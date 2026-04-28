# Load the packages needed to estimate fixed-effects count models, extract
# coefficients, run models in parallel and handle the panel data.
pacman::p_load(
  broom,
  fixest,
  furrr,
  future,
  here,
  lubridate,
  parallelly,
  tidyverse
)

# Increase the memory allocation available to R because the daily panel and
# fixed-effects model objects can be large.
mem.maxVSize(mem.maxVSize() * 4)

# LOAD PANEL DATA -----------------------------------------------------------

# Load the daily hexagon-level panel used by the main Notting Hill Carnival
# models.
nhc_panel <- read_rds(here("derived_data", "nhc_panel_350m_daily.rds"))

# Store the crime groups in the panel so the event-study models match the
# current panel definition.
nhc_crime_groups <- nhc_panel |>
  distinct(crime_group) |>
  arrange(crime_group) |>
  pull(crime_group) |>
  as.character()

# DEFINE EVENT-STUDY WINDOW -------------------------------------------------

# Store the day labels used in the event-study model, ordered from the Friday
# before Carnival to the Wednesday after Carnival.
nhc_event_day_levels <- c(
  "outside_window",
  "friday_before",
  "saturday_before",
  "carnival_sunday",
  "carnival_monday",
  "tuesday_after",
  "wednesday_after"
)

# Create one row for each day in the Friday-to-Wednesday window around each
# real Carnival Monday in the panel.
nhc_event_window_dates <- nhc_panel |>
  distinct(crime_date, nhc_monday) |>
  filter(nhc_monday == 1) |>
  transmute(
    nhc_year = year(crime_date),
    nhc_monday_date = crime_date
  ) |>
  expand_grid(event_day_offset = -3:2) |>
  transmute(
    crime_date = nhc_monday_date + days(event_day_offset),
    nhc_year = nhc_year,
    event_day_offset = event_day_offset,
    nhc_event_day = recode_values(
      event_day_offset,
      -3 ~ "friday_before",
      -2 ~ "saturday_before",
      -1 ~ "carnival_sunday",
      0 ~ "carnival_monday",
      1 ~ "tuesday_after",
      2 ~ "wednesday_after"
    )
  )

# Add the event-study day label to every row in the panel, with all other dates
# retained as the omitted outside-window category.
nhc_event_panel <- nhc_panel |>
  left_join(nhc_event_window_dates, by = join_by(crime_date)) |>
  mutate(
    nhc_event_day = replace_na(nhc_event_day, "outside_window"),
    nhc_event_day = factor(nhc_event_day, levels = nhc_event_day_levels)
  )

# Store the event-window days that are explicitly estimated in the model.
nhc_event_days <- nhc_event_day_levels |>
  setdiff("outside_window")

# SET UP PARALLEL MODELLING -------------------------------------------------

# Use a conservative number of workers because each fixed-effects model can use
# a large amount of memory while it is being estimated and saved.
nhc_event_study_workers <- min(3L, max(1L, as.integer(availableCores()) - 1L))

# Allow the parallel backend to register the large in-memory panel as a global.
options(future.globals.maxSize = 16 * 1024^3)

# MODEL FUNCTION ------------------------------------------------------------

# Estimate the event-study fixed-effects Poisson model for one crime type at a
# time.
fit_nhc_event_study_model <- function(crime_type) {
  # Keep the rows for the requested crime type so each model has the same
  # specification but a different outcome subset.
  model_data <- nhc_event_panel |>
    filter(crime_group == crime_type)

  # Estimate distance-band effects separately for each day in the event window,
  # using the 12-kilometre band and all outside-window dates as references.
  model <- fepois(
    crime_count ~
      i(dist, nhc_event_day, ref = "12km", ref2 = "outside_window") |
      hex_id + crime_date,
    data = model_data,
    cluster = ~hex_id
  )

  # Store the path for the saved event-study model so it follows the project
  # naming convention for reusable model outputs.
  model_path <- here(
    "derived_data",
    str_glue("nhc_model_event_study_{crime_type}.rds")
  )

  # Save the fitted model so it can be reused without re-estimating.
  write_rds(model, model_path, compress = "gz")

  # Clean up memory used by the model fit inside the worker process.
  gc()

  # Return a lightweight summary so the parent process does not have to hold all
  # fitted model objects in memory at the same time.
  tibble(crime_group = crime_type, model_path = model_path)
}

# RUN EVENT-STUDY MODELS ----------------------------------------------------

# Use forked workers when available so the large panel can be shared more
# efficiently between worker processes.
if (supportsMulticore()) {
  plan(multicore, workers = nhc_event_study_workers)
} else {
  plan(multisession, workers = nhc_event_study_workers)
}

# Estimate and save one event-study fixed-effects Poisson model for each crime
# group in the panel dataset.
nhc_event_study_model_files <- nhc_crime_groups |>
  future_map(
    fit_nhc_event_study_model,
    .options = furrr_options(seed = NULL)
  ) |>
  list_rbind()

# Return the future backend to sequential processing after the model run.
plan(sequential)

# EXTRACT COEFFICIENTS AND MODEL STATISTICS --------------------------------

# Extract coefficient and model-statistic summaries for each fitted event-study
# model.
nhc_event_study_results <- nhc_event_study_model_files |>
  mutate(
    model = map(model_path, read_rds),
    coefficients = map2(
      model,
      crime_group,
      \(saved_model, crime_group_name) {
        tidy(saved_model, conf.int = TRUE, conf.level = 0.95) |>
          mutate(crime_group = crime_group_name, .before = 1)
      }
    ),
    model_stats = map2(
      model,
      crime_group,
      \(saved_model, crime_group_name) {
        nhc_fit_stats <- fitstat(
          saved_model,
          c("ll", "aic", "bic", "rmse", "sq.cor", "pr2", "apr2")
        )

        nhc_event_day_test <- wald(saved_model, keep = "nhc_event_day")

        tibble(
          crime_group = crime_group_name,
          observations_original = saved_model$nobs_origin,
          observations_dropped = saved_model$nobs_origin - saved_model$nobs,
          observations_remaining = saved_model$nobs,
          dropped_terms = list(saved_model$collin.var),
          hexagon_fixed_effects = saved_model$fixef_sizes[["hex_id"]],
          date_fixed_effects = saved_model$fixef_sizes[["crime_date"]],
          log_likelihood = unname(nhc_fit_stats$ll),
          aic = unname(nhc_fit_stats$aic),
          bic = unname(nhc_fit_stats$bic),
          rmse = unname(nhc_fit_stats$rmse),
          squared_correlation = unname(nhc_fit_stats$sq.cor),
          pseudo_r2 = unname(nhc_fit_stats$pr2),
          adjusted_pseudo_r2 = unname(nhc_fit_stats$apr2),
          joint_event_window_test_statistic = unname(nhc_event_day_test$stat),
          joint_event_window_test_p_value = unname(nhc_event_day_test$p),
          joint_event_window_test_df1 = unname(nhc_event_day_test$df1),
          joint_event_window_test_df2 = unname(nhc_event_day_test$df2),
          joint_event_window_test_vcov = nhc_event_day_test$vcov
        )
      }
    )
  )

# Combine the event-study coefficient summaries into one tidy table.
nhc_event_study_coef <- nhc_event_study_results |>
  select(coefficients) |>
  unnest(coefficients) |>
  mutate(
    exponentiated_estimate = exp(estimate),
    exponentiated_conf_low = exp(conf.low),
    exponentiated_conf_high = exp(conf.high)
  ) |>
  arrange(crime_group, term)

# Combine the event-study model-statistic summaries into one tidy table.
nhc_event_study_stats <- nhc_event_study_results |>
  select(model_stats) |>
  unnest(model_stats) |>
  arrange(crime_group)

# Save the event-study coefficient summaries.
write_rds(
  nhc_event_study_coef,
  here("derived_data", "nhc_coef_event_study.rds"),
  compress = "gz"
)

# Save the event-study model-statistic summaries.
write_rds(
  nhc_event_study_stats,
  here("derived_data", "nhc_stats_event_study.rds"),
  compress = "gz"
)

# ESTIMATE EXTRA CRIMES BY EVENT DAY ----------------------------------------

# Estimate event-day counterfactual crime counts for each fitted event-study
# model and scale them to one annual instance of each event-window day.
nhc_event_study_extra <- nhc_event_study_results |>
  select(crime_group, model) |>
  mutate(
    extra_crimes = map2(
      crime_group,
      model,
      \(crime_group_name, saved_model) {
        # Keep only the event-window rows for the current crime group because
        # outside-window dates are the counterfactual baseline.
        event_day_data <- nhc_event_panel |>
          filter(
            crime_group == crime_group_name,
            nhc_event_day %in% nhc_event_days
          )

        # Create a counterfactual version of the event-window rows where the
        # same dates and hexagons are treated as outside the event window.
        no_event_day_data <- event_day_data |>
          mutate(
            nhc_event_day = factor(
              "outside_window",
              levels = nhc_event_day_levels
            )
          )

        # Store the model coefficients for use in the extra-crime calculations.
        nhc_model_coefs <- coef(saved_model)

        # Store the model variance-covariance matrix so the London-wide
        # intervals can account for covariance between distance-band terms.
        nhc_model_vcov <- vcov(saved_model)

        # Estimate additional crimes by event day and distance band.
        extra_by_band <- event_day_data |>
          mutate(
            predicted_event_day = predict(
              saved_model,
              newdata = event_day_data,
              type = "response"
            ),
            predicted_no_event_day = predict(
              saved_model,
              newdata = no_event_day_data,
              type = "response"
            ),
            predicted_event_day = replace_na(predicted_event_day, 0),
            predicted_no_event_day = replace_na(predicted_no_event_day, 0)
          ) |>
          summarise(
            annual_scale = 1 / n_distinct(crime_date),
            observed_event_day = sum(crime_count) * annual_scale,
            predicted_event_day = sum(predicted_event_day) * annual_scale,
            predicted_no_event_day = sum(predicted_no_event_day) * annual_scale,
            .by = c(event_day_offset, nhc_event_day, dist_km, dist)
          ) |>
          mutate(
            model_term = if_else(
              dist == "12km",
              NA_character_,
              str_glue("dist::{dist}:nhc_event_day::{nhc_event_day}")
            ),
            model_term_available = model_term %in% names(nhc_model_coefs),
            event_day_effect = if_else(
              is.na(model_term) | !model_term_available,
              0,
              nhc_model_coefs[model_term]
            ),
            event_day_effect_se = if_else(
              is.na(model_term) | !model_term_available,
              0,
              sqrt(diag(nhc_model_vcov)[model_term])
            ),
            extra_crimes = predicted_no_event_day *
              (exp(event_day_effect) - 1),
            extra_crimes_se = predicted_no_event_day *
              exp(event_day_effect) *
              event_day_effect_se,
            extra_crimes_conf_low = extra_crimes -
              qnorm(0.975) * extra_crimes_se,
            extra_crimes_conf_high = extra_crimes +
              qnorm(0.975) * extra_crimes_se
          )

        # Estimate London-wide additional crimes for each event day using the
        # full covariance matrix for the contributing distance-band terms.
        extra_all_london <- extra_by_band |>
          filter(model_term_available) |>
          summarise(
            model_terms = list(model_term),
            total_extra_crimes_gradient = list(
              predicted_no_event_day * exp(event_day_effect)
            ),
            dist_km = NA_real_,
            dist = "all_london",
            observed_event_day = sum(observed_event_day),
            predicted_event_day = sum(predicted_event_day),
            predicted_no_event_day = sum(predicted_no_event_day),
            extra_crimes = sum(extra_crimes),
            .by = c(event_day_offset, nhc_event_day)
          ) |>
          mutate(
            total_extra_crimes_se = map2_dbl(
              model_terms,
              total_extra_crimes_gradient,
              \(model_terms, total_extra_crimes_gradient) {
                sqrt(as.numeric(
                  t(total_extra_crimes_gradient) %*%
                    nhc_model_vcov[model_terms, model_terms, drop = FALSE] %*%
                    total_extra_crimes_gradient
                ))
              }
            ),
            extra_crimes_conf_low = extra_crimes -
              qnorm(0.975) * total_extra_crimes_se,
            extra_crimes_conf_high = extra_crimes +
              qnorm(0.975) * total_extra_crimes_se
          ) |>
          select(
            event_day_offset,
            nhc_event_day,
            dist_km,
            dist,
            observed_event_day,
            predicted_event_day,
            predicted_no_event_day,
            extra_crimes,
            extra_crimes_conf_low,
            extra_crimes_conf_high
          )

        # Combine distance-band and London-wide rows for the current crime
        # group.
        bind_rows(
          extra_by_band |>
            mutate(dist = as.character(dist)) |>
            select(
              event_day_offset,
              nhc_event_day,
              dist_km,
              dist,
              observed_event_day,
              predicted_event_day,
              predicted_no_event_day,
              extra_crimes,
              extra_crimes_conf_low,
              extra_crimes_conf_high
            ),
          extra_all_london
        ) |>
          mutate(crime_group = crime_group_name, .before = 1) |>
          arrange(event_day_offset, dist_km)
      }
    )
  ) |>
  select(extra_crimes) |>
  unnest(extra_crimes)

# Save the event-study extra-crime summaries.
write_rds(
  nhc_event_study_extra,
  here("derived_data", "nhc_extra_event_study.rds"),
  compress = "gz"
)
