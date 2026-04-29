# Load the packages needed to estimate fixed-effects count models, extract
# coefficients and handle the panel data.
pacman::p_load(
  broom,
  fixest,
  here,
  lubridate,
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

# Store the day labels used in the event-study model, ordered from the Thursday
# before Carnival to the Thursday after Carnival.
nhc_event_day_levels <- c(
  "outside_window",
  "thursday_before",
  "friday_before",
  "saturday_before",
  "carnival_sunday",
  "carnival_monday",
  "tuesday_after",
  "wednesday_after",
  "thursday_after"
)

# Create one row for each day in the Thursday-to-Thursday window around each
# real Carnival Monday in the panel.
nhc_event_window_dates <- nhc_panel |>
  distinct(crime_date, nhc_monday) |>
  filter(nhc_monday == 1) |>
  transmute(
    nhc_year = year(crime_date),
    nhc_monday_date = crime_date
  ) |>
  expand_grid(event_day_offset = -4:3) |>
  transmute(
    crime_date = nhc_monday_date + days(event_day_offset),
    nhc_year = nhc_year,
    event_day_offset = event_day_offset,
    nhc_event_day = recode_values(
      event_day_offset,
      -4 ~ "thursday_before",
      -3 ~ "friday_before",
      -2 ~ "saturday_before",
      -1 ~ "carnival_sunday",
      0 ~ "carnival_monday",
      1 ~ "tuesday_after",
      2 ~ "wednesday_after",
      3 ~ "thursday_after"
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

# MODEL FUNCTION ------------------------------------------------------------

# Estimate the event-study fixed-effects Poisson model for one crime type at a
# time and save it to the supplied model path.
fit_nhc_event_study_model <- function(crime_type, model_path) {
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

  # Save the fitted model so it can be reused without re-estimating.
  write_rds(model, model_path, compress = "gz")

  # Clean up memory used by the model fit before moving to the next crime group.
  gc()

  # Return a lightweight summary so the script does not hold all fitted model
  # objects in memory at the same time.
  tibble(crime_group = crime_type, model_path = model_path)
}

# RUN EVENT-STUDY MODELS ----------------------------------------------------

# Store the expected model path for each crime group so the script can resume
# after a failed or interrupted run without re-estimating completed models.
nhc_event_study_model_files <- tibble(
  crime_group = nhc_crime_groups,
  model_path = here(
    "derived_data",
    str_glue("nhc_model_event_study_{nhc_crime_groups}.rds")
  )
) |>
  mutate(model_exists = file.exists(model_path))

# Estimate and save only the event-study models that are not already present in
# the derived-data directory.
nhc_event_study_model_files |>
  filter(!model_exists) |>
  select(crime_group, model_path) |>
  pwalk(
    \(crime_group, model_path) {
      fit_nhc_event_study_model(crime_group, model_path)
    }
  ) |>
  invisible()

# Keep the completed model-file table used by the extraction steps.
nhc_event_study_model_files <- nhc_event_study_model_files |>
  select(crime_group, model_path)

# EXTRACT COEFFICIENTS AND MODEL STATISTICS --------------------------------

# Load and process one fitted model, then remove the model object from memory
# before returning compact summaries.
process_nhc_event_study_model <- function(crime_group_name, model_path) {
  # Load only the current crime-group model so previous and later models are not
  # kept in memory at the same time.
  saved_model <- read_rds(model_path)

  # Extract model coefficients with 95% confidence intervals.
  coef_results <- tidy(
    saved_model,
    conf.int = TRUE,
    conf.level = 0.95
  ) |>
    mutate(
      crime_group = crime_group_name,
      exponentiated_estimate = exp(estimate),
      exponentiated_conf_low = exp(conf.low),
      exponentiated_conf_high = exp(conf.high),
      .before = 1
    )

  # Store model fit statistics in a compact object for reporting.
  nhc_fit_stats <- fitstat(
    saved_model,
    c("ll", "aic", "bic", "rmse", "sq.cor", "pr2", "apr2")
  )

  # Test whether all estimated event-window coefficients are jointly equal to
  # zero.
  nhc_event_day_test <- wald(saved_model, keep = "nhc_event_day")

  # Store model-level statistics that are useful for manuscript tables and
  # appendix reporting.
  stats_results <- tibble(
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

  # Keep only the event-window rows for the current crime group because
  # outside-window dates are the counterfactual baseline.
  event_day_data <- nhc_event_panel |>
    filter(
      crime_group == crime_group_name,
      nhc_event_day %in% nhc_event_days
    )

  # Create a counterfactual version of the event-window rows where the same
  # dates and hexagons are treated as outside the event window.
  no_event_day_data <- event_day_data |>
    mutate(
      nhc_event_day = factor(
        "outside_window",
        levels = nhc_event_day_levels
      )
    )

  # Store the model coefficients for use in the extra-crime calculations.
  nhc_model_coefs <- coef(saved_model)

  # Store the model variance-covariance matrix so the London-wide intervals can
  # account for covariance between distance-band terms.
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
      extra_crimes = predicted_no_event_day * (exp(event_day_effect) - 1),
      extra_crimes_se = predicted_no_event_day *
        exp(event_day_effect) *
        event_day_effect_se,
      extra_crimes_conf_low = extra_crimes -
        qnorm(0.975) * extra_crimes_se,
      extra_crimes_conf_high = extra_crimes +
        qnorm(0.975) * extra_crimes_se
    )

  # Estimate London-wide additional crimes for each event day using the full
  # covariance matrix for the contributing distance-band terms.
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

  # Combine distance-band and London-wide extra-crime rows for the current
  # crime group.
  extra_results <- bind_rows(
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

  # Store compact outputs before removing the large model and prediction data.
  model_results <- tibble(
    crime_group = crime_group_name,
    coefficients = list(coef_results),
    model_stats = list(stats_results),
    extra_crimes = list(extra_results)
  )

  # Remove large objects from memory before the next model is loaded.
  rm(
    saved_model,
    event_day_data,
    no_event_day_data,
    extra_by_band,
    extra_all_london
  )
  gc()

  # Return only the compact summaries for the current crime group.
  model_results
}

# Load and process each event-study model in turn, retaining only compact
# summaries between iterations.
nhc_event_study_results <- nhc_event_study_model_files |>
  pmap(
    \(crime_group, model_path) {
      process_nhc_event_study_model(crime_group, model_path)
    }
  ) |>
  list_rbind()

# Combine the event-study coefficient summaries into one tidy table.
nhc_event_study_coef <- nhc_event_study_results |>
  select(coefficients) |>
  unnest(coefficients) |>
  arrange(crime_group, term)

# Combine the event-study model-statistic summaries into one tidy table.
nhc_event_study_stats <- nhc_event_study_results |>
  select(model_stats) |>
  unnest(model_stats) |>
  arrange(crime_group)

# Combine the event-day extra-crime summaries into one tidy table.
nhc_event_study_extra <- nhc_event_study_results |>
  select(extra_crimes) |>
  unnest(extra_crimes)

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

# Save the event-study extra-crime summaries.
write_rds(
  nhc_event_study_extra,
  here("derived_data", "nhc_extra_event_study.rds"),
  compress = "gz"
)
