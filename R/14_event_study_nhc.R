# Load the packages needed to estimate fixed-effects count models, extract
# coefficients and handle the panel data.
pacman::p_load(
  broom,
  cli,
  fixest,
  here,
  lubridate,
  tidyverse
)

# Increase the memory allocation available to R because the daily panel and
# fixed-effects model objects can be large.
mem.maxVSize(mem.maxVSize() * 6)

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

# Keep the requested crime group, limiting sexual-offence event-study models to
# the years with usable geographic co-ordinates.
filter_nhc_event_panel <- function(panel_data, crime_type) {
  panel_data |>
    filter(
      crime_group == crime_type,
      crime_type != "sexual_offences" |
        between(year(crime_date), 2013, 2019)
    )
}

# Store the event-window days that are explicitly estimated in the model.
nhc_event_days <- nhc_event_day_levels |>
  setdiff("outside_window")

# Count unexpected missing values in an object, including values stored inside
# list columns, so saved result files can be checked before large model files
# are deleted.
count_missing_values <- function(result_object) {
  # Ignore the expected missing distance value on the London-wide event-study
  # rows, because those rows aggregate across all distance bands.
  if (
    is.data.frame(result_object) &&
      all(c("dist", "dist_km") %in% names(result_object))
  ) {
    result_object <- result_object |>
      mutate(
        dist_km = if_else(dist == "all_london", 0, dist_km)
      )
  }

  # Recurse into list-like objects so nested vectors, matrices and data frames
  # are all included in the missing-value count.
  if (is.list(result_object)) {
    return(sum(map_int(result_object, count_missing_values)))
  }

  # Count missing values in atomic objects once all nested structures have been
  # traversed.
  sum(is.na(result_object))
}

# MODEL FUNCTION ------------------------------------------------------------

# Estimate the event-study fixed-effects Poisson model for one crime type at a
# time and save it to the supplied model path.
fit_nhc_event_study_model <- function(crime_type, model_path) {
  # Skip estimation when the fitted event-study model is already available on
  # disk.
  if (file.exists(model_path)) {
    return(tibble(crime_group = crime_type, model_path = model_path))
  }

  # Keep the rows for the requested crime type so each model has the same
  # specification but a different outcome subset.
  model_data <- nhc_event_panel |>
    filter_nhc_event_panel(crime_type)

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
  # Record whether each event-study model already exists on disk.
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

  # Test whether the estimated pre-Carnival coefficients are jointly equal to
  # zero, as a no-anticipation check for the event-study specification.
  nhc_no_anticipation_test <- wald(
    saved_model,
    keep = "thursday_before|friday_before|saturday_before"
  )

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
    joint_event_window_test_vcov = nhc_event_day_test$vcov,
    no_anticipation_test_statistic = unname(nhc_no_anticipation_test$stat),
    no_anticipation_test_p_value = unname(nhc_no_anticipation_test$p),
    no_anticipation_test_df1 = unname(nhc_no_anticipation_test$df1),
    no_anticipation_test_df2 = unname(nhc_no_anticipation_test$df2),
    no_anticipation_test_vcov = nhc_no_anticipation_test$vcov
  )

  # Keep only the event-window rows for the current crime group because
  # outside-window dates are the counterfactual baseline.
  event_day_data <- nhc_event_panel |>
    filter_nhc_event_panel(crime_group_name) |>
    filter(
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

# Store the event-study result paths so the same saved files can be checked
# before the large model files are deleted.
nhc_event_study_result_files <- tibble(
  result_file = c(
    "nhc_coef_event_study.rds",
    "nhc_stats_event_study.rds",
    "nhc_extra_event_study.rds"
  ),
  result_path = here("derived_data", result_file)
)

# Save the event-study coefficient summaries.
write_rds(
  nhc_event_study_coef,
  nhc_event_study_result_files |>
    filter(result_file == "nhc_coef_event_study.rds") |>
    pull(result_path),
  compress = "gz"
)

# Save the event-study model-statistic summaries.
write_rds(
  nhc_event_study_stats,
  nhc_event_study_result_files |>
    filter(result_file == "nhc_stats_event_study.rds") |>
    pull(result_path),
  compress = "gz"
)

# Save the event-study extra-crime summaries.
write_rds(
  nhc_event_study_extra,
  nhc_event_study_result_files |>
    filter(result_file == "nhc_extra_event_study.rds") |>
    pull(result_path),
  compress = "gz"
)

# Check each saved event-study result file for missing values before deciding
# whether it is safe to remove the reusable model objects.
nhc_event_study_missing_values <- nhc_event_study_result_files |>
  mutate(
    missing_values = map_int(
      result_path,
      \(result_path) {
        read_rds(result_path) |>
          count_missing_values()
      }
    )
  ) |>
  filter(missing_values > 0)

# Remove the saved event-study model objects only if every saved result file is
# complete, otherwise keep the models so the extraction can be investigated
# without re-estimating the expensive fixed-effects models.
if (nrow(nhc_event_study_missing_values) == 0) {
  unlink(nhc_event_study_model_files$model_path)
} else {
  cli_warn(c(
    "Event-study model files were kept because missing values were found in the saved result files.",
    set_names(
      str_c(
        nhc_event_study_missing_values$missing_values,
        " missing value",
        if_else(
          nhc_event_study_missing_values$missing_values == 1,
          "",
          "s"
        ),
        " in ",
        nhc_event_study_missing_values$result_file,
        "."
      ),
      rep("x", nrow(nhc_event_study_missing_values))
    )
  ))
}

# EXTENDED CRIMINAL-DAMAGE EVENT STUDY -------------------------------------

# Store the extended event-study labels from 10 days before Carnival Monday to
# three days after Carnival Monday.
nhc_criminal_damage_extended_event_days <- tibble(
  event_day_offset = -10:3,
  nhc_extended_event_day = c(
    "ten_days_before",
    "nine_days_before",
    "eight_days_before",
    "seven_days_before",
    "six_days_before",
    "five_days_before",
    "four_days_before",
    "three_days_before",
    "two_days_before",
    "carnival_sunday",
    "carnival_monday",
    "one_day_after",
    "two_days_after",
    "three_days_after"
  )
)

# Update the event-study factor levels for the extended criminal-damage model,
# keeping all dates outside the extended window as the omitted reference period.
nhc_event_day_levels <- c(
  "outside_window",
  nhc_criminal_damage_extended_event_days$nhc_extended_event_day
)

# Create one row for every date in the extended event-study window around each
# real Carnival Monday in the panel.
nhc_criminal_damage_extended_window_dates <- nhc_panel |>
  distinct(crime_date, nhc_monday) |>
  filter(nhc_monday == 1) |>
  transmute(
    nhc_year = year(crime_date),
    nhc_monday_date = crime_date
  ) |>
  expand_grid(nhc_criminal_damage_extended_event_days) |>
  transmute(
    crime_date = nhc_monday_date + days(event_day_offset),
    nhc_year = nhc_year,
    event_day_offset = event_day_offset,
    nhc_extended_event_day = nhc_extended_event_day
  )

# Add the extended event-study label to the criminal-damage panel using the same
# event-day column name as the main event-study helper functions.
nhc_event_panel <- nhc_panel |>
  filter(crime_group == "criminal_damage") |>
  left_join(
    nhc_criminal_damage_extended_window_dates,
    by = join_by(crime_date)
  ) |>
  mutate(
    nhc_event_day = replace_na(
      nhc_extended_event_day,
      "outside_window"
    ),
    nhc_event_day = factor(
      nhc_event_day,
      levels = nhc_event_day_levels
    )
  ) |>
  select(-nhc_extended_event_day)

# Store the extended event-window days that are explicitly estimated in the
# criminal-damage model.
nhc_event_days <- nhc_event_day_levels |>
  setdiff("outside_window")

# Store the path for the extended criminal-damage event-study model so the
# script can resume without re-estimating the model if it already exists.
nhc_criminal_damage_extended_model_file <- here(
  "derived_data",
  "nhc_model_event_study_criminal_damage_extended.rds"
)

# Estimate and save the extended criminal-damage model only if a saved copy is
# not already available.
if (!file.exists(nhc_criminal_damage_extended_model_file)) {
  fit_nhc_event_study_model(
    "criminal_damage",
    nhc_criminal_damage_extended_model_file
  )
}

# Load and process the extended criminal-damage model using the same extraction
# workflow as the main event-study models.
nhc_criminal_damage_extended_results <- process_nhc_event_study_model(
  "criminal_damage",
  nhc_criminal_damage_extended_model_file
)

# Unnest the extended criminal-damage coefficient summaries into one tidy table.
nhc_criminal_damage_extended_coef <- nhc_criminal_damage_extended_results |>
  select(coefficients) |>
  unnest(coefficients) |>
  arrange(term)

# Unnest the extended criminal-damage model-statistic summaries into one tidy
# table.
nhc_criminal_damage_extended_stats <- nhc_criminal_damage_extended_results |>
  select(model_stats) |>
  unnest(model_stats) |>
  mutate(
    event_window_start = -10,
    event_window_end = 3,
    .after = crime_group
  ) |>
  arrange(crime_group)

# Unnest the extended criminal-damage extra-crime summaries into one tidy table.
nhc_criminal_damage_extended_extra <- nhc_criminal_damage_extended_results |>
  select(extra_crimes) |>
  unnest(extra_crimes)

# Store the extended criminal-damage result paths so the saved files can be
# checked before removing the large model object.
nhc_criminal_damage_extended_result_files <- tibble(
  result_file = c(
    "nhc_coef_event_study_criminal_damage_extended.rds",
    "nhc_stats_event_study_criminal_damage_extended.rds",
    "nhc_extra_event_study_criminal_damage_extended.rds"
  ),
  result_path = here("derived_data", result_file)
)

# Save the extended criminal-damage coefficient summaries.
write_rds(
  nhc_criminal_damage_extended_coef,
  nhc_criminal_damage_extended_result_files |>
    filter(
      result_file == "nhc_coef_event_study_criminal_damage_extended.rds"
    ) |>
    pull(result_path),
  compress = "gz"
)

# Save the extended criminal-damage model-statistic summaries.
write_rds(
  nhc_criminal_damage_extended_stats,
  nhc_criminal_damage_extended_result_files |>
    filter(
      result_file == "nhc_stats_event_study_criminal_damage_extended.rds"
    ) |>
    pull(result_path),
  compress = "gz"
)

# Save the extended criminal-damage extra-crime summaries.
write_rds(
  nhc_criminal_damage_extended_extra,
  nhc_criminal_damage_extended_result_files |>
    filter(
      result_file == "nhc_extra_event_study_criminal_damage_extended.rds"
    ) |>
    pull(result_path),
  compress = "gz"
)

# Check the extended criminal-damage result files for missing values before
# deciding whether the reusable model object can be removed.
nhc_criminal_damage_extended_missing_values <-
  nhc_criminal_damage_extended_result_files |>
    mutate(
      missing_values = map_int(
        result_path,
        \(result_path) {
          read_rds(result_path) |>
            count_missing_values()
        }
      )
    ) |>
    filter(missing_values > 0)

# Remove the extended criminal-damage model object only if every saved result
# file is complete.
if (nrow(nhc_criminal_damage_extended_missing_values) == 0) {
  unlink(nhc_criminal_damage_extended_model_file)
} else {
  cli_warn(c(
    str_c(
      "Extended criminal-damage event-study model file was kept because",
      " missing values were found in the saved result files."
    ),
    set_names(
      str_c(
        nhc_criminal_damage_extended_missing_values$missing_values,
        " missing value",
        if_else(
          nhc_criminal_damage_extended_missing_values$missing_values == 1,
          "",
          "s"
        ),
        " in ",
        nhc_criminal_damage_extended_missing_values$result_file,
        "."
      ),
      rep("x", nrow(nhc_criminal_damage_extended_missing_values))
    )
  ))
}
