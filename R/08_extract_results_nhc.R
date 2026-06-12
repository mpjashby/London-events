# This code extracts coefficient and counterfactual extra-crime estimates from
# the saved Notting Hill Carnival crime-count models.

# This code might need a larger memory allocation.
mem.maxVSize(mem.maxVSize() * 4)

# Load the packages needed to load models, extract estimates and manipulate the
# panel data.
pacman::p_load(fixest, broom, here, janitor, tidyverse)

# LOAD PANEL DATA -----------------------------------------------------------

# Load the daily hexagon-level panel and standardise column names before any
# prediction operations.
nhc_panel <- read_rds(here("derived_data", "nhc_panel_350m_daily.rds"))

# LIST MODEL INPUTS ---------------------------------------------------------

# Store the current crime groups from the panel so that result extraction only
# uses models matching the latest panel definition.
nhc_crime_groups <- nhc_panel |>
  distinct(crime_group) |>
  arrange(crime_group) |>
  pull(crime_group) |>
  as.character()

# FILTER MODEL DATA ---------------------------------------------------------

# Keep the requested crime group, limiting sexual-offence model summaries and
# predictions to the years with usable geographic co-ordinates.
filter_nhc_model_panel <- function(panel_data, crime_group_name) {
  panel_data |>
    filter(
      crime_group == crime_group_name,
      crime_group_name != "sexual_offences" |
        between(year(crime_date), 2013, 2019)
    )
}

# EXTRACT RESULTS -----------------------------------------------------------

# Extract and save results for each fitted crime-type model.
nhc_crime_groups |>
  walk(function(crime_group_name) {
    # Load the fitted model for the current crime group.
    nhc_model <- read_rds(
      here("derived_data", str_glue("nhc_model_{crime_group_name}.rds"))
    )

    # Extract model coefficients with 95% confidence intervals.
    coef_results <- tidy(
      nhc_model,
      conf.int = TRUE,
      conf.level = 0.95
    ) |>
      # Add the crime-group label to the coefficient table.
      mutate(crime_group = crime_group_name, .before = 1)

    # Save the coefficient table for the current crime group.
    write_rds(
      coef_results,
      here("derived_data", str_glue("nhc_coef_{crime_group_name}.rds")),
      compress = "gz"
    )

    # Store model fit statistics in a compact object for reporting in the paper.
    nhc_fit_stats <- fitstat(
      nhc_model,
      c("ll", "aic", "bic", "rmse", "sq.cor", "pr2", "apr2")
    )

    # Test whether all estimated carnival distance-band coefficients are jointly
    # equal to zero.
    nhc_band_test <- wald(nhc_model, keep = "is_nhc")

    # Store the fixed-effect levels dropped by the model because they had only
    # zero outcomes.
    nhc_dropped_fixed_effects <- tibble(
      fixed_effect = names(nhc_model$fixef_removed),
      dropped_level = nhc_model$fixef_removed
    ) |>
      # Store all dropped fixed-effect levels as labels so different fixed-effect
      # types can be combined in one tidy table.
      mutate(dropped_level = map(dropped_level, as.character)) |>
      unnest_longer(dropped_level)

    # Store a compact count of rows removed at each model-selection step.
    nhc_dropped_row_summary <- enframe(
      nhc_model$obs_selection,
      name = "drop_step",
      value = "row_indices"
    ) |>
      mutate(rows_dropped = lengths(row_indices)) |>
      select(drop_step, rows_dropped)

    # Store model-level statistics that are useful for manuscript tables and
    # appendix reporting.
    stats_results <- tibble(
      crime_group = crime_group_name,
      observations_original = nhc_model$nobs_origin,
      observations_dropped = nhc_model$nobs_origin - nhc_model$nobs,
      observations_remaining = nhc_model$nobs,
      dropped_row_summary = list(nhc_dropped_row_summary),
      dropped_fixed_effects = list(nhc_dropped_fixed_effects),
      dropped_terms = list(nhc_model$collin.var),
      hexagon_fixed_effects = nhc_model$fixef_sizes[["hex_id"]],
      date_fixed_effects = nhc_model$fixef_sizes[["crime_date"]],
      log_likelihood = unname(nhc_fit_stats$ll),
      aic = unname(nhc_fit_stats$aic),
      bic = unname(nhc_fit_stats$bic),
      rmse = unname(nhc_fit_stats$rmse),
      squared_correlation = unname(nhc_fit_stats$sq.cor),
      pseudo_r2 = unname(nhc_fit_stats$pr2),
      adjusted_pseudo_r2 = unname(nhc_fit_stats$apr2),
      joint_band_test_statistic = unname(nhc_band_test$stat),
      joint_band_test_p_value = unname(nhc_band_test$p),
      joint_band_test_df1 = unname(nhc_band_test$df1),
      joint_band_test_df2 = unname(nhc_band_test$df2),
      joint_band_test_vcov = nhc_band_test$vcov
    )

    # Save the compact model statistics for the current crime group.
    write_rds(
      stats_results,
      here("derived_data", str_glue("nhc_stats_{crime_group_name}.rds")),
      compress = "gz"
    )

    # Keep the panel rows that correspond to the modelled years for the current
    # crime group.
    nhc_model_panel <- filter_nhc_model_panel(nhc_panel, crime_group_name)

    # Keep the carnival-day rows for the current crime group.
    nhc_day_data <- nhc_model_panel |>
      filter(is_nhc == 1)

    # Calculate the mean observed daily crime count on non-Carnival days in
    # each distance band for the current crime group.
    non_nhc_observed_results <- nhc_model_panel |>
      filter(is_nhc == 0) |>
      summarise(
        daily_crime_count = sum(crime_count),
        .by = c(crime_date, dist_km, dist)
      ) |>
      summarise(
        observed_non_nhc = mean(daily_crime_count),
        .by = c(dist_km, dist)
      )

    # Convert estimates across all carnival days in the data to estimates for a
    # single annual two-day occurrence of Carnival.
    carnival_annual_scale <- 2 / n_distinct(nhc_day_data$crime_date)

    # Create a counterfactual version of the carnival-day rows where carnival did
    # not occur on the same dates in the same hexagons.
    no_nhc_day_data <- nhc_day_data |>
      mutate(is_nhc = 0)

    # Store the model coefficients for use in the extra-crime calculations.
    nhc_model_coefs <- coef(nhc_model)

    # Store the model variance-covariance matrix so the London-wide interval can
    # account for covariance between distance-band coefficients.
    nhc_model_vcov <- vcov(nhc_model)

    # Estimate the additional crimes predicted on carnival days by distance band.
    extra_results <- nhc_day_data |>
      mutate(
        predicted_nhc = predict(
          nhc_model,
          newdata = nhc_day_data,
          type = "response"
        ),
        predicted_no_nhc = predict(
          nhc_model,
          newdata = no_nhc_day_data,
          type = "response"
        ),
        # Treat predictions for model-dropped all-zero fixed effects as zero.
        predicted_nhc = replace_na(predicted_nhc, 0),
        predicted_no_nhc = replace_na(predicted_no_nhc, 0)
      ) |>
      summarise(
        observed_nhc = sum(crime_count) * carnival_annual_scale,
        predicted_nhc = sum(predicted_nhc) * carnival_annual_scale,
        predicted_no_nhc = sum(predicted_no_nhc) * carnival_annual_scale,
        .by = c(dist_km, dist)
      ) |>
      # Link each distance band to the corresponding model coefficient.
      mutate(
        model_term = if_else(
          dist == "12km",
          NA_character_,
          str_glue("dist::{dist}:is_nhc")
        ),
        nhc_effect = if_else(
          is.na(model_term),
          0,
          nhc_model_coefs[model_term]
        ),
        nhc_effect_se = if_else(
          is.na(model_term),
          0,
          sqrt(diag(nhc_model_vcov)[model_term])
        ),
        extra_crimes = predicted_no_nhc * (exp(nhc_effect) - 1),
        extra_crimes_se = predicted_no_nhc * exp(nhc_effect) * nhc_effect_se,
        extra_crimes_conf_low = extra_crimes -
          qnorm(0.975) * extra_crimes_se,
        extra_crimes_conf_high = extra_crimes +
          qnorm(0.975) * extra_crimes_se
      ) |>
      # Add the observed non-Carnival daily crime count to match the observed
      # Carnival count already included in the results.
      left_join(
        non_nhc_observed_results,
        by = join_by(dist_km, dist),
        unmatched = "error"
      )

    # Store the model terms that contribute uncertainty to the London-wide total.
    total_model_terms <- extra_results |>
      filter(!is.na(model_term)) |>
      pull(model_term)

    # Store the derivative of the London-wide extra-crime estimate with respect
    # to each distance-band coefficient.
    total_extra_crimes_gradient <- extra_results |>
      filter(!is.na(model_term)) |>
      mutate(extra_crimes_gradient = predicted_no_nhc * exp(nhc_effect)) |>
      pull(extra_crimes_gradient)

    # Calculate the standard error for the London-wide total using the full
    # covariance matrix for the distance-band coefficients.
    total_extra_crimes_se <- sqrt(as.numeric(
      t(total_extra_crimes_gradient) %*%
        nhc_model_vcov[total_model_terms, total_model_terms, drop = FALSE] %*%
        total_extra_crimes_gradient
    ))

    # Add the crime-group label and the London-wide total row to the results.
    extra_results <- bind_rows(
      extra_results |>
        mutate(
          crime_group = crime_group_name,
          dist = as.character(dist),
          .before = 1
        ),
      extra_results |>
        summarise(
          crime_group = crime_group_name,
          dist_km = NA_real_,
          dist = "all_london",
          extra_crimes = sum(extra_crimes),
          extra_crimes_conf_low = extra_crimes -
            qnorm(0.975) * total_extra_crimes_se,
          extra_crimes_conf_high = extra_crimes +
            qnorm(0.975) * total_extra_crimes_se,
          observed_nhc = sum(observed_nhc),
          observed_non_nhc = sum(observed_non_nhc),
          predicted_nhc = sum(predicted_nhc),
          predicted_no_nhc = sum(predicted_no_nhc)
        )
    ) |>
      # Keep only the fields needed in the saved extra-crime results.
      select(
        crime_group,
        dist_km,
        dist,
        extra_crimes,
        extra_crimes_conf_low,
        extra_crimes_conf_high,
        observed_nhc,
        observed_non_nhc,
        predicted_nhc,
        predicted_no_nhc
      ) |>
      # Order distance bands from the footprint out, with the total row last.
      arrange(dist_km)

    # Save the extra-crime estimates for the current crime group.
    write_rds(
      extra_results,
      here("derived_data", str_glue("nhc_extra_{crime_group_name}.rds")),
      compress = "gz"
    )
  })

# EXTRACT PLACEBO CARNIVAL-DAY RESULTS -------------------------------------

# Store the last Sunday in July for each year in the panel as the placebo
# Carnival Sunday.
nhc_placebo_dates <- nhc_panel |>
  distinct(crime_date) |>
  filter(month(crime_date) == 7, wday(crime_date, week_start = 1) == 7) |>
  mutate(placebo_year = year(crime_date)) |>
  slice_max(crime_date, n = 1, by = placebo_year) |>
  transmute(
    placebo_year,
    placebo_sunday_date = crime_date,
    placebo_monday_date = crime_date + days(1)
  ) |>
  # Keep the placebo years aligned with the real Carnival treatment years by
  # excluding the years in which Carnival did not occur in its usual form.
  filter(!placebo_year %in% 2020:2021)

# Add placebo Carnival indicators to the panel for extracting counterfactual
# predictions from the saved placebo models.
nhc_panel <- nhc_panel |>
  mutate(
    nhc_placebo_sunday = as.numeric(crime_date %in% nhc_placebo_dates$placebo_sunday_date),
    nhc_placebo_monday = as.numeric(crime_date %in% nhc_placebo_dates$placebo_monday_date),
    is_nhc_placebo = as.numeric(nhc_placebo_sunday == 1 | nhc_placebo_monday == 1)
  )

# Extract and save results for each fitted placebo crime-type model.
nhc_crime_groups |>
  walk(function(crime_group_name) {
    # Load the fitted placebo model for the current crime group.
    nhc_placebo_model <- read_rds(
      here("derived_data", str_glue("nhc_model_placebo_{crime_group_name}.rds"))
    )

    # Extract placebo model coefficients with 95% confidence intervals.
    placebo_coef_results <- tidy(
      nhc_placebo_model,
      conf.int = TRUE,
      conf.level = 0.95
    ) |>
      # Add the crime-group label to the coefficient table.
      mutate(crime_group = crime_group_name, .before = 1)

    # Save the placebo coefficient table for the current crime group.
    write_rds(
      placebo_coef_results,
      here("derived_data", str_glue("nhc_coef_placebo_{crime_group_name}.rds")),
      compress = "gz"
    )

    # Store model fit statistics in a compact object for reporting in the paper.
    nhc_placebo_fit_stats <- fitstat(
      nhc_placebo_model,
      c("ll", "aic", "bic", "rmse", "sq.cor", "pr2", "apr2")
    )

    # Test whether all estimated placebo distance-band coefficients are jointly
    # equal to zero.
    nhc_placebo_band_test <- wald(nhc_placebo_model, keep = "is_nhc_placebo")

    # Store the fixed-effect levels dropped by the model because they had only
    # zero outcomes.
    nhc_placebo_dropped_fixed_effects <- tibble(
      fixed_effect = names(nhc_placebo_model$fixef_removed),
      dropped_level = nhc_placebo_model$fixef_removed
    ) |>
      # Store all dropped fixed-effect levels as labels so different fixed-effect
      # types can be combined in one tidy table.
      mutate(dropped_level = map(dropped_level, as.character)) |>
      unnest_longer(dropped_level)

    # Store a compact count of rows removed at each model-selection step.
    nhc_placebo_dropped_row_summary <- enframe(
      nhc_placebo_model$obs_selection,
      name = "drop_step",
      value = "row_indices"
    ) |>
      mutate(rows_dropped = lengths(row_indices)) |>
      select(drop_step, rows_dropped)

    # Store model-level statistics that are useful for manuscript tables and
    # appendix reporting.
    placebo_stats_results <- tibble(
      crime_group = crime_group_name,
      observations_original = nhc_placebo_model$nobs_origin,
      observations_dropped = nhc_placebo_model$nobs_origin - nhc_placebo_model$nobs,
      observations_remaining = nhc_placebo_model$nobs,
      dropped_row_summary = list(nhc_placebo_dropped_row_summary),
      dropped_fixed_effects = list(nhc_placebo_dropped_fixed_effects),
      dropped_terms = list(nhc_placebo_model$collin.var),
      hexagon_fixed_effects = nhc_placebo_model$fixef_sizes[["hex_id"]],
      date_fixed_effects = nhc_placebo_model$fixef_sizes[["crime_date"]],
      log_likelihood = unname(nhc_placebo_fit_stats$ll),
      aic = unname(nhc_placebo_fit_stats$aic),
      bic = unname(nhc_placebo_fit_stats$bic),
      rmse = unname(nhc_placebo_fit_stats$rmse),
      squared_correlation = unname(nhc_placebo_fit_stats$sq.cor),
      pseudo_r2 = unname(nhc_placebo_fit_stats$pr2),
      adjusted_pseudo_r2 = unname(nhc_placebo_fit_stats$apr2),
      joint_band_test_statistic = unname(nhc_placebo_band_test$stat),
      joint_band_test_p_value = unname(nhc_placebo_band_test$p),
      joint_band_test_df1 = unname(nhc_placebo_band_test$df1),
      joint_band_test_df2 = unname(nhc_placebo_band_test$df2),
      joint_band_test_vcov = nhc_placebo_band_test$vcov
    )

    # Save the compact placebo model statistics for the current crime group.
    write_rds(
      placebo_stats_results,
      here("derived_data", str_glue("nhc_stats_placebo_{crime_group_name}.rds")),
      compress = "gz"
    )

    # Keep the placebo Carnival-day rows for the current crime group.
    nhc_placebo_day_data <- nhc_panel |>
      filter(crime_group == crime_group_name, is_nhc_placebo == 1)

    # Convert estimates across all placebo days in the data to estimates for a
    # single annual two-day placebo occurrence.
    placebo_annual_scale <- 2 / n_distinct(nhc_placebo_day_data$crime_date)

    # Create a counterfactual version of the placebo-day rows where the placebo
    # event did not occur on the same dates in the same hexagons.
    no_nhc_placebo_day_data <- nhc_placebo_day_data |>
      mutate(is_nhc_placebo = 0)

    # Store the model coefficients for use in the placebo extra-crime
    # calculations.
    nhc_placebo_model_coefs <- coef(nhc_placebo_model)

    # Store the model variance-covariance matrix so the London-wide interval can
    # account for covariance between distance-band coefficients.
    nhc_placebo_model_vcov <- vcov(nhc_placebo_model)

    # Estimate the additional crimes predicted on placebo days by distance band.
    placebo_extra_results <- nhc_placebo_day_data |>
      mutate(
        predicted_nhc = predict(
          nhc_placebo_model,
          newdata = nhc_placebo_day_data,
          type = "response"
        ),
        predicted_no_nhc = predict(
          nhc_placebo_model,
          newdata = no_nhc_placebo_day_data,
          type = "response"
        ),
        # Treat predictions for model-dropped all-zero fixed effects as zero.
        predicted_nhc = replace_na(predicted_nhc, 0),
        predicted_no_nhc = replace_na(predicted_no_nhc, 0)
      ) |>
      summarise(
        observed_nhc = sum(crime_count) * placebo_annual_scale,
        predicted_nhc = sum(predicted_nhc) * placebo_annual_scale,
        predicted_no_nhc = sum(predicted_no_nhc) * placebo_annual_scale,
        .by = c(dist_km, dist)
      ) |>
      # Link each distance band to the corresponding placebo model coefficient.
      mutate(
        model_term = if_else(
          dist == "12km",
          NA_character_,
          str_glue("dist::{dist}:is_nhc_placebo")
        ),
        nhc_effect = if_else(
          is.na(model_term),
          0,
          nhc_placebo_model_coefs[model_term]
        ),
        nhc_effect_se = if_else(
          is.na(model_term),
          0,
          sqrt(diag(nhc_placebo_model_vcov)[model_term])
        ),
        extra_crimes = predicted_no_nhc * (exp(nhc_effect) - 1),
        extra_crimes_se = predicted_no_nhc * exp(nhc_effect) * nhc_effect_se,
        extra_crimes_conf_low = extra_crimes -
          qnorm(0.975) * extra_crimes_se,
        extra_crimes_conf_high = extra_crimes +
          qnorm(0.975) * extra_crimes_se
      )

    # Store the model terms that contribute uncertainty to the London-wide total.
    total_model_terms <- placebo_extra_results |>
      filter(!is.na(model_term)) |>
      pull(model_term)

    # Store the derivative of the London-wide placebo extra-crime estimate with
    # respect to each distance-band coefficient.
    total_extra_crimes_gradient <- placebo_extra_results |>
      filter(!is.na(model_term)) |>
      mutate(extra_crimes_gradient = predicted_no_nhc * exp(nhc_effect)) |>
      pull(extra_crimes_gradient)

    # Calculate the standard error for the London-wide total using the full
    # covariance matrix for the distance-band coefficients.
    total_extra_crimes_se <- sqrt(as.numeric(
      t(total_extra_crimes_gradient) %*%
        nhc_placebo_model_vcov[total_model_terms, total_model_terms, drop = FALSE] %*%
        total_extra_crimes_gradient
    ))

    # Add the crime-group label and the London-wide total row to the placebo
    # results.
    placebo_extra_results <- bind_rows(
      placebo_extra_results |>
        mutate(
          crime_group = crime_group_name,
          dist = as.character(dist),
          .before = 1
        ),
      placebo_extra_results |>
        summarise(
          crime_group = crime_group_name,
          dist_km = NA_real_,
          dist = "all_london",
          extra_crimes = sum(extra_crimes),
          extra_crimes_conf_low = extra_crimes -
            qnorm(0.975) * total_extra_crimes_se,
          extra_crimes_conf_high = extra_crimes +
            qnorm(0.975) * total_extra_crimes_se,
          observed_nhc = sum(observed_nhc),
          predicted_nhc = sum(predicted_nhc),
          predicted_no_nhc = sum(predicted_no_nhc)
        )
    ) |>
      # Keep only the fields needed in the saved placebo extra-crime results.
      select(
        crime_group,
        dist_km,
        dist,
        extra_crimes,
        extra_crimes_conf_low,
        extra_crimes_conf_high,
        observed_nhc,
        predicted_nhc,
        predicted_no_nhc
      ) |>
      # Order distance bands from the footprint out, with the total row last.
      arrange(dist_km)

    # Save the placebo extra-crime estimates for the current crime group.
    write_rds(
      placebo_extra_results,
      here("derived_data", str_glue("nhc_extra_placebo_{crime_group_name}.rds")),
      compress = "gz"
    )
  })
