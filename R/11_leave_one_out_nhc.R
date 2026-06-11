# This code runs leave-one-carnival-year-out robustness checks for the Notting
# Hill Carnival models. It only re-estimates crime groups whose main model found
# a statistically significant Carnival effect inside the operational footprint.

# This code might need a larger memory allocation.
mem.maxVSize(mem.maxVSize() * 4)

# Load the packages needed to estimate fixed-effects count models, extract
# coefficients and handle the panel data.
pacman::p_load(broom, fixest, furrr, future, here, janitor, parallelly, tidyverse)

# LOAD PANEL DATA -----------------------------------------------------------

# Load the daily hexagon-level panel used by the main NHC models.
nhc_panel <- read_rds(here("derived_data", "nhc_panel_350m_daily.rds"))

# Store the available crime groups from the current panel definition.
nhc_crime_groups <- nhc_panel |>
  distinct(crime_group) |>
  arrange(crime_group) |>
  pull(crime_group) |>
  as.character()

# FILTER MODEL DATA ---------------------------------------------------------

# Keep the requested crime group, limiting sexual-offence robustness checks to
# the years with usable geographic co-ordinates.
filter_nhc_model_panel <- function(panel_data, crime_type) {
  panel_data |>
    filter(
      crime_group == crime_type,
      crime_type != "sexual_offences" |
        between(year(crime_date), 2013, 2019)
    )
}

# IDENTIFY CRIME GROUPS TO CHECK -------------------------------------------

# Keep only the main combined-Carnival coefficient files, excluding the
# day-specific and placebo coefficient files.
nhc_coef_files <- here("derived_data", str_glue("nhc_coef_{nhc_crime_groups}.rds"))
names(nhc_coef_files) <- nhc_crime_groups
nhc_coef_files <- nhc_coef_files[file.exists(nhc_coef_files)]

# Select crime groups whose main model found a significant Carnival effect in
# the footprint band. This is the pre-specified subset for the leave-one-out
# robustness check.
nhc_leave_one_out_crime_groups <- imap_dfr(
  nhc_coef_files,
  \(coef_path, crime_group_name) {
    read_rds(coef_path) |>
      filter(term == "dist::0km:is_nhc") |>
      transmute(
        crime_group = crime_group_name,
        footprint_estimate = estimate,
        footprint_std_error = std.error,
        footprint_p_value = p.value,
        footprint_conf_low = conf.low,
        footprint_conf_high = conf.high,
        run_leave_one_out = p.value < 0.05
      )
  }
) |>
  filter(run_leave_one_out) |>
  arrange(crime_group)

write_rds(
  nhc_leave_one_out_crime_groups,
  here("derived_data", "nhc_leave_one_out_crime_groups.rds"),
  compress = "gz"
)

# Stop early with an informative message if none of the main footprint effects
# meet the significance criterion.
if (nrow(nhc_leave_one_out_crime_groups) == 0) {
  message("No crime groups had a significant 0 km Carnival effect; no leave-one-out models run.")
  quit(save = "no", status = 0)
}

# Use a conservative number of workers because each fixed-effects model can use
# a large amount of memory while it is being estimated and saved.
nhc_leave_one_out_workers <- min(3L, max(1L, as.integer(availableCores()) - 1L))

# Allow the parallel backend to register the large in-memory panel as a global.
options(future.globals.maxSize = 16 * 1024^3)

# Save leave-one-out model objects in their own subdirectory while results are
# being extracted. The model files are deleted after the summary files have been
# created so the directory should normally be empty after a successful run.
nhc_leave_one_out_model_dir <- here("derived_data", "nhc_leave_one_out_models")
dir.create(nhc_leave_one_out_model_dir, recursive = TRUE, showWarnings = FALSE)

# MODEL FUNCTION ------------------------------------------------------------

# Estimate the main NHC model for one crime group while omitting all rows from
# one Carnival year. This removes both Carnival days and comparison days from
# that calendar year, so each refit tests whether the main estimate is sensitive
# to the presence of that year's observations.
fit_nhc_leave_one_out_model <- function(crime_type, omitted_year) {
  # Store the path for the omitted-year model before estimating so existing
  # model files can be reused.
  model_path <- file.path(
    nhc_leave_one_out_model_dir,
    str_glue("nhc_model_leave_one_out_{crime_type}_omit_{omitted_year}.rds")
  )

  # Load the omitted-year model if it is already available on disk.
  if (file.exists(model_path)) {
    saved_model <- read_rds(model_path)
  } else {
    # Fit and save the omitted-year model only if it is not already available on
    # disk.
    model_data <- nhc_panel |>
      filter_nhc_model_panel(crime_type) |>
      filter(year(crime_date) != omitted_year)

    model <- fepois(
      crime_count ~ i(dist, is_nhc, ref = "12km") | hex_id + crime_date,
      data = model_data,
      cluster = ~hex_id
    )

    write_rds(model, model_path, compress = "gz")
    rm(model)
    saved_model <- read_rds(model_path)
  }

  coef_results <- tidy(
    saved_model,
    conf.int = TRUE,
    conf.level = 0.95
  ) |>
    mutate(
      crime_group = crime_type,
      omitted_year = omitted_year,
      model_path = model_path,
      .before = 1
    )

  nhc_fit_stats <- fitstat(
    saved_model,
    c("ll", "aic", "bic", "rmse", "sq.cor", "pr2", "apr2")
  )

  nhc_band_test <- wald(saved_model, keep = "is_nhc")

  stats_results <- tibble(
    crime_group = crime_type,
    omitted_year = omitted_year,
    model_path = model_path,
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
    joint_band_test_statistic = unname(nhc_band_test$stat),
    joint_band_test_p_value = unname(nhc_band_test$p),
    joint_band_test_df1 = unname(nhc_band_test$df1),
    joint_band_test_df2 = unname(nhc_band_test$df2),
    joint_band_test_vcov = nhc_band_test$vcov
  )

  # Keep compact term-level and model-level summaries. The temporary model file
  # is deleted only after the combined summary files have been written.
  tibble(
    crime_group = crime_type,
    omitted_year = omitted_year,
    coefficients = list(coef_results),
    model_stats = list(stats_results)
  )
}

# RUN LEAVE-ONE-OUT MODELS --------------------------------------------------

# Create a crime-specific omitted-year grid so sexual-offence models only omit
# Carnival years from 2013 to 2019.
nhc_leave_one_out_grid <- nhc_leave_one_out_crime_groups |>
  select(crime_group) |>
  mutate(
    omitted_year = map(
      crime_group,
      \(crime_group) {
        filter_nhc_model_panel(nhc_panel, crime_group) |>
          filter(is_nhc == 1) |>
          distinct(omitted_year = year(crime_date)) |>
          arrange(omitted_year) |>
          pull(omitted_year)
      }
    )
  ) |>
  unnest(omitted_year)

if (supportsMulticore()) {
  plan(multicore, workers = nhc_leave_one_out_workers)
} else {
  plan(multisession, workers = nhc_leave_one_out_workers)
}

nhc_leave_one_out_results <- nhc_leave_one_out_grid |>
  future_pmap(
    \(crime_group, omitted_year) {
      fit_nhc_leave_one_out_model(crime_group, omitted_year)
    },
    .options = furrr_options(seed = NULL)
  ) |>
  list_rbind()

plan(sequential)

# EXTRACT AND SAVE RESULTS --------------------------------------------------

nhc_leave_one_out_coef <- nhc_leave_one_out_results |>
  select(coefficients) |>
  unnest(coefficients) |>
  arrange(crime_group, omitted_year, term)

nhc_leave_one_out_stats <- nhc_leave_one_out_results |>
  select(model_stats) |>
  unnest(model_stats) |>
  arrange(crime_group, omitted_year)

# Keep a reporting-ready table for the footprint term that drove selection into
# the robustness check.
nhc_leave_one_out_footprint <- nhc_leave_one_out_coef |>
  filter(term == "dist::0km:is_nhc") |>
  left_join(
    nhc_leave_one_out_crime_groups |>
      select(
        crime_group,
        main_estimate = footprint_estimate,
        main_std_error = footprint_std_error,
        main_p_value = footprint_p_value,
        main_conf_low = footprint_conf_low,
        main_conf_high = footprint_conf_high
      ),
    by = join_by(crime_group)
  ) |>
  mutate(
    exponentiated_estimate = exp(estimate),
    exponentiated_conf_low = exp(conf.low),
    exponentiated_conf_high = exp(conf.high),
    main_exponentiated_estimate = exp(main_estimate),
    estimate_change_from_main = estimate - main_estimate,
    percent_change_from_main = 100 * (estimate - main_estimate) / abs(main_estimate),
    significant = p.value < 0.05,
    same_sign_as_main = sign(estimate) == sign(main_estimate)
  ) |>
  select(
    crime_group,
    omitted_year,
    term,
    estimate,
    std.error,
    statistic,
    p.value,
    conf.low,
    conf.high,
    exponentiated_estimate,
    exponentiated_conf_low,
    exponentiated_conf_high,
    main_estimate,
    main_std_error,
    main_p_value,
    main_conf_low,
    main_conf_high,
    main_exponentiated_estimate,
    estimate_change_from_main,
    percent_change_from_main,
    significant,
    same_sign_as_main
  ) |>
  arrange(crime_group, omitted_year)

write_rds(
  nhc_leave_one_out_coef,
  here("derived_data", "nhc_coef_leave_one_out.rds"),
  compress = "gz"
)

write_rds(
  nhc_leave_one_out_stats,
  here("derived_data", "nhc_stats_leave_one_out.rds"),
  compress = "gz"
)

write_rds(
  nhc_leave_one_out_footprint,
  here("derived_data", "nhc_footprint_leave_one_out.rds"),
  compress = "gz"
)

# Remove the saved leave-one-out model objects after all summary files have been
# extracted and written. The subdirectory is retained as an audit trail for where
# temporary model files were stored during the run.
nhc_leave_one_out_model_files <- list.files(
  nhc_leave_one_out_model_dir,
  pattern = "^nhc_model_leave_one_out_.*\\.rds$",
  full.names = TRUE
)
unlink(nhc_leave_one_out_model_files)
