# This code runs the Notting Hill Carnival crime-count models.

# This code might need a larger memory allocation
mem.maxVSize(mem.maxVSize() * 4)

# Load the packages needed to estimate fixed-effects count models and handle the
# panel data.
pacman::p_load(fixest, furrr, future, here, janitor, parallelly, tidyverse)

# LOAD PANEL DATA -----------------------------------------------------------

# Load the daily hexagon-level panel and standardise column names before any
# modelling operations.
nhc_panel <- read_rds(here("derived_data", "nhc_panel_350m_daily.rds"))

# Store the crime groups in the panel so the model script always matches the
# latest panel definition.
nhc_crime_groups <- nhc_panel |>
  distinct(crime_group) |>
  arrange(crime_group) |>
  pull(crime_group) |>
  as.character()

# Use a conservative number of workers because each fixed-effects model can use
# a large amount of memory while it is being estimated and saved.
nhc_model_workers <- min(3L, max(1L, as.integer(availableCores()) - 1L))

# Allow the parallel backend to register the large in-memory panel as a global.
# With forked workers, this avoids re-reading the panel for every model while
# still keeping the worker count conservative.
options(future.globals.maxSize = 16 * 1024^3)

# MODEL FUNCTION ------------------------------------------------------------

# Estimate the same fixed-effects Poisson model for one crime type at a time.
fit_nhc_model <- function(crime_type) {
  # Keep the rows for the requested crime type so each model has the same
  # specification but a different outcome subset.
  model <- fepois(
    crime_count ~ i(dist, is_nhc, ref = "12km") | hex_id + crime_date,
    data = filter(nhc_panel, crime_group == crime_type),
    cluster = ~hex_id
  )

  # Store the path for the saved model so the same value can be used for
  # writing the file and recording the run summary.
  model_path <- here("derived_data", str_glue("nhc_model_{crime_type}.rds"))

  # Save the fitted model so it can be reused without re-estimating.
  write_rds(model, model_path, compress = "gz")

  # Clean up memory used by the model fit inside the worker process.
  gc()

  # Return a lightweight summary so the parent process does not have to hold all
  # fitted model objects in memory at the same time.
  tibble(crime_group = crime_type, model_path = model_path)
}

# RUN MODELS ---------------------------------------------------------------

# Fit multiple crime-group models at the same time using a cautious worker
# count to reduce run time without overwhelming system memory.
if (supportsMulticore()) {
  # Use forked workers when available so the large panel can be shared more
  # efficiently between worker processes.
  plan(multicore, workers = nhc_model_workers)
} else {
  # Use socket workers when forked workers are not supported by the local R
  # session.
  plan(multisession, workers = nhc_model_workers)
}

# Estimate and save one fixed-effects Poisson model for each crime group in the
# panel dataset.
nhc_model_files <- nhc_crime_groups |>
  future_map(
    fit_nhc_model,
    .options = furrr_options(seed = NULL)
  ) |>
  list_rbind()

# Return the future backend to sequential processing after the model run.
plan(sequential)

# MODEL FUNCTION FOR SEPARATE CARNIVAL DAYS --------------------------------

# Estimate a fixed-effects Poisson model with separate distance-band effects
# for Carnival Sunday and Carnival Monday.
fit_nhc_day_model <- function(crime_type) {
  # Keep the rows for the requested crime type so the Sunday/Monday model uses
  # the same outcome subset as the combined Carnival-day model.
  model <- fepois(
    crime_count ~
      i(dist, nhc_sunday, ref = "12km") +
      i(dist, nhc_monday, ref = "12km") |
      hex_id + crime_date,
    data = filter(nhc_panel, crime_group == crime_type),
    cluster = ~hex_id
  )

  # Store the path for the saved day-specific model so it follows the project
  # naming convention for reusable model outputs.
  model_path <- here("derived_data", str_glue("nhc_model_day_{crime_type}.rds"))

  # Save the fitted model using compression because the fixed-effects model
  # objects are large.
  write_rds(model, model_path, compress = "gz")

  # Clean up memory used by the day-specific model fit inside the worker
  # process.
  gc()

  # Return a lightweight summary so the parent process does not collect all
  # fitted day-specific model objects in memory.
  tibble(crime_group = crime_type, model_path = model_path)
}

# RUN SEPARATE CARNIVAL-DAY MODELS -----------------------------------------

# Fit multiple Sunday/Monday models at the same time using the same cautious
# worker-count strategy as the combined Carnival-day models.
if (supportsMulticore()) {
  # Use forked workers when available so the large panel can be shared more
  # efficiently between worker processes.
  plan(multicore, workers = nhc_model_workers)
} else {
  # Use socket workers when forked workers are not supported by the local R
  # session.
  plan(multisession, workers = nhc_model_workers)
}

# Estimate and save one Sunday/Monday fixed-effects Poisson model for each
# crime group in the panel dataset.
nhc_day_model_files <- nhc_crime_groups |>
  future_map(
    fit_nhc_day_model,
    .options = furrr_options(seed = NULL)
  ) |>
  list_rbind()

# Return the future backend to sequential processing after the day-specific
# model run.
plan(sequential)

# PLACEBO CARNIVAL DATES ----------------------------------------------------

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

# Add placebo Carnival indicators to the panel for the last Sunday in July and
# the Monday immediately following it.
nhc_panel <- nhc_panel |>
  mutate(
    nhc_placebo_sunday = as.numeric(
      crime_date %in% nhc_placebo_dates$placebo_sunday_date
    ),
    nhc_placebo_monday = as.numeric(
      crime_date %in% nhc_placebo_dates$placebo_monday_date
    ),
    is_nhc_placebo = as.numeric(
      nhc_placebo_sunday == 1 | nhc_placebo_monday == 1
    )
  )

# MODEL FUNCTION FOR PLACEBO CARNIVAL DAYS ---------------------------------

# Estimate a fixed-effects Poisson model using the placebo Carnival date pair
# for one crime type at a time.
fit_nhc_placebo_model <- function(crime_type) {
  # Keep the rows for the requested crime type so the placebo model uses the
  # same outcome subset as the main Carnival-day model.
  model <- fepois(
    crime_count ~ i(dist, is_nhc_placebo, ref = "12km") | hex_id + crime_date,
    data = filter(nhc_panel, crime_group == crime_type),
    cluster = ~hex_id
  )

  # Store the path for the saved placebo model so it follows the project naming
  # convention for reusable model outputs.
  model_path <- here(
    "derived_data",
    str_glue("nhc_model_placebo_{crime_type}.rds")
  )

  # Save the fitted placebo model using compression because the fixed-effects
  # model objects are large.
  write_rds(model, model_path, compress = "gz")

  # Clean up memory used by the placebo model fit inside the worker process.
  gc()

  # Return a lightweight summary so the parent process does not collect all
  # fitted placebo model objects in memory.
  tibble(crime_group = crime_type, model_path = model_path)
}

# RUN PLACEBO CARNIVAL-DAY MODELS ------------------------------------------

# Fit multiple placebo models at the same time using the same cautious
# worker-count strategy as the main Carnival-day models.
if (supportsMulticore()) {
  # Use forked workers when available so the large panel can be shared more
  # efficiently between worker processes.
  plan(multicore, workers = nhc_model_workers)
} else {
  # Use socket workers when forked workers are not supported by the local R
  # session.
  plan(multisession, workers = nhc_model_workers)
}

# Estimate and save one placebo fixed-effects Poisson model for each crime group
# in the panel dataset.
nhc_placebo_model_files <- nhc_crime_groups |>
  future_map(
    fit_nhc_placebo_model,
    .options = furrr_options(seed = NULL)
  ) |>
  list_rbind()

# Return the future backend to sequential processing after the placebo model
# run.
plan(sequential)
