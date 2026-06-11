# This code extracts offence-method text for criminal-damage offences inside the
# Notting Hill Carnival operational footprint.

# Load the packages used to read data, handle dates, process spatial objects and
# manipulate the offence records.
pacman::p_load(
  here,
  janitor,
  lubridate,
  sf,
  tidyverse
)

# LOAD SHARED INPUTS --------------------------------------------------------

# Load the cleaned crime records created by the crime-build script.
crimes <- read_rds(here("derived_data", "crimes.rds"))

# Load the daily Notting Hill Carnival panel so the event dates match the
# event-study models.
nhc_panel <- read_rds(here("derived_data", "nhc_panel_350m_daily.rds"))

# Load the operational carnival footprint and dissolve it to a single geometry
# in the same British National Grid coordinate system as the crime records.
carnival_footprint <- here("original_data", "NHC operational footprint.kml") |>
  read_sf() |>
  clean_names() |>
  st_transform("EPSG:27700") |>
  st_union()

# DEFINE EXTENDED EVENT WINDOW ---------------------------------------------

# Store the extended event-study labels from 10 days before Carnival Monday to
# three days after Carnival Monday, matching the criminal-damage event study.
nhc_criminal_damage_extended_event_days <- tibble(
  event_day_offset = -10:3,
  nhc_event_day = c(
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

# Create one row for every date in the extended event-study window around each
# real Carnival Monday used in the panel.
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
    nhc_event_day = nhc_event_day
  )

# EXTRACT CRIMINAL-DAMAGE METHODS ------------------------------------------

# Keep criminal-damage offences with offence-method text, offence dates and
# valid coordinates, then convert them to points for filtering to the
# operational footprint.
nhc_criminal_damage_methods <- crimes |>
  mutate(
    crime_date = as_date(date_committed_from - hours(5)),
    method = na_if(method, "")
  ) |>
  filter(
    new_major_text == "Arson and Criminal Damage",
    !is.na(crime_date),
    !is.na(method),
    !is.na(easting),
    !is.na(northing)
  ) |>
  st_as_sf(coords = c("easting", "northing"), crs = "EPSG:27700") |>
  filter(lengths(st_within(geometry, carnival_footprint)) > 0) |>
  st_drop_geometry() |>
  left_join(
    nhc_criminal_damage_extended_window_dates,
    by = join_by(crime_date)
  ) |>
  mutate(
    nhc_event_day = replace_na(nhc_event_day, "outside_window"),
    nhc_event_day = factor(
      nhc_event_day,
      levels = c(
        "outside_window",
        nhc_criminal_damage_extended_event_days$nhc_event_day
      )
    )
  ) |>
  select(crime_date, nhc_event_day, method) |>
  arrange(crime_date, nhc_event_day, method)

# Save the compact offence-method extract for the criminal-damage method
# analysis.
write_rds(
  nhc_criminal_damage_methods,
  here("derived_data", "nhc_criminal_damage_methods.rds"),
  compress = "gz"
)
