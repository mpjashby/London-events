# This code loads the cleaned crime data and builds panel datasets for the
# different event studies used in the project.

# This code might need a larger memory allocation
mem.maxVSize(mem.maxVSize() * 4)

# Load the packages used throughout the panel-building workflow.
pacman::p_load(
  here,
  sf,
  janitor,
  lubridate,
  scales,
  timeDate,
  tsibble,
  tidyverse
)

# LOAD SHARED INPUTS --------------------------------------------------------

# Load the cleaned crime data and standardise column names before any other
# operations so that the rest of the script can rely on lower-case names.
crimes <- read_rds(here("derived_data", "crimes.rds")) |>
  # Standardise column names before filtering or converting the data.
  clean_names() |>
  # Convert the point-level crime records to an sf object in the British
  # National Grid so that the geometry aligns with the 250-metre grid definition
  filter(!is.na(easting), !is.na(northing)) |>
  st_as_sf(coords = c("easting", "northing"), crs = "EPSG:27700")

# Load dataset containing London boundary
if (!file.exists(here("derived_data/london_boundary.gpkg"))) {
  download.file(
    url = "https://data.london.gov.uk/download/20od9/114d1137-e339-4b50-b409-124c17f4b59a/gla.zip",
    destfile = london_boundary_file <- tempfile(fileext = ".zip")
  )
  unzip(london_boundary_file, exdir = tempdir())
  str_glue("{tempdir()}/gla/London_GLA_Boundary.dbf") |>
    # Standardise column names before saving the reusable boundary dataset.
    read_sf() |>
    clean_names() |>
    write_sf(here("derived_data/london_boundary.gpkg"))
}
london_boundary <- here("derived_data/london_boundary.gpkg") |>
  # Standardise column names before the boundary is used to make spatial data.
  read_sf() |>
  clean_names() |>
  select()


# NOTTING HILL CARNIVAL -----------------------------------------------------

# Load the operational carnival footprint and dissolve it to a single geometry
# so that distances are measured from the edge of the footprint area.
carnival_footprint <- here("original_data", "NHC operational footprint.kml") |>
  read_sf() |>
  clean_names() |>
  st_transform("EPSG:27700") |>
  st_union()

# Create a hexagonal grid covering the London crime extent.
nhc_grid <- st_make_grid(london_boundary, cellsize = 350, square = FALSE) |>
  st_as_sf() |>
  st_set_geometry("geometry") |>
  st_intersection(london_boundary) |>
  mutate(
    hex_id = row_number(),
    # Calculate each grid-cell centroid so that the footprint band is based on
    # whether the cell centre falls inside the operational footprint.
    centroid = st_centroid(geometry),
    # Identify cells whose centroid falls inside the carnival footprint.
    centroid_in_nhc = lengths(st_within(centroid, carnival_footprint)) > 0,
    # Measure the distance from each grid-cell centroid to the carnival
    # footprint.
    dist_to_nhc = as.numeric(st_distance(centroid, carnival_footprint)),
    # Store the distance band as the upper edge of each one-kilometre band for
    # plotting and simple trend checks.
    dist_km = if_else(centroid_in_nhc, 0, ceiling(dist_to_nhc / 1000)),
    # Cap all outer London cells at 12 km so that the sparse outer bands remain
    # a single comparison category.
    dist_km = oob_squish(dist_km, range = c(0, 12)),
    # Store the modelling distance band as an unordered factor with the furthest
    # distance band first so that it is the model reference category.
    dist = factor(
      str_glue("{dist_km}km"),
      levels = str_glue("{seq(12, 0, by = -1)}km")
    )
  ) |>
  select(-centroid, -centroid_in_nhc, -dist_to_nhc) |>
  write_sf(here("derived_data", "nhc_hex_grid_350m.gpkg"))

# Keep the crime categories used in the NHC analysis and harmonise them into
# the analysis groups used for separate models.
nhc_crime_counts <- crimes |>
  mutate(
    # Shift overnight crimes to the previous analysis day so that daily counts
    # run from 05:00 on each date to 04:59 on the following date.
    crime_date = as_date(date_committed_from - hours(5)),
    crime_group = case_when(
      new_minor_text %in% c("Violence with Injury", "Homicide") ~
        "violence_injury_homicide",
      new_minor_text == "Violence without Injury" ~
        "violence_no_injury",
      new_minor_text == "Robbery of Personal Property" ~
        "personal_robbery",
      new_minor_text == "Theft from Person" ~
        "theft_from_person",
      new_minor_text == "Bicycle Theft" ~
        "bicycle_theft",
      new_major_text == "Arson and Criminal Damage" ~ "criminal_damage",
      new_minor_text %in%
        c("Burglary - Residential", "Domestic Burglary") ~
        "residential_burglary",
      new_minor_text %in%
        c(
          "Burglary - Business and Community",
          "Burglary Business and Community"
        ) ~
        "commercial_business_burglary",
      new_major_text == "Sexual Offences" ~ "sexual_offences",
      new_major_text == "Vehicle Offences" ~ "vehicle_theft",
      TRUE ~ NA_character_
    )
  ) |>
  filter(
    !is.na(crime_date),
    # Keep only the crimes we are interested in
    !is.na(crime_group),
    # Keep only those crimes in July, August, September and October to avoid the
    # dataset becoming too large
    month(crime_date) %in% 7:10
  ) |>
  # Join crimes to grid cells so that each offence is assigned to one hexagon.
  st_join(nhc_grid) |>
  # Count the number of crimes in each cell on each day for the crime groups
  # used in the NHC analysis.
  st_drop_geometry() |>
  count(hex_id, crime_date, crime_group, name = "crime_count")

# Create a dataset specifying which days are carnival days.
nhc_dates <- holidayLONDON(year = 2013:2023) |>
  as_date() |>
  enframe(name = NULL, value = "date") |>
  mutate(is_holiday = TRUE) |>
  as_tsibble(index = date) |>
  fill_gaps(is_holiday = FALSE) |>
  mutate(
    # Identify the late-August bank-holiday Monday on which Carnival takes
    # place in normal years.
    bank_holiday_monday = is_holiday & month(date) == 8,
    # Identify the Sunday immediately before the late-August bank-holiday
    # Monday.
    nhc_sunday = lead(bank_holiday_monday),
    # Identify the late-August bank-holiday Monday itself.
    nhc_monday = bank_holiday_monday,
    # Exclude the years when Carnival did not occur in its usual form because
    # of the COVID-19 pandemic.
    across(
      c(nhc_sunday, nhc_monday),
      ~ .x & !year(date) %in% 2020:2021
    ),
    # Keep a combined Carnival-day indicator for the existing main model.
    is_nhc = nhc_sunday | nhc_monday
  ) |>
  as_tibble() |>
  select(date, nhc_sunday, nhc_monday, is_nhc) |>
  replace_na(list(nhc_sunday = FALSE, nhc_monday = FALSE, is_nhc = FALSE))

# Create a full hexagon-by-day panel so that dates with zero crimes are kept in
# the analysis dataset.
nhc_panel <- expand_grid(
  hex_id = pull(nhc_grid, hex_id),
  crime_date = seq(
    ymd("2013-07-01"),
    ymd("2023-10-31"),
    by = "day"
  ),
  crime_group = unique(pull(nhc_crime_counts, crime_group))
) |>
  filter(
    month(crime_date) %in% 7:10,
    # Data for September 2019 is missing, so remove the corresponding incorrect
    # counts from the panel
    !(year(crime_date) == 2019 & month(crime_date) == 9)
  ) |>
  # Join the crime counts
  left_join(nhc_crime_counts, by = join_by(hex_id, crime_date, crime_group)) |>
  # Replace missing counts (representing days with zero crimes in that cell)
  # with zeros
  replace_na(list(crime_count = 0)) |>
  # Join cell IDs
  left_join(st_drop_geometry(nhc_grid), by = join_by(hex_id)) |>
  # Join carnival dates
  left_join(nhc_dates, by = join_by(crime_date == date)) |>
  # Convert variables to factors to reduce dataset storage size and convert
  # logical Carnival-day flags to numeric indicators for modelling.
  mutate(
    across(where(is.character), factor),
    across(c(nhc_sunday, nhc_monday, is_nhc), as.numeric)
  )

# Write panel data to file
write_rds(
  nhc_panel,
  here("derived_data", "nhc_panel_350m_daily.rds"),
  compress = "gz"
)

# FOOTBALL ------------------------------------------------------------------

# Add the code for football panel construction in this section.

# OTHER LONDON EVENTS -------------------------------------------------------

# Add the code for other London event panels in this section.
