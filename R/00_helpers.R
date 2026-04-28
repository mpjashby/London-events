#!/usr/bin/env Rscript

# Stadium-centred spatio-temporal quasi-experiment helpers
# --------------------------------------------------------
# This file contains reusable functions shared by the fixtures, panel, and
# modelling scripts. The workflow is designed for large London crime datasets
# and emphasises explicit timezone handling, modularity, and inspectable
# intermediate objects.

pacman::p_load(fixest, here, readxl, sf, tidyverse)

# Keep these constants in one place so all scripts use the same timezone and
# coordinate systems when defining treatment timing and spatial exposure.
tz_london <- "Europe/London"
crs_wgs84 <- "EPSG:4326"
crs_bng <- "EPSG:27700"

# Return the first value unless it is NULL, empty, or all missing.
coalesce_or <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

# Create a directory if it does not already exist.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# Convert values to lower-case strings while preserving missing values.
safe_lower <- function(x) {
  str_to_lower(as.character(coalesce_or(x, NA_character_)))
}

# Parse mixed-format datetimes and return them in the London timezone.
as_london_datetime <- function(x, tz = tz_london) {
  # Inputs can arrive as text or POSIXct objects depending on source.
  # Standardising them here avoids inconsistent timezone handling later on.
  if (inherits(x, "POSIXt")) {
    return(with_tz(x, tzone = tz))
  }

  x_chr <- as.character(x)
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_

  parsed <- suppressWarnings(
    parse_date_time(
      x_chr,
      orders = c(
        "Ymd HMS", "Ymd HM", "YmdHMS", "YmdHM",
        "dmY HMS", "dmY HM", "dmYHMS", "dmYHM",
        "dmy HMS", "dmy HM", "dmyHMS", "dmyHM",
        "mdY HMS", "mdY HM", "mdYHMS", "mdYHM",
        "Y-m-d H:M:S", "Y-m-d H:M",
        "d/m/Y H:M:S", "d/m/Y H:M"
      ),
      tz = tz,
      exact = FALSE
    )
  )

  parsed
}

# Split comma-separated British National Grid references into coordinates.
parse_grid_ref <- function(x) {
  x_chr <- as.character(x)
  x_chr[x_chr %in% c("", "NA")] <- NA_character_
  mat <- str_match(x_chr, "^\\s*([0-9.]+)\\s*,\\s*([0-9.]+)\\s*$")
  tibble(
    easting = suppressWarnings(as.numeric(mat[, 2])),
    northing = suppressWarnings(as.numeric(mat[, 3]))
  )
}

# Convert a date into a football season label such as 2019/20.
make_season <- function(date_value) {
  yr <- year(date_value)
  mo <- month(date_value)
  season_start <- if_else(mo >= 7, yr, yr - 1)
  paste0(season_start, "/", str_sub(season_start + 1, 3, 4))
}

# Assign distances to labelled concentric stadium bands.
label_distance_band <- function(distance_meters, breaks) {
  labels <- paste0(head(breaks, -1), "_", tail(breaks, -1), "m")
  cut(
    distance_meters,
    breaks = breaks,
    labels = labels,
    right = FALSE,
    include.lowest = TRUE
  )
}

# Standardise competition labels for filtering and modelling.
normalise_competition <- function(x) {
  safe_lower(x) |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
}

# Standardise fixture kickoff and end times in the London timezone.
standardise_kickoff_datetimes <- function(
    fixtures,
    timezone = tz_london,
    match_duration_mins = 105L
) {
  # The event-study design is anchored to kickoff and final-whistle times.
  # This helper creates those anchors in one consistent timezone.
  fixtures |>
    mutate(
      kickoff_datetime = as_london_datetime(kickoff_datetime, tz = timezone),
      kickoff_datetime = force_tz(kickoff_datetime, tzone = timezone),
      match_date = as.Date(kickoff_datetime, tz = timezone),
      match_end_datetime = case_when(
        "match_end_datetime" %in% names(.) &
          !is.na(match_end_datetime) ~ as_london_datetime(match_end_datetime, tz = timezone),
        TRUE ~ kickoff_datetime + minutes(match_duration_mins)
      ),
      match_end_datetime = force_tz(match_end_datetime, tzone = timezone),
      season = coalesce(as.character(season), make_season(match_date))
    )
}

# Load a bank-holiday calendar if provided, otherwise return an empty table.
load_bank_holidays <- function(
    bank_holiday_path = NULL,
    timezone = tz_london
) {
  if (!is.null(bank_holiday_path) && file.exists(bank_holiday_path)) {
    holidays <- read_csv(bank_holiday_path, show_col_types = FALSE)
    return(
      holidays |>
        transmute(date = as.Date(date, tz = timezone), bank_holiday = TRUE) |>
        distinct()
    )
  }

  tibble(date = as.Date(character()), bank_holiday = logical())
}

# Recode raw offence labels into the three analysis crime categories.
recode_crime_category <- function(x, offence_map) {
  stopifnot(all(c("raw_value", "crime_category") %in% names(offence_map)))

  lookup <- offence_map |>
    mutate(raw_value = safe_lower(raw_value)) |>
    distinct(raw_value, crime_category)

  tibble(raw_value = safe_lower(x)) |>
    left_join(lookup, by = "raw_value") |>
    pull(crime_category)
}

# Load raw crime data, harmonise columns, parse datetimes, and return an sf object.
load_crime_data <- function(config) {
  message("Loading crime data...")

  # Read and stack all crime files before cleaning so the same rules are
  # applied to every incident record used in the panel.
  crime_tbl <- config$crime_input_paths |>
    map(function(path) {
      ext <- tools::file_ext(path)

      if (ext %in% c("csv", "txt")) {
        read_csv(path, show_col_types = FALSE, progress = FALSE)
      } else if (ext %in% c("rds", "RDS")) {
        read_rds(path)
      } else if (ext %in% c("xlsx", "xls")) {
        sheet_to_use <- coalesce_or(config$crime_excel_sheet, 1)
        read_excel(path, sheet = sheet_to_use)
      } else {
        stop(sprintf("Unsupported crime input extension: %s", ext), call. = FALSE)
      }
    }) |>
    bind_rows()

  crime_tbl <- crime_tbl |>
    mutate(
      crime_id = as.character(.data[[config$crime_id_col]]),
      offence_raw = as.character(.data[[config$offence_col]]),
      # The analysis is hourly, so incident timestamps are floored to the
      # start of the observed hour before any aggregation.
      crime_datetime = as_london_datetime(.data[[config$datetime_col]], tz = config$timezone),
      crime_datetime = floor_date(crime_datetime, unit = "hour")
    )

  if (!is.null(config$grid_ref_col)) {
    grid_ref_xy <- parse_grid_ref(crime_tbl[[config$grid_ref_col]])
    crime_tbl <- bind_cols(crime_tbl, grid_ref_xy)
  }

  if (!is.null(config$lon_col) && !is.null(config$lat_col)) {
    crime_sf <- crime_tbl |>
      filter(!is.na(.data[[config$lon_col]]), !is.na(.data[[config$lat_col]])) |>
      st_as_sf(coords = c(config$lon_col, config$lat_col), crs = config$crime_coord_crs_lonlat)
  } else {
    crime_sf <- crime_tbl |>
      filter(!is.na(easting), !is.na(northing)) |>
      st_as_sf(coords = c("easting", "northing"), crs = config$crime_coord_crs_projected)
  }

  crime_sf <- crime_sf |>
    st_transform(config$analysis_crs) |>
    mutate(
      crime_category = recode_crime_category(offence_raw, config$offence_map),
      date = as.Date(crime_datetime, tz = config$timezone),
      hour = hour(crime_datetime),
      season = make_season(date)
    ) |>
    filter(!is.na(crime_id), !is.na(crime_datetime), !is.na(crime_category)) |>
    distinct(crime_id, .keep_all = TRUE)

  crime_sf
}

# Return the time-bounded stadium-club analysis units used to build the panel.
get_analysis_units <- function(config) {
  config$analysis_units |>
    mutate(
      unit_start_date = as.Date(unit_start_date),
      unit_end_date = as.Date(unit_end_date)
    ) |>
    left_join(
      config$stadium_venues |>
        distinct(stadium_id, stadium_name),
      by = "stadium_id"
    )
}

# Build stadium point geometries, concentric rings, and overlap flags.
build_stadium_buffers <- function(config) {
  # These geometries define the spatial treatment zones used throughout the
  # analysis: a point for each venue, an outer catchment, and concentric rings.
  stadiums_sf <- config$stadium_venues |>
    distinct(stadium_id, stadium_name, longitude, latitude) |>
    st_as_sf(coords = c("longitude", "latitude"), crs = crs_wgs84) |>
    st_transform(config$analysis_crs)

  outer_buffer <- max(config$distance_breaks)
  max_buffer <- st_buffer(stadiums_sf, dist = outer_buffer)

  bands_tbl <- tibble(
    band_start = head(config$distance_breaks, -1),
    band_end = tail(config$distance_breaks, -1),
    distance_band = paste0(head(config$distance_breaks, -1), "_", tail(config$distance_breaks, -1), "m")
  )

  band_geometries <- map2(
    bands_tbl$band_start,
    bands_tbl$band_end,
    # Each ring is the area between two buffers rather than a cumulative disk.
    ~ st_difference(st_buffer(stadiums_sf, dist = .y), st_buffer(stadiums_sf, dist = .x))
  )

  rings_sf <- map2(
    band_geometries,
    seq_len(nrow(bands_tbl)),
    function(geom, idx) {
      st_as_sf(
        stadiums_sf |>
          st_drop_geometry() |>
          mutate(
            band_start = bands_tbl$band_start[idx],
            band_end = bands_tbl$band_end[idx],
            distance_band = bands_tbl$distance_band[idx]
          ),
        geometry = st_geometry(geom),
        crs = config$analysis_crs
      )
    }
  ) |>
    bind_rows()

  overlapping_flags <- st_intersects(max_buffer)
  stadium_overlap_tbl <- tibble(
    stadium_id = stadiums_sf$stadium_id,
    overlapping_catchment = lengths(overlapping_flags) > 1
  )

  list(
    stadiums_sf = stadiums_sf,
    max_buffer_sf = max_buffer,
    rings_sf = rings_sf,
    stadium_overlap_tbl = stadium_overlap_tbl
  )
}

# Attach each crime to its nearest stadium and assign a distance band.
assign_crimes_to_stadiums <- function(crime_sf, stadium_setup, config) {
  message("Assigning crimes to stadiums and distance bands...")

  # Each crime is linked to the nearest professional venue and then placed in a
  # distance band around that venue.
  stadiums_sf <- stadium_setup$stadiums_sf

  nearest_idx <- st_nearest_feature(crime_sf, stadiums_sf)
  nearest_stadium <- stadiums_sf |>
    st_drop_geometry() |>
    slice(nearest_idx)

  distance_m <- as.numeric(
    st_distance(crime_sf, stadiums_sf[nearest_idx, ], by_element = TRUE)
  )

  assigned <- crime_sf |>
    mutate(
      stadium_id = nearest_stadium$stadium_id,
      stadium_name = nearest_stadium$stadium_name,
      distance_to_stadium_m = distance_m,
      distance_band = label_distance_band(distance_m, config$distance_breaks),
      in_analysis_catchment = distance_to_stadium_m < max(config$distance_breaks)
    ) |>
    left_join(stadium_setup$stadium_overlap_tbl, by = "stadium_id") |>
    filter(in_analysis_catchment, !is.na(distance_band))

  assigned
}

# Expand fixtures into an hourly treatment panel at the stadium level.
create_match_hour_panel <- function(fixtures, config) {
  message("Creating stadium-hour match exposure panel...")

  # Home treatment includes ordinary home fixtures plus neutral-site London
  # matches intentionally coded as treated venue events, such as Wembley finals.
  stadium_lookup <- config$stadium_venues |>
    distinct(stadium_id, stadium_name)

  treated_fixtures <- fixtures |>
    filter(home_away %in% c("home", "neutral"))

  home_hours <- treated_fixtures |>
    pmap(function(
        match_id, club, opponent, home_away, stadium_id, stadium_name,
        match_date, kickoff_datetime, match_end_datetime, competition, season, ...) {
      build_single_match_hour_rows(
        match_id = match_id,
        club = club,
        stadium_id = stadium_id,
        stadium_name = stadium_name,
        competition = competition,
        season = season,
        kickoff_datetime = kickoff_datetime,
        match_end_datetime = match_end_datetime,
        home_away = "home",
        timezone = config$timezone
      )
    }) |>
    bind_rows()

  away_hours <- fixtures |>
    filter(home_away == "away") |>
    # Away fixtures remain linked to the club's home stadium so they act as the
    # "club away from home" comparison for that same local area.
    left_join(stadium_lookup, by = "stadium_id", suffix = c("", "_lookup")) |>
    mutate(stadium_name = coalesce(stadium_name, stadium_name_lookup)) |>
    pmap(function(
        match_id, club, opponent, home_away, stadium_id, stadium_name,
        match_date, kickoff_datetime, match_end_datetime, competition, season, ...) {
      build_single_match_hour_rows(
        match_id = match_id,
        club = club,
        stadium_id = stadium_id,
        stadium_name = stadium_name,
        competition = competition,
        season = season,
        kickoff_datetime = kickoff_datetime,
        match_end_datetime = match_end_datetime,
        home_away = "away",
        timezone = config$timezone
      )
    }) |>
    bind_rows()

  bind_rows(home_hours, away_hours) |>
    group_by(stadium_id, club, hour_start) |>
    summarise(
      home_pre_4_2 = as.integer(
        any(window_name == "pre_4_2" & home_away == "home")
      ),
      home_pre_2_0 = as.integer(
        any(window_name == "pre_2_0" & home_away == "home")
      ),
      home_during_match = as.integer(
        any(window_name == "during_match" & home_away == "home")
      ),
      home_post_0_2 = as.integer(
        any(window_name == "post_0_2" & home_away == "home")
      ),
      home_post_2_4 = as.integer(
        any(window_name == "post_2_4" & home_away == "home")
      ),
      away_pre_4_2 = as.integer(
        any(window_name == "pre_4_2" & home_away == "away")
      ),
      away_pre_2_0 = as.integer(
        any(window_name == "pre_2_0" & home_away == "away")
      ),
      away_during_match = as.integer(
        any(window_name == "during_match" & home_away == "away")
      ),
      away_post_0_2 = as.integer(
        any(window_name == "post_0_2" & home_away == "away")
      ),
      away_post_2_4 = as.integer(
        any(window_name == "post_2_4" & home_away == "away")
      ),
      home_match_any = as.integer(any(home_away == "home")),
      away_match_any = as.integer(any(home_away == "away")),
      competition_home = first(
        competition[home_away == "home"],
        default = NA_character_
      ),
      competition_away = first(
        competition[home_away == "away"],
        default = NA_character_
      ),
      season = first(season, default = NA_character_),
      .groups = "drop"
    ) |>
    mutate(non_match_hour = as.integer(home_match_any == 0 & away_match_any == 0))
}

# Generate hourly rows for each named event-time window around one match.
build_single_match_hour_rows <- function(
    match_id,
    club,
    stadium_id,
    stadium_name,
    competition,
    season,
    kickoff_datetime,
    match_end_datetime,
    home_away,
    timezone = tz_london
) {
  # These windows match the exact pre/during/post structure set out in the
  # research design.
  window_starts <- c(
    pre_4_2 = kickoff_datetime - hours(4),
    pre_2_0 = kickoff_datetime - hours(2),
    during_match = kickoff_datetime,
    post_0_2 = match_end_datetime,
    post_2_4 = match_end_datetime + hours(2)
  )
  window_ends <- c(
    pre_4_2 = kickoff_datetime - hours(2),
    pre_2_0 = kickoff_datetime,
    during_match = match_end_datetime,
    post_0_2 = match_end_datetime + hours(2),
    post_2_4 = match_end_datetime + hours(4)
  )

  names(window_starts) |>
    map(function(win) {
      # Convert each continuous interval into the set of whole hours that overlap
      # it, because the panel outcome is measured at hourly resolution.
      seq_hours <- seq(
        floor_date(window_starts[[win]], unit = "hour"),
        ceiling_date(window_ends[[win]], unit = "hour") - hours(1),
        by = "1 hour"
      )

      tibble(
        match_id = match_id,
        club = club,
        stadium_id = stadium_id,
        stadium_name = stadium_name,
        home_away = home_away,
        competition = competition,
        season = season,
        window_name = win,
        hour_start = force_tz(seq_hours, tzone = timezone)
      )
    }) |>
    bind_rows()
}

# Aggregate incident-level crimes into hourly stadium-distance-category counts.
aggregate_crime_counts <- function(crime_assigned) {
  message("Aggregating crime incidents to stadium-hour-distance-category panel...")

  # This produces the observed crime-count cells before zero-crime hours are
  # added back in by the full panel constructor.
  st_drop_geometry(crime_assigned) |>
    transmute(
      stadium_id,
      stadium_name,
      distance_band = as.character(distance_band),
      date = as.Date(crime_datetime, tz = tz_london),
      hour_start = crime_datetime,
      crime_category,
      crime_count = 1L,
      overlapping_catchment = overlapping_catchment
    ) |>
    group_by(stadium_id, stadium_name, distance_band, date, hour_start, crime_category) |>
    summarise(
      crime_count = sum(crime_count),
      overlapping_catchment = first(overlapping_catchment),
      .groups = "drop"
    )
}

# Build the zero-filled modelling panel and enrich it with time controls.
build_model_data <- function(
    crime_counts,
    match_hour_panel,
    config,
    stadium_setup,
    fixtures = NULL) {
  message("Building complete zero-filled model panel...")

  # Build the complete stadium-band-hour-category skeleton first, then join in
  # crime counts and match exposures so that zero-crime hours are explicit.
  crime_categories <- config$crime_categories

  analysis_hours <- seq(
    as.POSIXct(paste0(config$analysis_start_date, " 00:00:00"), tz = config$timezone),
    as.POSIXct(paste0(config$analysis_end_date, " 23:00:00"), tz = config$timezone),
    by = "1 hour"
  )

  analysis_units <- get_analysis_units(config)

  panel_skeleton <- crossing(
    analysis_units |>
      select(stadium_id, club, stadium_name, unit_start_date, unit_end_date),
    distance_band = paste0(head(config$distance_breaks, -1), "_", tail(config$distance_breaks, -1), "m"),
    hour_start = analysis_hours,
    crime_category = crime_categories
  ) |>
    mutate(
      date = as.Date(hour_start, tz = config$timezone)
    ) |>
    filter(date >= unit_start_date, date <= unit_end_date) |>
    select(-unit_start_date, -unit_end_date)

  bank_holidays <- load_bank_holidays(config$bank_holiday_path, timezone = config$timezone)

  panel <- panel_skeleton |>
    left_join(crime_counts, by = c(
      "stadium_id", "stadium_name", "distance_band",
      "date", "hour_start", "crime_category"
    )) |>
    mutate(
      crime_count = replace_na(crime_count, 0L),
      overlapping_catchment = replace_na(overlapping_catchment, FALSE)
    ) |>
    left_join(match_hour_panel, by = c("stadium_id", "club", "hour_start")) |>
    mutate(
      across(
        c(
          home_pre_4_2, home_pre_2_0, home_during_match, home_post_0_2, home_post_2_4,
          away_pre_4_2, away_pre_2_0, away_during_match, away_post_0_2, away_post_2_4,
          home_match_any, away_match_any, non_match_hour
        ),
        ~ replace_na(.x, 0L)
      ),
      competition_home = replace_na(competition_home, "none"),
      competition_away = replace_na(competition_away, "none"),
      distance_band = factor(
        distance_band,
        levels = paste0(head(config$distance_breaks, -1), "_", tail(config$distance_breaks, -1), "m")
      ),
      hour_of_day = hour(hour_start),
      day_of_week = wday(hour_start, week_start = 1, label = TRUE, abbr = TRUE),
      hour_of_week = factor(
        (wday(hour_start, week_start = 1) - 1) * 24 + hour_of_day
      ),
      month = month(date),
      year = year(date),
      season = make_season(date),
      weekend = as.integer(wday(date, week_start = 1) >= 6)
    ) |>
    left_join(bank_holidays, by = "date") |>
    mutate(
      # These are the fixed-effect identifiers and lower-dimensional calendar
      # controls used in the main and alternative specifications.
      bank_holiday = replace_na(bank_holiday, FALSE),
      stadium_distance_fe = interaction(stadium_id, distance_band, drop = TRUE, lex.order = TRUE),
      date_fe = factor(date),
      time_index_days = as.numeric(date - min(date)),
      time_index_days_sq = time_index_days^2,
      league_match_home = as.integer(
          str_detect(
          competition_home,
          "league|premier|championship|league_one|league_two"
        )
      ),
      kickoff_weekend = as.integer(weekend == 1)
    )

  if (!is.null(fixtures)) {
    overlap_summary <- fixtures |>
      filter(home_away == "home") |>
      count(stadium_id, kickoff_datetime) |>
      filter(n > 1)

    if (nrow(overlap_summary) > 0) {
      warning(
        "Some home fixtures overlap at the same stadium and kickoff time. Inspect fixtures_clean.rds."
      )
    }
  }

  panel
}

# Convert a fixest model object into a tidy coefficient table with confidence intervals.
tidy_fixest_model <- function(model, model_name, crime_category) {
  coef_est <- coef(model)
  conf_mat <- suppressWarnings(confint(model))

  tibble(
    term = names(coef_est),
    estimate = unname(coef_est),
    conf_low = conf_mat[, 1],
    conf_high = conf_mat[, 2],
    model = model_name,
    crime_category = crime_category
  )
}

# Fit the main PPML specifications for each crime category and save outputs.
fit_ppml_models <- function(model_data, config, output_dir) {
  ensure_dir(output_dir)

  message("Estimating PPML models...")

  # Model 1 is the main event-window specification. Model 2 collapses exposure
  # to any home or away match. The alternative-calendar model drops date fixed
  # effects in favour of richer calendar controls.
  home_window_terms <- c(
    "home_pre_4_2", "home_pre_2_0", "home_during_match", "home_post_0_2", "home_post_2_4"
  )
  away_window_terms <- c(
    "away_pre_4_2", "away_pre_2_0", "away_during_match", "away_post_0_2", "away_post_2_4"
  )

  event_formula <- as.formula(paste(
    "crime_count ~ 0 +",
    paste(sprintf("%s:distance_band", home_window_terms), collapse = " + "),
    "+",
    paste(sprintf("%s:distance_band", away_window_terms), collapse = " + "),
    "| stadium_distance_fe + hour_of_week + date_fe"
  ))

  robust_formula <- as.formula(
    "crime_count ~ 0 + home_match_any:distance_band + away_match_any:distance_band | stadium_distance_fe + hour_of_week + date_fe"
  )

  alt_calendar_formula <- as.formula(paste(
    "crime_count ~ 0 +",
    paste(sprintf("%s:distance_band", home_window_terms), collapse = " + "),
    "+",
    paste(sprintf("%s:distance_band", away_window_terms), collapse = " + "),
    "+ bank_holiday + weekend + i(month) + i(year) + time_index_days + time_index_days_sq",
    "| stadium_distance_fe + hour_of_week"
  ))

  results <- map(config$crime_categories, function(crime_cat) {
    dat <- model_data |>
      filter(crime_category == crime_cat)

    # Estimate one crime category at a time so each coefficient table is easy
    # to interpret and plot.
    model_1 <- fepois(
      formula = event_formula,
      data = dat,
      vcov = ~ stadium_id + date_fe
    )

    model_2 <- fepois(
      formula = robust_formula,
      data = dat,
      vcov = ~ stadium_id + date_fe
    )

    model_1_alt_calendar <- fepois(
      formula = alt_calendar_formula,
      data = dat,
      vcov = ~ stadium_id + date_fe
    )

    list(
      crime_category = crime_cat,
      model_1 = model_1,
      model_2 = model_2,
      model_1_alt_calendar = model_1_alt_calendar
    )
  })

  names(results) <- config$crime_categories

  coefficients_tbl <- bind_rows(
    map(results, ~ tidy_fixest_model(.x$model_1, "model_1", .x$crime_category)),
    map(results, ~ tidy_fixest_model(.x$model_2, "model_2", .x$crime_category)),
    map(
      results,
      ~ tidy_fixest_model(
        .x$model_1_alt_calendar,
        "model_1_alt_calendar",
        .x$crime_category
      )
    )
  ) |>
    mutate(
      exp_estimate = exp(estimate),
      exp_conf_low = exp(conf_low),
      exp_conf_high = exp(conf_high)
    )

  write_csv(
    coefficients_tbl,
    file.path(output_dir, "ppml_coefficients.csv")
  )

  for (crime_cat in names(results)) {
    sink(file.path(output_dir, paste0("summary_", crime_cat, ".txt")))
    print(summary(results[[crime_cat]]$model_1))
    cat("\n\n")
    print(summary(results[[crime_cat]]$model_2))
    cat("\n\n")
    print(summary(results[[crime_cat]]$model_1_alt_calendar))
    sink()
  }

  list(
    models = results,
      coefficients = coefficients_tbl
  )
}

# Extract home-match event-study coefficients into plot-friendly components.
extract_home_event_coefficients <- function(coefficients_tbl) {
  coefficients_tbl |>
    filter(str_detect(term, "^home_")) |>
    extract(
      term,
      into = c("window", "distance_band"),
      regex = "^(home_[a-z0-9_]+):distance_band(.+)$",
      remove = FALSE
    ) |>
    mutate(
      window = factor(
        window,
        levels = c(
          "home_pre_4_2", "home_pre_2_0", "home_during_match",
          "home_post_0_2", "home_post_2_4"
        ),
        labels = c(
          "Pre 4 to 2h",
          "Pre 2h to KO",
          "During match",
          "Post 0 to 2h",
          "Post 2 to 4h"
        )
      ),
      distance_band = as.character(distance_band)
    )
}

# Plot home-match event-study estimates by event window and distance band.
plot_event_study <- function(coefficients_tbl, crime_category, output_dir) {
  ensure_dir(output_dir)

  # The plotted values are incidence-rate ratios because PPML coefficients are
  # exponentiated before plotting.
  plot_tbl <- coefficients_tbl |>
    filter(crime_category == !!crime_category, model == "model_1") |>
    extract_home_event_coefficients()

  p <- ggplot(
    plot_tbl,
    aes(
      x = window,
      y = exp_estimate,
      ymin = exp_conf_low,
      ymax = exp_conf_high,
      color = distance_band,
      group = distance_band
    )
  ) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    labs(
      title = paste("Home-match event-study coefficients:", crime_category),
      x = NULL,
      y = "Incidence rate ratio",
      color = "Distance band"
    ) +
    theme_minimal(base_size = 12)

  ggsave(
    filename = file.path(output_dir, paste0("event_study_", crime_category, ".png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )

  invisible(p)
}

# Shift treatment windows forward in time for placebo tests.
run_placebo_shift <- function(model_data, shift_hours = 24L) {
  # Keep the observed outcomes fixed and move only treatment timing to test
  # whether the specification finds spurious effects.
  exposure_cols <- c(
    "home_pre_4_2", "home_pre_2_0", "home_during_match", "home_post_0_2", "home_post_2_4",
    "away_pre_4_2", "away_pre_2_0", "away_during_match", "away_post_0_2", "away_post_2_4",
    "home_match_any", "away_match_any", "non_match_hour"
  )

  shifted_exposure <- model_data |>
    select(stadium_id, club, hour_start, all_of(exposure_cols)) |>
    mutate(hour_start = hour_start + hours(shift_hours))

  model_data |>
    select(-all_of(exposure_cols)) |>
    left_join(shifted_exposure, by = c("stadium_id", "club", "hour_start")) |>
    mutate(across(all_of(exposure_cols), ~ replace_na(.x, 0L)))
}

# Reassign crimes to a user-specified alternative set of distance bands.
make_alternative_distance_assignment <- function(crime_sf, stadium_setup, new_breaks) {
  nearest_idx <- st_nearest_feature(crime_sf, stadium_setup$stadiums_sf)
  distance_m <- as.numeric(
    st_distance(
      crime_sf,
      stadium_setup$stadiums_sf[nearest_idx, ],
      by_element = TRUE
    )
  )

  crime_sf |>
    mutate(
      distance_to_stadium_m = distance_m,
      distance_band = label_distance_band(distance_m, new_breaks)
    ) |>
    filter(!is.na(distance_band))
}

# Estimate the preferred PPML specifications on useful robustness subsamples.
fit_optional_subsample_models <- function(model_data, config, output_dir) {
  ensure_dir(output_dir)

  # These subsamples correspond to the main robustness exercises described in
  # the analysis plan.
  subsamples <- list(
    weekday_only = model_data |> filter(weekend == 0),
    weekend_only = model_data |> filter(weekend == 1),
    non_overlapping = model_data |> filter(!overlapping_catchment),
    league_home_only = model_data |> filter(league_match_home == 1 | home_match_any == 0)
  )

  out <- imap(subsamples, function(dat, name) {
    if (nrow(dat) == 0) {
      return(NULL)
    }

    fit_ppml_models(dat, config, file.path(output_dir, name))
  })

  out
}
