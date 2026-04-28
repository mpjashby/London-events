#!/usr/bin/env Rscript

# Football fixtures builder from local CSV files
# ----------------------------------------------
# This script reads every match-level CSV in `original_data/football_matches/`
# and combines them into one tibble with one row per match.

source(here::here("R", "00_helpers.R"))

# USER INPUTS ---------------------------------------------------------------

config <- list(
  timezone = tz_london,
  football_matches_dir = here("original_data", "football_matches"),
  output_dir = here("derived_data")
)

# HELPERS -------------------------------------------------------------------

extract_competition_from_path <- function(path) {
  basename(path) |>
    str_remove("-matches-\\d{4}-to-\\d{4}-stats\\.csv$") |>
    str_replace_all("-", " ") |>
    str_squish()
}

extract_season_from_path <- function(path) {
  season_bits <- str_match(
    basename(path),
    "-matches-(\\d{4})-to-(\\d{4})-stats\\.csv$"
  )

  if_else(
    !is.na(season_bits[, 2]),
    paste0(season_bits[, 2], "/", str_sub(season_bits[, 3], 3, 4)),
    NA_character_
  )
}

parse_match_datetime <- function(x, tz = tz_london) {
  parse_date_time(
    x,
    orders = c("b d Y - I:M p", "b d Y - I:M%p"),
    tz = tz,
    exact = FALSE
  )
}

parse_attendance <- function(x) {
  x |>
    as.character() |>
    na_if("") |>
    na_if("N/A") |>
    str_replace_all(",", "") |>
    as.numeric()
}

parse_score <- function(x) {
  x |>
    as.character() |>
    na_if("") |>
    na_if("N/A") |>
    na_if("-1") |>
    as.integer()
}

read_match_file <- function(path, tz = tz_london) {
  message("Reading ", basename(path), "...")

  raw_matches <- read_csv(path, show_col_types = FALSE, progress = FALSE)

  match_datetime <- parse_match_datetime(raw_matches$date_GMT, tz = tz)

  raw_matches |>
    transmute(
      match_date = as.Date(match_datetime, tz = tz),
      start_time = format(match_datetime, "%H:%M"),
      venue = stadium_name,
      home_team = home_team_name,
      away_team = away_team_name,
      attendance = parse_attendance(attendance),
      home_goals = parse_score(home_team_goal_count),
      away_goals = parse_score(away_team_goal_count),
      competition = extract_competition_from_path(path),
      season = extract_season_from_path(path),
      source_file = basename(path)
    )
}

# LOAD AND COMBINE MATCH FILES ----------------------------------------------

match_files <- list.files(
  config$football_matches_dir,
  pattern = "matches-\\d{4}-to-\\d{4}-stats\\.csv$",
  full.names = TRUE
) |>
  sort()

fixtures_clean <- match_files |>
  map(read_match_file, tz = config$timezone) |>
  bind_rows() |>
  filter(
    !if_all(
      c(match_date, start_time, venue, home_team, away_team),
      is.na
    )
  ) |>
  distinct() |>
  arrange(match_date, start_time, home_team, away_team) |>
  mutate(venue = str_remove(venue, " \\(.+\\)$")) |>
  filter(
    venue %in%
      c(
        "Boleyn Ground",
        "Brentford Community Stadium",
        "Craven Cottage",
        "The Den",
        "Emirates Stadium",
        "Griffin Park",
        "Gtech Community Stadium",
        "Kiyan Prince Foundation Stadium",
        "London Stadium",
        "Loftus Road",
        "Loftus Road Stadium",
        "MATRADE Loftus Road",
        "Selhurst Park",
        "Stamford Bridge",
        "Tottenham Hotspur Stadium",
        "The Valley",
        "Wembley Stadium",
        "White Hart Lane"
      )
  ) |>
  mutate(
    venue = replace_values(
      venue,
      "Gtech Community Stadium" ~ "Brentford Community Stadium",
      c(
        "Kiyan Prince Foundation Stadium",
        "Loftus Road Stadium",
        "MATRADE Loftus Road"
      ) ~ "Loftus Road"
    )
  )

# SAVE OUTPUTS --------------------------------------------------------------

write_rds(
  fixtures_clean,
  file.path(config$output_dir, "fixtures_clean.rds"),
  compress = "gz"
)

write_csv(
  fixtures_clean,
  file.path(config$output_dir, "fixtures_clean.csv")
)
