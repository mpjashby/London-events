# This code loads data from the various files provided by the Met and outputs a
# single CSV file containing all the data.
#
# The file of data for 2010 to 2019 is provided in two files. The 'CRIS' file
# contains a unique ID for each incident and the incident co-ordinates. The
# 'TNO' file contains the date/time of the crime and the crime type. This code
# joins the co-ordinates to the offence records.
#
# The files for later years have all the fields in the same file, but the data
# split across multiple files by rows.

pacman::p_load(here, janitor, tidyverse)

# LOAD AND JOIN -------------------------------------------------------------

# Read the full TNO file, drop the unnamed export column, and standardise the
# remaining field names to snake case for downstream joins.
tno_2010 <- here(
  "original_data",
  "MPS All Crime",
  "2013-2019 (Yearly subsets)",
  "Data to use",
  "MPS-Merged-TNO-2019-10.csv"
) |>
  read_csv() |>
  clean_names() |>
  select(-x1, -date_cr_reported) |>
  rename(offences = number_of_total_notifiable_offences)

# Read the CRIS file and keep only the join key and coordinate columns
cris_2010 <- here(
  "original_data",
  "MPS All Crime",
  "2013-2019 (Yearly subsets)",
  "Data to use",
  "MPS-Merged-CRIS-2019-10.csv"
) |>
  read_csv() |>
  clean_names() |>
  select(unique_id, easting, northing) |>
  # This file contains some duplicate rows, so remove those
  slice_head(n = 1, by = unique_id)

# Keep every TNO field and add the CRIS easting and northing values using the
# shared unique crime identifier.
crimes_2010 <- left_join(tno_2010, cris_2010, by = join_by(unique_id))

# Load 2019 to 2020 data (seems to be missing September 2019)
crimes_2019 <- c(
  "original_data/MPS All Crime/2019-2022/UCL_Data_share_1920.csv",
  "original_data/MPS All Crime/2019-2022/UCL_Data_share_1_1920.csv",
  "original_data/MPS All Crime/2019-2022/UCL_Data_share_2_1920.csv",
  "original_data/MPS All Crime/2019-2022/UCL_Data_share_3_1920.csv"
) |>
  map(read_csv) |>
  bind_rows(.id = "file") |>
  clean_names() |>
  select(
    unique_id = cr_no,
    date_cr_recorded,
    date_committed_from,
    date_committed_to,
    new_major_text,
    new_minor_text,
    offences,
    easting = x,
    northing = y
  ) |>
  mutate(date_committed_to = as_datetime(date_committed_to))

# Load 2021 to 2022 data
crimes_2021 <- c(
  "original_data/MPS All Crime/2019-2022/UCL_Data._2122csv.csv",
  "original_data/MPS All Crime/2019-2022/UCL_Data__2122_1.csv",
  "original_data/MPS All Crime/2019-2022/UCL_Data__2122_2.csv",
  "original_data/MPS All Crime/2019-2022/UCL_Data__2122_3.csv"
) |>
  map(read_csv) |>
  bind_rows(.id = "file") |>
  clean_names() |>
  select(
    unique_id = cr_no,
    date_cr_recorded,
    date_committed_from,
    date_committed_to,
    new_major_text,
    new_minor_text,
    offences,
    easting = x,
    northing = y
  )

# Load 2022 to 2023 data
crimes_2022 <- c(
  "original_data/MPS All Crime/2022-2023/UCL_Data.csv",
  "original_data/MPS All Crime/2022-2023/UCL_Data_1.csv",
  "original_data/MPS All Crime/2022-2023/UCL_Data_2.csv",
  "original_data/MPS All Crime/2022-2023/UCL_Data_3.csv",
  "original_data/MPS All Crime/2022-2023/UCL_Data_4.csv"
) |>
  map(read_csv) |>
  bind_rows(.id = "file") |>
  clean_names() |>
  select(
    unique_id = cr_no,
    date_cr_recorded,
    date_committed_from,
    date_committed_to,
    new_major_text,
    new_minor_text,
    offences,
    easting = x,
    northing = y
  )

# Merge crimes
# The 2019-2020 and 2021-2022 datasets overlap, because the first dataset
# actually contains data to the end of 2021
crimes <- bind_rows(
  crimes_2010,
  filter(crimes_2019, as_date(date_cr_recorded) < ymd("2020-10-01")),
  filter(crimes_2021, as_date(date_cr_recorded) < ymd("2021-10-01")),
  crimes_2022
) |>
  arrange(date_cr_recorded) |>
  write_rds(here("derived_data/crimes.rds"), compress = "gz", version = 3)
