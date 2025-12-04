#==============================================================================
# Lake BT Detection Data - Initial Cleaning and Preparation
#==============================================================================
# Purpose: Clean and filter acoustic telemetry data from Lake BT
# Author: Marcus Michelangeli
#
# Input files:
#   - ./raw_tracking_data/raw/BT_B/results/animal/all.csv
#   - ./data/fish_size/biometric_data.csv
#   - ./data/fish_size/biometric_post_exp_data.csv
#   - ./raw_tracking_data/sunset_sunrise_data.csv
#   - ./data/lake_coords/lake_BT_polygon.gpkg
#
# Output:
#   - ./data/tracks_filtered/lake_BT/01_lake_BT_sub.rds
#==============================================================================

# Load required libraries ---------------------------------------------------
library(data.table)
library(tidyverse)
library(move)
library(move2)
library(mapedit)
library(sf)
library(lubridate)  # Explicit loading for date/time functions

# Define file paths ---------------------------------------------------------
save_filtered_data_path <- "./data/tracks_filtered/lake_BT/"
size_data_path <- "./data/fish_size/"
raw_data_path <- "./raw_tracking_data/"

cat("\014")  # Clear console

#==============================================================================
# 1. IMPORT AND INITIAL DATA PREPARATION
#==============================================================================

# Import raw detection data -------------------------------------------------
lake_BT <- fread(paste0(raw_data_path, "raw/BT_B/results/animal/all.csv"))
message("Imported ", nrow(lake_BT), " raw detections") #22967436 detections

# Fix swapped latitude/longitude columns ------------------------------------
lake_BT <- lake_BT %>%
  rename(Lat = Longitude, 
         Long = Latitude)

# Prepare timestamp columns -------------------------------------------------
# Original timestamps are in UTC; convert to local time (Europe/Stockholm)
lake_BT <- lake_BT %>%
  rename(timestamp_utc = Time) %>%
  mutate(
    timestamp_cest = with_tz(timestamp_utc, tzone = "Europe/Stockholm"),
    date_utc = as.Date(timestamp_utc),
    date_cest = as.Date(format(timestamp_cest, tz = "Europe/Stockholm"))
  )

# Verify timezone conversions
message("Non-aligned dates (UTC vs CEST): ", 
        sum(lake_BT$date_utc != lake_BT$date_cest)) #1755991

#==============================================================================
# 2. MERGE WITH BIOMETRIC DATA
#==============================================================================

# Import fish biometric data ------------------------------------------------
biometrics <- fread(paste0(size_data_path, "biometric_data.csv"))

# Extract individual ID from tag number and prepare for merging
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# Merge biometric data with detections --------------------------------------
lake_BT <- lake_BT %>%
  rename(individual_ID = Id) %>%
  left_join(biometrics, by = "individual_ID")

# Label reference tag and filter for Lake BT individuals -------------------
lake_BT <- lake_BT %>%
  mutate(individual_ID = if_else(FullId == 'H170-1802-65065', 
                                 'Reference', 
                                 as.character(individual_ID))) %>%
  filter(Lake == 'BT' | individual_ID == 'Reference')

message("After lake filter: ", nrow(lake_BT), " detections") #22957599 

#==============================================================================
# 3. TEMPORAL FILTERING AND DATA CLEANING
#==============================================================================

# Filter for post-introduction period (after pike introduction) -------------
pike_intro_date <- as.Date("2022-09-26")

lake_BT <- lake_BT %>%
  filter(date_cest >= pike_intro_date) %>%
  # Remove unnecessary columns to reduce memory footprint
  dplyr::select(-c(ID, Notes, Transmitter, `Tag Number`, `PIT Number`, `Tag Type`, 
            `Biologger Number`, Station, HPEm, TempData, DepthData, AccelData,
            RxDetected, nRxDetected, `Serial Number`, nRxUsed, RxUsed,
            timestamp_utc, date_utc, Date, Lake))

message("After temporal filter: ", nrow(lake_BT), " detections") #20832375

# Standardize individual IDs ------------------------------------------------
# Add 'F' prefix to distinguish fish individuals from reference
lake_BT <- lake_BT %>%
  mutate(individual_ID = paste0("F", individual_ID))

#==============================================================================
# 4. CONVERT TO MOVEMENT OBJECT AND SPATIAL FILTERING
#==============================================================================

# Convert to move2 object for spatial analysis -----------------------------
lake_BT_mv <- mt_as_move2(
  lake_BT,
  coords = c("Long", "Lat"),
  crs = "WGS84",
  time_column = "timestamp_cest",
  track_id_column = "individual_ID",
  na.fail = FALSE  # Allow missing coordinates
)

# Sort chronologically by individual and timestamp --------------------------
lake_BT_mv <- lake_BT_mv %>%
  arrange(individual_ID, timestamp_cest)

# Load lake boundary polygon ------------------------------------------------
lake_BT_poly <- st_read("./data/lake_coords/lake_BT_polygon.gpkg")

# Optional: Visualize reference tag to verify spatial accuracy
# ref_tag <- filter_track_data(lake_BT_mv, .track_id = "FReference")
# lake_BT_map <- mapview::mapView(mt_track_lines(ref_tag)$geometry)

# Note: To create or update the lake polygon, uncomment and run:
# BT_map <- mapedit::drawFeatures(map = lake_BT_map)
# st_write(BT_map, dsn="data/lake_coords/lake_BT_polygon.gpkg", 
#          driver="GPKG", delete_layer = TRUE)

# Filter detections within lake boundaries ---------------------------------
lake_BT_sub <- st_filter(lake_BT_mv, lake_BT_poly)
message("After spatial filter: ", nrow(lake_BT_sub), " detections") #20267565

# Extract coordinates back to standard columns ------------------------------
coords <- st_coordinates(lake_BT_sub$geometry)
lake_BT_sub$Long <- coords[, 1]
lake_BT_sub$Lat <- coords[, 2]

#==============================================================================
# 5. CREATE TEMPORAL CLASSIFICATION VARIABLES
#==============================================================================

# Classify experiment stage (Early vs Late) ---------------------------------
# Split dataset into two equal temporal periods
unique_dates <- unique(sort(lake_BT_sub$date_cest))
n_dates <- length(unique_dates)
early_dates <- unique_dates[1:ceiling(n_dates/2)]

lake_BT_sub <- lake_BT_sub %>%
  mutate(Stage = if_else(date_cest %in% early_dates, 'Early', 'Late'))

# Verify stage classification
message("Stage distribution:")
print(table(lake_BT_sub$Stage))
# Early     Late 
# 13299976  6967589 

# Classify time of day (Day vs Night) --------------------------------------
# Import sunrise/sunset data
sun_data <- fread(paste0(raw_data_path, "sunset_sunrise_data.csv"))

sun_data <- sun_data %>%
  dplyr::select(Sunrise_timestamp, Sunset_timestamp) %>%
  rename(sunrise = Sunrise_timestamp, sunset = Sunset_timestamp) %>%
  mutate(
    sunrise = dmy_hm(sunrise, tz = "Europe/Stockholm"),
    sunset = dmy_hm(sunset, tz = "Europe/Stockholm"),
    date_cest = as.Date(sunset)
  )

# Merge sun data and classify day/night -------------------------------------
lake_BT_sub <- lake_BT_sub %>%
  left_join(sun_data, by = "date_cest") %>%
  mutate(
    time_of_day = if_else(
      timestamp_cest >= sunrise & timestamp_cest <= sunset,
      "Day",
      "Night"
    )
  )

# Verify day/night classification
message("Time of day distribution:")
print(table(lake_BT_sub$time_of_day))
# Day    Night 
# 9615730 10651835 


#==============================================================================
# 6. MERGE POST-EXPERIMENT SURVIVAL DATA
#==============================================================================

# Import post-experiment biometric data -------------------------------------
post_biometrics <- fread(paste0(size_data_path, "biometric_post_exp_data.csv"))

post_size_cols <- post_biometrics %>%
  filter(Lake == 'BT') %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number))) %>%
  dplyr::select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found)

# Merge survival information ------------------------------------------------
lake_BT_sub <- lake_BT_sub %>%
  left_join(post_size_cols, by = "individual_ID")


#==============================================================================
# 7. SAVE CLEANED DATASET
#==============================================================================

# Save as RDS file ----------------------------------------------------------
output_file <- paste0(save_filtered_data_path, "01_lake_BT_sub.rds")
saveRDS(lake_BT_sub, file = output_file)

message("Final dataset: ", nrow(lake_BT_sub), " detections") #20267565 detections
message("Number of individuals: ", n_distinct(lake_BT_sub$individual_ID)) #67 individuals (including reference tag)
message("Date range: ", min(lake_BT_sub$date_cest), " to ", 
        max(lake_BT_sub$date_cest)) #Date range: 2022-09-26 to 2022-10-29
message("Saved to: ", output_file)
