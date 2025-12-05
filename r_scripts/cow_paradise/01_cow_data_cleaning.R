#==============================================================================
# Lake Cow Paradise Detection Data - Initial Cleaning and Preparation
#==============================================================================
# Purpose: Clean and filter acoustic telemetry data from Lake Cow Paradise
# Author: Marcus Michelangeli
# Input files:
#   - ./data/raw_tracking_data/cow_paradise_C/results/animal/all.csv
#   - ./data/fish_size/biometric_data.csv
#   - ./data/fish_size/biometric_post_exp_data.csv
#   - ./data/raw_tracking_data/sunset_sunrise_data.csv
#   - ./data/lake_coords/lake_cow_polygon.gpkg
#
# Output:
#   - ./data/tracks_filtered/lake_cow_paradise/01_lake_cow_sub.rds
#==============================================================================

# Load required libraries ---------------------------------------------------
library(data.table)
library(tidyverse)
library(move)
library(move2)
library(mapedit)
library(sf)
library(lubridate)  # Explicit loading for date/time functions

# Set timezone globally -----------------------------------------------------
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths ---------------------------------------------------------
raw_tracking_data_path <- "./data/raw_tracking_data/"
save_filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"
size_data_path <- "./data/fish_size/"

cat("\014")  # Clear console

#==============================================================================
# 1. IMPORT AND INITIAL DATA PREPARATION
#==============================================================================

# Import raw detection data -------------------------------------------------
lake_cow <- fread(paste0(raw_tracking_data_path, "cow_paradise_C/results/animal/all.csv"))
message("Imported ", nrow(lake_cow), " raw detections") #Imported 10368804 raw detections

# Fix swapped latitude/longitude columns ------------------------------------
lake_cow <- lake_cow %>%
  rename(Lat = Longitude,
         Long = Latitude)

# Prepare timestamp columns -------------------------------------------------
# Original timestamps are in UTC; convert to local time (Europe/Stockholm)
lake_cow <- lake_cow %>%
  rename(timestamp_utc = Time) %>%
  mutate(
    timestamp_cest = with_tz(timestamp_utc, tzone = "Europe/Stockholm"),
    date_utc = as.Date(timestamp_utc),
    date_cest = as.Date(format(timestamp_cest, tz = "Europe/Stockholm"))
  )

# Verify timezone conversions
message("Non-aligned dates (UTC vs CEST): ",
        sum(lake_cow$date_utc != lake_cow$date_cest))

#==============================================================================
# 2. MERGE WITH BIOMETRIC DATA
#==============================================================================

# Import fish biometric data ------------------------------------------------
biometrics <- fread(paste0(size_data_path, "biometric_data.csv"))

# Extract individual ID from tag number and prepare for merging
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# Merge biometric data with detections --------------------------------------
lake_cow <- lake_cow %>%
  rename(individual_ID = Id) %>%
  left_join(biometrics, by = "individual_ID")

# Label reference tag and filter for Lake Cow Paradise individuals ---------
lake_cow <- lake_cow %>%
  mutate(individual_ID = if_else(FullId == 'H170-1802-65064',
                                 'Reference',
                                 as.character(individual_ID))) %>%
  filter(Lake == 'Cow Paradise' | individual_ID == 'Reference')

message("After lake filter: ", nrow(lake_cow), " detections") #10364909 detections

#==============================================================================
# 3. TEMPORAL FILTERING AND DATA CLEANING
#==============================================================================

# Filter for post-introduction period (after pike introduction) -------------
pike_intro_date <- as.Date("2022-09-27")

lake_cow <- lake_cow %>%
  filter(date_cest >= pike_intro_date) %>%
  # Remove unnecessary columns to reduce memory footprint
  dplyr::select(-c(ID, Notes, Transmitter, `Tag Number`, `PIT Number`, `Tag Type`,
            `Biologger Number`, Station, HPEm, TempData, DepthData, AccelData,
            RxDetected, nRxDetected, `Serial Number`, nRxUsed, RxUsed,
            timestamp_utc, date_utc, Date, Lake))

message("After temporal filter: ", nrow(lake_cow), " detections") #9431255 detections

# Standardize individual IDs ------------------------------------------------
# Add 'F' prefix to distinguish fish individuals from reference
lake_cow <- lake_cow %>%
  mutate(individual_ID = paste0("F", individual_ID))

#==============================================================================
# 4. CONVERT TO MOVEMENT OBJECT AND SPATIAL FILTERING
#==============================================================================

# Convert to move2 object for spatial analysis -----------------------------
lake_cow_mv <- mt_as_move2(
  lake_cow,
  coords = c("Long", "Lat"),
  crs = "WGS84",
  time_column = "timestamp_cest",
  track_id_column = "individual_ID",
  na.fail = FALSE  # Allow missing coordinates
)

# Sort chronologically by individual and timestamp --------------------------
lake_cow_mv <- lake_cow_mv %>%
  arrange(individual_ID, timestamp_cest)

# Load lake boundary polygon ------------------------------------------------
lake_cow_poly <- st_read("./data/lake_coords/lake_cow_polygon.gpkg")

# # Optional: Visualize reference tag to verify spatial accuracy
# ref_tag <- filter_track_data(lake_cow_mv, .track_id = "F59818")
# lake_cow_map <- mapview::mapView(mt_track_lines(ref_tag)$geometry)

# Note: To create or update the lake polygon, uncomment and run:
# lake_cow_map <- mapedit::drawFeatures(map = lake_cow_map)
# st_write(lake_cow_map, dsn="data/lake_coords/lake_cow_polygon.gpkg",
#          driver="GPKG", delete_layer = TRUE)

# Filter detections within lake boundaries ---------------------------------
lake_cow_sub <- st_filter(lake_cow_mv, lake_cow_poly)
message("After spatial filter: ", nrow(lake_cow_sub), " detections") #9089970 detections

# Verify timezone preservation ----------------------------------------------
message("Timezone check: ", tz(lake_cow_sub$timestamp_cest))

# Extract coordinates back to standard columns ------------------------------
coords <- st_coordinates(lake_cow_sub$geometry)
lake_cow_sub$Long <- coords[, 1]
lake_cow_sub$Lat <- coords[, 2]

#==============================================================================
# 5. CREATE TEMPORAL CLASSIFICATION VARIABLES
#==============================================================================

# Classify experiment stage (Early vs Late) ---------------------------------
# Split dataset into two equal temporal periods
unique_dates <- unique(sort(lake_cow_sub$date_cest))
n_dates <- length(unique_dates)
early_dates <- unique_dates[1:ceiling(n_dates/2)]

message("\nTotal unique dates: ", n_dates) #Total unique dates: 35
message("Early period dates: ", length(early_dates)) #Early period dates: 18

lake_cow_sub <- lake_cow_sub %>%
  mutate(Stage = if_else(date_cest %in% early_dates, 'Early', 'Late'))

# Verify stage classification
message("\nStage distribution:")
print(table(lake_cow_sub$Stage))
# Early    Late 
# 6718978 2370992 


# Classify time of day (Day vs Night) --------------------------------------
# Import sunrise/sunset data
sun_data <- fread(paste0(raw_tracking_data_path, "sunset_sunrise_data.csv"))

sun_data <- sun_data %>%
  dplyr::select(Sunrise_timestamp, Sunset_timestamp) %>%
  rename(sunrise = Sunrise_timestamp, sunset = Sunset_timestamp) %>%
  mutate(
    sunrise = dmy_hm(sunrise, tz = "Europe/Stockholm"),
    sunset = dmy_hm(sunset, tz = "Europe/Stockholm"),
    date_cest = as.Date(sunset)
  )

# Merge sun data and classify day/night -------------------------------------
lake_cow_sub <- lake_cow_sub %>%
  left_join(sun_data, by = "date_cest") %>%
  mutate(
    time_of_day = if_else(
      timestamp_cest >= sunrise & timestamp_cest <= sunset,
      "Day",
      "Night"
    )
  )

# Verify day/night classification
message("\nTime of day distribution:")
print(table(lake_cow_sub$time_of_day))

# Day   Night 
# 4388787 4699445 


#==============================================================================
# 6. MERGE POST-EXPERIMENT SURVIVAL DATA
#==============================================================================

# Import post-experiment biometric data -------------------------------------
post_biometrics <- fread(paste0(size_data_path, "biometric_post_exp_data.csv"))

post_size_cols <- post_biometrics %>%
  filter(Lake == 'Cow Paradise') %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number))) %>%
  dplyr::select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found)

# Merge survival information ------------------------------------------------
lake_cow_sub <- lake_cow_sub %>%
  left_join(post_size_cols, by = "individual_ID")

#==============================================================================
# 7. SAVE CLEANED DATASET
#==============================================================================

# Save as RDS file ----------------------------------------------------------
output_file <- paste0(save_filtered_data_path, "01_lake_cow_sub.rds")
saveRDS(lake_cow_sub, file = output_file)

message("\n=== Data cleaning complete ===")
message("Final dataset: ", nrow(lake_cow_sub), " detections") #9089970 detections
message("Number of individuals: ", n_distinct(lake_cow_sub$individual_ID)) #Number of individuals: 67 (including the reference tag)
message("Date range: ", min(lake_cow_sub$date_cest), " to ",
        max(lake_cow_sub$date_cest)) #Date range: 2022-09-27 to 2022-10-31
message("Saved to: ", output_file)