#=============================================================================-
# Lake BT Detection Data - Initial Cleaning and Preparation ####
#=============================================================================-
# Purpose: Clean and filter acoustic telemetry data from Lake BT
# Author: Marcus Michelangeli
#
# Input files:
#   - ./data/raw_tracking_data/BT_B/results/animal/all.csv
#   - ./data/fish_biometrics/biometric_data.csv
#   - ./data/raw_tracking_data/sunset_sunrise_data.csv
#   - ./data/lake_coords/lake_BT_polygon.gpkg
#
# Output:
#   - ./data/tracks_filtered/BT/01_BT_sub.rds
#=============================================================================-

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
save_filtered_data_path <- "./data/tracks_filtered/BT/"
biometrics_data_path <- "./data/fish_biometrics/"

cat("\014")  # Clear console

#=============================================================================-
# 1. IMPORT AND INITIAL DATA PREPARATION ####
#=============================================================================-

# Import raw detection data -------------------------------------------------
BT <- fread(paste0(raw_tracking_data_path, "BT_B/results/animal/all.csv"))
message("Imported ", nrow(BT), " raw detections") #22967436

# Fix swapped latitude/longitude columns ------------------------------------
BT <- BT %>%
  rename(Lat = Longitude, 
         Long = Latitude)

# Prepare timestamp columns -------------------------------------------------
# Original timestamps are in UTC; convert to local time (Europe/Stockholm)
BT <- BT %>%
  rename(timestamp_utc = Time) %>%
  mutate(
    timestamp_cest = with_tz(timestamp_utc, tzone = "Europe/Stockholm"),
    date_utc = as.Date(timestamp_utc),
    date_cest = as.Date(format(timestamp_cest, tz = "Europe/Stockholm"))
  )

# Verify timezone conversions
message("Non-aligned dates (UTC vs CEST): ",
        sum(BT$date_utc != BT$date_cest)) #1755991

#=============================================================================-
# 2. MERGE WITH BIOMETRIC DATA ####
#=============================================================================-

# Import fish biometric data ------------------------------------------------
biometrics <- fread(paste0(biometrics_data_path, "biometric_data.csv"))

# Extract individual ID from tag number and prepare for merging
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# Merge biometric data with detections --------------------------------------
BT <- BT %>%
  rename(individual_ID = Id) %>%
  left_join(biometrics, by = "individual_ID")

# Label reference tag and filter for Lake BT individuals ------------
BT <- BT %>%
  mutate(individual_ID = if_else(FullId == 'H170-1802-65065',
                                 'Reference',
                                 as.character(individual_ID))) %>%
  filter(Lake == 'BT' | individual_ID == 'Reference')

message("After lake filter: ", nrow(BT), " detections") #22957599 detections

#=============================================================================-
# 3. TEMPORAL FILTERING AND DATA CLEANING ####
#=============================================================================-

# Filter for post-introduction period (after pike introduction) -------------
pike_intro_date <- as.Date("2022-09-26")

BT <- BT %>%
  filter(date_cest >= pike_intro_date) %>%
  # Remove unnecessary columns to reduce memory footprint
  dplyr::select(-c(ID, Notes, Transmitter, `Tag Number`, `PIT Number`, `Tag Type`,
                   `Biologger Number`, Station, HPEm, TempData, DepthData, AccelData,
                   RxDetected, nRxDetected, `Serial Number`, nRxUsed, RxUsed,
                   timestamp_utc, date_utc, Date, Lake))

message("After temporal filter: ", nrow(BT), " detections") #20832375 detections

# Standardize individual IDs ------------------------------------------------
# Add 'F' prefix to distinguish fish individuals from reference
BT <- BT %>%
  mutate(individual_ID = paste0("F", individual_ID))

#=============================================================================-
# 4. CONVERT TO MOVEMENT OBJECT AND SPATIAL FILTERING ####
#=============================================================================-

# Convert to move2 object for spatial analysis -----------------------------
BT_mv <- mt_as_move2(
  BT,
  coords = c("Long", "Lat"),
  crs = "WGS84",
  time_column = "timestamp_cest",
  track_id_column = "individual_ID",
  na.fail = FALSE  # Allow missing coordinates
)

# Sort chronologically by individual and timestamp --------------------------
BT_mv <- BT_mv %>%
  arrange(individual_ID, timestamp_cest)

# Load lake boundary polygon ------------------------------------------------
BT_polygon <- st_read("./data/lake_params/polygons/BT_polygon.gpkg")

# Optional: Visualize reference tag to verify spatial accuracy
# ref_tag <- filter_track_data(BT_mv, .track_id = "FReference")
# BT_map <- mapview::mapView(mt_track_lines(ref_tag)$geometry)

# Note: To create or update the lake polygon, uncomment and run:
# BT_map <- mapedit::drawFeatures(map = BT_map)
# st_write(BT_map, dsn="data/lake_coords/lake_BT_polygon.gpkg",
#          driver="GPKG", delete_layer = TRUE)

# Filter detections within lake boundaries ---------------------------------
BT_sub <- st_filter(BT_mv, BT_polygon)
message("After spatial filter: ", nrow(BT_sub), " detections") #20267565 detections

# Verify timezone preservation ----------------------------------------------
message("Timezone check: ", tz(BT_sub$timestamp_cest))

# Extract coordinates back to standard columns ------------------------------
coords <- st_coordinates(BT_sub$geometry)
BT_sub$Long <- coords[, 1]
BT_sub$Lat <- coords[, 2]

#=============================================================================-
# 5. CREATE TEMPORAL CLASSIFICATION VARIABLES ####
#=============================================================================-

# Classify experiment stage (Early vs Late) ---------------------------------
# Split dataset into two equal temporal periods
unique_dates <- unique(sort(BT_sub$date_cest))
n_dates <- length(unique_dates)
early_dates <- unique_dates[1:ceiling(n_dates/2)]

message("\nTotal unique dates: ", n_dates) #34
message("Early period dates: ", length(early_dates)) #17

BT_sub <- BT_sub %>%
  mutate(Stage = if_else(date_cest %in% early_dates, 'Early', 'Late'))

# Verify stage classification
message("\nStage distribution:")
print(table(BT_sub$Stage))
# Early     Late 
# 13299976  6967589  

# Classify time of day (Day vs Night) --------------------------------------
# Import sunrise/sunset data
sun_data <- fread(paste0("./data/lake_params/sunset_sunrise_data.csv"))

sun_data <- sun_data %>%
  dplyr::select(Sunrise_timestamp, Sunset_timestamp) %>%
  rename(sunrise = Sunrise_timestamp, sunset = Sunset_timestamp) %>%
  mutate(
    sunrise = dmy_hm(sunrise, tz = "Europe/Stockholm"),
    sunset = dmy_hm(sunset, tz = "Europe/Stockholm"),
    date_cest = as.Date(sunset)
  )

# Merge sun data and classify day/night -------------------------------------
BT_sub <- BT_sub %>%
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
print(table(BT_sub$time_of_day))
# Day    Night 
# 9615730 10651835 

#=============================================================================-
# 6. MERGE POST-EXPERIMENT SURVIVAL DATA ####
#=============================================================================-

# Import post-experiment biometric data -------------------------------------
post_biometrics <- fread(paste0(biometrics_data_path, "biometric_post_exp_data.csv"))

post_size_cols <- post_biometrics %>%
  filter(Lake == 'BT') %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number))) %>%
  dplyr::select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found,
         known_predated = Known_predated)

# Merge survival information ------------------------------------------------
BT_sub <- BT_sub %>%
  left_join(post_size_cols, by = "individual_ID")


#=============================================================================-
# 7. FILTER EARLY MORTALITY EVENTS ####
#=============================================================================-

# Note: Update individual ID if there are known early mortality in BT
# Remove individual that died on first day  -----------------
# 1 pike died on the first day
BT_sub <- BT_sub %>%
  filter(individual_ID != "F59889")


#=============================================================================-
# 8. SAVE CLEANED DATASET ####
#=============================================================================-

# Save as RDS file ----------------------------------------------------------
output_file <- paste0(save_filtered_data_path, "01_BT_sub.rds")
saveRDS(BT_sub, file = output_file)

message("Final dataset: ", nrow(BT_sub), " detections") #11098640 detections
message("Number of individuals: ", n_distinct(BT_sub$individual_ID)) #66 individuals (including reference tag)
message("Date range: ", min(BT_sub$date_cest), " to ",
        max(BT_sub$date_cest)) #Date range: 2022-09-25 to 2022-10-30 
message("Saved to: ", output_file)


#=============================================================================-
#9. EXTRACT TRUE LAKE BOUNDARY COORDINATES ####
#=============================================================================-

# Get coordinates from the polygon
coords <- st_coordinates(BT_polygon)

# Create a clean data frame with just longitude and latitude
lake_boundary_coords <- data.frame(
  longitude = coords[, "X"],
  latitude = coords[, "Y"]
)

# Remove duplicate closing point if it exists
if(nrow(lake_boundary_coords) > 1) {
  if(all(lake_boundary_coords[1, ] == lake_boundary_coords[nrow(lake_boundary_coords), ])) {
    lake_boundary_coords <- lake_boundary_coords[-nrow(lake_boundary_coords), ]
  }
}

message("Number of boundary points: ", nrow(lake_boundary_coords))

# Calculate lake properties
area_m2 <- st_area(BT_polygon)
area_hectares <- as.numeric(area_m2) / 10000
centroid <- st_centroid(BT_polygon)
centroid_coords <- st_coordinates(centroid)
bbox <- st_bbox(BT_polygon)

message("Lake area: ", round(area_hectares, 4), " hectares")
message("Centroid: Lon ", round(centroid_coords[1], 6), 
        ", Lat ", round(centroid_coords[2], 6))

# Save boundary coordinates as CSV
write.csv(lake_boundary_coords, "./data/lake_params/BT_boundary_coordinates.csv", 
          row.names = FALSE)

# Create summary object for use in other scripts
BT_lake_data <- list(
  boundary_coords = lake_boundary_coords,
  centroid = c(longitude = centroid_coords[1], latitude = centroid_coords[2]),
  area_hectares = as.numeric(area_hectares),
  area_m2 = as.numeric(area_m2),
  bbox = as.list(bbox),
  num_points = nrow(lake_boundary_coords)
)

# Save as RDS for easy loading in other scripts
saveRDS(BT_lake_data, "./data/lake_params/BT_coord_data.rds")
