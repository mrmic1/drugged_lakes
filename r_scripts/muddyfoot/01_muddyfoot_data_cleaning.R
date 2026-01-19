#==============================================================================
# Lake Muddyfoot Detection Data - Initial Cleaning and Preparation
#==============================================================================
# Purpose: Clean and filter acoustic telemetry data from Lake Muddyfoot
# Author: Marcus Michelangeli
#
# Input files:
#   - ./data/raw_tracking_data/muddyfoot_A/results/animal/all.csv
#   - ./data/fish_biometrics/biometric_data.csv
#   - ./data/raw_tracking_data/sunset_sunrise_data.csv
#   - ./data/lake_coords/lake_muddyfoot_polygon.gpkg
#
# Output:
#   - ./data/tracks_filtered/muddyfoot/01_muddyfoot_sub.rds
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
save_filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
biometrics_data_path <- "./data/fish_biometrics/"

cat("\014")  # Clear console

#==============================================================================
# 1. IMPORT AND INITIAL DATA PREPARATION
#==============================================================================

# Import raw detection data -------------------------------------------------
muddyfoot <- fread(paste0(raw_tracking_data_path, "muddyfoot_A/results/animal/all.csv"))
message("Imported ", nrow(muddyfoot), " raw detections") #13191157

# Standardize latitude/longitude column names -------------------------------
muddyfoot <- muddyfoot %>%
  rename(Lat = Latitude,
         Long = Longitude)

# Prepare timestamp columns -------------------------------------------------
# Original timestamps are in UTC; convert to local time (Europe/Stockholm)
muddyfoot <- muddyfoot %>%
  rename(timestamp_utc = Time) %>%
  mutate(
    timestamp_cest = with_tz(timestamp_utc, tzone = "Europe/Stockholm"),
    date_utc = as.Date(timestamp_utc),
    date_cest = as.Date(format(timestamp_cest, tz = "Europe/Stockholm"))
  )

# Verify timezone conversions
message("Non-aligned dates (UTC vs CEST): ",
        sum(muddyfoot$date_utc != muddyfoot$date_cest)) #1074230

#==============================================================================
# 2. MERGE WITH BIOMETRIC DATA
#==============================================================================

# Import fish biometric data ------------------------------------------------
biometrics <- fread(paste0(biometrics_data_path, "biometric_data.csv"))

# Extract individual ID from tag number and prepare for merging
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# Merge biometric data with detections --------------------------------------
muddyfoot <- muddyfoot %>%
  rename(individual_ID = Id) %>%
  left_join(biometrics, by = "individual_ID")

# Label reference tag and filter for Lake Muddyfoot individuals ------------
muddyfoot <- muddyfoot %>%
  mutate(individual_ID = if_else(FullId == 'H170-1802-65066',
                                 'Reference',
                                 as.character(individual_ID))) %>%
  filter(Lake == 'Muddyfoot' | individual_ID == 'Reference')

message("After lake filter: ", nrow(muddyfoot), " detections") #13182431

#==============================================================================
# 3. TEMPORAL FILTERING AND DATA CLEANING
#==============================================================================

# Filter for post-introduction period (after pike introduction) -------------
pike_intro_date <- as.Date("2022-09-25")

muddyfoot <- muddyfoot %>%
  filter(date_cest >= pike_intro_date) %>%
  # Remove unnecessary columns to reduce memory footprint
  dplyr::select(-c(ID, Notes, Transmitter, `Tag Number`, `PIT Number`, `Tag Type`,
            `Biologger Number`, Station, HPEm, TempData, DepthData, AccelData,
            RxDetected, nRxDetected, `Serial Number`, nRxUsed, RxUsed,
            timestamp_utc, date_utc, Date, Lake))

message("After temporal filter: ", nrow(muddyfoot), " detections") #11663216

# Standardize individual IDs ------------------------------------------------
# Add 'F' prefix to distinguish fish individuals from reference
muddyfoot <- muddyfoot %>%
  mutate(individual_ID = paste0("F", individual_ID))

#==============================================================================
# 4. CONVERT TO MOVEMENT OBJECT AND SPATIAL FILTERING
#==============================================================================

# Convert to move2 object for spatial analysis -----------------------------
muddyfoot_mv <- mt_as_move2(
  muddyfoot,
  coords = c("Long", "Lat"),
  crs = "WGS84",
  time_column = "timestamp_cest",
  track_id_column = "individual_ID",
  na.fail = FALSE  # Allow missing coordinates
)

# Sort chronologically by individual and timestamp --------------------------
muddyfoot_mv <- muddyfoot_mv %>%
  arrange(individual_ID, timestamp_cest)

# Load lake boundary polygon ------------------------------------------------
muddyfoot_polygon <- st_read("./data/lake_coords/polygons/muddyfoot_polygon.gpkg")

# Optional: Visualize reference tag to verify spatial accuracy
# ref_tag <- filter_track_data(muddyfoot_mv, .track_id = "FReference")
# muddyfoot_map <- mapview::mapView(mt_track_lines(ref_tag)$geometry)

# Note: To create or update the lake polygon, uncomment and run:
# muddyfoot_map <- mapedit::drawFeatures(map = muddyfoot_map)
# st_write(muddyfoot_map, dsn="data/lake_coords/lake_muddyfoot_polygon.gpkg",
#          driver="GPKG", delete_layer = TRUE)

# Filter detections within lake boundaries ---------------------------------
muddyfoot_sub <- st_filter(muddyfoot_mv, muddyfoot_polygon)
message("After spatial filter: ", nrow(muddyfoot_sub), " detections") #11098640

# Verify timezone preservation ----------------------------------------------
message("Timezone check: ", tz(muddyfoot_sub$timestamp_cest))

# Extract coordinates back to standard columns ------------------------------
coords <- st_coordinates(muddyfoot_sub$geometry)
muddyfoot_sub$Long <- coords[, 1]
muddyfoot_sub$Lat <- coords[, 2]

#==============================================================================
# 5. CREATE TEMPORAL CLASSIFICATION VARIABLES
#==============================================================================

# Classify experiment stage (Early vs Late) ---------------------------------
# Split dataset into two equal temporal periods
unique_dates <- unique(sort(muddyfoot_sub$date_cest))
n_dates <- length(unique_dates)
early_dates <- unique_dates[1:ceiling(n_dates/2)]

message("\nTotal unique dates: ", n_dates) #36
message("Early period dates: ", length(early_dates)) #18

muddyfoot_sub <- muddyfoot_sub %>%
  mutate(Stage = if_else(date_cest %in% early_dates, 'Early', 'Late'))

# Verify stage classification
message("\nStage distribution:")
print(table(muddyfoot_sub$Stage))
# Early    Late 
# 6356599 4742041 

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
muddyfoot_sub <- muddyfoot_sub %>%
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
print(table(muddyfoot_sub$time_of_day))
# Day   Night 
# 4925313 6173327

#==============================================================================
# 6. MERGE POST-EXPERIMENT SURVIVAL DATA
#==============================================================================

# Import post-experiment biometric data -------------------------------------
post_biometrics <- fread(paste0(biometrics_data_path, "biometric_post_exp_data.csv"))

post_size_cols <- post_biometrics %>%
  filter(Lake == 'Muddyfoot') %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number))) %>%
  dplyr::select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found,
         known_predated = Known_predated)

# Merge survival information ------------------------------------------------
muddyfoot_sub <- muddyfoot_sub %>%
  left_join(post_size_cols, by = "individual_ID")


#==============================================================================
# 7. FILTER MORTALITY EVENTS
#==============================================================================

# Note: Update individual ID if there are known early mortalities in Muddyfoot
# Remove individual that died on first day  -----------------
muddyfoot_sub <- muddyfoot_sub %>%
  filter(individual_ID != "F59707")
)

#==============================================================================
# 8. SAVE CLEANED DATASET
#==============================================================================

# Save as RDS file ----------------------------------------------------------
output_file <- paste0(filtered_data_path, "01_muddyfoot_sub.rds")
saveRDS(muddyfoot_sub, file = output_file)

message("Final dataset: ", nrow(muddyfoot_sub), " detections") #11098640 detections
message("Number of individuals: ", n_distinct(muddyfoot_sub$individual_ID)) #66 individuals (including reference tag)
message("Date range: ", min(muddyfoot_sub$date_cest), " to ",
        max(muddyfoot_sub$date_cest)) #Date range: 2022-09-25 to 2022-10-30 
message("Saved to: ", output_file)


#==============================================================================
#9. EXTRACT TRUE LAKE BOUNDARY COORDINATES
#==============================================================================

# Get coordinates from the polygon
coords <- st_coordinates(muddyfoot_polygon)

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
area_m2 <- st_area(muddyfoot_polygon)
area_hectares <- as.numeric(area_m2) / 10000
centroid <- st_centroid(muddyfoot_polygon)
centroid_coords <- st_coordinates(centroid)
bbox <- st_bbox(muddyfoot_polygon)

message("Lake area: ", round(area_hectares, 4), " hectares")
message("Centroid: Lon ", round(centroid_coords[1], 6), 
        ", Lat ", round(centroid_coords[2], 6))

# Save boundary coordinates as CSV
write.csv(lake_boundary_coords, "./data/lake_params/muddyfoot_boundary_coordinates.csv", 
          row.names = FALSE)

# Create summary object for use in other scripts
muddyfoot_lake_data <- list(
  boundary_coords = lake_boundary_coords,
  centroid = c(longitude = centroid_coords[1], latitude = centroid_coords[2]),
  area_hectares = as.numeric(area_hectares),
  area_m2 = as.numeric(area_m2),
  bbox = as.list(bbox),
  num_points = nrow(lake_boundary_coords)
)

# Save as RDS for easy loading in other scripts
saveRDS(muddyfoot_lake_data, "./data/lake_params/muddyfoot_coord_data.rds")
