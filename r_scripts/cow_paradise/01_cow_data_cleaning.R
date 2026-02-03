#=============================================================================-
# Lake Cow Paradise Detection Data - Initial Cleaning and Preparation ####
#=============================================================================-
# Purpose: Clean and filter acoustic telemetry data from Lake cow_paradise
# Author: Marcus Michelangeli
#
# Input files:
#   - ./data/raw_tracking_data/cow_paradise_B/results/animal/all.csv
#   - ./data/fish_biometrics/biometric_data.csv
#   - ./data/raw_tracking_data/sunset_sunrise_data.csv
#   - ./data/lake_coords/lake_cow_paradise_polygon.gpkg
#
# Output:
#   - ./data/tracks_filtered/cow_paradise/01_cow_paradise_sub.rds
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
save_filtered_data_path <- "./data/tracks_filtered/cow_paradise/"
biometrics_data_path <- "./data/fish_biometrics/"

cat("\014")  # Clear console

#=============================================================================-
# 1. IMPORT AND INITIAL DATA PREPARATION ####
#=============================================================================-

# Import raw detection data -------------------------------------------------
cow_paradise <- fread(paste0(raw_tracking_data_path, "cow_paradise_C/results/animal/all.csv"))
message("Imported ", nrow(cow_paradise), " raw detections") #10368804 raw detections

# Fix swapped latitude/longitude columns ------------------------------------
cow_paradise <- cow_paradise %>%
  rename(Lat = Longitude, 
         Long = Latitude)

# Prepare timestamp columns -------------------------------------------------
# Original timestamps are in UTC; convert to local time (Europe/Stockholm)
cow_paradise <- cow_paradise %>%
  rename(timestamp_utc = Time) %>%
  mutate(
    timestamp_cest = with_tz(timestamp_utc, tzone = "Europe/Stockholm"),
    date_utc = as.Date(timestamp_utc),
    date_cest = as.Date(format(timestamp_cest, tz = "Europe/Stockholm"))
  )

# Verify timezone conversions
message("Non-aligned dates (UTC vs CEST): ",
        sum(cow_paradise$date_utc != cow_paradise$date_cest)) #794986

#=============================================================================-
# 2. MERGE WITH BIOMETRIC DATA ####
#=============================================================================-

# Import fish biometric data ------------------------------------------------
biometrics <- fread(paste0(biometrics_data_path, "biometric_data.csv"))

# Extract individual ID from tag number and prepare for merging
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# Merge biometric data with detections --------------------------------------
cow_paradise <- cow_paradise %>%
  rename(individual_ID = Id) %>%
  left_join(biometrics, by = "individual_ID")

# Label reference tag and filter for Lake cow_paradise individuals ------------
cow_paradise <- cow_paradise %>%
  mutate(individual_ID = if_else(FullId == 'H170-1802-65064',
                                 'Reference',
                                 as.character(individual_ID))) %>%
  filter(Lake == 'Cow Paradise' | individual_ID == 'Reference')

message("After lake filter: ", nrow(cow_paradise), " detections") #10364909 detections

#=============================================================================-
# 3. TEMPORAL FILTERING AND DATA CLEANING ####
#=============================================================================-

# Filter for post-introduction period (after pike introduction) -------------
pike_intro_date <- as.Date("2022-09-27")

cow_paradise <- cow_paradise %>%
  filter(date_cest >= pike_intro_date) %>%
  # Remove unnecessary columns to reduce memory footprint
  dplyr::select(-c(ID, Notes, Transmitter, `Tag Number`, `PIT Number`, `Tag Type`,
                   `Biologger Number`, Station, HPEm, TempData, DepthData, AccelData,
                   RxDetected, nRxDetected, `Serial Number`, nRxUsed, RxUsed,
                   timestamp_utc, date_utc, Date, Lake))

message("After temporal filter: ", nrow(cow_paradise), " detections") #9431255 detections

# Standardize individual IDs ------------------------------------------------
# Add 'F' prefix to distinguish fish individuals from reference
cow_paradise <- cow_paradise %>%
  mutate(individual_ID = paste0("F", individual_ID))

#=============================================================================-
# 4. CONVERT TO MOVEMENT OBJECT AND SPATIAL FILTERING ####
#=============================================================================-

# Convert to move2 object for spatial analysis -----------------------------
cow_paradise_mv <- mt_as_move2(
  cow_paradise,
  coords = c("Long", "Lat"),
  crs = "WGS84",
  time_column = "timestamp_cest",
  track_id_column = "individual_ID",
  na.fail = FALSE  # Allow missing coordinates
)

# Sort chronologically by individual and timestamp --------------------------
cow_paradise_mv <- cow_paradise_mv %>%
  arrange(individual_ID, timestamp_cest)

# Load lake boundary polygon ------------------------------------------------
cow_paradise_polygon <- st_read("./data/lake_params/polygons/cow_polygon.gpkg")

# Optional: Visualize reference tag to verify spatial accuracy
# ref_tag <- filter_track_data(cow_paradise_mv, .track_id = "FReference")
# cow_paradise_map <- mapview::mapView(mt_track_lines(ref_tag)$geometry)

# Note: To create or update the lake polygon, uncomment and run:
# cow_paradise_map <- mapedit::drawFeatures(map = cow_paradise_map)
# st_write(cow_paradise_map, dsn="data/lake_coords/lake_cow_paradise_polygon.gpkg",
#          driver="GPKG", delete_layer = TRUE)

# Filter detections within lake boundaries ---------------------------------
cow_paradise_sub <- st_filter(cow_paradise_mv, cow_paradise_polygon)
message("After spatial filter: ", nrow(cow_paradise_sub), " detections") #9089970 detections

# Verify timezone preservation ----------------------------------------------
message("Timezone check: ", tz(cow_paradise_sub$timestamp_cest))

# Extract coordinates back to standard columns ------------------------------
coords <- st_coordinates(cow_paradise_sub$geometry)
cow_paradise_sub$Long <- coords[, 1]
cow_paradise_sub$Lat <- coords[, 2]

#=============================================================================-
# 5. CREATE TEMPORAL CLASSIFICATION VARIABLES ####
#=============================================================================-

# Classify experiment stage (Early vs Late) ---------------------------------
# Split dataset into two equal temporal periods
unique_dates <- unique(sort(cow_paradise_sub$date_cest))
n_dates <- length(unique_dates)
early_dates <- unique_dates[1:ceiling(n_dates/2)]

message("\nTotal unique dates: ", n_dates) #35
message("Early period dates: ", length(early_dates)) #18

cow_paradise_sub <- cow_paradise_sub %>%
  mutate(Stage = if_else(date_cest %in% early_dates, 'Early', 'Late'))

# Verify stage classification
message("\nStage distribution:")
print(table(cow_paradise_sub$Stage))
#  Early    Late 
# 6718978 2370992  

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
cow_paradise_sub <- cow_paradise_sub %>%
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
print(table(cow_paradise_sub$time_of_day))
# Day   Night 
# 4388787 4699445 

#=============================================================================-
# 6. MERGE POST-EXPERIMENT SURVIVAL DATA ####
#=============================================================================-

# Import post-experiment biometric data -------------------------------------
post_biometrics <- fread(paste0(biometrics_data_path, "biometric_post_exp_data.csv"))

post_size_cols <- post_biometrics %>%
  filter(Lake == 'cow_paradise') %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number))) %>%
  dplyr::select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found,
         known_predated = Known_predated)

# Merge survival information ------------------------------------------------
cow_paradise_sub <- cow_paradise_sub %>%
  left_join(post_size_cols, by = "individual_ID")


#=============================================================================-
# 7. FILTER EARLY MORTALITY EVENTS ####
#=============================================================================-

# Note: Update individual ID if there are known early mortality in cow_paradise
# Remove individual that died on first day  -----------------
# 1 pike died on the first day, 2 fish (1 perch and 1 roach) had very poor tracking so are being removed
cow_paradise_sub <- cow_paradise_sub %>%
  filter(!individual_ID %in% c("F59819", "F59873", "F59893"))


#=============================================================================-
# 8. SAVE CLEANED DATASET ####
#=============================================================================-

# Save as RDS file ----------------------------------------------------------
output_file <- paste0(save_filtered_data_path, "01_cow_sub.rds")
saveRDS(cow_paradise_sub, file = output_file)

message("Final dataset: ", nrow(cow_paradise_sub), " detections") #9012169 detections
message("Number of individuals: ", n_distinct(cow_paradise_sub$individual_ID)) #64 individuals (including reference tag)
message("Date range: ", min(cow_paradise_sub$date_cest), " to ",
        max(cow_paradise_sub$date_cest)) #Date range: 2022-09-27 to 2022-10-31 
message("Saved to: ", output_file)


#=============================================================================-
#9. EXTRACT TRUE LAKE BOUNDARY COORDINATES ####
#=============================================================================-

# Get coordinates from the polygon
coords <- st_coordinates(cow_paradise_polygon)

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
area_m2 <- st_area(cow_paradise_polygon)
area_hectares <- as.numeric(area_m2) / 10000
centroid <- st_centroid(cow_paradise_polygon)
centroid_coords <- st_coordinates(centroid)
bbox <- st_bbox(cow_paradise_polygon)

message("Lake area: ", round(area_hectares, 4), " hectares")
message("Centroid: Lon ", round(centroid_coords[1], 6), 
        ", Lat ", round(centroid_coords[2], 6))

# Save boundary coordinates as CSV
write.csv(lake_boundary_coords, "./data/lake_params/cow_paradise_boundary_coordinates.csv", 
          row.names = FALSE)

# Create summary object for use in other scripts
cow_paradise_lake_data <- list(
  boundary_coords = lake_boundary_coords,
  centroid = c(longitude = centroid_coords[1], latitude = centroid_coords[2]),
  area_hectares = as.numeric(area_hectares),
  area_m2 = as.numeric(area_m2),
  bbox = as.list(bbox),
  num_points = nrow(lake_boundary_coords)
)

# Save as RDS for easy loading in other scripts
saveRDS(cow_paradise_lake_data, "./data/lake_params/cow_paradise_coord_data.rds")
