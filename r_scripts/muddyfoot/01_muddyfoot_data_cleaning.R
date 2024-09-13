#-------------------------------------------------------------#
# Initial cleaning of the detection data in Lake Muddyfoot ####
#-------------------------------------------------------------#

# Load necessary libraries
library(data.table)
library(tidyverse)  # includes dplyr, ggplot2, and other utilities for data manipulation and visualization
library(move)       # useful for animal movement data analysis
library(move2)      # additional utilities for handling movement data
library(mapedit)    # for interactive map editing
library(sf)         # spatial data handling
library(ctmm)       # continuous-time movement models, useful for telemetry analysis

# Set the time zone environment to 'Europe/Stockholm' for consistent timestamp manipulation
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths
data_trans_path <- "./data/Transmitters/raw/"  # path to raw transmitter data
save_data_sub_path <- "./data/Transmitters/subsamples/"  # path to save subsampled data
save_data_clean_path <- "./data/Transmitters/clean/"  # path to save cleaned data


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Load detection data for Lake Muddyfoot
muddyfoot <- fread(paste0(data_trans_path, "muddyfoot_A/results/animal/all.csv"))

# Rename latitude and longitude columns for consistency with other datasets
names(muddyfoot)[names(muddyfoot) == "Longitude"] <- "Long"
names(muddyfoot)[names(muddyfoot) == "Latitude"] <- "Lat"

# Load biometric data
biometrics <- fread(paste0(data_trans_path, "biometric_data.csv"))

# Extract the last numeric part of 'Tag Number' to create a matching 'individual_ID' column
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# Rename the 'Id' column in muddyfoot to 'individual_ID' for consistency in merging
names(muddyfoot)[names(muddyfoot) == "Id"] <- "individual_ID"

# Join biometric data with the detection data on 'individual_ID'
muddyfoot <- left_join(muddyfoot, biometrics, by = "individual_ID")

# Convert muddyfoot to a data frame (optional but useful for further manipulation)
muddyfoot <- as.data.frame(muddyfoot)

# Identify and label the reference individual based on 'FullId' column
muddyfoot$individual_ID <- ifelse(muddyfoot$FullId == 'H170-1802-65066', 'Reference', muddyfoot$individual_ID)

# Filter data to only include individuals within Lake Muddyfoot or the reference individual
muddyfoot <- muddyfoot %>%
  dplyr::filter(Lake == 'Muddyfoot' | individual_ID == 'Reference')

# Filter for observations after pike were introduced (post 2022-09-25) and remove unnecessary columns
muddyfoot <- muddyfoot %>%
  filter(Time >= "2022-09-25") %>%
  dplyr::select(-ID, -Notes, -Transmitter, -`Tag Number`)  # removing redundant columns

# Set the timestamp to the correct time zone (converting from UTC to CEST)
muddyfoot$timestamp <- with_tz(muddyfoot$Time, "Europe/Stockholm")

# Convert muddyfoot to a move2 object for further spatial manipulation
muddyfoot_mv <- mt_as_move2(muddyfoot, 
                            coords = c("Long", "Lat"),  # specify coordinate columns
                            crs = "WGS84",  # use the WGS84 coordinate reference system
                            time_column = "timestamp",  # specify the timestamp column
                            track_id_column = "individual_ID",  # column identifying individual tracks
                            na.fail = FALSE)  # allows rows with missing coordinates

# Sort the data by individual ID and timestamp for chronological consistency
muddyfoot_mv <- muddyfoot_mv %>%
  dplyr::arrange(individual_ID, timestamp)

# Removing duplicates based on locations and timestamps (uncomment to activate if needed)
# muddyfoot_mv <- muddyfoot_mv %>% mt_filter_unique(criterion = "first")

#------------------------------------------------------------------------------#
# Handling spatial boundaries (polygon of Lake Muddyfoot) #
#------------------------------------------------------------------------------#

# Load the pre-drawn polygon representing Lake Muddyfoot's boundaries
muddyfoot_lake_poly <- sf::st_read("./data/Lakes/lake_muddyfoot_polygon.gpkg")

# Subset the data to include only points that fall within the lake's polygon
muddyfoot_sub <- st_intersection(muddyfoot_mv, muddyfoot_lake_poly)

# Save the subsampled data for future use
# saveRDS(muddyfoot_sub, file = paste0(save_data_sub_path, "muddyfoot_sub_spatialfilt.rds"))

# Optional: Thin the data to reduce size by filtering points that occur at least 25 seconds apart
# muddyfoot_sub <- muddyfoot_sub %>% mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#------------------------------------------------------------------------------#
# Time difference calculation and ID creation for daily tracking #
#------------------------------------------------------------------------------#

# Calculate the time difference (in seconds) between consecutive detections for each individual
muddyfoot_sub <- muddyfoot_sub %>%
  group_by(individual_ID) %>%
  mutate(time_diff = c(NA, diff(timestamp)))  # calculate time difference

# Round the time difference to 3 decimal places for clarity
muddyfoot_sub$time_diff <- round(as.numeric(muddyfoot_sub$time_diff), 3)

# Create additional columns to track the date, week, and day of the year
muddyfoot_sub$date <- strftime(muddyfoot_sub$timestamp, format = "%Y/%m/%d")
muddyfoot_sub$week <- strftime(muddyfoot_sub$timestamp, format = "%W")
muddyfoot_sub$day <- strftime(muddyfoot_sub$timestamp, format = "%j")

# Create a unique ID for each individual on each day
muddyfoot_sub <- muddyfoot_sub %>%
  mutate(individual_day = paste(individual_ID, day, sep = "_"))

#------------------------------------------------------------------------------#
# General dataframe cleaning #
#------------------------------------------------------------------------------#

# Identify cases where the number of receivers used differs from the number detected
inconsistent_rx <- muddyfoot_sub %>%
  filter(RxDetected != RxUsed)

# Remove unnecessary columns related to metadata
muddyfoot_sub <- muddyfoot_sub %>%
  dplyr::select(-Station, -HPEm, -TempData, -DepthData, -AccelData, 
                -Tag.Type, -Serial.Number, -PIT.Number, -Biologger.Number, 
                -RxDetected, -nRxDetected)

# Rename certain columns for clarity
names(muddyfoot_sub)[names(muddyfoot_sub) == "RxUsed"] <- "receivers_used"
names(muddyfoot_sub)[names(muddyfoot_sub) == "nRxUsed"] <- "num_receivers_used"

#------------------------------------------------------------------------------#
# Coordinate transformations (WGS84 to UTM) #
#------------------------------------------------------------------------------#

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(muddyfoot_sub$geometry)
muddyfoot_sub$Long <- coords[, 1]
muddyfoot_sub$Lat <- coords[, 2]

# Convert the dataframe to an sf object with the WGS84 CRS
data_sf <- st_as_sf(muddyfoot_sub, coords = c("Long", "Lat"), crs = 4326)

# Define the UTM CRS (using zone 34 for this specific location)
utm_crs <- st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs")

# Transform the coordinates to UTM
data_sf_utm <- st_transform(data_sf, crs = utm_crs)

# Extract the UTM coordinates
utm_coords <- st_coordinates(data_sf_utm$geometry)
muddyfoot_sub$Long_utm <- utm_coords[, 1]
muddyfoot_sub$Lat_utm <- utm_coords[, 2]

#------------------------------------------------------------------------------#
# Final save #
#------------------------------------------------------------------------------#

# Save the cleaned and filtered dataset
saveRDS(muddyfoot_sub, file = paste0(save_data_sub_path, "muddyfoot_sub.rds"))