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
raw_tracking_data_path <- "./data/raw_tracking_data/"  # path to raw transmitter data
save_filtered_data_path <- "./data/tracks_filtered/muddyfoot/"  # path to save subsampled data
size_data_path <- "./data/fish_size/" 

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Load detection data for Lake Muddyfoot
muddyfoot <- fread(paste0(raw_tracking_data_path, "muddyfoot_A/results/animal/all.csv"))

# Rename latitude and longitude columns for consistency with other datasets
names(muddyfoot)[names(muddyfoot) == "Longitude"] <- "Long"
names(muddyfoot)[names(muddyfoot) == "Latitude"] <- "Lat"

#check time column
tz(muddyfoot$Time) #timestamps are in UTC
#rename time column to time_utc
muddyfoot <- muddyfoot %>% 
  rename(timestamp_utc = Time)

#add column for time in CEST
muddyfoot$timestamp_cest <- with_tz(muddyfoot$timestamp_utc, tzone = "Europe/Stockholm" )
tz(muddyfoot$timestamp_cest)
head(muddyfoot %>% 
       dplyr::select(timestamp_utc, timestamp_cest), n = 20)

#date column for both utc and cest
muddyfoot$date_utc <- as.Date(muddyfoot$timestamp_utc)
muddyfoot$date_cest <- as.Date(format(muddyfoot$timestamp_cest, tz = "Europe/Stockholm"))

#find rows where date_utc and date_cest do not align
print(muddyfoot %>% 
        filter(date_utc != date_cest))


# Load fish biometric data
biometrics <- fread(paste0(size_data_path, "biometric_data.csv"))

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
#original rows: 13191157
#after filtering: 13182431 

# Filter for observations after pike were introduced (post 2022-09-25) and remove unnecessary columns
muddyfoot <- muddyfoot %>%
  filter(date_cest >= "2022-09-25") %>%
  dplyr::select(-ID, -Notes, -Transmitter, -`Tag Number`, -`PIT Number`, - `Tag Type`, - `Biologger Number`,
                -Station, -HPEm, -TempData, -DepthData, -AccelData, 
                -RxDetected, -nRxDetected, -`Serial Number`, -nRxUsed, -RxUsed,
                -timestamp_utc, -date_utc, -Date, -Lake) # removing redundant columns
#original rows: 13182431
#after filtering cest: 11663216 
#after filtering utc: 11663216


#Change individual_ID to include F at the front
#This will make it easier for analysis later
muddyfoot$individual_ID <- paste0("F", muddyfoot$individual_ID)


# Convert muddyfoot to a move2 object for further spatial manipulation
muddyfoot_mv <- mt_as_move2(muddyfoot, 
                            coords = c("Long", "Lat"),  # specify coordinate columns
                            crs = "WGS84",  # use the WGS84 coordinate reference system
                            time_column = "timestamp_cest",  # specify the timestamp column
                            track_id_column = "individual_ID",  # column identifying individual tracks
                            na.fail = FALSE)  # allows rows with missing coordinates

# Sort the data by individual ID and timestamp for chronological consistency
# It is also important because it place individual_ID in numerical order which will be beneficial for 
# later parts of the analysis
muddyfoot_mv <- muddyfoot_mv %>%
  dplyr::arrange(individual_ID, timestamp_cest)

# Removing duplicates based on locations and timestamps (uncomment to activate if needed)
# muddyfoot_mv <- muddyfoot_mv %>% mt_filter_unique(criterion = "first")

#------------------------------------------------------------------------------#
# Handling spatial boundaries (polygon of Lake Muddyfoot) #
#------------------------------------------------------------------------------#

# Load the pre-drawn polygon representing Lake Muddyfoot's boundaries
muddyfoot_lake_poly <- sf::st_read("./data/lake_coords/lake_muddyfoot_polygon.gpkg")

# Subset the data to include only points that fall within the lake's polygon
muddyfoot_sub <- st_intersection(muddyfoot_mv, muddyfoot_lake_poly)

# Save the subsampled data for future use
# saveRDS(muddyfoot_sub, file = paste0(save_data_sub_path, "muddyfoot_sub_spatialfilt.rds"))

# Optional: Thin the data to reduce size by filtering points that occur at least 25 seconds apart
# muddyfoot_sub <- muddyfoot_sub %>% mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe

muddyfoot_sub <- as.data.frame(muddyfoot_sub)

#check time
tz(muddyfoot_sub$timestamp_cest) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(muddyfoot_sub$geometry)
muddyfoot_sub$Long <- coords[, 1]
muddyfoot_sub$Lat <- coords[, 2]

# # Convert the dataframe to an sf object with the WGS84 CRS
# data_sf <- st_as_sf(muddyfoot_sub, coords = c("Long", "Lat"), crs = 4326)
# 
# # Define the UTM CRS (using zone 34 for this specific location)
# utm_crs <- st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs")
# 
# # Transform the coordinates to UTM
# data_sf_utm <- st_transform(data_sf, crs = utm_crs)
# 
# # Extract the UTM coordinates
# utm_coords <- st_coordinates(data_sf_utm$geometry)
# muddyfoot_sub$Long_utm <- utm_coords[, 1]
# muddyfoot_sub$Lat_utm <- utm_coords[, 2]

#------------------------------------------------------------------------------#
# Create temporal columns #
#-------------------------------------------------------------------------------

#Split data in half by date
unique(muddyfoot_sub$date_cest) #36 unique dates
class(muddyfoot_sub$date_cest)

#Create column called stage to split dataset into 'early' and 'late'
muddyfoot_sub <- 
  muddyfoot_sub %>%
  # Find the unique dates, arrange them, and assign 'Early' or 'Late'
  mutate(Stage = if_else(date_cest %in% unique(sort(date_cest))[1:18], 'Early', 'Late'))

#check outcome
muddyfoot_sub %>% 
  group_by(date_cest, Stage) %>% 
  summarise(n = n()) %>% 
  print(n = 36) #worked

#Create day and night columns
#import sunset and sunrise data
sun_data <- fread(paste0(raw_tracking_data_path, "sunset_sunrise_data.csv"))
sun_data <- sun_data %>% 
  dplyr::select(Sunrise_timestamp, Sunset_timestamp) %>% 
  rename(sunrise = Sunrise_timestamp,
         sunset = Sunset_timestamp)

# Convert Sunrise_timestamp and Sunset_timestamp to POSIXct format with timezone
sun_data[, sunrise := dmy_hm(sunrise, tz = "Europe/Stockholm")]
sun_data[, sunset := dmy_hm(sunset, tz = "Europe/Stockholm")]

# Convert the timestamps in sun_data to POSIXct with the same timezone
sun_data$sunrise <- as.POSIXct(sun_data$sunrise, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Stockholm")
sun_data$sunset <- as.POSIXct(sun_data$sunset, format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Stockholm")

# Verify the timezone is consistent
# attr(muddyfoot_sub$timestamp_cest, "tzone")
# attr(sun_data$sunset, "tzone")

#create date column in sun_data
sun_data$date_cest <- as.Date(sun_data$sunset)

# Merge the sun_data with muddyfoot_sub to get the sunrise and sunset times for each fish location timestamp
muddyfoot_sub <- merge(muddyfoot_sub, sun_data, by = "date_cest", all.x = TRUE)

# Create the Day_time column based on the comparison of timestamp with sunrise and sunset
muddyfoot_sub$time_of_day <- ifelse(
  muddyfoot_sub$timestamp_cest >= muddyfoot_sub$sunrise & muddyfoot_sub$timestamp_cest <= muddyfoot_sub$sunset,
  "Day", "Night"
)

# Check entries classified as 'Day'
# day_entries <- muddyfoot_sub[muddyfoot_sub$time_of_day == "Day", ]
# head(day_entries[, c("timestamp_cest", "sunrise", "sunset", "time_of_day")], 20)


# Check entries classified as 'Night'
# night_entries <- muddyfoot_sub[muddyfoot_sub$time_of_day == "Night", ]
# head(night_entries[, c("timestamp_cest", "sunrise", "sunset", "time_of_day")], 10)

#------------------------------------------------------------------------------#
# Final save #
#------------------------------------------------------------------------------#

# Save the cleaned and filtered dataset
saveRDS(muddyfoot_sub, file = paste0(save_filtered_data_path, "01_muddyfoot_sub.rds"))
