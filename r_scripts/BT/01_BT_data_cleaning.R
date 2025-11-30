#-------------------------------------------------------------#
# Initial cleaning of the detection data in Lake BT 
#-------------------------------------------------------------#

# Load necessary libraries
library(data.table)
library(tidyverse)  # includes dplyr, ggplot2, and other utilities for data manipulation and visualization
library(move)       # useful for animal movement data analysis
library(move2)      # additional utilities for handling movement data
library(mapedit)    # for interactive map editing
library(sf)         # spatial data handling

# Define file paths
raw_tracking_data_path <- "./data/raw_tracking_data/"  # path to raw transmitter data
save_filtered_data_path <- "./data/tracks_filtered/lake_BT/"  # path to save subsampled data
size_data_path <- "./data/fish_size/" 

cat("\014")

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

#import BT data
lake_BT = fread(paste0(raw_tracking_data_path, "BT_B/results/animal/all.csv"))

#noticed that lat and long are in the wrong columns, need to rename them
names(lake_BT)[names(lake_BT) == "Longitude"] <- "Lat"
names(lake_BT)[names(lake_BT) == "Latitude"] <- "Long"

#check time column
tz(lake_BT$Time) #timestamps are in UTC
#rename time column to time_utc
lake_BT <- lake_BT %>% 
  rename(timestamp_utc = Time)

#add column for time in CEST
lake_BT$timestamp_cest <- with_tz(lake_BT$timestamp_utc, tzone = "Europe/Stockholm" )
tz(lake_BT$timestamp_cest)

head(lake_BT %>% 
       dplyr::select(timestamp_utc, timestamp_cest), n = 20)

#date column for both utc and cest
lake_BT$date_utc <- as.Date(lake_BT$timestamp_utc)
#this next bit of code takes some time to run
lake_BT$date_cest <- as.Date(format(lake_BT$timestamp_cest, tz = "Europe/Stockholm"))


#find rows where date_utc and date_cest do not align
non_aligned_rows = lake_BT %>% 
        filter(date_utc != date_cest)
#1755991

#Load fish biometric data
biometrics <- fread(paste0(size_data_path, "biometric_data.csv"))

# Extract the last numeric part of 'Tag Number' to create a matching 'individual_ID' column
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# Rename the 'Id' column in lake_BT to 'individual_ID' for consistency in merging
names(lake_BT)[names(lake_BT) == "Id"] <- "individual_ID"

# Join biometric data with the detection data on 'individual_ID'
lake_BT <- left_join(lake_BT, biometrics, by = "individual_ID")

# Convert lake_BT to a data frame (optional but useful for further manipulation)
lake_BT <- as.data.frame(lake_BT)

# Identify and label the reference individual based on 'FullId' column
lake_BT$individual_ID <- ifelse(lake_BT$FullId == 'H170-1802-65065', 'Reference', lake_BT$individual_ID)

# Filter data to only include individuals within Lake lake_BT or the reference individual
lake_BT <- lake_BT %>%
  dplyr::filter(Lake == 'BT' | individual_ID == 'Reference')
#original rows: 22967436
#after filtering: 22957599 

# Filter for observations after pike were introduced (post 2022-09-26) and remove unnecessary columns
lake_BT <- lake_BT %>%
  filter(date_cest >= "2022-09-26") %>%
  dplyr::select(-ID, -Notes, -Transmitter, -`Tag Number`, -`PIT Number`, - `Tag Type`, - `Biologger Number`,
                -Station, -HPEm, -TempData, -DepthData, -AccelData, 
                -RxDetected, -nRxDetected, -`Serial Number`, -nRxUsed, -RxUsed,
                -timestamp_utc, -date_utc, -Date, -Lake) # removing redundant columns
#original rows: 22957599
#after filtering cest: 20832375 

#Change individual_ID to include F at the front
#This will make it easier for analysis later
lake_BT$individual_ID <- paste0("F", lake_BT$individual_ID)


# Convert lake_BT to a move2 object for further spatial manipulation
lake_BT_mv <- mt_as_move2(lake_BT, 
                            coords = c("Long", "Lat"),  # specify coordinate columns
                            crs = "WGS84",  # use the WGS84 coordinate reference system
                            time_column = "timestamp_cest",  # specify the timestamp column
                            track_id_column = "individual_ID",  # column identifying individual tracks
                            na.fail = FALSE)  # allows rows with missing coordinates

# Sort the data by individual ID and timestamp for chronological consistency
# It is also important because it place individual_ID in numerical order which will be beneficial for 
# later parts of the analysis
lake_BT_mv <- lake_BT_mv %>%
  dplyr::arrange(individual_ID, timestamp_cest)


##------------------------------------------------------------------------------#
# Handling spatial boundaries (polygon of Lake BT) #
#------------------------------------------------------------------------------#

fish1 <- filter_track_data(lake_BT_mv, .track_id = "FReference") #Reference tag
lake_BT_map <- mapview::mapView(mt_track_lines(fish1)$geometry) 
#good - for the reference tag no location fall outside of the lake, 
#and it also looks like there is not a lot of error

#drawing a polygon around the lake
#devtools::install_github("r-spatial/mapedit")
#library(mapedit)
BT_map <- mapedit::drawFeatures(map = lake_BT_map)
sf::st_write(BT_map, dsn="data/lake_coords/lake_BT_polygon.gpkg", driver="GPKG", delete_layer = TRUE) # for overwriting

# Load the pre-drawn polygon representing Lake lake_BT's boundaries
lake_BT_poly <- sf::st_read("./data/lake_coords/lake_BT_polygon.gpkg")

#check coordinates of bounding box
sf::st_bbox(lake_BT_poly)
#view polygon coordinates
sf::st_coordinates(lake_BT_poly)
sf::st_crs(lake_BT_poly)
sf::st_write(
  lake_BT_poly,
  "./data/lake_coords/lake_BT_polygon.kml",
  delete_dsn = TRUE
)

# Subset the data to include only points that fall within the lake's polygon
lake_BT_sub <- st_intersection(lake_BT_mv, lake_BT_poly)
#Original: 20832375
#New: 20709282


# Save the subsampled data for future use
saveRDS(lake_BT_sub, file = paste0(save_filtered_data_path, "lake_BT_sub_spatialfilt.rds"))
saveRDS(lake_BT_sub_sum, file = paste0(save_filtered_data_path, "lake_BT_sub_spatialfilt_sum.rds"))

# Calculate time differences between successive timestamps for each individual.
lake_BT_sub_sum <- 
  lake_BT_sub %>%
  group_by(individual_ID) %>% 
  mutate(time_diff = c(NA, diff(timestamp_cest)))  # First difference is NA.

# Round the time differences to 3 decimal places for clarity.
lake_BT_sub_sum$time_diff <- as.numeric(round(lake_BT_sub_sum$time_diff, digits = 3))

# Check the first 20 rows to ensure time differences are calculated correctly.
head(lake_BT_sub_sum %>% dplyr::select(timestamp_cest, time_diff), n = 20)

# Calculate mean and median time differences per individual.
lake_BT_sub_sum <- 
  lake_BT_sub_sum %>%
  group_by(individual_ID) %>% 
  mutate(mean_time_diff = mean(time_diff, na.rm = TRUE),
         median_time_diff = median(time_diff, na.rm = TRUE)) %>% 
  ungroup()

# Count the number of positions (i.e., tracking locations) for each individual.
lake_BT_sub_sum <- 
  lake_BT_sub_sum %>%
  group_by(individual_ID) %>% 
  mutate(n_positions = n()) %>% 
  ungroup()

# Calculate the number of unique days each individual was monitored.
lake_BT_sub_sum <- 
  lake_BT_sub_sum %>%
  group_by(individual_ID) %>% 
  mutate(n_days_tracked = length(unique(date_cest))) %>% 
  ungroup()

# Calculate the number of positions per day and median time differences between positions for each individual.
lake_BT_sub_sum <- 
  lake_BT_sub_sum %>%
  group_by(individual_ID, date_cest) %>% 
  mutate(n_positions_day = n(),
         median_day_time_diff = median(time_diff, na.rm = TRUE)) %>%  # Median time difference per day.
  ungroup()

# Calculate the number of positions per hour and per minute.
lake_BT_sub_sum <- 
  lake_BT_sub_sum %>%
  group_by(individual_ID, date_cest) %>% 
  mutate(n_positions_hourly = n_positions_day / 24,  # Average number of positions per hour.
         n_positions_per_min = n_positions_hourly / 60) %>%  # Average number of positions per minute.
  ungroup()

# Optional: Check the summary of number of days and positions per individual.
BT_filt_data %>% 
        dplyr::select(individual_ID, Treatment, n_positions, n_days_tracked) %>% 
        distinct()



# Optional: Thin the data to reduce size by filtering points that occur at least 25 seconds apart
lake_BT_sub <- lake_BT_sub %>% mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe

lake_BT_sub <- as.data.frame(lake_BT_sub)

#check time
tz(lake_BT_sub$timestamp_cest) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(lake_BT_sub$geometry)
lake_BT_sub$Long <- coords[, 1]
lake_BT_sub$Lat <- coords[, 2]

# # Convert the dataframe to an sf object with the WGS84 CRS
# data_sf <- st_as_sf(lake_BT_sub, coords = c("Long", "Lat"), crs = 4326)
# 
# # Define the UTM CRS (using zone 34 for this specific location)
# utm_crs <- st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs")
# 
# # Transform the coordinates to UTM
# data_sf_utm <- st_transform(data_sf, crs = utm_crs)
# 
# # Extract the UTM coordinates
# utm_coords <- st_coordinates(data_sf_utm$geometry)
# lake_BT_sub$Long_utm <- utm_coords[, 1]
# lake_BT_sub$Lat_utm <- utm_coords[, 2]

#------------------------------------------------------------------------------#
# Create temporal columns #
#-------------------------------------------------------------------------------

#Split data in half by date
unique(lake_BT_sub$date_cest) #34 unique dates
class(lake_BT_sub$date_cest)

#Create column called stage to split dataset into 'early' and 'late'
lake_BT_sub <- 
  lake_BT_sub %>%
  # Find the unique dates, arrange them, and assign 'Early' or 'Late'
  mutate(Stage = if_else(date_cest %in% unique(sort(date_cest))[1:17], 'Early', 'Late'))

#check outcome
lake_BT_sub %>% 
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
# attr(lake_BT_sub$timestamp_cest, "tzone")
# attr(sun_data$sunset, "tzone")

#create date column in sun_data
sun_data$date_cest <- as.Date(sun_data$sunset)

# Merge the sun_data with lake_BT_sub to get the sunrise and sunset times for each fish location timestamp
lake_BT_sub <- merge(lake_BT_sub, sun_data, by = "date_cest", all.x = TRUE)

# Create the Day_time column based on the comparison of timestamp with sunrise and sunset
lake_BT_sub$time_of_day <- 
  ifelse(
  lake_BT_sub$timestamp_cest >= lake_BT_sub$sunrise & lake_BT_sub$timestamp_cest <= lake_BT_sub$sunset,
  "Day", "Night"
)

# Check entries classified as 'Day'
# day_entries <- lake_BT_sub[lake_BT_sub$time_of_day == "Day", ]
# head(day_entries[, c("timestamp_cest", "sunrise", "sunset", "time_of_day")], 20)


# Check entries classified as 'Night'
# night_entries <- lake_BT_sub[lake_BT_sub$time_of_day == "Night", ]
# head(night_entries[, c("timestamp_cest", "sunrise", "sunset", "time_of_day")], 10)

#------------------------------------------------------------------------------#
# Final save #
#------------------------------------------------------------------------------#

# Save the cleaned and filtered dataset
saveRDS(lake_BT_sub, file = paste0(save_filtered_data_path, "01_lake_BT_sub.rds"))
