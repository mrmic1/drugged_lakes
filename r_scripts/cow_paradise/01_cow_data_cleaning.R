#-------------------------------------------------------------#
# Initial cleaning of the detection data in Lake Cow Paradise 
#-------------------------------------------------------------#

# Load necessary libraries
library(data.table)
library(tidyverse)  # includes dplyr, ggplot2, and other utilities for data manipulation and visualization
library(move)       # useful for animal movement data analysis
library(move2)      # additional utilities for handling movement data
library(mapedit)    # for interactive map editing
library(sf)         # spatial data handling


#set time zones to CEST (not CET because study was conducted in the summmer)
#helpful timezone link: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
Sys.setenv(TZ = 'Europe/Stockholm')


# Define file paths
raw_tracking_data_path <- "./data/raw_tracking_data/"  # path to raw transmitter data
save_filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"  # path to save subsampled data
size_data_path <- "./data/fish_size/" 

cat("\014")

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

#import BT data
lake_cow = fread(paste0(raw_tracking_data_path, "cow_paradise_C/results/animal/all.csv"))

#noticed that lat and long are in the wrong columns, need to rename them
names(lake_cow)[names(lake_cow) == "Longitude"] <- "Lat"
names(lake_cow)[names(lake_cow) == "Latitude"] <- "Long"

#check time column
tz(lake_cow$Time) #timestamps are in UTC
#rename time column to time_utc
lake_cow <- lake_cow %>% 
  rename(timestamp_utc = Time)

#add column for time in CEST
lake_cow$timestamp_cest <- with_tz(lake_cow$timestamp_utc, tzone = "Europe/Stockholm" )
tz(lake_cow$timestamp_cest)
head(lake_cow %>% 
       dplyr::select(timestamp_utc, timestamp_cest), n = 20)

#date column for both utc and cest
lake_cow$date_utc <- as.Date(lake_cow$timestamp_utc)
lake_cow$date_cest <- as.Date(format(lake_cow$timestamp_cest, tz = "Europe/Stockholm"))

#find rows where date_utc and date_cest do not align
print(lake_cow %>% 
        filter(date_utc != date_cest))


#Load fish biometric data
biometrics <- fread(paste0(size_data_path, "biometric_data.csv"))

# Extract the last numeric part of 'Tag Number' to create a matching 'individual_ID' column
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# Rename the 'Id' column in lake_cow to 'individual_ID' for consistency in merging
names(lake_cow)[names(lake_cow) == "Id"] <- "individual_ID"

# Join biometric data with the detection data on 'individual_ID'
lake_cow <- left_join(lake_cow, biometrics, by = "individual_ID")

# Convert lake_cow to a data frame (optional but useful for further manipulation)
lake_cow <- as.data.frame(lake_cow)

# Identify and label the reference individual based on 'FullId' column
lake_cow$individual_ID <- ifelse(lake_cow$FullId == 'H170-1802-65064', 'Reference', lake_cow$individual_ID)

# Filter data to only include individuals within Lake lake_cow or the reference individual
lake_cow <- lake_cow %>%
  dplyr::filter(Lake == 'Cow Paradise' | individual_ID == 'Reference')
#original rows: 10368804
#after filtering: 10364909 

# Filter for observations after pike were introduced (post 2022-09-26) and remove unnecessary columns
lake_cow <- lake_cow %>%
  filter(date_cest >= "2022-09-27") %>%
  dplyr::select(-ID, -Notes, -Transmitter, -`Tag Number`, -`PIT Number`, - `Tag Type`, - `Biologger Number`,
                -Station, -HPEm, -TempData, -DepthData, -AccelData, 
                -RxDetected, -nRxDetected, -`Serial Number`, -nRxUsed, -RxUsed,
                -timestamp_utc, -date_utc, -Date, -Lake) # removing redundant columns
#original rows: 10364909
#after filtering cest: 9431255 

#Change individual_ID to include F at the front
#This will make it easier for analysis later
lake_cow$individual_ID <- paste0("F", lake_cow$individual_ID)


# Convert lake_cow to a move2 object for further spatial manipulation
lake_cow_mv <- mt_as_move2(lake_cow, 
                          coords = c("Long", "Lat"),  # specify coordinate columns
                          crs = "WGS84",  # use the WGS84 coordinate reference system
                          time_column = "timestamp_cest",  # specify the timestamp column
                          track_id_column = "individual_ID",  # column identifying individual tracks
                          na.fail = FALSE)  # allows rows with missing coordinates

# Sort the data by individual ID and timestamp for chronological consistency
# It is also important because it place individual_ID in numerical order which will be beneficial for 
# later parts of the analysis
lake_cow_mv <- lake_cow_mv %>%
  dplyr::arrange(individual_ID, timestamp_cest)


#------------------------------------------------------------------------------#

##------------------------------------------------------------------------------#
# Handling spatial boundaries (polygon of Lake Muddyfoot) #
#------------------------------------------------------------------------------#

# Load the pre-drawn polygon representing Lake lake_cow's boundaries
lake_cow_poly <- sf::st_read("./data/lake_coords/lake_cow_polygon.gpkg")


# Subset the data to include only points that fall within the lake's polygon
lake_cow_sub <- st_intersection(lake_cow_mv, lake_cow_poly)

#before filtering: 9431255
#after filtering: 9089970

#Convert to dataframe
lake_cow_sub <- as.data.frame(lake_cow_sub)

#check time
tz(lake_cow_sub$timestamp_cest) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(lake_cow_sub$geometry)
lake_cow_sub$Long <- coords[, 1]
lake_cow_sub$Lat <- coords[, 2]

# # Convert the dataframe to an sf object with the WGS84 CRS
# data_sf <- st_as_sf(lake_cow_sub, coords = c("Long", "Lat"), crs = 4326)
# 
# # Define the UTM CRS (using zone 34 for this specific location)
# utm_crs <- st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs")
# 
# # Transform the coordinates to UTM
# data_sf_utm <- st_transform(data_sf, crs = utm_crs)
# 
# # Extract the UTM coordinates
# utm_coords <- st_coordinates(data_sf_utm$geometry)
# lake_cow_sub$Long_utm <- utm_coords[, 1]
# lake_cow_sub$Lat_utm <- utm_coords[, 2]



#------------------------------------------------------------------------------#
# Create temporal columns #
#-------------------------------------------------------------------------------

#Split data in half by date
unique(lake_cow_sub$date_cest) #35 unique dates
class(lake_cow_sub$date_cest)

#Create column called stage to split dataset into 'early' and 'late'
lake_cow_sub <- 
  lake_cow_sub %>%
  # Find the unique dates, arrange them, and assign 'Early' or 'Late'
  mutate(Stage = if_else(date_cest %in% unique(sort(date_cest))[1:17], 'Early', 'Late'))

#check outcome
lake_cow_sub %>% 
  group_by(date_cest, Stage) %>% 
  summarise(n = n()) %>% 
  print(n = 35) #worked

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
# attr(lake_cow_sub$timestamp_cest, "tzone")
# attr(sun_data$sunset, "tzone")

#create date column in sun_data
sun_data$date_cest <- as.Date(sun_data$sunset)

# Merge the sun_data with lake_cow_sub to get the sunrise and sunset times for each fish location timestamp
lake_cow_sub <- merge(lake_cow_sub, sun_data, by = "date_cest", all.x = TRUE)

# Create the Day_time column based on the comparison of timestamp with sunrise and sunset
lake_cow_sub$time_of_day <- 
  ifelse(
    lake_cow_sub$timestamp_cest >= lake_cow_sub$sunrise & lake_cow_sub$timestamp_cest <= lake_cow_sub$sunset,
    "Day", "Night"
  )

# Check entries classified as 'Day'
# day_entries <- lake_cow_sub[lake_cow_sub$time_of_day == "Day", ]
# head(day_entries[, c("timestamp_cest", "sunrise", "sunset", "time_of_day")], 20)


# Check entries classified as 'Night'
# night_entries <- lake_cow_sub[lake_cow_sub$time_of_day == "Night", ]
# head(night_entries[, c("timestamp_cest", "sunrise", "sunset", "time_of_day")], 10)

#------------------------------------------------------------------------------#
# Final save #
#------------------------------------------------------------------------------#

# Save the cleaned and filtered dataset
saveRDS(lake_cow_sub, file = paste0(save_filtered_data_path, "01_lake_cow_sub.rds"))
