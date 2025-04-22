### Exploring pike trajectories - Lake Cow Paradise ###

### SCRIPT DESCRIPITON ###

#In this script we explore the daily trajectories of pike in Lake Cow Paradise
#that have irregular tracking data or are suspected of dieing during
#the experiment

#--------------------------------------------------------------------------------------------#

### LIBRARIES ###

library(ctmm)
library(tidyverse)
library(move2)
library(sf)         # spatial data handling


### DIRECTORIES ###
# Define paths to various datasets and save locations for tables.
filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"  # Path for filtered tracking data.
telem_path <- "./data/telem_obj/cow_paradise/"
lake_polygon_path <- "./data/lake_coords/"

### Load data ###
cow_filt_data <- readRDS(paste0(filtered_data_path, '03_lake_cow_sub.rds'))

# Load the polygon representing the boundary of cow lake, in the form of a GeoPackage file
cow_polygon <- sf::st_read(paste0(lake_polygon_path, "lake_cow_polygon.gpkg"))

#--------------------------------------------------------------------------------------------------#

#filter for pike
cow_pike_dat <- cow_filt_data %>% 
  filter(Species == 'Pike') 

cow_pike_mv <- 
  mt_as_move2(cow_pike_dat, 
              coords = c("longitude", "latitude"),  # specify coordinate columns
              crs = "WGS84",  # use the WGS84 coordinate reference system
              time_column = "timestamp",  # specify the timestamp column
              track_id_column = "individual_ID",  # column identifying individual tracks
              na.fail = FALSE)  # allows rows with missing coordinates


cow_pike_mv <- cow_pike_mv %>%
  dplyr::arrange(individual_ID, timestamp)

cow_pike_mv <- 
  cow_pike_mv %>% 
  mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe
cow_pike_dat <- as.data.frame(cow_pike_mv)

#check time
tz(cow_pike_dat$timestamp) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(cow_pike_dat$geometry)
cow_pike_dat$Long <- coords[, 1]
cow_pike_dat$Lat <- coords[, 2]

# Add Date column
cow_pike_dat <- cow_pike_dat %>%
  mutate(Date = as.Date(timestamp))

#--------------------------------------------------------------------------------#

#### MORTALITY ####

#> F59893 ####

F59893_dat <- cow_pike_dat %>% 
  filter(individual_ID == 'F59893')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59893_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59893_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59893 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59896 ####

F59896_dat <- cow_pike_dat %>% 
  filter(individual_ID == 'F59896')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59896_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59896_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59896 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

