### Exploring roach trajectories - Lake Cow Paradise ###

### SCRIPT DESCRIPITON ###

#In this script we explore the daily trajectories of roach in Lake cow
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

#filter for roach
cow_roach_dat <- cow_filt_data %>% 
  filter(Species == 'Roach') 

cow_roach_mv <- 
  mt_as_move2(cow_roach_dat, 
              coords = c("longitude", "latitude"),  # specify coordinate columns
              crs = "WGS84",  # use the WGS84 coordinate reference system
              time_column = "timestamp",  # specify the timestamp column
              track_id_column = "individual_ID",  # column identifying individual tracks
              na.fail = FALSE)  # allows rows with missing coordinates


cow_roach_mv <- cow_roach_mv %>%
  dplyr::arrange(individual_ID, timestamp)

cow_roach_mv <- 
  cow_roach_mv %>% 
  mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe
cow_roach_dat <- as.data.frame(cow_roach_mv)

#check time
tz(cow_roach_dat$timestamp) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(cow_roach_dat$geometry)
cow_roach_dat$Long <- coords[, 1]
cow_roach_dat$Lat <- coords[, 2]

# Add Date column
cow_roach_dat <- cow_roach_dat %>%
  mutate(Date = as.Date(timestamp))

#--------------------------------------------------------------------------------#

#### KNOWN PREDATED ####
#> F59826 ####

F59826_dat <- cow_roach_dat %>% 
  filter(individual_ID == 'F59826')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59826_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59826_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59826 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

#### LIKELY PREDATED ####

#> F59828 ####

F59828_dat <- cow_roach_dat %>% 
  filter(individual_ID == 'F59828')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59828_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59828_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59828 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

#### LIKELY MORTALITY ####

#> F59817 ####

F59817_dat <- cow_roach_dat %>% 
  filter(individual_ID == 'F59817')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59817_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59817_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59817 in cow Lake",
    x = "Longitude", y = "Latitude"
  )


#### UNCLEAR ####

#> F59820 ####

F59820_dat <- cow_roach_dat %>% 
  filter(individual_ID == 'F59820')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59820_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59820_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59820 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59831 ####

F59831_dat <- cow_roach_dat %>% 
  filter(individual_ID == 'F59831')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59831_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59831_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59831 in cow Lake",
    x = "Longitude", y = "Latitude"
  )


#### REMOVE POOR TRACKING ####

#> F59819 ####

F59819_dat <- cow_roach_dat %>% 
  filter(individual_ID == 'F59819')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59819_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59819_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59819 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

#remove individual due to poor tracking