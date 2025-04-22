### Exploring perch trajectories ###

### SCRIPT DESCRIPITON ###

#In this script we explore the daily trajectories of perch in Lake BT
#that have irregular tracking data or are suspected of dieing during
#the experiment

#--------------------------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------------------#

#filter for perch
cow_perch_dat <- cow_filt_data %>% 
  filter(Species == 'Perch') 

cow_perch_mv <- 
  mt_as_move2(cow_perch_dat, 
              coords = c("longitude", "latitude"),  # specify coordinate columns
              crs = "WGS84",  # use the WGS84 coordinate reference system
              time_column = "timestamp",  # specify the timestamp column
              track_id_column = "individual_ID",  # column identifying individual tracks
              na.fail = FALSE)  # allows rows with missing coordinates


cow_perch_mv <- cow_perch_mv %>%
  dplyr::arrange(individual_ID, timestamp)

cow_perch_mv <- 
  cow_perch_mv %>% 
  mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe
cow_perch_dat <- as.data.frame(cow_perch_mv)

#check time
tz(cow_perch_dat$timestamp) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(cow_perch_dat$geometry)
cow_perch_dat$Long <- coords[, 1]
cow_perch_dat$Lat <- coords[, 2]

# Add Date column
cow_perch_dat <- cow_perch_dat %>%
  mutate(Date = as.Date(timestamp))

#---------------------------------------------------------------------------#

#### LIKELY MORTALITY ####

#> F59847 ####

F59847_dat <- cow_perch_dat %>% 
  filter(individual_ID == 'F59847')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59847_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59847_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59847 in cow Lake",
    x = "Longitude", y = "Latitude"
  )


#> F59849 ####

F59849_dat <- cow_perch_dat %>% 
  filter(individual_ID == 'F59849')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59849_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59849_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59849 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59852 ####

F59852_dat <- cow_perch_dat %>% 
  filter(individual_ID == 'F59852')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59852_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59852_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59852 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59853 ####

F59853_dat <- cow_perch_dat %>% 
  filter(individual_ID == 'F59853')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59853_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59853_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59853 in cow Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59878 ####

F59878_dat <- cow_perch_dat %>% 
  filter(individual_ID == 'F59878')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59878_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59878_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59878 in cow Lake",
    x = "Longitude", y = "Latitude"
  )


#### UNCLEAR ####

#### POOR TRACKING REMOVE ####

#> F59873 ####

F59873_dat <- cow_perch_dat %>% 
  filter(individual_ID == 'F59873')

# Plot
ggplot() +
  geom_sf(data = cow_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59873_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59873_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59873 in cow Lake",
    x = "Longitude", y = "Latitude"
  )




















