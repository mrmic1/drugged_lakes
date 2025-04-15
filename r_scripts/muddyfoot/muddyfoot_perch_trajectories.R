### Exploring perch trajectories - Lake muddyfoot ###

### SCRIPT DESCRIPITON ###

#In this script we explore the daily trajectories of perch in Lake muddyfoot
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
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"  # Path for filtered tracking data.
telem_path <- "./data/telem_obj/muddyfoot/"
lake_polygon_path <- "./data/lake_coords/"

### Load data ###
muddyfoot_filt_data <- readRDS(paste0(filtered_data_path, '03_muddyfoot_sub.rds'))

# Load the polygon representing the boundary of muddyfoot lake, in the form of a GeoPackage file
muddyfoot_polygon <- sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))

#--------------------------------------------------------------------------------------------------#

#filter for perch
muddyfoot_perch_dat <- muddyfoot_filt_data %>% 
  filter(Species == 'Perch') 

muddyfoot_perch_mv <- 
  mt_as_move2(muddyfoot_perch_dat, 
              coords = c("longitude", "latitude"),  # specify coordinate columns
              crs = "WGS84",  # use the WGS84 coordinate reference system
              time_column = "timestamp",  # specify the timestamp column
              track_id_column = "individual_ID",  # column identifying individual tracks
              na.fail = FALSE)  # allows rows with missing coordinates


muddyfoot_perch_mv <- muddyfoot_perch_mv %>%
  dplyr::arrange(individual_ID, timestamp)

muddyfoot_perch_mv <- 
  muddyfoot_perch_mv %>% 
  mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe
muddyfoot_perch_dat <- as.data.frame(muddyfoot_perch_mv)

#check time
tz(muddyfoot_perch_dat$timestamp) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(muddyfoot_perch_dat$geometry)
muddyfoot_perch_dat$Long <- coords[, 1]
muddyfoot_perch_dat$Lat <- coords[, 2]

# Add Date column
muddyfoot_perch_dat <- muddyfoot_perch_dat %>%
  mutate(Date = as.Date(timestamp))

#--------------------------------------------------------------------------------#


#### LIKELY PREDATED ####

#> F59712 ####

F59712_dat <- muddyfoot_perch_dat %>% 
  filter(individual_ID == 'F59712')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59712_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59712_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59712 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#### LIKELY MORTALITY ####

#> F59686 ####

F59686_dat <- muddyfoot_perch_dat %>% 
  filter(individual_ID == 'F59686')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59686_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59686_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59686 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#### UNCLEAR ####

#> F59682 ####

F59682_dat <- muddyfoot_perch_dat %>% 
  filter(individual_ID == 'F59682')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59682_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59682_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59682 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59689 ####

F59689_dat <- muddyfoot_perch_dat %>% 
  filter(individual_ID == 'F59689')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59689_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59689_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59689 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59698 ####

F59698_dat <- muddyfoot_perch_dat %>% 
  filter(individual_ID == 'F59698')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59698_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59698_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59698 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59720 ####

F59720_dat <- muddyfoot_perch_dat %>% 
  filter(individual_ID == 'F59720')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59720_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59720_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59720 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59735 ####

F59735_dat <- muddyfoot_perch_dat %>% 
  filter(individual_ID == 'F59735')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59735_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59735_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59735 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )
