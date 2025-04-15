### Exploring roach trajectories - Lake muddyfoot ###

### SCRIPT DESCRIPITON ###

#In this script we explore the daily trajectories of roach in Lake muddyfoot
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

#filter for roach
muddyfoot_roach_dat <- muddyfoot_filt_data %>% 
  filter(Species == 'Roach') 

muddyfoot_roach_mv <- 
  mt_as_move2(muddyfoot_roach_dat, 
              coords = c("longitude", "latitude"),  # specify coordinate columns
              crs = "WGS84",  # use the WGS84 coordinate reference system
              time_column = "timestamp",  # specify the timestamp column
              track_id_column = "individual_ID",  # column identifying individual tracks
              na.fail = FALSE)  # allows rows with missing coordinates


muddyfoot_roach_mv <- muddyfoot_roach_mv %>%
  dplyr::arrange(individual_ID, timestamp)

muddyfoot_roach_mv <- 
  muddyfoot_roach_mv %>% 
  mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe
muddyfoot_roach_dat <- as.data.frame(muddyfoot_roach_mv)

#check time
tz(muddyfoot_roach_dat$timestamp) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(muddyfoot_roach_dat$geometry)
muddyfoot_roach_dat$Long <- coords[, 1]
muddyfoot_roach_dat$Lat <- coords[, 2]

# Add Date column
muddyfoot_roach_dat <- muddyfoot_roach_dat %>%
  mutate(Date = as.Date(timestamp))

#--------------------------------------------------------------------------------#

### LIKELY PREDATED ###

#> F59701 ####

F59701_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59701')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59701_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59701_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59701 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59719 ####

F59719_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59719')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59719_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59719_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59719 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#movement looks like it changes around 03-10-22

#> F59729 ####

F59729_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59729')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59729_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59729_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59729 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#movement looks like it changes around 18-10-22

#> F59731 ####

F59731_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59731')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59731_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59731_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59731 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#> F59738 ####

F59738_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59738')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59738_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59738_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59738 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )


### LIKELY MORTALITY ###

#> F59683 ####

F59683_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59683')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59683_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59683_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59683 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#Movement clearly changed around 29-09-22

#> F59687 ####

F59687_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59687')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59687_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59687_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59687 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#Movement clearly changed around 28-09-22

#> F59707 ####

F59707_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59707')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59707_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59707_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59707 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#only has 2 days of movement
#died within the first days of the experiment (remove)

#> F59745 ####

#Fish was not found
#Likely predated

F59745_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59745')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59745_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59745_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59745 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#weird movement around the 19th Oct

### UNCLEAR ###

#> F59709 ####

F59709_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59709')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59709_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59709_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59709 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#not clear


#> F59742 ####

#Fish was not found
#Likely predated

F59742_dat <- muddyfoot_roach_dat %>% 
  filter(individual_ID == 'F59742')

# Plot
ggplot() +
  geom_sf(data = muddyfoot_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59742_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59742_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59742 in muddyfoot Lake",
    x = "Longitude", y = "Latitude"
  )

#weird movement around the 19th Oct
