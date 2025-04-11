### Exploring roach trajectories - Lake BT ###

### SCRIPT DESCRIPITON ###

#In this script we explore the daily trajectories of roach in Lake BT
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
filtered_data_path <- "./data/tracks_filtered/lake_BT/"  # Path for filtered tracking data.
telem_path <- "./data/telem_obj/BT/"
lake_polygon_path <- "./data/lake_coords/"

### Load data ###
BT_filt_data <- readRDS(paste0(filtered_data_path, '03_lake_BT_sub.rds'))

# Load the polygon representing the boundary of BT lake, in the form of a GeoPackage file
BT_polygon <- sf::st_read(paste0(lake_polygon_path, "lake_BT_polygon.gpkg"))

#--------------------------------------------------------------------------------------------------#

#filter for roach
BT_roach_dat <- BT_filt_data %>% 
  filter(Species == 'Roach') 

BT_roach_mv <- 
  mt_as_move2(BT_roach_dat, 
              coords = c("longitude", "latitude"),  # specify coordinate columns
              crs = "WGS84",  # use the WGS84 coordinate reference system
              time_column = "timestamp",  # specify the timestamp column
              track_id_column = "individual_ID",  # column identifying individual tracks
              na.fail = FALSE)  # allows rows with missing coordinates


BT_roach_mv <- BT_roach_mv %>%
  dplyr::arrange(individual_ID, timestamp)

BT_roach_mv <- 
  BT_roach_mv %>% 
  mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe
BT_roach_dat <- as.data.frame(BT_roach_mv)

#check time
tz(BT_roach_dat$timestamp) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(BT_roach_dat$geometry)
BT_roach_dat$Long <- coords[, 1]
BT_roach_dat$Lat <- coords[, 2]

# Add Date column
BT_roach_dat <- BT_roach_dat %>%
  mutate(Date = as.Date(timestamp))

#--------------------------------------------------------------------------------#

#> F59803 ####

F59803_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59803')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59803_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59803_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59803 in BT Lake",
    x = "Longitude", y = "Latitude"
  )



#> F59783 ####

F59783_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59783')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59783_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59783_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59783 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#Movement clearly changed around the 17-Oct or 18-Oct. Indicative of mortality

#> F59776 ####

F59776_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59776')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59776_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59776_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59776 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#not clear

#> F59777 ####

F59777_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59777')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59777_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59777_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59777 in BT Lake",
    x = "Longitude", y = "Latitude"
  )



#> F59807 ####

F59807_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59807')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59807_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59807_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59807 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#not clear


#> F59809 ####

#Fish was not found
#Likely predated

F59809_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59809')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59809_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59809_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59809 in BT Lake",
    x = "Longitude", y = "Latitude"
  )


#> F59810 ####

#Fish was not found
#Likely predated

F59810_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59810')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59810_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59810_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59810 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#Tracks pretty much stop after 10-04


#> F59812 ####

#Fish was not found
#Likely predated

F59812_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59812')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59812_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59812_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59812 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#Movement looks normal throughout the study.

#> F59815 ####

#Fish was not found
#Likely predated

F59815_dat <- BT_roach_dat %>% 
  filter(individual_ID == 'F59815')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59815_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59815_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59815 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#Movement looks normal throughout the study.
