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

#filter for pike
BT_pike_dat <- BT_filt_data %>% 
  filter(Species == 'Pike') 

BT_pike_mv <- 
  mt_as_move2(BT_pike_dat, 
              coords = c("longitude", "latitude"),  # specify coordinate columns
              crs = "WGS84",  # use the WGS84 coordinate reference system
              time_column = "timestamp",  # specify the timestamp column
              track_id_column = "individual_ID",  # column identifying individual tracks
              na.fail = FALSE)  # allows rows with missing coordinates


BT_pike_mv <- BT_pike_mv %>%
  dplyr::arrange(individual_ID, timestamp)

BT_pike_mv <- 
  BT_pike_mv %>% 
  mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe
BT_pike_dat <- as.data.frame(BT_pike_mv)

#check time
tz(BT_pike_dat$timestamp) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(BT_pike_dat$geometry)
BT_pike_dat$Long <- coords[, 1]
BT_pike_dat$Lat <- coords[, 2]

# Add Date column
BT_pike_dat <- BT_pike_dat %>%
  mutate(Date = as.Date(timestamp))

#--------------------------------------------------------------------------------#

#> F59889 ####

#Tags were found outside of lake
#May have been taken by an otter

F59889_dat <- BT_pike_dat %>% 
  filter(individual_ID == 'F59889')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59889_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59889_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59889 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#Died right at the beginning of the study


#> F59892 ####

#Tags were found outside of lake
#May have been taken by an otter

F59892_dat <- BT_pike_dat %>% 
  filter(individual_ID == 'F59892')


# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59892_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59892_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59892 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#Died around the 30 September


#> F59886 ####

#This individual was not found

F59886_dat <- BT_pike_dat %>% 
  filter(individual_ID == 'F59886')


# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59886_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59886_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59886 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#Diet around 02 October


#> F59891 ####

#This individual was found, but didn't move much

F59891_dat <- BT_pike_dat %>% 
  filter(individual_ID == 'F59891')


# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59891_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59891_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59891 in BT Lake",
    x = "Longitude", y = "Latitude"
  )
