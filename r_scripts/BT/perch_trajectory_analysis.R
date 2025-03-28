### Exploring perch trajectories ###

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
BT_filt_data <- readRDS(paste0(filtered_data_path, '04_lake_BT_sub.rds'))

# Load the polygon representing the boundary of BT lake, in the form of a GeoPackage file
BT_polygon <- sf::st_read(paste0(lake_polygon_path, "lake_BT_polygon.gpkg"))

#filter for perch
BT_perch_dat <- BT_filt_data %>% 
  filter(Species == 'Perch') 

BT_perch_mv <- 
  mt_as_move2(BT_perch_dat, 
              coords = c("longitude", "latitude"),  # specify coordinate columns
              crs = "WGS84",  # use the WGS84 coordinate reference system
              time_column = "timestamp",  # specify the timestamp column
              track_id_column = "individual_ID",  # column identifying individual tracks
              na.fail = FALSE)  # allows rows with missing coordinates


BT_perch_mv <- BT_perch_mv %>%
  dplyr::arrange(individual_ID, timestamp)

BT_perch_mv <- 
  BT_perch_mv %>% 
  mt_filter_per_interval(unit = "25 seconds", criterion = "first")

#Convert to dataframe
BT_perch_dat <- as.data.frame(BT_perch_mv)

#check time
tz(BT_perch_dat$timestamp) #Europe/Stockholm

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(BT_perch_dat$geometry)
BT_perch_dat$Long <- coords[, 1]
BT_perch_dat$Lat <- coords[, 2]

# Add Date column
BT_perch_dat <- BT_perch_dat %>%
  mutate(Date = as.Date(timestamp))



### F59801 ###

#Fish was not found
#Has a very small home range estimate relative to the population

F59801_dat <- BT_perch_dat %>% 
  filter(individual_ID == 'F59801')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59801_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59801_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59889 in BT Lake",
    x = "Longitude", y = "Latitude"
  )

#not clear

### F59764 ###

#Tags were found outside of lake
#May have been taken by an otter

F59764_dat <- BT_perch_dat %>% 
  filter(individual_ID == 'F59764')

# Plot
ggplot() +
  geom_sf(data = BT_polygon, fill = "lightblue", color = "black", alpha = 0.3) +  # lake polygon
  geom_path(data = F59764_dat, aes(x = Long, y = Lat, group = Date), color = "blue") +  # trajectory
  geom_point(data = F59764_dat, aes(x = Long, y = Lat), color = "black", size = 0.5) +  # points
  facet_wrap(~ Date) +  # one panel per date
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Daily Movement Trajectories of Individual F59889 in BT Lake",
    x = "Longitude", y = "Latitude"
  )


