# -----------------------------------# 
# Exploring interspecific interactions 
# -----------------------------------#

### LIBRARIES ###
library(dplyr)
library(move)
library(move2)


### DIRECTORIES ###
data_filter_path = "./data/tracks_filtered/"
lake_polygon_path = "./data/lake_coords/"
ctmm_path = "./data/ctmm_fits/"
telem_path = "./data/telem_obj/"


### LOAD DATA ###
muddyfoot_dat <- readRDS(paste0(data_filter_path, 'muddyfoot_final_filt_data.rds'))

### Extract individuals with irregular tracking ###
irreg_ind <- 
  muddyfoot_dat %>% 
  filter(irregular_sampling == 1)

#how many individuals
length(unique(irreg_ind$individual_id))
#19

#Explore F59707
#It was only tracked for 1 day
F59707_track = 
  irreg_ind %>% 
  filter(individual_id == 'F59707')


str(F59707_track)
class(F59707_track)

# now create a move2 object from a data.frame
F59707_track_mv <- 
  mt_as_move2(F59707_track, 
              coords = c("longitude","latitude"),
              crs = "WGS84",
              time_column = "timestamp",
              track_id_column = "individual_id",
              na.fail = F) # allows or not empty coordinates


#Plot tracks over outline of lake
#Load map polygon
muddyfoot_polygon = sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))

ggplot() + 
  coord_sf(crs = st_crs(F59707_track_mv)) + #set the coordinate system to that of the bats data
  geom_sf(data=muddyfoot_polygon) + 
  geom_sf(data=F59707_track_mv, col="red", alpha = 0.5)


#add scale bar
library(sf)
library(ggsn)
library(gganimate)

# Check the CRS of your lake polygon
st_crs(muddyfoot_polygon)

# Define the length of the scale bar in meters
scale_length <- 100  # 1 km

# Get the bounding box of the lake to determine where to place the scale bar
lake_bbox <- st_bbox(muddyfoot_polygon)

# Define the starting point of the scale bar near the bottom left corner of the map
x_start <- lake_bbox["xmin"] + (lake_bbox["xmax"] - lake_bbox["xmin"]) * 0.1
y_start <- lake_bbox["ymin"] + (lake_bbox["ymax"] - lake_bbox["ymin"]) * 0.1

# Define the ending point of the scale bar
x_end <- x_start + scale_length

p = ggplot() + 
  geom_sf(data=muddyfoot_polygon) + 
  geom_sf(data=F59707_track_mv, col="red", alpha = 0.5) + 
  ggsn::scalebar(
    data = muddyfoot_polygon, 
    dist = 5,  # Distance represented by the scale bar
    dist_unit = "m",  # Units of measurement
    transform = TRUE,  # Transform coordinates to projected CRS if needed
    model = 'WGS84',  # Model used for transformation, if needed
    location = "bottomright",  # Position of the scale bar
    st.dist = 0.05  # Distance of the text from the scale bar
  ) +
  theme_minimal() +
  transition_time(F59707_track_mv$timestamp) +  # Animate over the timestamp
  labs(title = "Fish Movement Over Time: {frame_time}") +  # Title that updates with time
  ease_aes('linear')

animate(p, nframes = 500, fps = 10, width = 800, height = 600)
 

#Later develop interactive plot to illustrate potential predation event. Plot predators for this day on top.  
#Also want to plot shelters on to map

#Extract pike data
pike_dat <- 
  muddyfoot_dat %>% 
  filter(Species == 'Pike' & Date == '2022-09-25')

#combine it with individual roach data
encounter_test_data = rbind(pike_dat, F59707_track)

#create move object
pred_encounter_track_mv <- 
  mt_as_move2(encounter_test_data, 
              coords = c("longitude","latitude"),
              crs = "WGS84",
              time_column = "timestamp",
              track_id_column = "individual_id",
              na.fail = F) # allows or not empty coordinates


test = encounter_test_data %>% 
  filter(individual_id == 'F59707' | individual_id == 'F59880')

p = ggplot() + 
  geom_sf(data = muddyfoot_polygon, fill = NA, color = "blue") +
  geom_point(data = test, 
             aes(x = longitude, 
                 y = latitude, 
                 group = individual_id,
                 color = individual_id), 
             size = 4, alpha = 0.3) +  # Plot the fish tracks with color by fish_id
  # geom_path(data = encounter_test_data, 
  #           aes(x = longitude, 
  #               y = latitude, 
  #               group = individual_id, 
  #               color = individual_id), 
  #           size = 1, 
  #           alpha = 0.5) +  # Add paths
  ggsn::scalebar(
    data = muddyfoot_polygon, 
    dist = 5,  # Distance represented by the scale bar
    dist_unit = "m",  # Units of measurement
    transform = TRUE,  # Transform coordinates to projected CRS if needed
    model = 'WGS84',  # Model used for transformation, if needed
    location = "bottomright",  # Position of the scale bar
    st.dist = 0.05  # Distance of the text from the scale bar
  ) +
  theme_minimal() +  
  labs(title = "Fish Movement Over Time: {frame_time}") +## Title that updates with time
  transition_reveal(along = timestamp) +  
  ease_aes('linear')

animate(p, nframes = 500, fps = 10, width = 800, height = 600)


#### Calculate distances and encounter probabilities 

#need CTMMs
pike_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
roach_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds")) 


#need telem object
pike_muddyfoot_tel = readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
roach_muddyfoot_tel = readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))
perch_muddyfoot_tel = readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))

#projections need to be the same to calculate pairwise speration distances
#current species specific telemetry objects have different projections 
#combine telemetry objects and reproject

test_telem <- list(pike_muddyfoot_tel$F59880, roach_muddyfoot_tel$F59707)
names(test_telem) <- c('F59880', 'F59707')

test_ctmm <- list(pike_muddyfoot_ctmm_fits$F59880, roach_muddyfoot_ctmm_fits$F59707)
names(test_ctmm) <- c('F59880', 'F59707')


#Center the projection on the geometric median of the data
ctmm::projection(test_telem) <- ctmm::median(test_telem)
ctmm::projection(test_telem)

new_proj <-  "+proj=aeqd +lat_0=63.7710945900604 +lon_0=20.0482803129897 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

ctmm::projection(test_ctmm$F59880) <- "+proj=aeqd +lat_0=63.7710945900604 +lon_0=20.0482803129897 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
ctmm::projection(test_ctmm$F59707)

#Pairwise separation distances
DISTS <- distances(c(test_telem$F59880, test_telem$F59707),
                   c(pike_muddyfoot_ctmm_fits$F59880, roach_muddyfoot_ctmm_fits$F59707))




pike_muddyfoot_ctmm_fits$F59880$

buffalo[c("Cilla","Mvubu")],
                   FITS[c("Cilla","Mvubu")])




#Map shelter and determine their use 

