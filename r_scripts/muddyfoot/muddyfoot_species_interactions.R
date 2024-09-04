# -----------------------------------------------# 
# Exploring interspecific interactions - Muddyfoot 
# ------------------------------------------------#

### LIBRARIES ###
library(dplyr)
library(move)
library(move2)
library(ctmm)
library(lubridate)


### DIRECTORIES ###
data_filter_path = "./data/tracks_filtered/"
lake_polygon_path = "./data/lake_coords/"
ctmm_path = "./data/ctmm_fits/"
telem_path = "./data/telem_obj/"
enc_path = "./data/encounters/"


#muddyfoot filtered dataframe
muddyfoot_filt_data <-  readRDS(paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))

#telemetry objects
pike_muddyfoot_tel = readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel = readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel = readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

#ctmm models
pike_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds")) 

#load in muddyfoot polygon
muddyfoot_polygon = sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))


#### Calculate distances and encounter probabilities 

#I want to calculate the distances between roach 59707 and all pike for shared timepoints
roach_sub_tel <- roach_muddyfoot_tel['F59707']
roach_sub_ctmm <- roach_muddyfoot_ctmm_fits['F59707']

#projections need to be the same to calculate pairwise seperation distances
#current species specific telemetry objects have different projections 
#combine telemetry objects and reproject

# Reproject roach objects to match pike projection
projection(roach_sub_tel) <- projection(pike_muddyfoot_tel)
projection(roach_sub_ctmm) <- projection(pike_muddyfoot_ctmm_fits)

# Initialize a list to store the predicted location differences
location_differences <- list()

# Loop over each pike (assuming the pike data starts from the second element)
for(i in 1:length(pike_muddyfoot_tel)) {
  # Create a list of telemetry data for the roach and the current pike
  combined_telemetry <- c(roach_sub_tel, pike_muddyfoot_tel[i])
  
  # Create a list of CTMM fits for the roach and the current pike
  combined_ctmm <- c(roach_sub_ctmm, pike_muddyfoot_ctmm_fits[i])
  
  # Predict the location difference between the roach and the current pike
  location_difference <- distances(combined_telemetry, combined_ctmm)
  
  # Extract the IDs of the roach and the current pike
  roach_id <- names(roach_sub_tel)[1]  # Use names of the list if available
  pike_id <- names(pike_muddyfoot_tel)[i]   # Use names of the list if available
  
  # Convert the result to a data frame and add the IDs as columns
  location_difference_df <- as.data.frame(location_difference)
  location_difference_df$Roach_ID <- roach_id
  location_difference_df$Pike_ID <- pike_id
  
  # Store the result in the list
  location_differences[[i]] <- location_difference_df
}

# Combine all the data frames into one
location_differences_df <- do.call(rbind, location_differences)

# Display the data frame with the predicted location differences and IDs
print(location_differences_df)

# Combine all the data frames into one
location_differences_df <- do.call(rbind, location_differences)

#timezone error/need to fix
location_differences_df$timestamp <- as.POSIXct(location_differences_df$timestamp, tz = "Europe/Stockholm")


#Visualise the separation distances
# Set up a 3x3 plotting layout
par(mfrow = c(2, 3))

# Loop over each set of distance data
for(i in 1:length(location_differences)) {
  # Extract the current distance data
  DISTS <- location_differences[[i]]
  
  # Create the plot for the current roach/pike combination
  plot(DISTS$est ~ DISTS$timestamp,
       type = "l",
       col = "#5e548e",
       main = paste("Roach:", DISTS$Roach_ID[1], "Pike:", DISTS$Pike_ID[1]),
       xlab = "Timestamp",
       ylab = "Distance (m)")
}

# Reset the plotting layout to default
par(mfrow = c(1, 1))



### Empirical encounters ###

location_differences_df$encounter <- ifelse(location_differences_df$est <= 1, 1, 0) 
#encounter less than 1 metre. 
#might want to consider the strike distance of pike as the threshold distance

#summary of possible encounters
location_differences_df %>%
  group_by(Roach_ID, Pike_ID) %>%
  summarize(Encounter_Count = sum(encounter), .groups = 'drop')
#looks only had encounters with F59880 and F59885

# Taking difference so that it only counts it as one encounter if it enters into that radius and stays there
unique_encounter <- 
  location_differences_df %>% 
  group_by(Roach_ID, Pike_ID) %>% 
  reframe(diff = diff(encounter)) %>% 
  ungroup() %>% 
  group_by(Roach_ID, Pike_ID) %>% 
  summarise(unique_encounter = length(which(diff == 1)))


#save location_differences_df
#saveRDS(location_differences_df, paste0(enc_path, "test_pred_prey.rds"))


# Get unique combinations of Roach_ID and Pike_ID
unique_pairs <- unique(location_differences_df[, c("Roach_ID", "Pike_ID")])


# Set up a 3x3 plotting layout
par(mfrow = c(2, 3))

# Loop over each unique Roach and Pike combination
for(i in 1:nrow(unique_pairs)) {
  # Extract the current Roach_ID and Pike_ID
  current_roach <- unique_pairs$Roach_ID[i]
  current_pike <- unique_pairs$Pike_ID[i]
  
  # Subset the dataframe for the current Roach and Pike combination
  DISTS <- subset(location_differences_df, Roach_ID == current_roach & Pike_ID == current_pike)
  
  # Create the plot for the current Roach/Pike combination
  plot(DISTS$encounter ~ DISTS$timestamp,
       col = "#5e548e",
       main = paste("Roach:", current_roach, "Pike:", current_pike),
       xlab = "Timestamp",
       ylab = "Distance (m)")
}
# Reset the plotting layout to default
par(mfrow = c(1, 1))








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

combined_telemetry
