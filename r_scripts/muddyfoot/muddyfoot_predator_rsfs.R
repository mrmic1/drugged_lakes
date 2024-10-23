# --------------------------------------------------------- #
# RESOURCE SELECTION FOR PREDATOR UDS - MUDDYFOOT ####
# --------------------------------------------------------- #

library(ctmm)
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(data.table)
library(parallel)    # Parallel computing support
library(foreach)     # Looping construct for parallel execution
library(doParallel)  # Parallel backend for foreach loops

### DIRECTORIES ###
polygon_path = "./data/lake_coords/"
ctmm_path = "./data/ctmm_fits/"
data_filter_path = "./data/tracks_filtered/"
telem_path = "./data/telem_obj/"
rec_data_path = "./data/lake_coords/reciever_and_habitat_locations/"
enc_path <- "./data/encounters/"                  # Directory for encounter data
akde_path <- "./data/akdes/"                      # Directory for AKDE (Autocorrelated Kernel Density Estimation) outputs
rsf_path <- "./data/rsfs/predators/"
save_ud_plots <- "./lakes_images_traces/ud_plots/"  # Directory for saving utilization distribution plots


### LOAD DATA ###
#receiver and habitat locations
mud_rec_locs_kml <- paste0(rec_data_path, "muddyfoot_rec_hab_locations.kml")
mud_rec_locs <- st_read(mud_rec_locs_kml)[1:5,]
mud_hab_locs <- st_read(mud_rec_locs_kml)[6:7,]#receiver and habitat locations

#fish telemetry data
pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

#fish ctmms 
pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds"))

#fish akdes
pike_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/pike_akdes_cg_list.rds"))
perch_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))
roach_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))

#pike population akde
pike_total_PKDE <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/population_akde/pike_total_PKDE.rds"))

#-----------------------------#
#> 1. Create a lake raster ####
#----------------------------#

#Load lake polygon
muddyfoot_polygon <- st_read(paste0(polygon_path, "lake_muddyfoot_polygon.gpkg"))

#Transform the CRS to UTM Zone 33N (EPSG:32633)
polyProj <- st_transform(muddyfoot_polygon, crs = "EPSG:32633")

# Create a raster object with the desired resolution and CRS (UTM Zone 33N)
resolution <- 0.5
lake_raster <- rast(ext(polyProj), 
                    res = resolution, 
                    crs = "EPSG:32633")

# Reproject the raster to WGS 84 (EPSG:4326)
lake_raster <- terra::project(lake_raster, "EPSG:4326")

# Ensure that the muddyfoot_polygon is also in WGS 84 (EPSG:4326)
muddyfoot_polygon <- st_transform(muddyfoot_polygon, crs = "EPSG:4326")

# Rasterize the lake polygon (masking cells outside the lake)
lake_raster <- rasterize(muddyfoot_polygon, lake_raster, background = NA)
plot(lake_raster, main = "Lake Raster")

#---------------------------------------#
#> 2. Create a raster of predator UD ####
#---------------------------------------#

predator_ud_raster <- raster(pike_total_PKDE)
plot(predator_ud_raster) #not bounded

#raster values
pred_raster_values <- values(predator_ud_raster)
value_matrix <- as.matrix(pred_raster_values) 

predator_ud_raster <- rast(predator_ud_raster, lake_raster, background = 0, res = 0.5)
plot(predator_ud_raster)
predator_ud_raster <- terra::project(predator_ud_raster, "EPSG:4326")


#To ensure that only the lake area is included in the habitat raster (and exclude areas outside the lake), 
#mask the raster using the lake polygon:
# Mask the habitat raster to the lake polygon
predator_ud_raster <- mask(predator_ud_raster, muddyfoot_polygon)
plot(predator_ud_raster, main = "Predator raster")

#convert to raster
#Need to convert it from terra to raster to work with rsf.fit()
predator_ud_raster <- raster::raster(predator_ud_raster)
plot(predator_ud_raster)

#save masked raster
#writeRaster(predator_ud_raster, paste0(rsf_path, "muddyfoot_predator_ud_raster.tif"), overwrite = TRUE)
)

ggsave(file = paste0(save_ud_plots, "predator_ud_raster.png"), 
       plot = predator_ud_raster, 
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


#---------------------------------------------#
#> 3. Split telemetry objects by treatment ####
#---------------------------------------------#

#first need to remove individuals that were predated and tracked while in stomack of pikes
# Load predation event data (pre-identified predation events)
mud_pred_events <- readRDS(paste0(enc_path, "muddyfoot_pred_events.rds"))

# Remove Roach individuals that were predated based on predation events
roach_ids_remove <- mud_pred_events %>%
  filter(Species == "Roach") %>%    # Filter for Roach species
  pull(individual_ID)               # Extract IDs of predated Roach

# Update the Roach AKDE list and telemetry data, removing predated individuals
roach_akdes_cg_list <- roach_akdes_cg_list[!(names(roach_akdes_cg_list) %in% roach_ids_remove)]
roach_muddyfoot_tel <- roach_muddyfoot_tel[!(names(roach_muddyfoot_tel) %in% roach_ids_remove)]

# Remove Perch individuals that were predated based on predation events
perch_ids_remove <- mud_pred_events %>%
  filter(Species == "Perch") %>%    # Filter for Perch species
  pull(individual_ID)               # Extract IDs of predated Perch

# Update the Perch AKDE list and telemetry data, removing predated individuals
perch_akdes_cg_list <- perch_akdes_cg_list[!(names(perch_akdes_cg_list) %in% perch_ids_remove)]
perch_muddyfoot_tel <- perch_muddyfoot_tel[!(names(perch_muddyfoot_tel) %in% perch_ids_remove)]


#use perch as an example
#telemetry objects
perch_control_tel <- perch_muddyfoot_tel[1:14]
perch_mix_tel <- perch_muddyfoot_tel[15:29]

#akdes
perch_control_akdes <- perch_akdes_cg_list[1:14]
perch_mix_akdes <- perch_akdes_cg_list[15:29]

#telemetry objects
roach_control_tel <- roach_muddyfoot_tel[1:11]
roach_mix_tel <- roach_muddyfoot_tel[12:24]

#akdes
roach_control_akdes <- roach_akdes_cg_list[1:11]
roach_mix_akdes <- roach_akdes_cg_list[12:24]



#-----------------------#
#> 4. Fit RSF models ####
#-----------------------#

#>>> 4.1 Perch ########################

### CONTROL ###

cl <- makeCluster(10)
doParallel::registerDoParallel(cl)
rsf_perch_control_list <- list()

# Assuming perch_control_tel and perch_control_akdes are lists of corresponding elements
rsf_perch_control_list <- foreach(i = seq_along(perch_control_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_control_model <- rsf.fit(
    perch_control_tel[[i]], 
    perch_control_akdes[[i]], 
    R=list(habitat1=predator_ud_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_control_model, 
          file = paste0(rsf_path, "muddyfoot_perch/", names(perch_control_tel)[i], "_predator_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_perch_control_list) <- names(perch_control_tel)

saveRDS(rsf_perch_control_list, paste0(rsf_path, "muddyfoot_perch/pred_rsf_perch_control_list.rds"))



### EXPOSED ###

cl <- makeCluster(10)
doParallel::registerDoParallel(cl)
rsf_perch_mix_list <- list()

# Assuming perch_control_tel and perch_control_akdes are lists of corresponding elements
rsf_perch_mix_list <- foreach(i = seq_along(perch_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    perch_mix_tel[[i]], 
    perch_mix_akdes[[i]], 
    R=list(habitat1=predator_ud_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "muddyfoot_perch/", names(perch_mix_tel)[i], "_predator_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_perch_mix_list) <- names(perch_mix_tel)

saveRDS(rsf_perch_mix_list, paste0(rsf_path, "muddyfoot_perch/pred_rsf_perch_mix_list.rds"))


#>>> 4.2 Roach ########################

### CONTROL ###

cl <- makeCluster(10)
doParallel::registerDoParallel(cl)
rsf_roach_control_list <- list()

# Assuming roach_control_tel and roach_control_akdes are lists of corresponding elements
rsf_roach_control_list <- foreach(i = seq_along(roach_control_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_control_model <- rsf.fit(
    roach_control_tel[[i]], 
    roach_control_akdes[[i]], 
    R=list(habitat1=predator_ud_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_control_model, 
          file = paste0(rsf_path, "muddyfoot_roach/", names(roach_control_tel)[i], "_predator_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_roach_control_list) <- names(roach_control_tel)

saveRDS(rsf_roach_control_list, paste0(rsf_path, "muddyfoot_roach/pred_rsf_roach_control_list.rds"))



### EXPOSED ###

cl <- makeCluster(10)
doParallel::registerDoParallel(cl)
rsf_roach_mix_list <- list()

# Assuming roach_control_tel and roach_control_akdes are lists of corresponding elements
rsf_roach_mix_list <- foreach(i = seq_along(roach_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    roach_mix_tel[[i]], 
    roach_mix_akdes[[i]], 
    R=list(habitat1=predator_ud_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "muddyfoot_roach/", names(roach_mix_tel)[i], "_predator_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_roach_mix_list) <- names(roach_mix_tel)

saveRDS(rsf_roach_mix_list, paste0(rsf_path, "muddyfoot_roach/pred_rsf_roach_mix_list.rds"))

#----------------------------------------------------------------------------------------------------------------#

#-----------------------------------------#
# > 5. Explore RSF model results ####
#-----------------------------------------#

#>>> 5.1. Perch ####

# Load the RSF (Resource Selection Function) results for control and exposed perch
rsf_perch_control_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/pred_rsf_perch_control_list.rds"))
rsf_perch_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/pred_rsf_perch_mix_list.rds"))

summary(rsf_perch_exposed_list) #weird value for F59736 (remove for the time being)
rsf_perch_exposed_list <- rsf_perch_exposed_list[-14]

# Calculate the mean of the RSF lists (assuming that each list is a collection of model objects)
rsf_perch_control_mean <- mean(rsf_perch_control_list)
rsf_perch_exposed_mean <- mean(rsf_perch_exposed_list)

# Extract coefficients' confidence intervals from the RSF summaries
# Convert the first confidence interval (CI) row to a data frame for each treatment
rsf_coef_control_perch <- as.data.frame(t(summary(rsf_perch_control_mean)$CI[1,]))
rsf_coef_exposed_perch <- as.data.frame(t(summary(rsf_perch_exposed_mean)$CI[1,]))

# Add a treatment column to each dataset to distinguish between Control and Exposed in the final plot
rsf_coef_control_perch$treatment <- "Control"
rsf_coef_exposed_perch$treatment <- "Exposed"

# Combine the coefficients for both treatments into one data frame
perch_predator_rsf_coefs <- rbind(rsf_coef_control_perch, rsf_coef_exposed_perch)

# Create the ggplot visualization of RSF coefficient estimates for perch habitats
(perch_predator_rsf_plot <- 
    ggplot(perch_predator_rsf_coefs, 
           aes(x = treatment, y = est)) +  # Aesthetic mapping: treatment on x-axis, estimate (est) on y-axis
    geom_errorbar(aes(ymin = low, ymax = high), 
                  width = 0.1,  # Width of error bars
                  size = 1,       # Thickness of error bars
                  color = "black") +  # Error bars in black
    geom_point(aes(shape = treatment, fill = treatment), 
               size = 4,  # Point size
               color = "black") +  # Outline of points in black
    scale_shape_manual(values = c(21, 21)) +  # Use the same shape (circle) for both treatments
    scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +  # White fill for Control, black fill for Exposed
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add a dashed horizontal line at y=0
    coord_cartesian(ylim = c(-2, 2)) +  # Set y-axis limits from -2 to 2
    labs(y = "Selection coefficient") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)

ggsave(file = paste0(save_ud_plots, "perch_predator_rsf_muddyfoot.png"), 
       plot = perch_predator_rsf_plot, 
       device = 'png',
       width = 8, 
       height = 8,
       units = 'cm',
       dpi = 300)



#>>> 5.2. Roach ####

# Load the RSF (Resource Selection Function) results for control and exposed roach
rsf_roach_control_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/pred_rsf_roach_control_list.rds"))
rsf_roach_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/pred_rsf_roach_mix_list.rds"))

# Calculate the mean of the RSF lists (assuming that each list is a collection of model objects)
rsf_roach_control_mean <- mean(rsf_roach_control_list)
rsf_roach_exposed_mean <- mean(rsf_roach_exposed_list)

# Extract coefficients' confidence intervals from the RSF summaries
# Convert the first confidence interval (CI) row to a data frame for each treatment
rsf_coef_control_roach <- as.data.frame(t(summary(rsf_roach_control_mean)$CI[1,]))
rsf_coef_exposed_roach <- as.data.frame(t(summary(rsf_roach_exposed_mean)$CI[1,]))

# Add a treatment column to each dataset to distinguish between Control and Exposed in the final plot
rsf_coef_control_roach$treatment <- "Control"
rsf_coef_exposed_roach$treatment <- "Exposed"

# Combine the coefficients for both treatments into one data frame
roach_predator_rsf_coefs <- rbind(rsf_coef_control_roach, rsf_coef_exposed_roach)

# Create the ggplot visualization of RSF coefficient estimates for roach habitats
(roach_predator_rsf_plot <- 
    ggplot(roach_predator_rsf_coefs, 
           aes(x = treatment, y = est)) +  # Aesthetic mapping: treatment on x-axis, estimate (est) on y-axis
    geom_errorbar(aes(ymin = low, ymax = high), 
                  width = 0.1,  # Width of error bars
                  size = 1,       # Thickness of error bars
                  color = "black") +  # Error bars in black
    geom_point(aes(shape = treatment, fill = treatment), 
               size = 4,  # Point size
               color = "black") +  # Outline of points in black
    scale_shape_manual(values = c(21, 21)) +  # Use the same shape (circle) for both treatments
    scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +  # White fill for Control, black fill for Exposed
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add a dashed horizontal line at y=0
    coord_cartesian(ylim = c(-2, 2)) +  # Set y-axis limits from -2 to 2
    labs(y = "Selection coefficient") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)

ggsave(file = paste0(save_ud_plots, "roach_predator_rsf_muddyfoot.png"), 
       plot = roach_predator_rsf_plot, 
       device = 'png',
       width = 8, 
       height = 8,
       units = 'cm',
       dpi = 300)
