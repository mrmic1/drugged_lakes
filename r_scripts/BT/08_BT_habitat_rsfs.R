# ---------------------------------------------------------#
# RESOURCE SELECTION FOR ARTIFICIAL HABITATS - LAKE BT ###
# ---------------------------------------------------------#

### LIBRARIES ###
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
filtered_data_path <- "./data/tracks_filtered/lake_BT/"
telem_path <- "./data/telem_obj/BT/" 
rec_data_path = "./data/lake_coords/reciever_and_habitat_locations/"
enc_path <- "./data/encounters/BT/"                  # Directory for encounter data
akde_path <- "./data/akdes/"                      # Directory for AKDE (Autocorrelated Kernel Density Estimation) outputs
rsf_path <- "./data/rsfs/habitats/"
save_plots <- "./plots/BT/"  # Directory for saving utilization distribution plots

### LOAD DATA ###
#receiver and habitat locations
BT_rec_locs_kml <- paste0(rec_data_path, "BT_rec_hab_locations.kml")
BT_locs <- st_read(BT_rec_locs_kml)
BT_rec_locs <- st_read(BT_rec_locs_kml)[1:5,]
BT_hab_locs <- st_read(BT_rec_locs_kml)[6:9,]#receiver and habitat locations

#fish telemetry data
pike_BT_tel <- readRDS(paste0(telem_path, 'pike_BT_tel.rds'))
perch_BT_tel <- readRDS(paste0(telem_path, 'perch_BT_tel.rds'))
roach_BT_tel <- readRDS(paste0(telem_path, 'roach_BT_tel.rds'))

#fish ctmms 
pike_BT_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_BT_pike_fits/lake_BT_pike_OUF_models.rds"))
perch_BT_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_BT_perch_fits/lake_BT_perch_OUF_models.rds"))
roach_BT_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_BT_roach_fits/lake_BT_roach_OUF_models.rds"))

#fish akdes
pike_akdes_cg_list <- readRDS(paste0(akde_path, "lake_BT_pike_akdes/akde_cg/lake_BT_pike_akdes_cg_list.rds"))
perch_akdes_cg_list <- readRDS(paste0(akde_path, "lake_BT_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))
roach_akdes_cg_list <- readRDS(paste0(akde_path, "lake_BT_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))

#load pkdes if necessary
perch_control_PKDE <- readRDS(paste0(akde_path, "lake_BT_perch_akdes/population_akde/perch_control_PKDE.rds"))
perch_mix_PKDE <- readRDS(paste0(akde_path, "lake_BT_perch_akdes/population_akde/perch_mix_PKDE.rds"))
roach_control_PKDE <- readRDS(paste0(akde_path, "lake_BT_roach_akdes/population_akde/roach_control_PKDE.rds"))
roach_mix_PKDE <- readRDS(paste0(akde_path, "lake_BT_roach_akdes/population_akde/roach_mix_PKDE.rds"))
pike_control_PKDE <- readRDS(paste0(akde_path, "lake_BT_pike_akdes/population_akde/lake_BT_pike_control_PKDE.rds"))
pike_mix_PKDE <- readRDS(paste0(akde_path, "lake_BT_pike_akdes/population_akde/lake_BT_pike_mix_PKDE.rds"))
pike_total_PKDE <- readRDS(paste0(akde_path, "lake_BT_pike_akdes/population_akde/lake_BT_pike_total_PKDE.rds"))


#Load lake polygon
BT_polygon <- st_read(paste0(polygon_path, "lake_BT_polygon.gpkg"))

#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#---------------------------------------------#
# 1. Plot population akdes with habitats ####
#---------------------------------------------#

# Function to generate the plot with habitats
generate_ud_plot_whabitats <- function(pkde_data, bbox, hab_locs) {
  ud_raster <- raster(pkde_data)
  masked_ud_raster <- mask(ud_raster, bbox)
  ud_df <- as.data.frame(as(masked_ud_raster, "SpatialPixelsDataFrame"))
  colnames(ud_df) <- c("value", "x", "y")
  
  ggplot() +
    geom_sf(data = bbox, color = "black") +
    geom_tile(data = ud_df, aes(x = x, y = y, fill = value), alpha = 0.6) +
    geom_sf(data = hab_locs, color = "green", size = 2, fill = NA, shape = 3, stroke = 1) + 
    scale_fill_viridis_c(na.value = 'transparent', option = 'magma') +
    coord_sf() +
    theme_classic() +
    labs(fill = "Utilization Distribution", x = "", y = '')+
    theme(legend.position = "none",
          axis.text = element_text(size = 8, color = 'black'))
}

# Function to generate the plot without habitats
generate_ud_plot <- function(pkde_data, bbox, hab_locs) {
  ud_raster <- raster(pkde_data)
  masked_ud_raster <- mask(ud_raster, bbox)
  ud_df <- as.data.frame(as(masked_ud_raster, "SpatialPixelsDataFrame"))
  colnames(ud_df) <- c("value", "x", "y")
  
  ggplot() +
    geom_sf(data = bbox, color = "black") +
    geom_tile(data = ud_df, aes(x = x, y = y, fill = value), alpha = 0.6) +
    scale_fill_viridis_c(na.value = 'transparent', option = 'magma') +
    coord_sf() +
    theme_classic() +
    labs(fill = "Utilization Distribution", x = "", y = '')+
    theme(legend.position = "none",
          axis.text = element_text(size = 8, color = 'black'))
}

# Transform bbox
BT_bbox_perch <- st_transform(BT_polygon, crs(raster(perch_control_PKDE)))
BT_bbox_roach <- st_transform(BT_polygon, crs(raster(roach_control_PKDE)))
BT_bbox_pike <- st_transform(BT_polygon, crs(raster(pike_total_PKDE)))

# Generate the plots

#> 1.1. Perch ####

## CONTROL ##
perch_control_plot <- generate_ud_plot(perch_control_PKDE, BT_bbox_perch, BT_hab_locs)
perch_control_plot_habitats <- generate_ud_plot_whabitats(perch_control_PKDE, BT_bbox_perch, BT_hab_locs)
ggsave(file = paste0(save_plots, "perch_control_UD_habitats_BT.png"), 
       plot = perch_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## EXPOSED ##
perch_exposed_plot <- generate_ud_plot(perch_mix_PKDE, BT_bbox_perch, BT_hab_locs)
perch_exposed_plot_habitats <- generate_ud_plot_whabitats(perch_mix_PKDE, BT_bbox_perch, BT_hab_locs)
ggsave(file = paste0(save_plots, "perch_exposed_UD_habitats_BT.png"), 
       plot = perch_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)



#> 1.2. Roach ####
## CONTROL ##
roach_control_plot <- generate_ud_plot(roach_control_PKDE, BT_bbox_roach, BT_hab_locs)
roach_control_plot_habitats <- generate_ud_plot_whabitats(roach_control_PKDE, BT_bbox_roach, BT_hab_locs)
ggsave(file = paste0(save_plots, "roach_control_UD_habitats_BT.png"), 
       plot = roach_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## EXPOSED ##
roach_exposed_plot <- generate_ud_plot(roach_mix_PKDE, BT_bbox_roach, BT_hab_locs)
roach_exposed_plot_habitats <- generate_ud_plot_whabitats(roach_mix_PKDE, BT_bbox_roach, BT_hab_locs)
ggsave(file = paste0(save_plots, "roach_exposed_UD_habitats_BT.png"), 
       plot = roach_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)



#> 1.3. Pike ####

## CONTROL ##
pike_control_plot <- generate_ud_plot(pike_control_PKDE, BT_bbox_pike, BT_hab_locs)
pike_control_plot_habitats <- generate_ud_plot_whabitats(pike_control_PKDE, BT_bbox_pike, BT_hab_locs)
ggsave(file = paste0(save_plots, "pike_control_UD_habitats_BT.png"), 
       plot = pike_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)



library(ctmm)

# Choose the individual to plot
pike_ud <- pike_akdes_cg_list[["F59892"]]  # Replace "F12345" with the actual name
# Or use index if unnamed: pike_ud <- pike_akdes_cg_list[[1]]

# Plot the UD
plot(pike_ud, col.level = "darkgreen")




## EXPOSED ##
pike_exposed_plot <- generate_ud_plot(pike_mix_PKDE, BT_bbox_pike, BT_hab_locs)
pike_exposed_plot_habitats <- generate_ud_plot_whabitats(pike_mix_PKDE, BT_bbox_pike, BT_hab_locs)
ggsave(file = paste0(save_plots, "pike_exposed_UD_habitats_BT.png"), 
       plot = pike_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## TOTAL ##
pike_total_plot <- generate_ud_plot(pike_total_PKDE, BT_bbox_pike, BT_hab_locs)
pike_total_plot_habitats <- generate_ud_plot_whabitats(pike_total_PKDE, BT_bbox_pike, BT_hab_locs)
ggsave(file = paste0(save_plots, "pike_total_UD_habitats_BT.png"), 
       plot = pike_total_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)

#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#-----------------------------------------#
# 2. Create a lake and habitat raster ####
#-----------------------------------------#

#Final raster is used in RSF models

#Transform the CRS to UTM Zone 33N (EPSG:32633)
polyProj <- st_transform(BT_polygon, crs = "EPSG:32633")

# Create a raster object with the desired resolution and CRS (UTM Zone 33N)
resolution <- 0.5
lake_raster <- rast(ext(polyProj), 
                    res = resolution, 
                    crs = "EPSG:32633")

# Reproject the raster to WGS 84 (EPSG:4326)
lake_raster <- terra::project(lake_raster, "EPSG:4326")

# Ensure that the BT_polygon is also in WGS 84 (EPSG:4326)
BT_polygon <- st_transform(BT_polygon, crs = "EPSG:4326")

# Rasterize the lake polygon (masking cells outside the lake)
lake_raster <- rasterize(BT_polygon, lake_raster, background = NA)
plot(lake_raster, main = "Lake Raster")

#Plot habitats onto lake raster

habitat_patches <- vect(BT_hab_locs)

# Assign a value of 1 to the habitat patches and leave other areas as 0 or NA
habitat_raster <- rasterize(habitat_patches, lake_raster, field = 1, background = 0)
plot(habitat_raster, main = "Habitat Raster")

#To ensure that only the lake area is included in the habitat raster (and exclude areas outside the lake), 
#mask the raster using the lake polygon:
# Mask the habitat raster to the lake polygon
habitat_raster <- mask(habitat_raster, BT_polygon)
plot(habitat_raster, main = "Habitat Raster Masked by Lake")

#convert to raster
#Need to convert it from terra to raster to work with rsf.fit()
habitat_raster <- raster::raster(habitat_raster)
plot(habitat_raster, main = "Habitat Raster Masked by Lake")

#save raster
writeRaster(habitat_raster, "./data/lake_coords/BT_habitat_raster.grd", overwrite = TRUE)
#open raster
habitat_raster <- rast("./data/lake_coords/BT_habitat_raster.grd")


#---------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------#
#3. Split telemetry objects by treatment ####
#---------------------------------------------#

# Separating into 'control' and 'mix'

#> 3.1 Perch ########################

# Check that treatments are in the right order
unique_ids <- names(perch_BT_tel)

# Create a data frame to store ID and unique Treatment information
result_table <- do.call(rbind, lapply(unique_ids, function(id) {
  # Extract the unique Treatment value for the current ID
  treatment_info <- unique(perch_BT_tel[[id]]$Treatment)
  
  # Ensure only one unique Treatment value is captured
  if (length(treatment_info) > 1) {
    warning(paste("Multiple treatments found for ID:", id, 
                  "using the first value. Treatments:", paste(treatment_info, collapse = ", ")))
    treatment_info <- treatment_info[1]
  }
  
  # Create a data frame for the ID and its treatment
  data.frame(ID = id, Treatment = treatment_info)
}))

# View the resulting table
print(result_table)

#use perch as an example
#telemetry objects
perch_control_tel <- perch_BT_tel[1:15]
perch_mix_tel <- perch_BT_tel[16:30]

#akdes
perch_control_akdes <- perch_akdes_cg_list[1:15]
perch_mix_akdes <- perch_akdes_cg_list[16:30]


#> 3.2 Roach ########################

# Check that treatments are in the right order
unique_ids <- names(roach_BT_tel)

# Create a data frame to store ID and unique Treatment information
result_table <- do.call(rbind, lapply(unique_ids, function(id) {
  # Extract the unique Treatment value for the current ID
  treatment_info <- unique(roach_BT_tel[[id]]$Treatment)
  
  # Ensure only one unique Treatment value is captured
  if (length(treatment_info) > 1) {
    warning(paste("Multiple treatments found for ID:", id, 
                  "using the first value. Treatments:", paste(treatment_info, collapse = ", ")))
    treatment_info <- treatment_info[1]
  }
  
  # Create a data frame for the ID and its treatment
  data.frame(ID = id, Treatment = treatment_info)
}))

print(result_table)



#telemetry objects
roach_control_tel <- roach_BT_tel[1:15]
roach_mix_tel <- roach_BT_tel[16:30]

#akdes
roach_control_akdes <- roach_akdes_cg_list[1:15]
roach_mix_akdes <- roach_akdes_cg_list[16:30]


#> 3.3 Pike ########################

# Separate telemetry objects
pike_control_tel <- pike_BT_tel[1:3]   # Control group
pike_mix_tel <- pike_BT_tel[4:6]       # Mixed group

# Separate AKDEs for the two groups
pike_control_akdes <- pike_akdes_cg_list[1:3]
pike_mix_akdes <- pike_akdes_cg_list[4:6]

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------#
# 4. Fit RSF models ####
#-----------------------#

#> 4.1 Perch ########################

### CONTROL ###

cl <- makeCluster(5)
doParallel::registerDoParallel(cl)
rsf_perch_control_list <- list()

# Assuming perch_control_tel and perch_control_akdes are lists of corresponding elements
rsf_perch_control_list <- foreach(i = seq_along(perch_control_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_control_model <- rsf.fit(
    perch_control_tel[[i]], 
    perch_control_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_control_model, 
          file = paste0(rsf_path, "BT_perch/", names(perch_control_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_perch_control_list) <- names(perch_control_tel)

#check
summary(rsf_perch_control_list$F59702)
summary(rsf_perch_control_list$F59697)

saveRDS(rsf_perch_control_list, paste0(rsf_path, "BT_perch/rsf_perch_control_list.rds"))


### EXPOSED ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_perch_mix_list <- list()

rsf_perch_mix_list <- foreach(i = seq_along(perch_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    perch_mix_tel[[i]], 
    perch_mix_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "BT_perch/", names(perch_mix_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_perch_mix_list) <- names(perch_mix_tel)

#check
summary(rsf_perch_mix_list$F59714)
summary(rsf_perch_mix_list$F59727)

saveRDS(rsf_perch_mix_list, paste0(rsf_path, "BT_perch/rsf_perch_mix_list.rds"))

#-----------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------#

#> 4.2 Roach ########################

### CONTROL ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_roach_control_list <- list()

# Assuming roach_control_tel and roach_control_akdes are lists of corresponding elements
rsf_roach_control_list <- foreach(i = seq_along(roach_control_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_control_model <- rsf.fit(
    roach_control_tel[[i]], 
    roach_control_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_control_model, 
          file = paste0(rsf_path, "BT_roach/", names(roach_control_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_roach_control_list) <- names(roach_control_tel)

#check
summary(rsf_roach_control_list$F59683)
summary(rsf_roach_control_list$F59684)

saveRDS(rsf_roach_control_list, paste0(rsf_path, "BT_roach/rsf_roach_control_list.rds"))



### EXPOSED ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_roach_mix_list <- list()

# Assuming roach_control_tel and roach_control_akdes are lists of corresponding elements
rsf_roach_mix_list <- foreach(i = seq_along(roach_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    roach_mix_tel[[i]], 
    roach_mix_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "BT_roach/", names(roach_mix_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_roach_mix_list) <- names(roach_mix_tel)

saveRDS(rsf_roach_mix_list, paste0(rsf_path, "BT_roach/rsf_roach_mix_list.rds"))

#-----------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------#

#> 4.3. Pike ########################

### CONTROL ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_pike_control_list <- list()

# Assuming pike_control_tel and pike_control_akdes are lists of corresponding elements
rsf_pike_control_list <- foreach(i = seq_along(pike_control_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_control_model <- rsf.fit(
    pike_control_tel[[i]], 
    pike_control_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_control_model, 
          file = paste0(rsf_path, "BT_pike/", names(pike_control_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_pike_control_list) <- names(pike_control_tel)


saveRDS(rsf_pike_control_list, paste0(rsf_path, "BT_pike/rsf_pike_control_list.rds"))


### EXPOSED ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_pike_mix_list <- list()

# Assuming pike_control_tel and pike_control_akdes are lists of corresponding elements
rsf_pike_mix_list <- foreach(i = seq_along(pike_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    pike_mix_tel[[i]], 
    pike_mix_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "BT_pike/", names(pike_mix_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_pike_mix_list) <- names(pike_mix_tel)

saveRDS(rsf_pike_mix_list, paste0(rsf_path, "BT_pike/rsf_pike_mix_list.rds"))

#-------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------#

#-----------------------------------------#
# 5. Explore RSF model results ####
#-----------------------------------------#

#> 5.1. Perch ####

# Load the RSF (Resource Selection Function) results for control and exposed perch
rsf_perch_control_list <- readRDS(paste0(rsf_path, "BT_perch/rsf_perch_control_list.rds"))
rsf_perch_exposed_list <- readRDS(paste0(rsf_path, "BT_perch/rsf_perch_mix_list.rds"))

# Calculate the mean of the RSF lists (assuming that each list is a collection of model objects)
# This assumes 'mean()' can handle these model objects, modify if you need more specific behavior
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
perch_habitat_rsf_coefs <- rbind(rsf_coef_control_perch, rsf_coef_exposed_perch)

# Create the ggplot visualization of RSF coefficient estimates for perch habitats
(perch_habitat_rsf_plot <- 
    ggplot(perch_habitat_rsf_coefs, 
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
    coord_cartesian(ylim = c(-1, 6)) +  # Set y-axis limits from -1 to 6
    labs(y = "Selection coefficient") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)


ggsave(file = paste0(save_ud_plots, "perch_habitats_rsf_BT.png"), 
       plot = perch_habitat_rsf_plot, 
       device = 'png',
       width = 8, 
       height = 8,
       units = 'cm',
       dpi = 300)

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#


#> 5.2 Roach ####

# Load the RSF (Resource Selection Function) results for control and exposed roach
rsf_roach_control_list <- readRDS(paste0(rsf_path, "BT_roach/rsf_roach_control_list.rds"))
rsf_roach_exposed_list <- readRDS(paste0(rsf_path, "BT_roach/rsf_roach_mix_list.rds"))

# Calculate the mean of the RSF lists (assuming that each list is a collection of model objects)
# This assumes 'mean()' can handle these model objects, modify if you need more specific behavior
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
roach_habitat_rsf_coefs <- rbind(rsf_coef_control_roach, rsf_coef_exposed_roach)

# Create the ggplot visualization of RSF coefficient estimates for roach habitats
(roach_habitat_rsf_plot <- 
    ggplot(roach_habitat_rsf_coefs, 
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
    coord_cartesian(ylim = c(-1, 6)) +  # Set y-axis limits from -1 to 6
    labs(y = "Selection coefficient") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)


ggsave(file = paste0(save_ud_plots, "roach_habitats_rsf_BT.png"), 
       plot = roach_habitat_rsf_plot, 
       device = 'png',
       width = 8, 
       height = 8,
       units = 'cm',
       dpi = 300)


#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

#> 5.3 Pike ####

# Load the RSF (Resource Selection Function) results for control and exposed pike
rsf_pike_control_list <- readRDS(paste0(rsf_path, "BT_pike/rsf_pike_control_list.rds"))
rsf_pike_exposed_list <- readRDS(paste0(rsf_path, "BT_pike/rsf_pike_mix_list.rds"))

# Calculate the mean of the RSF lists (assuming that each list is a collection of model objects)
# This assumes 'mean()' can handle these model objects, modify if you need more specific behavior
rsf_pike_control_mean <- mean(rsf_pike_control_list)
rsf_pike_exposed_mean <- mean(rsf_pike_exposed_list)

# Extract coefficients' confidence intervals from the RSF summaries
# Convert the first confidence interval (CI) row to a data frame for each treatment
rsf_coef_control_pike <- as.data.frame(t(summary(rsf_pike_control_mean)$CI[1,]))
rsf_coef_exposed_pike <- as.data.frame(t(summary(rsf_pike_exposed_mean)$CI[1,]))

# Add a treatment column to each dataset to distinguish between Control and Exposed in the final plot
rsf_coef_control_pike$treatment <- "Control"
rsf_coef_exposed_pike$treatment <- "Exposed"

# Combine the coefficients for both treatments into one data frame
pike_habitat_rsf_coefs <- rbind(rsf_coef_control_pike, rsf_coef_exposed_pike)

# Create the ggplot visualization of RSF coefficient estimates for pike habitats
(pike_habitat_rsf_plot <- 
    ggplot(pike_habitat_rsf_coefs, 
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
    coord_cartesian(ylim = c(-1, 6)) +  # Set y-axis limits from -1 to 6
    labs(y = "Selection coefficient") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)


ggsave(file = paste0(save_ud_plots, "pike_habitats_rsf_BT.png"), 
       plot = pike_habitat_rsf_plot, 
       device = 'png',
       width = 8, 
       height = 8,
       units = 'cm',
       dpi = 300)

#------------------------------------------------------------------------------------------------------------------------------#

