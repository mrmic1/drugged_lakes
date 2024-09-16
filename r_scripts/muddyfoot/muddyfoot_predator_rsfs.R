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
polygon_path = "./data/lake_coords/muddyfoot/"
ctmm_path = "./data/ctmm_fits/"
data_filter_path = "./data/tracks_filtered/"
telem_path = "./data/telem_obj/"
rec_data_path = "./data/lake_coords/reciever_and_habitat_locations/"
enc_path <- "./data/encounters/"                  # Directory for encounter data
akde_path <- "./data/akdes/"                      # Directory for AKDE (Autocorrelated Kernel Density Estimation) outputs
rsf_path <- "./data/rsfs/predators/"

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

predator_ud_raster <- rast(predator_ud_raster, lake_raster, background = 0, res = 0.1)
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
writeRaster(predator_ud_raster, paste0(rsf_path, "muddyfoot_predator_ud_raster.tif"), overwrite = TRUE)
)

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

# Separate telemetry objects
pike_control_tel <- pike_muddyfoot_tel[1:3]   # Control group
pike_mix_tel <- pike_muddyfoot_tel[4:6]       # Mixed group

# Separate AKDEs for the two groups
pike_control_akdes <- pike_akdes_cg_list[1:3]
pike_mix_akdes <- pike_akdes_cg_list[4:6]


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
          file = paste0(rsf_path, "predators/muddyfoot_perch/", names(perch_control_tel)[i], "_predator_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_perch_control_list) <- names(perch_control_tel)

#check
summary(rsf_perch_control_list$F59702)
summary(rsf_perch_control_list$F59697)

#saveRDS(rsf_perch_control_list, paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_list.rds"))



### EXPOSED ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_perch_mix_list <- list()

# Assuming perch_control_tel and perch_control_akdes are lists of corresponding elements
rsf_perch_mix_list <- foreach(i = seq_along(perch_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    perch_mix_tel[[i]], 
    perch_mix_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "predators/muddyfoot_perch/", names(perch_mix_tel)[i], "_predator_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_perch_mix_list) <- names(perch_mix_tel)

saveRDS(rsf_perch_mix_list, paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_list.rds"))


#>>> 4.2 Roach ########################

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
          file = paste0(rsf_path, "predators/muddyfoot_roach/", names(roach_control_tel)[i], "_predator_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_roach_control_list) <- names(roach_control_tel)

#check
summary(rsf_roach_control_list$F59702)
summary(rsf_roach_control_list$F59697)

saveRDS(rsf_roach_control_list, paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_list.rds"))



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
          file = paste0(rsf_path, "predators/muddyfoot_roach/", names(roach_mix_tel)[i], "_predator_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_roach_mix_list) <- names(roach_mix_tel)

saveRDS(rsf_roach_mix_list, paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_list.rds"))


#>>> 4.3. Pike ########################

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
          file = paste0(rsf_path, "muddyfoot_pike/", names(pike_control_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_pike_control_list) <- names(pike_control_tel)


saveRDS(rsf_pike_control_list, paste0(rsf_path, "muddyfoot_pike/rsf_pike_control_list.rds"))



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
          file = paste0(rsf_path, "muddyfoot_pike/", names(pike_mix_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_pike_mix_list) <- names(pike_mix_tel)

saveRDS(rsf_pike_mix_list, paste0(rsf_path, "muddyfoot_pike/rsf_pike_mix_list.rds"))







