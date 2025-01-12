#--------------------------------------------#
# HOME RANGE ANALYSIS - MUDDYFOOT
# -------------------------------------------#

# Load necessary libraries
library(ctmm)        # Continuous-Time Movement Models for animal telemetry data analysis
library(dplyr)       # Data manipulation and pipeline functions
library(parallel)    # Parallel computing support
library(foreach)     # Looping construct for parallel execution
library(doParallel)  # Parallel backend for foreach loops
library(ggplot2)     
library(sf)          
library(terra)       
library(raster)      

# Define paths to directories for loading/saving data
ctmm_path <- "./data/ctmm_fits/"                  # Directory for ctmm model fits
telem_path <- "./data/telem_obj/muddyfoot/"       # Directory for telemetry objects
akde_path <- "./data/akdes/"                      # Directory for AKDE (Autocorrelated Kernel Density Estimation) outputs
lake_polygon_path <- "./data/lake_coords/"        # Directory for lake polygon (boundary) data
rec_data_path <- "./data/lake_coords/reciever_and_habitat_locations/"  # Directory for receiver and habitat location data
save_ud_plots <- "./lakes_images_traces/ud_plots/"  # Directory for saving utilization distribution plots
enc_path <- "./data/encounters/"                  # Directory for encounter data
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"


### LOAD DATA ###

# Load telemetry objects for pike, perch, and roach species in Muddyfoot lake
pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

# Load ctmm model fits for pike, perch, and roach in Muddyfoot lake
# The models include continuous-time movement fits using the Ornstein-Uhlenbeck Foraging (OUF) model
pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds"))

# Load the polygon representing the boundary of Muddyfoot lake, in the form of a GeoPackage file
muddyfoot_polygon <- sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))

# Convert the simple features (sf) polygon data into a Spatial object for compatibility with other functions
muddyfoot_sp_data <- as(muddyfoot_polygon, "Spatial")

#-------------------------------------------------------------------------------------------------------#

#--------------------------#                 
# 1. AKDE ESTIMATION ####
#--------------------------#

# This script estimates the Autocorrelated Kernel Density Estimates (AKDEs)
# for three fish species: Pike, Perch, and Roach, from telemetry data. 
# We first estimate AKDEs on a consistent grid for each species, 
# using one individual's home range (HR) as a reference to standardize grid size.
#--------------------------#

#> 1.1. Pike ####

# Pike AKDE estimation with parallel processing
#Get reference akde for consistent grid
#Pick individuals with full dataset
unique(pike_muddyfoot_tel[[1]]$Date)
nrow(pike_muddyfoot_tel[[1]])
pike_akde_ref <- akde(pike_muddyfoot_tel[[1]], pike_muddyfoot_ctmm_fits[[1]], weights = FALSE)

# Initialize the list to store AKDEs with consistent grid
pike_akdes_cg <- list()

# Set up parallel processing using 3 cores
cl <- makeCluster(3)
doParallel::registerDoParallel(cl)

# Estimate AKDE for each individual pike using the reference grid
pike_akdes_cg <- foreach(i = 1:length(pike_muddyfoot_tel), .packages = 'ctmm') %dopar% {
  akde_fit_cg <- akde(
    pike_muddyfoot_tel[[i]],                  # Telemetry data for the i-th pike
    pike_muddyfoot_ctmm_fits[[i]],            # Corresponding CTMM fit
    weights = FALSE,                          # No weighting
    SP = muddyfoot_sp_data,                   # Spatial polygon for the lake boundary
    SP.in = TRUE,                             # Ensure AKDE confines the movement within the polygon
    grid = list(dr = pike_akde_ref$dr,             # Use reference grid resolution
                align.to.origin = TRUE))      # Align the grid to the origin
  # Save the AKDE result for each individual pike
  saveRDS(akde_fit_cg, file = paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/", names(pike_muddyfoot_tel)[i], "_akde_cg.rds"))
  akde_fit_cg
}

# Stop parallel cluster after computations
stopCluster(cl)

# Assign individual IDs to the AKDE list
names(pike_akdes_cg) <- names(pike_muddyfoot_tel)

# View summary of a specific individual's AKDE (e.g., F59880)
summary(pike_akdes_cg$F59880)

# Save the complete list of Pike AKDEs with consistent grid
saveRDS(pike_akdes_cg, paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/pike_akdes_cg_list.rds"))

#------------------------------------------------#

#> 1.2. Perch ####

# Perch AKDE estimation with parallel processing

#Get reference akde for consistent grid
#Pick individuals with full dataset
unique(perch_muddyfoot_tel[[9]]$Date)
nrow(perch_muddyfoot_tel[[9]])
perch_akde_ref <- akde(perch_muddyfoot_tel[[9]], perch_muddyfoot_ctmm_fits[[9]], weights = FALSE)

# Initialize the list to store AKDEs with consistent grid for Perch
perch_akdes_cg <- list()

# Set up parallel processing using 3 cores
cl <- makeCluster(3)
doParallel::registerDoParallel(cl)

# Estimate AKDE for each individual perch using the reference grid
perch_akdes_cg <- foreach(i = 1:length(perch_muddyfoot_tel), .packages = 'ctmm') %dopar% {
  akde_fit_cg <- akde(
    perch_muddyfoot_tel[[i]],                 # Telemetry data for the i-th perch
    perch_muddyfoot_ctmm_fits[[i]],           # Corresponding CTMM fit
    weights = FALSE,                          # No weighting
    SP = muddyfoot_sp_data,                   # Spatial polygon for the lake boundary
    SP.in = TRUE,                             # Ensure AKDE confines the movement within the polygon
    grid = list(dr = perch_akde_ref$dr,             # Use reference grid resolution
                align.to.origin = TRUE))      # Align the grid to the origin
  # Save the AKDE result for each individual perch
  saveRDS(akde_fit_cg, file = paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/", names(perch_muddyfoot_tel)[i], "_akde_cg.rds"))
  akde_fit_cg
}

# Stop parallel cluster after computations
stopCluster(cl)

# Assign individual IDs to the AKDE list
names(perch_akdes_cg) <- names(perch_muddyfoot_tel)

# View summary of all perch AKDEs
summary(perch_akdes_cg)

# Save the complete list of Perch AKDEs with consistent grid
saveRDS(perch_akdes_cg, paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))

#----------------------------------------------------------#

#> 1.3. Roach ####

# Roach AKDE estimation with parallel processing

#Get reference akde for consistent grid
#Pick individuals with full dataset
unique(roach_muddyfoot_tel[[9]]$Date)
nrow(roach_muddyfoot_tel[[9]])
roach_akde_ref <- akde(roach_muddyfoot_tel[[9]], roach_muddyfoot_ctmm_fits[[9]], weights = FALSE)

# Initialize the list to store AKDEs with consistent grid for Roach
roach_akdes_cg <- list()

# Set up parallel processing using 3 cores
cl <- makeCluster(10)
doParallel::registerDoParallel(cl)

# Estimate AKDE for each individual roach using the reference grid
roach_akdes_cg <- foreach(i = 1:length(roach_muddyfoot_tel), .packages = 'ctmm') %dopar% {
  akde_fit_cg <- akde(
    roach_muddyfoot_tel[[i]],                 # Telemetry data for the i-th roach
    roach_muddyfoot_ctmm_fits[[i]],           # Corresponding CTMM fit
    weights = FALSE,                          # No weighting
    SP = muddyfoot_sp_data,                   # Spatial polygon for the lake boundary
    SP.in = TRUE,                             # Ensure AKDE confines the movement within the polygon
    grid = list(dr = roach_akde_ref$dr,             # Use reference grid resolution
                align.to.origin = TRUE))      # Align the grid to the origin
  # Save the AKDE result for each individual roach
  saveRDS(akde_fit_cg, file = paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/", names(roach_muddyfoot_tel)[i], "_akde_cg.rds"))
  akde_fit_cg
}

# Stop parallel cluster after computations
stopCluster(cl)

# Assign individual IDs to the AKDE list
names(roach_akdes_cg) <- names(roach_muddyfoot_tel)

# View summary of all roach AKDEs
summary(roach_akdes_cg)

# Save the complete list of Roach AKDEs with consistent grid
saveRDS(roach_akdes_cg, paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))

#---------------------------------------------------------------------------------------------------------------#

#-------------------------------------#
# 2. POPULATION AKDE ESTIMATION ####
#-------------------------------------#

# This script estimates population-level Autocorrelated Kernel Density Estimates (AKDEs)
# for Pike, Perch, and Roach species, excluding individuals that were predated.

#> 2.1. Remove predated individuals from population AKDE estimation ####

# Load AKDEs for each species from previously saved files
pike_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/pike_akdes_cg_list.rds"))
perch_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))
roach_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))

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

#> 2.2. Pike  ####

# Separating Pike data into two groups: 'control' and 'mix'
# 'control' refers to the first three individuals, 'mix' refers to the next three

# Separate telemetry objects
pike_control_tel <- pike_muddyfoot_tel[1:3]   # Control group
pike_mix_tel <- pike_muddyfoot_tel[4:6]       # Mixed group

# Separate AKDEs for the two groups
pike_control_akdes <- pike_akdes_cg_list[1:3]
pike_mix_akdes <- pike_akdes_cg_list[4:6]

# Calculate population-level AKDE for the control group
pike_control_PKDE <- pkde(pike_control_tel,   # Telemetry data for the control group
                          pike_control_akdes, # AKDEs for the control group
                          SP = muddyfoot_sp_data,  # Spatial polygon for the lake boundary
                          SP.in = TRUE)       # Ensure the PKDE confines the movements within the polygon

# Save the population-level AKDE for the control group
saveRDS(pike_control_PKDE, paste0(akde_path, "muddyfoot_pike_akdes/population_akde/pike_control_PKDE.rds"))

# Calculate population-level AKDE for the mixed group
pike_mix_PKDE <- pkde(pike_mix_tel,           # Telemetry data for the mixed group
                      pike_mix_akdes,         # AKDEs for the mixed group
                      SP = muddyfoot_sp_data, # Spatial polygon for the lake boundary
                      SP.in = TRUE)           # Ensure the PKDE confines the movements within the polygon

saveRDS(pike_mix_PKDE, paste0(akde_path, "muddyfoot_pike_akdes/population_akde/pike_mix_PKDE.rds"))

# Calculate population-level AKDE for all pike - to be used as probability distribution raster
pike_total_PKDE <- pkde(pike_muddyfoot_tel,           # Telemetry data for individuals
                        pike_akdes_cg_list,         # AKDEs for all individuals
                        SP = muddyfoot_sp_data, # Spatial polygon for the lake boundary
                        SP.in = TRUE)           # Ensure the PKDE confines the movements within the polygon

# Save the population-level AKDE for the mixed group
saveRDS(pike_total_PKDE, paste0(akde_path, "muddyfoot_pike_akdes/population_akde/pike_total_PKDE.rds"))



#> 2.3. Perch ####

# Separating into 'control' and 'mix'
#telemetry objects
perch_control_tel <- perch_muddyfoot_tel[1:14]
perch_mix_tel <- perch_muddyfoot_tel[15:29]

#akdes
perch_control_akdes <- perch_akdes_cg_list[1:14]
perch_mix_akdes <- perch_akdes_cg_list[15:29]

#calculate population-level autocorrelated kernel density home range estimates for each treatment
perch_control_PKDE <- pkde(perch_control_tel,
                          perch_control_akdes, 
                          SP = muddyfoot_sp_data,
                          SP.in = TRUE)

#saveRDS(perch_control_PKDE, paste0(akde_path, "muddyfoot_perch_akdes/population_akde/perch_control_PKDE.rds"))

perch_mix_PKDE <- pkde(perch_mix_tel,
                      perch_mix_akdes, 
                      SP = muddyfoot_sp_data,
                      SP.in = TRUE)

saveRDS(perch_mix_PKDE, paste0(akde_path, "muddyfoot_perch_akdes/population_akde/perch_mix_PKDE.rds"))

# > 2.4. Roach ####

# Separating into 'control' and 'mix'
#telemetry objects
roach_control_tel <- roach_muddyfoot_tel[1:11]
roach_mix_tel <- roach_muddyfoot_tel[12:24]

#akdes
roach_control_akdes <- roach_akdes_cg_list[1:11]
roach_mix_akdes <- roach_akdes_cg_list[12:24]

#calculate population-level autocorrelated kernel density home range estimates
roach_control_PKDE <- pkde(roach_control_tel,
                          roach_control_akdes,
                          SP = muddyfoot_sp_data,
                          SP.in = TRUE)

saveRDS(roach_control_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/population_akde/roach_control_PKDE.rds"))

roach_mix_PKDE <- pkde(roach_mix_tel,
                      roach_mix_akdes, 
                      SP = muddyfoot_sp_data,
                      SP.in = TRUE)

saveRDS(roach_mix_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/population_akde/roach_mix_PKDE.rds"))


#-------------------------------------#
# >3. Compare home range size --------#
#-------------------------------------#

#> 3.1. Pike ####

# Combining 'control' and 'mix' back into a 'total' list
pike_akde_total <- list(Control = pike_control_akdes, Exposed = pike_mix_akdes)
pike_akde_meta_data <- ctmm::meta(pike_akde_total,col='black',sort=F, verbose = T, level.UD = 0.95)
pike_HR_control <- as.data.frame(t(pike_akde_meta_data$Control[1,]))
pike_HR_exposed <- as.data.frame(t(pike_akde_meta_data$Exposed[1,]))
  
# Add a treatment column to each dataset to distinguish between Control and Exposed in the final plot
pike_HR_control$treatment <- "Control"
pike_HR_exposed$treatment <- "Exposed"


# Combine the coefficients for both treatments into one data frame
pike_HR_coefs <- rbind(pike_HR_control, pike_HR_exposed)

# Create the ggplot visualization of RSF coefficient estimates for perch habitats
(pike_HR_coefs_plot <- 
    ggplot(pike_HR_coefs, 
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
    labs(y = "Home range size (m^2)") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)


#Plotting AKDEs by treatment

COL <- color(pike_akdes_cg_list, by='individual')
COL_control <- COL[1:3]
COL_mix <- COL[4:6]

#Calculate extent
EXT <- as.data.frame(t(as.matrix(extent(muddyfoot_sp_data))))

# plot AKDEs
plot(pike_control_akdes,
     col.UD = COL_control,
     col.DF = COL_control,
     col.level =COL_control, 
     col.grid = NA,
     SP = muddyfoot_sp_data,
     level.UD = 0.5,
     xlim = EXT$x,
     ylim = EXT$y)


plot(mix_AKDES,col.DF=COL_mix,col.level=COL_mix,col.grid=NA,level=NA,main="Mix AKDE", level.UD = 0.95,xlim=EXT$x,ylim=EXT$y)



