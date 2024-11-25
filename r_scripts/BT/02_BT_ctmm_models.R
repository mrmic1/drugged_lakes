#---------------------------------------------------------------------#
# Run CTMM models for each individual and species - Lake BT ####
# --------------------------------------------------------------------#

### LIBRARIES ###
# Load necessary libraries
library(data.table)   # For fast data manipulation
library(tidyverse)    # For data wrangling
library(ctmm)         # For continuous-time movement modeling
library(sf)           # For handling spatial data
library(parallel)     # For parallel processing
library(foreach)      # For parallel for loops
library(doParallel)   # For registering parallel backend

# Set the time zone to ensure consistent time handling
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths for reading and saving filtered telemetry and ctmm model results
filtered_data_path <- "./data/tracks_filtered/lake_BT/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"

#Load in the datasets
lake_BT_sub <- readRDS(paste0(filtered_data_path, '01_lake_BT_sub.rds'))

#-------------------------------------------------------------------------------#

#### 1. Northern pike #####

# Filter the data to isolate Northern Pike and the reference individual
# The reference individual is used for estimating location error (UERE)
pike_BT <- 
  lake_BT_sub %>%
  dplyr::filter(Species == 'Northern Pike' | individual_ID == 'FReference')
#2514194 rows


# Prepare a dataframe for conversion to a telemetry object
# This format is compatible with the movebank telemetry format
pike_movebank <- 
  with(pike_BT, 
       data.frame(
         "timestamp" = timestamp_cest,                        
         "location.long" = Long,                         
         "location.lat" = Lat, 
         "GPS.HDOP" = HPE,                               
         "individual-local-identifier" = individual_ID, 
         "Weight" = Weight,
         "Total_length" = Total_length,
         "Std_length" = Std_length,
         "Treatment" = Treatment,                        
         "Date" = date_cest,
         "Exp_Stage" = Stage,
         "Time_Of_Day" = time_of_day
       ))

# Convert the dataframe to a telemetry object using the ctmm package's as.telemetry function
# No projection is applied here, and the WGS84 datum is used (common geographic coordinate system)
pike_BT_tel_tpeqd <- as.telemetry(pike_movebank, 
                                         timezone = "Europe/Stockholm",   # Set time zone
                                         timeformat = "%Y-%m-%d %H:%M:%S",# Specify timestamp format
                                         projection = NULL,               # No projection
                                         datum = "WGS84",                 # World Geodetic System 1984
                                         keep = c("Weight","Total_length", "Std_length", "Treatment", 
                                                  "Date", "Exp_Stage", "Time_Of_Day"))  # Retain extra columns


# Check the projection of the telemetry object 
ctmm::projection(pike_BT_tel_tpeqd$F59886)
# Check the time zone of the timestamps (to confirm it's correctly set to Europe/Stockholm)
tz(pike_BT_tel_tpeqd$F59886$timestamp)

# Set a custom projection for the telemetry data
# This centers the data on the geometric median of the observed locations
# Using the median reduces the influence of extreme outliers on the projection
ctmm::projection(pike_BT_tel_tpeqd) <- ctmm::median(pike_BT_tel_tpeqd)

### Location Error Incorporation ###
# Fit UERE (user equivalent range error) for the reference individual
# UERE is used to model location error in the data, which helps improve movement model accuracy
UERE_tpeqd <- uere.fit(pike_BT_tel_tpeqd$FReference)

# Summarise the UERE to inspect the location error estimates
summary(UERE_tpeqd)

# Apply the UERE model to the entire pike telemetry dataset
# This incorporates the estimated location error into the telemetry object for each fish
uere(pike_BT_tel_tpeqd) <- UERE_tpeqd

# Remove the reference individual from the telemetry object, as we no longer need it for modeling
# Keep only the six actual pike individuals
pike_BT_tel <- pike_BT_tel_tpeqd[1:6]

### Outlier Detection and Removal ###
# Detect outliers based on unrealistic swimming speeds (greater than Ucrit speed)
# Ucrit (critical swimming speed) is a physiological limit beyond which fish are unlikely to sustain movement
out_pike <- outlie(pike_BT_tel, plot = FALSE)

# Inspect the first individual's outlier data (checking which speeds were detected as outliers)
head(out_pike[[1]])

# Summarize the number of speed outliers across all pike (speeds greater than 0.823 m/s)
sum(sapply(out_pike, function(x) sum(x$speed > 0.823)))
#11661

# Filter out observations with speeds greater than 0.823 m/s (Ucrit for Northern Pike)
# Create a logical vector indicating which speed values are below this threshold
which_lowSp <- lapply(out_pike, function(x) x$speed <= 0.823)

# Use Map to filter out the outliers for each fish
pike_BT_tel <- Map(function(x, y) x[y, ], pike_BT_tel, which_lowSp)

# Save the cleaned telemetry object (post-outlier removal)
saveRDS(pike_BT_tel, paste0(save_telem_path, "BT/pike_BT_tel.rds"))

#load object
pike_lake_BT_tel <- readRDS(paste0(save_telem_path, "pike_lake_BT_tel.rds"))

cl <- makeCluster(6)
doParallel::registerDoParallel(cl)
lake_BT_pike_select_fits <-  
  foreach(i = 1:length(pike_lake_BT_tel), .packages = 'ctmm') %dopar% {
    lake_BT_pike_guess <- ctmm.guess(pike_lake_BT_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.select(pike_lake_BT_tel[[i]], lake_BT_pike_guess, verbose = TRUE)
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_BT_pike_fits/", names(pike_lake_BT_tel)[i], ".rds"))
    model_fit
  }

stopCluster(cl)

# ###NEED TO LOOK INTO INDIVIDUAL F59889 ###
# #has 30,000 more positions than the next hightest individual.
# #I might subsample this pike in order to get the ctmm model to actually finish
# #first I will isolate the pike
# pike_F59889 <- pike_movebank %>% 
#   dplyr::filter(individual.local.identifier == 'F59889')
# library(move2)
# pike_F59889_mv <- mt_as_move2(pike_F59889, 
#                           coords = c("location.long","location.lat"),
#                           crs = "WGS84",
#                           time_column = "timestamp",
#                           track_id_column = "individual.local.identifier",
#                           na.fail = F) # allows or not empty coordinates
# 
# #Arrange by individual and datetime
# pike_F59889_mv <- pike_F59889_mv %>% 
#   dplyr::arrange(individual.local.identifier, timestamp)
# 
# pike_F59889_sub <- pike_F59889_mv %>% mt_filter_per_interval(unit = "4 seconds",  criterion="first")
# 
# pike_F59889_sub <- to_move(pike_F59889_sub)
# 
# pike_F59889_tel <- as.telemetry(pike_F59889_sub, 
#                                        timezone = "Europe/Stockholm", 
#                                        timeformat="%Y-%m-%d %H:%M:%S", 
#                                        projection= NULL,
#                                        datum="WGS84",
#                                        keep = c("treatment", "date", 
#                                                 "week", "individual_day")
# )
# 
# 
# ctmm::projection(pike_F59889_tel)
# #tpeqd projection
# tz(pike_F59889_tel$timestamp)
# #"Europe/Stockholm"
# 
# #Center the projection on the geometric median of the data
# ctmm::projection(pike_F59889_tel) <- ctmm::median(pike_F59889_tel)
# 
# #summary(UERE_utm)
# summary(UERE_tpeqd)
# #both are similar
# 
# # apply error model to data
# #uere(pike_lake_BT_tel_utm) <- UERE_utm
# uere(pike_F59889_tel) <- UERE_tpeqd
# #new column now called VAR.xy
# head(pike_F59889_tel)
# 
# #remove outliers based on speed
# out_pike <- outlie(pike_F59889_tel, plot = FALSE)
# head(out_pike[[1]])
# sum(out_pike$speed > 0.823)
# 
# #filter out speeds
# pike_F59889_tel <- pike_F59889_tel[out_pike$speed < 0.823, ]
# 
# 
# lake_BT_pike_guess <- ctmm.guess(pike_F59889_tel, CTMM=ctmm(error=TRUE), interactive = FALSE)
# F59889 <- ctmm.select(pike_F59889_tel, lake_BT_pike_guess, verbose = TRUE, trace = TRUE)
# saveRDS(F59889, file = paste0(save_ctmm_path, "lake_BT_pike_fits/", "F59889.rds"))
 
#---------------------------------------------------------------------------------#

#### 2. Perch ####

# Isolate data for Perch species and a reference individual
# We are filtering the dataset 'BT_sub' to only keep rows where the species is Perch or the individual ID is 'Reference'
perch_BT <- 
  lake_BT_sub %>% 
  dplyr::filter(Species == 'Perch' | individual_ID == 'FReference')

# 2. Prepare data for as.telemetry function (Movebank method)
# Create a data frame with the necessary column names for Movebank telemetry processing.
# These columns will be used by the `as.telemetry` function to convert the dataset into a telemetry object.
perch_movebank <-   with(perch_BT, 
                         data.frame(
                           "timestamp" = timestamp_cest,                        
                           "location.long" = Long,                         
                           "location.lat" = Lat, 
                           "GPS.HDOP" = HPE,                               
                           "individual-local-identifier" = individual_ID, 
                           "Weight" = Weight,
                           "Total_length" = Total_length,
                           "Std_length" = Std_length,
                           "Treatment" = Treatment,                        
                           "Date" = date_cest,
                           "Exp_Stage" = Stage,
                           "Time_Of_Day" = time_of_day
                         ))

perch_BT_tel_tpeqd <- as.telemetry(perch_movebank, 
                                          timezone = "Europe/Stockholm", 
                                          timeformat = "%Y-%m-%d %H:%M:%S", 
                                          projection = NULL,         
                                          datum = "WGS84",           
                                          keep = c("Weight","Total_length", "Std_length", "Treatment", 
                                                   "Date", "Exp_Stage", "Time_Of_Day"))


# Inspect the telemetry object parameters for a specific fish
ctmm::projection(perch_BT_tel_tpeqd$F59749)  # Check projection of the object
tz(perch_BT_tel_tpeqd$F59749$timestamp)      # Check timezone of timestamps

# Centering the projection of all telemetry data on the geometric median of the data to avoid projection bias.
ctmm::projection(perch_BT_tel_tpeqd) <- ctmm::median(perch_BT_tel_tpeqd)

# Fit UERE (User Equivalent Range Error) to reference individual
# Using reference data to estimate the error model and apply it to the telemetry data.
UERE_tpeqd <- uere.fit(perch_BT_tel_tpeqd$FReference)
summary(UERE_tpeqd)  # Review the fitted UERE model parameters
uere(perch_BT_tel_tpeqd) <- UERE_tpeqd  # Apply the error model to the telemetry data

# Keeping only the first 30 individuals (excluding the reference individual).
perch_BT_tel <- perch_BT_tel_tpeqd[1:30]
names(perch_BT_tel)

#Filter out unrealistic speeds (speeds higher than 0.977 m/s)
out_perch <- outlie(perch_BT_tel, plot = FALSE)
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))  # Count observations exceeding 0.977 m/s
#34023
which_lowSp <- lapply(out_perch, function(x) x$speed <= 0.977)
perch_BT_tel <- Map(function(x, y) x[y, ], perch_BT_tel, which_lowSp)

# Save the cleaned telemetry object
saveRDS(perch_BT_tel, paste0(save_telem_path, "BT/perch_BT_tel.rds"))

cl <- makeCluster(20)
doParallel::registerDoParallel(cl)
lake_BT_perch_select_fits <-  
  foreach(i = 1:length(perch_lake_BT_tel), .packages = 'ctmm') %dopar% {
    lake_BT_perch_guess <- ctmm.guess(perch_lake_BT_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.select(perch_lake_BT_tel[[i]], lake_BT_perch_guess, verbose = TRUE)
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_BT_perch_fits/", names(perch_lake_BT_tel)[i], ".rds"))
    model_fit
  }

stopCluster(cl)

# ###NEED TO LOOK INTO INDIVIDUAL F59789 ###
# #has 30,000 more positions than the next hightest individual.
# #I might subsample this perch in order to get the ctmm model to actually finish
# #first I will isolate the perch
# 
# print(perch_movebank %>%
#         group_by(individual.local.identifier) %>% 
#         summarise(count = n()), n = 30)
# #has a relativery low number of observations
# 
# perch_F59789 <- perch_movebank %>% 
#   dplyr::filter(individual.local.identifier == 'F59789')
# 
# 
# perch_F59789_tel <- as.telemetry(perch_F59789, 
#                                  timezone = "Europe/Stockholm", 
#                                  timeformat="%Y-%m-%d %H:%M:%S", 
#                                  projection= NULL,
#                                  datum="WGS84",
#                                  keep = c("treatment", "date", 
#                                           "week", "individual_day")
# )
# 
# 
# ctmm::projection(perch_F59789_tel)
# #tpeqd projection
# tz(perch_F59789_tel$timestamp)
# #"Europe/Stockholm"
# 
# #Center the projection on the geometric median of the data
# ctmm::projection(perch_F59789_tel) <- ctmm::median(perch_F59789_tel)
# 
# #summary(UERE_utm)
# summary(UERE_tpeqd)
# #both are similar
# 
# # apply error model to data
# #uere(perch_lake_BT_tel_utm) <- UERE_utm
# uere(perch_F59789_tel) <- UERE_tpeqd
# #new column now called VAR.xy
# head(perch_F59789_tel)
# 
# #remove outliers based on speed
# out_perch <- outlie(perch_F59789_tel, plot = FALSE)
# head(out_perch[[1]])
# sum(perch_F59789_tel$speed > 0.977)
# 
# #filter out speeds
# perch_F59789_tel <- perch_F59789_tel[out_perch$speed < 0.977, ]
# 
# 
# lake_BT_perch_guess <- ctmm.guess(perch_F59789_tel, CTMM=ctmm(error=TRUE), interactive = FALSE)
# F59789 <- ctmm.select(perch_F59789_tel, lake_BT_perch_guess, verbose = TRUE)
# saveRDS(F59789, file = paste0(save_ctmm_path, "lake_BT_perch_fits/", "F59789.rds"))

#------------------------------------------------------------------------------------------#

#### 3. Roach ####

# Isolate Roach data from the muddyfoot_sub dataset, including reference individuals
roach_BT <- 
  lake_BT_sub %>% 
  dplyr::filter(Species == 'Roach' | individual_ID == 'FReference')
#13596887

roach_movebank <- with(roach_BT, 
                       data.frame(
                         "timestamp" = timestamp_cest,                        
                         "location.long" = Long,                         
                         "location.lat" = Lat, 
                         "GPS.HDOP" = HPE,                               
                         "individual-local-identifier" = individual_ID, 
                         "Weight" = Weight,
                         "Total_length" = Total_length,
                         "Std_length" = Std_length,
                         "Treatment" = Treatment,                        
                         "Date" = date_cest,
                         "Exp_Stage" = Stage,
                         "Time_Of_Day" = time_of_day
                       ))

roach_BT_tel_tpeqd <- as.telemetry(
  roach_movebank, 
  timezone = "Europe/Stockholm", 
  timeformat="%Y-%m-%d %H:%M:%S", 
  projection= NULL,    
  datum="WGS84",       
  keep = c("Weight","Total_length", "Std_length", "Treatment", 
           "Date", "Exp_Stage", "Time_Of_Day"))


# Check the contents of the telemetry object.
ctmm::projection(roach_BT_tel_tpeqd$F59766)
tz(roach_BT_tel_tpeqd$F59766$timestamp)

# Center the projection on the geometric median of the data (for better alignment of coordinates).
ctmm::projection(roach_BT_tel_tpeqd) <- ctmm::median(roach_BT_tel_tpeqd)

# Fit error parameters using calibration data from reference tag.
UERE_tpeqd <- uere.fit(roach_BT_tel_tpeqd$FReference)

# Summarise the UERE model to inspect error estimates.
summary(UERE_tpeqd)

# Apply the error model to the entire telemetry dataset.
uere(roach_BT_tel_tpeqd) <- UERE_tpeqd

# Remove reference tag from the list, focusing only on tracked individuals.
roach_BT_tel <- roach_BT_tel_tpeqd[1:30]
names(roach_BT_tel)

# Identify potential outliers based on fish speed (m/s). Fish swimming speeds greater than 0.841 m/s are considered outliers.
out_roach <- outlie(roach_BT_tel, plot = FALSE)

# Count the number of speeds that exceed 0.841 m/s.
sum(sapply(out_roach, function(x) sum(x$speed > 0.841)))
# 114783 observations exceed the threshold of 0.841 m/s.

# Create a logical vector to identify observations with speeds below 0.841 m/s.
which_lowSp <- lapply(out_roach, function(x) x$speed <= 0.841)

# Filter out high-speed outliers from the telemetry data.
roach_BT_tel <- Map(function(x, y) x[y,], roach_BT_tel, which_lowSp)

### Saving the telemetry object for future use ###

# Save the cleaned telemetry object to a specified path.
saveRDS(roach_BT_tel, paste0(save_telem_path, "BT/roach_BT_tel.rds"))

### Model fitting using parallel processing ###

# Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(28)
doParallel::registerDoParallel(cl)

# Perform model fitting using ctmm for each individual in the telemetry data.
BT_roach_select_fits <- foreach(i = 1:length(roach_BT_tel), .packages = 'ctmm') %dopar% {
  BT_roach_guess <- ctmm.guess(roach_BT_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
  model_fit <- ctmm.select(roach_BT_tel[[i]], BT_roach_guess, verbose = TRUE)
  
  # Save individual model fits to a specified folder.
  saveRDS(model_fit, file = paste0(save_ctmm_path, "BT_roach_fits/", names(roach_BT_tel)[i], ".rds"))
  model_fit
}

# Stop the parallel cluster after the computation is complete.
stopCluster(cl)

# Assign names to the model fits for easier reference.
names(BT_roach_select_fits) <- names(roach_BT_tel)

# Save the full list of model fits to a file for later use.
saveRDS(BT_roach_select_fits, file = paste0(save_ctmm_path, "BT_roach_fits/", "BT_roach_ctmm_fits.rds"))

# To reload the saved models, you can use:
# mud_roach_fits <- readRDS(paste0(save_ctmm_path, "BT_roach_ctmm_fits.rds"))


#-----------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
#### 4. Make combined dataframe with outliers for each species removed ####
#-------------------------------------------------------------------------#

#If you need to reload the telemetry objects
pike_BT_tel <- readRDS(paste0(save_telem_path, "BT/pike_BT_tel.rds"))
perch_BT_tel <- readRDS(paste0(save_telem_path, "BT/perch_BT_tel.rds"))
roach_BT_tel <- readRDS(paste0(save_telem_path, "BT/roach_BT_tel.rds"))


# Create a function to extract the data and add the individual_id column
# This function merges telemetry data from different individuals into a single dataframe.
extract_data <- function(list_data) {
  data_combined <- do.call(rbind, lapply(names(list_data), function(id) {
    df <- list_data[[id]]
    df$individual_ID <- id  # Add individual ID to each dataframe.
    return(df)
  }))
  return(data_combined)
}

# Apply the function to each species' telemetry data and combine them into a single dataframe.
pike_data <- extract_data(pike_BT_tel)
perch_data <- extract_data(perch_BT_tel)
roach_data <- extract_data(roach_BT_tel)

# Combine all dataframes into one, with a species label added.
BT_telem_data <- 
  rbind(
    cbind(pike_data, Species = 'Pike'), 
    cbind(perch_data, Species = 'Perch'),
    cbind(roach_data, Species = 'Roach')
  )

saveRDS(BT_telem_data, paste0(filtered_data_path, "02_lake_BT_sub.rds"))
