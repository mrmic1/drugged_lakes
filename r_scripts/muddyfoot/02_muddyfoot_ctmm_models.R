#---------------------------------------------------------------------#
# Run CTMM models for each individual and species - Lake Muddyfoot ####
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
filtered_data_path <- "./data/tracks_filtered/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"

#-------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
# > 1. Northern Pike CTMM's ####
#-------------------------------------------------------------------------------#

# Load the dataset pre-processed in script 01_muddyfoot_data_cleaning.R
muddyfoot_sub <- readRDS(paste0(filtered_data_path, 'muddyfoot/01_muddyfoot_sub.rds'))

# Filter the data to isolate Northern Pike and the reference individual
# The reference individual is used for estimating location error (UERE)
pike_muddyfoot <- muddyfoot_sub %>%
  dplyr::filter(Species == 'Northern Pike' | individual_ID == 'FReference')
#1868158 rows


# Prepare a dataframe for conversion to a telemetry object
# This format is compatible with the movebank telemetry format
pike_movebank <- 
  with(pike_muddyfoot, 
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
pike_muddyfoot_tel_tpeqd <- as.telemetry(pike_movebank, 
                                         timezone = "Europe/Stockholm",   # Set time zone
                                         timeformat = "%Y-%m-%d %H:%M:%S",# Specify timestamp format
                                         projection = NULL,               # No projection
                                         datum = "WGS84",                 # World Geodetic System 1984
                                         keep = c("Weight","Total_length", "Std_length", "Treatment", 
                                                  "Date", "Exp_Stage", "Time_Of_Day"))  # Retain extra columns


# Check the structure of the telemetry object for one individual
# This helps ensure the conversion worked as expected
head(pike_muddyfoot_tel_tpeqd$F59880)

# Check the projection of the telemetry object 
ctmm::projection(pike_muddyfoot_tel_tpeqd$F59880)

# Check the time zone of the timestamps (to confirm it's correctly set to Europe/Stockholm)
tz(pike_muddyfoot_tel_tpeqd$F59880$timestamp)

# Set a custom projection for the telemetry data
# This centers the data on the geometric median of the observed locations
# Using the median reduces the influence of extreme outliers on the projection
ctmm::projection(pike_muddyfoot_tel_tpeqd) <- ctmm::median(pike_muddyfoot_tel_tpeqd)

### Location Error Incorporation ###
# Fit UERE (user equivalent range error) for the reference individual
# UERE is used to model location error in the data, which helps improve movement model accuracy
UERE_tpeqd <- uere.fit(pike_muddyfoot_tel_tpeqd$FReference)

# Summarize the UERE to inspect the location error estimates

summary(UERE_tpeqd)

# Apply the UERE model to the entire pike telemetry dataset
# This incorporates the estimated location error into the telemetry object for each fish
uere(pike_muddyfoot_tel_tpeqd) <- UERE_tpeqd

# Remove the reference individual from the telemetry object, as we no longer need it for modeling
# Keep only the six actual pike individuals
pike_muddyfoot_tel <- pike_muddyfoot_tel_tpeqd[1:6]

### Outlier Detection and Removal ###
# Detect outliers based on unrealistic swimming speeds (greater than Ucrit speed)
# Ucrit (critical swimming speed) is a physiological limit beyond which fish are unlikely to sustain movement
out_pike <- outlie(pike_muddyfoot_tel, plot = FALSE)

# Inspect the first individual's outlier data (checking which speeds were detected as outliers)
head(out_pike[[1]])

# Summarize the number of speed outliers across all pike (speeds greater than 0.823 m/s)
sum(sapply(out_pike, function(x) sum(x$speed > 0.823)))
#123469

# Filter out observations with speeds greater than 0.823 m/s (Ucrit for Northern Pike)
# Create a logical vector indicating which speed values are below this threshold
which_lowSp <- lapply(out_pike, function(x) x$speed <= 0.823)

# Use Map to filter out the outliers for each fish
pike_muddyfoot_tel <- Map(function(x, y) x[y, ], pike_muddyfoot_tel, which_lowSp)

# Save the cleaned telemetry object (post-outlier removal)
saveRDS(pike_muddyfoot_tel, paste0(save_telem_path, "pike_muddyfoot_tel.rds"))

### CTMM Model Fitting ###
# Run ctmm models for each pike using parallel processing
cl <- makeCluster(6)  # Initialize a parallel cluster with 6 cores (one per fish)
doParallel::registerDoParallel(cl)  # Register the parallel backend

# Fit ctmm models in parallel for each pike
muddyfoot_pike_ctmm_fits <- foreach(i = 1:length(pike_muddyfoot_tel), .packages = 'ctmm') %dopar% {
  # Generate an initial guess for the model parameters, incorporating location error
 muddyfoot_pike_guess <- ctmm.guess(pike_muddyfoot_tel[[i]], CTMM = ctmm(error = TRUE), interactive = FALSE)
  
  # Fit the ctmm model to the telemetry data using model selection
  model_fit <- ctmm.select(pike_muddyfoot_tel[[i]], muddyfoot_pike_guess, verbose = TRUE)
  
  # Save the fitted model for each fish individually
  saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_pike_fits/", names(pike_muddyfoot_tel)[i], ".rds"))
  
  model_fit  # Return the fitted model
}

# Stop the parallel cluster once model fitting is complete
stopCluster(cl)

# Assign names to the list of ctmm model fits (each element corresponds to a fish ID)
names(muddyfoot_pike_ctmm_fits) <- names(pike_muddyfoot_tel)

# Save the entire list of model fits to disk
saveRDS(muddyfoot_pike_ctmm_fits, paste0(save_ctmm_path, "muddyfoot_pike_ctmm_fits.rds"))

### Model Output Checking ###
# Loop through each fitted model and print a summary
# This provides a quick diagnostic of the fitted models for each pike
for (i in 1:length(muddyfoot_pike_ctmm_fits)) {
  cat("Summary for", names(muddyfoot_pike_ctmm_fits)[i], ":\n")
  print(summary(muddyfoot_pike_ctmm_fits[[i]]))
  cat("\n")
}

#-------------------------------------------------------------------------------#
# > 2. Perch CTMM's ####
#-------------------------------------------------------------------------------#

# Isolate data for Perch species and a reference individual
# We are filtering the dataset 'muddyfoot_sub' to only keep rows where the species is Perch or the individual ID is 'Reference'
perch_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Perch' | individual_ID == 'FReference')

# 2. Prepare data for as.telemetry function (Movebank method)
# Create a data frame with the necessary column names for Movebank telemetry processing.
# These columns will be used by the `as.telemetry` function to convert the dataset into a telemetry object.
perch_movebank <-   with(perch_muddyfoot, 
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

perch_muddyfoot_tel_tpeqd <- as.telemetry(perch_movebank, 
                                          timezone = "Europe/Stockholm", 
                                          timeformat = "%Y-%m-%d %H:%M:%S", 
                                          projection = NULL,         
                                          datum = "WGS84",           
                                          keep = c("Weight","Total_length", "Std_length", "Treatment", 
                                                   "Date", "Exp_Stage", "Time_Of_Day"))


# Inspect the telemetry object parameters for a specific fish
head(perch_muddyfoot_tel_tpeqd$F59682)
ctmm::projection(perch_muddyfoot_tel_tpeqd$F59682)  # Check projection of the object
tz(perch_muddyfoot_tel_tpeqd$F59682$timestamp)      # Check timezone of timestamps

# Centering the projection of all telemetry data on the geometric median of the data to avoid projection bias.
ctmm::projection(perch_muddyfoot_tel_tpeqd) <- ctmm::median(perch_muddyfoot_tel_tpeqd)

# Fit UERE (User Equivalent Range Error) to reference individual
# Using reference data to estimate the error model and apply it to the telemetry data.
UERE_tpeqd <- uere.fit(perch_muddyfoot_tel_tpeqd$FReference)
summary(UERE_tpeqd)  # Review the fitted UERE model parameters
uere(perch_muddyfoot_tel_tpeqd) <- UERE_tpeqd  # Apply the error model to the telemetry data

# Keeping only the first 30 individuals (excluding the reference individual).
perch_muddyfoot_tel <- perch_muddyfoot_tel_tpeqd[1:30]

#Filter out unrealistic speeds (speeds higher than 0.977 m/s)
out_perch <- outlie(perch_muddyfoot_tel, plot = FALSE)
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))  # Count observations exceeding 0.977 m/s
#156265
which_lowSp <- lapply(out_perch, function(x) x$speed <= 0.977)
perch_muddyfoot_tel <- Map(function(x, y) x[y, ], perch_muddyfoot_tel, which_lowSp)

# Save the cleaned telemetry object
saveRDS(perch_muddyfoot_tel, paste0(save_telem_path, "perch_muddyfoot_tel.rds")) 

### MODEL FITTING ###

# 13. Filter specific individuals based on their IDs
# Filter the telemetry data for individuals that we are interested in running the ctmm models on.
desired_identities <- c("F59682", "F59688", "F59689", "F59698", "F59699", "F59708",
                        "F59711", "F59712", "F59714", "F59717", "F59728", "F59730",
                        "F59733", "F59734")

object_names <- names(perch_muddyfoot_tel)
extracted_names <- sub(".*\\$", "", object_names)  # Extract IDs from object names
filtered_list <- perch_muddyfoot_tel[extracted_names %in% desired_identities]  # Filter by desired IDs

# 14. Run ctmm models in parallel for each perch individual
# Fitting OUF models for all filtered individuals in parallel.
cl <- makeCluster(14)  # Create a cluster with 14 cores
doParallel::registerDoParallel(cl)

muddyfoot_perch_ctmm_fits_OUF <- foreach(i = 1:length(filtered_list), .packages = 'ctmm') %dopar% {
  muddyfoot_perch_guess <- ctmm.guess(filtered_list[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
  model_fit <- ctmm.select(filtered_list[[i]], muddyfoot_perch_guess, verbose = TRUE)
  saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_perch_fits/", names(filtered_list)[i], "_OUF.rds"))
  model_fit
}

stopCluster(cl)  # Stop the parallel cluster

# 15. Add individual IDs to model results
# Assigning names to the model fits based on their respective individuals.
names(muddyfoot_perch_ctmm_fits_OUF) <- names(filtered_list)

# 16. Save the ctmm model fits
# Save the OUF model fits for future analysis.
saveRDS(muddyfoot_perch_ctmm_fits_OUF, paste0(save_ctmm_path, "muddyfoot_perch_ctmm_fits_OUF.rds")) 

# 17. Check model outputs
# Loop through each model fit and print a summary for review.
for (i in 1:length(muddyfoot_perch_select_fits)) {
  element_name <- names(muddyfoot_perch_select_fits)[i]
  cat("Summary for", element_name, ":\n")
  print(summary(muddyfoot_perch_select_fits[[i]]))
  cat("\n")
}

#-------------------------------------------------------------------------------#
# > 3. Roach CTMM's ####
#-------------------------------------------------------------------------------#

# Isolate Roach data from the muddyfoot_sub dataset, including reference individuals
roach_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Roach' | individual_ID == 'FReference')

roach_movebank <- with(roach_muddyfoot, 
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

roach_muddyfoot_tel_tpeqd <- as.telemetry(
  roach_movebank, 
  timezone = "Europe/Stockholm", 
  timeformat="%Y-%m-%d %H:%M:%S", 
  projection= NULL,    
  datum="WGS84",       
  keep = c("Weight","Total_length", "Std_length", "Treatment", 
           "Date", "Exp_Stage", "Time_Of_Day"))


# Check the contents of the telemetry object.
head(roach_muddyfoot_tel_tpeqd$F59683)
ctmm::projection(roach_muddyfoot_tel_tpeqd$F59683)
tz(roach_muddyfoot_tel_tpeqd$F59683$timestamp)

# Center the projection on the geometric median of the data (for better alignment of coordinates).
ctmm::projection(roach_muddyfoot_tel_tpeqd) <- ctmm::median(roach_muddyfoot_tel_tpeqd)

# Fit error parameters using calibration data from reference tag.
UERE_tpeqd <- uere.fit(roach_muddyfoot_tel_tpeqd$FReference)

# Summarise the UERE model to inspect error estimates.
summary(UERE_tpeqd)

# Apply the error model to the entire telemetry dataset.
uere(roach_muddyfoot_tel_tpeqd) <- UERE_tpeqd

# Remove reference tag from the list, focusing only on tracked individuals.
roach_muddyfoot_tel <- roach_muddyfoot_tel_tpeqd[1:29]

# Identify potential outliers based on fish speed (m/s). Fish swimming speeds greater than 0.841 m/s are considered outliers.
out_roach <- outlie(roach_muddyfoot_tel, plot = FALSE)

# Count the number of speeds that exceed 0.841 m/s.
sum(sapply(out_roach, function(x) sum(x$speed > 0.841)))
# 336192 observations exceed the threshold of 0.841 m/s.

# Create a logical vector to identify observations with speeds below 0.841 m/s.
which_lowSp <- lapply(out_roach, function(x) x$speed <= 0.841)

# Filter out high-speed outliers from the telemetry data.
roach_muddyfoot_tel <- Map(function(x, y) x[y,], roach_muddyfoot_tel, which_lowSp)

### Saving the telemetry object for future use ###

# Save the cleaned telemetry object to a specified path.
saveRDS(roach_muddyfoot_tel, paste0(save_telem_path, "roach_muddyfoot_tel.rds"))

### Model fitting using parallel processing ###

# Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(28)
doParallel::registerDoParallel(cl)

# Perform model fitting using ctmm for each individual in the telemetry data.
muddyfoot_roach_select_fits <- foreach(i = 1:length(roach_muddyfoot_tel), .packages = 'ctmm') %dopar% {
  muddyfoot_roach_guess <- ctmm.guess(roach_muddyfoot_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
  model_fit <- ctmm.select(roach_muddyfoot_tel[[i]], muddyfoot_roach_guess, verbose = TRUE)
  
  # Save individual model fits to a specified folder.
  saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_roach_fits/", names(roach_muddyfoot_tel)[i], ".rds"))
  model_fit
}

# Stop the parallel cluster after the computation is complete.
stopCluster(cl)

# Assign names to the model fits for easier reference.
names(muddyfoot_roach_select_fits) <- names(roach_muddyfoot_tel)

# Save the full list of model fits to a file for later use.
saveRDS(muddyfoot_roach_select_fits, file = paste0(save_ctmm_path, "muddyfoot_roach_fits/", "muddyfoot_roach_ctmm_fits.rds"))

# To reload the saved models, you can use:
# mud_roach_fits <- readRDS(paste0(save_ctmm_path, "muddyfoot_roach_ctmm_fits.rds"))


#---------------------------------------------------------------#####
# Make combined dataframe with outliers for each species removed ####
#---------------------------------------------------------------#####

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
pike_data <- extract_data(pike_muddyfoot_tel)
perch_data <- extract_data(perch_muddyfoot_tel)
roach_data <- extract_data(roach_muddyfoot_tel)

# Combine all dataframes into one, with a species label added.
muddyfoot_telem_data <- 
  rbind(
    cbind(pike_data, Species = 'Pike'), 
    cbind(perch_data, Species = 'Perch'),
    cbind(roach_data, Species = 'Roach')
  )

saveRDS(muddyfoot_telem_data, paste0(filtered_data_path, "muddyfoot/02_muddyfoot_sub.rds"))

