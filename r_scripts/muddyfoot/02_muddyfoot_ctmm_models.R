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
data_filter_path <- "./data/tracks_filtered/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"

#-------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
# >1. Northern Pike CTMM's ####
#-------------------------------------------------------------------------------#

# Load the dataset containing telemetry data for various species (previously filtered and pre-processed)
muddyfoot_sub <- readRDS(paste0(data_filter_path, 'muddyfoot_sub.rds'))

# Filter the data to isolate Northern Pike and the reference individual
# The reference individual is used for estimating location error (UERE)
pike_muddyfoot <- muddyfoot_sub %>%
  dplyr::filter(Species == 'Northern Pike' | individual_ID == 'Reference')

# Prepare a dataframe for conversion to a telemetry object
# This format is compatible with the movebank telemetry format
pike_movebank <- with(pike_muddyfoot, data.frame(
  "timestamp" = timestamp,                        # Timestamps for each recorded position
  "location.long" = Long,                         # Longitude coordinate
  "location.lat" = Lat,                           # Latitude coordinate
  "GPS.HDOP" = HPE,                               # Horizontal precision error from GPS
  "individual-local-identifier" = Fish_ID,        # Fish ID to identify individual fish
  "treatment" = Treatment,                        # Experimental treatment for each fish
  "date" = date,                                  # Date of observation
  "week" = week,                                  # Week number (useful for time-based grouping)
  "individual_day" = individual_day               # Unique day identifier for each individual
))

# Convert the dataframe to a telemetry object using the ctmm package's as.telemetry function
# No projection is applied here, and the WGS84 datum is used (common geographic coordinate system)
pike_muddyfoot_tel_tpeqd <- as.telemetry(pike_movebank, 
                                         timezone = "Europe/Stockholm",   # Set time zone
                                         timeformat = "%Y-%m-%d %H:%M:%S", # Specify timestamp format
                                         projection = NULL,               # No projection
                                         datum = "WGS84",                 # World Geodetic System 1984
                                         keep = c("treatment", "date", "week", "individual_day")  # Retain extra columns
)

# Check the structure of the telemetry object for one individual
# This helps ensure the conversion worked as expected
head(pike_muddyfoot_tel_tpeqd$F59880)

# Check the projection of the telemetry object (should be unprojected initially)
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
  lake_BT_pike_guess <- ctmm.guess(pike_muddyfoot_tel[[i]], CTMM = ctmm(error = TRUE), interactive = FALSE)
  
  # Fit the ctmm model to the telemetry data using model selection
  model_fit <- ctmm.select(pike_muddyfoot_tel[[i]], lake_BT_pike_guess, verbose = TRUE)
  
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
# >2. Perch CTMM's ####
#-------------------------------------------------------------------------------#

# 1. Isolate data for Perch species and a reference individual
# We are filtering the dataset 'muddyfoot_sub' to only keep rows where the species is Perch or the individual ID is 'Reference'
perch_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Perch' | individual_ID == 'Reference')

# 2. Prepare data for as.telemetry function (Movebank method)
# Create a data frame with the necessary column names for Movebank telemetry processing.
# These columns will be used by the `as.telemetry` function to convert the dataset into a telemetry object.
perch_movebank <- with(perch_muddyfoot, data.frame(
  "timestamp" = timestamp,                     # Timestamp of each observation
  "location.long" = Long,                      # Longitude of fish location
  "location.lat" = Lat,                        # Latitude of fish location
  "GPS.HDOP" = HPE,                            # Horizontal Dilution of Precision (accuracy)
  "individual-local-identifier" = Fish_ID,     # Unique identifier for each fish
  "treatment" = Treatment,                     # Treatment group
  "date" = date,                               # Date of observation
  "week" = week,                               # Week of observation
  "individual_day" = individual_day            # Individual fish day identifier
))

# 3. Convert the dataframe into a telemetry object (using equal-area projection)
# Using 'as.telemetry' without projection (TPEQD projection is set to NULL).
# Additional metadata such as 'treatment', 'date', etc., are retained.
perch_muddyfoot_tel_tpeqd <- as.telemetry(perch_movebank, 
                                          timezone = "Europe/Stockholm", 
                                          timeformat = "%Y-%m-%d %H:%M:%S", 
                                          projection = NULL,         # Projection is null since TPEQD projection is used
                                          datum = "WGS84",           # The datum is WGS84, a common standard for geographical data
                                          keep = c("treatment", "date", "week", "individual_day")  # Keep additional attributes
)

# 4. Inspect the telemetry object parameters for a specific fish
# Checking the telemetry object for the individual 'F59682'
head(perch_muddyfoot_tel_tpeqd$F59682)
ctmm::projection(perch_muddyfoot_tel_tpeqd$F59682)  # Check projection of the object
tz(perch_muddyfoot_tel_tpeqd$F59682$timestamp)      # Check timezone of timestamps

# 5. Set the projection to the geometric median of the data
# Centering the projection of all telemetry data on the geometric median of the data to avoid projection bias.
ctmm::projection(perch_muddyfoot_tel_tpeqd) <- ctmm::median(perch_muddyfoot_tel_tpeqd)

### INCORPORATING LOCATION ERROR ###

# 6. Fit UERE (User Equivalent Range Error) to reference individual
# Using reference data to estimate the error model and apply it to the telemetry data.
UERE_tpeqd <- uere.fit(perch_muddyfoot_tel_tpeqd$FReference)

# 7. Summarize and apply the UERE model
summary(UERE_tpeqd)  # Review the fitted UERE model parameters
uere(perch_muddyfoot_tel_tpeqd) <- UERE_tpeqd  # Apply the error model to the telemetry data

# 8. Remove the reference individual from the telemetry data
# Keeping only the first 30 individuals (excluding the reference individual).
perch_muddyfoot_tel <- perch_muddyfoot_tel_tpeqd[1:30]

# 9. Filter out unrealistic speeds (speeds higher than 0.977 m/s)
# 'outlie' detects possible outliers based on unrealistic speeds.
out_perch <- outlie(perch_muddyfoot_tel, plot = FALSE)
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))  # Count observations exceeding 0.977 m/s

# 10. Create a logical vector for filtering out these high-speed observations
# Retain only those observations where speed is below the critical speed threshold.
which_lowSp <- lapply(out_perch, function(x) x$speed <= 0.977)
perch_muddyfoot_tel <- Map(function(x, y) x[y, ], perch_muddyfoot_tel, which_lowSp)

# 11. Save the cleaned telemetry object
# Save the processed telemetry data to an RDS file for future use.
# saveRDS(perch_muddyfoot_tel, paste0(save_telem_path, "perch_muddyfoot_tel.rds")) 

# 12. Load telemetry object (if needed)
# perch_muddyfoot_tel <- readRDS(paste0(save_telem_path, "perch_muddyfoot_tel.rds"))

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
# >3. Roach CTMM's ####
#-------------------------------------------------------------------------------#

# Isolate Roach data from the muddyfoot_sub dataset, including reference individuals
roach_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Roach' | individual_ID == 'Reference')

# Eric's note: Movebank method. Providing as.telemetry with explicit column names makes it more efficient.
# There’s no need to coerce it into a move object first.
roach_movebank <- with(roach_muddyfoot, data.frame(
  "timestamp" = timestamp,            # Renaming to match the expected Movebank format
  "location.long" = Long,             # Longitude
  "location.lat" = Lat,               # Latitude
  "GPS.HDOP" = HPE,                   # Horizontal dilution of precision (error estimate)
  "individual-local-identifier" = Fish_ID,  # Fish ID
  "treatment" = Treatment,            # Experimental treatment
  "date" = date,                      # Date of observation
  "week" = week,                      # Week number
  "individual_day" = individual_day   # Day-specific identifier for the individual
))

# Transform the data into a telemetry object for analysis.
# Removed UTM projection for now, using TPEQD (Transverse Equidistant Projection).
roach_muddyfoot_tel_tpeqd <- as.telemetry(
  roach_movebank, 
  timezone = "Europe/Stockholm", 
  timeformat="%Y-%m-%d %H:%M:%S", 
  projection= NULL,    # No explicit projection
  datum="WGS84",       # World Geodetic System for consistency
  keep = c("treatment", "date", "week", "individual_day")  # Keeping additional columns for later use
)

# Check the contents of the telemetry object.
head(roach_muddyfoot_tel_tpeqd$F59880)   # Preview the data for Fish F59880

# Ensure the projection is properly defined (TPEQD).
ctmm::projection(roach_muddyfoot_tel_tpeqd$F59880)

# Check the timestamp timezone to ensure it's correctly set to Europe/Stockholm.
tz(roach_muddyfoot_tel_tpeqd$F59880$timestamp)

# Check column names of the telemetry object to confirm structure
names(roach_muddyfoot_tel_tpeqd)

# Center the projection on the geometric median of the data (for better alignment of coordinates).
ctmm::projection(roach_muddyfoot_tel_tpeqd) <- ctmm::median(roach_muddyfoot_tel_tpeqd)

### Incorporating location error ###

# Fit error parameters using calibration data from reference individuals.
UERE_tpeqd <- uere.fit(roach_muddyfoot_tel_tpeqd$FReference)   # FReference contains known error data

# Summarize the UERE model to inspect error estimates.
summary(UERE_tpeqd)

# Apply the error model to the entire telemetry dataset.
uere(roach_muddyfoot_tel_tpeqd) <- UERE_tpeqd

# New column `VAR.xy` added to represent error variance in x and y coordinates.
head(roach_muddyfoot_tel_tpeqd$F59684)  # Preview data after error incorporation
names(roach_muddyfoot_tel_tpeqd)

# Remove reference individuals from the list, focusing only on tracked individuals.
roach_muddyfoot_tel <- roach_muddyfoot_tel_tpeqd[1:29]

### Filtering out outliers based on unrealistic movement speeds ###

# Identify potential outliers based on fish speed (m/s). Fish swimming speeds greater than 0.841 m/s are considered outliers.
out_roach <- outlie(roach_muddyfoot_tel, plot = FALSE)

# View the speed distribution for the first individual.
head(out_roach[[1]])

# Count the number of speeds that exceed 0.841 m/s, using a review paper as reference for critical swimming speed.
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




