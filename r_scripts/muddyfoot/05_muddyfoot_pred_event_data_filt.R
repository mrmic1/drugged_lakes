#------------------------------------------------#
# Filter data to remove post-predation tracking  #
#------------------------------------------------#

#In this script we will remove data where individuals were tracked after being predated
#Predation events are identified in the 04_muddyfoot_pred_encounters script
#We will then re-run ctmm movement models for these individuals

### LIBRARIES ###
# Loading required packages
library(dplyr)
library(ctmm)
library(data.table)
library(sf)  # For spatial data handling
library(parallel)     # For parallel processing
library(foreach)      # For parallel for loops
library(doParallel)   # For registering parallel backend

# Set the time zone environment to 'Europe/Stockholm' for consistent timestamp manipulation
Sys.setenv(TZ = 'Europe/Stockholm')

### DIRECTORIES ###
# Paths to directories containing the necessary data files
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
telem_path <- "./data/telem_obj/"
enc_path <- "./data/encounters/"
save_ctmm_path = "./data/ctmm_fits/"

### DATA LOADING ###
muddyfoot_telem_data <-  readRDS(paste0(filtered_data_path, "03_muddyfoot_sub.rds"))
# Load the predation event dataframe 
mud_pred_events <- readRDS(paste0(enc_path, "muddyfoot_pred_events.rds"))

#-------------------------------------------------#
# 1. Filter out post-predation event data ####
#-------------------------------------------------#

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
pred_cols <- 
  mud_pred_events %>%
  filter(likely_predated == 1) %>% 
  dplyr::select(individual_ID, first_date_over_50)

# Filter out data after the predation event for each individual.
muddyfoot_filt_data <- 
  muddyfoot_telem_data %>%
  left_join(pred_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(first_date_over_50) | Date <= first_date_over_50)  # Keep only pre-predation data
#pre-filter rows: 10358641
#post-filter rows: 10067142
#291499 rows removed

print(paste0("Rows removed after predation filtering: ", nrow(muddyfoot_telem_data) - nrow(muddyfoot_filt_data)))

#save
saveRDS(muddyfoot_filt_data, paste0(filtered_data_path, "04_muddyfoot_sub.rds"))


#------------------------------------------------------------#
# 2. Rerun ctmm models for individuals that were predated ####
#------------------------------------------------------------#

#I need to pull the Reference tag data again 
#load dataframe that contains this information
get_ref_data <-  readRDS(paste0(filtered_data_path, "01_muddyfoot_sub.rds"))

ref_data <- get_ref_data %>% 
  filter(individual_ID == 'FReference') %>% 
  dplyr::select(individual_ID, HPE, timestamp_cest, Long, Lat)


ref_data <- with(ref_data, 
                 data.frame(
                   "timestamp" = timestamp_cest,                        
                   "location.long" = Long,                         
                   "location.lat" = Lat, 
                   "GPS.HDOP" = HPE,                               
                   "individual-local-identifier" = individual_ID
                 ))

#create telemetry object for the reference data
ref_tel <- as.telemetry(ref_data, 
                        timezone = "Europe/Stockholm",   
                        timeformat = "%Y-%m-%d %H:%M:%S",
                        projection = NULL,               
                        datum = "WGS84")

ctmm::projection(ref_tel) <- ctmm::median(ref_tel)
UERE <- uere.fit(ref_tel)
summary(UERE)


#> 2.1. Perch ####

#Extract predated individual from the data frame
predated_perch_ID <- 
  mud_pred_events %>%
  filter(likely_predated == 1 & Species == 'Perch') %>% 
  pull(individual_ID)

pred_perch_data <- 
  muddyfoot_filt_data %>% 
  filter(individual_ID == predated_perch_ID)

pred_perch_data <- 
  with(pred_perch_data, 
       data.frame(
         "timestamp" = timestamp,                        
         "location.long" = longitude,                         
         "location.lat" = latitude, 
         "GPS.HDOP" = HDOP,                               
         "individual-local-identifier" = individual_ID))

#create telemetry object for predated perch
pred_perch_tel <- as.telemetry(pred_perch_data, 
                        timezone = "Europe/Stockholm",   
                        timeformat = "%Y-%m-%d %H:%M:%S",
                        projection = NULL,               
                        datum = "WGS84")

ctmm::projection(pred_perch_tel) <- ctmm::median(pred_perch_tel)
uere(pred_perch_tel) <- UERE

pred_perch_guess <- ctmm.guess(pred_perch_tel, CTMM=ctmm(error=TRUE), interactive = FALSE)
pred_perch_fit <- ctmm.fit(pred_perch_tel, pred_perch_guess, method = 'ML', trace = TRUE)
saveRDS(pred_perch_fit, file = paste0(save_ctmm_path, "muddyfoot_perch_fits/", "F59692.rds"))

#> 2.2. Roach ####

#Extract predated individual from the data frame
predated_roach_IDs <- 
  mud_pred_events %>%
  filter(likely_predated == 1 & Species == 'Roach') %>% 
  pull(individual_ID)

pred_roach_data <- 
  muddyfoot_filt_data %>% 
  filter(individual_ID %in% predated_roach_IDs)

pred_roach_data <- 
  with(pred_roach_data, 
       data.frame(
         "timestamp" = timestamp,                        
         "location.long" = longitude,                         
         "location.lat" = latitude, 
         "GPS.HDOP" = HDOP,                               
         "individual-local-identifier" = individual_ID))

#create telemetry object for predated perch
pred_roach_tel <- as.telemetry(pred_roach_data, 
                               timezone = "Europe/Stockholm",   
                               timeformat = "%Y-%m-%d %H:%M:%S",
                               projection = NULL,               
                               datum = "WGS84")

ctmm::projection(pred_roach_tel) <- ctmm::median(pred_roach_tel)
uere(pred_roach_tel) <- UERE

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(5)
doParallel::registerDoParallel(cl)

# Perform model fitting using ctmm for each individual in the telemetry data.
pred_roach_fits <- 
  foreach(i = 1:length(pred_roach_tel), .packages = 'ctmm') %dopar% {
  pred_roach_guess <- ctmm.guess(pred_roach_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
  model_fit <- ctmm.fit(pred_roach_tel[[i]], pred_roach_guess, method = 'ML')
  # Save individual model fits to a specified folder.
  saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_roach_fits/", names(pred_roach_tel)[i], ".rds"))
  model_fit
}

