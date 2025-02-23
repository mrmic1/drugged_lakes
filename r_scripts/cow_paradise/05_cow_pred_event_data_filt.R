#-----------------------------------------------------------------------------#
# Filter data to remove post-predation and mortality tracking - Cow Paradise  #
#-----------------------------------------------------------------------------#

#In this script we will remove data where individuals were tracked after being predated or had 
#an identified mortality event
#Predation events are identified in the 04_cow_pred_encounters script
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
filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"
telem_path <- "./data/telem_obj/"
enc_path <- "./data/encounters/"
save_ctmm_path = "./data/ctmm_fits/"

### DATA LOADING ###
cow_telem_data <-  readRDS(paste0(filtered_data_path, "03_lake_cow_sub.rds"))
# Load the predation event dataframe 
cow_pred_mort_events <- readRDS(paste0(enc_path, "cow_paradise/cow_pred_encounter_summary_filtered.rds"))
# Load suspected dates for pike deaths
pike_deaths <- read.csv(paste0(enc_path, "pike_deaths.csv"))



#-------------------------------------------------#
# 1. Filter out post-predation event data ####
#-------------------------------------------------#

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
cow_pred_cols <- 
  cow_pred_mort_events %>%
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated') %>% 
  mutate(
    first_date_over_50 = as.Date(first_date_over_50),
    death_date = as.Date(death_date)) %>% 
  dplyr::select(individual_ID, Species, first_date_over_50, death_date)

# Merge into single death_date column
cow_pred_cols$death_date <- ifelse(!is.na(cow_pred_cols$first_date_over_50), cow_pred_cols$first_date_over_50, cow_pred_cols$death_date)

# Ensure the merged_date column is of Date type
cow_pred_cols$death_date <- as.Date(cow_pred_cols$death_date, origin = "1970-01-01")

# Filter out data after the predation or mortality event for each individual.
cow_telem_data_2 <- 
  cow_telem_data %>%
  left_join(cow_pred_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(death_date)| Date <= death_date)  # Keep only pre-predation data
#pre-filter rows: 8860993
#post-filter rows: 8824287
#36706 rows removed

print(paste0("Rows removed after  filtering: ", nrow(cow_telem_data) - nrow(cow_telem_data_2)))

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------#
# 2. Filter out pike post-mortality data ####
#-------------------------------------------------#

#Some pike died in Lake Cow Paradise. Need to filter out date post-suspected mortality date

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
pike_mort_cols <- 
  pike_deaths %>%
  filter(individual_ID %in% cow_telem_data_2$individual_ID) %>% 
  mutate(
    pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y"))


#check
str(pike_mort_cols)

# Filter out data after the predation or mortality event for each individual.
cow_filt_data <- 
  cow_telem_data_2 %>%
  left_join(pike_mort_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(pike_death_date)| Date <= pike_death_date)  # Keep only pre-predation data
#pre-filter rows: 8824287
#post-filter rows: 8793919

print(paste0("Rows removed after  filtering: ", nrow(cow_telem_data_2) - nrow(cow_filt_data)))
#30368

#save
saveRDS(cow_filt_data, paste0(filtered_data_path, "04_lake_cow_sub.rds"))
cow_filt_data <- readRDS(paste0(filtered_data_path, "04_lake_cow_sub.rds"))

#------------------------------------------------------------#
# 3. Rerun ctmm models for individuals that were predated ####
#------------------------------------------------------------#

#I need to pull the Reference tag data again 
#load dataframe that contains this information
get_ref_data <-  readRDS(paste0(filtered_data_path, "01_lake_cow_sub.rds"))

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
  cow_pred_cols %>%
  filter(Species == 'Perch') %>% 
  pull(individual_ID)

pred_perch_data <- 
  cow_filt_data %>% 
  filter(individual_ID %in% predated_perch_ID)

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

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(4)
doParallel::registerDoParallel(cl)

# Perform model fitting using ctmm for each individual in the telemetry data.
cow_pred_perch_fits <- 
  foreach(i = 1:length(pred_perch_tel), .packages = 'ctmm') %dopar% {
    pred_perch_guess <- ctmm.guess(pred_perch_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(pred_perch_tel[[i]], pred_perch_guess, method = 'ML')
    # Save individual model fits to a specified folder.
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_cow_perch_fits/", names(pred_perch_tel)[i],".rds"))
    model_fit
  }

#> 2.2. Roach ####

#Extract predated individual from the data frame
predated_roach_IDs <- 
  cow_pred_cols %>%
  filter(Species == 'Roach') %>% 
  pull(individual_ID)

pred_roach_data <- 
  cow_filt_data %>% 
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

pred_roach_guess <- ctmm.guess(pred_roach_tel, CTMM=ctmm(error=TRUE), interactive = FALSE)
F59828_ctmm_fit <-   ctmm.fit(pred_roach_tel, pred_roach_guess, method = 'ML')
saveRDS(F59828_ctmm_fit, file = paste0(save_ctmm_path, "lake_cow_roach_fits/", "F59828.rds"))



#> 2.3. Pike ####

pike_ids <- 
  pike_deaths %>%
  filter(individual_ID %in% cow_filt_data$individual_ID) %>% 
  mutate(
    pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y")) %>%  
  pull(individual_ID)

#only 1 pike, because 893 has no data as it died soon after release

pike_mort_data <- 
  cow_filt_data %>% 
  filter(individual_ID %in% pike_ids)

pike_mort_data <- 
  with(pike_mort_data, 
       data.frame(
         "timestamp" = timestamp,                        
         "location.long" = longitude,                         
         "location.lat" = latitude, 
         "GPS.HDOP" = HDOP,                               
         "individual-local-identifier" = individual_ID))

#create telemetry object for predated roach
pike_mort_tel <- as.telemetry(pike_mort_data, 
                              timezone = "Europe/Stockholm",   
                              timeformat = "%Y-%m-%d %H:%M:%S",
                              projection = NULL,               
                              datum = "WGS84")

ctmm::projection(pike_mort_tel) <- ctmm::median(pike_mort_tel)
uere(pike_mort_tel) <- UERE

mort_pike_guess <- ctmm.guess(pike_mort_tel, CTMM=ctmm(error=TRUE), interactive = FALSE)
mort_pike_fit <- ctmm.fit(pike_mort_tel, mort_pike_guess, method = 'ML', trace = TRUE)
saveRDS(mort_pike_fit, file = paste0(save_ctmm_path, "lake_cow_pike_fits/", "F59896.rds"))

#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#

### PIKE ###

#RE_EXTRACT OUF MODELS
#RE_CREATE CTMM LIST FOR COW PARADISE WITH UPDATED CTMM MODELS

###
open_ctmm_path = "./data/ctmm_fits/lake_cow_pike_fits/"
save_ctmm_path = "./data/ctmm_fits/lake_cow_pike_fits/"

# List all RDS files in the folder
rds_files <- list.files(path = open_ctmm_path, pattern = "\\.rds$", full.names = TRUE)

# Read all RDS files into a list
rds_list <- lapply(rds_files, readRDS)

# Read all RDS files into a list
names(rds_list) <- basename(rds_files)
names(rds_list) <- sub("\\.rds$", "", names(rds_list))

# Print the names of the loaded RDS files
print(names(rds_list))

#remove pred filtered ids 
#these were rerun using ctmm.fit so we do not need to extract OUF
#OUF is already extracted. 
#Use roach_pred_ids created earlier in the script

pike_pred_ids <- c("F59896")

lake_cow_pike_ctmm_fits <- rds_list[!names(rds_list) %in% pike_pred_ids]
print(names(lake_cow_pike_ctmm_fits))

# Create a list of remaining elements
remaining_list <- rds_list[names(rds_list) %in% pike_pred_ids]
print(names(remaining_list))

# ------------------ Extract OUF models ----------------------- #

lake_cow_pike_OUF_models <- list()

#Iterate over each object in muddyfoot_roach_ctmm_fits and extract the 'OUF anisotropic error' models
lake_cow_pike_OUF_models <- lapply(lake_cow_pike_ctmm_fits, function(x) x[['OUF anisotropic error']])
#check it worked
summary(lake_cow_pike_OUF_models$F59894)

#rejoin pred filtered models
lake_cow_pike_OUF_models <- c(lake_cow_pike_OUF_models, remaining_list)[names(rds_list)]
#check that ID are in chonological order
names(lake_cow_pike_OUF_models)

#save
saveRDS(lake_cow_pike_OUF_models, paste0(save_ctmm_path, "lake_cow_pike_OUF_models.rds"))


### perch ###

#RE_EXTRACT OUF MODELS
#RE_CREATE CTMM LIST FOR COW PARADISE WITH UPDATED CTMM MODELS

###
open_ctmm_path = "./data/ctmm_fits/lake_cow_perch_fits/"
save_ctmm_path = "./data/ctmm_fits/lake_cow_perch_fits/"

# List all RDS files in the folder
rds_files <- list.files(path = open_ctmm_path, pattern = "\\.rds$", full.names = TRUE)

# Read all RDS files into a list
rds_list <- lapply(rds_files, readRDS)

# Read all RDS files into a list
names(rds_list) <- basename(rds_files)
names(rds_list) <- sub("\\.rds$", "", names(rds_list))

# Print the names of the loaded RDS files
print(names(rds_list))

#remove pred filtered ids 
#these were rerun using ctmm.fit so we do not need to extract OUF
#OUF is already extracted. 
#Use roach_pred_ids created earlier in the script

perch_pred_ids <- c("F59852", "F59878", "F59853", "F59847")

lake_cow_perch_ctmm_fits <- rds_list[!names(rds_list) %in% perch_pred_ids]
print(names(lake_cow_perch_ctmm_fits))

# Create a list of remaining elements
remaining_list <- rds_list[names(rds_list) %in% perch_pred_ids]
print(names(remaining_list))

# ------------------ Extract OUF models ----------------------- #

lake_cow_perch_OUF_models <- list()

#Iterate over each object in muddyfoot_roach_ctmm_fits and extract the 'OUF anisotropic error' models
lake_cow_perch_OUF_models <- lapply(lake_cow_perch_ctmm_fits, function(x) x[['OUF anisotropic error']])
#check it worked
summary(lake_cow_perch_OUF_models$F59839)

#rejoin pred filtered models
lake_cow_perch_OUF_models <- c(lake_cow_perch_OUF_models, remaining_list)[names(rds_list)]
#check that ID are in chonological order
names(lake_cow_perch_OUF_models)

#save
saveRDS(lake_cow_perch_OUF_models, paste0(save_ctmm_path, "lake_cow_perch_OUF_models.rds"))

### Roach ###

#RE_EXTRACT OUF MODELS
#RE_CREATE CTMM LIST FOR COW PARADISE WITH UPDATED CTMM MODELS

###
open_ctmm_path = "./data/ctmm_fits/lake_cow_roach_fits/"
save_ctmm_path = "./data/ctmm_fits/lake_cow_roach_fits/"

# List all RDS files in the folder
rds_files <- list.files(path = open_ctmm_path, pattern = "\\.rds$", full.names = TRUE)

# Read all RDS files into a list
rds_list <- lapply(rds_files, readRDS)

# Read all RDS files into a list
names(rds_list) <- basename(rds_files)
names(rds_list) <- sub("\\.rds$", "", names(rds_list))

# Print the names of the loaded RDS files
print(names(rds_list))

#remove pred filtered ids 
#these were rerun using ctmm.fit so we do not need to extract OUF
#OUF is already extracted. 
#Use roach_pred_ids created earlier in the script

roach_pred_ids <- c("F59828", "F59819")

lake_cow_roach_ctmm_fits <- rds_list[!names(rds_list) %in% roach_pred_ids]
print(names(lake_cow_roach_ctmm_fits))

# Create a list of remaining elements
remaining_list <- rds_list[names(rds_list) %in% roach_pred_ids]
print(names(remaining_list))

# ------------------ Extract OUF models ----------------------- #

lake_cow_roach_OUF_models <- list()

#Iterate over each object in muddyfoot_roach_ctmm_fits and extract the 'OUF anisotropic error' models
lake_cow_roach_OUF_models <- lapply(lake_cow_roach_ctmm_fits, function(x) x[['OUF anisotropic error']])
#check it worked
summary(lake_cow_roach_OUF_models$F59832)

#rejoin pred filtered models
lake_cow_roach_OUF_models <- c(lake_cow_roach_OUF_models, remaining_list)[names(rds_list)]
#check that ID are in chonological order
names(lake_cow_roach_OUF_models)

summary(lake_cow_roach_OUF_models$F59819)

#save
saveRDS(lake_cow_roach_OUF_models, paste0(save_ctmm_path, "lake_cow_roach_OUF_models.rds"))


