#-------------------------------------------------------------------#
# Filter data to remove post-predation and mortality tracking - BT  #
#-------------------------------------------------------------------#

#SCRIPT DESCRIPTION

#In this script we will remove data where individuals were tracked after being predated or had 
#an identified mortality event
#Predation events are identified in the `04_BT_calculate_pred_prey_ints` script
#We will then re-run ctmm movement models for individuals that had their data filtered

#----------------------------------------------------------------------------------------------#

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
filtered_data_path <- "./data/tracks_filtered/lake_BT/"
telem_path <- "./data/telem_obj/BT/"
enc_path <- "./data/encounters/BT/"
save_ctmm_path = "./data/ctmm_fits/"

### DATA LOADING ###
BT_telem_data <-  readRDS(paste0(filtered_data_path, "03_lake_BT_sub.rds"))
# Load the predation event dataframe 
mortality_preds <- readxl::read_excel("./data/encounters/suspected_mortality_updated.xlsx")
# Load suspected dates for pike deaths
pike_deaths <- read.csv("./data/encounters/pike_deaths.csv")

#------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------#
# 1. Filter out post-predation and mortality tracking locations #####
#-------------------------------------------------------------------#

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
BT_pred_prey_cols <-
  mortality_preds %>%
  filter(lake == 'BT') %>% 
  filter(species == 'Roach'| species == 'Perch') %>%
  dplyr::select(individual_ID, species, revised_suspected_mortality, revised_likely_death_date) %>% 
  rename(death_date = revised_likely_death_date)

# Ensure the merged_date column is of Date type
BT_pred_prey_cols$death_date <- as.Date(BT_pred_prey_cols$death_date, origin = "1970-01-01")

# Filter out data after the predation or mortality event for each individual.
BT_telem_data_2 <- 
  BT_telem_data %>%
  left_join(BT_pred_prey_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(death_date)| Date < death_date)  # Keep locations only before the death date if one is recorded
#pre-filter rows: 20165124
#post-filter rows: 19936703

#how many rows removed

print(paste0("Rows removed after  filtering: ", nrow(BT_telem_data) - nrow(BT_telem_data_2)))
#228421

#check removed data
BT_removed_data <- 
  BT_telem_data %>%
  left_join(BT_pred_prey_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(!is.na(death_date) & Date >= death_date)

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------#
# 2. Filter out pike post-mortality data ####
#-------------------------------------------------#

#Some pike died in Lake BT. Need to filter out date post-suspected mortality date

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
pike_mort_cols <- 
  pike_deaths %>%
  filter(individual_ID %in% BT_telem_data$individual_ID) %>% 
  mutate(
    pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y"))

# Filter out data after the predation or mortality event for each individual.
BT_telem_data_3 <- 
  BT_telem_data_2 %>%
  left_join(pike_mort_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(pike_death_date)| Date < pike_death_date)  # Keep locations only before the death date if one is recorded
#pre-filter rows: 19936703
#post-filter rows: 18833395

print(paste0("Rows removed after  filtering: ", nrow(BT_telem_data_2) - nrow(BT_telem_data_3)))
#1103308

#save
saveRDS(BT_telem_data_3, paste0(filtered_data_path, "04_lake_BT_sub.rds"))

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#

#------------------------------------------------------------#
# 3. Rerun ctmm models for individuals that were predated ####
#------------------------------------------------------------#

#I need to pull the Reference tag data again 
#load dataframe that contains this information
get_ref_data <-  readRDS(paste0(filtered_data_path, "01_lake_BT_sub.rds"))

ref_data <- get_ref_data %>% 
  filter(individual_ID == 'FReference' | Species == 'Northern Pike') %>% 
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


ctmm::projection(ref_tel$FReference) <- ctmm::median(ref_tel)
BT_UERE <- uere.fit(ref_tel$FReference)
summary(BT_UERE)
# low      est      high
# all 0.37805 0.379299 0.3805479

#save BT UERE for future use
#saveRDS(UERE, paste0(filtered_data_path, 'BT_UERE.rds'))

#-------------------------------------------------------------------------------#

BT_filt_data <- readRDS(paste0(filtered_data_path, "04_lake_BT_sub.rds"))
#load UERE if needed
BT_UERE <- readRDS(paste0(filtered_data_path, "BT_UERE.rds"))


#> 2.1. Perch ####

perch_pred_ids <- 
  BT_pred_prey_cols %>%
  filter(species == 'Perch') %>% 
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated') %>% 
  pull(individual_ID)

#need to re-run ctmms for F59752, F59757 and F59792
#first double check differences between unfiltered and filtered dataset
# Count rows per ID in BT_filt_data
filt_counts <- BT_filt_data %>%
  filter(individual_ID %in% perch_pred_ids) %>%
  count(individual_ID, name = "BT_filt_data_rows")

# Count rows per ID in BT_telem_data
telem_counts <- BT_telem_data %>%
  filter(individual_ID %in% perch_pred_ids) %>%
  count(individual_ID, name = "BT_telem_data_rows")

# Combine results into one table
row_comparison <- full_join(filt_counts, telem_counts, by = "individual_ID")

# View result
print(row_comparison)

### Setup to re-run ctmms ###

pred_perch_data <- 
  BT_filt_data %>% 
  filter(individual_ID == perch_pred_ids)

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
uere(pred_perch_tel) <- BT_UERE

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(7)
doParallel::registerDoParallel(cl)

pred_perch_fits <- 
  foreach(i = 1:length(pred_perch_tel), .packages = 'ctmm') %dopar% {
    pred_perch_guess <- ctmm.guess(pred_perch_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(pred_perch_tel[[i]], pred_perch_guess, method = 'ML')
    # Save individual model fits to a specified folder.
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_BT_perch_fits/", names(pred_perch_tel)[i], ".rds"))
    model_fit
  }

saveRDS(pred_perch_fits, file = paste0(save_ctmm_path, "lake_BT_perch_fits/", "pred_perch_fits.rds"))

#-----------------------------------------------------------------------------------------------#

#> 2.2. Roach ####

roach_pred_ids <- 
  BT_pred_prey_cols %>%
  filter(species == 'Roach') %>% 
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated') %>% 
  pull(individual_ID)

#need to re-run ctmms for F59783, F59803, F59809, F59810
#need to also re-run ctmms for F59812, F59815, F59807, F59776

# Additional IDs to include
additional_ids <- c("F59812", "F59815", "F59807", "F59776")
# Combine and ensure uniqueness
roach_pred_ids <- union(roach_pred_ids, additional_ids)

#first double check differences between unfiltered and filtered dataset
# Count rows per ID in BT_filt_data
filt_counts <- BT_filt_data %>%
  filter(individual_ID %in% roach_pred_ids) %>%
  count(individual_ID, name = "BT_filt_data_rows")

# Count rows per ID in BT_telem_data
telem_counts <- BT_telem_data %>%
  filter(individual_ID %in% roach_pred_ids) %>%
  count(individual_ID, name = "BT_telem_data_rows")

# Combine results into one table
row_comparison <- full_join(filt_counts, telem_counts, by = "individual_ID")

print(row_comparison)


### Setup to re-run ctmms ###

pred_roach_data <- 
  BT_filt_data %>% 
  filter(individual_ID %in% roach_pred_ids)

pred_roach_data <- 
  with(pred_roach_data, 
       data.frame(
         "timestamp" = timestamp,                        
         "location.long" = longitude,                         
         "location.lat" = latitude, 
         "GPS.HDOP" = HDOP,                               
         "individual-local-identifier" = individual_ID))

#create telemetry object for predated roach
pred_roach_tel <- as.telemetry(pred_roach_data, 
                               timezone = "Europe/Stockholm",   
                               timeformat = "%Y-%m-%d %H:%M:%S",
                               projection = NULL,               
                               datum = "WGS84")

ctmm::projection(pred_roach_tel) <- ctmm::median(pred_roach_tel)
uere(pred_roach_tel) <- BT_UERE

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(8)
doParallel::registerDoParallel(cl)

pred_roach_fits <- 
  foreach(i = 1:length(pred_roach_tel), .packages = 'ctmm') %dopar% {
    pred_roach_guess <- ctmm.guess(pred_roach_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(pred_roach_tel[[i]], pred_roach_guess, method = 'ML')
    # Save individual model fits to a specified folder.
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_BT_roach_fits/", names(pred_roach_tel)[i], ".rds"))
    model_fit
  }

saveRDS(pred_roach_fits, file = paste0(save_ctmm_path, "lake_BT_roach_fits/", "pred_roach_fits.rds"))

#---------------------------------------------------------------------------------------------------------#

#> 2.3. Pike ####

pike_ids <- 
  pike_deaths %>%
  filter(individual_ID %in% BT_filt_data$individual_ID) %>% 
  mutate(
    pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y")) %>%  
  pull(individual_ID)


pike_mort_data <- 
  BT_filt_data %>% 
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

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(3)
doParallel::registerDoParallel(cl)

pike_mort_fits <- 
  foreach(i = 1:length(pike_mort_tel), .packages = 'ctmm') %dopar% {
    pike_mort_guess <- ctmm.guess(pike_mort_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(pike_mort_tel[[i]], pike_mort_guess, method = 'ML')
    # Save individual model fits to a specified folder.
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_BT_pike_fits/", names(pike_mort_tel)[i], ".rds"))
    model_fit
  }


#-----------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------#

#EXTRACT OUF MODELS
#CREATE CTMM LIST FOR BT ROACH

###
open_ctmm_path = "./data/ctmm_fits/lake_BT_pike_fits/"
save_ctmm_path = "./data/ctmm_fits/lake_BT_pike_fits/"

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

pike_pred_ids <- c("F59892", "F59886", "F59889")

lake_BT_pike_ctmm_fits <- rds_list[!names(rds_list) %in% pike_pred_ids]
print(names(lake_BT_pike_ctmm_fits))

# Create a list of remaining elements
remaining_list <- rds_list[names(rds_list) %in% pike_pred_ids]
print(names(remaining_list))

# ------------------ Extract OUF models ----------------------- #

lake_BT_pike_OUF_models <- list()

#Iterate over each object in muddyfoot_roach_ctmm_fits and extract the 'OUF anisotropic error' models
lake_BT_pike_OUF_models <- lapply(lake_BT_pike_ctmm_fits, function(x) x[['OUF anisotropic error']])
#check it worked
summary(lake_BT_pike_OUF_models$F59887)

#rejoin pred filtered models
lake_BT_pike_OUF_models <- c(lake_BT_pike_OUF_models, remaining_list)[names(rds_list)]
#check that ID are in chonological order
names(lake_BT_pike_OUF_models)

#save
saveRDS(lake_BT_pike_OUF_models, paste0(save_ctmm_path, "lake_BT_pike_OUF_models.rds"))

#----------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# 4. Rerun telemetry object for all individuals with filtered dataset ####
#------------------------------------------------------------------------#

BT_filt_data <- readRDS(paste0(filtered_data_path, "04_lake_BT_sub.rds"))

BT_movebank <- 
  with(BT_filt_data, 
       data.frame(
         "timestamp" = timestamp,                        
         "location.long" = longitude,                         
         "location.lat" = latitude, 
         "GPS.HDOP" = HDOP,                               
         "individual-local-identifier" = individual_ID, 
         "Species" = Species,
         "Treatment" = Treatment,                        
         "Date" = Date,
         "Exp_Stage" = Exp_Stage,
         "Time_Of_Day" = Time_Of_Day
       ))

BT_tel <- as.telemetry(BT_movebank, 
                        timezone = "Europe/Stockholm", 
                        timeformat="%Y-%m-%d %H:%M:%S", 
                        projection= NULL,
                        datum="WGS84",
                        keep = c("Species","Treatment", 
                                 "Date", "Exp_Stage", "Time_Of_Day"))

ctmm::projection(BT_tel) <- ctmm::median(BT_tel)

# Remove outliers based on species maximum swim speeds #

#Pike = 0.823
#Perch = 0.977
#Roach = 0.841

#first seperate telemetry object by species 
#check whether ids are in species order
BT_movebank %>%
  select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

#need to order them by id number
BT_tel <- BT_tel[order(names(BT_tel))]
names(BT_tel)

#mostly in order except for the first perch
pike_BT_tel <- BT_tel[61:66]
perch_BT_tel <- BT_tel[c(1:15, 31:44, 46)]
roach_BT_tel <- BT_tel[c(16:30, 45, 47:60)]

#remove outliers based on speed
#pike
out_pike <- outlie(pike_BT_tel, plot = FALSE)
sum(sapply(out_pike, function(x) sum(x$speed > 0.823)))

saveRDS(pike_BT_tel , paste0(telem_path, "BT/pike_BT_tel.rds")) 

#Perch
out_perch <- outlie(perch_BT_tel, plot = FALSE)
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))
#11

saveRDS(perch_BT_tel , paste0(telem_path, "BT/perch_BT_tel.rds"))

#Roach
out_roach <- outlie(roach_BT_tel, plot = FALSE)
sum(sapply(out_roach, function(x) sum(x$speed > 0.841)))

#save telemetry object
saveRDS(roach_BT_tel , paste0(telem_path, "BT/roach_BT_tel.rds"))



