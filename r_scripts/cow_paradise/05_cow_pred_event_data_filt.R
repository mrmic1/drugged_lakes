#-----------------------------------------------------------------------------#
# Filter data to remove post-predation and mortality tracking - Cow Paradise  #
#-----------------------------------------------------------------------------#

#SCRIPT DESCRIPTION

#In this script we will remove data where individuals were tracked after being predated or had 
#an identified mortality event
#Predation events are identified in the `04_cow_calculate_pred_prey_ints` script
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
filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"
telem_path <- "./data/telem_obj/cow_paradise/"
enc_path <- "./data/encounters/cow_paradise/"
save_ctmm_path = "./data/ctmm_fits/"

### DATA LOADING ###
cow_telem_data <-  readRDS(paste0(filtered_data_path, "03_lake_cow_sub.rds"))
# Load the predation event dataframe 
mortality_preds <- readxl::read_excel("./data/encounters/suspected_mortality_updated.xlsx")
# Load suspected dates for pike deaths
pike_deaths <- read.csv("./data/encounters/pike_deaths.csv")

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------#
# 1. Filter out post-predation and mortality tracking locations #####
#-------------------------------------------------------------------#

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
cow_pred_prey_cols <-
  mortality_preds %>%
  filter(lake == 'cow paradise') %>% 
  filter(species == 'Roach'| species == 'Perch') %>%
  dplyr::select(individual_ID, species, revised_suspected_mortality, revised_likely_death_date) %>% 
  rename(death_date = revised_likely_death_date)

# Ensure the merged_date column is of Date type
cow_pred_prey_cols$death_date <- as.Date(cow_pred_prey_cols$death_date, origin = "1970-01-01")

# Filter out data after the predation or mortality event for each individual.
cow_telem_data_2 <- 
  cow_telem_data %>%
  left_join(cow_pred_prey_cols, by = c("individual_ID" = "individual_ID")) %>%
  # Keep locations only before the death date if one is recorded
  filter(is.na(death_date)| Date < death_date) %>% 
  #Remove individuals with poor tracking
  filter(revised_suspected_mortality != "poor_tracking_remove" | is.na(revised_suspected_mortality))

#pre-filter rows: 8860993
#post-filter rows: 8686231

#how many rows removed

print(paste0("Rows removed after  filtering: ", nrow(cow_telem_data) - nrow(cow_telem_data_2)))
#174762

# #check removed data
# cow_removed_data <- 
#   cow_telem_data %>%
#   left_join(cow_pred_prey_cols, by = c("individual_ID" = "individual_ID")) %>%
#   filter(!is.na(death_date) & Date >= death_date)

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------#
# 2. Filter out pike post-mortality data ####
#-------------------------------------------------#

#Some pike died in Lake cow. Need to filter out date post-suspected mortality date

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
pike_mort_cols <- 
  pike_deaths %>%
  filter(individual_ID %in% cow_telem_data$individual_ID) %>% 
  mutate(
    pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y"))

# Filter out data after the predation or mortality event for each individual.
cow_telem_data_3 <- 
  cow_telem_data_2 %>%
  left_join(pike_mort_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(pike_death_date)| Date < pike_death_date)  # Keep locations only before the death date if one is recorded
#pre-filter rows: 8686231
#post-filter rows: 8649347

print(paste0("Rows removed after  filtering: ", nrow(cow_telem_data_2) - nrow(cow_telem_data_3)))
#36884

#save
saveRDS(cow_telem_data_3, paste0(filtered_data_path, "04_lake_cow_sub.rds"))

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#------------------------------------------------------------#
# 3. Rerun ctmm models for individuals that were predated ####
#------------------------------------------------------------#

#I need to pull the Reference tag data again 
#load dataframe that contains this information
get_ref_data <-  readRDS(paste0(filtered_data_path, "01_lake_cow_sub.rds"))

ref_data <- get_ref_data %>% 
  filter(individual_ID == 'FReference'| Species == 'Northern Pike') %>% 
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
cow_UERE <- uere.fit(ref_tel$FReference)
summary(cow_UERE)
# low       est      high
# all 0.8200953 0.8265717 0.8330474

#save BT UERE for future use
saveRDS(cow_UERE, paste0(filtered_data_path, 'cow_UERE.rds'))
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

cow_filt_data <- readRDS(paste0(filtered_data_path, "04_lake_cow_sub.rds"))
#load UERE if needed
cow_UERE <- readRDS(paste0(filtered_data_path, 'cow_UERE.rds'))

#> 2.1. Perch ####

perch_pred_ids <- 
  cow_pred_prey_cols %>%
  filter(species == 'Perch') %>% 
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated') %>% 
  pull(individual_ID)

#first double check differences between unfiltered and filtered dataset
# Count rows per ID in cow_filt_data
filt_counts <- cow_filt_data %>%
  filter(individual_ID %in% perch_pred_ids) %>%
  count(individual_ID, name = "cow_filt_data_rows")

# Count rows per ID in cow_telem_data
telem_counts <- cow_telem_data %>%
  filter(individual_ID %in% perch_pred_ids) %>%
  count(individual_ID, name = "cow_telem_data_rows")

# Combine results into one table
row_comparison <- full_join(filt_counts, telem_counts, by = "individual_ID")

# View result
print(row_comparison)

### Setup to re-run ctmms ###

pred_perch_data <- 
  cow_filt_data %>% 
  filter(individual_ID %in% perch_pred_ids)

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
uere(pred_perch_tel) <- cow_UERE

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(5)
doParallel::registerDoParallel(cl)

pred_perch_fits <- 
  foreach(i = 1:length(pred_perch_tel), .packages = 'ctmm') %dopar% {
    pred_perch_guess <- ctmm.guess(pred_perch_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(pred_perch_tel[[i]], pred_perch_guess, method = 'ML')
    # Save individual model fits to a specified folder.
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_cow_perch_fits/", names(pred_perch_tel)[i], ".rds"))
    model_fit
  }

saveRDS(pred_perch_fits, file = paste0(save_ctmm_path, "lake_cow_perch_fits/", "pred_perch_fits.rds"))

#F59849

F59849_tel <-pred_perch_tel$F59849

F59849_guess <- ctmm.guess(F59849_tel, CTMM=ctmm(error=TRUE), interactive = FALSE)
F59849_fit <- ctmm.fit(F59849_tel, F59849_guess, method = 'ML')
saveRDS(F59849_fit, file = paste0(save_ctmm_path, "lake_cow_perch_fits/", "F59849.rds"))



#-----------------------------------------------------------------------------------------------#
#> 2.2. Roach ####

roach_pred_ids <- 
  cow_pred_prey_cols %>%
  filter(species == 'Roach') %>% 
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated'| revised_suspected_mortality == 'known_predated') %>% 
  pull(individual_ID)


#first double check differences between unfiltered and filtered dataset
# Count rows per ID in cow_filt_data
filt_counts <- cow_filt_data %>%
  filter(individual_ID %in% roach_pred_ids) %>%
  count(individual_ID, name = "cow_filt_data_rows")

# Count rows per ID in cow_telem_data
telem_counts <- cow_telem_data %>%
  filter(individual_ID %in% roach_pred_ids) %>%
  count(individual_ID, name = "cow_telem_data_rows")

# Combine results into one table
row_comparison <- full_join(filt_counts, telem_counts, by = "individual_ID")

print(row_comparison)


### Setup to re-run ctmms ###

pred_roach_data <- 
  cow_filt_data %>% 
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
uere(pred_roach_tel) <- cow_UERE

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(3)
doParallel::registerDoParallel(cl)

pred_roach_fits <- 
  foreach(i = 1:length(pred_roach_tel), .packages = 'ctmm') %dopar% {
    pred_roach_guess <- ctmm.guess(pred_roach_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(pred_roach_tel[[i]], pred_roach_guess, method = 'ML')
    # Save individual model fits to a specified folder.
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_cow_roach_fits/", names(pred_roach_tel)[i], ".rds"))
    model_fit
  }

saveRDS(pred_roach_fits, file = paste0(save_ctmm_path, "lake_cow_roach_fits/", "pred_roach_fits.rds"))

#-------------------------------------------------------------------------------------------------------------------------#

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
uere(pike_mort_tel) <- cow_UERE

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

roach_pred_ids <- c("F59826", "F59828", "F59817")

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

summary(lake_cow_roach_OUF_models$F59817)

#save
saveRDS(lake_cow_roach_OUF_models, paste0(save_ctmm_path, "lake_cow_roach_OUF_models.rds"))

#----------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# 4. Rerun telemetry object for all individuals with filtered dataset ####
#------------------------------------------------------------------------#

cow_filt_data <- readRDS(paste0(filtered_data_path, "04_lake_cow_sub.rds"))

cow_movebank <- 
  with(cow_filt_data, 
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

cow_tel <- as.telemetry(cow_movebank, 
                       timezone = "Europe/Stockholm", 
                       timeformat="%Y-%m-%d %H:%M:%S", 
                       projection= NULL,
                       datum="WGS84",
                       keep = c("Species","Treatment", 
                                "Date", "Exp_Stage", "Time_Of_Day"))

ctmm::projection(cow_tel) <- ctmm::median(cow_tel)

# Remove outliers based on species maximum swim speeds #

#Pike = 0.823
#Perch = 0.977
#Roach = 0.841

#first seperate telemetry object by species 
#check whether ids are in species order
cow_movebank %>%
  select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

#need to order them by id number
cow_tel <- cow_tel[order(names(cow_tel))]
names(cow_tel)

#mostly in order except for the first perch
pike_cow_tel <- cow_tel[61:65]
perch_cow_tel <- cow_tel[c(1, 22:60)]
roach_cow_tel <- cow_tel[c(2:21)]

#remove outliers based on speed
#pike
out_pike <- outlie(pike_cow_tel, plot = FALSE)
sum(sapply(out_pike, function(x) sum(x$speed > 0.823)))

saveRDS(pike_cow_tel , paste0(telem_path, "cow_paradise/pike_cow_tel.rds")) 

#Perch
out_perch <- outlie(perch_cow_tel, plot = FALSE)
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))

saveRDS(perch_cow_tel , paste0(telem_path, "cow_paradise/perch_cow_tel.rds"))

#Roach
out_roach <- outlie(roach_cow_tel, plot = FALSE)
sum(sapply(out_roach, function(x) sum(x$speed > 0.841)))

#save telemetry object
saveRDS(roach_cow_tel , paste0(telem_path, "cow_paradise/roach_cow_tel.rds"))
