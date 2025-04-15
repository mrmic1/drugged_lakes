#--------------------------------------------------------------------------#
# Filter data to remove post-predation and mortality tracking - Muddyfoot  #
#--------------------------------------------------------------------------#

#SCRIPT DESCRIPTION

#In this script we will remove data where individuals were tracked after being predated or had 
#an identified mortality event
#Predation events are identified in the `04_muddyfoot_calculate_pred_prey_ints` script
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
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
telem_path <- "./data/telem_obj/muddyfoot/"
enc_path <- "./data/encounters/muddyfoot/"
save_ctmm_path = "./data/ctmm_fits/"

### DATA LOADING ###
muddyfoot_telem_data <-  readRDS(paste0(filtered_data_path, "03_muddyfoot_sub.rds"))
# Load the predation event dataframe 
mortality_preds <- readxl::read_excel("./data/encounters/suspected_mortality_updated.xlsx")

#------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------#
# 1. Filter out post-predation and mortality tracking locations #####
#-------------------------------------------------------------------#

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
muddyfoot_pred_prey_cols <-
  mortality_preds %>%
  filter(lake == 'muddyfoot') %>% 
  filter(species == 'Roach'| species == 'Perch') %>%
  dplyr::select(individual_ID, species, revised_suspected_mortality, revised_likely_death_date) %>% 
  rename(death_date = revised_likely_death_date)

# Ensure the merged_date column is of Date type
muddyfoot_pred_prey_cols$death_date <- as.Date(muddyfoot_pred_prey_cols$death_date, origin = "1970-01-01")

# Filter out data after the predation or mortality event for each individual.
muddyfoot_telem_data_2 <- 
  muddyfoot_telem_data %>%
  left_join(muddyfoot_pred_prey_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(death_date)| Date < death_date)  # Keep locations only before the death date if one is recorded
#pre-filter rows: 10358641
#post-filter rows: 9946860

#how many rows removed
print(paste0("Rows removed after  filtering: ", nrow(muddyfoot_telem_data) - nrow(muddyfoot_telem_data_2)))
#411781

#check removed data
muddyfoot_removed_data <- 
  muddyfoot_telem_data %>%
  left_join(muddyfoot_pred_prey_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(!is.na(death_date) & Date >= death_date)


#save
saveRDS(muddyfoot_telem_data_2, paste0(filtered_data_path, "04_muddyfoot_sub.rds"))

#------------------------------------------------------------------------------------------------------------#
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
muddyfoot_UERE <- uere.fit(ref_tel)
summary(muddyfoot_UERE)

#save muddyfoot UERE for future use
saveRDS(muddyfoot_UERE, paste0(filtered_data_path, 'muddyfoot_UERE.rds'))

#-------------------------------------------------------------------------------#

muddyfoot_filt_data <- readRDS(paste0(filtered_data_path, "04_muddyfoot_sub.rds"))
#load UERE if needed
muddyfoot_UERE <- readRDS(paste0(filtered_data_path, "muddyfoot_UERE.rds"))


#> 2.1. Perch ####

perch_pred_ids <- 
  muddyfoot_pred_prey_cols %>%
  filter(species == 'Perch') %>% 
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated') %>% 
  pull(individual_ID)

#need to re-run ctmms for F59752, F59757 and F59792
#first double check differences between unfiltered and filtered dataset
# Count rows per ID in muddyfoot_filt_data
filt_counts <- muddyfoot_filt_data %>%
  filter(individual_ID %in% perch_pred_ids) %>%
  count(individual_ID, name = "muddyfoot_filt_data_rows")

# Count rows per ID in muddyfoot_telem_data
telem_counts <- muddyfoot_telem_data %>%
  filter(individual_ID %in% perch_pred_ids) %>%
  count(individual_ID, name = "muddyfoot_telem_data_rows")

# Combine results into one table
row_comparison <- full_join(filt_counts, telem_counts, by = "individual_ID")

# View result
print(row_comparison)

### Setup to re-run ctmms ###

pred_perch_data <- 
  muddyfoot_filt_data %>% 
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
uere(pred_perch_tel) <- muddyfoot_UERE

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(2)
doParallel::registerDoParallel(cl)

pred_perch_fits <- 
  foreach(i = 1:length(pred_perch_tel), .packages = 'ctmm') %dopar% {
    pred_perch_guess <- ctmm.guess(pred_perch_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(pred_perch_tel[[i]], pred_perch_guess, method = 'ML')
    # Save individual model fits to a specified folder.
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_perch_fits/", names(pred_perch_tel)[i], ".rds"))
    model_fit
  }

saveRDS(pred_perch_fits, file = paste0(save_ctmm_path, "muddyfoot_perch_fits/", "pred_perch_fits.rds"))


#> 2.2. Roach ####

roach_pred_ids <- 
  mud_pred_mort_events %>%
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated'| n_missing_dates > 24) %>% 
  mutate(
    first_date_over_50 = as.Date(first_date_over_50),
    death_date = as.Date(death_date)) %>% 
  dplyr::select(individual_ID, Species, first_date_over_50, death_date) %>% 
  filter(Species == 'Roach') %>% 
  pull(individual_ID)

pred_roach_data <- 
  muddyfoot_filt_data %>% 
  filter(individual_ID == roach_pred_ids)

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
uere(pred_roach_tel) <- UERE

#Initialize a parallel cluster to speed up model fitting for each fish individual.
cl <- makeCluster(2)
doParallel::registerDoParallel(cl)

pred_roach_fits <- 
  foreach(i = 1:length(pred_roach_tel), .packages = 'ctmm') %dopar% {
    pred_roach_guess <- ctmm.guess(pred_roach_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(pred_roach_tel[[i]], pred_roach_guess, method = 'ML')
    # Save individual model fits to a specified folder.
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_roach_fits/", names(pred_roach_tel)[i], ".rds"))
    model_fit
  }

saveRDS(pred_roach_fits, file = paste0(save_ctmm_path, "muddyfoot_roach_fits/", "pred_roach_fits.rds"))





#----------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# 3. Rerun telemetry object for all individuals with filtered dataset ####
#------------------------------------------------------------------------#

muddyfoot_filt_data <- readRDS(paste0(filtered_data_path, "04_muddyfoot_sub.rds"))

mud_movebank <- 
  with(muddyfoot_filt_data, 
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

mud_tel <- as.telemetry(mud_movebank, 
                                   timezone = "Europe/Stockholm", 
                                   timeformat="%Y-%m-%d %H:%M:%S", 
                                   projection= NULL,
                                   datum="WGS84",
                                   keep = c("Species","Treatment", 
                                            "Date", "Exp_Stage", "Time_Of_Day"))


ctmm::projection(mud_tel) <- ctmm::median(mud_tel)

# Remove outliers based on species maximum swim speeds #

#Pike = 0.823
#Perch = 0.977
#Roach = 0.841

#first seperate telemetry object by species 
#check whether ids are in species order
mud_movebank %>%
  dplyr::select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

#need to order them by id number
mud_tel <- mud_tel[order(names(mud_tel))]
names(mud_tel)

#mostly in order except for the first perch
pike_mud_tel <- mud_tel[60:65]
perch_mud_tel <- mud_tel[c(1,4,5,7,8,10:12,15:17,20,25,28,29,31:33,35:37,39,40,41,43, 46:50)]
#For roach do not include F59707 - died at start of study
roach_mud_tel <- mud_tel[c(2,3,6,9,13,14,18,19,21:23,26,27,30,34,38,42,44,45,51:59)]

#remove outliers based on speed
#pike
out_pike <- outlie(pike_mud_tel, plot = FALSE)
sum(sapply(out_pike, function(x) sum(x$speed > 0.823)))
#14

saveRDS(pike_mud_tel , paste0(telem_path, "muddyfoot/pike_muddyfoot_tel.rds")) 

#Perch
out_perch <- outlie(perch_mud_tel, plot = FALSE)
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))
#11

saveRDS(perch_mud_tel , paste0(telem_path, "muddyfoot/perch_muddyfoot_tel.rds"))

#Roach
out_roach <- outlie(roach_mud_tel, plot = FALSE)
sum(sapply(out_roach, function(x) sum(x$speed > 0.841)))

saveRDS(roach_mud_tel , paste0(telem_path, "muddyfoot/roach_muddyfoot_tel.rds"))



