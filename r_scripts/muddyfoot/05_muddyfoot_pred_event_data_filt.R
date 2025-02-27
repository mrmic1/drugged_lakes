#--------------------------------------------------------------#
# Filter data to remove post-predation and mortality tracking  #
#--------------------------------------------------------------#

#In this script we will remove data where individuals were tracked after being predated or had 
#an identified mortality event
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
enc_path <- "./data/encounters/muddyfoot/"
save_ctmm_path = "./data/ctmm_fits/"

### DATA LOADING ###
muddyfoot_telem_data <-  readRDS(paste0(filtered_data_path, "03_muddyfoot_sub.rds"))
# Load the predation event dataframe 
mud_pred_mort_events <- readRDS(paste0(enc_path, "muddyfoot_pred_encounter_summary_filtered.rds"))

#-------------------------------------------------#
# 1. Filter out post-predation event data ####
#-------------------------------------------------#

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
pred_cols <- 
  mud_pred_mort_events %>%
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated'| n_missing_dates > 24) %>% 
  mutate(
    first_date_over_50 = as.Date(first_date_over_50),
    death_date = as.Date(death_date)) %>% 
  dplyr::select(individual_ID, first_date_over_50, death_date)

# Merge into single death_date column
pred_cols$death_date <- ifelse(!is.na(pred_cols$first_date_over_50), pred_cols$first_date_over_50, pred_cols$death_date)

# Ensure the merged_date column is of Date type
pred_cols$death_date <- as.Date(pred_cols$death_date, origin = "1970-01-01")

# Filter out data after the predation or mortality event for each individual.
muddyfoot_filt_data <- 
  muddyfoot_telem_data %>%
  left_join(pred_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(death_date)| Date <= death_date)  # Keep only pre-predation data
#pre-filter rows: 10358641
#post-filter rows: 10032122
#326519 rows removed

print(paste0("Rows removed after  filtering: ", nrow(muddyfoot_telem_data) - nrow(muddyfoot_filt_data)))

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

muddyfoot_filt_data <- readRDS(paste0(filtered_data_path, "04_muddyfoot_sub.rds"))

perch_pred_ids <- 
  mud_pred_mort_events %>%
  filter(revised_suspected_mortality == 'mortality' | revised_suspected_mortality == 'likely_predated'| n_missing_dates > 24) %>% 
  mutate(
    first_date_over_50 = as.Date(first_date_over_50),
    death_date = as.Date(death_date)) %>% 
  dplyr::select(individual_ID, Species, first_date_over_50, death_date) %>% 
  filter(Species == 'Perch') %>% 
  pull(individual_ID)

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
uere(pred_perch_tel) <- UERE

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


ctmm::projection(mud_tel_tpeqd$F59747)
tz(mud_tel_tpeqd$F59747$timestamp)
ctmm::projection(mud_tel_tpeqd) <- ctmm::median(mud_tel_tpeqd)

# Remove outliers based on species maximum swim speeds #

#Pike = 0.823
#Perch = 0.977
#Roach = 0.841

#first seperate telemetry object by species 
#check whether ids are in species order
mud_movebank %>%
  select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

#need to order them by id number
mud_tel <- mud_tel[order(names(mud_tel))]
names(mud_tel)

#mostly in order except for the first perch
pike_mud_tel <- mud_tel[60:65]
perch_mud_tel <- mud_tel[c(1,4,5,7,8,10:12,15:17,20,25,28,29,31:33,35:37,40,41,43, 46:50)]
roach_mud_tel <- mud_tel[c(2,3,6,9,13,14,18,19,21:24,26,27,30,34,38,42,44,45,51:59)]

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
#

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_roach, function(x) x$speed <= 0.841)
#Combining the lists and removing observations for which the logical vector was false
roach_mud_tel <- Map(function(x,y) x[y,], roach_mud_tel,which_lowSp)
#save telemetry object
saveRDS(roach_mud_tel , paste0(telem_path, "muddyfoot/roach_muddyfoot_tel.rds"))



