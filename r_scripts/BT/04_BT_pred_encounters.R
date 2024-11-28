# -----------------------------------------------# 
# Exploring predator-prey interactions - BT
# ------------------------------------------------#

# This script analyses telemetry data of prey (Roach, Perch) and predators (Pike)
# in BT, including pairwise distances and encounter rates.
# The goal is to identify potential predation events based on close proximity and to create 
# a dataframe that can be used for analysing potential treatment effects on predator encounters

### LIBRARIES ###
# Loading required packages
library(dplyr)
library(ctmm)
library(lubridate)
library(data.table)
library(sf)  # For spatial data handling
library(flextable)

# Set the time zone environment to 'Europe/Stockholm' for consistent timestamp manipulation
Sys.setenv(TZ = 'Europe/Stockholm')

### DIRECTORIES ###
# Paths to directories containing the necessary data files
filtered_data_path <- "./data/tracks_filtered/lake_BT/"
lake_polygon_path <- "./data/lake_coords/"
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/BT/"
enc_path <- "./data/encounters/BT/"
size_path <- "./data/fish_size/"
save_tables_path <- "./tables/BT/"  # Path to save summary tables.

# SET TABLE PLOTTING PARAMETERS
# Customize default settings for Flextable tables used later in the script.
set_flextable_defaults(
  font.color = "black",
  border.color = "black",
  font.family = 'Arial',
  line_spacing = 1
)


### DATA LOADING ###
BT_filt_data <-  readRDS(paste0(filtered_data_path, "03_lake_BT_sub.rds"))

pike_BT_tel <- readRDS(paste0(telem_path, 'pike_BT_tel.rds'))
perch_BT_tel <- readRDS(paste0(telem_path, 'perch_BT_tel.rds'))
roach_BT_tel <- readRDS(paste0(telem_path, 'roach_BT_tel.rds'))

# Load lake polygon data
BT_polygon <- sf::st_read(paste0(lake_polygon_path, "lake_BT_polygon.gpkg"))


#---------------------------------------------------------------------------------#
# 1. Calculate distances at shared timepoints for prey and predators #######
#---------------------------------------------------------------------------------#

# This section calculates pairwise distances between Roach/Perch (prey) and Pike (predator) at shared time points.

pike_BT_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_BT_pike_fits/lake_BT_pike_OUF_models.rds"))
perch_BT_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_BT_perch_fits/lake_BT_perch_OUF_models.rds"))
roach_BT_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_BT_roach_fits/lake_BT_roach_OUF_models.rds"))

#### > 1.1. Roach ####

# Reproject Roach telemetry objects to match the Pike projection
ctmm::projection(roach_BT_tel) <- ctmm::projection(pike_BT_tel)

#make sure pike ctmms are on the same trajectory
ctmm_pike <- pike_BT_ctmm_fits$F59886
ctmm::projection(pike_BT_ctmm_fits) <- ctmm::projection(ctmm_pike)
#now for roach
ctmm::projection(roach_BT_ctmm_fits) <- ctmm::projection(pike_BT_ctmm_fits)

#ctmm's with same projections
saveRDS(roach_BT_ctmm_fits, paste0(ctmm_path, "lake_BT_roach_fits/lake_BT_roach_OUF_models.rds"))
saveRDS(pike_BT_ctmm_fits, paste0(ctmm_path, "lake_BT_pike_fits/lake_BT_pike_OUF_models.rds"))

# Check projections to ensure they match
stopifnot(projection(roach_BT_tel) == projection(pike_BT_tel))

# Initialize an empty list to store distance calculations between Roach and Pike
roach_pike_distances <- list()

# Nested loop to calculate distances between each Roach and Pike telemetry pair
for(i in 1:length(roach_BT_tel)) {
  for(j in 1:length(pike_BT_tel)) {
    combined_telemetry <- c(roach_BT_tel[i], pike_BT_tel[j])
    combined_ctmm <- c(roach_BT_ctmm_fits[i], pike_BT_ctmm_fits[j])
    
    # Calculate pairwise distances
    location_difference <- distances(combined_telemetry, combined_ctmm)
    
    # Extract and store results with species IDs
    roach_id <- names(roach_BT_tel)[i]
    pike_id <- names(pike_BT_tel)[j]
    
    roach_pike_distances_df <- as.data.frame(location_difference)
    roach_pike_distances_df$Roach_ID <- roach_id
    roach_pike_distances_df$Pike_ID <- pike_id
    
    roach_pike_distances[[paste(roach_id, pike_id, sep = "_")]] <- roach_pike_distances_df
  }
}

# Combine all calculated distances into a single dataframe
roach_pike_distances_df <- do.call(rbind, roach_pike_distances)

# Convert timestamps to a consistent timezone if present
if ("timestamp" %in% colnames(roach_pike_distances_df)) {
  roach_pike_distances_df$timestamp <- as.POSIXct(roach_pike_distances_df$timestamp, tz = "Europe/Stockholm")
}

# View the first few rows of the distance data
head(roach_pike_distances_df)

# Save the distance data for future analysis
saveRDS(roach_pike_distances_df, paste0(enc_path, "BT/BT_pike_roach_distances_df.rds"))


#### > 1.2. Perch ####

# Reproject Perch telemetry objects to match Pike
projection(perch_BT_tel) <- projection(pike_BT_tel)
projection(perch_BT_ctmm_fits) <- projection(pike_BT_ctmm_fits)

saveRDS(perch_BT_ctmm_fits, paste0(ctmm_path, "lake_BT_perch_fits/lake_BT_perch_OUF_models.rds"))

# Ensure projections are correctly matched
stopifnot(projection(perch_BT_tel) == projection(pike_BT_tel))

# Initialize an empty list to store distance calculations between Perch and Pike
perch_pike_distances <- list()

# Nested loop for calculating Perch - Pike distances
for(i in 1:length(perch_BT_tel)) {
  for(j in 1:length(pike_BT_tel)) {
    combined_telemetry <- c(perch_BT_tel[i], pike_BT_tel[j])
    combined_ctmm <- c(perch_BT_ctmm_fits[i], pike_BT_ctmm_fits[j])
    
    location_difference <- distances(combined_telemetry, combined_ctmm)
    
    perch_id <- names(perch_BT_tel)[i]
    pike_id <- names(pike_BT_tel)[j]
    
    perch_pike_distances_df <- as.data.frame(location_difference)
    perch_pike_distances_df$Perch_ID <- perch_id
    perch_pike_distances_df$Pike_ID <- pike_id
    
    perch_pike_distances[[paste(perch_id, pike_id, sep = "_")]] <- perch_pike_distances_df
  }
}

# Combine Perch-Pike distance data into one dataframe
perch_pike_distances_df <- do.call(rbind, perch_pike_distances)

# Convert timestamps to consistent timezone
if ("timestamp" %in% colnames(perch_pike_distances_df)) {
  perch_pike_distances_df$timestamp <- as.POSIXct(perch_pike_distances_df$timestamp, tz = "Europe/Stockholm")
}

# View the first few rows
head(perch_pike_distances_df)

# Save the distance data for future analysis
# saveRDS(perch_pike_distances_df, paste0(enc_path, "BT_pike_perch_distances_df.rds"))


#---------------------------------------------------#
# 2. Calculate encounters #########
#--------------------------------------------------#

roach_pike_distances_df <- readRDS(paste0(enc_path, "BT_pike_roach_distances_df.rds"))
perch_pike_distances_df <- readRDS(paste0(enc_path, "BT_pike_perch_distances_df.rds"))

#### > 2.1. Roach - Pike encounters #### 

# Define encounters as those occurring within 0.45 meters (based on Pike's strike distance)
roach_pike_distances_df$encounter <- ifelse(roach_pike_distances_df$est <= 0.45, 1, 0)

# Add date column based on timestamp
# Note: this takes a while to run
roach_pike_distances_df <- 
  roach_pike_distances_df %>%
  mutate(Date = as.Date(format(roach_pike_distances_df$timestamp, tz = "Europe/Stockholm")))

#save with encounter and date column
saveRDS(roach_pike_distances_df, paste0(enc_path, "BT_pike_roach_distances_df.rds"))

# Summarise encounter data per day and calculate average daily distance for each
# unique roach and pike combination
roach_pike_encounter_sum <- 
  roach_pike_distances_df %>%
  group_by(Roach_ID, Pike_ID, Date) %>%
  summarise(encounter_count = sum(encounter),
            daily_avg_dist = mean(est, na.rm = TRUE))


#### ____ 2.1.1. Identify possible predation events ####

# This section identifies potential predation events by filtering for days when pike and roach
# had more than 50 encounters or when their average daily distance was less than 3 meters.
# First we load post-experiment biometric data to assess known deaths or survivals
# and add this information to our encounter summary

post_biometrics <- 
  fread(paste0(size_path, "biometric_post_exp_data.csv")) %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number)))

#Summarise number of survivors and known death for each species in BT
post_biometrics %>% 
  filter(Lake == 'BT') %>% 
  dplyr::select(individual_ID, Species, Treatment, Found) %>% 
  group_by(Species, Treatment) %>% 
  summarise(num_survive = sum(Found),
            num_died = length(unique(individual_ID)) - sum(Found),
            percent_survive = sum(Found)/length(unique(individual_ID)))


# Merge biometric data (e.g., whether the fish was found alive) with encounter summary
roach_post_size_cols <- 
  post_biometrics %>%
  filter(Lake == 'BT', Species == 'Roach') %>%
  dplyr::select(individual_ID, Found, Known_predated) %>% 
  rename(found_alive = Found)

roach_pike_encounter_sum <- 
  roach_pike_encounter_sum %>%
  left_join(roach_post_size_cols, by = c("Roach_ID" = "individual_ID"))
  

# Now we filter roach-pike encounters for days with more than 50 encounters or an avg daily distance < 3 meters.
# This our threshold for identifying potential predation events.
# This threshold was arbitrarily chosen, but seems conservative 

roach_over_50_encounters <- 
  roach_pike_encounter_sum %>%
  filter(encounter_count >= 50  | daily_avg_dist < 3 | Known_predated == 1)

#NOTE: we know that one individual F59774 was definitely predated as it was found in the stomach of a roach
#however interestingly it never has a day with over 50 encounters with a pike
#it is included in above calculation but ignore its encounter over 50 count. 

#how many unique roach individuals
unique(roach_over_50_encounters$Roach_ID)
#23 individuals

# Summarise roach-pike interactions: count days with > 50 encounters, consecutive days, and the first date
BT_roach_pike_pred_int <- 
  roach_over_50_encounters %>%
  group_by(Roach_ID, Pike_ID) %>%
  summarise(
    num_days_over_50 = n(),  # Count of days with over 50 encounters or daily avg distance < 3 meters
    encounter_dates = paste(format(Date, "%d/%m/%Y"), collapse = ", "),  # Combine encounter dates as a string
    consecutive_days = max(rle(cumsum(c(1, diff(Date) > 1)) == 1)$lengths),  # Calculate consecutive encounter days
    first_date_over_50 = min(Date)  # Identify first date with more than 50 encounters
  )

# Extract individual treatment and missing date information from filtered data
roach_sel_cols <- 
  BT_filt_data %>%
  filter(Species == 'Roach') %>%
  dplyr::select(individual_ID, Treatment, n_missing_dates) %>%
  distinct()

# Merge encounter summary with individual-level metadata
BT_roach_pike_pred_int <- 
  BT_roach_pike_pred_int %>%
  left_join(roach_sel_cols, by = c("Roach_ID" = "individual_ID"))

# Merge biometric data (e.g., whether the fish was found alive) with possible predation events
# Create column indicated instances where there might be a higher likelihood that an individual was predated
# specifically, when individuals 2 or more consecutive days where they had more than 50 encounters with a predator
BT_roach_pike_pred_int <- 
  BT_roach_pike_pred_int %>%
  left_join(roach_post_size_cols, by = c("Roach_ID" = "individual_ID")) %>%
  # Remove fish found alive at the end of the experiment
  filter(found_alive == '0') %>%
  # Flag likely predation events (i.e., when consecutive encounters are more than 2 days)
  mutate(likely_predated = ifelse(consecutive_days >= 2, 1, 0))


# Identify the date with the most encounters for each roach and pike pair
BT_roach_max_encounter_dates <- 
  roach_pike_encounter_sum %>%
  group_by(Roach_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(Roach_ID, Pike_ID, max_encounter_date = Date, max_encounter_count = encounter_count)

#Add max_encounter_date and max_encounter_count to potential predation events summary.
BT_roach_pike_pred_int <- 
  BT_roach_pike_pred_int %>%
  left_join(BT_roach_max_encounter_dates %>%
              dplyr::select(Roach_ID, Pike_ID, max_encounter_date, max_encounter_count),
            by = c("Roach_ID" = "Roach_ID", 
                   "Pike_ID" = "Pike_ID"))

# Save the summarised potential predation events for future use
saveRDS(BT_roach_pike_pred_int, paste0(enc_path, "BT_roach_possible_predation_events.rds"))
saveRDS(BT_roach_max_encounter_dates, paste0(enc_path, "BT_roach_max_encounters_for_each_ID.rds"))


#### > 2.2. Perch - Pike encounters #### 

# This section calculates encounters between Perch and Pike, specifically focusing on close interactions
# (within 0.45 meters, which is the max strike distance of a Pike, according to Harper 1991).

# Identify encounters where Perch and Pike are within 0.45 meters
perch_pike_distances_df$encounter <- ifelse(perch_pike_distances_df$est <= 0.45, 1, 0)

# Add a date column based on the timestamp for further daily analysis
# Note takes time to run
perch_pike_distances_df <- 
  perch_pike_distances_df %>%
  mutate(Date = as.Date(format(perch_pike_distances_df$timestamp, tz = "Europe/Stockholm")))

#save with encounter and date column
saveRDS(perch_pike_distances_df, paste0(enc_path, "BT_pike_perch_distances_df.rds"))

# Summarize the number of encounters and calculate the average daily distance for each Perch-Pike pair, excluding any outliers
perch_pike_encounter_sum <- 
  perch_pike_distances_df %>%
  group_by(Perch_ID, Pike_ID, Date) %>%
  summarise(encounter_count = sum(encounter),
            daily_avg_dist = mean(est, na.rm = TRUE))   

#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#

####____ 2.2.2. Identify possible predation events ####

# Merge biometric data (e.g., whether the fish was found alive) with perch encounter summary
perch_post_size_cols <- 
  post_biometrics %>%
  filter(Lake == 'BT', Species == 'Perch') %>%
  dplyr::select(individual_ID, Found, Known_predated) %>% 
  rename(found_alive = Found)

perch_pike_encounter_sum <- 
  perch_pike_encounter_sum %>%
  left_join(perch_post_size_cols, by = c("Perch_ID" = "individual_ID"))

# Filter perch-pike encounters for days with more than 50 encounters or an avg daily distance < 3 meters.
perch_over_50_encounters <- 
  perch_pike_encounter_sum %>%
  filter(encounter_count >= 50  | daily_avg_dist < 3 | Known_predated == 1)

#how many unique perch individuals
unique(perch_over_50_encounters$Perch_ID)
#2 individuals

# Summarise perch-pike interactions: count days with > 50 encounters, consecutive days, and the first date
BT_perch_pike_pred_int <- 
  perch_over_50_encounters %>%
  group_by(Perch_ID, Pike_ID) %>%
  summarise(
    num_days_over_50 = n(),  # Count of days with over 50 encounters or daily avg distance < 3 meters
    encounter_dates = paste(format(Date, "%d/%m/%Y"), collapse = ", "),  # Combine encounter dates as a string
    consecutive_days = max(rle(cumsum(c(1, diff(Date) > 1)) == 1)$lengths),  # Calculate consecutive encounter days
    first_date_over_50 = min(Date)  # Identify first date with more than 50 encounters
  )

# Extract individual treatment and missing date information from filtered data
perch_sel_cols <- 
  BT_filt_data %>%
  filter(Species == 'Perch') %>%
  dplyr::select(individual_ID, Treatment, n_missing_dates) %>%
  distinct()

# Merge encounter summary with individual-level metadata
BT_perch_pike_pred_int <- 
  BT_perch_pike_pred_int %>%
  left_join(perch_sel_cols, by = c("Perch_ID" = "individual_ID"))

# Merge biometric data (e.g., whether the fish was found alive) with possible predation events
# Create column indicated instances where there might be a higher likelihood that an individual was predated
# specifically, when individuals 2 or more consecutive days where they had more than 50 encounters with a predator
BT_perch_pike_pred_int <- 
  BT_perch_pike_pred_int %>%
  left_join(perch_post_size_cols, by = c("Perch_ID" = "individual_ID")) %>%
  # Remove fish found alive at the end of the experiment
  filter(found_alive == '0') %>%
  # Flag likely predation events (i.e., when consecutive encounters are more than 2 days)
  mutate(likely_predated = ifelse(consecutive_days >= 2, 1, 0))


# Identify the date with the most encounters for each perch and pike pair
BT_perch_max_encounter_dates <- 
  perch_pike_encounter_sum %>%
  group_by(Perch_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(Perch_ID, Pike_ID, max_encounter_date = Date, max_encounter_count = encounter_count)

#Add max_encounter_date and max_encounter_count to potential predation events summary.
BT_perch_pike_pred_int <- 
  BT_perch_pike_pred_int %>%
  left_join(BT_perch_max_encounter_dates %>%
              dplyr::select(Perch_ID, Pike_ID, max_encounter_date, max_encounter_count),
            by = c("Perch_ID" = "Perch_ID", 
                   "Pike_ID" = "Pike_ID"))

# Save the summarised potential predation events for future use
saveRDS(BT_perch_pike_pred_int, paste0(enc_path, "BT_perch_possible_predation_events.rds"))
saveRDS(BT_perch_max_encounter_dates, paste0(enc_path, "BT_perch_max_encounters_for_each_ID.rds"))


#---------------------------------------------------------#
# 3. Create dataframe of likely predation events in BT ####
#---------------------------------------------------------#

# Combine predation events into a single dataframe
BT_pred_events <- 
  rbind(
    BT_roach_pike_pred_int %>% rename(individual_ID = Roach_ID) %>% mutate(Species = 'Roach'),
    BT_perch_pike_pred_int %>% rename(individual_ID = Perch_ID) %>% mutate(Species = 'Perch')
  ) %>%
  dplyr::select(individual_ID, Species, Treatment, everything())

# Save the combined predation event data
saveRDS(BT_pred_events, paste0(enc_path, "BT_possible_pred_events.rds"))


#---------------------------------------------------------#
# 4. Create summary of encounter rate unfiltered ##########
#---------------------------------------------------------#

# The purpose of this section is to evaluate whether the number of encounters might be inflated due prey being
# eaten by a pike but still being tracked (e.g., high encounter rates post-predation). The analysis aims to determine
# when such high encounter rates started and remove data after potential predation.

# > 4.1. Roach ####

# >____ Total encounters ####

# Create a summary dataframe for the number of predator encounters each roach had per day
# This will be used to help try and identify the fate of individuals that were not found at the end of the experiment
# see below

BT_roach_daily_encounter_summary <-
  roach_pike_encounter_sum %>%
  group_by(Roach_ID, Date) %>%
  summarise(total_daily_encounter_count = sum(encounter_count),
            max_daily_encounter_count = max(encounter_count),
            pike_id_max_encounter_count = Pike_ID[which.max(encounter_count)],
            min_avg_dist_from_pred = min(daily_avg_dist),
            pike_id_min_dist = Pike_ID[which.min(daily_avg_dist)],
            found_alive = first(found_alive),
            found_predated = first(Known_predated))


# Total encounter summary
# Not unfiltered because we have not removed encounters calculated post-predation events
# this these encounter numbers are likely to be inflated for individuals that were predated

BT_roach_encounter_summary <-
  BT_roach_daily_encounter_summary %>%
  group_by(Roach_ID) %>%
  summarise(total_encounters = sum(total_daily_encounter_count),
            avg_minimum_daily_dist_from_pike = mean(mean(min_avg_dist_from_pred, na.rm = TRUE)))

# Merge encounter summary with metadata for individual Roach
BT_roach_encounter_summary <-
  BT_roach_encounter_summary %>%
  left_join(roach_sel_cols, by = c("Roach_ID" = "individual_ID"))

#also add number of poor tracking days info
#the loaded dataframe was created in script 03_BT_tracking_data_summary
BT_ID_irreg_sampling <- readRDS(paste0(filtered_data_path, "BT_ID_irreg_sampling.rds"))
roach_irreg_sampling_cols <- 
  BT_ID_irreg_sampling %>% 
  filter(Species == 'Roach') %>% 
  select(individual_ID, n_poor_tracking_days, found_alive)

BT_roach_encounter_summary <-
  BT_roach_encounter_summary %>%
  left_join(roach_irreg_sampling_cols, by = c("Roach_ID" = "individual_ID"))

BT_roach_encounter_summary <-
  BT_roach_encounter_summary %>% 
  mutate(days_tracked = 34 - n_missing_dates,  # Total tracking duration (max 34 days)
       avg_daily_pred_encounter_rate = round(total_encounters / days_tracked, 0))

head(BT_roach_encounter_summary)


# Save the unfiltered encounter summary
saveRDS(BT_roach_encounter_summary, paste0(enc_path, 'BT_roach_encounter_summmary_unfiltered.rds'))


# >____ Encounter summaries for roach not found at the end of the experiment ####

#Reduce summary data to only include individuals that were not found at the end of the experiment
daily_encounter_sum_roach_not_found <- 
  BT_roach_daily_encounter_summary %>% 
  filter(found_alive == 0)

#Add column to identify whether the date was irregularly sampled or dates were not tracked
#First I need to extract the required columns
roach_not_found_cols <- 
  BT_ID_irreg_sampling %>% 
  filter(Species == 'Roach') %>% 
  select(individual_ID, dates_irregular_positions, missing_dates)

#Then we need to do some data wrangling to extract the date information

daily_encounter_sum_roach_not_found$Date <- as.Date(daily_encounter_sum_roach_not_found$Date)
# Ensure dates_irregular_positions and missing_dates are properly formatted as lists of dates
roach_not_found_cols$dates_irregular_positions <- lapply(strsplit(roach_not_found_cols$dates_irregular_positions, ", "), as.Date)
roach_not_found_cols$missing_dates <- lapply(strsplit(roach_not_found_cols$missing_dates, ", "), as.Date)

# Merge the two datasets on the IDs
merged_df <- merge(daily_encounter_sum_roach_not_found, roach_not_found_cols, 
                   by.x = "Roach_ID", by.y = "individual_ID", all.x = TRUE)

# Create the two new columns
merged_df$poor_tracking_date <- mapply(function(date, irregular_dates) {
  if (!is.null(irregular_dates) && date %in% irregular_dates) 1 else 0
}, merged_df$Date, merged_df$dates_irregular_positions)

merged_df$no_tracking_date <- mapply(function(date, missing_dates) {
  if (!is.null(missing_dates) && date %in% missing_dates) 1 else 0
}, merged_df$Date, merged_df$missing_dates)

# Keep only the original columns from daily_encounter_sum_roach_found and the new ones
daily_encounter_sum_roach_not_found <- merged_df[, c("Roach_ID", "Date", "total_daily_encounter_count", 
                                                     "max_daily_encounter_count", "pike_id_max_encounter_count",
                                                     "min_avg_dist_from_pred", "pike_id_min_dist", "found_predated",
                                                     "poor_tracking_date", "no_tracking_date")]

# Save the unfiltered encounter summary
saveRDS(daily_encounter_sum_roach_not_found, paste0(enc_path, 'BT_daily_encounter_summaries_for_roach_not_found.rds'))


#Total encounter summary for roach not found
#This summarised dataframe will include information on
#total encounter count
#max encounter count
#date with max encounter count
#number of days with over 50 encounters
#first date with over 50 encounters

encounter_sum_roach_not_found <-
  daily_encounter_sum_roach_not_found %>% 
  group_by(Roach_ID) %>% 
  summarise(total_encounter_count = sum(total_daily_encounter_count),
            max_encounter_count = max(max_daily_encounter_count ),
            max_encounter_date = Date[which.max(max_daily_encounter_count )],
            max_encounter_pike = pike_id_max_encounter_count[which.max(max_daily_encounter_count)])

roach_not_found_cols <-
  BT_roach_pike_pred_int %>%
  dplyr::select(Roach_ID, num_days_over_50, first_date_over_50, consecutive_days) %>% 
  rename(consecutive_days_over_50 = consecutive_days)

encounter_sum_roach_not_found <-
  encounter_sum_roach_not_found %>%
  left_join(roach_not_found_cols, by = c("Roach_ID" = "Roach_ID"))

saveRDS(encounter_sum_roach_not_found, paste0(enc_path, 'BT_total_encounter_summary_for_roach_not_found.rds'))

#--------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------#

#> 4.2. Perch ####

# >____ Total encounters ####

BT_perch_daily_encounter_summary <-
  perch_pike_encounter_sum %>%
  group_by(Perch_ID, Date) %>%
  summarise(total_daily_encounter_count = sum(encounter_count),
            max_daily_encounter_count = max(encounter_count),
            pike_id_max_encounter_count = Pike_ID[which.max(encounter_count)],
            min_avg_dist_from_pred = min(daily_avg_dist),
            pike_id_min_dist = Pike_ID[which.min(daily_avg_dist)],
            found_alive = first(found_alive),
            found_predated = first(Known_predated))


# Total encounter summary
# Not unfiltered because we have not removed encounters calculated post-predation events
# this these encounter numbers are likely to be inflated for individuals that were predated

BT_perch_encounter_summary <-
  BT_perch_daily_encounter_summary %>%
  group_by(Perch_ID) %>%
  summarise(total_encounters = sum(total_daily_encounter_count),
            avg_minimum_daily_dist_from_pike = mean(mean(min_avg_dist_from_pred, na.rm = TRUE)))

# Merge encounter summary with metadata for individual perch
BT_perch_encounter_summary <-
  BT_perch_encounter_summary %>%
  left_join(perch_sel_cols, by = c("Perch_ID" = "individual_ID"))

perch_irreg_sampling_cols <- 
  BT_ID_irreg_sampling %>% 
  filter(Species == 'Perch') %>% 
  select(individual_ID, n_poor_tracking_days, found_alive)

BT_perch_encounter_summary <-
  BT_perch_encounter_summary %>%
  left_join(perch_irreg_sampling_cols, by = c("Perch_ID" = "individual_ID"))

BT_perch_encounter_summary <-
  BT_perch_encounter_summary %>% 
  mutate(days_tracked = 34 - n_missing_dates,  # Total tracking duration (max 34 days)
         avg_daily_pred_encounter_rate = round(total_encounters / days_tracked, 0))

head(BT_perch_encounter_summary)


# Save the unfiltered encounter summary
saveRDS(BT_perch_encounter_summary, paste0(enc_path, 'BT_perch_encounter_summmary_unfiltered.rds'))


# >____ Encounter summaries for perch not found at the end of the experiment ####

#Reduce summary data to only include individuals that were not found at the end of the experiment
daily_encounter_sum_perch_not_found <- 
  BT_perch_daily_encounter_summary %>% 
  filter(found_alive == 0)

#Add column to identify whether the date was irregularly sampled or dates were not tracked
#First I need to extract the required columns
perch_not_found_cols <- 
  BT_ID_irreg_sampling %>% 
  filter(Species == 'Perch') %>% 
  select(individual_ID, dates_irregular_positions, missing_dates)

#Then we need to do some data wrangling to extract the date information

daily_encounter_sum_perch_not_found$Date <- as.Date(daily_encounter_sum_perch_not_found$Date)
# Ensure dates_irregular_positions and missing_dates are properly formatted as lists of dates
perch_not_found_cols$dates_irregular_positions <- lapply(strsplit(perch_not_found_cols$dates_irregular_positions, ", "), as.Date)
perch_not_found_cols$missing_dates <- lapply(strsplit(perch_not_found_cols$missing_dates, ", "), as.Date)

# Merge the two datasets on the IDs
merged_df <- merge(daily_encounter_sum_perch_not_found, perch_not_found_cols, 
                   by.x = "Perch_ID", by.y = "individual_ID", all.x = TRUE)

# Create the two new columns
merged_df$poor_tracking_date <- mapply(function(date, irregular_dates) {
  if (!is.null(irregular_dates) && date %in% irregular_dates) 1 else 0
}, merged_df$Date, merged_df$dates_irregular_positions)

merged_df$no_tracking_date <- mapply(function(date, missing_dates) {
  if (!is.null(missing_dates) && date %in% missing_dates) 1 else 0
}, merged_df$Date, merged_df$missing_dates)

# Keep only the original columns from daily_encounter_sum_perch_found and the new ones
daily_encounter_sum_perch_not_found <- merged_df[, c("Perch_ID", "Date", "total_daily_encounter_count", 
                                                     "max_daily_encounter_count", "pike_id_max_encounter_count",
                                                     "min_avg_dist_from_pred", "pike_id_min_dist", "found_predated",
                                                     "poor_tracking_date", "no_tracking_date")]

# Save the unfiltered encounter summary
saveRDS(daily_encounter_sum_perch_not_found, paste0(enc_path, 'BT_daily_encounter_summaries_for_perch_not_found.rds'))


#Total encounter summary for perch not found
#This summarised dataframe will include information on
#total encounter count
#max encounter count
#date with max encounter count
#number of days with over 50 encounters
#first date with over 50 encounters

encounter_sum_perch_not_found <-
  daily_encounter_sum_perch_not_found %>% 
  group_by(Perch_ID) %>% 
  summarise(total_encounter_count = sum(total_daily_encounter_count),
            max_encounter_count = max(max_daily_encounter_count ),
            max_encounter_date = Date[which.max(max_daily_encounter_count )],
            max_encounter_pike = pike_id_max_encounter_count[which.max(max_daily_encounter_count)])

perch_not_found_cols <-
  BT_perch_pike_pred_int %>%
  dplyr::select(Perch_ID, num_days_over_50, first_date_over_50, consecutive_days) %>% 
  rename(consecutive_days_over_50 = consecutive_days)

encounter_sum_perch_not_found <-
  encounter_sum_perch_not_found %>%
  left_join(perch_not_found_cols, by = c("Perch_ID" = "Perch_ID"))

saveRDS(encounter_sum_perch_not_found, paste0(enc_path, 'BT_total_encounter_summary_for_perch_not_found.rds'))

#--------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------#

#Combine roach and perch encounter summaries for those individuals not found at the end of the study

# Modify the roach dataset
encounter_sum_roach_not_found <- 
  encounter_sum_roach_not_found %>%
  dplyr::rename(individual_ID = Roach_ID) %>%  # Rename roach_id to individual_ID
  dplyr::mutate(Species = "Roach")            # Add a Species column with "Roach"

# Modify the perch dataset
encounter_sum_perch_not_found <- 
  encounter_sum_perch_not_found %>%
  dplyr::rename(individual_ID = Perch_ID) %>% # Rename perch_ID to individual_ID
  dplyr::mutate(Species = "Perch")            # Add a Species column with "Perch"

# Combine the two datasets
encounters_sum_IDs_not_found <- dplyr::bind_rows(encounter_sum_roach_not_found, encounter_sum_perch_not_found)

#remove duplicate rows, retaining the row that has higest num_days_over_50
filtered_encounters <- 
  encounters_sum_IDs_not_found %>%
  group_by(individual_ID) %>%                               # Group by individual_ID
  arrange(desc(first_date_over_50),                        # Sort by most recent date
          desc(num_days_over_50),                          # Then by higher num_days_over_50
          .by_group = TRUE) %>%                            # Sort within each group
  slice(1) %>%                                             # Select the first row of each group
  ungroup()                                                # Ungroup the data


#merge with irregular sampling summary
#create in script 03_BT_tracking_data_summary
BT_ID_irreg_sampling_not_found <- 
  BT_ID_irreg_sampling %>% 
  filter(found_alive == 0)

#get required rows from filtered encounters
filtered_encounters_rows <- 
  filtered_encounters %>% 
  select(individual_ID, total_encounter_count, max_encounter_count, max_encounter_date, first_date_over_50, num_days_over_50, consecutive_days_over_50)

BT_ID_irreg_sampling_not_found <-
  BT_ID_irreg_sampling_not_found %>%
  left_join(filtered_encounters_rows, by = c("individual_ID" = "individual_ID")) %>% 
  select(-found_alive)

# Create a flextable summarizing the irregular and missing days for missing individuals
daily_irregular_sample_table <- 
  flextable(BT_ID_irreg_sampling_not_found) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("Species" = 'Species',
                    "individual_ID" = 'ID', 
                    "n_poor_tracking_days" = "N poor tracking days",
                    "total_missing_days" = "N days not tracked",
                    "total_irregular_or_missing_days" = "Total",
                    "dates_irregular_positions" = "Dates poor tracking",
                    "missing_dates" = 'Dates ID not tracked',
                    "Known_predated" = "Found predated (yes = 1)",
                    "total_encounter_count" = "Total pike encounters",
                    "max_encounter_count" = "Max encounters in day",
                    "max_encounter_date" = "Date with max encounters",
                    "first_date_over_50" = "First date with over 50 encounters",
                    "num_days_over_50" = "N days over 50 encounters",
                    "consecutive_days_over_50" = "consecutive days with 50 encounters"
                    )

save_as_docx(daily_irregular_sample_table, 
             path = paste0(save_tables_path, "BT_IDs_not_found_encounter and tracking_summary.docx"))

#--------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------#


# I will need to filter the encounter date once we decide who was predated/died and when they predated/died








# Create a dataframe that filters out potential post-predation encounters by removing rows after the first
# date with more than 50 encounters for likely predated Roach.

# Filter dataframe to remove post-predation encounters
pred_cols <-
  BT_roach_pike_pred_int %>%
  filter(likely_predated == '1') %>%
  dplyr::select(Roach_ID, likely_predated, first_date_over_50)

roach_pike_encounter_sum_filt <-
  roach_pike_encounter_sum %>%
  left_join(pred_cols, by = "Roach_ID") %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%
  filter(is.na(first_date_over_50) | Date <= first_date_over_50)  # Keep only pre-predation data

# Recalculate the most encounter-heavy date for each Roach-Pike pair after filtering
max_encounter_dates_filt <-
  roach_pike_encounter_sum_filt %>%
  group_by(Roach_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(Roach_ID, Pike_ID, max_encounter_date = Date, max_encounter_count = encounter_count)

# Recreate the summary table with filtered data
roach_pike_encounter_sum_filt <-
  roach_pike_encounter_sum_filt %>%
  group_by(Roach_ID, Pike_ID) %>%
  mutate(num_encounters_w_pike_id = sum(encounter_count)) %>%
  group_by(Roach_ID) %>%
  mutate(total_pike_encounters = sum(encounter_count),
         avg_daily_dist_from_pike = mean(daily_avg_dist, na.rm = TRUE)) %>%
  left_join(max_encounter_dates_filt, by = c("Roach_ID", "Pike_ID"))


# Extract individual treatment and missing date information from filtered data
sel_cols <- 
  BT_filt_data %>%
  filter(Species == 'Roach') %>%
  dplyr::select(individual_ID, Treatment, n_missing_dates) %>%
  distinct()

# Merge filtered encounter summary with metadata for individual Roach
roach_pike_encounter_sum_filt <-
  roach_pike_encounter_sum_filt %>%
  left_join(sel_cols, by = c("Roach_ID" = "individual_ID"))

# Create a smaller summary table for further analysis
roach_pike_encounter_sum_filt <- 
  roach_pike_encounter_sum_filt %>%
  dplyr::select(Roach_ID, Treatment, total_pike_encounters, avg_daily_dist_from_pike, n_missing_dates) %>% 
  distinct()  # Remove duplicate rows to ensure one entry per Perch

# Add information about whether each Perch was found alive at the end of the experiment
post_size_cols <- 
  post_biometrics %>%
  filter(Lake == 'BT', Species == 'Roach') %>%
  dplyr::select(individual_ID, Found, Known_predated) %>% 
  rename(found_alive = Found)

# Add information about whether the Perch were found alive at the end of the experiment
roach_pike_encounter_sum_filt <- 
  roach_pike_encounter_sum_filt %>%
  left_join(post_size_cols, by = c("Roach_ID" = "individual_ID"))

# Add information about likely predation events
pred_cols <- 
  BT_roach_pike_pred_int %>%
  filter(likely_predated == '1') %>%
  dplyr::select(Roach_ID, likely_predated)

# Merge predation information into the summary table
roach_pike_encounter_sum_filt <- 
  roach_pike_encounter_sum_filt %>%
  left_join(pred_cols, by = c("Roach_ID" = "Roach_ID")) %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%  # Replace NA values with 0 for non-predated fish
  mutate(likely_died = ifelse(n_missing_dates > 3 & found_alive == 0 | likely_predated == 1, 1, 0),  # Determine if fish likely died
         days_tracked = 36 - n_missing_dates,  # Calculate total days tracked (36 days max)
         avg_pred_daily_encounter_rate = round(total_pike_encounters / days_tracked, 0)) %>%  # Average daily encounter rate with Pike
  group_by(Roach_ID) %>%
  distinct()  # Remove duplicate rows


# Save the filtered summary table for Roach-Pike encounters
saveRDS(roach_pike_encounter_sum_filt, paste0(enc_path, 'roach_pike_encounter_sum_filt.rds'))




# > 4.2. Perch encounters ####         

## Filtered ##
# Filter dataframe to remove post-predation encounters after the first day with >50 encounters
pred_cols <- 
  BT_perch_pike_pred_int %>%
  filter(consecutive_days > 2) %>%
  dplyr::select(perch_ID, likely_predated, first_date_over_50)

perch_pike_encounter_sum_filt <- 
  perch_pike_encounter_sum %>%
  left_join(pred_cols, by = "perch_ID") %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%
  # Keep only rows before or on the first day with more than 50 encounters
  filter(is.na(first_date_over_50) | Date <= first_date_over_50)

# Identify the date with the most encounters for each Perch and Pike pair - filtered dataset
# This section calculates the day with the highest number of encounters between Perch and Pike after removing likely predation events.
max_encounter_dates_filt <- 
  perch_pike_encounter_sum_filt %>%
  group_by(perch_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%  # If there are multiple dates with the same max encounter count, select the first one
  dplyr::select(perch_ID, Pike_ID, max_encounter_date = Date, max_encounter_count = encounter_count)

# Summarize filtered encounter data for Perch-Pike interactions
perch_pike_encounter_sum_filt <- 
  perch_pike_encounter_sum_filt %>%
  group_by(perch_ID, Pike_ID) %>%
  mutate(num_encounters_w_pike_id = sum(encounter_count)) %>%  # Total number of encounters between each Perch and Pike
  group_by(perch_ID) %>%
  mutate(total_pike_encounters = sum(encounter_count),  # Total encounters for each Perch with all Pike
         avg_daily_dist_from_pike = mean(daily_avg_dist, na.rm = TRUE)) %>%  # Average distance between each Perch and Pike
  left_join(max_encounter_dates_filt, by = c("perch_ID", "Pike_ID")) %>%  # Add the date and count of the max encounters
  dplyr::select(-likely_predated, -first_date_over_50)  # Remove columns related to predation filtering

# Extract additional metadata such as treatment type and number of missing tracking dates for each Perch individual
sel_cols <- 
  BT_filt_data %>%
  filter(Species == 'Perch') %>%
  dplyr::select(individual_ID, Treatment, n_missing_dates) %>%
  distinct()

# Add metadata such as treatment type and number of missing tracking dates
perch_pike_encounter_sum_filt <- 
  perch_pike_encounter_sum_filt %>%
  left_join(sel_cols, by = c("perch_ID" = "individual_ID"))

# Create a smaller summary table for further analysis
perch_pike_encounter_sum_filt <- 
  perch_pike_encounter_sum_filt %>%
  dplyr::select(perch_ID, Treatment, total_pike_encounters, avg_daily_dist_from_pike, n_missing_dates) %>% 
  distinct()  # Remove duplicate rows to ensure one entry per Perch

# Add information about whether each Perch was found alive at the end of the experiment
post_size_cols <- 
  post_biometrics %>%
  filter(Lake == 'BT', Species == 'Perch') %>%
  dplyr::select(individual_ID, Found, Known_predated) %>% 
  rename(found_alive = Found)

# Add information about whether the Perch were found alive at the end of the experiment
perch_pike_encounter_sum_filt <- 
  perch_pike_encounter_sum_filt %>%
  left_join(post_size_cols, by = c("perch_ID" = "individual_ID"))

# Add information about likely predation events
pred_cols <- 
  BT_perch_pike_pred_int %>%
  filter(likely_predated == '1') %>%
  dplyr::select(perch_ID, likely_predated)

# Merge predation information into the summary table
perch_pike_encounter_sum_filt <- 
  perch_pike_encounter_sum_filt %>%
  left_join(pred_cols, by = c("perch_ID" = "perch_ID")) %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%  # Replace NA values with 0 for non-predated fish
  mutate(likely_died = ifelse(n_missing_dates > 3 & found_alive == 0 | likely_predated == 1, 1, 0),  # Determine if fish likely died
         days_tracked = 36 - n_missing_dates,  # Calculate total days tracked (36 days max)
         avg_pred_daily_encounter_rate = round(total_pike_encounters / days_tracked, 0)) %>%  # Average daily encounter rate with Pike
  group_by(perch_ID) %>%
  distinct()  # Remove duplicate rows

# Save the filtered summary table for Perch-Pike encounters
saveRDS(perch_pike_encounter_sum_filt, paste0(enc_path, 'perch_pike_encounter_sum_filt.rds'))


# # Create a summary dataframe with unfiltered data (i.e., all encounters)
# perch_pike_encounter_sum_unfilt <- 
#   perch_pike_encounter_sum %>%
#   group_by(perch_ID, Pike_ID) %>%
#   mutate(num_encounters_w_pike_id = sum(encounter_count)) %>%
#   group_by(perch_ID) %>%
#   mutate(total_pike_encounters = sum(num_encounters_w_pike_id),
#          avg_dist_from_pike = mean(daily_avg_dist, na.rm = TRUE)) %>%
#   left_join(max_encounter_dates, by = c("perch_ID", "Pike_ID"))
# 
# # Merge the encounter summary with metadata about each individual Perch
# perch_pike_encounter_sum_unfilt <- perch_pike_encounter_sum_unfilt %>%
#   left_join(sel_cols, by = c("perch_ID" = "individual_id"))
# 
# # Add information on whether each fish was found alive or likely predated
# perch_pike_encounter_sum_unfilt <- perch_pike_encounter_sum_unfilt %>%
#   left_join(post_size_cols, by = c("perch_ID" = "individual_id"))
# 
# # Merge information about likely predation events
# pred_cols <- BT_perch_pike_pred_int %>%
#   filter(likely_predated == '1') %>%
#   dplyr::select(perch_ID, likely_predated)
# 
# perch_pike_encounter_sum_unfilt <- perch_pike_encounter_sum_unfilt %>%
#   left_join(pred_cols, by = c("perch_ID" = "perch_ID")) %>%
#   tidyr::replace_na(list(likely_predated = 0)) %>%
#   # Flag individuals as likely dead if they were predated or had more than 3 missing tracking days
#   mutate(likely_died = ifelse(n_missing_dates > 3 & Found == 0 | likely_predated == 1, 1, 0),
#          days_tracked = 36 - n_missing_dates,  # Total tracking duration (max 36 days)
#          avg_pred_encounter_rate = round(total_pike_encounters / days_tracked, 0))
# 
# # Save the unfiltered encounter summary
# saveRDS(perch_pike_encounter_sum_unfilt, paste0(enc_path, 'perch_pike_encounter_sum_unfilt.rds'))


# > 4.3. BT encounter summary ####

# Combine Roach-Pike and Perch-Pike encounter data into a single summary table of predation events in BT lake.
BT_pred_encounter_summary <- 
  rbind(
    roach_pike_encounter_sum_filt %>% rename(individual_ID = Roach_ID) %>% mutate(Species = 'Roach'),
    perch_pike_encounter_sum_filt %>% rename(individual_ID = perch_ID) %>% mutate(Species = 'Perch')
  ) %>%
  dplyr::select(individual_ID, Species, Treatment, everything()) %>%
  tidyr::replace_na(list(Known_predated = 0))  # Replace NA values in Known_predated column

# Save the combined predation event summary
saveRDS(BT_pred_encounter_summary, paste0(enc_path, "BT_pred_encounter_summary.rds"))

#----------------------------------------------------------------------------------------------------------#
#---------------------------------------------------#
######### 5. Visualise predator encounters #########
#--------------------------------------------------#

#Filter IDs that have been identified as being predated
pred_id_dist <- BT_pike_roach_distances_df %>% 
  filter(Roach_ID %in% c("F59731", "F59719", "F59701", "F59738", "F59729"))

#make a date column
pred_id_dist <- pred_id_dist %>% 
  mutate(date = as.Date(timestamp, tz = "Europe/Stockholm"))  # Extract date from UTC+2 timestamp

#Add number of encounters per date
# Define encounters as those occurring within 0.45 meters (based on Pike's strike distance)
pred_id_dist$encounter <- ifelse(pred_id_dist$est <= 0.45, 1, 0)


#extract identified pred interactions
F59719_dist <- pred_id_dist %>% 
  filter(Roach_ID == 'F59719', Pike_ID == 'F59885') %>% 
  filter(date <= "2022-10-21") #for graphing purposes only

F59701_dist <- pred_id_dist %>% 
  filter(Roach_ID == 'F59701', Pike_ID == 'F59885') %>% 
  filter(date <= "2022-10-22")



#calculate average dist per date
F59719_dist_avg <- F59719_dist %>% 
  group_by(date) %>% 
  summarise(avg_dist = mean(est, na.rm = TRUE),
            encounters = sum(encounter))

F59701_dist_avg <- F59701_dist %>% 
  group_by(date) %>% 
  summarise(avg_dist = mean(est, na.rm = TRUE),
            encounters = sum(encounter))





library(ggplot2)

# Calculate scaling factor
scale_factor <- max(F59719_dist_avg$avg_dist, na.rm = TRUE) / max(F59719_dist_avg$encounters, na.rm = TRUE)


roach_pike_dist_plot <- ggplot(F59719_dist_avg, aes(x = date)) +
  # Bar chart for encounters, scaled and made semi-transparent
  geom_bar(aes(y = encounters * scale_factor), 
           stat = "identity", 
           fill = "grey", 
           alpha = 0.5) +
  
  # Line graph for average distance
  geom_line(aes(y = avg_dist), 
            color = "black", 
            size = 1) +
  geom_point(aes(y = avg_dist), 
             color = "black", 
             size = 2) +
  
  # Text labels for actual encounters on top of bars
  geom_text(aes(y = encounters * scale_factor, 
                label = encounters),
            vjust = -0.5, 
            size = 3) +
  
  # Primary y-axis for average distance
  scale_y_continuous(
    name = "Average distance from pike (m)",
    # Adjust the limits if necessary
    limits = c(0, max(c(F59719_dist_avg$avg_dist, F59719_dist_avg$encounters * scale_factor), na.rm = TRUE) * 1.1)
  ) +
  
  # Labels and theme
  labs(
    x = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1)
  )

ggsave(file = paste0(enc_path, "encounter_plots/roach_pike_dist_plot.png"), 
       plot = roach_pike_dist_plot, 
       device = 'png',
       width = 12.5, 
       height = 10,
       units = 'cm',
       dpi = 300)



### Visualise predator-prey interaction ###

#Extract data for individuals and dates
F59709_pred_dat <- BT_filt_data %>% 
  filter(individual_id %in% c('F59719','F59885')) %>%
  mutate(date_utc2 = as.Date(timestamp, tz = "Europe/Stockholm")) %>% 
  filter(date_utc2 <= "2022-10-22")





# now create a move2 object from a data.frame
F59709_pred_mv <- 
  mt_as_move2(F59709_pred_dat, 
              coords = c("longitude","latitude"),
              crs = "WGS84",
              time_column = "timestamp",
              track_id_column = "individual_id",
              na.fail = F) # allows or not empty coordinates


#Plot tracks over outline of lake
#Load map polygon
BT_polygon = sf::st_read(paste0(lake_polygon_path, "lake_BT_polygon.gpkg"))




#add scale bar
library(sf)
library(ggsn)
library(gganimate)

# Check the CRS of your lake polygon
st_crs(BT_polygon)

# Define the length of the scale bar in meters
scale_length <- 100  # 1 km

# Get the bounding box of the lake to determine where to place the scale bar
lake_bbox <- st_bbox(BT_polygon)

# Define the starting point of the scale bar near the bottom left corner of the map
x_start <- lake_bbox["xmin"] + (lake_bbox["xmax"] - lake_bbox["xmin"]) * 0.1
y_start <- lake_bbox["ymin"] + (lake_bbox["ymax"] - lake_bbox["ymin"]) * 0.1

# Define the ending point of the scale bar
x_end <- x_start + scale_length

#filter for plotting
F59709_pred_dat_filt <- F59709_pred_dat %>% 
  filter(date_utc2 > "2022-10-03")

#subset the data

# Define your ggplot object with the necessary layers
p <- ggplot() + 
  geom_sf(data = BT_polygon, fill = NA, color = "black", linewidth = 1) +
  geom_point(data = F59709_pred_dat_filt, 
             aes(x = longitude, 
                 y = latitude, 
                 group = individual_id,
                 color = individual_id), 
             size = 4) +  # Plot the fish tracks with color by individual_id
  ggsn::scalebar(
    data = BT_polygon, 
    dist = 5,  # Distance represented by the scale bar
    dist_unit = "m",  # Units of measurement
    transform = TRUE,  # Transform coordinates to projected CRS if needed
    model = 'WGS84',  # Model used for transformation, if needed
    location = "bottomright",  # Position of the scale bar
    st.dist = 0.05  # Distance of the text from the scale bar
  ) +
  theme_classic() + 
  transition_reveal(along = timestamp) +  # Reveal points along timestamp for all individuals
  ease_aes('linear') +
  labs(
    title = '{ifelse(format(as.POSIXct(frame_along), "%Y-%m-%d") == "2002-10-08", 
                   "Date: 2002-10-08 (Predation attempt occurred)", 
                   paste("Date:", format(as.POSIXct(frame_along), "%Y-%m-%d")))}',
    subtitle = '{ifelse(format(as.POSIXct(frame_along), "%Y-%m-%d") >= "2022-10-08" & 
                       format(as.POSIXct(frame_along), "%Y-%m-%d") <= "2022-10-21", 
                       "Roach predated!", 
                       "")}',
    x = "",
    y = ""
  ) +  # Show the date in the title and conditional subtitle
  scale_color_manual(values = c("Roach" = "blue", "Pike" = "red")) +  # Custom colors for Roach and Pike
  scale_color_discrete(name = "Species") +  # Adjust color legend for Species
  theme(
    plot.title = element_text(face = 'bold', size = 24),
    plot.subtitle = element_text(face = 'bold', color = 'red', size = 24, hjust = 0.5),  # Style for subtitle
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.position = 'none',
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )


p_animate <- animate(p, nframes = 500, fps = 12, width = 600, height = 600)

anim_save(paste0(enc_path, "roach_pred_animation.gif"), animation = p_animate)







?anim_save()





# #------------------------------------------------#
# ######## 5. Proximity calculation ################
# #------------------------------------------------#
# 
# #proximity function
# #Outputs a ratio estimate with confidence intervals
# #Value <1 indicate that the two individuals are closer on average than expected for independent movement
# #1 is consistent with independent movement
# #>1 individuals are farther from each other on average than expect for independent movement
# #TAKES A VERY LONG TIME TO RUN. 
# #Test proximity function
# #roach F59701 and  pike F59885
# 
# 
# # Reproject Roach telemetry objects to match the Pike projection
# projection(roach_BT_tel) <- projection(pike_BT_tel)
# projection(roach_BT_ctmm_fits) <- projection(pike_BT_ctmm_fits)
# 
# 
# combined_telemetry <- c(roach_BT_tel["F59701"], pike_BT_tel["F59885"])
# combined_ctmm <- c(roach_BT_ctmm_fits["F59701"], pike_BT_ctmm_fits["F59885"])
# 
# 
# #Pairwise separation distances
# PROX <- proximity(combined_telemetry,
#                   combined_ctmm)