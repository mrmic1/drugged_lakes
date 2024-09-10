# -----------------------------------------------# 
# Exploring interspecific interactions - Muddyfoot 
# ------------------------------------------------#

# This script analyses telemetry data of prey (Roach, Perch) and predators (Pike)
# in the Muddyfoot lake, including pairwise distances and encounter rates.
# The goal is to identify potential predation events based on close proximity and to create 
# a dataframe that can be used for analysing potential treatment effects on predator encounters

### LIBRARIES ###
# Loading required packages
`library(dplyr)
library(move)
library(move2)
library(ctmm)
library(lubridate)
library(data.table)
library(sf)  # For spatial data handling

####### 1. SETUP ########

### DIRECTORIES ###
# Paths to directories containing the necessary data files
data_filter_path <- "./data/tracks_filtered/"
lake_polygon_path <- "./data/lake_coords/"
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/"
enc_path <- "./data/encounters/"
size_path <- "./data/fish_size/"

### DATA LOADING ###
# Load preprocessed telemetry and model fit data for Roach, Perch, and Pike species
muddyfoot_filt_data <-  readRDS(paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))

pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds"))

# Load lake polygon data
muddyfoot_polygon <- sf::st_read(paste0(lake_polygon_path, "muddyfoot/lake_muddyfoot_polygon.gpkg"))

#------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------#
######## 2. Calculate distances at shared timepoints for prey and predators #######
#---------------------------------------------------------------------------------#

# This section calculates pairwise distances between Roach/Perch (prey) and Pike (predator) at shared time points.

# Reproject Roach telemetry objects to match the Pike projection
projection(roach_muddyfoot_tel) <- projection(pike_muddyfoot_tel)
projection(roach_muddyfoot_ctmm_fits) <- projection(pike_muddyfoot_ctmm_fits)

# Check projections to ensure they match
stopifnot(projection(roach_muddyfoot_tel) == projection(pike_muddyfoot_tel))

# Initialize an empty list to store distance calculations between Roach and Pike
roach_pike_distances <- list()

# Nested loop to calculate distances between each Roach and Pike telemetry pair
for(i in 1:length(roach_muddyfoot_tel)) {
  for(j in 1:length(pike_muddyfoot_tel)) {
    combined_telemetry <- c(roach_muddyfoot_tel[i], pike_muddyfoot_tel[j])
    combined_ctmm <- c(roach_muddyfoot_ctmm_fits[i], pike_muddyfoot_ctmm_fits[j])
    
    # Calculate pairwise distances
    location_difference <- distances(combined_telemetry, combined_ctmm)
    
    # Extract and store results with species IDs
    roach_id <- names(roach_muddyfoot_tel)[i]
    pike_id <- names(pike_muddyfoot_tel)[j]
    
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
saveRDS(roach_pike_distances_df, paste0(enc_path, "muddyfoot_pike_roach_distances_df.rds"))

# Repeat for Perch - Pike distances
# Reproject Perch telemetry objects to match Pike
projection(perch_muddyfoot_tel) <- projection(pike_muddyfoot_tel)
projection(perch_muddyfoot_ctmm_fits) <- projection(pike_muddyfoot_ctmm_fits)

# Ensure projections are correctly matched
stopifnot(projection(perch_muddyfoot_tel) == projection(pike_muddyfoot_tel))

# Initialize an empty list to store distance calculations between Perch and Pike
perch_pike_distances <- list()

# Nested loop for calculating Perch - Pike distances
for(i in 1:length(perch_muddyfoot_tel)) {
  for(j in 1:length(pike_muddyfoot_tel)) {
    combined_telemetry <- c(perch_muddyfoot_tel[i], pike_muddyfoot_tel[j])
    combined_ctmm <- c(perch_muddyfoot_ctmm_fits[i], pike_muddyfoot_ctmm_fits[j])
    
    location_difference <- distances(combined_telemetry, combined_ctmm)
    
    perch_id <- names(perch_muddyfoot_tel)[i]
    pike_id <- names(pike_muddyfoot_tel)[j]
    
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
saveRDS(perch_pike_distances_df, paste0(enc_path, "muddyfoot_pike_perch_distances_df.rds"))

#----------------------------------------------------------------------------------------------------#

#---------------------------------------------------#
######### 3. Calculate empirical encounters #########
#--------------------------------------------------#

#load distance df if not already
roach_pike_distances_df <- readRDS(paste0(enc_path, "muddyfoot_pike_roach_distances_df.rds"))
perch_pike_distances_df <- readRDS(paste0(enc_path, "muddyfoot_pike_perch_distances_df.rds"))

#> 3.1. Roach - Pike encounters #### 

# Define encounters as those occurring within 0.45 meters (based on Pike's strike distance)
roach_pike_distances_df$encounter <- ifelse(roach_pike_distances_df$est <= 0.45, 1, 0)

# Add date column based on timestamp
roach_pike_distances_df <- roach_pike_distances_df %>%
  mutate(date = as.Date(timestamp))

# Summarize encounter data per day and calculate average daily distance
roach_pike_encounter_sum <- roach_pike_distances_df %>%
  group_by(Roach_ID, Pike_ID, date) %>%
  summarize(encounter_count = sum(encounter),
            daily_avg_dist = mean(est, na.rm = TRUE)) %>%
  filter(date != '2022-09-24')  # Exclude outlier date


# >> 3.1.1. Identify possible predation events ####

# This section identifies potential predation events by filtering for days when Pike and Roach
# had more than 50 encounters or when their average daily distance was less than 3 meters.

# Filter Roach-Pike encounters for days with more than 50 encounters or an avg daily distance < 3 meters.
over_50_encounters <- roach_pike_encounter_sum %>%
  filter(encounter_count >= 50 | daily_avg_dist < 3)

# Summarize Roach-Pike interactions: count days with >50 encounters, consecutive days, and the first date
muddyfoot_roach_pike_pred_int <- over_50_encounters %>%
  group_by(Roach_ID, Pike_ID) %>%
  summarise(
    num_days_over_50 = n(),  # Count of days with over 50 encounters or daily avg distance < 3 meters
    encounter_dates = paste(format(date, "%d/%m/%Y"), collapse = ", "),  # Combine encounter dates as a string
    consecutive_days = max(rle(cumsum(c(1, diff(date) > 1)) == 1)$lengths),  # Calculate consecutive encounter days
    first_date_over_50 = min(date)  # Identify first date with more than 50 encounters
  )

# Extract individual treatment and missing date information from filtered data
sel_cols <- muddyfoot_filt_data %>%
  filter(Species == 'Roach') %>%
  dplyr::select(individual_id, treatment, n_missing_dates) %>%
  distinct()

# Merge encounter summary with individual-level metadata
muddyfoot_roach_pike_pred_int <- muddyfoot_roach_pike_pred_int %>%
  left_join(sel_cols, by = c("Roach_ID" = "individual_id"))

# Load post-experiment biometric data to assess known deaths or survivals
post_biometrics <- fread(paste0(size_path, "biometric_post_exp_data.csv")) %>%
  mutate(individual_id = paste0("F", sub(".*-", "", Tag_Number)))

# Merge biometric data (e.g., whether the fish was found alive) with encounter summary
post_size_cols <- post_biometrics %>%
  filter(Lake == 'Muddyfoot', Species == 'Roach') %>%
  dplyr::select(individual_id, Found)

muddyfoot_roach_pike_pred_int <- muddyfoot_roach_pike_pred_int %>%
  left_join(post_size_cols, by = c("Roach_ID" = "individual_id")) %>%
  # Remove fish found alive at the end of the experiment
  filter(Found == '0') %>%
  # Flag likely predation events (i.e., when consecutive encounters are more than 2 days)
  mutate(likely_predated = ifelse(consecutive_days > 2, 1, 0))

# Save the summarized potential predation events for future use
saveRDS(muddyfoot_roach_pike_pred_int, paste0(enc_path, "muddyfoot_roach_predation_events.rds"))


# >> 3.1.2. Create summary dataframe of encounters ####

# The purpose of this section is to evaluate whether the number of encounters might be inflated due to the Roach being
# eaten by a Pike but still being tracked (e.g., high encounter rates post-predation). The analysis aims to determine
# when such high encounter rates started and remove data after potential predation.

# >>>> Unfiltered data ####

# This creates a dataframe without filtering out potential predation events (i.e., no post-predation data is removed).

# Identify the date with the highest number of encounters for each Roach and Pike pair
max_encounter_dates <- roach_pike_encounter_sum %>%
  group_by(Roach_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(Roach_ID, Pike_ID, max_encounter_date = date, max_encounter_count = encounter_count)

# Create a summary dataframe for Roach-Pike encounters without filtering
roach_pike_encounter_sum_unfilt <- roach_pike_encounter_sum %>%
  group_by(Roach_ID, Pike_ID) %>%
  mutate(num_encounters_w_pike_id = sum(encounter_count)) %>%
  group_by(Roach_ID) %>%
  mutate(total_pike_encounters = sum(num_encounters_w_pike_id),
         avg_dist_from_pike = mean(daily_avg_dist, na.rm = TRUE)) %>%
  left_join(max_encounter_dates, by = c("Roach_ID", "Pike_ID"))

# Merge encounter summary with metadata for individual Roach
roach_pike_encounter_sum_unfilt <- roach_pike_encounter_sum_unfilt %>%
  left_join(sel_cols, by = c("Roach_ID" = "individual_id"))

# Add information about whether Roach were found alive or likely predated
roach_pike_encounter_sum_unfilt <- roach_pike_encounter_sum_unfilt %>%
  left_join(post_size_cols, by = c("Roach_ID" = "individual_id"))

# Merge likely predation information
pred_cols <- muddyfoot_roach_pike_pred_int %>%
  filter(likely_predated == '1') %>%
  dplyr::select(Roach_ID, likely_predated)

roach_pike_encounter_sum_unfilt <- roach_pike_encounter_sum_unfilt %>%
  left_join(pred_cols, by = c("Roach_ID" = "Roach_ID")) %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%
  # Flag individuals as likely dead if they were predated or had more than 3 missing tracking days
  mutate(likely_died = ifelse(n_missing_dates > 3 & Found == 0 | likely_predated == 1, 1, 0),
         days_tracked = 36 - n_missing_dates,  # Total tracking duration (max 36 days)
         avg_pred_encounter_rate = round(total_pike_encounters / days_tracked, 0))

# Save the unfiltered encounter summary
saveRDS(roach_pike_encounter_sum_unfilt, paste0(enc_path, 'roach_pike_encounter_sum_unfilt.rds'))

# >>>> Filtered data ####

# Create a dataframe that filters out potential post-predation encounters by removing rows after the first
# date with more than 50 encounters for likely predated Roach.

# Filter dataframe to remove post-predation encounters
pred_cols <- muddyfoot_roach_pike_pred_int %>%
  filter(likely_predated == '1') %>%
  dplyr::select(Roach_ID, likely_predated, first_date_over_50)

roach_pike_encounter_sum_filt <- roach_pike_encounter_sum %>%
  left_join(pred_cols, by = "Roach_ID") %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%
  filter(is.na(first_date_over_50) | date <= first_date_over_50)  # Keep only pre-predation data

# Recalculate the most encounter-heavy date for each Roach-Pike pair after filtering
max_encounter_dates_filt <- roach_pike_encounter_sum_filt %>%
  group_by(Roach_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(Roach_ID, Pike_ID, max_encounter_date = date, max_encounter_count = encounter_count)

# Recreate the summary table with filtered data
roach_pike_encounter_sum_filt <- roach_pike_encounter_sum_filt %>%
  group_by(Roach_ID, Pike_ID) %>%
  mutate(num_encounters_w_pike_id = sum(encounter_count)) %>%
  group_by(Roach_ID) %>%
  mutate(total_pike_encounters = sum(num_encounters_w_pike_id),
         avg_dist_from_pike = mean(daily_avg_dist, na.rm = TRUE)) %>%
  left_join(max_encounter_dates_filt, by = c("Roach_ID", "Pike_ID"))

# Merge filtered encounter summary with metadata for individual Roach
roach_pike_encounter_sum_filt <- roach_pike_encounter_sum_filt %>%
  left_join(sel_cols, by = c("Roach_ID" = "individual_id"))

# Add whether fish were found alive and merge with predation information
roach_pike_encounter_sum_filt <- roach_pike_encounter_sum_filt %>%
  left_join(post_size_cols, by = c("Roach_ID" = "individual_id")) %>%
  left_join(pred_cols, by = c("Roach_ID" = "Roach_ID")) %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%
  # Flag likely deaths based on missing tracking dates or predation
  mutate(likely_died = ifelse(n_missing_dates > 3 & Found == 0 | likely_predated == 1, 1, 0),
         days_tracked = 36 - n_missing_dates)  # Total tracking duration (max 36 days)
         
         
# > 3.2. Perch - Pike encounters #### 

# This section calculates encounters between Perch and Pike, specifically focusing on close interactions
# (within 0.45 meters, which is the max strike distance of a Pike, according to Harper 1991).

# Identify encounters where Perch and Pike are within 0.45 meters
perch_pike_distances_df$encounter <- ifelse(perch_pike_distances_df$est <= 0.45, 1, 0)

# Add a date column based on the timestamp for further daily analysis
perch_pike_distances_df <- perch_pike_distances_df %>%
  mutate(date = as.Date(timestamp))

# Summarize the number of encounters and calculate the average daily distance for each Perch-Pike pair, excluding any outliers
perch_pike_encounter_sum <- perch_pike_distances_df %>%
  group_by(perch_ID, Pike_ID, date) %>%
  summarize(encounter_count = sum(encounter),
            daily_avg_dist = mean(est, na.rm = TRUE)) %>%
  filter(date != '2022-09-24')  # Exclude any known problematic or outlier dates       
         
         
# >> 3.2.1. Identify possible predation events ####

# Filter for days when there were more than 50 encounters or when the average distance between Perch and Pike 
# was less than 3 meters, which may indicate possible predation events.
over_50_encounters <- perch_pike_encounter_sum %>%
  filter(encounter_count >= 50 | daily_avg_dist < 3)

# Summarize the Perch-Pike interactions with more than 50 encounters:
muddyfoot_perch_pike_pred_int <- over_50_encounters %>%
  group_by(perch_ID, Pike_ID) %>%
  summarise(
    num_days_over_50 = n(),  # Count of days with over 50 encounters or avg distance < 3 meters
    encounter_dates = paste(format(date, "%d/%m/%Y"), collapse = ", "),  # Concatenate dates of encounters into a string
    consecutive_days = max(rle(cumsum(c(1, diff(date) > 1)) == 1)$lengths),  # Calculate consecutive encounter days
    first_date_over_50 = min(date)  # Identify first date with more than 50 encounters
  )

# Extract additional metadata such as treatment type and number of missing tracking dates for each Perch individual
sel_cols <- muddyfoot_filt_data %>%
  filter(Species == 'Perch') %>%
  dplyr::select(individual_id, treatment, n_missing_dates) %>%
  distinct()

# Merge the summarized predation events with the individual-level metadata
muddyfoot_perch_pike_pred_int <- muddyfoot_perch_pike_pred_int %>%
  left_join(sel_cols, by = c("perch_ID" = "individual_id"))

# Load post-experiment biometric data to check for known/assumed deaths and survivals during the experiment
post_biometrics <- fread(paste0(size_path, "biometric_post_exp_data.csv")) %>%
  mutate(individual_id = paste0("F", sub(".*-", "", Tag_Number)))

# Add information about whether each Perch was found alive at the end of the experiment
post_size_cols <- post_biometrics %>%
  filter(Lake == 'Muddyfoot', Species == 'Perch') %>%
  dplyr::select(individual_id, Found, Known_predated)

muddyfoot_perch_pike_pred_int <- muddyfoot_perch_pike_pred_int %>%
  left_join(post_size_cols, by = c("perch_ID" = "individual_id")) %>%
  # Filter out fish that were found alive at the end of the experiment
  filter(Found == '0') %>%
  # Flag likely predation events based on consecutive days of encounters or if they were known to be predated
  mutate(likely_predated = ifelse(consecutive_days >= 2 | Known_predated == 1, 1, 0))

# Save the summarized predation events for Perch-Pike encounters
saveRDS(muddyfoot_perch_pike_pred_int, paste0(enc_path, "muddyfoot_perch_pike_pred_int.rds"))        
         

# >>>> Muddyfoot predation events ####

# Combine Roach-Pike and Perch-Pike predation events into a single dataframe for all predation events in Muddyfoot lake.
muddyfoot_pred_events <- rbind(
  muddyfoot_roach_pike_pred_int %>% rename(individual_ID = Roach_ID) %>% mutate(Species = 'Roach'),
  muddyfoot_perch_pike_pred_int %>% rename(individual_ID = perch_ID) %>% mutate(Species = 'Perch')
) %>%
  # Filter events with more than 1 consecutive day of encounters
  filter(consecutive_days > 1) %>%
  dplyr::select(individual_ID, Species, treatment, everything()) %>%
  dplyr::select(-likely_predated)

# Save the combined predation event data
saveRDS(muddyfoot_pred_events, paste0(enc_path, "muddyfoot_pred_events.rds"))


# >> 3.2.2. Create summary dataframe of encounters ####

# The next section investigates whether encounter rates might be inflated if a fish was eaten and still being tracked.
# This section aims to identify when the fish started having very high encounter rates and remove post-predation data.

# >>>> Unfiltered ####

# First, create a dataframe without removing encounters, including potential post-predation data.

# Identify the date with the most encounters for each Perch and Pike pair
max_encounter_dates <- perch_pike_encounter_sum %>%
  group_by(perch_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(perch_ID, Pike_ID, max_encounter_date = date, max_encounter_count = encounter_count)

# Create a summary dataframe with unfiltered data (i.e., all encounters)
perch_pike_encounter_sum_unfilt <- perch_pike_encounter_sum %>%
  group_by(perch_ID, Pike_ID) %>%
  mutate(num_encounters_w_pike_id = sum(encounter_count)) %>%
  group_by(perch_ID) %>%
  mutate(total_pike_encounters = sum(num_encounters_w_pike_id),
         avg_dist_from_pike = mean(daily_avg_dist, na.rm = TRUE)) %>%
  left_join(max_encounter_dates, by = c("perch_ID", "Pike_ID"))

# Merge the encounter summary with metadata about each individual Perch
perch_pike_encounter_sum_unfilt <- perch_pike_encounter_sum_unfilt %>%
  left_join(sel_cols, by = c("perch_ID" = "individual_id"))

# Add information on whether each fish was found alive or likely predated
perch_pike_encounter_sum_unfilt <- perch_pike_encounter_sum_unfilt %>%
  left_join(post_size_cols, by = c("perch_ID" = "individual_id"))

# Merge information about likely predation events
pred_cols <- muddyfoot_perch_pike_pred_int %>%
  filter(likely_predated == '1') %>%
  dplyr::select(perch_ID, likely_predated)

perch_pike_encounter_sum_unfilt <- perch_pike_encounter_sum_unfilt %>%
  left_join(pred_cols, by = c("perch_ID" = "perch_ID")) %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%
  # Flag individuals as likely dead if they were predated or had more than 3 missing tracking days
  mutate(likely_died = ifelse(n_missing_dates > 3 & Found == 0 | likely_predated == 1, 1, 0),
         days_tracked = 36 - n_missing_dates,  # Total tracking duration (max 36 days)
         avg_pred_encounter_rate = round(total_pike_encounters / days_tracked, 0))

# Save the unfiltered encounter summary
saveRDS(perch_pike_encounter_sum_unfilt, paste0(enc_path, 'perch_pike_encounter_sum_unfilt.rds'))


# >>>> Filtered ####

# Now, create a dataframe that removes encounters after likely predation events (i.e., post-predation encounters).

# Filter dataframe to remove post-predation encounters after the first day with >50 encounters
pred_cols <- muddyfoot_perch_pike_pred_int %>%
  filter(consecutive_days > 2) %>%
  dplyr::select(perch_ID, likely_predated, first_date_over_50)

perch_pike_encounter_sum_filt <- perch_pike_encounter_sum %>%
  left_join(pred_cols, by = "perch_ID") %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%
  # Keep only rows before or on the first day with more than 50 encounters
  filter(is.na(first_date_over_50) | date <= first_date_over_50)


# Identify the date with the most encounters for each Perch and Pike pair - filtered dataset
# This section calculates the day with the highest number of encounters between Perch and Pike after removing likely predation events.

max_encounter_dates_filt <- perch_pike_encounter_sum_filt %>%
  group_by(perch_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%  # If there are multiple dates with the same max encounter count, select the first one
  dplyr::select(perch_ID, Pike_ID, max_encounter_date = date, max_encounter_count = encounter_count)

# Summarize filtered encounter data for Perch-Pike interactions
perch_pike_encounter_sum_filt <- perch_pike_encounter_sum_filt %>%
  group_by(perch_ID, Pike_ID) %>%
  mutate(num_encounters_w_pike_id = sum(encounter_count)) %>%  # Total number of encounters between each Perch and Pike
  group_by(perch_ID) %>%
  mutate(total_pike_encounters = sum(num_encounters_w_pike_id),  # Total encounters for each Perch with all Pike
         avg_dist_from_pike = mean(daily_avg_dist, na.rm = TRUE)) %>%  # Average distance between each Perch and Pike
  left_join(max_encounter_dates_filt, by = c("perch_ID", "Pike_ID")) %>%  # Add the date and count of the max encounters
  dplyr::select(-likely_predated, -first_date_over_50)  # Remove columns related to predation filtering

# Add metadata such as treatment type and number of missing tracking dates
perch_pike_encounter_sum_filt <- perch_pike_encounter_sum_filt %>%
  left_join(sel_cols, by = c("perch_ID" = "individual_id"))

# Create a smaller summary table for further analysis
perch_pike_encounter_sum_filt <- perch_pike_encounter_sum_filt %>%
  dplyr::select(perch_ID, treatment, total_pike_encounters, avg_dist_from_pike, n_missing_dates) %>% 
  distinct()  # Remove duplicate rows to ensure one entry per Perch

# Add information about whether the Perch were found alive at the end of the experiment
perch_pike_encounter_sum_filt <- perch_pike_encounter_sum_filt %>%
  left_join(post_size_cols, by = c("perch_ID" = "individual_id"))

# Add information about likely predation events
pred_cols <- muddyfoot_perch_pike_pred_int %>%
  filter(likely_predated == '1') %>%
  dplyr::select(perch_ID, likely_predated)

# Merge predation information into the summary table
perch_pike_encounter_sum_filt <- perch_pike_encounter_sum_filt %>%
  left_join(pred_cols, by = c("perch_ID" = "perch_ID")) %>%
  tidyr::replace_na(list(likely_predated = 0)) %>%  # Replace NA values with 0 for non-predated fish
  mutate(likely_died = ifelse(n_missing_dates > 3 & Found == 0 | likely_predated == 1, 1, 0),  # Determine if fish likely died
         days_tracked = 36 - n_missing_dates,  # Calculate total days tracked (36 days max)
         avg_pred_encounter_rate = round(total_pike_encounters / days_tracked, 0)) %>%  # Average daily encounter rate with Pike
  group_by(perch_ID) %>%
  distinct()  # Remove duplicate rows

# Save the filtered summary table for Perch-Pike encounters
saveRDS(perch_pike_encounter_sum_filt, paste0(enc_path, 'perch_pike_encounter_sum_filt.rds'))


# >>>> Muddyfoot encounter summary ####

# Combine Roach-Pike and Perch-Pike encounter data into a single summary table of predation events in Muddyfoot lake.
muddyfoot_pred_encounter_summary <- rbind(
  roach_pike_encounter_sum_filt %>% rename(individual_ID = Roach_ID) %>% mutate(Species = 'Roach'),
  perch_pike_encounter_sum_filt %>% rename(individual_ID = perch_ID) %>% mutate(Species = 'Perch')
) %>%
  dplyr::select(individual_ID, Species, treatment, everything()) %>%
  tidyr::replace_na(list(Known_predated = 0))  # Replace NA values in Known_predated column

# Save the combined predation event summary
saveRDS(muddyfoot_pred_encounter_summary, paste0(enc_path, "muddyfoot_pred_encounter_summary.rds"))

#----------------------------------------------------------------------------------------------------------#

#---------------------------------------------------#
######### 4. Visualise predator encounters #########
#--------------------------------------------------#


#To be complete

### Extract individuals with irregular tracking ###
irreg_ind <- 
  muddyfoot_dat %>% 
  filter(irregular_sampling == 1)

#how many individuals
length(unique(irreg_ind$individual_id))
#19

#Explore F59707
#It was only tracked for 1 day
F59707_track = 
  irreg_ind %>% 
  filter(individual_id == 'F59707')


str(F59707_track)
class(F59707_track)

# now create a move2 object from a data.frame
F59707_track_mv <- 
  mt_as_move2(F59707_track, 
              coords = c("longitude","latitude"),
              crs = "WGS84",
              time_column = "timestamp",
              track_id_column = "individual_id",
              na.fail = F) # allows or not empty coordinates


#Plot tracks over outline of lake
#Load map polygon
muddyfoot_polygon = sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))

ggplot() + 
  coord_sf(crs = st_crs(F59707_track_mv)) + #set the coordinate system to that of the bats data
  geom_sf(data=muddyfoot_polygon) + 
  geom_sf(data=F59707_track_mv, col="red", alpha = 0.5)


#add scale bar
library(sf)
library(ggsn)
library(gganimate)

# Check the CRS of your lake polygon
st_crs(muddyfoot_polygon)

# Define the length of the scale bar in meters
scale_length <- 100  # 1 km

# Get the bounding box of the lake to determine where to place the scale bar
lake_bbox <- st_bbox(muddyfoot_polygon)

# Define the starting point of the scale bar near the bottom left corner of the map
x_start <- lake_bbox["xmin"] + (lake_bbox["xmax"] - lake_bbox["xmin"]) * 0.1
y_start <- lake_bbox["ymin"] + (lake_bbox["ymax"] - lake_bbox["ymin"]) * 0.1

# Define the ending point of the scale bar
x_end <- x_start + scale_length

p = ggplot() + 
  geom_sf(data=muddyfoot_polygon) + 
  geom_sf(data=F59707_track_mv, col="red", alpha = 0.5) + 
  ggsn::scalebar(
    data = muddyfoot_polygon, 
    dist = 5,  # Distance represented by the scale bar
    dist_unit = "m",  # Units of measurement
    transform = TRUE,  # Transform coordinates to projected CRS if needed
    model = 'WGS84',  # Model used for transformation, if needed
    location = "bottomright",  # Position of the scale bar
    st.dist = 0.05  # Distance of the text from the scale bar
  ) +
  theme_minimal() +
  transition_time(F59707_track_mv$timestamp) +  # Animate over the timestamp
  labs(title = "Fish Movement Over Time: {frame_time}") +  # Title that updates with time
  ease_aes('linear')

animate(p, nframes = 500, fps = 10, width = 800, height = 600)


#Later develop interactive plot to illustrate potential predation event. Plot predators for this day on top.  
#Also want to plot shelters on to map

#Extract pike data
pike_dat <- 
  muddyfoot_dat %>% 
  filter(Species == 'Pike' & Date == '2022-09-25')

#combine it with individual roach data
encounter_test_data = rbind(pike_dat, F59707_track)

#create move object
pred_encounter_track_mv <- 
  mt_as_move2(encounter_test_data, 
              coords = c("longitude","latitude"),
              crs = "WGS84",
              time_column = "timestamp",
              track_id_column = "individual_id",
              na.fail = F) # allows or not empty coordinates


test = encounter_test_data %>% 
  filter(individual_id == 'F59707' | individual_id == 'F59880')

p = ggplot() + 
  geom_sf(data = muddyfoot_polygon, fill = NA, color = "blue") +
  geom_point(data = test, 
             aes(x = longitude, 
                 y = latitude, 
                 group = individual_id,
                 color = individual_id), 
             size = 4, alpha = 0.3) +  # Plot the fish tracks with color by fish_id
  # geom_path(data = encounter_test_data, 
  #           aes(x = longitude, 
  #               y = latitude, 
  #               group = individual_id, 
  #               color = individual_id), 
  #           size = 1, 
  #           alpha = 0.5) +  # Add paths
  ggsn::scalebar(
    data = muddyfoot_polygon, 
    dist = 5,  # Distance represented by the scale bar
    dist_unit = "m",  # Units of measurement
    transform = TRUE,  # Transform coordinates to projected CRS if needed
    model = 'WGS84',  # Model used for transformation, if needed
    location = "bottomright",  # Position of the scale bar
    st.dist = 0.05  # Distance of the text from the scale bar
  ) +
  theme_minimal() +  
  labs(title = "Fish Movement Over Time: {frame_time}") +## Title that updates with time
  transition_reveal(along = timestamp) +  
  ease_aes('linear')

animate(p, nframes = 500, fps = 10, width = 800, height = 600)

combined_telemetry






#








### Extract individuals with irregular tracking ###
irreg_ind <- 
  muddyfoot_dat %>% 
  filter(irregular_sampling == 1)

#how many individuals
length(unique(irreg_ind$individual_id))
#19

#Explore F59707
#It was only tracked for 1 day
F59707_track = 
  irreg_ind %>% 
  filter(individual_id == 'F59707')


str(F59707_track)
class(F59707_track)

# now create a move2 object from a data.frame
F59707_track_mv <- 
  mt_as_move2(F59707_track, 
              coords = c("longitude","latitude"),
              crs = "WGS84",
              time_column = "timestamp",
              track_id_column = "individual_id",
              na.fail = F) # allows or not empty coordinates


#Plot tracks over outline of lake
#Load map polygon
muddyfoot_polygon = sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))

ggplot() + 
  coord_sf(crs = st_crs(F59707_track_mv)) + #set the coordinate system to that of the bats data
  geom_sf(data=muddyfoot_polygon) + 
  geom_sf(data=F59707_track_mv, col="red", alpha = 0.5)


#add scale bar
library(sf)
library(ggsn)
library(gganimate)

# Check the CRS of your lake polygon
st_crs(muddyfoot_polygon)

# Define the length of the scale bar in meters
scale_length <- 100  # 1 km

# Get the bounding box of the lake to determine where to place the scale bar
lake_bbox <- st_bbox(muddyfoot_polygon)

# Define the starting point of the scale bar near the bottom left corner of the map
x_start <- lake_bbox["xmin"] + (lake_bbox["xmax"] - lake_bbox["xmin"]) * 0.1
y_start <- lake_bbox["ymin"] + (lake_bbox["ymax"] - lake_bbox["ymin"]) * 0.1

# Define the ending point of the scale bar
x_end <- x_start + scale_length

p = ggplot() + 
  geom_sf(data=muddyfoot_polygon) + 
  geom_sf(data=F59707_track_mv, col="red", alpha = 0.5) + 
  ggsn::scalebar(
    data = muddyfoot_polygon, 
    dist = 5,  # Distance represented by the scale bar
    dist_unit = "m",  # Units of measurement
    transform = TRUE,  # Transform coordinates to projected CRS if needed
    model = 'WGS84',  # Model used for transformation, if needed
    location = "bottomright",  # Position of the scale bar
    st.dist = 0.05  # Distance of the text from the scale bar
  ) +
  theme_minimal() +
  transition_time(F59707_track_mv$timestamp) +  # Animate over the timestamp
  labs(title = "Fish Movement Over Time: {frame_time}") +  # Title that updates with time
  ease_aes('linear')

animate(p, nframes = 500, fps = 10, width = 800, height = 600)


#Later develop interactive plot to illustrate potential predation event. Plot predators for this day on top.  
#Also want to plot shelters on to map

#Extract pike data
pike_dat <- 
  muddyfoot_dat %>% 
  filter(Species == 'Pike' & Date == '2022-09-25')

#combine it with individual roach data
encounter_test_data = rbind(pike_dat, F59707_track)

#create move object
pred_encounter_track_mv <- 
  mt_as_move2(encounter_test_data, 
              coords = c("longitude","latitude"),
              crs = "WGS84",
              time_column = "timestamp",
              track_id_column = "individual_id",
              na.fail = F) # allows or not empty coordinates


test = encounter_test_data %>% 
  filter(individual_id == 'F59707' | individual_id == 'F59880')

p = ggplot() + 
  geom_sf(data = muddyfoot_polygon, fill = NA, color = "blue") +
  geom_point(data = test, 
             aes(x = longitude, 
                 y = latitude, 
                 group = individual_id,
                 color = individual_id), 
             size = 4, alpha = 0.3) +  # Plot the fish tracks with color by fish_id
  # geom_path(data = encounter_test_data, 
  #           aes(x = longitude, 
  #               y = latitude, 
  #               group = individual_id, 
  #               color = individual_id), 
  #           size = 1, 
  #           alpha = 0.5) +  # Add paths
  ggsn::scalebar(
    data = muddyfoot_polygon, 
    dist = 5,  # Distance represented by the scale bar
    dist_unit = "m",  # Units of measurement
    transform = TRUE,  # Transform coordinates to projected CRS if needed
    model = 'WGS84',  # Model used for transformation, if needed
    location = "bottomright",  # Position of the scale bar
    st.dist = 0.05  # Distance of the text from the scale bar
  ) +
  theme_minimal() +  
  labs(title = "Fish Movement Over Time: {frame_time}") +## Title that updates with time
  transition_reveal(along = timestamp) +  
  ease_aes('linear')

animate(p, nframes = 500, fps = 10, width = 800, height = 600)

combined_telemetry
