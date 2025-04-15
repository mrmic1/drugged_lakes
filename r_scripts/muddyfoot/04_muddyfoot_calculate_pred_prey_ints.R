# ------------------------------------------------# 
# Exploring predator-prey interactions - Muddyfoot
# ------------------------------------------------#

### SCRIPT DESCRIPTION ###

# The goal of this script is to identify potential predation events based on close proximity and track alignment and to create 
# a dataframe that can be used for analysing potential treatment effects on predator encounter rates

#------------------------------------------------------------------------------------------------------------#

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
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
lake_polygon_path <- "./data/lake_coords/"
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/muddyfoot/"
enc_path <- "./data/encounters/muddyfoot/"
size_path <- "./data/fish_size/"
save_tables_path <- "./tables/muddyfoot/"  # Path to save summary tables.


### DATA LOADING ###
muddyfoot_filt_data <-  readRDS(paste0(filtered_data_path, "03_muddyfoot_sub.rds"))

#--------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------#
# 1. Calculate distances at shared timepoints for prey and predators #######
#---------------------------------------------------------------------------------#

#This part of the script takes a long time to run
#Outcomes are stored in 
# roach_pike_distances_df <- readRDS(paste0(enc_path, "muddyfoot_pike_roach_distances_df.rds"))
# perch_pike_distances_df <- readRDS(paste0(enc_path, "muddyfoot_perch_distances_df.rds"))


pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds"))

pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

# This section calculates pairwise distances between Roach/Perch (prey) and Pike (predator) at shared time points.

#### > 1.1. Roach ####

# Reproject Roach telemetry objects to match the Pike projection
ctmm::projection(roach_muddyfoot_tel) <- ctmm::projection(pike_muddyfoot_tel)
ctmm::projection(roach_muddyfoot_ctmm_fits) <- ctmm::projection(pike_muddyfoot_ctmm_fits)

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


#### > 1.2. Perch ####

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
# saveRDS(perch_pike_distances_df, paste0(enc_path, "muddyfoot_pike_perch_distances_df.rds"))

#---------------------------------------------------#
# 2. Calculate empirical encounters #########
#--------------------------------------------------#

roach_pike_distances_df <- readRDS(paste0(enc_path, "muddyfoot_pike_roach_distances_df.rds"))
perch_pike_distances_df <- readRDS(paste0(enc_path, "muddyfoot_pike_perch_distances_df.rds"))

#--------------------------------------------------------------------------------------------#

# > 2.1. Total encounters #### 

# This section calculates encounters between prey and pike (predator), specifically focusing on close interactions
# (within 0.45 meters, which is the max strike distance of a Pike, according to Harper 1991).

# Define encounters as those occurring within 0.45 meters (based on Pike's strike distance)
# roach_pike_distances_df$encounter <- ifelse(roach_pike_distances_df$est <= 0.45, 1, 0)
# perch_pike_distances_df$encounter <- ifelse(perch_pike_distances_df$est <= 0.45, 1, 0)

# Add a date column based on the timestamp for further daily analysis
# Note takes time to run - do not run if you don't have to
# time column has been saved in the provided files

# perch_pike_distances_df <- 
#    perch_pike_distances_df %>%
#    mutate(Date = as.Date(format(perch_pike_distances_df$timestamp, tz = "Europe/Stockholm")))
# 
# roach_pike_distances_df <-
#   roach_pike_distances_df %>%
#   mutate(Date = as.Date(format(roach_pike_distances_df$timestamp, tz = "Europe/Stockholm")))

# Summarise the number of encounters and calculate the 
# average daily distance for each perch-pike pair,
# and minimum daily distance

roach_pike_encounter_sum <- 
  roach_pike_distances_df %>%
  group_by(Roach_ID, Pike_ID, Date) %>%
  summarise(encounter_count = sum(encounter),
            daily_avg_dist_from_pike = mean(est, na.rm = TRUE),
            daily_min_dist_from_pike = min(est, na.rm = TRUE))
#5286 rows

perch_pike_encounter_sum <- 
  perch_pike_distances_df %>%
  group_by(perch_ID, Pike_ID, Date) %>%
  summarise(encounter_count = sum(encounter),
            daily_avg_dist_from_pike = mean(est, na.rm = TRUE),
            daily_min_dist_from_pike = min(est, na.rm = TRUE))
#5500 rows

# #save roach and perch predator encounter summaries
saveRDS(roach_pike_encounter_sum, paste0(enc_path, "muddyfoot_roach_pike_encounter_summary_unfiltered.rds"))
saveRDS(perch_pike_encounter_sum, paste0(enc_path, "muddyfoot_perch_pike_encounter_summary_unfiltered.rds"))


# Renaming the columns to have a common name before binding
roach_pike_encounter_sum <- roach_pike_encounter_sum %>%
  rename(Prey_ID = Roach_ID)

perch_pike_encounter_sum <- perch_pike_encounter_sum %>%
  rename(Prey_ID = perch_ID)

# Combine the data frames into a single data frame
prey_pike_encounter_sum <- rbind(cbind(roach_pike_encounter_sum, Species = 'Roach'), 
                                 cbind(perch_pike_encounter_sum, Species = 'Perch'))
#10836 rows

# Check the structure of the combined data frame
str(prey_pike_encounter_sum)

#save prey predator encounter summary
saveRDS(prey_pike_encounter_sum, paste0(enc_path, "muddyfoot_prey_pike_encounter_summary_unfiltered.rds"))

#----------------------------------------------------------------------------------------------------#

#> 2.2. Suspected predation events ####

# This section identifies potential predation events by filtering for days when pike and prey
# had more than 25 encounters or when their average daily distance was less than 0.2 meters.

post_biometrics <- 
  fread(paste0(size_path, "biometric_post_exp_data.csv")) %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number)))

#Summarise number of survivors and known death for each species in muddyfoot
post_biometrics %>% 
  filter(Lake == 'Muddyfoot') %>% 
  dplyr::select(individual_ID, Species, Treatment, Found) %>% 
  group_by(Species, Treatment) %>% 
  summarise(num_survive = sum(Found),
            num_died = length(unique(individual_ID)) - sum(Found),
            percent_survive = sum(Found)/length(unique(individual_ID)))


# Merge biometric data (e.g., whether the fish was found alive) with encounter summary
prey_post_size_cols <- 
  post_biometrics %>%
  filter(Lake == 'Muddyfoot', Species == 'Roach'| Species == 'Perch') %>%
  dplyr::select(individual_ID, Found, Known_predated) %>% 
  rename(found_alive = Found,
         found_predated = Known_predated)

prey_pike_encounter_sum <- 
  prey_pike_encounter_sum %>%
  left_join(prey_post_size_cols, by = c("Prey_ID" = "individual_ID"))

# Now we filter prey-pike encounters for days with more than 25 encounters or an avg daily distance < 0.2 meters.
# This our threshold for identifying potential predation events.
# This threshold was arbitrarily chosen, but seems conservative 

prey_over_25_encounters <- 
  prey_pike_encounter_sum %>%
  filter(encounter_count >= 25  | daily_min_dist_from_pike < 0.2)

#how many unique roach individuals
unique(prey_over_25_encounters$Prey_ID)
#42 unique individuals

# Summarise prey-pike interactions: count days with > 25 encounters, consecutive days, and the first date
muddyfoot_prey_pike_interactions <- 
  prey_over_25_encounters %>%
  group_by(Prey_ID, Pike_ID, Species) %>%
  summarise(
    # Count of days with at least 25 encounters
    num_days_25_encounters = sum(encounter_count >= 25),
    # Count of days with at least 100 encounters
    num_days_100_encounters = sum(encounter_count >= 100),
    # Count of days where min dist from pike was less than 0.2 m
    num_days_min_dist_less_0.2m = sum(daily_min_dist_from_pike < 0.2),
    # Combine encounter dates as a string
    encounter_dates = paste(format(Date, "%d/%m/%Y"), collapse = ", "),
    # Identify first date with 25+ encounters
    first_date_over_25 = if (any(encounter_count >= 25)) min(Date[encounter_count >= 25]) else as.Date(NA),
    # Identify first date with 100+ encounters
    first_date_over_100 = if (any(encounter_count >= 100)) min(Date[encounter_count >= 100]) else as.Date(NA),
    # Calculate max consecutive days with 25+ encounters
    consecutive_days_25 = {
      dates_25 <- Date[encounter_count >= 25]
      if (length(dates_25) == 0) 0 else max(rle(c(1, diff(sort(dates_25)) == 1))$lengths)
    },
    # Calculate max consecutive days with 100+ encounters
    consecutive_days_100 = {
      dates_100 <- Date[encounter_count >= 100]
      if (length(dates_100) == 0) 0 else max(rle(c(1, diff(sort(dates_100)) == 1))$lengths)
    }
  )
#69 rows

# Extract individual treatment and missing dates information from filtered data
prey_sel_cols <- 
  muddyfoot_filt_data %>%
  filter(Species == 'Roach'| Species == 'Perch') %>%
  dplyr::select(individual_ID, Treatment, n_missing_dates) %>%
  distinct()

# Merge encounter summary with individual-level metadata
muddyfoot_prey_pike_interactions <- 
  muddyfoot_prey_pike_interactions %>%
  left_join(prey_sel_cols, by = c("Prey_ID" = "individual_ID"))

# Merge bio metric data (e.g., whether the fish was found alive) with possible predator events
# Create column indicating instances where there might be a higher likelihood that an individual was predated
# specifically, when individuals 2 or more consecutive days where they had more than 25 encounters with a predator
muddyfoot_prey_pike_interactions <- 
  muddyfoot_prey_pike_interactions %>%
  left_join(prey_post_size_cols, by = c("Prey_ID" = "individual_ID")) %>%
  # Remove fish found alive at the end of the experiment
  filter(found_alive == '0') %>%
  # Flag likely predation events (i.e., when consecutive days of over 100 encounters are more than 2 days)
  mutate(likely_predated = ifelse(consecutive_days_100 >= 2, 1, 0))

# Identify the date with the most encounters for each roach and pike pair
muddyfoot_prey_max_encounter_dates <- 
  prey_pike_encounter_sum %>%
  group_by(Prey_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(Prey_ID, Pike_ID, Species, max_encounter_date = Date, max_encounter_count = encounter_count)

#Add max_encounter_date and max_encounter_count to potential predation events summary.
muddyfoot_prey_pike_interactions <- 
  muddyfoot_prey_pike_interactions %>%
  left_join(muddyfoot_prey_max_encounter_dates %>%
              dplyr::select(Prey_ID, Pike_ID, max_encounter_date, max_encounter_count),
            by = c("Prey_ID" = "Prey_ID", 
                   "Pike_ID" = "Pike_ID"))

# Save the summarised potential predation events for future use
saveRDS(muddyfoot_prey_max_encounter_dates, paste0(enc_path, "muddyfoot_prey_daily_max_pred_encounters.rds"))

# Save the suspected predation event data
saveRDS(muddyfoot_prey_pike_interactions, paste0(enc_path, "muddyfoot_suspected_predation_events.rds"))

#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------#
# 3. Encounter rate summaries - unfiltered #######################
#----------------------------------------------------------------#

# The purpose of this section is to evaluate whether the number of encounters might be inflated due prey being
# eaten by a pike but still being tracked (e.g., high encounter rates post-predation). The analysis aims to determine
# when such high encounter rates started and remove data after potential predation.

# > 3.1. Daily and total encounter summary tables ####

# Create a summary dataframe for the daily number of predator encounters each prey had 
# This will be used to help try and identify the fate of individuals that were not found at the end of the experiment
# see below

muddyfoot_prey_daily_encounter_summary <-
  prey_pike_encounter_sum %>%
  group_by(Prey_ID, Date) %>%
  summarise(Species = first(Species),
            daily_avg_dist_from_pike = mean(daily_avg_dist_from_pike),
            total_daily_encounter_count = sum(encounter_count),
            max_daily_encounter_count = max(encounter_count),
            pike_id_max_encounter_count = Pike_ID[which.max(encounter_count)],
            daily_min_dist_from_pike = min(daily_min_dist_from_pike),
            pike_id_min_dist = Pike_ID[which.min(daily_min_dist_from_pike)],
            found_alive = first(found_alive),
            found_predated = first(found_predated))


#Total experiment encounter summary

muddyfoot_prey_total_encounter_summary <-
  muddyfoot_prey_daily_encounter_summary %>%
  group_by(Prey_ID) %>%
  summarise(Species = first(Species),
            total_encounters = sum(total_daily_encounter_count),
            min_daily_dist_from_pike = min(daily_min_dist_from_pike),
            avg_daily_min_dist_from_pike = min(daily_avg_dist_from_pike, na.rm = TRUE))


# Merge total encounter summary with metadata for each individual.
muddyfoot_prey_total_encounter_summary <-
  muddyfoot_prey_total_encounter_summary %>%
  left_join(prey_sel_cols, by = c("Prey_ID" = "individual_ID"))

#also add the number of poor tracking days info
#the loaded dataframe was created in script 03_muddyfoot_tracking_data_summary
muddyfoot_ID_irreg_sampling <- readRDS(paste0(filtered_data_path, "muddyfoot_ID_irreg_sampling.rds"))
prey_irreg_sampling_cols <- 
  muddyfoot_ID_irreg_sampling %>% 
  filter(Species == 'Roach'| Species == 'Perch') %>% 
  dplyr::select(individual_ID, n_poor_tracking_days, found_alive)

muddyfoot_prey_total_encounter_summary <-
  muddyfoot_prey_total_encounter_summary %>%
  left_join(prey_irreg_sampling_cols, by = c("Prey_ID" = "individual_ID"))

muddyfoot_prey_total_encounter_summary <-
  muddyfoot_prey_total_encounter_summary %>% 
  mutate(days_tracked = 36 - n_missing_dates,  # Total tracking duration (max 34 days)
         avg_daily_pred_encounter_rate = round(total_encounters / days_tracked, 0))

head(muddyfoot_prey_total_encounter_summary)


# Save the unfiltered encounter summary
saveRDS(muddyfoot_prey_total_encounter_summary, paste0(enc_path, 'muddyfoot_prey_total_encounter_summary_unfiltered.rds'))
saveRDS(muddyfoot_prey_daily_encounter_summary, paste0(enc_path, 'muddyfoot_prey_daily_encounter_summary_unfiltered.rds'))

#------------------------------------------------------------------------------------------------------#


# > 3.3. Daily and total encounter summary tables for prey not found ####

#Reduce summary data to only include individuals that were not found at the end of the experiment
muddyfoot_daily_encounter_summary_prey_not_found <- 
  muddyfoot_prey_daily_encounter_summary %>% 
  filter(found_alive == 0)
#381 rows


#Add column to identify whether the date was irregularly sampled or dates were not tracked
#First I need to extract the required columns
prey_not_found_cols <- 
  muddyfoot_ID_irreg_sampling %>% 
  filter(Species == 'Roach'| Species == 'Perch') %>% 
  dplyr::select(individual_ID, dates_irregular_positions, missing_dates)

#Then we need to do some data wrangling to extract the date information
muddyfoot_daily_encounter_summary_prey_not_found$Date <- as.Date(muddyfoot_daily_encounter_summary_prey_not_found$Date)

# Ensure dates_irregular_positions and missing_dates are properly formatted as lists of dates
prey_not_found_cols$dates_irregular_positions <- lapply(strsplit(prey_not_found_cols$dates_irregular_positions, ", "), as.Date)
prey_not_found_cols$missing_dates <- lapply(strsplit(prey_not_found_cols$missing_dates, ", "), as.Date)

# Merge the two datasets on the IDs
merged_df <- merge(muddyfoot_daily_encounter_summary_prey_not_found, prey_not_found_cols, 
                   by.x = "Prey_ID", by.y = "individual_ID", all.x = TRUE)

# Create the two new columns
merged_df$poor_tracking_date <- mapply(function(date, irregular_dates) {
  if (!is.null(irregular_dates) && date %in% irregular_dates) 1 else 0
}, merged_df$Date, merged_df$dates_irregular_positions)

merged_df$no_tracking_date <- mapply(function(date, missing_dates) {
  if (!is.null(missing_dates) && date %in% missing_dates) 1 else 0
}, merged_df$Date, merged_df$missing_dates)

# Keep only the original columns from daily_encounter_sum_roach_found and the new ones
muddyfoot_daily_encounter_summary_prey_not_found <- merged_df[, c("Prey_ID", "Date", "Species", "daily_avg_dist_from_pike","total_daily_encounter_count", 
                                                           "max_daily_encounter_count", "pike_id_max_encounter_count",
                                                           "daily_min_dist_from_pike", "pike_id_min_dist", "found_predated",
                                                           "poor_tracking_date", "no_tracking_date")]

# Save the unfiltered encounter summary
saveRDS(muddyfoot_daily_encounter_summary_prey_not_found, paste0(enc_path, 'muddyfoot_daily_encounter_summaries_for_prey_not_found_unfiltered.rds'))

length(unique(muddyfoot_daily_encounter_summary_prey_not_found$Prey_ID))
#19

#Total encounter summary for prey not found
#This summarised dataframe will include information on
#total encounter count
#max encounter count
#date with max encounter count
#number of days with over 50 encounters
#first date with over 50 encounters

muddyfoot_total_encounter_summary_prey_not_found <-
  muddyfoot_daily_encounter_summary_prey_not_found %>% 
  group_by(Prey_ID) %>% 
  summarise(Species = first(Species),
            total_encounter_count = sum(total_daily_encounter_count),
            max_encounter_count = max(max_daily_encounter_count ),
            max_encounter_date = Date[which.max(max_daily_encounter_count )],
            max_encounter_pike = pike_id_max_encounter_count[which.max(max_daily_encounter_count)],
            min_dist_from_pike = min(daily_min_dist_from_pike),
            min_dist_date = Date[which.min(daily_min_dist_from_pike)],
            min_dist_pike = pike_id_min_dist[which.min(daily_min_dist_from_pike)])

prey_not_found_cols <-
  muddyfoot_prey_pike_interactions %>%
  dplyr::select(Prey_ID, num_days_25_encounters, num_days_100_encounters, num_days_min_dist_less_0.2m,  
                first_date_over_25,consecutive_days_25, first_date_over_100, consecutive_days_100) %>% 
  group_by(Prey_ID) %>% 
  arrange(desc(num_days_100_encounters),
          desc(num_days_100_encounters),
          .by_group = TRUE) %>% 
  slice(1) %>% 
  ungroup()

muddyfoot_total_encounter_summary_prey_not_found <-
  muddyfoot_total_encounter_summary_prey_not_found %>%
  left_join(prey_not_found_cols, by = c("Prey_ID" = "Prey_ID")) %>% 
  dplyr::select(-Pike_ID)


#merge with irregular sampling summary
#create in script 03_muddyfoot_tracking_data_summary
muddyfoot_ID_irreg_sampling_not_found <- 
  muddyfoot_ID_irreg_sampling %>% 
  filter(found_alive == 0)

#get required rows from filtered encounters
filtered_encounters_rows <- 
  muddyfoot_ID_irreg_sampling_not_found %>% 
  dplyr::select(individual_ID, n_poor_tracking_days, total_missing_days)

muddyfoot_total_encounter_summary_prey_not_found <-
  muddyfoot_total_encounter_summary_prey_not_found %>%
  left_join(filtered_encounters_rows, by = c("Prey_ID" = "individual_ID")) %>% 
  dplyr::select(Prey_ID, Species, everything())

# Create a flextable summarizing the irregular and missing days for missing individuals
muddyfoot_daily_irregular_sample_table <- 
  flextable(muddyfoot_total_encounter_summary_prey_not_found) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("Prey_ID" = 'ID',
                    "Species" = 'Species',
                    "total_encounter_count" = "Total encounters",
                    "max_encounter_count" = "Max encounters in day",
                    "max_encounter_date" = "Date w/max encounters",
                    "max_encounter_pike" = "Pike ID w/max encounters",
                    "min_dist_from_pike" = "Minimium recorded dist from pike",
                    "min_dist_date" = "Date w/min  dist",
                    "min_dist_pike" = "Pike ID w/min dist",
                    "num_days_25_encounters" = "Days >= 25 encounters", 
                    "num_days_100_encounters" = "Days >= 100 encounters",
                    "num_days_min_dist_less_0.2m" = "Days < 0.2m distance from pike",
                    "first_date_over_25" = "First date >= 25",
                    "consecutive_days_25" = "consecutive days > 25",
                    "first_date_over_100" = "First date >= 100",
                    "consecutive_days_100" = "consecutive days > 100",
                    "n_poor_tracking_days" = "N poor tracking days",
                    "total_missing_days" = "N days not tracked")

save_as_docx(muddyfoot_daily_irregular_sample_table, 
             path = paste0(save_tables_path, "muddyfoot IDs not found encounter and tracking summary.docx"))

saveRDS(muddyfoot_total_encounter_summary_prey_not_found, paste0(enc_path, 'muddyfoot_total_encounter_summary_for_prey_not_found_unfiltered.rds'))

# Load the openxlsx package
if (!require(openxlsx)) install.packages("openxlsx")
library(openxlsx)

# Save the dataframe as an Excel file
write.xlsx(muddyfoot_total_encounter_summary_prey_not_found, file = paste0(save_tables_path, "muddyfoot_IDs_not_found_total_encounter_summary.xlsx"))

#create seperate table for dates

muddyfoot_dates_table <- 
  muddyfoot_ID_irreg_sampling_not_found %>% 
  dplyr::select(individual_ID, Species, n_poor_tracking_days, dates_irregular_positions, 
                total_missing_days, missing_dates) %>% 
  flextable() %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_ID" = 'ID',
                    "Species" = 'Species',
                    "n_poor_tracking_days" = "N poor tracking days",
                    "dates_irregular_positions" = "Dates poor tracking",
                    "total_missing_days" = "N days not tracked",
                    "missing_dates" = 'Dates not tracked')

save_as_docx(muddyfoot_dates_table, 
             path = paste0(save_tables_path, "muddyfoot_IDs_not_found_tracking_dates.docx"))

#--------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------#

#------------------------------#
# 4. Filter encounter data #### 
#------------------------------#

prey_encounters_sum_IDs_not_found <- readRDS(paste0(enc_path, "muddyfoot_total_encounter_summary_for_prey_not_found_unfiltered.rds"))
mortality_preds <- readxl::read_excel("./data/encounters/suspected_mortality_updated.xlsx")
prey_encounter_sum <- readRDS(paste0(enc_path, "muddyfoot_prey_pike_encounter_summary_unfiltered.rds"))

#-------------------------------------------------------------------------------------------------------------------#

#> 4.1. Filter out encounters for dead prey ####

#Extracted individuals that we believe to be predated
muddyfoot_pred_prey_cols <-
  mortality_preds %>%
  filter(lake == 'muddyfoot') %>% 
  filter(species == 'Roach'| species == 'Perch') %>%
  dplyr::select(individual_ID, revised_suspected_mortality, revised_likely_death_date) %>% 
  rename(death_date = revised_likely_death_date)

prey_encounter_sum_filtered <-
  prey_encounter_sum %>%
  left_join(muddyfoot_pred_prey_cols, c("Prey_ID" = "individual_ID")) %>% 
  filter(is.na(death_date)| Date <= death_date)  # Keep only pre-predation data
#original number of rows: 10836
#new: 10092

# Recalculate the most encounter-heavy date 
max_encounter_dates_filtered <-
  prey_encounter_sum_filtered %>%
  group_by(Prey_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(Prey_ID, Pike_ID, max_encounter_date = Date, max_encounter_count = encounter_count)

# Recreate the summary table with filtered data
prey_encounter_sum_filtered <-
  prey_encounter_sum_filtered %>%
  group_by(Prey_ID, Pike_ID) %>%
  mutate(num_encounters_w_pike_id = sum(encounter_count),
         avg_dist_from_pike = mean(daily_avg_dist_from_pike, na.rm = TRUE),
         min_dist_from_pike = min(daily_min_dist_from_pike, na.rm = TRUE)) %>%
  group_by(Prey_ID) %>%
  mutate(total_encounters = sum(encounter_count),
         avg_dist_from_pred = mean(daily_avg_dist_from_pike, na.rm = TRUE),
         min_dist_from_pred = min(daily_min_dist_from_pike, na.rm = TRUE)) %>%
  left_join(max_encounter_dates_filtered, by = c("Prey_ID", "Pike_ID"))


# Extract individual treatment and missing date information from filtered data
prey_sel_cols <- 
  muddyfoot_filt_data %>%
  filter(Species == 'Roach'| Species == 'Perch') %>%
  dplyr::select(individual_ID, Treatment, n_missing_dates) %>%
  distinct()

# Merge filtered encounter summary with metadata for prey
prey_encounter_sum_filtered <-
  prey_encounter_sum_filtered %>%
  left_join(prey_sel_cols, by = c("Prey_ID" = "individual_ID"))

# Create a smaller summary table for further analysis
prey_total_encounter_summary_filtered <- 
  prey_encounter_sum_filtered %>%
  dplyr::select(Prey_ID, Treatment, total_encounters, avg_dist_from_pred, min_dist_from_pred, 
                n_missing_dates, revised_suspected_mortality, death_date) %>% 
  distinct() 

# Add information about whether each roach was found alive at the end of the experiment
post_size_cols <- 
  post_biometrics %>%
  filter(Lake == 'muddyfoot', Species == 'Roach'| Species == 'Perch') %>%
  dplyr::select(individual_ID, Found, Known_predated) %>% 
  rename(found_alive = Found,
         found_predated = Known_predated)

# Add information about whether the Perch were found alive at the end of the experiment
prey_total_encounter_summary_filtered <- 
  prey_total_encounter_summary_filtered %>%
  left_join(post_size_cols, by = c("Prey_ID" = "individual_ID"))

# Save the combined predation event summary
saveRDS(prey_total_encounter_summary_filtered, paste0(enc_path, "muddyfoot_pred_prey_encounter_summary_filtered.rds"))

#----------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------#

#---------------------------------------------------#
######### 6. Visualise predator encounters #########
#--------------------------------------------------#

#Filter IDs that have been identified as being predated
pred_id_dist <- muddyfoot_pike_roach_distances_df %>% 
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
F59709_pred_dat <- muddyfoot_filt_data %>% 
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
muddyfoot_polygon = sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))




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

#filter for plotting
F59709_pred_dat_filt <- F59709_pred_dat %>% 
  filter(date_utc2 > "2022-10-03")

#subset the data

# Define your ggplot object with the necessary layers
p <- ggplot() + 
  geom_sf(data = muddyfoot_polygon, fill = NA, color = "black", linewidth = 1) +
  geom_point(data = F59709_pred_dat_filt, 
             aes(x = longitude, 
                 y = latitude, 
                 group = individual_id,
                 color = individual_id), 
             size = 4) +  # Plot the fish tracks with color by individual_id
  ggsn::scalebar(
    data = muddyfoot_polygon, 
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
# projection(roach_muddyfoot_tel) <- projection(pike_muddyfoot_tel)
# projection(roach_muddyfoot_ctmm_fits) <- projection(pike_muddyfoot_ctmm_fits)
# 
# 
# combined_telemetry <- c(roach_muddyfoot_tel["F59701"], pike_muddyfoot_tel["F59885"])
# combined_ctmm <- c(roach_muddyfoot_ctmm_fits["F59701"], pike_muddyfoot_ctmm_fits["F59885"])
# 
# 
# #Pairwise separation distances
# PROX <- proximity(combined_telemetry,
#                   combined_ctmm)









