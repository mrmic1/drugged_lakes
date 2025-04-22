# ----------------------------------------------------------## 
# Exploring predator-prey interactions - Lake Cow Paradise ###
# ----------------------------------------------------------##

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
filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"
lake_polygon_path <- "./data/lake_coords/"
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/cow_paradise/"
enc_path <- "./data/encounters/cow_paradise/"
size_path <- "./data/fish_size/"
save_tables_path <- "./tables/cow_paradise/"  # Path to save summary tables.

### DATA LOADING ###
cow_filt_data <-  readRDS(paste0(filtered_data_path, "03_lake_cow_sub.rds"))

#--------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------#
# 1. Calculate distances at shared timepoints for prey and predators ##############
#---------------------------------------------------------------------------------#

# This section calculates pairwise distances between Roach/Perch (prey) and Pike (predator) at shared time points.
# Need the ctmms and telemetry objects for this part.

pike_cow_tel <- readRDS(paste0(telem_path, 'pike_lake_cow_tel.rds'))
perch_cow_tel <- readRDS(paste0(telem_path, 'perch_lake_cow_tel.rds'))
roach_cow_tel <- readRDS(paste0(telem_path, 'roach_lake_cow_tel.rds'))

pike_cow_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_cow_pike_fits/lake_cow_pike_OUF_models.rds"))
perch_cow_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_cow_perch_fits/lake_cow_perch_OUF_models.rds"))
roach_cow_ctmm_fits <- readRDS(paste0(ctmm_path, "lake_cow_roach_fits/lake_cow_roach_OUF_models.rds"))


# This section calculates pairwise distances between Roach/Perch (prey) and Pike (predator) at shared time points.

#### > 1.1. Roach ####

# Reproject Roach telemetry objects to match the Pike projection
ctmm::projection(roach_cow_tel) <- ctmm::projection(pike_cow_tel)
ctmm::projection(roach_cow_ctmm_fits) <- ctmm::projection(pike_cow_ctmm_fits)

#ctmm's with same projections
saveRDS(roach_cow_ctmm_fits, paste0(ctmm_path, "lake_cow_roach_fits/lake_cow_roach_OUF_models.rds"))
saveRDS(pike_cow_ctmm_fits, paste0(ctmm_path, "lake_cow_pike_fits/lake_cow_pike_OUF_models.rds"))

# Check projections to ensure they match
stopifnot(projection(roach_cow_tel) == projection(pike_cow_tel))

# Initialize an empty list to store distance calculations between Roach and Pike
roach_pike_distances <- list()

# Nested loop to calculate distances between each Roach and Pike telemetry pair
for(i in 1:length(roach_cow_tel)) {
  for(j in 1:length(pike_cow_tel)) {
    combined_telemetry <- c(roach_cow_tel[i], pike_cow_tel[j])
    combined_ctmm <- c(roach_cow_ctmm_fits[i], pike_cow_ctmm_fits[j])
    
    # Calculate pairwise distances
    location_difference <- distances(combined_telemetry, combined_ctmm)
    
    # Extract and store results with species IDs
    roach_id <- names(roach_cow_tel)[i]
    pike_id <- names(pike_cow_tel)[j]
    
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
saveRDS(roach_pike_distances_df, paste0(enc_path, "cow/cow_pike_roach_distances_df.rds"))


#### > 1.2. Perch ####

# Reproject Perch telemetry objects to match Pike
projection(perch_cow_tel) <- projection(pike_cow_tel)
projection(perch_cow_ctmm_fits) <- projection(pike_cow_ctmm_fits)

saveRDS(perch_cow_ctmm_fits, paste0(ctmm_path, "lake_cow_perch_fits/lake_cow_perch_OUF_models.rds"))

# Ensure projections are correctly matched
stopifnot(projection(perch_cow_tel) == projection(pike_cow_tel))

# Initialize an empty list to store distance calculations between Perch and Pike
perch_pike_distances <- list()

# Nested loop for calculating Perch - Pike distances
for(i in 1:length(perch_cow_tel)) {
  for(j in 1:length(pike_cow_tel)) {
    combined_telemetry <- c(perch_cow_tel[i], pike_cow_tel[j])
    combined_ctmm <- c(perch_cow_ctmm_fits[i], pike_cow_ctmm_fits[j])
    
    location_difference <- distances(combined_telemetry, combined_ctmm)
    
    perch_id <- names(perch_cow_tel)[i]
    pike_id <- names(pike_cow_tel)[j]
    
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
# saveRDS(perch_pike_distances_df, paste0(enc_path, "cow_pike_perch_distances_df.rds"))

#-------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------#
# 2. Calculate empirical encounters between prey and predators #########
#----------------------------------------------------------------------#

roach_pike_distances_df <- readRDS(paste0(enc_path, "cow_pike_roach_distances_df.rds"))
perch_pike_distances_df <- readRDS(paste0(enc_path, "cow_pike_perch_distances_df.rds"))

#-------------------------------------------------------------------------------------------#
# > 2.1. Total encounters #### 

# This section calculates encounters between prey and pike (predator), specifically focusing on close interactions
# (within 0.45 meters, which is the max strike distance of a Pike, according to Harper 1991).

# Define encounters as those occurring within 0.45 meters (based on Pike's strike distance)
roach_pike_distances_df$encounter <- ifelse(roach_pike_distances_df$est <= 0.45, 1, 0)
perch_pike_distances_df$encounter <- ifelse(perch_pike_distances_df$est <= 0.45, 1, 0)

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
#3992 rows

perch_pike_encounter_sum <- 
  perch_pike_distances_df %>%
  group_by(Perch_ID, Pike_ID, Date) %>%
  summarise(encounter_count = sum(encounter),
            daily_avg_dist_from_pike = mean(est, na.rm = TRUE),
            daily_min_dist_from_pike = min(est, na.rm = TRUE))
#7627 rows

# #save roach and perch predator encounter summaries
saveRDS(roach_pike_encounter_sum, paste0(enc_path, "cow_roach_pike_encounter_summary_unfiltered.rds"))
saveRDS(perch_pike_encounter_sum, paste0(enc_path, "cow_perch_pike_encounter_summary_unfiltered.rds"))


# Renaming the columns to have a common name before binding
roach_pike_encounter_sum <- roach_pike_encounter_sum %>%
  rename(Prey_ID = Roach_ID)

perch_pike_encounter_sum <- perch_pike_encounter_sum %>%
  rename(Prey_ID = Perch_ID)

# Combine the data frames into a single data frame
prey_pike_encounter_sum <- rbind(cbind(roach_pike_encounter_sum, Species = 'Roach'), 
                                 cbind(perch_pike_encounter_sum, Species = 'Perch'))
#11619 rows

#Need to filter out encounters with identified dead pike
#Dead pike track and mortality event identified by looking daily trajectories for pike that
#were not retrieved at the end of the experiment

pike_deaths <- read.csv("./data/encounters/pike_deaths.csv")

# Filter out encounters from dead pike #

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
pike_mort_cols <- 
  pike_deaths %>%
  filter(individual_ID %in% cow_filt_data$individual_ID) %>% 
  mutate(
    pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y")) %>%
  dplyr::select(-likely_death_date)

prey_pike_encounter_sum <-
  prey_pike_encounter_sum %>%
  left_join(pike_mort_cols, c("Pike_ID" = "individual_ID")) %>% 
  filter(is.na(pike_death_date)| Date <= pike_death_date) %>% 
  filter(Pike_ID != "F59889") #This pike died right at the beginning of the study

#original number of rows: 11619
#new: 8901

#save prey predator encounter summary
#saveRDS(prey_pike_encounter_sum, paste0(enc_path, "cow_prey_pike_encounter_summary_unfiltered.rds"))
#prey_pike_encounter_sum <- readRDS(paste0(enc_path, "cow_prey_pike_encounter_summary_unfiltered.rds"))

#----------------------------------------------------------------------------------------------------#

#> 2.2. Suspected predation events ####

# This section identifies potential predation events by filtering for days when pike and prey
# had more than 25 encounters or when their average daily distance was less than 0.2 meters.

post_biometrics <- 
  fread(paste0(size_path, "biometric_post_exp_data.csv")) %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number)))

#Summarise number of survivors and known death for each species in cow
post_biometrics %>% 
  filter(Lake == 'Cow Paradise') %>% 
  dplyr::select(individual_ID, Species, Treatment, Found) %>% 
  group_by(Species, Treatment) %>% 
  summarise(num_survive = sum(Found),
            num_died = length(unique(individual_ID)) - sum(Found),
            percent_survive = sum(Found)/length(unique(individual_ID)))


# Merge biometric data (e.g., whether the fish was found alive) with encounter summary
prey_post_size_cols <- 
  post_biometrics %>%
  filter(Lake == 'Cow Paradise', Species == 'Roach'| Species == 'Perch') %>%
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
#21 unique individuals

# Summarise prey-pike interactions: count days with > 25 encounters, consecutive days, and the first date
cow_prey_pike_interactions <- 
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
#27 rows

# Extract individual treatment and missing dates information from filtered data
prey_sel_cols <- 
  cow_filt_data %>%
  filter(Species == 'Roach'| Species == 'Perch') %>%
  dplyr::select(individual_ID, Treatment, n_missing_dates) %>%
  distinct()

# Merge encounter summary with individual-level metadata
cow_prey_pike_interactions <- 
  cow_prey_pike_interactions %>%
  left_join(prey_sel_cols, by = c("Prey_ID" = "individual_ID"))

# Merge bio metric data (e.g., whether the fish was found alive) with possible predator events
# Create column indicating instances where there might be a higher likelihood that an individual was predated
# specifically, when individuals 2 or more consecutive days where they had more than 25 encounters with a predator
cow_prey_pike_interactions <- 
  cow_prey_pike_interactions %>%
  left_join(prey_post_size_cols, by = c("Prey_ID" = "individual_ID")) %>%
  # Remove fish found alive at the end of the experiment
  filter(found_alive == '0') %>%
  # Flag likely predation events (i.e., when consecutive days of over 100 encounters are more than 2 days)
  mutate(likely_predated = ifelse(consecutive_days_100 >= 2, 1, 0))

# Identify the date with the most encounters for each roach and pike pair
cow_prey_max_encounter_dates <- 
  prey_pike_encounter_sum %>%
  group_by(Prey_ID, Pike_ID) %>%
  filter(encounter_count == max(encounter_count)) %>%
  slice(1) %>%
  dplyr::select(Prey_ID, Pike_ID, Species, max_encounter_date = Date, max_encounter_count = encounter_count)

#Add max_encounter_date and max_encounter_count to potential predation events summary.
cow_prey_pike_interactions <- 
  cow_prey_pike_interactions %>%
  left_join(cow_prey_max_encounter_dates %>%
              dplyr::select(Prey_ID, Pike_ID, max_encounter_date, max_encounter_count),
            by = c("Prey_ID" = "Prey_ID", 
                   "Pike_ID" = "Pike_ID"))

# Save the summarised potential predation events for future use
saveRDS(cow_prey_max_encounter_dates, paste0(enc_path, "cow_prey_daily_max_pred_encounters.rds"))

# Save the suspected predation event data
saveRDS(cow_prey_pike_interactions, paste0(enc_path, "cow_suspected_predation_events.rds"))

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

cow_prey_daily_encounter_summary <-
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

cow_prey_total_encounter_summary <-
  cow_prey_daily_encounter_summary %>%
  group_by(Prey_ID) %>%
  summarise(Species = first(Species),
            total_encounters = sum(total_daily_encounter_count),
            min_daily_dist_from_pike = min(daily_min_dist_from_pike),
            avg_daily_min_dist_from_pike = min(daily_avg_dist_from_pike, na.rm = TRUE))


# Merge total encounter summary with metadata for each individual.
cow_prey_total_encounter_summary <-
  cow_prey_total_encounter_summary %>%
  left_join(prey_sel_cols, by = c("Prey_ID" = "individual_ID"))

#also add the number of poor tracking days info
#the loaded dataframe was created in script 03_cow_tracking_data_summary
cow_ID_irreg_sampling <- readRDS(paste0(filtered_data_path, "cow_ID_irreg_sampling.rds"))
prey_irreg_sampling_cols <- 
  cow_ID_irreg_sampling %>% 
  filter(Species == 'Roach'| Species == 'Perch') %>% 
  dplyr::select(individual_ID, n_poor_tracking_days, found_alive)

cow_prey_total_encounter_summary <-
  cow_prey_total_encounter_summary %>%
  left_join(prey_irreg_sampling_cols, by = c("Prey_ID" = "individual_ID"))

cow_prey_total_encounter_summary <-
  cow_prey_total_encounter_summary %>% 
  mutate(days_tracked = 35 - n_missing_dates,  # Total tracking duration (max 35 days)
         avg_daily_pred_encounter_rate = round(total_encounters / days_tracked, 0))

head(cow_prey_total_encounter_summary)


# Save the unfiltered encounter summary
saveRDS(cow_prey_total_encounter_summary, paste0(enc_path, 'cow_prey_total_encounter_summary_unfiltered.rds'))
saveRDS(cow_prey_daily_encounter_summary, paste0(enc_path, 'cow_prey_daily_encounter_summary_unfiltered.rds'))

#------------------------------------------------------------------------------------------------------#


# > 3.3. Daily and total encounter summary tables for prey not found ####

#Reduce summary data to only include individuals that were not found at the end of the experiment
cow_daily_encounter_summary_prey_not_found <- 
  cow_prey_daily_encounter_summary %>% 
  filter(found_alive == 0)
#326 rows


#Add column to identify whether the date was irregularly sampled or dates were not tracked
#First I need to extract the required columns
prey_not_found_cols <- 
  cow_ID_irreg_sampling %>% 
  filter(Species == 'Roach'| Species == 'Perch') %>% 
  dplyr::select(individual_ID, dates_irregular_positions, missing_dates)

#Then we need to do some data wrangling to extract the date information
cow_daily_encounter_summary_prey_not_found$Date <- as.Date(cow_daily_encounter_summary_prey_not_found$Date)

# Ensure dates_irregular_positions and missing_dates are properly formatted as lists of dates
prey_not_found_cols$dates_irregular_positions <- lapply(strsplit(prey_not_found_cols$dates_irregular_positions, ", "), as.Date)
prey_not_found_cols$missing_dates <- lapply(strsplit(prey_not_found_cols$missing_dates, ", "), as.Date)

# Merge the two datasets on the IDs
merged_df <- merge(cow_daily_encounter_summary_prey_not_found, prey_not_found_cols, 
                   by.x = "Prey_ID", by.y = "individual_ID", all.x = TRUE)

# Create the two new columns
merged_df$poor_tracking_date <- mapply(function(date, irregular_dates) {
  if (!is.null(irregular_dates) && date %in% irregular_dates) 1 else 0
}, merged_df$Date, merged_df$dates_irregular_positions)

merged_df$no_tracking_date <- mapply(function(date, missing_dates) {
  if (!is.null(missing_dates) && date %in% missing_dates) 1 else 0
}, merged_df$Date, merged_df$missing_dates)

# Keep only the original columns from daily_encounter_sum_roach_found and the new ones
cow_daily_encounter_summary_prey_not_found <- merged_df[, c("Prey_ID", "Date", "Species", "daily_avg_dist_from_pike","total_daily_encounter_count", 
                                                           "max_daily_encounter_count", "pike_id_max_encounter_count",
                                                           "daily_min_dist_from_pike", "pike_id_min_dist", "found_predated",
                                                           "poor_tracking_date", "no_tracking_date")]

# Save the unfiltered encounter summary
saveRDS(cow_daily_encounter_summary_prey_not_found, paste0(enc_path, 'cow_daily_encounter_summaries_for_prey_not_found_unfiltered.rds'))

length(unique(cow_daily_encounter_summary_prey_not_found$Prey_ID))
#12

#Total encounter summary for prey not found
#This summarised dataframe will include information on
#total encounter count
#max encounter count
#date with max encounter count
#number of days with over 50 encounters
#first date with over 50 encounters

cow_total_encounter_summary_prey_not_found <-
  cow_daily_encounter_summary_prey_not_found %>% 
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
  cow_prey_pike_interactions %>%
  dplyr::select(Prey_ID, num_days_25_encounters, num_days_100_encounters, num_days_min_dist_less_0.2m,  
                first_date_over_25,consecutive_days_25, first_date_over_100, consecutive_days_100) %>% 
  group_by(Prey_ID) %>% 
  arrange(desc(num_days_100_encounters),
          desc(num_days_100_encounters),
          .by_group = TRUE) %>% 
  slice(1) %>% 
  ungroup()

cow_total_encounter_summary_prey_not_found <-
  cow_total_encounter_summary_prey_not_found %>%
  left_join(prey_not_found_cols, by = c("Prey_ID" = "Prey_ID")) %>% 
  dplyr::select(-Pike_ID)


#merge with irregular sampling summary
#create in script 03_cow_tracking_data_summary
cow_ID_irreg_sampling_not_found <- 
  cow_ID_irreg_sampling %>% 
  filter(found_alive == 0)

#get required rows from filtered encounters
filtered_encounters_rows <- 
  cow_ID_irreg_sampling_not_found %>% 
  dplyr::select(individual_ID, n_poor_tracking_days, total_missing_days)

cow_total_encounter_summary_prey_not_found <-
  cow_total_encounter_summary_prey_not_found %>%
  left_join(filtered_encounters_rows, by = c("Prey_ID" = "individual_ID")) %>% 
  dplyr::select(Prey_ID, Species, everything())

# Create a flextable summarizing the irregular and missing days for missing individuals
cow_daily_irregular_sample_table <- 
  flextable(cow_total_encounter_summary_prey_not_found) %>% 
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

save_as_docx(cow_daily_irregular_sample_table, 
             path = paste0(save_tables_path, "cow IDs not found encounter and tracking summary.docx"))

saveRDS(cow_total_encounter_summary_prey_not_found, paste0(enc_path, 'cow_total_encounter_summary_for_prey_not_found_unfiltered.rds'))

# Load the openxlsx package
if (!require(openxlsx)) install.packages("openxlsx")
library(openxlsx)

# Save the dataframe as an Excel file
write.xlsx(cow_total_encounter_summary_prey_not_found, file = paste0(save_tables_path, "cow_IDs_not_found_total_encounter_summary.xlsx"))

#create seperate table for dates

cow_dates_table <- 
  cow_ID_irreg_sampling_not_found %>% 
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

save_as_docx(cow_dates_table, 
             path = paste0(save_tables_path, "cow_IDs_not_found_tracking_dates.docx"))

#--------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------#

#------------------------------#
# 4. Filter encounter data #### 
#------------------------------#

prey_encounters_sum_IDs_not_found <- readRDS(paste0(enc_path, "cow_total_encounter_summary_for_prey_not_found_unfiltered.rds"))
mortality_preds <- readxl::read_excel("./data/encounters/suspected_mortality_updated.xlsx")
prey_encounter_sum <- readRDS(paste0(enc_path, "cow_prey_pike_encounter_summary_unfiltered.rds"))

#-------------------------------------------------------------------------------------------------------------------#

#> 4.1. Filter out encounters for dead prey ####

#Extracted individuals that we believe to be predated
cow_pred_prey_cols <-
  mortality_preds %>%
  filter(lake == 'cow paradise') %>% 
  filter(species == 'Roach'| species == 'Perch') %>%
  dplyr::select(individual_ID, revised_suspected_mortality, revised_likely_death_date) %>% 
  rename(death_date = revised_likely_death_date)

prey_encounter_sum_filtered <-
  prey_encounter_sum %>%
  left_join(cow_pred_prey_cols, c("Prey_ID" = "individual_ID")) %>% 
  # Keep only pre-predation data
  filter(is.na(death_date)| Date <= death_date)

#original number of rows: 8901
#new: 8531

#Filter out individuals with poor tracking

prey_encounter_sum_filtered <- 
  prey_encounter_sum_filtered %>% 
  filter(revised_suspected_mortality != "poor_tracking_remove" | is.na(revised_suspected_mortality))

#new: 8396

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
  cow_filt_data %>%
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
  filter(Lake == 'Cow Paradise', Species == 'Roach'| Species == 'Perch') %>%
  dplyr::select(individual_ID, Found, Known_predated) %>% 
  rename(found_alive = Found,
         found_predated = Known_predated)

# Add information about whether the Perch were found alive at the end of the experiment
prey_total_encounter_summary_filtered <- 
  prey_total_encounter_summary_filtered %>%
  left_join(post_size_cols, by = c("Prey_ID" = "individual_ID"))

# Save the combined predation event summary
saveRDS(prey_total_encounter_summary_filtered, paste0(enc_path, "cow_pred_prey_encounter_summary_filtered.rds"))
