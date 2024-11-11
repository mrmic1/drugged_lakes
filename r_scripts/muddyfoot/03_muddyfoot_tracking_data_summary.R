#-------------------------------------------#
### EXTRACTING TRACKING SUMMARY STATISTICS
#-------------------------------------------#

#This script calculates tracking summary statistics prior to any filtering based on predation events 
#It also adds data columns for missing data and irregular sampling which will be used to help to identify predation events
#This script create summary tables for tracking data that could go into a supplementary section. 

# LIBRARIES
# Loading essential libraries for data manipulation, spatial analysis, and formatting tables.
library(tidyverse)  # For data manipulation and visualization.
library(sp)  # For handling spatial data.
library(sf)  # For handling spatial vector data using simple features.
library(flextable)  # For creating and customizing tables.
library(kableExtra)  # Additional table formatting options.
library(officer)  # For exporting tables to Word documents.
library(ctmm)  # For continuous-time movement modeling (ctmm) package.
library(lubridate)  # For handling date-time data.

# Set the time zone to ensure consistent time handling
Sys.setenv(TZ = 'Europe/Stockholm')

# SET TABLE PLOTTING PARAMETERS
# Customize default settings for Flextable tables used later in the script.
set_flextable_defaults(
  font.color = "black",
  border.color = "black",
  font.family = 'Arial',
  line_spacing = 1
)

### DIRECTORIES ###
# Define paths to various datasets and save locations for tables.
ctmm_path <- "./data/ctmm_fits/"  # Path for ctmm model fits.
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"  # Path for filtered tracking data.
telem_path <- "./data/telem_obj/"  # Path for telemetry object files.
save_tables_path <- "./tables/muddyfoot/"  # Path to save summary tables.
enc_path <- "./data/encounters/"

### LOAD DATA ###

muddyfoot_filt_data <-  readRDS(paste0(filtered_data_path, "02_muddyfoot_sub.rds"))

# Load telemetry objects for different species.
pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

# # Load ctmm model fits for each species.
# pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
# perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
# roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds"))

#--------------------------------------------------------#
# 1. Calculate tracking summary metrics #################
#--------------------------------------------------------#

#Below we calculate the following metrics
#- mean and median time difference between timestamps
#- number of locations per individual
#- number of days individuals were tracked for.
#- number of locations per day, hour and min
#- effective sample size for each individual used for calculating akdes
#- daily location information

# > 1.1 Calculate time difference between positions ####

# Calculate time differences between successive timestamps for each individual.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_ID) %>% 
  mutate(time_diff = c(NA, diff(timestamp)))  # First difference is NA.

# Round the time differences to 3 decimal places for clarity.
muddyfoot_filt_data$time_diff <- as.numeric(round(muddyfoot_filt_data$time_diff, digits = 3))

# Check the first 20 rows to ensure time differences are calculated correctly.
head(muddyfoot_filt_data %>% dplyr::select(timestamp, time_diff), n = 20)

# Calculate mean and median time differences per individual.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_ID) %>% 
  mutate(mean_time_diff = mean(time_diff, na.rm = TRUE),
         median_time_diff = median(time_diff, na.rm = TRUE)) %>% 
  ungroup()

# > 1.2 Calculate number of positions per individual ####

# Count the number of positions (i.e., tracking locations) for each individual.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_ID) %>% 
  mutate(n_positions = n()) %>% 
  ungroup()

# Print a summary of unique individuals and their corresponding number of positions.
print(muddyfoot_filt_data %>% 
        dplyr::select(individual_ID, n_positions) %>% 
        distinct(), 
      n = 65)

# > 1.3 Calculate number of days individuals were tracked ####

# Calculate the number of unique days each individual was monitored.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_ID) %>% 
  mutate(n_days_tracked = length(unique(Date))) %>% 
  ungroup()

# Optional: Check the summary of number of days and positions per individual.
print(muddyfoot_filt_data %>% 
        dplyr::select(individual_ID, Treatment, n_positions, n_days_tracked) %>% 
        distinct(), 
      n = 65)


# Create a summary table of the calculated metrics.
positions_sum <- 
  muddyfoot_filt_data %>% 
  dplyr::select(individual_ID, Treatment, Species,
                n_positions, 
                n_days_tracked, 
                mean_time_diff,
                median_time_diff) %>% 
  distinct()


#> 1.4. Calculate daily positions sum stats ####

# This section calculates daily tracking metrics for each individual, 
# including the number of positions per day and time differences between positions.

# Calculate the number of positions per day and median time differences between positions for each individual.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_ID, Date) %>% 
  mutate(n_positions_day = n(),
         median_day_time_diff = median(time_diff, na.rm = TRUE)) %>%  # Median time difference per day.
  ungroup()

# Calculate the number of positions per hour and per minute.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_ID, Date) %>% 
  mutate(n_positions_hourly = n_positions_day / 24,  # Average number of positions per hour.
         n_positions_per_min = n_positions_hourly / 60) %>%  # Average number of positions per minute.
  ungroup()

# Optional: Check a specific individual for daily tracking stats.
print(
  muddyfoot_filt_data %>% 
    filter(individual_ID == 'F59704') %>%  # Example individual ID.
    dplyr::select(individual_ID, Species, Date, n_positions_day, median_day_time_diff, n_positions_hourly, n_positions_per_min) %>% 
    distinct(), n = 36
)


#------------------------------------------------------------------------------------------------#

#-------------------------------------------#
# 2. Find periods of irregular sampling ####
#-------------------------------------------#

# This section identifies days with irregular sampling based on the number of positions and time differences.
# Summarise average, standard deviation, minimum, and maximum time differences between positions for each individual per day.
daily_sampling_sum <- 
  muddyfoot_filt_data %>% 
  group_by(individual_ID, Treatment, Date) %>%
  summarise(
    avg_time_diff = mean(time_diff, na.rm = TRUE),
    stdev_time_diff = sd(time_diff, na.rm = TRUE),
    min_time_diff = min(time_diff, na.rm = TRUE),
    max_time_diff = max(time_diff, na.rm = TRUE),
    n_days_tracked = first(n_days_tracked)  # Number of days monitored.
  )

# > 2.1 Dates that individuals were not tracked ####

# Identify dates where an individual was not tracked at all by comparing against a full date range.
full_date_range <- seq(as.Date("2022-09-25"), as.Date("2022-10-30"), by = "day")  # Define full date range.

# Create a table of all possible combinations of individual_id and dates.
all_combinations <- 
  expand.grid(
  individual_ID = unique(muddyfoot_filt_data$individual_ID),
  Date = full_date_range
)

# Merge with original tracking data to include species and treatment information.
all_combinations <- 
  all_combinations %>%
  left_join(muddyfoot_filt_data %>% 
              dplyr::select(individual_ID, Species, Treatment) %>% 
              distinct(), 
            by = "individual_ID")


# Perform an anti-join to identify missing tracking dates.
missing_dates <- 
  all_combinations %>%
  anti_join(muddyfoot_filt_data, by = c("individual_ID", "Date")) %>%
  arrange(individual_ID, Date)


# Create a summary of missing dates for each individual.
summary_missing_dates <- 
  missing_dates %>%
  group_by(individual_ID) %>%
  summarise(Species = first(Species),
            Treatment = first(Treatment),
            missing_dates = paste(Date, collapse = ", "),  # List missing dates.
            n = n(),  # Count of missing dates.
            n_breaks = if_else(n() == 1, 0, sum(diff(Date) != 1)))  # Number of breaks in consecutive dates.

# Create a table of missing dates using flextable.
missing_dates_table <- 
  flextable(summary_missing_dates) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_ID" = 'ID', 
                    "missing_dates" = 'Dates ID was not tracked') %>% 
  width(j = 4, 11, unit = 'cm')

# Optional: Save the missing dates table as a Word document.
save_as_docx(missing_dates_table, path = paste0(save_tables_path, "muddyfoot_IDs_missing_dates.docx"))


#> 2.2 Add information about missing dates to positions summary ####

# Combine missing dates information with the positions summary.
positions_sum <- 
  positions_sum %>% 
  left_join(summary_missing_dates[, c('individual_ID', 'n', 'n_breaks')], by = "individual_ID") %>%
  mutate(n_missing_dates = ifelse(is.na(n), 0, n),  # If no missing dates, set to 0.
         n_breaks = ifelse(is.na(n_breaks), 0, n_breaks)) %>% 
  mutate_if(is.numeric, round, 1)  # Round all numeric columns.

# Create a summary table of positions using flextable.
positions_sum_table <- 
  flextable(positions_sum %>% 
              dplyr::select(-n, -mean_time_diff)) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_ID" = 'ID',
                    "treatment" = 'Treatment',
                    "n_positions" = 'Number of locations',
                    "n_days_tracked" = 'Number of days tracked',
                    "median_time_diff" = 'Median location frequency',
                    "n_breaks" = "Date sequence breaks in tracking",
                    "n_missing_dates" = "Number of untracked days") 

# Optional: Save the positions summary table as a Word document.
save_as_docx(positions_sum_table, path = paste0(save_tables_path, "muddyfoot_ID_loc_sum_pred_events_unfilt.docx"))

# Update the main dataset by adding effective sample size and missing dates.
muddyfoot_filt_data <- merge(muddyfoot_filt_data, 
                             positions_sum[, c("individual_ID", "n_breaks", "n_missing_dates")],
                             by = "individual_ID", 
                             all.x = TRUE)


# > 2.3 Extract individual IDs that had irregularities in tracking ####

# This section identifies individuals with irregular tracking patterns, possibly due to predation, 
# equipment failure, or battery issues. These individuals may need to have adjusted home range estimates.

# Summarise the range of locations (n_positions) to determine thresholds for irregular sampling.
positions_sum %>% 
  dplyr::select(n_positions) %>% 
  summarise(mean = mean(n_positions, na.rm = TRUE),
            median = median(n_positions, na.rm = TRUE),
            stdev = sd(n_positions, na.rm = TRUE),
            min = min(n_positions, na.rm = TRUE),
            max = max(n_positions, na.rm = TRUE))

# Set a threshold for irregular sampling as -1 standard deviation from the mean number of positions.
# mean: 159364 - sd: 95511 = 63853
irreg_indiv <- positions_sum %>% 
  filter(n_positions < 63853 | n_missing_dates >= 10 | n_breaks >= 2)

# Add a flag for individuals with irregular sampling in the main dataset.
irregular_ids <- irreg_indiv$individual_ID

muddyfoot_filt_data <- 
  muddyfoot_filt_data %>% 
  mutate(irregular_sampling = ifelse(individual_ID %in% irregular_ids, 1, 0))

### Identify irregular days ###

daily_pos_irreg_ind <- 
  muddyfoot_filt_data %>% 
  group_by(Species) %>% 
  mutate(avg_median_time_diff = mean(median_time_diff, na.rm = TRUE),
         std_median_time_diff = sd(median_time_diff, na.rm = TRUE),
         avg_daily_positions = mean(n_positions_day, na.rm = TRUE),
         std_daily_positions = sd(n_positions_day, na.rm = TRUE),
         avg_daily_hour_positions = mean(n_positions_hourly, na.rm = TRUE),
         std_daily_hour_positions = sd(n_positions_hourly, na.rm = TRUE)
  ) %>% 
  # Identify days where tracking is below thresholds for irregular sampling.
  filter((median_day_time_diff <= (avg_median_time_diff - 3 * std_median_time_diff)) |
           (n_positions_day <= (avg_daily_positions - 1 * std_daily_positions)) |
           (n_positions_hourly <= (avg_daily_hour_positions - 1 * std_daily_hour_positions))) %>%
  dplyr::select(Species, individual_ID, Date, median_day_time_diff, n_positions_day, n_positions_hourly) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  filter((n_positions_day <= 1000)) %>%  
  filter((median_day_time_diff > 5)) %>%  
  filter((n_positions_hourly < 30)) %>% 
  distinct(individual_ID, Date, .keep_all = TRUE) %>% 
  ungroup()

# Flag irregular sampling days in the main dataset.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  mutate(irreg_sample_day = ifelse(paste(individual_ID, Date) %in% paste(daily_pos_irreg_ind$individual_ID, 
                                                                         daily_pos_irreg_ind$Date), 1, 0))

# Check the flagged data for irregular days.
check <- print(
  muddyfoot_filt_data %>% 
    filter(irreg_sample_day == 1) %>% 
    dplyr::select(individual_ID, Date, n_positions_day, median_day_time_diff, n_positions_hourly) %>% 
    mutate_if(is.numeric, round, 1) %>% 
    distinct()
)

# Create a table of identified irregular sampling days.
individual_irregular_sample_table <- 
  flextable(check) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_ID" = 'ID', 
                    "Date" = "Date",
                    "n_positions_day" = "Number of locations",
                    "median_day_time_diff" = "Median location frequency (s)",
                    "n_positions_hourly" = "Number of location per hour") %>% 
  width(j = c(2:5), width = 3, unit = 'cm')

# Optional: Save the irregular sampling table as a Word document.
save_as_docx(individual_irregular_sample_table, path = paste0(save_tables_path, "muddyfoot_ID_day_locs_w_irregular_sampling.docx"))


# > 2.4. Summarise irregular days and missing data per individual ####

# Reduce the data to one row per individual per date, summarising irregular and missing days.
reduced_data <- 
  muddyfoot_filt_data %>%
  group_by(Species, individual_ID, Date) %>%
  summarise(
    irreg_sample_day = max(irreg_sample_day),  # Maximum indicates any irregular sampling.
    n_missing_dates = max(n_missing_dates),  # Maximum for missing dates.
    .groups = 'drop'
  )

# Summarise the number of irregular and missing days per individual.
summary_data <- 
  reduced_data %>%
  group_by(Species, individual_ID) %>%
  summarise(
    n_irregular_days = sum(irreg_sample_day == 1),  # Count irregular days.
    total_missing_days = max(n_missing_dates),  # Count missing days.
    total_irregular_or_missing_days = n_irregular_days + total_missing_days,  # Combine irregular and missing days.
    dates_irregular_positions = paste(Date[irreg_sample_day == 1], collapse = ", "),  # List irregular dates.
    .groups = 'drop'
  )

# Merge the missing dates information.
missing_dates_sub <- 
  summary_missing_dates %>% 
  dplyr::select(individual_ID, missing_dates)

summary_data_test <- 
  summary_data %>% 
  left_join(missing_dates_sub, by = 'individual_ID')

# Create a flextable summarizing the irregular and missing days.
daily_irregular_sample_table <- 
  flextable(summary_data_test) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_ID" = 'ID', 
                    "n_irregular_days" = "Number of days with irregular tracking",
                    "total_missing_days" = "Number of days not tracked",
                    "total_irregular_or_missing_days" = "Total number of poor tracking days",
                    "dates_irregular_positions" = "Dates with irregular sampling",
                    "missing_dates" = 'Dates ID was not tracked') %>% 
  width(j = c(6,7), 13, unit = 'cm')

# Optional: Save the daily irregular sample table as a Word document.
save_as_docx(daily_irregular_sample_table, 
             path = paste0(save_tables_path, "muddyfoot_ID_irregular_sampling_summary.docx"))


#--------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------#

#SAVE
saveRDS(summary_data_test, paste0(filtered_data_path, "muddyfoot_ID_irreg_sampling.rds"))
saveRDS(muddyfoot_filt_data, paste0(filtered_data_path, "03_muddyfoot_sub.rds"))

