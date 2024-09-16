#-------------------------------------------#
### EXTRACTING TRACKING SUMMARY STATISTICS
#-------------------------------------------#

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
data_filter_path <- "./data/tracks_filtered/muddyfoot/"  # Path for filtered tracking data.
telem_path <- "./data/telem_obj/"  # Path for telemetry object files.
save_tables_path <- "./data/tracks_filtered/sum_tables/"  # Path to save summary tables.
enc_path <- "./data/encounters/"

### LOAD DATA ###
# Load telemetry objects for different species.
pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

# Load ctmm model fits for each species.
pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds"))

#--------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------#
#> Make combined dataframe to extract summary statistics ####
#---------------------------------------------------------------#

# Create a function to extract the data and add the individual_id column
# This function merges telemetry data from different individuals into a single dataframe.
extract_data <- function(list_data) {
  data_combined <- do.call(rbind, lapply(names(list_data), function(id) {
    df <- list_data[[id]]
    df$individual_id <- id  # Add individual ID to each dataframe.
    return(df)
  }))
  return(data_combined)
}

# Apply the function to each species' telemetry data and combine them into a single dataframe.
pike_data <- extract_data(pike_muddyfoot_tel)
perch_data <- extract_data(perch_muddyfoot_tel)
roach_data <- extract_data(roach_muddyfoot_tel)

# Combine all dataframes into one, with a species label added.
muddyfoot_telem_data <- 
  rbind(
    cbind(pike_data, Species = 'Pike'), 
    cbind(perch_data, Species = 'Perch'),
    cbind(roach_data, Species = 'Roach')
  )

#saveRDS(muddyfoot_telem_data, paste0(data_filter_path, "02_muddyfoot_sub.rds"))

#------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------#
# 2. Filter out post-predation event data ####
#-------------------------------------------------#

# This section filters tracking data to remove data from prey individuals post-predation events. 
# Predation events were identified in a separate script and are used to clean the data.

# Load the predation event dataframe (identified in the 'muddyfoot_species_interactions.R' script).
mud_pred_events <- readRDS(paste0(enc_path, "muddyfoot_pred_events.rds"))

# Select relevant columns from predation event data to identify the first date the prey was tracked post-predation.
pred_cols <- mud_pred_events %>%
  dplyr::select(individual_ID, first_date_over_50)

# Filter out data after the predation event for each individual.
muddyfoot_filt_data <- 
  muddyfoot_telem_data %>%
  left_join(pred_cols, by = c("individual_id" = "individual_ID")) %>%
  filter(is.na(first_date_over_50) | date <= first_date_over_50)  # Keep only pre-predation data

# Check how many rows were removed during filtering (optional).
# For example, 293,214 rows were removed, as indicated in the comment.

print(paste0("Rows removed after predation filtering: ", nrow(muddyfoot_telem_data) - nrow(muddyfoot_filt_data)))

# #check filtering worked
# # Filter data for individual 'F59709'
# filtered_test <- test %>% filter(individual_id == 'F59709')
# filtered_muddyfoot <- muddyfoot_filt_data %>% filter(individual_id == 'F59709')
# # Calculate the difference in the number of rows
# abs(nrow(filtered_test) - nrow(filtered_muddyfoot))

#saveRDS(muddyfoot_filt_data, paste0(data_filter_path, "03_muddyfoot_sub.rds"))

#-------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------#
# 3. Calculate tracking summary metrics #################
#--------------------------------------------------------#

#Below we calculate the following metrics
#- mean and median time difference between timestamps
#- number of locations per individual
#- number of days individuals were tracked for.
#- number of locations per day, hour and min
#- effective sample size for each individual used for calculating akdes
#- daily location information

# >>> 3.1 Calculate time difference between positions ####

# Calculate time differences between successive timestamps for each individual.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(time_diff = c(NA, diff(timestamp)))  # First difference is NA.

# Round the time differences to 3 decimal places for clarity.
muddyfoot_filt_data$time_diff <- as.numeric(round(muddyfoot_filt_data$time_diff, digits = 3))

# Optional: Check the first 20 rows to ensure time differences are calculated correctly.
# head(muddyfoot_filt_data %>% dplyr::select(timestamp, time_diff), n = 20)

# Calculate mean and median time differences per individual.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(mean_time_diff = mean(time_diff, na.rm = TRUE),
         median_time_diff = median(time_diff, na.rm = TRUE)) %>% 
  ungroup()

# >>> 3.2 Calculate number of positions per individual ####

# Count the number of positions (i.e., tracking locations) for each individual.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(n_positions = n()) %>% 
  ungroup()

# Print a summary of unique individuals and their corresponding number of positions.
print(muddyfoot_filt_data %>% 
        dplyr::select(individual_id, n_positions) %>% 
        distinct(), 
      n = 65)

# >>> 3.3 Calculate number of days individuals were tracked ####

#create date column in correct timezone

muddyfoot_filt_data$Date <- format(with_tz(ymd_hms(muddyfoot_filt_data$timestamp), 
                                           tzone = "Europe/Stockholm"), "%Y/%m/%d")
muddyfoot_filt_data$Date <- as.POSIXct(muddyfoot_filt_data$Date, tz = "Europe/Stockholm")

#check time
tz(muddyfoot_filt_data$Date) #Europe/Stockholm
tz(muddyfoot_filt_data$date) #UTC
tz(muddyfoot_filt_data$timestamp) #Europe/Stockholm

#to avoid issues I am going to remove the date column
muddyfoot_filt_data <- muddyfoot_filt_data %>% 
  dplyr::select(-date)

# Calculate the number of unique days each individual was monitored.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(n_days_mon = length(unique(Date))) %>% 
  ungroup()

# Optional: Check the summary of number of days and positions per individual.
print(muddyfoot_filt_data %>% 
        dplyr::select(individual_id, treatment, n_positions, n_days_mon) %>% 
        distinct(), 
      n = 65)


# Create a summary table of the calculated metrics.
positions_sum <- 
  muddyfoot_filt_data %>% 
  dplyr::select(individual_id, treatment, Species,
                n_positions, 
                n_days_mon, 
                mean_time_diff) %>% 
  distinct()


# >>> 3.4 Extract effective sample sizes for akde estimation ####

# The effective sample size is important for akde (Autocorrelated Kernel Density Estimation) 
# and helps to account for autocorrelation in the tracking data.

# Combine all ctmm model fits into a list for easier access.
muddyfoot_all_ctmm_fits <- list(Pike_fits = pike_muddyfoot_ctmm_fits, 
                                Perch_fits = perch_muddyfoot_ctmm_fits, 
                                Roach_fits = roach_muddyfoot_ctmm_fits)

# Initialize an empty dataframe to store summary results.
summary_df <- data.frame(
  species = character(),
  ID = character(),
  effective_n = numeric(),
  stringsAsFactors = FALSE
)

# List of species to loop through.
species_list <- c("Pike_fits", "Perch_fits", "Roach_fits")

# Loop through each species and extract the effective sample size (DOF - "area").
for (species in species_list) {
  individuals <- names(muddyfoot_all_ctmm_fits[[species]])  # Get individual IDs.
  
  # For each individual, extract the effective sample size from the ctmm model summary.
  for (ID in individuals) {
    summary_obj <- summary(muddyfoot_all_ctmm_fits[[species]][[ID]])  # Extract model summary.
    effective_n <- summary_obj$DOF["area"]  # Extract the effective sample size (DOF).
    
    # Append the result to the summary dataframe.
    summary_df <- rbind(summary_df, data.frame(
      species = gsub("_fits", "", species),  # Remove '_fits' from species name.
      ID = ID,
      effective_n = effective_n
    ))
  }
}

# Print the summary dataframe of effective sample sizes.
print(summary_df)

# Replace effective sample sizes with NA for individuals that were predated.
summary_df <- summary_df %>% 
  mutate(effective_n = ifelse(ID %in% mud_pred_events$individual_ID, NA, effective_n))

# Store the final effective sample size summary.
effective_n_sum <- summary_df

# Combine the effective sample sizes with individual tracking positions summary.
positions_sum <- merge(positions_sum, 
                       effective_n_sum[, c("ID", "effective_n")], 
                       by.x = "individual_id", 
                       by.y = "ID", 
                       all.x = TRUE)

# Check the correlation between the number of positions and the effective sample size.
cor(positions_sum$n_positions, positions_sum$effective_n, method = 'spearman', use = 'complete.obs')
# Example output: 0.73 - indicating a strong positive correlation.


#>>> 3.5. Calculate daily positions sum stats ####

# This section calculates daily tracking metrics for each individual, 
# including the number of positions per day and time differences between positions.

# Calculate the number of positions per day and median time differences between positions for each individual.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id, Date) %>% 
  mutate(n_positions_day = n(),
         median_day_time_diff = median(time_diff, na.rm = TRUE)) %>%  # Median time difference per day.
  ungroup()

# Calculate the number of positions per hour and per minute.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id, Date) %>% 
  mutate(n_positions_hourly = n_positions_day / 24,  # Average number of positions per hour.
         n_positions_per_min = n_positions_hourly / 60) %>%  # Average number of positions per minute.
  ungroup()

# Optional: Check a specific individual for daily tracking stats.
print(
  muddyfoot_filt_data %>% 
    filter(individual_id == 'F59704') %>%  # Example individual ID.
    dplyr::select(individual_id, Species, Date, n_positions_day, median_day_time_diff, n_positions_hourly, n_positions_per_min) %>% 
    distinct(), n = 36
)


#------------------------------------------------------------------------------------------------#

#-------------------------------------------#
# 4. Find periods of irregular sampling ####
#-------------------------------------------#

# This section identifies days with irregular sampling based on the number of positions and time differences.
# It excludes data after predation events, as those individuals are already filtered out.

# Summarize average, standard deviation, minimum, and maximum time differences between positions for each individual per day.
daily_sampling_sum <- 
  muddyfoot_filt_data %>% 
  group_by(individual_id, treatment, Date) %>%
  summarise(
    avg_time_diff = mean(time_diff, na.rm = TRUE),
    stdev_time_diff = sd(time_diff, na.rm = TRUE),
    min_time_diff = min(time_diff, na.rm = TRUE),
    max_time_diff = max(time_diff, na.rm = TRUE),
    n_days_mon = first(n_days_mon)  # Number of days monitored.
  )

# >>> 4.1 Dates that individuals were not tracked ####

# Identify dates where an individual was not tracked at all by comparing against a full date range.
full_date_range <- seq(as.Date("2022-09-25"), as.Date("2022-10-30"), by = "day")  # Define full date range.

# Create a table of all possible combinations of individual_id and dates.
all_combinations <- expand.grid(
  individual_id = unique(muddyfoot_filt_data$individual_id),
  Date = full_date_range
)

# Merge with original tracking data to include species and treatment information.
all_combinations <- all_combinations %>%
  left_join(muddyfoot_filt_data %>% 
              dplyr::select(individual_id, Species, treatment) %>% 
              distinct(), 
            by = "individual_id")



# Perform an anti-join to identify missing tracking dates.
missing_dates <- 
  all_combinations %>%
  anti_join(muddyfoot_filt_data, by = c("individual_id", "Date")) %>%
  arrange(individual_id, Date)

# Create a summary of missing dates for each individual.
summary_missing_dates <- missing_dates %>%
  group_by(individual_id) %>%
  summarise(Species = first(Species),
            Treatment = first(treatment),
            missing_dates = paste(Date, collapse = ", "),  # List missing dates.
            n = n(),  # Count of missing dates.
            n_breaks = if_else(n() == 1, 0, sum(diff(Date) != 1)))  # Number of breaks in consecutive dates.

# Create a table of missing dates using flextable for visualization.
missing_dates_table <- 
  flextable(summary_missing_dates) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_id" = 'Fish ID', 
                    "missing_dates" = 'Dates fish was not tracked') %>% 
  width(j = 4, 11, unit = 'cm')

# Optional: Save the missing dates table as a Word document.
save_as_docx(missing_dates_table, path = paste0(save_tables_path, "muddyfoot_fish_w_missing_dates.docx"))


#>>> 4.2 Add information about missing dates to positions summary ####

# Combine missing dates information with the positions summary.
positions_sum <- 
  positions_sum %>% 
  left_join(summary_missing_dates[, c('individual_id', 'n', 'n_breaks')], by = "individual_id") %>%
  mutate(n_missing_dates = ifelse(is.na(n), 0, n),  # If no missing dates, set to 0.
         n_breaks = ifelse(is.na(n_breaks), 0, n_breaks)) %>% 
  mutate_if(is.numeric, round, 1)  # Round all numeric columns.

# Create a summary table of positions using flextable.
positions_sum_table <- 
  flextable(positions_sum) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_id" = 'Fish ID',
                    "treatment" = 'Treatment',
                    "n_positions" = 'Number of locations',
                    "n_days_mon" = 'Number of days tracked',
                    "mean_time_diff" = 'Mean location frequency (s)',
                    "effective_n" = "Effective sample size",
                    "n_breaks" = "Date sequence breaks in tracking",
                    "n_missing_dates" = "Number of untracked days") 

# Optional: Save the positions summary table as a Word document.
save_as_docx(positions_sum_table, path = paste0(save_tables_path, "muddyfoot_fish_locations_sum.docx"))

# Update the main dataset by adding effective sample size and missing dates.
muddyfoot_filt_data <- merge(muddyfoot_filt_data, 
                             positions_sum[, c("individual_id", "effective_n", "n_breaks", "n_missing_dates")],
                             by = "individual_id", 
                             all.x = TRUE)


# >>> 4.3 Extract individual IDs that had irregularities in tracking ####

# This section identifies individuals with irregular tracking patterns, possibly due to predation, 
# equipment failure, or battery issues. These individuals may need adjusted home range estimates.

# Summarize the range of locations (n_positions) to determine thresholds for irregular sampling.
positions_sum %>% 
  dplyr::select(n_positions) %>% 
  summarise(mean = mean(n_positions, na.rm = TRUE),
            median = median(n_positions, na.rm = TRUE),
            stdev = sd(n_positions, na.rm = TRUE),
            min = min(n_positions, na.rm = TRUE),
            max = max(n_positions, na.rm = TRUE))

# Set a threshold for irregular sampling as -1 standard deviation from the mean number of positions.
# This threshold will be used to identify individuals with fewer than 58,574 locations.
irreg_indiv <- positions_sum %>% 
  filter(n_positions < 58574.04 | n_missing_dates >= 10 | n_breaks >= 2)

# Add a flag for individuals with irregular sampling in the main dataset.
irregular_ids <- irreg_indiv$individual_id

muddyfoot_filt_data <- 
  muddyfoot_filt_data %>% 
  mutate(irregular_sampling = ifelse(individual_id %in% irregular_ids, 1, 0))

### Identify irregular days ###

# Calculate summary statistics (mean, standard deviation) for daily time_diff and positions per species.
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
  dplyr::select(Species, individual_id, Date, median_day_time_diff, n_positions_day, n_positions_hourly) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  filter((n_positions_day <= 1000)) %>%  # Remove individuals with very low positions.
  filter((median_day_time_diff > 5)) %>%  # Remove very short time gaps.
  filter((n_positions_hourly < 30)) %>% 
  distinct(individual_id, Date, .keep_all = TRUE) %>% 
  ungroup()

# Flag irregular sampling days in the main dataset.
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  mutate(irreg_sample_day = ifelse(paste(individual_id, Date) %in% paste(daily_pos_irreg_ind$individual_id, 
                                                                         daily_pos_irreg_ind$Date), 1, 0))

# Check the flagged data for irregular days.
check <- print(
  muddyfoot_filt_data %>% 
    filter(irreg_sample_day == 1) %>% 
    dplyr::select(individual_id, Date, n_positions_day, median_day_time_diff, n_positions_hourly) %>% 
    mutate_if(is.numeric, round, 1) %>% 
    distinct()
)

# Create a table of identified irregular sampling days using flextable.
individual_irregular_sample_table <- 
  flextable(check) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_id" = 'Fish ID', 
                    "Date" = "Date",
                    "n_positions_day" = "Number of positions",
                    "median_day_time_diff" = "Median location frequency (s)",
                    "n_positions_hourly" = "Hourly positions") %>% 
  width(j = c(3:5), width = 4, unit = 'cm')

# Optional: Save the irregular sampling table as a Word document.
save_as_docx(individual_irregular_sample_table, path = paste0(save_tables_path, "muddyfoot_individual_irregular_sampling.docx"))


#----------------------------------------------------------------------#
# >>> 4.4. Summarize irregular days and missing data per individual ####
#-----------------------------------------------------------------------#

# Reduce the data to one row per individual per date, summarizing irregular and missing days.
reduced_data <- muddyfoot_filt_data %>%
  group_by(Species, individual_id, Date) %>%
  summarise(
    irreg_sample_day = max(irreg_sample_day),  # Maximum indicates any irregular sampling.
    n_missing_dates = max(n_missing_dates),  # Maximum for missing dates.
    .groups = 'drop'
  )

# Summarize the number of irregular and missing days per individual.
summary_data <- reduced_data %>%
  group_by(Species, individual_id) %>%
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
  dplyr::select(individual_id, missing_dates)

summary_data_test <- 
  summary_data %>% 
  left_join(missing_dates_sub, by = 'individual_id')

# Create a flextable summarizing the irregular and missing days.
daily_irregular_sample_table <- 
  flextable(summary_data_test) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_id" = 'Fish ID', 
                    "n_irregular_days" = "Number of days with irregular sampling",
                    "total_missing_days" = "Number of days not tracked",
                    "total_irregular_or_missing_days" = "Total poor tracking days",
                    "dates_irregular_positions" = "Dates with irregular sampling",
                    "missing_dates" = 'Dates fish were not tracked') %>% 
  width(j = c(6,7), 13, unit = 'cm')

# Optional: Save the daily irregular sample table as a Word document.
save_as_docx(daily_irregular_sample_table, 
             path = paste0(save_tables_path, "muddyfoot_irregular_sampling.docx"))

#----------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------#
# 5. Create a summary table for each species and treatment ####
#---------------------------------------------------------#

# Create a summary table showing key metrics for each species and treatment combination.

# Set the species factor order for the summary table.
muddyfoot_filt_data$Species <- factor(muddyfoot_filt_data$Species, levels = c("Roach", "Perch", "Pike"))

# Summarize median frequency, positions, tracking duration, and effective sample size.
summary_table <- 
  muddyfoot_filt_data %>%
  group_by(Species, treatment) %>%
  summarise(
    median_freq = median(time_diff, na.rm = TRUE),  # Median time between positions.
    mean_positions = mean(n_positions, na.rm = TRUE),  # Average number of positions.
    min_positions = min(n_positions, na.rm = TRUE),  # Minimum number of positions.
    max_positions = max(n_positions, na.rm = TRUE),  # Maximum number of positions.
    avg_days_mon = mean(n_days_mon, na.rm = TRUE),  # Average number of days tracked.
    min_days = min(n_days_mon, na.rm = TRUE),  # Minimum days tracked.
    max_days = max(n_days_mon, na.rm = TRUE),  # Maximum days tracked.
    mean_eff_n = mean(effective_n, na.rm = TRUE),  # Average effective sample size.
    min_eff_n = min(effective_n, na.rm = TRUE),  # Minimum effective sample size.
    max_eff_n = max(effective_n, na.rm = TRUE)  # Maximum effective sample size.
  ) %>% 
  mutate_if(is.numeric, round, 0)  # Round all numeric values to whole numbers.

# Format the summary table by combining metrics into readable strings.
summary_table <- 
  summary_table %>% 
  mutate(median_freq = median_freq,
         positions = paste(mean_positions, "(", min_positions, "-", max_positions, ")", sep = ""),
         days_mon = paste(avg_days_mon, "(", min_days, "-", max_days, ")", sep = ""),
         effective_n = paste(mean_eff_n, "(", min_eff_n, "-", max_eff_n, ")", sep = "")) %>% 
  dplyr::select(Species, treatment, median_freq, positions, days_mon, effective_n)

# Create the summary table using flextable for easy viewing and export.
muddyfoot_species_positions_sum <- 
  flextable(summary_table) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("individual_id" = 'Fish ID',
                    "treatment" = 'Treatment',
                    "median_freq" = 'Frequency (s)',
                    "positions" = 'Locations',
                    "days_mon" = 'Duration (days)',
                    "effective_n" = 'Effective sample size') %>% 
  width(width = 3, unit = 'cm') %>% 
  width(j = c(1,2), width = 2, unit = 'cm')

# Save the final summary table as a Word document.
save_as_docx(muddyfoot_species_positions_sum, 
             path = paste0(save_tables_path, "muddyfoot_species_location_summary_MS.docx"))

# Optional: Save the filtered dataset for future use.
saveRDS(muddyfoot_filt_data, paste0(data_filter_path, "04_muddyfoot_sub.rds"))
saveRDS(positions_sum, paste0(data_filter_path, "muddyfoot_daily_location_sum.rds"))

