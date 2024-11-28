#---------------------------------------------------#
### EXTRACTING TRACKING SUMMARY STATISTICS - LAKE BT
#---------------------------------------------------#

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
library(data.table)


#SET TABLE PLOTTING PARAMETERS
set_flextable_defaults(
  font.color = "black",
  border.color = "black",
  font.family = 'Arial',
  line_spacing = 1
)

### DIRECTORIES ###
# Define paths to various datasets and save locations for tables.
ctmm_path <- "./data/ctmm_fits/"  # Path for ctmm model fits.
filtered_data_path <- "./data/tracks_filtered/lake_BT/"  # Path for filtered tracking data.
telem_path <- "./data/telem_obj/BT/"  # Path for telemetry object files.
save_tables_path <- "./tables/BT/"  # Path to save summary tables.
enc_path <- "./data/encounters/"
size_path <- "./data/fish_size/"


#---------------------##
#### Load data ####
#----------------------#

### Load pike lake_BT dataset and ctmm ###
BT_filt_data <- readRDS(paste0(filtered_data_path, '02_lake_BT_sub.rds'))

#--------------------------------------------------------#
# 1. Calculate tracking summary metrics ##################
#--------------------------------------------------------#

#Below we calculate the following metrics
#- mean and median time difference between timestamps
#- number of locations per individual
#- number of days individuals were tracked for.
#- number of locations per day, hour and min
#- effective sample size for each individual used for calculating akdes
#- daily location information

#tracking date range
date_range <- range(BT_filt_data$Date, na.rm = FALSE)
date_range
#"2022-09-26" "2022-10-29"

number_of_days <- as.integer(difftime(date_range[2], date_range[1], units = "days")) + 1
number_of_days
#34

# > 1.1 Calculate time difference between positions ####

# Calculate time differences between successive timestamps for each individual.
BT_filt_data <- 
  BT_filt_data %>%
  group_by(individual_ID) %>% 
  mutate(time_diff = c(NA, diff(timestamp)))  # First difference is NA.

# Round the time differences to 3 decimal places for clarity.
BT_filt_data$time_diff <- as.numeric(round(BT_filt_data$time_diff, digits = 3))

# Check the first 20 rows to ensure time differences are calculated correctly.
head(BT_filt_data %>% dplyr::select(timestamp, time_diff), n = 20)

# Calculate mean and median time differences per individual.
BT_filt_data <- 
  BT_filt_data %>%
  group_by(individual_ID) %>% 
  mutate(mean_time_diff = mean(time_diff, na.rm = TRUE),
         median_time_diff = median(time_diff, na.rm = TRUE)) %>% 
  ungroup()

# > 1.2 Calculate number of positions per individual ####

# Count the number of positions (i.e., tracking locations) for each individual.
BT_filt_data <- 
  BT_filt_data %>%
  group_by(individual_ID) %>% 
  mutate(n_positions = n()) %>% 
  ungroup()

# Print a summary of unique individuals and their corresponding number of positions.
print(BT_filt_data %>% 
        dplyr::select(individual_ID, n_positions) %>% 
        distinct(), 
      n = 66)


# > 1.3 Calculate number of days individuals were tracked ####

# Calculate the number of unique days each individual was monitored.
BT_filt_data <- 
  BT_filt_data %>%
  group_by(individual_ID) %>% 
  mutate(n_days_tracked = length(unique(Date))) %>% 
  ungroup()

# Optional: Check the summary of number of days and positions per individual.
print(BT_filt_data %>% 
        dplyr::select(individual_ID, Treatment, n_positions, n_days_tracked) %>% 
        distinct(), 
      n = 65)


#> 1.4. Calculate daily positions sum stats ####

# This section calculates daily tracking metrics for each individual, 
# including the number of positions per day and time differences between positions.

# Calculate the number of positions per day and median time differences between positions for each individual.
BT_filt_data <- 
  BT_filt_data %>%
  group_by(individual_ID, Date) %>% 
  mutate(n_positions_day = n(),
         median_day_time_diff = median(time_diff, na.rm = TRUE)) %>%  # Median time difference per day.
  ungroup()

# Calculate the number of positions per hour and per minute.
BT_filt_data <- 
  BT_filt_data %>%
  group_by(individual_ID, Date) %>% 
  mutate(n_positions_hourly = n_positions_day / 24,  # Average number of positions per hour.
         n_positions_per_min = n_positions_hourly / 60) %>%  # Average number of positions per minute.
  ungroup()


#----------------------------------------------------#
# 2. Find dates that individuals were not tracked ####
#----------------------------------------------------#

# Identify dates where an individual was not tracked at all by comparing against a full date range.
full_date_range <- seq(as.Date("2022-09-26"), as.Date("2022-10-29"), by = "day")  # Define full date range.

# Create a table of all possible combinations of individual_id and dates.
all_combinations <- 
  expand.grid(
    individual_ID = unique(BT_filt_data$individual_ID),
    Date = full_date_range
  )

# Merge with original tracking data to include species and treatment information.
all_combinations <- 
  all_combinations %>%
  left_join(BT_filt_data %>% 
              dplyr::select(individual_ID, Species, Treatment) %>% 
              distinct(), 
            by = "individual_ID")


# Perform an anti-join to identify missing tracking dates.
missing_dates <- 
  all_combinations %>%
  anti_join(BT_filt_data, by = c("individual_ID", "Date")) %>%
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
save_as_docx(missing_dates_table, path = paste0(save_tables_path, "BT_IDs_missing_dates.docx"))


#> 2.1. Add information about missing dates to positions summary ####

# Create a summary table of the calculated metrics.
positions_sum <- 
  BT_filt_data %>% 
  dplyr::select(individual_ID, Treatment, Species,
                n_positions, 
                n_days_tracked, 
                mean_time_diff,
                median_time_diff) %>% 
  distinct()

# Combine missing dates information with the positions summary.
positions_sum <- 
  positions_sum %>% 
  left_join(summary_missing_dates[, c('individual_ID', 'n', 'n_breaks')], by = "individual_ID") %>%
  mutate(n_missing_dates = ifelse(is.na(n), 0, n),  # If no missing dates, set to 0.
         n_breaks = ifelse(is.na(n_breaks), 0, n_breaks)) %>% 
  mutate_if(is.numeric, round, 1)  # Round all numeric columns.


#Add information as whether individual was found dead or alive at the end of the experiment
#Load post-experiment biometric data to assess known deaths or survivals

post_biometrics <- 
  fread(paste0(size_path, "biometric_post_exp_data.csv")) %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number)))

# Merge biometric data (e.g., whether the fish was found alive) with encounter summary
post_size_cols <- 
  post_biometrics %>%
  filter(Lake == 'BT') %>%
  dplyr::select(individual_ID, Found, Known_predated) %>% 
  rename(found_alive = Found)


positions_sum <- 
  positions_sum %>%
  left_join(post_size_cols, by = c("individual_ID" = "individual_ID"))


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
                    "n_missing_dates" = "Number of untracked days",
                    "found_alive" = "Found alive (yes = 1, no = 0)",
                    "Known_predated" = "Found predated (yes = 1, no = 0)") 

# Optional: Save the positions summary table as a Word document.
save_as_docx(positions_sum_table, path = paste0(save_tables_path, "BT_ID_loc_sum_pred_events_unfilt.docx"))

# Update the main dataset by adding missing dates and survival information.
BT_filt_data <- merge(BT_filt_data, 
                       positions_sum[, c("individual_ID", "n_breaks", "n_missing_dates", "found_alive", "Known_predated")],
                       by = "individual_ID", 
                       all.x = TRUE)



#------------------------------------------------------#
# 3. Find dates that individuals had poor tracking ####
#------------------------------------------------------#

# This section identifies individuals with irregular tracking patterns, possibly due to predation, 
# equipment failure, or battery issues. These individuals may need to have adjusted home range estimates
# or we may need to re-run their ctmms depending on the identified issues.

# > 3.1.Summarise the daily positional data for each individual ####

daily_sampling_sum <- 
  BT_filt_data %>% 
  group_by(individual_ID, Date) %>%
  summarise(
    Species = first(Species),
    Treatment = first(Treatment),
    n_positions_day = n(), 
    avg_time_diff = mean(time_diff, na.rm = TRUE),
    median_time_diff = median(time_diff, na.rm = TRUE),
    stdev_time_diff = sd(time_diff, na.rm = TRUE),
    min_time_diff = min(time_diff, na.rm = TRUE),
    max_time_diff = max(time_diff, na.rm = TRUE),
    n_days_tracked = first(n_days_tracked)  # Number of days monitored.
  )



# > 3.2. Summary statistics of the positional data for each species

#For each day what is the mean, median, IQR, min and max number of locations
species_positions_stats <- 
  daily_sampling_sum %>%
  dplyr::select(Species, n_positions_day) %>% 
  group_by(Species) %>% 
  summarise(mean = mean(n_positions_day, na.rm = TRUE),
            median = median(n_positions_day, na.rm = TRUE),
            Q1_n_positions = quantile(n_positions_day, 0.25, na.rm = TRUE),
            Q3_n_positions = quantile(n_positions_day, 0.75, na.rm = TRUE),
            min = min(n_positions_day, na.rm = TRUE),
            max = max(n_positions_day, na.rm = TRUE),
            )


#Histogram of daily locations for each Species
species_positions_histogram =
  ggplot(daily_sampling_sum, aes(x = n_positions_day, fill = Species)) +
  geom_histogram(binwidth = 500,
                 position = "dodge", 
                 color = "black",
                 fill = 'lightgrey',
                 alpha = 0.8) +
  facet_wrap(~Species, scales = "free_y") +
  geom_vline(data = species_positions_stats, 
             aes(xintercept = mean, color = "Mean"), 
             linetype = "solid", 
             size = 1.5) +
  geom_vline(data = species_positions_stats, 
             aes(xintercept = median, color = "Median"), 
             linetype = "dashed", 
             size = 1) +
  geom_vline(data = species_positions_stats, 
             aes(xintercept = Q1_n_positions, color = "25% IQR"), 
             linetype = "dotted", 
             size = 1.5) +
  geom_vline(data = species_positions_stats, 
             aes(xintercept = Q3_n_positions, color = "75% IQR"), 
             linetype = "dotted", 
             size = 1.5) +
  # Add text labels
  geom_text(
    data = species_positions_stats,
    aes(x = mean, y = 10, label = "Mean"), color = "black", angle = 90, vjust = -0.5, hjust = 0
  ) +
  geom_text(
    data = species_positions_stats,
    aes(x = median, y = 10, label = "Median"), color = "black", angle = 90, vjust = -0.5, hjust = 0
  ) +
  geom_text(
    data = species_positions_stats,
    aes(x = Q1_n_positions, y = 10, label = "Q1"), color = "black", angle = 90, vjust = -0.5, hjust = 0
  ) +
  geom_text(
    data = species_positions_stats,
    aes(x = Q3_n_positions, y = 10, label = "Q3"), color = "black", angle = 90, vjust = -0.5, hjust = 0
  ) +
  scale_color_manual(
    values = c("Mean" = "blue", "Median" = "red", "25% IQR" = "purple", "75% IQR" = "purple"),
    name = "Statistics"
  ) +
  labs(
    x = "Number of Positions per Day",
    y = "Frequency"
  ) +
  coord_cartesian(xlim = c(0, NA)) + # Adjust x-axis dynamically
  facet_wrap(~Species, scales = "free") +
  theme_bw() +
  theme(legend.position = 'none')


#Species daily position stats
#Look at the positional stats for each species per date
#Expected number of positions may depend on the date, for example dates towards
#the end of the study may have fewer locations due to battery life
#all days with poor weather or receiver malfunction will have fewer locations.
#extract for each species and date the mean, median and lower 25% positions for that date

daily_species_positions_stats <- 
  daily_sampling_sum %>%
  dplyr::select(Species, Date, n_positions_day) %>% 
  group_by(Species, Date) %>% 
  summarise(mean = mean(n_positions_day, na.rm = TRUE),
            median = median(n_positions_day, na.rm = TRUE),
            lower_25_quartile = quantile(n_positions_day, 0.25, na.rm = TRUE),
  )


# Bring in Q1_positions to daily sampling summary
daily_sampling_sum <- 
  daily_sampling_sum %>%
  left_join(daily_species_positions_stats %>%
              select(Date, Species, lower_25_quartile),
            by = c("Date", "Species"))

# Add a new column identifying individuals with n_positions <= lower_25_quartile
daily_sampling_sum <- 
  daily_sampling_sum %>%
  mutate(poor_tracking_day = ifelse(n_positions_day <= lower_25_quartile, 1, 0))

# Count the number of dates below Q1 and the maximum consecutive dates below Q1
daily_sampling_sum <- 
  daily_sampling_sum %>%
  group_by(individual_ID) %>% # Group by individual ID
  mutate(
    # Count the number of dates below Q1
    n_poor_tracking_days = sum(poor_tracking_day, na.rm = TRUE),
    
    # Calculate the maximum number of consecutive dates below Q1
    max_consecutive_poor_tracking_days = {
      # Identify streaks of 1s (equivalent to TRUE values)
      streaks <- rle(poor_tracking_day == 1)
      # Extract the maximum streak of 1s
      max_streak <- ifelse(any(streaks$values), max(streaks$lengths[streaks$values == TRUE]), 0)
      max_streak
    }
  ) %>%
  ungroup()


#what is the median number of poor tracking days
daily_sampling_sum %>% 
  summarise(median = median(n_poor_tracking_days, na.rm = TRUE),
            mean = mean(n_poor_tracking_days, na.rm = TRUE))


#extract irregular individuals
irreg_indiv <- daily_sampling_sum %>% 
  #n_poor_tracking day threshold is based on the median
  filter(poor_tracking_day == 1 & n_poor_tracking_days > 6)


# Flag irregular sampling days in the main dataset.
BT_filt_data <- 
  BT_filt_data %>%
  mutate(poor_tracking_day = ifelse(paste(individual_ID, Date) %in% paste(irreg_indiv$individual_ID, 
                                                                          irreg_indiv$Date), 1, 0))
# Check the flagged data for irregular days.
check <- print(
  BT_filt_data %>% 
    filter(poor_tracking_day == 1) %>% 
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
                    "n_positions_hourly" = "Number of locations per hour") %>% 
  width(j = c(2:5), width = 3, unit = 'cm')

# Optional: Save the irregular sampling table as a Word document.
save_as_docx(individual_irregular_sample_table, 
             path = paste0(save_tables_path, "BT_ID_day_locs_w_irregular_sampling.docx"))


#--------------------------------------------------------------------#
# 4. Summarise irregular sampling and missing data per individual ####
#--------------------------------------------------------------------#

# Reduce the data to one row per individual per date, summarising irregular and missing days.
reduced_data <- 
  BT_filt_data %>%
  group_by(Species, individual_ID, Date) %>%
  summarise(
    poor_tracking_day = max(poor_tracking_day), 
    n_missing_dates = max(n_missing_dates),  
    .groups = 'drop'
  )

# Summarise the number of irregular and missing days per individual.
summary_data <- 
  reduced_data %>%
  group_by(Species, individual_ID) %>%
  summarise(
    n_poor_tracking_days = sum(poor_tracking_day == 1),  # Count irregular days.
    total_missing_days = max(n_missing_dates),  # Count missing days.
    total_irregular_or_missing_days = n_poor_tracking_days + total_missing_days,  # Combine irregular and missing days.
    dates_irregular_positions = paste(Date[poor_tracking_day == 1], collapse = ", "),  # List irregular dates.
    .groups = 'drop'
  )

# Merge the missing dates information.
missing_dates_sub <- 
  summary_missing_dates %>% 
  dplyr::select(individual_ID, missing_dates)

summary_data <- 
  summary_data %>% 
  left_join(missing_dates_sub, by = 'individual_ID')

#add info on known survival
summary_data <- 
  summary_data %>%
  left_join(post_size_cols, by = c("individual_ID" = "individual_ID"))

# Create a flextable summarizing the irregular and missing days.
daily_irregular_sample_table <- 
  flextable(BT_ID_irreg_sampling) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("Species" = 'Species',
                    "individual_ID" = 'ID', 
                    "n_poor_tracking_days" = "Number of poor tracking days",
                    "total_missing_days" = "Number of days not tracked",
                    "total_irregular_or_missing_days" = "Total number of irregular tracking days",
                    "dates_irregular_positions" = "Dates with poor tracking",
                    "missing_dates" = 'Dates ID was not tracked',
                    "found_alive" = 'Found alive (yes = 1, no = 0)',
                    "Known_predated" = "Found predated (yes = 1, no = 0)") %>% 
  width(j = c(6,7), 13, unit = 'cm')

# Optional: Save the daily irregular sample table as a Word document.
save_as_docx(daily_irregular_sample_table, 
             path = paste0(save_tables_path, "BT_ID_poor_tracking_summary.docx"))


#--------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------#

#SAVE
saveRDS(summary_data, paste0(filtered_data_path, "BT_ID_irreg_sampling.rds"))
saveRDS(BT_filt_data, paste0(filtered_data_path, "03_lake_BT_sub.rds"))
saveRDS(daily_sampling_sum, paste0(filtered_data_path, "BT_daily_tracking_ID_stats"))

