################################################################################
### TRACKING DATA QUALITY ANALYSIS - LAKE MUDDYFOOT
### Purpose: Extract and analyze tracking summary statistics, identify data
###          quality issues, and flag irregular sampling patterns
### Author: Marcus Michelangeli
################################################################################

# =============================================================================
# SETUP
# =============================================================================

## 1.1 Load Required Libraries ----
# Core data manipulation
library(tidyverse)      # Data wrangling and visualization
library(data.table)     # Fast data manipulation
library(lubridate)      # Date-time handling

# Spatial analysis
library(sp)             # Spatial data structures
library(sf)             # Simple features for spatial vectors
library(ctmm)           # Continuous-time movement models

# Table formatting and export
library(flextable)      # Professional table creation
library(kableExtra)     # Additional table formatting
library(officer)        # Export to Word documents

## 1.2 Configure Table Defaults ----
set_flextable_defaults(
  font.color = "black",
  border.color = "black",
  font.family = "Arial",
  line_spacing = 1
)

## 1.3 Define File Paths ----
# Input directories
ctmm_path <- "./data/ctmm_fits/"
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
telem_path <- "./data/telem_obj/muddyfoot/"
enc_path <- "./data/encounters/"
size_path <- "./data/fish_size/"

# Output directories
save_tables_path <- "./tables/muddyfoot/"
save_figures_path <- "./figures/muddyfoot/"

## 1.4 Load Data ----
muddyfoot_filt_data <- readRDS(paste0(filtered_data_path, "03_muddyfoot_sub.rds"))

## 1.5 Define Study Period ----
date_range <- range(muddyfoot_filt_data$Date, na.rm = TRUE)
number_of_days <- as.integer(difftime(date_range[2], date_range[1], units = "days")) + 1

cat("Study Period:\n")
cat("  Start:", as.character(date_range[1]), "\n")
cat("  End:", as.character(date_range[2]), "\n")
cat("  Total Days:", number_of_days, "\n\n") #36 days

# =============================================================================
# SECTION 1: CALCULATE BASIC TRACKING METRICS
# =============================================================================

## 1.1 Time Intervals Between Positions ----
# Calculate time difference (in seconds) between consecutive GPS fixes
muddyfoot_filt_data <- muddyfoot_filt_data %>%
  group_by(individual_ID) %>%
  mutate(
    time_diff = c(NA, diff(timestamp)),  # First fix has no previous fix
    time_diff = as.numeric(round(time_diff, digits = 3))
  ) %>%
  ungroup()

# Preview time intervals
cat("Sample of time intervals between positions:\n")
print(head(muddyfoot_filt_data %>% select(timestamp, time_diff), n = 20))

## 1.2 Summary Statistics per Individual ----
# Calculate mean and median time intervals for each fish
muddyfoot_filt_data <- muddyfoot_filt_data %>%
  group_by(individual_ID) %>%
  mutate(
    mean_time_diff = mean(time_diff, na.rm = TRUE),
    median_time_diff = median(time_diff, na.rm = TRUE),
    n_positions = n()  # Total number of GPS fixes
  ) %>%
  ungroup()

cat("\nIndividual tracking summary:\n")
print(
  muddyfoot_filt_data %>%
    select(individual_ID, n_positions, mean_time_diff, median_time_diff) %>%
    distinct() %>%
    head(10)
)

## 1.3 Calculate Tracking Duration ----
# Determine how many unique days each fish was tracked
muddyfoot_filt_data <- muddyfoot_filt_data %>%
  group_by(individual_ID) %>%
  mutate(n_days_tracked = n_distinct(Date)) %>%
  ungroup()

## 1.4 Daily Position Metrics ----
# Calculate positions per day and location frequency
muddyfoot_filt_data <- muddyfoot_filt_data %>%
  group_by(individual_ID, Date) %>%
  mutate(
    n_positions_day = n(),
    median_day_time_diff = median(time_diff, na.rm = TRUE),
    n_positions_hourly = n_positions_day / 24,
    n_positions_per_min = n_positions_hourly / 60
  ) %>%
  ungroup()


# =============================================================================
# SECTION 2: IDENTIFY MISSING TRACKING DAYS
# =============================================================================

## 2.1 Create Complete Date Sequence ----
# Generate all possible individual-date combinations for the study period
full_date_range <- seq(as.Date(date_range[1]), as.Date(date_range[2]), by = "day")

all_combinations <- expand.grid(
  individual_ID = unique(muddyfoot_filt_data$individual_ID),
  Date = full_date_range,
  stringsAsFactors = FALSE
)

# Add species and treatment information
all_combinations <- all_combinations %>%
  left_join(
    muddyfoot_filt_data %>%
      select(individual_ID, Species, Treatment) %>%
      distinct(),
    by = "individual_ID"
  )

## 2.2 Find Missing Dates ----
# Identify dates where individuals were not tracked at all
missing_dates <- all_combinations %>%
  anti_join(muddyfoot_filt_data, by = c("individual_ID", "Date")) %>%
  arrange(individual_ID, Date)

## 2.3 Summarize Missing Data Patterns ----
summary_missing_dates <- missing_dates %>%
  group_by(individual_ID) %>%
  summarise(
    Species = first(Species),
    Treatment = first(Treatment),
    missing_dates = paste(Date, collapse = ", "),
    n = n(),  # Total missing days
    n_breaks = if_else(
      n() == 1,
      0,
      sum(diff(Date) != 1)  # Number of gaps in tracking
    ),
    .groups = "drop"
  )

## 2.4 Create Missing Dates Table ----
missing_dates_table <- flextable(summary_missing_dates) %>%
  fontsize(part = "all", size = 11) %>%
  bold(part = "header") %>%
  set_header_labels(
    individual_ID = "ID",
    missing_dates = "Dates ID was not tracked"
  ) %>%
  width(j = 4, width = 11, unit = "cm")

# Export table
save_as_docx(
  missing_dates_table,
  path = paste0(save_tables_path, "muddyfoot_IDs_missing_dates.docx")
)

## 2.5 Merge Missing Data Info with Main Dataset ----
# Create summary of tracking quality metrics
positions_sum <- muddyfoot_filt_data %>%
  select(
    individual_ID, Treatment, Species,
    n_positions, n_days_tracked,
    mean_time_diff, median_time_diff
  ) %>%
  distinct()

# Add missing dates information
positions_sum <- positions_sum %>%
  left_join(
    summary_missing_dates %>% select(individual_ID, n, n_breaks),
    by = "individual_ID"
  ) %>%
  mutate(
    n_missing_dates = ifelse(is.na(n), 0, n),
    n_breaks = ifelse(is.na(n_breaks), 0, n_breaks)
  ) %>%
  select(-n) %>%
  mutate(across(where(is.numeric), ~ round(., 1)))

## 2.6 Create Comprehensive Positions Summary Table ----
positions_sum_table <- flextable(positions_sum %>% select(-mean_time_diff)) %>%
  fontsize(part = "all", size = 11) %>%
  bold(part = "header") %>%
  set_header_labels(
    individual_ID = "ID",
    Treatment = "Treatment",
    n_positions = "Number of locations",
    n_days_tracked = "Number of days tracked",
    median_time_diff = "Median location frequency (s)",
    n_breaks = "Date sequence breaks in tracking",
    n_missing_dates = "Number of untracked days",
    found_alive = "Found alive (1=yes, 0=no)",
    known_predated = "Found predated (1=yes, 0=no)"
  )

# Export table
save_as_docx(
  positions_sum_table,
  path = paste0(save_tables_path, "muddyfoot_ID_posititions_sum_unfilt.docx")
)

## 2.8 Update Main Dataset ----
muddyfoot_filt_data <- muddyfoot_filt_data %>%
  left_join(
    positions_sum %>%
      select(individual_ID, n_missing_dates),
    by = "individual_ID"
  )


# =============================================================================
# SECTION 3: IDENTIFY DAYS WITH POOR TRACKING QUALITY
# =============================================================================

## 3.1 Create Daily Sampling Summary ----
# Aggregate position data by individual and date
daily_sampling_sum <- muddyfoot_filt_data %>%
  group_by(individual_ID, Date) %>%
  summarise(
    Species = first(Species),
    Treatment = first(Treatment),
    n_positions_day = n(),
    avg_time_diff = mean(time_diff, na.rm = TRUE),
    median_day_time_diff = first(median_day_time_diff),
    stdev_time_diff = sd(time_diff, na.rm = TRUE),
    min_time_diff = min(time_diff, na.rm = TRUE),
    max_time_diff = max(time_diff, na.rm = TRUE),
    n_days_tracked = first(n_days_tracked),
    .groups = "drop"
  )

## 3.2 Calculate Species-Level Position Statistics ----
# Overall statistics across all dates
species_day_positions_stats <- daily_sampling_sum %>%
  group_by(Species) %>%
  summarise(
    mean = mean(n_positions_day, na.rm = TRUE),
    median = median(n_positions_day, na.rm = TRUE),
    Q1_n_positions = quantile(n_positions_day, 0.25, na.rm = TRUE),
    Q3_n_positions = quantile(n_positions_day, 0.75, na.rm = TRUE),
    min = min(n_positions_day, na.rm = TRUE),
    max = max(n_positions_day, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nSpecies-level daily position statistics:\n")
print(species_positions_stats)

## 3.3 Visualize Daily Position Distribution ----
species_positions_histogram <- ggplot(
  daily_sampling_sum,
  aes(x = n_positions_day, fill = Species)
) +
  geom_histogram(
    binwidth = 500,
    position = "dodge",
    color = "black",
    fill = "lightgrey",
    alpha = 0.8
  ) +
  facet_wrap(~Species, scales = "free") +
  # Add reference lines
  geom_vline(
    data = species_positions_stats,
    aes(xintercept = mean, color = "Mean"),
    linetype = "solid",
    linewidth = 1.5
  ) +
  geom_vline(
    data = species_positions_stats,
    aes(xintercept = median, color = "Median"),
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_vline(
    data = species_positions_stats,
    aes(xintercept = Q1_n_positions, color = "25% IQR"),
    linetype = "dotted",
    linewidth = 1.5
  ) +
  geom_vline(
    data = species_positions_stats,
    aes(xintercept = Q3_n_positions, color = "75% IQR"),
    linetype = "dotted",
    linewidth = 1.5
  ) +
  scale_color_manual(
    values = c(
      "Mean" = "blue",
      "Median" = "red",
      "25% IQR" = "purple",
      "75% IQR" = "purple"
    ),
    name = "Statistics"
  ) +
  labs(
    x = "Number of positions per day",
    y = "Frequency"
  ) +
  coord_cartesian(xlim = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none")

print(species_positions_histogram)
ggsave(
  filename = "muddyfoot_species_daily_positions_histogram.jpg",
  plot = species_positions_histogram,
  path = save_figures_path,
  width = 210,        # A4 width in mm
  height = 99,        # ~1/3 of A4 height (297mm / 3)
  units = "mm",
  dpi = 300,          # High quality for publication
  device = "jpg"
)

## 3.4 Calculate Date-Specific Position Statistics ----
# Account for temporal variation (e.g., battery degradation, weather)
daily_species_positions_stats <- 
  daily_sampling_sum %>%
  group_by(Species, Date) %>%
  summarise(
    mean = mean(n_positions_day, na.rm = TRUE),
    median = median(n_positions_day, na.rm = TRUE),
    lower_25_quartile = quantile(n_positions_day, 0.25, na.rm = TRUE),
    .groups = "drop"
  )

## 3.5 Flag Poor Tracking Days ----
# Days with positions below the 25th percentile for that species/date
daily_sampling_sum <- daily_sampling_sum %>%
  left_join(
    daily_species_positions_stats %>%
      select(Date, Species, lower_25_quartile),
    by = c("Date", "Species")
  ) %>%
  mutate(
    poor_tracking_day = if_else(n_positions_day <= lower_25_quartile, 1, 0)
  )

## 3.6 Identify Individuals with Persistent Poor Tracking ----
daily_sampling_sum <- daily_sampling_sum %>%
  group_by(individual_ID) %>%
  mutate(
    # Count total number of poor tracking days
    n_poor_tracking_days = sum(poor_tracking_day, na.rm = TRUE),
    
    # Find longest consecutive streak of poor tracking
    max_consecutive_poor_tracking_days = {
      streaks <- rle(poor_tracking_day == 1)
      if (any(streaks$values)) {
        max(streaks$lengths[streaks$values == TRUE])
      } else {
        0
      }
    }
  ) %>%
  ungroup()

# Assess overall poor tracking prevalence
poor_tracking_summary <- daily_sampling_sum %>%
  summarise(
    median = median(n_poor_tracking_days, na.rm = TRUE),
    mean = mean(n_poor_tracking_days, na.rm = TRUE)
  )

cat("\nPoor tracking days summary:\n")
print(poor_tracking_summary)

## 3.7 Extract Irregular Individuals ----
# Threshold based on median + buffer
POOR_TRACKING_THRESHOLD <- 9  # Days with poor tracking

irreg_indiv <- daily_sampling_sum %>%
  filter(poor_tracking_day == 1 & n_poor_tracking_days > POOR_TRACKING_THRESHOLD)

cat(
  "\nIndividuals with >", POOR_TRACKING_THRESHOLD,
  "poor tracking days:", n_distinct(irreg_indiv$individual_ID), "\n"
)

## 3.8 Flag Poor Tracking Days in Main Dataset ----
muddyfoot_filt_data <- muddyfoot_filt_data %>%
  mutate(
    poor_tracking_day = if_else(
      paste(individual_ID, Date) %in% paste(irreg_indiv$individual_ID, irreg_indiv$Date),
      1,
      0
    )
  )

# =============================================================================
# SECTION 4: COMPREHENSIVE INDIVIDUAL-LEVEL SUMMARY
# =============================================================================

## 4.1 Aggregate Data Quality Issues ----
# Reduce to one row per individual per date
reduced_data <- 
  muddyfoot_filt_data %>%
  group_by(Species, individual_ID, Date) %>%
  summarise(
    poor_tracking_day = max(poor_tracking_day),
    n_missing_dates = max(n_missing_dates),
    .groups = "drop"
  )

## 4.2 Create Individual-Level Summary ----
summary_data <- reduced_data %>%
  group_by(Species, individual_ID) %>%
  summarise(
    n_poor_tracking_days = sum(poor_tracking_day == 1),
    total_missing_days = max(n_missing_dates),
    total_irregular_or_missing_days = n_poor_tracking_days + total_missing_days,
    dates_irregular_positions = paste(Date[poor_tracking_day == 1], collapse = ", "),
    .groups = "drop"
  )

## 4.3 Add Missing Dates Details ----
missing_dates_sub <- summary_missing_dates %>%
  select(individual_ID, missing_dates)

summary_data <- summary_data %>%
  left_join(missing_dates_sub, by = "individual_ID")

## 4.4 Add Survival Information ----
post_size_cols <- muddyfoot_filt_data %>%
  select(individual_ID, found_alive, known_predated) %>%
  distinct()

summary_data <- summary_data %>%
  left_join(post_size_cols, by = "individual_ID")

## 4.5 Create Final Summary Table ----
daily_irregular_sample_table <- flextable(summary_data) %>%
  fontsize(part = "all", size = 11) %>%
  bold(part = "header") %>%
  set_header_labels(
    Species = "Species",
    individual_ID = "ID",
    n_poor_tracking_days = "Number of poor tracking days",
    total_missing_days = "Number of days not tracked",
    total_irregular_or_missing_days = "Total number of irregular tracking days",
    dates_irregular_positions = "Dates with poor tracking",
    missing_dates = "Dates ID was not tracked",
    found_alive = "Found alive (1=yes, 0=no)",
    known_predated = "Found predated (1=yes, 0=no)"
  ) %>%
  width(j = c(6, 7), width = 13, unit = "cm")

# Export table
save_as_docx(
  daily_irregular_sample_table,
  path = paste0(save_tables_path, "muddyfoot_ID_irregular_tracking_summary.docx")
)

################################################################################
# END OF SCRIPT
################################################################################