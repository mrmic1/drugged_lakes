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
# saveRDS(muddyfoot_filt_data, paste0(data_filter_path, "04_muddyfoot_sub.rds"))
# saveRDS(positions_sum, paste0(data_filter_path, "muddyfoot_daily_location_sum.rds"))


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


# > 3.4 Extract effective sample sizes for akde estimation ####

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