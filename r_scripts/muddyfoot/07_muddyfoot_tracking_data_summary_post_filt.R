#----------------------------------------------------------#
### Create a summary table for each species and treatment ##
#----------------------------------------------------------#

# LIBRARIES
# Loading essential libraries for data manipulation, spatial analysis, and formatting tables.
library(tidyverse)  # For data manipulation and visualization.
library(sp)  # For handling spatial data.
library(sf)  # For handling spatial vector data using simple features.
library(flextable)  # For creating and customizing tables.
library(kableExtra)  # Additional table formatting options.
library(officer)  # For exporting tables to Word documents.
library(ctmm)  # For continuous-time movement modeling (ctmm) package.

### DIRECTORIES ###
# Define paths to various datasets and save locations for tables.
ctmm_path <- "./data/ctmm_fits/"  # Path for ctmm model fits.
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"  # Path for filtered tracking data.
telem_path <- "./data/telem_obj/muddyfoot/"  # Path for telemetry object files.
save_tables_path <- "./tables/muddyfoot/"  # Path to save summary tables.
enc_path <- "./data/encounters/muddyfoot/"
size_path <- "./data/fish_size/"
akde_path <- "./data/akdes/"      


### Load data ###
muddyfoot_filt_data <- readRDS(paste0(filtered_data_path, '04_muddyfoot_sub.rds'))

#Load akdes to extract effective sample sizes
pike_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/pike_akdes_cg_list.rds"))
perch_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))
roach_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))

#Load telemetry objects for pike, perch, and roach species in Muddyfoot lake
pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

# Check predation event data (pre-identified predation events) for individual with many missing dates
mud_pred_mort_events <- readRDS(paste0(enc_path, "muddyfoot_pred_encounter_summary_filtered.rds"))

# Remove Roach individuals that were predated based on predation events
roach_ids_remove <- mud_pred_mort_events %>%
  filter(Species == "Roach" & revised_suspected_mortality == 'poor_tracking') %>%    # Filter for Roach species
  pull(individual_ID)               # Extract IDs of predated Roach

# Update the Roach AKDE list and telemetry data, removing predated individuals
roach_akdes_cg_list <- roach_akdes_cg_list[!(names(roach_akdes_cg_list) %in% roach_ids_remove)]
roach_muddyfoot_tel <- roach_muddyfoot_tel[!(names(roach_muddyfoot_tel) %in% roach_ids_remove)]
muddyfoot_filt_data <- muddyfoot_filt_data %>% 
  filter(!individual_ID == roach_ids_remove)

saveRDS(roach_muddyfoot_tel, paste0(telem_path, "roach_muddyfoot_tel.rds"))

#------------------------------------------------------------#
# 1. Extract effective sample size information ###############
#------------------------------------------------------------#

#effective sample size
roach_akdes_cg_list$F59683$DOF.area[1]

# Function to extract ID anDOF.area# Function to extract ID and DOF.H from a given list
extract_dofh <- function(akdes_list) {
  tibble(
    individual_ID = names(akdes_list),
    effective_n = map_dbl(akdes_list, ~ .x$DOF.area[1])
  )
}

# Extract for each species
roach_df <- extract_dofh(roach_akdes_cg_list)
#check
head(roach_df)
summary(roach_akdes_cg_list$F59684)

perch_df <- extract_dofh(perch_akdes_cg_list)
pike_df <- extract_dofh(pike_akdes_cg_list)

# Combine into a single dataframe
akde_effective_n_df <- bind_rows(roach_df, perch_df, pike_df)

# View the result
print(akde_effective_n_df)

### Combine to full dataset ###

muddyfoot_filt_data <- merge(muddyfoot_filt_data, akde_effective_n_df, by = "individual_ID", all.x = TRUE)

#save dataframe modified dataframe
saveRDS(muddyfoot_filt_data, paste0(filtered_data_path, "05_muddyfoot_sub.rds"))

#--------------------------------------------------------#
# 2. Calculate tracking summary metrics ##################
#--------------------------------------------------------#

# Create a summary table showing key metrics for each species and treatment combination.

# Set the species factor order for the summary table.
muddyfoot_filt_data$Species <- factor(muddyfoot_filt_data$Species, levels = c("Roach", "Perch", "Pike"))

# Summarise median frequency, positions, tracking duration, and effective sample size
# for each specie and treatment

summary_table <- 
  muddyfoot_filt_data %>%
  group_by(Species, Treatment) %>%
  summarise(
    median_freq = median(time_diff, na.rm = TRUE),  # Median time between positions.
    mean_positions = mean(n_positions, na.rm = TRUE),  # Average number of positions.
    min_positions = min(n_positions, na.rm = TRUE),  # Minimum number of positions.
    max_positions = max(n_positions, na.rm = TRUE),  # Maximum number of positions.
    avg_days_mon = mean(n_days_tracked, na.rm = TRUE),  # Average number of days tracked.
    min_days = min(n_days_tracked, na.rm = TRUE),  # Minimum days tracked.
    max_days = max(n_days_tracked, na.rm = TRUE),  # Maximum days tracked.
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
  dplyr::select(Species, Treatment, median_freq, positions, days_mon, effective_n)

# Create the summary table using flextable for easy viewing and export.
muddyfoot_species_positions_sum <- 
  flextable(summary_table) %>% 
  fontsize(part = "all", size = 11) %>% 
  bold(part = 'header') %>% 
  set_header_labels("Species" = 'Species',
                    "Treatment" = 'Treatment',
                    "median_freq" = 'Location frequency (s)',
                    "positions" = 'Number of locations',
                    "days_mon" = 'Days tracked',
                    "effective_n" = 'Effective n') %>% 
  width(width = 3, unit = 'cm') %>% 
  width(j = c(1,2), width = 2, unit = 'cm')

# Save the final summary table as a Word document.
save_as_docx(muddyfoot_species_positions_sum, 
             path = paste0(save_tables_path, "muddyfoot_species_location_summary_MS.docx"))

# Optional: Save the filtered dataset for future use.
# saveRDS(muddyfoot_filt_data, paste0(data_filter_path, "04_muddyfoot_sub.rds"))
# saveRDS(positions_sum, paste0(data_filter_path, "muddyfoot_daily_location_sum.rds"))
