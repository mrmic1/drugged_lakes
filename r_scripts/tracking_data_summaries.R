#### EXTRACTING TRACKING SUMMARY STATISTICS

#LIBRARIES 
library(tidyverse)
library(sp)
library(sf)
library(flextable)
library(kableExtra)
library(officer)

#SET TABLE PLOTTING PARAMETERS
set_flextable_defaults(
  font.color = "black",
  border.color = "black",
  font.family = 'Arial',
  line_spacing = 1
)

### DIRECTORIES ###
ctmm_path = "./data/ctmm_fits/"
data_filter_path = "./data/tracks_filtered/"
telem_path = "./data/telem_obj/"
save_tables_path = "./data/tracks_filtered/sum_tables/" 

#### MUDDYFOOT ####

### Load pike muddyfoot dataset and ctmm ###

muddyfoot_sub <- readRDS(paste0(data_filter_path, 'muddyfoot_sub.rds'))

#isolate pike
pike_mud <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Northern Pike')

#telemetry objects
pike_muddyfoot_tel = readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel = readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel = readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

#ctmm lists
pike_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds")) 

#> 1. Make combined dataframe to extract summary statistics ####
# Create a function to extract the data and add the individual_id column
extract_data <- function(list_data) {
  data_combined <- do.call(rbind, lapply(names(list_data), function(id) {
    df <- list_data[[id]]
    df$individual_id <- id
    return(df)
  }))
  return(data_combined)
}

# Apply the function to each list and combine them
pike_data <- extract_data(pike_muddyfoot_tel)
perch_data <- extract_data(perch_muddyfoot_tel)
roach_data <- extract_data(roach_muddyfoot_tel)

# Combine all dataframes into one
muddyfoot_filt_data <- 
  rbind(
    cbind(pike_data, Species = 'Pike'), 
    cbind(perch_data, Species = 'Perch'),
    cbind(roach_data, Species = 'Roach'))

#saveRDS(muddyfoot_filt_data, paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))
muddyfoot_filt_data <-  readRDS(paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))
#this is the filtered data including the removal of outliers using the out() function in ctmm

#---------------------------------------------------------------------------------#

#> 2. Calculate time difference between positions ####
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(time_diff = c(NA, diff(timestamp)))

muddyfoot_filt_data$time_diff <- as.numeric(round(muddyfoot_filt_data$time_diff, digits = 3))

#check time_diff
head(muddyfoot_filt_data %>% 
       dplyr::select(timestamp, time_diff), n = 20)

#average time_diff per individual
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(mean_time_diff = mean(time_diff, na.rm = TRUE),
         median_time_diff = median(time_diff, na.rm = TRUE)) %>% 
  ungroup()


#> 3. Calculate number of positions per individual ####
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(n_positions = n()) %>% 
  ungroup()

print(muddyfoot_filt_data %>% 
       dplyr::select(individual_id,n_positions) %>% 
       distinct(), 
     n = 65)

#> 4. Calculate number of days individuals were tracked or monitored (n_days_mon) ####
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(n_days_mon = length(unique(date))) %>% 
  ungroup()

print(muddyfoot_filt_data %>% 
        dplyr::select(individual_id, treatment, n_positions, n_days_mon) %>% 
        distinct(), 
      n = 65)

#calculate number of positions per day
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id) %>% 
  mutate(n_positions_per_day = n_positions/n_days_mon,
         n_positions_per_hour = n_positions_per_day/24,
         n_positions_per_min = n_positions_per_hour/60) %>% 
  ungroup()

positions_sum = 
  muddyfoot_filt_data %>% 
  dplyr::select(individual_id, treatment, Species,
                n_positions, 
                n_days_mon, 
                n_positions_per_day,
                n_positions_per_hour,
                n_positions_per_min,
                mean_time_diff) %>% 
  distinct() 


#> 5. Calculate number of individuals for each species within each treatment ####
muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(Species, treatment) %>% 
  mutate(n_individuals = length(unique(individual_id))) %>% 
  ungroup()

print(muddyfoot_filt_data %>% 
        dplyr::select(Species, treatment, n_individuals) %>% 
        distinct(), 
      n = 6)

#> 6. Extract effective sample sizes for home range estimation ####

#first combine all ctmm fits
muddyfoot_all_ctmm_fits <- list(Pike_fits = pike_muddyfoot_ctmm_fits, 
                                Perch_fits = perch_muddyfoot_ctmm_fits, 
                                Roach_fits = roach_muddyfoot_ctmm_fits)


# Initialize an empty data frame to store the results
summary_df <- data.frame(
  species = character(),
  ID = character(),
  effective_n = numeric(),
  stringsAsFactors = FALSE
)

# List the species in your dataset
species_list <- c("Pike_fits", "Perch_fits", "Roach_fits")

# Loop through each species in the list
for (species in species_list) {
  # Get the list of individuals for the species
  individuals <- names(muddyfoot_all_ctmm_fits[[species]])
  
  # Loop through each individual in the species
  for (ID in individuals) {
    # Extract the summary object
    summary_obj <- summary(muddyfoot_all_ctmm_fits[[species]][[ID]])
    
    # Extract the effective_n value using the correct indexing
    effective_n <- summary_obj$DOF["area"]
    
    # Append the results to the summary data frame
    summary_df <- rbind(summary_df, data.frame(
      species = gsub("_fits", "", species),  # Remove the '_fits' suffix for species name
      ID = ID,
      effective_n = effective_n
    ))
  }
}

# Print the summary data frame
print(summary_df)

effective_n_sum = summary_df


#> 7. Combine effective sample sizes with individual positions data ####

positions_sum <- merge(positions_sum, 
              effective_n_sum[, c("ID", "effective_n")], 
              by.x = "individual_id", 
              by.y = "ID", 
              all_x = T)

#What is the correlation between effective sample size and the number of locations
cor(positions_sum$n_positions, positions_sum$effective_n, method = 'spearman')
#0.71 - pretty high

#> 8. Calculate daily positions sum stats ####

muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id, date) %>% 
  mutate(n_day_positions = n(),
         median_day_time_diff = median(time_diff, na.rm = T)) %>% 
  ungroup()

muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  group_by(individual_id, date) %>% 
  mutate(n_day_positions_per_hour = n_day_positions/24,
         n_day_positions_per_min = n_day_positions/60) %>% 
  ungroup()

#check
print(
  muddyfoot_filt_data %>% 
  filter(individual_id == 'F59704') %>% 
  select(individual_id, Species, date, n_day_positions, median_day_time_diff, n_day_positions_per_hour, n_day_positions_per_min) %>% 
  distinct(), n = 36)


#------------------------------------------------------------------------------------------------#

#> 9. Find periods of irregular sampling ####
#For each individual and unique date calculate the average, stdev, min and max time difference between positions

daily_sampling_sum = 
  muddyfoot_filt_data %>% 
  group_by(individual_id, treatment, date) %>%
  summarise(
    avg_time_diff = mean(time_diff, na.rm = TRUE),
    stdev_time_diff = sd(time_diff, na.rm = TRUE),
    min_time_diff = min(time_diff, na.rm = TRUE),
    max_time_diff = max(time_diff, na.rm = TRUE),
    n_days_mon = first(n_days_mon))


# >>> Dates that individuals were not tracked ####

#I also list the dates that an individual was tracked and was not tracked
# Define the full sequence of dates
full_date_range  <- seq(as.Date("2022-09-25"), as.Date("2022-10-30"), by = "day")

# Create a table of all possible individual_id and dates
all_combinations <- expand.grid(
  individual_id = unique(muddyfoot_filt_data$individual_id),
  Date = full_date_range
)

# Merge with the original data to include 'Species' and 'Treatment'
all_combinations <- all_combinations %>%
  left_join(muddyfoot_filt_data %>% 
              select(individual_id, Species, treatment) %>% 
              distinct(), 
            by = "individual_id")

#need to make column Date
muddyfoot_filt_data$Date = format(with_tz(ymd_hms(muddyfoot_filt_data$timestamp), 
                                          tzone = "Europe/Stockholm"), "
                                  %Y/%m/%d")

muddyfoot_filt_data$Date = as.Date(muddyfoot_filt_data$Date)

# Identify the missing dates by performing an anti-join with the original tracking data
missing_dates <- all_combinations %>%
  anti_join(muddyfoot_filt_data, by = c("individual_id", "Date")) %>%
  arrange(individual_id, Date)

# Create summary table
summary_missing_dates <- missing_dates %>%
  group_by(individual_id) %>%
  summarise(Species = first(Species),
            Treatment = first(treatment),
            missing_dates = paste(Date, collapse = ", "),
            n = n(),
            n_breaks = if_else(n() == 1, 0, sum(diff(Date) != 1)))


(missing_dates_table <- 
    flextable(summary_missing_dates) %>% 
    fontsize(part = "all", size = 11) %>% 
    bold(part = 'header')  %>% 
    set_header_labels("individual_id" = 'Fish ID', 
                      "missing_dates" = 'Dates fish was not tracked') %>% 
    width(j = 4, 11, unit = 'cm'))



#save table
#save_as_docx(missing_dates_table, 
#              path = paste0(save_tables_path, "muddyfoot_fish_w_missing_dates.docx"))


#>>> Add information about individual missing dates to positions summary ####

positions_sum <- 
  positions_sum %>% 
  left_join(summary_missing_dates[, c('individual_id', 'n', 'n_breaks')], by = "individual_id") %>%
  mutate(n_missing_dates = ifelse(is.na(n), 0, n),
         n_breaks = ifelse(is.na(n_breaks), 0, n_breaks)) %>%
  select(-n, -n_positions_per_day, -n_positions_per_hour, -n_positions_per_min) %>% 
  mutate_if(is.numeric, round, 1) 
  
#might need to add median tracking frequency.

(positions_sum_table <- 
    flextable(positions_sum) %>% 
    fontsize(part = "all", size = 11) %>% 
    bold(part = 'header')  %>% 
    set_header_labels("individual_id" = 'Fish ID',
                      "treatment" = 'Treatment',
                      "n_positions" = 'Number of locations',
                      "n_days_mon" = 'Number of days tracked',
                      "mean_time_diff" = 'Mean location frequency (s)',
                      "effective_n" = "Effective sample size",
                      "n_breaks" = "Date sequence breaks in tracking",
                      "n_missing_dates" = "Number of untracked days")) 

save_as_docx(positions_sum_table, 
                          path = paste0(save_tables_path, "muddyfoot_fish_locations_sum.docx"))

#29/65 individuals had some sort of missed tracking

muddyfoot_filt_data <- merge(muddyfoot_filt_data, 
              positions_sum[, c("individual_id", "effective_n", "n_breaks", "n_missing_dates")],
              by = "individual_id", 
              all.x = TRUE)

#------------------------------------------------------------------------------------------------------------------#

# > 10. Extract individual IDs that had irregularities in tracking ####

#These individuals likely died or had battery malfunctions. 
#If they died it could be because of predation - good candidates to assess this
#Irregularly sampled individuals may also need to have their home range estimates weighted.

#range of locations
positions_sum %>% 
  select(n_positions) %>% 
  summarise(mean = mean(n_positions, na.rm = TRUE),
            median = median(n_positions, na.rm = TRUE),
            stdev = sd(n_positions, na.rm = TRUE),
            min = min(n_positions, na.rm = TRUE),
            max = max(n_positions, na.rm = TRUE))

# -1 std from the mean
# Look for individuals that had less than 63852.95 locations

irreg_indiv = positions_sum %>% 
  filter(n_positions < 63852.95 | n_missing_dates >= 10 | n_breaks >= 2)

#add column to muddyfoot_filt_data flagging irregular individuals
irregular_ids = irreg_indiv$individual_id

muddyfoot_filt_data <- 
  muddyfoot_filt_data %>% 
  mutate(irregular_sampling = ifelse(individual_id %in% irregular_ids, 1, 0))

saveRDS(muddyfoot_filt_data, paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))


#Identify irregular days 

#What is the median daily time_diff and number of positions per species
daily_pos_irreg_ind = 
  muddyfoot_filt_data %>% 
  group_by(Species) %>% 
  mutate(avg_median_time_diff = mean(median_time_diff, na.rm = T),
         std_median_time_diff = sd(median_time_diff, na.rm = T),
         avg_daily_positions = mean(n_day_positions, na.rm = T),
         std_daily_positions = sd(n_day_positions, na.rm = T),
         avg_daily_hour_positions = mean(n_day_positions_per_hour, na.rm = T),
         std_daily_hour_positions = sd(n_day_positions_per_hour, na.rm = T)
  ) %>% 
  filter((median_day_time_diff <= (avg_median_time_diff - 3 * std_median_time_diff))|
           (n_day_positions <= (avg_daily_positions - 1 * std_daily_positions))|
           (n_day_positions_per_hour <= (avg_daily_hour_positions - 1 * std_daily_hour_positions))) %>%
  select(Species, individual_id, date, median_day_time_diff, n_day_positions, n_day_positions_per_hour) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  filter((n_day_positions <= 1000)) %>%
  filter((median_day_time_diff > 5)) %>%
  filter((n_day_positions_per_hour < 30)) %>%
  distinct(individual_id, date, .keep_all = T) %>% 
  ungroup()

muddyfoot_filt_data <- 
  muddyfoot_filt_data %>%
  mutate(irreg_sample_day = ifelse(paste(individual_id, date) %in% paste(daily_pos_irreg_ind$individual_id, 
                                                                         daily_pos_irreg_ind$date), 1, 0))

#check
check = print(
  muddyfoot_filt_data %>% 
    filter(irreg_sample_day == '1') %>% 
    select(individual_id, date, n_day_positions, median_day_time_diff, n_day_positions_per_hour) %>% 
    mutate_if(is.numeric, round, 1) %>% 
    distinct())

#saveRDS(muddyfoot_filt_data, paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))
#muddyfoot_filt_data <-  readRDS(paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))
#this is the filtered data including the removal of outliers using the out() function in ctmm


(individual_irregular_sample_table <- 
    flextable(check) %>% 
    fontsize(part = "all", size = 11) %>% 
    bold(part = 'header')  %>% 
    set_header_labels("individual_id" = 'Fish ID', 
                      "date" = "Date",
                      "n_day_positions" = "Number of positions",
                      "median_day_time_diff" = "Median location frequency (s)",
                      "n_day_positions_per_hour" = "Hourly positions") %>% 
    width(j = c(3:5), width = 4, unit = 'cm'))



#save table
save_as_docx(individual_irregular_sample_table, 
             path = paste0(save_tables_path, "muddyfoot_individual_irregular_sampling.docx"))



#how many days per individuals had irregular sampling or had no locations
# First, reduce the data to one row per individual per date
reduced_data <- muddyfoot_filt_data %>%
  group_by(Species,individual_id, date) %>%
  summarise(
    irreg_sample_day = max(irreg_sample_day),
    n_missing_dates = max(n_missing_dates),
    .groups = 'drop'
  )

# Then, summarize the number of irregular days and add the missing days information
summary_data <- reduced_data %>%
  group_by(Species, individual_id) %>%
  summarise(
    n_irregular_days = sum(irreg_sample_day == 1),
    total_missing_days = max(n_missing_dates),
    total_irregular_or_missing_days = n_irregular_days + total_missing_days,
    dates_irregular_positions  = paste(date[irreg_sample_day == 1], collapse = ", "),
    .groups = 'drop'
  )

#merge missing dates
missing_dates_sub <- 
  summary_missing_dates %>% 
  select(individual_id, missing_dates)

summary_data_test = 
  summary_data %>% 
  left_join(missing_dates_sub, by = 'individual_id')


(daily_irregular_sample_table <- 
    flextable(summary_data_test) %>% 
    fontsize(part = "all", size = 11) %>% 
    bold(part = 'header')  %>% 
    set_header_labels("individual_id" = 'Fish ID', 
                      "n_irregular_days" = "Number of days with irregular sampling",
                      "total_missing_days" = "Number of days not tracked",
                      "total_irregular_or_missing_days" = "Total poor tracking days",
                      "dates_irregular_positions" = "Dates with irregular sampling",
                      "missing_dates" = 'Dates fish were not tracked') %>% 
    width(j = c(6,7), 13, unit = 'cm'))



#save table
save_as_docx(daily_irregular_sample_table, 
              path = paste0(save_tables_path, "muddyfoot_irregular_sampling.docx"))


#----------------------------------------------------------------------------------------------------------------#

#Need to create a summary table for each species and treatment (i.e regarding sample size and location related statistics)

# Summary table for manuscript

#change species order for summary table
muddyfoot_filt_data$Species <- factor(muddyfoot_filt_data$Species, levels = c("Roach", "Perch", "Pike"))

summary_table <- 
  muddyfoot_filt_data %>%
  group_by(Species, treatment) %>%
  summarise(
    median_freq = median(time_diff, na.rm = TRUE),
    mean_positions = mean(n_positions, na.rm = TRUE),
    min_positions = min(n_positions, na.rm = TRUE),
    max_positions = max(n_positions, na.rm = TRUE),
    avg_days_mon = mean(n_days_mon, na.rm = TRUE),
    min_days = min(n_days_mon, na.rm = TRUE),
    max_days = max(n_days_mon, na.rm = TRUE),
    mean_eff_n = mean(effective_n, na.rm = TRUE),
    min_eff_n = min(effective_n, na.rm = TRUE),
    max_eff_n = max(effective_n, na.rm = TRUE)
  ) %>% 
  mutate_if(is.numeric, round, 0)

(summary_table <- 
  summary_table %>% 
  mutate(median_freq = median_freq,
         positions = paste(mean_positions, "(", min_positions, "-", max_positions, ")", sep = ""),
         days_mon = paste(avg_days_mon, "(", min_days, "-", max_days, ")", sep = ""),
         effective_n = paste(mean_eff_n, "(", min_eff_n, "-", max_eff_n, ")", sep = "")) %>% 
  select(Species, treatment, median_freq, positions, days_mon, effective_n))
  

(muddyfoot_species_positions_sum <- 
    flextable(summary_table) %>% 
    fontsize(part = "all", size = 11) %>% 
    bold(part = 'header')  %>% 
    set_header_labels("individual_id" = 'Fish ID',
                      "treatment" = 'Treatment',
                      "median_freq" = 'Frequency (s)',
                      "positions" = 'Locations',
                      "days_mon" = 'Duration (days)',
                      "effective_n" = 'N')) %>% 
  width(width = 3, unit = 'cm') %>% 
  width(j = c(1,2), width = 2, unit = 'cm')
  
save_as_docx(muddyfoot_species_positions_sum, 
                           path = paste0(save_tables_path, "muddyfoot_species_location_summary_MS.docx"))
             