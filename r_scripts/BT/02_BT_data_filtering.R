### LIBRARIES ###
# Load necessary libraries
library(data.table)   # For fast data manipulation
library(tidyverse)    # For data wrangling
library(ctmm)         # For continuous-time movement modeling
library(sf)           # For handling spatial data
library(parallel)     # For parallel processing
library(foreach)      # For parallel for loops
library(doParallel)   # For registering parallel backend

# Set the time zone to ensure consistent time handling
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths for reading and saving filtered telemetry and ctmm model results
filtered_data_path <- "./data/tracks_filtered/lake_BT/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"

#Load in the datasets
lake_BT_sub <- readRDS(paste0(filtered_data_path, '01_lake_BT_sub.rds'))

#-------------------------------------------------------------------------------#
# Early data exploration
#------------------------------------------------------------------------------#

#Create data subset
set.seed(123)  # reproducible
lake_BT_sub_dt <- as.data.table(lake_BT_sub)
selected_ids <- lake_BT_sub_dt[
  Species != "FReference",
  .(individual_ID = sample(unique(individual_ID), min(2, .N))),
  by = Species
]$individual_ID

selected_ids

ref_id <- lake_BT_sub_dt[individual_ID == "FReference", unique(individual_ID)]
selected_ids <- c(selected_ids, ref_id)

subset_ids <- lake_BT_sub_dt[individual_ID %in% selected_ids]

# subset to only include dates between
unique(lake_BT_sub_dt$date_cest)
keep_dates <- as.Date(c(
  "2022-09-28", "2022-09-29", "2022-09-30"
))

subset_ids <- subset_ids[date_cest %in% keep_dates]


subset_ids[, timestamp_cest := as.POSIXct(timestamp_cest, tz = "Europe/Stockholm")]
setorder(subset_ids, individual_ID, timestamp_cest)


# Prepare a dataframe for conversion to a telemetry object
# This format is compatible with the movebank telemetry format
subet_ids_movebank <- 
  with(subset_ids, 
       data.frame(
         "timestamp" = timestamp_cest,                        
         "location.long" = Long,                         
         "location.lat" = Lat, 
         "GPS.HDOP" = HPE,                               
         "individual-local-identifier" = individual_ID,
         "Species" = Species,
         "Weight" = Weight,
         "Total_length" = Total_length,
         "Std_length" = Std_length,
         "Treatment" = Treatment,                        
         "Date" = date_cest,
         "Exp_Stage" = Stage,
         "Time_Of_Day" = time_of_day
       ))


# Convert the dataframe to a telemetry object using the ctmm package's as.telemetry function
# No projection is applied here, and the WGS84 datum is used (common geographic coordinate system)
subset_id_tel_list <- as.telemetry(subet_ids_movebank, 
                                  timezone = "Europe/Stockholm",   # Set time zone
                                  timeformat = "%Y-%m-%d %H:%M:%S",# Specify timestamp format
                                  projection = NULL,               # No projection
                                  datum = "WGS84",                 # World Geodetic System 1984
                                  keep = c("Species", "Weight","Total_length", "Std_length", "Treatment", 
                                           "Date", "Exp_Stage", "Time_Of_Day"))  # Retain extra columns


plot(subset_id_tel_list[[1]], error = FALSE)
plot(subset_id_tel_list[[5]], error = FALSE)
plot(subset_id_tel_list[[7]], error = FALSE) #Reference tag

#create variograms
subset_varios <- lapply(subset_id_tel_list, function(tel) {
  v <- variogram(tel, fast = TRUE)
  plot(v)   # will pop up a plot per individual
  v
})

#plotting parameters for variograms
level <- c(0.95)
xlim <- c(0, 5 %#% "hour")

#Roach example semivariance
plot(subset_varios[[1]], xlim = xlim, level = level)
plot(subset_varios[[1]], level = level)
guess.roach = ctmm.guess(subset_id_tel_list[[1]], name = "roach.guess", interactive = TRUE)

#Perch example semivariance
plot(subset_varios[[3]], xlim = xlim, level = level)
plot(subset_varios[[3]], level = level)
guess.perch = ctmm.guess(subset_id_tel_list[[1]], name = "perch.guess", interactive = TRUE)

#Pike example semivariance
plot(subset_varios[[5]], xlim = xlim, level = level)
plot(subset_varios[[5]], level = level)
guess.pike = ctmm.guess(subset_id_tel_list[[6]], name = "pike.guess", interactive = TRUE)


#run some exploratory movement models
#first create seperate species lists
roach_subset_id_tel_list <- subset_id_tel_list[c(1,2)]
perch_subset_id_tel_list <- subset_id_tel_list[c(3,4)]
pike_subset_id_tel_list <- subset_id_tel_list[c(5,6)]


#roach sample ctmms
cl <- makeCluster(2)
doParallel::registerDoParallel(cl)
BT_sample_roach_ctmm_fits <-  
  foreach(i = 1:length(roach_subset_id_tel_list), .packages = 'ctmm') %dopar% {
    model_fit <- ctmm.select(roach_subset_id_tel_list[[i]], roach.guess, verbose = TRUE)
    model_fit
  }

stopCluster(cl)

summary(BT_sample_roach_ctmm_fits[[1]])
summary(BT_sample_roach_ctmm_fits[[2]])

#perch sample ctmms

BT_sample_perch_ctmm_fits <-  
  foreach(i = 1:length(perch_subset_id_tel_list), .packages = 'ctmm') %dopar% {
    model_fit <- ctmm.select(perch_subset_id_tel_list[[i]], perch.guess, verbose = TRUE)
    model_fit
  }





#convert to data table
setDT(lake_BT_sub)
  
# 1. Compute time between successive detections and per-individual summaries
lake_BT_sub <- lake_BT_sub %>% 
  arrange(individual_ID, timestamp_cest) %>%         # ensure correct temporal order within individuals
  group_by(individual_ID) %>% 
  mutate(
    time_diff        = c(NA, diff(timestamp_cest)),  # lagged time difference
    time_diff        = round(as.numeric(time_diff), 3),
    mean_time_diff   = mean(time_diff, na.rm = TRUE),
    median_time_diff = median(time_diff, na.rm = TRUE),
    n_positions      = n(),
    n_days_tracked   = n_distinct(date_cest)
  ) %>% 
  ungroup()

# 2. Create a compact BT tracking summary table
BT_track_sum <- lake_BT_sub %>% 
    group_by(individual_ID) %>% 
    dplyr::select(
      individual_ID, Treatment, Species,
      mean_time_diff, median_time_diff,
      n_positions, n_days_tracked
    ) %>% 
    mutate(
      mean_time_diff   = round(as.numeric(mean_time_diff), 2),
      median_time_diff = round(as.numeric(median_time_diff), 2)) %>% 
    distinct()
  
BT_track_sum_table <- 
    flextable(BT_track_sum) %>% 
    fontsize(part = "all", size = 11) %>% 
    bold(part = 'header') %>% 
    set_header_labels("individual_ID" = 'ID',
                      "mean_time_diff" = 'Mean lag',
                      "median_time_diff" = 'Median lag',
                      "n_positions" = 'Locations',
                      "n_days_tracked" = 'Days')
  
  
  # Optional: Save BT tracking summary.
  save_as_docx(BT_track_sum_table, path = paste0(save_filtered_data_path, "BT_track_sum_table_unfilt.docx"))
  
  #summary per species
  BT_track_sum %>% 
    group_by(Species) %>% 
    summarise(mean_lag = mean(mean_time_diff),
              median_lag = mean(median_time_diff),
              mean_positions = mean(n_positions))
  
  # Species       mean_lag median_lag mean_positions
  # <chr>            <dbl>      <dbl>          <dbl>
  # 1 Northern Pike    19.0        1.12        410999.
  # 2 Perch            23.0        2.74        150847.
  # 3 Roach             6.68       2.12        454310.
  
  
  # 3. Compute per-day sampling intensity and median time step for each individual
  lake_BT_sub <- lake_BT_sub %>% 
    group_by(individual_ID, date_cest) %>% 
    mutate(
      n_positions_day      = n(),
      median_day_time_diff = median(time_diff, na.rm = TRUE),
      n_positions_hourly   = n_positions_day / 24,
      n_positions_per_min  = n_positions_hourly / 60
    ) %>% 
    ungroup()
  
  
  
  # Create a compact per-ID / per-day tracking summary table
  # 1. Species-level average median_time_diff -------------------------------
  species_lag <- BT_track_sum %>% 
    group_by(Species) %>% 
    summarise(
      species_median_lag = mean(median_time_diff, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2. Per-ID, per-date summary with new columns ---------------------------
  ID_track_sum <- lake_BT_sub %>% 
    dplyr::select(
      individual_ID, Treatment, Species, date_cest,
      n_positions, n_days_tracked,
      median_day_time_diff, n_positions_hourly, n_positions_per_min
    ) %>% 
    distinct() %>% 
    left_join(species_lag, by = "Species") %>% 
    group_by(individual_ID, Species) %>% 
    mutate(
      # NEW COLUMN 1: Is this date ABOVE the species’ median lag?
      above_species_median = median_day_time_diff > species_median_lag,
      
      # NEW COLUMN 2: How many dates for this individual are above the species median?
      n_days_above_species_median = sum(above_species_median, na.rm = TRUE)
    ) %>% 
    ungroup()
  
  
  
  saveRDS(lake_BT_sub_sum, file = paste0(save_filtered_data_path, "lake_BT_sub_spatialfilt_sum.rds"))
  
  
  
  
  
  
  
  
  
  #Convert to dataframe
  
  lake_BT_sub <- as.data.frame(lake_BT_sub)
  
  
  
  