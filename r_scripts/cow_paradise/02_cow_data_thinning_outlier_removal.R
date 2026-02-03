#==============================================================================
# Lake Cow Paradise Detection Data - Outlier Removal and Temporal Thinning
#==============================================================================
# Purpose: Remove outliers, thin data, and filter mortality events
# Author: Marcus Michelangeli
#
# Input files:
#   - ./data/tracks_filtered/cow_paradise/01_cow_paradise_sub.rds
#   - ./data/lake_coords/lake_cow_polygon.gpkg
#   - ./data/pike_deaths.csv
#
# Output:
#   - ./data/telem_obj/cow_paradise/cow_paradise_UERE.rds
#   - ./data/telem_obj/cow_paradise/[species]_cow_paradise_tel_unthinned.rds (3 files)
#   - ./data/tracks_filtered/cow_paradise/02_cow_paradise_sub.rds
#   - ./data/tracks_filtered/cow_paradise/03_cow_paradise_sub.rds
#   - ./data/lake_coords/cow_paradise_boundary_coordinates.csv (NEW)
#   - ./data/lake_coords/cow_paradise_lake_data.rds (NEW)
#   - ./daily_trajectory_plots/cow_paradise/[individual]_daily_trajectory_plot.png
#==============================================================================

# Load required libraries ---------------------------------------------------
library(data.table)
library(tidyverse)
library(ctmm)
library(sf)
library(parallel)
library(foreach)
library(doParallel)
library(geosphere)
library(move2)

# Set timezone globally -----------------------------------------------------
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths ---------------------------------------------------------
filtered_data_path <- "./data/tracks_filtered/cow_paradise/"
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/cow_paradise/"
polygon_path <- "./data/lake_params/polygons/"
rec_loc_path <- "./data/lake_params/reciever_and_habitat_locations/"
plot_output_path <- "./daily_trajectory_plots/cow_paradise/"

# Import data ---------------------------------------------------------------
cow_paradise_sub <- readRDS(paste0(filtered_data_path, '01_cow_sub.rds'))
cow_paradise_sub_dt <- as.data.table(cow_paradise_sub)
cow_paradise_polygon <- st_read(paste0(polygon_path, "cow_polygon.gpkg"))

message("Imported ", nrow(cow_paradise_sub_dt), " detections for ", 
        n_distinct(cow_paradise_sub_dt$individual_ID), " individuals")

#==============================================================================
# 1. ESTIMATE LOCATION ERROR FROM REFERENCE TAG
#==============================================================================

# Extract reference tag detections ------------------------------------------
ref <- cow_paradise_sub_dt[individual_ID == "FReference"]
message("\n=== Reference Tag Analysis ===")
message("Reference tag detections: ", nrow(ref)) #15646

# Known fixed coordinates of reference tag ---------------------------------
# Reference receiver number
unique(ref$FullId)
#"H170-1802-65064"

# Get the coordinates for the location of the reference receiver
rec_locs <- fread(paste0(rec_loc_path, "cow_paradise_rec_hab_locations.csv"))

#The location of the reference tag
ref_tag_loc <- rec_locs %>% filter(ID == 'ref')

#The first location recorded of the reference tag
ref_first_location_coords <- rec_locs %>% filter(ID == 'ref2')

#The average coordinates recorded from the reference tag
ref_coords_avg <- ref %>%
  summarise(
    avg_lat = mean(Lat, na.rm = TRUE),
    avg_long = mean(Long, na.rm = TRUE),
    n_detections = n()
  )

# Use these average coordinates as the true reference position
true_lat <- ref_coords_avg$avg_lat
true_long <- ref_coords_avg$avg_long

# Function to parse DMS string and convert to decimal degrees
parse_dms_to_decimal <- function(dms_string) {
  # Remove direction letters (N, S, E, W) and clean the string
  dms_clean <- gsub("[NSEW]", "", dms_string)
  
  # Extract degrees, minutes, and seconds using regex
  # Pattern matches: number followed by degree symbol (or �), then minutes with ', then seconds with "
  pattern <- "([0-9]+)[°�]([0-9]+)'([0-9.]+)\""
  matches <- regmatches(dms_clean, regexec(pattern, dms_clean))[[1]]
  
  if(length(matches) == 4) {
    degrees <- as.numeric(matches[2])
    minutes <- as.numeric(matches[3])
    seconds <- as.numeric(matches[4])
    
    decimal <- degrees + minutes/60 + seconds/3600
    
    # Apply negative sign if South or West
    if(grepl("[SW]", dms_string)) {
      decimal <- -decimal
    }
    
    return(decimal)
  } else {
    stop("Could not parse DMS string: ", dms_string)
  }
}

# Extract coordinates from ref_coords dataframe
true_lat_ref_loc <- parse_dms_to_decimal(ref_tag_loc$Latitude)
true_long_ref_loc <- parse_dms_to_decimal(ref_tag_loc$Longitude)

message("Reference tag coordinates:")
message("  Latitude: ", round(true_lat, 7), "°N (", ref_coords$Latitude, ")")
message("  Longitude: ", round(true_long, 7), "°E (", ref_coords$Longitude, ")")

# Calculate distance between each detection and true position (m) -----------
ref[, error_dist_m := distHaversine(
  cbind(Long, Lat),
  cbind(true_long, true_lat)
)]

ref[, error_dist_m_ref_loc := distHaversine(
  cbind(Long, Lat),
  cbind(true_long_ref_loc, true_lat_ref_loc)
)]

# Summarize positional error ------------------------------------------------
message("\n=== Reference Tag Location Error (meters) ===")
print(summary(ref$error_dist_m))
print(summary(ref$error_dist_m_ref_loc))

sigma_error_med <- median(ref$error_dist_m, na.rm = TRUE)
message("Median error: ", round(sigma_error_med, 4), " m")
#Median error: 0.9551 m based on average coordinates

sigma_error_med_ref_loc <- median(ref$error_dist_m_ref_loc, na.rm = TRUE)
message("Median error: ", round(sigma_error_med_ref_loc, 4), " m")
#Median error: 1.0569 m based on the marked location of the reference tag

#In summary, the error is about 1 m

#==============================================================================
# 2. ANALYZE MOVEMENT SCALE VS TIME LAG
#==============================================================================
# Purpose: Determine optimal temporal thinning interval where movement
# is reliably distinguishable from location error

# Select representative individuals (2 per species) -------------------------
set.seed(123)  # For reproducibility
selected_ids <- cow_paradise_sub_dt[
  Species != "FReference",
  .(individual_ID = sample(unique(individual_ID), min(2, .N))),
  by = Species
]$individual_ID

subset_ids <- cow_paradise_sub_dt[individual_ID %in% selected_ids]
setorder(subset_ids, individual_ID, timestamp_cest)

# Calculate step metrics ----------------------------------------------------
# Time difference between successive fixes (seconds)
subset_ids[, dt_sec := as.numeric(
  difftime(timestamp_cest, data.table::shift(timestamp_cest), units = "secs")
), by = individual_ID]

# Step length between successive fixes (m) using Haversine distance
subset_ids[, step_dist_m := distHaversine(
  cbind(data.table::shift(Long), data.table::shift(Lat)),
  cbind(Long, Lat)
), by = individual_ID]

# Bin step lengths by time lag ----------------------------------------------
# Define time lag bins (seconds)
breaks <- c(0, 3, 7, 12, 18, 25, 40, 70, Inf)
labels <- c("~2", "~5", "~10", "~15", "~20", "~30", ">30", ">70")

subset_ids[!is.na(step_dist_m),
           dt_bin := cut(dt_sec,
                         breaks = breaks,
                         labels = labels,
                         include.lowest = TRUE)]

# Summarize movement by time lag bin ----------------------------------------
step_summary <- subset_ids[!is.na(step_dist_m),
                           .(
                             n_steps = .N,
                             median_step = median(step_dist_m),
                             mean_step = mean(step_dist_m),
                             q25_step = quantile(step_dist_m, 0.25),
                             q75_step = quantile(step_dist_m, 0.75)
                           ),
                           by = dt_bin]

step_summary[, dt_bin := factor(dt_bin, levels = labels)]
setorder(step_summary, dt_bin)

# Calculate movement-to-error ratio -----------------------------------------
step_summary[, movement_to_error := median_step / lake_cow_UERE$UERE[,1]]

message("\n=== Movement vs Time Lag Analysis ===")
print(step_summary)

message("\nConclusion: Movement becomes most distinguishable from error at 30s+")

#==============================================================================
# 3. PREPARE DATA FOR OUTLIER DETECTION
#==============================================================================

# Convert to Movebank format for ctmm package ------------------------------
cow_paradise_movebank <- with(
  cow_paradise_sub,
  data.frame(
    "timestamp" = timestamp_cest,
    "location.long" = Long,
    "location.lat" = Lat,
    "GPS.HDOP" = HPE,
    "individual-local-identifier" = individual_ID,
    "species" = Species,
    "weight" = Weight,
    "total_length" = Total_length,
    "std_length" = Std_length,
    "treatment" = Treatment,
    "date" = date_cest,
    "exp_stage" = Stage,
    "time_of_day" = time_of_day,
    "found_alive" = found_alive,
    "known_predated" = known_predated
  )
)

# Convert to telemetry object -----------------------------------------------
cow_paradise_tels <- as.telemetry(
  cow_paradise_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,
  datum = "WGS84",
  keep = c("species", "weight", "total_length", "std_length", "treatment",
           "date", "exp_stage", "time_of_day", "found_alive", "known_predated")
)

# Center projection on geometric median -------------------------------------
projection(cow_paradise_tels) <- ctmm::median(cow_paradise_tels)

#==============================================================================
# 4. INCORPORATE LOCATION ERROR MODEL
#==============================================================================

# Fit error parameters using reference tag ----------------------------------
cow_paradise_UERE <- uere.fit(cow_paradise_tels$FReference)
print(summary(cow_paradise_UERE))
# , , horizontal
# 
# low       est      high
# all 0.8200953 0.8265717 0.8330474

# Apply error model to all telemetry objects -------------------------------
uere(cow_paradise_tels) <- cow_paradise_UERE

# Remove reference tag from analysis ---------------------------------------
# Note: Adjust the index range based on actual number of individuals
n_individuals <- length(cow_paradise_tels) - 1
cow_paradise_tels <- cow_paradise_tels[1:n_individuals]
message("Removed reference tag. Remaining individuals: ", length(cow_paradise_tels))

# Save UERE model -----------------------------------------------------------
saveRDS(cow_paradise_UERE, paste0(telem_path, "cow_paradise_UERE.rds"))

#==============================================================================
# 5. REMOVE SPEED-BASED OUTLIERS BY SPECIES
#==============================================================================
# Maximum sustained swimming speeds (Ucrit in m/s) from literature
# Reference: "Key factors explaining critical swimming speed in freshwater fish:
# A review and statistical analysis using Iberian species"

speed_thresholds <- list(
  Pike = 0.823,
  Perch = 0.977,
  Roach = 0.841
)

# Organize telemetry objects by species ------------------------------------
# First, verify species distribution and sort by individual ID
species_order <- cow_paradise_movebank %>%
  select(species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

message("\n=== Species Distribution (Before Filtering) ===")
print(table(species_order$species, useNA = "ifany"))

# Filter out reference tag and any individuals with NA species
species_order_filtered <- species_order %>%
  filter(!is.na(species))

message("\n=== Species Distribution (After Filtering) ===")
print(table(species_order_filtered$species))

# Add index column to match telemetry list order
species_order_filtered$telem_index <- 1:nrow(species_order_filtered)

# Display the mapping table for verification
message("\n=== Individual ID to Species Mapping (Filtered) ===")
print(species_order_filtered)

# Automatically extract indices for each species
pike_indices <- species_order_filtered$telem_index[species_order_filtered$species == "Northern Pike"]
perch_indices <- species_order_filtered$telem_index[species_order_filtered$species == "Perch"]
roach_indices <- species_order_filtered$telem_index[species_order_filtered$species == "Roach"]

# Extract telemetry objects by species
pike_cow_paradise_tel <- cow_paradise_tels[pike_indices]
perch_cow_paradise_tel <- cow_paradise_tels[perch_indices]
roach_cow_paradise_tel <- cow_paradise_tels[roach_indices]

message("\n=== Species Telemetry Objects Created ===")
message("Pike: ", length(pike_cow_paradise_tel), " individuals (indices: ", 
        paste(pike_indices, collapse = ", "), ")")
message("Perch: ", length(perch_cow_paradise_tel), " individuals (indices: ", 
        paste(perch_indices, collapse = ", "), ")")
message("Roach: ", length(roach_cow_paradise_tel), " individuals (indices: ", 
        paste(roach_indices, collapse = ", "), ")")

# Verify totals match
total_individuals <- length(pike_cow_paradise_tel) + length(perch_cow_paradise_tel) + length(roach_cow_paradise_tel)
message("\nTotal individuals across all species: ", total_individuals)
message("Expected total (excluding reference tag): ", length(cow_paradise_tels))

if(total_individuals != length(cow_paradise_tels)) {
  warning("Mismatch in total individuals! Check species assignment.")
}

# Function to remove speed outliers -----------------------------------------
remove_speed_outliers <- function(telem_list, species_name, max_speed) {
  
  message("\n--- Processing ", species_name, " ---")
  message("Max speed threshold: ", max_speed, " m/s")
  
  # Calculate speeds
  out <- outlie(telem_list, plot = FALSE)
  
  # Count outliers
  n_outliers <- sum(sapply(out, function(x) sum(x$speed > max_speed)))
  message("Outliers detected: ", n_outliers)
  
  # Filter to retain only realistic speeds
  which_lowSp <- lapply(out, function(x) x$speed <= max_speed)
  filtered_telem <- Map(function(x, y) x[y, ], telem_list, which_lowSp)
  
  message("Filtering complete")
  return(filtered_telem)
}

# Apply outlier removal to each species ------------------------------------
pike_cow_paradise_tel <- remove_speed_outliers(
  pike_cow_paradise_tel, "Pike", speed_thresholds$Pike
) #outliers detected: 9882
saveRDS(pike_cow_paradise_tel, paste0(telem_path, "pike_cow_paradise_tel_unthinned.rds"))

perch_cow_paradise_tel <- remove_speed_outliers(
  perch_cow_paradise_tel, "Perch", speed_thresholds$Perch
) #outliers detected: 33878
saveRDS(perch_cow_paradise_tel, paste0(telem_path, "perch_cow_paradise_tel_unthinned.rds"))

roach_cow_paradise_tel <- remove_speed_outliers(
  roach_cow_paradise_tel, "Roach", speed_thresholds$Roach
) #outliers detected: 339565
saveRDS(roach_cow_paradise_tel, paste0(telem_path, "roach_cow_paradise_tel_unthinned.rds"))

#==============================================================================
# 6. COMBINE FILTERED DATA FROM ALL SPECIES
#==============================================================================

# Function to extract and combine telemetry data ---------------------------
extract_data <- function(list_data) {
  data_combined <- do.call(rbind, lapply(names(list_data), function(id) {
    df <- list_data[[id]]
    df$individual_ID <- id
    return(df)
  }))
  return(data_combined)
}

# Extract and label data by species ----------------------------------------
pike_data <- extract_data(pike_cow_paradise_tel)
perch_data <- extract_data(perch_cow_paradise_tel)
roach_data <- extract_data(roach_cow_paradise_tel)

# Combine into single dataframe --------------------------------------------
cow_paradise_telem_data <- rbind(
  cbind(pike_data, Species = 'Pike'),
  cbind(perch_data, Species = 'Perch'),
  cbind(roach_data, Species = 'Roach')
)

message("\n=== Data After Outlier Removal ===")
message("Total detections: ", nrow(cow_paradise_telem_data)) #Total detections: 8783485
message("Outliers removed: ", nrow(cow_paradise_sub_dt) - nrow(ref) - nrow(cow_paradise_telem_data)) #Outliers removed: 213038

# Save outlier-filtered data ------------------------------------------------
saveRDS(cow_paradise_telem_data, paste0(filtered_data_path, "02_cow_paradise_sub.rds"))

# Clean up large objects ----------------------------------------------------
rm(cow_paradise_movebank, cow_paradise_sub, cow_paradise_sub_dt, cow_paradise_tels,
   roach_data, perch_data, pike_data,
   roach_cow_paradise_tel, pike_cow_paradise_tel, perch_cow_paradise_tel,
   subset_ids, ref)
gc()

#==============================================================================
# 7. TEMPORAL THINNING
#==============================================================================

cow_paradise_telem_data <- readRDS(paste0(filtered_data_path, "02_cow_paradise_sub.rds"))

# Convert to move2 object ---------------------------------------------------
cow_paradise_mv <- mt_as_move2(
  cow_paradise_telem_data,
  coords = c("longitude", "latitude"),
  crs = "WGS84",
  time_column = "timestamp",
  track_id_column = "individual_ID",
  na.fail = FALSE
)

cow_paradise_mv <- cow_paradise_mv %>%
  arrange(individual_ID, timestamp)

# Apply temporal thinning ---------------------------------------------------

#15 seconds
# thinning_interval <- "15 seconds"
# message("\n=== Temporal Thinning ===")
# message("Applying interval: ", thinning_interval)
# 
# cow_paradise_thin_15_data <- cow_paradise_mv %>%
#   mt_filter_per_interval(unit = thinning_interval, criterion = "first")
# 
# message("Before thinning: ", nrow(cow_paradise_mv), " detections") #10459343 detections
# message("After thinning: ", nrow(cow_paradise_thin_data), " detections") #4952567 detections
# message("Reduction: ", round(100 * (1 - nrow(cow_paradise_thin_data) / nrow(cow_paradise_mv)), 1), "%") #Reduction: 52.6%

#30 seconds
thinning_interval <- "30 seconds"
message("\n=== Temporal Thinning ===")
message("Applying interval: ", thinning_interval)

cow_paradise_thin_data <- cow_paradise_mv %>%
  mt_filter_per_interval(unit = thinning_interval, criterion = "first")

message("Before thinning: ", nrow(cow_paradise_mv), " detections") #8783485  detections
message("After thinning: ", nrow(cow_paradise_thin_data), " detections") #After thinning: 1963882  detections
message("Reduction: ", round(100 * (1 - nrow(cow_paradise_thin_data) / nrow(cow_paradise_mv)), 1), "%") #Reduction: 77.6%%

# Verify thinning intervals -------------------------------------------------
cow_paradise_thin_track_sum <- cow_paradise_thin_data %>%
  group_by(individual_ID) %>%
  arrange(timestamp) %>%
  mutate(dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
  summarise(
    min_dt = min(dt, na.rm = TRUE),
    median_dt = median(dt, na.rm = TRUE),
    max_dt = max(dt, na.rm = TRUE)
  )

message("\n=== Temporal Resolution Summary ===")
print(summary(cow_paradise_thin_track_sum[, c("min_dt", "median_dt", "max_dt")]))
# min_dt         median_dt         max_dt                 geometry 
# Min.   :0.7997   Min.   :30.26   Min.   :  51396   MULTIPOINT   :63  
# 1st Qu.:1.8001   1st Qu.:30.66   1st Qu.:  62953   epsg:4326    : 0  
# Median :1.8007   Median :31.48   Median :  80700   +proj=long...: 0  
# Mean   :1.7222   Mean   :32.27   Mean   : 138793                     
# 3rd Qu.:1.8009   3rd Qu.:33.49   3rd Qu.: 159461                     
# Max.   :1.8254   Max.   :39.15   Max.   :1174976    

#==============================================================================
# 8. IDENTIFY TEMPORAL GAPS IN TRACKING
#==============================================================================

# Convert to data.table and calculate time gaps ----------------------------
cow_paradise_thin_data <- as.data.table(cow_paradise_thin_data)
setorder(cow_paradise_thin_data, individual_ID, timestamp)

# Calculate time gaps between successive fixes (minutes)
cow_paradise_thin_data[, dt_mins := as.numeric(
  difftime(timestamp, data.table::shift(timestamp), units = "mins")
), by = individual_ID]

# Store previous timestamp and date (start of gap)
cow_paradise_thin_data[, `:=`(
  prev_time = data.table::shift(timestamp),
  prev_date = data.table::shift(date)
), by = individual_ID]

# Find maximum gap for each individual --------------------------------------
max_gaps <- cow_paradise_thin_data[
  !is.na(dt_mins),
  .SD[which.max(dt_mins)],
  by = individual_ID
][, .(
  individual_ID,
  gap_start_time = prev_time,
  gap_start_date = prev_date,
  gap_end_time = timestamp,
  gap_length_mins = dt_mins
)]

message("\n=== Largest Temporal Gaps by Individual ===")
print(max_gaps[order(-gap_length_mins)])

# Save temporal gaps to Excel -----------------------------------------------
library(openxlsx)

# Create formatted Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "Temporal Gaps Summary")

# Prepare data for export with formatted columns
gaps_export <- max_gaps[order(-gap_length_mins)]
gaps_export[, `:=`(
  gap_start_time = format(gap_start_time, "%Y-%m-%d %H:%M:%S"),
  gap_end_time = format(gap_end_time, "%Y-%m-%d %H:%M:%S"),
  gap_length_hours = round(gap_length_mins / 60, 2),
  gap_length_days = round(gap_length_mins / 1440, 2)
)]

# Reorder and rename columns for clarity
gaps_export <- gaps_export[, .(
  Individual_ID = individual_ID,
  Gap_Start_Date = gap_start_date,
  Gap_Start_Time = gap_start_time,
  Gap_End_Time = gap_end_time,
  Gap_Length_Minutes = round(gap_length_mins, 1),
  Gap_Length_Hours = gap_length_hours,
  Gap_Length_Days = gap_length_days
)]

# Write data to worksheet
writeData(wb, "Temporal Gaps Summary", gaps_export, startRow = 1)

# Format header row
headerStyle <- createStyle(
  fontSize = 12,
  fontColour = "#FFFFFF",
  halign = "center",
  fgFill = "#4F81BD",
  border = "TopBottomLeftRight",
  borderColour = "#4F81BD",
  textDecoration = "bold"
)

addStyle(wb, "Temporal Gaps Summary", headerStyle, rows = 1, cols = 1:7, gridExpand = TRUE)

# Format data cells
dataStyle <- createStyle(
  halign = "left",
  border = "TopBottomLeftRight",
  borderColour = "#CCCCCC"
)

addStyle(wb, "Temporal Gaps Summary", dataStyle, 
         rows = 2:(nrow(gaps_export) + 1), 
         cols = 1:7, 
         gridExpand = TRUE)

# Format numeric columns
numStyle <- createStyle(
  halign = "right",
  border = "TopBottomLeftRight",
  borderColour = "#CCCCCC"
)

addStyle(wb, "Temporal Gaps Summary", numStyle,
         rows = 2:(nrow(gaps_export) + 1),
         cols = 5:7,
         gridExpand = TRUE)

# Highlight large gaps (> 24 hours)
largeGapStyle <- createStyle(
  fgFill = "#FFC7CE",
  fontColour = "#9C0006"
)

large_gap_rows <- which(gaps_export$Gap_Length_Hours > 24) + 1
if(length(large_gap_rows) > 0) {
  addStyle(wb, "Temporal Gaps Summary", largeGapStyle,
           rows = large_gap_rows,
           cols = 1:7,
           gridExpand = TRUE)
}

# Set column widths
setColWidths(wb, "Temporal Gaps Summary", cols = 1:7, 
             widths = c(15, 15, 20, 20, 18, 16, 15))

# Freeze header row
freezePane(wb, "Temporal Gaps Summary", firstRow = TRUE)

# Save the workbook
output_file <- paste0(filtered_data_path, "cow_paradise_temporal_gaps_summary.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

message("\nTemporal gaps summary saved to Excel: ", output_file)
message("Rows with gaps > 24 hours are highlighted in red")

# Extract coordinates from geometry -----------------------------------------
coords_data <- st_coordinates(cow_paradise_thin_data$geometry)
cow_paradise_thin_data$Long <- coords_data[, 1]
cow_paradise_thin_data$Lat <- coords_data[, 2]

#Remove unneeded columns
cow_paradise_thin_data <- cow_paradise_thin_data %>% 
  dplyr::select(-Species, -dt_mins, -prev_time, -prev_date)

#==============================================================================
# 9. VISUALIZE DAILY MOVEMENT TRAJECTORIES
#==============================================================================

# Function to plot daily trajectories for one individual -------------------
plot_daily_traj_individual <- function(dat, lake_poly) {
  
  this_id <- unique(dat$individual_ID)
  
  ggplot() +
    geom_sf(data = lake_poly, fill = "lightblue", color = "black", alpha = 0.3) +
    geom_path(
      data = dat,
      aes(x = Long, y = Lat, group = date),
      linewidth = 0.3
    ) +
    geom_point(
      data = dat,
      aes(x = Long, y = Lat),
      size = 0.4
    ) +
    facet_wrap(~ date) +
    coord_sf() +
    theme_minimal() +
    labs(
      title = paste("Daily Movement Trajectories -", this_id),
      x = "Longitude",
      y = "Latitude"
    )
}

# Generate plots for all individuals ----------------------------------------
message("\n=== Generating Daily Trajectory Plots ===")

traj_plots <- cow_paradise_thin_data %>%
  group_by(individual_ID) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$individual_ID))) %>%
  map(~ plot_daily_traj_individual(.x, cow_paradise_polygon))

# Save plots to files -------------------------------------------------------
walk(names(traj_plots), function(id) {
  filename <- paste0(plot_output_path, id, "_daily_trajectory_plot.png")
  
  ggsave(
    filename = filename,
    plot = traj_plots[[id]],
    width = 10,
    height = 8,
    dpi = 300
  )
  
  message("Saved: ", filename)
})


#Identified some weird tracking going on particularly during later dates.
#I am going to calculate some positional information 
#We may need to remove some dates where tracking was really poor
# Count positions per date
positions_per_date <- cow_paradise_thin_data %>%
  group_by(date) %>%
  summarise(n_positions = n(), .groups = 'drop')

# View the counts
print(positions_per_date, n = 35)
#the number of posistion falls by about 50% from the 19-10-2022

#remove date 31-10-2022
cow_paradise_thin_data <- cow_paradise_thin_data %>%
  filter(date != as.Date("2022-10-31"))


#The tracking from the 19th seems odd for most individuals
# Filter data from 2022-10-19 onwards
filtered_data <- cow_paradise_thin_data %>%
  filter(date >= as.Date("2022-10-19"))

# Calculate positions per individual per day and compare to daily median
individual_daily_positions <- filtered_data %>%
  group_by(date, individual_ID) %>%
  summarise(n_positions = n(), .groups = 'drop') %>%
  group_by(date) %>%
  mutate(
    median_positions = median(n_positions),
    below_median = n_positions < median_positions
  ) %>%
  ungroup()

# View individuals below median
below_median_individuals <- individual_daily_positions %>%
  filter(below_median == TRUE)

print(below_median_individuals)

# Optional: Summary of how many individuals are below median each day
summary_below_median <- below_median_individuals %>%
  group_by(date) %>%
  summarise(
    n_individuals_below_median = n(),
    .groups = 'drop'
  )

print(summary_below_median)
# Plot showing all individuals with those below median highlighted
ggplot(individual_daily_positions, aes(x = date, y = n_positions, color = below_median)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(y = median_positions), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                     labels = c("At/Above Median", "Below Median")) +
  labs(title = "Daily Positions per Individual",
       subtitle = "Red points indicate individuals below median for that day",
       x = "Date", 
       y = "Number of Positions",
       color = "") +
  theme_minimal()



#Individuals identified as having particularly poor tracking
# F59851
# F59853



# Calculate positions per individual per day and compare to daily median
species_daily_positions <- cow_paradise_thin_track_sum %>%
  group_by(date, species) %>%
  summarise(n_positions = n(), .groups = 'drop') %>%
  group_by(date) %>%
  mutate(
    median_positions = median(n_positions),
  ) %>%
  ungroup()



# Save final filtered dataset -----------------------------------------------
saveRDS(cow_paradise_thin_data, paste0(filtered_data_path, "03_cow_sub.rds"))

message("\n=== Data Processing Complete ===")
message("Final dataset: ", nrow(cow_paradise_thin_data), " detections") #Final dataset: 1963300 detections
message("Number of individuals: ", n_distinct(cow_paradise_thin_data$individual_ID)) #63 individuals
message("Date range: ", min(cow_paradise_thin_data$date), " to ", 
        max(cow_paradise_thin_data$date)) #Date range: 2022-09-27 to 2022-10-30
message("Saved to: ", paste0(filtered_data_path, "03_cow_paradise_sub.rds"))

#==============================================================================
# SUMMARY OF OUTPUT FILES
#==============================================================================
message("\n=== Output Files Generated ===")
message("1. Lake boundary coordinates: ", paste0(lake_polygon_path, "cow_paradise_boundary_coordinates.csv"))
message("2. Lake data summary: ", paste0(lake_polygon_path, "cow_paradise_lake_data.rds"))
message("3. UERE model: ", paste0(save_telem_path, "cow_paradise/cow_paradise_UERE.rds"))
message("4. Outlier-filtered data: ", paste0(filtered_data_path, "02_cow_paradise_sub.rds"))
message("5. Final thinned data: ", paste0(filtered_data_path, "03_cow_paradise_sub.rds"))
message("6. Daily trajectory plots: ", plot_output_path)

















#==============================================================================
# Lake cow_paradise Detection Data - Outlier Removal and Temporal Thinning
#==============================================================================
# Purpose: Remove outliers, thin data, and filter mortality events
# Author: Marcus Michelangeli
#
# Input files:
#   - ./data/tracks_filtered/lake_cow_paradise/01_lake_cow_paradise_sub.rds
#   - ./data/lake_coords/lake_cow_paradise_polygon.gpkg
#   - ./data/pike_deaths.csv
#
# Output:
#   - ./data/telem_obj/cow_paradise/lake_cow_paradise_UERE.rds
#   - ./data/telem_obj/cow_paradise/[species]_lake_cow_paradise_tel_unthinned.rds (3 files)
#   - ./data/tracks_filtered/lake_cow_paradise/02_lake_cow_paradise_sub.rds
#   - ./data/tracks_filtered/lake_cow_paradise/03_lake_cow_paradise_sub.rds
#   - ./daily_trajectory_plots/lake_cow_paradise/[individual]_daily_trajectory_plot.png
#==============================================================================

# Load required libraries ---------------------------------------------------
library(data.table)
library(tidyverse)
library(ctmm)
library(sf)
library(parallel)
library(foreach)
library(doParallel)
library(geosphere)
library(move2)

# Set timezone globally -----------------------------------------------------
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths ---------------------------------------------------------
filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"
lake_polygon_path <- "./data/lake_coords/"
plot_output_path <- "./daily_trajectory_plots/lake_cow_paradise/"

# Import data ---------------------------------------------------------------
lake_cow_paradise_sub <- readRDS(paste0(filtered_data_path, '01_lake_cow_paradise_sub.rds'))
lake_cow_paradise_sub_dt <- as.data.table(lake_cow_paradise_sub)
cow_paradise_polygon <- st_read(paste0(lake_polygon_path, "lake_cow_paradise_polygon.gpkg"))

message("Imported ", nrow(lake_cow_paradise_sub_dt), " detections for ", 
        n_distinct(lake_cow_paradise_sub_dt$individual_ID), " individuals")

#==============================================================================
# 1. ESTIMATE LOCATION ERROR FROM REFERENCE TAG
#==============================================================================

# Extract reference tag detections ------------------------------------------
ref <- lake_cow_paradise_sub_dt[individual_ID == "FReference"]
message("Reference tag detections: ", nrow(ref)) #88576

# Known fixed coordinates of reference tag ---------------------------------
true_long <- 20.0472989
true_lat <- 63.7707882

# Calculate distance between each detection and true position (m) -----------
ref[, error_dist_m := distHaversine(
  cbind(Long, Lat),
  cbind(true_long, true_lat)
)]

# Summarize positional error ------------------------------------------------
message("\n=== Reference Tag Location Error (meters) ===")
print(summary(ref$error_dist_m))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003683 0.0810723 0.1092805 0.1117890 0.1363092 3.8638857


sigma_error_med <- median(ref$error_dist_m, na.rm = TRUE)
message("Median error: ", round(sigma_error_med, 4), " m")
# Median error: 0.1093 m


#==============================================================================
# 2. ANALYZE MOVEMENT SCALE VS TIME LAG
#==============================================================================
# Purpose: Determine optimal temporal thinning interval where movement
# is reliably distinguishable from location error

# Select representative individuals (2 per species) -------------------------
set.seed(123)  # For reproducibility
selected_ids <- lake_cow_paradise_sub_dt[
  Species != "FReference",
  .(individual_ID = sample(unique(individual_ID), min(2, .N))),
  by = Species
]$individual_ID

subset_ids <- lake_cow_paradise_sub_dt[individual_ID %in% selected_ids]
setorder(subset_ids, individual_ID, timestamp_cest)

# Calculate step metrics ----------------------------------------------------
# Time difference between successive fixes (seconds)
subset_ids[, dt_sec := as.numeric(
  difftime(timestamp_cest, shift(timestamp_cest), units = "secs")
), by = individual_ID]

# Step length between successive fixes (m) using Haversine distance
subset_ids[, step_dist_m := distHaversine(
  cbind(shift(Long), shift(Lat)),
  cbind(Long, Lat)
), by = individual_ID]

# Bin step lengths by time lag ----------------------------------------------
# Define time lag bins (seconds)
breaks <- c(0, 3, 7, 12, 18, 25, 40, 70, Inf)
labels <- c("~2", "~5", "~10", "~15", "~20", "~30", ">30", ">70")

subset_ids[!is.na(step_dist_m),
           dt_bin := cut(dt_sec,
                         breaks = breaks,
                         labels = labels,
                         include.lowest = TRUE)]

# Summarize movement by time lag bin ----------------------------------------
step_summary <- subset_ids[!is.na(step_dist_m),
                           .(
                             n_steps = .N,
                             median_step = median(step_dist_m),
                             mean_step = mean(step_dist_m),
                             q25_step = quantile(step_dist_m, 0.25),
                             q75_step = quantile(step_dist_m, 0.75)
                           ),
                           by = dt_bin]

step_summary[, dt_bin := factor(dt_bin, levels = labels)]
setorder(step_summary, dt_bin)

# Calculate movement-to-error ratio -----------------------------------------
step_summary[, movement_to_error := median_step / sigma_error_med]

message("\n=== Movement vs Time Lag Analysis ===")
print(step_summary)

message("\nConclusion: Movement becomes reliably distinguishable from error at 10s+")
message("Reliability stabilizes at 15-20s. No meaningful improvement beyond 20-30s")
message("Selected thinning interval: 10 seconds")

#==============================================================================
# 3. PREPARE DATA FOR OUTLIER DETECTION
#==============================================================================

# Convert to Movebank format for ctmm package ------------------------------
lake_cow_paradise_movebank <- with(
  lake_cow_paradise_sub,
  data.frame(
    "timestamp" = timestamp_cest,
    "location.long" = Long,
    "location.lat" = Lat,
    "GPS.HDOP" = HPE,
    "individual-local-identifier" = individual_ID,
    "species" = Species,
    "weight" = Weight,
    "total_length" = Total_length,
    "std_length" = Std_length,
    "treatment" = Treatment,
    "date" = date_cest,
    "exp_stage" = Stage,
    "time_of_day" = time_of_day,
    "found_alive" = found_alive,
    "known_predated" = Known_predated
  )
)

# Convert to telemetry object -----------------------------------------------
lake_cow_paradise_tels <- as.telemetry(
  lake_cow_paradise_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,
  datum = "WGS84",
  keep = c("Species", "Weight", "Total_length", "Std_length", "Treatment",
           "Date", "Exp_Stage", "Time_Of_Day", "found_alive", "known_predated")
)

# Center projection on geometric median -------------------------------------
projection(lake_cow_paradise_tels) <- ctmm::median(lake_cow_paradise_tels)

#==============================================================================
# 4. INCORPORATE LOCATION ERROR MODEL
#==============================================================================

# Fit error parameters using reference tag ----------------------------------
lake_cow_paradise_UERE <- uere.fit(lake_cow_paradise_tels$FReference)
message("\n=== UERE Model Summary ===")
print(summary(lake_cow_paradise_UERE))
# low      est      high
# all 0.37805 0.379299 0.3805479

# Apply error model to all telemetry objects -------------------------------
uere(lake_cow_paradise_tels) <- lake_cow_paradise_UERE

# Remove reference tag from analysis ---------------------------------------
lake_cow_paradise_tels <- lake_cow_paradise_tels[1:66]
message("Removed reference tag. Remaining individuals: ", length(lake_cow_paradise_tels))

# Save UERE model -----------------------------------------------------------
saveRDS(lake_cow_paradise_UERE, paste0(save_telem_path, "cow_paradise/lake_cow_paradise_UERE.rds"))

#==============================================================================
# 5. REMOVE SPEED-BASED OUTLIERS BY SPECIES
#==============================================================================
# Maximum sustained swimming speeds (Ucrit in m/s) from literature
# Reference: "Key factors explaining critical swimming speed in freshwater fish:
# A review and statistical analysis using Iberian species"

speed_thresholds <- list(
  Pike = 0.823,
  Perch = 0.977,
  Roach = 0.841
)

# Organize telemetry objects by species ------------------------------------
# First, verify and sort by individual ID

# Assign individuals to species groups
pike_lake_cow_paradise_tel <- lake_cow_paradise_tels[61:66]
perch_lake_cow_paradise_tel <- lake_cow_paradise_tels[c(1:15, 31:44, 46)]
roach_lake_cow_paradise_tel <- lake_cow_paradise_tels[c(16:30, 45, 47:60)]

message("\nSpecies distribution:")
message("Pike: ", length(pike_lake_cow_paradise_tel), " individuals")
message("Perch: ", length(perch_lake_cow_paradise_tel), " individuals")
message("Roach: ", length(roach_lake_cow_paradise_tel), " individuals")

# Function to remove speed outliers -----------------------------------------
remove_speed_outliers <- function(telem_list, species_name, max_speed) {
  
  message("\n--- Processing ", species_name, " ---")
  message("Max speed threshold: ", max_speed, " m/s")
  
  # Calculate speeds
  out <- outlie(telem_list, plot = FALSE)
  
  # Count outliers
  n_outliers <- sum(sapply(out, function(x) sum(x$speed > max_speed)))
  message("Outliers detected: ", n_outliers)
  
  # Filter to retain only realistic speeds
  which_lowSp <- lapply(out, function(x) x$speed <= max_speed)
  filtered_telem <- Map(function(x, y) x[y, ], telem_list, which_lowSp)
  
  message("Filtering complete")
  return(filtered_telem)
}

# Apply outlier removal to each species ------------------------------------
pike_lake_cow_paradise_tel <- remove_speed_outliers(
  pike_lake_cow_paradise_tel, "Pike", speed_thresholds$Pike
) #outliers removed: 11538
saveRDS(pike_lake_cow_paradise_tel, paste0(save_telem_path, "cow_paradise/pike_lake_cow_paradise_tel_unthinned.rds"))

perch_lake_cow_paradise_tel <- remove_speed_outliers(
  perch_lake_cow_paradise_tel, "Perch", speed_thresholds$Perch
) #outliers removed: 33878
saveRDS(perch_lake_cow_paradise_tel, paste0(save_telem_path, "cow_paradise/perch_lake_cow_paradise_tel_unthinned.rds"))

roach_lake_cow_paradise_tel <- remove_speed_outliers(
  roach_lake_cow_paradise_tel, "Roach", speed_thresholds$Roach
) #outlier removed: 114260
saveRDS(roach_lake_cow_paradise_tel, paste0(save_telem_path, "cow_paradise/roach_lake_cow_paradise_tel_unthinned.rds"))

#==============================================================================
# 6. COMBINE FILTERED DATA FROM ALL SPECIES
#==============================================================================

# Function to extract and combine telemetry data ---------------------------
extract_data <- function(list_data) {
  data_combined <- do.call(rbind, lapply(names(list_data), function(id) {
    df <- list_data[[id]]
    df$individual_ID <- id
    return(df)
  }))
  return(data_combined)
}

# Extract and label data by species ----------------------------------------
pike_data <- extract_data(pike_lake_cow_paradise_tel)
perch_data <- extract_data(perch_lake_cow_paradise_tel)
roach_data <- extract_data(roach_lake_cow_paradise_tel)

# Combine into single dataframe --------------------------------------------
cow_paradise_telem_data <- rbind(
  cbind(pike_data, Species = 'Pike'),
  cbind(perch_data, Species = 'Perch'),
  cbind(roach_data, Species = 'Roach')
)

message("\n=== Data After Outlier Removal ===")
message("Total detections: ", nrow(cow_paradise_telem_data)) #20019313
message("Outliers removed: ", nrow(lake_cow_paradise_sub_dt) - nrow(ref) - nrow(cow_paradise_telem_data)) #159676

# Save outlier-filtered data ------------------------------------------------
saveRDS(cow_paradise_telem_data, paste0(filtered_data_path, "02_lake_cow_paradise_sub.rds"))

# Clean up large objects ----------------------------------------------------
rm(lake_cow_paradise_movebank, lake_cow_paradise_sub, lake_cow_paradise_sub_dt, lake_cow_paradise_tels,
   roach_data, perch_data, pike_data,
   roach_lake_cow_paradise_tel, pike_lake_cow_paradise_tel, perch_lake_cow_paradise_tel,
   subset_ids, ref)
gc()

#==============================================================================
# 7. TEMPORAL THINNING
#==============================================================================

# Convert to move2 object ---------------------------------------------------
lake_cow_paradise_mv <- mt_as_move2(
  cow_paradise_telem_data,
  coords = c("longitude", "latitude"),
  crs = "WGS84",
  time_column = "timestamp",
  track_id_column = "individual_ID",
  na.fail = FALSE
)

lake_cow_paradise_mv <- lake_cow_paradise_mv %>%
  arrange(individual_ID, timestamp)

# Apply temporal thinning ---------------------------------------------------
thinning_interval <- "10 seconds"
message("\nApplying temporal thinning: ", thinning_interval)

lake_cow_paradise_thin_data <- lake_cow_paradise_mv %>%
  mt_filter_per_interval(unit = thinning_interval, criterion = "first")

message("Before thinning: ", nrow(lake_cow_paradise_mv), " detections") #20019313 detections
message("After thinning: ", nrow(lake_cow_paradise_thin_data), " detections") #7339440 detections
message("Reduction: ", round(100 * (1 - nrow(lake_cow_paradise_thin_data) / nrow(lake_cow_paradise_mv)), 1), "%") #Reduction: 63.3%

# Verify thinning intervals -------------------------------------------------
lake_cow_paradise_thin_track_sum <- lake_cow_paradise_thin_data %>%
  group_by(individual_ID) %>%
  arrange(timestamp) %>%
  mutate(dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
  summarise(
    min_dt = min(dt, na.rm = TRUE),
    median_dt = median(dt, na.rm = TRUE),
    max_dt = max(dt, na.rm = TRUE)
  )

message("\n=== Temporal Resolution Summary ===")
print(summary(lake_cow_paradise_thin_track_sum[, c("min_dt", "median_dt", "max_dt")]))

# min_dt         median_dt         max_dt                 geometry 
# Min.   :0.7993   Min.   :10.05   Min.   :   1249   MULTIPOINT   :66  
# 1st Qu.:1.7992   1st Qu.:10.18   1st Qu.:  29077   epsg:4326    : 0  
# Median :1.7999   Median :10.28   Median :  47226   +proj=long...: 0  
# Mean   :1.7089   Mean   :10.49   Mean   : 101684                     
# 3rd Qu.:1.8004   3rd Qu.:10.47   3rd Qu.:  86707                     
# Max.   :1.8009   Max.   :12.35   Max.   :1300537      

#==============================================================================
# 8. IDENTIFY TEMPORAL GAPS IN TRACKING
#==============================================================================

# Convert to data.table and calculate time gaps ----------------------------
lake_cow_paradise_thin_data <- as.data.table(lake_cow_paradise_thin_data)
setorder(lake_cow_paradise_thin_data, individual_ID, timestamp)

# Calculate time gaps between successive fixes (minutes)
lake_cow_paradise_thin_data[, dt_mins := as.numeric(
  difftime(timestamp, shift(timestamp), units = "mins")
), by = individual_ID]

# Store previous timestamp and date (start of gap)
lake_cow_paradise_thin_data[, `:=`(
  prev_time = shift(timestamp),
  prev_date = shift(Date)
), by = individual_ID]

# Find maximum gap for each individual --------------------------------------
max_gaps <- lake_cow_paradise_thin_data[
  !is.na(dt_mins),
  .SD[which.max(dt_mins)],
  by = individual_ID
][, .(
  individual_ID,
  gap_start_time = prev_time,
  gap_start_date = prev_date,
  gap_end_time = timestamp,
  gap_length_mins = dt_mins
)]

message("\n=== Largest Temporal Gaps by Individual ===")
print(max_gaps[order(-gap_length_mins)])

# Save temporal gaps to Excel -----------------------------------------------
library(openxlsx)

# Create formatted Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "Temporal Gaps Summary")

# Prepare data for export with formatted columns
gaps_export <- max_gaps[order(-gap_length_mins)]
gaps_export[, `:=`(
  gap_start_time = format(gap_start_time, "%Y-%m-%d %H:%M:%S"),
  gap_end_time = format(gap_end_time, "%Y-%m-%d %H:%M:%S"),
  gap_length_hours = round(gap_length_mins / 60, 2),
  gap_length_days = round(gap_length_mins / 1440, 2)
)]

# Reorder and rename columns for clarity
gaps_export <- gaps_export[, .(
  Individual_ID = individual_ID,
  Gap_Start_Date = gap_start_date,
  Gap_Start_Time = gap_start_time,
  Gap_End_Time = gap_end_time,
  Gap_Length_Minutes = round(gap_length_mins, 1),
  Gap_Length_Hours = gap_length_hours,
  Gap_Length_Days = gap_length_days
)]

# Write data to worksheet
writeData(wb, "Temporal Gaps Summary", gaps_export, startRow = 1)

# Format header row
headerStyle <- createStyle(
  fontSize = 12,
  fontColour = "#FFFFFF",
  halign = "center",
  fgFill = "#4F81BD",
  border = "TopBottomLeftRight",
  borderColour = "#4F81BD",
  textDecoration = "bold"
)

addStyle(wb, "Temporal Gaps Summary", headerStyle, rows = 1, cols = 1:7, gridExpand = TRUE)

# Format data cells
dataStyle <- createStyle(
  halign = "left",
  border = "TopBottomLeftRight",
  borderColour = "#CCCCCC"
)

addStyle(wb, "Temporal Gaps Summary", dataStyle, 
         rows = 2:(nrow(gaps_export) + 1), 
         cols = 1:7, 
         gridExpand = TRUE)

# Format numeric columns
numStyle <- createStyle(
  halign = "right",
  border = "TopBottomLeftRight",
  borderColour = "#CCCCCC"
)

addStyle(wb, "Temporal Gaps Summary", numStyle,
         rows = 2:(nrow(gaps_export) + 1),
         cols = 5:7,
         gridExpand = TRUE)

# Highlight large gaps (> 24 hours)
largeGapStyle <- createStyle(
  fgFill = "#FFC7CE",
  fontColour = "#9C0006"
)

large_gap_rows <- which(gaps_export$Gap_Length_Hours > 24) + 1
if(length(large_gap_rows) > 0) {
  addStyle(wb, "Temporal Gaps Summary", largeGapStyle,
           rows = large_gap_rows,
           cols = 1:7,
           gridExpand = TRUE)
}

# Set column widths
setColWidths(wb, "Temporal Gaps Summary", cols = 1:7, 
             widths = c(15, 15, 20, 20, 18, 16, 15))

# Freeze header row
freezePane(wb, "Temporal Gaps Summary", firstRow = TRUE)

# Save the workbook
output_file <- paste0(filtered_data_path, "lake_cow_paradise_temporal_gaps_summary.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

message("\nTemporal gaps summary saved to Excel: ", output_file)
message("Rows with gaps > 24 hours are highlighted in red")

# Extract coordinates from geometry -----------------------------------------
coords_data <- st_coordinates(lake_cow_paradise_thin_data$geometry)
lake_cow_paradise_thin_data$Long <- coords_data[, 1]
lake_cow_paradise_thin_data$Lat <- coords_data[, 2]

#==============================================================================
# 9. VISUALIZE DAILY MOVEMENT TRAJECTORIES
#==============================================================================

# Function to plot daily trajectories for one individual -------------------
plot_daily_traj_individual <- function(dat, lake_poly) {
  
  this_id <- unique(dat$individual_ID)
  
  ggplot() +
    geom_sf(data = lake_poly, fill = "lightblue", color = "black", alpha = 0.3) +
    geom_path(
      data = dat,
      aes(x = Long, y = Lat, group = Date),
      linewidth = 0.3
    ) +
    geom_point(
      data = dat,
      aes(x = Long, y = Lat),
      size = 0.4
    ) +
    facet_wrap(~ Date) +
    coord_sf() +
    theme_minimal() +
    labs(
      title = paste("Daily Movement Trajectories -", this_id),
      x = "Longitude",
      y = "Latitude"
    )
}

# Generate plots for all individuals ----------------------------------------
message("\n=== Generating Daily Trajectory Plots ===")

traj_plots <- lake_cow_paradise_thin_data %>%
  group_by(individual_ID) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$individual_ID))) %>%
  map(~ plot_daily_traj_individual(.x, cow_paradise_polygon))

# Save plots to files -------------------------------------------------------
walk(names(traj_plots), function(id) {
  filename <- paste0(plot_output_path, id, "_daily_trajectory_plot.png")
  
  ggsave(
    filename = filename,
    plot = traj_plots[[id]],
    width = 10,
    height = 8,
    dpi = 300
  )
  
  message("Saved: ", filename)
})

#==============================================================================
# 10. FILTER MORTALITY EVENTS
#==============================================================================

# Remove individual that died on first day ----------------------------------
n_before <- n_distinct(lake_cow_paradise_thin_data$individual_ID)
lake_cow_paradise_thin_data <- lake_cow_paradise_thin_data %>%
  filter(individual_ID != "F59889")

message("\nRemoved F59889 (died on Day 1)")
message("Individuals remaining: ", n_distinct(lake_cow_paradise_thin_data$individual_ID))

# Load pike mortality information -------------------------------------------
pike_deaths <- read.csv("./data/pike_deaths.csv")

pike_mort_cols <- pike_deaths %>%
  filter(individual_ID %in% lake_cow_paradise_thin_data$individual_ID) %>%
  mutate(pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y"))

message("\nPike mortality records loaded: ", nrow(pike_mort_cols))

# Filter out detections after mortality events ------------------------------
n_rows_before <- nrow(lake_cow_paradise_thin_data)

lake_cow_paradise_thin_data <- lake_cow_paradise_thin_data %>%
  left_join(pike_mort_cols, by = "individual_ID") %>%
  filter(is.na(pike_death_date) | Date < pike_death_date)

n_rows_after <- nrow(lake_cow_paradise_thin_data)
message("Detections removed after mortality events: ", n_rows_before - n_rows_after) #48037

# Save final filtered dataset -----------------------------------------------
saveRDS(lake_cow_paradise_thin_data, paste0(filtered_data_path, "03_lake_cow_paradise_sub.rds"))

message("Final dataset: ", nrow(lake_cow_paradise_thin_data), " detections") #7143462 detections
message("Number of individuals: ", n_distinct(lake_cow_paradise_thin_data$individual_ID)) #Number of individuals: 65
message("Date range: ", min(lake_cow_paradise_thin_data$Date), " to ", 
        max(lake_cow_paradise_thin_data$Date))
message("Saved to: ", paste0(filtered_data_path, "03_lake_cow_paradise_sub.rds"))





#==============================================================================
# Lake Cow Paradise Detection Data - Outlier Removal and Temporal Thinning
#==============================================================================
# Purpose: Remove outliers, thin data, and filter mortality events
# Author: Marcus Michelangeli
# Input files:
#   - ./data/tracks_filtered/lake_cow_paradise/01_lake_cow_sub.rds
#   - ./data/lake_coords/lake_cow_polygon.gpkg
#   - ./data/pike_deaths.csv
#
# Output:
#   - ./data/telem_obj/cow_paradise/lake_cow_UERE.rds
#   - ./data/telem_obj/cow_paradise/[species]_lake_cow_tel_unthinned.rds (3 files)
#   - ./data/tracks_filtered/lake_cow_paradise/02_lake_cow_sub.rds
#   - ./data/tracks_filtered/lake_cow_paradise/03_lake_cow_sub.rds
#   - ./daily_trajectory_plots/lake_cow_paradise/[individual]_daily_trajectory_plot.png
#==============================================================================

# Load required libraries ---------------------------------------------------
library(data.table)
library(tidyverse)
library(ctmm)
library(sf)
library(parallel)
library(foreach)
library(doParallel)
library(geosphere)
library(move2)

# Set timezone globally -----------------------------------------------------
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths ---------------------------------------------------------
filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"
lake_polygon_path <- "./data/lake_coords/"
plot_output_path <- "./daily_trajectory_plots/lake_cow_paradise/"

# Import data ---------------------------------------------------------------
lake_cow_sub <- readRDS(paste0(filtered_data_path, '01_lake_cow_sub.rds'))
lake_cow_sub_dt <- as.data.table(lake_cow_sub)
cow_polygon <- st_read(paste0(lake_polygon_path, "lake_cow_polygon.gpkg"))

message("Imported ", nrow(lake_cow_sub_dt), " detections for ", 
        n_distinct(lake_cow_sub_dt$individual_ID), " individuals")

#==============================================================================
# 1. ESTIMATE LOCATION ERROR FROM REFERENCE TAG
#==============================================================================

# Extract reference tag detections ------------------------------------------
ref <- lake_cow_sub_dt[individual_ID == "FReference"]
message("Reference tag detections: ", nrow(ref)) #15646

# Known fixed coordinates of reference tag ---------------------------------
# Update these coordinates to match Lake Cow Paradise reference tag location
true_long <- 20.06347865833  
true_lat <- 63.77199213611   

# Calculate distance between each detection and true position (m) -----------
ref[, error_dist_m := distHaversine(
  cbind(Long, Lat),
  cbind(true_long, true_lat)
)]

# Summarize positional error ------------------------------------------------
message("\n=== Reference Tag Location Error (meters) ===")
print(summary(ref$error_dist_m))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3908  0.8034  1.0569  1.4183  2.1318 13.6876 

sigma_error_med <- median(ref$error_dist_m, na.rm = TRUE)
message("Median error: ", round(sigma_error_med, 4), " m")
#Median error: 1.0569 m


#==============================================================================
# 2. ANALYZE MOVEMENT SCALE VS TIME LAG
#==============================================================================
# Purpose: Determine optimal temporal thinning interval where movement
# is reliably distinguishable from location error

# Select representative individuals (2 per species) -------------------------
set.seed(123)  # For reproducibility
selected_ids <- lake_cow_sub_dt[
  Species != "FReference",
  .(individual_ID = sample(unique(individual_ID), min(2, .N))),
  by = Species
]$individual_ID

subset_ids <- lake_cow_sub_dt[individual_ID %in% selected_ids]
setorder(subset_ids, individual_ID, timestamp_cest)

# Calculate step metrics ----------------------------------------------------
# Time difference between successive fixes (seconds)
subset_ids[, dt_sec := as.numeric(
  difftime(timestamp_cest, shift(timestamp_cest), units = "secs")
), by = individual_ID]

# Step length between successive fixes (m) using Haversine distance
subset_ids[, step_dist_m := distHaversine(
  cbind(shift(Long), shift(Lat)),
  cbind(Long, Lat)
), by = individual_ID]

# Bin step lengths by time lag ----------------------------------------------
# Define time lag bins (seconds)
breaks <- c(0, 3, 7, 12, 18, 25, 40, 70, Inf)
labels <- c("~2", "~5", "~10", "~15", "~20", "~30", ">30", ">70")

subset_ids[!is.na(step_dist_m),
           dt_bin := cut(dt_sec,
                         breaks = breaks,
                         labels = labels,
                         include.lowest = TRUE)]

# Summarize movement by time lag bin ----------------------------------------
step_summary <- subset_ids[!is.na(step_dist_m),
                           .(
                             n_steps = .N,
                             median_step = median(step_dist_m),
                             mean_step = mean(step_dist_m),
                             q25_step = quantile(step_dist_m, 0.25),
                             q75_step = quantile(step_dist_m, 0.75)
                           ),
                           by = dt_bin]

step_summary[, dt_bin := factor(dt_bin, levels = labels)]
setorder(step_summary, dt_bin)

# Calculate movement-to-error ratio -----------------------------------------
step_summary[, movement_to_error := median_step / sigma_error_med]

message("\n=== Movement vs Time Lag Analysis ===")
print(step_summary)


#==============================================================================
# 3. PREPARE DATA FOR OUTLIER DETECTION
#==============================================================================

# Convert to Movebank format for ctmm package ------------------------------
lake_cow_movebank <- with(
  lake_cow_sub,
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
    "Time_Of_Day" = time_of_day,
    "found_alive" = found_alive,
    "known_predated" = Known_predated
  )
)

# Convert to telemetry object -----------------------------------------------
lake_cow_tels <- as.telemetry(
  lake_cow_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,
  datum = "WGS84",
  keep = c("Species", "Weight", "Total_length", "Std_length", "Treatment",
           "Date", "Exp_Stage", "Time_Of_Day", "found_alive", "known_predated")
)

# Center projection on geometric median -------------------------------------
projection(lake_cow_tels) <- ctmm::median(lake_cow_tels)

#==============================================================================
# 4. INCORPORATE LOCATION ERROR MODEL
#==============================================================================

# Fit error parameters using reference tag ----------------------------------
lake_cow_UERE <- uere.fit(lake_cow_tels$FReference)
message("\n=== UERE Model Summary ===")
print(summary(lake_cow_UERE))
# low       est      high
# all 0.8200953 0.8265717 0.8330474

# Apply error model to all telemetry objects -------------------------------
uere(lake_cow_tels) <- lake_cow_UERE

# Remove reference tag from analysis ---------------------------------------
# Adjust the index based on actual number of individuals
n_individuals <- length(lake_cow_tels) - 1
lake_cow_tels <- lake_cow_tels[1:n_individuals]
message("Removed reference tag. Remaining individuals: ", length(lake_cow_tels))

# Save UERE model -----------------------------------------------------------
saveRDS(lake_cow_UERE, paste0(save_telem_path, "cow_paradise/lake_cow_UERE.rds"))

#==============================================================================
# 5. REMOVE SPEED-BASED OUTLIERS BY SPECIES
#==============================================================================
# Maximum sustained swimming speeds (Ucrit in m/s) from literature
# Reference: "Key factors explaining critical swimming speed in freshwater fish:
# A review and statistical analysis using Iberian species"

speed_thresholds <- list(
  Pike = 0.823,
  Perch = 0.977,
  Roach = 0.841
)

# Organize telemetry objects by species ------------------------------------
# First, verify species distribution and sort by individual ID
species_order <- lake_cow_movebank %>%
  select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

message("\n=== Species Distribution (Before Filtering) ===")
print(table(species_order$Species, useNA = "ifany"))

# Filter out reference tag and any individuals with NA species
species_order_filtered <- species_order %>%
  filter(!is.na(Species))

message("\n=== Species Distribution (After Filtering) ===")
print(table(species_order_filtered$Species))

# Add index column to match telemetry list order
species_order_filtered$telem_index <- 1:nrow(species_order_filtered)

# Display the mapping table for verification
message("\n=== Individual ID to Species Mapping (Filtered) ===")
print(species_order_filtered)

# Automatically extract indices for each species
pike_indices <- species_order_filtered$telem_index[species_order_filtered$Species == "Northern Pike"]
perch_indices <- species_order_filtered$telem_index[species_order_filtered$Species == "Perch"]
roach_indices <- species_order_filtered$telem_index[species_order_filtered$Species == "Roach"]

# Extract telemetry objects by species
pike_lake_cow_tel <- lake_cow_tels[pike_indices]
perch_lake_cow_tel <- lake_cow_tels[perch_indices]
roach_lake_cow_tel <- lake_cow_tels[roach_indices]

message("\nSpecies distribution:")
message("Pike: ", length(pike_lake_cow_tel), " individuals")
message("Perch: ", length(perch_lake_cow_tel), " individuals")
message("Roach: ", length(roach_lake_cow_tel), " individuals")

# Function to remove speed outliers -----------------------------------------
remove_speed_outliers <- function(telem_list, species_name, max_speed) {
  
  message("\n--- Processing ", species_name, " ---")
  message("Max speed threshold: ", max_speed, " m/s")
  
  # Calculate speeds
  out <- outlie(telem_list, plot = FALSE)
  
  # Count outliers
  n_outliers <- sum(sapply(out, function(x) sum(x$speed > max_speed)))
  message("Outliers detected: ", n_outliers)
  
  # Filter to retain only realistic speeds
  which_lowSp <- lapply(out, function(x) x$speed <= max_speed)
  filtered_telem <- Map(function(x, y) x[y, ], telem_list, which_lowSp)
  
  message("Filtering complete")
  return(filtered_telem)
}

# Apply outlier removal to each species ------------------------------------
pike_lake_cow_tel <- remove_speed_outliers(
  pike_lake_cow_tel, "Pike", speed_thresholds$Pike
) #Outliers detected: 26812
saveRDS(pike_lake_cow_tel, paste0(save_telem_path, "cow_paradise/pike_lake_cow_tel_unthinned.rds"))

perch_lake_cow_tel <- remove_speed_outliers(
  perch_lake_cow_tel, "Perch", speed_thresholds$Perch
) #Outliers detected: 58011
saveRDS(perch_lake_cow_tel, paste0(save_telem_path, "cow_paradise/perch_lake_cow_tel_unthinned.rds"))

roach_lake_cow_tel <- remove_speed_outliers(
  roach_lake_cow_tel, "Roach", speed_thresholds$Roach
) #Outliers detected: 128508
saveRDS(roach_lake_cow_tel, paste0(save_telem_path, "cow_paradise/roach_lake_cow_tel_unthinned.rds"))

#==============================================================================
# 6. COMBINE FILTERED DATA FROM ALL SPECIES
#==============================================================================

# Function to extract and combine telemetry data ---------------------------
extract_data <- function(list_data) {
  data_combined <- do.call(rbind, lapply(names(list_data), function(id) {
    df <- list_data[[id]]
    df$individual_ID <- id
    return(df)
  }))
  return(data_combined)
}

# Extract and label data by species ----------------------------------------
pike_data <- extract_data(pike_lake_cow_tel)
perch_data <- extract_data(perch_lake_cow_tel)
roach_data <- extract_data(roach_lake_cow_tel)

# Combine into single dataframe --------------------------------------------
cow_telem_data <- rbind(
  cbind(pike_data, Species = 'Pike'),
  cbind(perch_data, Species = 'Perch'),
  cbind(roach_data, Species = 'Roach')
)

message("\n=== Data After Outlier Removal ===")
message("Total detections: ", nrow(cow_telem_data)) #8860993
message("Outliers removed: ", nrow(lake_cow_sub_dt) - nrow(ref) - nrow(cow_telem_data)) #213331

# Save outlier-filtered data ------------------------------------------------
saveRDS(cow_telem_data, paste0(filtered_data_path, "02_lake_cow_sub.rds"))

# Clean up large objects ----------------------------------------------------
rm(lake_cow_movebank, lake_cow_sub, lake_cow_sub_dt, lake_cow_tels,
   roach_data, perch_data, pike_data,
   roach_lake_cow_tel, pike_lake_cow_tel, perch_lake_cow_tel,
   subset_ids, ref)
gc()

#==============================================================================
# 7. TEMPORAL THINNING
#==============================================================================

# Convert to move2 object ---------------------------------------------------
lake_cow_mv <- mt_as_move2(
  cow_telem_data,
  coords = c("longitude", "latitude"),
  crs = "WGS84",
  time_column = "timestamp",
  track_id_column = "individual_ID",
  na.fail = FALSE
)

lake_cow_mv <- lake_cow_mv %>%
  arrange(individual_ID, timestamp)

# Apply temporal thinning ---------------------------------------------------
thinning_interval <- "10 seconds"
message("\nApplying temporal thinning: ", thinning_interval)

lake_cow_thin_data <- lake_cow_mv %>%
  mt_filter_per_interval(unit = thinning_interval, criterion = "first")

message("Before thinning: ", nrow(lake_cow_mv), " detections") #8860993 detections
message("After thinning: ", nrow(lake_cow_thin_data), " detections") #4073867 detections
message("Reduction: ", round(100 * (1 - nrow(lake_cow_thin_data) / nrow(lake_cow_mv)), 1), "%") #Reduction: 54%

# Verify thinning intervals -------------------------------------------------
lake_cow_thin_track_sum <- lake_cow_thin_data %>%
  group_by(individual_ID) %>%
  arrange(timestamp) %>%
  mutate(dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
  summarise(
    min_dt = min(dt, na.rm = TRUE),
    median_dt = median(dt, na.rm = TRUE),
    max_dt = max(dt, na.rm = TRUE)
  )

message("\n=== Temporal Resolution Summary ===")
print(summary(lake_cow_thin_track_sum[, c("min_dt", "median_dt", "max_dt")]))

# min_dt         median_dt         max_dt                 geometry 
# Min.   :0.7992   Min.   :10.30   Min.   :  21801   MULTIPOINT   :66  
# 1st Qu.:1.7947   1st Qu.:11.27   1st Qu.:  61072   epsg:4326    : 0  
# Median :1.7985   Median :11.67   Median :  80453   +proj=long...: 0  
# Mean   :1.7072   Mean   :12.06   Mean   : 146064                     
# 3rd Qu.:1.7996   3rd Qu.:12.62   3rd Qu.: 160698                     
# Max.   :1.8008   Max.   :16.61   Max.   :1174976        


#==============================================================================
# 8. IDENTIFY TEMPORAL GAPS IN TRACKING
#==============================================================================

# Convert to data.table and calculate time gaps ----------------------------
lake_cow_thin_data <- as.data.table(lake_cow_thin_data)
setorder(lake_cow_thin_data, individual_ID, timestamp)

# Calculate time gaps between successive fixes (minutes)
lake_cow_thin_data[, dt_mins := as.numeric(
  difftime(timestamp, shift(timestamp), units = "mins")
), by = individual_ID]

# Store previous timestamp and date (start of gap)
lake_cow_thin_data[, `:=`(
  prev_time = shift(timestamp),
  prev_date = shift(Date)
), by = individual_ID]

# Find maximum gap for each individual --------------------------------------
max_gaps <- lake_cow_thin_data[
  !is.na(dt_mins),
  .SD[which.max(dt_mins)],
  by = individual_ID
][, .(
  individual_ID,
  gap_start_time = prev_time,
  gap_start_date = prev_date,
  gap_end_time = timestamp,
  gap_length_mins = dt_mins
)]

message("\n=== Largest Temporal Gaps by Individual ===")
print(max_gaps[order(-gap_length_mins)])

# Save temporal gaps to Excel -----------------------------------------------
library(openxlsx)

# Create formatted Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "Temporal Gaps Summary")

# Prepare data for export with formatted columns
gaps_export <- max_gaps[order(-gap_length_mins)]
gaps_export[, `:=`(
  gap_start_time = format(gap_start_time, "%Y-%m-%d %H:%M:%S"),
  gap_end_time = format(gap_end_time, "%Y-%m-%d %H:%M:%S"),
  gap_length_hours = round(gap_length_mins / 60, 2),
  gap_length_days = round(gap_length_mins / 1440, 2)
)]

# Reorder and rename columns for clarity
gaps_export <- gaps_export[, .(
  Individual_ID = individual_ID,
  Gap_Start_Date = gap_start_date,
  Gap_Start_Time = gap_start_time,
  Gap_End_Time = gap_end_time,
  Gap_Length_Minutes = round(gap_length_mins, 1),
  Gap_Length_Hours = gap_length_hours,
  Gap_Length_Days = gap_length_days
)]

# Write data to worksheet
writeData(wb, "Temporal Gaps Summary", gaps_export, startRow = 1)

# Format header row
headerStyle <- createStyle(
  fontSize = 12,
  fontColour = "#FFFFFF",
  halign = "center",
  fgFill = "#4F81BD",
  border = "TopBottomLeftRight",
  borderColour = "#4F81BD",
  textDecoration = "bold"
)

addStyle(wb, "Temporal Gaps Summary", headerStyle, rows = 1, cols = 1:7, gridExpand = TRUE)

# Format data cells
dataStyle <- createStyle(
  halign = "left",
  border = "TopBottomLeftRight",
  borderColour = "#CCCCCC"
)

addStyle(wb, "Temporal Gaps Summary", dataStyle, 
         rows = 2:(nrow(gaps_export) + 1), 
         cols = 1:7, 
         gridExpand = TRUE)

# Format numeric columns
numStyle <- createStyle(
  halign = "right",
  border = "TopBottomLeftRight",
  borderColour = "#CCCCCC"
)

addStyle(wb, "Temporal Gaps Summary", numStyle,
         rows = 2:(nrow(gaps_export) + 1),
         cols = 5:7,
         gridExpand = TRUE)

# Highlight large gaps (> 24 hours)
largeGapStyle <- createStyle(
  fgFill = "#FFC7CE",
  fontColour = "#9C0006"
)

large_gap_rows <- which(gaps_export$Gap_Length_Hours > 24) + 1
if(length(large_gap_rows) > 0) {
  addStyle(wb, "Temporal Gaps Summary", largeGapStyle,
           rows = large_gap_rows,
           cols = 1:7,
           gridExpand = TRUE)
}

# Set column widths
setColWidths(wb, "Temporal Gaps Summary", cols = 1:7, 
             widths = c(15, 15, 20, 20, 18, 16, 15))

# Freeze header row
freezePane(wb, "Temporal Gaps Summary", firstRow = TRUE)

# Save the workbook
output_file <- paste0(filtered_data_path, "lake_cow_temporal_gaps_summary.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

message("\nTemporal gaps summary saved to Excel: ", output_file)
message("Rows with gaps > 24 hours are highlighted in red")

# Extract coordinates from geometry -----------------------------------------
coords_data <- st_coordinates(lake_cow_thin_data$geometry)
lake_cow_thin_data$Long <- coords_data[, 1]
lake_cow_thin_data$Lat <- coords_data[, 2]

#==============================================================================
# 9. VISUALIZE DAILY MOVEMENT TRAJECTORIES
#==============================================================================

# Function to plot daily trajectories for one individual -------------------
plot_daily_traj_individual <- function(dat, lake_poly) {
  
  this_id <- unique(dat$individual_ID)
  
  ggplot() +
    geom_sf(data = lake_poly, fill = "lightblue", color = "black", alpha = 0.3) +
    geom_path(
      data = dat,
      aes(x = Long, y = Lat, group = Date),
      linewidth = 0.3
    ) +
    geom_point(
      data = dat,
      aes(x = Long, y = Lat),
      size = 0.4
    ) +
    facet_wrap(~ Date) +
    coord_sf() +
    theme_minimal() +
    labs(
      title = paste("Daily Movement Trajectories -", this_id),
      x = "Longitude",
      y = "Latitude"
    )
}

# Generate plots for all individuals ----------------------------------------
message("\n=== Generating Daily Trajectory Plots ===")

traj_plots <- lake_cow_thin_data %>%
  group_by(individual_ID) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$individual_ID))) %>%
  map(~ plot_daily_traj_individual(.x, cow_polygon))

# Save plots to files -------------------------------------------------------
walk(names(traj_plots), function(id) {
  filename <- paste0(plot_output_path, id, "_daily_trajectory_plot.png")
  
  ggsave(
    filename = filename,
    plot = traj_plots[[id]],
    width = 10,
    height = 8,
    dpi = 300
  )
  
  message("Saved: ", filename)
})

#individuals with very little tracking identified from the plots (need further exploration)
#F59819
#F59873
#F59893


#==============================================================================
# 10. FILTER MORTALITY EVENTS
#==============================================================================

# Note: Update individual ID if there are known early mortalities
# Remove individual that died on first day (if applicable) or had very few fixes -----------------
early_mortality_ids <- c("F59893", "F59873", "F59819")

lake_cow_thin_data <- lake_cow_thin_data %>%
  filter(!individual_ID %in% early_mortality_ids) 

message("Removed ", length(early_mortality_ids), " individuals with early mortality")
message("Individuals remaining: ", n_distinct(lake_cow_thin_data$individual_ID))

# Load pike mortality information -------------------------------------------
pike_deaths <- read.csv("./data/pike_deaths.csv")

pike_mort_cols <- pike_deaths %>%
  filter(individual_ID %in% lake_cow_thin_data$individual_ID) %>%
  mutate(pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y"))

message("\nPike mortality records loaded: ", nrow(pike_mort_cols))

# Filter out detections after mortality events ------------------------------
n_rows_before <- nrow(lake_cow_thin_data)

lake_cow_thin_data <- lake_cow_thin_data %>%
  left_join(pike_mort_cols, by = "individual_ID") %>%
  filter(is.na(pike_death_date) | Date < pike_death_date)

n_rows_after <- nrow(lake_cow_thin_data)
message("Detections removed after mortality events: ", n_rows_before - n_rows_after) #12657

# Save final filtered dataset -----------------------------------------------
saveRDS(lake_cow_thin_data, paste0(filtered_data_path, "03_lake_cow_sub.rds"))

message("\n=== Data Processing Complete ===")
message("Final dataset: ", nrow(lake_cow_thin_data), " detections") #4018236 detections
message("Number of individuals: ", n_distinct(lake_cow_thin_data$individual_ID)) #63 individuals
message("Date range: ", min(lake_cow_thin_data$Date), " to ", 
        max(lake_cow_thin_data$Date)) #2022-09-27 to 2022-10-31
message("Saved to: ", paste0(filtered_data_path, "03_lake_cow_sub.rds"))