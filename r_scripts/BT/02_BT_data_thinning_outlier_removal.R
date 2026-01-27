#==============================================================================
# Lake BT Detection Data - Outlier Removal and Temporal Thinning
#==============================================================================
# Purpose: Remove outliers, thin data, and filter mortality events
# Author: Marcus Michelangeli
#
# Input files:
#   - ./data/tracks_filtered/BT/01_BT_sub.rds
#   - ./data/lake_coords/lake_BT_polygon.gpkg
#   - ./data/pike_deaths.csv
#
# Output:
#   - ./data/telem_obj/BT/BT_UERE.rds
#   - ./data/telem_obj/BT/[species]_BT_tel_unthinned.rds (3 files)
#   - ./data/tracks_filtered/BT/02_BT_sub.rds
#   - ./data/tracks_filtered/BT/03_BT_sub.rds
#   - ./data/lake_coords/BT_boundary_coordinates.csv (NEW)
#   - ./data/lake_coords/BT_lake_data.rds (NEW)
#   - ./daily_trajectory_plots/BT/[individual]_daily_trajectory_plot.png
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
filtered_data_path <- "./data/tracks_filtered/BT/"
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/BT/"
polygon_path <- "./data/lake_params/polygons/"
rec_loc_path <- "./data/lake_params/reciever_and_habitat_locations/"
plot_output_path <- "./daily_trajectory_plots/BT/"

# Import data ---------------------------------------------------------------
BT_sub <- readRDS(paste0(filtered_data_path, '01_BT_sub.rds'))
BT_sub_dt <- as.data.table(BT_sub)
BT_polygon <- st_read(paste0(polygon_path, "BT_polygon.gpkg"))

message("Imported ", nrow(BT_sub_dt), " detections for ", 
        n_distinct(BT_sub_dt$individual_ID), " individuals")

#==============================================================================
# 1. ESTIMATE LOCATION ERROR FROM REFERENCE TAG
#==============================================================================

# Extract reference tag detections ------------------------------------------
ref <- BT_sub_dt[individual_ID == "FReference"]
message("\n=== Reference Tag Analysis ===")
message("Reference tag detections: ", nrow(ref)) #88576

# Known fixed coordinates of reference tag ---------------------------------
# Reference receiver number
unique(ref$FullId)
#"H170-1802-65065"

# Get the coordinates for the location of the reference receiver
rec_locs <- fread(paste0(rec_loc_path, "BT_rec_hab_locations.csv"))
ref_coords <- rec_locs %>% filter(ID == 'ref')

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
true_lat <- parse_dms_to_decimal(ref_coords$Latitude)
true_long <- parse_dms_to_decimal(ref_coords$Longitude)

message("Reference tag coordinates:")
message("  Latitude: ", round(true_lat, 7), "°N (", ref_coords$Latitude, ")")
message("  Longitude: ", round(true_long, 7), "°E (", ref_coords$Longitude, ")")

# Calculate distance between each detection and true position (m) -----------
ref[, error_dist_m := distHaversine(
  cbind(Long, Lat),
  cbind(true_long, true_lat)
)]

# Summarize positional error ------------------------------------------------
message("\n=== Reference Tag Location Error (meters) ===")
print(summary(ref$error_dist_m))

sigma_error_med <- median(ref$error_dist_m, na.rm = TRUE)
message("Median error: ", round(sigma_error_med, 4), " m")
#Median error: 0.1091 m

#==============================================================================
# 2. ANALYZE MOVEMENT SCALE VS TIME LAG
#==============================================================================
# Purpose: Determine optimal temporal thinning interval where movement
# is reliably distinguishable from location error

# Select representative individuals (2 per species) -------------------------
set.seed(123)  # For reproducibility
selected_ids <- BT_sub_dt[
  Species != "FReference",
  .(individual_ID = sample(unique(individual_ID), min(2, .N))),
  by = Species
]$individual_ID

subset_ids <- BT_sub_dt[individual_ID %in% selected_ids]
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
step_summary[, movement_to_error := median_step / BT_UERE$UERE[,1]]

message("\n=== Movement vs Time Lag Analysis ===")
print(step_summary)

message("\nConclusion: Movement becomes most distinguishable from error at 30s+")

#==============================================================================
# 3. PREPARE DATA FOR OUTLIER DETECTION
#==============================================================================

# Convert to Movebank format for ctmm package ------------------------------
BT_movebank <- with(
  BT_sub,
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
BT_tels <- as.telemetry(
  BT_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,
  datum = "WGS84",
  keep = c("species", "weight", "total_length", "std_length", "treatment",
           "date", "exp_stage", "time_of_day", "found_alive", "known_predated")
)

# Center projection on geometric median -------------------------------------
projection(BT_tels) <- ctmm::median(BT_tels)

#==============================================================================
# 4. INCORPORATE LOCATION ERROR MODEL
#==============================================================================

# Fit error parameters using reference tag ----------------------------------
BT_UERE <- uere.fit(BT_tels$FReference)
print(summary(BT_UERE))
# low      est      high
# all 0.37805 0.379299 0.3805479

# Apply error model to all telemetry objects -------------------------------
uere(BT_tels) <- BT_UERE

# Remove reference tag from analysis ---------------------------------------
# Note: Adjust the index range based on actual number of individuals
n_individuals <- length(BT_tels) - 1
BT_tels <- BT_tels[1:n_individuals]
message("Removed reference tag. Remaining individuals: ", length(BT_tels))

# Save UERE model -----------------------------------------------------------
saveRDS(BT_UERE, paste0(telem_path, "BT_UERE.rds"))

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
species_order <- BT_movebank %>%
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
pike_BT_tel <- BT_tels[pike_indices]
perch_BT_tel <- BT_tels[perch_indices]
roach_BT_tel <- BT_tels[roach_indices]

message("\n=== Species Telemetry Objects Created ===")
message("Pike: ", length(pike_BT_tel), " individuals (indices: ", 
        paste(pike_indices, collapse = ", "), ")")
message("Perch: ", length(perch_BT_tel), " individuals (indices: ", 
        paste(perch_indices, collapse = ", "), ")")
message("Roach: ", length(roach_BT_tel), " individuals (indices: ", 
        paste(roach_indices, collapse = ", "), ")")

# Verify totals match
total_individuals <- length(pike_BT_tel) + length(perch_BT_tel) + length(roach_BT_tel)
message("\nTotal individuals across all species: ", total_individuals)
message("Expected total (excluding reference tag): ", length(BT_tels))

if(total_individuals != length(BT_tels)) {
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
pike_BT_tel <- remove_speed_outliers(
  pike_BT_tel, "Pike", speed_thresholds$Pike
) #outliers detected: 9882
saveRDS(pike_BT_tel, paste0(telem_path, "pike_BT_tel_unthinned.rds"))

perch_BT_tel <- remove_speed_outliers(
  perch_BT_tel, "Perch", speed_thresholds$Perch
) #outliers detected: 33878
saveRDS(perch_BT_tel, paste0(telem_path, "perch_BT_tel_unthinned.rds"))

roach_BT_tel <- remove_speed_outliers(
  roach_BT_tel, "Roach", speed_thresholds$Roach
) #outliers detected: 339565
saveRDS(roach_BT_tel, paste0(telem_path, "roach_BT_tel_unthinned.rds"))

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
pike_data <- extract_data(pike_BT_tel)
perch_data <- extract_data(perch_BT_tel)
roach_data <- extract_data(roach_BT_tel)

# Combine into single dataframe --------------------------------------------
BT_telem_data <- rbind(
  cbind(pike_data, Species = 'Pike'),
  cbind(perch_data, Species = 'Perch'),
  cbind(roach_data, Species = 'Roach')
)

message("\n=== Data After Outlier Removal ===")
message("Total detections: ", nrow(BT_telem_data)) #Total detections: 19165902
message("Outliers removed: ", nrow(BT_sub_dt) - nrow(ref) - nrow(BT_telem_data)) #Outliers removed: 158020

# Save outlier-filtered data ------------------------------------------------
saveRDS(BT_telem_data, paste0(filtered_data_path, "02_BT_sub.rds"))

# Clean up large objects ----------------------------------------------------
rm(BT_movebank, BT_sub, BT_sub_dt, BT_tels,
   roach_data, perch_data, pike_data,
   roach_BT_tel, pike_BT_tel, perch_BT_tel,
   subset_ids, ref)
gc()

#==============================================================================
# 7. TEMPORAL THINNING
#==============================================================================

BT_telem_data <- readRDS(paste0(filtered_data_path, "02_BT_sub.rds"))

# Convert to move2 object ---------------------------------------------------
BT_mv <- mt_as_move2(
  BT_telem_data,
  coords = c("longitude", "latitude"),
  crs = "WGS84",
  time_column = "timestamp",
  track_id_column = "individual_ID",
  na.fail = FALSE
)

BT_mv <- BT_mv %>%
  arrange(individual_ID, timestamp)

# Apply temporal thinning ---------------------------------------------------

#15 seconds
# thinning_interval <- "15 seconds"
# message("\n=== Temporal Thinning ===")
# message("Applying interval: ", thinning_interval)
# 
# BT_thin_15_data <- BT_mv %>%
#   mt_filter_per_interval(unit = thinning_interval, criterion = "first")
# 
# message("Before thinning: ", nrow(BT_mv), " detections") #10459343 detections
# message("After thinning: ", nrow(BT_thin_data), " detections") #4952567 detections
# message("Reduction: ", round(100 * (1 - nrow(BT_thin_data) / nrow(BT_mv)), 1), "%") #Reduction: 52.6%

#30 seconds
thinning_interval <- "30 seconds"
message("\n=== Temporal Thinning ===")
message("Applying interval: ", thinning_interval)

BT_thin_data <- BT_mv %>%
  mt_filter_per_interval(unit = thinning_interval, criterion = "first")

message("Before thinning: ", nrow(BT_mv), " detections") #19165902 detections
message("After thinning: ", nrow(BT_thin_data), " detections") #After thinning: 3027526 detections
message("Reduction: ", round(100 * (1 - nrow(BT_thin_data) / nrow(BT_mv)), 1), "%") #Reduction: 84.2%%

# Verify thinning intervals -------------------------------------------------
BT_thin_track_sum <- BT_thin_data %>%
  group_by(individual_ID) %>%
  arrange(timestamp) %>%
  mutate(dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
  summarise(
    min_dt = min(dt, na.rm = TRUE),
    median_dt = median(dt, na.rm = TRUE),
    max_dt = max(dt, na.rm = TRUE)
  )

message("\n=== Temporal Resolution Summary ===")
print(summary(BT_thin_track_sum[, c("min_dt", "median_dt", "max_dt")]))
# min_dt         median_dt         max_dt                 geometry 
# Min.   :0.7998   Min.   :30.05   Min.   :   1249   MULTIPOINT   :65  
# 1st Qu.:1.8003   1st Qu.:30.17   1st Qu.:  29614   epsg:4326    : 0  
# Median :1.8007   Median :30.39   Median :  48367   +proj=long...: 0  
# Mean   :1.7267   Mean   :30.54   Mean   : 103135                     
# 3rd Qu.:1.8009   3rd Qu.:30.69   3rd Qu.:  87910                     
# Max.   :1.8254   Max.   :32.46   Max.   :1300537     

#==============================================================================
# 8. IDENTIFY TEMPORAL GAPS IN TRACKING
#==============================================================================

# Convert to data.table and calculate time gaps ----------------------------
BT_thin_data <- as.data.table(BT_thin_data)
setorder(BT_thin_data, individual_ID, timestamp)

# Calculate time gaps between successive fixes (minutes)
BT_thin_data[, dt_mins := as.numeric(
  difftime(timestamp, data.table::shift(timestamp), units = "mins")
), by = individual_ID]

# Store previous timestamp and date (start of gap)
BT_thin_data[, `:=`(
  prev_time = data.table::shift(timestamp),
  prev_date = data.table::shift(date)
), by = individual_ID]

# Find maximum gap for each individual --------------------------------------
max_gaps <- BT_thin_data[
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
output_file <- paste0(filtered_data_path, "BT_temporal_gaps_summary.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

message("\nTemporal gaps summary saved to Excel: ", output_file)
message("Rows with gaps > 24 hours are highlighted in red")

# Extract coordinates from geometry -----------------------------------------
coords_data <- st_coordinates(BT_thin_data$geometry)
BT_thin_data$Long <- coords_data[, 1]
BT_thin_data$Lat <- coords_data[, 2]

#Remove unneeded columns
BT_thin_data <- BT_thin_data %>% 
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

traj_plots <- BT_thin_data %>%
  group_by(individual_ID) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$individual_ID))) %>%
  map(~ plot_daily_traj_individual(.x, BT_polygon))

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


# Save final filtered dataset -----------------------------------------------
saveRDS(BT_thin_data, paste0(filtered_data_path, "03_BT_sub.rds"))

message("\n=== Data Processing Complete ===")
message("Final dataset: ", nrow(BT_thin_data), " detections") #Final dataset: 3027526 detections
message("Number of individuals: ", n_distinct(BT_thin_data$individual_ID)) 
message("Date range: ", min(BT_thin_data$date), " to ", 
        max(BT_thin_data$date)) #Date range: 2022-09-26 to 2022-10-29
message("Saved to: ", paste0(filtered_data_path, "03_BT_sub.rds"))

#==============================================================================
# SUMMARY OF OUTPUT FILES
#==============================================================================
message("\n=== Output Files Generated ===")
message("1. Lake boundary coordinates: ", paste0(lake_polygon_path, "BT_boundary_coordinates.csv"))
message("2. Lake data summary: ", paste0(lake_polygon_path, "BT_lake_data.rds"))
message("3. UERE model: ", paste0(save_telem_path, "BT/BT_UERE.rds"))
message("4. Outlier-filtered data: ", paste0(filtered_data_path, "02_BT_sub.rds"))
message("5. Final thinned data: ", paste0(filtered_data_path, "03_BT_sub.rds"))
message("6. Daily trajectory plots: ", plot_output_path)

















#==============================================================================
# Lake BT Detection Data - Outlier Removal and Temporal Thinning
#==============================================================================
# Purpose: Remove outliers, thin data, and filter mortality events
# Author: Marcus Michelangeli
#
# Input files:
#   - ./data/tracks_filtered/lake_BT/01_lake_BT_sub.rds
#   - ./data/lake_coords/lake_BT_polygon.gpkg
#   - ./data/pike_deaths.csv
#
# Output:
#   - ./data/telem_obj/BT/lake_BT_UERE.rds
#   - ./data/telem_obj/BT/[species]_lake_BT_tel_unthinned.rds (3 files)
#   - ./data/tracks_filtered/lake_BT/02_lake_BT_sub.rds
#   - ./data/tracks_filtered/lake_BT/03_lake_BT_sub.rds
#   - ./daily_trajectory_plots/lake_BT/[individual]_daily_trajectory_plot.png
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
filtered_data_path <- "./data/tracks_filtered/lake_BT/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"
lake_polygon_path <- "./data/lake_coords/"
plot_output_path <- "./daily_trajectory_plots/lake_BT/"

# Import data ---------------------------------------------------------------
lake_BT_sub <- readRDS(paste0(filtered_data_path, '01_lake_BT_sub.rds'))
lake_BT_sub_dt <- as.data.table(lake_BT_sub)
BT_polygon <- st_read(paste0(lake_polygon_path, "lake_BT_polygon.gpkg"))

message("Imported ", nrow(lake_BT_sub_dt), " detections for ", 
        n_distinct(lake_BT_sub_dt$individual_ID), " individuals")

#==============================================================================
# 1. ESTIMATE LOCATION ERROR FROM REFERENCE TAG
#==============================================================================

# Extract reference tag detections ------------------------------------------
ref <- lake_BT_sub_dt[individual_ID == "FReference"]
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
selected_ids <- lake_BT_sub_dt[
  Species != "FReference",
  .(individual_ID = sample(unique(individual_ID), min(2, .N))),
  by = Species
]$individual_ID

subset_ids <- lake_BT_sub_dt[individual_ID %in% selected_ids]
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
lake_BT_movebank <- with(
  lake_BT_sub,
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
lake_BT_tels <- as.telemetry(
  lake_BT_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,
  datum = "WGS84",
  keep = c("Species", "Weight", "Total_length", "Std_length", "Treatment",
           "Date", "Exp_Stage", "Time_Of_Day", "found_alive", "known_predated")
)

# Center projection on geometric median -------------------------------------
projection(lake_BT_tels) <- ctmm::median(lake_BT_tels)

#==============================================================================
# 4. INCORPORATE LOCATION ERROR MODEL
#==============================================================================

# Fit error parameters using reference tag ----------------------------------
lake_BT_UERE <- uere.fit(lake_BT_tels$FReference)
message("\n=== UERE Model Summary ===")
print(summary(lake_BT_UERE))
# low      est      high
# all 0.37805 0.379299 0.3805479

# Apply error model to all telemetry objects -------------------------------
uere(lake_BT_tels) <- lake_BT_UERE

# Remove reference tag from analysis ---------------------------------------
lake_BT_tels <- lake_BT_tels[1:66]
message("Removed reference tag. Remaining individuals: ", length(lake_BT_tels))

# Save UERE model -----------------------------------------------------------
saveRDS(lake_BT_UERE, paste0(save_telem_path, "BT/lake_BT_UERE.rds"))

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
pike_lake_BT_tel <- lake_BT_tels[61:66]
perch_lake_BT_tel <- lake_BT_tels[c(1:15, 31:44, 46)]
roach_lake_BT_tel <- lake_BT_tels[c(16:30, 45, 47:60)]

message("\nSpecies distribution:")
message("Pike: ", length(pike_lake_BT_tel), " individuals")
message("Perch: ", length(perch_lake_BT_tel), " individuals")
message("Roach: ", length(roach_lake_BT_tel), " individuals")

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
pike_lake_BT_tel <- remove_speed_outliers(
  pike_lake_BT_tel, "Pike", speed_thresholds$Pike
) #outliers removed: 11538
saveRDS(pike_lake_BT_tel, paste0(save_telem_path, "BT/pike_lake_BT_tel_unthinned.rds"))

perch_lake_BT_tel <- remove_speed_outliers(
  perch_lake_BT_tel, "Perch", speed_thresholds$Perch
) #outliers removed: 33878
saveRDS(perch_lake_BT_tel, paste0(save_telem_path, "BT/perch_lake_BT_tel_unthinned.rds"))

roach_lake_BT_tel <- remove_speed_outliers(
  roach_lake_BT_tel, "Roach", speed_thresholds$Roach
) #outlier removed: 114260
saveRDS(roach_lake_BT_tel, paste0(save_telem_path, "BT/roach_lake_BT_tel_unthinned.rds"))

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
pike_data <- extract_data(pike_lake_BT_tel)
perch_data <- extract_data(perch_lake_BT_tel)
roach_data <- extract_data(roach_lake_BT_tel)

# Combine into single dataframe --------------------------------------------
BT_telem_data <- rbind(
  cbind(pike_data, Species = 'Pike'),
  cbind(perch_data, Species = 'Perch'),
  cbind(roach_data, Species = 'Roach')
)

message("\n=== Data After Outlier Removal ===")
message("Total detections: ", nrow(BT_telem_data)) #20019313
message("Outliers removed: ", nrow(lake_BT_sub_dt) - nrow(ref) - nrow(BT_telem_data)) #159676

# Save outlier-filtered data ------------------------------------------------
saveRDS(BT_telem_data, paste0(filtered_data_path, "02_lake_BT_sub.rds"))

# Clean up large objects ----------------------------------------------------
rm(lake_BT_movebank, lake_BT_sub, lake_BT_sub_dt, lake_BT_tels,
   roach_data, perch_data, pike_data,
   roach_lake_BT_tel, pike_lake_BT_tel, perch_lake_BT_tel,
   subset_ids, ref)
gc()

#==============================================================================
# 7. TEMPORAL THINNING
#==============================================================================

# Convert to move2 object ---------------------------------------------------
lake_BT_mv <- mt_as_move2(
  BT_telem_data,
  coords = c("longitude", "latitude"),
  crs = "WGS84",
  time_column = "timestamp",
  track_id_column = "individual_ID",
  na.fail = FALSE
)

lake_BT_mv <- lake_BT_mv %>%
  arrange(individual_ID, timestamp)

# Apply temporal thinning ---------------------------------------------------
thinning_interval <- "10 seconds"
message("\nApplying temporal thinning: ", thinning_interval)

lake_BT_thin_data <- lake_BT_mv %>%
  mt_filter_per_interval(unit = thinning_interval, criterion = "first")

message("Before thinning: ", nrow(lake_BT_mv), " detections") #20019313 detections
message("After thinning: ", nrow(lake_BT_thin_data), " detections") #7339440 detections
message("Reduction: ", round(100 * (1 - nrow(lake_BT_thin_data) / nrow(lake_BT_mv)), 1), "%") #Reduction: 63.3%

# Verify thinning intervals -------------------------------------------------
lake_BT_thin_track_sum <- lake_BT_thin_data %>%
  group_by(individual_ID) %>%
  arrange(timestamp) %>%
  mutate(dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
  summarise(
    min_dt = min(dt, na.rm = TRUE),
    median_dt = median(dt, na.rm = TRUE),
    max_dt = max(dt, na.rm = TRUE)
  )

message("\n=== Temporal Resolution Summary ===")
print(summary(lake_BT_thin_track_sum[, c("min_dt", "median_dt", "max_dt")]))

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
lake_BT_thin_data <- as.data.table(lake_BT_thin_data)
setorder(lake_BT_thin_data, individual_ID, timestamp)

# Calculate time gaps between successive fixes (minutes)
lake_BT_thin_data[, dt_mins := as.numeric(
  difftime(timestamp, shift(timestamp), units = "mins")
), by = individual_ID]

# Store previous timestamp and date (start of gap)
lake_BT_thin_data[, `:=`(
  prev_time = shift(timestamp),
  prev_date = shift(Date)
), by = individual_ID]

# Find maximum gap for each individual --------------------------------------
max_gaps <- lake_BT_thin_data[
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
output_file <- paste0(filtered_data_path, "lake_BT_temporal_gaps_summary.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

message("\nTemporal gaps summary saved to Excel: ", output_file)
message("Rows with gaps > 24 hours are highlighted in red")

# Extract coordinates from geometry -----------------------------------------
coords_data <- st_coordinates(lake_BT_thin_data$geometry)
lake_BT_thin_data$Long <- coords_data[, 1]
lake_BT_thin_data$Lat <- coords_data[, 2]

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

traj_plots <- lake_BT_thin_data %>%
  group_by(individual_ID) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$individual_ID))) %>%
  map(~ plot_daily_traj_individual(.x, BT_polygon))

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
n_before <- n_distinct(lake_BT_thin_data$individual_ID)
lake_BT_thin_data <- lake_BT_thin_data %>%
  filter(individual_ID != "F59889")

message("\nRemoved F59889 (died on Day 1)")
message("Individuals remaining: ", n_distinct(lake_BT_thin_data$individual_ID))

# Load pike mortality information -------------------------------------------
pike_deaths <- read.csv("./data/pike_deaths.csv")

pike_mort_cols <- pike_deaths %>%
  filter(individual_ID %in% lake_BT_thin_data$individual_ID) %>%
  mutate(pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y"))

message("\nPike mortality records loaded: ", nrow(pike_mort_cols))

# Filter out detections after mortality events ------------------------------
n_rows_before <- nrow(lake_BT_thin_data)

lake_BT_thin_data <- lake_BT_thin_data %>%
  left_join(pike_mort_cols, by = "individual_ID") %>%
  filter(is.na(pike_death_date) | Date < pike_death_date)

n_rows_after <- nrow(lake_BT_thin_data)
message("Detections removed after mortality events: ", n_rows_before - n_rows_after) #48037

# Save final filtered dataset -----------------------------------------------
saveRDS(lake_BT_thin_data, paste0(filtered_data_path, "03_lake_BT_sub.rds"))

message("Final dataset: ", nrow(lake_BT_thin_data), " detections") #7143462 detections
message("Number of individuals: ", n_distinct(lake_BT_thin_data$individual_ID)) #Number of individuals: 65
message("Date range: ", min(lake_BT_thin_data$Date), " to ", 
        max(lake_BT_thin_data$Date))
message("Saved to: ", paste0(filtered_data_path, "03_lake_BT_sub.rds"))
