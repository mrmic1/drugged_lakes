#==============================================================================
# Lake Muddyfoot Detection Data - Outlier Removal and Temporal Thinning
#==============================================================================
# Purpose: Remove outliers, thin data, and filter mortality events
# Author: Marcus Michelangeli
#
# Input files:
#   - ./data/tracks_filtered/muddyfoot/01_muddyfoot_sub.rds
#   - ./data/lake_coords/lake_muddyfoot_polygon.gpkg
#   - ./data/pike_deaths.csv
#
# Output:
#   - ./data/telem_obj/muddyfoot/muddyfoot_UERE.rds
#   - ./data/telem_obj/muddyfoot/[species]_muddyfoot_tel_unthinned.rds (3 files)
#   - ./data/tracks_filtered/muddyfoot/02_muddyfoot_sub.rds
#   - ./data/tracks_filtered/muddyfoot/03_muddyfoot_sub.rds
#   - ./data/lake_coords/muddyfoot_boundary_coordinates.csv (NEW)
#   - ./data/lake_coords/muddyfoot_lake_data.rds (NEW)
#   - ./daily_trajectory_plots/muddyfoot/[individual]_daily_trajectory_plot.png
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
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/muddyfoot/"
polygon_path <- "./data/lake_params/polygons/"
rec_loc_path <- "./data/lake_params/reciever_and_habitat_locations/"
plot_output_path <- "./daily_trajectory_plots/muddyfoot/"

# Import data ---------------------------------------------------------------
muddyfoot_sub <- readRDS(paste0(filtered_data_path, '01_muddyfoot_sub.rds'))
muddyfoot_sub_dt <- as.data.table(muddyfoot_sub)
muddyfoot_polygon <- st_read(paste0(polygon_path, "muddyfoot_polygon.gpkg"))

message("Imported ", nrow(muddyfoot_sub_dt), " detections for ", 
        n_distinct(muddyfoot_sub_dt$individual_ID), " individuals")

#==============================================================================
# 1. ESTIMATE LOCATION ERROR FROM REFERENCE TAG
#==============================================================================

# Extract reference tag detections ------------------------------------------
ref <- muddyfoot_sub_dt[individual_ID == "FReference"]
message("\n=== Reference Tag Analysis ===")
message("Reference tag detections: ", nrow(ref)) #22467

# Known fixed coordinates of reference tag ---------------------------------
# Reference receiver number
unique(ref$FullId)
#"H170-1802-65066"

# Get the coordinates for the location of the reference receiver
rec_locs <- fread(paste0(rec_loc_path, "muddyfoot_rec_hab_locations.csv"))
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


#==============================================================================
# 2. ANALYZE MOVEMENT SCALE VS TIME LAG
#==============================================================================
# Purpose: Determine optimal temporal thinning interval where movement
# is reliably distinguishable from location error

# Select representative individuals (2 per species) -------------------------
set.seed(123)  # For reproducibility
selected_ids <- muddyfoot_sub_dt[
  Species != "FReference",
  .(individual_ID = sample(unique(individual_ID), min(2, .N))),
  by = Species
]$individual_ID

subset_ids <- muddyfoot_sub_dt[individual_ID %in% selected_ids]
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
step_summary[, movement_to_error := median_step / sigma_error_med]

message("\n=== Movement vs Time Lag Analysis ===")
print(step_summary)

message("\nConclusion: Movement becomes most distinguishable from error at 30s+")

#==============================================================================
# 3. PREPARE DATA FOR OUTLIER DETECTION
#==============================================================================

# Convert to Movebank format for ctmm package ------------------------------
muddyfoot_movebank <- with(
  muddyfoot_sub,
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
muddyfoot_tels <- as.telemetry(
  muddyfoot_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,
  datum = "WGS84",
  keep = c("species", "weight", "total_length", "std_length", "treatment",
           "date", "exp_stage", "time_of_day", "found_alive", "known_predated")
)

# Center projection on geometric median -------------------------------------
projection(muddyfoot_tels) <- ctmm::median(muddyfoot_tels)

#==============================================================================
# 4. INCORPORATE LOCATION ERROR MODEL
#==============================================================================

# Fit error parameters using reference tag ----------------------------------
muddyfoot_UERE <- uere.fit(muddyfoot_tels$FReference)
print(summary(muddyfoot_UERE))
# low       est      high
# all 0.5133307 0.5167092 0.5200873

# Apply error model to all telemetry objects -------------------------------
uere(muddyfoot_tels) <- muddyfoot_UERE

# Remove reference tag from analysis ---------------------------------------
# Note: Adjust the index range based on actual number of individuals
n_individuals <- length(muddyfoot_tels) - 1
muddyfoot_tels <- muddyfoot_tels[1:n_individuals]
message("Removed reference tag. Remaining individuals: ", length(muddyfoot_tels))

# Save UERE model -----------------------------------------------------------
saveRDS(muddyfoot_UERE, paste0(telem_path, "muddyfoot_UERE.rds"))

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
species_order <- muddyfoot_movebank %>%
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
pike_muddyfoot_tel <- muddyfoot_tels[pike_indices]
perch_muddyfoot_tel <- muddyfoot_tels[perch_indices]
roach_muddyfoot_tel <- muddyfoot_tels[roach_indices]

message("\n=== Species Telemetry Objects Created ===")
message("Pike: ", length(pike_muddyfoot_tel), " individuals (indices: ", 
        paste(pike_indices, collapse = ", "), ")")
message("Perch: ", length(perch_muddyfoot_tel), " individuals (indices: ", 
        paste(perch_indices, collapse = ", "), ")")
message("Roach: ", length(roach_muddyfoot_tel), " individuals (indices: ", 
        paste(roach_indices, collapse = ", "), ")")

# Verify totals match
total_individuals <- length(pike_muddyfoot_tel) + length(perch_muddyfoot_tel) + length(roach_muddyfoot_tel)
message("\nTotal individuals across all species: ", total_individuals)
message("Expected total (excluding reference tag): ", length(muddyfoot_tels))

if(total_individuals != length(muddyfoot_tels)) {
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
pike_muddyfoot_tel <- remove_speed_outliers(
  pike_muddyfoot_tel, "Pike", speed_thresholds$Pike
) #outliers detected: 122566
saveRDS(pike_muddyfoot_tel, paste0(telem_path, "pike_muddyfoot_tel_unthinned.rds"))

perch_muddyfoot_tel <- remove_speed_outliers(
  perch_muddyfoot_tel, "Perch", speed_thresholds$Perch
) #outliers detected: 154699
saveRDS(perch_muddyfoot_tel, paste0(telem_path, "perch_muddyfoot_tel_unthinned.rds"))

roach_muddyfoot_tel <- remove_speed_outliers(
  roach_muddyfoot_tel, "Roach", speed_thresholds$Roach
) #outliers detected: 339565
saveRDS(roach_muddyfoot_tel, paste0(telem_path, "roach_muddyfoot_tel_unthinned.rds"))

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
pike_data <- extract_data(pike_muddyfoot_tel)
perch_data <- extract_data(perch_muddyfoot_tel)
roach_data <- extract_data(roach_muddyfoot_tel)

# Combine into single dataframe --------------------------------------------
muddyfoot_telem_data <- rbind(
  cbind(pike_data, Species = 'Pike'),
  cbind(perch_data, Species = 'Perch'),
  cbind(roach_data, Species = 'Roach')
)

message("\n=== Data After Outlier Removal ===")
message("Total detections: ", nrow(muddyfoot_telem_data)) #10459343
message("Outliers removed: ", nrow(muddyfoot_sub_dt) - nrow(ref) - nrow(muddyfoot_telem_data)) #616830

# Save outlier-filtered data ------------------------------------------------
saveRDS(muddyfoot_telem_data, paste0(filtered_data_path, "02_muddyfoot_sub.rds"))

# Clean up large objects ----------------------------------------------------
rm(muddyfoot_movebank, muddyfoot_sub, muddyfoot_sub_dt, muddyfoot_tels,
   roach_data, perch_data, pike_data,
   roach_muddyfoot_tel, pike_muddyfoot_tel, perch_muddyfoot_tel,
   subset_ids, ref)
gc()

#==============================================================================
# 7. TEMPORAL THINNING
#==============================================================================

muddyfoot_telem_data <- readRDS(paste0(filtered_data_path, "02_muddyfoot_sub.rds"))

# Convert to move2 object ---------------------------------------------------
muddyfoot_mv <- mt_as_move2(
  muddyfoot_telem_data,
  coords = c("longitude", "latitude"),
  crs = "WGS84",
  time_column = "timestamp",
  track_id_column = "individual_ID",
  na.fail = FALSE
)

muddyfoot_mv <- muddyfoot_mv %>%
  arrange(individual_ID, timestamp)

# Apply temporal thinning ---------------------------------------------------

#15 seconds
# thinning_interval <- "15 seconds"
# message("\n=== Temporal Thinning ===")
# message("Applying interval: ", thinning_interval)
# 
# muddyfoot_thin_15_data <- muddyfoot_mv %>%
#   mt_filter_per_interval(unit = thinning_interval, criterion = "first")
# 
# message("Before thinning: ", nrow(muddyfoot_mv), " detections") #10459343 detections
# message("After thinning: ", nrow(muddyfoot_thin_data), " detections") #4952567 detections
# message("Reduction: ", round(100 * (1 - nrow(muddyfoot_thin_data) / nrow(muddyfoot_mv)), 1), "%") #Reduction: 52.6%

#30 seconds
thinning_interval <- "30 seconds"
message("\n=== Temporal Thinning ===")
message("Applying interval: ", thinning_interval)

muddyfoot_thin_data <- muddyfoot_mv %>%
  mt_filter_per_interval(unit = thinning_interval, criterion = "first")

message("Before thinning: ", nrow(muddyfoot_mv), " detections") #10459343 detections
message("After thinning: ", nrow(muddyfoot_thin_data), " detections") #2500625 detections
message("Reduction: ", round(100 * (1 - nrow(muddyfoot_thin_data) / nrow(muddyfoot_mv)), 1), "%") #Reduction: 76.1%

# Verify thinning intervals -------------------------------------------------
muddyfoot_thin_track_sum <- muddyfoot_thin_data %>%
  group_by(individual_ID) %>%
  arrange(timestamp) %>%
  mutate(dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
  summarise(
    min_dt = min(dt, na.rm = TRUE),
    median_dt = median(dt, na.rm = TRUE),
    max_dt = max(dt, na.rm = TRUE)
  )

message("\n=== Temporal Resolution Summary ===")
print(summary(muddyfoot_thin_track_sum[, c("min_dt", "median_dt", "max_dt")]))
# min_dt         median_dt         max_dt                geometry 
# Min.   :0.7998   Min.   :30.29   Min.   :  3700   MULTIPOINT   :65  
# 1st Qu.:1.8003   1st Qu.:31.16   1st Qu.: 32061   epsg:4326    : 0  
# Median :1.8007   Median :31.58   Median : 50413   +proj=long...: 0  
# Mean   :1.7086   Mean   :32.03   Mean   :100902                     
# 3rd Qu.:1.8009   3rd Qu.:32.51   3rd Qu.:120338                     
# Max.   :1.8254   Max.   :37.55   Max.   :698753   

#==============================================================================
# 8. IDENTIFY TEMPORAL GAPS IN TRACKING
#==============================================================================

# Convert to data.table and calculate time gaps ----------------------------
muddyfoot_thin_data <- as.data.table(muddyfoot_thin_data)
setorder(muddyfoot_thin_data, individual_ID, timestamp)

# Calculate time gaps between successive fixes (minutes)
muddyfoot_thin_data[, dt_mins := as.numeric(
  difftime(timestamp, data.table::shift(timestamp), units = "mins")
), by = individual_ID]

# Store previous timestamp and date (start of gap)
muddyfoot_thin_data[, `:=`(
  prev_time = data.table::shift(timestamp),
  prev_date = data.table::shift(date)
), by = individual_ID]

# Find maximum gap for each individual --------------------------------------
max_gaps <- muddyfoot_thin_data[
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
output_file <- paste0(filtered_data_path, "muddyfoot_temporal_gaps_summary.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

message("\nTemporal gaps summary saved to Excel: ", output_file)
message("Rows with gaps > 24 hours are highlighted in red")

# Extract coordinates from geometry -----------------------------------------
coords_data <- st_coordinates(muddyfoot_thin_data$geometry)
muddyfoot_thin_data$Long <- coords_data[, 1]
muddyfoot_thin_data$Lat <- coords_data[, 2]

#Remove unneeded columns
muddyfoot_thin_data <- muddyfoot_thin_data %>% 
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

traj_plots <- muddyfoot_thin_data %>%
  group_by(individual_ID) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$individual_ID))) %>%
  map(~ plot_daily_traj_individual(.x, muddyfoot_polygon))

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
saveRDS(muddyfoot_thin_data, paste0(filtered_data_path, "03_muddyfoot_sub.rds"))

message("\n=== Data Processing Complete ===")
message("Final dataset: ", nrow(muddyfoot_thin_data), " detections") 
message("Number of individuals: ", n_distinct(muddyfoot_thin_data$individual_ID)) 
message("Date range: ", min(muddyfoot_thin_data$Date), " to ", 
        max(muddyfoot_thin_data$Date)) #Date range: 2022-09-25 to 2022-10-30
message("Saved to: ", paste0(filtered_data_path, "03_muddyfoot_sub.rds"))

#==============================================================================
# SUMMARY OF OUTPUT FILES
#==============================================================================
message("\n=== Output Files Generated ===")
message("1. Lake boundary coordinates: ", paste0(lake_polygon_path, "muddyfoot_boundary_coordinates.csv"))
message("2. Lake data summary: ", paste0(lake_polygon_path, "muddyfoot_lake_data.rds"))
message("3. UERE model: ", paste0(save_telem_path, "muddyfoot/muddyfoot_UERE.rds"))
message("4. Outlier-filtered data: ", paste0(filtered_data_path, "02_muddyfoot_sub.rds"))
message("5. Final thinned data: ", paste0(filtered_data_path, "03_muddyfoot_sub.rds"))
message("6. Daily trajectory plots: ", plot_output_path)