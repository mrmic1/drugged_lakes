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
