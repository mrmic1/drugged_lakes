### LIBRARIES ###
# Load necessary libraries
library(data.table)   # For fast data manipulation
library(tidyverse)    # For data wrangling
library(ctmm)         # For continuous-time movement modeling
library(sf)           # For handling spatial data
library(parallel)     # For parallel processing
library(foreach)      # For parallel for loops
library(doParallel)   # For registering parallel backend
library(geosphere)
library(move2)

# Set the time zone to ensure consistent time handling
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths for reading and saving filtered telemetry and ctmm model results
filtered_data_path <- "./data/tracks_filtered/lake_BT/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"
lake_polygon_path <- "./data/lake_coords/"


#Load in the datasets
lake_BT_sub <- readRDS(paste0(filtered_data_path, '01_lake_BT_sub.rds'))
lake_BT_sub_dt <- as.data.table(lake_BT_sub)
BT_polygon <- sf::st_read(paste0(lake_polygon_path, "lake_BT_polygon.gpkg"))


#----------------------------------------------------------------------------------------#

# Get an error scale from the reference tag


ref <- lake_BT_sub_dt[individual_ID == "FReference"]

#Known fixed coordinates of the reference tag
true_long <- 20.0472989
true_lat <-  63.7707882

# Distance between each detection and the true position (m)
ref[, error_dist_m := distHaversine(
  cbind(Long, Lat),
  cbind(true_long, true_lat)
)]

summary(ref$error_dist_m)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003683 0.0810723 0.1092805 0.1117890 0.1363092 3.8638857

sigma_error_med  <- median(ref$error_dist_m, na.rm = TRUE)

#-----------------------------------------------------------------------------------#

# Look at how step lengths grow with time lag for fish


selected_ids <- lake_BT_sub_dt[
  Species != "FReference",
  .(individual_ID = sample(unique(individual_ID), min(2, .N))),
  by = Species
]$individual_ID

subset_ids <- lake_BT_sub_dt[individual_ID %in% selected_ids]

# Ensure sorted
setorder(subset_ids, individual_ID, timestamp_cest)

# Time difference between successive fixes (seconds)
subset_ids[, dt_sec := as.numeric(
  difftime(timestamp_cest, shift(timestamp_cest), units = "secs")
), by = individual_ID]

# Step length between successive fixes (m) using Haversine distance
subset_ids[, step_dist_m := distHaversine(
  cbind(shift(Long), shift(Lat)),
  cbind(Long, Lat)
), by = individual_ID]

#----------------------------------------------------------------------------------------#

# Bin by time lage and compare movement to error

# Define Δt bins (seconds)
breaks <- c(0, 3, 7, 12, 18, 25, 40, 70, Inf)
labels <- c("~2", "~5", "~10", "~15", "~20", "~30", ">30", ">70")

subset_ids[!is.na(step_dist_m),
           dt_bin := cut(dt_sec,
                         breaks = breaks,
                         labels = labels,
                         include.lowest = TRUE)]

step_summary <- subset_ids[!is.na(step_dist_m),
                         .(
                           n_steps      = .N,
                           median_step  = median(step_dist_m),
                           mean_step    = mean(step_dist_m),
                           q25_step     = quantile(step_dist_m, 0.25),
                           q75_step     = quantile(step_dist_m, 0.75)
                         ),
                         by = dt_bin
]

# enforce factor level order
step_summary[, dt_bin := factor(dt_bin, levels = labels)]

# order the summary table by dt_bin
setorder(step_summary, dt_bin)

step_summary

step_summary[, movement_to_error := median_step / sigma_error_med]
step_summary

#Movement becomes reliably distinguisable from erro at 10 seconds or over
#Relaibility stabilises at around 15-20 seconds
#No meaningful improvement beyond 20-30s

#-----------------------------------------------------------------------------------------------#

#Based on analysis on above, I have chosen to thin the dataset by 15s
#Before thinning, I will remove outliers for each species based on excessive speeds

lake_BT_movebank <- 
  with(lake_BT_sub, 
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
       ))


# Convert the dataframe to a telemetry object using the ctmm package's as.telemetry function
# No projection is applied here, and the WGS84 datum is used (common geographic coordinate system)
lake_BT_tels <- as.telemetry(lake_BT_movebank, 
                                   timezone = "Europe/Stockholm",   # Set time zone
                                   timeformat = "%Y-%m-%d %H:%M:%S",# Specify timestamp format
                                   projection = NULL,               # No projection
                                   datum = "WGS84",                 # World Geodetic System 1984
                                   keep = c("Species", "Weight","Total_length", "Std_length", "Treatment", 
                                            "Date", "Exp_Stage", "Time_Of_Day", "found_alive", "known_predated"))  # Retain extra columns




# Incoporate location error

names(lake_BT_tels)

#Center the projection on the geometric median of the data
ctmm::projection(lake_BT_tels) <- ctmm::median(lake_BT_tels)

### INCORPORATING LOCATION ERROR
# fit error parameters to calibration data
lake_BT_UERE <- uere.fit(lake_BT_tels$FReference)
summary(lake_BT_UERE)
#        low      est      high
#all 0.37805 0.379299 0.3805479

# apply error model to data
uere(lake_BT_tels) <- lake_BT_UERE
#new column now called VAR.xy

#remove reference list
names(lake_BT_tels)
lake_BT_tels <- lake_BT_tels[1:66]

#save UERE
saveRDS(lake_BT_UERE , paste0(save_telem_path, "BT/lake_BT_UERE.rds"))

#-----------------------------------------------------------------------------#

# Remove outliers based on species maximum swim speeds ####

#Pike = 0.823
#Perch = 0.977
#Roach = 0.841
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES)

#first seperate telemetry object by species 
#check whether ids are in species order
lake_BT_movebank %>%
  select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

#need to order them by id number
lake_BT_tels <- lake_BT_tels[order(names(lake_BT_tels))]
names(lake_BT_tels)

#mostly in order except for the first perch
pike_lake_BT_tel <- lake_BT_tels[61:66]
perch_lake_BT_tel <- lake_BT_tels[c(1:15, 31:44, 46)]
roach_lake_BT_tel <- lake_BT_tels[c(16:30, 45, 47:60)]

#remove outliers based on speed
#PIKE
out_pike <- outlie(pike_lake_BT_tel, plot = FALSE)
sum(sapply(out_pike, function(x) sum(x$speed > 0.823)))
#11661
#Need to filter out unrealistic speeds
which_lowSp <- lapply(out_pike, function(x) x$speed <= 0.823)
#Combining the lists and removing observations for which the logical vector was false
pike_lake_BT_tel <- Map(function(x,y) x[y,], pike_lake_BT_tel,which_lowSp)
#save telemetry object
saveRDS(pike_lake_BT_tel, paste0(save_telem_path, "BT/pike_lake_BT_tel_unthinned.rds")) 

#PERCH
out_perch <- outlie(perch_lake_BT_tel, plot = FALSE)
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))
# 34023 
#Need to filter out unrealistic speeds
which_lowSp <- lapply(out_perch, function(x) x$speed <= 0.977)
#Combining the lists and removing observations for which the logical vector was false
perch_lake_BT_tel <- Map(function(x,y) x[y,], perch_lake_BT_tel,which_lowSp)
#save telemetry object
saveRDS(perch_lake_BT_tel , paste0(save_telem_path, "BT/perch_lake_BT_tel_unthinned.rds")) 

#ROACH
out_roach <- outlie(roach_lake_BT_tel, plot = FALSE)
sum(sapply(out_roach, function(x) sum(x$speed > 0.841)))
#114783 
#Need to filter out unrealistic speeds
which_lowSp <- lapply(out_roach, function(x) x$speed <= 0.841)
#Combining the lists and removing observations for which the logical vector was false
roach_lake_BT_tel <- Map(function(x,y) x[y,], roach_lake_BT_tel,which_lowSp)
#save telemetry object
saveRDS(roach_lake_BT_tel , paste0(save_telem_path, "BT/roach_lake_BT_tel_unthinned.rds"))

#------------------------------------------------------------------------------------------------------#

#Make combined dataframe with outliers for each species removed ####

# Create a function to extract the data and add the individual_id column
# This function merges telemetry data from different individuals into a single dataframe.
extract_data <- function(list_data) {
  data_combined <- do.call(rbind, lapply(names(list_data), function(id) {
    df <- list_data[[id]]
    df$individual_ID <- id  # Add individual ID to each dataframe.
    return(df)
  }))
  return(data_combined)
}

# Apply the function to each species' telemetry data and combine them into a single dataframe.
pike_data <- extract_data(pike_lake_BT_tel)
perch_data <- extract_data(perch_lake_BT_tel)
roach_data <- extract_data(roach_lake_BT_tel)

# Combine all dataframes into one, with a species label added.
BT_telem_data <- 
  rbind(
    cbind(pike_data, Species = 'Pike'), 
    cbind(perch_data, Species = 'Perch'),
    cbind(roach_data, Species = 'Roach')
  )
#Original: 20414167
#New: 20165124
#Difference: 249043 (which checks out with removal of reference tag data and outliers)

saveRDS(BT_telem_data, paste0(filtered_data_path, "02_lake_BT_sub.rds"))

#------------------------------------------------------------------------------#

#remove unrequired large files
rm(
  lake_BT_movebank, lake_BT_sub, lake_BT_sub_dt, lake_BT_tels,
  roach_data, perch_data, pike_data,
  roach_lake_BT_tel, pike_lake_BT_tel, perch_lake_BT_tel,
  out_perch, out_pike, out_roach, subset_ids, which_lowSp
)

#free up memory
gc()

#Now I want to thin the data based on  15 second intervals
#I need to first convert dataframe into a move object
lake_BT_mv <- mt_as_move2(BT_telem_data, 
                          coords = c("longitude", "latitude"),  # specify coordinate columns
                          crs = "WGS84",  # use the WGS84 coordinate reference system
                          time_column = "timestamp",  # specify the timestamp column
                          track_id_column = "individual_ID",  # column identifying individual tracks
                          na.fail = FALSE)  # allows rows with missing coordinates

# Sort the data by individual ID and timestamp for chronological consistency
lake_BT_mv <- lake_BT_mv %>%
  dplyr::arrange(individual_ID, timestamp)

lake_BT_thin_data <- lake_BT_mv %>% mt_filter_per_interval(unit = "10 seconds", criterion = "first")
#OLD: 20165124
#NEW: 7400067

lake_BT_thin_track_sum <-
  lake_BT_thin_data %>%
  group_by(individual_ID) %>%
  arrange(timestamp) %>%
  mutate(dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
  summarise(
    min_dt = min(dt, na.rm = TRUE),
    median_dt = median(dt, na.rm = TRUE),
    max_dt = max(dt, na.rm = TRUE)
  )

#find when the largest time gaps occur for each individual
# If not already a data.table:
lake_BT_thin_data <- as.data.table(lake_BT_thin_data)

# Sort by individual and time (important!)
setorder(lake_BT_thin_data, individual_ID, timestamp)

# Time gap in seconds between successive fixes for each individual
lake_BT_thin_data[, dt_mins := as.numeric(
  difftime(timestamp, shift(timestamp), units = "mins")
), by = individual_ID]

# Also store the previous timestamp and its date (start of the gap)
lake_BT_thin_data[, prev_time := shift(timestamp), by = individual_ID]
lake_BT_thin_data[, prev_date := shift(Date), by = individual_ID]

max_gaps <- lake_BT_thin_data[!is.na(dt_sec),
                              .SD[which.max(dt_sec)],
                              by = individual_ID
][
  ,
  .(
    individual_ID,
    gap_start_time = prev_time,
    gap_start_date = prev_date,
    gap_end_time   = timestamp,
    gap_length_mins = dt_sec
  )
]

max_gaps

# Extract the longitude and latitude columns from the geometry column
coords <- st_coordinates(lake_BT_thin_data$geometry)
lake_BT_thin_data$Long <- coords[, 1]
lake_BT_thin_data$Lat <- coords[, 2]

#------------------------------------------------------------------------------------------#

#Explore daily movement trajectories
#Remove any individuals that clearly that died early in the experiment
#Identify any problematic individuals that might need some further exploration

#Function that takes data for a single individual and returns a ggplot
plot_daily_traj_individual <- function(dat, lake_poly) {
  
  this_id <- unique(dat$individual_ID)
  
  ggplot() +
    # lake polygon
    geom_sf(data = lake_poly, fill = "lightblue", color = "black", alpha = 0.3) +
    
    # trajectories (grouped by date so each day's path is separate)
    geom_path(
      data = dat,
      aes(x = Long, y = Lat, group = Date),
      linewidth = 0.3
    ) +
    
    # points
    geom_point(
      data = dat,
      aes(x = Long, y = Lat),
      size = 0.4
    ) +
    
    facet_wrap(~ Date) +
    coord_sf() +
    theme_minimal() +
    labs(
      title = paste("Daily Movement Trajectories of", this_id, "in BT Lake"),
      x = "Longitude",
      y = "Latitude"
    )
}

#Split data by individual and build a named list of plots
traj_plots <- lake_BT_sub %>%
  group_by(individual_ID) %>%
  group_split() %>%                                   # list of data frames, one per individual
  set_names(map_chr(., ~ unique(.x$individual_ID))) %>%  # name each list element by the ID
  map(~ plot_daily_traj_individual(.x, BT_polygon))  


walk(names(traj_plots), function(id) {
  
  filename <- paste0("./daily_trajectory_plots/lake_BT/", id, "_daily_trajectory_plot.png")
  
  ggsave(
    filename = filename,
    plot = traj_plots[[id]],
    width = 10,
    height = 8,
    dpi = 300
  )
  
  message("Saved: ", filename)
})

#------------------------------------------------------------------------------#

#Filter out dead pike days

#Filter out individual F59889 because it died on the first day of tracking
lake_BT_sub <- lake_BT_sub %>% 
  filter(individual_ID != "F59889")

#load information about pike estimated death days
pike_deaths <- read.csv("./data/pike_deaths.csv")

pike_mort_cols <- 
  pike_deaths %>%
  filter(individual_ID %in% lake_BT_sub$individual_ID) %>% 
  mutate(
    pike_death_date = as.Date(likely_death_date, format = "%d/%m/%Y"))

# Filter out data after the predation or mortality event for each individual.
lake_BT_sub <- 
  lake_BT_sub %>%
  left_join(pike_mort_cols, by = c("individual_ID" = "individual_ID")) %>%
  filter(is.na(pike_death_date)| Date < pike_death_date)  # Keep locations only before the death date if one is recorded
#pre-filter rows: 7252120
#post-filter rows: 7204217

print(paste0("Rows removed after  filtering: ", nrow(lake_BT_sub) - nrow(test)))
#47903

saveRDS(lake_BT_sub, paste0(filtered_data_path, "03_lake_BT_sub.rds"))

