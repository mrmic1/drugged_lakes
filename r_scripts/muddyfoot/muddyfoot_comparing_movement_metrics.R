#---------------------------------------#
# Comparing movement metrics - Muddyfoot  
#---------------------------------------#

#LIBRARIES
library(ctmm)
library(tidyverse)
library(move2)
library(tidybayes)
library(nlme)
library(lme4)
library(performance)
library(modelbased)
library(emmeans)

# Set the time zone environment to 'Europe/Stockholm' for consistent timestamp manipulation
Sys.setenv(TZ = 'Europe/Stockholm')

### DIRECTORIES ###
# Define paths to directories for loading/saving data
ctmm_path <- "./data/ctmm_fits/"                  # Directory for ctmm model fits
data_filter_path <- "./data/tracks_filtered/"     # Directory for filtered telemetry data
telem_path <- "./data/telem_obj/"                 # Directory for telemetry objects
akde_path <- "./data/akdes/"                      # Directory for AKDE (Autocorrelated Kernel Density Estimation) outputs
lake_polygon_path <- "./data/lake_coords/"        # Directory for lake polygon (boundary) data
enc_path <- "./data/encounters/"                  # Directory for encounter data


### LOADING DATA ###

# Load pre-filter muddyfoot tracking data
muddyfoot_track_data <-  readRDS(paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))

# Load telemetry objects for pike, perch, and roach species in Muddyfoot lake
pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

# Load ctmm model fits for pike, perch, and roach in Muddyfoot lake
# The models include continuous-time movement fits using the Ornstein-Uhlenbeck Foraging (OUF) model
pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds"))

#Load akdes for pike, perch, and roach in muddyfoot lake
pike_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/pike_akdes_cg_list.rds"))
perch_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))
roach_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))

# Load the polygon representing the boundary of Muddyfoot lake, in the form of a GeoPackage file
muddyfoot_polygon <- sf::st_read(paste0(lake_polygon_path, "muddyfoot/lake_muddyfoot_polygon.gpkg"))

#---------------------------------------------------------------------------------------------------------------#

# Remove predated individuals ####

# Load predation event data (pre-identified predation events)
mud_pred_events <- readRDS(paste0(enc_path, "muddyfoot_pred_events.rds"))

# Remove Roach individuals that were predated based on predation events
roach_ids_remove <- mud_pred_events %>%
  filter(Species == "Roach") %>%    # Filter for Roach species
  pull(individual_ID)               # Extract IDs of predated Roach

# Update the Roach AKDE list and telemetry data, removing predated individuals
roach_akdes_cg_list <- roach_akdes_cg_list[!(names(roach_akdes_cg_list) %in% roach_ids_remove)]
roach_muddyfoot_tel <- roach_muddyfoot_tel[!(names(roach_muddyfoot_tel) %in% roach_ids_remove)]

# Remove Perch individuals that were predated based on predation events
perch_ids_remove <- mud_pred_events %>%
  filter(Species == "Perch") %>%    # Filter for Perch species
  pull(individual_ID)               # Extract IDs of predated Perch

# Update the Perch AKDE list and telemetry data, removing predated individuals
perch_akdes_cg_list <- perch_akdes_cg_list[!(names(perch_akdes_cg_list) %in% perch_ids_remove)]
perch_muddyfoot_tel <- perch_muddyfoot_tel[!(names(perch_muddyfoot_tel) %in% perch_ids_remove)]

#----------------------------------------------------------------------------------------------------------------#

#----------------------------#
#> 1. Trajectory analysis ####
#---------------------------#

#create a move2 object from a data.frame
muddyfoot_mv <- mt_as_move2(muddyfoot_track_data,
                            coords = c("longitude","latitude"),
                            crs = "WGS84",
                            time_column = "timestamp",
                            track_id_column = "individual_id",
                            na.fail = F) # allows or not empty coordinates



# Annotate speed, azimuth and turning angle to the trajectory.
muddyfoot_mv <- muddyfoot_mv %>% 
  mutate(
    azimuth = mt_azimuth(muddyfoot_mv), 
    speed = mt_speed(muddyfoot_mv), 
    turnangle = mt_turnangle(muddyfoot_mv),
    distance = mt_distance(muddyfoot_mv))


# Change units
muddyfoot_mv$speed <- units::set_units(muddyfoot_mv$speed, cm/s)
muddyfoot_mv$distance <- units::set_units(muddyfoot_mv$distance, cm)
str(muddyfoot_mv$speed)

#convert back to dataframe
muddyfoot_track_data <- as.data.frame(muddyfoot_mv)

#Get the date

muddyfoot_track_data$new_date <- strftime(muddyfoot_track_data$timestamp, format="%Y/%m/%d")

tz(muddyfoot_track_data$new_date)
tz(muddyfoot_track_data$timestamp)
tz(muddyfoot_track_data$Date)

#Summarise per day
muddyfoot_daily_tracks <- 
  muddyfoot_track_data %>%
  dplyr::group_by(individual_id, Date) %>%
  dplyr::mutate(
    avg_daily_speed = mean(speed),
    total_distance = sum(distance)) %>% 
  dplyr::distinct(individual_id, Date, .keep_all = TRUE) %>%
  ungroup()

#-----------------------------------------------------------------------------------------#


analyze_movement <- function(species_name, response_var, muddyfoot_daily_tracks) {
  
  # Filter the dataset for the given species
  species_move_daily <- muddyfoot_daily_tracks %>% filter(Species == species_name)
  
  # Create a factor for date
  species_move_daily$Date_fact <- as.factor(species_move_daily$Date)
  
  # Plot for the given response variable
  ggplot(species_move_daily, 
         aes(x = factor(treatment), 
             y = !!sym(response_var), 
             fill = factor(treatment))) +
    stat_halfeye(adjust = 0.75, justification = -0.1, .width = 0, point_colour = NA) +
    geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.5) +
    theme_bw() +
    labs(x = "Treatment", y = paste("Mean", gsub("_", " ", response_var))) +
    coord_flip() +
    theme(legend.position = "none")
  
  # Linear mixed model
  model <- lmer(as.formula(paste(response_var, "~ treatment + (1|Date_fact) + (1|individual_id)")), 
                data = species_move_daily)
  
  
  # Estimated means and contrasts
  means <- estimate_means(model)
  contrasts <- estimate_contrasts(model, contrast = "treatment", p_adjust = "tukey")
  
  return(list(means = means, contrasts = contrasts))
}


# Run for each species and response variables
species_list <- c("Perch", "Roach", "Pike")
response_vars <- c("avg_daily_speed", "total_distance")

results <- list()

for (species in species_list) {
  for (response in response_vars) {
    results[[paste(species, response, sep = "_")]] <- analyze_movement(species, response, muddyfoot_daily_tracks)
  }
}

results$Perch_total_distance










#> 2.1. Perch ####

perch_move_daily <- 
  muddyfoot_daily_tracks %>% 
  filter(Species == 'Perch')

#Plot
ggplot(perch_move_daily, 
       aes(x = factor(treatment), 
           y = avg_daily_speed, 
           fill = factor(treatment))) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.75,
    # move to the right
    justification = -0.1,
    # remove the slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "Mean daily speed",
  ) +
  coord_flip()+
  theme(legend.position="none")

#Linear mixed model
head(perch_move_daily)
perch_move_daily$Date_fact <- as.factor(perch_move_daily$Date)

perch_speed <- lmer(avg_daily_speed ~ treatment + (1|Date_fact) + (1|individual_id), 
                    data = perch_move_daily)

performance::check_model(perch_speed)

summary(perch_speed)

# 2. Obtain estimated means
perch_speed_means <- estimate_means(perch_speed)


# 3. Obtain estimated contrasts
perch_speed_contrast <- estimate_contrasts(perch_speed, contrast = "treatment", p_adjust = "tukey")

perch_speed_means
perch_speed_contrast




















#----------------------------------#
# > 1. AKDE Meta-analysis ##########
#----------------------------------#


#> 1.2. Pike AKDE treatment comparison ####

meta(pike_akdes_cg_list,col='black',sort=F, verbose = T, level.UD = 0.95)

# plot AKDEs with consistent grids
COL <- color(pike_akdes_cg_list, by='individual')
plot(pike_akdes_cg_list,
     col.UD=COL_cg,
     col.level=COL_cg,
     col.grid=NA,
     level=NA,
     main="Muddyfoot pike AKDEs")

#Seperating treatment groups 
pike_akde_control <- pike_akdes_cg_list[1:3]
pike_akde_mix <- pike_akdes_cg_list[4:6]

#Combine back into list based on treatment
pike_akde_total <- list(pike_akde_control = pike_akde_control, 
                        pike_akde_mix = pike_akde_mix)

#Estimate population-level mean parameters
meta(pike_akde_total,col='black',sort=F, verbose = T, level.UD = 0.95)




OVER <- overlap(pike_akdes_cg_list)



#Estimating the CDE
CDE <- encounter(pike_akdes_list[c("F59880", "F59884")])
#inconsistent grid resolution
CDE_cg <- encounter(pike_akdes_cg_list[c("F59880", "F59884")])
summary(CDE_cg)
CDE_cg

#Plot the data and HR estimates
ctmm::plot(pike_muddyfoot_tel[c("F59880", "F59884")],
     UD = pike_akdes_cg_list[c("F59880", "F59884")],
     col = NA,
     col.DF= c("#f4a261", "#2a9d8f"),
     col.grid = NA)
#this works


#Visualise the CDE
ctmm::plot(pike_muddyfoot_tel[c("F59880", "F59884")],
           UD = CDE_cg,
           col.UD="red",
           #col=c("#e76f51", "#264653"),
           error = FALSE,
           col.grid = NA)
#Error: unable to find an inherited method for function ‘extent’ for signature ‘x = "overlap"’
#the same error message pops up when running the same code in ctmm_interactions tutorial.


