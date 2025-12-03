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

#-----------------------------------------------------------------------------#

#Load in the datasets
lake_BT_sub <- readRDS(paste0(filtered_data_path, '03_lake_BT_sub.rds'))

#Need to create telemetry objects again
lake_BT_movebank <- 
  with(lake_BT_sub, 
       data.frame(
         "timestamp" = timestamp,                        
         "location.long" = Long,                         
         "location.lat" = Lat, 
         "GPS.HDOP" = HDOP,                               
         "individual-local-identifier" = individual_ID,
         "Species" = Species,
         "Weight" = Weight,
         "Total_length" = Total_length,
         "Std_length" = Std_length,
         "Treatment" = Treatment,                        
         "Date" = Date,
         "Exp_Stage" = Exp_Stage,
         "Time_Of_Day" = Time_Of_Day,
         "found_alive" = found_alive,
         "known_predated" = known_predated
       ))

rm(lake_BT_sub)

# Convert the dataframe to a telemetry object using the ctmm package's as.telemetry function
# No projection is applied here, and the WGS84 datum is used (common geographic coordinate system)
lake_BT_tels <- as.telemetry(lake_BT_movebank, 
                             timezone = "Europe/Stockholm",   # Set time zone
                             timeformat = "%Y-%m-%d %H:%M:%S",# Specify timestamp format
                             projection = NULL,               # No projection
                             datum = "WGS84",                 # World Geodetic System 1984
                             keep = c("Species", "Weight","Total_length", "Std_length", "Treatment", 
                                      "Date", "Exp_Stage", "Time_Of_Day", "found_alive", "known_predated"))


#columns x and y are very different from 1. I think it has something to do with the projection
ctmm::projection(lake_BT_tels$F59749)
#tpeqd projection
tz(lake_BT_tels$F59749$timestamp)
#"Europe/Stockholm"

#Split telemetry object by species

lake_BT_movebank %>%
  select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

#mostly in order except for the first perch
pike_lake_BT_tel <- lake_BT_tels[61:65]
perch_lake_BT_tel <- lake_BT_tels[c(1:15, 31:44, 46)]
roach_lake_BT_tel <- lake_BT_tels[c(16:30, 45, 47:60)]

#save telemetry objects
saveRDS(pike_lake_BT_tel, paste0(save_telem_path, "BT/pike_lake_BT_tel_thinned.rds"))
saveRDS(perch_lake_BT_tel, paste0(save_telem_path, "BT/perch_lake_BT_tel_thinned.rds"))
saveRDS(roach_lake_BT_tel, paste0(save_telem_path, "BT/roach_lake_BT_tel_thinned.rds"))

#---------------------------------------------------------------------------------#

#### RUN CTMMS ####

#load telemetry objects if needed
pike_lake_BT_tel <- readRDS(paste0(save_telem_path, "BT/pike_lake_BT_tel_thinned.rds"))
perch_lake_BT_tel <- readRDS(paste0(save_telem_path, "BT/perch_lake_BT_tel_thinned.rds"))
roach_lake_BT_tel <- readRDS(paste0(save_telem_path, "BT/roach_lake_BT_tel_thinned.rds"))

### >>> Pike ctmms ####
#Run ctmm models for five pike
names(pike_lake_BT_tel)

#look at what ctmm guess predicts the best model for each individual will be 
for (i in seq_along(pike_lake_BT_tel)) {
  
  tel_i <- pike_lake_BT_tel[[i]]
  id_i  <- names(pike_lake_BT_tel)[i]
  
  # Run ctmm.guess
  guess_i <- ctmm.guess(
    tel_i,
    CTMM        = ctmm(error = TRUE),
    interactive = FALSE
  )
  
  # Determine model type by tau length
  tau_len <- length(guess_i$tau)
  
  model_type <- dplyr::case_when(
    tau_len == 0 ~ "IID (uncorrelated)",
    tau_len == 1 ~ "OU (range resident, no velocity autocorrelation)",
    tau_len == 2 ~ "OUF (directional persistence + range residency)",
    TRUE         ~ "Unknown"
  )
  
  # Print summary
  message("ID: ", id_i,
          "   |   tau length = ", tau_len,
          "   |   guessed model = ", model_type)
}

# ID: F59886   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59887   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59888   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59891   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59892   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)

#Prepare list for ctmm fit
lake_BT_pike_ctmm_fits <- vector("list", length(pike_lake_BT_tel))
names(lake_BT_pike_ctmm_fits) <- names(pike_lake_BT_tel)

# Loop over individuals
for (i in seq_along(pike_lake_BT_tel)) {
  
  # Current individual's telemetry data
  tel_i   <- pike_lake_BT_tel[[i]]
  id_i    <- names(pike_lake_BT_tel)[i]
  
  message("Fitting ctmm for pike: ", id_i, " (", i, "/", length(pike_lake_BT_tel), ")")
  
  # 1. Get an initial guess model (error-informed)
  lake_BT_pike_guess <- ctmm.guess(
    tel_i,
    CTMM       = ctmm(error = TRUE),
    interactive = FALSE
  )
  
  # 2. Fit the model using maximum likelihood
  model_fit <- ctmm.fit(
    data   = tel_i,
    CTMM   = lake_BT_pike_guess,
    method = "ML"
  )
  
  # 3. Save the fitted model to disk
  saveRDS(
    model_fit,
    file = file.path(
      save_ctmm_path,
      "lake_BT_pike_fits",
      paste0(id_i, "_ctmm_fit.rds")
    )
  )
  
  # 4. Store in list
  lake_BT_pike_ctmm_fits[[i]] <- model_fit
}

#add ID to lists
names(lake_BT_pike_ctmm_fits ) <- names(pike_lake_BT_tel)
summary(lake_BT_pike_ctmm_fits[[1]])
#check fits
plot(pike_lake_BT_tel[[1]], lake_BT_pike_ctmm_fits[[1]], error = FALSE)

#save fit list
saveRDS(lake_BT_pike_ctmm_fits, paste0(save_ctmm_path, "lake_BT_pike_fits/lake_BT_pike_ctmm_fits"))

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

### >>> Perch ctmms ####
#Run ctmm models for 30 perch
names(perch_lake_BT_tel)

#look at what ctmm guess predicts the best model for each individual will be 
for (i in seq_along(perch_lake_BT_tel)) {
  
  tel_i <- perch_lake_BT_tel[[i]]
  id_i  <- names(perch_lake_BT_tel)[i]
  
  # Run ctmm.guess
  guess_i <- ctmm.guess(
    tel_i,
    CTMM        = ctmm(error = TRUE),
    interactive = FALSE
  )
  
  # Determine model type by tau length
  tau_len <- length(guess_i$tau)
  
  model_type <- dplyr::case_when(
    tau_len == 0 ~ "IID (uncorrelated)",
    tau_len == 1 ~ "OU (range resident, no velocity autocorrelation)",
    tau_len == 2 ~ "OUF (directional persistence + range residency)",
    TRUE         ~ "Unknown"
  )
  
  # Print summary
  message("ID: ", id_i,
          "   |   tau length = ", tau_len,
          "   |   guessed model = ", model_type)
}

# ID: F59749   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59750   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59751   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59752   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59753   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59754   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59755   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59756   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59757   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59758   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59760   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59762   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59763   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59764   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59765   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59784   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59785   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59786   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59787   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59789   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59790   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59791   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59792   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59793   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59794   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59795   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59796   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59797   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59799   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59801   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)


#Prepare list for ctmm fits
lake_BT_perch_ctmm_fits <- vector("list", length(perch_lake_BT_tel))
names(lake_BT_perch_ctmm_fits) <- names(perch_lake_BT_tel)


# Loop over individuals
for (i in seq_along(perch_lake_BT_tel)) {
  
  # Current individual's telemetry data
  tel_i   <- perch_lake_BT_tel[[i]]
  id_i    <- names(perch_lake_BT_tel)[i]
  
  message("Fitting ctmm for perch: ", id_i, " (", i, "/", length(perch_lake_BT_tel), ")")
  
  # 1. Get an initial guess model (error-informed)
  lake_BT_perch_guess <- ctmm.guess(
    tel_i,
    CTMM       = ctmm(error = TRUE),
    interactive = FALSE
  )
  
  # 2. Fit the model using maximum likelihood
  model_fit <- ctmm.fit(
    data   = tel_i,
    CTMM   = lake_BT_perch_guess,
    method = "ML"
  )
  
  # 3. Save the fitted model to disk
  saveRDS(
    model_fit,
    file = file.path(
      save_ctmm_path,
      "lake_BT_perch_fits",
      paste0(id_i, "_ctmm_fit.rds")
    )
  )
  
  # 4. Store in list
  lake_BT_perch_ctmm_fits[[i]] <- model_fit
}

#add ID to lists
names(lake_BT_perch_ctmm_fits ) <- names(perch_lake_BT_tel)
summary(lake_BT_perch_ctmm_fits[[1]])
#check fits
plot(perch_lake_BT_tel[[1]], lake_BT_perch_ctmm_fits[[1]], error = FALSE)

#save fit list
saveRDS(lake_BT_perch_ctmm_fits, paste0(save_ctmm_path, "lake_BT_perch_fits/lake_BT_perch_ctmm_fits"))

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

### >>> Roach ctmms ####

#Run ctmm models for 30 roach
names(roach_lake_BT_tel)

#look at what ctmm guess predicts the best model for each individual will be 
for (i in seq_along(roach_lake_BT_tel)) {
  
  tel_i <- roach_lake_BT_tel[[i]]
  id_i  <- names(roach_lake_BT_tel)[i]
  
  # Run ctmm.guess
  guess_i <- ctmm.guess(
    tel_i,
    CTMM        = ctmm(error = TRUE),
    interactive = FALSE
  )
  
  # Determine model type by tau length
  tau_len <- length(guess_i$tau)
  
  model_type <- dplyr::case_when(
    tau_len == 0 ~ "IID (uncorrelated)",
    tau_len == 1 ~ "OU (range resident, no velocity autocorrelation)",
    tau_len == 2 ~ "OUF (directional persistence + range residency)",
    TRUE         ~ "Unknown"
  )
  
  # Print summary
  message("ID: ", id_i,
          "   |   tau length = ", tau_len,
          "   |   guessed model = ", model_type)
}


# ID: F59766   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59767   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59768   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59770   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59772   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59773   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59774   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59775   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59776   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59777   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59779   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59780   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59781   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59782   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59783   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59800   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59802   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59803   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59804   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59805   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59807   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59808   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59809   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59810   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59811   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59812   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59813   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59814   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59815   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)
# ID: F59816   |   tau length = 2   |   guessed model = OUF (directional persistence + range residency)

#Prepare list for ctmm fits
lake_BT_roach_ctmm_fits <- vector("list", length(roach_lake_BT_tel))
names(lake_BT_roach_ctmm_fits) <- names(roach_lake_BT_tel)


# Loop over individuals
for (i in seq_along(roach_lake_BT_tel)) {
  
  # Current individual's telemetry data
  tel_i   <- roach_lake_BT_tel[[i]]
  id_i    <- names(roach_lake_BT_tel)[i]
  
  message("Fitting ctmm for roach: ", id_i, " (", i, "/", length(roach_lake_BT_tel), ")")
  
  # 1. Get an initial guess model (error-informed)
  lake_BT_roach_guess <- ctmm.guess(
    tel_i,
    CTMM       = ctmm(error = TRUE),
    interactive = FALSE
  )
  
  # 2. Fit the model using maximum likelihood
  model_fit <- ctmm.fit(
    data   = tel_i,
    CTMM   = lake_BT_roach_guess,
    method = "ML"
  )
  
  # 3. Save the fitted model to disk
  saveRDS(
    model_fit,
    file = file.path(
      save_ctmm_path,
      "lake_BT_roach_fits",
      paste0(id_i, "_ctmm_fit.rds")
    )
  )
  
  # 4. Store in list
  lake_BT_roach_ctmm_fits[[i]] <- model_fit
}

#add ID to lists
names(lake_BT_roach_ctmm_fits ) <- names(roach_lake_BT_tel)
summary(lake_BT_roach_ctmm_fits[[1]])
#check fits
plot(roach_lake_BT_tel[[1]], lake_BT_roach_ctmm_fits[[1]], error = FALSE)

#save fit list
saveRDS(lake_BT_roach_ctmm_fits, paste0(save_ctmm_path, "lake_BT_roach_fits/lake_BT_roach_ctmm_fits"))