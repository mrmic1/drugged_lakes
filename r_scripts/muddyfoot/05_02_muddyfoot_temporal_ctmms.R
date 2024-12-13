#------------------------------------------------------#
# Create temporal dataframes and CTMM objects #
#------------------------------------------------------#

#The goal of this script is to create two dataframes that split the data
#into two timepoints. These will be used to compare behavioral differences early and later 
#in the experiment. I will create individual movement models for each timepoint for each individual

###  LIBRARIES ###

library(data.table)
library(tidyverse)
library(ctmm)
library(sf)
#for parallel processing
library(parallel)
library(foreach)
library(doParallel)


# Paths to directories containing the necessary data files
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"


#Load data where data has been filtered rows outside of a species critical speed limit
#Post-predation tracks have also been filtered
mud_track_data <- readRDS(paste0(filtered_data_path, "04_muddyfoot_sub.rds"))

#I need to pull the Reference tag data again 
#load dataframe that contains this information
get_ref_data <-  readRDS(paste0(filtered_data_path, "01_muddyfoot_sub.rds"))

ref_data <- get_ref_data %>% 
  filter(individual_ID == 'FReference') %>% 
  dplyr::select(individual_ID, HPE, timestamp_cest, Long, Lat)


ref_data <- with(ref_data, 
                 data.frame(
                   "timestamp" = timestamp_cest,                        
                   "location.long" = Long,                         
                   "location.lat" = Lat, 
                   "GPS.HDOP" = HPE,                               
                   "individual-local-identifier" = individual_ID
                 ))

#create telemetry object for the reference data
ref_tel <- as.telemetry(ref_data, 
                        timezone = "Europe/Stockholm",   
                        timeformat = "%Y-%m-%d %H:%M:%S",
                        projection = NULL,               
                        datum = "WGS84")

ctmm::projection(ref_tel) <- ctmm::median(ref_tel)
UERE <- uere.fit(ref_tel)
summary(UERE)

# > 1. PIKE ####

#make seperate telemetry object for pike early and late experiment 

pike_data <- mud_track_data %>% 
  filter(Species == 'Pike')

pike_data <- with(pike_data, 
                  data.frame(
                    "timestamp" = timestamp,                       
                    "location.long" = longitude,                         
                    "location.lat" = latitude,                           
                    "GPS.HDOP" = HDOP,                               
                    "individual-local-identifier" = individual_ID,        
                    "Treatment" = Treatment,                        
                    "Date" = Date,                                  
                    "Exp_Stage" = Exp_Stage,
                    "Time_Of_Day" = Time_Of_Day))


#seperate into early and late to make seperate telemetry objects
pike_data_early <- pike_data %>% 
  filter(Exp_Stage == 'Early')

pike_data_late <- pike_data %>% 
  filter(Exp_Stage == 'Late')

#create telemetry objects
pike_early_tel <- as.telemetry(pike_data_early, 
                               timezone = "Europe/Stockholm",   
                               timeformat = "%Y-%m-%d %H:%M:%S", 
                               projection = NULL,               
                               datum = "WGS84",                 
                               keep = c("Treatment", "Date", "Exp_Stage", "Time_Of_Day"))  


pike_late_tel <- as.telemetry(pike_data_late, 
                              timezone = "Europe/Stockholm",   
                              timeformat = "%Y-%m-%d %H:%M:%S", 
                              projection = NULL,               
                              datum = "WGS84",                 
                              keep = c("Treatment", "Date", "Exp_Stage", "Time_Of_Day")  
)  


# Centers the data on the geometric median of the observed locations
ctmm::projection(pike_early_tel) <- ctmm::median(pike_early_tel)
ctmm::projection(pike_late_tel) <- ctmm::median(pike_late_tel)

#Add error model
uere(pike_early_tel) <- UERE
uere(pike_late_tel) <- UERE

#save telemetry objects
saveRDS(pike_early_tel, paste0(save_telem_path, "muddyfoot/stage/pike_mud_early_tel.rds"))
saveRDS(pike_late_tel, paste0(save_telem_path, "muddyfoot/stage/pike_mud_late_tel.rds"))


### CTMM Model Fitting ###

#load telemetry objects
pike_early_tel <- readRDS(paste0(save_telem_path, "muddyfoot/stage/pike_mud_early_tel.rds"))
pike_late_tel <- readRDS(paste0(save_telem_path, "muddyfoot/stage/pike_mud_late_tel.rds"))

### EARLY ###
cl <- makeCluster(6)  
doParallel::registerDoParallel(cl)

# Fit ctmm models in parallel for each pike
pike_mud_early_ctmm_fits <- 
  foreach(i = 1:length(pike_early_tel), .packages = 'ctmm') %dopar% {
  # Generate an initial guess for the model parameters, incorporating location error
  pike_early_guess <- ctmm.guess(pike_early_tel[[i]], CTMM = ctmm(error = TRUE), interactive = FALSE)
  
  # Fit the ctmm model to the telemetry data using model selection
  model_fit <- ctmm.fit(pike_early_tel[[i]], pike_early_guess, method = 'ML')
  
  # Save the fitted model for each fish individually
  saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_pike_fits/stage/", names(pike_early_tel)[i], ".rds"))
  
  model_fit  # Return the fitted model
}

# Stop the parallel cluster once model fitting is complete
stopCluster(cl)


### LATE ###
cl <- makeCluster(6)  
doParallel::registerDoParallel(cl)

pike_mud_late_ctmm_fits <- 
  foreach(i = 1:length(pike_late_tel), .packages = 'ctmm') %dopar% {
    # Generate an initial guess for the model parameters, incorporating location error
    pike_late_guess <- ctmm.guess(pike_late_tel[[i]], CTMM = ctmm(error = TRUE), interactive = FALSE)
    
    # Fit the ctmm model to the telemetry data using model selection
    model_fit <- ctmm.fit(pike_late_tel[[i]], pike_late_guess, method = 'ML')
    
    # Save the fitted model for each fish individually
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_pike_fits/stage/", names(pike_late_tel)[i], ".rds"))
    
    model_fit  # Return the fitted model
  }

# Stop the parallel cluster once model fitting is complete
stopCluster(cl)


#----------------------------------------------------------------------------------------#

# > 2. PERCH ####

#make seperate telemetry object for perch early and late experiment 

perch_data <- mud_track_data %>% 
  filter(Species == 'Perch')

perch_data <- with(perch_data, 
                  data.frame(
                    "timestamp" = timestamp,                       
                    "location.long" = longitude,                         
                    "location.lat" = latitude,                           
                    "GPS.HDOP" = HDOP,                               
                    "individual-local-identifier" = individual_ID,        
                    "Treatment" = Treatment,                        
                    "Date" = Date,                                  
                    "Exp_Stage" = Exp_Stage,
                    "Time_Of_Day" = Time_Of_Day))


#seperate into early and late to make seperate telemetry objects
perch_data_early <- perch_data %>% 
  filter(Exp_Stage == 'Early')

perch_data_late <- perch_data %>% 
  filter(Exp_Stage == 'Late')

#create telemetry objects
perch_early_tel <- as.telemetry(perch_data_early, 
                               timezone = "Europe/Stockholm",   
                               timeformat = "%Y-%m-%d %H:%M:%S", 
                               projection = NULL,               
                               datum = "WGS84",                 
                               keep = c("Treatment", "Date", "Exp_Stage", "Time_Of_Day"))  


perch_late_tel <- as.telemetry(perch_data_late, 
                              timezone = "Europe/Stockholm",   
                              timeformat = "%Y-%m-%d %H:%M:%S", 
                              projection = NULL,               
                              datum = "WGS84",                 
                              keep = c("Treatment", "Date", "Exp_Stage", "Time_Of_Day")  
)  


# Centers the data on the geometric median of the observed locations
ctmm::projection(perch_early_tel) <- ctmm::median(perch_early_tel)
ctmm::projection(perch_late_tel) <- ctmm::median(perch_late_tel)

#Add error model
uere(perch_early_tel) <- UERE
uere(perch_late_tel) <- UERE

#save telemetry objects
saveRDS(perch_early_tel, paste0(save_telem_path, "muddyfoot/stage/perch_mud_early_tel.rds"))
saveRDS(perch_late_tel, paste0(save_telem_path, "muddyfoot/stage/perch_mud_late_tel.rds"))


### CTMM Model Fitting ###
#load telemetry models
perch_early_tel <- readRDS(paste0(save_telem_path, "muddyfoot/stage/perch_mud_early_tel.rds"))
perch_late_tel <- readRDS(paste0(save_telem_path, "muddyfoot/stage/perch_mud_late_tel.rds"))


### EARLY ###
cl <- makeCluster(15)  
doParallel::registerDoParallel(cl)

# Fit ctmm models in parallel for each perch
perch_mud_early_ctmm_fits <- 
  foreach(i = 1:length(perch_early_tel), .packages = 'ctmm') %dopar% {
    # Generate an initial guess for the model parameters, incorporating location error
    perch_early_guess <- ctmm.guess(perch_early_tel[[i]], CTMM = ctmm(error = TRUE), interactive = FALSE)
    
    # Fit the ctmm model to the telemetry data using model selection
    model_fit <- ctmm.fit(perch_early_tel[[i]], perch_early_guess, method = 'ML')
    
    # Save the fitted model for each fish individually
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_perch_fits/stage/", names(perch_early_tel)[i], ".rds"))
    
    model_fit  # Return the fitted model
  }

# Stop the parallel cluster once model fitting is complete
stopCluster(cl)


### LATE ###
cl <- makeCluster(15)  
doParallel::registerDoParallel(cl)

perch_mud_late_ctmm_fits <- 
  foreach(i = 1:length(perch_late_tel), .packages = 'ctmm') %dopar% {
    # Generate an initial guess for the model parameters, incorporating location error
    perch_late_guess <- ctmm.guess(perch_late_tel[[i]], CTMM = ctmm(error = TRUE), interactive = FALSE)
    
    # Fit the ctmm model to the telemetry data using model selection
    model_fit <- ctmm.fit(perch_late_tel[[i]], perch_late_guess, method = 'ML')
    
    # Save the fitted model for each fish individually
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_perch_fits/stage/", names(perch_late_tel)[i], ".rds"))
    
    model_fit  # Return the fitted model
  }

# Stop the parallel cluster once model fitting is complete
stopCluster(cl)


#----------------------------------------------------------------------------------------------------#

# > 3. Roach ####

#make seperate telemetry object for roach early and late experiment 

roach_data <- mud_track_data %>% 
  filter(Species == 'Roach')

roach_data <- with(roach_data, 
                   data.frame(
                     "timestamp" = timestamp,                       
                     "location.long" = longitude,                         
                     "location.lat" = latitude,                           
                     "GPS.HDOP" = HDOP,                               
                     "individual-local-identifier" = individual_ID,        
                     "Treatment" = Treatment,                        
                     "Date" = Date,                                  
                     "Exp_Stage" = Exp_Stage,
                     "Time_Of_Day" = Time_Of_Day))


#seperate into early and late to make seperate telemetry objects
roach_data_early <- roach_data %>% 
  filter(Exp_Stage == 'Early')

roach_data_late <- roach_data %>% 
  filter(Exp_Stage == 'Late')

#create telemetry objects
roach_early_tel <- as.telemetry(roach_data_early, 
                                timezone = "Europe/Stockholm",   
                                timeformat = "%Y-%m-%d %H:%M:%S", 
                                projection = NULL,               
                                datum = "WGS84",                 
                                keep = c("Treatment", "Date", "Exp_Stage", "Time_Of_Day"))  


roach_late_tel <- as.telemetry(roach_data_late, 
                               timezone = "Europe/Stockholm",   
                               timeformat = "%Y-%m-%d %H:%M:%S", 
                               projection = NULL,               
                               datum = "WGS84",                 
                               keep = c("Treatment", "Date", "Exp_Stage", "Time_Of_Day")  
)  

# Centers the data on the geometric median of the observed locations
ctmm::projection(roach_early_tel) <- ctmm::median(roach_early_tel)
ctmm::projection(roach_late_tel) <- ctmm::median(roach_late_tel)

#Add error model
uere(roach_early_tel) <- UERE
uere(roach_late_tel) <- UERE

#save telemetry objects
saveRDS(roach_early_tel, paste0(save_telem_path, "muddyfoot/stage/roach_mud_early_tel.rds"))
saveRDS(roach_late_tel, paste0(save_telem_path, "muddyfoot/stage/roach_mud_late_tel.rds"))

### CTMM Model Fitting ###

#load telemetry models
roach_early_tel <- readRDS(paste0(save_telem_path, "muddyfoot/stage/roach_mud_early_tel.rds"))
roach_late_tel <- readRDS(paste0(save_telem_path, "muddyfoot/stage/roach_mud_late_tel.rds"))

### EARLY ###
cl <- makeCluster(15)  
doParallel::registerDoParallel(cl)

# Fit ctmm models in parallel for each roach
roach_mud_early_ctmm_fits <- 
  foreach(i = 1:length(roach_early_tel), .packages = 'ctmm') %dopar% {
    # Generate an initial guess for the model parameters, incorporating location error
    roach_early_guess <- ctmm.guess(roach_early_tel[[i]], CTMM = ctmm(error = TRUE), interactive = FALSE)
    
    # Fit the ctmm model to the telemetry data using model selection
    model_fit <- ctmm.fit(roach_early_tel[[i]], roach_early_guess, method = 'ML')
    
    # Save the fitted model for each fish individually
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_roach_fits/stage/", names(roach_early_tel)[i], ".rds"))
    
    model_fit  # Return the fitted model
  }

# Stop the parallel cluster once model fitting is complete
stopCluster(cl)


### LATE ###
cl <- makeCluster(15)  
doParallel::registerDoParallel(cl)

roach_mud_late_ctmm_fits <- 
  foreach(i = 1:length(roach_late_tel), .packages = 'ctmm') %dopar% {
    # Generate an initial guess for the model parameters, incorporating location error
    roach_late_guess <- ctmm.guess(roach_late_tel[[i]], CTMM = ctmm(error = TRUE), interactive = FALSE)
    
    # Fit the ctmm model to the telemetry data using model selection
    model_fit <- ctmm.fit(roach_late_tel[[i]], roach_late_guess, method = 'ML')
    
    # Save the fitted model for each fish individually
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_roach_fits/stage/", names(roach_late_tel)[i], ".rds"))
    
    model_fit  # Return the fitted model
  }

# Stop the parallel cluster once model fitting is complete
stopCluster(cl)