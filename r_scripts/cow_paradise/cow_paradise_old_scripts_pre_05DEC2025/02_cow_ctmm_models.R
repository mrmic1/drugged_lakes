#--------------------------------------#
# Run CTMM models for Lake Cow Paradise #
#--------------------------------------#

# > 1. Setup ####

###  LIBRARIES ###

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
filtered_data_path <- "./data/tracks_filtered/lake_cow_paradise/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"

#-------------------------------------------------------------------------------#

#> 2. Create telemetry objects

#Load in the datasets
lake_cow_sub <- readRDS(paste0(filtered_data_path, '01_lake_cow_sub.rds'))

#Eric - movebank method. Provide as.telemetry with columns names it works better with
#No need to coerce it into a move object first.
lake_cow_movebank <- 
  with(lake_cow_sub, 
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

lake_cow_tel_tpeqd <- as.telemetry(lake_cow_movebank, 
                                   timezone = "Europe/Stockholm", 
                                   timeformat="%Y-%m-%d %H:%M:%S", 
                                   projection= NULL,
                                   datum="WGS84",
                                   keep = c("Species","Weight","Total_length", "Std_length", "Treatment", 
                                            "Date", "Exp_Stage", "Time_Of_Day"))


#Check some of the parameters
head(lake_cow_tel_tpeqd$F59747)
#columns x and y are very different from 1. I think it has something to do with the projection
ctmm::projection(lake_cow_tel_tpeqd$F59747)
#tpeqd projection
tz(lake_cow_tel_tpeqd$F59747$timestamp)
#"Europe/Stockholm"

#--------------------------------------------------------------------------------#

#> 3. Incoporate location error ####

names(lake_cow_tel_tpeqd)

#Center the projection on the geometric median of the data
ctmm::projection(lake_cow_tel_tpeqd) <- ctmm::median(lake_cow_tel_tpeqd)

### INCORPORATING LOCATION ERROR
# fit error parameters to calibration data
#UERE_utm <- uere.fit(pike_lake_cow_tel_utm$FReference)
UERE_tpeqd <- uere.fit(lake_cow_tel_tpeqd$FReference)
# do not run uere.fit on tracking data

#summary(UERE_utm)
summary(UERE_tpeqd)
#both are similar

# apply error model to data
#uere(pike_lake_cow_tel_utm) <- UERE_utm
uere(lake_cow_tel_tpeqd) <- UERE_tpeqd
#new column now called VAR.xy

#remove reference list
names(lake_cow_tel_tpeqd)
lake_cow_tel <- lake_cow_tel_tpeqd[1:66]

#-----------------------------------------------------------------------------#

#> 4. Remove outliers based on species maximum swim speeds ####

#Pike = 0.823
#Perch = 0.977
#Roach = 0.841
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES)

#first seperate telemetry object by species 
#check whether ids are in species order
lake_cow_movebank %>%
  select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

#need to order them by id number
lake_cow_tel <- lake_cow_tel[order(names(lake_cow_tel))]
names(lake_cow_tel)

#mostly in order except for the first perch
pike_lake_cow_tel <- lake_cow_tel[61:66]
perch_lake_cow_tel <- lake_cow_tel[c(1,22:60)]
roach_lake_cow_tel <- lake_cow_tel[2:21]

#remove outliers based on speed
#pike
out_pike <- outlie(pike_lake_cow_tel, plot = FALSE)
sum(sapply(out_pike, function(x) sum(x$speed > 0.823)))
#26812

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_pike, function(x) x$speed <= 0.823)
#Combining the lists and removing observations for which the logical vector was false
pike_lake_cow_tel <- Map(function(x,y) x[y,], pike_lake_cow_tel,which_lowSp)
#save telemetry object
saveRDS(pike_lake_cow_tel , paste0(save_telem_path, "cow_paradise/pike_lake_cow_tel.rds")) 

#Perch
out_perch <- outlie(perch_lake_cow_tel, plot = FALSE)
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))
#Function to get range in speeds (m/s) for each individual. 
#58011 are > than 0.823 m/s 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_perch, function(x) x$speed <= 0.977)
#Combining the lists and removing observations for which the logical vector was false
perch_lake_cow_tel <- Map(function(x,y) x[y,], perch_lake_cow_tel,which_lowSp)
#save telemetry object
saveRDS(perch_lake_cow_tel , paste0(save_telem_path, "cow_paradise/perch_lake_cow_tel.rds")) 

#Roach
out_roach <- outlie(roach_lake_cow_tel, plot = FALSE)
head(out_roach[[1]])
sum(sapply(out_roach, function(x) sum(x$speed > 0.841)))
#Function to get range in speeds (m/s) for each individual. 
#128508 are > than 0.841 m/s 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_roach, function(x) x$speed <= 0.841)
#Combining the lists and removing observations for which the logical vector was false
roach_lake_cow_tel <- Map(function(x,y) x[y,], roach_lake_cow_tel,which_lowSp)
#save telemetry object
saveRDS(roach_lake_cow_tel , paste0(save_telem_path, "cow_paradise/roach_lake_cow_tel.rds"))



#> 5. Make combined dataframe with outliers for each species removed ####


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
pike_data <- extract_data(pike_lake_cow_tel)
perch_data <- extract_data(perch_lake_cow_tel)
roach_data <- extract_data(roach_lake_cow_tel)

# Combine all dataframes into one, with a species label added.
cow_telem_data <- 
  rbind(
    cbind(pike_data, Species = 'Pike'), 
    cbind(perch_data, Species = 'Perch'),
    cbind(roach_data, Species = 'Roach')
  )

saveRDS(cow_telem_data, paste0(filtered_data_path, "02_lake_cow_sub.rds"))

#------------------------------------------------------------------------------#

#> 5. Run CTMM models ####

### >>> Pike ####
#Run ctmm models for six pike in lake_cow
names(pike_lake_cow_tel)

cl <- makeCluster(6)
doParallel::registerDoParallel(cl)
lake_cow_pike_ctmm_fits <-  
  foreach(i = 1:length(pike_lake_cow_tel), .packages = 'ctmm') %dopar% {
    lake_BT_pike_guess <- ctmm.guess(pike_lake_cow_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.select(pike_lake_cow_tel[[i]], lake_BT_pike_guess, verbose = TRUE, method = 'ML')
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_cow_pike_fits/", names(pike_lake_cow_tel)[i], ".rds"))
    model_fit
  }

stopCluster(cl)


#add ID to lists
names(lake_cow_pike_ctmm_fits ) <- names(pike_lake_cow_tel)


#Did not finish F59893
#Need to run seperately
F59893_guess <- ctmm.guess(pike_lake_cow_tel[[1]], CTMM=ctmm(error=TRUE), interactive = FALSE)
F59893_fit <- ctmm.fit(pike_lake_cow_tel[[1]], F59893_guess, method = 'ML', trace = TRUE)
summary(F59893_fit)
lake_cow_sub %>%
  filter(individual_ID == '59893') %>% 
  nrow()
#6086 - low number of positions
lake_cow_sub %>%
  filter(individual_ID == '59893') %>% 
  summarise(n_dates = length(unique(date)))
#tracked for 12 days

#save fits
saveRDS(F59893_fit , paste0(save_ctmm_path, "lake_cow_pike_fits/F59893.rds")) 

#Check model outputs
for (i in 1:length(lake_cow_pike_select_fits)) {
  element_name <- names(lake_cow_pike_select_fits)[i]
  cat("Summary for", element_name, ":\n")
  print(summary(lake_cow_pike_select_fits[[i]]))
  cat("\n")
}



### >>> Perch ###############

#Run ctmm models for six perch in lake_cow
names(perch_lake_cow_tel) #40 individuals

cl <- makeCluster(20)
doParallel::registerDoParallel(cl)
lake_cow_perch_ctmm_fits <- list()
lake_cow_perch_ctmm_fits <-  
  foreach(i = 1:length(perch_lake_cow_tel), .packages = 'ctmm') %dopar% {
    lake_BT_perch_guess <- ctmm.guess(perch_lake_cow_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.select(perch_lake_cow_tel[[i]], lake_BT_perch_guess, verbose = TRUE, method = 'ML')
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_cow_perch_fits/", names(perch_lake_cow_tel)[i], ".rds"))
    model_fit
  }

stopCluster(cl)


#add ID to lists
names(lake_cow_perch_ctmm_fits ) <- names(perch_lake_cow_tel)

#save fits
saveRDS(lake_cow_perch_ctmm_fits , paste0(save_ctmm_path, "lake_cow_perch_fits/lake_cow_perch_ctmm_fits.rds")) 

###>>> Roach ###############

names(roach_lake_cow_tel) #40 individuals

cl <- makeCluster(15)
doParallel::registerDoParallel(cl)
lake_cow_roach_ctmm_fits <-  
  foreach(i = 1:length(roach_lake_cow_tel), .packages = 'ctmm') %dopar% {
    lake_BT_roach_guess <- ctmm.guess(roach_lake_cow_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.select(roach_lake_cow_tel[[i]], lake_BT_roach_guess, verbose = TRUE, method = 'ML')
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_cow_roach_fits/", names(roach_lake_cow_tel)[i], ".rds"))
    model_fit
  }

stopCluster(cl)

names(roach_lake_cow_tel)

names(roach_lake_cow_tel) #20 individuals

cl <- makeCluster(15)
doParallel::registerDoParallel(cl)
lake_cow_roach_ctmm_fits <-  
  foreach(i = 1:length(roach_lake_cow_tel), .packages = 'ctmm') %dopar% {
    lake_BT_roach_guess <- ctmm.guess(roach_lake_cow_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.select(roach_lake_cow_tel[[i]], lake_BT_roach_guess, verbose = TRUE, method = 'ML')
    saveRDS(model_fit, file = paste0(save_ctmm_path, "lake_cow_roach_fits/", names(roach_lake_cow_tel)[i], ".rds"))
    model_fit
  }

stopCluster(cl)

names(roach_lake_cow_tel)

#Did not finish F59819
#Need to run seperately
F59819_guess <- ctmm.guess(roach_lake_cow_tel[[3]], CTMM=ctmm(error=TRUE), interactive = FALSE)
F59819_fit <- ctmm.fit(roach_lake_cow_tel[[3]], F59893_guess, method = 'ML', trace = TRUE)
summary(F59819_fit)
lake_cow_sub %>%
  filter(individual_ID == '59818') %>% 
  nrow()
#325678 
lake_cow_sub %>%
  filter(individual_ID == '59818') %>% 
  summarise(n_dates = length(unique(date)))
#tracked for 35 days

#save fits
saveRDS(F59819_fit , paste0(save_ctmm_path, "lake_cow_roach_fits/F59819.rds")) 

#Check model outputs
for (i in 1:length(lake_cow_pike_select_fits)) {
  element_name <- names(lake_cow_pike_select_fits)[i]
  cat("Summary for", element_name, ":\n")
  print(summary(lake_cow_pike_select_fits[[i]]))
  cat("\n")
}