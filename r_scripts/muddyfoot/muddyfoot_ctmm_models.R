# Muddyfoot 
# Run CTMM models for each species 

###  LIBRARIES ###

library(data.table)
library(tidyverse)
library(ctmm)
library(sf)
#for parallel processing
library(parallel)
library(foreach)
library(doParallel)

#set time zones
Sys.setenv(TZ = 'Europe/Stockholm')

#create data path for transmitter data
data_filter_path = "./data/tracks_filtered/"
save_ctmm_path = "./data/ctmm_fits/"
save_telem_path = "./data/telem_obj/"

#-------------------------------------------------------------------------------#

#######################-
#### Northern pike #####
#######################-


#Load in the datasets
muddyfoot_sub <- readRDS(paste0(data_filter_path, 'muddyfoot_sub.rds'))

#isolate pike
pike_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Northern Pike'| individual_ID == 'Reference')


#Eric - movebank method. Provide as.telemetry with columns names it works better with
#No need to coerce it into a move object first.
pike_movebank <- with(pike_muddyfoot, data.frame("timestamp" = timestamp, "location.long" = Long,
                                                 "location.lat" = Lat, "GPS.HDOP" = HPE,
                                                 "individual-local-identifier" = Fish_ID,
                                                 "treatment" = Treatment,
                                                 "date" = date,
                                                 "week" = week,
                                                 "individual_day" = individual_day))

# pike_muddyfoot_tel_utm <- as.telemetry(pike_movebank, 
#                                        timezone = "Europe/Stockholm", 
#                                        timeformat="%Y-%m-%d %H:%M:%S", 
#                                        projection= "+init=epsg:32634",
#                                        datum="WGS84",
#                                        keep = c("treatment", "date", 
#                                                 "week", "individual_day"))

pike_muddyfoot_tel_tpeqd <- as.telemetry(pike_movebank, 
                                         timezone = "Europe/Stockholm", 
                                         timeformat="%Y-%m-%d %H:%M:%S", 
                                         projection= NULL,
                                         datum="WGS84",
                                         keep = c("treatment", "date", 
                                                  "week", "individual_day")
)


#Check some of the parameters

head(pike_muddyfoot_tel_tpeqd$F59880)
#columns x and y are very different from 1. I think it has something to do with the projection
ctmm::projection(pike_muddyfoot_tel_tpeqd$F59880)
#tpeqd projection
tz(pike_muddyfoot_tel_tpeqd$F59880$timestamp)
#"Europe/Stockholm"

#--------------------------------------------------------------------------------#

names(pike_muddyfoot_tel_tpeqd)

#Center the projection on the geometric median of the data
ctmm::projection(pike_muddyfoot_tel_tpeqd) <- ctmm::median(pike_muddyfoot_tel_tpeqd)

### INCORPORATING LOCATION ERROR
# fit error parameters to calibration data
#UERE_utm <- uere.fit(pike_muddyfoot_tel_utm$FReference)
UERE_tpeqd <- uere.fit(pike_muddyfoot_tel_tpeqd$FReference)
# do not run uere.fit on tracking data

#summary(UERE_utm)
summary(UERE_tpeqd)
#both are similar

# apply error model to data
#uere(pike_muddyfoot_tel_utm) <- UERE_utm
uere(pike_muddyfoot_tel_tpeqd) <- UERE_tpeqd
#new column now called VAR.xy

#remove reference list
pike_muddyfoot_tel <- pike_muddyfoot_tel_tpeqd[1:6]

#remove outliers based on speed
out_pike <- outlie(pike_muddyfoot_tel, plot = FALSE)
head(out_pike[[1]])
sum(sapply(out_pike, function(x) sum(x$speed > 0.823)))
#Function to get range in speeds (m/s) for each individual. 
#123469 are > than 0.823 m/s 
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_pike, function(x) x$speed <= 0.823)
#Combining the lists and removing observations for which the logical vector was false
pike_muddyfoot_tel <- Map(function(x,y) x[y,], pike_muddyfoot_tel,which_lowSp)


#save telemetry object
saveRDS(pike_muddyfoot_tel , paste0(save_telem_path, "pike_muddyfoot_tel.rds")) 

#load object
#pike_muddyfoot_tel <- readRDS(paste0(save_telem_path, "pike_muddyfoot_tel.rds"))

#Run ctmm models for six pike in muddyfoot
cl <- makeCluster(6)
doParallel::registerDoParallel(cl)
muddyfoot_pike_ctmm_fits <-  
  foreach(i = 1:length(pike_muddyfoot_tel), .packages = 'ctmm') %dopar% {
    lake_BT_pike_guess <- ctmm.guess(pike_muddyfoot_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.select(pike_muddyfoot_tel[[i]], lake_BT_pike_guess, verbose = TRUE)
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_pike_fits/", names(pike_muddyfoot_tel)[i], ".rds"))
    model_fit
  }

stopCluster(cl)

#add ID to lists
names(muddyfoot_pike_ctmm_fits ) <- names(pike_muddyfoot_tel)

#save fits
saveRDS(muddyfoot_pike_ctmm_fits , paste0(save_ctmm_path, "muddyfoot_pike_ctmm_fits.rds")) 

#Check model outputs
for (i in 1:length(muddyfoot_pike_select_fits)) {
  element_name <- names(muddyfoot_pike_select_fits)[i]
  cat("Summary for", element_name, ":\n")
  print(summary(muddyfoot_pike_select_fits[[i]]))
  cat("\n")
}

#############################-
####### Perch ###############
############################-

#isolate perch
perch_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Perch'| individual_ID == 'Reference')


#Eric - movebank method. Provide as.telemetry with columns names it works better with
#No need to coerce it into a move object first.
perch_movebank <- with(perch_muddyfoot, data.frame("timestamp" = timestamp, "location.long" = Long,
                                                 "location.lat" = Lat, "GPS.HDOP" = HPE,
                                                 "individual-local-identifier" = Fish_ID,
                                                 "treatment" = Treatment,
                                                 "date" = date,
                                                 "week" = week,
                                                 "individual_day" = individual_day))

# perch_muddyfoot_tel_utm <- as.telemetry(perch_movebank, 
#                                        timezone = "Europe/Stockholm", 
#                                        timeformat="%Y-%m-%d %H:%M:%S", 
#                                        projection= "+init=epsg:32634",
#                                        datum="WGS84",
#                                        keep = c("treatment", "date", 
#                                                 "week", "individual_day"))

perch_muddyfoot_tel_tpeqd <- as.telemetry(perch_movebank, 
                                         timezone = "Europe/Stockholm", 
                                         timeformat="%Y-%m-%d %H:%M:%S", 
                                         projection= NULL,
                                         datum="WGS84",
                                         keep = c("treatment", "date", 
                                                  "week", "individual_day")
)


#Check some of the parameters

head(perch_muddyfoot_tel_tpeqd$F59682)
ctmm::projection(perch_muddyfoot_tel_tpeqd$F59682)
tz(perch_muddyfoot_tel_tpeqd$F59682$timestamp)

#--------------------------------------------------------------------------------#

names(perch_muddyfoot_tel_tpeqd)

#Center the projection on the geometric median of the data
ctmm::projection(perch_muddyfoot_tel_tpeqd) <- ctmm::median(perch_muddyfoot_tel_tpeqd)

### INCORPORATING LOCATION ERROR
# fit error parameters to calibration data
#UERE_utm <- uere.fit(perch_muddyfoot_tel_utm$FReference)
UERE_tpeqd <- uere.fit(perch_muddyfoot_tel_tpeqd$FReference)
# do not run uere.fit on tracking data

#summary(UERE_utm)
summary(UERE_tpeqd)
#both are similar

# apply error model to data
#uere(perch_muddyfoot_tel_utm) <- UERE_utm
uere(perch_muddyfoot_tel_tpeqd) <- UERE_tpeqd
#new column now called VAR.xy

#remove reference list
perch_muddyfoot_tel <- perch_muddyfoot_tel_tpeqd[1:30]

#remove outliers based on speed
out_perch <- outlie(perch_muddyfoot_tel, plot = FALSE)
head(out_perch[[1]])
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))
#Function to get range in speeds (m/s) for each individual. 
#156265 are > than 0.977 m/s 
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_perch, function(x) x$speed <= 0.977)
#Combining the lists and removing observations for which the logical vector was false
perch_muddyfoot_tel <- Map(function(x,y) x[y,], perch_muddyfoot_tel,which_lowSp)

#save telemetry object
#saveRDS(perch_muddyfoot_tel , paste0(save_telem_path, "perch_muddyfoot_tel.rds")) 
perch_muddyfoot_tel <- readRDS(paste0(save_telem_path, "perch_muddyfoot_tel.rds"))

#load object
#perch_muddyfoot_tel <- readRDS(paste0(save_telem_path, "perch_muddyfoot_tel.rds"))

#some individuals did not get OUF. Rerun for these individuals
# List of desired identities
desired_identities <- c("F59682", "F59688", "F59689", "F59698", "F59699", "F59708",
                        "F59711", "F59712", "F59714", "F59717", "F59728", "F59730",
                        "F59733", "F59734")

# Extract names from the list
object_names <- names(perch_muddyfoot_tel)

# Use a regular expression to extract the desired parts of the names
extracted_names <- sub(".*\\$", "", object_names)

# Filter the list based on the extracted names
filtered_list <- perch_muddyfoot_tel[extracted_names %in% desired_identities]

# Print the filtered list
print(filtered_list)

#Run ctmm models for 14 perch in muddyfoot
cl <- makeCluster(14)
doParallel::registerDoParallel(cl)
muddyfoot_perch_ctmm_fits_OUF <-  
  foreach(i = 1:length(filtered_list), .packages = 'ctmm') %dopar% {
    muddyfoot_perch_guess <- ctmm.guess(filtered_list[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model_fit <- ctmm.select(filtered_list[[i]], muddyfoot_perch_guess, verbose = TRUE)
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_perch_fits/", names(filtered_list)[i], "_OUF.rds"))
    model_fit
  }

stopCluster(cl)

#add ID to lists
names(muddyfoot_perch_ctmm_fits_OUF ) <- names(filtered_list)

#save fits
saveRDS(muddyfoot_perch_ctmm_fits_OUF , paste0(save_ctmm_path, "muddyfoot_perch_ctmm_fits_OUF.rds")) 

#Check model outputs
for (i in 1:length(muddyfoot_perch_select_fits)) {
  element_name <- names(muddyfoot_perch_select_fits)[i]
  cat("Summary for", element_name, ":\n")
  print(summary(muddyfoot_perch_select_fits[[i]]))
  cat("\n")
}

#############################-
####### Roach ###############
############################-

#isolate roach
roach_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Roach'| individual_ID == 'Reference')


#Eric - movebank method. Provide as.telemetry with columns names it works better with
#No need to coerce it into a move object first.
roach_movebank <- with(roach_muddyfoot, data.frame("timestamp" = timestamp, "location.long" = Long,
                                                   "location.lat" = Lat, "GPS.HDOP" = HPE,
                                                   "individual-local-identifier" = Fish_ID,
                                                   "treatment" = Treatment,
                                                   "date" = date,
                                                   "week" = week,
                                                   "individual_day" = individual_day))

# roach_muddyfoot_tel_utm <- as.telemetry(roach_movebank, 
#                                        timezone = "Europe/Stockholm", 
#                                        timeformat="%Y-%m-%d %H:%M:%S", 
#                                        projection= "+init=epsg:32634",
#                                        datum="WGS84",
#                                        keep = c("treatment", "date", 
#                                                 "week", "individual_day"))

roach_muddyfoot_tel_tpeqd <- as.telemetry(roach_movebank, 
                                          timezone = "Europe/Stockholm", 
                                          timeformat="%Y-%m-%d %H:%M:%S", 
                                          projection= NULL,
                                          datum="WGS84",
                                          keep = c("treatment", "date", 
                                                   "week", "individual_day")
)


# #check our telemetery objects
# head(roach_muddyfoot_tel_utm$F59880)
# ctmm::projection(roach_muddyfoot_tel_utm) 
# #utm projection
# tz(roach_muddyfoot_tel_utm$F59880$timestamp)
# #"Europe/Stockholm"

head(roach_muddyfoot_tel_tpeqd$F59880)
#columns x and y are very different from 1. I think it has something to do with the projection
ctmm::projection(roach_muddyfoot_tel_tpeqd$F59880)
#tpeqd projection
tz(roach_muddyfoot_tel_tpeqd$F59880$timestamp)
#"Europe/Stockholm"

#--------------------------------------------------------------------------------#

names(roach_muddyfoot_tel_tpeqd)

#Center the projection on the geometric median of the data
ctmm::projection(roach_muddyfoot_tel_tpeqd) <- ctmm::median(roach_muddyfoot_tel_tpeqd)

### INCORPORATING LOCATION ERROR
# fit error parameters to calibration data
#UERE_utm <- uere.fit(roach_muddyfoot_tel_utm$FReference)
UERE_tpeqd <- uere.fit(roach_muddyfoot_tel_tpeqd$FReference)
# do not run uere.fit on tracking data

#summary(UERE_utm)
summary(UERE_tpeqd)
#both are similar

# apply error model to data
#uere(roach_muddyfoot_tel_utm) <- UERE_utm
uere(roach_muddyfoot_tel_tpeqd) <- UERE_tpeqd
#new column now called VAR.xy
head(roach_muddyfoot_tel_tpeqd$F59684)
names(roach_muddyfoot_tel_tpeqd)

#remove reference list
roach_muddyfoot_tel <- roach_muddyfoot_tel_tpeqd[1:29]

#remove outliers based on speed
out_roach <- outlie(roach_muddyfoot_tel, plot = FALSE)
head(out_roach[[1]])
sum(sapply(out_roach, function(x) sum(x$speed > 0.841)))
#Function to get range in speeds (m/s) for each individual. 
#336192 are > than 0.841 m/s 
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_roach, function(x) x$speed <= 0.841)
#Combining the lists and removing observations for which the logical vector was false
roach_muddyfoot_tel <- Map(function(x,y) x[y,], roach_muddyfoot_tel,which_lowSp)

#save telem object
saveRDS(roach_muddyfoot_tel , paste0(save_telem_path, "roach_muddyfoot_tel.rds"))


cl <- makeCluster(28)
doParallel::registerDoParallel(cl)
muddyfoot_roach_select_fits <-  
  foreach(i = 1:length(roach_muddyfoot_tel), .packages = 'ctmm') %dopar% {
    muddyfoot_roach_guess <- ctmm.guess(roach_muddyfoot_tel[[i]], CTMM=ctmm(error=TRUE), interactive = FALSE)
    model.fit <- ctmm.select(roach_muddyfoot_tel[[i]], muddyfoot_roach_guess, verbose = TRUE)
    saveRDS(model_fit, file = paste0(save_ctmm_path, "muddyfoot_roach_fits/", names(roach_muddyfoot_tel)[i], ".rds"))
    model_fit
  }

stopCluster(cl)

names(muddyfoot_roach_select_fits ) <- names(roach_muddyfoot_tel)
saveRDS(muddyfoot_roach_select_fits , file = paste0(save_ctmm_path, "muddyfoot_roach_fits/", "muddyfoot_roach_ctmm_fits.rds")) 
#mud_roach_fits <- readRDS(paste0(save_ctmm_path, "muddyfoot_roach_ctmm_fits.rds"))




