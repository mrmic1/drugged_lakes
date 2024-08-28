#Script for cleaning the detection data in Lake BT - adapted from Jack Brand script
#Author: Marcus Michelangeli

library(data.table)
library(tidyverse)
library(move)
library(move2)
library(mapedit)
library(sf)
library(ctmm)


#set time zones to CEST (not CET because study was conducted in the summmer)
#helpful timezone link: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
Sys.setenv(TZ = 'Europe/Stockholm')


#create data path for transmitter data
data_trans_path = "./data/Transmitters/raw/"
save_data_sub_path = "./data/Transmitters/subsamples/" #where to save sub-sampled data or test data
save_data_clean_path = "./data/Transmitters/clean/" #where to save clean data files used in processing and analysis

cat("\014")

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

#import BT data
lake_BT = fread(paste0(data_trans_path, "BT_B/results/animal/all.csv"))

#noticed that lat and long are in the wrong columns, need to rename them
names(lake_BT)[names(lake_BT) == "Longitude"] <- "Lat"
names(lake_BT)[names(lake_BT) == "Latitude"] <- "Long"

#import biometic data
biometrics = fread(paste0(data_trans_path, "biometric_data.csv"))

#need to create matching column to merge datasets by ID
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

#For some reason below rename function kept crashing R
#lake_BT <- lake_BT %>%
 # dplyr::rename(indiv_id = Id)

#base function does not crash R
names(lake_BT)[names(lake_BT) == "Id"] <- "individual_ID"


#Join biometric and detection information
lake_BT <- left_join(lake_BT, biometrics, 
                       by = "individual_ID")

lake_BT <- as.data.frame(lake_BT)

# Identify reference individual
lake_BT$individual_ID <- ifelse(lake_BT$FullId == 'H170-1802-65065',
                                'Reference', lake_BT$individual_ID)

#Filter data to  only inlcude fish found within the lake
lake_BT <- lake_BT %>%
  dplyr::filter(Lake == 'BT' | individual_ID == 'Reference')

#Select observations when pike were present in the lake (i.e. after they were introduced)
lake_BT <- lake_BT %>% 
  filter(Time >= "2022-09-26") %>% 
  dplyr::select(-ID, -Notes, -Transmitter, -'Tag Number') #remove unnecessary columns (just a bit of cleaning)

#Change timezone so that it is accurate
lake_BT$timestamp <- with_tz(lake_BT$Time, 'Europe/Stockholm')
lake_BT$timestamp[1] #CEST
lake_BT$Time[1] #UTC

#save
#saveRDS(lake_BT, paste0(save_data_clean_path, "lake_BT.rds"))

#------------------------------------------------------------------------------#

#Reload lake_BT file here

#create a move2 object from a data.frame
#you will need library(move2)
lake_BT_mv <- mt_as_move2(lake_BT, 
                          coords = c("Long","Lat"),
                          crs = "WGS84",
                          time_column = "timestamp",
                          track_id_column = "individual_ID",
                          na.fail = F) # allows or not empty coordinates

#Arrange by individual and datetime
lake_BT_mv <- lake_BT_mv %>% 
  dplyr::arrange(individual_ID, timestamp)

#-------THINNING--------#

#Thinning data to reduce size
lake_BT_2s_mv <- lake_BT_mv %>% mt_filter_per_interval(unit = "2 seconds",  criterion="first")# removed just over 1 million points. Not a necessary step, done to make the dataset more manageable. 
lake_BT_1min_mv <- lake_BT_mv %>% mt_filter_per_interval(unit = "1 minute",  criterion="first")# removed around 19 million points

# --------- SPATIAL FILTERING USING BOUNDING POLYGONS --------- #

#Create a polygon or spatial bounds of the lake to determine how many points fall out of this polygon
fish1 <- filter_track_data(lake_BT_mv, .track_id = "59804") #536959 detections for this fish
#Creating map of the lake
lake_BT_map <- mapview::mapView(mt_track_lines(fish1)$geometry) 
#Drawing a polygon around the lake 
lake_BT_polygon <- mapedit::drawFeatures(map = lake_BT_map)
#save this drawing 
sf::st_write(lake_BT_polygon, dsn="./data/Lakes/lake_BT_polygon.gpkg", driver="GPKG", delete_layer = TRUE) # for overwriting
#load this draqing
lake_BT_polygon <- sf::st_read("./data/Lakes/lake_BT_polygon.gpkg")

#subsampling points that only fall within the polygon that we drew
#this step takes time to process
lake_BT_sub <- st_intersection(lake_BT_mv, lake_BT_polygon) #full dataset: reduced by only ~400,000
lake_BT_2s_sub <- st_intersection(lake_BT_2s_mv, lake_BT_polygon) #2s subset: reduced by ~400,000
lake_BT_1min_sub <- st_intersection(lake_BT_1min_mv, lake_BT_polygon) #1min subset´: reduced by ~100,000

#Save sub-sampled data
saveRDS(lake_BT_sub, file= paste0(save_data_sub_path, "lake_BT_sub_full.rds"))
saveRDS(lake_BT_2s_sub, file= paste0(save_data_sub_path, "lake_BT_sub_2s.rds"))
saveRDS(lake_BT_1min_sub, file= paste0(save_data_sub_path, "lake_BT_sub_1min.rds"))

#Compare spatially filtered data to unfiltered data
fish_unfilt <- filter_track_data(lake_BT_mv, .track_id = "59804")
fish_filt <- filter_track_data(lake_BT_sub, .track_id = "59804")

#Plot below takes a fair while
ggplot() + 
  geom_sf(data=lake_BT_polygon) + 
  #geom_sf(data=lake_BT) +
  geom_sf(data=fish_unfilt, col="red", alpha = 0.3) +
  geom_sf(data=fish_filt, col="blue", alpha = 0.5)

#looks like the filtering worked

#------------------------------------------------------------------#

#Read in rds file
#lake_BT_sub <- readRDS(file= paste0(save_data_path, "lake_BT_sub_full.rds"))

#Creat new id column for each individual on each day
lake_BT_sub$date <- strftime(lake_BT_sub$timestamp, format="%Y/%m/%d")
lake_BT_sub$week <- strftime(lake_BT_sub$timestamp, format="%W")
lake_BT_sub <- 
  lake_BT_sub %>% 
  unite(individual_week, c("individual_ID", "week"), remove = FALSE)

#How many observations days for each individual
#This code takes some time
print.data.frame(lake_BT_sub %>% 
                   dplyr::group_by(individual_ID) %>%
                   dplyr::summarise(n_days_mon = n_distinct(date),
                                    percet_days_mon =  n_distinct(date)/34,
                                    Species = Species) %>%
                   ungroup())
#Fish were monitored for 34 total days


#Split dataset by species
perch_lake_BT <- lake_BT_sub %>%
  filter(Species == 'Perch')
roach_lake_BT <- lake_BT_sub %>% 
  filter(Species == 'Roach')
pike_lake_BT <- lake_BT_sub %>% 
  filter(Species == 'Northern Pike')


# --------- Perch lake_BT ------------- #

#Create a move2 object
perch_lake_BT_mv <- mt_as_move2(perch_lake_BT, time_column = 'timestamp', track_id_column = 'individual_ID')
#Convert to a move bject
perch_lake_BT_mv <- to_move(perch_lake_BT_mv)
perch_lake_BT_tel <- as.telemetry(perch_lake_BT_mv,
                                    timeformat="auto",
                                    timezone="Europe/Stockholm",
                                    projection=NULL,
                                    datum="WGS84",
                                    timeout=Inf,
                                    na.rm="row",
                                    mark.rm=FALSE,
                                    keep=TRUE,
                                    drop=FALSE)
#Outliers using ctmm
out_perch <- outlie(perch_lake_BT_tel)
#how many times does a location exceed speed of 0.977
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))
#Function to get range in speeds (m/s) for each individual. 
#17281 are > than 0.977 m/s -----MUCH HIGHER THAN MUDDYFOOT
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_perch, function(x) x$speed <= 0.977)
#Combining the lists and removing observations for which the logical vector was false
perch_lake_BT_tel <- mapply(function(x,y){x <- x[y,]},
                              x=perch_lake_BT_tel, y=which_lowSp)
#Large matrix format


# --------- Roach lake_BT ------------- #

#Create a move2 object
roach_lake_BT_mv <- mt_as_move2(roach_lake_BT, time_column = 'timestamp', track_id_column = 'individual_ID')
#Convert to a move bject
roach_lake_BT_mv <- to_move(roach_lake_BT_mv)
roach_lake_BT_tel <- as.telemetry(roach_lake_BT_mv,
                                    timeformat="auto",
                                    timezone="Europe/Stockholm",
                                    projection=NULL,
                                    datum="WGS84",
                                    timeout=Inf,
                                    na.rm="row",
                                    mark.rm=FALSE,
                                    keep=TRUE,
                                    drop=FALSE)
#Outliers using ctmm
roach_mud_out <- outlie(roach_lake_BT_tel)
#how many times does a location exceed speed of 0.841 m/s
sum(sapply(roach_mud_out, function(x) sum(x$speed > 0.841)))
#Function to get range in speeds (m/s) for each individual. 
#51356 are > than 0.841 m/s 
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(roach_mud_out, function(x) x$speed <= 0.841)
#Combining the lists and removing observations for which the logical vector was false
roach_lake_BT_tel <- mapply(function(x,y){x <- x[y,]},
                              x=roach_lake_BT_tel, y=which_lowSp)



# --------- Pike lake_BT ------------- #

#Create a move2 object
pike_lake_BT_mv <- mt_as_move2(pike_lake_BT, time_column = 'timestamp', track_id_column = 'individual_ID')
#Convert to a move bject
pike_lake_BT_mv <- to_move(pike_lake_BT_mv)
pike_lake_BT_tel <- as.telemetry(pike_lake_BT_mv,
                                   timeformat="auto",
                                   timezone="Europe/Stockholm",
                                   projection=NULL,
                                   datum="WGS84",
                                   timeout=Inf,
                                   na.rm="row",
                                   mark.rm=FALSE,
                                   keep=TRUE,
                                   drop=FALSE)
#Outliers using ctmm
pike_mud_out <- outlie(pike_lake_BT_tel)
#how many times does a location exceed speed of 0.823 m/s m/s
sum(sapply(pike_mud_out, function(x) sum(x$speed > 0.823)))
#Function to get range in speeds (m/s) for each individual. 
#4070 are > than 0.823 m/s 
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(pike_mud_out, function(x) x$speed <= 0.823)
#Combining the lists and removing observations for which the logical vector was false
pike_lake_BT_tel <- mapply(function(x,y){x <- x[y,]},
                             x=pike_lake_BT_tel, y=which_lowSp)

# --------------- Combine datasets and save -------------------------#

#not sure the below code is useful yet
#lake_BT_telem <- c(perch_lake_BT_tel, roach_lake_BT_tel, pike_lake_BT_tel)  
#save file
#saveRDS(lake_BT_telem, 
#        file=paste0(save_data_path, "lake_BT_telem.rds"))
#load file
#lake_BT_telem <- readRDS(file=paste0(save_data_path, "lake_BT_telem.rds"))


#instead save individual telem files 
saveRDS(perch_lake_BT_tel, file=paste0(save_data_sub_path, "perch_lake_BT_telem.rds"))
saveRDS(roach_lake_BT_tel, file=paste0(save_data_sub_path, "roach_lake_BT_telem.rds"))
saveRDS(pike_lake_BT_tel, file=paste0(save_data_sub_path, "pike_lake_BT_telem.rds"))
