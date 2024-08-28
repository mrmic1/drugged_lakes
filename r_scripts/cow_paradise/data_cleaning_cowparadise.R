#Script for preliminary cleaning the detection data in Lake Cow Paradise - adapted from Jack Brand script
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

#import cow data
lake_cow = fread(paste0(data_trans_path, "cow_paradise_C/results/animal/all.csv"))

#noticed that lat and long are in the wrong columns, need to rename them
names(lake_cow)[names(lake_cow) == "Longitude"] <- "Lat"
names(lake_cow)[names(lake_cow) == "Latitude"] <- "Long"

#import biometic data
biometrics = fread(paste0(data_trans_path, "biometric_data.csv"))

#need to create matching column to merge datasets by ID
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

#For some reason below rename function kept crashing R
#lake_cow <- lake_cow %>%
# dplyr::rename(indiv_id = Id)

#base function does not crash R
names(lake_cow)[names(lake_cow) == "Id"] <- "individual_ID"


#Join biometric and detection information
lake_cow <- left_join(lake_cow, biometrics, 
                     by = "individual_ID")

lake_cow <- as.data.frame(lake_cow)

# Identify reference individual
lake_cow$individual_ID <- ifelse(lake_cow$FullId == 'H170-1802-65064',
                                'Reference', lake_cow$individual_ID)

#Filter data to  only inlcude fish found within the lake
lake_cow <- lake_cow %>%
  dplyr::filter(Lake == 'Cow Paradise' | individual_ID == 'Reference')

#Select observations when pike were present in the lake (i.e. after they were introduced)
lake_cow <- lake_cow %>% 
  filter(Time >= "2022-09-27") %>% 
  dplyr::select(-ID, -Notes, -Transmitter, -'Tag Number') #remove unnecessary columns (just a bit of cleaning)

#Change timezone so that it is accurate
lake_cow$timestamp <- with_tz(lake_cow$Time, 'Europe/Stockholm')
lake_cow$timestamp[1] #CEST
lake_cow$Time[1] #UTC

#save
#saveRDS(lake_cow, paste0(save_data_clean_path, "lake_cow.rds"))

#------------------------------------------------------------------------------#

#Reload lake_cow file here


#create a move2 object from a data.frame
#you will need library(move2)
lake_cow_mv <- mt_as_move2(lake_cow, 
                          coords = c("Long","Lat"),
                          crs = "WGS84",
                          time_column = "timestamp",
                          track_id_column = "individual_ID",
                          na.fail = F) # allows or not empty coordinates

#Arrange by individual and datetime
lake_cow_mv <- lake_cow_mv %>% 
  dplyr::arrange(individual_ID, timestamp)

#-------THINNING--------#

#Thinning data to reduce size
lake_cow_2s_mv <- lake_cow_mv %>% mt_filter_per_interval(unit = "2 seconds",  criterion="first")# removed around 400,000 points. Not a necessary step, done to make the dataset more manageable. 
lake_cow_1min_mv <- lake_cow_mv %>% mt_filter_per_interval(unit = "1 minute",  criterion="first")# removed around 8 million points

# --------- SPATIAL FILTERING USING BOUNDING POLYGONS --------- #

#Create a polygon or spatial bounds of the lake to determine how many points fall out of this polygon
fish1 <- filter_track_data(lake_cow_mv, .track_id = "59831") #284981 detections for this fish
#Creating map of the lake
lake_cow_map <- mapview::mapView(mt_track_lines(fish1)$geometry) 
#Drawing a polygon around the lake 
lake_cow_polygon <- mapedit::drawFeatures(map = lake_cow_map)
#save this drawing 
sf::st_write(lake_cow_polygon, dsn="./data/Lakes/lake_cow_polygon.gpkg", driver="GPKG", delete_layer = TRUE) # for overwriting
#load this drawing
#lake_cow_polygon <- sf::st_read("./data/Lakes/lake_cow_polygon.gpkg")

#subsampling points that only fall within the polygon that we drew
#this step takes time to process
lake_cow_sub <- st_intersection(lake_cow_mv, lake_cow_polygon) #full dataset: reduced by only ~400,000
lake_cow_2s_sub <- st_intersection(lake_cow_2s_mv, lake_cow_polygon) #2s subset: reduced by ~300,000
lake_cow_1min_sub <- st_intersection(lake_cow_1min_mv, lake_cow_polygon) #1min subset´: reduced by ~100,000

#Save sub-sampled data
saveRDS(lake_cow_sub, file= paste0(save_data_sub_path, "lake_cow_sub_full.rds"))
saveRDS(lake_cow_2s_sub, file= paste0(save_data_sub_path, "lake_cow_sub_2s.rds"))
saveRDS(lake_cow_1min_sub, file= paste0(save_data_sub_path, "lake_cow_sub_1min.rds"))

#Compare spatially filtered data to unfiltered data
fish_unfilt <- filter_track_data(lake_cow_mv, .track_id = "59831")
fish_filt <- filter_track_data(lake_cow_sub, .track_id = "59831")

#Plot below takes a fair while
ggplot() + 
  geom_sf(data=lake_cow_polygon) + 
  #geom_sf(data=lake_cow) +
  geom_sf(data=fish_unfilt, col="red", alpha = 0.3) #+
  #geom_sf(data=fish_filt, col="blue", alpha = 0.5)

#looks like the filtering worked

#------------------------------------------------------------------#

#Read in rds file
#lake_cow_sub <- readRDS(file= paste0(save_data_path, "lake_cow_sub_full.rds"))

#Creat new id column for each individual on each day
lake_cow_sub$date <- strftime(lake_cow_sub$timestamp, format="%Y/%m/%d")
lake_cow_sub$week <- strftime(lake_cow_sub$timestamp, format="%W")
lake_cow_sub <- 
  lake_cow_sub %>% 
  unite(individual_week, c("individual_ID", "week"), remove = FALSE)

#How many observations days for each individual
#This code takes some time
print.data.frame(lake_cow_sub %>% 
                   dplyr::group_by(individual_ID) %>%
                   dplyr::summarise(n_days_mon = n_distinct(date)) %>%
                   ungroup())
#Fish were monitored for 35 total days


#Split dataset by species
perch_lake_cow <- lake_cow_sub %>%
  filter(Species == 'Perch')
roach_lake_cow <- lake_cow_sub %>% 
  filter(Species == 'Roach')
pike_lake_cow <- lake_cow_sub %>% 
  filter(Species == 'Northern Pike')


# --------- Perch lake_cow ------------- #

#Create a move2 object
perch_lake_cow_mv <- mt_as_move2(perch_lake_cow, time_column = 'timestamp', track_id_column = 'individual_ID')
#Convert to a move bject
perch_lake_cow_mv <- to_move(perch_lake_cow_mv)
perch_lake_cow_tel <- as.telemetry(perch_lake_cow_mv,
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
out_perch <- outlie(perch_lake_cow_tel)
#how many times does a location exceed speed of 0.977
sum(sapply(out_perch, function(x) sum(x$speed > 0.977)))
#Function to get range in speeds (m/s) for each individual. 
#10100 are > than 0.977 m/s -----MUCH HIGHER THAN MUDDYFOOT
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(out_perch, function(x) x$speed <= 0.977)
#Combining the lists and removing observations for which the logical vector was false
perch_lake_cow_tel <- mapply(function(x,y){x <- x[y,]},
                            x=perch_lake_cow_tel, y=which_lowSp)
#Large matrix format


# --------- Roach lake_cow ------------- #

#Create a move2 object
roach_lake_cow_mv <- mt_as_move2(roach_lake_cow, time_column = 'timestamp', track_id_column = 'individual_ID')
#Convert to a move bject
roach_lake_cow_mv <- to_move(roach_lake_cow_mv)
roach_lake_cow_tel <- as.telemetry(roach_lake_cow_mv,
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
roach_mud_out <- outlie(roach_lake_cow_tel)
#how many times does a location exceed speed of 0.841 m/s
sum(sapply(roach_mud_out, function(x) sum(x$speed > 0.841)))
#Function to get range in speeds (m/s) for each individual. 
#16506 are > than 0.841 m/s 
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(roach_mud_out, function(x) x$speed <= 0.841)
#Combining the lists and removing observations for which the logical vector was false
roach_lake_cow_tel <- mapply(function(x,y){x <- x[y,]},
                            x=roach_lake_cow_tel, y=which_lowSp)



# --------- Pike lake_cow ------------- #

#Create a move2 object
pike_lake_cow_mv <- mt_as_move2(pike_lake_cow, time_column = 'timestamp', track_id_column = 'individual_ID')
#Convert to a move bject
pike_lake_cow_mv <- to_move(pike_lake_cow_mv)
pike_lake_cow_tel <- as.telemetry(pike_lake_cow_mv,
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
pike_mud_out <- outlie(pike_lake_cow_tel)
#how many times does a location exceed speed of 0.823 m/s m/s
sum(sapply(pike_mud_out, function(x) sum(x$speed > 0.823)))
#Function to get range in speeds (m/s) for each individual. 
#1231 are > than 0.823 m/s 
#Ucrit speeds taken from 
#KEY FACTORS EXPLAINING CRITICAL SWIMMING SPEED IN FRESHWATER FISH:  
#A REVIEW AND STATISTICAL ANALYSIS USING IBERIAN SPECIES), 

#Need to filter out unrealistic speeds
#Making a logical vector
which_lowSp <- lapply(pike_mud_out, function(x) x$speed <= 0.823)
#Combining the lists and removing observations for which the logical vector was false
pike_lake_cow_tel <- mapply(function(x,y){x <- x[y,]},
                           x=pike_lake_cow_tel, y=which_lowSp)

# --------------- Combine datasets and save -------------------------#

#not sure the below code is useful yet
#lake_cow_telem <- c(perch_lake_cow_tel, roach_lake_cow_tel, pike_lake_cow_tel)  
#save file
#saveRDS(lake_cow_telem, 
#        file=paste0(save_data_path, "lake_cow_telem.rds"))
#load file
#lake_cow_telem <- readRDS(file=paste0(save_data_path, "lake_cow_telem.rds"))


#instead save individual telem files 
saveRDS(perch_lake_cow_tel, file=paste0(save_data_sub_path, "perch_lake_cow_telem.rds"))
saveRDS(roach_lake_cow_tel, file=paste0(save_data_sub_path, "roach_lake_cow_telem.rds"))
saveRDS(pike_lake_cow_tel, file=paste0(save_data_sub_path, "pike_lake_cow_telem.rds"))
