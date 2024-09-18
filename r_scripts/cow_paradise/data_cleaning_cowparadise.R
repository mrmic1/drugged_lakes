#-------------------------------------------------#
# Initial cleaning of lake cow paradise data 
#-------------------------------------------------#

#--------------#
#> 1. SETUP ####
#--------------#

#LIBRARIES

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
data_trans_path = "./raw_tracking_data/raw/"
save_data_sub_path = "./raw_tracking_data/subsamples/" #where to save sub-sampled data or test data
save_data_clean_path = "./data/tracks_filtered/" #where to save clean data files used in processing and analysis

cat("\014")

#import cow data
lake_cow = fread(paste0(data_trans_path, "cow_paradise_C/results/animal/all.csv"))


#--------------------------#
#> 2. Initial cleaning ####
#-------------------------#

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

#Filter data to  only include fish found within the lake
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

#---------------------------------------------------#
#> 3. Spatial filtering using boundary polygons #####
#---------------------------------------------------#

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

#Thinning data to reduce size
#lake_cow_2s_mv <- lake_cow_mv %>% mt_filter_per_interval(unit = "2 seconds",  criterion="first")# removed around 400,000 points. Not a necessary step, done to make the dataset more manageable. 
#lake_cow_1min_mv <- lake_cow_mv %>% mt_filter_per_interval(unit = "1 minute",  criterion="first")# removed around 8 million points

#Create a polygon or spatial bounds of the lake to determine how many points fall out of this polygon
fish1 <- filter_track_data(lake_cow_mv, .track_id = "59831") #284981 detections for this fish

#Creating map of the lake
lake_cow_map <- mapview::mapView(mt_track_lines(fish1)$geometry) 

#Drawing a polygon around the lake 
lake_cow_polygon <- mapedit::drawFeatures(map = lake_cow_map)

#save this drawing 
sf::st_write(lake_cow_polygon, dsn="./data/lake_coords/lake_cow_polygon.gpkg", driver="GPKG", delete_layer = TRUE) # for overwriting

#load this drawing
#lake_cow_polygon <- sf::st_read("./data/Lakes/lake_cow_polygon.gpkg")

#subsampling points that only fall within the polygon that we drew
#this step takes time to process
lake_cow_sub <- st_intersection(lake_cow_mv, lake_cow_polygon) #full dataset: reduced by only ~400,000
#lake_cow_2s_sub <- st_intersection(lake_cow_2s_mv, lake_cow_polygon) #2s subset: reduced by ~300,000
#lake_cow_1min_sub <- st_intersection(lake_cow_1min_mv, lake_cow_polygon) #1min subset´: reduced by ~100,000

#Save sub-sampled data
saveRDS(lake_cow_sub, file= paste0(save_data_sub_path, "lake_cow_sub.rds")) #removed about 300,000 locations
#saveRDS(lake_cow_2s_sub, file= paste0(save_data_sub_path, "lake_cow_sub_2s.rds"))
#saveRDS(lake_cow_1min_sub, file= paste0(save_data_sub_path, "lake_cow_sub_1min.rds"))

#Compare spatially filtered data to unfiltered data
fish_unfilt <- filter_track_data(lake_cow_mv, .track_id = "59831")
fish_filt <- filter_track_data(lake_cow_sub, .track_id = "59831")

#Plot below takes a fair while
ggplot() + 
  geom_sf(data=lake_cow_polygon) + 
  #geom_sf(data=lake_cow) +
  #geom_sf(data=fish_unfilt, col="red", alpha = 0.3) #+
  geom_sf(data=fish_filt, col="blue", alpha = 0.5)

#looks like the filtering worked

#------------------------------------------------------------------#

#--------------------------------------------#
#> 4. Adding useful metrics to dataframe ####
#--------------------------------------------#

#Read in rds file
#lake_cow_sub <- readRDS(file= paste0(save_data_path, "lake_cow_sub_full.rds"))

#Create new id column for each individual on each day
lake_cow_sub$date <- strftime(lake_cow_sub$timestamp, format="%Y/%m/%d")
lake_cow_sub$week <- strftime(lake_cow_sub$timestamp, format="%W")
lake_cow_sub <- 
  lake_cow_sub %>% 
  unite(individual_week, c("individual_ID", "week"), remove = FALSE)

#How many observations days for each individual
#This code takes some time
#print.data.frame(lake_cow_sub %>% 
#                   dplyr::group_by(individual_ID) %>%
#                   dplyr::summarise(n_days_mon = n_distinct(date)) %>%
#                   ungroup())
#Fish were monitored for 35 total days

## GENERAL DATAFRAME CLEANING ##

#receiver detection data. Do detected columns match used columns
lake_cow_sub %>% 
  filter(RxDetected != RxUsed) %>% 
  nrow()
#211,190

#remove unnecessary columns
lake_cow_sub <- lake_cow_sub %>% 
  dplyr::select(-Station, -HPEm, -TempData, -DepthData, -AccelData, 
                -Tag.Type, -Serial.Number, -PIT.Number, -Biologger.Number, -RxDetected, -nRxDetected)


#rename some columns
names(lake_cow_sub)[names(lake_cow_sub) == "RxUsed"] <- "recievers_used"
names(lake_cow_sub)[names(lake_cow_sub) == "nRxUsed"] <- "num_recievers_used"

#arrange columns more logically
lake_cow_sub <- lake_cow_sub %>% 
  dplyr::select(21,3, 24:25, 15:20, everything())

#retrospectively also delete Date column - as refers to dates that fish were tagged
lake_cow_sub <- lake_cow_sub %>% 
  dplyr::select(-Date)

#add long and lat columns back in the dataframe
coords <- st_coordinates(lake_cow_sub$geometry)
lake_cow_sub$Long <- coords[,1] 
lake_cow_sub$Lat <- coords[, 2]

#Add columns for coordinaate converted to UTM

#Convert data to sf object with WGS84 CRS
data_sf <- st_as_sf(lake_cow_sub, coords = c("Long", "Lat"), crs = 4326) 
#note that crs = 4326 is equivalent to WGS84

#Define the UTM CRS
utm_crs <- st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs")
#zone 34 is the zone that Djacknabole is located within

#Transform coordinates to utm
data_sf_utm <- st_transform(data_sf, crs = utm_crs)

#Extract utm x and y coordinates
utm_coords <- st_coordinates(data_sf_utm$geometry)
lake_cow_sub$Long_utm <- utm_coords[,1]
lake_cow_sub$Lat_utm <- utm_coords[,2]

#Add ID column to work better with ctmm
lake_cow_sub$Fish_ID = paste("F", lake_cow_sub$individual_ID, sep = "")

#Save filtered data
saveRDS(lake_cow_sub, file= paste0(save_data_clean_path, "lake_cow_sub.rds"))
