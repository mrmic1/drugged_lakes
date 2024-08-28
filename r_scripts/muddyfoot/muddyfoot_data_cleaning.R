#Script for cleaning the detection data in Lake Muddyfoot - adapted from Jack B script
#Author: Marcus Michelangeli

library(data.table)
library(tidyverse)
library(move)
library(move2)
library(mapedit)
library(sf)
library(ctmm)


#set time zones
Sys.setenv(TZ = 'Europe/Stockholm')

#create data path for transmitter data
data_trans_path = "./data/Transmitters/raw/"
save_data_sub_path = "./data/Transmitters/subsamples/" #where to save sub-sampled data or test data
save_data_clean_path = "./data/Transmitters/clean/" #where to save clean data files used in processing and analysis

cat("\014")

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#import muddyfoot data
muddyfoot = fread(paste0(data_trans_path, "muddyfoot_A/results/animal/all.csv"))


#renaming lat and long columns to align with other lakes (see other scripts for reason for change)
names(muddyfoot)[names(muddyfoot) == "Longitude"] <- "Long"
names(muddyfoot)[names(muddyfoot) == "Latitude"] <- "Lat"

#import biometic data
biometrics = fread(paste0(data_trans_path, "biometric_data.csv"))

#need to create matching column to merge datasets by ID
biometrics <- biometrics %>%
  mutate(individual_ID = as.numeric(sub(".*-.*-(\\d+)", "\\1", `Tag Number`)))

# muddyfoot <- muddyfoot %>%
#   dplyr::rename(individual_ID = Id)

names(muddyfoot)[names(muddyfoot) == "Id"] <- "individual_ID"


#Join biometric and detection information
muddyfoot <- left_join(muddyfoot, biometrics, 
                       by = "individual_ID")

muddyfoot <- as.data.frame(muddyfoot)

# Identify reference individual
muddyfoot$individual_ID <- ifelse(muddyfoot$FullId == 'H170-1802-65066', 
                                  'Reference', muddyfoot$individual_ID)

#Filter data to  only include fish found within the lake
muddyfoot <- muddyfoot %>%
  dplyr::filter(Lake == 'Muddyfoot' | individual_ID == 'Reference')

#Select observations when pike were present in the lake (i.e. after they were introduced)
muddyfoot <- muddyfoot %>% 
  filter(Time >= "2022-09-25") %>% 
  dplyr::select(-ID, -Notes, -Transmitter, -'Tag Number') #remove unnecessary columns (just a bit of initial cleaning)

#Change timezone so that it is accurate
muddyfoot$timestamp <- with_tz(muddyfoot$Time, "Europe/Stockholm")
muddyfoot$timestamp[1] #CEST
muddyfoot$Time[1] #UTC which is CEST - 2hrs

#create a move2 object from a data.frame
#you will need library(move2)
muddyfoot_mv <- mt_as_move2(muddyfoot, 
                            coords = c("Long","Lat"),
                            crs = "WGS84",
                            time_column = "timestamp",
                            track_id_column = "individual_ID",
                            na.fail = F) # allows or not empty coordinates

#Arrange by individual and datetime
muddyfoot_mv <- muddyfoot_mv %>% 
  dplyr::arrange(individual_ID, timestamp)

#Removing duplicate locations and timestamps
#muddyfoot_mv <- muddyfoot_mv %>% mt_filter_unique(criterion = "first") 
#does not remove anything from the dataset


# #Create a polygon or spatial bounds of the lake to determine how many points fall out of this polygon
# fish1 <- filter_track_data(muddyfoot_mv, .track_id = "59687") #259930 detections for this fish
# #Creating map of the lake
# muddyfoot_map <- mapview::mapView(mt_track_lines(fish1)$geometry) 
# #Drawing a polygon around the lake 
# muddyfoot_lake_poly <- mapedit::drawFeatures(map = muddyfoot_map)
# #save this drawing 
# #sf::st_write(muddyfoot_lake_poly, dsn="./data/Lakes/lake_muddyfoot_polygon.gpkg", driver="GPKG", delete_layer = TRUE) # for overwriting

#load this draqing
muddyfoot_lake_poly <- sf::st_read("./data/Lakes/lake_muddyfoot_polygon.gpkg")

#subsampling points that only fall within the polygon that we drew
#this step takes time to process
muddyfoot_sub <- st_intersection(muddyfoot_mv, muddyfoot_lake_poly) 

#Save filtered data
#saveRDS(muddyfoot_sub, file= paste0(save_data_sub_path, "muddyfoot_sub_spatialfilt.rds"))

#Thinning data to reduce size
#muddyfoot_sub <- muddyfoot_sub %>% mt_filter_per_interval(unit = "25 seconds",  criterion="first")
#muddyfoot_sub_2 <- muddyfoot_sub_2 %>% mt_filter_per_interval(unit = "25 seconds",  criterion="first")
#removed ~5 million points

# #check that all location fall within map polygon
# #they do, no need to re-run code
# ggplot() + 
#   geom_sf(data=muddyfoot_lake) + 
#   #geom_sf(data=muddyfoot) +
#   geom_sf(data=muddyfoot_sub, col="red", alpha = 0.5)

### TIME DIFFERENCES BETWEEN LOCATIONS ###

#Add column calculating the time difference between timestamps
muddyfoot_sub <- 
  muddyfoot_sub %>%
  group_by(individual_ID) %>% 
  mutate(time_diff = c(NA, diff(timestamp)))
  

muddyfoot_sub$time_diff <- as.numeric(round(muddyfoot_sub$time_diff, digits = 3))

#check time_diff
head(muddyfoot_sub %>% 
       dplyr::select(timestamp, time_diff), n = 20)

#Creat new id column for each individual on each day
muddyfoot_sub$date <- strftime(muddyfoot_sub$timestamp, format="%Y/%m/%d")
muddyfoot_sub$week <- strftime(muddyfoot_sub$timestamp, format="%W")
muddyfoot_sub$day <- strftime(muddyfoot_sub$timestamp, format="%j")
muddyfoot_sub <- 
  muddyfoot_sub %>% 
  mutate(individual_day = paste(individual_ID, day, sep = "_"))

# #How many observations days for each individual
# print.data.frame(muddyfoot_sub %>% 
#                    dplyr::group_by(individual_ID) %>%
#                    dplyr::summarise(n_days_mon = n_distinct(date),
#                                     percet_days_mon =  n_distinct(date)/36) %>%
#                    ungroup())
# #Fish were monitored for 36 total days

# #one fish was only detected for one day. lets check
# odd_fish = muddyfoot_sub %>% 
#   filter(individual_ID == '59707')

 

# time_diff_stats_muddy =
#   muddyfoot_sub %>% 
#   group_by(Species, date) %>% 
#   summarise(mean_time_diff = mean(time_diff, na.rm = TRUE),
#             sd_time_diff = sd(time_diff, na.rm = TRUE),
#             median_time_diff = median(time_diff, na.rm = TRUE),
#             min_time_diff = min(time_diff, na.rm = TRUE),
#             max_time_diff = max(time_diff, na.rm = TRUE))

muddyfoot_sub %>%
  group_by(date) %>%
  summarise(mean_time_diff = mean(time_diff, na.rm = TRUE),
            sd_time_diff = sd(time_diff, na.rm = TRUE),
            median_time_diff = median(time_diff, na.rm = TRUE),
            min_time_diff = min(time_diff, na.rm = TRUE),
            max_time_diff = max(time_diff, na.rm = TRUE))



## GENERAL DATAFRAME CLEANING ##

#receiver detection data. Do detected columns match used columns
muddyfoot_sub %>% 
  filter(RxDetected != RxUsed) %>% 
  nrow()

#remove unnecessary columns
muddyfoot_sub <- muddyfoot_sub %>% 
  dplyr::select(-Station, -HPEm, -TempData, -DepthData, -AccelData, 
         -Tag.Type, -Serial.Number, PIT.Number, Biologger.Number, -RxDetected, -nRxDetected)



#rename some columns
names(muddyfoot_sub)[names(muddyfoot_sub) == "RxUsed"] <- "recievers_used"
names(muddyfoot_sub)[names(muddyfoot_sub) == "nRxUsed"] <- "num_recievers_used"

#arrange columns more logically
muddyfoot_sub <- muddyfoot_sub %>% 
  dplyr::select(2,25, 20, 15:19, everything())

#retrospectively also delete Date column - as refers to dates that fish were tagged
muddyfoot_sub <- muddyfoot_sub %>% 
  dplyr::select(-Date)

#add long and lat columns back in the dataframe
coords <- st_coordinates(muddyfoot_sub$geometry)
muddyfoot_sub$Long <- coords[,1] 
muddyfoot_sub$Lat <- coords[, 2]

#Add columns for coordinaate converted to UTM

#Convert data to sf object with WGS84 CRS
data_sf <- st_as_sf(muddyfoot_sub, coords = c("Long", "Lat"), crs = 4326) 
#note that crs = 4326 is equivalent to WGS84

#Define the UTM CRS
utm_crs <- st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m +no_defs")
#zone 34 is the zone that Djacknabole is located within

#Transform coordinates to utm
data_sf_utm <- st_transform(data_sf, crs = utm_crs)

#Extract utm x and y coordinates
utm_coords <- st_coordinates(data_sf_utm$geometry)
muddyfoot_sub$Long_utm <- utm_coords[,1]
muddyfoot_sub$Lat_utm <- utm_coords[,2]

#Add ID column to work better with ctmm
muddyfoot_sub$Fish_ID = paste("F", muddyfoot_sub$individual_ID, sep = "")

#Save filtered data
saveRDS(muddyfoot_sub, file= paste0(save_data_sub_path, "muddyfoot_sub.rds"))



#--------------------------------------------------------------#
#### 1) Fish biometrics summary --------------------------------
#--------------------------------------------------------------#

## Fish standard lengths
muddyfoot %>% 
  filter(!individual_ID == 'Reference') %>% 
  dplyr::distinct(individual_ID, .keep_all = TRUE) %>%
  dplyr::group_by(Species, Treatment) %>%
  dplyr::summarise(mean = mean(Std_length, na.rm = T),
                   sd = sd(Std_length, na.rm = T),
                   max = max(Std_length, na.rm = T),
                   min = min(Std_length, na.rm = T),
                   num_unique_ID = length(unique(individual_ID)),
                   num_data_pts = sum(!is.na(Std_length)),
                   mean_datapts_ID = num_data_pts/num_unique_ID) %>%
  ungroup()


## Fish weights
muddyfoot %>% 
  filter(!individual_ID == 'Reference') %>% 
  dplyr::distinct(individual_ID, .keep_all = TRUE) %>%
  dplyr::group_by(Species,Treatment) %>%
  dplyr::summarise(mean = mean(Weight, na.rm = T),
                   sd = sd(Weight, na.rm = T),
                   max = max(Weight, na.rm = T),
                   min = min(Weight, na.rm = T),
                   num_unique_ID = length(unique(individual_ID)),
                   num_data_pts = sum(!is.na(Weight)),
                   mean_datapts_ID = num_data_pts/num_unique_ID) %>%
  ungroup()




