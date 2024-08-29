# _________________________ #
# HOME RANGE ANALYSIS - BT
# __________________________#

#LIBRARIES
library(ctmm)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)


#DIRECTORIES
ctmm_path = "./data/ctmm_fits/"
data_filter_path = "./data/tracks_filtered/"
telem_path = "./data/telem_obj/"
save_tables_path = "./data/tracks_filtered/sum_tables/" 
save_akde_path = "./data/akdes/" 
lake_polygon_path = "./data/lake_coords/"


#____________________#
#### 1. lake_BT ####
#____________________#

#>>> 1. Setup ####

### LOAD DATA ###

#lake_BT filtered dataframe
lake_BT_filt_data <-  readRDS(paste0(data_filter_path, "lake_BT_final_filt_data.rds"))

#telemetry objects
pike_lake_BT_tel = readRDS(paste0(telem_path, 'pike_lake_BT_tel.rds'))
perch_lake_BT_tel = readRDS(paste0(telem_path, 'perch_lake_BT_tel.rds'))
roach_lake_BT_tel = readRDS(paste0(telem_path, 'roach_lake_BT_tel.rds'))

#ctmm models
pike_lake_BT_ctmm_fits = readRDS(paste0(ctmm_path, "lake_BT_pike_fits/lake_BT_pike_OUF_models.rds"))
perch_lake_BT_ctmm_fits = readRDS(paste0(ctmm_path, "lake_BT_perch_fits/lake_BT_perch_OUF_models.rds"))
roach_lake_BT_ctmm_fits = readRDS(paste0(ctmm_path, "lake_BT_roach_fits/lake_BT_roach_OUF_models.rds")) 

#load in lake_BT polygon
lake_BT_polygon = sf::st_read(paste0(lake_polygon_path, "lake_BT_polygon.gpkg"))


#>>> 2. AKDE ESTIMATION ####

#Pike
pike_akdes <- list()
cl <- makeCluster(6)
doParallel::registerDoParallel(cl)
pike_akdes <- 
  foreach(i = 1:length(pike_lake_BT_tel), .packages = 'ctmm') %dopar% {
    akde_fit = akde(pike_lake_BT_tel[[i]],pike_lake_BT_ctmm_fits[[i]], weights = FALSE, sp = lake_BT_polygon)
    saveRDS(akde_fit, file = paste0(save_akde_path, "lake_BT_pike_akdes/", names(pike_lake_BT_tel)[i], "_akde.rds"))
    akde_fit
  }


#add ID to lists
names(pike_akdes ) <- names(pike_lake_BT_tel)
summary(pike_akdes$F59880)

#Perch
perch_akdes <- list()
cl <- makeCluster(20)
doParallel::registerDoParallel(cl)
perch_akdes <- 
  foreach(i = 1:length(perch_lake_BT_tel), .packages = 'ctmm') %dopar% {
    akde_fit = akde(perch_lake_BT_tel[[i]], 
                    perch_lake_BT_ctmm_fits[[i]], 
                    weights = FALSE, 
                    sp = lake_BT_polygon)
    saveRDS(akde_fit, 
            file = paste0(save_akde_path, 
                          "lake_BT_perch_akdes/", 
                          names(perch_lake_BT_tel)[i], 
                          "_akde.rds"))
    akde_fit
  }

stopCluster(cl)

#add ID to lists
names(perch_akdes ) <- names(perch_lake_BT_tel)
summary(perch_akdes$F59682)  

#save
saveRDS(perch_akdes, paste0(save_akde_path, "lake_BT_perch_akdes/perch_akdes_list.rds"))

#Roach
roach_akdes <- list()
cl <- makeCluster(20)
doParallel::registerDoParallel(cl)
roach_akdes <- 
  foreach(i = 1:length(roach_lake_BT_tel), .packages = 'ctmm') %dopar% {
    akde_fit = akde(roach_lake_BT_tel[[i]], 
                    roach_lake_BT_ctmm_fits[[i]], 
                    weights = FALSE, 
                    sp = lake_BT_polygon)
    saveRDS(akde_fit, 
            file = paste0(save_akde_path, 
                          "lake_BT_roach_akdes/", 
                          names(roach_lake_BT_tel)[i], 
                          "_akde.rds"))
    akde_fit
  }

stopCluster(cl)

#add ID to lists
names(roach_akdes ) <- names(roach_lake_BT_tel)

#save
saveRDS(roach_akdes, paste0(save_akde_path, "lake_BT_roach_akdes/roach_akdes_list.rds"))
