# _________________________ #
# HOME RANGE ANALYSIS
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
#### 1. MUDDYFOOT ####
#____________________#

#>>> 1. Setup ####

### LOAD DATA ###

#muddyfoot filtered dataframe
muddyfoot_filt_data <-  readRDS(paste0(data_filter_path, "muddyfoot_final_filt_data.rds"))

#telemetry objects
pike_muddyfoot_tel = readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel = readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel = readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

#ctmm models
pike_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds")) 

#load in muddyfoot polygon
muddyfoot_polygon = sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))

                 

#>>> 2. AKDE ESTIMATION ####

#Pike
pike_akdes <- list()
cl <- makeCluster(6)
doParallel::registerDoParallel(cl)
pike_akdes <- 
  foreach(i = 1:length(pike_muddyfoot_tel), .packages = 'ctmm') %dopar% {
    akde_fit = akde(pike_muddyfoot_tel[[i]],pike_muddyfoot_ctmm_fits[[i]], weights = FALSE, sp = muddyfoot_polygon)
    saveRDS(akde_fit, file = paste0(save_akde_path, "muddyfoot_pike_akdes/", names(pike_muddyfoot_tel)[i], "_akde.rds"))
    akde_fit
  }


#add ID to lists
names(pike_akdes ) <- names(pike_muddyfoot_tel)
summary(pike_akdes$F59880)

#Perch
perch_akdes <- list()
cl <- makeCluster(20)
doParallel::registerDoParallel(cl)
perch_akdes <- 
  foreach(i = 1:length(perch_muddyfoot_tel), .packages = 'ctmm') %dopar% {
    akde_fit = akde(perch_muddyfoot_tel[[i]], 
                    perch_muddyfoot_ctmm_fits[[i]], 
                    weights = FALSE, 
                    sp = muddyfoot_polygon)
    saveRDS(akde_fit, 
            file = paste0(save_akde_path, 
                          "muddyfoot_perch_akdes/", 
                          names(perch_muddyfoot_tel)[i], 
                          "_akde.rds"))
    akde_fit
  }

stopCluster(cl)

#add ID to lists
names(perch_akdes ) <- names(perch_muddyfoot_tel)
summary(perch_akdes$F59682)  

#save
saveRDS(perch_akdes, paste0(save_akde_path, "muddyfoot_perch_akdes/perch_akdes_list.rds"))

#Roach
roach_akdes <- list()
cl <- makeCluster(20)
doParallel::registerDoParallel(cl)
roach_akdes <- 
  foreach(i = 1:length(roach_muddyfoot_tel), .packages = 'ctmm') %dopar% {
    akde_fit = akde(roach_muddyfoot_tel[[i]], 
                    roach_muddyfoot_ctmm_fits[[i]], 
                    weights = FALSE, 
                    sp = muddyfoot_polygon)
    saveRDS(akde_fit, 
            file = paste0(save_akde_path, 
                          "muddyfoot_roach_akdes/", 
                          names(roach_muddyfoot_tel)[i], 
                          "_akde.rds"))
    akde_fit
  }

stopCluster(cl)

#add ID to lists
names(roach_akdes ) <- names(roach_muddyfoot_tel)

#save
saveRDS(roach_akdes, paste0(save_akde_path, "muddyfoot_roach_akdes/roach_akdes_list.rds"))

#try weighted akde - takes much longer, potentialy at least 24 hrs to run. 
#Roach
roach_akdes_W <- list()
cl <- makeCluster(20)
doParallel::registerDoParallel(cl)
roach_akdes_W <- 
  foreach(i = 1:length(roach_muddyfoot_tel), .packages = 'ctmm') %dopar% {
    akde_fit = akde(roach_muddyfoot_tel[[i]], 
                    roach_muddyfoot_ctmm_fits[[i]], 
                    weights = TRUE, 
                    sp = muddyfoot_polygon)
    saveRDS(akde_fit, 
            file = paste0(save_akde_path, 
                          "muddyfoot_roach_akdes/akdeW/", 
                          names(roach_muddyfoot_tel)[i], 
                          "_akde.rds"))
    akde_fit
  }

stopCluster(cl)