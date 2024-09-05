#--------------------------------------------#
# HOME RANGE ANALYSIS - MUDDYFOOT
# -------------------------------------------#

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
akde_path = "./data/akdes/" 
lake_polygon_path = "./data/lake_coords/"


#Need to make akde_cgs
#population akdes
#change projection of ctmms objects
#Add LFS to github repository
#Create raster maps

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

#-------------------------------------------------------------------------------------------------------

#--------------------------#                 
#>>> 2. AKDE ESTIMATION ####
#--------------------------#

### Pike HR analysis ###

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


# calculate AKDES on a consistent grid
# first estimate aKDE for one individual whose HR is representative enough to use as a reference for the grid size
# load pike akde list
pike_akdes_list <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_akdes_list.rds"))
names(pike_akdes_list)
akde_ref <- pike_akdes_list[[1]]

pike_akdes_cg <- list()
cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
pike_akdes_cg <- foreach(i = 1:length(pike_muddyfoot_tel), .packages = 'ctmm') %dopar% {
  akde_fit_cg = akde(
    pike_muddyfoot_tel[[i]],
    pike_muddyfoot_ctmm_fits[[i]], 
    weights = FALSE,
    grid = list(dr = akde_ref$dr,
                align.to.origin = T))
  saveRDS(akde_fit_cg, file = paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/", names(pike_muddyfoot_tel)[i], "_akde_cg.rds"))
  akde_fit_cg
}

stopCluster(cl)

#add ID to lists
names(pike_akdes_cg) <- names(pike_muddyfoot_tel)
summary(pike_akdes_cg)  

#save
saveRDS(pike_akdes_cg, paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/pike_akdes_cg_list.rds"))

#------------------------------------------------#

### Perch HR analysis ###

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


# calculate AKDES on a consistent grid
# first estimate aKDE for one individual whose HR is representative enough to use as a reference for the grid size
# load perch akde list
perch_akdes_list <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/perch_akdes_list.rds"))
names(perch_akdes_list)
akde_ref <- perch_akdes_list[[9]]

perch_akdes_cg <- list()
cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
perch_akdes_cg <- foreach(i = 1:length(perch_muddyfoot_tel), .packages = 'ctmm') %dopar% {
  akde_fit_cg = akde(
    perch_muddyfoot_tel[[i]],
    perch_muddyfoot_ctmm_fits[[i]], 
    weights = FALSE,
    grid = list(dr = akde_ref$dr,
                align.to.origin = T))
  saveRDS(akde_fit_cg, file = paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/", names(perch_muddyfoot_tel)[i], "_akde_cg.rds"))
  akde_fit_cg
}

stopCluster(cl)

#add ID to lists
names(perch_akdes_cg) <- names(perch_muddyfoot_tel)
summary(perch_akdes_cg)  

#save
#saveRDS(perch_akdes_cg, paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))

#----------------------------------------------------------#

### Roach HR analysis ###

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

# calculate AKDES on a consistent grid
# first estimate aKDE for one individual whose HR is representative enough to use as a reference for the grid size
# load perch akde list
roach_akdes_list <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/roach_akdes_list.rds"))
names(roach_akdes_list)
akde_ref <- roach_akdes_list[[7]]

roach_akdes_cg <- list()
cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
roach_akdes_cg <- foreach(i = 1:length(roach_muddyfoot_tel), .packages = 'ctmm') %dopar% {
  akde_fit_cg = akde(
    roach_muddyfoot_tel[[i]],
    roach_muddyfoot_ctmm_fits[[i]], 
    weights = FALSE,
    grid = list(dr = akde_ref$dr,
                align.to.origin = T))
  saveRDS(akde_fit_cg, file = paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/", names(roach_muddyfoot_tel)[i], "_akde_cg.rds"))
  akde_fit_cg
}

stopCluster(cl)

#add ID to lists
names(roach_akdes_cg) <- names(roach_muddyfoot_tel)
summary(roach_akdes_cg)  

#save
#saveRDS(roach_akdes_cg, paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))

#------------------------------------------------------------------------------------------------#

#-------------------------------------#
#>>> 3. POPULATION AKDE ESTIMATION ####
#-------------------------------------#

#make sure muddyfoot polygon is loaded

### Pike population akdes ###

# Separating into 'control' and 'mix'
#telemetry objects
pike_control_tel <- pike_muddyfoot_tel[1:3]
pike_mix_tel <- pike_muddyfoot_tel[4:6]

#ctmms
pike_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/pike_akdes_cg_list.rds"))
pike_control_akdes <- pike_akdes_cg_list[1:3]
pike_mix_akdes <- pike_akdes_cg_list[4:6]

#calculate population-level autocorrelated kernel desnity home range estimates
pike_control_PKDE <- pkde(pike_control_tel,
                          pike_control_akdes, 
                          sp = muddyfoot_polygon)

saveRDS(pike_control_PKDE, paste0(akde_path, "muddyfoot_pike_akdes/population_akde/pike_control_PKDE.rds"))

pike_mix_PKDE <- pkde(pike_mix_tel,
                      pike_mix_akdes, 
                      sp = muddyfoot_polygon)

saveRDS(pike_mix_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/population_akde/pike_mix_PKDE.rds"))

### Perch population akdes ###

# Separating into 'control' and 'mix'
#telemetry objects
perch_control_tel <- perch_muddyfoot_tel[1:3]
perch_mix_tel <- perch_muddyfoot_tel[4:6]

#ctmms
perch_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))
perch_control_akdes <- perch_akdes_cg_list[1:3]
perch_mix_akdes <- perch_akdes_cg_list[4:6]

#calculate population-level autocorrelated kernel desnity home range estimates
perch_control_PKDE <- pkde(perch_control_tel,
                          perch_control_akdes, 
                          sp = muddyfoot_polygon)

saveRDS(perch_control_PKDE, paste0(akde_path, "muddyfoot_perch_akdes/population_akde/perch_control_PKDE.rds"))

perch_mix_PKDE <- pkde(perch_mix_tel,
                      perch_mix_akdes, 
                      sp = muddyfoot_polygon)

saveRDS(perch_mix_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/population_akde/perch_mix_PKDE.rds"))

### roach population akdes ###

# Separating into 'control' and 'mix'
#telemetry objects
roach_control_tel <- roach_muddyfoot_tel[1:3]
roach_mix_tel <- roach_muddyfoot_tel[4:6]

#ctmms
roach_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))
roach_control_akdes <- roach_akdes_cg_list[1:3]
roach_mix_akdes <- roach_akdes_cg_list[4:6]

#calculate population-level autocorrelated kernel desnity home range estimates
roach_control_PKDE <- pkde(roach_control_tel,
                          roach_control_akdes, 
                          sp = muddyfoot_polygon)

saveRDS(roach_control_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/population_akde/roach_control_PKDE.rds"))

roach_mix_PKDE <- pkde(roach_mix_tel,
                      roach_mix_akdes, 
                      sp = muddyfoot_polygon)

#saveRDS(roach_mix_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/population_akde/roach_mix_PKDE.rds"))


