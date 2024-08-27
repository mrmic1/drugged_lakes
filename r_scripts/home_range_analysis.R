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
save_akde_path = "./data/akdes/muddyfoot/" 
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

muddyfoot_filt_data %>% 
  filter(Species == 'Pike') %>% 
  select(individual_id) %>% 
  distinct()


#create a sample dataset to test parallel processing
test_dat = 
  muddyfoot_filt_data %>% 
  filter(individual_id == 'F59880' | individual_id == 'F59881') %>% 
  filter(date == '2022/10/01')
  
test_tel  <- as.telemetry(test_dat, 
                          timezone = "Europe/Stockholm", 
                          timeformat="%Y-%m-%d %H:%M:%S", 
                          projection= NULL,
                          datum="WGS84")
test_fits <- pike_muddyfoot_ctmm_fits[1:2]


#Pike
test_akdes <- list()
cl <- makeCluster(2)
doParallel::registerDoParallel(cl)
test_akdes <- 
  foreach(i = 1:length(test_tel), .packages = 'ctmm') %dopar% {
    akde_fit = akde(test_tel[[i]],test_fits[[i]], weights = TRUE, sp = muddyfoot_polygon)
    saveRDS(akde_fit, file = paste0(save_akde_path, "muddyfoot_pike_akdes/", names(test_tel)[i], "_akde.rds"))
    akde_fit
  }

    
#add ID to lists
names(test_akdes ) <- names(test_tel)
summary(test_akdes$F59880)


#Seperate telemetry objects based on irregular sampling



# View the summary
print(summary_data)

  



