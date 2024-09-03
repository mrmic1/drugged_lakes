#----------------------------------------#
# Comparing akdes #######################
#---------------------------------------#

#LIBRARIES
library(ctmm)

#DIRECTORIES
pike_akde_path = "./data/akdes/muddyfoot_pike_akdes/" 
pike_ctmm_path = "./data/ctmm_fits/muddyfoot_pike_fits"
pike_telem_path = "./data/telem_obj/"

#akdes without consistent grid
pike_akdes_list = readRDS(paste0(pike_akde_path, "pike_akdes_list.rds"))

#akdes on consistent grid
pike_akdes_cg_list = readRDS(paste0(pike_akde_path, "akdes_cg/pike_akdes_cg_list.rds"))


#Compare summaries
summary(pike_akdes_list$F59885)
#estimate = 433.51

summary(pike_akdes_cg_list$F59885)
#estimate = 433.55

#Almost identical

#For plotting
COL <- color(pike_akdes_list, by='individual')
#code does not work
COL_cg <- color(pike_akdes_cg_list, by='individual')
#code does work

# plot AKDEs with consistent grids
plot(pike_akdes_cg_list,
     col.UD=COL_cg,
     col.level=COL_cg,
     col.grid=NA,
     level=NA,
     main="Muddyfoot pike AKDEs")

#Meta-analysis results
meta(pike_akdes_list,col='black',sort=F, verbose = T, level.UD = 0.95)
meta(pike_akdes_cg_list,col='black',sort=F, verbose = T, level.UD = 0.95)
#basically the same


#Separating into 'control' and 'mix'
pike_akde_control <- pike_akdes_list[1:3]
pike_akde_mix <- pike_akdes_list[4:6]
pike_akde_total <- list(pike_akde_control = pike_akde_control, pike_akde_mix = pike_akde_mix)
meta(pike_akde_total,col='black',sort=F, verbose = T, level.UD = 0.95)

pike_akde_control_cg <- pike_akdes_cg_list[1:3]
pike_akde_mix_cg <- pike_akdes_cg_list[4:6]
pike_akde_total_cg <- list(pike_akde_control_cg = pike_akde_control_cg, pike_akde_mix_cg = pike_akde_mix_cg)
meta(pike_akde_total_cg,col='black',sort=F, verbose = T, level.UD = 0.95)
#Identical

#What about for plotting an encounter distribution
over <- ctmm::overlap(pike_akdes_list[c("F59880", "F59884")])
OVER <- overlap(pike_akdes_list)
#can't do because of inconsistent grid resolution
over_cg <- ctmm::overlap(pike_akdes_cg_list[c("F59880", "F59884")])
OVER_cg <- overlap(pike_akdes_cg_list)
#this works
over_cg
OVER_cg

#Estimating the CDE
CDE <- encounter(pike_akdes_list[c("F59880", "F59884")])
#inconsistent grid resolution
CDE_cg <- encounter(pike_akdes_cg_list[c("F59880", "F59884")])
summary(CDE_cg)
CDE_cg

#Plot the data and HR estimates
ctmm::plot(pike_muddyfoot_tel[c("F59880", "F59884")],
     UD = pike_akdes_cg_list[c("F59880", "F59884")],
     col = NA,
     col.DF= c("#f4a261", "#2a9d8f"),
     col.grid = NA)
#this works


#Visualise the CDE
ctmm::plot(pike_muddyfoot_tel[c("F59880", "F59884")],
           UD = CDE_cg,
           col.UD="red",
           #col=c("#e76f51", "#264653"),
           error = FALSE,
           col.grid = NA)
#Error: unable to find an inherited method for function ‘extent’ for signature ‘x = "overlap"’
#the same error message pops up when running the same code in ctmm_interactions tutorial.


