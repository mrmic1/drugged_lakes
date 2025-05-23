#-----------------------------------#
# Combine cow paradise perch ctmms
#-----------------------------------#
  
  ctmm_path = "./data/ctmm_fits/lake_cow_perch_fits/"
  
  # List all RDS files in the folder
  rds_files <- list.files(path = ctmm_path, pattern = "\\.rds$", full.names = TRUE)
  
  # Read all RDS files into a list
  rds_list <- lapply(rds_files, readRDS)
  
  # Read all RDS files into a list
  names(rds_list) <- basename(rds_files)
  names(rds_list) <- sub("\\.rds$", "", names(rds_list))
  
  #Need to remove last element which contain a list of the ctmms
  names(rds_list)
  rds_list <- rds_list[-41]
  
  
  
  ## Now I need to extract the OUF models ##
  lake_cow_perch_OUF_models <- list()
  
  #Iterate over each object in rds_list and extract the 'OUF anisotropic error' models
  lake_cow_perch_OUF_models <- lapply(rds_list, function(x) x[['OUF anisotropic error']])
  
  # Set the names of the new list to the identities of the original objects (if needed)
  names(lake_cow_perch_OUF_models) <- names(rds_list)
  print(names(lake_cow_perch_OUF_models))
  
  #Check model outputs
  for (i in 1:length(lake_cow_perch_OUF_models)) {
    element_name <- names(lake_cow_perch_OUF_models)[i]
    cat("Summary for", element_name, ":\n")
    print(summary(lake_cow_perch_OUF_models[[i]]))
    cat("\n")
  }
  

  names(lake_cow_perch_OUF_models)
  
  
  #save list
  saveRDS(lake_cow_perch_OUF_models, paste0(ctmm_path, "lake_cow_perch_OUF_models.rds"))
  