###
open_ctmm_path = "./data/ctmm_fits/lake_cow_roach_fits/"
save_ctmm_path = "./data/ctmm_fits/lake_cow_roach_fits/"

# List all RDS files in the folder
rds_files <- list.files(path = open_ctmm_path, pattern = "\\.rds$", full.names = TRUE)

# Read all RDS files into a list
rds_list <- lapply(rds_files, readRDS)

# Read all RDS files into a list
names(rds_list) <- basename(rds_files)
names(rds_list) <- sub("\\.rds$", "", names(rds_list))

# Print the names of the loaded RDS files
print(names(rds_list))

#remove pred filtered ids 
#these were rerun using ctmm.fit so we do not need to extract OUF
#OUF is already extracted. 

exclude_ids <- c("F59701", "F59738", "F59719", "F59709", "F59731", "F59687",
                 "F59745", "F59729", "F59683")

lake_cow_roach_ctmm_fits <- rds_list[!names(rds_list) %in% exclude_ids]

# Create a list of remaining elements
remaining_list <- rds_list[names(rds_list) %in% exclude_ids]

# #Isolate the files than need to be renamed
# # Create a new list to store renamed objects
# renamed_list <- list()
# 
# for (i in 1:8) {
#   # Extract the identity value
#   identity_value <- rds_list[[as.character(i)]]$`OUF anisotropic error`$ISO@info$identity
#   # Rename the object with its identity value
#   renamed_list[[identity_value]] <- rds_list[[as.character(i)]]
# }
# 
# # Function to save each list as an .rds file
# save_lists_as_rds <- function(fit_list, fit_name) {
#   saveRDS(fit_list, file = paste0(save_ctmm_path, fit_name, ".rds"))
# }
# 
# # Iterate over each element in renamed_list
# for (fit_name in names(renamed_list)) {
#   fit_list <- renamed_list[[fit_name]]
#   save_lists_as_rds(fit_list, fit_name)
# }
# 
# #save combined list of roach lake_cow ctmm fits with correct names
# lake_cow_roach_ctmm_fits <- rds_list
# saveRDS(lake_cow_roach_ctmm_fits, paste0(save_ctmm_path, "lake_cow_roach_ctmm_fits.rds"))

# ------------------ Extract OUF models ----------------------- #

lake_cow_roach_OUF_models <- list()

#Iterate over each object in lake_cow_roach_ctmm_fits and extract the 'OUF anisotropic error' models
lake_cow_roach_OUF_models <- lapply(lake_cow_roach_ctmm_fits, function(x) x[['OUF anisotropic error']])
#check it worked
summary(lake_cow_roach_OUF_models$F59710)

#rejoin pred filtered models
lake_cow_roach_OUF_models <- c(lake_cow_roach_OUF_models, remaining_list)[names(rds_list)]
#check that ID are in chonological order
names(lake_cow_roach_OUF_models)

#save
saveRDS(lake_cow_roach_OUF_models, paste0(save_ctmm_path, "lake_cow_roach_OUF_models.rds"))
