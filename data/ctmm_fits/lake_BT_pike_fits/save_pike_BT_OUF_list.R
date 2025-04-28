#-----------------------------#
# Combine BT pike ctmms
#-----------------------------#

ctmm_path = "./data/ctmm_fits/lake_BT_pike_fits/"

# List all RDS files in the folder
rds_files <- list.files(path = ctmm_path, pattern = "\\.rds$", full.names = TRUE)

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

exclude_ids <- c("F59886", "F59892")

BT_pike_ctmm_fits <- rds_list[!names(rds_list) %in% exclude_ids]

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
# #save combined list of roach muddyfoot ctmm fits with correct names
# muddyfoot_roach_ctmm_fits <- rds_list
# saveRDS(muddyfoot_roach_ctmm_fits, paste0(save_ctmm_path, "muddyfoot_roach_ctmm_fits.rds"))

# ------------------ Extract OUF models ----------------------- #

BT_pike_OUF_models <- list()

#Iterate over each object in muddyfoot_roach_ctmm_fits and extract the 'OUF anisotropic error' models
BT_pike_OUF_models <- lapply(BT_pike_ctmm_fits, function(x) x[['OUF anisotropic error']])
#check it worked
summary(BT_pike_OUF_models$F59887)

#rejoin pred filtered models
BT_pike_OUF_models <- c(BT_pike_OUF_models, remaining_list)[names(rds_list)]
#check that ID are in chonological order
names(BT_pike_OUF_models)
summary(BT_pike_OUF_models$F59886)


#save
saveRDS(BT_pike_OUF_models, paste0(ctmm_path, "BT_pike_OUF_models.rds"))
