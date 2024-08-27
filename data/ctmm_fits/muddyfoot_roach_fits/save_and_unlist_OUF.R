###
open_ctmm_path = "./data/ctmm_fits/muddyfoot_roach_fits/"
save_ctmm_path = "./data/ctmm_fits/muddyfoot_roach_fits/"

# List all RDS files in the folder
rds_files <- list.files(path = open_ctmm_path, pattern = "\\.rds$", full.names = TRUE)

# Read all RDS files into a list
rds_list <- lapply(rds_files, readRDS)

# Read all RDS files into a list
names(rds_list) <- basename(rds_files)
names(rds_list) <- sub("\\.rds$", "", names(rds_list))

# Print the names of the loaded RDS files
print(names(rds_list))

#remove the the last list
rds_list <- rds_list[1:29]
names(rds_list)

#Isolate the files than need to be renamed
# Create a new list to store renamed objects
renamed_list <- list()

for (i in 1:8) {
  # Extract the identity value
  identity_value <- rds_list[[as.character(i)]]$`OUF anisotropic error`$ISO@info$identity
  # Rename the object with its identity value
  renamed_list[[identity_value]] <- rds_list[[as.character(i)]]
}

# Function to save each list as an .rds file
save_lists_as_rds <- function(fit_list, fit_name) {
  saveRDS(fit_list, file = paste0(save_ctmm_path, fit_name, ".rds"))
}

# Iterate over each element in renamed_list
for (fit_name in names(renamed_list)) {
  fit_list <- renamed_list[[fit_name]]
  save_lists_as_rds(fit_list, fit_name)
}

#save combined list of roach muddyfoot ctmm fits with correct names
muddyfoot_roach_ctmm_fits <- rds_list
saveRDS(muddyfoot_roach_ctmm_fits, paste0(save_ctmm_path, "muddyfoot_roach_ctmm_fits.rds"))

# ------------------ Extract OUF models ----------------------- #

muddyfoot_roach_OUF_models <- list()

#Iterate over each object in muddyfoot_roach_ctmm_fits and extract the 'OUF anisotropic error' models
muddyfoot_roach_OUF_models <- lapply(muddyfoot_roach_ctmm_fits, function(x) x[['OUF anisotropic error']])

# Set the names of the new list to the identities of the original objects (if needed)
names(muddyfoot_roach_OUF_models) <- names(muddyfoot_roach_ctmm_fits)
saveRDS(muddyfoot_roach_OUF_models, paste0(save_ctmm_path, "muddyfoot_roach_OUF_models.rds"))
