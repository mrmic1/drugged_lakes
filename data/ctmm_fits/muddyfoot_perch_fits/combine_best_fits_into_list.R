# Set the folder path
folder_path <- "./data/ctmm_fits/muddyfoot_perch_fits/"

# Get all .rds files that match the pattern
file_list <- list.files(
  path = folder_path,
  pattern = "_ctmm_best_fit(_refitted)?\\.rds$",
  full.names = TRUE
)

# Read all files into a list
muddyfoot_perch_best_ctmm_model_fits <- lapply(file_list, readRDS)




# Extract ID numbers from filenames for naming the list elements
id_names <- basename(file_list)
id_names <- gsub("_ctmm_best_fit(_refitted)?\\.rds$", "", id_names)

# Name the list elements with the ID numbers
names(muddyfoot_perch_best_ctmm_model_fits) <- id_names

# Save the combined list
saveRDS(muddyfoot_perch_best_ctmm_model_fits, 
        file = paste0(folder_path, "muddyfoot_perch_best_ctmm_model_fits.rds"))

# Print summary
cat("Combined", length(muddyfoot_perch_best_ctmm_model_fits), "ctmm model fits\n")
cat("IDs included:", paste(names(muddyfoot_perch_best_ctmm_model_fits), collapse = ", "), "\n")
