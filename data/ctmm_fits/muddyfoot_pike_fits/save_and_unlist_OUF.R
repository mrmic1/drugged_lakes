###
open_ctmm_path = "./data/ctmm_fits/muddyfoot_pike_fits/"
save_ctmm_path = "./data/ctmm_fits/muddyfoot_pike_fits/"

# ------------------ Extract OUF models ----------------------- #

muddyfoot_pike_ctmm_fits <- readRDS(paste0(open_ctmm_path, "muddyfoot_pike_ctmm_fits.rds"))

muddyfoot_pike_OUF_models <- list()

#Iterate over each object in muddyfoot_pike_ctmm_fits and extract the 'OUF anisotropic error' models
muddyfoot_pike_OUF_models <- lapply(muddyfoot_pike_ctmm_fits, function(x) x[['OUF anisotropic error']])

# Set the names of the new list to the identities of the original objects (if needed)
names(muddyfoot_pike_OUF_models) <- names(muddyfoot_pike_ctmm_fits)

#Check model outputs
for (i in 1:length(muddyfoot_pike_OUF_models)) {
  element_name <- names(muddyfoot_pike_OUF_models)[i]
  cat("Summary for", element_name, ":\n")
  print(summary(muddyfoot_pike_OUF_models[[i]]))
  cat("\n")
}


saveRDS(muddyfoot_pike_OUF_models, paste0(save_ctmm_path, "muddyfoot_pike_OUF_models.rds"))
