#Check select fits
#Seperate fits into seperate files if needed

muddyfoot_perch_select_fits <- readRDS("~/Desktop/Git_Repo/pharma-LOF-mac/data/ctmm_fits/muddyfoot_perch_fits/muddyfoot_perch_select_fits.rds")

#Check model outputs
for (i in 1:length(muddyfoot_perch_select_fits)) {
  element_name <- names(muddyfoot_perch_select_fits)[i]
  cat("Summary for", element_name, ":\n")
  print(summary(muddyfoot_perch_select_fits[[i]]))
  cat("\n")
}

# Assuming muddyfoot_perch_select_fits is a list of lists, where each sublist contains fits
# We will use lapply to iterate over each element and generate a summary

summaries <- lapply(muddyfoot_perch_select_fits, function(fit) {
  summary(fit)$name
})

# Display summaries
summaries

#14/30 = "OU ansiotropic error"

save_ctmm_path = "./data/ctmm_fits/muddyfoot_perch_fits/"

# Function to save each list as an .rds file
save_lists_as_rds <- function(fit_list, fit_name) {
  saveRDS(fit_list, file = paste0(save_ctmm_path, fit_name, ".rds"))
}

# Iterate over each element in muddyfoot_perch_select_fits
for (fit_name in names(muddyfoot_perch_select_fits)) {
  fit_list <- muddyfoot_perch_select_fits[[fit_name]]
  save_lists_as_rds(fit_list, fit_name)
}


#----------------------------------------------------------------------------#

# Set the directory where your .rds files are stored
setwd("C:/Users/MarcusMichelangeli/Documents/Git_repo/pharma-LOF-griffith/pharma-LOF/data/ctmm_fits/muddyfoot_perch_fits")

# List all files in the directory that have "_OUF" in the name and end with .rds
files <- list.files(pattern = "_OUF.*\\.rds$")

# Create an empty list to store your models
muddyfoot_perch_ctmm_fits <- list()

# Loop through each file and read the .rds into the list
for (file in files) {
  # Extract the identifier from the filename (assuming it is the part before the first underscore)
  identifier <- sub("_.*", "", file)
  
  # Read the .rds file and store it in the list with the identifier as the name
  muddyfoot_perch_ctmm_fits[[identifier]] <- readRDS(file)
}

# Now, muddyfoot_perch_ctmm_fits contains all the models with their identifiers as names

## Now I need to extract the OUF models ##
muddyfoot_perch_OUF_models_2 <- list()

#Iterate over each object in muddyfoot_pike_ctmm_fits and extract the 'OUF anisotropic error' models
muddyfoot_perch_OUF_models_2 <- lapply(muddyfoot_perch_ctmm_fits, function(x) x[['OUF anisotropic error']])

# Set the names of the new list to the identities of the original objects (if needed)
names(muddyfoot_perch_OUF_models_2) <- names(muddyfoot_perch_ctmm_fits)

muddyfoot_perch_OUF_models_2 = muddyfoot_perch_OUF_models_2[1:14]

#Check model outputs
for (i in 1:length(muddyfoot_perch_OUF_models_2)) {
  element_name <- names(muddyfoot_perch_OUF_models_2)[i]
  cat("Summary for", element_name, ":\n")
  print(summary(muddyfoot_perch_OUF_models_2[[i]]))
  cat("\n")
}

#save each file

save_ctmm_path = "./data/ctmm_fits/muddyfoot_perch_fits/"

# Function to save each list as an .rds file
save_lists_as_rds <- function(fit_list, fit_name) {
  saveRDS(fit_list, file = paste0(save_ctmm_path, fit_name, ".rds"))
}

# Iterate over each element in muddyfoot_perch_OUF_models_2
for (fit_name in names(muddyfoot_perch_OUF_models_2)) {
  fit_list <- muddyfoot_perch_OUF_models_2[[fit_name]]
  save_lists_as_rds(fit_list, fit_name)
}


#------------------------------------------------------------------------------------------#

#Combine all OUF models into a single list

# Set the working directory to the folder containing the .rds files
setwd("C:/Users/MarcusMichelangeli/Documents/Git_repo/pharma-LOF-griffith/pharma-LOF/data/ctmm_fits/muddyfoot_perch_fits")

# List all .rds files in the directory
files <- list.files(pattern = "\\.rds$")

# Initialize an empty list to store the models
muddyfoot_perch_OUF_models <- list()

# Loop through each file, read it, and add it to the list
for (file in files) {
  # Extract the identifier (e.g., F59682) from the file name, remove any .rds or suffixes
  identifier <- sub("\\.rds$", "", sub("_.*", "", file))
  
  # Read the .rds file
  model <- readRDS(file)
  
  # Add the model to the list with the identifier as the name
  muddyfoot_perch_OUF_models[[identifier]] <- model
}


saveRDS(muddyfoot_perch_OUF_models, file = "muddyfoot_perch_OUF_models.rds")

