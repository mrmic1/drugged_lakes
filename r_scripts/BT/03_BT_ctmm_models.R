#==============================================================================
# Lake BT Detection Data - Continuous-Time Movement Modeling (CTMM)
#==============================================================================
# Purpose: Fit continuous-time movement models to individual fish trajectories
# Author: Marcus Michelangeli
#
# Models fitted:
#   - Ornstein-Uhlenbeck Foraging (OUF) models with location error
#   - Separate models for each individual fish
#   - Species-specific processing (Pike, Perch, Roach)
#
# Input files:
#   - ./data/tracks_filtered/lake_BT/03_lake_BT_sub.rds
#
# Output files:
#   - ./data/telem_obj/BT/[species]_lake_BT_tel_thinned.rds (3 files)
#   - ./data/ctmm_fits/lake_BT_[species]_fits/[individual]_ctmm_fit.rds
#   - ./data/ctmm_fits/lake_BT_[species]_fits/lake_BT_[species]_ctmm_fits (3 files)
#==============================================================================

# Load required libraries ---------------------------------------------------
library(data.table)
library(tidyverse)
library(ctmm)
library(sf)
library(parallel)
library(foreach)
library(doParallel)

# Set timezone globally -----------------------------------------------------
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths ---------------------------------------------------------
filtered_data_path <- "./data/tracks_filtered/lake_BT/"
save_ctmm_path <- "./data/ctmm_fits/"
save_telem_path <- "./data/telem_obj/"

#==============================================================================
# 1. LOAD AND PREPARE DATA
#==============================================================================

# Load filtered tracking data -----------------------------------------------
lake_BT_sub <- readRDS(paste0(filtered_data_path, '03_lake_BT_sub.rds'))
message("Loaded ", nrow(lake_BT_sub), " detections for ", 
        n_distinct(lake_BT_sub$individual_ID), " individuals")

# Convert to Movebank format for ctmm package ------------------------------
lake_BT_movebank <- with(
  lake_BT_sub,
  data.frame(
    "timestamp" = timestamp,
    "location.long" = Long,
    "location.lat" = Lat,
    "GPS.HDOP" = HDOP,
    "individual-local-identifier" = individual_ID,
    "Species" = Species,
    "Weight" = Weight,
    "Total_length" = Total_length,
    "Std_length" = Std_length,
    "Treatment" = Treatment,
    "Date" = Date,
    "Exp_Stage" = Exp_Stage,
    "Time_Of_Day" = Time_Of_Day,
    "found_alive" = found_alive,
    "known_predated" = known_predated
  )
)

# Free up memory ------------------------------------------------------------
rm(lake_BT_sub)
gc()

# Convert to telemetry object -----------------------------------------------
message("\nConverting to ctmm telemetry object...")
lake_BT_tels <- as.telemetry(
  lake_BT_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,  # Will be set to geometric median automatically
  datum = "WGS84",
  keep = c("Species", "Weight", "Total_length", "Std_length", "Treatment",
           "Date", "Exp_Stage", "Time_Of_Day", "found_alive", "known_predated")
)

message("Telemetry object created with ", length(lake_BT_tels), " individuals")
message("Projection: ", projection(lake_BT_tels[[1]]))
message("Timezone: ", tz(lake_BT_tels[[1]]$timestamp))

#==============================================================================
# 2. ORGANIZE INDIVIDUALS BY SPECIES
#==============================================================================

#==============================================================================
# 2. ORGANIZE INDIVIDUALS BY SPECIES
#==============================================================================

# Verify individual order ---------------------------------------------------
species_order <- lake_BT_movebank %>%
  select(Species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

message("\n=== Species Distribution (Before Filtering) ===")
print(table(species_order$Species, useNA = "ifany"))

# Filter out reference tag and any individuals with NA species --------------
species_order_filtered <- species_order %>%
  filter(!is.na(Species))

message("\n=== Species Distribution (After Filtering) ===")
print(table(species_order_filtered$Species))

# Add index column to match telemetry list order ----------------------------
species_order_filtered$telem_index <- 1:nrow(species_order_filtered)

# Display the mapping table for verification --------------------------------
message("\n=== Individual ID to Species Mapping (Filtered) ===")
print(species_order_filtered)

# Automatically extract indices for each species ----------------------------
pike_indices <- species_order_filtered$telem_index[species_order_filtered$Species == "Northern Pike"]
perch_indices <- species_order_filtered$telem_index[species_order_filtered$Species == "Perch"]
roach_indices <- species_order_filtered$telem_index[species_order_filtered$Species == "Roach"]

# Split telemetry objects by species ----------------------------------------
pike_lake_BT_tel <- lake_BT_tels[pike_indices]
perch_lake_BT_tel <- lake_BT_tels[perch_indices]
roach_lake_BT_tel <- lake_BT_tels[roach_indices]

message("\n=== Species Telemetry Objects Created ===")
message("Pike: ", length(pike_lake_BT_tel), " individuals (indices: ", 
        paste(pike_indices, collapse = ", "), ")")
message("Perch: ", length(perch_lake_BT_tel), " individuals (indices: ", 
        paste(perch_indices, collapse = ", "), ")")
message("Roach: ", length(roach_lake_BT_tel), " individuals (indices: ", 
        paste(roach_indices, collapse = ", "), ")")

# Verify totals match --------------------------------------------------------
total_individuals <- length(pike_lake_BT_tel) + length(perch_lake_BT_tel) + length(roach_lake_BT_tel)
message("\nTotal individuals across all species: ", total_individuals)
message("Expected total (excluding reference tag): ", length(lake_BT_tels))

if(total_individuals != length(lake_BT_tels)) {
  warning("Mismatch in total individuals! Check species assignment.")
}

# Save species-specific telemetry objects -----------------------------------
saveRDS(pike_lake_BT_tel, paste0(save_telem_path, "BT/pike_lake_BT_tel_thinned.rds"))
saveRDS(perch_lake_BT_tel, paste0(save_telem_path, "BT/perch_lake_BT_tel_thinned.rds"))
saveRDS(roach_lake_BT_tel, paste0(save_telem_path, "BT/roach_lake_BT_tel_thinned.rds"))

message("\nTelemetry objects saved to:", save_telem_path, "BT/")

#==============================================================================
# 3. HELPER FUNCTIONS
#==============================================================================

# Function to safely determine number of cores to use ----------------------
get_safe_cores <- function(max_cores = NULL, reserve_cores = 5) {
  available <- detectCores()
  
  if (is.null(max_cores)) {
    # Use 60-70% of available cores, leaving some for system
    safe_cores <- max(1, floor(available * 0.65))
  } else {
    safe_cores <- min(max_cores, available - reserve_cores)
  }
  
  safe_cores <- max(1, safe_cores)  # Ensure at least 1 core
  
  message(sprintf("Using %d cores out of %d available", safe_cores, available))
  return(safe_cores)
}

# Function to assess initial model guess ------------------------------------
assess_ctmm_guess <- function(telem_list) {
  
  message("\n=== Assessing Model Selection for Each Individual ===")
  
  for (i in seq_along(telem_list)) {
    tel_i <- telem_list[[i]]
    id_i <- names(telem_list)[i]
    
    # Run ctmm.guess
    guess_i <- ctmm.guess(
      tel_i,
      CTMM = ctmm(error = TRUE),
      interactive = FALSE
    )
    
    # Determine model type by tau length
    tau_len <- length(guess_i$tau)
    
    model_type <- case_when(
      tau_len == 0 ~ "IID (uncorrelated)",
      tau_len == 1 ~ "OU (range resident, no velocity autocorrelation)",
      tau_len == 2 ~ "OUF (directional persistence + range residency)",
      TRUE ~ "Unknown"
    )
    
    # Print summary
    message(sprintf("ID: %-10s | tau length = %d | guessed model = %s",
                    id_i, tau_len, model_type))
  }
}

# Parallel fitting function -------------------------------------------------
fit_ctmm_species_parallel <- function(telem_list, species_name, 
                                      max_cores = NULL, 
                                      save_individual_fits = TRUE) {
  
  message("\n", strrep("=", 80))
  message("=== FITTING CTMM MODELS FOR ", toupper(species_name), " ===")
  message(strrep("=", 80))
  
  # Assess expected models
  assess_ctmm_guess(telem_list)
  
  # Create output directory
  output_dir <- file.path(save_ctmm_path, paste0("lake_BT_", species_name, "_fits"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine number of cores
  n_cores <- get_safe_cores(max_cores = max_cores)
  
  # Set up parallel cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to cluster
  clusterExport(cl, c("telem_list", "species_name", "output_dir", 
                      "save_individual_fits", "save_ctmm_path"),
                envir = environment())
  
  message("\n=== Starting Parallel Model Fitting ===")
  message("Processing ", length(telem_list), " individuals using ", n_cores, " cores")
  start_time <- Sys.time()
  
  # Parallel fitting with error handling
  results <- foreach(
    i = 1:length(telem_list),
    .packages = c('ctmm'),
    .errorhandling = 'pass',  # Continue even if one fails
    .verbose = FALSE
  ) %dopar% {
    
    tryCatch({
      tel_i <- telem_list[[i]]
      id_i <- names(telem_list)[i]
      
      ind_start <- Sys.time()
      
      # Get initial guess
      guess_model <- ctmm.guess(
        tel_i,
        CTMM = ctmm(error = TRUE),
        interactive = FALSE
      )
      
      # Fit the model
      model_fit <- ctmm.fit(
        data = tel_i,
        CTMM = guess_model,
        method = "ML"
      )
      
      ind_elapsed <- as.numeric(difftime(Sys.time(), ind_start, units = "secs"))
      
      # Save individual fit if requested
      if (save_individual_fits) {
        output_file <- file.path(output_dir, paste0(id_i, "_ctmm_fit.rds"))
        saveRDS(model_fit, file = output_file)
      }
      
      # Clean up memory
      gc(verbose = FALSE)
      
      list(
        id = id_i,
        fit = model_fit,
        time = ind_elapsed,
        success = TRUE,
        error = NULL
      )
      
    }, error = function(e) {
      list(
        id = names(telem_list)[i],
        fit = NULL,
        time = NA,
        success = FALSE,
        error = as.character(e)
      )
    })
  }
  
  # Stop cluster
  stopCluster(cl)
  
  total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  # Process results
  successful_fits <- sapply(results, function(x) x$success)
  n_success <- sum(successful_fits)
  n_failed <- sum(!successful_fits)
  
  message("\n=== Fitting Complete ===")
  message(sprintf("Total time: %.1f minutes", total_elapsed))
  message(sprintf("Successful fits: %d/%d", n_success, length(telem_list)))
  
  if (n_failed > 0) {
    message(sprintf("\nFailed fits: %d", n_failed))
    failed_ids <- sapply(results[!successful_fits], function(x) x$id)
    message("Failed IDs: ", paste(failed_ids, collapse = ", "))
    
    # Print error messages
    for (i in which(!successful_fits)) {
      message(sprintf("\n%s error: %s", results[[i]]$id, results[[i]]$error))
    }
  }
  
  # Extract successful fits
  ctmm_fits <- lapply(results[successful_fits], function(x) x$fit)
  names(ctmm_fits) <- sapply(results[successful_fits], function(x) x$id)
  
  # Calculate timing statistics
  fit_times <- sapply(results[successful_fits], function(x) x$time)
  message(sprintf("\nTiming stats (seconds):"))
  message(sprintf("  Mean: %.1f  |  Median: %.1f  |  Max: %.1f", 
                  mean(fit_times), median(fit_times), max(fit_times)))
  
  # Save combined fit list
  if (length(ctmm_fits) > 0) {
    output_list_file <- file.path(output_dir, paste0("lake_BT_", species_name, "_ctmm_fits.rds"))
    saveRDS(ctmm_fits, output_list_file)
    message("\nCombined fit list saved to: ", output_list_file)
  }
  
  # Return both fits and diagnostic info
  return(list(
    fits = ctmm_fits,
    diagnostics = list(
      total_time = total_elapsed,
      n_success = n_success,
      n_failed = n_failed,
      failed_ids = if(n_failed > 0) sapply(results[!successful_fits], function(x) x$id) else NULL,
      fit_times = fit_times
    )
  ))
}

# Sequential fitting function (backup if parallel causes issues) -----------
fit_ctmm_species_sequential <- function(telem_list, species_name) {
  
  message("\n", strrep("=", 80))
  message("=== FITTING CTMM MODELS FOR ", toupper(species_name), " (SEQUENTIAL) ===")
  message(strrep("=", 80))
  
  assess_ctmm_guess(telem_list)
  
  ctmm_fits <- vector("list", length(telem_list))
  names(ctmm_fits) <- names(telem_list)
  
  output_dir <- file.path(save_ctmm_path, paste0("lake_BT_", species_name, "_fits"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("\n=== Starting Sequential Model Fitting ===")
  start_time <- Sys.time()
  fit_times <- numeric(length(telem_list))
  
  for (i in seq_along(telem_list)) {
    tel_i <- telem_list[[i]]
    id_i <- names(telem_list)[i]
    
    # Estimate time remaining
    if (i > 1) {
      avg_time <- mean(fit_times[1:(i-1)])
      remaining <- (length(telem_list) - i + 1) * avg_time / 60
      message(sprintf("\n[%d/%d] %s | Est. remaining: %.1f min", 
                      i, length(telem_list), id_i, remaining))
    } else {
      message(sprintf("\n[%d/%d] %s", i, length(telem_list), id_i))
    }
    
    ind_start <- Sys.time()
    
    guess_model <- ctmm.guess(tel_i, CTMM = ctmm(error = TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(data = tel_i, CTMM = guess_model, method = "ML")
    
    fit_times[i] <- as.numeric(difftime(Sys.time(), ind_start, units = "secs"))
    message(sprintf("  Fitted in %.1f seconds", fit_times[i]))
    
    output_file <- file.path(output_dir, paste0(id_i, "_ctmm_fit.rds"))
    saveRDS(model_fit, file = output_file)
    
    ctmm_fits[[i]] <- model_fit
    gc()
  }
  
  total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf("\n=== Completed in %.1f minutes ===", total_elapsed))
  
  output_list_file <- file.path(output_dir, paste0("lake_BT_", species_name, "_ctmm_fits.rds"))
  saveRDS(ctmm_fits, output_list_file)
  
  return(ctmm_fits)
}

# Function to verify model fits ---------------------------------------------
verify_fits <- function(telem_list, fit_list, species_name, n_check = 1) {
  
  message("\n=== Verifying ", species_name, " Model Fits ===")
  
  # Print summary of first model
  message("\nSummary of first individual (", names(fit_list)[1], "):")
  print(summary(fit_list[[1]]))
  
  # Optional: create diagnostic plot for first individual
  if (n_check > 0) {
    message("\nGenerating diagnostic plot for ", names(fit_list)[1])
    plot(telem_list[[1]], fit_list[[1]], error = FALSE)
  }
}

#==============================================================================
# 4. FIT CTMM MODELS - PIKE
#==============================================================================

# Option to reload if needed ------------------------------------------------
# pike_lake_BT_tel <- readRDS(paste0(save_telem_path, "BT/pike_lake_BT_tel_thinned.rds"))

# Fit models using parallel processing -------------------------------------
# Use max_cores=3 for pike (small group) to avoid overloading PC
pike_results <- fit_ctmm_species_parallel(
  pike_lake_BT_tel, 
  "pike",
  max_cores = 3  # Conservative setting to reduce noise/heat
)

lake_BT_pike_ctmm_fits <- pike_results$fits

# Verify fits ---------------------------------------------------------------
verify_fits(pike_lake_BT_tel, lake_BT_pike_ctmm_fits, "Pike")

#==============================================================================
# 5. FIT CTMM MODELS - PERCH
#==============================================================================

# Option to reload if needed ------------------------------------------------
# perch_lake_BT_tel <- readRDS(paste0(save_telem_path, "BT/perch_lake_BT_tel_thinned.rds"))

# Fit models using parallel processing -------------------------------------
# Use more cores for perch (larger group, 30 individuals)
perch_results <- fit_ctmm_species_parallel(
  perch_lake_BT_tel,
  "perch",
  max_cores = 6  # Can use more cores for larger groups
)

lake_BT_perch_ctmm_fits <- perch_results$fits

# Verify fits ---------------------------------------------------------------
verify_fits(perch_lake_BT_tel, lake_BT_perch_ctmm_fits, "Perch")

#==============================================================================
# 6. FIT CTMM MODELS - ROACH
#==============================================================================

# Option to reload if needed ------------------------------------------------
# roach_lake_BT_tel <- readRDS(paste0(save_telem_path, "BT/roach_lake_BT_tel_thinned.rds"))

# Fit models using parallel processing -------------------------------------
roach_results <- fit_ctmm_species_parallel(
  roach_lake_BT_tel,
  "roach",
  max_cores = 6  # Can use more cores for larger groups
)

lake_BT_roach_ctmm_fits <- roach_results$fits

# Verify fits ---------------------------------------------------------------
verify_fits(roach_lake_BT_tel, lake_BT_roach_ctmm_fits, "Roach")

#==============================================================================
# 7. FINAL SUMMARY
#==============================================================================

message("\n", strrep("=", 80))
message("=== ALL CTMM MODELS COMPLETED ===")
message(strrep("=", 80))

message("\nModels fitted:")
message("  Pike: ", length(lake_BT_pike_ctmm_fits), " individuals")
message("  Perch: ", length(lake_BT_perch_ctmm_fits), " individuals")
message("  Roach: ", length(lake_BT_roach_ctmm_fits), " individuals")
message("  Total: ", length(lake_BT_pike_ctmm_fits) + 
          length(lake_BT_perch_ctmm_fits) + 
          length(lake_BT_roach_ctmm_fits), " individuals")

message("\nTotal processing times:")
message("  Pike: ", sprintf("%.1f", pike_results$diagnostics$total_time), " minutes")
message("  Perch: ", sprintf("%.1f", perch_results$diagnostics$total_time), " minutes")
message("  Roach: ", sprintf("%.1f", roach_results$diagnostics$total_time), " minutes")

message("\nOutput locations:")
message("  Individual fits: ", save_ctmm_path, "lake_BT_[species]_fits/")
message("  Combined lists: ", save_ctmm_path, "lake_BT_[species]_fits/lake_BT_[species]_ctmm_fits.rds")
message("  Telemetry objects: ", save_telem_path, "BT/")

message("\nTo reload fitted models:")
message("  pike_fits <- readRDS('", save_ctmm_path, "lake_BT_pike_fits/lake_BT_pike_ctmm_fits.rds')")
message("  perch_fits <- readRDS('", save_ctmm_path, "lake_BT_perch_fits/lake_BT_perch_ctmm_fits.rds')")
message("  roach_fits <- readRDS('", save_ctmm_path, "lake_BT_roach_fits/lake_BT_roach_ctmm_fits.rds')")

message("\n", strrep("=", 80))

#==============================================================================
# NOTES ON USAGE
#==============================================================================

# If your PC gets too loud or hot during processing:
# 1. Reduce max_cores to 2-3 for all species
# 2. Or use the sequential version:
#    lake_BT_pike_ctmm_fits <- fit_ctmm_species_sequential(pike_lake_BT_tel, "pike")
#
# To adjust core usage dynamically:
# - max_cores = NULL will use ~65% of available cores
# - max_cores = 3 is conservative and quiet
# - max_cores = 6 is good for larger groups on modern PCs
#
# The parallel version will save significant time on the 30-individual spe