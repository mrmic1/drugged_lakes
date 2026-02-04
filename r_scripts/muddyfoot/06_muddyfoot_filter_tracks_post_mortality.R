#--------------------------------------------------------------------------#
# Filter data to remove post-predation and mortality tracking - Muddyfoot  #
#--------------------------------------------------------------------------#

# SCRIPT DESCRIPTION
# In this script we will:
# 1. Remove data where individuals were tracked after being predated or had a mortality event
# 2. Re-run ctmm model selection for individuals that had their data filtered
# 3. Use ctmm.select to find the best candidate movement model for each individual

#----------------------------------------------------------------------------------------------#

### LIBRARIES ###
library(dplyr)
library(ctmm)
library(data.table)
library(sf)
library(parallel)
library(foreach)
library(doParallel)

# Set the time zone environment to 'Europe/Stockholm' for consistent timestamp manipulation
Sys.setenv(TZ = 'Europe/Stockholm')

### DIRECTORIES ###
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
save_telem_path <- "./data/telem_obj/"
save_ctmm_path <- "./data/ctmm_fits/"
enc_path <- "./data/encounters/muddyfoot/"

#==============================================================================
# 1. FILTER OUT POST-PREDATION AND MORTALITY TRACKING LOCATIONS
#==============================================================================

# Load data
muddyfoot_telem_data <- readRDS(paste0(filtered_data_path, "04_muddyfoot_sub.rds"))
mortality_preds <- readxl::read_excel("./data/encounters/suspected_mortality_updated.xlsx")

message("\nOriginal dataset: ", nrow(muddyfoot_telem_data), " rows")

# Prepare mortality data
muddyfoot_pred_prey_cols <- mortality_preds %>%
  filter(lake == 'muddyfoot') %>% 
  filter(species == 'Roach' | species == 'Perch') %>%
  dplyr::select(individual_ID, species, revised_suspected_mortality, revised_likely_death_date) %>% 
  rename(death_date = revised_likely_death_date)

# Ensure death_date is Date type
muddyfoot_pred_prey_cols$death_date <- as.Date(muddyfoot_pred_prey_cols$death_date, origin = "1970-01-01")

message("\nIndividuals with mortality events:")
# Get species info from original data for summary
mortality_summary <- muddyfoot_telem_data %>%
  filter(individual_ID %in% muddyfoot_pred_prey_cols$individual_ID) %>%
  distinct(individual_ID, species) %>%
  left_join(muddyfoot_pred_prey_cols, by = "individual_ID")

print(table(mortality_summary$species, mortality_summary$revised_suspected_mortality))

# Filter out data after the predation or mortality event
muddyfoot_telem_data_2 <- muddyfoot_telem_data %>%
  left_join(muddyfoot_pred_prey_cols %>% select(individual_ID, death_date), 
            by = "individual_ID") %>%
  filter(is.na(death_date) | date < death_date) %>%
  select(-death_date)  # Remove the death_date column after filtering

rows_removed <- nrow(muddyfoot_telem_data) - nrow(muddyfoot_telem_data_2)
message("\nRows removed after filtering: ", format(rows_removed, big.mark = ","))
message("Filtered dataset: ", nrow(muddyfoot_telem_data_2), " rows")

# Save filtered data
saveRDS(muddyfoot_telem_data_2, paste0(filtered_data_path, "05_muddyfoot_sub.rds"))

#==============================================================================
# 2. HELPER FUNCTIONS
#==============================================================================

# Function to safely determine number of cores to use ----------------------
get_safe_cores <- function(max_cores = NULL, reserve_cores = 2) {
  available <- detectCores()
  
  if (is.null(max_cores)) {
    safe_cores <- max(1, floor(available * 0.65))
  } else {
    safe_cores <- min(max_cores, available - reserve_cores)
  }
  
  safe_cores <- max(1, safe_cores)
  message(sprintf("Using %d cores out of %d available", safe_cores, available))
  return(safe_cores)
}

# Function to assess initial model guess ------------------------------------
assess_ctmm_guess <- function(telem_list) {
  message("\n=== Assessing Model Selection for Each Individual ===")
  
  for (i in seq_along(telem_list)) {
    tel_i <- telem_list[[i]]
    id_i <- names(telem_list)[i]
    
    guess_i <- ctmm.guess(tel_i, CTMM = ctmm(error = TRUE), interactive = FALSE)
    tau_len <- length(guess_i$tau)
    
    model_type <- case_when(
      tau_len == 0 ~ "IID (uncorrelated)",
      tau_len == 1 ~ "OU (range resident, no velocity autocorrelation)",
      tau_len == 2 ~ "OUF (directional persistence + range residency)",
      TRUE ~ "Unknown"
    )
    
    message(sprintf("ID: %-10s | tau length = %d | guessed model = %s",
                    id_i, tau_len, model_type))
  }
}

# Parallel model selection function ----------------------------------------
refit_ctmm_models_parallel <- function(telem_list, species_name, 
                                       affected_ids,
                                       lake_name = "muddyfoot",
                                       max_cores = NULL, 
                                       ic = "AICc") {
  
  message("\n", strrep("=", 80))
  message("=== RE-FITTING CTMM MODELS FOR ", toupper(species_name), " ===")
  message("=== (Post-Mortality Filtering) ===")
  message(strrep("=", 80))
  
  # Filter to only affected individuals
  telem_subset <- telem_list[names(telem_list) %in% affected_ids]
  
  if (length(telem_subset) == 0) {
    message("\nNo individuals to refit for ", species_name)
    return(list(best_models = list(), diagnostics = list(n_success = 0)))
  }
  
  message("\nRefitting models for ", length(telem_subset), " individuals:")
  message(paste(names(telem_subset), collapse = ", "))
  
  # Assess expected models
  assess_ctmm_guess(telem_subset)
  
  # Output directory
  output_dir <- file.path(save_ctmm_path, paste0(lake_name, "_", species_name, "_fits"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine number of cores
  n_cores <- get_safe_cores(max_cores = max_cores)
  
  # Set up parallel cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export necessary objects
  clusterExport(cl, c("telem_subset", "species_name", "lake_name", "output_dir", 
                      "save_ctmm_path", "ic"),
                envir = environment())
  
  message("\n=== Starting Parallel Model Selection ===")
  message("Processing ", length(telem_subset), " individuals using ", n_cores, " cores")
  message("Information criterion: ", ic)
  start_time <- Sys.time()
  
  # Parallel processing
  results <- foreach(
    i = 1:length(telem_subset),
    .packages = c('ctmm'),
    .errorhandling = 'pass',
    .verbose = FALSE
  ) %dopar% {
    
    tryCatch({
      tel_i <- telem_subset[[i]]
      id_i <- names(telem_subset)[i]
      
      ind_start <- Sys.time()
      
      # Get initial guess
      guess_model <- ctmm.guess(tel_i, CTMM = ctmm(error = TRUE), interactive = FALSE)
      
      # Use ctmm.select to find best model
      model_selection <- ctmm.select(
        data = tel_i,
        CTMM = guess_model,
        method = "ML",
        IC = ic,
        verbose = TRUE
      )
      
      # Extract best model
      best_model <- model_selection[[1]]
      
      ind_elapsed <- as.numeric(difftime(Sys.time(), ind_start, units = "secs"))
      
      # Save results
      output_file_selection <- file.path(output_dir, paste0(id_i, "_ctmm_selection_refitted.rds"))
      output_file_best <- file.path(output_dir, paste0(id_i, "_ctmm_best_fit_refitted.rds"))
      
      saveRDS(model_selection, file = output_file_selection)
      saveRDS(best_model, file = output_file_best)
      
      # Also overwrite the original files
      saveRDS(model_selection, file.path(output_dir, paste0(id_i, "_ctmm_selection.rds")))
      saveRDS(best_model, file.path(output_dir, paste0(id_i, "_ctmm_best_fit.rds")))
      
      gc(verbose = FALSE)
      
      list(
        id = id_i,
        selection = model_selection,
        best_fit = best_model,
        time = ind_elapsed,
        success = TRUE,
        error = NULL
      )
      
    }, error = function(e) {
      list(
        id = names(telem_subset)[i],
        selection = NULL,
        best_fit = NULL,
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
  successful <- sapply(results, function(x) x$success)
  n_success <- sum(successful)
  n_failed <- sum(!successful)
  
  message("\n=== Model Selection Complete ===")
  message(sprintf("Total time: %.1f minutes", total_elapsed))
  message(sprintf("Successful selections: %d/%d", n_success, length(telem_subset)))
  
  if (n_failed > 0) {
    message(sprintf("\nFailed selections: %d", n_failed))
    failed_ids <- sapply(results[!successful], function(x) x$id)
    message("Failed IDs: ", paste(failed_ids, collapse = ", "))
    
    for (i in which(!successful)) {
      message(sprintf("\n%s error: %s", results[[i]]$id, results[[i]]$error))
    }
  }
  
  # Extract results
  best_models_list <- lapply(results[successful], function(x) x$best_fit)
  names(best_models_list) <- sapply(results[successful], function(x) x$id)
  
  selection_list <- lapply(results[successful], function(x) x$selection)
  names(selection_list) <- sapply(results[successful], function(x) x$id)
  
  # Timing statistics
  fit_times <- sapply(results[successful], function(x) x$time)
  message(sprintf("\nTiming stats (seconds):"))
  message(sprintf("  Mean: %.1f  |  Median: %.1f  |  Max: %.1f", 
                  mean(fit_times), median(fit_times), max(fit_times)))
  
  # Model summary
  message("\n=== Model Selection Summary ===")
  selected_models <- sapply(best_models_list, function(x) summary(x)$name)
  model_table <- table(selected_models)
  message("Selected models across individuals:")
  print(model_table)
  
  # Save refitted models list
  if (length(best_models_list) > 0) {
    output_refitted_file <- file.path(output_dir, 
                                      paste0(lake_name, "_", species_name, "_refitted_best_models.rds"))
    saveRDS(best_models_list, output_refitted_file)
    message("\nRefitted models saved to: ", output_refitted_file)
  }
  
  message("\n", strrep("=", 80))
  
  return(list(
    best_models = best_models_list,
    selection_results = selection_list,
    all_results = results,
    diagnostics = list(
      total_time = total_elapsed,
      n_success = n_success,
      n_failed = n_failed,
      failed_ids = if(n_failed > 0) sapply(results[!successful], function(x) x$id) else NULL,
      fit_times = fit_times,
      model_counts = as.list(model_table)
    )
  ))
}

#==============================================================================
# 3. PREPARE DATA FOR REFITTING
#==============================================================================


# Load filtered data
muddyfoot_filt_data <- readRDS(paste0(filtered_data_path, "05_muddyfoot_sub.rds"))
muddyfoot_UERE <- readRDS(paste0(save_telem_path, "muddyfoot/muddyfoot_UERE.rds"))

# Prepare mortality data for identifying affected individuals
muddyfoot_pred_prey_cols <- mortality_preds %>%
  filter(lake == 'muddyfoot') %>% 
  filter(species == 'Roach' | species == 'Perch') %>%
  dplyr::select(individual_ID, species, revised_suspected_mortality, revised_likely_death_date) %>% 
  rename(death_date = revised_likely_death_date)

muddyfoot_pred_prey_cols$death_date <- as.Date(muddyfoot_pred_prey_cols$death_date, origin = "1970-01-01")

# Identify affected individuals by species
perch_pred_ids <- muddyfoot_pred_prey_cols %>%
  filter(species == 'Perch') %>% 
  filter(revised_suspected_mortality %in% c('confirmed_mortality', 
                                            'confirmed_predated',
                                            'likely_predated', 
                                            'known_predated')) %>% 
  pull(individual_ID)

roach_pred_ids <- muddyfoot_pred_prey_cols %>%
  filter(species == 'Roach') %>% 
  filter(revised_suspected_mortality %in% c('confirmed_mortality', 
                                            'confirmed_predated',
                                            'likely_predated', 
                                            'known_predated')) %>% 
  filter(individual_ID != 'F59707') %>%  # Exclude this individual
  pull(individual_ID)

message("\nAffected individuals:")
message("  Perch: ", length(perch_pred_ids), " individuals")
if (length(perch_pred_ids) > 0) message("    IDs: ", paste(perch_pred_ids, collapse = ", "))
message("  Roach: ", length(roach_pred_ids), " individuals")
if (length(roach_pred_ids) > 0) message("    IDs: ", paste(roach_pred_ids, collapse = ", "))

# Verify data differences and identify individuals with no data remaining
message("\n=== Verifying Data Changes ===")

# Load the BEFORE filtering data for comparison
muddyfoot_telem_data_before <- readRDS(paste0(filtered_data_path, "04_muddyfoot_sub.rds"))

# Initialize vectors to track individuals with no data
perch_no_data <- character(0)
roach_no_data <- character(0)

if (length(perch_pred_ids) > 0) {
  message("\nPerch - Row counts before/after filtering:")
  
  filt_counts <- muddyfoot_filt_data %>%
    filter(individual_ID %in% perch_pred_ids) %>%
    count(individual_ID, name = "after_filtering")
  
  telem_counts <- muddyfoot_telem_data_before %>%
    filter(individual_ID %in% perch_pred_ids) %>%
    count(individual_ID, name = "before_filtering")
  
  row_comparison <- full_join(telem_counts, filt_counts, by = "individual_ID") %>%
    mutate(after_filtering = ifelse(is.na(after_filtering), 0, after_filtering),
           rows_removed = before_filtering - after_filtering)
  
  print(row_comparison)
  
  # Identify individuals with no data remaining
  perch_no_data <- row_comparison %>%
    filter(after_filtering == 0) %>%
    pull(individual_ID)
  
  if (length(perch_no_data) > 0) {
    message("\n*** WARNING: Perch individuals with NO data remaining after filtering:")
    message("    ", paste(perch_no_data, collapse = ", "))
    message("    These will be excluded from model refitting.")
  }
}

if (length(roach_pred_ids) > 0) {
  message("\nRoach - Row counts before/after filtering:")
  
  filt_counts <- muddyfoot_filt_data %>%
    filter(individual_ID %in% roach_pred_ids) %>%
    count(individual_ID, name = "after_filtering")
  
  telem_counts <- muddyfoot_telem_data_before %>%
    filter(individual_ID %in% roach_pred_ids) %>%
    count(individual_ID, name = "before_filtering")
  
  row_comparison <- full_join(telem_counts, filt_counts, by = "individual_ID") %>%
    mutate(after_filtering = ifelse(is.na(after_filtering), 0, after_filtering),
           rows_removed = before_filtering - after_filtering)
  
  print(row_comparison)
  
  # Identify individuals with no data remaining
  roach_no_data <- row_comparison %>%
    filter(after_filtering == 0) %>%
    pull(individual_ID)
  
  if (length(roach_no_data) > 0) {
    message("\n*** WARNING: Roach individuals with NO data remaining after filtering:")
    message("    ", paste(roach_no_data, collapse = ", "))
    message("    These will be excluded from model refitting.")
  }
}

# Remove individuals with no data from the ID lists
if (length(perch_no_data) > 0) {
  perch_pred_ids <- setdiff(perch_pred_ids, perch_no_data)
  message("\nPerch individuals to refit after exclusions: ", length(perch_pred_ids))
  if (length(perch_pred_ids) > 0) {
    message("  IDs: ", paste(perch_pred_ids, collapse = ", "))
  }
}

if (length(roach_no_data) > 0) {
  roach_pred_ids <- setdiff(roach_pred_ids, roach_no_data)
  message("\nRoach individuals to refit after exclusions: ", length(roach_pred_ids))
  if (length(roach_pred_ids) > 0) {
    message("  IDs: ", paste(roach_pred_ids, collapse = ", "))
  }
}

#==============================================================================
# 5. REFIT PERCH MODELS
#==============================================================================

if (length(perch_pred_ids) > 0) {
  
  message("\n", strrep("=", 80))
  message("=== REFITTING PERCH MODELS ===")
  message(strrep("=", 80))
  
  # Prepare perch data (only for individuals with data remaining)
  pred_perch_data <- muddyfoot_filt_data %>% 
    filter(individual_ID %in% perch_pred_ids)
  
  # Verify we have data
  if (nrow(pred_perch_data) == 0) {
    message("\n*** No perch data available for refitting after filtering ***")
    muddyfoot_perch_refitted <- list()
  } else {
    
    pred_perch_data <- with(pred_perch_data, 
                            data.frame(
                              "timestamp" = timestamp,                        
                              "location.long" = Long,                         
                              "location.lat" = Lat, 
                              "GPS.HDOP" = HDOP,                               
                              "individual-local-identifier" = individual_ID))
    
    # Create telemetry object
    pred_perch_tel <- as.telemetry(pred_perch_data, 
                                   timezone = "Europe/Stockholm",   
                                   timeformat = "%Y-%m-%d %H:%M:%S",
                                   projection = NULL,               
                                   datum = "WGS84")
    
    ctmm::projection(pred_perch_tel) <- ctmm::median(pred_perch_tel)
    uere(pred_perch_tel) <- muddyfoot_UERE
    
    # Refit models
    perch_refit_results <- refit_ctmm_models_parallel(
      telem_list = pred_perch_tel,
      species_name = "perch",
      affected_ids = perch_pred_ids,
      lake_name = "muddyfoot",
      max_cores = 4,
      ic = "AICc"
    )
    
    muddyfoot_perch_refitted <- perch_refit_results$best_models
    
    message("\n=== Perch Refitting Summary ===")
    message("Successfully refitted: ", perch_refit_results$diagnostics$n_success, " individuals")
    if (perch_refit_results$diagnostics$n_success > 0) {
      message("Model distribution:")
      print(perch_refit_results$diagnostics$model_counts)
    }
  }
  
} else {
  message("\n*** No perch individuals require refitting (all filtered out or no data) ***")
  muddyfoot_perch_refitted <- list()
}

#==============================================================================
# 6. REFIT ROACH MODELS
#==============================================================================

roach_ids <- c("F59710", "F59719")

if (length(roach_ids) > 0) {
  
  message("\n", strrep("=", 80))
  message("=== REFITTING ROACH MODELS ===")
  message(strrep("=", 80))
  
  # Prepare roach data (only for individuals with data remaining)
  pred_roach_data <- muddyfoot_filt_data %>% 
    filter(individual_ID %in% roach_ids)
  
  # Verify we have data
  if (nrow(pred_roach_data) == 0) {
    message("\n*** No roach data available for refitting after filtering ***")
    muddyfoot_roach_refitted <- list()
  } else {
    
    pred_roach_data <- with(pred_roach_data, 
                            data.frame(
                              "timestamp" = timestamp,                        
                              "location.long" = Long,                         
                              "location.lat" = Lat, 
                              "GPS.HDOP" = HDOP,                               
                              "individual-local-identifier" = individual_ID))
    
    # Create telemetry object
    pred_roach_tel <- as.telemetry(pred_roach_data, 
                                   timezone = "Europe/Stockholm",   
                                   timeformat = "%Y-%m-%d %H:%M:%S",
                                   projection = NULL,               
                                   datum = "WGS84")
    
    ctmm::projection(pred_roach_tel) <- ctmm::median(pred_roach_tel)
    uere(pred_roach_tel) <- muddyfoot_UERE
    
    # Refit models
    roach_refit_results <- refit_ctmm_models_parallel(
      telem_list = pred_roach_tel,
      species_name = "roach",
      affected_ids = roach_ids,
      lake_name = "muddyfoot",
      max_cores = 2,
      ic = "AICc"
    )
    
    muddyfoot_roach_refitted <- roach_refit_results$best_models
    
    message("\n=== Roach Refitting Summary ===")
    message("Successfully refitted: ", roach_refit_results$diagnostics$n_success, " individuals")
    if (roach_refit_results$diagnostics$n_success > 0) {
      message("Model distribution:")
      print(roach_refit_results$diagnostics$model_counts)
    }
  }
  
} else {
  message("\n*** No roach individuals require refitting (all filtered out or no data) ***")
  muddyfoot_roach_refitted <- list()
}

#==============================================================================
# 7. UPDATE MASTER MODEL LISTS AND TELEMETRY OBJECTS
#==============================================================================

message("\n", strrep("=", 80))
message("=== UPDATING MASTER MODEL LISTS AND TELEMETRY OBJECTS ===")
message(strrep("=", 80))

# Function to update model lists and remove individuals with no data
update_model_list <- function(species_name, refitted_models, no_data_ids) {
  
  # Load existing models
  existing_best_file <- file.path(save_ctmm_path, 
                                  paste0("muddyfoot_", species_name, "_fits"),
                                  paste0("muddyfoot_", species_name, "_best_models.rds"))
  
  if (file.exists(existing_best_file)) {
    existing_models <- readRDS(existing_best_file)
    
    # Remove individuals with no data remaining
    if (length(no_data_ids) > 0) {
      message("\nRemoving ", species_name, " individuals with no data after filtering:")
      message("  IDs: ", paste(no_data_ids, collapse = ", "))
      existing_models <- existing_models[!names(existing_models) %in% no_data_ids]
    }
    
    # Update with refitted models
    if (length(refitted_models) > 0) {
      for (id in names(refitted_models)) {
        existing_models[[id]] <- refitted_models[[id]]
      }
      message("\nUpdated ", species_name, " master model list with ", 
              length(refitted_models), " refitted models")
    }
    
    # Save updated list
    saveRDS(existing_models, existing_best_file)
    message("Saved to: ", existing_best_file)
    message("Final ", species_name, " model count: ", length(existing_models))
    
  } else {
    warning("\nCould not find existing model file: ", existing_best_file)
    message("Refitted models saved separately but not merged into master list")
  }
}

# Function to update telemetry object lists
update_telemetry_list <- function(species_name, refitted_tel, no_data_ids) {
  
  # Load existing telemetry objects
  existing_tel_file <- file.path(save_telem_path, "muddyfoot",
                                 paste0(species_name, "_muddyfoot_tel_thinned.rds"))
  
  if (file.exists(existing_tel_file)) {
    existing_tel <- readRDS(existing_tel_file)
    
    message("\n--- Updating ", species_name, " telemetry objects ---")
    message("Original telemetry list: ", length(existing_tel), " individuals")
    
    # Remove individuals with no data remaining
    if (length(no_data_ids) > 0) {
      message("Removing individuals with no data: ", paste(no_data_ids, collapse = ", "))
      existing_tel <- existing_tel[!names(existing_tel) %in% no_data_ids]
    }
    
    # Update with refitted telemetry objects
    if (length(refitted_tel) > 0) {
      message("Updating with refitted telemetry: ", length(refitted_tel), " individuals")
      for (id in names(refitted_tel)) {
        existing_tel[[id]] <- refitted_tel[[id]]
      }
    }
    
    message("Final telemetry list: ", length(existing_tel), " individuals")
    
    # Save updated telemetry list
    output_tel_file <- file.path(save_telem_path, "muddyfoot",
                                 paste0(species_name, "_muddyfoot_tel_thinned_and_filtered.rds"))
    saveRDS(existing_tel, output_tel_file)
    message("Saved to: ", output_tel_file)
    
    return(existing_tel)
    
  } else {
    warning("\nCould not find existing telemetry file: ", existing_tel_file)
    message("Refitted telemetry saved separately but not merged into master list")
    return(NULL)
  }
}

# Update master model lists
if (exists("muddyfoot_perch_refitted")) {
  update_model_list("perch", muddyfoot_perch_refitted, perch_no_data)
}

if (exists("muddyfoot_roach_refitted")) {
  update_model_list("roach", muddyfoot_roach_refitted, roach_no_data)
}

# Update and save telemetry object lists
message("\n", strrep("-", 80))
message("--- UPDATING TELEMETRY OBJECT LISTS ---")
message(strrep("-", 80))

# Update perch telemetry
if (exists("pred_perch_tel") && length(perch_pred_ids) > 0) {
  perch_tel_updated <- update_telemetry_list("perch", pred_perch_tel, perch_no_data)
} else {
  # Still create filtered version even if no refitting occurred
  if (length(perch_no_data) > 0) {
    existing_tel_file <- file.path(save_telem_path, "muddyfoot", "perch_muddyfoot_tel_thinned.rds")
    if (file.exists(existing_tel_file)) {
      existing_tel <- readRDS(existing_tel_file)
      message("\n--- Updating perch telemetry objects ---")
      message("Original telemetry list: ", length(existing_tel), " individuals")
      message("Removing individuals with no data: ", paste(perch_no_data, collapse = ", "))
      existing_tel <- existing_tel[!names(existing_tel) %in% perch_no_data]
      message("Final telemetry list: ", length(existing_tel), " individuals")
      
      output_tel_file <- file.path(save_telem_path, "muddyfoot", 
                                   "perch_muddyfoot_tel_thinned_and_filtered.rds")
      saveRDS(existing_tel, output_tel_file)
      message("Saved to: ", output_tel_file)
      perch_tel_updated <- existing_tel
    }
  } else {
    # No changes needed, just copy the original
    existing_tel_file <- file.path(save_telem_path, "muddyfoot", "perch_muddyfoot_tel_thinned.rds")
    if (file.exists(existing_tel_file)) {
      existing_tel <- readRDS(existing_tel_file)
      output_tel_file <- file.path(save_telem_path, "muddyfoot", 
                                   "perch_muddyfoot_tel_thinned_and_filtered.rds")
      saveRDS(existing_tel, output_tel_file)
      message("\nNo perch individuals affected - copied original telemetry to filtered version")
      message("Saved to: ", output_tel_file)
      perch_tel_updated <- existing_tel
    }
  }
}

# Update roach telemetry
if (exists("pred_roach_tel") && length(roach_pred_ids) > 0) {
  roach_tel_updated <- update_telemetry_list("roach", pred_roach_tel, roach_no_data)
} else {
  # Still create filtered version even if no refitting occurred
  if (length(roach_no_data) > 0) {
    existing_tel_file <- file.path(save_telem_path, "muddyfoot", "roach_muddyfoot_tel_thinned.rds")
    if (file.exists(existing_tel_file)) {
      existing_tel <- readRDS(existing_tel_file)
      message("\n--- Updating roach telemetry objects ---")
      message("Original telemetry list: ", length(existing_tel), " individuals")
      message("Removing individuals with no data: ", paste(roach_no_data, collapse = ", "))
      existing_tel <- existing_tel[!names(existing_tel) %in% roach_no_data]
      message("Final telemetry list: ", length(existing_tel), " individuals")
      
      output_tel_file <- file.path(save_telem_path, "muddyfoot", 
                                   "roach_muddyfoot_tel_thinned_and_filtered.rds")
      saveRDS(existing_tel, output_tel_file)
      message("Saved to: ", output_tel_file)
      roach_tel_updated <- existing_tel
    }
  } else {
    # No changes needed, just copy the original
    existing_tel_file <- file.path(save_telem_path, "muddyfoot", "roach_muddyfoot_tel_thinned.rds")
    if (file.exists(existing_tel_file)) {
      existing_tel <- readRDS(existing_tel_file)
      output_tel_file <- file.path(save_telem_path, "muddyfoot", 
                                   "roach_muddyfoot_tel_thinned_and_filtered.rds")
      saveRDS(existing_tel, output_tel_file)
      message("\nNo roach individuals affected - copied original telemetry to filtered version")
      message("Saved to: ", output_tel_file)
      roach_tel_updated <- existing_tel
    }
  }
}

#==============================================================================
# 8. FINAL SUMMARY
#==============================================================================

message("\n", strrep("=", 80))
message("=== POST-MORTALITY FILTERING AND REFITTING COMPLETE ===")
message(strrep("=", 80))

message("\nData Filtering Summary:")
message("  Total rows removed: ", format(rows_removed, big.mark = ","))
message("  Filtered data saved to: ", paste0(filtered_data_path, "05_muddyfoot_sub.rds"))

# Summary of individuals with no data
if (length(perch_no_data) > 0 || length(roach_no_data) > 0) {
  message("\nIndividuals with NO data remaining after filtering:")
  if (length(perch_no_data) > 0) {
    message("  Perch (", length(perch_no_data), "): ", paste(perch_no_data, collapse = ", "))
  }
  if (length(roach_no_data) > 0) {
    message("  Roach (", length(roach_no_data), "): ", paste(roach_no_data, collapse = ", "))
  }
  message("  These individuals have been REMOVED from model and telemetry lists")
}

message("\nModel Refitting Summary:")

if (exists("perch_refit_results")) {
  message("\n  Perch:")
  message("    Individuals refitted: ", perch_refit_results$diagnostics$n_success)
  message("    Processing time: ", sprintf("%.1f", perch_refit_results$diagnostics$total_time), " minutes")
  if (perch_refit_results$diagnostics$n_success > 0) {
    message("    Models selected:")
    for (model_name in names(perch_refit_results$diagnostics$model_counts)) {
      message("      ", model_name, ": ", perch_refit_results$diagnostics$model_counts[[model_name]])
    }
  }
}

if (exists("roach_refit_results")) {
  message("\n  Roach:")
  message("    Individuals refitted: ", roach_refit_results$diagnostics$n_success)
  message("    Processing time: ", sprintf("%.1f", roach_refit_results$diagnostics$total_time), " minutes")
  if (roach_refit_results$diagnostics$n_success > 0) {
    message("    Models selected:")
    for (model_name in names(roach_refit_results$diagnostics$model_counts)) {
      message("      ", model_name, ": ", roach_refit_results$diagnostics$model_counts[[model_name]])
    }
  }
}

message("\nTelemetry Object Summary:")
if (exists("perch_tel_updated")) {
  message("  Perch telemetry: ", length(perch_tel_updated), " individuals")
}
if (exists("roach_tel_updated")) {
  message("  Roach telemetry: ", length(roach_tel_updated), " individuals")
}

message("\nOutput Files:")
message("  Updated models:")
message("    ", save_ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_best_models.rds")
message("    ", save_ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_best_models.rds")
message("  Updated telemetry objects:")
message("    ", save_telem_path, "muddyfoot/perch_muddyfoot_tel_thinned_and_filtered.rds")
message("    ", save_telem_path, "muddyfoot/roach_muddyfoot_tel_thinned_and_filtered.rds")

message("\nTo reload updated data:")
message("  # Models")
message("  perch_models <- readRDS('", save_ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_best_models.rds')")
message("  roach_models <- readRDS('", save_ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_best_models.rds')")
message("\n  # Telemetry objects")
message("  perch_tel <- readRDS('", save_telem_path, "muddyfoot/perch_muddyfoot_tel_thinned_and_filtered.rds')")
message("  roach_tel <- readRDS('", save_telem_path, "muddyfoot/roach_muddyfoot_tel_thinned_and_filtered.rds')")

message("\n", strrep("=", 80))