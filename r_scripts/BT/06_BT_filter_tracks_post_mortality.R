#--------------------------------------------------------------------------#
# Filter data to remove post-predation and mortality tracking - BT  #
#--------------------------------------------------------------------------#

# SCRIPT DESCRIPTION
# In this script we will:
# 1. Remove data where individuals were tracked after being predated or had a mortality event
# 2. Re-run ctmm model selection for all individuals witht the filtered dataset
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
filtered_data_path <- "./data/tracks_filtered/BT/"
save_telem_path <- "./data/telem_obj/"
save_ctmm_path <- "./data/ctmm_fits/"
enc_path <- "./data/encounters/BT/"

#==============================================================================-
# 1. FILTER OUT POST-PREDATION AND MORTALITY TRACKING LOCATIONS ####
#==============================================================================-

# Load data
BT_telem_data <- readRDS(paste0(filtered_data_path, "04_BT_sub.rds"))
mortality_preds <- readxl::read_excel("./data/encounters/suspected_mortality_updated.xlsx")
message("\nOriginal dataset: ", nrow(BT_telem_data), " rows")
pike_mortality <- read.csv("./data/pike_deaths.csv")

# Prepare mortality data
BT_pred_prey_cols <- mortality_preds %>%
  filter(lake == 'BT') %>% 
  filter(species == 'Roach' | species == 'Perch') %>%
  dplyr::select(individual_ID, revised_suspected_mortality, revised_likely_death_date) %>%  # Remove species here
  rename(death_date = revised_likely_death_date)

# Ensure death_date is Date type
BT_pred_prey_cols$death_date <- as.Date(BT_pred_prey_cols$death_date, origin = "1970-01-01")

message("\nIndividuals with mortality events:")
# Get species info from original data for summary
mortality_summary <- 
  BT_telem_data %>%
  filter(individual_ID %in% BT_pred_prey_cols$individual_ID) %>%
  distinct(individual_ID, species) %>%
  left_join(BT_pred_prey_cols, by = "individual_ID")

print(table(mortality_summary$species, mortality_summary$revised_suspected_mortality))

# Filter out data after the predation or mortality event
BT_telem_data_2 <- BT_telem_data %>%
  left_join(BT_pred_prey_cols %>% select(individual_ID, death_date), 
            by = "individual_ID") %>%
  filter(is.na(death_date) | date < death_date) %>%
  select(-death_date)  # Remove the death_date column after filtering

rows_removed <- nrow(BT_telem_data) - nrow(BT_telem_data_2)
message("\nRows removed after filtering: ", format(rows_removed, big.mark = ",")) #Rows removed after filtering: 168,308
message("Filtered dataset: ", nrow(BT_telem_data_2), " rows") #2859218 rows

#Now I need to filter out pike mortality
BT_pike_mort <- pike_mortality %>%
  filter(lake == 'BT') %>% 
  dplyr::select(individual_ID, likely_death_date) %>%  # Remove species here
  rename(death_date = likely_death_date)

# Ensure death_date is Date type
BT_pike_mort$death_date <- as.Date(BT_pike_mort$death_date, origin = "1970-01-01")

# Filter out data after the predation or mortality event
BT_telem_data_3 <- BT_telem_data_2 %>%
  left_join(BT_pike_mort %>% select(individual_ID, death_date), 
            by = "individual_ID") %>%
  filter(is.na(death_date) | date < death_date) %>%
  select(-death_date)  # Remove the death_date column after filtering

rows_removed <- nrow(BT_telem_data_2) - nrow(BT_telem_data_3)
message("\nRows removed after filtering: ", format(rows_removed, big.mark = ",")) #Rows removed after filtering: 31,886
message("Filtered dataset: ", nrow(BT_telem_data_3), " rows") #2827332 rows


# Save filtered data
saveRDS(BT_telem_data_2, paste0(filtered_data_path, "05_BT_sub.rds"))

#==============================================================================-
# 2. HELPER FUNCTIONS ####
#==============================================================================-

# Function to safely determine number of cores to use ----------------------
get_safe_cores <- function(max_cores = NULL, reserve_cores = 2) {
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

# Parallel model selection function (UPDATED) ------------------------------
fit_ctmm_species_parallel <- function(telem_list, species_name, lake_name = "BT",
                                      max_cores = NULL, 
                                      save_individual_fits = TRUE,
                                      ic = "AICc") {
  
  message("\n", strrep("=", 80))
  message("=== SELECTING BEST CTMM MODELS FOR ", toupper(species_name), " ===")
  message(strrep("=", 80))
  
  # Assess expected models
  assess_ctmm_guess(telem_list)
  
  # Create output directory
  output_dir <- file.path(save_ctmm_path, paste0(lake_name, "_", species_name, "_fits"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine number of cores
  n_cores <- get_safe_cores(max_cores = max_cores)
  
  # Set up parallel cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to cluster
  clusterExport(cl, c("telem_list", "species_name", "lake_name", "output_dir", 
                      "save_individual_fits", "save_ctmm_path", "ic"),
                envir = environment())
  
  message("\n=== Starting Parallel Model Selection ===")
  message("Processing ", length(telem_list), " individuals using ", n_cores, " cores")
  message("Information criterion: ", ic)
  start_time <- Sys.time()
  
  # Parallel fitting with error handling
  results <- foreach(
    i = 1:length(telem_list),
    .packages = c('ctmm'),
    .errorhandling = 'pass',
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
      
      # Use ctmm.select to find best model
      model_selection <- ctmm.select(
        data = tel_i,
        CTMM = guess_model,
        method = "ML",
        IC = ic,
        verbose = TRUE
      )
      
      # Extract the best model (first in the list)
      best_model <- model_selection[[1]]
      
      ind_elapsed <- as.numeric(difftime(Sys.time(), ind_start, units = "secs"))
      
      # Save individual fit if requested
      if (save_individual_fits) {
        # Save both the selection results and best model
        output_file_selection <- file.path(output_dir, paste0(id_i, "_ctmm_selection.rds"))
        output_file_best <- file.path(output_dir, paste0(id_i, "_ctmm_best_fit.rds"))
        
        saveRDS(model_selection, file = output_file_selection)
        saveRDS(best_model, file = output_file_best)
      }
      
      # Clean up memory
      gc(verbose = FALSE)
      
      list(
        id = id_i,
        selection = model_selection,  # Full selection results
        best_fit = best_model,        # Best model only
        time = ind_elapsed,
        success = TRUE,
        error = NULL
      )
      
    }, error = function(e) {
      list(
        id = names(telem_list)[i],
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
  successful_fits <- sapply(results, function(x) x$success)
  n_success <- sum(successful_fits)
  n_failed <- sum(!successful_fits)
  
  message("\n=== Model Selection Complete ===")
  message(sprintf("Total time: %.1f minutes", total_elapsed))
  message(sprintf("Successful selections: %d/%d", n_success, length(telem_list)))
  
  if (n_failed > 0) {
    message(sprintf("\nFailed selections: %d", n_failed))
    failed_ids <- sapply(results[!successful_fits], function(x) x$id)
    message("Failed IDs: ", paste(failed_ids, collapse = ", "))
    
    # Print error messages
    for (i in which(!successful_fits)) {
      message(sprintf("\n%s error: %s", results[[i]]$id, results[[i]]$error))
    }
  }
  
  # Extract best models list
  best_models_list <- lapply(results[successful_fits], function(x) x$best_fit)
  names(best_models_list) <- sapply(results[successful_fits], function(x) x$id)
  
  # Extract full selection results list
  selection_list <- lapply(results[successful_fits], function(x) x$selection)
  names(selection_list) <- sapply(results[successful_fits], function(x) x$id)
  
  # Calculate timing statistics
  fit_times <- sapply(results[successful_fits], function(x) x$time)
  message(sprintf("\nTiming stats (seconds):"))
  message(sprintf("  Mean: %.1f  |  Median: %.1f  |  Max: %.1f", 
                  mean(fit_times), median(fit_times), max(fit_times)))
  
  # Print model summary
  message("\n=== Model Selection Summary ===")
  selected_models <- sapply(best_models_list, function(x) summary(x)$name)
  model_table <- table(selected_models)
  message("Selected models across individuals:")
  print(model_table)
  
  # Save combined lists
  if (length(best_models_list) > 0) {
    output_best_file <- file.path(output_dir, paste0(lake_name, "_", species_name, "_best_models.rds"))
    output_selection_file <- file.path(output_dir, paste0(lake_name, "_", species_name, "_all_selections.rds"))
    
    saveRDS(best_models_list, output_best_file)
    saveRDS(selection_list, output_selection_file)
    
    message("\nBest models list saved to: ", output_best_file)
    message("Full selection results saved to: ", output_selection_file)
  }
  
  message("\n", strrep("=", 80))
  
  # Return comprehensive results
  return(list(
    best_models = best_models_list,      # List of best ctmm objects
    selection_results = selection_list,  # Full selection results for each individual
    all_results = results,
    diagnostics = list(
      total_time = total_elapsed,
      n_success = n_success,
      n_failed = n_failed,
      failed_ids = if(n_failed > 0) sapply(results[!successful_fits], function(x) x$id) else NULL,
      fit_times = fit_times,
      model_counts = as.list(model_table)
    )
  ))
}

# Sequential model selection function (backup) -----------------------------
fit_ctmm_species_sequential <- function(telem_list, species_name, lake_name = "BT",
                                        ic = "AICc") {
  
  message("\n", strrep("=", 80))
  message("=== SELECTING BEST CTMM MODELS FOR ", toupper(species_name), " (SEQUENTIAL) ===")
  message(strrep("=", 80))
  
  assess_ctmm_guess(telem_list)
  
  best_models_list <- vector("list", length(telem_list))
  selection_list <- vector("list", length(telem_list))
  names(best_models_list) <- names(telem_list)
  names(selection_list) <- names(telem_list)
  
  output_dir <- file.path(save_ctmm_path, paste0(lake_name, "_", species_name, "_fits"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("\n=== Starting Sequential Model Selection ===")
  message("Information criterion: ", ic)
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
    
    # Use ctmm.select
    model_selection <- ctmm.select(
      data = tel_i,
      CTMM = guess_model,
      method = "ML",
      IC = ic,
      verbose = TRUE
    )
    
    best_model <- model_selection[[1]]
    
    fit_times[i] <- as.numeric(difftime(Sys.time(), ind_start, units = "secs"))
    message(sprintf("  Selected best model in %.1f seconds: %s", 
                    fit_times[i], summary(best_model)$name))
    
    # Save files
    output_file_selection <- file.path(output_dir, paste0(id_i, "_ctmm_selection.rds"))
    output_file_best <- file.path(output_dir, paste0(id_i, "_ctmm_best_fit.rds"))
    
    saveRDS(model_selection, file = output_file_selection)
    saveRDS(best_model, file = output_file_best)
    
    best_models_list[[i]] <- best_model
    selection_list[[i]] <- model_selection
    gc()
  }
  
  total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf("\n=== Completed in %.1f minutes ===", total_elapsed))
  
  # Print model summary
  message("\n=== Model Selection Summary ===")
  selected_models <- sapply(best_models_list, function(x) summary(x)$name)
  model_table <- table(selected_models)
  message("Selected models across individuals:")
  print(model_table)
  
  # Save combined lists
  output_best_file <- file.path(output_dir, paste0(lake_name, "_", species_name, "_best_models.rds"))
  output_selection_file <- file.path(output_dir, paste0(lake_name, "_", species_name, "_all_selections.rds"))
  
  saveRDS(best_models_list, output_best_file)
  saveRDS(selection_list, output_selection_file)
  
  return(list(
    best_models = best_models_list,
    selection_results = selection_list
  ))
}

# Function to verify model fits ---------------------------------------------
verify_fits <- function(telem_list, fit_list, species_name, n_check = NULL) {
  
  message("\n=== Verifying ", species_name, " Model Fits ===")
  
  if (length(fit_list) == 0) {
    message("No fits to verify!")
    return(invisible())
  }
  
  # If n_check is NULL, check all individuals
  if (is.null(n_check)) {
    n_check <- length(fit_list)
  }
  
  # Ensure n_check doesn't exceed available fits
  n_check <- min(n_check, length(fit_list))
  
  message(sprintf("\nChecking %d out of %d individuals", n_check, length(fit_list)))
  
  # Loop through individuals to check
  for (i in 1:n_check) {
    id_i <- names(fit_list)[i]
    
    message("\n", strrep("-", 80))
    message(sprintf("Individual %d/%d: %s", i, n_check, id_i))
    message(strrep("-", 80))
    
    # Print summary
    print(summary(fit_list[[i]]))
    
    # Create diagnostic plot
    message("\nGenerating diagnostic plot...")
    plot(telem_list[[id_i]], fit_list[[i]], error = FALSE)
    
    # Pause for user to review (except for last one)
    if (i < n_check) {
      message("\nPress [Enter] to continue to next individual, or type 'q' to quit verification...")
      user_input <- readline()
      if (tolower(trimws(user_input)) == 'q') {
        message("Verification stopped by user.")
        break
      }
    }
  }
  
  message("\n", strrep("=", 80))
  message("Verification complete!")
  message(strrep("=", 80))
}

# Quick verification function (summaries only, no plots) -------------------
verify_fits_quick <- function(fit_list, species_name) {
  
  message("\n=== Quick Summary of ", species_name, " Model Fits ===")
  
  if (length(fit_list) == 0) {
    message("No fits to verify!")
    return(invisible())
  }
  
  # Create a summary dataframe
  summary_df <- data.frame(
    ID = names(fit_list),
    Model = character(length(fit_list)),
    AIC = numeric(length(fit_list)),
    DOF_area = numeric(length(fit_list)),
    DOF_speed = numeric(length(fit_list)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(fit_list)) {
    model_sum <- summary(fit_list[[i]])
    summary_df$Model[i] <- model_sum$name
    summary_df$AIC[i] <- model_sum$IC
    
    # Extract DOF information if available
    if ("DOF" %in% names(model_sum)) {
      summary_df$DOF_area[i] <- model_sum$DOF["area"]
      summary_df$DOF_speed[i] <- ifelse("speed" %in% names(model_sum$DOF), 
                                        model_sum$DOF["speed"], NA)
    }
  }
  
  print(summary_df)
  
  message("\nModel type distribution:")
  print(table(summary_df$Model))
  
  return(invisible(summary_df))
}


#=============================================================================-
# 3. PREPARE DATA FOR REFITTING ####
#=============================================================================-


# Load filtered data
BT_filt_data <- readRDS(paste0(filtered_data_path, "05_BT_sub.rds"))


# Convert to Movebank format for ctmm package ------------------------------
BT_movebank <- with(
  BT_filt_data,
  data.frame(
    "timestamp" = timestamp,
    "location.long" = Long,
    "location.lat" = Lat,
    "GPS.HDOP" = HDOP,
    "individual-local-identifier" = individual_ID,
    "species" = species,
    "weight" = weight,
    "total_length" = total_length,
    "std_length" = std_length,
    "treatment" = treatment,
    "date" = date,
    "exp_stage" = exp_stage,
    "time_of_day" = time_of_day,
    "found_alive" = found_alive,
    "known_predated" = known_predated
  )
)

# Free up memory ------------------------------------------------------------
rm(BT_filt_data)
gc()

# Convert to telemetry object -----------------------------------------------
message("\nConverting to ctmm telemetry object...")
BT_tels <- as.telemetry(
  BT_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,  # Will be set to geometric median automatically
  datum = "WGS84",
  keep = c("species", "weight", "total_length", "std_length", "treatment",
           "date", "exp_stage", "time_of_day", "found_alive", "known_predated")
)

message("Telemetry object created with ", length(BT_tels), " individuals")
message("Projection: ", ctmm::projection(BT_tels[[1]]))
message("Timezone: ", tz(BT_tels[[1]]$timestamp))


# Incorporate UERE error into telemetry objects ------------------------------
#load UERE
BT_UERE <- readRDS(paste0(save_telem_path, "BT/BT_UERE.rds"))
print(summary(BT_UERE))
uere(BT_tels) <- BT_UERE


#==============================================================================-
# 4. ORGANIZE INDIVIDUALS BY SPECIES ####
#==============================================================================-

# Verify individual order ---------------------------------------------------
species_order <- BT_movebank %>%
  dplyr::select(species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

print(table(species_order$species))

# Display the ordering to help with indexing
print(species_order)

# Split telemetry objects by species ---------------------------------------
pike_ids <- species_order %>% filter(species == "Northern Pike") %>% pull(individual.local.identifier)
perch_ids <- species_order %>% filter(species == "Perch") %>% pull(individual.local.identifier)
roach_ids <- species_order %>% filter(species == "Roach") %>% pull(individual.local.identifier)

pike_BT_tel <- BT_tels[names(BT_tels) %in% pike_ids]
perch_BT_tel <- BT_tels[names(BT_tels) %in% perch_ids]
roach_BT_tel <- BT_tels[names(BT_tels) %in% roach_ids]

message("\nSpecies groups created:")
message("Pike: ", length(pike_BT_tel), " individuals")
message("Perch: ", length(perch_BT_tel), " individuals")
message("Roach: ", length(roach_BT_tel), " individuals")

# Save species-specific telemetry objects -----------------------------------
saveRDS(pike_BT_tel, paste0(save_telem_path, "BT/pike_BT_tel_thinned_final.rds"))
saveRDS(perch_BT_tel, paste0(save_telem_path, "BT/perch_BT_tel_thinned_final.rds"))
saveRDS(roach_BT_tel, paste0(save_telem_path, "BT/roach_BT_tel_thinned_final.rds"))


#==============================================================================#
# 5. FIT CTMM MODELS - PIKE ####
#==============================================================================#

# Option to reload if needed ------------------------------------------------
pike_BT_tel <- readRDS(paste0(save_telem_path, "BT/pike_BT_tel_thinned_final.rds"))

if (length(pike_BT_tel) > 0) {
  # Fit models using parallel processing -----------------------------------
  pike_results <- fit_ctmm_species_parallel(
    pike_BT_tel, 
    "pike",
    lake_name = "BT",
    max_cores = 3,
    ic = "AICc"
  )
  
  BT_pike_best_models <- pike_results$best_models
  BT_pike_selections <- pike_results$selection_results
  
  # Verify fits -----------------------------------------------------------
  # To reload if needed:
  # BT_pike_best_models <- readRDS(paste0(save_ctmm_path, "BT_pike_fits/BT_pike_best_models.rds"))
  
  verify_fits(pike_BT_tel, BT_pike_best_models, "Pike")
} else {
  message("\n*** No pike individuals found in BT dataset ***")
  BT_pike_best_models <- list()
}

#==============================================================================-
# 6. FIT CTMM MODELS - PERCH ####
#==============================================================================-

# Option to reload if needed ------------------------------------------------
perch_BT_tel <- readRDS(paste0(save_telem_path, "BT/perch_BT_tel_thinned_final.rds"))

if (length(perch_BT_tel) > 0) {
  # Fit models using parallel processing -----------------------------------
  perch_cores <- ifelse(length(perch_BT_tel) > 15, 15, 3)
  
  perch_results <- fit_ctmm_species_parallel(
    perch_BT_tel,
    "perch",
    lake_name = "BT",
    max_cores = perch_cores,
    ic = "AICc"
  )
  
  BT_perch_best_models <- perch_results$best_models
  BT_perch_selections <- perch_results$selection_results
  
  # Verify fits -----------------------------------------------------------
  # To reload if needed:
  # BT_perch_best_models <- readRDS(paste0(save_ctmm_path, "BT_perch_fits/BT_perch_best_models.rds"))
  
  verify_fits(perch_BT_tel, BT_perch_best_models, "Perch")
} else {
  message("\n*** No perch individuals found in BT dataset ***")
  BT_perch_best_models <- list()
}

#==============================================================================-
# 6. FIT CTMM MODELS - ROACH ####
#==============================================================================-

# Option to reload if needed ------------------------------------------------
roach_BT_tel <- readRDS(paste0(save_telem_path, "BT/roach_BT_tel_thinned_final.rds"))

if (length(roach_BT_tel) > 0) {
  # Fit models using parallel processing -----------------------------------
  roach_cores <- ifelse(length(roach_BT_tel) > 15, 15, 3)
  
  roach_results <- fit_ctmm_species_parallel(
    roach_BT_tel,
    "roach",
    lake_name = "BT",
    max_cores = roach_cores,
    ic = "AICc"
  )
  
  BT_roach_best_models <- roach_results$best_models
  BT_roach_selections <- roach_results$selection_results
  
  # Verify fits -----------------------------------------------------------
  verify_fits(roach_BT_tel, BT_roach_best_models, "Roach")
} else {
  message("\n*** No roach individuals found in BT dataset ***")
  BT_roach_best_models <- list()
}

#=============================================================================-
# 7. FINAL SUMMARY ####
#=============================================================================-

message("\n", strrep("=", 80))
message("=== ALL CTMM MODEL SELECTIONS COMPLETED - LAKE BT ===")
message(strrep("=", 80))

message("\nBest models selected:")
message("  Pike: ", length(BT_pike_best_models), " individuals")
message("  Perch: ", length(BT_perch_best_models), " individuals")
message("  Roach: ", length(BT_roach_best_models), " individuals")
message("  Total: ", length(BT_pike_best_models) + 
          length(BT_perch_best_models) + 
          length(BT_roach_best_models), " individuals")

if (exists("pike_results") && length(pike_results$best_models) > 0) {
  message("\nTotal processing times:")
  message("  Pike: ", sprintf("%.1f", pike_results$diagnostics$total_time), " minutes")
  message("  Pike model distribution:")
  print(pike_results$diagnostics$model_counts)
}
if (exists("perch_results") && length(perch_results$best_models) > 0) {
  message("  Perch: ", sprintf("%.1f", perch_results$diagnostics$total_time), " minutes")
  message("  Perch model distribution:")
  print(perch_results$diagnostics$model_counts)
}
if (exists("roach_results") && length(roach_results$best_models) > 0) {
  message("  Roach: ", sprintf("%.1f", roach_results$diagnostics$total_time), " minutes")
  message("  Roach model distribution:")
  print(roach_results$diagnostics$model_counts)
}

message("\nOutput locations:")
message("  Individual selections: ", save_ctmm_path, "BT_[species]_fits/[ID]_ctmm_selection.rds")
message("  Individual best models: ", save_ctmm_path, "BT_[species]_fits/[ID]_ctmm_best_fit.rds")
message("  Combined best models: ", save_ctmm_path, "BT_[species]_fits/BT_[species]_best_models.rds")
message("  All selections: ", save_ctmm_path, "BT_[species]_fits/BT_[species]_all_selections.rds")

message("\nTo reload best models:")
message("  pike_models <- readRDS('", save_ctmm_path, "BT_pike_fits/BT_pike_best_models.rds')")
message("  perch_models <- readRDS('", save_ctmm_path, "BT_perch_fits/BT_perch_best_models.rds')")
message("  roach_models <- readRDS('", save_ctmm_path, "BT_roach_fits/BT_roach_best_models.rds')")

message("\n", strrep("=", 80))

#==============================================================================
# NOTES ON USAGE
#==============================================================================

# The updated workflow now:
# 1. Uses ctmm.select to compare multiple candidate models per individual
# 2. Saves both the full selection results and the best model for each individual
# 3. Returns a list of best ctmm objects that can be used directly for downstream analyses
# 4. Provides model selection summaries showing which models were most commonly selected
#
# Access the results:
# - pike_results$best_models: List of best ctmm objects for each pike
# - pike_results$selection_results: Full ctmm.select output for each pike (all candidate models)
# - pike_results$diagnostics$model_counts: Summary of selected model types
#
# If your PC gets too loud or hot during processing:
# 1. Reduce max_cores to 2-3 for all species
# 2. Or use the sequential version:
#    pike_results <- fit_ctmm_species_sequential(pike_BT_tel, "pike", "BT", ic = "AICc")

