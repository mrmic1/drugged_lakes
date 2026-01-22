# ===================================================================
# Predator-Prey Interaction Analysis - Muddyfoot
# ===================================================================
# 
# ANALYSIS OVERVIEW:
# This script identifies potential predation events by analyzing movement 
# tracking data for prey (roach/perch) and predators (pike). Enhanced to use
# ctmm::proximity() to detect non-independent movement patterns.
#
# >>> CHANGE: Now uses TIERED ENCOUNTER THRESHOLDS based on telemetry error
# Predation is inferred from:
#   1. Close proximity using tiered thresholds:
#      - High-confidence: ≤1.5m (within GPS error)
#      - Probable: 1.5-3m (with temporal confirmation)
#      - Possible: 3-5m (shared space use)
#   2. Sustained encounters over consecutive days
#   3. Movement pattern changes suggesting predation (proximity analysis)
#   4. Non-independent movement detected by ctmm::proximity()
#
# Author: Marcus Michelangeli
# ===================================================================

# ===================================================================
# 1. SETUP & CONFIGURATION
# ===================================================================

## 1.1 Load Required Libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ctmm)
  library(lubridate)
  library(data.table)
  library(sf)
  library(flextable)
  library(parallel)
  library(future.apply)
  library(progressr)
  library(ggplot2)
  library(openxlsx)
})

# Configure time zone
Sys.setenv(TZ = 'Europe/Stockholm')

## 1.2 Define Directory Paths ----
paths <- list(
  filtered_data = "./data/tracks_filtered/muddyfoot/",
  lake_polygon  = "./data/lake_coords/",
  ctmm          = "./data/ctmm_fits/",
  telem         = "./data/telem_obj/muddyfoot/",
  encounters    = "./data/encounters/muddyfoot/",
  proximity     = "./data/proximity/muddyfoot/",
  size          = "./data/fish_size/",
  tables        = "./tables/muddyfoot/",
  figures       = "./figures/muddyfoot/" 
)

# Create proximity directory if it doesn't exist
if (!dir.exists(paths$proximity)) {
  dir.create(paths$proximity, recursive = TRUE)
}

## 1.3 Define Analysis Parameters ----
# >>> CHANGE: Added tiered encounter thresholds based on GPS error analysis
params <- list(
  # >>> CHANGE: Tiered distance thresholds
  strike_distance = 0.45,              # Pike max strike distance (m)
  high_confidence_threshold = 1.5,     # ≤1.5m: High confidence (within GPS error ~1.4m)
  probable_threshold = 3.0,            # ≤3m: Probable encounter (with confirmation)
  possible_threshold = 5.0,            # ≤5m: Possible encounter (shared space)
  
  # Encounter count thresholds (kept for backward compatibility)
  encounter_threshold_low = 20,        # Minimum encounters for suspected predation
  encounter_threshold_high = 100,      # High encounter count threshold
  
  # >>> CHANGE: Updated to use high-confidence threshold
  min_distance_threshold = 1.5,        # Minimum distance for high-confidence encounter
  consecutive_days_threshold = 2,      # Consecutive days for likely predation
  proximity_encounter_threshold = 100  # Total encounters threshold for proximity analysis
)

# >>> CHANGE: Print encounter strategy
message("\n=== ENCOUNTER DETECTION STRATEGY ===")
message("Tiered approach based on GPS error (σ_combined ≈ 0.72m):")
message(sprintf("  High-confidence: ≤%.1fm (within 95%% GPS error ~1.4m)", 
                params$high_confidence_threshold))
message(sprintf("  Probable: %.1f-%.1fm (requires temporal confirmation)", 
                params$high_confidence_threshold, params$probable_threshold))
message(sprintf("  Possible: %.1f-%.1fm (shared space use)", 
                params$probable_threshold, params$possible_threshold))
message("=====================================\n")

# ===================================================================
# 2. CUSTOM FUNCTIONS
# ===================================================================

#' Calculate Pairwise Distances Between Prey and Predator
#' 
#' >>> CHANGE: Now flags encounters using tiered thresholds
#' 
#' @param prey_tel Telemetry object for prey
#' @param pred_tel Telemetry object for predator
#' @param prey_fits CTMM fits for prey
#' @param pred_fits CTMM fits for predator
#' @param prey_species Species name for prey (e.g., "Roach", "Perch")
#' @param strike_dist Maximum strike distance for encounter definition (deprecated, use thresholds)
#' @param high_conf_threshold High-confidence encounter threshold (m)
#' @param probable_threshold Probable encounter threshold (m)
#' @param possible_threshold Possible encounter threshold (m)
#' @param parallel Logical, use parallel processing?
#' @param n_cores Number of cores for parallel processing
#' 
#' @return Data frame with pairwise distances and tiered encounter flags
calculate_pairwise_distances <- function(prey_tel, pred_tel, 
                                         prey_fits, pred_fits,
                                         prey_species = "Prey",
                                         strike_dist = 0.45,  # Kept for compatibility
                                         high_conf_threshold = 1.5,
                                         probable_threshold = 3.0,
                                         possible_threshold = 5.0,
                                         parallel = TRUE,
                                         n_cores = NULL) {
  
  # Ensure matching projections
  ctmm::projection(prey_tel) <- ctmm::projection(pred_tel)
  ctmm::projection(prey_fits) <- ctmm::projection(pred_fits)
  
  stopifnot(projection(prey_tel) == projection(pred_tel))
  
  # Prepare combinations
  n_prey <- length(prey_tel)
  n_pred <- length(pred_tel)
  combinations <- expand.grid(prey_idx = 1:n_prey, pred_idx = 1:n_pred)
  
  message(sprintf("Calculating distances for %d %s x %d Pike = %d combinations",
                  n_prey, prey_species, n_pred, nrow(combinations)))
  
  # Define distance calculation function
  calc_dist <- function(idx) {
    i <- combinations$prey_idx[idx]
    j <- combinations$pred_idx[idx]
    
    tryCatch({
      combined_telemetry <- c(prey_tel[i], pred_tel[j])
      combined_ctmm <- c(prey_fits[i], pred_fits[j])
      
      # Calculate pairwise distances
      location_difference <- ctmm::distances(combined_telemetry, combined_ctmm)
      
      # Extract IDs
      prey_id <- names(prey_tel)[i]
      pred_id <- names(pred_tel)[j]
      
      # Convert to data frame
      df <- as.data.frame(location_difference)
      df$Prey_ID <- prey_id
      df$Pred_ID <- pred_id
      df$Species <- prey_species
      
      # >>> CHANGE: Add tiered encounter flags
      df$encounter <- ifelse(df$est <= strike_dist, 1, 0)  # Legacy flag
      df$encounter_high_conf <- ifelse(df$est <= high_conf_threshold, 1, 0)
      df$encounter_probable <- ifelse(df$est > high_conf_threshold & 
                                        df$est <= probable_threshold, 1, 0)
      df$encounter_possible <- ifelse(df$est > probable_threshold & 
                                        df$est <= possible_threshold, 1, 0)
      
      # >>> CHANGE: Add categorical encounter type
      df$encounter_type <- case_when(
        df$est <= high_conf_threshold ~ "high_confidence",
        df$est <= probable_threshold ~ "probable",
        df$est <= possible_threshold ~ "possible",
        TRUE ~ "no_encounter"
      )
      
      # Add date column
      if ("timestamp" %in% colnames(df)) {
        df$timestamp <- as.POSIXct(df$timestamp, tz = "Europe/Stockholm")
        df$Date <- as.Date(df$timestamp, tz = "Europe/Stockholm")
      }
      
      return(df)
    }, error = function(e) {
      warning(sprintf("Error for %s-%s: %s", prey_id, pred_id, e$message))
      return(NULL)
    })
  }
  
  # Execute calculations (parallel or sequential)
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
    message(sprintf("Using %d cores for parallel processing", n_cores))
    
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(cl, c("prey_tel", "pred_tel", "prey_fits", 
                                  "pred_fits", "strike_dist", "prey_species",
                                  "high_conf_threshold", "probable_threshold",
                                  "possible_threshold", "combinations"),
                            envir = environment())
    
    parallel::clusterEvalQ(cl, {
      library(ctmm)
      library(dplyr)
    })
    
    results <- parallel::parLapply(cl, 1:nrow(combinations), calc_dist)
  } else {
    results <- lapply(1:nrow(combinations), function(idx) {
      if (idx %% 10 == 0) message(sprintf("Processing %d/%d", idx, nrow(combinations)))
      calc_dist(idx)
    })
  }
  
  # Combine results
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) {
    stop("No valid distance calculations completed")
  }
  
  distances_df <- data.table::rbindlist(results)
  
  message(sprintf("Completed: %d distance calculations", nrow(distances_df)))
  
  return(as.data.frame(distances_df))
}


#' Summarize Daily Encounters with Tiered Approach
#' 
#' >>> CHANGE: Now calculates encounter counts for each tier
#' 
#' @param distances_df Data frame with distance calculations
#' @return Daily encounter summary with tiered counts
summarize_daily_encounters <- function(distances_df) {
  daily_summary <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      Species = first(Species),
      
      # >>> CHANGE: Tiered encounter counts
      encounter_count_high_conf = sum(encounter_high_conf, na.rm = TRUE),
      encounter_count_probable = sum(encounter_probable, na.rm = TRUE),
      encounter_count_possible = sum(encounter_possible, na.rm = TRUE),
      encounter_count_legacy = sum(encounter, na.rm = TRUE),  # Keep for compatibility
      
      daily_avg_dist = mean(est, na.rm = TRUE),
      daily_min_dist = min(est, na.rm = TRUE),
      
      # >>> CHANGE: Flag if there are consecutive probable encounters (temporal confirmation)
      has_consecutive_probable = any(encounter_probable == 1) & 
        any(lag(encounter_probable, default = 0) == 1),
      
      .groups = "drop"
    )
  
  return(daily_summary)
}


#' Calculate Total Encounters Per Prey-Predator Pair
#' 
#' >>> CHANGE: Now includes tiered encounter totals
#' 
#' @param distances_df Data frame with distance calculations
#' @return Summary with tiered total encounters
calculate_total_encounters <- function(distances_df) {
  total_encounters <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Species) %>%
    summarise(
      # >>> CHANGE: Tiered total encounters
      total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
      total_encounters_probable = sum(encounter_probable, na.rm = TRUE),
      total_encounters_possible = sum(encounter_possible, na.rm = TRUE),
      total_encounters_legacy = sum(encounter, na.rm = TRUE),  # Keep for compatibility
      
      min_distance = min(est, na.rm = TRUE),
      mean_distance = mean(est, na.rm = TRUE),
      n_days = n_distinct(Date),
      
      # >>> CHANGE: Count days with different encounter types
      n_days_high_conf = n_distinct(Date[encounter_high_conf == 1]),
      n_days_probable = n_distinct(Date[encounter_probable == 1]),
      n_days_possible = n_distinct(Date[encounter_possible == 1]),
      
      .groups = "drop"
    ) %>%
    arrange(desc(total_encounters_high_conf))
  
  return(total_encounters)
}


#' Perform Proximity Analysis on Prey-Predator Pair
#' 
#' @param prey_tel Telemetry object for single prey
#' @param pred_tel Telemetry object for single predator
#' @param prey_fit CTMM fit for prey
#' @param pred_fit CTMM fit for predator
#' @param prey_id ID of prey individual
#' @param pred_id ID of predator individual
#' 
#' @return List containing proximity results and statistics
perform_proximity_analysis <- function(prey_tel, pred_tel, 
                                       prey_fit, pred_fit,
                                       prey_id, pred_id) {
  
  tryCatch({
    # Ensure matching projections
    projection(prey_tel) <- projection(pred_tel)
    projection(prey_fit) <- projection(pred_fit)
    
    # Combine telemetry and fits
    combined_tel <- c(prey_tel, pred_tel)
    combined_fits <- c(prey_fit, pred_fit)
    
    # Run proximity analysis
    message(sprintf("  Running proximity for %s - %s", prey_id, pred_id))
    proximity_result <- ctmm::proximity(combined_tel, combined_fits)
    
    # Extract key statistics
    proximity_stats <- list(
      prey_id = prey_id,
      pred_id = pred_id,
      z_statistic = proximity_result$CI[1, "est"],  # Z-statistic
      z_lower = proximity_result$CI[1, "low"],
      z_upper = proximity_result$CI[1, "high"],
      p_value = proximity_result$P,  # P-value for independence test
      independent = proximity_result$P > 0.05,  # Movement is independent if p > 0.05
      proximity_object = proximity_result
    )
    
    return(proximity_stats)
    
  }, error = function(e) {
    warning(sprintf("Proximity analysis failed for %s - %s: %s", 
                    prey_id, pred_id, e$message))
    return(list(
      prey_id = prey_id,
      pred_id = pred_id,
      z_statistic = NA,
      z_lower = NA,
      z_upper = NA,
      p_value = NA,
      independent = NA,
      error = e$message
    ))
  })
}


#' Calculate Maximum Consecutive Days
#' 
#' @param dates Vector of dates
#' @return Maximum number of consecutive days
calc_max_consecutive <- function(dates) {
  if (length(dates) == 0) return(0)
  sorted_dates <- sort(dates)
  consecutive <- rle(c(1, diff(sorted_dates) == 1))$lengths
  return(max(consecutive, na.rm = TRUE))
}


#' Identify Suspected Predation Events with Tiered Approach
#' 
#' >>> CHANGE: Complete rewrite to use tiered thresholds
#' 
#' @param distances_df Data frame with pairwise distances
#' @param enc_threshold_low Minimum high-confidence encounters per day
#' @param enc_threshold_high High encounter threshold
#' @param min_dist_threshold High-confidence distance threshold (m)
#' @param consec_days_threshold Consecutive days for likely predation
#' 
#' @return Data frame with suspected predation events using tiered classification
identify_predation_events <- function(distances_df,
                                      enc_threshold_low = 10,  # >>> CHANGE: Lower threshold for high-confidence
                                      enc_threshold_high = 50,  # >>> CHANGE: Adjusted for high-confidence
                                      min_dist_threshold = 1.5,  # >>> CHANGE: High-confidence threshold
                                      consec_days_threshold = 2) {
  
  # >>> CHANGE: Daily summary with tiered encounters
  daily_encounters <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      # Tiered encounter counts
      encounter_count_high_conf = sum(encounter_high_conf, na.rm = TRUE),
      encounter_count_probable = sum(encounter_probable, na.rm = TRUE),
      encounter_count_probable_confirmed = sum(encounter_probable == 1 & 
                                                 lag(encounter_probable, default = 0) == 1, 
                                               na.rm = TRUE),
      daily_min_dist = min(est, na.rm = TRUE),
      .groups = "drop"
    )
  
  # >>> CHANGE: Filter for potential predation days using tiered approach
  high_risk_days <- daily_encounters %>%
    filter(
      encounter_count_high_conf >= enc_threshold_low |  # High-confidence encounters
        encounter_count_probable_confirmed >= enc_threshold_low |  # Confirmed probable
        daily_min_dist < min_dist_threshold  # Very close approach
    )
  
  # >>> CHANGE: Summarize predation patterns with tiered metrics
  predation_summary <- high_risk_days %>%
    group_by(Prey_ID, Pred_ID) %>%
    summarise(
      # High-confidence metrics
      num_days_high_conf_10 = sum(encounter_count_high_conf >= 10),
      num_days_high_conf_50 = sum(encounter_count_high_conf >= 50),
      
      # Probable encounter metrics (with temporal confirmation)
      num_days_probable_confirmed_10 = sum(encounter_count_probable_confirmed >= 10),
      
      # Distance-based metrics
      num_days_min_dist_less_1.5m = sum(daily_min_dist < 1.5),
      num_days_min_dist_less_0.5m = sum(daily_min_dist < 0.5),
      
      # First occurrence dates
      first_date_high_conf = if(any(encounter_count_high_conf >= enc_threshold_low)) 
        min(Date[encounter_count_high_conf >= enc_threshold_low]) else as.Date(NA),
      first_date_probable = if(any(encounter_count_probable_confirmed >= enc_threshold_low)) 
        min(Date[encounter_count_probable_confirmed >= enc_threshold_low]) else as.Date(NA),
      
      # Consecutive days patterns
      consecutive_days_high_conf = calc_max_consecutive(
        Date[encounter_count_high_conf >= enc_threshold_low]
      ),
      consecutive_days_probable = calc_max_consecutive(
        Date[encounter_count_probable_confirmed >= enc_threshold_low]
      ),
      
      # >>> CHANGE: Store encounter dates by type
      encounter_dates_high_conf = paste(
        format(Date[encounter_count_high_conf >= 10], "%d/%m/%Y"), 
        collapse = ", "
      ),
      encounter_dates_probable = paste(
        format(Date[encounter_count_probable_confirmed >= 10], "%d/%m/%Y"), 
        collapse = ", "
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      # >>> CHANGE: New tiered predation likelihood classification
      likely_predated = case_when(
        consecutive_days_high_conf >= consec_days_threshold ~ 1,  # Strong evidence
        num_days_high_conf_50 >= 2 ~ 1,  # Multiple days with many encounters
        num_days_min_dist_less_0.5m >= 2 ~ 1,  # Very close approaches
        TRUE ~ 0
      )
    )
  
  return(predation_summary)
}


# ===================================================================
# 3. DATA LOADING
# ===================================================================

## 3.1 Load Telemetry Objects ----
pike_tel <- readRDS(file.path(paths$telem, 'pike_muddyfoot_tel_thinned.rds'))
perch_tel <- readRDS(file.path(paths$telem, 'perch_muddyfoot_tel_thinned.rds'))
roach_tel <- readRDS(file.path(paths$telem, 'roach_muddyfoot_tel_thinned.rds'))

message(sprintf("Loaded telemetry: %d pike, %d perch, %d roach", 
                length(pike_tel), length(perch_tel), length(roach_tel)))

## 3.2 Load CTMM Fits ----
pike_fits <- readRDS(file.path(paths$ctmm, "muddyfoot_pike_fits/muddyfoot_pike_ctmm_fits.rds"))
perch_fits <- readRDS(file.path(paths$ctmm, "muddyfoot_perch_fits/muddyfoot_perch_ctmm_fits.rds"))
roach_fits <- readRDS(file.path(paths$ctmm, "muddyfoot_roach_fits/muddyfoot_roach_ctmm_fits.rds"))

## 3.3 Load Metadata ----
muddyfoot_filt_data <- readRDS(file.path(paths$filtered_data, "03_muddyfoot_sub.rds"))
post_biometrics <- fread(file.path(paths$size, "biometric_post_exp_data.csv")) %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number)))

# Create post-biometric columns for filtering
post_biometric_cols <- post_biometrics %>%
  filter(Lake == 'Muddyfoot', Species %in% c('Roach', 'Perch')) %>%
  select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found, found_predated = Known_predated)


# ===================================================================
# 4. DISTANCE CALCULATIONS
# ===================================================================

message("\n=== DISTANCE CALCULATIONS WITH TIERED THRESHOLDS ===")

## 4.1 Roach-Pike Distances ----
roach_pike_file <- file.path(paths$encounters, "muddyfoot_pike_roach_distances_tiered_df.rds")

if (file.exists(roach_pike_file)) {
  message("Loading existing Roach-Pike distances...")
  roach_pike_distances_df <- readRDS(roach_pike_file)
} else {
  message("Calculating Roach-Pike distances with tiered thresholds...")
  roach_pike_distances_df <- calculate_pairwise_distances(
    prey_tel = roach_tel,
    pred_tel = pike_tel,
    prey_fits = roach_fits,
    pred_fits = pike_fits,
    prey_species = "Roach",
    strike_dist = params$strike_distance,
    # >>> CHANGE: Pass tiered thresholds
    high_conf_threshold = params$high_confidence_threshold,
    probable_threshold = params$probable_threshold,
    possible_threshold = params$possible_threshold,
    parallel = TRUE,
    n_cores = 16
  )
  saveRDS(roach_pike_distances_df, roach_pike_file)
}

## 4.2 Perch-Pike Distances ----
perch_pike_file <- file.path(paths$encounters, "muddyfoot_pike_perch_distances_tiered_df.rds")

if (file.exists(perch_pike_file)) {
  message("Loading existing Perch-Pike distances...")
  perch_pike_distances_df <- readRDS(perch_pike_file)
} else {
  message("Calculating Perch-Pike distances with tiered thresholds...")
  perch_pike_distances_df <- calculate_pairwise_distances(
    prey_tel = perch_tel,
    pred_tel = pike_tel,
    prey_fits = perch_fits,
    pred_fits = pike_fits,
    prey_species = "Perch",
    strike_dist = params$strike_distance,
    # >>> CHANGE: Pass tiered thresholds
    high_conf_threshold = params$high_confidence_threshold,
    probable_threshold = params$probable_threshold,
    possible_threshold = params$possible_threshold,
    parallel = TRUE,
    n_cores = 16
  )
  saveRDS(perch_pike_distances_df, perch_pike_file)
}

message(sprintf("Loaded distances: %d roach-pike, %d perch-pike observations",
                nrow(roach_pike_distances_df), nrow(perch_pike_distances_df)))

# >>> CHANGE: Print summary of encounter types
message("\n=== ENCOUNTER TYPE DISTRIBUTION ===")
combined_distances <- bind_rows(roach_pike_distances_df, perch_pike_distances_df)
encounter_summary <- combined_distances %>%
  group_by(Species, encounter_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Species) %>%
  mutate(pct = 100 * n / sum(n))

print(encounter_summary)


# ===================================================================
# 5. IDENTIFY CANDIDATE PAIRS FOR PROXIMITY ANALYSIS
# ===================================================================

message("\n=== IDENTIFYING CANDIDATES FOR PROXIMITY ANALYSIS ===")

## 5.1 Calculate Total Encounters (with tiered metrics) ----
roach_total <- calculate_total_encounters(roach_pike_distances_df)
perch_total <- calculate_total_encounters(perch_pike_distances_df)

all_encounters <- bind_rows(roach_total, perch_total)

## 5.2 Filter for High-Risk Pairs ----
# >>> CHANGE: Use high-confidence encounters for proximity analysis selection
candidate_pairs <- all_encounters %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
  filter(
    found_alive == 0,  # Not recovered alive
    # >>> CHANGE: Use high-confidence or probable encounters
    (total_encounters_high_conf > 50 | 
       total_encounters_probable > params$proximity_encounter_threshold)
  ) %>%
  arrange(desc(total_encounters_high_conf))

message(sprintf("Identified %d candidate prey-predator pairs for proximity analysis",
                nrow(candidate_pairs)))

# >>> CHANGE: Print breakdown by encounter type
message("\nCandidate pair encounter breakdown:")
message(sprintf("  Mean high-confidence encounters: %.1f", 
                mean(candidate_pairs$total_encounters_high_conf)))
message(sprintf("  Mean probable encounters: %.1f", 
                mean(candidate_pairs$total_encounters_probable)))

# Save candidate pairs
saveRDS(candidate_pairs, 
        file.path(paths$proximity, "muddyfoot_proximity_candidate_pairs_tiered.rds"))


# ===================================================================
# 6. PROXIMITY ANALYSIS - TEST EXAMPLE
# ===================================================================

message("\n=== TESTING PROXIMITY ANALYSIS (TOP 3 PAIRS) ===")

## 6.1 Select Test Cases ----
test_pairs <- candidate_pairs %>%
  slice_head(n = 3)

message("Test pairs:")
print(test_pairs %>% select(Prey_ID, Pred_ID, Species, 
                            total_encounters_high_conf, 
                            total_encounters_probable))

## 6.2 Run Proximity Analysis on Test Cases ----
test_proximity_results <- list()

for (i in 1:nrow(test_pairs)) {
  prey_id <- test_pairs$Prey_ID[i]
  pred_id <- test_pairs$Pred_ID[i]
  species <- test_pairs$Species[i]
  
  message(sprintf("\nTest %d/%d: %s (%s) vs %s", 
                  i, nrow(test_pairs), prey_id, species, pred_id))
  
  # Get appropriate telemetry objects
  if (species == "Roach") {
    prey_tel_single <- roach_tel[[prey_id]]
    prey_fit_single <- roach_fits[[prey_id]]
  } else {
    prey_tel_single <- perch_tel[[prey_id]]
    prey_fit_single <- perch_fits[[prey_id]]
  }
  
  pred_tel_single <- pike_tel[[pred_id]]
  pred_fit_single <- pike_fits[[pred_id]]
  
  # Run proximity analysis
  prox_result <- perform_proximity_analysis(
    prey_tel_single, pred_tel_single,
    prey_fit_single, pred_fit_single,
    prey_id, pred_id
  )
  
  test_proximity_results[[i]] <- prox_result
  
  # Print results
  if (!is.na(prox_result$p_value)) {
    message(sprintf("  Z-statistic: %.3f [%.3f, %.3f]",
                    prox_result$z_statistic,
                    prox_result$z_lower,
                    prox_result$z_upper))
    message(sprintf("  P-value: %.4f", prox_result$p_value))
    message(sprintf("  Movement independent: %s", prox_result$independent))
    message(sprintf("  Interpretation: %s",
                    ifelse(prox_result$independent,
                           "No evidence of predation (movements independent)",
                           "POTENTIAL PREDATION (movements non-independent)")))
  }
}

# Convert test results to data frame
test_proximity_df <- do.call(rbind, lapply(test_proximity_results, function(x) {
  data.frame(
    Prey_ID = x$prey_id,
    Pred_ID = x$pred_id,
    z_statistic = x$z_statistic,
    z_lower = x$z_lower,
    z_upper = x$z_upper,
    p_value = x$p_value,
    independent = x$independent,
    stringsAsFactors = FALSE
  )
}))

# Save test results
saveRDS(test_proximity_results, 
        file.path(paths$proximity, "test_proximity_results_detailed_tiered.rds"))
write.csv(test_proximity_df, 
          file.path(paths$proximity, "test_proximity_results_summary_tiered.csv"),
          row.names = FALSE)

message("\nTest proximity analysis complete!")


# ===================================================================
# 7. FULL PROXIMITY ANALYSIS (ALL CANDIDATE PAIRS)
# ===================================================================

message("\n=== RUNNING FULL PROXIMITY ANALYSIS ===")
message(sprintf("Analyzing %d candidate pairs", nrow(candidate_pairs)))

all_proximity_results <- list()

for (i in 1:nrow(candidate_pairs)) {
  prey_id <- candidate_pairs$Prey_ID[i]
  pred_id <- candidate_pairs$Pred_ID[i]
  species <- candidate_pairs$Species[i]
  
  if (i %% 10 == 0) {
    message(sprintf("Progress: %d/%d (%.1f%%)", 
                    i, nrow(candidate_pairs), 100*i/nrow(candidate_pairs)))
  }
  
  # Get telemetry objects
  if (species == "Roach") {
    prey_tel_single <- roach_tel[[prey_id]]
    prey_fit_single <- roach_fits[[prey_id]]
  } else {
    prey_tel_single <- perch_tel[[prey_id]]
    prey_fit_single <- perch_fits[[prey_id]]
  }
  
  pred_tel_single <- pike_tel[[pred_id]]
  pred_fit_single <- pike_fits[[pred_id]]
  
  # Run proximity analysis
  prox_result <- perform_proximity_analysis(
    prey_tel_single, pred_tel_single,
    prey_fit_single, pred_fit_single,
    prey_id, pred_id
  )
  
  all_proximity_results[[i]] <- prox_result
}

## 7.2 Compile Full Results ----
full_proximity_df <- do.call(rbind, lapply(all_proximity_results, function(x) {
  data.frame(
    Prey_ID = x$prey_id,
    Pred_ID = x$pred_id,
    z_statistic = x$z_statistic,
    z_lower = x$z_lower,
    z_upper = x$z_upper,
    p_value = x$p_value,
    independent = x$independent,
    stringsAsFactors = FALSE
  )
}))

# Merge with encounter data
proximity_with_encounters <- full_proximity_df %>%
  left_join(candidate_pairs, by = c("Prey_ID", "Pred_ID"))

# Save full results
saveRDS(all_proximity_results, 
        file.path(paths$proximity, "full_proximity_results_detailed_tiered.rds"))
saveRDS(proximity_with_encounters,
        file.path(paths$proximity, "full_proximity_results_with_encounters_tiered.rds"))
write.csv(proximity_with_encounters,
          file.path(paths$proximity, "full_proximity_results_tiered.csv"),
          row.names = FALSE)

message("\nFull proximity analysis complete!")


# ===================================================================
# 8. INTEGRATE PROXIMITY RESULTS WITH PREDATION CLASSIFICATION
# ===================================================================

message("\n=== INTEGRATING PROXIMITY WITH PREDATION EVENTS ===")

## 8.1 Daily Encounter Summaries ----
roach_daily <- summarize_daily_encounters(roach_pike_distances_df) %>%
  rename(Pike_ID = Pred_ID)

perch_daily <- summarize_daily_encounters(perch_pike_distances_df) %>%
  rename(Pike_ID = Pred_ID)

prey_daily_encounters <- bind_rows(roach_daily, perch_daily)

## 8.2 Identify Suspected Predation (Distance-Based with Tiered Approach) ----
roach_predation <- identify_predation_events(
  roach_pike_distances_df,
  # >>> CHANGE: Updated thresholds for high-confidence encounters
  enc_threshold_low = 10,
  enc_threshold_high = 50,
  min_dist_threshold = params$high_confidence_threshold,
  consec_days_threshold = params$consecutive_days_threshold
) %>% rename(Pike_ID = Pred_ID) %>% mutate(Species = "Roach")

perch_predation <- identify_predation_events(
  perch_pike_distances_df,
  enc_threshold_low = 10,
  enc_threshold_high = 50,
  min_dist_threshold = params$high_confidence_threshold,
  consec_days_threshold = params$consecutive_days_threshold
) %>% rename(Pike_ID = Pred_ID) %>% mutate(Species = "Perch")

suspected_predation_distance <- bind_rows(roach_predation, perch_predation)

## 8.3 Merge with Proximity Results ----
# >>> CHANGE: Enhanced classification using tiered encounter data
suspected_predation_enhanced <- suspected_predation_distance %>%
  left_join(
    proximity_with_encounters %>% 
      select(Prey_ID, Pred_ID, z_statistic, p_value, independent,
             total_encounters_high_conf, total_encounters_probable),
    by = c("Prey_ID", "Pike_ID" = "Pred_ID")
  ) %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
  filter(found_alive == 0) %>%
  mutate(
    # Proximity suggests predation
    proximity_suggests_predation = !independent & !is.na(independent),
    
    # >>> CHANGE: Enhanced scoring system using tiered encounters
    combined_predation_score = case_when(
      # Very strong evidence: distance + proximity + high-confidence encounters
      likely_predated == 1 & proximity_suggests_predation & 
        num_days_high_conf_50 >= 2 ~ 4,
      
      # Strong evidence: distance + proximity OR many high-confidence encounters
      likely_predated == 1 & proximity_suggests_predation ~ 3,
      likely_predated == 1 & num_days_high_conf_50 >= 3 ~ 3,
      
      # Moderate evidence: one strong signal
      likely_predated == 1 ~ 2,
      proximity_suggests_predation ~ 2,
      num_days_high_conf_10 >= 3 ~ 2,
      
      # Weak evidence
      TRUE ~ 1
    ),
    
    # >>> CHANGE: Updated classification with more nuance
    classification = case_when(
      combined_predation_score == 4 ~ "Very high confidence predation",
      combined_predation_score == 3 ~ "High confidence predation",
      combined_predation_score == 2 ~ "Moderate confidence predation",
      TRUE ~ "Low confidence / uncertain"
    ),
    
    # >>> CHANGE: Add encounter tier summary
    primary_encounter_tier = case_when(
      num_days_high_conf_50 >= 2 ~ "High-confidence (≤1.5m)",
      num_days_high_conf_10 >= 2 ~ "High-confidence (≤1.5m)",
      num_days_probable_confirmed_10 >= 2 ~ "Probable (1.5-3m, confirmed)",
      TRUE ~ "Uncertain"
    )
  )

# Save enhanced results
saveRDS(suspected_predation_enhanced,
        file.path(paths$encounters, "muddyfoot_suspected_predation_enhanced_tiered.rds"))

message(sprintf("Enhanced predation classification complete: %d events analyzed",
                nrow(suspected_predation_enhanced)))

# >>> CHANGE: Updated summary by classification and tier
classification_summary <- suspected_predation_enhanced %>%
  group_by(Species, classification, primary_encounter_tier) %>%
  summarise(
    n_events = n(),
    mean_high_conf_encounters = mean(num_days_high_conf_10, na.rm = TRUE),
    mean_probable_encounters = mean(num_days_probable_confirmed_10, na.rm = TRUE),
    mean_p_value = mean(p_value, na.rm = TRUE),
    .groups = "drop"
  )

print("Classification Summary (Tiered Approach):")
print(classification_summary)


# ===================================================================
# 9. VISUALIZATION OF TIERED PROXIMITY RESULTS
# ===================================================================

message("\n=== CREATING VISUALIZATIONS ===")

## >>> CHANGE: 9.1 Encounter Type Distribution by Species ----
p0 <- ggplot(combined_distances %>% filter(encounter_type != "no_encounter"), 
             aes(x = encounter_type, fill = Species)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Distribution of Encounter Types",
    subtitle = "High-confidence (≤1.5m), Probable (1.5-3m), Possible (3-5m)",
    x = "Encounter Type",
    y = "Count",
    fill = "Prey Species"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(paths$figures, "encounter_type_distribution_tiered.png"),
       p0, width = 10, height = 6, dpi = 300)

## 9.2 Proximity P-value Distribution ----
p1 <- ggplot(proximity_with_encounters, aes(x = p_value)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Distribution of Proximity P-values",
    subtitle = "Red line = 0.05 significance threshold",
    x = "P-value (independence test)",
    y = "Count"
  ) +
  theme_classic()

ggsave(file.path(paths$proximity, "proximity_pvalue_distribution_tiered.png"),
       p1, width = 8, height = 6, dpi = 300)

## >>> CHANGE: 9.3 Z-statistic vs High-Confidence Encounters ----
p2 <- ggplot(proximity_with_encounters, 
             aes(x = total_encounters_high_conf, y = z_statistic)) +
  geom_point(aes(color = independent), size = 3, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("TRUE" = "darkgreen", "FALSE" = "darkred"),
    labels = c("TRUE" = "Independent", "FALSE" = "Non-independent"),
    name = "Movement Pattern"
  ) +
  labs(
    title = "Proximity Analysis: Z-statistic vs High-Confidence Encounters",
    subtitle = "High-confidence encounters: ≤1.5m (within GPS error)",
    x = "Total High-Confidence Encounters",
    y = "Z-statistic"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(paths$proximity, "proximity_z_vs_high_conf_encounters_tiered.png"),
       p2, width = 10, height = 7, dpi = 300)

## >>> CHANGE: 9.4 Classification Comparison with Encounter Tiers ----
p3 <- ggplot(suspected_predation_enhanced, 
             aes(x = classification, fill = primary_encounter_tier)) +
  geom_bar(position = "stack", alpha = 0.8) +
  facet_wrap(~Species) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(
    title = "Predation Events by Confidence Level and Encounter Tier",
    x = "Classification",
    y = "Number of Events",
    fill = "Primary Encounter Type"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(paths$proximity, "predation_classification_by_tier.png"),
       p3, width = 12, height = 7, dpi = 300)

## >>> CHANGE: 9.5 Distance Distribution by Encounter Type ----
p4 <- ggplot(combined_distances %>% 
               filter(encounter_type != "no_encounter") %>%
               sample_n(min(10000, n())),  # Sample for performance
             aes(x = est, fill = encounter_type)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = c(1.5, 3.0), linetype = "dashed", color = "red") +
  scale_fill_manual(
    values = c("high_confidence" = "#d73027", 
               "probable" = "#fee090", 
               "possible" = "#91bfdb"),
    labels = c("High-confidence (≤1.5m)", 
               "Probable (1.5-3m)", 
               "Possible (3-5m)")
  ) +
  labs(
    title = "Distance Distribution by Encounter Type",
    subtitle = "Vertical lines show tier thresholds (1.5m, 3m)",
    x = "Distance (m)",
    y = "Count",
    fill = "Encounter Type"
  ) +
  xlim(0, 5) +
  theme_classic()

ggsave(file.path(paths$figures, "distance_distribution_by_tier.png"),
       p4, width = 10, height = 6, dpi = 300)


# ===================================================================
# 10. EXPORT COMPREHENSIVE RESULTS
# ===================================================================

message("\n=== EXPORTING RESULTS ===")

## 10.1 Create Excel Workbook ----
wb <- createWorkbook()

# >>> CHANGE: Add tiered encounter summary sheet
addWorksheet(wb, "Tiered Encounter Summary")
tier_summary <- combined_distances %>%
  group_by(Species, encounter_type) %>%
  summarise(
    n_observations = n(),
    median_distance = median(est),
    mean_distance = mean(est),
    min_distance = min(est),
    max_distance = max(est),
    .groups = "drop"
  )
writeData(wb, "Tiered Encounter Summary", tier_summary)

# Proximity test results
addWorksheet(wb, "Test Proximity")
writeData(wb, "Test Proximity", test_proximity_df)

# Full proximity results
addWorksheet(wb, "Full Proximity")
writeData(wb, "Full Proximity", proximity_with_encounters)

# Enhanced predation classification
addWorksheet(wb, "Enhanced Predation")
writeData(wb, "Enhanced Predation", suspected_predation_enhanced)

# Classification summary
addWorksheet(wb, "Classification Summary")
writeData(wb, "Classification Summary", classification_summary)

# Candidate pairs
addWorksheet(wb, "Candidate Pairs")
writeData(wb, "Candidate Pairs", candidate_pairs)

# >>> CHANGE: Add daily encounters with tiers
addWorksheet(wb, "Daily Encounters Tiered")
writeData(wb, "Daily Encounters Tiered", prey_daily_encounters)

# Save workbook
saveWorkbook(wb, 
             file.path(paths$tables, "muddyfoot_proximity_analysis_tiered_complete.xlsx"),
             overwrite = TRUE)

message("Excel workbook saved successfully!")


# ===================================================================
# 11. SUMMARY REPORT WITH TIERED METRICS
# ===================================================================

message("\n" %+% paste(rep("=", 70), collapse = ""))
message("ANALYSIS COMPLETE - TIERED APPROACH SUMMARY")
message(paste(rep("=", 70), collapse = ""))

# >>> CHANGE: Report tiered encounter statistics
total_observations <- nrow(combined_distances)
message(sprintf("\nTotal observations: %,d", total_observations))

message("\nEncounter distribution:")
enc_dist <- table(combined_distances$encounter_type)
for (type in names(enc_dist)) {
  pct <- 100 * enc_dist[type] / total_observations
  message(sprintf("  %s: %,d (%.2f%%)", type, enc_dist[type], pct))
}

message(sprintf("\nTotal prey-predator pairs analyzed: %d", nrow(candidate_pairs)))
message(sprintf("Pairs with non-independent movement (p < 0.05): %d", 
                sum(!proximity_with_encounters$independent, na.rm = TRUE)))
message(sprintf("Percentage showing non-independence: %.1f%%",
                100 * mean(!proximity_with_encounters$independent, na.rm = TRUE)))

message("\nPredation confidence classification:")
print(table(suspected_predation_enhanced$classification))

message("\nPrimary encounter tier for suspected predation:")
print(table(suspected_predation_enhanced$primary_encounter_tier))

message("\nOutput locations:")
message(sprintf("  Proximity results: %s", paths$proximity))
message(sprintf("  Encounter data: %s", paths$encounters))
message(sprintf("  Tables: %s", paths$tables))
message(sprintf("  Figures: %s", paths$figures))

message("\n" %+% paste(rep("=", 70), collapse = ""))


# ===================================================================
# 12. SESSION INFO
# ===================================================================

sink(file.path(paths$tables, "session_info_proximity_tiered.txt"))
cat("Proximity Analysis Session Info (Tiered Approach)\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")
cat("Analysis Parameters:\n")
cat(sprintf("  High-confidence threshold: %.1fm\n", params$high_confidence_threshold))
cat(sprintf("  Probable threshold: %.1f-%.1fm\n", 
            params$high_confidence_threshold, params$probable_threshold))
cat(sprintf("  Possible threshold: %.1f-%.1fm\n", 
            params$probable_threshold, params$possible_threshold))
cat("\n")
print(sessionInfo())
sink()

message("\nAll analyses complete! Check output directories for results.")
message("Key files to review:")
message("  - muddyfoot_proximity_analysis_tiered_complete.xlsx")
message("  - encounter_type_distribution_tiered.png")
message("  - predation_classification_by_tier.png")