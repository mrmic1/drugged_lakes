# ===================================================================
# Predator-Prey Interaction Analysis - Muddyfoot
# ===================================================================
# 
# ANALYSIS OVERVIEW:
# This script identifies potential predation events by analyzing movement 
# tracking data for prey (roach/perch) and predators (pike). Enhanced to use
# ctmm::proximity() to detect non-independent movement patterns.
#
# Predation is inferred from:
#   1. Close proximity (<0.45m - pike strike distance)
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
  proximity     = "./data/proximity/muddyfoot/",  # New directory for proximity results
  size          = "./data/fish_size/",
  tables        = "./tables/muddyfoot/",
  figures       = "./figures/muddyfoot/" 
)

# Create proximity directory if it doesn't exist
if (!dir.exists(paths$proximity)) {
  dir.create(paths$proximity, recursive = TRUE)
}

## 1.3 Define Analysis Parameters ----
params <- list(
  strike_distance = 0.45,      # Pike max strike distance (m)
  encounter_threshold_low = 20,  # Minimum encounters for suspected predation
  encounter_threshold_high = 100, # High encounter count threshold
  min_distance_threshold = 0.45,   # Minimum distance threshold (m)
  consecutive_days_threshold = 2,  # Consecutive days for likely predation
  proximity_encounter_threshold = 100  # Total encounters threshold for proximity analysis
)

# ===================================================================
# 2. CUSTOM FUNCTIONS
# ===================================================================

#' Calculate Pairwise Distances Between Prey and Predator
#' 
#' @param prey_tel Telemetry object for prey
#' @param pred_tel Telemetry object for predator
#' @param prey_fits CTMM fits for prey
#' @param pred_fits CTMM fits for predator
#' @param prey_species Species name for prey (e.g., "Roach", "Perch")
#' @param strike_dist Maximum strike distance for encounter definition
#' @param parallel Logical, use parallel processing?
#' @param n_cores Number of cores for parallel processing
#' 
#' @return Data frame with pairwise distances and encounter flags
calculate_pairwise_distances <- function(prey_tel, pred_tel, 
                                         prey_fits, pred_fits,
                                         prey_species = "Prey",
                                         strike_dist = 0.45,
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
      
      # Add encounter flag
      df$encounter <- ifelse(df$est <= strike_dist, 1, 0)
      
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
                                  "combinations"),
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


#' Summarize Daily Encounters
#' 
#' @param distances_df Data frame with distance calculations
#' @return Daily encounter summary
summarize_daily_encounters <- function(distances_df) {
  daily_summary <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      Species = first(Species),
      encounter_count = sum(encounter, na.rm = TRUE),
      daily_avg_dist = mean(est, na.rm = TRUE),
      daily_min_dist = min(est, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(daily_summary)
}


#' Calculate Total Encounters Per Prey-Predator Pair
#' 
#' @param distances_df Data frame with distance calculations
#' @return Summary with total encounters
calculate_total_encounters <- function(distances_df) {
  total_encounters <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Species) %>%
    summarise(
      total_encounters = sum(encounter, na.rm = TRUE),
      min_distance = min(est, na.rm = TRUE),
      mean_distance = mean(est, na.rm = TRUE),
      n_days = n_distinct(Date),
      .groups = "drop"
    ) %>%
    arrange(desc(total_encounters))
  
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


#' Identify Suspected Predation Events
#' 
#' @param distances_df Data frame with pairwise distances
#' @param enc_threshold_low Minimum encounters per day
#' @param enc_threshold_high High encounter threshold
#' @param min_dist_threshold Minimum distance threshold (m)
#' @param consec_days_threshold Consecutive days for likely predation
#' 
#' @return Data frame with suspected predation events
identify_predation_events <- function(distances_df,
                                      enc_threshold_low = 25,
                                      enc_threshold_high = 100,
                                      min_dist_threshold = 0.2,
                                      consec_days_threshold = 2) {
  
  # Filter for potential predation days
  high_encounter_days <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      encounter_count = sum(encounter, na.rm = TRUE),
      daily_min_dist = min(est, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(encounter_count >= enc_threshold_low | daily_min_dist < min_dist_threshold)
  
  # Summarize predation patterns
  predation_summary <- high_encounter_days %>%
    group_by(Prey_ID, Pred_ID) %>%
    summarise(
      num_days_25_encounters = sum(encounter_count >= enc_threshold_low),
      num_days_100_encounters = sum(encounter_count >= enc_threshold_high),
      num_days_min_dist_less_0.2m = sum(daily_min_dist < min_dist_threshold),
      
      first_date_over_25 = if(any(encounter_count >= enc_threshold_low)) 
        min(Date[encounter_count >= enc_threshold_low]) else as.Date(NA),
      first_date_over_100 = if(any(encounter_count >= enc_threshold_high)) 
        min(Date[encounter_count >= enc_threshold_high]) else as.Date(NA),
      
      consecutive_days_25 = calc_max_consecutive(Date[encounter_count >= enc_threshold_low]),
      consecutive_days_100 = calc_max_consecutive(Date[encounter_count >= enc_threshold_high]),
      
      encounter_dates = paste(format(Date, "%d/%m/%Y"), collapse = ", "),
      
      .groups = "drop"
    ) %>%
    mutate(likely_predated = ifelse(consecutive_days_100 >= consec_days_threshold, 1, 0))
  
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
# 4. DISTANCE CALCULATIONS (Already Completed)
# ===================================================================

message("\n=== LOADING DISTANCE CALCULATIONS ===")

## 4.1 Load Pre-calculated Distances ----
roach_pike_file <- file.path(paths$encounters, "muddyfoot_pike_roach_distances_df.rds")
perch_pike_file <- file.path(paths$encounters, "muddyfoot_pike_perch_distances_df.rds")

if (!file.exists(roach_pike_file) || !file.exists(perch_pike_file)) {
  stop("Distance files not found. Please run distance calculations first.")
}

## 4.1 Roach-Pike Distances ----
# Check if already calculated
roach_pike_file <- file.path(paths$encounters, "muddyfoot_pike_roach_distances_df.rds")

if (file.exists(roach_pike_file)) {
  message("Loading existing Roach-Pike distances...")
  roach_pike_distances_df <- readRDS(roach_pike_file)
} else {
  message("Calculating Roach-Pike distances (this may take time)...")
  roach_pike_distances_df <- calculate_pairwise_distances(
    prey_tel = roach_tel,
    pred_tel = pike_tel,
    prey_fits = roach_fits,
    pred_fits = pike_fits,
    prey_species = "Roach",
    strike_dist = params$strike_distance,
    parallel = TRUE,
    n_cores = 16  # Auto-detect
  )
  saveRDS(roach_pike_distances_df, roach_pike_file)
}

## 4.2 Perch-Pike Distances ----
perch_pike_file <- file.path(paths$encounters, "muddyfoot_pike_perch_distances_df.rds")

if (file.exists(perch_pike_file)) {
  message("Loading existing Perch-Pike distances...")
  perch_pike_distances_df <- readRDS(perch_pike_file)
} else {
  message("Calculating Perch-Pike distances (this may take time)...")
  perch_pike_distances_df <- calculate_pairwise_distances(
    prey_tel = perch_tel,
    pred_tel = pike_tel,
    prey_fits = perch_fits,
    pred_fits = pike_fits,
    prey_species = "Perch",
    strike_dist = params$strike_distance,
    parallel = TRUE,
    n_cores = 16
  )
  saveRDS(perch_pike_distances_df, perch_pike_file)
}

roach_pike_distances_df <- readRDS(roach_pike_file)
perch_pike_distances_df <- readRDS(perch_pike_file)

message(sprintf("Loaded distances: %d roach-pike, %d perch-pike observations",
                nrow(roach_pike_distances_df), nrow(perch_pike_distances_df)))


# ===================================================================
# 5. IDENTIFY CANDIDATE PAIRS FOR PROXIMITY ANALYSIS
# ===================================================================

message("\n=== IDENTIFYING CANDIDATES FOR PROXIMITY ANALYSIS ===")

## 5.1 Calculate Total Encounters ----
roach_total <- calculate_total_encounters(roach_pike_distances_df)
perch_total <- calculate_total_encounters(perch_pike_distances_df)

all_encounters <- bind_rows(roach_total, perch_total)

## 5.2 Filter for High-Risk Pairs ----
# Criteria: Not found alive AND high encounter rate (>100 total encounters)
candidate_pairs <- all_encounters %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
  filter(
    found_alive == 0,  # Not recovered alive
    total_encounters > params$proximity_encounter_threshold  # High encounter rate
  ) %>%
  arrange(desc(total_encounters))

message(sprintf("Identified %d candidate prey-predator pairs for proximity analysis",
                nrow(candidate_pairs)))

# Save candidate pairs
saveRDS(candidate_pairs, 
        file.path(paths$proximity, "muddyfoot_proximity_candidate_pairs.rds"))


# ===================================================================
# 6. PROXIMITY ANALYSIS - TEST EXAMPLE
# ===================================================================

message("\n=== TESTING PROXIMITY ANALYSIS (TOP 3 PAIRS) ===")

## 6.1 Select Test Cases ----
# Take top 3 pairs with highest encounters for testing
test_pairs <- candidate_pairs %>%
  slice_head(n = 3)

message("Test pairs:")
print(test_pairs %>% select(Prey_ID, Pred_ID, Species, total_encounters))

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
        file.path(paths$proximity, "test_proximity_results_detailed.rds"))
write.csv(test_proximity_df, 
          file.path(paths$proximity, "test_proximity_results_summary.csv"),
          row.names = FALSE)

message("\nTest proximity analysis complete!")
message(sprintf("Results saved to: %s", paths$proximity))


# ===================================================================
# 7. FULL PROXIMITY ANALYSIS (ALL CANDIDATE PAIRS)
# ===================================================================

message("\n=== RUNNING FULL PROXIMITY ANALYSIS ===")
message("This may take considerable time depending on number of pairs...")
message(sprintf("Analyzing %d candidate pairs", nrow(candidate_pairs)))

## 7.1 Run Proximity on All Candidates ----
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
        file.path(paths$proximity, "full_proximity_results_detailed.rds"))
saveRDS(proximity_with_encounters,
        file.path(paths$proximity, "full_proximity_results_with_encounters.rds"))
write.csv(proximity_with_encounters,
          file.path(paths$proximity, "full_proximity_results.csv"),
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

## 8.2 Identify Suspected Predation (Distance-Based) ----
roach_predation <- identify_predation_events(
  roach_pike_distances_df,
  enc_threshold_low = params$encounter_threshold_low,
  enc_threshold_high = params$encounter_threshold_high,
  min_dist_threshold = params$min_distance_threshold,
  consec_days_threshold = params$consecutive_days_threshold
) %>% rename(Pike_ID = Pred_ID) %>% mutate(Species = "Roach")

perch_predation <- identify_predation_events(
  perch_pike_distances_df,
  enc_threshold_low = params$encounter_threshold_low,
  enc_threshold_high = params$encounter_threshold_high,
  min_dist_threshold = params$min_distance_threshold,
  consec_days_threshold = params$consecutive_days_threshold
) %>% rename(Pike_ID = Pred_ID) %>% mutate(Species = "Perch")

suspected_predation_distance <- bind_rows(roach_predation, perch_predation)

## 8.3 Merge with Proximity Results ----
suspected_predation_enhanced <- suspected_predation_distance %>%
  left_join(
    proximity_with_encounters %>% 
      select(Prey_ID, Pred_ID = Pike_ID, z_statistic, p_value, independent),
    by = c("Prey_ID", "Pike_ID" = "Pred_ID")
  ) %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
  filter(found_alive == 0) %>%
  mutate(
    # Combined predation classification
    proximity_suggests_predation = !independent & !is.na(independent),
    combined_predation_score = case_when(
      likely_predated == 1 & proximity_suggests_predation ~ 3,  # Strong evidence
      likely_predated == 1 & is.na(proximity_suggests_predation) ~ 2,  # Distance only
      proximity_suggests_predation ~ 2,  # Proximity only
      TRUE ~ 1  # Weak evidence
    ),
    classification = case_when(
      combined_predation_score == 3 ~ "High confidence predation",
      combined_predation_score == 2 ~ "Moderate confidence predation",
      TRUE ~ "Low confidence / uncertain"
    )
  )

# Save enhanced results
saveRDS(suspected_predation_enhanced,
        file.path(paths$encounters, "muddyfoot_suspected_predation_enhanced.rds"))

message(sprintf("Enhanced predation classification complete: %d events analyzed",
                nrow(suspected_predation_enhanced)))

# Summary by classification
classification_summary <- suspected_predation_enhanced %>%
  group_by(Species, classification) %>%
  summarise(
    n_events = n(),
    mean_encounters_100 = mean(num_days_100_encounters, na.rm = TRUE),
    mean_p_value = mean(p_value, na.rm = TRUE),
    .groups = "drop"
  )

print("Classification Summary:")
print(classification_summary)


# ===================================================================
# 9. VISUALIZATION OF PROXIMITY RESULTS
# ===================================================================

message("\n=== CREATING VISUALIZATIONS ===")

## 9.1 Proximity P-value Distribution ----
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

ggsave(file.path(paths$proximity, "proximity_pvalue_distribution.png"),
       p1, width = 8, height = 6, dpi = 300)

## 9.2 Z-statistic vs Total Encounters ----
p2 <- ggplot(proximity_with_encounters, aes(x = total_encounters, y = z_statistic)) +
  geom_point(aes(color = independent), size = 3, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("TRUE" = "darkgreen", "FALSE" = "darkred"),
    labels = c("TRUE" = "Independent", "FALSE" = "Non-independent"),
    name = "Movement Pattern"
  ) +
  labs(
    title = "Proximity Analysis: Z-statistic vs Encounter Rate",
    x = "Total Encounters (within strike distance)",
    y = "Z-statistic"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(paths$proximity, "proximity_z_vs_encounters.png"),
       p2, width = 10, height = 7, dpi = 300)

## 9.3 Classification Comparison ----
p3 <- ggplot(suspected_predation_enhanced, 
             aes(x = classification, fill = Species)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  labs(
    title = "Predation Events by Confidence Level",
    x = "Classification",
    y = "Number of Events",
    fill = "Prey Species"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(paths$proximity, "predation_classification_summary.png"),
       p3, width = 10, height = 6, dpi = 300)


# ===================================================================
# 10. EXPORT COMPREHENSIVE RESULTS
# ===================================================================

message("\n=== EXPORTING RESULTS ===")

## 10.1 Create Excel Workbook ----
wb <- createWorkbook()

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

# Save workbook
saveWorkbook(wb, 
             file.path(paths$tables, "muddyfoot_proximity_analysis_complete.xlsx"),
             overwrite = TRUE)

message("Excel workbook saved successfully!")


# ===================================================================
# 11. SUMMARY REPORT
# ===================================================================

message("\n" %+% paste(rep("=", 70), collapse = ""))
message("ANALYSIS COMPLETE - SUMMARY")
message(paste(rep("=", 70), collapse = ""))

message(sprintf("\nTotal prey-predator pairs analyzed: %d", nrow(candidate_pairs)))
message(sprintf("Pairs with non-independent movement (p < 0.05): %d", 
                sum(!proximity_with_encounters$independent, na.rm = TRUE)))
message(sprintf("Percentage showing non-independence: %.1f%%",
                100 * mean(!proximity_with_encounters$independent, na.rm = TRUE)))

message("\nPredation confidence classification:")
print(table(suspected_predation_enhanced$classification))

message("\nOutput locations:")
message(sprintf("  Proximity results: %s", paths$proximity))
message(sprintf("  Encounter data: %s", paths$encounters))
message(sprintf("  Tables: %s", paths$tables))
message(sprintf("  Figures: %s", paths$figures))

message("\n" %+% paste(rep("=", 70), collapse = ""))


# ===================================================================
# 12. SESSION INFO
# ===================================================================

sink(file.path(paths$tables, "session_info_proximity.txt"))
cat("Proximity Analysis Session Info\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")
print(sessionInfo())
sink()

message("\nAll analyses complete! Check output directories for results.")
