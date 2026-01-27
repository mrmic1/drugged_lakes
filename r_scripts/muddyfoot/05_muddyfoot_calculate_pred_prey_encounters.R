# =================================================================-
# Predator-Prey Interaction Analysis - Muddyfoot ####
# =================================================================-
# 
# ANALYSIS OVERVIEW:
# This script identifies potential predation events by analyzing movement 
# tracking data for prey (roach/perch) and predators (pike).
#
# METHODOLOGY:
# 1. TIERED ENCOUNTER DETECTION based on GPS error (0.5m)
#    - High-confidence: ≤1.4m (within 95% combined GPS error)
#    - Probable: 1.4-2.8m (requires temporal confirmation)
#    - Possible: 2.8-4.2m (shared space use)
#
# 2. MOVEMENT CESSATION DETECTION accounting for GPS noise
#    - GPS error produces ~38m/day apparent movement (random walk)
#    - Mortality threshold: ~95m/day (2.5x random walk)
#    - Requires 2+ consecutive days below threshold
#
# 3. OPTIONAL: Proximity analysis (ctmm::proximity) to detect
#    non-independent movement patterns between prey-predator pairs
#
# Author: Marcus Michelangeli
# Date: 2024
# =================================================================-

# =================================================================-
# 1. SETUP & CONFIGURATION ####
# =================================================================-

## 1.1 Load Required Libraries ----
suppressPackageStartupMessages({
  library(dplyr)        # Data manipulation
  library(ctmm)         # Movement modeling and distance calculations
  library(lubridate)    # Date/time handling
  library(data.table)   # Fast data operations
  library(sf)           # Spatial data
  library(flextable)    # Table formatting
  library(parallel)     # Parallel processing
  library(future.apply) # Parallel operations
  library(progressr)    # Progress bars
  library(ggplot2)      # Plotting
  library(openxlsx)     # Excel export
  library(patchwork)    # Combining plots
})

# Configure time zone
Sys.setenv(TZ = 'Europe/Stockholm')

## 1.2 Define Directory Paths ----
paths <- list(
  filtered_data = "./data/tracks_filtered/muddyfoot/",
  lake_polygon  = "./data/lake_params/polygons/",
  ctmm          = "./data/ctmm_fits/",
  telem         = "./data/telem_obj/muddyfoot/",
  encounters    = "./data/encounters/muddyfoot/",
  proximity     = "./data/proximity/muddyfoot/",
  size          = "./data/fish_biometrics/",
  tables        = "./tables/muddyfoot/",
  figures       = "./figures/muddyfoot/",
  distance_plots = "./figures/muddyfoot/distance_plots/"
)

# Create directories if they don't exist
for (path in paths) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

## 1.3 Define Analysis Parameters ----
# Combined GPS error for prey+predator: √(0.5² + 0.5²) = 0.71m
# 95% CI: ~1.4m
params <- list(
  # Tiered distance thresholds (meters)
  strike_distance = 0.45,              # Pike max strike distance
  high_confidence_threshold = 1.4,     # ≤1.4m: Within 95% combined GPS error
  probable_threshold = 2.8,            # ≤2.8m: 2x GPS error
  possible_threshold = 4.2,            # ≤4.2m: 3x GPS error
  
  # Encounter count thresholds
  encounter_threshold_low = 10,        # Min high-confidence encounters/day
  encounter_threshold_high = 50,       # High encounter threshold
  
  # Distance-based thresholds
  min_distance_threshold = 1.4,        # Min distance for high-confidence
  consecutive_days_threshold = 2,      # Consecutive days for predation
  
  # Proximity analysis settings (OPTIONAL - computationally expensive)
  run_proximity_analysis = FALSE,      # Toggle proximity analysis
  proximity_encounter_threshold = 100, # Min encounters for proximity
  
  # Distance plot settings
  create_distance_plots = TRUE,        # Create time series plots
  distance_plot_min_encounters = 50,   # Min encounters to generate plot
  
  # Mortality detection settings
  gps_error_sd = 0.5,                  # Single fish GPS error (meters)
  fixes_per_day = 2880,                # 30-second fix intervals
  mortality_cessation_threshold = 0.15, # 15% of baseline movement
  mortality_consecutive_days = 2       # Consecutive low-movement days
)

# Print analysis strategy
message("\n=== ENCOUNTER DETECTION STRATEGY ===")
message("Tiered approach based on GPS error (σ = 0.5m per fish):")
message(sprintf("  High-confidence: ≤%.1fm (within 95%% combined GPS error)",
                params$high_confidence_threshold))
message(sprintf("  Probable: %.1f-%.1fm (2x GPS error)",
                params$high_confidence_threshold, params$probable_threshold))
message(sprintf("  Possible: %.1f-%.1fm (3x GPS error)",
                params$probable_threshold, params$possible_threshold))
message(sprintf("\nProximity analysis: %s",
                ifelse(params$run_proximity_analysis, "ENABLED", "DISABLED")))
message(sprintf("Distance plots: %s",
                ifelse(params$create_distance_plots, "ENABLED", "DISABLED")))
message("=====================================\n")

# =================================================================-
# 2. CUSTOM FUNCTIONS ####
# =================================================================-

#' Calculate Expected Apparent Movement Due to GPS Error Alone
#' 
#' Estimates how much a stationary object would appear to move due to GPS error.
#' Uses random walk model: error accumulates as sqrt(n) not linearly.
#' 
#' @param gps_error Standard deviation of GPS error (meters)
#' @param fixes_per_day Number of GPS fixes per day
#' 
#' @return Data frame with apparent movement estimates
calculate_gps_noise_threshold <- function(gps_error = 0.5, fixes_per_day = 2880) {
  
  # Error in distance between two consecutive positions
  # Combines errors from both positions: √(σ₁² + σ₂²)
  distance_error_sd <- sqrt(2) * gps_error
  distance_error_95ci <- 2 * distance_error_sd
  
  # Expected apparent movement per fix
  movement_per_fix <- distance_error_sd
  
  # Daily apparent movement (random walk model)
  # Errors partially cancel: use sqrt(n) scaling, not n
  daily_movement_random_walk <- movement_per_fix * sqrt(fixes_per_day)
  
  # Realistic maximum accounting for temporal correlation
  daily_movement_realistic_max <- daily_movement_random_walk * 2
  
  # Theoretical maximum (all errors same direction - impossible)
  daily_movement_theoretical_max <- movement_per_fix * fixes_per_day
  
  results <- data.frame(
    gps_error_sd = gps_error,
    fixes_per_day = fixes_per_day,
    fixes_per_hour = fixes_per_day / 24,
    distance_error_sd = distance_error_sd,
    distance_error_95ci = distance_error_95ci,
    movement_per_fix = movement_per_fix,
    daily_movement_random_walk = daily_movement_random_walk,
    daily_movement_realistic_max = daily_movement_realistic_max,
    daily_movement_theoretical_max = daily_movement_theoretical_max
  )
  
  return(results)
}


#' Calculate Pairwise Distances Between Prey and Predator
#' 
#' Uses ctmm::distances() to calculate separation distances accounting for
#' GPS error and movement uncertainty. Flags encounters using tiered thresholds.
#' 
#' @param prey_tel Telemetry object for prey
#' @param pred_tel Telemetry object for predator
#' @param prey_fits CTMM fits for prey
#' @param pred_fits CTMM fits for predator
#' @param prey_species Species name (e.g., "Roach", "Perch")
#' @param strike_dist Legacy strike distance (kept for compatibility)
#' @param high_conf_threshold High-confidence threshold (m)
#' @param probable_threshold Probable encounter threshold (m)
#' @param possible_threshold Possible encounter threshold (m)
#' @param parallel Use parallel processing?
#' @param n_cores Number of cores (NULL = auto-detect)
#' 
#' @return Data frame with distances and tiered encounter flags
calculate_pairwise_distances <- function(prey_tel, pred_tel,
                                         prey_fits, pred_fits,
                                         prey_species = "Prey",
                                         strike_dist = 0.45,
                                         high_conf_threshold = 1.4,
                                         probable_threshold = 2.8,
                                         possible_threshold = 4.2,
                                         parallel = TRUE,
                                         n_cores = NULL) {
  
  # Ensure matching projections (required by ctmm)
  ctmm::projection(prey_tel) <- ctmm::projection(pred_tel)
  ctmm::projection(prey_fits) <- ctmm::projection(pred_fits)
  
  stopifnot(projection(prey_tel) == projection(pred_tel))
  
  # Prepare all prey-predator combinations
  n_prey <- length(prey_tel)
  n_pred <- length(pred_tel)
  combinations <- expand.grid(prey_idx = 1:n_prey, pred_idx = 1:n_pred)
  
  message(sprintf("Calculating distances for %d %s x %d Pike = %d combinations",
                  n_prey, prey_species, n_pred, nrow(combinations)))
  
  # Define distance calculation function for single pair
  calc_dist <- function(idx) {
    i <- combinations$prey_idx[idx]
    j <- combinations$pred_idx[idx]
    
    tryCatch({
      # Combine telemetry and fits for this pair
      combined_telemetry <- c(prey_tel[i], pred_tel[j])
      combined_ctmm <- c(prey_fits[i], pred_fits[j])
      
      # Calculate pairwise distances using ctmm
      # Returns data frame with: timestamp, est (distance), low, high
      location_difference <- ctmm::distances(combined_telemetry, combined_ctmm)
      
      # Extract IDs
      prey_id <- names(prey_tel)[i]
      pred_id <- names(pred_tel)[j]
      
      # Convert to data frame and add metadata
      df <- as.data.frame(location_difference)
      df$Prey_ID <- prey_id
      df$Pred_ID <- pred_id
      df$Species <- prey_species
      
      # Add tiered encounter flags
      # Legacy flag (for backward compatibility)
      df$encounter <- ifelse(df$est <= strike_dist, 1, 0)
      
      # Tiered flags based on GPS error
      df$encounter_high_conf <- ifelse(df$est <= high_conf_threshold, 1, 0)
      df$encounter_probable <- ifelse(df$est > high_conf_threshold &
                                        df$est <= probable_threshold, 1, 0)
      df$encounter_possible <- ifelse(df$est > probable_threshold &
                                        df$est <= possible_threshold, 1, 0)
      
      # Categorical encounter type
      df$encounter_type <- dplyr::case_when(
        df$est <= high_conf_threshold ~ "high_confidence",
        df$est <= probable_threshold ~ "probable",
        df$est <= possible_threshold ~ "possible",
        TRUE ~ "no_encounter"
      )
      
      # Add date column (for daily aggregation)
      if ("timestamp" %in% colnames(df)) {
        df$timestamp <- as.POSIXct(df$timestamp, tz = "Europe/Stockholm")
        df$Date <- as.Date(df$timestamp, tz = "Europe/Stockholm")
      }
      
      return(df)
    }, error = function(e) {
      warning(sprintf("Error for %s-%s: %s",
                      names(prey_tel)[i], names(pred_tel)[j], e$message))
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
    
    # Export required objects to cluster
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
  
  # Combine results and remove NULLs
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
#' Aggregates encounter data by prey-predator pair and date.
#' Counts encounters in each tier and calculates daily distance metrics.
#' 
#' @param distances_df Data frame from calculate_pairwise_distances()
#' @return Daily encounter summary with tiered counts
summarize_daily_encounters <- function(distances_df) {
  daily_summary <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      Species = first(Species),
      
      # Count encounters by tier
      encounter_count_high_conf = sum(encounter_high_conf, na.rm = TRUE),
      encounter_count_probable = sum(encounter_probable, na.rm = TRUE),
      encounter_count_possible = sum(encounter_possible, na.rm = TRUE),
      encounter_count_legacy = sum(encounter, na.rm = TRUE),
      
      # Daily distance metrics
      daily_avg_dist = mean(est, na.rm = TRUE),
      daily_min_dist = min(est, na.rm = TRUE),
      
      # Check for consecutive probable encounters (temporal confirmation)
      has_consecutive_probable = any(encounter_probable == 1) &
        any(lag(encounter_probable, default = 0) == 1),
      
      .groups = "drop"
    )
  
  return(daily_summary)
}


#' Calculate Total Encounters Per Prey-Predator Pair
#' 
#' Aggregates encounter counts across all observations for each pair.
#' Provides summary statistics useful for identifying high-risk pairs.
#' 
#' @param distances_df Data frame from calculate_pairwise_distances()
#' @return Summary with tiered total encounters
calculate_total_encounters <- function(distances_df) {
  total_encounters <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Species) %>%
    summarise(
      # Total encounters by tier
      total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
      total_encounters_probable = sum(encounter_probable, na.rm = TRUE),
      total_encounters_possible = sum(encounter_possible, na.rm = TRUE),
      total_encounters_legacy = sum(encounter, na.rm = TRUE),
      
      # Distance metrics
      min_distance = min(est, na.rm = TRUE),
      mean_distance = mean(est, na.rm = TRUE),
      
      # Temporal metrics
      n_days = n_distinct(Date),
      n_days_high_conf = n_distinct(Date[encounter_high_conf == 1]),
      n_days_probable = n_distinct(Date[encounter_probable == 1]),
      n_days_possible = n_distinct(Date[encounter_possible == 1]),
      
      .groups = "drop"
    ) %>%
    arrange(desc(total_encounters_high_conf))
  
  return(total_encounters)
}


#' Calculate Maximum Consecutive Days
#' 
#' Finds the longest run of consecutive days in a date vector.
#' Useful for identifying sustained encounter periods.
#' 
#' @param dates Vector of dates
#' @return Maximum number of consecutive days
calc_max_consecutive <- function(dates) {
  # Handle edge cases
  if (length(dates) == 0 || all(is.na(dates))) return(0)
  
  # Clean and sort dates
  dates <- as.Date(dates[!is.na(dates)])
  if (length(dates) == 0) return(0)
  if (length(dates) == 1) return(1)
  
  sorted_dates <- sort(unique(dates))
  if (length(sorted_dates) == 1) return(1)
  
  # Use run-length encoding to find consecutive days
  # diff() == 1 means consecutive days
  consecutive <- rle(diff(sorted_dates) == 1)
  
  # Get lengths of TRUE runs (consecutive days)
  true_runs <- consecutive$lengths[consecutive$values]
  
  # If no consecutive days, return 1
  if (length(true_runs) == 0) return(1)
  
  # Add 1 because we're counting days, not gaps
  max_run <- max(true_runs, na.rm = TRUE)
  if (is.infinite(max_run) || is.na(max_run)) return(1)
  
  return(max_run + 1)
}


#' Perform Proximity Analysis on Prey-Predator Pair
#' 
#' Uses ctmm::proximity() to test if two individuals' movements are
#' statistically independent. Non-independence suggests predation.
#' 
#' Proximity ratio < 1: closer than expected (potential predation)
#' Proximity ratio > 1: farther than expected (avoidance)
#' CI contains 1: movements are independent
#' 
#' @param prey_tel Single prey telemetry object
#' @param pred_tel Single predator telemetry object
#' @param prey_fit CTMM fit for prey
#' @param pred_fit CTMM fit for predator
#' @param prey_id Prey individual ID
#' @param pred_id Predator individual ID
#' 
#' @return List with proximity results and statistics
perform_proximity_analysis <- function(prey_tel, pred_tel,
                                       prey_fit, pred_fit,
                                       prey_id, pred_id) {
  
  tryCatch({
    # Ensure matching projections
    target_proj <- projection(pred_tel)
    projection(prey_tel) <- target_proj
    projection(pred_tel) <- target_proj
    projection(prey_fit) <- target_proj
    projection(pred_fit) <- target_proj
    
    # Create named lists (required format for ctmm::proximity)
    combined_tel <- list(prey_tel, pred_tel)
    names(combined_tel) <- c(prey_id, pred_id)
    
    combined_fits <- list(prey_fit, pred_fit)
    names(combined_fits) <- c(prey_id, pred_id)
    
    # Verify structure
    stopifnot(
      length(combined_tel) == 2,
      length(combined_fits) == 2,
      all(names(combined_tel) == names(combined_fits))
    )
    
    # Run proximity analysis
    message(sprintf("  Running proximity for %s - %s", prey_id, pred_id))
    proximity_result <- ctmm::proximity(combined_tel, combined_fits,
                                        GUESS = ctmm(error = FALSE))
    
    # Extract statistics
    # proximity() returns named vector: 'low', 'est', 'high'
    proximity_stats <- list(
      prey_id = prey_id,
      pred_id = pred_id,
      proximity_ratio_est = unname(proximity_result["est"]),
      proximity_ratio_low = unname(proximity_result["low"]),
      proximity_ratio_high = unname(proximity_result["high"]),
      # Independent if CI contains 1
      independent = proximity_result["low"] < 1 & proximity_result["high"] > 1,
      proximity_object = proximity_result
    )
    
    return(proximity_stats)
    
  }, error = function(e) {
    warning(sprintf("Proximity analysis failed for %s - %s: %s",
                    prey_id, pred_id, e$message))
    return(list(
      prey_id = prey_id,
      pred_id = pred_id,
      proximity_ratio_est = NA,
      proximity_ratio_low = NA,
      proximity_ratio_high = NA,
      independent = NA,
      error = e$message
    ))
  })
}


#' Identify Suspected Predation Events with Tiered Approach
#' 
#' Identifies predation candidates based on:
#' 1. High frequency of close encounters
#' 2. Sustained encounters over consecutive days
#' 3. Very close minimum approach distances
#' 
#' @param distances_df Data frame with pairwise distances
#' @param enc_threshold_low Min encounters per day
#' @param enc_threshold_high High encounter threshold
#' @param min_dist_threshold Distance threshold (m)
#' @param consec_days_threshold Consecutive days required
#' 
#' @return Data frame with suspected predation events
identify_predation_events <- function(distances_df,
                                      enc_threshold_low = 10,
                                      enc_threshold_high = 50,
                                      min_dist_threshold = 1.4,
                                      consec_days_threshold = 2) {
  
  # Aggregate by day
  daily_encounters <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      encounter_count_high_conf = sum(encounter_high_conf, na.rm = TRUE),
      encounter_count_probable = sum(encounter_probable, na.rm = TRUE),
      encounter_count_probable_confirmed = sum(
        encounter_probable == 1 & lag(encounter_probable, default = 0) == 1,
        na.rm = TRUE
      ),
      daily_min_dist = min(est, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Filter for high-risk days
  high_risk_days <- daily_encounters %>%
    filter(
      encounter_count_high_conf >= enc_threshold_low |
        encounter_count_probable_confirmed >= enc_threshold_low |
        daily_min_dist < min_dist_threshold
    )
  
  # Summarize predation patterns
  predation_summary <- high_risk_days %>%
    group_by(Prey_ID, Pred_ID) %>%
    summarise(
      # Count days with different encounter intensities
      num_days_high_conf_10 = sum(encounter_count_high_conf >= 10),
      num_days_high_conf_50 = sum(encounter_count_high_conf >= 50),
      num_days_probable_confirmed_10 = sum(encounter_count_probable_confirmed >= 10),
      
      # Distance-based metrics
      num_days_min_dist_less_1.5m = sum(daily_min_dist < 1.5),
      num_days_min_dist_less_0.5m = sum(daily_min_dist < 0.5),
      
      # First occurrence dates
      first_date_high_conf = if(any(encounter_count_high_conf >= enc_threshold_low)) {
        min(Date[encounter_count_high_conf >= enc_threshold_low])
      } else {
        as.Date(NA)
      },
      first_date_probable = if(any(encounter_count_probable_confirmed >= enc_threshold_low)) {
        min(Date[encounter_count_probable_confirmed >= enc_threshold_low])
      } else {
        as.Date(NA)
      },
      
      # Consecutive days patterns
      consecutive_days_high_conf = calc_max_consecutive(
        Date[encounter_count_high_conf >= enc_threshold_low]
      ),
      consecutive_days_probable = calc_max_consecutive(
        Date[encounter_count_probable_confirmed >= enc_threshold_low]
      ),
      
      # Store encounter dates
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
      # Binary predation flag based on strong evidence
      likely_predated = case_when(
        consecutive_days_high_conf >= consec_days_threshold ~ 1,
        num_days_high_conf_50 >= 2 ~ 1,
        num_days_min_dist_less_0.5m >= 2 ~ 1,
        TRUE ~ 0
      )
    )
  
  return(predation_summary)
}


#' Plot Distance Time Series for Prey-Predator Pair
#' 
#' Creates time series plot of separation distance with annotations for:
#' - Encounter thresholds
#' - First date with 50+ high-confidence encounters
#' - Date with maximum encounters (if total > 100)
#' 
#' @param distances_df Data frame with distance calculations
#' @param prey_id Prey individual ID
#' @param pred_id Predator individual ID
#' @param species Species name
#' @param output_dir Directory to save plot
#' @param threshold_lines Add threshold lines?
#' 
#' @return ggplot object
plot_distance_timeseries <- function(distances_df,
                                     prey_id,
                                     pred_id,
                                     species,
                                     output_dir,
                                     threshold_lines = TRUE) {
  
  # Filter for this pair
  pair_data <- distances_df %>%
    filter(Prey_ID == prey_id, Pred_ID == pred_id) %>%
    arrange(timestamp)
  
  if (nrow(pair_data) == 0) {
    warning(sprintf("No data found for %s - %s", prey_id, pred_id))
    return(NULL)
  }
  
  # Calculate summary statistics
  n_encounters_high <- sum(pair_data$encounter_high_conf)
  n_encounters_prob <- sum(pair_data$encounter_probable)
  mean_dist <- mean(pair_data$est, na.rm = TRUE)
  min_dist <- min(pair_data$est, na.rm = TRUE)
  
  # Calculate daily high-confidence encounter counts
  daily_high_conf <- pair_data %>%
    group_by(Date) %>%
    summarise(
      daily_high_conf_count = sum(encounter_high_conf, na.rm = TRUE),
      timestamp = first(timestamp),
      .groups = "drop"
    )
  
  # Find first date with 50+ encounters
  first_50_date <- daily_high_conf %>%
    filter(daily_high_conf_count >= 50) %>%
    slice_min(Date, n = 1)
  
  # Find date with most encounters (only if total > 100)
  max_encounter_date <- NULL
  if (n_encounters_high > 100) {
    max_encounter_date <- daily_high_conf %>%
      slice_max(daily_high_conf_count, n = 1)
  }
  
  # Set up x-axis breaks (every 2 days)
  date_range <- range(pair_data$Date, na.rm = TRUE)
  all_dates <- seq(from = date_range[1], to = date_range[2], by = "day")
  date_breaks <- all_dates[seq(1, length(all_dates), by = 2)]
  
  # Create base plot
  p <- ggplot(pair_data, aes(x = timestamp, y = est)) +
    geom_line(color = "#5e548e", alpha = 0.7, linewidth = 0.5) +
    scale_x_datetime(
      breaks = as.POSIXct(paste(date_breaks, "12:00:00"), tz = "Europe/Stockholm"),
      date_labels = "%d %b",
      expand = c(0.02, 0)
    ) +
    labs(
      title = sprintf("%s: %s vs Pike %s", species, prey_id, pred_id),
      subtitle = sprintf(
        "High-conf encounters: %d | Probable: %d | Mean dist: %.1fm | Min dist: %.2fm",
        n_encounters_high, n_encounters_prob, mean_dist, min_dist
      ),
      x = "Date",
      y = "Separation distance (m)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "gray30"),
      axis.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
  # Add threshold lines
  if (threshold_lines) {
    p <- p +
      geom_hline(yintercept = 1.4, linetype = "dashed",
                 color = "#d73027", alpha = 0.7, linewidth = 0.8) +
      geom_hline(yintercept = 2.8, linetype = "dashed",
                 color = "#fee090", alpha = 0.7, linewidth = 0.8) +
      annotate("text", x = min(pair_data$timestamp), y = 1.4,
               label = "High-conf (1.4m)", hjust = 0, vjust = -0.5,
               size = 3, color = "#d73027") +
      annotate("text", x = min(pair_data$timestamp), y = 2.8,
               label = "Probable (2.8m)", hjust = 0, vjust = -0.5,
               size = 3, color = "#fee090")
  }
  
  # Add first 50+ encounter date
  if (nrow(first_50_date) > 0) {
    p <- p +
      geom_vline(xintercept = as.numeric(first_50_date$timestamp[1]),
                 linetype = "dotted", color = "#e63946", linewidth = 1) +
      annotate("text",
               x = first_50_date$timestamp[1],
               y = max(pair_data$est, na.rm = TRUE) * 0.95,
               label = sprintf("First 50+ encounters\n%s (%d)",
                               format(first_50_date$Date[1], "%d %b"),
                               first_50_date$daily_high_conf_count[1]),
               hjust = -0.05,
               vjust = 1,
               size = 3,
               color = "#e63946",
               fontface = "bold")
  }
  
  # Add max encounter date (if different and total > 100)
  if (!is.null(max_encounter_date) && nrow(max_encounter_date) > 0) {
    is_different <- nrow(first_50_date) == 0 ||
      max_encounter_date$Date[1] != first_50_date$Date[1]
    
    if (is_different) {
      p <- p +
        geom_vline(xintercept = as.numeric(max_encounter_date$timestamp[1]),
                   linetype = "dotdash", color = "#f77f00", linewidth = 1) +
        annotate("text",
                 x = max_encounter_date$timestamp[1],
                 y = max(pair_data$est, na.rm = TRUE) * 0.85,
                 label = sprintf("Max encounters\n%s (%d)",
                                 format(max_encounter_date$Date[1], "%d %b"),
                                 max_encounter_date$daily_high_conf_count[1]),
                 hjust = -0.05,
                 vjust = 1,
                 size = 3,
                 color = "#f77f00",
                 fontface = "bold")
    }
  }
  
  # Save plot
  filename <- sprintf("%s_%s_%s_distance_timeseries.png",
                      species, prey_id, pred_id)
  ggsave(
    filename = file.path(output_dir, filename),
    plot = p,
    width = 12,
    height = 6,
    dpi = 300
  )
  
  message(sprintf("  Saved distance plot: %s", filename))
  
  return(p)
}


#' Plot Movement Metrics for Missing Individuals with Mortality Detection
#' 
#' Creates plots showing daily distance moved and speed for individuals
#' not recovered at study end. Detects potential mortality dates by
#' identifying sustained periods of reduced movement, accounting for
#' GPS-induced apparent movement.
#' 
#' MORTALITY DETECTION LOGIC:
#' 1. Calculate GPS noise threshold (2.5x random walk estimate)
#' 2. Set mortality threshold as MAX of:
#'    - 15% of baseline movement (relative)
#'    - GPS noise threshold (absolute)
#' 3. Flag days with movement below threshold
#' 4. Identify first occurrence of 2+ consecutive low-movement days
#' 
#' @param telemetry_obj Single telemetry object
#' @param individual_id Individual ID
#' @param species Species name
#' @param output_dir Directory to save plots
#' @param gps_error GPS position error SD (meters)
#' @param cessation_threshold Proportion of baseline for relative threshold
#' @param consecutive_days Consecutive low-movement days required
#' 
#' @return List with plots, metrics, and mortality info
plot_movement_metrics <- function(telemetry_obj,
                                  individual_id,
                                  species,
                                  output_dir,
                                  gps_error = 0.5,
                                  cessation_threshold = 0.15,
                                  consecutive_days = 2) {
  
  tryCatch({
    # Validate input
    if (!inherits(telemetry_obj, "telemetry")) {
      stop("Input is not a telemetry object")
    }
    
    # Extract coordinates and timestamps
    tel_df <- data.frame(
      timestamp = telemetry_obj$timestamp,
      x = telemetry_obj$x,
      y = telemetry_obj$y,
      stringsAsFactors = FALSE
    )
    
    if (nrow(tel_df) < 2) {
      warning(sprintf("Insufficient data for %s", individual_id))
      return(NULL)
    }
    
    # Format timestamps
    tel_df$timestamp <- as.POSIXct(tel_df$timestamp,
                                   origin = "1970-01-01",
                                   tz = "Europe/Stockholm")
    tel_df$Date <- as.Date(tel_df$timestamp, tz = "Europe/Stockholm")
    
    # Calculate step distances and speeds
    tel_df <- tel_df %>%
      arrange(timestamp) %>%
      mutate(
        # Euclidean distance between consecutive points
        step_distance = sqrt((x - lag(x))^2 + (y - lag(y))^2),
        # Time interval in hours
        time_interval = as.numeric(difftime(timestamp, lag(timestamp), units = "hours")),
        # Speed (m/hour)
        speed = step_distance / time_interval
      )
    
    # Aggregate to daily metrics
    daily_metrics <- tel_df %>%
      group_by(Date) %>%
      summarise(
        total_distance = sum(step_distance, na.rm = TRUE),
        mean_speed = mean(speed, na.rm = TRUE),
        median_speed = median(speed, na.rm = TRUE),
        max_speed = max(speed, na.rm = TRUE),
        n_observations = n(),
        .groups = "drop"
      ) %>%
      filter(!is.na(Date)) %>%
      arrange(Date)
    
    if (nrow(daily_metrics) == 0) {
      warning(sprintf("No daily metrics calculated for %s", individual_id))
      return(NULL)
    }
    
    # === GPS NOISE CALCULATION ===
    
    mean_fixes_per_day <- mean(daily_metrics$n_observations, na.rm = TRUE)
    
    gps_noise_estimate <- calculate_gps_noise_threshold(
      gps_error = gps_error,
      fixes_per_day = mean_fixes_per_day
    )
    
    # Set absolute threshold as 2.5x random walk
    # Accounts for GPS noise plus some temporal correlation
    gps_noise_threshold <- gps_noise_estimate$daily_movement_random_walk * 2.5
    
    # === MORTALITY DETECTION ===
    
    # Calculate baseline (first 50% or first 7 days)
    n_baseline_days <- max(7, floor(nrow(daily_metrics) * 0.5))
    baseline_days <- min(n_baseline_days, nrow(daily_metrics))
    
    baseline_distance <- daily_metrics %>%
      slice_head(n = baseline_days) %>%
      pull(total_distance) %>%
      mean(na.rm = TRUE)
    
    baseline_speed <- daily_metrics %>%
      slice_head(n = baseline_days) %>%
      pull(mean_speed) %>%
      mean(na.rm = TRUE)
    
    # Mortality threshold = MAX(relative, absolute)
    cessation_distance_threshold <- max(
      baseline_distance * cessation_threshold,
      gps_noise_threshold
    )
    
    cessation_speed_threshold <- max(
      baseline_speed * cessation_threshold,
      gps_noise_threshold / 24
    )
    
    # Flag low-movement days
    daily_metrics <- daily_metrics %>%
      mutate(
        low_movement = (total_distance < cessation_distance_threshold) |
          (mean_speed < cessation_speed_threshold),
        consecutive_low = NA
      )
    
    # Find consecutive low-movement periods
    if (any(daily_metrics$low_movement, na.rm = TRUE)) {
      rle_result <- rle(daily_metrics$low_movement)
      
      run_ends <- cumsum(rle_result$lengths)
      run_starts <- c(1, run_ends[-length(run_ends)] + 1)
      
      long_low_periods <- which(rle_result$values == TRUE &
                                  rle_result$lengths >= consecutive_days)
      
      if (length(long_low_periods) > 0) {
        for (period_idx in long_low_periods) {
          start_idx <- run_starts[period_idx]
          end_idx <- run_ends[period_idx]
          daily_metrics$consecutive_low[start_idx:end_idx] <- TRUE
        }
      }
    }
    
    # Identify mortality date
    mortality_date <- NULL
    mortality_info <- NULL
    
    if (any(daily_metrics$consecutive_low == TRUE, na.rm = TRUE)) {
      mortality_idx <- which(daily_metrics$consecutive_low == TRUE)[1]
      mortality_date <- daily_metrics$Date[mortality_idx]
      
      post_mortality <- daily_metrics %>%
        filter(Date >= mortality_date)
      
      mortality_info <- list(
        date = mortality_date,
        days_tracked_after = nrow(post_mortality),
        mean_distance_after = mean(post_mortality$total_distance, na.rm = TRUE),
        mean_speed_after = mean(post_mortality$mean_speed, na.rm = TRUE),
        baseline_distance = baseline_distance,
        baseline_speed = baseline_speed,
        gps_noise_threshold = gps_noise_threshold,
        reduction_distance = 100 * (1 - mean(post_mortality$total_distance, na.rm = TRUE) / baseline_distance),
        reduction_speed = 100 * (1 - mean(post_mortality$mean_speed, na.rm = TRUE) / baseline_speed),
        # KEY METRIC: Movement relative to GPS noise
        # <1.0 = very likely dead, 1-2 = possibly dead, >2 = uncertain
        movement_vs_gps_noise = mean(post_mortality$total_distance, na.rm = TRUE) / gps_noise_threshold
      )
    }
    
    # Summary statistics
    total_days <- nrow(daily_metrics)
    mean_daily_dist <- mean(daily_metrics$total_distance, na.rm = TRUE)
    mean_daily_speed <- mean(daily_metrics$mean_speed, na.rm = TRUE)
    last_date <- max(daily_metrics$Date, na.rm = TRUE)
    
    # X-axis breaks (every 2 days)
    date_range <- range(daily_metrics$Date, na.rm = TRUE)
    all_dates <- seq(from = date_range[1], to = date_range[2], by = "day")
    date_breaks <- all_dates[seq(1, length(all_dates), by = 2)]
    
    # === PLOT 1: Daily Total Distance ===
    p_distance <- ggplot(daily_metrics, aes(x = Date, y = total_distance)) +
      geom_line(color = "#264653", linewidth = 0.8, alpha = 0.8) +
      geom_point(color = "#264653", size = 2, alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE, color = "#e76f51",
                  fill = "#e76f51", alpha = 0.2, linewidth = 0.6)
    
    # Add GPS noise reference line
    p_distance <- p_distance +
      geom_hline(yintercept = gps_noise_threshold,
                 linetype = "dotted", color = "gray40", alpha = 0.7, linewidth = 0.8) +
      annotate("text",
               x = min(daily_metrics$Date),
               y = gps_noise_threshold,
               label = sprintf("GPS noise (~%.0fm)", gps_noise_threshold),
               hjust = 0, vjust = -0.5, size = 2.5, color = "gray40")
    
    # Add mortality indicator
    if (!is.null(mortality_date)) {
      p_distance <- p_distance +
        geom_vline(xintercept = as.numeric(mortality_date),
                   linetype = "dashed", color = "#d62828", linewidth = 1.2) +
        geom_hline(yintercept = cessation_distance_threshold,
                   linetype = "dotted", color = "#d62828", alpha = 0.5, linewidth = 0.8) +
        annotate("text",
                 x = mortality_date,
                 y = max(daily_metrics$total_distance, na.rm = TRUE) * 0.95,
                 label = sprintf("Potential mortality\n%s\n(%.0f%% ↓, %.1fx GPS noise)",
                                 format(mortality_date, "%d %b"),
                                 mortality_info$reduction_distance,
                                 mortality_info$movement_vs_gps_noise),
                 hjust = -0.05,
                 vjust = 1,
                 size = 3.5,
                 color = "#d62828",
                 fontface = "bold")
    }
    
    p_distance <- p_distance +
      scale_x_date(
        breaks = date_breaks,
        date_labels = "%d %b",
        expand = c(0.02, 0)
      ) +
      labs(
        title = sprintf("%s %s: Daily Total Distance Moved", species, individual_id),
        subtitle = sprintf(
          "Mean: %.1fm/day | Baseline: %.1fm/day | GPS noise: ~%.0fm | Days: %d | Last: %s",
          mean_daily_dist, baseline_distance, gps_noise_threshold, total_days, format(last_date, "%d %b")
        ),
        x = "Date",
        y = "Total distance moved (m)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, color = "gray30"),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )
    
    # === PLOT 2: Daily Average Speed ===
    speed_gps_threshold <- gps_noise_threshold / 24
    
    p_speed <- ggplot(daily_metrics, aes(x = Date, y = mean_speed)) +
      geom_line(color = "#2a9d8f", linewidth = 0.8, alpha = 0.8) +
      geom_point(color = "#2a9d8f", size = 2, alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE, color = "#f4a261",
                  fill = "#f4a261", alpha = 0.2, linewidth = 0.6)
    
    # Add GPS noise reference line
    p_speed <- p_speed +
      geom_hline(yintercept = speed_gps_threshold,
                 linetype = "dotted", color = "gray40", alpha = 0.7, linewidth = 0.8) +
      annotate("text",
               x = min(daily_metrics$Date),
               y = speed_gps_threshold,
               label = sprintf("GPS noise (~%.1fm/hr)", speed_gps_threshold),
               hjust = 0, vjust = -0.5, size = 2.5, color = "gray40")
    
    # Add mortality indicator
    if (!is.null(mortality_date)) {
      p_speed <- p_speed +
        geom_vline(xintercept = as.numeric(mortality_date),
                   linetype = "dashed", color = "#d62828", linewidth = 1.2) +
        geom_hline(yintercept = cessation_speed_threshold,
                   linetype = "dotted", color = "#d62828", alpha = 0.5, linewidth = 0.8) +
        annotate("text",
                 x = mortality_date,
                 y = max(daily_metrics$mean_speed, na.rm = TRUE) * 0.95,
                 label = sprintf("Potential mortality\n%s\n(%.0f%% ↓ speed)",
                                 format(mortality_date, "%d %b"),
                                 mortality_info$reduction_speed),
                 hjust = -0.05,
                 vjust = 1,
                 size = 3.5,
                 color = "#d62828",
                 fontface = "bold")
    }
    
    p_speed <- p_speed +
      scale_x_date(
        breaks = date_breaks,
        date_labels = "%d %b",
        expand = c(0.02, 0)
      ) +
      labs(
        title = sprintf("%s %s: Daily Average Speed", species, individual_id),
        subtitle = sprintf(
          "Mean: %.1fm/hr | Baseline: %.1fm/hr | GPS noise: ~%.1fm/hr | Max: %.1fm/hr",
          mean_daily_speed, baseline_speed, speed_gps_threshold,
          max(daily_metrics$mean_speed, na.rm = TRUE)
        ),
        x = "Date",
        y = "Average speed (m/hour)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, color = "gray30"),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )
    
    # === COMBINED PLOT ===
    combined_title <- sprintf("%s %s: Movement Summary (Not Recovered)",
                              species, individual_id)
    if (!is.null(mortality_date)) {
      combined_title <- sprintf(
        "%s %s: Movement Summary - Potential Mortality: %s (%.1fx GPS noise)",
        species, individual_id, format(mortality_date, "%d %b %Y"),
        mortality_info$movement_vs_gps_noise
      )
    }
    
    p_combined <- p_distance / p_speed +
      plot_annotation(
        title = combined_title,
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    # Save plots
    ggsave(
      filename = file.path(output_dir, sprintf("%s_%s_daily_distance.png",
                                               species, individual_id)),
      plot = p_distance,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(output_dir, sprintf("%s_%s_daily_speed.png",
                                               species, individual_id)),
      plot = p_speed,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(output_dir, sprintf("%s_%s_movement_summary.png",
                                               species, individual_id)),
      plot = p_combined,
      width = 10,
      height = 10,
      dpi = 300
    )
    
    # Print summary
    if (!is.null(mortality_date)) {
      message(sprintf("  %s: Potential mortality on %s (%.0f%% ↓ dist, %.0f%% ↓ speed, %.1fx GPS noise)",
                      individual_id,
                      format(mortality_date, "%d %b %Y"),
                      mortality_info$reduction_distance,
                      mortality_info$reduction_speed,
                      mortality_info$movement_vs_gps_noise))
    } else {
      message(sprintf("  %s: No mortality signal (movement > %.0fm/day threshold)",
                      individual_id, cessation_distance_threshold))
    }
    
    return(list(
      distance = p_distance,
      speed = p_speed,
      combined = p_combined,
      metrics = daily_metrics,
      mortality_info = mortality_info,
      baseline_distance = baseline_distance,
      baseline_speed = baseline_speed,
      gps_noise_threshold = gps_noise_threshold,
      cessation_threshold = cessation_distance_threshold
    ))
    
  }, error = function(e) {
    warning(sprintf("Failed to create movement plots for %s: %s",
                    individual_id, e$message))
    return(NULL)
  })
}


# =================================================================-
# 3. DATA LOADING ####
# =================================================================-

## 3.1 Load Telemetry Objects ----
pike_tel <- readRDS(file.path(paths$telem, 'pike_muddyfoot_tel_thinned.rds'))
perch_tel <- readRDS(file.path(paths$telem, 'perch_muddyfoot_tel_thinned.rds'))
roach_tel <- readRDS(file.path(paths$telem, 'roach_muddyfoot_tel_thinned.rds'))

message(sprintf("Loaded telemetry: %d pike, %d perch, %d roach",
                length(pike_tel), length(perch_tel), length(roach_tel)))

# Remove problematic individual
roach_tel <- roach_tel[names(roach_tel) != "F59707"]

## 3.2 Load CTMM Fits ----
pike_fits <- readRDS(file.path(paths$ctmm, "muddyfoot_pike_fits/muddyfoot_pike_best_models.rds"))
perch_fits <- readRDS(file.path(paths$ctmm, "muddyfoot_perch_fits/muddyfoot_perch_best_models.rds"))
roach_fits <- readRDS(file.path(paths$ctmm, "muddyfoot_roach_fits/muddyfoot_roach_best_models.rds"))

message(sprintf("Loaded ctmms: %d pike, %d perch, %d roach",
                length(pike_fits), length(perch_fits), length(roach_fits)))

## 3.3 Load Metadata ----
muddyfoot_filt_data <- readRDS(file.path(paths$filtered_data, "04_muddyfoot_sub.rds"))
post_biometrics <- fread(file.path(paths$size, "biometric_post_exp_data.csv")) %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number)))

# Create biometric filtering columns
post_biometric_cols <- post_biometrics %>%
  filter(Lake == 'Muddyfoot', Species %in% c('Roach', 'Perch')) %>%
  select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found, found_predated = Known_predated)


# =================================================================-
# 4. DISTANCE CALCULATIONS ####
# =================================================================-
# 
# Uses ctmm::distances() to calculate separation distances between
# all prey-predator pairs. Accounts for GPS error and movement
# uncertainty in distance estimates.
# =================================================================-

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

## 4.3 Print Encounter Summary ----
message("\n=== ENCOUNTER TYPE DISTRIBUTION ===")
combined_distances <- bind_rows(roach_pike_distances_df, perch_pike_distances_df)
encounter_summary <- combined_distances %>%
  group_by(Species, encounter_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Species) %>%
  mutate(pct = 100 * n / sum(n))

print(encounter_summary)


# =================================================================-
# 5. DAILY & TOTAL ENCOUNTER SUMMARIES ####
# =================================================================-
# 
# Aggregates encounter data to daily and total summaries.
# These summaries are used to identify high-risk prey-predator pairs.
# =================================================================-

message("\n=== CALCULATING ENCOUNTER SUMMARIES ===")

## 5.1 Daily Encounter Summaries ----
roach_daily <- summarize_daily_encounters(roach_pike_distances_df) %>%
  rename(Pike_ID = Pred_ID)

perch_daily <- summarize_daily_encounters(perch_pike_distances_df) %>%
  rename(Pike_ID = Pred_ID)

prey_daily_encounters <- bind_rows(roach_daily, perch_daily)

message(sprintf("Daily encounters summarized: %d days", nrow(prey_daily_encounters)))

## 5.2 Calculate Total Encounters ----
roach_total <- calculate_total_encounters(roach_pike_distances_df)
perch_total <- calculate_total_encounters(perch_pike_distances_df)

all_encounters <- bind_rows(roach_total, perch_total)

message(sprintf("Total encounter pairs: %d", nrow(all_encounters)))


# =================================================================-
# 6. IDENTIFY SUSPECTED PREDATION EVENTS ####
# =================================================================-
# 
# Identifies predation candidates based on:
# 1. High frequency of close encounters (50+ high-confidence/day)
# 2. Sustained encounters over consecutive days (2+)
# 3. Very close minimum approaches (<0.5m)
# =================================================================-

message("\n=== IDENTIFYING SUSPECTED PREDATION EVENTS ===")

## 6.1 Roach Predation Events ----
roach_predation <- identify_predation_events(
  roach_pike_distances_df,
  enc_threshold_low = params$encounter_threshold_low,
  enc_threshold_high = params$encounter_threshold_high,
  min_dist_threshold = params$min_distance_threshold,
  consec_days_threshold = params$consecutive_days_threshold
) %>%
  rename(Pike_ID = Pred_ID) %>%
  mutate(Species = "Roach")

## 6.2 Perch Predation Events ----
perch_predation <- identify_predation_events(
  perch_pike_distances_df,
  enc_threshold_low = params$encounter_threshold_low,
  enc_threshold_high = params$encounter_threshold_high,
  min_dist_threshold = params$min_distance_threshold,
  consec_days_threshold = params$consecutive_days_threshold
) %>%
  rename(Pike_ID = Pred_ID) %>%
  mutate(Species = "Perch")

## 6.3 Combine and Filter ----
# Filter for individuals not found alive AND with substantial encounters
suspected_predation_distance <- bind_rows(roach_predation, perch_predation) %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
  filter(
    found_alive == 0,           # Not recovered alive
    num_days_high_conf_50 > 1   # At least 2 days with 50+ encounters
  )

message(sprintf("Suspected predation events identified: %d",
                nrow(suspected_predation_distance)))

# Print summary
predation_summary_table <- suspected_predation_distance %>%
  group_by(Species, likely_predated) %>%
  summarise(
    n_cases = n(),
    mean_high_conf_days = mean(num_days_high_conf_10, na.rm = TRUE),
    mean_consecutive_days = mean(consecutive_days_high_conf, na.rm = TRUE),
    .groups = "drop"
  )

message("\nPredation Summary:")
print(predation_summary_table)


# =================================================================-
# 7. DISTANCE TIME SERIES PLOTS (OPTIONAL) ####
# =================================================================-
# 
# Creates time series plots showing prey-predator separation distances
# over time. Annotates key dates (first 50+ encounters, max encounters).
# =================================================================-

if (params$create_distance_plots) {
  message("\n=== CREATING DISTANCE TIME SERIES PLOTS ===")
  
  # Get pairs with sufficient encounters for plotting
  high_encounter_pairs <- suspected_predation_distance %>%
    filter(
      num_days_high_conf_10 >= 2 |
        num_days_high_conf_50 >= 1
    ) %>%
    select(Prey_ID, Pike_ID, Species) %>%
    rename(Pred_ID = Pike_ID)
  
  message(sprintf("Creating distance plots for %d pairs", nrow(high_encounter_pairs)))
  
  # Create plots
  plot_list <- list()
  
  for (i in 1:nrow(high_encounter_pairs)) {
    prey_id <- high_encounter_pairs$Prey_ID[i]
    pred_id <- high_encounter_pairs$Pred_ID[i]
    species <- high_encounter_pairs$Species[i]
    
    # Get appropriate distance data
    if (species == "Roach") {
      dist_data <- roach_pike_distances_df
    } else {
      dist_data <- perch_pike_distances_df
    }
    
    # Create plot
    p <- plot_distance_timeseries(
      distances_df = dist_data,
      prey_id = prey_id,
      pred_id = pred_id,
      species = species,
      output_dir = paths$distance_plots,
      threshold_lines = TRUE
    )
    
    plot_list[[i]] <- p
  }
  
  message(sprintf("Distance plots complete: %d plots created", length(plot_list)))
}


# =================================================================-
# 7B. MOVEMENT METRICS PLOTS FOR MISSING INDIVIDUALS ####
# =================================================================-
# 
# Creates plots showing daily movement patterns for individuals not
# recovered at study end. Detects potential mortality dates by
# identifying sustained periods of reduced movement that exceed
# what would be expected from GPS error alone.
#
# GPS ERROR CONTEXT (30-second fixes, σ = 0.5m):
# - Random walk apparent movement: ~38 m/day
# - Mortality threshold: ~95 m/day (2.5x random walk)
# - Movement < GPS noise = very likely dead
# - Movement 1-2x GPS noise = possibly dead
# - Movement > 2x GPS noise = uncertain
# =================================================================-

message("\n=== CREATING MOVEMENT METRICS PLOTS FOR MISSING INDIVIDUALS ===")

# Create directory
movement_plots_dir <- file.path(paths$figures, "movement_metrics")
if (!dir.exists(movement_plots_dir)) {
  dir.create(movement_plots_dir, recursive = TRUE)
}

# Get individuals not found alive
missing_individuals <- post_biometric_cols %>%
  filter(found_alive == 0)

message(sprintf("Creating movement plots for %d missing individuals",
                nrow(missing_individuals)))

# Calculate GPS noise reference
gps_reference <- calculate_gps_noise_threshold(
  gps_error = params$gps_error_sd,
  fixes_per_day = params$fixes_per_day
)

message(sprintf("\nGPS Error Context (σ = %.2fm, %d fixes/day):",
                params$gps_error_sd, params$fixes_per_day))
message(sprintf("  Random walk apparent movement: ~%.0fm/day",
                gps_reference$daily_movement_random_walk))
message(sprintf("  Mortality detection threshold: ~%.0fm/day (2.5x random walk)",
                gps_reference$daily_movement_random_walk * 2.5))

# Initialize storage
movement_plot_results <- list()

## Process Roach ----
missing_roach <- missing_individuals %>%
  filter(individual_ID %in% names(roach_tel))

if (nrow(missing_roach) > 0) {
  message(sprintf("\nProcessing %d missing Roach individuals...", nrow(missing_roach)))
  
  for (i in 1:nrow(missing_roach)) {
    ind_id <- missing_roach$individual_ID[i]
    
    if (ind_id %in% names(roach_tel)) {
      result <- plot_movement_metrics(
        telemetry_obj = roach_tel[[ind_id]],
        individual_id = ind_id,
        species = "Roach",
        output_dir = movement_plots_dir,
        gps_error = params$gps_error_sd,
        cessation_threshold = params$mortality_cessation_threshold,
        consecutive_days = params$mortality_consecutive_days
      )
      
      movement_plot_results[[ind_id]] <- result
    }
  }
}

## Process Perch ----
missing_perch <- missing_individuals %>%
  filter(individual_ID %in% names(perch_tel))

if (nrow(missing_perch) > 0) {
  message(sprintf("\nProcessing %d missing Perch individuals...", nrow(missing_perch)))
  
  for (i in 1:nrow(missing_perch)) {
    ind_id <- missing_perch$individual_ID[i]
    
    if (ind_id %in% names(perch_tel)) {
      result <- plot_movement_metrics(
        telemetry_obj = perch_tel[[ind_id]],
        individual_id = ind_id,
        species = "Perch",
        output_dir = movement_plots_dir,
        gps_error = params$gps_error_sd,
        cessation_threshold = params$mortality_cessation_threshold,
        consecutive_days = params$mortality_consecutive_days
      )
      
      movement_plot_results[[ind_id]] <- result
    }
  }
}

## Compile Results ----

# Compile daily metrics
all_daily_metrics <- bind_rows(
  lapply(names(movement_plot_results), function(ind_id) {
    result <- movement_plot_results[[ind_id]]
    if (!is.null(result) && !is.null(result$metrics)) {
      result$metrics %>%
        mutate(
          individual_ID = ind_id,
          species = ifelse(ind_id %in% names(roach_tel), "Roach", "Perch")
        )
    }
  })
)

# Compile mortality information
mortality_summary <- bind_rows(
  lapply(names(movement_plot_results), function(ind_id) {
    result <- movement_plot_results[[ind_id]]
    species <- ifelse(ind_id %in% names(roach_tel), "Roach", "Perch")
    
    if (!is.null(result)) {
      if (!is.null(result$mortality_info)) {
        data.frame(
          individual_ID = ind_id,
          species = species,
          potential_mortality_date = result$mortality_info$date,
          days_tracked_after = result$mortality_info$days_tracked_after,
          baseline_distance = result$baseline_distance,
          baseline_speed = result$baseline_speed,
          mean_distance_after = result$mortality_info$mean_distance_after,
          mean_speed_after = result$mortality_info$mean_speed_after,
          reduction_distance_pct = result$mortality_info$reduction_distance,
          reduction_speed_pct = result$mortality_info$reduction_speed,
          gps_noise_threshold = result$gps_noise_threshold,
          movement_vs_gps_noise = result$mortality_info$movement_vs_gps_noise,
          mortality_threshold = result$cessation_threshold,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          individual_ID = ind_id,
          species = species,
          potential_mortality_date = as.Date(NA),
          days_tracked_after = NA,
          baseline_distance = result$baseline_distance,
          baseline_speed = result$baseline_speed,
          mean_distance_after = NA,
          mean_speed_after = NA,
          reduction_distance_pct = NA,
          reduction_speed_pct = NA,
          gps_noise_threshold = result$gps_noise_threshold,
          movement_vs_gps_noise = NA,
          mortality_threshold = result$cessation_threshold,
          stringsAsFactors = FALSE
        )
      }
    }
  })
)

# Save compiled metrics
if (nrow(all_daily_metrics) > 0) {
  write.csv(all_daily_metrics,
            file.path(paths$tables, "missing_individuals_daily_movement_metrics.csv"),
            row.names = FALSE)
  
  write.csv(mortality_summary,
            file.path(paths$tables, "missing_individuals_mortality_dates.csv"),
            row.names = FALSE)
  
  message(sprintf("\nMovement metrics complete: %d individuals processed",
                  length(movement_plot_results)))
  
  # Summary statistics
  summary_stats <- all_daily_metrics %>%
    group_by(individual_ID, species) %>%
    summarise(
      n_days = n(),
      mean_daily_distance = mean(total_distance, na.rm = TRUE),
      sd_daily_distance = sd(total_distance, na.rm = TRUE),
      median_daily_distance = median(total_distance, na.rm = TRUE),
      mean_speed = mean(mean_speed, na.rm = TRUE),
      sd_speed = sd(mean_speed, na.rm = TRUE),
      first_observation = min(Date, na.rm = TRUE),
      last_observation = max(Date, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      mortality_summary %>%
        select(individual_ID, potential_mortality_date, movement_vs_gps_noise,
               gps_noise_threshold, mortality_threshold),
      by = "individual_ID"
    )
  
  write.csv(summary_stats,
            file.path(paths$tables, "missing_individuals_movement_summary.csv"),
            row.names = FALSE)
  
  message("\n=== Movement Summary Statistics ===")
  print(summary_stats %>%
          select(individual_ID, species, n_days, mean_daily_distance,
                 potential_mortality_date, movement_vs_gps_noise))
  
  # Mortality detection summary
  message("\n=== Mortality Detection Summary ===")
  mortality_detected <- sum(!is.na(mortality_summary$potential_mortality_date))
  message(sprintf("Individuals with detected mortality signal: %d / %d (%.1f%%)",
                  mortality_detected,
                  nrow(mortality_summary),
                  100 * mortality_detected / nrow(mortality_summary)))
  
  if (mortality_detected > 0) {
    message("\nDetected mortality events:")
    mortality_table <- mortality_summary %>%
      filter(!is.na(potential_mortality_date)) %>%
      arrange(potential_mortality_date) %>%
      select(individual_ID, species, potential_mortality_date,
             reduction_distance_pct, movement_vs_gps_noise) %>%
      mutate(
        interpretation = case_when(
          movement_vs_gps_noise < 1.0 ~ "Very likely dead (< GPS noise)",
          movement_vs_gps_noise < 2.0 ~ "Possibly dead (< 2x GPS noise)",
          TRUE ~ "Uncertain (> 2x GPS noise)"
        )
      )
    
    print(mortality_table)
    
    # Summary by interpretation
    message("\nMortality confidence breakdown:")
    interp_summary <- mortality_table %>%
      group_by(interpretation) %>%
      summarise(
        n = n(),
        mean_reduction = mean(reduction_distance_pct, na.rm = TRUE),
        .groups = "drop"
      )
    print(interp_summary)
  } else {
    message("No mortality signals detected above GPS noise threshold")
  }
  
  # Species comparison
  message("\n=== Species Comparison ===")
  species_summary <- mortality_summary %>%
    group_by(species) %>%
    summarise(
      n_individuals = n(),
      n_mortality_detected = sum(!is.na(potential_mortality_date)),
      pct_mortality_detected = 100 * mean(!is.na(potential_mortality_date)),
      mean_baseline_distance = mean(baseline_distance, na.rm = TRUE),
      mean_gps_threshold = mean(gps_noise_threshold, na.rm = TRUE),
      .groups = "drop"
    )
  print(species_summary)
  
  # Temporal pattern
  if (mortality_detected > 0) {
    message("\n=== Temporal Pattern of Mortality Events ===")
    mortality_timeline <- mortality_summary %>%
      filter(!is.na(potential_mortality_date)) %>%
      arrange(potential_mortality_date) %>%
      mutate(
        week = floor(as.numeric(potential_mortality_date - min(potential_mortality_date)) / 7) + 1
      ) %>%
      group_by(week, species) %>%
      summarise(n_events = n(), .groups = "drop")
    
    print(mortality_timeline)
  }
  
  # Save detailed mortality summary
  if (mortality_detected > 0) {
    detailed_mortality <- mortality_summary %>%
      filter(!is.na(potential_mortality_date)) %>%
      arrange(potential_mortality_date) %>%
      mutate(
        interpretation = case_when(
          movement_vs_gps_noise < 1.0 ~ "Very likely dead",
          movement_vs_gps_noise < 2.0 ~ "Possibly dead",
          TRUE ~ "Uncertain"
        ),
        days_since_start = as.numeric(potential_mortality_date -
                                        min(potential_mortality_date))
      )
    
    write.csv(detailed_mortality,
              file.path(paths$tables, "detailed_mortality_events.csv"),
              row.names = FALSE)
  }
  
} else {
  message("\nWarning: No movement metrics were successfully calculated")
}

message("\n=== Movement Metrics Analysis Complete ===")
message(sprintf("Output directory: %s", movement_plots_dir))
message("Generated files:")
message("  - Individual movement plots (distance and speed)")
message("  - missing_individuals_daily_movement_metrics.csv")
message("  - missing_individuals_mortality_dates.csv")
message("  - missing_individuals_movement_summary.csv")
if (exists("detailed_mortality") && nrow(detailed_mortality) > 0) {
  message("  - detailed_mortality_events.csv")
}
message("========================================\n")


# =================================================================-
# 8. PROXIMITY ANALYSIS (OPTIONAL) ####
# =================================================================-
# 
# OPTIONAL: Uses ctmm::proximity() to test if prey-predator movements
# are statistically independent. Computationally expensive, so disabled
# by default. Enable by setting run_proximity_analysis = TRUE in params.
#
# Proximity ratio interpretation:
# - < 1: Animals closer than expected (potential predation/attraction)
# - = 1: Independent movement (null hypothesis)
# - > 1: Animals farther than expected (avoidance)
# - CI contains 1: Cannot reject independence
# =================================================================-

#If I want to conduct proximity analysis, I need to identify the stretch of time that 
#I believe a individual was predated, because the whole track could be non-independent (due to the size of the lakes)
#Focus only on individuals that were identified as being likely predated.
#This needs to be thought through as how to deal with the tracking data.

if (params$run_proximity_analysis) {
  message("\n=== PROXIMITY ANALYSIS (OPTIONAL COMPONENT) ===")
  
  ## 8.1 Identify Candidate Pairs ----
  candidate_pairs <- all_encounters %>%
    left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
    filter(
      found_alive == 0,
      (total_encounters_high_conf > 50 |
         total_encounters_probable > params$proximity_encounter_threshold)
    ) %>%
    arrange(desc(total_encounters_high_conf))
  
  message(sprintf("Identified %d candidate pairs for proximity analysis",
                  nrow(candidate_pairs)))
  
  if (nrow(candidate_pairs) > 0) {
    ## 8.2 Run Proximity Analysis ----
    all_proximity_results <- list()
    
    for (i in 1:nrow(candidate_pairs)) {
      prey_id <- candidate_pairs$Prey_ID[i]
      pred_id <- candidate_pairs$Pred_ID[i]
      species <- candidate_pairs$Species[i]
      
      if (i %% 10 == 0 || i == 1) {
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
    
    ## 8.3 Compile Proximity Results ----
    proximity_df <- do.call(rbind, lapply(all_proximity_results, function(x) {
      if (!is.null(x) && !is.na(x$proximity_ratio_est)) {
        data.frame(
          Prey_ID = x$prey_id,
          Pred_ID = x$pred_id,
          proximity_ratio_est = x$proximity_ratio_est,
          proximity_ratio_low = x$proximity_ratio_low,
          proximity_ratio_high = x$proximity_ratio_high,
          independent = x$independent,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))
    
    if (!is.null(proximity_df) && nrow(proximity_df) > 0) {
      # Merge with encounter data
      proximity_with_encounters <- proximity_df %>%
        left_join(candidate_pairs, by = c("Prey_ID", "Pred_ID"))
      
      # Save results
      saveRDS(all_proximity_results,
              file.path(paths$proximity, "proximity_results_detailed_tiered.rds"))
      saveRDS(proximity_with_encounters,
              file.path(paths$proximity, "proximity_results_with_encounters_tiered.rds"))
      write.csv(proximity_with_encounters,
                file.path(paths$proximity, "proximity_results_tiered.csv"),
                row.names = FALSE)
      
      message(sprintf("Proximity analysis complete: %d pairs analyzed",
                      nrow(proximity_df)))
      
      ## 8.4 Integrate Proximity with Predation Classification ----
      suspected_predation_enhanced <- suspected_predation_distance %>%
        left_join(
          proximity_with_encounters %>%
            select(Prey_ID, Pred_ID, proximity_ratio_est, proximity_ratio_low,
                   proximity_ratio_high, independent),
          by = c("Prey_ID", "Pike_ID" = "Pred_ID")
        ) %>%
        mutate(
          # Proximity suggests predation if ratio < 1 (closer than expected)
          proximity_suggests_predation = !is.na(independent) & !independent,
          
          # Enhanced scoring system
          combined_predation_score = case_when(
            # Very strong: distance + proximity + high encounters
            likely_predated == 1 & proximity_suggests_predation &
              num_days_high_conf_50 >= 2 ~ 4,
            
            # Strong evidence
            likely_predated == 1 & proximity_suggests_predation ~ 3,
            likely_predated == 1 & num_days_high_conf_50 >= 3 ~ 3,
            
            # Moderate evidence
            likely_predated == 1 ~ 2,
            proximity_suggests_predation ~ 2,
            num_days_high_conf_10 >= 3 ~ 2,
            
            # Weak evidence
            TRUE ~ 1
          ),
          
          # Classification
          classification = case_when(
            combined_predation_score == 4 ~ "Very high confidence",
            combined_predation_score == 3 ~ "High confidence",
            combined_predation_score == 2 ~ "Moderate confidence",
            TRUE ~ "Low confidence"
          ),
          
          # Encounter tier summary
          primary_encounter_tier = case_when(
            num_days_high_conf_50 >= 2 ~ "High-confidence (≤1.4m)",
            num_days_high_conf_10 >= 2 ~ "High-confidence (≤1.4m)",
            num_days_probable_confirmed_10 >= 2 ~ "Probable (1.4-2.8m)",
            TRUE ~ "Uncertain"
          )
        )
    } else {
      message("No successful proximity analyses completed")
      suspected_predation_enhanced <- suspected_predation_distance %>%
        mutate(
          proximity_suggests_predation = NA,
          combined_predation_score = case_when(
            likely_predated == 1 & num_days_high_conf_50 >= 3 ~ 3,
            likely_predated == 1 ~ 2,
            num_days_high_conf_10 >= 3 ~ 2,
            TRUE ~ 1
          ),
          classification = case_when(
            combined_predation_score == 3 ~ "High confidence",
            combined_predation_score == 2 ~ "Moderate confidence",
            TRUE ~ "Low confidence"
          ),
          primary_encounter_tier = case_when(
            num_days_high_conf_50 >= 2 ~ "High-confidence (≤1.4m)",
            num_days_high_conf_10 >= 2 ~ "High-confidence (≤1.4m)",
            num_days_probable_confirmed_10 >= 2 ~ "Probable (1.4-2.8m)",
            TRUE ~ "Uncertain"
          )
        )
    }
  } else {
    message("No candidate pairs for proximity analysis")
    suspected_predation_enhanced <- suspected_predation_distance %>%
      mutate(
        proximity_suggests_predation = NA,
        combined_predation_score = case_when(
          likely_predated == 1 & num_days_high_conf_50 >= 3 ~ 3,
          likely_predated == 1 ~ 2,
          num_days_high_conf_10 >= 3 ~ 2,
          TRUE ~ 1
        ),
        classification = case_when(
          combined_predation_score == 3 ~ "High confidence",
          combined_predation_score == 2 ~ "Moderate confidence",
          TRUE ~ "Low confidence"
        ),
        primary_encounter_tier = case_when(
          num_days_high_conf_50 >= 2 ~ "High-confidence (≤1.4m)",
          num_days_high_conf_10 >= 2 ~ "High-confidence (≤1.4m)",
          num_days_probable_confirmed_10 >= 2 ~ "Probable (1.4-2.8m)",
          TRUE ~ "Uncertain"
        )
      )
  }
  
} else {
  message("\n=== PROXIMITY ANALYSIS SKIPPED (disabled in params) ===")
  
  # Classification without proximity
  suspected_predation_enhanced <- suspected_predation_distance %>%
    mutate(
      proximity_suggests_predation = NA,
      combined_predation_score = case_when(
        likely_predated == 1 & num_days_high_conf_50 >= 3 ~ 3,
        likely_predated == 1 ~ 2,
        num_days_high_conf_10 >= 3 ~ 2,
        TRUE ~ 1
      ),
      classification = case_when(
        combined_predation_score == 3 ~ "High confidence",
        combined_predation_score == 2 ~ "Moderate confidence",
        TRUE ~ "Low confidence"
      ),
      primary_encounter_tier = case_when(
        num_days_high_conf_50 >= 2 ~ "High-confidence (≤1.4m)",
        num_days_high_conf_10 >= 2 ~ "High-confidence (≤1.4m)",
        num_days_probable_confirmed_10 >= 2 ~ "Probable (1.4-2.8m)",
        TRUE ~ "Uncertain"
      )
    )
}

# Save enhanced results
saveRDS(suspected_predation_enhanced,
        file.path(paths$encounters, "suspected_predation_enhanced_tiered.rds"))

# Classification summary
classification_summary <- suspected_predation_enhanced %>%
  group_by(Species, classification, primary_encounter_tier) %>%
  summarise(
    n_events = n(),
    mean_high_conf_days = mean(num_days_high_conf_10, na.rm = TRUE),
    mean_consecutive_days = mean(consecutive_days_high_conf, na.rm = TRUE),
    .groups = "drop"
  )

message("\nFinal Classification Summary:")
print(classification_summary)


# =================================================================-
# 9. VISUALIZATIONS ####
# =================================================================-

message("\n=== CREATING VISUALIZATIONS ===")

## 9.1 Encounter Type Distribution ----
p1 <- ggplot(combined_distances %>% filter(encounter_type != "no_encounter"),
             aes(x = encounter_type, fill = Species)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Distribution of Encounter Types",
    subtitle = sprintf("High-confidence (≤%.1fm), Probable (%.1f-%.1fm), Possible (%.1f-%.1fm)",
                       params$high_confidence_threshold,
                       params$high_confidence_threshold, params$probable_threshold,
                       params$probable_threshold, params$possible_threshold),
    x = "Encounter Type",
    y = "Count",
    fill = "Prey Species"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(paths$figures, "encounter_type_distribution.png"),
       p1, width = 10, height = 6, dpi = 300)

## 9.2 Predation Classification ----
p2 <- ggplot(suspected_predation_enhanced,
             aes(x = classification, fill = primary_encounter_tier)) +
  geom_bar(position = "stack", alpha = 0.8) +
  facet_wrap(~Species) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(
    title = "Predation Events by Confidence Level",
    x = "Classification",
    y = "Number of Events",
    fill = "Primary Encounter Type"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(paths$figures, "predation_classification.png"),
       p2, width = 12, height = 7, dpi = 300)

## 9.3 Distance Distribution ----
p3 <- ggplot(combined_distances %>%
               filter(encounter_type != "no_encounter") %>%
               sample_n(min(50000, n())),
             aes(x = est, fill = encounter_type)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = params$high_confidence_threshold,
             linetype = "dashed", color = "red") +
  geom_vline(xintercept = params$probable_threshold,
             linetype = "dashed", color = "orange") +
  scale_fill_manual(
    values = c("high_confidence" = "#d73027",
               "probable" = "#fee090",
               "possible" = "#91bfdb"),
    labels = c("High-confidence", "Probable", "Possible")
  ) +
  labs(
    title = "Distance Distribution by Encounter Type",
    subtitle = sprintf("Thresholds: %.1fm (high-conf), %.1fm (probable)",
                       params$high_confidence_threshold, params$probable_threshold),
    x = "Distance (m)",
    y = "Count",
    fill = "Encounter Type"
  ) +
  xlim(0, 5) +
  theme_classic()

ggsave(file.path(paths$figures, "distance_distribution_by_tier.png"),
       p3, width = 10, height = 6, dpi = 300)

## 9.4 Proximity Results (if available) ----
if (params$run_proximity_analysis && exists("proximity_with_encounters")) {
  
  p4 <- ggplot(proximity_with_encounters,
               aes(x = total_encounters_high_conf, y = proximity_ratio_est)) +
    geom_point(aes(color = independent), size = 3, alpha = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    scale_color_manual(
      values = c("TRUE" = "darkgreen", "FALSE" = "darkred"),
      labels = c("TRUE" = "Independent", "FALSE" = "Non-independent"),
      name = "Movement Pattern"
    ) +
    labs(
      title = "Proximity Ratio vs High-Confidence Encounters",
      subtitle = "Ratio < 1: closer than expected | Ratio > 1: farther than expected",
      x = "Total High-Confidence Encounters",
      y = "Proximity Ratio"
    ) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(paths$proximity, "proximity_ratio_vs_encounters.png"),
         p4, width = 10, height = 7, dpi = 300)
}

message("Visualizations complete!")


# =================================================================-
# 10. EXPORT RESULTS ####
# =================================================================-

message("\n=== EXPORTING RESULTS ===")

## 10.1 Create Excel Workbook ----
wb <- createWorkbook()

# Encounter summary
addWorksheet(wb, "Encounter Summary")
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
writeData(wb, "Encounter Summary", tier_summary)

# Daily encounters
addWorksheet(wb, "Daily Encounters")
writeData(wb, "Daily Encounters", prey_daily_encounters)

# Total encounters
addWorksheet(wb, "Total Encounters")
writeData(wb, "Total Encounters", all_encounters)

# Suspected predation
addWorksheet(wb, "Suspected Predation")
writeData(wb, "Suspected Predation", suspected_predation_enhanced)

# Classification summary
addWorksheet(wb, "Classification Summary")
writeData(wb, "Classification Summary", classification_summary)

# Proximity results (if available)
if (params$run_proximity_analysis && exists("proximity_with_encounters")) {
  addWorksheet(wb, "Proximity Results")
  writeData(wb, "Proximity Results", proximity_with_encounters)
}

#Save workbook
saveWorkbook(wb,
             file.path(paths$tables, "muddyfoot_predation_analysis_complete.xlsx"),
             overwrite = TRUE)
message("Excel workbook saved!")

# 10.2 Save Individual CSV Files ----
  write.csv(suspected_predation_enhanced,
            file.path(paths$tables, "suspected_predation_events.csv"),
            row.names = FALSE)
write.csv(all_encounters,
          file.path(paths$tables, "all_encounter_totals.csv"),
          row.names = FALSE)
write.csv(prey_daily_encounters,
          file.path(paths$tables, "daily_encounters.csv"),
          row.names = FALSE)

#=================================================================-
# 11. SUMMARY REPORT ####
#=================================================================-

  message("\n", paste(rep("=", 70), collapse = ""))
message("ANALYSIS COMPLETE - SUMMARY")
message(paste(rep("=", 70), collapse = ""))
Total observations
total_obs <- nrow(combined_distances)
message(sprintf("\nTotal observations: %s", format(total_obs, big.mark = ",")))

#Encounter distribution
message("\nEncounter distribution:")
enc_dist <- table(combined_distances$encounter_type)
for (type in names(enc_dist)) {
  pct <- 100 * enc_dist[type] / total_obs
  message(sprintf("  %s: %s (%.2f%%)",
                  type, format(enc_dist[type], big.mark = ","), pct))
}

#Suspected predation
message(sprintf("\nSuspected predation events: %d",
                nrow(suspected_predation_enhanced)))
message("\nClassification breakdown:")
class_table <- table(suspected_predation_enhanced$classification)
for (class in names(class_table)) {
  message(sprintf("  %s: %d", class, class_table[class]))
}

#Proximity analysis
if (params$run_proximity_analysis && exists("proximity_with_encounters")) {
  n_non_indep <- sum(!proximity_with_encounters$independent, na.rm = TRUE)
  pct_non_indep <- 100 * mean(!proximity_with_encounters$independent, na.rm = TRUE)
  message(sprintf("\nProximity analysis:"))
  message(sprintf("  Pairs analyzed: %d", nrow(proximity_with_encounters)))
  message(sprintf("  Non-independent movement: %d (%.1f%%)",
                  n_non_indep, pct_non_indep))
} else {
  message("\nProximity analysis: Not performed")
}

#Output locations
message("\nOutput locations:")
message(sprintf("  Figures: %s", paths$figures))
message(sprintf("  Tables: %s", paths$tables))
message(sprintf("  Encounters: %s", paths$encounters))
if (params$create_distance_plots) {
  message(sprintf("  Distance plots: %s", paths$distance_plots))
}
if (params$run_proximity_analysis) {
  message(sprintf("  Proximity: %s", paths$proximity))
}
message("\n", paste(rep("=", 70), collapse = ""))

#=================================================================-
# 12. SESSION INFO ####
#=================================================================-

sink(file.path(paths$tables, "session_info.txt"))
cat("Predation Analysis Session Info\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")
cat("Analysis Parameters:\n")
cat(sprintf("  GPS error: %.2fm\n", params$gps_error_sd))
cat(sprintf("  Fixes per day: %d\n", params$fixes_per_day))
cat(sprintf("  High-confidence threshold: %.1fm\n", params$high_confidence_threshold))
cat(sprintf("  Probable threshold: %.1fm\n", params$probable_threshold))
cat(sprintf("  Possible threshold: %.1fm\n", params$possible_threshold))
cat(sprintf("  Proximity analysis: %s\n",
            ifelse(params$run_proximity_analysis, "Enabled", "Disabled")))
cat(sprintf("  Distance plots: %s\n",
            ifelse(params$create_distance_plots, "Enabled", "Disabled")))
cat("\n")
print(sessionInfo())
sink()
message("\nAll analyses complete!")
message("\nKey files to review:")
message("  - muddyfoot_predation_analysis_complete.xlsx")
message("  - suspected_predation_events.csv")
message("  - encounter_type_distribution.png")
message("  - predation_classification.png")
if (params$create_distance_plots) {
  message("  - Individual distance plots in: distance_plots/")
}
if (exists("detailed_mortality") && nrow(detailed_mortality) > 0) {
  message("  - detailed_mortality_events.csv")
  message("  - Individual movement plots in: movement_metrics/")
}