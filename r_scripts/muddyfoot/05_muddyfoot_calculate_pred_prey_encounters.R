# =================================================================-
# Predator-Prey Interactions and encounter rate - Muddyfoot ####
# =================================================================-
# 
# ANALYSIS OVERVIEW:
# This script identifies potential predation events by analyzing movement 
# tracking data for prey (roach/perch) and predators (pike), and calculates
# encounter summaries used by the downstream encounter rate analysis script.

# METHODOLOGY:
# 1. TIERED ENCOUNTER DETECTION based on GPS error (0.5m)
#    - High-confidence: ≤1.4m (within 95% combined GPS error)
#    - Probable:        1.4-2.8m (requires temporal confirmation)
#
# 2. MOVEMENT CESSATION DETECTION
#
## OUTPUTS USED DOWNSTREAM:
#   data/encounters/muddyfoot/muddyfoot_pike_roach_distances_tiered_df.rds
#   data/encounters/muddyfoot/muddyfoot_pike_perch_distances_tiered_df.rds
#   -> Read by 07_muddyfoot_encounter_rate_analysis.R
#
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
  filtered_data  = "./data/tracks_filtered/muddyfoot/",
  lake_polygon   = "./data/lake_params/polygons/",
  ctmm           = "./data/ctmm_fits/",
  telem          = "./data/telem_obj/muddyfoot/",
  encounters     = "./data/encounters/muddyfoot/",
  size           = "./data/fish_biometrics/",
  tables         = "./tables/muddyfoot/",
  figures        = "./figures/muddyfoot/",
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
  
  # Encounter count thresholds (passed to identify_predation_events)
  encounter_threshold_low = 10,        # Min high-conf encounters/day to flag a day
  encounter_threshold_high = 50,       # High encounter threshold (stricter flag)
  consecutive_days_threshold = 2,      # Consecutive days required for predation signal
  
  # Distance plot settings
  create_distance_plots = TRUE,        # Create time series plots
  distance_plot_min_encounters = 50,   # Min encounters to generate plot
  
  # Mortality detection settings
  gps_error_sd = 0.5,                  # Single fish GPS error (meters)
  fixes_per_day = 2880,                # 30-second fix intervals
  mortality_cessation_threshold = 0.15, # 15% of baseline movement
  mortality_consecutive_days = 2       # Consecutive low-movement days
)

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
      
      # Categorical encounter type
      df$encounter_type <- dplyr::case_when(
        df$est <= high_conf_threshold ~ "high_confidence",
        df$est <= probable_threshold ~ "probable",
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
                                  "high_conf_threshold", "probable_threshold", "combinations"),
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




#' Identify Suspected Predation Events
#'
#' Detects sustained contact periods using multi-criteria approach.
#' Thresholds scaled from BT (sigma=0.379m) to Muddyfoot (sigma=0.5m):
#'   Scale factor = combined_error_muddyfoot / combined_error_BT
#'                = (sqrt(2)*0.5) / (sqrt(2)*0.379) = 0.707 / 0.536 = 1.32
#'
#' Sustained-contact distance cutoffs (* 1.32 from BT):
#'   BT <0.5m  -> Muddyfoot <0.66m  |  BT <0.1m  -> Muddyfoot <0.13m
#'   BT <1.0m  -> Muddyfoot <1.32m  |  BT <1.5m  -> Muddyfoot <2.0m
#'   BT <2.0m  -> Muddyfoot <2.64m
#' "At-zero" threshold = combined GPS error = 0.707m
#'
#' @param distances_df Data frame from calculate_pairwise_distances()
#' @param enc_threshold_low Min encounters per day for high-conf flag
#' @param enc_threshold_high High encounter threshold
#' @param min_dist_threshold Distance threshold (m)
#' @param consec_days_threshold Consecutive days required
#'
#' @return Data frame with predation likelihood per prey-predator pair
identify_predation_events <- function(distances_df,
                                      enc_threshold_low      = 10,
                                      enc_threshold_high     = 50,
                                      min_dist_threshold     = 0.66,
                                      consec_days_threshold  = 2,
                                      probable_threshold_low = 50,
                                      probable_threshold_high = 100,
                                      final_days_window      = 2) {
  
  # Daily summaries
  daily_summary <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      n_encounters_high_conf = sum(encounter_type == "high_confidence", na.rm = TRUE),
      n_encounters_probable  = sum(encounter_type == "probable",        na.rm = TRUE),
      min_distance           = min(est, na.rm = TRUE),
      mean_distance          = mean(est, na.rm = TRUE),
      # "at zero" = within combined GPS error for Muddyfoot (0.707m)
      n_at_zero              = sum(est < 0.707, na.rm = TRUE),
      pct_at_zero            = n_at_zero / n() * 100,
      .groups = "drop"
    ) %>%
    arrange(Prey_ID, Pred_ID, Date)
  
  # Sustained contact flagging (thresholds scaled * 1.32 from BT)
  predation_windows <- daily_summary %>%
    group_by(Prey_ID, Pred_ID) %>%
    arrange(Date) %>%
    mutate(
      # A day qualifies as sustained contact only if distances are CONSISTENTLY
      # small — not just briefly. mean_distance < 4.0m rules out days where the
      # prey-predator distance swings widely (low min but high mean), which
      # indicate the prey is still moving freely and has not been predated.
      # This prevents brief close approaches within an otherwise large-distance
      # day from triggering a false sustained-contact window.
      sustained_contact = (
        (min_distance < 0.66 & mean_distance < 4.0 & n_encounters_high_conf > 10) |
          (min_distance < 0.13 & mean_distance < 4.0 & n_encounters_high_conf > 5)  |
          (pct_at_zero > 40    & mean_distance < 4.0 & n_encounters_high_conf > 5)  |
          (min_distance < 1.32 & mean_distance < 4.0 & n_encounters_high_conf > 50) |
          (min_distance < 2.0  & mean_distance < 4.0 & n_encounters_high_conf > 100)|
          (min_distance < 2.0  & mean_distance < 4.0 & n_encounters_high_conf > 20  & n_encounters_probable > 50) |
          (min_distance < 2.64 & mean_distance < 4.0 & n_encounters_probable > 150)
      ),
      
      day_diff = as.numeric(difftime(Date, lag(Date), units = "days")),
      is_consecutive_contact = (day_diff == 1 & sustained_contact & lag(sustained_contact)),
      contact_group = cumsum(!is_consecutive_contact | is.na(is_consecutive_contact))
    ) %>%
    group_by(Prey_ID, Pred_ID, contact_group) %>%
    mutate(
      consecutive_contact_days = if_else(sustained_contact, n(), 0L),
      contact_period_start     = if_else(sustained_contact, min(Date), as.Date(NA)),
      contact_period_end       = if_else(sustained_contact, max(Date), as.Date(NA))
    ) %>%
    ungroup()
  
  # Consecutive encounter days
  consecutive_info <- predation_windows %>%
    group_by(Prey_ID, Pred_ID) %>%
    arrange(Date) %>%
    mutate(
      high_conf_flag      = n_encounters_high_conf >= enc_threshold_low,
      is_consecutive      = (day_diff == 1 & high_conf_flag & lag(high_conf_flag)),
      consec_group        = cumsum(!is_consecutive | is.na(is_consecutive)),
      probable_flag       = n_encounters_probable >= probable_threshold_low,
      is_consecutive_prob = (day_diff == 1 & probable_flag & lag(probable_flag)),
      consec_group_prob   = cumsum(!is_consecutive_prob | is.na(is_consecutive_prob))
    ) %>%
    group_by(Prey_ID, Pred_ID, consec_group) %>%
    mutate(consecutive_days_high_conf = if_else(high_conf_flag, n(), 0L)) %>%
    group_by(Prey_ID, Pred_ID, consec_group_prob) %>%
    mutate(consecutive_days_probable  = if_else(probable_flag,  n(), 0L)) %>%
    ungroup()
  
  # Aggregate by pair
  consecutive_info %>%
    group_by(Prey_ID, Pred_ID) %>%
    summarise(
      total_encounters_high_conf        = sum(n_encounters_high_conf, na.rm = TRUE),
      num_days_high_conf_10             = sum(n_encounters_high_conf >= enc_threshold_low,  na.rm = TRUE),
      num_days_high_conf_50             = sum(n_encounters_high_conf >= enc_threshold_high, na.rm = TRUE),
      num_days_high_conf_100            = sum(n_encounters_high_conf >= 100, na.rm = TRUE),
      max_consecutive_days_high_conf    = max(consecutive_days_high_conf, na.rm = TRUE),
      total_encounters_probable         = sum(n_encounters_probable,  na.rm = TRUE),
      num_days_probable_50              = sum(n_encounters_probable >= probable_threshold_low,  na.rm = TRUE),
      num_days_probable_100             = sum(n_encounters_probable >= probable_threshold_high, na.rm = TRUE),
      num_days_probable_200             = sum(n_encounters_probable >= 200, na.rm = TRUE),
      max_consecutive_days_probable     = max(consecutive_days_probable, na.rm = TRUE),
      num_days_sustained_contact        = sum(sustained_contact, na.rm = TRUE),
      max_consecutive_contact_days      = max(consecutive_contact_days, na.rm = TRUE),
      sustained_contact_start           = { v <- contact_period_start[!is.na(contact_period_start)]; if (length(v)) min(v) else as.Date(NA) },
      sustained_contact_end             = { v <- contact_period_end[!is.na(contact_period_end)];     if (length(v)) max(v) else as.Date(NA) },
      # "at zero" = within combined GPS error (0.707m for Muddyfoot)
      total_points_at_zero              = sum(n_at_zero, na.rm = TRUE),
      overall_min_distance              = min(min_distance,  na.rm = TRUE),
      mean_of_daily_min_distance        = mean(min_distance, na.rm = TRUE),
      # Distance thresholds scaled * 1.32 from BT
      num_days_under_0.66m              = sum(min_distance < 0.66, na.rm = TRUE),
      num_days_under_1.32m              = sum(min_distance < 1.32, na.rm = TRUE),
      num_days_under_2.0m               = sum(min_distance < 2.0,  na.rm = TRUE),
      first_encounter_date              = min(Date, na.rm = TRUE),
      last_encounter_date               = max(Date, na.rm = TRUE),
      tracking_duration_days            = as.numeric(difftime(max(Date), min(Date), units = "days")) + 1,
      total_days_with_encounters        = n_distinct(Date),
      encounters_final_2days_high_conf  = sum(n_encounters_high_conf[Date >= (max(Date) - final_days_window)], na.rm = TRUE),
      encounters_final_2days_probable   = sum(n_encounters_probable[ Date >= (max(Date) - final_days_window)], na.rm = TRUE),
      high_conf_intensity               = total_encounters_high_conf / tracking_duration_days,
      probable_intensity                = total_encounters_probable  / tracking_duration_days,
      .groups = "drop"
    ) %>%
    mutate(
      combined_encounter_score = (total_encounters_high_conf * 1.0) +
        (total_encounters_probable * 0.5),
      
      # Predation likelihood — mirrors BT/Cow logic with Muddyfoot distance values
      likely_predated = case_when(
        (num_days_sustained_contact >= 3 & max_consecutive_contact_days >= 2) ~ "very_likely",
        (max_consecutive_contact_days >= 5)                                   ~ "very_likely",
        (total_points_at_zero > 100 & num_days_sustained_contact >= 2)       ~ "very_likely",
        (num_days_sustained_contact >= 2 & overall_min_distance < 0.66)      ~ "likely",
        (num_days_high_conf_50 > 1      & overall_min_distance < 0.66)       ~ "likely",
        (num_days_probable_100 > 1      & overall_min_distance < 0.99)       ~ "likely",
        num_days_sustained_contact >= 2                                        ~ "likely",
        (num_days_high_conf_10 > 2      & total_encounters_probable > 250)   ~ "possible",
        (num_days_probable_50 > 1       & overall_min_distance < 1.32)       ~ "possible",
        (total_encounters_high_conf > 200 & num_days_under_1.32m > 0)        ~ "possible",
        num_days_sustained_contact >= 1                                        ~ "possible",
        (total_encounters_high_conf > 150 | total_encounters_probable > 400) ~ "weak_signal",
        TRUE                                                                   ~ "unlikely"
      )
    )
}


#' Plot Distance Time Series for Prey-Predator Pair
#'
#' Creates time series of separation distance with:
#'   - Tiered threshold lines (high-confidence and probable only)
#'   - Highlighted sustained-contact periods (potential predation windows)
#' This mirrors the BT and Cow Paradise plot style.
#'
#' @param distances_df Data frame from calculate_pairwise_distances()
#' @param prey_id Prey individual ID
#' @param pred_id Predator individual ID
#' @param species Species name
#' @param output_dir Directory to save plot (NULL = no save)
#' @param threshold_lines Add threshold lines?
#' @param predation_summary Output from identify_predation_events() for this species
#'
#' @return ggplot object
plot_distance_timeseries <- function(distances_df, prey_id, pred_id, species,
                                     output_dir = NULL, threshold_lines = TRUE,
                                     predation_summary = NULL) {
  
  pair_data <- distances_df %>%
    filter(Prey_ID == prey_id, Pred_ID == pred_id) %>%
    arrange(timestamp)
  
  if (nrow(pair_data) == 0) {
    warning(sprintf("No data found for %s - %s", prey_id, pred_id))
    return(NULL)
  }
  
  # Look up predation info for this pair
  predation_info <- NULL
  if (!is.null(predation_summary)) {
    col_check <- if ("Pike_ID" %in% colnames(predation_summary)) "Pike_ID" else "Pred_ID"
    predation_info <- predation_summary %>%
      filter(Prey_ID == prey_id, .data[[col_check]] == pred_id)
  }
  
  n_encounters_high <- sum(pair_data$encounter_high_conf)
  n_encounters_prob <- sum(pair_data$encounter_probable)
  mean_dist         <- mean(pair_data$est, na.rm = TRUE)
  min_dist          <- min(pair_data$est,  na.rm = TRUE)
  
  # X-axis breaks (every 2 days)
  date_range  <- range(pair_data$Date, na.rm = TRUE)
  date_breaks <- seq(from = date_range[1], to = date_range[2], by = "day")
  date_breaks <- date_breaks[seq(1, length(date_breaks), by = 2)]
  
  p <- ggplot(pair_data, aes(x = timestamp, y = est)) +
    geom_line(color = "gray50", alpha = 0.5) +
    geom_point(aes(color = encounter_type), size = 1, alpha = 0.6) +
    scale_color_manual(
      values = c("high_confidence" = "#d62728", "probable" = "#ff7f0e",
                 "no_encounter" = "gray80"),
      labels = c("high_confidence" = "High Confidence (<=1.4m)",
                 "probable"        = "Probable (1.4-2.8m)",
                 "no_encounter"    = "No Encounter"),
      name = "Encounter Type"
    ) +
    labs(
      title    = sprintf("%s: %s vs Pike %s", species, prey_id, pred_id),
      subtitle = if (!is.null(predation_info) && nrow(predation_info) > 0) {
        sprintf("Predation likelihood: %s | Sustained contact days: %d | High-conf: %d | Probable: %d",
                predation_info$likely_predated,
                predation_info$num_days_sustained_contact,
                n_encounters_high, n_encounters_prob)
      } else {
        sprintf("High-conf encounters: %d | Probable: %d | Mean dist: %.1fm | Min dist: %.2fm",
                n_encounters_high, n_encounters_prob, mean_dist, min_dist)
      },
      x = "Date",
      y = "Separation distance (m)"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "gray30"),
      axis.text     = element_text(size = 9),
      axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
  if (threshold_lines) {
    p <- p +
      geom_hline(yintercept = 0.45, linetype = "dashed", color = "#d62728", alpha = 0.7) +
      geom_hline(yintercept = 1.4,  linetype = "dashed", color = "#ff7f0e", alpha = 0.7) +
      annotate("text", x = min(pair_data$timestamp), y = 0.45,
               label = "0.45m (strike dist)", hjust = 0, vjust = -0.5, size = 3, color = "#d62728") +
      annotate("text", x = min(pair_data$timestamp), y = 1.4,
               label = "1.4m (high-conf limit)", hjust = 0, vjust = -0.5, size = 3, color = "#ff7f0e")
  }
  
  # Highlight sustained contact period (potential predation window)
  if (!is.null(predation_info) && nrow(predation_info) > 0 &&
      !is.na(predation_info$sustained_contact_start) &&
      !is.na(predation_info$sustained_contact_end) &&
      !is.infinite(predation_info$sustained_contact_start) &&
      !is.infinite(predation_info$sustained_contact_end) &&
      predation_info$num_days_sustained_contact > 0) {
    
    p <- p +
      annotate("rect",
               xmin = as.POSIXct(predation_info$sustained_contact_start),
               xmax = as.POSIXct(predation_info$sustained_contact_end + 1),
               ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15) +
      annotate("text",
               x     = as.POSIXct(predation_info$sustained_contact_end),
               y     = max(pair_data$est, na.rm = TRUE) * 0.9,
               label = "Sustained\nContact\n(Likely Predation)",
               hjust = 1, vjust = 1, size = 3.5, fontface = "bold", color = "#d62728")
  }
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename <- sprintf("%s_%s_vs_%s_distances.png", species, prey_id, pred_id)
    ggsave(file.path(output_dir, filename), plot = p, width = 10, height = 6, dpi = 300)
    message(sprintf("  Saved distance plot: %s", filename))
  }
  
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
    parallel = TRUE,
    n_cores = 16
  )
  saveRDS(perch_pike_distances_df, perch_pike_file)
}

message(sprintf("Loaded distances: %d roach-pike, %d perch-pike observations",
                nrow(roach_pike_distances_df), nrow(perch_pike_distances_df)))


perch_pike_distances_df <- readRDS(paste0(paths$encounters, "muddyfoot_pike_perch_distances_tiered_df.rds"))
roach_pike_distances_df <- readRDS(paste0(paths$encounters, "muddyfoot_pike_roach_distances_tiered_df.rds"))

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
# Groups by prey-predator pair and date, counting encounters in each tier
# and calculating daily distance metrics for each pair x day combination.
roach_daily <- roach_pike_distances_df %>%
  group_by(Prey_ID, Pred_ID, Date) %>%
  summarise(
    Species                   = first(Species),
    encounter_count_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    encounter_count_probable  = sum(encounter_probable,  na.rm = TRUE),
    daily_avg_dist            = mean(est, na.rm = TRUE),
    daily_min_dist            = min(est,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Pike_ID = Pred_ID)

perch_daily <- perch_pike_distances_df %>%
  group_by(Prey_ID, Pred_ID, Date) %>%
  summarise(
    Species                   = first(Species),
    encounter_count_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    encounter_count_probable  = sum(encounter_probable,  na.rm = TRUE),
    daily_avg_dist            = mean(est, na.rm = TRUE),
    daily_min_dist            = min(est,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Pike_ID = Pred_ID)

#combine
prey_daily_encounters <- bind_rows(roach_daily, perch_daily)


## 5.2 Calculate Total Encounters ----
# Aggregates across all observations for each prey-predator pair, giving
# total encounter counts, distance metrics, and temporal spread (n_days).
roach_total <- roach_pike_distances_df %>%
  group_by(Prey_ID, Pred_ID, Species) %>%
  summarise(
    total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    total_encounters_probable  = sum(encounter_probable,  na.rm = TRUE),
    total_encounters_legacy    = sum(encounter,           na.rm = TRUE),
    min_distance               = min(est, na.rm = TRUE),
    mean_distance              = mean(est, na.rm = TRUE),
    n_days                     = n_distinct(Date),
    n_days_high_conf           = n_distinct(Date[encounter_high_conf == 1]),
    n_days_probable            = n_distinct(Date[encounter_probable  == 1]),
    .groups = "drop"
  ) %>%
  arrange(desc(total_encounters_high_conf))

perch_total <- perch_pike_distances_df %>%
  group_by(Prey_ID, Pred_ID, Species) %>%
  summarise(
    total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    total_encounters_probable  = sum(encounter_probable,  na.rm = TRUE),
    total_encounters_legacy    = sum(encounter,           na.rm = TRUE),
    min_distance               = min(est, na.rm = TRUE),
    mean_distance              = mean(est, na.rm = TRUE),
    n_days                     = n_distinct(Date),
    n_days_high_conf           = n_distinct(Date[encounter_high_conf == 1]),
    n_days_probable            = n_distinct(Date[encounter_probable  == 1]),
    .groups = "drop"
  ) %>%
  arrange(desc(total_encounters_high_conf))

#combine
all_encounters <- bind_rows(roach_total, perch_total)

total_encounters_ID <- all_encounters %>% 
  group_by(Prey_ID, Species) %>%
  summarise(
    total_encounters_high_conf = sum(total_encounters_high_conf),
    total_encounters_probable = sum(total_encounters_probable),
    n_predators = n_distinct(Pred_ID),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_encounters_high_conf))

# Add treatment information
total_encounters_ID <- total_encounters_ID %>%
  left_join(
    muddyfoot_filt_data %>% 
      select(individual_ID, treatment) %>% 
      distinct(),
    by = c("Prey_ID" = "individual_ID")
  )

saveRDS(total_encounters_ID, paste0(paths$encounters, "muddyfoot_total_encounters_ID.rds"))

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
  enc_threshold_low     = params$encounter_threshold_low,
  enc_threshold_high    = params$encounter_threshold_high,
  consec_days_threshold = params$consecutive_days_threshold
) %>% rename(Pike_ID = Pred_ID) %>% mutate(Species = "Roach")

## 6.2 Perch Predation Events ----
perch_predation <- identify_predation_events(
  perch_pike_distances_df,
  enc_threshold_low     = params$encounter_threshold_low,
  enc_threshold_high    = params$encounter_threshold_high,
  consec_days_threshold = params$consecutive_days_threshold
) %>% rename(Pike_ID = Pred_ID) %>% mutate(Species = "Perch")

## 6.3 Combine and Filter ----
# Less restrictive filter to capture brief predation events where a pike may have
# regurgitated the prey tag shortly after ingestion — these may appear as only a
# short window of near-zero distances rather than sustained multi-day contact.
suspected_predation_distance <- bind_rows(roach_predation, perch_predation) %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
  filter(
    found_alive == 0,
    (num_days_high_conf_50 > 1 | num_days_probable_50 > 1 | num_days_probable_100 > 0 |
       total_encounters_high_conf > 200 | total_encounters_probable > 500 |
       combined_encounter_score > 400   | overall_min_distance < 0.99 |
       num_days_under_1.32m > 2         | encounters_final_2days_high_conf > 100 |
       encounters_final_2days_probable > 200 |
       (num_days_high_conf_10 > 2 & total_encounters_probable > 250) |
       likely_predated %in% c("very_likely", "likely", "possible"))
  ) %>%
  arrange(desc(likely_predated), desc(combined_encounter_score))

# Print summary
predation_summary_table <- suspected_predation_distance %>%
  group_by(Species, likely_predated) %>%
  summarise(
    n_cases                    = n(),
    mean_high_conf_encounters  = mean(total_encounters_high_conf, na.rm = TRUE),
    mean_probable_encounters   = mean(total_encounters_probable,  na.rm = TRUE),
    mean_min_distance          = mean(overall_min_distance,       na.rm = TRUE),
    mean_sustained_days        = mean(num_days_sustained_contact, na.rm = TRUE),
    .groups = "drop"
  )

message("\nPredation Summary:")
print(predation_summary_table)


# =================================================================-
# 7. DISTANCE TIME SERIES PLOTS  ####
# =================================================================-
# 
# Creates time series plots showing prey-predator separation distances
# over time. Annotates key dates (first 50+ encounters, max encounters).
# =================================================================-

if (params$create_distance_plots) {
  message("\n=== CREATING DISTANCE TIME SERIES PLOTS ===")
  
  high_encounter_pairs <- suspected_predation_distance %>%
    filter(num_days_high_conf_10 >= 2 | num_days_high_conf_50 >= 1 |
             num_days_sustained_contact >= 1) %>%
    select(Prey_ID, Pike_ID, Species)
  
  message(sprintf("Creating distance plots for %d pairs", nrow(high_encounter_pairs)))
  
  for (i in 1:nrow(high_encounter_pairs)) {
    prey_id  <- high_encounter_pairs$Prey_ID[i]
    pred_id  <- high_encounter_pairs$Pike_ID[i]
    species  <- high_encounter_pairs$Species[i]
    
    dist_data <- if (species == "Roach") roach_pike_distances_df else perch_pike_distances_df
    pred_data <- if (species == "Roach") roach_predation         else perch_predation
    
    plot_distance_timeseries(
      distances_df      = dist_data,
      prey_id           = prey_id,
      pred_id           = pred_id,
      species           = species,
      output_dir        = paths$distance_plots,
      threshold_lines   = TRUE,
      predation_summary = pred_data
    )
  }
  
  message(sprintf("Distance plots complete: %d plots created", nrow(high_encounter_pairs)))
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