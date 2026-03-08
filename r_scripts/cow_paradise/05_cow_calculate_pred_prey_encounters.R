# =================================================================-
# Predator-Prey Interaction Analysis - COW PARADISE ####
# =================================================================-
# Identifies potential predation events using GPS tracking data
# Key approach: Sustained contact detection (distances ~0m for multiple days)
#
# IMPORTANT DIFFERENCE FROM BT LAKE:
# GPS error is substantially higher in Cow Paradise (σ = 0.827m vs 0.379m in BT)
# This means:
#   - Combined pairwise GPS error: √(0.827² + 0.827²) = 1.17m (vs 0.536m in BT)
#   - Distance thresholds are scaled up accordingly (~2.2x BT values)
#   - High-confidence: ≤3.0m (vs ≤1.4m in BT)
#   - Probable:        3.0–6.0m (vs 1.4–2.8m in BT)
#   - Possible:        6.0–9.0m (vs 2.8–4.2m in BT)
#   - GPS noise (random walk): ~80m/day (vs ~30m/day in BT)
#   - Mortality threshold: ~200m/day (vs ~75m/day in BT)
#   - Sustained contact criteria are relaxed proportionally
#
# Author: Marcus Michelangeli | Date: 2024
# =================================================================-


# =================================================================-
# 1. SETUP ####
# =================================================================-

## Libraries ----
suppressPackageStartupMessages({
  library(dplyr); library(ctmm); library(lubridate); library(data.table)
  library(sf); library(parallel); library(ggplot2); library(openxlsx); library(patchwork)
})
Sys.setenv(TZ = 'Europe/Stockholm')

## Paths ----
paths <- list(
  filtered_data  = "./data/tracks_filtered/cow_paradise/",
  ctmm           = "./data/ctmm_fits/",
  telem          = "./data/telem_obj/cow_paradise/",
  encounters     = "./data/encounters/cow_paradise/",
  size           = "./data/fish_biometrics/",
  tables         = "./tables/cow_paradise/",
  figures        = "./figures/cow_paradise/",
  distance_plots = "./figures/cow_paradise/distance_plots/"
)
for (path in paths) if (!dir.exists(path)) dir.create(path, recursive = TRUE)

## Parameters ----
# GPS error for Cow Paradise is 0.827m (vs 0.379m for BT).
# Combined pairwise error = √2 × 0.827 = 1.170m.
# All distance thresholds are scaled by the ratio of combined errors:
#   scale_factor ≈ 1.170 / 0.536 ≈ 2.18
# Thresholds are rounded to ecologically sensible values.
params <- list(
  gps_error_sd                  = 0.827,   # GPS error SD (meters) — Cow Paradise specific
  fixes_per_day                 = 2880,    # 30-second fixes
  
  # Tiered distance thresholds (scaled from BT by ~2.18x)
  high_confidence_threshold     = 3.0,     # ≤3.0m: Within 95% combined GPS error
  probable_threshold            = 6.0,     # ≤6.0m: 2x combined GPS error
  possible_threshold            = 9.0,     # ≤9.0m: 3x combined GPS error
  strike_distance               = 0.45,    # Pike strike distance (unchanged)
  
  # Encounter count thresholds (same logic as BT — counts are scale-independent)
  encounter_threshold_low       = 10,      # Min encounters/day
  encounter_threshold_high      = 50,      # High encounter threshold
  
  # Distance-based thresholds (scaled from BT)
  min_distance_threshold        = 3.0,     # Min distance for high-confidence flag
  
  # Sustained contact: same consecutive-day logic as BT
  consecutive_days_threshold    = 2,       # Days required for predation signal
  
  # Mortality detection (scaled from BT: ~2.18x higher noise floor)
  mortality_cessation_threshold = 0.15,    # 15% of baseline movement
  mortality_consecutive_days    = 2,       # Consecutive low-movement days required
  
  create_distance_plots         = TRUE     # Generate time series plots
)

message("\n=== ENCOUNTER DETECTION STRATEGY (COW PARADISE) ===")
message(sprintf("GPS error: %.3fm | Combined pairwise: %.3fm",
                params$gps_error_sd, sqrt(2) * params$gps_error_sd))
message(sprintf("High-confidence: ≤%.1fm | Probable: %.1f-%.1fm | Possible: %.1f-%.1fm",
                params$high_confidence_threshold,
                params$high_confidence_threshold, params$probable_threshold,
                params$probable_threshold, params$possible_threshold))
message(sprintf("(BT equivalents were: ≤1.4m | 1.4-2.8m | 2.8-4.2m)"))
gps_combined <- sqrt(2) * params$gps_error_sd
gps_noise_daily <- gps_combined * sqrt(params$fixes_per_day)
message(sprintf("Expected GPS noise: ~%.0fm/day | Mortality threshold: ~%.0fm/day",
                gps_noise_daily, gps_noise_daily * 2.5))
message("====================================================\n")


# =================================================================-
# 2. FUNCTIONS ####
# =================================================================-

## GPS Noise Threshold ----
# Calculates apparent movement from GPS error alone (random walk model).
# A stationary fish will appear to move ~sqrt(n) × combined_error per day.
calculate_gps_noise_threshold <- function(gps_error = 0.827, fixes_per_day = 2880) {
  distance_error_sd            <- sqrt(2) * gps_error
  movement_per_fix             <- distance_error_sd
  daily_movement_random_walk   <- movement_per_fix * sqrt(fixes_per_day)
  daily_movement_realistic_max <- daily_movement_random_walk * 2
  data.frame(
    gps_error_sd                  = gps_error,
    fixes_per_day                 = fixes_per_day,
    distance_error_sd             = distance_error_sd,
    daily_movement_random_walk    = daily_movement_random_walk,
    daily_movement_realistic_max  = daily_movement_realistic_max
  )
}

## Pairwise Distance Calculation ----
# Uses ctmm::distances() to calculate prey-predator separations with GPS error
# accounting. Identical logic to BT script; thresholds passed as arguments.
calculate_pairwise_distances <- function(prey_tel, pred_tel, prey_fits, pred_fits,
                                         prey_species = "Prey", strike_dist = 0.45,
                                         high_conf_threshold = 3.0,
                                         probable_threshold   = 6.0,
                                         possible_threshold   = 9.0,
                                         parallel = TRUE, n_cores = NULL) {
  
  ctmm::projection(prey_tel)   <- ctmm::projection(pred_tel)
  ctmm::projection(prey_fits)  <- ctmm::projection(pred_fits)
  stopifnot(projection(prey_tel) == projection(pred_tel))
  
  n_prey       <- length(prey_tel)
  n_pred       <- length(pred_tel)
  combinations <- expand.grid(prey_idx = 1:n_prey, pred_idx = 1:n_pred)
  
  message(sprintf("Calculating distances for %d %s × %d Pike = %d combinations",
                  n_prey, prey_species, n_pred, nrow(combinations)))
  
  calc_dist <- function(idx) {
    i <- combinations$prey_idx[idx]
    j <- combinations$pred_idx[idx]
    tryCatch({
      combined_telemetry    <- c(prey_tel[i], pred_tel[j])
      combined_ctmm         <- c(prey_fits[i], pred_fits[j])
      location_difference   <- ctmm::distances(combined_telemetry, combined_ctmm)
      
      df           <- as.data.frame(location_difference)
      df$Prey_ID   <- names(prey_tel)[i]
      df$Pred_ID   <- names(pred_tel)[j]
      df$Species   <- prey_species
      
      # Legacy strike flag
      df$encounter <- ifelse(df$est <= strike_dist, 1, 0)
      
      # Tiered encounter flags (Cow Paradise thresholds)
      df$encounter_high_conf <- ifelse(df$est <= high_conf_threshold, 1, 0)
      df$encounter_probable  <- ifelse(df$est > high_conf_threshold &
                                         df$est <= probable_threshold, 1, 0)
      df$encounter_possible  <- ifelse(df$est > probable_threshold &
                                         df$est <= possible_threshold, 1, 0)
      
      df$encounter_type <- dplyr::case_when(
        df$est <= high_conf_threshold ~ "high_confidence",
        df$est <= probable_threshold  ~ "probable",
        df$est <= possible_threshold  ~ "possible",
        TRUE                          ~ "no_encounter"
      )
      
      if ("timestamp" %in% colnames(df)) {
        df$timestamp <- as.POSIXct(df$timestamp, tz = "Europe/Stockholm")
        df$Date      <- as.Date(df$timestamp, tz = "Europe/Stockholm")
      }
      return(df)
    }, error = function(e) NULL)
  }
  
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl,
                            c("prey_tel", "pred_tel", "prey_fits", "pred_fits",
                              "strike_dist", "prey_species", "high_conf_threshold",
                              "probable_threshold", "possible_threshold", "combinations"),
                            envir = environment())
    parallel::clusterEvalQ(cl, { library(ctmm); library(dplyr) })
    results <- parallel::parLapply(cl, 1:nrow(combinations), calc_dist)
  } else {
    results <- lapply(1:nrow(combinations), calc_dist)
  }
  
  results <- results[!sapply(results, is.null)]
  data.table::rbindlist(results) %>% as.data.frame()
}

## Daily Encounter Summaries ----
summarize_daily_encounters <- function(distances_df) {
  distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      Species                       = first(Species),
      encounter_count_high_conf     = sum(encounter_high_conf, na.rm = TRUE),
      encounter_count_probable      = sum(encounter_probable,  na.rm = TRUE),
      encounter_count_possible      = sum(encounter_possible,  na.rm = TRUE),
      daily_avg_dist                = mean(est, na.rm = TRUE),
      daily_min_dist                = min(est,  na.rm = TRUE),
      .groups = "drop"
    )
}

## Total Encounters Per Pair ----
calculate_total_encounters <- function(distances_df) {
  distances_df %>%
    group_by(Prey_ID, Pred_ID, Species) %>%
    summarise(
      total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
      total_encounters_probable  = sum(encounter_probable,  na.rm = TRUE),
      total_encounters_possible  = sum(encounter_possible,  na.rm = TRUE),
      min_distance               = min(est, na.rm = TRUE),
      mean_distance              = mean(est, na.rm = TRUE),
      n_days                     = n_distinct(Date),
      n_days_high_conf           = n_distinct(Date[encounter_high_conf == 1]),
      n_days_probable            = n_distinct(Date[encounter_probable  == 1]),
      .groups = "drop"
    ) %>%
    arrange(desc(total_encounters_high_conf))
}

## Identify Predation Events ----
# Ports BT's sophisticated multi-criteria detection, with all distance and
# encounter thresholds rescaled for Cow Paradise's higher GPS error (~2.18x BT).
#
# SCALING RATIONALE:
#   BT combined error = 0.536m → Cow Paradise combined error = 1.170m
#   Scale factor ≈ 2.18
#
#   Sustained-contact distance cut-offs:
#     BT:  <0.5m, <0.1m, <1.0m, <1.5m, <2.0m
#     Cow: <1.1m, <0.2m, <2.2m, <3.3m, <4.4m  (each × 2.18, rounded)
#
#   Encounter count thresholds remain the same (counts are scale-independent).
#
identify_predation_events <- function(distances_df,
                                      enc_threshold_low         = 10,
                                      enc_threshold_high        = 50,
                                      min_dist_threshold        = 3.0,
                                      consec_days_threshold     = 2,
                                      probable_threshold_low    = 50,
                                      probable_threshold_high   = 100,
                                      final_days_window         = 2) {
  
  # --- Daily summaries ---
  daily_summary <- distances_df %>%
    group_by(Prey_ID, Pred_ID, Date) %>%
    summarise(
      n_encounters_high_conf = sum(encounter_type == "high_confidence", na.rm = TRUE),
      n_encounters_probable  = sum(encounter_type == "probable",        na.rm = TRUE),
      n_encounters_possible  = sum(encounter_type == "possible",        na.rm = TRUE),
      min_distance           = min(est, na.rm = TRUE),
      mean_distance          = mean(est, na.rm = TRUE),
      # "at zero" = within combined GPS error (1.17m for Cow Paradise)
      n_at_zero              = sum(est < 1.17, na.rm = TRUE),
      pct_at_zero            = n_at_zero / n() * 100,
      .groups = "drop"
    ) %>%
    arrange(Prey_ID, Pred_ID, Date)
  
  # --- Sustained contact flagging (thresholds scaled × 2.18 from BT) ---
  predation_windows <- daily_summary %>%
    group_by(Prey_ID, Pred_ID) %>%
    arrange(Date) %>%
    mutate(
      sustained_contact = (
        # Very close sustained contact (BT <0.5m → Cow <1.1m; BT 10 enc → same)
        (min_distance < 1.1  & n_encounters_high_conf > 10) |
          # Near-zero contact  (BT <0.1m → Cow <0.2m)
          (min_distance < 0.2  & n_encounters_high_conf > 5) |
          # High pct at combined-error distance (threshold unchanged: 40%)
          (pct_at_zero > 40    & n_encounters_high_conf > 5) |
          
          # Moderate distance but high encounter density
          # (BT <1.0m/50enc → Cow <2.2m/50enc; BT <1.5m/100enc → Cow <3.3m/100enc)
          (min_distance < 2.2  & n_encounters_high_conf > 50)  |
          (min_distance < 3.3  & n_encounters_high_conf > 100) |
          
          # Mix of high-conf and probable
          # (BT <1.5m/20hc/50prob → Cow <3.3m/20hc/50prob)
          (min_distance < 3.3  & n_encounters_high_conf > 20 &
             n_encounters_probable > 50) |
          
          # Very high probable encounters alone
          # (BT <2.0m/150prob → Cow <4.4m/150prob)
          (min_distance < 4.4  & n_encounters_probable > 150)
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
  
  # --- Consecutive encounter days ---
  consecutive_info <- predation_windows %>%
    group_by(Prey_ID, Pred_ID) %>%
    arrange(Date) %>%
    mutate(
      high_conf_flag       = n_encounters_high_conf >= enc_threshold_low,
      is_consecutive       = (day_diff == 1 & high_conf_flag & lag(high_conf_flag)),
      consec_group         = cumsum(!is_consecutive | is.na(is_consecutive)),
      probable_flag        = n_encounters_probable >= probable_threshold_low,
      is_consecutive_prob  = (day_diff == 1 & probable_flag & lag(probable_flag)),
      consec_group_prob    = cumsum(!is_consecutive_prob | is.na(is_consecutive_prob))
    ) %>%
    group_by(Prey_ID, Pred_ID, consec_group) %>%
    mutate(consecutive_days_high_conf = if_else(high_conf_flag, n(), 0L)) %>%
    group_by(Prey_ID, Pred_ID, consec_group_prob) %>%
    mutate(consecutive_days_probable  = if_else(probable_flag, n(), 0L)) %>%
    ungroup()
  
  # --- Aggregate by pair ---
  consecutive_info %>%
    group_by(Prey_ID, Pred_ID) %>%
    summarise(
      total_encounters_high_conf        = sum(n_encounters_high_conf, na.rm = TRUE),
      num_days_high_conf_10             = sum(n_encounters_high_conf >= enc_threshold_low,  na.rm = TRUE),
      num_days_high_conf_50             = sum(n_encounters_high_conf >= enc_threshold_high, na.rm = TRUE),
      num_days_high_conf_100            = sum(n_encounters_high_conf >= 100, na.rm = TRUE),
      max_consecutive_days_high_conf    = { v <- consecutive_days_high_conf[!is.na(consecutive_days_high_conf)]; if (length(v)) max(v) else NA_integer_ },
      total_encounters_probable         = sum(n_encounters_probable,  na.rm = TRUE),
      num_days_probable_50              = sum(n_encounters_probable >= probable_threshold_low,  na.rm = TRUE),
      num_days_probable_100             = sum(n_encounters_probable >= probable_threshold_high, na.rm = TRUE),
      num_days_probable_200             = sum(n_encounters_probable >= 200, na.rm = TRUE),
      max_consecutive_days_probable     = { v <- consecutive_days_probable[!is.na(consecutive_days_probable)]; if (length(v)) max(v) else NA_integer_ },
      total_encounters_possible         = sum(n_encounters_possible, na.rm = TRUE),
      num_days_sustained_contact        = sum(sustained_contact, na.rm = TRUE),
      max_consecutive_contact_days      = { v <- consecutive_contact_days[!is.na(consecutive_contact_days)]; if (length(v)) max(v) else NA_integer_ },
      sustained_contact_start           = { v <- contact_period_start[!is.na(contact_period_start)]; if (length(v)) min(v) else as.Date(NA) },
      sustained_contact_end             = { v <- contact_period_end[!is.na(contact_period_end)];   if (length(v)) max(v) else as.Date(NA) },
      # "at zero" = within combined GPS error for Cow Paradise (1.17m)
      total_points_at_zero              = sum(n_at_zero, na.rm = TRUE),
      overall_min_distance              = min(min_distance,  na.rm = TRUE),
      mean_of_daily_min_distance        = mean(min_distance, na.rm = TRUE),
      # Distance thresholds scaled × 2.18 from BT (0.5m → 1.1m; 1.0m → 2.2m; 1.5m → 3.3m)
      num_days_under_1.1m               = sum(min_distance < 1.1, na.rm = TRUE),
      num_days_under_2.2m               = sum(min_distance < 2.2, na.rm = TRUE),
      num_days_under_3.3m               = sum(min_distance < 3.3, na.rm = TRUE),
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
        (total_encounters_probable  * 0.5) +
        (total_encounters_possible  * 0.2),
      
      # Predation likelihood classification.
      # Logic mirrors BT exactly; distance comparisons use Cow Paradise values.
      # "at zero" threshold uses 1.17m (combined GPS error for Cow Paradise).
      likely_predated = case_when(
        (num_days_sustained_contact >= 3 & max_consecutive_contact_days >= 2)  ~ "very_likely",
        (max_consecutive_contact_days >= 5)                                     ~ "very_likely",
        (total_points_at_zero > 100 & num_days_sustained_contact >= 2)         ~ "very_likely",
        (num_days_sustained_contact >= 2 & overall_min_distance < 1.1)         ~ "likely",
        (num_days_high_conf_50 > 1      & overall_min_distance < 1.1)          ~ "likely",
        (num_days_probable_100 > 1      & overall_min_distance < 1.65)         ~ "likely",  # BT 0.75m × 2.18
        num_days_sustained_contact >= 2                                         ~ "likely",
        (num_days_high_conf_10 > 2      & total_encounters_probable > 250)      ~ "possible",
        (num_days_probable_50 > 1       & overall_min_distance < 2.2)          ~ "possible",  # BT 1.0m × 2.18
        (total_encounters_high_conf > 200 & num_days_under_2.2m > 0)           ~ "possible",
        num_days_sustained_contact >= 1                                         ~ "possible",
        (total_encounters_high_conf > 150 | total_encounters_probable > 400)   ~ "weak_signal",
        TRUE                                                                    ~ "unlikely"
      )
    )
}

## Distance Time Series Plot ----
# Creates plots showing prey-predator separation over time.
# Threshold lines drawn at Cow Paradise-specific values (3.0m / 6.0m).
plot_distance_timeseries <- function(distances_df, prey_id, pred_id, species,
                                     output_dir = NULL, threshold_lines = TRUE,
                                     predation_summary = NULL) {
  
  pair_data <- distances_df %>% filter(Prey_ID == prey_id, Pred_ID == pred_id) %>%
    arrange(timestamp)
  if (nrow(pair_data) == 0) return(NULL)
  
  predation_info <- NULL
  if (!is.null(predation_summary)) {
    col_check <- if ("Pike_ID" %in% colnames(predation_summary)) "Pike_ID" else "Pred_ID"
    predation_info <- predation_summary %>%
      filter(Prey_ID == prey_id, .data[[col_check]] == pred_id)
  }
  
  p <- ggplot(pair_data, aes(x = timestamp, y = est)) +
    geom_line(color = "gray50", alpha = 0.5) +
    geom_point(aes(color = encounter_type), size = 1, alpha = 0.6) +
    scale_color_manual(
      values = c("high_confidence" = "#d62728", "probable" = "#ff7f0e",
                 "possible" = "#1f77b4", "no_encounter" = "gray80"),
      labels = c("high_confidence" = "High Confidence (≤3.0m)",
                 "probable"        = "Probable (3.0–6.0m)",
                 "possible"        = "Possible (6.0–9.0m)",
                 "no_encounter"    = "No Encounter"),
      name = "Encounter Type"
    ) +
    labs(
      title    = sprintf("%s: %s vs Pike %s", species, prey_id, pred_id),
      subtitle = if (!is.null(predation_info) && nrow(predation_info) > 0) {
        sprintf("Predation likelihood: %s | Sustained contact days: %d",
                predation_info$likely_predated, predation_info$num_days_sustained_contact)
      } else NULL,
      x = "Date", y = "Distance (m)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, color = "gray30"))
  
  if (threshold_lines) {
    p <- p +
      # Cow Paradise thresholds: 1.1m (near-zero), 3.0m (high-conf), 6.0m (probable)
      geom_hline(yintercept = 1.1, linetype = "dashed", color = "#8B0000", alpha = 0.6) +
      geom_hline(yintercept = 3.0, linetype = "dashed", color = "#d62728", alpha = 0.5) +
      geom_hline(yintercept = 6.0, linetype = "dashed", color = "#ff7f0e", alpha = 0.5) +
      geom_hline(yintercept = 9.0, linetype = "dashed", color = "#1f77b4", alpha = 0.5) +
      annotate("text", x = min(pair_data$timestamp), y = 1.1,
               label = "1.1m (≈0 contact)", hjust = 0, vjust = -0.5, size = 3, color = "#8B0000") +
      annotate("text", x = min(pair_data$timestamp), y = 3.0,
               label = "3.0m (high-conf)",  hjust = 0, vjust = -0.5, size = 3, color = "#d62728") +
      annotate("text", x = min(pair_data$timestamp), y = 6.0,
               label = "6.0m (probable)",   hjust = 0, vjust = -0.5, size = 3, color = "#ff7f0e") +
      annotate("text", x = min(pair_data$timestamp), y = 9.0,
               label = "9.0m (possible)",   hjust = 0, vjust = -0.5, size = 3, color = "#1f77b4")
  }
  
  # Highlight sustained contact period
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
               ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) +
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
  }
  
  return(p)
}

## Movement Metrics for Mortality Detection ----
# GPS noise for Cow Paradise (~80m/day random walk, ~200m/day threshold)
# is ~2.18x higher than BT (~37m/day, ~92m/day). Same proportional logic applies.
plot_movement_metrics <- function(telemetry_obj, individual_id, species, output_dir,
                                  gps_error = 0.827, cessation_threshold = 0.15,
                                  consecutive_days = 2) {
  
  tryCatch({
    if (!inherits(telemetry_obj, "telemetry")) stop("Input is not a telemetry object")
    
    tel_df <- data.frame(
      timestamp = telemetry_obj$timestamp,
      x = telemetry_obj$x,
      y = telemetry_obj$y,
      stringsAsFactors = FALSE
    )
    if (nrow(tel_df) < 2) return(NULL)
    
    tel_df$timestamp <- as.POSIXct(tel_df$timestamp, origin = "1970-01-01", tz = "Europe/Stockholm")
    tel_df$Date      <- as.Date(tel_df$timestamp, tz = "Europe/Stockholm")
    
    tel_df <- tel_df %>%
      arrange(timestamp) %>%
      mutate(
        step_distance = sqrt((x - lag(x))^2 + (y - lag(y))^2),
        time_interval = as.numeric(difftime(timestamp, lag(timestamp), units = "hours")),
        speed         = step_distance / time_interval
      )
    
    daily_metrics <- tel_df %>%
      group_by(Date) %>%
      summarise(
        total_distance = sum(step_distance, na.rm = TRUE),
        mean_speed     = mean(speed,         na.rm = TRUE),
        n_observations = n(),
        .groups = "drop"
      ) %>%
      filter(!is.na(Date)) %>%
      arrange(Date)
    
    if (nrow(daily_metrics) == 0) return(NULL)
    
    # GPS noise (Cow Paradise: ~80m/day random walk; threshold = 2.5× = ~200m/day)
    mean_fixes_per_day <- mean(daily_metrics$n_observations, na.rm = TRUE)
    gps_noise_estimate <- calculate_gps_noise_threshold(
      gps_error = gps_error, fixes_per_day = mean_fixes_per_day)
    gps_noise_threshold <- gps_noise_estimate$daily_movement_random_walk * 2.5
    
    # Baseline: first 50% of tracking or first 7 days
    n_baseline_days  <- max(7, floor(nrow(daily_metrics) * 0.5))
    baseline_days    <- min(n_baseline_days, nrow(daily_metrics))
    baseline_distance <- daily_metrics %>% slice_head(n = baseline_days) %>%
      pull(total_distance) %>% mean(na.rm = TRUE)
    
    # Mortality threshold = MAX(relative, absolute GPS noise)
    cessation_distance_threshold <- max(baseline_distance * cessation_threshold,
                                        gps_noise_threshold)
    
    daily_metrics <- daily_metrics %>%
      mutate(low_movement   = total_distance < cessation_distance_threshold,
             consecutive_low = NA)
    
    # Identify consecutive low-movement periods
    mortality_date <- NULL
    mortality_info <- NULL
    
    if (any(daily_metrics$low_movement, na.rm = TRUE)) {
      rle_result  <- rle(daily_metrics$low_movement)
      run_ends    <- cumsum(rle_result$lengths)
      run_starts  <- c(1, run_ends[-length(run_ends)] + 1)
      long_low    <- which(rle_result$values == TRUE & rle_result$lengths >= consecutive_days)
      
      if (length(long_low) > 0) {
        for (period_idx in long_low) {
          s <- run_starts[period_idx]; e <- run_ends[period_idx]
          daily_metrics$consecutive_low[s:e] <- TRUE
        }
      }
    }
    
    if (any(daily_metrics$consecutive_low == TRUE, na.rm = TRUE)) {
      mortality_idx  <- which(daily_metrics$consecutive_low == TRUE)[1]
      mortality_date <- daily_metrics$Date[mortality_idx]
      post_mortality <- daily_metrics %>% filter(Date >= mortality_date)
      mortality_info <- list(
        date                  = mortality_date,
        mean_distance_after   = mean(post_mortality$total_distance, na.rm = TRUE),
        baseline_distance     = baseline_distance,
        movement_vs_gps_noise = mean(post_mortality$total_distance, na.rm = TRUE) / gps_noise_threshold
      )
    }
    
    date_breaks <- seq(from = min(daily_metrics$Date), to = max(daily_metrics$Date), by = "2 days")
    
    p <- ggplot(daily_metrics, aes(x = Date, y = total_distance)) +
      geom_line(color = "#264653", linewidth = 0.8, alpha = 0.8) +
      geom_point(color = "#264653", size = 2, alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE, color = "#e76f51", fill = "#e76f51",
                  alpha = 0.2, linewidth = 0.6) +
      geom_hline(yintercept = gps_noise_threshold, linetype = "dotted",
                 color = "gray40", alpha = 0.7, linewidth = 0.8) +
      annotate("text", x = min(daily_metrics$Date), y = gps_noise_threshold,
               label = sprintf("GPS noise (~%.0fm)", gps_noise_threshold),
               hjust = 0, vjust = -0.5, size = 2.5, color = "gray40")
    
    if (!is.null(mortality_date)) {
      p <- p +
        geom_vline(xintercept = as.numeric(mortality_date), linetype = "dashed",
                   color = "#d62828", linewidth = 1.2) +
        annotate("text", x = mortality_date,
                 y = max(daily_metrics$total_distance, na.rm = TRUE) * 0.95,
                 label = sprintf("Potential mortality\n%s\n(%.1fx GPS noise)",
                                 format(mortality_date, "%d %b"),
                                 mortality_info$movement_vs_gps_noise),
                 hjust = -0.05, vjust = 1, size = 3.5, color = "#d62828", fontface = "bold")
    }
    
    p <- p +
      scale_x_date(breaks = date_breaks, date_labels = "%d %b", expand = c(0.02, 0)) +
      labs(
        title    = sprintf("%s %s: Daily Total Distance Moved", species, individual_id),
        subtitle = sprintf("Mean: %.1fm/day | Baseline: %.1fm/day | GPS noise: ~%.0fm | (σ_GPS = %.3fm)",
                           mean(daily_metrics$total_distance, na.rm = TRUE),
                           baseline_distance, gps_noise_threshold, gps_error),
        x = "Date", y = "Total distance moved (m)"
      ) +
      theme_classic() +
      theme(plot.title    = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 8, color = "gray30"),
            axis.text     = element_text(size = 9),
            axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1))
    
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      filename <- sprintf("%s_%s_movement_metrics.png", species, individual_id)
      ggsave(file.path(output_dir, filename), plot = p, width = 10, height = 6, dpi = 300)
    }
    
    return(list(plot = p, metrics = daily_metrics, mortality_info = mortality_info,
                baseline_distance = baseline_distance, gps_noise_threshold = gps_noise_threshold))
    
  }, error = function(e) {
    warning(sprintf("Failed to create movement plots for %s: %s", individual_id, e$message))
    return(NULL)
  })
}


# =================================================================-
# 3. DATA LOADING ####
# =================================================================-

pike_tel   <- readRDS(file.path(paths$telem, 'pike_cow_paradise_tel_thinned.rds'))
perch_tel  <- readRDS(file.path(paths$telem, 'perch_cow_paradise_tel_thinned.rds'))
roach_tel  <- readRDS(file.path(paths$telem, 'roach_cow_paradise_tel_thinned.rds'))

pike_fits  <- readRDS(file.path(paths$ctmm, "cow_paradise_pike_fits/cow_paradise_pike_best_models.rds"))
perch_fits <- readRDS(file.path(paths$ctmm, "cow_paradise_perch_fits/cow_paradise_perch_best_models.rds"))
roach_fits <- readRDS(file.path(paths$ctmm, "cow_paradise_roach_fits/cow_paradise_roach_best_models.rds"))

cow_filt_data  <- readRDS(file.path(paths$filtered_data, "04_cow_sub.rds"))
post_biometrics <- fread(file.path(paths$size, "biometric_post_exp_data.csv")) %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number)))

post_biometric_cols <- post_biometrics %>%
  filter(Lake == 'Cow Paradise', Species %in% c('Roach', 'Perch')) %>%
  select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found, found_predated = Known_predated)

message(sprintf("Loaded: %d pike, %d perch, %d roach",
                length(pike_tel), length(perch_tel), length(roach_tel)))


# =================================================================-
# 4. DISTANCE CALCULATIONS ####
# =================================================================-

message("\n=== CALCULATING DISTANCES ===")

## Roach-Pike ----
roach_pike_file <- file.path(paths$encounters, "cow_pike_roach_distances_tiered_df.rds")
if (file.exists(roach_pike_file)) {
  roach_pike_distances_df <- readRDS(roach_pike_file)
} else {
  roach_pike_distances_df <- calculate_pairwise_distances(
    prey_tel = roach_tel, pred_tel = pike_tel,
    prey_fits = roach_fits, pred_fits = pike_fits,
    prey_species = "Roach", strike_dist = params$strike_distance,
    high_conf_threshold = params$high_confidence_threshold,
    probable_threshold  = params$probable_threshold,
    possible_threshold  = params$possible_threshold,
    parallel = TRUE, n_cores = 16
  )
  saveRDS(roach_pike_distances_df, roach_pike_file)
}

## Perch-Pike ----
perch_pike_file <- file.path(paths$encounters, "cow_pike_perch_distances_tiered_df.rds")
if (file.exists(perch_pike_file)) {
  perch_pike_distances_df <- readRDS(perch_pike_file)
} else {
  perch_pike_distances_df <- calculate_pairwise_distances(
    prey_tel = perch_tel, pred_tel = pike_tel,
    prey_fits = perch_fits, pred_fits = pike_fits,
    prey_species = "Perch", strike_dist = params$strike_distance,
    high_conf_threshold = params$high_confidence_threshold,
    probable_threshold  = params$probable_threshold,
    possible_threshold  = params$possible_threshold,
    parallel = TRUE, n_cores = 16
  )
  saveRDS(perch_pike_distances_df, perch_pike_file)
}

message(sprintf("Distances: %d roach-pike, %d perch-pike observations",
                nrow(roach_pike_distances_df), nrow(perch_pike_distances_df)))


# =================================================================-
# 5. ENCOUNTER SUMMARIES ####
# =================================================================-

message("\n=== CALCULATING ENCOUNTER SUMMARIES ===")

## Daily summaries ----
roach_daily <- summarize_daily_encounters(roach_pike_distances_df) %>% rename(Pike_ID = Pred_ID)
perch_daily <- summarize_daily_encounters(perch_pike_distances_df) %>% rename(Pike_ID = Pred_ID)
prey_daily_encounters <- bind_rows(roach_daily, perch_daily)

## Total encounters ----
roach_total <- calculate_total_encounters(roach_pike_distances_df)
perch_total <- calculate_total_encounters(perch_pike_distances_df)
all_encounters <- bind_rows(roach_total, perch_total)

## Encounters by prey ID ----
total_encounters_ID <- all_encounters %>%
  group_by(Prey_ID, Species) %>%
  summarise(
    total_encounters_high_conf = sum(total_encounters_high_conf),
    total_encounters_probable  = sum(total_encounters_probable),
    total_encounters_possible  = sum(total_encounters_possible),
    n_predators                = n_distinct(Pred_ID),
    .groups = 'drop'
  ) %>%
  left_join(cow_filt_data %>% select(individual_ID, treatment) %>% distinct(),
            by = c("Prey_ID" = "individual_ID")) %>%
  arrange(desc(total_encounters_high_conf))

saveRDS(total_encounters_ID, file.path(paths$encounters, "cow_total_encounters_ID.rds"))


# =================================================================-
# 6. IDENTIFY PREDATION EVENTS ####
# =================================================================-

message("\n=== IDENTIFYING PREDATION EVENTS ===")

roach_predation <- identify_predation_events(
  roach_pike_distances_df,
  enc_threshold_low      = params$encounter_threshold_low,
  enc_threshold_high     = params$encounter_threshold_high,
  min_dist_threshold     = params$min_distance_threshold,
  consec_days_threshold  = params$consecutive_days_threshold
) %>% rename(Pike_ID = Pred_ID) %>% mutate(Species = "Roach")

perch_predation <- identify_predation_events(
  perch_pike_distances_df,
  enc_threshold_low      = params$encounter_threshold_low,
  enc_threshold_high     = params$encounter_threshold_high,
  min_dist_threshold     = params$min_distance_threshold,
  consec_days_threshold  = params$consecutive_days_threshold
) %>% rename(Pike_ID = Pred_ID) %>% mutate(Species = "Perch")

## Filter for suspected predation ----
# Thresholds for filtering are scaled for Cow Paradise:
#   BT overall_min_distance < 0.75m → Cow < 1.65m  (× 2.18)
#   BT num_days_under_1m > 2        → Cow num_days_under_2.2m > 2
#   BT encounters_final_2days_high_conf > 100 → same (count-based)
suspected_predation_distance <- bind_rows(roach_predation, perch_predation) %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
  filter(
    found_alive == 0,
    (num_days_high_conf_50 > 1 | num_days_probable_50 > 1 | num_days_probable_100 > 0 |
       total_encounters_high_conf > 200 | total_encounters_probable > 500 |
       combined_encounter_score > 400   | overall_min_distance < 1.65 |
       num_days_under_2.2m > 2          | encounters_final_2days_high_conf > 100 |
       encounters_final_2days_probable > 200 |
       (num_days_high_conf_10 > 2 & total_encounters_probable > 250) |
       likely_predated %in% c("very_likely", "likely", "possible"))
  ) %>%
  arrange(desc(likely_predated), desc(combined_encounter_score))

message(sprintf("Suspected predation events: %d", nrow(suspected_predation_distance)))

## Summary ----
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
print(predation_summary_table)


# =================================================================-
# 7. DISTANCE PLOTS ####
# =================================================================-

if (params$create_distance_plots) {
  message("\n=== CREATING DISTANCE PLOTS ===")
  
  high_encounter_pairs <- suspected_predation_distance %>%
    filter(num_days_high_conf_10 >= 2 | num_days_high_conf_50 >= 1 |
             num_days_sustained_contact >= 1)
  
  message(sprintf("Creating %d distance plots...", nrow(high_encounter_pairs)))
  
  for (i in 1:nrow(high_encounter_pairs)) {
    prey_id  <- high_encounter_pairs$Prey_ID[i]
    pred_id  <- high_encounter_pairs$Pike_ID[i]
    species  <- high_encounter_pairs$Species[i]
    
    dist_data <- if (species == "Roach") roach_pike_distances_df else perch_pike_distances_df
    pred_data <- if (species == "Roach") roach_predation         else perch_predation
    
    plot_distance_timeseries(
      distances_df = dist_data, prey_id = prey_id, pred_id = pred_id, species = species,
      output_dir = paths$distance_plots, threshold_lines = TRUE, predation_summary = pred_data
    )
  }
  
  message(sprintf("Created %d distance plots", nrow(high_encounter_pairs)))
}


# =================================================================-
# 8. MOVEMENT METRICS FOR MORTALITY DETECTION ####
# =================================================================-

message("\n=== CREATING MOVEMENT METRICS PLOTS ===")

movement_plots_dir <- file.path(paths$figures, "movement_metrics")
if (!dir.exists(movement_plots_dir)) dir.create(movement_plots_dir, recursive = TRUE)

# Print GPS context for Cow Paradise
gps_ref <- calculate_gps_noise_threshold(params$gps_error_sd, params$fixes_per_day)
message(sprintf("GPS noise context (Cow Paradise, σ = %.3fm):", params$gps_error_sd))
message(sprintf("  Random walk apparent movement:  ~%.0fm/day", gps_ref$daily_movement_random_walk))
message(sprintf("  Mortality detection threshold:  ~%.0fm/day (2.5× random walk)",
                gps_ref$daily_movement_random_walk * 2.5))
message(sprintf("  (BT equivalents: ~%.0fm/day and ~%.0fm/day)",
                calculate_gps_noise_threshold(0.379, 2880)$daily_movement_random_walk,
                calculate_gps_noise_threshold(0.379, 2880)$daily_movement_random_walk * 2.5))

missing_individuals <- post_biometric_cols %>% filter(found_alive == 0)
message(sprintf("Creating movement plots for %d missing individuals", nrow(missing_individuals)))

movement_plot_results <- list()

## Process Roach ----
missing_roach <- missing_individuals %>% filter(individual_ID %in% names(roach_tel))
if (nrow(missing_roach) > 0) {
  message(sprintf("Processing %d missing Roach...", nrow(missing_roach)))
  for (i in 1:nrow(missing_roach)) {
    ind_id <- missing_roach$individual_ID[i]
    if (ind_id %in% names(roach_tel)) {
      movement_plot_results[[ind_id]] <- plot_movement_metrics(
        telemetry_obj      = roach_tel[[ind_id]], individual_id = ind_id, species = "Roach",
        output_dir         = movement_plots_dir,  gps_error = params$gps_error_sd,
        cessation_threshold = params$mortality_cessation_threshold,
        consecutive_days    = params$mortality_consecutive_days
      )
    }
  }
}

## Process Perch ----
missing_perch <- missing_individuals %>% filter(individual_ID %in% names(perch_tel))
if (nrow(missing_perch) > 0) {
  message(sprintf("Processing %d missing Perch...", nrow(missing_perch)))
  for (i in 1:nrow(missing_perch)) {
    ind_id <- missing_perch$individual_ID[i]
    if (ind_id %in% names(perch_tel)) {
      movement_plot_results[[ind_id]] <- plot_movement_metrics(
        telemetry_obj      = perch_tel[[ind_id]], individual_id = ind_id, species = "Perch",
        output_dir         = movement_plots_dir,  gps_error = params$gps_error_sd,
        cessation_threshold = params$mortality_cessation_threshold,
        consecutive_days    = params$mortality_consecutive_days
      )
    }
  }
}

## Compile mortality summary ----
mortality_summary <- bind_rows(
  lapply(names(movement_plot_results), function(ind_id) {
    result  <- movement_plot_results[[ind_id]]
    species <- ifelse(ind_id %in% names(roach_tel), "Roach", "Perch")
    if (!is.null(result)) {
      if (!is.null(result$mortality_info)) {
        data.frame(
          individual_ID         = ind_id, species = species,
          potential_mortality_date = result$mortality_info$date,
          baseline_distance     = result$baseline_distance,
          mean_distance_after   = result$mortality_info$mean_distance_after,
          gps_noise_threshold   = result$gps_noise_threshold,
          movement_vs_gps_noise = result$mortality_info$movement_vs_gps_noise,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          individual_ID         = ind_id, species = species,
          potential_mortality_date = as.Date(NA),
          baseline_distance     = result$baseline_distance,
          mean_distance_after   = NA,
          gps_noise_threshold   = result$gps_noise_threshold,
          movement_vs_gps_noise = NA,
          stringsAsFactors = FALSE
        )
      }
    }
  })
)

if (!is.null(mortality_summary) && nrow(mortality_summary) > 0) {
  write.csv(mortality_summary,
            file.path(paths$tables, "mortality_detection_summary.csv"), row.names = FALSE)
  message(sprintf("\nMortality signals detected: %d / %d individuals",
                  sum(!is.na(mortality_summary$potential_mortality_date)),
                  nrow(mortality_summary)))
}


# =================================================================-
# 9. VISUALISATIONS ####
# =================================================================-

message("\n=== CREATING VISUALIZATIONS ===")

combined_distances <- bind_rows(roach_pike_distances_df, perch_pike_distances_df)

## Encounter distribution ----
p1 <- combined_distances %>%
  filter(encounter_type %in% c("high_confidence", "probable")) %>%
  left_join(cow_filt_data %>% select(individual_ID, treatment) %>% distinct(),
            by = c("Prey_ID" = "individual_ID")) %>%
  filter(!is.na(treatment)) %>%
  mutate(encounter_type = factor(encounter_type,
                                 levels = c("high_confidence", "probable"),
                                 labels = c("High Confidence (≤3.0m)", "Probable (3.0–6.0m)"))) %>%
  ggplot(aes(x = encounter_type, fill = treatment)) +
  geom_bar(position = "dodge", alpha = 0.8, color = "black", linewidth = 0.5) +
  facet_wrap(~Species, scales = "free_y") +
  scale_fill_brewer(palette = "Set1", name = "Treatment") +
  labs(y = "Count",
       caption = sprintf("GPS error: σ = %.3fm | Thresholds scaled × 2.18 from BT lake",
                         params$gps_error_sd)) +
  theme_classic() +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x  = element_blank(),
        strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 1),
        strip.text    = element_text(face = "bold"),
        panel.border  = element_rect(color = "black", fill = NA, linewidth = 1))

print(p1)
ggsave(file.path(paths$figures, "encounter_type_distribution.png"),
       p1, width = 12, height = 6, dpi = 300)


# =================================================================-
# 10. EXPORT RESULTS ####
# =================================================================-

message("\n=== EXPORTING RESULTS ===")

wb <- createWorkbook()
addWorksheet(wb, "Daily Encounters");    writeData(wb, "Daily Encounters",    prey_daily_encounters)
addWorksheet(wb, "Total Encounters");    writeData(wb, "Total Encounters",    all_encounters)
addWorksheet(wb, "Suspected Predation"); writeData(wb, "Suspected Predation", suspected_predation_distance)
addWorksheet(wb, "Summary");             writeData(wb, "Summary",             predation_summary_table)
saveWorkbook(wb, file.path(paths$tables, "cow_predation_analysis.xlsx"), overwrite = TRUE)

write.csv(suspected_predation_distance,
          file.path(paths$tables, "suspected_predation_events.csv"),      row.names = FALSE)
write.csv(all_encounters,
          file.path(paths$tables, "all_encounter_totals.csv"),            row.names = FALSE)
write.csv(prey_daily_encounters,
          file.path(paths$tables, "daily_encounters.csv"),                row.names = FALSE)

message("\n=== ANALYSIS COMPLETE ===")
message(sprintf("Suspected predation events: %d", nrow(suspected_predation_distance)))
message(sprintf("Events with sustained contact: %d",
                sum(suspected_predation_distance$num_days_sustained_contact > 0)))
message(sprintf("Output: %s", paths$tables))