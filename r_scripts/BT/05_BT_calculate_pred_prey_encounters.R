# ===================================================================
# Predator-Prey Interaction Analysis - BT
# ===================================================================
# 
# ANALYSIS OVERVIEW:
# This script identifies potential predation events by analyzing movement 
# tracking data for prey (roach/perch) and predators (pike). Predation is 
# inferred from:
#   1. Close proximity (<0.45m - pike strike distance)
#   2. Sustained encounters over consecutive days
#   3. Movement pattern changes suggesting predation
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
  library(parallel)      # For parallel processing
  library(future.apply)  # For parallelized apply functions
  library(progressr)     # For progress bars
})

# Configure time zone
Sys.setenv(TZ = 'Europe/Stockholm')

## 1.2 Define Directory Paths ----
paths <- list(
  filtered_data = "./data/tracks_filtered/BT/",
  lake_polygon  = "./data/lake_coords/",
  ctmm          = "./data/ctmm_fits/",
  telem         = "./data/telem_obj/BT/",
  encounters    = "./data/encounters/BT/",
  size          = "./data/fish_size/",
  tables        = "./tables/BT/",
  figures       = "./figures/BT/" 
)

## 1.3 Define Analysis Parameters ----
params <- list(
  strike_distance = 0.45,      # Pike max strike distance (m)
  encounter_threshold_low = 25,  # Minimum encounters for suspected predation
  encounter_threshold_high = 100, # High encounter count threshold
  min_distance_threshold = 0.2,   # Minimum distance threshold (m)
  consecutive_days_threshold = 2  # Consecutive days for likely predation
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
#' 
#' @details This function performs the computationally intensive distance
#' calculation. Parallelization significantly reduces runtime for large datasets.
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
    
    # Set up parallel backend
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Export necessary objects to cluster
    parallel::clusterExport(cl, c("prey_tel", "pred_tel", "prey_fits", 
                                  "pred_fits", "strike_dist", "prey_species",
                                  "combinations"),
                            envir = environment())
    
    # Load required packages on each worker
    parallel::clusterEvalQ(cl, {
      library(ctmm)
      library(dplyr)
    })
    
    # Run parallel computation
    results <- parallel::parLapply(cl, 1:nrow(combinations), calc_dist)
  } else {
    # Sequential processing with progress bar
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
      # Count days exceeding thresholds
      num_days_25_encounters = sum(encounter_count >= enc_threshold_low),
      num_days_100_encounters = sum(encounter_count >= enc_threshold_high),
      num_days_min_dist_less_0.2m = sum(daily_min_dist < min_dist_threshold),
      
      # Identify key dates
      first_date_over_25 = if(any(encounter_count >= enc_threshold_low)) 
        min(Date[encounter_count >= enc_threshold_low]) else as.Date(NA),
      first_date_over_100 = if(any(encounter_count >= enc_threshold_high)) 
        min(Date[encounter_count >= enc_threshold_high]) else as.Date(NA),
      
      # Calculate consecutive days
      consecutive_days_25 = calc_max_consecutive(Date[encounter_count >= enc_threshold_low]),
      consecutive_days_100 = calc_max_consecutive(Date[encounter_count >= enc_threshold_high]),
      
      # Compile encounter dates
      encounter_dates = paste(format(Date, "%d/%m/%Y"), collapse = ", "),
      
      .groups = "drop"
    ) %>%
    # Flag likely predation events
    mutate(likely_predated = ifelse(consecutive_days_100 >= consec_days_threshold, 1, 0))
  
  return(predation_summary)
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


#' Filter Post-Predation Data
#' 
#' @param encounter_df Encounter summary data frame
#' @param mortality_data Data frame with mortality dates
#' @return Filtered data frame excluding post-predation encounters
filter_post_predation <- function(encounter_df, mortality_data) {
  
  encounter_filtered <- encounter_df %>%
    left_join(
      mortality_data %>% select(individual_ID, death_date = revised_likely_death_date),
      by = c("Prey_ID" = "individual_ID")
    ) %>%
    filter(is.na(death_date) | Date <= death_date)
  
  message(sprintf("Filtered from %d to %d rows (removed post-predation data)",
                  nrow(encounter_df), nrow(encounter_filtered)))
  
  return(encounter_filtered)
}


#' Create Encounter Summary Statistics
#' 
#' @param daily_encounters Daily encounter data
#' @param metadata Individual metadata (treatment, tracking quality)
#' @return Summary statistics by individual
create_encounter_summary <- function(daily_encounters, metadata) {
  
  summary_stats <- daily_encounters %>%
    group_by(Prey_ID) %>%
    summarise(
      Species = first(Species),
      total_encounters = sum(encounter_count, na.rm = TRUE),
      mean_daily_encounters = mean(encounter_count, na.rm = TRUE),
      max_daily_encounters = max(encounter_count, na.rm = TRUE),
      max_encounter_date = Date[which.max(encounter_count)],
      min_distance_overall = min(daily_min_dist, na.rm = TRUE),
      mean_distance = mean(daily_avg_dist, na.rm = TRUE),
      n_tracking_days = n_distinct(Date),
      .groups = "drop"
    ) %>%
    left_join(metadata, by = c("Prey_ID" = "individual_ID")) %>%
    mutate(
      days_tracked = 36 - n_missing_dates,
      avg_daily_encounter_rate = round(total_encounters / days_tracked, 2)
    )
  
  return(summary_stats)
}

# ===================================================================
# 3. DATA LOADING
# ===================================================================

## 3.1 Load Telemetry Objects ----
pike_tel <- readRDS(file.path(paths$telem, 'pike_BT_tel_thinned.rds'))
perch_tel <- readRDS(file.path(paths$telem, 'perch_BT_tel_thinned.rds'))
roach_tel <- readRDS(file.path(paths$telem, 'roach_BT_tel_thinned.rds'))

## 3.2 Load CTMM Fits ----
pike_fits <- readRDS(file.path(paths$ctmm, "BT_pike_fits/BT_pike_ctmm_fits.rds"))
perch_fits <- readRDS(file.path(paths$ctmm, "BT_perch_fits/BT_perch_ctmm_fits.rds"))
roach_fits <- readRDS(file.path(paths$ctmm, "BT_roach_fits/BT_roach_ctmm_fits.rds"))

## 3.3 Load Metadata ----
BT_filt_data <- readRDS(file.path(paths$filtered_data, "03_BT_sub.rds"))
post_biometrics <- fread(file.path(paths$size, "biometric_post_exp_data.csv")) %>%
  mutate(individual_ID = paste0("F", sub(".*-", "", Tag_Number)))

# ===================================================================
# 4. DISTANCE CALCULATIONS
# ===================================================================

## 4.1 Roach-Pike Distances ----
# Check if already calculated
roach_pike_file <- file.path(paths$encounters, "BT_pike_roach_distances_df.rds")

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
    n_cores = 20  # Auto-detect
  )
  saveRDS(roach_pike_distances_df, roach_pike_file)
}

## 4.2 Perch-Pike Distances ----
perch_pike_file <- file.path(paths$encounters, "BT_pike_perch_distances_df.rds")

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
    n_cores = 20
  )
  saveRDS(perch_pike_distances_df, perch_pike_file)
}


# ===================================================================
# 5. ENCOUNTER ANALYSIS
# ===================================================================

message("Analyzing encounters...")

## 5.1 Daily Encounter Summaries ----
roach_daily <- summarize_daily_encounters(roach_pike_distances_df) %>%
  rename(Pike_ID = Pred_ID)

perch_daily <- summarize_daily_encounters(perch_pike_distances_df) %>%
  rename(Pike_ID = Pred_ID)

# Combine prey species
prey_daily_encounters <- bind_rows(roach_daily, perch_daily)

## 5.2 Identify Suspected Predation Events ----
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

# Combine suspected predation events
suspected_predation <- bind_rows(roach_predation, perch_predation)

# Add biometric data
post_biometric_cols <- post_biometrics %>%
  filter(Lake == 'BT', Species %in% c('Roach', 'Perch')) %>%
  select(individual_ID, Found, Known_predated) %>%
  rename(found_alive = Found, found_predated = Known_predated)

suspected_predation <- suspected_predation %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID")) %>%
  filter(found_alive == 0)  # Focus on individuals not recovered alive

# Save results
saveRDS(suspected_predation, 
        file.path(paths$encounters, "BT_suspected_predation_events.rds"))

message(sprintf("Identified %d suspected predation events", nrow(suspected_predation)))


## 5.3 Total Encounter Summary ----
# Prepare metadata
prey_metadata <- BT_filt_data %>%
  filter(Species %in% c('Roach', 'Perch')) %>%
  select(individual_ID, Treatment, n_missing_dates) %>%
  distinct()

# Create comprehensive summary
prey_encounter_summary <- create_encounter_summary(prey_daily_encounters, prey_metadata)

# Add survival data
prey_encounter_summary <- prey_encounter_summary %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID"))

# Save
saveRDS(prey_encounter_summary,
        file.path(paths$encounters, "BT_prey_total_encounter_summary.rds"))


# ===================================================================
# 6. FILTER POST-PREDATION DATA
# ===================================================================

message("Filtering post-predation encounters...")

# Load mortality predictions
mortality_preds <- readxl::read_excel("./data/encounters/suspected_mortality_updated.xlsx") %>%
  filter(lake == 'BT', species %in% c('Roach', 'Perch'))

# Filter encounters
prey_encounters_filtered <- filter_post_predation(prey_daily_encounters, mortality_preds)

# Recalculate summary with filtered data
prey_summary_filtered <- create_encounter_summary(prey_encounters_filtered, prey_metadata) %>%
  left_join(post_biometric_cols, by = c("Prey_ID" = "individual_ID"))

# Save filtered results
saveRDS(prey_encounters_filtered,
        file.path(paths$encounters, "BT_prey_encounters_filtered.rds"))
saveRDS(prey_summary_filtered,
        file.path(paths$encounters, "BT_prey_encounter_summary_filtered.rds"))


# ===================================================================
# 7. STATISTICAL COMPARISON BY TREATMENT
# ===================================================================

message("Comparing encounter rates by treatment...")

## 7.1 Treatment-Level Summary ----
treatment_comparison <- prey_summary_filtered %>%
  group_by(Treatment, Species) %>%
  summarise(
    n_individuals = n(),
    mean_encounter_rate = mean(avg_daily_encounter_rate, na.rm = TRUE),
    sd_encounter_rate = sd(avg_daily_encounter_rate, na.rm = TRUE),
    median_encounter_rate = median(avg_daily_encounter_rate, na.rm = TRUE),
    min_encounter_rate = min(avg_daily_encounter_rate, na.rm = TRUE),
    max_encounter_rate = max(avg_daily_encounter_rate, na.rm = TRUE),
    .groups = "drop"
  )

print(treatment_comparison)

## 7.2 Statistical Tests ----
# Example: Wilcoxon test for treatment differences
if ("Treatment" %in% names(prey_summary_filtered)) {
  treatments <- unique(prey_summary_filtered$Treatment)
  
  if (length(treatments) == 2) {
    for (sp in c("Roach", "Perch")) {
      species_data <- prey_summary_filtered %>% filter(Species == sp)
      
      if (nrow(species_data) > 0) {
        test_result <- wilcox.test(
          avg_daily_encounter_rate ~ Treatment,
          data = species_data
        )
        
        message(sprintf("\n%s - Wilcoxon test p-value: %.4f", sp, test_result$p.value))
      }
    }
  }
}


# ===================================================================
# 8. CREATE SUMMARY TABLES
# ===================================================================

message("Creating summary tables...")

## 8.1 Suspected Predation Table ----
predation_table <- suspected_predation %>%
  select(Prey_ID, Species, Pike_ID, num_days_100_encounters, 
         consecutive_days_100, first_date_over_100, likely_predated) %>%
  arrange(desc(consecutive_days_100)) %>%
  flextable() %>%
  set_header_labels(
    Prey_ID = "Prey ID",
    Species = "Species",
    Pike_ID = "Pike ID",
    num_days_100_encounters = "Days ≥100 Encounters",
    consecutive_days_100 = "Consecutive Days",
    first_date_over_100 = "First Date",
    likely_predated = "Likely Predated"
  ) %>%
  theme_booktabs() %>%
  autofit()

# Save table
save_as_docx(predation_table,
             path = file.path(paths$tables, "suspected_predation_events.docx"))


## 8.2 Treatment Comparison Table ----
treatment_table <- treatment_comparison %>%
  flextable() %>%
  set_header_labels(
    Treatment = "Treatment",
    Species = "Species",
    n_individuals = "N",
    mean_encounter_rate = "Mean Daily Encounters",
    sd_encounter_rate = "SD",
    median_encounter_rate = "Median",
    min_encounter_rate = "Min",
    max_encounter_rate = "Max"
  ) %>%
  colformat_double(j = c("mean_encounter_rate", "sd_encounter_rate"), digits = 2) %>%
  theme_booktabs() %>%
  autofit()

save_as_docx(treatment_table,
             path = file.path(paths$tables, "treatment_encounter_comparison.docx"))


# ===================================================================
# 9. VISUALIZATION
# ===================================================================

message("Creating visualizations...")

## 9.1 Encounter Rate Distribution ----
library(ggplot2)

encounter_dist_plot <- ggplot(prey_summary_filtered, 
                              aes(x = Treatment, y = avg_daily_encounter_rate, 
                                  fill = Species)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~Species) +
  labs(
    title = "Daily Predator Encounter Rate by Treatment",
    x = "Treatment",
    y = "Average Daily Encounters",
    fill = "Species"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "lightgray"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(
  filename = file.path(paths$encounters, "encounter_plots/treatment_comparison.png"),
  plot = encounter_dist_plot,
  width = 10, height = 6, dpi = 300
)


## 9.2 Time Series of High-Risk Individual ----
# Example for one suspected predation event
if (nrow(suspected_predation) > 0) {
  # Select individual with highest consecutive encounters
  top_prey <- suspected_predation %>%
    arrange(desc(consecutive_days_100)) %>%
    slice(1)
  
  # Extract their daily data
  prey_id <- top_prey$Prey_ID
  pike_id <- top_prey$Pike_ID
  
  time_series_data <- prey_daily_encounters %>%
    filter(Prey_ID == prey_id, Pike_ID == pike_id)
  
  # Calculate scaling factor for dual y-axis
  scale_factor <- max(time_series_data$daily_avg_dist, na.rm = TRUE) / 
    max(time_series_data$encounter_count, na.rm = TRUE)
  
  ts_plot <- ggplot(time_series_data, aes(x = Date)) +
    geom_bar(aes(y = encounter_count * scale_factor), 
             stat = "identity", fill = "grey", alpha = 0.5) +
    geom_line(aes(y = daily_avg_dist), color = "black", linewidth = 1) +
    geom_point(aes(y = daily_avg_dist), color = "black", size = 2) +
    geom_text(aes(y = encounter_count * scale_factor, label = encounter_count),
              vjust = -0.5, size = 3) +
    scale_y_continuous(
      name = "Average Distance from Pike (m)",
      limits = c(0, max(c(time_series_data$daily_avg_dist, 
                          time_series_data$encounter_count * scale_factor), 
                        na.rm = TRUE) * 1.1)
    ) +
    labs(
      title = sprintf("Suspected Predation: %s by %s", prey_id, pike_id),
      subtitle = sprintf("Species: %s | Consecutive days >100 encounters: %d",
                         top_prey$Species, top_prey$consecutive_days_100),
      x = "Date"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  ggsave(
    filename = file.path(paths$encounters, 
                         sprintf("encounter_plots/%s_%s_timeline.png", prey_id, pike_id)),
    plot = ts_plot,
    width = 12, height = 8, units = "cm", dpi = 300
  )
}


# ===================================================================
# 10. EXPORT RESULTS
# ===================================================================

message("Exporting final results...")

# Export key datasets to Excel for easy review
library(openxlsx)

# Create workbook
wb <- createWorkbook()

# Add sheets
addWorksheet(wb, "Encounter Summary")
writeData(wb, "Encounter Summary", prey_summary_filtered)

addWorksheet(wb, "Suspected Predation")
writeData(wb, "Suspected Predation", suspected_predation)

addWorksheet(wb, "Treatment Comparison")
writeData(wb, "Treatment Comparison", treatment_comparison)

# Save workbook
saveWorkbook(wb, 
             file = file.path(paths$tables, "BT_predation_analysis_results.xlsx"),
             overwrite = TRUE)

message("Analysis complete! Check output directories for results.")


# ===================================================================
# 11. SESSION INFO
# ===================================================================

# Document R session for reproducibility
sink(file.path(paths$tables, "session_info.txt"))
print(sessionInfo())
sink()