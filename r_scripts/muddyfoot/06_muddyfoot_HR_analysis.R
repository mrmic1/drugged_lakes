#===============================================================================
# MUDDYFOOT LAKE - FISH HOME RANGE ANALYSIS WITH AKDE AND BOUNDARY CORRECTION
#===============================================================================
# 
# PURPOSE:
# Estimate individual and population-level home ranges for pike, perch, and 
# roach using autocorrelated kernel density estimation (AKDE) with proper
# boundary enforcement for lake shorelines.
#
# KEY FEATURES:
# - Data-driven reference individual selection based on effective sample size
# - Local boundary correction (Hollins et al. 2025)
# - Population-level inference with meta-analysis (Fleming et al. 2022)
# - Treatment group comparisons (Control vs Exposed)
# - Parallel processing with error handling
#
# REFERENCES:
# Fleming et al. (2022) Methods Ecol Evol - Population-level home range inference
# Silva et al. (2022) Methods Ecol Evol - AKDE methodology review
# Hollins et al. (2025) Methods Ecol Evol - Boundary spillover correction
#
# AUTHOR: [Your Name]
# DATE: February 2026
#===============================================================================

# Load required libraries ---------------------------------------------------
library(ctmm)
library(sf)
library(ggplot2)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

#===============================================================================
# SETUP - DEFINE PATHS
#===============================================================================

# Input paths
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/muddyfoot/"
lake_polygon_path <- "./data/lake_params/polygons/"

# Output paths
akde_path <- "./data/akdes/"
figure_path <- "./figures/muddyfoot/"
save_tables_path <- "./tables/muddyfoot/"

# Create output directories if they don't exist
dir.create(akde_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
dir.create(save_tables_path, recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# LOAD DATA
#===============================================================================

message("=== Loading Data ===\n")

# Load telemetry objects
pike_muddyfoot_tel <- readRDS(paste0(telem_path, "pike_muddyfoot_tel_thinned_final.rds"))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, "perch_muddyfoot_tel_thinned_final.rds"))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, "roach_muddyfoot_tel_thinned_final.rds"))

message("Pike: ", length(pike_muddyfoot_tel), " individuals")
message("Perch: ", length(perch_muddyfoot_tel), " individuals")
message("Roach: ", length(roach_muddyfoot_tel), " individuals")

# Load movement model fits
pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_best_ctmm_model_fits.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_best_ctmm_model_fits.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_best_ctmm_model_fits.rds"))

# Load lake boundary polygon
muddyfoot_polygon <- st_read(paste0(lake_polygon_path, "muddyfoot_polygon_tpeqd.gpkg"), 
                             quiet = TRUE)
muddyfoot_sp_data <- as(muddyfoot_polygon, "Spatial")

message("\nLake boundary polygon loaded")
message("Polygon CRS: ", st_crs(muddyfoot_polygon)$input, "\n")

#===============================================================================
# FUNCTION: DIAGNOSE MOVEMENT MODELS
#===============================================================================
# 
# PURPOSE:
# Evaluate quality of movement model fits by calculating effective sample size
# (N̂ = T/τ) for each individual. Identifies range-resident vs nomadic individuals.
#
# PARAMETERS:
# - tel_list: List of telemetry objects
# - fits_list: List of fitted ctmm models
#
# RETURNS:
# Data frame with diagnostic metrics for each individual
#===============================================================================

diagnose_movement_models <- function(tel_list, fits_list) {
  
  n_individuals <- length(tel_list)
  
  diagnostics <- data.frame(
    Individual = names(tel_list),
    n_locations = sapply(tel_list, nrow),
    duration_days = sapply(tel_list, function(tel) {
      as.numeric(diff(range(tel$timestamp)), units = "days")
    }),
    model_type = sapply(fits_list, function(fit) {
      summary(fit)$name
    }),
    tau_hours = NA_real_,
    tau_days = NA_real_,
    N_effective = NA_real_,
    temporal_coverage_days = sapply(tel_list, function(tel) {
      length(unique(as.Date(tel$timestamp)))
    }),
    stringsAsFactors = FALSE
  )
  
  # Calculate effective sample size with proper unit handling
  for(i in 1:n_individuals) {
    fit <- fits_list[[i]]
    model_name <- summary(fit)$name
    
    # For IID models (no autocorrelation), N ≈ number of locations
    if(grepl("IID", model_name, ignore.case = TRUE)) {
      diagnostics$N_effective[i] <- diagnostics$n_locations[i]
      diagnostics$tau_hours[i] <- NA
      diagnostics$tau_days[i] <- NA
      next
    }
    
    # For OU models, extract tau[position] and convert to days
    if("position" %in% names(fit$tau)) {
      tau_pos <- fit$tau["position"]
      
      # Convert to hours for display
      tau_hours <- as.numeric(tau_pos, units = "hours")
      diagnostics$tau_hours[i] <- tau_hours
      
      # Convert to days for N calculation
      tau_days <- as.numeric(tau_pos, units = "days")
      diagnostics$tau_days[i] <- tau_days
      
      # Check if tau is valid (finite and positive)
      if(is.finite(tau_days) && tau_days > 0) {
        # N̂ = T / τ
        diagnostics$N_effective[i] <- diagnostics$duration_days[i] / tau_days
      } else {
        diagnostics$N_effective[i] <- 0
      }
    } else {
      # No position timescale (shouldn't happen for OU models)
      diagnostics$N_effective[i] <- 0
    }
  }
  
  # Calculate quality score: N̂ × temporal coverage
  diagnostics$quality_score <- diagnostics$N_effective * diagnostics$temporal_coverage_days
  
  return(diagnostics)
}

#===============================================================================
# FUNCTION: SELECT REFERENCE INDIVIDUAL
#===============================================================================
#
# PURPOSE:
# Select the best individual to use as reference for grid spacing (dr parameter)
# in AKDE calculations. Selection is based on effective sample size (N̂) and
# temporal coverage, following Fleming et al. (2022) recommendations.
#
# RATIONALE:
# - High N̂ (T/τ) indicates individual crossed its home range many times
# - High temporal coverage indicates consistent tracking over many days
# - Product of these metrics identifies highest-quality individual
#
# PARAMETERS:
# - tel_list: List of telemetry objects
# - fits_list: List of fitted ctmm models
# - min_N: Minimum effective sample size for range residency (default = 3)
#
# RETURNS:
# Index of selected reference individual
#===============================================================================

select_reference_individual_final <- function(tel_list, fits_list, min_N = 0.5) {
  
  n_individuals <- length(tel_list)
  
  diagnostics <- data.frame(
    Individual = names(tel_list),
    n_locations = sapply(tel_list, nrow),
    duration_days = sapply(tel_list, function(tel) {
      as.numeric(diff(range(tel$timestamp)), units = "days")
    }),
    model_type = sapply(fits_list, function(fit) {
      summary(fit)$name
    }),
    tau_hours = NA_real_,
    tau_days = NA_real_,
    N_effective = NA_real_,
    temporal_coverage_days = sapply(tel_list, function(tel) {
      length(unique(as.Date(tel$timestamp)))
    }),
    stringsAsFactors = FALSE
  )
  
  # Extract tau from summary CI table (most reliable)
  for(i in 1:n_individuals) {
    fit <- fits_list[[i]]
    sum_fit <- summary(fit)
    model_name <- sum_fit$name
    
    # For IID models, N ≈ number of locations
    if(grepl("IID", model_name, ignore.case = TRUE)) {
      diagnostics$N_effective[i] <- diagnostics$n_locations[i]
      next
    }
    
    # For OU models, extract from CI table
    ci_table <- sum_fit$CI
    
    # Look for tau[position] row with different time units
    tau_row_patterns <- list(
      hours = c("τ[position] (hours)", "tau[position] (hours)"),
      minutes = c("τ[position] (minutes)", "tau[position] (minutes)"),
      days = c("τ[position] (days)", "tau[position] (days)"),
      seconds = c("τ[position] (seconds)", "tau[position] (seconds)")
    )
    
    tau_found <- FALSE
    tau_hours <- NA
    
    # Try each unit type
    for(unit_type in names(tau_row_patterns)) {
      for(pattern in tau_row_patterns[[unit_type]]) {
        if(pattern %in% rownames(ci_table)) {
          tau_value <- ci_table[pattern, "est"]
          
          # Convert to hours
          tau_hours <- switch(unit_type,
                              hours = tau_value,
                              minutes = tau_value / 60,
                              days = tau_value * 24,
                              seconds = tau_value / 3600)
          
          tau_found <- TRUE
          break
        }
      }
      if(tau_found) break
    }
    
    if(tau_found && !is.na(tau_hours)) {
      tau_days <- tau_hours / 24
      
      diagnostics$tau_hours[i] <- tau_hours
      diagnostics$tau_days[i] <- tau_days
      
      if(is.finite(tau_days) && tau_days > 0) {
        diagnostics$N_effective[i] <- diagnostics$duration_days[i] / tau_days
      } else {
        diagnostics$N_effective[i] <- 0
      }
    } else {
      # Couldn't find tau
      diagnostics$N_effective[i] <- 0
    }
  }
  
  # Calculate quality score
  diagnostics$quality_score <- diagnostics$N_effective * diagnostics$temporal_coverage_days
  
  # Print diagnostics
  message("\n========================================")
  message("    MOVEMENT MODEL DIAGNOSTICS")
  message("========================================\n")
  print(diagnostics)
  
  # Statistics
  max_N <- max(diagnostics$N_effective, na.rm = TRUE)
  
  message("\n----------------------------------------")
  message(sprintf("Maximum N̂ observed: %.2f", max_N))
  message(sprintf("Range-resident threshold: N̂ ≥ %.1f", min_N))
  
  range_resident <- !is.na(diagnostics$N_effective) & 
    diagnostics$N_effective >= min_N
  
  message(sprintf("Range-resident individuals: %d / %d",
                  sum(range_resident), length(tel_list)))
  message("----------------------------------------\n")
  
  # Select best individual
  if(sum(range_resident) == 0) {
    warning("\nNo individuals meet the N̂ >= ", min_N, " threshold.")
    message("\n📊 INTERPRETATION:")
    message("   Fish have very large home ranges relative to tracking period")
    message("   Average τ ~ ", round(mean(diagnostics$tau_days, na.rm=TRUE), 1), 
            " days vs tracking duration ~ ", round(mean(diagnostics$duration_days), 1), " days")
    message("\n   This is biologically plausible for large predatory fish in lakes")
    message("   AKDE can still estimate home ranges, but with higher uncertainty")
    message("\n   Proceeding with individual with highest N̂...\n")
    
    best_idx <- which.max(diagnostics$N_effective)
  } else {
    valid_indices <- which(range_resident)
    best_idx <- valid_indices[which.max(diagnostics$quality_score[valid_indices])]
  }
  
  # Print selection
  message("========================================")
  message("    SELECTED REFERENCE INDIVIDUAL")
  message("========================================")
  message(sprintf("Individual:         %s", diagnostics$Individual[best_idx]))
  message(sprintf("Model:              %s", diagnostics$model_type[best_idx]))
  message(sprintf("Locations:          %d", diagnostics$n_locations[best_idx]))
  message(sprintf("Duration:           %.1f days", diagnostics$duration_days[best_idx]))
  message(sprintf("Temporal coverage:  %d days", diagnostics$temporal_coverage_days[best_idx]))
  message(sprintf("τ[position]:        %.2f hours (%.3f days)", 
                  diagnostics$tau_hours[best_idx],
                  diagnostics$tau_days[best_idx]))
  message(sprintf("N̂ (effective):      %.1f home range crossings", 
                  diagnostics$N_effective[best_idx]))
  message(sprintf("Quality score:      %.1f", diagnostics$quality_score[best_idx]))
  
  # Interpretation
  N_val <- diagnostics$N_effective[best_idx]
  if(N_val >= 10) {
    message("\n✓ EXCELLENT: N̂ ≥ 10 - High quality home range estimate expected")
  } else if(N_val >= 3) {
    message("\n✓ GOOD: N̂ ≥ 3 - Reliable home range estimate")
  } else if(N_val >= 1) {
    message("\n⚠️  CAUTION: 1 ≤ N̂ < 3 - Moderate uncertainty in home range estimate")
  } else if(N_val >= 0.5) {
    message("\n⚠️  N̂ < 1: Fish did not complete one full home range crossing")
    message("   - Home range estimate will be uncertain and possibly underestimated")
  } else {
    message("\n⚠️  N̂ < 0.5: Fish completed less than half a home range crossing")
    message("   - Home range is much larger than area covered during tracking")
    message("   - Estimates represent 'occurrence distribution' not true home range")
  }
  
  message("========================================\n")
  
  return(best_idx)
}

#===============================================================================
# PIKE ANALYSIS
#===============================================================================

message("\n")
message("###############################################################################")
message("#                         PIKE HOME RANGE ANALYSIS                            #")
message("###############################################################################\n")

# Separate by treatment group -----------------------------------------------

# Extract treatment information from telemetry objects
pike_treatments <- sapply(pike_muddyfoot_tel, function(tel) {
  unique(tel$treatment)
})
print(pike_treatments)

pike_control_tel <- pike_muddyfoot_tel[1:3]  # Update indices based on your data
pike_mix_tel <- pike_muddyfoot_tel[4:6]      # Update indices based on your data

pike_control_fits <- pike_muddyfoot_ctmm_fits[1:3]
pike_mix_fits <- pike_muddyfoot_ctmm_fits[4:6]

message("Control group: ", length(pike_control_tel), " individuals")
message("Exposed group: ", length(pike_mix_tel), " individuals\n")

# Select reference individual -----------------------------------------------
pike_ref_idx <- select_reference_individual_final(pike_muddyfoot_tel, 
                                            pike_muddyfoot_ctmm_fits,
                                            min_N = 3)

pike_ref_tel <- pike_muddyfoot_tel[[pike_ref_idx]]
pike_ref_fit <- pike_muddyfoot_ctmm_fits[[pike_ref_idx]]
pike_ref_name <- names(pike_muddyfoot_tel)[pike_ref_idx]

# Calculate reference AKDE with boundary ------------------------------------
message("Calculating reference AKDE with boundary correction...")

pike_ref_akde_bounded <- akde(pike_ref_tel, 
                              pike_ref_fit,
                              SP = muddyfoot_sp_data,
                              SP.in = TRUE)

pike_ref_akde_unbounded <- akde(pike_ref_tel, 
                                pike_ref_fit)

# Calculate spillover percentage
pike_spillover_pct <- (summary(pike_ref_akde_unbounded)$CI[2] - 
                         summary(pike_ref_akde_bounded)$CI[2]) / 
  summary(pike_ref_akde_unbounded)$CI[2] * 100

message(sprintf("Reference AKDE complete"))
message(sprintf("Bounded area:   %.2f m² (95%% CI: %.2f - %.2f)",
                summary(pike_ref_akde_bounded)$CI[2],
                summary(pike_ref_akde_bounded)$CI[1],
                summary(pike_ref_akde_bounded)$CI[3]))
message(sprintf("Unbounded area: %.2f m²", summary(pike_ref_akde_unbounded)$CI[2]))
message(sprintf("Spillover:      %.1f%%\n", pike_spillover_pct))

# Extract grid spacing from reference AKDE ----------------------------------
pike_dr <- pike_ref_akde_bounded$dr

message(sprintf("Grid spacing (dr): %.3f m\n", pike_dr))

# Calculate individual AKDEs in parallel ------------------------------------
message("Calculating individual AKDEs with boundary correction...")
message("Using parallel processing with ", detectCores() - 1, " cores\n")

# Setup parallel processing
cl <- makeCluster(detectCores() - 3)
registerDoParallel(cl)

# CRITICAL: Export all necessary objects to workers
clusterExport(cl, c("pike_muddyfoot_tel", 
                    "pike_muddyfoot_ctmm_fits", 
                    "muddyfoot_sp_data"))

# Now run the foreach loop
pike_akdes <- foreach(i = 1:length(pike_muddyfoot_tel),
                      .packages = c('ctmm'),
                      .errorhandling = 'pass') %dopar% {
                        
                        akde(pike_muddyfoot_tel[[i]],
                             pike_muddyfoot_ctmm_fits[[i]],
                             SP = muddyfoot_sp_data,
                             SP.in = TRUE)
                      }
                      

stopCluster(cl)

# Check for failures
names(pike_akdes) <- names(pike_muddyfoot_tel)
failures <- sapply(pike_akdes, function(x) inherits(x, "error"))

if(any(failures)) {
  warning(sprintf("%d AKDE calculations failed:", sum(failures)))
  warning(paste(names(pike_akdes)[failures], collapse = ", "))
  
  # Remove failures
  pike_akdes <- pike_akdes[!failures]
  pike_muddyfoot_tel <- pike_muddyfoot_tel[!failures]
  pike_muddyfoot_ctmm_fits <- pike_muddyfoot_ctmm_fits[!failures]
}

message(sprintf("Successfully calculated %d individual AKDEs\n", length(pike_akdes)))

# Save individual AKDEs -----------------------------------------------------
saveRDS(pike_akdes, paste0(akde_path, "muddyfoot_pike_akdes/pike_muddyfoot_akdes.rds"))
message("Individual AKDEs saved\n")

# Separate by treatment group (after removing failures) ---------------------
# NOTE: Update these indices if any calculations failed
pike_control_akdes <- pike_akdes[1:3]
pike_mix_akdes <- pike_akdes[4:6]


# Population-level estimates ------------------------------------------------
message("Calculating population-level home ranges (PKDE)...\n")

# Control group PKDE
pike_control_PKDE <- pkde(pike_control_tel,
                          pike_control_akdes,
                          SP = muddyfoot_sp_data,
                          SP.in = TRUE)

# Exposed group PKDE
pike_mix_PKDE <- pkde(pike_mix_tel,
                      pike_mix_akdes,
                      SP = muddyfoot_sp_data,
                      SP.in = TRUE)

# Overall population PKDE
pike_total_PKDE <- pkde(pike_muddyfoot_tel,
                        pike_akdes,
                        SP = muddyfoot_sp_data,
                        SP.in = TRUE)

message("Control group PKDE: ", 
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)",
                summary(pike_control_PKDE)$CI[2],
                summary(pike_control_PKDE)$CI[1],
                summary(pike_control_PKDE)$CI[3]))

message("Exposed group PKDE: ",
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)",
                summary(pike_mix_PKDE)$CI[2],
                summary(pike_mix_PKDE)$CI[1],
                summary(pike_mix_PKDE)$CI[3]))

message("Total population PKDE: ",
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)\n",
                summary(pike_total_PKDE)$CI[2],
                summary(pike_total_PKDE)$CI[1],
                summary(pike_total_PKDE)$CI[3]))

saveRDS(pike_control_PKDE, paste0(akde_path, "muddyfoot_pike_akdes/pike_control_PKDE.rds"))
saveRDS(pike_mix_PKDE, paste0(akde_path, "muddyfoot_pike_akdes/pike_mix_PKDE.rds"))
saveRDS(pike_total_PKDE, paste0(akde_path, "muddyfoot_pike_akdes/pike_total_PKDE.rds"))


# Meta-analysis for treatment comparison ------------------------------------
message("Performing meta-analysis for treatment comparison...\n")

pike_akde_total <- list(Control = pike_control_akdes, 
                        Exposed = pike_mix_akdes)

pike_akde_meta_data <- ctmm::meta(pike_akde_total, 
                            sort = FALSE, 
                            level = 0.95)

print(pike_akde_meta_data)
summary(pike_akde_meta_data)

#Get individual group estimates
pike_meta_control <- ctmm::meta(pike_control_akdes, level = 0.95)
pike_meta_exposed <- ctmm::meta(pike_mix_akdes, level = 0.95)

print(pike_meta_control)
print(pike_meta_exposed)

message("\n=== Pike Analysis Complete ===\n")

#===============================================================================
# PERCH ANALYSIS
#===============================================================================

message("\n")
message("###############################################################################")
message("#                        PERCH HOME RANGE ANALYSIS                            #")
message("###############################################################################\n")

# Separate by treatment group -----------------------------------------------

# Extract treatment information from telemetry objects
perch_treatments <- sapply(perch_muddyfoot_tel, function(tel) {
  unique(tel$treatment)
})
print(perch_treatments)


# Update these indices based on your data structure
perch_control_tel <- perch_muddyfoot_tel[1:15]
perch_mix_tel <- perch_muddyfoot_tel[16:30]

perch_control_fits <- perch_muddyfoot_ctmm_fits[1:15]
perch_mix_fits <- perch_muddyfoot_ctmm_fits[16:30]

message("Control group: ", length(perch_control_tel), " individuals")
message("Exposed group: ", length(perch_mix_tel), " individuals\n")

# Select reference individual -----------------------------------------------
perch_ref_idx <- select_reference_individual_final(perch_muddyfoot_tel, 
                                             perch_muddyfoot_ctmm_fits,
                                             min_N = 3)

perch_ref_tel <- perch_muddyfoot_tel[[perch_ref_idx]]
perch_ref_fit <- perch_muddyfoot_ctmm_fits[[perch_ref_idx]]
perch_ref_name <- names(perch_muddyfoot_tel)[perch_ref_idx]

# Calculate reference AKDE with boundary ------------------------------------
message("Calculating reference AKDE with boundary correction...")

perch_ref_akde_bounded <- akde(perch_ref_tel, 
                               perch_ref_fit,
                               SP = muddyfoot_sp_data,
                               SP.in = TRUE)

perch_ref_akde_unbounded <- akde(perch_ref_tel, 
                                 perch_ref_fit)

# Calculate spillover percentage
perch_spillover_pct <- (summary(perch_ref_akde_unbounded)$CI[2] - 
                          summary(perch_ref_akde_bounded)$CI[2]) / 
  summary(perch_ref_akde_unbounded)$CI[2] * 100

message(sprintf("Reference AKDE complete"))
message(sprintf("Bounded area:   %.2f m² (95%% CI: %.2f - %.2f)",
                summary(perch_ref_akde_bounded)$CI[2],
                summary(perch_ref_akde_bounded)$CI[1],
                summary(perch_ref_akde_bounded)$CI[3]))
message(sprintf("Unbounded area: %.2f m²", summary(perch_ref_akde_unbounded)$CI[2]))
message(sprintf("Spillover:      %.1f%%\n", perch_spillover_pct))

# Extract grid spacing
perch_dr <- perch_ref_akde_bounded$dr
message(sprintf("Grid spacing (dr): %.3f m\n", perch_dr))

# Calculate individual AKDEs in parallel ------------------------------------
message("Calculating individual AKDEs with boundary correction...")
message("Using parallel processing with ", detectCores() - 1, " cores\n")

cl <- makeCluster(detectCores() - 3)
registerDoParallel(cl)

perch_akdes <- foreach(i = 1:length(perch_muddyfoot_tel),
                       .packages = c('ctmm'),
                       .errorhandling = 'pass') %dopar% {
                         
                         akde(perch_muddyfoot_tel[[i]],
                              perch_muddyfoot_ctmm_fits[[i]],
                              SP = muddyfoot_sp_data,
                              SP.in = TRUE)
                       }

stopCluster(cl)

# Check for failures
names(perch_akdes) <- names(perch_muddyfoot_tel)
failures <- sapply(perch_akdes, function(x) inherits(x, "error"))

if(any(failures)) {
  warning(sprintf("%d AKDE calculations failed:", sum(failures)))
  warning(paste(names(perch_akdes)[failures], collapse = ", "))
  
  perch_akdes <- perch_akdes[!failures]
  perch_muddyfoot_tel <- perch_muddyfoot_tel[!failures]
  perch_muddyfoot_ctmm_fits <- perch_muddyfoot_ctmm_fits[!failures]
}

message(sprintf("Successfully calculated %d individual AKDEs\n", length(perch_akdes)))

# Save individual AKDEs
saveRDS(perch_akdes, paste0(akde_path, "muddyfoot_perch_akdes/perch_muddyfoot_akdes.rds"))
message("Individual AKDEs saved\n")

# Separate by treatment group (update if failures occurred)
perch_control_akdes <- perch_akdes[1:15]
perch_mix_akdes <- perch_akdes[16:30]

perch_control_tel <- perch_muddyfoot_tel[1:15]
perch_mix_tel <- perch_muddyfoot_tel[16:30]

# Population-level estimates ------------------------------------------------
message("Calculating population-level home ranges (PKDE)...\n")

perch_control_PKDE <- pkde(perch_control_tel,
                           perch_control_akdes,
                           SP = muddyfoot_sp_data,
                           SP.in = TRUE)

perch_mix_PKDE <- pkde(perch_mix_tel,
                       perch_mix_akdes,
                       SP = muddyfoot_sp_data,
                       SP.in = TRUE)

perch_total_PKDE <- pkde(perch_muddyfoot_tel,
                         perch_akdes,
                         SP = muddyfoot_sp_data,
                         SP.in = TRUE)

message("Control group PKDE: ", 
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)",
                summary(perch_control_PKDE)$CI[2],
                summary(perch_control_PKDE)$CI[1],
                summary(perch_control_PKDE)$CI[3]))

message("Exposed group PKDE: ",
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)",
                summary(perch_mix_PKDE)$CI[2],
                summary(perch_mix_PKDE)$CI[1],
                summary(perch_mix_PKDE)$CI[3]))

message("Total population PKDE: ",
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)\n",
                summary(perch_total_PKDE)$CI[2],
                summary(perch_total_PKDE)$CI[1],
                summary(perch_total_PKDE)$CI[3]))

saveRDS(perch_control_PKDE, paste0(akde_path, "muddyfoot_perch_akdes/perch_control_PKDE.rds"))
saveRDS(perch_mix_PKDE, paste0(akde_path, "muddyfoot_perch_akdes/perch_mix_PKDE.rds"))
saveRDS(perch_total_PKDE, paste0(akde_path, "muddyfoot_perch_akdes/perch_total_PKDE.rds"))

# Meta-analysis
message("Performing meta-analysis for treatment comparison...\n")

perch_akde_total <- list(Control = perch_control_akdes, 
                         Exposed = perch_mix_akdes)

perch_akde_meta_data <- ctmm::meta(perch_akde_total, 
                             sort = FALSE, 
                             level = 0.95)

print(perch_akde_meta_data)

message("\n=== Perch Analysis Complete ===\n")

#===============================================================================
# ROACH ANALYSIS
#===============================================================================

message("\n")
message("###############################################################################")
message("#                        ROACH HOME RANGE ANALYSIS                            #")
message("###############################################################################\n")

# Separate by treatment group -----------------------------------------------
# Extract treatment information from telemetry objects
roach_treatments <- sapply(roach_muddyfoot_tel, function(tel) {
  unique(tel$treatment)
})
print(roach_treatments)


# Update these indices based on your data structure
roach_control_tel <- roach_muddyfoot_tel[1:13]
roach_mix_tel <- roach_muddyfoot_tel[14:26]

roach_control_fits <- roach_muddyfoot_ctmm_fits[1:13]
roach_mix_fits <- roach_muddyfoot_ctmm_fits[14:26]

message("Control group: ", length(roach_control_tel), " individuals")
message("Exposed group: ", length(roach_mix_tel), " individuals\n")

# Select reference individual -----------------------------------------------
roach_ref_idx <- select_reference_individual_final(roach_muddyfoot_tel, 
                                             roach_muddyfoot_ctmm_fits,
                                             min_N = 3)

roach_ref_tel <- roach_muddyfoot_tel[[roach_ref_idx]]
roach_ref_fit <- roach_muddyfoot_ctmm_fits[[roach_ref_idx]]
roach_ref_name <- names(roach_muddyfoot_tel)[roach_ref_idx]

# Calculate reference AKDE with boundary ------------------------------------
message("Calculating reference AKDE with boundary correction...")

roach_ref_akde_bounded <- akde(roach_ref_tel, 
                               roach_ref_fit,
                               SP = muddyfoot_sp_data,
                               SP.in = TRUE)

roach_ref_akde_unbounded <- akde(roach_ref_tel, 
                                 roach_ref_fit)

# Calculate spillover percentage
roach_spillover_pct <- (summary(roach_ref_akde_unbounded)$CI[2] - 
                          summary(roach_ref_akde_bounded)$CI[2]) / 
  summary(roach_ref_akde_unbounded)$CI[2] * 100

message(sprintf("Reference AKDE complete"))
message(sprintf("Bounded area:   %.2f m² (95%% CI: %.2f - %.2f)",
                summary(roach_ref_akde_bounded)$CI[2],
                summary(roach_ref_akde_bounded)$CI[1],
                summary(roach_ref_akde_bounded)$CI[3]))
message(sprintf("Unbounded area: %.2f m²", summary(roach_ref_akde_unbounded)$CI[2]))
message(sprintf("Spillover:      %.1f%%\n", roach_spillover_pct))

# Extract grid spacing
roach_dr <- roach_ref_akde_bounded$dr
message(sprintf("Grid spacing (dr): %.3f m\n", roach_dr))

# Calculate individual AKDEs in parallel ------------------------------------
message("Calculating individual AKDEs with boundary correction...")
message("Using parallel processing with ", detectCores() - 3, " cores\n")

cl <- makeCluster(detectCores() - 10)
registerDoParallel(cl)
# Export all necessary objects to parallel workers
clusterExport(cl, varlist = c("roach_muddyfoot_tel", 
                              "roach_muddyfoot_ctmm_fits", 
                              "muddyfoot_sp_data"))

roach_akdes <- foreach(i = 1:length(roach_muddyfoot_tel),
                       .packages = c('ctmm')) %dopar% {
                         
                         akde(roach_muddyfoot_tel[[i]],
                              roach_muddyfoot_ctmm_fits[[i]],
                              SP = muddyfoot_sp_data,
                              SP.in = TRUE)
                       }

stopCluster(cl)

# Check for failures
names(roach_akdes) <- names(roach_muddyfoot_tel)
failures <- sapply(roach_akdes, function(x) inherits(x, "error"))

if(any(failures)) {
  warning(sprintf("%d AKDE calculations failed:", sum(failures)))
  warning(paste(names(roach_akdes)[failures], collapse = ", "))
  
  roach_akdes <- roach_akdes[!failures]
  roach_muddyfoot_tel <- roach_muddyfoot_tel[!failures]
  roach_muddyfoot_ctmm_fits <- roach_muddyfoot_ctmm_fits[!failures]
}

message(sprintf("Successfully calculated %d individual AKDEs\n", length(roach_akdes)))

# Save individual AKDEs
saveRDS(roach_akdes, paste0(akde_path, "roach_muddyfoot_akdes.rds"))
message("Individual AKDEs saved\n")

# Separate by treatment group (update if failures occurred)
roach_control_akdes <- roach_akdes[1:13]
roach_mix_akdes <- roach_akdes[14:26]

# Population-level estimates ------------------------------------------------
message("Calculating population-level home ranges (PKDE)...\n")

roach_control_PKDE <- pkde(roach_control_tel,
                           roach_control_akdes,
                           SP = muddyfoot_sp_data,
                           SP.in = TRUE)

roach_mix_PKDE <- pkde(roach_mix_tel,
                       roach_mix_akdes,
                       SP = muddyfoot_sp_data,
                       SP.in = TRUE)

roach_total_PKDE <- pkde(roach_muddyfoot_tel,
                         roach_akdes,
                         SP = muddyfoot_sp_data,
                         SP.in = TRUE)

message("Control group PKDE: ", 
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)",
                summary(roach_control_PKDE)$CI[2],
                summary(roach_control_PKDE)$CI[1],
                summary(roach_control_PKDE)$CI[3]))

message("Exposed group PKDE: ",
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)",
                summary(roach_mix_PKDE)$CI[2],
                summary(roach_mix_PKDE)$CI[1],
                summary(roach_mix_PKDE)$CI[3]))

message("Total population PKDE: ",
        sprintf("%.2f m² (95%% CI: %.2f - %.2f)\n",
                summary(roach_total_PKDE)$CI[2],
                summary(roach_total_PKDE)$CI[1],
                summary(roach_total_PKDE)$CI[3]))

saveRDS(roach_control_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/roach_control_PKDE.rds"))
saveRDS(roach_mix_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/roach_mix_PKDE.rds"))
saveRDS(roach_total_PKDE, paste0(akde_path, "muddyfoot_roach_akdes/roach_total_PKDE.rds"))

# Meta-analysis
message("Performing meta-analysis for treatment comparison...\n")

roach_akde_total <- list(Control = roach_control_akdes, 
                         Exposed = roach_mix_akdes)

roach_akde_meta_data <- meta(roach_akde_total, 
                             sort = FALSE, 
                             level = 0.95)

print(roach_akde_meta_data)

message("\n=== Roach Analysis Complete ===\n")

#===============================================================================
# SUMMARY
#===============================================================================

message("\n")
message("###############################################################################")
message("#                          ANALYSIS SUMMARY                                   #")
message("###############################################################################\n")

# Create summary table
summary_data <- data.frame(
  Species = rep(c("Pike", "Perch", "Roach"), each = 3),
  Group = rep(c("Control", "Exposed", "Total"), 3),
  Mean_Area_m2 = c(
    summary(pike_control_PKDE)$CI[2],
    summary(pike_mix_PKDE)$CI[2],
    summary(pike_total_PKDE)$CI[2],
    summary(perch_control_PKDE)$CI[2],
    summary(perch_mix_PKDE)$CI[2],
    summary(perch_total_PKDE)$CI[2],
    summary(roach_control_PKDE)$CI[2],
    summary(roach_mix_PKDE)$CI[2],
    summary(roach_total_PKDE)$CI[2]
  ),
  CI_Lower = c(
    summary(pike_control_PKDE)$CI[1],
    summary(pike_mix_PKDE)$CI[1],
    summary(pike_total_PKDE)$CI[1],
    summary(perch_control_PKDE)$CI[1],
    summary(perch_mix_PKDE)$CI[1],
    summary(perch_total_PKDE)$CI[1],
    summary(roach_control_PKDE)$CI[1],
    summary(roach_mix_PKDE)$CI[1],
    summary(roach_total_PKDE)$CI[1]
  ),
  CI_Upper = c(
    summary(pike_control_PKDE)$CI[3],
    summary(pike_mix_PKDE)$CI[3],
    summary(pike_total_PKDE)$CI[3],
    summary(perch_control_PKDE)$CI[3],
    summary(perch_mix_PKDE)$CI[3],
    summary(perch_total_PKDE)$CI[3],
    summary(roach_control_PKDE)$CI[3],
    summary(roach_mix_PKDE)$CI[3],
    summary(roach_total_PKDE)$CI[3]
))

print(summary_data)

# Save summary table
write.csv(summary_data, 
          paste0(save_tables_path, "home_range_summary.csv"),
          row.names = FALSE)

message("\nSummary table saved to: ", save_tables_path, "home_range_summary.csv")

message("\n###############################################################################")
message("#                        ANALYSIS COMPLETE                                    #")
message("###############################################################################\n")