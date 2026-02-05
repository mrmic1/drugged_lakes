#===============================================================================
# DIAGNOSTIC SCRIPT - INVESTIGATE TAU EXTRACTION
#===============================================================================

# Test with first pike individual
fit <- pike_muddyfoot_ctmm_fits[[1]]

cat("=== RAW TAU OBJECT ===\n")
print(fit$tau)
cat("\n")

cat("=== TAU POSITION EXTRACTION ===\n")
tau_pos <- fit$tau["position"]
print(tau_pos)
cat("Class:", class(tau_pos), "\n")
cat("Units attribute:", attr(tau_pos, "units"), "\n\n")

cat("=== UNIT CONVERSIONS ===\n")
cat("as.numeric(tau_pos):", as.numeric(tau_pos), "\n")
cat("as.numeric(tau_pos, units='hours'):", as.numeric(tau_pos, units='hours'), "\n")
cat("as.numeric(tau_pos, units='days'):", as.numeric(tau_pos, units='days'), "\n")
cat("as.numeric(tau_pos, units='seconds'):", as.numeric(tau_pos, units='seconds'), "\n\n")

cat("=== SUMMARY METHOD ===\n")
print(summary(fit)$CI)
cat("\n")

cat("=== TELEMETRY TIME RANGE ===\n")
tel <- pike_muddyfoot_tel[[1]]
time_diff <- diff(range(tel$timestamp))
print(time_diff)
cat("as.numeric in days:", as.numeric(time_diff, units="days"), "\n")
cat("as.numeric in hours:", as.numeric(time_diff, units="hours"), "\n\n")

cat("=== EFFECTIVE SAMPLE SIZE CALCULATION ===\n")
tau_hours <- as.numeric(tau_pos, units = "hours")
tau_days <- as.numeric(tau_pos, units = "days")
duration_days <- as.numeric(time_diff, units = "days")

cat("Duration (days):", duration_days, "\n")
cat("Tau (hours):", tau_hours, "\n")
cat("Tau (days):", tau_days, "\n")
cat("N_eff = duration/tau:", duration_days / tau_days, "\n\n")

# Check if there's a conversion issue
cat("=== CHECKING FOR UNIT MISMATCH ===\n")
# The correct tau should be around 1.49 hours based on summary
expected_tau_hours <- 1.49
cat("Expected tau (from summary):", expected_tau_hours, "hours\n")
cat("Calculated tau:", tau_hours, "hours\n")
cat("Ratio:", tau_hours / expected_tau_hours, "\n")


#===============================================================================
# ROBUST REFERENCE INDIVIDUAL SELECTION - HANDLES TAU EXTRACTION PROPERLY
#===============================================================================

diagnose_movement_models_robust <- function(tel_list, fits_list) {
  
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
  
  # Calculate effective sample size with robust unit handling
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
    
    # For OU models, extract tau[position] using summary (most reliable)
    if("position" %in% rownames(summary(fit)$CI)) {
      
      # METHOD 1: Use summary CI table (most reliable)
      # The summary table shows tau in the units specified in the model
      tau_from_summary <- summary(fit)$CI["τ[position]", "est"]
      
      # Check the DOF table to understand the time units
      # DOF table typically uses hours for tau
      tau_hours <- tau_from_summary
      tau_days <- tau_hours / 24
      
      diagnostics$tau_hours[i] <- tau_hours
      diagnostics$tau_days[i] <- tau_days
      
      # Calculate N_effective
      if(is.finite(tau_days) && tau_days > 0) {
        diagnostics$N_effective[i] <- diagnostics$duration_days[i] / tau_days
      } else {
        diagnostics$N_effective[i] <- 0
      }
      
    } else if("position" %in% names(fit$tau)) {
      
      # METHOD 2: Direct extraction from fit object
      tau_pos <- fit$tau["position"]
      
      # Try to extract the numeric value properly
      # The issue is that fit$tau may store values in different internal units
      tau_hours <- as.numeric(tau_pos, units = "hours")
      tau_days <- tau_hours / 24
      
      # Sanity check: tau should be reasonable relative to tracking duration
      # If tau > tracking duration, there may be a unit issue
      if(tau_days > diagnostics$duration_days[i]) {
        warning(sprintf("Individual %s: tau (%.1f days) > tracking duration (%.1f days). Check units!",
                        diagnostics$Individual[i], tau_days, diagnostics$duration_days[i]))
      }
      
      diagnostics$tau_hours[i] <- tau_hours
      diagnostics$tau_days[i] <- tau_days
      
      # Calculate N_effective
      if(is.finite(tau_days) && tau_days > 0) {
        diagnostics$N_effective[i] <- diagnostics$duration_days[i] / tau_days
      } else {
        diagnostics$N_effective[i] <- 0
      }
      
    } else {
      # No position timescale available
      diagnostics$N_effective[i] <- 0
    }
  }
  
  # Calculate quality score
  diagnostics$quality_score <- diagnostics$N_effective * diagnostics$temporal_coverage_days
  
  return(diagnostics)
}

select_reference_individual_robust <- function(tel_list, fits_list, min_N = 3) {
  
  # Run diagnostics
  diag <- diagnose_movement_models_robust(tel_list, fits_list)
  
  # Print diagnostic table
  message("\n========================================")
  message("    MOVEMENT MODEL DIAGNOSTICS")
  message("========================================\n")
  print(diag)
  
  # Check for range residency
  range_resident <- !is.na(diag$N_effective) & 
    diag$N_effective >= min_N
  
  message("\n----------------------------------------")
  message(sprintf("Range-resident individuals (N̂ ≥ %.1f): %d / %d",
                  min_N, sum(range_resident), length(tel_list)))
  message("----------------------------------------\n")
  
  # If no individuals meet min_N, show what we have and use best available
  if(sum(range_resident) == 0) {
    warning("\nNo individuals meet the N̂ >= ", min_N, " threshold.")
    message("\nThis suggests:")
    message("  1. Pike have very large home ranges (τ >> tracking duration)")
    message("  2. Pike are nomadic/migratory during this period")
    message("  3. Short tracking duration relative to movement scale")
    message("\nProceeding with best available individual...")
    message("Note: Home range estimates will have HIGH UNCERTAINTY\n")
    
    # Use individual with highest N_effective, even if < min_N
    best_idx <- which.max(diag$N_effective)
    
  } else {
    # Warn about low sample sizes
    low_N_count <- sum(range_resident & diag$N_effective < 10)
    if(low_N_count > 0) {
      warning(sprintf("\n%d individuals have N̂ < 10. Home range estimates will have high uncertainty.",
                      low_N_count))
    }
    
    # Select individual with highest quality score
    valid_indices <- which(range_resident)
    best_idx <- valid_indices[which.max(diag$quality_score[valid_indices])]
  }
  
  # Print selection summary
  message("========================================")
  message("    SELECTED REFERENCE INDIVIDUAL")
  message("========================================")
  message(sprintf("Individual:         %s", diag$Individual[best_idx]))
  message(sprintf("Model:              %s", diag$model_type[best_idx]))
  message(sprintf("Locations:          %d", diag$n_locations[best_idx]))
  message(sprintf("Duration:           %.1f days", diag$duration_days[best_idx]))
  message(sprintf("Temporal coverage:  %d days", diag$temporal_coverage_days[best_idx]))
  message(sprintf("τ[position]:        %.2f hours (%.3f days)", 
                  diag$tau_hours[best_idx],
                  diag$tau_days[best_idx]))
  message(sprintf("N̂ (effective):      %.2f home range crossings", 
                  diag$N_effective[best_idx]))
  message(sprintf("Quality score:      %.1f", diag$quality_score[best_idx]))
  
  # Add interpretation
  if(diag$N_effective[best_idx] < 1) {
    message("\n⚠️  WARNING: N̂ < 1 indicates fish did not complete one full home range crossing")
    message("   Home range size estimates will be HIGHLY UNCERTAIN and likely underestimated")
  } else if(diag$N_effective[best_idx] < 3) {
    message("\n⚠️  CAUTION: N̂ < 3 indicates limited home range crossings")
    message("   Home range estimates will have considerable uncertainty")
  }
  
  message("========================================\n")
  
  return(best_idx)
}

#===============================================================================
# TEST WITH PIKE DATA
#===============================================================================

message("Testing with pike data...\n")

pike_ref_idx <- select_reference_individual_robust(pike_muddyfoot_tel, 
                                                   pike_muddyfoot_ctmm_fits,
                                                   min_N = 3)

# Show the summary for the selected individual
message("\n=== Detailed Summary of Selected Individual ===\n")
print(summary(pike_muddyfoot_ctmm_fits[[pike_ref_idx]]))
