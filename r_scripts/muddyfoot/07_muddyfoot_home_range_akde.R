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
# - Pre-flight alignment checks (projection + ID matching)
#
# REFERENCES:
# Fleming et al. (2022) Methods Ecol Evol - Population-level home range inference
# Silva et al. (2022) Methods Ecol Evol - AKDE methodology review
# Hollins et al. (2025) Methods Ecol Evol - Boundary spillover correction
#
# AUTHOR: [Your Name]
# DATE: February 2026
#===============================================================================

#===============================================================================
# LIBRARIES
#===============================================================================

library(ctmm)
library(sf)
library(ggplot2)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

#===============================================================================
# PATHS
#===============================================================================

ctmm_path    <- "./data/ctmm_fits/"
telem_path   <- "./data/telem_obj/muddyfoot/"
polygon_path <- "./data/lake_params/polygons/"
akde_path    <- "./data/akdes/"
figure_path  <- "./figures/muddyfoot/"
tables_path  <- "./tables/muddyfoot/"

dir.create(akde_path,  recursive = TRUE, showWarnings = FALSE)
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_path, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(akde_path, "muddyfoot_pike_akdes/"),  showWarnings = FALSE)
dir.create(paste0(akde_path, "muddyfoot_perch_akdes/"), showWarnings = FALSE)
dir.create(paste0(akde_path, "muddyfoot_roach_akdes/"), showWarnings = FALSE)

#===============================================================================
# LOAD DATA
#===============================================================================

message("=== Loading Data ===\n")

# Telemetry
pike_muddyfoot_tel  <- readRDS(paste0(telem_path, "pike_muddyfoot_tel_thinned_final.rds"))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, "perch_muddyfoot_tel_thinned_final.rds"))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, "roach_muddyfoot_tel_thinned_final.rds"))

message("Pike:  ", length(pike_muddyfoot_tel),  " individuals")
message("Perch: ", length(perch_muddyfoot_tel), " individuals")
message("Roach: ", length(roach_muddyfoot_tel), " individuals")

# ctmm fits
pike_muddyfoot_ctmm_fits  <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_best_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_best_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_best_models.rds"))

# Lake boundary
muddyfoot_polygon <- st_read(paste0(polygon_path, "muddyfoot_polygon_tpeqd.gpkg"), quiet = TRUE)
muddyfoot_sp_data <- as(muddyfoot_polygon, "Spatial")

message("\nPolygon CRS: ", st_crs(muddyfoot_polygon)$input, "\n")

#===============================================================================
# PRE-FLIGHT CHECKS: TELEMETRY vs CTMM ALIGNMENT
#===============================================================================

message("=== Pre-flight Alignment Checks ===\n")

for(species in c("perch", "roach", "pike")) {
  tel  <- get(paste0(species, "_muddyfoot_tel"))
  fits <- get(paste0(species, "_muddyfoot_ctmm_fits"))
  
  cat("Species:", species, "\n")
  
  # Check 1: IDs match
  tel_ids  <- names(tel)
  fit_ids  <- names(fits)
  ids_match <- identical(tel_ids, fit_ids)
  cat("  IDs match:", ids_match, "\n")
  if(!ids_match) {
    cat("  !!! Mismatched IDs:", setdiff(tel_ids, fit_ids), "\n")
  }
  
  # Check 2: Projections match
  proj_mismatches <- 0
  for(i in seq_along(tel)) {
    tel_proj <- ctmm::projection(tel[[i]])
    fit_proj <- ctmm::projection(fits[[i]])
    if(!identical(tel_proj, fit_proj)) {
      proj_mismatches <- proj_mismatches + 1
      cat("  !!! Projection mismatch at", names(tel)[i], ":\n")
      cat("      Tel:", tel_proj, "\n")
      cat("      Fit:", fit_proj, "\n")
    }
  }
  if(proj_mismatches == 0) cat("  Projections match: TRUE\n")
  cat("\n")
}

#===============================================================================
# FUNCTION: SELECT REFERENCE INDIVIDUAL
#===============================================================================
#
# Selects the individual with the highest effective sample size (N̂ = T/τ)
# weighted by temporal coverage, for use as the AKDE grid reference.

select_reference_individual <- function(tel_list, fits_list, min_N = 0.5) {
  
  n <- length(tel_list)
  
  diag <- data.frame(
    Individual         = names(tel_list),
    n_locations        = sapply(tel_list, nrow),
    duration_days      = sapply(tel_list, function(x) as.numeric(diff(range(x$timestamp)), units = "days")),
    model_type         = sapply(fits_list, function(x) summary(x)$name),
    tau_hours          = NA_real_,
    tau_days           = NA_real_,
    N_effective        = NA_real_,
    temporal_coverage  = sapply(tel_list, function(x) length(unique(as.Date(x$timestamp)))),
    stringsAsFactors   = FALSE
  )
  
  tau_patterns <- list(
    hours   = c("τ[position] (hours)",   "tau[position] (hours)"),
    minutes = c("τ[position] (minutes)", "tau[position] (minutes)"),
    days    = c("τ[position] (days)",    "tau[position] (days)"),
    seconds = c("τ[position] (seconds)", "tau[position] (seconds)")
  )
  
  for(i in seq_len(n)) {
    fit <- fits_list[[i]]
    sum_fit <- summary(fit)
    
    if(grepl("IID", sum_fit$name, ignore.case = TRUE)) {
      diag$N_effective[i] <- diag$n_locations[i]
      next
    }
    
    ci <- sum_fit$CI
    tau_found <- FALSE
    
    for(unit in names(tau_patterns)) {
      for(pat in tau_patterns[[unit]]) {
        if(pat %in% rownames(ci)) {
          tau_val <- ci[pat, "est"]
          tau_h   <- switch(unit,
                            hours   = tau_val,
                            minutes = tau_val / 60,
                            days    = tau_val * 24,
                            seconds = tau_val / 3600)
          diag$tau_hours[i] <- tau_h
          diag$tau_days[i]  <- tau_h / 24
          if(is.finite(tau_h / 24) && tau_h / 24 > 0)
            diag$N_effective[i] <- diag$duration_days[i] / (tau_h / 24)
          else
            diag$N_effective[i] <- 0
          tau_found <- TRUE
          break
        }
      }
      if(tau_found) break
    }
    
    if(!tau_found) diag$N_effective[i] <- 0
  }
  
  diag$quality_score <- diag$N_effective * diag$temporal_coverage
  
  message("\n========================================")
  message("    MOVEMENT MODEL DIAGNOSTICS")
  message("========================================\n")
  print(diag)
  
  range_resident <- !is.na(diag$N_effective) & diag$N_effective >= min_N
  message(sprintf("\nRange-resident (N̂ ≥ %.1f): %d / %d", min_N, sum(range_resident), n))
  
  if(sum(range_resident) == 0) {
    warning("No individuals meet N̂ >= ", min_N, " threshold. Using highest N̂.")
    best_idx <- which.max(diag$N_effective)
  } else {
    valid <- which(range_resident)
    best_idx <- valid[which.max(diag$quality_score[valid])]
  }
  
  message("\n========================================")
  message("    SELECTED REFERENCE INDIVIDUAL")
  message("========================================")
  message(sprintf("Individual:  %s",       diag$Individual[best_idx]))
  message(sprintf("Model:       %s",       diag$model_type[best_idx]))
  message(sprintf("N locations: %d",       diag$n_locations[best_idx]))
  message(sprintf("Duration:    %.1f days",diag$duration_days[best_idx]))
  message(sprintf("τ[position]: %.2f hrs", diag$tau_hours[best_idx]))
  message(sprintf("N̂:           %.1f",     diag$N_effective[best_idx]))
  message("========================================\n")
  
  return(best_idx)
}

#===============================================================================
# FUNCTION: RUN SPECIES AKDE ANALYSIS
#===============================================================================
#
# Runs the full AKDE pipeline for one species:
# - Reference individual selection and spillover calculation
# - Parallel individual AKDEs with boundary correction
# - Population-level PKDEs by treatment group
# - Meta-analysis for treatment comparison

run_species_akde <- function(tel_list,
                             fits_list,
                             sp_data,
                             control_idx,
                             exposed_idx,
                             species_name,
                             akde_save_path,
                             min_N = 3) {
  
  message("\n")
  message("###############################################################################")
  message(sprintf("#  %s HOME RANGE ANALYSIS", toupper(species_name)))
  message("###############################################################################\n")
  
  # Split by treatment
  control_tel  <- tel_list[control_idx]
  exposed_tel  <- tel_list[exposed_idx]
  control_fits <- fits_list[control_idx]
  exposed_fits <- fits_list[exposed_idx]
  
  message("Control group: ", length(control_tel),  " individuals")
  message("Exposed group: ", length(exposed_tel), " individuals\n")
  
  # Print treatment assignments
  treatments <- sapply(tel_list, function(x) unique(x$treatment))
  message("Treatment assignments:")
  print(treatments)
  cat("\n")
  
  # Reference individual ---------------------------------------------------
  ref_idx  <- select_reference_individual(tel_list, fits_list, min_N = min_N)
  ref_tel  <- tel_list[[ref_idx]]
  ref_fit  <- fits_list[[ref_idx]]
  ref_name <- names(tel_list)[ref_idx]
  
  message("Calculating reference AKDE with boundary correction...")
  
  ref_akde_bounded   <- akde(ref_tel, ref_fit, SP = sp_data, SP.in = TRUE)
  ref_akde_unbounded <- akde(ref_tel, ref_fit)
  
  spillover_pct <- (summary(ref_akde_unbounded)$CI[2] - summary(ref_akde_bounded)$CI[2]) /
    summary(ref_akde_unbounded)$CI[2] * 100
  
  message(sprintf("Bounded area:   %.2f m²", summary(ref_akde_bounded)$CI[2]))
  message(sprintf("Unbounded area: %.2f m²", summary(ref_akde_unbounded)$CI[2]))
  message(sprintf("Spillover:      %.1f%%\n", spillover_pct))
  
  # Individual AKDEs -------------------------------------------------------
  message("Calculating individual AKDEs...")
  message("Using ", detectCores() - 1, " cores\n")
  
  n_workers <- length(tel_list)
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  
  akde_list <- foreach(i = seq_along(tel_list),
                       .packages   = "ctmm",
                       .errorhandling = "pass") %dopar% {
                         akde(tel_list[[i]], fits_list[[i]], SP = sp_data, SP.in = TRUE)
                       }
  
  stopCluster(cl)
  
  names(akde_list) <- names(tel_list)
  
  # Check failures
  failures <- sapply(akde_list, inherits, "error")
  if(any(failures)) {
    warning(sprintf("%d AKDE(s) failed: %s",
                    sum(failures),
                    paste(names(akde_list)[failures], collapse = ", ")))
    akde_list <- akde_list[!failures]
    tel_list  <- tel_list[!failures]
    fits_list <- fits_list[!failures]
  }
  
  message(sprintf("Successfully calculated %d individual AKDEs\n", length(akde_list)))
  
  # Save individual AKDEs
  saveRDS(akde_list, paste0(akde_save_path, tolower(species_name), "_muddyfoot_akdes.rds"))
  
  # Re-split after potential failure removal
  # (uses names to ensure correct assignment even if some failed)
  control_names <- names(tel_list[control_idx])
  exposed_names <- names(tel_list[exposed_idx])
  
  control_akdes <- akde_list[names(akde_list) %in% control_names]
  exposed_akdes <- akde_list[names(akde_list) %in% exposed_names]
  control_tel   <- tel_list[names(tel_list) %in% control_names]
  exposed_tel   <- tel_list[names(tel_list) %in% exposed_names]
  
  # Population-level PKDEs -------------------------------------------------
  message("Calculating population-level PKDEs...\n")
  
  control_PKDE <- pkde(control_tel, control_akdes, SP = sp_data, SP.in = TRUE)
  exposed_PKDE <- pkde(exposed_tel, exposed_akdes, SP = sp_data, SP.in = TRUE)
  total_PKDE   <- pkde(tel_list,    akde_list,     SP = sp_data, SP.in = TRUE)
  
  message(sprintf("Control PKDE: %.2f m² (%.2f - %.2f)",
                  summary(control_PKDE)$CI[2],
                  summary(control_PKDE)$CI[1],
                  summary(control_PKDE)$CI[3]))
  message(sprintf("Exposed PKDE: %.2f m² (%.2f - %.2f)",
                  summary(exposed_PKDE)$CI[2],
                  summary(exposed_PKDE)$CI[1],
                  summary(exposed_PKDE)$CI[3]))
  message(sprintf("Total   PKDE: %.2f m² (%.2f - %.2f)\n",
                  summary(total_PKDE)$CI[2],
                  summary(total_PKDE)$CI[1],
                  summary(total_PKDE)$CI[3]))
  
  saveRDS(control_PKDE, paste0(akde_save_path, tolower(species_name), "_control_PKDE.rds"))
  saveRDS(exposed_PKDE, paste0(akde_save_path, tolower(species_name), "_mix_PKDE.rds"))
  saveRDS(total_PKDE,   paste0(akde_save_path, tolower(species_name), "_total_PKDE.rds"))
  
  # Meta-analysis ----------------------------------------------------------
  message("Performing meta-analysis for treatment comparison...\n")
  
  meta_list <- list(Control = control_akdes, Exposed = exposed_akdes)
  meta_result <- ctmm::meta(meta_list, sort = FALSE, level = 0.95)
  print(meta_result)
  
  message(sprintf("\n=== %s Analysis Complete ===\n", species_name))
  
  # Return results
  list(
    akdes         = akde_list,
    control_akdes = control_akdes,
    exposed_akdes = exposed_akdes,
    control_PKDE  = control_PKDE,
    exposed_PKDE  = exposed_PKDE,
    total_PKDE    = total_PKDE,
    meta          = meta_result,
    spillover_pct = spillover_pct
  )
}

#===============================================================================
# PIKE ANALYSIS
#===============================================================================

pike_results <- run_species_akde(
  tel_list      = pike_muddyfoot_tel,
  fits_list     = pike_muddyfoot_ctmm_fits,
  sp_data       = muddyfoot_sp_data,
  control_idx   = 1:3,
  exposed_idx   = 4:6,
  species_name  = "Pike",
  akde_save_path = paste0(akde_path, "muddyfoot_pike_akdes/"),
  min_N         = 3
)

pike_muddyfoot_akdes  <- pike_results$akdes
pike_control_akdes    <- pike_results$control_akdes
pike_mix_akdes        <- pike_results$exposed_akdes
pike_control_PKDE     <- pike_results$control_PKDE
pike_mix_PKDE         <- pike_results$exposed_PKDE
pike_total_PKDE       <- pike_results$total_PKDE

#===============================================================================
# PERCH ANALYSIS
#===============================================================================

perch_results <- run_species_akde(
  tel_list      = perch_muddyfoot_tel,
  fits_list     = perch_muddyfoot_ctmm_fits,
  sp_data       = muddyfoot_sp_data,
  control_idx   = 1:15,
  exposed_idx   = 16:30,
  species_name  = "Perch",
  akde_save_path = paste0(akde_path, "muddyfoot_perch_akdes/"),
  min_N         = 3
)

perch_muddyfoot_akdes <- perch_results$akdes
perch_control_akdes   <- perch_results$control_akdes
perch_mix_akdes       <- perch_results$exposed_akdes
perch_control_PKDE    <- perch_results$control_PKDE
perch_mix_PKDE        <- perch_results$exposed_PKDE
perch_total_PKDE      <- perch_results$total_PKDE

#===============================================================================
# ROACH ANALYSIS
#===============================================================================

roach_results <- run_species_akde(
  tel_list      = roach_muddyfoot_tel,
  fits_list     = roach_muddyfoot_ctmm_fits,
  sp_data       = muddyfoot_sp_data,
  control_idx   = 1:13,
  exposed_idx   = 14:26,
  species_name  = "Roach",
  akde_save_path = paste0(akde_path, "muddyfoot_roach_akdes/"),
  min_N         = 3
)

roach_muddyfoot_akdes <- roach_results$akdes
roach_control_akdes   <- roach_results$control_akdes
roach_mix_akdes       <- roach_results$exposed_akdes
roach_control_PKDE    <- roach_results$control_PKDE
roach_mix_PKDE        <- roach_results$exposed_PKDE
roach_total_PKDE      <- roach_results$total_PKDE

#===============================================================================
# SUMMARY TABLE
#===============================================================================

message("\n")
message("###############################################################################")
message("#                          ANALYSIS SUMMARY                                   #")
message("###############################################################################\n")

get_ci <- function(pkde_obj, idx) summary(pkde_obj)$CI[idx]

summary_data <- data.frame(
  Species      = rep(c("Pike", "Perch", "Roach"), each = 3),
  Group        = rep(c("Control", "Exposed", "Total"), 3),
  Mean_m2      = c(
    get_ci(pike_control_PKDE, 2),  get_ci(pike_mix_PKDE, 2),  get_ci(pike_total_PKDE, 2),
    get_ci(perch_control_PKDE, 2), get_ci(perch_mix_PKDE, 2), get_ci(perch_total_PKDE, 2),
    get_ci(roach_control_PKDE, 2), get_ci(roach_mix_PKDE, 2), get_ci(roach_total_PKDE, 2)
  ),
  CI_Lower     = c(
    get_ci(pike_control_PKDE, 1),  get_ci(pike_mix_PKDE, 1),  get_ci(pike_total_PKDE, 1),
    get_ci(perch_control_PKDE, 1), get_ci(perch_mix_PKDE, 1), get_ci(perch_total_PKDE, 1),
    get_ci(roach_control_PKDE, 1), get_ci(roach_mix_PKDE, 1), get_ci(roach_total_PKDE, 1)
  ),
  CI_Upper     = c(
    get_ci(pike_control_PKDE, 3),  get_ci(pike_mix_PKDE, 3),  get_ci(pike_total_PKDE, 3),
    get_ci(perch_control_PKDE, 3), get_ci(perch_mix_PKDE, 3), get_ci(perch_total_PKDE, 3),
    get_ci(roach_control_PKDE, 3), get_ci(roach_mix_PKDE, 3), get_ci(roach_total_PKDE, 3)
  ),
  Spillover_pct = c(
    pike_results$spillover_pct,  NA, NA,
    perch_results$spillover_pct, NA, NA,
    roach_results$spillover_pct, NA, NA
  )
)

print(summary_data)

write.csv(summary_data,
          paste0(tables_path, "home_range_summary.csv"),
          row.names = FALSE)

message("\nSummary saved to: ", tables_path, "home_range_summary.csv")
message("\n###############################################################################")
message("#                        ANALYSIS COMPLETE                                    #")
message("###############################################################################\n")