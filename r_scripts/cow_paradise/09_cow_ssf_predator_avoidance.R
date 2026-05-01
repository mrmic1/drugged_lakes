# ================================================================ #
# DYNAMIC SSF PREDATOR AVOIDANCE ANALYSIS - COW PARADISE             #
# ================================================================ #
# Based on Schlägel et al. (2019) Methods in Ecology and Evolution
#
# Approach: Use rolling occurrence distributions (ODs) of pike
# as dynamic covariates in prey step-selection functions (SSFs).
# A negative SSF coefficient = prey actively moves away from
# areas of recent pike space use = predator avoidance.
# Compare control vs. exposed treatments within each prey species.
# ================================================================ #

library(ctmm)
library(tidyverse)
library(amt)
library(raster)
library(sf)
library(parallel)
library(foreach)
library(doParallel)
library(terra)

# -------------------- #
# 1. DIRECTORIES ####
# -------------------- #

telem_path   <- "./data/telem_obj/cow_paradise/"
akde_path    <- "./data/akdes/"
ssf_path     <- "./data/ssfs/predators/"
fig_path     <- "./figures/cow_paradise/"

# Create output directories
dir.create(paste0(ssf_path, "cow_perch/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(ssf_path, "cow_roach/"), recursive = TRUE, showWarnings = FALSE)

# -------------------- #
# 2. LOAD DATA ####
# -------------------- #

# Telemetry objects (ctmm format)
pike_cow_tel   <- readRDS(paste0(telem_path, "pike_cow_paradise_tel_thinned_final.rds"))
perch_cow_tel  <- readRDS(paste0(telem_path, "perch_cow_paradise_tel_thinned_final.rds"))
roach_cow_tel  <- readRDS(paste0(telem_path, "roach_cow_paradise_tel_thinned_final.rds"))

# ------------------------------------------ #
# 3. CONVERT TELEMETRY TO AMT TRACKS ####
# ------------------------------------------ #
# Following the same approach as previous RSF analysis:
# convert WGS84 coordinates to UTM 33N before making tracks.
# ctmm telemetry objects store raw coordinates as longitude/latitude
# in the 'longitude' and 'latitude' columns (not x/y which are tpeqd)

# Check column names available in telemetry objects
cat("Pike telemetry columns:", names(pike_cow_tel[[1]]), "\n")

# Function to convert ctmm telemetry to UTM 33N data frame
ctmm_tel_to_utm <- function(tel_obj, id_label, extra_cols = NULL) {
  
  df <- data.frame(
    x_wgs = tel_obj$longitude,
    y_wgs = tel_obj$latitude,
    t     = as.POSIXct(tel_obj$timestamp),
    id    = id_label
  )
  
  if (!is.null(extra_cols)) {
    for (col in names(extra_cols)) {
      df[[col]] <- extra_cols[[col]]
    }
  }
  
  df_utm <- df %>%
    sf::st_as_sf(coords = c("x_wgs", "y_wgs"), crs = 4326) %>%
    sf::st_transform(32633) %>%
    dplyr::mutate(
      x_ = sf::st_coordinates(geometry)[, 1],
      y_ = sf::st_coordinates(geometry)[, 2]
    ) %>%
    sf::st_drop_geometry()
  
  df_utm
}

# Convert pike
pike_df <- bind_rows(
  lapply(names(pike_cow_tel), function(id) {
    ctmm_tel_to_utm(pike_cow_tel[[id]], id_label = id)
  })
)

# Convert perch (with treatment column)
perch_df <- bind_rows(
  lapply(names(perch_cow_tel), function(id) {
    tel <- perch_cow_tel[[id]]
    ctmm_tel_to_utm(tel, id_label = id,
                    extra_cols = list(treatment = unique(tel$treatment)))
  })
)

# Convert roach (with treatment column)
roach_df <- bind_rows(
  lapply(names(roach_cow_tel), function(id) {
    tel <- roach_cow_tel[[id]]
    ctmm_tel_to_utm(tel, id_label = id,
                    extra_cols = list(treatment = unique(tel$treatment)))
  })
)

# Verify coordinates look sensible (UTM 33N northings ~7,000,000 for Sweden)
cat("Pike x range (UTM easting):",  range(pike_df$x_), "\n")
cat("Pike y range (UTM northing):", range(pike_df$y_), "\n")

# Create amt track objects following the same pattern as previous RSF analysis
utm33_crs <- sf::st_crs(32633)

pike_track <- pike_df %>%
  transmute(x_ = x_, y_ = y_, t_ = t, id = id) %>%
  { class(.) <- c("track_xyt", "track_xy", "tbl_df", "tbl", "data.frame"); . } %>%
  { attr(., "crs_") <- utm33_crs; . }

perch_track <- perch_df %>%
  transmute(x_ = x_, y_ = y_, t_ = t, id = id, treatment = treatment) %>%
  { class(.) <- c("track_xyt", "track_xy", "tbl_df", "tbl", "data.frame"); . } %>%
  { attr(., "crs_") <- utm33_crs; . }

roach_track <- roach_df %>%
  transmute(x_ = x_, y_ = y_, t_ = t, id = id, treatment = treatment) %>%
  { class(.) <- c("track_xyt", "track_xy", "tbl_df", "tbl", "data.frame"); . } %>%
  { attr(., "crs_") <- utm33_crs; . }

# Verify track structure
print(head(pike_track))
cat("Track CRS:", sf::st_crs(attr(pike_track, "crs_"))$input, "\n")
print(attr(pike_track, "crs_"))


# Also verify the track structure looks correct
cat("\nTrack class:", class(pike_track), "\n")
cat("Track dimensions:", nrow(pike_track), "rows,", ncol(pike_track), "cols\n")
cat("Column names:", names(pike_track), "\n")
cat("\nHead of pike track:\n")
print(head(pike_track))

# Verify coordinate ranges make sense for UTM 33N (Sweden)
cat("\nCoordinate ranges:\n")
cat("x (easting):  ", range(pike_track$x_), "\n")
cat("y (northing): ", range(pike_track$y_), "\n")


# Update template raster to UTM 33N
trast <- raster::raster(
  raster::extent(
    min(pike_df$x_) - 20,
    max(pike_df$x_) + 20,
    min(pike_df$y_) - 20,
    max(pike_df$y_) + 20
  ),
  res = 1
)
raster::crs(trast) <- sp::CRS("EPSG:32633")
trast[] <- 1
plot(trast, main = "Template raster - UTM 33N")

# Confirm template raster dimensions
cat("Template raster dimensions:", nrow(trast), "rows x", ncol(trast), "cols\n")
cat("Template raster extent:", as.vector(raster::extent(trast)), "\n")
cat("Total cells:", raster::ncell(trast), "\n")

# ------------------------------------------ #
# 4. BUILD TEMPLATE RASTER ####
# ------------------------------------------ #
# amt 0.3.1.0 uses terra SpatRaster, not raster RasterLayer

trast_terra <- terra::rast(
  xmin = min(pike_df$x_) - 20,
  xmax = max(pike_df$x_) + 20,
  ymin = min(pike_df$y_) - 20,
  ymax = max(pike_df$y_) + 20,
  resolution = 1,
  crs = "EPSG:32633"
)

terra::values(trast_terra) <- 1

# Verify
print(trast_terra)
terra::crs(trast_terra)
plot(trast_terra, main = "Template raster - UTM 33N (terra)")

# ------------------------------------------ #
# 5. COMPUTE ROLLING PIKE ODs ####
# ------------------------------------------ #

# Thin telemetry to 5-minute fixes to reduce computation time.
# 30-second fixes are temporally redundant for estimating pike
# occurrence distributions at the scale of this small lake.
thin_interval_mins <- 5

pike_df_thinned <- pike_df %>%
  dplyr::group_by(id) %>%                    # process each pike separately
  dplyr::arrange(t) %>%                      # ensure chronological order
  dplyr::filter(
    as.numeric(difftime(t,                   # time difference between current
                        dplyr::lag(t,        # and previous fix
                                   default = t[1]),     # first fix always kept
                        units = "mins")) >= thin_interval_mins |
      dplyr::row_number() == 1              # always keep first location
  ) %>%
  dplyr::ungroup()

cat("Thinned location counts per pike:\n")
pike_df_thinned %>% dplyr::count(id) %>% print()

# n_points_window: number of consecutive locations used to estimate
# each occurrence distribution. 12 x 5-min fixes = 60-min window,
# representing the timescale over which prey perceive pike presence.
n_points_window <- 12

# resolution: spatial grain of the occurrence distribution raster.
# 2m is sufficient for a ~40x30m lake and reduces computation time
# compared to 1m (4x fewer cells).
resolution <- 1

# Output directory for GeoTIFF files — one per pike individual.
# Saving as GeoTIFF avoids terra pointer issues that occur when
# SpatRaster objects are passed across parallel workers or saved
# with saveRDS().
od_dir <- paste0(ssf_path, "pike_ods/cow/")
dir.create(od_dir, recursive = TRUE, showWarnings = FALSE)

pike_ids <- unique(pike_df_thinned$id)

for (pid in pike_ids) {
  
  cat("Computing rolling OD for:", pid, "—", format(Sys.time(), "%H:%M:%S"), "\n")
  
  # Subset to single pike and rename t to t_ as required by amt
  pike_sub <- pike_df_thinned %>%
    dplyr::filter(id == pid) %>%
    dplyr::select(x_, y_, t) %>%
    dplyr::rename(t_ = t)
  
  # make_track(): converts data frame to amt track_xyt object.
  # Required input format for all amt movement functions.
  pike_track <- amt::make_track(
    pike_sub, .x = x_, .y = y_, .t = t_,
    crs = sf::st_crs(32633)             # UTM 33N coordinate system
  )
  
  # Build the template raster that defines the spatial grid on which
  # occurrence distributions are computed. Extent covers the full lake
  # plus a 20m buffer. All values set to 1 as a placeholder — the
  # actual OD values are computed by rolling_od().
  trast_worker <- terra::rast(
    xmin = min(pike_df_thinned$x_) - 20,
    xmax = max(pike_df_thinned$x_) + 20,
    ymin = min(pike_df_thinned$y_) - 20,
    ymax = max(pike_df_thinned$y_) + 20,
    resolution = resolution,
    crs = "EPSG:32633"
  )
  terra::values(trast_worker) <- 1
  
  # rolling_od(): core function from Schlägel et al. (2019) via amt.
  # For each time step k, uses locations k to k+n_points_window to
  # compute a kriged occurrence distribution (OD) via ctmm.
  # Returns a SpatRaster stack with one layer per time step,
  # where high values indicate areas the pike most likely occupied
  # during that time window. This becomes the dynamic predation
  # risk covariate in the prey SSF.
  result <- tryCatch(
    amt::rolling_od(pike_track, trast_worker,
                    n.points = n_points_window,  # locations per OD window
                    show.progress = TRUE),        # progress bar per pike
    error = function(e) {
      message("FAILED: ", pid, " — ", e$message)
      return(NULL)
    }
  )
  
  # writeRaster(): saves the SpatRaster stack to a GeoTIFF file.
  # This is done immediately inside the loop so results are preserved
  # even if a later pike fails. GeoTIFF is the correct format for
  # SpatRaster objects — unlike saveRDS(), it stores the actual raster
  # data on disk rather than relying on C++ memory pointers.
  if (!is.null(result)) {
    out_file <- paste0(od_dir, pid, "_rolling_od.tif")
    terra::writeRaster(result, out_file, overwrite = TRUE)
    cat("  Saved:", out_file, "\n")
    rm(result)   # remove from memory immediately to free RAM
    gc()         # trigger garbage collection
  }
}

# ---- Reload from GeoTIFF files ----
# terra::rast() loads each GeoTIFF back into R as a fresh SpatRaster
# with valid pointers. This is always safe regardless of how the
# files were originally created.
cat("\nReloading pike ODs from GeoTIFF files...\n")

od_files <- list.files(od_dir, pattern = "_rolling_od.tif$",
                       full.names = TRUE)
print(od_files)

pike_rolling_ods <- lapply(od_files, terra::rast)

# Extract pike IDs from filenames to use as list names
names(pike_rolling_ods) <- gsub("_rolling_od.tif", "",
                                basename(od_files))

# Verify pointers are valid after reload
cat("\nVerification:\n")
for (pid in names(pike_rolling_ods)) {
  cat(sprintf("  %-10s  layers: %d  pointer: %s\n",
              pid,
              terra::nlyr(pike_rolling_ods[[pid]]),
              tryCatch({ terra::nlyr(pike_rolling_ods[[pid]]); "VALID" },
                       error = function(e) "INVALID")))
}


# Check temporal coverage of the OD stack
od <- pike_rolling_ods$F59894

cat("Total layers:", terra::nlyr(od), "\n")

# Layer names contain the timestamps
layer_names <- names(od)
cat("First layer:", layer_names[1], "\n")
cat("Last layer: ", layer_names[terra::nlyr(od)], "\n")

# Or check via the thinned track directly
pike_df_thinned %>%
  dplyr::filter(id == "F59894") %>%
  dplyr::summarise(
    start = min(t),
    end   = max(t),
    n_locs = dplyr::n(),
    days  = round(difftime(max(t), min(t), units = "days"), 1)
  ) %>%
  print()


# ------------------------------------------ #
# 6. FIT PREY SSFs ####
# ------------------------------------------ #
# For each prey individual, fit an SSF where the covariate
# is the mean rolling OD across all pike (averaged to give
# a single "pike presence" surface at each time step).
#
# Steps:
# 1. Generate observed steps from track
# 2. Generate random available steps (n = 10 per observed)
# 3. Extract the pike OD value at each step endpoint
# 4. Fit conditional logistic regression

# Extract pike_od value for each pike individually, then average
# across pike. This preserves time attributes and handles different
# layer counts correctly.

fit_prey_ssf <- function(steps_xyt, pike_od_list, max_time_mins = 60) {
  
  tryCatch({
    
    pike_od_extractions <- lapply(names(pike_od_list), function(pid) {
      
      od <- pike_od_list[[pid]]
      
      if (is.null(terra::time(od))) {
        stop("No time attribute on OD for pike: ", pid)
      }
      
      amt::extract_covariates_var_time(
        steps_xyt,
        od,
        max_time   = lubridate::minutes(max_time_mins),
        when       = "any",
        name_covar = paste0("pike_od_", pid)
      )[[paste0("pike_od_", pid)]]
    })
    
    od_matrix <- do.call(cbind, pike_od_extractions)
    
    steps_xyt$pike_od <- -rowMeans(od_matrix, na.rm = TRUE)
    
    cat("    Non-NA pike_od values:",
        sum(!is.na(steps_xyt$pike_od)), "/", nrow(steps_xyt), "\n")
    
    steps_xyt %>%
      dplyr::filter(!is.na(pike_od)) %>%
      amt::fit_ssf(case_ ~ pike_od + strata(step_id_))
    
  }, error = function(e) {
    message("SSF failed: ", e$message)
    return(NULL)
  })
}

# ---- Rebuild perch steps and test on one individual ----
pid_test  <- perch_ids[1]
prey_sub  <- dplyr::filter(perch_track, id == pid_test)

steps_test <- prey_sub %>%
  amt::steps() %>%
  amt::random_steps(n_control = 10) %>%
  dplyr::mutate(
    t1_ = lubridate::with_tz(t1_, "UTC"),
    t2_ = lubridate::with_tz(t2_, "UTC")
  )

cat("Testing SSF for:", pid_test, "\n")
model_test <- fit_prey_ssf(steps_test, pike_rolling_ods)

if (!is.null(model_test)) {
  cat("SUCCESS\n")
  print(summary(model_test$model)$coefficients)
} else {
  cat("FAILED\n")
}

# ---- 6.1 Perch ####

perch_ids       <- unique(perch_track$id)
perch_treatments <- perch_df %>% distinct(id, treatment)

perch_ssf_results <- lapply(perch_ids, function(pid) {
  
  cat("Fitting SSF for perch:", pid, "\n")
  
  prey_sub  <- dplyr::filter(perch_track, id == pid)
  treatment <- perch_treatments$treatment[perch_treatments$id == pid]
  
  steps_xyt <- prey_sub %>%
    amt::steps() %>%
    amt::random_steps(n_control = 10) %>%
    dplyr::mutate(
      t1_ = lubridate::with_tz(t1_, "UTC"),
      t2_ = lubridate::with_tz(t2_, "UTC")
    )
  
  model <- fit_prey_ssf(steps_xyt, pike_rolling_ods)
  
  if (!is.null(model)) {
    saveRDS(model, paste0(ssf_path, "cow_perch/", pid, "_ssf.rds"))
  }
  
  list(id = pid, treatment = treatment, model = model)
})

names(perch_ssf_results) <- perch_ids
saveRDS(perch_ssf_results, paste0(ssf_path, "cow_perch/perch_ssf_results.rds"))



# ---- 6.2 Roach ####

roach_ids        <- unique(roach_track$id)
roach_treatments <- roach_df %>% distinct(id, treatment)

roach_ssf_results <- lapply(roach_ids, function(rid) {
  
  cat("Fitting SSF for roach:", rid, "\n")
  
  prey_sub  <- dplyr::filter(roach_track, id == rid)
  treatment <- roach_treatments$treatment[roach_treatments$id == rid]
  
  # Convert to steps first — this creates the case_ column
  steps_xyt <- prey_sub %>%
    amt::steps() %>%
    amt::random_steps(n_control = 10) %>%
    dplyr::mutate(
      t1_ = lubridate::with_tz(t1_, "UTC"),
      t2_ = lubridate::with_tz(t2_, "UTC")
    )
  
  model <- fit_prey_ssf(steps_xyt, pike_rolling_ods)
  
  if (!is.null(model)) {
    saveRDS(model, paste0(ssf_path, "cow_roach/", rid, "_ssf.rds"))
  }
  
  list(id = rid, treatment = treatment, model = model)
})

names(roach_ssf_results) <- roach_ids
saveRDS(roach_ssf_results, paste0(ssf_path, "cow_roach/roach_ssf_results.rds"))

# ------------------------------------------ #
# 7. EXTRACT SSF COEFFICIENTS ####
# ------------------------------------------ #
# Extract the pike_od selection coefficient for each individual.
# Negative = avoidance of pike space use
# Positive = attraction to pike space use

extract_ssf_coefs <- function(ssf_results, species_label) {
  
  bind_rows(lapply(ssf_results, function(res) {
    
    if (is.null(res$model)) return(NULL)
    
    coef_table <- tryCatch(
      as.data.frame(summary(res$model$model)$coefficients),
      error = function(e) NULL
    )
    
    if (is.null(coef_table)) return(NULL)
    
    # Get CI from model
    ci <- tryCatch(
      confint(res$model$model),
      error = function(e) NULL
    )
    
    data.frame(
      id        = res$id,
      species   = species_label,
      treatment = res$treatment,
      est       = coef_table["pike_od", "coef"],
      se        = coef_table["pike_od", "se(coef)"],
      low       = if (!is.null(ci)) ci["pike_od", 1] else NA,
      high      = if (!is.null(ci)) ci["pike_od", 2] else NA,
      p         = coef_table["pike_od", "Pr(>|z|)"]
    )
  }))
}

# Load if re-running from saved files
# perch_ssf_results <- readRDS(paste0(ssf_path, "cow_perch/perch_ssf_results.rds"))
# roach_ssf_results <- readRDS(paste0(ssf_path, "cow_roach/roach_ssf_results.rds"))

perch_ssf_coefs <- extract_ssf_coefs(perch_ssf_results, "Perch")
roach_ssf_coefs <- extract_ssf_coefs(roach_ssf_results, "Roach")

# Tidy treatment labels
perch_ssf_coefs <- perch_ssf_coefs %>%
  mutate(treatment = ifelse(treatment == "Control", "Control", "Exposed"))
roach_ssf_coefs <- roach_ssf_coefs %>%
  mutate(treatment = ifelse(treatment == "Control", "Control", "Exposed"))

print(perch_ssf_coefs)
print(roach_ssf_coefs)

# Save
write_csv(perch_ssf_coefs, paste0(ssf_path, "cow_perch/perch_ssf_coefs.csv"))
write_csv(roach_ssf_coefs, paste0(ssf_path, "cow_roach/roach_ssf_coefs.csv"))

# ------------------------------------------ #
# 8. PLOT SSF COEFFICIENTS ####
# ------------------------------------------ #
# Same plot style as the RSF analysis for consistency.
# Individual points + population mean (simple mean of
# individual estimates) with SE.

plot_ssf_coefs <- function(coef_df, species_label, y_limits = c(-5, 5)) {
  
  pop_summary <- coef_df %>%
    group_by(treatment) %>%
    summarise(
      mean_est = mean(est, na.rm = TRUE),
      se_est   = sd(est, na.rm = TRUE) / sqrt(n()),
      .groups  = "drop"
    ) %>%
    mutate(
      low  = mean_est - 1.96 * se_est,
      high = mean_est + 1.96 * se_est
    )
  
  ggplot() +
    geom_jitter(data = coef_df,
                aes(x = treatment, y = est, color = treatment),
                width = 0.08, size = 1.8, alpha = 0.5, shape = 16) +
    geom_errorbar(data = pop_summary,
                  aes(x = treatment, ymin = low, ymax = high),
                  width = 0.12, linewidth = 0.9, color = "black") +
    geom_point(data = pop_summary,
               aes(x = treatment, y = mean_est, fill = treatment),
               size = 4, shape = 21, color = "black") +
    scale_fill_manual(values  = c("Control" = "white", "Exposed" = "black")) +
    scale_color_manual(values = c("Control" = "grey60", "Exposed" = "grey20")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    coord_cartesian(ylim = y_limits) +
    labs(y     = "Selection coefficient (pike OD)",
         title = species_label) +
    theme_classic() +
    theme(legend.position  = "none",
          axis.title.x     = element_blank(),
          axis.title.y     = element_text(face = "bold", size = 14, margin = margin(r = 10)),
          axis.text        = element_text(size = 12, color = "black"),
          plot.title       = element_text(face = "italic", size = 13, hjust = 0.5),
          panel.border     = element_rect(color = "black", fill = NA, linewidth = 1))
}

perch_ssf_plot <- plot_ssf_coefs(perch_ssf_coefs, "Perch")
roach_ssf_plot <- plot_ssf_coefs(roach_ssf_coefs, "Roach")

print(perch_ssf_plot)
print(roach_ssf_plot)

ggsave(paste0(fig_path, "perch_ssf_pike_od_cow.png"), perch_ssf_plot,
       width = 8, height = 8, units = "cm", dpi = 300, device = "png")
ggsave(paste0(fig_path, "roach_ssf_pike_od_cow.png"), roach_ssf_plot,
       width = 8, height = 8, units = "cm", dpi = 300, device = "png")


# ---- Statistical comparison: Control vs Exposed ----
# Two-sample t-test on individual coefficients
# This is the standard "two-stage" approach in movement ecology:
# Stage 1 = fit individual models (done)
# Stage 2 = compare coefficients across groups

compare_treatments <- function(coef_df, species_label) {
  
  cat("\n=============================\n")
  cat(species_label, "— Treatment comparison\n")
  cat("=============================\n")
  
  control <- coef_df$est[coef_df$treatment == "Control"]
  exposed <- coef_df$est[coef_df$treatment == "Exposed"]
  
  cat("Control n:", length(control), "\n")
  cat("Exposed n:", length(exposed), "\n\n")
  
  # Descriptive statistics
  cat("Control  mean:", round(mean(control), 3),
      " SD:", round(sd(control), 3), "\n")
  cat("Exposed  mean:", round(mean(exposed), 3),
      " SD:", round(sd(exposed), 3), "\n\n")
  
  # Test if each group differs from zero (is there avoidance?)
  cat("--- One-sample t-tests vs zero ---\n")
  t_control <- t.test(control, mu = 0)
  t_exposed <- t.test(exposed, mu = 0)
  cat("Control vs 0:  t =", round(t_control$statistic, 3),
      " p =", round(t_control$p.value, 4), "\n")
  cat("Exposed vs 0:  t =", round(t_exposed$statistic, 3),
      " p =", round(t_exposed$p.value, 4), "\n\n")
  
  # Test if groups differ from each other
  cat("--- Two-sample t-test: Control vs Exposed ---\n")
  t_between <- t.test(exposed, control, var.equal = FALSE)  # Welch t-test
  cat("t =", round(t_between$statistic, 3),
      " df =", round(t_between$parameter, 1),
      " p =", round(t_between$p.value, 4), "\n")
  cat("Mean difference (Exposed - Control):",
      round(mean(exposed) - mean(control), 3), "\n")
  cat("95% CI of difference: [",
      round(t_between$conf.int[1], 3), ",",
      round(t_between$conf.int[2], 3), "]\n\n")
  
  # Non-parametric alternative (use if small n or non-normal)
  cat("--- Wilcoxon rank-sum test (non-parametric) ---\n")
  w_test <- wilcox.test(exposed, control, exact = FALSE)
  cat("W =", w_test$statistic,
      " p =", round(w_test$p.value, 4), "\n")
  
  # Return results as data frame
  data.frame(
    species          = species_label,
    n_control        = length(control),
    n_exposed        = length(exposed),
    mean_control     = mean(control),
    mean_exposed     = mean(exposed),
    t_control_vs0    = t_control$statistic,
    p_control_vs0    = t_control$p.value,
    t_exposed_vs0    = t_exposed$statistic,
    p_exposed_vs0    = t_exposed$p.value,
    t_between        = t_between$statistic,
    p_between        = t_between$p.value,
    mean_diff        = mean(exposed) - mean(control),
    ci_low           = t_between$conf.int[1],
    ci_high          = t_between$conf.int[2],
    p_wilcoxon       = w_test$p.value
  )
}

perch_stats <- compare_treatments(perch_ssf_coefs, "Perch")
roach_stats  <- compare_treatments(roach_ssf_coefs,  "Roach")

# Combined results table
all_stats <- dplyr::bind_rows(perch_stats, roach_stats)
print(all_stats)

write.csv(all_stats, paste0(ssf_path, "treatment_comparison_stats.csv"),
          row.names = FALSE)

# -----------------