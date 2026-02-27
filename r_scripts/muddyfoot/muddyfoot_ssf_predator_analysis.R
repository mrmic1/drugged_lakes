# ================================================================ #
# DYNAMIC SSF PREDATOR AVOIDANCE ANALYSIS - MUDDYFOOT             #
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

telem_path   <- "./data/telem_obj/muddyfoot/"
akde_path    <- "./data/akdes/"
ssf_path     <- "./data/ssfs/predators/"
fig_path     <- "./figures/muddyfoot/"

# Create output directories
dir.create(paste0(ssf_path, "muddyfoot_perch/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(ssf_path, "muddyfoot_roach/"), recursive = TRUE, showWarnings = FALSE)

# -------------------- #
# 2. LOAD DATA ####
# -------------------- #

# Telemetry objects (ctmm format)
pike_muddyfoot_tel   <- readRDS(paste0(telem_path, "pike_muddyfoot_tel_thinned_final.rds"))
perch_muddyfoot_tel  <- readRDS(paste0(telem_path, "perch_muddyfoot_tel_thinned_final.rds"))
roach_muddyfoot_tel  <- readRDS(paste0(telem_path, "roach_muddyfoot_tel_thinned_final.rds"))

# ------------------------------------------ #
# 3. CONVERT TELEMETRY TO AMT TRACKS ####
# ------------------------------------------ #
# Following the same approach as previous RSF analysis:
# convert WGS84 coordinates to UTM 33N before making tracks.
# ctmm telemetry objects store raw coordinates as longitude/latitude
# in the 'longitude' and 'latitude' columns (not x/y which are tpeqd)

# Check column names available in telemetry objects
cat("Pike telemetry columns:", names(pike_muddyfoot_tel[[1]]), "\n")

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
  lapply(names(pike_muddyfoot_tel), function(id) {
    ctmm_tel_to_utm(pike_muddyfoot_tel[[id]], id_label = id)
  })
)

# Convert perch (with treatment column)
perch_df <- bind_rows(
  lapply(names(perch_muddyfoot_tel), function(id) {
    tel <- perch_muddyfoot_tel[[id]]
    ctmm_tel_to_utm(tel, id_label = id,
                    extra_cols = list(treatment = unique(tel$treatment)))
  })
)

# Convert roach (with treatment column)
roach_df <- bind_rows(
  lapply(names(roach_muddyfoot_tel), function(id) {
    tel <- roach_muddyfoot_tel[[id]]
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

# Verify median fix interval for pike
median_interval_mins <- median(diff(as.numeric(pike_muddyfoot_tel[[1]]$timestamp))) / 60
cat("Median fix interval:", round(median_interval_mins, 2), "minutes\n")

# With ~30 second (0.5 min) fixes:
# 120 points = 60 minutes  (1 hour)
# 60  points = 30 minutes
# 240 points = 120 minutes (2 hours)
#
# For a small lake (~40 x 30 m), pike can traverse the entire
# lake in minutes, so a 1-hour window is ecologically reasonable
# as a timescale over which prey could perceive cumulative pike
# space use through direct encounters or chemical/visual cues.
# This matches the 4-hour window used in Schlägel et al. (2019)
# scaled to the much finer fix rate here.

n_points_window <- 120  # ~1 hour at 30-second fix rate
cat("Window size:", n_points_window, "points =",
    round(n_points_window * median_interval_mins, 0), "minutes\n\n")

pike_ids <- unique(pike_df$id)
cat("Computing rolling ODs for", length(pike_ids), "pike individuals\n")

pike_rolling_ods <- lapply(pike_ids, function(pid) {
  cat("Computing rolling OD for pike:", pid, "\n")
  
  # Check column names first
  cat("  Columns in pike_df:", names(pike_df), "\n")
  
  # Select using actual column names from pike_df
  pike_sub_df <- pike_df %>%
    dplyr::filter(id == pid) %>%
    dplyr::select(x_, y_, t)  # t not t_
  
  # Rename t to t_ as required by amt track format
  names(pike_sub_df)[names(pike_sub_df) == "t"] <- "t_"
  
  # Rebuild as fresh minimal track
  pike_sub_track <- structure(
    pike_sub_df,
    class = c("track_xyt", "track_xy", "tbl_df", "tbl", "data.frame"),
    crs_  = sf::st_crs(32633)
  )
  
  cat("  Locations:", nrow(pike_sub_track), "\n")
  cat("  Columns:", names(pike_sub_track), "\n")
  
  tryCatch(
    rolling_od(pike_sub_track, trast_terra,
               n.points      = n_points_window,
               show.progress = TRUE),
    error = function(e) {
      message("  rolling_od FAILED: ", e$message)
      return(NULL)
    }
  )
})

names(pike_rolling_ods) <- pike_ids

names(pike_rolling_ods) <- pike_ids

# Summary
cat("\n--- Rolling OD summary ---\n")
for (pid in pike_ids) {
  od <- pike_rolling_ods[[pid]]
  if (is.null(od)) {
    cat(sprintf("  %-10s  FAILED\n", pid))
  } else {
    cat(sprintf("  %-10s  SUCCESS — %d layers\n", pid, raster::nlayers(od)))
  }
}

pike_rolling_ods <- Filter(Negate(is.null), pike_rolling_ods)
cat("\nSuccessful ODs:", length(pike_rolling_ods), "/", length(pike_ids), "\n")

saveRDS(pike_rolling_ods, paste0(ssf_path, "pike_rolling_ods.rds"))
# pike_rolling_ods <- readRDS(paste0(ssf_path, "pike_rolling_ods.rds"))

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

fit_prey_ssf <- function(prey_track_sub, pike_od_list, max_time_mins = 30) {
  
  # Average pike ODs across individuals to get a single
  # population-level pike occurrence surface per time step
  # (simpler than fitting one covariate per pike individual,
  # and more appropriate for a predation risk question)
  if (length(pike_od_list) > 1) {
    # Stack and average across pike ODs
    # (raster stacks must have same dimensions - check trast alignment)
    pike_od_mean <- Reduce("+", pike_od_list) / length(pike_od_list)
  } else {
    pike_od_mean <- pike_od_list[[1]]
  }
  
  tryCatch({
    prey_track_sub %>%
      steps() %>%
      random_steps(n_control = 10) %>%
      extract_covariates_var_time(
        pike_od_mean,
        max_time  = minutes(max_time_mins),  # max time lag for covariate matching
        when      = "end",                   # extract at step endpoint
        name_covar = "pike_od"
      ) %>%
      filter(!is.na(pike_od)) %>%
      fit_ssf(case_ ~ pike_od + strata(step_id_))
  },
  error = function(e) {
    message("SSF failed: ", e$message)
    return(NULL)
  })
}

# ---- 6.1 Perch ####

perch_ids       <- unique(perch_track$id)
perch_treatments <- perch_df %>% distinct(id, treatment)

perch_ssf_results <- lapply(perch_ids, function(pid) {
  
  cat("Fitting SSF for perch:", pid, "\n")
  
  prey_sub  <- filter(perch_track, id == pid)
  treatment <- perch_treatments$treatment[perch_treatments$id == pid]
  
  model <- fit_prey_ssf(prey_sub, pike_rolling_ods)
  
  if (!is.null(model)) {
    saveRDS(model, paste0(ssf_path, "muddyfoot_perch/", pid, "_ssf.rds"))
  }
  
  list(id = pid, treatment = treatment, model = model)
})

names(perch_ssf_results) <- perch_ids
saveRDS(perch_ssf_results, paste0(ssf_path, "muddyfoot_perch/perch_ssf_results.rds"))

# ---- 6.2 Roach ####

roach_ids        <- unique(roach_track$id)
roach_treatments <- roach_df %>% distinct(id, treatment)

roach_ssf_results <- lapply(roach_ids, function(rid) {
  
  cat("Fitting SSF for roach:", rid, "\n")
  
  prey_sub  <- filter(roach_track, id == rid)
  treatment <- roach_treatments$treatment[roach_treatments$id == rid]
  
  model <- fit_prey_ssf(prey_sub, pike_rolling_ods)
  
  if (!is.null(model)) {
    saveRDS(model, paste0(ssf_path, "muddyfoot_roach/", rid, "_ssf.rds"))
  }
  
  list(id = rid, treatment = treatment, model = model)
})

names(roach_ssf_results) <- roach_ids
saveRDS(roach_ssf_results, paste0(ssf_path, "muddyfoot_roach/roach_ssf_results.rds"))

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
# perch_ssf_results <- readRDS(paste0(ssf_path, "muddyfoot_perch/perch_ssf_results.rds"))
# roach_ssf_results <- readRDS(paste0(ssf_path, "muddyfoot_roach/roach_ssf_results.rds"))

perch_ssf_coefs <- extract_ssf_coefs(perch_ssf_results, "Perch")
roach_ssf_coefs <- extract_ssf_coefs(roach_ssf_results, "Roach")

# Tidy treatment labels
perch_ssf_coefs <- perch_ssf_coefs %>%
  mutate(treatment = ifelse(treatment == "control", "Control", "Exposed"))
roach_ssf_coefs <- roach_ssf_coefs %>%
  mutate(treatment = ifelse(treatment == "control", "Control", "Exposed"))

print(perch_ssf_coefs)
print(roach_ssf_coefs)

# Save
write_csv(perch_ssf_coefs, paste0(ssf_path, "muddyfoot_perch/perch_ssf_coefs.csv"))
write_csv(roach_ssf_coefs, paste0(ssf_path, "muddyfoot_roach/roach_ssf_coefs.csv"))

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

ggsave(paste0(fig_path, "perch_ssf_pike_od_muddyfoot.png"), perch_ssf_plot,
       width = 8, height = 8, units = "cm", dpi = 300, device = "png")
ggsave(paste0(fig_path, "roach_ssf_pike_od_muddyfoot.png"), roach_ssf_plot,
       width = 8, height = 8, units = "cm", dpi = 300, device = "png")
