# ============================================================ #
# PREDATOR SPACE USE AVOIDANCE - RSF & HR OVERLAP - MUDDYFOOT #
# ============================================================ #
# Goal: Determine whether prey fish (perch, roach) avoid areas
# of high pike space use, and whether pharmaceutical exposure
# (treatment) reduces this avoidance behaviour.
#
# Approach:
#   1. Use pike population PKDE as a predation risk raster
#   2. Fit individual RSFs with predator UD as covariate
#   3. Compare selection coefficients between control vs. exposed
#   4. Complement with home range overlap (prey AKDE ~ pike PKDE)
# ============================================================ #

library(ctmm)
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(parallel)
library(foreach)
library(doParallel)

# -------------------- #
# 1. DIRECTORIES ####
# -------------------- #

telem_path    <- "./data/telem_obj/muddyfoot/"
ctmm_path     <- "./data/ctmm_fits/"
akde_path     <- "./data/akdes/"
rsf_path      <- "./data/rsfs/predators/"
polygon_path  <- "./data/lake_params/polygons/"
fig_path      <- "./figures/muddyfoot/"
rec_data_path     <- "./data/lake_params/reciever_and_habitat_locations/"

# -------------------- #
# 2. LOAD DATA ####
# -------------------- #

# Telemetry objects
perch_muddyfoot_tel <- readRDS(paste0(telem_path, "perch_muddyfoot_tel_thinned_final.rds"))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, "roach_muddyfoot_tel_thinned_final.rds"))

# CTMM fits
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_best_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_best_models.rds"))

# Individual AKDEs
perch_muddyfoot_akdes <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/perch_muddyfoot_akdes.rds"))
roach_muddyfoot_akdes <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/roach_muddyfoot_akdes.rds"))

# Pike population AKDE (PKDE) - used as predation risk surface
pike_total_PKDE <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_total_PKDE.rds"))

# Lake polygon for masking (used below)
muddyfoot_polygon <- st_read(paste0(polygon_path, "muddyfoot_polygon.gpkg"))  # adjust filename as needed

# Receiver and habitat locations
mud_rec_locs_kml <- paste0(rec_data_path, "muddyfoot_rec_hab_locations.kml")
mud_rec_locs     <- st_read(mud_rec_locs_kml)[1:5,]
mud_hab_locs     <- st_read(mud_rec_locs_kml)[6:7,]


# ------------------------------------------ #
# 3. BUILD PREDATOR RISK RASTER ####
# ------------------------------------------ #
# Convert the pike PKDE to a predation risk raster.
# The PKDE default output is the CDF, so we explicitly request
# DF="PDF" so that high values = high pike use (core areas).
# The raster remains in the native tpeqd CRS which matches
# the telemetry objects - no reprojection needed.

# Step 1: Extract PDF raster directly from PKDE
predator_ud_raster <- raster::raster(pike_total_PKDE, DF = "PDF")
plot(predator_ud_raster, main = "Pike UD - PDF (high = high pike use)")

# Step 2: Mask raster to lake polygon to exclude areas outside lake
predator_ud_terra <- terra::rast(predator_ud_raster)
muddyfoot_polygon_tpeqd <- st_transform(muddyfoot_polygon,
                                        crs = st_crs(predator_ud_raster))
predator_ud_terra_masked <- terra::mask(predator_ud_terra,
                                        terra::vect(muddyfoot_polygon_tpeqd))

# Step 3: Convert back to {raster} for ctmm::rsf.fit() compatibility
predator_ud_raster <- raster::raster(predator_ud_terra_masked)
plot(predator_ud_raster, main = "Pike UD - masked")

# Step 4: Build data frame for ggplot
predator_ud_df <- as.data.frame(predator_ud_terra_masked, xy = TRUE) %>%
  rename(pred_ud = 3)

# Step 5: Inspect value distribution
cat("Non-NA cells:", sum(!is.na(predator_ud_df$pred_ud)), "\n")
cat("Value range:", range(predator_ud_df$pred_ud, na.rm = TRUE), "\n")
print(quantile(predator_ud_df$pred_ud,
               probs = c(0.5, 0.75, 0.9, 0.95, 0.99),
               na.rm = TRUE))

# Transform habitat locations to match raster CRS
mud_hab_locs_tpeqd <- st_transform(mud_hab_locs, crs = st_crs(predator_ud_raster))
mud_rec_locs_tpeqd <- st_transform(mud_rec_locs, crs = st_crs(predator_ud_raster))

# Step 6: Plot predation risk surface with habitat locations
pred_value_plot <- ggplot(predator_ud_df, aes(x = x, y = y)) +
  geom_raster(aes(fill = pred_ud)) +
  scale_fill_gradientn(
    name    = "Predation\nrisk",
    colours = c("#2d004b", "#3b0f70", "#8c2981", "#de4968", "#fe9f6d", "#fcfdbf"),
    limits  = c(0, quantile(predator_ud_df$pred_ud, 0.995, na.rm = TRUE)),
    oob     = scales::squish,
    na.value = NA
  ) +
  geom_sf(data = muddyfoot_polygon_tpeqd, fill = NA,
          color = "black", linewidth = 0.5,
          inherit.aes = FALSE) +
  geom_sf(data = mud_hab_locs_tpeqd, 
          color = "green", size = 3, shape = 3, stroke = 2,
          inherit.aes = FALSE) +
  coord_sf() +
  labs(title = "Predation Risk Surface",
       x = "X (m)", y = "Y (m)") +
  theme_classic() +
  theme(axis.text    = element_text(size = 9, color = "black"),
        axis.title   = element_text(face = "bold", size = 11),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

print(pred_value_plot)

# Step 7: Save raster and plot
writeRaster(predator_ud_raster,
            paste0(rsf_path, "muddyfoot_predator_ud_raster.tif"),
            overwrite = TRUE)

ggsave(filename = paste0(fig_path, "predator_ud_raster.png"),
       plot     = pred_value_plot,
       device   = "png",
       width    = 12, height = 7, units = "cm", dpi = 300)

# ------------------------------------------ #
# 4. SPLIT TELEMETRY & AKDEs BY TREATMENT ####
# ------------------------------------------ #
# Use treatment labels stored in telemetry objects rather than
# hard-coded indices, which are fragile if list order changes.

# Perch
perch_treatments <- sapply(perch_muddyfoot_tel, function(tel) unique(tel$treatment))
print(perch_treatments)

perch_control_tel   <- perch_muddyfoot_tel[perch_treatments == "Control"]
perch_exposed_tel   <- perch_muddyfoot_tel[perch_treatments == "Mix"]
perch_control_akdes <- perch_muddyfoot_akdes[perch_treatments == "Control"]
perch_exposed_akdes <- perch_muddyfoot_akdes[perch_treatments == "Mix"]

# Roach
roach_treatments <- sapply(roach_muddyfoot_tel, function(tel) unique(tel$treatment))
print(roach_treatments)

roach_control_tel   <- roach_muddyfoot_tel[roach_treatments == "Control"]
roach_exposed_tel   <- roach_muddyfoot_tel[roach_treatments == "Mix"]
roach_control_akdes <- roach_muddyfoot_akdes[roach_treatments == "Control"]
roach_exposed_akdes <- roach_muddyfoot_akdes[roach_treatments == "Mix"]

# ------------------------------------------ #
# 5. FIT RSF MODELS ####
# ------------------------------------------ #
# Fit individual RSFs with the predator UD raster as the single
# covariate. A negative coefficient = avoidance of pike-use areas.
# Models are saved individually AND as a named list.

# Step 1: Verify inputs before fitting
cat("Number of telemetry objects:", length(perch_control_tel), "\n")
cat("Number of AKDE objects:", length(perch_control_akdes), "\n")
cat("Individual IDs:\n")
print(names(perch_control_tel))

# Step 2: Check raster values at observed locations for each individual
cat("\n--- Raster value check ---\n")
for (i in seq_along(perch_control_tel)) {
  id     <- names(perch_control_tel)[i]
  coords <- data.frame(x = perch_control_tel[[i]]$x,
                       y = perch_control_tel[[i]]$y)
  vals   <- raster::extract(predator_ud_raster, coords)
  cat(sprintf("%-12s  n = %4d  NAs = %3d  mean = %.6f  sd = %.6f  range = [%.6f, %.6f]\n",
              id, nrow(coords), sum(is.na(vals)),
              mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
              min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
}

# Step 3: Visual check - plot first individual's locations on raster
test_id     <- names(perch_control_tel)[1]
test_coords <- data.frame(x = perch_control_tel[[test_id]]$x,
                          y = perch_control_tel[[test_id]]$y)
plot(predator_ud_raster, main = paste("Predator UD —", test_id))
points(test_coords$x, test_coords$y, pch = 16, cex = 0.3, col = "red")

# Step 4: Fit RSF for each individual sequentially
# (sequential so we can see exactly which one fails and why)
rsf_perch_control_list <- list()

for (i in seq_along(perch_control_tel)) {
  
  id <- names(perch_control_tel)[i]
  cat("\nFitting RSF for:", id, "(", i, "of", length(perch_control_tel), ")\n")
  
  rsf_perch_control_list[[id]] <- tryCatch({
    
    model <- rsf.fit(
      perch_control_tel[[i]],
      perch_control_akdes[[i]],
      R = list(predator_ud = predator_ud_raster)
    )
    
    cat("  SUCCESS — est =", round(summary(model)$CI[1, "est"], 4),
        "  CI = [", round(summary(model)$CI[1, "low"], 4), ",",
        round(summary(model)$CI[1, "high"], 4), "]\n")
    
    saveRDS(model, paste0(rsf_path, "muddyfoot_perch/", id, "_predator_rsf.rds"))
    
    model
    
  }, error = function(e) {
    cat("  FAILED —", e$message, "\n")
    return(NULL)
  })
}

# Step 5: Summary of what worked
cat("\n--- Fitting summary ---\n")
succeeded <- names(rsf_perch_control_list)[!sapply(rsf_perch_control_list, is.null)]
failed    <- names(rsf_perch_control_list)[ sapply(rsf_perch_control_list, is.null)]
cat("Succeeded:", length(succeeded), "—", paste(succeeded, collapse = ", "), "\n")
cat("Failed:   ", length(failed),    "—", paste(failed,    collapse = ", "), "\n")


# ---- Verify RSF models ----

verify_rsf_list <- function(rsf_list, group_label) {
  
  cat("\n=============================\n")
  cat("Verification:", group_label, "\n")
  cat("=============================\n")
  
  # Total number of models
  cat("Total individuals:", length(rsf_list), "\n")
  
  # Check which are NULL (failed completely)
  is_null <- sapply(rsf_list, is.null)
  if (any(is_null)) {
    cat("FAILED (NULL) models:", sum(is_null), "\n")
    cat("  Failed IDs:", paste(names(rsf_list)[is_null], collapse = ", "), "\n")
  } else {
    cat("No NULL models - all individuals returned a result\n")
  }
  
  # Check which returned an error object instead of a model
  is_error <- sapply(rsf_list, function(x) inherits(x, "error") | inherits(x, "try-error"))
  if (any(is_error)) {
    cat("ERROR models:", sum(is_error), "\n")
    cat("  Error IDs:", paste(names(rsf_list)[is_error], collapse = ", "), "\n")
  }
  
  # For models that exist, try to extract a summary
  valid_ids   <- names(rsf_list)[!is_null & !is_error]
  failed_summ <- c()
  
  for (id in valid_ids) {
    ci <- tryCatch(
      summary(rsf_list[[id]])$CI,
      error = function(e) NULL
    )
    if (is.null(ci)) {
      failed_summ <- c(failed_summ, id)
    }
  }
  
  if (length(failed_summ) > 0) {
    cat("Models that exist but summary FAILED:", length(failed_summ), "\n")
    cat("  IDs:", paste(failed_summ, collapse = ", "), "\n")
  } else {
    cat("All valid models passed summary extraction\n")
  }
  
  # Print coefficient estimates for all valid models
  cat("\nIndividual coefficients:\n")
  for (id in valid_ids) {
    ci <- tryCatch(summary(rsf_list[[id]])$CI, error = function(e) NULL)
    if (!is.null(ci)) {
      cat(sprintf("  %-12s  est = %7.4f  [%7.4f, %7.4f]\n",
                  id, ci[1, "est"], ci[1, "low"], ci[1, "high"]))
    }
  }
  
  cat("\nSuccessful models:", length(valid_ids) - length(failed_summ), 
      "/", length(rsf_list), "\n")
}


# >>> 5.1 Perch ####
rsf_perch_control_list <- fit_rsf_parallel(
  tel_list    = perch_control_tel,
  akde_list   = perch_control_akdes,
  risk_raster = predator_ud_raster,
  save_dir    = paste0(rsf_path, "muddyfoot_perch/")
)

verify_rsf_list(rsf_perch_control_list, "Perch - Control")

saveRDS(rsf_perch_control_list, paste0(rsf_path, "muddyfoot_perch/pred_rsf_perch_control_list.rds"))




rsf_perch_exposed_list <- fit_rsf_parallel(
  tel_list    = perch_exposed_tel,
  akde_list   = perch_exposed_akdes,
  risk_raster = predator_ud_raster,
  save_dir    = paste0(rsf_path, "muddyfoot_perch/")
)
saveRDS(rsf_perch_exposed_list, paste0(rsf_path, "muddyfoot_perch/pred_rsf_perch_exposed_list.rds"))

# >>> 5.2 Roach ####
rsf_roach_control_list <- fit_rsf_parallel(
  tel_list    = roach_control_tel,
  akde_list   = roach_control_akdes,
  risk_raster = predator_ud_raster,
  save_dir    = paste0(rsf_path, "muddyfoot_roach/")
)
saveRDS(rsf_roach_control_list, paste0(rsf_path, "muddyfoot_roach/pred_rsf_roach_control_list.rds"))

rsf_roach_exposed_list <- fit_rsf_parallel(
  tel_list    = roach_exposed_tel,
  akde_list   = roach_exposed_akdes,
  risk_raster = predator_ud_raster,
  save_dir    = paste0(rsf_path, "muddyfoot_roach/")
)
saveRDS(rsf_roach_exposed_list, paste0(rsf_path, "muddyfoot_roach/pred_rsf_roach_exposed_list.rds"))

# ------------------------------------------ #
# 6. EXTRACT RSF COEFFICIENTS ####
# ------------------------------------------ #
# Extract individual-level coefficients first (important for
# checking distributions and identifying outliers transparently),
# then summarise at population level using ctmm::mean().

extract_rsf_coefs <- function(rsf_list, treatment_label, species_label) {
  
  # Individual-level coefficients
  ind_coefs <- lapply(names(rsf_list), function(id) {
    ci <- tryCatch(summary(rsf_list[[id]])$CI, error = function(e) NULL)
    if (is.null(ci)) return(NULL)
    data.frame(
      id        = id,
      species   = species_label,
      treatment = treatment_label,
      low       = ci[1, "low"],
      est       = ci[1, "est"],
      high      = ci[1, "high"]
    )
  })
  
  bind_rows(Filter(Negate(is.null), ind_coefs))
}

# Load if re-running from saved files (comment out if running in sequence)
# rsf_perch_control_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/pred_rsf_perch_control_list.rds"))
# rsf_perch_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/pred_rsf_perch_exposed_list.rds"))
# rsf_roach_control_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/pred_rsf_roach_control_list.rds"))
# rsf_roach_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/pred_rsf_roach_exposed_list.rds"))

# Extract individual-level coefficients
perch_rsf_coefs_ind <- bind_rows(
  extract_rsf_coefs(rsf_perch_control_list, "Control", "Perch"),
  extract_rsf_coefs(rsf_perch_exposed_list, "Exposed", "Perch")
)

roach_rsf_coefs_ind <- bind_rows(
  extract_rsf_coefs(rsf_roach_control_list, "Control", "Roach"),
  extract_rsf_coefs(rsf_roach_exposed_list, "Exposed", "Roach")
)

# Inspect for outliers before computing population means
print(perch_rsf_coefs_ind)
print(roach_rsf_coefs_ind)

# Extract population-mean coefficients using ctmm::mean()
# (this is preferable to manually averaging; ctmm propagates uncertainty correctly)
extract_pop_mean_coef <- function(rsf_list, treatment_label, species_label) {
  pop_mean <- mean(rsf_list)
  ci <- summary(pop_mean)$CI
  data.frame(
    species   = species_label,
    treatment = treatment_label,
    low       = ci[1, "low"],
    est       = ci[1, "est"],
    high      = ci[1, "high"]
  )
}

perch_rsf_coefs_pop <- bind_rows(
  extract_pop_mean_coef(rsf_perch_control_list, "Control", "Perch"),
  extract_pop_mean_coef(rsf_perch_exposed_list, "Exposed", "Perch")
)

roach_rsf_coefs_pop <- bind_rows(
  extract_pop_mean_coef(rsf_roach_control_list, "Control", "Roach"),
  extract_pop_mean_coef(rsf_roach_exposed_list, "Exposed", "Roach")
)

# ------------------------------------------ #
# 7. PLOT RSF COEFFICIENTS ####
# ------------------------------------------ #
# Plot both individual-level (jittered points) and population-mean
# (larger point + CI) on the same panel for transparency.

plot_rsf_coefs <- function(ind_coefs, pop_coefs, species_label, y_limits = c(-3, 3)) {
  
  ggplot() +
    # Individual-level estimates (jittered for visibility)
    geom_jitter(data = ind_coefs,
                aes(x = treatment, y = est, color = treatment),
                width = 0.08, size = 1.8, alpha = 0.5, shape = 16) +
    # Population-mean CI
    geom_errorbar(data = pop_coefs,
                  aes(x = treatment, ymin = low, ymax = high),
                  width = 0.12, linewidth = 0.9, color = "black") +
    # Population-mean point
    geom_point(data = pop_coefs,
               aes(x = treatment, y = est, fill = treatment),
               size = 4, shape = 21, color = "black") +
    scale_fill_manual(values  = c("Control" = "white", "Exposed" = "black")) +
    scale_color_manual(values = c("Control" = "grey60", "Exposed" = "grey20")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    coord_cartesian(ylim = y_limits) +
    labs(y     = "Selection coefficient",
         title = species_label) +
    theme_classic() +
    theme(legend.position  = "none",
          axis.title.x     = element_blank(),
          axis.title.y     = element_text(face = "bold", size = 14, margin = margin(r = 10)),
          axis.text        = element_text(size = 12, color = "black"),
          plot.title       = element_text(face = "italic", size = 13, hjust = 0.5),
          panel.border     = element_rect(color = "black", fill = NA, linewidth = 1))
}

perch_rsf_plot <- plot_rsf_coefs(perch_rsf_coefs_ind, perch_rsf_coefs_pop, "Perch")
roach_rsf_plot <- plot_rsf_coefs(roach_rsf_coefs_ind, roach_rsf_coefs_pop, "Roach")

print(perch_rsf_plot)
print(roach_rsf_plot)

ggsave(paste0(fig_path, "perch_predator_rsf_muddyfoot.png"), perch_rsf_plot,
       width = 8, height = 8, units = "cm", dpi = 300, device = "png")

ggsave(paste0(fig_path, "roach_predator_rsf_muddyfoot.png"), roach_rsf_plot,
       width = 8, height = 8, units = "cm", dpi = 300, device = "png")

# ------------------------------------------ #
# 8. HOME RANGE OVERLAP WITH PIKE PKDE ####
# ------------------------------------------ #
# Compute Bhattacharyya Affinity (BA) overlap between each prey
# individual's AKDE and the pike population PKDE.
# Higher BA = more spatial overlap with pike space use.
# We expect: Exposed > Control if pharma reduces avoidance.
#
# ctmm::overlap() returns BA with CIs. We call it pairwise:
# overlap(list(prey_akde, pike_pkde)) and extract the [1,2] element.

compute_overlap_with_pkde <- function(akde_list, pkde, treatment_label, species_label) {
  
  results <- lapply(names(akde_list), function(id) {
    
    ov <- tryCatch(
      ctmm::overlap(list(akde_list[[id]], pkde))$CI,
      error = function(e) {
        message("Overlap failed for ", id, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(ov)) return(NULL)
    
    # overlap() returns a 2x2x3 array; [1,2,] gives prey~pike overlap
    data.frame(
      id        = id,
      species   = species_label,
      treatment = treatment_label,
      low       = ov[1, 2, "low"],
      est       = ov[1, 2, "est"],
      high      = ov[1, 2, "high"]
    )
  })
  
  bind_rows(Filter(Negate(is.null), results))
}

# Compute for all groups (sequential — overlap() is fast enough)
overlap_perch <- bind_rows(
  compute_overlap_with_pkde(perch_control_akdes, pike_total_PKDE, "Control", "Perch"),
  compute_overlap_with_pkde(perch_exposed_akdes, pike_total_PKDE, "Exposed", "Perch")
)

overlap_roach <- bind_rows(
  compute_overlap_with_pkde(roach_control_akdes, pike_total_PKDE, "Control", "Roach"),
  compute_overlap_with_pkde(roach_exposed_akdes, pike_total_PKDE, "Exposed", "Roach")
)

# Save overlap tables
write_csv(overlap_perch, paste0(rsf_path, "muddyfoot_perch/perch_pkde_overlap.csv"))
write_csv(overlap_roach, paste0(rsf_path, "muddyfoot_roach/roach_pkde_overlap.csv"))

print(overlap_perch)
print(overlap_roach)

# -----------------