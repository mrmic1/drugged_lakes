# ---------------------------------------------------------#
# RESOURCE SELECTION FOR ARTIFICIAL HABITATS - Cow Paradise ###
# ---------------------------------------------------------#
# Date: 2026-02-17

### LIBRARIES ###
library(ctmm)
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(data.table)
library(parallel)
library(foreach)
library(doParallel)

### DIRECTORIES ###
polygon_path      <- "./data/lake_params/polygons/"
ctmm_path         <- "./data/ctmm_fits/"
filtered_data_path <- "./data/tracks_filtered/cow_paradise/"
telem_path        <- "./data/telem_obj/cow_paradise/"
rec_data_path     <- "./data/lake_params/reciever_and_habitat_locations/"
enc_path          <- "./data/encounters/cow_paradise/"
akde_path         <- "./data/akdes/"
rsf_path          <- "./data/rsfs/habitats/"
figure_path       <- "./figures/cow_paradise/"

### LOAD DATA ###

# Receiver and habitat locations
cow_rec_locs_kml <- paste0(rec_data_path, "cow_paradise_rec_hab_locations.kml")
cow_rec_locs <- st_read(cow_rec_locs_kml)
cow_rec_locs     <- st_read(cow_rec_locs_kml)[c(1,4,6,8,9),]
cow_hab_locs     <- st_read(cow_rec_locs_kml)[c(2,3,5,7),]

# Fish telemetry data
pike_cow_tel  <- readRDS(paste0(telem_path, 'pike_cow_paradise_tel_thinned_final.rds'))
perch_cow_tel <- readRDS(paste0(telem_path, 'perch_cow_paradise_tel_thinned_final.rds'))
roach_cow_tel <- readRDS(paste0(telem_path, 'roach_cow_paradise_tel_thinned_final.rds'))

# Fish ctmm fits (freshly rerun)
pike_cow_ctmm_fits  <- readRDS(paste0(ctmm_path, "cow_paradise_pike_fits/cow_paradise_pike_best_models.rds"))
perch_cow_ctmm_fits <- readRDS(paste0(ctmm_path, "cow_paradise_perch_fits/cow_paradise_perch_best_models.rds"))
roach_cow_ctmm_fits <- readRDS(paste0(ctmm_path, "cow_paradise_roach_fits/cow_paradise_roach_best_models.rds"))

# Fish AKDEs (freshly rerun)
pike_cow_akdes  <- readRDS(paste0(akde_path, "cow_paradise_pike_akdes/pike_cow_paradise_akdes.rds"))
perch_cow_akdes <- readRDS(paste0(akde_path, "cow_paradise_perch_akdes/perch_cow_paradise_akdes.rds"))
roach_cow_akdes <- readRDS(paste0(akde_path, "cow_paradise_roach_akdes/roach_cow_paradise_akdes.rds"))

# PKDEs
perch_control_PKDE <- readRDS(paste0(akde_path, "cow_paradise_perch_akdes/perch_control_PKDE.rds"))
perch_mix_PKDE     <- readRDS(paste0(akde_path, "cow_paradise_perch_akdes/perch_mix_PKDE.rds"))
roach_control_PKDE <- readRDS(paste0(akde_path, "cow_paradise_roach_akdes/roach_control_PKDE.rds"))
roach_mix_PKDE     <- readRDS(paste0(akde_path, "cow_paradise_roach_akdes/roach_mix_PKDE.rds"))
pike_control_PKDE  <- readRDS(paste0(akde_path, "cow_paradise_pike_akdes/pike_control_PKDE.rds"))
pike_mix_PKDE      <- readRDS(paste0(akde_path, "cow_paradise_pike_akdes/pike_mix_PKDE.rds"))
pike_total_PKDE    <- readRDS(paste0(akde_path, "cow_paradise_pike_akdes/pike_total_PKDE.rds"))

# Lake polygon
cow_polygon <- st_read(paste0(polygon_path, "cow_polygon.gpkg"))

#-----------------------------------------#
# 1. Plot population AKDEs with habitats ####
#-----------------------------------------#

generate_ud_plot_whabitats <- function(pkde_data, bbox, hab_locs, title = "") {
  ud_raster        <- raster(pkde_data, DF = "CDF")
  masked_ud_raster <- mask(ud_raster, bbox)
  ud_df            <- as.data.frame(masked_ud_raster, xy = TRUE, na.rm = TRUE)
  colnames(ud_df)  <- c("x", "y", "value")
  ud_df$value      <- 1 - ud_df$value  # Invert: high = core use area
  
  ggplot() +
    geom_sf(data = bbox, color = "black") +
    geom_tile(data = ud_df, aes(x = x, y = y, fill = value), alpha = 0.6) +
    geom_sf(data = hab_locs, color = "green", size = 3, fill = NA, shape = 3, stroke = 2) +
    scale_fill_viridis_c(na.value = 'transparent', option = 'magma', direction = -1) +
    coord_sf() +
    theme_classic() +
    labs(fill = "Utilization Distribution", title = title, x = "", y = '') +
    theme(
      legend.position  = "bottom",
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      axis.text        = element_blank(),
      axis.title       = element_blank()
    )
}

# Transform bbox to match PKDE projection
cow_bbox_perch <- st_transform(cow_polygon, crs(raster(perch_control_PKDE)))
cow_bbox_roach <- st_transform(cow_polygon, crs(raster(roach_control_PKDE)))
cow_bbox_pike  <- st_transform(cow_polygon, crs(raster(pike_total_PKDE)))


# Perch plots
perch_control_plot_habitats <- generate_ud_plot_whabitats(
  perch_control_PKDE, cow_bbox_perch, cow_hab_locs, "Perch Control")
perch_exposed_plot_habitats <- generate_ud_plot_whabitats(
  perch_mix_PKDE, cow_bbox_perch, cow_hab_locs, "Perch Exposed")

# Roach plots
roach_control_plot_habitats <- generate_ud_plot_whabitats(
  roach_control_PKDE, cow_bbox_roach, cow_hab_locs, "Roach Control")
roach_exposed_plot_habitats <- generate_ud_plot_whabitats(
  roach_mix_PKDE, cow_bbox_roach, cow_hab_locs, "Roach Exposed")

#Pike plot
pike_total_plot_habitats <- generate_ud_plot_whabitats(
  pike_total_PKDE, cow_bbox_pike, cow_hab_locs, "Pike Total")


#View plots
print(perch_control_plot_habitats)
print(perch_exposed_plot_habitats)
print(roach_control_plot_habitats)
print(roach_exposed_plot_habitats)
print(pike_total_plot_habitats)


#Save plots
ggsave(paste0(figure_path, "UD_plots/perch_control_UD_habitats_cow.png"),
       perch_control_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)
ggsave(paste0(figure_path, "UD_plots/perch_exposed_UD_habitats_cow.png"),
       perch_exposed_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)

ggsave(paste0(figure_path, "UD_plots/roach_control_UD_habitats_cow.png"),
       roach_control_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)
ggsave(paste0(figure_path, "UD_plots/roach_exposed_UD_habitats_cow.png"),
       roach_exposed_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)

ggsave(paste0(figure_path, "UD_plots/pike_total_UD_habitats_cow.png"),
       pike_total_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)

#-----------------------------------------#
# 2. Create habitat raster ####
#-----------------------------------------#

# NOTE: rsf.fit() uses longitude/latitude columns from telemetry internally,
# so the raster must be in WGS84 (EPSG:4326). We build in UTM for accuracy,
# then transform to WGS84.

# Helper function to create square buffers
create_square_buffer <- function(points, size = 2) {
  squares <- lapply(1:nrow(points), function(i) {
    center <- st_coordinates(points[i,])
    bbox   <- st_bbox(c(xmin = center[1] - size/2,
                        xmax = center[1] + size/2,
                        ymin = center[2] - size/2,
                        ymax = center[2] + size/2),
                      crs = st_crs(points))
    st_as_sfc(bbox)
  })
  st_sf(geometry = do.call(c, squares), crs = st_crs(points))
}

# Transform polygon to UTM and add 2m buffer to catch edge cases
polyProj_buffered <- st_buffer(st_transform(cow_polygon, "EPSG:32633"), dist = 2)

# Create base raster in UTM
lake_raster <- rast(ext(polyProj_buffered), res = 0.25, crs = "EPSG:32633")
lake_raster <- rasterize(polyProj_buffered, lake_raster, background = NA)

# Create 2m × 2m square habitat patches in UTM
cow_hab_locs_utm      <- st_transform(cow_hab_locs, crs = "EPSG:32633")
habitat_patches      <- create_square_buffer(cow_hab_locs_utm, size = 2.0)
habitat_patches_vect <- vect(habitat_patches)

# Rasterize habitats
habitat_raster_utm <- rasterize(habitat_patches_vect, lake_raster, 
                                field = 1, background = 0)
habitat_raster_utm <- mask(habitat_raster_utm, polyProj_buffered)

# Transform to WGS84 (required by rsf.fit)
habitat_raster_wgs84 <- project(habitat_raster_utm, "EPSG:4326", method = "near")

# Force binary values after reprojection
habitat_raster_wgs84[habitat_raster_wgs84 > 0.5] <- 1
habitat_raster_wgs84[habitat_raster_wgs84 <= 0.5 & !is.na(habitat_raster_wgs84)] <- 0

# Convert to raster package format for ctmm
habitat_raster_final <- raster::raster(habitat_raster_wgs84)

plot(habitat_raster_final)

# Verify
cat("Raster extent (WGS84):\n")
print(raster::extent(habitat_raster_final))
cat("Unique values:", unique(raster::values(habitat_raster_final)), "\n")
plot(habitat_raster_final, main = "cow Habitat Raster (WGS84, 2m buffer)")

# Save
raster::writeRaster(habitat_raster_final,
                    "./data/lake_params/cow_habitat_raster_wgs84_buffered.grd",
                    overwrite = TRUE)



#-----------------------------------------#
# 3. Split telemetry and AKDEs by treatment ####
#-----------------------------------------#

#> 3.1 Perch ####
perch_treatments <- sapply(perch_cow_tel, function(tel) unique(tel$treatment))
print(perch_treatments)

perch_control_tel  <- perch_cow_tel[c(1, 21:37)]
perch_mix_tel      <- perch_cow_tel[2:20]
perch_control_akdes <- perch_cow_akdes[c(1, 21:37)]
perch_mix_akdes     <- perch_cow_akdes[2:20]

#> 3.2 Roach ####
roach_treatments <- sapply(roach_cow_tel, function(tel) unique(tel$treatment))
print(roach_treatments)

roach_control_tel  <- roach_cow_tel[1:9]
roach_mix_tel      <- roach_cow_tel[10:19]
roach_control_akdes <- roach_cow_akdes[1:9]
roach_mix_akdes     <- roach_cow_akdes[10:19]


#-----------------------#
# 4. Fit RSF models ####
#-----------------------#

# Helper function to check results and report failures
check_rsf_results <- function(rsf_list, tel_list, label) {
  failed <- sapply(rsf_list, function(x) {
    !inherits(x, "ctmm") || is.null(x) || inherits(x, "error")
  })
  cat("\n", label, "\n")
  cat("  Succeeded:", sum(!failed), "/", length(tel_list), "\n")
  if (any(failed)) {
    cat("  Failed IDs:", names(tel_list)[failed], "\n")
    for (i in which(failed)) {
      if (inherits(rsf_list[[i]], "error")) {
        cat("   ", names(tel_list)[i], ":", rsf_list[[i]]$message, "\n")
      }
    }
  }
  return(failed)
}

#> 4.1 Perch ####
cat("\n=== Fitting Perch RSFs ===\n")

# Control
cl <- makeCluster(length(perch_control_tel))
registerDoParallel(cl)
rsf_perch_control_list <- foreach(i = seq_along(perch_control_tel),
                                  .packages = c("ctmm", "raster"),
                                  .errorhandling = "pass") %dopar% {
                                    rsf_model <- rsf.fit(perch_control_tel[[i]],
                                                         perch_control_akdes[[i]],
                                                         R = list(habitat1 = habitat_raster_final))
                                    saveRDS(rsf_model, paste0(rsf_path, "cow_perch/",
                                                              names(perch_control_tel)[i], "_habitat_rsf.rds"))
                                    rsf_model
                                  }
stopCluster(cl)

failed_perch_control <- check_rsf_results(rsf_perch_control_list, perch_control_tel, "Perch Control")
rsf_perch_control_list <- rsf_perch_control_list[!failed_perch_control]
names(rsf_perch_control_list) <- names(perch_control_tel)[!failed_perch_control]
saveRDS(rsf_perch_control_list, paste0(rsf_path, "cow_perch/rsf_perch_control_list.rds"))

# Exposed
cl <- makeCluster(length(perch_mix_tel))
registerDoParallel(cl)
rsf_perch_mix_list <- foreach(i = seq_along(perch_mix_tel),
                              .packages = c("ctmm", "raster"),
                              .errorhandling = "pass") %dopar% {
                                rsf_model <- rsf.fit(perch_mix_tel[[i]],
                                                     perch_mix_akdes[[i]],
                                                     R = list(habitat1 = habitat_raster_final))
                                saveRDS(rsf_model, paste0(rsf_path, "cow_perch/",
                                                          names(perch_mix_tel)[i], "_habitat_rsf.rds"))
                                rsf_model
                              }
stopCluster(cl)

failed_perch_mix <- check_rsf_results(rsf_perch_mix_list, perch_mix_tel, "Perch Exposed")
rsf_perch_mix_list <- rsf_perch_mix_list[!failed_perch_mix]
names(rsf_perch_mix_list) <- names(perch_mix_tel)[!failed_perch_mix]
saveRDS(rsf_perch_mix_list, paste0(rsf_path, "cow_perch/rsf_perch_mix_list.rds"))

#> 4.2 Roach ####
cat("\n=== Fitting Roach RSFs ===\n")

# Control
cl <- makeCluster(length(roach_control_tel))
registerDoParallel(cl)
rsf_roach_control_list <- foreach(i = seq_along(roach_control_tel),
                                  .packages = c("ctmm", "raster"),
                                  .errorhandling = "pass") %dopar% {
                                    rsf_model <- rsf.fit(roach_control_tel[[i]],
                                                         roach_control_akdes[[i]],
                                                         R = list(habitat1 = habitat_raster_final))
                                    saveRDS(rsf_model, paste0(rsf_path, "cow_roach/",
                                                              names(roach_control_tel)[i], "_habitat_rsf.rds"))
                                    rsf_model
                                  }
stopCluster(cl)

failed_roach_control <- check_rsf_results(rsf_roach_control_list, roach_control_tel, "Roach Control")
rsf_roach_control_list <- rsf_roach_control_list[!failed_roach_control]
names(rsf_roach_control_list) <- names(roach_control_tel)[!failed_roach_control]
saveRDS(rsf_roach_control_list, paste0(rsf_path, "cow_roach/rsf_roach_control_list.rds"))

# Exposed
cl <- makeCluster(length(roach_mix_tel))
registerDoParallel(cl)
rsf_roach_mix_list <- foreach(i = seq_along(roach_mix_tel),
                              .packages = c("ctmm", "raster"),
                              .errorhandling = "pass") %dopar% {
                                rsf_model <- rsf.fit(roach_mix_tel[[i]],
                                                     roach_mix_akdes[[i]],
                                                     R = list(habitat1 = habitat_raster_final))
                                saveRDS(rsf_model, paste0(rsf_path, "cow_roach/",
                                                          names(roach_mix_tel)[i], "_habitat_rsf.rds"))
                                rsf_model
                              }
stopCluster(cl)

failed_roach_mix <- check_rsf_results(rsf_roach_mix_list, roach_mix_tel, "Roach Exposed")
rsf_roach_mix_list <- rsf_roach_mix_list[!failed_roach_mix]
names(rsf_roach_mix_list) <- names(roach_mix_tel)[!failed_roach_mix]
saveRDS(rsf_roach_mix_list, paste0(rsf_path, "cow_roach/rsf_roach_mix_list.rds"))

#> 4.3 Pike ####
cat("\n=== Fitting Pike RSFs ===\n")

# All pike together (too few for treatment split)
cl <- makeCluster(length(pike_cow_tel))
registerDoParallel(cl)
rsf_pike_list <- foreach(i = seq_along(pike_cow_tel),
                         .packages = c("ctmm", "raster"),
                         .errorhandling = "pass") %dopar% {
                           rsf_model <- rsf.fit(pike_cow_tel[[i]],
                                                pike_cow_akdes[[i]],
                                                R = list(habitat1 = habitat_raster_final))
                           saveRDS(rsf_model, paste0(rsf_path, "cow_pike/",
                                                     names(pike_cow_tel)[i], "_habitat_rsf.rds"))
                           rsf_model
                         }
stopCluster(cl)

failed_pike <- check_rsf_results(rsf_pike_list, pike_cow_tel, "Pike")
rsf_pike_list <- rsf_pike_list[!failed_pike]
names(rsf_pike_list) <- names(pike_cow_tel)[!failed_pike]
saveRDS(rsf_pike_list, paste0(rsf_path, "cow_pike/rsf_pike_list.rds"))

cat("\n=== RSF Fitting Complete ===\n")
cat("Perch control:", length(rsf_perch_control_list), "models\n")
cat("Perch exposed:", length(rsf_perch_mix_list), "models\n")
cat("Roach control:", length(rsf_roach_control_list), "models\n")
cat("Roach exposed:", length(rsf_roach_mix_list), "models\n")
cat("Pike:", length(rsf_pike_list), "models\n")

#-----------------------------------------#
# 5. Explore RSF results ####
#-----------------------------------------#

#> 5.1 Perch ####
cat("\n=== PERCH ANALYSIS ===\n")

rsf_perch_control_list <- readRDS(paste0(rsf_path, "cow_perch/rsf_perch_control_list.rds"))
rsf_perch_mix_list     <- readRDS(paste0(rsf_path, "cow_perch/rsf_perch_mix_list.rds"))

# Extract individual coefficients - CONTROL
perch_control_coefs <- data.frame(
  id = character(), est = numeric(), low = numeric(), high = numeric()
)

for(i in seq_along(rsf_perch_control_list)) {
  summ <- summary(rsf_perch_control_list[[i]])
  coef <- summ$CI["habitat1 (1/habitat1)",]
  perch_control_coefs <- rbind(perch_control_coefs, 
                               data.frame(id = names(rsf_perch_control_list)[i],
                                          est = coef["est"],
                                          low = coef["low"],
                                          high = coef["high"]))
}

# Manual population mean - CONTROL
mean_est <- mean(perch_control_coefs$est)
se_est <- sd(perch_control_coefs$est) / sqrt(nrow(perch_control_coefs))

rsf_coef_control_perch <- data.frame(
  low = mean_est - 1.96 * se_est,
  est = mean_est,
  high = mean_est + 1.96 * se_est,
  treatment = "Control"
)

# Extract individual coefficients - EXPOSED
perch_exposed_coefs <- data.frame(
  id = character(), est = numeric(), low = numeric(), high = numeric()
)

for(i in seq_along(rsf_perch_mix_list)) {
  summ <- summary(rsf_perch_mix_list[[i]])
  coef <- summ$CI["habitat1 (1/habitat1)",]
  perch_exposed_coefs <- rbind(perch_exposed_coefs, 
                               data.frame(id = names(rsf_perch_mix_list)[i],
                                          est = coef["est"],
                                          low = coef["low"],
                                          high = coef["high"]))
}

# Manual population mean - EXPOSED
mean_est_exp <- mean(perch_exposed_coefs$est)
se_est_exp <- sd(perch_exposed_coefs$est) / sqrt(nrow(perch_exposed_coefs))

rsf_coef_exposed_perch <- data.frame(
  low = mean_est_exp - 1.96 * se_est_exp,
  est = mean_est_exp,
  high = mean_est_exp + 1.96 * se_est_exp,
  treatment = "Exposed"
)

perch_habitat_rsf_coefs <- rbind(rsf_coef_control_perch, rsf_coef_exposed_perch)

# Combine individual coefficients and add treatment label
perch_control_coefs$treatment <- "Control"
perch_exposed_coefs$treatment <- "Exposed"
ind_coefs_perch <- rbind(perch_control_coefs, perch_exposed_coefs)

cat("Perch Control CI:\n"); print(rsf_coef_control_perch)
cat("Perch Exposed CI:\n"); print(rsf_coef_exposed_perch)

# Plot
perch_habitat_rsf_plot <- ggplot(perch_habitat_rsf_coefs, aes(x = treatment, y = est)) +
  geom_jitter(data = ind_coefs_perch,
              aes(x = treatment, y = est, color = treatment),
              width = 0.08, size = 1.8, alpha = 0.5, shape = 16) +
  scale_color_manual(values = c("Control" = "grey50", "Exposed" = "grey50")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, linewidth = 1, color = "black") +
  geom_point(aes(shape = treatment, fill = treatment), size = 4, color = "black") +
  scale_shape_manual(values = c(21, 21)) +
  scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +
  coord_cartesian(ylim = c(-2.5, 6)) +
  labs(y = "Selection coefficient") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),
        axis.text = element_text(size = 12, color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))

print(perch_habitat_rsf_plot)
ggsave(paste0(figure_path, "rsf_plots/perch_habitats_rsf_cow.png"),
       perch_habitat_rsf_plot, width = 8, height = 8, units = 'cm', dpi = 300)

#> 5.2 Roach ####
cat("\n=== ROACH ANALYSIS ===\n")

rsf_roach_control_list <- readRDS(paste0(rsf_path, "cow_roach/rsf_roach_control_list.rds"))
rsf_roach_mix_list     <- readRDS(paste0(rsf_path, "cow_roach/rsf_roach_mix_list.rds"))

# Extract individual coefficients - CONTROL
roach_control_coefs <- data.frame(
  id = character(), est = numeric(), low = numeric(), high = numeric()
)

for(i in seq_along(rsf_roach_control_list)) {
  summ <- summary(rsf_roach_control_list[[i]])
  coef <- summ$CI["habitat1 (1/habitat1)",]
  roach_control_coefs <- rbind(roach_control_coefs, 
                               data.frame(id = names(rsf_roach_control_list)[i],
                                          est = coef["est"],
                                          low = coef["low"],
                                          high = coef["high"]))
}

# Manual population mean - CONTROL
mean_est <- mean(roach_control_coefs$est)
se_est <- sd(roach_control_coefs$est) / sqrt(nrow(roach_control_coefs))

rsf_coef_control_roach <- data.frame(
  low = mean_est - 1.96 * se_est,
  est = mean_est,
  high = mean_est + 1.96 * se_est,
  treatment = "Control"
)

# Extract individual coefficients - EXPOSED
roach_exposed_coefs <- data.frame(
  id = character(), est = numeric(), low = numeric(), high = numeric()
)

for(i in seq_along(rsf_roach_mix_list)) {
  summ <- summary(rsf_roach_mix_list[[i]])
  coef <- summ$CI["habitat1 (1/habitat1)",]
  roach_exposed_coefs <- rbind(roach_exposed_coefs, 
                               data.frame(id = names(rsf_roach_mix_list)[i],
                                          est = coef["est"],
                                          low = coef["low"],
                                          high = coef["high"]))
}

# Manual population mean - EXPOSED
mean_est_exp <- mean(roach_exposed_coefs$est)
se_est_exp <- sd(roach_exposed_coefs$est) / sqrt(nrow(roach_exposed_coefs))

rsf_coef_exposed_roach <- data.frame(
  low = mean_est_exp - 1.96 * se_est_exp,
  est = mean_est_exp,
  high = mean_est_exp + 1.96 * se_est_exp,
  treatment = "Exposed"
)

roach_habitat_rsf_coefs <- rbind(rsf_coef_control_roach, rsf_coef_exposed_roach)

# Combine individual coefficients and add treatment label
roach_control_coefs$treatment <- "Control"
roach_exposed_coefs$treatment <- "Exposed"
ind_coefs_roach <- rbind(roach_control_coefs, roach_exposed_coefs)

cat("Roach Control CI:\n"); print(rsf_coef_control_roach)
cat("Roach Exposed CI:\n"); print(rsf_coef_exposed_roach)

# Plot
roach_habitat_rsf_plot <- ggplot(roach_habitat_rsf_coefs, aes(x = treatment, y = est)) +
  geom_jitter(data = ind_coefs_roach,
              aes(x = treatment, y = est, color = treatment),
              width = 0.08, size = 1.8, alpha = 0.5, shape = 16) +
  scale_color_manual(values = c("Control" = "grey50", "Exposed" = "grey50")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, linewidth = 1, color = "black") +
  geom_point(aes(shape = treatment, fill = treatment), size = 4, color = "black") +
  scale_shape_manual(values = c(21, 21)) +
  scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +
  coord_cartesian(ylim = c(-2, 4)) +
  labs(y = "Selection coefficient") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),
        axis.text = element_text(size = 12, color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))

print(roach_habitat_rsf_plot)
ggsave(paste0(figure_path, "roach_habitats_rsf_cow.png"),
       roach_habitat_rsf_plot, width = 8, height = 8, units = 'cm', dpi = 300)

#> 5.3 Pike ####
cat("\n=== PIKE ANALYSIS (Overall) ===\n")

rsf_pike_list <- readRDS(paste0(rsf_path, "cow_pike/rsf_pike_list.rds"))

# Extract individual coefficients - ALL PIKE
pike_all_coefs <- data.frame(
  id = character(), est = numeric(), low = numeric(), high = numeric()
)

for(i in seq_along(rsf_pike_list)) {
  summ <- summary(rsf_pike_list[[i]])
  coef <- summ$CI["habitat1 (1/habitat1)",]
  pike_all_coefs <- rbind(pike_all_coefs, 
                          data.frame(id = names(rsf_pike_list)[i],
                                     est = coef["est"],
                                     low = coef["low"],
                                     high = coef["high"]))
}

# Manual population mean - ALL PIKE
mean_est <- mean(pike_all_coefs$est)
se_est <- sd(pike_all_coefs$est) / sqrt(nrow(pike_all_coefs))

rsf_coef_pike <- data.frame(
  low = mean_est - 1.96 * se_est,
  est = mean_est,
  high = mean_est + 1.96 * se_est,
  species = "Pike"
)

cat("Pike Overall CI:\n"); print(rsf_coef_pike)
cat("Individual pike (n =", nrow(pike_all_coefs), "):\n")
print(pike_all_coefs)

# Plot (single point)
pike_habitat_rsf_plot <- ggplot(rsf_coef_pike, aes(x = species, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, linewidth = 1, color = "black") +
  geom_point(size = 4, color = "black", fill = "gray", shape = 21) +
  coord_cartesian(ylim = c(-2, 6)) +
  labs(y = "Selection coefficient") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),
        axis.text = element_text(size = 12, color = 'black'),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))

print(pike_habitat_rsf_plot)
ggsave(paste0(figure_path, "pike_habitats_rsf_cow.png"),
       pike_habitat_rsf_plot, width = 8, height = 8, units = 'cm', dpi = 300)

cat("\n=== ANALYSIS COMPLETE ===\n")
