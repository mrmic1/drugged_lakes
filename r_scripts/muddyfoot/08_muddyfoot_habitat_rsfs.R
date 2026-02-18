# ---------------------------------------------------------#
# RESOURCE SELECTION FOR ARTIFICIAL HABITATS - MUDDYFOOT ###
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
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
telem_path        <- "./data/telem_obj/muddyfoot/"
rec_data_path     <- "./data/lake_params/reciever_and_habitat_locations/"
enc_path          <- "./data/encounters/muddyfoot/"
akde_path         <- "./data/akdes/"
rsf_path          <- "./data/rsfs/habitats/"
figure_path       <- "./figures/muddyfoot/"

### LOAD DATA ###

# Receiver and habitat locations
mud_rec_locs_kml <- paste0(rec_data_path, "muddyfoot_rec_hab_locations.kml")
mud_rec_locs     <- st_read(mud_rec_locs_kml)[1:5,]
mud_hab_locs     <- st_read(mud_rec_locs_kml)[6:7,]

# Fish telemetry data
pike_muddyfoot_tel  <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel_thinned_final.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel_thinned_final.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel_thinned_final.rds'))

# Fish ctmm fits (freshly rerun)
pike_muddyfoot_ctmm_fits  <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_best_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_best_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_best_models.rds"))

# Fish AKDEs (freshly rerun)
pike_muddyfoot_akdes  <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_muddyfoot_akdes.rds"))
perch_muddyfoot_akdes <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/perch_muddyfoot_akdes.rds"))
roach_muddyfoot_akdes <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/roach_muddyfoot_akdes.rds"))

# PKDEs
perch_control_PKDE <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/perch_control_PKDE.rds"))
perch_mix_PKDE     <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/perch_mix_PKDE.rds"))
roach_control_PKDE <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/roach_control_PKDE.rds"))
roach_mix_PKDE     <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/roach_mix_PKDE.rds"))
pike_control_PKDE  <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_control_PKDE.rds"))
pike_mix_PKDE      <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_mix_PKDE.rds"))
pike_total_PKDE    <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_total_PKDE.rds"))

# Lake polygon
muddyfoot_polygon <- st_read(paste0(polygon_path, "muddyfoot_polygon.gpkg"))

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
muddyfoot_bbox_perch <- st_transform(muddyfoot_polygon, crs(raster(perch_control_PKDE)))
muddyfoot_bbox_roach <- st_transform(muddyfoot_polygon, crs(raster(roach_control_PKDE)))
muddyfoot_bbox_pike  <- st_transform(muddyfoot_polygon, crs(raster(pike_total_PKDE)))

#> 1.1 Perch ####
perch_control_plot_habitats <- generate_ud_plot_whabitats(
  perch_control_PKDE, muddyfoot_bbox_perch, mud_hab_locs, "Perch Control")
perch_exposed_plot_habitats <- generate_ud_plot_whabitats(
  perch_mix_PKDE, muddyfoot_bbox_perch, mud_hab_locs, "Perch Exposed")

ggsave(paste0(figure_path, "UD_plots/perch_control_UD_habitats_muddyfoot.png"),
       perch_control_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)
ggsave(paste0(figure_path, "UD_plots/perch_exposed_UD_habitats_muddyfoot.png"),
       perch_exposed_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)

#> 1.2 Roach ####
roach_control_plot_habitats <- generate_ud_plot_whabitats(
  roach_control_PKDE, muddyfoot_bbox_roach, mud_hab_locs, "Roach Control")
roach_exposed_plot_habitats <- generate_ud_plot_whabitats(
  roach_mix_PKDE, muddyfoot_bbox_roach, mud_hab_locs, "Roach Exposed")

ggsave(paste0(figure_path, "UD_plots/roach_control_UD_habitats_muddyfoot.png"),
       roach_control_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)
ggsave(paste0(figure_path, "UD_plots/roach_exposed_UD_habitats_muddyfoot.png"),
       roach_exposed_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)

#> 1.3 Pike ####
pike_total_plot_habitats <- generate_ud_plot_whabitats(
  pike_total_PKDE, muddyfoot_bbox_pike, mud_hab_locs, "Pike Total")

ggsave(paste0(figure_path, "UD_plots/pike_total_UD_habitats_muddyfoot.png"),
       pike_total_plot_habitats, width = 9, height = 6.5, units = 'cm', dpi = 300)

#-----------------------------------------#
# 2. Create habitat raster ####
#-----------------------------------------#

# NOTE: rsf.fit() uses longitude/latitude columns from telemetry internally,
# so the raster must be in WGS84 (EPSG:4326)

# Step 1: Create raster in UTM for accurate square geometry
polyProj <- st_transform(muddyfoot_polygon, crs = "EPSG:32633")

lake_raster <- rast(ext(polyProj), res = 0.25, crs = "EPSG:32633")
lake_raster <- rasterize(polyProj, lake_raster, background = NA)

# Step 2: Create 2m x 2m square habitat patches in UTM
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

mud_hab_locs_utm      <- st_transform(mud_hab_locs, crs = "EPSG:32633")
habitat_patches_buffered <- create_square_buffer(mud_hab_locs_utm, size = 2.0)
habitat_patches_vect  <- vect(habitat_patches_buffered)

habitat_raster_utm <- rasterize(habitat_patches_vect, lake_raster, field = 1, background = 0)
habitat_raster_utm <- mask(habitat_raster_utm, polyProj)

# Step 3: Transform to WGS84 (required by rsf.fit which uses lon/lat internally)
habitat_raster_wgs84_terra <- project(habitat_raster_utm, "EPSG:4326", method = "near")

# Force binary values after reprojection
habitat_raster_wgs84_terra[habitat_raster_wgs84_terra > 0.5]  <- 1
habitat_raster_wgs84_terra[habitat_raster_wgs84_terra <= 0.5 & !is.na(habitat_raster_wgs84_terra)] <- 0

# Convert to raster package format for ctmm
habitat_raster_final <- raster::raster(habitat_raster_wgs84_terra)

cat("Raster extent (WGS84):\n")
print(raster::extent(habitat_raster_final))
cat("Unique values:", unique(raster::values(habitat_raster_final)), "\n")

# Save
raster::writeRaster(habitat_raster_final,
                    "./data/lake_coords/muddyfoot_habitat_raster_wgs84.grd",
                    overwrite = TRUE)

#-----------------------------------------#
# 3. Split telemetry and AKDEs by treatment ####
#-----------------------------------------#

#> 3.1 Perch ####
perch_treatments <- sapply(perch_muddyfoot_tel, function(tel) unique(tel$treatment))
print(perch_treatments)

perch_control_tel  <- perch_muddyfoot_tel[1:15]
perch_mix_tel      <- perch_muddyfoot_tel[16:30]
perch_control_akdes <- perch_muddyfoot_akdes[1:15]
perch_mix_akdes     <- perch_muddyfoot_akdes[16:30]

#> 3.2 Roach ####
roach_treatments <- sapply(roach_muddyfoot_tel, function(tel) unique(tel$treatment))
print(roach_treatments)

roach_control_tel  <- roach_muddyfoot_tel[1:13]
roach_mix_tel      <- roach_muddyfoot_tel[14:26]
roach_control_akdes <- roach_muddyfoot_akdes[1:13]
roach_mix_akdes     <- roach_muddyfoot_akdes[14:26]

#> 3.3 Pike ####
pike_control_tel  <- pike_muddyfoot_tel[1:3]
pike_mix_tel      <- pike_muddyfoot_tel[4:6]
pike_control_akdes <- pike_muddyfoot_akdes[1:3]
pike_mix_akdes     <- pike_muddyfoot_akdes[4:6]

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
                                    saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_perch/",
                                                              names(perch_control_tel)[i], "_habitat_rsf.rds"))
                                    rsf_model
                                  }
stopCluster(cl)

failed_perch_control <- check_rsf_results(rsf_perch_control_list, perch_control_tel, "Perch Control")
rsf_perch_control_list <- rsf_perch_control_list[!failed_perch_control]
names(rsf_perch_control_list) <- names(perch_control_tel)[!failed_perch_control]
saveRDS(rsf_perch_control_list, paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_list.rds"))

# Exposed
cl <- makeCluster(length(perch_mix_tel))
registerDoParallel(cl)
rsf_perch_mix_list <- foreach(i = seq_along(perch_mix_tel),
                              .packages = c("ctmm", "raster"),
                              .errorhandling = "pass") %dopar% {
                                rsf_model <- rsf.fit(perch_mix_tel[[i]],
                                                     perch_mix_akdes[[i]],
                                                     R = list(habitat1 = habitat_raster_final))
                                saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_perch/",
                                                          names(perch_mix_tel)[i], "_habitat_rsf.rds"))
                                rsf_model
                              }
stopCluster(cl)

failed_perch_mix <- check_rsf_results(rsf_perch_mix_list, perch_mix_tel, "Perch Exposed")
rsf_perch_mix_list <- rsf_perch_mix_list[!failed_perch_mix]
names(rsf_perch_mix_list) <- names(perch_mix_tel)[!failed_perch_mix]
saveRDS(rsf_perch_mix_list, paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_list.rds"))

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
                                    saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_roach/",
                                                              names(roach_control_tel)[i], "_habitat_rsf.rds"))
                                    rsf_model
                                  }
stopCluster(cl)

failed_roach_control <- check_rsf_results(rsf_roach_control_list, roach_control_tel, "Roach Control")
rsf_roach_control_list <- rsf_roach_control_list[!failed_roach_control]
names(rsf_roach_control_list) <- names(roach_control_tel)[!failed_roach_control]
saveRDS(rsf_roach_control_list, paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_list.rds"))

# Exposed
cl <- makeCluster(length(roach_mix_tel))
registerDoParallel(cl)
rsf_roach_mix_list <- foreach(i = seq_along(roach_mix_tel),
                              .packages = c("ctmm", "raster"),
                              .errorhandling = "pass") %dopar% {
                                rsf_model <- rsf.fit(roach_mix_tel[[i]],
                                                     roach_mix_akdes[[i]],
                                                     R = list(habitat1 = habitat_raster_final))
                                saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_roach/",
                                                          names(roach_mix_tel)[i], "_habitat_rsf.rds"))
                                rsf_model
                              }
stopCluster(cl)

failed_roach_mix <- check_rsf_results(rsf_roach_mix_list, roach_mix_tel, "Roach Exposed")
rsf_roach_mix_list <- rsf_roach_mix_list[!failed_roach_mix]
names(rsf_roach_mix_list) <- names(roach_mix_tel)[!failed_roach_mix]
saveRDS(rsf_roach_mix_list, paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_list.rds"))

#> 4.3 Pike ####
cat("\n=== Fitting Pike RSFs ===\n")

# All pike together (too few for treatment split)
cl <- makeCluster(length(pike_muddyfoot_tel))
registerDoParallel(cl)
rsf_pike_list <- foreach(i = seq_along(pike_muddyfoot_tel),
                         .packages = c("ctmm", "raster"),
                         .errorhandling = "pass") %dopar% {
                           rsf_model <- rsf.fit(pike_muddyfoot_tel[[i]],
                                                pike_muddyfoot_akdes[[i]],
                                                R = list(habitat1 = habitat_raster_final))
                           saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_pike/",
                                                     names(pike_muddyfoot_tel)[i], "_habitat_rsf.rds"))
                           rsf_model
                         }
stopCluster(cl)

failed_pike <- check_rsf_results(rsf_pike_list, pike_muddyfoot_tel, "Pike")
rsf_pike_list <- rsf_pike_list[!failed_pike]
names(rsf_pike_list) <- names(pike_muddyfoot_tel)[!failed_pike]
saveRDS(rsf_pike_list, paste0(rsf_path, "muddyfoot_pike/rsf_pike_list.rds"))

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

rsf_perch_control_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_list.rds"))
rsf_perch_mix_list     <- readRDS(paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_list.rds"))

rsf_perch_control_mean <- mean(rsf_perch_control_list)
rsf_perch_exposed_mean <- mean(rsf_perch_mix_list)

rsf_coef_control_perch <- as.data.frame(t(summary(rsf_perch_control_mean)$CI[1,]))
rsf_coef_exposed_perch <- as.data.frame(t(summary(rsf_perch_exposed_mean)$CI[1,]))
rsf_coef_control_perch$treatment <- "Control"
rsf_coef_exposed_perch$treatment <- "Exposed"

perch_habitat_rsf_coefs <- rbind(rsf_coef_control_perch, rsf_coef_exposed_perch)

cat("Perch Control CI:\n"); print(rsf_coef_control_perch)
cat("Perch Exposed CI:\n"); print(rsf_coef_exposed_perch)

perch_habitat_rsf_plot <- ggplot(perch_habitat_rsf_coefs, aes(x = treatment, y = est)) +
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

print(perch_habitat_rsf_plot)
ggsave(paste0(figure_path, "perch_habitats_rsf_muddyfoot.png"),
       perch_habitat_rsf_plot, width = 8, height = 8, units = 'cm', dpi = 300)

#> 5.2 Roach ####
cat("\n=== ROACH ANALYSIS ===\n")

rsf_roach_control_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_list.rds"))
rsf_roach_mix_list     <- readRDS(paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_list.rds"))

rsf_roach_control_mean <- mean(rsf_roach_control_list)
rsf_roach_exposed_mean <- mean(rsf_roach_mix_list)

rsf_coef_control_roach <- as.data.frame(t(summary(rsf_roach_control_mean)$CI[1,]))
rsf_coef_exposed_roach <- as.data.frame(t(summary(rsf_roach_exposed_mean)$CI[1,]))
rsf_coef_control_roach$treatment <- "Control"
rsf_coef_exposed_roach$treatment <- "Exposed"

roach_habitat_rsf_coefs <- rbind(rsf_coef_control_roach, rsf_coef_exposed_roach)

cat("Roach Control CI:\n"); print(rsf_coef_control_roach)
cat("Roach Exposed CI:\n"); print(rsf_coef_exposed_roach)

roach_habitat_rsf_plot <- ggplot(roach_habitat_rsf_coefs, aes(x = treatment, y = est)) +
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
ggsave(paste0(figure_path, "roach_habitats_rsf_muddyfoot.png"),
       roach_habitat_rsf_plot, width = 8, height = 8, units = 'cm', dpi = 300)

#> 5.3 Pike ####
cat("\n=== PIKE ANALYSIS ===\n")

rsf_pike_list <- readRDS(paste0(rsf_path, "muddyfoot_pike/rsf_pike_list.rds"))

# Split into control and exposed for plotting
rsf_pike_control_list <- rsf_pike_list[1:3]
rsf_pike_mix_list     <- rsf_pike_list[4:6]

rsf_pike_control_mean <- mean(rsf_pike_control_list)
rsf_pike_exposed_mean <- mean(rsf_pike_mix_list)

rsf_coef_control_pike <- as.data.frame(t(summary(rsf_pike_control_mean)$CI[1,]))
rsf_coef_exposed_pike <- as.data.frame(t(summary(rsf_pike_exposed_mean)$CI[1,]))
rsf_coef_control_pike$treatment <- "Control"
rsf_coef_exposed_pike$treatment <- "Exposed"

pike_habitat_rsf_coefs <- rbind(rsf_coef_control_pike, rsf_coef_exposed_pike)

cat("Pike Control CI:\n"); print(rsf_coef_control_pike)
cat("Pike Exposed CI:\n"); print(rsf_coef_exposed_pike)

pike_habitat_rsf_plot <- ggplot(pike_habitat_rsf_coefs, aes(x = treatment, y = est)) +
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

print(pike_habitat_rsf_plot)
ggsave(paste0(figure_path, "pike_habitats_rsf_muddyfoot.png"),
       pike_habitat_rsf_plot, width = 8, height = 8, units = 'cm', dpi = 300)

cat("\n=== ANALYSIS COMPLETE ===\n")
