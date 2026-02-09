# ---------------------------------------------------------#
# RESOURCE SELECTION FOR ARTIFICIAL HABITATS - MUDDYFOOT ###
# ---------------------------------------------------------#
# Date: 2026-02-06

### LIBRARIES ###
library(ctmm)
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(data.table)
library(parallel)    # Parallel computing support
library(foreach)     # Looping construct for parallel execution
library(doParallel)  # Parallel backend for foreach loops

### DIRECTORIES ###
polygon_path <- "./data/lake_params/polygons/"
ctmm_path <- "./data/ctmm_fits/"
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
telem_path <- "./data/telem_obj/muddyfoot/" 
rec_data_path <- "./data/lake_params/reciever_and_habitat_locations/"
enc_path <- "./data/encounters/muddyfoot/"
akde_path <- "./data/akdes/"
rsf_path <- "./data/rsfs/habitats/"
figure_path <- "./figures/muddyfoot/"

### LOAD DATA ###

# Receiver and habitat locations
mud_rec_locs_kml <- paste0(rec_data_path, "muddyfoot_rec_hab_locations.kml")
mud_rec_locs <- st_read(mud_rec_locs_kml)[1:5,]
mud_hab_locs <- st_read(mud_rec_locs_kml)[6:7,]

# Fish telemetry data
pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel_thinned_final.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel_thinned_final.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel_thinned_final.rds'))

# Fish ctmms 
pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_best_ctmm_model_fits.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_best_ctmm_model_fits.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_best_ctmm_model_fits.rds"))

# Fish akdes
pike_muddyfoot_akdes <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_muddyfoot_akdes.rds"))
perch_muddyfoot_akdes <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/perch_muddyfoot_akdes.rds"))
roach_muddyfoot_akdes <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/roach_muddyfoot_akdes.rds"))

# Load PKDEs
perch_control_PKDE <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/perch_control_PKDE.rds"))
perch_mix_PKDE <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/perch_mix_PKDE.rds"))
roach_control_PKDE <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/roach_control_PKDE.rds"))
roach_mix_PKDE <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/roach_mix_PKDE.rds"))
pike_control_PKDE <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_control_PKDE.rds"))
pike_mix_PKDE <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_mix_PKDE.rds"))
pike_total_PKDE <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/pike_total_PKDE.rds"))

# Load lake polygon
muddyfoot_polygon <- st_read(paste0(polygon_path, "muddyfoot_polygon.gpkg"))

#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#---------------------------------------------#
# 1. Plot population AKDEs with habitats ####
#---------------------------------------------#

# Function to generate the plot with habitats
generate_ud_plot_whabitats <- function(pkde_data, bbox, hab_locs, title = "") {
  # Extract the CDF (cumulative distribution function) - default
  ud_raster <- raster(pkde_data, DF = "CDF")  # or just raster(pkde_data)
  
  masked_ud_raster <- mask(ud_raster, bbox)
  
  # Convert to dataframe
  ud_df <- as.data.frame(masked_ud_raster, xy = TRUE, na.rm = TRUE)
  colnames(ud_df) <- c("x", "y", "value")
  
  # Invert the CDF so high values = high use
  # CDF: low values = core areas, high values = periphery
  # Inverted: high values = core areas, low values = periphery
  ud_df$value <- 1 - ud_df$value
  
  ggplot() +
    geom_sf(data = bbox, color = "black") +
    geom_tile(data = ud_df, aes(x = x, y = y, fill = value), alpha = 0.6) +
    geom_sf(data = hab_locs, color = "green", size = 3, fill = NA, shape = 3, stroke = 2) + 
    scale_fill_viridis_c(na.value = 'transparent', option = 'magma', direction = -1) +
    coord_sf() +
    theme_classic() +
    labs(fill = "Utilization Distribution", title = title, x = "", y = '') +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}


# Function to generate the plot without habitats
generate_ud_plot <- function(pkde_data, bbox, hab_locs, title = "") {
  ud_raster <- raster(pkde_data)
  masked_ud_raster <- mask(ud_raster, bbox)
  ud_df <- as.data.frame(as(masked_ud_raster, "SpatialPixelsDataFrame"))
  colnames(ud_df) <- c("value", "x", "y")
  
  ggplot() +
    geom_sf(data = bbox, color = "black") +
    geom_tile(data = ud_df, aes(x = x, y = y, fill = value), alpha = 0.6) +
    scale_fill_viridis_c(na.value = 'transparent', option = 'magma') +
    coord_sf() +
    theme_classic() +
    labs(fill = "Utilization Distribution", title = title, x = "", y = '') +
    theme(legend.position = "none",
          axis.text = element_text(size = 8, color = 'black'),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

# Transform bbox
muddyfoot_bbox_perch <- st_transform(muddyfoot_polygon, crs(raster(perch_control_PKDE)))
muddyfoot_bbox_roach <- st_transform(muddyfoot_polygon, crs(raster(roach_control_PKDE)))
muddyfoot_bbox_pike <- st_transform(muddyfoot_polygon, crs(raster(pike_total_PKDE)))

#> 1.1. Perch ####

## CONTROL ##
perch_control_plot <- generate_ud_plot(perch_control_PKDE, muddyfoot_bbox_perch, 
                                       mud_hab_locs, "Perch Control")

perch_control_plot_habitats <- generate_ud_plot_whabitats(perch_control_PKDE, 
                                                          muddyfoot_bbox_perch, 
                                                          mud_hab_locs, 
                                                          "Perch Control")

print(perch_control_plot_habitats)

ggsave(file = paste0(figure_path, "UD_plots/perch_control_UD_habitats_muddyfoot.png"), 
       plot = perch_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)

## EXPOSED ##
perch_exposed_plot <- generate_ud_plot(perch_mix_PKDE, muddyfoot_bbox_perch, 
                                       mud_hab_locs, "Perch Exposed")
perch_exposed_plot_habitats <- generate_ud_plot_whabitats(perch_mix_PKDE, 
                                                          muddyfoot_bbox_perch, 
                                                          mud_hab_locs, 
                                                          "Perch Exposed")

print(perch_exposed_plot_habitats)

ggsave(file = paste0(figure_path, "UD_plots/perch_exposed_UD_habitats_muddyfoot.png"), 
       plot = perch_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


#> 1.2. Roach ####
## CONTROL ##
roach_control_plot <- generate_ud_plot(roach_control_PKDE, muddyfoot_bbox_roach, 
                                       mud_hab_locs, "Roach Control")
roach_control_plot_habitats <- generate_ud_plot_whabitats(roach_control_PKDE, 
                                                          muddyfoot_bbox_roach, 
                                                          mud_hab_locs, 
                                                          "Roach Control")

print(roach_control_plot_habitats)

ggsave(file = paste0(figure_path, "UD_plots/roach_control_UD_habitats_muddyfoot.png"), 
       plot = roach_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)

## EXPOSED ##
roach_exposed_plot <- generate_ud_plot(roach_mix_PKDE, muddyfoot_bbox_roach, 
                                       mud_hab_locs, "Roach Exposed")
roach_exposed_plot_habitats <- generate_ud_plot_whabitats(roach_mix_PKDE, 
                                                          muddyfoot_bbox_roach, 
                                                          mud_hab_locs, 
                                                          "Roach Exposed")

print(roach_exposed_plot_habitats)

ggsave(file = paste0(figure_path, "UD_plots/roach_exposed_UD_habitats_muddyfoot.png"), 
       plot = roach_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


#> 1.3. Pike ####

## CONTROL ##
pike_control_plot <- generate_ud_plot(pike_control_PKDE, muddyfoot_bbox_pike, 
                                      mud_hab_locs, "Pike Control")
pike_control_plot_habitats <- generate_ud_plot_whabitats(pike_control_PKDE, 
                                                         muddyfoot_bbox_pike, 
                                                         mud_hab_locs, 
                                                         "Pike Control")

print(pike_control_plot_habitats)

ggsave(file = paste0(figure_path, "UD_plots/pike_control_UD_habitats_muddyfoot.png"), 
       plot = pike_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## EXPOSED ##
pike_exposed_plot <- generate_ud_plot(pike_mix_PKDE, muddyfoot_bbox_pike, 
                                      mud_hab_locs, "Pike Exposed")

pike_exposed_plot_habitats <- generate_ud_plot_whabitats(pike_mix_PKDE, 
                                                         muddyfoot_bbox_pike, 
                                                         mud_hab_locs, 
                                                         "Pike Exposed")

print(pike_exposed_plot_habitats)


ggsave(file = paste0(figure_path, "UD_plots/pike_exposed_UD_habitats_muddyfoot.png"), 
       plot = pike_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## TOTAL ##
pike_total_plot <- generate_ud_plot(pike_total_PKDE, muddyfoot_bbox_pike, 
                                    mud_hab_locs, "Pike Total")

pike_total_plot_habitats <- generate_ud_plot_whabitats(pike_total_PKDE, 
                                                       muddyfoot_bbox_pike, 
                                                       mud_hab_locs, 
                                                       "Pike Total")

print(pike_total_plot_habitats)


ggsave(file = paste0(figure_path, "UD_plots/pike_total_UD_habitats_muddyfoot.png"), 
       plot = pike_total_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)

#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#-----------------------------------------#
# 2. Create lake and habitat raster - IMPROVED ####
#-----------------------------------------#

# Transform the CRS to UTM Zone 33N (EPSG:32633)
polyProj <- st_transform(muddyfoot_polygon, crs = "EPSG:32633")

# Resolution in meters
resolution <- 0.25  # meters
lake_raster <- rast(ext(polyProj), 
                    res = resolution, 
                    crs = "EPSG:32633")

# Rasterize the lake polygon
lake_raster <- rasterize(polyProj, lake_raster, background = NA)

cat("Lake raster created.\n")

# Transform habitats to UTM
mud_hab_locs_utm <- st_transform(mud_hab_locs, crs = "EPSG:32633")

# Custom function to create perfect 2m × 2m squares
create_square_buffer <- function(points, size = 2) {
  squares <- lapply(1:nrow(points), function(i) {
    center <- st_coordinates(points[i,])
    bbox <- st_bbox(c(xmin = center[1] - size/2, 
                      xmax = center[1] + size/2,
                      ymin = center[2] - size/2, 
                      ymax = center[2] + size/2),
                    crs = st_crs(points))
    st_as_sfc(bbox)
  })
  st_sf(geometry = do.call(c, squares), crs = st_crs(points))
}

# Create exact 2m × 2m squares in UTM
habitat_patches_buffered <- create_square_buffer(mud_hab_locs_utm, size = 2.0)

# Convert to SpatVector
habitat_patches_vect <- vect(habitat_patches_buffered)

# Rasterize with buffered habitats IN UTM
habitat_raster <- rasterize(habitat_patches_vect, 
                            lake_raster, 
                            field = 1, 
                            background = 0)

# Mask to lake polygon
habitat_raster <- mask(habitat_raster, polyProj)

# Visualize in UTM (perfect squares!)
plot(habitat_raster, main = "Habitat Raster with 2m × 2m Square Structures (UTM)")
plot(st_geometry(polyProj), add = TRUE, border = "black", lwd = 2)
plot(st_geometry(mud_hab_locs_utm), add = TRUE, col = "red", pch = 4, cex = 2)
plot(st_geometry(habitat_patches_buffered), add = TRUE, border = "green", lwd = 2)

cat("Habitat raster created with 2m × 2m square structures.\n")

# Convert to raster package format for ctmm (KEEP IN UTM!)
habitat_raster_final <- raster::raster(habitat_raster)

plot(habitat_raster_final)


#IF we need transform the raster to WGS84
# Convert raster package object to terra
habitat_raster_terra <- rast(habitat_raster_final)

# Check current CRS
cat("Current CRS:\n")
print(crs(habitat_raster_terra))

# Transform to WGS84
habitat_raster_wgs84_terra <- project(habitat_raster_terra, 
                                      "EPSG:4326",  # WGS84
                                      method = "near")  # nearest neighbor

# Force binary values
habitat_raster_wgs84_terra[habitat_raster_wgs84_terra > 0.5] <- 1
habitat_raster_wgs84_terra[habitat_raster_wgs84_terra <= 0.5 & !is.na(habitat_raster_wgs84_terra)] <- 0

# Convert back to raster package format for ctmm
habitat_raster_wgs84 <- raster::raster(habitat_raster_wgs84_terra)

# Verify the transformation
cat("\nOriginal raster (UTM) extent:\n")
print(raster::extent(habitat_raster_final))

cat("\nTransformed raster (WGS84) extent:\n")
print(raster::extent(habitat_raster_wgs84))

cat("\nTelemetry extent:\n")
cat("  lon:", range(perch_control_tel[[1]]$longitude), "\n")
cat("  lat:", range(perch_control_tel[[1]]$latitude), "\n")

# Test extraction with WGS84 raster
coords_test <- data.frame(x = perch_control_tel[[1]]$longitude[1:100],
                          y = perch_control_tel[[1]]$latitude[1:100])
hab_values_test <- raster::extract(habitat_raster_wgs84, coords_test)

cat("\nTest extraction with WGS84 raster:\n")
cat("  Total points:", length(hab_values_test), "\n")
cat("  NA values:", sum(is.na(hab_values_test)), "\n")
cat("  Unique values:", unique(hab_values_test), "\n")


plot(habitat_raster_wgs84)

# If this works, save it
raster::writeRaster(habitat_raster_wgs84, 
                    "./data/lake_coords/muddyfoot_habitat_raster_2m_squares_WGS84.grd", 
                    overwrite = TRUE)

cat("\n✓ Habitat raster transformed to WGS84 and saved\n")


#---------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------#
# 3. Split telemetry objects by treatment ####
#---------------------------------------------#

#> 3.1 Perch ########################

# Split telemetry objects based on verified treatment assignment
perch_control_tel <- perch_muddyfoot_tel[1:15]
perch_mix_tel <- perch_muddyfoot_tel[16:30]

# Split AKDEs
perch_control_akdes <- perch_muddyfoot_akdes[1:15]
perch_mix_akdes <- perch_muddyfoot_akdes[16:30]

#> 3.2 Roach ########################

# Split telemetry objects
roach_control_tel <- roach_muddyfoot_tel[1:13]
roach_mix_tel <- roach_muddyfoot_tel[14:26]

# Split AKDEs
roach_control_akdes <- roach_muddyfoot_akdes[1:13]
roach_mix_akdes <- roach_muddyfoot_akdes[14:26]

#> 3.3 Pike ########################

# Split telemetry objects
pike_control_tel <- pike_muddyfoot_tel[1:3]
pike_mix_tel <- pike_muddyfoot_tel[4:6]

# Split AKDEs
pike_control_akdes <- pike_muddyfoot_akdes[1:3]
pike_mix_akdes <- pike_muddyfoot_akdes[4:6]

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------

#-----------------------------------------#
# 3. TEST RSF on single individual ####
#-----------------------------------------#

cat("\n=== TESTING RSF ===\n")

cat("Testing on first perch control individual...\n")

# Test 1: Check CRS compatibility
cat("\nTest 1: Checking CRS...\n")
cat("  Telemetry CRS:", projection(perch_control_tel[[1]]), "\n")
cat("  Raster CRS:", crs(habitat_raster_final), "\n")

# Test 2: Check raster values
cat("\nTest 2: Checking raster values...\n")
raster_values <- unique(values(habitat_raster_final))
cat("  Unique raster values:", raster_values[!is.na(raster_values)], "\n")
cat("  Number of habitat cells:", sum(values(habitat_raster) == 1, na.rm = TRUE), "\n")
cat("  Number of non-habitat cells:", sum(values(habitat_raster) == 0, na.rm = TRUE), "\n")

# Test 3: Check if fish locations overlap with raster
cat("\nTest 3: Checking spatial overlap...\n")
tel_extent <- c(range(perch_control_tel[[1]]$longitude), 
                range(perch_control_tel[[1]]$latitude))
raster_extent <- as.vector(extent(habitat_raster_final))
cat("  Telemetry extent (lon):", tel_extent[1:2], "\n")
cat("  Telemetry extent (lat):", tel_extent[3:4], "\n")
cat("  Raster extent:", raster_extent, "\n")

# Test 4: Try fitting RSF on single individual
cat("\nTest 4: Fitting test RSF...\n")
test_rsf <- tryCatch({
  rsf.fit(perch_control_tel[[1]], 
          perch_control_akdes[[1]], 
          R = list(habitat1 = habitat_raster_final))
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  return(NULL)
})

if(!is.null(test_rsf)) {
  cat("  ✓ Test RSF fitted successfully!\n")
  cat("\nTest RSF Summary:\n")
  print(summary(test_rsf))
  
  # Save test result
  saveRDS(test_rsf, paste0(rsf_path, "test_rsf_result.rds"))
  
  cat("\n=== Tests passed! Proceeding with full analysis ===\n")
  
} else {
  cat("\n  ✗ Test RSF FAILED!\n")
  cat("\nPossible issues:\n")
  cat("  1. CRS mismatch between telemetry and raster\n")
  cat("  2. Raster values not binary (0/1)\n")
  cat("  3. No spatial overlap between data and raster\n")
  cat("  4. AKDE object incompatible with telemetry\n")
  cat("\nStopping analysis. Fix issues above before proceeding.\n")
  stop("RSF test failed. See error messages above.")
}

summary(test_rsf)


#-----------------------#
# 4. Fit RSF models ####
#-----------------------#

#> 4.1 Perch ####

cat("\n=== Fitting Perch RSFs ===\n")

# Control
cat("Perch Control...\n")
n_workers <- length(perch_control_tel)
cl <- makeCluster(n_workers)
registerDoParallel(cl)

cat("Running RSF on", length(perch_control_tel), "fish with", n_workers, "workers\n")

# Run RSF
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

# Check results
cat("\nTotal models returned:", length(rsf_perch_control_list), "\n")

# Check for failures
failed <- sapply(rsf_perch_control_list, function(x) {
  !inherits(x, "ctmm") || is.null(x) || inherits(x, "error")
})

cat("Failed:", sum(failed), "\n")
cat("Succeeded:", sum(!failed), "\n")

if(any(failed)) {
  cat("\nFailed IDs:", names(perch_control_tel)[failed], "\n")
  
  # Show error messages for failed models
  for(i in which(failed)) {
    cat("\nModel", i, "-", names(perch_control_tel)[i], ":\n")
    if(inherits(rsf_perch_control_list[[i]], "error")) {
      cat("  Error:", rsf_perch_control_list[[i]]$message, "\n")
    }
  }
}

# Remove failed models and fix names
rsf_perch_control_list <- rsf_perch_control_list[!failed]
names(rsf_perch_control_list) <- names(perch_control_tel)[!failed]

# Save successful models
saveRDS(rsf_perch_control_list, paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_list.rds"))

cat("\nFinal: Saved", length(rsf_perch_control_list), "successful models\n")



# Exposed
cat("Perch Exposed...\n")
cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_perch_mix_list <- foreach(i = seq_along(perch_mix_tel), 
                              .packages = c("ctmm", "raster"),
                              .errorhandling = "pass") %dopar% {
                                rsf_model <- rsf.fit(perch_mix_tel[[i]], 
                                                     perch_mix_akdes[[i]], 
                                                     R = list(habitat1 = habitat_raster))
                                saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_perch/", 
                                                          names(perch_mix_tel)[i], "_habitat_rsf.rds"))
                                rsf_model
                              }
stopCluster(cl)

failed_mix <- sapply(rsf_perch_mix_list, function(x) inherits(x, "error"))
if(any(failed_mix)) {
  cat("  Warning:", sum(failed_mix), "models failed\n")
  rsf_perch_mix_list <- rsf_perch_mix_list[!failed_mix]
}

names(rsf_perch_mix_list) <- names(perch_mix_tel)[1:length(rsf_perch_mix_list)]
saveRDS(rsf_perch_mix_list, paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_list.rds"))
cat("  Completed:", length(rsf_perch_mix_list), "models\n")

#> 4.2 Roach ####

cat("\n=== Fitting Roach RSFs ===\n")

# Control
cat("Roach Control...\n")
cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_roach_control_list <- foreach(i = seq_along(roach_control_tel), 
                                  .packages = c("ctmm", "raster"),
                                  .errorhandling = "pass") %dopar% {
                                    rsf_model <- rsf.fit(roach_control_tel[[i]], 
                                                         roach_control_akdes[[i]], 
                                                         R = list(habitat1 = habitat_raster))
                                    saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_roach/", 
                                                              names(roach_control_tel)[i], "_habitat_rsf.rds"))
                                    rsf_model
                                  }
stopCluster(cl)

failed_control <- sapply(rsf_roach_control_list, function(x) inherits(x, "error"))
if(any(failed_control)) {
  cat("  Warning:", sum(failed_control), "models failed\n")
  rsf_roach_control_list <- rsf_roach_control_list[!failed_control]
}

names(rsf_roach_control_list) <- names(roach_control_tel)[1:length(rsf_roach_control_list)]
saveRDS(rsf_roach_control_list, paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_list.rds"))
cat("  Completed:", length(rsf_roach_control_list), "models\n")

# Exposed
cat("Roach Exposed...\n")
cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_roach_mix_list <- foreach(i = seq_along(roach_mix_tel), 
                              .packages = c("ctmm", "raster"),
                              .errorhandling = "pass") %dopar% {
                                rsf_model <- rsf.fit(roach_mix_tel[[i]], 
                                                     roach_mix_akdes[[i]], 
                                                     R = list(habitat1 = habitat_raster))
                                saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_roach/", 
                                                          names(roach_mix_tel)[i], "_habitat_rsf.rds"))
                                rsf_model
                              }
stopCluster(cl)

failed_mix <- sapply(rsf_roach_mix_list, function(x) inherits(x, "error"))
if(any(failed_mix)) {
  cat("  Warning:", sum(failed_mix), "models failed\n")
  rsf_roach_mix_list <- rsf_roach_mix_list[!failed_mix]
}

names(rsf_roach_mix_list) <- names(roach_mix_tel)[1:length(rsf_roach_mix_list)]
saveRDS(rsf_roach_mix_list, paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_list.rds"))
cat("  Completed:", length(rsf_roach_mix_list), "models\n")

#> 4.3 Pike ####

cat("\n=== Fitting Pike RSFs ===\n")

# Control
cat("Pike Control...\n")
cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_pike_control_list <- foreach(i = seq_along(pike_control_tel), 
                                 .packages = c("ctmm", "raster"),
                                 .errorhandling = "pass") %dopar% {
                                   rsf_model <- rsf.fit(pike_control_tel[[i]], 
                                                        pike_control_akdes[[i]], 
                                                        R = list(habitat1 = habitat_raster))
                                   saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_pike/", 
                                                             names(pike_control_tel)[i], "_habitat_rsf.rds"))
                                   rsf_model
                                 }
stopCluster(cl)

failed_control <- sapply(rsf_pike_control_list, function(x) inherits(x, "error"))
if(any(failed_control)) {
  cat("  Warning:", sum(failed_control), "models failed\n")
  rsf_pike_control_list <- rsf_pike_control_list[!failed_control]
}

names(rsf_pike_control_list) <- names(pike_control_tel)[1:length(rsf_pike_control_list)]
saveRDS(rsf_pike_control_list, paste0(rsf_path, "muddyfoot_pike/rsf_pike_control_list.rds"))
cat("  Completed:", length(rsf_pike_control_list), "models\n")

# Exposed
cat("Pike Exposed...\n")
cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_pike_mix_list <- foreach(i = seq_along(pike_mix_tel), 
                             .packages = c("ctmm", "raster"),
                             .errorhandling = "pass") %dopar% {
                               rsf_model <- rsf.fit(pike_mix_tel[[i]], 
                                                    pike_mix_akdes[[i]], 
                                                    R = list(habitat1 = habitat_raster))
                               saveRDS(rsf_model, paste0(rsf_path, "muddyfoot_pike/", 
                                                         names(pike_mix_tel)[i], "_habitat_rsf.rds"))
                               rsf_model
                             }
stopCluster(cl)

failed_mix <- sapply(rsf_pike_mix_list, function(x) inherits(x, "error"))
if(any(failed_mix)) {
  cat("  Warning:", sum(failed_mix), "models failed\n")
  rsf_pike_mix_list <- rsf_pike_mix_list[!failed_mix]
}

names(rsf_pike_mix_list) <- names(pike_mix_tel)[1:length(rsf_pike_mix_list)]
saveRDS(rsf_pike_mix_list, paste0(rsf_path, "muddyfoot_pike/rsf_pike_mix_list.rds"))
cat("  Completed:", length(rsf_pike_mix_list), "models\n")

cat("\n=== RSF Analysis Complete ===\n")


















#-----------------------#
# 4. Fit RSF models - IMPROVED ####
#-----------------------#

# Enhanced RSF fitting function with diagnostics
fit_rsf_with_diagnostics <- function(tel_obj, akde_obj, raster_obj, id_name, save_path) {
  result <- list(
    id = id_name,
    success = FALSE,
    model = NULL,
    diagnostics = NULL,
    error = NULL
  )
  
  tryCatch({
    # Fit RSF
    rsf_model <- rsf.fit(
      tel_obj, 
      akde_obj, 
      R = list(habitat1 = raster_obj)
    )
    
    # Extract diagnostics
    summ <- summary(rsf_model)
    
    diagnostics <- list(
      coefficient = summ$CI["habitat1", "est"],
      se = (summ$CI["habitat1", "high"] - summ$CI["habitat1", "low"]) / (2 * 1.96),
      ci_low = summ$CI["habitat1", "low"],
      ci_high = summ$CI["habitat1", "high"],
      significant = !(summ$CI["habitat1", "low"] <= 0 & summ$CI["habitat1", "high"] >= 0),
      effective_n = summ$DOF["area"] / 2,
      ci_width = summ$CI["habitat1", "high"] - summ$CI["habitat1", "low"]
    )
    
    # Save individual model
    saveRDS(rsf_model, file = paste0(save_path, id_name, "_habitat_rsf.rds"))
    
    result$success <- TRUE
    result$model <- rsf_model
    result$diagnostics <- diagnostics
    
    cat("  ✓", id_name, "completed successfully\n")
    
  }, error = function(e) {
    warning(paste("RSF fitting failed for", id_name, ":", e$message))
    result$error <- e$message
    cat("  ✗", id_name, "FAILED:", e$message, "\n")
  })
  
  return(result)
}

#> 4.1 Perch - IMPROVED ####

cat("\n=== Fitting Perch RSF Models ===\n")

### CONTROL ###
cat("\nFitting Perch Control RSFs...\n")

cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_perch_control_results <- foreach(
  i = seq_along(perch_control_tel), 
  .packages = c("ctmm", "raster"),
  .errorhandling = "pass"
) %dopar% {
  fit_rsf_with_diagnostics(
    perch_control_tel[[i]], 
    perch_control_akdes[[i]], 
    habitat_raster,
    names(perch_control_tel)[i],
    paste0(rsf_path, "muddyfoot_perch/")
  )
}

stopCluster(cl)

# Extract successful models
rsf_perch_control_list <- lapply(rsf_perch_control_results, function(x) x$model)
rsf_perch_control_list <- rsf_perch_control_list[!sapply(rsf_perch_control_list, is.null)]
names(rsf_perch_control_list) <- names(perch_control_tel)[1:length(rsf_perch_control_list)]

# Compile diagnostics
perch_control_diag <- do.call(rbind, lapply(rsf_perch_control_results, function(x) {
  if(x$success) data.frame(x$diagnostics, stringsAsFactors = FALSE) else NULL
}))
perch_control_diag$id <- names(perch_control_tel)[1:nrow(perch_control_diag)]
perch_control_diag$treatment <- "Control"

cat("\nPerch Control RSF Summary:\n")
cat("  Successfully fit:", length(rsf_perch_control_list), "/", length(perch_control_tel), "\n")
cat("  Significant selection:", sum(perch_control_diag$significant), "individuals\n")
cat("  Mean coefficient:", round(mean(perch_control_diag$coefficient), 3), "\n")

# Save
saveRDS(rsf_perch_control_list, 
        paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_list.rds"))
saveRDS(perch_control_diag,
        paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_diagnostics.rds"))

### EXPOSED ###
cat("\nFitting Perch Exposed RSFs...\n")

cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_perch_mix_results <- foreach(
  i = seq_along(perch_mix_tel), 
  .packages = c("ctmm", "raster"),
  .errorhandling = "pass"
) %dopar% {
  fit_rsf_with_diagnostics(
    perch_mix_tel[[i]], 
    perch_mix_akdes[[i]], 
    habitat_raster,
    names(perch_mix_tel)[i],
    paste0(rsf_path, "muddyfoot_perch/")
  )
}

stopCluster(cl)

# Extract and compile
rsf_perch_mix_list <- lapply(rsf_perch_mix_results, function(x) x$model)
rsf_perch_mix_list <- rsf_perch_mix_list[!sapply(rsf_perch_mix_list, is.null)]
names(rsf_perch_mix_list) <- names(perch_mix_tel)[1:length(rsf_perch_mix_list)]

perch_mix_diag <- do.call(rbind, lapply(rsf_perch_mix_results, function(x) {
  if(x$success) data.frame(x$diagnostics, stringsAsFactors = FALSE) else NULL
}))
perch_mix_diag$id <- names(perch_mix_tel)[1:nrow(perch_mix_diag)]
perch_mix_diag$treatment <- "Exposed"

cat("\nPerch Exposed RSF Summary:\n")
cat("  Successfully fit:", length(rsf_perch_mix_list), "/", length(perch_mix_tel), "\n")
cat("  Significant selection:", sum(perch_mix_diag$significant), "individuals\n")
cat("  Mean coefficient:", round(mean(perch_mix_diag$coefficient), 3), "\n")

saveRDS(rsf_perch_mix_list, 
        paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_list.rds"))
saveRDS(perch_mix_diag,
        paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_diagnostics.rds"))

#-----------------------------------------------------------------------------------------#

#> 4.2 Roach - IMPROVED ####

cat("\n=== Fitting Roach RSF Models ===\n")

### CONTROL ###
cat("\nFitting Roach Control RSFs...\n")

cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_roach_control_results <- foreach(
  i = seq_along(roach_control_tel), 
  .packages = c("ctmm", "raster"),
  .errorhandling = "pass"
) %dopar% {
  fit_rsf_with_diagnostics(
    roach_control_tel[[i]], 
    roach_control_akdes[[i]], 
    habitat_raster,
    names(roach_control_tel)[i],
    paste0(rsf_path, "muddyfoot_roach/")
  )
}

stopCluster(cl)

rsf_roach_control_list <- lapply(rsf_roach_control_results, function(x) x$model)
rsf_roach_control_list <- rsf_roach_control_list[!sapply(rsf_roach_control_list, is.null)]
names(rsf_roach_control_list) <- names(roach_control_tel)[1:length(rsf_roach_control_list)]

roach_control_diag <- do.call(rbind, lapply(rsf_roach_control_results, function(x) {
  if(x$success) data.frame(x$diagnostics, stringsAsFactors = FALSE) else NULL
}))
roach_control_diag$id <- names(roach_control_tel)[1:nrow(roach_control_diag)]
roach_control_diag$treatment <- "Control"

cat("\nRoach Control RSF Summary:\n")
cat("  Successfully fit:", length(rsf_roach_control_list), "/", length(roach_control_tel), "\n")
cat("  Significant selection:", sum(roach_control_diag$significant), "individuals\n")
cat("  Mean coefficient:", round(mean(roach_control_diag$coefficient), 3), "\n")

saveRDS(rsf_roach_control_list, 
        paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_list.rds"))
saveRDS(roach_control_diag,
        paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_diagnostics.rds"))

### EXPOSED ###
cat("\nFitting Roach Exposed RSFs...\n")

cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_roach_mix_results <- foreach(
  i = seq_along(roach_mix_tel), 
  .packages = c("ctmm", "raster"),
  .errorhandling = "pass"
) %dopar% {
  fit_rsf_with_diagnostics(
    roach_mix_tel[[i]], 
    roach_mix_akdes[[i]], 
    habitat_raster,
    names(roach_mix_tel)[i],
    paste0(rsf_path, "muddyfoot_roach/")
  )
}

stopCluster(cl)

rsf_roach_mix_list <- lapply(rsf_roach_mix_results, function(x) x$model)
rsf_roach_mix_list <- rsf_roach_mix_list[!sapply(rsf_roach_mix_list, is.null)]
names(rsf_roach_mix_list) <- names(roach_mix_tel)[1:length(rsf_roach_mix_list)]

roach_mix_diag <- do.call(rbind, lapply(rsf_roach_mix_results, function(x) {
  if(x$success) data.frame(x$diagnostics, stringsAsFactors = FALSE) else NULL
}))
roach_mix_diag$id <- names(roach_mix_tel)[1:nrow(roach_mix_diag)]
roach_mix_diag$treatment <- "Exposed"

cat("\nRoach Exposed RSF Summary:\n")
cat("  Successfully fit:", length(rsf_roach_mix_list), "/", length(roach_mix_tel), "\n")
cat("  Significant selection:", sum(roach_mix_diag$significant), "individuals\n")
cat("  Mean coefficient:", round(mean(roach_mix_diag$coefficient), 3), "\n")

saveRDS(rsf_roach_mix_list, 
        paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_list.rds"))
saveRDS(roach_mix_diag,
        paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_diagnostics.rds"))

#-----------------------------------------------------------------------------------------------------------#

#> 4.3 Pike - IMPROVED ####

cat("\n=== Fitting Pike RSF Models ===\n")

### CONTROL ###
cat("\nFitting Pike Control RSFs...\n")

cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_pike_control_results <- foreach(
  i = seq_along(pike_control_tel), 
  .packages = c("ctmm", "raster"),
  .errorhandling = "pass"
) %dopar% {
  fit_rsf_with_diagnostics(
    pike_control_tel[[i]], 
    pike_control_akdes[[i]], 
    habitat_raster,
    names(pike_control_tel)[i],
    paste0(rsf_path, "muddyfoot_pike/")
  )
}

stopCluster(cl)

rsf_pike_control_list <- lapply(rsf_pike_control_results, function(x) x$model)
rsf_pike_control_list <- rsf_pike_control_list[!sapply(rsf_pike_control_list, is.null)]
names(rsf_pike_control_list) <- names(pike_control_tel)[1:length(rsf_pike_control_list)]

pike_control_diag <- do.call(rbind, lapply(rsf_pike_control_results, function(x) {
  if(x$success) data.frame(x$diagnostics, stringsAsFactors = FALSE) else NULL
}))
pike_control_diag$id <- names(pike_control_tel)[1:nrow(pike_control_diag)]
pike_control_diag$treatment <- "Control"

cat("\nPike Control RSF Summary:\n")
cat("  Successfully fit:", length(rsf_pike_control_list), "/", length(pike_control_tel), "\n")
cat("  Significant selection:", sum(pike_control_diag$significant), "individuals\n")
cat("  Mean coefficient:", round(mean(pike_control_diag$coefficient), 3), "\n")

saveRDS(rsf_pike_control_list, 
        paste0(rsf_path, "muddyfoot_pike/rsf_pike_control_list.rds"))
saveRDS(pike_control_diag,
        paste0(rsf_path, "muddyfoot_pike/rsf_pike_control_diagnostics.rds"))

### EXPOSED ###
cat("\nFitting Pike Exposed RSFs...\n")

cl <- makeCluster(min(3, detectCores() - 1))
registerDoParallel(cl)

rsf_pike_mix_results <- foreach(
  i = seq_along(pike_mix_tel), 
  .packages = c("ctmm", "raster"),
  .errorhandling = "pass"
) %dopar% {
  fit_rsf_with_diagnostics(
    pike_mix_tel[[i]], 
    pike_mix_akdes[[i]], 
    habitat_raster,
    names(pike_mix_tel)[i],
    paste0(rsf_path, "muddyfoot_pike/")
  )
}

stopCluster(cl)

rsf_pike_mix_list <- lapply(rsf_pike_mix_results, function(x) x$model)
rsf_pike_mix_list <- rsf_pike_mix_list[!sapply(rsf_pike_mix_list, is.null)]
names(rsf_pike_mix_list) <- names(pike_mix_tel)[1:length(rsf_pike_mix_list)]

pike_mix_diag <- do.call(rbind, lapply(rsf_pike_mix_results, function(x) {
  if(x$success) data.frame(x$diagnostics, stringsAsFactors = FALSE) else NULL
}))
pike_mix_diag$id <- names(pike_mix_tel)[1:nrow(pike_mix_diag)]
pike_mix_diag$treatment <- "Exposed"

cat("\nPike Exposed RSF Summary:\n")
cat("  Successfully fit:", length(rsf_pike_mix_list), "/", length(pike_mix_tel), "\n")
cat("  Significant selection:", sum(pike_mix_diag$significant), "individuals\n")
cat("  Mean coefficient:", round(mean(pike_mix_diag$coefficient), 3), "\n")

saveRDS(rsf_pike_mix_list, 
        paste0(rsf_path, "muddyfoot_pike/rsf_pike_mix_list.rds"))
saveRDS(pike_mix_diag,
        paste0(rsf_path, "muddyfoot_pike/rsf_pike_mix_diagnostics.rds"))

cat("\n=== All RSF models fitted successfully ===\n")

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------#
# 5. Explore RSF model results - IMPROVED ####
#-----------------------------------------#

# Function to compare treatments statistically
compare_treatments <- function(control_est, control_se, exposed_est, exposed_se) {
  diff <- exposed_est - control_est
  se_diff <- sqrt(control_se^2 + exposed_se^2)
  z_stat <- diff / se_diff
  p_value <- 2 * pnorm(-abs(z_stat))
  
  list(
    difference = diff,
    se_difference = se_diff,
    z_statistic = z_stat,
    p_value = p_value,
    significant = p_value < 0.05
  )
}

# Function to create comprehensive summary tables
create_rsf_summary_table <- function(species_name, control_diag, exposed_diag, 
                                     control_mean, exposed_mean, comparison) {
  
  # Individual-level summary
  ind_summary <- data.frame(
    Treatment = c("Control", "Exposed"),
    N = c(nrow(control_diag), nrow(exposed_diag)),
    N_significant = c(sum(control_diag$significant), sum(exposed_diag$significant)),
    Pct_significant = c(
      round(100 * sum(control_diag$significant) / nrow(control_diag), 1),
      round(100 * sum(exposed_diag$significant) / nrow(exposed_diag), 1)
    ),
    Mean_coef = c(mean(control_diag$coefficient), mean(exposed_diag$coefficient)),
    SD_coef = c(sd(control_diag$coefficient), sd(exposed_diag$coefficient)),
    Mean_eff_N = c(mean(control_diag$effective_n), mean(exposed_diag$effective_n)),
    Median_CI_width = c(median(control_diag$ci_width), median(exposed_diag$ci_width))
  )
  
  # Population-level summary
  control_summ <- summary(control_mean)$CI[1,]
  exposed_summ <- summary(exposed_mean)$CI[1,]
  
  pop_summary <- data.frame(
    Treatment = c("Control", "Exposed"),
    Population_est = c(control_summ["est"], exposed_summ["est"]),
    CI_low = c(control_summ["low"], exposed_summ["low"]),
    CI_high = c(control_summ["high"], exposed_summ["high"]),
    Significant = c(
      !(control_summ["low"] <= 0 & control_summ["high"] >= 0),
      !(exposed_summ["low"] <= 0 & exposed_summ["high"] >= 0)
    )
  )
  
  cat("\n========================================\n")
  cat(species_name, "RSF Summary\n")
  cat("========================================\n\n")
  cat("Individual-level summary:\n")
  print(ind_summary, row.names = FALSE)
  cat("\nPopulation-level summary:\n")
  print(pop_summary, row.names = FALSE)
  cat("\nTreatment comparison:\n")
  cat("  Difference (Exposed - Control):", round(comparison$difference, 3), "\n")
  cat("  Standard Error:", round(comparison$se_difference, 3), "\n")
  cat("  Z-statistic:", round(comparison$z_statistic, 3), "\n")
  cat("  P-value:", format.pval(comparison$p_value, digits = 3), "\n")
  cat("  Interpretation:", 
      ifelse(comparison$significant, 
             "*** Treatments differ significantly ***", 
             "No significant difference between treatments"), "\n")
  
  return(list(
    individual = ind_summary, 
    population = pop_summary, 
    comparison = comparison
  ))
}

#> 5.1. Perch - IMPROVED ####

cat("\n\n=== PERCH ANALYSIS ===\n")

# Load RSF results
rsf_perch_control_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_list.rds"))
rsf_perch_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_list.rds"))
perch_control_diag <- readRDS(paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_diagnostics.rds"))
perch_mix_diag <- readRDS(paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_diagnostics.rds"))

# Calculate population-level means
rsf_perch_control_mean <- mean(rsf_perch_control_list)
rsf_perch_exposed_mean <- mean(rsf_perch_exposed_list)

# Extract coefficients
rsf_coef_control_perch <- as.data.frame(t(summary(rsf_perch_control_mean)$CI[1,]))
rsf_coef_exposed_perch <- as.data.frame(t(summary(rsf_perch_exposed_mean)$CI[1,]))

rsf_coef_control_perch$treatment <- "Control"
rsf_coef_exposed_perch$treatment <- "Exposed"

perch_habitat_rsf_coefs <- rbind(rsf_coef_control_perch, rsf_coef_exposed_perch)

# Calculate SEs from CIs
control_se <- (rsf_coef_control_perch$high - rsf_coef_control_perch$low) / (2 * 1.96)
exposed_se <- (rsf_coef_exposed_perch$high - rsf_coef_exposed_perch$low) / (2 * 1.96)

# Statistical comparison
perch_comparison <- compare_treatments(
  rsf_coef_control_perch$est, control_se,
  rsf_coef_exposed_perch$est, exposed_se
)

# Create comprehensive summary
perch_summary <- create_rsf_summary_table(
  "Perch",
  perch_control_diag,
  perch_mix_diag,
  rsf_perch_control_mean,
  rsf_perch_exposed_mean,
  perch_comparison
)

# Save summary
saveRDS(perch_summary, paste0(rsf_path, "muddyfoot_perch/perch_rsf_summary.rds"))

# Enhanced population-level plot
perch_habitat_rsf_plot <- ggplot(perch_habitat_rsf_coefs, 
                                 aes(x = treatment, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = low, ymax = high), 
                width = 0.1, size = 1, color = "black") +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, color = "black", stroke = 1) +
  scale_shape_manual(values = c(21, 21)) +
  scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +
  coord_cartesian(ylim = c(min(perch_habitat_rsf_coefs$low) - 0.5, 
                           max(perch_habitat_rsf_coefs$high) + 0.5)) +
  labs(y = "Selection coefficient",
       title = "Perch Habitat Selection",
       subtitle = paste0("Δ = ", round(perch_comparison$difference, 2),
                         ", p = ", format.pval(perch_comparison$p_value, digits = 2))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

print(perch_habitat_rsf_plot)

ggsave(file = paste0(save_plots, "perch_habitats_rsf_muddyfoot.png"), 
       plot = perch_habitat_rsf_plot, 
       device = 'png',
       width = 10, 
       height = 10,
       units = 'cm',
       dpi = 300)

# Individual-level variation plot
perch_all_diag <- rbind(perch_control_diag, perch_mix_diag)

perch_individual_plot <- ggplot(perch_all_diag, 
                                aes(x = treatment, y = coefficient, color = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, fill = NA, linewidth = 0.8) +
  scale_color_manual(values = c("Control" = "#0072B2", "Exposed" = "#D55E00")) +
  labs(y = "Selection coefficient", 
       title = "Individual variation in perch habitat selection",
       subtitle = paste0("Control: n=", nrow(perch_control_diag), 
                         ", Exposed: n=", nrow(perch_mix_diag))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 'bold', size = 14),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(perch_individual_plot)

ggsave(file = paste0(save_plots, "perch_individual_variation_rsf.png"), 
       plot = perch_individual_plot, 
       device = 'png',
       width = 12, 
       height = 10,
       units = 'cm',
       dpi = 300)

#-------------------------------------------------------------------------------#

#> 5.2 Roach - IMPROVED ####

cat("\n\n=== ROACH ANALYSIS ===\n")

# Load RSF results
rsf_roach_control_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_list.rds"))
rsf_roach_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_list.rds"))
roach_control_diag <- readRDS(paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_diagnostics.rds"))
roach_mix_diag <- readRDS(paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_diagnostics.rds"))

# Calculate means
rsf_roach_control_mean <- mean(rsf_roach_control_list)
rsf_roach_exposed_mean <- mean(rsf_roach_exposed_list)

# Extract coefficients
rsf_coef_control_roach <- as.data.frame(t(summary(rsf_roach_control_mean)$CI[1,]))
rsf_coef_exposed_roach <- as.data.frame(t(summary(rsf_roach_exposed_mean)$CI[1,]))

rsf_coef_control_roach$treatment <- "Control"
rsf_coef_exposed_roach$treatment <- "Exposed"

roach_habitat_rsf_coefs <- rbind(rsf_coef_control_roach, rsf_coef_exposed_roach)

# Calculate SEs
control_se <- (rsf_coef_control_roach$high - rsf_coef_control_roach$low) / (2 * 1.96)
exposed_se <- (rsf_coef_exposed_roach$high - rsf_coef_exposed_roach$low) / (2 * 1.96)

# Statistical comparison
roach_comparison <- compare_treatments(
  rsf_coef_control_roach$est, control_se,
  rsf_coef_exposed_roach$est, exposed_se
)

# Create summary
roach_summary <- create_rsf_summary_table(
  "Roach",
  roach_control_diag,
  roach_mix_diag,
  rsf_roach_control_mean,
  rsf_roach_exposed_mean,
  roach_comparison
)

saveRDS(roach_summary, paste0(rsf_path, "muddyfoot_roach/roach_rsf_summary.rds"))

# Population-level plot
roach_habitat_rsf_plot <- ggplot(roach_habitat_rsf_coefs, 
                                 aes(x = treatment, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = low, ymax = high), 
                width = 0.1, size = 1, color = "black") +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, color = "black", stroke = 1) +
  scale_shape_manual(values = c(21, 21)) +
  scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +
  coord_cartesian(ylim = c(min(roach_habitat_rsf_coefs$low) - 0.5, 
                           max(roach_habitat_rsf_coefs$high) + 0.5)) +
  labs(y = "Selection coefficient",
       title = "Roach Habitat Selection",
       subtitle = paste0("Δ = ", round(roach_comparison$difference, 2),
                         ", p = ", format.pval(roach_comparison$p_value, digits = 2))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

print(roach_habitat_rsf_plot)

ggsave(file = paste0(save_plots, "roach_habitats_rsf_muddyfoot.png"), 
       plot = roach_habitat_rsf_plot, 
       device = 'png',
       width = 10, 
       height = 10,
       units = 'cm',
       dpi = 300)

# Individual variation plot
roach_all_diag <- rbind(roach_control_diag, roach_mix_diag)

roach_individual_plot <- ggplot(roach_all_diag, 
                                aes(x = treatment, y = coefficient, color = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, fill = NA, linewidth = 0.8) +
  scale_color_manual(values = c("Control" = "#0072B2", "Exposed" = "#D55E00")) +
  labs(y = "Selection coefficient", 
       title = "Individual variation in roach habitat selection",
       subtitle = paste0("Control: n=", nrow(roach_control_diag), 
                         ", Exposed: n=", nrow(roach_mix_diag))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 'bold', size = 14),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(roach_individual_plot)

ggsave(file = paste0(save_plots, "roach_individual_variation_rsf.png"), 
       plot = roach_individual_plot, 
       device = 'png',
       width = 12, 
       height = 10,
       units = 'cm',
       dpi = 300)

#-------------------------------------------------------------------------------#

#> 5.3 Pike - IMPROVED ####

cat("\n\n=== PIKE ANALYSIS ===\n")

# Load RSF results
rsf_pike_control_list <- readRDS(paste0(rsf_path, "muddyfoot_pike/rsf_pike_control_list.rds"))
rsf_pike_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_pike/rsf_pike_mix_list.rds"))
pike_control_diag <- readRDS(paste0(rsf_path, "muddyfoot_pike/rsf_pike_control_diagnostics.rds"))
pike_mix_diag <- readRDS(paste0(rsf_path, "muddyfoot_pike/rsf_pike_mix_diagnostics.rds"))

# Calculate means
rsf_pike_control_mean <- mean(rsf_pike_control_list)
rsf_pike_exposed_mean <- mean(rsf_pike_exposed_list)

# Extract coefficients
rsf_coef_control_pike <- as.data.frame(t(summary(rsf_pike_control_mean)$CI[1,]))
rsf_coef_exposed_pike <- as.data.frame(t(summary(rsf_pike_exposed_mean)$CI[1,]))

rsf_coef_control_pike$treatment <- "Control"
rsf_coef_exposed_pike$treatment <- "Exposed"

pike_habitat_rsf_coefs <- rbind(rsf_coef_control_pike, rsf_coef_exposed_pike)

# Calculate SEs
control_se <- (rsf_coef_control_pike$high - rsf_coef_control_pike$low) / (2 * 1.96)
exposed_se <- (rsf_coef_exposed_pike$high - rsf_coef_exposed_pike$low) / (2 * 1.96)

# Statistical comparison
pike_comparison <- compare_treatments(
  rsf_coef_control_pike$est, control_se,
  rsf_coef_exposed_pike$est, exposed_se
)

# Create summary
pike_summary <- create_rsf_summary_table(
  "Pike",
  pike_control_diag,
  pike_mix_diag,
  rsf_pike_control_mean,
  rsf_pike_exposed_mean,
  pike_comparison
)

saveRDS(pike_summary, paste0(rsf_path, "muddyfoot_pike/pike_rsf_summary.rds"))

# Population-level plot
pike_habitat_rsf_plot <- ggplot(pike_habitat_rsf_coefs, 
                                aes(x = treatment, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = low, ymax = high), 
                width = 0.1, size = 1, color = "black") +
  geom_point(aes(shape = treatment, fill = treatment), 
             size = 4, color = "black", stroke = 1) +
  scale_shape_manual(values = c(21, 21)) +
  scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +
  coord_cartesian(ylim = c(min(pike_habitat_rsf_coefs$low) - 0.5, 
                           max(pike_habitat_rsf_coefs$high) + 0.5)) +
  labs(y = "Selection coefficient",
       title = "Pike Habitat Selection",
       subtitle = paste0("Δ = ", round(pike_comparison$difference, 2),
                         ", p = ", format.pval(pike_comparison$p_value, digits = 2))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

print(pike_habitat_rsf_plot)

ggsave(file = paste0(save_plots, "pike_habitats_rsf_muddyfoot.png"), 
       plot = pike_habitat_rsf_plot, 
       device = 'png',
       width = 10, 
       height = 10,
       units = 'cm',
       dpi = 300)

# Individual variation plot
pike_all_diag <- rbind(pike_control_diag, pike_mix_diag)

pike_individual_plot <- ggplot(pike_all_diag, 
                               aes(x = treatment, y = coefficient, color = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, fill = NA, linewidth = 0.8) +
  scale_color_manual(values = c("Control" = "#0072B2", "Exposed" = "#D55E00")) +
  labs(y = "Selection coefficient", 
       title = "Individual variation in pike habitat selection",
       subtitle = paste0("Control: n=", nrow(pike_control_diag), 
                         ", Exposed: n=", nrow(pike_mix_diag))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 'bold', size = 14),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(pike_individual_plot)

ggsave(file = paste0(save_plots, "pike_individual_variation_rsf.png"), 
       plot = pike_individual_plot, 
       device = 'png',
       width = 12, 
       height = 10,
       units = 'cm',
       dpi = 300)

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------#
# 6. Create multi-species comparison ####
#-----------------------------------------#

cat("\n\n=== MULTI-SPECIES COMPARISON ===\n")

# Compile all results
all_pop_results <- rbind(
  data.frame(species = "Perch", perch_habitat_rsf_coefs),
  data.frame(species = "Roach", roach_habitat_rsf_coefs),
  data.frame(species = "Pike", pike_habitat_rsf_coefs)
)

# Multi-species plot
multi_species_plot <- ggplot(all_pop_results, 
                             aes(x = species, y = est, shape = treatment, fill = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = low, ymax = high), 
                position = position_dodge(width = 0.5),
                width = 0.2, size = 1, color = "black") +
  geom_point(position = position_dodge(width = 0.5),
             size = 4, color = "black", stroke = 1) +
  scale_shape_manual(values = c(21, 21)) +
  scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +
  labs(y = "Selection coefficient",
       title = "Habitat Selection Across Species",
       fill = "Treatment",
       shape = "Treatment") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 'bold', size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = 'bold', size = 12)
  )

print(multi_species_plot)

ggsave(file = paste0(save_plots, "multi_species_habitat_rsf.png"), 
       plot = multi_species_plot, 
       device = 'png',
       width = 16, 
       height = 12,
       units = 'cm',
       dpi = 300)

# Create summary table
comparison_summary <- data.frame(
  Species = c("Perch", "Roach", "Pike"),
  N_control = c(nrow(perch_control_diag), nrow(roach_control_diag), nrow(pike_control_diag)),
  N_exposed = c(nrow(perch_mix_diag), nrow(roach_mix_diag), nrow(pike_mix_diag)),
  Control_est = c(rsf_coef_control_perch$est, rsf_coef_control_roach$est, rsf_coef_control_pike$est),
  Exposed_est = c(rsf_coef_exposed_perch$est, rsf_coef_exposed_roach$est, rsf_coef_exposed_pike$est),
  Difference = c(perch_comparison$difference, roach_comparison$difference, pike_comparison$difference),
  P_value = c(perch_comparison$p_value, roach_comparison$p_value, pike_comparison$p_value),
  Significant = c(perch_comparison$significant, roach_comparison$significant, pike_comparison$significant)
)

cat("\nMulti-species comparison summary:\n")
print(comparison_summary, row.names = FALSE)

saveRDS(comparison_summary, paste0(rsf_path, "multi_species_comparison_summary.rds"))
write.csv(comparison_summary, paste0(rsf_path, "multi_species_comparison_summary.csv"), 
          row.names = FALSE)

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#

cat("\n\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")
cat("Summary of results:\n")
cat("- All models fitted successfully\n")
cat("- Diagnostic plots created\n")
cat("- Statistical comparisons completed\n")
cat("- Results saved to:", rsf_path, "\n")
cat("- Plots saved to:", save_plots, "\n\n")

# Print final summary
cat("Key findings:\n\n")
cat("PERCH:\n")
cat("  Treatment effect: p =", format.pval(perch_comparison$p_value, digits = 3), "\n")
cat("  Control coefficient:", round(rsf_coef_control_perch$est, 3), "\n")
cat("  Exposed coefficient:", round(rsf_coef_exposed_perch$est, 3), "\n\n")

cat("ROACH:\n")
cat("  Treatment effect: p =", format.pval(roach_comparison$p_value, digits = 3), "\n")
cat("  Control coefficient:", round(rsf_coef_control_roach$est, 3), "\n")
cat("  Exposed coefficient:", round(rsf_coef_exposed_roach$est, 3), "\n\n")

cat("PIKE:\n")
cat("  Treatment effect: p =", format.pval(pike_comparison$p_value, digits = 3), "\n")
cat("  Control coefficient:", round(rsf_coef_control_pike$est, 3), "\n")
cat("  Exposed coefficient:", round(rsf_coef_exposed_pike$est, 3), "\n\n")

cat("========================================\n\n")


























# ---------------------------------------------------------#
# RESOURCE SELECTION FOR ARTIFICIAL HABITATS - MUDDYFOOT ###
# ---------------------------------------------------------#

### LIBRARIES ###
library(ctmm)
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(data.table)
library(parallel)    # Parallel computing support
library(foreach)     # Looping construct for parallel execution
library(doParallel)  # Parallel backend for foreach loops

### DIRECTORIES ###
polygon_path = "./data/lake_coords/"
ctmm_path = "./data/ctmm_fits/"
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
telem_path <- "./data/telem_obj/muddyfoot/" 
rec_data_path = "./data/lake_coords/reciever_and_habitat_locations/"
enc_path <- "./data/encounters/muddyfoot/"                  # Directory for encounter data
akde_path <- "./data/akdes/"                      # Directory for AKDE (Autocorrelated Kernel Density Estimation) outputs
rsf_path <- "./data/rsfs/habitats/"
save_plots <- "./plots/muddyfoot/"  # Directory for saving utilization distribution plots

### LOAD DATA ###
#receiver and habitat locations
mud_rec_locs_kml <- paste0(rec_data_path, "muddyfoot_rec_hab_locations.kml")
mud_rec_locs <- st_read(mud_rec_locs_kml)[1:5,]
mud_hab_locs <- st_read(mud_rec_locs_kml)[6:7,]#receiver and habitat locations

#fish telemetry data
pike_muddyfoot_tel <- readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

#fish ctmms 
pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds"))

#fish akdes
pike_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/akde_cg/pike_akdes_cg_list.rds"))
perch_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/akde_cg/perch_akdes_cg_list.rds"))
roach_akdes_cg_list <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/akde_cg/roach_akdes_cg_list.rds"))

#load pkdes if necessary
perch_control_PKDE <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/population_akde/perch_control_PKDE.rds"))
perch_mix_PKDE <- readRDS(paste0(akde_path, "muddyfoot_perch_akdes/population_akde/perch_mix_PKDE.rds"))
roach_control_PKDE <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/population_akde/roach_control_PKDE.rds"))
roach_mix_PKDE <- readRDS(paste0(akde_path, "muddyfoot_roach_akdes/population_akde/roach_mix_PKDE.rds"))
pike_control_PKDE <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/population_akde/pike_control_PKDE.rds"))
pike_mix_PKDE <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/population_akde/pike_mix_PKDE.rds"))
pike_total_PKDE <- readRDS(paste0(akde_path, "muddyfoot_pike_akdes/population_akde/pike_total_PKDE.rds"))


#Load lake polygon
muddyfoot_polygon <- st_read(paste0(polygon_path, "lake_muddyfoot_polygon.gpkg"))

#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#---------------------------------------------#
# 1. Plot population akdes with habitats ####
#---------------------------------------------#

# Function to generate the plot with habitats
generate_ud_plot_whabitats <- function(pkde_data, bbox, hab_locs) {
  ud_raster <- raster(pkde_data)
  masked_ud_raster <- mask(ud_raster, bbox)
  ud_df <- as.data.frame(as(masked_ud_raster, "SpatialPixelsDataFrame"))
  colnames(ud_df) <- c("value", "x", "y")
  
  ggplot() +
    geom_sf(data = bbox, color = "black") +
    geom_tile(data = ud_df, aes(x = x, y = y, fill = value), alpha = 0.6) +
    geom_sf(data = hab_locs, color = "green", size = 2, fill = NA, shape = 3, stroke = 1) + 
    scale_fill_viridis_c(na.value = 'transparent', option = 'magma') +
    coord_sf() +
    theme_classic() +
    labs(fill = "Utilization Distribution", x = "", y = '')+
    theme(legend.position = "none",
          axis.text = element_text(size = 8, color = 'black'))
}

# Function to generate the plot without habitats
generate_ud_plot <- function(pkde_data, bbox, hab_locs) {
  ud_raster <- raster(pkde_data)
  masked_ud_raster <- mask(ud_raster, bbox)
  ud_df <- as.data.frame(as(masked_ud_raster, "SpatialPixelsDataFrame"))
  colnames(ud_df) <- c("value", "x", "y")
  
  ggplot() +
    geom_sf(data = bbox, color = "black") +
    geom_tile(data = ud_df, aes(x = x, y = y, fill = value), alpha = 0.6) +
    scale_fill_viridis_c(na.value = 'transparent', option = 'magma') +
    coord_sf() +
    theme_classic() +
    labs(fill = "Utilization Distribution", x = "", y = '')+
    theme(legend.position = "none",
          axis.text = element_text(size = 8, color = 'black'))
}

# Transform bbox
muddyfoot_bbox_perch <- st_transform(muddyfoot_polygon, crs(raster(perch_control_PKDE)))
muddyfoot_bbox_roach <- st_transform(muddyfoot_polygon, crs(raster(roach_control_PKDE)))
muddyfoot_bbox_pike <- st_transform(muddyfoot_polygon, crs(raster(pike_total_PKDE)))

# Generate the plots

#> 1.1. Perch ####

## CONTROL ##
perch_control_plot <- generate_ud_plot(perch_control_PKDE, muddyfoot_bbox_perch, mud_hab_locs)
perch_control_plot_habitats <- generate_ud_plot_whabitats(perch_control_PKDE, muddyfoot_bbox_perch, mud_hab_locs)
ggsave(file = paste0(save_plots, "perch_control_UD_habitats_muddyfoot.png"), 
       plot = perch_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## EXPOSED ##
perch_exposed_plot <- generate_ud_plot(perch_mix_PKDE, muddyfoot_bbox_perch, mud_hab_locs)
perch_exposed_plot_habitats <- generate_ud_plot_whabitats(perch_mix_PKDE, muddyfoot_bbox_perch, mud_hab_locs)
ggsave(file = paste0(save_plots, "perch_exposed_UD_habitats_muddyfoot.png"), 
       plot = perch_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)



#> 1.2. Roach ####
## CONTROL ##
roach_control_plot <- generate_ud_plot(roach_control_PKDE, muddyfoot_bbox_roach, mud_hab_locs)
roach_control_plot_habitats <- generate_ud_plot_whabitats(roach_control_PKDE, muddyfoot_bbox_roach, mud_hab_locs)
ggsave(file = paste0(save_plots, "roach_control_UD_habitats_muddyfoot.png"), 
       plot = roach_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## EXPOSED ##
roach_exposed_plot <- generate_ud_plot(roach_mix_PKDE, muddyfoot_bbox_roach, mud_hab_locs)
roach_exposed_plot_habitats <- generate_ud_plot_whabitats(roach_mix_PKDE, muddyfoot_bbox_roach, mud_hab_locs)
ggsave(file = paste0(save_plots, "roach_exposed_UD_habitats_muddyfoot.png"), 
       plot = roach_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)



#> 1.3. Pike ####

## CONTROL ##
pike_control_plot <- generate_ud_plot(pike_control_PKDE, muddyfoot_bbox_pike, mud_hab_locs)
pike_control_plot_habitats <- generate_ud_plot_whabitats(pike_control_PKDE, muddyfoot_bbox_pike, mud_hab_locs)
ggsave(file = paste0(save_plots, "pike_control_UD_habitats_muddyfoot.png"), 
       plot = pike_control_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## EXPOSED ##
pike_exposed_plot <- generate_ud_plot(pike_mix_PKDE, muddyfoot_bbox_pike, mud_hab_locs)
pike_exposed_plot_habitats <- generate_ud_plot_whabitats(pike_mix_PKDE, muddyfoot_bbox_pike, mud_hab_locs)
ggsave(file = paste0(save_plots, "pike_exposed_UD_habitats_muddyfoot.png"), 
       plot = pike_exposed_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)


## TOTAL ##
pike_total_plot <- generate_ud_plot(pike_total_PKDE, muddyfoot_bbox_pike, mud_hab_locs)
pike_total_plot_habitats <- generate_ud_plot_whabitats(pike_total_PKDE, muddyfoot_bbox_pike, mud_hab_locs)
ggsave(file = paste0(save_plots, "pike_total_UD_habitats_muddyfoot.png"), 
       plot = pike_total_plot_habitats, 
       device = 'png',
       width = 9, 
       height = 6.5,
       units = 'cm',
       dpi = 300)

#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#-----------------------------------------#
# 2. Create a lake and habitat raster ####
#-----------------------------------------#

#Final raster is used in RSF models

#Transform the CRS to UTM Zone 33N (EPSG:32633)
polyProj <- st_transform(muddyfoot_polygon, crs = "EPSG:32633")

# Create a raster object with the desired resolution and CRS (UTM Zone 33N)
resolution <- 0.5
lake_raster <- rast(ext(polyProj), 
                    res = resolution, 
                    crs = "EPSG:32633")

# Reproject the raster to WGS 84 (EPSG:4326)
lake_raster <- terra::project(lake_raster, "EPSG:4326")

# Ensure that the muddyfoot_polygon is also in WGS 84 (EPSG:4326)
muddyfoot_polygon <- st_transform(muddyfoot_polygon, crs = "EPSG:4326")

# Rasterize the lake polygon (masking cells outside the lake)
lake_raster <- rasterize(muddyfoot_polygon, lake_raster, background = NA)
plot(lake_raster, main = "Lake Raster")

#Plot habitats onto lake raster

habitat_patches <- vect(mud_hab_locs)

# Assign a value of 1 to the habitat patches and leave other areas as 0 or NA
habitat_raster <- rasterize(habitat_patches, lake_raster, field = 1, background = 0)
plot(habitat_raster, main = "Habitat Raster")

#To ensure that only the lake area is included in the habitat raster (and exclude areas outside the lake), 
#mask the raster using the lake polygon:
# Mask the habitat raster to the lake polygon
habitat_raster <- mask(habitat_raster, muddyfoot_polygon)
plot(habitat_raster, main = "Habitat Raster Masked by Lake")

#convert to raster
#Need to convert it from terra to raster to work with rsf.fit()
habitat_raster <- raster::raster(habitat_raster)

#save raster
writeRaster(habitat_raster, "./data/lake_coords/muddyfoot_habitat_raster.grd", overwrite = TRUE)
#open raster
habitat_raster <- rast("./data/lake_coords/muddyfoot_habitat_raster.grd")


#---------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------#
#3. Split telemetry objects by treatment ####
#---------------------------------------------#

# Separating into 'control' and 'mix'

#> 3.1 Perch ########################

# Check that treatments are in the right order
unique_ids <- names(perch_muddyfoot_tel)

# Create a data frame to store ID and unique Treatment information
result_table <- do.call(rbind, lapply(unique_ids, function(id) {
  # Extract the unique Treatment value for the current ID
  treatment_info <- unique(perch_muddyfoot_tel[[id]]$Treatment)
  
  # Ensure only one unique Treatment value is captured
  if (length(treatment_info) > 1) {
    warning(paste("Multiple treatments found for ID:", id, 
                  "using the first value. Treatments:", paste(treatment_info, collapse = ", ")))
    treatment_info <- treatment_info[1]
  }
  
  # Create a data frame for the ID and its treatment
  data.frame(ID = id, Treatment = treatment_info)
}))

# View the resulting table
print(result_table)

#use perch as an example
#telemetry objects
perch_control_tel <- perch_muddyfoot_tel[1:15]
perch_mix_tel <- perch_muddyfoot_tel[16:30]

#akdes
perch_control_akdes <- perch_akdes_cg_list[1:15]
perch_mix_akdes <- perch_akdes_cg_list[16:30]


#> 3.2 Roach ########################

# Check that treatments are in the right order
unique_ids <- names(roach_muddyfoot_tel)

# Create a data frame to store ID and unique Treatment information
result_table <- do.call(rbind, lapply(unique_ids, function(id) {
  # Extract the unique Treatment value for the current ID
  treatment_info <- unique(roach_muddyfoot_tel[[id]]$Treatment)
  
  # Ensure only one unique Treatment value is captured
  if (length(treatment_info) > 1) {
    warning(paste("Multiple treatments found for ID:", id, 
                  "using the first value. Treatments:", paste(treatment_info, collapse = ", ")))
    treatment_info <- treatment_info[1]
  }
  
  # Create a data frame for the ID and its treatment
  data.frame(ID = id, Treatment = treatment_info)
}))

print(result_table)



#telemetry objects
roach_control_tel <- roach_muddyfoot_tel[1:14]
roach_mix_tel <- roach_muddyfoot_tel[15:28]

#akdes
roach_control_akdes <- roach_akdes_cg_list[1:14]
roach_mix_akdes <- roach_akdes_cg_list[15:28]


#> 3.3 Pike ########################

# Separate telemetry objects
pike_control_tel <- pike_muddyfoot_tel[1:3]   # Control group
pike_mix_tel <- pike_muddyfoot_tel[4:6]       # Mixed group

# Separate AKDEs for the two groups
pike_control_akdes <- pike_akdes_cg_list[1:3]
pike_mix_akdes <- pike_akdes_cg_list[4:6]

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------#
# 4. Fit RSF models ####
#-----------------------#

#> 4.1 Perch ########################

### CONTROL ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_perch_control_list <- list()

# Assuming perch_control_tel and perch_control_akdes are lists of corresponding elements
rsf_perch_control_list <- foreach(i = seq_along(perch_control_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_control_model <- rsf.fit(
    perch_control_tel[[i]], 
    perch_control_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_control_model, 
          file = paste0(rsf_path, "muddyfoot_perch/", names(perch_control_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_perch_control_list) <- names(perch_control_tel)

#check
summary(rsf_perch_control_list$F59702)
summary(rsf_perch_control_list$F59697)

saveRDS(rsf_perch_control_list, paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_list.rds"))


### EXPOSED ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_perch_mix_list <- list()

rsf_perch_mix_list <- foreach(i = seq_along(perch_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    perch_mix_tel[[i]], 
    perch_mix_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "muddyfoot_perch/", names(perch_mix_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_perch_mix_list) <- names(perch_mix_tel)

#check
summary(rsf_perch_mix_list$F59714)
summary(rsf_perch_mix_list$F59727)

saveRDS(rsf_perch_mix_list, paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_list.rds"))

#-----------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------#

#> 4.2 Roach ########################

### CONTROL ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_roach_control_list <- list()

# Assuming roach_control_tel and roach_control_akdes are lists of corresponding elements
rsf_roach_control_list <- foreach(i = seq_along(roach_control_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_control_model <- rsf.fit(
    roach_control_tel[[i]], 
    roach_control_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_control_model, 
          file = paste0(rsf_path, "muddyfoot_roach/", names(roach_control_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_roach_control_list) <- names(roach_control_tel)

#check
summary(rsf_roach_control_list$F59683)
summary(rsf_roach_control_list$F59684)

saveRDS(rsf_roach_control_list, paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_list.rds"))



### EXPOSED ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_roach_mix_list <- list()

# Assuming roach_control_tel and roach_control_akdes are lists of corresponding elements
rsf_roach_mix_list <- foreach(i = seq_along(roach_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    roach_mix_tel[[i]], 
    roach_mix_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "muddyfoot_roach/", names(roach_mix_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_roach_mix_list) <- names(roach_mix_tel)

saveRDS(rsf_roach_mix_list, paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_list.rds"))

#-----------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------#

#> 4.3. Pike ########################

### CONTROL ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_pike_control_list <- list()

# Assuming pike_control_tel and pike_control_akdes are lists of corresponding elements
rsf_pike_control_list <- foreach(i = seq_along(pike_control_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_control_model <- rsf.fit(
    pike_control_tel[[i]], 
    pike_control_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_control_model, 
          file = paste0(rsf_path, "muddyfoot_pike/", names(pike_control_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_control_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_pike_control_list) <- names(pike_control_tel)


saveRDS(rsf_pike_control_list, paste0(rsf_path, "muddyfoot_pike/rsf_pike_control_list.rds"))


### EXPOSED ###

cl <- makeCluster(3)
doParallel::registerDoParallel(cl)
rsf_pike_mix_list <- list()

# Assuming pike_control_tel and pike_control_akdes are lists of corresponding elements
rsf_pike_mix_list <- foreach(i = seq_along(pike_mix_tel), .packages = "ctmm") %dopar% {
  
  # Run the rsf.select model
  rsf_mix_model <- rsf.fit(
    pike_mix_tel[[i]], 
    pike_mix_akdes[[i]], 
    R=list(habitat1=habitat_raster)
  )
  
  # Save the model to the 'rsfs' folder with an appropriate name
  saveRDS(rsf_mix_model, 
          file = paste0(rsf_path, "muddyfoot_pike/", names(pike_mix_tel)[i], "_habitat_rsf.rds"))
  
  # Return the model in case you want to store it in a list
  rsf_mix_model
}

stopCluster(cl)

# Assign individual IDs to the AKDE list
names(rsf_pike_mix_list) <- names(pike_mix_tel)

saveRDS(rsf_pike_mix_list, paste0(rsf_path, "muddyfoot_pike/rsf_pike_mix_list.rds"))

#-------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------#

#-----------------------------------------#
# 5. Explore RSF model results ####
#-----------------------------------------#

#> 5.1. Perch ####

# Load the RSF (Resource Selection Function) results for control and exposed perch
rsf_perch_control_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/rsf_perch_control_list.rds"))
rsf_perch_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_perch/rsf_perch_mix_list.rds"))

# Calculate the mean of the RSF lists (assuming that each list is a collection of model objects)
# This assumes 'mean()' can handle these model objects, modify if you need more specific behavior
rsf_perch_control_mean <- mean(rsf_perch_control_list)
rsf_perch_exposed_mean <- mean(rsf_perch_exposed_list)

# Extract coefficients' confidence intervals from the RSF summaries
# Convert the first confidence interval (CI) row to a data frame for each treatment
rsf_coef_control_perch <- as.data.frame(t(summary(rsf_perch_control_mean)$CI[1,]))
rsf_coef_exposed_perch <- as.data.frame(t(summary(rsf_perch_exposed_mean)$CI[1,]))

# Add a treatment column to each dataset to distinguish between Control and Exposed in the final plot
rsf_coef_control_perch$treatment <- "Control"
rsf_coef_exposed_perch$treatment <- "Exposed"

# Combine the coefficients for both treatments into one data frame
perch_habitat_rsf_coefs <- rbind(rsf_coef_control_perch, rsf_coef_exposed_perch)

# Create the ggplot visualization of RSF coefficient estimates for perch habitats
(perch_habitat_rsf_plot <- 
    ggplot(perch_habitat_rsf_coefs, 
           aes(x = treatment, y = est)) +  # Aesthetic mapping: treatment on x-axis, estimate (est) on y-axis
    geom_errorbar(aes(ymin = low, ymax = high), 
                  width = 0.1,  # Width of error bars
                  size = 1,       # Thickness of error bars
                  color = "black") +  # Error bars in black
    geom_point(aes(shape = treatment, fill = treatment), 
               size = 4,  # Point size
               color = "black") +  # Outline of points in black
    scale_shape_manual(values = c(21, 21)) +  # Use the same shape (circle) for both treatments
    scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +  # White fill for Control, black fill for Exposed
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add a dashed horizontal line at y=0
    coord_cartesian(ylim = c(-1, 6)) +  # Set y-axis limits from -1 to 6
    labs(y = "Selection coefficient") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)


ggsave(file = paste0(save_ud_plots, "perch_habitats_rsf_muddyfoot.png"), 
       plot = perch_habitat_rsf_plot, 
       device = 'png',
       width = 8, 
       height = 8,
       units = 'cm',
       dpi = 300)

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#


#> 5.2 Roach ####

# Load the RSF (Resource Selection Function) results for control and exposed roach
rsf_roach_control_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/rsf_roach_control_list.rds"))
rsf_roach_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_roach/rsf_roach_mix_list.rds"))

# Calculate the mean of the RSF lists (assuming that each list is a collection of model objects)
# This assumes 'mean()' can handle these model objects, modify if you need more specific behavior
rsf_roach_control_mean <- mean(rsf_roach_control_list)
rsf_roach_exposed_mean <- mean(rsf_roach_exposed_list)

# Extract coefficients' confidence intervals from the RSF summaries
# Convert the first confidence interval (CI) row to a data frame for each treatment
rsf_coef_control_roach <- as.data.frame(t(summary(rsf_roach_control_mean)$CI[1,]))
rsf_coef_exposed_roach <- as.data.frame(t(summary(rsf_roach_exposed_mean)$CI[1,]))

# Add a treatment column to each dataset to distinguish between Control and Exposed in the final plot
rsf_coef_control_roach$treatment <- "Control"
rsf_coef_exposed_roach$treatment <- "Exposed"

# Combine the coefficients for both treatments into one data frame
roach_habitat_rsf_coefs <- rbind(rsf_coef_control_roach, rsf_coef_exposed_roach)

# Create the ggplot visualization of RSF coefficient estimates for roach habitats
(roach_habitat_rsf_plot <- 
    ggplot(roach_habitat_rsf_coefs, 
           aes(x = treatment, y = est)) +  # Aesthetic mapping: treatment on x-axis, estimate (est) on y-axis
    geom_errorbar(aes(ymin = low, ymax = high), 
                  width = 0.1,  # Width of error bars
                  size = 1,       # Thickness of error bars
                  color = "black") +  # Error bars in black
    geom_point(aes(shape = treatment, fill = treatment), 
               size = 4,  # Point size
               color = "black") +  # Outline of points in black
    scale_shape_manual(values = c(21, 21)) +  # Use the same shape (circle) for both treatments
    scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +  # White fill for Control, black fill for Exposed
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add a dashed horizontal line at y=0
    coord_cartesian(ylim = c(-1, 6)) +  # Set y-axis limits from -1 to 6
    labs(y = "Selection coefficient") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)


ggsave(file = paste0(save_ud_plots, "roach_habitats_rsf_muddyfoot.png"), 
       plot = roach_habitat_rsf_plot, 
       device = 'png',
       width = 8, 
       height = 8,
       units = 'cm',
       dpi = 300)


#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

#> 5.3 Pike ####

# Load the RSF (Resource Selection Function) results for control and exposed pike
rsf_pike_control_list <- readRDS(paste0(rsf_path, "muddyfoot_pike/rsf_pike_control_list.rds"))
rsf_pike_exposed_list <- readRDS(paste0(rsf_path, "muddyfoot_pike/rsf_pike_mix_list.rds"))

# Calculate the mean of the RSF lists (assuming that each list is a collection of model objects)
# This assumes 'mean()' can handle these model objects, modify if you need more specific behavior
rsf_pike_control_mean <- mean(rsf_pike_control_list)
rsf_pike_exposed_mean <- mean(rsf_pike_exposed_list)

# Extract coefficients' confidence intervals from the RSF summaries
# Convert the first confidence interval (CI) row to a data frame for each treatment
rsf_coef_control_pike <- as.data.frame(t(summary(rsf_pike_control_mean)$CI[1,]))
rsf_coef_exposed_pike <- as.data.frame(t(summary(rsf_pike_exposed_mean)$CI[1,]))

# Add a treatment column to each dataset to distinguish between Control and Exposed in the final plot
rsf_coef_control_pike$treatment <- "Control"
rsf_coef_exposed_pike$treatment <- "Exposed"

# Combine the coefficients for both treatments into one data frame
pike_habitat_rsf_coefs <- rbind(rsf_coef_control_pike, rsf_coef_exposed_pike)

# Create the ggplot visualization of RSF coefficient estimates for pike habitats
(pike_habitat_rsf_plot <- 
    ggplot(pike_habitat_rsf_coefs, 
           aes(x = treatment, y = est)) +  # Aesthetic mapping: treatment on x-axis, estimate (est) on y-axis
    geom_errorbar(aes(ymin = low, ymax = high), 
                  width = 0.1,  # Width of error bars
                  size = 1,       # Thickness of error bars
                  color = "black") +  # Error bars in black
    geom_point(aes(shape = treatment, fill = treatment), 
               size = 4,  # Point size
               color = "black") +  # Outline of points in black
    scale_shape_manual(values = c(21, 21)) +  # Use the same shape (circle) for both treatments
    scale_fill_manual(values = c("Control" = "white", "Exposed" = "black")) +  # White fill for Control, black fill for Exposed
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add a dashed horizontal line at y=0
    coord_cartesian(ylim = c(-1, 6)) +  # Set y-axis limits from -1 to 6
    labs(y = "Selection coefficient") +  # Label for the y-axis
    theme_classic() +  # Use a clean classic theme for the plot
    theme(legend.position = "none",  # Remove the legend as it's not necessary
          axis.title.x = element_blank(),  # Remove the x-axis title for simplicity
          axis.title.y = element_text(face = 'bold', size = 16, margin = margin(r = 10)),  # Bold y-axis title with larger font
          axis.text = element_text(size = 12, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1))  # Set font size for axis labels
)


ggsave(file = paste0(save_ud_plots, "pike_habitats_rsf_muddyfoot.png"), 
       plot = pike_habitat_rsf_plot, 
       device = 'png',
       width = 8, 
       height = 8,
       units = 'cm',
       dpi = 300)

#------------------------------------------------------------------------------------------------------------------------------#

