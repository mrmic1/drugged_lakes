# Resource Selection Function Analysis for Fish Tracking Data
# Analyzes habitat selection patterns for multiple fish species
# Author: [Your name]
# Date: 2026-01-21

# Load packages ----------------------------------------------------------------
library(amt)      # Animal movement tools
library(sf)       # Spatial data handling
library(terra)    # Raster operations
library(dplyr)    # Data manipulation
library(purrr)    # Functional programming
library(ggplot2)  # Visualization
library(broom)    # Model tidying
library(tidyr)    # Data reshaping

# Set paths --------------------------------------------------------------------
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
rec_data_path <- "./data/lake_params/reciever_and_habitat_locations/"
polygon_path <- "./data/lake_params/polygons/"
rsf_output_path <- "./output/muddyfoot/rsf/"


#==============================================================================
# PART 1: PREPARE TRACKING DATA
#==============================================================================

# Load and prepare tracking data ----------------------------------------------
muddyfoot_filt_data <- readRDS(paste0(filtered_data_path, "06_muddyfoot_sub.rds"))

cat("\n=== DATA SUMMARY ===\n")
cat("Timestamp range:", format(range(muddyfoot_filt_data$timestamp)), "\n")
cat("Number of individuals:", n_distinct(muddyfoot_filt_data$individual_ID), "\n")
cat("Species:", paste(unique(muddyfoot_filt_data$species), collapse = ", "), "\n")
cat("Total observations:", nrow(muddyfoot_filt_data), "\n")

# Convert to track format and project to UTM 33N -------------------------------
muddyfoot_trk <- muddyfoot_filt_data |>
  transmute(
    x_ = Long,
    y_ = Lat,
    t_ = timestamp,
    id = individual_ID,
    species = species,
    treatment = treatment
  ) |>
  # Convert to spatial points and transform to UTM
  st_as_sf(coords = c("x_", "y_"), crs = 4326) |>
  st_transform(32633) |>
  mutate(
   x_ = st_coordinates(geometry)[, 1],
   y_ = st_coordinates(geometry)[, 2]
  ) |>
  st_drop_geometry()

# Set track class and CRS
class(muddyfoot_trk) <- c("track_xyt", "track_xy", "tbl_df", "tbl", "data.frame")
attr(muddyfoot_trk, "crs_") <- st_crs(32633)

# Clean up
rm(muddyfoot_filt_data)
gc()

# Resample to regular intervals ------------------------------------------------
cat("\n=== RESAMPLING TRACKS ===\n")

muddyfoot_trk_resampled <- muddyfoot_trk |>
  nest(data = -c(id, species, treatment)) |>
  mutate(
    n_original = map_int(data, nrow),
    
    # Resample to 30-second intervals
    data_resampled = map(data, ~ track_resample(
      .x, 
      rate = seconds(30),
      tolerance = seconds(30)
    )),
    
    # Filter out short bursts
    data_filtered = map(data_resampled, ~ filter_min_n_burst(.x, min_n = 3)),
    
    n_final = map_int(data_filtered, nrow),
    pct_retained = round((n_final / n_original) * 100, 1)
  )

# Summary
cat("Data retention:\n")
print(muddyfoot_trk_resampled |> 
        select(id, species, treatment, n_original, n_final, pct_retained))

#==============================================================================
# PART 2: PREPARE SPATIAL DATA
#==============================================================================

# Load spatial data -----------------------------------------------------------
cat("\n=== LOADING SPATIAL DATA ===\n")

mud_rec_locs_kml <- paste0(rec_data_path, "muddyfoot_rec_hab_locations.kml")
mud_hab_locs <- st_read(mud_rec_locs_kml, quiet = TRUE)[6:7, ]
muddyfoot_polygon <- st_read(paste0(polygon_path, "muddyfoot_polygon.gpkg"), quiet = TRUE)

# Transform to UTM 33N
mud_hab_locs <- st_transform(mud_hab_locs, crs = 32633)
muddyfoot_polygon <- st_transform(muddyfoot_polygon, crs = 32633)

# Buffer habitats to account for telemetry error ------------------------------
cat("Buffering habitat patches by 0.5m for telemetry error...\n")

mud_hab_locs_buffered <- st_buffer(mud_hab_locs, dist = 0.5)

# Create habitat raster -------------------------------------------------------
cat("\n=== CREATING HABITAT RASTER ===\n")

# Use finer resolution to capture the 0.5m buffer accurately
resolution <- 0.5  # Changed from 1m to 0.5m
lake_extent <- ext(muddyfoot_polygon)

# Create raster template
lake_raster <- rast(lake_extent, res = resolution, crs = "EPSG:32633")
lake_raster <- rasterize(vect(muddyfoot_polygon), lake_raster, field = 1, background = NA)

# Create binary habitat raster (buffered)
habitat_raster <- rasterize(
  vect(mud_hab_locs_buffered),
  lake_raster,
  field = 1,
  background = 0
)
habitat_raster <- mask(habitat_raster, vect(muddyfoot_polygon))

# Better visualization to show buffer vs original
par(mfrow = c(1, 2))

# Plot 1: Original vs Buffered polygons
plot(st_geometry(muddyfoot_polygon), main = "Original (blue) vs Buffered (red)", 
     col = "lightgray", border = "black")
plot(st_geometry(mud_hab_locs), col = "blue", add = TRUE, lwd = 2)
plot(st_geometry(mud_hab_locs_buffered), col = "red", add = TRUE, lwd = 2, border = "red")
legend("topright", legend = c("Original 1.5x1.5m", "Buffered +0.5m"), 
       col = c("blue", "red"), lwd = 2)

# Plot 2: Rasterized habitat
plot(habitat_raster, main = "Habitat Raster (0.5m resolution)", 
     col = c("white", "red"), legend = FALSE)
plot(st_geometry(muddyfoot_polygon), add = TRUE, lwd = 2)
plot(st_geometry(mud_hab_locs), col = "blue", add = TRUE, lwd = 2)
text(st_coordinates(st_centroid(mud_hab_locs))[,1], 
     st_coordinates(st_centroid(mud_hab_locs))[,2], 
     labels = c("Habitat 1", "Habitat 2"), cex = 0.8)

par(mfrow = c(1, 1))

# Verify buffer dimensions
cat("\nHabitat dimensions:\n")
original_bbox <- st_bbox(mud_hab_locs)
buffered_bbox <- st_bbox(mud_hab_locs_buffered)
cat("Original width: ~1.5m, Buffered width: ~2.5m (should be +1m total)\n")
cat("Original bbox:", paste(round(original_bbox, 2), collapse = ", "), "\n")
cat("Buffered bbox:", paste(round(buffered_bbox, 2), collapse = ", "), "\n")

# Count habitat cells
n_habitat_cells <- global(habitat_raster == 1, "sum", na.rm = TRUE)
habitat_area_m2 <- n_habitat_cells * (resolution^2)
cat("Habitat area from raster:", round(habitat_area_m2, 2), "m²\n")
cat("Expected buffered area: ~12-13 m² (2 patches × ~6-6.5 m² each)\n")

# Save raster
writeRaster(habitat_raster, paste0(rsf_output_path, "habitat_raster_buffered.tif"), 
            overwrite = TRUE)

# Save visualization
png(paste0(rsf_output_path, "habitat_buffer_check.png"), width = 1200, height = 600)
par(mfrow = c(1, 2))
plot(st_geometry(muddyfoot_polygon), main = "Original (blue) vs Buffered (red)", 
     col = "lightgray", border = "black")
plot(st_geometry(mud_hab_locs), col = "blue", add = TRUE, lwd = 2)
plot(st_geometry(mud_hab_locs_buffered), col = "red", add = TRUE, lwd = 2, border = "red")
legend("topright", legend = c("Original 1.5x1.5m", "Buffered +0.5m"), 
       col = c("blue", "red"), lwd = 2)

plot(habitat_raster, main = "Habitat Raster (0.5m resolution)", 
     col = c("white", "red"), legend = FALSE)
plot(st_geometry(muddyfoot_polygon), add = TRUE, lwd = 2)
plot(st_geometry(mud_hab_locs), col = "blue", add = TRUE, lwd = 2)
dev.off()
par(mfrow = c(1, 1))

#==============================================================================
# PART 3: RSF ANALYSIS
#==============================================================================

# Prepare location data -------------------------------------------------------
cat("\n=== PREPARING RSF DATA ===\n")

muddyfoot_trk_locs <- muddyfoot_trk_resampled |>
  select(id, species, treatment, data_filtered) |>
  unnest(cols = data_filtered) |>
  select(x_, y_, t_, id, species, treatment)

# Restore track class
class(muddyfoot_trk_locs) <- c("track_xyt", "track_xy", "tbl_df", "tbl", "data.frame")
attr(muddyfoot_trk_locs, "crs_") <- st_crs(32633)

cat("Total locations:", nrow(muddyfoot_trk_locs), "\n")
cat("Individuals:", n_distinct(muddyfoot_trk_locs$id), "\n")

# Generate random available points and extract habitat ------------------------
prepare_rsf_data <- function(trk_data, polygon, hab_raster) {
  
  # Restore track attributes
  class(trk_data) <- c("track_xyt", "track_xy", "tbl_df", "tbl", "data.frame")
  attr(trk_data, "crs_") <- st_crs(32633)
  
  # Generate random available points (10x the number of used points)
  rsf_data <- random_points(
    x = polygon,
    n = nrow(trk_data) * 10,
    presence = trk_data
  )
  
  # Extract habitat values
  coords_matrix <- as.matrix(rsf_data[, c("x_", "y_")])
  rsf_data$habitat <- terra::extract(hab_raster, coords_matrix)[, 1]
  
  return(as.data.frame(rsf_data))
}

# Nest by individual and prepare RSF data
muddyfoot_rsf_nested <- muddyfoot_trk_locs |>
  nest(data = -c(id, species, treatment)) |>
  mutate(
    rsf_data = map(
      data, 
      ~ prepare_rsf_data(.x, muddyfoot_polygon, habitat_raster)
    ),
    n_used = map_int(data, nrow),
    n_available = map_int(rsf_data, ~ sum(.x$case_ == FALSE))
  )

cat("RSF data prepared for", nrow(muddyfoot_rsf_nested), "individuals\n")
print(muddyfoot_rsf_nested |> select(id, species, treatment, n_used, n_available))

# Fit RSF models --------------------------------------------------------------
cat("\n=== FITTING RSF MODELS ===\n")

rsf_results <- muddyfoot_rsf_nested |>
  mutate(
    # Fit logistic regression: Pr(used) ~ habitat presence
    model = map(rsf_data, ~ fit_rsf(.x, case_ ~ habitat)),
    
    # Extract the actual GLM model object from fit_rsf output
    glm_model = map(model, ~ .$model),
    
    # Extract model summaries
    aic = map_dbl(glm_model, AIC),
    coefs = map(glm_model, ~ broom::tidy(.x))
  )

# Extract coefficients --------------------------------------------------------
rsf_coefs <- rsf_results |>
  select(id, species, treatment, coefs) |>
  unnest(coefs) |>
  filter(term == "habitat")  # Only keep habitat coefficient

cat("\n=== RSF RESULTS ===\n")
print(rsf_coefs |> select(id, species, treatment, estimate, std.error, p.value))

# Summarize by species and treatment ------------------------------------------
species_treatment_summary <- rsf_coefs |>
  group_by(species, treatment) |>
  summarise(
    n = n(),
    mean_estimate = round(mean(estimate), 4),
    se_mean = round(sd(estimate) / sqrt(n()), 4),
    n_positive = sum(estimate > 0),
    n_significant_positive = sum(estimate > 0 & p.value < 0.05),
    n_significant_negative = sum(estimate < 0 & p.value < 0.05),
    .groups = "drop"
  )

cat("\n=== SUMMARY BY SPECIES AND TREATMENT ===\n")
print(species_treatment_summary)

# Visualizations --------------------------------------------------------------
cat("\n=== CREATING VISUALIZATIONS ===\n")

# 1. Coefficient estimates by species and treatment
p1 <- rsf_coefs |>
  ggplot(aes(x = treatment, y = estimate, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ species) +
  labs(
    title = "Habitat Selection Coefficients by Species and Treatment",
    x = "Treatment",
    y = "Coefficient Estimate",
    caption = "Positive = selection for habitat, Negative = avoidance"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)
ggsave(paste0(rsf_output_path, "rsf_coefficients_by_species_treatment.png"), 
       p1, width = 10, height = 6)

# Save results ----------------------------------------------------------------
saveRDS(rsf_results, paste0(rsf_output_path, "rsf_models.rds"))
write.csv(rsf_coefs, paste0(rsf_output_path, "rsf_all_coefficients.csv"), 
          row.names = FALSE)
write.csv(species_treatment_summary, paste0(rsf_output_path, "rsf_species_treatment_summary.csv"), 
          row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", rsf_output_path, "\n")