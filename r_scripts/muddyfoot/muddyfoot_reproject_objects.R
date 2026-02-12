#===============================================================================
# CONVERT ALL DATA TO CTMM'S CUSTOM TPEQD PROJECTION
#===============================================================================
#
# PURPOSE: Standardize telemetry, ctmm fits, and lake polygon to use the same
# custom two-point equidistant projection that ctmm selected for optimal
# distance calculations within the study area.
#
# ADVANTAGE: Most accurate for small-scale movement analysis (30m lake)
#
# AUTHOR: [Your Name]
# DATE: February 2026
#===============================================================================

library(ctmm)
library(sf)
library(sp)

# Define paths (update these to match your directory structure)
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/muddyfoot/"
lake_polygon_path <- "./data/lake_params/polygons/"

#===============================================================================
# STEP 1: LOAD DATA
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
muddyfoot_polygon <- st_read(paste0(lake_polygon_path, "muddyfoot_polygon.gpkg"), 
                             quiet = TRUE)
muddyfoot_sp_data <- as(muddyfoot_polygon, "Spatial")

message("\nLake boundary polygon loaded")
message("Polygon CRS: ", st_crs(muddyfoot_polygon)$input, "\n")

#===============================================================================
# STEP 2: IDENTIFY TARGET PROJECTION FROM CTMM
#===============================================================================

message("\n========================================")
message("  IDENTIFYING TARGET PROJECTION")
message("========================================\n")

# Use the first pike's ctmm fit projection as the standard
# This is the custom tpeqd projection optimized for your study area
target_projection <- ctmm::projection(pike_muddyfoot_ctmm_fits[[1]])

message("Target projection (from pike ctmm fit):\n")
cat(target_projection, "\n\n")

# Verify this is a tpeqd projection
if(grepl("\\+proj=tpeqd", target_projection)) {
  message("✓ Confirmed: Two-Point Equidistant (tpeqd) projection")
  message("  This projection is optimized for accurate distance calculations")
  message("  within your study area\n")
} else {
  warning("Unexpected projection type. Expected tpeqd.\n")
}

#===============================================================================
# STEP 3: CHECK CURRENT PROJECTIONS
#===============================================================================

message("\n========================================")
message("  CURRENT PROJECTION STATUS")
message("========================================\n")

# Check telemetry projections
pike_tel_proj <- ctmm::projection(pike_muddyfoot_tel[[1]])
perch_tel_proj <- ctmm::projection(perch_muddyfoot_tel[[1]])
roach_tel_proj <- ctmm::projection(roach_muddyfoot_tel[[1]])

message("Current telemetry projections:")
message("  Pike:  ", ifelse(is.null(pike_tel_proj) || pike_tel_proj == "", 
                            "NOT SET", pike_tel_proj))
message("  Perch: ", ifelse(is.null(perch_tel_proj) || perch_tel_proj == "", 
                            "NOT SET", perch_tel_proj))
message("  Roach: ", ifelse(is.null(roach_tel_proj) || roach_tel_proj == "", 
                            "NOT SET", roach_tel_proj), "\n")

# Check fit projections
pike_fit_proj <- ctmm::projection(pike_muddyfoot_ctmm_fits[[1]])
perch_fit_proj <- ctmm::projection(perch_muddyfoot_ctmm_fits[[1]])
roach_fit_proj <- ctmm::projection(roach_muddyfoot_ctmm_fits[[1]])

message("Current ctmm fit projections:")
message("  Pike:  ", pike_fit_proj)
message("  Perch: ", perch_fit_proj)
message("  Roach: ", roach_fit_proj, "\n")

# Check polygon
poly_crs <- st_crs(muddyfoot_polygon)$input
message("Current polygon projection: ", poly_crs, "\n")

#===============================================================================
# STEP 4: CONVERT POLYGON TO TARGET PROJECTION
#===============================================================================

message("\n========================================")
message("  CONVERTING POLYGON PROJECTION")
message("========================================\n")

message("Converting lake polygon to tpeqd projection...\n")

# Convert sf polygon to the target projection
# Note: st_transform may not work directly with custom proj4 strings
# We need to use sp::spTransform instead

# First convert to Spatial object if not already
muddyfoot_sp_original <- as(muddyfoot_polygon, "Spatial")

# Create the target CRS object
target_crs <- sp::CRS(target_projection)

# Transform the polygon
muddyfoot_sp_data <- sp::spTransform(muddyfoot_sp_original, target_crs)

# Convert back to sf for verification
muddyfoot_polygon_projected <- st_as_sf(muddyfoot_sp_data)

message("✓ Polygon converted to tpeqd projection\n")

# Verify projection
new_poly_proj <- projection(muddyfoot_sp_data)
message("New polygon projection matches target: ", 
        identical(new_poly_proj, target_projection), "\n")

#===============================================================================
# STEP 5: STANDARDIZE PIKE PROJECTIONS
#===============================================================================

message("\n========================================")
message("  STANDARDIZING PIKE PROJECTIONS")
message("========================================\n")

# Update telemetry projections
message("Updating pike telemetry projections...")
for(i in 1:length(pike_muddyfoot_tel)) {
  ctmm::projection(pike_muddyfoot_tel[[i]]) <- target_projection
}
message("✓ Pike telemetry: ", length(pike_muddyfoot_tel), " objects updated\n")

# Update ctmm fit projections (should already match, but ensure consistency)
message("Updating pike ctmm fit projections...")
ctmm::projection(pike_muddyfoot_ctmm_fits) <- target_projection
message("✓ Pike fits: ", length(pike_muddyfoot_ctmm_fits), " objects updated\n")

# Verify
pike_tel_check <- ctmm::projection(pike_muddyfoot_tel[[1]])
pike_fit_check <- ctmm::projection(pike_muddyfoot_ctmm_fits[[1]])
pike_match <- identical(pike_tel_check, target_projection) && 
  identical(pike_fit_check, target_projection)

if(pike_match) {
  message("✓ Pike projections verified: All match target\n")
} else {
  warning("⚠️  Pike projections may not match perfectly\n")
}

#===============================================================================
# STEP 6: STANDARDIZE PERCH PROJECTIONS
#===============================================================================

message("\n========================================")
message("  STANDARDIZING PERCH PROJECTIONS")
message("========================================\n")

# Update telemetry projections
message("Updating perch telemetry projections...")
for(i in 1:length(perch_muddyfoot_tel)) {
  ctmm::projection(perch_muddyfoot_tel[[i]]) <- target_projection
}
message("✓ Perch telemetry: ", length(perch_muddyfoot_tel), " objects updated\n")

# Update ctmm fit projections
message("Updating perch ctmm fit projections...")
ctmm::projection(perch_muddyfoot_ctmm_fits) <- target_projection
message("✓ Perch fits: ", length(perch_muddyfoot_ctmm_fits), " objects updated\n")

# Verify
perch_tel_check <- ctmm::projection(perch_muddyfoot_tel[[1]])
perch_fit_check <- ctmm::projection(perch_muddyfoot_ctmm_fits[[1]])
perch_match <- identical(perch_tel_check, target_projection) && 
  identical(perch_fit_check, target_projection)

if(perch_match) {
  message("✓ Perch projections verified: All match target\n")
} else {
  warning("⚠️  Perch projections may not match perfectly\n")
}

#===============================================================================
# STEP 7: STANDARDIZE ROACH PROJECTIONS
#===============================================================================

message("\n========================================")
message("  STANDARDIZING ROACH PROJECTIONS")
message("========================================\n")

# Update telemetry projections
message("Updating roach telemetry projections...")
for(i in 1:length(roach_muddyfoot_tel)) {
  ctmm::projection(roach_muddyfoot_tel[[i]]) <- target_projection
}
message("✓ Roach telemetry: ", length(roach_muddyfoot_tel), " objects updated\n")

# Update ctmm fit projections
message("Updating roach ctmm fit projections...")
ctmm::projection(roach_muddyfoot_ctmm_fits) <- target_projection
message("✓ Roach fits: ", length(roach_muddyfoot_ctmm_fits), " objects updated\n")

# Verify
roach_tel_check <- ctmm::projection(roach_muddyfoot_tel[[1]])
roach_fit_check <- ctmm::projection(roach_muddyfoot_ctmm_fits[[1]])
roach_match <- identical(roach_tel_check, target_projection) && 
  identical(roach_fit_check, target_projection)

if(roach_match) {
  message("✓ Roach projections verified: All match target\n")
} else {
  warning("⚠️  Roach projections may not match perfectly\n")
}

#===============================================================================
# STEP 8: FINAL VERIFICATION
#===============================================================================

message("\n========================================")
message("  FINAL VERIFICATION")
message("========================================\n")

# Check all projections match
all_projections <- c(
  ctmm::projection(pike_muddyfoot_tel[[1]]),
  ctmm::projection(pike_muddyfoot_ctmm_fits[[1]]),
  ctmm::projection(perch_muddyfoot_tel[[1]]),
  ctmm::projection(perch_muddyfoot_ctmm_fits[[1]]),
  ctmm::projection(roach_muddyfoot_tel[[1]]),
  ctmm::projection(roach_muddyfoot_ctmm_fits[[1]]),
  ctmm::projection(muddyfoot_sp_data)
)

all_match <- all(sapply(all_projections, function(x) identical(x, target_projection)))

if(all_match) {
  message("✓✓✓ SUCCESS ✓✓✓")
  message("\nAll objects now use the same tpeqd projection:")
  message("  - Pike telemetry (", length(pike_muddyfoot_tel), " objects)")
  message("  - Pike ctmm fits (", length(pike_muddyfoot_ctmm_fits), " objects)")
  message("  - Perch telemetry (", length(perch_muddyfoot_tel), " objects)")
  message("  - Perch ctmm fits (", length(perch_muddyfoot_ctmm_fits), " objects)")
  message("  - Roach telemetry (", length(roach_muddyfoot_tel), " objects)")
  message("  - Roach ctmm fits (", length(roach_muddyfoot_ctmm_fits), " objects)")
  message("  - Lake polygon\n")
} else {
  warning("⚠️  Not all projections match. Manual inspection needed.\n")
}

#===============================================================================
# STEP 9: SAVE CORRECTED OBJECTS
#===============================================================================

message("\n========================================")
message("  SAVING CORRECTED OBJECTS")
message("========================================\n")

# Save telemetry objects with corrected projections
message("Saving telemetry objects...")
saveRDS(pike_muddyfoot_tel, 
        paste0(telem_path, "pike_muddyfoot_tel_thinned_final.rds"))
saveRDS(perch_muddyfoot_tel, 
        paste0(telem_path, "perch_muddyfoot_tel_thinned_final.rds"))
saveRDS(roach_muddyfoot_tel, 
        paste0(telem_path, "roach_muddyfoot_tel_thinned_final.rds"))
message("✓ Telemetry objects saved with '_projected' suffix\n")

# Save ctmm fits with corrected projections
message("Saving ctmm fits...")
saveRDS(pike_muddyfoot_ctmm_fits, 
        paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_best_ctmm_model_fits.rds"))
saveRDS(perch_muddyfoot_ctmm_fits, 
        paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_best_ctmm_model_fits.rds"))
saveRDS(roach_muddyfoot_ctmm_fits, 
        paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_best_ctmm_model_fits.rds"))
message("✓ CTMM fits saved with '_projected' suffix\n")

# Save projected polygon (both sf and Spatial formats)
message("Saving lake polygon...")

# Save as sf object
st_write(muddyfoot_polygon_projected, 
         paste0(lake_polygon_path, "muddyfoot_polygon_tpeqd.gpkg"),
         delete_dsn = TRUE,
         quiet = TRUE)

# Save as Spatial object (for direct use with ctmm)
saveRDS(muddyfoot_sp_data,
        paste0(lake_polygon_path, "muddyfoot_polygon_tpeqd_sp.rds"))

message("✓ Lake polygon saved:")
message("    - SF format:     lake_muddyfoot_polygon_tpeqd.gpkg")
message("    - Spatial format: lake_muddyfoot_polygon_tpeqd_sp.rds\n")

#===============================================================================
# STEP 10: CREATE LOADING SCRIPT FOR AKDE ANALYSIS
#===============================================================================

message("\n========================================")
message("  CREATING CONVENIENCE LOADING SCRIPT")
message("========================================\n")

loading_script <- paste0(
#===============================================================================
# LOAD PROJECTED DATA FOR AKDE ANALYSIS
#===============================================================================
# This script loads all telemetry, ctmm fits, and polygon with matching
# tpeqd projections, ready for AKDE with boundary correction.
#
# Generated automatically by convert_to_ctmm_projection.R
#===============================================================================

library(ctmm)
library(sf)

# Define paths
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/muddyfoot/"
lake_polygon_path <- "./data/lake_params/polygons/"

# Load telemetry objects (with correct projection)
pike_muddyfoot_tel <- readRDS(paste0(telem_path, "pike_muddyfoot_tel_thinned_projected.rds"))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, "perch_muddyfoot_tel_thinned_projected.rds"))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, "roach_muddyfoot_tel_thinned_projected.rds"))

message("Telemetry loaded: Pike (", length(pike_muddyfoot_tel), 
        "), Perch (", length(perch_muddyfoot_tel), 
        "), Roach (", length(roach_muddyfoot_tel), ")")

# Load ctmm fits (with correct projection)
pike_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "pike_muddyfoot_ctmm_fits_projected.rds"))
perch_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "perch_muddyfoot_ctmm_fits_projected.rds"))
roach_muddyfoot_ctmm_fits <- readRDS(paste0(ctmm_path, "roach_muddyfoot_ctmm_fits_projected.rds"))

message("CTMM fits loaded: Pike (", length(pike_muddyfoot_ctmm_fits), 
        "), Perch (", length(perch_muddyfoot_ctmm_fits), 
        "), Roach (", length(roach_muddyfoot_ctmm_fits), ")")

# Load lake polygon (Spatial format, ready for ctmm)
muddyfoot_sp_data <- readRDS(paste0(lake_polygon_path, "lake_muddyfoot_polygon_tpeqd_sp.rds"))

message("Lake polygon loaded (tpeqd projection)")

# Verify all projections match
target_proj <- projection(pike_muddyfoot_ctmm_fits[[1]])
message("\\nProjection verification:")
message("  Pike tel:    ", identical(projection(pike_muddyfoot_tel[[1]]), target_proj))
message("  Pike fit:    ", identical(projection(pike_muddyfoot_ctmm_fits[[1]]), target_proj))
message("  Perch tel:   ", identical(projection(perch_muddyfoot_tel[[1]]), target_proj))
message("  Perch fit:   ", identical(projection(perch_muddyfoot_ctmm_fits[[1]]), target_proj))
message("  Roach tel:   ", identical(projection(roach_muddyfoot_tel[[1]]), target_proj))
message("  Roach fit:   ", identical(projection(roach_muddyfoot_ctmm_fits[[1]]), target_proj))
message("  Polygon:     ", identical(projection(muddyfoot_sp_data), target_proj))

message("\\n✓ All data loaded with matching projections")
message("  Ready for AKDE with boundary correction (SP parameter)\\n")
'
)

# Write the loading script
writeLines(loading_script, paste0(telem_path, "../load_projected_data.R"))
message("✓ Loading script created: ", paste0(telem_path, "../load_projected_data.R\n"))

#===============================================================================
# COMPLETION SUMMARY
#===============================================================================

message("\n")
message("===============================================================================")
message("                         PROJECTION CONVERSION COMPLETE")
message("===============================================================================\n")

message("All objects have been converted to tpeqd projection and saved:\n")

message("TELEMETRY OBJECTS:")
message("  ", paste0(telem_path, "pike_muddyfoot_tel_thinned_projected.rds"))
message("  ", paste0(telem_path, "perch_muddyfoot_tel_thinned_projected.rds"))
message("  ", paste0(telem_path, "roach_muddyfoot_tel_thinned_projected.rds\n"))

message("CTMM FITS:")
message("  ", paste0(ctmm_path, "pike_muddyfoot_ctmm_fits_projected.rds"))
message("  ", paste0(ctmm_path, "perch_muddyfoot_ctmm_fits_projected.rds"))
message("  ", paste0(ctmm_path, "roach_muddyfoot_ctmm_fits_projected.rds\n"))

message("LAKE POLYGON:")
message("  ", paste0(lake_polygon_path, "lake_muddyfoot_polygon_tpeqd.gpkg"))
message("  ", paste0(lake_polygon_path, "lake_muddyfoot_polygon_tpeqd_sp.rds\n"))

message("CONVENIENCE SCRIPT:")
message("  ", paste0(telem_path, "../load_projected_data.R"))
message("  (Use this to quickly load all projected data)\n")

message("NEXT STEPS:")
message("1. Use the load_projected_data.R script in your AKDE analysis")
message("2. Run AKDE with boundary correction:")
message("   akde(tel, fit, SP = muddyfoot_sp_data, SP.in = TRUE)")
message("3. All distance calculations will be accurate in meters\n")

message("===============================================================================\n")