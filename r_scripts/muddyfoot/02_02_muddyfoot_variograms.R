#==============================================================================
# Prepare Telemetry Objects and Explore Movement Patterns - Lake Muddyfoot
#==============================================================================
# This script:
# 1. Loads and prepares filtered tracking data
# 2. Converts data to ctmm telemetry objects
# 3. Incorporates UERE (User Equivalent Range Error)
# 4. Explores movement patterns using variograms
# 5. Creates diagnostic plots and HTML reports
# 6. Saves species-specific telemetry objects for modeling
#==============================================================================

### LIBRARIES ###
library(data.table)   # For fast data manipulation
library(tidyverse)    # For data wrangling
library(ctmm)         # For continuous-time movement modeling
library(sf)           # For handling spatial data

# Set the time zone to ensure consistent time handling
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
telem_path <- "./data/telem_obj/muddyfoot/"
figure_path <- "./figures/muddyfoot/variograms/"

#==============================================================================
# 1. LOAD AND PREPARE DATA
#==============================================================================

# Load filtered tracking data -----------------------------------------------
muddyfoot_sub <- readRDS(paste0(filtered_data_path, '03_muddyfoot_sub.rds'))
message("\nLoaded ", nrow(muddyfoot_sub), " detections for ", 
        n_distinct(muddyfoot_sub$individual_ID), " individuals")

# Convert to Movebank format for ctmm package ------------------------------
muddyfoot_movebank <- with(
  muddyfoot_sub,
  data.frame(
    "timestamp" = timestamp,
    "location.long" = Long,
    "location.lat" = Lat,
    "GPS.HDOP" = HDOP,
    "individual-local-identifier" = individual_ID,
    "species" = species,
    "weight" = weight,
    "total_length" = total_length,
    "std_length" = std_length,
    "treatment" = treatment,
    "date" = date,
    "exp_stage" = exp_stage,
    "time_of_day" = time_of_day,
    "found_alive" = found_alive,
    "known_predated" = known_predated
  )
)

# Free up memory ------------------------------------------------------------
rm(muddyfoot_sub)
gc()

# Convert to telemetry object -----------------------------------------------
muddyfoot_tels <- as.telemetry(
  muddyfoot_movebank,
  timezone = "Europe/Stockholm",
  timeformat = "%Y-%m-%d %H:%M:%S",
  projection = NULL,  # Will be set to geometric median automatically
  datum = "WGS84",
  keep = c("species", "weight", "total_length", "std_length", "treatment",
           "date", "exp_stage", "time_of_day", "found_alive", "known_predated")
)

message("Telemetry object created with ", length(muddyfoot_tels), " individuals")
message("Projection: ", projection(muddyfoot_tels[[1]]))
message("Timezone: ", tz(muddyfoot_tels[[1]]$timestamp))

# Incorporate UERE error into telemetry objects ------------------------------
muddyfoot_UERE <- readRDS(paste0(telem_path, "muddyfoot_UERE.rds"))
print(summary(muddyfoot_UERE))
uere(muddyfoot_tels) <- muddyfoot_UERE

#==============================================================================
# 2. ORGANIZE INDIVIDUALS BY SPECIES
#==============================================================================

# Verify individual order ---------------------------------------------------
species_order <- muddyfoot_movebank %>%
  select(species, individual.local.identifier) %>%
  distinct() %>%
  arrange(individual.local.identifier)

message("\nSpecies distribution:")
print(table(species_order$species))

# Display the ordering to help with indexing
message("\nIndividual-Species mapping:")
print(species_order)

# Split telemetry objects by species ---------------------------------------
pike_ids <- species_order %>% 
  filter(species == "Northern Pike") %>% 
  pull(individual.local.identifier)
perch_ids <- species_order %>% 
  filter(species == "Perch") %>% 
  pull(individual.local.identifier)
roach_ids <- species_order %>% 
  filter(species == "Roach") %>% 
  pull(individual.local.identifier)

pike_muddyfoot_tel <- muddyfoot_tels[names(muddyfoot_tels) %in% pike_ids]
perch_muddyfoot_tel <- muddyfoot_tels[names(muddyfoot_tels) %in% perch_ids]
roach_muddyfoot_tel <- muddyfoot_tels[names(muddyfoot_tels) %in% roach_ids]

message("\nSpecies groups created:")
message("  Pike: ", length(pike_muddyfoot_tel), " individuals")
message("  Perch: ", length(perch_muddyfoot_tel), " individuals")
message("  Roach: ", length(roach_muddyfoot_tel), " individuals")

#==============================================================================
# 3. VARIOGRAM ANALYSIS FUNCTIONS
#==============================================================================

# Function to plot variograms for all individuals --------------------------
plot_variograms <- function(tel_list, species_name, output_dir) {
  
  n_individuals <- length(tel_list)
  message("\n", strrep("-", 70))
  message("Processing ", species_name, ": ", n_individuals, " individuals")
  message(strrep("-", 70))
  
  # Calculate grid dimensions for plotting
  n_cols <- ceiling(sqrt(n_individuals))
  n_rows <- ceiling(n_individuals / n_cols)
  
  # Create output filename
  pdf_filename <- paste0(output_dir, species_name, "_variograms_all_individuals.pdf")
  
  # Open PDF device
  pdf(pdf_filename, width = 4 * n_cols, height = 4 * n_rows)
  
  # Set up multi-panel plot
  par(mfrow = c(n_rows, n_cols))
  
  # Loop through each individual
  for (i in 1:n_individuals) {
    ind_name <- names(tel_list)[i]
    tel_data <- tel_list[[i]]
    
    message("  [", i, "/", n_individuals, "] ", ind_name)
    
    # Calculate and plot variogram
    tryCatch({
      vario <- variogram(tel_data)
      plot(vario, main = ind_name, cex.main = 1.2)
      
    }, error = function(e) {
      message("    ERROR: ", e$message)
      # Create empty plot with error message
      plot.new()
      text(0.5, 0.5, paste0("Error:\n", ind_name), cex = 1.2)
    })
  }
  
  dev.off()
  
  message("\nSaved: ", pdf_filename)
  
  return(invisible(NULL))
}

# Function to analyze variogram characteristics ----------------------------
analyze_variograms <- function(tel_list, species_name, output_dir) {
  
  n_individuals <- length(tel_list)
  
  # Create data frame to store variogram characteristics
  vario_summary <- data.frame(
    species = character(),
    individual = character(),
    n_points = integer(),
    timespan_days = numeric(),
    has_asymptote = logical(),
    range_estimate_m = numeric(),
    tau_position = numeric(),
    tau_velocity = numeric(),
    irregularities_detected = character(),
    stringsAsFactors = FALSE
  )
  
  message("\n", strrep("-", 70))
  message("Analyzing variograms for ", species_name, "...")
  message(strrep("-", 70))
  
  for (i in 1:n_individuals) {
    ind_name <- names(tel_list)[i]
    tel_data <- tel_list[[i]]
    
    tryCatch({
      # Calculate variogram
      vario <- variogram(tel_data)
      
      # Extract basic info
      n_points <- nrow(tel_data)
      timespan_days <- as.numeric(diff(range(tel_data$timestamp)), units = "days")
      
      # Check for asymptote (range residency indicator)
      # If SVF plateaus, animal is range resident
      has_asymptote <- !is.null(vario$SVF) && 
        length(vario$SVF) > 10 && 
        sd(tail(vario$SVF, 5)) < mean(tail(vario$SVF, 5)) * 0.1
      
      # Try to estimate range size (if asymptote exists)
      range_estimate <- if(has_asymptote && !is.null(vario$SVF)) {
        max(vario$SVF, na.rm = TRUE)  # Maximum semi-variance as range proxy
      } else {
        NA_real_
      }
      
      # Extract timescales if available
      tau_pos <- if(!is.null(vario$tau) && "position" %in% names(vario$tau)) {
        vario$tau["position"]
      } else {
        NA_real_
      }
      
      tau_vel <- if(!is.null(vario$tau) && "velocity" %in% names(vario$tau)) {
        vario$tau["velocity"]
      } else {
        NA_real_
      }
      
      # Detect irregularities
      irregularities <- c()
      
      # Check for gaps in data
      time_diffs <- diff(as.numeric(tel_data$timestamp))
      median_interval <- median(time_diffs)
      large_gaps <- sum(time_diffs > 5 * median_interval)
      if (large_gaps > 0) {
        irregularities <- c(irregularities, paste0("Large gaps: ", large_gaps))
      }
      
      # Check for very short tracking duration
      if (timespan_days < 7) {
        irregularities <- c(irregularities, "Short duration (<7 days)")
      }
      
      # Check for insufficient data
      if (n_points < 50) {
        irregularities <- c(irregularities, "Few relocations (<50)")
      }
      
      irregularities_text <- if(length(irregularities) > 0) {
        paste(irregularities, collapse = "; ")
      } else {
        "None detected"
      }
      
      # Add to summary
      vario_summary <- rbind(vario_summary, data.frame(
        species = species_name,
        individual = ind_name,
        n_points = n_points,
        timespan_days = round(timespan_days, 1),
        has_asymptote = has_asymptote,
        range_estimate_m = round(range_estimate, 0),
        tau_position = round(tau_pos, 2),
        tau_velocity = round(tau_vel, 2),
        irregularities_detected = irregularities_text,
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      message("  ERROR for ", ind_name, ": ", e$message)
    })
  }
  
  # Save summary
  csv_filename <- paste0(output_dir, species_name, "_variogram_summary.csv")
  write.csv(vario_summary, csv_filename, row.names = FALSE)
  message("\nSaved summary: ", csv_filename)
  
  return(vario_summary)
}

# Function to create interactive HTML report -------------------------------
create_html_report <- function(all_summaries, output_dir) {
  
  # Combine all summaries
  combined <- do.call(rbind, all_summaries)
  
  # Create HTML content
  html_content <- paste0(
    "<!DOCTYPE html>
<html>
<head>
  <title>Variogram Analysis Report</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    h1 { color: #2c3e50; }
    h2 { color: #34495e; margin-top: 30px; }
    table { border-collapse: collapse; width: 100%; margin-top: 20px; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #3498db; color: white; }
    tr:nth-child(even) { background-color: #f2f2f2; }
    .warning { background-color: #fff3cd; }
    .good { background-color: #d4edda; }
    .summary-box { 
      background-color: #ecf0f1; 
      padding: 15px; 
      border-radius: 5px; 
      margin: 20px 0; 
    }
  </style>
</head>
<body>
  <h1>Movement Data Variogram Analysis</h1>
  <p>Generated: ", Sys.time(), "</p>
  
  <div class='summary-box'>
    <h2>Summary Statistics</h2>
    <p><strong>Total Individuals:</strong> ", nrow(combined), "</p>
    <p><strong>Species:</strong> ", paste(unique(combined$species), collapse = ", "), "</p>
    <p><strong>Range Resident Individuals:</strong> ", 
    sum(combined$has_asymptote, na.rm = TRUE), " (", 
    round(100 * sum(combined$has_asymptote, na.rm = TRUE) / nrow(combined), 1), "%)</p>
    <p><strong>Individuals with Irregularities:</strong> ", 
    sum(combined$irregularities_detected != "None detected"), "</p>
  </div>
  
  <h2>Interpretation Guide</h2>
  <ul>
    <li><strong>Has Asymptote:</strong> TRUE indicates range residency (animal has a defined home range)</li>
    <li><strong>Range Estimate:</strong> Approximate size of home range in meters (when asymptote present)</li>
    <li><strong>Tau Position:</strong> Timescale of position autocorrelation</li>
    <li><strong>Tau Velocity:</strong> Timescale of velocity autocorrelation</li>
    <li><strong>Irregularities:</strong> Data quality issues that may affect analysis</li>
  </ul>
  
  <h2>Detailed Results by Species</h2>
  "
  )
  
  # Add tables for each species
  for (sp in unique(combined$species)) {
    sp_data <- combined[combined$species == sp, ]
    
    html_content <- paste0(html_content, "
    <h3>", sp, "</h3>
    <table>
      <tr>
        <th>Individual</th>
        <th>N Points</th>
        <th>Duration (days)</th>
        <th>Range Resident</th>
        <th>Range Est. (m)</th>
        <th>Irregularities</th>
      </tr>
    ")
    
    for (i in 1:nrow(sp_data)) {
      row_class <- if(sp_data$irregularities_detected[i] != "None detected") {
        "warning"
      } else if(sp_data$has_asymptote[i]) {
        "good"
      } else {
        ""
      }
      
      html_content <- paste0(html_content, "
      <tr class='", row_class, "'>
        <td>", sp_data$individual[i], "</td>
        <td>", sp_data$n_points[i], "</td>
        <td>", sp_data$timespan_days[i], "</td>
        <td>", ifelse(sp_data$has_asymptote[i], "YES", "NO"), "</td>
        <td>", ifelse(is.na(sp_data$range_estimate_m[i]), "-", 
                      format(sp_data$range_estimate_m[i], big.mark = ",")), "</td>
        <td>", sp_data$irregularities_detected[i], "</td>
      </tr>
      ")
    }
    
    html_content <- paste0(html_content, "
    </table>
    ")
  }
  
  html_content <- paste0(html_content, "
</body>
</html>
  ")
  
  # Save HTML report
  html_filename <- paste0(output_dir, "variogram_analysis_report.html")
  writeLines(html_content, html_filename)
  message("\nHTML report saved: ", html_filename)
}

#==============================================================================
# 4. RUN VARIOGRAM ANALYSIS
#==============================================================================

#load telemetry objects if necessary
pike_muddyfoot_tel <- readRDS(paste0(telem_path, "pike_muddyfoot_tel_thinned.rds"))
perch_muddyfoot_tel <- readRDS(paste0(telem_path, "perch_muddyfoot_tel_thinned.rds"))
roach_muddyfoot_tel <- readRDS(paste0(telem_path, "roach_muddyfoot_tel_thinned.rds"))


# Plot variograms
if (length(pike_muddyfoot_tel) > 0) {
  plot_variograms(pike_muddyfoot_tel, "Northern_Pike", figure_path)
}
if (length(perch_muddyfoot_tel) > 0) {
  plot_variograms(perch_muddyfoot_tel, "Perch", figure_path)
}
if (length(roach_muddyfoot_tel) > 0) {
  plot_variograms(roach_muddyfoot_tel, "Roach", figure_path)
}

# Analyze variograms
pike_summary <- if (length(pike_muddyfoot_tel) > 0) {
  analyze_variograms(pike_muddyfoot_tel, "Northern_Pike", figure_path)
} else {
  NULL
}

perch_summary <- if (length(perch_muddyfoot_tel) > 0) {
  analyze_variograms(perch_muddyfoot_tel, "Perch", figure_path)
} else {
  NULL
}

roach_summary <- if (length(roach_muddyfoot_tel) > 0) {
  analyze_variograms(roach_muddyfoot_tel, "Roach", figure_path)
} else {
  NULL
}

# Create HTML report
all_summaries <- Filter(Negate(is.null), list(pike_summary, perch_summary, roach_summary))
if (length(all_summaries) > 0) {
  create_html_report(all_summaries, figure_path)
}

#==============================================================================
# 5. ADDITIONAL EXPLORATORY ANALYSES (recommended by ctmm workflow)
#==============================================================================

message("\n", strrep("=", 80))
message("ADDITIONAL EXPLORATORY ANALYSES")
message(strrep("=", 80))

# Function to create summary statistics ------------------------------------
create_tracking_summary <- function(tel_list, species_name) {
  
  message("\nSummary for ", species_name, ":")
  
  summary_df <- data.frame(
    individual = names(tel_list),
    n_locations = sapply(tel_list, nrow),
    duration_days = sapply(tel_list, function(x) {
      as.numeric(diff(range(x$timestamp)), units = "days")
    }),
    median_interval_hours = sapply(tel_list, function(x) {
      median(diff(as.numeric(x$timestamp))) / 3600
    }),
    min_interval_hours = sapply(tel_list, function(x) {
      min(diff(as.numeric(x$timestamp))) / 3600
    }),
    max_interval_hours = sapply(tel_list, function(x) {
      max(diff(as.numeric(x$timestamp))) / 3600
    })
  )
  
  print(summary_df)
  
  # Save summary
  csv_file <- paste0(figure_path, species_name, "_tracking_summary.csv")
  write.csv(summary_df, csv_file, row.names = FALSE)
  message("Saved: ", csv_file)
  
  return(summary_df)
}

# Create summaries for each species
if (length(pike_muddyfoot_tel) > 0) {
  pike_tracking_summary <- create_tracking_summary(pike_muddyfoot_tel, "Northern_Pike")
}
if (length(perch_muddyfoot_tel) > 0) {
  perch_tracking_summary <- create_tracking_summary(perch_muddyfoot_tel, "Perch")
}
if (length(roach_muddyfoot_tel) > 0) {
  roach_tracking_summary <- create_tracking_summary(roach_muddyfoot_tel, "Roach")
}

#==============================================================================
# 6. SAVE TELEMETRY OBJECTS
#==============================================================================

message("\n", strrep("=", 80))
message("SAVING TELEMETRY OBJECTS")
message(strrep("=", 80))

# Save species-specific telemetry objects
if (length(pike_muddyfoot_tel) > 0) {
  saveRDS(pike_muddyfoot_tel, paste0(telem_path, "pike_muddyfoot_tel_thinned.rds"))
  message("Saved: ", telem_path, "pike_muddyfoot_tel_thinned.rds")
}

if (length(perch_muddyfoot_tel) > 0) {
  saveRDS(perch_muddyfoot_tel, paste0(telem_path, "perch_muddyfoot_tel_thinned.rds"))
  message("Saved: ", telem_path, "perch_muddyfoot_tel_thinned.rds")
}

if (length(roach_muddyfoot_tel) > 0) {
  saveRDS(roach_muddyfoot_tel, paste0(telem_path, "roach_muddyfoot_tel_thinned.rds"))
  message("Saved: ", telem_path, "roach_muddyfoot_tel_thinned.rds")
}

#==============================================================================
# 7. FINAL SUMMARY
#==============================================================================

message("\n", strrep("=", 80))
message("EXPLORATION COMPLETE")
message(strrep("=", 80))

message("\nOutputs created:")
message("  Telemetry objects: ", telem_path)
message("  Variogram plots: ", figure_path, "*_variograms_all_individuals.pdf")
message("  Variogram summaries: ", figure_path, "*_variogram_summary.csv")
message("  Tracking summaries: ", figure_path, "*_tracking_summary.csv")
message("  HTML report: ", figure_path, "variogram_analysis_report.html")

message("\nNext steps:")
message("  1. Review the variogram plots and HTML report")
message("  2. Identify any problematic individuals or data quality issues")
message("  3. Run script 02_fit_ctmm_models.R to fit movement models")

message("\n", strrep("=", 80))
