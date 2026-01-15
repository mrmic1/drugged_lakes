#==============================================================================
# Run CTMM models for each individual and species - Lake Muddyfoot
#==============================================================================

### LIBRARIES ###
library(data.table)   # For fast data manipulation
library(tidyverse)    # For data wrangling
library(ctmm)         # For continuous-time movement modeling
library(sf)           # For handling spatial data
library(parallel)     # For parallel processing
library(foreach)      # For parallel for loops
library(doParallel)   # For registering parallel backend
library(gridExtra)

# Set the time zone to ensure consistent time handling
Sys.setenv(TZ = 'Europe/Stockholm')

# Define file paths for reading and saving filtered telemetry and ctmm model results
filtered_data_path <- "./data/tracks_filtered/muddyfoot/"
ctmm_path <- "./data/ctmm_fits/"
telem_path <- "./data/telem_obj/muddyfoot/"
figure_path <- "./figures"

#==============================================================================
# 1. LOAD AND PREPARE DATA
#==============================================================================

# Load filtered tracking data -----------------------------------------------
muddyfoot_sub <- readRDS(paste0(filtered_data_path, '03_muddyfoot_sub.rds'))
message("Loaded ", nrow(muddyfoot_sub), " detections for ", 
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
message("\nConverting to ctmm telemetry object...")
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
#load UERE
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

print(table(species_order$species))

# Display the ordering to help with indexing
print(species_order)

# Split telemetry objects by species ---------------------------------------
# Note: You'll need to adjust these indices based on your actual data
# Check the species_order output above to determine correct indices

# For now, let's create a flexible approach:
pike_ids <- species_order %>% filter(species == "Northern Pike") %>% pull(individual.local.identifier)
perch_ids <- species_order %>% filter(species == "Perch") %>% pull(individual.local.identifier)
roach_ids <- species_order %>% filter(species == "Roach") %>% pull(individual.local.identifier)

pike_muddyfoot_tel <- muddyfoot_tels[names(muddyfoot_tels) %in% pike_ids]
perch_muddyfoot_tel <- muddyfoot_tels[names(muddyfoot_tels) %in% perch_ids]
roach_muddyfoot_tel <- muddyfoot_tels[names(muddyfoot_tels) %in% roach_ids]

message("\nSpecies groups created:")
message("Pike: ", length(pike_muddyfoot_tel), " individuals")
message("Perch: ", length(perch_muddyfoot_tel), " individuals")
message("Roach: ", length(roach_muddyfoot_tel), " individuals")

# Save species-specific telemetry objects -----------------------------------
saveRDS(pike_muddyfoot_tel, paste0(telem_path, "pike_muddyfoot_tel_thinned.rds"))
saveRDS(perch_muddyfoot_tel, paste0(telem_path, "perch_muddyfoot_tel_thinned.rds"))
saveRDS(roach_muddyfoot_tel, paste0(telem_path, "roach_muddyfoot_tel_thinned.rds"))

plot_variograms <- function(tel_list, species_name, output_dir) {
  
  n_individuals <- length(tel_list)
  message(paste0("\nProcessing ", species_name, ": ", n_individuals, " individuals"))
  
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
    
    message(paste0("  Processing individual ", i, "/", n_individuals, ": ", ind_name))
    
    # Calculate and plot variogram
    tryCatch({
      vario <- variogram(tel_data)
      
      # Plot directly (ctmm's plot function handles the plotting)
      plot(vario, main = ind_name, cex.main = 1.2)
      
      # Check for range residency indicators
      # Asymptote suggests range residency
      # No asymptote suggests migration/nomadism
      
    }, error = function(e) {
      message(paste0("    ERROR for ", ind_name, ": ", e$message))
      # Create empty plot with error message
      plot.new()
      text(0.5, 0.5, paste0("Error:\n", ind_name), cex = 1.2)
    })
  }
  
  dev.off()
  
  message(paste0("  Saved: ", pdf_filename))
  
  return(invisible(NULL))
}

# Function to create detailed variogram analysis ---------------------------
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
  
  message(paste0("\nAnalyzing variograms for ", species_name, "..."))
  
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
      message(paste0("  ERROR for ", ind_name, ": ", e$message))
    })
  }
  
  # Save summary
  csv_filename <- paste0(output_dir, species_name, "_variogram_summary.csv")
  write.csv(vario_summary, csv_filename, row.names = FALSE)
  message(paste0("  Saved summary: ", csv_filename))
  
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
  message(paste0("\nHTML report saved: ", html_filename))
}

# Run analysis for all species ---------------------------------------------
message(paste(rep("=", 70), collapse = ""))
message("VARIOGRAM ANALYSIS")
message(paste(rep("=", 70), collapse = ""))

# Plot variograms
plot_variograms(pike_muddyfoot_tel, "Northern_Pike", figure_path)
plot_variograms(perch_muddyfoot_tel, "Perch", output_path)
plot_variograms(roach_muddyfoot_tel, "Roach", output_path)

# Analyze variograms
pike_summary <- analyze_variograms(pike_muddyfoot_tel, "Northern_Pike", output_path)
perch_summary <- analyze_variograms(perch_muddyfoot_tel, "Perch", output_path)
roach_summary <- analyze_variograms(roach_muddyfoot_tel, "Roach", output_path)

# Create HTML report
create_html_report(
  list(pike_summary, perch_summary, roach_summary), 
  output_path
)

message("\n", paste(rep("=", 70), collapse = ""))
message("ANALYSIS COMPLETE")
message(paste(rep("=", 70), collapse = ""))
message("\nOutputs saved to: ", output_path)
message("\nFiles created:")
message("  1. PDF plots: [Species]_variograms_all_individuals.pdf")
message("  2. CSV summaries: [Species]_variogram_summary.csv")
message("  3. HTML report: variogram_analysis_report.html")
message("\nInterpretation:")
message("  - Variograms with asymptotes indicate range residency")
message("  - Variograms that continue rising suggest migration/nomadism")
message("  - Irregular patterns may indicate data quality issues")
message("  - Check the HTML report for a comprehensive overview")




















#==============================================================================
# 3. HELPER FUNCTIONS
#==============================================================================

# Function to safely determine number of cores to use ----------------------
get_safe_cores <- function(max_cores = NULL, reserve_cores = 2) {
  available <- detectCores()
  
  if (is.null(max_cores)) {
    # Use 60-70% of available cores, leaving some for system
    safe_cores <- max(1, floor(available * 0.65))
  } else {
    safe_cores <- min(max_cores, available - reserve_cores)
  }
  
  safe_cores <- max(1, safe_cores)  # Ensure at least 1 core
  
  message(sprintf("Using %d cores out of %d available", safe_cores, available))
  return(safe_cores)
}

# Function to assess initial model guess ------------------------------------
assess_ctmm_guess <- function(telem_list) {
  
  message("\n=== Assessing Model Selection for Each Individual ===")
  
  for (i in seq_along(telem_list)) {
    tel_i <- telem_list[[i]]
    id_i <- names(telem_list)[i]
    
    # Run ctmm.guess
    guess_i <- ctmm.guess(
      tel_i,
      CTMM = ctmm(error = TRUE),
      interactive = FALSE
    )
    
    # Determine model type by tau length
    tau_len <- length(guess_i$tau)
    
    model_type <- case_when(
      tau_len == 0 ~ "IID (uncorrelated)",
      tau_len == 1 ~ "OU (range resident, no velocity autocorrelation)",
      tau_len == 2 ~ "OUF (directional persistence + range residency)",
      TRUE ~ "Unknown"
    )
    
    # Print summary
    message(sprintf("ID: %-10s | tau length = %d | guessed model = %s",
                    id_i, tau_len, model_type))
  }
}

# Parallel fitting function -------------------------------------------------
fit_ctmm_species_parallel <- function(telem_list, species_name, lake_name = "muddyfoot",
                                      max_cores = NULL, 
                                      save_individual_fits = TRUE) {
  
  message("\n", strrep("=", 80))
  message("=== FITTING CTMM MODELS FOR ", toupper(species_name), " ===")
  message(strrep("=", 80))
  
  # Assess expected models
  assess_ctmm_guess(telem_list)
  
  # Create output directory
  output_dir <- file.path(save_ctmm_path, paste0(lake_name, "_", species_name, "_fits"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine number of cores
  n_cores <- get_safe_cores(max_cores = max_cores)
  
  # Set up parallel cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to cluster
  clusterExport(cl, c("telem_list", "species_name", "lake_name", "output_dir", 
                      "save_individual_fits", "save_ctmm_path"),
                envir = environment())
  
  message("\n=== Starting Parallel Model Fitting ===")
  message("Processing ", length(telem_list), " individuals using ", n_cores, " cores")
  start_time <- Sys.time()
  
  # Parallel fitting with error handling
  results <- foreach(
    i = 1:length(telem_list),
    .packages = c('ctmm'),
    .errorhandling = 'pass',  # Continue even if one fails
    .verbose = FALSE
  ) %dopar% {
    
    tryCatch({
      tel_i <- telem_list[[i]]
      id_i <- names(telem_list)[i]
      
      ind_start <- Sys.time()
      
      # Get initial guess
      guess_model <- ctmm.guess(
        tel_i,
        CTMM = ctmm(error = TRUE),
        interactive = FALSE
      )
      
      # Fit the model
      model_fit <- ctmm.fit(
        data = tel_i,
        CTMM = guess_model,
        method = "ML"
      )
      
      ind_elapsed <- as.numeric(difftime(Sys.time(), ind_start, units = "secs"))
      
      # Save individual fit if requested
      if (save_individual_fits) {
        output_file <- file.path(output_dir, paste0(id_i, "_ctmm_fit.rds"))
        saveRDS(model_fit, file = output_file)
      }
      
      # Clean up memory
      gc(verbose = FALSE)
      
      list(
        id = id_i,
        fit = model_fit,
        time = ind_elapsed,
        success = TRUE,
        error = NULL
      )
      
    }, error = function(e) {
      list(
        id = names(telem_list)[i],
        fit = NULL,
        time = NA,
        success = FALSE,
        error = as.character(e)
      )
    })
  }
  
  # Stop cluster
  stopCluster(cl)
  
  total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  # Process results
  successful_fits <- sapply(results, function(x) x$success)
  n_success <- sum(successful_fits)
  n_failed <- sum(!successful_fits)
  
  message("\n=== Fitting Complete ===")
  message(sprintf("Total time: %.1f minutes", total_elapsed))
  message(sprintf("Successful fits: %d/%d", n_success, length(telem_list)))
  
  if (n_failed > 0) {
    message(sprintf("\nFailed fits: %d", n_failed))
    failed_ids <- sapply(results[!successful_fits], function(x) x$id)
    message("Failed IDs: ", paste(failed_ids, collapse = ", "))
    
    # Print error messages
    for (i in which(!successful_fits)) {
      message(sprintf("\n%s error: %s", results[[i]]$id, results[[i]]$error))
    }
  }
  
  # Extract successful fits
  ctmm_fits <- lapply(results[successful_fits], function(x) x$fit)
  names(ctmm_fits) <- sapply(results[successful_fits], function(x) x$id)
  
  # Calculate timing statistics
  fit_times <- sapply(results[successful_fits], function(x) x$time)
  message(sprintf("\nTiming stats (seconds):"))
  message(sprintf("  Mean: %.1f  |  Median: %.1f  |  Max: %.1f", 
                  mean(fit_times), median(fit_times), max(fit_times)))
  
  # Save combined fit list
  if (length(ctmm_fits) > 0) {
    output_list_file <- file.path(output_dir, paste0("lake_", lake_name, "_", species_name, "_ctmm_fits.rds"))
    saveRDS(ctmm_fits, output_list_file)
    message("\nCombined fit list saved to: ", output_list_file)
  }
  
  # Return both fits and diagnostic info
  return(list(
    fits = ctmm_fits,
    diagnostics = list(
      total_time = total_elapsed,
      n_success = n_success,
      n_failed = n_failed,
      failed_ids = if(n_failed > 0) sapply(results[!successful_fits], function(x) x$id) else NULL,
      fit_times = fit_times
    )
  ))
}

# Sequential fitting function (backup if parallel causes issues) -----------
fit_ctmm_species_sequential <- function(telem_list, species_name, lake_name = "muddyfoot") {
  
  message("\n", strrep("=", 80))
  message("=== FITTING CTMM MODELS FOR ", toupper(species_name), " (SEQUENTIAL) ===")
  message(strrep("=", 80))
  
  assess_ctmm_guess(telem_list)
  
  ctmm_fits <- vector("list", length(telem_list))
  names(ctmm_fits) <- names(telem_list)
  
  output_dir <- file.path(save_ctmm_path, paste0(lake_name, "_", species_name, "_fits"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("\n=== Starting Sequential Model Fitting ===")
  start_time <- Sys.time()
  fit_times <- numeric(length(telem_list))
  
  for (i in seq_along(telem_list)) {
    tel_i <- telem_list[[i]]
    id_i <- names(telem_list)[i]
    
    # Estimate time remaining
    if (i > 1) {
      avg_time <- mean(fit_times[1:(i-1)])
      remaining <- (length(telem_list) - i + 1) * avg_time / 60
      message(sprintf("\n[%d/%d] %s | Est. remaining: %.1f min", 
                      i, length(telem_list), id_i, remaining))
    } else {
      message(sprintf("\n[%d/%d] %s", i, length(telem_list), id_i))
    }
    
    ind_start <- Sys.time()
    
    guess_model <- ctmm.guess(tel_i, CTMM = ctmm(error = TRUE), interactive = FALSE)
    model_fit <- ctmm.fit(data = tel_i, CTMM = guess_model, method = "ML")
    
    fit_times[i] <- as.numeric(difftime(Sys.time(), ind_start, units = "secs"))
    message(sprintf("  Fitted in %.1f seconds", fit_times[i]))
    
    output_file <- file.path(output_dir, paste0(id_i, "_ctmm_fit.rds"))
    saveRDS(model_fit, file = output_file)
    
    ctmm_fits[[i]] <- model_fit
    gc()
  }
  
  total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf("\n=== Completed in %.1f minutes ===", total_elapsed))
  
  output_list_file <- file.path(output_dir, paste0("lake_", lake_name, "_", species_name, "_ctmm_fits.rds"))
  saveRDS(ctmm_fits, output_list_file)
  
  return(ctmm_fits)
}

# Function to verify model fits ---------------------------------------------
verify_fits <- function(telem_list, fit_list, species_name, n_check = NULL) {
  
  message("\n=== Verifying ", species_name, " Model Fits ===")
  
  if (length(fit_list) == 0) {
    message("No fits to verify!")
    return(invisible())
  }
  
  # If n_check is NULL, check all individuals
  if (is.null(n_check)) {
    n_check <- length(fit_list)
  }
  
  # Ensure n_check doesn't exceed available fits
  n_check <- min(n_check, length(fit_list))
  
  message(sprintf("\nChecking %d out of %d individuals", n_check, length(fit_list)))
  
  # Loop through individuals to check
  for (i in 1:n_check) {
    id_i <- names(fit_list)[i]
    
    message("\n", strrep("-", 80))
    message(sprintf("Individual %d/%d: %s", i, n_check, id_i))
    message(strrep("-", 80))
    
    # Print summary
    print(summary(fit_list[[i]]))
    
    # Create diagnostic plot
    message("\nGenerating diagnostic plot...")
    plot(telem_list[[id_i]], fit_list[[i]], error = FALSE)
    
    # Pause for user to review (except for last one)
    if (i < n_check) {
      message("\nPress [Enter] to continue to next individual, or type 'q' to quit verification...")
      user_input <- readline()
      if (tolower(trimws(user_input)) == 'q') {
        message("Verification stopped by user.")
        break
      }
    }
  }
  
  message("\n", strrep("=", 80))
  message("Verification complete!")
  message(strrep("=", 80))
}

# Quick verification function (summaries only, no plots) -------------------
verify_fits_quick <- function(fit_list, species_name) {
  
  message("\n=== Quick Summary of ", species_name, " Model Fits ===")
  
  if (length(fit_list) == 0) {
    message("No fits to verify!")
    return(invisible())
  }
  
  # Create a summary dataframe
  summary_df <- data.frame(
    ID = names(fit_list),
    Model = character(length(fit_list)),
    AIC = numeric(length(fit_list)),
    DOF_area = numeric(length(fit_list)),
    DOF_speed = numeric(length(fit_list)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(fit_list)) {
    model_sum <- summary(fit_list[[i]])
    summary_df$Model[i] <- names(fit_list[[i]]$tau)[1]
    summary_df$AIC[i] <- model_sum$AIC
    
    # Extract DOF information if available
    if ("DOF" %in% names(model_sum)) {
      summary_df$DOF_area[i] <- model_sum$DOF["area"]
      summary_df$DOF_speed[i] <- ifelse("speed" %in% names(model_sum$DOF), 
                                        model_sum$DOF["speed"], NA)
    }
  }
  
  print(summary_df)
  
  message("\nModel type distribution:")
  print(table(summary_df$Model))
  
  return(invisible(summary_df))
}

#==============================================================================#
#### 4. FIT CTMM MODELS - PIKE ####
#==============================================================================#

# Option to reload if needed ------------------------------------------------
pike_muddyfoot_tel <- readRDS(paste0(save_telem_path, "muddyfoot/pike_muddyfoot_tel_thinned.rds"))

if (length(pike_muddyfoot_tel) > 0) {
  # Fit models using parallel processing -----------------------------------
  # Use max_cores=3 for small groups to avoid overloading PC
  pike_results <- fit_ctmm_species_parallel(
    pike_muddyfoot_tel, 
    "pike",
    lake_name = "muddyfoot",
    max_cores = 3  # Conservative setting to reduce noise/heat
  )
  
  muddyfoot_pike_ctmm_fits <- pike_results$fits
  
# Verify fits ----------------------------------------------------------- #
  #if you need to reload
  muddyfoot_pike_ctmm_fits <-  readRDS(paste0(save_ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_ctmm_fits.rds"))
  
verify_fits(pike_muddyfoot_tel, muddyfoot_pike_ctmm_fits, "Pike")
} else {
  message("\n*** No pike individuals found in Muddyfoot dataset ***")
  muddyfoot_pike_ctmm_fits <- list()
}

#==============================================================================
# 5. FIT CTMM MODELS - PERCH ####
#==============================================================================

# Option to reload if needed ------------------------------------------------
perch_muddyfoot_tel <- readRDS(paste0(save_telem_path, "muddyfoot/perch_muddyfoot_tel_thinned.rds"))

if (length(perch_muddyfoot_tel) > 0) {
  # Fit models using parallel processing -----------------------------------
  # Adjust max_cores based on number of individuals
  perch_cores <- ifelse(length(perch_muddyfoot_tel) > 15, 6, 3)
  
  perch_results <- fit_ctmm_species_parallel(
    perch_muddyfoot_tel,
    "perch",
    lake_name = "muddyfoot",
    max_cores = perch_cores
  )
  
  muddyfoot_perch_ctmm_fits <- perch_results$fits

  
  
  
# Verify fits ----------------------------------------------------------- #
#if you need to reload
muddyfoot_perch_ctmm_fits <-  readRDS(paste0(save_ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_ctmm_fits.rds"))
verify_fits(perch_muddyfoot_tel, muddyfoot_perch_ctmm_fits, "Perch")
} else {
  message("\n*** No perch individuals found in Muddyfoot dataset ***")
  muddyfoot_perch_ctmm_fits <- list()
}

#==============================================================================
# 6. FIT CTMM MODELS - ROACH
#==============================================================================

# Option to reload if needed ------------------------------------------------
roach_muddyfoot_tel <- readRDS(paste0(save_telem_path, "muddyfoot/roach_muddyfoot_tel_thinned.rds"))

if (length(roach_muddyfoot_tel) > 0) {
  # Fit models using parallel processing -----------------------------------
  # Adjust max_cores based on number of individuals
  roach_cores <- ifelse(length(roach_muddyfoot_tel) > 15, 6, 3)
  
  roach_results <- fit_ctmm_species_parallel(
    roach_muddyfoot_tel,
    "roach",
    lake_name = "muddyfoot",
    max_cores = roach_cores
  )
  
  muddyfoot_roach_ctmm_fits <- roach_results$fits
  
  # Verify fits -----------------------------------------------------------
  verify_fits(roach_muddyfoot_tel, muddyfoot_roach_ctmm_fits, "Roach")
} else {
  message("\n*** No roach individuals found in Muddyfoot dataset ***")
  muddyfoot_roach_ctmm_fits <- list()
}

#==============================================================================
# 7. FINAL SUMMARY
#==============================================================================

message("\n", strrep("=", 80))
message("=== ALL CTMM MODELS COMPLETED - LAKE MUDDYFOOT ===")
message(strrep("=", 80))

message("\nModels fitted:")
message("  Pike: ", length(muddyfoot_pike_ctmm_fits), " individuals")
message("  Perch: ", length(muddyfoot_perch_ctmm_fits), " individuals")
message("  Roach: ", length(muddyfoot_roach_ctmm_fits), " individuals")
message("  Total: ", length(muddyfoot_pike_ctmm_fits) + 
          length(muddyfoot_perch_ctmm_fits) + 
          length(muddyfoot_roach_ctmm_fits), " individuals")

if (exists("pike_results") && length(pike_results$fits) > 0) {
  message("\nTotal processing times:")
  message("  Pike: ", sprintf("%.1f", pike_results$diagnostics$total_time), " minutes")
}
if (exists("perch_results") && length(perch_results$fits) > 0) {
  message("  Perch: ", sprintf("%.1f", perch_results$diagnostics$total_time), " minutes")
}
if (exists("roach_results") && length(roach_results$fits) > 0) {
  message("  Roach: ", sprintf("%.1f", roach_results$diagnostics$total_time), " minutes")
}

message("\nOutput locations:")
message("  Individual fits: ", save_ctmm_path, "muddyfoot_[species]_fits/")
message("  Combined lists: ", save_ctmm_path, "muddyfoot_[species]_fits/muddyfoot_[species]_ctmm_fits.rds")
message("  Telemetry objects: ", save_telem_path, "muddyfoot/")

message("\nTo reload fitted models:")
message("  pike_fits <- readRDS('", save_ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_ctmm_fits.rds')")
message("  perch_fits <- readRDS('", save_ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_ctmm_fits.rds')")
message("  roach_fits <- readRDS('", save_ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_ctmm_fits.rds')")

message("\n", strrep("=", 80))

#==============================================================================
# NOTES ON USAGE
#==============================================================================

# If your PC gets too loud or hot during processing:
# 1. Reduce max_cores to 2-3 for all species
# 2. Or use the sequential version:
#    muddyfoot_pike_ctmm_fits <- fit_ctmm_species_sequential(pike_muddyfoot_tel, "pike", "muddyfoot")
#
# To adjust core usage dynamically:
# - max_cores = NULL will use ~65% of available cores
# - max_cores = 3 is conservative and quiet
# - max_cores = 6 is good for larger groups on modern PCs
#
# The parallel version will save significant time on larger groups!