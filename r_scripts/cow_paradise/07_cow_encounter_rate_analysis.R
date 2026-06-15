# =================================================================-
# Encounter Rate Analysis - Cow Paradise ####
# =================================================================-
#
# ANALYSIS OVERVIEW:
# Analyses encounter rates between prey fish (roach and perch) and
# pike, using encounter data produced by script 05. Addresses three
# questions:
#
# Q1: Do treatment groups (Control vs. Exposed) differ in their
#     overall encounter rate with pike?
#
# Q2: When during the day do prey most regularly encounter predators?
# Q3: Does the diel pattern of encounters differ between treatments?
#
# Analyses conducted on the post-predation-filtered dataset. Encounter
# rates in the pre-filtering dataset are inflated by prey fish that
# were predated.
#
# NOTE — COW PARADISE GPS ERROR:
# GPS error is higher here (σ = 0.827m) than in BT or Muddyfoot, but
# encounter thresholds are set at ≤2.0m (high-confidence) to account for
# the higher GPS error at this lake:
#   High-confidence: ≤2.0m | Probable: 2.0–4.0m
#
#=================================================================-
# =================================================================-
# 1. SETUP ####
# =================================================================-

## 1.1 Load Required Libraries ----
suppressPackageStartupMessages({
  library(dplyr)           # Data manipulation
  library(lubridate)       # Date/time handling
  library(data.table)      # Fast I/O
  library(tidyr)           # Data reshaping (expand/nesting for predictions)
  library(purrr)           # map_df for per-individual ACF
  library(ggplot2)         # Plotting
  library(ggthemes)        # Plot themes (theme_few)
  library(ggpubr)          # Multi-panel figure assembly (ggarrange)
  library(mgcv)            # GAMMs (bam)
  library(gratia)          # GAMM diagnostics (appraise, draw)
  library(itsadug)         # GAMM residual checks (check_resid)
  library(marginaleffects) # Treatment contrasts on response scale
  library(glmmTMB)         # Negative binomial GLMM
  library(DHARMa)          # Residual diagnostics for glmmTMB
  library(emmeans)         # Marginal means and pairwise contrasts
  library(suncalc)         # Sunrise/sunset times for diel annotation
})

Sys.setenv(TZ = "Europe/Stockholm")

## 1.2 Define Directory Paths ----
paths <- list(
  encounters = "./data/encounters/cow_paradise/",
  filtered   = "./data/tracks_filtered/cow_paradise/",
  mortality  = "./data/encounters/suspected_mortality_updated.xlsx",
  figures    = "./figures/cow_paradise/encounter_analysis/",
  tables     = "./tables/cow_paradise/encounter_analysis/"
)

for (p in paths) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

## 1.3 Diel Reference Times ----
# Calculated for mid-experiment at lake coordinates (lat 63.77, lon 20.06)
# using suncalc::getSunlightTimes(). BT, Muddyfoot, and Cow Paradise are
# all in the same region of northern Sweden, so shared diel times apply.
nauticalDawn <- 5.494649   # ~05:30
sunrise      <- 7.189919   # ~07:11
sunset       <- 17.75258   # ~17:45
nauticalDusk <- 19.44784   # ~19:27

## 1.4 Model output helper functions ----
print_glmm_results <- function(model, lake, question, species, response_var) {
  header <- sprintf("=== %s | %s | %s NB-GLMM (%s) ===",
                    lake, question, species, response_var)
  sep  <- paste(rep("=", nchar(header)), collapse = "")
  dash <- paste(rep("-", nchar(header)), collapse = "")
  message("\n", sep); message(header); message(sep)
  message("  Model     : glmmTMB | nbinom2 (log link) | offset = log(n_days)")
  message(sprintf("  Response  : %s", response_var))
  message(sprintf("  AIC       : %.2f", AIC(model)))
  message(sprintf("  Log-Lik   : %.2f", as.numeric(logLik(model))))
  message(sprintf("  Dispersion: theta = %.3f", sigma(model)))
  message("\n  Fixed Effects [log scale | IRR = exp(Estimate)]:")
  s     <- summary(model)
  coefs <- as.data.frame(s$coefficients$cond)
  out   <- data.frame(
    Estimate = round(coefs[["Estimate"]],   4),
    SE       = round(coefs[["Std. Error"]], 4),
    z        = round(coefs[["z value"]],    3),
    p        = round(coefs[["Pr(>|z|)"]],   4),
    IRR      = round(exp(coefs[["Estimate"]]), 3)
  )
  rownames(out) <- rownames(coefs)
  print(out)
  message(dash)
}

print_gamm_results <- function(model, lake, species) {
  header <- sprintf("=== %s | Q2/Q3 | %s GAMM (diel encounter pattern) ===",
                    lake, species)
  sep  <- paste(rep("=", nchar(header)), collapse = "")
  dash <- paste(rep("-", nchar(header)), collapse = "")
  message("\n", sep); message(header); message(sep)
  s <- summary(model)
  message(sprintf("  Family    : %s (log link)", s$family[1]))
  message(sprintf("  R-sq(adj) : %.3f",          s$r.sq))
  message(sprintf("  Dev. expl : %.1f%%",         s$dev.expl * 100))
  message(sprintf("  AIC       : %.2f",           AIC(model)))
  message(sprintf("  n (obs)   : %d",             s$n))
  message("\n  Parametric Coefficients:")
  print(round(s$p.table, 4))
  message("\n  Smooth Terms (approximate significance):")
  print(round(s$s.table, 4))
  message(dash)
}


# =================================================================-
# 2. DATA LOADING ####
# =================================================================-

## 2.1 Load encounter distance data (produced by script 05) ----
# encounter_high_conf = ≤2.0m | encounter_probable = 2.0–4.0m
# (Cow Paradise threshold accounts for higher GPS error: σ = 0.827m)
roach_pike_distances_df <- readRDS(
  paste0(paths$encounters, "cow_pike_roach_distances_tiered_df.rds")
)
perch_pike_distances_df <- readRDS(
  paste0(paths$encounters, "cow_pike_perch_distances_tiered_df.rds")
)

## 2.2 Load metadata ----
cow_filt_data <- readRDS(
  paste0(paths$filtered, "05_cow_sub.rds")
)

# Extracting ID, treatment and weight information
cow_meta <- cow_filt_data %>%
  select(individual_ID, treatment, weight) %>%
  distinct()

# =================================================================-
# 3. FILTER OUT POST-PREDATION ENCOUNTER RECORDS ####
# =================================================================-

## 3.1 Load mortality reference ----
mortality_preds <- readxl::read_excel(paths$mortality)

cow_death_dates <- mortality_preds %>%
  filter(lake == "cow paradise") %>%
  filter(species %in% c("Roach", "Perch")) %>%
  dplyr::select(individual_ID, species, revised_suspected_mortality,
                revised_likely_death_date) %>%
  rename(death_date = revised_likely_death_date) %>%
  mutate(death_date = as.Date(death_date, origin = "1970-01-01"))

print(table(cow_death_dates$species,
            cow_death_dates$revised_suspected_mortality))


## 3.2 Filter roach-pike distance data ----
roach_pike_distances_df <- roach_pike_distances_df %>%
  left_join(
    cow_death_dates %>% select(individual_ID, death_date),
    by = c("Prey_ID" = "individual_ID")
  ) %>%
  filter(is.na(death_date) | Date < death_date) %>%
  select(-death_date)

#before: 6819358
#after: 6704862

## 3.3 Filter perch-pike distance data ----
perch_pike_distances_df <- perch_pike_distances_df %>%
  left_join(
    cow_death_dates %>% select(individual_ID, death_date),
    by = c("Prey_ID" = "individual_ID")
  ) %>%
  filter(is.na(death_date) | Date < death_date) %>%
  select(-death_date)

#before: 8402633
#after: 8043628


# =================================================================-
# 4. BUILD ANALYSIS DATASETS ####
# =================================================================-

## 4.1 Individual-total encounter dataset (for GLMM, Q1) ----
# One row per individual: total high-confidence encounter count and the
# number of days tracked (used as an offset to model a daily rate).
# weight_z is scaled within each species separately.

# Calculate total encounters per individual per category
enc_roach_total <- roach_pike_distances_df %>%
  rename(individual_ID = Prey_ID) %>%
  group_by(individual_ID) %>%
  summarise(
    total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    total_encounters_probable  = sum(encounter_probable,  na.rm = TRUE),
    total_encounters_combined  = sum(encounter_high_conf, na.rm = TRUE) +
      sum(encounter_probable,  na.rm = TRUE),
    n_days                     = n_distinct(Date),
    .groups = "drop"
  ) %>%
  left_join(cow_meta, by = "individual_ID") %>%
  filter(!is.na(treatment)) %>%
  mutate(
    Species       = "Roach",
    treatment     = factor(treatment, levels = c("Control", "Mix")),
    individual_ID = factor(individual_ID),
    weight_z      = as.numeric(scale(weight))
  )

enc_perch_total <- perch_pike_distances_df %>%
  rename(individual_ID = Prey_ID) %>%
  group_by(individual_ID) %>%
  summarise(
    total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    total_encounters_probable  = sum(encounter_probable,  na.rm = TRUE),
    total_encounters_combined  = sum(encounter_high_conf, na.rm = TRUE) +
      sum(encounter_probable,  na.rm = TRUE),
    n_days                     = n_distinct(Date),
    .groups = "drop"
  ) %>%
  left_join(cow_meta, by = "individual_ID") %>%
  filter(!is.na(treatment)) %>%
  mutate(
    Species       = "Perch",
    treatment     = factor(treatment, levels = c("Control", "Mix")),
    individual_ID = factor(individual_ID),
    weight_z      = as.numeric(scale(weight))
  )

## 4.1a Summary statistics — encounter rates by species and treatment ----
enc_rate_summary <- bind_rows(enc_roach_total, enc_perch_total) %>%
  mutate(
    rate_high_conf = total_encounters_high_conf / n_days,
    rate_combined  = total_encounters_combined  / n_days
  ) %>%
  group_by(Species, treatment) %>%
  summarise(
    n                = n(),
    mean_rate_hc     = mean(rate_high_conf,            na.rm = TRUE),
    sd_rate_hc       = sd(rate_high_conf,              na.rm = TRUE),
    median_rate_hc   = median(rate_high_conf,          na.rm = TRUE),
    iqr_rate_hc      = IQR(rate_high_conf,             na.rm = TRUE),
    mean_rate_comb   = mean(rate_combined,             na.rm = TRUE),
    sd_rate_comb     = sd(rate_combined,               na.rm = TRUE),
    median_rate_comb = median(rate_combined,           na.rm = TRUE),
    iqr_rate_comb    = IQR(rate_combined,              na.rm = TRUE),
    total_enc_hc     = sum(total_encounters_high_conf, na.rm = TRUE),
    total_enc_comb   = sum(total_encounters_combined,  na.rm = TRUE),
    .groups = "drop"
  )

message("\nEncounter rate summary by species and treatment:")
print(enc_rate_summary)

message(sprintf(
  "Encounter totals: %d roach, %d perch individuals",
  nrow(enc_roach_total), nrow(enc_perch_total)
))

## 4.2 Hourly encounter dataset (for GAMM, Q2 & Q3) ----
roach_hourly <- roach_pike_distances_df %>%
  mutate(
    timestamp = as.POSIXct(timestamp, tz = "Europe/Stockholm"),
    Date      = as.Date(timestamp, tz = "Europe/Stockholm"),
    hour      = as.integer(format(timestamp, "%H"))
  ) %>%
  group_by(Prey_ID, Date, hour) %>%
  summarise(
    n_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    n_encounters_probable  = sum(encounter_probable,  na.rm = TRUE),
    n_encounters_combined  = sum(encounter_high_conf, na.rm = TRUE) +
      sum(encounter_probable,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(individual_ID = Prey_ID) %>%
  left_join(cow_meta, by = "individual_ID") %>%
  filter(!is.na(treatment)) %>%
  mutate(Species = "Roach")

perch_hourly <- perch_pike_distances_df %>%
  mutate(
    timestamp = as.POSIXct(timestamp, tz = "Europe/Stockholm"),
    Date      = as.Date(timestamp, tz = "Europe/Stockholm"),
    hour      = as.integer(format(timestamp, "%H"))
  ) %>%
  group_by(Prey_ID, Date, hour) %>%
  summarise(
    n_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    n_encounters_probable  = sum(encounter_probable,  na.rm = TRUE),
    n_encounters_combined  = sum(encounter_high_conf, na.rm = TRUE) +
      sum(encounter_probable,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(individual_ID = Prey_ID) %>%
  left_join(cow_meta, by = "individual_ID") %>%
  filter(!is.na(treatment)) %>%
  mutate(Species = "Perch")


# =================================================================-
# 4.5 Q1 PRE-FILTERING DATASET & MODELS ####
# =================================================================-
# Encounter rates estimated from the UNFILTERED dataset (before removal
# of post-predation/mortality tracks). High-confidence encounters only.
# Note: high-confidence = ≤2.0m (Cow Paradise threshold).

message("\n=== Q1 PRE-FILTERING ANALYSIS (COW PARADISE) ===")

## 4.5.1 Reload raw (unfiltered) encounter distances ----
roach_pike_raw <- readRDS(
  paste0(paths$encounters, "cow_pike_roach_distances_tiered_df.rds")
)
perch_pike_raw <- readRDS(
  paste0(paths$encounters, "cow_pike_perch_distances_tiered_df.rds")
)

## 4.5.2 Load pre-filtered metadata (04_cow_sub.rds) ----
cow_pre_data <- readRDS(
  paste0(paths$filtered, "04_cow_sub.rds")
)

cow_pre_meta <- cow_pre_data %>%
  select(individual_ID, treatment, weight) %>%
  distinct()

## 4.5.3 Build individual-level pre-filter totals ----
enc_roach_pre <- roach_pike_raw %>%
  rename(individual_ID = Prey_ID) %>%
  group_by(individual_ID) %>%
  summarise(
    total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    n_days                     = n_distinct(Date),
    .groups = "drop"
  ) %>%
  left_join(cow_pre_meta, by = "individual_ID") %>%
  filter(!is.na(treatment)) %>%
  mutate(
    Species        = "Roach",
    rate_high_conf = total_encounters_high_conf / n_days,
    treatment      = factor(treatment, levels = c("Control", "Mix")),
    individual_ID  = factor(individual_ID),
    weight_z       = as.numeric(scale(weight))
  )

enc_perch_pre <- perch_pike_raw %>%
  rename(individual_ID = Prey_ID) %>%
  group_by(individual_ID) %>%
  summarise(
    total_encounters_high_conf = sum(encounter_high_conf, na.rm = TRUE),
    n_days                     = n_distinct(Date),
    .groups = "drop"
  ) %>%
  left_join(cow_pre_meta, by = "individual_ID") %>%
  filter(!is.na(treatment)) %>%
  mutate(
    Species        = "Perch",
    rate_high_conf = total_encounters_high_conf / n_days,
    treatment      = factor(treatment, levels = c("Control", "Mix")),
    individual_ID  = factor(individual_ID),
    weight_z       = as.numeric(scale(weight))
  )

## 4.5.4 Summary statistics — pre-filter ----
enc_rate_summary_pre <- bind_rows(enc_roach_pre, enc_perch_pre) %>%
  group_by(Species, treatment) %>%
  summarise(
    n              = n(),
    mean_rate_hc   = mean(rate_high_conf,  na.rm = TRUE),
    sd_rate_hc     = sd(rate_high_conf,    na.rm = TRUE),
    median_rate_hc = median(rate_high_conf, na.rm = TRUE),
    iqr_rate_hc    = IQR(rate_high_conf,   na.rm = TRUE),
    total_enc_hc   = sum(total_encounters_high_conf, na.rm = TRUE),
    .groups = "drop"
  )

message("\nPre-filter encounter rate summary — Cow Paradise:")
print(enc_rate_summary_pre)

## 4.5.5 Roach — pre-filter high-confidence GLMM ----
enc_mod_roach_hc_pre <- glmmTMB(
  total_encounters_high_conf ~
    treatment +
    weight_z,
  data   = enc_roach_pre,
  family = nbinom2(link = "log"),
  offset = log(n_days)
)

print(summary(enc_mod_roach_hc_pre))

sim_res_roach_hc_pre <- simulateResiduals(enc_mod_roach_hc_pre, n = 500)
png(file.path(paths$figures, "glmm_dharma_roach_hc_pre.png"),
    width = 1200, height = 600)
plot(sim_res_roach_hc_pre)
dev.off()

emm_roach_hc_pre <- emmeans(enc_mod_roach_hc_pre,
                            specs  = ~ treatment,
                            type   = "response",
                            offset = log(1))
print(emm_roach_hc_pre)
print(pairs(emm_roach_hc_pre))

print_glmm_results(enc_mod_roach_hc_pre,
                   lake         = "COW PARADISE",
                   question     = "Q1 Pre-filter",
                   species      = "Roach",
                   response_var = "total_encounters_high_conf (<=2.0m)")
message("  Marginal means (encounters/day):")
print(emm_roach_hc_pre)
message("  Pairwise contrast (Mix vs Control):")
print(pairs(emm_roach_hc_pre))

## 4.5.6 Perch — pre-filter high-confidence GLMM ----
enc_mod_perch_hc_pre <- glmmTMB(
  total_encounters_high_conf ~
    treatment +
    weight_z,
  data   = enc_perch_pre,
  family = nbinom2(link = "log"),
  offset = log(n_days)
)

print(summary(enc_mod_perch_hc_pre))

sim_res_perch_hc_pre <- simulateResiduals(enc_mod_perch_hc_pre, n = 500)
png(file.path(paths$figures, "glmm_dharma_perch_hc_pre.png"),
    width = 1200, height = 600)
plot(sim_res_perch_hc_pre)
dev.off()

emm_perch_hc_pre <- emmeans(enc_mod_perch_hc_pre,
                            specs  = ~ treatment,
                            type   = "response",
                            offset = log(1))
print(emm_perch_hc_pre)
print(pairs(emm_perch_hc_pre))

print_glmm_results(enc_mod_perch_hc_pre,
                   lake         = "COW PARADISE",
                   question     = "Q1 Pre-filter",
                   species      = "Perch",
                   response_var = "total_encounters_high_conf (<=2.0m)")
message("  Marginal means (encounters/day):")
print(emm_perch_hc_pre)
message("  Pairwise contrast (Mix vs Control):")
print(pairs(emm_perch_hc_pre))


# =================================================================-
# 5. Q1: OVERALL ENCOUNTER RATE — NEGATIVE BINOMIAL GLMM ####
# =================================================================-
# Models run separately per species and for two encounter definitions:
#
#   (A) HIGH-CONFIDENCE ONLY: encounters <= 2.0m (Cow Paradise threshold).
#
#   (B) COMBINED (high-confidence + probable): encounters <= 4.0m.
#       Both definitions are run to assess robustness of treatment effects
#       to threshold choice.
#
# Single-lake analysis (Cow Paradise); individual as random intercept.

# ---- A: HIGH-CONFIDENCE ENCOUNTERS ONLY (<=2.0m) ----

## 5.1 Roach — high-confidence GLMM ----
enc_mod_roach_hc <- glmmTMB(
  total_encounters_high_conf ~
    treatment +
    weight_z,
  data   = enc_roach_total,
  family = nbinom2(link = "log"),
  offset = log(n_days)
)

print(summary(enc_mod_roach_hc))

sim_res_roach_hc <- simulateResiduals(enc_mod_roach_hc, n = 500)
png(file.path(paths$figures, "glmm_dharma_roach_high_conf.png"),
    width = 1200, height = 600)
plot(sim_res_roach_hc)
dev.off()

emm_roach_hc <- emmeans(enc_mod_roach_hc,
                        specs  = ~ treatment,
                        type   = "response",
                        offset = log(1))
print(emm_roach_hc)

print(pairs(emm_roach_hc))

print_glmm_results(enc_mod_roach_hc,
                   lake         = "COW PARADISE",
                   question     = "Q1 Post-filter",
                   species      = "Roach",
                   response_var = "total_encounters_high_conf (<=2.0m)")
message("  Marginal means (encounters/day):")
print(emm_roach_hc)
message("  Pairwise contrast (Mix vs Control):")
print(pairs(emm_roach_hc))

## 5.2 Perch — high-confidence GLMM ----

enc_mod_perch_hc <- glmmTMB(
  total_encounters_high_conf ~
    treatment +
    weight_z,
  data   = enc_perch_total,
  family = nbinom2(link = "log"),
  offset = log(n_days)
)

print(summary(enc_mod_perch_hc))

sim_res_perch_hc <- simulateResiduals(enc_mod_perch_hc, n = 500)
png(file.path(paths$figures, "glmm_dharma_perch_high_conf.png"),
    width = 1200, height = 600)
plot(sim_res_perch_hc)
dev.off()

emm_perch_hc <- emmeans(enc_mod_perch_hc,
                        specs  = ~ treatment,
                        type   = "response",
                        offset = log(1))

print(emm_perch_hc)

print(pairs(emm_perch_hc))

print_glmm_results(enc_mod_perch_hc,
                   lake         = "COW PARADISE",
                   question     = "Q1 Post-filter",
                   species      = "Perch",
                   response_var = "total_encounters_high_conf (<=2.0m)")
message("  Marginal means (encounters/day):")
print(emm_perch_hc)
message("  Pairwise contrast (Mix vs Control):")
print(pairs(emm_perch_hc))

# ---- B: COMBINED ENCOUNTERS (high-confidence + probable, <=4.0m) ----


## 5.3 Roach — combined GLMM ----

enc_mod_roach_comb <- glmmTMB(
  total_encounters_combined ~
    treatment +
    weight_z,
  data   = enc_roach_total,
  family = nbinom2(link = "log"),
  offset = log(n_days)
)

print(summary(enc_mod_roach_comb))

sim_res_roach_comb <- simulateResiduals(enc_mod_roach_comb, n = 500)
png(file.path(paths$figures, "glmm_dharma_roach_combined.png"),
    width = 1200, height = 600)
plot(sim_res_roach_comb)
dev.off()

emm_roach_comb <- emmeans(enc_mod_roach_comb,
                          specs  = ~ treatment,
                          type   = "response",
                          offset = log(1))

print(emm_roach_comb)
print(pairs(emm_roach_comb))

print_glmm_results(enc_mod_roach_comb,
                   lake         = "COW PARADISE",
                   question     = "Q1 Post-filter",
                   species      = "Roach",
                   response_var = "total_encounters_combined (<=4.0m)")
message("  Marginal means (encounters/day):")
print(emm_roach_comb)
message("  Pairwise contrast (Mix vs Control):")
print(pairs(emm_roach_comb))

## 5.4 Perch — combined GLMM ----

enc_mod_perch_comb <- glmmTMB(
  total_encounters_combined ~
    treatment +
    weight_z,  
  data   = enc_perch_total,
  family = nbinom2(link = "log"),
  offset = log(n_days)
)

print(summary(enc_mod_perch_comb))

sim_res_perch_comb <- simulateResiduals(enc_mod_perch_comb, n = 500)
png(file.path(paths$figures, "glmm_dharma_perch_combined.png"),
    width = 1200, height = 600)
plot(sim_res_perch_comb)
dev.off()

emm_perch_comb <- emmeans(enc_mod_perch_comb,
                          specs  = ~ treatment,
                          type   = "response",
                          offset = log(1))
print(emm_perch_comb)
print(pairs(emm_perch_comb))

print_glmm_results(enc_mod_perch_comb,
                   lake         = "COW PARADISE",
                   question     = "Q1 Post-filter",
                   species      = "Perch",
                   response_var = "total_encounters_combined (<=4.0m)")
message("  Marginal means (encounters/day):")
print(emm_perch_comb)
message("  Pairwise contrast (Mix vs Control):")
print(pairs(emm_perch_comb))

## 5.5 Encounter rate plots ----
enc_colours <- c("Control" = "#004D40", "Mix" = "#D81B60")

p_enc_hc <- ggplot(
  bind_rows(enc_roach_total, enc_perch_total) %>%
    mutate(rate = total_encounters_high_conf / n_days),
  aes(x = treatment, y = rate, fill = treatment)
) +
  geom_boxplot(alpha = 0.75, outlier.size = 1) +
  scale_fill_manual(values = enc_colours, name = "Treatment") +
  facet_wrap(~Species, scales = "free_y") +
  labs(x = NULL, y = "Encounters per day") +
  theme_classic() +
  theme(
    legend.position  = "none",
    axis.title.y     = element_text(size = 12),
    axis.title.x     = element_text(size = 12),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 10),
    strip.text       = element_text(size = 12),
    strip.background = element_rect(fill = "grey95", colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA)
  )

print(p_enc_hc)

ggsave(file.path(paths$figures, "encounter_rates_high_conf.pdf"),
       p_enc_hc, width = 8, height = 5, dpi = 300)

p_enc_comb <- ggplot(
  bind_rows(enc_roach_total, enc_perch_total) %>%
    mutate(rate = total_encounters_combined / n_days),
  aes(x = treatment, y = rate, fill = treatment)
) +
  geom_boxplot(alpha = 0.75, outlier.size = 1) +
  scale_fill_manual(values = enc_colours) +
  facet_wrap(~Species, scales = "free_y") +
  labs(x = NULL, y = "Encounters per day") +
  theme_classic() +
  theme(
    legend.position  = "none",
    axis.title.y     = element_text(size = 12),
    axis.title.x     = element_text(size = 12),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 10),
    strip.text       = element_text(size = 12),
    strip.background = element_rect(fill = "grey95", colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA)
  )

ggsave(file.path(paths$figures, "cow_encounter_rates_combined.png"),
       p_enc_comb, width = 8, height = 5, dpi = 300)


# =================================================================-
# 6. Q2 & Q3: DIEL ENCOUNTER PATTERNS — GAMMs ####
# =================================================================-
# Same model structure as BT and Muddyfoot. Encounter thresholds
# (2.0m / 4.0m) are already embedded in the
# encounter_high_conf / encounter_probable columns.
# Single-lake analysis (Cow Paradise).

## 6.1 Prepare hourly datasets ----
# Experiment start date for Cow Paradise: 2022-09-27.
roach_h <- roach_hourly %>%
  mutate(
    exp_day = as.integer(
      as.duration(as.POSIXct(paste(Date, "00:00:00")) -
                    as.POSIXct("2022-09-27")) / lubridate::ddays(1)
    ),
    treatment     = factor(treatment, levels = c("Control", "Mix")),
    individual_ID = factor(individual_ID),
    hour          = as.integer(hour),
    exp_day_z     = as.numeric(scale(exp_day)),
    weight_z      = as.numeric(scale(weight))
  ) %>%
  filter(!is.na(exp_day)) %>%
  arrange(individual_ID, exp_day, hour)

perch_h <- perch_hourly %>%
  mutate(
    exp_day = as.integer(
      as.duration(as.POSIXct(paste(Date, "00:00:00")) -
                    as.POSIXct("2022-09-27")) / lubridate::ddays(1)
    ),
    treatment     = factor(treatment, levels = c("Control", "Mix")),
    individual_ID = factor(individual_ID),
    hour          = as.integer(hour),
    exp_day_z     = as.numeric(scale(exp_day)),
    weight_z      = as.numeric(scale(weight))
  ) %>%
  filter(!is.na(exp_day)) %>%
  arrange(individual_ID, exp_day, hour)

## 6.2 Helper: add AR1 start flag ----
add_ar_start <- function(data) {
  data %>%
    arrange(individual_ID, exp_day, hour) %>%
    mutate(start = is.na(lag(individual_ID)) | individual_ID != lag(individual_ID))
}

## 6.3 Helper: estimate lag-1 autocorrelation (rho) ----
estimate_rho <- function(data, mod) {
  data$residuals <- resid(mod, type = "pearson")
  data %>%
    group_split(individual_ID) %>%
    map_df(~{
      acf_res <- acf(.x$residuals, plot = FALSE)
      tibble(individual_ID = unique(.x$individual_ID),
             lag1_value    = acf_res$acf[2])
    }) %>%
    summarise(mean_rho = mean(lag1_value, na.rm = TRUE)) %>%
    pull(mean_rho)
}

## 6.4 Shared ggplot elements for diel plots ----
enc_x_scale <- scale_x_continuous(
  breaks = c(3, 6, 9, 12, 15, 18, 21),
  labels = c("03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00")
)

add_diel_lines <- function(p) {
  p +
    geom_vline(xintercept = nauticalDawn, linetype = "dashed", alpha = 0.45) +
    geom_vline(xintercept = sunrise,      linetype = "dashed", alpha = 0.45) +
    geom_vline(xintercept = sunset,       linetype = "dashed", alpha = 0.45) +
    geom_vline(xintercept = nauticalDusk, linetype = "dashed", alpha = 0.45)
}

plot_diel_enc <- function(pred_summary, species) {
  p <- ggplot(pred_summary, aes(x = hour, y = encounters,
                                group = treatment, colour = treatment)) +
    geom_ribbon(aes(x = hour, ymin = lower, ymax = upper, group = treatment),
                fill = "black", alpha = 0.1, inherit.aes = FALSE) +
    geom_line(linewidth = 1.8) +
    scale_colour_manual(values = enc_colours, name = "Treatment") +
    enc_x_scale +
    labs(
         y     = "Predicted encounters\nwith pike per hour\n",
         x     = "\nTime of day") +
    theme_few() +
    theme(legend.position = "none")
  add_diel_lines(p)
}

plot_enc_contrast <- function(contrasts_df, species) {
  p <- ggplot(contrasts_df, aes(x = hour, y = estimate,
                                ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(fill = "grey80", alpha = 0.5) +
    geom_line(linewidth = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    enc_x_scale +
    labs(title = species,
         y     = "Predicted difference\nbetween treatments\n(Mix – Control)\n",
         x     = "\nTime of day") +
    theme_few()
  add_diel_lines(p)
}

# ----------------------------------------------------------------
# ROACH GAMM
# ----------------------------------------------------------------

## 6.5 Roach — first-pass model to estimate rho ----
roach_h <- roach_h %>% add_ar_start()

roach_gamm_init <- bam(
  n_encounters_high_conf ~
    treatment + exp_day_z + weight_z +
    s(hour, by = treatment, k = 15, bs = "cc") +
    s(hour, individual_ID, bs = "fs", xt = list(bs = "cc")),
  data     = roach_h,
  method   = "fREML",
  discrete = TRUE,
  knots    = list(hour = c(0, 24)),
  family   = nb(link = "log")
)

roach_rho <- estimate_rho(roach_h, roach_gamm_init)

## 6.6 Roach — AR-corrected GAMM ----
roach_mod <- bam(
  n_encounters_high_conf ~
    treatment + exp_day_z + weight_z +
    s(hour, k = 15, bs = "cc") +          # single shared smooth
    s(hour, individual_ID, bs = "fs",
      xt = list(bs = "cc")),
  data     = roach_h,
  method   = "fREML",
  rho      = round(roach_rho, 1),
  AR.start = roach_h$start,
  discrete = TRUE,
  knots    = list(hour = c(0, 24)),
  family   = nb(link = "log")
)

## 6.7 Roach — diagnostics, summary, anova ----
check_resid(roach_mod, ask = FALSE)
print(appraise(roach_mod))
summary(roach_mod)
anova.gam(roach_mod)

print_gamm_results(roach_mod, lake = "COW PARADISE", species = "Roach")

## 6.8 Roach — predictions ----
newdata_roach <- tidyr::expand(roach_h, nesting(individual_ID, treatment),
                               hour = seq(0, 23, length.out = 96)) %>%
  mutate(exp_day_z = 0, weight_z = 0)

ilink_roach <- family(roach_mod)$linkinv

roach_pred_summary <- bind_cols(
  newdata_roach,
  as.data.frame(predict(roach_mod, newdata = newdata_roach,
                        type = "link", se.fit = TRUE))
) %>%
  transform(fitted = ilink_roach(fit),
            upper  = ilink_roach(fit + 2 * se.fit),
            lower  = ilink_roach(fit - 2 * se.fit)) %>%
  group_by(treatment, hour) %>%
  summarise(encounters = mean(fitted), upper = mean(upper),
            lower = mean(lower), .groups = "drop")

## 6.9 Roach — treatment contrasts ----
roach_comp <- comparisons(roach_mod, newdata = newdata_roach,
                          variables = "treatment", type = "response",
                          vcov = TRUE, by = "hour")
print(
  roach_comp %>%
    mutate(hour_int = round(hour, 0)) %>%
    group_by(hour_int) %>%
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high), .groups = "drop") %>%
    mutate(abs_diff = abs(estimate)) %>%
    filter(abs_diff == max(abs_diff) | abs_diff == min(abs_diff))
)

## 6.9b Roach — diel period peak analysis ----
newdata_roach_periods <- newdata_roach %>%
  mutate(
    diel_period = factor(
      case_when(
        hour >= nauticalDawn & hour < sunrise      ~ "Dawn",
        hour >= sunrise      & hour < sunset       ~ "Day",
        hour >= sunset       & hour < nauticalDusk ~ "Dusk",
        TRUE                                        ~ "Night"
      ),
      levels = c("Night", "Dawn", "Day", "Dusk")
    )
  )

roach_period_raw <- bind_cols(
  newdata_roach_periods,
  as.data.frame(predict(roach_mod, newdata = newdata_roach_periods,
                        type = "link", se.fit = TRUE))
) %>%
  mutate(
    fitted = ilink_roach(fit),
    upper  = ilink_roach(fit + 1.96 * se.fit),
    lower  = ilink_roach(fit - 1.96 * se.fit)
  )

roach_period_summary <- roach_period_raw %>%
  group_by(diel_period) %>%
  summarise(mean_encounters = mean(fitted),
            lower_95ci      = mean(lower),
            upper_95ci      = mean(upper), .groups = "drop") %>%
  arrange(desc(mean_encounters))

print(roach_period_summary)

roach_period_trt_summary <- roach_period_raw %>%
  group_by(treatment, diel_period) %>%
  summarise(mean_encounters = mean(fitted),
            lower_95ci      = mean(lower),
            upper_95ci      = mean(upper), .groups = "drop") %>%
  arrange(treatment, desc(mean_encounters))

message("\nRoach — mean predicted encounters by diel period and treatment:")
print(roach_period_trt_summary)

roach_period_preds <- avg_predictions(
  roach_mod, newdata = newdata_roach_periods,
  by = "diel_period", type = "response"
)

print(roach_period_preds)

roach_period_pairs <- hypotheses(
  roach_period_preds,
  hypothesis = c(
    "b1 - b2 = 0",
    "b1 - b3 = 0",
    "b1 - b4 = 0",
    "b2 - b3 = 0",
    "b2 - b4 = 0",
    "b3 - b4 = 0"
  )
)

print(roach_period_pairs)

## 6.10 Roach — diel encounter figure ----
roach_diel_plot     <- plot_diel_enc(roach_pred_summary, "Roach")

print(roach_diel_plot)

ggsave(file.path(paths$figures, "cow_roach_diel_encounter_figure.pdf"),
       roach_diel_plot, width = 12, height = 5, dpi = 300)

# ----------------------------------------------------------------
# PERCH GAMM
# ----------------------------------------------------------------


## 6.11 Perch — first-pass model to estimate rho ----
perch_h <- perch_h %>% add_ar_start()

perch_gamm_init <- bam(
  n_encounters_high_conf ~
    treatment + exp_day_z + weight_z +
    s(hour, by = treatment, k = 15, bs = "cc") +
    s(hour, individual_ID, bs = "fs", xt = list(bs = "cc")),
  data     = perch_h,
  method   = "fREML",
  discrete = TRUE,
  knots    = list(hour = c(0, 24)),
  family   = nb(link = "log")
)

perch_rho <- estimate_rho(perch_h, perch_gamm_init)
message(sprintf("Perch lag-1 autocorrelation (rho): %.2f", perch_rho))

## 6.12 Perch — AR-corrected GAMM ----
perch_mod <- bam(
  n_encounters_high_conf ~
    treatment + exp_day_z + weight_z +
    s(hour, by = treatment, k = 15, bs = "cc") +
    s(hour, individual_ID, bs = "fs", xt = list(bs = "cc")),
  data     = perch_h,
  method   = "fREML",
  rho      = round(perch_rho, 1),
  AR.start = perch_h$start,
  discrete = TRUE,
  knots    = list(hour = c(0, 24)),
  family   = nb(link = "log")
)

## 6.13 Perch — diagnostics, summary, anova ----
check_resid(perch_mod, ask = FALSE)
print(appraise(perch_mod))
summary(perch_mod)
anova.gam(perch_mod)

print_gamm_results(perch_mod, lake = "COW PARADISE", species = "Perch")

## 6.14 Perch — predictions ----
newdata_perch <- tidyr::expand(perch_h, nesting(individual_ID, treatment),
                               hour = seq(0, 23, length.out = 96)) %>%
  mutate(exp_day_z = 0, weight_z = 0)

ilink_perch <- family(perch_mod)$linkinv

perch_pred_summary <- bind_cols(
  newdata_perch,
  as.data.frame(predict(perch_mod, newdata = newdata_perch,
                        type = "link", se.fit = TRUE))
) %>%
  transform(fitted = ilink_perch(fit),
            upper  = ilink_perch(fit + 2 * se.fit),
            lower  = ilink_perch(fit - 2 * se.fit)) %>%
  group_by(treatment, hour) %>%
  summarise(encounters = mean(fitted), upper = mean(upper),
            lower = mean(lower), .groups = "drop")

## 6.15 Perch — treatment contrasts ----
perch_comp <- comparisons(perch_mod, newdata = newdata_perch,
                          variables = "treatment", type = "response",
                          vcov = TRUE, by = "hour")

message("\nPerch — peak/trough treatment difference:")
print(
  perch_comp %>%
    mutate(hour_int = round(hour, 0)) %>%
    group_by(hour_int) %>%
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high), .groups = "drop") %>%
    mutate(abs_diff = abs(estimate)) %>%
    filter(abs_diff == max(abs_diff) | abs_diff == min(abs_diff))
)

## 6.15b Perch — diel period peak analysis ----
newdata_perch_periods <- newdata_perch %>%
  mutate(
    diel_period = factor(
      case_when(
        hour >= nauticalDawn & hour < sunrise      ~ "Dawn",
        hour >= sunrise      & hour < sunset       ~ "Day",
        hour >= sunset       & hour < nauticalDusk ~ "Dusk",
        TRUE                                        ~ "Night"
      ),
      levels = c("Night", "Dawn", "Day", "Dusk")
    )
  )

perch_period_raw <- bind_cols(
  newdata_perch_periods,
  as.data.frame(predict(perch_mod, newdata = newdata_perch_periods,
                        type = "link", se.fit = TRUE))
) %>%
  mutate(
    fitted = ilink_perch(fit),
    upper  = ilink_perch(fit + 1.96 * se.fit),
    lower  = ilink_perch(fit - 1.96 * se.fit)
  )

perch_period_summary <- perch_period_raw %>%
  group_by(diel_period) %>%
  summarise(mean_encounters = mean(fitted),
            lower_95ci      = mean(lower),
            upper_95ci      = mean(upper), .groups = "drop") %>%
  arrange(desc(mean_encounters))

print(perch_period_summary)

perch_period_trt_summary <- perch_period_raw %>%
  group_by(treatment, diel_period) %>%
  summarise(mean_encounters = mean(fitted),
            lower_95ci      = mean(lower),
            upper_95ci      = mean(upper), .groups = "drop") %>%
  arrange(treatment, desc(mean_encounters))

message("\nPerch — mean predicted encounters by diel period and treatment:")
print(perch_period_trt_summary)

perch_period_preds <- avg_predictions(
  perch_mod, newdata = newdata_perch_periods,
  by = "diel_period", type = "response"
)

message("\nPerch — marginal mean encounters per diel period:")
print(perch_period_preds)

perch_period_pairs <- hypotheses(
  perch_period_preds,
  hypothesis = c(
    "b1 - b2 = 0",
    "b1 - b3 = 0",
    "b1 - b4 = 0",
    "b2 - b3 = 0",
    "b2 - b4 = 0",
    "b3 - b4 = 0"
  )
)

message("\nPerch — pairwise diel period contrasts:")
print(perch_period_pairs)

## 6.16 Perch — diel encounter figure ----
perch_diel_plot     <- plot_diel_enc(perch_pred_summary, "Perch")

plot(perch_diel_plot)

ggsave(file.path(paths$figures, "cow_perch_diel_encounter_figure.pdf"),
       perch_diel_plot, width = 12, height = 5, dpi = 300)

## 6.17 Combined species overview ----
enc_diel_summary <- bind_rows(
  roach_pred_summary %>% mutate(Species = "Roach"),
  perch_pred_summary %>% mutate(Species = "Perch")
)

p_enc_diel_overview <- ggplot(enc_diel_summary,
                              aes(x = hour, y = encounters,
                                  colour = treatment, group = treatment)) +
  geom_ribbon(aes(x = hour, ymin = lower, ymax = upper, group = treatment),
              fill = "black", alpha = 0.08, inherit.aes = FALSE) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = enc_colours, name = "Treatment") +
  facet_wrap(~Species, scales = "free_y") +
  enc_x_scale +
  geom_vline(xintercept = nauticalDawn, linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = sunrise,      linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = sunset,       linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = nauticalDusk, linetype = "dashed", alpha = 0.4) +
  labs(
    y = "Predicted encounters with pike per hour\n",
    x = "\nTime of day"
  ) +
  theme_few()

ggsave(file.path(paths$figures, "cow_combined_diel_encounter_overview.png"),
       p_enc_diel_overview, width = 12, height = 5, dpi = 300)
message("Saved: cow_combined_diel_encounter_overview.png")


# =================================================================-
# 7. SAVE OUTPUTS ####
# =================================================================-

# Pre-filter models
saveRDS(enc_mod_roach_hc_pre, file.path(paths$tables, "glmm_roach_hc_pre.rds"))
saveRDS(enc_mod_perch_hc_pre, file.path(paths$tables, "glmm_perch_hc_pre.rds"))

# Post-filter models
saveRDS(enc_mod_roach_hc,   file.path(paths$tables, "glmm_roach_high_conf.rds"))
saveRDS(enc_mod_perch_hc,   file.path(paths$tables, "glmm_perch_high_conf.rds"))
saveRDS(enc_mod_roach_comb, file.path(paths$tables, "glmm_roach_combined.rds"))
saveRDS(enc_mod_perch_comb, file.path(paths$tables, "glmm_perch_combined.rds"))
saveRDS(roach_mod,          file.path(paths$tables, "gamm_roach_encounter_diel.rds"))
saveRDS(perch_mod,          file.path(paths$tables, "gamm_perch_encounter_diel.rds"))

# Pre-filter summary data
write.csv(enc_rate_summary_pre,
          file.path(paths$tables, "encounter_rate_summary_pre_filter.csv"),
          row.names = FALSE)
write.csv(bind_rows(enc_roach_pre, enc_perch_pre),
          file.path(paths$tables, "individual_total_encounters_pre_filter.csv"),
          row.names = FALSE)

# Post-filter summary data
write.csv(roach_pred_summary,
          file.path(paths$tables, "roach_diel_encounter_predictions.csv"),
          row.names = FALSE)
write.csv(perch_pred_summary,
          file.path(paths$tables, "perch_diel_encounter_predictions.csv"),
          row.names = FALSE)
write.csv(bind_rows(enc_roach_total, enc_perch_total),
          file.path(paths$tables, "individual_total_encounters_with_treatment.csv"),
          row.names = FALSE)

message("\n=== ANALYSIS COMPLETE ===")
message(sprintf("Figures saved to: %s", paths$figures))
message(sprintf("Tables saved to:  %s", paths$tables))


# =================================================================-
# 8. SESSION INFO ####
# =================================================================-

sink(file.path(paths$tables, "session_info_encounter_analysis.txt"))
cat("Encounter Rate Analysis — Cow Paradise — Session Info\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")
print(sessionInfo())
sink()
