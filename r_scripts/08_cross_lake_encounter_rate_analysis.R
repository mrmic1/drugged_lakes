# =================================================================-
# Cross-Lake Encounter Rate Analysis ####
# =================================================================-
#
# ANALYSIS OVERVIEW:
# Combines encounter rate data from all three lakes (BT, Muddyfoot,
# Cow Paradise) to test for overarching treatment effects on predator
# encounter rates across lakes. Lake is included as a random effect
# to account for between-lake variation in baseline encounter rates,
# GPS error magnitude, and environmental conditions.
#
# Addresses:
# Q1: Is there a consistent cross-lake effect of treatment (Control vs.
#     Mix) on overall encounter rate with pike?
#
# Q2: Is there a consistent cross-lake diel pattern of encounters?
# Q3: Does the diel pattern of encounters differ between treatments,
#     consistently across lakes?

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
  library(tidyr)           # Data reshaping
  library(purrr)           # map_df for ACF
  library(ggplot2)         # Plotting
  library(ggthemes)        # Plot themes
  library(ggpubr)          # Multi-panel figures
  library(mgcv)            # GAMMs (bam)
  library(gratia)          # GAMM diagnostics
  library(itsadug)         # GAMM residual checks
  library(marginaleffects) # Treatment contrasts
  library(glmmTMB)         # Negative binomial GLMM
  library(DHARMa)          # Residual diagnostics
  library(emmeans)         # Marginal means and contrasts
  library(suncalc)         # Diel annotation
})

Sys.setenv(TZ = "Europe/Stockholm")

## 1.2 Define Directory Paths ----
paths <- list(
  enc_bt  = "./data/encounters/BT/",
  enc_mf  = "./data/encounters/muddyfoot/",
  enc_cow = "./data/encounters/cow_paradise/",
  filt_bt  = "./data/tracks_filtered/BT/",
  filt_mf  = "./data/tracks_filtered/muddyfoot/",
  filt_cow = "./data/tracks_filtered/cow_paradise/",
  mortality = "./data/encounters/suspected_mortality_updated.xlsx",
  figures   = "./figures/cross_lake/encounter_analysis/",
  tables    = "./tables/cross_lake/encounter_analysis/"
)

for (p in c(paths$figures, paths$tables)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

## 1.3 Diel Reference Times ----
# All three lakes are in northern Sweden (~lat 63.77, lon 20.05).
nauticalDawn <- 5.494649   # ~05:30
sunrise      <- 7.189919   # ~07:11
sunset       <- 17.75258   # ~17:45
nauticalDusk <- 19.44784   # ~19:27

## 1.4 Plot colours ----
enc_colours <- c("Control" = "#004D40", "Mix" = "#D81B60")


# =================================================================-
# 2. DATA LOADING ####
# =================================================================-

## 2.1 Load encounter distance data ----
message("Loading encounter distance data...")

# BT
roach_bt  <- readRDS(file.path(paths$enc_bt, "BT_pike_roach_distances_tiered_df.rds"))
perch_bt  <- readRDS(file.path(paths$enc_bt, "BT_pike_perch_distances_tiered_df.rds"))

# Muddyfoot
roach_mf  <- readRDS(file.path(paths$enc_mf, "muddyfoot_pike_roach_distances_tiered_df.rds"))
perch_mf  <- readRDS(file.path(paths$enc_mf, "muddyfoot_pike_perch_distances_tiered_df.rds"))

# Cow Paradise
roach_cow <- readRDS(file.path(paths$enc_cow, "cow_pike_roach_distances_tiered_df.rds"))
perch_cow <- readRDS(file.path(paths$enc_cow, "cow_pike_perch_distances_tiered_df.rds"))

## 2.2 Load metadata ----
meta_bt  <- readRDS(file.path(paths$filt_bt,  "05_BT_sub.rds")) %>%
  select(individual_ID, treatment, weight) %>% distinct()

meta_mf  <- readRDS(file.path(paths$filt_mf,  "05_muddyfoot_sub.rds")) %>%
  select(individual_ID, treatment, weight) %>% distinct()

meta_cow <- readRDS(file.path(paths$filt_cow, "05_cow_sub.rds")) %>%
  select(individual_ID, treatment, weight) %>% distinct()


# =================================================================-
# 3. FILTER OUT POST-PREDATION ENCOUNTER RECORDS ####
# =================================================================-

## 3.1 Load mortality reference ----
mortality_preds <- readxl::read_excel(paths$mortality)

get_death_dates <- function(lake_name) {
  mortality_preds %>%
    filter(lake == lake_name,
           species %in% c("Roach", "Perch")) %>%
    select(individual_ID, revised_likely_death_date) %>%
    rename(death_date = revised_likely_death_date) %>%
    mutate(death_date = as.Date(death_date, origin = "1970-01-01"))
}

filter_post_predation <- function(df, death_dates) {
  df %>%
    left_join(death_dates, by = c("Prey_ID" = "individual_ID")) %>%
    filter(is.na(death_date) | Date < death_date) %>%
    select(-death_date)
}

death_bt  <- get_death_dates("BT")
death_mf  <- get_death_dates("muddyfoot")
death_cow <- get_death_dates("cow paradise")

roach_bt  <- filter_post_predation(roach_bt,  death_bt)
perch_bt  <- filter_post_predation(perch_bt,  death_bt)
roach_mf  <- filter_post_predation(roach_mf,  death_mf)
perch_mf  <- filter_post_predation(perch_mf,  death_mf)
roach_cow <- filter_post_predation(roach_cow, death_cow)
perch_cow <- filter_post_predation(perch_cow, death_cow)

message("Post-predation filtering complete.")


# =================================================================-
# 4. BUILD COMBINED ANALYSIS DATASETS ####
# =================================================================-

## 4.1 Helper: individual-total encounter summary ----
make_enc_total <- function(distances_df, meta, lake_name, species_name) {
  distances_df %>%
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
    left_join(meta, by = "individual_ID") %>%
    filter(!is.na(treatment)) %>%
    mutate(
      Species       = species_name,
      Lake          = lake_name,
      treatment     = factor(treatment, levels = c("Control", "Mix")),
      individual_ID = factor(individual_ID),
      weight_z      = as.numeric(scale(weight))
    )
}

## 4.2 Build per-lake totals and combine ----
enc_roach_total <- bind_rows(
  make_enc_total(roach_bt,  meta_bt,  "BT",           "Roach"),
  make_enc_total(roach_mf,  meta_mf,  "Muddyfoot",    "Roach"),
  make_enc_total(roach_cow, meta_cow, "Cow Paradise",  "Roach")
) %>% mutate(Lake = factor(Lake))

enc_perch_total <- bind_rows(
  make_enc_total(perch_bt,  meta_bt,  "BT",           "Perch"),
  make_enc_total(perch_mf,  meta_mf,  "Muddyfoot",    "Perch"),
  make_enc_total(perch_cow, meta_cow, "Cow Paradise",  "Perch")
) %>% mutate(Lake = factor(Lake))

message(sprintf("Combined totals: %d roach, %d perch individuals across lakes",
                nrow(enc_roach_total), nrow(enc_perch_total)))
message("\nRoach by lake:")
print(table(enc_roach_total$Lake, enc_roach_total$treatment))
message("\nPerch by lake:")
print(table(enc_perch_total$Lake, enc_perch_total$treatment))

## 4.3 Summary statistics — encounter rates by lake, species, treatment ----
enc_rate_summary <- bind_rows(enc_roach_total, enc_perch_total) %>%
  mutate(
    rate_high_conf = total_encounters_high_conf / n_days,
    rate_combined  = total_encounters_combined  / n_days
  ) %>%
  group_by(Lake, Species, treatment) %>%
  summarise(
    n                = n(),
    mean_rate_hc     = mean(rate_high_conf,  na.rm = TRUE),
    sd_rate_hc       = sd(rate_high_conf,    na.rm = TRUE),
    median_rate_hc   = median(rate_high_conf, na.rm = TRUE),
    mean_rate_comb   = mean(rate_combined,   na.rm = TRUE),
    sd_rate_comb     = sd(rate_combined,     na.rm = TRUE),
    median_rate_comb = median(rate_combined,  na.rm = TRUE),
    .groups = "drop"
  )

message("\nCross-lake encounter rate summary:")
print(enc_rate_summary)

## 4.4 Helper: hourly encounter dataset ----
# Experiment start dates: BT = 2022-09-25, Muddyfoot = 2022-09-25, Cow = 2022-09-27
exp_starts <- c(BT = "2022-09-25", Muddyfoot = "2022-09-25", "Cow Paradise" = "2022-09-27")

make_hourly <- function(distances_df, meta, lake_name, species_name) {
  distances_df %>%
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
        sum(encounter_probable, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(individual_ID = Prey_ID) %>%
    left_join(meta, by = "individual_ID") %>%
    filter(!is.na(treatment)) %>%
    mutate(
      Species  = species_name,
      Lake     = lake_name,
      exp_day  = as.integer(
        as.duration(as.POSIXct(paste(Date, "00:00:00")) -
                      as.POSIXct(exp_starts[lake_name])) / lubridate::ddays(1)
      )
    )
}

roach_hourly_all <- bind_rows(
  make_hourly(roach_bt,  meta_bt,  "BT",           "Roach"),
  make_hourly(roach_mf,  meta_mf,  "Muddyfoot",    "Roach"),
  make_hourly(roach_cow, meta_cow, "Cow Paradise",  "Roach")
)

perch_hourly_all <- bind_rows(
  make_hourly(perch_bt,  meta_bt,  "BT",           "Perch"),
  make_hourly(perch_mf,  meta_mf,  "Muddyfoot",    "Perch"),
  make_hourly(perch_cow, meta_cow, "Cow Paradise",  "Perch")
)


# =================================================================-
# 5. Q1: CROSS-LAKE OVERALL ENCOUNTER RATE — NEGATIVE BINOMIAL GLMM ####
# =================================================================-
# Models run separately per species and for two encounter definitions.
#
# Random effects structure:
#   (1 | Lake / individual_ID)  — individuals nested within lakes.
#   This accounts for:
#     (a) between-lake variation in baseline encounter rates
#     (b) between-individual variation within each lake
#
# The offset log(n_days) converts raw counts to daily rates, enabling
# comparison across individuals and lakes tracked for different durations.

message("\n=== Q1: CROSS-LAKE ENCOUNTER RATE GLMMs ===")

# ---- A: HIGH-CONFIDENCE ENCOUNTERS ONLY ----

message("\n--- A: High-confidence encounters only ---")
message("(≤1.4m for BT/Muddyfoot; ≤2.0m for Cow Paradise)")

## 5.1 Roach — high-confidence ----
message("\nRoach:")

enc_mod_roach_hc <- glmmTMB(
  total_encounters_high_conf ~
    treatment +
    weight_z +
    (1 | Lake / individual_ID),
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
message("\nRoach marginal means — high-confidence (encounters / day):")
print(emm_roach_hc)
message("\nRoach contrast — high-confidence (Mix vs Control):")
print(pairs(emm_roach_hc))

## 5.2 Perch — high-confidence ----
message("\nPerch:")

enc_mod_perch_hc <- glmmTMB(
  total_encounters_high_conf ~
    treatment +
    weight_z +
    (1 | Lake / individual_ID),
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
message("\nPerch marginal means — high-confidence (encounters / day):")
print(emm_perch_hc)
message("\nPerch contrast — high-confidence (Mix vs Control):")
print(pairs(emm_perch_hc))

# ---- B: COMBINED ENCOUNTERS ----

message("\n--- B: Combined encounters (high-confidence + probable) ---")
message("(≤2.8m for BT/Muddyfoot; ≤4.0m for Cow Paradise)")

## 5.3 Roach — combined ----
message("\nRoach:")

enc_mod_roach_comb <- glmmTMB(
  total_encounters_combined ~
    treatment +
    weight_z +
    (1 | Lake / individual_ID),
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
message("\nRoach marginal means — combined (encounters / day):")
print(emm_roach_comb)
message("\nRoach contrast — combined (Mix vs Control):")
print(pairs(emm_roach_comb))

## 5.4 Perch — combined ----
message("\nPerch:")

enc_mod_perch_comb <- glmmTMB(
  total_encounters_combined ~
    treatment +
    weight_z +
    (1 | Lake / individual_ID),
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
message("\nPerch marginal means — combined (encounters / day):")
print(emm_perch_comb)
message("\nPerch contrast — combined (Mix vs Control):")
print(pairs(emm_perch_comb))

## 5.5 Cross-lake encounter rate plot ----
# Faceted by Lake x Species to show between-lake variation.
p_enc_hc_cross <- bind_rows(enc_roach_total, enc_perch_total) %>%
  mutate(rate = total_encounters_high_conf / n_days) %>%
  ggplot(aes(x = treatment, y = rate, fill = treatment)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1) +
  scale_fill_manual(values = enc_colours, name = "Treatment") +
  facet_grid(Species ~ Lake, scales = "free_y") +
  labs(x = NULL, y = "Encounters per day") +
  theme_classic(base_family = "Times New Roman") +
  theme(
    legend.position  = "none",
    axis.title.y     = element_text(size = 12),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 10),
    strip.text       = element_text(size = 11, family = "Times New Roman"),
    strip.background = element_rect(fill = "grey95", colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA)
  )

ggsave(file.path(paths$figures, "cross_lake_encounter_rates_high_conf.png"),
       p_enc_hc_cross, width = 12, height = 7, dpi = 300)


# =================================================================-
# 6. Q2 & Q3: CROSS-LAKE DIEL ENCOUNTER PATTERNS — GAMMs ####
# =================================================================-
# Individual random smooth allows each individual to have its own
# diel encounter shape. Lake is included as a parametric factor
# (fixed effect for the mean encounter level) and as a random smooth
# via s(hour, Lake, bs="fs") to allow lake-specific diel shapes,
# while s(hour, by=treatment) tests the population-level treatment
# difference in diel patterns.

message("\n=== Q2/Q3: CROSS-LAKE DIEL ENCOUNTER PATTERN GAMMs ===")

## 6.1 Prepare hourly datasets ----
prepare_hourly_gamm <- function(hourly_df) {
  hourly_df %>%
    mutate(
      treatment     = factor(treatment, levels = c("Control", "Mix")),
      individual_ID = factor(individual_ID),
      Lake          = factor(Lake),
      hour          = as.integer(hour),
      exp_day_z     = as.numeric(scale(exp_day)),
      weight_z      = as.numeric(scale(weight))
    ) %>%
    filter(!is.na(exp_day)) %>%
    arrange(Lake, individual_ID, exp_day, hour)
}

roach_h <- prepare_hourly_gamm(roach_hourly_all)
perch_h <- prepare_hourly_gamm(perch_hourly_all)

## 6.2 Helper: add AR1 start flag ----
# Reset AR1 at the start of each individual's time series within each lake.
add_ar_start <- function(data) {
  data %>%
    arrange(Lake, individual_ID, exp_day, hour) %>%
    mutate(start = is.na(lag(individual_ID)) |
             individual_ID != lag(individual_ID) |
             Lake != lag(Lake))
}

## 6.3 Helper: estimate rho ----
estimate_rho <- function(data, mod) {
  data$residuals <- resid(mod, type = "pearson")
  data %>%
    group_split(Lake, individual_ID) %>%
    map_df(~{
      acf_res <- acf(.x$residuals, plot = FALSE)
      tibble(Lake          = unique(.x$Lake),
             individual_ID = unique(.x$individual_ID),
             lag1_value    = acf_res$acf[2])
    }) %>%
    summarise(mean_rho = mean(lag1_value, na.rm = TRUE)) %>%
    pull(mean_rho)
}

## 6.4 Shared plot elements ----
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

plot_diel_enc_cross <- function(pred_summary, species) {
  p <- ggplot(pred_summary,
              aes(x = hour, y = encounters, colour = treatment, group = treatment)) +
    geom_ribbon(aes(x = hour, ymin = lower, ymax = upper, group = treatment),
                fill = "black", alpha = 0.08, inherit.aes = FALSE) +
    geom_line(linewidth = 1.5) +
    scale_colour_manual(values = enc_colours, name = "Treatment") +
    facet_wrap(~Lake, scales = "free_y") +
    enc_x_scale +
    labs(title = species,
         y     = "Predicted encounters\nwith pike per hour\n",
         x     = "\nTime of day") +
    theme_few()
  add_diel_lines(p)
}

plot_enc_contrast_cross <- function(contrasts_df, species) {
  p <- ggplot(contrasts_df,
              aes(x = hour, y = estimate, ymin = conf.low, ymax = conf.high)) +
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
# ROACH GAMM — CROSS-LAKE
# ----------------------------------------------------------------

message("\n--- Q2/Q3: Roach cross-lake diel encounter GAMM ---")

## 6.5 Roach — first-pass model ----
roach_h <- roach_h %>% add_ar_start()

roach_gamm_init <- bam(
  n_encounters_high_conf ~
    treatment + Lake + exp_day_z + weight_z +
    s(hour, by = treatment, k = 15, bs = "cc") +
    s(hour, Lake, bs = "fs", xt = list(bs = "cc")) +
    s(hour, individual_ID, bs = "fs", xt = list(bs = "cc")),
  data     = roach_h,
  method   = "fREML",
  discrete = TRUE,
  knots    = list(hour = c(0, 24)),
  family   = nb(link = "log")
)

roach_rho <- estimate_rho(roach_h, roach_gamm_init)
message(sprintf("Roach lag-1 autocorrelation (rho): %.2f", roach_rho))

## 6.6 Roach — AR-corrected GAMM ----
# Model structure:
#   treatment            — fixed effect: overall treatment difference (Q1/Q3)
#   Lake                 — fixed effect: mean encounter level per lake
#   exp_day_z            — fixed effect: linear trend over experiment
#   weight_z             — fixed effect: body size covariate
#   s(hour, by=treatment)— separate diel curves per treatment (Q2/Q3)
#   s(hour, Lake, bs=fs) — lake-specific deviation in diel shape (random smooth)
#   s(hour, individual_ID, bs=fs) — individual-level diel variation
roach_mod <- bam(
  n_encounters_high_conf ~
    treatment + Lake + exp_day_z + weight_z +
    s(hour, by = treatment, k = 15, bs = "cc") +
    s(hour, Lake, bs = "fs", xt = list(bs = "cc")) +
    s(hour, individual_ID, bs = "fs", xt = list(bs = "cc")),
  data     = roach_h,
  method   = "fREML",
  rho      = round(roach_rho, 1),
  AR.start = roach_h$start,
  discrete = TRUE,
  knots    = list(hour = c(0, 24)),
  family   = nb(link = "log")
)

## 6.7 Roach — diagnostics ----
check_resid(roach_mod, ask = FALSE)
print(appraise(roach_mod))
summary(roach_mod)
anova.gam(roach_mod)

## 6.8 Roach — predictions ----
newdata_roach <- tidyr::expand(
  roach_h,
  nesting(Lake, individual_ID, treatment),
  hour = seq(0, 23, length.out = 96)
) %>%
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
  group_by(Lake, treatment, hour) %>%
  summarise(encounters = mean(fitted), upper = mean(upper),
            lower = mean(lower), .groups = "drop")

## 6.9 Roach — treatment contrasts (pooled across lakes) ----
newdata_roach_pool <- newdata_roach %>% mutate(exp_day_z = 0, weight_z = 0)

roach_comp <- comparisons(roach_mod, newdata = newdata_roach_pool,
                          variables = "treatment", type = "response",
                          vcov = TRUE, by = "hour")

message("\nRoach — peak/trough cross-lake treatment difference:")
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

## 6.9b Roach — diel period peak analysis (cross-lake) ----
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
  mutate(fitted = ilink_roach(fit),
         upper  = ilink_roach(fit + 1.96 * se.fit),
         lower  = ilink_roach(fit - 1.96 * se.fit))

roach_period_summary <- roach_period_raw %>%
  group_by(diel_period) %>%
  summarise(mean_encounters = mean(fitted),
            lower_95ci      = mean(lower),
            upper_95ci      = mean(upper), .groups = "drop") %>%
  arrange(desc(mean_encounters))

message("\nRoach — mean predicted encounters by diel period (cross-lake, pooled):")
print(roach_period_summary)

roach_period_preds <- avg_predictions(
  roach_mod, newdata = newdata_roach_periods,
  by = "diel_period", type = "response"
)

roach_period_pairs <- hypotheses(
  roach_period_preds,
  hypothesis = c("b1-b2=0", "b1-b3=0", "b1-b4=0",
                 "b2-b3=0", "b2-b4=0", "b3-b4=0")
)

message("\nRoach — pairwise diel period contrasts (cross-lake):")
print(roach_period_pairs)

## 6.10 Roach — figures ----
# Diel curves faceted by lake + pooled contrast panel
roach_diel_plot     <- plot_diel_enc_cross(roach_pred_summary, "Roach")

# Pooled-across-lakes contrast
roach_comp_pooled <- roach_comp %>%
  mutate(hour_int = round(hour, 2)) %>%
  group_by(hour_int) %>%
  summarise(estimate  = mean(estimate),
            conf.low  = mean(conf.low),
            conf.high = mean(conf.high), .groups = "drop") %>%
  rename(hour = hour_int)

roach_contrast_plot <- plot_enc_contrast_cross(roach_comp_pooled, "Roach (pooled across lakes)")

roach_figure <- ggarrange(roach_diel_plot, roach_contrast_plot,
                          labels = c("A", "B"), ncol = 2, widths = c(2, 1))

ggsave(file.path(paths$figures, "cross_lake_roach_diel_encounter_figure.png"),
       roach_figure, width = 16, height = 6, dpi = 300)
message("Saved: cross_lake_roach_diel_encounter_figure.png")

# ----------------------------------------------------------------
# PERCH GAMM — CROSS-LAKE
# ----------------------------------------------------------------

message("\n--- Q2/Q3: Perch cross-lake diel encounter GAMM ---")

## 6.11 Perch — first-pass model ----
perch_h <- perch_h %>% add_ar_start()

perch_gamm_init <- bam(
  n_encounters_high_conf ~
    treatment + Lake + exp_day_z + weight_z +
    s(hour, by = treatment, k = 15, bs = "cc") +
    s(hour, Lake, bs = "fs", xt = list(bs = "cc")) +
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
    treatment + Lake + exp_day_z + weight_z +
    s(hour, by = treatment, k = 15, bs = "cc") +
    s(hour, Lake, bs = "fs", xt = list(bs = "cc")) +
    s(hour, individual_ID, bs = "fs", xt = list(bs = "cc")),
  data     = perch_h,
  method   = "fREML",
  rho      = round(perch_rho, 1),
  AR.start = perch_h$start,
  discrete = TRUE,
  knots    = list(hour = c(0, 24)),
  family   = nb(link = "log")
)

## 6.13 Perch — diagnostics ----
check_resid(perch_mod, ask = FALSE)
print(appraise(perch_mod))
summary(perch_mod)
anova.gam(perch_mod)

## 6.14 Perch — predictions ----
newdata_perch <- tidyr::expand(
  perch_h,
  nesting(Lake, individual_ID, treatment),
  hour = seq(0, 23, length.out = 96)
) %>%
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
  group_by(Lake, treatment, hour) %>%
  summarise(encounters = mean(fitted), upper = mean(upper),
            lower = mean(lower), .groups = "drop")

## 6.15 Perch — treatment contrasts ----
perch_comp <- comparisons(perch_mod, newdata = newdata_perch,
                          variables = "treatment", type = "response",
                          vcov = TRUE, by = "hour")

message("\nPerch — peak/trough cross-lake treatment difference:")
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

## 6.15b Perch — diel period peak analysis (cross-lake) ----
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
  mutate(fitted = ilink_perch(fit),
         upper  = ilink_perch(fit + 1.96 * se.fit),
         lower  = ilink_perch(fit - 1.96 * se.fit))

perch_period_summary <- perch_period_raw %>%
  group_by(diel_period) %>%
  summarise(mean_encounters = mean(fitted),
            lower_95ci      = mean(lower),
            upper_95ci      = mean(upper), .groups = "drop") %>%
  arrange(desc(mean_encounters))

message("\nPerch — mean predicted encounters by diel period (cross-lake, pooled):")
print(perch_period_summary)

perch_period_preds <- avg_predictions(
  perch_mod, newdata = newdata_perch_periods,
  by = "diel_period", type = "response"
)

perch_period_pairs <- hypotheses(
  perch_period_preds,
  hypothesis = c("b1-b2=0", "b1-b3=0", "b1-b4=0",
                 "b2-b3=0", "b2-b4=0", "b3-b4=0")
)

message("\nPerch — pairwise diel period contrasts (cross-lake):")
print(perch_period_pairs)

## 6.16 Perch — figures ----
perch_diel_plot <- plot_diel_enc_cross(perch_pred_summary, "Perch")

perch_comp_pooled <- perch_comp %>%
  mutate(hour_int = round(hour, 2)) %>%
  group_by(hour_int) %>%
  summarise(estimate  = mean(estimate),
            conf.low  = mean(conf.low),
            conf.high = mean(conf.high), .groups = "drop") %>%
  rename(hour = hour_int)

perch_contrast_plot <- plot_enc_contrast_cross(perch_comp_pooled, "Perch (pooled across lakes)")

perch_figure <- ggarrange(perch_diel_plot, perch_contrast_plot,
                          labels = c("A", "B"), ncol = 2, widths = c(2, 1))

ggsave(file.path(paths$figures, "cross_lake_perch_diel_encounter_figure.png"),
       perch_figure, width = 16, height = 6, dpi = 300)
message("Saved: cross_lake_perch_diel_encounter_figure.png")

## 6.17 Combined cross-lake diel overview ----
enc_diel_summary <- bind_rows(
  roach_pred_summary %>% mutate(Species = "Roach"),
  perch_pred_summary %>% mutate(Species = "Perch")
)

p_diel_cross <- ggplot(enc_diel_summary,
                       aes(x = hour, y = encounters,
                           colour = treatment, group = treatment)) +
  geom_ribbon(aes(x = hour, ymin = lower, ymax = upper, group = treatment),
              fill = "black", alpha = 0.08, inherit.aes = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = enc_colours, name = "Treatment") +
  facet_grid(Species ~ Lake, scales = "free_y") +
  enc_x_scale +
  geom_vline(xintercept = nauticalDawn, linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = sunrise,      linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = sunset,       linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = nauticalDusk, linetype = "dashed", alpha = 0.4) +
  labs(
    subtitle = "Dashed lines: nautical dawn / sunrise / sunset / nautical dusk",
    y        = "Predicted encounters with pike per hour\n",
    x        = "\nTime of day"
  ) +
  theme_few()

ggsave(file.path(paths$figures, "cross_lake_combined_diel_overview.png"),
       p_diel_cross, width = 16, height = 8, dpi = 300)
message("Saved: cross_lake_combined_diel_overview.png")


# =================================================================-
# 7. SAVE OUTPUTS ####
# =================================================================-

saveRDS(enc_mod_roach_hc,   file.path(paths$tables, "glmm_roach_high_conf.rds"))
saveRDS(enc_mod_perch_hc,   file.path(paths$tables, "glmm_perch_high_conf.rds"))
saveRDS(enc_mod_roach_comb, file.path(paths$tables, "glmm_roach_combined.rds"))
saveRDS(enc_mod_perch_comb, file.path(paths$tables, "glmm_perch_combined.rds"))
saveRDS(roach_mod,          file.path(paths$tables, "gamm_roach_encounter_diel_crosslake.rds"))
saveRDS(perch_mod,          file.path(paths$tables, "gamm_perch_encounter_diel_crosslake.rds"))

write.csv(roach_pred_summary,
          file.path(paths$tables, "roach_diel_encounter_predictions_crosslake.csv"),
          row.names = FALSE)
write.csv(perch_pred_summary,
          file.path(paths$tables, "perch_diel_encounter_predictions_crosslake.csv"),
          row.names = FALSE)
write.csv(enc_rate_summary,
          file.path(paths$tables, "cross_lake_encounter_rate_summary.csv"),
          row.names = FALSE)
write.csv(bind_rows(enc_roach_total, enc_perch_total),
          file.path(paths$tables, "individual_total_encounters_all_lakes.csv"),
          row.names = FALSE)

message("\n=== CROSS-LAKE ANALYSIS COMPLETE ===")
message(sprintf("Figures saved to: %s", paths$figures))
message(sprintf("Tables saved to:  %s", paths$tables))


# =================================================================-
# 8. SESSION INFO ####
# =================================================================-

sink(file.path(paths$tables, "session_info_cross_lake_encounter_analysis.txt"))
cat("Cross-Lake Encounter Rate Analysis — Session Info\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")
print(sessionInfo())
sink()
