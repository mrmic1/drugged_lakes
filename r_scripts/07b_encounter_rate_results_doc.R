# =================================================================-
# Encounter Rate Analysis — Results Word Document ####
# =================================================================-
#
# Reads all saved model objects and data frames produced by:
#   07_muddyfoot_encounter_rate_analysis.R
#   07_BT_encounter_rate_analysis.R
#   07_cow_encounter_rate_analysis.R
#   08_cross_lake_encounter_rate_analysis.R
#
# Compiles results into a single .docx file using officer + flextable,
# organised as:
#   1. Muddyfoot
#   2. BT
#   3. Cow Paradise
#   4. Combined (all lakes)
#
# Within each section:
#   (A) Pre-filtered Q1: encounter rate by treatment (high-conf only)
#   (B) Post-filtered Q1: encounter rate by treatment (HC + combined)
#   (C) Post-filtered Q2/Q3: diel pattern GAMM
#
# =================================================================-

# =================================================================-
# 1. SETUP ####
# =================================================================-

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(officer)
  library(flextable)
  library(emmeans)
  library(glmmTMB)
  library(mgcv)
})

# Output path
out_dir  <- "./tables/"
out_file <- file.path(out_dir, "encounter_rate_analysis_results.docx")

# =================================================================-
# 2. TABLE PATHS ####
# =================================================================-

paths <- list(
  mf  = "./tables/muddyfoot/encounter_analysis/",
  bt  = "./tables/BT/encounter_analysis/",
  cow = "./tables/cow_paradise/encounter_analysis/",
  cl  = "./tables/cross_lake/encounter_analysis/"
)

# =================================================================-
# 3. HELPER FUNCTIONS ####
# =================================================================-

## 3.1 GLMM fixed-effects table ----
# Returns a data frame of fixed effects with IRR and p-value stars.
glmm_coef_df <- function(model) {
  s     <- summary(model)
  coefs <- as.data.frame(s$coefficients$cond)
  data.frame(
    Parameter = rownames(coefs),
    Estimate  = round(coefs[["Estimate"]],   3),
    SE        = round(coefs[["Std. Error"]], 3),
    z         = round(coefs[["z value"]],    2),
    p         = round(coefs[["Pr(>|z|)"]],   4),
    IRR       = round(exp(coefs[["Estimate"]]), 3),
    stringsAsFactors = FALSE
  ) %>%
    mutate(Sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE       ~ ""
    ))
}

## 3.2 GLMM fit-statistics one-liner ----
glmm_fit_str <- function(model) {
  sprintf("AIC = %.1f | Log-Lik = %.1f | theta (dispersion) = %.3f",
          AIC(model), as.numeric(logLik(model)), sigma(model))
}

## 3.3 Emmeans table (treatment contrasts, response scale) ----
# data must be passed explicitly so emmeans can reconstruct the reference
# grid when the model was fitted in a different session and loaded from .rds.
# CI column names vary by model type (asymp.LCL vs lower.CL), so we detect
# them dynamically rather than hardcoding.

emm_df <- function(model, label = "", data = NULL) {
  emm <- emmeans(model, specs = ~ treatment, type = "response",
                 offset = log(1), data = data)
  df  <- as.data.frame(emm)

  # Detect response and CI column names
  resp_col <- intersect(c("response", "emmean", "rate"), names(df))[1]
  lcl_col  <- intersect(c("asymp.LCL", "lower.CL", "lwr"), names(df))[1]
  ucl_col  <- intersect(c("asymp.UCL", "upper.CL", "upr"), names(df))[1]

  out <- df %>%
    mutate(Model = label) %>%
    rename(Treatment   = treatment,
           Rate_per_day = !!resp_col)

  if (!is.na(lcl_col))
    out <- out %>% rename(Lower_95CI = !!lcl_col, Upper_95CI = !!ucl_col)

  out %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    select(any_of(c("Model", "Treatment", "Rate_per_day", "SE",
                    "Lower_95CI", "Upper_95CI")))
}

emm_contrast_df <- function(model, label = "", data = NULL) {
  emm  <- emmeans(model, specs = ~ treatment, type = "response",
                  offset = log(1), data = data)
  # confint() adds CI columns; as.data.frame(pairs()) gives p-values
  cont_ci <- as.data.frame(confint(pairs(emm)))
  cont_p  <- as.data.frame(pairs(emm)) %>% select(contrast, p.value)
  cont    <- left_join(cont_ci, cont_p, by = "contrast")

  # Detect ratio and CI column names
  ratio_col <- intersect(c("ratio", "odds.ratio"), names(cont))[1]
  lcl_col   <- intersect(c("asymp.LCL", "lower.CL", "lwr"), names(cont))[1]
  ucl_col   <- intersect(c("asymp.UCL", "upper.CL", "upr"), names(cont))[1]

  out <- cont %>%
    mutate(Model = label) %>%
    rename(Contrast = contrast)

  if (!is.na(ratio_col))
    out <- out %>% rename(Ratio = !!ratio_col)
  if (!is.na(lcl_col))
    out <- out %>% rename(Lower_95CI = !!lcl_col, Upper_95CI = !!ucl_col)

  out %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    select(any_of(c("Model", "Contrast", "Ratio", "SE",
                    "Lower_95CI", "Upper_95CI", "p.value")))
}

## 3.4 GAMM parametric-coefficients table ----
gamm_param_df <- function(model) {
  s <- summary(model)
  as.data.frame(s$p.table) %>%
    tibble::rownames_to_column("Parameter") %>%
    rename(
      Estimate = Estimate,
      SE       = `Std. Error`,
      t        = `t value`,
      p        = `Pr(>|t|)`
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    mutate(Sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE       ~ ""
    ))
}

## 3.5 GAMM smooth-terms table ----
gamm_smooth_df <- function(model) {
  s <- summary(model)
  as.data.frame(s$s.table) %>%
    tibble::rownames_to_column("Smooth") %>%
    rename(edf = edf, Ref.df = Ref.df, F = F, p = `p-value`) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    mutate(Sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE       ~ ""
    ))
}

## 3.6 GAMM fit string ----
gamm_fit_str <- function(model) {
  s <- summary(model)
  sprintf("R²(adj) = %.3f | Dev. explained = %.1f%% | AIC = %.1f | n = %d",
          s$r.sq, s$dev.expl * 100, AIC(model), s$n)
}

## 3.7 flextable theme ----
ft_style <- function(ft, caption = NULL) {
  ft <- ft %>%
    theme_booktabs() %>%
    autofit() %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    align(j = 1, align = "left",   part = "body") %>%
    align(j = -1, align = "center", part = "body")
  ft
}

## 3.8 Officer heading helpers ----
h1 <- function(doc, text) body_add_par(doc, text, style = "heading 1")
h2 <- function(doc, text) body_add_par(doc, text, style = "heading 2")
h3 <- function(doc, text) body_add_par(doc, text, style = "heading 3")
cap <- function(doc, text) body_add_par(doc, text, style = "Normal")
gap <- function(doc)       body_add_par(doc, "", style = "Normal")

add_ft <- function(doc, ft) {
  doc <- body_add_flextable(doc, ft)
  gap(doc)
}

## 3.9 Build one lake's Q1 block (pre or post) ----
add_q1_block <- function(doc, roach_mod, perch_mod,
                          summary_csv, ind_csv,
                          lake, filter_label, hc_label) {

  h3(doc, sprintf("%s — %s | Q1 Encounter Rate (Treatment Effect)", lake, filter_label))

  # Load individual-level data to pass to emmeans (needed when model was
  # loaded from .rds and original data frame is no longer in environment)
  ind_data <- tryCatch(
    read.csv(ind_csv) %>%
      mutate(
        treatment     = factor(treatment, levels = c("Control", "Mix")),
        individual_ID = factor(individual_ID)
      ),
    error = function(e) NULL
  )

  roach_data <- if (!is.null(ind_data)) filter(ind_data, Species == "Roach") else NULL
  perch_data <- if (!is.null(ind_data)) filter(ind_data, Species == "Perch") else NULL

  # --- Summary statistics ---
  doc <- gap(doc)
  doc <- body_add_par(doc, "Summary statistics", style = "Normal")

  summ <- tryCatch(read.csv(summary_csv), error = function(e) NULL)
  if (!is.null(summ)) {
    ft <- summ %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      flextable() %>%
      ft_style()
    doc <- body_add_flextable(doc, ft)
    doc <- cap(doc, sprintf(
      "Table: Mean (SD) and median (IQR) high-confidence encounter rate (encounters/day) by species and treatment group. %s, %s.",
      lake, filter_label))
  }

  doc <- gap(doc)

  # --- Roach GLMM ---
  doc <- body_add_par(doc, sprintf("Roach — Negative Binomial GLMM (%s)", hc_label),
                      style = "Normal")
  doc <- body_add_par(doc, glmm_fit_str(roach_mod), style = "Normal")

  ft_coef <- glmm_coef_df(roach_mod) %>% flextable() %>% ft_style()
  doc <- body_add_flextable(doc, ft_coef)
  doc <- cap(doc, sprintf(
    "Table: Fixed-effect coefficients from negative binomial GLMM. Estimate and SE are on the log scale; IRR = exp(Estimate) is the incidence rate ratio relative to the reference level (Control). %s, %s, Roach.",
    lake, filter_label))
  doc <- gap(doc)

  ft_emm <- emm_df(roach_mod,
                   label = sprintf("%s %s Roach", lake, filter_label),
                   data  = roach_data) %>%
    flextable() %>% ft_style()
  doc <- body_add_flextable(doc, ft_emm)
  doc <- cap(doc, sprintf(
    "Table: Estimated marginal means (encounters/day, response scale) by treatment group. %s, %s, Roach.",
    lake, filter_label))
  doc <- gap(doc)

  ft_cont <- emm_contrast_df(roach_mod,
                              label = sprintf("%s %s Roach", lake, filter_label),
                              data  = roach_data) %>%
    flextable() %>% ft_style()
  doc <- body_add_flextable(doc, ft_cont)
  doc <- cap(doc, sprintf(
    "Table: Pairwise treatment contrast (ratio of marginal means; Mix vs Control). %s, %s, Roach.",
    lake, filter_label))
  doc <- gap(doc)

  # --- Perch GLMM ---
  doc <- body_add_par(doc, sprintf("Perch — Negative Binomial GLMM (%s)", hc_label),
                      style = "Normal")
  doc <- body_add_par(doc, glmm_fit_str(perch_mod), style = "Normal")

  ft_coef_p <- glmm_coef_df(perch_mod) %>% flextable() %>% ft_style()
  doc <- body_add_flextable(doc, ft_coef_p)
  doc <- cap(doc, sprintf(
    "Table: Fixed-effect coefficients from negative binomial GLMM. %s, %s, Perch.",
    lake, filter_label))
  doc <- gap(doc)

  ft_emm_p <- emm_df(perch_mod,
                     label = sprintf("%s %s Perch", lake, filter_label),
                     data  = perch_data) %>%
    flextable() %>% ft_style()
  doc <- body_add_flextable(doc, ft_emm_p)
  doc <- cap(doc, sprintf(
    "Table: Estimated marginal means (encounters/day) by treatment group. %s, %s, Perch.",
    lake, filter_label))
  doc <- gap(doc)

  ft_cont_p <- emm_contrast_df(perch_mod,
                                label = sprintf("%s %s Perch", lake, filter_label),
                                data  = perch_data) %>%
    flextable() %>% ft_style()
  doc <- body_add_flextable(doc, ft_cont_p)
  doc <- cap(doc, sprintf(
    "Table: Pairwise treatment contrast. %s, %s, Perch.",
    lake, filter_label))
  doc <- gap(doc)

  doc
}

## 3.10 Build one lake's GAMM block ----
add_gamm_block <- function(doc, roach_gamm, perch_gamm,
                            roach_pred_csv, perch_pred_csv,
                            lake) {

  h3(doc, sprintf("%s — Q2 & Q3 Diel Encounter Pattern (GAMM)", lake))

  for (sp in c("Roach", "Perch")) {
    mod      <- if (sp == "Roach") roach_gamm else perch_gamm
    pred_csv <- if (sp == "Roach") roach_pred_csv else perch_pred_csv

    doc <- body_add_par(doc, sprintf("%s — AR-corrected GAMM", sp), style = "Normal")
    doc <- body_add_par(doc, gamm_fit_str(mod), style = "Normal")
    doc <- gap(doc)

    # Parametric coefficients
    ft_param <- gamm_param_df(mod) %>% flextable() %>% ft_style()
    doc <- body_add_flextable(doc, ft_param)
    doc <- cap(doc, sprintf(
      "Table: Parametric coefficients from the AR-corrected GAMM. Estimates are on the log scale. %s, %s.",
      lake, sp))
    doc <- gap(doc)

    # Smooth terms
    ft_smooth <- gamm_smooth_df(mod) %>% flextable() %>% ft_style()
    doc <- body_add_flextable(doc, ft_smooth)
    doc <- cap(doc, sprintf(
      "Table: Approximate significance of smooth terms. edf = estimated degrees of freedom. %s, %s.",
      lake, sp))
    doc <- gap(doc)

    # Diel predictions summary (averaged across individuals)
    pred <- tryCatch(read.csv(pred_csv), error = function(e) NULL)
    if (!is.null(pred)) {
      pred_summ <- pred %>%
        group_by(treatment) %>%
        summarise(
          peak_hour     = hour[which.max(encounters)],
          peak_enc_rate = round(max(encounters), 4),
          mean_enc_rate = round(mean(encounters), 4),
          .groups = "drop"
        )
      ft_pred <- pred_summ %>% flextable() %>% ft_style()
      doc <- body_add_flextable(doc, ft_pred)
      doc <- cap(doc, sprintf(
        "Table: Summary of predicted diel encounter curve — peak hour and peak/mean predicted encounter rate per treatment group (model predictions averaged across individuals, at mean covariate values). %s, %s.",
        lake, sp))
      doc <- gap(doc)
    }
  }

  doc
}


# =================================================================-
# 4. BUILD DOCUMENT ####
# =================================================================-

doc <- read_docx()

# Title page
doc <- body_add_par(doc,
  "Predator–Prey Encounter Rate Analysis: Model Results",
  style = "heading 1")
doc <- body_add_par(doc,
  "Muddyfoot | BT | Cow Paradise | Combined (all lakes)",
  style = "heading 2")
doc <- body_add_par(doc,
  paste("Generated:", format(Sys.time(), "%d %B %Y %H:%M")),
  style = "Normal")
doc <- gap(doc)
doc <- body_add_par(doc,
  paste(
    "This document summarises results from the encounter rate analyses",
    "conducted for each lake separately and across all lakes combined.",
    "Results are presented in three subsections per lake:",
    "(A) Pre-filtered Q1 models (encounter rates before removal of",
    "post-predation tracks — rates for predated individuals will be inflated);",
    "(B) Post-filtered Q1 models (encounter rates after removal of",
    "post-predation tracks); and",
    "(C) Q2/Q3 diel pattern GAMMs (post-filtered data only).",
    "BT and Muddyfoot use high-confidence = ≤1.4m, combined = ≤2.8m.",
    "Cow Paradise uses high-confidence = ≤2.0m, combined = ≤4.0m, to account",
    "for its higher GPS error (σ = 0.827m vs 0.379–0.500m at BT/Muddyfoot)."
  ),
  style = "Normal")
doc <- gap(doc)

# =================================================================-
# SECTION 1: MUDDYFOOT
# =================================================================-

mf_path <- paths$mf

doc <- h1(doc, "1. Muddyfoot")

# ---- 1A: Pre-filter ----
doc <- h2(doc, "1A. Pre-Filtered Encounter Rate Analysis (Q1)")
doc <- body_add_par(doc,
  paste("Models fitted to the full (unfiltered) dataset.",
        "Encounter rates for predated individuals are expected to be inflated.",
        "Only high-confidence encounters (≤1.4m) are modelled here."),
  style = "Normal")
doc <- gap(doc)

mf_roach_pre <- readRDS(file.path(mf_path, "glmm_roach_hc_pre.rds"))
mf_perch_pre <- readRDS(file.path(mf_path, "glmm_perch_hc_pre.rds"))

doc <- add_q1_block(
  doc,
  roach_mod    = mf_roach_pre,
  perch_mod    = mf_perch_pre,
  summary_csv  = file.path(mf_path, "encounter_rate_summary_pre_filter.csv"),
  ind_csv      = file.path(mf_path, "individual_total_encounters_pre_filter.csv"),
  lake         = "Muddyfoot",
  filter_label = "Pre-filter",
  hc_label     = "≤2.0m"
)

# ---- 1B: Post-filter ----
doc <- h2(doc, "1B. Post-Filtered Encounter Rate Analysis (Q1)")
doc <- body_add_par(doc,
  paste("Models fitted after removing tracking records following confirmed",
        "or suspected predation/mortality events."),
  style = "Normal")
doc <- gap(doc)

mf_roach_hc   <- readRDS(file.path(mf_path, "glmm_roach_high_conf.rds"))
mf_perch_hc   <- readRDS(file.path(mf_path, "glmm_perch_high_conf.rds"))
mf_roach_comb <- readRDS(file.path(mf_path, "glmm_roach_combined.rds"))
mf_perch_comb <- readRDS(file.path(mf_path, "glmm_perch_combined.rds"))

doc <- add_q1_block(
  doc,
  roach_mod    = mf_roach_hc,
  perch_mod    = mf_perch_hc,
  summary_csv  = file.path(mf_path, "individual_total_encounters_with_treatment.csv"),
  ind_csv      = file.path(mf_path, "individual_total_encounters_with_treatment.csv"),
  lake         = "Muddyfoot",
  filter_label = "Post-filter (high-confidence)",
  hc_label     = "≤2.0m"
)

doc <- h3(doc, "Muddyfoot — Post-Filter | Combined Encounters (≤2.8m)")

mf_post_data <- read.csv(file.path(mf_path, "individual_total_encounters_with_treatment.csv")) %>%
  mutate(treatment = factor(treatment, levels = c("Control", "Mix")),
         individual_ID = factor(individual_ID))
mf_post_roach <- filter(mf_post_data, Species == "Roach")
mf_post_perch <- filter(mf_post_data, Species == "Perch")

for (sp_label in c("Roach", "Perch")) {
  mod      <- if (sp_label == "Roach") mf_roach_comb else mf_perch_comb
  sp_data  <- if (sp_label == "Roach") mf_post_roach else mf_post_perch
  lbl      <- sprintf("Muddyfoot Post-filter %s Combined", sp_label)
  thresh   <- "≤2.8m"
  doc <- body_add_par(doc, sprintf("%s — Combined GLMM (%s)", sp_label, thresh), style = "Normal")
  doc <- body_add_par(doc, glmm_fit_str(mod), style = "Normal")
  doc <- body_add_flextable(doc, glmm_coef_df(mod) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Fixed-effect coefficients, combined encounter GLMM (%s). Muddyfoot, Post-filter, %s.", thresh, sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_df(mod, lbl, data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Estimated marginal means, combined encounters. Muddyfoot, Post-filter, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_contrast_df(mod, lbl, data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Pairwise treatment contrast, combined encounters. Muddyfoot, Post-filter, %s.", sp_label))
  doc <- gap(doc)
}

# ---- 1C: GAMM ----
doc <- h2(doc, "1C. Diel Encounter Pattern — GAMMs (Q2 & Q3)")

mf_roach_gamm <- readRDS(file.path(mf_path, "gamm_roach_encounter_diel.rds"))
mf_perch_gamm <- readRDS(file.path(mf_path, "gamm_perch_encounter_diel.rds"))

doc <- add_gamm_block(
  doc,
  roach_gamm     = mf_roach_gamm,
  perch_gamm     = mf_perch_gamm,
  roach_pred_csv = file.path(mf_path, "roach_diel_encounter_predictions.csv"),
  perch_pred_csv = file.path(mf_path, "perch_diel_encounter_predictions.csv"),
  lake           = "Muddyfoot"
)


# =================================================================-
# SECTION 2: BT
# =================================================================-

bt_path <- paths$bt

doc <- h1(doc, "2. BT")

# ---- 2A: Pre-filter ----
doc <- h2(doc, "2A. Pre-Filtered Encounter Rate Analysis (Q1)")
doc <- body_add_par(doc,
  "Models fitted to the full (unfiltered) dataset. High-confidence encounters (≤1.4m) only.",
  style = "Normal")
doc <- gap(doc)

bt_roach_pre <- readRDS(file.path(bt_path, "glmm_roach_hc_pre.rds"))
bt_perch_pre <- readRDS(file.path(bt_path, "glmm_perch_hc_pre.rds"))

doc <- add_q1_block(
  doc,
  roach_mod    = bt_roach_pre,
  perch_mod    = bt_perch_pre,
  summary_csv  = file.path(bt_path, "encounter_rate_summary_pre_filter.csv"),
  ind_csv      = file.path(bt_path, "individual_total_encounters_pre_filter.csv"),
  lake         = "BT",
  filter_label = "Pre-filter",
  hc_label     = "≤2.0m"
)

# ---- 2B: Post-filter ----
doc <- h2(doc, "2B. Post-Filtered Encounter Rate Analysis (Q1)")
doc <- gap(doc)

bt_roach_hc   <- readRDS(file.path(bt_path, "glmm_roach_high_conf.rds"))
bt_perch_hc   <- readRDS(file.path(bt_path, "glmm_perch_high_conf.rds"))
bt_roach_comb <- readRDS(file.path(bt_path, "glmm_roach_combined.rds"))
bt_perch_comb <- readRDS(file.path(bt_path, "glmm_perch_combined.rds"))

doc <- add_q1_block(
  doc,
  roach_mod    = bt_roach_hc,
  perch_mod    = bt_perch_hc,
  summary_csv  = file.path(bt_path, "individual_total_encounters_with_treatment.csv"),
  ind_csv      = file.path(bt_path, "individual_total_encounters_with_treatment.csv"),
  lake         = "BT",
  filter_label = "Post-filter (high-confidence)",
  hc_label     = "≤2.0m"
)

doc <- h3(doc, "BT — Post-Filter | Combined Encounters (≤2.8m)")

bt_post_data <- read.csv(file.path(bt_path, "individual_total_encounters_with_treatment.csv")) %>%
  mutate(treatment = factor(treatment, levels = c("Control", "Mix")),
         individual_ID = factor(individual_ID))

for (sp_label in c("Roach", "Perch")) {
  mod     <- if (sp_label == "Roach") bt_roach_comb else bt_perch_comb
  sp_data <- filter(bt_post_data, Species == sp_label)
  lbl     <- sprintf("BT Post-filter %s Combined", sp_label)
  doc <- body_add_par(doc, sprintf("%s — Combined GLMM (≤4.0m)", sp_label), style = "Normal")
  doc <- body_add_par(doc, glmm_fit_str(mod), style = "Normal")
  doc <- body_add_flextable(doc, glmm_coef_df(mod) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Fixed-effect coefficients, combined GLMM (≤2.8m). BT, Post-filter, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_df(mod, lbl, data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Estimated marginal means, combined encounters. BT, Post-filter, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_contrast_df(mod, lbl, data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Pairwise treatment contrast, combined encounters. BT, Post-filter, %s.", sp_label))
  doc <- gap(doc)
}

# ---- 2C: GAMM ----
doc <- h2(doc, "2C. Diel Encounter Pattern — GAMMs (Q2 & Q3)")

bt_roach_gamm <- readRDS(file.path(bt_path, "gamm_roach_encounter_diel.rds"))
bt_perch_gamm <- readRDS(file.path(bt_path, "gamm_perch_encounter_diel.rds"))

doc <- add_gamm_block(
  doc,
  roach_gamm     = bt_roach_gamm,
  perch_gamm     = bt_perch_gamm,
  roach_pred_csv = file.path(bt_path, "roach_diel_encounter_predictions.csv"),
  perch_pred_csv = file.path(bt_path, "perch_diel_encounter_predictions.csv"),
  lake           = "BT"
)


# =================================================================-
# SECTION 3: COW PARADISE
# =================================================================-

cow_path <- paths$cow

doc <- h1(doc, "3. Cow Paradise")
doc <- body_add_par(doc,
  paste("Note: GPS error at Cow Paradise (σ = 0.827m) is higher than at BT or Muddyfoot.",
        "Encounter thresholds are set accordingly:",
        "high-confidence = ≤2.0m, combined (high-confidence + probable) = ≤4.0m."),
  style = "Normal")
doc <- gap(doc)

# ---- 3A: Pre-filter ----
doc <- h2(doc, "3A. Pre-Filtered Encounter Rate Analysis (Q1)")
doc <- gap(doc)

cow_roach_pre <- readRDS(file.path(cow_path, "glmm_roach_hc_pre.rds"))
cow_perch_pre <- readRDS(file.path(cow_path, "glmm_perch_hc_pre.rds"))

doc <- add_q1_block(
  doc,
  roach_mod    = cow_roach_pre,
  perch_mod    = cow_perch_pre,
  summary_csv  = file.path(cow_path, "encounter_rate_summary_pre_filter.csv"),
  ind_csv      = file.path(cow_path, "individual_total_encounters_pre_filter.csv"),
  lake         = "Cow Paradise",
  filter_label = "Pre-filter",
  hc_label     = "≤2.0m"
)

# ---- 3B: Post-filter ----
doc <- h2(doc, "3B. Post-Filtered Encounter Rate Analysis (Q1)")
doc <- gap(doc)

cow_roach_hc   <- readRDS(file.path(cow_path, "glmm_roach_high_conf.rds"))
cow_perch_hc   <- readRDS(file.path(cow_path, "glmm_perch_high_conf.rds"))
cow_roach_comb <- readRDS(file.path(cow_path, "glmm_roach_combined.rds"))
cow_perch_comb <- readRDS(file.path(cow_path, "glmm_perch_combined.rds"))

doc <- add_q1_block(
  doc,
  roach_mod    = cow_roach_hc,
  perch_mod    = cow_perch_hc,
  summary_csv  = file.path(cow_path, "individual_total_encounters_with_treatment.csv"),
  ind_csv      = file.path(cow_path, "individual_total_encounters_with_treatment.csv"),
  lake         = "Cow Paradise",
  filter_label = "Post-filter (high-confidence)",
  hc_label     = "≤2.0m"
)

doc <- h3(doc, "Cow Paradise — Post-Filter | Combined Encounters (≤4.0m)")

cow_post_data <- read.csv(file.path(cow_path, "individual_total_encounters_with_treatment.csv")) %>%
  mutate(treatment = factor(treatment, levels = c("Control", "Mix")),
         individual_ID = factor(individual_ID))

for (sp_label in c("Roach", "Perch")) {
  mod     <- if (sp_label == "Roach") cow_roach_comb else cow_perch_comb
  sp_data <- filter(cow_post_data, Species == sp_label)
  lbl     <- sprintf("Cow Paradise Post-filter %s Combined", sp_label)
  doc <- body_add_par(doc, sprintf("%s — Combined GLMM (≤4.0m)", sp_label), style = "Normal")
  doc <- body_add_par(doc, glmm_fit_str(mod), style = "Normal")
  doc <- body_add_flextable(doc, glmm_coef_df(mod) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Fixed-effect coefficients, combined GLMM (≤6.0m). Cow Paradise, Post-filter, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_df(mod, lbl, data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Estimated marginal means, combined encounters. Cow Paradise, Post-filter, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_contrast_df(mod, lbl, data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Pairwise treatment contrast. Cow Paradise, Post-filter, %s.", sp_label))
  doc <- gap(doc)
}

# ---- 3C: GAMM ----
doc <- h2(doc, "3C. Diel Encounter Pattern — GAMMs (Q2 & Q3)")

cow_roach_gamm <- readRDS(file.path(cow_path, "gamm_roach_encounter_diel.rds"))
cow_perch_gamm <- readRDS(file.path(cow_path, "gamm_perch_encounter_diel.rds"))

doc <- add_gamm_block(
  doc,
  roach_gamm     = cow_roach_gamm,
  perch_gamm     = cow_perch_gamm,
  roach_pred_csv = file.path(cow_path, "roach_diel_encounter_predictions.csv"),
  perch_pred_csv = file.path(cow_path, "perch_diel_encounter_predictions.csv"),
  lake           = "Cow Paradise"
)


# =================================================================-
# SECTION 4: COMBINED (ALL LAKES)
# =================================================================-

cl_path <- paths$cl

doc <- h1(doc, "4. Combined Analysis — All Lakes")
doc <- body_add_par(doc,
  paste("Cross-lake models include lake as a fixed effect (parametric) and as a",
        "random smooth s(hour, Lake, bs='fs') in the GAMMs. Individual is also",
        "included as a random smooth. GLMMs use (1 | Lake / individual_ID).",
        "Cow Paradise encounter thresholds differ from BT/Muddyfoot",
        "(see lake-specific note above); the offset-based daily rate",
        "makes encounter counts comparable across lakes."),
  style = "Normal")
doc <- gap(doc)

# ---- 4A: Combined Q1 ----
doc <- h2(doc, "4A. Cross-Lake Encounter Rate Analysis (Q1)")

cl_roach_hc   <- readRDS(file.path(cl_path, "glmm_roach_high_conf.rds"))
cl_perch_hc   <- readRDS(file.path(cl_path, "glmm_perch_high_conf.rds"))
cl_roach_comb <- readRDS(file.path(cl_path, "glmm_roach_combined.rds"))
cl_perch_comb <- readRDS(file.path(cl_path, "glmm_perch_combined.rds"))

# Load individual-level cross-lake data for passing to emmeans
cl_ind_data <- read.csv(file.path(cl_path, "individual_total_encounters_all_lakes.csv")) %>%
  mutate(
    treatment     = factor(treatment, levels = c("Control", "Mix")),
    individual_ID = factor(individual_ID),
    Lake          = factor(Lake)
  )

# Cross-lake summary CSV (from 08 script)
cl_summ <- tryCatch(
  read.csv(file.path(cl_path, "cross_lake_encounter_rate_summary.csv")),
  error = function(e) NULL
)
if (!is.null(cl_summ)) {
  ft_cl <- cl_summ %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    flextable() %>% ft_style()
  doc <- body_add_flextable(doc, ft_cl)
  doc <- cap(doc,
    "Table: Mean (SD) and median (IQR) high-confidence encounter rate (encounters/day) by lake, species, and treatment group.")
  doc <- gap(doc)
}

for (sp_label in c("Roach", "Perch")) {
  hc_mod   <- if (sp_label == "Roach") cl_roach_hc   else cl_perch_hc
  comb_mod <- if (sp_label == "Roach") cl_roach_comb else cl_perch_comb
  sp_data  <- filter(cl_ind_data, Species == sp_label)

  doc <- h3(doc, sprintf("All Lakes — %s (high-confidence)", sp_label))
  doc <- body_add_par(doc, glmm_fit_str(hc_mod), style = "Normal")
  doc <- body_add_flextable(doc, glmm_coef_df(hc_mod) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Fixed-effect coefficients, cross-lake high-confidence GLMM. All lakes, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_df(hc_mod, sprintf("All Lakes %s HC", sp_label), data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Estimated marginal means (encounters/day, response scale), cross-lake. All lakes, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_contrast_df(hc_mod, sprintf("All Lakes %s HC", sp_label), data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Pairwise treatment contrast, cross-lake high-confidence GLMM. All lakes, %s.", sp_label))
  doc <- gap(doc)

  doc <- h3(doc, sprintf("All Lakes — %s (combined)", sp_label))
  doc <- body_add_par(doc, glmm_fit_str(comb_mod), style = "Normal")
  doc <- body_add_flextable(doc, glmm_coef_df(comb_mod) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Fixed-effect coefficients, cross-lake combined GLMM. All lakes, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_df(comb_mod, sprintf("All Lakes %s Comb", sp_label), data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Estimated marginal means, combined encounters. All lakes, %s.", sp_label))
  doc <- gap(doc)
  doc <- body_add_flextable(doc, emm_contrast_df(comb_mod, sprintf("All Lakes %s Comb", sp_label), data = sp_data) %>% flextable() %>% ft_style())
  doc <- cap(doc, sprintf("Table: Pairwise treatment contrast, combined encounters. All lakes, %s.", sp_label))
  doc <- gap(doc)
}

# ---- 4B: Cross-lake GAMM ----
doc <- h2(doc, "4B. Cross-Lake Diel Encounter Pattern — GAMMs (Q2 & Q3)")

cl_roach_gamm <- readRDS(file.path(cl_path, "gamm_roach_encounter_diel_crosslake.rds"))
cl_perch_gamm <- readRDS(file.path(cl_path, "gamm_perch_encounter_diel_crosslake.rds"))

doc <- add_gamm_block(
  doc,
  roach_gamm     = cl_roach_gamm,
  perch_gamm     = cl_perch_gamm,
  roach_pred_csv = file.path(cl_path, "roach_diel_encounter_predictions_crosslake.csv"),
  perch_pred_csv = file.path(cl_path, "perch_diel_encounter_predictions_crosslake.csv"),
  lake           = "All Lakes (combined)"
)


# =================================================================-
# 5. SAVE DOCUMENT ####
# =================================================================-

print(doc, target = out_file)
message(sprintf("\nDocument saved to: %s", out_file))
