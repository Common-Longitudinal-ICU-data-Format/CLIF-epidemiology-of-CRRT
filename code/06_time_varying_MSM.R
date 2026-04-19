####### Marginal Structural Model (MSM) For First 24h CRRT Dose ###############

# ================================ #
# ---- 0. SETUP ----
# ================================ #

## ---- A. Set up R environment ----

### ---- Clear existing environment ----
# Clear global environment
rm(list = ls())

# Close all open plotting devices
while (!is.null(dev.list())) dev.off()

# Clear console
cat("\014")

cat("Environment and plots cleared.\n")

# Set working directory to project root
# Set working directory to project root using config.json project_root
.find_config <- function() {
  candidates <- c("../config/config.json", "config/config.json")
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(sub("^--file=", "", file_arg))
    candidates <- c(file.path(script_dir, "..", "config", "config.json"), candidates)
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
    candidates <- c(file.path(script_dir, "..", "config", "config.json"), candidates)
  }
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  stop("Cannot find config/config.json.")
}
.config_path <- .find_config()
.config <- jsonlite::fromJSON(.config_path)
if (!is.null(.config$project_root) && nchar(.config$project_root) > 0) {
  setwd(.config$project_root)
} else {
  setwd(dirname(dirname(.config_path)))
}

cat("Working directory set to:", getwd(), "\n\n")

## ---- B. Load required packages, installing if needed ----

# Set CRAN mirror if not already configured (needed when --no-init-file is used)
if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

required_packages <- c("tidyverse", "readr", "arrow", "gtsummary", "cmprsk",
                       "survival", "jsonlite", "MatchIt", "WeightIt", "broom",
                       "cobalt", "EValue", "SuperLearner", "randomForest",
                       "xgboost","gam","survminer", "survey", "mice")
new_packages <- required_packages[!(required_packages %in% installed.packages()
                                    [,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, require, character.only = TRUE)

## ---- C. Create output directory if it doesn't exist ----
output_dir <- "output/final/msm/time_varying"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

## ---- D. Load configuration ----
config_path <- "config/config.json"
if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path)
}

config <- jsonlite::fromJSON(config_path)
SITE_NAME <- config$site_name

cat("Site:", SITE_NAME, "\n")
cat("Timezone:", config$timezone, "\n\n")

## ---- E. Load data ----
data_path <- "output/intermediate/msm_competing_risk_df.parquet"
if (!file.exists(data_path)) {
  stop("File '", data_path, "' not found.")
}
df <- arrow::read_parquet(data_path)

cat("Loaded data:", nrow(df), "rows\n")

# Define CCI component variable names (reused throughout script)
cci_vars <- c(
  "cci_myocardial_infarction", "cci_congestive_heart_failure",
  "cci_peripheral_vascular_disease", "cci_cerebrovascular_disease",
  "cci_dementia", "cci_chronic_pulmonary_disease",
  "cci_connective_tissue_disease", "cci_peptic_ulcer_disease",
  "cci_mild_liver_disease", "cci_diabetes_uncomplicated",
  "cci_diabetes_with_complications", "cci_hemiplegia",
  "cci_renal_disease", "cci_cancer",
  "cci_moderate_severe_liver_disease", "cci_metastatic_solid_tumor",
  "cci_aids"
)

# Pretty labels for CCI components (used in tables and balance plots)
cci_labels <- c(
  cci_myocardial_infarction         = "Myocardial Infarction",
  cci_congestive_heart_failure      = "Congestive Heart Failure",
  cci_peripheral_vascular_disease   = "Peripheral Vascular Disease",
  cci_cerebrovascular_disease       = "Cerebrovascular Disease",
  cci_dementia                      = "Dementia",
  cci_chronic_pulmonary_disease     = "Chronic Pulmonary Disease",
  cci_connective_tissue_disease     = "Connective Tissue Disease",
  cci_peptic_ulcer_disease          = "Peptic Ulcer Disease",
  cci_mild_liver_disease            = "Mild Liver Disease",
  cci_diabetes_uncomplicated        = "Diabetes Uncomplicated",
  cci_diabetes_with_complications   = "Diabetes w/ Complications",
  cci_hemiplegia                    = "Hemiplegia",
  cci_renal_disease                 = "Renal Disease",
  cci_cancer                        = "Cancer",
  cci_moderate_severe_liver_disease = "Moderate/Severe Liver Disease",
  cci_metastatic_solid_tumor        = "Metastatic Solid Tumor",
  cci_aids                          = "AIDS"
)

# Make sure outcome and covariates exist
required_vars <- c(
  "encounter_block",
  "time_to_event_90d", "outcome",
  "age_at_admission", "sex_category", "race_category", "ethnicity_category",
  "crrt_mode_category", "weight_kg",
  "crrt_dose_0_12", "crrt_dose_12_24",
  # Labs at t=0 and t=12
  "lactate_0", "bicarbonate_0", "potassium_0",
  "lactate_12", "bicarbonate_12", "potassium_12",
  # SOFA (kept in data for descriptive use; NOT used in models)
  "sofa_total_0", "sofa_total_12",
  # New time-varying covariates at t=0 and t=12
  "oxygenation_index_0", "norepinephrine_equivalent_0", "imv_status_0",
  "oxygenation_index_12", "norepinephrine_equivalent_12", "imv_status_12",
  # CCI components (baseline only)
  cci_vars
)

if (!all(required_vars %in% names(df))) {
  missing <- setdiff(required_vars, names(df))
  stop("Missing variables from data frame: ", paste(missing, collapse = ", "))
}

## ---- F. Collapse race into 3 categories ----
# Normalize to lowercase first so sites with different capitalization map correctly
df$race_category <- tolower(as.character(df$race_category))
df$race_category <- forcats::fct_collapse(
  factor(df$race_category),
  "White"                  = "white",
  "Black/African American" = "black or african american",
  other_level = "Other"
)
df <- droplevels(df)

# ================================ #
# ---- 1. DATA CLEANING ----
# ================================ #

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Generating Sample Characteristics\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

## ---- A. Define model variables for complete case analysis ----
model_vars <- required_vars

## ---- B. CRRT Dose Cutoff ----
# =================================== #

# CRRT dose cutoff (mL/kg/hr)
dose_cutoff <- 30

# Generate a label for filenames and titles
dose_label <- paste0("CRRT_", dose_cutoff, "cutoff")

cat("Using CRRT dose cutoff of:", dose_cutoff, "mL/kg/hr\n\n")

## ---- C. CRRT Dose Outlier Exclusion ----
# Adjustable cutoff — rows exceeding this are excluded                                                                                                                
CRRT_DOSE_MAX <- 300  # mL/kg/hr                                                                                                                                      

dose_cols <- c("crrt_dose_ml_kg_hr_0", "crrt_dose_0_12", "crrt_dose_12_24")                                                                                           
dose_cols <- dose_cols[dose_cols %in% names(df)]

outlier_mask <- rowSums(
  sapply(dose_cols, function(col) {
    !is.na(df[[col]]) & df[[col]] > CRRT_DOSE_MAX
  })
) > 0

n_outliers <- sum(outlier_mask)
if (n_outliers > 0) {
  cat("Excluding", n_outliers, "rows with CRRT dose >", CRRT_DOSE_MAX, "mL/kg/hr\n")
  cat("  Breakdown:\n")
  for (col in dose_cols) {
    n <- sum(!is.na(df[[col]]) & df[[col]] > CRRT_DOSE_MAX)
    if (n > 0) cat("   ", col, ":", n, "rows\n")
  }
  df <- df[!outlier_mask, ]
  cat("  Remaining rows:", nrow(df), "\n\n")
} else {
  cat("No CRRT dose outliers found (threshold:", CRRT_DOSE_MAX, "mL/kg/hr)\n\n")
}

## ---- D. Propensity Score Cutoffs for Truncation in SL Step 4A ----
PS_FLOOR <- 0.00                                                                                                                                                    
PS_CEIL  <- 1.00
#^Currently not truncating at all since it did not help

## ---- E. MICE Imputation and Missingness - df_complete ----

### ---- I. Report missingness before imputation ----
cat("Missing values per variable:\n")
miss_counts <- colSums(is.na(df[, model_vars]))
print(miss_counts[miss_counts > 0])
cat("\n")

# NOTE: 24h exclusion (died or off CRRT within 24h) now handled in Python step 04.

# MICE imputation for remaining missingness (oxygenation index only)
miss_counts <- colSums(is.na(df[, model_vars]))
vars_with_na <- names(miss_counts[miss_counts > 0])
n_incomplete <- sum(!complete.cases(df[, model_vars]))
if (n_incomplete > 0) {
  cat("Imputing", n_incomplete, "incomplete rows via MICE (pmm, m=5) …\n")
  library(mice)
  
  predictor_vars <- c("age_at_admission", "sex_category", "lactate_0",
                      "bicarbonate_0", "potassium_0", "norepinephrine_equivalent_0",
                      "imv_status_0", "crrt_dose_ml_kg_hr_0")
  mice_vars <- unique(c(vars_with_na, intersect(predictor_vars, model_vars)))
  mice_vars <- mice_vars[mice_vars %in% names(df)]
  
  imp <- mice(df[, mice_vars], m = 5, method = "pmm", seed = 42, maxit = 10, printFlag = FALSE)
  imputed_block <- complete(imp, 1)
  df_complete <- df
  df_complete[, mice_vars] <- imputed_block
  cat("MICE imputed", sum(miss_counts[miss_counts > 0]), "total missing values\n")
} else {
  cat("No missing values — skipping MICE\n")
  df_complete <- df
}

# Safety net: coerce NaN/Inf to NA in numeric columns (MICE doesn't handle these)
sanitize_numeric <- function(d) {
  num_cols <- names(d)[sapply(d, is.numeric)]
  for (col in num_cols) {
    bad <- is.nan(d[[col]]) | is.infinite(d[[col]])
    if (any(bad)) {
      cat("  Warning:", sum(bad), "NaN/Inf values in", col, "— converted to NA\n")
      d[[col]][bad] <- NA
    }
  }
  d
}
cat("Checking for NaN/Inf values after imputation...\n")
df_complete <- sanitize_numeric(df_complete)
cat("NaN/Inf check complete.\n")

# Drop any rows still incomplete in model_vars after MICE + sanitize
n_still_incomplete <- sum(!complete.cases(df_complete[, model_vars]))
if (n_still_incomplete > 0) {
  cat("Warning:", n_still_incomplete, "rows still incomplete after MICE + sanitize — dropping\n")
  df_complete <- df_complete[complete.cases(df_complete[, model_vars]), ]
  df_complete <- droplevels(df_complete)
}

cat("Analysis sample:", nrow(df_complete), "\n\n")

# Check for sufficient data
if (nrow(df_complete) < 50) {
  stop("Insufficient cases for modeling (n = ", nrow(df_complete), ")")
}

# ================================ #
# ---- 2. DATA SUMMARIZATION ----
# ================================ #

## ---- A. Calculate summary statistics ----
sample_chars <- data.frame(
  site_name = SITE_NAME,
  total_n = nrow(df_complete),

  # Age
  age_mean = mean(df_complete$age_at_admission, na.rm = TRUE),
  age_sd = sd(df_complete$age_at_admission, na.rm = TRUE),
  age_median = median(df_complete$age_at_admission, na.rm = TRUE),
  age_q25 = quantile(df_complete$age_at_admission, 0.25, na.rm = TRUE),
  age_q75 = quantile(df_complete$age_at_admission, 0.75, na.rm = TRUE),

  # Sex
  sex_male_n = sum(df_complete$sex_category == "male", na.rm = TRUE),
  sex_male_pct = mean(df_complete$sex_category == "male", na.rm = TRUE) * 100,
  
  # Weight (kg)                                                                                                                                                         
  weight_mean = mean(df_complete$weight_kg, na.rm = TRUE),                                                                                                              
  weight_sd = sd(df_complete$weight_kg, na.rm = TRUE),                                                                                                                  
  weight_median = median(df_complete$weight_kg, na.rm = TRUE),                                                                                                          
  weight_q25 = quantile(df_complete$weight_kg, 0.25, na.rm = TRUE),                                                                                                     
  weight_q75 = quantile(df_complete$weight_kg, 0.75, na.rm = TRUE),
  weight_min = min(df_complete$weight_kg, na.rm = TRUE),
  weight_max = max(df_complete$weight_kg, na.rm = TRUE),

  # SOFA score (descriptive only)
  sofa_mean = mean(df_complete$sofa_total_0, na.rm = TRUE),
  sofa_sd = sd(df_complete$sofa_total_0, na.rm = TRUE),
  sofa_median = median(df_complete$sofa_total_0, na.rm = TRUE),
  sofa_q25 = quantile(df_complete$sofa_total_0, 0.25, na.rm = TRUE),
  sofa_q75 = quantile(df_complete$sofa_total_0, 0.75, na.rm = TRUE),

  # CRRT dose at initiation (t=0)                                                                                                                                     
  crrt_dose_t0_mean = mean(df_complete$crrt_dose_ml_kg_hr_0, na.rm = TRUE),                                                                                           
  crrt_dose_t0_sd = sd(df_complete$crrt_dose_ml_kg_hr_0, na.rm = TRUE),                                                                                               
  crrt_dose_t0_median = median(df_complete$crrt_dose_ml_kg_hr_0, na.rm = TRUE),                                                                                       
  crrt_dose_t0_q25 = quantile(df_complete$crrt_dose_ml_kg_hr_0, 0.25, na.rm = TRUE),                                                                                  
  crrt_dose_t0_q75 = quantile(df_complete$crrt_dose_ml_kg_hr_0, 0.75, na.rm = TRUE),
  crrt_dose_t0_min = min(df_complete$crrt_dose_ml_kg_hr_0, na.rm = TRUE),
  crrt_dose_t0_max = max(df_complete$crrt_dose_ml_kg_hr_0, na.rm = TRUE),
  
  # CRRT dose mean 0-12h
  crrt_dose_0_12_mean = mean(df_complete$crrt_dose_0_12, na.rm = TRUE),
  crrt_dose_0_12_sd = sd(df_complete$crrt_dose_0_12, na.rm = TRUE),
  crrt_dose_0_12_median = median(df_complete$crrt_dose_0_12, na.rm = TRUE),
  crrt_dose_0_12_q25 = quantile(df_complete$crrt_dose_0_12, 0.25, na.rm = TRUE),
  crrt_dose_0_12_q75 = quantile(df_complete$crrt_dose_0_12, 0.75, na.rm = TRUE),
  crrt_dose_0_12_min = min(df_complete$crrt_dose_0_12, na.rm = TRUE),
  crrt_dose_0_12_max = max(df_complete$crrt_dose_0_12, na.rm = TRUE),
  
  # CRRT dose mean 12-24h
  crrt_dose_12_24_mean = mean(df_complete$crrt_dose_12_24, na.rm = TRUE),
  crrt_dose_12_24_sd = sd(df_complete$crrt_dose_12_24, na.rm = TRUE),
  crrt_dose_12_24_median = median(df_complete$crrt_dose_12_24, na.rm = TRUE),
  crrt_dose_12_24_q25 = quantile(df_complete$crrt_dose_12_24, 0.25, na.rm = TRUE),
  crrt_dose_12_24_q75 = quantile(df_complete$crrt_dose_12_24, 0.75, na.rm = TRUE),
  crrt_dose_12_24_min = min(df_complete$crrt_dose_12_24, na.rm = TRUE),
  crrt_dose_12_24_max = max(df_complete$crrt_dose_12_24, na.rm = TRUE),

  # Oxygenation Index (time t=0)
  oxygenation_index_mean   = mean(df_complete$oxygenation_index_0, na.rm = TRUE),
  oxygenation_index_sd     = sd(df_complete$oxygenation_index_0, na.rm = TRUE),
  oxygenation_index_median = median(df_complete$oxygenation_index_0, na.rm = TRUE),
  oxygenation_index_q25    = quantile(df_complete$oxygenation_index_0, 0.25, na.rm = TRUE),
  oxygenation_index_q75    = quantile(df_complete$oxygenation_index_0, 0.75, na.rm = TRUE),

  # Norepinephrine Equivalent (time t=0)
  norepi_eq_mean   = mean(df_complete$norepinephrine_equivalent_0, na.rm = TRUE),
  norepi_eq_sd     = sd(df_complete$norepinephrine_equivalent_0, na.rm = TRUE),
  norepi_eq_median = median(df_complete$norepinephrine_equivalent_0, na.rm = TRUE),
  norepi_eq_q25    = quantile(df_complete$norepinephrine_equivalent_0, 0.25, na.rm = TRUE),
  norepi_eq_q75    = quantile(df_complete$norepinephrine_equivalent_0, 0.75, na.rm = TRUE),

  # IMV Status (time t=0)
  imv_status_n   = sum(df_complete$imv_status_0 == 1, na.rm = TRUE),
  imv_status_pct = mean(df_complete$imv_status_0 == 1, na.rm = TRUE) * 100,

  # Lactate (time t=0)
  lactate_mean   = mean(df_complete$lactate_0, na.rm = TRUE),
  lactate_sd     = sd(df_complete$lactate_0, na.rm = TRUE),
  lactate_median = median(df_complete$lactate_0, na.rm = TRUE),
  lactate_q25    = quantile(df_complete$lactate_0, 0.25, na.rm = TRUE),
  lactate_q75    = quantile(df_complete$lactate_0, 0.75, na.rm = TRUE),

  # Bicarbonate (time t=0)
  bicarbonate_mean   = mean(df_complete$bicarbonate_0, na.rm = TRUE),
  bicarbonate_sd     = sd(df_complete$bicarbonate_0, na.rm = TRUE),
  bicarbonate_median = median(df_complete$bicarbonate_0, na.rm = TRUE),
  bicarbonate_q25    = quantile(df_complete$bicarbonate_0, 0.25, na.rm = TRUE),
  bicarbonate_q75    = quantile(df_complete$bicarbonate_0, 0.75, na.rm = TRUE),

  # Potassium (time t=0)
  potassium_mean   = mean(df_complete$potassium_0, na.rm = TRUE),
  potassium_sd     = sd(df_complete$potassium_0, na.rm = TRUE),
  potassium_median = median(df_complete$potassium_0, na.rm = TRUE),
  potassium_q25    = quantile(df_complete$potassium_0, 0.25, na.rm = TRUE),
  potassium_q75    = quantile(df_complete$potassium_0, 0.75, na.rm = TRUE),

  # Outcomes
  outcome_censored_n = sum(df_complete$outcome == 0, na.rm = TRUE),
  outcome_censored_pct = mean(df_complete$outcome == 0, na.rm = TRUE) * 100,
  outcome_discharge_n = sum(df_complete$outcome == 1, na.rm = TRUE),
  outcome_discharge_pct = mean(df_complete$outcome == 1, na.rm = TRUE) * 100,
  outcome_death_n = sum(df_complete$outcome == 2, na.rm = TRUE),
  outcome_death_pct = mean(df_complete$outcome == 2, na.rm = TRUE) * 100,

  # CRRT duration
  crrt_duration_days_mean   = mean(df_complete$crrt_duration_days, na.rm = TRUE),
  crrt_duration_days_sd     = sd(df_complete$crrt_duration_days, na.rm = TRUE),
  crrt_duration_days_median = median(df_complete$crrt_duration_days, na.rm = TRUE),
  crrt_duration_days_q25    = quantile(df_complete$crrt_duration_days, 0.25, na.rm = TRUE),
  crrt_duration_days_q75    = quantile(df_complete$crrt_duration_days, 0.75, na.rm = TRUE),

  # IMV duration
  imv_duration_days_mean   = mean(df_complete$imv_duration_days, na.rm = TRUE),
  imv_duration_days_sd     = sd(df_complete$imv_duration_days, na.rm = TRUE),
  imv_duration_days_median = median(df_complete$imv_duration_days, na.rm = TRUE),
  imv_duration_days_q25    = quantile(df_complete$imv_duration_days, 0.25, na.rm = TRUE),
  imv_duration_days_q75    = quantile(df_complete$imv_duration_days, 0.75, na.rm = TRUE),

  stringsAsFactors = FALSE
)

# Save sample characteristics
write.csv(sample_chars,
          file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                       "_sample_characteristics.csv")),
          row.names = FALSE)

cat("Sample characteristics saved.\n\n")

## ---- B. Histogram CRRT dose ----

### ---- I. Separate Histograms of CRRT dose ----
# Helper function to avoid repeating ggplot code
plot_dose_hist <- function(data, dose_col, title_label) {                                                                                                             
  ggplot(data, aes(x = .data[[dose_col]])) +
    geom_histogram(binwidth = 5, color = "black", fill = "skyblue") +
    geom_vline(xintercept = dose_cutoff, linetype = "dashed", linewidth = 1) +
    labs(
      title = paste0("Distribution of ", title_label,
                     " (Cutoff = ", dose_cutoff, " mL/kg/hr)"),
      x = "CRRT Dose (mL/kg/hr)",
      y = "Count"
    ) +
    theme_bw(base_size = 12)
}

# 1. Initial dose (t=0)
crrt_dose_hist_t0 <- plot_dose_hist(
  df_complete, "crrt_dose_ml_kg_hr_0", "Initial CRRT Dose (t=0)")
crrt_dose_hist_t0

ggsave(file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                    "_hist_crrt_dose_t0.png")),
       crrt_dose_hist_t0, width = 6, height = 4)

# 2. Mean dose 0-12h
crrt_dose_hist_0_12 <- plot_dose_hist(
  df_complete, "crrt_dose_0_12", "Mean CRRT Dose (0-12h)")
crrt_dose_hist_0_12

ggsave(file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                    "_hist_crrt_dose_0_12.png")),
       crrt_dose_hist_0_12, width = 6, height = 4)

# 3. Mean dose 12-24h
crrt_dose_hist_12_24 <- plot_dose_hist(
  df_complete, "crrt_dose_12_24", "Mean CRRT Dose (12-24h)")
crrt_dose_hist_12_24

ggsave(file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                    "_hist_crrt_dose_12_24.png")),
       crrt_dose_hist_12_24, width = 6, height = 4)

cat("Saved all three dose histogram PNGs\n")

### ---- II. Combined Histogram of CRRT dose ----
# Pivot to long format for overlaid plotting                                                                                                                          
dose_long <- df_complete %>%                                                                                                                                          
  select(encounter_block, crrt_dose_ml_kg_hr_0, crrt_dose_0_12, crrt_dose_12_24) %>%                                                                                  
  pivot_longer(                                                                                                                                                       
    cols = c(crrt_dose_ml_kg_hr_0, crrt_dose_0_12, crrt_dose_12_24),                                                                                                  
    names_to = "window",
    values_to = "dose"
  ) %>%
  mutate(
    window = factor(window,
                    levels = c("crrt_dose_ml_kg_hr_0", "crrt_dose_0_12", "crrt_dose_12_24"),
                    labels = c("At Initiation", "Mean 0-12h", "Mean 12-24h")
    )
  )

crrt_dose_hist_combined <- ggplot(dose_long, aes(x = dose, fill = window)) +
  geom_histogram(binwidth = 5, alpha = 0.5, position = "identity",
                 color = "black", linewidth = 0.2) +
  geom_vline(xintercept = dose_cutoff, linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c(
    "At Initiation" = "#E69F00",
    "Mean 0-12h"    = "#56B4E9",
    "Mean 12-24h"   = "#009E73"
  )) +
  labs(
    title = paste0("CRRT Dose Distribution Over First 24h",
                   " (Cutoff = ", dose_cutoff, " mL/kg/hr)"),
    x = "CRRT Dose (mL/kg/hr)",
    y = "Count",
    fill = "Time Window"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

crrt_dose_hist_combined

ggsave(file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                    "_hist_crrt_dose_combined.png")),
       crrt_dose_hist_combined, width = 8, height = 5)
cat("Saved combined dose histogram PNG\n")

## ---- C. Designate High vs Low CRRT Dose ----
cat("Defining treatment groups using cutoff =", dose_cutoff, "mL/kg/hr\n")

# Dynamic complete-case variable list: all covariates that models will touch
complete_case_vars <- c(
  "age_at_admission", "sex_category", "race_category", "weight_kg",
  "sofa_total_0", "crrt_dose_ml_kg_hr_0",
  "lactate_0", "bicarbonate_0", "potassium_0",
  "oxygenation_index_0", "norepinephrine_equivalent_0", "imv_status_0",
  cci_vars,
  "time_to_event_90d", "outcome"
)
# Keep only vars that exist in the data (sites without CCI columns, etc.)
complete_case_vars <- intersect(complete_case_vars, names(df_complete))

df_tte_bin <- df_complete %>%
  mutate(
    crrt_high = ifelse(crrt_dose_ml_kg_hr_0 >= dose_cutoff, 1L, 0L),
    crrt_high = factor(
      crrt_high,
      levels = c(0,1),
      labels = c(paste0("<", dose_cutoff), paste0(">=", dose_cutoff))
    )
  ) %>%
  drop_na(all_of(complete_case_vars)) %>%
  droplevels()

cat("After drop_na on all model covariates:", nrow(df_tte_bin), "rows\n")

# Print Distribution
df_tte_bin %>%
  count(crrt_high) %>%
  mutate(prop = n/sum(n))

# Save df_tte_bin treatment distribution
bin_dist <- df_tte_bin %>%
  count(crrt_high) %>%
  mutate(prop = n/sum(n))
write.csv(bin_dist, file.path(output_dir,
                              paste0(SITE_NAME, "_", dose_label,
                                     "_crrt_bin_distribution.csv")),
          row.names = FALSE)
cat("Saved bin distribution\n")


## ---- D. Table 1 ----

### ---- Modify df to include all three outcomes in one column
df_tte_table1 <- df_tte_bin %>%
  mutate(
    #### ---- Rename CRRT groups
    crrt_group = factor(
      crrt_high,
      levels = c(paste0("<", dose_cutoff),
                 paste0(">=", dose_cutoff)),
      labels = c(
        paste0("Low CRRT dose (<", dose_cutoff, " mL/kg/hr)"),
        paste0("High CRRT dose (>=", dose_cutoff, " mL/kg/hr)")
      )
    ),

    #### ---- Ensure sex is a factor
    sex_category = factor(sex_category),

    #### ---- Clean race labels
    race_category = factor(race_category),
    
    #### ---- Clean CRRT modality labels
    crrt_mode_category = factor(crrt_mode_category),

    #### ---- Three-category outcome
    outcome_3cat = factor(
      case_when(
        outcome == 1 ~ "Discharged",
        outcome == 2 ~ "Died",
        TRUE ~ "Censored"
      ),
      levels = c("Discharged",
                 "Died",
                 "Censored")
    )
  )

### ---- Build Table 1

#### ---- Helper function to collapse categorical variables like Sex ----
collapse_binary_in_gtsummary <- function(tbl, variable, keep_level, label) {

  modify_table_body(tbl, function(body) {

    # Locate rows belonging to this variable
    idx <- which(body$variable == variable)

    # Split into target variable rows vs all others
    sub <- body[idx, ]
    rest <- body[-idx, ]

    # Detect p-value column
    pcol <- intersect(c("stat_pvalue", "p.value", "stat_p", "pvalue"), names(sub))
    if (length(pcol) == 0) stop("No p-value column found.")
    pcol <- pcol[1]

    # Extract p-value from the original header row
    pval <- sub[[pcol]][sub$row_type == "label"][1]

    # Build the replacement 1-row block
    sub2 <- sub %>%
      # remove header and unwanted levels
      dplyr::filter(!(row_type == "label"),
                    !(row_type == "level" & label != keep_level)) %>%
      # promote remaining level row to label row
      dplyr::mutate(
        row_type = "label",
        label = !!label,
        !!pcol := pval
      )

    # Insert the updated row back in the original location
    out <- dplyr::bind_rows(
      rest[seq_len(min(idx) - 1), ],   # before original female position
      sub2,                            # updated female row
      rest[seq(min(idx), nrow(rest)),] # remaining rows
    )

    out
  })
}

#### Include covariates
vars_table1 <- c(
  "age_at_admission",
  "sex_category",
  "weight_kg",
  "race_category",
  "sofa_total_0",
  "creatinine_0",
  "oxygenation_index_0",
  "norepinephrine_equivalent_0",
  "imv_status_0",
  "lactate_0",
  "bicarbonate_0",
  "potassium_0",
  "crrt_mode_category",
  "crrt_duration_days",
  "imv_duration_days",
  "crrt_dose_ml_kg_hr_0",
  cci_vars,
  "outcome_3cat"
)

# Build type list (CCI vars auto-detected as dichotomous from 0/1 integers)
table1_type <- list(
  age_at_admission            ~ "continuous",
  sex_category                ~ "categorical",
  weight_kg                   ~ "continuous",
  race_category               ~ "categorical",
  sofa_total_0                ~ "continuous",
  creatinine_0                ~ "continuous",
  oxygenation_index_0         ~ "continuous",
  norepinephrine_equivalent_0 ~ "continuous",
  imv_status_0                ~ "dichotomous",
  lactate_0                   ~ "continuous",
  bicarbonate_0               ~ "continuous",
  potassium_0                 ~ "continuous",
  crrt_mode_category          ~ "categorical",
  crrt_duration_days          ~ "continuous",
  imv_duration_days           ~ "continuous",
  crrt_dose_ml_kg_hr_0        ~ "continuous"
)
# Add CCI vars explicitly as dichotomous
for (v in cci_vars) {
  table1_type[[length(table1_type) + 1]] <- as.formula(
    paste0(v, ' ~ "dichotomous"'))
}

# Build label list
table1_label <- list(
  age_at_admission             ~ "Age at Admission (years)",
  sex_category                 ~ "Female (%)",
  weight_kg                    ~ "Weight (kg)",
  race_category                ~ "Race",
  sofa_total_0                 ~ "SOFA Score",
  creatinine_0                 ~ "Creatinine at CRRT Start (mg/dL)",
  oxygenation_index_0          ~ "Oxygenation Index (P/F or S/F)",
  norepinephrine_equivalent_0  ~ "NE Equivalent (mcg/kg/min)",
  imv_status_0                 ~ "On IMV (%)",
  lactate_0                    ~ "Lactate at CRRT Start (mmol/L)",
  bicarbonate_0                ~ "Bicarbonate at CRRT Start (mEq/L)",
  potassium_0                  ~ "Potassium at CRRT Start (mEq/L)",
  crrt_mode_category           ~ "CRRT Modality",
  crrt_duration_days           ~ "Duration of CRRT (Days)",
  imv_duration_days            ~ "Duration of IMV (Days)",
  crrt_dose_ml_kg_hr_0         ~ "Initial CRRT Dose (mL/kg/hr)",
  outcome_3cat                 ~ "90-day Outcome"
)
# Add CCI labels
for (v in cci_vars) {
  table1_label[[length(table1_label) + 1]] <- as.formula(
    paste0(v, ' ~ "', cci_labels[v], '"'))
}

# Auto-detect all dichotomous variables and set value=1 for gtsummary
# Drop dichotomous vars that are all-zero (no events) — gtsummary crashes on value=1 if 1 doesn't exist
table1_value <- list()
drop_vars <- character(0)
for (i in seq_along(table1_type)) {
  if (grepl("dichotomous", deparse(table1_type[[i]]))) {
    vname <- all.vars(table1_type[[i]])[1]
    if (vname %in% names(df_tte_table1) && !any(df_tte_table1[[vname]] == 1L, na.rm = TRUE)) {
      drop_vars <- c(drop_vars, vname)
    } else {
      table1_value[[length(table1_value) + 1]] <- as.formula(paste0(vname, " ~ 1L"))
    }
  }
}
if (length(drop_vars) > 0) {
  cat("  Dropping all-zero dichotomous vars from Table 1:", paste(drop_vars, collapse = ", "), "\n")
  vars_table1 <- setdiff(vars_table1, drop_vars)
  table1_type <- table1_type[!sapply(table1_type, function(f) all.vars(f)[1] %in% drop_vars)]
  table1_label <- table1_label[!sapply(table1_label, function(f) all.vars(f)[1] %in% drop_vars)]
}

table1 <- df_tte_table1 %>%
  select(crrt_group, all_of(vars_table1)
  ) %>%
  tbl_summary(
    by = crrt_group,
    type = table1_type,
    value = table1_value,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = table1_label
  ) %>%
  add_overall() %>%
  {
    tryCatch(
      add_p(., test = list(
        race_category ~ "chisq.test.no.correct",
        outcome_3cat  ~ "chisq.test.no.correct"
      )),
      error = function(e) {
        cat("  Warning: add_p() failed (", conditionMessage(e), "), skipping p-values\n")
        .
      }
    )
  } %>%
  bold_labels() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(c("stat_1","stat_2") ~ "**CRRT dose group**") %>%
  modify_caption(
    "**Table 1. Baseline Characteristics of the Cohort by Initial CRRT Dose Group**")

table1 <- table1 %>%
  collapse_binary_in_gtsummary(
    variable   = "sex_category",
    keep_level = "female",
    label      = "Female (%)"
  )

# Save Table 1 as HTML
gt::gtsave(
  as_gt(table1),
  filename = file.path(output_dir,
                       paste0(SITE_NAME, "_",
                              dose_label, "_Table1_unadjusted.html"))
)
write.csv(as_tibble(table1),
          file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_Table1_unadjusted.csv")),
          row.names = FALSE)


# ============================================================= #
# ---- 3. CREATING WIDE DATASET ----
# ============================================================= #

## ---- A. Making and visualizing wide MSM datasets ----

# Add A_0 and A_12 labels
df_msm_wide <- df_complete %>%
  mutate(
    A_0  = as.integer(crrt_dose_0_12  >= dose_cutoff),
      # treatment for interval 0-12h
    A_12 = as.integer(crrt_dose_12_24 >= dose_cutoff)
      # treatment for interval 12-24h
  )

# keep a labeled factor version for plots/tables
df_msm_wide <- df_msm_wide %>%
  mutate(
    A_0_f  = factor(A_0,  levels=c(0,1),
                    labels=c(paste0("<",dose_cutoff), paste0(">=",dose_cutoff))),
    A_12_f = factor(A_12, levels=c(0,1),
                    labels=c(paste0("<",dose_cutoff), paste0(">=",dose_cutoff)))
  )

## ---- B. Pivot to long dataset ----

id_var <- "encounter_block"

df_msm_long <- bind_rows(
  df_msm_wide %>%
    transmute(
      encounter_block = .data[[id_var]],
      tpt = 0L,
      A = A_0,
      # baseline covariates (carried on both rows)
      age_at_admission, sex_category, race_category, ethnicity_category,
      crrt_mode_category, weight_kg,
      # time-varying covariates at decision time
      lactate = lactate_0,
      bicarbonate = bicarbonate_0,
      potassium = potassium_0,
      oxygenation_index = oxygenation_index_0,
      norepinephrine_equivalent = norepinephrine_equivalent_0,
      imv_status = imv_status_0,
      # SOFA (descriptive only, not used in models)
      sofa_total = sofa_total_0,
      # CCI components (baseline only, same on both rows)
      across(all_of(cci_vars)),
      # outcome (patient-level, repeated here for convenience)
      time_to_event_90d, outcome, censored_at_90d
    ),
  df_msm_wide %>%
    transmute(
      encounter_block = .data[[id_var]],
      tpt = 12L,
      A = A_12,
      age_at_admission, sex_category, race_category, ethnicity_category,
      crrt_mode_category, weight_kg,
      lactate = lactate_12,
      bicarbonate = bicarbonate_12,
      potassium = potassium_12,
      oxygenation_index = oxygenation_index_12,
      norepinephrine_equivalent = norepinephrine_equivalent_12,
      imv_status = imv_status_12,
      sofa_total = sofa_total_12,
      across(all_of(cci_vars)),
      time_to_event_90d, outcome, censored_at_90d
    )
) %>%
  arrange(encounter_block, tpt) %>%
  group_by(encounter_block) %>%
  mutate(
    A_lag = lag(A),
    A_lag = if_else(is.na(A_lag), 0L, A_lag)
  ) %>%
  ungroup()

stopifnot(all(df_msm_long$A %in% c(0L, 1L)))

# ============================================================= #
# ---- 4. STABILIZED WEIGHTS WITH SUPERLEARNER ----
# ============================================================= #

## ---- A. Propensity Weighting with SuperLearner  ----
set.seed(42)

# Baseline terms: demographics + CCI components
# NOTE: sofa_total deliberately excluded from model covariates
baseline_terms <- c(
  "age_at_admission", "sex_category", "race_category", "weight_kg",
  cci_vars
)

# Time-varying terms: labs + new acute severity markers
# NOTE: sofa_total deliberately excluded from model covariates
tv_terms <- c(
  "lactate", "bicarbonate", "potassium",
  "oxygenation_index", "norepinephrine_equivalent", "imv_status"
)

# Dynamic filter: drop baseline terms with <2 unique values at this site
# (e.g., CCI components with zero prevalence)
baseline_terms <- baseline_terms[
  sapply(baseline_terms, function(v) {
    length(unique(df_msm_long[[v]])) >= 2
  })
]

cat("Baseline terms after dynamic filter (", length(baseline_terms), "):\n")
cat(" ", paste(baseline_terms, collapse = ", "), "\n\n")

sl_lib <- c("SL.glm","SL.gam","SL.randomForest")

sl_fit_prob <- function(d, y, xvars, sl_lib) {
  # Complete-case filter on y + xvars to prevent model.matrix row-dropping mismatch
  keep <- complete.cases(d[, c(y, xvars)])
  n_drop <- sum(!keep)
  if (n_drop > 0) cat("  SL: dropping", n_drop, "incomplete rows\n")
  d_cc <- d[keep, ]

  X <- d_cc %>% select(all_of(xvars))
  mm <- model.matrix(~ . , data = X)[, -1, drop = FALSE]
  colnames(mm) <- make.names(colnames(mm))

  preds_cc <- tryCatch({
    fit <- SuperLearner::SuperLearner(
      Y = d_cc[[y]],
      X = as.data.frame(mm),
      family = binomial(),
      SL.library = sl_lib
    )
    as.vector(fit$SL.predict)
  }, error = function(e) {
    cat("  WARNING: SuperLearner failed:", conditionMessage(e), "\n")
    cat("  Falling back to GLM\n")
    glm_fit <- glm(
      reformulate(colnames(mm), response = y),
      data = cbind(d_cc[, y, drop = FALSE], as.data.frame(mm)),
      family = binomial()
    )
    predict(glm_fit, type = "response")
  })

  # Expand back to full length — NA for dropped rows
  out <- rep(NA_real_, nrow(d))
  out[keep] <- preds_cc
  out
}

df_w <- df_msm_long %>%
  mutate(p_denom = NA_real_, p_num = NA_real_, ps_denom = NA_real_)

for (tt in c(0, 12)) {

  dft <- df_w %>% filter(tpt == tt)

  # Denominator: P(A_t | baseline + current L_t + prior A)
  # Numerator:   P(A_t | baseline + prior A)  (stabilization)
  if (tt == 0) {
    denom_vars <- c(baseline_terms, tv_terms)
    num_vars   <- c(baseline_terms)
  } else {
    denom_vars <- c(baseline_terms, tv_terms, "A_lag")
    num_vars   <- c(baseline_terms, "A_lag")
  }

  p1_denom <- sl_fit_prob(dft, y = "A", xvars = denom_vars, sl_lib = sl_lib)
  p1_num   <- sl_fit_prob(dft, y = "A", xvars = num_vars,   sl_lib = sl_lib)

  # Fill NA predictions (from dropped incomplete rows) with marginal P(A=1)
  # so those patients get uninformative weights (~1)
  marginal_p <- mean(dft$A, na.rm = TRUE)
  p1_denom[is.na(p1_denom)] <- marginal_p
  p1_num[is.na(p1_num)]     <- marginal_p

  ### ---- I. Truncate extreme propensity scores to enforce positivity  ----
  # If unwanted, add #s
  p1_denom <- pmin(pmax(p1_denom, PS_FLOOR), PS_CEIL)

  # Probability of the observed treatment A
  p_obs_denom <- ifelse(dft$A == 1, p1_denom, 1 - p1_denom)
  p_obs_num   <- ifelse(dft$A == 1, p1_num,   1 - p1_num)

  idx <- which(df_w$tpt == tt)
  df_w$p_denom[idx]  <- p_obs_denom
  df_w$p_num[idx]    <- p_obs_num
  df_w$ps_denom[idx] <- p1_denom
}

# Cumulative stabilized weight through the treatment window (0-24h)
df_w <- df_w %>%
  group_by(encounter_block) %>%
  arrange(tpt, .by_group = TRUE) %>%
  mutate(w_msm = cumprod(p_num / p_denom)) %>%
  ungroup()

# Final weight for outcome model = weight after t=12 decision (covers 0-12 and 12-24)
df_w_final <- df_w %>%
  filter(tpt == 12) %>%
  select(encounter_block, w_msm)

## ---- B. Weight and Overlap Diagnostics ----

### ---- I. Weight summaries (final MSM weights at t=12) ----
summary(df_w_final$w_msm)
quantile(df_w_final$w_msm, c(0, .01, .05, .1, .5, .9, .95, .99, 1), na.rm = TRUE)

ESS <- (sum(df_w_final$w_msm)^2) / sum(df_w_final$w_msm^2)
ESS_prop <- ESS / nrow(df_w_final)

write.csv(
  tibble(N = nrow(df_w_final), ESS = ESS, ESS_prop = ESS_prop),
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_MSM_IPTW_ESS.csv"))
)

cat("ESS:", ESS, "\n")
cat("ESS_prop:", ESS_prop, "\n")
cat("Output path:", file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_MSM_IPTW_ESS.csv")), "\n")
cat("Output dir exists:", dir.exists(output_dir), "\n")

### ---- II. Weight Overlap Diagnostics (Positivity per Timepoint) ----

# Safety check: ps_denom must exist (created in the sequential weighting loop)
if (!("ps_denom" %in% names(df_w))) {
  stop("df_w is missing 'ps_denom'.
       Add ps_denom storage in the MSM weighting loop (ps_denom = p1_denom).")
}

plot_overlap_timepoint <- function(
    df_w, tt, output_dir, SITE_NAME, dose_label, dose_cutoff) {

  dft <- df_w %>%
    dplyr::filter(tpt == tt) %>%
    dplyr::mutate(
      trt = factor(A, levels = c(0, 1),
                   labels = c(paste0("<", dose_cutoff), paste0(">=", dose_cutoff)))
    )

  plot_overlap_timepoint_gg <- ggplot(dft, aes(x = ps_denom, fill = trt)) +
    geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
    theme_bw(base_size = 10) +
    labs(
      title = paste0("Denominator Propensity Score Overlap at t=", tt, "h (MSM)"),
      subtitle = "Distribution of P(High dose | history) by observed dose group",
      x = "Predicted P(High dose | history) [denominator model]",
      y = "Count",
      fill = "Observed dose"
    )

  out_path <- file.path(output_dir,
                        paste0(SITE_NAME, "_", dose_label,
                               "_MSM_PS_overlap_t", tt, ".png"))
  ggsave(out_path, plot_overlap_timepoint_gg, width = 6, height = 4, dpi = 300)

  cat("Saved MSM PS overlap plot for t=", tt, "h\n", sep = "")
  return(plot_overlap_timepoint_gg)
}

# Save overlap plots for each decision time
plot_overlap_t0  <- plot_overlap_timepoint(
  df_w, 0,  output_dir, SITE_NAME, dose_label, dose_cutoff)
plot_overlap_t12 <- plot_overlap_timepoint(
  df_w, 12, output_dir, SITE_NAME, dose_label, dose_cutoff)

plot_overlap_t0
plot_overlap_t12

# Optional numeric screen for extreme PS (positivity red flags)
extreme_summary <- df_w %>%
  dplyr::group_by(tpt) %>%
  dplyr::summarise(
    n = dplyr::n(),
    ps_lt_0.01 = mean(ps_denom < 0.01, na.rm = TRUE),
    ps_gt_0.99 = mean(ps_denom > 0.99, na.rm = TRUE),
    ps_lt_0.05 = mean(ps_denom < 0.05, na.rm = TRUE),
    ps_gt_0.95 = mean(ps_denom > 0.95, na.rm = TRUE)
  )

readr::write_csv(
  extreme_summary,
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_MSM_PS_extremes.csv"))
)

print(extreme_summary)

### ---- III. Balance Diagnostics (including Love Plots) ----

# Pretty variable names for balance plots
pretty_names_msm <- c(
  age_at_admission          = "Age at Admission",
  weight_kg                 = "Weight (kg)",
  sex_category              = "Sex",
  race_category             = "Race",
  ethnicity_category        = "Ethnicity",
  crrt_mode_category        = "CRRT Modality",
  lactate                   = "Lactate",
  bicarbonate               = "Bicarbonate",
  potassium                 = "Potassium",
  oxygenation_index         = "Oxygenation Index",
  norepinephrine_equivalent = "NE Equivalent",
  imv_status                = "IMV Status",
  A_lag                     = "Prior Dose Group",
  cci_labels  # named vector: cci_myocardial_infarction = "Myocardial Infarction", etc.
)

# Balance check function
run_balance_timepoint <- function(df_w, tt, output_dir, SITE_NAME,
                                  dose_label, dose_cutoff, pretty_names,
                                  baseline_terms, tv_terms) {

  dft <- df_w %>%
    dplyr::filter(tpt == tt) %>%
    dplyr::mutate(
      trt = factor(A, levels = c(0,1),
                   labels = c(paste0("<", dose_cutoff), paste0(">=", dose_cutoff)))
    )

  # Covariates to assess balance on at that timepoint
  covs <- c(baseline_terms, tv_terms)
  if (tt != 0) covs <- c(covs, "A_lag")

  bal <- cobalt::bal.tab(
    x = dft[, covs, drop = FALSE],
    treat = dft$trt,
    weights = dft$w_msm,
    method = "weighting",
    s.d.denom = "pooled",
    un = TRUE
  )

  bal_df <- as.data.frame(bal$Balance)

  # Save balance as CSV
  csv_path <- file.path(output_dir, paste0(
    SITE_NAME, "_", dose_label, "_MSM_balance_t", tt, ".csv"))
  write.csv(bal_df, csv_path, row.names = TRUE)

  # Love plot
  plot_loveplot_msm <- cobalt::love.plot(
    bal,
    stats = "mean.diffs",
    abs = TRUE,
    var.names = pretty_names,
    thresholds = c(m = 0.1),
    var.order = "unadjusted",
    line.size = 0.8,
    point.size = 3,
    sample.names = c("Unweighted", "Weighted"),
    title = paste0("Covariate Balance: MSM IPTW at t=", tt, "h"),
    subtitle = "Standardized Mean Differences Before and After Weighting",
    grid = TRUE
  )

  print(plot_loveplot_msm)

  # Save Love plot
  png_path <- file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                           "_MSM_LovePlot_t", tt, ".png"))
  png(png_path, width = 8, height = 7, units = "in", res = 300)
  print(plot_loveplot_msm)
  dev.off()

  # ESS at this timepoint
  ESS <- (sum(dft$w_msm)^2) / sum(dft$w_msm^2)
  ESS_prop <- ESS / nrow(dft)

  ess_path <- file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                           "_MSM_ESS_t", tt, ".csv"))
  readr::write_csv(
    tibble::tibble(tpt = tt, N = nrow(dft), ESS = ESS, ESS_prop = ESS_prop),
    ess_path
  )

  cat("Saved balance + love plot + ESS for t=", tt, "h\n", sep = "")
  invisible(list(bal = bal, bal_df = bal_df, love_plot = plot_loveplot_msm))
}

### ---- IV. Run balance checks at specific decision times ----
bal_t0  <- run_balance_timepoint(
  df_w, tt = 0,  output_dir, SITE_NAME, dose_label,
  dose_cutoff, pretty_names_msm, baseline_terms, tv_terms
)
bal_t12 <- run_balance_timepoint(
  df_w, tt = 12, output_dir, SITE_NAME, dose_label,
  dose_cutoff, pretty_names_msm, baseline_terms, tv_terms
)

# ============================================================= #
# ---- 5. MSM OUTCOME MODELS ----
# ============================================================= #

## ---- A. Exposure history summary: number of high-dose intervals (0,1,2) ----
A_hist <- df_msm_wide %>%
  transmute(
    encounter_block,
    high_count = A_0 + A_12,
    always_high = as.integer(high_count == 2),
    always_low  = as.integer(high_count == 0)
  )

df_out_msm <- df_msm_wide %>%
  left_join(df_w_final, by = "encounter_block") %>%
  left_join(A_hist,     by = "encounter_block") %>%
  mutate(
    death_event = as.integer(outcome == 2),
    disch_event = as.integer(outcome == 1)
  )

## ---- B. Death cause-specific hazard ----
# Doubly robust, adding weight, lactate, and bibcarb due to residual confounding
fit_msm_cs_death <- coxph(
  Surv(time_to_event_90d, death_event) ~ high_count + 
    weight_kg + lactate_0 + bicarbonate_0,
  data = df_out_msm,
  weights = w_msm,
  robust = TRUE,
  cluster = encounter_block
)

## ---- C. Discharge cause-specific hazard ----
# Doubly robust, adding weight, lactate, and bibcarb due to residual confounding
fit_msm_cs_disch <- coxph(
  Surv(time_to_event_90d, disch_event) ~ high_count + 
    weight_kg + lactate_0 + bicarbonate_0,
  data = df_out_msm,
  weights = w_msm,
  robust = TRUE,
  cluster = encounter_block
)

summary(fit_msm_cs_death)
summary(fit_msm_cs_disch)

## ---- D. Extract MSM Cox results + save CSV ----

extract_cox_robust <- function(fit, model_label) {
  s <- summary(fit)
  co <- as.data.frame(s$coefficients)
  co$variable <- rownames(co)
  rownames(co) <- NULL

  ci <- as.data.frame(confint(fit))
  ci$variable <- rownames(ci)
  rownames(ci) <- NULL
  names(ci)[1:2] <- c("ci_lower", "ci_upper")

  # Detect the z, p-value, and SE column names (vary across survival versions)
  z_col <- intersect(c("robust z", "z"), names(co))[1]
  p_col <- intersect(c("Pr(>|z|)", "Robust Pr(>|z|)"), names(co))[1]
  se_col <- intersect(c("robust se", "se(coef)"), names(co))[1]

  out <- co %>%
    dplyr::left_join(ci, by = "variable") %>%
    dplyr::mutate(
      model = model_label,
      HR = exp(coef),
      HR_lower = exp(ci_lower),
      HR_upper = exp(ci_upper),
      se_log_hr = .data[[se_col]],
      z = .data[[z_col]],
      p_value = .data[[p_col]]
    ) %>%
    dplyr::select(model, variable, HR, HR_lower, HR_upper, se_log_hr, z, p_value)

  out
}

msm_results <- dplyr::bind_rows(
  extract_cox_robust(fit_msm_cs_death, "MSM IPTW Cause-Specific Cox - Death"),
  extract_cox_robust(fit_msm_cs_disch, "MSM IPTW Cause-Specific Cox - Discharge")
)

write.csv(
  msm_results,
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_MSM_IPTW_CauseSpecificCox_results.csv")),
  row.names = FALSE
)

cat("Saved MSM Cox results CSV.\n")

## ---- E. Table S1 ----
# Weighted cohort characteristics in the MSM pseudopopulation

# New collapse function without p-values
collapse_binary_no_p <- function(tbl, variable, keep_level, label) {

  modify_table_body(tbl, function(body) {

    # Identify rows for that variable
    idx <- which(body$variable == variable)
    sub <- body[idx, ]
    rest <- body[-idx, ]

    # Build 1-row replacement without expecting a p-value column
    sub2 <- sub %>%
      dplyr::filter(!(row_type == "label"),
                    !(row_type == "level" & label != keep_level)) %>%
      dplyr::mutate(
        row_type = "label",
        label = !!label
      )

    # Reinsert in place
    out <- dplyr::bind_rows(
      rest[seq_len(min(idx) - 1), ],
      sub2,
      rest[seq(min(idx), nrow(rest)), ]
    )

    out
  })
}

# Ensure df_out_msm contains: w_msm, high_count, A_0, A_12
stopifnot(all(c("w_msm","high_count") %in% names(df_out_msm)))

# Prepare MSM Table S1 dataframe (mirror Table 1 preprocessing)
df_msm_tableS1 <- df_out_msm %>%
  mutate(
    # Regime summary groups (4 directional groups matching Section 6)
    crrt_group = factor(
      case_when(
        A_0 == 0 & A_12 == 0 ~ "low_low",
        A_0 == 0 & A_12 == 1 ~ "low_high",
        A_0 == 1 & A_12 == 0 ~ "high_low",
        A_0 == 1 & A_12 == 1 ~ "high_high"
      ),
      levels = c("low_low", "low_high", "high_low", "high_high"),
      labels = c(
        paste0("Low -> Low (<", dose_cutoff, " both)"),
        paste0("Low -> High (escalated to >=", dose_cutoff, ")"),
        paste0("High -> Low (de-escalated from >=", dose_cutoff, ")"),
        paste0("High -> High (>=", dose_cutoff, " both)")
      )
    ),

    # Ensure sex is factor
    sex_category = factor(sex_category),

    # Clean race labels (same as Table 1)
    race_category = factor(race_category),

    # Clean CRRT modality labels (same as Table 1)
    crrt_mode_category = factor(crrt_mode_category),

    # Three-category outcome (same as Table 1)
    outcome_3cat = factor(
      case_when(
        outcome == 1 ~ "Discharged",
        outcome == 2 ~ "Died",
        TRUE         ~ "Censored"
      ),
      levels = c("Discharged", "Died", "Censored")
    )
  ) %>%
  # Keep only variables for the table
  select(
    crrt_group,
    age_at_admission, sex_category, weight_kg, race_category, ethnicity_category,
    sofa_total_0,
    oxygenation_index_0, norepinephrine_equivalent_0, imv_status_0,
    lactate_0, bicarbonate_0, potassium_0,
    crrt_mode_category,
    all_of(cci_vars),
    outcome_3cat,
    w_msm
  )

# Survey design using MSM IPTW weights
design_msm <- survey::svydesign(
  ids = ~1,
  weights = ~w_msm,
  data = df_msm_tableS1
)

# Build type list for Table S1
tableS1_type <- list(
  age_at_admission            ~ "continuous",
  sex_category                ~ "categorical",
  weight_kg                   ~ "continuous",
  race_category               ~ "categorical",
  ethnicity_category          ~ "categorical",
  sofa_total_0                ~ "continuous",
  oxygenation_index_0         ~ "continuous",
  norepinephrine_equivalent_0 ~ "continuous",
  imv_status_0                ~ "dichotomous",
  lactate_0                   ~ "continuous",
  bicarbonate_0               ~ "continuous",
  potassium_0                 ~ "continuous",
  crrt_mode_category          ~ "categorical",
  outcome_3cat                ~ "categorical"
)
for (v in cci_vars) {
  tableS1_type[[length(tableS1_type) + 1]] <- as.formula(
    paste0(v, ' ~ "dichotomous"'))
}

# Auto-detect dichotomous vars and set value=1; drop all-zero vars
tableS1_value <- list()
tableS1_drop <- character(0)
for (i in seq_along(tableS1_type)) {
  if (grepl("dichotomous", deparse(tableS1_type[[i]]))) {
    vname <- all.vars(tableS1_type[[i]])[1]
    if (vname %in% names(df_out_msm) && !any(df_out_msm[[vname]] == 1L, na.rm = TRUE)) {
      tableS1_drop <- c(tableS1_drop, vname)
    } else {
      tableS1_value[[length(tableS1_value) + 1]] <- as.formula(paste0(vname, " ~ 1L"))
    }
  }
}
if (length(tableS1_drop) > 0) {
  cat("  Dropping all-zero dichotomous vars from Table S1:", paste(tableS1_drop, collapse = ", "), "\n")
  tableS1_type <- tableS1_type[!sapply(tableS1_type, function(f) all.vars(f)[1] %in% tableS1_drop)]
}

# Build label list for Table S1
tableS1_label <- list(
  age_at_admission             ~ "Age at Admission (years)",
  sex_category                 ~ "Female (%)",
  weight_kg                    ~ "Weight (kg)",
  race_category                ~ "Race",
  ethnicity_category           ~ "Ethnicity",
  sofa_total_0                 ~ "SOFA Score at CRRT Start",
  oxygenation_index_0          ~ "Oxygenation Index (P/F or S/F)",
  norepinephrine_equivalent_0  ~ "NE Equivalent (mcg/kg/min)",
  imv_status_0                 ~ "On IMV (%)",
  lactate_0                    ~ "Lactate at CRRT Start (mmol/L)",
  bicarbonate_0                ~ "Bicarbonate at CRRT Start (mEq/L)",
  potassium_0                  ~ "Potassium at CRRT Start (mEq/L)",
  crrt_mode_category           ~ "CRRT Modality",
  outcome_3cat                 ~ "90-day Outcome"
)
for (v in cci_vars) {
  tableS1_label[[length(tableS1_label) + 1]] <- as.formula(
    paste0(v, ' ~ "', cci_labels[v], '"'))
}

# Also drop from label list
if (length(tableS1_drop) > 0) {
  tableS1_label <- tableS1_label[!sapply(tableS1_label, function(f) all.vars(f)[1] %in% tableS1_drop)]
}

# Build weighted Table S1
tableS1_msm <- tbl_svysummary(
  design_msm,
  by = crrt_group,
  type = tableS1_type,
  value = tableS1_value,
  statistic = list(
    all_continuous()  ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)"
  ),
  label = tableS1_label
) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(
    c("stat_1","stat_2","stat_3","stat_4")
    ~ "**CRRT Dose Strategy in the First 24h**") %>%
  modify_caption(
    "**Table S1. Weighted Cohort Characteristics After MSM IPTW (SuperLearner)**")

# Collapse Sex row
tableS1_msm <- tableS1_msm %>%
  collapse_binary_no_p(
    variable   = "sex_category",
    keep_level = "female",
    label      = "Female (%)"
  )

# Save Table S1
gt::gtsave(
  as_gt(tableS1_msm),
  filename = file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_TableS1_MSM_IPTW.html"))
)
write.csv(as_tibble(tableS1_msm),
          file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_TableS1_MSM_IPTW.csv")),
          row.names = FALSE)

cat("Saved MSM Table S1.\n")


# ============================================================= #
# ---- 6. MSM-STANDARDIZED COMPETING RISK CIF CURVES ----
# (model-based standardized CIFs from weighted cause-specific Cox)
# ============================================================= #

## ---- A. Define the MSM strategy grouping ----

### ---- I. Choose between Option 1, 2, or 3 ---- 

# ___________________________________________________________________
#  OPTION 1 - 3 groups - 'Always low' / 'Mixed' / 'Always high'
#df_plot_msm <- df_out_msm %>%
#  mutate(
#    crrt_strategy = factor(
#      high_count,
#      levels = c(0, 1, 2),
#      labels = c(
#        paste0("Always low (<", dose_cutoff, ")"),
#        paste0("Mixed (1/2 intervals high)"),
#        paste0("Always high (>=", dose_cutoff, ")")
#      )
#    )
#  )

# OPTION 2 - Two groups - 'Always low' / 'Always high'
#df_plot_msm <- df_out_msm %>%
#  filter(high_count %in% c(0, 2)) %>%
#  mutate(
#    crrt_strategy = factor(
#      high_count,
#      levels = c(0, 2),
#      labels = c(
#        paste0("Always low (<", dose_cutoff, ")"),
#        paste0("Always high (>=", dose_cutoff, ")")
#      )
#    )
#  )

# OPTION 3 - 4 groups - directional dose strategies
df_plot_msm <- df_out_msm %>%
  mutate(
    crrt_strategy = factor(
      case_when(
        A_0 == 0 & A_12 == 0 ~ "low_low",
        A_0 == 0 & A_12 == 1 ~ "low_high",
        A_0 == 1 & A_12 == 0 ~ "high_low",
        A_0 == 1 & A_12 == 1 ~ "high_high"
      ),
      levels = c("low_low", "low_high", "high_low", "high_high"),
      labels = c(
        paste0("Low -> Low (<", dose_cutoff, " both)"),
        paste0("Low -> High (escalated to >=", dose_cutoff, ")"),
        paste0("High -> Low (de-escalated from >=", dose_cutoff, ")"),
        paste0("High -> High (>=", dose_cutoff, " both)")
      )
    )
  )
# ___________________________________________________________________

### ---- II. Strategy group counts and outcomes ----                                                                                                                                
df_4group_summary <- df_plot_msm %>%
  group_by(crrt_strategy) %>%
  summarise(
    n_patients = n(),
    deaths_90d = sum(outcome == 2),
    discharged_90d = sum(outcome == 1),
    censored_90d = sum(outcome == 0),
    mortality_pct = round(mean(outcome == 2) * 100, 1),
    discharge_pct = round(mean(outcome == 1) * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(
    group_pct = round(n_patients / sum(n_patients) * 100, 1),
    dose_cutoff_mlkghr = dose_cutoff
  )

write.csv(
  df_4group_summary,
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                               "_MSM_4group_strategy_summary.csv")),
  row.names = FALSE
)
cat("\n4-group strategy summary:\n")
print(df_4group_summary)

## ---- B. Fit weighted cause-specific Cox models ----

# Death cause-specific hazard
fit_cs_death_msm <- coxph(
  Surv(time_to_event_90d, outcome == 2) ~ crrt_strategy + 
    weight_kg + lactate_0 + bicarbonate_0,
  data = df_plot_msm,
  weights = w_msm,
  robust = TRUE,
  cluster = encounter_block
)

# Discharge cause-specific hazard
fit_cs_disch_msm <- coxph(
  Surv(time_to_event_90d, outcome == 1) ~ crrt_strategy + 
    weight_kg + lactate_0 + bicarbonate_0,
  data = df_plot_msm,
  weights = w_msm,
  robust = TRUE,
  cluster = encounter_block
)

cif_cox_results <- dplyr::bind_rows(
  extract_cox_robust(fit_cs_death_msm, "MSM CIF Cause-Specific Cox - Death (4-group)"),
  extract_cox_robust(fit_cs_disch_msm, "MSM CIF Cause-Specific Cox - Discharge (4-group)")
)

write.csv(
  cif_cox_results,
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_MSM_CIF_CauseSpecificCox_4group_results.csv")),
  row.names = FALSE
)

cat("Saved CIF cause-specific Cox (4-group strategy) results CSV.\n")

## ---- C. Helpers to compute standardized CIFs from cause-specific hazards ----

# Carry-forward cumulative baseline hazard to a common grid
cumhaz_at_grid <- function(basehaz_df, grid_time) {
  idx <- findInterval(grid_time, basehaz_df$time)
  out <- numeric(length(grid_time))
  out[idx > 0] <- basehaz_df$hazard[idx[idx > 0]]
  out
}

# Aalen-Johansen recursion for 2 competing causes:
# Cause 1 = death, Cause 2 = discharge
compute_cifs_two_causes <- function(mult_death, mult_disch, tgrid, dH0_death, dH0_disch) {
  dH_death <- dH0_death * mult_death
  dH_disch <- dH0_disch * mult_disch

  S <- numeric(length(tgrid))
  CIF_death <- numeric(length(tgrid))
  CIF_disch <- numeric(length(tgrid))

  S_prev <- 1
  cd_prev <- 0
  cc_prev <- 0

  for (i in seq_along(tgrid)) {
    cd_prev <- cd_prev + S_prev * dH_death[i]
    cc_prev <- cc_prev + S_prev * dH_disch[i]
    S_prev  <- S_prev * exp(-(dH_death[i] + dH_disch[i]))

    CIF_death[i] <- cd_prev
    CIF_disch[i] <- cc_prev
    S[i] <- S_prev
  }

  tibble::tibble(time = tgrid, cif_death = CIF_death, cif_disch = CIF_disch, surv = S)
}

# Build CIF curves for each strategy level
build_standardized_cifs <- function(fit_death, fit_disch, strategy_levels, newdata) {
  
  # Baseline cumulative hazards
  bh_death <- basehaz(fit_death, centered = FALSE)
  bh_disch <- basehaz(fit_disch, centered = FALSE)
  
  # Common time grid
  tgrid <- sort(unique(c(bh_death$time, bh_disch$time)))
  tgrid <- tgrid[tgrid >= 0]
  
  # Baseline cumhaz on grid + increments
  H0_death <- cumhaz_at_grid(bh_death, tgrid)
  H0_disch <- cumhaz_at_grid(bh_disch, tgrid)
  dH0_death <- c(H0_death[1], diff(H0_death))
  dH0_disch <- c(H0_disch[1], diff(H0_disch))
  
  nt <- length(tgrid)
  
  # For each strategy level, compute population-averaged CIFs (g-computation)
  out <- lapply(strategy_levels, function(lv) {
    
    # Counterfactual: assign everyone to this strategy, keep real covariates
    cf_data <- newdata
    cf_data$crrt_strategy <- factor(lv, levels = strategy_levels)
    
    # Individual linear predictors under this counterfactual strategy
    lp_death <- predict(fit_death, newdata = cf_data, type = "lp")
    lp_disch <- predict(fit_disch, newdata = cf_data, type = "lp")
    
    mult_death <- exp(lp_death)  # n-vector
    mult_disch <- exp(lp_disch)  # n-vector
    n <- length(mult_death)
    
    # Vectorized Aalen-Johansen across patients, sequential over time
    S_prev <- rep(1, n)
    cif_death_accum <- rep(0, n)
    cif_disch_accum <- rep(0, n)
    
    avg_cif_death <- numeric(nt)
    avg_cif_disch <- numeric(nt)
    avg_surv      <- numeric(nt)
    
    for (j in seq_len(nt)) {
      dH_death_j <- dH0_death[j] * mult_death
      dH_disch_j <- dH0_disch[j] * mult_disch
      
      cif_death_accum <- cif_death_accum + S_prev * dH_death_j
      cif_disch_accum <- cif_disch_accum + S_prev * dH_disch_j
      S_prev <- S_prev * exp(-(dH_death_j + dH_disch_j))
      
      avg_cif_death[j] <- mean(cif_death_accum)
      avg_cif_disch[j] <- mean(cif_disch_accum)
      avg_surv[j]      <- mean(S_prev)
    }
    
    tibble::tibble(
      time = tgrid,
      cif_death = avg_cif_death,
      cif_disch = avg_cif_disch,
      surv = avg_surv,
      crrt_strategy = lv
    )
  }) %>%
    dplyr::bind_rows() %>%
    mutate(crrt_strategy = factor(crrt_strategy, levels = strategy_levels))
  
  out
}

## ---- D. Generate CIF dataframe ----
strategy_levels <- levels(df_plot_msm$crrt_strategy)

cif_df <- build_standardized_cifs(
  fit_death = fit_cs_death_msm,
  fit_disch = fit_cs_disch_msm,
  strategy_levels = strategy_levels,
  newdata = df_plot_msm
)

## ---- E. Bootstrap CIs for standardized CIFs ----
N_BOOT <- 200  # Adjustable
set.seed(42)

cat("Bootstrapping CIF confidence intervals (", N_BOOT, "replicates) …\n")

# Common time grid from the point estimate
tgrid_common <- sort(unique(cif_df$time))

# Storage: list of CIF matrices per strategy
boot_death <- list()
boot_disch <- list()
for (lv in strategy_levels) {
  boot_death[[lv]] <- matrix(NA, nrow = N_BOOT, ncol = length(tgrid_common))
  boot_disch[[lv]] <- matrix(NA, nrow = N_BOOT, ncol = length(tgrid_common))
}

unique_ids <- unique(df_plot_msm$encounter_block)
n_ids <- length(unique_ids)

for (b in seq_len(N_BOOT)) {
  # Resample encounter_blocks with replacement
  boot_ids <- sample(unique_ids, n_ids, replace = TRUE)
  
  # Build bootstrap dataset (handles duplicate IDs)
  boot_df <- data.frame(encounter_block = boot_ids, .boot_idx = seq_along(boot_ids)) %>%
    left_join(df_plot_msm, by = "encounter_block", relationship = "many-to-many")
  
  # Refit cause-specific Cox models
  fit_b <- tryCatch({
    fit_d <- coxph(
      Surv(time_to_event_90d, outcome == 2) ~ crrt_strategy +
        weight_kg + lactate_0 + bicarbonate_0,
      data = boot_df, weights = w_msm, robust = FALSE
    )
    fit_c <- coxph(
      Surv(time_to_event_90d, outcome == 1) ~ crrt_strategy +
        weight_kg + lactate_0 + bicarbonate_0,
      data = boot_df, weights = w_msm, robust = FALSE
    )
    list(death = fit_d, disch = fit_c)
  }, error = function(e) NULL)
  
  if (is.null(fit_b)) next
  
  # Compute CIFs on the common grid
  boot_cif <- tryCatch(
    build_standardized_cifs(fit_b$death, fit_b$disch, strategy_levels, 
                            newdata = boot_df),
    error = function(e) NULL
  )
  
  if (is.null(boot_cif)) next
  
  # Interpolate to common grid and store
  for (lv in strategy_levels) {
    lv_df <- boot_cif %>% filter(crrt_strategy == lv)
    if (nrow(lv_df) == 0) next
    boot_death[[lv]][b, ] <- approx(lv_df$time, lv_df$cif_death,
                                    xout = tgrid_common, rule = 2)$y
    boot_disch[[lv]][b, ] <- approx(lv_df$time, lv_df$cif_disch,
                                    xout = tgrid_common, rule = 2)$y
  }
  
  if (b %% 50 == 0) cat("  ", b, "/", N_BOOT, "\n")
}

# Compute pointwise 2.5% and 97.5% percentiles
ci_list <- lapply(strategy_levels, function(lv) {
  tibble(
    time = tgrid_common,
    crrt_strategy = lv,
    cif_death_lower = apply(boot_death[[lv]], 2, quantile, 0.025, na.rm = TRUE),
    cif_death_upper = apply(boot_death[[lv]], 2, quantile, 0.975, na.rm = TRUE),
    cif_disch_lower = apply(boot_disch[[lv]], 2, quantile, 0.025, na.rm = TRUE),
    cif_disch_upper = apply(boot_disch[[lv]], 2, quantile, 0.975, na.rm = TRUE)
  )
})
ci_df <- bind_rows(ci_list) %>%
  mutate(crrt_strategy = factor(crrt_strategy, levels = strategy_levels))

# Merge CIs into the point estimate dataframe
cif_df <- cif_df %>%
  left_join(ci_df, by = c("time", "crrt_strategy"))

cat("Bootstrap complete.\n")

## ---- F. Plot CIFs (Death + Discharge) ----

### ---- I. Death CIF ----
plot_cif_death <- ggplot(cif_df, aes(x = time, y = cif_death, color = crrt_strategy,
                                     fill = crrt_strategy)) +
  geom_ribbon(aes(ymin = cif_death_lower, ymax = cif_death_upper),
              alpha = 0.2, linewidth = 0) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 7)) +
  guides(color = guide_legend(nrow = 2),
         fill = guide_legend(nrow = 2)) +
  labs(
    title = "MSM Standardized CIF: Death",
    x = "Time from CRRT Initiation (Days)",
    y = "Cumulative Incidence (Death)",
    color = "CRRT Strategy (0-24h)",
    fill = "CRRT Strategy (0-24h)"
  )

plot_cif_death

ggsave(
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_MSM_CIF_Death.png")),
  plot_cif_death, width = 6, height = 4, dpi = 300
)

### ---- II. Discharge CIF ----
plot_cif_disch <- ggplot(cif_df, aes(x = time, y = cif_disch, color = crrt_strategy,
                                     fill = crrt_strategy)) +
  geom_ribbon(aes(ymin = cif_disch_lower, ymax = cif_disch_upper),
              alpha = 0.2, linewidth = 0) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 7)) +
  guides(color = guide_legend(nrow = 2),
         fill = guide_legend(nrow = 2)) +
  labs(
    title = "MSM Standardized CIF: Discharge",
    x = "Time from CRRT Initiation (Days)",
    y = "Cumulative Incidence (Discharge)",
    color = "CRRT Strategy (0-24h)",
    fill = "CRRT Strategy (0-24h)"
  )

plot_cif_disch

ggsave(
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_MSM_CIF_Discharge.png")),
  plot_cif_disch, width = 6, height = 4, dpi = 300
)

# Export MSM CIF data as CSV
write.csv(cif_df,
          file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_MSM_CIF_data.csv")),
          row.names = FALSE)
cat("MSM standardized CIF plots with CIs saved as PNGs + CSV.\n")

# ============================================================= #
# ---- 8. FINISH! ----
# ============================================================= #
cat("Code run complete!\n")
