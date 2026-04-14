####### Propensity Score Matching (PSM) Model #######################
####### Inverse Probability of Treatment Weighting (IPTW) Model #####

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
# Works both interactively (RStudio) and via Rscript
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  # Running in RStudio - use active document path
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  root_dir <- dirname(script_dir)
  setwd(root_dir)
} else {
  # Running via Rscript or not in RStudio
  # Try to get script location
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    script_dir <- dirname(script_path)
    root_dir <- dirname(script_dir)
    setwd(root_dir)
  } else {
    # Assume already in project root or code directory
    if (basename(getwd()) == "code") {
      setwd("..")
    }
    # If current directory has 'code' subdirectory, we're in root
    if (!dir.exists("code")) {
      stop("Please run this script from the project root directory or use Rscript")
    }
  }
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
output_dir <- "output/final/psm_iptw"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

## ---- D. Rubin's rules helper for manual pooling ----
# Used for models that mice::pool() doesn't support (crr, weighted coxph)
rubins_pool <- function(betas, ses) {
  # betas: vector of m coefficient estimates
  # ses:   vector of m standard errors
  # Returns: list(beta, se, hr, hr_lower, hr_upper, p_value, df)
  m <- length(betas)
  beta_bar <- mean(betas)
  W <- mean(ses^2)               # Within-imputation variance
  B <- var(betas)                 # Between-imputation variance
  T_var <- W + (1 + 1/m) * B     # Total variance
  se_pooled <- sqrt(T_var)

  # Barnard-Rubin degrees of freedom
  r <- (1 + 1/m) * B / W
  df_old <- (m - 1) * (1 + 1/r)^2
  # Simplified: use df_old (sufficient for m=5)
  df <- df_old

  t_stat <- beta_bar / se_pooled
  p_value <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)

  list(
    beta     = beta_bar,
    se       = se_pooled,
    hr       = exp(beta_bar),
    hr_lower = exp(beta_bar - 1.96 * se_pooled),
    hr_upper = exp(beta_bar + 1.96 * se_pooled),
    p_value  = p_value,
    df       = df,
    fmi      = (B + B/m) / T_var  # Fraction of missing information
  )
}

## ---- E. Load configuration ----
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
  "time_to_event_30d", "outcome",
  "age_at_admission", "sex_category", "race_category", "ethnicity_category",
  "crrt_mode_category", "weight_kg",
  "crrt_dose_ml_kg_hr_0",
  # Labs at t=0
  "lactate_0", "bicarbonate_0", "potassium_0",
  # SOFA (kept in data for descriptive use; NOT used in models)
  "sofa_total_0",
  # New covariates at t=0
  "oxygenation_index_0", "norepinephrine_equivalent_0", "imv_status_0",
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
  "White" = "white",
  "Black" = "black or african american",
  other_level = "Other"
)

# ================================ #
# ---- 1. DATA CLEANING ----
# ================================ #

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Generating Sample Characteristics\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

## ---- A. Define model variables for complete case analysis ----
model_vars <- required_vars

# Report missingness before imputation
cat("Missing values per variable:\n")
miss_counts <- colSums(is.na(df[, model_vars]))
print(miss_counts[miss_counts > 0])
cat("\n")

# MICE imputation — only on variables with missing values, using nearby predictors
n_incomplete <- sum(!complete.cases(df[, model_vars]))
if (n_incomplete > 0) {
  cat("Imputing", n_incomplete, "incomplete rows via MICE (pmm, m=5) …\n")
  library(mice)

  # Only impute vars that actually have missingness, plus a small set of predictors
  vars_with_na <- names(miss_counts[miss_counts > 0])
  predictor_vars <- c("age_at_admission", "sex_category", "lactate_0",
                       "bicarbonate_0", "potassium_0", "norepinephrine_equivalent_0",
                       "imv_status_0", "crrt_dose_ml_kg_hr_0")
  mice_vars <- unique(c(vars_with_na, intersect(predictor_vars, model_vars)))
  mice_vars <- mice_vars[mice_vars %in% names(df)]

  imp <- mice(df[, mice_vars], m = 5, method = "pmm", seed = 42, maxit = 10, printFlag = FALSE)

  # Retain all 5 imputed datasets for Rubin's rules pooling
  N_IMP <- 5
  imp_list <- lapply(seq_len(N_IMP), function(m) {
    d <- df
    d[, mice_vars] <- complete(imp, m)
    d
  })
  df_complete <- imp_list[[1]]  # Primary dataset for PSM, SL, bootstrap, etc.
  cat("MICE imputed", sum(miss_counts[miss_counts > 0]),
      "total missing values (", N_IMP, "imputations retained)\n")
} else {
  cat("No missing values — skipping MICE\n")
  df_complete <- df
  N_IMP <- 1
  imp_list <- list(df_complete)
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
imp_list <- lapply(imp_list, sanitize_numeric)
cat("NaN/Inf check complete.\n\n")

cat("Analysis sample:", nrow(df_complete), "of", nrow(df), "\n\n")

# Check for sufficient data
if (nrow(df_complete) < 50) {
  stop("Insufficient cases for modeling (n = ", nrow(df_complete), ")")
}

## ---- B. Build model covariate set and formulas dynamically ----

# All potential covariates for propensity / outcome models
# NOTE: sofa_total_0 deliberately excluded from models
model_covariates <- c(
  "age_at_admission", "sex_category", "race_category", "weight_kg",
  "lactate_0", "bicarbonate_0", "potassium_0",
  "oxygenation_index_0", "norepinephrine_equivalent_0", "imv_status_0",
  cci_vars
)

# Dynamic filter: drop covariates with <2 unique values at this site
# (e.g., CCI components with zero prevalence, single-level categoricals)
model_covariates <- model_covariates[
  sapply(model_covariates, function(v) {
    length(unique(df_complete[[v]])) >= 2
  })
]

cat("Model covariates after dynamic filter (", length(model_covariates), "):\n")
cat(" ", paste(model_covariates, collapse = ", "), "\n\n")

# Build formulas programmatically (reused across PSM, IPTW, Cox)
psm_formula <- as.formula(
  paste("crrt_high ~", paste(model_covariates, collapse = " + "))
)

cat("PSM/IPTW formula:\n")
print(psm_formula)
cat("\n")

# Complete-case variable list: all model covariates + outcome/time variables
# Used by all drop_na() calls to ensure data is complete for every formula variable
complete_case_vars <- c(model_covariates, "sofa_total_0",
                        "crrt_dose_ml_kg_hr_0", "time_to_event_30d", "outcome")

## ---- C. CRRT Dose Cutoff ----
# =================================== #

# CRRT dose cutoff (mL/kg/hr)
dose_cutoff <- 30


cat("Using CRRT dose cutoff of:", dose_cutoff, "mL/kg/hr\n\n")
# =================================== #

## ---- D. Calculate summary statistics ----
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

  # SOFA score (descriptive only)
  sofa_mean = mean(df_complete$sofa_total_0, na.rm = TRUE),
  sofa_sd = sd(df_complete$sofa_total_0, na.rm = TRUE),
  sofa_median = median(df_complete$sofa_total_0, na.rm = TRUE),
  sofa_q25 = quantile(df_complete$sofa_total_0, 0.25, na.rm = TRUE),
  sofa_q75 = quantile(df_complete$sofa_total_0, 0.75, na.rm = TRUE),

  # CRRT dose (modality independent, time t=0)
  crrt_dose_mean = mean(df_complete$crrt_dose_ml_kg_hr_0, na.rm = TRUE),
  crrt_dose_sd = sd(df_complete$crrt_dose_ml_kg_hr_0, na.rm = TRUE),
  crrt_dose_median = median(df_complete$crrt_dose_ml_kg_hr_0, na.rm = TRUE),
  crrt_dose_q25 = quantile(df_complete$crrt_dose_ml_kg_hr_0, 0.25, na.rm = TRUE),
  crrt_dose_q75 = quantile(df_complete$crrt_dose_ml_kg_hr_0, 0.75, na.rm = TRUE),

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
          file.path(output_dir, paste0(SITE_NAME, "_sample_characteristics.csv")),
          row.names = FALSE)

cat("Sample characteristics saved.\n\n")

## ---- E. Histogram CRRT dose ----
crrt_dose_hist <- ggplot(df_complete, aes(x = crrt_dose_ml_kg_hr_0)) +
  geom_histogram(binwidth = 5, color = "black", fill = "skyblue") +
  geom_vline(xintercept = dose_cutoff, linetype = "dashed", linewidth = 1) +
      # Line at cutoff
  labs(
    title = paste0("Distribution of Initial CRRT Dose (Cutoff = 
                   ", dose_cutoff, " mL/kg/hr)"),
    x = "CRRT Dose (mL/kg/hr)",
    y = "Count"
  ) +
  theme_bw(base_size = 12)
crrt_dose_hist

# Save histogram
ggsave(file.path(output_dir, paste0(SITE_NAME, "_hist_crrt_dose.png")),
       crrt_dose_hist, width=6, height=4)
cat("Saved histogram PNG\n")

## ---- F. Designate High vs Low CRRT Dose ----
cat("Defining treatment groups using cutoff =", dose_cutoff, "mL/kg/hr\n")

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
                              paste0(SITE_NAME, "_crrt_bin_distribution.csv")),
          row.names = FALSE)
cat("Saved bin distribution\n")


## ---- G. TABLE 1 ----
# =================================== #

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
    crrt_mode_category = factor(toupper(as.character(crrt_mode_category))),

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
  outcome_3cat                 ~ "30-day Outcome"
)
for (v in cci_vars) {
  table1_label[[length(table1_label) + 1]] <- as.formula(
    paste0(v, ' ~ "', cci_labels[v], '"'))
}

# Auto-detect all dichotomous variables and set value=1 for gtsummary
table1_value <- list()
for (i in seq_along(table1_type)) {
  if (grepl("dichotomous", deparse(table1_type[[i]]))) {
    vname <- all.vars(table1_type[[i]])[1]
    table1_value[[length(table1_value) + 1]] <- as.formula(paste0(vname, " ~ 1L"))
  }
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
    # add_p can fail when categorical variables have <2 levels — fall back gracefully
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
    "**Table 1. Baseline Characteristics of the Cohort by CRRT Dose Group**")

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
                       paste0(SITE_NAME, "_Table1_unadjusted.html"))
)
write.csv(as_tibble(table1),
          file.path(output_dir, paste0(SITE_NAME, "_Table1_unadjusted.csv")),
          row.names = FALSE)

# ============================================ #
# ---- 2. PROPENSITY SCORE MATCHING BRANCH ----
# ============================================ #

## ---- A. Making PSM dataset ----
### ---- I. PSM ----
m.out <- matchit(
  psm_formula,
  data = df_tte_bin,
  method = "nearest",
  distance = "logit",
  ratio = 1,
  caliper = 0.2
)

# Balance Summary
summary(m.out)
# Save PSM summary (both balance and counts)
psm_bal <- summary(m.out)$sum.matched
psm_counts <- summary(m.out)$nn %>%
  as.data.frame() %>%
  rownames_to_column(var = "Category")
write.csv(psm_bal, file.path(output_dir, paste0(SITE_NAME, "_psm_balance_summary.csv")),
          row.names = FALSE)
write.csv(psm_counts, file.path(output_dir, paste0(SITE_NAME, "_psm_counts_summary.csv")),
          row.names = FALSE)
cat("Saved PSM summaries\n")

### ---- II. Save Matched Dataset ----
df_match <- match.data(m.out)

# Pretty variable names for Love plot
pretty_names_loveplot <- c(
  age_at_admission              = "Age",
  sex_category                  = "Sex",
  race_category                 = "Race",
  weight_kg                     = "Weight (kg)",
  lactate_0                     = "Lactate",
  bicarbonate_0                 = "Bicarbonate",
  potassium_0                   = "Potassium",
  oxygenation_index_0           = "Oxygenation Index",
  norepinephrine_equivalent_0   = "NE Equivalent",
  imv_status_0                  = "IMV Status",

  # For IPTW SL (ps column)
  prop.score                    = "Propensity Score",

  # For PSM (MatchIt internal name)
  distance                      = "Propensity Score",

  # CCI components
  cci_labels  # named vector appended here
)

## ---- B. Love Plot (SMD + variance ratio) ----
plot_loveplot_psm <- love.plot(
  m.out,
  stat = c("m", "v"),
  grid = TRUE,
  var.names = pretty_names_loveplot,
  title = "Covariate Balance: PSM",
  threshold = c(m = .25, v = 1.25)
)

print(plot_loveplot_psm)

png(file.path(output_dir, paste0(SITE_NAME, "_psm_loveplot.png")),
    width = 8, height = 7, units = "in", res = 300)
print(plot_loveplot_psm)
dev.off()

cat("Saved Love plot for PSM as PNG\n")

## ---- C. TABLE S1 (PSM) ----

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

# Prepare matched data (mirror the Table 1 preprocessing)
df_tte_tableS1 <- df_match %>%
  mutate(
    # ---- Rename CRRT groups using dose_cutoff
    crrt_group = factor(
      crrt_high,
      levels = c(paste0("<", dose_cutoff),
                 paste0(">=", dose_cutoff)),
      labels = c(
        paste0("Low CRRT dose (<", dose_cutoff, " mL/kg/hr)"),
        paste0("High CRRT dose (>=", dose_cutoff, " mL/kg/hr)")
      )
    ),

    # ---- Ensure sex is a factor
    sex_category = factor(sex_category),

    # ---- Clean race labels
    race_category = factor(race_category),

    # ---- Clean CRRT modality labels
    crrt_mode_category = forcats::fct_recode(
      crrt_mode_category,
      "CVVH"   = "cvvh",
      "CVVHD"  = "cvvhd",
      "CVVHDF" = "cvvhdf",
      "SCUF"   = "scuf",
      "AVVH"   = "avvh"
    ),

    # ---- Three-category outcome
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


# Build TABLE S1 (same specification as Table 1)
tableS1 <- df_tte_tableS1 %>%
  select(crrt_group, all_of(vars_table1)) %>%
  tbl_summary(
    by = crrt_group,
    type = table1_type,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = table1_label
  ) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(
    c("stat_1","stat_2") ~ "**CRRT dose group**"
  ) %>%
  modify_caption(
    "**Table S1.
    Baseline Characteristics of the Propensity Score-Matched Cohort**"
  )

# Collapse Sex row
tableS1 <- tableS1 %>%
  collapse_binary_no_p(
    variable   = "sex_category",
    keep_level = "female",
    label      = "Female (%)"
  )

# Save Table S1 as HTML
gt::gtsave(
  as_gt(tableS1),
  filename = file.path(
    output_dir,
    paste0(SITE_NAME, "_TableS1_matched.html")
  )
)
write.csv(as_tibble(tableS1),
          file.path(output_dir, paste0(SITE_NAME, "_TableS1_matched.csv")),
          row.names = FALSE)

## ---- D. Analysis of PSM ----
### ---- I. Fine-Gray Analysis on Matched Patients ----

# Matrix for treatment only (reference = low dose)
X_trt <- model.matrix(~ crrt_high, df_match)[, -1, drop = FALSE]

fg_death_psm <- cmprsk::crr(
  ftime   = df_match$time_to_event_30d,
  fstatus = df_match$outcome,
  cov1    = X_trt,
  failcode = 2,
  cencode  = 0
)

fg_disch_psm <- cmprsk::crr(
  ftime   = df_match$time_to_event_30d,
  fstatus = df_match$outcome,
  cov1    = X_trt,
  failcode = 1,
  cencode  = 0
)

summary(fg_death_psm)
summary(fg_disch_psm)

# Add doubly robust model (treatment + covariates)
dr_formula <- as.formula(
  paste("~", paste(c("crrt_high", model_covariates), collapse = " + "))
)
X_dr <- model.matrix(dr_formula, data = df_match)[, -1, drop = FALSE]
colnames(X_dr) <- make.names(colnames(X_dr))  # SL.gam compatibility

fg_death_psm_dr <- cmprsk::crr(
  ftime   = df_match$time_to_event_30d,
  fstatus = df_match$outcome,
  cov1    = X_dr,
  failcode = 2,
  cencode  = 0
)

fg_disch_psm_dr <- cmprsk::crr(
  ftime   = df_match$time_to_event_30d,
  fstatus = df_match$outcome,
  cov1    = X_dr,
  failcode = 1,
  cencode  = 0
)
summary(fg_death_psm_dr)
summary(fg_disch_psm_dr)

# Save doubly robust Fine Gray summary
extract_fg <- function(fg_model, outcome_label) {
  s <- summary(fg_model)$coef
  data.frame(
    outcome = outcome_label,
    variable = rownames(s),
    coef = s[, "coef"],
    SHR = s[, "exp(coef)"],
    se = s[, "se(coef)"],
    p_value = s[, "p-value"],
    SHR_lower = s[, "exp(coef)"] / exp(1.96*s[,"se(coef)"]),
    SHR_upper = s[, "exp(coef)"] * exp(1.96*s[,"se(coef)"]),
    stringsAsFactors = FALSE
  )
}

fg_results <- bind_rows(
  extract_fg(fg_death_psm_dr, "Death"),
  extract_fg(fg_disch_psm_dr, "Discharge")
)

write.csv(fg_results, file.path(output_dir, paste0(SITE_NAME, "_fg_psm_dr_results.csv")),
          row.names = FALSE)
cat("Saved Fine Gray DR model results\n")

### ---- I-b. Rubin's Rules Pooled Fine-Gray Models ----

cat("\nPooling Fine-Gray DR models across", N_IMP, "MICE imputations...\n")

# Use the same PSM matched set (from imputation 1) but substitute imputed
# covariate values from each imputation. The matched set (row membership)
# stays fixed — only the imputed values change.
matched_ids <- as.integer(rownames(df_match))

pool_fg_treatment <- function(failcode_val, outcome_label) {
  betas <- numeric(N_IMP)
  ses   <- numeric(N_IMP)

  for (m in seq_len(N_IMP)) {
    d_imp <- imp_list[[m]]
    # Apply same data prep as df_tte_bin
    d_imp$crrt_high <- ifelse(d_imp$crrt_dose_ml_kg_hr_0 >= dose_cutoff, 1L, 0L)
    d_imp$crrt_high <- factor(d_imp$crrt_high, levels = c(0, 1),
                              labels = c(paste0("<", dose_cutoff),
                                         paste0(">=", dose_cutoff)))
    d_imp <- tidyr::drop_na(d_imp, dplyr::all_of(complete_case_vars))
    d_imp <- droplevels(d_imp)
    # Subset to matched rows
    d_match_m <- d_imp[rownames(d_imp) %in% as.character(matched_ids), ]

    X_dr_m <- model.matrix(dr_formula, data = d_match_m)[, -1, drop = FALSE]
    colnames(X_dr_m) <- make.names(colnames(X_dr_m))

    fg_m <- tryCatch(
      cmprsk::crr(
        ftime    = d_match_m$time_to_event_30d,
        fstatus  = d_match_m$outcome,
        cov1     = X_dr_m,
        failcode = failcode_val,
        cencode  = 0
      ),
      error = function(e) NULL
    )
    if (is.null(fg_m)) next

    s <- summary(fg_m)$coef
    trt_row <- grep("crrt_high", rownames(s))[1]
    betas[m] <- s[trt_row, "coef"]
    ses[m]   <- s[trt_row, "se(coef)"]
  }

  pooled <- rubins_pool(betas, ses)
  data.frame(
    model    = paste0("PSM FG (pooled) - ", outcome_label),
    HR_type  = "SHR",
    HR       = pooled$hr,
    HR_lower = pooled$hr_lower,
    HR_upper = pooled$hr_upper,
    se_log_hr = pooled$se,
    p_value  = pooled$p_value,
    fmi      = pooled$fmi,
    stringsAsFactors = FALSE
  )
}

pooled_fg_death <- pool_fg_treatment(2, "Death")
pooled_fg_disch <- pool_fg_treatment(1, "Discharge")
pooled_fg_results <- bind_rows(pooled_fg_death, pooled_fg_disch)

cat("Pooled Fine-Gray results:\n")
for (i in seq_len(nrow(pooled_fg_results))) {
  r <- pooled_fg_results[i, ]
  cat(sprintf("  %-35s SHR=%.2f (%.2f-%.2f) p=%.4f FMI=%.3f\n",
              r$model, r$HR, r$HR_lower, r$HR_upper, r$p_value, r$fmi))
}
cat("\n")

write.csv(pooled_fg_results,
          file.path(output_dir, paste0(SITE_NAME, "_fg_psm_pooled_results.csv")),
          row.names = FALSE)
cat("Saved pooled Fine-Gray results CSV\n")

### ---- II. CIF Curves ----

# Helper: tidy a cuminc object into a data.frame          
tidy_cuminc <- function(ci_obj) {                                                                                                                                     
  bind_rows(lapply(names(ci_obj), function(name) {                                                                                                                    
    if (!grepl(" ", name)) return(NULL)
    data.frame(
      group   = sub(" .*", "", name),
      outcome = sub(".* ", "", name),
      time    = ci_obj[[name]]$time,
      est     = ci_obj[[name]]$est
    )
  }))
}

# Point-estimate CIF
ci <- cuminc(
  ftime   = df_match$time_to_event_30d,
  fstatus = df_match$outcome,
  group   = df_match$crrt_high,
  cencode = 0
)
tidy_ci <- tidy_cuminc(ci)

# Common time grid for interpolation
tgrid_psm <- sort(unique(tidy_ci$time))
grp_levels <- levels(df_match$crrt_high)

# Bootstrap CIs (200 replicates)
N_BOOT_PSM <- 500
set.seed(42)
cat("Bootstrapping PSM CIF confidence intervals (", N_BOOT_PSM, "replicates) …\n")

boot_death_psm <- boot_disch_psm <- list()
for (lv in grp_levels) {
  boot_death_psm[[lv]] <- matrix(NA, nrow = N_BOOT_PSM, ncol = length(tgrid_psm))
  boot_disch_psm[[lv]] <- matrix(NA, nrow = N_BOOT_PSM, ncol = length(tgrid_psm))
}

for (b in seq_len(N_BOOT_PSM)) {
  boot_idx <- sample(seq_len(nrow(df_match)), nrow(df_match), replace = TRUE)
  boot_df  <- df_match[boot_idx, ]
  
  boot_ci <- tryCatch(
    cuminc(ftime = boot_df$time_to_event_30d,
           fstatus = boot_df$outcome,
           group = boot_df$crrt_high,
           cencode = 0),
    error = function(e) NULL
  )
  if (is.null(boot_ci)) next
  
  boot_tidy <- tidy_cuminc(boot_ci)
  
  for (lv in grp_levels) {
    # Death (outcome == 2)
    lv_death <- boot_tidy %>% filter(group == lv, outcome == "2")
    if (nrow(lv_death) > 0) {
      boot_death_psm[[lv]][b, ] <- approx(lv_death$time, lv_death$est,
                                          xout = tgrid_psm, method = "constant",
                                          rule = 2, f = 0)$y
    }
    # Discharge (outcome == 1)
    lv_disch <- boot_tidy %>% filter(group == lv, outcome == "1")
    if (nrow(lv_disch) > 0) {
      boot_disch_psm[[lv]][b, ] <- approx(lv_disch$time, lv_disch$est,
                                          xout = tgrid_psm, method = "constant",
                                          rule = 2, f = 0)$y
    }
  }
  if (b %% 100 == 0) cat("  ", b, "/", N_BOOT_PSM, "\n")
}

# Compute pointwise 95% CIs and merge with point estimates
ci_psm_list <- lapply(grp_levels, function(lv) {
  tibble(
    group = lv, time = tgrid_psm,
    ci_death_lower = apply(boot_death_psm[[lv]], 2, quantile, 0.025, na.rm = TRUE),
    ci_death_upper = apply(boot_death_psm[[lv]], 2, quantile, 0.975, na.rm = TRUE),
    ci_disch_lower = apply(boot_disch_psm[[lv]], 2, quantile, 0.025, na.rm = TRUE),
    ci_disch_upper = apply(boot_disch_psm[[lv]], 2, quantile, 0.975, na.rm = TRUE)
  )
})
ci_psm_df <- bind_rows(ci_psm_list)

# Join CIs onto point-estimate data
tidy_ci_death <- tidy_ci %>%
  filter(outcome == "2") %>%
  left_join(ci_psm_df %>% select(group, time, ci_death_lower, ci_death_upper),
            by = c("group", "time"))

tidy_ci_discharge <- tidy_ci %>%
  filter(outcome == "1") %>%
  left_join(ci_psm_df %>% select(group, time, ci_disch_lower, ci_disch_upper),
            by = c("group", "time"))

cat("PSM bootstrap complete.\n")

# Plot CIF for death
plot_cif_death <- ggplot(tidy_ci_death, aes(x = time, y = est,
                                            color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_death_lower, ymax = ci_death_upper),
              alpha = 0.2, linewidth = 0) +
  geom_step(linewidth = 1) +
  labs(
    title = "PSM Cumulative Incidence: Death",
    x = "Time from CRRT Initiation (Days)",
    y = "Cumulative Incidence (Death)",
    color = "CRRT Dose", fill = "CRRT Dose"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
plot_cif_death

ggsave(file.path(output_dir, paste0(SITE_NAME, "_cif_death.png")),
       plot_cif_death, width = 6, height = 4, dpi = 300)
cat("Saved CIF for death as PNG\n")

# Export PSM CIF data as CSV (death + discharge combined)
psm_cif_export <- bind_rows(
  tidy_ci_death %>% mutate(cif_lower = ci_death_lower, cif_upper = ci_death_upper) %>%
    select(group, outcome, time, est, cif_lower, cif_upper),
  tidy_ci_discharge %>% mutate(cif_lower = ci_disch_lower, cif_upper = ci_disch_upper) %>%
    select(group, outcome, time, est, cif_lower, cif_upper)
) %>% rename(cif = est)
write.csv(psm_cif_export,
          file.path(output_dir, paste0(SITE_NAME, "_PSM_CIF_data.csv")),
          row.names = FALSE)
cat("Saved PSM CIF data as CSV\n")

# Plot CIF for discharge
plot_cif_discharge <- ggplot(tidy_ci_discharge, aes(x = time, y = est,
                                                    color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_disch_lower, ymax = ci_disch_upper),
              alpha = 0.2, linewidth = 0) +
  geom_step(linewidth = 1) +
  labs(
    title = "PSM Cumulative Incidence: Discharge",
    x = "Time from CRRT Initiation (Days)",
    y = "Cumulative Incidence (Discharge)",
    color = "CRRT Dose", fill = "CRRT Dose"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
plot_cif_discharge

ggsave(file.path(output_dir, paste0(SITE_NAME, "_cif_discharge.png")),
       plot_cif_discharge, width = 6, height = 4, dpi = 300)
cat("Saved CIF for discharge as PNG\n")

# ============================================================= #
# ---- 3. SUPER LEARNER PROPENSITY WEIGHTING (IPTW BRANCH) ----
# ============================================================= #

## ---- A. Making IPTW dataset ----

### Super Learner Propensity Weighting ###
set.seed(42)

# Make SL copy of df_tte_bin to avoid overwriting the PSM path
df_tte_sl <- df_tte_bin

# IPTW using Super Learner for propensity score estimation
# Falls back to GLM if SuperLearner fails (e.g., variable length errors at some sites)
w_sl <- tryCatch(
  weightit(
    psm_formula,
    data = df_tte_sl,
    method = "super",
    SL.library = c("SL.glm", "SL.gam", "SL.randomForest"),
    estimand = "ATE",
    stabilize = TRUE
  ),
  error = function(e) {
    cat("\n*** WARNING: SuperLearner failed:", conditionMessage(e), "\n")
    cat("*** Falling back to GLM propensity scores\n\n")
    weightit(
      psm_formula,
      data = df_tte_sl,
      method = "glm",
      estimand = "ATE",
      stabilize = TRUE
    )
  }
)
cat("Propensity score method used:", ifelse(inherits(w_sl, "weightit"),
    w_sl$method, "unknown"), "\n")

# attach SL weights + PS
df_tte_sl$w  <- w_sl$weights
df_tte_sl$ps <- w_sl$ps

## ---- B. Overlap Diagnostics ----
plot_sl_overlap <- ggplot(df_tte_sl, aes(x = ps, fill = crrt_high)) +
  geom_histogram(alpha=0.6, position="identity", bins=30) +
  theme_bw(base_size=10) +
  labs(title="Propensity Score Distribution by CRRT Dose Arm (SuperLearner)",
       x="Propensity Score", y="Count", fill="CRRT Dose")
plot_sl_overlap
ggsave(file.path(output_dir, paste0(SITE_NAME, "_SL_PS_overlap.png")),
       plot_sl_overlap, width=6, height=4)

cat("Saved SL PS overlap plot as PNG\n")

# Positivity diagnostics: extreme propensity score proportions
ps_thresholds <- c(0.05, 0.10, 0.90, 0.95)
ps_extremes <- data.frame(
  threshold = ps_thresholds,
  direction = c("< 0.05", "< 0.10", "> 0.90", "> 0.95"),
  n = sapply(ps_thresholds, function(t) {
    if (t <= 0.5) sum(df_tte_sl$ps < t) else sum(df_tte_sl$ps > t)
  }),
  pct = sapply(ps_thresholds, function(t) {
    if (t <= 0.5) mean(df_tte_sl$ps < t) * 100 else mean(df_tte_sl$ps > t) * 100
  })
)

cat("\nPositivity diagnostics (extreme PS proportions):\n")
cat(sprintf("  %-10s  %5s  %6s\n", "Threshold", "N", "Pct"))
for (i in seq_len(nrow(ps_extremes))) {
  cat(sprintf("  %-10s  %5d  %5.1f%%\n",
              ps_extremes$direction[i], ps_extremes$n[i], ps_extremes$pct[i]))
}
cat("\n")

write.csv(ps_extremes,
          file.path(output_dir, paste0(SITE_NAME, "_SL_PS_extremes.csv")),
          row.names = FALSE)
cat("Saved PS extremes CSV\n\n")

# Weighted balance check BEFORE trimming (Look for Diff.Adj being all <0.1)
# Effective Sampling (aim for close to 1.0)
bal_sl <- cobalt::bal.tab(w_sl, un = TRUE)
bal_df <- as.data.frame(bal_sl$Balance)
bal_sl
write.csv(bal_df, file.path(output_dir,
                            paste0(SITE_NAME, "_SL_IPTW_balance.csv")),
          row.names = TRUE)

ESS <- (sum(df_tte_sl$w)^2) / sum(df_tte_sl$w^2)
ESS_prop <- ESS / nrow(df_tte_sl)
ESS_prop
write_csv(
  tibble(N = nrow(df_tte_sl), ESS = ESS, ESS_prop = ESS_prop),
  file.path(output_dir, paste0(SITE_NAME, "_SL_IPTW_ESS.csv"))
)

### ---- I. Love Plot (SMDs of covariates between groups) ----
# Demonstrates balance of chosen covariates

plot_loveplot_sl <- love.plot(
  bal_sl,
  stats = "mean.diffs",
  abs = TRUE,
  var.names = pretty_names_loveplot,
  thresholds = c(m = .1),
  var.order = "unadjusted",
  line.size = 0.8,
  point.size = 3,
  sample.names = c("Unweighted", "Weighted"),
  title = "Covariate Balance: IPTW via SuperLearner",
  subtitle = "Standardized Mean Differences Before and After Weighting",
  grid = TRUE
)

print(plot_loveplot_sl)

# Save Love Plot
png(file.path(output_dir,
              paste0(SITE_NAME, "_SL_LovePlot.png")),
    width = 8, height = 7, units = "in", res = 300)
print(plot_loveplot_sl)
dev.off()

cat("Saved Love Plot PNG\n")

### ---- II. TABLE S2 (IPTW) ----
# Like Table 1 but median/IQR/% are weighted for the pseudopopulation

# Prepare IPTW df (mirror Table 1 preprocessing)
df_tte_tableS2 <- df_tte_sl %>%
  mutate(
    #### ---- Rename CRRT groups using dose_cutoff
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
    crrt_mode_category = forcats::fct_recode(
      crrt_mode_category,
      "CVVH"   = "cvvh",
      "CVVHD"  = "cvvhd",
      "CVVHDF" = "cvvhdf",
      "SCUF"   = "scuf",
      "AVVH"   = "avvh"
    ),

    #### ---- Three-category outcome
    outcome_3cat = factor(
      case_when(
        outcome == 1 ~ "Discharged",
        outcome == 2 ~ "Died",
        TRUE         ~ "Censored"
      ),
      levels = c("Discharged", "Died", "Censored")
    )
  ) %>%
  # IMPORTANT: Keep only the Table 1 variables + weight
  select(crrt_group, all_of(vars_table1), w)

### ---- Survey design using IPTW weights
design_iptw <- survey::svydesign(
  ids = ~1,
  weights = ~w,
  data = df_tte_tableS2
)


### ---- Build WEIGHTED Table S2
tableS2 <- tbl_svysummary(
  design_iptw,
  by = crrt_group,
  type = table1_type,
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)"
  ),
  label = table1_label
) %>%
  add_overall() %>%
  bold_labels() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(
    c("stat_1","stat_2") ~ "**CRRT dose group**"
  ) %>%
  modify_caption(
    "**Table S2.
    Weighted Baseline Characteristics After IPTW (Super Learner PS)**")

  # Collapse Sex row
  tableS2 <- tableS2 %>%
    collapse_binary_no_p(
      variable   = "sex_category",
      keep_level = "female",
      label      = "Female (%)"
    )

### ---- Save Table S2
gt::gtsave(
  as_gt(tableS2),
  filename = file.path(
    output_dir,
    paste0(SITE_NAME, "_TableS2_IPTW.html")
  )
)
write.csv(as_tibble(tableS2),
          file.path(output_dir, paste0(SITE_NAME, "_TableS2_IPTW.csv")),
          row.names = FALSE)

# ============================================================= #
### ---- OPTIONAL TRIMMING ----
## ---- trimming for non-overlap ---- ONLY DO IF WEIGHTED BAL Diff.Adj > 0.1
#keep <- df_tte_sl$ps >= .05 & df_tte_sl$ps <= .95
#df_tte_sl_trim <- df_tte_sl[keep,]

## winsorize weights (1st and 99th percentiles)
#q <- quantile(df_tte_sl_trim$w, c(.01,.99))
#df_tte_sl_trim$w_trim <- pmin(pmax(df_tte_sl_trim$w, q[1]), q[2])

#cat("SL IPTW: original n =", nrow(df_tte_sl), " trimmed n =",
#nrow(df_tte_sl_trim), "\n")

# Save trimmed SL dataset
#write.csv(df_tte_sl_trim,
#          file.path(output_dir, paste0(SITE_NAME, "_SL_IPTW_trimmed_dataset.csv")),
#          row.names = FALSE)
#cat("Saved SL trimmed dataset\n\n")
# ============================================================= #

## ---- C. Analysis of IPTW ----

### ---- I. Marginal Structural Cause-Specific Cox Models (with IPTW) ----

# References and RowIDs for sandwich SEs
# Dynamically set the reference to the lower-dose group based on cutoff
ref_label <- paste0("<", dose_cutoff)
df_tte_sl$crrt_high <- relevel(df_tte_sl$crrt_high, ref = ref_label)

# Robust sandwich SEs (cluster on row id to be safe)
df_tte_sl$id <- seq_len(nrow(df_tte_sl))

# Build Cox formula dynamically (doubly robust: treatment + covariates)
cox_rhs <- paste(c("crrt_high", model_covariates), collapse = " + ")

# Death cause-specific hazard (event code 2)
fit_cs_death <- coxph(
  as.formula(paste("Surv(time_to_event_30d, outcome == 2) ~", cox_rhs)),
  data = df_tte_sl,
  weights = w,          # IPTW (ATE) from SuperLearner
  robust = TRUE,
  cluster = id
)
summary(fit_cs_death)

# Discharge cause-specific hazard (event code 1)
fit_cs_disch <- coxph(
  as.formula(paste("Surv(time_to_event_30d, outcome == 1) ~", cox_rhs)),
  data = df_tte_sl,
  weights = w,
  robust = TRUE,
  cluster = id
)
summary(fit_cs_disch)

#### ---- II. Extract IPTW Cox cause-specific HR results (full covariates) ----

extract_iptw_cox <- function(fit, label){
  s <- summary(fit)
  co <- s$coefficients
  ci <- confint(fit)    # CI already on coef scale

  # detect correct column names
  z_col <- intersect(colnames(co), c("robust z","z"))[1]
  p_col <- intersect(colnames(co), c("Pr(>|z|)","Robust Pr(>|z|)","p"))[1]
  se_col <- intersect(colnames(co), c("robust se","se(coef)"))[1]

  data.frame(
    outcome = label,
    variable = rownames(co),
    HR = exp(co[, "coef"]),
    HR_lower = exp(ci[,1]),
    HR_upper = exp(ci[,2]),
    se_log_hr = co[, se_col],
    z = co[, z_col],
    p_value = co[, p_col],
    row.names = NULL,
    check.names = FALSE
  )
}
# Bind and write csv
iptw_results <- bind_rows(
  extract_iptw_cox(fit_cs_death,  "Death (IPTW Cause-Specific Cox)"),
  extract_iptw_cox(fit_cs_disch,  "Discharge (IPTW Cause-Specific Cox)")
)
write.csv(
  iptw_results,
  file.path(output_dir, paste0(SITE_NAME, "_IPTW_CauseSpecificCox_FULLresults.csv")),
  row.names = FALSE
)
cat("Saved full IPTW Cause-Specific Cox results.\n")

### ---- II-b. Rubin's Rules Pooled IPTW Cox Models (across 5 imputations) ----

cat("\nPooling IPTW Cox models across", N_IMP, "MICE imputations...\n")

# Build imputed SL datasets: same data prep as df_tte_sl but on each imputation
# Reuse same SL weights (fitted on imputation 1) — standard pragmatic approach
imp_sl_list <- lapply(seq_len(N_IMP), function(m) {
  d <- imp_list[[m]]
  d$crrt_high <- ifelse(d$crrt_dose_ml_kg_hr_0 >= dose_cutoff, 1L, 0L)
  d$crrt_high <- factor(d$crrt_high, levels = c(0, 1),
                        labels = c(paste0("<", dose_cutoff),
                                   paste0(">=", dose_cutoff)))
  d <- tidyr::drop_na(d, dplyr::all_of(complete_case_vars))
  d <- droplevels(d)
  d$crrt_high <- relevel(d$crrt_high, ref = ref_label)
  d$w  <- w_sl$weights
  d$ps <- w_sl$ps
  d$id <- seq_len(nrow(d))
  d
})

# Fit Cox death + discharge on each imputed dataset
cox_formula_death <- as.formula(
  paste("Surv(time_to_event_30d, outcome == 2) ~", cox_rhs))
cox_formula_disch <- as.formula(
  paste("Surv(time_to_event_30d, outcome == 1) ~", cox_rhs))

imp_cox_fits <- lapply(seq_len(N_IMP), function(m) {
  d <- imp_sl_list[[m]]
  fit_d <- coxph(cox_formula_death, data = d, weights = w,
                 robust = TRUE, cluster = id)
  fit_c <- coxph(cox_formula_disch, data = d, weights = w,
                 robust = TRUE, cluster = id)
  list(death = fit_d, disch = fit_c)
})

# Extract treatment coefficient from each fit and pool
p_col_detect <- function(co) {
  intersect(colnames(co), c("Pr(>|z|)", "Robust Pr(>|z|)"))[1]
}

pool_treatment_hr <- function(fits, outcome_label) {
  trt_name <- grep("crrt_high", names(coef(fits[[1]])), value = TRUE)[1]
  betas <- sapply(fits, function(f) coef(f)[trt_name])
  ses   <- sapply(fits, function(f) {
    co <- summary(f)$coefficients
    co[trt_name, "robust se"]
  })
  pooled <- rubins_pool(betas, ses)
  data.frame(
    model     = paste0("IPTW Cox (pooled) - ", outcome_label),
    HR_type   = "CauseSpecific HR",
    HR        = pooled$hr,
    HR_lower  = pooled$hr_lower,
    HR_upper  = pooled$hr_upper,
    se_log_hr = pooled$se,
    p_value   = pooled$p_value,
    fmi       = pooled$fmi,
    stringsAsFactors = FALSE
  )
}

pooled_iptw_death <- pool_treatment_hr(
  lapply(imp_cox_fits, `[[`, "death"), "Death")
pooled_iptw_disch <- pool_treatment_hr(
  lapply(imp_cox_fits, `[[`, "disch"), "Discharge")

pooled_iptw_results <- bind_rows(pooled_iptw_death, pooled_iptw_disch)

cat("Pooled IPTW Cox results:\n")
for (i in seq_len(nrow(pooled_iptw_results))) {
  r <- pooled_iptw_results[i, ]
  cat(sprintf("  %-35s HR=%.2f (%.2f-%.2f) p=%.4f FMI=%.3f\n",
              r$model, r$HR, r$HR_lower, r$HR_upper, r$p_value, r$fmi))
}
cat("\n")

write.csv(pooled_iptw_results,
          file.path(output_dir, paste0(SITE_NAME, "_IPTW_pooled_results.csv")),
          row.names = FALSE)
cat("Saved pooled IPTW results CSV\n")

### ---- III. Standardized CIF Curves (competing risks) ----

# Helper: carry-forward cumulative baseline hazard to a common grid
cumhaz_at_grid <- function(basehaz_df, grid_time) {
  idx <- findInterval(grid_time, basehaz_df$time)
  out <- numeric(length(grid_time))
  out[idx > 0] <- basehaz_df$hazard[idx[idx > 0]]
  out
}

# Population-averaged standardized CIFs via g-computation
build_standardized_cifs <- function(fit_death, fit_disch, trt_var, trt_levels, newdata) {
  bh_death <- basehaz(fit_death, centered = FALSE)
  bh_disch <- basehaz(fit_disch, centered = FALSE)
  tgrid <- sort(unique(c(bh_death$time, bh_disch$time)))
  tgrid <- tgrid[tgrid >= 0]
  H0_death <- cumhaz_at_grid(bh_death, tgrid)
  H0_disch <- cumhaz_at_grid(bh_disch, tgrid)
  dH0_death <- c(H0_death[1], diff(H0_death))
  dH0_disch <- c(H0_disch[1], diff(H0_disch))
  nt <- length(tgrid)
  
  out <- lapply(trt_levels, function(lv) {
    cf_data <- newdata
    cf_data[[trt_var]] <- factor(lv, levels = trt_levels)
    lp_death <- predict(fit_death, newdata = cf_data, type = "lp")
    lp_disch <- predict(fit_disch, newdata = cf_data, type = "lp")
    mult_death <- exp(lp_death)
    mult_disch <- exp(lp_disch)
    n <- length(mult_death)
    S_prev <- rep(1, n)
    cif_death_accum <- rep(0, n)
    cif_disch_accum <- rep(0, n)
    avg_cif_death <- numeric(nt)
    avg_cif_disch <- numeric(nt)
    avg_surv <- numeric(nt)
    for (j in seq_len(nt)) {
      dH_death_j <- dH0_death[j] * mult_death
      dH_disch_j <- dH0_disch[j] * mult_disch
      cif_death_accum <- cif_death_accum + S_prev * dH_death_j
      cif_disch_accum <- cif_disch_accum + S_prev * dH_disch_j
      S_prev <- S_prev * exp(-(dH_death_j + dH_disch_j))
      avg_cif_death[j] <- mean(cif_death_accum)
      avg_cif_disch[j] <- mean(cif_disch_accum)
      avg_surv[j] <- mean(S_prev)
    }
    tibble::tibble(time = tgrid, cif_death = avg_cif_death,
                   cif_disch = avg_cif_disch, surv = avg_surv,
                   trt_group = lv)
  }) %>%
    dplyr::bind_rows() %>%
    mutate(trt_group = factor(trt_group, levels = trt_levels))
  
  out
}

# Generate point-estimate CIFs
  trt_levels <- levels(df_tte_sl$crrt_high)

cif_df <- build_standardized_cifs(
  fit_death = fit_cs_death,
  fit_disch = fit_cs_disch,
  trt_var = "crrt_high",
  trt_levels = trt_levels,
  newdata = df_tte_sl
)

# Bootstrap CIs (200 replicates)
N_BOOT <- 500
set.seed(42)
cat("Bootstrapping IPTW CIF confidence intervals (", N_BOOT, "replicates) …\n")

tgrid_common <- sort(unique(cif_df$time))
boot_death <- boot_disch <- list()
for (lv in trt_levels) {
  boot_death[[lv]] <- matrix(NA, nrow = N_BOOT, ncol = length(tgrid_common))
  boot_disch[[lv]] <- matrix(NA, nrow = N_BOOT, ncol = length(tgrid_common))
}

unique_ids <- seq_len(nrow(df_tte_sl))  # row-level resampling (no cluster ID)

for (b in seq_len(N_BOOT)) {
  boot_idx <- sample(unique_ids, length(unique_ids), replace = TRUE)
  boot_df <- df_tte_sl[boot_idx, ]
  boot_df$id <- seq_len(nrow(boot_df))
  
  fit_b <- tryCatch({
    cox_rhs_boot <- paste(c("crrt_high", model_covariates), collapse = " + ")
    fit_d <- coxph(
      as.formula(paste("Surv(time_to_event_30d, outcome == 2) ~", cox_rhs_boot)),
      data = boot_df, weights = w, robust = FALSE
    )
    fit_c <- coxph(
      as.formula(paste("Surv(time_to_event_30d, outcome == 1) ~", cox_rhs_boot)),
      data = boot_df, weights = w, robust = FALSE
    )
    list(death = fit_d, disch = fit_c)
  }, error = function(e) NULL)
  
  if (is.null(fit_b)) next
  
  boot_cif <- tryCatch(
    build_standardized_cifs(fit_b$death, fit_b$disch, "crrt_high",
                            trt_levels, newdata = boot_df),
    error = function(e) NULL
  )
  if (is.null(boot_cif)) next
  
  for (lv in trt_levels) {
    lv_df <- boot_cif %>% filter(trt_group == lv)
    if (nrow(lv_df) == 0) next
    boot_death[[lv]][b, ] <- approx(lv_df$time, lv_df$cif_death,
                                    xout = tgrid_common, rule = 2)$y
    boot_disch[[lv]][b, ] <- approx(lv_df$time, lv_df$cif_disch,
                                    xout = tgrid_common, rule = 2)$y
  }
  if (b %% 100 == 0) cat("  ", b, "/", N_BOOT, "\n")
}

# Compute 95% CIs
ci_list <- lapply(trt_levels, function(lv) {
  tibble(time = tgrid_common, trt_group = lv,
         cif_death_lower = apply(boot_death[[lv]], 2, quantile, 0.025, na.rm = TRUE),
         cif_death_upper = apply(boot_death[[lv]], 2, quantile, 0.975, na.rm = TRUE),
         cif_disch_lower = apply(boot_disch[[lv]], 2, quantile, 0.025, na.rm = TRUE),
         cif_disch_upper = apply(boot_disch[[lv]], 2, quantile, 0.975, na.rm = TRUE))
})
ci_df <- bind_rows(ci_list) %>%
  mutate(trt_group = factor(trt_group, levels = trt_levels))

cif_df <- cif_df %>% left_join(ci_df, by = c("time", "trt_group"))
cat("Bootstrap complete.\n")

# Plot CIF Death
plot_cif_death <- ggplot(cif_df, aes(x = time, y = cif_death,
                                     color = trt_group, fill = trt_group)) +
  geom_ribbon(aes(ymin = cif_death_lower, ymax = cif_death_upper),
              alpha = 0.2, linewidth = 0) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "IPTW Standardized CIF: Death",
       x = "Time from CRRT Initiation (Days)",
       y = "Cumulative Incidence (Death)",
       color = "CRRT Dose", fill = "CRRT Dose")
plot_cif_death
ggsave(file.path(output_dir, paste0(SITE_NAME, "_IPTW_CIF_Death.png")),
       plot_cif_death, width = 6, height = 4, dpi = 300)

# Plot CIF Discharge
plot_cif_disch <- ggplot(cif_df, aes(x = time, y = cif_disch,
                                     color = trt_group, fill = trt_group)) +
  geom_ribbon(aes(ymin = cif_disch_lower, ymax = cif_disch_upper),
              alpha = 0.2, linewidth = 0) +
  geom_line(linewidth = 1) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "IPTW Standardized CIF: Discharge",
       x = "Time from CRRT Initiation (Days)",
       y = "Cumulative Incidence (Discharge)",
       color = "CRRT Dose", fill = "CRRT Dose")
plot_cif_disch
ggsave(file.path(output_dir, paste0(SITE_NAME, "_IPTW_CIF_Discharge.png")),
       plot_cif_disch, width = 6, height = 4, dpi = 300)

# Export IPTW CIF data as CSV
write.csv(cif_df,
          file.path(output_dir, paste0(SITE_NAME, "_IPTW_CIF_data.csv")),
          row.names = FALSE)
cat("IPTW standardized CIF plots with CIs saved as PNGs + CSV.\n")

# ============================================================= #
# ----  4. MODEL COMPARISON - PSM FG vs IPTW Cox ----
# ============================================================= #

# Extract FG SHR treatment rows
extract_fg_trt <- function(fg, label){
  s <- summary(fg)$coef
  tr <- s[grep("crrt_high", rownames(s)), , drop=FALSE]
  data.frame(
    model = label,
    HR_type = "SHR",
    HR = tr[,"exp(coef)"],
    HR_lower = tr[,"exp(coef)"]/exp(1.96*tr[,"se(coef)"]),
    HR_upper = tr[,"exp(coef)"]*exp(1.96*tr[,"se(coef)"]),
    se_log_hr = tr[,"se(coef)"],
    p_value = tr[,"p-value"],
    row.names = NULL
  )
}

fg_death_row  <- extract_fg_trt(fg_death_psm_dr,  "PSM FG - Death")
fg_disch_row  <- extract_fg_trt(fg_disch_psm_dr,  "PSM FG - Discharge")

# Extract IPTW Cox treatment rows
extract_iptw_trt <- function(fit, label){
  s <- summary(fit)
  co <- s$coefficients
  ci <- confint(fit)
  idx <- grep("crrt_high", rownames(co))

  # Detect column names dynamically
  p_col <- intersect(colnames(co), c("Pr(>|z|)", "Robust Pr(>|z|)"))[1]
  se_col <- intersect(colnames(co), c("robust se", "se(coef)"))[1]

  data.frame(
    model = label,
    HR_type = "CauseSpecific HR",
    HR = exp(co[idx,"coef"]),
    HR_lower = exp(ci[idx,1]),
    HR_upper = exp(ci[idx,2]),
    se_log_hr = co[idx, se_col],
    p_value = co[idx, p_col],
    row.names=NULL
  )
}

iptw_death_row <- extract_iptw_trt(fit_cs_death, "IPTW Cox - Death")
iptw_disch_row <- extract_iptw_trt(fit_cs_disch, "IPTW Cox - Discharge")

comparison_table <- bind_rows(
  fg_death_row,
  fg_disch_row,
  iptw_death_row,
  iptw_disch_row,
  pooled_fg_results %>% select(model, HR_type, HR, HR_lower, HR_upper, se_log_hr, p_value),
  pooled_iptw_results %>% select(model, HR_type, HR, HR_lower, HR_upper, se_log_hr, p_value)
)

comparison_table

write.csv(comparison_table,
          file.path(output_dir, paste0(SITE_NAME, "_ModelComparison_PSMvsIPTW.csv")),
          row.names = FALSE)

cat("Saved model comparison CSV\n")

# ============================================================= #
# ---- 5. MODEL COMPARISON COMPLETE ----
# ============================================================= #
cat("PSM and IPTW analyses complete.\n\n")

# ============================================================= #
# ---- 5b. E-VALUE SENSITIVITY ANALYSIS ----
# ============================================================= #

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Section 5b: E-Value Sensitivity Analysis for Unmeasured Confounding\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# E-values quantify the minimum strength of association that an unmeasured
# confounder would need to have with both treatment AND outcome to fully
# explain away the observed effect. Higher E-values = more robust to
# unmeasured confounding.

# Extract treatment-row HRs from comparison_table (already computed)
# comparison_table has: model, HR_type, HR, HR_lower, HR_upper, p_value

compute_evalue_row <- function(model_label, hr, hr_lower, hr_upper) {
  # evalues.HR expects HR >= 1 for the point estimate direction;
  # it handles HR < 1 internally by computing 1/HR
  ev <- evalues.HR(est = hr, lo = hr_lower, hi = hr_upper, rare = FALSE)

  # ev is a matrix; row 1 = point estimate, row 2 = CI closest to null
  data.frame(
    model       = model_label,
    HR          = hr,
    HR_lower    = hr_lower,
    HR_upper    = hr_upper,
    evalue_point = ev["E-values", "point"],
    evalue_ci    = ev["E-values", "lower"],
    stringsAsFactors = FALSE
  )
}

# Compute E-values for all models in the comparison table
evalue_results <- bind_rows(lapply(seq_len(nrow(comparison_table)), function(i) {
  compute_evalue_row(comparison_table$model[i],
                     comparison_table$HR[i],
                     comparison_table$HR_lower[i],
                     comparison_table$HR_upper[i])
}))

cat("E-value results:\n")
cat(sprintf("  %-25s  HR=%5.2f (%4.2f-%4.2f)  E-value: point=%.2f, CI=%.2f\n",
            evalue_results$model,
            evalue_results$HR,
            evalue_results$HR_lower,
            evalue_results$HR_upper,
            evalue_results$evalue_point,
            evalue_results$evalue_ci))
cat("\n")

write.csv(evalue_results,
          file.path(output_dir, paste0(SITE_NAME, "_evalue_sensitivity.csv")),
          row.names = FALSE)
cat("Saved E-value sensitivity CSV\n\n")

# ============================================================= #
# ---- 6. HYPOTHESIS-GENERATING SUBGROUP ANALYSIS ----
# ============================================================= #

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Section 6: Subgroup Analysis (Hypothesis-Generating)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

## ---- A. Define subgroup variables with clinical cutpoints ----

# Compute CCI total (sum of 17 binary components)
df_tte_sl$cci_total <- rowSums(df_tte_sl[, cci_vars], na.rm = TRUE)

# Named list of subgroup definitions
# Each entry: var_expr (evaluated on df), label, reference level, factor levels
subgroup_defs <- list(
  age_group = list(
    var_expr = quote(ifelse(age_at_admission >= 65, ">=65", "<65")),
    label = "Age (years)",
    levels = c("<65", ">=65")
  ),
  sofa_group = list(
    var_expr = quote(ifelse(sofa_total_0 >= 12, ">=12", "<12")),
    label = "SOFA Score",
    levels = c("<12", ">=12")
  ),
  lactate_group = list(
    var_expr = quote(ifelse(lactate_0 >= 4, ">=4", "<4")),
    label = "Lactate (mmol/L)",
    levels = c("<4", ">=4")
  ),
  sex_group = list(
    var_expr = quote(ifelse(sex_category == "female", "Female", "Male")),
    label = "Sex",
    levels = c("Male", "Female")
  ),
  imv_group = list(
    var_expr = quote(ifelse(imv_status_0 == 1, "On IMV", "No IMV")),
    label = "Mechanical Ventilation",
    levels = c("No IMV", "On IMV")
  ),
  vasopressor_group = list(
    var_expr = quote(ifelse(norepinephrine_equivalent_0 > 0.1, "High", "Low")),
    label = "Vasopressor (NEE >0.1)",
    levels = c("Low", "High")
  ),
  crrt_mode_group = list(
    var_expr = quote(ifelse(crrt_mode_category == "cvvhdf", "CVVHDF", "Other")),
    label = "CRRT Modality",
    levels = c("Other", "CVVHDF")
  ),
  cci_group = list(
    var_expr = quote(ifelse(cci_total >= 4, ">=4", "<4")),
    label = "CCI Total",
    levels = c("<4", ">=4")
  ),
  weight_group = list(
    var_expr = quote(ifelse(weight_kg >= 80, ">=80 kg", "<80 kg")),
    label = "Body Weight",
    levels = c("<80 kg", ">=80 kg")
  )
)

# Create subgroup columns on the IPTW analysis dataset
cat("Subgroup distributions:\n")
for (sg_name in names(subgroup_defs)) {
  sg <- subgroup_defs[[sg_name]]
  df_tte_sl[[sg_name]] <- factor(
    eval(sg$var_expr, envir = df_tte_sl),
    levels = sg$levels
  )
  tbl <- table(df_tte_sl[[sg_name]])
  cat("  ", sg$label, ": ",
      paste(paste0(names(tbl), "=", tbl), collapse = ", "), "\n")
}
cat("\n")

## ---- B. Fit interaction models (pooled across imputations) ----

cat("Fitting interaction models for death outcome (pooled across",
    N_IMP, "imputations)...\n")
subgroup_results <- list()

# Create subgroup columns on each imputed dataset
for (m in seq_len(N_IMP)) {
  d <- imp_sl_list[[m]]
  d$cci_total <- rowSums(d[, cci_vars], na.rm = TRUE)
  for (sg_name in names(subgroup_defs)) {
    sg <- subgroup_defs[[sg_name]]
    d[[sg_name]] <- factor(
      eval(sg$var_expr, envir = d),
      levels = sg$levels
    )
  }
  imp_sl_list[[m]] <- d
}

for (sg_name in names(subgroup_defs)) {
  sg <- subgroup_defs[[sg_name]]

  # Skip if either subgroup level has < 20 death events (check on imp 1)
  event_counts <- tapply(df_tte_sl$outcome == 2, df_tte_sl[[sg_name]], sum)
  if (any(event_counts < 20, na.rm = TRUE)) {
    cat("  Skipping ", sg$label,
        " -- insufficient events (", paste(event_counts, collapse = "/"), ")\n")
    next
  }

  # Formula: Surv(...) ~ crrt_high * subgroup_var + covariates
  int_formula <- as.formula(paste(
    "Surv(time_to_event_30d, outcome == 2) ~",
    "crrt_high *", sg_name, "+",
    paste(model_covariates, collapse = " + ")
  ))

  # Fit on each imputed dataset and collect coefficients
  betas_ref <- numeric(N_IMP)     # Treatment effect at reference level
  ses_ref   <- numeric(N_IMP)
  betas_alt <- numeric(N_IMP)     # Treatment effect at non-reference level
  ses_alt   <- numeric(N_IMP)
  betas_int <- numeric(N_IMP)     # Interaction coefficient
  ses_int   <- numeric(N_IMP)
  fit_ok    <- logical(N_IMP)

  for (m in seq_len(N_IMP)) {
    fit_int <- tryCatch(
      coxph(int_formula,
            data = imp_sl_list[[m]],
            weights = w,
            robust = TRUE,
            cluster = id),
      error = function(e) NULL
    )
    if (is.null(fit_int)) next
    fit_ok[m] <- TRUE

    beta <- coef(fit_int)
    V <- vcov(fit_int)

    trt_idx <- grep("^crrt_high", names(beta))[1]
    int_rows <- grep(
      paste0("crrt_high.*:", sg_name, "|", sg_name, ":.*crrt_high"),
      names(beta)
    )

    # Reference level: treatment main effect
    betas_ref[m] <- beta[trt_idx]
    ses_ref[m]   <- sqrt(V[trt_idx, trt_idx])

    # Non-reference level: treatment + interaction
    if (length(int_rows) > 0) {
      int_idx <- int_rows[1]
      betas_int[m] <- beta[int_idx]
      ses_int[m]   <- sqrt(V[int_idx, int_idx])

      betas_alt[m] <- beta[trt_idx] + beta[int_idx]
      ses_alt[m]   <- sqrt(V[trt_idx, trt_idx] + V[int_idx, int_idx] +
                             2 * V[trt_idx, int_idx])
    }
  }

  if (sum(fit_ok) < 3) {
    cat("  Skipping ", sg$label, " -- too few successful fits\n")
    next
  }

  # Pool using Rubin's rules
  pooled_ref <- rubins_pool(betas_ref[fit_ok], ses_ref[fit_ok])
  pooled_alt <- rubins_pool(betas_alt[fit_ok], ses_alt[fit_ok])
  pooled_int <- rubins_pool(betas_int[fit_ok], ses_int[fit_ok])

  # Count n and death events per subgroup level (from imputation 1)
  n_ref <- sum(df_tte_sl[[sg_name]] == sg$levels[1], na.rm = TRUE)
  n_alt <- sum(df_tte_sl[[sg_name]] == sg$levels[2], na.rm = TRUE)
  ev_ref <- sum(df_tte_sl[[sg_name]] == sg$levels[1] &
                  df_tte_sl$outcome == 2, na.rm = TRUE)
  ev_alt <- sum(df_tte_sl[[sg_name]] == sg$levels[2] &
                  df_tte_sl$outcome == 2, na.rm = TRUE)

  subgroup_results[[sg_name]] <- tibble(
    subgroup    = sg$label,
    level       = c(sg$levels[1], sg$levels[2]),
    n           = c(n_ref, n_alt),
    events      = c(ev_ref, ev_alt),
    HR          = c(pooled_ref$hr, pooled_alt$hr),
    HR_lower    = c(pooled_ref$hr_lower, pooled_alt$hr_lower),
    HR_upper    = c(pooled_ref$hr_upper, pooled_alt$hr_upper),
    p_interaction = pooled_int$p_value
  )

  cat("  ", sg$label, ": p-interaction =",
      sprintf("%.4f", pooled_int$p_value), "\n")
}
cat("\n")

## ---- C. Assemble results and compute FDR q-values ----

sg_results_df <- bind_rows(subgroup_results)

# Add overall HR from the pooled IPTW model for context
overall_row <- tibble(
  subgroup      = "Overall",
  level         = "All patients",
  n             = nrow(df_tte_sl),
  events        = sum(df_tte_sl$outcome == 2),
  HR            = pooled_iptw_death$HR,
  HR_lower      = pooled_iptw_death$HR_lower,
  HR_upper      = pooled_iptw_death$HR_upper,
  p_interaction = NA
)

sg_results_df <- bind_rows(overall_row, sg_results_df)

# Compute BH-adjusted FDR q-values on the unique interaction p-values
unique_p <- sg_results_df %>%
  filter(!is.na(p_interaction)) %>%
  distinct(subgroup, .keep_all = TRUE) %>%
  mutate(q_interaction = p.adjust(p_interaction, method = "BH"))

sg_results_df <- sg_results_df %>%
  left_join(unique_p %>% select(subgroup, q_interaction), by = "subgroup")

# Format for display
sg_results_df <- sg_results_df %>%
  mutate(
    HR_label = sprintf("%.2f (%.2f-%.2f)", HR, HR_lower, HR_upper),
    p_int_label = ifelse(is.na(p_interaction), "",
                         sprintf("%.3f", p_interaction)),
    q_int_label = ifelse(is.na(q_interaction), "",
                         sprintf("%.3f", q_interaction))
  )

# Save CSV
write.csv(sg_results_df,
          file.path(output_dir, paste0(SITE_NAME, "_subgroup_analysis_results.csv")),
          row.names = FALSE)
cat("Saved subgroup analysis results CSV\n")

# Print summary
cat("\nSubgroup Analysis Summary (Death Outcome):\n")
cat(sprintf("%-25s %-12s %5s %6s   %-22s  %s\n",
            "Subgroup", "Level", "N", "Events", "HR (95% CI)", "p-int"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in seq_len(nrow(sg_results_df))) {
  r <- sg_results_df[i, ]
  cat(sprintf("%-25s %-12s %5d %6d   %-22s  %s\n",
              r$subgroup, r$level, r$n, r$events,
              r$HR_label, r$p_int_label))
}
cat("\n")

## ---- D. Forest plot of subgroup-specific HRs ----

# Filter out subgroup levels with NA HR (e.g., empty subgroups like CVVHDF=0)
sg_plot_data <- sg_results_df %>%
  filter(subgroup != "Overall") %>%
  filter(!is.na(HR))

# Also drop subgroups that lost a level after NA filtering
sg_with_both <- sg_plot_data %>%
  count(subgroup) %>%
  filter(n == 2) %>%
  pull(subgroup)
sg_plot_data <- sg_plot_data %>% filter(subgroup %in% sg_with_both)

# Build plot data with row ordering:
#   - Overall at the top
#   - Then each subgroup header followed by its two levels
plot_rows <- list()
row_counter <- 0

# Overall row
row_counter <- row_counter + 1
plot_rows[[row_counter]] <- tibble(
  row_label = "Overall", is_header = FALSE, is_overall = TRUE,
  HR = overall_row$HR, HR_lower = overall_row$HR_lower,
  HR_upper = overall_row$HR_upper,
  p_interaction = NA, row_order = row_counter
)

# Spacer after overall
row_counter <- row_counter + 1
plot_rows[[row_counter]] <- tibble(
  row_label = "", is_header = FALSE, is_overall = FALSE,
  HR = NA, HR_lower = NA, HR_upper = NA,
  p_interaction = NA, row_order = row_counter
)

# Each subgroup: header + two levels (only subgroups with both levels)
for (sg_label in unique(sg_plot_data$subgroup)) {
  sg_data <- sg_plot_data %>% filter(subgroup == sg_label)

  # Subgroup header row
  row_counter <- row_counter + 1
  plot_rows[[row_counter]] <- tibble(
    row_label = sg_label, is_header = TRUE, is_overall = FALSE,
    HR = NA, HR_lower = NA, HR_upper = NA,
    p_interaction = sg_data$p_interaction[1], row_order = row_counter
  )

  # Level rows
  for (j in seq_len(nrow(sg_data))) {
    row_counter <- row_counter + 1
    r <- sg_data[j, ]
    plot_rows[[row_counter]] <- tibble(
      row_label = paste0("    ", r$level),
      is_header = FALSE, is_overall = FALSE,
      HR = r$HR, HR_lower = r$HR_lower, HR_upper = r$HR_upper,
      p_interaction = NA, row_order = row_counter
    )
  }
}

plot_df <- bind_rows(plot_rows)

# Reverse y so Overall is at the top
plot_df$y_pos <- max(plot_df$row_order) - plot_df$row_order + 1

# Determine x-axis range for annotation placement
x_max_ci <- max(plot_df$HR_upper, na.rm = TRUE)
x_annot <- x_max_ci * 1.5

# Data for points (non-header, non-spacer rows with valid HR)
point_df <- plot_df %>%
  filter(!is.na(HR)) %>%
  mutate(pt_color = ifelse(is_overall, "black", "steelblue"))

# Build forest plot
p_forest <- ggplot(plot_df, aes(y = y_pos)) +
  # Reference line at HR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  # Point + CI
  geom_pointrange(
    data = point_df,
    aes(x = HR, xmin = HR_lower, xmax = HR_upper, color = pt_color),
    size = 0.4, fatten = 3, show.legend = FALSE
  ) +
  scale_color_identity() +
  # P-interaction annotations for header rows
  geom_text(
    data = plot_df %>% filter(is_header & !is.na(p_interaction)),
    aes(x = x_annot, y = y_pos,
        label = paste0("p-int = ", sprintf("%.3f", p_interaction))),
    hjust = 0, size = 3, fontface = "italic", color = "grey30"
  ) +
  # Log scale x-axis
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4")
  ) +
  # Y-axis labels
  scale_y_continuous(
    breaks = plot_df$y_pos,
    labels = plot_df$row_label,
    expand = expansion(add = 0.8)
  ) +
  coord_cartesian(xlim = c(
    min(0.2, min(plot_df$HR_lower, na.rm = TRUE) * 0.8),
    x_annot * 2
  )) +
  labs(
    title = "Subgroup Analysis: CRRT Dose and 30-Day Mortality",
    subtitle = paste0("IPTW Cause-Specific HR (High vs Low Dose, cutoff = ",
                      dose_cutoff, " mL/kg/hr)"),
    x = "Hazard Ratio (log scale)",
    y = NULL,
    caption = paste0(
      "Hypothesis-generating analysis. P-interaction values are unadjusted.\n",
      "Subgroups with <20 death events in either level were excluded."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(hjust = 0, size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    plot.caption = element_text(hjust = 0, face = "italic", size = 8),
    plot.margin = ggplot2::margin(10, 40, 10, 10, unit = "pt")
  )

print(p_forest)

# Save PNG and PDF
ggsave(file.path(output_dir, paste0(SITE_NAME, "_subgroup_forest_plot.png")),
       p_forest, width = 10, height = 8, dpi = 300)

ggsave(file.path(output_dir, paste0(SITE_NAME, "_subgroup_forest_plot.pdf")),
       p_forest, width = 10, height = 8)

cat("Saved subgroup forest plot (PNG + PDF)\n")

# ============================================================= #
# ---- 7. FINISH! ----
# ============================================================= #
cat("\nAll analyses complete!\n")
