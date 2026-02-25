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
  imputed_block <- complete(imp, 1)
  df_complete <- df
  df_complete[, mice_vars] <- imputed_block
  cat("MICE imputed", sum(miss_counts[miss_counts > 0]), "total missing values\n")
} else {
  cat("No missing values — skipping MICE\n")
  df_complete <- df
}

cat("Analysis sample:", nrow(df_complete), "of", nrow(df), "\n\n")

# Check for sufficient data
if (nrow(df_complete) < 50) {
  stop("Insufficient cases for modeling (n = ", nrow(df_complete), ")")
}

## ---- A2. Build model covariate set and formulas dynamically ----

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

## ---- B. CRRT Dose Cutoff ----
# =================================== #

# CRRT dose cutoff (mL/kg/hr)
dose_cutoff <- 25

# Generate a label for filenames and titles
dose_label <- paste0("CRRT_", dose_cutoff, "cutoff")

cat("Using CRRT dose cutoff of:", dose_cutoff, "mL/kg/hr\n\n")
# =================================== #

## ---- C. Calculate summary statistics ----
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
          file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                       "_sample_characteristics.csv")),
          row.names = FALSE)

cat("Sample characteristics saved.\n\n")

## ---- D. Histogram CRRT dose ----
crrt_dose_hist <- ggplot(df_complete, aes(x = crrt_dose_ml_kg_hr_0)) +
  geom_histogram(binwidth = 5, color = "black", fill = "skyblue") +
  geom_vline(xintercept = dose_cutoff, linetype = "dashed", linewidth = 1) +
      # Line at cutoff
  labs(
    title = paste0("Distribution of CRRT Dose (Cutoff = ", dose_cutoff, " mL/kg/hr)"),
    x = "CRRT Dose (mL/kg/hr)",
    y = "Count"
  ) +
  theme_bw(base_size = 12)
crrt_dose_hist

# Save histogram
ggsave(file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                    "_hist_crrt_dose.png")),
       crrt_dose_hist, width=6, height=4)
cat("Saved histogram PNG\n")

## ---- E. Designate High vs Low CRRT Dose ----
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
  drop_na(
    crrt_high, age_at_admission, sex_category, sofa_total_0,
    lactate_0, bicarbonate_0, potassium_0,
    oxygenation_index_0, norepinephrine_equivalent_0, imv_status_0,
    time_to_event_90d, outcome
  )

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


## ---- TABLE 1 ----
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
    race_category = forcats::fct_recode(
      race_category,
      "Black/African American" = "black or african american",
      "White" = "white",
      "Asian" = "asian",
      "Native American" = "american indian or alaska native",
      "Pacific Islander" = "native hawaiian or other pacific islander",
      "Unknown" = "unknown",
      "Other" = "other"
    ),

    #### ---- Clean CRRT modality labels
    crrt_mode_category = forcats::fct_recode(
      crrt_mode_category,
      "CVVH" = "cvvh",
      "CVVHD" = "cvvhd",
      "CVVHDF" = "cvvhdf",
      "SCUF" = "scuf",
      "AVVH" = "avvh"
    ),

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
  outcome_3cat                 ~ "90-day Outcome"
)
for (v in cci_vars) {
  table1_label[[length(table1_label) + 1]] <- as.formula(
    paste0(v, ' ~ "', cci_labels[v], '"'))
}

table1 <- df_tte_table1 %>%
  select(crrt_group, all_of(vars_table1)
         ) %>%
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
  add_p(
    test = list(
      race_category ~ "chisq.test.no.correct",
      outcome_3cat  ~ "chisq.test.no.correct"
    )
  ) %>%
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
                       paste0(SITE_NAME, "_",
                              dose_label, "_Table1_unadjusted.html"))
)

# ============================================ #
# ---- 2. PROPENSITY SCORE MATCHING BRANCH ----
# ============================================ #

## ---- A. Making and visualizing PSM dataset ----
### ---- PSM ----
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
write.csv(psm_bal, file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                                "_psm_balance_summary.csv")),
          row.names = FALSE)
write.csv(psm_counts, file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                                "_psm_counts_summary.csv")),
          row.names = FALSE)
cat("Saved PSM summaries\n")

### ---- Save Matched Dataset ----
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

### ---- Love Plot (SMD + variance ratio) ----
plot_loveplot_psm <- love.plot(
  m.out,
  stat = c("m", "v"),
  grid = TRUE,
  var.names = pretty_names_loveplot,
  title = "Covariate Balance: PSM",
  threshold = c(m = .25, v = 1.25)
)

print(plot_loveplot_psm)

png(file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_psm_loveplot.png")),
    width = 8, height = 7, units = "in", res = 300)
print(plot_loveplot_psm)
dev.off()

cat("Saved Love plot for PSM as PNG\n")

### ---- TABLE S1 (PSM) ----

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
    race_category = forcats::fct_recode(
      race_category,
      "Black/African American"  = "black or african american",
      "White"                   = "white",
      "Asian"                   = "asian",
      "Native American"         = "american indian or alaska native",
      "Pacific Islander"        = "native hawaiian or other pacific islander",
      "Unknown"                 = "unknown",
      "Other"                   = "other"
    ),

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
    paste0(SITE_NAME, "_", dose_label, "_TableS1_matched.html")
  )
)

## ---- B. Analysis of PSM ----
### ---- Fine-Gray Analysis on Matched Patients ----

# Matrix for treatment only (reference = low dose)
X_trt <- model.matrix(~ crrt_high, df_match)[, -1, drop = FALSE]

fg_death_psm <- cmprsk::crr(
  ftime   = df_match$time_to_event_90d,
  fstatus = df_match$outcome,
  cov1    = X_trt,
  failcode = 2,
  cencode  = 0
)

fg_disch_psm <- cmprsk::crr(
  ftime   = df_match$time_to_event_90d,
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
  ftime   = df_match$time_to_event_90d,
  fstatus = df_match$outcome,
  cov1    = X_dr,
  failcode = 2,
  cencode  = 0
)

fg_disch_psm_dr <- cmprsk::crr(
  ftime   = df_match$time_to_event_90d,
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

write.csv(fg_results, file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                                   "_fg_psm_dr_results.csv")),
          row.names = FALSE)
cat("Saved Fine Gray DR model results\n")

### ---- CIF Curves ----

# Fine Gray cumulative incidence estimate by group (death outcome failcode=2)
ci <- cuminc(
  ftime = df_match$time_to_event_90d,
  fstatus = df_match$outcome,
  group = df_match$crrt_high,
  cencode = 0
)

# Convert cuminc object to tidy df for ggplot
tidy_ci <- bind_rows(lapply(names(ci), function(name) {
  if(!grepl(" ", name)) return(NULL) # skip variance elements etc

  data.frame(
    group = sub(" .*", "", name),
    outcome = sub(".* ", "", name),         # should be "2" or "1"
    time = ci[[name]]$time,
    est = ci[[name]]$est
  )
}), .id = NULL)

# Filter only death (failcode 2)
tidy_ci_death <- tidy_ci %>% filter(outcome == "2")

# Plot CIF for death
plot_cif_death <- ggplot(tidy_ci_death, aes(x = time, y = est, color = group)) +
  geom_step(linewidth = 1) +
  labs(
    title = paste0("Cumulative Incidence of Death by CRRT Dose Group"),
    x = "Time from CRRT Initiation (Days)",
    y = "Cumulative Incidence",
    color = "CRRT Dose"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
plot_cif_death

# Save CIF for death
ggsave(file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                    "_cif_death.png")),
       plot_cif_death, width=6, height=4)
cat("Saved CIF for death as PNG\n")

# Filter only discharge (failcode 1)
tidy_ci_discharge <- tidy_ci %>% filter(outcome == "1")

# Plot CIF for discharge
plot_cif_discharge <- ggplot(tidy_ci_discharge, aes(x = time, y = est,
                                                    color = group)) +
  geom_step(linewidth = 1) +
  labs(
    title = paste0("Cumulative Incidence of Discharge by CRRT Dose Group"),
    x = "Time from CRRT Initiation (Days)",
    y = "Cumulative Incidence",
    color = "CRRT Dose"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
plot_cif_discharge

# Save CIF for discharge
ggsave(file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                    "_cif_discharge.png")),
       plot_cif_discharge, width=6, height=4)
cat("Saved CIF for discharge as PNG\n")

# ============================================================= #
# ---- 3. SUPER LEARNER PROPENSITY WEIGHTING (IPTW BRANCH) ----
# ============================================================= #

## ---- A. Making and visualizing IPTW dataset ----

### Super Learner Propensity Weighting ###
set.seed(42)

# Make SL copy of df_tte_bin to avoid overwriting the PSM path
df_tte_sl <- df_tte_bin

### ---- IPTW using Super Learner for propensity score estimation ----
# NOTE: If SL.gam fails due to spaces in race_category level names,
#       either remove SL.gam from library or apply make.names() to data
w_sl <- weightit(
  psm_formula,
  data = df_tte_sl,
  method = "super",
  SL.library = c("SL.glm","SL.gam","SL.randomForest"),
  estimand = "ATE",
  stabilize = TRUE
)

# attach SL weights + PS
df_tte_sl$w  <- w_sl$weights
df_tte_sl$ps <- w_sl$ps

### ---- Overlap Diagnostics ----
plot_sl_overlap <- ggplot(df_tte_sl, aes(x = ps, fill = crrt_high)) +
  geom_histogram(alpha=0.6, position="identity", bins=30) +
  theme_bw(base_size=10) +
  labs(title="Propensity Score Distribution by CRRT Dose Arm (SuperLearner)",
       x="Propensity Score", y="Count", fill="CRRT Dose")
plot_sl_overlap
ggsave(file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                    "_SL_PS_overlap.png")),
       plot_sl_overlap, width=6, height=4)

cat("Saved SL PS overlap plot as PNG\n")

# Weighted balance check BEFORE trimming (Look for Diff.Adj being all <0.1)
# Effective Sampling (aim for close to 1.0)
bal_sl <- cobalt::bal.tab(w_sl, un = TRUE)
bal_df <- as.data.frame(bal_sl$Balance)
bal_sl
write.csv(bal_df, file.path(output_dir,
                            paste0(SITE_NAME, "_", dose_label,
                                   "_SL_IPTW_balance.csv")),
          row.names = TRUE)

ESS <- (sum(df_tte_sl$w)^2) / sum(df_tte_sl$w^2)
ESS_prop <- ESS / nrow(df_tte_sl)
ESS_prop
write_csv(
  tibble(N = nrow(df_tte_sl), ESS = ESS, ESS_prop = ESS_prop),
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                               "_SL_IPTW_ESS.csv"))
)

### ---- Love Plot (SMDs of covariates between groups) ----
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
              paste0(SITE_NAME, "_", dose_label, "_SL_LovePlot.png")),
    width = 8, height = 7, units = "in", res = 300)
print(plot_loveplot_sl)
dev.off()

cat("Saved Love Plot PNG\n")

### ---- TABLE S2 (IPTW) ----
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
    race_category = forcats::fct_recode(
      race_category,
      "Black/African American" = "black or african american",
      "White"                  = "white",
      "Asian"                  = "asian",
      "Native American"        = "american indian or alaska native",
      "Pacific Islander"       = "native hawaiian or other pacific islander",
      "Unknown"                = "unknown",
      "Other"                  = "other"
    ),

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
    paste0(SITE_NAME, "_", dose_label, "_TableS2_IPTW.html")
  )
)

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

## ---- B. Analysis of IPTW ----

### ---- Marginal Structural Cause-Specific Cox Models (with IPTW) ----

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
  as.formula(paste("Surv(time_to_event_90d, outcome == 2) ~", cox_rhs)),
  data = df_tte_sl,
  weights = w,          # IPTW (ATE) from SuperLearner
  robust = TRUE,
  cluster = id
)
summary(fit_cs_death)

# Discharge cause-specific hazard (event code 1)
fit_cs_disch <- coxph(
  as.formula(paste("Surv(time_to_event_90d, outcome == 1) ~", cox_rhs)),
  data = df_tte_sl,
  weights = w,
  robust = TRUE,
  cluster = id
)
summary(fit_cs_disch)

#### ---- Extract IPTW Cox cause-specific HR results (full covariates) ----

extract_iptw_cox <- function(fit, label){
  s <- summary(fit)
  co <- s$coefficients
  ci <- confint(fit)    # CI already on coef scale

  # detect correct z column name
  z_col <- intersect(colnames(co), c("robust z","z"))[1]
  # detect correct p column name
  p_col <- intersect(colnames(co), c("Pr(>|z|)","Robust Pr(>|z|)","p"))[1]

  data.frame(
    outcome = label,
    variable = rownames(co),
    HR = exp(co[, "coef"]),
    HR_lower = exp(ci[,1]),
    HR_upper = exp(ci[,2]),
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
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                               "_IPTW_CauseSpecificCox_FULLresults.csv")),
  row.names = FALSE
)
cat("Saved full IPTW Cause-Specific Cox results.\n")

### ---- Weighted KM Curves ----

# Weighted KM for Death
fit_km_death <- survfit(
  Surv(time_to_event_90d, outcome == 2) ~ crrt_high,
  data = df_tte_sl,
  weights = w
)

plot_km_death <- ggsurvplot(
  fit_km_death,
  data = df_tte_sl,
  fun = "event",  # plots cumulative incidence-like 1 - S(t)
  conf.int = FALSE,
  legend.title = "CRRT Dose",
  legend.labs = levels(df_tte_sl$crrt_high),
  ggtheme = theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5)
    ),
  title = paste0("Weighted Kaplan-Meier: Death"),
  xlab = "Time from CRRT Initiation (Days)",
  ylab = "Cumulative Incidence"
)
plot_km_death

# Save as PNG
ggsave(file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                    "_IPTW_KM_Death.png")),
       plot_km_death$plot, width=6, height=4)

# Weighted KM for Discharge
fit_km_disch <- survfit(
  Surv(time_to_event_90d, outcome == 1) ~ crrt_high,
  data = df_tte_sl,
  weights = w
)

plot_km_disch <- ggsurvplot(
  fit_km_disch,
  data = df_tte_sl,
  fun = "event",
  conf.int = FALSE,
  legend.title = "CRRT Dose",
  legend.labs = levels(df_tte_sl$crrt_high),
  ggtheme = theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5)
    ),
  title = paste0("Weighted Kaplan-Meier: Discharge"),
  xlab = "Time from CRRT Initiation (Days)",
  ylab = "Cumulative Incidence"
)
plot_km_disch
ggsave(file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                    "_IPTW_KM_Discharge.png")),
       plot_km_disch$plot, width=6, height=4)

cat("Weighted KM plots saved as PNGs.\n")

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

  data.frame(
    model = label,
    HR_type = "CauseSpecific HR",
    HR = exp(co[idx,"coef"]),
    HR_lower = exp(ci[idx,1]),
    HR_upper = exp(ci[idx,2]),
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
  iptw_disch_row
)

comparison_table

write.csv(comparison_table,
          file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                       "_ModelComparison_PSMvsIPTW.csv")),
          row.names = FALSE)

cat("Saved model comparison CSV\n")

# ============================================================= #
# ---- 5. FINISH! ----
# ============================================================= #
cat("Code run complete!\n")
