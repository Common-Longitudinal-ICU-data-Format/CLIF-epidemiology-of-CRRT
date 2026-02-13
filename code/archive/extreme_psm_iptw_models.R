# =============================================== #
# SECONDARY ANALYSIS
# PSM and IPTW on Extreme CRRT Doses
# =============================================== #

# ================================ #
# ---- 0. SETUP ----
# ================================ #

## ---- A. Set up R environment ----

### ---- Clear existing environment ----
  # Clear global environment
  rm(list = ls())
  # Close all open plotting devices safely
  while (!is.null(dev.list())) {
    dev.off()
  }
  # Clear console
  cat("\014")
  # Clear warnings
  assign("last.warning", NULL, envir = baseenv())

cat("Environment, plots, and console cleared.\n")

# Set working directory to project root
# Works both interactively (RStudio) and via Rscript
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
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
                       "xgboost","gam","survminer", "survey")
new_packages <- required_packages[!(required_packages %in% installed.packages()
                                    [,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, require, character.only = TRUE)

## ---- C. Create output directory if it doesn't exist ----
output_dir <- "output/final/extreme"
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
data_path <- "output/intermediate/competing_risk_final.parquet"
if (!file.exists(data_path)) {
  stop("File '", data_path, "' not found.")
}
df <- arrow::read_parquet(data_path)

cat("Loaded data:", nrow(df), "rows\n")

# Make sure outcome and covariates exist
required_vars <- c("time_to_event_90d", "outcome",
                   "crrt_dose_ml_kg_hr_full", "age_at_admission",
                   "sex_category", "sofa_total", "lactate_peri_crrt", 
                   "bicarbonate_peri_crrt", "potassium_peri_crrt")
if (!all(required_vars %in% names(df))) {
  stop("One or more required variables are missing from the data frame.")
}

# ================================ #
# ---- 1. DATA CLEANING ----
# ================================ #

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Generating Sample Characteristics\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

## ---- A. Define model variables for complete case analysis ----
model_vars <- c("time_to_event_90d", "outcome", "crrt_dose_ml_kg_hr_full",
                "age_at_admission", "sex_category", "sofa_total", 
                "lactate_peri_crrt", "bicarbonate_peri_crrt", 
                "potassium_peri_crrt")

# Filter to complete cases
df_complete <- df[complete.cases(df[, model_vars]), ]

cat("Complete cases:", nrow(df_complete), "of", nrow(df), "\n")
cat("Excluded due to missing data:", nrow(df) - nrow(df_complete),
    "(", round((nrow(df) - nrow(df_complete)) / nrow(df) * 100, 1), "%)\n\n")

# Check for sufficient data
if (nrow(df_complete) < 50) {
  stop("Insufficient complete cases for modeling (n = ", nrow(df_complete), ")")
}

## ---- B. Extreme CRRT Dose Cutoffs ----
# ======================================== #

dose_low_cutoff_extreme  <- 20
dose_high_cutoff_extreme <- 35

# Label for filenames, tables, plots
dose_label_extreme <- paste0(
  "CRRT_low",  dose_low_cutoff_extreme,
  "_high", dose_high_cutoff_extreme
)

cat(
  "Using EXTREME CRRT dose definition:\n",
  "  Low CRRT dose (<", dose_low_cutoff_extreme, " mL/kg/hr) vs ",
  "High CRRT dose (>", dose_high_cutoff_extreme, " mL/kg/hr)\n\n",
  sep = ""
)

# ======================================== #

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
  
  # SOFA score
  sofa_mean = mean(df_complete$sofa_total, na.rm = TRUE),
  sofa_sd = sd(df_complete$sofa_total, na.rm = TRUE),
  sofa_median = median(df_complete$sofa_total, na.rm = TRUE),
  sofa_q25 = quantile(df_complete$sofa_total, 0.25, na.rm = TRUE),
  sofa_q75 = quantile(df_complete$sofa_total, 0.75, na.rm = TRUE),
  
  # CRRT dose (modality independent)
  crrt_dose_mean = mean(df_complete$crrt_dose_ml_kg_hr_full, na.rm = TRUE),
  crrt_dose_sd = sd(df_complete$crrt_dose_ml_kg_hr_full, na.rm = TRUE),
  crrt_dose_median = median(df_complete$crrt_dose_ml_kg_hr_full, na.rm = TRUE),
  crrt_dose_q25 = quantile(df_complete$crrt_dose_ml_kg_hr_full, 0.25, na.rm = TRUE),
  crrt_dose_q75 = quantile(df_complete$crrt_dose_ml_kg_hr_full, 0.75, na.rm = TRUE),
  
  # Lactate (peri-CRRT)
  lactate_mean   = mean(df_complete$lactate_peri_crrt, na.rm = TRUE),
  lactate_sd     = sd(df_complete$lactate_peri_crrt, na.rm = TRUE),
  lactate_median = median(df_complete$lactate_peri_crrt, na.rm = TRUE),
  lactate_q25    = quantile(df_complete$lactate_peri_crrt, 0.25, na.rm = TRUE),
  lactate_q75    = quantile(df_complete$lactate_peri_crrt, 0.75, na.rm = TRUE),
  
  # Bicarbonate (peri-CRRT)
  bicarbonate_mean   = mean(df_complete$bicarbonate_peri_crrt, na.rm = TRUE),
  bicarbonate_sd     = sd(df_complete$bicarbonate_peri_crrt, na.rm = TRUE),
  bicarbonate_median = median(df_complete$bicarbonate_peri_crrt, na.rm = TRUE),
  bicarbonate_q25    = quantile(df_complete$bicarbonate_peri_crrt, 0.25, na.rm = TRUE),
  bicarbonate_q75    = quantile(df_complete$bicarbonate_peri_crrt, 0.75, na.rm = TRUE),
  
  # Potassium (peri-CRRT)
  potassium_mean   = mean(df_complete$potassium_peri_crrt, na.rm = TRUE),
  potassium_sd     = sd(df_complete$potassium_peri_crrt, na.rm = TRUE),
  potassium_median = median(df_complete$potassium_peri_crrt, na.rm = TRUE),
  potassium_q25    = quantile(df_complete$potassium_peri_crrt, 0.25, na.rm = TRUE),
  potassium_q75    = quantile(df_complete$potassium_peri_crrt, 0.75, na.rm = TRUE),
  
  # Outcomes
  outcome_censored_n = sum(df_complete$outcome == 0, na.rm = TRUE),
  outcome_censored_pct = mean(df_complete$outcome == 0, na.rm = TRUE) * 100,
  outcome_discharge_n = sum(df_complete$outcome == 1, na.rm = TRUE),
  outcome_discharge_pct = mean(df_complete$outcome == 1, na.rm = TRUE) * 100,
  outcome_death_n = sum(df_complete$outcome == 2, na.rm = TRUE),
  outcome_death_pct = mean(df_complete$outcome == 2, na.rm = TRUE) * 100,
  
  stringsAsFactors = FALSE
)

# Save sample characteristics
write.csv(sample_chars,
          file.path(output_dir, paste0(SITE_NAME,"_", dose_label_extreme, 
                                       "_sample_characteristics.csv")),
          row.names = FALSE)

cat("Sample characteristics saved.\n\n")

## ---- D. Histogram CRRT dose ----
crrt_dose_hist_extreme <- ggplot(
  df_complete,
  aes(x = crrt_dose_ml_kg_hr_full)
) +
  geom_histogram(binwidth = 5, color = "black", fill = "skyblue") +
  geom_vline(
    xintercept = dose_low_cutoff_extreme,
    linetype = "dashed",
    linewidth = 1,
    color = "darkred"
  ) +
  geom_vline(
    xintercept = dose_high_cutoff_extreme,
    linetype = "dashed",
    linewidth = 1,
    color = "darkgreen"
  ) +
  labs(
    title = paste0(
      "Distribution of CRRT Dose\n",
      "Low: <", dose_low_cutoff_extreme,
      " vs High: >", dose_high_cutoff_extreme, " mL/kg/hr"
    ),
    x = "CRRT Dose (mL/kg/hr)",
    y = "Count"
  ) +
  theme_bw(base_size = 12)

crrt_dose_hist_extreme

# Save histogram
ggsave(
  file.path(
    output_dir,
    paste0(
      SITE_NAME, "_", dose_label_extreme,
      "_hist_crrt_dose_extreme.png"
    )
  ),
  crrt_dose_hist_extreme,
  width = 6,
  height = 4
)
cat("Saved EXTREME histogram PNG\n\n")

## ---- E. Designate Extreme High vs Low CRRT Dose ----
cat(
  "Defining EXTREME treatment groups: <",
  dose_low_cutoff_extreme, " vs >",
  dose_high_cutoff_extreme, " mL/kg/hr\n",
  sep = ""
)

df_tte_bin_extreme <- df_complete %>%
  mutate(
    # 0 = Low (<20), 1 = High (>35), middle becomes NA
    crrt_high_extreme = case_when(
      crrt_dose_ml_kg_hr_full < dose_low_cutoff_extreme  ~ 0L,
      crrt_dose_ml_kg_hr_full > dose_high_cutoff_extreme ~ 1L,
      TRUE                                               ~ NA_integer_
    )
  ) %>%
  # Drop the middle (20–35) explicitly
  filter(!is.na(crrt_high_extreme)) %>%
  mutate(
    crrt_high_extreme = factor(
      crrt_high_extreme,
      levels = c(0L, 1L),
      labels = c(
        paste0("<", dose_low_cutoff_extreme),
        paste0(">", dose_high_cutoff_extreme)
      )
    )
  ) %>%
  drop_na(
    crrt_high_extreme,
    age_at_admission, sex_category, sofa_total,
    lactate_peri_crrt, bicarbonate_peri_crrt, potassium_peri_crrt,
    time_to_event_90d, outcome
  )

# Distribution of extreme groups
bin_dist_extreme <- df_tte_bin_extreme %>%
  count(crrt_high_extreme) %>%
  mutate(prop = n / sum(n))

bin_dist_extreme

# Save distribution
write.csv(
  bin_dist_extreme,
  file.path(
    output_dir,
    paste0(
      SITE_NAME, "_", dose_label_extreme,
      "_crrt_bin_distribution_extreme.csv"
    )
  ),
  row.names = FALSE
)
cat("Saved EXTREME bin distribution\n\n")

## ---- TABLE 1 ----
# =================================== #

### ---- Modify df to include all three outcomes in one column
df_tte_table1_extreme <- df_tte_bin_extreme %>%
  mutate(
    #### ---- Rename CRRT groups (EXTREME labels)
    crrt_group_extreme = factor(
      crrt_high_extreme,
      levels = c(
        paste0("<", dose_low_cutoff_extreme),
        paste0(">", dose_high_cutoff_extreme)
      ),
      labels = c(
        paste0(
          "Low CRRT dose (<", dose_low_cutoff_extreme, " mL/kg/hr)"
        ),
        paste0(
          "High CRRT dose (>", dose_high_cutoff_extreme, " mL/kg/hr)"
        )
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
  )

### ---- Covariates for Table 1 (EXTREME)
vars_table1_extreme <- c(
  "age_at_admission",
  "sex_category",
  "weight_kg",
  "race_category",
  "sofa_total",
  "creatinine_peri_crrt",
  "lactate_peri_crrt",
  "bicarbonate_peri_crrt",
  "potassium_peri_crrt",
  "crrt_mode_category",
  "crrt_dose_ml_kg_hr_full",
  "outcome_3cat"
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

table1_extreme <- df_tte_table1_extreme %>%
  select(crrt_group_extreme, all_of(vars_table1_extreme)) %>%
  tbl_summary(
    by = crrt_group_extreme,
    
    type = list(
      age_at_admission        ~ "continuous",
      sex_category            ~ "categorical",
      weight_kg               ~ "continuous",
      race_category           ~ "categorical",
      sofa_total              ~ "continuous",
      creatinine_peri_crrt    ~ "continuous",
      lactate_peri_crrt       ~ "continuous",
      bicarbonate_peri_crrt   ~ "continuous",
      potassium_peri_crrt     ~ "continuous",
      crrt_mode_category      ~ "categorical",
      crrt_dose_ml_kg_hr_full ~ "continuous"
    ),
    
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    label = list(
      age_at_admission        ~ "Age at Admission (years)",
      sex_category            ~ "Female (%)",
      weight_kg               ~ "Weight (kg)",
      race_category           ~ "Race",
      sofa_total              ~ "SOFA Score",
      creatinine_peri_crrt    ~ "Creatinine at CRRT Start (mg/dL)",
      lactate_peri_crrt       ~ "Lactate at CRRT Start (mmol/L)",
      bicarbonate_peri_crrt   ~ "Bicarbonate at CRRT Start (mEq/L)",
      potassium_peri_crrt     ~ "Potassium at CRRT Start (mEq/L)",
      crrt_mode_category      ~ "CRRT Modality",
      crrt_dose_ml_kg_hr_full ~ "Initial CRRT Dose (mL/kg/hr)",
      outcome_3cat            ~ "90-day Outcome"
    )
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
  modify_spanning_header(
    c("stat_1", "stat_2") ~ "**CRRT dose group (EXTREME definition)**"
  ) %>%
  modify_caption(
    "**Table 1 (EXTREME). Baseline Characteristics by Extreme CRRT Dose Group**"
  )

# Collapse Sex row using your helper
table1_extreme <- table1_extreme %>%
  collapse_binary_in_gtsummary(
    variable   = "sex_category",
    keep_level = "female",
    label      = "Female (%)"
  )

# Save Table 1 EXTREME as HTML
gt::gtsave(
  as_gt(table1_extreme),
  filename = file.path(
    output_dir,
    paste0(
      SITE_NAME, "_",
      dose_label_extreme,
      "_Table1_unadjusted_extreme.html"
    )
  )
)

# ============================================ #
# ---- 2. PROPENSITY SCORE MATCHING BRANCH ----
# ============================================ #

## ---- A. Making and visualizing PSM dataset ----
### ---- PSM ----
m.out_extreme <- matchit(
  crrt_high_extreme ~ age_at_admission + sex_category + sofa_total +
    lactate_peri_crrt + bicarbonate_peri_crrt + potassium_peri_crrt,
  data   = df_tte_bin_extreme,
  method = "nearest",
  distance = "logit",
  ratio  = 1,
  caliper = 0.2
)

# Balance summary for PSM
summary(m.out_extreme)

# Save PSM balance and counts
psm_bal_extreme <- summary(m.out_extreme)$sum.matched
psm_counts_extreme <- summary(m.out_extreme)$nn %>%
  as.data.frame() %>%
  rownames_to_column(var = "Category")

write.csv(
  psm_bal_extreme,
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme, "_psm_balance_summary_extreme.csv")
  ),
  row.names = FALSE
)

write.csv(
  psm_counts_extreme,
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme, "_psm_counts_summary_extreme.csv")
  ),
  row.names = FALSE
)

cat("Saved EXTREME PSM summaries\n\n")

### ---- Save Matched Dataset ----
df_match_extreme <- match.data(m.out_extreme)

cat("Matched EXTREME dataset created: n = ",
    nrow(df_match_extreme), "\n\n", sep = "")

### ---- Love Plot (EXTREME PSM) ----

pretty_names_loveplot_extreme <- c(
  age_at_admission      = "Age",
  sex_category          = "Sex",
  sofa_total            = "SOFA score",
  lactate_peri_crrt     = "Lactate",
  bicarbonate_peri_crrt = "Bicarbonate",
  potassium_peri_crrt   = "Potassium",
  prop.score            = "Propensity Score",
  distance              = "Propensity Score"
)

plot_loveplot_psm_extreme <- love.plot(
  m.out_extreme,
  stat = c("m", "v"),
  grid = TRUE,
  var.names = pretty_names_loveplot_extreme,
  title = "Covariate Balance: PSM on Extremes",
  threshold = c(m = .25, v = 1.25)
)

plot_loveplot_psm_extreme

ggsave(
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme, "_psm_loveplot_extreme.png")
  ),
  plot_loveplot_psm_extreme,
  width = 6,
  height = 4
)

cat("Saved EXTREME PSM Love plot\n\n")

### ---- TABLE S1 Extreme ----

df_tte_tableS1_extreme <- df_match_extreme %>%
  mutate(
    crrt_group_extreme = factor(
      crrt_high_extreme,
      levels = c(
        paste0("<", dose_low_cutoff_extreme),
        paste0(">", dose_high_cutoff_extreme)
      ),
      labels = c(
        paste0("Low CRRT dose (<", dose_low_cutoff_extreme, " mL/kg/hr)"),
        paste0("High CRRT dose (>", dose_high_cutoff_extreme, " mL/kg/hr)")
      )
    ),
    
    sex_category = factor(sex_category),
    
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
    
    crrt_mode_category = forcats::fct_recode(
      crrt_mode_category,
      "CVVH" = "cvvh",
      "CVVHD" = "cvvhd",
      "CVVHDF" = "cvvhdf",
      "SCUF" = "scuf",
      "AVVH" = "avvh"
    ),
    
    outcome_3cat = factor(
      case_when(
        outcome == 1 ~ "Discharged",
        outcome == 2 ~ "Died",
        TRUE ~ "Censored"
      ),
      levels = c("Discharged", "Died", "Censored")
    )
  )


# Build TABLE S1 (same specification as Table 1)
tableS1_extreme <- df_tte_tableS1_extreme %>%
  select(crrt_group_extreme, all_of(vars_table1_extreme)) %>%
  tbl_summary(
    by = crrt_group_extreme,
    type = list(
      age_at_admission        ~ "continuous",
      sex_category            ~ "categorical",
      weight_kg               ~ "continuous",
      race_category           ~ "categorical",
      sofa_total              ~ "continuous",
      creatinine_peri_crrt    ~ "continuous",
      lactate_peri_crrt       ~ "continuous",
      bicarbonate_peri_crrt   ~ "continuous",
      potassium_peri_crrt     ~ "continuous",
      crrt_mode_category      ~ "categorical",
      crrt_dose_ml_kg_hr_full ~ "continuous"
    ),
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age_at_admission        ~ "Age at Admission (years)",
      sex_category            ~ "Female (%)",
      weight_kg               ~ "Weight (kg)",
      race_category           ~ "Race",
      sofa_total              ~ "SOFA Score",
      creatinine_peri_crrt    ~ "Creatinine at CRRT Start (mg/dL)",
      lactate_peri_crrt       ~ "Lactate at CRRT Start (mmol/L)",
      bicarbonate_peri_crrt   ~ "Bicarbonate at CRRT Start (mEq/L)",
      potassium_peri_crrt     ~ "Potassium at CRRT Start (mEq/L)",
      crrt_mode_category      ~ "CRRT Modality",
      crrt_dose_ml_kg_hr_full ~ "Initial CRRT Dose (mL/kg/hr)",
      outcome_3cat            ~ "90-day Outcome"
    )
  ) %>%
  add_overall() %>%
  add_p(
    test = list(
      race_category ~ "chisq.test.no.correct",
      outcome_3cat ~ "chisq.test.no.correct"
    )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(
    c("stat_1", "stat_2") ~
      "**CRRT dose group (EXTREME matched sample)**"
  ) %>%
  modify_caption(
    "**Table S1 (EXTREME). Baseline Characteristics of the PSM-Matched Cohort**"
  )

# Collapse sex
tableS1_extreme <- tableS1_extreme %>%
  collapse_binary_in_gtsummary(
    variable   = "sex_category",
    keep_level = "female",
    label      = "Female (%)"
  )

# Save
gt::gtsave(
  as_gt(tableS1_extreme),
  filename = file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme,
           "_TableS1_matched_extreme.html")
  )
)

cat("Saved EXTREME Table S1 (PSM)\n\n")


## ---- B. Analysis of PSM ----
### ---- Fine-Gray Analysis on Matched Patients ----

# Matrix for treatment only (reference "<25")
X_trt_extreme <- model.matrix(~ crrt_high_extreme, df_match_extreme)[, -1]

# Death (event = 2)
fg_death_psm_extreme <- cmprsk::crr(
  ftime   = df_match_extreme$time_to_event_90d,
  fstatus = df_match_extreme$outcome,
  cov1    = X_trt_extreme,
  failcode = 2,
  cencode = 0
)

# Discharge (event = 1)
fg_disch_psm_extreme <- cmprsk::crr(
  ftime   = df_match_extreme$time_to_event_90d,
  fstatus = df_match_extreme$outcome,
  cov1    = X_trt_extreme,
  failcode = 1,
  cencode = 0
)

summary(fg_death_psm_extreme)
summary(fg_disch_psm_extreme)

# Add doubly robust model
X_dr_extreme <- model.matrix(
  ~ crrt_high_extreme + age_at_admission + sex_category +
    sofa_total + lactate_peri_crrt + bicarbonate_peri_crrt +
    potassium_peri_crrt,
  data = df_match_extreme
)[, -1]

fg_death_psm_dr_extreme <- cmprsk::crr(
  ftime   = df_match_extreme$time_to_event_90d,
  fstatus = df_match_extreme$outcome,
  cov1    = X_dr_extreme,
  failcode = 2,
  cencode = 0
)

fg_disch_psm_dr_extreme <- cmprsk::crr(
  ftime   = df_match_extreme$time_to_event_90d,
  fstatus = df_match_extreme$outcome,
  cov1    = X_dr_extreme,
  failcode = 1,
  cencode = 0
)

summary(fg_death_psm_dr_extreme)
summary(fg_disch_psm_dr_extreme)

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

fg_results_extreme <- bind_rows(
  extract_fg(fg_death_psm_dr_extreme, "Death"),
  extract_fg(fg_disch_psm_dr_extreme, "Discharge")
)

write.csv(
  fg_results_extreme,
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme,
           "_fg_psm_dr_results_extreme.csv")
  ),
  row.names = FALSE
)

cat("Saved EXTREME FG DR model results\n\n")


### ---- CIF Curves ----

# Fine Gray cumulative incidence estimate by group (death outcome failcode=2)
ci_extreme <- cuminc(
  ftime = df_match_extreme$time_to_event_90d,
  fstatus = df_match_extreme$outcome,
  group = df_match_extreme$crrt_high_extreme,
  cencode = 0
)

# Convert cuminc object to tidy df for ggplot
tidy_ci_extreme <- bind_rows(lapply(names(ci_extreme), function(nm) {
  if (!grepl(" ", nm)) return(NULL)
  data.frame(
    group = sub(" .*", "", nm),
    outcome = sub(".* ", "", nm),
    time = ci_extreme[[nm]]$time,
    est  = ci_extreme[[nm]]$est
  )
}), .id = NULL)

# Filter only death (failcode 2)
tidy_ci_death_extreme <- tidy_ci_extreme %>% filter(outcome == "2")

# Plot CIF for death
plot_cif_death_extreme <- ggplot(
  tidy_ci_death_extreme,
  aes(x = time, y = est, color = group)
) +
  geom_step(linewidth = 1) +
  labs(
    title = paste0(
      "Cumulative Incidence of Death \n",
      "Low CRRT dose (<", dose_low_cutoff_extreme, " mL/kg/hr) vs ",
      "High CRRT dose (>", dose_high_cutoff_extreme, " mL/kg/hr)"
    ),
    x = "Days since CRRT start",
    y = "Cumulative incidence",
    color = "CRRT Dose Group"
  ) +
  theme_bw(base_size = 10)

plot_cif_death_extreme

# Save CIF for death
ggsave(
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label_extreme,
                               "_cif_death_extreme.png")),
  plot_cif_death_extreme, width = 6, height = 4
)

cat("Saved EXTREME CIF for death\n\n")

# Filter only discharge (failcode 1)
tidy_ci_discharge_extreme <- tidy_ci_extreme %>%
  filter(outcome == "1")

# Plot CIF for discharge
plot_cif_discharge_extreme <- ggplot(
  tidy_ci_discharge_extreme,
  aes(x = time, y = est, color = group)
) +
  geom_step(linewidth = 1) +
  labs(
    title = paste0(
      "Cumulative Incidence of Discharge \n",
      "Low CRRT dose (<", dose_low_cutoff_extreme, " mL/kg/hr) vs ",
      "High CRRT dose (>", dose_high_cutoff_extreme, " mL/kg/hr)"
    ),
    x = "Days since CRRT initiation",
    y = "Cumulative incidence",
    color = "CRRT Dose Group"
  ) +
  theme_bw(base_size = 10)

plot_cif_discharge_extreme

# Save CIF for discharge
ggsave(
  file.path(
    output_dir,
    paste0(
      SITE_NAME, "_", dose_label_extreme, "_cif_discharge_extreme.png"
    )
  ),
  plot_cif_discharge_extreme,
  width = 6,
  height = 4
)

cat("Saved EXTREME CIF for discharge\n\n")

# ============================================================= #
# ---- 3. SUPER LEARNER PROPENSITY WEIGHTING (IPTW BRANCH) ---- 
# ============================================================= #

## ---- A. Making and visualizing IPTW dataset ----

### Super Learner Propensity Weighting ###
set.seed(42)

# Make SL copy of df_tte_bin to avoid overwriting the PSM path
df_tte_sl_extreme <- df_tte_bin_extreme

### ---- IPTW using Super Learner for propensity score estimation ----
w_sl_extreme <- weightit(
  crrt_high_extreme ~ age_at_admission + sex_category + sofa_total +
    lactate_peri_crrt + bicarbonate_peri_crrt + potassium_peri_crrt,
  data     = df_tte_sl_extreme,
  method   = "super",
  SL.library = c("SL.glm", "SL.gam", "SL.randomForest", "SL.xgboost"),
  estimand = "ATE",
  stabilize = TRUE
)

# attach SL weights + PS
df_tte_sl_extreme$w  <- w_sl_extreme$weights
df_tte_sl_extreme$ps <- w_sl_extreme$ps

### ---- Overlap Diagnostics ----
plot_sl_overlap_extreme <- ggplot(
  df_tte_sl_extreme,
  aes(x = ps, fill = crrt_high_extreme)
) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_bw(base_size = 10) +
  labs(
    title = 
      "Propensity Score Distribution by Extreme CRRT Dose Arm (SuperLearner)",
    x = "Propensity Score",
    y = "Count",
    fill = "CRRT Dose Group"
  )

plot_sl_overlap_extreme

ggsave(
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme, "_SL_PS_overlap_extreme.png")
  ),
  plot_sl_overlap_extreme,
  width = 6,
  height = 4
)

cat("Saved EXTREME SL overlap plot\n\n")

# Weighted balance check BEFORE trimming (Look for Diff.Adj being all <0.1)
# Effective Sampling (aim for close to 1.0)
bal_sl_extreme <- cobalt::bal.tab(w_sl_extreme, un = TRUE)
bal_sl_extreme_df <- as.data.frame(bal_sl_extreme$Balance)

write.csv(
  bal_sl_extreme_df,
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme, "_SL_IPTW_balance_extreme.csv")
  ),
  row.names = TRUE
)

ESS_extreme <- (sum(df_tte_sl_extreme$w)^2) /
  sum(df_tte_sl_extreme$w^2)
ESS_prop_extreme <- ESS_extreme / nrow(df_tte_sl_extreme)

write_csv(
  tibble(
    N = nrow(df_tte_sl_extreme),
    ESS = ESS_extreme,
    ESS_prop = ESS_prop_extreme
  ),
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme, "_SL_IPTW_ESS_extreme.csv")
  )
)

cat("Saved EXTREME balance + ESS\n\n")

### ---- Love Plot (SMDs of covariates between groups) ----
# Demonstrates balance of chosen covariates

plot_loveplot_sl_extreme <- love.plot(
  bal_sl_extreme,
  stats = "mean.diffs",
  abs = TRUE,
  var.names = pretty_names_loveplot_extreme,
  thresholds = c(m = .1),
  sample.names = c("Unweighted", "Weighted"),
  var.order = "unadjusted",
  title = "Covariate Balance: IPTW via SuperLearner, Extreme",
  subtitle = "Standardized Mean Differences Before and After Weighting",
  grid = TRUE,
)

plot_loveplot_sl_extreme

# Save Love Plot
ggsave(
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme, "_SL_LovePlot_extreme.png")
  ),
  plot_loveplot_sl_extreme,
  width = 6,
  height = 5
)

cat("Saved EXTREME Love plot\n\n")

### ---- TABLE S2 (IPTW) ----
# Like Table 1 but median/IQR/%/p-values are weighted for the pseudopopulation

# Prepare IPTW df (mirror Table 1 preprocessing)
df_tte_tableS2_extreme <- df_tte_sl_extreme %>%
  mutate(
    #### ---- Rename CRRT groups using dose_cutoff
    crrt_group_extreme = factor(
      crrt_high_extreme,
      levels = c(
        paste0("<", dose_low_cutoff_extreme),
        paste0(">", dose_high_cutoff_extreme)
      ),
      labels = c(
        paste0("Low CRRT dose (<", dose_low_cutoff_extreme, " mL/kg/hr)"),
        paste0("High CRRT dose (>", dose_high_cutoff_extreme, " mL/kg/hr)")
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
  # IMPORTANT: Keep only the Table 1 variables
  select(crrt_group_extreme, all_of(vars_table1_extreme), w)

### ---- Survey design using IPTW weights
design_iptw_extreme <- survey::svydesign(
  ids = ~1,
  weights = ~w,
  data = df_tte_tableS2_extreme
)


### ---- Build WEIGHTED Table S2
tableS2_extreme <- tbl_svysummary(
  design_iptw_extreme,        # positional argument (required for your version)
  by = crrt_group_extreme,
  
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{p}%"
  ),
  
  type = list(
    age_at_admission        ~ "continuous",
    sex_category            ~ "categorical",
    weight_kg               ~ "continuous",
    race_category           ~ "categorical",
    sofa_total              ~ "continuous",
    creatinine_peri_crrt    ~ "continuous",
    lactate_peri_crrt       ~ "continuous",
    bicarbonate_peri_crrt   ~ "continuous",
    potassium_peri_crrt     ~ "continuous",
    crrt_mode_category      ~ "categorical",
    crrt_dose_ml_kg_hr_full ~ "continuous"
  ),
  
  label = list(
    age_at_admission        ~ "Age at Admission (years)",
    sex_category            ~ "Female (%)",
    weight_kg               ~ "Weight (kg)",
    race_category           ~ "Race",
    sofa_total              ~ "SOFA Score",
    creatinine_peri_crrt    ~ "Creatinine at CRRT Start (mg/dL)",
    lactate_peri_crrt       ~ "Lactate at CRRT Start (mmol/L)",
    bicarbonate_peri_crrt   ~ "Bicarbonate at CRRT Start (mEq/L)",
    potassium_peri_crrt     ~ "Potassium at CRRT Start (mEq/L)",
    crrt_mode_category      ~ "CRRT Modality",
    crrt_dose_ml_kg_hr_full ~ "Initial CRRT Dose (mL/kg/hr)",
    outcome_3cat            ~ "90-day Outcome",
    w                       ~ "IPTW Weight"
  )
) %>%
  add_overall() %>%
  add_p() %>%                 # weighted p-values are OK
  bold_labels() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(
    c("stat_1", "stat_2") ~ "**CRRT dose group (IPTW weighted)**"
  ) %>%
  modify_caption("**Table S2. Weighted Baseline Characteristics After IPTW (Super Learner PS)**")

# Collapse Sex row
tableS2_extreme <- tableS2_extreme %>%
  collapse_binary_in_gtsummary(
    variable   = "sex_category",
    keep_level = "female",
    label      = "Female (%)"
  )

### ---- Save Table S2
gt::gtsave(
  as_gt(tableS2_extreme),
  filename = file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme, "_TableS2_IPTW_extreme.html")
  )
)
cat("Saved EXTREME Table S2 (IPTW)\n\n")

# ============================================================= #
### ---- OPTIONAL TRIMMING ----
## ---- trimming for non-overlap ---- ONLY DO IF WEIGHTED BAL Diff.Adj > 0.1
#keep <- df_tte_sl_extreme$ps >= .05 & df_tte_sl_extreme$ps <= .95
#df_tte_sl_trim_extreme <- df_tte_sl_extreme[keep,]

## winsorize weights (1st and 99th percentiles)
#q <- quantile(df_tte_sl_trim_extreme$w, c(.01,.99))
#df_tte_sl_trim_extreme$w_trim <- pmin(pmax(df_tte_sl_trim_extreme$w, q[1]), q[2])

#cat("SL IPTW: original n =", nrow(df_tte_sl_extreme), " trimmed n =", 
#nrow(df_tte_sl_trim_extreme), "\n")

# Save trimmed SL dataset
#write.csv(df_tte_sl_trim_extreme, 
#          file.path(output_dir, paste0(SITE_NAME, 
#          "_SL_IPTW_trimmed_extreme_dataset.csv")),
#          row.names = FALSE)
#cat("Saved SL trimmed dataset\n\n")
# ============================================================= #

## ---- B. Analysis of IPTW ----

### ---- Marginal Structural Cause-Specific Cox Models (with IPTW) ----

# References and RowIDs for sandwich SEs
# Dynamically set the reference to the lower-dose group based on cutoff
ref_label_extreme <- paste0("<", dose_low_cutoff_extreme)
df_tte_sl_extreme$crrt_high_extreme <- relevel(
  df_tte_sl_extreme$crrt_high_extreme,
  ref = ref_label_extreme
)

# Robust sandwich SEs (cluster on row id to be safe)
df_tte_sl_extreme$id_extreme <- seq_len(nrow(df_tte_sl_extreme))

# Death cause-specific hazard (event code 2)
fit_cs_death_extreme <- coxph(
  Surv(time_to_event_90d, outcome == 2) ~
    crrt_high_extreme + age_at_admission + sex_category + sofa_total +
    lactate_peri_crrt + bicarbonate_peri_crrt + potassium_peri_crrt,
  data = df_tte_sl_extreme,
  weights = w,
  robust = TRUE,
  cluster = id_extreme
)
summary(fit_cs_death_extreme)

# Discharge cause-specific hazard (event code 1)
fit_cs_disch_extreme <- coxph(
  Surv(time_to_event_90d, outcome == 1) ~
    crrt_high_extreme + age_at_admission + sex_category + sofa_total +
    lactate_peri_crrt + bicarbonate_peri_crrt + potassium_peri_crrt,
  data = df_tte_sl_extreme,
  weights = w,
  robust = TRUE,
  cluster = id_extreme
)

summary(fit_cs_disch_extreme)

#### ---- Extract MSM Cox IPTW cause-specific HR results (full covariates) ----

extract_msm <- function(fit, label){
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
msm_results_extreme <- bind_rows(
  extract_msm(fit_cs_death_extreme,  "Death (EXTREME IPTW MSM Cox)"),
  extract_msm(fit_cs_disch_extreme,  "Discharge (EXTREME IPTW MSM Cox)")
)

write.csv(
  msm_results_extreme,
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme,
           "_IPTW_MSM_CauseSpecificCox_FULLresults_extreme.csv")
  ),
  row.names = FALSE
)

cat("Saved EXTREME MSM Cox results\n\n")

### ---- Weighted KM Curves ----

# Weighted KM for Death
fit_km_death_extreme <- survfit(
  Surv(time_to_event_90d, outcome == 2) ~ crrt_high_extreme,
  data = df_tte_sl_extreme,
  weights = w
)

plot_km_death_extreme <- ggsurvplot(
  fit_km_death_extreme,
  data = df_tte_sl_extreme,
  fun = "event",
  conf.int = FALSE,
  legend.title = "CRRT Dose Cutoffs (mL/kg/hr)",
  legend.labs = levels(df_tte_sl_extreme$crrt_high_extreme),
  ggtheme = theme_bw(base_size = 12),
  title = paste0(
    "Weighted Kaplain-Meier (IPTW) – Death"
  )
)

plot_km_death_extreme

# Save as PNG
ggsave(
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme,
           "_IPTW_KM_Death_extreme.png")
  ),
  plot_km_death_extreme$plot,
  width = 6,
  height = 4
)

# Weighted KM for Discharge 
fit_km_disch_extreme <- survfit(
  Surv(time_to_event_90d, outcome == 1) ~ crrt_high_extreme,
  data = df_tte_sl_extreme,
  weights = w
)

plot_km_disch_extreme <- ggsurvplot(
  fit_km_disch_extreme,
  data = df_tte_sl_extreme,
  fun = "event",
  conf.int = FALSE,
  legend.title = "CRRT Dose Cutoffs (mL/kg/hr)",
  legend.labs = levels(df_tte_sl_extreme$crrt_high_extreme),
  ggtheme = theme_bw(base_size = 12),
  title = paste0(
    "Weighted Kaplain-Meier (IPTW) – Discharge"
  )
)

plot_km_disch_extreme

# Save as PNG
ggsave(
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme,
           "_IPTW_KM_Discharge_extreme.png")
  ),
  plot_km_disch_extreme$plot,
  width = 6,
  height = 4
)

cat("Saved EXTREME weighted KM plots\n\n")

# ============================================================= #
# ----  4. MODEL COMPARISON - PSM FG vs IPTW MSM ----
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

fg_death_row_extreme <- extract_fg_trt(
  fg_death_psm_dr_extreme,
  "PSM FG EXTREME – Death"
)

fg_disch_row_extreme <- extract_fg_trt(
  fg_disch_psm_dr_extreme,
  "PSM FG EXTREME – Discharge"
)

# Extract IPTW MSM treatment rows
extract_msm_trt <- function(fit,label){
  s <- summary(fit)
  co <- s$coefficients
  ci <- confint(fit)
  idx <- grep("crrt_high", rownames(co))
  data.frame(
    model = label,
    HR_type = "CauseSpecific HR",
    HR = exp(co[idx,"coef"]),
    HR_lower = exp(ci[idx,1]),
    HR_upper = exp(ci[idx,2]),
    p_value = co[idx,"Pr(>|z|)"],
    row.names=NULL
  )
}

msm_death_row_extreme <- extract_msm_trt(
  fit_cs_death_extreme,
  "IPTW MSM EXTREME – Death"
)

msm_disch_row_extreme <- extract_msm_trt(
  fit_cs_disch_extreme,
  "IPTW MSM EXTREME – Discharge"
)

comparison_table_extreme <- bind_rows(
  fg_death_row_extreme,
  fg_disch_row_extreme,
  msm_death_row_extreme,
  msm_disch_row_extreme
)

write.csv(
  comparison_table_extreme,
  file.path(
    output_dir,
    paste0(SITE_NAME, "_", dose_label_extreme,
           "_ModelComparison_PSMvsIPTW_extreme.csv")
  ),
  row.names = FALSE
)

cat("Saved EXTREME model comparison\n\n")

# ============================================================= #
# ---- 5. FINISH! ----
# ============================================================= #
cat("Code run complete!\n")
