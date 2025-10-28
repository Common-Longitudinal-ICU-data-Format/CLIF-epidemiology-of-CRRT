# Set up R environment 

# Load required packages, installing if needed
required_packages <- c("readr", "arrow", "cmprsk", "survival")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, require, character.only = TRUE)

# Load data
if (!file.exists("output/intermediate/competing_risk_final.parquet")) {
  stop("File 'output/intermediate/competing_risk_final.parquet' not found.")
}
df <- arrow::read_parquet("output/intermediate/competing_risk_final.parquet")

cat("Loaded data:", nrow(df), "rows\n")

# Make sure outcome and covariates exist
required_vars <- c("time_to_event_90d", "outcome",
                   "crrt_dose_ml_kg_hr", "age_at_admission",
                   "sex_category", "sofa_total")
if (!all(required_vars %in% names(df))) {
  stop("One or more required variables are missing from the data frame.")
}

# ============================================================================
# Cumulative Incidence Curves (All Data)
# ============================================================================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Cumulative Incidence Analysis (All Data)\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

cr_fit <- cmprsk::cuminc(
  ftime = df$time_to_event_90d,
  fstatus = df$outcome,
  cencode = 0
)

plot(cr_fit,
     xlab = "Days from CRRT Initiation",
     ylab = "Cumulative Incidence",
     main = "Competing Risks: Discharge vs Death")

# ============================================================================
# Complete Case Analysis for Fine-Gray Models
# ============================================================================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Fine-Gray Regression Models (Complete Case Analysis)\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Define model variables
model_vars <- c("time_to_event_90d", "outcome", "crrt_dose_ml_kg_hr",
                "age_at_admission", "sex_category", "sofa_total")

# Filter to complete cases
df_complete <- df[complete.cases(df[, model_vars]), ]

cat("Complete cases:", nrow(df_complete), "of", nrow(df), "\n")
cat("Excluded due to missing data:", nrow(df) - nrow(df_complete),
    "(", round((nrow(df) - nrow(df_complete)) / nrow(df) * 100, 1), "%)\n\n")

# Check for sufficient data
if (nrow(df_complete) < 50) {
  stop("Insufficient complete cases for modeling (n = ", nrow(df_complete), ")")
}

# Create model matrix from complete cases
covariates <- model.matrix(~ crrt_dose_ml_kg_hr + age_at_admission +
                             sex_category + sofa_total,
                           data = df_complete)[, -1, drop=FALSE]

cat("Model covariates:\n")
print(colnames(covariates))
cat("\n")

# ============================================================================
# Model 1: Fine-Gray for Death (Outcome = 2)
# ============================================================================
cat(paste(rep("-", 80), collapse=""), "\n")
cat("Model 1: Fine-Gray for Death (Outcome = 2)\n")
cat(paste(rep("-", 80), collapse=""), "\n\n")

fg_death <- cmprsk::crr(
  ftime = df_complete$time_to_event_90d,
  fstatus = df_complete$outcome,
  cov1 = covariates,
  failcode = 2,  # Death
  cencode = 0
)

print(summary(fg_death))

# ============================================================================
# Model 2: Fine-Gray for Discharge Alive (Outcome = 1)
# ============================================================================
cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("Model 2: Fine-Gray for Discharge Alive (Outcome = 1)\n")
cat(paste(rep("-", 80), collapse=""), "\n\n")

fg_discharge <- cmprsk::crr(
  ftime = df_complete$time_to_event_90d,
  fstatus = df_complete$outcome,
  cov1 = covariates,
  failcode = 1,  # Discharged alive
  cencode = 0
)

print(summary(fg_discharge))

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Analysis Complete\n")
cat(paste(rep("=", 80), collapse=""), "\n")