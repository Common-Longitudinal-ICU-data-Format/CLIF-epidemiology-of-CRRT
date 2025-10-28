# ============================================================================
# COMPETING RISK ANALYSIS - CRRT EPIDEMIOLOGY
# ============================================================================
# This script generates competing risk analysis results including:
# - Sample characteristics (Table 1)
# - Cumulative incidence at key time points
# - Fine-Gray model results (subdistribution hazard ratios)
# - Model diagnostics
# - Cumulative incidence plots
# ============================================================================

# Set up R environment
# ----------------------------------------------------------------------------

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

# Load required packages, installing if needed
required_packages <- c("readr", "arrow", "cmprsk", "survival", "jsonlite")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, require, character.only = TRUE)

# Create output directory if it doesn't exist
output_dir <- "output/final"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load configuration
# ----------------------------------------------------------------------------
config_path <- "config/config.json"
if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path)
}

config <- jsonlite::fromJSON(config_path)
SITE_NAME <- config$site_name

cat("Site:", SITE_NAME, "\n")
cat("Timezone:", config$timezone, "\n\n")

# Load data
# ----------------------------------------------------------------------------
data_path <- "output/intermediate/competing_risk_final.parquet"
if (!file.exists(data_path)) {
  stop("File '", data_path, "' not found.")
}
df <- arrow::read_parquet(data_path)

cat("Loaded data:", nrow(df), "rows\n")

# Make sure outcome and covariates exist
required_vars <- c("time_to_event_90d", "outcome",
                   "crrt_dose_ml_kg_hr", "age_at_admission",
                   "sex_category", "sofa_total")
if (!all(required_vars %in% names(df))) {
  stop("One or more required variables are missing from the data frame.")
}

# ============================================================================
# 1. SAMPLE CHARACTERISTICS (Table 1)
# ============================================================================
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Generating Sample Characteristics\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Define model variables for complete case analysis
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

# Calculate summary statistics
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

  # CRRT dose
  crrt_dose_mean = mean(df_complete$crrt_dose_ml_kg_hr, na.rm = TRUE),
  crrt_dose_sd = sd(df_complete$crrt_dose_ml_kg_hr, na.rm = TRUE),
  crrt_dose_median = median(df_complete$crrt_dose_ml_kg_hr, na.rm = TRUE),
  crrt_dose_q25 = quantile(df_complete$crrt_dose_ml_kg_hr, 0.25, na.rm = TRUE),
  crrt_dose_q75 = quantile(df_complete$crrt_dose_ml_kg_hr, 0.75, na.rm = TRUE),

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
          file.path(output_dir, paste0(SITE_NAME, "_sample_characteristics.csv")),
          row.names = FALSE)

cat("Sample characteristics saved.\n\n")

# ============================================================================
# 2. CUMULATIVE INCIDENCE DATA AT KEY TIME POINTS
# ============================================================================
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Calculating Cumulative Incidence at Key Time Points\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Fit cumulative incidence
cr_fit <- cmprsk::cuminc(
  ftime = df_complete$time_to_event_90d,
  fstatus = df_complete$outcome,
  cencode = 0
)

# Debug: Check what's in the cuminc object
cat("Cumulative incidence object structure:\n")
cat("Available outcomes in cuminc object:\n")
print(names(cr_fit))
cat("\n")

# Extract cumulative incidence at specific time points
time_points <- c(7, 14, 21, 28, 60, 90)

# Function to extract CIF estimates at specific times
extract_cif_at_times <- function(cuminc_obj, times, outcome_code, outcome_label) {
  # Get the appropriate element from cuminc object
  # The names in cuminc object are formatted as "group outcome"
  # Try multiple possible naming patterns

  possible_names <- c(
    paste0("0 ", outcome_code),  # "0 1" or "0 2"
    paste0("1 ", outcome_code),  # "1 1" or "1 2" (when no outcome=0)
    as.character(outcome_code)   # "1" or "2" (no grouping)
  )

  cif_data <- NULL
  for (name in possible_names) {
    if (name %in% names(cuminc_obj)) {
      cif_data <- cuminc_obj[[name]]
      cat("Using outcome:", name, "for", outcome_label, "\n")
      break
    }
  }

  # Check if the outcome exists in the data
  if (is.null(cif_data)) {
    cat("Warning: No data for outcome code", outcome_code, "\n")
    cat("Available outcomes in cuminc object:\n")
    print(names(cuminc_obj))
    return(data.frame())
  }

  results <- data.frame()
  for (t in times) {
    # Find the closest time point
    idx <- which.min(abs(cif_data$time - t))

    # Extract estimate
    est <- cif_data$est[idx]

    # Extract variance - handle different possible structures
    if (!is.null(cif_data$var)) {
      var_val <- cif_data$var[idx]
      # Check if var is numeric
      if (is.numeric(var_val) && !is.na(var_val)) {
        se <- sqrt(var_val)
      } else {
        se <- NA
      }
    } else {
      se <- NA
    }

    # Calculate 95% CI
    if (!is.na(se)) {
      ci_lower <- est - 1.96 * se
      ci_upper <- est + 1.96 * se
    } else {
      ci_lower <- NA
      ci_upper <- NA
    }

    # Calculate number at risk and events
    n_at_risk <- nrow(df_complete) - sum(df_complete$time_to_event_90d <= cif_data$time[idx])
    n_events <- sum(df_complete$outcome == outcome_code &
                      df_complete$time_to_event_90d <= cif_data$time[idx])

    results <- rbind(results, data.frame(
      site_name = SITE_NAME,
      time_days = t,
      outcome = outcome_label,
      cif_estimate = est,
      cif_se = ifelse(is.na(se), NA, se),
      cif_lower = ifelse(is.na(ci_lower), NA, ci_lower),
      cif_upper = ifelse(is.na(ci_upper), NA, ci_upper),
      n_at_risk = n_at_risk,
      n_events = n_events,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

# Extract for both outcomes
cif_death <- extract_cif_at_times(cr_fit, time_points, 2, "death")
cif_discharge <- extract_cif_at_times(cr_fit, time_points, 1, "discharge")

# Combine
cif_data <- rbind(cif_death, cif_discharge)

# Save cumulative incidence data
write.csv(cif_data,
          file.path(output_dir, paste0(SITE_NAME, "_cumulative_incidence.csv")),
          row.names = FALSE)

cat("Cumulative incidence data saved.\n\n")

# ============================================================================
# 3. FINE-GRAY MODEL RESULTS
# ============================================================================
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Fitting Fine-Gray Models\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Create model matrix from complete cases
covariates <- model.matrix(~ crrt_dose_ml_kg_hr + age_at_admission +
                             sex_category + sofa_total,
                           data = df_complete)[, -1, drop=FALSE]

cat("Model covariates:\n")
print(colnames(covariates))
cat("\n")

# Function to extract model results
extract_model_results <- function(fg_model, outcome_label, site_name) {
  # Get coefficient summary
  coef_summary <- summary(fg_model)$coef

  # Extract variable names
  var_names <- rownames(coef_summary)

  results <- data.frame(
    site_name = site_name,
    outcome = outcome_label,
    variable = var_names,
    coefficient = coef_summary[, "coef"],
    se = coef_summary[, "se(coef)"],
    z_value = coef_summary[, "z"],
    p_value = coef_summary[, "p-value"],
    shr = coef_summary[, "exp(coef)"],
    shr_lower = coef_summary[, "exp(coef)"] / exp(1.96 * coef_summary[, "se(coef)"]),
    shr_upper = coef_summary[, "exp(coef)"] * exp(1.96 * coef_summary[, "se(coef)"]),
    n_patients = fg_model$n,
    n_events = sum(fg_model$uftime != 0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(results)
}

# Model 1: Death (failcode = 2)
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

death_results <- extract_model_results(fg_death, "death", SITE_NAME)

# Model 2: Discharge Alive (failcode = 1)
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

discharge_results <- extract_model_results(fg_discharge, "discharge", SITE_NAME)

# Combine model results
model_results <- rbind(death_results, discharge_results)

# Save model results
write.csv(model_results,
          file.path(output_dir, paste0(SITE_NAME, "_model_results.csv")),
          row.names = FALSE)

cat("\nModel results saved.\n\n")

# ============================================================================
# 4. MODEL DIAGNOSTICS
# ============================================================================
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Saving Model Diagnostics\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

model_diagnostics <- data.frame(
  site_name = SITE_NAME,
  outcome = c("death", "discharge"),
  n_complete_cases = c(nrow(df_complete), nrow(df_complete)),
  n_excluded = c(nrow(df) - nrow(df_complete), nrow(df) - nrow(df_complete)),
  pct_excluded = c((nrow(df) - nrow(df_complete)) / nrow(df) * 100,
                   (nrow(df) - nrow(df_complete)) / nrow(df) * 100),
  pseudo_loglik = c(fg_death$loglik, fg_discharge$loglik),
  convergence = c(fg_death$converged, fg_discharge$converged),
  stringsAsFactors = FALSE
)

write.csv(model_diagnostics,
          file.path(output_dir, paste0(SITE_NAME, "_model_diagnostics.csv")),
          row.names = FALSE)

cat("Model diagnostics saved.\n\n")

# ============================================================================
# 5. GENERATE CUMULATIVE INCIDENCE PLOT
# ============================================================================
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Generating Cumulative Incidence Plot\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Save plot as PDF
pdf(file.path(output_dir, paste0(SITE_NAME, "_cumulative_incidence_plot.pdf")),
    width = 8, height = 6)

plot(cr_fit,
     xlab = "Days from CRRT Initiation",
     ylab = "Cumulative Incidence",
     main = paste("Competing Risks:", SITE_NAME),
     lty = c(1, 2),
     col = c("blue", "red"),
     lwd = 2)

legend("right",
       legend = c("Discharge Alive", "Death"),
       lty = c(1, 2),
       col = c("blue", "red"),
       lwd = 2,
       bty = "n")

dev.off()

cat("Plot saved.\n\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat(paste(rep("=", 80), collapse=""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

cat("The following files have been saved to", output_dir, ":\n\n")
cat("1.", paste0(SITE_NAME, "_sample_characteristics.csv"), "\n")
cat("   - Sample size, demographics, outcomes distribution\n\n")
cat("2.", paste0(SITE_NAME, "_cumulative_incidence.csv"), "\n")
cat("   - CIF estimates at key time points (7, 14, 21, 28, 60, 90 days)\n\n")
cat("3.", paste0(SITE_NAME, "_model_results.csv"), "\n")
cat("   - Fine-Gray model coefficients, SHRs, CIs, p-values\n\n")
cat("4.", paste0(SITE_NAME, "_model_diagnostics.csv"), "\n")
cat("   - Model fit statistics and diagnostics\n\n")
cat("5.", paste0(SITE_NAME, "_cumulative_incidence_plot.pdf"), "\n")
cat("   - Visual representation of competing risks\n\n")

cat("These aggregated results can be used for multi-site analysis.\n\n")

cat(paste(rep("=", 80), collapse=""), "\n")