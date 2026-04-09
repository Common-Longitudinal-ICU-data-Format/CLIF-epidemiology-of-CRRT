#!/usr/bin/env Rscript
# ===========================================================================
# 05b_dose_response_analysis.R
#
# Continuous dose-response analysis for CRRT dose on 30-day mortality.
# Three complementary approaches:
#   Section 1: Dose decile analysis (descriptive, unadjusted)
#   Section 2: Restricted cubic spline Cox model (regression-adjusted)
#   Section 3: Generalized propensity score (propensity-weighted, causal)
#   Section 4: Combined three-panel figure
#   Section 5: Target trial emulation specification table
#
# Input:  output/intermediate/msm_competing_risk_df.parquet
# Output: output/final/psm_iptw/{SITE}_dose_*.{csv,png,pdf}
#
# Usage:
#   cd /path/to/CLIF-epidemiology-of-CRRT
#   Rscript code/05b_dose_response_analysis.R
# ===========================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("05b: Dose-Response Analysis (Decile + RCS + GPS)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ===================================================================
# 0. SETUP
# ===================================================================

## ---- A. Working directory ----
if (basename(getwd()) == "code") setwd("..")
cat("Working directory:", getwd(), "\n")

## ---- B. Load packages ----
required_packages <- c(
  "tidyverse", "arrow", "survival", "cmprsk", "jsonlite",
  "WeightIt", "cobalt", "survey", "mice", "patchwork"
)
library(splines)  # base R — natural splines for dose-response
new_packages <- required_packages[!(required_packages %in%
                                      installed.packages()[, "Package"])]
if (length(new_packages)) {
  cat("Installing:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "https://cran.r-project.org")
}
lapply(required_packages, require, character.only = TRUE)

## ---- C. Output directory ----
output_dir <- "output/final/psm_iptw"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ---- D. Load configuration ----
config <- jsonlite::fromJSON("config/config.json")
SITE_NAME <- config$site_name
cat("Site:", SITE_NAME, "\n\n")

## ---- E. Load data ----
data_path <- "output/intermediate/msm_competing_risk_df.parquet"
if (!file.exists(data_path)) stop("File not found: ", data_path)
df <- arrow::read_parquet(data_path)
cat("Loaded:", nrow(df), "rows x", ncol(df), "columns\n")

## ---- F. CCI variables ----
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

## ---- G. Race collapse (match script 05) ----
df$race_category <- forcats::fct_collapse(
  df$race_category,
  "White" = "white",
  "Black" = "black or african american",
  other_level = "Other"
)

## ---- H. MICE imputation (same as script 05) ----
model_covariates_full <- c(
  "age_at_admission", "sex_category", "race_category", "weight_kg",
  "lactate_0", "bicarbonate_0", "potassium_0",
  "oxygenation_index_0", "norepinephrine_equivalent_0", "imv_status_0",
  cci_vars
)

required_vars <- c(
  "time_to_event_30d", "outcome", "crrt_dose_ml_kg_hr_0",
  model_covariates_full
)
miss_counts <- colSums(is.na(df[, required_vars]))
if (any(miss_counts > 0)) {
  cat("Imputing", sum(!complete.cases(df[, required_vars])),
      "incomplete rows via MICE (pmm, m=5)...\n")
  vars_with_na <- names(miss_counts[miss_counts > 0])
  predictor_vars <- c("age_at_admission", "sex_category", "lactate_0",
                       "bicarbonate_0", "potassium_0",
                       "norepinephrine_equivalent_0",
                       "imv_status_0", "crrt_dose_ml_kg_hr_0")
  mice_vars <- unique(c(vars_with_na, intersect(predictor_vars, names(df))))
  imp <- mice(df[, mice_vars], m = 5, method = "pmm", seed = 42,
              maxit = 10, printFlag = FALSE)
  df[, mice_vars] <- complete(imp, 1)
  cat("  Imputation complete.\n")
} else {
  cat("No missing values — skipping MICE.\n")
}

## ---- I. Dynamic covariate filter ----
model_covariates <- model_covariates_full[
  sapply(model_covariates_full, function(v) length(unique(df[[v]])) >= 2)
]
cat("Model covariates (", length(model_covariates), "):",
    paste(model_covariates[1:min(5, length(model_covariates))], collapse = ", "),
    "...\n\n")

## ---- J. Constants ----
DOSE_VAR   <- "crrt_dose_ml_kg_hr_0"
DOSE_CUTOFF <- 30
dose_vals  <- df[[DOSE_VAR]]
overall_mort <- mean(df$outcome == 2)

cat("Dose range:", round(min(dose_vals), 1), "-",
    round(max(dose_vals), 1), "mL/kg/hr\n")
cat("Median dose:", round(median(dose_vals), 1), "mL/kg/hr\n")
cat("Overall 30-day mortality:", round(overall_mort * 100, 1), "%\n\n")


# ===================================================================
# 1. DOSE DECILE ANALYSIS (Descriptive / Unadjusted)
# ===================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Section 1: Dose Decile Analysis (Unadjusted)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

df$dose_decile <- ntile(dose_vals, 10)

decile_summary <- df %>%
  group_by(dose_decile) %>%
  summarise(
    n               = n(),
    dose_median     = median(.data[[DOSE_VAR]]),
    dose_q25        = quantile(.data[[DOSE_VAR]], 0.25),
    dose_q75        = quantile(.data[[DOSE_VAR]], 0.75),
    dose_min        = min(.data[[DOSE_VAR]]),
    dose_max        = max(.data[[DOSE_VAR]]),
    n_died          = sum(outcome == 2),
    n_discharged    = sum(outcome == 1),
    n_censored      = sum(outcome == 0),
    mortality_rate  = mean(outcome == 2),
    discharge_rate  = mean(outcome == 1),
    .groups = "drop"
  ) %>%
  mutate(
    # Wilson 95% CI for mortality
    mort_ci_lower = mortality_rate -
      1.96 * sqrt(mortality_rate * (1 - mortality_rate) / n),
    mort_ci_upper = mortality_rate +
      1.96 * sqrt(mortality_rate * (1 - mortality_rate) / n),
    mort_ci_lower = pmax(mort_ci_lower, 0),
    mort_ci_upper = pmin(mort_ci_upper, 1)
  )

print(decile_summary)
write.csv(decile_summary,
          file.path(output_dir, paste0(SITE_NAME, "_dose_decile_mortality.csv")),
          row.names = FALSE)

# Decile plot
p_decile <- ggplot(decile_summary, aes(x = dose_median, y = mortality_rate)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = mort_ci_lower, ymax = mort_ci_upper),
                width = 0.8, color = "steelblue", linewidth = 0.6) +
  geom_line(color = "steelblue", linewidth = 0.5, alpha = 0.5) +
  geom_hline(yintercept = overall_mort, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, NA)) +
  labs(
    title = "30-Day Mortality by CRRT Dose Decile",
    subtitle = "Unadjusted mortality rates with 95% CI",
    x = "Initial CRRT Dose (mL/kg/hr)",
    y = "30-Day Mortality Rate"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(output_dir, paste0(SITE_NAME, "_dose_decile_plot.png")),
       p_decile, width = 6, height = 4, dpi = 300)
ggsave(file.path(output_dir, paste0(SITE_NAME, "_dose_decile_plot.pdf")),
       p_decile, width = 6, height = 4)
cat("Saved dose decile plot.\n\n")


# ===================================================================
# 2. NATURAL SPLINE COX MODEL (Adjusted)
# ===================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Section 2: Natural Spline Cox Model (Adjusted)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Helper: fit spline Cox model & extract dose-response predictions
fit_spline_cox <- function(df, dose_var, covariates, spline_df = 4,
                           dose_grid = NULL, weights = NULL) {
  # Build formula: ns(dose, df) + covariates
  spline_term <- paste0("ns(", dose_var, ", df = ", spline_df, ")")
  rhs <- paste(c(spline_term, covariates), collapse = " + ")
  fml <- as.formula(paste("Surv(time_to_event_30d, outcome == 2) ~", rhs))

  if (is.null(weights)) {
    fit <- coxph(fml, data = df)
  } else {
    fit <- coxph(fml, data = df, weights = weights)
  }

  if (is.null(dose_grid)) {
    dose_grid <- seq(quantile(df[[dose_var]], 0.01),
                     quantile(df[[dose_var]], 0.99),
                     length.out = 200)
  }

  # Build reference row from the first observation, then set covariates to
  # median (numeric) or most common level (factor)
  ref_row <- df[1, , drop = FALSE]
  all_vars <- c(covariates, dose_var)
  for (v in all_vars) {
    if (v %in% names(df)) {
      if (is.factor(df[[v]]) || is.character(df[[v]])) {
        tbl <- sort(table(df[[v]]), decreasing = TRUE)
        ref_row[[v]] <- names(tbl)[1]
      } else {
        ref_row[[v]] <- median(df[[v]], na.rm = TRUE)
      }
    }
  }
  # Ensure factors have correct levels
  for (v in all_vars) {
    if (is.factor(df[[v]])) {
      ref_row[[v]] <- factor(ref_row[[v]], levels = levels(df[[v]]))
    }
  }

  # Replicate for dose grid
  pred_data <- ref_row[rep(1, length(dose_grid)), , drop = FALSE]
  rownames(pred_data) <- NULL
  pred_data[[dose_var]] <- dose_grid

  # Reference point at median dose
  ref_point <- ref_row
  ref_point[[dose_var]] <- median(df[[dose_var]])

  lp_grid <- predict(fit, newdata = pred_data, type = "lp", se.fit = TRUE)
  lp_ref  <- predict(fit, newdata = ref_point, type = "lp")

  log_hr <- lp_grid$fit - lp_ref
  se     <- lp_grid$se.fit

  data.frame(
    dose     = dose_grid,
    log_hr   = log_hr,
    hr       = exp(log_hr),
    hr_lower = exp(log_hr - 1.96 * se),
    hr_upper = exp(log_hr + 1.96 * se)
  )
}

# Fit primary model (df=4, ~5 effective knots)
cat("Fitting natural spline Cox model (df=4)...\n")
dose_grid <- seq(quantile(dose_vals, 0.01),
                 quantile(dose_vals, 0.99),
                 length.out = 200)

spline_fml <- as.formula(
  paste0("Surv(time_to_event_30d, outcome == 2) ~ ns(",
         DOSE_VAR, ", df = 4) + ",
         paste(model_covariates, collapse = " + "))
)
ns_fit <- coxph(spline_fml, data = df)

# Nonlinearity test: compare spline model to linear dose model
linear_fml <- as.formula(
  paste0("Surv(time_to_event_30d, outcome == 2) ~ ",
         DOSE_VAR, " + ",
         paste(model_covariates, collapse = " + "))
)
lin_fit <- coxph(linear_fml, data = df)
nonlin_test <- anova(lin_fit, ns_fit)
cat("\nNonlinearity test (linear vs spline):\n")
print(nonlin_test)
nonlin_p <- nonlin_test[2, "Pr(>|Chi|)"]
cat("Nonlinearity p-value:", format.pval(nonlin_p, digits = 3), "\n")

# Save nonlinearity test
anova_out <- data.frame(
  model = c("Linear dose", "Spline dose (df=4)"),
  loglik = c(logLik(lin_fit), logLik(ns_fit)),
  df = c(length(coef(lin_fit)), length(coef(ns_fit))),
  chi_sq = c(NA, nonlin_test[2, "Chisq"]),
  p_value = c(NA, nonlin_p)
)
write.csv(anova_out,
          file.path(output_dir,
                    paste0(SITE_NAME, "_dose_response_rcs_anova.csv")),
          row.names = FALSE)

# Generate predictions
pred_df <- fit_spline_cox(df, DOSE_VAR, model_covariates,
                          spline_df = 4, dose_grid = dose_grid)
names(pred_df)[1] <- "crrt_dose_ml_kg_hr_0"

write.csv(pred_df,
          file.path(output_dir, paste0(SITE_NAME, "_dose_response_rcs.csv")),
          row.names = FALSE)

# Plot
p_rcs <- ggplot(pred_df, aes(x = crrt_dose_ml_kg_hr_0)) +
  geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = hr), color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = DOSE_CUTOFF, linetype = "dashed",
             color = "red3", linewidth = 0.5) +
  geom_vline(xintercept = median(dose_vals), linetype = "dotted",
             color = "grey60", linewidth = 0.5) +
  scale_y_continuous(trans = "log2",
                     breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2)) +
  coord_cartesian(ylim = c(0.5, 2)) +
  labs(
    title = "Adjusted Hazard Ratio for 30-Day Mortality",
    subtitle = paste0("Covariate-adjusted Cox model with natural spline (nonlinearity p=",
                      format.pval(nonlin_p, digits = 2), ")"),
    x = "Initial CRRT Dose (mL/kg/hr)",
    y = "Hazard Ratio (log scale)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(output_dir, paste0(SITE_NAME, "_dose_response_rcs_plot.png")),
       p_rcs, width = 6, height = 4, dpi = 300)
ggsave(file.path(output_dir, paste0(SITE_NAME, "_dose_response_rcs_plot.pdf")),
       p_rcs, width = 6, height = 4)
cat("Saved adjusted dose-response plot.\n")

# Sensitivity: df=3 and df=5
for (sdf in c(3, 5)) {
  pred_k <- fit_spline_cox(df, DOSE_VAR, model_covariates,
                           spline_df = sdf, dose_grid = dose_grid)
  cat("  Sensitivity (df=", sdf, "): HR range",
      round(min(pred_k$hr), 2), "-", round(max(pred_k$hr), 2), "\n")
}
cat("\n")


# ===================================================================
# 3. GENERALIZED PROPENSITY SCORE (Causal)
# ===================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Section 3: Generalized Propensity Score (CBPS)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Build GPS formula (continuous treatment)
gps_formula <- as.formula(
  paste(DOSE_VAR, "~", paste(model_covariates, collapse = " + "))
)

cat("Estimating generalized propensity scores via CBPS...\n")
W <- tryCatch(
  weightit(gps_formula, data = df, method = "cbps", over = FALSE),
  error = function(e) {
    cat("  CBPS failed:", conditionMessage(e), "\n")
    cat("  Falling back to GLM-based GPS...\n")
    weightit(gps_formula, data = df, method = "glm")
  }
)

# Weight diagnostics
cat("\nGPS weight summary:\n")
w_summary <- data.frame(
  metric = c("min", "Q1", "median", "mean", "Q3", "max",
             "ESS", "max_weight_ratio"),
  value = c(
    min(W$weights), quantile(W$weights, 0.25),
    median(W$weights), mean(W$weights),
    quantile(W$weights, 0.75), max(W$weights),
    ESS(W$weights),
    max(W$weights) / min(W$weights)
  )
)
print(w_summary)
write.csv(w_summary,
          file.path(output_dir, paste0(SITE_NAME, "_gps_weights_summary.csv")),
          row.names = FALSE)

# Covariate balance
bal <- bal.tab(W, stats = c("m", "ks"), un = TRUE)
cat("\nCovariate balance (weighted vs unweighted):\n")
print(bal)

bal_df <- as.data.frame(bal$Balance)
bal_df$variable <- rownames(bal_df)
write.csv(bal_df,
          file.path(output_dir, paste0(SITE_NAME, "_gps_balance.csv")),
          row.names = FALSE)

# Pretty variable names (matches script 05 convention)
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
pretty_names <- c(
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
  cci_labels
)

# Love plot: covariate balance before/after GPS weighting
# CBPS = Covariate Balancing Propensity Score (Imai & Ratkovic 2014).
# Extends propensity scores to continuous treatments by simultaneously
# optimizing covariate balance and treatment model fit. Unlike standard
# GPS (which models treatment | covariates then inverts), CBPS directly
# targets the moment conditions that define balance.
plot_loveplot_gps <- love.plot(
  W,
  stats = "correlations",
  abs = TRUE,
  var.names = pretty_names,
  thresholds = c(correlations = .1),
  var.order = "unadjusted",
  line.size = 0.8,
  point.size = 3,
  sample.names = c("Unweighted", "CBPS-Weighted"),
  title = "Covariate Balance: Generalized Propensity Score (CBPS)",
  subtitle = paste0("Absolute Treatment-Covariate Correlations | ",
                    "Dashed line = 0.1 threshold | ESS = ",
                    round(ESS(W$weights), 0), " / ", nrow(df)),
  grid = TRUE
)

print(plot_loveplot_gps)

png(file.path(output_dir, paste0(SITE_NAME, "_gps_loveplot.png")),
    width = 8, height = 7, units = "in", res = 300)
print(plot_loveplot_gps)
dev.off()

pdf(file.path(output_dir, paste0(SITE_NAME, "_gps_loveplot.pdf")),
    width = 8, height = 7)
print(plot_loveplot_gps)
dev.off()

cat("Saved GPS love plot.\n")

# Weight distribution plot (histogram + by dose quartile)
df$gps_weights_tmp <- W$weights
df$dose_quartile <- factor(ntile(df[[DOSE_VAR]], 4),
                           labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"))

p_wt_hist <- ggplot(df, aes(x = gps_weights_tmp)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "GPS Weight Distribution",
       subtitle = paste0("ESS = ", round(ESS(W$weights), 0),
                         " / ", nrow(df), " | Max weight = ",
                         round(max(W$weights), 1)),
       x = "GPS Weight", y = "Count") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p_wt_dose <- ggplot(df, aes(x = dose_quartile, y = gps_weights_tmp)) +
  geom_boxplot(fill = "lightblue", outlier.color = "red3",
               outlier.size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(title = "GPS Weights by Dose Quartile",
       x = "Dose Quartile", y = "GPS Weight") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p_wt_combined <- p_wt_hist + p_wt_dose +
  plot_layout(widths = c(1, 1))

ggsave(file.path(output_dir, paste0(SITE_NAME, "_gps_weight_diagnostics.png")),
       p_wt_combined, width = 12, height = 5, dpi = 300)
ggsave(file.path(output_dir, paste0(SITE_NAME, "_gps_weight_diagnostics.pdf")),
       p_wt_combined, width = 12, height = 5)
cat("Saved GPS weight diagnostics.\n")

df$gps_weights_tmp <- NULL
df$dose_quartile <- NULL

# Doubly-robust weighted Cox model: GPS weights + covariate adjustment
# Consistent if EITHER the GPS model OR the outcome model is correct.
# This matches the approach in script 05 (IPTW + covariates in Cox).
cat("\nFitting doubly-robust GPS-weighted Cox model (ns, df=3 + covariates)...\n")
df$gps_weights <- W$weights

gps_pred_df <- fit_spline_cox(df, DOSE_VAR, covariates = model_covariates,
                              spline_df = 3, dose_grid = dose_grid,
                              weights = df$gps_weights)
names(gps_pred_df)[1] <- "crrt_dose_ml_kg_hr_0"
cat("Doubly-robust GPS model fit complete.\n")

write.csv(gps_pred_df,
          file.path(output_dir, paste0(SITE_NAME, "_gps_dose_response.csv")),
          row.names = FALSE)

# GPS dose-response plot
p_gps <- ggplot(gps_pred_df, aes(x = crrt_dose_ml_kg_hr_0)) +
  geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper),
              fill = "forestgreen", alpha = 0.2) +
  geom_line(aes(y = hr), color = "forestgreen", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = DOSE_CUTOFF, linetype = "dashed",
             color = "red3", linewidth = 0.5) +
  geom_vline(xintercept = median(dose_vals), linetype = "dotted",
             color = "grey60", linewidth = 0.5) +
  scale_y_continuous(trans = "log2",
                     breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2)) +
  coord_cartesian(ylim = c(0.5, 2)) +
  labs(
    title = "Doubly-Robust Hazard Ratio for 30-Day Mortality",
    subtitle = "CBPS-weighted and covariate-adjusted Cox model",
    x = "Initial CRRT Dose (mL/kg/hr)",
    y = "Hazard Ratio (log scale)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(output_dir, paste0(SITE_NAME, "_gps_dose_response_plot.png")),
       p_gps, width = 6, height = 4, dpi = 300)
ggsave(file.path(output_dir, paste0(SITE_NAME, "_gps_dose_response_plot.pdf")),
       p_gps, width = 6, height = 4)
cat("Saved GPS dose-response plot.\n\n")


# ===================================================================
# 4. COMBINED THREE-PANEL FIGURE
# ===================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Section 4: Combined Figure\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Panel A: Decile (reformat for consistent y-axis)
p_A <- p_decile +
  labs(title = "A. Unadjusted (Dose Deciles)",
       subtitle = NULL) +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

# Panel B: RCS
p_B <- p_rcs +
  labs(title = "B. Regression-Adjusted (Natural Spline Cox)",
       subtitle = NULL) +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

# Panel C: GPS
p_C <- p_gps +
  labs(title = "C. Doubly-Robust (GPS-CBPS + Covariates)",
       subtitle = NULL) +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

# Combine
library(patchwork)
p_combined <- p_A / p_B / p_C +
  plot_annotation(
    title = paste0("Dose-Response Analysis \u2014 ", SITE_NAME),
    subtitle = "Initial CRRT dose vs 30-day mortality across three analytic approaches",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5)
    )
  )

ggsave(file.path(output_dir, paste0(SITE_NAME, "_dose_response_combined.png")),
       p_combined, width = 9, height = 12, dpi = 300)
ggsave(file.path(output_dir, paste0(SITE_NAME, "_dose_response_combined.pdf")),
       p_combined, width = 9, height = 12)
cat("Saved combined three-panel figure.\n\n")


# ===================================================================
# 5. SENSITIVITY: INCLUDING PATIENTS WITH CRRT < 24h
# ===================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Section 5: Sensitivity — Including Patients Excluded by 24h Criterion\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# The primary analysis excludes 503 patients who died or stopped CRRT
# within 24h. This conditions on a post-treatment variable, which could
# bias results if high-dose CRRT causes early death. This sensitivity
# analysis re-fits the linear dose model on the full descriptive cohort.

# Load the full cohort (N=2,136) from tableone_analysis_df
tbl_path <- "output/intermediate/tableone_analysis_df.parquet"
if (!file.exists(tbl_path)) {
  cat("  SKIPPED: tableone_analysis_df.parquet not found.\n\n")
} else {
  df_full <- arrow::read_parquet(tbl_path)
  cat("Full descriptive cohort:", nrow(df_full), "rows\n")

  # Compute time-to-event (days from CRRT initiation to death or censoring)
  df_full$crrt_initiation_time <- as.POSIXct(df_full$crrt_initiation_time)
  df_full$death_dttm <- as.POSIXct(df_full$death_dttm)
  df_full$final_outcome_dttm <- as.POSIXct(df_full$final_outcome_dttm)

  df_full$time_to_event_days <- as.numeric(difftime(
    ifelse(!is.na(df_full$death_dttm), df_full$death_dttm, df_full$final_outcome_dttm),
    df_full$crrt_initiation_time, units = "days"
  ))
  # Cap at 30 days
  df_full$event_30d <- ifelse(
    df_full$death_30d == 1 & df_full$time_to_event_days <= 30, 1, 0
  )
  df_full$time_30d <- pmin(df_full$time_to_event_days, 30)
  df_full$time_30d <- pmax(df_full$time_30d, 0.01)  # avoid zero times

  # Rename covariates to match the primary analysis
  df_full <- df_full %>%
    rename(
      lactate_0 = lactate_baseline,
      bicarbonate_0 = bicarbonate_baseline,
      potassium_0 = potassium_baseline,
      sofa_total_0 = sofa_total
    )

  # Dose from index_crrt (already in tableone_analysis_df)
  df_full$crrt_dose <- df_full$crrt_dose_ml_kg_hr

  # Identify which patients are in the primary vs excluded
  primary_ebs <- df$encounter_block
  df_full$in_primary <- df_full$encounter_block %in% primary_ebs

  n_primary <- sum(df_full$in_primary)
  n_excluded <- sum(!df_full$in_primary)
  n_excluded_with_dose <- sum(!df_full$in_primary & !is.na(df_full$crrt_dose))
  n_excluded_no_dose <- sum(!df_full$in_primary & is.na(df_full$crrt_dose))
  n_sensitivity <- sum(!is.na(df_full$crrt_dose))

  cat("\n--- Cohort Breakdown ---\n")
  cat("  Primary analysis cohort:           ", n_primary, "\n")
  cat("  Excluded by 24h criterion:         ", n_excluded, "\n")
  cat("    - with computable dose:          ", n_excluded_with_dose, "\n")
  cat("    - without dose (excluded again): ", n_excluded_no_dose, "\n")
  cat("  Sensitivity analysis cohort:       ", n_sensitivity, "\n")
  cat("  Patients added back:               ", n_sensitivity - n_primary, "\n\n")

  # Mortality comparison
  mort_primary <- mean(df_full$event_30d[df_full$in_primary], na.rm = TRUE)
  mort_excluded <- mean(df_full$event_30d[!df_full$in_primary & !is.na(df_full$crrt_dose)],
                        na.rm = TRUE)
  cat("30-day mortality:\n")
  cat("  Primary cohort:   ", round(mort_primary * 100, 1), "%\n")
  cat("  Added-back group: ", round(mort_excluded * 100, 1), "%\n\n")

  # Fit linear dose Cox model on full cohort with available covariates
  # (oxygenation_index and NEE not available in tableone_analysis_df)
  df_sens <- df_full %>% filter(!is.na(crrt_dose))

  # Race collapse
  df_sens$race_category <- forcats::fct_collapse(
    df_sens$race_category,
    "White" = "white", "Black" = "black or african american",
    other_level = "Other"
  )

  # Available covariates (subset of primary model)
  sens_covariates <- c(
    "age_at_admission", "sex_category", "race_category",
    "lactate_0", "bicarbonate_0", "potassium_0", "sofa_total_0"
  )
  # Filter to covariates with >=2 unique values
  sens_covariates <- sens_covariates[
    sapply(sens_covariates, function(v) {
      v %in% names(df_sens) && length(unique(df_sens[[v]])) >= 2
    })
  ]

  cat("Sensitivity model covariates:", paste(sens_covariates, collapse = ", "), "\n")
  cat("(Note: oxygenation index, NEE, IMV status, and CCI not available ",
      "for the full cohort in tableone_analysis_df)\n\n")

  sens_fml <- as.formula(
    paste0("Surv(time_30d, event_30d) ~ crrt_dose + ",
           paste(sens_covariates, collapse = " + "))
  )

  sens_fit <- coxph(sens_fml, data = df_sens)
  sens_summary <- summary(sens_fit)
  dose_coef <- sens_summary$coefficients["crrt_dose", ]

  cat("--- Sensitivity Result ---\n")
  cat(sprintf("  Linear dose HR: %.4f (95%% CI: %.4f - %.4f), p = %.3f\n",
              exp(dose_coef["coef"]),
              exp(dose_coef["coef"] - 1.96 * dose_coef["se(coef)"]),
              exp(dose_coef["coef"] + 1.96 * dose_coef["se(coef)"]),
              dose_coef["Pr(>|z|)"]))

  # Compare to primary analysis
  primary_dose_coef <- summary(lin_fit)$coefficients["crrt_dose_ml_kg_hr_0", ]
  cat(sprintf("  Primary analysis HR: %.4f (95%% CI: %.4f - %.4f), p = %.3f\n\n",
              exp(primary_dose_coef["coef"]),
              exp(primary_dose_coef["coef"] - 1.96 * primary_dose_coef["se(coef)"]),
              exp(primary_dose_coef["coef"] + 1.96 * primary_dose_coef["se(coef)"]),
              primary_dose_coef["Pr(>|z|)"]))

  # Save results
  sens_results <- data.frame(
    analysis = c("Primary (N=1,633)", paste0("Sensitivity (N=", n_sensitivity, ")")),
    n = c(n_primary, n_sensitivity),
    n_deaths = c(sum(df$outcome == 2),
                 sum(df_sens$event_30d, na.rm = TRUE)),
    hr = c(exp(primary_dose_coef["coef"]), exp(dose_coef["coef"])),
    hr_lower = c(
      exp(primary_dose_coef["coef"] - 1.96 * primary_dose_coef["se(coef)"]),
      exp(dose_coef["coef"] - 1.96 * dose_coef["se(coef)"])
    ),
    hr_upper = c(
      exp(primary_dose_coef["coef"] + 1.96 * primary_dose_coef["se(coef)"]),
      exp(dose_coef["coef"] + 1.96 * dose_coef["se(coef)"])
    ),
    se_log_hr = c(primary_dose_coef["se(coef)"], dose_coef["se(coef)"]),
    p_value = c(primary_dose_coef["Pr(>|z|)"], dose_coef["Pr(>|z|)"]),
    patients_added_back = c(0, n_sensitivity - n_primary),
    mortality_rate = c(
      round(mort_primary * 100, 1),
      round(mean(df_sens$event_30d, na.rm = TRUE) * 100, 1)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(sens_results,
            file.path(output_dir, paste0(SITE_NAME, "_sensitivity_24h_exclusion.csv")),
            row.names = FALSE)
  cat("Saved sensitivity results.\n\n")
}


# ===================================================================
# DONE
# ===================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("05b complete. Outputs in:", output_dir, "\n")

outputs <- list.files(output_dir,
                      pattern = paste0(SITE_NAME, "_(dose_|gps_|sensitivity_)"),
                      full.names = FALSE)
cat("  Files generated:\n")
for (f in sort(outputs)) cat("    ", f, "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
