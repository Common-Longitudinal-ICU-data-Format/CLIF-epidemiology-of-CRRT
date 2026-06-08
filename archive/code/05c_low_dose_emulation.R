####### Low-Dose CRRT Trial Emulation (NCT06021288) ##################
####### Standard-of-care 25–30 mL/kg/hr vs intervention 10–15 ml/kg/hr
####### Descriptive: apply trial exclusions, histogram CRRT dose,
####### count patients in the 10-15 mL/kg/hr intervention range.

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

# Set working directory to project root using config.json project_root
# Try multiple config paths: relative to code/, relative to cwd, or via --file=
.find_config <- function() {
  candidates <- c("../config/config.json", "config/config.json")
  # Also try relative to script location
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(sub("^--file=", "", file_arg))
    candidates <- c(file.path(script_dir, "..", "config", "config.json"), candidates)
  }
  # RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
    candidates <- c(file.path(script_dir, "..", "config", "config.json"), candidates)
  }
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p))
  }
  stop("Cannot find config/config.json. Please run from the project root or code/ directory.")
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

if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

required_packages <- c("tidyverse", "arrow", "jsonlite", "mice",
                       "survival", "cmprsk", "broom", "survminer", "car")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
lapply(required_packages, require, character.only = TRUE)

## ---- C. Create output directory if it doesn't exist ----
output_dir <- "output/final/low_dose"
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

## ---- E. Load primary data ----
data_path <- "output/intermediate/msm_competing_risk_df.parquet"
if (!file.exists(data_path)) {
  stop("File '", data_path, "' not found.")
}
df <- arrow::read_parquet(data_path)

cat("Loaded msm_competing_risk_df:", nrow(df), "rows\n")

## ---- F. Pull baseline pH from auxiliary tableone parquet ----
# msm_competing_risk_df does not carry pH; tableone_analysis_df does
# (windowed -12h to +3h, last value, forward-filled — same window as the
# other baseline labs).
tableone_path <- "output/intermediate/tableone_analysis_df.parquet"
if (!file.exists(tableone_path)) {
  stop("File '", tableone_path, "' not found (needed for baseline pH).")
}
tableone_cols <- names(arrow::read_parquet(tableone_path, as_data_frame = FALSE))
if (!"ph_arterial_baseline" %in% tableone_cols) {
  stop("Column 'ph_arterial_baseline' not found in tableone_analysis_df.parquet — ",
       "cannot apply NCT06021288 acidosis exclusion at this site.")
}
ph_df <- arrow::read_parquet(
  tableone_path,
  col_select = all_of(c("encounter_block", "ph_arterial_baseline"))
) %>%
  rename(ph_0 = ph_arterial_baseline)

df <- df %>% left_join(ph_df, by = "encounter_block")
cat("Merged baseline pH (n with pH measured:",
    sum(!is.na(df$ph_0)), "/", nrow(df), ")\n\n")

# CCI binary components — same set 05 uses (lines 134-144 in 05).
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

# Full covariate set for the Cox / Fine-Gray models in Section 6 — matches
# 05_PSM_IPTW_CRRT_dose.R lines 277-282 (sofa_total_0 deliberately excluded).
# Column-name fallback: sites that haven't yet rerun step 04 under the
# pf_sf_ratio rename still expose the column as oxygenation_index_0;
# rename in-place so downstream code is uniform.
if (!("pf_sf_ratio_0" %in% names(df)) &&
      "oxygenation_index_0" %in% names(df)) {
  df <- df %>% rename(pf_sf_ratio_0 = oxygenation_index_0)
  cat("Renamed legacy column oxygenation_index_0 -> pf_sf_ratio_0\n")
}

model_covariates <- c(
  "age_at_admission", "sex_category", "race_category", "weight_kg",
  "lactate_0", "bicarbonate_0", "potassium_0",
  "pf_sf_ratio_0", "norepinephrine_equivalent_0", "imv_status_0",
  cci_vars
)

# Required variables for this script (exclusions + outcome + Cox covariates)
required_vars <- unique(c(
  "encounter_block",
  "crrt_dose_median_3h",
  "ph_0",                              # added by tableone merge above
  "time_to_event_30d", "outcome",
  model_covariates,
  "cci_moderate_severe_liver_disease"  # also used as exclusion proxy
))

if (!all(required_vars %in% names(df))) {
  missing <- setdiff(required_vars, names(df))
  stop("Missing variables from data frame: ", paste(missing, collapse = ", "))
}

# Collapse race into 3 categories (matches 05 lines 189-197)
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

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Preparing complete-case dataframe (mirrors 05_PSM_IPTW_CRRT_dose.R)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

## ---- A. Define model variables for complete-case analysis ----
model_vars <- required_vars

# Report missingness before imputation
cat("Missing values per variable:\n")
miss_counts <- colSums(is.na(df[, model_vars]))
print(miss_counts[miss_counts > 0])
cat("\n")

## ---- Pre-MICE missingness visualization ----
# Horizontal bar chart of missingness per model variable so the reviewer can
# see how much of the cohort relies on imputed values for each exclusion.
miss_df <- data.frame(
  variable = names(miss_counts),
  n_missing = as.integer(miss_counts),
  n_total = nrow(df),
  stringsAsFactors = FALSE
)
miss_df$pct_missing <- 100 * miss_df$n_missing / miss_df$n_total
miss_df <- miss_df[order(miss_df$n_missing, decreasing = TRUE), ]
miss_df$variable <- factor(miss_df$variable, levels = rev(miss_df$variable))

miss_plot <- ggplot(miss_df, aes(x = variable, y = n_missing)) +
  geom_col(fill = "#fb801b", color = "black", width = 0.7) +
  geom_text(
    aes(label = sprintf("%d (%.1f%%)", n_missing, pct_missing)),
    hjust = -0.1, size = 3.3
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = paste0("Pre-MICE missingness in NCT06021288 emulation cohort (",
                   SITE_NAME, ")"),
    subtitle = paste0("Denominator n = ", nrow(df),
                      " encounters; variables imputed via PMM (m=5)"),
    x = NULL,
    y = "Missing count"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

ggsave(
  file.path(output_dir, paste0(SITE_NAME, "_low_dose_missingness.png")),
  miss_plot, width = 7, height = 4
)
write.csv(
  miss_df,
  file.path(output_dir, paste0(SITE_NAME, "_low_dose_missingness.csv")),
  row.names = FALSE
)
cat("Saved pre-MICE missingness plot + CSV.\n\n")

# MICE imputation — mirror 05's pattern (pmm, m=5, seed=42, maxit=10).
# This descriptive script only uses the first imputation; Rubin's rules
# pooling is not needed because we are not estimating effects.
n_incomplete <- sum(!complete.cases(df[, model_vars]))
if (n_incomplete > 0) {
  cat("Imputing", n_incomplete, "incomplete rows via MICE (pmm, m=5) …\n")
  library(mice)

  vars_with_na <- names(miss_counts[miss_counts > 0])
  # Predictor set mirrors 05 (lines 224-226) + ph_0 + full covariate panel so
  # the Cox/Fine-Gray fit later is on the same imputed data the exclusions
  # were applied to.
  predictor_vars <- c("age_at_admission", "sex_category", "lactate_0",
                      "bicarbonate_0", "potassium_0",
                      "norepinephrine_equivalent_0",
                      "imv_status_0", "crrt_dose_median_3h", "ph_0",
                      model_covariates)
  mice_vars <- unique(c(vars_with_na, intersect(predictor_vars, model_vars)))
  mice_vars <- mice_vars[mice_vars %in% names(df)]

  imp <- mice(df[, mice_vars], m = 5, method = "pmm",
              seed = 42, maxit = 10, printFlag = FALSE)

  df_complete <- df
  df_complete[, mice_vars] <- complete(imp, 1)
  cat("MICE imputed", sum(miss_counts[miss_counts > 0]),
      "total missing values (first imputation retained for descriptive use)\n")
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
cat("NaN/Inf check complete.\n\n")

# Drop rows still missing any exclusion-relevant variable (defensive — MICE
# normally fills these, but a site without arterial pH at all would leave
# ph_0 entirely NA and we should not silently fail-open the acidosis filter).
exclusion_vars <- c("ph_0", "potassium_0", "cci_moderate_severe_liver_disease",
                    "crrt_dose_median_3h")
df_complete <- df_complete %>% drop_na(all_of(exclusion_vars))

cat("Complete-case analysis sample (post-MICE, post drop_na on exclusion ",
    "variables):", nrow(df_complete), "rows\n\n", sep = "")

if (nrow(df_complete) < 50) {
  stop("Insufficient complete-case sample for emulation (n = ",
       nrow(df_complete), ")")
}

# ================================ #
# ---- 2. NCT06021288 EXCLUSIONS ----
# ================================ #
# Trial:           NCT06021288 (low-dose CRRT vs standard 25-30 mL/kg/hr)
# Inclusion:       Adults, AKI in need of CRRT (already met by our cohort)
# Trial exclusions implemented here:
#   - Severe acidosis:      pH < 7.20
#   - Severe hyperkalemia:  K  > 6.0 mmol/L
#   - Severe liver cirrhosis (Child-Pugh C, proxied by CCI
#     moderate/severe liver disease == 1)
# Not implemented:
#   - Chronic dialysis dependency: already excluded upstream via ESRD ICD codes
#   - CKD with eGFR < 30: no baseline outpatient creatinine available
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Applying NCT06021288 trial exclusion criteria\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Sequential CONSORT-style log
log_rows <- list()
log_step <- function(label, df_in, n_excluded = NA_integer_) {
  log_rows[[length(log_rows) + 1L]] <<- data.frame(
    site_name   = SITE_NAME,
    step        = label,
    n_excluded  = n_excluded,
    n_remaining = nrow(df_in),
    stringsAsFactors = FALSE
  )
}

df_cohort <- df_complete
log_step("Complete-case cohort (post-MICE)", df_cohort)
cat(sprintf("Start (complete-case post-MICE):                %6d\n",
            nrow(df_cohort)))

# 1) Severe acidosis: pH < 7.20
n_before <- nrow(df_cohort)
df_cohort <- df_cohort %>% filter(!(ph_0 < 7.20))
n_excl <- n_before - nrow(df_cohort)
log_step("Exclude severe acidosis (pH < 7.20)", df_cohort, n_excl)
cat(sprintf("- Excluded for pH < 7.20:                       %6d  (remaining: %d)\n",
            n_excl, nrow(df_cohort)))

# 2) Severe hyperkalemia: K > 6.0 mmol/L
n_before <- nrow(df_cohort)
df_cohort <- df_cohort %>% filter(!(potassium_0 > 6.0))
n_excl <- n_before - nrow(df_cohort)
log_step("Exclude severe hyperkalemia (K > 6.0 mmol/L)", df_cohort, n_excl)
cat(sprintf("- Excluded for K > 6.0:                         %6d  (remaining: %d)\n",
            n_excl, nrow(df_cohort)))

# 3) Severe liver cirrhosis proxy: CCI moderate/severe liver == 1 (weight = 3
#    in the Charlson scoring, which we use as a Child-Pugh C proxy per
#    the user's specification).
n_before <- nrow(df_cohort)
df_cohort <- df_cohort %>% filter(cci_moderate_severe_liver_disease == 0)
n_excl <- n_before - nrow(df_cohort)
log_step("Exclude severe liver cirrhosis (CCI moderate/severe liver == 1)",
         df_cohort, n_excl)
cat(sprintf("- Excluded for moderate/severe liver disease:   %6d  (remaining: %d)\n",
            n_excl, nrow(df_cohort)))

cat(sprintf("\nFinal NCT06021288 emulated cohort:              %6d\n\n",
            nrow(df_cohort)))

# Save the CONSORT log
exclusion_log <- do.call(rbind, log_rows)
write.csv(exclusion_log,
          file.path(output_dir,
                    paste0(SITE_NAME, "_low_dose_emulation_exclusion_log.csv")),
          row.names = FALSE)
cat("Saved exclusion log to ",
    paste0(SITE_NAME, "_low_dose_emulation_exclusion_log.csv"), "\n\n",
    sep = "")

## ---- B. CONSORT-style flow diagram ----
# ggplot-rendered CONSORT mirroring the visual style of the existing causal
# CONSORT in 04_build_msm_competing_risk_df.py (rounded boxes, main column +
# excluded-count side branches). Boxes are placed in [0, 1] x [0, 1] space.

consort_main <- data.frame(
  step = c(
    "Start: complete-case cohort\n(post-MICE)",
    "After acidosis exclusion\n(pH ≥ 7.20)",
    "After hyperkalemia exclusion\n(K ≤ 6.0 mmol/L)",
    "Final NCT06021288-emulated cohort\n(no moderate/severe liver disease)"
  ),
  n_remaining = exclusion_log$n_remaining,
  stringsAsFactors = FALSE
)
# One side box per transition between consecutive main boxes
# (3 transitions for 4 main boxes).
consort_excl <- data.frame(
  step = c(
    paste0("Excluded: pH < 7.20\nn = ",
           exclusion_log$n_excluded[2]),
    paste0("Excluded: K > 6.0 mmol/L\nn = ",
           exclusion_log$n_excluded[3]),
    paste0("Excluded: moderate/severe\nliver disease (CCI)\nn = ",
           exclusion_log$n_excluded[4])
  ),
  stringsAsFactors = FALSE
)

n_steps <- nrow(consort_main)
box_w <- 0.42
box_h <- 0.11
x_main <- 0.05
x_excl <- 0.55
y_top <- 0.86
y_bottom <- 0.08
v_spacing <- (y_top - y_bottom) / (n_steps - 1)

main_df <- data.frame(
  xmin = x_main,
  xmax = x_main + box_w,
  ymid = y_top - (seq_len(n_steps) - 1) * v_spacing,
  label = sprintf("%s\nn = %s",
                  consort_main$step,
                  formatC(consort_main$n_remaining, big.mark = ",",
                          format = "d"))
)
main_df$ymin <- main_df$ymid - box_h / 2
main_df$ymax <- main_df$ymid + box_h / 2

n_excl <- nrow(consort_excl)
excl_df <- data.frame(
  xmin = x_excl,
  xmax = x_excl + box_w,
  ymid = head(main_df$ymid, n_excl) - v_spacing / 2,
  label = consort_excl$step,
  stringsAsFactors = FALSE
)
excl_df$ymin <- excl_df$ymid - box_h / 2
excl_df$ymax <- excl_df$ymid + box_h / 2

# Main vertical arrows between consecutive main boxes
main_arrows <- data.frame(
  x = x_main + box_w / 2,
  xend = x_main + box_w / 2,
  y = head(main_df$ymin, -1),
  yend = tail(main_df$ymax, -1)
)

# Side arrows from the midpoint of each transition to the excluded box
side_arrows <- data.frame(
  x = x_main + box_w / 2,
  xend = excl_df$xmin - 0.005,
  y = excl_df$ymid,
  yend = excl_df$ymid
)

consort_plot <- ggplot() +
  geom_rect(
    data = main_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white", color = "black", linewidth = 0.6
  ) +
  geom_text(
    data = main_df,
    aes(x = (xmin + xmax) / 2, y = ymid, label = label),
    size = 3.4, lineheight = 1.0
  ) +
  geom_rect(
    data = excl_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white", color = "black", linewidth = 0.5
  ) +
  geom_text(
    data = excl_df,
    aes(x = (xmin + xmax) / 2, y = ymid, label = label),
    size = 3.2, lineheight = 1.0
  ) +
  geom_segment(
    data = main_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    linewidth = 0.6
  ) +
  geom_segment(
    data = side_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    linewidth = 0.5
  ) +
  annotate(
    "text", x = 0.5, y = 0.97,
    label = paste0("NCT06021288 emulation — CONSORT (", SITE_NAME, ")"),
    fontface = "bold", size = 4.5
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void()

ggsave(
  file.path(output_dir, paste0(SITE_NAME, "_low_dose_consort.png")),
  consort_plot, width = 9, height = 6, dpi = 200
)
cat("Saved CONSORT PNG to ",
    paste0(SITE_NAME, "_low_dose_consort.png"), "\n\n", sep = "")

# ================================ #
# ---- 3. CRRT DOSE HISTOGRAM ----
# ================================ #

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("CRRT dose distribution in NCT06021288-emulated cohort\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Trial dose ranges
intervention_lo  <- 10
intervention_hi  <- 15
standard_lo      <- 25
standard_hi      <- 30

dose_hist <- ggplot(df_cohort, aes(x = crrt_dose_median_3h)) +
  geom_histogram(binwidth = 5, color = "black", fill = "skyblue") +
  geom_vline(xintercept = c(intervention_lo, intervention_hi),
             linetype = "dashed", linewidth = 0.8, color = "#1e417c") +
  geom_vline(xintercept = c(standard_lo, standard_hi),
             linetype = "dashed", linewidth = 0.8, color = "#fb801b") +
  labs(
    title    = paste0("CRRT dose distribution — NCT06021288 emulation (",
                      SITE_NAME, ")"),
    subtitle = paste0("Median first 3h; blue dashed: intervention arm ",
                      intervention_lo, "-", intervention_hi,
                      " mL/kg/hr; orange dashed: standard arm ",
                      standard_lo, "-", standard_hi, " mL/kg/hr"),
    x = "CRRT Dose (mL/kg/hr)",
    y = "Count"
  ) +
  theme_bw(base_size = 12)

print(dose_hist)

ggsave(file.path(output_dir,
                 paste0(SITE_NAME, "_hist_crrt_dose_low_dose_emulation.png")),
       dose_hist, width = 7, height = 4.5)
cat("Saved histogram PNG\n\n")

# ================================ #
# ---- 4. RANGE COUNTS ----
# ================================ #

# Patients in each trial arm range (inclusive boundaries — closer match to
# how a trial would enroll dose-band-targeted patients).
n_total          <- nrow(df_cohort)
n_intervention   <- sum(df_cohort$crrt_dose_median_3h >= intervention_lo &
                          df_cohort$crrt_dose_median_3h <= intervention_hi)
n_standard       <- sum(df_cohort$crrt_dose_median_3h >= standard_lo &
                          df_cohort$crrt_dose_median_3h <= standard_hi)
n_below_inter    <- sum(df_cohort$crrt_dose_median_3h < intervention_lo)
n_between        <- sum(df_cohort$crrt_dose_median_3h > intervention_hi &
                          df_cohort$crrt_dose_median_3h < standard_lo)
n_above_standard <- sum(df_cohort$crrt_dose_median_3h > standard_hi)

cat("Dose-range distribution in emulated cohort (n = ", n_total, "):\n", sep = "")
cat(sprintf("  Below intervention range  (<%d):              %6d  (%.1f%%)\n",
            intervention_lo, n_below_inter, 100 * n_below_inter / n_total))
cat(sprintf("  Intervention range        (%d-%d mL/kg/hr):  %6d  (%.1f%%)\n",
            intervention_lo, intervention_hi,
            n_intervention, 100 * n_intervention / n_total))
cat(sprintf("  Between arms              (>%d, <%d):        %6d  (%.1f%%)\n",
            intervention_hi, standard_lo,
            n_between, 100 * n_between / n_total))
cat(sprintf("  Standard-of-care range    (%d-%d mL/kg/hr):  %6d  (%.1f%%)\n",
            standard_lo, standard_hi,
            n_standard, 100 * n_standard / n_total))
cat(sprintf("  Above standard range      (>%d):              %6d  (%.1f%%)\n",
            standard_hi, n_above_standard,
            100 * n_above_standard / n_total))
cat("\n")

dose_summary <- data.frame(
  site_name             = SITE_NAME,
  n_total               = n_total,
  n_below_intervention  = n_below_inter,
  n_intervention_10_15  = n_intervention,
  n_between_15_25       = n_between,
  n_standard_25_30      = n_standard,
  n_above_standard_30   = n_above_standard,
  pct_intervention_10_15 = round(100 * n_intervention / n_total, 2),
  pct_standard_25_30    = round(100 * n_standard / n_total, 2),
  dose_median           = median(df_cohort$crrt_dose_median_3h, na.rm = TRUE),
  dose_q25              = quantile(df_cohort$crrt_dose_median_3h, 0.25,
                                    na.rm = TRUE),
  dose_q75              = quantile(df_cohort$crrt_dose_median_3h, 0.75,
                                    na.rm = TRUE),
  dose_min              = min(df_cohort$crrt_dose_median_3h, na.rm = TRUE),
  dose_max              = max(df_cohort$crrt_dose_median_3h, na.rm = TRUE),
  stringsAsFactors      = FALSE
)

write.csv(dose_summary,
          file.path(output_dir,
                    paste0(SITE_NAME,
                           "_low_dose_emulation_dose_distribution.csv")),
          row.names = FALSE)
cat("Saved dose-distribution summary to ",
    paste0(SITE_NAME, "_low_dose_emulation_dose_distribution.csv"), "\n\n",
    sep = "")

# ================================ #
# ---- 5. 30-DAY MORTALITY BY ARM ----
# ================================ #
# Descriptive outcomes for the two trial-band groups. No formal modeling
# here — Cox model deferred to a follow-up step. Reports:
#   - n in each arm
#   - 30-day deaths (outcome == 2) with %, Wilson 95% CI
#   - median (IQR) follow-up time among each arm (time_to_event_30d)
#   - median (IQR) follow-up time among deaths only (when n_deaths >= 1)
# The MSM dataframe codes `outcome`: 0 = censored, 1 = discharge, 2 = death,
# matching 05's convention (05_PSM_IPTW_CRRT_dose.R:389-394).

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("30-day mortality by trial arm (10-15 vs 25-30 mL/kg/hr)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Wilson 95% CI for a binomial proportion
wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, center - half), min(1, center + half))
}

fmt_iqr <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return("—")
  q <- quantile(v, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

mort_arms <- list(
  intervention = df_cohort %>%
    filter(crrt_dose_median_3h >= intervention_lo &
             crrt_dose_median_3h <= intervention_hi),
  standard = df_cohort %>%
    filter(crrt_dose_median_3h >= standard_lo &
             crrt_dose_median_3h <= standard_hi)
)

mort_rows <- lapply(names(mort_arms), function(arm) {
  d <- mort_arms[[arm]]
  n <- nrow(d)
  n_deaths <- sum(d$outcome == 2, na.rm = TRUE)
  ci <- wilson_ci(n_deaths, n)
  data.frame(
    site_name = SITE_NAME,
    arm = arm,
    dose_range = if (arm == "intervention")
                   sprintf("%d-%d", intervention_lo, intervention_hi)
                 else sprintf("%d-%d", standard_lo, standard_hi),
    n_total = n,
    n_deaths_30d = n_deaths,
    pct_mortality_30d = if (n > 0) round(100 * n_deaths / n, 1)
                        else NA_real_,
    mortality_ci_lo = if (n > 0) round(100 * ci[1], 1) else NA_real_,
    mortality_ci_hi = if (n > 0) round(100 * ci[2], 1) else NA_real_,
    followup_days_median_iqr = fmt_iqr(d$time_to_event_30d),
    followup_days_deaths_median_iqr =
      fmt_iqr(d$time_to_event_30d[d$outcome == 2]),
    stringsAsFactors = FALSE
  )
})
mortality_summary <- do.call(rbind, mort_rows)

# Pretty-print to console
for (i in seq_len(nrow(mortality_summary))) {
  r <- mortality_summary[i, ]
  cat(sprintf("Arm: %s (%s mL/kg/hr)\n", r$arm, r$dose_range))
  cat(sprintf("  n total:                       %d\n", r$n_total))
  if (r$n_total > 0) {
    cat(sprintf(
      "  30-day deaths:                 %d  (%.1f%%, 95%% CI %.1f-%.1f)\n",
      r$n_deaths_30d, r$pct_mortality_30d,
      r$mortality_ci_lo, r$mortality_ci_hi
    ))
    cat(sprintf("  Follow-up days median [IQR]:   %s\n",
                r$followup_days_median_iqr))
    cat(sprintf("  Among deaths, days median [IQR]: %s\n\n",
                r$followup_days_deaths_median_iqr))
  } else {
    cat("  (no patients in this arm)\n\n")
  }
}

write.csv(
  mortality_summary,
  file.path(output_dir,
            paste0(SITE_NAME, "_low_dose_mortality_by_arm.csv")),
  row.names = FALSE
)
cat("Saved mortality-by-arm summary to ",
    paste0(SITE_NAME, "_low_dose_mortality_by_arm.csv"), "\n\n", sep = "")

# ================================ #
# ---- 6. COX + FINE-GRAY MODELS ----
# ================================ #
# Survival comparison: 10-15 mL/kg/hr (intervention) vs 25-30 mL/kg/hr
# (standard of care, reference). No propensity scores per user direction —
# this is a regression-adjusted competing-risk analysis matching the model
# structure of 05's IPTW Cox + Fine-Gray (lines 1432-1453, 889-931) without
# the IPTW weights.
#
# Four models fit (mirrors 05's parallel csHR + SHR pattern):
#   1. cs-Cox death     — coxph(Surv(t, outcome==2) ~ arm + covariates)
#   2. cs-Cox discharge — coxph(Surv(t, outcome==1) ~ arm + covariates)
#   3. FG death         — cmprsk::crr(failcode=2, cencode=0)
#   4. FG discharge     — cmprsk::crr(failcode=1, cencode=0)
#
# Sample size caveat: with n_intervention=23 and ~11-15 events expected,
# this model is below the conventional EPV>=10 threshold. Treatment HR
# and its CI are reported but should be interpreted as exploratory until
# pooled across sites.

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Cox + Fine-Gray models: 10-15 vs 25-30 mL/kg/hr\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Build two-arm subset with binary treatment indicator (ref = standard).
df_cox <- df_cohort %>%
  mutate(
    crrt_arm = case_when(
      crrt_dose_median_3h >= intervention_lo &
        crrt_dose_median_3h <= intervention_hi ~ "10-15",
      crrt_dose_median_3h >= standard_lo &
        crrt_dose_median_3h <= standard_hi ~ "25-30",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(crrt_arm)) %>%
  mutate(crrt_arm = factor(crrt_arm, levels = c("25-30", "10-15")))

# Drop_na on the full covariate set (matches 05's complete_case_vars approach)
complete_case_vars <- c(model_covariates, "time_to_event_30d", "outcome",
                        "crrt_dose_median_3h")
df_cox <- df_cox %>% drop_na(all_of(complete_case_vars)) %>% droplevels()

cat("Two-arm Cox sample:", nrow(df_cox), "rows\n")
cat("  Standard (25-30, ref): ",
    sum(df_cox$crrt_arm == "25-30"), "\n", sep = "")
cat("  Intervention (10-15):  ",
    sum(df_cox$crrt_arm == "10-15"), "\n", sep = "")
n_deaths_total <- sum(df_cox$outcome == 2, na.rm = TRUE)
n_disch_total  <- sum(df_cox$outcome == 1, na.rm = TRUE)
cat("  Total deaths (event=2):    ", n_deaths_total, "\n", sep = "")
cat("  Total discharges (event=1):", n_disch_total, "\n\n", sep = "")

# Dynamic filter — tighter than 05's <2-level filter to clear small-n
# convergence warnings on sparse binary covariates (e.g., AIDS,
# metastatic_solid_tumor). Rule: keep continuous covariates with any
# variance; for binary/factor, require min cross-tab cell >= 5 with the
# treatment indicator. Drops covariates with perfect/near-perfect
# separation that drive "coefficient may be infinite" warnings.
keep_cox_covariate <- function(v, d) {
  x <- d[[v]]
  nonmiss <- x[!is.na(x)]
  if (length(unique(nonmiss)) < 2) return(FALSE)
  # Treat as binary indicator if values are 0/1 only
  is_binary <- is.numeric(x) && all(nonmiss %in% c(0, 1))
  if (is.factor(x) || is.character(x) || is_binary) {
    tab <- table(d$crrt_arm, x)
    return(min(tab) >= 5)
  }
  TRUE  # continuous, retain
}
keep_flags <- sapply(model_covariates,
                     function(v) keep_cox_covariate(v, df_cox))
cox_covariates <- model_covariates[keep_flags]
dropped_covariates <- model_covariates[!keep_flags]

cat("Covariates after dynamic filter (", length(cox_covariates),
    "/", length(model_covariates), "):\n", sep = "")
cat("  ", paste(cox_covariates, collapse = ", "), "\n", sep = "")
if (length(dropped_covariates) > 0) {
  cat("Dropped (insufficient cell counts in cross-tab with treatment):\n  ",
      paste(dropped_covariates, collapse = ", "), "\n", sep = "")
}
cat("\n")

# EPV diagnostic
n_params_est <- length(model.matrix(
  as.formula(paste("~ crrt_arm +", paste(cox_covariates, collapse = " + "))),
  data = df_cox
)[1, ]) - 1  # subtract intercept
epv_death <- n_deaths_total / n_params_est
epv_disch <- n_disch_total / n_params_est
cat(sprintf(
  "Events-per-variable: %d params; deaths EPV = %.2f; discharge EPV = %.2f\n",
  n_params_est, epv_death, epv_disch
))
if (epv_death < 10) {
  cat("  WARNING: deaths EPV < 10 — treatment HR CIs will be wide; ",
      "interpret as exploratory single-site estimate.\n", sep = "")
}
cat("\n")

cox_rhs <- paste(c("crrt_arm", cox_covariates), collapse = " + ")

## ---- A. Cause-specific Cox (death and discharge) ----

fit_cs_death <- coxph(
  as.formula(paste("Surv(time_to_event_30d, outcome == 2) ~", cox_rhs)),
  data = df_cox
)
fit_cs_disch <- coxph(
  as.formula(paste("Surv(time_to_event_30d, outcome == 1) ~", cox_rhs)),
  data = df_cox
)

extract_cox <- function(fit, label) {
  s <- summary(fit)
  co <- s$coefficients
  ci <- confint(fit)
  data.frame(
    outcome = label,
    variable = rownames(co),
    HR = exp(co[, "coef"]),
    HR_lower = exp(ci[, 1]),
    HR_upper = exp(ci[, 2]),
    se_log_hr = co[, "se(coef)"],
    z = co[, "z"],
    p_value = co[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
cs_results <- bind_rows(
  extract_cox(fit_cs_death, "Death (cause-specific Cox)"),
  extract_cox(fit_cs_disch, "Discharge (cause-specific Cox)")
)
write.csv(
  cs_results,
  file.path(output_dir,
            paste0(SITE_NAME, "_low_dose_CauseSpecificCox_FULL.csv")),
  row.names = FALSE
)

## ---- B. Fine-Gray subdistribution hazards (death and discharge) ----

dr_formula <- as.formula(
  paste("~", paste(c("crrt_arm", cox_covariates), collapse = " + "))
)
X_dr <- model.matrix(dr_formula, data = df_cox)[, -1, drop = FALSE]
colnames(X_dr) <- make.names(colnames(X_dr))

fit_fg_death <- tryCatch(
  cmprsk::crr(
    ftime = df_cox$time_to_event_30d,
    fstatus = df_cox$outcome,
    cov1 = X_dr,
    failcode = 2,
    cencode = 0
  ),
  error = function(e) {
    cat("  FG death model failed to converge:", conditionMessage(e), "\n")
    NULL
  }
)
fit_fg_disch <- tryCatch(
  cmprsk::crr(
    ftime = df_cox$time_to_event_30d,
    fstatus = df_cox$outcome,
    cov1 = X_dr,
    failcode = 1,
    cencode = 0
  ),
  error = function(e) {
    cat("  FG discharge model failed to converge:",
        conditionMessage(e), "\n")
    NULL
  }
)

extract_fg <- function(fit, label) {
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coef
  data.frame(
    outcome = label,
    variable = rownames(s),
    SHR = s[, "exp(coef)"],
    SHR_lower = exp(s[, "coef"] - 1.96 * s[, "se(coef)"]),
    SHR_upper = exp(s[, "coef"] + 1.96 * s[, "se(coef)"]),
    se_log_shr = s[, "se(coef)"],
    z = s[, "z"],
    p_value = s[, "p-value"],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
fg_results <- bind_rows(
  extract_fg(fit_fg_death, "Death (Fine-Gray SHR)"),
  extract_fg(fit_fg_disch, "Discharge (Fine-Gray SHR)")
)
if (!is.null(fg_results) && nrow(fg_results) > 0) {
  write.csv(
    fg_results,
    file.path(output_dir,
              paste0(SITE_NAME, "_low_dose_FineGray_FULL.csv")),
    row.names = FALSE
  )
}

## ---- C. Pull treatment-only summary across all 4 models ----

# The treatment coefficient name differs between coxph (factor expansion) and
# model.matrix (sanitized via make.names) — handle both.
trt_row_cox <- "crrt_arm10-15"
trt_row_fg <- make.names("crrt_arm10-15")

trt_rows <- list()

if (trt_row_cox %in% cs_results$variable) {
  d <- cs_results[cs_results$variable == trt_row_cox, ]
  trt_rows[[length(trt_rows) + 1L]] <- data.frame(
    site_name = SITE_NAME,
    model = "Cause-specific Cox",
    outcome = ifelse(grepl("Death", d$outcome), "Death", "Discharge"),
    estimand = "csHR",
    HR = d$HR, HR_lower = d$HR_lower, HR_upper = d$HR_upper,
    p_value = d$p_value,
    n = nrow(df_cox), n_events = c(n_deaths_total, n_disch_total)[
      match(ifelse(grepl("Death", d$outcome), "Death", "Discharge"),
            c("Death", "Discharge"))],
    stringsAsFactors = FALSE
  )
}

if (!is.null(fg_results) && trt_row_fg %in% fg_results$variable) {
  d <- fg_results[fg_results$variable == trt_row_fg, ]
  trt_rows[[length(trt_rows) + 1L]] <- data.frame(
    site_name = SITE_NAME,
    model = "Fine-Gray",
    outcome = ifelse(grepl("Death", d$outcome), "Death", "Discharge"),
    estimand = "SHR",
    HR = d$SHR, HR_lower = d$SHR_lower, HR_upper = d$SHR_upper,
    p_value = d$p_value,
    n = nrow(df_cox), n_events = c(n_deaths_total, n_disch_total)[
      match(ifelse(grepl("Death", d$outcome), "Death", "Discharge"),
            c("Death", "Discharge"))],
    stringsAsFactors = FALSE
  )}

trt_summary <- do.call(rbind, trt_rows)
write.csv(
  trt_summary,
  file.path(output_dir,
            paste0(SITE_NAME, "_low_dose_treatment_HRs.csv")),
  row.names = FALSE
)

cat("Treatment HRs (10-15 vs 25-30 mL/kg/hr, reference = 25-30):\n")
for (i in seq_len(nrow(trt_summary))) {
  r <- trt_summary[i, ]
  cat(sprintf("  %s (%s) — %s: %.2f (95%% CI %.2f-%.2f, p=%.3f, ",
              r$outcome, r$estimand, r$model,
              r$HR, r$HR_lower, r$HR_upper, r$p_value))
  cat(sprintf("n=%d, events=%d)\n", r$n, r$n_events))
}
cat("\nSaved Cox + Fine-Gray outputs:\n",
    "  ", SITE_NAME, "_low_dose_CauseSpecificCox_FULL.csv\n",
    "  ", SITE_NAME, "_low_dose_FineGray_FULL.csv\n",
    "  ", SITE_NAME, "_low_dose_treatment_HRs.csv\n", sep = "")

# ================================ #
# ---- 7. COX DIAGNOSTICS ----
# ================================ #
# Five-check diagnostic battery for the Section 6 cause-specific Cox
# models. See output/final/low_dose/cox_diagnostics_explainer.md for the
# conceptual walk-through.
#   A. Proportional hazards (Schoenfeld residuals via cox.zph)
#   B. Influential observations (DFBETA on the treatment coefficient)
#   C. Functional form for continuous covariates (Martingale residuals)
#   D. Discrimination (Harrell's C-statistic)
#   E. Multicollinearity (Generalized VIF)
# Pass/fail summary lands in {SITE}_cox_diag_summary.csv.

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Cox model diagnostics\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

## ---- A. Proportional hazards (Schoenfeld) ----
cat("A. Proportional hazards (cox.zph)...\n")

zph_death <- cox.zph(fit_cs_death)
zph_disch <- cox.zph(fit_cs_disch)

zph_to_df <- function(zph, model_label) {
  tab <- as.data.frame(zph$table)
  tab$variable <- rownames(tab)
  tab$model <- model_label
  tab
}
zph_combined <- bind_rows(
  zph_to_df(zph_death, "cs-Cox Death"),
  zph_to_df(zph_disch, "cs-Cox Discharge")
)
write.csv(zph_combined,
          file.path(output_dir,
                    paste0(SITE_NAME, "_cox_diag_zph.csv")),
          row.names = FALSE)

ph_global_death <- zph_death$table["GLOBAL", "p"]
ph_global_disch <- zph_disch$table["GLOBAL", "p"]

# cox.zph collapses 2-level factors back to the base variable name
# (e.g., "crrt_arm" rather than "crrt_arm10-15"); try both forms.
get_trt_ph <- function(zph) {
  candidates <- c("crrt_arm", "crrt_arm10-15")
  hit <- intersect(candidates, rownames(zph$table))
  if (length(hit) == 0) return(NA_real_)
  zph$table[hit[1], "p"]
}
ph_trt_death <- get_trt_ph(zph_death)
ph_trt_disch <- get_trt_ph(zph_disch)

cat(sprintf("  Global PH p (death):     %.3f  %s\n",
            ph_global_death,
            ifelse(ph_global_death < 0.05,
                   "[FAIL - PH violated]", "[pass]")))
cat(sprintf("  Global PH p (discharge): %.3f  %s\n",
            ph_global_disch,
            ifelse(ph_global_disch < 0.05,
                   "[FAIL - PH violated]", "[pass]")))
cat(sprintf("  Treatment PH p (death):     %.3f\n", ph_trt_death))
cat(sprintf("  Treatment PH p (discharge): %.3f\n", ph_trt_disch))

# Plot top-3 worst-p covariates (death model)
worst3_tab <- zph_death$table[
  rownames(zph_death$table) != "GLOBAL", , drop = FALSE
]
worst3_names <- rownames(worst3_tab)[order(worst3_tab[, "p"])][
  seq_len(min(3, nrow(worst3_tab)))
]

if (length(worst3_names) > 0) {
  zph_plots <- survminer::ggcoxzph(zph_death, var = worst3_names)
  png(file.path(output_dir,
                paste0(SITE_NAME, "_cox_diag_zph_death.png")),
      width = 1400, height = 450, res = 120)
  print(zph_plots)
  dev.off()
  cat("  Saved Schoenfeld plot for top-3 worst-p covariates (death)\n")
}

## ---- B. Influential observations (DFBETA on treatment) ----
cat("\nB. Influential observations (DFBETA)...\n")

extract_dfbeta_trt <- function(fit) {
  dfb <- residuals(fit, type = "dfbeta")
  cn <- names(coef(fit))
  if (is.null(dim(dfb))) {
    # Single-coefficient case (shouldn't happen here)
    return(as.numeric(dfb))
  }
  colnames(dfb) <- cn
  trt_idx <- grep("^crrt_arm", cn)[1]
  as.numeric(dfb[, trt_idx])
}

dfbeta_death <- extract_dfbeta_trt(fit_cs_death)
dfbeta_disch <- extract_dfbeta_trt(fit_cs_disch)

dfbeta_threshold <- 2 / sqrt(nrow(df_cox))
n_inf_death <- sum(abs(dfbeta_death) > dfbeta_threshold)
n_inf_disch <- sum(abs(dfbeta_disch) > dfbeta_threshold)
max_abs_dfb_death <- max(abs(dfbeta_death))
max_abs_dfb_disch <- max(abs(dfbeta_disch))

cat(sprintf("  Threshold (2/sqrt(n)): %.4f\n", dfbeta_threshold))
cat(sprintf("  Death: max |DFBETA| = %.4f; n above threshold = %d / %d\n",
            max_abs_dfb_death, n_inf_death, nrow(df_cox)))
cat(sprintf("  Discharge: max |DFBETA| = %.4f; n above threshold = %d / %d\n",
            max_abs_dfb_disch, n_inf_disch, nrow(df_cox)))

dfbeta_df <- bind_rows(
  data.frame(idx = seq_along(dfbeta_death), dfbeta = dfbeta_death,
             arm = df_cox$crrt_arm, model = "Death", died = df_cox$outcome == 2,
             stringsAsFactors = FALSE),
  data.frame(idx = seq_along(dfbeta_disch), dfbeta = dfbeta_disch,
             arm = df_cox$crrt_arm, model = "Discharge",
             died = df_cox$outcome == 2, stringsAsFactors = FALSE)
)

dfbeta_plot <- ggplot(dfbeta_df, aes(x = idx, y = dfbeta, color = arm)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_hline(yintercept = c(-dfbeta_threshold, dfbeta_threshold),
             linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, color = "gray60", linewidth = 0.3) +
  facet_wrap(~ model, scales = "free_y") +
  scale_color_manual(values = c("25-30" = "#fb801b", "10-15" = "#1e417c")) +
  labs(
    title = paste0("DFBETA for treatment coefficient (", SITE_NAME, ")"),
    subtitle = sprintf(paste0("Threshold |DFBETA| > 2/sqrt(n) = %.3f; ",
                              "patients above threshold visibly drive ",
                              "the treatment HR"), dfbeta_threshold),
    x = "Patient index", y = "DFBETA (change in log-HR if dropped)",
    color = "Arm"
  ) +
  theme_bw(base_size = 11)
ggsave(file.path(output_dir,
                 paste0(SITE_NAME, "_cox_diag_dfbeta.png")),
       dfbeta_plot, width = 10, height = 4.5)
cat("  Saved DFBETA plot\n")

## ---- C. Functional form (Martingale residuals) ----
cat("\nC. Functional form (Martingale residuals)...\n")

continuous_covars <- c("age_at_admission", "weight_kg",
                       "lactate_0", "bicarbonate_0", "potassium_0",
                       "pf_sf_ratio_0", "norepinephrine_equivalent_0")
continuous_covars <- intersect(continuous_covars, cox_covariates)

mart_death <- residuals(fit_cs_death, type = "martingale")
mart_panels <- lapply(continuous_covars, function(v) {
  data.frame(
    covariate = v,
    value = df_cox[[v]],
    martingale = mart_death,
    stringsAsFactors = FALSE
  )
})
mart_long <- do.call(rbind, mart_panels)

mart_plot <- ggplot(mart_long, aes(x = value, y = martingale)) +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_smooth(method = "loess", se = FALSE,
              color = "#1e417c", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  facet_wrap(~ covariate, scales = "free_x", ncol = 3) +
  labs(
    title = paste0("Martingale residuals vs continuous covariates (",
                   SITE_NAME, ", cs-Cox Death)"),
    subtitle = paste0("Loess curve should be flat near 0; ",
                      "curvature suggests nonlinear effect, ",
                      "consider rcs() transform"),
    x = "Covariate value", y = "Martingale residual"
  ) +
  theme_bw(base_size = 10)
ggsave(file.path(output_dir,
                 paste0(SITE_NAME, "_cox_diag_martingale.png")),
       mart_plot, width = 10, height = 6)
cat("  Saved martingale residual plot (",
    length(continuous_covars), " continuous covariates)\n", sep = "")

## ---- D. Discrimination (C-statistic) ----
cat("\nD. Discrimination (Harrell's C)...\n")

c_death <- summary(fit_cs_death)$concordance
c_disch <- summary(fit_cs_disch)$concordance

cat(sprintf("  C-statistic (death):     %.3f (SE %.3f) %s\n",
            c_death["C"], c_death["se(C)"],
            ifelse(c_death["C"] > 0.75, "[strong]",
                   ifelse(c_death["C"] > 0.65, "[reasonable]",
                          "[weak]"))))
cat(sprintf("  C-statistic (discharge): %.3f (SE %.3f) %s\n",
            c_disch["C"], c_disch["se(C)"],
            ifelse(c_disch["C"] > 0.75, "[strong]",
                   ifelse(c_disch["C"] > 0.65, "[reasonable]",
                          "[weak]"))))

## ---- E. Multicollinearity (VIF) ----
cat("\nE. Multicollinearity (car::vif)...\n")

compute_max_vif <- function(fit) {
  v <- tryCatch(car::vif(fit), error = function(e) NULL)
  if (is.null(v)) return(NA_real_)
  # Factor covariates produce a 3-col matrix with GVIF^(1/(2*Df))
  if (is.matrix(v)) {
    vif_vals <- v[, "GVIF^(1/(2*Df))"]^2
  } else {
    vif_vals <- v
  }
  max(vif_vals, na.rm = TRUE)
}
max_vif_death <- compute_max_vif(fit_cs_death)
max_vif_disch <- compute_max_vif(fit_cs_disch)

vif_status <- function(v) {
  if (is.na(v)) return("[unavailable]")
  if (v > 10) return("[FAIL >10]")
  if (v > 5) return("[caution >5]")
  "[pass]"
}
cat(sprintf("  Max VIF (death):     %.2f  %s\n",
            max_vif_death, vif_status(max_vif_death)))
cat(sprintf("  Max VIF (discharge): %.2f  %s\n",
            max_vif_disch, vif_status(max_vif_disch)))

# Save per-covariate VIF for the death model (the primary one)
vif_full_death <- tryCatch(car::vif(fit_cs_death), error = function(e) NULL)
if (!is.null(vif_full_death)) {
  if (is.matrix(vif_full_death)) {
    vif_df <- data.frame(
      variable = rownames(vif_full_death),
      GVIF = vif_full_death[, "GVIF"],
      Df = vif_full_death[, "Df"],
      VIF_equivalent = vif_full_death[, "GVIF^(1/(2*Df))"]^2,
      stringsAsFactors = FALSE
    )
  } else {
    vif_df <- data.frame(
      variable = names(vif_full_death),
      VIF = as.numeric(vif_full_death),
      stringsAsFactors = FALSE
    )
  }
  write.csv(vif_df,
            file.path(output_dir,
                      paste0(SITE_NAME, "_cox_diag_vif_death.csv")),
            row.names = FALSE)
}

## ---- F. Pass/fail summary ----
diag_summary <- data.frame(
  site_name = SITE_NAME,
  model = c("cs-Cox Death", "cs-Cox Discharge"),
  n = nrow(df_cox),
  n_events = c(n_deaths_total, n_disch_total),
  n_params = n_params_est,
  epv = c(n_deaths_total / n_params_est,
          n_disch_total / n_params_est),
  c_statistic = c(c_death["C"], c_disch["C"]),
  c_se = c(c_death["se(C)"], c_disch["se(C)"]),
  global_ph_p = c(ph_global_death, ph_global_disch),
  trt_ph_p = c(ph_trt_death, ph_trt_disch),
  max_vif = c(max_vif_death, max_vif_disch),
  max_abs_dfbeta_trt = c(max_abs_dfb_death, max_abs_dfb_disch),
  n_influential_trt = c(n_inf_death, n_inf_disch),
  dfbeta_threshold = dfbeta_threshold,
  stringsAsFactors = FALSE
)
write.csv(diag_summary,
          file.path(output_dir,
                    paste0(SITE_NAME, "_cox_diag_summary.csv")),
          row.names = FALSE)

cat("\nDiagnostic summary saved to ",
    paste0(SITE_NAME, "_cox_diag_summary.csv"), "\n", sep = "")
cat("See cox_diagnostics_explainer.md for interpretation.\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("05c complete.\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
