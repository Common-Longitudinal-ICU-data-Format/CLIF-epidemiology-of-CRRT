###################################################################
############# Propensity Score Matching (PSM) Model ###############
##### Inverse Probability of Treatment Weighting (IPTW) Model #####
###################################################################

################
### 0. SETUP ###
################

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
# ----------------------------------------------------------------------------
required_packages <- c("tidyverse", "readr", "arrow", "cmprsk", "survival", 
                       "jsonlite", "MatchIt", "WeightIt", "cobalt", "EValue", 
                       "SuperLearner", "randomForest","xgboost","gam",
                       "survminer")
new_packages <- required_packages[!(required_packages %in% installed.packages()
                                    [,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, require, character.only = TRUE)

# Create output directory if it doesn't exist
# ----------------------------------------------------------------------------
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
                   "crrt_dose_ml_kg_hr_full", "age_at_admission",
                   "sex_category", "sofa_total", "lactate_peri_crrt", 
                   "bicarbonate_peri_crrt", "potassium_peri_crrt")
if (!all(required_vars %in% names(df))) {
  stop("One or more required variables are missing from the data frame.")
}

########################
### 1. DATA CLEANING ###
########################

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Generating Sample Characteristics\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Define model variables for complete case analysis
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


#### DEBUG THE PATIENTS DROPPED HERE
df_dropped <- df[!complete.cases(df[, model_vars]), ]

nrow(df_dropped)

# Restrict to the variables used in the model
dropped_model <- df_dropped[, model_vars]

# Count of missing by variable
missing_counts <- colSums(is.na(dropped_model))

# Percentage of missing by variable (among dropped patients)
missing_perc <- round(colMeans(is.na(dropped_model)) * 100, 1)

missing_summary <- data.frame(
  variable = names(missing_counts),
  missing_n = as.integer(missing_counts),
  missing_pct_among_dropped = missing_perc,
  row.names = NULL
)

missing_summary


###############################################
########## DEFINING CRRT DOSE CUTOFF ##########
###############################################

# CRRT dose cutoff (mL/kg/hr)
dose_cutoff <- 30

# Generate a label for filenames and titles
dose_label <- paste0("CRRT_", dose_cutoff, "cutoff")

cat("Using CRRT dose cutoff of:", dose_cutoff, "mL/kg/hr\n\n")

###############################################

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
          file.path(output_dir, paste0(SITE_NAME,"_", dose_label, 
                                       "_sample_characteristics.csv")),
          row.names = FALSE)

cat("Sample characteristics saved.\n\n")

# Quick histogram check of CRRT doses to assess bins
crrt_dose_hist <- ggplot(df_complete, aes(x = crrt_dose_ml_kg_hr_full)) +
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

# Designate High vs Low CRRT Dose
cat("Defining treatment groups using cutoff =", dose_cutoff, "mL/kg/hr\n")

df_tte_bin <- df_complete %>% 
  mutate(
    crrt_high = ifelse(crrt_dose_ml_kg_hr_full >= dose_cutoff, 1L, 0L),
    crrt_high = factor(
      crrt_high,
      levels = c(0,1),
      labels = c(paste0("<", dose_cutoff), paste0("≥", dose_cutoff))
    )
  ) %>%
  drop_na(
    crrt_high, age_at_admission, sex_category, sofa_total, 
    lactate_peri_crrt, bicarbonate_peri_crrt, potassium_peri_crrt,
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

############################################
### 2A. PROPENSITY SCORE MATCHING BRANCH ###
############################################

# PSM
m.out <- matchit(
  crrt_high ~ age_at_admission + sex_category + sofa_total +
    lactate_peri_crrt + bicarbonate_peri_crrt + potassium_peri_crrt,
  data = df_tte_bin,
  method = "nearest",
  distance = "logit",
  ratio = 1,
  caliper = 0.2
)

# Balance Summary
summary(m.out)
# Save PSM summary
psm_sum <- broom::tidy(summary(m.out)$sum.matched)
write.csv(psm_sum, file.path(output_dir, paste0(SITE_NAME, "_", dose_label, 
                                                "_psm_balance_summary.csv")), 
          row.names = FALSE)
cat("Saved PSM summary\n")

# Save Matched Dataset
df_match <- match.data(m.out)  

# Pretty variable names for Love plot
pretty_names_loveplot <- c(
  age_at_admission      = "Age",
  sex_category          = "Male sex",
  sofa_total            = "SOFA score",
  lactate_peri_crrt     = "Lactate",
  bicarbonate_peri_crrt = "Bicarbonate",
  potassium_peri_crrt   = "Potassium",
  
  # For IPTW SL (ps column)
  prop.score            = "Propensity Score",
  
  # For PSM (MatchIt internal name)
  distance              = "Propensity Score"
)

# Love Plot (SMD + variance ratio)
plot_loveplot_psm <- love.plot(
  m.out,
  stat = c("m", "v"),
  grid = TRUE,
  var.names = pretty_names_loveplot,
  title = "Covariate Balance: PSM",
  threshold = c(m = .25, v = 1.25)
)

plot_loveplot_psm

# Save Love plot
png(file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_psm_loveplot.png")),
    width = 6, height = 4, units = "in", res = 300)
print(plot_loveplot_psm)
dev.off()

cat("Saved Love plot for PSM as PNG\n")


### Fine=Gray Analysis on Matched Patients ###

# Matrix for treatment only (reference "<25")
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

# Add doubly robust model
X_dr <- model.matrix(~ crrt_high + age_at_admission + sex_category + sofa_total +
                       lactate_peri_crrt + bicarbonate_peri_crrt + potassium_peri_crrt,
                     data = df_match)[, -1, drop = FALSE]
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

### CIF Curves ###

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
    title = paste0("Cumulative Incidence of Death by CRRT Dose Group 
                   (Cutoff = ", dose_cutoff, " mL/kg/hr)"),
    x = "Days since CRRT initiation",
    y = "Cumulative incidence",
    color = "CRRT Dose"
  ) +
  theme_bw(base_size = 10)
plot_cif_death
# Save CIF for death
ggsave(file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                    "_cif_death.png")), 
       plot_cif_death, width=6, height=4)
cat("Saved CIF for death as PNG\n")

# Filter only discharge (failcode 1)
tidy_ci_discharge <- tidy_ci %>% filter(outcome == "1")

# Plot CIF for discharge
plot_cif_discharge <- ggplot(tidy_ci_discharge, aes(x = time, y = est, color = group)) +
  geom_step(linewidth = 1) +
  labs(
    title = paste0("Cumulative Incidence of Discharge by CRRT Dose Group 
                   (Cutoff = ", dose_cutoff, " mL/kg/hr)"),
    x = "Days since CRRT initiation",
    y = "Cumulative incidence",
    color = "CRRT Dose"
  ) +
  theme_bw(base_size = 10)
plot_cif_discharge
# Save CIF for discharge
ggsave(file.path(output_dir, paste0(SITE_NAME, "_", dose_label,
                                    "_cif_discharge.png")), 
       plot_cif_discharge, width=6, height=4)
cat("Saved CIF for discharge as PNG\n")

############################################################
### 2B. SUPER LEARNER PROPENSITY WEIGHTING (IPTW BRANCH) ###
############################################################

### Super Learner Propensity Weighting ###
set.seed(42)

# Make SL copy of df_tte_bin to avoid overwriting the PSM path
df_tte_sl <- df_tte_bin

# IPTW using Super Learner for propensity score estimation
w_sl <- weightit(
  crrt_high ~ age_at_admission + sex_category + sofa_total +
    lactate_peri_crrt + bicarbonate_peri_crrt + potassium_peri_crrt,
  data = df_tte_sl,
  method = "super",
  SL.library = c("SL.glm","SL.gam","SL.randomForest","SL.xgboost"),
  estimand = "ATE",
  stabilize = TRUE
)

# attach SL weights + PS
df_tte_sl$w  <- w_sl$weights
df_tte_sl$ps <- w_sl$ps

# ---- overlap diagnostics ----
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

# Love Plot (SMDs of covariates between groups)
# Demonstrates balance of chosen covariates

plot_loveplot_sl <- love.plot(
  bal_sl,
  stats = "mean.diffs",
  abs = TRUE,
  var.names = pretty_names_loveplot,       # <-- rename here
  thresholds = c(m = .1),
  var.order = "unadjusted",
  line.size = 0.8,
  point.size = 3,
  sample.names = c("Unweighted", "Weighted"),
  title = "Covariate Balance: IPTW via SuperLearner",
  subtitle = "Standardized Mean Differences Before and After Weighting",
  grid = TRUE
)

plot_loveplot_sl

# Save Love Plot
png(file.path(output_dir, paste0(SITE_NAME, "_", dose_label, "_SL_LovePlot.png")),
    width = 6, height = 5, units = "in", res = 300)
print(plot_loveplot_sl)
dev.off()

cat("Saved Love Plot PNG\n")

################## OPTIONAL TRIMMING ######################
## ---- trimming for non-overlap ---- ONLY DO IF WEIGHTED BAL Diff.Adj > 0.1
#keep <- df_tte_sl$ps >= .05 & df_tte_sl$ps <= .95
#df_tte_sl_trim <- df_tte_sl[keep,]

## winsorize weights (1st and 99th percentiles)
#q <- quantile(df_tte_sl_trim$w, c(.01,.99))
#df_tte_sl_trim$w_trim <- pmin(pmax(df_tte_sl_trim$w, q[1]), q[2])

#cat("SL IPTW: original n =", nrow(df_tte_sl), " trimmed n =", nrow(df_tte_sl_trim), "\n")

# Save trimmed SL dataset
#write.csv(df_tte_sl_trim, 
#          file.path(output_dir, paste0(SITE_NAME, "_SL_IPTW_trimmed_dataset.csv")),
#          row.names = FALSE)
#cat("Saved SL trimmed dataset\n\n")
###########################################################

### Marginal Structural Cause-Specific Cox Models (with IPTW) ###

# References and RowIDs for sandwich SEs
# Dynamically set the reference to the lower-dose group based on cutoff
ref_label <- paste0("<", dose_cutoff)
df_tte_sl$crrt_high <- relevel(df_tte_sl$crrt_high, ref = ref_label)

# Robust sandwich SEs (cluster on row id to be safe)
df_tte_sl$id <- seq_len(nrow(df_tte_sl))

# Death cause-specific hazard (event code 2)
fit_cs_death <- coxph(
  Surv(time_to_event_90d, outcome == 2) ~ crrt_high + age_at_admission +
    sex_category + sofa_total + lactate_peri_crrt + bicarbonate_peri_crrt +
    potassium_peri_crrt,
  data = df_tte_sl,
  weights = w,          # IPTW (ATE) from SuperLearner
  robust = TRUE,
  cluster = id
)
summary(fit_cs_death)

# Discharge cause-specific hazard (event code 1)
fit_cs_disch <- coxph(
  Surv(time_to_event_90d, outcome == 1) ~ crrt_high + age_at_admission +
    sex_category + sofa_total + lactate_peri_crrt + bicarbonate_peri_crrt +
    potassium_peri_crrt,
  data = df_tte_sl,
  weights = w,
  robust = TRUE,
  cluster = id
)
summary(fit_cs_disch)

### Extract MSM Cox IPTW cause-specific HR results (full covariates) ###

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
msm_results <- bind_rows(
  extract_msm(fit_cs_death,  "Death (IPTW MSM Cause-Specific Cox)"),
  extract_msm(fit_cs_disch,  "Discharge (IPTW MSM Cause-Specific Cox)")
)
write.csv(
  msm_results,
  file.path(output_dir, paste0(SITE_NAME, "_", dose_label, 
                               "_IPTW_MSM_CauseSpecificCox_FULLresults.csv")),
  row.names = FALSE
)
cat("Saved full MSM IPTW Cause-Specific Cox results.\n")

### Weighted KM Curves ###

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
  ggtheme = theme_bw(base_size = 12),
  title = paste0("Weighted KM (IPTW) – Death (Cutoff = ",
                 dose_cutoff, " mL/kg/hr)")
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
  ggtheme = theme_bw(base_size = 12),
  title = paste0("Weighted KM (IPTW) – Discharge (Cutoff = ", 
                 dose_cutoff, " mL/kg/hr)")
)
plot_km_disch
ggsave(file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                    "_IPTW_KM_Discharge.png")),
       plot_km_disch$plot, width=6, height=4)

cat("Weighted KM plots saved as PNGs.\n")

################################################
### 3. MODEL COMPARISON - PSM FG vs IPTW MSM ###
################################################

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

msm_death_row <- extract_msm_trt(fit_cs_death, "IPTW MSM Cox - Death")
msm_disch_row <- extract_msm_trt(fit_cs_disch, "IPTW MSM Cox - Discharge")

comparison_table <- bind_rows(
  fg_death_row,
  fg_disch_row,
  msm_death_row,
  msm_disch_row
)

comparison_table

write.csv(comparison_table,
          file.path(output_dir, paste0(SITE_NAME,"_", dose_label,
                                       "_ModelComparison_PSMvsIPTW.csv")),
          row.names = FALSE)

cat("Saved model comparison CSV\n")
cat("Code run complete!\n")
