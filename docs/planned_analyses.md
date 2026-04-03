# Planned Analyses: Target Trial Emulation, Sensitivity Analyses, and Dose-Response

*Created 2026-04-03. Three methodological improvements to strengthen the causal inference component of the CRRT dose study.*

---

## 1. Explicit Target Trial Emulation Framework

### Motivation

The current PSM/IPTW analysis implicitly emulates a trial but never states the target trial specification. Making it explicit accomplishes three things:

- **Clarity**: Readers can evaluate whether the observational design credibly answers the question by comparing it to the hypothetical RCT.
- **Rigor**: Forces us to articulate assumptions about eligibility, treatment assignment, and follow-up that are currently implicit in the code.
- **Convention**: Target trial emulation is increasingly the expected framework for observational causal inference in critical care (Hernan & Robins, *Causal Inference* textbook; Hernan 2016 *Am J Epidemiol*).

### Target Trial Specification

| Component | Specification |
|-----------|--------------|
| **Eligibility** | Adults (age ≥18) with AKI requiring CRRT, without pre-existing ESRD, with documented weight, surviving ≥24h on CRRT, with dialysate flow (excludes SCUF-only). N = 1,633 at UCMC. |
| **Treatment strategies** | Strategy A: Initial CRRT dose ≥30 mL/kg/hr. Strategy B: Initial CRRT dose <30 mL/kg/hr. Dose defined as total effluent flow rate (dialysate + replacement + ultrafiltration) divided by body weight at CRRT initiation. |
| **Treatment assignment** | Patients are assigned at CRRT initiation based on the first recorded dose. In the target trial, this would be randomized; in the emulation, we adjust for confounders via propensity score methods. |
| **Follow-up start (time zero)** | CRRT initiation time (first recorded CRRT therapy event per encounter block). |
| **Outcome** | 30-day all-cause mortality from time zero, with hospital discharge alive as a competing event. |
| **Causal contrast** | Intention-to-treat analogue: effect of initial dose assignment, regardless of subsequent dose changes. |
| **Analysis plan** | (a) PSM with nearest-neighbor matching + Fine-Gray competing risk model. (b) IPTW via Super Learner propensity scores + cause-specific Cox model. Both estimate the effect of initial dose assignment on 30-day mortality. |

### Deviations from an Ideal Trial

| Issue | Impact | How Addressed |
|-------|--------|---------------|
| No randomization | Unmeasured confounding | Propensity score adjustment; E-values for sensitivity |
| Observational dose assignment | Confounding by indication | PS adjusts for SOFA, lactate, oxygenation, vasopressors, IMV |
| Delivered dose ≠ prescribed dose | Misclassification of treatment | Using delivered (recorded) dose; may bias toward null |
| 24h survival requirement | Conditions on post-treatment variable | Necessary to define meaningful exposure; sensitivity analysis planned (see Section 2) |
| Weight measurement in ICU | Dose misclassification (fluid overload) | Non-differential error; biases toward null |
| Time-varying dose not captured | ITT-like estimand, not per-protocol | MSM attempted but positivity violations preclude reliable estimation |

### Implementation Plan

**What to create**: A methods documentation file (`output/final/psm_iptw/target_trial_specification.md`) and a formatted table for the manuscript.

**Code changes**: None required for the analysis itself — the current PSM/IPTW pipeline already implements this design. The target trial specification is a documentation/framing exercise that maps existing code to the framework.

**Optional code addition**: Add the target trial table to the combined dashboard or as a standalone HTML output.

---

## 2. Sensitivity Analysis: Including Patients Excluded by 24h Criterion

### Motivation

Step 04 (`04_build_msm_competing_risk_df.py`) excludes 503 patients:
- Died within 24h of CRRT initiation, OR
- CRRT duration < 24h (discontinued early)

This exclusion is clinically justified — patients who die or stop CRRT within hours haven't received enough therapy for dose to plausibly matter. However, conditioning on a post-treatment variable (surviving 24h on CRRT) can introduce bias:

- **If high-dose CRRT causes hemodynamic instability**: Some patients assigned high dose may die within 24h *because* of the high dose. Excluding them removes dose-related deaths from the high-dose group, biasing toward null.
- **If low-dose patients are sicker**: Sicker patients may be started on lower doses (hemodynamic instability limits flow rates) and die quickly. Excluding them removes the sickest low-dose patients, potentially biasing in either direction.

Showing that the null finding holds even without this exclusion strengthens the overall conclusion.

### Analysis Plan

**Approach**: Re-run the PSM/IPTW analysis on the full 2,136-patient descriptive cohort (or as close as possible — patients without a computable dose at time 0 still need exclusion).

**Specific steps**:

1. **Modify `04_build_msm_competing_risk_df.py`** to optionally skip the 24h exclusion:
   - Add a config flag or command-line argument: `exclude_short_crrt` (default: True)
   - When False, keep all patients but still compute dose, outcomes, and covariates
   - Patients who died within 24h have `time_to_event_30d` values near 0 — these are valid data points
   - SCUF-only patients (no dialysate) still need exclusion since they have no meaningful dose

2. **Create a sensitivity output directory**: `output/final/psm_iptw_sensitivity_24h/`

3. **Run `05_PSM_IPTW_CRRT_dose.R`** on the sensitivity dataset:
   - Same PS model specification
   - Same matching/weighting approach
   - Compare the treatment effect estimate to the primary analysis
   - If estimates are similar → 24h exclusion doesn't drive the null
   - If estimates differ → the exclusion matters and should be discussed

4. **Report**: Side-by-side table of primary vs. sensitivity HRs, CIs, and p-values.

### Expected Challenges

- **Dose availability**: Patients who died within 24h may not have a recorded dose in the first measurement window. Need to check what fraction of the 503 excluded patients have `crrt_dose_ml_kg_hr` available at time 0.
- **Covariate completeness**: Some baseline covariates (SOFA, labs) may be missing for patients who died quickly.
- **Very short follow-up**: Patients dying at hour 2 contribute minimal information to the Cox model. The proportional hazards assumption may be violated.
- **Sample size**: The 503 additional patients include many early deaths, which could shift the mortality rate substantially.

### What Success Looks Like

- If the sensitivity analysis shows HR ≈ 0.98 (similar to primary): Strong evidence that the 24h exclusion doesn't bias the result. Report as: "Results were robust to inclusion of patients with CRRT duration <24h."
- If the sensitivity analysis shows a different HR: Important finding — discuss the direction and possible mechanisms. This would suggest that early CRRT mortality is non-random with respect to dose.

---

## 3. Dose-Response Analysis (Continuous Dose)

### Motivation

Dichotomizing dose at 30 mL/kg/hr discards information and introduces several problems:

- **Arbitrary boundary**: A patient at 29.5 mL/kg/hr is "low dose" and one at 30.5 is "high dose" — clinically identical, statistically different groups.
- **Loss of power**: Continuous analysis uses all available dose variation, not just above/below a cutoff.
- **Assumes no dose-response shape**: The binary analysis assumes the effect is a step function at 30. In reality, there could be:
  - No effect at any dose (flat curve)
  - Linear dose-response (more dose → more/less mortality)
  - U-shaped curve (optimal dose in the middle, harm at extremes)
  - Threshold effect (benefit above some dose, but not necessarily 30)
- **Weight measurement error**: If ICU weight is inaccurate by ±10%, patients near the cutoff are frequently misclassified. A continuous analysis is more robust to this because it doesn't depend on a single threshold.

### Analysis Plan

**Approach A: Descriptive dose-response curve**

1. **Dose decile analysis**: Divide the cohort into deciles by initial CRRT dose. Plot 30-day mortality by decile with 95% CIs. This is model-free and immediately interpretable.
   - Implementation: Simple pandas groupby + binomial CI calculation
   - Output: Figure + CSV in `output/final/psm_iptw/`

2. **Restricted cubic spline (RCS)**: Fit a Cox model with dose as an RCS (4-5 knots at quantile-based positions). Plot the HR as a function of dose, with the median dose as the reference.
   - This reveals nonlinear relationships: U-shapes, thresholds, plateaus
   - Implementation: In R using `rms::cph()` with `rcs()` term, or in Python using `patsy` spline basis + `lifelines`
   - Adjust for the same covariates used in the IPTW model

**Approach B: Continuous treatment effect (generalized propensity score)**

1. **Generalized propensity score (GPS)**: Extend the binary PS framework to continuous treatment.
   - Model dose as a continuous outcome of baseline confounders (e.g., linear regression or random forest)
   - Use the GPS to weight or stratify the outcome analysis
   - Implementation: R package `CBPS` or `WeightIt` with continuous treatment support
   - This is more methodologically complex and may not add much over the RCS approach

**Recommendation**: Start with Approach A (decile plot + RCS). It's simpler, more interpretable, and answers the key question: "Is there any dose level where mortality differs?" If the RCS curve is flat, it's strong evidence of no dose-response at any level — more convincing than a single binary comparison.

### Specific Implementation

**New R script**: `code/05b_dose_response_analysis.R`

```
Section 1: Dose decile analysis
  - Compute deciles of crrt_dose_ml_kg_hr_0
  - Calculate mortality, discharge, and censoring rates per decile
  - Plot with error bars
  - Output: dose_decile_mortality.csv, dose_decile_plot.png

Section 2: Restricted cubic spline
  - Fit Cox PH model: Surv(time, death) ~ rcs(dose, 5) + covariates
  - Plot HR vs. dose with 95% CI band
  - Test for nonlinearity (Wald test on spline terms)
  - Output: dose_response_rcs.png, dose_response_rcs_table.csv

Section 3: Comparison to binary analysis
  - Overlay the binary HR (at cutoff=30) on the continuous curve
  - Report whether 30 appears to be a meaningful threshold
```

**Input**: Same `msm_competing_risk_df.parquet` used by script 05.

**Dependencies**: `rms` package (for RCS), already available `survival` package.

### What Success Looks Like

- **Flat RCS curve**: No dose-response at any level. This would be the strongest possible null result — not just "no effect at the 30 cutoff" but "no effect anywhere in the dose range." Highly publishable and clinically important.
- **Threshold effect**: Mortality drops above some dose (maybe 25, maybe 35). This would suggest the 30 cutoff is approximately right (or wrong) and that dose does matter — just not at the specific cutoff tested.
- **U-shaped curve**: Harm at both extremes. This would be a novel finding suggesting an optimal dose range, consistent with the hypothesis that too-low dose provides insufficient clearance while too-high dose causes hemodynamic stress.

---

## Priority Order

| # | Analysis | Effort | Value |
|---|----------|--------|-------|
| 1 | Target trial specification | Low (documentation only) | High — frames the entire study |
| 2 | Dose-response curve (decile + RCS) | Medium (new R script) | High — could be the most important figure |
| 3 | 24h exclusion sensitivity | Medium (modify step 04 + re-run 05) | Medium — confirms robustness |

---

## References

- Hernan MA, Robins JM. Using big data to emulate a target trial when a randomized trial is not available. *Am J Epidemiol*. 2016;183(8):758-764.
- Hernan MA, Robins JM. *Causal Inference: What If*. Chapman & Hall/CRC, 2020.
- ATN Study Group. Intensity of renal support in critically ill patients with acute kidney injury. *NEJM*. 2008;359(1):7-20.
- RENAL Replacement Therapy Study Investigators. Intensity of continuous renal-replacement therapy in critically ill patients. *NEJM*. 2009;361(17):1627-1638.
- Hirano K, Imbens GW. The propensity score with continuous treatments. In: Gelman A, Meng XL, eds. *Applied Bayesian Modeling and Causal Inference from Incomplete-Data Perspectives*. Wiley; 2004:73-84.
