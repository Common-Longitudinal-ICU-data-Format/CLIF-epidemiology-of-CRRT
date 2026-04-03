# Multi-Site Preparation: Issues & Proposed Solutions

*Prepared 2026-04-03. Reference document for preparing the CRRT epidemiology pipeline for deployment across CLIF consortium sites.*

---

## Current State

The pipeline has solid multi-site foundations:
- **Config system** is fully site-agnostic (`config_template.json`)
- **All CSVs include a `site` column** for downstream aggregation
- **`07_combine_site_results.py`** auto-discovers sites in `all_site_data/{site}/final/` and builds a combined HTML dashboard
- **`08_severity_analysis.py`** generates cross-site severity comparisons, forest plots, and meta-regressions
- **Table 1 + visualization pooling** works end-to-end for descriptive analyses

The issues below are gaps that could cause problems when sites beyond UCMC run the pipeline.

---

## Issue 1: No Multi-Site Workflow Documentation

**Problem:** There are no instructions explaining how sites should run the pipeline, what outputs to share, or how the central coordinator combines results. Sites receiving the code would not know the expected workflow.

**Affected files:** `README.md` (or new doc)

**Proposed solution:** Add a "Multi-Site Workflow" section covering:

1. **Per-site setup:**
   - Clone repo, copy `config/config_template.json` to `config/config.json`, fill in site-specific values
   - Install dependencies with `uv sync`
   - Run `bash run_pipeline.sh` (steps 00-03 for descriptive; 04-06b for causal inference)

2. **What to send to coordinator:**
   - The entire `output/final/` directory (three subdirectories: `crrt_epi/`, `msm/`, `psm_iptw/`)
   - `config/outlier_config.json` (for comparability documentation)
   - Confirmation of `has_crrt_settings` value

3. **Coordinator workflow:**
   - Create `all_site_data/` at repo root
   - Place each site's `output/final/` at `all_site_data/{site_name}/final/`
   - Run `uv run python code/07_combine_site_results.py`
   - Run `uv run python code/08_severity_analysis.py`
   - Review `all_site_data/combined_dashboard.html`

4. **Directory layout example:**
   ```
   all_site_data/
     UCMC/
       final/
         crrt_epi/
         msm/
         psm_iptw/
     MIMIC/
       final/
         crrt_epi/
         msm/
         psm_iptw/
   ```

**Effort:** Low (documentation only)

---

## Issue 2: `has_crrt_settings=false` Incomplete Coverage

**Problem:** The `has_crrt_settings` flag correctly controls whether dose computation runs in `00_cohort.py` and whether CRRT columns are loaded in `01_create_wide_df.py`. However, `02_construct_crrt_tableone.py` still renders CRRT settings rows (flow rate, mode, dose) in Table 1 even when the flag is `false`. When pooling sites with mixed settings availability, these rows will have data for some sites and NaN/missing for others.

**Affected files:**
- `code/02_construct_crrt_tableone.py` (~line 743-749, mode-specific rows)
- `code/07_combine_site_results.py` (pooling logic doesn't account for partial CRRT data)

**Proposed solution:**
1. In `02_construct_crrt_tableone.py`: wrap CRRT settings rows (flow rate, mode, dose statistics) in a conditional block checking `config["has_crrt_settings"]`. When `false`, omit these rows entirely from both `table1_crrt.csv` and `table1_crrt_long.csv`.
2. In `07_combine_site_results.py`: when pooling Table 1 across sites, note which sites contributed CRRT settings data. Add a footnote to the combined dashboard indicating partial availability.
3. For dose visualizations: skip overlay generation for dose-related graphs when fewer than 2 sites have CRRT settings.

**Effort:** Low

---

## Issue 3: Causal Analysis Outputs Not Aggregated

**Problem:** `07_combine_site_results.py` currently aggregates only descriptive outputs (Table 1, STROBE counts, visualization CSVs). It does **not** handle any outputs from the causal inference scripts (05, 06, 06b):

- PSM/IPTW treatment effect estimates (`_ModelComparison_PSMvsIPTW.csv`)
- E-value sensitivity (`_evalue_sensitivity.csv`)
- Subgroup analysis results (`_subgroup_analysis_results.csv`)
- MSM time-varying estimates
- Pooled Rubin's rules results (`_IPTW_pooled_results.csv`, `_fg_psm_pooled_results.csv`)

Without this, each site's causal findings exist in isolation with no cross-site comparison.

**Affected files:** `code/07_combine_site_results.py`

**Proposed solution:**
1. **Add PSM/IPTW aggregation section** to script 07:
   - Load each site's `_ModelComparison_PSMvsIPTW.csv` and combine into a forest plot of site-level treatment effects
   - Pool E-values into a summary table
   - Compare subgroup interaction p-values across sites (especially CCI and body weight signals)

2. **Add MSM aggregation section:**
   - Load each site's MSM treatment effect estimates
   - Generate cross-site comparison of positivity diagnostics (% with PS > 0.95)
   - Flag sites where MSM may be unreliable due to positivity violations

3. **Meta-analysis consideration:**
   - For pooling treatment effects across sites, use inverse-variance weighted random-effects meta-analysis
   - Document that per-site MICE imputation means global Rubin's rules are not applied — this is a limitation
   - Consider using the `meta` R package or Python's `statsmodels` for formal meta-analysis

4. **Dashboard integration:**
   - Add a "Causal Inference" tab to the combined dashboard with forest plots of per-site HRs and pooled estimates

**Effort:** Medium-High (new aggregation logic + meta-analysis)

---

## Issue 4: Hardcoded CRRT Mode Names

**Problem:** CRRT mode names are hardcoded in multiple scripts:

- `00_cohort.py`: `DOSE_MODES = {'cvvh', 'cvvhd', 'cvvhdf'}` (line ~1994)
- `02_construct_crrt_tableone.py`: iterates over mode-specific rows
- `07_combine_site_results.py`: regex drops `last_crrt_mode` rows by pattern (line ~126)

If a site uses different mode naming (capitalization, localized terms, additional modes), dose calculations and Table 1 construction break silently.

**Affected files:**
- `code/00_cohort.py`
- `code/02_construct_crrt_tableone.py`
- `code/07_combine_site_results.py`

**Proposed solution:**

Option A (minimal): Add a note in the multi-site documentation that `crrt_mode_category` values must be lowercase and match `{scuf, cvvh, cvvhd, cvvhdf}`. This is already the CLIF 2.1.0 standard, so sites using `clifpy` for validation should comply.

Option B (robust): Move mode lists to config:
```json
{
  "crrt_modes": {
    "valid": ["scuf", "cvvh", "cvvhd", "cvvhdf"],
    "dose_eligible": ["cvvh", "cvvhd", "cvvhdf"]
  }
}
```
Then reference `config["crrt_modes"]["dose_eligible"]` instead of hardcoded sets.

**Recommendation:** Option A is sufficient if all sites use `clifpy` for CLIF table validation (which standardizes category values). Option B is safer but adds configuration complexity.

**Effort:** Low (Option A) / Medium (Option B)

---

## Issue 5: No Data Completeness Reporting

**Problem:** When pooling across sites, some sites may lack certain labs, vitals, or variables. The current pipeline handles missing data gracefully (NaN propagation), but there's no visibility into what's missing. A combined Table 1 might show a variable with N=500 when the total cohort is N=2000 — without explanation that 3 sites didn't have that variable.

**Affected files:** `code/07_combine_site_results.py`, `code/08_severity_analysis.py`

**Proposed solution:**
1. Add a **data completeness matrix** to the combined dashboard:
   - Rows = variables (labs, vitals, CRRT settings, etc.)
   - Columns = sites
   - Cells = % non-missing or "Not collected"

2. Source this from each site's `table1_crrt_long.csv` — variables with 0 N or entirely NaN indicate unavailability.

3. Add a footnote system to the combined Table 1 marking variables with < 100% site coverage.

**Effort:** Medium

---

## Issue 6: Outlier Config Not Portable or Compared

**Problem:** `config/outlier_config.json` defines physiological range thresholds for vitals, labs, and CRRT flow rates. These are applied during cohort construction (`00_cohort.py` via `utils.handle_crrt_outliers()`). If sites use different thresholds — for example, because their CRRT machines have different flow rate ranges — the pooled data is not directly comparable.

**Affected files:**
- `config/outlier_config.json`
- `code/00_cohort.py`

**Proposed solution:**
1. Include `outlier_config.json` in the site output package (alongside `output/final/`)
2. In `07_combine_site_results.py`: load each site's outlier config and generate a comparison table
3. Flag any variables where sites differ by > 20% in their thresholds
4. Add the comparison to the combined dashboard as a "Methods Consistency" tab

**Alternative:** Standardize a single `outlier_config.json` across all sites (simpler, but may not account for legitimate site differences).

**Effort:** Low

---

## Issue 7: Cross-Site Validation Checks

**Problem:** After loading site CSVs, scripts 07 and 08 don't validate consistency:
- Do N totals match across stratifications (e.g., dead + alive = total)?
- Do percentages sum correctly?
- Are STROBE exclusion counts internally consistent?

Silent inconsistencies could propagate into the combined dashboard.

**Affected files:** `code/07_combine_site_results.py`

**Proposed solution:**
1. Add a validation step after site discovery that checks:
   - STROBE: `included = total - sum(exclusions)`
   - Table 1: `N(dead) + N(alive) = N(overall)` for each variable
   - Percentages: categorical variable percentages sum to ~100% within each subgroup
2. Log warnings for any inconsistencies; include a validation report in the dashboard
3. Do not block execution — treat as warnings, not errors (sites may have legitimate rounding differences)

**Effort:** Medium

---

## Priority Order

| Priority | Issue | Effort | Rationale |
|----------|-------|--------|-----------|
| 1 | Workflow documentation (#1) | Low | Prerequisite for any site to run the pipeline |
| 2 | `has_crrt_settings` gaps (#2) | Low | Prevents errors at sites without flow rate data |
| 3 | Hardcoded mode names (#4) | Low | Documentation fix prevents silent failures |
| 4 | Outlier config portability (#6) | Low | Simple to include in output package |
| 5 | Data completeness reporting (#5) | Medium | Important for interpreting pooled results |
| 6 | Cross-site validation (#7) | Medium | Catches errors before they enter the dashboard |
| 7 | Causal analysis aggregation (#3) | Medium-High | Most complex; may require meta-analysis expertise |

---

## Notes

- **CLIF version compatibility:** All scripts assume CLIF 2.1.0 column names. Sites on other versions would need to verify schema compatibility. Consider adding a CLIF version field to `config.json`.
- **Causal inference pooling strategy:** Per-site MICE imputation means we cannot apply global Rubin's rules across sites. The recommended approach is per-site causal analysis with cross-site meta-analysis of point estimates.
- **Site name uniqueness:** `07_combine_site_results.py` uses directory names for site identification. Ensure no two sites use the same `site_name` in their config.
- **UCMC-specific finding:** UCMC has 0 CVVHDF patients — the `crrt_modality` subgroup was excluded. Other sites may have this mode, so the subgroup analysis code handles it generically.
