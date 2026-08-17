"""
Build standalone manuscript artifacts (figures, tables, pooled CSVs) from
the per-site outputs of scripts 04-06b.

Consumes per-site CSVs in `output/multi_site/<site>/final/` and produces:

  output/multi_site/pooled/         7 pooled meta-analysis CSVs
  output/multi_site/figures/        5 manuscript figures (PDF + PNG)
  output/multi_site/tables/         3 manuscript tables (CSV + standalone HTML)
                                      T1 baseline by dose band (full cohort, 02)
                                      T2 unadjusted balance (causal cohort, 05)
                                      T3 pooled point-treatment HRs

Independent of `07_build_dashboard.py` (the dashboard builder). Run order does
not matter between the two — they consume the same per-site CSVs and write to
disjoint output directories. The dashboard does its own per-site meta-analysis
and does NOT read the pooled CSVs that this script emits.

Shared infrastructure (config, IO, color palette, meta-analysis primitives,
combined-table compute, between-site heterogeneity) lives in `report_core.py`,
imported normally below — the dashboard builder imports the same module.

Run standalone:
    uv run python code/08_manuscript_artifacts.py

Or iterate on a single artifact via the `# %% Render Figure N` cells at
the bottom of this file (VS Code Jupyter Interactive Window workflow).
"""
from __future__ import annotations

import html
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import site_validation

from report_core import (
    ROOT,
    SITE_LABELS,
    discover_sites,
    _site_final_dir,
    _load_csv,
    meta_analyze,
    _collect_iptw,
    _collect_psm_fg,
    _se_from_ci,
    _compute_combined_table_df,
    _render_gtsummary_table_html,
    _read_site_causal_n,
    BALANCE_TABLE_PATTERNS,
    compute_epi_heterogeneity,
    EPI_HET_COLS,
    compute_dose_ibw_pooled,
    DOSE_IBW_POOLED_COLS,
    _weighted_median,
    collect_smr,
    collect_smr_calibration,
)


# ── Pooled CSV Exports for Manuscript ─────────────────────────────────────

def _bh_qvalues(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values (BH-FDR). Returns array same shape as input."""
    p = np.asarray(p_values, dtype=float)
    k = len(p)
    if k == 0:
        return p
    order = np.argsort(p)
    ranks = np.empty(k, dtype=float)
    ranks[order] = np.arange(1, k + 1)
    q_raw = p * k / ranks
    # Step-up: enforce monotonicity from largest p downward
    sorted_idx = order
    q_sorted = q_raw[sorted_idx]
    q_sorted_rev = q_sorted[::-1]
    cum_min_rev = np.minimum.accumulate(q_sorted_rev)
    q_sorted = cum_min_rev[::-1]
    q_out = np.empty(k, dtype=float)
    q_out[sorted_idx] = q_sorted
    return np.minimum(q_out, 1.0)


def _read_table1_psm_counts(site_dir: Path) -> dict | None:
    """Parse the unadjusted balance table to recover counts on the PSM-eligible
    cohort. (Manuscript Table 2; the file is named "_table2_unadjusted_balance"
    by 05.R, with the legacy "_Table1_unadjusted" tried as a fallback until all
    sites re-run 05.)

    The balance table is built downstream of the PSM `drop_na` (all model
    covariates + treatment + outcome non-missing), so its column headers and
    outcome-row cells give counts restricted to the same population MatchIt
    sees. This is the canonical denominator for `per_site_summary.csv` — it
    matches `psm_counts_summary.csv` exactly and the manuscript figures derived
    from it.

    Returns a dict with: n_total, n_low_dose, n_high_dose,
    deaths_30d, discharge_n, censored_n. Returns None if the file is missing
    or unparsable.
    """
    d = _site_final_dir(site_dir) / "psm_iptw"
    fp = (next(d.glob("*_table2_unadjusted_balance.csv"), None)
          or next(d.glob("*_Table1_unadjusted.csv"), None))
    if fp is None or not fp.exists():
        return None
    try:
        df = pd.read_csv(fp)
    except Exception:
        return None
    if df.empty or len(df.columns) < 4:
        return None

    # Strip markdown bold + carriage returns from column headers; pull N from each.
    def _n_from_header(col: str) -> int | None:
        m = re.search(r"N\s*=\s*([\d,]+)", str(col))
        return int(m.group(1).replace(",", "")) if m else None

    cols = list(df.columns)
    n_total = _n_from_header(cols[1])
    n_low = _n_from_header(cols[2])
    n_high = _n_from_header(cols[3])
    if n_total is None:
        return None

    # Outcome rows live under a "30-day Outcome" header. Match the leading
    # cell against {Discharged, Died, Censored} (case-insensitive, allowing
    # surrounding whitespace).
    char_col = df.columns[0]
    overall_col = df.columns[1]

    def _n_from_cell(label_pattern: str) -> int | None:
        mask = df[char_col].astype(str).str.strip().str.fullmatch(label_pattern, case=False, na=False)
        if not mask.any():
            return None
        cell = str(df.loc[mask, overall_col].iloc[0])
        m = re.match(r"\s*([\d,]+)\s*\(", cell)
        return int(m.group(1).replace(",", "")) if m else None

    return {
        "n_total": n_total,
        "n_low_dose": n_low,
        "n_high_dose": n_high,
        "deaths_30d": _n_from_cell("Died"),
        "discharge_n": _n_from_cell("Discharged"),
        "censored_n": _n_from_cell("Censored"),
    }


def export_pooled_csvs(sites, labels, out_dir: Path, all_sites=None) -> None:
    """Write pooled meta-analysis numbers as CSVs for downstream manuscript use.
    `sites` is the CAUSAL-ELIGIBLE set: pooled causal estimates (primary
    analyses, subgroup contrasts, dose-response) are computed from it. Per-site
    DIAGNOSTIC outputs are computed from `all_sites` instead, defaulting to
    `sites`. That split is deliberate: per_site_summary.csv carries the arm sizes
    the eligibility rule is applied to, and the cox-diagnostics exports are the
    per-site evidence that adjustment behaved, so restricting either would hide
    the basis for the exclusion (S-Methods 4, site eligibility).

    Outputs (in out_dir):
      per_site_summary.csv                  - cohort N, deaths, mortality, dose stats per site
      pooled_primary_analyses.csv           - PSM FG / IPTW / MSM pooled HRs (FE + RE)
      pooled_primary_per_site.csv           - per-site contributions to each pooled estimate
      pooled_subgroup_levels.csv            - per-(subgroup, level) pooled HR
      pooled_subgroup_interactions.csv      - per-subgroup pooled p-interaction + BH q
      pooled_dose_response_linear.csv       - per-site + FE/RE pooled linear-dose Cox HR
                                              (per mL/kg/hr and per 10 mL/kg/hr)
      pooled_dose_response_nonlinearity.csv - per-site LRT-for-nonlinearity p +
                                              Fisher's combined p
      pooled_cox_diagnostics.csv            - per-site C-statistic + global &
                                              term-level Schoenfeld PH p for the
                                              IPTW cs-Cox and dose-response Cox
                                              models (Supplement table source)
      pooled_cox_diagnostics_per_covariate.csv
                                            - long-format per-(site, model,
                                              covariate) Schoenfeld PH p; used
                                              to identify which specific
                                              covariates have time-varying
                                              effects when global PH fails
    """
    all_sites = list(all_sites) if all_sites is not None else list(sites)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nExporting pooled CSVs to {out_dir}/")

    # ── 1. Per-site summary (cohort + outcomes + dose) ──────────────────
    # Cohort N + outcome counts come from the unadjusted balance table
    # (table2_unadjusted_balance.csv, legacy Table1_unadjusted.csv) which is
    # built post-PSM-eligibility filter so it matches the figures' analytic
    # cohort). Medians/IQRs/dose stats come from sample_characteristics.csv,
    # which is computed pre-eligibility filter on a slightly larger cohort
    # (≤6 patients larger pooled, all with missing initial-dose classification);
    # the median/IQR drift from those few patients is negligible.
    rows = []
    for sd in all_sites:
        sc = _load_csv(sd, "psm_iptw", "sample_characteristics")
        psm = _load_csv(sd, "psm_iptw", "psm_counts_summary")
        if sc is None or sc.empty:
            continue
        s = sc.iloc[0]

        # Canonical cohort counts: prefer Table 1 (PSM-eligible cohort).
        # Fall back to psm_counts_summary for n_total/dose-group counts and
        # to sample_characteristics outcome columns only if Table 1 is missing.
        t1 = _read_table1_psm_counts(sd)
        n_total = None
        n_low = n_high = None
        deaths_30d = discharge_n = censored_n = None
        if t1 is not None:
            n_total = t1["n_total"]
            n_low = t1["n_low_dose"]
            n_high = t1["n_high_dose"]
            deaths_30d = t1["deaths_30d"]
            discharge_n = t1["discharge_n"]
            censored_n = t1["censored_n"]
        if (n_total is None or n_low is None or n_high is None) and psm is not None \
                and not psm.empty and "Category" in psm.columns:
            mask = psm["Category"].astype(str).str.fullmatch("All", case=False, na=False)
            if not mask.any():
                mask = psm["Category"].astype(str).str.contains("All", case=False, na=False) & \
                       ~psm["Category"].astype(str).str.contains("ESS", case=False, na=False)
            if mask.any():
                allrow = psm[mask].iloc[0]
                if n_low is None and pd.notna(allrow.get("Control")):
                    n_low = int(allrow["Control"])
                if n_high is None and pd.notna(allrow.get("Treated")):
                    n_high = int(allrow["Treated"])
                if n_total is None and n_low is not None and n_high is not None:
                    n_total = n_low + n_high
        if n_total is None:
            n_total = int(s.get("total_n", 0)) if pd.notna(s.get("total_n", float("nan"))) else 0

        def _f(col):
            v = s.get(col, float("nan"))
            return float(v) if pd.notna(v) else None

        def _i(col):
            v = s.get(col, float("nan"))
            return int(v) if pd.notna(v) else None

        if deaths_30d is None:
            deaths_30d = _i("outcome_death_n")
        if discharge_n is None:
            discharge_n = _i("outcome_discharge_n")
        if censored_n is None:
            censored_n = _i("outcome_censored_n")

        rows.append({
            "site_id": sd.name,
            "site_label": labels.get(sd.name, sd.name),
            "n_total": n_total,
            "n_high_dose": n_high,
            "n_low_dose": n_low,
            "high_dose_pct": (n_high / n_total * 100) if (n_high is not None and n_total) else None,
            "deaths_30d": deaths_30d,
            "mortality_pct": (deaths_30d / n_total * 100) if (deaths_30d is not None and n_total) else None,
            "discharge_n": discharge_n,
            "discharge_pct": (discharge_n / n_total * 100) if (discharge_n is not None and n_total) else None,
            "censored_n": censored_n,
            "age_median": _f("age_median"),
            "age_q25": _f("age_q25"),
            "age_q75": _f("age_q75"),
            "male_pct": _f("sex_male_pct"),
            "sofa_median": _f("sofa_median"),
            "sofa_q25": _f("sofa_q25"),
            "sofa_q75": _f("sofa_q75"),
            "crrt_dose_median": _f("crrt_dose_median"),
            "crrt_dose_q25": _f("crrt_dose_q25"),
            "crrt_dose_q75": _f("crrt_dose_q75"),
            "lactate_median": _f("lactate_median"),
            "imv_pct": _f("imv_status_pct"),
            "norepi_eq_median": _f("norepi_eq_median"),
            "crrt_duration_days_median": _f("crrt_duration_days_median"),
        })

    if rows:
        df_sites = pd.DataFrame(rows)
        # Append a pooled total row
        pooled = {
            "site_id": "POOLED",
            "site_label": "Pooled (all sites)",
            "n_total": int(df_sites["n_total"].sum()),
            "n_high_dose": int(df_sites["n_high_dose"].fillna(0).sum())
                           if df_sites["n_high_dose"].notna().any() else None,
            "n_low_dose": int(df_sites["n_low_dose"].fillna(0).sum())
                          if df_sites["n_low_dose"].notna().any() else None,
            "deaths_30d": int(df_sites["deaths_30d"].fillna(0).sum())
                          if df_sites["deaths_30d"].notna().any() else None,
            "discharge_n": int(df_sites["discharge_n"].fillna(0).sum())
                           if df_sites["discharge_n"].notna().any() else None,
            "censored_n": int(df_sites["censored_n"].fillna(0).sum())
                          if df_sites["censored_n"].notna().any() else None,
        }
        n = pooled["n_total"]
        if n and pooled["deaths_30d"] is not None:
            pooled["mortality_pct"] = pooled["deaths_30d"] / n * 100
        if n and pooled["discharge_n"] is not None:
            pooled["discharge_pct"] = pooled["discharge_n"] / n * 100
        if n and pooled["n_high_dose"] is not None:
            pooled["high_dose_pct"] = pooled["n_high_dose"] / n * 100
        df_sites = pd.concat([df_sites, pd.DataFrame([pooled])], ignore_index=True)
        df_sites.to_csv(out_dir / "per_site_summary.csv", index=False)
        print(f"  wrote per_site_summary.csv ({len(rows)} sites + pooled row)")

    # ── 2. Pooled primary analyses (PSM FG / IPTW) ────────────────
    # ASCII analysis names so downstream tools handle them cleanly.
    primary_analyses = [
        ("Point-Treatment IPTW - Death",        _collect_iptw,    dict(outcome_kw="Death")),
        ("Point-Treatment IPTW - Discharge",    _collect_iptw,    dict(outcome_kw="Discharge")),
        ("Point-Treatment PSM FG - Death",      _collect_psm_fg,  dict(outcome_kw="Death")),
        ("Point-Treatment PSM FG - Discharge",  _collect_psm_fg,  dict(outcome_kw="Discharge")),
    ]

    pooled_rows, per_site_rows = [], []
    for title, collector, kwargs in primary_analyses:
        s_lbls, s_ids, log_h, se_h, h, lo, hi = collector(sites, labels, **kwargs)
        for sl, sid, hr, l, u in zip(s_lbls, s_ids, h, lo, hi):
            per_site_rows.append({
                "analysis": title, "site_id": sid, "site_label": sl,
                "hr": hr, "hr_lower": l, "hr_upper": u,
            })
        row = {
            "analysis": title,
            "n_sites": len(log_h),
            "site_ids": ";".join(s_ids),
        }
        if len(log_h) >= 2:
            r = meta_analyze(log_h, se_h)
            row.update({
                "fe_hr": r["fe_hr"], "fe_lo": r["fe_lo"], "fe_hi": r["fe_hi"], "fe_p": r["fe_p"],
                "re_hr": r["re_hr"], "re_lo": r["re_lo"], "re_hi": r["re_hi"], "re_p": r["re_p"],
                "I2_pct": r["I2"], "tau2": r["tau2"], "Q": r["Q"], "Q_p": r["Q_p"],
            })
        pooled_rows.append(row)

    pd.DataFrame(pooled_rows).to_csv(out_dir / "pooled_primary_analyses.csv", index=False)
    print(f"  wrote pooled_primary_analyses.csv ({len(pooled_rows)} analyses)")
    if per_site_rows:
        pd.DataFrame(per_site_rows).to_csv(out_dir / "pooled_primary_per_site.csv", index=False)
        print(f"  wrote pooled_primary_per_site.csv ({len(per_site_rows)} rows)")

    # ── 3. Subgroup-level pooled HRs + interactions ─────────────────────
    site_subgroup = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "subgroup_analysis_results")
        if df is None:
            continue
        site_subgroup.append((sd.name, labels.get(sd.name, sd.name), df))

    sub_rows, int_rows = [], []
    if site_subgroup:
        subgroup_order = []
        for _, _, df in site_subgroup:
            for s in df["subgroup"].unique():
                if s not in subgroup_order:
                    subgroup_order.append(s)

        for subgroup in subgroup_order:
            all_levels = []
            for _, _, df in site_subgroup:
                sub = df[df["subgroup"] == subgroup]
                for lv in sub["level"].unique():
                    if lv not in all_levels:
                        all_levels.append(lv)

            for level in all_levels:
                log_hrs, se_hrs, sites_in = [], [], []
                n_sum = events_sum = 0
                for sid, _slbl, df in site_subgroup:
                    sub = df[(df["subgroup"] == subgroup) & (df["level"] == level)]
                    if sub.empty:
                        continue
                    r = sub.iloc[0]
                    if pd.notna(r.get("n")):
                        n_sum += int(r["n"])
                    if pd.notna(r.get("events")):
                        events_sum += int(r["events"])
                    hr, lo, hi = r.get("HR"), r.get("HR_lower"), r.get("HR_upper")
                    if pd.isna(hr) or hr <= 0 or pd.isna(lo) or lo <= 0 or pd.isna(hi) or hi <= 0:
                        continue
                    se = _se_from_ci(hr, lo, hi)
                    if np.isnan(se) or se <= 0:
                        continue
                    log_hrs.append(np.log(hr))
                    se_hrs.append(se)
                    sites_in.append(sid)
                row = {
                    "subgroup": subgroup,
                    "level": level,
                    "n_sites_with_hr": len(log_hrs),
                    "n_total": n_sum,
                    "events_total": events_sum,
                    "site_ids": ";".join(sites_in),
                }
                if len(log_hrs) >= 2:
                    res = meta_analyze(np.array(log_hrs), np.array(se_hrs))
                    row.update({
                        "fe_hr": res["fe_hr"], "fe_lo": res["fe_lo"], "fe_hi": res["fe_hi"], "fe_p": res["fe_p"],
                        "re_hr": res["re_hr"], "re_lo": res["re_lo"], "re_hi": res["re_hi"], "re_p": res["re_p"],
                        "I2_pct": res["I2"], "tau2": res["tau2"], "Q": res["Q"], "Q_p": res["Q_p"],
                    })
                sub_rows.append(row)

            if len(all_levels) >= 2 and subgroup.lower() != "overall":
                lv_a, lv_b = all_levels[0], all_levels[1]
                diff_log_hrs, diff_ses, sites_in = [], [], []
                for sid, _slbl, df in site_subgroup:
                    sub = df[df["subgroup"] == subgroup]
                    ra = sub[sub["level"] == lv_a]
                    rb = sub[sub["level"] == lv_b]
                    if ra.empty or rb.empty:
                        continue
                    ra, rb = ra.iloc[0], rb.iloc[0]
                    if pd.isna(ra["HR"]) or pd.isna(rb["HR"]) or ra["HR"] <= 0 or rb["HR"] <= 0:
                        continue
                    se_a = _se_from_ci(ra["HR"], ra["HR_lower"], ra["HR_upper"])
                    se_b = _se_from_ci(rb["HR"], rb["HR_lower"], rb["HR_upper"])
                    if np.isnan(se_a) or np.isnan(se_b):
                        continue
                    diff_log_hrs.append(np.log(rb["HR"]) - np.log(ra["HR"]))
                    diff_ses.append(np.sqrt(se_a**2 + se_b**2))
                    sites_in.append(sid)
                irow = {
                    "subgroup": subgroup,
                    "level_a": lv_a,
                    "level_b": lv_b,
                    "n_sites": len(diff_log_hrs),
                    "site_ids": ";".join(sites_in),
                }
                if len(diff_log_hrs) >= 2:
                    dm = meta_analyze(np.array(diff_log_hrs), np.array(diff_ses))
                    irow.update({
                        "fe_log_diff": float(np.log(dm["fe_hr"])),
                        "re_log_diff": float(np.log(dm["re_hr"])),
                        "fe_p_int": dm["fe_p"],
                        "re_p_int": dm["re_p"],
                        "I2_pct": dm["I2"],
                        "Q_p": dm["Q_p"],
                    })
                int_rows.append(irow)

    if sub_rows:
        pd.DataFrame(sub_rows).to_csv(out_dir / "pooled_subgroup_levels.csv", index=False)
        print(f"  wrote pooled_subgroup_levels.csv ({len(sub_rows)} rows)")

    if int_rows:
        df_int = pd.DataFrame(int_rows)
        if "re_p_int" in df_int.columns and df_int["re_p_int"].notna().any():
            valid_mask = df_int["re_p_int"].notna()
            df_int.loc[valid_mask, "BH_q_re_p_int"] = _bh_qvalues(
                df_int.loc[valid_mask, "re_p_int"].to_numpy()
            )
        df_int.to_csv(out_dir / "pooled_subgroup_interactions.csv", index=False)
        print(f"  wrote pooled_subgroup_interactions.csv ({len(int_rows)} rows)")

    # ── 5. Dose-response meta-analysis (linear-dose HR + LRT for nonlinearity) ──
    # Per-site linear-dose Cox HR per mL/kg/hr lives in *_dose_response_linear.csv;
    # LRT for nonlinearity (spline vs linear, df=3) lives in
    # *_dose_response_rcs_anova.csv (Spline dose row).
    #
    # BACKWARD COMPATIBILITY (2026-08-04): before the 05b cleanup this estimate was
    # emitted only as the comparison arm inside *_sensitivity_24h_exclusion.csv
    # ("Primary" row), whose 24h-exclusion sensitivity analysis has since been
    # removed as uninterpretable (it varied cohort, covariate set, imputation and
    # winsorization simultaneously). Sites that have not yet re-run still ship only
    # the old file, so we prefer the new one and fall back to the old. The numeric
    # value is identical: both are the coefficient of the same linear-dose Cox
    # model (verified locally 2026-08-04, hr agreed to 9 decimal places).
    # Drop the fallback once every site has re-run.
    def _dose_linear_row(sd):
        new = _load_csv(sd, "psm_iptw", "dose_response_linear")
        if new is not None and not new.empty:
            return new.iloc[0], "dose_response_linear"
        old = _load_csv(sd, "psm_iptw", "sensitivity_24h_exclusion")
        if old is not None and not old.empty and "analysis" in old.columns:
            mask = old["analysis"].astype(str).str.startswith("Primary")
            if mask.any():
                return old[mask].iloc[0], "sensitivity_24h_exclusion (legacy)"
        return None, None

    dr_linear_rows = []
    dr_lrt_rows = []
    _legacy_dr_sites = []
    for sd in sites:
        site_label = labels.get(sd.name, sd.name)

        p_row, _src = _dose_linear_row(sd)
        if p_row is not None:
            if _src and _src.endswith("(legacy)"):
                _legacy_dr_sites.append(sd.name)
            if pd.notna(p_row.get("hr")) and pd.notna(p_row.get("se_log_hr")):
                dr_linear_rows.append({
                    "site": sd.name,
                    "site_label": site_label,
                    "n": int(p_row.get("n", 0)) if pd.notna(p_row.get("n")) else None,
                    "n_deaths": int(p_row.get("n_deaths", 0)) if pd.notna(p_row.get("n_deaths")) else None,
                    "hr": float(p_row["hr"]),
                    "hr_lo": float(p_row["hr_lower"]) if pd.notna(p_row.get("hr_lower")) else None,
                    "hr_hi": float(p_row["hr_upper"]) if pd.notna(p_row.get("hr_upper")) else None,
                    "log_hr": float(np.log(p_row["hr"])),
                    "se_log_hr": float(p_row["se_log_hr"]),
                    "p": float(p_row.get("p_value")) if pd.notna(p_row.get("p_value")) else None,
                    "source": _src,
                    "row_type": "per_site",
                })

        anova = _load_csv(sd, "psm_iptw", "dose_response_rcs_anova")
        if anova is not None and not anova.empty and "model" in anova.columns:
            spl_mask = anova["model"].astype(str).str.startswith("Spline")
            if spl_mask.any():
                sp = anova[spl_mask].iloc[0]
                if pd.notna(sp.get("p_value")):
                    dr_lrt_rows.append({
                        "site": sd.name,
                        "site_label": site_label,
                        "chi_sq": float(sp.get("chi_sq")) if pd.notna(sp.get("chi_sq")) else None,
                        "df_diff": 3,
                        "p_lrt_nonlinearity": float(sp["p_value"]),
                        "row_type": "per_site",
                    })

    if _legacy_dr_sites:
        print(f"  NOTE: linear dose-response read from the LEGACY "
              f"sensitivity_24h_exclusion.csv for {len(_legacy_dr_sites)} site(s): "
              f"{', '.join(sorted(_legacy_dr_sites))}. These sites predate the 05b "
              f"cleanup; the value is the same model's coefficient. Remove the "
              f"fallback in _dose_linear_row() once all sites have re-run.")

    if len(dr_linear_rows) >= 2:
        df_dr = pd.DataFrame(dr_linear_rows)
        m = meta_analyze(df_dr["log_hr"].to_numpy(), df_dr["se_log_hr"].to_numpy())
        fe_log = float(np.log(m["fe_hr"]))
        re_log = float(np.log(m["re_hr"]))
        z_crit = 1.959964
        # Recover SE on log-HR scale from the FE/RE intervals (meta_analyze does
        # not expose SE directly, only HR and bounds).
        fe_se_log = float((np.log(m["fe_hi"]) - np.log(m["fe_lo"])) / (2 * z_crit))
        re_se_log = float((np.log(m["re_hi"]) - np.log(m["re_lo"])) / (2 * z_crit))
        pooled_linear = [
            {"site": "POOLED_FE", "site_label": "Pooled FE per mL/kg/hr",
             "hr": m["fe_hr"], "hr_lo": m["fe_lo"], "hr_hi": m["fe_hi"],
             "log_hr": fe_log, "se_log_hr": fe_se_log, "p": m["fe_p"],
             "Q": m["Q"], "Q_p": m["Q_p"], "I2_pct": m["I2"], "tau2": m["tau2"],
             "row_type": "pooled_fe_per_ml_kg_hr"},
            {"site": "POOLED_RE", "site_label": "Pooled RE per mL/kg/hr",
             "hr": m["re_hr"], "hr_lo": m["re_lo"], "hr_hi": m["re_hi"],
             "log_hr": re_log, "se_log_hr": re_se_log, "p": m["re_p"],
             "Q": m["Q"], "Q_p": m["Q_p"], "I2_pct": m["I2"], "tau2": m["tau2"],
             "row_type": "pooled_re_per_ml_kg_hr"},
            {"site": "POOLED_FE_10", "site_label": "Pooled FE per 10 mL/kg/hr",
             "hr": float(np.exp(fe_log * 10)),
             "hr_lo": float(np.exp((fe_log - z_crit * fe_se_log) * 10)),
             "hr_hi": float(np.exp((fe_log + z_crit * fe_se_log) * 10)),
             "log_hr": fe_log * 10, "se_log_hr": fe_se_log * 10, "p": m["fe_p"],
             "Q": m["Q"], "Q_p": m["Q_p"], "I2_pct": m["I2"], "tau2": m["tau2"],
             "row_type": "pooled_fe_per_10_ml_kg_hr"},
            {"site": "POOLED_RE_10", "site_label": "Pooled RE per 10 mL/kg/hr",
             "hr": float(np.exp(re_log * 10)),
             "hr_lo": float(np.exp((re_log - z_crit * re_se_log) * 10)),
             "hr_hi": float(np.exp((re_log + z_crit * re_se_log) * 10)),
             "log_hr": re_log * 10, "se_log_hr": re_se_log * 10, "p": m["re_p"],
             "Q": m["Q"], "Q_p": m["Q_p"], "I2_pct": m["I2"], "tau2": m["tau2"],
             "row_type": "pooled_re_per_10_ml_kg_hr"},
        ]
        out = pd.concat([df_dr, pd.DataFrame(pooled_linear)], ignore_index=True)
        out.to_csv(out_dir / "pooled_dose_response_linear.csv", index=False)
        print(f"  wrote pooled_dose_response_linear.csv ({len(out)} rows)")

    if len(dr_lrt_rows) >= 2:
        df_lrt = pd.DataFrame(dr_lrt_rows)
        pvals = df_lrt["p_lrt_nonlinearity"].to_numpy()
        # Fisher's combined p:  -2 * sum(ln p_i)  ~  chi2(2k)
        fisher_chi2 = float(-2.0 * np.log(pvals).sum())
        fisher_df = 2 * len(pvals)
        fisher_p = float(stats.chi2.sf(fisher_chi2, fisher_df))
        # NOTE (removed 2026-08-04): a Stouffer combination used to be written here
        # alongside Fisher's. It was invalid as implemented and is deliberately gone.
        #
        # Stouffer's method requires SIGNED z from one-sided p-values. The
        # nonlinearity LRT is an omnibus chi-square test with no direction, so the
        # code converted a two-sided p to |z| = Phi^-1(1 - p/2) and summed those.
        # Under the null each |z| has mean sqrt(2/pi) ~ 0.798 rather than 0, so the
        # statistic sum(|z|)/sqrt(k) has null mean ~0.798*sqrt(k) -- about 2.53 at
        # k=10 -- and drifts upward with the number of sites no matter what the data
        # say. Scoring it against N(0,1) therefore reported p = 0.059, reading as
        # borderline evidence of nonlinearity, when the observed value sat at the
        # 14th percentile of its OWN null (correct p ~ 0.86), i.e. LESS extreme than
        # chance. Simulated 2026-08-04 to confirm.
        #
        # Fisher's method is the appropriate combination for independent omnibus
        # tests and is retained. Do not reinstate Stouffer here without a signed
        # per-site statistic to combine, which an LRT does not provide.
        pooled_lrt = [
            {"site": "POOLED_FISHER", "site_label": "Fisher's combined p (LRT for nonlinearity)",
             "chi_sq": fisher_chi2, "df_diff": fisher_df,
             "p_lrt_nonlinearity": fisher_p, "row_type": "fishers_combined"},
        ]
        out = pd.concat([df_lrt, pd.DataFrame(pooled_lrt)], ignore_index=True)
        out.to_csv(out_dir / "pooled_dose_response_nonlinearity.csv", index=False)
        print(f"  wrote pooled_dose_response_nonlinearity.csv ({len(out)} rows)")

    # ── 8. Cox PH and discrimination diagnostics ────────────────────────
    # Per-site C-statistic, global Schoenfeld PH p, and treatment/dose PH p
    # for: (a) the IPTW cs-Cox models (death, discharge) from script 05,
    # and (b) the linear + spline dose-response Cox models from script 05b.
    # Provides a Supplement-table view of model adequacy across sites for
    # the manuscript's primary point-treatment and continuous-dose models.
    cox_diag_rows = []
    for sd in all_sites:
        lbl = labels.get(sd.name, sd.name)
        iptw = _load_csv(sd, "psm_iptw", "iptw_cox_diagnostics")
        if iptw is not None and not iptw.empty:
            for _, row in iptw.iterrows():
                cox_diag_rows.append({
                    "site": sd.name,
                    "site_label": lbl,
                    "script": "05",
                    "model": row.get("model", ""),
                    "n": row.get("n"),
                    "n_events": row.get("n_events"),
                    "c_statistic": row.get("c_statistic"),
                    "c_se": row.get("c_se"),
                    "global_ph_p": row.get("global_ph_p"),
                    "primary_term_ph_p": row.get("trt_ph_p"),
                    "primary_term": "treatment (crrt_high)",
                })
        dose = _load_csv(sd, "psm_iptw", "dose_response_cox_diagnostics")
        if dose is not None and not dose.empty:
            for _, row in dose.iterrows():
                cox_diag_rows.append({
                    "site": sd.name,
                    "site_label": lbl,
                    "script": "05b",
                    "model": row.get("model", ""),
                    "n": row.get("n"),
                    "n_events": row.get("n_events"),
                    "c_statistic": row.get("c_statistic"),
                    "c_se": row.get("c_se"),
                    "global_ph_p": row.get("global_ph_p"),
                    "primary_term_ph_p": row.get("dose_ph_p"),
                    "primary_term": "dose",
                })
    if cox_diag_rows:
        cox_diag_df = pd.DataFrame(cox_diag_rows)
        cox_diag_df.to_csv(out_dir / "pooled_cox_diagnostics.csv", index=False)
        print(
            f"  wrote pooled_cox_diagnostics.csv ({len(cox_diag_df)} rows "
            f"across {cox_diag_df['site'].nunique()} sites, "
            f"{cox_diag_df['model'].nunique()} models)"
        )

    # ── 8b. Per-covariate Schoenfeld PH (long format) ───────────────────
    # Granular drill-down for sites/models where global PH fails. One row
    # per (site, model, covariate); GLOBAL row already lives in the
    # summary CSV above.
    per_cov_rows = []
    for sd in all_sites:
        lbl = labels.get(sd.name, sd.name)
        iptw = _load_csv(sd, "psm_iptw",
                         "iptw_cox_diagnostics_per_covariate")
        if iptw is not None and not iptw.empty:
            for _, row in iptw.iterrows():
                per_cov_rows.append({
                    "site": sd.name,
                    "site_label": lbl,
                    "script": "05",
                    "model": row.get("model", ""),
                    "term": row.get("term", ""),
                    "chi_sq": row.get("chi_sq"),
                    "df": row.get("df"),
                    "p_value": row.get("p_value"),
                })
        dose = _load_csv(sd, "psm_iptw",
                         "dose_response_cox_diagnostics_per_covariate")
        if dose is not None and not dose.empty:
            for _, row in dose.iterrows():
                per_cov_rows.append({
                    "site": sd.name,
                    "site_label": lbl,
                    "script": "05b",
                    "model": row.get("model", ""),
                    "term": row.get("term", ""),
                    "chi_sq": row.get("chi_sq"),
                    "df": row.get("df"),
                    "p_value": row.get("p_value"),
                })
    if per_cov_rows:
        per_cov_df = pd.DataFrame(per_cov_rows)
        per_cov_df.to_csv(
            out_dir / "pooled_cox_diagnostics_per_covariate.csv", index=False
        )
        print(
            f"  wrote pooled_cox_diagnostics_per_covariate.csv "
            f"({len(per_cov_df)} rows across "
            f"{per_cov_df['site'].nunique()} sites, "
            f"{per_cov_df['term'].nunique()} unique terms)"
        )


# ── Manuscript Figures ─────────────────────────────────────────────────────

def build_cohort_flow_figure(sites, out_dir: Path, causal_excluded=None) -> dict:
    """Build the pooled cohort exclusion flow diagram.

    Placement in the manuscript is editorial and is tracked in docs/display_items.csv;
    do not encode a figure number here.

    `causal_excluded` is the output of split_causal_eligible(). When supplied and
    non-empty, a final SITE-LEVEL step is appended showing the sites and patients
    dropped from the pooled causal analysis by the minimum-arm rule. That step is
    drawn only on the causal branch: those sites remain in Table 1, the SMR and
    every descriptive analysis, so the box is labelled as a site-level exclusion
    from the pooled causal estimate rather than as patients leaving the study.

    Replicates the layout of build_combined_consort() but renders to
    standalone PDF + PNG for manuscript inclusion. Same data sources:
    per-site `*_strobe_counts.csv` (steps 1b-6) and per-site
    `*_psm_counts_summary.csv` (causal cohort N).

    Sites missing psm_counts_summary.csv are excluded from the pooled
    diagram (they would otherwise make the final-box count wrong).
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    out_dir.mkdir(parents=True, exist_ok=True)

    site_data = []
    for sd in sites:
        df = _load_csv(sd, "crrt_epi", "strobe_counts")
        if df is None or df.empty or "counter" not in df.columns or "value" not in df.columns:
            continue
        counts = dict(zip(df["counter"].astype(str), df["value"]))
        n_causal = _read_site_causal_n(sd)
        site_data.append((sd.name, counts, n_causal))

    full = [(sid, c, nc) for sid, c, nc in site_data if nc is not None]
    if len(full) < 2:
        print(f"  [cohort flow] need psm_counts_summary.csv from >=2 sites; skipping")
        return {}

    def pool(key: str) -> int:
        return int(sum(c.get(key, 0) or 0 for _, c, _ in full))

    pooled_counts = {
        "1b_after_stitching": pool("1b_after_stitching"),
        "2_crrt_blocks": pool("2_crrt_blocks"),
        "3_encounter_blocks_without_esrd": pool("3_encounter_blocks_without_esrd"),
        "4_encounter_blocks_with_weight": pool("4_encounter_blocks_with_weight"),
        "5_encounter_blocks_with_crrt_settings": pool("5_encounter_blocks_with_crrt_settings"),
        "6_encounter_blocks_with_required_labs": pool("6_encounter_blocks_with_required_labs"),
        "causal": int(sum(nc for _, _, nc in full)),
    }
    start_n = pooled_counts["1b_after_stitching"]
    n_sites = len(full)

    step_defs = [
        (pooled_counts["2_crrt_blocks"],
         "CRRT encounter blocks", "Excluded: No CRRT"),
        (pooled_counts["3_encounter_blocks_without_esrd"],
         "After ESRD exclusion", "Excluded: ESRD diagnosis"),
        (pooled_counts["4_encounter_blocks_with_weight"],
         "With documented weight and sex", "Excluded: Missing weight or sex"),
        (pooled_counts["5_encounter_blocks_with_crrt_settings"],
         "With CRRT settings", "Excluded: Missing CRRT settings"),
        (pooled_counts["6_encounter_blocks_with_required_labs"],
         "Descriptive analytic cohort\n(Table 1 + descriptive analyses)", "Excluded: Missing required labs"),
        (pooled_counts["causal"],
         "Causal analysis cohort\n(complete-case, dichotomized dose)",
         "Excluded: Died or off CRRT ≤24h /\nSCUF-only"),
    ]

    # Site-level step: sites whose higher-dose arm is too small to support the
    # contrast are dropped from POOLED CAUSAL estimates only. Patient counts here
    # are whole sites leaving the pooled estimate, not individually excluded
    # patients, and the labels say so.
    n_excl_sites = len(causal_excluded or [])
    if n_excl_sites:
        n_excl_pts = sum(int(e.get("n_total") or 0) for e in causal_excluded)
        step_defs.append((
            max(pooled_counts["causal"] - n_excl_pts, 0),
            f"Pooled causal analysis\n({n_sites - n_excl_sites} eligible sites)",
            f"Excluded at site level: {n_excl_sites} site(s)\n"
            f"with high-dose arm <{MIN_HIGH_DOSE_ARM}\n"
            f"(retained for Table 1, SMR and\nall descriptive analyses)",
        ))

    steps = []
    parent = start_n
    for remaining, remaining_lbl, excl_lbl in step_defs:
        excluded = max(parent - remaining, 0)
        if excl_lbl.endswith("CRRT settings") and excluded == 0:
            continue
        steps.append({
            "remaining_n": remaining,
            "remaining_label": remaining_lbl,
            "excluded_n": excluded,
            "excluded_label": excl_lbl,
        })
        parent = remaining

    n_boxes = len(steps) + 1
    fig_h = max(8.5, 1.5 * n_boxes)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box_h = 0.08
    box_w = 0.40
    x_main_start = 0.05
    x_excl_start = 0.55
    excl_arrow_gap = 0.015
    v_spacing = 0.86 / n_boxes

    def draw_box(x, y, w, h, text, fontsize=12, weight="normal", facecolor="white"):
        rect = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.01",
            linewidth=2, edgecolor="black", facecolor=facecolor,
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text,
                ha="center", va="center", fontsize=fontsize,
                fontweight=weight, wrap=True)
        return x + w / 2, y

    ax.text(0.5, 0.98,
            f"Pooled Cohort Selection ({n_sites} CLIF Consortium Sites)",
            ha="center", va="center", fontsize=15, fontweight="bold")
    arrow_style = dict(arrowstyle="->", lw=2, color="black")

    top_y = 0.90 - box_h
    draw_box(
        x_main_start, top_y, box_w, box_h,
        f"All adult hospitalizations\n(pooled across {n_sites} sites)\nn = {start_n:,}",
    )

    for i, step in enumerate(steps):
        y_parent = top_y - (i * v_spacing)
        current_y = top_y - ((i + 1) * v_spacing)
        box_center_y_top = y_parent + box_h / 2
        box_center_y_bottom = current_y + box_h / 2
        arrow_vertical_center = (box_center_y_top + box_center_y_bottom) / 2

        is_final = (i == len(steps) - 1)
        # Highlight BOTH result cohorts: the descriptive analytic cohort
        # (second-to-last) and the causal cohort (last).
        is_cohort = (i >= len(steps) - 2)
        remain_text = (
            f"{step['remaining_label']}\nn = {step['remaining_n']:,}"
            if is_cohort else
            f"Remaining hospitalizations\n{step['remaining_label']}\n"
            f"n = {step['remaining_n']:,}"
        )
        facecolor = "#f1f5f9" if is_cohort else "white"
        weight = "bold" if is_cohort else "normal"
        remain_center_x, _ = draw_box(
            x_main_start, current_y, box_w, box_h, remain_text,
            facecolor=facecolor, weight=weight,
        )

        if step["excluded_n"] > 0:
            excl_text = f"{step['excluded_label']}\nn = {step['excluded_n']:,}"
            draw_box(
                x_excl_start,
                arrow_vertical_center - box_h / 2,
                box_w, box_h, excl_text,
            )

        ax.annotate("", xy=(remain_center_x, current_y + box_h),
                    xytext=(remain_center_x, y_parent),
                    arrowprops=arrow_style)
        if step["excluded_n"] > 0:
            ax.annotate(
                "", xy=(x_excl_start - excl_arrow_gap, arrow_vertical_center),
                xytext=(remain_center_x, arrow_vertical_center),
                arrowprops=arrow_style, annotation_clip=False,
            )

    pdf_path = out_dir / "pooled_cohort_flow_diagram.pdf"
    png_path = out_dir / "pooled_cohort_flow_diagram.png"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"  wrote {pdf_path.name}")
    print(f"  wrote {png_path.name}")

    return {"pdf": pdf_path, "png": png_path}


def build_iptw_forest_figure(sites, out_dir: Path, pooled_dir: Path) -> dict:
    """Build the site-specific and pooled IPTW csHR forest plot
    for 30-day mortality.

    Reads `pooled_primary_per_site.csv` and `pooled_primary_analyses.csv`
    from `pooled_dir`. Sites are anonymized as "Site 1, ..., Site N" in the
    figure (per the global anonymization rule); ordering follows the site
    order in the per-site CSV (alphabetical by site_id from
    export_pooled_csvs).
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

    out_dir.mkdir(parents=True, exist_ok=True)

    per_site_path = pooled_dir / "pooled_primary_per_site.csv"
    pooled_path = pooled_dir / "pooled_primary_analyses.csv"
    if not per_site_path.exists() or not pooled_path.exists():
        print(f"  [iptw forest] missing pooled CSVs in {pooled_dir}; skipping")
        return {}

    per_site_df = pd.read_csv(per_site_path)
    pooled_df = pd.read_csv(pooled_path)

    site_rows = per_site_df[
        per_site_df["analysis"] == "Point-Treatment IPTW - Death"
    ].reset_index(drop=True)
    pooled_match = pooled_df[
        pooled_df["analysis"] == "Point-Treatment IPTW - Death"
    ]
    if site_rows.empty or pooled_match.empty:
        print(f"  [iptw forest] no IPTW Death analysis row; skipping")
        return {}
    pooled_row = pooled_match.iloc[0]

    n_sites = len(site_rows)
    site_anon_labels = [f"Site {i+1}" for i in range(n_sites)]

    fig_h = max(4.0, 0.55 * (n_sites + 2) + 1.5)
    fig, (ax, ax_txt) = plt.subplots(
        1, 2, figsize=(10, fig_h),
        gridspec_kw={"width_ratios": [3, 1.3]},
    )
    ax_txt.axis("off")

    for i, (_, r) in enumerate(site_rows.iterrows()):
        y = n_sites - i  # Site 1 at top
        hr = r["hr"]
        lo = r["hr_lower"]
        hi = r["hr_upper"]
        ax.errorbar(
            hr, y, xerr=[[hr - lo], [hi - hr]],
            fmt="s", color="#374151", markersize=8, capsize=4,
            linewidth=1.8, zorder=2,
        )
        ax_txt.text(0.0, y, site_anon_labels[i],
                    va="center", ha="left", fontsize=11, fontweight="bold")
        ax_txt.text(0.45, y, f"{hr:.2f} ({lo:.2f}–{hi:.2f})",
                    va="center", ha="left", fontsize=11)

    # Pooled random-effects diamond only exists with >= 2 sites (meta_analyze is
    # written into the pooled row only then, 08 ~L331). With a single causal-
    # eligible site there is nothing to pool, so skip the diamond and label it
    # rather than KeyError'ing the whole manuscript-figure build.
    y_re = 0.0  # one full row below Site N for breathing room
    if pd.notna(pooled_row.get("re_hr")):
        re_hr = pooled_row["re_hr"]
        re_lo = pooled_row["re_lo"]
        re_hi = pooled_row["re_hi"]
        i2 = pooled_row["I2_pct"]
        diamond_x = [re_lo, re_hr, re_hi, re_hr]
        diamond_y = [y_re, y_re - 0.30, y_re, y_re + 0.30]
        ax.fill(diamond_x, diamond_y, color="#1e417c", alpha=0.90, zorder=3)
        ax.plot([re_lo, re_hi], [y_re, y_re], color="#1e417c", linewidth=2, zorder=3)
        ax_txt.text(0.0, y_re,
                    f"Pooled (RE)\nI² = {i2:.0f}%",
                    va="center", ha="left", fontsize=11, fontweight="bold", color="#1e417c")
        ax_txt.text(0.45, y_re, f"{re_hr:.2f} ({re_lo:.2f}–{re_hi:.2f})",
                    va="center", ha="left", fontsize=11, fontweight="bold", color="#1e417c")
    else:
        ax_txt.text(0.0, y_re, "Pooled\n(n/a: 1 site)",
                    va="center", ha="left", fontsize=11, fontweight="bold", color="#6b7280")

    ax.axvline(x=1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.65, 1.4)
    ax.xaxis.set_major_locator(FixedLocator([0.7, 0.85, 1.0, 1.2, 1.4]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.7", "0.85", "1.0", "1.2", "1.4"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("IPTW Cause-Specific Hazard Ratio (95% CI) — log scale", fontsize=11)
    ax.set_yticks([])
    ax.set_ylim(-0.7, n_sites + 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", labelsize=10)

    ax_txt.set_xlim(0, 1)
    ax_txt.set_ylim(-0.7, n_sites + 1)
    ax_txt.text(0.0, n_sites + 0.6, "Site",
                va="center", ha="left", fontsize=11, fontweight="bold")
    ax_txt.text(0.45, n_sites + 0.6, "HR (95% CI)",
                va="center", ha="left", fontsize=11, fontweight="bold")

    fig.suptitle(
        "IPTW Cause-Specific Hazard Ratio for 30-Day Mortality (High vs. Low CRRT Dose)",
        fontsize=12.5, fontweight="bold", y=0.97,
    )
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.93])

    pdf_path = out_dir / "pooled_iptw_forest.pdf"
    png_path = out_dir / "pooled_iptw_forest.png"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"  wrote {pdf_path.name}")
    print(f"  wrote {png_path.name}")

    return {"pdf": pdf_path, "png": png_path}


def build_psm_cif_figure(sites, out_dir: Path) -> dict:
    """Build the pooled PSM-matched cumulative incidence figure
    of 30-day mortality by initial CRRT dose group across CLIF consortium sites.

    Reads per-site `*_PSM_CIF_data.csv` (long format with columns: group,
    outcome, time, cif, cif_lower, cif_upper). Filters to outcome == 2
    (Death). Each site's CIF is interpolated onto a common 30-day time grid
    via step interpolation; site-specific curves are pooled to a mean curve
    with a 95% CI ribbon (t-distribution across sites) for each dose group.
    Both groups are plotted on a single panel with overlaid ribbons.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from scipy.stats import t as _t_dist

    out_dir.mkdir(parents=True, exist_ok=True)

    TIME_MAX = 30.0
    n_grid = 300
    tgrid = np.linspace(0, TIME_MAX, n_grid)

    # Group definitions: (raw label in CIF csv, display label, color)
    # ASN brand palette: deep blue + orange.
    groups = [
        ("<30",  "Low dose (<30 mL/kg/hr)",  "#1e417c"),
        (">=30", "High dose (≥30 mL/kg/hr)", "#fb801b"),
    ]

    per_group_curves: dict[str, list[np.ndarray]] = {gid: [] for gid, _, _ in groups}
    n_sites_with_data = 0
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "PSM_CIF_data")
        if df is None or df.empty:
            continue
        if not all(c in df.columns for c in ("group", "outcome", "time", "cif")):
            continue
        sub = df[df["outcome"] == 2]
        if sub.empty:
            continue
        site_added = False
        for gid, _, _ in groups:
            gdf = sub[sub["group"].astype(str) == gid].sort_values("time")
            if gdf.empty:
                continue
            interp = np.interp(
                tgrid,
                gdf["time"].to_numpy(),
                gdf["cif"].to_numpy(),
                left=0.0,
                right=np.nan,
            )
            per_group_curves[gid].append(interp)
            site_added = True
        if site_added:
            n_sites_with_data += 1

    if not all(per_group_curves[gid] for gid, _, _ in groups):
        print("  [psm cif] insufficient per-site PSM CIF data; skipping")
        return {}

    def _pooled_ci(arr: np.ndarray, conf: float = 0.95):
        mean = np.nanmean(arr, axis=0)
        sd = np.nanstd(arr, axis=0, ddof=1)
        n = (~np.isnan(arr)).sum(axis=0)
        valid = n >= 2
        df_ = np.where(valid, n - 1, 1)
        t_crit = np.where(valid, _t_dist.ppf((1 + conf) / 2, df=df_), np.nan)
        se = np.where(valid, sd / np.sqrt(np.maximum(n, 1)), np.nan)
        half = t_crit * se
        return mean, mean - half, mean + half

    fig, ax = plt.subplots(figsize=(8, 5.2))
    legend_handles: list = []

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="Mean of empty slice")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="Degrees of freedom <= 0 for slice")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="invalid value encountered")
        for gid, label, color in groups:
            arr = np.vstack(per_group_curves[gid])
            mean, lo, hi = _pooled_ci(arr)
            valid_mean = ~np.isnan(mean)
            valid_band = valid_mean & ~np.isnan(lo) & ~np.isnan(hi)
            ax.fill_between(
                tgrid[valid_band],
                np.clip(lo[valid_band], 0, 1),
                np.clip(hi[valid_band], 0, 1),
                color=color, alpha=0.18, zorder=1,
            )
            ax.plot(tgrid[valid_mean], mean[valid_mean],
                    color=color, linewidth=2.5, zorder=3)
            legend_handles.append(Line2D([0], [0], color=color, linewidth=2.5,
                                         label=label))
    legend_handles.append(Patch(facecolor="gray", alpha=0.30,
                                label="95% CI across sites"))

    ax.set_xlim(0, TIME_MAX)
    ax.set_ylim(0, None)
    ax.set_xlabel("Time from CRRT Initiation (Days)", fontsize=11)
    ax.set_ylabel("Cumulative Incidence of Death", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25)
    ax.legend(handles=legend_handles, loc="upper left",
              frameon=False, fontsize=10)

    fig.suptitle(
        "PSM-Matched Cumulative Incidence of 30-Day Mortality by Initial CRRT Dose",
        fontsize=12.5, fontweight="bold", y=0.97,
    )
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.94])

    pdf_path = out_dir / "pooled_psm_cif_death.pdf"
    png_path = out_dir / "pooled_psm_cif_death.png"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"  wrote {pdf_path.name}")
    print(f"  wrote {png_path.name}")

    return {"pdf": pdf_path, "png": png_path}


def build_dose_response_figure(sites, out_dir: Path) -> dict:
    """Build the 3-panel continuous dose-response figure.

    Layout: 1 x 3 horizontal.
      Panel (a) - Unadjusted dose-decile mortality
      Panel (b) - Covariate-adjusted natural cubic spline Cox HR
      Panel (c) - Doubly-robust GPS-CBPS HR

    Each panel shows per-site curves overlaid (anonymized site colors,
    Site 1 ... Site N) and a pooled mean curve with a 95% CI ribbon
    across sites. Curves are interpolated to a common dose grid for the
    pooled summary; per-site curves are drawn at their native dose values.
    Panels (b)/(c) smoothed in log-HR space (sigma=8 grid pts ≈ 3 mL/kg/hr)
    so site entry/exit doesn't show as slope discontinuities.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from scipy.stats import t as _t_dist
    from scipy.ndimage import gaussian_filter1d as _gaussian_smooth

    out_dir.mkdir(parents=True, exist_ok=True)

    # DOSE_FILTER bounds the per-site data we accept (full physiologic range);
    # DOSE_XLIM is the visible plot range, which trims the choppy ribbon edges
    # at dose extremes (<4 sites contributing). Matplotlib clips drawn curves
    # to the axis range automatically.
    DOSE_FILTER = (0, 80)
    DOSE_XLIM = (10, 60)
    grid = np.linspace(DOSE_FILTER[0], DOSE_FILTER[1], 200)
    # Suppress np.nanmean/std warnings on grid columns where every site is NaN
    # or only one site contributes (legitimate at dose-grid extremes; valid
    # masks filter them out before plotting).
    _ignore_allnan = warnings.catch_warnings()
    _ignore_allnan.__enter__()
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN slice encountered")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="Degrees of freedom <= 0 for slice")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered")

    # Show ribbon only at grid points with >= N_MIN_RIBBON sites; with t-CI,
    # n=2 gives t_crit=12.7 and n=3 gives 4.30, which produce noisy spikes
    # at dose-grid extremes. n>=4 (t_crit<=3.18) is a clean compromise.
    N_MIN_RIBBON = 4

    def _pooled_ci_linear(arr, conf: float = 0.95, n_min: int = N_MIN_RIBBON):
        """Per-grid-point mean and `conf` CI across sites on linear scale (t-distribution).

        Returns (mean, lo, hi) arrays the same length as the grid; ribbon
        entries are NaN where fewer than `n_min` sites contribute. The mean
        itself is preserved wherever any site contributes.
        """
        mean = np.nanmean(arr, axis=0)
        sd = np.nanstd(arr, axis=0, ddof=1)
        n = (~np.isnan(arr)).sum(axis=0)
        valid = n >= n_min
        df = np.where(valid, n - 1, 1)
        t_crit = np.where(valid, _t_dist.ppf((1 + conf) / 2, df=df), np.nan)
        se = np.where(valid, sd / np.sqrt(np.maximum(n, 1)), np.nan)
        half = t_crit * se
        return mean, mean - half, mean + half

    def _pooled_ci_log(arr, conf: float = 0.95, n_min: int = N_MIN_RIBBON,
                       smooth_sigma_grid: float | None = None):
        """Same as `_pooled_ci_linear` but in log space; back-transformed to HR scale.

        If `smooth_sigma_grid` is provided, applies a 1D Gaussian smoother
        (in grid units; sigma is in *grid points*, not mL/kg/hr) to the
        log-mean and log-CI-bounds within their contiguous non-NaN segment
        before exponentiation. Smoothing in log space keeps HR>0.
        """
        with np.errstate(invalid="ignore", divide="ignore"):
            log_arr = np.log(arr)
        log_mean, log_lo, log_hi = _pooled_ci_linear(log_arr, conf=conf, n_min=n_min)
        if smooth_sigma_grid is not None:
            log_mean = _smooth_nan_segment(log_mean, smooth_sigma_grid)
            log_lo = _smooth_nan_segment(log_lo, smooth_sigma_grid)
            log_hi = _smooth_nan_segment(log_hi, smooth_sigma_grid)
        return np.exp(log_mean), np.exp(log_lo), np.exp(log_hi)

    def _smooth_nan_segment(y: np.ndarray, sigma: float) -> np.ndarray:
        """Apply 1D Gaussian smoothing to the contiguous non-NaN portion of y.

        NaN entries are preserved; values inside the contiguous valid range
        are replaced with their Gaussian-smoothed counterparts (mode="nearest"
        edge handling). Caller passes `sigma` in *grid points*.
        """
        valid = ~np.isnan(y)
        if valid.sum() < 3:
            return y
        out = y.copy()
        # Smooth only the contiguous valid run; if there are gaps, smooth each run
        # independently so a NaN gap doesn't bleed across the smoother.
        idx = np.where(valid)[0]
        runs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
        for run in runs:
            if len(run) >= 3:
                out[run] = _gaussian_smooth(y[run], sigma=sigma, mode="nearest")
        return out

    site_labels = [f"Site {i+1}" for i in range(len(sites))]
    cmap = plt.colormaps["tab10"]
    site_colors = [cmap(i % 10) for i in range(len(sites))]

    def _interp(df, dose_col, value_col):
        d = df[dose_col].to_numpy()
        v = df[value_col].to_numpy()
        order = np.argsort(d)
        return np.interp(grid, d[order], v[order], left=np.nan, right=np.nan)

    decile_curves = []
    spline_curves = []
    gps_curves = []
    for i, sd in enumerate(sites):
        lbl = site_labels[i]
        color = site_colors[i]
        d = _load_csv(sd, "psm_iptw", "dose_decile_mortality")
        if d is not None and not d.empty and "dose_median" in d.columns and "mortality_rate" in d.columns:
            d = d[(d["dose_median"] >= DOSE_FILTER[0]) & (d["dose_median"] <= DOSE_FILTER[1])]
            if not d.empty:
                decile_curves.append((lbl, color, d))
        d = _load_csv(sd, "psm_iptw", "dose_response_rcs")
        if d is not None and not d.empty:
            dose_col = next((c for c in ("crrt_dose_median_3h", "crrt_dose_ml_kg_hr_0")
                             if c in d.columns), None)
            if dose_col:
                d = d[(d[dose_col] >= DOSE_FILTER[0]) & (d[dose_col] <= DOSE_FILTER[1])]
                if not d.empty:
                    spline_curves.append((lbl, color, d, dose_col))
        d = _load_csv(sd, "psm_iptw", "gps_dose_response")
        if d is not None and not d.empty:
            dose_col = next((c for c in ("crrt_dose_median_3h", "crrt_dose_ml_kg_hr_0")
                             if c in d.columns), None)
            if dose_col:
                d = d[(d[dose_col] >= DOSE_FILTER[0]) & (d[dose_col] <= DOSE_FILTER[1])]
                if not d.empty:
                    gps_curves.append((lbl, color, d, dose_col))

    if not (decile_curves and spline_curves and gps_curves):
        print(f"  [dose response] insufficient per-site dose-response data; skipping")
        return {}

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # ── Panel (a): Decile mortality ────────────────────────────────────────
    ax = axes[0]
    decile_grid_values = []
    for lbl, color, d in decile_curves:
        ax.plot(d["dose_median"], d["mortality_rate"], "o-",
                color=color, linewidth=1.2, markersize=4, alpha=0.35, zorder=2)
        decile_grid_values.append(_interp(d, "dose_median", "mortality_rate"))
    if decile_grid_values:
        arr = np.vstack(decile_grid_values)
        mean, lo, hi = _pooled_ci_linear(arr)
        valid = ~np.isnan(mean) & ~np.isnan(lo) & ~np.isnan(hi)
        ax.fill_between(grid[valid], lo[valid], hi[valid],
                        color="black", alpha=0.15, zorder=1)
        plot_mean = ~np.isnan(mean)
        ax.plot(grid[plot_mean], mean[plot_mean],
                color="black", linewidth=2.5, zorder=3)
    ax.set_xlim(*DOSE_XLIM)
    ax.set_ylim(0, 1)
    ax.set_xlabel("CRRT Dose (mL/kg/hr)", fontsize=11)
    ax.set_ylabel("30-Day Mortality Rate", fontsize=11)
    ax.set_title("(a) Unadjusted Dose-Decile Mortality", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Panel (b): Spline Cox ───────────────────────────────────────────────
    # Smooth the pooled mean + CI bounds in log-HR space (sigma=8 grid points
    # ≈ 3.2 mL/kg/hr) so site entry/exit doesn't show as slope discontinuities.
    # Mask the pooled mean to the same support as the ribbon so the line
    # doesn't extend past the band visibly.
    ax = axes[1]
    spline_grid_values = []
    for lbl, color, d, dose_col in spline_curves:
        ax.plot(d[dose_col], d["hr"], color=color, linewidth=1.2, alpha=0.35, zorder=2)
        spline_grid_values.append(_interp(d, dose_col, "hr"))
    if spline_grid_values:
        arr = np.vstack(spline_grid_values)
        mean, lo, hi = _pooled_ci_log(arr, smooth_sigma_grid=8)
        valid = ~np.isnan(mean) & ~np.isnan(lo) & ~np.isnan(hi)
        ax.fill_between(grid[valid], lo[valid], hi[valid],
                        color="black", alpha=0.15, zorder=1)
        ax.plot(grid[valid], mean[valid],
                color="black", linewidth=2.5, zorder=3)
    ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, zorder=1)
    ax.set_xlim(*DOSE_XLIM)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 2.0)
    ax.yaxis.set_major_locator(FixedLocator([0.5, 0.7, 1.0, 1.5, 2.0]))
    ax.yaxis.set_major_formatter(FixedFormatter(["0.5", "0.7", "1.0", "1.5", "2.0"]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("CRRT Dose (mL/kg/hr)", fontsize=11)
    ax.set_ylabel("Hazard Ratio (log scale)", fontsize=11)
    ax.set_title("(b) Covariate-Adjusted Natural Cubic Spline", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Panel (c): GPS ──────────────────────────────────────────────────────
    ax = axes[2]
    gps_grid_values = []
    for lbl, color, d, dose_col in gps_curves:
        ax.plot(d[dose_col], d["hr"], color=color, linewidth=1.2, alpha=0.35, zorder=2)
        gps_grid_values.append(_interp(d, dose_col, "hr"))
    if gps_grid_values:
        arr = np.vstack(gps_grid_values)
        mean, lo, hi = _pooled_ci_log(arr, smooth_sigma_grid=8)
        valid = ~np.isnan(mean) & ~np.isnan(lo) & ~np.isnan(hi)
        ax.fill_between(grid[valid], lo[valid], hi[valid],
                        color="black", alpha=0.15, zorder=1)
        ax.plot(grid[valid], mean[valid],
                color="black", linewidth=2.5, zorder=3)
    ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, zorder=1)
    ax.set_xlim(*DOSE_XLIM)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 2.0)
    ax.yaxis.set_major_locator(FixedLocator([0.5, 0.7, 1.0, 1.5, 2.0]))
    ax.yaxis.set_major_formatter(FixedFormatter(["0.5", "0.7", "1.0", "1.5", "2.0"]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("CRRT Dose (mL/kg/hr)", fontsize=11)
    ax.set_ylabel("Hazard Ratio (log scale)", fontsize=11)
    ax.set_title("(c) Doubly-Robust GPS-CBPS", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        Line2D([0], [0], color=site_colors[i], linewidth=2, label=site_labels[i])
        for i in range(len(sites))
    ]
    legend_handles.append(Line2D([0], [0], color="black", linewidth=2.5, label="Pooled mean"))
    legend_handles.append(Patch(facecolor="black", alpha=0.15, label="95% CI across sites"))
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=min(8, len(legend_handles)), fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, -0.01))

    fig.suptitle(
        "Continuous Dose-Response: CRRT Dose and 30-Day Mortality Across Six CLIF Consortium Sites",
        fontsize=12.5, fontweight="bold", y=0.9,
    )
    plt.tight_layout(rect=[0.02, 0.07, 0.98, 0.93])

    pdf_path = out_dir / "pooled_dose_response.pdf"
    png_path = out_dir / "pooled_dose_response.png"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    _ignore_allnan.__exit__(None, None, None)

    print(f"  wrote {pdf_path.name}")
    print(f"  wrote {png_path.name}")

    return {"pdf": pdf_path, "png": png_path}


def build_subgroup_forest_figure(sites, out_dir: Path, pooled_dir: Path) -> dict:
    """Build the multi-site subgroup forest plot in a 4x2 landscape grid.

    8 testable subgroups across 2 rows x 4 columns. Each panel shows:
      - 6 site-specific HRs as small jittered dots (no per-site CIs)
      - Pooled random-effects HRs as diamonds with 95% CI whiskers (one per stratum)
      - Vertical reference line at HR = 1
      - Panel header: subgroup name + pooled p-int + BH q-value

    Reads pooled CSVs from `pooled_dir` (assumes export_pooled_csvs has run).
    Writes PDF + PNG to `out_dir`. Sites are not labeled on the figure (per
    the manuscript anonymization rule); per-site detail with labels and CIs
    belongs in a Supplement figure.

    Returns: dict with PDF/PNG output paths, or empty dict if inputs missing.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

    out_dir.mkdir(parents=True, exist_ok=True)

    levels_path = pooled_dir / "pooled_subgroup_levels.csv"
    interactions_path = pooled_dir / "pooled_subgroup_interactions.csv"
    if not levels_path.exists() or not interactions_path.exists():
        print(f"  [subgroup forest] missing pooled CSVs in {pooled_dir}; skipping")
        return {}
    levels = pd.read_csv(levels_path)
    interactions = pd.read_csv(interactions_path)

    pi = interactions[
        (interactions["n_sites"] >= 2) & (interactions["subgroup"] != "Overall")
    ].copy()
    pi = pi[pi["re_p_int"].notna()]
    poolable = pi["subgroup"].tolist()
    if not poolable:
        print(f"  [subgroup forest] no poolable subgroups; skipping")
        return {}

    per_site_dfs = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "subgroup_analysis_results")
        if df is None or df.empty:
            continue
        df = df.copy()
        df["site_id"] = sd.name
        per_site_dfs.append(df)
    if not per_site_dfs:
        print(f"  [subgroup forest] no per-site subgroup data found; skipping")
        return {}
    per_site = pd.concat(per_site_dfs, ignore_index=True)
    per_site["HR_num"] = pd.to_numeric(per_site["HR"], errors="coerce")

    n_panels = min(8, len(poolable))
    n_cols = 4
    n_rows = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 5.5), sharex=True)
    axes = axes.flatten()

    # ASN brand palette: deep blue + orange.
    color_a = "#1e417c"  # blue: lower stratum
    color_b = "#fb801b"  # orange: upper stratum
    rng = np.random.default_rng(42)

    for i in range(n_panels):
        ax = axes[i]
        subgroup = poolable[i]

        # Preserve the level order from the source CSV (Low/High, Male/Female,
        # No IMV/On IMV are not alphabetical — sorting would invert them).
        sg_levels = levels[levels["subgroup"] == subgroup].copy().reset_index(drop=True)
        if len(sg_levels) < 2:
            ax.set_visible(False)
            continue

        sg_int = pi[pi["subgroup"] == subgroup].iloc[0]
        p_int = sg_int["re_p_int"]
        q = sg_int.get("BH_q_re_p_int", float("nan"))

        for j, lvl_row in sg_levels.iterrows():
            level = lvl_row["level"]
            color = color_a if j == 0 else color_b
            # Each stratum gets a vertical "lane": dots in the upper half,
            # pooled diamond in the lower half, so they don't visually overlap.
            y_center = (1 - j) * 1.0     # j=0 -> y=1 (top stratum), j=1 -> y=0
            y_dots = y_center + 0.22
            y_pooled = y_center - 0.22

            site_hrs = per_site[
                (per_site["subgroup"] == subgroup) &
                (per_site["level"] == level)
            ]["HR_num"].dropna().tolist()
            for site_hr in site_hrs:
                jitter = rng.uniform(-0.07, 0.07)
                ax.scatter(
                    site_hr, y_dots + jitter,
                    s=30, color=color, alpha=0.78,
                    edgecolors="white", linewidths=0.6, zorder=3,
                )

            pooled_hr = lvl_row["re_hr"]
            pooled_lo = lvl_row["re_lo"]
            pooled_hi = lvl_row["re_hi"]
            ax.errorbar(
                pooled_hr, y_pooled,
                xerr=[[max(0.0, pooled_hr - pooled_lo)],
                      [max(0.0, pooled_hi - pooled_hr)]],
                fmt="D", color=color, markersize=8, linewidth=1.4,
                capsize=3, capthick=1.0, zorder=2,
                markeredgecolor="black", markeredgewidth=0.4,
            )

            # Stratum label sits in the inter-panel gutter (x just outside the
            # left axis edge in axes-fraction coords) so it doesn't collide
            # with the dots/diamond when the HR range is tight.
            ax.text(
                -0.04, y_center, str(level),
                va="center", ha="right", fontsize=8.5, color=color,
                transform=ax.get_yaxis_transform(),
                fontweight="bold",
                clip_on=False,
            )

        ax.axvline(1.0, color="gray", linestyle="--",
                   linewidth=0.8, alpha=0.7, zorder=1)
        ax.set_xscale("log")
        ax.set_xlim(0.6, 1.5)
        # Override matplotlib's default log-scale formatter, which produces
        # scientific notation like "4 × 10^-1" that breaks readability.
        ax.xaxis.set_major_locator(FixedLocator([0.6, 0.8, 1.0, 1.25, 1.5]))
        ax.xaxis.set_major_formatter(FixedFormatter(["0.6", "0.8", "1.0", "1.25", "1.5"]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.tick_params(axis="x", labelsize=9)
        ax.set_yticks([])
        ax.set_ylim(-0.5, 1.5)
        ax.set_title(
            f"{subgroup}\np-int = {p_int:.2f}   BH q = {q:.2f}",
            fontsize=10, pad=1,
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", alpha=0.2, zorder=0)

    for i in range(n_panels, n_rows * n_cols):
        axes[i].set_visible(False)

    fig.suptitle(
        "Subgroup-Stratified Pooled IPTW csHR for 30-Day Mortality (High vs. Low CRRT Dose)",
        fontsize=12, fontweight="bold", y=0.90,
    )
    fig.text(
        0.5, 0.055,
        "x-axis (log scale): cause-specific hazard ratio for high-dose vs. low-dose CRRT, by stratum.",
        ha="center", fontsize=8.5,
    )
    fig.text(
        0.5, 0.020,
        "Each stratum shows site-specific HRs as dots, and pooled random-effects HRs with 95% CI as a diamond and whiskers."
        "Dashed line: HR=1.  No q-value reached q < 0.05.",
        ha="center", fontsize=8, style="italic",
    )

    plt.tight_layout(rect=[0.02, 0.09, 0.98, 0.94])

    pdf_path = out_dir / "pooled_subgroup_forest.pdf"
    png_path = out_dir / "pooled_subgroup_forest.png"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"  wrote {pdf_path.name}")
    print(f"  wrote {png_path.name}")

    return {"pdf": pdf_path, "png": png_path}


# ── Manuscript Tables ──────────────────────────────────────────────────────

_MANUSCRIPT_TABLE_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    color: #0f172a;
    max-width: 1100px;
    margin: 2rem auto;
    padding: 0 1.5rem;
  }}
  h1 {{ font-size: 1.15rem; margin: 0 0 0.25rem; }}
  p.caption {{ color: #475569; font-size: 0.9rem; margin: 0 0 1rem; }}
  table.results-table {{
    border-collapse: collapse;
    width: auto;
    margin: 0;
    font-size: 0.92rem;
  }}
  table.results-table th,
  table.results-table td {{
    padding: 0.45rem 0.85rem;
    border-bottom: 1px solid #e2e8f0;
    text-align: left;
    vertical-align: top;
  }}
  table.results-table thead th {{
    background: #f1f5f9;
    border-bottom: 2px solid #cbd5e1;
    font-weight: 600;
    white-space: pre-line;
  }}
  table.results-table tbody tr:hover {{ background: #f8fafc; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="caption">{caption}</p>
{table}
</body>
</html>
"""


def _build_combined_manuscript_table(
    sites, labels, out_dir: Path, *, subdir, pattern, out_stem, title,
    caption_lead, log_tag, overall_only: bool = False,
    include_site_range: bool = False,
) -> dict:
    """Pool a per-site gtsummary-shaped table across sites and render it as a
    standalone manuscript table (CSV + HTML).

    Shared by manuscript Table 1 (full-cohort baseline by dose band, from 02's
    `crrt_epi/table1_crrt`) and Table 2 (unadjusted balance on the causal
    cohort, from 05's `psm_iptw/table2_unadjusted_balance`). Both pool through
    `_compute_combined_table_df()` (weighted medians for continuous; pooled
    counts with the full pooled column N as denominator for categorical) and
    write the raw combined DataFrame (headers preserve embedded Ns) plus an
    HTML using the manuscript-table chrome. `caption_lead` is prepended to the
    standard per-site-contribution caption tail.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    result_df, site_names, site_overall_ns = _compute_combined_table_df(
        sites, labels, subdir, pattern, include_site_range=include_site_range
    )
    if result_df is None:
        print(f"  [{log_tag}] no per-site files found ({subdir}/{pattern}); skipping")
        return {}

    if overall_only:
        keep = [result_df.columns[0]] + [
            c for c in result_df.columns
            if "combined" in str(c).lower() or "between-site range" in str(c).lower()
        ]
        result_df = result_df[keep]

    csv_path = out_dir / f"{out_stem}.csv"
    html_path = out_dir / f"{out_stem}.html"

    result_df.to_csv(csv_path, index=False)

    n_sites = len(site_names)
    total_n = sum(site_overall_ns)
    site_info = ", ".join(
        f"{name} (N={n:,})" for name, n in zip(site_names, site_overall_ns)
    )
    has_smd = any("smd" in str(c).lower() for c in result_df.columns)
    smd_note = (
        " SMD, standardized mean difference; |SMD| > 0.10 conventionally "
        "indicates a meaningful between-group difference."
        if has_smd else ""
    )
    has_range = any(
        "between-site range" in str(c).lower() for c in result_df.columns
    )
    range_note = (
        " Between-site range is the span from the lowest to the highest "
        "per-site summary value across contributing sites (each site's median "
        "for continuous variables, each site's percentage for categorical "
        "variables); it is not the minimum and maximum of individual patient "
        "values."
        if has_range else ""
    )
    caption = (
        f"{caption_lead} Combined from {n_sites} CLIF consortium "
        f"sites (N = {total_n:,}). Per-site contributions: {site_info}. "
        "Continuous variables reported as weighted median (IQR); "
        f"categorical variables as n (% of pooled cohort).{range_note}{smd_note}"
    )
    html_doc = _MANUSCRIPT_TABLE_HTML_TEMPLATE.format(
        title=title,
        caption=caption,
        table=_render_gtsummary_table_html(result_df),
    )
    html_path.write_text(html_doc, encoding="utf-8")

    print(f"  wrote {csv_path.name}")
    print(f"  wrote {html_path.name}")
    return {"csv": csv_path, "html": html_path}


def build_baseline_table(sites, labels, out_dir: Path) -> dict:
    """Manuscript Table 1: pooled full-cohort baseline characteristics of the
    overall CRRT analytic cohort (single Combined column, not stratified by
    dose band), sourced from 02's per-site `crrt_epi/table1_crrt.csv`. Dose-band
    variation is presented as a figure, not as Table 1 strata; the per-site
    dose-band CSVs remain available as a supplement."""
    return _build_combined_manuscript_table(
        sites, labels, out_dir,
        subdir="crrt_epi", pattern="table1_crrt",
        out_stem="pooled_baseline_overall",
        title="Table 1 — Baseline Characteristics of the CRRT Cohort",
        caption_lead=(
            "Baseline characteristics of the full CRRT analytic cohort at "
            "CRRT initiation."
        ),
        log_tag="Table 1",
        overall_only=True,
        include_site_range=True,
    )


def build_unadjusted_balance_table(sites, labels, out_dir: Path) -> dict:
    """Manuscript Table 2: pooled unadjusted balance between dose groups on the
    causal (PSM-eligible) cohort, sourced from 05's per-site
    `psm_iptw/table2_unadjusted_balance.csv` (legacy `Table1_unadjusted.csv`
    fallback). Shows the confounding-by-indication pattern that motivates the
    causal analyses (high- vs low-intensity, < 30 vs >= 30 mL/kg/hr)."""
    return _build_combined_manuscript_table(
        sites, labels, out_dir,
        subdir="psm_iptw", pattern=BALANCE_TABLE_PATTERNS,
        out_stem="pooled_unadjusted_balance",
        title=("Table 2 — Unadjusted Balance Between Dose Groups "
               "(Causal Cohort)"),
        caption_lead=(
            "Unadjusted baseline balance between high-intensity (>= 30 mL/kg/hr) "
            "and low-intensity (< 30 mL/kg/hr) initial CRRT dose on the "
            "causal (PSM-eligible) cohort."
        ),
        log_tag="Table 2",
    )


def build_causal_hrs_table(sites, labels, out_dir: Path, pooled_dir: Path) -> dict:
    """Build the pooled point-treatment HR table (PSM Fine-Gray
    SHR + IPTW csHR) for 30-day mortality and discharge alive across CLIF
    consortium sites.

    Reads `pooled_primary_analyses.csv` (produced by `export_pooled_csvs()`)
    and renders four pooled rows: IPTW csHR (Death, Discharge) and PSM FG
    SHR (Death, Discharge). Each row carries the random-effects pooled
    estimate (95% CI), I², and Cochran's Q p-value. Per-site detail is in
    Figure 2 (Death) and the underlying source CSVs.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    pooled_path = pooled_dir / "pooled_primary_analyses.csv"
    if not pooled_path.exists():
        print(f"  [causal hrs] missing {pooled_path}; skipping")
        return {}

    df = pd.read_csv(pooled_path)
    # Cohort N: derive from sites using the same path that drives all the
    # other manuscript denominators (per-site Table1 headers Control+Treated).
    n_total = 0
    for sd in sites:
        t1 = _read_table1_psm_counts(sd)
        if t1 is not None and t1.get("n_total"):
            n_total += t1["n_total"]
    n_sites = sum(
        1 for sd in sites
        if _read_table1_psm_counts(sd) is not None
    )

    # Order matches the Results-section narrative: Death first (primary),
    # Discharge second (competing event); IPTW (csHR) before PSM (SHR) to
    # mirror Figure 2's primary estimand.
    rows_spec = [
        ("Point-Treatment IPTW - Death",
         "IPTW cause-specific HR — 30-day mortality"),
        ("Point-Treatment IPTW - Discharge",
         "IPTW cause-specific HR — discharge alive (competing event)"),
        ("Point-Treatment PSM FG - Death",
         "PSM Fine-Gray subdistribution HR — 30-day mortality"),
        ("Point-Treatment PSM FG - Discharge",
         "PSM Fine-Gray subdistribution HR — discharge alive (competing event)"),
    ]

    out_rows = []
    for analysis_key, display_label in rows_spec:
        match = df[df["analysis"] == analysis_key]
        if match.empty:
            print(f"  [causal hrs] missing pooled row '{analysis_key}'; skipping it")
            continue
        r = match.iloc[0]
        re_hr = r["re_hr"]
        re_lo = r["re_lo"]
        re_hi = r["re_hi"]
        re_p = r["re_p"]
        i2 = r["I2_pct"]
        q_p = r["Q_p"]
        out_rows.append({
            "Estimand": display_label,
            "Random-effects pooled HR (95% CI)":
                f"{re_hr:.2f} ({re_lo:.2f}–{re_hi:.2f})",
            "p-value": f"{re_p:.2f}" if re_p >= 0.01 else f"{re_p:.3f}",
            "I² (%)": f"{i2:.0f}",
            "Cochran's Q p": f"{q_p:.2f}" if q_p >= 0.01 else f"{q_p:.3f}",
        })

    if not out_rows:
        print("  [causal hrs] no rows matched in pooled_primary_analyses.csv; skipping")
        return {}

    out_df = pd.DataFrame(out_rows)

    csv_path = out_dir / "pooled_causal_hrs.csv"
    html_path = out_dir / "pooled_causal_hrs.html"
    out_df.to_csv(csv_path, index=False)

    caption = (
        f"Pooled point-treatment causal estimates from {n_sites} CLIF consortium "
        f"sites (N = {n_total:,}). Random-effects (DerSimonian-Laird) pooled hazard "
        "ratios for high-intensity (≥ 30 mL/kg/hr) versus low-intensity "
        "(< 30 mL/kg/hr) initial CRRT dose. PSM = propensity score matching "
        "(Fine-Gray subdistribution model); IPTW = inverse probability of "
        "treatment weighting (covariate-adjusted cause-specific Cox). I² and "
        "Cochran's Q p quantify between-site heterogeneity."
    )
    # Flat-style render (no gtsummary-style sub-level indentation).
    table_rows = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(v), quote=False)}</td>"
                         for v in row) + "</tr>"
        for row in out_df.itertuples(index=False, name=None)
    )
    flat_table_html = (
        '<table class="results-table" border="0">'
        "<thead><tr>"
        + "".join(f"<th>{html.escape(c, quote=False)}</th>" for c in out_df.columns)
        + "</tr></thead><tbody>"
        + table_rows
        + "</tbody></table>"
    )
    html_doc = _MANUSCRIPT_TABLE_HTML_TEMPLATE.format(
        title=("Table 3 — Pooled Causal Estimates: PSM Fine-Gray SHR and IPTW "
               "Cause-Specific HR for the Effect of Initial CRRT Dose on 30-Day "
               "Outcomes"),
        caption=caption,
        table=flat_table_html,
    )
    html_path.write_text(html_doc, encoding="utf-8")

    print(f"  wrote {csv_path.name}")
    print(f"  wrote {html_path.name}")
    return {"csv": csv_path, "html": html_path}


# ── Descriptive epidemiology pooled CSVs (from scripts 07/08) ───────────────

def export_descriptive_pooled_csvs(sites, labels, out_dir: Path) -> None:
    """Pool the per-site descriptive CSVs from 07/08 into manuscript pooled CSVs:
       pooled_crrt_incidence.csv, pooled_crrt_practice_quality.csv,
       pooled_low_dose_counts.csv. Each carries per-site rows plus a POOLED row.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- CRRT incidence (sum numerators/denominators across sites per stratum) ---
    inc_rows = []
    for sd in sites:
        df = _load_csv(sd, "crrt_epi", "crrt_incidence")
        if df is None or df.empty:
            continue
        for _, r in df.iterrows():
            inc_rows.append({
                "site_id": sd.name, "site_label": labels.get(sd.name, sd.name),
                "stratum": r["stratum"], "n_denominator": int(r["n_denominator"]),
                "n_crrt": int(r["n_crrt"]), "incidence_pct": r["incidence_pct"],
            })
    if inc_rows:
        inc = pd.DataFrame(inc_rows)
        pooled = inc.groupby("stratum", as_index=False).agg(
            n_denominator=("n_denominator", "sum"), n_crrt=("n_crrt", "sum"))
        pooled["incidence_pct"] = (100 * pooled["n_crrt"] / pooled["n_denominator"]).round(2)
        pooled["site_id"] = "POOLED"
        pooled["site_label"] = "Pooled"
        out = pd.concat([inc, pooled.reindex(columns=inc.columns)], ignore_index=True)
        out.to_csv(out_dir / "pooled_crrt_incidence.csv", index=False)
        print(f"  wrote pooled_crrt_incidence.csv ({inc['site_id'].nunique()} sites)")

    # --- Practice / quality (stack per-site; pool count-based metrics) ---
    pq_list = []
    for sd in sites:
        df = _load_csv(sd, "crrt_epi", "crrt_practice_quality")
        if df is None or df.empty:
            continue
        df = df.copy()
        df["site_id"] = sd.name
        df["site_label"] = labels.get(sd.name, sd.name)
        pq_list.append(df)
    if pq_list:
        pq = pd.concat(pq_list, ignore_index=True)
        cnt = pq[pq["denominator"].notna()].copy()
        cnt["n"] = pd.to_numeric(cnt["n"], errors="coerce")
        cnt["denominator"] = pd.to_numeric(cnt["denominator"], errors="coerce")
        pooled = cnt.groupby(["variable", "stat"], as_index=False).agg(
            n=("n", "sum"), denominator=("denominator", "sum"))
        pooled["value"] = (100 * pooled["n"] / pooled["denominator"]).round(1)
        pooled["site_id"] = "POOLED"
        pooled["site_label"] = "Pooled"
        out = pd.concat([pq, pooled.reindex(columns=pq.columns)], ignore_index=True)
        out.to_csv(out_dir / "pooled_crrt_practice_quality.csv", index=False)
        print(f"  wrote pooled_crrt_practice_quality.csv ({pq['site_id'].nunique()} sites)")

    # --- Low-dose dose-band counts (sum counts across sites per band) ---
    ld_list = []
    for sd in sites:
        df = _load_csv(sd, "low_dose", "low_dose_counts")
        if df is None or df.empty:
            continue
        df = df.copy()
        df["site_id"] = sd.name
        df["site_label"] = labels.get(sd.name, sd.name)
        ld_list.append(df)
    if ld_list:
        ld = pd.concat(ld_list, ignore_index=True)
        bands = ld[ld["band"] != "Total with dose"]
        total = int(ld.loc[ld["band"] == "Total with dose", "n"].sum())
        pooled = bands.groupby("band", as_index=False)["n"].sum()
        pooled["pct_of_dosed"] = (100 * pooled["n"] / total).round(1) if total else np.nan
        pooled["site_id"] = "POOLED"
        pooled["site_label"] = "Pooled"
        pooled = pd.concat([pooled, pd.DataFrame([{
            "band": "Total with dose", "n": total, "pct_of_dosed": 100.0,
            "site_id": "POOLED", "site_label": "Pooled"}])], ignore_index=True)
        out = pd.concat([ld, pooled.reindex(columns=ld.columns)], ignore_index=True)
        out.to_csv(out_dir / "pooled_low_dose_counts.csv", index=False)
        print(f"  wrote pooled_low_dose_counts.csv ({ld['site_id'].nunique()} sites)")


def export_pooled_low_dose_baseline(sites, labels, out_dir: Path) -> None:
    """Pool the per-site very-low-dose (10-15 mL/kg/hr) vs rest baseline
    comparison across sites -> pooled_low_dose_baseline.csv (tidy: per-site rows
    + POOLED rows) and pooled_low_dose_baseline_table.csv (readable Very-low vs
    Rest, one column each).

    POOLING METHOD = crude / N-weighted, IDENTICAL to the main Table 1
    (_compute_combined_table_df): continuous -> N-weighted median of the per-site
    medians (weights = per-site subgroup n), IQR -> N-weighted median of the
    per-site q25/q75; categorical -> summed n / summed total. WHY: this is a
    descriptive portrait of the pooled patient population, so the patient-
    weighted (crude) estimand is the right one and stays consistent with Table 1,
    the other baseline table a reader compares it to. It is deliberately NOT the
    inverse-variance random-effects method used for pooled_manuscript_numbers /
    pooled_epi_heterogeneity (that answers "how much do sites disagree" and
    carries a CI + I2; a baseline table needs neither). Every row is labeled with
    `pool_method` so the two conventions never get mixed.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    ld_list = []
    for sd in sites:
        df = _load_csv(sd, "low_dose", "low_dose_long")
        if df is None or df.empty:
            continue
        df = df.copy()
        df["site_id"] = sd.name
        df["site_label"] = labels.get(sd.name, sd.name)
        ld_list.append(df)
    if not ld_list:
        return
    ld = pd.concat(ld_list, ignore_index=True)
    for c in ("n", "total", "median", "q25", "q75"):
        ld[c] = pd.to_numeric(ld.get(c), errors="coerce")

    pooled_rows = []
    # continuous -> N-weighted median of site medians (weights = site subgroup n)
    cont = ld[ld["stat_type"] == "continuous"]
    for (var, sub), g in cont.groupby(["variable", "subgroup"], sort=False):
        gg = g[g["median"].notna() & g["n"].notna() & (g["n"] > 0)]
        if gg.empty:
            continue
        w = gg["n"].tolist()
        has_iqr = bool(gg["q25"].notna().all() and gg["q75"].notna().all())
        pooled_rows.append({
            "site_id": "POOLED", "site_label": "Pooled", "variable": var,
            "level": "", "subgroup": sub, "stat_type": "continuous",
            "n": int(gg["n"].sum()), "total": int(gg["n"].sum()), "pct": np.nan,
            "median": round(_weighted_median(gg["median"].tolist(), w), 3),
            "q25": round(_weighted_median(gg["q25"].tolist(), w), 3) if has_iqr else np.nan,
            "q75": round(_weighted_median(gg["q75"].tolist(), w), 3) if has_iqr else np.nan,
            "pval": np.nan, "k_sites": int(gg["site_id"].nunique()),
            "pool_method": "N-weighted median of site medians (crude; = Table 1)",
        })
    # categorical -> summed n / summed total
    cat = ld[ld["stat_type"] == "categorical"]
    for (var, lvl, sub), g in cat.groupby(["variable", "level", "subgroup"], sort=False):
        gg = g[g["n"].notna() & g["total"].notna()]
        if gg.empty:
            continue
        nsum, tsum = int(gg["n"].sum()), int(gg["total"].sum())
        pooled_rows.append({
            "site_id": "POOLED", "site_label": "Pooled", "variable": var,
            "level": lvl, "subgroup": sub, "stat_type": "categorical",
            "n": nsum, "total": tsum,
            "pct": round(100 * nsum / tsum, 1) if tsum else np.nan,
            "median": np.nan, "q25": np.nan, "q75": np.nan, "pval": np.nan,
            "k_sites": int(gg["site_id"].nunique()),
            "pool_method": "summed n / summed total (crude; = Table 1)",
        })

    cols = ["site_id", "site_label", "variable", "level", "subgroup", "stat_type",
            "n", "total", "pct", "median", "q25", "q75", "pval", "k_sites", "pool_method"]
    ld_out = ld.copy()
    ld_out["pct"] = np.where(
        (ld_out["stat_type"] == "categorical") & ld_out["total"].gt(0),
        (100 * ld_out["n"] / ld_out["total"]).round(1), np.nan)
    ld_out["k_sites"] = np.nan
    ld_out["pool_method"] = ""
    tidy = pd.concat([ld_out.reindex(columns=cols),
                      pd.DataFrame(pooled_rows).reindex(columns=cols)], ignore_index=True)
    tidy.to_csv(out_dir / "pooled_low_dose_baseline.csv", index=False)

    # readable wide table: Characteristic | Very low (10-15) | Rest | k sites
    def _cell(r):
        if r["stat_type"] == "continuous":
            if pd.isna(r["median"]):
                return ""
            if pd.notna(r["q25"]) and pd.notna(r["q75"]):
                return f'{r["median"]:g} ({r["q25"]:g}, {r["q75"]:g})'
            return f'{r["median"]:g}'
        return f'{int(r["n"])} ({r["pct"]:g}%)' if pd.notna(r["pct"]) else f'{int(r["n"])}'
    lut, kmap = {}, {}
    for r in pooled_rows:
        key = (r["variable"], r["level"] or "")
        lut[(key, r["subgroup"])] = _cell(r)
        kmap[key] = max(kmap.get(key, 0), r["k_sites"])
    seen, order = set(), []
    for _, rr in ld.iterrows():
        key = (rr["variable"], "" if rr["stat_type"] == "continuous" else rr["level"])
        if key not in seen:
            seen.add(key)
            order.append((key, rr["stat_type"]))
    wide_rows = []
    for key, st in order:
        label = key[0] if st == "continuous" else f"{key[0]}: {key[1]}"
        wide_rows.append({
            "Characteristic": label,
            "Very low (10-15)": lut.get((key, "Very low (10-15)"), ""),
            "Rest": lut.get((key, "Rest"), ""),
            "k sites": kmap.get(key, ""),
        })
    pd.DataFrame(wide_rows).to_csv(out_dir / "pooled_low_dose_baseline_table.csv", index=False)
    print(f"  wrote pooled_low_dose_baseline.csv + _table.csv ({ld['site_id'].nunique()} sites)")


# ── Between-site heterogeneity of the descriptive epi metrics ────────────────
# Federated comparison: each site ships only summaries (counts, means, SDs); the
# pooling + between-site tests live in 07 (compute_epi_heterogeneity, the single
# source of truth, bridged above) so the dashboard can render the same numbers
# inline without reading this CSV. This thin wrapper just persists them.

def export_epi_heterogeneity_csv(sites, labels, out_dir: Path) -> None:
    """Write pooled_epi_heterogeneity.csv (one row per metric) from the shared
    07 computation. labels is accepted for call-site symmetry (unused: metrics
    are site-agnostic)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = compute_epi_heterogeneity(sites)
    if rows:
        df = pd.DataFrame(rows)[EPI_HET_COLS]
        df.to_csv(out_dir / "pooled_epi_heterogeneity.csv", index=False)
        print(f"  wrote pooled_epi_heterogeneity.csv ({len(df)} metrics; "
              f"max k={int(df['k_sites'].max())} sites)")
    else:
        print("  [epi heterogeneity] no per-site epi summaries found; skipping")


def build_supplement_causal_diagnostics(sites, labels, out_dir: Path, pooled_dir: Path = None) -> None:
    """Supplemental causal-model diagnostics, PER SITE, anonymized (Site 1..N) — NOT
    pooled (pool the estimate via the forests; report the diagnostics per site).
    Reads the existing per-site `psm_iptw/` outputs (no re-run needed) and emits
    compact CSV + HTML tables:
      supp_causal_balance      - IPTW covariate balance: max |SMD| unadjusted vs
                                 weighted + # covariates still >0.10 (from SL_IPTW_balance;
                                 the numeric substitute for love plots)
      supp_causal_positivity   - treated arm n/%, effective sample size, PS-tail shares
                                 (from psm_counts_summary + SL_IPTW_ESS + SL_PS_extremes)
      supp_causal_model_comparison - PSM Fine-Gray SHR vs IPTW cause-specific HR (death),
                                 per site (shows the two estimators agree)
      supp_causal_evalue       - E-value (IPTW cause-specific, death) per site, for the
                                 one-line prose sensitivity sentence
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    anon = {sd.name: f"Site {i + 1}" for i, sd in enumerate(sites)}

    def _num(x):
        return pd.to_numeric(x, errors="coerce")

    bal, pos, mc, ev = [], [], [], []
    for sd in sites:
        site = anon[sd.name]

        # Balance for BOTH methods: one shared unadjusted reference + after-PSM
        # (matching) + after-IPTW (weighting). Unadjusted from SL_IPTW_balance
        # Diff.Un; after-IPTW from Diff.Adj; after-PSM from psm_balance_summary's
        # "Std. Mean Diff." (matched-set SMD). PSM typically balances better than
        # IPTW here, so reporting both is more complete (and fairer).
        row = {"Site": site}
        biptw = _load_csv(sd, "psm_iptw", "SL_IPTW_balance")
        if biptw is not None and not biptw.empty:
            cov = biptw.columns[0]
            bb = biptw[~biptw[cov].astype(str).str.contains("prop.score|Distance", case=False, regex=True, na=False)]
            un, adj = _num(bb.get("Diff.Un")).abs(), _num(bb.get("Diff.Adj")).abs()
            if un.notna().any():
                row["Covariates (n)"] = int(un.notna().sum())
                row["Max |SMD| unadjusted"] = round(un.max(), 3)
            if adj.notna().any():
                row["Max |SMD| after IPTW"] = round(adj.max(), 3)
                row["N >0.10 (IPTW)"] = int((adj > 0.10).sum())
        bpsm = _load_csv(sd, "psm_iptw", "psm_balance_summary")
        if bpsm is not None and not bpsm.empty and "Std. Mean Diff." in bpsm.columns:
            smd = _num(bpsm["Std. Mean Diff."]).abs()
            if smd.notna().any():
                row["Max |SMD| after PSM"] = round(smd.max(), 3)
                row["N >0.10 (PSM)"] = int((smd > 0.10).sum())
        if len(row) > 1:
            # order columns: unadjusted, then PSM, then IPTW
            bal.append({k: row.get(k) for k in
                        ["Site", "Covariates (n)", "Max |SMD| unadjusted",
                         "Max |SMD| after PSM", "N >0.10 (PSM)",
                         "Max |SMD| after IPTW", "N >0.10 (IPTW)"] if k in row})

        # Positivity: treated arm, effective sample size, PS tails
        cnt = _load_csv(sd, "psm_iptw", "psm_counts_summary")
        ess = _load_csv(sd, "psm_iptw", "SL_IPTW_ESS")
        ext = _load_csv(sd, "psm_iptw", "SL_PS_extremes")
        pr = {"Site": site}
        if cnt is not None and not cnt.empty:
            a = cnt[cnt.iloc[:, 0].astype(str).str.strip().str.lower() == "all"]
            if len(a):
                t, c = _num(a.iloc[0].get("Treated")), _num(a.iloc[0].get("Control"))
                if pd.notna(t) and pd.notna(c) and (t + c) > 0:
                    pr["High-dose n"] = int(t); pr["High-dose %"] = round(100 * t / (t + c), 1)
        if ess is not None and not ess.empty:
            e = _num(ess.iloc[0].get("ESS"))
            pr["Effective N (weighted)"] = int(round(e)) if pd.notna(e) else None
            pr["ESS proportion"] = round(_num(ess.iloc[0].get("ESS_prop")), 2)
        if ext is not None and not ext.empty:
            def _pct(d):
                r = ext[ext["direction"].astype(str).str.replace(" ", "") == d]
                return round(_num(r.iloc[0]["pct"]), 1) if len(r) else None
            pr["% PS<0.10"], pr["% PS>0.90"] = _pct("<0.10"), _pct(">0.90")
        if len(pr) > 1:
            pos.append(pr)

        # Model comparison: PSM Fine-Gray SHR vs IPTW cause-specific HR (death)
        m = _load_csv(sd, "psm_iptw", "ModelComparison")
        if m is not None and not m.empty:
            def _hr(key):
                r = m[m["model"].astype(str).str.contains(key, case=False, na=False)]
                if not len(r):
                    return "—"
                rr = r.iloc[0]
                return f"{_num(rr['HR']):.2f} ({_num(rr['HR_lower']):.2f}–{_num(rr['HR_upper']):.2f})"
            mc.append({"Site": site, "PSM Fine-Gray SHR (death)": _hr("PSM FG - Death"),
                       "IPTW cause-specific HR (death)": _hr("IPTW Cox - Death")})

        # E-value (primary = IPTW cause-specific, death)
        e = _load_csv(sd, "psm_iptw", "evalue_sensitivity")
        if e is not None and not e.empty:
            r = e[e["model"].astype(str).str.contains("IPTW Cox - Death", case=False, na=False)]
            if len(r) and pd.notna(_num(r.iloc[0].get("evalue_point"))):
                ev.append({"Site": site, "HR (death)": round(_num(r.iloc[0]["HR"]), 2),
                           "E-value (point)": round(_num(r.iloc[0]["evalue_point"]), 2)})

    def _write(rows, stem, caption):
        if not rows:
            print(f"  [{stem}] no per-site data found; skipping")
            return
        df = pd.DataFrame(rows)
        # Integer-valued columns render as whole numbers (no trailing ".0"); a
        # column that is int at some sites but NaN at others (e.g. IPTW counts
        # where a site lacks the SL_IPTW files) is upcast to float by pandas, so
        # coerce whole-number columns to int-string with an em-dash for NaN.
        for c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            if s.notna().any() and (s.dropna() % 1 == 0).all():
                df[c] = s.map(lambda v: "—" if pd.isna(v) else str(int(v)))
        df.to_csv(out_dir / f"{stem}.csv", index=False)
        html = (f'<h3>{caption}</h3>'
                + df.to_html(index=False, na_rep="—", border=0,
                             classes="results-table", justify="center"))
        (out_dir / f"{stem}.html").write_text(html)
        print(f"  wrote {stem}.csv/.html ({len(df)} sites)")

    # Pooled random-effects "diamond" row for the model-comparison table, read
    # from the SAME pooled CSV the main-text Table 2 / forests use (so the numbers
    # match exactly), not re-meta-analyzed here.
    if mc and pooled_dir is not None and (pooled_dir / "pooled_primary_analyses.csv").exists():
        ppa = pd.read_csv(pooled_dir / "pooled_primary_analyses.csv")

        def _re(label):
            r = ppa[ppa["analysis"].astype(str).str.strip() == label]
            if not len(r):
                return "—"
            rr = r.iloc[0]
            return f"{_num(rr['re_hr']):.2f} ({_num(rr['re_lo']):.2f}–{_num(rr['re_hi']):.2f})"
        mc.append({"Site": "Pooled (random-effects)",
                   "PSM Fine-Gray SHR (death)": _re("Point-Treatment PSM FG - Death"),
                   "IPTW cause-specific HR (death)": _re("Point-Treatment IPTW - Death")})

    print("\nSupplement causal diagnostics (per-site, anonymized):")
    _write(bal, "supp_causal_balance", "Covariate balance (IPTW), per site — anonymized")
    _write(pos, "supp_causal_positivity", "Positivity diagnostics, per site — anonymized")
    _write(mc, "supp_causal_model_comparison", "PSM vs IPTW primary estimate (death), per site — anonymized")
    _write(ev, "supp_causal_evalue", "E-value (IPTW cause-specific, death), per site — anonymized")


def build_descriptive_table(sites, labels, out_dir: Path) -> dict:
    """Rojas-style descriptive Table: CRRT incidence + CRRT dose + very-low-dose
    subcohort, one column per site plus a Pooled column. CSV + standalone HTML."""
    out_dir.mkdir(parents=True, exist_ok=True)

    per_site, pool = {}, {"inc": {}, "low": 0, "dosed": 0, "mort": {}}
    site_order = []
    for sd in sites:
        inc = _load_csv(sd, "crrt_epi", "crrt_incidence")
        pq = _load_csv(sd, "crrt_epi", "crrt_practice_quality")
        ld = _load_csv(sd, "low_dose", "low_dose_counts")
        if inc is None or inc.empty:
            continue
        sid = sd.name
        site_order.append(sid)
        per_site[sid] = {"inc": inc, "pq": pq, "ld": ld}
        for strat in ["Overall (adult ICU)", "On invasive ventilation", "On vasopressors"]:
            r = inc[inc["stratum"] == strat]
            if len(r):
                d = pool["inc"].setdefault(strat, {"n": 0, "d": 0})
                d["n"] += int(r.iloc[0]["n_crrt"]); d["d"] += int(r.iloc[0]["n_denominator"])
        if ld is not None and not ld.empty:
            lo = ld[ld["band"] == "10-15 (very low)"]
            tot = ld[ld["band"] == "Total with dose"]
            pool["low"] += int(lo.iloc[0]["n"]) if len(lo) else 0
            pool["dosed"] += int(tot.iloc[0]["n"]) if len(tot) else 0
        if pq is not None and not pq.empty:
            for mvar in ["30-day mortality", "In-hospital mortality"]:
                mr = pq[(pq["variable"] == mvar) & (pq["stat"] == "pct")]
                if len(mr):
                    mm = pool["mort"].setdefault(mvar, {"n": 0, "d": 0})
                    mm["n"] += int(mr.iloc[0]["n"]); mm["d"] += int(mr.iloc[0]["denominator"])
    if not site_order:
        print("  [descriptive table] no per-site CSVs found; skipping")
        return {}

    def inc_cell(d, strat):
        r = d["inc"][d["inc"]["stratum"] == strat]
        if not len(r):
            return "—"
        return f"{r.iloc[0]['incidence_pct']:.1f} ({int(r.iloc[0]['n_crrt']):,}/{int(r.iloc[0]['n_denominator']):,})"

    def pq_cell(d, var, stat):
        if d["pq"] is None:
            return "—"
        r = d["pq"][(d["pq"]["variable"] == var) & (d["pq"]["stat"] == stat)]
        return str(r.iloc[0]["value"]) if len(r) else "—"

    def low_cell(d):
        if d["ld"] is None:
            return "—"
        r = d["ld"][d["ld"]["band"] == "10-15 (very low)"]
        return f"{int(r.iloc[0]['n'])} ({r.iloc[0]['pct_of_dosed']:.1f}%)" if len(r) else "—"

    def pooled_inc(strat):
        s = pool["inc"].get(strat)
        return f"{100*s['n']/s['d']:.1f} ({s['n']:,}/{s['d']:,})" if s and s["d"] else "—"

    def mort_cell(d, var):
        if d["pq"] is None:
            return "—"
        r = d["pq"][(d["pq"]["variable"] == var) & (d["pq"]["stat"] == "pct")]
        if not len(r):
            return "—"
        rr = r.iloc[0]
        return f"{rr['value']}% ({int(rr['n']):,}/{int(rr['denominator']):,})"

    def pooled_mort(var):
        s = pool["mort"].get(var)
        return f"{100*s['n']/s['d']:.1f}% ({s['n']:,}/{s['d']:,})" if s and s["d"] else "—"

    cols = ["Characteristic"] + [labels.get(s, s) for s in site_order]
    multi = len(site_order) > 1
    if multi:
        cols.append("Pooled")

    def row(label, value_fn, pooled_val=""):
        vals = [value_fn(per_site[s]) for s in site_order]
        return [label] + vals + ([pooled_val] if multi else [])

    rows = [
        ["__CRRT incidence, % (n/N)__"] + [""] * (len(cols) - 1),
        row("Overall (adult ICU)", lambda d: inc_cell(d, "Overall (adult ICU)"), pooled_inc("Overall (adult ICU)")),
        row("Among ventilated", lambda d: inc_cell(d, "On invasive ventilation"), pooled_inc("On invasive ventilation")),
        row("Among vasopressor-exposed", lambda d: inc_cell(d, "On vasopressors"), pooled_inc("On vasopressors")),
        ["__CRRT dose (analytic cohort)__"] + [""] * (len(cols) - 1),
        row("Median [IQR], mL/kg/hr", lambda d: pq_cell(d, "CRRT dose (mL/kg/hr)", "median_iqr")),
        row("% <20 mL/kg/hr", lambda d: pq_cell(d, "Dose band", "pct_<20")),
        row("% 20-30 mL/kg/hr", lambda d: pq_cell(d, "Dose band", "pct_20-30 (KDIGO)")),
        row("% >30 mL/kg/hr", lambda d: pq_cell(d, "Dose band", "pct_>30")),
        ["__Very-low-dose subcohort__"] + [""] * (len(cols) - 1),
        row("10-15 mL/kg/hr, n (%)", low_cell,
            f"{pool['low']} ({100*pool['low']/pool['dosed']:.1f}%)" if pool["dosed"] else "—"),
        ["__Outcomes__"] + [""] * (len(cols) - 1),
        row("30-day mortality, % (n/N)", lambda d: mort_cell(d, "30-day mortality"),
            pooled_mort("30-day mortality")),
        row("In-hospital mortality, % (n/N)", lambda d: mort_cell(d, "In-hospital mortality"),
            pooled_mort("In-hospital mortality")),
    ]
    result_df = pd.DataFrame(rows, columns=cols)

    csv_path = out_dir / "pooled_descriptive_epi.csv"
    html_path = out_dir / "pooled_descriptive_epi.html"
    result_df.to_csv(csv_path, index=False)
    caption = (f"CRRT incidence, CRRT dose, and very-low-dose subcohort across "
               f"{len(site_order)} site(s). Incidence is computed at the adult ICU "
               f"hospitalization level: for each stratum it is the number of "
               f"hospitalizations that received CRRT at any point during the "
               f"admission (numerator) divided by all hospitalizations in that "
               f"stratum (denominator). The ventilated and vasopressor strata are "
               f"therefore the proportion of ventilated, or vasopressor-exposed, "
               f"ICU hospitalizations that went on to receive CRRT, not the "
               f"proportion of CRRT recipients who were ventilated or on "
               f"vasopressors. Denominators are hospitalizations (not unique "
               f"patients) and overlap across strata, so they do not sum to 100%. "
               f"Dose metrics are from the analytic cohort.")
    html_doc = _MANUSCRIPT_TABLE_HTML_TEMPLATE.format(
        title="Table — CRRT Epidemiology Across the CLIF Consortium",
        caption=caption,
        table=_render_gtsummary_table_html(result_df),
    )
    html_path.write_text(html_doc, encoding="utf-8")
    print(f"  wrote {csv_path.name}")
    print(f"  wrote {html_path.name}")
    return {"csv": csv_path, "html": html_path}


def build_descriptive_figures(sites, labels, out_dir: Path, *, anonymize=True) -> dict:
    """Cross-site descriptive figures (practice variation). PNG + PDF.
      - pooled_descriptive_incidence: CRRT incidence by site (overall/vent/vaso)
      - pooled_descriptive_dosebands: dose-band composition by site (very-low hi-lit)

    When ``anonymize`` (default) sites are labeled "Site 1..N" per the manuscript
    anonymization rule and the canonical stems are written (consumed by the
    manuscript and the anonymized dashboard). When ``anonymize=False`` the real
    site display names (``labels``) are used and a "_named" variant is written,
    consumed only by the non-anonymized dashboard.
    """
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    BLUE, ORANGE, GREEN = "#1e417c", "#fb801b", "#4CAF50"
    suffix = "" if anonymize else "_named"

    inc_data, band_data = [], []
    inc_names, band_names = [], []  # real site keys, aligned to the filtered data
    for sd in sites:
        inc = _load_csv(sd, "crrt_epi", "crrt_incidence")
        if inc is not None and not inc.empty:
            def g(strat):
                r = inc[inc["stratum"] == strat]
                return float(r.iloc[0]["incidence_pct"]) if len(r) else np.nan
            inc_data.append((g("Overall (adult ICU)"), g("On invasive ventilation"), g("On vasopressors")))
            inc_names.append(sd.name)
        ld = _load_csv(sd, "low_dose", "low_dose_counts")
        if ld is not None and not ld.empty:
            tr = ld.loc[ld["band"] == "Total with dose", "n"]
            tot = float(tr.iloc[0]) if len(tr) else np.nan
            bd = {r["band"]: (100 * float(r["n"]) / tot if tot else np.nan)
                  for _, r in ld.iterrows() if r["band"] != "Total with dose"}
            band_data.append(bd)
            band_names.append(sd.name)

    figs = {}

    # (1) CRRT incidence by site — grouped bars
    if inc_data:
        xlabels = ([f"Site {i+1}" for i in range(len(inc_data))] if anonymize
                   else [labels.get(n, n) for n in inc_names])
        x = np.arange(len(xlabels)); w = 0.26
        fig, ax = plt.subplots(figsize=(max(7, 1.5 * len(xlabels) + 3), 5))
        ax.bar(x - w, [d[0] for d in inc_data], w, label="Overall (adult ICU)", color=BLUE)
        ax.bar(x,     [d[1] for d in inc_data], w, label="Among ventilated", color=ORANGE)
        ax.bar(x + w, [d[2] for d in inc_data], w, label="Among vasopressor-exposed", color=GREEN)
        ax.set_xticks(x); ax.set_xticklabels(xlabels)
        ax.set_ylabel("% of hospitalizations receiving CRRT")
        ax.set_title("CRRT incidence across CLIF consortium sites")
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / f"pooled_descriptive_incidence{suffix}.pdf", bbox_inches="tight")
        fig.savefig(out_dir / f"pooled_descriptive_incidence{suffix}.png", dpi=300, bbox_inches="tight")
        plt.close(fig); figs["incidence"] = True

    # (2) Dose-band composition by site — 100% stacked horizontal bars
    if band_data:
        order = ["<10", "10-15 (very low)", ">15-20", ">20-25", ">25-30", ">30"]
        # Sequential Blues, light (lowest dose) -> dark (highest dose); thin white
        # separators keep the six steps crisp without any single band standing out.
        bcol = {"<10": "#c6dbef", "10-15 (very low)": "#9ecae1", ">15-20": "#6baed6",
                ">20-25": "#4292c6", ">25-30": "#2171b5", ">30": "#08306b"}
        ylabels = ([f"Site {i+1}" for i in range(len(band_data))] if anonymize
                   else [labels.get(n, n) for n in band_names])
        y = np.arange(len(ylabels))
        fig, ax = plt.subplots(figsize=(9, max(3, 0.8 * len(ylabels) + 2)))
        left = np.zeros(len(ylabels))
        for b in order:
            vals = np.array([d.get(b, 0) or 0 for d in band_data])
            ax.barh(y, vals, left=left, color=bcol[b], label=b,
                    edgecolor="white", linewidth=0.5)
            left += vals
        ax.set_yticks(y); ax.set_yticklabels(ylabels)
        ax.set_xlabel("% of dosed encounters")
        ax.set_xlim(0, 100)
        ax.set_title("CRRT dose-band composition across CLIF consortium sites")
        ax.legend(ncol=1, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  title="Dose band (mL/kg/hr)")
        fig.tight_layout()
        fig.savefig(out_dir / f"pooled_descriptive_dosebands{suffix}.pdf", bbox_inches="tight")
        fig.savefig(out_dir / f"pooled_descriptive_dosebands{suffix}.png", dpi=300, bbox_inches="tight")
        plt.close(fig); figs["dosebands"] = True

    print(f"  wrote descriptive figures: {list(figs)}")
    return figs


def build_pooled_dose_distribution(sites, labels, out_dir: Path, *, anonymize=True) -> dict:
    """Cross-site POOLED CRRT dose distribution. The per-site dose_distribution
    exchange CSVs carry counts over shared fixed 2.5 mL/kg/hr bins, so the pooled
    distribution is the summed bin counts (no patient data leaves a site). The
    pooled distribution is drawn as filled bars (% of the pooled cohort); each
    site is overlaid as a thin density outline, anonymized "Site 1..N", so
    between-site practice variation is visible against the consortium pooled
    shape. KDIGO 20-30 band shaded; pooled grouped-data median + pooled band
    percentages annotated. PNG + PDF."""
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    BLUE, GREEN = "#1e417c", "#4CAF50"
    SITE_CYCLE = ["#fb801b", "#6a51a3", "#238b45", "#cb181d", "#2171b5", "#d94801",
                  "#54278f", "#006d2c", "#a50f15", "#08519c"]
    suffix = "" if anonymize else "_named"

    # Per-site: parse the fixed-bin histogram + band counts from the exchange CSV.
    per_site = []          # (bin_pct vector aligned to canonical lefts, n)
    per_site_names = []    # real site keys, aligned to per_site
    bands = {"<20 (low)": 0, "20-30 (guideline)": 0, ">30 (high)": 0}
    lefts = rights = None
    for sd in sites:
        dd = _load_csv(sd, "crrt_epi/graphs", "dose_distribution")
        if dd is None or dd.empty or "variable" not in dd.columns:
            continue
        dd = dd[dd["variable"] == "initial_crrt_dose_ml_kg_hr"]
        hist = dd[dd["kind"] == "hist"]
        # Main bins only: keys are "lo-hi"; skip "<lo"/">hi" under/overflow rows.
        edges_counts = {}
        for _, r in hist.iterrows():
            key = str(r["key"])
            if key[0] in "<>" or "-" not in key:
                continue
            lo, hi = key.split("-")
            edges_counts[(float(lo), float(hi))] = float(r["value"])
        if not edges_counts:
            continue
        order = sorted(edges_counts)
        site_lefts = np.array([lo for lo, _ in order])
        site_rights = np.array([hi for _, hi in order])
        counts = np.array([edges_counts[e] for e in order])
        if lefts is None:
            lefts, rights = site_lefts, site_rights
        elif len(site_lefts) != len(lefts) or not np.allclose(site_lefts, lefts):
            # Shared edges expected; if a site diverges, skip it (don't misalign).
            print(f"    [pooled dose] {sd.name} bin edges differ; skipped")
            continue
        summ = dd[dd["kind"] == "summary"].set_index("key")["value"]
        n = int(summ["n"]) if "n" in summ.index else int(counts.sum())
        per_site.append((counts / n * 100 if n else counts * 0, n))
        per_site_names.append(sd.name)
        bd = dd[dd["kind"] == "band"].set_index("key")["value"]
        for b in bands:
            if b in bd.index:
                bands[b] += float(bd[b])

    if not per_site or lefts is None:
        print("  [pooled dose] no per-site dose_distribution CSVs found; skipping")
        return {}

    width = float(rights[0] - lefts[0])
    edges = np.append(lefts, rights[-1])
    pooled_counts = np.sum([p * n / 100 for p, n in per_site], axis=0)
    pooled_n = int(sum(n for _, n in per_site))
    pooled_pct = pooled_counts / pooled_n * 100 if pooled_n else pooled_counts

    # Pooled median by grouped-data interpolation on the summed histogram.
    cum = np.cumsum(pooled_counts)
    half = cum[-1] / 2.0
    k = min(int(np.searchsorted(cum, half, side="left")), len(pooled_counts) - 1)
    cprev = cum[k - 1] if k > 0 else 0.0
    med = (lefts[k] + (half - cprev) / pooled_counts[k] * width
           if pooled_counts[k] > 0 else float(lefts[k]))
    band_pct = {b: (100 * v / pooled_n if pooled_n else np.nan) for b, v in bands.items()}

    fig, ax = plt.subplots(figsize=(8.5, 5))
    h_span = ax.axvspan(20, 30, color=GREEN, alpha=0.18, zorder=1,
                        label=f"Guideline 20–30: pooled {band_pct['20-30 (guideline)']:.1f}%")
    h_bar = ax.bar(lefts, pooled_pct, width=width, align="edge", color=BLUE, alpha=0.80,
                   edgecolor="black", linewidth=0.4, zorder=2,
                   label=f"Pooled ({len(per_site)} site(s), N={pooled_n:,})")
    # Per-site density outlines, drawn over the pooled bars.
    site_out_labels = ([f"Site {i + 1}" for i in range(len(per_site))] if anonymize
                       else [labels.get(n, n) for n in per_site_names])
    h_sites = [ax.stairs(pct, edges, color=SITE_CYCLE[i % len(SITE_CYCLE)],
                         linewidth=1.3, alpha=0.85, zorder=3, label=site_out_labels[i])
               for i, (pct, _n) in enumerate(per_site)]
    h_med = ax.axvline(med, color="black", linestyle="--", linewidth=1.2, zorder=4,
                       label=f"Pooled median {med:.1f} mL/kg/hr")
    h_lo = ax.plot([], [], " ", label=f"<20: pooled {band_pct['<20 (low)']:.1f}%")[0]
    h_hi = ax.plot([], [], " ", label=f">30: pooled {band_pct['>30 (high)']:.1f}%")[0]
    ax.set_xlim(0, 80)
    ax.set_xlabel("Initial CRRT Dose (mL/kg/hr)")
    ax.set_ylabel("% of encounters")
    ax.set_title("CRRT dose distribution across CLIF consortium sites")
    # Lead with the pooled summary, then per-site outlines.
    ax.legend(handles=[h_bar, h_med, h_span, h_lo, h_hi, *h_sites],
              fontsize=8, ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    fig.savefig(out_dir / f"pooled_descriptive_dose_distribution{suffix}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"pooled_descriptive_dose_distribution{suffix}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote pooled_descriptive_dose_distribution (pooled N={pooled_n:,}, "
          f"{len(per_site)} site(s); pooled median {med:.1f})")
    return {"dose_distribution": True}


def _pool_fixed_hist(sites, pattern: str, variable: str):
    """Sum the fixed-edge histogram bins for `variable` across sites from the
    `pattern` exchange CSV (crrt_epi/graphs). Main bins only (skips '<lo'/'>hi'
    under/overflow rows). Returns (lefts, rights, pooled_counts, pooled_n); the
    first three are None / n is 0 if no site carries it. Sites whose bin edges
    diverge from the first are skipped (never misaligned)."""
    lefts = rights = pooled_counts = None
    pooled_n = 0
    for sd in sites:
        dd = _load_csv(sd, "crrt_epi/graphs", pattern)
        if dd is None or dd.empty or "variable" not in dd.columns:
            continue
        dd = dd[dd["variable"] == variable]
        edges_counts = {}
        for _, r in dd[dd["kind"] == "hist"].iterrows():
            key = str(r["key"])
            if key[0] in "<>" or "-" not in key:
                continue
            lo, hi = key.split("-")
            edges_counts[(float(lo), float(hi))] = float(r["value"])
        if not edges_counts:
            continue
        order = sorted(edges_counts)
        sl = np.array([lo for lo, _ in order])
        sr = np.array([hi for _, hi in order])
        counts = np.array([edges_counts[e] for e in order])
        if lefts is None:
            lefts, rights, pooled_counts = sl, sr, counts.copy()
        elif len(sl) == len(lefts) and np.allclose(sl, lefts):
            pooled_counts = pooled_counts + counts
        else:
            print(f"    [pooled dose-ibw] {sd.name} {variable} bin edges differ; skipped")
            continue
        summ = dd[dd["kind"] == "summary"].set_index("key")["value"]
        pooled_n += int(summ["n"]) if "n" in summ.index else int(counts.sum())
    return lefts, rights, pooled_counts, pooled_n


def build_pooled_dose_ibw(sites, labels, out_dir: Path) -> dict:
    """Cross-site POOLED CRRT dose distribution under actual vs Devine ideal body
    weight, on the identical height-available paired subset (03 exports both
    histograms on that subset with shared fixed edges, so pooled = summed bins).
    Overlays the two pooled distributions (% of the pooled paired cohort), KDIGO
    20-30 shaded, with an annotation box carrying the pooled reclassification
    across the 30 mL/kg/hr threshold from compute_dose_ibw_pooled (single source
    of truth, also written to pooled_dose_ibw_comparison.csv). Median lines/values
    come from that same pooler, so every number matches the CSV + dashboard table.
    Anonymized; safe at k==1 (pooled == the single contributing site). PNG + PDF."""
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    BLUE, ORANGE, GREEN = "#1e417c", "#fb801b", "#4CAF50"

    la, ra, ca, na = _pool_fixed_hist(
        sites, "dose_distribution_actual_paired",
        "initial_crrt_dose_actual_paired_ml_kg_hr")
    lb, _rb, cb, nb = _pool_fixed_hist(
        sites, "dose_distribution_ibw", "initial_crrt_dose_ibw_ml_kg_hr")
    pooled = compute_dose_ibw_pooled(sites)
    if la is None or lb is None or pooled is None:
        print("  [pooled dose-ibw] no per-site dose-by-IBW exchange CSVs; skipping")
        return {}

    width = float(ra[0] - la[0])
    ka = pooled["k_sites"]
    pa = ca / na * 100 if na else ca * 0
    pb = cb / nb * 100 if nb else cb * 0

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axvspan(20, 30, color=GREEN, alpha=0.18, zorder=1,
               label="KDIGO recommendation (20–30)")
    ax.bar(la, pa, width=width, align="edge", color=BLUE, alpha=0.55,
           edgecolor="black", linewidth=0.4, zorder=2,
           label=f"Actual body weight (median {pooled['median_dose_actual']:.1f})")
    ax.bar(lb, pb, width=width, align="edge", color=ORANGE, alpha=0.55,
           edgecolor="black", linewidth=0.4, zorder=2,
           label=f"Ideal body weight, Devine (median {pooled['median_dose_ibw']:.1f})")
    ax.axvline(pooled["median_dose_actual"], color=BLUE, ls=":", lw=1.2, zorder=3)
    ax.axvline(pooled["median_dose_ibw"], color=ORANGE, ls=":", lw=1.2, zorder=3)
    ax.axvline(30, color="black", ls="--", lw=1.0, zorder=4,
               label="High/low cutoff (30)")
    ax.set_xlim(0, 80)
    ax.set_xlabel("CRRT Dose (mL/kg/hr)")
    ax.set_ylabel("% of encounters")
    ax.set_title("CRRT dose: actual vs ideal body weight across CLIF consortium sites")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.text(0.97, 0.55,
            (f"Pooled paired n = {pooled['n_paired']:,} ({ka} site(s))\n"
             f"Median actual/IBW weight ratio: {pooled['median_weight_ibw_ratio']:.2f}\n"
             f"Heavier than ideal: {pooled['pct_heavier_than_ideal']:.0f}%\n"
             f"≥30 mL/kg/hr: actual {pooled['pct_ge30_actual']:.0f}% → "
             f"IBW {pooled['pct_ge30_ibw']:.0f}%\n"
             f"In KDIGO 20–30: actual {pooled['pct_kdigo_actual']:.0f}% → "
             f"IBW {pooled['pct_kdigo_ibw']:.0f}%"),
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#cccccc"))
    fig.tight_layout()
    fig.savefig(out_dir / "pooled_descriptive_dose_ibw.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "pooled_descriptive_dose_ibw.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote pooled_descriptive_dose_ibw (pooled paired N={pooled['n_paired']:,}, "
          f"{ka} site(s); median {pooled['median_dose_actual']:.1f}→"
          f"{pooled['median_dose_ibw']:.1f})")
    return {"dose_ibw": True}


def export_pooled_dose_ibw_csv(sites, labels, out_dir: Path) -> None:
    """Write pooled_dose_ibw_comparison.csv (one pooled row) from
    compute_dose_ibw_pooled — the poolable numerical details behind the
    actual-vs-IBW dose figure (median shift, weight ratio, reclassification %s,
    between-site heterogeneity of % >=30 on IBW). Single source of truth shared
    with the 07 dashboard table."""
    out_dir.mkdir(parents=True, exist_ok=True)
    pooled = compute_dose_ibw_pooled(sites)
    if pooled is None:
        print("  [pooled dose-ibw csv] no per-site comparison CSVs; skipping")
        return
    row = {c: pooled.get(c) for c in DOSE_IBW_POOLED_COLS}
    pd.DataFrame([row]).to_csv(out_dir / "pooled_dose_ibw_comparison.csv", index=False)
    print(f"  wrote pooled_dose_ibw_comparison.csv (k={pooled['k_sites']}, "
          f"n_paired={pooled['n_paired']:,})")


# ── Cross-site pooled descriptive figures: families C (net-UF mortality),
#    D (state over time), B (30-day trajectories) ─────────────────────────────
# Each consumes the per-site exchange CSVs 03 writes under crrt_epi/graphs, pools
# by the method the data supports (summed fixed bins for counts; N-weighted median
# of per-site medians for trajectories), and renders an anonymized consortium
# figure matching the per-site aesthetics. All degrade gracefully at k=1 (pooled
# == the single contributing site).
_PAL_BLUE, _PAL_ORANGE, _PAL_GREEN = "#1e417c", "#fb801b", "#4CAF50"
_TRAJ_MIN_N = 10   # mask pooled time bins with fewer than this many encounters


def build_pooled_uf_mortality(sites, labels, out_dir: Path) -> dict:
    """Family C: pooled crude 30-day mortality vs first-72h net-UF intensity.
    The per-site uf_mortality CSVs share fixed intensity bins, so pooled mortality
    per bin = sum(n_deaths)/sum(n). Marker size scaled to pooled n; Murugan middle
    band (1.01-1.75) shaded. PNG + PDF."""
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    acc, k_sites = {}, 0
    for sd in sites:
        df = _load_csv(sd, "crrt_epi/graphs", "uf_mortality")
        if df is None or df.empty:
            continue
        k_sites += 1
        for _, r in df.iterrows():
            a = acc.setdefault(float(r["bin_mid"]), {"n": 0, "deaths": 0})
            a["n"] += int(r["n"]); a["deaths"] += int(r["n_deaths"])
    if not acc:
        print("  [pooled uf-mortality] no per-site uf_mortality CSVs found; skipping")
        return {}
    mids = np.array(sorted(acc))
    n = np.array([acc[m]["n"] for m in mids], float)
    deaths = np.array([acc[m]["deaths"] for m in mids], float)
    show = n >= 15
    if show.sum() < 3:
        print(f"  [pooled uf-mortality] only {int(show.sum())} bins with pooled n>=15; skipping")
        return {}
    x, mort, nn = mids[show], 100 * deaths[show] / n[show], n[show]
    sizes = 30 + 200 * (nn / nn.max())
    pooled_n = int(n.sum())
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axvspan(1.01, 1.75, color=_PAL_GREEN, alpha=0.20, label="Murugan Middle (1.01–1.75)")
    ax.plot(x, mort, color=_PAL_BLUE, linewidth=1.5, zorder=2)
    ax.scatter(x, mort, s=sizes, color=_PAL_BLUE, alpha=0.85, zorder=3,
               label=f"Pooled mortality ({k_sites} site(s), N={pooled_n:,}; marker size scaled to n)")
    ax.set_xlabel("Net Ultrafiltration Intensity, First 72 h (mL/kg/hr)")
    ax.set_ylabel("Crude 30-Day Mortality (%)")
    ax.set_title("Crude 30-day mortality vs net UF intensity across CLIF consortium sites")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "pooled_descriptive_uf_mortality.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "pooled_descriptive_uf_mortality.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote pooled_descriptive_uf_mortality (pooled N={pooled_n:,}, {k_sites} site(s))")
    return {"uf_mortality": True}


def _pool_state(sites, stem, ncols):
    """Sum per-hour state counts across sites -> pooled proportions per hour.
    ncols: list of n-column names. Returns (hours, {ncol: pct_array}, k_sites)."""
    acc, k = {}, 0
    for sd in sites:
        df = _load_csv(sd, "crrt_epi/graphs", stem)
        if df is None or df.empty:
            continue
        k += 1
        for _, r in df.iterrows():
            a = acc.setdefault(int(r["hour"]), {"tot": 0.0})
            a["tot"] += float(r.get("n_total", 0) or 0)
            for c in ncols:
                a[c] = a.get(c, 0.0) + float(r.get(c, 0) or 0)
    if not acc:
        return None, None, 0
    hours = np.array(sorted(acc))
    pct = {c: np.array([100 * acc[h].get(c, 0) / acc[h]["tot"] if acc[h]["tot"] else 0
                        for h in hours]) for c in ncols}
    return hours, pct, k


def build_pooled_state_figures(sites, labels, out_dir: Path) -> dict:
    """Family D: pooled stacked-area patient-state census over 30 days (IMV state,
    CRRT state). Per-hour state counts summed across sites -> pooled proportions;
    exact per-site stacking order + colors (per-state alpha via RGBA). PNG + PDF."""
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    out_dir.mkdir(parents=True, exist_ok=True)
    specs = [
        ("imv_state_over_crrt", "pooled_descriptive_imv_state",
         "Invasive mechanical ventilation state over 30 days across CLIF consortium sites",
         [("n_imv", "On IMV", _PAL_ORANGE, 0.75),
          ("n_off_imv", "Off IMV", "#9ecae1", 0.85),
          ("n_discharged", "Discharged alive", _PAL_GREEN, 0.60),
          ("n_dead", "Dead", "#9e9e9e", 0.75)]),
        ("crrt_state_over_crrt", "pooled_descriptive_crrt_state",
         "CRRT state over 30 days across CLIF consortium sites",
         [("n_on_crrt", "On CRRT", _PAL_BLUE, 0.80),
          ("n_off_crrt", "Off CRRT (alive; iHD or no RRT)", "#9ecae1", 0.85),
          ("n_discharged", "Discharged alive", _PAL_GREEN, 0.60),
          ("n_dead", "Dead", "#9e9e9e", 0.75)]),
    ]
    figs = {}
    for stem, outstem, title, states in specs:
        hours, pct, k = _pool_state(sites, stem, [c for c, _, _, _ in states])
        if hours is None:
            continue
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.stackplot(hours / 24.0, *[pct[c] for c, _, _, _ in states],
                     colors=[mcolors.to_rgba(col, a) for _, _, col, a in states],
                     labels=[lbl for _, lbl, _, _ in states])
        ax.set_xlim(0, 30); ax.set_ylim(0, 100)
        ax.set_xlabel("Days from CRRT Initiation")
        ax.set_ylabel("Proportion of Patients (%)")
        ax.set_title(title)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / f"{outstem}.pdf", bbox_inches="tight")
        fig.savefig(out_dir / f"{outstem}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        figs[outstem] = k
    print(f"  wrote pooled state figures: {list(figs)}")
    return figs


def _pool_traj(sites, stem, group_col=None):
    """Pool a per-bin trajectory CSV across sites. Per time bin, the pooled central
    line is the N-weighted median of per-site medians (weights = per-site bin n),
    and the ribbon is the N-weighted q25/q75. Returns {group_value: DataFrame}
    with columns hour, p_median, p_q25, p_q75, p_n, k_sites (group_value is None
    when group_col is None)."""
    frames = [df for sd in sites
              if (df := _load_csv(sd, "crrt_epi/graphs", stem)) is not None and not df.empty]
    if not frames:
        return {}
    alldf = pd.concat(frames, ignore_index=True)
    # Fold the older "lactate" and current "lab_lactate" naming to one canonical
    # lab key so every site pools into the same trajectory (older site CSVs would
    # otherwise be silently dropped by the "lab_"-prefixed matching below).
    if group_col == "lab" and group_col in alldf.columns:
        alldf[group_col] = alldf[group_col].astype(str).str.replace(r"^lab_", "", regex=True)
    groups = list(dict.fromkeys(alldf[group_col].dropna())) if group_col else [None]
    out = {}
    for g in groups:
        sub = alldf if g is None else alldf[alldf[group_col] == g]
        rows = []
        for hour, h in sub.groupby("hour"):
            med = pd.to_numeric(h["median"], errors="coerce")
            nn = pd.to_numeric(h["n"], errors="coerce")
            q25 = pd.to_numeric(h["q25"], errors="coerce")
            q75 = pd.to_numeric(h["q75"], errors="coerce")
            ok = med.notna() & nn.notna() & (nn > 0)
            if not ok.any():
                continue
            w = nn[ok].astype(float).tolist()
            rows.append({
                "hour": float(hour),
                "p_median": _weighted_median(med[ok].astype(float).tolist(), w),
                "p_q25": _weighted_median(q25[ok].fillna(med[ok]).astype(float).tolist(), w),
                "p_q75": _weighted_median(q75[ok].fillna(med[ok]).astype(float).tolist(), w),
                "p_n": int(sum(w)), "k_sites": int(ok.sum())})
        if rows:
            out[g] = pd.DataFrame(rows).sort_values("hour").reset_index(drop=True)
    return out


def _pool_census(sites):
    """Pool the per-day 'alive and on CRRT' census across sites into a consortium
    number-at-risk context line (mirrors the per-site grey dashed line in 03).
    Each site exports per-day n (encounters still alive & on CRRT) plus its
    constant baseline n_base; the pooled percentage is (Sum n(day)) / (Sum n_base)
    across sites. Returns a DataFrame [day, pct] or None."""
    frames = [df for sd in sites
              if (df := _load_csv(sd, "crrt_epi/graphs", "crrt_census")) is not None
              and not df.empty and {"day", "n", "n_base"} <= set(df.columns)]
    if not frames:
        return None
    # One n_base per site (constant within a site); total baseline = sum over sites.
    n_base_total = float(sum(pd.to_numeric(f["n_base"], errors="coerce").max() for f in frames))
    if not n_base_total:
        return None
    alldf = pd.concat(frames, ignore_index=True)
    per = (alldf.assign(n=lambda x: pd.to_numeric(x["n"], errors="coerce"))
           .groupby("day", as_index=False)["n"].sum())
    per["pct"] = per["n"] / n_base_total * 100.0
    return per.sort_values("day").reset_index(drop=True)


def _plot_pooled_traj(ax, d, color, ylabel, title, band=None, band_label=None,
                      census=None, census_ylabel=True):
    """Median line + N-weighted IQR ribbon over 30 days; bins with pooled n below
    _TRAJ_MIN_N are dropped. Optional horizontal reference band. When *census* is
    provided (pooled per-day % alive and on CRRT), overlay it as a grey dashed
    number-at-risk context line on a secondary 0-100% axis, mirroring the per-site
    figures (03)."""
    d = d[d["p_n"] >= _TRAJ_MIN_N]
    x = d["hour"] / 24.0
    if band is not None:
        ax.axhspan(band[0], band[1], color=_PAL_GREEN, alpha=0.18, label=band_label)
    ax.fill_between(x, d["p_q25"], d["p_q75"], color=color, alpha=0.18)
    ax.plot(x, d["p_median"], color=color, linewidth=1.8, zorder=3, label="Pooled median [IQR]")
    ax.set_xlim(0, 30)
    ax.set_xlabel("Days from CRRT Initiation")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)
    if census is not None and not census.empty:
        GREY = "#888888"
        ax2 = ax.twinx()
        ax2.plot(census["day"], census["pct"], color=GREY, ls="--", lw=1.4, alpha=0.85,
                 zorder=2, label="Alive and on CRRT (%, Right Axis)")
        ax2.set_ylim(0, 100)
        ax2.tick_params(axis="y", labelcolor=GREY, labelsize=8)
        ax2.set_ylabel("% of Cohort Alive and on CRRT" if census_ylabel else "",
                       color=GREY, fontsize=9)
        if not census_ylabel:
            ax2.set_yticklabels([])
        return ax2
    return None


def _pooled_at_risk_row(ax, census, days=(0, 5, 10, 15, 20, 25, 30)):
    """KM-style number-at-risk row beneath the x-axis: pooled N (and %) alive & on
    CRRT at each day. Mirrors the per-site _at_risk_row (03) so the pooled figures
    read the same. Draw on the PRIMARY axis (data-x = days)."""
    if census is None or census.empty:
        return
    GREY = "#888888"
    idx = census.set_index("day")
    tr = ax.get_xaxis_transform()  # x in data (days), y in axes fraction
    ax.text(0.0, -0.155, "Alive and on CRRT, N (%):", transform=ax.transAxes,
            ha="left", va="top", fontsize=8, color=GREY, fontweight="bold")
    for d in days:
        if d in idx.index:
            t = f"{int(idx.loc[d, 'n'])}\n({idx.loc[d, 'pct']:.0f}%)"
        else:
            t = "0\n(0%)"
        ax.text(d, -0.215, t, transform=tr, ha="center", va="top", fontsize=7.5, color=GREY)


def build_pooled_trajectories(sites, labels, out_dir: Path) -> dict:
    """Family B: pooled 30-day trajectories (CRRT dose, net-UF rate + cumulative,
    NEE single figures + a 5-lab panel). Per-bin N-weighted median of per-site
    medians; N-weighted IQR ribbon; reference bands where applicable. PNG + PDF."""
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)
    figs = {}
    census = _pool_census(sites)  # pooled % alive-and-on-CRRT number-at-risk line
    single = [
        ("crrt_dose_hourly", "pooled_descriptive_traj_dose", _PAL_BLUE,
         "CRRT Dose (mL/kg/hr)", "Pooled CRRT dose over 30 days across CLIF consortium sites",
         (20, 30), "KDIGO Recommendation (20–30)"),
        ("net_uf_rate_hourly", "pooled_descriptive_traj_uf_rate", _PAL_BLUE,
         "Net Ultrafiltration Rate (mL/kg/hr)",
         "Pooled net UF rate over 30 days across CLIF consortium sites",
         (1.01, 1.75), "Murugan Middle (1.01–1.75)"),
        ("net_uf_cumulative", "pooled_descriptive_traj_uf_cumulative", _PAL_BLUE,
         "Cumulative Net Ultrafiltration (mL/kg)",
         "Pooled cumulative net UF over 30 days across CLIF consortium sites", None, None),
        ("nee_over_crrt", "pooled_descriptive_traj_nee", _PAL_ORANGE,
         "Norepinephrine Equivalents (mcg/kg/min)",
         "Pooled vasopressor (NEE) over 30 days across CLIF consortium sites", None, None),
    ]
    for stem, outstem, color, ylabel, title, band, band_label in single:
        d = _pool_traj(sites, stem).get(None)
        if d is None or d.empty:
            continue
        fig, ax = plt.subplots(figsize=(11, 5.5))
        ax2 = _plot_pooled_traj(ax, d, color, ylabel, title, band, band_label,
                                census=census)
        handles, labs = ax.get_legend_handles_labels()
        if ax2 is not None:
            h2, l2 = ax2.get_legend_handles_labels()
            handles += h2; labs += l2
        ax.legend(handles, labs, fontsize=8, loc="upper right")
        have_census = census is not None and not census.empty
        if have_census:
            _pooled_at_risk_row(ax, census)
        # Leave room at the bottom for the number-at-risk row when it is drawn.
        fig.tight_layout(rect=[0, 0.08, 1, 1] if have_census else None)
        fig.savefig(out_dir / f"{outstem}.pdf", bbox_inches="tight")
        fig.savefig(out_dir / f"{outstem}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        figs[outstem] = True

    lab_pool = _pool_traj(sites, "lab_distributions_over_crrt", group_col="lab")
    lab_meta = [("lactate", "Lactate (mmol/L)"), ("ph_arterial", "Arterial pH"),
                ("bicarbonate", "Bicarbonate (mEq/L)"), ("potassium", "Potassium (mEq/L)"),
                ("phosphate", "Phosphate (mg/dL)")]
    present = [(c, lbl) for c, lbl in lab_meta if c in lab_pool and not lab_pool[c].empty]
    if present:
        fig, axes = plt.subplots(3, 2, figsize=(14, 13))
        axes = axes.flatten()
        for i, (c, lbl) in enumerate(present):
            _plot_pooled_traj(axes[i], lab_pool[c], _PAL_GREEN, lbl, lbl,
                              census=census, census_ylabel=False)
        # 6th panel = pooled census reference (mirrors the per-site bottom-right
        # panel; slot 5 is reserved for it even when fewer than 5 labs are present).
        axc = axes[5]
        if census is not None and not census.empty:
            axc.plot(census["day"], census["pct"], color="#888888", ls="--", lw=1.8)
            axc.set_xlim(0, 30); axc.set_ylim(0, 100)
            axc.set_xlabel("Days from CRRT Initiation")
            axc.set_ylabel("% of Cohort"); axc.set_title("Alive and on CRRT (%)")
            axc.grid(axis="y", alpha=0.25)
            _pooled_at_risk_row(axc, census)
        else:
            axc.set_visible(False)
        for j in range(len(present), 5):  # hide gap panels, keep slot 5 = census
            axes[j].set_visible(False)
        have_census = census is not None and not census.empty
        fig.suptitle("Pooled lab trajectories over 30 days across CLIF consortium sites"
                     + ("\n(Dashed Grey: % of Cohort Alive and on CRRT; See Bottom-Right Panel)"
                        if have_census else ""),
                     fontsize=14, y=1.0)
        fig.tight_layout()
        fig.savefig(out_dir / "pooled_descriptive_traj_labs.pdf", bbox_inches="tight")
        fig.savefig(out_dir / "pooled_descriptive_traj_labs.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        figs["traj_labs"] = True

    print(f"  wrote pooled trajectory figures: {list(figs)}")
    return figs


# ── Risk-standardized mortality (SMR) ──────────────────────────────────────
def export_pooled_smr_csvs(sites, labels, out_dir: Path) -> None:
    """Pooled SMR CSVs: per-site O/E rows + pooled ΣO/ΣE (Byar CI)."""
    per_site, pooled = collect_smr(sites, labels)
    if not per_site:
        print("  (no SMR data; skipping pooled SMR CSVs)")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(per_site).to_csv(out_dir / "pooled_smr_per_site.csv", index=False)
    if pooled:
        pd.DataFrame([pooled]).to_csv(out_dir / "pooled_smr.csv", index=False)
    print(f"  wrote pooled_smr.csv + pooled_smr_per_site.csv ({len(per_site)} site(s))")


# ---------------------------------------------------------------------------
# Site eligibility for POOLED CAUSAL estimates (S-Methods 4, site eligibility)
# ---------------------------------------------------------------------------
# A site is excluded from the pooled CAUSAL analysis if fewer than this many of
# its patients received the higher dose. Two grounds, kept distinct: the arm is
# too small to fit the propensity, outcome and diagnostic models; and there is
# clinical suspicion that higher-dose CRRT is effectively never delivered there,
# i.e. a structural positivity violation.
#
# SCOPE, AND THIS IS THE EASY THING TO GET WRONG: the rule applies ONLY to
# pooled causal estimates. Excluded sites remain in Table 1, the SMR, incidence,
# practice variation, the low-dose subcohort, and every per-site diagnostic
# table. Dropping them from the diagnostic tables would hide the very evidence
# the rule rests on.
MIN_HIGH_DOSE_ARM = 100


def split_causal_eligible(sites, labels):
    """(eligible_sites, excluded_rows) under the MIN_HIGH_DOSE_ARM rule.

    Arm sizes come from the same per-site parse that feeds per_site_summary.csv,
    so the rule is applied to exactly the numbers reported (_read_table1_psm_counts).
    A site whose arm size
    cannot be determined is RETAINED and flagged, because silently dropping a
    site for a missing file would be an unannounced exclusion.
    """
    eligible, excluded = [], []
    for sd in sites:
        n_high = None
        t1 = _read_table1_psm_counts(sd)
        if t1:
            n_high = t1.get("n_high_dose")
        if n_high is None:
            print(f"  [eligibility] {sd.name}: arm size undetermined; RETAINED")
            eligible.append(sd)
            continue
        if int(n_high) < MIN_HIGH_DOSE_ARM:
            excluded.append({"site": sd.name, "label": labels.get(sd.name, sd.name),
                             "n_high_dose": int(n_high),
                             "n_total": t1.get("n_total")})
        else:
            eligible.append(sd)
    return eligible, excluded


def _write_supp_table(df, stem: str, caption: str, out_dir: Path):
    """Shared writer for the supplemental tables rendered straight from collected
    per-site CSVs. Mirrors the supp_causal_* convention: CSV plus a standalone
    HTML fragment with the same caption."""
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"{stem}.csv", index=False)
    html = (f"<h3>{caption}</h3>"
            + df.to_html(index=False, na_rep="\u2014", border=0,
                         classes="results-table", justify="center"))
    (out_dir / f"{stem}.html").write_text(html)
    print(f"  wrote {stem}.csv/.html ({len(df)} rows)")


# Covariates reported in the missingness table: the model covariates that are
# ever incomplete. Age, sex, weight, IMV and vasopressor dose are complete at
# every site and are stated as such in the caption rather than given empty rows.
_MISSINGNESS_VARS = [
    ("lactate_0", "Lactate"),
    ("crrt_dose_median_3h", "CRRT dose"),
    ("bicarbonate_0", "Bicarbonate"),
    ("potassium_0", "Potassium"),
    ("pf_sf_ratio_0", "P/F or S/F ratio"),
]


def build_missingness_table(sites, labels, out_dir: Path) -> dict:
    """Missingness of the model covariates -> tables/supp_missingness.{csv,html}.

    Pooled percent missing is sum(n_missing)/sum(n_total), i.e. patient-weighted,
    NOT a mean of per-site percentages. The between-site range is the span of
    per-site percentages, matching how every other between-site range in the
    paper is defined.
    """
    frames = []
    for sd in sites:
        df = _load_csv(sd, "diagnostics", "missingness_summary")
        if df is not None and not df.empty and "variable" in df.columns:
            frames.append(df.assign(_site=sd.name))
    if not frames:
        print("  (no missingness_summary.csv collected; skipping)")
        return {}
    allm = pd.concat(frames, ignore_index=True)

    rows = []
    for var, label in _MISSINGNESS_VARS:
        sub = allm[allm["variable"] == var]
        if sub.empty:
            continue
        tot, miss = sub["n_total"].sum(), sub["n_missing"].sum()
        pct = sub["pct_missing"].astype(float)
        rows.append({
            "Covariate": label,
            "Sites reporting": int(sub["_site"].nunique()),
            "Pooled % missing": f"{100.0 * miss / tot:.1f}" if tot else "\u2014",
            "Between-site range (%)": f"{pct.min():.1f}\u2013{pct.max():.1f}",
        })
    if not rows:
        print("  (missingness_summary.csv present but no expected covariates; skipping)")
        return {}
    _write_supp_table(pd.DataFrame(rows), "supp_missingness",
                      "Missingness of model covariates, pooled across sites. Pooled percent is "
                      "patient-weighted (summed missing over summed denominator); the range is the span "
                      "of per-site percentages. Age, sex, body weight, mechanical ventilation status and "
                      "vasopressor dose were complete at every site and are omitted.",
                      out_dir)
    return {"n_rows": len(rows)}


def build_modality_table(sites, labels, out_dir: Path) -> dict:
    """CRRT modality mix and mapping completeness per site ->
    tables/supp_modality.{csv,html}.

    Reads the long-format practice-quality export, where modality appears as
    variable="Modality" with stat="pct_<LEVEL>". Levels are NOT fixed: a site
    that never records a modality simply has no row for it, so the column set is
    discovered from the data rather than assumed.

    The Unknown column is a DATA COMPLETENESS measure, not a clinical category:
    it is the share of encounters whose modality is null at the initiation
    record. Reported per site because it varies by an order of magnitude and a
    pooled figure would hide that.
    """
    per_site, levels = [], []
    for i, sd in enumerate(sites):
        pq = _load_csv(sd, "crrt_epi", "crrt_practice_quality")
        if pq is None or pq.empty or "variable" not in pq.columns:
            continue
        mod = pq[pq["variable"].astype(str).str.strip().str.lower() == "modality"]
        if mod.empty:
            continue
        denom = pd.to_numeric(mod["denominator"], errors="coerce").max()
        vals = {}
        for _, r in mod.iterrows():
            lvl = str(r["stat"]).replace("pct_", "").strip().upper()
            vals[lvl] = pd.to_numeric(r["value"], errors="coerce")
            if lvl not in levels:
                levels.append(lvl)
        per_site.append({"Site": f"Site {i + 1}",
                         "N": int(denom) if pd.notna(denom) else None, "_v": vals})
    if not per_site:
        print("  (no crrt_practice_quality.csv modality rows; skipping)")
        return {}

    # Known modalities first in clinical order, then anything else, Unknown last.
    order = [l for l in ["CVVHD", "CVVH", "CVVHDF", "SCUF"] if l in levels]
    order += [l for l in levels if l not in order and l != "UNKNOWN"]
    if "UNKNOWN" in levels:
        order.append("UNKNOWN")

    rows = []
    for ps in per_site:
        row = {"Site": ps["Site"], "N with a modality record": ps["N"]}
        for lvl in order:
            v = ps["_v"].get(lvl)
            key = "Unknown (%)" if lvl == "UNKNOWN" else f"{lvl} (%)"
            row[key] = "\u2014" if v is None or pd.isna(v) else f"{float(v):.1f}"
        rows.append(row)
    _write_supp_table(pd.DataFrame(rows), "supp_modality",
                      "CRRT modality distribution and mapping completeness by site \u2014 anonymized. "
                      "Percentages are of encounters with a modality record at initiation; a dash means "
                      "the site recorded no encounters of that modality. The Unknown column is a measure "
                      "of mapping completeness rather than a clinical category.",
                      out_dir)
    return {"n_sites": len(rows)}


# Minimum cell count below which a per-site figure is suppressed in the supplement.
_SUPP_MIN_N = 10


def build_crrt_settings_by_mode_table(sites, labels, out_dir: Path) -> dict:
    """S-Table 1 — per-site CRRT settings at initiation by modality ->
    tables/supp_crrt_settings_by_mode.{csv,html}.

    The per-site source (diagnostics/{SITE}_crrt_initial_settings_by_mode.csv)
    already carries median [IQR] strings per mode, so this is a per-site side-by-
    side table (sites anonymized Site 1..N), NOT a meta-analytic pool — raw values
    never leave a site, so pooled medians cannot be recomputed centrally. Modes
    with fewer than _SUPP_MIN_N patients are small-cell suppressed: the mode is
    kept (so the reader sees it exists) but its N shows "<N" and its settings blank.
    """
    SRC = [("blood_flow_rate", "Blood Flow Rate (mL/min)"),
           ("dialysate_flow_rate", "Dialysate (mL/hr)"),
           ("pre_filter_replacement_fluid_rate", "Pre-filter Replacement (mL/hr)"),
           ("post_filter_replacement_fluid_rate", "Post-filter Replacement (mL/hr)"),
           ("ultrafiltration_out", "Ultrafiltration (mL/hr)"),
           ("crrt_dose_ml_kg_hr", "Delivered Dose (mL/kg/hr)")]
    rows = []
    for i, sd in enumerate(sites):
        df = _load_csv(sd, "diagnostics", "crrt_initial_settings_by_mode")
        if df is None or df.empty or "Mode" not in df.columns:
            continue
        site = f"Site {i + 1}"
        for _, r in df.iterrows():
            n = pd.to_numeric(r.get("N_patients"), errors="coerce")
            small = pd.notna(n) and n < _SUPP_MIN_N
            row = {"Site": site, "Modality": str(r["Mode"]).upper(),
                   "N": f"<{_SUPP_MIN_N}" if small else (int(n) if pd.notna(n) else "—")}
            for col, lbl in SRC:
                row[lbl] = "—" if (small or col not in df.columns) else r.get(col, "—")
            rows.append(row)
    if not rows:
        print("  (no crrt_initial_settings_by_mode.csv collected; skipping)")
        return {}
    _write_supp_table(pd.DataFrame(rows), "supp_crrt_settings_by_mode",
                      "CRRT settings at initiation by modality, per site (anonymized). Each cell is the "
                      "median [IQR] of the first recorded setting per patient. Sites are shown side by side "
                      "rather than pooled because only per-site summaries are shared. Modes with fewer than "
                      f"{_SUPP_MIN_N} patients are shown with settings suppressed.",
                      out_dir)
    return {"n_rows": len(rows)}


_SEVERITY_METRIC_ORDER = ["potassium_t1", "lactate_t1", "bicarbonate_t1",
                          "creatinine_t1", "bun_t1"]


def build_severity_thresholds_table(sites, labels, out_dir: Path) -> dict:
    """S-Table 4 — baseline severity thresholds by site ->
    tables/supp_severity_thresholds.{csv,html}.

    One row per severe-derangement threshold (K>6, lactate>4, bicarbonate<15,
    creatinine>4, BUN>80). Shows each site's percent meeting the threshold
    (anonymized), a patient-weighted pooled percent (Σ n_meeting / Σ n_with_lab,
    NOT a mean of site percentages), and the between-site range. Per-site cells
    with a numerator below _SUPP_MIN_N are suppressed to "<N".
    """
    per_site = []           # list of (label, {metric_col: (n_meet, n_lab, pct, threshold, metric_label)})
    for i, sd in enumerate(sites):
        df = _load_csv(sd, "crrt_epi", "baseline_severity_thresholds")
        if df is None or df.empty or "column" not in df.columns:
            continue
        d = {}
        for _, r in df.iterrows():
            d[str(r["column"])] = (
                pd.to_numeric(r.get("n_meeting"), errors="coerce"),
                pd.to_numeric(r.get("n_with_lab"), errors="coerce"),
                pd.to_numeric(r.get("pct_meeting"), errors="coerce"),
                str(r.get("threshold", "")), str(r.get("metric", "")))
        per_site.append((f"Site {i + 1}", d))
    if not per_site:
        print("  (no baseline_severity_thresholds.csv collected; skipping)")
        return {}

    # metric order: canonical list first, then any extras present in the data
    seen = [m for m in _SEVERITY_METRIC_ORDER if any(m in d for _, d in per_site)]
    for _, d in per_site:
        for m in d:
            if m not in seen:
                seen.append(m)
    site_cols = [lbl for lbl, _ in per_site]
    rows = []
    for m in seen:
        label = threshold = ""
        num_tot = den_tot = 0.0
        pcts = []
        rec = {}
        for lbl, d in per_site:
            if m not in d:
                rec[lbl] = "—"
                continue
            n_meet, n_lab, pct, thr, met = d[m]
            label = label or met
            threshold = threshold or thr
            if pd.notna(n_meet) and pd.notna(n_lab):
                num_tot += float(n_meet); den_tot += float(n_lab)
            if pd.notna(pct):
                pcts.append(float(pct))
            rec[lbl] = (f"<{_SUPP_MIN_N}" if (pd.notna(n_meet) and n_meet < _SUPP_MIN_N)
                        else (f"{float(pct):.1f}" if pd.notna(pct) else "—"))
        row = {"Severity threshold": label, "Threshold": threshold}
        for lbl in site_cols:
            row[f"{lbl} (%)"] = rec.get(lbl, "—")
        row["Pooled % (patient-weighted)"] = f"{100.0 * num_tot / den_tot:.1f}" if den_tot else "—"
        row["Between-site range (%)"] = f"{min(pcts):.1f}–{max(pcts):.1f}" if pcts else "—"
        rows.append(row)
    _write_supp_table(pd.DataFrame(rows), "supp_severity_thresholds",
                      "Prevalence of severe baseline derangements at CRRT initiation, by site (anonymized). "
                      "Each cell is the percent of the site's patients with a recorded value meeting the "
                      "threshold. The pooled column is patient-weighted (summed numerator over summed "
                      "denominator); the range is the span of per-site percentages. Captures the severe tail "
                      "that the Table 1 medians hide.",
                      out_dir)
    return {"n_rows": len(rows)}


def build_ph_diagnostics_table(sites, labels, out_dir: Path, pooled_dir: Path) -> dict:
    """Proportional-hazards and discrimination diagnostics ->
    tables/supp_ph_diagnostics.{csv,html}.

    Rendered from the already-pooled cox-diagnostics CSVs rather than re-collected.
    Reports, per model, how many sites fail the global test and how many fail for
    the PRIMARY TERM (treatment or dose). That distinction is the whole point of
    the table: global failures are driven by covariates and do not threaten the
    reported effect, whereas a primary-term failure would.
    """
    src = pooled_dir / "pooled_cox_diagnostics.csv"
    if not src.exists():
        print("  (pooled_cox_diagnostics.csv missing; skipping PH table)")
        return {}
    d = pd.read_csv(src)
    rows = []
    for model, g in d.groupby("model", sort=False):
        k = len(g)
        gl = (g["global_ph_p"] < 0.05).sum()
        pt = (g["primary_term_ph_p"] < 0.05).sum()
        rows.append({
            "Model": model,
            "Sites": k,
            "Global PH violated": f"{gl} / {k}",
            "Primary term": str(g["primary_term"].iloc[0]),
            "Primary-term PH violated": f"{pt} / {int(g['primary_term_ph_p'].notna().sum())}",
            "Min primary-term p": f"{g['primary_term_ph_p'].min():.3f}",
            "C-statistic, median (range)":
                f"{g['c_statistic'].median():.3f} ({g['c_statistic'].min():.3f}\u2013{g['c_statistic'].max():.3f})",
        })
    _write_supp_table(pd.DataFrame(rows), "supp_ph_diagnostics",
                      "Proportional-hazards and discrimination diagnostics by model, across sites. "
                      "PH assessed by scaled Schoenfeld residual tests at \u03b1 = 0.05. The primary term "
                      "is treatment for the weighted cause-specific models and dose for the "
                      "dose-response models; global tests additionally pool departures across all "
                      "covariates. Per-covariate results are in pooled_cox_diagnostics_per_covariate.csv.",
                      out_dir)
    return {"n_models": len(rows)}


def build_iptw_cif_figure(sites, out_dir: Path) -> dict:
    """Pooled IPTW-standardized CIFs, death AND discharge, two panels ->
    figures/pooled_iptw_cif.{png,pdf}.  (Supplemental figure.)

    Companion to the PSM cumulative-incidence figure. Same pooling contract, and
    it must be described the same way: each site's standardized curve is
    interpolated onto a common 30-day grid, the pooled curve is the UNWEIGHTED
    arithmetic mean across sites, and the band is a t-interval on the between-site
    standard deviation. The band is therefore between-site dispersion, NOT a
    patient-level confidence interval, and sites contribute equally regardless of
    size. Per-site bootstrap intervals exist in the source CSVs and are
    deliberately not propagated.

    Both competing events are shown because reporting only death would leave a
    reader unable to tell whether a dose effect on death was offset by an effect
    on discharge.
    """
    import matplotlib.pyplot as plt
    from scipy.stats import t as _t_dist

    TIME_MAX, n_grid = 30.0, 300
    tgrid = np.linspace(0, TIME_MAX, n_grid)
    arms = [("<30", "Lower dose (<30 mL/kg/hr)", "#1e417c"),
            (">=30", "Higher dose (>=30 mL/kg/hr)", "#b2182b")]
    outcomes = [("cif_death", "Death"), ("cif_disch", "Discharge alive")]
    curves = {(o, a): [] for o, _ in outcomes for a, _, _ in arms}

    n_sites = 0
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "IPTW_CIF_data")
        if df is None or df.empty or "trt_group" not in df.columns:
            continue
        n_sites += 1
        for arm, _, _ in arms:
            sub = df[df["trt_group"].astype(str).str.strip() == arm].sort_values("time")
            if sub.empty:
                continue
            for col, _ in outcomes:
                if col not in sub.columns:
                    continue
                curves[(col, arm)].append(
                    np.interp(tgrid, sub["time"].to_numpy(), sub[col].to_numpy(),
                              left=0.0, right=np.nan))
    if n_sites < 2:
        print(f"  [iptw cif] need >=2 sites with IPTW_CIF_data.csv (found {n_sites}); skipping")
        return {}

    def _pooled_ci(arr, conf=0.95):
        mean = np.nanmean(arr, axis=0)
        sd = np.nanstd(arr, axis=0, ddof=1)
        n = (~np.isnan(arr)).sum(axis=0)
        valid = n > 1
        t_crit = np.where(valid, _t_dist.ppf((1 + conf) / 2, df=np.where(valid, n - 1, 1)), np.nan)
        se = np.where(valid, sd / np.sqrt(np.maximum(n, 1)), np.nan)
        return mean, mean - t_crit * se, mean + t_crit * se

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=False)
    for ax, (col, oname) in zip(axes, outcomes):
        for arm, alabel, colour in arms:
            stack = curves[(col, arm)]
            if not stack:
                continue
            mean, lo, hi = _pooled_ci(np.vstack(stack))
            ok = ~np.isnan(mean)
            band = ok & ~np.isnan(lo) & ~np.isnan(hi)
            ax.fill_between(tgrid[band], np.clip(lo[band], 0, 1), np.clip(hi[band], 0, 1),
                            color=colour, alpha=0.18, linewidth=0)
            ax.plot(tgrid[ok], mean[ok], color=colour, lw=2, label=alabel)
        ax.set_title(f"Cumulative incidence: {oname}")
        ax.set_xlabel("Days from CRRT initiation")
        ax.set_ylabel("Cumulative incidence")
        ax.set_xlim(0, TIME_MAX)
        ax.set_ylim(0, None)
        ax.grid(alpha=0.25)
    axes[0].legend(loc="upper left", frameon=False, fontsize=9)
    fig.suptitle(f"IPTW-standardized cumulative incidence by initial CRRT dose "
                 f"(pooled across {n_sites} sites)", y=1.02)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "pooled_iptw_cif.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "pooled_iptw_cif.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote pooled_iptw_cif.{{png,pdf}} ({n_sites} sites, both competing events)")
    return {"n_sites": n_sites}


def build_smr_per_site_table(sites, labels, out_dir: Path) -> dict:
    """Per-site SMR provenance table -> tables/supp_smr_per_site.{csv,html}.

    Numeric backing for the main-text SMR forest (F5). That figure prints each
    site's SMR and CI but not the counts the ratio is built from, so a reader
    cannot check O/E, judge whether a wide interval reflects few deaths or
    genuine noise, or see how unequally sites contribute.

    Site labels and row order are taken from the SAME collect_smr() call the
    forest uses, so "Site 1" means the same institution in both. Do not sort
    this table independently.
    """
    per_site, pooled = collect_smr(sites, labels)
    if not per_site:
        print("  (no SMR data; skipping per-site SMR table)")
        return {}
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for i, p in enumerate(per_site):
        auc = p.get("auc")
        crude = p.get("crude_mort_pct")
        rows.append({
            "Site": f"Site {i + 1}",
            "N": int(p["n"]),
            "Observed deaths": int(p["observed"]),
            "Expected deaths": f"{p['expected']:.0f}",
            "SMR (95% CI)": f"{p['smr']:.2f} ({p['lo']:.2f}–{p['hi']:.2f})",
            "AUC": "—" if auc is None or pd.isna(auc) else f"{auc:.3f}",
            "Crude mortality (%)": "—" if crude is None or pd.isna(crude) else f"{crude:.1f}",
        })
    if pooled and pooled.get("k", 0) > 1:
        sN = sum(int(p["n"]) for p in per_site)
        sO = sum(int(p["observed"]) for p in per_site)
        sE = sum(float(p["expected"]) for p in per_site)
        rows.append({
            "Site": f"Pooled (k={pooled['k']})",
            "N": sN,
            "Observed deaths": sO,
            "Expected deaths": f"{sE:.0f}",
            "SMR (95% CI)": f"{pooled['smr']:.2f} ({pooled['lo']:.2f}–{pooled['hi']:.2f})",
            # AUC is deliberately blank: discrimination does not pool by summing,
            # and an unweighted mean of site AUCs corresponds to no quantity.
            "AUC": "—",
            "Crude mortality (%)": f"{100.0 * sO / sN:.1f}" if sN else "—",
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "supp_smr_per_site.csv", index=False)
    caption = ("Risk-standardized 30-day mortality by site — anonymized. "
               "Expected deaths are from the frozen external reference model; "
               "AUC is that model's discrimination when applied at each site. "
               "Site numbering matches the SMR forest figure.")
    html = (f"<h3>{caption}</h3>"
            + df.to_html(index=False, na_rep="—", border=0,
                         classes="results-table", justify="center"))
    (out_dir / "supp_smr_per_site.html").write_text(html)
    print(f"  wrote supp_smr_per_site.csv/.html ({len(per_site)} sites"
          f"{' + pooled' if len(df) > len(per_site) else ''})")
    return {"n_sites": len(per_site)}


def build_smr_forest_figure(sites, labels, out_dir: Path) -> dict:
    """Anonymized (Site 1..N) SMR forest + pooled diamond; reference line at 1.0."""
    import matplotlib.pyplot as plt
    per_site, pooled = collect_smr(sites, labels)
    if not per_site:
        print("  (no SMR data; skipping SMR forest)")
        return {}
    out_dir.mkdir(parents=True, exist_ok=True)
    smr = [p["smr"] for p in per_site]; lo = [p["lo"] for p in per_site]; hi = [p["hi"] for p in per_site]
    show_pool = bool(pooled) and pooled["k"] > 1
    m = len(per_site) + (1 if show_pool else 0)
    yy = list(range(m))[::-1]
    ylabels = [f"Site {i+1}" for i in range(len(per_site))] + (["Pooled"] if show_pool else [])
    fig, ax = plt.subplots(figsize=(8, max(2.4, 0.6 * m + 1.4)))
    for i in range(len(per_site)):
        ax.errorbar(smr[i], yy[i], xerr=[[smr[i] - lo[i]], [hi[i] - smr[i]]], fmt="o",
                    color="#1e417c", ecolor="#1e417c", capsize=4, markersize=7)
    xmax = max(hi + ([pooled["hi"]] if show_pool else []))
    if show_pool:
        yp = yy[-1]
        ax.fill([pooled["lo"], pooled["smr"], pooled["hi"], pooled["smr"]],
                [yp, yp - 0.3, yp, yp + 0.3], color="#b2182b", alpha=0.9)
    ax.axvline(1.0, ls="--", color="black", lw=1)
    for i in range(len(per_site)):
        ax.text(xmax * 1.02, yy[i], f"{smr[i]:.2f} ({lo[i]:.2f}–{hi[i]:.2f})", va="center", fontsize=9)
    if show_pool:
        ax.text(xmax * 1.02, yy[-1], f"{pooled['smr']:.2f} ({pooled['lo']:.2f}–{pooled['hi']:.2f})",
                va="center", fontsize=9, fontweight="bold")
    ax.set_yticks(yy); ax.set_yticklabels(ylabels)
    ax.set_xlim(min(lo) - 0.03, xmax * 1.20)
    ax.set_xlabel("Standardized Mortality Ratio (O/E), 95% CI")
    ax.set_title("Risk-Standardized 30-Day Mortality by Site")
    fig.tight_layout()
    fig.savefig(out_dir / "pooled_descriptive_smr_forest.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "pooled_descriptive_smr_forest.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote pooled_descriptive_smr_forest.{png,pdf}")
    return {"forest": True}


def build_smr_calibration_figure(sites, labels, out_dir: Path) -> dict:
    """Anonymized (Site 1..N) SMR transfer-calibration figure (obs vs exp by decile)."""
    import matplotlib.pyplot as plt
    cal = collect_smr_calibration(sites)
    if not cal:
        print("  (no SMR calibration data; skipping)")
        return {}
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 6))
    mx = 1.0
    for i, (sid, df) in enumerate(cal.items()):
        obs = df["observed"] / df["n"] * 100
        exp = df["expected"] / df["n"] * 100
        ax.plot(exp, obs, "o-", alpha=0.85, label=f"Site {i+1}", markersize=5)
        mx = max(mx, float(exp.max()), float(obs.max()))
    lim = min(100.0, mx * 1.08)
    ax.plot([0, lim], [0, lim], ls="--", color="grey", lw=1)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Expected mortality (%)"); ax.set_ylabel("Observed mortality (%)")
    ax.set_title("SMR Transfer Calibration by Risk Decile")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "pooled_descriptive_smr_calibration.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "pooled_descriptive_smr_calibration.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote pooled_descriptive_smr_calibration.{png,pdf}")
    return {"calibration": True}


# ── Pipeline entrypoint ────────────────────────────────────────────────────

# ── Consolidated pooled-numbers reference (writing companion) ───────────────

def _num_or_dash(v, dp=2, pct=False, dash="—"):
    """Format a scalar for the reference table; '—' for missing/non-finite."""
    if v is None:
        return dash
    try:
        f = float(v)
    except (TypeError, ValueError):
        return dash
    if not np.isfinite(f):
        return dash
    return f"{f:.{dp}f}" + ("%" if pct else "")


def _range_or_dash(lo, hi, dp=2, pct=False):
    """Format a 'lo–hi' range (en dash); '—' if either endpoint is missing."""
    a, b = _num_or_dash(lo, dp, pct), _num_or_dash(hi, dp, pct)
    if a == "—" or b == "—":
        return "—"
    return a if a == b else f"{a}–{b}"


def build_manuscript_numbers_reference(sites, labels, pooled_dir: Path,
                                       out_dir: Path) -> dict:
    """One consolidated 'manuscript numbers' sheet spanning descriptive epi + SMR +
    causal: for every pooled estimate, its 95% CI, the between-site range (min–max
    of per-site point estimates), I², a homogeneity p, and N. A WRITING COMPANION
    (not a manuscript display item) so all the numbers a writer needs live in one
    place. Descriptive + SMR rows are computed live from per-site CSVs; causal rows
    are read from the pooled CSVs the rest of main() has already written into
    `pooled_dir`. Robust to k==1 and to any missing causal input."""
    out_dir.mkdir(parents=True, exist_ok=True)
    COLS = ["Section", "Metric", "k", "Pooled", "95% CI", "Between-site range",
            "I² (%)", "Homogeneity p", "N", "Scale / estimand"]
    rows: list[dict] = []

    def add(section, metric, k, pooled, ci, rng, i2, homo_p, n, scale):
        rows.append(dict(zip(COLS, [section, metric, k, pooled, ci, rng, i2,
                                    homo_p, n, scale])))

    # ── Descriptive epidemiology (incidence, dose, bands, net UF, blood flow,
    #    time-to-CRRT, duration, modality, crude mortality) ──
    for r in compute_epi_heterogeneity(sites):
        prop = r.get("scale") == "proportion (%)"
        dp = 1 if prop else 2
        # Proportions report the CRUDE additive pool (Σnum/Σden) as the headline, not
        # the random-effects pool (which a 0%/100% boundary site can dominate). I²
        # is still shown from the RE fit as a heterogeneity flag.
        pooled_val = r.get("pooled_crude") if prop else r.get("pooled")
        ci_lo = r.get("pooled_crude_lo") if prop else r.get("pooled_lo")
        ci_hi = r.get("pooled_crude_hi") if prop else r.get("pooled_hi")
        scale_lbl = (r.get("scale", "") + ", crude pool (Σnum/Σden)") if prop else r.get("scale", "")
        add(r.get("group", ""), r.get("metric", ""), r.get("k_sites"),
            _num_or_dash(pooled_val, dp, pct=prop),
            _range_or_dash(ci_lo, ci_hi, dp, pct=prop),
            _range_or_dash(r.get("site_min"), r.get("site_max"), dp, pct=prop),
            _num_or_dash(r.get("I2_pct"), 1),
            _num_or_dash(r.get("between_site_p"), 4),
            int(r["n_total"]) if pd.notna(r.get("n_total")) else "",
            scale_lbl)

    # ── Risk-standardized mortality ──
    per_site_smr, pooled_smr = collect_smr(sites, labels)
    if pooled_smr:
        smrs = [p["smr"] for p in per_site_smr
                if pd.notna(p.get("smr")) and np.isfinite(p["smr"])]
        add("Risk-standardized mortality", "SMR (observed/expected)",
            pooled_smr.get("k"), _num_or_dash(pooled_smr.get("smr"), 3),
            _range_or_dash(pooled_smr.get("lo"), pooled_smr.get("hi"), 3),
            _range_or_dash(min(smrs), max(smrs), 3) if smrs else "—", "—",
            _num_or_dash(pooled_smr.get("het_p"), 4), pooled_smr.get("n"),
            "SMR = ΣO/ΣE, Byar 95% CI; homogeneity = Breslow-Day")

    # ── Causal: point-treatment HRs (high vs low dose) ──
    prim = pooled_dir / "pooled_primary_analyses.csv"
    if prim.exists():
        dfp = pd.read_csv(prim)
        ps_path = pooled_dir / "pooled_primary_per_site.csv"
        dfps = pd.read_csv(ps_path) if ps_path.exists() else None
        for _, r in dfp.iterrows():
            rng = "—"
            if dfps is not None and "analysis" in dfps.columns:
                hrs = pd.to_numeric(
                    dfps.loc[dfps["analysis"] == r["analysis"], "hr"],
                    errors="coerce").dropna()
                if len(hrs):
                    rng = _range_or_dash(hrs.min(), hrs.max(), 2)
            add("Causal — dose effect", str(r["analysis"]), r.get("n_sites"),
                _num_or_dash(r.get("re_hr"), 2),
                _range_or_dash(r.get("re_lo"), r.get("re_hi"), 2), rng,
                _num_or_dash(r.get("I2_pct"), 1), _num_or_dash(r.get("Q_p"), 4),
                "", "HR, random-effects (≥ 30 vs < 30 mL/kg/hr)")

    # ── Causal: continuous dose-response (per 10 mL/kg/hr) ──
    drp = pooled_dir / "pooled_dose_response_linear.csv"
    if drp.exists():
        dfd = pd.read_csv(drp)
        m = dfd[dfd["site"].astype(str) == "POOLED_RE_10"]
        if len(m):
            rr = m.iloc[0]
            add("Causal — dose effect", "Continuous dose, HR per 10 mL/kg/hr", None,
                _num_or_dash(rr.get("hr"), 2),
                _range_or_dash(rr.get("hr_lo"), rr.get("hr_hi"), 2), "—",
                _num_or_dash(rr.get("I2_pct"), 1), _num_or_dash(rr.get("Q_p"), 4),
                "", "HR per 10 mL/kg/hr, random-effects")

    # ── Causal: subgroup interactions (summary across pre-specified subgroups) ──
    sgp = pooled_dir / "pooled_subgroup_interactions.csv"
    if sgp.exists():
        dfs = pd.read_csv(sgp)
        if "re_p_int" in dfs.columns and dfs["re_p_int"].notna().any():
            pv = pd.to_numeric(dfs["re_p_int"], errors="coerce").dropna()
            qv = (pd.to_numeric(dfs["BH_q_re_p_int"], errors="coerce").dropna()
                  if "BH_q_re_p_int" in dfs.columns else pd.Series(dtype=float))
            add("Causal — subgroups",
                f"Interaction p-value ({len(pv)} pre-specified subgroups)", len(pv),
                _range_or_dash(pv.min(), pv.max(), 3), "—", "—", "—",
                _num_or_dash(qv.min(), 3) if len(qv) else "—", "",
                "range of pooled interaction p; Homogeneity-p col = min BH q")

    ref_df = pd.DataFrame(rows, columns=COLS)
    csv_path = out_dir / "pooled_manuscript_numbers.csv"
    html_path = out_dir / "pooled_manuscript_numbers.html"
    ref_df.to_csv(csv_path, index=False)

    caption = (
        f"Consolidated pooled estimates with between-site variation across "
        f"{len(sites)} site(s). Writing companion, not a manuscript display item. "
        "Descriptive rows are pooled by random-effects meta-analysis (proportions "
        "on the logit scale, continuous metrics on the mean); the between-site "
        "range is the min-max of per-site point estimates, and I-squared with the "
        "homogeneity p quantify between-site heterogeneity. SMR is pooled as "
        "sum(O)/sum(E) with a Byar 95% CI and a Breslow-Day homogeneity test. "
        "Causal rows are random-effects hazard ratios for high- versus low-dose "
        "CRRT. Note that I-squared is inflated for high-N metrics (labs, "
        "incidence): with tens of thousands of observations the sampling error is "
        "tiny, so even a clinically trivial between-site difference registers as "
        "high heterogeneity. For those rows the between-site range is the "
        "interpretable measure of variation, not I-squared. All pooled values are "
        "provisional and will change as sites re-run and additional sites join."
    )
    table_rows = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(v), quote=False)}</td>"
                         for v in row) + "</tr>"
        for row in ref_df.itertuples(index=False, name=None)
    )
    flat_html = (
        '<table class="results-table" border="0"><thead><tr>'
        + "".join(f"<th>{html.escape(c, quote=False)}</th>" for c in ref_df.columns)
        + "</tr></thead><tbody>" + table_rows + "</tbody></table>"
    )
    html_doc = _MANUSCRIPT_TABLE_HTML_TEMPLATE.format(
        title="Consolidated Pooled Numbers — CRRT Manuscript (writing companion)",
        caption=caption,
        table=flat_html,
    )
    html_path.write_text(html_doc, encoding="utf-8")
    print(f"  wrote {csv_path.name} ({len(ref_df)} metrics)")
    print(f"  wrote {html_path.name}")
    return {"csv": csv_path, "html": html_path}


def main():
    sites = discover_sites()
    if not sites:
        print("No sites found in", ROOT)
        return
    print(f"Found {len(sites)} sites: {[s.name for s in sites]}")

    # Causal pooling is restricted; descriptive pooling and the SMR are not.
    causal_sites, causal_excluded = split_causal_eligible(sites, SITE_LABELS)
    print(f"\nCausal-pooling eligibility (high-dose arm >= {MIN_HIGH_DOSE_ARM}):")
    if causal_excluded:
        for e in causal_excluded:
            print(f"  EXCLUDED {e['label']}: high-dose arm n={e['n_high_dose']}"
                  f" of {e['n_total']}")
        _kept = sum(int(e["n_total"] or 0) for e in causal_excluded)
        print(f"  {len(causal_sites)} of {len(sites)} sites eligible; "
              f"{len(causal_excluded)} excluded ({_kept} patients)")
    else:
        print("  all sites eligible")
    print("  (descriptive pooling, SMR and per-site diagnostics use ALL sites)")

    # Refresh the gate verdicts before pooling anything, so site_validation.json
    # and the review tracker always describe the site set these artifacts were
    # built from. (Documentation only — pooling does not yet filter on it.)
    print("\nRefreshing site validation (gates -> site_validation.json + tracker)...")
    site_validation.run(root=ROOT, verbose=False)

    # Pooled CSVs (figures/tables read from these).
    export_pooled_csvs(causal_sites, SITE_LABELS, ROOT / "pooled", all_sites=sites)

    # Descriptive epidemiology pooled CSVs (from scripts 07/08).
    print("\nExporting descriptive epidemiology pooled CSVs...")
    export_descriptive_pooled_csvs(sites, SITE_LABELS, ROOT / "pooled")
    export_pooled_low_dose_baseline(sites, SITE_LABELS, ROOT / "pooled")
    export_epi_heterogeneity_csv(sites, SITE_LABELS, ROOT / "pooled")
    export_pooled_dose_ibw_csv(sites, SITE_LABELS, ROOT / "pooled")

    # Manuscript figures (1 → 5 in numerical order).
    print("\nBuilding manuscript Figure 1 (cohort exclusion flow diagram)...")
    build_cohort_flow_figure(sites, ROOT / "figures", causal_excluded=causal_excluded)

    print("\nBuilding manuscript Figure 2 (IPTW forest plot)...")
    build_iptw_forest_figure(causal_sites, ROOT / "figures", ROOT / "pooled")

    print("\nBuilding manuscript Figure 3 (PSM CIF death)...")
    build_psm_cif_figure(causal_sites, ROOT / "figures")

    print("\nBuilding manuscript Figure 4 (dose-response)...")
    build_dose_response_figure(causal_sites, ROOT / "figures")

    print("\nBuilding manuscript Figure 5 (subgroup forest)...")
    build_subgroup_forest_figure(causal_sites, ROOT / "figures", ROOT / "pooled")

    # Manuscript tables (T1: pooled baseline by dose; T2: pooled HR table).
    print("\nBuilding manuscript Table 1 (full-cohort baseline by dose band)...")
    build_baseline_table(sites, SITE_LABELS, ROOT / "tables")

    print("\nBuilding manuscript Table 2 (unadjusted balance, causal cohort)...")
    build_unadjusted_balance_table(sites, SITE_LABELS, ROOT / "tables")

    print("\nBuilding manuscript Table 3 (pooled HRs)...")
    build_causal_hrs_table(causal_sites, SITE_LABELS, ROOT / "tables", ROOT / "pooled")

    print("\nBuilding descriptive epidemiology table...")
    build_descriptive_table(sites, SITE_LABELS, ROOT / "tables")

    # Supplemental causal-model diagnostics (per-site, anonymized; inputs already on disk)
    build_supplement_causal_diagnostics(sites, SITE_LABELS, ROOT / "tables", ROOT / "pooled")

    print("\nBuilding descriptive epidemiology figures...")
    build_descriptive_figures(sites, SITE_LABELS, ROOT / "figures")
    build_pooled_dose_distribution(sites, SITE_LABELS, ROOT / "figures")
    build_pooled_dose_ibw(sites, SITE_LABELS, ROOT / "figures")
    # Named "_named" variants of the site-identifying descriptive figures, consumed
    # only by the non-anonymized dashboard (07). The manuscript + anonymized
    # dashboard keep using the canonical anonymized stems above.
    build_descriptive_figures(sites, SITE_LABELS, ROOT / "figures", anonymize=False)
    build_pooled_dose_distribution(sites, SITE_LABELS, ROOT / "figures", anonymize=False)
    build_pooled_uf_mortality(sites, SITE_LABELS, ROOT / "figures")
    build_pooled_state_figures(sites, SITE_LABELS, ROOT / "figures")
    build_pooled_trajectories(sites, SITE_LABELS, ROOT / "figures")

    print("\nBuilding SMR artifacts...")
    export_pooled_smr_csvs(sites, SITE_LABELS, ROOT / "pooled")
    build_smr_forest_figure(sites, SITE_LABELS, ROOT / "figures")
    build_smr_calibration_figure(sites, SITE_LABELS, ROOT / "figures")
    build_smr_per_site_table(sites, SITE_LABELS, ROOT / "tables")
    build_iptw_cif_figure(causal_sites, ROOT / "figures")

    print("\nBuilding remaining supplement tables...")
    build_crrt_settings_by_mode_table(sites, SITE_LABELS, ROOT / "tables")   # S-Table 1
    build_severity_thresholds_table(sites, SITE_LABELS, ROOT / "tables")     # S-Table 4
    build_missingness_table(sites, SITE_LABELS, ROOT / "tables")
    build_modality_table(sites, SITE_LABELS, ROOT / "tables")
    build_ph_diagnostics_table(sites, SITE_LABELS, ROOT / "tables", ROOT / "pooled")

    # Consolidated pooled-numbers reference (writing companion; reads the pooled
    # CSVs written above, so it runs last).
    print("\nBuilding consolidated pooled-numbers reference...")
    build_manuscript_numbers_reference(sites, SITE_LABELS, ROOT / "pooled",
                                       ROOT / "tables")


# %% Pipeline entrypoint — full manuscript-artifact rebuild
if __name__ == "__main__":
    main()


# %% Render Figure 1 only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_cohort_flow_figure(_sites, ROOT / "figures")


# %% Render Figure 2 only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_iptw_forest_figure(_sites, ROOT / "figures", ROOT / "pooled")


# %% Render Figure 3 only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_psm_cif_figure(_sites, ROOT / "figures")


# %% Render Figure 4 only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_dose_response_figure(_sites, ROOT / "figures")


# %% Render Figure 5 only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_subgroup_forest_figure(_sites, ROOT / "figures", ROOT / "pooled")


# %% Render Table 1 only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_baseline_table(_sites, SITE_LABELS, ROOT / "tables")


# %% Render Table 2 only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_unadjusted_balance_table(_sites, SITE_LABELS, ROOT / "tables")


# %% Render Table 3 only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_causal_hrs_table(_sites, SITE_LABELS, ROOT / "tables", ROOT / "pooled")


# %% Export descriptive epidemiology pooled CSVs only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    export_descriptive_pooled_csvs(_sites, SITE_LABELS, ROOT / "pooled")
    export_pooled_low_dose_baseline(_sites, SITE_LABELS, ROOT / "pooled")
    export_epi_heterogeneity_csv(_sites, SITE_LABELS, ROOT / "pooled")
    export_pooled_dose_ibw_csv(_sites, SITE_LABELS, ROOT / "pooled")


# %% Render descriptive epidemiology table only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_descriptive_table(_sites, SITE_LABELS, ROOT / "tables")


# %% Render descriptive epidemiology figures only (interactive cell)
if "ipykernel" in sys.modules:
    _sites = discover_sites()
    build_descriptive_figures(_sites, SITE_LABELS, ROOT / "figures")
    build_pooled_dose_distribution(_sites, SITE_LABELS, ROOT / "figures")
    build_pooled_dose_ibw(_sites, SITE_LABELS, ROOT / "figures")
    build_pooled_uf_mortality(_sites, SITE_LABELS, ROOT / "figures")
    build_pooled_state_figures(_sites, SITE_LABELS, ROOT / "figures")
    build_pooled_trajectories(_sites, SITE_LABELS, ROOT / "figures")
