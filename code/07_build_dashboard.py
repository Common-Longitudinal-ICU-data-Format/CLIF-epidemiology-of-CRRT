"""
07_build_dashboard.py
Combines per-site results (descriptive CRRT epidemiology + PSM/IPTW / dose-
response causal inference) into a single tabbed HTML dashboard. (The time-
varying MSM analysis was cut from the manuscript and archived, tag msm-v1.)

Shared infrastructure (config, IO, color/HTML rendering, meta-analysis,
combined-table compute, between-site heterogeneity) lives in `report_core.py`,
imported below — the manuscript-artifact builder (08_manuscript_artifacts.py)
imports the same module.

Input:  output/multi_site/{SITE_NAME}/final/ — one folder per site,
        mirroring the output/final/ structure (psm_iptw/, crrt_epi/, etc.)
Output: output/multi_site/crrt_dashboard.html
"""

import base64
import html
import io
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Presentation-ready defaults: larger fonts + bold titles/labels for
# readability when figures are projected or embedded in slides. Tight layout
# pad reduced globally via PAD_DEFAULT below.
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 13,
    "axes.titlesize": 16,
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "axes.labelweight": "bold",
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "legend.title_fontsize": 13,
    "figure.titlesize": 16,
})
PAD_DEFAULT = 0.5  # tight_layout pad — smaller = less whitespace
import numpy as np
import pandas as pd
from scipy import stats

import site_validation

# ── Config ──────────────────────────────────────────────────────────────────

from report_core import (
    ROOT,
    CONSUMED,
    SITE_LABELS,
    SITE_COLORS,
    _EXTRA_COLORS,
    BALANCE_TABLE_PATTERNS,
    get_color,
    discover_sites,
    _site_final_dir,
    fig_to_data_uri,
    img_to_data_uri,
    _find_csv,
    _find_png,
    _load_csv,
    _missing,
    _styled_table,
    _se_from_ci,
    meta_analyze,
    _collect_iptw,
    _collect_psm_fg,
    _parse_gtsummary_csv,
    _extract_n_from_header,
    _weighted_median,
    _parse_stat_cell,
    _decimal_places_for_label,
    _CCI_LABELS,
    _MODALITY_MISSING_LABELS,
    _uppercase_modality_rows,
    _render_gtsummary_table_html,
    _compute_combined_table_df,
    _combine_table_across_sites,
    _read_site_causal_n,
    _logit,
    _expit,
    _het_pool,
    compute_epi_heterogeneity,
    EPI_HET_COLS,
    compute_dose_ibw_pooled,
    collect_smr,
    collect_smr_calibration,
    load_smr_model,
)
OUTPUT = ROOT / "crrt_dashboard.html"

# Pooled figures embedded in the "Cross-Site Pooled Figures" section. These PNGs
# are BUILT BY 08_manuscript_artifacts.py — 07 only embeds them, so 08 must run
# first after any change to the site set (`_warn_if_pooled_figures_stale`
# enforces this by checking mtimes against the collected site data). Single
# source of truth for both the staleness check and the render loop.
_POOLED_FIGURE_SPECS = [
    ("pooled_descriptive_incidence",
     "CRRT incidence across consortium sites (overall, ventilated, vasopressor-exposed)"),
    ("pooled_descriptive_dose_distribution",
     "Pooled CRRT dose distribution (histogram bins summed across sites; per-site density outlines)"),
    ("pooled_descriptive_dose_ibw",
     "Pooled CRRT dose: actual vs ideal body weight (Devine), on the height-available paired subset"),
    ("pooled_descriptive_dosebands",
     "CRRT dose-band composition by site"),
    ("pooled_descriptive_traj_dose",
     "Pooled CRRT dose over 30 days (N-weighted median, IQR ribbon)"),
    ("pooled_descriptive_traj_uf_rate",
     "Pooled net ultrafiltration rate over 30 days (N-weighted median, IQR ribbon)"),
    ("pooled_descriptive_traj_uf_cumulative",
     "Pooled cumulative net ultrafiltration over 30 days (N-weighted median, IQR ribbon)"),
    ("pooled_descriptive_traj_nee",
     "Pooled vasopressor (norepinephrine-equivalent) over 30 days (N-weighted median, IQR ribbon)"),
    ("pooled_descriptive_traj_labs",
     "Pooled lab trajectories over 30 days (lactate, pH, bicarbonate, potassium, phosphate)"),
    ("pooled_descriptive_uf_mortality",
     "Pooled crude 30-day mortality vs first-72h net UF intensity (counts summed per fixed bin)"),
    ("pooled_descriptive_imv_state",
     "Pooled invasive mechanical ventilation state over 30 days (per-hour counts summed)"),
    ("pooled_descriptive_crrt_state",
     "Pooled CRRT state over 30 days (per-hour counts summed)"),
]




# Per-site unadjusted balance table (causal cohort, by dose group) = manuscript
# Table 2. 05.R now emits "{SITE}_table2_unadjusted_balance.csv"; the legacy
# "{SITE}_Table1_unadjusted.csv" is kept as a fallback (tried second) until
# every site re-runs 05. Resolvers below try these patterns in order and
# prefer the new name, so a stale legacy orphan never shadows the new file.
# Remove the legacy entry once all sites have re-run 05.

# ── Expected Output Files (for completeness checklist) ────────────────────
# Each entry: (subdir, pattern, label, extension, script_hint)

EXPECTED_FILES = {
    "Descriptive Epidemiology (00 / 02 / 03)": [
        # Cohort + Table 1
        ("crrt_epi", "strobe_counts", "STROBE counts", "csv", "00"),
        ("crrt_epi", "table1_crrt.csv", "Table 1 (baseline by dose band)", "csv", "02"),
        # Incidence + practice/quality (feed pooled incidence figure + descriptive table)
        ("crrt_epi", "crrt_incidence", "CRRT incidence (context + ICU subtype)", "csv", "03"),
        ("crrt_epi", "crrt_practice_quality", "CRRT practice / quality metrics", "csv", "03"),
        # Distribution + net-UF intensity summaries
        ("crrt_epi/graphs", "dose_distribution.csv", "Dose distribution", "csv", "03"),
        ("crrt_epi/graphs", "net_uf_intensity_summary", "Net UF intensity summary (first 72h)", "csv", "03"),
        # 30-day trajectory + state CSVs
        ("crrt_epi/graphs", "crrt_dose_hourly", "Dose over time (30d)", "csv", "03"),
        ("crrt_epi/graphs", "lab_distributions_over_crrt", "Lab trajectories (30d)", "csv", "03"),
        ("crrt_epi/graphs", "nee_over_crrt", "Vasopressor (NEE) trajectory (30d)", "csv", "03"),
        ("crrt_epi/graphs", "net_uf_rate_hourly", "Net UF rate trajectory (30d)", "csv", "03"),
        ("crrt_epi/graphs", "net_uf_cumulative", "Net UF cumulative (30d)", "csv", "03"),
        ("crrt_epi/graphs", "uf_mortality", "Net UF intensity vs mortality", "csv", "03"),
        ("crrt_epi/graphs", "imv_state_over_crrt", "IMV state (30d)", "csv", "03"),
        ("crrt_epi/graphs", "crrt_state_over_crrt", "CRRT state (30d)", "csv", "03"),
        # Internal QC / data-completeness diagnostics (00 / 04) -> diagnostics/
        ("diagnostics", "missingness_summary", "Missingness summary", "csv", "04"),
        ("diagnostics", "crrt_initial_field_completeness", "CRRT field completeness at start (first-3h)", "csv", "00"),
        ("diagnostics", "crrt_column_missingness", "CRRT column missingness", "csv", "00"),
        ("diagnostics", "crrt_settings_distribution", "CRRT settings distribution", "csv", "00"),
        ("diagnostics", "crrt_settings_summary", "CRRT settings summary", "csv", "00"),
        # Representative PNGs
        ("crrt_epi/graphs", "crrt_incidence_by_context", "Incidence by context", "png", "03"),
        ("crrt_epi/graphs", "dose_over_time", "Dose over time", "png", "03"),
        ("crrt_epi/graphs", "net_uf_rate_over_time", "Net UF rate over time", "png", "03"),
        ("crrt_epi/graphs", "imv_state", "IMV state (30d)", "png", "03"),
        ("crrt_epi/graphs", "crrt_state", "CRRT state (30d)", "png", "03"),
    ],
    "Risk-Standardized Mortality (Script 03b)": [
        ("crrt_epi", "smr.csv", "SMR (observed/expected, Byar CI)", "csv", "03b"),
        ("crrt_epi", "smr_calibration", "SMR transfer calibration (by decile)", "csv", "03b"),
    ],
    "Low-Dose Characterization (Script 06)": [
        ("low_dose", "low_dose_counts", "Dose-band counts (incl. very-low 10-15)", "csv", "06"),
        ("low_dose", "low_dose_long", "Very-low-vs-rest baseline (long)", "csv", "06"),
        ("low_dose", "low_dose_table", "Very-low-vs-rest table", "csv", "06"),
    ],
    "PSM / IPTW (Script 05)": [
        ("psm_iptw", "sample_characteristics", "Sample characteristics", "csv", "05"),
        ("psm_iptw", "IPTW_pooled_results", "IPTW pooled results", "csv", "05"),
        ("psm_iptw", "fg_psm_pooled_results", "Fine-Gray PSM pooled results", "csv", "05"),
        ("psm_iptw", "PSM_CIF_data", "PSM CIF data", "csv", "05"),
        ("psm_iptw", "IPTW_CIF_data", "IPTW CIF data", "csv", "05"),
        ("psm_iptw", "evalue_sensitivity", "E-value sensitivity", "csv", "05"),
        ("psm_iptw", "SL_IPTW_ESS", "SL IPTW ESS", "csv", "05"),
        ("psm_iptw", "SL_PS_extremes", "PS extremes", "csv", "05"),
        ("psm_iptw", BALANCE_TABLE_PATTERNS, "Table 2 (unadjusted balance, causal cohort)", "csv", "05"),
        ("psm_iptw", "TableS1_matched", "Table S1 matched", "csv", "05"),
        ("psm_iptw", "TableS2_IPTW", "Table S2 IPTW", "csv", "05"),
        ("psm_iptw", "subgroup_analysis_results", "Subgroup analysis", "csv", "05"),
        ("psm_iptw", "ModelComparison", "Model comparison", "csv", "05"),
        ("psm_iptw", "causal_consort_diagram", "Causal cohort flow", "png", "04"),
        ("psm_iptw", "psm_loveplot", "PSM love plot", "png", "05"),
        ("psm_iptw", "SL_LovePlot", "SL love plot", "png", "05"),
        ("psm_iptw", "SL_PS_overlap", "SL PS overlap", "png", "05"),
        ("psm_iptw", "IPTW_CIF_Death", "IPTW CIF Death", "png", "05"),
        ("psm_iptw", "IPTW_CIF_Discharge", "IPTW CIF Discharge", "png", "05"),
        ("psm_iptw", "subgroup_forest_plot", "Subgroup forest plot", "png", "05"),
    ],
    "Dose-Response (Script 05b)": [
        ("psm_iptw", "dose_decile_mortality", "Dose-decile mortality (Fig 4a)", "csv", "05b"),
        ("psm_iptw", "dose_response_rcs.csv", "RCS dose-response (Fig 4b)", "csv", "05b"),
        ("psm_iptw", "dose_response_rcs_anova", "RCS ANOVA (nonlinearity LRT)", "csv", "05b"),
        ("psm_iptw", "gps_dose_response.csv", "GPS dose-response (Fig 4c)", "csv", "05b"),
        ("psm_iptw", "dose_response_combined", "Combined dose-response", "png", "05b"),
        ("psm_iptw", "dose_response_rcs_plot", "RCS dose-response plot", "png", "05b"),
        ("psm_iptw", "gps_dose_response_plot", "GPS dose-response plot", "png", "05b"),
        ("psm_iptw", "gps_loveplot", "GPS love plot", "png", "05b"),
    ],
}


# Recognized per-site outputs that legitimately do NOT feed the dashboard, so
# the completeness flag must not mislabel them as stale/orphan. Matched as
# case-insensitive substrings of the (site-prefixed) filename. Two families:
#   - Manuscript (08) federated-exchange inputs: script 08 builds the POOLED
#     figures from these per-site CSVs; the dashboard embeds 08's finished PNGs,
#     never these raw inputs.
#   - Deep causal diagnostics (05 / 05b): written for provenance/QC, not shown.
# Files under diagnostics/ are handled structurally (a QC bucket) and need not
# be listed here. Keep this in sync when 03 / 05 / 05b add a non-dashboard CSV.
NON_DASHBOARD_OUTPUTS = (
    # manuscript-08 pooled-figure inputs (script 03)
    "crrt_census", "crrt_course_time_to_event",
    "dose_distribution_ibw", "dose_distribution_actual_paired",
    # deep causal diagnostics (05 / 05b) — "cox_diagnostics" also covers the
    # dose_response_* and iptw_* _per_covariate variants
    "crrt_bin_distribution", "cox_diagnostics", "fg_psm_dr_results",
    "gps_balance", "gps_weights_summary", "iptw_causespecificcox_fullresults",
    "psm_balance_summary", "sl_iptw_balance", "winsorization_audit",
)


# ── Helpers ─────────────────────────────────────────────────────────────────

















def _build_combined_cif(
    site_cifs: list[tuple[str, str, pd.DataFrame]],
    time_col: str, cif_col: str, group_col: str,
    title: str, ylabel: str,
    max_time: float = 30.0,
    n_grid: int = 200,
) -> str:
    """Build a combined CIF plot across sites with mean line + min/max ribbon.

    site_cifs: list of (site_label, site_id, dataframe) tuples
    Returns HTML string with embedded plot.
    """
    if not site_cifs:
        return _missing(f"{title}: CIF data not available.")

    # Common time grid
    tgrid = np.linspace(0, max_time, n_grid)

    # Collect all groups across sites
    all_groups = []
    for _, _, df in site_cifs:
        all_groups.extend(df[group_col].unique())
    groups = sorted(set(all_groups), key=str)

    if not groups:
        return _missing(f"{title}: no groups found in CIF data.")

    # For each group, interpolate each site's CIF onto tgrid, then compute mean + range
    fig, ax = plt.subplots(figsize=(10, 5))
    group_colors = {}
    color_cycle = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                   "#56B4E9", "#F0E442", "#882255"]

    for gi, grp in enumerate(groups):
        color = color_cycle[gi % len(color_cycle)]
        group_colors[grp] = color

        # Interpolate CIF for each site
        site_curves = []
        for _site_lbl, _site_id, df in site_cifs:
            gdf = df[df[group_col] == grp].sort_values(time_col)
            if gdf.empty:
                continue
            # Step-function interpolation (CIF is a step function)
            interp = np.interp(tgrid, gdf[time_col].values, gdf[cif_col].values,
                               left=0.0)
            site_curves.append(interp)

        if not site_curves:
            continue

        curves = np.array(site_curves)
        mean_curve = curves.mean(axis=0)
        min_curve = curves.min(axis=0)
        max_curve = curves.max(axis=0)

        # Shorten group label for legend
        short_grp = grp.split("(")[0].strip() if "(" in grp else grp
        ax.plot(tgrid, mean_curve, color=color, linewidth=2, label=short_grp)
        ax.fill_between(tgrid, min_curve, max_curve, color=color, alpha=0.15)

        # If only 1 site, also plot per-site bootstrap CI as lighter ribbon
        if len(site_curves) == 1:
            _, _, df = site_cifs[0]
            gdf = df[df[group_col] == grp].sort_values(time_col)
            # Try standard naming patterns
            for lo_try, hi_try in [
                (f"{cif_col}_lower", f"{cif_col}_upper"),
                ("cif_lower", "cif_upper"),
                (f"ci_{cif_col.split('_')[-1]}_lower", f"ci_{cif_col.split('_')[-1]}_upper"),
            ]:
                if lo_try in gdf.columns and hi_try in gdf.columns:
                    lo_interp = np.interp(tgrid, gdf[time_col].values,
                                          gdf[lo_try].values, left=0.0)
                    hi_interp = np.interp(tgrid, gdf[time_col].values,
                                          gdf[hi_try].values, left=0.0)
                    ax.fill_between(tgrid, lo_interp, hi_interp, color=color, alpha=0.08)
                    break

    n_sites = len(site_cifs)
    ribbon_label = "min\u2013max across sites" if n_sites > 1 else "bootstrap 95% CI"
    ax.set_xlabel("Time from CRRT Initiation (Days)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\n({n_sites} site{'s' if n_sites != 1 else ''}, "
                 f"ribbon = {ribbon_label})")
    ax.legend(loc="upper left")
    ax.set_xlim(0, max_time)
    ax.set_ylim(0, None)
    ax.grid(alpha=0.3)
    fig.tight_layout(pad=PAD_DEFAULT)

    return (
        f'<div class="figure-block"><h3>{title}</h3>'
        f'<img src="{fig_to_data_uri(fig)}" alt="{title}"></div>'
    )


def build_combined_cif_psm(sites, labels) -> str:
    """Combined CIF curves from PSM analysis across sites."""
    html_parts = []
    for outcome_code, outcome_name in [("2", "Death"), ("1", "Discharge")]:
        site_cifs = []
        for sd in sites:
            df = _load_csv(sd, "psm_iptw", "PSM_CIF_data")
            if df is None:
                continue
            sub = df[df["outcome"].astype(str) == outcome_code]
            if sub.empty:
                continue
            site_cifs.append((labels.get(sd.name, sd.name), sd.name, sub))
        html_parts.append(_build_combined_cif(
            site_cifs, time_col="time", cif_col="cif", group_col="group",
            title=f"PSM Cumulative Incidence: {outcome_name}",
            ylabel=f"Cumulative Incidence ({outcome_name})",
        ))
    return "\n".join(html_parts) if html_parts else ""


def build_combined_cif_iptw(sites, labels) -> str:
    """Combined CIF curves from IPTW analysis across sites."""
    html_parts = []
    for cif_col, outcome_name in [("cif_death", "Death"), ("cif_disch", "Discharge")]:
        site_cifs = []
        for sd in sites:
            df = _load_csv(sd, "psm_iptw", "IPTW_CIF_data")
            if df is None:
                continue
            site_cifs.append((labels.get(sd.name, sd.name), sd.name, df))
        html_parts.append(_build_combined_cif(
            site_cifs, time_col="time", cif_col=cif_col, group_col="trt_group",
            title=f"IPTW Standardized CIF: {outcome_name}",
            ylabel=f"Cumulative Incidence ({outcome_name})",
        ))
    return "\n".join(html_parts) if html_parts else ""


def _fmt_hr(row) -> str:
    return f"{row['HR']:.2f} ({row['HR_lower']:.2f}\u2013{row['HR_upper']:.2f})"


def _fmt_p(p) -> str:
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


# ── Section Builders ────────────────────────────────────────────────────────

def build_sample_characteristics(sites, labels) -> str:
    """Table of sample characteristics across sites."""
    rows = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "sample_characteristics")
        if df is None:
            continue
        r = df.iloc[0]
        rows.append({
            "site_id": sd.name,
            "Site": labels.get(sd.name, sd.name),
            "N": int(r.get("total_n", 0)),
            "Age (mean)": f"{r.get('age_mean', 0):.1f}",
            "Male %": f"{r.get('sex_male_pct', r.get('sex_pct_male', 0)):.1f}",
            "SOFA (median)": f"{r.get('sofa_median', 0):.1f}",
            "Dose (median)": f"{r.get('crrt_dose_median', 0):.1f}",
            "Mortality %": f"{r.get('outcome_death_pct', 0):.1f}",
        })
    if not rows:
        return _missing("Sample characteristics not available.")
    tbl = pd.DataFrame(rows)
    return (
        '<div class="section"><h2>Sample Characteristics</h2>'
        '<p class="fig-caption">Key demographics and clinical features at CRRT initiation for each site.</p>'
        + _styled_table(tbl)
        + '</div>'
    )


def build_forest_plot_iptw(sites, labels) -> str:
    """Forest plot of pooled IPTW Cox HRs across sites."""
    rows = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "IPTW_pooled_results")
        if df is None:
            continue
        for _, r in df.iterrows():
            rows.append({
                "site": labels.get(sd.name, sd.name),
                "site_id": sd.name,
                "model": r["model"],
                "HR": r["HR"],
                "HR_lower": r["HR_lower"],
                "HR_upper": r["HR_upper"],
                "p_value": r["p_value"],
            })
    if not rows:
        return _missing("IPTW pooled results not available.")

    all_df = pd.DataFrame(rows)
    html_parts = []

    for outcome in ["Death", "Discharge"]:
        sub = all_df[all_df["model"].str.contains(outcome, case=False)]
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(8, max(2, len(sub) * 0.6 + 1)))
        y_pos = list(range(len(sub)))
        for i, (_, r) in enumerate(sub.iterrows()):
            color = get_color(r["site_id"])
            ax.errorbarx = ax.errorbar(
                r["HR"], i,
                xerr=[[r["HR"] - r["HR_lower"]], [r["HR_upper"] - r["HR"]]],
                fmt="o", color=color, markersize=8, capsize=4, linewidth=2,
                label=r["site"]
            )
        ax.axvline(x=1, color="grey", linestyle="--", linewidth=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(sub["site"].tolist())
        ax.set_xlabel("Hazard Ratio (95% CI)")
        ax.set_title(f"IPTW Cox — {outcome} (High vs Low Dose)")
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
        fig.tight_layout(pad=PAD_DEFAULT)
        html_parts.append(
            f'<div class="figure-block"><h3>IPTW Forest Plot — {outcome}</h3>'
            f'<img src="{fig_to_data_uri(fig)}" alt="Forest plot {outcome}"></div>'
        )

    # Results table
    tbl_rows = []
    for _, r in all_df.iterrows():
        tbl_rows.append({
            "site_id": r["site_id"],
            "Site": r["site"],
            "Outcome": ("Death" if "Death" in r["model"] else
                        "Discharge" if "Discharge" in r["model"] else
                        r["model"]),
            "HR (95% CI)": _fmt_hr(r),
            "p": _fmt_p(r["p_value"]),
        })
    tbl = pd.DataFrame(tbl_rows)
    html_parts.append(_styled_table(tbl))

    return ('<div class="section"><h2>Point-Treatment IPTW Cox Results</h2>'
            '<p class="fig-caption">Cause-specific Cox hazard ratios for high (&ge;30 mL/kg/hr) vs. low (&lt;30) initial CRRT dose, estimated via Super Learner IPTW.</p>'
            + "\n".join(html_parts) + '</div>')


def build_forest_plot_psm_fg(sites, labels) -> str:
    """Forest plot of PSM Fine-Gray SHRs across sites."""
    rows = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "fg_psm_pooled_results")
        if df is None:
            continue
        for _, r in df.iterrows():
            rows.append({
                "site": labels.get(sd.name, sd.name),
                "site_id": sd.name,
                "model": r["model"],
                "HR": r["HR"],
                "HR_lower": r["HR_lower"],
                "HR_upper": r["HR_upper"],
                "p_value": r["p_value"],
            })
    if not rows:
        return _missing("PSM Fine-Gray results not available.")

    all_df = pd.DataFrame(rows)
    html_parts = []

    for outcome in ["Death", "Discharge"]:
        sub = all_df[all_df["model"].str.contains(outcome, case=False)]
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(8, max(2, len(sub) * 0.6 + 1)))
        y_pos = list(range(len(sub)))
        for i, (_, r) in enumerate(sub.iterrows()):
            color = get_color(r["site_id"])
            ax.errorbar(
                r["HR"], i,
                xerr=[[r["HR"] - r["HR_lower"]], [r["HR_upper"] - r["HR"]]],
                fmt="o", color=color, markersize=8, capsize=4, linewidth=2,
                label=r["site"]
            )
        ax.axvline(x=1, color="grey", linestyle="--", linewidth=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(sub["site"].tolist())
        ax.set_xlabel("Subdistribution Hazard Ratio (95% CI)")
        ax.set_title(f"PSM Fine-Gray — {outcome} (High vs Low Dose)")
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
        fig.tight_layout(pad=PAD_DEFAULT)
        html_parts.append(
            f'<div class="figure-block"><h3>PSM Fine-Gray Forest Plot — {outcome}</h3>'
            f'<img src="{fig_to_data_uri(fig)}" alt="PSM FG forest {outcome}"></div>'
        )

    tbl_rows = []
    for _, r in all_df.iterrows():
        tbl_rows.append({
            "site_id": r["site_id"],
            "Site": r["site"],
            "Outcome": ("Death" if "Death" in r["model"] else
                        "Discharge" if "Discharge" in r["model"] else
                        r["model"]),
            "SHR (95% CI)": _fmt_hr(r),
            "p": _fmt_p(r["p_value"]),
        })
    tbl = pd.DataFrame(tbl_rows)
    html_parts.append(_styled_table(tbl))

    return ('<div class="section"><h2>Point-Treatment PSM Fine-Gray Results</h2>'
            '<p class="fig-caption">Subdistribution hazard ratios (Fine-Gray competing risks) for high (&ge;30 mL/kg/hr) vs. low (&lt;30) initial CRRT dose, estimated on the propensity-score-matched cohort.</p>'
            + "\n".join(html_parts) + '</div>')


def build_model_comparison(sites, labels) -> str:
    """Model comparison table (PSM FG vs IPTW Cox) across sites."""
    rows = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "ModelComparison")
        if df is None:
            continue
        # Drop the MICE-pooled duplicate rows; FMI ~0 makes them
        # statistically identical to the single-imputation rows and the
        # "(pooled)" label is easily confused with cross-site pooling.
        df = df[~df["model"].str.contains(r"\(pooled\)", regex=True, na=False)]
        for _, r in df.iterrows():
            rows.append({
                "site_id": sd.name,
                "Site": labels.get(sd.name, sd.name),
                "Model": r["model"],
                "HR Type": r.get("HR_type", ""),
                "HR (95% CI)": _fmt_hr(r),
                "p": _fmt_p(r["p_value"]),
            })
    if not rows:
        return _missing("Model comparison not available.")
    tbl = pd.DataFrame(rows)
    return (
        '<div class="section"><h2>Model Comparison (PSM vs IPTW)</h2>'
        '<p class="fig-caption">Concordance between propensity score matching (Fine-Gray SHR) and IPTW (cause-specific Cox HR) across sites.</p>'
        + _styled_table(tbl)
        + '</div>'
    )


def build_dose_response(sites, labels) -> str:
    """Overlay dose-response curves (RCS + GPS) across sites."""
    # Display window for dose-response overlays. Clinical CRRT doses span
    # ~15-45 mL/kg/hr; trials go to 35. Beyond 80 is extrapolation and at
    # one site reflects known data-quality outliers (Penn: max=264).
    DOSE_XLIM = (0, 80)
    html_parts = []

    for _analysis, pattern, title in [
        ("RCS", "dose_response_rcs", "Adjusted Dose-Response (Restricted Cubic Spline)"),
        ("GPS", "gps_dose_response", "Generalized Propensity Score Dose-Response"),
    ]:
        fig, ax = plt.subplots(figsize=(10, 5))
        any_data = False
        for sd in sites:
            df = _load_csv(sd, "psm_iptw", pattern)
            if df is None:
                continue
            dose_col = next((c for c in ("crrt_dose_median_3h",
                                          "crrt_dose_ml_kg_hr_0")
                             if c in df.columns), None)
            if dose_col is None:
                continue
            df = df[df[dose_col].between(*DOSE_XLIM)]
            if df.empty:
                continue
            lbl = labels.get(sd.name, sd.name)
            color = get_color(sd.name)
            ax.plot(df[dose_col], df["hr"], label=lbl,
                    color=color, linewidth=2)
            ax.fill_between(df[dose_col],
                            df["hr_lower"], df["hr_upper"],
                            alpha=0.15, color=color)
            any_data = True

        if any_data:
            ax.axhline(y=1, color="grey", linestyle="--", linewidth=1)
            ax.set_xlim(*DOSE_XLIM)
            ax.set_xlabel(f"CRRT Dose (mL/kg/hr) — display clipped to {DOSE_XLIM[0]}-{DOSE_XLIM[1]}")
            ax.set_ylabel("Hazard Ratio")
            ax.set_title(title)
            ax.legend()
            ax.grid(alpha=0.3)
            fig.tight_layout(pad=PAD_DEFAULT)
            html_parts.append(
                f'<div class="figure-block"><h3>{title}</h3>'
                f'<img src="{fig_to_data_uri(fig)}" alt="{title}"></div>'
            )
        else:
            plt.close(fig)
            html_parts.append(_missing(f"{title}: data not available."))

    # Dose decile mortality overlay
    fig, ax = plt.subplots(figsize=(10, 5))
    any_data = False
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "dose_decile_mortality")
        if df is None:
            continue
        lbl = labels.get(sd.name, sd.name)
        color = get_color(sd.name)
        ax.plot(df["dose_median"], df["mortality_rate"], "o-",
                label=lbl, color=color, linewidth=2, markersize=5)
        any_data = True
    if any_data:
        ax.set_xlabel("CRRT Dose (mL/kg/hr) — Decile Median")
        ax.set_ylabel("30-day Mortality Rate")
        ax.set_title("Dose-Decile Mortality")
        ax.legend()
        ax.grid(alpha=0.3)
        fig.tight_layout(pad=PAD_DEFAULT)
        html_parts.append(
            '<div class="figure-block"><h3>Dose-Decile Mortality</h3>'
            f'<img src="{fig_to_data_uri(fig)}" alt="Dose decile"></div>'
        )
    else:
        plt.close(fig)

    if not html_parts:
        return _missing("Dose-response data not available.")
    return ('<div class="section"><h2>Dose-Response Analysis</h2>'
            '<p class="fig-caption">Continuous relationship between CRRT dose and 30-day mortality, assessed via three complementary approaches.</p>'
            + "\n".join(html_parts) + '</div>')




def _subgroup_forest(subgroup_name, lv_a, lv_b,
                     site_labels, site_ids,
                     hrs_a, los_a, his_a,
                     hrs_b, los_b, his_b,
                     pooled_a, pooled_b) -> str:
    """Compact side-by-side forest plot for a two-level subgroup.

    Shows site-level HRs for both levels as colored markers (circles vs squares)
    with RE pooled diamonds at the bottom.
    """
    k = len(site_labels)
    fig, ax = plt.subplots(figsize=(8.5, k * 0.85 + 2.8))

    color_a = "#f9776d"  # coral — matches CIF high-dose curve
    color_b = "#2cbfc3"  # teal — matches CIF low-dose curve
    offset = 0.15  # vertical offset to separate the two markers per site

    for i, (_lbl, _sid) in enumerate(zip(site_labels, site_ids)):
        # Level A (circle, slightly above)
        if i < len(hrs_a) and not np.isnan(hrs_a[i]):
            ax.errorbar(hrs_a[i], i - offset,
                        xerr=[[hrs_a[i] - los_a[i]], [his_a[i] - hrs_a[i]]],
                        fmt="o", color=color_a, markersize=8, capsize=4, linewidth=1.8)
        # Level B (square, slightly below)
        if i < len(hrs_b) and not np.isnan(hrs_b[i]):
            ax.errorbar(hrs_b[i], i + offset,
                        xerr=[[hrs_b[i] - los_b[i]], [his_b[i] - hrs_b[i]]],
                        fmt="s", color=color_b, markersize=8, capsize=4, linewidth=1.8)

    # Pooled diamonds
    y_base = k + 0.4
    for j, (pool, color, _level_name) in enumerate([
        (pooled_a, color_a, lv_a), (pooled_b, color_b, lv_b)
    ]):
        if pool is None:
            continue
        y = y_base + j * 0.7
        dx = [pool["re_lo"], pool["re_hr"], pool["re_hi"], pool["re_hr"]]
        dy = [y, y - 0.2, y, y + 0.2]
        ax.fill(dx, dy, color=color, alpha=0.8)
        ax.text(max(pool["re_hi"] + 0.02, 1.05), y,
                f'{pool["re_hr"]:.2f} [{pool["re_lo"]:.2f}\u2013{pool["re_hi"]:.2f}]',
                va="center", fontsize=11, color=color, fontweight="bold")

    ax.axvline(1, color="grey", linestyle="--", linewidth=0.8)
    yticks = list(range(k)) + [y_base, y_base + 0.7]
    ylabels = list(site_labels) + [f"Pooled: {lv_a}", f"Pooled: {lv_b}"]
    # Trim if a pooled level is missing
    if pooled_b is None:
        yticks = yticks[:-1]
        ylabels = ylabels[:-1]
    if pooled_a is None:
        yticks = [y for y in yticks if y != y_base]
        ylabels = [l for l in ylabels if not l.startswith("Pooled:")]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.invert_yaxis()
    ax.set_xlabel("HR (95% CI)")
    ax.set_title(subgroup_name)
    ax.grid(axis="x", alpha=0.2)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color=color_a, label=lv_a,
               markersize=8, linestyle="None"),
        Line2D([0], [0], marker="s", color=color_b, label=lv_b,
               markersize=8, linestyle="None"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=11,
              framealpha=0.9, edgecolor="#e2e8f0")

    fig.tight_layout(pad=PAD_DEFAULT)
    fig.subplots_adjust(left=0.28, right=0.96)
    uri = fig_to_data_uri(fig)
    return f'<img src="{uri}" style="max-width:100%; height:auto;" alt="Forest: {subgroup_name}">'


def build_subgroup_analysis(sites, labels) -> str:
    """Subgroup interaction analysis grouped by subgroup, pivoted by level.

    Each subgroup becomes a side-by-side section: table on the left with
    site-specific and pooled HRs, forest plot on the right.
    """
    site_data = []  # (site_label, site_id, dataframe)
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "subgroup_analysis_results")
        if df is None:
            continue
        site_data.append((labels.get(sd.name, sd.name), sd.name, df))

    if not site_data:
        return _missing("Subgroup analysis not available.")

    # Determine ordered list of subgroups (preserve original order from first site)
    subgroup_order = []
    for _, _, df in site_data:
        for s in df["subgroup"].unique():
            if s not in subgroup_order:
                subgroup_order.append(s)

    def _pool_level(subgroup_name, level_name):
        """Collect HRs for one subgroup level across sites, return arrays + meta result."""
        log_hrs, se_hrs, hrs, los, his = [], [], [], [], []
        s_labels, s_ids = [], []
        for lbl, sid, df in site_data:
            sub = df[(df["subgroup"] == subgroup_name) & (df["level"] == level_name)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            hr, lo, hi = r["HR"], r["HR_lower"], r["HR_upper"]
            if pd.isna(hr) or hr <= 0 or lo <= 0 or hi <= 0:
                continue
            se = _se_from_ci(hr, lo, hi)
            if np.isnan(se) or se <= 0:
                continue
            log_hrs.append(np.log(hr))
            se_hrs.append(se)
            hrs.append(hr)
            los.append(lo)
            his.append(hi)
            s_labels.append(lbl)
            s_ids.append(sid)
        result = None
        if len(log_hrs) >= 2:
            result = meta_analyze(np.array(log_hrs), np.array(se_hrs))
        return s_labels, s_ids, np.array(hrs), np.array(los), np.array(his), result

    # ── Overall section (single-level) ──────────────────────────────────
    overall_rows = []
    ov_log_hrs, ov_se_hrs = [], []
    for lbl, sid, df in site_data:
        ov = df[df["subgroup"].str.lower() == "overall"]
        if ov.empty:
            continue
        r = ov.iloc[0]
        color = get_color(sid)
        overall_rows.append(
            f'<tr style="border-left:4px solid {color};">'
            f'<td>{lbl}</td>'
            f'<td>{int(r["n"])}</td><td>{int(r["events"])}</td>'
            f'<td>{_fmt_hr(r)}</td>'
            f'</tr>'
        )
        hr, lo, hi = r["HR"], r["HR_lower"], r["HR_upper"]
        if not pd.isna(hr) and hr > 0 and lo > 0 and hi > 0:
            se = _se_from_ci(hr, lo, hi)
            if not np.isnan(se) and se > 0:
                ov_log_hrs.append(np.log(hr))
                ov_se_hrs.append(se)

    # Pooled overall row
    if len(ov_log_hrs) >= 2:
        ov_meta = meta_analyze(np.array(ov_log_hrs), np.array(ov_se_hrs))
        pooled_hr_str = (f'{ov_meta["re_hr"]:.2f} '
                         f'({ov_meta["re_lo"]:.2f}\u2013{ov_meta["re_hi"]:.2f})')
        overall_rows.append(
            f'<tr style="background:#f0fdfa; font-weight:600;">'
            f'<td>Pooled (RE)</td>'
            f'<td></td><td></td>'
            f'<td>{pooled_hr_str}</td>'
            f'</tr>'
        )

    overall_html = (
        '<table class="results-table" border="0">'
        '<thead><tr><th>Site</th><th>N</th><th>Events</th>'
        '<th>HR (95% CI)</th></tr></thead>'
        '<tbody>' + "\n".join(overall_rows) + '</tbody></table>'
    )

    # ── Two-level subgroup sections ─────────────────────────────────────
    subgroup_sections = []
    for subgroup in subgroup_order:
        if subgroup.lower() == "overall":
            continue
        # Collect levels across all sites to determine column order
        all_levels = []
        for _, _, df in site_data:
            sub = df[df["subgroup"] == subgroup]
            for lv in sub["level"].unique():
                if lv not in all_levels:
                    all_levels.append(lv)
        if len(all_levels) < 2:
            continue
        lv_a, lv_b = all_levels[0], all_levels[1]

        # Skip subgroups where fewer than 2 sites have valid HRs in BOTH levels
        # (e.g., CRRT Modality when most sites use a single modality)
        sites_with_both = 0
        for _, _, df_chk in site_data:
            sub_chk = df_chk[df_chk["subgroup"] == subgroup]
            ra_chk = sub_chk[sub_chk["level"] == lv_a]
            rb_chk = sub_chk[sub_chk["level"] == lv_b]
            has_a = not ra_chk.empty and not pd.isna(ra_chk.iloc[0]["HR"]) and ra_chk.iloc[0]["n"] > 0
            has_b = not rb_chk.empty and not pd.isna(rb_chk.iloc[0]["HR"]) and rb_chk.iloc[0]["n"] > 0
            if has_a and has_b:
                sites_with_both += 1
        if sites_with_both < 2:
            continue

        # Pool each level (only the pooled estimate is consumed here; the per-site
        # arrays are recomputed below from the raw data with NaN padding for sites
        # missing one of the two levels — see all_hrs_a / all_hrs_b further down).
        *_, pooled_a = _pool_level(subgroup, lv_a)
        *_, pooled_b = _pool_level(subgroup, lv_b)

        # Pooled interaction p-value (meta-analyze the log(HR) difference across sites)
        pooled_p_int = ""
        diff_log_hrs, diff_ses = [], []
        for _, _, df in site_data:
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
        if len(diff_log_hrs) >= 2:
            diff_meta = meta_analyze(np.array(diff_log_hrs), np.array(diff_ses))
            pooled_p_int = _fmt_p(diff_meta["re_p"])

        # Build header
        hdr = (
            f'<th>Site</th>'
            f'<th style="text-align:right;">N</th><th>HR ({lv_a})</th>'
            f'<th style="text-align:right;">N</th><th>HR ({lv_b})</th>'
            f'<th>p-int</th><th>q-int (BH)</th>'
        )
        rows_html = []
        for lbl, sid, df in site_data:
            sub = df[df["subgroup"] == subgroup]
            if sub.empty:
                rows_html.append(
                    f'<tr style="border-left:4px solid {get_color(sid)};">'
                    f'<td>{lbl}</td>' + '<td></td>' * 6 + '</tr>'
                )
                continue
            color = get_color(sid)
            ra = sub[sub["level"] == lv_a]
            rb = sub[sub["level"] == lv_b]
            n_a = str(int(ra.iloc[0]["n"])) if not ra.empty else ""
            n_b = str(int(rb.iloc[0]["n"])) if not rb.empty else ""
            hr_a = _fmt_hr(ra.iloc[0]) if not ra.empty else ""
            hr_b = _fmt_hr(rb.iloc[0]) if not rb.empty else ""
            p_int = ""
            q_int = ""
            for piece in [ra, rb]:
                if not piece.empty:
                    if pd.notna(piece.iloc[0].get("p_interaction")):
                        p_int = _fmt_p(piece.iloc[0]["p_interaction"])
                    if pd.notna(piece.iloc[0].get("q_interaction")):
                        q_int = _fmt_p(piece.iloc[0]["q_interaction"])
                    if p_int:
                        break
            rows_html.append(
                f'<tr style="border-left:4px solid {color};">'
                f'<td>{lbl}</td>'
                f'<td style="text-align:right;">{n_a}</td><td>{hr_a}</td>'
                f'<td style="text-align:right;">{n_b}</td><td>{hr_b}</td>'
                f'<td>{p_int}</td><td>{q_int}</td>'
                f'</tr>'
            )

        # Pooled row
        pooled_hr_a = ""
        pooled_hr_b = ""
        if pooled_a:
            pooled_hr_a = (f'{pooled_a["re_hr"]:.2f} '
                           f'({pooled_a["re_lo"]:.2f}\u2013{pooled_a["re_hi"]:.2f})')
        if pooled_b:
            pooled_hr_b = (f'{pooled_b["re_hr"]:.2f} '
                           f'({pooled_b["re_lo"]:.2f}\u2013{pooled_b["re_hi"]:.2f})')
        rows_html.append(
            f'<tr style="background:#f0fdfa; font-weight:600;">'
            f'<td>Pooled (RE)</td>'
            f'<td></td><td>{pooled_hr_a}</td>'
            f'<td></td><td>{pooled_hr_b}</td>'
            f'<td>{pooled_p_int}</td><td></td>'
            f'</tr>'
        )

        table_html = (
            f'<table class="results-table" border="0">'
            f'<thead><tr>{hdr}</tr></thead>'
            f'<tbody>' + "\n".join(rows_html) + '</tbody></table>'
        )

        # Forest plot — need aligned arrays across all sites (NaN for missing)
        all_hrs_a, all_los_a, all_his_a = [], [], []
        all_hrs_b, all_los_b, all_his_b = [], [], []
        all_labels, all_sids = [], []
        for lbl, sid, df in site_data:
            all_labels.append(lbl)
            all_sids.append(sid)
            sub = df[df["subgroup"] == subgroup]
            ra = sub[sub["level"] == lv_a]
            rb = sub[sub["level"] == lv_b]
            if not ra.empty and not pd.isna(ra.iloc[0]["HR"]):
                all_hrs_a.append(ra.iloc[0]["HR"])
                all_los_a.append(ra.iloc[0]["HR_lower"])
                all_his_a.append(ra.iloc[0]["HR_upper"])
            else:
                all_hrs_a.append(float("nan"))
                all_los_a.append(float("nan"))
                all_his_a.append(float("nan"))
            if not rb.empty and not pd.isna(rb.iloc[0]["HR"]):
                all_hrs_b.append(rb.iloc[0]["HR"])
                all_los_b.append(rb.iloc[0]["HR_lower"])
                all_his_b.append(rb.iloc[0]["HR_upper"])
            else:
                all_hrs_b.append(float("nan"))
                all_los_b.append(float("nan"))
                all_his_b.append(float("nan"))

        forest_html = _subgroup_forest(
            subgroup, lv_a, lv_b,
            all_labels, all_sids,
            np.array(all_hrs_a), np.array(all_los_a), np.array(all_his_a),
            np.array(all_hrs_b), np.array(all_los_b), np.array(all_his_b),
            pooled_a, pooled_b,
        )

        subgroup_sections.append(
            f'<h3 style="margin:20px 0 6px; font-size:14px; color:#334155;">'
            f'{subgroup}</h3>'
            f'<div style="display:flex; gap:40px; align-items:flex-start; flex-wrap:wrap;">'
            f'<div style="flex:0 0 auto;">{table_html}</div>'
            f'<div style="flex:0 0 auto; max-width:480px;">{forest_html}</div>'
            f'</div>'
        )

    return (
        '<div class="section"><h2>Subgroup Analysis</h2>'
        '<p class="fig-caption">Hypothesis-generating interaction analysis '
        'testing whether high CRRT dose (\u226530 mL/kg/hr) affects 30-day '
        'mortality differently across pre-specified clinical subgroups. '
        'HRs reflect the effect of high vs. low dose within each '
        'subpopulation; p-interaction tests whether the dose effect varies '
        'by subgroup level. q-values are Benjamini\u2013Hochberg adjusted '
        'to account for multiple comparisons across subgroups. '
        'Pooled estimates use DerSimonian\u2013Laird random-effects meta-analysis; '
        'pooled p-int tests the cross-site pooled interaction.</p>'
        '<h3 style="margin:16px 0 6px; font-size:14px; color:#334155;">Overall</h3>'
        + overall_html
        + "".join(subgroup_sections)
        + '</div>'
    )


def build_evalue(sites, labels) -> str:
    """E-value sensitivity analysis table."""
    rows = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "evalue_sensitivity")
        if df is None:
            continue
        # Drop MICE-pooled duplicate rows (see build_model_comparison).
        if "model" in df.columns:
            df = df[~df["model"].astype(str).str.contains(
                r"\(pooled\)", regex=True, na=False
            )]
        lbl = labels.get(sd.name, sd.name)
        for _, r in df.iterrows():
            rows.append({
                "site_id": sd.name,
                "Site": lbl,
                "Model": r.get("model", ""),
                "HR": f"{r['HR']:.2f}",
                "E-value (point)": f"{r['evalue_point']:.2f}",
                "E-value (CI)": f"{r['evalue_ci']:.2f}",
            })
    if not rows:
        return _missing("E-value results not available.")
    tbl = pd.DataFrame(rows)
    return (
        '<div class="section"><h2>E-Value Sensitivity Analysis</h2>'
        '<p class="fig-caption">Minimum strength of unmeasured confounding (on the risk ratio scale) needed to explain away the observed result.</p>'
        + _styled_table(tbl)
        + '</div>'
    )


def build_diagnostics_iptw(sites, labels) -> str:
    """PSM/IPTW diagnostics (t=0 point-treatment): ESS, PS overlap, love plots, PS extremes."""
    html_parts = []

    # ESS table
    ess_rows = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "SL_IPTW_ESS")
        if df is not None:
            r = df.iloc[0]
            ess_rows.append({
                "site_id": sd.name,
                "Site": labels.get(sd.name, sd.name),
                "N": int(r.get("N", 0)),
                "ESS": f"{r.get('ESS', 0):.0f}",
                "ESS %": f"{r.get('ESS_prop', 0) * 100:.1f}",
            })
    if ess_rows:
        tbl = pd.DataFrame(ess_rows)
        html_parts.append(
            '<h3>Effective Sample Size</h3>'
            + _styled_table(tbl)
        )

    # PS extremes table
    ext_rows = []
    for sd in sites:
        df = _load_csv(sd, "psm_iptw", "PS_extremes")
        if df is None:
            continue
        lbl = labels.get(sd.name, sd.name)
        for _, r in df.iterrows():
            direction = str(r.get("direction", ""))
            ext_rows.append({
                "site_id": sd.name,
                "Site": lbl,
                "Threshold": direction,
                "N": int(r.get("n", r.get("N", 0))),
                "%": f"{r.get('percentage', r.get('pct', 0)):.2f}",
            })
    if ext_rows:
        tbl = pd.DataFrame(ext_rows)
        html_parts.append(
            '<h3>Propensity Score Extremes</h3>'
            + _styled_table(tbl)
        )

    # PS overlap and love plots
    for pattern, title in [
        ("SL_PS_overlap", "Propensity Score Overlap"),
        ("SL_LovePlot", "Covariate Balance (Love Plot)"),
    ]:
        imgs = []
        for sd in sites:
            png = _find_png(sd, "psm_iptw", pattern)
            if png is None:
                continue
            lbl = labels.get(sd.name, sd.name)
            imgs.append(
                f'<div style="display:inline-block; margin:10px;">'
                f'<h4>{lbl}</h4>'
                f'<img src="{img_to_data_uri(png)}" style="max-width:450px;" alt="{title} {lbl}"></div>'
            )
        if imgs:
            html_parts.append(
                f'<div class="figure-block"><h3>{title}</h3>' + "".join(imgs) + '</div>'
            )

    if not html_parts:
        return _missing("IPTW diagnostics not available.")
    return (
        '<div class="section"><h2>PSM/IPTW Diagnostics (Point-Treatment, t=0)</h2>'
        '<p><em>These diagnostics assess balance at a single time point (CRRT initiation). '
        'ESS retention is high and PS overlap is generally adequate.</em></p>'
        + "\n".join(html_parts) + '</div>'
    )


# ── Meta-Analysis ───────────────────────────────────────────────────────────



def _forest_plot(site_labels, hrs, lowers, uppers, site_ids, meta_result, title) -> str:
    """Draw a forest plot with site estimates + pooled diamonds, return HTML."""
    k = len(site_labels)
    fig, (ax, ax_txt) = plt.subplots(1, 2, figsize=(15, k * 1.0 + 3.0),
                                      gridspec_kw={"width_ratios": [3, 1]})
    ax_txt.axis("off")  # right panel is text-only

    # Individual sites
    for i, (_lbl, hr, lo, hi, sid) in enumerate(zip(site_labels, hrs, lowers, uppers, site_ids)):
        color = get_color(sid)
        ax.errorbar(hr, i, xerr=[[hr - lo], [hi - hr]],
                    fmt="s", color=color, markersize=10, capsize=6, linewidth=2.2)
        ax_txt.text(0.05, i, f"{hr:.2f} [{lo:.2f}, {hi:.2f}]",
                    va="center", fontsize=12, transform=ax_txt.get_yaxis_transform())

    # Random-effects pooled diamond
    re = meta_result
    y_re = k + 0.6
    diamond_x = [re["re_lo"], re["re_hr"], re["re_hi"], re["re_hr"]]
    diamond_y = [y_re, y_re - 0.25, y_re, y_re + 0.25]
    ax.fill(diamond_x, diamond_y, color="#D32F2F", alpha=0.8)
    ax_txt.text(0.05, y_re,
                f"RE: {re['re_hr']:.2f} [{re['re_lo']:.2f}, {re['re_hi']:.2f}] p={_fmt_p(re['re_p'])}",
                va="center", fontsize=12, fontweight="bold",
                transform=ax_txt.get_yaxis_transform())

    # Fixed-effect pooled diamond
    y_fe = k + 1.3
    diamond_x = [re["fe_lo"], re["fe_hr"], re["fe_hi"], re["fe_hr"]]
    diamond_y = [y_fe, y_fe - 0.25, y_fe, y_fe + 0.25]
    ax.fill(diamond_x, diamond_y, color="#1976D2", alpha=0.8)
    ax_txt.text(0.05, y_fe,
                f"FE: {re['fe_hr']:.2f} [{re['fe_lo']:.2f}, {re['fe_hi']:.2f}] p={_fmt_p(re['fe_p'])}",
                va="center", fontsize=12, fontweight="bold",
                transform=ax_txt.get_yaxis_transform())

    ax.axvline(1, color="grey", linestyle="--", linewidth=1)
    all_yticks = list(range(k)) + [y_re, y_fe]
    ax.set_yticks(all_yticks)
    ax.set_yticklabels(list(site_labels) + ["RE Pooled", "FE Pooled"])
    ax.invert_yaxis()
    ax.set_xlabel("Hazard Ratio (95% CI)")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.3)

    # Sync y-axis limits on text panel
    ax_txt.set_ylim(ax.get_ylim())

    # Heterogeneity annotation
    ax.text(0.02, 0.02,
            f"I\u00b2 = {re['I2']:.1f}%    \u03c4\u00b2 = {re['tau2']:.4f}    Q p = {_fmt_p(re['Q_p'])}",
            transform=ax.transAxes, fontsize=11, color="#666")
    fig.tight_layout(pad=PAD_DEFAULT)
    uri = fig_to_data_uri(fig)
    return (
        f'<div class="figure-block"><h3>{title}</h3>'
        f'<img src="{uri}" alt="{title}"></div>'
    )






def build_meta_analysis(sites, labels) -> str:
    """Meta-analysis tab: pool site-level HRs across all causal analyses."""
    html_parts = []
    summary_rows = []

    # Define all analyses to pool
    analyses = [
        # (title, collector_func, kwargs)
        # MSM analyses removed: the time-varying MSM was cut from the
        # manuscript (positivity violations) and archived (tag msm-v1).
        ("Point-Treatment IPTW \u2014 Death",
         _collect_iptw, dict(outcome_kw="Death")),
        ("Point-Treatment IPTW \u2014 Discharge",
         _collect_iptw, dict(outcome_kw="Discharge")),

        ("Point-Treatment PSM Fine-Gray \u2014 Death",
         _collect_psm_fg, dict(outcome_kw="Death")),
        ("Point-Treatment PSM Fine-Gray \u2014 Discharge",
         _collect_psm_fg, dict(outcome_kw="Discharge")),
    ]

    for title, collector, kwargs in analyses:
        s_lbls, s_ids, log_h, se_h, h, lo, hi = collector(sites, labels, **kwargs)
        if len(log_h) < 2:
            continue
        result = meta_analyze(log_h, se_h)
        html_parts.append(_forest_plot(s_lbls, h, lo, hi, s_ids, result, title))
        summary_rows.append({
            "Analysis": title,
            "Sites": len(log_h),
            "FE HR (95% CI)": f"{result['fe_hr']:.2f} ({result['fe_lo']:.2f}\u2013{result['fe_hi']:.2f})",
            "FE p": _fmt_p(result["fe_p"]),
            "RE HR (95% CI)": f"{result['re_hr']:.2f} ({result['re_lo']:.2f}\u2013{result['re_hi']:.2f})",
            "RE p": _fmt_p(result["re_p"]),
            "I\u00b2 %": f"{result['I2']:.1f}",
            "Q p": _fmt_p(result["Q_p"]),
        })

    # Summary table — grouped by analysis type
    if summary_rows:
        # Classify each row into a (group_label, short_name) pair. The ordering
        # below — Point-Treatment IPTW, then Point-Treatment PSM Fine-Gray —
        # also defines the display order in the rendered table.
        def _classify(title: str) -> tuple[str, str]:
            """Return (group_label, short_name) for a summary row."""
            if title.startswith("Point-Treatment PSM Fine-Gray"):
                outcome = title.split("\u2014")[-1].strip()
                return ("Point-Treatment PSM Fine-Gray", outcome)
            if title.startswith("Point-Treatment IPTW"):
                outcome = title.split("\u2014")[-1].strip()
                return ("Point-Treatment IPTW", outcome)
            return ("Other", title)

        cols = ["Outcome", "Sites", "FE HR (95% CI)", "FE p",
                "RE HR (95% CI)", "RE p", "I\u00b2 %", "Q p"]
        ncols = len(cols)
        header = "".join(f"<th>{c}</th>" for c in cols)

        rows_html = []
        prev_group = None
        for row in summary_rows:
            group, short_name = _classify(row["Analysis"])
            if group != prev_group:
                rows_html.append(
                    f'<tr><td colspan="{ncols}" style="background:#f0fdfa; '
                    f'font-weight:600; color:#0f766e; padding:10px 12px; '
                    f'border-bottom:2px solid #ccfbf1; border-top:1px solid #e2e8f0; '
                    f'font-size:13px; letter-spacing:0.2px;">{group}</td></tr>'
                )
                prev_group = group
            tds = (
                f'<td style="padding-left:24px;">{short_name}</td>'
                + f'<td>{row["Sites"]}</td>'
                + f'<td>{row["FE HR (95% CI)"]}</td>'
                + f'<td>{row["FE p"]}</td>'
                + f'<td>{row["RE HR (95% CI)"]}</td>'
                + f'<td>{row["RE p"]}</td>'
                + f'<td>{row["I² %"]}</td>'
                + f'<td>{row["Q p"]}</td>'
            )
            rows_html.append(f"<tr>{tds}</tr>")

        summary_html = (
            '<h3>Summary of All Pooled Estimates</h3>'
            f'<table class="results-table" border="0">'
            f'<thead><tr>{header}</tr></thead>'
            f'<tbody>{"".join(rows_html)}</tbody></table>'
        )
        html_parts.insert(0, summary_html)

    if not html_parts:
        return _missing("Meta-analysis not available (need >=2 sites with se_log_hr).")
    return ('<div class="section"><h2>Meta-Analysis (Inverse-Variance Weighted)</h2>'
            '<p class="fig-caption">Fixed-effect and DerSimonian-Laird random-effects pooled estimates across all participating sites.</p>'
            + "\n".join(html_parts) + '</div>')


# ── Combined Tables ────────────────────────────────────────────────────────

























def build_combined_consort(sites, labels) -> str:
    """Pooled causal cohort exclusion flow diagram across all sites.

    Mirrors the per-site causal cohort flow emitted by step 04: strobe-based
    exclusions (CRRT → ESRD → weight → settings → labs) feed into one
    additional causal step (died/off-CRRT-≤24h + SCUF-only) to yield the
    final causal analysis cohort.

    Data sources per site:
      - final/crrt_epi/{site}_strobe_counts.csv   (descriptive steps 1-6)
      - final/psm_iptw/{site}_psm_counts_summary.csv  (n_causal)

    Sites missing psm_counts_summary are shown in the per-site breakdown
    table with '—' in the causal column but are excluded from the pooled
    diagram and pooled-total row (otherwise the final box would be wrong)."""
    from matplotlib.patches import FancyBboxPatch

    site_data: list[tuple[str, str, dict, int | None]] = []
    for sd in sites:
        sid = sd.name
        lbl = labels.get(sid, sid)
        df = _load_csv(sd, "crrt_epi", "strobe_counts")
        if df is None or df.empty or "counter" not in df.columns or "value" not in df.columns:
            continue
        counts = dict(zip(df["counter"].astype(str), df["value"]))
        n_causal = _read_site_causal_n(sd)
        site_data.append((sid, lbl, counts, n_causal))

    if len(site_data) < 2:
        return _missing("Combined cohort flow diagram requires strobe_counts.csv from at least 2 sites.")

    # Sites contributing to the pooled diagram = those with a causal N.
    full = [(sid, lbl, c, nc) for sid, lbl, c, nc in site_data if nc is not None]
    if len(full) < 2:
        return _missing(
            "Combined causal cohort flow diagram requires psm_counts_summary.csv from at "
            "least 2 sites. Ask sites to rerun step 05 to populate this file."
        )

    def pool(key: str) -> int:
        return int(sum(c.get(key, 0) or 0 for _, _, c, _ in full))

    pooled_counts = {
        "1b_after_stitching": pool("1b_after_stitching"),
        "2_crrt_blocks": pool("2_crrt_blocks"),
        "3_encounter_blocks_without_esrd": pool("3_encounter_blocks_without_esrd"),
        "4_encounter_blocks_with_weight": pool("4_encounter_blocks_with_weight"),
        "5_encounter_blocks_with_crrt_settings": pool("5_encounter_blocks_with_crrt_settings"),
        "6_encounter_blocks_with_required_labs": pool("6_encounter_blocks_with_required_labs"),
        "causal": int(sum(nc for _, _, _, nc in full)),
    }
    start_n = pooled_counts["1b_after_stitching"]

    # Step plan: (remaining_n, remaining_label, excluded_label)
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

    steps = []
    parent = start_n
    for remaining, remaining_lbl, excl_lbl in step_defs:
        excluded = max(parent - remaining, 0)
        # Skip "Settings" step if no exclusions (matches single-site behavior).
        if excl_lbl.endswith("CRRT settings") and excluded == 0:
            continue
        steps.append({
            "remaining_n": remaining,
            "remaining_label": remaining_lbl,
            "excluded_n": excluded,
            "excluded_label": excl_lbl,
        })
        parent = remaining

    # ── Render diagram ──────────────────────────────────────────────────────
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
    v_spacing = 0.86 / n_boxes  # fit all boxes between y=0 and y=0.90

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

    n_sites = len(full)
    ax.text(0.5, 0.98,
            f"Pooled Cohort Selection ({n_sites} CLIF sites)",
            ha="center", va="center", fontsize=18, fontweight="bold")
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

    diagram_uri = fig_to_data_uri(fig)

    # ── Per-site breakdown table ────────────────────────────────────────────
    display_labels = [
        ("1b_after_stitching", "Adult hosp."),
        ("2_crrt_blocks", "CRRT blocks"),
        ("3_encounter_blocks_without_esrd", "No ESRD"),
        ("4_encounter_blocks_with_weight", "With weight & sex"),
        ("5_encounter_blocks_with_crrt_settings", "With CRRT settings"),
        ("6_encounter_blocks_with_required_labs", "With required labs"),
    ]

    # Include every site (not just sites with causal N) in the breakdown.
    rows = []
    for sid, lbl, counts, nc in site_data:
        row = {"site_id": sid, "Site": lbl}
        for k, col_lbl in display_labels:
            v = counts.get(k)
            row[col_lbl] = f"{int(v):,}" if v is not None and not pd.isna(v) else "—"
        row["Causal cohort"] = f"{nc:,}" if nc is not None else "—"
        rows.append(row)

    total_row = {"site_id": "", "Site": "Pooled total"}
    for k, col_lbl in display_labels:
        total_row[col_lbl] = f"{pooled_counts[k]:,}"
    total_row["Causal cohort"] = f"{pooled_counts['causal']:,}"
    rows.append(total_row)

    display_cols = ["Site"] + [lbl for _, lbl in display_labels] + ["Causal cohort"]
    header = "".join(f"<th>{c}</th>" for c in display_cols)
    tr_rows = []
    for r in rows:
        sid = r["site_id"]
        is_total = (r["Site"] == "Pooled total")
        if is_total:
            tds = "".join(f"<td><strong>{r[c]}</strong></td>" for c in display_cols)
            tr_rows.append(
                f'<tr style="background:#f1f5f9; border-top:2px solid #0f172a;">{tds}</tr>'
            )
        else:
            color = get_color(sid) if sid else "#cbd5e1"
            tds = "".join(f"<td>{r[c]}</td>" for c in display_cols)
            tr_rows.append(f'<tr style="border-left:4px solid {color};">{tds}</tr>')

    breakdown_table = (
        '<table class="results-table" border="0">'
        f'<thead><tr>{header}</tr></thead>'
        '<tbody>' + "\n".join(tr_rows) + '</tbody></table>'
    )

    missing_note = ""
    missing_sites = [lbl for _sid, lbl, _, nc in site_data if nc is None]
    if missing_sites:
        missing_note = (
            f'<p style="color:#b45309; font-size:13px;"><strong>Note:</strong> '
            f'{", ".join(missing_sites)} missing psm_counts_summary.csv — '
            f'excluded from the pooled diagram but shown in the breakdown below. '
            f'Ask these sites to rerun script 05 to populate the causal N.</p>'
        )

    return (
        '<div class="section"><h2>Combined Causal Cohort Exclusion Flow Diagram</h2>'
        '<p><em>Patient flow across all participating sites ending at the '
        'causal analysis cohort. Counts are at the encounter-block level. '
        'The final step combines the died/off-CRRT-≤24h and SCUF-only '
        'exclusions applied in script 04.</em></p>'
        + missing_note +
        f'<div style="text-align:center; margin:16px 0;">'
        f'<img src="{diagram_uri}" style="max-width:900px; width:100%;" '
        'alt="Pooled causal cohort exclusion flow diagram"/></div>'
        '<h3 style="margin-top:24px;">Per-Site Breakdown</h3>'
        '<p><em>Encounter-block counts at each selection step. The final '
        'column is the causal cohort N from each site\'s psm_counts_summary.csv. '
        'The bottom row is the pooled total used in the diagram above.</em></p>'
        + breakdown_table +
        '</div>'
    )


def build_combined_tables(sites, labels) -> str:
    """Combined Table 1, S1, S2 across sites using weighted medians."""
    html_parts = []

    # Table 1 — single project baseline table: full CRRT cohort, OVERALL
    # (not dose-banded), produced by script 02 (crrt_epi/table1_crrt). The
    # dose-band variation is shown as a figure, not as Table 1 strata.
    # (Sites still on the old survival-format table1_crrt are skipped by the
    # parser until they re-run 02 — no "Overall" column to anchor on.)
    html_parts.append(
        _combine_table_across_sites(
            sites, labels, "crrt_epi", "table1_crrt",
            "Table 1 — Baseline Characteristics of the CRRT Cohort",
            overall_only=True, include_site_range=True)
    )

    # Causal-cohort balance tables (script 05): Table 2 = unadjusted balance
    # between dose groups; Tables S1/S2 = post-adjustment (PSM / IPTW) balance.
    for pattern, title in [
        (BALANCE_TABLE_PATTERNS, "Table 2 — Unadjusted Balance Between Dose Groups (Causal Cohort)"),
        ("TableS1_matched", "Table S1 — PSM-Matched Cohort"),
        ("TableS2_IPTW", "Table S2 — IPTW-Weighted Cohort"),
    ]:
        html_parts.append(
            _combine_table_across_sites(sites, labels, "psm_iptw", pattern, title)
        )
    # (MSM combined tables removed — the time-varying MSM analysis was archived.)

    consort_html = build_combined_consort(sites, labels)

    tables_body = ""
    if any(p and "not available" not in p for p in html_parts):
        tables_body = (
            '<div class="section"><h2>Combined Tables</h2>'
            '<p><em>Tables below combine summary statistics from all sites. '
            'Continuous variables use weighted medians (weighted by site N). '
            'Categorical variables are pooled counts with recomputed percentages.</em></p>'
            + "\n".join(html_parts)
            + '</div>'
        )

    if not consort_html and not tables_body:
        return _missing("Combined tables not available.")
    return consort_html + tables_body


# Pattern → (title, one-line description) for per-site PNG figures.
# Patterns are matched case-insensitively against the filename stem (no extension).
# Order matters: first match wins.
_FIGURE_CATALOG = [
    # ── PSM / IPTW ──
    ("causal_consort_diagram",
     "Causal Cohort Exclusion Flow Diagram",
     "Patient flow from the descriptive cohort through causal analysis exclusion criteria."),
    ("subgroup_forest_plot",
     "Subgroup Interaction Forest Plot",
     "Hazard ratios for high vs. low CRRT dose within pre-specified clinical subgroups."),
    ("SL_LovePlot",
     "Covariate Balance After IPTW (Love Plot)",
     "Standardized mean differences before and after inverse probability weighting via Super Learner."),
    ("SL_PS_overlap",
     "Propensity Score Overlap",
     "Distribution of estimated propensity scores by treatment group, assessing positivity."),
    ("IPTW_CIF_Death",
     "IPTW Cumulative Incidence of Death",
     "Weighted cumulative incidence curves for mortality, comparing high vs. low CRRT dose."),
    ("IPTW_CIF_Discharge",
     "IPTW Cumulative Incidence of Discharge",
     "Weighted cumulative incidence curves for live discharge, comparing high vs. low CRRT dose."),
    ("cif_death",
     "PSM Cumulative Incidence of Death",
     "Cumulative incidence of mortality in the propensity-score matched cohort."),
    ("cif_discharge",
     "PSM Cumulative Incidence of Discharge",
     "Cumulative incidence of live discharge in the propensity-score matched cohort."),
    ("dose_response_combined",
     "Dose-Response Analysis (Combined Panel)",
     "Three complementary dose-response approaches: decile, spline Cox, and generalized propensity score."),
    ("dose_response_rcs",
     "Dose-Response Curve (Restricted Cubic Spline)",
     "Covariate-adjusted hazard ratio as a smooth function of continuous CRRT dose."),
    ("gps_dose_response",
     "Generalized Propensity Score Dose-Response",
     "Doubly-robust causal dose-response curve using covariate-balancing propensity scores."),
    ("gps_loveplot",
     "GPS Covariate Balance (Love Plot)",
     "Correlation between covariates and CRRT dose before and after GPS weighting."),
    ("gps_weight_diagnostics",
     "GPS Weight Diagnostics",
     "Distribution and extremes of generalized propensity score weights."),
    ("dose_decile",
     "Dose-Decile Mortality Plot",
     "Unadjusted 30-day mortality rate by CRRT dose decile."),
    ("hist_crrt_dose",
     "CRRT Dose Distribution",
     "Histogram of initial CRRT dose (mL/kg/hr) in the analysis cohort."),
    # ── Descriptive Epidemiology (03 / 07 / 08) ──
    ("consort_diagram_straight",
     "Descriptive Cohort Exclusion Flow Diagram",
     "Patient flow into the descriptive CRRT analytic cohort."),
    ("crrt_incidence_by_context",
     "CRRT Incidence by Clinical Context",
     "Proportion of adult ICU hospitalizations receiving CRRT, overall and among "
     "ventilated / vasopressor-exposed patients."),
    ("crrt_incidence_by_icu_subtype",
     "CRRT Incidence by ICU Subtype",
     "Proportion of adult ICU hospitalizations receiving CRRT by ICU location type."),
    ("dose_distribution",
     "CRRT Dose Distribution (KDIGO 20-30 Band)",
     "Distribution of prescribed CRRT dose (mL/kg/hr) against the 20-30 band."),
    ("dose_by_ibw",
     "CRRT Dose: Actual vs Ideal Body Weight",
     "Overlaid dose histograms normalized by actual weight vs Devine ideal body "
     "weight (the ARDSnet tidal-volume basis); shows how the dose distribution and "
     "high/low classification shift when the denominator is lean rather than actual mass."),
    ("net_uf_rate_over_time",
     "Net Ultrafiltration Rate Trajectory (30 Days)",
     "Per-patient net ultrafiltration rate (mL/kg/hr) over the 30 days following "
     "initiation, against the Murugan 1.01-1.75 band. Removal ramps up from near "
     "zero at initiation toward a sustained rate."),
    ("net_uf_cumulative_over_time",
     "Cumulative Net Ultrafiltration (30 Days)",
     "Median cumulative net ultrafiltration volume (mL/kg) over the course, among "
     "encounters still on CRRT each day."),
    ("uf_mortality",
     "Crude Mortality vs Net Ultrafiltration Intensity (First 72 h)",
     "Unadjusted 30-day mortality across first-72h net ultrafiltration "
     "intensity (mL/kg/hr) bins (Murugan / 2025 literature definition)."),
    ("crrt_dose_over_time",
     "CRRT Dose Trajectory (30 Days)",
     "Per-patient CRRT dose over the 30 days following initiation."),
    ("lab_distributions_over_crrt",
     "Laboratory Trajectories (30 Days)",
     "Per-patient lab trajectories over the 30 days following CRRT initiation."),
    ("nee_over_crrt",
     "Vasopressor (NEE) Trajectory (30 Days)",
     "Per-patient norepinephrine-equivalent dose over 30 days following initiation."),
    ("imv_state_over_crrt",
     "IMV State Over Time (30 Days)",
     "Share of patients on invasive mechanical ventilation vs discharged / "
     "deceased over 30 days."),
    ("crrt_state_over_crrt",
     "CRRT State Over Time (30 Days)",
     "Share of patients on/off CRRT vs discharged / deceased over 30 days."),
    ("dose_band_distribution",
     "Dose-Band Distribution (Low-Dose Subcohort)",
     "Patients across CRRT dose bands, highlighting the very-low-dose subcohort."),
    # (MSM figure-catalog entries removed — the time-varying MSM analysis was
    # cut from the manuscript and archived, tag msm-v1.)
    ("hist_crrt_dose_0_12",
     "CRRT Dose Distribution (0\u201312h)",
     "Histogram of CRRT dose during the first 12-hour interval."),
    ("hist_crrt_dose_12_24",
     "CRRT Dose Distribution (12\u201324h)",
     "Histogram of CRRT dose during the second 12-hour interval."),
    ("hist_crrt_dose_0_24",
     "CRRT Dose Distribution (0\u201324h)",
     "Histogram of CRRT dose during the first 24-hour interval."),
    ("hist_crrt_dose_24_48",
     "CRRT Dose Distribution (24\u201348h)",
     "Histogram of CRRT dose during the second 24-hour interval."),
    ("hist_crrt_dose_combined",
     "CRRT Dose Distribution (Combined Intervals)",
     "Side-by-side histograms of CRRT dose across all time intervals."),
    ("hist_crrt_dose_t0",
     "CRRT Dose Distribution at Baseline (t=0)",
     "Histogram of CRRT dose at the initial measurement."),
]


def _figure_title_desc(filename_stem: str, site_name: str | None = None) -> tuple[str, str]:
    """Look up a human-readable title and description for a figure filename.
    Falls back to cleaned-up filename if no catalog entry matches.

    `site_name`, if given, is stripped from the front of the stem first so a
    per-site figure title never embeds the real site name (the fallback's
    `^[A-Z]{2,}_` rule only catches ALL-CAPS prefixes, so mixed-case dir names
    like `Penn_`/`Michigan_` would otherwise leak into the anonymized build)."""
    stem = filename_stem
    if site_name:
        stem = re.sub(rf"^{re.escape(site_name)}[_-]", "", stem, flags=re.IGNORECASE)
    for pattern, title, desc in _FIGURE_CATALOG:
        if pattern.lower() in stem.lower():
            return title, desc
    # Fallback: strip any residual (all-caps) site prefix and clean up
    clean = re.sub(r"^[A-Z]{2,}_(?:CRRT_\d+cutoff_)?", "", stem)
    clean = clean.replace("_", " ").strip()
    return clean, ""


def _build_site_baseline_tables_html(site_dir: Path) -> str:
    """Per-site baseline tables for the site tab, matching the manuscript
    numbering: Table 1 = full-cohort baseline by dose band (from 02
    crrt_epi/table1_crrt); Table 2 = unadjusted balance on the causal cohort
    (from 05 table2_unadjusted_balance, legacy fallback). Each is rendered only
    when its CSV is present (sites that have not re-run 02's dose-banded format
    show only Table 2)."""
    specs = [
        ("crrt_epi", "table1_crrt",
         "Table 1 — Baseline Characteristics by CRRT Dose Band",
         "Full analytic cohort at CRRT initiation, by initial dose band."),
        ("psm_iptw", BALANCE_TABLE_PATTERNS,
         "Table 2 — Unadjusted Balance Between Dose Groups (Causal Cohort)",
         "Complete-case causal cohort at CRRT initiation, by dose group "
         "(low- vs high-intensity, 30 mL/kg/hr cutoff)."),
    ]
    out = []
    for subdir, pattern, title, caption in specs:
        csv_path = _find_csv(site_dir, subdir, pattern)
        if csv_path is None:
            continue
        tbl = _render_gtsummary_table_html(_parse_gtsummary_csv(csv_path))
        if not tbl:
            continue
        out.append(
            f'<div class="section"><h2>{title}</h2>'
            f'<p class="fig-caption">{caption}</p>{tbl}</div>'
        )
    return "\n".join(out)


# Per-site tab figure layout: ordered functional sections so each figure sits
# under a clear "what is this for" header rather than alphabetical filename
# order. Each entry is (header, subdir-under-<site>/final, [filename substrings
# in display order]). Figures matched by an earlier section are not re-used;
# any PNG not claimed by a section falls into a final "Additional Figures"
# group so nothing is silently dropped.
_SITE_FIGURE_SECTIONS = [
    # ----- Descriptive epidemiology (leads, matching the epi-first dashboard) -----
    ("Descriptive Epidemiology — Cohort Flow", "crrt_epi/graphs",
     ["consort_diagram_straight"]),
    ("Descriptive Epidemiology — Incidence & Utilization", "crrt_epi/graphs",
     ["crrt_incidence_by_context", "crrt_incidence_by_icu_subtype"]),
    ("Descriptive Epidemiology — Practice & Quality", "crrt_epi/graphs",
     ["dose_distribution", "dose_by_ibw", "uf_mortality"]),
    ("Descriptive Epidemiology — Longitudinal Trajectories", "crrt_epi/graphs",
     ["crrt_dose_over_time", "net_uf_rate_over_time", "net_uf_cumulative_over_time",
      "lab_distributions_over_crrt", "nee_over_crrt",
      "imv_state_over_crrt", "crrt_state_over_crrt"]),
    ("Low-Dose Subcohort", "low_dose/graphs", ["dose_band_distribution"]),
    # ----- Causal analysis (point-treatment + continuous dose) -----
    ("Causal Analysis — Cohort Flow", "psm_iptw", ["causal_consort_diagram"]),
    ("Causal Analysis — Exposure (CRRT Dose)", "psm_iptw", ["hist_crrt_dose"]),
    ("Causal Analysis — Covariate Balance", "psm_iptw",
     ["psm_loveplot", "SL_LovePlot"]),
    ("Causal Analysis — Propensity Score Overlap & Positivity", "psm_iptw",
     ["SL_PS_overlap"]),
    ("Causal Analysis — Primary Outcome: Cumulative Incidence", "psm_iptw",
     ["cif_death", "cif_discharge", "IPTW_CIF_Death", "IPTW_CIF_Discharge"]),
    ("Causal Analysis — Continuous Dose-Response", "psm_iptw",
     ["dose_decile_plot", "dose_response_rcs_plot", "gps_dose_response_plot",
      "dose_response_combined", "gps_loveplot", "gps_weight_diagnostics"]),
    ("Causal Analysis — Subgroups", "psm_iptw", ["subgroup_forest_plot"]),
]


def _figure_block_html(p: Path, site_name: str | None = None) -> str:
    title, desc = _figure_title_desc(p.stem, site_name)
    caption = f'<p class="fig-caption">{desc}</p>' if desc else ""
    # Cap per-site figures at a smaller width than the container so the long
    # per-site tabs scroll less (figures otherwise render near full width).
    return (
        f'<div class="figure-block"><h3>{title}</h3>{caption}'
        f'<img src="{img_to_data_uri(p)}" alt="{title}" '
        f'style="max-width:1000px; width:100%; height:auto;"></div>'
    )


def _build_site_initial_settings_html(site_dir: Path) -> str:
    """Per-site CRRT settings at initiation (first 3h): median [IQR] AND % missing
    for each setting, from step-00 `{SITE}_crrt_initial_field_completeness.csv`
    (which now carries both the overall first-3h median/IQR and the completeness).
    One clean per-site view of what each CRRT setting looks like at CRRT start and
    how often it is recorded. Mode-dependent fields (dialysate/replacement/dose)
    are legitimately absent for modes that don't use them, so high missingness
    there is often expected — noted in the caption."""
    df = _load_csv(site_dir, "diagnostics", "crrt_initial_field_completeness")
    if df is None or df.empty or "pct_missing" not in df.columns:
        return ""
    have_median = "median" in df.columns
    by = {str(r["field"]).lower(): r for _, r in df.iterrows()}
    ORDER = [("blood_flow_rate", "Blood flow rate (mL/min)"),
             ("dialysate_flow_rate", "Dialysate flow rate (mL/hr)"),
             ("pre_filter_replacement_fluid_rate", "Pre-filter replacement (mL/hr)"),
             ("post_filter_replacement_fluid_rate", "Post-filter replacement (mL/hr)"),
             ("crrt_dose_ml_kg_hr", "Calculated dose (mL/kg/hr)"),
             ("ultrafiltration_out", "Net ultrafiltration (mL/hr)")]

    def _fmt_med(r, dose):
        if not have_median or pd.isna(r.get("median")):
            return "&mdash;"
        d = 1 if dose else 0
        return f'{float(r["median"]):.{d}f} [{float(r["q25"]):.{d}f}&ndash;{float(r["q75"]):.{d}f}]'

    def _miss_cell(p):
        p = pd.to_numeric(p, errors="coerce")
        if pd.isna(p):
            return '<td style="text-align:right;color:#94a3b8;">&mdash;</td>'
        col = ("#991b1b" if p >= 50 else "#9a3412" if p >= 30
               else "#854d0e" if p >= 10 else "#166534")
        return f'<td style="text-align:right;color:{col};font-weight:600;">{p:.0f}%</td>'

    rows = []
    for key, lbl in ORDER:
        r = by.get(key)
        if r is None:
            continue
        med = _fmt_med(r, dose=(key == "crrt_dose_ml_kg_hr"))
        rows.append(
            f'<tr><td style="font-size:13px;">{lbl}</td>'
            f'<td style="white-space:nowrap;">{med}</td>'
            f'{_miss_cell(r.get("pct_missing"))}</tr>')
    if not rows:
        return ""
    n_enc = ""
    try:
        n_enc = f" (n = {int(pd.to_numeric(df['n_encounters'].iloc[0])):,} encounters)"
    except Exception:
        pass
    return (
        '<div class="section"><h2>CRRT Settings at Initiation (first 3h)</h2>'
        f'<p style="font-size:13px;color:#64748b;">Median [IQR] of each setting over '
        'the first 3&nbsp;h after CRRT start (among encounters that have it), with the '
        f'% of encounters missing that setting at start{n_enc}. Mode-dependent fields '
        '(dialysate, replacement, dose) are legitimately absent for modes that don&rsquo;t '
        'use them, so high missingness there can be expected; blood flow / net UF should '
        'be recorded for most starts.</p>'
        '<table class="results-table" border="0" style="width:auto;">'
        '<thead><tr><th style="text-align:left;">Setting</th>'
        '<th>Median [IQR]</th><th>% missing</th></tr></thead>'
        '<tbody>' + "\n".join(rows) + '</tbody></table></div>'
    )


def _build_site_missingness_html(site_dir: Path) -> str:
    """Per-site data-completeness panel from step-04 `{SITE}_missingness_summary.csv`
    (missingness of the causal-df analysis variables: labs, dose, SOFA, oxygenation,
    vasopressors, CCI, ...). Rendered as a variable list sorted by %missing with a
    color-graded inline bar, so high-missingness variables stand out at a glance —
    a data-quality/ETL smell test (a variable that is far more missing than at peer
    sites, or unexpectedly missing, is worth a look). Variable names are generic
    (phosphate_0, lactate_0, ...) so this is anonymization-safe on the per-site tab."""
    df = _load_csv(site_dir, "diagnostics", "missingness_summary")
    if df is None or df.empty or "pct_missing" not in df.columns:
        return ""
    df = df.copy()
    df["pct_missing"] = pd.to_numeric(df["pct_missing"], errors="coerce")
    df = df.dropna(subset=["pct_missing"]).sort_values("pct_missing", ascending=False)

    def _band(p):  # (bar color, row emphasis)
        if p >= 50:   return "#dc2626"
        if p >= 30:   return "#ea580c"
        if p >= 10:   return "#ca8a04"
        return "#16a34a"

    rows = []
    for _, r in df.iterrows():
        p = float(r["pct_missing"]); col = _band(p)
        npres = r.get("n_present", ""); ntot = r.get("n_total", "")
        bar = (f'<div style="background:#f1f5f9;border-radius:3px;height:14px;'
               f'width:180px;overflow:hidden;display:inline-block;vertical-align:middle;">'
               f'<div style="background:{col};height:100%;width:{min(p,100):.0f}%;"></div></div>')
        rows.append(
            f'<tr>'
            f'<td style="font-family:monospace;font-size:12px;">{html.escape(str(r["variable"]), quote=False)}</td>'
            f'<td style="text-align:right;color:{col};font-weight:600;white-space:nowrap;">{p:.1f}%</td>'
            f'<td>{bar}</td>'
            f'<td style="font-size:12px;color:#64748b;white-space:nowrap;">{npres}/{ntot} present</td>'
            f'</tr>')
    n_hi = int((df["pct_missing"] >= 50).sum())
    hi_note = (f' <strong style="color:#dc2626;">{n_hi} variable(s) &ge;50% missing.</strong>'
               if n_hi else '')
    return (
        '<div class="section"><h2>Data Completeness (causal-df variables)</h2>'
        '<p style="font-size:13px;color:#64748b;">Missingness of each analysis variable '
        '(step 04), sorted highest-first. Bars: '
        '<span style="color:#16a34a;">&lt;10%</span> / '
        '<span style="color:#ca8a04;">10&ndash;30%</span> / '
        '<span style="color:#ea580c;">30&ndash;50%</span> / '
        '<span style="color:#dc2626;">&ge;50%</span>. Some missingness is real/expected '
        '(e.g. a lab drawn less often in the window); a variable far more missing than at '
        'peer sites is worth checking for an ETL/mapping issue.' + hi_note + '</p>'
        '<table class="results-table" border="0" style="width:auto;">'
        '<thead><tr><th>Variable</th><th>% missing</th><th></th><th>Present</th></tr></thead>'
        '<tbody>' + "\n".join(rows) + '</tbody></table></div>'
    )


def build_site_tab(site_dir: Path) -> str:
    """Per-site tab: baseline tables (T1 descriptive, T2 balance) followed by
    every per-site figure grouped under a clear functional header — causal
    analysis (balance / positivity / primary outcome / dose-response /
    subgroups) and descriptive epidemiology.

    Only figures named in `_SITE_FIGURE_SECTIONS` are rendered: this is a
    curated layout, not a directory dump, so stale or cut figures left on disk
    at sites that have not re-run (e.g. the removed MAP / respiratory
    trajectories, QC heatmaps) are NOT resurfaced. A genuinely new figure is
    added to the section list deliberately."""
    html_parts = []
    tables = _build_site_baseline_tables_html(site_dir)
    if tables:
        html_parts.append(tables)

    init_settings = _build_site_initial_settings_html(site_dir)
    if init_settings:
        html_parts.append(init_settings)

    missingness = _build_site_missingness_html(site_dir)
    if missingness:
        html_parts.append(missingness)

    claimed: set[Path] = set()
    for header, subdir, patterns in _SITE_FIGURE_SECTIONS:
        d = _site_final_dir(site_dir) / subdir
        if not d.exists():
            continue
        pngs = sorted(d.glob("*.png"))
        picked: list[Path] = []
        for pat in patterns:
            for p in pngs:
                if p not in claimed and pat.lower() in p.name.lower():
                    picked.append(p)
                    claimed.add(p)
        if picked:
            html_parts.append(
                f'<div class="section"><h2>{header}</h2>'
                + "\n".join(_figure_block_html(p, site_dir.name) for p in picked)
                + '</div>'
            )

    return "\n".join(html_parts) if html_parts else _missing("No results for this site.")


# ── Completeness Checklist ─────────────────────────────────────────────────

def _check_file(site_dir: Path, subdir: str, pattern, ext: str) -> bool:
    """Check if a file matching pattern exists in site_dir/final/subdir/.
    `pattern` may be a single substring or a sequence of substrings (present if
    any one matches)."""
    d = _site_final_dir(site_dir) / subdir
    if not d.exists():
        return False
    patterns = [pattern] if isinstance(pattern, str) else list(pattern)
    glob_ext = f"*.{ext}"
    for f in d.glob(glob_ext):
        if any(pat.lower() in f.name.lower() for pat in patterns):
            return True
    return False


def _expected_matchers() -> list[tuple[str, list[str], str]]:
    """Flatten EXPECTED_FILES into (subdir, patterns, ext) matcher tuples used to
    DEMOTE an unconsumed file from stale to expected-but-unrendered.

    For a multi-pattern entry (e.g. BALANCE_TABLE_PATTERNS = new-name, legacy),
    only the PREFERRED (first) pattern is recognized. A not-yet-rerun site's
    legacy file is always consumed by `_find_csv`'s fallback (so it never reaches
    classification); the only way a legacy name arrives here unconsumed is as a
    leftover beside the consumed new file — a genuine shadowing orphan that must
    stay in the stale bucket, not be demoted."""
    out = []
    for files in EXPECTED_FILES.values():
        for subdir, pattern, _label, ext, _script in files:
            preferred = pattern if isinstance(pattern, str) else pattern[0]
            out.append((subdir, [preferred.lower()], ext.lower()))
    return out


def _matches_expected(rel_subdir: str, name: str, ext: str,
                      matchers: list[tuple[str, list[str], str]]) -> bool:
    """True if a file (subdir relative to the site's final dir, filename, ext)
    matches any EXPECTED_FILES entry — i.e. it is a recognized pipeline output,
    just not rendered by this build."""
    nl = name.lower()
    for m_subdir, m_pats, m_ext in matchers:
        if m_subdir == rel_subdir and m_ext == ext and any(p in nl for p in m_pats):
            return True
    return False


def _classify_uncontributed(rel_subdir: str, name: str, ext: str,
                            matchers: list[tuple[str, list[str], str]]) -> str:
    """Classify a file present on disk but not consumed by this build into one
    of four buckets (checked in precedence order):
      - "qc"       : anything under diagnostics/ — internal QC, never displayed.
      - "expected" : a recognized EXPECTED_FILES output a tab just didn't render.
      - "other"    : a recognized non-dashboard output (manuscript-08 input or a
                     deep causal diagnostic) — contributes elsewhere, not here.
      - "stale"    : matches none of the above — a genuine renamed/cut/stale
                     orphan that can shadow current output. The only red flag.
    """
    if rel_subdir == "diagnostics" or rel_subdir.startswith("diagnostics/"):
        return "qc"
    if _matches_expected(rel_subdir, name, ext, matchers):
        return "expected"
    nl = name.lower()
    if any(s in nl for s in NON_DASHBOARD_OUTPUTS):
        return "other"
    return "stale"


def _build_not_contributing_html(sites, labels) -> str:
    """Flag files present on disk but not consumed by this dashboard build.

    Four buckets (see `_classify_uncontributed`). Only STALE is a genuine
    concern — a renamed/cut/stale file that can shadow current output (recall the
    `_find_csv` shorter-name preference). QC, recognized non-dashboard outputs,
    and expected-but-unrendered files are benign and shown collapsed.

    Only `.csv` and `.png` are scanned — the two extensions the dashboard
    consumes; sibling `.pdf`/`.parquet` print/data artifacts are ignored.
    """
    matchers = _expected_matchers()
    BUCKETS = (
        ("stale", "Stale / orphan", "#991b1b",
         "not consumed, and not a recognized output &mdash; likely renamed, cut, "
         "or stale; can shadow current results"),
        ("qc", "QC / not displayed", "#64748b",
         "internal QC under diagnostics/, intentionally not shown"),
        ("other", "Non-dashboard output", "#64748b",
         "recognized manuscript-08 input or deep causal diagnostic; contributes "
         "to a pooled figure or provenance, not this dashboard"),
        ("expected", "Expected but not rendered", "#94a3b8",
         "a recognized output a tab simply did not surface this build "
         "(section skipped or k&lt;threshold)"),
    )
    cards = []
    total_stale = 0
    for sd in sites:
        sid = sd.name
        lbl = labels.get(sid, sid)
        color = get_color(sid)
        final = _site_final_dir(sd)
        found: dict[str, list[str]] = {k: [] for k, *_ in BUCKETS}
        if final.is_dir():
            files = sorted(final.rglob("*.csv")) + sorted(final.rglob("*.png"))
            for f in files:
                if f.resolve() in CONSUMED:
                    continue
                rel = f.relative_to(final)
                rel_subdir = str(rel.parent) if str(rel.parent) != "." else ""
                ext = f.suffix.lstrip(".").lower()
                found[_classify_uncontributed(rel_subdir, f.name, ext, matchers)].append(str(rel))
        total_stale += len(found["stale"])

        def _list(items, col):
            return "".join(
                f'<div style="padding:2px 0;font-family:monospace;font-size:12px;'
                f'color:{col};">{html.escape(p, quote=False)}</div>'
                for p in sorted(items)
            )

        if not any(found.values()):
            body = ('<div style="color:#16a34a;font-size:13px;padding:4px 0;">'
                    '&#10003; Every file on disk contributed to a figure or table.</div>')
        else:
            body = ""
            for key, title, col, desc in BUCKETS:
                items = found[key]
                if not items:
                    continue
                warn = "&#9888; " if key == "stale" else ""
                body += (
                    f'<details style="margin:6px 0;"{" open" if key == "stale" else ""}>'
                    f'<summary style="cursor:pointer;font-weight:600;color:{col};'
                    f'font-size:13px;">{warn}{title} ({len(items)})'
                    f'<span style="font-weight:400;color:#94a3b8;"> &mdash; {desc}</span>'
                    f'</summary>{_list(items, col)}</details>'
                )

        n_stale = len(found["stale"])
        n_other = sum(len(v) for k, v in found.items() if k != "stale")
        cards.append(
            f'<details style="border:1px solid #e2e8f0;border-left:4px solid {color};'
            f'border-radius:6px;padding:10px 16px;margin-bottom:10px;"'
            f'{" open" if n_stale else ""}>'
            f'<summary style="cursor:pointer;font-weight:600;color:#334155;">'
            f'{html.escape(lbl, quote=False)} '
            f'<span style="color:{"#991b1b" if n_stale else "#16a34a"};font-weight:600;">'
            f'{n_stale} stale</span>'
            f'<span style="color:#94a3b8;font-weight:400;"> '
            f'&middot; {n_other} other non-contributing</span></summary>'
            f'<div style="padding:6px 0 2px;">{body}</div></details>'
        )

    banner = (
        f'<p style="color:#991b1b;font-weight:600;font-size:13px;">'
        f'&#9888; {total_stale} stale / orphan file(s) across all sites &mdash; '
        f'review before deletion.</p>' if total_stale else
        '<p style="color:#16a34a;font-weight:600;font-size:13px;">'
        '&#10003; No stale / orphan files: every unrecognized csv/png on disk is '
        'accounted for (QC, manuscript input, or an unrendered expected output).</p>'
    )
    return (
        '<div class="section"><h2>Present but Not Contributing</h2>'
        '<p style="color:#64748b;font-size:13px;">Files on disk that this dashboard '
        'build never opened, grouped by why. Only <strong style="color:#991b1b;">'
        'stale / orphan</strong> files warrant review &mdash; they are unrecognized '
        'outputs (renamed, cut) that can shadow current results. QC diagnostics, '
        'recognized manuscript-08 inputs, and expected-but-unrendered outputs are '
        'benign and shown collapsed.</p>' + banner + "".join(cards) + '</div>'
    )


# ── Pooling-eligibility ledger (from site_validation.py) ────────────────────
# Surfaces the three-gate validation verdicts (output/multi_site/site_validation.json)
# in the Completeness tab so it is transparent WHICH sites carry issues that would
# make them ineligible for WHICH pooled analysis, and WHY. This is a documentation
# layer only — it does NOT filter any pooled figure/table. Inclusion decisions are
# a separate task. Run `code/site_validation.py` to (re)generate the JSON.

# machine domain key → human label (order = display order)
_ELIG_DOMAINS = [
    ("descriptive_baseline",  "Descriptive baseline (Table 1)"),
    ("descriptive_mortality", "Crude mortality"),
    ("dose_distribution",     "Dose distribution"),
    ("modality_mix",          "Modality mix"),
    ("SMR",                   "SMR (risk-standardized mortality)"),
    ("causal_dose_contrast",  "Causal dose contrast (≥30 vs <30)"),
    ("causal_outcome",        "Causal outcome (mortality)"),
]
_POOL_COLOR = {"CLEARED": ("#16a34a", "#dcfce7"),
               "CONDITIONAL": ("#ca8a04", "#fef9c3"),
               "EXCLUDED": ("#dc2626", "#fee2e2")}


def _load_site_validation():
    """Load output/multi_site/site_validation.json, or None if absent."""
    p = ROOT / "site_validation.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (ValueError, OSError):
        return None


def _blocking_reasons(site_rec: dict) -> dict:
    """Map domain → list of human reasons for a site, from its gate checks.
    Gate 1/2 FAIL blocks every domain (hard exclusion); Gate 3 HOLD blocks only
    the domains each check scopes to."""
    reasons: dict = {}
    checks = site_rec.get("checks", {})
    # Gate 1/2 hard failures apply to all domains
    for gk, glabel in (("gate1", "Integrity"), ("gate2", "Consistency")):
        for c in checks.get(gk, []):
            if c.get("status") == "FAIL":
                for dk, _ in _ELIG_DOMAINS:
                    reasons.setdefault(dk, []).append(f"Gate 2 {glabel}: {c['msg']}"
                                                      if gk == "gate2" else f"Gate 1: {c['msg']}")
    # Gate 3 scoped holds
    for c in checks.get("gate3", []):
        if c.get("status") == "HOLD":
            for dk in c.get("blocks", []):
                reasons.setdefault(dk, []).append(f"Gate 3 {c['name']}: {c['msg']}")
    return reasons


def build_pooling_eligibility_html(sites, labels) -> str:
    """Ledger of per-site validation verdicts + which pooled analyses each site's
    issues would affect and why. Documentation only; filters nothing."""
    val = _load_site_validation()
    if val is None:
        return (
            '<div class="section"><h2>Pooling Eligibility &amp; Site Issues</h2>'
            '<p class="missing">No <code>site_validation.json</code> found. Run '
            '<code>uv run python code/site_validation.py</code> from the repo root '
            'to generate per-site gate verdicts, then rebuild.</p></div>'
        )
    recs = val.get("sites", {})
    site_ids = [sd.name for sd in sites if sd.name in recs]

    # ── (1) Gate scorecard: one row per site ──
    hdr = ("<th>Site</th><th>Gate 1<br>Integrity</th><th>Gate 2<br>Consistency</th>"
           "<th>Gate 3<br>Plausibility</th><th>Pooling status</th><th>Headline issue</th>")
    rows = []
    _sev = {"CLEARED": 0, "CONDITIONAL": 1, "EXCLUDED": 2}
    for sid in sorted(site_ids, key=lambda s: (-_sev.get(recs[s].get("pooling"), 0), s)):
        r = recs[sid]
        pool = r.get("pooling", "?")
        fg, bg = _POOL_COLOR.get(pool, ("#334155", "#f1f5f9"))
        def _cell(st):
            c = {"PASS": "#16a34a", "HOLD": "#ca8a04", "FAIL": "#dc2626"}.get(st, "#64748b")
            return f'<td style="text-align:center;color:{c};font-weight:600;">{st}</td>'
        # headline = first FAIL (gate2/1) else first HOLD msg
        headline = "&mdash;"
        ch = r.get("checks", {})
        fails = [c for gk in ("gate2", "gate1") for c in ch.get(gk, []) if c.get("status") == "FAIL"]
        holds = [c for c in ch.get("gate3", []) if c.get("status") == "HOLD"]
        if fails:
            headline = html.escape(fails[0]["msg"], quote=False)
        elif holds:
            headline = html.escape("; ".join(c["name"] for c in holds), quote=False)
        rows.append(
            f'<tr style="border-left:4px solid {get_color(sid)};">'
            f'<td><strong>{html.escape(labels.get(sid, sid), quote=False)}</strong></td>'
            f'{_cell(r.get("gate1"))}{_cell(r.get("gate2"))}{_cell(r.get("gate3"))}'
            f'<td style="text-align:center;"><span style="background:{bg};color:{fg};'
            f'padding:2px 10px;border-radius:12px;font-weight:600;font-size:12px;">{pool}</span></td>'
            f'<td style="font-size:12px;color:#475569;">{headline}</td></tr>'
        )
    scorecard = (
        '<table class="results-table" border="0"><thead><tr>'
        + hdr + '</tr></thead><tbody>' + "\n".join(rows) + '</tbody></table>'
    )

    # ── (2) Domain × site matrix: ✓ eligible / ⚠ flagged / ✗ excluded ──
    mhdr = "<th>Analysis domain</th>" + "".join(
        f'<th style="border-bottom:3px solid {get_color(s)};font-size:11px;">'
        f'{html.escape(labels.get(s, s), quote=False)}</th>' for s in site_ids)
    mrows = []
    for dk, dlabel in _ELIG_DOMAINS:
        tds = f'<td style="font-size:12px;"><strong>{dlabel}</strong></td>'
        n_ok = 0
        for sid in site_ids:
            r = recs[sid]
            blocked = dk in r.get("blocked_domains", [])
            if not blocked:
                n_ok += 1; mark, col = "&#10003;", "#16a34a"
            elif r.get("pooling") == "EXCLUDED":
                mark, col = "&#10007;", "#dc2626"
            else:
                mark, col = "&#9888;", "#ca8a04"
            tds += f'<td style="text-align:center;color:{col};font-weight:700;">{mark}</td>'
        tds += f'<td style="text-align:center;font-weight:600;color:#334155;">{n_ok}/{len(site_ids)}</td>'
        mrows.append(f"<tr>{tds}</tr>")
    matrix = (
        '<table class="results-table" border="0"><thead><tr>'
        + mhdr + '<th>Eligible</th></tr></thead><tbody>' + "\n".join(mrows)
        + '</tbody></table>'
        '<p style="font-size:12px;color:#64748b;">&#10003; eligible &middot; '
        '<span style="color:#ca8a04;">&#9888;</span> flagged (conditional &mdash; pending '
        'adjudication) &middot; <span style="color:#dc2626;">&#10007;</span> excluded '
        '(Gate&nbsp;1/2 failure). Counts document which sites carry issues per analysis; '
        'they do <strong>not</strong> filter any figure &mdash; inclusion is a separate decision.</p>'
    )

    # ── (3) Per-domain "why" breakdown ──
    why_blocks = []
    for dk, dlabel in _ELIG_DOMAINS:
        entries = []
        for sid in site_ids:
            r = recs[sid]
            if dk in r.get("blocked_domains", []):
                rs = _blocking_reasons(r).get(dk, [])
                reason = "; ".join(rs) if rs else "(see gate detail)"
                col = "#dc2626" if r.get("pooling") == "EXCLUDED" else "#ca8a04"
                entries.append(
                    f'<div style="padding:2px 0;font-size:12px;">'
                    f'<span style="color:{col};font-weight:600;">'
                    f'{html.escape(labels.get(sid, sid), quote=False)}</span>'
                    f' <span style="color:#475569;">&mdash; {html.escape(reason, quote=False)}</span></div>')
        if entries:
            why_blocks.append(
                f'<details style="margin-bottom:6px;"><summary style="cursor:pointer;'
                f'font-weight:600;color:#334155;font-size:13px;">{dlabel} '
                f'<span style="color:#94a3b8;font-weight:400;">({len(entries)} flagged)</span>'
                f'</summary><div style="padding:4px 0 6px 20px;">{"".join(entries)}</div></details>')
    why = ("".join(why_blocks) if why_blocks
           else '<p style="color:#16a34a;">No site issues recorded &mdash; all sites clear every domain.</p>')

    gen = val.get("generated", "")
    return (
        '<div class="section"><h2>Pooling Eligibility &amp; Site Issues</h2>'
        f'<p style="color:#64748b;font-size:13px;">Three-gate validation verdicts from '
        f'<code>site_validation.py</code> (generated {html.escape(str(gen)[:19], quote=False)}). '
        'This documents which sites carry data-quality issues and which pooled analyses each issue '
        'affects. <strong>It does not drop any site from any figure</strong> &mdash; whether to '
        'include a flagged site is a separate, explicit decision (tracked in '
        '<code>results_review_tracker.md</code> &sect;5b).</p>'
        '<h3 style="margin-top:16px;">Gate scorecard</h3>' + scorecard
        + '<h3 style="margin-top:20px;">Eligibility by analysis domain</h3>' + matrix
        + '<h3 style="margin-top:20px;">Why sites are flagged (per domain)</h3>' + why
        + '</div>'
    )


def build_crrt_initial_completeness_html(sites, labels) -> str:
    """Cross-site comparison of PATIENT-LEVEL first-3h CRRT field completeness
    (step-00 `{SITE}_crrt_initial_field_completeness.csv`): for each CRRT field,
    the % of encounters with NO value recorded in the first 3h after CRRT start —
    a clean, un-confounded ETL/mapping signal (unlike the raw record-level
    crrt_column_missingness, which reflects long-format charting structure). A
    fields × sites heatmap so a field far more missing at one site than its peers
    stands out. Completeness tab (internal QC; omitted from anon). Transitional:
    only sites that have re-run step 00 under this change appear.

    Split by interpretation: mode / blood flow / net UF should be recorded for
    ~all CRRT starts, so high missingness there is a real ETL flag; dialysate,
    replacement, and dose are MODE-dependent (dose only for cvvh/cvvhd/cvvhdf;
    dialysate absent for CVVH, replacement for CVVHD), so read those with the
    modality mix."""
    CORE = [("crrt_mode_category", "CRRT mode"),
            ("blood_flow_rate", "Blood flow rate"),
            ("ultrafiltration_out", "Net ultrafiltration")]
    MODE_DEP = [("dialysate_flow_rate", "Dialysate flow rate"),
                ("pre_filter_replacement_fluid_rate", "Pre-filter replacement"),
                ("post_filter_replacement_fluid_rate", "Post-filter replacement"),
                ("crrt_dose_ml_kg_hr", "Calculated dose")]
    site_pct, present = {}, []
    for sd in sites:
        df = _load_csv(sd, "diagnostics", "crrt_initial_field_completeness")
        if df is None or df.empty or "pct_missing" not in df.columns:
            continue
        site_pct[sd.name] = {str(k).lower(): pd.to_numeric(v, errors="coerce")
                             for k, v in zip(df["field"], df["pct_missing"])}
        present.append(sd)
    if not present:
        return ('<div class="section"><h2>CRRT Field Completeness at Start (cross-site)</h2>'
                '<p class="missing">No site has re-run step 00 under the first-3h field-completeness '
                'change yet. Populates on re-run '
                '(<code>crrt_initial_field_completeness.csv</code>).</p></div>')

    def _cell(pct):
        if pct is None or (isinstance(pct, float) and pct != pct):
            return '<td style="text-align:center;color:#94a3b8;background:#f8fafc;">&mdash;</td>'
        if pct >= 50:   bg, fg = "#fee2e2", "#991b1b"
        elif pct >= 30: bg, fg = "#ffedd5", "#9a3412"
        elif pct >= 10: bg, fg = "#fef9c3", "#854d0e"
        else:           bg, fg = "#dcfce7", "#166534"
        return f'<td style="text-align:center;background:{bg};color:{fg};font-weight:600;">{pct:.0f}%</td>'

    hdr = "<th style='text-align:left;'>CRRT field</th>" + "".join(
        f'<th style="border-bottom:3px solid {get_color(sd.name)};font-size:11px;">'
        f'{html.escape(labels.get(sd.name, sd.name), quote=False)}</th>' for sd in present)

    def _rows(fields):
        out = []
        for key, lbl in fields:
            tds = f'<td style="font-size:12px;">{lbl}</td>'
            for sd in present:
                v = site_pct.get(sd.name, {}).get(key)
                tds += _cell(None if (v is None or (isinstance(v, float) and v != v)) else float(v))
            out.append(f"<tr>{tds}</tr>")
        return "\n".join(out)

    return (
        '<div class="section"><h2>CRRT Field Completeness at Start (cross-site)</h2>'
        '<p style="color:#64748b;font-size:13px;">% of encounters with <strong>no value '
        'recorded in the first 3&nbsp;h after CRRT start</strong> for each field '
        '(patient-level, step 00 &mdash; the same 3h-median window as the dose exposure). '
        'Scan a row across sites: a field far more missing at one site than peers is a '
        'likely <strong>ETL / mapping issue</strong>. Cells: '
        '<span style="color:#166534;">&lt;10%</span> / '
        '<span style="color:#854d0e;">10&ndash;30%</span> / '
        '<span style="color:#9a3412;">30&ndash;50%</span> / '
        '<span style="color:#991b1b;">&ge;50%</span>.</p>'
        '<h3 style="font-size:14px;margin:10px 0 4px;">Core fields '
        '<span style="font-weight:400;color:#64748b;">&mdash; recorded for ~all CRRT starts; '
        'high missingness = real ETL flag</span></h3>'
        '<table class="results-table" border="0" style="width:auto;"><thead><tr>'
        + hdr + '</tr></thead><tbody>' + _rows(CORE) + '</tbody></table>'
        '<h3 style="font-size:14px;margin:16px 0 4px;">Mode-dependent fields '
        '<span style="font-weight:400;color:#64748b;">&mdash; legitimately absent for modes that '
        'do not use them (dialysate for CVVH, replacement for CVVHD, dose for SCUF/AVVH); read '
        'with the modality mix</span></h3>'
        '<table class="results-table" border="0" style="width:auto;"><thead><tr>'
        + hdr + '</tr></thead><tbody>' + _rows(MODE_DEP) + '</tbody></table></div>'
    )


def build_modality_dq_html(sites, labels) -> str:
    """Modality data-quality diagnostic (step-00 outputs): (1) which non-CRRT modes
    are DROPPED at the valid-mode filter per site (excluded_modes_breakdown), and
    (2) the fluid-pattern composition of the KEPT null-mode 'Unknown' rows
    (unknown_modality_composition) — inferable modes (mapping gap) vs genuinely
    uncharacterized. Documents the two distinct modality problems for site ETL
    follow-up. Transitional: only sites that have re-run step 00 appear. Lives in
    the Completeness tab (internal QC; omitted from the anonymized build)."""
    excl_rows, null_rows, ran_sites = [], [], set()
    for sd in sites:
        lbl = labels.get(sd.name, sd.name)
        ex = _load_csv(sd, "diagnostics", "excluded_modes_breakdown")
        if ex is not None:
            ran_sites.add(sd.name)
        if ex is not None and not ex.empty:
            for _, r in ex.iterrows():
                excl_rows.append({
                    "site_id": sd.name, "Site": lbl,
                    "Mode category": r.get("crrt_mode_category", ""),
                    "Raw mode name": r.get("crrt_mode_name", "—"),
                    "Records": r.get("n_records", ""),
                    "Encounters": r.get("n_encounter_blocks", "")})
        nl = _load_csv(sd, "diagnostics", "unknown_modality_composition")
        if nl is not None and not nl.empty:
            for _, r in nl.iterrows():
                null_rows.append({
                    "site_id": sd.name, "Site": lbl,
                    "Inferred mode": r.get("inferred_mode", ""),
                    "Raw mode name": r.get("crrt_mode_name", "—"),
                    "Records": r.get("n_records", ""),
                    "Encounters": r.get("n_encounter_blocks", "")})
    if not excl_rows and not null_rows:
        if ran_sites:
            return (
                '<div class="section"><h2>Modality Data Quality</h2>'
                f'<p style="color:#16a34a;">&#10003; {len(ran_sites)} site(s) have re-run '
                'step 00 under the modality diagnostic and show <strong>0 dropped modes '
                'and 0 Unknown (null) rows</strong> &mdash; clean modality mapping. The '
                'affected sites (UCSF/UMN/OHSU/NU/Emory) populate this section when they '
                're-run.</p></div>'
            )
        return (
            '<div class="section"><h2>Modality Data Quality</h2>'
            '<p class="missing">No site has re-run step 00 under the modality-diagnostic '
            'change yet. Populates on re-run (<code>excluded_modes_breakdown</code> / '
            '<code>unknown_modality_composition</code>).</p></div>'
        )
    excl_html = (_styled_table(pd.DataFrame(excl_rows),
                 display_cols=["Site", "Mode category", "Raw mode name", "Records", "Encounters"])
                 if excl_rows else "<p style='color:#16a34a;'>No modes dropped at any re-run site.</p>")
    null_html = (_styled_table(pd.DataFrame(null_rows),
                 display_cols=["Site", "Inferred mode", "Raw mode name", "Records", "Encounters"])
                 if null_rows else "<p style='color:#16a34a;'>No null-mode (Unknown) rows at any re-run site.</p>")
    return (
        '<div class="section"><h2>Modality Data Quality</h2>'
        '<p style="color:#64748b;font-size:13px;">Two distinct modality issues from '
        'step 00. Transitional &mdash; only sites that have re-run appear.</p>'
        '<h3>1. Dropped at the valid-mode filter '
        '(<code>excluded_non_crrt_modes</code>)</h3>'
        '<p style="font-size:13px;">Rows removed because the mode is <em>present but '
        'not a recognized continuous CRRT mode</em> (IHD, PIRRT, or an unmapped '
        'variant that should map to a valid mode). These episodes are excluded from '
        'the cohort. Per site: the mode category, the raw <code>crrt_mode_name</code> '
        'text, and record/encounter counts.</p>' + excl_html +
        '<h3 style="margin-top:20px;">2. Composition of KEPT null-mode '
        '(&ldquo;Unknown&rdquo;) rows</h3>'
        '<p style="font-size:13px;">Null <code>crrt_mode_category</code> is <em>not</em> '
        'dropped &mdash; it survives as Unknown. Mode is <strong>inferred from the fluid '
        'pattern</strong> (dialysate&rarr;CVVHD, replacement&rarr;CVVH, both&rarr;CVVHDF, '
        'UF-only&rarr;SCUF): inferable rows are a <strong>labeling/ETL gap</strong> (the '
        'therapy is fully characterized, only the category is blank &mdash; fixable at the '
        'site), while &ldquo;uncharacterized&rdquo; rows have no flows. Record-level over '
        'all CRRT records (pre-cohort-restriction). The raw <code>crrt_mode_name</code> is '
        'the text the site&rsquo;s ETL did not map.</p>' + null_html + '</div>'
    )


def build_completeness_checklist(sites, labels) -> str:
    """Build a completeness checklist showing expected vs present files per site."""

    # Gather results: {site_name: {section: [(label, found, script)]}}
    results = {}
    for sd in sites:
        sid = sd.name
        results[sid] = {}
        for section, files in EXPECTED_FILES.items():
            checks = []
            for subdir, pattern, label, ext, script in files:
                found = _check_file(sd, subdir, pattern, ext)
                checks.append((label, found, script))
            results[sid][section] = checks

    sections = list(EXPECTED_FILES.keys())

    # ── Summary Matrix ──
    # Header row: Section | Site1 | Site2 | ...
    site_ids = [sd.name for sd in sites]
    header = "<th>Section</th>" + "".join(
        f'<th style="border-bottom:3px solid {get_color(sid)};">{labels.get(sid, sid)}</th>'
        for sid in site_ids
    )
    rows = []
    # Per-section row
    for section in sections:
        tds = f"<td><strong>{section}</strong></td>"
        for sid in site_ids:
            checks = results[sid][section]
            total = len(checks)
            found = sum(1 for _, ok, _ in checks if ok)
            pct = round(100 * found / total) if total else 0
            if pct >= 90:
                bg, fg = "#dcfce7", "#166534"
            elif pct >= 50:
                bg, fg = "#fef9c3", "#854d0e"
            else:
                bg, fg = "#fee2e2", "#991b1b"
            tds += (
                f'<td style="background:{bg};color:{fg};text-align:center;'
                f'font-weight:600;border-radius:4px;">'
                f'{found}/{total} ({pct}%)</td>'
            )
        rows.append(f"<tr>{tds}</tr>")

    # Overall row
    overall_tds = "<td><strong>Overall</strong></td>"
    for sid in site_ids:
        total = sum(len(checks) for checks in results[sid].values())
        found = sum(
            sum(1 for _, ok, _ in checks if ok)
            for checks in results[sid].values()
        )
        pct = round(100 * found / total) if total else 0
        if pct >= 90:
            bg, fg = "#dcfce7", "#166534"
        elif pct >= 50:
            bg, fg = "#fef9c3", "#854d0e"
        else:
            bg, fg = "#fee2e2", "#991b1b"
        overall_tds += (
            f'<td style="background:{bg};color:{fg};text-align:center;'
            f'font-weight:700;border-radius:4px;font-size:15px;">'
            f'{found}/{total} ({pct}%)</td>'
        )
    rows.append(f'<tr style="border-top:2px solid #cbd5e1;">{overall_tds}</tr>')

    matrix_html = (
        '<div class="section"><h2>Output Completeness Summary</h2>'
        '<p style="color:#64748b;font-size:13px;">Check that each site has uploaded '
        'all required output files before running the combined dashboard.</p>'
        '<table class="results-table" border="0">'
        f'<thead><tr>{header}</tr></thead>'
        '<tbody>' + "\n".join(rows) + '</tbody></table></div>'
    )

    # ── Present-but-not-contributing flag (two-tier) ──
    # The reverse of the checklist above: files that ARE on disk but that this
    # dashboard build never opened. Because this tab is built last, the shared
    # `CONSUMED` set (populated by report_core's `_find_csv` / `img_to_data_uri`
    # as every other tab rendered) is a complete record of what contributed.
    notcontrib_html = _build_not_contributing_html(sites, labels)

    # ── Per-Site Detail Cards ──
    cards_html = ""
    for sd in sites:
        sid = sd.name
        lbl = labels.get(sid, sid)
        color = get_color(sid)

        total_all = sum(len(checks) for checks in results[sid].values())
        found_all = sum(
            sum(1 for _, ok, _ in checks if ok)
            for checks in results[sid].values()
        )
        pct_all = round(100 * found_all / total_all) if total_all else 0
        if pct_all >= 90:
            badge_bg, badge_fg = "#dcfce7", "#166534"
        elif pct_all >= 50:
            badge_bg, badge_fg = "#fef9c3", "#854d0e"
        else:
            badge_bg, badge_fg = "#fee2e2", "#991b1b"

        sections_html = ""
        for section in sections:
            checks = results[sid][section]
            s_total = len(checks)
            s_found = sum(1 for _, ok, _ in checks if ok)

            items = ""
            for label, found, script in checks:
                if found:
                    icon = '<span style="color:#16a34a;">&#10003;</span>'
                    text = f'<span>{label}</span>'
                else:
                    icon = '<span style="color:#dc2626;">&#10007;</span>'
                    text = (
                        f'<span style="color:#dc2626;">{label}</span>'
                        f' <span style="color:#94a3b8;font-size:12px;">'
                        f'(run script {script})</span>'
                    )
                items += f'<div style="display:flex;gap:8px;align-items:center;padding:3px 0;">{icon}{text}</div>'

            s_pct = round(100 * s_found / s_total) if s_total else 0
            if s_pct >= 90:
                s_color = "#16a34a"
            elif s_pct >= 50:
                s_color = "#ca8a04"
            else:
                s_color = "#dc2626"

            sections_html += (
                f'<details style="margin-bottom:8px;">'
                f'<summary style="cursor:pointer;padding:6px 0;font-weight:600;'
                f'color:#334155;font-size:13px;">'
                f'{section} '
                f'<span style="color:{s_color};font-weight:500;">'
                f'{s_found}/{s_total}</span></summary>'
                f'<div style="padding:4px 0 8px 24px;font-size:13px;">{items}</div>'
                f'</details>'
            )

        cards_html += (
            f'<div style="border:1px solid #e2e8f0;border-left:4px solid {color};'
            f'border-radius:6px;padding:16px 20px;margin-bottom:16px;">'
            f'<div style="display:flex;align-items:center;gap:12px;margin-bottom:12px;">'
            f'<h3 style="margin:0;">{lbl}</h3>'
            f'<span style="background:{badge_bg};color:{badge_fg};padding:2px 10px;'
            f'border-radius:12px;font-size:13px;font-weight:600;">'
            f'{found_all}/{total_all} ({pct_all}%)</span></div>'
            f'{sections_html}</div>'
        )

    detail_html = f'<div class="section"><h2>Per-Site Details</h2>{cards_html}</div>'

    eligibility_html = build_pooling_eligibility_html(sites, labels)
    initcomp_html = build_crrt_initial_completeness_html(sites, labels)
    modality_dq_html = build_modality_dq_html(sites, labels)

    return (matrix_html + eligibility_html + initcomp_html
            + modality_dq_html + notcontrib_html + detail_html)


# ── Dashboard Assembly ──────────────────────────────────────────────────────

def _slugify(s: str) -> str:
    """DOM-safe slug from a label, e.g. 'Site 8' → 'site8'. Used for anonymized
    per-site tab ids so real site names never appear in DOM ids / switchTab / JS."""
    return re.sub(r"[^a-z0-9]+", "", s.lower()) or "site"


def _tab_btn(tab_id: str, label: str, active: bool = False) -> str:
    cls = "tab-btn active" if active else "tab-btn"
    return f'<button class="{cls}" onclick="switchTab(\'{tab_id}\')" id="btn-{tab_id}">{label}</button>'


def _tab_group(group_label: str, buttons: list[str]) -> str:
    """Wrap a set of tab buttons in a labeled group."""
    btns = "\n".join(buttons)
    return (
        f'<div class="tab-group">'
        f'<span class="tab-group-label">{group_label}</span>'
        f'<div class="tab-group-buttons">{btns}</div>'
        f'</div>'
    )


# ── Between-site heterogeneity of the descriptive epi metrics ────────────────
# Federated comparison: each site ships only summaries (counts, means, SDs); we
# pool + test for between-site heterogeneity here. Proportions are pooled on the
# logit scale (DerSimonian-Laird RE) with a k×2 chi-square test of homogeneity as
# the direct between-site difference test; continuous metrics are pooled by
# inverse-variance RE on per-site means (SE = SD/sqrt(n)) with Cochran's Q. Q / I²
# quantify the magnitude of heterogeneity either way. This is the single source of
# truth for the computation: the dashboard renders it inline (build_crrt_epidemiology)
# and 08's export_epi_heterogeneity_csv bridges compute_epi_heterogeneity() to
# write pooled_epi_heterogeneity.csv.









# Column order for pooled_epi_heterogeneity.csv (08 writer) + dashboard render.


def _smr_forest_fig(per_site, pooled):
    """Forest of per-site SMRs (O/E) with exact-Poisson CIs + pooled, ref line at 1."""
    import matplotlib.pyplot as plt
    rows = list(per_site)
    show_pool = bool(pooled) and pooled["k"] > 1
    labels_y = [p["label"] for p in rows] + (["Pooled"] if show_pool else [])
    smr = [p["smr"] for p in rows] + ([pooled["smr"]] if show_pool else [])
    lo = [p["lo"] for p in rows] + ([pooled["lo"]] if show_pool else [])
    hi = [p["hi"] for p in rows] + ([pooled["hi"]] if show_pool else [])
    colors = [get_color(p["site_id"]) for p in rows] + (["#1e417c"] if show_pool else [])
    m = len(labels_y)
    yy = list(range(m))[::-1]
    fig, ax = plt.subplots(figsize=(8, max(2.2, 0.62 * m + 1.2)))
    for i in range(m):
        ax.errorbar(smr[i], yy[i], xerr=[[smr[i] - lo[i]], [hi[i] - smr[i]]], fmt="o",
                    color=colors[i], ecolor=colors[i], capsize=4,
                    markersize=10 if labels_y[i] == "Pooled" else 7,
                    elinewidth=2 if labels_y[i] == "Pooled" else 1.3)
    ax.axvline(1.0, ls="--", color="black", lw=1)
    xmax = max(hi) if hi else 1.2
    for i in range(m):
        ax.text(xmax * 1.02, yy[i], f"{smr[i]:.2f} ({lo[i]:.2f}–{hi[i]:.2f})",
                va="center", fontsize=9)
    ax.set_yticks(yy); ax.set_yticklabels(labels_y)
    ax.set_xlim(min(lo) - 0.03 if lo else 0.9, xmax * 1.20)
    ax.set_xlabel("Standardized Mortality Ratio (O/E), 95% CI")
    ax.set_title("Risk-Standardized 30-Day Mortality")
    fig.tight_layout()
    return fig


def _smr_calibration_fig(cal, labels):
    """Transfer calibration: observed vs expected mortality by risk decile per site."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 6))
    mx = 1.0
    for sid, df in cal.items():
        obs = df["observed"] / df["n"] * 100
        exp = df["expected"] / df["n"] * 100
        ax.plot(exp, obs, "o-", color=get_color(sid), alpha=0.85,
                label=labels.get(sid, sid), markersize=5)
        mx = max(mx, float(exp.max()), float(obs.max()))
    lim = min(100.0, mx * 1.08)
    ax.plot([0, lim], [0, lim], ls="--", color="grey", lw=1)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Expected mortality (%)"); ax.set_ylabel("Observed mortality (%)")
    ax.set_title("Transfer Calibration by Risk Decile")
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def _smr_model_card(model) -> str:
    """Reference-model diagnostics card (development coefficients, AUC, VIF)."""
    covs = model.get("covariates", [])
    coef = model.get("coef", {})
    vif = model.get("vif", {})
    rows = [{"site_id": "", "Covariate": html.escape(c, quote=False),
             "Coefficient": f"{coef.get(c, float('nan')):+.4f}",
             "VIF": (f"{vif[c]:.2f}" if c in vif else "—")} for c in covs]
    tbl = _styled_table(pd.DataFrame(rows), display_cols=["Covariate", "Coefficient", "VIF"])
    dev_n = model.get("dev_n", "?")
    meta = (f"Reference model developed on <strong>{html.escape(str(model.get('dev_site','?')), quote=False)}</strong> "
            f"(n={dev_n:,}); pre-specified parsimonious linear logistic. "
            f"Intercept {model.get('intercept', 0):+.3f}. "
            f"Discrimination: c-statistic {model.get('dev_auc','?')} "
            f"(optimism-corrected {model.get('dev_auc_optimism_corrected','?')}). "
            f"Missing predictors: frozen development-median single imputation. "
            f"All VIF near 1 indicates no collinearity.")
    return ('<div class="section"><h2>Reference Model (diagnostics)</h2>'
            f'<p class="fig-caption">{meta}</p>' + tbl + '</div>')


def build_smr_tab(sites, labels) -> str:
    """Risk-Standardized Mortality (SMR) tab: per-site O/E vs an external
    (MIMIC-developed) case-mix reference, with forest, transfer calibration, and
    the reference-model diagnostics. See docs/smr_addition_plan.md."""
    per_site, pooled = collect_smr(sites, labels)
    head = (
        '<div class="section"><h2>Risk-Standardized 30-Day Mortality (SMR)</h2>'
        '<p class="fig-caption">Each site\'s observed 30-day deaths divided by the number '
        '<em>expected</em> from a case-mix model (age, sex, SOFA, baseline lactate, Charlson '
        'index) developed on an independent CLIF dataset (MIMIC-IV) and applied unchanged at '
        'every site. SMR &gt; 1 means more deaths than the case-mix predicts; a 95% CI crossing '
        '1.0 indicates no difference from the reference after adjustment. Because the reference '
        'is external, no participating site is privileged. This is a descriptive benchmark, not '
        'a site quality ranking. The SMR is computed on the same analytic CRRT cohort as Table 1 '
        '(adult, continuous CRRT, non-ESRD, with required baseline data); no minimum '
        'CRRT-duration filter is applied, so early deaths are retained.</p></div>')
    if not per_site:
        return head + _missing(
            "No SMR results yet. Run <code>03b_crrt_epi_smr.py</code> at each site "
            "(after shipping <code>config/smr_reference_model.json</code>).")

    rows = []
    for p in per_site:
        rows.append({"site_id": p["site_id"], "Site": p["label"], "N": f"{p['n']:,}",
                     "Observed": p["observed"], "Expected": f"{p['expected']:.0f}",
                     "SMR (95% CI)": f"{p['smr']:.2f} ({p['lo']:.2f}–{p['hi']:.2f})",
                     "Crude mortality": f"{p['crude_mort_pct']:.1f}%",
                     "c-stat": f"{p['auc']:.3f}"})
    if pooled and pooled["k"] > 1:
        rows.append({"site_id": "", "Site": "<strong>Pooled</strong>", "N": f"{pooled['n']:,}",
                     "Observed": pooled["observed"], "Expected": f"{pooled['expected']:.0f}",
                     "SMR (95% CI)": f"<strong>{pooled['smr']:.2f} ({pooled['lo']:.2f}–{pooled['hi']:.2f})</strong>",
                     "Crude mortality": "", "c-stat": ""})
    cols = ["Site", "N", "Observed", "Expected", "SMR (95% CI)", "Crude mortality", "c-stat"]
    if pooled and pooled.get("het_p") is not None:
        pstr = "&lt;0.001" if pooled["het_p"] < 0.001 else f"{pooled['het_p']:.3f}"
        homo_html = (
            '<p class="fig-caption"><strong>Between-site homogeneity of SMRs:</strong> '
            f'&chi;&sup2; = {pooled["het_Q"]:.2f} (df {pooled["het_df"]}), '
            f'I&sup2; = {pooled["het_I2_pct"]:.0f}%, p = {pstr} '
            '(Poisson / Breslow-Day test, as in CRRTnet). A non-significant p means sites do '
            'not differ in risk-standardized mortality beyond chance, despite practice variation.</p>')
    else:
        k = pooled["k"] if pooled else 0
        homo_html = (
            '<p class="fig-caption"><strong>Between-site homogeneity of SMRs:</strong> '
            f'requires at least 2 sites (currently k = {k}); the &chi;&sup2;, I&sup2;, and p '
            'will populate once additional sites run 03b.</p>')
    table = ('<div class="section"><h2>Per-Site SMR</h2>'
             '<p class="fig-caption">Observed vs expected 30-day deaths; pooled SMR is &Sigma;O/&Sigma;E '
             'with an exact-Poisson (Byar) interval. c-statistic is the model\'s discrimination at '
             'that site (secondary; calibration drives SMR validity).</p>'
             + _styled_table(pd.DataFrame(rows), display_cols=cols)
             + homo_html + '</div>')

    forest = ('<div class="section"><h2>Forest Plot</h2>'
              f'<div class="figure-block"><img src="{fig_to_data_uri(_smr_forest_fig(per_site, pooled))}" '
              'alt="SMR forest plot"></div></div>')

    cal = collect_smr_calibration(sites)
    calib = ""
    if cal:
        calib = ('<div class="section"><h2>Transfer Calibration (diagnostics)</h2>'
                 '<p class="fig-caption">Observed vs expected mortality within predicted-risk deciles '
                 'at each site. Points near the dashed identity line mean the MIMIC-developed model is '
                 'well-calibrated at that site; systematic departure flags case-mix transportability '
                 'limits (the constant-risk caveat) and is the key check for SMR validity.</p>'
                 f'<div class="figure-block"><img src="{fig_to_data_uri(_smr_calibration_fig(cal, labels))}" '
                 'alt="SMR calibration plot"></div></div>')

    model = load_smr_model()
    card = _smr_model_card(model) if model else ""
    return head + table + forest + calib + card


def build_modality_mapping_html(sites, labels) -> str:
    """Cross-site CRRT modality-mapping completeness: % of encounters with a blank
    `crrt_mode_category` ("Unknown") at CRRT initiation, per site, from the existing
    `crrt_practice_quality.csv` Modality rows (patient-level; no re-run needed).
    A high Unknown share is a likely ETL mapping gap (raw crrt_mode_name not mapping
    to a standard category) or a charting-timing artifact (mode not on the very
    first CRRT record). Sorted worst-first, color-graded."""
    rows = []
    for sd in sites:
        pq = _load_csv(sd, "crrt_epi", "crrt_practice_quality")
        if pq is None or pq.empty:
            continue
        mod = pq[pq["variable"] == "Modality"]
        if mod.empty:
            continue
        u = mod[mod["stat"].astype(str).str.contains("UNKNOWN", case=False, na=False)]
        denom = None
        try:
            denom = int(pd.to_numeric(mod["denominator"], errors="coerce").dropna().iloc[0])
        except Exception:
            pass
        if len(u):
            n = pd.to_numeric(u.iloc[0]["n"], errors="coerce")
            pct = pd.to_numeric(u.iloc[0].get("value"), errors="coerce")
            if pd.isna(pct) and pd.notna(n) and denom:
                pct = 100 * n / denom
        else:
            n, pct = 0, 0.0  # no Unknown row = fully mapped
        rows.append({"site_id": sd.name, "Site": labels.get(sd.name, sd.name),
                     "pct": float(pct) if pd.notna(pct) else float("nan"),
                     "n": n, "denom": denom})
    if not rows:
        return ""
    rows.sort(key=lambda r: (-(r["pct"] if r["pct"] == r["pct"] else -1)))

    def _cell(p):
        if p != p:  # NaN
            return '<td style="text-align:right;color:#94a3b8;">&mdash;</td>'
        col = ("#991b1b" if p >= 15 else "#9a3412" if p >= 5 else "#166534")
        return f'<td style="text-align:right;color:{col};font-weight:600;">{p:.1f}%</td>'

    body = []
    for r in rows:
        nd = (f'{"<5" if (isinstance(r["n"], float) and r["n"]!=r["n"]) else int(r["n"])}'
              f'/{r["denom"]}' if r["denom"] else "&mdash;")
        body.append(
            f'<tr style="border-left:4px solid {get_color(r["site_id"])};">'
            f'<td>{html.escape(r["Site"], quote=False)}</td>'
            f'{_cell(r["pct"])}'
            f'<td style="font-size:12px;color:#64748b;">{nd}</td></tr>')
    n_flag = sum(1 for r in rows if r["pct"] == r["pct"] and r["pct"] >= 15)
    return (
        '<div class="section"><h2>CRRT Modality Mapping Completeness</h2>'
        '<p style="font-size:13px;color:#64748b;">% of encounters whose CRRT mode '
        '(<code>crrt_mode_category</code>) is <strong>blank ("Unknown") at initiation</strong> '
        '&mdash; the mode on the CRRT record at the initiation timestamp (patient-level). '
        'A high share points to an <strong>ETL mapping gap</strong> (raw <code>crrt_mode_name</code> '
        'not mapping to a standard category) or a charting-timing artifact (mode absent on the '
        'first record though present moments later). '
        + (f'<strong style="color:#991b1b;">{n_flag} site(s) &ge;15%.</strong>' if n_flag else '')
        + '</p>'
        '<table class="results-table" border="0" style="width:auto;">'
        '<thead><tr><th style="text-align:left;">Site</th><th>Unknown modality</th>'
        '<th>n / cohort</th></tr></thead><tbody>' + "\n".join(body) + '</tbody></table></div>'
    )


def build_initial_settings_table(sites, labels) -> str:
    """Cross-site table of FIRST-3h INITIAL CRRT settings, per site × per mode
    (from step-00 `{SITE}_crrt_initial_settings_by_mode.csv`): BFR, dialysate,
    pre/post replacement, calculated dose, and net UF as median [IQR], aggregated
    one value per patient (first-3h median) then across patients. Purpose:
    cross-site BFR variation, to inform local anticoagulation-practice (e.g.
    citrate) discussions. Sites that have not re-run step 00 under this change are
    silently skipped (transitional)."""
    _COLS = [
        ("blood_flow_rate", "BFR (mL/min)"),
        ("dialysate_flow_rate", "Dialysate (mL/hr)"),
        ("pre_filter_replacement_fluid_rate", "Pre-filter repl. (mL/hr)"),
        ("post_filter_replacement_fluid_rate", "Post-filter repl. (mL/hr)"),
        ("crrt_dose_ml_kg_hr", "Dose (mL/kg/hr)"),
        ("ultrafiltration_out", "Net UF (mL/hr)"),
    ]
    rows = []
    for sd in sites:
        df = _load_csv(sd, "diagnostics", "crrt_initial_settings_by_mode")
        if df is None or df.empty:
            continue
        for _, r in df.iterrows():
            row = {"site_id": sd.name, "Site": labels.get(sd.name, sd.name),
                   "Mode": r.get("Mode", ""), "N": r.get("N_patients", "")}
            for key, _lbl in _COLS:
                row[_lbl] = r.get(key, "—")
            rows.append(row)
    if not rows:
        return (
            '<div class="section"><h2>Initial CRRT Settings (first 3h) by Mode</h2>'
            '<p class="missing">No site has re-run step 00 under the initial-settings '
            'change yet, so this table is empty. It populates as sites re-run '
            '(<code>{SITE}_crrt_initial_settings_by_mode.csv</code>).</p></div>'
        )
    n_sites = len({r["site_id"] for r in rows})
    display = ["Site", "Mode", "N"] + [lbl for _k, lbl in _COLS]
    tbl = _styled_table(pd.DataFrame(rows), display_cols=display)
    return (
        '<div class="section"><h2>Initial CRRT Settings (first 3h) by Mode</h2>'
        f'<p>Median [IQR] of each setting over the first 3&nbsp;h after CRRT '
        'initiation, aggregated <strong>one value per patient then across '
        'patients</strong> (the same window as the dose exposure) &mdash; a true '
        '<em>initial</em> view, not the whole-course record-level distribution. '
        'Rows are per site &times; per mode. <strong>BFR variation</strong> across '
        'sites can reflect local anticoagulation practice (e.g. lower blood-flow '
        'targets with regional citrate); &ldquo;No data&rdquo; means that mode '
        'does not use that fluid (dialysate vs replacement) or the site does not '
        f'chart it. Dose is defined for CVVH/CVVHD/CVVHDF only. Showing {n_sites} '
        f'site(s) that have re-run step 00.</p>' + tbl + '</div>'
    )


def build_crrt_epidemiology(sites, labels, *, anonymized=False) -> str:
    """CRRT descriptive epidemiology tab (from per-site outputs of scripts 07/08):
    incidence/utilization, practice variation + Tier-A quality, low-dose subcohort.
    Sites-as-rows tables; cross-site variation is the headline."""

    # ── Incidence / utilization ──
    def _inc(df, stratum, col):
        r = df[df["stratum"] == stratum]
        return r.iloc[0][col] if len(r) else None

    inc_rows, inc_pool = [], {"den": 0, "num": 0, "vd": 0, "vn": 0, "pd": 0, "pn": 0}
    for sd in sites:
        df = _load_csv(sd, "crrt_epi", "crrt_incidence")
        if df is None or df.empty:
            continue
        den = _inc(df, "Overall (adult ICU)", "n_denominator")
        num = _inc(df, "Overall (adult ICU)", "n_crrt")
        vd = _inc(df, "On invasive ventilation", "n_denominator")
        vn = _inc(df, "On invasive ventilation", "n_crrt")
        pd_, pn = _inc(df, "On vasopressors", "n_denominator"), _inc(df, "On vasopressors", "n_crrt")
        inc_rows.append({
            "site_id": sd.name, "Site": labels.get(sd.name, sd.name),
            "Adult ICU (N)": f"{int(den):,}" if den else "—",
            "CRRT (n)": f"{int(num):,}" if num else "—",
            "Overall %": _inc(df, "Overall (adult ICU)", "incidence_pct"),
            "Ventilated %": _inc(df, "On invasive ventilation", "incidence_pct"),
            "Vasopressors %": _inc(df, "On vasopressors", "incidence_pct"),
        })
        for k, v in [("den", den), ("num", num), ("vd", vd), ("vn", vn), ("pd", pd_), ("pn", pn)]:
            inc_pool[k] += int(v) if v else 0
    if inc_rows and len(inc_rows) > 1:
        p = inc_pool
        inc_rows.append({
            "site_id": "", "Site": "Pooled", "Adult ICU (N)": f"{p['den']:,}", "CRRT (n)": f"{p['num']:,}",
            "Overall %": round(100 * p["num"] / p["den"], 2) if p["den"] else None,
            "Ventilated %": round(100 * p["vn"] / p["vd"], 2) if p["vd"] else None,
            "Vasopressors %": round(100 * p["pn"] / p["pd"], 2) if p["pd"] else None,
        })
    inc_cols = ["Site", "Adult ICU (N)", "CRRT (n)", "Overall %", "Ventilated %", "Vasopressors %"]
    inc_html = (_styled_table(pd.DataFrame(inc_rows), display_cols=inc_cols)
                if inc_rows else "<p>No incidence data available yet.</p>")

    # ── Practice variation + Tier-A quality ──
    def _pq(df, variable, stat):
        r = df[(df["variable"] == variable) & (df["stat"] == stat)]
        return r.iloc[0]["value"] if len(r) else "—"

    def _uf_intensity(sd):
        # Course-average net-UF intensity (mL/kg/hr): median [IQR] from the
        # per-site summary written by 03 build_uf_trajectories().
        s = _load_csv(sd, "crrt_epi/graphs", "net_uf_intensity_summary")
        if s is None or s.empty:
            return "—"
        r = s.iloc[0]
        return f"{r['median']:.2f} [{r['p25']:.2f}, {r['p75']:.2f}]"

    pq_rows = []
    for sd in sites:
        df = _load_csv(sd, "crrt_epi", "crrt_practice_quality")
        if df is None or df.empty:
            continue
        pq_rows.append({
            "site_id": sd.name, "Site": labels.get(sd.name, sd.name),
            "Dose (mL/kg/hr)": _pq(df, "CRRT dose (mL/kg/hr)", "median_iqr"),
            "% <20": _pq(df, "Dose band", "pct_<20"),
            "% 20-30": _pq(df, "Dose band", "pct_20-30 (KDIGO)"),
            "% >30": _pq(df, "Dose band", "pct_>30"),
            "Net UF intensity (mL/kg/hr, 72h)": _uf_intensity(sd),
            "Time to CRRT (h)": _pq(df, "Time to CRRT from first vital (h)", "median_iqr"),
            "Duration (d)": _pq(df, "CRRT duration (days)", "median_iqr"),
            "30d mortality %": _pq(df, "30-day mortality", "pct"),
        })
    pq_cols = ["Site", "Dose (mL/kg/hr)", "% <20", "% 20-30", "% >30",
               "Net UF intensity (mL/kg/hr, 72h)", "Time to CRRT (h)", "Duration (d)",
               "30d mortality %"]
    pq_html = (_styled_table(pd.DataFrame(pq_rows), display_cols=pq_cols)
               if pq_rows else "<p>No practice/quality data available yet.</p>")

    # ── Between-site heterogeneity ──
    # One row per descriptive epi metric: random-effects pooled estimate (95% CI),
    # the number of contributing sites, I², and the direct between-site test p.
    # Single source of truth = compute_epi_heterogeneity() (also written to
    # pooled_epi_heterogeneity.csv by 08). Q/I²/p are undefined at k=1.
    het = compute_epi_heterogeneity(sites)
    het_rows = []
    for g in dict.fromkeys(r["group"] for r in het):  # preserve first-seen order
        het_rows.append({"site_id": "", "Metric": f"<strong>{html.escape(g)}</strong>",
                         "Sites (k)": "", "Pooled (95% CI)": "", "Range (sites)": "", "I²": "",
                         "Between-site p": ""})
        for r in (x for x in het if x["group"] == g):
            is_prop = r["scale"] == "proportion (%)"
            lo, hi = r["pooled_lo"], r["pooled_hi"]
            suffix = "%" if is_prop else ""
            ci = f"{r['pooled']:g} ({lo:g}, {hi:g}){suffix}"
            smin, smax = r.get("site_min"), r.get("site_max")
            if pd.notna(smin) and pd.notna(smax):
                rng = f"{smin:g}{suffix}" if smin == smax else f"{smin:g}-{smax:g}{suffix}"
            else:
                rng = "—"
            k = int(r["k_sites"])
            i2 = f"{r['I2_pct']:g}%" if k > 1 and pd.notna(r["I2_pct"]) else "—"
            bp = r["between_site_p"]
            bp_disp = ("&lt;0.001" if (k > 1 and pd.notna(bp) and bp < 0.001)
                       else f"{bp:g}" if (k > 1 and pd.notna(bp)) else "—")
            metric = str(r["metric"])
            if metric.startswith("ICU subtype: "):  # burn_icu -> Burn ICU
                metric = (metric.replace("ICU subtype: ", "").replace("_", " ")
                          .strip().title().replace("Icu", "ICU"))
            het_rows.append({
                "site_id": "", "Metric": " " + html.escape(metric),
                "Sites (k)": k, "Pooled (95% CI)": ci, "Range (sites)": rng, "I²": i2,
                "Between-site p": bp_disp})
    het_cols = ["Metric", "Sites (k)", "Pooled (95% CI)", "Range (sites)", "I²",
                "Between-site p"]
    het_html = (_styled_table(pd.DataFrame(het_rows), display_cols=het_cols)
                if het_rows else "<p>No heterogeneity data available yet.</p>")

    # ── CRRT dose: actual vs ideal body weight (height-available paired subset) ──
    # Single source of truth = compute_dose_ibw_pooled() (also written to
    # pooled_dose_ibw_comparison.csv by 08). Shows how normalizing delivered dose
    # to Devine IBW instead of actual weight reclassifies encounters across the
    # 30 mL/kg/hr cutoff. Percentages/medians are pooled; k=1 == single site.
    dib = compute_dose_ibw_pooled(sites)
    if dib is None:
        dib_html = "<p>No dose-by-IBW comparison data available yet.</p>"
    else:
        _p = lambda x: f"{x:.0f}%" if pd.notna(x) else "—"
        _m = lambda x: f"{x:.1f}" if pd.notna(x) else "—"
        dib_rows = [
            {"site_id": "", "Metric": "Median dose (mL/kg/hr)",
             "Actual body weight": _m(dib["median_dose_actual"]),
             "Ideal body weight (Devine)": _m(dib["median_dose_ibw"])},
            {"site_id": "", "Metric": "&ge;30 mL/kg/hr (high-intensity)",
             "Actual body weight": _p(dib["pct_ge30_actual"]),
             "Ideal body weight (Devine)": _p(dib["pct_ge30_ibw"])},
            {"site_id": "", "Metric": "Within KDIGO 20&ndash;30 band",
             "Actual body weight": _p(dib["pct_kdigo_actual"]),
             "Ideal body weight (Devine)": _p(dib["pct_kdigo_ibw"])},
            {"site_id": "", "Metric": "&lt;15 mL/kg/hr (very low)",
             "Actual body weight": _p(dib["pct_lt15_actual"]),
             "Ideal body weight (Devine)": _p(dib["pct_lt15_ibw"])},
        ]
        dib_cols = ["Metric", "Actual body weight", "Ideal body weight (Devine)"]
        k = dib["k_sites"]
        i2 = f"{dib['ge30_ibw_I2']:g}%" if k > 1 and pd.notna(dib["ge30_ibw_I2"]) else "—"
        net30 = dib["counts"]["n_ibw_ge30"] - dib["counts"]["n_actual_ge30"]
        dib_html = (
            _styled_table(pd.DataFrame(dib_rows), display_cols=dib_cols)
            + f'<p style="font-size:12px;color:#64748b;">Paired (height-available) '
            f'encounters n = {dib["n_paired"]:,} across {k} site(s). Median '
            f'actual/IBW weight ratio {dib["median_weight_ibw_ratio"]:.2f}; '
            f'{_p(dib["pct_heavier_than_ideal"])} heavier than ideal weight. '
            f'Normalizing delivered dose to ideal body weight moves the '
            f'high-intensity (&ge;30) share from {_p(dib["pct_ge30_actual"])} to '
            f'{_p(dib["pct_ge30_ibw"])} (net {net30:+,} encounters across the 30 '
            f'mL/kg/hr cutoff); {dib["n_left_kdigo_ibw"]:,} left and '
            f'{dib["n_entered_kdigo_ibw"]:,} entered the KDIGO 20&ndash;30 band. '
            f'Between-site I&sup2; for the &ge;30 IBW share: {i2}.</p>')

    # ── Low-dose (10-15 mL/kg/hr) subcohort ──
    def _band(df, band, col):
        r = df[df["band"] == band]
        return r.iloc[0][col] if len(r) else None

    ld_rows, ld_pool = [], {"low": 0, "tot": 0}
    for sd in sites:
        df = _load_csv(sd, "low_dose", "low_dose_counts")
        if df is None or df.empty:
            continue
        low = _band(df, "10-15 (very low)", "n")
        tot = _band(df, "Total with dose", "n")
        ld_rows.append({
            "site_id": sd.name, "Site": labels.get(sd.name, sd.name),
            "Dosed (N)": f"{int(tot):,}" if tot else "—",
            "Very low 10-15 (n)": int(low) if low is not None else "—",
            "Very low %": _band(df, "10-15 (very low)", "pct_of_dosed"),
            "<10 (n)": int(_band(df, "<10", "n") or 0),
        })
        ld_pool["low"] += int(low) if low else 0
        ld_pool["tot"] += int(tot) if tot else 0
    if ld_rows and len(ld_rows) > 1:
        ld_rows.append({
            "site_id": "", "Site": "Pooled", "Dosed (N)": f"{ld_pool['tot']:,}",
            "Very low 10-15 (n)": ld_pool["low"],
            "Very low %": round(100 * ld_pool["low"] / ld_pool["tot"], 1) if ld_pool["tot"] else None,
            "<10 (n)": "",
        })
    ld_cols = ["Site", "Dosed (N)", "Very low 10-15 (n)", "Very low %", "<10 (n)"]
    ld_html = (_styled_table(pd.DataFrame(ld_rows), display_cols=ld_cols)
               if ld_rows else "<p>No low-dose data available yet.</p>")

    # ── Per-site figure gallery (07 epidemiology + 08 low-dose figures) ──
    gallery_specs = [
        ("crrt_epi/graphs", "crrt_incidence_by_context", "CRRT incidence by clinical context"),
        ("crrt_epi/graphs", "crrt_incidence_by_icu_subtype", "CRRT incidence by ICU subtype"),
        ("crrt_epi/graphs", "dose_distribution", "CRRT dose distribution (KDIGO 20-30 band)"),
        ("crrt_epi/graphs", "dose_by_ibw", "CRRT dose: actual vs ideal body weight (Devine)"),
        ("crrt_epi/graphs", "uf_mortality", "Crude 30-day mortality vs first-72h net ultrafiltration intensity"),
        ("crrt_epi/graphs", "crrt_dose_over_time", "CRRT dose over 30 days"),
        ("crrt_epi/graphs", "net_uf_rate_over_time", "Net ultrafiltration rate over 30 days"),
        ("crrt_epi/graphs", "net_uf_cumulative_over_time", "Cumulative net ultrafiltration over 30 days"),
        ("crrt_epi/graphs", "lab_distributions_over_crrt", "Lab trajectories over 30 days"),
        ("crrt_epi/graphs", "nee_over_crrt", "Vasopressor (NEE) over 30 days"),
        ("crrt_epi/graphs", "imv_state_over_crrt", "Invasive mechanical ventilation state over 30 days (IMV / discharge / death)"),
        ("crrt_epi/graphs", "crrt_state_over_crrt", "CRRT state over 30 days (on/off CRRT, discharged, dead)"),
    ]
    gallery_blocks = []
    for sd in sites:
        cards = []
        for subdir, stem, caption in gallery_specs:
            p = _find_png(sd, subdir, stem)
            if p is None:
                continue
            uri = img_to_data_uri(p)
            cards.append(
                '<figure style="margin:0;flex:1 1 380px;max-width:540px;">'
                f'<img src="{uri}" alt="{html.escape(caption)}" '
                'style="width:100%;border:1px solid #e2e8f0;border-radius:6px;"/>'
                '<figcaption style="font-size:12px;color:#64748b;margin-top:4px;">'
                f'{html.escape(caption)}</figcaption></figure>'
            )
        if cards:
            site_lbl = str(labels.get(sd.name, sd.name))
            color = get_color(sd.name)
            gallery_blocks.append(
                f'<h3 style="margin:18px 0 8px;color:#334155;font-size:15px;'
                f'border-left:4px solid {color};padding-left:8px;">{html.escape(site_lbl)}</h3>'
                f'<div style="display:flex;flex-wrap:wrap;gap:16px;">{"".join(cards)}</div>'
            )
    gallery_html = ("".join(gallery_blocks) if gallery_blocks
                    else "<p>No epidemiology figures available yet.</p>")

    # ── Cross-site pooled figures (08 pooled_descriptive_*, anonymized Site 1..N) ──
    pooled_specs = _POOLED_FIGURE_SPECS
    pooled_cards = []
    for stem, caption in pooled_specs:
        p = ROOT / "figures" / f"{stem}.png"
        # In the non-anonymized dashboard, prefer the real-site-name "_named"
        # variant (emitted by 08 for the site-identifying figures); the pooled-only
        # figures have no _named variant and fall through to the canonical stem.
        if not anonymized:
            named_p = ROOT / "figures" / f"{stem}_named.png"
            if named_p.exists():
                p = named_p
        if not p.exists():
            continue
        uri = img_to_data_uri(p)
        pooled_cards.append(
            '<figure style="margin:0;flex:1 1 380px;max-width:540px;">'
            f'<img src="{uri}" alt="{html.escape(caption)}" '
            'style="width:100%;border:1px solid #e2e8f0;border-radius:6px;"/>'
            '<figcaption style="font-size:12px;color:#64748b;margin-top:4px;">'
            f'{html.escape(caption)}</figcaption></figure>'
        )
    pooled_fig_html = ('<div style="display:flex;flex-wrap:wrap;gap:16px;">'
                       + "".join(pooled_cards) + '</div>' if pooled_cards
                       else "<p>No cross-site pooled figures available yet.</p>")

    return (
        '<div class="section"><h2>CRRT Incidence &amp; Utilization</h2>'
        '<p>Each percentage is the fraction of that column&rsquo;s denominator that '
        'received CRRT at any point during the hospitalization (numerator = CRRT '
        'recipients, denominator = the stratum). For example, Ventilated % = CRRT '
        'recipients &divide; all ventilated adult ICU hospitalizations, <em>not</em> '
        'the share of CRRT patients who were ventilated. Denominators are '
        'hospitalization-level (not unique patients) and overlap across strata. '
        'Cross-site variation reflects differing case-mix and practice.</p>'
        + inc_html + '</div>'
        '<div class="section"><h2>Practice Variation &amp; Quality (Tier A)</h2>'
        '<p>CRRT dose (prescribed, as recorded in CLIF) vs the displayed 20-30 '
        'mL/kg/hr band, which spans prescribed targets that achieve the '
        'KDIGO-recommended 20-25 delivered dose; net ultrafiltration rate vs the '
        'Murugan 1.01-1.75 mL/kg/hr band; timing, duration, and 30-day mortality '
        'in the analytic CRRT cohort.</p>'
        + pq_html + '</div>'
        + build_modality_mapping_html(sites, labels)
        + build_initial_settings_table(sites, labels)
        + '<div class="section"><h2>Between-Site Heterogeneity</h2>'
        '<p>Random-effects pooled estimate (95% CI) for each descriptive metric '
        'with the number of contributing sites, I&sup2; (share of variance from '
        'between-site heterogeneity), and the direct between-site test p '
        '(chi-square homogeneity for proportions, Cochran&rsquo;s Q for means). '
        'I&sup2; and the between-site p require at least two sites and read '
        '&ldquo;&mdash;&rdquo; until additional sites contribute.</p>'
        + het_html + '</div>'
        '<div class="section"><h2>Very-Low-Dose Subcohort (10–15 mL/kg/hr)</h2>'
        '<p>Size of the very-low-dose subcohort targeted by upcoming RCTs '
        '(e.g. NCT06021288), informing trial feasibility across the consortium.</p>'
        + ld_html + '</div>'
        '<div class="section"><h2>CRRT Dose: Actual vs Ideal Body Weight</h2>'
        '<p>Delivered CRRT dose is conventionally normalized to actual body '
        'weight. Because most CRRT patients exceed their ideal body weight, '
        'normalizing to Devine ideal body weight (IBW) shifts the dose '
        'distribution rightward and reclassifies a substantial share of '
        'encounters across the 30 mL/kg/hr high/low cutoff used in the causal '
        'analysis. Computed on the height-available paired subset, so both '
        'columns share the same denominator.</p>'
        + dib_html + '</div>'
        '<div class="section"><h2>Cross-Site Pooled Figures</h2>'
        '<p>Consortium-level descriptive figures pooled across sites (script 08), '
        + ('with sites anonymized as Site 1&hellip;N. ' if anonymized
           else 'with sites shown by name. ')
        + 'The dose distribution sums fixed-width '
        'histogram bins across sites; per-site outlines show practice variation '
        'against the pooled shape.</p>'
        + pooled_fig_html + '</div>'
        '<div class="section"><h2>Per-Site Figures</h2>'
        '<p>Descriptive epidemiology and low-dose figures rendered per site '
        '(scripts 07 / 08). Cross-site variation is the headline.</p>'
        + gallery_html + '</div>'
    )


def build_dashboard(sites, labels, output_path, *, anonymized=False):
    tab_panels = []

    # ── INITIAL POINT-TREATMENT GROUP ──
    print("  Building Initial Point-Treatment tab...")
    overview = (
        build_sample_characteristics(sites, labels)
        + build_forest_plot_iptw(sites, labels)
        + build_forest_plot_psm_fg(sites, labels)
        + build_model_comparison(sites, labels)
        + build_combined_cif_psm(sites, labels)
        + build_combined_cif_iptw(sites, labels)
        + build_dose_response(sites, labels)
        + build_evalue(sites, labels)
    )
    tab_panels.append(f'<div class="tab-panel" id="panel-overview" style="display:none;">{overview}</div>')

    print("  Building Subgroups (Initial) tab...")
    subgroup = build_subgroup_analysis(sites, labels)
    tab_panels.append(f'<div class="tab-panel" id="panel-subgroup" style="display:none;">{subgroup}</div>')

    print("  Building IPTW Diagnostics tab...")
    diag_iptw = build_diagnostics_iptw(sites, labels)
    tab_panels.append(f'<div class="tab-panel" id="panel-diag_iptw" style="display:none;">{diag_iptw}</div>')

    pt_group = _tab_group("Point Treatment", [
        _tab_btn("overview", "Results"),
        _tab_btn("subgroup", "Subgroups"),
        _tab_btn("diag_iptw", "Diagnostics"),
    ])

    # ── POOLED / TABLES GROUP ──
    print("  Building Meta-Analysis tab...")
    meta = build_meta_analysis(sites, labels)
    tab_panels.append(f'<div class="tab-panel" id="panel-meta" style="display:none;">{meta}</div>')

    print("  Building Combined Tables tab...")
    combined_tables = build_combined_tables(sites, labels)
    tab_panels.append(f'<div class="tab-panel" id="panel-tables" style="display:none;">{combined_tables}</div>')

    # The Completeness / Eligibility / not-contributing QC tabs are internal
    # coordinating-center review tools whose file listings inherently show real
    # site-named paths — omit them entirely from the anonymized build.
    pooled_btns = [_tab_btn("meta", "Meta-Analysis"), _tab_btn("tables", "Combined Tables")]
    if not anonymized:
        pooled_btns.append(_tab_btn("completeness", "Completeness"))
    pooled_group = _tab_group("Pooled", pooled_btns)

    # ── DESCRIPTIVE EPIDEMIOLOGY GROUP ──
    print("  Building CRRT Epidemiology tab...")
    epi = build_crrt_epidemiology(sites, labels, anonymized=anonymized)
    tab_panels.append(f'<div class="tab-panel" id="panel-epi" style="display:block;">{epi}</div>')

    epi_group = _tab_group("Descriptive", [
        _tab_btn("epi", "CRRT Epidemiology", active=True),
    ])

    print("  Building SMR tab...")
    smr = build_smr_tab(sites, labels)
    tab_panels.append(f'<div class="tab-panel" id="panel-smr" style="display:none;">{smr}</div>')
    smr_group = _tab_group("Risk-Standardized", [_tab_btn("smr", "SMR")])

    # ── Per-site tabs ──
    # In the anonymized build the tab id (used in DOM ids, switchTab(), and the
    # JS id array) must NOT be the real site name — derive a slug from the anon
    # label ("Site 8" → "site8") so the HTML is source-clean, not just visually.
    def _tab_id(sd):
        return _slugify(labels.get(sd.name, sd.name)) if anonymized else sd.name
    site_btns = []
    site_tab_ids = []
    for sd in sites:
        sid = sd.name
        lbl = labels.get(sid, sid)
        tid = _tab_id(sd)
        site_tab_ids.append(tid)
        print(f"  Building {lbl} tab...")
        content = build_site_tab(sd)
        site_btns.append(_tab_btn(tid, lbl))
        tab_panels.append(f'<div class="tab-panel" id="panel-{tid}" style="display:none;">{content}</div>')

    site_group = _tab_group("Per-Site", site_btns)

    # Completeness is built LAST, after every other tab, so the build-time
    # `CONSUMED` set (populated by report_core's resolvers as each tab renders)
    # is complete before the checklist diffs it against the files on disk to
    # flag present-but-not-contributing outputs. Omitted from the anonymized
    # build (internal QC tool with site-named file paths).
    if not anonymized:
        print("  Building Completeness Checklist tab...")
        completeness = build_completeness_checklist(sites, labels)
        tab_panels.append(f'<div class="tab-panel" id="panel-completeness" style="display:none;">{completeness}</div>')

    # Descriptive epidemiology leads (epi tab is the default); causal/pooled follow.
    tabs_bar = "\n".join([epi_group, smr_group, pt_group, pooled_group, site_group])
    panels = "\n".join(tab_panels)
    all_ids = ["epi", "smr", "overview", "meta", "tables", "subgroup", "diag_iptw"]
    if not anonymized:
        all_ids.append("completeness")
    all_ids += site_tab_ids
    site_ids_js = ", ".join(f'"{t}"' for t in all_ids)

    from datetime import datetime
    gen_time = datetime.now().strftime("%B %d, %Y at %I:%M %p")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>CRRT Epidemiology and Dose Causal Inference Dashboard</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: Inter, -apple-system, 'Segoe UI', system-ui, sans-serif;
            margin: 0; padding: 24px;
            background: #f8f9fa;
            color: #1e293b;
            line-height: 1.55;
            font-size: 14px;
        }}
        .container {{
            max-width: 1600px; margin: 0 auto;
            background: white;
            padding: 40px 48px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.06), 0 4px 16px rgba(0,0,0,0.04);
            border-radius: 8px;
        }}

        /* ── Header ── */
        .dashboard-header {{
            border-bottom: 2px solid #e2e8f0;
            padding-bottom: 16px;
            margin-bottom: 28px;
        }}
        .dashboard-header h1 {{
            color: #0f172a;
            font-size: 26px;
            font-weight: 700;
            margin: 0 0 4px 0;
            letter-spacing: -0.3px;
        }}
        .dashboard-header .subtitle {{
            color: #64748b;
            font-size: 15px;
            font-weight: 400;
            margin: 0;
        }}
        .dashboard-header .gen-time {{
            color: #94a3b8;
            font-size: 12px;
            margin-top: 6px;
        }}

        /* ── Tab Bar ── */
        .tab-nav {{
            display: flex;
            gap: 24px;
            margin-bottom: 28px;
            flex-wrap: wrap;
        }}
        .tab-group {{
            display: flex;
            flex-direction: column;
            gap: 6px;
        }}
        .tab-group-label {{
            font-size: 10px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.8px;
            color: #94a3b8;
            padding-left: 4px;
        }}
        .tab-group-buttons {{
            display: flex;
            gap: 3px;
        }}
        .tab-btn {{
            padding: 8px 16px;
            border: 1px solid #e2e8f0;
            background: #f8fafc;
            cursor: pointer;
            font-size: 13px;
            font-weight: 500;
            border-radius: 6px;
            color: #475569;
            transition: all 0.15s ease;
        }}
        .tab-btn:hover {{
            background: #f1f5f9;
            border-color: #cbd5e1;
            color: #1e293b;
        }}
        .tab-btn.active {{
            background: #0f766e;
            color: white;
            border-color: #0f766e;
            box-shadow: 0 1px 3px rgba(15,118,110,0.25);
        }}

        /* ── Sections & Cards ── */
        .section {{
            margin-bottom: 48px;
            padding: 28px;
            background: #ffffff;
            border: 1px solid #f1f5f9;
            border-radius: 8px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.03);
        }}
        .section h2 {{
            color: #0f172a;
            font-size: 19px;
            font-weight: 600;
            border-bottom: 1px solid #e2e8f0;
            padding-bottom: 8px;
            margin-top: 0;
            margin-bottom: 20px;
        }}

        /* ── Tables ── */
        table {{
            width: auto;
            border-collapse: collapse;
            margin-top: 16px;
            font-size: 13px;
        }}
        thead {{ border-bottom: 2px solid #334155; }}
        th {{
            background: #f8fafc;
            color: #334155;
            font-weight: 600;
            padding: 10px 12px;
            text-align: left;
            border: none;
            border-bottom: 2px solid #334155;
        }}
        td {{
            padding: 9px 12px;
            border: none;
            border-bottom: 1px solid #f1f5f9;
            vertical-align: top;
            text-align: left;
        }}
        tr:nth-child(even) {{ background: #fafbfc; }}
        tr:hover {{ background: #f1f5f9; }}
        tbody tr:last-child td {{ border-bottom: 1px solid #cbd5e1; }}

        /* ── Figures ── */
        .figure-block {{
            margin-bottom: 36px;
            text-align: center;
        }}
        .figure-block h3 {{
            color: #334155;
            font-size: 15px;
            font-weight: 600;
            margin-bottom: 12px;
            text-align: left;
        }}
        .figure-block img {{
            max-width: 100%;
            height: auto;
            border-radius: 6px;
            box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        }}
        .fig-caption {{
            color: #64748b;
            font-size: 13px;
            margin: 4px 0 10px 0;
            text-align: left;
        }}

        /* ── Missing Data ── */
        .missing {{
            background: #fffbeb;
            border: 1px solid #fde68a;
            border-radius: 6px;
            padding: 12px 16px;
            color: #92400e;
            font-size: 13px;
            font-style: normal;
        }}
        .missing::before {{
            content: "— ";
            font-weight: 600;
        }}

        /* ── Miscellaneous ── */
        em {{ color: #64748b; }}
        p {{ margin: 8px 0; }}
        h3 {{ color: #334155; }}
        h4 {{ color: #475569; font-weight: 600; font-size: 14px; }}
    </style>
</head>
<body>
<div class="container">
    <div class="dashboard-header">
        <h1>CRRT Epidemiology and Dose Causal Inference Dashboard</h1>
        <p class="subtitle">CLIF Consortium — Multi-Site Results</p>
        <p class="gen-time">Generated {gen_time}</p>
    </div>
    <div class="tab-nav">{tabs_bar}</div>
    {panels}
</div>
<script>
    const TABS = [{site_ids_js}];
    function switchTab(id) {{
        TABS.forEach(function(t) {{
            document.getElementById('panel-' + t).style.display = 'none';
            document.getElementById('btn-' + t).classList.remove('active');
        }});
        document.getElementById('panel-' + id).style.display = 'block';
        document.getElementById('btn-' + id).classList.add('active');
    }}
</script>
</body>
</html>"""

    output_path.write_text(html, encoding="utf-8")
    print(f"  Written to {output_path}")
    print(f"  File size: {output_path.stat().st_size / 1024 / 1024:.1f} MB")


OUTPUT_ANON = ROOT / "crrt_dashboard_anonymized.html"


def _warn_if_pooled_figures_stale(sites) -> bool:
    """Warn when the pooled figures this dashboard embeds predate the site data.

    The dashboard's "Cross-Site Pooled Figures" section embeds finished PNGs
    written by `08_manuscript_artifacts.py`; 07 never regenerates them. So a
    site drop followed by 07-alone silently ships last build's pooled figures —
    a wrong dashboard that looks fine. This compares every embedded stem against
    the newest collected per-site file and says so out loud.

    Warning, not a hard failure: rebuilding 07 alone is legitimate when only the
    dashboard's own layout/tables changed. Returns True when something is stale.
    """
    newest, newest_src = 0.0, None
    for sd in sites:
        for p in Path(_site_final_dir(sd)).rglob("*"):
            if p.is_file() and p.suffix.lower() in (".csv", ".png", ".pdf", ".html"):
                mt = p.stat().st_mtime
                if mt > newest:
                    newest, newest_src = mt, p

    figs = ROOT / "figures"
    stale, missing = [], []
    for stem, _ in _POOLED_FIGURE_SPECS:
        candidates = [figs / f"{stem}.png", figs / f"{stem}_named.png"]
        present = [p for p in candidates if p.exists()]
        if not present:
            missing.append(stem)
        else:
            stale += [p.name for p in present if p.stat().st_mtime < newest]

    if not (stale or missing):
        return False
    print("\n" + "!"*78)
    print("  STALE POOLED FIGURES — run 08_manuscript_artifacts.py, then re-run 07.")
    print(f"  Newest collected site file: {newest_src.relative_to(ROOT) if newest_src else '?'}"
          f"  ({datetime.fromtimestamp(newest):%Y-%m-%d %H:%M})")
    if stale:
        print(f"  {len(stale)} pooled figure(s) older than that: {', '.join(sorted(stale)[:6])}"
              + (f", +{len(stale)-6} more" if len(stale) > 6 else ""))
    if missing:
        print(f"  {len(missing)} never built: {', '.join(missing)}")
    print("  The dashboard will still build, but its pooled figures describe an older site set.")
    print("!"*78)
    return True


def main():
    sites = discover_sites()
    if not sites:
        print("No sites found in", ROOT)
        return
    print(f"Found {len(sites)} sites: {[s.name for s in sites]}")

    # 08 owns the pooled figures this dashboard embeds; shout if they predate
    # the collected site data rather than shipping a stale-figure dashboard.
    _warn_if_pooled_figures_stale(sites)

    # Refresh the gate verdicts FIRST: build_pooling_eligibility_html reads
    # site_validation.json, so running the gates after the build would embed a
    # stale ledger. Also rewrites the AUTO blocks of the review tracker.
    print("\nRefreshing site validation (gates -> site_validation.json + tracker)...")
    site_validation.run(root=ROOT, verbose=False)

    # Manuscript artifacts (pooled CSVs, figures, tables) live in
    # `code/08_manuscript_artifacts.py` — run that script separately to
    # rebuild them. They are independent of the dashboard.

    print("\nBuilding causal inference dashboard...")
    build_dashboard(sites, SITE_LABELS, OUTPUT)

    # Anonymized version: Site 1, Site 2, ...
    anon_labels = {sd.name: f"Site {i+1}" for i, sd in enumerate(sites)}
    anon_colors = {sd.name: _EXTRA_COLORS[i % len(_EXTRA_COLORS)] for i, sd in enumerate(sites)}
    # Temporarily override SITE_COLORS for anonymized run
    saved_colors = SITE_COLORS.copy()
    SITE_COLORS.clear()
    SITE_COLORS.update(anon_colors)
    if hasattr(get_color, "_cache"):
        get_color._cache.clear()

    print("\nBuilding anonymized dashboard...")
    build_dashboard(sites, anon_labels, OUTPUT_ANON, anonymized=True)

    # Restore original colors
    SITE_COLORS.clear()
    SITE_COLORS.update(saved_colors)
    if hasattr(get_color, "_cache"):
        get_color._cache.clear()


if __name__ == "__main__":
    main()
