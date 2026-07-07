"""
06_low_dose_characterization.py — very-low-dose CRRT subcohort (descriptive).

Characterizes the very-low-dose CRRT subcohort (dose 10-15 mL/kg/hr,
the range targeted by upcoming low-dose RCTs, e.g. NCT06021288) across the
consortium, to inform trial feasibility and external validity. This is the
DESCRIPTIVE companion to the causal emulation in 05c_low_dose_emulation.R — it
makes no causal claim; it answers "how many such patients exist at each site and
how do they differ from the rest of the CRRT cohort?"

Emits AGGREGATE per-site CSVs only (no patient-level data). Pooling + rendering
happen in 10 (manuscript). Three products, all keyed by site:

  {SITE}_low_dose_counts.csv   dose-band tally (the RCT-feasibility numbers)
  {SITE}_low_dose_long.csv     long-format Very-low(10-15) vs Rest comparison
                               (gtsummary-shaped)
  {SITE}_low_dose_table.csv    standalone readable two-column table + p-values

Requires has_crrt_settings=True (dose-defined subcohort). Skips gracefully otherwise.

Run standalone (from repo root):
    uv run python code/06_low_dose_characterization.py
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from pipeline_helpers import load_config

warnings.filterwarnings("ignore")

config = load_config()  # honors CLIF_CONFIG; defaults to config/config.json

SITE_NAME = config["site_name"]
HAS_CRRT_SETTINGS = config.get("has_crrt_settings", True)

INTER = Path("../output/intermediate_phi")
OUT = Path("../output/final_no_phi/low_dose")
OUT.mkdir(parents=True, exist_ok=True)

LOW_LO, LOW_HI = 10.0, 15.0  # very-low-dose band (mL/kg/hr), inclusive

t_start = time.time()
print(f"=== 06 low-dose CRRT characterization — {SITE_NAME} ===")

if not HAS_CRRT_SETTINGS:
    print("  has_crrt_settings=False — dose-defined subcohort not computable; skipping.")
    raise SystemExit(0)

cohort = pd.read_parquet(INTER / "tableone_analysis_df.parquet")
dose = pd.to_numeric(cohort["crrt_dose_ml_kg_hr"], errors="coerce")

# ── Dose-band tally (RCT feasibility) ───────────────────────────────────────
bands = [
    ("<10", dose < 10),
    ("10-15 (very low)", (dose >= 10) & (dose <= 15)),
    (">15-20", (dose > 15) & (dose <= 20)),
    (">20-25", (dose > 20) & (dose <= 25)),
    (">25-30", (dose > 25) & (dose <= 30)),
    (">30", dose > 30),
]
n_dose = int(dose.notna().sum())
count_rows = [{"site": SITE_NAME, "band": label, "n": int(m.sum()),
               "pct_of_dosed": round(100 * m.sum() / n_dose, 1) if n_dose else np.nan}
              for label, m in bands]
count_rows.append({"site": SITE_NAME, "band": "Total with dose", "n": n_dose, "pct_of_dosed": 100.0})
counts = pd.DataFrame(count_rows)
counts.to_csv(OUT / f"{SITE_NAME}_low_dose_counts.csv", index=False)
n_low = int(((dose >= LOW_LO) & (dose <= LOW_HI)).sum())
print(f"  very-low-dose (10-15 mL/kg/hr) subcohort: n={n_low}  of {n_dose} dosed encounters")

# NOTE: the dose-band distribution bar chart was retired — it duplicated the
# continuous initial-dose histogram in 03 (which shows shape + the very-low /
# KDIGO reference bands). Per-band proportions live in the practice-quality CSV
# (03) and the descriptive table (10). The RCT-feasibility counts CSV above is
# retained for the dashboard's very-low (10-15) lookup.

# ── Very-low vs Rest comparison ─────────────────────────────────────────────
grp = pd.Series(np.where((dose >= LOW_LO) & (dose <= LOW_HI), "Very low (10-15)",
                         np.where(dose.notna(), "Rest", None)), index=cohort.index)
low = cohort[grp == "Very low (10-15)"]
rest = cohort[grp == "Rest"]
GROUPS = [("Very low (10-15)", low), ("Rest", rest)]
print(f"  comparison groups: Very low n={len(low)}  |  Rest n={len(rest)}")

long_rows: list = []
table_rows: list = []  # readable: [variable, very_low_str, rest_str, pval]


def _fmt_mi(s: pd.Series) -> str:
    v = pd.to_numeric(s, errors="coerce").dropna()
    if v.empty:
        return "—"
    return f"{v.median():.1f} [{v.quantile(.25):.1f}, {v.quantile(.75):.1f}]"


def add_continuous(label: str, col: str):
    if col not in cohort.columns:
        return
    a = pd.to_numeric(low[col], errors="coerce").dropna()
    b = pd.to_numeric(rest[col], errors="coerce").dropna()
    p = np.nan
    if len(a) >= 2 and len(b) >= 2:
        try:
            p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            pass
    for gname, gdf in GROUPS:
        v = pd.to_numeric(gdf[col], errors="coerce").dropna()
        long_rows.append({
            "variable": label, "level": "", "subgroup": gname, "stat_type": "continuous",
            "n": len(v), "total": len(gdf),
            "median": float(v.median()) if len(v) else None,
            "q25": float(v.quantile(.25)) if len(v) else None,
            "q75": float(v.quantile(.75)) if len(v) else None,
            "pval": p,
        })
    table_rows.append([label, _fmt_mi(low[col]), _fmt_mi(rest[col]), p])


def add_categorical(label: str, col: str, levels=None):
    if col not in cohort.columns:
        return
    s_all = cohort[col]
    if levels is None:
        levels = [lvl for lvl in pd.Series(s_all.dropna().unique()).tolist()]
    # contingency for p (level x group)
    try:
        ct = pd.crosstab(s_all[grp.isin(["Very low (10-15)", "Rest"])],
                         grp[grp.isin(["Very low (10-15)", "Rest"])])
        if ct.shape == (2, 2):
            p = stats.fisher_exact(ct.values)[1]
        elif ct.size and (ct.values.sum() > 0):
            p = stats.chi2_contingency(ct.values)[1]
        else:
            p = np.nan
    except ValueError:
        p = np.nan
    table_rows.append([label, "", "", p])  # header row carries the p-value
    for lvl in levels:
        for gname, gdf in GROUPS:
            col_s = gdf[col]
            n_lvl = int((col_s == lvl).sum())
            tot = int(col_s.notna().sum())
            long_rows.append({
                "variable": label, "level": str(lvl), "subgroup": gname, "stat_type": "categorical",
                "n": n_lvl, "total": tot, "median": None, "q25": None, "q75": None, "pval": p,
            })
        lo_n = int((low[col] == lvl).sum()); lo_t = int(low[col].notna().sum())
        re_n = int((rest[col] == lvl).sum()); re_t = int(rest[col].notna().sum())
        table_rows.append([
            f"  {lvl}",
            f"{lo_n} ({100*lo_n/lo_t:.0f}%)" if lo_t else "—",
            f"{re_n} ({100*re_n/re_t:.0f}%)" if re_t else "—",
            np.nan,
        ])


# Confirm dose separation
add_continuous("CRRT dose (mL/kg/hr)", "crrt_dose_ml_kg_hr")
# Demographics
add_continuous("Age (years)", "age_at_admission")
add_categorical("Sex", "sex_category")
# Severity
add_continuous("SOFA total (baseline)", "sofa_total")
add_continuous("SOFA renal (baseline)", "sofa_renal")
add_continuous("SOFA respiratory (baseline)", "sofa_resp")
# Baseline labs
for col, label in [
    ("creatinine_baseline", "Creatinine (baseline)"),
    ("bun_baseline", "BUN (baseline)"),
    ("lactate_baseline", "Lactate (baseline)"),
    ("ph_arterial_baseline", "Arterial pH (baseline)"),
    ("bicarbonate_baseline", "Bicarbonate (baseline)"),
    ("potassium_baseline", "Potassium (baseline)"),
]:
    add_continuous(label, col)
# Support / context
add_continuous("NEE (mcg/kg/min, baseline)", "nee_baseline")
add_categorical("Sepsis within window", "sepsis_within_window")
add_categorical("CRRT modality", "crrt_mode_category")
# Outcomes
add_continuous("CRRT duration (days)", "duration_days")
add_continuous("ICU LOS (days)", "icu_los_days")
add_categorical("30-day mortality", "death_30d", levels=[1])
add_categorical("In-hospital mortality", "in_hosp_death", levels=[1])

long_df = pd.DataFrame(long_rows)
long_df["site"] = SITE_NAME
long_df.to_csv(OUT / f"{SITE_NAME}_low_dose_long.csv", index=False)

def _fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ""
    return "<0.001" if p < 0.001 else f"{p:.3f}"


table_df = pd.DataFrame(table_rows, columns=[
    "Variable", f"Very low 10-15 (N={len(low)})", f"Rest (N={len(rest)})", "p"])
table_df["p"] = table_df["p"].map(_fmt_p)
table_df.to_csv(OUT / f"{SITE_NAME}_low_dose_table.csv", index=False)

print(f"  wrote counts / long / table CSVs to {OUT}")
print(f"=== 06 complete in {time.time() - t_start:.1f}s ===")
