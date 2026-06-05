"""
07_crrt_descriptive_epi.py — CRRT descriptive epidemiology + practice variation.

Per-site descriptive companion to the dose causal analysis (04-06b). Produces the
"epidemiology of CRRT across the consortium" worked example for the reframed,
Rojas-style manuscript (value of CLIF to ICU nephrologists). Emits AGGREGATE
per-site CSVs only — no patient-level data leaves the site. Multi-site pooling
and rendering happen in 09 (dashboard) / 10 (manuscript).

Two independent products:

  A. CRRT INCIDENCE / UTILIZATION (hospitalization level)
     Denominator = adult ICU hospitalizations (age>=18, admission-year filter,
     any ICU ADT stay). Numerator = raw CRRT receipt (any crrt_therapy record).
     Stratified: overall, among ventilated (any IMV), among vasopressor-exposed,
     and by ICU location_type (subtype). This is a UTILIZATION rate and is
     deliberately distinct from the post-exclusion analytic cohort.
        -> output/final/crrt_epi/{SITE}_crrt_incidence.csv

  B. PRACTICE VARIATION + TIER-A QUALITY (analytic CRRT cohort)
     From the analytic cohort (index_crrt_df + tableone_analysis_df): delivered
     dose vs KDIGO band, net ultrafiltration rate vs Murugan band, modality mix,
     blood flow, time-to-initiation, acid-base correction, duration, mortality.
        -> output/final/crrt_epi/{SITE}_crrt_practice_quality.csv  (long format)

Quality-metric scope is Tier A (computable from data the pipeline already loads).
Anticoagulation, delivered:prescribed ratio, downtime, and filter life are NOT in
CLIF and are out of scope. RRT-dependence-at-discharge / renal recovery are
deferred (need discharge-time RRT/creatinine not in the analytic frame).

Run standalone (from repo root):
    uv run python code/07_crrt_descriptive_epi.py
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline_helpers import validate_config, safe_load_clif_table

warnings.filterwarnings("ignore")

# ── Config (direct load, mirrors 00_cohort.py; cwd is code/) ────────────────
with open("../config/config.json", "r") as f:
    config = json.load(f)
config = validate_config(config)

SITE_NAME = config["site_name"]
TABLES_PATH = config["tables_path"]
HAS_CRRT_SETTINGS = config.get("has_crrt_settings", True)
YEAR_START = config.get("admission_year_start", 2018)
YEAR_END = config.get("admission_year_end", None)

INTER = Path("../output/intermediate")
OUT = Path("../output/final/crrt_epi")
OUT.mkdir(parents=True, exist_ok=True)

# Pure vasopressors (exclude inotropes dobutamine/milrinone/isoproterenol) for the
# "on vasopressors" denominator stratum.
VASO_CATEGORIES = [
    "norepinephrine", "epinephrine", "phenylephrine",
    "vasopressin", "dopamine", "angiotensin",
]

t_start = time.time()
print(f"=== 07 CRRT descriptive epidemiology — {SITE_NAME} ===")


# ── Part A: adult-ICU denominator (hospitalization level) ───────────────────
def build_incidence() -> pd.DataFrame:
    from clifpy.clif_orchestrator import ClifOrchestrator

    clif = ClifOrchestrator(
        data_directory=TABLES_PATH,
        filetype=config["file_type"],
        timezone=config["timezone"],
    )

    # Core tables, full population
    for tbl in ["hospitalization", "adt"]:
        safe_load_clif_table(clif, tbl, tables_path=TABLES_PATH)
    hosp = clif.hospitalization.df
    adt = clif.adt.df

    # Adult + admission-year filter -> denominator hospitalization set
    adult = hosp[(hosp["age_at_admission"] >= 18) & hosp["age_at_admission"].notna()].copy()
    yr = adult["admission_dttm"].dt.year
    mask = yr >= YEAR_START
    if YEAR_END is not None:
        mask &= yr <= YEAR_END
    adult = adult[mask]
    adult_ids = set(adult["hospitalization_id"].unique())
    print(f"  adult hospitalizations (>=18, yr>={YEAR_START}): {len(adult_ids):,}")

    # ICU stays among adults; capture ICU location_types visited per hospitalization
    adt_adult = adt[adt["hospitalization_id"].isin(adult_ids)].copy()
    icu_rows = adt_adult[adt_adult["location_category"].astype("string").str.lower() == "icu"]
    icu_ids = set(icu_rows["hospitalization_id"].unique())
    print(f"  adult ICU hospitalizations (denominator): {len(icu_ids):,}")

    # location_type strata: hosp ids per ICU subtype (a hosp may appear in several)
    subtype_ids: dict[str, set] = {}
    for lt, grp in icu_rows.groupby(icu_rows["location_type"].astype("string").fillna("unknown")):
        subtype_ids[str(lt)] = set(grp["hospitalization_id"].unique())

    # Ventilated stratum: any IMV device among ICU hospitalizations
    safe_load_clif_table(
        clif, "respiratory_support", tables_path=TABLES_PATH,
        columns=["hospitalization_id", "recorded_dttm", "device_category"],
        filters={"hospitalization_id": list(icu_ids)},
    )
    rs = clif.respiratory_support.df
    vent_ids = set(
        rs.loc[rs["device_category"].astype("string").str.lower() == "imv", "hospitalization_id"].unique()
    )
    print(f"  ICU hospitalizations ventilated (IMV): {len(vent_ids):,}")

    # Vasopressor stratum: any continuous vasopressor among ICU hospitalizations
    safe_load_clif_table(
        clif, "medication_admin_continuous", tables_path=TABLES_PATH,
        columns=["hospitalization_id", "admin_dttm", "med_category"],
        filters={"hospitalization_id": list(icu_ids), "med_category": VASO_CATEGORIES},
    )
    mc = clif.medication_admin_continuous.df
    vaso_ids = set(mc["hospitalization_id"].unique())
    print(f"  ICU hospitalizations on vasopressors: {len(vaso_ids):,}")

    # Numerator: raw CRRT receipt (any crrt_therapy record) among ICU hospitalizations
    safe_load_clif_table(
        clif, "crrt_therapy", tables_path=TABLES_PATH,
        filters={"hospitalization_id": list(icu_ids)},
    )
    crrt = clif.crrt_therapy.df
    crrt_ids = set(crrt["hospitalization_id"].unique())
    print(f"  ICU hospitalizations with CRRT receipt (numerator): {len(crrt_ids):,}")

    # Assemble incidence rows
    def row(stratum: str, denom: set) -> dict:
        n_d = len(denom)
        n_c = len(denom & crrt_ids)
        return {
            "stratum": stratum,
            "n_denominator": n_d,
            "n_crrt": n_c,
            "incidence_pct": round(100 * n_c / n_d, 2) if n_d else np.nan,
        }

    rows = [
        row("Overall (adult ICU)", icu_ids),
        row("On invasive ventilation", vent_ids & icu_ids),
        row("On vasopressors", vaso_ids & icu_ids),
    ]
    for lt in sorted(subtype_ids):
        rows.append(row(f"ICU subtype: {lt}", subtype_ids[lt]))

    df = pd.DataFrame(rows)
    df.insert(0, "site", SITE_NAME)
    return df


# ── Part B: practice variation + Tier-A quality (analytic cohort) ───────────
def _q(series: pd.Series) -> tuple[float, float, float, int]:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return (np.nan, np.nan, np.nan, 0)
    return (s.median(), s.quantile(0.25), s.quantile(0.75), len(s))


def _long(rows: list, variable, stat, value, n=None, denom=None):
    rows.append({"site": SITE_NAME, "variable": variable, "stat": stat,
                 "value": value, "n": n, "denominator": denom})


def build_practice_quality() -> pd.DataFrame:
    idx = pd.read_parquet(INTER / "index_crrt_df.parquet")
    t1 = pd.read_parquet(INTER / "tableone_analysis_df.parquet")
    # Analytic cohort = encounter_blocks in tableone frame; bring weight_kg from index
    w = idx[["encounter_block", "weight_kg"]].drop_duplicates("encounter_block")
    coh = t1.merge(w, on="encounter_block", how="left")
    n_cohort = len(coh)
    print(f"  analytic CRRT cohort: {n_cohort:,} encounters")

    rows: list = []
    _long(rows, "Analytic cohort size", "n", n_cohort)

    if HAS_CRRT_SETTINGS:
        # Delivered dose (median-first-3h) + KDIGO bands
        dose = pd.to_numeric(coh["crrt_dose_ml_kg_hr"], errors="coerce")
        med, q1, q3, n = _q(dose)
        _long(rows, "Delivered dose (mL/kg/hr)", "median_iqr", f"{med:.1f} [{q1:.1f}, {q3:.1f}]", n)
        d = dose.dropna()
        bands = {
            "<15 (very low)": (d < 15),
            "15-20": (d >= 15) & (d < 20),
            "20-25 (KDIGO)": (d >= 20) & (d <= 25),
            ">25": (d > 25),
        }
        for label, m in bands.items():
            _long(rows, "Dose band", f"pct_{label}", round(100 * m.sum() / len(d), 1) if len(d) else np.nan,
                  int(m.sum()), len(d))

        # Net ultrafiltration rate = ultrafiltration_out / weight  (mL/kg/hr)
        uf = pd.to_numeric(coh["ultrafiltration_out"], errors="coerce")
        wt = pd.to_numeric(coh["weight_kg"], errors="coerce")
        uf_rate = uf / wt
        uf_rate = uf_rate.where(uf_rate.between(0, 12))  # drop implausible (>~1000 mL/hr)
        med, q1, q3, n = _q(uf_rate)
        _long(rows, "Net UF rate (mL/kg/hr)", "median_iqr", f"{med:.2f} [{q1:.2f}, {q3:.2f}]", n)
        r = uf_rate.dropna()
        in_band = (r >= 1.0) & (r <= 1.75)
        _long(rows, "Net UF rate", "pct_in_1.0-1.75_band",
              round(100 * in_band.sum() / len(r), 1) if len(r) else np.nan, int(in_band.sum()), len(r))

        # Modality mix
        modes = coh["crrt_mode_category"].astype("string").str.upper().fillna("UNKNOWN")
        for mode, cnt in modes.value_counts().items():
            _long(rows, "Modality", f"pct_{mode}", round(100 * cnt / len(modes), 1), int(cnt), len(modes))

        # Blood flow rate
        med, q1, q3, n = _q(coh["blood_flow_rate"])
        _long(rows, "Blood flow rate (mL/min)", "median_iqr", f"{med:.0f} [{q1:.0f}, {q3:.0f}]", n)

        # Acid-base delivery adequacy: among acidotic at baseline, % corrected by post window
        ph_b = pd.to_numeric(coh["ph_arterial_baseline"], errors="coerce")
        hco3_b = pd.to_numeric(coh["bicarbonate_baseline"], errors="coerce")
        ph_p = pd.to_numeric(coh["ph_arterial_post_crrt"], errors="coerce")
        hco3_p = pd.to_numeric(coh["bicarbonate_post_crrt"], errors="coerce")
        acidotic = (ph_b < 7.30) | (hco3_b < 22)
        corrected = ((ph_p >= 7.30) | (hco3_p >= 22)) & acidotic
        n_ac = int(acidotic.sum())
        _long(rows, "Acidosis corrected (baseline->post)", "pct",
              round(100 * corrected.sum() / n_ac, 1) if n_ac else np.nan, int(corrected.sum()), n_ac)

    # Time to CRRT initiation from first vital (proxy for presentation->CRRT), hours
    init = pd.to_datetime(coh["crrt_initiation_time"], errors="coerce")
    fv = pd.to_datetime(coh["first_vital_dttm"], errors="coerce")
    tti = (init - fv).dt.total_seconds() / 3600.0
    tti = tti.where(tti >= 0)
    med, q1, q3, n = _q(tti)
    _long(rows, "Time to CRRT from first vital (h)", "median_iqr", f"{med:.1f} [{q1:.1f}, {q3:.1f}]", n)

    # CRRT duration (days)
    med, q1, q3, n = _q(coh["duration_days"])
    _long(rows, "CRRT duration (days)", "median_iqr", f"{med:.1f} [{q1:.1f}, {q3:.1f}]", n)

    # Outcomes
    for col, label in [("death_30d", "30-day mortality"), ("in_hosp_death", "In-hospital mortality")]:
        s = pd.to_numeric(coh[col], errors="coerce").dropna()
        _long(rows, label, "pct", round(100 * s.mean(), 1) if len(s) else np.nan, int(s.sum()), len(s))

    return pd.DataFrame(rows)


# ── Run ─────────────────────────────────────────────────────────────────────
print("\n[A] CRRT incidence / utilization ...")
inc = build_incidence()
inc_path = OUT / f"{SITE_NAME}_crrt_incidence.csv"
inc.to_csv(inc_path, index=False)
print(f"  wrote {inc_path}  ({len(inc)} strata)")

print("\n[B] Practice variation + Tier-A quality ...")
pq = build_practice_quality()
pq_path = OUT / f"{SITE_NAME}_crrt_practice_quality.csv"
pq.to_csv(pq_path, index=False)
print(f"  wrote {pq_path}  ({len(pq)} metric rows)")

print(f"\n=== 07 complete in {time.time() - t_start:.1f}s ===")
