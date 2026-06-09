"""
03_crrt_epidemiology.py — comprehensive CRRT descriptive epidemiology.

Per-site descriptive epidemiology for the Rojas-style manuscript (value of CLIF
to ICU nephrologists). Consolidates the former 07 (incidence + practice/quality
+ distribution figures) with the reworked longitudinal trajectories from the
former 03_crrt_visualizations.py. Emits AGGREGATE per-site CSVs/figures only —
no patient-level data leaves the site. Multi-site pooling and rendering happen
in 09 (dashboard) / 10 (manuscript).

Products:

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

  C. DISTRIBUTION FIGURES (analytic cohort): incidence by context / ICU subtype,
     CRRT dose distribution, net-UF distribution (Murugan groups), and crude
     30-day mortality vs net UF.

  D. LONGITUDINAL TRAJECTORIES over 30 days (x-axis in days), each carrying a
     "% alive and on CRRT" census line: CRRT dose, lab trajectories, vasopressor
     (NEE), and patient state (On IMV / Off IMV / Discharged / Dead). MAP and
     respiratory (FiO2/IMV) trajectories were cut.

NOTE (t=0 / limitation): t=0 is each encounter's CRRT initiation time. ~2% of
encounters have a recorded death at or before t=0 (death_dttm <= crrt start) --
a death/CRRT-start timestamp ordering artifact -- so they appear "Dead" at day 0
in BOTH state-over-time figures (patient-state and CRRT-state), which therefore
start near 98% rather than 100%. Left as-is to show the data faithfully; flag in
the manuscript Limitations.

Quality-metric scope is Tier A (computable from data the pipeline already loads).
Anticoagulation, delivered:prescribed ratio, downtime, and filter life are NOT in
CLIF and are out of scope. RRT-dependence-at-discharge / renal recovery are
deferred (need discharge-time RRT/creatinine not in the analytic frame).

Run standalone (from repo root):
    uv run python code/03_crrt_epidemiology.py
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Shared manuscript visual language (matches 09/10): Arial sans-serif + size.
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 13,
})

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
GRAPHS = OUT / "graphs"
GRAPHS.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, GREEN = "#1e417c", "#fb801b", "#4CAF50"  # ASN palette + KDIGO-band green

# Pure vasopressors (exclude inotropes dobutamine/milrinone/isoproterenol) for the
# "on vasopressors" denominator stratum.
VASO_CATEGORIES = [
    "norepinephrine", "epinephrine", "phenylephrine",
    "vasopressin", "dopamine", "angiotensin",
]

t_start = time.time()
print(f"=== 03 CRRT epidemiology — {SITE_NAME} ===")


# ── Part A: adult-ICU denominator (hospitalization level) ───────────────────
def build_incidence() -> pd.DataFrame:
    from clifpy.clif_orchestrator import ClifOrchestrator

    clif = ClifOrchestrator(
        data_directory=TABLES_PATH,
        filetype=config["file_type"],
        timezone=config["timezone"],
        output_directory="../output",  # else clifpy logs to <cwd>/output (= code/output)
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
        # CRRT dose (median-first-3h; prescribed machine dose) + KDIGO bands
        dose = pd.to_numeric(coh["crrt_dose_ml_kg_hr"], errors="coerce")
        med, q1, q3, n = _q(dose)
        _long(rows, "CRRT dose (mL/kg/hr)", "median_iqr", f"{med:.1f} [{q1:.1f}, {q3:.1f}]", n)
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


# ── Part C: per-site figures (output/final/crrt_epi/graphs/) ────────────────
def _pretty_icu(stratum: str) -> str:
    """ICU subtype stratum label -> human-readable (medical_icu -> Medical ICU)."""
    s = stratum.replace("ICU subtype: ", "").replace("_", " ").strip().title()
    return s.replace("Icu", "ICU")


def _barh_incidence(d: pd.DataFrame, title: str, fname: str, prettify: bool) -> None:
    """Horizontal incidence bars with %, n/N labels. One homogeneous color per
    figure (color no longer encodes a grouping the reader cannot infer)."""
    if d.empty:
        return
    d = d.sort_values("incidence_pct")
    labels = [_pretty_icu(s) for s in d["stratum"]] if prettify else list(d["stratum"])
    xmax = float(d["incidence_pct"].max()) if len(d) else 1.0
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(labels, d["incidence_pct"], color=BLUE)
    for y, (v, n, N) in enumerate(zip(d["incidence_pct"], d["n_crrt"], d["n_denominator"])):
        ax.text(v + xmax * 0.01, y, f"{v:.1f}% ({int(n):,}/{int(N):,})",
                va="center", fontsize=9)
    ax.set_xlabel("% of Hospitalizations Receiving CRRT at Any Point")
    ax.set_title(title)
    # Headroom so the "% (n/N)" labels never collide with the right spine.
    ax.set_xlim(0, xmax * 1.42)
    fig.tight_layout()
    fig.savefig(GRAPHS / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)


def build_figures(inc_df: pd.DataFrame) -> None:
    idx = pd.read_parquet(INTER / "index_crrt_df.parquet")
    t1 = pd.read_parquet(INTER / "tableone_analysis_df.parquet")

    # (1) CRRT incidence — split into two figures so a single homogeneous color
    #     suffices in each (the prior single figure used orange/navy to encode a
    #     grouping the reader could not infer):
    #       (1a) by clinical context: overall + co-indications (IMV, vasopressors)
    #       (1b) by ICU subtype: one bar per location_type (labels prettified)
    context = {"Overall (adult ICU)", "On invasive ventilation", "On vasopressors"}
    ctx = inc_df[inc_df["stratum"].isin(context)].copy()
    _barh_incidence(ctx, f"CRRT Incidence by Clinical Context: {SITE_NAME}",
                    f"{SITE_NAME}_crrt_incidence_by_context.png", prettify=False)
    sub = inc_df[inc_df["stratum"].str.startswith("ICU subtype: ")].copy()
    _barh_incidence(sub, f"CRRT Incidence by ICU Subtype: {SITE_NAME}",
                    f"{SITE_NAME}_crrt_incidence_by_icu_subtype.png", prettify=True)

    if not HAS_CRRT_SETTINGS:
        return

    # (2) Initial CRRT dose distribution (initial dose = median of the first 3 h).
    #     The 20-30 band spans prescribed targets achieving KDIGO 20-25 delivered.
    #     Median computed on the FULL vector; the [0, 80] clip governs display only.
    dose_all = pd.to_numeric(t1["crrt_dose_ml_kg_hr"], errors="coerce").dropna()
    if len(dose_all):
        med = dose_all.median()
        pct_vlow = 100 * dose_all.between(10, 15).sum() / len(dose_all)
        pct_kdigo = 100 * dose_all.between(20, 30).sum() / len(dose_all)
        dose_disp = dose_all[dose_all.between(0, 80)]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(dose_disp, bins=40, color=BLUE, alpha=0.85,
                edgecolor="black", linewidth=0.5)
        ax.axvspan(10, 15, color=ORANGE, alpha=0.30,
                   label=f"Very Low (10–15 mL/kg/hr): {pct_vlow:.1f}%")
        ax.axvspan(20, 30, color=GREEN, alpha=0.30,
                   label=f"KDIGO Recommendation (20–30 mL/kg/hr): {pct_kdigo:.1f}%")
        ax.axvline(med, color="black", linestyle="--", linewidth=1.2,
                   label=f"Median {med:.1f} mL/kg/hr")
        ax.set_xlabel("Initial CRRT Dose (mL/kg/hr)")
        ax.set_ylabel("Encounters")
        ax.set_title(f"Initial CRRT Dose Distribution (Median First 3 h): {SITE_NAME}")
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(GRAPHS / f"{SITE_NAME}_dose_distribution.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # Per-encounter net UF rate (mL/kg/hr), shared by figures (3) and (4).
    w = idx[["encounter_block", "weight_kg"]].drop_duplicates("encounter_block")
    m = t1[["encounter_block", "ultrafiltration_out", "death_30d"]].merge(
        w, on="encounter_block", how="left")
    uf_rate = (pd.to_numeric(m["ultrafiltration_out"], errors="coerce")
               / pd.to_numeric(m["weight_kg"], errors="coerce"))

    # (3) Net ultrafiltration rate distribution with the three Murugan 2019 groups
    #     (low <1.01, middle 1.01-1.75, high >1.75). All percentages are of the
    #     total analytic cohort (incl. the UF=0 subset, which sits within Low).
    uf = uf_rate[uf_rate.between(0, 6)].dropna()
    if len(uf):
        n = len(uf)
        pct_low = 100 * (uf < 1.01).sum() / n
        pct_mid = 100 * uf.between(1.01, 1.75).sum() / n
        pct_high = 100 * (uf > 1.75).sum() / n
        pct_zero = 100 * (uf == 0).sum() / n
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(uf, bins=40, color=BLUE, alpha=0.85, edgecolor="black", linewidth=0.5)
        ax.axvspan(1.01, 1.75, color=GREEN, alpha=0.30, label="Middle (1.01–1.75)")
        ax.axvline(1.01, color="#555555", linestyle=":", linewidth=1.0)
        ax.axvline(1.75, color="#555555", linestyle=":", linewidth=1.0)
        ax.axvline(uf.median(), color="black", linestyle="--", linewidth=1.2,
                   label=f"Median {uf.median():.2f} mL/kg/hr")
        ax.set_xlabel("Net Ultrafiltration Rate (mL/kg/hr)")
        ax.set_ylabel("Encounters")
        ax.set_title(f"Net Ultrafiltration Rate Distribution: {SITE_NAME}")
        ax.legend(fontsize=9, loc="center right")
        txt = (f"Murugan 2019 groups (% of all encounters)\n"
               f"Low (<1.01): {pct_low:.1f}%\n"
               f"Middle (1.01–1.75): {pct_mid:.1f}%\n"
               f"High (>1.75): {pct_high:.1f}%\n"
               f"Zero net UF: {pct_zero:.1f}%")
        ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top",
                fontsize=9, bbox=dict(boxstyle="round", facecolor="white",
                                      alpha=0.85, edgecolor="#cccccc"))
        fig.tight_layout()
        fig.savefig(GRAPHS / f"{SITE_NAME}_net_uf_distribution.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # (4) Crude 30-day mortality vs net UF rate (Murugan 2019 Fig e2 style; binned,
    #     marker size ∝ n). Descriptive only — shows the nonlinear (J-shaped)
    #     association; deliberately NOT the Gray time-varying-coefficient survival
    #     model used in that paper (beyond scope here).
    mm = pd.DataFrame({"uf": uf_rate,
                       "death": pd.to_numeric(m["death_30d"], errors="coerce")})
    mm = mm[mm["uf"].between(0, 4)].dropna()
    if len(mm) >= 50:
        edges = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 4.0]
        mm["bin"] = pd.cut(mm["uf"], bins=edges, include_lowest=True)
        g = mm.groupby("bin", observed=True)["death"].agg(["mean", "size"])
        g = g[g["size"] >= 15]
        if len(g) >= 3:
            x = [iv.mid for iv in g.index]
            y = (g["mean"] * 100).to_numpy()
            smax = float(g["size"].max())
            sizes = (30 + 200 * g["size"] / smax).to_numpy()
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.axvspan(1.01, 1.75, color=GREEN, alpha=0.20, label="Murugan middle (1.01–1.75)")
            ax.plot(x, y, color=BLUE, linewidth=1.5, zorder=2)
            ax.scatter(x, y, s=sizes, color=BLUE, alpha=0.85, zorder=3,
                       label="Crude 30-day mortality (marker size scaled to n)")
            ax.set_xlabel("Net Ultrafiltration Rate (mL/kg/hr)")
            ax.set_ylabel("Crude 30-Day Mortality (%)")
            ax.set_title(f"Crude 30-Day Mortality vs Net Ultrafiltration Rate: {SITE_NAME}")
            ax.legend(fontsize=9)
            fig.tight_layout()
            fig.savefig(GRAPHS / f"{SITE_NAME}_uf_mortality.png", dpi=150, bbox_inches="tight")
            plt.close(fig)


# ── Part D: longitudinal trajectories over the CRRT course ──────────────────
# Reworked from the legacy 03_crrt_visualizations.py. Each trajectory is shown
# over 30 days (720 h) with a 3-day (72 h) inset zoom, since most change is
# early. MAP and respiratory (FiO2/IMV) trajectories were cut. Shared ASN palette.
import pyarrow.parquet as _pq

MAX_HOURS_TRAJ = 720   # 30 days
INSET_HOURS = 72       # 3-day zoom
BIN_H = 12             # 12-hour bins for labs / vasopressor
GREY = "#555555"


def _crrt_dose_per_row(df: pd.DataFrame) -> np.ndarray:
    """Mode-specific total CRRT flow / weight -> dose (mL/kg/hr). Mirrors the
    formula in 00_cohort.py (kept the single shared definition for trajectories)."""
    mode = df["crrt_mode_category"]
    dia = pd.to_numeric(df.get("crrt_dialysate_flow_rate"), errors="coerce")
    pre = pd.to_numeric(df.get("crrt_pre_filter_replacement_fluid_rate"), errors="coerce")
    post = pd.to_numeric(df.get("crrt_post_filter_replacement_fluid_rate"), errors="coerce")
    total = np.select(
        [mode == "cvvhd", mode == "cvvh", mode == "cvvhdf"],
        [dia, pre.fillna(0) + post.fillna(0), dia.fillna(0) + pre.fillna(0) + post.fillna(0)],
        default=np.nan,
    )
    wt = pd.to_numeric(df["weight_kg"], errors="coerce")
    return np.where((wt > 0) & (total > 0), total / wt, np.nan)


def _load_traj_wide() -> pd.DataFrame:
    """wide_df (needed columns only) joined to encounter_block + crrt_initiation
    time + weight, with hours_from_crrt computed."""
    crrt_init = pd.read_parquet(INTER / "crrt_initiation.parquet")
    idx = pd.read_parquet(INTER / "index_crrt_df.parquet")
    cohort = pd.read_parquet(INTER / "cohort_df.parquet",
                             columns=["hospitalization_id", "encounter_block"])
    crrt_cols = (["crrt_dialysate_flow_rate", "crrt_pre_filter_replacement_fluid_rate",
                  "crrt_post_filter_replacement_fluid_rate", "crrt_mode_category"]
                 if HAS_CRRT_SETTINGS else [])
    lab_cols = ["lab_lactate", "lab_ph_arterial", "lab_bicarbonate", "lab_potassium", "lab_phosphate"]
    other = ["med_cont_nee", "resp_device_category"]  # resp_device_category drives patient-state IMV
    avail = {f.name for f in _pq.read_schema(INTER / "wide_df.parquet")}
    need = ["hospitalization_id", "event_dttm"] + [c for c in crrt_cols + lab_cols + other if c in avail]
    w = pd.read_parquet(INTER / "wide_df.parquet", columns=need)
    w = w.merge(cohort, on="hospitalization_id", how="inner")
    w = w.merge(crrt_init[["encounter_block", "crrt_initiation_time"]], on="encounter_block", how="inner")
    w = w.merge(idx[["encounter_block", "weight_kg"]].drop_duplicates("encounter_block"),
                on="encounter_block", how="left")
    w["hours_from_crrt"] = (w["event_dttm"] - w["crrt_initiation_time"]).dt.total_seconds() / 3600.0
    return w


def _median_iqr_by_bin(df: pd.DataFrame, col: str, bin_h: int = BIN_H,
                       max_h: int = MAX_HOURS_TRAJ) -> pd.DataFrame:
    d = df[(df["hours_from_crrt"] >= 0) & (df["hours_from_crrt"] <= max_h)].copy()
    d["_v"] = pd.to_numeric(d[col], errors="coerce")
    d = d.dropna(subset=["_v"])
    if d.empty:
        return pd.DataFrame(columns=["hour", "median", "q25", "q75", "n", "site"])
    d["hour"] = (d["hours_from_crrt"] // bin_h).astype(int) * bin_h + bin_h / 2.0
    g = (d.groupby("hour")["_v"]
         .agg(median="median", q25=lambda x: x.quantile(0.25),
              q75=lambda x: x.quantile(0.75), n="count")
         .reset_index())
    g["site"] = SITE_NAME
    return g


MAX_DAYS_TRAJ = MAX_HOURS_TRAJ / 24.0  # 30 days (x-axis is in days)


def _crrt_census(w: pd.DataFrame) -> pd.DataFrame | None:
    """Per integer day, the proportion of the baseline cohort with a CRRT record
    that day (i.e., still alive and on CRRT). Denominator = distinct encounter
    blocks in the cohort. Used as a number-at-risk style context line."""
    if "crrt_mode_category" not in w.columns:
        return None
    n_base = w["encounter_block"].nunique()
    c = w[(w["hours_from_crrt"] >= 0) & (w["hours_from_crrt"] <= MAX_HOURS_TRAJ)
          & w["crrt_mode_category"].notna()].copy()
    if c.empty or not n_base:
        return None
    c["day"] = (c["hours_from_crrt"] // 24).astype(int)
    per = c.groupby("day")["encounter_block"].nunique().reset_index(name="n")
    per["pct"] = per["n"] / n_base * 100.0
    return per


def _add_census(ax, census, ylabel=True):
    """Overlay the % alive-and-on-CRRT census on a secondary (0-100%) axis.
    Returns the line handle (for a combined legend)."""
    if census is None or census.empty:
        return None
    ax2 = ax.twinx()
    (ln,) = ax2.plot(census["day"], census["pct"], color=GREY, ls="--", lw=1.4, alpha=0.85)
    ax2.set_ylim(0, 100)
    ax2.tick_params(axis="y", labelcolor=GREY, labelsize=8)
    ax2.set_ylabel("% of Cohort Alive and on CRRT" if ylabel else "", color=GREY, fontsize=9)
    if not ylabel:
        ax2.set_yticklabels([])
    return ln


def _at_risk_row(ax, census, days=(0, 5, 10, 15, 20, 25, 30)):
    """KM-style number-at-risk row beneath the x-axis: N (and %) alive & on CRRT."""
    if census is None or census.empty:
        return
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


def _draw_traj(ax, g, color, band=None, band_label=None, line_label=None):
    if band is not None:
        ax.axhspan(band[0], band[1], color=GREEN, alpha=0.18, label=band_label)
    x = g["hour"] / 24.0  # hours -> days
    ax.plot(x, g["median"], color=color, lw=1.8, zorder=3, label=line_label)
    ax.fill_between(x, g["q25"], g["q75"], color=color, alpha=0.18, zorder=2)
    ax.grid(axis="y", alpha=0.25)


def _traj_axis(ax, g, ylabel, color, band=None, band_label=None, line_label=None,
               census=None, census_ylabel=True, legend=False, atrisk=False):
    """Draw a 30-day trajectory (x in days) on ax, with optional census line,
    combined legend (incl. the grey census line), and number-at-risk row."""
    _draw_traj(ax, g, color, band=band, band_label=band_label, line_label=line_label)
    ax.set_xlim(0, MAX_DAYS_TRAJ)
    ax.set_xlabel("Days from CRRT Initiation")
    ax.set_ylabel(ylabel)
    cens = _add_census(ax, census, ylabel=census_ylabel) if census is not None else None
    if legend:
        h, lab = ax.get_legend_handles_labels()
        if cens is not None:
            h.append(cens)
            lab.append("Alive and on CRRT (%, Right Axis)")
        if h:
            ax.legend(h, lab, fontsize=8, loc="upper right", framealpha=0.9)
    if atrisk:
        _at_risk_row(ax, census)


def build_trajectories(w: pd.DataFrame) -> None:
    # Number-at-risk style context line shared across the "over time" figures:
    # proportion of the cohort still alive and on CRRT each day.
    census = _crrt_census(w)

    # (D1) CRRT dose over time
    if HAS_CRRT_SETTINGS and "crrt_mode_category" in w.columns:
        dr = w[(w["hours_from_crrt"] >= 0) & (w["crrt_mode_category"].notna())].copy()
        dr["dose"] = _crrt_dose_per_row(dr)
        g = _median_iqr_by_bin(dr, "dose", bin_h=6)
        if not g.empty:
            g.to_csv(GRAPHS / f"{SITE_NAME}_crrt_dose_hourly.csv", index=False)
            fig, ax = plt.subplots(figsize=(11, 5.5))
            _traj_axis(ax, g, "CRRT Dose (mL/kg/hr)", BLUE, band=(20, 30),
                       band_label="KDIGO Recommendation (20–30)", line_label="Median CRRT Dose",
                       census=census, legend=True, atrisk=True)
            ax.set_title(f"CRRT Dose Over 30 Days: {SITE_NAME}")
            fig.tight_layout(rect=[0, 0.08, 1, 1])
            fig.savefig(GRAPHS / f"{SITE_NAME}_crrt_dose_over_time.png", dpi=150, bbox_inches="tight")
            plt.close(fig)

    # (D2) Lab trajectories (5-panel grid)
    lab_info = [
        ("lab_lactate", "Lactate (mmol/L)"),
        ("lab_ph_arterial", "Arterial pH"),
        ("lab_bicarbonate", "Bicarbonate (mEq/L)"),
        ("lab_potassium", "Potassium (mEq/L)"),
        ("lab_phosphate", "Phosphate (mg/dL)"),
    ]
    present = [(c, lbl) for c, lbl in lab_info if c in w.columns]
    if present:
        all_lab = []
        fig, axes = plt.subplots(3, 2, figsize=(14, 13))
        axes = axes.flatten()
        for i, (col, lbl) in enumerate(present):
            g = _median_iqr_by_bin(w, col)
            if g.empty:
                axes[i].set_visible(False)
                continue
            all_lab.append(g.assign(lab=col))
            _traj_axis(axes[i], g, lbl, GREEN, census=census, census_ylabel=False)
            axes[i].set_title(lbl)
        # 6th panel = census reference (legend for the dashed grey line), with at-risk row
        axc = axes[5]
        if census is not None and not census.empty:
            axc.plot(census["day"], census["pct"], color=GREY, ls="--", lw=1.8)
            axc.set_xlim(0, MAX_DAYS_TRAJ)
            axc.set_ylim(0, 100)
            axc.set_xlabel("Days from CRRT Initiation")
            axc.set_ylabel("% of Cohort")
            axc.set_title("Alive and on CRRT (%)")
            axc.grid(axis="y", alpha=0.25)
            _at_risk_row(axc, census)
        else:
            axc.set_visible(False)
        for j in range(len(present), 5):  # hide any gap panels (keep slot 5 = census)
            axes[j].set_visible(False)
        if all_lab:
            pd.concat(all_lab, ignore_index=True).to_csv(
                GRAPHS / f"{SITE_NAME}_lab_distributions_over_crrt.csv", index=False)
        fig.suptitle(f"Lab Trajectories Over 30 Days Post-CRRT: {SITE_NAME}\n"
                     "(Dashed Grey: % of Cohort Alive and on CRRT; See Bottom-Right Panel)",
                     fontsize=13, y=1.0)
        fig.tight_layout()
        fig.savefig(GRAPHS / f"{SITE_NAME}_lab_distributions_over_crrt.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # (D3) Vasopressor (norepinephrine-equivalent) trajectory
    if "med_cont_nee" in w.columns:
        g = _median_iqr_by_bin(w, "med_cont_nee")
        if not g.empty:
            g.to_csv(GRAPHS / f"{SITE_NAME}_nee_over_crrt.csv", index=False)
            fig, ax = plt.subplots(figsize=(11, 5.5))
            _traj_axis(ax, g, "Norepinephrine Equivalents (NEE, mcg/kg/min)", ORANGE,
                       line_label="Median NEE", census=census, legend=True, atrisk=True)
            ax.set_title(f"Vasopressor (Norepinephrine-Equivalent) Over 30 Days: {SITE_NAME}")
            fig.tight_layout(rect=[0, 0.08, 1, 1])
            fig.savefig(GRAPHS / f"{SITE_NAME}_nee_over_crrt.png", dpi=150, bbox_inches="tight")
            plt.close(fig)


def build_patient_course(w: pd.DataFrame) -> None:
    """(D4) Invasive mechanical ventilation (IMV) state over the CRRT course
    (On IMV / Off IMV / Discharged alive / Dead) as a stacked area over 30 days,
    + a time-to-event CSV. Ported from legacy 03; the outcome-trajectory
    companion to the practice/quality duration & mortality summaries (Part B)."""
    import matplotlib.colors as mcolors
    MAXH = MAX_HOURS_TRAJ
    crrt_init = pd.read_parquet(INTER / "crrt_initiation.parquet")
    oc = pd.read_parquet(INTER / "outcomes_df.parquet")
    cols = ["encounter_block", "in_hosp_death", "died", "death_dttm", "last_vital_dttm"]
    if "discharge_dttm" in oc.columns:
        cols.append("discharge_dttm")
    ot = oc[cols].merge(crrt_init[["encounter_block", "crrt_initiation_time"]],
                        on="encounter_block", how="inner")
    ot["hours_to_death"] = np.where(
        ot["died"] == 1,
        (ot["death_dttm"].fillna(ot["last_vital_dttm"]) - ot["crrt_initiation_time"]).dt.total_seconds() / 3600,
        np.nan)
    disc_col = "discharge_dttm" if "discharge_dttm" in ot.columns else "last_vital_dttm"
    ot["hours_to_discharge"] = np.where(
        ot["died"] == 0,
        (ot[disc_col] - ot["crrt_initiation_time"]).dt.total_seconds() / 3600,
        np.nan)

    # Hourly IMV status (last value per encounter per hour)
    if "resp_device_category" in w.columns:
        rr = w[w["resp_device_category"].notna()
               & (w["hours_from_crrt"] >= 0) & (w["hours_from_crrt"] <= MAXH)].copy()
        rr["hour"] = rr["hours_from_crrt"].astype(int)
        rr["is_imv"] = (rr["resp_device_category"].astype("string").str.lower() == "imv").astype(int)
        hourly_imv = (rr.sort_values(["encounter_block", "hours_from_crrt"])
                      .groupby(["encounter_block", "hour"])["is_imv"].last())
    else:
        hourly_imv = pd.Series(dtype=int)

    all_patients = crrt_init["encounter_block"].unique()
    n_pat = len(all_patients)
    dead_ebs = ot.loc[ot["hours_to_death"].notna(), ["encounter_block", "hours_to_death"]]
    disc_ebs = ot.loc[ot["hours_to_discharge"].notna(), ["encounter_block", "hours_to_discharge"]]

    rows = []
    for h in range(0, MAXH + 1):
        dead_by = set(dead_ebs.loc[dead_ebs["hours_to_death"] <= h, "encounter_block"])
        disc_by = set(disc_ebs.loc[disc_ebs["hours_to_discharge"] <= h, "encounter_block"])
        imv_at = (hourly_imv.xs(h, level="hour") if (len(hourly_imv) and h in hourly_imv.index.get_level_values("hour"))
                  else pd.Series(dtype=int))
        n_dead = n_disc = n_imv = n_off = 0
        for eb in all_patients:
            if eb in dead_by:
                n_dead += 1
            elif eb in disc_by:
                n_disc += 1
            elif eb in imv_at.index:
                if imv_at.loc[eb] == 1:
                    n_imv += 1
                else:
                    n_off += 1
            else:
                n_off += 1
        rows.append({"hour": h, "n_total": n_pat, "n_dead": n_dead,
                     "n_discharged": n_disc, "n_imv": n_imv, "n_off_imv": n_off})
    sdf = pd.DataFrame(rows)
    for k in ["imv", "off_imv", "discharged", "dead"]:
        sdf[f"prop_{k}"] = sdf[f"n_{k}"] / sdf["n_total"] * 100
    sdf["site"] = SITE_NAME
    sdf.to_csv(GRAPHS / f"{SITE_NAME}_imv_state_over_crrt.csv", index=False)

    colors = {
        "On IMV": mcolors.to_rgba(ORANGE, 0.75),
        "Off IMV": mcolors.to_rgba("#9ecae1", 0.85),
        "Discharged alive": mcolors.to_rgba(GREEN, 0.60),
        "Dead": mcolors.to_rgba("#9e9e9e", 0.75),
    }
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.stackplot(sdf["hour"] / 24.0, sdf["prop_imv"], sdf["prop_off_imv"],
                 sdf["prop_discharged"], sdf["prop_dead"],
                 labels=list(colors), colors=list(colors.values()))
    ax.set_xlim(0, MAXH / 24.0); ax.set_ylim(0, 100)
    ax.set_xlabel("Days from CRRT Initiation"); ax.set_ylabel("Proportion of Patients (%)")
    ax.set_title(f"Invasive Mechanical Ventilation State Over 30 Days Post-CRRT: {SITE_NAME}")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=9)
    fig.tight_layout()
    fig.savefig(GRAPHS / f"{SITE_NAME}_imv_state_over_crrt.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Time-to-event medians (companion to Part B duration & mortality)
    ttd = ot.loc[ot["hours_to_death"].notna() & (ot["hours_to_death"] > 0), "hours_to_death"]
    ttc = ot.loc[ot["hours_to_discharge"].notna() & (ot["hours_to_discharge"] > 0), "hours_to_discharge"]

    def _tte(name, s, unit="hours"):
        return {"metric": name, "n": int(s.count()),
                "median": round(s.median(), 1) if len(s) else np.nan,
                "q25": round(s.quantile(0.25), 1) if len(s) else np.nan,
                "q75": round(s.quantile(0.75), 1) if len(s) else np.nan,
                "unit": unit, "site": SITE_NAME}
    pd.DataFrame([
        _tte("Time to death from CRRT", ttd),
        _tte("Time to discharge alive from CRRT", ttc),
        _tte("Time to death (days)", ttd / 24, "days"),
        _tte("Time to discharge alive (days)", ttc / 24, "days"),
    ]).to_csv(GRAPHS / f"{SITE_NAME}_crrt_course_time_to_event.csv", index=False)


def build_crrt_state(w: pd.DataFrame) -> None:
    """(D5) CRRT-state stacked area over 30 days: On CRRT / Off CRRT (alive) /
    Discharged alive / Dead. "On CRRT" spans t=0 to each encounter's last CRRT
    record; "Off CRRT" = alive and still hospitalized but no longer on CRRT
    (either iHD or no renal replacement — indistinguishable in this CLIF version)."""
    import matplotlib.colors as mcolors
    MAXH = MAX_HOURS_TRAJ
    crrt_init = pd.read_parquet(INTER / "crrt_initiation.parquet")
    oc = pd.read_parquet(INTER / "outcomes_df.parquet")
    cols = ["encounter_block", "died", "death_dttm", "last_vital_dttm"]
    if "discharge_dttm" in oc.columns:
        cols.append("discharge_dttm")
    ot = oc[cols].merge(crrt_init[["encounter_block", "crrt_initiation_time"]],
                        on="encounter_block", how="inner")
    ot["hours_to_death"] = np.where(
        ot["died"] == 1,
        (ot["death_dttm"].fillna(ot["last_vital_dttm"]) - ot["crrt_initiation_time"]).dt.total_seconds() / 3600,
        np.nan)
    disc_col = "discharge_dttm" if "discharge_dttm" in ot.columns else "last_vital_dttm"
    ot["hours_to_discharge"] = np.where(
        ot["died"] == 0,
        (ot[disc_col] - ot["crrt_initiation_time"]).dt.total_seconds() / 3600,
        np.nan)

    # Last CRRT-record hour per encounter (defines the on-CRRT span from t=0).
    if "crrt_mode_category" in w.columns:
        cr = w[(w["crrt_mode_category"].notna())
               & (w["hours_from_crrt"] >= 0) & (w["hours_from_crrt"] <= MAXH)]
        last_crrt_h = cr.groupby("encounter_block")["hours_from_crrt"].max()
    else:
        last_crrt_h = pd.Series(dtype=float)

    all_patients = crrt_init["encounter_block"].unique()
    n_pat = len(all_patients)
    dead_ebs = ot.loc[ot["hours_to_death"].notna(), ["encounter_block", "hours_to_death"]]
    disc_ebs = ot.loc[ot["hours_to_discharge"].notna(), ["encounter_block", "hours_to_discharge"]]

    rows = []
    for h in range(0, MAXH + 1):
        dead_by = set(dead_ebs.loc[dead_ebs["hours_to_death"] <= h, "encounter_block"])
        disc_by = set(disc_ebs.loc[disc_ebs["hours_to_discharge"] <= h, "encounter_block"])
        n_dead = n_disc = n_on = n_off = 0
        for eb in all_patients:
            if eb in dead_by:
                n_dead += 1
            elif eb in disc_by:
                n_disc += 1
            elif eb in last_crrt_h.index and h <= last_crrt_h[eb]:
                n_on += 1
            else:
                n_off += 1
        rows.append({"hour": h, "n_total": n_pat, "n_on_crrt": n_on,
                     "n_off_crrt": n_off, "n_discharged": n_disc, "n_dead": n_dead})
    sdf = pd.DataFrame(rows)
    for k in ["on_crrt", "off_crrt", "discharged", "dead"]:
        sdf[f"prop_{k}"] = sdf[f"n_{k}"] / sdf["n_total"] * 100
    sdf["site"] = SITE_NAME
    sdf.to_csv(GRAPHS / f"{SITE_NAME}_crrt_state_over_crrt.csv", index=False)

    colors = {
        "On CRRT": mcolors.to_rgba(BLUE, 0.80),
        "Off CRRT (alive; iHD or no RRT)": mcolors.to_rgba("#9ecae1", 0.85),
        "Discharged alive": mcolors.to_rgba(GREEN, 0.60),
        "Dead": mcolors.to_rgba("#9e9e9e", 0.75),
    }
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.stackplot(sdf["hour"] / 24.0, sdf["prop_on_crrt"], sdf["prop_off_crrt"],
                 sdf["prop_discharged"], sdf["prop_dead"],
                 labels=list(colors), colors=list(colors.values()))
    ax.set_xlim(0, MAXH / 24.0)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Days from CRRT Initiation")
    ax.set_ylabel("Proportion of Patients (%)")
    ax.set_title(f"CRRT State Over 30 Days Post-CRRT: {SITE_NAME}")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=9)
    fig.tight_layout()
    fig.savefig(GRAPHS / f"{SITE_NAME}_crrt_state_over_crrt.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


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

print("\n[C] Distribution figures ...")
build_figures(inc)
print(f"  wrote figures to {GRAPHS}")

print("\n[D] Longitudinal trajectories (30 d, 3 d inset) ...")
_traj_w = _load_traj_wide()
build_trajectories(_traj_w)
build_patient_course(_traj_w)
build_crrt_state(_traj_w)
print(f"  wrote trajectory figures to {GRAPHS}")

print(f"\n=== 03 CRRT epidemiology complete in {time.time() - t_start:.1f}s ===")
