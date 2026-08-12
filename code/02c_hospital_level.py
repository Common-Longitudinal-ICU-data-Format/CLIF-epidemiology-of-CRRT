"""
02c_hospital_level.py — hospital-level descriptive aggregates.

Stratifies the CRRT cohort by `hospital_id_at_init` (the facility where CRRT was
started; descriptive-at-initiation only — transfers make it ambiguous over the
course). Produces three per-site CSVs:
  1. {site}_hospital_counts.csv          — encounters per hospital
  2. {site}_table1_by_hospital.csv       — the full Table 1, one column per hospital
  3. {site}_crrt_settings_by_hospital.csv — 3h-median CRRT settings + mode mix per hospital

SHARING: these are shareable as written (changed 2026-08-11). Two protections,
    both required — either alone is insufficient:
      1. SMALL-CELL SUPPRESSION at n < SUPPRESS_MIN_N: hospitals below the
         threshold are dropped entirely, and categorical cells with a numerator
         below it render as "<N" rather than a count.
      2. PSEUDONYMOUS HOSPITAL LABELS: the raw hospital_id never leaves the site.
         Hospitals are relabelled {SITE}_H1..{SITE}_Hk in DESCENDING order of
         encounter count. The crosswalk back to the real ids is written to
         intermediate_phi/ (local, git-ignored) so the site can audit its own
         output; it is never part of the shareable set.
    The consortium builders rely on the {SITE}_H<k> form — report_core's
    collect_hospital_dose parses the trailing ordinal for its figure labels.

Run from code/ :  uv run python 02c_hospital_level.py
"""
import json  # noqa: F401
from pipeline_helpers import load_config, load_intermediate, get_output_root  # noqa: E402
import pandas as pd
import numpy as np

# ---- small-cell suppression -------------------------------------------------
# Consortium threshold, set 2026-08-11: drop hospitals with < 10 encounters and
# render categorical cells with a numerator < 10 as "<10". Setting this to None
# reverts to unsuppressed output, which must then NOT leave the site.
#
# Distinct from the ANALYSIS threshold applied downstream: the pooled
# hospital-dose figure additionally drops hospitals under
# report_core.HOSP_DOSE_MIN_N (20) as too small for a stable median. That is a
# statistical judgement, not a disclosure control; this one is the disclosure
# control. Do not collapse them into a single number.
SUPPRESS_MIN_N = 10

config = load_config()
OUTPUT_ROOT = get_output_root(config)
INTERMEDIATE_DIR = OUTPUT_ROOT / "intermediate_phi"
OUT_DIR = OUTPUT_ROOT / "final_no_phi" / "hospital_level"   # shareable: suppressed + pseudonymous
OUT_DIR.mkdir(parents=True, exist_ok=True)
CROSSWALK_DIR = INTERMEDIATE_DIR / "hospital_level"         # local only: label -> real hospital_id
CROSSWALK_DIR.mkdir(parents=True, exist_ok=True)
SITE_NAME = config["site_name"]
HAS_CRRT_SETTINGS = config.get("has_crrt_settings", True)
HCOL = "hospital_id_at_init"

print("=== 02c: hospital-level aggregates ===")

t1 = load_intermediate(INTERMEDIATE_DIR / "tableone_full_df.parquet")
idx = load_intermediate(INTERMEDIATE_DIR / "index_crrt_df.parquet")

# hospital label: raw id as string; missing -> "unknown"
for _df in (t1, idx):
    if HCOL not in _df.columns:
        _df[HCOL] = pd.NA
    _df[HCOL] = _df[HCOL].astype("string").fillna("unknown")

# ---- pseudonymize hospital ids ------------------------------------------
# Applied HERE, at the single point where HCOL is normalized, so every
# downstream group/header/count uses the pseudonym and no later code has to
# remember to. Ordered by DESCENDING encounter count, so H1 is the largest.
# "unknown" is a category, not a hospital, and passes through unchanged.
_raw_counts = t1[HCOL].value_counts(dropna=False)
_real = [h for h in _raw_counts.index if h != "unknown"]
_pseudo = {h: f"{SITE_NAME}_H{i + 1}" for i, h in enumerate(_real)}
_pseudo["unknown"] = "unknown"
pd.DataFrame({"label": [_pseudo[h] for h in _raw_counts.index],
              "hospital_id_at_init": list(_raw_counts.index),
              "n_encounters": _raw_counts.values}).to_csv(
    CROSSWALK_DIR / f"{SITE_NAME}_hospital_id_crosswalk.csv", index=False)
for _df in (t1, idx):
    if HCOL in _df.columns:
        _df[HCOL] = _df[HCOL].map(lambda v: _pseudo.get(v, v)).astype("string")
print(f"  pseudonymized {len(_real)} hospital(s) -> {SITE_NAME}_H1..H{len(_real)}; "
      f"crosswalk (LOCAL ONLY) -> {CROSSWALK_DIR}")

_counts = t1[HCOL].value_counts(dropna=False)
# groups: Overall + each hospital (unknown last), optionally suppressing small ones
_hosps = [h for h in _counts.index if h != "unknown"] + (["unknown"] if "unknown" in _counts.index else [])
if SUPPRESS_MIN_N:
    _dropped = [h for h in _hosps if _counts[h] < SUPPRESS_MIN_N]
    _hosps = [h for h in _hosps if _counts[h] >= SUPPRESS_MIN_N]
    if _dropped:
        print(f"  suppression: dropped {len(_dropped)} hospital(s) with n < {SUPPRESS_MIN_N}: "
              + ", ".join(f"{h} (n={_counts[h]})" for h in _dropped))
GROUPS = ["Overall"] + _hosps
print(f"  {SITE_NAME}: {t1[HCOL].nunique()} distinct hospital_id_at_init | {len(t1):,} encounters")


# =====================================================================
# 0. Site-level hospital SUMMARY  (for pooled reporting)
# ---------------------------------------------------------------------
# The consortium needs "k sites comprising N hospitals". That total CANNOT be
# recovered by counting rows of hospital_counts.csv, which lists only the
# hospitals surviving suppression — counting those rows undercounts every site
# where suppression fired. So the totals are recorded explicitly, here, before
# any filtering. Counts of hospitals are not patient data; this is shareable.
# =====================================================================
_n_unknown = int(_raw_counts.get("unknown", 0))
_n_reported = len([h for h in _hosps if h != "unknown"])
pd.DataFrame([{
    "site": SITE_NAME,
    "n_hospitals_total": len(_real),               # distinct real hospital_id, pre-suppression
    "n_hospitals_reported": _n_reported,           # surviving n >= SUPPRESS_MIN_N
    "n_hospitals_suppressed": len(_real) - _n_reported,
    "n_encounters_total": int(len(t1)),
    "n_encounters_unknown_hospital": _n_unknown,
    "suppress_min_n": SUPPRESS_MIN_N if SUPPRESS_MIN_N else 0,
}]).to_csv(OUT_DIR / f"{SITE_NAME}_hospital_summary.csv", index=False)
print(f"  hospital summary: {len(_real)} total, {_n_reported} reported, "
      f"{len(_real) - _n_reported} suppressed, {_n_unknown} encounter(s) with unknown hospital_id")


# =====================================================================
# 0b. Hospital TYPE (academic / community)
# ---------------------------------------------------------------------
# hospital_type is a required ADT column, captured at CRRT initiation by 00
# alongside hospital_id. Emitted as its own small file so the consortium can
# describe practice by hospital type without any site shipping a hospital name.
# Restricted to the surviving groups, like hospital_counts. A site whose ADT
# leaves the column empty simply produces no file.
# =====================================================================
TCOL = "hospital_type_at_init"
_type_src = next((f for f in (idx, t1) if TCOL in f.columns), None)
if _type_src is not None and _type_src[TCOL].notna().any():
    _ht = (_type_src[[HCOL, TCOL]].dropna(subset=[TCOL])
           .astype({HCOL: "string", TCOL: "string"}))
    _ht[TCOL] = _ht[TCOL].str.strip().str.lower()
    # One type per hospital: take the most frequent, and say so if a hospital
    # ever disagrees with itself (a data-quality signal, not something to hide).
    _agg = (_ht.groupby(HCOL)[TCOL]
            .agg(hospital_type=lambda s: s.mode().iloc[0] if not s.mode().empty else pd.NA,
                 n_distinct_types="nunique")
            .reset_index().rename(columns={HCOL: "hospital"}))
    _mixed = _agg[_agg["n_distinct_types"] > 1]
    if len(_mixed):
        print(f"  WARNING: {len(_mixed)} hospital(s) report >1 hospital_type; "
              f"taking the most frequent: {', '.join(_mixed['hospital'])}")
    _agg = _agg[_agg["hospital"].isin(_hosps)]
    _n = _type_src.groupby(HCOL).size().rename("n_encounters")
    _agg = _agg.merge(_n, left_on="hospital", right_index=True, how="left")
    _agg[["hospital", "hospital_type", "n_encounters"]].to_csv(
        OUT_DIR / f"{SITE_NAME}_hospital_type.csv", index=False)
    print(f"  hospital_type: {_agg['hospital_type'].value_counts().to_dict()}")
else:
    print("  hospital_type: not available (ADT column absent or empty) — file not written")


# =====================================================================
# 1. Hospital counts
# =====================================================================
# Restricted to the SURVIVING groups. Built from the raw value_counts this
# would have re-published every suppressed hospital with an exact count,
# defeating the suppression applied to the tables beside it.
_cnt = (t1[HCOL].value_counts(dropna=False).rename_axis("hospital_id").reset_index(name="n_encounters"))
_cnt = _cnt[_cnt["hospital_id"].isin(_hosps)].copy()
# Denominator stays the FULL site cohort, so the percentages remain
# interpretable and visibly fail to sum to 100 when anything was suppressed.
_cnt["pct_of_site"] = (_cnt["n_encounters"] / len(t1) * 100).round(1)
_cnt.to_csv(OUT_DIR / f"{SITE_NAME}_hospital_counts.csv", index=False)


# =====================================================================
# 2. Full Table 1 by hospital  (mirrors 02's variable spec, no p-value)
# =====================================================================
gframe = {"Overall": t1, **{h: t1[t1[HCOL] == h] for h in _hosps}}
gn = {g: len(gframe[g]) for g in GROUPS}
COL = "Characteristic"
HDR = {g: (f"Overall (N={gn[g]:,})" if g == "Overall" else f"{g} (N={gn[g]:,})") for g in GROUPS}
_columns = [COL] + [HDR[g] for g in GROUPS]


def _suppress_cell(n):
    return SUPPRESS_MIN_N is not None and 0 < n < SUPPRESS_MIN_N


def _iqr(s, dec):
    s = pd.to_numeric(s, errors="coerce").dropna()
    if s.empty:
        return "NA"
    return f"{s.median():.{dec}f} ({s.quantile(0.25):.{dec}f}, {s.quantile(0.75):.{dec}f})"


def _npct(n, d):
    if _suppress_cell(n):
        return f"<{SUPPRESS_MIN_N}"
    return "NA" if d == 0 else f"{n} ({n / d * 100:.0f}%)"


_rows = []


def _row_cont(label, col, dec):
    if col not in t1.columns:
        return
    r = {COL: f"__{label}__"}
    for g in GROUPS:
        r[HDR[g]] = _iqr(gframe[g][col], dec)
    _rows.append(r)


def _row_binary(label, boolcol):
    if boolcol not in t1.columns:
        return
    r = {COL: f"__{label}__"}
    for g in GROUPS:
        f = gframe[g]
        r[HDR[g]] = _npct(int(pd.Series(f[boolcol]).fillna(False).astype(bool).sum()), len(f))
    _rows.append(r)


def _row_multi(label, col, level_order):
    if col not in t1.columns:
        return
    _rows.append({COL: f"__{label}__", **{HDR[g]: "" for g in GROUPS}})
    for lv, disp in level_order:
        r = {COL: disp}
        for g in GROUPS:
            f = gframe[g]
            r[HDR[g]] = _npct(int((f[col] == lv).sum()), len(f))
        _rows.append(r)


# Hospital type first: it describes the COLUMN (the hospital), so a reader knows
# which stratum is academic and which is community before reading any patient
# characteristic, without joining {SITE}_hospital_type.csv.
#
# In the per-hospital columns this is 100%/0% by construction — a hospital has
# one type — and that doubles as a self-check. The OVERALL column is the
# informative one: it gives the share of the site's CRRT delivered at academic
# vs community hospitals, which nothing else in this file reports.
#
# Kept as a ROW, not folded into the column header: report_core's
# collect_hospital_dose parses the trailing "_H<k>" off the header to label the
# pooled caterpillar, and appending ", academic" to the label would break it.
if TCOL in t1.columns and t1[TCOL].notna().any():
    _types = sorted(t1[TCOL].dropna().astype("string").str.strip().str.lower().unique())
    t1[TCOL] = t1[TCOL].astype("string").str.strip().str.lower()
    for _g in GROUPS:                       # gframe holds slices taken before this normalization
        gframe[_g] = gframe[_g].assign(**{TCOL: gframe[_g][TCOL].astype("string").str.strip().str.lower()})
    _row_multi("Hospital Type", TCOL, [(t, t.capitalize()) for t in _types])

# Demographics
_row_cont("Age at Admission (years)", "age_at_admission", 0)
_row_binary("Female (%)", "_female")
_row_multi("Race", "_race_grp", [("Black", "Black"), ("White", "White"), ("Other", "Other")])
_row_cont("Charlson Comorbidity Index (total)", "cci_score", 0)
# Severity and labs at CRRT initiation
_row_cont("SOFA Score", "sofa_total", 0)
_row_cont("SOFA Score, Non-Renal", "sofa_nonrenal", 0)
_row_cont("Creatinine (mg/dL)", "creatinine_t1", 2)
_row_cont("BUN (mg/dL)", "bun_t1", 0)
_row_cont("Lactate (mmol/L)", "lactate_t1", 1)
_row_cont("Bicarbonate (mEq/L)", "bicarbonate_t1", 1)
_row_cont("Potassium (mEq/L)", "potassium_t1", 2)
_row_cont("Phosphate (mg/dL)", "phosphate_t1", 1)
_row_cont("Arterial pH", "ph_arterial_t1", 2)
_row_cont("NE Equivalent (mcg/kg/min)", "nee_baseline", 2)
_row_binary("On IMV (%)", "_imv")
# CRRT practice descriptors
_row_cont("Initial CRRT Dose (mL/kg/hr)", "crrt_dose_ml_kg_hr", 1)
if HAS_CRRT_SETTINGS and "crrt_mode_category" in t1.columns:
    _modes = list(t1["crrt_mode_category"].dropna().astype("string").str.lower().value_counts().index)
    if _modes:
        _row_multi("CRRT Modality", "crrt_mode_category", [(m, str(m).upper()) for m in _modes])
    # Net UF Intensity omitted (UF sign convention unreliable across sites) — mirrors 02.
# Outcome
# Labels match 02's site-wide Table 1: both are IN-HOSPITAL (death comes from
# the discharge disposition, so post-discharge death is unobservable) and the
# 90-day row is in_hosp_death, which un-counts deaths >90d after CRRT start.
_row_binary("30-Day In-Hospital Mortality (%)", "_death30")
if "in_hosp_death" in t1.columns:
    _row_binary("90-Day In-Hospital Mortality (%)", "in_hosp_death")

pd.DataFrame(_rows, columns=_columns).to_csv(OUT_DIR / f"{SITE_NAME}_table1_by_hospital.csv", index=False)


# =====================================================================
# 3. CRRT 3h-median settings by hospital
# =====================================================================
_SET = [("blood_flow_rate", "Blood Flow Rate (mL/min)", 0),
        ("dialysate_flow_rate", "Dialysate Flow Rate (mL/hr)", 0),
        ("pre_filter_replacement_fluid_rate", "Pre-filter Replacement (mL/hr)", 0),
        ("post_filter_replacement_fluid_rate", "Post-filter Replacement (mL/hr)", 0),
        ("ultrafiltration_out", "Ultrafiltration Out (mL/hr)", 0),
        ("crrt_dose_ml_kg_hr", "Delivered Dose (mL/kg/hr)", 1)]
if HCOL not in idx.columns:
    idx[HCOL] = "unknown"
_ig = {"Overall": idx, **{h: idx[idx[HCOL] == h] for h in _hosps}}
_srows = []
for col, label, dec in _SET:
    if col not in idx.columns:
        continue
    _srows.append({COL: f"__{label}__ (median [IQR])",
                   **{HDR[g]: _iqr(_ig[g][col], dec) for g in GROUPS}})
# mode mix
if HAS_CRRT_SETTINGS and "crrt_mode_category" in idx.columns:
    _m = idx["crrt_mode_category"].dropna().astype("string").str.lower()
    _srows.append({COL: "__CRRT Modality (n %)__", **{HDR[g]: "" for g in GROUPS}})
    for mode in list(_m.value_counts().index):
        r = {COL: str(mode).upper()}
        for g in GROUPS:
            f = _ig[g]
            n = int((f["crrt_mode_category"].astype("string").str.lower() == mode).sum())
            r[HDR[g]] = _npct(n, len(f))
        _srows.append(r)
pd.DataFrame(_srows, columns=_columns).to_csv(OUT_DIR / f"{SITE_NAME}_crrt_settings_by_hospital.csv", index=False)

print(f"  wrote 4 files -> {OUT_DIR}")
print(f"  (pseudonymous labels + n<{SUPPRESS_MIN_N} suppression -> shareable)"
      if SUPPRESS_MIN_N else
      "  (SUPPRESS_MIN_N is None -> UNSUPPRESSED; do NOT share this output)")
print("Done!")
