"""
02c_hospital_level.py — hospital-level descriptive aggregates.

Stratifies the CRRT cohort by `hospital_id_at_init` (the facility where CRRT was
started; descriptive-at-initiation only — transfers make it ambiguous over the
course). Produces four per-site CSVs:
  1. {site}_hospital_counts.csv          — encounters per hospital
  2. {site}_table1_by_hospital.csv       — the full Table 1, one column per hospital
  3. {site}_crrt_settings_by_hospital.csv — 3h-median CRRT settings + mode mix per hospital
  4. {site}_hospital_type.csv            — anon label -> academic/community (for caterpillar)

SHARING: outputs are de-identified for the consortium — hospital labels are
    anonymized to {site}_H1..Hn (ordered by descending encounter count; no
    crosswalk leaves the site), and small cells are suppressed two ways:
      * whole hospitals with < SUPPRESS_MIN_N encounters are dropped entirely, and
      * any categorical cell with numerator < SUPPRESS_MIN_N is blanked to "<N".
    Written to final_no_phi/hospital_level/. Set SUPPRESS_MIN_N = None to disable
    suppression for LOCAL inspection only (never share an unsuppressed run).

Run from code/ :  uv run python 02c_hospital_level.py
"""
import json  # noqa: F401
from pipeline_helpers import load_config, load_intermediate, get_output_root  # noqa: E402
import pandas as pd
import numpy as np

# ---- small-cell suppression (consortium default; None = local inspection only) ----
SUPPRESS_MIN_N = 10     # drop hospitals with < 10 encounters and blank categorical
                        # cells with numerator < 10. Set None to disable (do NOT share).

config = load_config()
OUTPUT_ROOT = get_output_root(config)
INTERMEDIATE_DIR = OUTPUT_ROOT / "intermediate_phi"
OUT_DIR = OUTPUT_ROOT / "final_no_phi" / "hospital_level"   # de-identified + suppressed -> shareable
OUT_DIR.mkdir(parents=True, exist_ok=True)
SITE_NAME = config["site_name"]
HAS_CRRT_SETTINGS = config.get("has_crrt_settings", True)
HCOL = "hospital_id_at_init"
TCOL = "hospital_type_at_init"

print("=== 02c: hospital-level aggregates ===")

t1 = load_intermediate(INTERMEDIATE_DIR / "tableone_full_df.parquet")
idx = load_intermediate(INTERMEDIATE_DIR / "index_crrt_df.parquet")

# hospital label: raw id as string; missing -> "unknown"
for _df in (t1, idx):
    if HCOL not in _df.columns:
        _df[HCOL] = pd.NA
    _df[HCOL] = _df[HCOL].astype("string").fillna("unknown")

_counts = t1[HCOL].value_counts(dropna=False)
# groups: Overall + each hospital (unknown last), optionally suppressing small ones
_hosps = [h for h in _counts.index if h != "unknown"] + (["unknown"] if "unknown" in _counts.index else [])
if SUPPRESS_MIN_N:
    _dropped = [h for h in _hosps if _counts[h] < SUPPRESS_MIN_N]
    _hosps = [h for h in _hosps if _counts[h] >= SUPPRESS_MIN_N]
    if _dropped:
        print(f"  suppression: dropped {len(_dropped)} hospital(s) with n < {SUPPRESS_MIN_N}")
GROUPS = ["Overall"] + _hosps
# anonymize surviving hospitals -> {site}_H1..Hn (descending encounter count, which is
# _hosps' order); "unknown" and "Overall" pass through. The raw->anon map never leaves here.
LABEL = {"Overall": "Overall", "unknown": "unknown"}
_h = 0
for _hid in _hosps:
    if _hid == "unknown":
        continue
    _h += 1
    LABEL[_hid] = f"{SITE_NAME}_H{_h}"
print(f"  {SITE_NAME}: {t1[HCOL].nunique()} distinct hospital_id_at_init | {len(t1):,} encounters")
if SUPPRESS_MIN_N:
    print(f"  suppression ON (min n = {SUPPRESS_MIN_N}); {len([h for h in _hosps if h != 'unknown'])} hospital(s) retained + anonymized")


# =====================================================================
# 1. Hospital counts
# =====================================================================
# only retained (unsuppressed) hospitals, anonymized; pcts are of the full site total
_cnt = pd.DataFrame({"hospital_id": [LABEL[h] for h in _hosps],
                     "n_encounters": [int(_counts[h]) for h in _hosps]})
_cnt["pct_of_site"] = (_cnt["n_encounters"] / len(t1) * 100).round(1)
_cnt.to_csv(OUT_DIR / f"{SITE_NAME}_hospital_counts.csv", index=False)


# =====================================================================
# 2. Full Table 1 by hospital  (mirrors 02's variable spec, no p-value)
# =====================================================================
gframe = {"Overall": t1, **{h: t1[t1[HCOL] == h] for h in _hosps}}
gn = {g: len(gframe[g]) for g in GROUPS}
COL = "Characteristic"
HDR = {g: f"{LABEL[g]} (N={gn[g]:,})" for g in GROUPS}
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


def _row_dose_bins(label, col):
    # N(%) of dosed patients per KDIGO-style band; denominator = patients with a
    # non-missing dose (so bands sum to 100% of the dosed), edges: <20, 20-30, >30.
    if col not in t1.columns:
        return
    bands = [("< 20 mL/kg/hr", lambda s: s < 20),
             ("20–30 mL/kg/hr", lambda s: (s >= 20) & (s <= 30)),
             ("> 30 mL/kg/hr", lambda s: s > 30)]
    _rows.append({COL: f"__{label}__", **{HDR[g]: "" for g in GROUPS}})
    for disp, pred in bands:
        r = {COL: disp}
        for g in GROUPS:
            s = pd.to_numeric(gframe[g][col], errors="coerce")
            r[HDR[g]] = _npct(int(pred(s).sum()), int(s.notna().sum()))
        _rows.append(r)


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
_row_dose_bins("Initial CRRT Dose, banded (% of dosed)", "crrt_dose_ml_kg_hr")
if HAS_CRRT_SETTINGS and "crrt_mode_category" in t1.columns:
    _modes = list(t1["crrt_mode_category"].dropna().astype("string").str.lower().value_counts().index)
    if _modes:
        _row_multi("CRRT Modality", "crrt_mode_category", [(m, str(m).upper()) for m in _modes])
    # Net UF Intensity omitted (UF sign convention unreliable across sites) — mirrors 02.
# Outcome
_row_binary("30-Day Mortality (%)", "_death30")

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


# =====================================================================
# 4. Hospital type (academic / community) — anon label -> type, for the caterpillar
# =====================================================================
def _first_type(s):
    v = s.dropna()
    return str(v.iloc[0]).lower() if len(v) else "unknown"


if TCOL in idx.columns:
    _typ = idx.groupby(HCOL)[TCOL].agg(_first_type)
else:
    _typ = pd.Series(dtype="object")   # site did not populate ADT hospital_type
_trows = [{"hospital": LABEL[h],
           "hospital_type": _typ.get(h, "unknown"),
           "n_encounters": int(_counts[h])}
          for h in _hosps if h != "unknown"]
pd.DataFrame(_trows, columns=["hospital", "hospital_type", "n_encounters"]).to_csv(
    OUT_DIR / f"{SITE_NAME}_hospital_type.csv", index=False)

print(f"  wrote 4 files -> {OUT_DIR}")
print(f"  (anonymized labels + small-cell suppression @ n<{SUPPRESS_MIN_N} -> shareable)"
      if SUPPRESS_MIN_N else "  (SUPPRESSION OFF -> local inspection only, DO NOT share)")
print("Done!")
