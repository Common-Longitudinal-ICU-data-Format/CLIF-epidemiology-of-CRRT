"""03b — Risk-standardized mortality (SMR) for CRRT, per site.

Part of the descriptive-epi family (runs after 02). Site-agnostic: the data
path/site come from the selected config, so identical covariate definitions are
used at the MIMIC development site and every consortium application site
(Option E in docs/smr_addition_plan.md).

Every site runs this:
    uv run python 03b_crrt_epi_smr.py                       # UChicago / any site
    CLIF_CONFIG=config/config_mimic.json uv run python 03b_crrt_epi_smr.py  # MIMIC

COHORT (harmonized 2026-06-29): the SMR uses the SAME analytic CRRT cohort as
Table 1 — `output[_root]/intermediate/tableone_analysis_df.parquet`, built by
00->01->02 (adult >=18, continuous CRRT mode, ESRD-excluded, with required
baseline data). There is **no minimum-CRRT-duration filter**: a mortality-
standardization model must not exclude patients who died early in CRRT (that
would bias the observed count). Cohort key = encounter_block, so the SMR N
equals the Table-1 N exactly. (Previously 03b re-derived a hospitalization-level
cohort with a >=24h span filter; that drift from plan §4.2 is removed.)

It (1) assembles the per-site SMR cohort (covariates + outcome) and writes
output[_root]/smr/{SITE}_smr_cohort.parquet (gitignored, row-level), then (2) if
the frozen reference model config/smr_reference_model.json is present, applies it
to compute observed/expected, the SMR with exact-Poisson (Byar) 95% CI, and AUC,
and writes the per-site exchange CSV output[_root]/final/crrt_epi/{SITE}_smr.csv
(which 07/08 pool into the anonymized cross-site forest + dashboard).

Covariates (5, parsimonious linear): age, female, SOFA total, baseline lactate,
CCI. age / female / sofa_total / lactate (nearest-measured {lab}_t1) / 30-day
mortality are read straight from the pipeline (tz-correct, identical to Table 1).
CCI is the only covariate computed here (clifpy ICD-10 + a Quan-2005 ICD-9 map),
keyed to the cohort's encounter_blocks via the hospitalization map; it has no
datetime/tz dependence.

The reference model is fit ONCE on MIMIC by the coordinating site via
smr_fit_reference.py; consortium sites only need the shipped coefficient JSON,
NOT MIMIC data.
"""
import warnings
warnings.filterwarnings("ignore")
import json
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pyarrow.parquet as pq
import clifpy

from pipeline_helpers import load_config, get_output_root

config = load_config()
SITE = config["site_name"]
TABLES_PATH = config["tables_path"]
FILE_TYPE = config["file_type"]
print(f"=== 03b CRRT SMR | site={SITE} ===")
print(f"    tables_path: {TABLES_PATH}")

OUTPUT_ROOT = get_output_root(config)            # honors config['output_dir']
INTER = OUTPUT_ROOT / "intermediate"
OUT = OUTPUT_ROOT / "smr"
OUT.mkdir(parents=True, exist_ok=True)
EPI_OUT = OUTPUT_ROOT / "final" / "crrt_epi"
EPI_OUT.mkdir(parents=True, exist_ok=True)
MODEL_PATH = Path("../config/smr_reference_model.json")

def P(name): return f"{TABLES_PATH}/clif_{name}.{FILE_TYPE}"
def pct(n, d): return f"{100*n/d:.1f}%" if d else "n/a"

def byar_ci(O, E):
    """Exact-Poisson (Byar) 95% CI for the SMR = O/E."""
    O = int(O)
    if O == 0:
        return (0.0, 3.689 / E)
    lo = (O / E) * (1 - 1/(9*O) - 1.96/(3*np.sqrt(O)))**3
    hi = ((O + 1) / E) * (1 - 1/(9*(O+1)) + 1.96/(3*np.sqrt(O+1)))**3
    return lo, hi

# ── Cohort + covariates from the pipeline (the Table-1 analytic cohort) ──────
TBL_PATH = INTER / "tableone_analysis_df.parquet"
if not TBL_PATH.exists():
    raise FileNotFoundError(
        f"{TBL_PATH} not found. Run the pipeline first (00 -> 01 -> 02):\n"
        f"  uv run python 00_cohort.py && uv run python 01_create_wide_df.py && "
        f"uv run python 02_construct_crrt_tableone.py")
NEED = ["encounter_block", "age_at_admission", "sex_category",
        "sofa_total", "lactate_t1", "death_30d"]
_avail = set(pq.read_schema(TBL_PATH).names)
_missing = [c for c in NEED if c not in _avail]
if _missing:
    raise KeyError(
        f"tableone_analysis_df.parquet is missing {_missing}. Re-run "
        f"02_construct_crrt_tableone.py (it persists nearest-measured lactate_t1).")
tbl = pd.read_parquet(TBL_PATH, columns=NEED)
N = len(tbl)

base = pd.DataFrame({
    "encounter_block": tbl["encounter_block"],
    "age": pd.to_numeric(tbl["age_at_admission"], errors="coerce"),
    "female": tbl["sex_category"].astype("string").str.lower().isin(["female", "f"]).astype("int8"),
    "sofa_total": pd.to_numeric(tbl["sofa_total"], errors="coerce"),
    "lactate": pd.to_numeric(tbl["lactate_t1"], errors="coerce"),
    "died_30d": (pd.to_numeric(tbl["death_30d"], errors="coerce") == 1).astype("int8"),
})
print(f"\n-- cohort (== Table-1 analytic CRRT cohort; no duration filter) --")
print(f"  N encounter_blocks: {N:,}")
print(f"  died within 30d: {int(base['died_30d'].sum()):,} ({pct(int(base['died_30d'].sum()), N)})")
print(f"  completeness: age {pct(base['age'].notna().sum(), N)}, "
      f"SOFA {pct(base['sofa_total'].notna().sum(), N)}, "
      f"lactate {pct(base['lactate'].notna().sum(), N)}")

# ── CCI: ICD-10 via clifpy + ICD-9 via Quan 2005 mapping (same clifpy scale) ─
# Computed here (the only raw-data covariate; no datetime/tz). Keyed to the
# cohort's encounter_blocks via the hospitalization map (cohort_df), OR-ing the
# condition flags across a block's hospitalizations before scoring.
CCI_WEIGHTS = {
    "myocardial_infarction": 0, "congestive_heart_failure": 2,
    "peripheral_vascular_disease": 0, "cerebrovascular_disease": 0,
    "dementia": 2, "chronic_pulmonary_disease": 1, "connective_tissue_disease": 1,
    "peptic_ulcer_disease": 0, "mild_liver_disease": 2, "diabetes_uncomplicated": 0,
    "diabetes_with_complications": 1, "hemiplegia": 2, "renal_disease": 1,
    "cancer": 2, "moderate_severe_liver_disease": 4, "metastatic_solid_tumor": 6,
    "aids": 4,
}
CONDS = list(CCI_WEIGHTS.keys())
HIER = [("diabetes_with_complications", "diabetes_uncomplicated"),
        ("moderate_severe_liver_disease", "mild_liver_disease"),
        ("metastatic_solid_tumor", "cancer")]
ICD9 = {
    "myocardial_infarction": ["410", "412"],
    "congestive_heart_failure": ["39891","40201","40211","40291","40401","40403",
        "40411","40413","40491","40493","4254","4255","4256","4257","4258","4259","428"],
    "peripheral_vascular_disease": ["0930","4373","440","441","4431","4432","4433",
        "4434","4435","4436","4437","4438","4439","4471","5571","5579","v434"],
    "cerebrovascular_disease": ["36234","430","431","432","433","434","435","436","437","438"],
    "dementia": ["290","2941","3312"],
    "chronic_pulmonary_disease": ["4168","4169","490","491","492","493","494","495",
        "496","500","501","502","503","504","505","5064","5081","5088"],
    "connective_tissue_disease": ["4465","7100","7101","7102","7103","7104","7140",
        "7141","7142","7148","725"],
    "peptic_ulcer_disease": ["531","532","533","534"],
    "mild_liver_disease": ["07022","07023","07032","07033","07044","07054","0706",
        "0709","570","571","5733","5734","5738","5739","v427"],
    "diabetes_uncomplicated": ["2500","2501","2502","2503","2508","2509"],
    "diabetes_with_complications": ["2504","2505","2506","2507"],
    "hemiplegia": ["3341","342","343","3440","3441","3442","3443","3444","3445","3446","3449"],
    "renal_disease": ["40301","40311","40391","40402","40403","40412","40413","40492",
        "40493","582","5830","5831","5832","5834","5835","5836","5837","585","586",
        "5880","v420","v451","v56"],
    "cancer": [str(x) for x in range(140,173)] + [str(x) for x in range(174,196)]
               + [str(x) for x in range(200,209)] + ["2386"],
    "moderate_severe_liver_disease": ["4560","4561","4562","5722","5723","5724",
        "5725","5726","5727","5728"],
    "metastatic_solid_tumor": ["196","197","198","199"],
    "aids": ["042","043","044"],
}

# encounter_block <-> hospitalization_id map (pipeline), restricted to the cohort
eb_map = pd.read_parquet(INTER / "cohort_df.parquet",
                         columns=["hospitalization_id", "encounter_block"])
eb_map = eb_map[eb_map["encounter_block"].isin(set(base["encounter_block"]))].copy()
eb_map["hospitalization_id"] = eb_map["hospitalization_id"].astype(str)
cohort_hosp_ids = eb_map["hospitalization_id"].unique().tolist()

diag_raw = pl.read_parquet(P("hospital_diagnosis"),
    columns=["hospitalization_id", "diagnosis_code", "diagnosis_code_format", "poa_present"])
diag = (diag_raw.with_columns(pl.col("hospitalization_id").cast(pl.Utf8))
        .filter(pl.col("hospitalization_id").is_in(cohort_hosp_ids)).to_pandas())
diag["diagnosis_code"] = diag["diagnosis_code"].astype(str).str.replace(".", "", regex=False).str.lower()

def _flags_from_clifpy_icd10(d10):
    res = clifpy.calculate_cci(
        d10[["hospitalization_id", "diagnosis_code", "diagnosis_code_format"]],
        hierarchy=False)
    if isinstance(res, pl.DataFrame):
        res = res.to_pandas()
    keep = ["hospitalization_id"] + [c for c in CONDS if c in res.columns]
    return pl.from_pandas(res[keep]).with_columns([pl.col(c).cast(pl.Int8) for c in keep[1:]])

def _flags_from_icd9(d9):
    rows = pl.from_pandas(d9[["hospitalization_id", "diagnosis_code"]])
    exprs = [pl.col("diagnosis_code").str.contains("^(?:" + "|".join(ICD9[c]) + ")").alias(c)
             for c in CONDS]
    return (rows.with_columns(exprs).group_by("hospitalization_id")
            .agg([pl.col(c).any().cast(pl.Int8) for c in CONDS]))

d10 = diag[diag["diagnosis_code_format"].astype(str).str.upper() == "ICD10CM"]
d9 = diag[diag["diagnosis_code_format"].astype(str).str.upper() == "ICD9CM"]
f10 = _flags_from_clifpy_icd10(d10) if len(d10) else None
f9 = _flags_from_icd9(d9) if len(d9) else None

if f10 is not None and f9 is not None:
    comb = f10.join(f9.rename({c: f"{c}__9" for c in CONDS}),
                    on="hospitalization_id", how="full", coalesce=True)
    for c in CONDS:
        comb = comb.with_columns(
            pl.max_horizontal(pl.col(c).fill_null(0), pl.col(f"{c}__9").fill_null(0)).alias(c)
        ).drop(f"{c}__9")
else:
    comb = f10 if f10 is not None else f9

cdf = comb.to_pandas()
cdf["hospitalization_id"] = cdf["hospitalization_id"].astype(str)
for c in CONDS:
    cdf[c] = cdf[c].fillna(0).astype(int)
# Map hospitalization-level flags to encounter_block, OR across the block's
# hospitalizations, THEN apply the Charlson hierarchy + weights per block.
cdf = cdf.merge(eb_map, on="hospitalization_id", how="inner")
blk = cdf.groupby("encounter_block")[CONDS].max().reset_index()
for hi, lo in HIER:
    blk.loc[blk[hi] == 1, lo] = 0
blk["cci_score"] = sum(blk[c] * w for c, w in CCI_WEIGHTS.items())
cci = blk[["encounter_block", "cci_score"]]
print(f"-- CCI (ICD-10 clifpy + ICD-9 Quan) computed for {pct(len(cci), N)} of cohort; "
      f"median {cci['cci_score'].median():.1f}")

# ── Assemble + write cohort (encounter_block-keyed) ─────────────────────────
cohort = base.merge(cci, on="encounter_block", how="left")
# Blocks with no diagnosis records -> cci_score NA; left to the frozen-median
# imputation in the apply step (consistent with how other missing covariates
# are handled, and with the reference fitter).
cohort_path = OUT / f"{SITE}_smr_cohort.parquet"
cohort.to_parquet(cohort_path, index=False)
print(f"\nWROTE {cohort_path}  (n={len(cohort):,}, 30d mort {pct(int(cohort['died_30d'].sum()), len(cohort))})")

# ── Apply the frozen reference model (if present) → per-site exchange CSV ────
if not MODEL_PATH.exists():
    print(f"\nNo reference model at {MODEL_PATH}. Cohort built only.")
    print("  Run smr_fit_reference.py on the development site (MIMIC) to create it,")
    print("  then re-run 03b at each site to emit output/final/crrt_epi/{SITE}_smr.csv.")
else:
    from sklearn.metrics import roc_auc_score
    model = json.loads(MODEL_PATH.read_text())
    COVARS = model["covariates"]
    coef, intercept = model["coef"], model["intercept"]
    pdf = cohort.copy()
    # frozen development-median imputation (outcome-free; identical at every site)
    for cv, med in model.get("impute_medians", {}).items():
        if cv in pdf.columns:
            pdf[cv] = pdf[cv].fillna(med)
    n_pre = len(pdf)
    pdf = pdf.dropna(subset=COVARS + ["died_30d"])
    # parsimonious linear model: P(death) = sigmoid(intercept + sum coef*covar)
    lp = intercept + sum(coef[c] * pdf[c] for c in COVARS)
    p = 1.0 / (1.0 + np.exp(-lp))
    y = pdf["died_30d"].to_numpy(int)
    O = int(y.sum()); E = float(p.sum()); smr = O / E
    lo, hi = byar_ci(O, E)
    auc = roc_auc_score(y, p) if pdf["died_30d"].nunique() > 1 else float("nan")
    row = {"site": SITE, "n": int(len(pdf)), "n_dropped_missing": int(n_pre - len(pdf)),
           "observed": O, "expected": round(E, 2), "smr": round(smr, 3),
           "lo": round(lo, 3), "hi": round(hi, 3), "auc": round(float(auc), 3),
           "crude_mort_pct": round(100 * O / len(pdf), 1),
           "reference_dev_site": model.get("dev_site", "")}
    pd.DataFrame([row]).to_csv(EPI_OUT / f"{SITE}_smr.csv", index=False)
    # Transfer calibration: observed vs expected by predicted-risk decile
    cal = (pd.DataFrame({"y": y, "p": p})
           .assign(decile=lambda t: pd.qcut(t["p"], 10, labels=False, duplicates="drop"))
           .groupby("decile").agg(n=("p", "size"), observed=("y", "sum"),
                                  expected=("p", "sum")).reset_index())
    cal.insert(0, "site", SITE)
    cal["expected"] = cal["expected"].round(2)
    cal.to_csv(EPI_OUT / f"{SITE}_smr_calibration.csv", index=False)
    print(f"\n-- applied reference model (dev={model.get('dev_site','?')}, "
          f"parsimonious linear) --")
    print(f"  SMR {smr:.3f} ({lo:.3f}-{hi:.3f}); O={O} E={E:.1f}; AUC {auc:.3f}; "
          f"crude {row['crude_mort_pct']}%")
    print(f"WROTE {EPI_OUT / f'{SITE}_smr.csv'} + _smr_calibration.csv")
