"""Build a site's SMR cohort (covariates + outcome) on one consistent definition.

Site-agnostic: the data path/site come from the selected config, so the SAME
covariate definitions are used at the MIMIC development site and at every
consortium application site (Option E in docs/smr_addition_plan.md).

    # MIMIC (reference / development cohort)
    CLIF_CONFIG=config/config_mimic.json uv run python smr_build_mimic_dev.py
    # UChicago (or any consortium site; default config)
    uv run python smr_build_mimic_dev.py

Produces ONE artifact, output/smr/{SITE}_smr_cohort.parquet, with the locked SMR
covariates + outcome:
    hospitalization_id, age, female, sofa_total, lactate, cci_score, died_30d

Self-contained: writes only under output/smr/, so it cannot clobber the main
single-site pipeline outputs. Cohort is hospitalization-level (no encounter-block
stitching) — a deliberate approximation for a case-mix standardization model
(see docs/smr_addition_plan.md §4.9). Reuses the pipeline's SOFA calculator and
clifpy CCI so covariate definitions match the rest of the pipeline.
"""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
from datetime import timedelta

import numpy as np
import polars as pl
import pyarrow.parquet as pq
import clifpy

from pipeline_helpers import load_config
from sofa_calculator import compute_sofa_polars

config = load_config()  # expects CLIF_CONFIG=config/config_mimic.json
SITE = config["site_name"]
TABLES_PATH = config["tables_path"]
FILE_TYPE = config["file_type"]
TIMEZONE = config["timezone"]
print(f"=== SMR dev-cohort builder | site={SITE} | tz={TIMEZONE} ===")
print(f"    tables_path: {TABLES_PATH}")

OUT = Path("../output/smr")
OUT.mkdir(parents=True, exist_ok=True)

def P(name): return f"{TABLES_PATH}/clif_{name}.{FILE_TYPE}"
def pct(n, d): return f"{100*n/d:.1f}%" if d else "n/a"

# ── Load tables (selected columns) ──────────────────────────────────────────
crrt = pl.read_parquet(P("crrt_therapy"),
    columns=["hospitalization_id", "recorded_dttm", "crrt_mode_category"])
hosp = pl.read_parquet(P("hospitalization"),
    columns=["hospitalization_id", "patient_id", "age_at_admission",
             "discharge_category", "discharge_dttm"])
pat = pl.read_parquet(P("patient"),
    columns=["patient_id", "sex_category", "death_dttm"])

# ── Cohort funnel (hospitalization-level; mirrors 00 rules) ─────────────────
print("\n-- cohort funnel --")
n0 = crrt["hospitalization_id"].n_unique()
print(f"  0. any CRRT record:                       {n0:,}")

adult_ids = hosp.filter(pl.col("age_at_admission") >= 18)["hospitalization_id"]
crrt = crrt.filter(pl.col("hospitalization_id").is_in(adult_ids.to_list()))
print(f"  1. adult (>=18):                          {crrt['hospitalization_id'].n_unique():,}")

valid_modes = ["cvvhdf", "cvvhd", "cvvh"]
has_valid = (crrt.with_columns(pl.col("crrt_mode_category").str.to_lowercase().alias("m"))
             .filter(pl.col("m").is_in(valid_modes))["hospitalization_id"].unique())
crrt = crrt.filter(pl.col("hospitalization_id").is_in(has_valid.to_list()))
print(f"  2. continuous mode (cvvhdf/cvvhd/cvvh):   {crrt['hospitalization_id'].n_unique():,}")

span = (crrt.group_by("hospitalization_id")
        .agg((pl.col("recorded_dttm").max() - pl.col("recorded_dttm").min()).alias("span")))
ge24 = span.filter(pl.col("span") >= timedelta(hours=24))["hospitalization_id"]
crrt = crrt.filter(pl.col("hospitalization_id").is_in(ge24.to_list()))
print(f"  3. CRRT span >= 24h:                      {crrt['hospitalization_id'].n_unique():,}")

# ESRD exclusion (codes + poa in {1, NA}); normalize codes
esrd_codes = ['z992','z9115','i120','n186','i132','z91158','i1311','z4931','z4901',
              'i272','5856','40391','40311','v4511','v4512']
diag_raw = pl.read_parquet(P("hospital_diagnosis"),
    columns=["hospitalization_id", "diagnosis_code", "diagnosis_code_format", "poa_present"])
diag_norm = diag_raw.with_columns(
    pl.col("diagnosis_code").str.to_lowercase().str.replace_all(r"\.", "").alias("code"))
esrd_ids = set(diag_norm.filter(
    pl.col("code").is_in(esrd_codes) &
    ((pl.col("poa_present") == 1) | (pl.col("poa_present").is_null()))
)["hospitalization_id"].unique().to_list())
cohort_ids = [h for h in crrt["hospitalization_id"].unique().to_list() if h not in esrd_ids]
N = len(cohort_ids)
print(f"  4. after ESRD exclusion == DEV COHORT:    {N:,}")

# CRRT initiation = first record per hospitalization
init = (crrt.filter(pl.col("hospitalization_id").is_in(cohort_ids))
        .group_by("hospitalization_id")
        .agg(pl.col("recorded_dttm").min().alias("crrt_init")))

# ── Demographics + 30-day mortality ─────────────────────────────────────────
# Death definition mirrors 00_cohort.py: died = discharge_category in
# {expired, hospice}; death timing = death_dttm if present else discharge_dttm
# (00 falls back to last_vital; discharge_dttm is the standalone-builder proxy).
# died_30d = died AND death timing within 30 days of CRRT initiation.
base = (hosp.filter(pl.col("hospitalization_id").is_in(cohort_ids))
        .join(pat, on="patient_id", how="left")
        .join(init, on="hospitalization_id", how="left")
        .with_columns([
            pl.col("age_at_admission").alias("age"),
            (pl.col("sex_category").str.to_lowercase().is_in(["female", "f"]))
                .cast(pl.Int8).alias("female"),
            pl.col("discharge_category").str.to_lowercase().is_in(["expired", "hospice"])
                .alias("died_flag"),
            pl.coalesce([pl.col("death_dttm"), pl.col("discharge_dttm")]).alias("death_time"),
        ])
        .with_columns(
            (pl.col("died_flag")
             & pl.col("death_time").is_not_null()
             & ((pl.col("death_time") - pl.col("crrt_init")) >= timedelta(0))
             & ((pl.col("death_time") - pl.col("crrt_init")) <= timedelta(days=30)))
                .cast(pl.Int8).alias("died_30d")
        ))
_died_any = int(base["died_flag"].sum())
print(f"\n-- mortality (discharge_category in expired/hospice) --")
print(f"  died in hospital (any time): {_died_any:,} ({pct(_died_any, N)})")
print(f"  died within 30d of CRRT init: {int(base['died_30d'].sum()):,} "
      f"({pct(int(base['died_30d'].sum()), N)})")

# ── SOFA total at baseline (-12h..+3h), via the pipeline calculator ─────────
print("\n-- computing SOFA (compute_sofa_polars, -12h..+3h) --")
sofa_cohort = init.with_columns([
    (pl.col("crrt_init") - pl.duration(hours=12)).alias("start_dttm"),
    (pl.col("crrt_init") + pl.duration(hours=3)).alias("end_dttm"),
]).select(["hospitalization_id", "start_dttm", "end_dttm"])
sofa = compute_sofa_polars(
    data_directory=TABLES_PATH, cohort_df=sofa_cohort, filetype=FILE_TYPE,
    id_name="hospitalization_id", extremal_type="worst",
    fill_na_scores_with_zero=False, remove_outliers=True, timezone=TIMEZONE,
)
sofa = sofa.select(["hospitalization_id", "sofa_total"])
print(f"  SOFA computed for {sofa.height:,} hosps; sofa_total non-null "
      f"{pct(sofa['sofa_total'].drop_nulls().len(), N)}")

# ── Baseline lactate: nearest measured in [-12h,+3h] ────────────────────────
print("\n-- baseline lactate (nearest in -12h..+3h) --")
# Prefer lab_value_numeric (clean double); fall back to parsing lab_value string
labs_cols = pq.read_schema(P("labs")).names
val_col = "lab_value_numeric" if "lab_value_numeric" in labs_cols else "lab_value"
lac = (pl.scan_parquet(P("labs"))
       .select(["hospitalization_id", "lab_category", "lab_result_dttm", val_col])
       .filter((pl.col("hospitalization_id").is_in(cohort_ids))
               & (pl.col("lab_category") == "lactate")).collect()
       .with_columns(pl.col(val_col).cast(pl.Float64, strict=False).alias("lactate"))
       .filter(pl.col("lactate").is_not_null())
       .join(init, on="hospitalization_id", how="inner")
       .with_columns((pl.col("lab_result_dttm") - pl.col("crrt_init")).alias("dt"))
       .filter((pl.col("dt") >= timedelta(hours=-12)) & (pl.col("dt") <= timedelta(hours=3)))
       .with_columns(pl.col("dt").abs().alias("absdt"))
       .sort(["hospitalization_id", "absdt"])
       .group_by("hospitalization_id").first()
       .select(["hospitalization_id", "lactate"]))
print(f"  (lactate source column: {val_col})")
print(f"  lactate available: {pct(lac['hospitalization_id'].n_unique(), N)}")
if lac.height:
    print(f"  lactate median (UNIT CHECK: ~1-5=mmol/L good, ~10-40=mg/dL mismatch): "
          f"{lac['lactate'].median():.2f}")

# ── CCI: ICD-10 via clifpy + ICD-9 via Quan 2005 mapping ────────────────────
# clifpy maps ICD10CM only; MIMIC has a large ICD-9 share (older encounters),
# so ICD-9-only hospitalizations would be dropped (~33% missing, non-random by
# era). We add a Quan(2005) ICD-9-CM Charlson mapping, OR condition flags from
# both formats at the condition level, then apply clifpy's hierarchy + weights
# (read from cci.yaml) so ICD-9 scores land on the exact clifpy scale.
print("\n-- CCI (ICD-10 via clifpy + ICD-9 via Quan 2005 mapping) --")
# clifpy's Quan weights + hierarchy (from cci.yaml; matched exactly)
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
# Quan (2005) enhanced ICD-9-CM code prefixes (dot-stripped, lowercased)
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

diag = diag_raw.filter(pl.col("hospitalization_id").is_in(cohort_ids)).to_pandas()
diag["diagnosis_code"] = diag["diagnosis_code"].astype(str).str.replace(".", "", regex=False).str.lower()
fmt_counts = diag["diagnosis_code_format"].value_counts(dropna=False)
print(f"  diagnosis_code_format mix:\n{fmt_counts.to_string()}")
print("  (POA not used for CCI: MIMIC poa_present is sparse)")

def _flags_from_clifpy_icd10(d10):
    res = clifpy.calculate_cci(
        d10[["hospitalization_id", "diagnosis_code", "diagnosis_code_format"]],
        hierarchy=False)  # raw flags; hierarchy applied after OR-ing formats
    if isinstance(res, pl.DataFrame):
        res = res.to_pandas()
    keep = ["hospitalization_id"] + [c for c in CONDS if c in res.columns]
    return pl.from_pandas(res[keep]).with_columns(
        [pl.col(c).cast(pl.Int8) for c in keep[1:]])

def _flags_from_icd9(d9):
    rows = pl.from_pandas(d9[["hospitalization_id", "diagnosis_code"]])
    exprs = [pl.col("diagnosis_code").str.contains("^(?:" + "|".join(ICD9[c]) + ")")
             .alias(c) for c in CONDS]
    return (rows.with_columns(exprs)
            .group_by("hospitalization_id")
            .agg([pl.col(c).any().cast(pl.Int8) for c in CONDS]))

d10 = diag[diag["diagnosis_code_format"].astype(str).str.upper() == "ICD10CM"]
d9 = diag[diag["diagnosis_code_format"].astype(str).str.upper() == "ICD9CM"]
f10 = _flags_from_clifpy_icd10(d10) if len(d10) else None
f9 = _flags_from_icd9(d9) if len(d9) else None

# OR condition flags across formats (full outer join, max per condition)
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
for c in CONDS:
    cdf[c] = cdf[c].fillna(0).astype(int)
for hi, lo in HIER:                       # assign0 hierarchy
    cdf.loc[cdf[hi] == 1, lo] = 0
cdf["cci_score"] = sum(cdf[c] * w for c, w in CCI_WEIGHTS.items())

# Scale-consistency check: my score must reproduce clifpy's on ICD-10-only hosps
chk = clifpy.calculate_cci(
    d10[["hospitalization_id", "diagnosis_code", "diagnosis_code_format"]], hierarchy=True)
if isinstance(chk, pl.DataFrame):
    chk = chk.to_pandas()
icd9_ids = set(f9["hospitalization_id"].to_list()) if f9 is not None else set()
only10 = cdf[~cdf["hospitalization_id"].isin(icd9_ids)]
m = only10.merge(chk[["hospitalization_id", "cci_score"]], on="hospitalization_id",
                 suffixes=("_mine", "_clifpy"))
agree = (m["cci_score_mine"] == m["cci_score_clifpy"]).mean() if len(m) else float("nan")
print(f"  scale check (my score == clifpy on ICD-10-only, n={len(m):,}): {agree*100:.1f}% agree")

cci = pl.from_pandas(cdf[["hospitalization_id", "cci_score"]]).with_columns(
    pl.col("hospitalization_id").cast(pl.Utf8))
print(f"  cci_score available: {pct(cci['hospitalization_id'].n_unique(), N)}; "
      f"median {cci['cci_score'].median():.1f}")

# ── Assemble + write ────────────────────────────────────────────────────────
dev = (base.select(["hospitalization_id", "age", "female", "died_30d"])
       .join(sofa, on="hospitalization_id", how="left")
       .join(lac, on="hospitalization_id", how="left")
       .join(cci, on="hospitalization_id", how="left"))

print("\n-- dev cohort assembled --")
print(f"  n = {dev.height:,}")
print(f"  died_30d rate: {pct(dev['died_30d'].sum(), dev.height)}")
for c in ["age", "female", "sofa_total", "lactate", "cci_score"]:
    print(f"  {c:12s} non-null {pct(dev[c].drop_nulls().len(), dev.height)}")

out_path = OUT / f"{SITE}_smr_cohort.parquet"
dev.write_parquet(out_path)
print(f"\nWROTE {out_path.resolve()}")
