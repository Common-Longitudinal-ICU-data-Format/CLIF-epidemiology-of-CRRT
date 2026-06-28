"""Fit the SMR reference (case-mix) model ON THE DEVELOPMENT SITE (MIMIC).

DEV/COORDINATING SITE ONLY — consortium application sites do NOT run this; they
receive the frozen coefficient JSON and apply it via 03b_crrt_epi_smr.py. This
is the one place MIMIC-CLIF is required (Option E, docs/smr_addition_plan.md §4.4).

    # build the MIMIC cohort first (writes output/smr/MIMIC_smr_cohort.parquet):
    CLIF_CONFIG=config/config_mimic.json uv run python 03b_crrt_epi_smr.py
    # then fit + freeze the reference model:
    CLIF_CONFIG=config/config_mimic.json uv run python smr_fit_reference.py

Reads  output/smr/{DEV_SITE}_smr_cohort.parquet   (DEV_SITE = config site_name)
Writes config/smr_reference_model.json            (TRACKED; shipped to all sites)

Unpenalized logistic MLE -> in-sample sum(p)=sum(y), so the development site's
own SMR is ~1.0 by construction (a calibration sanity check, not a comparison).
"""
import json
from pathlib import Path

import numpy as np
import polars as pl
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from pipeline_helpers import load_config

config = load_config()
DEV_SITE = config["site_name"]
COVARS = ["age", "female", "sofa_total", "lactate", "cci_score"]

OUT = Path("../output/smr")
cohort_path = OUT / f"{DEV_SITE}_smr_cohort.parquet"
if not cohort_path.exists():
    raise FileNotFoundError(
        f"{cohort_path} not found. Build the development cohort first:\n"
        f"  CLIF_CONFIG=config/config_mimic.json uv run python 03b_crrt_epi_smr.py")

print(f"=== SMR reference-model fit | development site = {DEV_SITE} ===")
dev = pl.read_parquet(cohort_path).to_pandas()
lac_median = float(np.nanmedian(dev["lactate"]))
dev["lactate"] = dev["lactate"].fillna(lac_median)
n_pre = len(dev)
dev = dev.dropna(subset=COVARS + ["died_30d"])
print(f"  n={len(dev):,} (dropped {n_pre - len(dev)} for missing covariates); "
      f"lactate impute median={lac_median:.2f}")

X = dev[COVARS].to_numpy(float)
y = dev["died_30d"].to_numpy(int)
clf = LogisticRegression(penalty=None, max_iter=5000)
clf.fit(X, y)
intercept = float(clf.intercept_[0])
coef = dict(zip(COVARS, clf.coef_[0].astype(float).tolist()))

lp = intercept + sum(coef[c] * dev[c] for c in COVARS)
p = 1.0 / (1.0 + np.exp(-lp))
auc = roc_auc_score(y, p)

model = {"dev_site": DEV_SITE, "covariates": COVARS, "intercept": intercept,
         "coef": coef, "lactate_impute_median": lac_median,
         "dev_n": int(len(dev)), "dev_auc": round(float(auc), 3)}
MODEL_PATH = Path("../config/smr_reference_model.json")
MODEL_PATH.write_text(json.dumps(model, indent=2))

print(f"  intercept {intercept:+.3f}; " + "; ".join(f"{c} {coef[c]:+.3f}" for c in COVARS))
print(f"  in-sample AUC {auc:.3f}")
# Calibration by risk decile (observed vs expected)
cal = (pd.DataFrame({"y": y, "p": p})
       .assign(decile=lambda d: pd.qcut(d["p"], 10, labels=False, duplicates="drop"))
       .groupby("decile").agg(obs=("y", "mean"), exp=("p", "mean"), n=("p", "size")))
print("  calibration (decile: observed vs expected):")
for _, r in cal.iterrows():
    print(f"    n={int(r['n']):4d}  obs {r['obs']*100:5.1f}%  exp {r['exp']*100:5.1f}%")
print(f"\nWROTE {MODEL_PATH}  (tracked; ship to consortium sites)")
