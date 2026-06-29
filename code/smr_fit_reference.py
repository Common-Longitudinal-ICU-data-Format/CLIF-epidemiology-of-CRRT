"""Fit the SMR reference (case-mix) model ON THE DEVELOPMENT SITE (MIMIC).

DEV/COORDINATING SITE ONLY — consortium application sites do NOT run this; they
receive the frozen coefficient JSON and apply it via 03b_crrt_epi_smr.py. This
is the one place MIMIC-CLIF is required (Option E, docs/smr_addition_plan.md §4.4).

    CLIF_CONFIG=config/config_mimic.json uv run python 03b_crrt_epi_smr.py   # build cohort
    CLIF_CONFIG=config/config_mimic.json uv run python smr_fit_reference.py  # fit + freeze

Reads  output[_root]/smr/{DEV_SITE}_smr_cohort.parquet
Writes config/smr_reference_model.json   (TRACKED; shipped to all sites)

Cohort (harmonized 2026-06-29): 03b builds {DEV_SITE}_smr_cohort.parquet from the
SAME analytic CRRT cohort as Table 1 (00->02; adult, continuous CRRT, ESRD-
excluded, with required baseline data) with NO minimum-CRRT-duration filter, so
the reference model is fit on the same cohort definition the consortium applies.
At MIMIC this requires the full pipeline run into its isolated output_mimic/ tree
(CLIF_CONFIG=config/config_mimic.json on 00 -> 01 -> 02 -> 03b).

Model: pre-specified PARSIMONIOUS LINEAR logistic — age, female, sofa_total,
lactate, cci_score (all orthogonal to SOFA components; no double-counting). We
tested enriching with bicarbonate + potassium and a non-linearity (restricted
cubic spline) check: only bicarbonate was non-linear (p=0.015) yet enrichment
left discrimination unchanged (AUC 0.652 vs 0.654) and the SMR unchanged
(UChicago 1.052 vs 1.053), so for portability/parsimony the deployed model is
linear in the five covariates (rationale in plan §4.3.3). Imputation = frozen
development medians (NOT MICE — transported prediction model; plan §4.3.2).
Diagnostics: VIF, in-sample + bootstrap optimism-corrected AUC, decile
calibration. Unpenalized MLE -> in-sample sum(p)=sum(y), so the development
site's own SMR is ~1.0 by construction.
"""
import json
from pathlib import Path

import numpy as np
import polars as pl
import pandas as pd
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score

from pipeline_helpers import load_config, get_output_root

config = load_config()
DEV_SITE = config["site_name"]
COVARS = ["age", "female", "sofa_total", "lactate", "cci_score"]

OUT = get_output_root(config) / "smr"  # honors config['output_dir'] (dev-site isolation)
cohort_path = OUT / f"{DEV_SITE}_smr_cohort.parquet"
if not cohort_path.exists():
    raise FileNotFoundError(
        f"{cohort_path} not found. Build the development cohort first:\n"
        f"  CLIF_CONFIG=config/config_mimic.json uv run python 03b_crrt_epi_smr.py")

print(f"=== SMR reference-model fit | development site = {DEV_SITE} ===")
dev = pl.read_parquet(cohort_path).to_pandas()
dev = dev.dropna(subset=["died_30d"])
y = dev["died_30d"].to_numpy(int)

# Frozen single imputation: development medians for every covariate (plan §4.3.2)
impute_medians = {cv: float(np.nanmedian(dev[cv])) for cv in COVARS}
miss = {cv: float(dev[cv].isna().mean()) for cv in COVARS}
for cv in COVARS:
    dev[cv] = dev[cv].fillna(impute_medians[cv])
print(f"  n={len(dev):,}; missingness pre-impute: "
      + ", ".join(f"{cv} {miss[cv]*100:.0f}%" for cv in COVARS if miss[cv] > 0.001))

# ── Fit parsimonious linear logistic ────────────────────────────────────────
X = dev[COVARS].to_numpy(float)
clf = LogisticRegression(penalty=None, max_iter=5000).fit(X, y)
coef = dict(zip(COVARS, clf.coef_[0].astype(float).tolist()))
intercept = float(clf.intercept_[0])
p_in = clf.predict_proba(X)[:, 1]
auc = float(roc_auc_score(y, p_in))

# ── Collinearity (VIF) ──────────────────────────────────────────────────────
vif = {}
for i, c in enumerate(COVARS):
    others = [j for j in range(len(COVARS)) if j != i]
    r2 = LinearRegression().fit(X[:, others], X[:, i]).score(X[:, others], X[:, i])
    vif[c] = float(1 / (1 - r2)) if r2 < 1 else float("inf")

# ── Bootstrap optimism-corrected AUC ────────────────────────────────────────
rng = np.random.default_rng(20260628)
B, n, opt = 200, len(dev), []
for _ in range(B):
    idx = rng.integers(0, n, n)
    clfb = LogisticRegression(penalty=None, max_iter=5000).fit(X[idx], y[idx])
    opt.append(roc_auc_score(y[idx], clfb.predict_proba(X[idx])[:, 1])
               - roc_auc_score(y, clfb.predict_proba(X)[:, 1]))
auc_corr = float(auc - np.mean(opt))

# ── Write frozen model ──────────────────────────────────────────────────────
model = {"dev_site": DEV_SITE, "covariates": COVARS,
         "intercept": intercept, "coef": coef, "impute_medians": impute_medians,
         "dev_n": int(len(dev)), "dev_auc": round(auc, 3),
         "dev_auc_optimism_corrected": round(auc_corr, 3),
         "vif": {k: round(v, 2) for k, v in vif.items()}}
MODEL_PATH = Path("../config/smr_reference_model.json")
MODEL_PATH.write_text(json.dumps(model, indent=2))

# ── Report ──────────────────────────────────────────────────────────────────
print(f"\n-- model (intercept {intercept:+.3f}) --")
for c in COVARS:
    print(f"  {c:12s} {coef[c]:+.4f}")
print(f"\n  AUC in-sample {auc:.3f}; optimism-corrected {auc_corr:.3f}")
print("  VIF: " + ", ".join(f"{c} {vif[c]:.2f}" for c in COVARS))
cal = (pd.DataFrame({"y": y, "p": p_in})
       .assign(d=lambda t: pd.qcut(t["p"], 10, labels=False, duplicates="drop"))
       .groupby("d").agg(obs=("y", "mean"), exp=("p", "mean"), n=("p", "size")))
print("  calibration (decile obs vs exp):")
for _, r in cal.iterrows():
    print(f"    n={int(r['n']):4d}  obs {r['obs']*100:5.1f}%  exp {r['exp']*100:5.1f}%")
print(f"\nWROTE {MODEL_PATH}  (tracked; ship to consortium sites)")
