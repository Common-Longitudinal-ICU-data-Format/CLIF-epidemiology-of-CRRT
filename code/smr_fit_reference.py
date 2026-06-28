"""Fit the SMR reference (case-mix) model ON THE DEVELOPMENT SITE (MIMIC).

DEV/COORDINATING SITE ONLY — consortium application sites do NOT run this; they
receive the frozen coefficient JSON and apply it via 03b_crrt_epi_smr.py. This
is the one place MIMIC-CLIF is required (Option E, docs/smr_addition_plan.md §4.4).

    CLIF_CONFIG=config/config_mimic.json uv run python 03b_crrt_epi_smr.py   # build cohort
    CLIF_CONFIG=config/config_mimic.json uv run python smr_fit_reference.py  # fit + freeze

Reads  output/smr/{DEV_SITE}_smr_cohort.parquet
Writes config/smr_reference_model.json   (TRACKED; shipped to all sites)

Covariates: age, female, sofa_total, lactate, bicarbonate, potassium, cci_score
(all orthogonal to SOFA components — no double-counting). Continuous covariates
get a non-linearity LR test (4-knot restricted cubic spline vs linear); a spline
is adopted only where significant (knots stored, so every site applies the same
basis). Imputation = frozen development medians (NOT MICE — this is a transported
prediction model; see plan §4.3.2). Diagnostics: VIF, in-sample + bootstrap
optimism-corrected AUC, decile calibration. Unpenalized MLE -> in-sample sum(p)=
sum(y) so the development site's own SMR is ~1.0 by construction.
"""
import json
from pathlib import Path

import numpy as np
import polars as pl
import pandas as pd
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score
from scipy.stats import chi2

from pipeline_helpers import load_config
from smr_common import design_matrix

config = load_config()
DEV_SITE = config["site_name"]
COVARS = ["age", "female", "sofa_total", "lactate", "bicarbonate", "potassium", "cci_score"]
CONTINUOUS = ["age", "sofa_total", "lactate", "bicarbonate", "potassium"]  # test these
KNOT_Q = [0.05, 0.35, 0.65, 0.95]   # Harrell 4-knot placement (-> 2 nonlinear terms)
ALPHA = 0.05

OUT = Path("../output/smr")
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


def fit_logit(X):
    return LogisticRegression(penalty=None, max_iter=5000).fit(X, y)


def loglik(clf, X):
    p = np.clip(clf.predict_proba(X)[:, 1], 1e-12, 1 - 1e-12)
    return float(np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)))


# ── Non-linearity LR test per continuous covariate (spline vs linear) ───────
base_spec = {"covariates": COVARS, "splines": {}}
X0, _ = design_matrix(dev, base_spec)
ll0 = loglik(fit_logit(X0), X0)

nl_p, knots_for, splines = {}, {}, {}
print("\n-- non-linearity test (4-knot RCS vs linear, 2 df) --")
for cv in CONTINUOUS:
    kn = list(np.quantile(dev[cv].to_numpy(float), KNOT_Q))
    if len(set(np.round(kn, 6))) < len(kn):       # tied knots -> can't spline
        nl_p[cv] = float("nan")
        print(f"  {cv:12s} knots tied; skipped")
        continue
    spec = {"covariates": COVARS, "splines": {cv: kn}}
    Xs, _ = design_matrix(dev, spec)
    stat = 2 * (loglik(fit_logit(Xs), Xs) - ll0)
    p = float(chi2.sf(stat, len(kn) - 2))
    nl_p[cv] = p
    knots_for[cv] = kn
    flag = "  -> SPLINE" if p < ALPHA else ""
    print(f"  {cv:12s} p={p:.4f}{flag}")
    if p < ALPHA:
        splines[cv] = kn

# ── Final model ─────────────────────────────────────────────────────────────
spec = {"covariates": COVARS, "splines": splines}
X, names = design_matrix(dev, spec)
clf = fit_logit(X)
coef = dict(zip(names, clf.coef_[0].astype(float).tolist()))
intercept = float(clf.intercept_[0])
p_in = clf.predict_proba(X)[:, 1]
auc = float(roc_auc_score(y, p_in))

# ── Collinearity (VIF on the base covariates) ───────────────────────────────
M = dev[COVARS].to_numpy(float)
vif = {}
for i, c in enumerate(COVARS):
    others = [j for j in range(len(COVARS)) if j != i]
    r2 = LinearRegression().fit(M[:, others], M[:, i]).score(M[:, others], M[:, i])
    vif[c] = float(1 / (1 - r2)) if r2 < 1 else float("inf")

# ── Bootstrap optimism-corrected AUC (fixed final spec) ─────────────────────
rng = np.random.default_rng(20260628)
B, n, opt = 200, len(dev), []
for _ in range(B):
    idx = rng.integers(0, n, n)
    Xb, yb = X[idx], y[idx]
    clfb = LogisticRegression(penalty=None, max_iter=5000).fit(Xb, yb)
    opt.append(roc_auc_score(yb, clfb.predict_proba(Xb)[:, 1])
               - roc_auc_score(y, clfb.predict_proba(X)[:, 1]))
auc_corr = float(auc - np.mean(opt))

# ── Write frozen model ──────────────────────────────────────────────────────
model = {"dev_site": DEV_SITE, "covariates": COVARS, "splines": splines,
         "intercept": intercept, "coef": coef, "impute_medians": impute_medians,
         "dev_n": int(len(dev)), "dev_auc": round(auc, 3),
         "dev_auc_optimism_corrected": round(auc_corr, 3),
         "nonlinearity_p": {k: (round(v, 4) if v == v else None) for k, v in nl_p.items()},
         "vif": {k: round(v, 2) for k, v in vif.items()}}
MODEL_PATH = Path("../config/smr_reference_model.json")
MODEL_PATH.write_text(json.dumps(model, indent=2))

# ── Report ──────────────────────────────────────────────────────────────────
print(f"\n-- model (intercept {intercept:+.3f}) --")
for nm in names:
    print(f"  {nm:18s} {coef[nm]:+.4f}")
print(f"\n  splined: {list(splines) or 'none'}")
print(f"  AUC in-sample {auc:.3f}; optimism-corrected {auc_corr:.3f}")
print("  VIF: " + ", ".join(f"{c} {vif[c]:.2f}" for c in COVARS))
cal = (pd.DataFrame({"y": y, "p": p_in})
       .assign(d=lambda t: pd.qcut(t["p"], 10, labels=False, duplicates="drop"))
       .groupby("d").agg(obs=("y", "mean"), exp=("p", "mean"), n=("p", "size")))
print("  calibration (decile obs vs exp):")
for _, r in cal.iterrows():
    print(f"    n={int(r['n']):4d}  obs {r['obs']*100:5.1f}%  exp {r['exp']*100:5.1f}%")
print(f"\nWROTE {MODEL_PATH}  (tracked; ship to consortium sites)")
