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
lactate, cci_score (all orthogonal to SOFA components; no double-counting).
A covariate-enrichment sensitivity test (adding bicarbonate + potassium, plus a
restricted-cubic-spline non-linearity check) was run BEFORE the 2026-06-29
cohort harmonization and has NOT been repeated since. Its conclusion —
enrichment left discrimination and the SMR essentially unchanged, so parsimony
wins for portability (plan §4.3.3) — may well still hold, but THE NUMBERS IT
REPORTED (AUC 0.652 vs 0.654; dev-site SMR 1.052 vs 1.053) DESCRIBE A COHORT AND
A MODEL THAT NO LONGER EXIST and must not be quoted. Re-run the comparison
against the current cohort before citing it. Imputation = frozen development
medians (NOT MICE — transported prediction model; plan §4.3.2).
Diagnostics: VIF, Wald SEs / ORs / 95% CIs, in-sample + bootstrap
optimism-corrected AUC, decile calibration. Unpenalized MLE -> in-sample
sum(p)=sum(y), so the development site's own SMR is ~1.0 by construction (and
so its in-sample calibration slope is 1 and intercept 0 by construction, which
is why neither is reported here — calibration is only informative on transfer,
where 03b measures it).

FINGERPRINT SAFETY: 03b hashes ONLY {covariates, coef, intercept,
impute_medians} to fingerprint the deployed model, and refuses to pool sites
whose fingerprints disagree. Diagnostics may therefore be added to the JSON
freely, but the fit itself must not be re-specified casually: any change to the
four hashed keys invalidates every site that has already re-run. This script
prints UNCHANGED/CHANGED against the model it is overwriting so that is never a
surprise.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import polars as pl
import pandas as pd
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score

from pipeline_helpers import load_config, get_output_root

config = load_config()
DEV_SITE = config["site_name"]
COVARS = ["age", "female", "sofa_total", "lactate", "cci_score"]

# The dev cohort is row-level PHI, so 03b writes it to intermediate_phi/
# (03b:22,227-228). This was "smr/" — the pre-rename layout — and was missed
# by the July output rename, so the fitter could not find a cohort that had
# just been built. Honors config['output_dir'] for dev-site isolation.
OUT = get_output_root(config) / "intermediate_phi"
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

# ── Wald standard errors, odds ratios, 95% CIs ──────────────────────────────
# Computed analytically from the observed information at the fitted solution
# rather than by refitting under statsmodels: statsmodels is not in the pinned
# compute stack, and a second optimizer could perturb the coefficients in their
# last decimals — which would change the fingerprint and invalidate every site
# that has already re-run. This reads the sklearn fit, it does not redo it.
#   cov(beta_hat) = (Xd' W Xd)^-1,  W = diag(p(1-p)),  Xd = [1, X]
# For an unpenalized logistic MLE this is the same covariance statsmodels would
# report. With no statsmodels available to cross-check against, the closed form
# is verified below against a FINITE-DIFFERENCE Hessian of the log-likelihood,
# which shares no algebra with it — a genuinely independent check, unlike
# re-deriving the same expression a second way.
Xd = np.column_stack([np.ones(len(X)), X])
W = p_in * (1 - p_in)
fisher = Xd.T @ (W[:, None] * Xd)
cov = np.linalg.inv(fisher)
se_all = np.sqrt(np.diag(cov))
beta_all = np.concatenate([[intercept], clf.coef_[0]])
z_all = beta_all / se_all
p_all = 2 * norm.sf(np.abs(z_all))
_names = ["(Intercept)"] + COVARS
se = {n: float(s) for n, s in zip(_names, se_all)}
# ORs are undefined/uninterpretable for the intercept, so covariates only.
odds_ratio = {c: float(np.exp(b)) for c, b in zip(COVARS, clf.coef_[0])}
or_ci = {c: [float(np.exp(b - 1.96 * s)), float(np.exp(b + 1.96 * s))]
         for c, b, s in zip(COVARS, clf.coef_[0], se_all[1:])}
wald_p = {n: float(pv) for n, pv in zip(_names, p_all)}


def _neg_loglik(b):
    eta = Xd @ b
    # log(1+exp(eta)) computed stably for large |eta|
    return float(np.sum(np.logaddexp(0.0, eta) - y * eta))


# Central-difference Hessian of the negative log-likelihood at beta_hat. At the
# MLE this equals the observed information, so its inverse must equal `cov`.
_h = 1e-5
_H = np.empty((len(beta_all), len(beta_all)))
for _i in range(len(beta_all)):
    for _j in range(len(beta_all)):
        _ei, _ej = np.zeros(len(beta_all)), np.zeros(len(beta_all))
        _ei[_i], _ej[_j] = _h, _h
        _H[_i, _j] = (_neg_loglik(beta_all + _ei + _ej) - _neg_loglik(beta_all + _ei - _ej)
                      - _neg_loglik(beta_all - _ei + _ej) + _neg_loglik(beta_all - _ei - _ej)
                      ) / (4 * _h * _h)
_se_fd = np.sqrt(np.diag(np.linalg.inv(_H)))
_se_dev = float(np.max(np.abs(_se_fd - se_all) / se_all))
if _se_dev > 1e-3:
    raise RuntimeError(
        f"Analytic and finite-difference standard errors disagree by {_se_dev:.2%}. "
        "The closed-form covariance is wrong, or the fit did not converge — do NOT "
        "ship these standard errors.")

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
         "vif": {k: round(v, 2) for k, v in vif.items()},
         # ── Diagnostics for the supplement's model card ─────────────────────
         # NONE of these are read by 03b or enter the fingerprint; they exist so
         # the deployed model can be reported and reproduced without re-fitting.
         "se": {k: round(v, 5) for k, v in se.items()},
         "odds_ratio": {k: round(v, 4) for k, v in odds_ratio.items()},
         "or_ci95": {k: [round(v[0], 4), round(v[1], 4)] for k, v in or_ci.items()},
         "wald_p": wald_p,
         "dev_deaths": int(y.sum()),
         "dev_events_per_variable": round(float(y.sum()) / len(COVARS), 1),
         # Pre-imputation missingness in the DEVELOPMENT cohort, as a fraction.
         # Reported so a reader can see how much of the frozen-median imputation
         # was actually exercised where the coefficients were estimated.
         "dev_missingness": {k: round(v, 4) for k, v in miss.items()},
         "auc_bootstrap_reps": B,
         # Which outcome this model predicts. 03b compares the site's OBSERVED
         # deaths against EXPECTED from this model, which is only meaningful when
         # both count the same thing; it refuses to run on a mismatch and warns
         # when this key is absent. death_30d was re-anchored from ~admission to
         # CRRT initiation on 2026-08-11 (5fda0cc), and a model fitted before
         # that predicts a different outcome.
         "outcome_definition": "death_30d_from_crrt_initiation"}
# Resolved against THIS FILE, not the cwd — as "../config/..." it silently wrote
# the frozen model to a path that only existed when run from code/.
MODEL_PATH = Path(__file__).resolve().parent.parent / "config" / "smr_reference_model.json"


def _fingerprint(m):
    """The 12-char id 03b computes. MUST stay in sync with 03b:274-280 — it
    hashes the four keys that change predictions, and nothing else."""
    return hashlib.sha256(json.dumps(
        {"covariates": m.get("covariates"), "coef": m.get("coef"),
         "intercept": m.get("intercept"), "impute_medians": m.get("impute_medians")},
        sort_keys=True).encode()).hexdigest()[:12]


# A refit that changes the fingerprint invalidates every site that has already
# run 03b against the old one: their SMRs were computed under different
# coefficients and pooling them together is meaningless. That is sometimes the
# intent (a genuine refit) and sometimes an accident (a stray optimizer change
# while adding diagnostics), so say which happened, loudly, every time.
_new_fp = _fingerprint(model)
_old_fp = _fingerprint(json.loads(MODEL_PATH.read_text())) if MODEL_PATH.exists() else None
MODEL_PATH.write_text(json.dumps(model, indent=2))
if _old_fp is None:
    print(f"\n  FINGERPRINT {_new_fp} (no previous model on disk)")
elif _old_fp == _new_fp:
    print(f"\n  FINGERPRINT {_new_fp} UNCHANGED — predictions are identical to the "
          "deployed model; sites that have already re-run stay valid.")
else:
    print(f"\n  *** FINGERPRINT CHANGED {_old_fp} -> {_new_fp} ***\n"
          "  This is a genuine refit. EVERY site must re-run 03b before its SMR can\n"
          "  be pooled; results from the old fingerprint are not comparable.")

# ── Report ──────────────────────────────────────────────────────────────────
print(f"\n-- model card (n={len(dev):,}; deaths={int(y.sum()):,}; "
      f"EPV={y.sum()/len(COVARS):.0f}) --")
print(f"  {'term':13s} {'beta':>9s} {'SE':>8s} {'OR':>7s} {'95% CI':>16s} "
      f"{'p':>9s} {'VIF':>5s}")
print(f"  {'(Intercept)':13s} {intercept:+9.4f} {se['(Intercept)']:8.4f} "
      f"{'—':>7s} {'—':>16s} {wald_p['(Intercept)']:9.2g} {'—':>5s}")
for c in COVARS:
    print(f"  {c:13s} {coef[c]:+9.4f} {se[c]:8.4f} {odds_ratio[c]:7.3f} "
          f"{f'{or_ci[c][0]:.3f}-{or_ci[c][1]:.3f}':>16s} {wald_p[c]:9.2g} {vif[c]:5.2f}")
print(f"\n  standard errors verified against a finite-difference Hessian "
      f"(max deviation {_se_dev:.2e})")
print(f"  AUC in-sample {auc:.3f}; optimism-corrected {auc_corr:.3f} ({B} bootstrap reps)")
cal = (pd.DataFrame({"y": y, "p": p_in})
       .assign(d=lambda t: pd.qcut(t["p"], 10, labels=False, duplicates="drop"))
       .groupby("d").agg(obs=("y", "mean"), exp=("p", "mean"), n=("p", "size")))
print("  calibration (decile obs vs exp):")
for _, r in cal.iterrows():
    print(f"    n={int(r['n']):4d}  obs {r['obs']*100:5.1f}%  exp {r['exp']*100:5.1f}%")
print(f"\nWROTE {MODEL_PATH}  (tracked; ship to consortium sites)")
