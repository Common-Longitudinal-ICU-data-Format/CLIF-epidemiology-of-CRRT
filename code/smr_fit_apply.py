"""Fit the SMR reference (case-mix) model on the development cohort (Option E:
MIMIC), freeze the coefficients, and apply them to every site's SMR cohort to
get O/E, the standardized mortality ratio with exact-Poisson (Byar) 95% CIs,
AUC, and a calibration check. Renders an (anonymized) forest plot.

Run from code/ (after smr_build_cohort.py has produced the site cohorts):
    uv run python smr_fit_apply.py

Reads  output/smr/{SITE}_smr_cohort.parquet
Writes output/smr/smr_reference_model.json  (frozen coefficients; aggregate, safe)
       output/smr/smr_results.csv           (per-site O/E/SMR/CI/AUC; aggregate)
       output/smr/smr_forest.png            (anonymized forest)

Reference model is fit ONLY on the development site (MIMIC), so every consortium
site -- including its own home site -- is evaluated against an external standard
(no site is privileged; see docs/smr_addition_plan.md s4.4). Unpenalized logistic
MLE => in-sample sum(p)=sum(y), so the development site's SMR is ~1.0 by
construction (a calibration sanity check, not a real comparison point).
"""
import json
import numpy as np
import polars as pl
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

OUT = Path("../output/smr")
DEV_SITE = "MIMIC"
COVARS = ["age", "female", "sofa_total", "lactate", "cci_score"]


def byar_ci(O, E):
    """Exact-Poisson (Byar) 95% CI for O/E."""
    O = int(O)
    if O == 0:
        return (0.0, 3.689 / E)
    lo = (O / E) * (1 - 1/(9*O) - 1.96/(3*np.sqrt(O)))**3
    hi = ((O + 1) / E) * (1 - 1/(9*(O+1)) + 1.96/(3*np.sqrt(O+1)))**3
    return lo, hi


def load_site(site):
    return pl.read_parquet(OUT / f"{site}_smr_cohort.parquet").to_pandas()


# ── Fit the reference model on the development cohort ────────────────────────
dev = load_site(DEV_SITE)
lac_median = float(np.nanmedian(dev["lactate"]))


def prep(d):
    d = d.copy()
    d["lactate"] = d["lactate"].fillna(lac_median)
    return d.dropna(subset=COVARS + ["died_30d"])


devp = prep(dev)
X = devp[COVARS].to_numpy(float)
y = devp["died_30d"].to_numpy(int)
clf = LogisticRegression(penalty=None, max_iter=5000)
clf.fit(X, y)
intercept = float(clf.intercept_[0])
coef = dict(zip(COVARS, clf.coef_[0].astype(float).tolist()))


def predict(d):
    lp = intercept + sum(coef[c] * d[c] for c in COVARS)
    return 1.0 / (1.0 + np.exp(-lp))


dev_auc = roc_auc_score(y, predict(devp))
model = {"dev_site": DEV_SITE, "covariates": COVARS, "intercept": intercept,
         "coef": coef, "lactate_impute_median": lac_median,
         "dev_n": int(len(devp)), "dev_auc": round(float(dev_auc), 3)}
(OUT / "smr_reference_model.json").write_text(json.dumps(model, indent=2))
print(f"Reference model (development = {DEV_SITE}, n={len(devp):,}):")
print(f"  intercept {intercept:+.3f}; " + "; ".join(f"{c} {coef[c]:+.3f}" for c in COVARS))
print(f"  in-sample AUC = {dev_auc:.3f}; lactate impute median = {lac_median:.2f}")

# Calibration on the development cohort (observed vs expected by risk decile)
dp = devp.assign(p=predict(devp))
dp["decile"] = pd.qcut(dp["p"], 10, labels=False, duplicates="drop")
cal = dp.groupby("decile").agg(obs=("died_30d", "mean"), exp=("p", "mean"), n=("p", "size"))
print("\nDevelopment-cohort calibration (decile: observed vs expected):")
for _, r in cal.iterrows():
    print(f"  n={int(r['n']):4d}  observed {r['obs']*100:5.1f}%   expected {r['exp']*100:5.1f}%")

# ── Apply frozen coefficients to every site ─────────────────────────────────
rows = []
for f in sorted(OUT.glob("*_smr_cohort.parquet")):
    site = f.stem.replace("_smr_cohort", "")
    d = prep(load_site(site))
    p = predict(d)
    O = int(d["died_30d"].sum())
    E = float(p.sum())
    lo, hi = byar_ci(O, E)
    auc = roc_auc_score(d["died_30d"], p) if d["died_30d"].nunique() > 1 else np.nan
    rows.append({"site": site, "n": int(len(d)), "observed": O, "expected": round(E, 1),
                 "smr": round(O / E, 3), "lo": round(lo, 3), "hi": round(hi, 3),
                 "auc": round(float(auc), 3), "crude_mort_pct": round(100 * O / len(d), 1),
                 "is_dev": site == DEV_SITE})
res = pd.DataFrame(rows)
res.to_csv(OUT / "smr_results.csv", index=False)
print("\nPer-site SMR (O/E vs the MIMIC-developed reference):")
print(res.drop(columns="is_dev").to_string(index=False))

# ── Forest plot (anonymized consortium sites; dev site marked) ──────────────
res_plot = res.sort_values("is_dev")  # consortium sites first, dev last
labels, anon = [], 0
for _, r in res_plot.iterrows():
    if r["is_dev"]:
        labels.append(f"{r['site']} (reference)")
    else:
        anon += 1
        labels.append(f"Site {anon}")
yy = np.arange(len(res_plot))[::-1]
fig, ax = plt.subplots(figsize=(7, 0.7 * len(res_plot) + 1.5))
xerr = np.array([(res_plot["smr"] - res_plot["lo"]).values,
                 (res_plot["hi"] - res_plot["smr"]).values])
colors = ["#999999" if dv else "#1e417c" for dv in res_plot["is_dev"]]
for i, (yi, c) in enumerate(zip(yy, colors)):
    ax.errorbar(res_plot["smr"].iloc[i], yi, xerr=xerr[:, i:i+1], fmt="o",
                color=c, ecolor=c, capsize=4, markersize=8)
ax.axvline(1.0, ls="--", color="black", lw=1)
ax.set_yticks(yy); ax.set_yticklabels(labels)
ax.set_xlabel("Standardized Mortality Ratio (O/E), 95% CI\nreference model developed in MIMIC")
ax.set_title("Risk-Standardized 30-Day Mortality (SMR)")
for i, yi in enumerate(yy):
    r = res_plot.iloc[i]
    ax.text(ax.get_xlim()[1], yi, f"  {r['smr']:.2f} ({r['lo']:.2f}-{r['hi']:.2f})",
            va="center", fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "smr_forest.png", dpi=200, bbox_inches="tight")
print(f"\nWROTE {(OUT/'smr_reference_model.json')}, smr_results.csv, smr_forest.png")
