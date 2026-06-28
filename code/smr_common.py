"""Shared SMR helpers, imported by both the dev fitter (smr_fit_reference.py)
and the per-site apply (03b_crrt_epi_smr.py) so the design matrix, spline basis,
imputation, and prediction are byte-for-byte identical at fit and deployment.

Frozen-model JSON schema (config/smr_reference_model.json):
    dev_site, dev_n, dev_auc, dev_auc_optimism_corrected
    covariates : ordered list of base covariates
    splines    : {covar: [knots]} for covariates adopted as restricted cubic splines
    intercept  : float
    coef       : {term_name: value}; term = covar (linear) or "covar__s{j}" (spline)
    impute_medians : {covar: development median}  (frozen single imputation)
    nonlinearity_p : {covar: LR-test p}   vif : {covar: value}   (diagnostics)
"""
import numpy as np


def byar_ci(O, E):
    """Exact-Poisson (Byar) 95% CI for the SMR = O/E."""
    O = int(O)
    if O == 0:
        return (0.0, 3.689 / E)
    lo = (O / E) * (1 - 1/(9*O) - 1.96/(3*np.sqrt(O)))**3
    hi = ((O + 1) / E) * (1 - 1/(9*(O+1)) + 1.96/(3*np.sqrt(O+1)))**3
    return lo, hi


def rcs_terms(x, knots):
    """Harrell restricted cubic spline nonlinear basis: (k-2) columns for k knots,
    normalized by (t_last - t_first)^2. The raw linear term is added separately by
    design_matrix(); these are only the nonlinear terms."""
    x = np.asarray(x, float)
    t = np.asarray(knots, float)
    k = len(t)
    denom = (t[-1] - t[0]) ** 2
    def cub(u):
        return np.where(u > 0, u, 0.0) ** 3
    cols = []
    for j in range(k - 2):
        cols.append((cub(x - t[j])
                     - cub(x - t[k-2]) * (t[k-1] - t[j]) / (t[k-1] - t[k-2])
                     + cub(x - t[k-1]) * (t[k-2] - t[j]) / (t[k-1] - t[k-2])) / denom)
    return np.column_stack(cols) if cols else np.empty((len(x), 0))


def design_matrix(df, spec):
    """(X, names) for a model spec. spec['covariates'] is the ordered base list;
    spec['splines'] = {covar: knots}. A spline covar contributes its raw linear
    column plus (k-2) nonlinear columns named covar__s1.. ."""
    splines = spec.get("splines", {}) or {}
    cols, names = [], []
    for cv in spec["covariates"]:
        cols.append(np.asarray(df[cv], float))
        names.append(cv)
        if cv in splines:
            st = rcs_terms(df[cv], splines[cv])
            for j in range(st.shape[1]):
                cols.append(st[:, j])
                names.append(f"{cv}__s{j+1}")
    return np.column_stack(cols), names


def impute(df, model):
    """Frozen single imputation: fill each covariate's NA with its development
    median stored in the model (outcome-free, identical at every site)."""
    df = df.copy()
    for cv, med in model.get("impute_medians", {}).items():
        if cv in df.columns:
            df[cv] = df[cv].fillna(med)
    return df


def predict(df, model):
    """Predicted P(death) under the frozen reference model."""
    X, names = design_matrix(df, model)
    beta = np.array([model["coef"][n] for n in names], float)
    lp = model["intercept"] + X @ beta
    return 1.0 / (1.0 + np.exp(-lp))
