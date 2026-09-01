#!/usr/bin/env python
"""
Diagnostic: why is SOFA missing for so many CRRT patients at this site?

Runs against a site's own CLIF tables and reports three things:

  A. Which SOFA COMPONENT is missing (renal / coag / liver / resp / cns / cv),
     by null rate on the CRRT cohort. The component with a high null rate is the
     one blocking the total SOFA. This is VERSION-INDEPENDENT: each component
     score is null when its input is missing, even on polars versions where the
     summed total silently treats the missing component as 0. (This step
     recomputes SOFA and may take a few minutes.)

  B. Coverage of each SOFA INPUT in the source tables (labs / patient_assessments
     / vitals): is the field present, under what exact category name/casing, and
     for how many hospitalizations does it actually carry a value.

  C. Environment (polars / clifpy version, git commit) so the coordinating center
     can reproduce.

The console output is AGGREGATE ONLY (category names + counts + rates) — no
patient-level data — so it is safe to copy back to the coordinating center.

Usage (from the repo root):
    # point CLIF_CONFIG at your site config, or have config/config.json present
    export CLIF_CONFIG=config.json          # or: config/config_<site>.json
    python code/debug_sofa_coverage.py

On Windows PowerShell:
    $env:CLIF_CONFIG="config.json"; python code/debug_sofa_coverage.py
"""
import os
import sys
import json
import subprocess
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # so `from sofa_calculator import ...` works

# SOFA component column -> (input field it needs, source table key or None)
SOFA_INPUTS = {
    "sofa_renal": ("creatinine",      "labs"),
    "sofa_coag":  ("platelet_count",  "labs"),
    "sofa_liver": ("bilirubin_total", "labs"),
    "sofa_cns":   ("gcs_total",       "patient_assessments"),
    "sofa_cv_97": ("map",             "vitals"),
    "sofa_resp":  ("po2_arterial/spo2/fio2_set", None),  # multi-table; not a single-field check
}
# source table key -> (category column, candidate numeric-value columns in priority order)
# Value-column names vary by site/ETL; patient_assessments in particular stores the
# number in `numerical_value` (the SOFA loader coalesces numerical_value+categorical_value).
_TABLE_META = {
    "labs":                ("lab_category",        ["lab_value_numeric", "lab_value"]),
    "patient_assessments": ("assessment_category", ["numerical_value", "categorical_value", "assessment_value"]),
    "vitals":              ("vital_category",      ["vital_value", "vital_value_numeric"]),
}


def _section(title):
    print("\n" + "=" * 72 + f"\n{title}\n" + "=" * 72)


def load_config():
    cfg_path = os.environ.get("CLIF_CONFIG") or str(HERE.parent / "config" / "config.json")
    with open(cfg_path) as f:
        cfg = json.load(f)
    print(f"  config: {cfg_path}")
    return cfg


def read_table(tables_path, file_type, name, columns=None):
    """Read a CLIF table, trying both `clif_<name>` and `<name>` filenames.

    `columns` limits the read to those columns (fast on huge tables); silently
    falls back to a full read if the column subset is unavailable.
    """
    for stem in (f"clif_{name}", name):
        p = Path(tables_path) / f"{stem}.{file_type}"
        if not p.exists():
            continue
        if file_type == "parquet":
            try:
                return pd.read_parquet(p, columns=columns)
            except Exception:
                return pd.read_parquet(p)
        else:
            try:
                return pd.read_csv(p, usecols=columns)
            except Exception:
                return pd.read_csv(p)
    return None


def _nonnull_value_mask(df, valcols):
    """True where ANY of the candidate value columns is non-null."""
    mask = pd.Series(False, index=df.index)
    for vc in valcols:
        if vc in df.columns:
            mask = mask | df[vc].notna()
    return mask


def part_c_environment():
    _section("C. Environment")
    try:
        import polars as pl
        print("  polars :", pl.__version__)
    except Exception as e:  # pragma: no cover
        print("  polars : (not importable)", e)
    try:
        import clifpy
        print("  clifpy :", getattr(clifpy, "__version__", "?"))
    except Exception as e:
        print("  clifpy : (not importable)", e)
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(HERE.parent), "log", "-1", "--oneline"]
        ).decode().strip()
        print("  commit :", commit)
    except Exception:
        pass


def part_a_component_nulls(tables_path, file_type, timezone):
    _section("A. SOFA component null rates on the CRRT cohort")
    print("  (the component with a HIGH null rate is what blocks the total SOFA;")
    print("   this recomputes SOFA and may take a few minutes)")
    try:
        import polars as pl
        from sofa_calculator import compute_sofa_polars

        crrt = read_table(tables_path, file_type, "crrt_therapy",
                          columns=["hospitalization_id", "recorded_dttm"])
        if crrt is None:
            raise FileNotFoundError("crrt_therapy table not found")
        crrt = crrt.dropna(subset=["recorded_dttm"]).copy()
        crrt["hospitalization_id"] = crrt["hospitalization_id"].astype(str)
        crrt["recorded_dttm"] = pd.to_datetime(crrt["recorded_dttm"])
        init = crrt.groupby("hospitalization_id", as_index=False)["recorded_dttm"].min()
        init["start_dttm"] = init["recorded_dttm"] - pd.Timedelta(hours=12)
        init["end_dttm"] = init["recorded_dttm"] + pd.Timedelta(hours=3)
        cohort = init[["hospitalization_id", "start_dttm", "end_dttm"]]
        print(f"  CRRT hospitalizations in cohort: {len(cohort):,}")

        sofa = compute_sofa_polars(
            data_directory=tables_path,
            cohort_df=pl.from_pandas(cohort),
            filetype=file_type,
            id_name="hospitalization_id",
            fill_na_scores_with_zero=False,   # keep component nulls VISIBLE
            remove_outliers=True,
            timezone=timezone,
        ).to_pandas()

        comp = [c for c in
                ["sofa_renal", "sofa_coag", "sofa_liver", "sofa_cns", "sofa_cv_97", "sofa_resp"]
                if c in sofa.columns]
        if not comp:
            print(f"  !! no sofa component columns found. Columns: {list(sofa.columns)}")
            return
        rates = sofa[comp].isna().mean().sort_values(ascending=False)
        print(f"  SOFA rows computed: {len(sofa):,}")
        print("  component null rate (fraction of cohort with that component missing):")
        for c, r in rates.items():
            field = SOFA_INPUTS.get(c, ("?",))[0]
            flag = "   <== LIKELY CULPRIT" if r >= 0.30 else ""
            print(f"    {c:11s} {r:6.1%}   needs: {field}{flag}")
    except Exception as e:
        print(f"  [Part A could not run: {type(e).__name__}: {e}]")
        print("  -> rely on Part B below to find the missing input.")


def part_b_source_coverage(tables_path, file_type):
    _section("B. Source-table coverage of SOFA inputs")
    tables = {}
    for tbl, (catcol, valcols) in _TABLE_META.items():
        df = read_table(tables_path, file_type, tbl,
                        columns=["hospitalization_id", catcol] + valcols)
        tables[tbl] = (df, catcol, valcols)
        if df is None:
            print(f"\n  [{tbl}]  TABLE NOT FOUND")
            continue
        nh = df["hospitalization_id"].nunique() if "hospitalization_id" in df.columns else "?"
        print(f"\n  [{tbl}]  rows={len(df):,}  hospitalizations={nh}")
        if catcol not in df.columns:
            print(f"    !! expected column '{catcol}' missing. Columns present: {list(df.columns)}")
            continue
        print(f"    distinct {catcol} values (top 60 by count):")
        for k, v in df[catcol].value_counts().head(60).items():
            print(f"       {repr(k):42s} {v:,}")

    print("\n  SOFA input coverage — hospitalizations with >=1 NON-NULL value:")
    for comp, (field, tbl) in SOFA_INPUTS.items():
        if tbl is None:
            continue  # respiratory uses several tables; skip the single-field check
        df, catcol, valcols = tables.get(tbl, (None, None, None))
        if df is None or catcol not in df.columns:
            print(f"    {field:16s} ({tbl}): table/column unavailable")
            continue
        sub = df[df[catcol].astype(str).str.strip().str.lower() == field]
        nn = sub.loc[_nonnull_value_mask(sub, valcols), "hospitalization_id"].nunique() if len(sub) else 0
        tot = df["hospitalization_id"].nunique() or 1
        print(f"    {field:16s} ({tbl}): {nn:,}/{tot:,} hosp = {nn/tot:5.1%}")


def main():
    print("SOFA coverage diagnostic")
    cfg = load_config()
    tables_path = cfg["tables_path"]
    file_type = cfg.get("file_type", "parquet")
    timezone = cfg.get("timezone", "UTC")
    print(f"  tables_path: {tables_path}")
    print(f"  file_type:   {file_type}")

    part_c_environment()
    part_b_source_coverage(tables_path, file_type)   # fast: pinpoints the missing input
    part_a_component_nulls(tables_path, file_type, timezone)  # slower: confirms on the SOFA window

    _section("DONE — copy this ENTIRE console output back to the coordinating center")


if __name__ == "__main__":
    main()
