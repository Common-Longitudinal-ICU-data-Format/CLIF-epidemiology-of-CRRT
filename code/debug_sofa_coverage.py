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


def _to_utc(s, tz):
    """Return a tz-aware UTC pandas datetime Series (localize naive to `tz`)."""
    s = pd.to_datetime(s, errors="coerce")
    if getattr(s.dtype, "tz", None) is None:
        s = s.dt.tz_localize(tz, ambiguous="NaT", nonexistent="NaT")
    return s.dt.tz_convert("UTC")


def part_d_window_alignment(cfg, tables_path, file_type, timezone):
    """Confirm WHETHER script 02's SOFA window (anchored on crrt_initiation_time)
    actually overlaps the data. Reads the intermediate crrt_initiation.parquet from
    the last pipeline run and compares its anchor to the first CRRT record, then
    checks how many encounter_blocks have vitals in each window.
    """
    _section("D. SOFA window alignment (why script 02 drops encounters)")
    try:
        # locate the intermediate outputs from the last run
        proj = Path(cfg.get("project_root", HERE.parent))
        outdir = Path(cfg.get("output_dir", "output"))
        if not outdir.is_absolute():
            outdir = proj / outdir
        interm = outdir / "intermediate_phi"
        init_f = interm / "crrt_initiation.parquet"
        xwalk_f = interm / "cohort_df.parquet"
        if not init_f.exists() or not xwalk_f.exists():
            print(f"  [need {init_f} and {xwalk_f} from the last run — not found; skipping]")
            return

        init = pd.read_parquet(init_f)  # [encounter_block, crrt_initiation_time]
        init["crrt_initiation_time"] = _to_utc(init["crrt_initiation_time"], timezone)
        xwalk = pd.read_parquet(xwalk_f, columns=["hospitalization_id", "encounter_block"])
        xwalk["hospitalization_id"] = xwalk["hospitalization_id"].astype(str)
        n_blocks = init["encounter_block"].nunique()
        print(f"  SOFA cohort encounter_blocks: {n_blocks:,}")

        # first CRRT record per encounter_block (the debug script's anchor)
        crrt = read_table(tables_path, file_type, "crrt_therapy",
                          columns=["hospitalization_id", "recorded_dttm"])
        crrt["hospitalization_id"] = crrt["hospitalization_id"].astype(str)
        crrt["recorded_dttm"] = _to_utc(crrt["recorded_dttm"], timezone)
        crrt = crrt.merge(xwalk, on="hospitalization_id", how="inner")
        first_crrt = crrt.groupby("encounter_block", as_index=False)["recorded_dttm"].min() \
                         .rename(columns={"recorded_dttm": "first_crrt_time"})

        m = init.merge(first_crrt, on="encounter_block", how="left")
        m["gap_hours"] = (m["crrt_initiation_time"] - m["first_crrt_time"]).dt.total_seconds() / 3600.0
        g = m["gap_hours"].dropna()
        print("\n  gap = crrt_initiation_time - first_CRRT_record (hours):")
        if len(g):
            print(f"    median {g.median():.2f}h | p25 {g.quantile(.25):.2f} | p75 {g.quantile(.75):.2f} "
                  f"| min {g.min():.1f} | max {g.max():.1f}")
            for thr in (3, 15, 24, 72):
                n = int((g.abs() > thr).sum())
                print(f"    |gap| > {thr:>2}h : {n:,} ({n/len(g):.1%})   "
                      f"{'<== these fall outside the -12h/+3h SOFA window' if thr==15 else ''}")

        # Does each window overlap ANY vitals? (vitals are the densest signal)
        vit = read_table(tables_path, file_type, "vitals",
                         columns=["hospitalization_id", "recorded_dttm"])
        vit["hospitalization_id"] = vit["hospitalization_id"].astype(str)
        vit["recorded_dttm"] = _to_utc(vit["recorded_dttm"], timezone)
        vit = vit.merge(xwalk, on="hospitalization_id", how="inner")[["encounter_block", "recorded_dttm"]]

        def _covered(anchor_col):
            w = m.merge(vit, on="encounter_block", how="left")
            lo = w[anchor_col] - pd.Timedelta(hours=12)
            hi = w[anchor_col] + pd.Timedelta(hours=3)
            inwin = w["recorded_dttm"].between(lo, hi)
            return w.loc[inwin, "encounter_block"].nunique()

        cov_init = _covered("crrt_initiation_time")
        cov_first = _covered("first_crrt_time")
        print("\n  encounter_blocks with >=1 vital inside the -12h/+3h window:")
        print(f"    anchored on crrt_initiation_time (what 02 uses): {cov_init:,}/{n_blocks:,} "
              f"= {cov_init/n_blocks:.1%}   <-- should ~match the SOFA row count")
        print(f"    anchored on first CRRT record (what debug uses):  {cov_first:,}/{n_blocks:,} "
              f"= {cov_first/n_blocks:.1%}")
        print("\n  INTERPRETATION:")
        print("    - If the two rows differ a lot -> crrt_initiation_time is mis-anchored;")
        print("      the window misses data the first-CRRT window finds (the bug).")
        print("    - If both are high but SOFA still dropped rows -> the drop is a timezone/")
        print("      precision mismatch in the pipeline's window filter, not the anchor.")
    except Exception as e:
        print(f"  [Part D could not run: {type(e).__name__}: {e}]")


def _resolve_table_path(tables_path, file_type, name):
    for stem in (f"clif_{name}", name):
        p = Path(tables_path) / f"{stem}.{file_type}"
        if p.exists():
            return p
    return None


def part_e_pipeline_repro(cfg, tables_path, file_type, timezone):
    """Reproduce the pipeline's SOFA drop with the ACTUAL loader code path
    (_load_vitals) on the 02-style cohort, and show the datetime dtypes the
    window filter actually compares. This pinpoints the drop that Part D proved
    is NOT a data-availability problem.
    """
    _section("E. Reproduce the pipeline SOFA drop (faithful _load_vitals)")
    try:
        import polars as pl
        from sofa_calculator import _load_vitals, ensure_timezone_lazy

        proj = Path(cfg.get("project_root", HERE.parent))
        outdir = Path(cfg.get("output_dir", "output"))
        if not outdir.is_absolute():
            outdir = proj / outdir
        interm = outdir / "intermediate_phi"
        init = pd.read_parquet(interm / "crrt_initiation.parquet")
        eb = pd.read_parquet(interm / "cohort_df.parquet",
                             columns=["hospitalization_id", "encounter_block"])
        eb["hospitalization_id"] = eb["hospitalization_id"].astype(str)

        # build the SOFA cohort EXACTLY like script 02
        sc = init.merge(eb, on="encounter_block")
        sc["crrt_initiation_time"] = pd.to_datetime(sc["crrt_initiation_time"])
        sc["start_dttm"] = sc["crrt_initiation_time"] - pd.Timedelta(hours=12)
        sc["end_dttm"] = sc["crrt_initiation_time"] + pd.Timedelta(hours=3)
        cohort = pl.from_pandas(sc[["hospitalization_id", "encounter_block", "start_dttm", "end_dttm"]]) \
                   .with_columns(pl.col("hospitalization_id").cast(pl.Utf8))
        n_blocks = cohort["encounter_block"].n_unique()
        hosp_ids = cohort["hospitalization_id"].unique().to_list()

        # faithful pipeline call (cohort window as script 02 passes it)
        vit = _load_vitals(tables_path, file_type, hosp_ids, cohort, timezone).collect()
        surv = vit["encounter_block"].n_unique() if "encounter_block" in vit.columns else 0
        print(f"  cohort encounter_blocks:                         {n_blocks:,}")
        print(f"  _load_vitals AS-IS -> blocks with vitals:        {surv:,}"
              f"  ({surv/n_blocks:.1%})   <-- SOFA row count tracks this")

        # candidate FIX: cast the cohort window to microseconds so both sides of the
        # >= / <= filter share the same time unit as the data.
        cohort_us = cohort.with_columns([
            pl.col("start_dttm").dt.cast_time_unit("us"),
            pl.col("end_dttm").dt.cast_time_unit("us"),
        ])
        vit2 = _load_vitals(tables_path, file_type, hosp_ids, cohort_us, timezone).collect()
        surv2 = vit2["encounter_block"].n_unique() if "encounter_block" in vit2.columns else 0
        print(f"  _load_vitals with us-cast window -> blocks:      {surv2:,}"
              f"  ({surv2/n_blocks:.1%})   <-- if this jumps to ~100%, THAT is the fix")
        print(f"  (Part D, correct handling, found {n_blocks:,} have vitals in-window)")

        # show the two dtypes the window filter compares
        vf = _resolve_table_path(tables_path, file_type, "vitals")
        if vf is not None:
            samp = pl.scan_parquet(vf).select("recorded_dttm").head(2000)
            raw_dt = samp.collect_schema()["recorded_dttm"]
            post_dt = samp.with_columns(ensure_timezone_lazy(pl.col("recorded_dttm"), timezone)) \
                          .collect_schema()["recorded_dttm"]
            print(f"\n  recorded_dttm dtype: raw {raw_dt}")
            print(f"                       after ensure_timezone_lazy -> {post_dt}")
            print(f"  cohort start_dttm dtype:                     {cohort['start_dttm'].dtype}")
            print("  (a timezone/unit mismatch between the last two silently breaks the >= / <= filter)")
    except Exception as e:
        print(f"  [Part E could not run: {type(e).__name__}: {e}]")


def main():
    print("SOFA coverage diagnostic")
    cfg = load_config()
    tables_path = cfg["tables_path"]
    file_type = cfg.get("file_type", "parquet")
    timezone = cfg.get("timezone", "UTC")
    print(f"  tables_path: {tables_path}")
    print(f"  file_type:   {file_type}")

    # `--window` runs ONLY the environment + window-alignment check (fast: reuses
    # the last run's intermediates, no SOFA recompute). Default runs everything.
    window_only = "--window" in sys.argv

    part_c_environment()
    part_d_window_alignment(cfg, tables_path, file_type, timezone)  # pinpoints the window/anchor drop
    part_e_pipeline_repro(cfg, tables_path, file_type, timezone)    # reproduces drop via real loader
    if not window_only:
        part_b_source_coverage(tables_path, file_type)       # pinpoints a missing input
        part_a_component_nulls(tables_path, file_type, timezone)  # confirms on the SOFA window

    _section("DONE — copy this ENTIRE console output back to the coordinating center")


if __name__ == "__main__":
    main()
