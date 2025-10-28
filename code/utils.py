import pandas as pd
import numpy as np
from typing import Union, Dict, Set, Tuple
from pathlib import Path
import json


def filter_adult_encounters(
    all_encounters: pd.DataFrame,
    age_col: str = "age_at_admission",
    admission_col: str = "admission_dttm",
    min_age: int = 18,
    start_year: int = 2018,
    end_year: int = 2024
) -> pd.DataFrame:
    """
    Filter encounters for adult patients within specified date range.

    Parameters
    ----------
    all_encounters : pd.DataFrame
        DataFrame with hospitalization data
    age_col : str
        Column name for age at admission
    admission_col : str
        Column name for admission datetime
    min_age : int
        Minimum age for inclusion (default 18)
    start_year : int
        Start year for inclusion (default 2018)
    end_year : int
        End year for inclusion (default 2024)

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with adult encounters
    """
    adult_encounters = all_encounters[
        (all_encounters[age_col] >= min_age) &
        (all_encounters[age_col].notna()) &
        (all_encounters[admission_col].dt.year >= start_year) &
        (all_encounters[admission_col].dt.year <= end_year)
    ].copy()

    return adult_encounters


def identify_esrd_hospitalizations(
    hospital_diagnosis_df: pd.DataFrame,
    esrd_codes: list = None
) -> Tuple[Set, Set]:
    """
    Identify hospitalizations with ESRD diagnosis present on admission.

    Parameters
    ----------
    hospital_diagnosis_df : pd.DataFrame
        Hospital diagnosis DataFrame with 'diagnosis_code' and 'poa_present' columns
    esrd_codes : list, optional
        List of ESRD ICD codes. If None, uses default codes.

    Returns
    -------
    Tuple[Set, Set]
        (hospitalization_ids_with_esrd, encounter_blocks_with_esrd)
    """
    if esrd_codes is None:
        esrd_codes = [
            'z992', 'z9115', 'i120', 'n186', 'i132', 'z991158', 'i1311',
            '5856', '40391', '40311', 'v4511', 'v4512'
        ]

    # Filter for ESRD codes present on admission or unknown
    esrd_mask = (
        hospital_diagnosis_df['diagnosis_code'].isin(esrd_codes) &
        ((hospital_diagnosis_df['poa_present'] == 1) |
         (hospital_diagnosis_df['poa_present'].isna()))
    )

    hosp_ids_with_esrd = set(hospital_diagnosis_df[esrd_mask]['hospitalization_id'].unique())
    blocks_with_esrd = set(hospital_diagnosis_df[esrd_mask]['encounter_block'].unique())

    return hosp_ids_with_esrd, blocks_with_esrd


def identify_aki_hospitalizations(
    hospital_diagnosis_df: pd.DataFrame,
    cohort_blocks: Set,
    aki_codes: list = None
) -> Tuple[Set, float]:
    """
    Identify hospitalizations with AKI diagnosis and calculate percentage.

    Parameters
    ----------
    hospital_diagnosis_df : pd.DataFrame
        Hospital diagnosis DataFrame
    cohort_blocks : Set
        Set of encounter blocks in the cohort
    aki_codes : list, optional
        List of AKI ICD codes. If None, uses default codes.

    Returns
    -------
    Tuple[Set, float]
        (encounter_blocks_with_aki, percentage_with_aki)
    """
    if aki_codes is None:
        aki_codes = [
            'n170', 'n171', 'n172', 'n178', 'n179',
            'r34', 'n990', 't795', '5845', '5849', '5848'
        ]

    # Filter to cohort encounters
    cohort_diagnoses = hospital_diagnosis_df[
        hospital_diagnosis_df['encounter_block'].isin(cohort_blocks)
    ]

    # Identify AKI
    aki_mask = cohort_diagnoses['diagnosis_code'].isin(aki_codes)
    blocks_with_aki = set(cohort_diagnoses[aki_mask]['encounter_block'].unique())

    # Calculate percentage
    total_blocks = len(cohort_blocks)
    aki_percentage = (len(blocks_with_aki) / total_blocks * 100) if total_blocks > 0 else 0

    return blocks_with_aki, aki_percentage


def calculate_icu_percentage(
    adt_df: pd.DataFrame,
    cohort_blocks: Set
) -> float:
    """
    Calculate percentage of encounters with ICU stay.

    Parameters
    ----------
    adt_df : pd.DataFrame
        ADT DataFrame with location_category column
    cohort_blocks : Set
        Set of encounter blocks in the cohort

    Returns
    -------
    float
        Percentage of encounters with at least one ICU stay
    """
    cohort_adt = adt_df[adt_df['encounter_block'].isin(cohort_blocks)]
    icu_blocks = cohort_adt[cohort_adt['location_category'] == 'icu']['encounter_block'].unique()

    total_blocks = len(cohort_blocks)
    icu_percentage = (len(icu_blocks) / total_blocks * 100) if total_blocks > 0 else 0

    return icu_percentage


def save_intermediate_output(
    data: Union[pd.DataFrame, dict],
    filename: str,
    output_dir: Union[str, Path] = "output/intermediate"
) -> Path:
    """
    Save intermediate outputs (DataFrames or dicts) to the intermediate directory.

    Parameters
    ----------
    data : pd.DataFrame or dict
        Data to save
    filename : str
        Output filename (with extension)
    output_dir : str or Path
        Output directory path

    Returns
    -------
    Path
        Path to the saved file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / filename

    if isinstance(data, pd.DataFrame):
        if filename.endswith('.parquet'):
            data.to_parquet(output_path, index=False)
        elif filename.endswith('.csv'):
            data.to_csv(output_path, index=False)
        else:
            raise ValueError(f"Unsupported file extension for DataFrame: {filename}")
    elif isinstance(data, dict):
        if filename.endswith('.json'):
            with open(output_path, 'w') as f:
                json.dump(data, f, indent=2, default=str)
        else:
            raise ValueError(f"Unsupported file extension for dict: {filename}")
    else:
        raise ValueError(f"Unsupported data type: {type(data)}")

    return output_path


def create_consort_diagram(
    strobe_counts: Dict,
    output_dir: Union[str, Path] = "../output/final/graphs"
) -> Path:
    """
    Create CONSORT flow diagram for CRRT cohort selection.

    Parameters
    ----------
    strobe_counts : dict
        Dictionary containing counts for each step of cohort selection.
        Expected keys:
        - '0_total_hospitalizations': Total hospitalizations
        - '1_adult_hospitalizations': Adult hospitalizations
        - '1b_after_stitching': Encounter blocks after stitching
        - '2_crrt_blocks': CRRT encounter blocks
        - '3_encounter_blocks_with_esrd': Blocks with ESRD
        - '3_encounter_blocks_without_esrd': Blocks without ESRD
        - '4_encounter_blocks_with_weight': Blocks with weight data
        - '5_encounter_blocks_with_crrt_settings': Blocks with CRRT settings
        - '6_encounter_blocks_with_AKI_no_esrd': Blocks with AKI (for display)
        - '6_Percentage_non_ESRD_encounter_blocks_with_AKI_codes': AKI %
        - '6_number_hosp_without_ICU_stay': Blocks without ICU (for ICU % calc)
    output_dir : str or Path
        Directory to save the diagram (default: ../output/final/graphs)

    Returns
    -------
    Path
        Path to the saved diagram file
    """
    from matplotlib.patches import FancyBboxPatch
    import matplotlib.pyplot as plt

    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Calculate percentages from strobe_counts
    aki_percentage = strobe_counts.get('6_Percentage_non_ESRD_encounter_blocks_with_AKI_codes', 0)
    
    non_esrd_blocks = strobe_counts.get('3_encounter_blocks_without_esrd', 0)
    no_icu_blocks = strobe_counts.get('6_number_hosp_without_ICU_stay', 0)
    icu_percentage = 0
    if non_esrd_blocks > 0:
        icu_percentage = ((non_esrd_blocks - no_icu_blocks) / non_esrd_blocks) * 100

    # Calculate exclusion counts dynamically
    excluded_esrd = strobe_counts.get('2_crrt_blocks', 0) - strobe_counts.get('3_encounter_blocks_without_esrd', 0)
    excluded_no_weight = strobe_counts.get('3_encounter_blocks_without_esrd', 0) - strobe_counts.get('4_encounter_blocks_with_weight', 0)
    excluded_no_crrt_settings = strobe_counts.get('4_encounter_blocks_with_weight', 0) - strobe_counts.get('5_encounter_blocks_with_crrt_settings', 0)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 20))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.5, 0.97, 'CONSORT Flow Diagram', fontsize=16, fontweight='bold', ha='center')
    ax.text(0.5, 0.94, 'CLIF-CRRT Cohort Selection', fontsize=12, ha='center')

    # Define box positions
    BOX_W_MAIN = 0.35
    BOX_H = 0.055
    BOX_H_FINAL = 0.07
    X_MAIN = 0.325
    X_EXCL = 0.75

    # Y positions (properly spaced to avoid overlap)
    y1 = 0.87   # Total hospitalizations
    y2 = 0.77   # Adult hospitalizations
    y3 = 0.67   # After stitching
    y4 = 0.57   # CRRT hospitalizations
    y5 = 0.47   # After ESRD exclusion
    y6 = 0.37   # After weight exclusion
    y7 = 0.27   # With CRRT settings
    y8 = 0.15   # AKI/ICU info boxes
    yF = 0.04   # Final cohort

    # Helper function to draw boxes
    def draw_box(x, y, w, h, text, style='main'):
        if style == 'main':
            fc = 'lightblue'
            ec = 'navy'
        elif style == 'exclude':
            fc = 'mistyrose'
            ec = 'darkred'
        elif style == 'info':
            fc = 'lightyellow'
            ec = 'darkgoldenrod'
        else:  # 'final'
            fc = 'lightgreen'
            ec = 'darkgreen'

        box = FancyBboxPatch((x, y), w, h,
                             boxstyle="round,pad=0.01",
                             linewidth=2, edgecolor=ec, facecolor=fc)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center',
                fontsize=9, wrap=True)

    # Draw main flow boxes
    draw_box(X_MAIN, y1, BOX_W_MAIN, BOX_H,
             f"All Hospitalizations\nn = {strobe_counts.get('0_total_hospitalizations', 0):,}",
             'main')

    draw_box(X_MAIN, y2, BOX_W_MAIN, BOX_H,
             f"Adult Hospitalizations\n(Age ≥ 18, Years 2018-2024)\nn = {strobe_counts.get('1_adult_hospitalizations', 0):,}",
             'main')

    draw_box(X_MAIN, y3, BOX_W_MAIN, BOX_H,
             f"After Encounter Stitching\nEncounter blocks = {strobe_counts.get('1b_after_stitching', 0):,}",
             'main')

    draw_box(X_MAIN, y4, BOX_W_MAIN, BOX_H,
             f"CRRT Encounter Blocks\nn = {strobe_counts.get('2_crrt_blocks', 0):,}",
             'main')

    draw_box(X_MAIN, y5, BOX_W_MAIN, BOX_H,
             f"After ESRD Exclusion\nEncounter blocks = {strobe_counts.get('3_encounter_blocks_without_esrd', 0):,}",
             'main')

    draw_box(X_MAIN, y6, BOX_W_MAIN, BOX_H,
             f"With Weight Data\nEncounter blocks = {strobe_counts.get('4_encounter_blocks_with_weight', 0):,}",
             'main')

    draw_box(X_MAIN, y7, BOX_W_MAIN, BOX_H,
             f"With CRRT Settings\nEncounter blocks = {strobe_counts.get('5_encounter_blocks_with_crrt_settings', 0):,}",
             'main')

    # AKI and ICU info boxes (side by side)
    if aki_percentage > 0:
        draw_box(X_MAIN - 0.02, y8, 0.17, BOX_H,
                 f"With AKI Diagnosis\n{aki_percentage:.1f}%\n" +
                 f"n = {strobe_counts.get('6_encounter_blocks_with_AKI_no_esrd', 0):,}",
                 'info')

    if icu_percentage > 0:
        draw_box(X_MAIN + 0.19, y8, 0.17, BOX_H,
                 f"With ICU Stay\n{icu_percentage:.1f}%",
                 'info')

    draw_box(X_MAIN, yF, BOX_W_MAIN, BOX_H_FINAL,
             f"FINAL COHORT\n" +
             f"Encounter blocks = {strobe_counts.get('5_encounter_blocks_with_crrt_settings', 0):,}",
             'final')

    # Draw exclusion boxes
    draw_box(X_EXCL, y5, 0.2, BOX_H,
             f"Excluded:\nESRD diagnosis\nn = {excluded_esrd:,}",
             'exclude')

    draw_box(X_EXCL, y6, 0.2, BOX_H,
             f"Excluded:\nNo weight data\nn = {excluded_no_weight:,}",
             'exclude')

    draw_box(X_EXCL, y7, 0.2, BOX_H,
             f"Excluded:\nNo CRRT settings\nn = {excluded_no_crrt_settings:,}",
             'exclude')

    # Draw arrows - main flow
    arrow_props = dict(arrowstyle='->', connectionstyle='arc3', lw=2, color='black')

    # Main flow arrows
    ax.annotate('', xy=(X_MAIN + BOX_W_MAIN/2, y2 + BOX_H),
                xytext=(X_MAIN + BOX_W_MAIN/2, y1),
                arrowprops=arrow_props)

    ax.annotate('', xy=(X_MAIN + BOX_W_MAIN/2, y3 + BOX_H),
                xytext=(X_MAIN + BOX_W_MAIN/2, y2),
                arrowprops=arrow_props)

    ax.annotate('', xy=(X_MAIN + BOX_W_MAIN/2, y4 + BOX_H),
                xytext=(X_MAIN + BOX_W_MAIN/2, y3),
                arrowprops=arrow_props)

    ax.annotate('', xy=(X_MAIN + BOX_W_MAIN/2, y5 + BOX_H),
                xytext=(X_MAIN + BOX_W_MAIN/2, y4),
                arrowprops=arrow_props)

    ax.annotate('', xy=(X_MAIN + BOX_W_MAIN/2, y6 + BOX_H),
                xytext=(X_MAIN + BOX_W_MAIN/2, y5),
                arrowprops=arrow_props)

    ax.annotate('', xy=(X_MAIN + BOX_W_MAIN/2, y7 + BOX_H),
                xytext=(X_MAIN + BOX_W_MAIN/2, y6),
                arrowprops=arrow_props)

    ax.annotate('', xy=(X_MAIN + BOX_W_MAIN/2, yF + BOX_H_FINAL),
                xytext=(X_MAIN + BOX_W_MAIN/2, y7),
                arrowprops=arrow_props)

    # Exclusion arrows
    ax.annotate('', xy=(X_EXCL, y5 + BOX_H/2),
                xytext=(X_MAIN + BOX_W_MAIN, y5 + BOX_H/2),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='darkred'))

    ax.annotate('', xy=(X_EXCL, y6 + BOX_H/2),
                xytext=(X_MAIN + BOX_W_MAIN, y6 + BOX_H/2),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='darkred'))

    ax.annotate('', xy=(X_EXCL, y7 + BOX_H/2),
                xytext=(X_MAIN + BOX_W_MAIN, y7 + BOX_H/2),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='darkred'))

    # Save figure
    consort_file = output_dir / "consort_diagram_crrt.png"
    plt.savefig(consort_file, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close()

    return consort_file


def process_crrt_waterfall(
    crrt: pd.DataFrame,
    *,
    id_col: str = "hospitalization_id",
    gap_thresh: Union[str, pd.Timedelta] = "2h",
    infer_modes: bool = True,          # infer missing mode from numeric pattern
    flag_missing_bfr: bool = True,     # add QC flag if blood-flow still NaN
    wipe_unused: bool = True,          # null parameters not used by the mode
    fix_islands: bool = True,          # relabel single-row SCUF islands
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Clean + episode-aware forward-fill for the CLIF `crrt_therapy` table.
    Episode-aware clean-up and forward-fill of the CLIF `crrt_therapy` table.

    The function mirrors the respiratory-support “waterfall” logic but adapts it to
    the quirks of Continuous Renal Replacement Therapy (CRRT):

    1. **Episode detection** - a new `crrt_episode_id` starts whenever  
       • `crrt_mode_category` changes **OR**  
       • the gap between successive rows exceeds *gap_thresh* (default 2 h).
    2. **Numeric forward-fill inside an episode** - fills *only* the parameters
       that are clinically relevant for the active mode.
    3. **Mode-specific wiping** after filling, parameters that are **not used**
       in the current mode (e.g. `dialysate_flow_rate` in SCUF) are nulled so
       stale data never bleed across modes.
    4. **Deduplication & ordering** guarantees exactly **one row per
       `(id_col, recorded_dttm)`**, chronologically ordered.

    Parameters
    ----------
    crrt : pd.DataFrame
        Raw `crrt_therapy` table **in UTC**. Must contain the schema columns
        defined on the CLIF website (see docstring footer).
    id_col : str, default ``"hospitalization_id"``
        Encounter-level identifier.
    gap_thresh : str or pd.Timedelta, default ``"2h"``
        Maximum tolerated gap **inside** an episode before a new episode is
        forced. Accepts any pandas-parsable offset string (``"90min"``, ``"3h"``,
        …) or a ``pd.Timedelta``.
    verbose : bool, default ``True``
        If *True* prints progress banners.

    Returns
    -------
    pd.DataFrame
        Processed CRRT DataFrame with

        * ``crrt_episode_id`` (int32) - sequential per encounter,
        * forward-filled numeric parameters **within** each episode,
        * unused parameters blanked per mode,
        * unique, ordered rows ``id_col, recorded_dttm``.

    Add-ons v2.0
    ------------
    • Optional numeric-pattern inference of `crrt_mode_category`.
    • Flags rows that *should* have blood-flow but don't.
    • Optional fix for single-row modality islands (sandwiched rows).
    • Optional wipe vs. keep of parameters not used by the active mode.

    Key steps
    ----------
    0.  Lower-case strings, coerce numerics, **infer** mode when blank.
    1.  **Relabel single-row SCUF islands** (if *fix_islands*).
    2.  Detect `crrt_episode_id` (mode change or >gap_thresh).
    3.  Forward-fill numeric parameters *within* an episode.
    4.  QC flag → `blood_flow_missing_after_ffill` (optional).
    5.  Wipe / flag parameters not valid for the mode (configurable).
    6.  Deduplicate & order ⇒ one row per ``(id_col, recorded_dttm)``.
    """
    p = print if verbose else (lambda *_, **__: None)
    gap_thresh = pd.Timedelta(gap_thresh)

    # ───────────── Phase 0 — prep, numeric coercion, optional inference
    p("✦ Phase 0: prep & numeric coercion (+optional mode inference)")
    df = crrt.copy()

    df["crrt_mode_category"] = df["crrt_mode_category"].str.lower()
    # save original dialysate_flow_rate values
    df["_orig_df"] = df["dialysate_flow_rate"]

    # 0a) RAW SCUF DF‐OUT sanity check
    # look for rows that are already labeled “scuf”
    # and that have a non‐zero dialysate_flow_rate in the raw data
    raw_scuf = df["crrt_mode_category"].str.lower() == "scuf"
    raw_df_positive = df["_orig_df"].fillna(0) > 0

    n_bad = (raw_scuf & raw_df_positive).sum()
    if n_bad:
        print(f"!!!  Found {n_bad} raw SCUF rows with dialysate_flow_rate > 0 (should be 0 or NA)")
        print(" Converting these mode category to NA, keep recorded numerical values as the ground truth")
        df.loc[raw_df_positive, "crrt_mode_category"] = np.nan
    else:
        print("!!! No raw SCUF rows had dialysate_flow_rate > 0")

    NUM_COLS = [
        "blood_flow_rate",
        "pre_filter_replacement_fluid_rate",
        "post_filter_replacement_fluid_rate",
        "dialysate_flow_rate",
        "ultrafiltration_out",
    ]
    NUM_COLS = [c for c in NUM_COLS if c in df.columns]
    df[NUM_COLS] = df[NUM_COLS].apply(pd.to_numeric, errors="coerce")

    #  any row whose original ultrafiltration_out was >0 must never be SCUF
    def drop_scuf_on_positive_df(df, p):
        bad_df  = df["_orig_df"].fillna(0) > 0
        scuf_now = df["crrt_mode_category"] == "scuf"
        n = (bad_df & scuf_now).sum()
        if n:
            p(f"→ Removing {n:,} SCUF labels on rows with DF>0")
            df.loc[bad_df & scuf_now, "crrt_mode_category"] = np.nan
            

    if infer_modes:
        miss = df["crrt_mode_category"].isna()
        pre  = df["pre_filter_replacement_fluid_rate"].notna()
        post = df["post_filter_replacement_fluid_rate"].notna()
        dial = df["dialysate_flow_rate"].notna()
        bf   = df["blood_flow_rate"].notna()
        uf   = df["ultrafiltration_out"].notna()
        all_num_present = df[NUM_COLS].notna().all(axis=1)

        df.loc[miss & all_num_present,                       "crrt_mode_category"] = "cvvhdf"
        df.loc[miss & (~dial) & pre & post,                  "crrt_mode_category"] = "cvvh"
        df.loc[miss & dial & (~pre) & (~post),               "crrt_mode_category"] = "cvvhd"
        df.loc[miss & (~dial) & (~pre) & (~post) & bf & uf,  "crrt_mode_category"] = "scuf"

        filled = (miss & df["crrt_mode_category"].notna()).sum()
        p(f"  • numeric-pattern inference filled {filled:,} missing modes")
        drop_scuf_on_positive_df(df, p)

    # ───────────── Phase 1 — sort and *fix islands before episodes*
    p("✦ Phase 1: sort + SCUF-island fix")
    df = df.sort_values([id_col, "recorded_dttm"]).reset_index(drop=True)

    if fix_islands:
        # after sorting, BEFORE episode detection
        prev_mode = df.groupby(id_col)["crrt_mode_category"].shift()
        next_mode = df.groupby(id_col)["crrt_mode_category"].shift(-1)

        scuf_island = (
            (df["crrt_mode_category"] == "scuf") &
            (prev_mode.notna()) & (next_mode.notna()) &     # ensure we have neighbours
            (prev_mode == next_mode)                        # both neighbours agree
        )

        df.loc[scuf_island, "crrt_mode_category"] = prev_mode[scuf_island]
        n_fixed = scuf_island.sum()
        p(f"  • relabelled {n_fixed:,} SCUF-island rows")
        drop_scuf_on_positive_df(df, p)


    # ───────────── Phase 2 — episode detection (now with fixed modes)
    p("✦ Phase 2: derive `crrt_episode_id`")
    mode_change = (
        df.groupby(id_col)["crrt_mode_category"]
          .apply(lambda s: s != s.shift())
          .reset_index(level=0, drop=True)
    )
    time_gap = df.groupby(id_col)["recorded_dttm"].diff().gt(gap_thresh).fillna(False)
    df["crrt_episode_id"] = ((mode_change | time_gap)
                              .groupby(df[id_col]).cumsum()
                              .astype("int32"))

    # ───────────── Phase 3 — forward-fill numerics inside episodes
    p("✦ Phase 3: forward-fill numeric vars inside episodes")
    tqdm.pandas(disable=not verbose, desc="ffill per episode")
    df[NUM_COLS] = (
        df.groupby([id_col, "crrt_episode_id"], sort=False, group_keys=False)[NUM_COLS]
          .progress_apply(lambda g: g.ffill())
    )

    # QC: blood-flow still missing?
    if flag_missing_bfr and "blood_flow_rate" in NUM_COLS:
        need_bfr = df["crrt_mode_category"].isin(["scuf", "cvvh", "cvvhd", "cvvhdf"])
        df["blood_flow_missing_after_ffill"] = need_bfr & df["blood_flow_rate"].isna()
        p(f"  • blood-flow still missing where required: "
          f"{df['blood_flow_missing_after_ffill'].mean():.1%}")
        
    # Bridge tiny episodes
    
    single_row_ep = (
        df.groupby([id_col, "crrt_episode_id"]).size() == 1
    ).reset_index(name="n").query("n == 1")
    print("Bridging single row episodes")

    rows_to_bridge = df.merge(single_row_ep[[id_col, "crrt_episode_id"]],
                            on=[id_col, "crrt_episode_id"]).index
    
    CAT_COLS = [c for c in ["crrt_mode_category"] if c in df.columns]

    # Combine with the numeric columns we already had
    BRIDGE_COLS = NUM_COLS + CAT_COLS

    # Forward-fill (and back-fill just in case the island is the first row of the encounter)
    df.loc[rows_to_bridge, BRIDGE_COLS] = (
        df.loc[rows_to_bridge, BRIDGE_COLS]
        .groupby(df.loc[rows_to_bridge, id_col])      # keep encounter boundaries
        .apply(lambda g: g.ffill())          
        .reset_index(level=0, drop=True)
    )
    drop_scuf_on_positive_df(df, p)
    # ───────────── Phase 4 — wipe / flag unused parameters
    p("✦ Phase 4: handle parameters not valid for the mode")
    MODE_PARAM_MAP = {
        "scuf":   {"blood_flow_rate", "ultrafiltration_out"},
        "cvvh":   {"blood_flow_rate", "pre_filter_replacement_fluid_rate",
                   "post_filter_replacement_fluid_rate", "ultrafiltration_out"},
        "cvvhd":  {"blood_flow_rate", "dialysate_flow_rate", "ultrafiltration_out"},
        "cvvhdf": {"blood_flow_rate", "pre_filter_replacement_fluid_rate","post_filter_replacement_fluid_rate",
                   "dialysate_flow_rate", "ultrafiltration_out"},
    }

    wiped_totals = {c: 0 for c in NUM_COLS}
    for mode, keep in MODE_PARAM_MAP.items():
        mask = df["crrt_mode_category"] == mode
        drop_cols = list(set(NUM_COLS) - keep)
        if wipe_unused:
            for col in drop_cols:
                wiped_totals[col] += df.loc[mask, col].notna().sum()
            df.loc[mask, drop_cols] = np.nan
        else:
            for col in drop_cols:
                df.loc[mask & df[col].notna(), f"{col}_unexpected"] = True

    if verbose and wipe_unused:
        p("  • cells set → NA by wipe:")
        for col, n in wiped_totals.items():
            p(f"    {col:<35} {n:>8,}")
    # ───────────── Phase 4a — SCUF‐specific sanity check
    if "dialysate_flow_rate" in df.columns:
        # only consider rows that were originally SCUF mode
        # and whose original _orig_df was non‐zero/non‐NA
        scuf_rows = df["crrt_mode_category"] == "scuf"
        orig_bad = df["_orig_df"].fillna(0) > 0

        # these are rows where the *original* data had UF>0 despite SCUF
        bad_orig_scuf = scuf_rows & orig_bad

        n_bad_orig = bad_orig_scuf.sum()
        if n_bad_orig:
            p(f"!!! {n_bad_orig} rows originally labeled SCUF had DF>0 (raw data); forcing DF→NA for those")
            df.loc[bad_orig_scuf, "dialysate_flow_rate"] = np.nan
        else:
            p("!!! No SCUF rows with DF>0")

    # then drop the helper column
    df = df.drop(columns="_orig_df")

    # ───────────── Phase 5 — deduplicate & order
    p("✦ Phase 5: deduplicate & order")
    pre = len(df)
    df = (
        df.drop_duplicates(subset=[id_col, "recorded_dttm"])
          .sort_values([id_col, "recorded_dttm"])
          .reset_index(drop=True)
    )
    p(f"  • dropped {pre - len(df):,} duplicate rows")

    if verbose:
        sparse = df[NUM_COLS].isna().all(axis=1).mean()
        p(f"  • rows with all NUM_COLS missing: {sparse:.1%}")

    p("[OK] CRRT waterfall complete.")
    return df


def handle_crrt_outliers(
    crrt_df: pd.DataFrame,
    config_path: str = "config/outlier_config.json"
) -> Tuple[pd.DataFrame, Dict]:
    """
    Remove CRRT parameter outliers based on physiologically plausible ranges.

    Parameters
    ----------
    crrt_df : pd.DataFrame
        CRRT data with flow rate parameters
    config_path : str
        Path to outlier_config.json (default: config/outlier_config.json)

    Returns
    -------
    tuple[pd.DataFrame, dict]
        (cleaned_df, outlier_summary)
        - cleaned_df: DataFrame with outliers set to NaN
        - outlier_summary: Dict with outlier counts per parameter
    """
    # Load config
    config_file = Path(config_path)
    if not config_file.exists():
        print(f"⚠️  Config file not found: {config_path}")
        print("   Skipping outlier removal")
        return crrt_df.copy(), {}

    with open(config_file, 'r') as f:
        config = json.load(f)

    cleaned_df = crrt_df.copy()
    outlier_summary = {}

    # CRRT-specific parameters to check
    crrt_params = [
        'blood_flow_rate',
        'pre_filter_replacement_fluid_rate',
        'post_filter_replacement_fluid_rate',
        'dialysate_flow_rate',
        'ultrafiltration_out'
    ]

    print("\n" + "=" * 80)
    print("Removing CRRT Parameter Outliers")
    print("=" * 80)

    for param in crrt_params:
        if param not in cleaned_df.columns:
            continue

        if param not in config:
            print(f"   ⚠️  {param} not in config, skipping")
            continue

        min_val, max_val = config[param]

        # Count outliers
        below_min = (cleaned_df[param] < min_val).sum()
        above_max = (cleaned_df[param] > max_val).sum()
        total_outliers = below_min + above_max
        total_values = cleaned_df[param].notna().sum()

        if total_outliers > 0:
            # Set outliers to NaN
            cleaned_df.loc[cleaned_df[param] < min_val, param] = np.nan
            cleaned_df.loc[cleaned_df[param] > max_val, param] = np.nan

            outlier_summary[param] = {
                'below_min': int(below_min),
                'above_max': int(above_max),
                'total_outliers': int(total_outliers),
                'total_values': int(total_values),
                'percent_outliers': float((total_outliers / total_values * 100) if total_values > 0 else 0),
                'range': [min_val, max_val]
            }

            print(f"\n   {param}:")
            print(f"     Range: {min_val}-{max_val}")
            print(f"     Below min: {below_min:,}")
            print(f"     Above max: {above_max:,}")
            print(f"     Total outliers: {total_outliers:,}/{total_values:,} ({total_outliers/total_values*100:.1f}%)")

    if not outlier_summary:
        print("\n   ✓ No outliers detected")
    else:
        total_removed = sum(s['total_outliers'] for s in outlier_summary.values())
        print(f"\n   Total parameter values set to NaN: {total_removed:,}")

    print("=" * 80)

    return cleaned_df, outlier_summary


def summarize_by_crrt_mode(
    crrt_df: pd.DataFrame,
    outcomes_df: pd.DataFrame = None,
    output_path: str = None
) -> pd.DataFrame:
    """
    Generate summary statistics by CRRT mode category.

    Parameters
    ----------
    crrt_df : pd.DataFrame
        CRRT analysis data with mode, dose, flow rates
    outcomes_df : pd.DataFrame, optional
        Outcomes data with mortality information
    output_path : str, optional
        Path to save CSV output

    Returns
    -------
    pd.DataFrame
        Summary table with statistics by CRRT mode
    """
    print("\n" + "=" * 80)
    print("CRRT Mode-Specific Summary Statistics")
    print("=" * 80)

    if 'crrt_mode_category' not in crrt_df.columns:
        print("⚠️  crrt_mode_category column not found")
        return pd.DataFrame()

    # Merge with outcomes if provided
    if outcomes_df is not None:
        analysis_df = crrt_df.merge(
            outcomes_df[['encounter_block', 'in_hosp_death', 'death_30d']],
            on='encounter_block',
            how='left'
        )
    else:
        analysis_df = crrt_df.copy()

    # Standardize mode names
    analysis_df['mode'] = analysis_df['crrt_mode_category'].str.upper()

    summary_rows = []

    for mode in sorted(analysis_df['mode'].dropna().unique()):
        mode_data = analysis_df[analysis_df['mode'] == mode]

        row = {
            'CRRT Mode': mode,
            'N': len(mode_data),
            'Percent': f"{len(mode_data) / len(analysis_df) * 100:.1f}%"
        }

        # CRRT Dose
        if 'crrt_dose_ml_kg_hr' in mode_data.columns:
            dose = mode_data['crrt_dose_ml_kg_hr'].dropna()
            if len(dose) > 0:
                row['Dose (mL/kg/hr) - Mean±SD'] = f"{dose.mean():.1f}±{dose.std():.1f}"
                row['Dose - Median [IQR]'] = f"{dose.median():.1f} [{dose.quantile(0.25):.1f}-{dose.quantile(0.75):.1f}]"

        # Flow rates
        flow_params = {
            'Blood Flow (mL/min)': 'blood_flow_rate_filled',
            'Dialysate (mL/hr)': 'dialysate_flow_rate',
            'Pre-filter (mL/hr)': 'pre_filter_replacement_fluid_rate',
            'Post-filter (mL/hr)': 'post_filter_replacement_fluid_rate',
            'UF Out (mL/hr)': 'ultrafiltration_out'
        }

        for label, col in flow_params.items():
            if col in mode_data.columns:
                values = mode_data[col].dropna()
                if len(values) > 0:
                    row[f'{label} - Median [IQR]'] = f"{values.median():.0f} [{values.quantile(0.25):.0f}-{values.quantile(0.75):.0f}]"

        # Mortality
        if 'in_hosp_death' in mode_data.columns:
            deaths = mode_data['in_hosp_death'].sum()
            row['In-Hospital Mortality'] = f"{deaths} ({deaths/len(mode_data)*100:.1f}%)"

        if 'death_30d' in mode_data.columns:
            deaths_30d = mode_data['death_30d'].sum()
            row['30-Day Mortality'] = f"{deaths_30d} ({deaths_30d/len(mode_data)*100:.1f}%)"

        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)

    # Display summary
    print("\n" + summary_df.to_string(index=False))

    # Save if path provided
    if output_path:
        summary_df.to_csv(output_path, index=False)
        print(f"\n✓ Summary saved to: {output_path}")

    print("\n" + "=" * 80)

    return summary_df

def create_table_one_competing_risk(
    df: 'pd.DataFrame',
    stratify_by: str = 'outcome',
    output_path = None,
    outcome_labels = None
) -> 'pd.DataFrame':
    """
    Create Table 1 for competing risk analysis with baseline characteristics.

    Parameters
    ----------
    df : pd.DataFrame
        Competing risk dataset with all variables
    stratify_by : str, optional
        Column to stratify by. Default is 'outcome' (0=Censored, 1=Discharged, 2=Died).
        Can also use None for overall only, or 'crrt_mode_category' for mode comparison.
    output_path : str, optional
        Path to save CSV output
    outcome_labels : dict, optional
        Custom labels for outcome categories. Default: {0: 'Censored', 1: 'Discharged', 2: 'Died'}

    Returns
    -------
    pd.DataFrame
        Table 1 with characteristics by strata

    Example
    -------
    >>> table1 = create_table_one_competing_risk(
    ...     competing_risk_final,
    ...     stratify_by='outcome',
    ...     output_path='output/final/table1_by_outcome.csv'
    ... )
    """

    if outcome_labels is None:
        outcome_labels = {0: 'Censored (>90d)', 1: 'Discharged Alive', 2: 'Died'}

    # Define variable groups
    continuous_vars = {
        'Demographics': ['age_at_admission'],
        'CRRT Parameters': [
            'crrt_dose_ml_kg_hr', 'blood_flow_rate_filled',
            'dialysate_flow_rate', 'pre_filter_replacement_fluid_rate',
            'post_filter_replacement_fluid_rate', 'ultrafiltration_out',
            'total_flow_rate', 'weight_kg'
        ],
        'Labs (Peri-CRRT)': [
            'ph_arterial_peri_crrt', 'lactate_peri_crrt', 'bicarbonate_peri_crrt',
            'potassium_peri_crrt', 'sodium_peri_crrt', 'creatinine_peri_crrt',
            'bun_peri_crrt', 'hemoglobin_peri_crrt', 'glucose_serum_peri_crrt'
        ],
        'SOFA Scores': [
            'sofa_cv_97', 'sofa_coag', 'sofa_liver',
            'sofa_resp', 'sofa_cns', 'sofa_renal', 'sofa_total'
        ],
        'Length of Stay': ['icu_los_days', 'hosp_los_days']
    }

    categorical_vars = {
        'Demographics': ['sex_category', 'race_category', 'ethnicity_category'],
        'CRRT Parameters': ['crrt_mode_category'],
        'Data Quality': ['has_any_lab', 'analysis_ready']
    }

    # Helper function for continuous variables
    def format_continuous(series):
        """Format as Mean ± SD, Median [IQR], (n/N with data)"""
        valid = series.dropna()
        n_valid = len(valid)
        n_total = len(series)

        if n_valid == 0:
            return "—"

        mean = valid.mean()
        sd = valid.std()
        median = valid.median()
        q25 = valid.quantile(0.25)
        q75 = valid.quantile(0.75)

        return f"{mean:.1f} ± {sd:.1f}, {median:.1f} [{q25:.1f}-{q75:.1f}] ({n_valid}/{n_total})"

    # Helper function for categorical variables
    def format_categorical(series):
        """Format as n (%), with missingness noted"""
        counts = series.value_counts()
        total = len(series)
        missing = series.isna().sum()

        result = []
        for cat, count in counts.items():
            pct = count / total * 100
            result.append(f"  {cat}: {count} ({pct:.1f}%)")

        if missing > 0:
            result.append(f"  Missing: {missing} ({missing/total*100:.1f}%)")

        return "\n".join(result)

    # Build table
    table_rows = []

    # Overall N
    if stratify_by is not None:
        strata = df[stratify_by].unique()
        strata = sorted([s for s in strata if pd.notna(s)])

        # N row
        n_row = {'Variable': 'N', 'Overall': str(len(df))}
        for stratum in strata:
            stratum_df = df[df[stratify_by] == stratum]
            label = (
                outcome_labels.get(stratum, str(stratum))
                if stratify_by == 'outcome'
                else str(stratum)
            )
            n_row[label] = str(len(stratum_df))
        table_rows.append(n_row)
    else:
        table_rows.append({'Variable': 'N', 'Overall': str(len(df))})

    # Process continuous variables
    for group_name, var_list in continuous_vars.items():
        # Group header
        table_rows.append({'Variable': f'**{group_name}**', 'Overall': ''})

        for var in var_list:
            if var not in df.columns:
                continue

            # Clean variable name
            var_label = var.replace('_', ' ').title()
            if 'peri_crrt' in var:
                var_label = var_label.replace(' Peri Crrt', '')

            row = {'Variable': f"  {var_label}"}

            # Overall
            row['Overall'] = format_continuous(df[var])

            # By strata
            if stratify_by is not None:
                for stratum in strata:
                    stratum_df = df[df[stratify_by] == stratum]
                    label = (
                        outcome_labels.get(stratum, str(stratum))
                        if stratify_by == 'outcome'
                        else str(stratum)
                    )
                    row[label] = format_continuous(stratum_df[var])

            table_rows.append(row)

    # Process categorical variables
    for group_name, var_list in categorical_vars.items():
        # Group header
        group_header = {'Variable': f'**{group_name}**', 'Overall': ''}
        if stratify_by is not None:
            for stratum in strata:
                label = (
                    outcome_labels.get(stratum, str(stratum))
                    if stratify_by == 'outcome'
                    else str(stratum)
                )
                group_header[label] = ''
        table_rows.append(group_header)

        for var in var_list:
            if var not in df.columns:
                continue

            # Clean variable name
            var_label = var.replace('_', ' ').title()

            # Add variable name as subheader
            var_header = {'Variable': f"  {var_label}", 'Overall': ''}
            if stratify_by is not None:
                for stratum in strata:
                    label = (
                        outcome_labels.get(stratum, str(stratum))
                        if stratify_by == 'outcome'
                        else str(stratum)
                    )
                    var_header[label] = ''
            table_rows.append(var_header)

            # Get all unique categories for this variable
            all_categories = df[var].dropna().unique()

            # Add row for each category
            for category in sorted(all_categories, key=str):
                cat_row = {'Variable': f"    {category}"}

                # Overall count
                overall_count = (df[var] == category).sum()
                overall_pct = overall_count / len(df) * 100
                cat_row['Overall'] = f"{overall_count} ({overall_pct:.1f}%)"

                # By strata
                if stratify_by is not None:
                    for stratum in strata:
                        stratum_df = df[df[stratify_by] == stratum]
                        strat_count = (stratum_df[var] == category).sum()
                        strat_pct = strat_count / len(stratum_df) * 100 if len(stratum_df) > 0 else 0
                        label = (
                            outcome_labels.get(stratum, str(stratum))
                            if stratify_by == 'outcome'
                            else str(stratum)
                        )
                        cat_row[label] = f"{strat_count} ({strat_pct:.1f}%)"

                table_rows.append(cat_row)

            # Add missing row if there are missing values
            missing_count = df[var].isna().sum()
            if missing_count > 0:
                miss_row = {'Variable': '    Missing'}
                miss_pct = missing_count / len(df) * 100
                miss_row['Overall'] = f"{missing_count} ({miss_pct:.1f}%)"

                if stratify_by is not None:
                    for stratum in strata:
                        stratum_df = df[df[stratify_by] == stratum]
                        strat_miss = stratum_df[var].isna().sum()
                        strat_miss_pct = strat_miss / len(stratum_df) * 100 if len(stratum_df) > 0 else 0
                        label = (
                            outcome_labels.get(stratum, str(stratum))
                            if stratify_by == 'outcome'
                            else str(stratum)
                        )
                        miss_row[label] = f"{strat_miss} ({strat_miss_pct:.1f}%)"

                table_rows.append(miss_row)

    # Time-to-event (only if stratifying by outcome)
    if stratify_by == 'outcome':
        table_rows.append({'Variable': '**Time-to-Event**', 'Overall': ''})

        row = {'Variable': '  Time to Event (days)'}
        row['Overall'] = format_continuous(df['time_to_event_90d'])

        for stratum in strata:
            stratum_df = df[df[stratify_by] == stratum]
            label = outcome_labels.get(stratum, str(stratum))
            row[label] = format_continuous(stratum_df['time_to_event_90d'])

        table_rows.append(row)

    # Convert to DataFrame
    table1_df = pd.DataFrame(table_rows)

    # Save if path provided
    if output_path is not None:
        table1_df.to_csv(output_path, index=False)
        print(f"Table 1 saved to: {output_path}")

        # Also save as HTML
        html_path = output_path.replace('.csv', '.html')

        # Create styled HTML
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Table 1 - Baseline Characteristics</title>
    <style>
        body {{
            font-family: 'Arial', sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #333;
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 13px;
        }}
        th {{
            background-color: #4CAF50;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
            border: 1px solid #ddd;
        }}
        td {{
            padding: 10px;
            border: 1px solid #ddd;
            vertical-align: top;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        tr:hover {{
            background-color: #f0f0f0;
        }}
        .group-header {{
            font-weight: bold;
            background-color: #e8f5e9 !important;
            font-size: 14px;
        }}
        .variable-name {{
            padding-left: 20px;
        }}
        .footer {{
            margin-top: 20px;
            padding-top: 10px;
            border-top: 1px solid #ddd;
            font-size: 11px;
            color: #666;
        }}
        pre {{
            white-space: pre-wrap;
            word-wrap: break-word;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Table 1: Baseline Characteristics - Competing Risk Analysis</h1>
        <p><strong>Stratified by:</strong> {stratify_by if stratify_by else 'Overall (no stratification)'}</p>
        <p><em>Continuous variables shown as: Mean ± SD, Median [IQR] (n/N)</em></p>
        <p><em>Categorical variables shown as: n (%)</em></p>

        {table1_df.to_html(index=False, escape=False, classes='table')}

        <div class="footer">
            <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>CRRT Competing Risk Analysis Pipeline</p>
        </div>
    </div>
</body>
</html>
"""

        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"Table 1 HTML saved to: {html_path}")

    return table1_df


# Optional: Simplified version for quick viewing
def print_table_one_summary(df: pd.DataFrame) -> None:
    """
    Quick summary statistics for competing risk dataset.

    Parameters
    ----------
    df : pd.DataFrame
        Competing risk dataset
    """
    print("=" * 80)
    print("Table 1 Summary - Competing Risk Analysis")
    print("=" * 80)

    print("\n" + "-" * 80)
    print("OUTCOME DISTRIBUTION")
    print("-" * 80)
    outcome_labels = {0: 'Censored (>90d)', 1: 'Discharged Alive', 2: 'Died'}
    for outcome_val in sorted(df['outcome'].unique()):
        count = (df['outcome'] == outcome_val).sum()
        pct = count / len(df) * 100
        print(f"  {outcome_labels.get(outcome_val, outcome_val)}: {count:,} ({pct:.1f}%)")

    print("\n" + "-" * 80)
    print("KEY DEMOGRAPHICS")
    print("-" * 80)
    print(f"Age (years): {df['age_at_admission'].mean():.1f} ± {df['age_at_admission'].std():.1f}")
    print(f"  Median [IQR]: {df['age_at_admission'].median():.1f} [{df['age_at_admission'].quantile(0.25):.1f}-{df['age_at_admission'].quantile(0.75):.1f}]")

    if 'sex_category' in df.columns:
        print("\nSex:")
        print(df['sex_category'].value_counts().to_string())

    print("\n" + "-" * 80)
    print("CRRT PARAMETERS")
    print("-" * 80)
    print("\nCRRT Mode:")
    print(df['crrt_mode_category'].value_counts().to_string())

    print(f"\nCRRT Dose (mL/kg/hr):")
    dose_valid = df['crrt_dose_ml_kg_hr'].dropna()
    print(f"  Mean ± SD: {dose_valid.mean():.1f} ± {dose_valid.std():.1f}")
    print(f"  Median [IQR]: {dose_valid.median():.1f} [{dose_valid.quantile(0.25):.1f}-{dose_valid.quantile(0.75):.1f}]")
    print(f"  Available: {len(dose_valid)}/{len(df)} ({len(dose_valid)/len(df)*100:.1f}%)")

    print("\n" + "-" * 80)
    print("SOFA SCORES")
    print("-" * 80)
    if 'sofa_total' in df.columns:
        sofa_valid = df['sofa_total'].dropna()
        print(f"Total SOFA: {sofa_valid.mean():.1f} ± {sofa_valid.std():.1f}")
        print(f"  Median [IQR]: {sofa_valid.median():.1f} [{sofa_valid.quantile(0.25):.1f}-{sofa_valid.quantile(0.75):.1f}]")
        print(f"  Available: {len(sofa_valid)}/{len(df)} ({len(sofa_valid)/len(df)*100:.1f}%)")

    print("\n" + "-" * 80)
    print("DATA COMPLETENESS")
    print("-" * 80)
    print(f"Has any lab: {df['has_any_lab'].sum():,} ({df['has_any_lab'].sum()/len(df)*100:.1f}%)")
    print(f"Analysis ready: {df['analysis_ready'].sum():,} ({df['analysis_ready'].sum()/len(df)*100:.1f}%)")

    print("\n" + "=" * 80)

