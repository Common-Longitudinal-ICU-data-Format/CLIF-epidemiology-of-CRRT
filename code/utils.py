import pandas as pd
import numpy as np
from typing import Dict, Tuple
from pathlib import Path
import json


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
