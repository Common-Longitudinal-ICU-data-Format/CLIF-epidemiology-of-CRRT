"""
05_severity_analysis.py
Cross-site severity analysis for CRRT epidemiology.

Compares baseline severity metrics (SOFA, NEE, lactate, mortality) across sites.
Generates forest plots and heatmaps.

Output: all_site_data/severity_analysis.html (standalone report)
"""

import base64
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ── Config ──────────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent.parent / "all_site_data"
OUTPUT = ROOT / "severity_analysis.html"

SITE_LABELS = {
    "mimic_iv": "MIMIC-IV",
    "MIMIC": "MIMIC-IV",
    "nu": "Northwestern",
    "NU": "Northwestern",
    "rush": "Rush",
    "Rush": "Rush",
    "ucmc": "UChicago",
    "UCMC": "UChicago",
    "umich": "UMichigan",
    "UMich": "UMichigan",
    "upenn": "UPenn",
    "penn": "UPenn",
    "ohsu": "OHSU",
    "OHSU": "OHSU",
}

SITE_COLORS = {
    "mimic_iv": "#0072B2",
    "nu": "#E69F00",
    "rush": "#009E73",
    "ucmc": "#CC79A7",
    "umich": "#882255",
    "upenn": "#332288",
    "ohsu": "#44AA99",
}

# Variable names as they appear in table1_crrt_long.csv
VAR_SOFA_TOTAL = "sofa_total, median [Q1,Q3]"
VAR_SOFA_CV = "sofa_cv_97, median [Q1,Q3]"
VAR_SOFA_COAG = "sofa_coag, median [Q1,Q3]"
VAR_SOFA_RENAL = "sofa_renal, median [Q1,Q3]"
VAR_SOFA_LIVER = "sofa_liver, median [Q1,Q3]"
VAR_SOFA_RESP = "sofa_resp, median [Q1,Q3]"
VAR_SOFA_CNS = "sofa_cns, median [Q1,Q3]"
VAR_NEE = "NEE (mcg/kg/min), median [Q1,Q3]"
VAR_LACTATE = "lactate, median [Q1,Q3]"
VAR_VASO_DUR = "Duration of vasopressor support (hours), median [Q1,Q3]"
VAR_MORTALITY = "Mortality overall"

SOFA_COMPONENTS = {
    VAR_SOFA_CV: "CV",
    VAR_SOFA_COAG: "Coag",
    VAR_SOFA_RENAL: "Renal",
    VAR_SOFA_LIVER: "Liver",
    VAR_SOFA_RESP: "Resp",
    VAR_SOFA_CNS: "CNS",
}


# ── Data Loading ────────────────────────────────────────────────────────────

def discover_sites() -> list[Path]:
    """Return sorted list of site dirs that contain final/ outputs."""
    return sorted(
        d for d in ROOT.iterdir()
        if d.is_dir() and (d / "final").is_dir()
    )


def load_all_table1(sites: list[Path]) -> pd.DataFrame:
    """Load and concatenate table1_crrt_long.csv from all sites."""
    frames = []
    for site_dir in sites:
        csv_path = site_dir / "final" / "table1_crrt_long.csv"
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            df["site_dir"] = site_dir.name  # directory name for colors
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def normalize_site_id(site_col: str, site_dir: str) -> str:
    """Map CSV site column value to canonical site_dir name."""
    mapping = {
        "MIMIC": "mimic_iv", "mimic_iv": "mimic_iv",
        "NU": "nu", "nu": "nu",
        "Rush": "rush", "rush": "rush",
        "UCMC": "ucmc", "ucmc": "ucmc",
        "UMich": "umich", "umich": "umich",
        "penn": "upenn", "upenn": "upenn", "UPenn": "upenn",
        "OHSU": "ohsu", "ohsu": "ohsu",
    }
    return mapping.get(site_col, site_dir)


def extract_severity_data(df: pd.DataFrame) -> pd.DataFrame:
    """Extract baseline severity metrics per site into a summary DataFrame."""
    baseline = df[df["subgroup"] == "Baseline"].copy()

    # Get unique sites by directory name
    site_dirs = baseline["site_dir"].unique()
    rows = []

    for sd in site_dirs:
        site_df = baseline[baseline["site_dir"] == sd]
        site_name = site_df["site"].iloc[0] if "site" in site_df.columns else sd
        sid = normalize_site_id(site_name, sd)

        row = {"site_dir": sid, "site_label": SITE_LABELS.get(sid, sid)}

        # Patient count
        n_row = site_df[(site_df["variable"] == "n_patients") & (site_df["stat_type"] == "count")]
        row["n_patients"] = int(n_row["n"].iloc[0]) if len(n_row) > 0 else np.nan

        # Mortality rate
        mort = site_df[(site_df["variable"] == VAR_MORTALITY) & (site_df["level"].astype(str) == "1")]
        if len(mort) > 0:
            row["mortality_n"] = int(mort["n"].iloc[0])
            row["mortality_total"] = int(mort["total"].iloc[0])
            row["mortality_rate"] = row["mortality_n"] / row["mortality_total"]
        else:
            row["mortality_n"] = row["mortality_total"] = np.nan
            row["mortality_rate"] = np.nan

        # Continuous variables: median, q25, q75
        for var_name, col_prefix in [
            (VAR_SOFA_TOTAL, "sofa_total"),
            (VAR_NEE, "nee"),
            (VAR_LACTATE, "lactate"),
            (VAR_VASO_DUR, "vaso_dur"),
        ]:
            v = site_df[(site_df["variable"] == var_name) & (site_df["stat_type"] == "continuous")]
            if len(v) > 0:
                row[f"{col_prefix}_median"] = v["median"].iloc[0]
                row[f"{col_prefix}_q25"] = v["q25"].iloc[0]
                row[f"{col_prefix}_q75"] = v["q75"].iloc[0]
                row[f"{col_prefix}_n"] = int(v["n"].iloc[0])
            else:
                row[f"{col_prefix}_median"] = row[f"{col_prefix}_q25"] = row[f"{col_prefix}_q75"] = np.nan
                row[f"{col_prefix}_n"] = 0

        # SOFA components
        for var_name, comp in SOFA_COMPONENTS.items():
            v = site_df[(site_df["variable"] == var_name) & (site_df["stat_type"] == "continuous")]
            if len(v) > 0:
                row[f"sofa_{comp}_median"] = v["median"].iloc[0]
            else:
                row[f"sofa_{comp}_median"] = np.nan

        rows.append(row)

    return pd.DataFrame(rows).sort_values("site_label").reset_index(drop=True)


# ── Figure Helpers ──────────────────────────────────────────────────────────

def fig_to_data_uri(fig) -> str:
    """Render a matplotlib figure to a base64 PNG data URI."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def get_color(sid: str) -> str:
    return SITE_COLORS.get(sid, "#666666")


# ── Figures ─────────────────────────────────────────────────────────────────

def build_severity_table_html(sev: pd.DataFrame) -> str:
    """Build an HTML summary table of severity metrics by site."""
    header = (
        "<table>"
        "<thead><tr>"
        "<th>Site</th><th>N</th><th>Mortality (%)</th>"
        "<th>SOFA Total</th><th>NEE (mcg/kg/min)</th>"
        "<th>Lactate</th><th>Vasopressor Duration (h)</th>"
        "</tr></thead><tbody>"
    )
    rows = []
    for _, r in sev.iterrows():
        mort_pct = f"{r['mortality_rate']*100:.1f}%" if pd.notna(r["mortality_rate"]) else "—"
        sofa = f"{r['sofa_total_median']:.0f} [{r['sofa_total_q25']:.0f}, {r['sofa_total_q75']:.0f}]" if pd.notna(r["sofa_total_median"]) else "—"
        nee = f"{r['nee_median']:.3f} [{r['nee_q25']:.3f}, {r['nee_q75']:.3f}]" if pd.notna(r["nee_median"]) else "—"
        lac = f"{r['lactate_median']:.1f} [{r['lactate_q25']:.1f}, {r['lactate_q75']:.1f}]" if pd.notna(r["lactate_median"]) else "—"
        vd = f"{r['vaso_dur_median']:.1f} [{r['vaso_dur_q25']:.1f}, {r['vaso_dur_q75']:.1f}]" if pd.notna(r["vaso_dur_median"]) else "—"
        n = f"{int(r['n_patients']):,}" if pd.notna(r["n_patients"]) else "—"
        rows.append(f"<tr><td><b>{r['site_label']}</b></td><td>{n}</td><td>{mort_pct}</td><td>{sofa}</td><td>{nee}</td><td>{lac}</td><td>{vd}</td></tr>")
    return header + "\n".join(rows) + "</tbody></table>"


def build_forest_plots(sev: pd.DataFrame, labels: dict | None = None) -> str:
    """4-panel forest plot: mortality, SOFA, NEE, lactate by site."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    def get_label(row):
        if labels:
            return labels.get(row["site_dir"], row["site_label"])
        return row["site_label"]

    # Panel 1: Mortality rate with 95% CI
    ax = axes[0][0]
    sev_sorted = sev.dropna(subset=["mortality_rate"]).sort_values("mortality_rate")
    y_pos = range(len(sev_sorted))
    for i, (_, r) in enumerate(sev_sorted.iterrows()):
        ci_low, ci_high = stats.binom.interval(0.95, int(r["mortality_total"]), r["mortality_rate"])
        ci_low /= r["mortality_total"]
        ci_high /= r["mortality_total"]
        color = get_color(r["site_dir"])
        ax.errorbar(r["mortality_rate"] * 100, i, xerr=[[r["mortality_rate"]*100 - ci_low*100], [ci_high*100 - r["mortality_rate"]*100]],
                     fmt="o", color=color, markersize=8, capsize=4, linewidth=2)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels([get_label(r) for _, r in sev_sorted.iterrows()])
    ax.set_xlabel("Mortality Rate (%)")
    ax.set_title("In-Hospital Mortality by Site")
    ax.grid(True, alpha=0.3, axis="x")

    # Panel 2: SOFA total (median with IQR bars)
    ax = axes[0][1]
    sev_sorted = sev.dropna(subset=["sofa_total_median"]).sort_values("sofa_total_median")
    y_pos = range(len(sev_sorted))
    for i, (_, r) in enumerate(sev_sorted.iterrows()):
        color = get_color(r["site_dir"])
        ax.errorbar(r["sofa_total_median"], i,
                     xerr=[[r["sofa_total_median"] - r["sofa_total_q25"]], [r["sofa_total_q75"] - r["sofa_total_median"]]],
                     fmt="o", color=color, markersize=8, capsize=4, linewidth=2)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels([get_label(r) for _, r in sev_sorted.iterrows()])
    ax.set_xlabel("SOFA Score")
    ax.set_title("Baseline SOFA Total by Site")
    ax.grid(True, alpha=0.3, axis="x")

    # Panel 3: NEE (median with IQR bars)
    ax = axes[1][0]
    sev_sorted = sev.dropna(subset=["nee_median"]).sort_values("nee_median")
    y_pos = range(len(sev_sorted))
    for i, (_, r) in enumerate(sev_sorted.iterrows()):
        color = get_color(r["site_dir"])
        ax.errorbar(r["nee_median"], i,
                     xerr=[[r["nee_median"] - r["nee_q25"]], [r["nee_q75"] - r["nee_median"]]],
                     fmt="o", color=color, markersize=8, capsize=4, linewidth=2)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels([get_label(r) for _, r in sev_sorted.iterrows()])
    ax.set_xlabel("NEE (mcg/kg/min)")
    ax.set_title("Baseline NEE by Site")
    ax.grid(True, alpha=0.3, axis="x")

    # Panel 4: Lactate (median with IQR bars)
    ax = axes[1][1]
    sev_sorted = sev.dropna(subset=["lactate_median"]).sort_values("lactate_median")
    y_pos = range(len(sev_sorted))
    for i, (_, r) in enumerate(sev_sorted.iterrows()):
        color = get_color(r["site_dir"])
        ax.errorbar(r["lactate_median"], i,
                     xerr=[[r["lactate_median"] - r["lactate_q25"]], [r["lactate_q75"] - r["lactate_median"]]],
                     fmt="o", color=color, markersize=8, capsize=4, linewidth=2)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels([get_label(r) for _, r in sev_sorted.iterrows()])
    ax.set_xlabel("Lactate (mmol/L)")
    ax.set_title("Baseline Lactate by Site")
    ax.grid(True, alpha=0.3, axis="x")

    fig.suptitle("Cross-Site Severity Comparison (Baseline)", fontsize=14, y=1.02)
    fig.tight_layout()
    return fig_to_data_uri(fig)


def build_sofa_heatmap(sev: pd.DataFrame, labels: dict | None = None) -> str:
    """Heatmap: sites as rows, 6 SOFA components as columns."""
    components = ["CV", "Coag", "Renal", "Liver", "Resp", "CNS"]
    cols = [f"sofa_{c}_median" for c in components]

    valid = sev.dropna(subset=cols)
    if valid.empty:
        return ""

    site_labels = []
    for _, r in valid.iterrows():
        if labels:
            site_labels.append(labels.get(r["site_dir"], r["site_label"]))
        else:
            site_labels.append(r["site_label"])

    data = valid[cols].values

    fig, ax = plt.subplots(figsize=(10, max(4, len(valid) * 0.8 + 1)))
    im = ax.imshow(data, cmap="YlOrRd", aspect="auto", vmin=0, vmax=4)

    ax.set_xticks(range(len(components)))
    ax.set_xticklabels(components, fontsize=12)
    ax.set_yticks(range(len(site_labels)))
    ax.set_yticklabels(site_labels, fontsize=12)

    # Annotate cells
    for i in range(len(site_labels)):
        for j in range(len(components)):
            val = data[i, j]
            text_color = "white" if val > 2.5 else "black"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=13, fontweight="bold", color=text_color)

    ax.set_title("Baseline SOFA Component Scores by Site", fontsize=14, pad=12)
    fig.colorbar(im, ax=ax, label="SOFA Score", shrink=0.8)
    fig.tight_layout()
    return fig_to_data_uri(fig)


# ── HTML Report ─────────────────────────────────────────────────────────────

def generate_severity_html(sev: pd.DataFrame, labels: dict | None = None) -> str:
    """Generate the full severity analysis HTML content block."""
    table_html = build_severity_table_html(sev)
    forest_uri = build_forest_plots(sev, labels=labels)
    heatmap_uri = build_sofa_heatmap(sev, labels=labels)

    return f"""
    <div class="section">
        <h2>Severity Comparison Table</h2>
        <p><em>Baseline values at CRRT initiation. Continuous variables: median [Q1, Q3].</em></p>
        {table_html}
    </div>
    <div class="section">
        <h2>Forest Plots</h2>
        <p><em>Mortality: point estimate with 95% binomial CI. Others: median with IQR bars.</em></p>
        <div class="figure-block">
            <img src="{forest_uri}" alt="Forest Plots">
        </div>
    </div>
    <div class="section">
        <h2>SOFA Component Heatmap</h2>
        <p><em>Baseline median SOFA component scores by site. Scale: 0-4.</em></p>
        <div class="figure-block">
            <img src="{heatmap_uri}" alt="SOFA Heatmap">
        </div>
    </div>
    """


def build_standalone_report(sev: pd.DataFrame) -> str:
    """Build a complete standalone HTML report."""
    content = generate_severity_html(sev)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>CRRT Severity Analysis — Cross-Site Comparison</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: 'Arial', sans-serif;
            margin: 0; padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1400px; margin: 0 auto;
            background: white; padding: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #333;
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        .section {{ margin-bottom: 40px; }}
        .section h2 {{
            color: #333;
            border-bottom: 2px solid #e0e0e0;
            padding-bottom: 6px;
        }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 12px; font-size: 13px; }}
        th {{ background: #4CAF50; color: white; padding: 12px; text-align: left; border: 1px solid #ddd; }}
        td {{ padding: 10px; border: 1px solid #ddd; vertical-align: top; }}
        tr:nth-child(even) {{ background: #f9f9f9; }}
        tr:hover {{ background: #f0f0f0; }}
        .figure-block {{ margin-bottom: 32px; }}
        .figure-block img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }}
    </style>
</head>
<body>
<div class="container">
    <h1>CRRT Severity Analysis — Cross-Site Comparison</h1>
    {content}
</div>
</body>
</html>"""


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    sites = discover_sites()
    if not sites:
        print("No sites found in", ROOT)
        return

    print(f"Found {len(sites)} sites: {[s.name for s in sites]}")

    print("Loading table1 data...")
    df = load_all_table1(sites)
    if df.empty:
        print("No table1 data found.")
        return

    print("Extracting severity metrics...")
    sev = extract_severity_data(df)
    print("\nSeverity summary:")
    print(sev[["site_label", "n_patients", "mortality_rate", "sofa_total_median", "nee_median", "lactate_median"]].to_string(index=False))

    print("\nGenerating report...")
    html = build_standalone_report(sev)
    OUTPUT.write_text(html, encoding="utf-8")
    print(f"Written to {OUTPUT}")
    print(f"File size: {OUTPUT.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
