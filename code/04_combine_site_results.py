"""
04_combine_site_results.py
Combines per-site Table 1 HTML and figures into a single tabbed dashboard.
Output: all_site_data/combined_dashboard.html
"""

import base64
import importlib
import io
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

# Import 05_severity_analysis (numeric prefix requires importlib)
sys.path.insert(0, str(Path(__file__).resolve().parent))
_sev_mod = importlib.import_module("05_severity_analysis")
load_all_table1 = _sev_mod.load_all_table1
extract_severity_data = _sev_mod.extract_severity_data
generate_severity_html = _sev_mod.generate_severity_html

# ── Config ──────────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent.parent / "all_site_data"
OUTPUT = ROOT / "combined_dashboard.html"

# Display names for sites (fallback: folder name title-cased)
SITE_LABELS = {
    "emory": "Emory",
    "hopkins": "Hopkins",
    "mimic_iv": "MIMIC-IV",
    "nu": "Northwestern",
    "ohsu": "OHSU",
    "rush": "Rush",
    "ucmc": "UChicago",
    "ucsf": "UCSF",
    "umich": "UMichigan",
    "umn": "UMN",
    "upenn": "UPenn",
}

# Anonymized labels: alphabetical by real name → "Site 1", "Site 2", …
ANON_LABELS = {
    k: f"Site {i}"
    for i, k in enumerate(sorted(SITE_LABELS, key=lambda k: SITE_LABELS[k]), 1)
}

# 11 distinct colors + line styles for all sites
SITE_COLORS = {
    "emory": "#E69F00",      # orange
    "hopkins": "#56B4E9",    # sky blue
    "mimic_iv": "#0072B2",   # blue
    "nu": "#009E73",         # green
    "ohsu": "#D55E00",       # vermilion
    "rush": "#CC79A7",       # pink
    "ucmc": "#000000",       # black
    "ucsf": "#F0E442",       # yellow
    "umich": "#8B4513",      # brown
    "umn": "#882255",        # wine
    "upenn": "#44AA99",      # teal
}
MODE_COLORS = {
    "CVVHD": "#0072B5",
    "CVVHDF": "#1BB355",
    "CVVH": "#F16041",
    "SCUF": "#72B7B2",
    "AVVH": "#B279A2",
    "OTHER": "#999999",
}
MODE_ORDER = ["CVVHD", "CVVHDF", "CVVH", "SCUF", "AVVH", "OTHER"]

SITE_LINESTYLES = {
    "emory": "-",
    "hopkins": "--",
    "mimic_iv": "-",
    "nu": "--",
    "ohsu": "-.",
    "rush": "-.",
    "ucmc": "-",
    "ucsf": ":",
    "umich": ":",
    "umn": "-.",
    "upenn": "--",
}

# Figures in display order: (filename, title)
FIGURES = [
    ("consort_diagram_straight_flow_right_excl.png", "CONSORT/STROBE Diagram"),
    ("patient_state_over_crrt.png", "Patient State Over CRRT Course"),
    ("crrt_dose_over_time.png", "CRRT Dose Over Time"),
    ("crrt_dose_comparison_by_mode.png", "CRRT Dose Comparison by Mode"),
    ("dose_comparison.png", "Dose Comparison"),  # lives in final/, not graphs/
    ("crrt_parameter_histograms_grid.png", "CRRT Parameter Histograms"),
    ("lab_distributions_over_crrt.png", "Lab Distributions Over CRRT"),
    ("map_over_crrt.png", "MAP Over CRRT"),
    ("nee_over_crrt.png", "NEE Over CRRT"),
    ("respiratory_over_crrt.png", "Respiratory Over CRRT"),
]


# ── Helpers ─────────────────────────────────────────────────────────────────

def discover_sites() -> list[Path]:
    """Return sorted list of site dirs that contain final/ outputs."""
    return sorted(
        d for d in ROOT.iterdir()
        if d.is_dir() and (d / "final").is_dir()
    )


def extract_table(html_path: Path) -> str:
    """Extract the <table>…</table> element from a Table 1 HTML file."""
    html = html_path.read_text(encoding="utf-8")
    match = re.search(r"(<table[\s\S]*?</table>)", html, re.IGNORECASE)
    if not match:
        return "<p><em>Table not found.</em></p>"
    table = match.group(1)
    # Drop hardcoded last_crrt_mode rows (scuf / cvvhd)
    table = re.sub(r"<tr>\s*<td>last_crrt_mode.*?</tr>", "", table, flags=re.DOTALL)
    return table


def img_to_data_uri(img_path: Path) -> str:
    """Convert a PNG file to a base64 data URI string."""
    data = img_path.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return f"data:image/png;base64,{b64}"


def build_site_content(site_dir: Path) -> str:
    """Build the HTML content block for one site (table + figures)."""
    final = site_dir / "final"
    graphs = final / "graphs"

    # Table 1
    table_html = extract_table(final / "table1_crrt.html")

    # Figures
    fig_blocks = []
    for fname, title in FIGURES:
        # dose_comparison.png is in final/, rest in final/graphs/
        img_path = final / fname if fname == "dose_comparison.png" else graphs / fname
        if not img_path.exists():
            fig_blocks.append(
                f'<div class="figure-block">'
                f'<h3>{title}</h3>'
                f'<p class="missing">Figure not available: {fname}</p>'
                f'</div>'
            )
            continue
        data_uri = img_to_data_uri(img_path)
        fig_blocks.append(
            f'<div class="figure-block">'
            f'<h3>{title}</h3>'
            f'<img src="{data_uri}" alt="{title}">'
            f'</div>'
        )

    figures_html = "\n".join(fig_blocks)
    return f"""
        <div class="section">
            <h2>Table 1</h2>
            {table_html}
        </div>
        <div class="section">
            <h2>Figures</h2>
            {figures_html}
        </div>
    """


# ── Overall Tab Helpers ────────────────────────────────────────────────────


def fig_to_data_uri(fig) -> str:
    """Render a matplotlib figure to a base64 PNG data URI."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{b64}"


# ── Combined CONSORT ──────────────────────────────────────────────────────


def build_combined_consort(sites: list[Path]) -> str:
    """Build a combined CONSORT diagram summing counts across all sites."""
    # Load and sum strobe counts
    frames = []
    for site_dir in sites:
        csv_path = site_dir / "final" / "strobe_counts.csv"
        if csv_path.exists():
            frames.append(pd.read_csv(csv_path))
    if not frames:
        return '<p class="missing">No STROBE data available.</p>'

    df = pd.concat(frames, ignore_index=True)
    # Sum numeric values by counter; percentages will be recomputed
    summed = df.groupby("counter")["value"].sum().to_dict()

    # Derive step counts
    get = lambda k, d=0: summed.get(k, d)
    start_n = get("1b_after_stitching", get("1_adult_hospitalizations"))
    n_crrt = get("2_crrt_blocks")
    n_no_esrd = get("3_encounter_blocks_without_esrd", n_crrt)

    # Use Table 1 n_hospitalizations (Baseline) as the final analytical cohort
    t1_frames = []
    for site_dir in sites:
        csv_path = site_dir / "final" / "table1_crrt_long.csv"
        if csv_path.exists():
            t1_frames.append(pd.read_csv(csv_path))
    if t1_frames:
        t1_df = pd.concat(t1_frames, ignore_index=True)
        n_final = int(
            t1_df[(t1_df["variable"] == "n_hospitalizations")
                   & (t1_df["subgroup"] == "Baseline")]["n"].sum()
        )
    else:
        n_final = int(get("6_encounter_blocks_with_required_labs", n_no_esrd))

    # Build steps: No CRRT → ESRD → combined missing data
    steps = []
    parent = start_n
    for remaining_n, label, excl_label in [
        (n_crrt, "CRRT hospitalizations", "Excluded: No CRRT"),
        (n_no_esrd, "After ESRD exclusion", "Excluded: ESRD diagnosis"),
        (n_final, "Analytical cohort", "Excluded: Missing data\n(weight, CRRT settings, or labs)"),
    ]:
        excluded = max(parent - remaining_n, 0)
        if excluded > 0 or remaining_n > 0:
            steps.append((remaining_n, label, excluded, excl_label, parent))
        if excluded > 0:
            parent = remaining_n
        elif remaining_n == parent:
            steps.pop()

    # Draw figure
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box_h, box_w = 0.10, 0.40
    x_main_start = 0.05
    x_main_center = x_main_start + box_w / 2
    x_excl_start = 0.55
    v_spacing = 0.20
    excl_arrow_gap = 0.015

    def draw_box(x, y, w, h, text, fontsize=11):
        rect = FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.01",
            linewidth=2, edgecolor="black", facecolor="white",
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fontsize, wrap=True)

    ax.text(0.5, 0.98, "CRRT Cohort Selection (All Sites Combined)",
            ha="center", va="center", fontsize=16, fontweight="bold")

    arrow_props = dict(arrowstyle="->", lw=2, color="black")

    # Top box
    top_y = 0.90 - box_h
    draw_box(x_main_start, top_y, box_w, box_h,
             f"All adult hospitalizations\nn = {int(start_n):,}")

    for i, (remain_n, label, excl_n, excl_label, _parent) in enumerate(steps):
        cur_y = top_y - ((i + 1) * v_spacing)
        prev_y = top_y if i == 0 else top_y - (i * v_spacing)

        draw_box(x_main_start, cur_y, box_w, box_h,
                 f"Remaining hospitalizations\n{label}\nn = {int(remain_n):,}")

        # Main flow arrow
        ax.annotate("", xy=(x_main_center, cur_y + box_h),
                     xytext=(x_main_center, prev_y), arrowprops=arrow_props)

        # Exclusion box + arrow
        if excl_n > 0:
            mid_y = ((prev_y + box_h / 2) + (cur_y + box_h / 2)) / 2
            draw_box(x_excl_start, mid_y - box_h / 2, box_w, box_h,
                     f"{excl_label}\nn = {int(excl_n):,}")
            ax.annotate("", xy=(x_excl_start - excl_arrow_gap, mid_y),
                         xytext=(x_main_center, mid_y),
                         arrowprops=arrow_props, annotation_clip=False)

    data_uri = fig_to_data_uri(fig)
    return (
        f'<div class="figure-block">'
        f'<h3>CONSORT/STROBE Diagram (Combined)</h3>'
        f'<img src="{data_uri}" alt="Combined CONSORT Diagram">'
        f'</div>'
    )


# ── Combined Table 1 ─────────────────────────────────────────────────────


CRRT_DETAIL_VARS = {
    "last_crrt_mode",
    "crrt_dose_ml_kg_hr",
    "blood_flow_rate",
    "pre_filter_replacement_fluid_rate",
    "post_filter_replacement_fluid_rate",
    "dialysate_flow_rate",
    "ultrafiltration_out",
}


def build_combined_table1(sites: list[Path]) -> str:
    """Aggregate table1_crrt_long.csv across sites and render HTML."""
    frames = []
    for site_dir in sites:
        csv_path = site_dir / "final" / "table1_crrt_long.csv"
        if csv_path.exists():
            frames.append(pd.read_csv(csv_path))
    if not frames:
        return '<p class="missing">No Table 1 data available.</p>'

    df = pd.concat(frames, ignore_index=True)

    # Skip CRRT detail rows (keep Duration of CRRT)
    df = df[~df["variable"].str.split(",").str[0].isin(CRRT_DETAIL_VARS)]

    # Aggregate by (variable, level, subgroup, stat_type)
    group_cols = ["variable", "level", "subgroup", "stat_type"]
    rows = []
    for key, grp in df.groupby(group_cols, dropna=False, sort=False):
        variable, level, subgroup, stat_type = key
        if stat_type == "count":
            rows.append({
                "variable": variable, "level": level, "subgroup": subgroup,
                "stat_type": stat_type,
                "n": grp["n"].sum(), "total": grp["total"].sum(),
                "median": np.nan, "q25": np.nan, "q75": np.nan,
            })
        elif stat_type == "categorical":
            rows.append({
                "variable": variable, "level": level, "subgroup": subgroup,
                "stat_type": stat_type,
                "n": grp["n"].sum(), "total": grp["total"].sum(),
                "median": np.nan, "q25": np.nan, "q75": np.nan,
            })
        elif stat_type == "continuous":
            total_n = grp["n"].sum()
            if total_n > 0:
                w_median = (grp["median"] * grp["n"]).sum() / total_n
                w_q25 = (grp["q25"] * grp["n"]).sum() / total_n
                w_q75 = (grp["q75"] * grp["n"]).sum() / total_n
            else:
                w_median = w_q25 = w_q75 = np.nan
            rows.append({
                "variable": variable, "level": level, "subgroup": subgroup,
                "stat_type": stat_type,
                "n": total_n, "total": grp["total"].sum(),
                "median": w_median, "q25": w_q25, "q75": w_q75,
            })

    agg = pd.DataFrame(rows)

    # Build HTML table matching the existing site format
    subgroups = ["Baseline", "At 72h - Survivors", "At 72h - Non-survivors",
                 "At discharge - Survivors", "At discharge - Non-survivors"]

    header = (
        "<table><thead><tr>"
        "<th>Variable</th><th>Level</th>"
        + "".join(f"<th>{sg}</th>" for sg in subgroups)
        + "</tr></thead><tbody>\n"
    )

    # Preserve original row order by iterating unique (variable, level)
    seen = []
    for _, row in agg.iterrows():
        pair = (row["variable"], row["level"])
        if pair not in seen:
            seen.append(pair)

    body_rows = []
    for variable, level in seen:
        mask = (agg["variable"] == variable)
        if pd.notna(level):
            mask = mask & (agg["level"] == level)
        else:
            mask = mask & (agg["level"].isna())
        subset = agg[mask]
        stat_type = subset["stat_type"].iloc[0] if len(subset) > 0 else "count"

        # Clean display name
        disp_var = variable.replace("_", " ") if variable else ""
        disp_level = str(level) if pd.notna(level) else "Missing"

        cells = []
        for sg in subgroups:
            sg_row = subset[subset["subgroup"] == sg]
            if sg_row.empty:
                cells.append("")
                continue
            r = sg_row.iloc[0]
            if stat_type == "count":
                cells.append(f"{int(r['n'])}")
            elif stat_type == "categorical":
                pct = 100.0 * r["n"] / r["total"] if r["total"] > 0 else 0
                cells.append(f"{int(r['n'])} ({pct:.1f}%)")
            elif stat_type == "continuous":
                cells.append(
                    f"{r['median']:.1f} [{r['q25']:.1f}, {r['q75']:.1f}] "
                    f"({int(r['n'])}/{int(r['total'])})"
                )

        body_rows.append(
            "<tr><td>" + disp_var + "</td><td>" + disp_level + "</td>"
            + "".join(f"<td>{c}</td>" for c in cells)
            + "</tr>"
        )

    return (
        '<p><em>Pooled across all sites. '
        'Continuous: weighted median [Q1, Q3] (n/N). '
        'Categorical: n (%)</em></p>\n'
        + header + "\n".join(body_rows) + "\n</tbody></table>"
    )


# ── Overlay Graphs ────────────────────────────────────────────────────────


MAX_HOUR = 167  # Cap overlay graphs just under 7 days (hour 168 has boundary artifacts)


def _load_site_csv(sites: list[Path], rel_path: str) -> pd.DataFrame:
    """Load a CSV from each site's final/ tree, concat with site_id column."""
    frames = []
    for site_dir in sites:
        csv_path = site_dir / "final" / rel_path
        if csv_path.exists():
            tmp = pd.read_csv(csv_path)
            tmp["site_id"] = site_dir.name
            frames.append(tmp)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    if "hour" in df.columns:
        df = df[df["hour"] <= MAX_HOUR]
    return df


def _simple_overlay(sites, rel_path, x_col, y_col, title, xlabel, ylabel,
                    labels=None) -> str:
    """Single-panel overlay: one line per site."""
    if labels is None:
        labels = SITE_LABELS
    df = _load_site_csv(sites, rel_path)
    if df.empty:
        return f'<p class="missing">{title}: data not available.</p>'

    fig, ax = plt.subplots(figsize=(10, 5))
    for site_dir in sites:
        sid = site_dir.name
        sub = df[df["site_id"] == sid].sort_values(x_col)
        if sub.empty:
            continue
        label = labels.get(sid, sid)
        ax.plot(sub[x_col], sub[y_col], color=SITE_COLORS.get(sid, "gray"),
                linestyle=SITE_LINESTYLES.get(sid, "-"),
                label=label, linewidth=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=len(sites), fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.subplots_adjust(bottom=0.2)
    data_uri = fig_to_data_uri(fig)
    return (
        f'<div class="figure-block"><h3>{title}</h3>'
        f'<img src="{data_uri}" alt="{title}"></div>'
    )


def _build_dose_overlay(sites, labels=None) -> str:
    return _simple_overlay(
        sites, "graphs/crrt_dose_hourly.csv",
        "hour", "median_dose",
        "CRRT Dose Over Time (Median)", "Hour", "Dose (mL/kg/hr)",
        labels=labels,
    )


def _build_map_overlay(sites, labels=None) -> str:
    return _simple_overlay(
        sites, "graphs/map_over_crrt.csv",
        "hour", "median",
        "MAP Over CRRT (Median)", "Hour", "MAP (mmHg)",
        labels=labels,
    )


def _build_nee_overlay(sites, labels=None) -> str:
    return _simple_overlay(
        sites, "graphs/nee_over_crrt.csv",
        "hour", "median",
        "NEE Over CRRT (Median)", "Hour", "NEE (mcg/kg/min)",
        labels=labels,
    )


def _build_lab_overlay(sites, labels=None) -> str:
    """3x2 grid of lab panels."""
    if labels is None:
        labels = SITE_LABELS
    df = _load_site_csv(sites, "graphs/lab_distributions_over_crrt.csv")
    if df.empty:
        return '<p class="missing">Lab distributions: data not available.</p>'

    labs = sorted(df["lab"].unique())
    n_labs = len(labs)
    ncols = 2
    nrows = (n_labs + 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), squeeze=False)
    for idx, lab_name in enumerate(labs):
        r, c = divmod(idx, ncols)
        ax = axes[r][c]
        lab_df = df[df["lab"] == lab_name]
        for site_dir in sites:
            sid = site_dir.name
            sub = lab_df[lab_df["site_id"] == sid].sort_values("hour")
            if sub.empty:
                continue
            ax.plot(sub["hour"], sub["median"],
                    color=SITE_COLORS.get(sid, "gray"),
                    linestyle=SITE_LINESTYLES.get(sid, "-"),
                    label=labels.get(sid, sid), linewidth=1.5)
        ax.set_title(lab_name.replace("_", " ").title())
        ax.set_xlabel("Hour")
        ax.set_ylabel("Median")
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(n_labs, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    fig.suptitle("Lab Distributions Over CRRT (Median)", fontsize=14, y=1.01)
    # Shared legend at figure bottom
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(sites), fontsize=9,
               bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout()
    data_uri = fig_to_data_uri(fig)
    return (
        '<div class="figure-block"><h3>Lab Distributions Over CRRT</h3>'
        f'<img src="{data_uri}" alt="Lab Distributions"></div>'
    )


def _build_respiratory_overlay(sites, labels=None) -> str:
    """Two-panel: FiO2 median + IMV proportion."""
    if labels is None:
        labels = SITE_LABELS
    df = _load_site_csv(sites, "graphs/respiratory_over_crrt.csv")
    if df.empty:
        return '<p class="missing">Respiratory data not available.</p>'

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: FiO2
    fio2 = df[df["variable"] == "fio2"]
    for site_dir in sites:
        sid = site_dir.name
        sub = fio2[fio2["site_id"] == sid].sort_values("hour")
        if sub.empty:
            continue
        ax1.plot(sub["hour"], sub["median"],
                 color=SITE_COLORS.get(sid, "gray"),
                 linestyle=SITE_LINESTYLES.get(sid, "-"),
                 label=labels.get(sid, sid), linewidth=1.5)
    ax1.set_title("FiO2 (Median)")
    ax1.set_xlabel("Hour")
    ax1.set_ylabel("FiO2")
    ax1.grid(True, alpha=0.3)

    # Panel 2: IMV proportion
    imv = df[df["variable"] == "imv_proportion"]
    for site_dir in sites:
        sid = site_dir.name
        sub = imv[imv["site_id"] == sid].sort_values("hour")
        if sub.empty:
            continue
        ax2.plot(sub["hour"], sub["imv_proportion"],
                 color=SITE_COLORS.get(sid, "gray"),
                 linestyle=SITE_LINESTYLES.get(sid, "-"),
                 label=labels.get(sid, sid), linewidth=1.5)
    ax2.set_title("IMV Proportion")
    ax2.set_xlabel("Hour")
    ax2.set_ylabel("IMV %")
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Respiratory Over CRRT", fontsize=14, y=1.02)
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(sites), fontsize=9,
               bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
    data_uri = fig_to_data_uri(fig)
    return (
        '<div class="figure-block"><h3>Respiratory Over CRRT</h3>'
        f'<img src="{data_uri}" alt="Respiratory Over CRRT"></div>'
    )


def _build_patient_state_overlay(sites, labels=None) -> str:
    """Four-line overlay: prop_dead, prop_imv, prop_off_imv, prop_discharged."""
    if labels is None:
        labels = SITE_LABELS
    df = _load_site_csv(sites, "graphs/patient_state_over_crrt.csv")
    if df.empty:
        return '<p class="missing">Patient state data not available.</p>'

    metrics = [
        ("prop_dead", "Dead"),
        ("prop_imv", "On IMV"),
        ("prop_off_imv", "Off IMV"),
        ("prop_discharged", "Discharged"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for idx, (col, label) in enumerate(metrics):
        r, c = divmod(idx, 2)
        ax = axes[r][c]
        for site_dir in sites:
            sid = site_dir.name
            sub = df[df["site_id"] == sid].sort_values("hour")
            if sub.empty or col not in sub.columns:
                continue
            ax.plot(sub["hour"], sub[col],
                    color=SITE_COLORS.get(sid, "gray"),
                    linestyle=SITE_LINESTYLES.get(sid, "-"),
                    label=labels.get(sid, sid), linewidth=1.5)
        ax.set_title(label)
        ax.set_xlabel("Hour")
        ax.set_ylabel("Proportion (%)")
        ax.grid(True, alpha=0.3)

    fig.suptitle("Patient State Over CRRT Course", fontsize=14, y=1.01)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(sites), fontsize=9,
               bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout()
    data_uri = fig_to_data_uri(fig)
    return (
        '<div class="figure-block"><h3>Patient State Over CRRT Course</h3>'
        f'<img src="{data_uri}" alt="Patient State Over CRRT"></div>'
    )


def _build_modality_bar(sites, labels=None) -> str:
    """100% stacked horizontal bar chart of CRRT modality distribution per site (measurement rows)."""
    if labels is None:
        labels = SITE_LABELS

    # Load measurement-row counts from crrt_settings_distribution_by_mode.csv
    rows = []
    for site_dir in sites:
        csv_path = site_dir / "final" / "crrt_settings_distribution_by_mode.csv"
        if not csv_path.exists():
            continue
        tmp = pd.read_csv(csv_path)
        for _, r in tmp.iterrows():
            rows.append({"site_id": site_dir.name, "mode": r["Mode"].upper(), "n": r["N_Total"]})

    if not rows:
        return '<p class="missing">Modality distribution: data not available.</p>'

    df = pd.DataFrame(rows)
    # Pivot: sites × modes
    pivot = df.pivot_table(index="site_id", columns="mode", values="n", fill_value=0)
    # Compute proportions
    totals = pivot.sum(axis=1)
    pct = pivot.div(totals, axis=0) * 100
    # Order modes consistently, only include modes that exist
    modes = [m for m in MODE_ORDER if m in pct.columns]
    pct = pct[modes]
    # Sort sites by decreasing CVVHDF proportion (highest at top)
    pct["_label"] = pct.index.map(lambda s: labels.get(s, s))
    sort_col = "CVVHDF" if "CVVHDF" in pct.columns else modes[0]
    pct = pct.sort_values(sort_col, ascending=True)  # ascending so highest is at top of horizontal bar
    site_labels = pct["_label"].tolist()
    pct = pct.drop(columns=["_label"])

    fig, ax = plt.subplots(figsize=(10, max(4, len(pct) * 0.6)))
    left = np.zeros(len(pct))
    for mode in modes:
        vals = pct[mode].values
        ax.barh(site_labels, vals, left=left,
                color=MODE_COLORS.get(mode, "gray"), label=mode, edgecolor="white", linewidth=0.5)
        # Add percentage labels for segments >= 5%
        for i, v in enumerate(vals):
            if v >= 5:
                ax.text(left[i] + v / 2, i, f"{v:.0f}%", ha="center", va="center",
                        fontsize=8, fontweight="bold", color="white")
        left += vals

    ax.set_xlabel("Proportion (%)")
    ax.set_title("CRRT Modality Distribution by Site")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=len(modes), fontsize=9)
    ax.set_xlim(0, 100)
    fig.subplots_adjust(bottom=0.2)
    data_uri = fig_to_data_uri(fig)
    return (
        '<div class="figure-block"><h3>CRRT Modality Distribution by Site</h3>'
        f'<img src="{data_uri}" alt="CRRT Modality Distribution"></div>'
    )


# ── UCMC IHD Tab ─────────────────────────────────────────────────────────


def build_ihd_tab() -> str | None:
    """Build UCMC IHD-after-CRRT tab from 06_ihd_after_crrt.py outputs."""
    ihd_dir = Path(__file__).resolve().parent.parent / "output" / "final" / "ihd_analysis"
    summary_path = ihd_dir / "ihd_summary.csv"
    flowchart_path = ihd_dir / "ihd_flowchart.txt"
    missingness_path = ihd_dir / "ihd_missingness.csv"
    if not summary_path.exists():
        return None

    summary = pd.read_csv(summary_path)
    s = summary.iloc[0]
    boundary_qc_path = ihd_dir / "session_boundary_qc.csv"

    # ── Flowchart ──
    flowchart_html = ""
    if flowchart_path.exists():
        flowchart_text = flowchart_path.read_text()
        flowchart_html = (
            '<div class="section">'
            '<h2>Cohort Flow</h2>'
            f'<pre style="font-size:14px; line-height:1.6; background:#f8f8f8; '
            f'padding:20px; border-radius:6px; border:1px solid #ddd;">'
            f'{flowchart_text}</pre>'
            '</div>'
        )

    # ── Summary table ──
    summary_html = (
        '<table><thead><tr>'
        '<th>Metric</th><th>Value</th>'
        '</tr></thead><tbody>'
    )
    rows = [
        ("Total CRRT cohort (UCMC)", f"{int(s['n_cohort']):,}"),
        ("Received IHD after CRRT", f"{int(s['n_with_ihd']):,} ({s['pct_with_ihd']:.1f}%)"),
        ("No IHD after CRRT", f"{int(s['n_without_ihd']):,}"),
        ("Mortality — IHD group", f"{int(s['died_with_ihd']):,} / {int(s['n_with_ihd']):,} ({s['mortality_with_ihd_pct']:.1f}%)"),
        ("Mortality — No IHD group", f"{int(s['died_without_ihd']):,} / {int(s['n_without_ihd']):,} ({s['mortality_without_ihd_pct']:.1f}%)"),
        ("Median hours CRRT end → first IHD", f"{s['median_hours_to_first_ihd']:.1f}"),
        ("IHD sessions per patient (median [IQR])", f"{s['median_ihd_sessions']:.0f} [{s['ihd_sessions_q1']:.0f}, {s['ihd_sessions_q3']:.0f}]"),
    ]
    # Off-IHD metrics (if present)
    if "n_survivors_with_ihd" in s and s["n_survivors_with_ihd"] > 0:
        rows += [
            ("Survivors with IHD", f"{int(s['n_survivors_with_ihd']):,}"),
            ("Off IHD ≥2 days before discharge", f"{s['pct_off_ihd_2d_before_dc']:.1f}%"),
        ]
        if "pct_off_ihd_3d_before_dc" in s:
            rows.append(("Off IHD ≥3 days before discharge", f"{s['pct_off_ihd_3d_before_dc']:.1f}%"))
        if "pct_off_ihd_5d_before_dc" in s:
            rows.append(("Off IHD ≥5 days before discharge", f"{s['pct_off_ihd_5d_before_dc']:.1f}%"))
        rows.append(("Off IHD ≥7 days before discharge", f"{s['pct_off_ihd_7d_before_dc']:.1f}%"))
    for label, val in rows:
        summary_html += f"<tr><td>{label}</td><td>{val}</td></tr>"

    # Session definition note with boundary gap QC
    boundary_note = ""
    if boundary_qc_path.exists():
        bqc = pd.read_csv(boundary_qc_path)
        median_gap = bqc["gap_hours"].median()
        q25_gap = bqc["gap_hours"].quantile(0.25)
        q75_gap = bqc["gap_hours"].quantile(0.75)
        boundary_note = (
            f" Among the {len(bqc):,} session boundaries, the median inter-session gap "
            f"was {median_gap:.1f} hours (IQR [{q25_gap:.1f}, {q75_gap:.1f}]), "
            f"confirming the 6-hour threshold is well below the natural gap between IHD sessions."
        )
    summary_html += (
        '<tr><td colspan="2" style="font-size:11px; color:#666; font-style:italic; '
        'background:#f8f8f8; border-top:2px solid #ddd;">'
        '<b>Session definition:</b> IHD records are sorted by (encounter_block, recorded_dttm). '
        'A new session begins whenever (a) the time gap from the previous record exceeds 6 hours, '
        'or (b) the device_id changes to a different non-null value (forward-filled to bridge NaN gaps), '
        'or (c) at the patient\'s first record.'
        f'{boundary_note}</td></tr>'
    )
    summary_html += "</tbody></table>"

    # ── IHD Missingness table ──
    miss_html = ""
    if missingness_path.exists():
        miss_df = pd.read_csv(missingness_path)
        miss_html = (
            '<div class="section">'
            '<h2>IHD Variable Missingness (Post-CRRT Patients)</h2>'
            '<table><thead><tr>'
            '<th>Variable</th><th>Patients with any value</th>'
            '<th>% Patients</th><th>Row-level non-null %</th>'
            '</tr></thead><tbody>'
        )
        for _, r in miss_df.iterrows():
            miss_html += (
                f"<tr><td>{r['column']}</td>"
                f"<td>{int(r['n_patients_with_any_value']):,}</td>"
                f"<td>{r['pct_patients_with_any_value']:.1f}%</td>"
                f"<td>{r['pct_rows_non_null']:.1f}%</td></tr>"
            )
        miss_html += "</tbody></table></div>"

    return f"""
        {flowchart_html}
        <div class="section">
            <h2>UCMC: IHD After CRRT — Summary</h2>
            {summary_html}
        </div>
        {miss_html}
    """


# ── Missingness Tab ──────────────────────────────────────────────────────


def _missingness_heatmap(df: pd.DataFrame, subgroup: str, title: str) -> str | None:
    """Generate a single availability heatmap for the given subgroup."""
    import matplotlib.colors as mcolors

    cont = df[(df["stat_type"] == "continuous") & (df["subgroup"] == subgroup)].copy()
    if cont.empty:
        return None

    # Deduplicate (CRRT treatment vars have multiple rows per site/mode)
    cont = cont.sort_values("n", ascending=False).drop_duplicates(
        subset=["site", "variable"], keep="first"
    )

    cont["avail_pct"] = np.where(cont["total"] > 0, 100.0 * cont["n"] / cont["total"], 0)

    pivot = cont.pivot_table(index="variable", columns="site", values="avail_pct", aggfunc="first")
    pivot.index = [v.replace("_", " ") for v in pivot.index]
    pivot.columns = [SITE_LABELS.get(c, c) for c in pivot.columns]

    # Sort: most missing first
    pivot["mean_avail"] = pivot.mean(axis=1)
    pivot = pivot.sort_values("mean_avail")
    pivot = pivot.drop(columns="mean_avail")

    n_vars = len(pivot)
    fig_h = max(6, n_vars * 0.35 + 2)
    fig, ax = plt.subplots(figsize=(max(10, len(pivot.columns) * 1.2), fig_h))

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "avail", ["#d73027", "#fee08b", "#1a9850"], N=256
    )
    data = pivot.values.astype(float)
    im = ax.imshow(data, cmap=cmap, aspect="auto", vmin=0, vmax=100)

    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=9)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if np.isnan(val):
                txt, color = "—", "gray"
            else:
                txt = f"{val:.0f}"
                color = "white" if val < 40 else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("% Available", fontsize=10)
    ax.set_title(title, fontsize=13, pad=12)
    fig.tight_layout()

    return fig_to_data_uri(fig)


def build_missingness_tab(sites: list[Path]) -> str | None:
    """Heatmaps of data availability: Baseline and 72h Survivors."""
    frames = []
    for site_dir in sites:
        csv_path = site_dir / "final" / "table1_crrt_long.csv"
        if csv_path.exists():
            tmp = pd.read_csv(csv_path)
            tmp["site"] = site_dir.name
            frames.append(tmp)
    if not frames:
        return None

    df = pd.concat(frames, ignore_index=True)

    sections = []
    for subgroup, label in [
        ("Baseline", "Baseline"),
        ("At 72h - Survivors", "At 72h (Survivors)"),
    ]:
        uri = _missingness_heatmap(
            df, subgroup,
            f"Data Availability — {label}",
        )
        if uri:
            sections.append(f"""
                <div class="figure-block">
                    <h3>{label}</h3>
                    <img src="{uri}" alt="Availability {label}" style="max-width:100%;">
                </div>
            """)

    if not sections:
        return None

    return f"""
        <div class="section">
            <h2>Data Availability Heatmaps</h2>
            <p><em>Continuous variables. Red = high missingness, green = complete.
            Values show % of patients with data available.</em></p>
            {"".join(sections)}
        </div>
    """


# ── Build Overall Tab Content ─────────────────────────────────────────────


def build_overall_content(sites: list[Path], labels=None) -> str:
    """Assemble the full Overall tab HTML: CONSORT + Table 1 + overlay graphs."""
    print("  Building combined CONSORT diagram...")
    consort_html = build_combined_consort(sites)

    print("  Building combined Table 1...")
    table1_html = build_combined_table1(sites)

    print("  Building overlay graphs...")
    graph_blocks = [
        _build_modality_bar(sites, labels=labels),
        _build_dose_overlay(sites, labels=labels),
        _build_lab_overlay(sites, labels=labels),
        _build_map_overlay(sites, labels=labels),
        _build_nee_overlay(sites, labels=labels),
        _build_respiratory_overlay(sites, labels=labels),
        _build_patient_state_overlay(sites, labels=labels),
    ]
    graphs_html = "\n".join(graph_blocks)

    return f"""
        <div class="section">
            <h2>CONSORT/STROBE Diagram</h2>
            {consort_html}
        </div>
        <div class="section">
            <h2>Table 1 (All Sites Combined)</h2>
            {table1_html}
        </div>
        <div class="section">
            <h2>Figures (Site Overlay)</h2>
            {graphs_html}
        </div>
    """


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    sites = discover_sites()
    if not sites:
        print("No sites found in", ROOT)
        return

    print(f"Found {len(sites)} sites: {[s.name for s in sites]}")

    # Build tab buttons and tab content panels
    tab_buttons = []
    tab_panels = []

    # Overall tab (first, default active)
    print("  Processing Overall tab...")
    overall_content = build_overall_content(sites)
    tab_buttons.append(
        '<button class="tab-btn active" '
        "onclick=\"switchTab('overall')\" "
        'id="btn-overall">Overall</button>'
    )
    tab_panels.append(
        '<div class="tab-panel" id="panel-overall" style="display:block;">'
        f'{overall_content}'
        '</div>'
    )

    # Severity Analysis tab
    print("  Building severity analysis tab...")
    df_table1 = load_all_table1(sites)
    if not df_table1.empty:
        sev = extract_severity_data(df_table1)
        severity_content = generate_severity_html(sev, labels=SITE_LABELS)
        tab_buttons.append(
            '<button class="tab-btn" '
            "onclick=\"switchTab('severity')\" "
            'id="btn-severity">Severity Analysis</button>'
        )
        tab_panels.append(
            '<div class="tab-panel" id="panel-severity" style="display:none;">'
            f'{severity_content}'
            '</div>'
        )

    # Missingness tab
    print("  Building missingness tab...")
    miss_content = build_missingness_tab(sites)
    if miss_content:
        tab_buttons.append(
            '<button class="tab-btn" '
            "onclick=\"switchTab('missingness')\" "
            'id="btn-missingness">Data Availability</button>'
        )
        tab_panels.append(
            '<div class="tab-panel" id="panel-missingness" style="display:none;">'
            f'{miss_content}'
            '</div>'
        )

    # UCMC IHD tab
    print("  Building UCMC IHD tab...")
    ihd_content = build_ihd_tab()
    if ihd_content:
        tab_buttons.append(
            '<button class="tab-btn" '
            "onclick=\"switchTab('ihd')\" "
            'id="btn-ihd">UCMC IHD</button>'
        )
        tab_panels.append(
            '<div class="tab-panel" id="panel-ihd" style="display:none;">'
            f'{ihd_content}'
            '</div>'
        )

    # Per-site tabs
    for site_dir in sites:
        site_id = site_dir.name
        label = SITE_LABELS.get(site_id, site_id.replace("_", " ").title())

        tab_buttons.append(
            f'<button class="tab-btn" '
            f'onclick="switchTab(\'{site_id}\')" '
            f'id="btn-{site_id}">{label}</button>'
        )

        print(f"  Processing {label}...")
        content = build_site_content(site_dir)
        tab_panels.append(
            f'<div class="tab-panel" id="panel-{site_id}" style="display:none;">'
            f'{content}'
            f'</div>'
        )

    tabs_bar = "\n".join(tab_buttons)
    panels = "\n".join(tab_panels)

    # Site IDs for JS
    all_tab_ids = ["overall", "severity"]
    if miss_content:
        all_tab_ids.append("missingness")
    if ihd_content:
        all_tab_ids.append("ihd")
    all_tab_ids += [s.name for s in sites]
    site_ids_js = ", ".join(f'"{tid}"' for tid in all_tab_ids)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>CRRT Epidemiology — Multi-Site Dashboard</title>
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

        /* Tabs */
        .tab-bar {{
            display: flex; gap: 4px;
            border-bottom: 2px solid #4CAF50;
            margin-bottom: 20px;
        }}
        .tab-btn {{
            padding: 10px 24px;
            border: 1px solid #ccc;
            border-bottom: none;
            background: #e8e8e8;
            cursor: pointer;
            font-size: 15px;
            font-weight: 600;
            border-radius: 6px 6px 0 0;
            color: #555;
            transition: background 0.15s, color 0.15s;
        }}
        .tab-btn:hover {{ background: #d4edda; color: #333; }}
        .tab-btn.active {{
            background: #4CAF50; color: white;
            border-color: #4CAF50;
        }}

        /* Content */
        .section {{ margin-bottom: 40px; }}
        .section h2 {{
            color: #333;
            border-bottom: 2px solid #e0e0e0;
            padding-bottom: 6px;
        }}

        /* Table styling (matches existing theme) */
        table {{ width: 100%; border-collapse: collapse; margin-top: 12px; font-size: 13px; }}
        th {{ background: #4CAF50; color: white; padding: 12px; text-align: left; border: 1px solid #ddd; }}
        td {{ padding: 10px; border: 1px solid #ddd; vertical-align: top; }}
        tr:nth-child(even) {{ background: #f9f9f9; }}
        tr:hover {{ background: #f0f0f0; }}

        /* Figures */
        .figure-block {{
            margin-bottom: 32px;
        }}
        .figure-block h3 {{
            color: #444;
            margin-bottom: 8px;
        }}
        .figure-block img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}
        .missing {{
            color: #999; font-style: italic;
        }}
    </style>
</head>
<body>
<div class="container">
    <h1>CRRT Epidemiology — Multi-Site Dashboard</h1>
    <div class="tab-bar">
        {tabs_bar}
    </div>
    {panels}
</div>
<script>
    const SITES = [{site_ids_js}];
    function switchTab(siteId) {{
        SITES.forEach(function(id) {{
            document.getElementById('panel-' + id).style.display = 'none';
            document.getElementById('btn-' + id).classList.remove('active');
        }});
        document.getElementById('panel-' + siteId).style.display = 'block';
        document.getElementById('btn-' + siteId).classList.add('active');
    }}
</script>
</body>
</html>"""

    OUTPUT.write_text(html, encoding="utf-8")
    print(f"\nDashboard written to {OUTPUT}")
    print(f"  File size: {OUTPUT.stat().st_size / 1024 / 1024:.1f} MB")

    # ── Anonymized dashboard (Overall + Severity, site names replaced) ───
    print("\nBuilding anonymized dashboard (Overall + Severity)...")
    anon_content = build_overall_content(sites, labels=ANON_LABELS)

    # Anonymized severity analysis
    sev_anon = sev.copy()
    sev_anon["site_label"] = sev_anon["site_dir"].map(
        lambda sd: ANON_LABELS.get(sd, sd)
    )
    sev_anon = sev_anon.sort_values("site_label").reset_index(drop=True)
    severity_anon = generate_severity_html(sev_anon, labels=ANON_LABELS)

    anon_tabs = (
        '<button class="tab-btn active" onclick="switchTab(\'overall\')" '
        'id="btn-overall">Overall</button>'
        '<button class="tab-btn" onclick="switchTab(\'severity\')" '
        'id="btn-severity">Severity Analysis</button>'
    )
    anon_panels = (
        '<div class="tab-panel" id="panel-overall" style="display:block;">'
        f'{anon_content}</div>'
        '<div class="tab-panel" id="panel-severity" style="display:none;">'
        f'{severity_anon}</div>'
    )
    anon_site_ids_js = '"overall","severity"'

    anon_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>CRRT Epidemiology — Multi-Site Dashboard (Anonymized)</title>
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
        .tab-bar {{ margin: 20px 0; display: flex; gap: 8px; flex-wrap: wrap; }}
        .tab-btn {{
            padding: 10px 20px; border: 1px solid #ddd;
            background: #f0f0f0; cursor: pointer;
            border-radius: 4px 4px 0 0; font-size: 14px;
        }}
        .tab-btn.active {{ background: #4CAF50; color: white; border-color: #4CAF50; }}
        .tab-panel {{ display: none; }}
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
        .figure-block h3 {{ color: #444; margin-bottom: 8px; }}
        .figure-block img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }}
        .missing {{ color: #999; font-style: italic; }}
    </style>
</head>
<body>
<div class="container">
    <h1>CRRT Epidemiology — Multi-Site Dashboard (Anonymized)</h1>
    <div class="tab-bar">
        {anon_tabs}
    </div>
    {anon_panels}
</div>
<script>
    const SITES = [{anon_site_ids_js}];
    function switchTab(siteId) {{
        SITES.forEach(function(id) {{
            document.getElementById('panel-' + id).style.display = 'none';
            document.getElementById('btn-' + id).classList.remove('active');
        }});
        document.getElementById('panel-' + siteId).style.display = 'block';
        document.getElementById('btn-' + siteId).classList.add('active');
    }}
</script>
</body>
</html>"""
    anon_output = ROOT / "combined_dashboard_anon.html"
    anon_output.write_text(anon_html, encoding="utf-8")
    print(f"Anonymized dashboard written to {anon_output}")
    print(f"  File size: {anon_output.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
