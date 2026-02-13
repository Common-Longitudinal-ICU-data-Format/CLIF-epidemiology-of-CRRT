"""
04_combine_site_results.py
Combines per-site Table 1 HTML and figures into a single tabbed dashboard.
Output: all_site_data/combined_dashboard.html
"""

import base64
import io
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent.parent / "all_site_data"
OUTPUT = ROOT / "combined_dashboard.html"

# Display names for sites (fallback: folder name title-cased)
SITE_LABELS = {
    "mimic_iv": "MIMIC-IV",
    "nu": "Northwestern",
    "rush": "Rush",
    "ucmc": "UChicago",
}

# Colors for overlay graphs (keyed by site folder name)
SITE_COLORS = {
    "mimic_iv": "#1f77b4",   # blue
    "nu": "#ff7f0e",         # orange
    "rush": "#2ca02c",       # green
    "ucmc": "#d62728",       # red
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

    # Derive step counts (same logic as create_consort_diagram_straight_flow)
    get = lambda k, d=0: summed.get(k, d)
    start_n = get("1b_after_stitching", get("1_adult_hospitalizations"))
    n_crrt = get("2_crrt_blocks")
    n_no_esrd = get("3_encounter_blocks_without_esrd", n_crrt)
    n_with_weight = get("4_encounter_blocks_with_weight", n_no_esrd)
    n_with_settings = get("5_encounter_blocks_with_crrt_settings", n_with_weight)
    n_with_labs = get("6_encounter_blocks_with_required_labs", n_with_settings)

    # Build steps
    steps = []
    parent = start_n
    for remaining_n, label, excl_label in [
        (n_crrt, "CRRT hospitalizations", "Excluded: No CRRT"),
        (n_no_esrd, "After ESRD exclusion", "Excluded: ESRD diagnosis"),
        (n_with_weight, "With documented weight", "Excluded: Missing weight"),
        (n_with_settings, "With CRRT settings", "Excluded: Missing CRRT settings"),
        (n_with_labs, "With required labs", "Excluded: Missing required labs"),
    ]:
        excluded = max(parent - remaining_n, 0)
        if excluded > 0 or remaining_n > 0:
            steps.append((remaining_n, label, excluded, excl_label, parent))
        if excluded > 0:
            parent = remaining_n
        # If excluded == 0, skip the step but keep parent unchanged
        elif remaining_n == parent:
            steps.pop()  # remove the just-added zero-exclusion step

    # Draw figure
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box_h, box_w = 0.08, 0.40
    x_main_start = 0.05
    x_main_center = x_main_start + box_w / 2
    x_excl_start = 0.55
    v_spacing = 0.14
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

    # Analytical notes
    n_aki = get("6_encounter_blocks_with_AKI_no_esrd", 0)
    n_hosp_no_icu = get("6_number_hosp_without_ICU_stay", 0)
    if n_with_labs > 0:
        pct_aki = 100.0 * n_aki / n_with_labs
        n_with_icu = int(n_with_labs - n_hosp_no_icu)
        pct_icu = 100.0 * n_with_icu / n_with_labs
        note = (
            f"AKI codes present (non-ESRD): {int(n_aki):,} / {int(n_with_labs):,} ({pct_aki:.1f}%)\n"
            f"CRRT hospitalizations with ICU admission: {n_with_icu:,} / {int(n_with_labs):,} ({pct_icu:.1f}%)"
        )
        note_rect = FancyBboxPatch(
            (0.05, 0.005), 0.89, 0.085, boxstyle="round,pad=0.01",
            linewidth=1, edgecolor="#bbbbbb", facecolor="#f6f6ee", zorder=0,
        )
        ax.add_patch(note_rect)
        ax.text(0.5, 0.05, note, ha="center", va="center", fontsize=11)

    data_uri = fig_to_data_uri(fig)
    return (
        f'<div class="figure-block">'
        f'<h3>CONSORT/STROBE Diagram (Combined)</h3>'
        f'<img src="{data_uri}" alt="Combined CONSORT Diagram">'
        f'</div>'
    )


# ── Combined Table 1 ─────────────────────────────────────────────────────


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

    # Aggregate by (variable, level, subgroup, stat_type)
    group_cols = ["variable", "level", "subgroup", "stat_type"]
    rows = []
    for key, grp in df.groupby(group_cols, dropna=False):
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
        disp_level = str(level) if pd.notna(level) else ""

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


def _simple_overlay(sites, rel_path, x_col, y_col, title, xlabel, ylabel) -> str:
    """Single-panel overlay: one line per site."""
    df = _load_site_csv(sites, rel_path)
    if df.empty:
        return f'<p class="missing">{title}: data not available.</p>'

    fig, ax = plt.subplots(figsize=(10, 5))
    for site_dir in sites:
        sid = site_dir.name
        sub = df[df["site_id"] == sid].sort_values(x_col)
        if sub.empty:
            continue
        label = SITE_LABELS.get(sid, sid)
        ax.plot(sub[x_col], sub[y_col], color=SITE_COLORS.get(sid, "gray"),
                label=label, linewidth=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    data_uri = fig_to_data_uri(fig)
    return (
        f'<div class="figure-block"><h3>{title}</h3>'
        f'<img src="{data_uri}" alt="{title}"></div>'
    )


def _build_dose_overlay(sites) -> str:
    return _simple_overlay(
        sites, "graphs/crrt_dose_hourly.csv",
        "hour", "median_dose",
        "CRRT Dose Over Time (Median)", "Hour", "Dose (mL/kg/hr)",
    )


def _build_map_overlay(sites) -> str:
    return _simple_overlay(
        sites, "graphs/map_over_crrt.csv",
        "hour", "median",
        "MAP Over CRRT (Median)", "Hour", "MAP (mmHg)",
    )


def _build_nee_overlay(sites) -> str:
    return _simple_overlay(
        sites, "graphs/nee_over_crrt.csv",
        "hour", "median",
        "NEE Over CRRT (Median)", "Hour", "NEE (mcg/kg/min)",
    )


def _build_lab_overlay(sites) -> str:
    """3x2 grid of lab panels."""
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
                    label=SITE_LABELS.get(sid, sid), linewidth=1.5)
        ax.set_title(lab_name.replace("_", " ").title())
        ax.set_xlabel("Hour")
        ax.set_ylabel("Median")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(n_labs, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    fig.suptitle("Lab Distributions Over CRRT (Median)", fontsize=14, y=1.01)
    fig.tight_layout()
    data_uri = fig_to_data_uri(fig)
    return (
        '<div class="figure-block"><h3>Lab Distributions Over CRRT</h3>'
        f'<img src="{data_uri}" alt="Lab Distributions"></div>'
    )


def _build_respiratory_overlay(sites) -> str:
    """Two-panel: FiO2 median + IMV proportion."""
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
                 label=SITE_LABELS.get(sid, sid), linewidth=1.5)
    ax1.set_title("FiO2 (Median)")
    ax1.set_xlabel("Hour")
    ax1.set_ylabel("FiO2")
    ax1.legend(fontsize=8)
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
                 label=SITE_LABELS.get(sid, sid), linewidth=1.5)
    ax2.set_title("IMV Proportion")
    ax2.set_xlabel("Hour")
    ax2.set_ylabel("IMV %")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Respiratory Over CRRT", fontsize=14, y=1.02)
    fig.tight_layout()
    data_uri = fig_to_data_uri(fig)
    return (
        '<div class="figure-block"><h3>Respiratory Over CRRT</h3>'
        f'<img src="{data_uri}" alt="Respiratory Over CRRT"></div>'
    )


def _build_patient_state_overlay(sites) -> str:
    """Four-line overlay: prop_dead, prop_imv, prop_off_imv, prop_discharged."""
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
                    label=SITE_LABELS.get(sid, sid), linewidth=1.5)
        ax.set_title(label)
        ax.set_xlabel("Hour")
        ax.set_ylabel("Proportion (%)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Patient State Over CRRT Course", fontsize=14, y=1.01)
    fig.tight_layout()
    data_uri = fig_to_data_uri(fig)
    return (
        '<div class="figure-block"><h3>Patient State Over CRRT Course</h3>'
        f'<img src="{data_uri}" alt="Patient State Over CRRT"></div>'
    )


# ── Build Overall Tab Content ─────────────────────────────────────────────


def build_overall_content(sites: list[Path]) -> str:
    """Assemble the full Overall tab HTML: CONSORT + Table 1 + overlay graphs."""
    print("  Building combined CONSORT diagram...")
    consort_html = build_combined_consort(sites)

    print("  Building overlay graphs...")
    graph_blocks = [
        _build_dose_overlay(sites),
        _build_lab_overlay(sites),
        _build_map_overlay(sites),
        _build_nee_overlay(sites),
        _build_respiratory_overlay(sites),
        _build_patient_state_overlay(sites),
    ]
    graphs_html = "\n".join(graph_blocks)

    return f"""
        <div class="section">
            <h2>CONSORT/STROBE Diagram</h2>
            {consort_html}
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

    # Site IDs for JS (include "overall")
    all_tab_ids = ["overall"] + [s.name for s in sites]
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


if __name__ == "__main__":
    main()
