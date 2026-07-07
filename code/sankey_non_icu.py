"""
Sankey of ADT location trajectories for the CRRT encounter_blocks flagged as
"no ICU stay" (00_cohort.py). Reads the pipeline's saved
`{output}/intermediate/adt_df_non_icu_hosps.csv` — aggregate location flows only.

Usage:  uv run python code/sankey_non_icu.py output_nu    # or output_ucmc
"""
import collections
import sys
from pathlib import Path

import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")  # headless for CLI; when imported (00_cohort.py) the caller owns the backend
import matplotlib.path
import matplotlib.patches
import matplotlib.patheffects
import matplotlib.pyplot as plt
import pandas as pd

# Categorical node colors — CVD-safe slots from the dataviz palette (fixed order).
COLORS = {
    "ED":         "#2a78d6",  # blue
    "Procedural": "#1baf7a",  # aqua
    "Stepdown":   "#eda100",  # yellow
    "Ward":       "#008300",  # green
    "ICU":        "#e34948",  # red — flags the anomalous ICU records
    "Other":      "#9aa0a6",  # neutral
    "Discharged": "#cbccc0",  # neutral terminal
}
ORDER_ALL = ["ED", "Procedural", "Stepdown", "Ward", "ICU", "Other", "Discharged"]


def simplify_location_category(loc):
    if pd.isna(loc):
        return "Other"
    l = str(loc).strip().lower()
    return {"icu": "ICU", "ward": "Ward", "ed": "ED",
            "procedural": "Procedural", "stepdown": "Stepdown"}.get(l, "Other")


def prepare_wide(adt, max_locations=6):
    a = adt.copy()
    a["in_dttm"] = pd.to_datetime(a["in_dttm"], errors="coerce", utc=True)
    a = a.sort_values(["encounter_block", "in_dttm"])
    a["loc"] = a["location_category"].apply(simplify_location_category)
    # Drop reverse transitions back to ED (stitching artifacts).
    a["prev"] = a.groupby("encounter_block")["loc"].shift(1)
    rev = (a["loc"] == "ED") & (a["prev"].isin(["ICU", "Ward", "Procedural", "Stepdown"]))
    if rev.sum():
        print(f"  dropping {int(rev.sum())} reverse->ED transitions (stitching artifacts)")
    a = a[~rev].drop(columns=["prev"])
    a["rank"] = a.groupby("encounter_block").cumcount() + 1
    a = a[a["rank"] <= max_locations]
    wide = a.pivot_table(index="encounter_block", columns="rank",
                         values="loc", aggfunc="first")
    cols = list(wide.columns)
    wide = wide.reset_index()
    for c in cols:
        wide[c] = wide[c].fillna("Discharged")
    return wide, cols


def admission_hist(adt, site, gdir):
    """Histogram of admission timing (earliest ADT in_dttm per encounter_block)
    for the no-ICU CRRT encounters, with the COVID window shaded — to see whether
    these are concentrated in 2020–2021."""
    a = adt.copy()
    a["in_dttm"] = pd.to_datetime(a["in_dttm"], errors="coerce", utc=True)
    adm = a.dropna(subset=["in_dttm"]).groupby("encounter_block")["in_dttm"].min()
    adm = adm.dt.tz_convert(None)  # tz-naive for plotting
    n = len(adm)
    if n == 0:
        print("  no admission times available — skipping histogram")
        return

    lo = adm.min().to_period("M").to_timestamp()
    hi = adm.max().to_period("M").to_timestamp() + pd.offsets.MonthBegin(1)
    bins = pd.date_range(lo, hi, freq="MS")

    fig, ax = plt.subplots(figsize=(10, 4.2))
    covid_lo, covid_hi = pd.Timestamp("2020-03-01"), pd.Timestamp("2021-12-31")
    ax.axvspan(covid_lo, covid_hi, color="#eda100", alpha=0.16, label="COVID (Mar 2020–2021)")
    ax.hist(adm, bins=bins, color="#2a78d6", edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Admission (earliest ADT in_dttm, monthly bins)")
    ax.set_ylabel("no-ICU CRRT encounters")
    ax.set_title(f"{site}: admission timing of the {n} CRRT encounters with no ICU stay",
                 loc="left", fontsize=12)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    out = gdir / f"{site}_non_icu_admission_hist.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    yr = adm.dt.year.value_counts().sort_index()
    covid_n = int(((adm >= covid_lo) & (adm <= covid_hi)).sum())
    print(f"  admission by year:  " + "  ".join(f"{y}:{c}" for y, c in yr.items()))
    print(f"  COVID (Mar 2020–Dec 2021): {covid_n}/{n} ({covid_n/n*100:.0f}%)  -> {out}")


# ---- the provided Sankey class (get_distinct_colors dependency removed; colors required) ----
class Sankey:
    def __init__(self, df, plot_width=9, plot_height=7, gap=0.12, alpha=0.35,
                 block_width=0.3, block_fontsize=11, order=None, colors=None,
                 flow_color_func=None, ax=None):
        self.df = df
        self.plot_width, self.plot_height = plot_width, plot_height
        self.gap, self.block_width, self.block_fontsize = gap, block_width, block_fontsize
        self.alpha, self.order, self.flow_color_func = alpha, order, flow_color_func
        self.init_figure(ax)
        self.init_flows()
        self.init_nodes(order)
        self.init_widths()
        self.resolution = (plot_height - gap * (len(order) - 1)) / df.shape[0]
        self.colors = colors
        self.init_offsets()

    def init_figure(self, ax):
        if ax is None:
            self.fig = plt.figure(figsize=(self.plot_width, self.plot_height))
            self.ax = plt.Axes(self.fig, [0, 0, 1, 1]); self.fig.add_axes(self.ax)
        else:
            self.fig, self.ax = ax.figure, ax

    def init_flows(self):
        self.flows = []
        for i in range(self.df.columns.size - 1):
            x, y = self.df.iloc[:, i], self.df.iloc[:, i + 1]
            self.flows.append(collections.Counter(zip(x, y)))

    def init_nodes(self, order):
        self.nodes = []
        for i in range(self.df.columns.size):
            col = collections.OrderedDict()
            counts = self.df.iloc[:, i].value_counts()
            for item in order:
                col[item] = int(counts[item]) if item in counts else 0
            self.nodes.append(col)

    def init_widths(self):
        self.left_stop = self.block_width
        self.right_stop = self.plot_width - self.block_width
        self.stops = []
        n = self.df.columns.size
        self.flow_width = (self.plot_width - self.block_width * (n - 2)) / (n - 1)
        for i in range(1, n):
            s1 = self.block_width * i + self.flow_width * (i - 1) + self.flow_width * 7 / 20
            s2 = self.block_width * i + self.flow_width * (i - 1) + self.flow_width * 13 / 20
            self.stops.append((s1, s2))

    def init_offsets(self):
        self.offsets = []
        for col in self.nodes:
            off, offs = 0, collections.OrderedDict()
            for name, size in col.items():
                offs[name] = off
                off += size * self.resolution + self.gap
            self.offsets.append(offs)

    def draw_flow(self, x, left, right, flow, nol, nor):
        P = matplotlib.path.Path
        ly = self.offsets[x][left] + nol[left]
        ry = self.offsets[x + 1][right] + nor[right]
        flow *= self.resolution
        nol[left] += flow; nor[right] += flow
        color = self.colors[left]
        lx = self.flow_width * x + self.block_width * (x + 1)
        rx = lx + self.flow_width
        pd_ = [(P.MOVETO, (lx, -ly)), (P.LINETO, (lx, -ly - flow)),
               (P.CURVE4, (self.stops[x][0], -ly - flow)), (P.CURVE4, (self.stops[x][1], -ry - flow)),
               (P.CURVE4, (rx, -ry - flow)), (P.LINETO, (rx, -ry)),
               (P.CURVE4, (self.stops[x][1], -ry)), (P.CURVE4, (self.stops[x][0], -ly)),
               (P.CURVE4, (lx, -ly)), (P.CLOSEPOLY, (lx, -ly))]
        codes, verts = zip(*pd_)
        self.ax.add_patch(matplotlib.patches.PathPatch(
            P(verts, codes), facecolor=color, alpha=self.alpha, edgecolor="none"))

    def draw_node(self, x, y, size, name):
        if size <= 0:
            return
        yy = -list(self.offsets[x].values())[y] - size * self.resolution
        xx = self.flow_width * x + self.block_width * x
        self.ax.add_patch(matplotlib.patches.Rectangle(
            (xx, yy), self.block_width, size * self.resolution,
            facecolor=self.colors[name], edgecolor="none"))
        self.ax.text(xx + self.block_width / 2, yy + size * self.resolution / 2,
                     f"{name}\n{size}", color="black", va="center", ha="center",
                     size=self.block_fontsize,
                     path_effects=[matplotlib.patheffects.Stroke(linewidth=2, foreground="white"),
                                   matplotlib.patheffects.Normal()])

    def draw(self):
        for x, col in enumerate(self.nodes):
            for y, (name, size) in enumerate(col.items()):
                self.draw_node(x, y, size, name)
        for x, flows in enumerate(self.flows):
            nol, nor = collections.Counter(), collections.Counter()
            for (l, r), fl in sorted(flows.items(),
                                     key=lambda t: (self.order.index(t[0][0]), self.order.index(t[0][1]))):
                self.draw_flow(x, l, r, fl, nol, nor)
        self.ax.set_ylim(-self.resolution * self.df.shape[0] - self.gap * (len(self.order) - 1), 0)
        self.ax.set_xlim(0, self.block_width * self.df.shape[1] + self.flow_width * (self.df.shape[1] - 1))
        self.ax.get_xaxis().set_visible(False); self.ax.get_yaxis().set_visible(False)
        for k in self.ax.spines: self.ax.spines[k].set_visible(False)


def generate(adt, site, gdir):
    """Build the no-ICU CRRT location Sankey + admission-timing histogram from an
    in-memory ADT frame (needs columns: encounter_block, in_dttm, location_category)
    and save both PNGs to gdir. Called by the pipeline (00_cohort.py) and the CLI."""
    gdir = Path(gdir)
    gdir.mkdir(parents=True, exist_ok=True)
    wide, cols = prepare_wide(adt)
    n = len(wide)
    if n < 3:
        print(f"  [{site}] only {n} no-ICU encounter(s) — skipping Sankey")
    else:
        order = [c for c in ORDER_ALL if c in set(pd.unique(wide[cols].values.ravel()))]
        sk = Sankey(wide[cols], order=order, colors=COLORS, plot_height=max(6, n / 12))
        sk.draw()
        sk.ax.set_title(f"{site}: location trajectory of the {n} CRRT encounters with NO ICU stay",
                        fontsize=12, pad=10, loc="left")
        outpng = gdir / f"{site}_non_icu_sankey.png"
        sk.fig.savefig(outpng, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(sk.fig)
        print(f"  [{site}] non-ICU Sankey ({n} encounters) -> {outpng}")
    admission_hist(adt, site, gdir)


def main():
    OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("output")
    if not OUT.is_absolute():
        OUT = Path(__file__).resolve().parent.parent / OUT
    logs = sorted(OUT.glob("final/*_pipeline_*.log"))
    site = logs[-1].name.split("_pipeline_")[0] if logs else OUT.name.replace("output_", "").upper()
    f = OUT / "intermediate" / "adt_df_non_icu_hosps.csv"
    if not f.exists():
        print(f"[{site}] {f} not found — no no-ICU CRRT encounters to plot.")
        return
    generate(pd.read_csv(f), site, OUT / "final" / "diagnostics" / "graphs")


if __name__ == "__main__":
    main()
