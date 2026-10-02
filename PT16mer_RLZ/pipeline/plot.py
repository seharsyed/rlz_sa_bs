#!/usr/bin/env python3
"""The experiment's plots, from collect.sh's two CSV files.

    ./plot.py [--results DIR] [--out DIR] [DATASET ...]

DIR holds per_file.csv and summary.csv (collect.sh; default: config.sh's
RESULTS). Per dataset (default: every dataset in the files), into --out
(default: DIR/plots):

  <dataset>_throughput_boxplot.{pdf,png}  plot 1: per variant, one box per
        reference, each over the parse throughput of the collection's files
        (input size / parse time, MB/s with MB = 10^6 bytes; 1 thread)
  <dataset>_thread_sweep.{pdf,png}        plot 2: the whole collection's parse
        (wall time) against the thread count, one line per variant (the
        dataset's sweep reference), each with its ideal linear scaling
        (dashed: 1-thread time / threads)
  <dataset>_medians.txt                   the plotted numbers, for a look

and, with several datasets, both plots with the datasets side by side (shared
y-axis): all_throughput_boxplot.{pdf,png}, all_thread_sweep.{pdf,png}, and the
thread sweep as throughput (MB/s of the collection): all_thread_sweep_throughput

Needs matplotlib and pandas.
"""

import argparse
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

# A LaTeX-like serif font (matplotlib's own Computer Modern).
# Figures are drawn at their printed size (A4 portrait, text width: about
# 16 cm = 6.3 in), so these font sizes are the sizes on the page.
TEXT_WIDTH = 6.3
plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "font.family": "serif",
    "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,
    "axes.unicode_minus": False,
})

# ggplot2-like panels (the boxplots): grey background, white grid, a grey
# header strip with the panel's name, no frame.
PANEL_BACKGROUND = "#EBEBEB"
STRIP_BACKGROUND = "#D9D9D9"


def ggplot_panel(ax, title):
    ax.set_facecolor(PANEL_BACKGROUND)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(axis="y", which="major", color="white", linewidth=1.0)
    ax.grid(axis="y", which="minor", color="white", linewidth=0.5)
    ax.grid(axis="x", which="major", color="white", linewidth=1.0)
    ax.set_axisbelow(True)
    ax.tick_params(which="both", length=3, color="#333333")
    ax.tick_params(which="minor", length=0)
    # The header strip: a grey band above the panel, as wide as the panel.
    ax.add_patch(Rectangle((0, 1.0), 1, 0.085, transform=ax.transAxes,
                           facecolor=STRIP_BACKGROUND, edgecolor="none",
                           clip_on=False))
    ax.text(0.5, 1.0425, title, transform=ax.transAxes, ha="center", va="center",
            fontsize=9)


def tex_text(text):
    """Text in mathtext's roman font, so that "_" shows (Computer Modern's
    text font has no underscore)."""
    return "$\\mathrm{" + text.replace("_", "\\_") + "}$"


def plain_log_ticks(ax):
    """Log-axis labels as plain numbers (10, 100, 1000) instead of powers."""
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{v:g}" if v < 1e4 else f"{v / 1e3:g}k"))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())

# The variants in the order of config.sh (fallback: as they appear).
# The baselines (SA, LRF, Varki), then Powered, then the PT16 tables, ending
# with PT16-Pow.
VARIANT_ORDER = ["sa-binary-search", "lrf-ms", "varki", "powered-fwd-escape",
                 "pt16-mlr", "pt16", "pt16-v2-mlr", "pt16-v2", "sassy-mlr", "sassy",
                 "powered-pt16-fwd-escape"]


# The names in the paper (the data keeps the programs' names). The PT16-SA
# variants with and without the mlr binary search get the same name: a run
# has one of them (the results folder tells which).
LABELS = {
    "sa-binary-search": "SA",
    "lrf-ms": "LRF",
    "pt16-mlr": "PT16-SA",
    "pt16-v2-mlr": "PT16-SA-v2",
    "sassy-mlr": "PT16-SA",
    "pt16": "PT16-SA",
    "pt16-v2": "PT16-SA-v2",
    "sassy": "PT16-SA",
    "varki": "Varki",
    "powered-fwd-escape": "Powered",
    "powered-pt16-fwd-escape": "PT16-Pow",
}


# The variants in the plots: of the three PT16 tables on the suffix-array side
# only sassy (the fastest, never significantly slower than the other two),
# named PT16-SA. The others stay in the data.
SHOWN = {"sa-binary-search", "lrf-ms", "varki", "powered-fwd-escape", "sassy-mlr",
         "sassy", "powered-pt16-fwd-escape"}


def label(variant):
    """A variant's name in the plots."""
    return LABELS.get(variant, variant)


def default_results():
    here = Path(__file__).resolve().parent
    out = subprocess.run(["bash", "-c", f'source "{here}/config.sh"; echo "$RESULTS"'],
                         check=True, capture_output=True, text=True).stdout.strip()
    return Path(out)


def short(reference):
    """A reference's label: its accession (GCA_903819205.2_..._genomic ->
    GCA_903819205.2), else its name."""
    parts = reference.split("_")
    if reference.startswith("GC") and len(parts) >= 2:
        return parts[0] + "_" + parts[1]
    return reference


def ordered(values):
    """The shown variants among `values`, in VARIANT_ORDER."""
    values = {v for v in values if v in SHOWN}
    known = [v for v in VARIANT_ORDER if v in values]
    return known + sorted(v for v in values if v not in VARIANT_ORDER)


def box_data(dataset, per_file):
    """The 1-thread per-file rows of a dataset, with their throughput."""
    data = per_file[(per_file.dataset == dataset) & (per_file.threads == 1)].copy()
    # Throughput per file: bytes per millisecond / 1000 = MB/s (10^6 bytes).
    data["mb_per_s"] = data.input_bytes / data.parse_ms / 1000.0
    return data


def draw_boxes(ax, dataset, data, variants, ylabel=True):
    """Plot 1 into `ax`, in ggplot2 style: per variant, one white box per
    reference (1, 2, 3 from left to right within each variant, labelled below
    the boxes)."""
    refs = sorted(set(data.reference))
    width = 0.22
    spacing = width * 1.1  # the boxes of a variant almost touch
    line = {"color": "#333333", "linewidth": 0.8}
    for k, ref in enumerate(refs):
        values = [data[(data.parser == v) & (data.reference == ref)].mb_per_s.values
                  for v in variants]
        positions = [i + (k - (len(refs) - 1) / 2) * spacing for i in range(len(variants))]
        ax.boxplot(values, positions=positions, widths=width, patch_artist=True,
                   boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.8},
                   whiskerprops=line, capprops={"linewidth": 0},
                   medianprops={"color": "#333333", "linewidth": 1.6},
                   flierprops={"marker": "o", "markersize": 2.5, "markerfacecolor": "#333333",
                               "markeredgecolor": "#333333"})
        for x in positions:
            ax.text(x, 0.012, f"{k + 1}", transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=6.5, color="#666666")
    # The references are numbered 1, 2, ... (in name order; a table in the
    # paper names them).
    ax.set_xticks(range(len(variants)))
    ax.set_xticklabels([label(v) for v in variants], rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_xlim(-0.5, len(variants) - 0.5)
    ax.set_yscale("log")
    plain_log_ticks(ax)
    # Short labels: details belong in the caption (1 thread, one box per
    # reference over the files).
    if ylabel:
        ax.set_ylabel("throughput per file (MB/s)")
    ggplot_panel(ax, dataset)


def reference_legend(ax, count):
    """The legend of the numbers under the boxes, in `ax`'s empty top-left
    corner."""
    numbers = ", ".join(str(k + 1) for k in range(count))
    ax.text(0.03, 0.97, f"boxes {numbers}: reference {numbers}", transform=ax.transAxes,
            ha="left", va="top", fontsize=7,
            bbox={"boxstyle": "square,pad=0.3", "facecolor": "white", "edgecolor": "none"})


def plot_boxes(dataset, per_file, out, notes):
    data = box_data(dataset, per_file)
    if data.empty:
        print(f"  {dataset}: no 1-thread per-file rows")
        return
    variants = ordered(set(data.parser))

    # Compact: the three datasets' plots go side by side.
    fig, ax = plt.subplots(figsize=(TEXT_WIDTH / 2, 2.8))  # half the text width
    draw_boxes(ax, dataset, data, variants)
    reference_legend(ax, data.reference.nunique())
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{dataset}_throughput_boxplot.{ext}", dpi=150, bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)

    medians = data.groupby(["parser", "reference"]).mb_per_s.median().unstack()
    medians = medians.reindex(variants)
    medians.index = [label(v) for v in medians.index]
    medians.columns = [short(c) for c in medians.columns]
    notes.append(f"{dataset}: median parse throughput per file (MB/s), 1 thread\n"
                 f"{medians.round(1).to_string()}\n")


def sweep_data(dataset, summary):
    return summary[(summary.dataset == dataset) & (summary.sweep_reference == 1)]


def draw_sweep(ax, dataset, data, variants, ylabel=True, legend=True,
               throughput=False):
    """Plot 2 into `ax`: wall time (or, with throughput, the collection's
    throughput in MB/s) against threads, one line per variant, each with its
    ideal linear scaling (dashed)."""
    reference = data.reference.iloc[0]
    for variant in variants:
        rows = data[data.parser == variant].sort_values("threads")
        if rows.empty:
            continue
        if throughput:
            # The whole collection: bytes per millisecond / 1000 = MB/s.
            values = rows.input_bytes / rows.wall_ms / 1000.0
        else:
            values = rows.wall_ms / 1000.0
        line, = ax.plot(rows.threads, values, marker="o", markersize=4,
                        label=label(variant))
        # Ideal (linear) scaling from 1 thread: time / threads, throughput * threads.
        one = values[rows.threads == 1]
        if len(one):
            ideal = one.iloc[0] * rows.threads if throughput else one.iloc[0] / rows.threads
            ax.plot(rows.threads, ideal, linestyle="--", linewidth=1,
                    color=line.get_color(), alpha=0.6)
    # "ideal": the dashed lines, 1-thread time / threads (see the caption).
    ax.plot([], [], linestyle="--", linewidth=1, color="gray", label="ideal")
    threads = sorted(data.threads.unique())
    ax.set_xscale("log", base=2)
    ax.set_xticks(threads)
    ax.set_xticklabels([str(t) for t in threads])
    ax.set_yscale("log")
    ax.set_xlabel("threads")
    ax.minorticks_off()
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    if ylabel:
        ax.set_ylabel("throughput (MB/s)" if throughput else "parse time of the collection (s)")
    plain_log_ticks(ax)
    if legend:
        ax.legend(fontsize=7)
    # Same look as the boxplots; the sweep's reference (the dataset's first)
    # belongs in the caption.
    ggplot_panel(ax, dataset)
    ax.grid(axis="y", which="minor", color="white", linewidth=0.5)


def plot_sweep(dataset, summary, out, notes):
    data = sweep_data(dataset, summary)
    if data.empty:
        print(f"  {dataset}: no sweep rows")
        return
    variants = ordered(set(data.parser))
    reference = data.reference.iloc[0]

    fig, ax = plt.subplots(figsize=(TEXT_WIDTH / 2, 2.8))
    draw_sweep(ax, dataset, data, variants)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{dataset}_thread_sweep.{ext}", dpi=150, bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)

    table = data.pivot_table(index="parser", columns="threads", values="wall_ms") / 1000.0
    table = table.reindex(variants)
    table.index = [label(v) for v in table.index]
    notes.append(f"{dataset}: wall time of the collection (s), reference {short(reference)}\n"
                 f"{table.round(2).to_string()}\n")


def plot_combined(datasets, per_file, summary, out):
    """The datasets side by side, sharing the y-axis: all_throughput_boxplot
    (plot 1) and all_thread_sweep (plot 2)."""
    boxes = [(d, box_data(d, per_file)) for d in datasets]
    boxes = [(d, b) for d, b in boxes if not b.empty]
    if boxes:
        variants = ordered(set().union(*(set(b.parser) for _, b in boxes)))
        fig, axes = plt.subplots(1, len(boxes), sharey=True, figsize=(TEXT_WIDTH, 2.9))
        axes = [axes] if len(boxes) == 1 else list(axes)
        for k, (ax, (dataset, data)) in enumerate(zip(axes, boxes)):
            draw_boxes(ax, dataset, data, variants, ylabel=(k == 0))
        reference_legend(axes[0], boxes[0][1].reference.nunique())
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.04)  # panels close together, as ggplot facets
        for ext in ("pdf", "png"):
            fig.savefig(out / f"all_throughput_boxplot.{ext}", dpi=150, bbox_inches="tight", pad_inches=0.01)
        plt.close(fig)

    sweeps = [(d, sweep_data(d, summary)) for d in datasets]
    sweeps = [(d, s) for d, s in sweeps if not s.empty]
    if sweeps:
        variants = ordered(set().union(*(set(s.parser) for _, s in sweeps)))
        for throughput, name in ((False, "all_thread_sweep"),
                                 (True, "all_thread_sweep_throughput")):
            fig, axes = plt.subplots(1, len(sweeps), sharey=True,
                                     figsize=(TEXT_WIDTH, 2.3))
            axes = [axes] if len(sweeps) == 1 else list(axes)
            for k, (ax, (dataset, data)) in enumerate(zip(axes, sweeps)):
                draw_sweep(ax, dataset, data, variants, ylabel=(k == 0), legend=False,
                           throughput=throughput)
                ax.tick_params(axis="y", which="both", left=True)
            # One legend for all panels, below them.
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.06),
                       ncol=len(labels), frameon=False, handlelength=1.6,
                       columnspacing=1.0, handletextpad=0.4)
            fig.tight_layout()
            fig.subplots_adjust(wspace=0.06)
            for ext in ("pdf", "png"):
                fig.savefig(out / f"{name}.{ext}", dpi=150, bbox_inches="tight",
                            pad_inches=0.01)
            plt.close(fig)


def dataset_order(names, per_file):
    """The datasets by increasing size (the collection's total input bytes,
    for one reference and parser at 1 thread)."""
    one = per_file[per_file.threads == 1]
    size = {}
    for d in names:
        rows = one[one.dataset == d]
        if rows.empty:
            size[d] = float("inf")
            continue
        first = rows[(rows.reference == rows.reference.iloc[0]) &
                     (rows.parser == rows.parser.iloc[0])]
        size[d] = first.input_bytes.sum()
    return sorted(names, key=lambda d: (size[d], d))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("datasets", nargs="*")
    args = parser.parse_args()

    results = args.results or default_results()
    out = args.out or results / "plots"
    out.mkdir(parents=True, exist_ok=True)
    per_file = pd.read_csv(results / "per_file.csv")
    summary = pd.read_csv(results / "summary.csv")

    datasets = args.datasets or dataset_order(set(per_file.dataset) | set(summary.dataset),
                                              per_file)
    for dataset in datasets:
        notes = []
        plot_boxes(dataset, per_file, out, notes)
        plot_sweep(dataset, summary, out, notes)
        (out / f"{dataset}_medians.txt").write_text("\n".join(notes))
        print(f"  {dataset}: plots and numbers in {out}")
    if len(datasets) > 1:
        plot_combined(datasets, per_file, summary, out)
        print(f"  all: all_throughput_boxplot, all_thread_sweep in {out}")


if __name__ == "__main__":
    main()
