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
y-axis): all_throughput_boxplot.{pdf,png}, all_thread_sweep.{pdf,png}

Needs matplotlib and pandas.
"""

import argparse
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

# The variants in the order of config.sh (fallback: as they appear).
VARIANT_ORDER = ["sa-binary-search", "lrf-ms", "pt16-mlr", "pt16-v2-mlr", "sassy-mlr",
                 "pt16", "pt16-v2", "sassy", "varki", "powered-fwd-escape",
                 "powered-pt16-fwd-escape"]


# The names in the paper (the data keeps the programs' names).
# The PT16-SA variants are the mlr versions (range hits narrowed by the mlr
# binary search); the older versions narrowing one character at a time are
# marked "(char)".
LABELS = {
    "sa-binary-search": "SA",
    "lrf-ms": "LRF",
    "pt16-mlr": "PT16-SA",
    "pt16-v2-mlr": "PT16-SA-v2",
    "sassy-mlr": "PT16-SA-sassy",
    "pt16": "PT16-SA (char)",
    "pt16-v2": "PT16-SA-v2 (char)",
    "sassy": "PT16-SA-sassy (char)",
    "varki": "Varki",
    "powered-fwd-escape": "Powered",
    "powered-pt16-fwd-escape": "PT16-Pow",
}


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
    known = [v for v in VARIANT_ORDER if v in values]
    return known + sorted(v for v in values if v not in VARIANT_ORDER)


def box_data(dataset, per_file):
    """The 1-thread per-file rows of a dataset, with their throughput."""
    data = per_file[(per_file.dataset == dataset) & (per_file.threads == 1)].copy()
    # Throughput per file: bytes per millisecond / 1000 = MB/s (10^6 bytes).
    data["mb_per_s"] = data.input_bytes / data.parse_ms / 1000.0
    return data


def draw_boxes(ax, dataset, data, variants, ylabel=True):
    """Plot 1 into `ax`: per variant, one box per reference. One colour for
    every box; the references are told apart by position (R1, R2, ... from
    left to right within each variant, labelled below the boxes and keyed in
    the corner)."""
    refs = sorted(set(data.reference))
    width = 0.22
    spacing = width * 1.1  # the boxes of a variant almost touch
    for k, ref in enumerate(refs):
        values = [data[(data.parser == v) & (data.reference == ref)].mb_per_s.values
                  for v in variants]
        positions = [i + (k - (len(refs) - 1) / 2) * spacing for i in range(len(variants))]
        box = ax.boxplot(values, positions=positions, widths=width,
                         patch_artist=True, flierprops={"markersize": 2},
                         medianprops={"color": "black"})
        for patch in box["boxes"]:
            patch.set_facecolor("#4C72B0")
            patch.set_alpha(0.75)
        for x in positions:
            ax.text(x, 0.01, f"R{k + 1}", transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=6, color="dimgray")
    key = "\n".join(f"R{k + 1}: {short(ref)}" for k, ref in enumerate(refs))
    ax.text(0.01, 0.98, key, transform=ax.transAxes, ha="left", va="top", fontsize=7,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "lightgray"})
    ax.set_xticks(range(len(variants)))
    ax.set_xticklabels([label(v) for v in variants], rotation=35, ha="right",
                       rotation_mode="anchor")
    ax.set_xlim(-0.5, len(variants) - 0.5)
    ax.set_yscale("log")
    # Short labels: details belong in the caption (1 thread, one box per
    # reference over the files).
    if ylabel:
        ax.set_ylabel("throughput per file (MB/s)")
    ax.set_title(dataset)
    ax.grid(axis="y", which="both", alpha=0.3)


def plot_boxes(dataset, per_file, out, notes):
    data = box_data(dataset, per_file)
    if data.empty:
        print(f"  {dataset}: no 1-thread per-file rows")
        return
    variants = ordered(set(data.parser))

    # Compact: the three datasets' plots go side by side.
    fig, ax = plt.subplots(figsize=(0.8 * len(variants) + 1.2, 4.5))
    draw_boxes(ax, dataset, data, variants)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{dataset}_throughput_boxplot.{ext}", dpi=150)
    plt.close(fig)

    medians = data.groupby(["parser", "reference"]).mb_per_s.median().unstack()
    medians = medians.reindex(variants)
    medians.index = [label(v) for v in medians.index]
    medians.columns = [short(c) for c in medians.columns]
    notes.append(f"{dataset}: median parse throughput per file (MB/s), 1 thread\n"
                 f"{medians.round(1).to_string()}\n")


def sweep_data(dataset, summary):
    return summary[(summary.dataset == dataset) & (summary.sweep_reference == 1)]


def draw_sweep(ax, dataset, data, variants, ylabel=True, legend=True):
    """Plot 2 into `ax`: wall time against threads, one line per variant,
    each with its ideal linear scaling (dashed: 1-thread time / threads)."""
    reference = data.reference.iloc[0]
    for variant in variants:
        rows = data[data.parser == variant].sort_values("threads")
        if rows.empty:
            continue
        line, = ax.plot(rows.threads, rows.wall_ms / 1000.0, marker="o",
                        markersize=4, label=label(variant))
        # Ideal (linear) scaling from the 1-thread time: wall(1) / threads.
        one = rows[rows.threads == 1]
        if len(one):
            ax.plot(rows.threads, one.wall_ms.iloc[0] / 1000.0 / rows.threads,
                    linestyle="--", linewidth=1, color=line.get_color(), alpha=0.6)
    ax.plot([], [], linestyle="--", linewidth=1, color="gray",
            label="ideal (1-thread time / threads)")
    threads = sorted(data.threads.unique())
    ax.set_xscale("log", base=2)
    ax.set_xticks(threads)
    ax.set_xticklabels([str(t) for t in threads])
    ax.set_yscale("log")
    ax.set_xlabel("threads")
    if ylabel:
        ax.set_ylabel("parse time of the collection (s)")
    ax.set_title(f"{dataset} (reference {short(reference)})", fontsize=10)
    if legend:
        ax.legend(fontsize=7)
    ax.grid(which="both", alpha=0.3)


def plot_sweep(dataset, summary, out, notes):
    data = sweep_data(dataset, summary)
    if data.empty:
        print(f"  {dataset}: no sweep rows")
        return
    variants = ordered(set(data.parser))
    reference = data.reference.iloc[0]

    fig, ax = plt.subplots(figsize=(7, 5))
    draw_sweep(ax, dataset, data, variants)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{dataset}_thread_sweep.{ext}", dpi=150)
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
        fig, axes = plt.subplots(1, len(boxes), sharey=True,
                                 figsize=(0.62 * len(variants) * len(boxes) + 1.2, 4.2))
        axes = [axes] if len(boxes) == 1 else list(axes)
        for k, (ax, (dataset, data)) in enumerate(zip(axes, boxes)):
            draw_boxes(ax, dataset, data, variants, ylabel=(k == 0))
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.03)  # panels close together (shared y-axis)
        for ext in ("pdf", "png"):
            fig.savefig(out / f"all_throughput_boxplot.{ext}", dpi=150)
        plt.close(fig)

    sweeps = [(d, sweep_data(d, summary)) for d in datasets]
    sweeps = [(d, s) for d, s in sweeps if not s.empty]
    if sweeps:
        variants = ordered(set().union(*(set(s.parser) for _, s in sweeps)))
        fig, axes = plt.subplots(1, len(sweeps), sharey=True,
                                 figsize=(4.0 * len(sweeps) + 2.2, 4.0))
        axes = [axes] if len(sweeps) == 1 else list(axes)
        for k, (ax, (dataset, data)) in enumerate(zip(axes, sweeps)):
            draw_sweep(ax, dataset, data, variants, ylabel=(k == 0), legend=False)
            ax.tick_params(axis="y", which="both", left=True)
        # One legend for all panels, to the right.
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.04)
        for ext in ("pdf", "png"):
            fig.savefig(out / f"all_thread_sweep.{ext}", dpi=150, bbox_inches="tight")
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
