# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8"]
# ///
"""Figures for the README and the performance page, drawn from the suite's `results.tsv`.

    uv run benchmarks/suite/figures.py [benchmark_results/latest] [--out DIR]

Two kinds:

- **Accuracy:** missed variants (false negatives) on x against false calls (false positives) on y,
  one line per caller through its depths, both axes logarithmic. Down and to the left is better.
  Counts are vcfeval's with the genotype required to be right.
- **Performance:** wall time on x against peak memory on y, one line per caller through its
  configurations (depths, or numbers of samples), both axes logarithmic. Down and to the left is
  faster and smaller.

Only `results.tsv` is read, so the figures can be redrawn from any run. Static SVG and PNG, on the
light chart surface, for a markdown page.
"""

import argparse
import collections
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

# The reference palette's first three categorical slots (validated all-pairs, light mode). pvc's
# two ways of running share its hue and differ in line style and marker.
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e6e5e1"
STYLE = {
    "pvc (defaults)": dict(color="#2a78d6", linestyle="-", marker="o", label="pvc, one step"),
    "pvc (fitted)": dict(color="#2a78d6", linestyle="--", marker="s", label="pvc, fitted"),
    "GATK": dict(color="#eb6834", linestyle="-", marker="^", label="GATK"),
    "freebayes": dict(color="#1baf7a", linestyle="-", marker="D", label="freebayes"),
    # No fourth hue separates from these three for every reader (the palette validator fails every
    # candidate on all pairs), so DeepVariant, the accuracy reference, is drawn in neutral ink with
    # its own marker; against the three hues it passes both separation checks.
    "DeepVariant": dict(color="#3a3936", linestyle="-", marker="P", label="DeepVariant"),
}
CALLERS = list(STYLE)
DEPTH_ORDER = ["5x", "10x", "15x", "20x", "30x", "50x", "full", "300x"]


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def depth_key(d):
    return DEPTH_ORDER.index(d) if d in DEPTH_ORDER else 99


def depth_label(d):
    return {"full": "~60×", "300x": "~300×"}.get(d, d.replace("x", "×"))


def style_axes(ax, xlabel, ylabel, title):
    ax.set_facecolor(SURFACE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="major", color=GRID, linewidth=0.8)
    ax.grid(False, which="minor")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_2)
    ax.tick_params(colors=TEXT_2, labelsize=9)
    fmt = FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}")
    # Between decades, label only the 2s and 5s, in the same plain numbers as the decades.
    minor = FuncFormatter(lambda v, _: (f"{v:,.0f}" if v >= 1 else f"{v:g}")
                          if f"{v:.0e}"[0] in "25" else "")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(fmt)
        axis.set_minor_formatter(minor)
    ax.set_xlabel(xlabel, color=TEXT_2, fontsize=10)
    ax.set_ylabel(ylabel, color=TEXT_2, fontsize=10)
    ax.set_title(title, color=TEXT, fontsize=11, loc="left")


def draw_lines(ax, series, point_labels_for="pvc (defaults)", styles=None, order=None):
    """`series`: line -> list of (x, y, point label) in order."""
    styles = styles or STYLE
    for caller in order or CALLERS:
        pts = series.get(caller)
        if not pts:
            continue
        st = styles[caller]
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        ax.plot(xs, ys, color=st["color"], linestyle=st["linestyle"], linewidth=2,
                marker=st["marker"], markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.5,
                label=st["label"], zorder=3)
        # A direct label at the line's last point, in text ink, beside the colored mark.
        ax.annotate(st["label"], (xs[-1], ys[-1]), xytext=(6, 0), textcoords="offset points",
                    fontsize=8.5, color=TEXT, va="center", zorder=4, gid="caller")
        if caller == point_labels_for:
            for x, y, lab in pts:
                ax.annotate(lab, (x, y), xytext=(-7, -13), textcoords="offset points", fontsize=8,
                            color=TEXT_2, ha="right", zorder=4, gid="point")


def declutter(fig):
    """Move a text label the shortest way, up or down, until it overlaps no label already placed."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        placed = []
        # Callers' names are placed first; the points' labels give way to them.
        for text in sorted(ax.texts, key=lambda t: t.get_gid() != "caller"):
            x0, y0 = text.xyann
            for shift in (0, 10, -10, 20, -20, 30, -30, 40, -40):
                text.xyann = (x0, y0 + shift)
                fig.canvas.draw()
                box = text.get_window_extent(renderer).expanded(1.02, 1.1)
                if not any(box.overlaps(other) for other in placed):
                    break
            placed.append(text.get_window_extent(renderer))


def save(fig, out, name):
    declutter(fig)
    os.makedirs(out, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), facecolor=SURFACE, dpi=160,
                    bbox_inches="tight")
    plt.close(fig)
    return os.path.join(out, f"{name}.svg")


def accuracy_series(rows, benchmark, cls, aggregate_datasets=False):
    """caller -> [(missed, false calls, depth label)] for one benchmark and variant class."""
    acc = collections.defaultdict(lambda: collections.defaultdict(float))
    for r in rows:
        if (r["benchmark"] != benchmark or r["stage"] != "genotype match"
                or r["variant_class"] != cls or r["caller"] not in STYLE
                or r["metric"] not in ("missed", "false_calls")):
            continue
        key = (r["caller"], r["depth"]) if aggregate_datasets else (r["caller"], r["depth"], r["dataset"])
        acc[key][r["metric"]] += float(r["value"])
    series = collections.defaultdict(list)
    for key, m in acc.items():
        series[key[0]].append((key[1], m["missed"], m["false_calls"]))
    out = {}
    for caller, pts in series.items():
        pts.sort(key=lambda p: depth_key(p[0]))
        # A zero cannot sit on a log axis; it is drawn at 0.8, just under the 1 gridline.
        out[caller] = [(max(fn, 0.8), max(fp, 0.8), depth_label(d)) for d, fn, fp in pts]
    return out


def accuracy_figure(rows, benchmark, title, subtitle, out, name, aggregate=False):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), facecolor=SURFACE)
    for ax, cls in zip(axes, ("SNP", "indel")):
        series = accuracy_series(rows, benchmark, cls, aggregate)
        if not series:
            ax.set_visible(False)
            continue
        style_axes(ax, "missed variants (false negatives)", "false calls (false positives)",
                   f"{cls}s")
        draw_lines(ax, series)
    fig.suptitle(title, x=0.01, ha="left", color=TEXT, fontsize=13, y=1.04)
    fig.text(0.01, 0.97, subtitle, ha="left", color=TEXT_2, fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.06), labelcolor=TEXT)
    return save(fig, out, name)


# The performance figures' lines: (line, caller in the table, its stages in order of preference).
# pvc is three lines, its three ways of spending time: in one step from the reads, building the
# psps from the reads, and calling from the psps. Its parameter fit is optional, made once for a
# cohort, and is in no line (the summary reports it on its own).
PERF_LINES = [
    ("pvc, one step", "pvc (defaults)", ("call-from-alignments",)),
    ("pvc, building psps", "pvc (fitted)", ("generate-psps",)),
    ("pvc, calling from psps", "pvc (fitted)", ("call-from-psps",)),
    ("GATK", "GATK", ("whole pipeline", "HaplotypeCaller")),
    ("freebayes", "freebayes", ("freebayes",)),
    ("DeepVariant", "DeepVariant", ("DeepVariant + GLnexus", "DeepVariant")),
]
PERF_STYLE = dict(STYLE)
PERF_STYLE["pvc, one step"] = dict(STYLE["pvc (defaults)"])
PERF_STYLE["pvc, building psps"] = dict(color="#2a78d6", linestyle="--", marker="s",
                                        label="pvc, building psps")
PERF_STYLE["pvc, calling from psps"] = dict(color="#2a78d6", linestyle=":", marker="v",
                                            label="pvc, calling from psps")


def performance_series(rows, benchmark, label_of):
    vals = collections.defaultdict(dict)
    for line, caller, stages in PERF_LINES:
        for r in rows:
            if r["benchmark"] != benchmark or r["caller"] != caller or r["stage"] not in stages:
                continue
            if r["metric"] not in ("wall_seconds", "peak_rss_mb_all_processes"):
                continue
            key = (line, r["samples"], r["depth"])
            rank = stages.index(r["stage"])
            prior = vals[key].get(("rank", r["metric"]))
            if prior is not None and prior < rank:
                continue
            vals[key][r["metric"]] = float(r["value"])
            vals[key][("rank", r["metric"])] = rank
    series = collections.defaultdict(list)
    for (line, samples, depth), m in vals.items():
        if "wall_seconds" in m and "peak_rss_mb_all_processes" in m:
            series[line].append((label_of(samples, depth), m["wall_seconds"],
                                 m["peak_rss_mb_all_processes"] / 1024))
    return series


def performance_figure(rows, benchmark, order_key, label_of, title, subtitle, out, name):
    raw = performance_series(rows, benchmark, label_of)
    series = {c: [(w, mem, depth_label(lab) if lab in DEPTH_ORDER else lab)
                  for lab, w, mem in sorted(pts, key=lambda p: order_key(p[0]))]
              for c, pts in raw.items()}
    fig, ax = plt.subplots(figsize=(7.5, 5), facecolor=SURFACE)
    style_axes(ax, "wall time (seconds)", "peak memory (GB)", "")
    ax.margins(x=0.08, y=0.12)
    draw_lines(ax, series, point_labels_for="pvc, one step", styles=PERF_STYLE,
               order=[line for line, _, _ in PERF_LINES])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}" if f"{v:.0e}"[0] in "25" else ""))
    fig.suptitle(title, x=0.01, ha="left", color=TEXT, fontsize=13, y=1.02)
    fig.text(0.01, 0.955, subtitle, ha="left", color=TEXT_2, fontsize=9)
    ax.legend(loc="upper left", frameon=False, fontsize=9, labelcolor=TEXT)
    return save(fig, out, name)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run", nargs="?", default="benchmark_results/latest")
    p.add_argument("--out", default=None)
    args = p.parse_args()
    rows = read_rows(os.path.join(args.run, "results.tsv"))
    out = args.out or os.path.join(args.run, "figures")
    made = [
        accuracy_figure(rows, "giab_chr20", "Accuracy — HG002, chromosome 20",
                        "GIAB v4.2.1, confident regions; vcfeval, genotype must match; "
                        "points from 5× to ~60×, lower left is better", out, "accuracy_giab_chr20"),
        accuracy_figure(rows, "giab_panel", "Accuracy — GIAB HG002, HG003, HG004 panel",
                        "100 confident regions per sample, summed over the three; points from 5× "
                        "to ~300×; a count of zero is drawn just below 1", out,
                        "accuracy_giab_panel", aggregate=True),
        performance_figure(
            rows, "giab_chr20", depth_key, lambda s, d: d,
            # depth_key orders the raw depths; their labels are made readable below.
            "Speed and memory — HG002, chromosome 20, one sample",
            "6 threads; points from 5× to ~60×; pvc's optional parameter fit is in no line",
            out,
            "performance_giab_chr20"),
        performance_figure(
            rows, "speed", lambda s: int(s), lambda s, d: s,
            "Speed and memory — tomato cohort, 8 Mb at ~3 reads a position",
            "6 threads; points are 1, 8 and 62 accessions; pvc's optional parameter fit is in no "
            "line", out, "performance_tomato"),
    ]
    for m in made:
        print(m)


if __name__ == "__main__":
    main()
