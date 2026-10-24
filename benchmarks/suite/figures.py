# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib>=3.8"]
# ///
"""Figures for the README and the performance page, drawn from the suite's `results.tsv`.

    uv run benchmarks/suite/figures.py [benchmark_results/latest] [--out DIR]

Two kinds:

- **Errors:** a dot plot of the three error types — wrong GT (truth ref), wrong GT (truth var), no
  call (./.) — one row per caller and depth on one logarithmic count axis, one panel per class of
  variant (the section "accuracy: the three error types" below). Counts are pvc_eval's for the
  GIAB benchmarks and the simulated-STR scorer's for the stutter grid.
- **Performance:** wall time on x against peak memory on y, one line per caller through its
  configurations (depths, or numbers of samples), both axes logarithmic. Down and to the left is
  faster and smaller.

Only `results.tsv` is read, so the figures can be redrawn from any run. Static SVG and PNG, on the
light chart surface, for a markdown page.
"""

import argparse
import collections
import csv
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter, NullFormatter  # noqa: E402

# The reference palette's first three categorical slots (validated all-pairs, light mode). pvc's
# two ways of running share its hue and differ in line style and marker.
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e6e5e1"
STYLE = {
    "pvc (defaults)": dict(color="#2a78d6", linestyle="-", marker="o", label="pvc, defaults"),
    "pvc (fitted)": dict(color="#2a78d6", linestyle="--", marker="s", label="pvc, fitted"),
    "GATK": dict(color="#eb6834", linestyle="-", marker="^", label="GATK"),
    "freebayes": dict(color="#1baf7a", linestyle="-", marker="D", label="freebayes"),
    # No fourth hue separates from these three for every reader (the palette validator fails every
    # candidate on all pairs), so DeepVariant, the accuracy reference, is drawn in neutral ink with
    # its own marker; against the three hues it passes both separation checks.
    "DeepVariant": dict(color="#3a3936", linestyle="-", marker="P", label="DeepVariant"),
    # HipSTR, only in the repeat benchmarks, takes a lighter neutral on the same reasoning; its two
    # runs differ by line style and marker, and every line carries a direct label.
    "HipSTR": dict(color="#8a8883", linestyle="-", marker="X", label="HipSTR"),
    "HipSTR (fixed stutter)": dict(color="#8a8883", linestyle="--", marker="v",
                                   label="HipSTR, fixed stutter"),
}
CALLERS = list(STYLE)
DEPTH_ORDER = ["5x", "10x", "15x", "20x", "30x", "50x", "60x", "full", "300x"]


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def depth_key(d):
    return DEPTH_ORDER.index(d) if d in DEPTH_ORDER else 99


def depth_label(d):
    return {"full": "~60×", "300x": "~300×"}.get(d, d.replace("x", "×"))


def is_two_or_five(v):
    """Whether a tick is a 2 or a 5 of its decade (2,000, 50, 0.2), and not, say, 16,000."""
    if v <= 0:
        return False
    mantissa = v / 10 ** math.floor(math.log10(v) + 1e-9)
    return any(abs(mantissa - m) < 1e-6 for m in (2, 5))


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
                          if is_two_or_five(v) else "")
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


# ---- accuracy: the three error types ---------------------------------------------------------
#
# Every accuracy figure is a dot plot (designed with the owner, 24 October 2026). One row per
# caller, grouped by depth; one shared count axis, logarithmic, with a short linear stretch from 0
# to 1 where some count is zero; one panel per class of variant. A row's three error counts are
# three marks on its line, nested by size on one lane so tied counts show as a target, and drawn
# in the caller's colour so a colour means the same caller in every figure:
#
# | mark                | error type           | what it counts                                   |
# |---------------------|----------------------|--------------------------------------------------|
# | medium filled cross | wrong GT (truth ref) | a variant called where the truth is reference    |
# | large open ring     | wrong GT (truth var) | the truth has a variant and the call is not it   |
# | small light dot     | no call (./.)        | the genotype left missing, the lesser error      |
#
# The line joins all three, zeros included: a short line near 0 is the good result.

ERROR_MARKS = [
    ("wrong_gt_truth_ref", "wrong GT (truth ref)",
     dict(marker="X", fill=True, size=9.5, z=4, edge=1.2, alpha=1)),
    ("wrong_gt_truth_var", "wrong GT (truth var)",
     dict(marker="o", fill=False, size=11.5, z=3, edge=1.8, alpha=1)),
    ("no_call", "no call (./.)",
     dict(marker="o", fill=True, size=6, z=5, edge=1.0, alpha=0.8)),
]
ERROR_SUBTITLE = ("wrong GT (truth ref): a variant called where the truth is reference; "
                  "wrong GT (truth var): the truth has a variant and the call is not its genotype; "
                  "no call (./.): genotype left missing")


def error_counts(rows, benchmark, kind, aggregate_datasets=False, exclude=()):
    """{(depth, caller): {error type: count}} from pvc_eval's rows, for one class of variant.

    The loci scored right are counted too, under "right", so that a class with no locus at all
    (no STR in the trio's regions) can be told from one with no error, and left out."""
    acc = collections.defaultdict(lambda: collections.defaultdict(float))
    wanted = {metric for metric, _, _ in ERROR_MARKS} | {"right"}
    datasets = {r["dataset"] for r in rows if r["benchmark"] == benchmark}
    if len(datasets) > 1 and not aggregate_datasets:
        raise ValueError(f"{benchmark} has {len(datasets)} datasets; sum them or pick one")
    for r in rows:
        if (r["benchmark"] == benchmark and r["stage"] == "error type"
                and r["variant_class"] == kind and r["caller"] in STYLE
                and r["caller"] not in exclude and r["metric"] in wanted):
            acc[(r["depth"], r["caller"])][r["metric"]] += float(r["value"])
    if not any(v for counts in acc.values() for v in counts.values()):
        return {}
    return {k: {m: v for m, v in c.items() if m != "right"} for k, c in acc.items()}


def row_layout(panels):
    """One row position for every (depth, caller) any panel has, shared by every panel, so a row
    reads across them even where a panel lacks a caller (HipSTR and the STR bundles)."""
    keys = set().union(*(counts.keys() for _, counts in panels))
    depths = sorted({d for d, _ in keys}, key=depth_key)
    y, place, groups = 0, {}, []
    for depth in depths:
        top = y
        for caller in CALLERS:
            if (depth, caller) in keys:
                place[(depth, caller)] = y
                y += 1
        groups.append((depth_label(depth), top))
        y += 0.8
    return place, groups, y - 0.8


def dot_panel(ax, counts, layout, show_labels):
    place, groups, bottom = layout
    for key, y in place.items():
        if key not in counts:
            continue
        colour = STYLE[key[1]]["color"]
        values = [counts[key].get(metric, 0) for metric, _, _ in ERROR_MARKS]
        ax.plot([min(values), max(values)], [y, y], color=colour, linewidth=1.2, alpha=0.6,
                zorder=2)
        for (metric, _, mark), value in zip(ERROR_MARKS, values):
            ax.plot(value, y, linestyle="none", marker=mark["marker"], markersize=mark["size"],
                    markerfacecolor=colour if mark["fill"] else "none",
                    markeredgecolor=SURFACE if mark["fill"] else colour,
                    markeredgewidth=mark["edge"], alpha=mark["alpha"], zorder=mark["z"])
    ax.set_facecolor(SURFACE)
    values = [v for c in counts.values() for v in c.values()] or [1]
    if min(values) == 0:
        # A zero has no place on a log axis: a short linear stretch from 0 to 1 holds it.
        ax.set_xscale("symlog", linthresh=1, linscale=0.4)
        ax.set_xlim(left=-0.3)
    else:
        ax.set_xscale("log")
        ax.set_xlim(left=10 ** math.floor(math.log10(min(values))))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else "0"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(bottom + 0.5, -0.7)
    ax.set_yticks(list(place.values()))
    ax.set_yticklabels([STYLE[caller]["label"] for _, caller in place])
    ax.tick_params(axis="y", length=0, labelsize=8.5, labelcolor=TEXT, labelleft=show_labels)
    ax.tick_params(axis="x", colors=TEXT_2, labelsize=9)
    ax.grid(True, axis="x", which="major", color=GRID, linewidth=0.8)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(TEXT_2)
    for label, top in groups:
        if top:
            ax.axhline(top - 0.9, color=GRID, linewidth=0.8, zorder=1)
        if show_labels:
            ax.annotate(label, (0, top), xycoords=("axes fraction", "data"), xytext=(-125, 0),
                        textcoords="offset points", ha="left", va="center", fontsize=9,
                        color=TEXT, fontweight="bold")


def error_figure(panels, title, subtitle, out, name):
    """`panels`: [(panel title, {(depth, caller): {error type: count}})]; empty panels dropped."""
    panels = [(t, c) for t, c in panels if c]
    if not panels:
        return None
    layout = row_layout(panels)
    fig, axes = plt.subplots(1, len(panels), figsize=(4.4 * len(panels) + 1.6,
                                                      1.4 + 0.27 * len(layout[0])),
                             facecolor=SURFACE, sharey=True, squeeze=False)
    for i, ((panel, counts), ax) in enumerate(zip(panels, axes[0])):
        dot_panel(ax, counts, layout, show_labels=(i == 0))
        ax.set_title(panel, loc="left", color=TEXT, fontsize=11)
        ax.set_xlabel("errors (log scale)", color=TEXT_2, fontsize=9)
    fig.tight_layout(rect=(0.04, 0.05, 1, 1))
    handles = [Line2D([], [], linestyle="none", marker=m["marker"], markersize=m["size"],
                      markeredgewidth=m["edge"], alpha=m["alpha"],
                      markerfacecolor=TEXT_2 if m["fill"] else "none",
                      markeredgecolor=SURFACE if m["fill"] else TEXT_2, label=label)
               for _, label, m in ERROR_MARKS]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.01), labelcolor=TEXT)
    pt = 1 / 72 / fig.get_figheight()
    fig.text(0.01, 1 + 4 * pt, ERROR_SUBTITLE, ha="left", va="bottom", color=TEXT_2, fontsize=8)
    fig.text(0.01, 1 + 17 * pt, subtitle, ha="left", va="bottom", color=TEXT_2, fontsize=9)
    fig.text(0.01, 1 + 32 * pt, title, ha="left", va="bottom", color=TEXT, fontsize=13)
    os.makedirs(out, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), facecolor=SURFACE, dpi=160,
                    bbox_inches="tight")
    plt.close(fig)
    return os.path.join(out, f"{name}.svg")


# pvc_eval's classes of variant, one panel each.
KIND_PANELS = [("SNPs", "SNP"), ("indels", "indel"), ("STRs", "STR"), ("complex", "complex")]


def vcfeval_benchmark_figure(rows, benchmark, title, subtitle, out, name, aggregate=False,
                             exclude=()):
    panels = [(title_, error_counts(rows, benchmark, kind, aggregate, exclude))
              for title_, kind in KIND_PANELS]
    return error_figure(panels, title, subtitle, out, name)


# The simulated stuttering benchmark scores with its own scorer (scripts/ng_str_stutter_score.py):
# its panels, the scorer's classes summed into each, and its outcomes gathered into the error types.
STUTTER_PANELS = [
    ("homopolymers", ["STR segment, motif length 1"]),
    ("motifs of 2 to 6 bases", [f"STR segment, motif length {n}" for n in range(2, 7)]),
    ("STR bundles", ["STR bundle"]),
]
STUTTER_ERRORS = {
    "wrong_gt_truth_ref": ("false_call",),
    "wrong_gt_truth_var": ("wrong", "length_right", "called_reference"),
    "no_call": ("no_call",),
}


def stutter_counts(rows, dataset, samples, classes):
    acc = collections.defaultdict(lambda: collections.defaultdict(float))
    for r in rows:
        if (r["benchmark"] == "stutter_sim" and r["dataset"] == dataset
                and r["samples"] == samples and r["variant_class"] in classes
                and r["caller"] in STYLE):
            acc[(r["depth"], r["caller"])][r["metric"]] += float(r["value"])
    return {key: {error: sum(m.get(o, 0) for o in outcomes)
                  for error, outcomes in STUTTER_ERRORS.items()}
            for key, m in acc.items()}


def stutter_figures(rows, out):
    made = []
    configs = sorted({(r["dataset"], r["samples"]) for r in rows
                      if r["benchmark"] == "stutter_sim" and r["metric"] == "variant_pairs"},
                     key=lambda c: (c[0], int(c[1])))
    for dataset, samples in configs:
        panels = [(panel, stutter_counts(rows, dataset, samples, classes))
                  for panel, classes in STUTTER_PANELS]
        design = dataset.split(", ")[-1]
        plural = "sample" if samples == "1" else "samples"
        slug = design.replace(" stutter", "").replace("-", "_").replace(" ", "_")
        made.append(error_figure(
            panels, f"Errors — simulated STRs, {design}, {samples} {plural}",
            "Drawn 4 Mb genome with exact truth; the unit is one sample at one STR; "
            "HipSTR is not given the bundles", out, f"errors_stutter_{slug}_{samples}s"))
    return made


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
PERF_STYLE["pvc, one step"] = dict(STYLE["pvc (defaults)"], label="pvc, one step")
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
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}" if is_two_or_five(v) else ""))
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
        vcfeval_benchmark_figure(rows, "giab_chr20", "Errors — HG002, chromosome 20",
                                 "GIAB v4.2.1 confident regions; scored by pvc_eval; depths from "
                                 "5× to ~60×", out, "errors_giab_chr20"),
        vcfeval_benchmark_figure(rows, "giab_panel", "Errors — GIAB HG002, HG003, HG004 panel",
                                 "100 confident regions per sample, each sample called alone, "
                                 "summed over the three; scored by pvc_eval", out,
                                 "errors_giab_panel", aggregate=True),
        vcfeval_benchmark_figure(rows, "giab_trio", "Errors — GIAB trio, called jointly",
                                 "HG002, HG003 and HG004 called together over 100 common regions, "
                                 "summed over the three; scored by pvc_eval", out,
                                 "errors_giab_trio", aggregate=True),
        vcfeval_benchmark_figure(rows, "hg002_5mb", "Errors — HG002 over 5 Mb",
                                 "1,000 confident intervals (5.08 Mb); scored by pvc_eval", out,
                                 "errors_hg002_5mb"),
        vcfeval_benchmark_figure(rows, "hg002_tr", "Errors — HG002, GIAB tandem-repeat regions",
                                 "GIAB TR benchmark v1.0.1, 50,000 regions (6.1 Mb); scored by "
                                 "pvc_eval; HipSTR left out: it calls repeat lengths only, and "
                                 "its long records make pvc_eval class its loci as indels and "
                                 "complex, not STRs", out, "errors_hg002_tr",
                                 exclude=("HipSTR",)),
        *stutter_figures(rows, out),
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
        if m:
            print(m)


if __name__ == "__main__":
    main()
