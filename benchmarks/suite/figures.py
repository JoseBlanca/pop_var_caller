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
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter  # noqa: E402

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
DEPTH_ORDER = ["1x", "2x", "3x", "5x", "10x", "15x", "20x", "30x", "50x", "60x", "full", "300x"]


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
#
# **pvc is drawn once, as fitted, named plain "pvc"** (owner, 25 October 2026): its run with the
# defaults is usually a little worse and two pvc rows crowded every depth group. The table keeps
# both.
ERROR_FIGURE_CALLERS = [c for c in CALLERS if c != "pvc (defaults)"]
ERROR_FIGURE_LABEL = {"pvc (fitted)": "pvc"}


def error_label(caller):
    return ERROR_FIGURE_LABEL.get(caller, STYLE[caller]["label"])


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
                and r["variant_class"] == kind and r["caller"] in ERROR_FIGURE_CALLERS
                and r["caller"] not in exclude and r["metric"] in wanted):
            acc[(r["depth"], r["caller"])][r["metric"]] += float(r["value"])
    if not any(v for counts in acc.values() for v in counts.values()):
        return {}
    return {k: {m: v for m, v in c.items() if m != "right"} for k, c in acc.items()}


def row_layout(panels, callers=None):
    """One row position for every (depth, caller) any panel has, shared by every panel, so a row
    reads across them even where a panel lacks a caller (HipSTR and the STR bundles)."""
    callers = callers or ERROR_FIGURE_CALLERS
    keys = set().union(*(counts.keys() for _, counts in panels))
    depths = sorted({d for d, _ in keys}, key=depth_key)
    y, place, groups = 0, {}, []
    for depth in depths:
        top = y
        for caller in callers:
            if (depth, caller) in keys:
                place[(depth, caller)] = y
                y += 1
        groups.append((depth_label(depth), top))
        y += 0.8
    return place, groups, y - 0.8


def dot_panel(ax, counts, layout, show_labels, label=error_label):
    place, groups, bottom = layout
    for key, y in place.items():
        if key not in counts:
            continue
        # pvc's other runs (the inbreeding benchmark's stated F) share pvc's blue.
        colour = STYLE.get(key[1], STYLE["pvc (fitted)"])["color"]
        values = [counts[key].get(metric, 0) for metric, _, _ in ERROR_MARKS]
        ax.plot([min(values), max(values)], [y, y], color=colour, linewidth=1.2, alpha=0.6,
                zorder=2)
        for (metric, _, mark), value in zip(ERROR_MARKS, values):
            ax.plot(value, y, linestyle="none", marker=mark["marker"], markersize=mark["size"],
                    markerfacecolor=colour if mark["fill"] else "none",
                    markeredgecolor=SURFACE if mark["fill"] else colour,
                    markeredgewidth=mark["edge"], alpha=mark["alpha"], zorder=mark["z"])
    ax.set_facecolor(SURFACE)
    # Every axis starts at 0 (owner, 7 November 2026). A zero has no place on a log axis: a short
    # linear stretch from 0 to 1 holds it.
    ax.set_xscale("symlog", linthresh=1, linscale=0.4)
    ax.set_xlim(left=-0.3)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else "0"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(bottom + 0.5, -0.7)
    ax.set_yticks(list(place.values()))
    ax.set_yticklabels([label(caller) for _, caller in place])
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


def error_figure(panels, title, subtitle, out, name, ncols=None, callers=None, label=error_label):
    """`panels`: [(panel title, {(depth, caller): {error type: count}})]; empty panels dropped.

    The panels sit in one row unless `ncols` says how many go across. `callers` orders the rows
    and `label` names them; both default to the accuracy figures' callers."""
    panels = [(t, c) for t, c in panels if c]
    if not panels:
        return None
    layout = row_layout(panels, callers)
    ncols = min(ncols or len(panels), len(panels))
    nrows = math.ceil(len(panels) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols + 1.6,
                                                    (1.4 + 0.27 * len(layout[0])) * nrows),
                             facecolor=SURFACE, sharey=True, squeeze=False)
    for i, ((panel, counts), ax) in enumerate(zip(panels, axes.flat)):
        dot_panel(ax, counts, layout, show_labels=(i % ncols == 0), label=label)
        ax.set_title(panel, loc="left", color=TEXT, fontsize=11)
        ax.set_xlabel("errors (log scale)", color=TEXT_2, fontsize=9)
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    fig.tight_layout(rect=(0.04, 0.05, 1, 1))
    handles = [Line2D([], [], linestyle="none", marker=m["marker"], markersize=m["size"],
                      markeredgewidth=m["edge"], alpha=m["alpha"],
                      markerfacecolor=TEXT_2 if m["fill"] else "none",
                      markeredgecolor=SURFACE if m["fill"] else TEXT_2, label=label)
               for _, label, m in ERROR_MARKS]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.01), labelcolor=TEXT)
    # The heading's lines, stacked upwards from the top of the panels and wrapped to the figure's
    # width (a character is about half its font size wide), so a narrow figure is not widened.
    pt = 1 / 72 / fig.get_figheight()
    y = 4
    for text, size, colour in ((ERROR_SUBTITLE, 8, TEXT_2), (subtitle, 9, TEXT_2),
                               (title, 13, TEXT)):
        wrapped = textwrap.fill(text, int(fig.get_figwidth() * 72 / (0.5 * size)))
        fig.text(0.01, 1 + y * pt, wrapped, ha="left", va="bottom", color=colour, fontsize=size,
                 linespacing=1.3)
        y += (wrapped.count("\n") + 1) * size * 1.3 + 2
    os.makedirs(out, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), facecolor=SURFACE, dpi=160,
                    bbox_inches="tight")
    plt.close(fig)
    return os.path.join(out, f"{name}.svg")


# pvc_eval's classes of variant, one panel each, in two figures of two panels: each figure is one
# column, so a panel is as wide as the page it is shown on.
KIND_FIGURES = [
    ("snps_indels", "SNPs and indels", [("SNPs", "SNP"), ("indels", "indel")]),
    ("strs_complex", "STRs and complex variants", [("STRs", "STR"), ("complex", "complex")]),
]


def vcfeval_benchmark_figure(rows, benchmark, title, subtitle, out, name, aggregate=False,
                             exclude=()):
    """Two figures, `<name>_snps_indels` and `<name>_strs_complex`; returns their paths."""
    made = []
    for suffix, kinds_title, kinds in KIND_FIGURES:
        panels = [(title_, error_counts(rows, benchmark, kind, aggregate, exclude))
                  for title_, kind in kinds]
        made.append(error_figure(panels, f"{title}: {kinds_title}", subtitle, out,
                                 f"{name}_{suffix}", ncols=1))
    return made


# The refinement benchmark (bench/refine_sim.py): one panel per cohort, at the depth where the
# refinement matters most, genotypes with and without reads summed. SNPs only, because STITCH calls
# only SNPs; the figure with every class of variant and no STITCH told the same story (owner,
# 6 November 2026).
# GATK alone fills nothing and is in the page's tables, not the figure, which it crowded (owner,
# 7 November 2026). GATK + Beagle, which refines GATK's calls, takes GATK's colour and mark, and
# STITCH the light neutral HipSTR has in the repeat benchmarks, where STITCH never appears.
REFINE_PANELS = [
    # (panel, dataset, {the benchmark's depth group: the depth it is drawn as})
    ("F2, 20 plants at 1×", "simulated F2", {"F2 1x": "1x"}),
    ("outbred, 100 samples a depth", "simulated natural population, outbred",
     {d: d for d in ("2x", "5x", "10x", "30x")}),
    ("selfing F = 0.9, 100 samples a depth", "simulated natural population, selfing F = 0.9",
     {d: d for d in ("2x", "5x", "10x", "30x")}),
]
REFINE_STYLE = {"pvc (refinement on)": dict(STYLE["pvc (fitted)"], label="pvc"),
                "GATK + Beagle 4.1": dict(STYLE["GATK"], label="GATK + Beagle"),
                "STITCH": dict(STYLE["HipSTR"], label="STITCH")}
REFINE_CALLERS = ["pvc (refinement on)", "pvc (refinement off)", "GATK + Beagle 4.1", "STITCH"]


def refine_figures(rows, out):
    """`errors_refine_sim_snps`: false positives against false negatives in each cohort, SNPs only
    (STITCH calls only SNPs), genotypes with and without reads summed, the genotypes left missing
    in a panel of their own under each."""
    wanted = {metric for metric, _, _ in ERROR_MARKS}
    classes = {"SNP, with reads", "SNP, without reads"}
    panels = []
    for panel, dataset, groups in REFINE_PANELS:
        acc = collections.defaultdict(lambda: collections.defaultdict(float))
        for r in rows:
            if (r["benchmark"] == "refine_sim" and r["dataset"] == dataset
                    and r["depth"] in groups and r["stage"] == "error type"
                    and r["variant_class"] in classes and r["caller"] in REFINE_CALLERS):
                # Every scored point is kept, a caller with no error at all included.
                counts = acc[(groups[r["depth"]], r["caller"])]
                if r["metric"] in wanted:
                    counts[r["metric"]] += float(r["value"])
        # A count the table has no row for is a zero (Beagle leaves no genotype missing).
        panels.append((panel, {k: {m: v.get(m, 0.0) for m in wanted} for k, v in acc.items()}))
    if not all(c for _, c in panels):
        return []
    # pvc with the refinement off, the baseline the refinement is judged against (owner,
    # 7 November 2026): pvc's blue, paler, and a mark of its own.
    style = dict(REFINE_STYLE, **{"pvc (refinement off)": dict(
        STYLE["pvc (fitted)"], color=tint(STYLE["pvc (fitted)"]["color"], 0.45), marker="o",
        label="pvc, refinement off")})
    saved = {c: STYLE.get(c) for c in style}
    STYLE.update(style)  # scatter_figure styles its callers by STYLE
    try:
        return [scatter_figure(panels, "Errors — refined genotypes, simulated cohorts: SNPs", None,
                               out, "errors_refine_sim_snps", callers=REFINE_CALLERS,
                               missing_panel=True)]
    finally:
        for c, st in saved.items():
            if st is None:
                STYLE.pop(c, None)
            else:
                STYLE[c] = st


# pvc with the refinement off and on, each kind of error in a panel of its own (owner, 7 November
# 2026): the counts run from tens to tens of thousands, so each panel has its own linear scale from
# 0, and the cohorts sit side by side within it.
# The lowest depth of each cohort only, where the refinement matters most (owner, 7 November 2026).
REFINE_BAR_GROUPS = [
    ("F2, 1×", "simulated F2", "F2 1x"),
    ("outbred, 2×", "simulated natural population, outbred", "2x"),
    ("selfing, 2×", "simulated natural population, selfing F = 0.9", "2x"),
]
REFINE_BAR_PANELS = [
    ("wrong_gt_truth_var", "missed: the truth has a variant and the call is not its genotype"),
    ("wrong_gt_truth_ref", "false: a variant called where the truth is reference"),
    ("no_call", "no call: the genotype left missing (./.)"),
]


def refine_on_off_figure(rows, out):
    """`errors_refine_on_off`: pvc's SNP errors with the refinement off and on, genotypes with and
    without reads summed; one panel per kind of error, two bars per cohort and depth."""
    callers = [("pvc (refinement off)", "refinement off"), ("pvc (refinement on)", "refinement on")]
    counts = collections.defaultdict(float)
    for r in rows:
        if (r["benchmark"] == "refine_sim" and r["stage"] == "error type"
                and r["variant_class"] in ("SNP, with reads", "SNP, without reads")
                and r["caller"] in dict(callers) and r["value"] != ""):
            counts[(r["dataset"], r["depth"], r["caller"], r["metric"])] += float(r["value"])
    if not counts:
        return None
    blue = STYLE["pvc (fitted)"]["color"]
    colours = {"pvc (refinement off)": tint(blue, 0.45), "pvc (refinement on)": blue}
    fig, axes = plt.subplots(len(REFINE_BAR_PANELS), 1, figsize=(5.2, 2.8 * len(REFINE_BAR_PANELS)),
                             facecolor=SURFACE)
    width = 0.38
    for ax, (metric, panel) in zip(axes, REFINE_BAR_PANELS):
        ax.set_facecolor(SURFACE)
        top = max(counts[(ds, d, c, metric)] for _, ds, d in REFINE_BAR_GROUPS for c, _ in callers)
        for i, (caller, _) in enumerate(callers):
            xs = [g + (i - 0.5) * (width + 0.02) for g in range(len(REFINE_BAR_GROUPS))]
            ys = [counts[(ds, d, caller, metric)] for _, ds, d in REFINE_BAR_GROUPS]
            ax.bar(xs, ys, width=width, color=colours[caller], zorder=2)
            for x, y in zip(xs, ys):
                ax.annotate(f"{y:,.0f}", (x, y), xytext=(0, 2), textcoords="offset points",
                            ha="center", va="bottom", fontsize=7, color=TEXT_2)
        ax.set_ylim(0, top * 1.15)
        ax.set_xticks(range(len(REFINE_BAR_GROUPS)))
        ax.set_xticklabels([g for g, _, _ in REFINE_BAR_GROUPS], fontsize=9, color=TEXT)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.grid(True, axis="y", color=GRID, linewidth=0.8, zorder=0)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(TEXT_2)
        ax.tick_params(colors=TEXT_2, labelsize=9)
        ax.set_title(panel, color=TEXT, fontsize=10, loc="left")
        ax.set_ylabel("SNP genotypes", color=TEXT_2, fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    handles = [Line2D([], [], linestyle="none", marker="s", markersize=10,
                      markerfacecolor=colours[c], markeredgecolor=colours[c], label=label)
               for c, label in callers]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0), labelcolor=TEXT)
    fig.text(0.01, 0.985, "pvc, refinement off and on — SNPs",
             ha="left", va="top", color=TEXT, fontsize=13)
    return save(fig, out, "errors_refine_on_off")


# The simulated inbreeding benchmark (bench/inbreeding_sim.py) counts each sample's genotype against
# the truth by true and called genotype; these sums gather them into the error types. One figure per
# group of samples with the same F, at one depth. pvc is drawn only with F estimated (owner,
# 1 November 2026); its runs with a wrong F are in the page's tables.
INBREEDING_CALLERS = ["pvc (fitted)", "GATK", "DeepVariant"]
INBREEDING_LABEL = {"pvc (fitted)": "pvc, F estimated"}
INBREEDING_ERRORS = {
    # A planted SNP's homozygous reference called a variant, and a variant called where nothing
    # was planted (the scorer's false_variant_genotypes).
    "wrong_gt_truth_ref": tuple(f"true hom_ref, called {c}" for c in ("het", "hom_alt",
                                                                      "other_allele"))
                          + ("false_variant_genotypes",),
    "wrong_gt_truth_var": tuple(f"true het, called {c}" for c in ("hom_ref", "hom_alt",
                                                                  "other_allele"))
                          + tuple(f"true hom_alt, called {c}" for c in ("hom_ref", "het",
                                                                        "other_allele")),
    "no_call": tuple(f"true {t}, called no_call" for t in ("hom_ref", "het", "hom_alt")),
}
INBREEDING_FIGURES = [("F = 1.00", "1", "f1"), ("F = 0.75", "0.75", "f075")]


def inbreeding_figures(rows, out, depth="15x"):
    made = []
    for group, f, slug in INBREEDING_FIGURES:
        acc = collections.defaultdict(lambda: collections.defaultdict(float))
        for r in rows:
            if (r["benchmark"] == "inbreeding_sim" and r["depth"] == depth
                    and r["variant_class"] == group and r["caller"] in INBREEDING_CALLERS
                    and r["value"] != ""):
                acc[(r["depth"], r["caller"])][r["metric"]] += float(r["value"])
        counts = {key: {error: sum(m.get(o, 0) for o in outcomes)
                        for error, outcomes in INBREEDING_ERRORS.items()}
                  for key, m in acc.items()}
        made.append(error_figure(
            [("SNPs", counts)], f"Errors — simulated inbreeding, the 10 samples at F = {f}",
            f"50 simulated tomato samples called together, 10 at each F from 0 to 1; "
            f"{depth_label(depth)}; summed over the 10", out, f"errors_inbreeding_{slug}",
            callers=INBREEDING_CALLERS, label=lambda c: INBREEDING_LABEL[c] if c in INBREEDING_LABEL else error_label(c)))
    return made


# The scatter's two kinds of error, from the dot plot's three: a genotype called wrong (either way)
# and a genotype left missing.
INBREEDING_KINDS = [
    ("wrong", "wrong genotype", INBREEDING_ERRORS["wrong_gt_truth_ref"]
     + INBREEDING_ERRORS["wrong_gt_truth_var"]),
    ("missing", "missing genotype", INBREEDING_ERRORS["no_call"]),
]


def tint(colour, share=0.5):
    """`colour` mixed with the page's surface: the same hue, subdued."""
    rgb = matplotlib.colors.to_rgb(colour)
    base = matplotlib.colors.to_rgb(SURFACE)
    return tuple(c * (1 - share) + b * share for c, b in zip(rgb, base))


def inbreeding_totals(rows, depth):
    """{(caller, kind): [(F, count over the group's 10 samples)]}, one point per group."""
    groups = sorted({r["variant_class"] for r in rows if r["benchmark"] == "inbreeding_sim"
                     and r["variant_class"].startswith("F = ")})
    series = {}
    for group in groups:
        for caller in INBREEDING_CALLERS:
            m = {r["metric"]: float(r["value"]) for r in rows
                 if r["benchmark"] == "inbreeding_sim" and r["depth"] == depth
                 and r["variant_class"] == group and r["caller"] == caller and r["value"] != ""}
            if not m:
                continue
            for kind, _, outcomes in INBREEDING_KINDS:
                series.setdefault((caller, kind), []).append(
                    (float(group.split("=")[1]), sum(m.get(o, 0) for o in outcomes)))
    return series


def inbreeding_scatter_figure(rows, out, depths=("15x", "30x")):
    """Errors in a group of samples against the group's F, one panel per depth on its own scale:
    for each caller a line of wrong genotypes and a paler, dashed line of missing ones."""
    panels = [(d, inbreeding_totals(rows, d)) for d in depths]
    panels = [(d, s) for d, s in panels if s]
    if not panels:
        return None
    fig, axes = plt.subplots(len(panels), 1, figsize=(6.4, 4.6 * len(panels)), facecolor=SURFACE,
                             squeeze=False)
    for (depth, series), ax in zip(panels, axes[:, 0]):
        style_axes(ax, "inbreeding coefficient F of the group's 10 samples",
                   "errors in the group's 10 samples\n(log scale)", depth_label(depth))
        ax.set_xscale("linear")
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_xlim(-0.05, 1.32)
        counts = [y for pts in series.values() for _, y in pts]
        if min(counts) == 0:
            # A zero has no place on a log axis: a short linear stretch from 0 to 1 holds it.
            ax.set_yscale("symlog", linthresh=1, linscale=0.4)
            low = -0.3
        else:
            ax.set_yscale("log")
            low = min(counts) / 1.3
        high = max(counts) * 1.3
        ax.set_ylim(low, high)
        ax.set_yticks([t for t in (0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000,
                                   20000) if low <= t <= high])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        for caller in INBREEDING_CALLERS:
            st = STYLE[caller]
            label = SCATTER_LABEL.get(caller, st["label"])
            for kind, kind_label, _ in INBREEDING_KINDS:
                pts = series.get((caller, kind))
                if not pts:
                    continue
                missing = kind == "missing"
                colour = tint(st["color"]) if missing else st["color"]
                xs, ys = [p[0] for p in pts], [p[1] for p in pts]
                ax.plot(xs, ys, color=colour, linestyle="--" if missing else "-", linewidth=2,
                        marker=st["marker"], markersize=7,
                        markerfacecolor=SURFACE if missing else colour,
                        markeredgecolor=colour if missing else SURFACE, markeredgewidth=1.5,
                        zorder=3)
                ax.annotate(f"{label}, {'missing' if missing else 'wrong'}", (xs[-1], ys[-1]),
                            xytext=(6, 0), textcoords="offset points", fontsize=8.5,
                            color=TEXT if not missing else TEXT_2, va="center", zorder=4,
                            gid="caller")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    handles = [Line2D([], [], color=TEXT_2, linestyle="-", linewidth=2, marker="o",
                      markerfacecolor=TEXT_2, markeredgecolor=SURFACE, label="wrong genotype"),
               Line2D([], [], color=tint(TEXT_2), linestyle="--", linewidth=2, marker="o",
                      markerfacecolor=SURFACE, markeredgecolor=tint(TEXT_2),
                      label="missing genotype (./.)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.01), labelcolor=TEXT)
    pt = 1 / 72 / fig.get_figheight()
    y = 4
    for text, size, colour in (
            ("50 samples called together, 10 at each F; each point counts one group's 10 samples' "
             "genotypes at the planted SNPs, plus variants called where none was planted; pvc "
             "estimates each sample's F, GATK and DeepVariant take none", 9, TEXT_2),
            ("Errors against inbreeding — simulated tomato samples", 13, TEXT)):
        wrapped = textwrap.fill(text, int(fig.get_figwidth() * 72 / (0.5 * size)))
        fig.text(0.01, 1 + y * pt, wrapped, ha="left", va="bottom", color=colour, fontsize=size,
                 linespacing=1.3)
        y += (wrapped.count("\n") + 1) * size * 1.3 + 2
    return save(fig, out, "errors_inbreeding_by_f")


# ---- accuracy: false positives against false negatives ---------------------------------------
#
# For a benchmark where no caller leaves a genotype missing (one sample called alone: every
# caller writes a genotype for each record it emits), the dot plot's third mark is always 0 and
# two numbers are left. A scatter shows them better: false positives across, the variants the
# truth has and the call does not give down, one line per caller through its depths, so the
# corner nearest 0, 0 is the good one. Both axes are logarithmic with a short linear stretch from
# 0 to 1, which holds the zeros.

SCATTER_LABEL = {"pvc (fitted)": "pvc"}
# Marker sizes for the depths, smallest first: depth is the marker's size, so every caller's points
# say their depth without a label each.
SCATTER_SIZES = [5, 8.5, 12, 15.5, 19]


# SNPs, and every other class summed: the split most readers expect (owner, 7 November 2026), used
# by doc/user/accuracy.md's chromosome 20 and trio figures. The STRs alone are on its STR section.
SNPS_OTHER_PANELS = [("SNPs", ("SNP",)), ("indels, STRs and complex variants", ("indel", "STR", "complex"))]


def snps_other_panels(rows, benchmark, aggregate=False):
    """[(panel title, {(depth, caller): error counts})]: SNPs, then the other classes summed."""
    panels = []
    for panel, kinds in SNPS_OTHER_PANELS:
        summed = collections.defaultdict(lambda: collections.defaultdict(float))
        for kind in kinds:
            for key, counts in error_counts(rows, benchmark, kind, aggregate).items():
                for metric, value in counts.items():
                    summed[key][metric] += value
        if summed:
            panels.append((panel, {k: dict(c) for k, c in summed.items()}))
    return panels


def snps_other_scatter_figure(rows, benchmark, title, out, name, aggregate=False,
                              missing_panel=False, from_data=()):
    """`<name>_scatter_snps_other`: false positives against false negatives, two panels; the
    genotypes left missing, of both kinds together, in one panel under them (owner, 8 November
    2026)."""
    return scatter_figure(snps_other_panels(rows, benchmark, aggregate), title, None, out,
                          f"{name}_scatter_snps_other", missing_panel=missing_panel,
                          sum_missing=True, from_data=from_data)


def zero_log_axes(ax, xs, ys):
    """Both axes logarithmic from 0, through a short linear stretch from 0 to 1 that holds the
    zeros (owner, 7 November 2026); `xs`, `ys` the values each axis must reach."""
    for setter, limits, values in ((ax.set_xscale, ax.set_xlim, xs), (ax.set_yscale, ax.set_ylim, ys)):
        if values is None:
            continue
        setter("symlog", linthresh=1, linscale=0.4)
        limits(-0.3, max(max(values), 1) * 1.6)
    fmt = FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else "0")
    minor = FuncFormatter(lambda v, _: f"{v:,.0f}" if is_two_or_five(v) else "")
    for axis, values in ((ax.xaxis, xs), (ax.yaxis, ys)):
        if values is None:
            continue
        low, high = axis.get_view_interval()
        axis.set_major_formatter(fmt)
        # Over more than three decades the 2s and 5s crowd the decades' own labels.
        crowded = high / max(low, 1) > 1e3
        if not crowded:
            # The zero stretch's axis has no ticks between decades of its own.
            axis.set_minor_locator(FixedLocator(
                [m * 10 ** k for k in range(5) for m in (2, 5) if m * 10 ** k < high]))
        axis.set_minor_formatter(NullFormatter() if crowded else minor)


def scatter_figure(panels, title, subtitle, out, name, callers=ERROR_FIGURE_CALLERS,
                   shared_axes=False, missing_panel=False, sum_missing=False, from_data=()):
    """One panel a row, false positives against false negatives, one line a caller through its
    depths, a marker's size its depth. `panels`: (panel title, {(depth, caller): error counts}).
    With `shared_axes`, every panel has the same scales, so the panels compare at a glance. With
    `missing_panel`, and some genotype left missing, each panel has a third as tall a panel below
    it: the genotypes left missing (./.) at each depth, in the same lines and marks (owner,
    8 November 2026). With `sum_missing` as well, there is one such panel, under them all, its
    counts the panels' summed. Every axis starts at 0 except those of the panels named in
    `from_data`."""
    styles = {c: dict(STYLE[c], label=SCATTER_LABEL.get(c, STYLE[c]["label"]),
                      linestyle="-") for c in callers}
    depths = sorted({d for _, c in panels for d, _ in c}, key=depth_key)
    size = {d: SCATTER_SIZES[i] for i, d in enumerate(depths)}
    missing = [v.get("no_call", 0) for _, cs in panels for v in cs.values()]
    missing_panel = missing_panel and any(missing)
    if missing_panel and sum_missing:
        fig, grid = plt.subplots(len(panels) + 1, 1, figsize=(6.4, 5.6 * len(panels) + 2.2),
                                 facecolor=SURFACE, squeeze=False,
                                 gridspec_kw={"height_ratios": [3] * len(panels) + [1]})
        axes, missing_axes = grid[:-1, 0], [None] * len(panels)
    elif missing_panel:
        fig, grid = plt.subplots(2 * len(panels), 1, figsize=(6.4, 7.6 * len(panels)),
                                 facecolor=SURFACE, squeeze=False,
                                 gridspec_kw={"height_ratios": [3, 1] * len(panels)})
        axes, missing_axes = grid[0::2, 0], grid[1::2, 0]
    else:
        fig, grid = plt.subplots(len(panels), 1, figsize=(6.4, 5.6 * len(panels)),
                                 facecolor=SURFACE, squeeze=False)
        axes, missing_axes = grid[:, 0], [None] * len(panels)
    for (panel, counts), ax, mx in zip(panels, axes, missing_axes):
        series = {}
        for depth, caller in sorted(counts, key=lambda k: depth_key(k[0])):
            c = counts[(depth, caller)]
            series.setdefault(caller, []).append(
                (c.get("wrong_gt_truth_ref", 0), c.get("wrong_gt_truth_var", 0),
                 depth))
        style_axes(ax, "false positives: a variant called where the truth is reference",
                   "false negatives: the truth has a variant and the call\n"
                   "is not its genotype", panel)
        for caller in callers:
            pts = series.get(caller)
            if not pts:
                continue
            st = styles[caller]
            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            ax.plot(xs, ys, color=st["color"], linewidth=1.6, alpha=0.7, zorder=2)
            # Hollow, so markers that overlap where the callers crowd still show.
            for x, y, depth in pts:
                ax.plot(x, y, linestyle="none", marker=st["marker"],
                        markersize=size[depth], markerfacecolor="none",
                        markeredgecolor=st["color"], markeredgewidth=1.8, zorder=3)
            ax.annotate(st["label"], (xs[-1], ys[-1]), xytext=(10, 0),
                        textcoords="offset points", fontsize=8.5, color=TEXT, va="center",
                        zorder=4, gid="caller")
        points = [p for pts in series.values() for p in pts]
        if shared_axes:
            points = [(c.get("wrong_gt_truth_ref", 0), c.get("wrong_gt_truth_var", 0))
                      for _, cs in panels for c in cs.values()]
        if panel in from_data:
            # The panels named in `from_data` are on log axes from just under their data, not from
            # 0 (owner, 8 November 2026: the chromosome 20 SNPs, whose counts are all far above 0).
            for limits, values in ((ax.set_xlim, [p[0] for p in points]),
                                   (ax.set_ylim, [p[1] for p in points])):
                limits(min(values) / 1.6, max(values) * 1.6)
        else:
            zero_log_axes(ax, [p[0] for p in points], [p[1] for p in points])
        if mx is not None:
            missing_by_depth(mx, counts, callers, styles,
                             sorted({d for d, _ in counts}, key=depth_key),
                             missing if shared_axes
                             else [v.get("no_call", 0) for v in counts.values()])
    if missing_panel and sum_missing:
        summed = collections.defaultdict(lambda: {"no_call": 0})
        for _, counts in panels:
            for key, c in counts.items():
                summed[key]["no_call"] += c.get("no_call", 0)
        missing_by_depth(grid[-1, 0], summed, callers, styles, depths,
                         [c["no_call"] for c in summed.values()],
                         "genotypes left missing, all kinds of variant")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=SCATTER_SIZES[i],
                      markerfacecolor="none", markeredgecolor=TEXT_2, markeredgewidth=1.8,
                      label=depth_label(d)) for i, d in enumerate(depths)]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.01), labelcolor=TEXT,
               title="read depth", title_fontsize=9)
    pt = 1 / 72 / fig.get_figheight()
    y = 4
    for text, size, colour in ((subtitle, 9, TEXT_2), (title, 13, TEXT)):
        if not text:
            continue
        wrapped = textwrap.fill(text, int(fig.get_figwidth() * 72 / (0.5 * size)))
        fig.text(0.01, 1 + y * pt, wrapped, ha="left", va="bottom", color=colour,
                 fontsize=size, linespacing=1.3)
        y += (wrapped.count("\n") + 1) * size * 1.3 + 2
    return save(fig, out, name)


def missing_by_depth(ax, counts, callers, styles, depths, reach, title=""):
    """The genotypes each caller left missing (./.) against read depth: the panel's depths evenly
    spaced across, the counts up from 0 to `reach`'s most."""
    drawn = [c for c in callers if any(k[1] == c for k in counts)]
    if len(depths) == 1:
        # One depth has no line to draw: a column a caller, in its colour (owner, 8 November 2026).
        style_axes(ax, f"caller, at {depth_label(depths[0])}", "genotypes left\nmissing (./.)",
                   title)
        ax.set_xscale("linear")
        for i, caller in enumerate(drawn):
            value = counts[(depths[0], caller)].get("no_call", 0)
            ax.bar(i, value, width=0.6, color=styles[caller]["color"], zorder=2)
            ax.annotate(f"{value:,.0f}", (i, value), xytext=(0, 2), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7.5, color=TEXT_2)
        ax.set_xticks(range(len(drawn)))
        ax.set_xticklabels([styles[c]["label"] for c in drawn], fontsize=8)
        ax.set_xlim(-0.6, len(drawn) - 0.4)
        ax.grid(False, axis="x")
        zero_log_axes(ax, None, [v * 2 for v in reach])
        return
    style_axes(ax, "read depth", "genotypes left\nmissing (./.)", title)
    ax.set_xscale("linear")
    for caller in drawn:
        # Each caller a little to the side, so callers with the same count, often none, all show.
        shift = 0.07 * (drawn.index(caller) - (len(drawn) - 1) / 2)
        pts = [(depths.index(d) + shift, counts[(d, c)].get("no_call", 0))
               for d, c in sorted(counts, key=lambda k: depth_key(k[0])) if c == caller]
        if not pts:
            continue
        st = styles[caller]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=st["color"], linewidth=1.6,
                alpha=0.7, marker=st["marker"], markersize=7, markerfacecolor="none",
                markeredgecolor=st["color"], markeredgewidth=1.6, zorder=3)
    ax.set_xticks(range(len(depths)))
    ax.set_xticklabels([depth_label(d) for d in depths])
    ax.set_xlim(-0.4, len(depths) - 0.6)
    zero_log_axes(ax, None, list(reach))


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
                and r["caller"] in ERROR_FIGURE_CALLERS):
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


# The homopolymer figure: the 50-sample cohort, one panel per library, drawn like the chromosome 20
# scatter, with the genotypes left missing in a panel of their own under each (owner, 7 and
# 8 November 2026): folded into the false negatives, DeepVariant's empty genotypes at reference homopolymers
# read as missed variants, and the trio's dot plot drew a caller with none as a long line.
HOMOPOLYMER_CALLERS = ["pvc (fitted)", "GATK", "freebayes", "DeepVariant", "HipSTR"]


def stutter_homopolymer_figure(rows, out, samples="50"):
    panels = []
    for title, dataset in (("PCR stutter", "simulated 4 Mb, pcr stutter"),
                           ("PCR-free", "simulated 4 Mb, pcr-free stutter")):
        counts = stutter_counts(rows, dataset, samples, ["STR segment, motif length 1"])
        panels.append((title, {key: c for key, c in counts.items() if key[1] in HOMOPOLYMER_CALLERS}))
    if not all(c for _, c in panels):
        return None
    return scatter_figure(panels, f"Errors — simulated homopolymers, {samples} samples", None,
                          out, f"errors_stutter_homopolymers_{samples}s",
                          callers=HOMOPOLYMER_CALLERS, shared_axes=True, missing_panel=True)


# The performance figures' lines: (line, caller in the table, its stages in order of preference).
# pvc is three lines, its three ways of spending time: in one step from the reads, building the
# psps from the reads, and calling from the psps. Its parameter fit is optional, made once for a
# cohort, and is in no line (the summary reports it on its own).
PERF_LINES = [
    ("pvc, one step", "pvc (defaults)", ("call-from-alignments",)),
    ("pvc, one step, 3 segments", "pvc (defaults, 3 segments)", ("whole pipeline",)),
    ("pvc, building psps", "pvc (fitted)", ("generate-psps",)),
    ("pvc, calling from psps", "pvc (fitted)", ("call-from-psps",)),
    ("GATK", "GATK", ("whole pipeline", "HaplotypeCaller")),
    ("freebayes", "freebayes", ("freebayes",)),
    ("DeepVariant", "DeepVariant", ("DeepVariant + GLnexus", "DeepVariant")),
]
PERF_STYLE = dict(STYLE)
PERF_STYLE["pvc, one step"] = dict(STYLE["pvc (defaults)"], label="pvc, one step")
PERF_STYLE["pvc, one step, 3 segments"] = dict(color="#2a78d6", linestyle="-.", marker="*",
                                               label="pvc, one step, 3 segments")
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
            if r.get("outcome"):  # did not finish: a lower bound, not a point on the curve
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
        *vcfeval_benchmark_figure(rows, "giab_chr20", "Errors — HG002, chromosome 20",
                                 "GIAB v4.2.1 confident regions; scored by pvc_eval; depths from "
                                 "5× to ~60×", out, "errors_giab_chr20"),
        # No subtitle: doc/user/accuracy.md explains the figure in its caption.
        snps_other_scatter_figure(rows, "giab_chr20", "Errors — HG002, chromosome 20", out,
                                  "errors_giab_chr20", from_data=("SNPs",)),
        *vcfeval_benchmark_figure(rows, "giab_panel", "Errors — GIAB HG002, HG003, HG004 panel",
                                 "100 confident regions per sample, each sample called alone, "
                                 "summed over the three; scored by pvc_eval", out,
                                 "errors_giab_panel", aggregate=True),
        # The trio's empty genotypes in the marks' size, as the homopolymer figure's.
        snps_other_scatter_figure(rows, "giab_trio", "Errors — GIAB trio, called jointly", out,
                                  "errors_giab_trio", aggregate=True, missing_panel=True),
        *vcfeval_benchmark_figure(rows, "hg002_5mb", "Errors — HG002 over 5 Mb",
                                 "1,000 confident intervals (5.08 Mb); scored by pvc_eval", out,
                                 "errors_hg002_5mb"),
        *vcfeval_benchmark_figure(rows, "hg002_tr", "Errors — HG002, GIAB tandem-repeat regions",
                                 "GIAB TR benchmark v1.0.1, 50,000 regions (6.1 Mb); scored by "
                                 "pvc_eval; HipSTR left out: it calls repeat lengths only, and "
                                 "its long records make pvc_eval class its loci as indels and "
                                 "complex, not STRs", out, "errors_hg002_tr",
                                 exclude=("HipSTR",)),
        *stutter_figures(rows, out),
        *refine_figures(rows, out),
        refine_on_off_figure(rows, out),
        stutter_homopolymer_figure(rows, out),
        *inbreeding_figures(rows, out),
        inbreeding_scatter_figure(rows, out),
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
