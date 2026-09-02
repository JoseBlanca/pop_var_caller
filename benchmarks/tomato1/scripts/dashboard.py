# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
#     "matplotlib-venn",
# ]
# ///

# Marimo dashboard comparing the SNP callers — pop_var_caller ("ours"),
# freebayes, GATK, and ng — on the tomato cohort test set. Toggle between the
# single-sample and cohort variants of the test, and pick which callers to
# compare.
#
# Recommended invocation (no global marimo install needed):
#
#   uvx marimo edit --sandbox benchmarks/tomato1/scripts/dashboard.py
#
# To serve as a read-only app instead of opening the editor:
#
#   uvx marimo run --sandbox benchmarks/tomato1/scripts/dashboard.py
#
# There is no truth set for tomato. This dashboard is about **agreement**, not
# accuracy: which calls the callers share, and how their confidences relate.
# Accuracy lives in the GIAB benchmarks.
#
# Measures:
#   0. Inputs — each VCF's size, record count and sample count. **Read the
#      sample counts before anything else**: a VCF built from a different
#      number of accessions is not comparable with the others, and a stale one
#      is the likeliest reason two callers seem to disagree wildly.
#   1. Variant agreement — per-caller totals, pairwise overlap, and the
#      shared/exclusive split, at (CHROM, POS, REF, ALT) granularity. No
#      PASS / QUAL filter applied; multi-allelic records are split per ALT;
#      `*` (gVCF spanning deletion) is skipped. NB: indel calls aren't
#      normalised, so the same indel may show as discordant if callers
#      anchored it on different REF bases. Add `bcftools norm` upstream if
#      that matters for a measure. A Venn is drawn when two or three callers
#      are selected.
#   2. Pairwise QUAL hexbins on shared variants — one 2D density plot per
#      caller pair, log colour scale, y=x reference line; bounded by the
#      QUAL-cap slider.
#   3. Per-caller QUAL distributions — one histogram per selected caller,
#      sharing the same x-axis (0..cap) and bin edges; independent y-axes
#      (callers can differ by an order of magnitude in record count);
#      y-scale linear or log via radio toggle.
#
# One thing to hold against ng's record count: it does not build loci inside
# tandem repeats yet — about 5 bases in every 100 of this benchmark's ground —
# so calls the others make there have no ng counterpart to agree with. Its run
# report (results/ng/cohort.log) says how much ground that was.

import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import gzip
    from pathlib import Path

    import marimo as mo

    return Path, gzip, mo


@app.cell
def _(mo):
    mo.md("""
    # SNP caller comparison — tomato cohort test
    """)
    return


@app.cell
def _(mo):
    mode = mo.ui.radio(
        options=["single", "cohort"],
        value="cohort",
        label="Test set",
    )
    mode
    return (mode,)


@app.cell
def _(Path, mode):
    # benchmarks/tomato1/scripts/dashboard.py -> benchmarks/tomato1/
    test_dir = Path(__file__).resolve().parent.parent
    results = test_dir / "results"
    sample = "SRR17274057"
    if mode.value == "single":
        all_vcf_paths = {
            "ours": results / "ours" / f"single_{sample}.vcf",
            "freebayes": results / "freebayes" / f"single_{sample}.vcf",
            "gatk": results / "gatk" / f"single_{sample}.vcf",
            "ng": results / "ng" / f"single_{sample}.vcf",
        }
    else:
        all_vcf_paths = {
            "ours": results / "ours" / "cohort" / "cohort.vcf",
            "freebayes": results / "freebayes" / "cohort.vcf",
            "gatk": results / "gatk" / "cohort" / "cohort.vcf",
            "ng": results / "ng" / "cohort.vcf",
        }
    present_paths = {n: p for n, p in all_vcf_paths.items() if p.exists()}
    return all_vcf_paths, present_paths


@app.cell
def _(gzip):
    def vcf_shape(path):
        """(record count, sample count) of a VCF, without loading it whole.

        The sample count is what makes two result files comparable or not, and
        it is invisible in the filename — a cohort VCF left over from a run
        with fewer accessions looks exactly like a current one.
        """
        opener = gzip.open if str(path).endswith(".gz") else open
        records = 0
        samples = 0
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = max(len(line.rstrip("\n").split("\t")) - 9, 0)
                    continue
                if line.strip():
                    records += 1
        return records, samples

    return (vcf_shape,)


@app.cell
def _(all_vcf_paths, mo, present_paths, vcf_shape):
    shapes = {n: vcf_shape(p) for n, p in present_paths.items()}
    _absent = [n for n in all_vcf_paths if n not in present_paths]

    _rows = "\n".join(
        f"| {n} | `{present_paths[n].name}` | "
        f"{present_paths[n].stat().st_size / 1e6:.1f} MB | "
        f"{shapes[n][0]:,} | {shapes[n][1]} |"
        for n in present_paths
    )
    _sample_counts = {s for _, s in shapes.values()}
    _warn = ""
    if len(_sample_counts) > 1:
        _warn = (
            "\n\n> **These VCFs do not hold the same accessions.** "
            f"The sample counts present are {sorted(_sample_counts)}. A file "
            "built from fewer accessions carries fewer segregating sites, so "
            "the agreement counts below are measuring the sample sets as much "
            "as the callers. Re-run the short caller on the full set before "
            "reading anything into its exclusive calls."
        )
    _miss = ""
    if _absent:
        _miss = "\n\nNot present (run its `benchmarks/lib/run_*.sh` first): " + ", ".join(
            f"`{n}`" for n in _absent
        )
    inputs_view = mo.md(
        "### 0. Inputs\n\n"
        "| caller | file | size | records | samples |\n|---|---|---:|---:|---:|\n"
        + _rows
        + _warn
        + _miss
    )
    inputs_view
    return inputs_view, shapes


@app.cell
def _(mo, present_paths):
    # Which callers to compare. Defaults to everything present, but the Venn
    # below only draws for two or three — deselect to get one.
    caller_sel = mo.ui.multiselect(
        options=list(present_paths.keys()),
        value=list(present_paths.keys()),
        label="callers to compare",
    )
    caller_sel
    return (caller_sel,)


@app.cell
def _(gzip):
    def variant_keys(path):
        """Read a VCF and return a set of (CHROM, POS, REF, ALT) tuples.

        Multi-allelic records are split into one entry per ALT.
        `*` (gVCF spanning deletion marker) is skipped — it is not a
        real variant.  No PASS / QUAL filter is applied here; callers
        of this helper layer filtering on top if needed.
        """
        opener = gzip.open if str(path).endswith(".gz") else open
        keys: set[tuple[str, int, str, str]] = set()
        with opener(path, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alts = fields[4]
                for alt in alts.split(","):
                    if alt == "*" or alt == ".":
                        continue
                    keys.add((chrom, pos, ref, alt))
        return keys

    def variant_quals(path):
        """Return a dict mapping (CHROM, POS, REF, ALT) -> float QUAL.

        QUAL is a *site*-level statistic, so every ALT of a
        multi-allelic record inherits the same value.  Records with
        QUAL = `.` or an unparsable value are skipped — the hexbin
        plots cannot consume missing/non-numeric coordinates.
        """
        opener = gzip.open if str(path).endswith(".gz") else open
        result: dict[tuple[str, int, str, str], float] = {}
        with opener(path, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                qual_str = fields[5]
                if qual_str in (".", ""):
                    continue
                try:
                    qual = float(qual_str)
                except ValueError:
                    continue
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alts = fields[4]
                for alt in alts.split(","):
                    if alt == "*" or alt == ".":
                        continue
                    result[(chrom, pos, ref, alt)] = qual
        return result

    return variant_keys, variant_quals


@app.cell
def _(caller_sel, present_paths, variant_keys):
    sets = {n: variant_keys(present_paths[n]) for n in caller_sel.value}
    return (sets,)


@app.cell
def _(mo, sets: dict[str, set]):
    # Three things, in the order they answer a reader's questions: how many
    # calls each caller made, how much of each pair is shared, and how the
    # selected callers split into shared / exclusive.
    _names = list(sets)
    if len(_names) < 2:
        counts_view = mo.md("_Select at least two callers._")
    else:
        _totals = "\n".join(f"| {n} | {len(sets[n]):,} |" for n in _names)

        _pairs = []
        for _i, _a in enumerate(_names):
            for _b in _names[_i + 1:]:
                _shared = len(sets[_a] & sets[_b])
                _union = len(sets[_a] | sets[_b])
                _pairs.append(
                    f"| {_a} ∩ {_b} | {_shared:,} | "
                    f"{_shared / len(sets[_a]):.1%} | "
                    f"{_shared / len(sets[_b]):.1%} | "
                    f"{_shared / _union:.1%} |"
                )

        _all_shared = set.intersection(*(sets[n] for n in _names))
        _excl = []
        for _n in _names:
            _others = set.union(*(sets[m] for m in _names if m != _n))
            _excl.append(f"| only {_n} | {len(sets[_n] - _others):,} |")

        counts_view = mo.md(
            "### 1. Variant agreement\n\n"
            "_(key: CHROM, POS, REF, ALT; multi-allelic split per ALT; no "
            "PASS/QUAL filter beyond what each runner already applied)_\n\n"
            "**Calls made**\n\n| caller | variants |\n|---|---:|\n" + _totals + "\n\n"
            "**Pairwise overlap** — the same shared count read as a share of "
            "each caller's own set, then of their union (Jaccard).\n\n"
            "| pair | shared | of left | of right | of union |\n"
            "|---|---:|---:|---:|---:|\n" + "\n".join(_pairs) + "\n\n"
            "**Shared and exclusive across the whole selection**\n\n"
            "| set | variants |\n|---|---:|\n"
            f"| in all {len(_names)} | {len(_all_shared):,} |\n" + "\n".join(_excl)
        )
    counts_view
    return (counts_view,)


@app.cell
def _(mo, sets: dict[str, set]):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    _names = list(sets)
    if len(_names) == 3:
        from matplotlib_venn import venn3

        fig, ax = plt.subplots(figsize=(8, 6))
        v = venn3([sets[n] for n in _names], set_labels=tuple(_names), ax=ax)
        ax.set_title("Variant agreement — (CHROM, POS, REF, ALT)")
        # A colour-keyed legend built from the venn patches, because the set
        # labels matplotlib_venn draws beside the circles are hard to
        # associate at a glance.
        _handles = [
            Patch(
                facecolor=v.get_patch_by_id(pid).get_facecolor(),
                edgecolor="black",
                label=lbl,
            )
            for pid, lbl in zip(("100", "010", "001"), _names)
            if v.get_patch_by_id(pid) is not None
        ]
        ax.legend(handles=_handles, loc="upper left", bbox_to_anchor=(0.0, 1.0))
        venn_view = fig
    elif len(_names) == 2:
        from matplotlib_venn import venn2

        fig, ax = plt.subplots(figsize=(7, 5))
        venn2([sets[n] for n in _names], set_labels=tuple(_names), ax=ax)
        ax.set_title("Variant agreement — (CHROM, POS, REF, ALT)")
        venn_view = fig
    else:
        venn_view = mo.md(
            "_A Venn is drawn for two or three callers; the table above covers "
            "any number._"
        )
    venn_view
    return (plt, venn_view)


@app.cell
def _(caller_sel, present_paths, variant_quals):
    # Per-caller QUAL maps. Same shape as `sets` upstream but mapping the
    # variant key to its QUAL float rather than just membership.
    quals = {n: variant_quals(present_paths[n]) for n in caller_sel.value}
    return (quals,)


@app.cell
def _(mo, quals):
    # Slider to cap the QUAL axis range on the hexbin plots, so the
    # dense low/mid-QUAL bulk is readable instead of getting squashed
    # by a handful of very-high-QUAL outliers. Default value is the
    # 95th percentile of pooled QUALs — focuses on the bulk while
    # still showing some tail. Re-drawn live as the slider moves.
    all_quals = sorted(q for caller in quals.values() for q in caller.values())
    if all_quals:
        data_max = float(all_quals[-1])
        p95 = float(all_quals[int(len(all_quals) * 0.95)])
    else:
        # Inert defaults so the widget still renders when nothing is selected.
        data_max = 100.0
        p95 = 100.0
    # ~500 slider positions across the data range — gives smooth
    # sliding for both small (data_max ~ 100) and large
    # (data_max ~ 10 000) value ranges.
    step = max(1.0, data_max / 500)
    qual_cap = mo.ui.slider(
        start=1.0,
        stop=data_max,
        step=step,
        value=p95,
        label="Max QUAL on hexbin axes",
        show_value=True,
        full_width=True,
    )
    qual_cap
    return all_quals, data_max, p95, qual_cap, step


@app.cell
def _(mo, plt, qual_cap, quals):
    # Pairwise hexbins of QUAL on the variants both callers called, one panel
    # per pair. Hexbin (vs scatter) handles the 10⁵-record cohort output
    # cleanly; log colour scale so dense ridges and sparse tails are both
    # readable. A red y=x reference line shows where the two callers would
    # agree on confidence — each caller defines QUAL as −10·log10 P over a
    # slightly different event, so the bulk rarely sits on the diagonal and it
    # is the *slope* and *spread* that carry the information.
    #
    # `plt` is dependency-injected from the Venn cell, which already owns the
    # matplotlib import — marimo enforces single-definition for every name.
    _names = list(quals)
    _pairs = [
        (a, b) for i, a in enumerate(_names) for b in _names[i + 1:]
    ]
    if not _pairs:
        hex_view = mo.md("_Select at least two callers._")
    else:
        cap = float(qual_cap.value)
        _ncol = min(3, len(_pairs))
        _nrow = (len(_pairs) + _ncol - 1) // _ncol
        hex_fig, hex_axes = plt.subplots(
            _nrow, _ncol, figsize=(6 * _ncol, 5.5 * _nrow), squeeze=False
        )
        _flat = [hex_axes[r][c] for r in range(_nrow) for c in range(_ncol)]
        for hex_ax, (a, b) in zip(_flat, _pairs):
            qa = quals[a]
            qb = quals[b]
            all_shared = set(qa) & set(qb)
            # Apply the cap symmetrically: both axes need to fall in
            # range, otherwise we'd be looking at variants where one
            # caller is extreme and the other isn't — interesting but
            # for a different view.
            in_cap = [(qa[k], qb[k]) for k in all_shared if qa[k] <= cap and qb[k] <= cap]
            if not in_cap:
                hex_ax.set_title(f"{a} vs {b} — no shared variants ≤ cap")
                hex_ax.set_axis_off()
                continue
            xs = [p[0] for p in in_cap]
            ys = [p[1] for p in in_cap]
            hb = hex_ax.hexbin(
                xs, ys,
                gridsize=40,
                mincnt=1,
                cmap="viridis",
                bins="log",
                extent=(0, cap, 0, cap),
            )
            hex_ax.plot([0, cap], [0, cap], "r--", alpha=0.5, linewidth=1, label="y = x")
            hex_ax.set_xlim(0, cap)
            hex_ax.set_ylim(0, cap)
            hex_ax.set_xlabel(f"{a} QUAL")
            hex_ax.set_ylabel(f"{b} QUAL")
            hex_ax.set_title(
                f"{a} vs {b}  ({len(in_cap)} of {len(all_shared)} ≤ cap)"
            )
            hex_fig.colorbar(hb, ax=hex_ax, label="count (log)")
            hex_ax.legend(loc="upper left", fontsize=8)
        for _spare in _flat[len(_pairs):]:
            _spare.set_axis_off()
        hex_fig.suptitle(
            f"2. Pairwise QUAL agreement on shared variants  (cap = {cap:.0f})",
            y=1.02,
            fontsize=13,
        )
        hex_fig.tight_layout()
        hex_view = hex_fig
    hex_view
    return (hex_view,)


@app.cell
def _(mo):
    # Y-axis scale for the per-caller QUAL distributions below.
    # Default linear; flip to log when a heavy low-QUAL tail dominates
    # the bulk and squashes shape detail.
    qual_yscale = mo.ui.radio(
        options=["linear", "log"],
        value="linear",
        label="QUAL distribution y-scale",
    )
    qual_yscale
    return (qual_yscale,)


@app.cell
def _(mo, plt, qual_cap, qual_yscale, quals):
    # Per-caller QUAL distributions. One panel per caller stacked vertically,
    # all sharing the same x-axis (0..cap) and bin edges so bar positions line
    # up across rows; y-axes are independent because callers can differ by an
    # order of magnitude in record count (freebayes typically emits many more
    # low-QUAL sites than the others) and a shared y would squash the smaller
    # panels.
    _callers = list(quals)
    if not _callers:
        dist_view = mo.md("_Select at least one caller._")
    else:
        dist_cap = float(qual_cap.value)
        n_bins = 50
        bin_edges = [dist_cap * i / n_bins for i in range(n_bins + 1)]
        dist_fig, dist_axes = plt.subplots(
            len(_callers), 1,
            figsize=(10, 2.7 * len(_callers)),
            sharex=True,
            sharey=False,
            squeeze=False,
        )
        for _row, _name in zip(dist_axes, _callers):
            dist_ax = _row[0]
            _values = [q for q in quals[_name].values() if q <= dist_cap]
            _total = len(quals[_name])
            dist_ax.hist(_values, bins=bin_edges, color="steelblue", edgecolor="white")
            dist_ax.set_yscale(qual_yscale.value)
            dist_ax.set_ylabel(f"{_name}\ncount")
            dist_ax.annotate(
                f"n = {len(_values):,} of {_total:,} ≤ cap",
                xy=(0.99, 0.92),
                xycoords="axes fraction",
                ha="right",
                va="top",
                fontsize=9,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7),
            )
        dist_axes[-1][0].set_xlabel("QUAL")
        dist_axes[-1][0].set_xlim(0.0, dist_cap)
        dist_fig.suptitle(
            f"3. Per-caller QUAL distributions  "
            f"(cap = {dist_cap:.0f}, {n_bins} bins, y={qual_yscale.value})",
            y=0.995,
            fontsize=13,
        )
        dist_fig.tight_layout()
        dist_view = dist_fig
    dist_view
    return (dist_view,)


if __name__ == "__main__":
    app.run()
