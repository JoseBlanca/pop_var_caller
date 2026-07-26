# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "numpy", "matplotlib"]
# ///
"""ng STR delimiter bake-off — the "shape of the data" across period × tract-length.

The empirical input to ng_proposal.md §step-3's open question: *how much stuttering
exists for different repeat sizes and lengths, and is there a difference between the
aligners?* Three tract delimiters are compared on the same loci:

  * flat_gap    = algorithm 3, SsrFlatGapAligner   (one flat per-base gap inside the
                  tract; the production-parity port)
  * unit_slip   = algorithm 4, SsrUnitSlipAligner  (whole-unit slips priced apart,
                  direction-asymmetric; the former default)
  * unit_robust = algorithm 4u, SsrUnitRobustAligner (algorithm 4 plus a narrow
                  junction guard and an evidence-based capped anchor test; the
                  current default — the better ruler)

For every microsatellite tract the ng STR locus generator emits, each read is
delimited by ALL THREE aligners and its observation tagged complete / partial /
no-border. This notebook bins those observations by motif period (1-6) and reference
tract length (bp) and shows, per cell: the delimiter-outcome mix, the stutter
distribution, and where two chosen aligners disagree — by default unit_slip vs
unit_robust, the old ruler against the new one.

Data comes from the paired dump tool (one region-typing walk, all three aligners):

    ./scripts/dev.sh cargo run --release --example ng_ssr_aligner_bakeoff -- \\
        benchmarks/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \\
        benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam \\
        chr20 chr21 chr22 > tmp/bakeoff/chr20_22_30x.tsv

Point the notebook at a different dump with PVC_BAKEOFF_TSV=<path>.

Run:  uv run marimo run  benchmarks/ssr_hg002/scripts/ng_ssr_aligner_bakeoff_dashboard.py
Edit: uv run marimo edit benchmarks/ssr_hg002/scripts/ng_ssr_aligner_bakeoff_dashboard.py
"""
import marimo

app = marimo.App(width="medium")


@app.cell
def _():
    import os
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm

    # dataviz palette — the first three categorical slots (validated all-pairs, both modes),
    # assigned so the headline comparison (old default vs new) gets the strongest separation.
    ALIGNER_C = {
        "unit_robust": "#2a78d6",  # slot 1 (blue)   — the current default, algorithm 4u
        "unit_slip": "#eb6834",  # slot 2 (orange) — the former default, algorithm 4
        "flat_gap": "#1baf7a",  # slot 3 (aqua)   — the production-parity port, algorithm 3
    }
    ALIGNERS = list(ALIGNER_C)
    ALIGNER_LABEL = {
        "unit_robust": "unit_robust (algo 4u, default)",
        "unit_slip": "unit_slip (algo 4, former default)",
        "flat_gap": "flat_gap (algo 3, port)",
    }
    # Outcome categories (complete / partial / no-border): categorical slots 1,2 + recessive grey.
    OUTCOME_C = {"complete": "#2a78d6", "partial": "#eb6834", "no_border": "#b8b6ad"}
    # Cells below this many complete reads (grid 2) / this many loci (grid 3) are too sparse to
    # characterise are left blank rather than shown as noise. The chr20–22 slice needs this most;
    # the whole-genome dump fills the penta/hexa long-tract cells that the slice leaves thin.
    MIN_COMPLETE_READS = 8
    MIN_LOCI = 3
    return (
        ALIGNER_C,
        ALIGNER_LABEL,
        ALIGNERS,
        MIN_COMPLETE_READS,
        MIN_LOCI,
        Normalize,
        OUTCOME_C,
        Path,
        TwoSlopeNorm,
        mo,
        np,
        os,
        pd,
        plt,
    )


@app.cell
def _(mo):
    mo.md(
        r"""
        # ng STR delimiter bake-off — the shape of the data (period × tract length)

        Three tract delimiters, the **same loci**, one region-typing walk. Each read is
        delimited by all three and tagged **complete** (both borders seen — an exact
        length), **partial** (one border, ran off its own end — a *lower bound*), or
        **no-border** (reached the aligner, anchored nothing).

        | aligner | algorithm | tract-gap model |
        |---|---|---|
        | **unit_robust** | 4u (`SsrUnitRobustAligner`) | algorithm 4 plus a narrow junction guard (a sequencing error can no longer slide the tract boundary) and an evidence-based capped anchor test (a *complete* is reported only when the flank was actually matched) — **the current default** |
        | **unit_slip** | 4 (`SsrUnitSlipAligner`) | whole-unit slips priced apart, direction-asymmetric — the former default |
        | **flat_gap** | 3 (`SsrFlatGapAligner`) | one flat per-base gap inside the tract — the production-parity port |

        **Stutter** here = `len(observed) − len(reference tract)` in bp: negative is a
        contraction, positive an expansion. The routing-frontier question is *where in
        (period, length) does the stutter/partial/divergence load make the STR path pay
        for itself* — this is its empirical picture.

        The old ruler **fabricates completes**: on a tract long enough that the read runs
        off inside it, `unit_slip` still reports a short exact length. Those enter the
        stutter distribution as spurious short alleles. `unit_robust` reports them as the
        lower bounds they are, so its complete-read stutter shape should be cleaner where
        tracts are long — and its complete→partial frontier should sit earlier.
        """
    )
    return


@app.cell
def _(ALIGNERS, Path, mo, os, pd):
    _here = Path(__file__).resolve()
    _durable = _here.parent.parent / "results" / "ng_aligner_bakeoff" / "chr20_22_30x.tsv"
    _scratch = _here.parents[3] / "tmp" / "bakeoff" / "chr20_22_30x.tsv"
    _default = _durable if _durable.exists() else _scratch
    tsv_path = Path(os.environ.get("PVC_BAKEOFF_TSV") or _default)
    mo.stop(
        not tsv_path.is_file(),
        mo.md(
            f"**Missing** `{tsv_path}`. Generate it (≈110 s for chr20–22, ≈28 min genome-wide):\n\n"
            "```\n./scripts/dev.sh cargo run --release --example ng_ssr_aligner_bakeoff -- \\\n"
            "  benchmarks/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \\\n"
            "  benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam \\\n"
            "  chr20 chr21 chr22 > tmp/bakeoff/chr20_22_30x.tsv\n```\n\n"
            "Drop the trailing contig names for a whole-genome dump. Point the notebook at any "
            "dump with `PVC_BAKEOFF_TSV`."
        ),
    )

    # The `#`-prefixed header carries the authoritative run-level counts, one line per aligner.
    header = {}
    with open(tsv_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            parts = line[1:].split()
            kv = dict(p.split("=", 1) for p in parts if "=" in p)
            header[kv["aligner"]] = {k: int(v) for k, v in kv.items() if k != "aligner"}
    run_counts = pd.DataFrame(header).T.rename_axis("aligner").reset_index()

    # A dump predating unit_robust would otherwise fail deep inside a comparison with a bare
    # KeyError; say what is missing and how to refresh it instead.
    _absent = [a for a in ALIGNERS if a not in header]
    mo.stop(
        bool(_absent),
        mo.md(
            f"`{tsv_path.name}` carries only **{', '.join(header) or 'no'}** — this notebook needs "
            f"**{', '.join(_absent)}** as well. Regenerate the dump with the current "
            "`examples/ng_ssr_aligner_bakeoff.rs`, which runs every delimiter in one walk."
        ),
    )

    df = pd.read_csv(tsv_path, sep="\t", comment="#", dtype={"observed": str})
    df["observed"] = df["observed"].fillna("")
    df["period"] = df["motif"].str.len()
    df["ref_len"] = df["ref_tract"].str.len()
    df["obs_len"] = df["observed"].str.len()
    df["stutter"] = df["obs_len"] - df["ref_len"]
    return df, run_counts, tsv_path


@app.cell
def _(np, pd):
    # Shared axes for every cell grid: motif period (rows) and reference tract-length band (cols).
    PERIOD_NAME = {1: "mono", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
    PERIODS = [1, 2, 3, 4, 5, 6]
    # bp bands chosen against 148 bp reads: spanning degrades as the tract approaches read length.
    LEN_EDGES = [0, 10, 15, 20, 30, 40, 60, 90, 10**9]
    LEN_LABELS = ["<10", "10-14", "15-19", "20-29", "30-39", "40-59", "60-89", "90+"]

    def len_band(series):
        return pd.cut(series, bins=LEN_EDGES, labels=LEN_LABELS, right=False, ordered=True)

    def period_label(p):
        return f"{PERIOD_NAME.get(p, f'{p}-mer')} ({p})"

    def empty_grid():
        return np.full((len(PERIODS), len(LEN_LABELS)), np.nan)

    def annotate_grid(ax, grid, cmap, norm, fmt, ncnt=None):
        """Write each finite cell's value, picking white/dark text by the cell fill's luminance so
        the number stays legible on both pale and saturated cells (the dataviz contrast rule)."""
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                v = grid[i, j]
                if not np.isfinite(v):
                    continue
                r, g, b, _ = cmap(norm(v))
                lum = 0.299 * r + 0.587 * g + 0.114 * b
                txt = fmt.format(v)
                if ncnt is not None and np.isfinite(ncnt[i, j]):
                    txt = f"{txt}\nn={int(ncnt[i, j])}"
                ax.text(
                    j, i, txt, ha="center", va="center", fontsize=6.2,
                    color="white" if lum < 0.5 else "#222",
                )

    def grid_axes(ax):
        ax.set_xticks(range(len(LEN_LABELS)))
        ax.set_xticklabels(LEN_LABELS, rotation=45, ha="right", fontsize=7.5)
        ax.set_yticks(range(len(PERIODS)))
        ax.set_yticklabels([period_label(p) for p in PERIODS], fontsize=8)
        ax.set_xlabel("reference tract length (bp)")

    return (
        LEN_LABELS,
        PERIODS,
        PERIOD_NAME,
        annotate_grid,
        empty_grid,
        grid_axes,
        len_band,
        period_label,
    )


@app.cell
def _(df, len_band):
    # One tidy frame the plot cells share: every row tagged with its period/length cell.
    cells = df.copy()
    cells["len_band"] = len_band(cells["ref_len"])
    cells["cls"] = cells["coverage"].map(
        {
            "complete": "complete",
            "partial_left": "partial",
            "partial_right": "partial",
            "no_border": "no_border",
            "capped": "capped",
        }
    )
    return (cells,)


@app.cell
def _(ALIGNERS, cells, df, mo, run_counts):
    # Headline accounting + the one-number old-ruler-vs-new-ruler result.
    def _fmt_counts(rc):
        r = rc.copy()
        r["on_ref_frac"] = r["aligner"].map(_on_ref_frac(cells))
        cols = [
            "aligner",
            "obs_complete",
            "obs_partial",
            "reads_without_observation",
            "reads_fetched",
            "on_ref_frac",
        ]
        r = r[cols].rename(columns={"reads_without_observation": "no_border"})
        r["on_ref_frac"] = r["on_ref_frac"].map(lambda x: f"{x:.3f}")
        # Table order = the palette order: current default first.
        return r.set_index("aligner").reindex(ALIGNERS).reset_index()

    def _on_ref_frac(cells):
        out = {}
        comp = cells[cells["cls"] == "complete"]
        for al, g in comp.groupby("aligner"):
            out[al] = float(g.loc[g["stutter"] == 0, "reads"].sum()) / float(g["reads"].sum())
        return out

    # Per-locus complete-length evidence, per aligner, compared old default vs new — the whole
    # histogram (which any withdrawn read changes) and the mode (the length a caller would read
    # off the locus). The two answer different questions and are reported as such.
    def _len_hist(g):
        return tuple(sorted((int(l), int(r)) for l, r in zip(g["obs_len"], g["reads"])))

    def _mode_len(g):
        return int(sorted(zip(-g["reads"], g["obs_len"]))[0][1])

    _comp = cells[cells["cls"] == "complete"]

    def _per_locus(fn):
        h = (
            _comp.groupby(["contig", "start", "end", "aligner"])
            .apply(fn, include_groups=False)
            .unstack("aligner")[["unit_slip", "unit_robust"]]
            .dropna()
        )
        return int((h["unit_slip"] != h["unit_robust"]).sum()), len(h)

    _diverge, _n = _per_locus(_len_hist)
    _mode_diverge, _ = _per_locus(_mode_len)

    _rc = run_counts.set_index("aligner")
    _d_complete = int(_rc.loc["unit_robust", "obs_complete"] - _rc.loc["unit_slip", "obs_complete"])
    _d_partial = int(_rc.loc["unit_robust", "obs_partial"] - _rc.loc["unit_slip", "obs_partial"])

    n_loci = df[["contig", "start", "end"]].drop_duplicates().shape[0]
    # Name the slice from the data, not from a hardcoded label — the same notebook serves the
    # chr20–22 dump and the whole-genome one.
    _contigs = sorted(df["contig"].unique())
    scope = ", ".join(_contigs) if len(_contigs) <= 4 else f"{len(_contigs)} contigs"
    mo.vstack(
        [
            mo.md(
                f"**{n_loci:,} covered loci** ({scope}, HG002 30×). Swapping the old ruler for "
                f"the new one moves **{_d_complete:+,} complete** reads and **{_d_partial:+,} "
                f"partial** — reads `unit_slip` gave an exact length that `unit_robust` will only "
                f"bound. Of the **{_n:,}** loci with a complete observation under both, "
                f"**{100 * _diverge / max(_n, 1):.0f}%** get a *different* complete-length "
                f"histogram, but only **{_mode_diverge:,} "
                f"({100 * _mode_diverge / max(_n, 1):.1f}%)** change their **modal** length. "
                f"The ruler swap re-weighs the evidence at nearly every locus; it re-calls the "
                f"length at very few. Run-level accounting "
                f"(`reads_fetched = complete + partial + no_border`):"
            ),
            mo.ui.table(_fmt_counts(run_counts), selection=None, pagination=False),
        ]
    )
    return n_loci, scope


@app.cell
def _(PERIODS, PERIOD_NAME, cells, mo, plt, run_counts):
    # SECTION 1 — the catalog itself: how many STRs of each period the walk actually covered. Every
    # later grid is a slice of this, so it belongs first; it is aligner-independent (all three
    # delimiters see the same locus set), hence one bar per period rather than a series per ruler.
    def catalog_figure():
        loci = cells[["contig", "start", "end", "period"]].drop_duplicates()
        counts = [int((loci["period"] == p).sum()) for p in PERIODS]
        total = sum(counts)

        fig, ax = plt.subplots(figsize=(9, 3.4))
        x = range(len(PERIODS))
        # A single measure with no aligner identity — recessive ink, not a categorical hue.
        ax.bar(x, counts, color="#8a8880", width=0.7)
        for xi, c in zip(x, counts):
            ax.text(
                xi,
                c + total * 0.012,
                f"{c:,}\n{c / total:.0%}",
                ha="center",
                va="bottom",
                fontsize=8,
                color="#52514e",
            )
        ax.set_xticks(list(x))
        ax.set_xticklabels([f"{PERIOD_NAME[p]}\n({p} bp)" for p in PERIODS], fontsize=8.5)
        ax.set_xlabel("motif period")
        ax.set_ylabel("covered STR loci")
        ax.set_ylim(0, max(counts) * 1.18)
        ax.grid(True, axis="y", alpha=0.25)
        ax.set_axisbelow(True)
        fig.suptitle(
            f"The covered catalog — {total:,} STR loci by motif period",
            fontweight="bold",
            fontsize=11,
        )
        fig.tight_layout()
        return fig, counts, total

    catalog_fig, _counts, _total = catalog_figure()
    _mono, _di = _counts[0], _counts[1]
    # The walk types far more tracts than the reads ever reach; the header carries both totals.
    _walked = int(run_counts["ssr_loci"].iloc[0])
    _zero = int(run_counts["zero_coverage"].iloc[0])
    mo.vstack(
        [
            mo.md(
                "### 1 · How many STRs of each period\n"
                "The locus set every later grid slices. **Homopolymers dominate it** — mono and "
                f"di alone are {(_mono + _di) / _total:.0%} of the {_total:,} covered loci — so "
                "any figure pooled over all periods is mostly a statement about them. That is "
                "why the rest of this notebook is binned by period rather than pooled, and why "
                "the penta/hexa cells stay thin. The tail is also not monotonic: **tetra "
                "outnumbers tri**."
            ),
            catalog_fig,
            mo.md(
                f"*Covered loci only — a locus no read reaches emits no rows and is absent here. "
                f"That is a severe filter, and it is what shapes the mix above: the walk typed "
                f"{_walked:,} STR tracts and {_zero:,} ({_zero / _walked:.1%}) of them carry no "
                f"reads at all, because this HG002 BAM is sliced to the GIAB tandem-repeat "
                f"benchmark's 50,000-region Tier subset. Read the counts as the composition of "
                f"the benchmark's covered regions, not of the genome. They are "
                f"aligner-independent — one region-typing walk feeds all three delimiters, so "
                f"they see the identical locus set.*"
            ),
        ]
    )
    return


@app.cell
def _(ALIGNER_C, ALIGNER_LABEL, ALIGNERS, mo):
    aligner_sel = mo.ui.radio(
        {ALIGNER_LABEL[a]: a for a in ALIGNERS},
        value=ALIGNER_LABEL["unit_robust"],
        label="Aligner (grids 2–3 below)",
    )
    mo.md(
        f"The two grids below are shown for **one** aligner at a time — switch here. "
        f"<span style='color:{ALIGNER_C['unit_robust']}'>■</span> unit_robust is the default; "
        f"flip to <span style='color:{ALIGNER_C['unit_slip']}'>■</span> unit_slip to see the "
        f"same cell under the old ruler."
    )
    return (aligner_sel,)


@app.cell
def _(aligner_sel):
    aligner_sel
    return


@app.cell
def _(
    LEN_LABELS,
    OUTCOME_C,
    PERIODS,
    aligner_sel,
    cells,
    mo,
    period_label,
    plt,
):
    # GRID 2 — delimiter-outcome mix per cell: what fraction of a cell's fetched reads land as
    # complete vs partial vs no-border. This is where the tract outgrows the read.
    def outcome_figure(aligner):
        sub = cells[(cells["aligner"] == aligner) & (cells["cls"] != "capped")]
        fig, axes = plt.subplots(2, 3, figsize=(12, 6.6), sharey=True)
        order = ["complete", "partial", "no_border"]
        for ax, period in zip(axes.flat, PERIODS):
            g = sub[sub["period"] == period]
            reads = (
                g.groupby(["len_band", "cls"], observed=False)["reads"]
                .sum()
                .unstack("cls")
                .reindex(columns=order, fill_value=0)
                .reindex(LEN_LABELS)
            )
            totals = reads.sum(axis=1)
            frac = reads.div(totals.replace(0, 1), axis=0)
            bottom = 0
            x = range(len(LEN_LABELS))
            for cls in order:
                ax.bar(
                    x,
                    frac[cls].values,
                    bottom=bottom,
                    color=OUTCOME_C[cls],
                    width=0.82,
                    label=cls,
                    edgecolor="white",
                    linewidth=0.6,
                )
                bottom = bottom + frac[cls].values
            for xi, t in zip(x, totals.values):
                if t > 0:
                    ax.text(xi, 1.02, f"{int(t)}", ha="center", va="bottom", fontsize=6.5, color="#52514e")
            ax.set_title(period_label(period), fontsize=10)
            ax.set_xticks(list(x))
            ax.set_xticklabels(LEN_LABELS, rotation=45, ha="right", fontsize=7.5)
            ax.set_ylim(0, 1.12)
            ax.margins(x=0.02)
        for ax in axes[:, 0]:
            ax.set_ylabel("fraction of fetched reads")
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.005))
        fig.suptitle(
            f"Delimiter outcomes by period × tract length — {aligner}  (n above each bar = fetched reads)",
            fontweight="bold",
            y=1.05,
        )
        fig.supxlabel("reference tract length (bp)", y=-0.02)
        fig.tight_layout()
        return fig

    mo.vstack(
        [
            mo.md(
                "### 2 · Delimiter outcomes — complete / partial / no-border\n"
                "As a tract lengthens toward the 148 bp read, fewer reads span it: **complete** "
                "gives way to **partial** (a censored lower bound) and eventually **no-border**. "
                "The frontier where complete collapses is period-dependent."
            ),
            outcome_figure(aligner_sel.value),
        ]
    )
    return


@app.cell
def _(
    LEN_LABELS,
    MIN_COMPLETE_READS,
    Normalize,
    PERIODS,
    TwoSlopeNorm,
    aligner_sel,
    annotate_grid,
    cells,
    empty_grid,
    grid_axes,
    mo,
    np,
    plt,
):
    # GRID 3 — stutter shape over complete reads: how much (off-reference fraction) and which way
    # (mean signed bp). Two heatmaps sharing the period × length grid. Cells with too few complete
    # reads to characterise (MIN_COMPLETE_READS) are left blank rather than shown as noise.
    def stutter_figure(aligner):
        comp = cells[(cells["aligner"] == aligner) & (cells["cls"] == "complete")]
        off = empty_grid()
        signed = empty_grid()
        n = empty_grid()
        for i, period in enumerate(PERIODS):
            for j, band in enumerate(LEN_LABELS):
                g = comp[(comp["period"] == period) & (comp["len_band"] == band)]
                w = g["reads"].to_numpy()
                if w.sum() < MIN_COMPLETE_READS:
                    continue
                s = g["stutter"].to_numpy()
                n[i, j] = w.sum()
                off[i, j] = w[s != 0].sum() / w.sum()
                signed[i, j] = float((s * w).sum() / w.sum())

        fig, axes = plt.subplots(1, 2, figsize=(13, 4.3))

        cmap0 = plt.get_cmap("Blues")
        norm0 = Normalize(vmin=0, vmax=max(np.nanmax(off) if np.isfinite(off).any() else 0.05, 0.05))
        im0 = axes[0].imshow(off, cmap=cmap0, norm=norm0, aspect="auto")
        axes[0].set_title("off-reference read fraction (stutter intensity)", fontsize=10)
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.03)
        # n = complete reads behind the cell, so a one-locus cell can be discounted on sight.
        annotate_grid(axes[0], off, cmap0, norm0, "{:.2f}", ncnt=n)

        # Robust cap so a couple of tiny-n outliers don't wash out the bulk contraction-bias signal.
        _absv = np.abs(signed[np.isfinite(signed)])
        _cap = max(float(np.nanpercentile(_absv, 90)) if _absv.size else 2.0, 2.0)
        # RdBu_r, not RdBu: the reversed ramp is the one that puts BLUE at the negative
        # (contraction) end, which is what the panel title and the section text both claim.
        cmap1 = plt.get_cmap("RdBu_r")
        norm1 = TwoSlopeNorm(vcenter=0, vmin=-_cap, vmax=_cap)
        im1 = axes[1].imshow(np.clip(signed, -_cap, _cap), cmap=cmap1, norm=norm1, aspect="auto")
        axes[1].set_title(
            f"mean signed stutter (bp), clipped ±{_cap:.0f} — blue = contraction bias", fontsize=9.5
        )
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.03)
        # Annotate with the TRUE value (not the clipped one), coloured by the clipped fill.
        annotate_grid(axes[1], signed, cmap1, norm1, "{:+.1f}")

        for ax in axes:
            grid_axes(ax)
        fig.suptitle(f"Stutter shape over complete reads — {aligner}", fontweight="bold")
        fig.tight_layout()
        return fig

    mo.vstack(
        [
            mo.md(
                "### 3 · Stutter shape — how much, and which way\n"
                "Left: fraction of complete reads whose length ≠ reference (stutter load). Right: "
                "mean signed stutter in bp — **blue = net contraction** (losing units, the known "
                "PCR-stutter asymmetry), red = net expansion. Reads-weighted; the diverging scale is "
                f"clipped at a robust cap so a few tiny-n outliers don't dominate. Cells with fewer "
                f"than **{MIN_COMPLETE_READS}** complete reads are left blank."
            ),
            stutter_figure(aligner_sel.value),
        ]
    )
    return


@app.cell
def _(ALIGNER_LABEL, ALIGNERS, mo):
    baseline_sel = mo.ui.dropdown(
        {ALIGNER_LABEL[a]: a for a in ALIGNERS},
        value=ALIGNER_LABEL["unit_slip"],
        label="baseline (A)",
    )
    challenger_sel = mo.ui.dropdown(
        {ALIGNER_LABEL[a]: a for a in ALIGNERS},
        value=ALIGNER_LABEL["unit_robust"],
        label="challenger (B)",
    )
    mo.md(
        "### 4 · Where two rulers diverge — old default vs new\n"
        "The pair defaults to **unit_slip → unit_robust**: the ruler ng used against the ruler "
        "ng uses now. Every panel below is **B − A**, so red means *the challenger does more of "
        "it*. Pick another pair to recover the original algorithm-3-vs-4 comparison."
    )
    return baseline_sel, challenger_sel


@app.cell
def _(baseline_sel, challenger_sel, mo):
    mo.hstack([baseline_sel, challenger_sel], justify="start", gap=2)
    return


@app.cell
def _(
    LEN_LABELS,
    MIN_LOCI,
    Normalize,
    PERIODS,
    TwoSlopeNorm,
    annotate_grid,
    baseline_sel,
    cells,
    challenger_sel,
    empty_grid,
    grid_axes,
    mo,
    np,
    pd,
    plt,
):
    # GRID 4 — where two chosen aligners diverge, per cell. Three heatmaps: (a) fraction of loci
    # whose MODAL complete length differs, (b) B − A off-reference fraction among complete reads
    # (is the challenger's surviving stutter shape cleaner?), (c) B − A complete-read fraction —
    # the complete→partial frontier shift, i.e. how many exact lengths the challenger declines to
    # claim. Cells with fewer than MIN_LOCI loci are left blank in (a)/(b) so a single-locus cell
    # can't read as "100% divergence".
    #
    # Panel (a) uses the mode, not the whole histogram, precisely because (c) is non-zero: once a
    # ruler withdraws reads to partial, the full histogram differs at almost every locus with
    # enough depth, which measures the frontier move a second time instead of the thing a caller
    # would notice — the length it reads off the locus.
    def _mode_len(g):
        # Heaviest complete length; ties broken by the shorter length so the answer is stable.
        best = sorted(zip(-g["reads"], g["obs_len"]))
        return int(best[0][1])

    def divergence_grids(a, b):
        comp = cells[cells["cls"] == "complete"]
        # per-locus modal complete length per aligner
        h = (
            comp.groupby(["contig", "start", "end", "period", "len_band", "aligner"], observed=True)
            .apply(_mode_len, include_groups=False)
            .unstack("aligner")
        )
        h = h.reset_index()
        for al in (a, b):
            if al not in h.columns:
                h[al] = None
        # Only loci with a complete observation under BOTH rulers have two modes to compare; where
        # one ruler emits nothing complete the difference is a frontier move, which panel (c)
        # measures instead.
        h = h.dropna(subset=[a, b])
        h["diverge"] = h[a] != h[b]

        diverge = empty_grid()
        nloci = empty_grid()
        for i, period in enumerate(PERIODS):
            for j, band in enumerate(LEN_LABELS):
                sub = h[(h["period"] == period) & (h["len_band"] == band)]
                if len(sub) < MIN_LOCI:
                    continue
                nloci[i, j] = len(sub)
                diverge[i, j] = sub["diverge"].mean()

        # B − A off-reference fraction among complete reads, per cell (reads-weighted)
        def off_ref(al):
            g = comp[comp["aligner"] == al]
            return g.groupby(["period", "len_band"], observed=True).apply(
                lambda x: pd.Series(
                    {"off": x.loc[x["stutter"] != 0, "reads"].sum(), "tot": x["reads"].sum()}
                ),
                include_groups=False,
            )

        ga, gb = off_ref(a), off_ref(b)
        delta = empty_grid()
        for i, period in enumerate(PERIODS):
            for j, band in enumerate(LEN_LABELS):
                key = (period, band)
                if key in ga.index and key in gb.index and ga.loc[key, "tot"] and gb.loc[key, "tot"]:
                    delta[i, j] = (
                        gb.loc[key, "off"] / gb.loc[key, "tot"]
                        - ga.loc[key, "off"] / ga.loc[key, "tot"]
                    )
        # blank the delta panel where the divergence panel is too sparse, so the two agree on scope
        delta[~np.isfinite(diverge)] = np.nan

        # B − A complete-read fraction over ALL fetched reads — the frontier shift. Denominator is
        # every non-capped read in the cell, which both rulers share, so the difference is purely
        # reclassification.
        def complete_frac(al):
            g = cells[(cells["aligner"] == al) & (cells["cls"] != "capped")]
            return g.groupby(["period", "len_band"], observed=True).apply(
                lambda x: pd.Series(
                    {
                        "comp": x.loc[x["cls"] == "complete", "reads"].sum(),
                        "tot": x["reads"].sum(),
                    }
                ),
                include_groups=False,
            )

        fa, fb = complete_frac(a), complete_frac(b)
        frontier = empty_grid()
        for i, period in enumerate(PERIODS):
            for j, band in enumerate(LEN_LABELS):
                key = (period, band)
                if key in fa.index and key in fb.index and fa.loc[key, "tot"] and fb.loc[key, "tot"]:
                    frontier[i, j] = (
                        fb.loc[key, "comp"] / fb.loc[key, "tot"]
                        - fa.loc[key, "comp"] / fa.loc[key, "tot"]
                    )
        return diverge, nloci, delta, frontier

    def divergence_figure(a, b):
        diverge, nloci, delta, frontier = divergence_grids(a, b)
        fig, axes = plt.subplots(1, 3, figsize=(17, 4.3))

        cmap0 = plt.get_cmap("Purples")
        norm0 = Normalize(vmin=0, vmax=1.0)
        im0 = axes[0].imshow(diverge, cmap=cmap0, norm=norm0, aspect="auto")
        axes[0].set_title("loci whose modal complete length differs", fontsize=9.5)
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.03)
        annotate_grid(axes[0], diverge, cmap0, norm0, "{:.0%}", ncnt=nloci)

        _lim = max(float(np.nanmax(np.abs(delta))) if np.isfinite(delta).any() else 0.05, 0.05)
        cmap1 = plt.get_cmap("RdBu_r")
        norm1 = TwoSlopeNorm(vcenter=0, vmin=-_lim, vmax=_lim)
        im1 = axes[1].imshow(delta, cmap=cmap1, norm=norm1, aspect="auto")
        axes[1].set_title(
            f"off-ref frac among completes: {b} − {a}\n(red = B keeps more reads off REF)",
            fontsize=9.5,
        )
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.03)
        annotate_grid(axes[1], delta, cmap1, norm1, "{:+.2f}")

        _flim = max(
            float(np.nanmax(np.abs(frontier))) if np.isfinite(frontier).any() else 0.05, 0.05
        )
        norm2 = TwoSlopeNorm(vcenter=0, vmin=-_flim, vmax=_flim)
        im2 = axes[2].imshow(frontier, cmap=cmap1, norm=norm2, aspect="auto")
        axes[2].set_title(
            f"complete-read frac: {b} − {a}\n(blue = B claims fewer exact lengths)", fontsize=9.5
        )
        fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.03)
        annotate_grid(axes[2], frontier, cmap1, norm2, "{:+.2f}")

        for ax in axes:
            grid_axes(ax)
        fig.suptitle(f"{a} → {b}: divergence by period × tract length", fontweight="bold")
        fig.tight_layout()
        return fig

    mo.vstack(
        [
            mo.md(
                "**Left**: fraction of covered loci (with a complete read under both) whose "
                "**modal complete length differs** — the cells where the ruler changes the length "
                "you would read off the locus. (The mode, not the whole histogram: once a ruler "
                "withdraws reads to partial, the full histogram differs nearly everywhere, which "
                "just re-measures the right-hand panel.) "
                "**Middle**: off-reference fraction *among the reads each ruler still calls "
                "complete*; red = B keeps more of them off REF, blue = B's surviving completes sit "
                "on the reference more often (spurious off-REF alleles withdrawn). **Right**: the "
                "**complete→partial frontier shift** — blue = B declines to claim an exact length "
                "where A did, reporting a lower bound instead. Cells with fewer than "
                f"**{MIN_LOCI}** loci are blank in the first two panels."
            ),
            divergence_figure(baseline_sel.value, challenger_sel.value),
        ]
    )
    return


@app.cell
def _(LEN_LABELS, PERIODS, mo, period_label):
    detail_period = mo.ui.dropdown(
        {period_label(p): p for p in PERIODS}, value=period_label(2), label="period"
    )
    detail_band = mo.ui.dropdown(LEN_LABELS, value="20-29", label="tract length (bp)")
    mo.md("### 5 · Drill into one cell — the allele-length distribution")
    return detail_band, detail_period


@app.cell
def _(detail_band, detail_period, mo):
    mo.hstack([detail_period, detail_band], justify="start", gap=2)
    return


@app.cell
def _(
    ALIGNER_C,
    ALIGNERS,
    cells,
    detail_band,
    detail_period,
    mo,
    np,
    plt,
):
    # The per-cell "shape": reads-weighted complete-stutter distributions, all three rulers
    # overlaid, plus the partial / no-border tallies the frontier grids summarise.
    def detail_figure(period, band):
        sub = cells[(cells["period"] == period) & (cells["len_band"] == band)]
        comp = sub[sub["cls"] == "complete"]
        if comp.empty:
            fig, ax = plt.subplots(figsize=(9, 3))
            ax.text(0.5, 0.5, "no complete reads in this cell", ha="center", va="center")
            ax.axis("off")
            return fig, "—"

        lo = int(comp["stutter"].min())
        hi = int(comp["stutter"].max())
        lo, hi = min(lo, -1), max(hi, 1)
        bins = np.arange(lo - 0.5, hi + 1.5, 1)
        centers = np.arange(lo, hi + 1)

        # Two views of the same curves. Linear carries the headline shape (the reference spike
        # dominates); log carries the tail, which is where the rulers actually differ and which a
        # linear axis flattens to nothing. Only the default is filled — three overlapping fills
        # would be mud.
        fig, axes = plt.subplots(1, 2, figsize=(13, 3.8))
        curves = {}
        for al in ALIGNERS:
            g = comp[comp["aligner"] == al]
            hist, _ = np.histogram(g["stutter"], bins=bins, weights=g["reads"])
            curves[al] = hist / hist.sum() if hist.sum() else hist
        for ax in axes:
            for al in ALIGNERS:
                ax.step(centers, curves[al], where="mid", color=ALIGNER_C[al], linewidth=2, label=al)
            ax.fill_between(
                centers, curves[ALIGNERS[0]], step="mid", color=ALIGNER_C[ALIGNERS[0]], alpha=0.12
            )
            ax.axvline(0, color="#b8b6ad", linewidth=1, zorder=0)
            ax.set_xlabel("stutter Δ = observed − reference (bp)")
            ax.grid(True, axis="y", alpha=0.25)
        axes[0].set_ylabel("fraction of complete reads")
        axes[0].set_title("linear — the bulk", fontsize=10)
        axes[1].set_yscale("log")
        axes[1].set_title("log — the off-reference tail", fontsize=10)
        axes[0].legend(frameon=False)
        fig.suptitle(
            f"Complete-read stutter — period {period}, tract {band} bp", fontweight="bold"
        )
        fig.tight_layout()

        # outcome tally caption
        n_loci = sub[["contig", "start", "end"]].drop_duplicates().shape[0]
        lines = [
            f"**{n_loci} loci** in this cell. Reads by outcome, and the on-reference share of "
            "each ruler's surviving completes:",
            "",
            "| ruler | complete | partial | no-border | on-REF |",
            "|---|---:|---:|---:|---:|",
        ]
        for al in ALIGNERS:
            g = sub[sub["aligner"] == al]
            counts = [
                int(g.loc[g["cls"] == c, "reads"].sum())
                for c in ("complete", "partial", "no_border")
            ]
            gc = g[g["cls"] == "complete"]
            on_ref = gc.loc[gc["stutter"] == 0, "reads"].sum() / max(gc["reads"].sum(), 1)
            lines.append(
                f"| `{al}` | " + " | ".join(f"{c:,}" for c in counts) + f" | {on_ref:.1%} |"
            )
        return fig, "\n".join(lines)

    _fig, _cap = detail_figure(detail_period.value, detail_band.value)
    mo.vstack([_fig, mo.md(_cap)])
    return


@app.cell
def _(ALIGNER_C, PERIODS, cells, np, period_label, plt):
    # Which complete reads did the new ruler withdraw? A fabricated complete is a read that ran
    # off INSIDE the tract, so it reads short — if that is what unit_robust takes back, the
    # withdrawn reads must be enriched for contractions relative to the pool they came from.
    # Lollipop: base contraction rate → contraction share of the withdrawn, one row per period.
    def withdrawal_figure():
        comp = cells[cells["cls"] == "complete"]

        def classes(al, period):
            g = comp[(comp["aligner"] == al) & (comp["period"] == period)]
            return np.array(
                [
                    g.loc[g["stutter"] < 0, "reads"].sum(),
                    g.loc[g["stutter"] == 0, "reads"].sum(),
                    g.loc[g["stutter"] > 0, "reads"].sum(),
                ],
                dtype=float,
            )

        rows = []
        for period in PERIODS:
            slip, robust = classes("unit_slip", period), classes("unit_robust", period)
            withdrawn = slip - robust
            if slip.sum() <= 0 or withdrawn.sum() <= 0:
                continue
            rows.append(
                (period, slip[0] / slip.sum(), withdrawn[0] / withdrawn.sum(), withdrawn.sum())
            )

        fig, ax = plt.subplots(figsize=(9, 3.4))
        y = range(len(rows))
        for i, (_, base, got, _n) in zip(y, rows):
            ax.plot([base, got], [i, i], color="#d6d5cd", linewidth=2.5, zorder=1, solid_capstyle="round")
            ax.scatter([base], [i], s=70, color="#8a8880", zorder=2)
            ax.scatter([got], [i], s=70, color=ALIGNER_C["unit_robust"], zorder=3)
            ax.text(
                got + 0.012, i, f"{got / base:.1f}×", va="center", fontsize=8, color="#52514e"
            )
        ax.set_yticks(list(y))
        ax.set_yticklabels([period_label(p) for p, *_ in rows], fontsize=8.5)
        ax.invert_yaxis()
        ax.set_xlabel("contraction share of complete reads")
        ax.set_xlim(0, max(g for _, _, g, _ in rows) * 1.22)
        ax.xaxis.set_major_formatter(lambda v, _pos: f"{v:.0%}")
        ax.grid(True, axis="x", alpha=0.25)
        ax.scatter([], [], s=70, color="#8a8880", label="all unit_slip completes (base rate)")
        ax.scatter(
            [], [], s=70, color=ALIGNER_C["unit_robust"], label="the completes unit_robust withdrew"
        )
        # Legend above the plot: inside the axes it lands on the bottom rows' marks and labels.
        ax.legend(
            frameon=False,
            fontsize=8.5,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.0),
            ncol=2,
        )
        fig.suptitle(
            "The withdrawn reads are the short ones — and increasingly so with period",
            fontsize=10.5,
            fontweight="bold",
        )
        fig.tight_layout()
        return fig

    withdrawal_fig = withdrawal_figure()
    return (withdrawal_fig,)


@app.cell
def _(cells, mo, pd, withdrawal_fig):
    # SECTION 6 — the read of grids 3 and 4, with its numbers computed from the loaded dump
    # rather than written down, so a re-run (or a whole-genome dump) restates itself.
    def _shape(al):
        comp = cells[(cells["aligner"] == al) & (cells["cls"] == "complete")]
        allr = cells[(cells["aligner"] == al) & (cells["cls"] != "capped")]
        tot = comp["reads"].sum()
        return {
            "complete reads": int(tot),
            "of fetched reads": comp["reads"].sum() / allr["reads"].sum(),
            "on reference": comp.loc[comp["stutter"] == 0, "reads"].sum() / tot,
            "contraction": comp.loc[comp["stutter"] < 0, "reads"].sum() / tot,
            "expansion": comp.loc[comp["stutter"] > 0, "reads"].sum() / tot,
        }

    _t = pd.DataFrame({al: _shape(al) for al in ("unit_slip", "unit_robust")}).T
    _fmt = _t.copy()
    for _c in ("of fetched reads", "on reference", "contraction", "expansion"):
        _fmt[_c] = _t[_c].map(lambda x: f"{x:.2%}")
    _fmt["complete reads"] = _t["complete reads"].map(lambda x: f"{int(x):,}")

    # How lopsided is the withdrawal? Fabricated completes are reads that ran off INSIDE the
    # tract, so they read short — the contraction column should shrink harder than the expansion.
    _c_drop = 1 - _t.loc["unit_robust", "contraction"] / _t.loc["unit_slip", "contraction"]
    _e_drop = 1 - _t.loc["unit_robust", "expansion"] / _t.loc["unit_slip", "expansion"]

    mo.vstack(
        [
            mo.md(
                "### 6 · What the better ruler changed\n"
                "Pooled over every complete read in the dump, old ruler against new:"
            ),
            mo.ui.table(_fmt.rename_axis("ruler").reset_index(), selection=None, pagination=False),
            mo.md(
                f"The withdrawal is **asymmetric**, which is the signature the fix predicts: the "
                f"contraction share of complete reads falls **{_c_drop:.0%}** while the expansion "
                f"share falls only **{_e_drop:.0%}**. A read that runs off *inside* a tract reads "
                f"**short**, so a ruler that fabricates completes manufactures contractions "
                f"specifically — and that is the part `unit_robust` takes back. The on-reference "
                f"share of what survives rises accordingly: the stutter shape left over is "
                f"cleaner, not merely smaller."
            ),
            withdrawal_fig,
            mo.md(
                "Read the lollipop as a test of that claim. If the withdrawn reads were a random "
                "sample of the completes, both dots would sit on top of each other. They do not: "
                "the withdrawn set is several times richer in contractions than the pool it came "
                "from, and the multiple **grows with period** — at mono/di the new ruler is mostly "
                "withdrawing on-reference reads (the anchor-evidence test firing everywhere), "
                "while from tri upward the reads it declines to call are predominantly the short "
                "ones. Homopolymer and dinucleotide contractions largely *survive*, which is the "
                "right outcome: that stutter is real."
            ),
        ]
    )
    return


@app.cell
def _(mo, n_loci, scope, tsv_path):
    mo.md(
        f"""
        ---
        *Source: `{tsv_path.name}` ({scope}) — {n_loci:,} covered loci, all three delimiters, one region-typing
        walk (`examples/ng_ssr_aligner_bakeoff.rs`). Periods 1–6 (region typing emits homopolymers
        too, which §4.2 of the alignment spec deliberately keeps in the comparison). Partials are
        recorded but never fed to a length model — they are censored lower bounds.*
        """
    )
    return


if __name__ == "__main__":
    app.run()
