# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "numpy", "matplotlib"]
# ///
"""ng STR delimiter bake-off — the "shape of the data" across period × tract-length.

The empirical input to ng_proposal.md §step-3's open question: *how much stuttering
exists for different repeat sizes and lengths, and is there a difference between the
aligners?* Two tract delimiters are compared on the same loci:

  * flat_gap  = algorithm 3, SsrFlatGapAligner  (one flat per-base gap inside the
                tract; the production-parity port)
  * unit_slip = algorithm 4, SsrUnitSlipAligner  (whole-unit slips priced apart,
                direction-asymmetric; the recommended default)

For every microsatellite tract the ng STR locus generator emits, each read is
delimited by BOTH aligners and its observation tagged complete / partial / no-border.
This notebook bins those observations by motif period (1-6) and reference tract length
(bp) and shows, per cell: the delimiter-outcome mix, the stutter distribution, and
where the two aligners disagree.

Data comes from the paired dump tool (one region-typing walk, both aligners):

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

    # dataviz palette — categorical slots (CVD-safe pair), reserved for the two aligners.
    UNIT_SLIP_C = "#2a78d6"  # slot 1 (blue)   — the recommended default
    FLAT_GAP_C = "#eb6834"  # slot 2 (orange) — the production-parity port
    # Outcome categories (complete / partial / no-border): categorical slots 1,2 + recessive grey.
    OUTCOME_C = {"complete": "#2a78d6", "partial": "#eb6834", "no_border": "#b8b6ad"}
    # Cells below this many complete reads (grid 2) / this many loci (grid 3) are too sparse to
    # characterise on a chr20–22 slice, so they are left blank rather than shown as noise.
    MIN_COMPLETE_READS = 8
    MIN_LOCI = 3
    return (
        FLAT_GAP_C,
        MIN_COMPLETE_READS,
        MIN_LOCI,
        Normalize,
        OUTCOME_C,
        Path,
        TwoSlopeNorm,
        UNIT_SLIP_C,
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

        Two tract delimiters, the **same loci**, one region-typing walk. Each read is
        delimited by both and tagged **complete** (both borders seen — an exact length),
        **partial** (one border, ran off its own end — a *lower bound*), or **no-border**
        (reached the aligner, anchored nothing).

        | aligner | algorithm | tract-gap model |
        |---|---|---|
        | **unit_slip** | 4 (`SsrUnitSlipAligner`) | whole-unit slips priced apart, direction-asymmetric — the recommended default |
        | **flat_gap** | 3 (`SsrFlatGapAligner`) | one flat per-base gap inside the tract — the production-parity port |

        **Stutter** here = `len(observed) − len(reference tract)` in bp: negative is a
        contraction, positive an expansion. The routing-frontier question is *where in
        (period, length) does the stutter/partial/divergence load make the STR path pay
        for itself* — this is its empirical picture.
        """
    )
    return


@app.cell
def _(Path, mo, os, pd):
    _here = Path(__file__).resolve()
    _durable = _here.parent.parent / "results" / "ng_aligner_bakeoff" / "chr20_22_30x.tsv"
    _scratch = _here.parents[3] / "tmp" / "bakeoff" / "chr20_22_30x.tsv"
    _default = _durable if _durable.exists() else _scratch
    tsv_path = Path(os.environ.get("PVC_BAKEOFF_TSV") or _default)
    mo.stop(
        not tsv_path.is_file(),
        mo.md(
            f"**Missing** `{tsv_path}`. Generate it (≈70 s for chr20–22):\n\n"
            "```\n./scripts/dev.sh cargo run --release --example ng_ssr_aligner_bakeoff -- \\\n"
            "  benchmarks/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \\\n"
            "  benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam \\\n"
            "  chr20 chr21 chr22 > tmp/bakeoff/chr20_22_30x.tsv\n```\n\n"
            "or set `PVC_BAKEOFF_TSV` to another dump."
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
def _(cells, df, mo, run_counts):
    # Headline accounting + the one-number divergence result.
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
        return r

    def _on_ref_frac(cells):
        out = {}
        comp = cells[cells["cls"] == "complete"]
        for al, g in comp.groupby("aligner"):
            out[al] = float(g.loc[g["stutter"] == 0, "reads"].sum()) / float(g["reads"].sum())
        return out

    # Per-locus complete-length histogram, per aligner, compared flat_gap vs unit_slip.
    def _len_hist(g):
        return tuple(sorted((int(l), int(r)) for l, r in zip(g["obs_len"], g["reads"])))

    _comp = cells[cells["cls"] == "complete"]
    _h = (
        _comp.groupby(["contig", "start", "end", "aligner"])
        .apply(_len_hist, include_groups=False)
        .unstack("aligner")
    )
    _both = _h.dropna()
    _diverge = int((_both["flat_gap"] != _both["unit_slip"]).sum())
    _n = len(_both)

    n_loci = df[["contig", "start", "end"]].drop_duplicates().shape[0]
    mo.vstack(
        [
            mo.md(
                f"**{n_loci:,} covered loci** (chr20–22, HG002 30×). Of the **{_n:,}** with a "
                f"complete observation under both aligners, **{_diverge:,} "
                f"({100 * _diverge / max(_n, 1):.1f}%)** get a *different* complete-length "
                f"histogram — the aligner choice is not cosmetic. Run-level accounting "
                f"(`reads_fetched = complete + partial + no_border`):"
            ),
            mo.ui.table(_fmt_counts(run_counts), selection=None, pagination=False),
        ]
    )
    return (n_loci,)


@app.cell
def _(UNIT_SLIP_C, mo):
    aligner_sel = mo.ui.radio(
        {"unit_slip (algo 4, default)": "unit_slip", "flat_gap (algo 3, port)": "flat_gap"},
        value="unit_slip (algo 4, default)",
        label="Aligner (grids below)",
    )
    mo.md(
        f"The three grids below are shown for **one** aligner at a time — switch here. "
        f"<span style='color:{UNIT_SLIP_C}'>■</span> unit_slip is the default."
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
    # GRID 1 — delimiter-outcome mix per cell: what fraction of a cell's fetched reads land as
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
                "### 1 · Delimiter outcomes — complete / partial / no-border\n"
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
    # GRID 2 — stutter shape over complete reads: how much (off-reference fraction) and which way
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
        annotate_grid(axes[0], off, cmap0, norm0, "{:.2f}")

        # Robust cap so a couple of tiny-n outliers don't wash out the bulk contraction-bias signal.
        _absv = np.abs(signed[np.isfinite(signed)])
        _cap = max(float(np.nanpercentile(_absv, 90)) if _absv.size else 2.0, 2.0)
        cmap1 = plt.get_cmap("RdBu")
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
                "### 2 · Stutter shape — how much, and which way\n"
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
def _(
    LEN_LABELS,
    MIN_LOCI,
    Normalize,
    PERIODS,
    TwoSlopeNorm,
    annotate_grid,
    cells,
    empty_grid,
    grid_axes,
    mo,
    np,
    pd,
    plt,
):
    # GRID 3 — where the two aligners diverge, per cell. Two heatmaps: (a) fraction of loci whose
    # complete-length histogram differs, (b) unit_slip minus flat_gap off-reference fraction (does
    # unit_slip keep more stutter reads off the reference, as its design intends?). Cells with fewer
    # than MIN_LOCI loci are left blank so a single-locus cell can't read as "100% divergence".
    def _len_hist(g):
        return tuple(sorted((int(l), int(r)) for l, r in zip(g["obs_len"], g["reads"])))

    def divergence_grids():
        comp = cells[cells["cls"] == "complete"]
        # per-locus complete-length histogram per aligner
        h = (
            comp.groupby(["contig", "start", "end", "period", "len_band", "aligner"], observed=True)
            .apply(_len_hist, include_groups=False)
            .unstack("aligner")
        )
        h = h.reset_index()
        for al in ("flat_gap", "unit_slip"):
            if al not in h.columns:
                h[al] = None
        h["diverge"] = h["flat_gap"] != h["unit_slip"]

        diverge = empty_grid()
        nloci = empty_grid()
        for i, period in enumerate(PERIODS):
            for j, band in enumerate(LEN_LABELS):
                sub = h[(h["period"] == period) & (h["len_band"] == band)]
                if len(sub) < MIN_LOCI:
                    continue
                nloci[i, j] = len(sub)
                diverge[i, j] = sub["diverge"].mean()

        # unit_slip - flat_gap off-reference fraction per cell (reads-weighted)
        def off_ref(al):
            g = comp[comp["aligner"] == al]
            return g.groupby(["period", "len_band"], observed=True).apply(
                lambda x: pd.Series(
                    {"off": x.loc[x["stutter"] != 0, "reads"].sum(), "tot": x["reads"].sum()}
                ),
                include_groups=False,
            )

        u = off_ref("unit_slip")
        f = off_ref("flat_gap")
        delta = empty_grid()
        for i, period in enumerate(PERIODS):
            for j, band in enumerate(LEN_LABELS):
                key = (period, band)
                if key in u.index and key in f.index and u.loc[key, "tot"] and f.loc[key, "tot"]:
                    delta[i, j] = (
                        u.loc[key, "off"] / u.loc[key, "tot"]
                        - f.loc[key, "off"] / f.loc[key, "tot"]
                    )
        # blank the delta panel where the divergence panel is too sparse, so the two agree on scope
        delta[~np.isfinite(diverge)] = np.nan
        return diverge, nloci, delta

    def divergence_figure():
        diverge, nloci, delta = divergence_grids()
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.3))

        cmap0 = plt.get_cmap("Purples")
        norm0 = Normalize(vmin=0, vmax=1.0)
        im0 = axes[0].imshow(diverge, cmap=cmap0, norm=norm0, aspect="auto")
        axes[0].set_title("loci with a different complete-length histogram", fontsize=10)
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.03)
        annotate_grid(axes[0], diverge, cmap0, norm0, "{:.0%}", ncnt=nloci)

        _lim = max(float(np.nanmax(np.abs(delta))) if np.isfinite(delta).any() else 0.05, 0.05)
        cmap1 = plt.get_cmap("RdBu_r")
        norm1 = TwoSlopeNorm(vcenter=0, vmin=-_lim, vmax=_lim)
        im1 = axes[1].imshow(delta, cmap=cmap1, norm=norm1, aspect="auto")
        axes[1].set_title(
            "off-ref frac: unit_slip − flat_gap  (red = unit_slip keeps more stutter)", fontsize=9.5
        )
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.03)
        annotate_grid(axes[1], delta, cmap1, norm1, "{:+.2f}")

        for ax in axes:
            grid_axes(ax)
        fig.suptitle("Aligner divergence by period × tract length", fontweight="bold")
        fig.tight_layout()
        return fig

    mo.vstack(
        [
            mo.md(
                "### 3 · Where the aligners diverge\n"
                "Left: fraction of covered loci (with a complete read under both) whose complete "
                "**length histogram differs** between the two delimiters — the cells where the "
                "choice changes calls. Right: unit_slip **−** flat_gap off-reference fraction; "
                "**red = unit_slip keeps more reads off the reference** (its cheaper-slip design "
                "declining to collapse stutter onto REF), blue = the reverse. Cells with fewer than "
                f"**{MIN_LOCI}** loci are blank."
            ),
            divergence_figure(),
        ]
    )
    return


@app.cell
def _(LEN_LABELS, PERIODS, mo, period_label):
    detail_period = mo.ui.dropdown(
        {period_label(p): p for p in PERIODS}, value=period_label(2), label="period"
    )
    detail_band = mo.ui.dropdown(LEN_LABELS, value="20-29", label="tract length (bp)")
    mo.md("### 4 · Drill into one cell — the allele-length distribution")
    return detail_band, detail_period


@app.cell
def _(detail_band, detail_period, mo):
    mo.hstack([detail_period, detail_band], justify="start", gap=2)
    return


@app.cell
def _(
    FLAT_GAP_C,
    UNIT_SLIP_C,
    cells,
    detail_band,
    detail_period,
    mo,
    np,
    plt,
):
    # The per-cell "shape": reads-weighted complete-stutter distributions, both aligners overlaid,
    # plus the partial / no-border tallies the frontier grids summarise.
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

        fig, ax = plt.subplots(figsize=(10, 3.6))
        for al, color in (("unit_slip", UNIT_SLIP_C), ("flat_gap", FLAT_GAP_C)):
            g = comp[comp["aligner"] == al]
            hist, _ = np.histogram(g["stutter"], bins=bins, weights=g["reads"])
            frac = hist / hist.sum() if hist.sum() else hist
            ax.step(centers, frac, where="mid", color=color, linewidth=2, label=al)
            ax.fill_between(centers, frac, step="mid", color=color, alpha=0.12)
        ax.axvline(0, color="#b8b6ad", linewidth=1, zorder=0)
        ax.set_xlabel("stutter Δ = observed − reference (bp)")
        ax.set_ylabel("fraction of complete reads")
        ax.set_title(f"Complete-read stutter — period {period}, tract {band} bp", fontsize=11)
        ax.legend(frameon=False)
        ax.grid(True, axis="y", alpha=0.25)
        fig.tight_layout()

        # outcome tally caption
        tally = {}
        for al in ("unit_slip", "flat_gap"):
            g = sub[sub["aligner"] == al]
            tally[al] = {
                c: int(g.loc[g["cls"] == c, "reads"].sum())
                for c in ("complete", "partial", "no_border")
            }
        n_loci = sub[["contig", "start", "end"]].drop_duplicates().shape[0]
        cap = (
            f"**{n_loci} loci** in this cell. Reads by outcome — "
            f"unit_slip: {tally['unit_slip']}, flat_gap: {tally['flat_gap']}."
        )
        return fig, cap

    _fig, _cap = detail_figure(detail_period.value, detail_band.value)
    mo.vstack([_fig, mo.md(_cap)])
    return


@app.cell
def _(mo, n_loci, tsv_path):
    mo.md(
        f"""
        ---
        *Source: `{tsv_path.name}` — {n_loci:,} covered loci, both delimiters, one region-typing
        walk (`examples/ng_ssr_aligner_bakeoff.rs`). Periods 1–6 (region typing emits homopolymers
        too, which §4.2 of the alignment spec deliberately keeps in the comparison). Partials are
        recorded but never fed to a length model — they are censored lower bounds.*
        """
    )
    return


if __name__ == "__main__":
    app.run()
