# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "numpy", "matplotlib"]
# ///
"""Per-sample STR stutter across the tomato cohort — how strong, and from what tract length.

Two questions, both per sample:

  1. **At which tract length does stuttering start to matter?** Stutter grows with the length of
     the repeat, and does so at a different rate for each motif period, so the onset is read off a
     length axis faceted by period — not from a single pooled number.
  2. **How strong is the stuttering in each sample?** Library chemistry differs between samples;
     PCR amplification before sequencing adds slippage that a PCR-free library does not have. If
     that is present in this cohort, samples should separate on stutter strength.

The measure that answers both is **off-mode** stutter, not off-reference. A tomato accession is not
the reference: a read whose tract differs from the reference may simply be carrying a real allele.
So for each (sample, locus) the sample's own **modal complete length** stands in for its genotype,
and stutter is the spread of its reads around that mode. Off-reference is kept alongside as the
naive comparison, and the gap between the two is itself informative.

Caveat that cannot be engineered away here: at a heterozygous locus the second allele's reads count
as off-mode. These accessions are largely inbred, so most loci are homozygous, but the off-mode rate
is an upper bound on stutter.

Data comes from the cohort dump (one region-typing walk, every sample, ng's default delimiter):

    ./scripts/dev.sh cargo run --release --example ng_ssr_cohort_stutter -- \\
        --regions benchmarks/ssr_tomato1/ssr_regions.bed \\
        $HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \\
        benchmarks/ssr_tomato1/crams/*.bench.cram > tmp/tomato_stutter/cohort.tsv

plus a sample table (`sample_meta.tsv`) carrying each sample's duplicate rate, which is used as an
independent, data-side proxy for PCR amplification — the CRAM headers carry no library-prep flag.

Override the inputs with PVC_STUTTER_TSV / PVC_STUTTER_META.

Run:  uv run marimo run  benchmarks/ssr_tomato1/scripts/ng_ssr_cohort_stutter_dashboard.py
Edit: uv run marimo edit benchmarks/ssr_tomato1/scripts/ng_ssr_cohort_stutter_dashboard.py
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
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    # dataviz — one sequential ramp (the documented blue scale) because the sample axis is encoded
    # by a *continuous* covariate, the duplicate rate. 51 samples cannot take categorical hues:
    # the rule is a fixed order of at most eight, never cycled, so a per-sample palette is not on
    # the table. Recessive grey carries "a sample, unhighlighted".
    BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
    DUP_CMAP = LinearSegmentedColormap.from_list("dup", BLUES)
    GREY = "#b8b6ad"
    INK = "#52514e"
    ACCENT = "#eb6834"  # categorical slot 2 — reserved for the cohort aggregate line

    # A locus needs enough complete reads for its modal length to mean anything; below this the
    # "mode" is one read and every other read looks like stutter.
    MIN_LOCUS_DEPTH = 6
    # Cells thinner than this are left out of the per-sample curves rather than drawn as noise.
    MIN_CELL_READS = 30
    return (
        ACCENT,
        DUP_CMAP,
        GREY,
        INK,
        MIN_CELL_READS,
        MIN_LOCUS_DEPTH,
        Normalize,
        Path,
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
        # Tomato cohort — per-sample STR stutter

        Every sample delimited at the **same** microsatellite tracts, by ng's default delimiter
        (`SsrUnitRobustAligner`), in one region-typing walk. Two questions:

        1. **From what tract length does stuttering start to matter**, per motif period?
        2. **How strong is each sample's stuttering** — and do samples separate the way a
           PCR / PCR-free split would predict?

        **Stutter is measured off-mode, not off-reference.** These accessions are not the
        reference, so a tract that differs from the reference may be a real allele. For each
        (sample, locus) the sample's own modal complete length stands in for its genotype, and
        stutter is the spread of its reads around that mode. At a heterozygous locus the second
        allele counts as off-mode, so this is an **upper bound** — these accessions are largely
        inbred, which keeps the bound tight.
        """
    )
    return


@app.cell
def _(Path, mo, os, pd):
    _here = Path(__file__).resolve()
    _root = _here.parents[3]
    _default_tsv = _root / "tmp" / "tomato_stutter" / "cohort.tsv"
    _default_meta = _root / "tmp" / "tomato_stutter" / "sample_meta.tsv"
    tsv_path = Path(os.environ.get("PVC_STUTTER_TSV") or _default_tsv)
    meta_path = Path(os.environ.get("PVC_STUTTER_META") or _default_meta)

    mo.stop(
        not tsv_path.is_file(),
        mo.md(
            f"**Missing** `{tsv_path}`. Generate it:\n\n"
            "```\n./scripts/dev.sh cargo run --release --example ng_ssr_cohort_stutter -- \\\n"
            "  --regions benchmarks/ssr_tomato1/ssr_regions.bed \\\n"
            "  $HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \\\n"
            "  benchmarks/ssr_tomato1/crams/*.bench.cram > tmp/tomato_stutter/cohort.tsv\n```\n\n"
            "or set `PVC_STUTTER_TSV`."
        ),
    )

    # A full-cohort dump is millions of rows, and every analysis below is length-based — the tract
    # sequences are needed only to measure them. So the strings are turned into lengths and then
    # dropped, and the repeated identifiers are held as categories; without this the frame is
    # several gigabytes of text that is never read again.
    # The `#rg` header is the run's read-group table: the rows carry only a numeric id, and this is
    # what resolves it to sample / library / experiment. Reading it first also means the dump can
    # gain read groups without the loader below caring.
    _rg_rows = []
    with open(tsv_path) as _fh:
        for _line in _fh:
            if not _line.startswith("#"):
                break
            if _line.startswith("#rg\t"):
                _rg_rows.append(_line.rstrip("\n").split("\t")[1:])
    rg_table = pd.DataFrame(
        _rg_rows,
        columns=[
            "read_group", "rg_id", "sample", "library", "library_origin",
            "experiment", "experiment_origin", "platform", "file",
        ],
    )
    if not rg_table.empty:
        rg_table["read_group"] = rg_table["read_group"].astype("int32")

    raw = pd.read_csv(
        tsv_path,
        sep="\t",
        comment="#",
        dtype={
            "observed": str,
            "motif": str,
            "ref_tract": str,
            "sample": "category",
            "contig": "category",
            "coverage": "category",
        },
    )
    raw["observed"] = raw["observed"].fillna("")
    raw["period"] = raw["motif"].str.len().astype("int16")
    raw["ref_len"] = raw["ref_tract"].str.len().astype("int32")
    raw["obs_len"] = raw["observed"].str.len().astype("int32")
    raw = raw.drop(columns=["motif", "ref_tract", "observed"])

    meta = (
        pd.read_csv(meta_path, sep="\t")
        if meta_path.is_file()
        else pd.DataFrame(columns=["sample", "run", "project", "reads", "duplicates", "dup_pct"])
    )
    return meta, meta_path, raw, rg_table, tsv_path


@app.cell
def _(mo, rg_table):
    # The grain to analyse at. Chemistry is a property of the library preparation, so the sample is
    # the wrong default whenever one sample holds several libraries — and in this archive one holds
    # sixteen. The dump stores the finest grain (the read group) precisely so this stays a choice.
    _choices = {
        "library — one preparation (recommended)": "library",
        "experiment — preparation + sequencing config": "experiment",
        "read group — one @RG, usually a lane": "read_group",
        "sample — the individual (merges libraries)": "sample",
    }
    unit_sel = mo.ui.radio(
        _choices, value="library — one preparation (recommended)", label="Analysis unit"
    )
    _have_rg = not rg_table.empty
    mo.vstack(
        [
            mo.md(
                "Chemistry belongs to the **library**, not to the individual. Where a sample holds "
                "more than one library, folding to the sample averages across preparations — which "
                "is the very thing a stutter comparison is trying to separate. Pick the grain:"
                if _have_rg
                else "*This dump predates the read-group table, so only the sample grain is "
                "available. Regenerate it with the current `ng_ssr_cohort_stutter` to fold by "
                "library or experiment.*"
            ),
            unit_sel if _have_rg else mo.md(""),
        ]
    )
    return (unit_sel,)


@app.cell
def _(MIN_LOCUS_DEPTH, np, pd, raw, rg_table, unit_sel):
    # The one derived frame everything below shares: complete reads only, each tagged with its
    # unit's modal length at that locus, and with both stutter measures.
    #
    # Only complete reads carry a length at all — a partial is a censored lower bound and cannot
    # enter a spread. Loci with fewer than MIN_LOCUS_DEPTH complete reads are dropped because their
    # "mode" would be a single read, against which every other read scores as stutter.
    comp = raw[raw["coverage"] == "complete"].copy()

    # `unit` is the analysis grain. It is resolved per row from the read-group table, so the modal
    # length, the depth filter and every statistic below are computed within one chemistry rather
    # than across a mixture of them.
    if rg_table.empty or unit_sel.value == "sample":
        comp["unit"] = comp["sample"].astype(str)
    else:
        _map = rg_table.set_index("read_group")[unit_sel.value]
        comp["unit"] = comp["read_group"].map(_map).astype(str)
    comp["unit"] = comp["unit"].astype("category")

    _locus = ["unit", "contig", "start", "end"]

    _depth = comp.groupby(_locus, observed=True)["reads"].transform("sum")
    comp = comp[_depth >= MIN_LOCUS_DEPTH].copy()
    comp["locus_depth"] = comp.groupby(_locus, observed=True)["reads"].transform("sum")

    # Reads sharing a length must be pooled before anything asks "which length is heaviest". The
    # dump emits one row per distinct observed *sequence*, and two sequences can have the same
    # length — an interrupted tract and a pure one of equal size — so a length can arrive split
    # across rows. Keying on rows instead of lengths silently under-counts the true allele and
    # inflates stutter.
    comp["len_reads"] = comp.groupby([*_locus, "obs_len"], observed=True)["reads"].transform("sum")

    # Modal complete length per (sample, locus): the heaviest length, ties broken by the shorter
    # one so the answer does not depend on row order.
    _ordered = comp.sort_values(
        [*_locus, "len_reads", "obs_len"], ascending=[True] * 4 + [False, True]
    )
    _modal = (
        _ordered.groupby(_locus, observed=True)
        .first()
        .rename(columns={"obs_len": "modal_len"})[["modal_len"]]
    )
    comp = comp.merge(_modal, left_on=_locus, right_index=True, how="left")

    comp["off_mode"] = comp["obs_len"] - comp["modal_len"]
    comp["off_ref"] = comp["obs_len"] - comp["ref_len"]

    # **Leave-one-out off-mode.** The plain off-mode rate is biased *downward*, and the bias
    # depends on depth: the mode is estimated from the very reads it is then compared against, so
    # each read votes for its own agreement. At 7 reads a locus that bias is large, and because
    # depth varies between samples it would manufacture exactly the between-sample differences
    # this notebook exists to measure.
    #
    # The fix is to score every read against the mode of the *other* reads at its locus. On the
    # aggregated counts that is exact and cheap: a read of length L at a locus whose length counts
    # are `n` is off-mode iff L is not the argmax of `n` with one L removed. Only the reads of the
    # modal length can change the answer, and only when the mode's lead is a single read.
    # Done with a per-locus top-two rather than a Python function per group: at full cohort scope
    # this runs over millions of loci, where a groupby-apply would take longer than the dump did.
    #
    # Removing one read of length L changes only L's own count, so the winner among the *other*
    # lengths is the locus's top length — or its runner-up, when L is the top one. A read agrees
    # with the leave-one-out mode iff its own reduced count still beats that, ties going to the
    # shorter length (the same rule the modal length above uses).
    _lens = (
        comp[[*_locus, "obs_len", "len_reads"]]
        .drop_duplicates([*_locus, "obs_len"])
        .sort_values([*_locus, "len_reads", "obs_len"], ascending=[True] * 4 + [False, True])
    )
    _lens["_rank"] = _lens.groupby(_locus, observed=True).cumcount()
    _top1 = _lens[_lens["_rank"] == 0].rename(
        columns={"obs_len": "t1_len", "len_reads": "t1_n"}
    )[[*_locus, "t1_len", "t1_n"]]
    _top2 = _lens[_lens["_rank"] == 1].rename(
        columns={"obs_len": "t2_len", "len_reads": "t2_n"}
    )[[*_locus, "t2_len", "t2_n"]]
    comp = comp.merge(_top1, on=_locus, how="left").merge(_top2, on=_locus, how="left")

    _is_top = comp["obs_len"] == comp["t1_len"]
    # A locus with one distinct length has no runner-up: -1 loses to any count, so its reads agree.
    _other_n = np.where(_is_top, comp["t2_n"].fillna(-1), comp["t1_n"])
    _other_len = np.where(_is_top, comp["t2_len"].fillna(np.inf), comp["t1_len"])
    _mine = comp["len_reads"] - 1
    _agrees = (_mine > _other_n) | ((_mine == _other_n) & (comp["obs_len"] < _other_len))
    comp["loo_off_reads"] = np.where(_agrees, 0, comp["reads"]).astype(int)
    comp = comp.drop(columns=["t1_len", "t1_n", "t2_len", "t2_n"])

    # Length bands against ~150 bp reads, matching the human bake-off dashboard's axis.
    LEN_EDGES = [0, 10, 15, 20, 30, 40, 60, 10**9]
    LEN_LABELS = ["<10", "10-14", "15-19", "20-29", "30-39", "40-59", "60+"]
    comp["len_band"] = pd.cut(
        comp["ref_len"], bins=LEN_EDGES, labels=LEN_LABELS, right=False, ordered=True
    )
    # The same tract measured in repeat UNITS rather than bases. Slippage is a per-unit process, so
    # copy number is the axis the mechanism would predict, while base length is the axis the read
    # geometry cares about — and the two orderings differ, because 20 bp is ten dinucleotide copies
    # but only three hexamer ones. Which of them stutter actually tracks is the question.
    comp["ref_copies"] = comp["ref_len"] / comp["period"]
    COPY_EDGES = [3, 4, 5, 6, 7, 9, 12, 16, 25, 10**9]
    COPY_LABELS = ["3", "4", "5", "6", "7-8", "9-11", "12-15", "16-24", "25+"]
    comp["copy_band"] = pd.cut(
        comp["ref_copies"], bins=COPY_EDGES, labels=COPY_LABELS, right=False, ordered=True
    )
    PERIOD_NAME = {1: "mono", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
    PERIODS = [1, 2, 3, 4, 5, 6]
    return COPY_LABELS, LEN_LABELS, PERIOD_NAME, PERIODS, comp


@app.cell
def _(LEN_LABELS, PERIODS, np, plt):
    # A period × tract-length grid, the axis this whole question lives on. Shared with the human
    # bake-off dashboard deliberately: the same two dimensions, read the same way, so a number seen
    # in one is comparable with a number seen in the other.
    def empty_grid(cols=None):
        return np.full((len(PERIODS), len(cols if cols is not None else LEN_LABELS)), np.nan)

    def annotate_grid(ax, grid, cmap, norm, fmt, ncnt=None):
        """Write each finite cell's value, picking white or dark text by the cell's luminance so
        the number stays legible on pale and saturated fills alike."""
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                v = grid[i, j]
                if not np.isfinite(v):
                    continue
                r, g, b, _ = cmap(norm(v))
                lum = 0.299 * r + 0.587 * g + 0.114 * b
                txt = fmt.format(v)
                if ncnt is not None and np.isfinite(ncnt[i, j]):
                    txt = f"{txt}\nn={int(ncnt[i, j]):,}"
                ax.text(
                    j, i, txt, ha="center", va="center", fontsize=6.2,
                    color="white" if lum < 0.5 else "#222",
                )

    def grid_axes(ax, period_label, labels=None, xlabel="reference tract length (bp)"):
        labels = labels if labels is not None else LEN_LABELS
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7.5)
        ax.set_yticks(range(len(PERIODS)))
        ax.set_yticklabels([period_label(p) for p in PERIODS], fontsize=8)
        ax.set_xlabel(xlabel)

    return annotate_grid, empty_grid, grid_axes


@app.cell
def _(MIN_LOCUS_DEPTH, comp, meta, mo, pd, raw):
    # Per-unit headline: how much data each unit contributes and how much it stutters overall.
    def _per_sample():
        g = comp.groupby("unit", observed=True)
        out = pd.DataFrame(
            {
                "loci": g.apply(
                    lambda x: x[["contig", "start", "end"]].drop_duplicates().shape[0],
                    include_groups=False,
                ),
                "complete_reads": g["reads"].sum(),
                "mean_depth": g.apply(
                    lambda x: x.groupby(["contig", "start", "end"], observed=True)["reads"]
                    .sum()
                    .mean(),
                    include_groups=False,
                ),
                "off_mode": g.apply(
                    lambda x: x.loc[x["off_mode"] != 0, "reads"].sum() / x["reads"].sum(),
                    include_groups=False,
                ),
                "off_mode_loo": g.apply(
                    lambda x: x["loo_off_reads"].sum() / x["reads"].sum(),
                    include_groups=False,
                ),
                "off_ref": g.apply(
                    lambda x: x.loc[x["off_ref"] != 0, "reads"].sum() / x["reads"].sum(),
                    include_groups=False,
                ),
            }
        ).reset_index()
        # The duplicate-rate metadata is keyed by sample. A unit belongs to exactly one sample
        # (a library is prepared from one individual), so carrying that sample across is a lookup,
        # never an aggregation.
        owner = comp.groupby("unit", observed=True)["sample"].first().astype(str)
        out["sample"] = out["unit"].map(owner)
        if not meta.empty:
            out = out.merge(meta[["sample", "run", "dup_pct"]], on="sample", how="left")
        return out.sort_values("off_mode", ascending=False)

    per_sample = _per_sample()

    _n = per_sample.shape[0]
    _lo, _hi = per_sample["off_mode"].min(), per_sample["off_mode"].max()
    _show = per_sample.copy()
    for _c in ("off_mode", "off_mode_loo", "off_ref"):
        _show[_c] = _show[_c].map(lambda v: f"{v:.2%}")
    _show["complete_reads"] = _show["complete_reads"].map(lambda v: f"{int(v):,}")
    _show["loci"] = _show["loci"].map(lambda v: f"{int(v):,}")
    _show["mean_depth"] = _show["mean_depth"].map(lambda v: f"{v:.1f}")

    mo.vstack(
        [
            mo.md(
                f"**{_n} samples**, {raw['sample'].nunique()} in the dump, over "
                f"{comp[['contig', 'start', 'end']].drop_duplicates().shape[0]:,} distinct loci "
                f"with at least {MIN_LOCUS_DEPTH} complete reads. Off-mode stutter spans "
                f"**{_lo:.2%} to {_hi:.2%}** across samples — a "
                f"**{_hi / max(_lo, 1e-9):.1f}×** range. `dup_pct` is the duplicate rate, the "
                f"data-side proxy for PCR amplification (the CRAM headers carry no library-prep "
                f"flag).\n\n"
                f"`off_mode` is biased **downward** at low depth — the mode is estimated from the "
                f"same reads it is compared against, so each read votes for its own agreement, and "
                f"the bias shrinks as depth grows. Since depth differs between samples that bias "
                f"could fake a between-sample difference, so `off_mode_loo` scores every read "
                f"against the mode of the *other* reads at its locus, which removes it. Read "
                f"`off_mode_loo` against `mean_depth`: if the two track each other, depth is still "
                f"talking. (At a heterozygous locus the leave-one-out measure over-counts, since "
                f"each allele's reads disagree with the other's — an upper bound, tightened by "
                f"these accessions being largely inbred.)"
            ),
            mo.ui.table(_show, selection=None, pagination=True, page_size=12),
        ]
    )
    return (per_sample,)


@app.cell
def _(ACCENT, GREY, INK, mo, np, per_sample, plt):
    # SECTION 1 — how strong, per sample. A ranked dot plot: 51 samples is far past what a
    # categorical palette can carry, and the question is "where does this sample sit in the
    # cohort", which rank + position answers directly.
    def strength_figure():
        d = per_sample.sort_values("off_mode")
        y = np.arange(len(d))
        fig, ax = plt.subplots(figsize=(9, max(3.5, 0.17 * len(d))))
        ax.hlines(y, 0, d["off_mode"], color=GREY, linewidth=1.2, zorder=1)
        ax.scatter(d["off_mode"], y, s=34, color="#2a78d6", zorder=2, label="off-mode (stutter)")
        ax.scatter(
            d["off_ref"], y, s=22, facecolors="none", edgecolors=ACCENT, linewidth=1.2, zorder=3,
            label="off-reference (stutter + real alleles)",
        )
        ax.set_yticks(y)
        ax.set_yticklabels(d["unit"], fontsize=6)
        ax.set_xlabel("fraction of complete reads away from the sample's modal length")
        ax.xaxis.set_major_formatter(lambda v, _p: f"{v:.0%}")
        ax.grid(True, axis="x", alpha=0.25)
        ax.set_axisbelow(True)
        ax.margins(y=0.01)
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        med = float(d["off_mode"].median())
        ax.axvline(med, color=INK, linewidth=1, linestyle="--", zorder=0)
        ax.text(med, len(d) - 0.5, f" cohort median {med:.1%}", fontsize=7.5, color=INK, va="top")
        fig.suptitle("Stutter strength per sample", fontweight="bold")
        fig.tight_layout()
        return fig

    mo.vstack(
        [
            mo.md(
                "## 1 · How strong is the stuttering, per sample\n"
                "Filled = **off-mode**, the stutter measure. Hollow = off-reference, which adds "
                "the sample's genuine differences from the reference genome and is therefore "
                "always larger; the gap between them is roughly how far that accession sits from "
                "the reference. A PCR-amplified library should sit high on the *filled* axis."
            ),
            strength_figure(),
        ]
    )
    return


@app.cell
def _(comp, mo, np, pd):
    # Is the spread in that chart real? Each sample contributes k off-mode reads out of n, so if
    # every sample shared one true rate the k's would be binomial(n, p). A range means nothing on
    # its own — with a few thousand reads per sample the extremes move a long way on noise alone.
    # This is the test that separates "samples differ" from "small n".
    def overdispersion(sub):
        g = sub.groupby("unit", observed=True)
        tab = pd.DataFrame(
            {
                "n": g["reads"].sum(),
                "k": g.apply(
                    lambda x: x.loc[x["off_mode"] != 0, "reads"].sum(), include_groups=False
                ),
            }
        )
        tab = tab[tab["n"] > 0]
        if len(tab) < 3:
            return None
        p = tab["k"].sum() / tab["n"].sum()
        rate = tab["k"] / tab["n"]
        expected = tab["n"] * p
        chi2 = float((((tab["k"] - expected) ** 2) / (expected * (1 - p))).sum())
        dof = len(tab) - 1
        sigma = np.sqrt(p * (1 - p) / tab["n"])
        return {
            "samples": len(tab),
            "median n": f"{tab['n'].median():,.0f}",
            "pooled rate": f"{p:.2%}",
            "observed spread": f"{rate.min():.2%} – {rate.max():.2%}",
            "χ²/dof": f"{chi2 / max(dof, 1):.2f}",
            "variance vs binomial": f"{rate.var(ddof=1) / (p * (1 - p) / tab['n']).mean():.1f}×",
            "outside ±3σ": f"{int((np.abs(rate - p) > 3 * sigma).sum())} / {len(tab)}",
        }

    _strata = {
        "all loci": comp,
        "mononucleotide": comp[comp["period"] == 1],
        "dinucleotide ≥15 bp": comp[(comp["period"] == 2) & (comp["ref_len"] >= 15)],
    }
    _rows = [{"stratum": k, **v} for k, v in
             ((k, overdispersion(v)) for k, v in _strata.items()) if v]

    mo.vstack(
        [
            mo.md(
                "### Are those differences real, or counting noise?\n"
                "**χ²/dof = 1.0 is what pure noise looks like.** Above it, samples genuinely "
                "differ; at it, the range in the chart above is an artefact of how few reads each "
                "sample contributes. Read each stratum separately — a stratum can be real overall "
                "and powerless once split by period and length."
            ),
            mo.ui.table(pd.DataFrame(_rows), selection=None, pagination=False),
        ]
    )
    return


@app.cell
def _(
    COPY_LABELS,
    LEN_LABELS,
    MIN_CELL_READS,
    Normalize,
    PERIODS,
    PERIOD_NAME,
    annotate_grid,
    comp,
    empty_grid,
    grid_axes,
    mo,
    np,
    plt,
):
    # SECTION 2 — the onset, as a period × size grid, drawn on both size axes. This is the shape of
    # the answer: a single number per cell, cohort-pooled, read straight off two axes. The
    # per-sample curves that follow answer a different question (do samples agree?) and are much
    # harder to read for this one.
    #
    # Both axes share this one implementation, and — this is the point of drawing both — share one
    # colour scale, so a cell's shade means the same thing in each and the two fronts are directly
    # comparable rather than each normalised to its own maximum.
    def onset_grid_figure(column, labels, xlabel, title, norm=None):
        off = empty_grid(labels)
        n = empty_grid(labels)
        for i, period in enumerate(PERIODS):
            for j, band in enumerate(labels):
                g = comp[(comp["period"] == period) & (comp[column] == band)]
                tot = g["reads"].sum()
                if tot < MIN_CELL_READS:
                    continue
                n[i, j] = tot
                off[i, j] = g.loc[g["off_mode"] != 0, "reads"].sum() / tot

        fig, ax = plt.subplots(figsize=(max(9.5, 1.15 * len(labels)), 4.2))
        cmap = plt.get_cmap("Blues")
        if norm is None:
            # Capped at the 95th percentile so one saturated cell does not flatten the gradient the
            # onset is read from; the annotations carry the true value regardless.
            finite = off[np.isfinite(off)]
            vmax = max(float(np.nanpercentile(finite, 95)) if finite.size else 0.2, 0.05)
            norm = Normalize(vmin=0, vmax=vmax)
        im = ax.imshow(np.clip(off, 0, norm.vmax), cmap=cmap, norm=norm, aspect="auto")
        fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02).set_label(
            "off-mode read fraction", fontsize=8.5
        )
        annotate_grid(ax, off, cmap, norm, "{:.1%}", ncnt=n)
        grid_axes(ax, lambda p: f"{PERIOD_NAME.get(p, p)} ({p})", labels=labels, xlabel=xlabel)
        ax.set_title(title, fontweight="bold", fontsize=11)
        fig.tight_layout()
        return fig, norm

    _len_fig, _shared_norm = onset_grid_figure(
        "len_band",
        LEN_LABELS,
        "reference tract length (bp)",
        "Stutter by period × tract LENGTH — cohort-pooled off-mode fraction",
    )
    _copy_fig, _ = onset_grid_figure(
        "copy_band",
        COPY_LABELS,
        "reference tract length (repeat units)",
        "Stutter by period × REPEAT COUNT — same reads, same colour scale",
        norm=_shared_norm,
    )

    mo.vstack(
        [
            mo.md(
                "## 2 · From what size does stuttering start to matter\n"
                "One cell per (period, size), pooled over the cohort: the fraction of complete "
                "reads sitting away from their own unit's modal length. Read the **onset** as the "
                "column where a row stops being pale. Cells with fewer than "
                f"**{MIN_CELL_READS}** complete reads are left blank.\n\n"
                "**Two size axes, because they are different claims about the mechanism.** "
                "Base length is what the *read* cares about — a tract competes with the read for "
                "room to be spanned. Repeat count is what the *polymerase* cares about — slippage "
                "happens per unit, so a run of ten copies offers ten chances to slip whether the "
                "unit is 2 bp or 6 bp. If stutter tracked base length the front would fall in the "
                "same column for every period; if it tracks copy number it should straighten out "
                "on the second grid."
            ),
            _len_fig,
            _copy_fig,
            mo.md(
                "*The two grids hold the same reads and share one colour scale, so shades are "
                "comparable between them. The empty lower-left of the copy grid is not missing "
                "data: the catalog's per-period copy floors (6 for mono, 4 for di and tri, 3 "
                "above) mean short-copy tracts of some periods are never admitted in the first "
                "place.*"
            ),
        ]
    )
    return


@app.cell
def _(
    DUP_CMAP,
    GREY,
    LEN_LABELS,
    MIN_CELL_READS,
    Normalize,
    PERIOD_NAME,
    comp,
    meta,
    mo,
    np,
    per_sample,
    plt,
):
    # SECTION 2b — the same axis, one line per sample. Lines are coloured by duplicate rate so the
    # PCR question is answered by looking: if amplification drives stutter, the dark lines separate
    # upward.
    def onset_curves_figure():
        dup = (
            meta.set_index("sample")["dup_pct"].to_dict() if not meta.empty else {}
        )
        norm = Normalize(
            vmin=min(dup.values()) if dup else 0, vmax=max(dup.values()) if dup else 1
        )
        periods = sorted(comp["period"].unique())[:6]
        fig, axes = plt.subplots(2, 3, figsize=(13, 6.4), sharey=True, sharex=True)
        for ax, period in zip(axes.flat, periods):
            sub = comp[comp["period"] == period]
            for sample, g in sub.groupby("unit", observed=True):
                xs, ys = [], []
                for j, band in enumerate(LEN_LABELS):
                    cell = g[g["len_band"] == band]
                    tot = cell["reads"].sum()
                    if tot < MIN_CELL_READS:
                        continue
                    xs.append(j)
                    ys.append(cell.loc[cell["off_mode"] != 0, "reads"].sum() / tot)
                if len(xs) < 2:
                    continue
                colour = DUP_CMAP(norm(dup[sample])) if sample in dup else GREY
                ax.plot(xs, ys, color=colour, linewidth=1.1, alpha=0.85)
            ax.set_title(f"{PERIOD_NAME.get(period, period)} ({period} bp)", fontsize=10)
            ax.set_xticks(range(len(LEN_LABELS)))
            ax.set_xticklabels(LEN_LABELS, rotation=45, ha="right", fontsize=7.5)
            ax.grid(True, axis="y", alpha=0.25)
            ax.set_axisbelow(True)
        for ax in axes[:, 0]:
            ax.set_ylabel("off-mode read fraction")
        for ax in axes.flat[len(periods) :]:
            ax.axis("off")
        if dup:
            sm = plt.cm.ScalarMappable(cmap=DUP_CMAP, norm=norm)
            cb = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
            cb.set_label("duplicate rate (%) — the PCR proxy", fontsize=8.5)
        fig.suptitle(
            "Stutter onset: off-mode fraction vs reference tract length, one line per sample",
            fontweight="bold",
        )
        return fig

    mo.vstack(
        [
            mo.md(
                "### 2b · The same thing per sample\n"
                "One line per sample, faceted by motif period, coloured by duplicate rate. The grid "
                "above answers *where* stutter starts; this answers *whether samples agree about "
                "it*. Fifty-one lines is a thicket by design — what to read is whether the dark "
                f"(high-duplicate) lines sit above the pale ones. Cells with fewer than "
                f"**{MIN_CELL_READS}** complete reads are skipped."
            ),
            onset_curves_figure(),
        ]
    )
    return


@app.cell
def _(LEN_LABELS, MIN_CELL_READS, PERIOD_NAME, comp, mo, np, pd, plt):
    # SECTION 3 — the onset as a number: the shortest length band whose off-mode fraction crosses a
    # threshold, per sample × period. A distribution of onsets across samples, not a single value,
    # because the whole question is whether samples differ.
    def onset_table(threshold):
        rows = []
        for (sample, period), g in comp.groupby(["unit", "period"], observed=True):
            onset = None
            for band in LEN_LABELS:
                cell = g[g["len_band"] == band]
                tot = cell["reads"].sum()
                if tot < MIN_CELL_READS:
                    continue
                frac = cell.loc[cell["off_mode"] != 0, "reads"].sum() / tot
                if frac >= threshold:
                    onset = band
                    break
            rows.append({"sample": sample, "period": period, "onset": onset})
        return pd.DataFrame(rows)

    def onset_band_figure(threshold):
        t = onset_table(threshold)
        periods = sorted(comp["period"].unique())[:6]
        counts = np.zeros((len(periods), len(LEN_LABELS) + 1))
        for i, period in enumerate(periods):
            sub = t[t["period"] == period]
            for j, band in enumerate(LEN_LABELS):
                counts[i, j] = (sub["onset"] == band).sum()
            counts[i, -1] = sub["onset"].isna().sum()

        labels = [*LEN_LABELS, "never"]
        fig, ax = plt.subplots(figsize=(11, 3.8))
        bottom = np.zeros(len(periods))
        x = np.arange(len(periods))
        ramp = plt.get_cmap("Blues")
        for j, band in enumerate(labels):
            colour = "#d6d5cd" if band == "never" else ramp(0.25 + 0.65 * j / len(LEN_LABELS))
            ax.bar(
                x, counts[:, j], bottom=bottom, color=colour, width=0.72, label=band,
                edgecolor="white", linewidth=0.6,
            )
            bottom += counts[:, j]
        ax.set_xticks(x)
        ax.set_xticklabels([f"{PERIOD_NAME.get(p, p)}\n({p} bp)" for p in periods], fontsize=8.5)
        ax.set_ylabel("samples")
        ax.grid(True, axis="y", alpha=0.25)
        ax.set_axisbelow(True)
        ax.legend(
            frameon=False, fontsize=8, ncol=len(labels), loc="lower center",
            bbox_to_anchor=(0.5, 1.0), title=None,
        )
        fig.suptitle(
            f"Onset length: shortest band whose off-mode fraction reaches {threshold:.0%}",
            fontweight="bold",
            y=1.12,
        )
        fig.tight_layout()
        return fig

    onset_threshold = mo.ui.slider(
        start=0.02, stop=0.30, step=0.02, value=0.10, label="onset threshold"
    )
    return onset_band_figure, onset_threshold


@app.cell
def _(mo, onset_threshold):
    mo.vstack(
        [
            mo.md(
                "## 3 · The onset as a number\n"
                "For each sample and period, the shortest length band whose off-mode fraction "
                "reaches the threshold. A bar is the cohort's distribution of onsets — if samples "
                "all shared one chemistry the bars would be nearly single-coloured; spread means "
                "samples start stuttering at genuinely different lengths."
            ),
            onset_threshold,
        ]
    )
    return


@app.cell
def _(mo, onset_band_figure, onset_threshold):
    mo.vstack([onset_band_figure(onset_threshold.value)])
    return


@app.cell
def _(DUP_CMAP, INK, Normalize, meta, mo, np, per_sample, plt):
    # SECTION 4 — the PCR hypothesis, tested rather than assumed. If amplification drives stutter,
    # duplicate rate and stutter strength should rise together.
    def pcr_figure():
        if meta.empty or "dup_pct" not in per_sample.columns:
            fig, ax = plt.subplots(figsize=(7, 3))
            ax.text(0.5, 0.5, "no sample metadata loaded", ha="center", va="center")
            ax.axis("off")
            return fig, "—"
        d = per_sample.dropna(subset=["dup_pct"])
        x = d["dup_pct"].to_numpy(dtype=float)
        y = d["off_mode"].to_numpy(dtype=float)
        norm = Normalize(vmin=x.min(), vmax=x.max())

        fig, ax = plt.subplots(figsize=(8, 4.2))
        ax.scatter(x, y, s=52, c=[DUP_CMAP(norm(v)) for v in x], edgecolors="white", linewidth=0.8)
        # Spearman, computed from ranks — the relationship need not be linear, and four
        # high-duplicate outliers would dominate a Pearson r.
        rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
        rho = float(np.corrcoef(rx, ry)[0, 1]) if len(x) > 2 else float("nan")
        if len(x) > 2:
            fit = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 50)
            ax.plot(xs, np.polyval(fit, xs), color=INK, linewidth=1.2, linestyle="--", zorder=0)
        ax.set_xlabel("duplicate rate (%) — PCR proxy")
        ax.set_ylabel("off-mode stutter fraction")
        ax.yaxis.set_major_formatter(lambda v, _p: f"{v:.1%}")
        ax.grid(True, alpha=0.25)
        ax.set_axisbelow(True)
        fig.suptitle("Does the PCR proxy predict stutter?", fontweight="bold")
        fig.tight_layout()
        verdict = (
            f"Spearman **ρ = {rho:+.2f}** over {len(x)} samples. "
            + (
                "A positive rank correlation is what the PCR hypothesis predicts — but duplicate "
                "rate also rises with sequencing depth and falls with library complexity, so this "
                "is consistent with, not proof of, an amplification split."
                if rho > 0.3
                else "That is too weak to support an amplification split on this evidence; "
                "whatever separates these samples on stutter, the duplicate rate does not track "
                "it. A real library-prep label would settle it."
            )
        )
        return fig, verdict

    _fig, _verdict = pcr_figure()
    mo.vstack(
        [
            mo.md(
                "## 4 · Is it PCR?\n"
                "The CRAM headers carry no library-prep flag, so the hypothesis is tested against "
                "a data-side proxy: PCR-amplified libraries duplicate more."
            ),
            _fig,
            mo.md(_verdict),
        ]
    )
    return


@app.cell
def _(GREY, INK, PERIOD_NAME, comp, mo, np, plt):
    # SECTION 5 — the shape of the stutter, which is what a read model actually needs. Two claims
    # the STR model rests on are testable here: that slippage moves the tract by WHOLE motif units,
    # and that the step distribution falls away from ±1.
    #
    # The **0 bar is included**, so each panel is the whole distribution of where reads land rather
    # than the off-allele part alone: how much stutter there IS is then read off the same picture as
    # what shape it has. That forces a log axis — agreement is ~99% and the far steps are near
    # 1-in-10,000, four orders of magnitude a linear axis would flatten to a single spike and a row
    # of nothing.
    def shape_figure():
        # Denominator: reads that either match the allele or sit a whole number of units from it.
        # Reads at a non-unit offset are excluded because they are not slippage at all — they are
        # the residue the panel above measures, and mixing them in would make this neither a
        # slippage distribution nor a complete one.
        d = comp.copy()
        d["whole_unit"] = (d["off_mode"] == 0) | (d["off_mode"] % d["period"] == 0)
        d = d[d["whole_unit"]].copy()
        d["units"] = (d["off_mode"] // d["period"]).astype(int)

        periods = [p for p in sorted(d["period"].unique()) if p <= 6]
        fig, axes = plt.subplots(2, 3, figsize=(13, 6.4), sharex=True, sharey=True)
        steps = list(range(-4, 5))
        for ax, period in zip(axes.flat, periods):
            g = d[d["period"] == period]
            tot = g["reads"].sum()
            if tot <= 0:
                ax.axis("off")
                continue
            frac = [g.loc[g["units"] == u, "reads"].sum() / tot for u in steps]
            # Grey for "matches the allele", blue for contractions, orange for expansions — the same
            # sign convention the rest of this work uses, so a reader never re-learns which is which.
            colours = [GREY if u == 0 else ("#2a78d6" if u < 0 else "#eb6834") for u in steps]
            ax.bar(range(len(steps)), frac, color=colours, width=0.74)
            on_allele = frac[steps.index(0)]
            ax.set_title(
                f"{PERIOD_NAME.get(period, period)} ({period} bp) — "
                f"{1 - on_allele:.2%} of reads stutter",
                fontsize=9.5,
            )
            ax.set_xticks(range(len(steps)))
            ax.set_xticklabels([("0" if u == 0 else f"{u:+d}") for u in steps], fontsize=8)
            ax.grid(True, axis="y", alpha=0.25)
            ax.set_axisbelow(True)
            # Top-left: the far-contraction bars are ~1e-4 while the axis runs to 2, so that corner
            # is empty, whereas the bottom-left is exactly where those bars sit.
            ax.text(
                0.02, 0.95, f"n={int(tot):,}", transform=ax.transAxes, fontsize=7.5, color=INK,
                va="top",
            )
        for ax in axes.flat[len(periods) :]:
            ax.axis("off")
        for ax in axes.flat[: len(periods)]:
            ax.set_yscale("log")
            ax.set_ylim(1e-5, 2.0)
        for ax in axes[:, 0]:
            ax.set_ylabel("fraction of reads\n(log scale)")
        fig.supxlabel("read's distance from the allele, in whole motif units (− = shorter)", y=-0.01)
        fig.suptitle(
            "Where reads land relative to the allele — the 0 bar is agreement",
            fontweight="bold",
        )
        fig.tight_layout()
        return fig

    # The headline of this section as its own chart rather than a percentage buried in six panel
    # titles: whether slippage moves the tract by whole motif units is a claim the read model rests
    # on, so it deserves to be readable at a glance.
    def whole_unit_figure():
        off = comp[comp["off_mode"] != 0].copy()
        off["whole_unit"] = off["off_mode"] % off["period"] == 0
        periods, fracs, counts = [], [], []
        for period in sorted(off["period"].unique()):
            if period > 6:
                continue
            g = off[off["period"] == period]
            tot = g["reads"].sum()
            if tot < 100:
                continue
            periods.append(period)
            fracs.append(g.loc[g["whole_unit"], "reads"].sum() / tot)
            counts.append(int(tot))

        fig, ax = plt.subplots(figsize=(9, 3.4))
        x = np.arange(len(periods))
        # Period 1 is grey, not blue: it is 1.0 by construction and must not read as a measurement.
        colours = [GREY if p == 1 else "#2a78d6" for p in periods]
        ax.bar(x, fracs, color=colours, width=0.7)
        for xi, (f, c, p) in enumerate(zip(fracs, counts, periods)):
            label = f"{f:.0%}\nn={c:,}" + ("\nby construction" if p == 1 else "")
            ax.text(xi, f + 0.02, label, ha="center", va="bottom", fontsize=7.5, color=INK)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{PERIOD_NAME.get(p, p)}\n({p} bp)" for p in periods], fontsize=8.5)
        ax.set_ylabel("of off-mode reads")
        ax.set_ylim(0, 1.28)
        ax.yaxis.set_major_formatter(lambda v, _p: f"{v:.0%}")
        ax.grid(True, axis="y", alpha=0.25)
        ax.set_axisbelow(True)
        ax.set_title(
            "Stutter moves the tract by whole motif units", fontweight="bold", fontsize=11
        )
        fig.tight_layout()
        return fig

    mo.vstack(
        [
            mo.md(
                "## 5 · Is stutter whole-unit, and what shape is it?\n"
                "The read model prices slippage in **whole motif units**, and treats the chance of "
                "slipping *k* units as falling off geometrically with *k*. Both are assumptions "
                "worth testing rather than asserting. A read whose length difference is *not* a "
                "multiple of the period is not slippage at all — it is an indel, an interruption, "
                "or a mis-delimited tract. **Period 1 is 1.0 by construction** — every integer is "
                "a multiple of one — so it is drawn in grey and carries no evidence either way."
            ),
            whole_unit_figure(),
            mo.md(
                "The model's premise holds where it can be tested: di and tri are 95–98% whole-unit "
                "and hexa 93%. Tetra and penta fall short, and their non-unit residue sits at ±1 **bp** "
                "— single-base indels rather than slippage, which is what the remainder *should* look "
                "like if the model is right about the rest. Both cells are thin, so treat the shortfall "
                "as a flag for the synthetic validation to settle, not as a measured rate.\n\n"
                "Given whole-unit steps, this is where reads actually land:"
            ),
            shape_figure(),
            mo.md(
                "The **0 bar is the reads that agree with the allele**, so each panel carries both "
                "halves of the question at once: how much stutter there is — the height of "
                "everything that is *not* 0, quoted in each title — and what shape it takes. The "
                "axis is logarithmic because those two live four orders of magnitude apart; on a "
                "linear axis the 0 bar would be the only thing visible.\n\n"
                "Three things to read off it. Agreement dominates: **over 99% of reads sit exactly "
                "on the allele** at every period, so stutter is a small perturbation rather than a "
                "pervasive one — but it is not small where it matters, since the grids above show "
                "it concentrating in the long-tract cells this average hides. The fall-off from ±1 "
                "is steep and roughly straight on a log axis, which is what a geometric looks like. "
                "And it is **asymmetric**, contractions outnumbering expansions by more and more as "
                "the period grows. That asymmetry is not an artefact of what we can see — a "
                "censoring explanation (long alleles being harder to span) would make it *grow* "
                "with tract length, and it does not: it is already there at the shortest tracts."
            ),
        ]
    )
    return


@app.cell
def _(comp, meta_path, mo, tsv_path):
    mo.md(
        f"""
        ---
        *Source: `{tsv_path.name}` — {comp['unit'].nunique()} units,
        {comp[['contig', 'start', 'end']].drop_duplicates().shape[0]:,} loci, ng's default STR
        delimiter, one region-typing walk (`examples/ng_ssr_cohort_stutter.rs`). Sample metadata:
        `{meta_path.name}`. Only complete reads carry a length, so partials — censored lower
        bounds — are excluded throughout.*
        """
    )
    return


if __name__ == "__main__":
    app.run()
