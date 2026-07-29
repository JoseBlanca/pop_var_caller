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
            "read_group", "sample", "library", "library_origin",
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
    PERIOD_NAME = {1: "mono", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
    return LEN_LABELS, PERIOD_NAME, comp


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
    # SECTION 2 — the onset. Off-mode fraction against tract length, one line per sample, faceted
    # by period. Lines are coloured by duplicate rate so the PCR question is answered by looking:
    # if amplification drives stutter, the dark lines separate upward.
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
                "## 2 · From what length does stuttering start to matter\n"
                "One line per sample, faceted by motif period, coloured by duplicate rate. Read "
                "the **onset** as the length band where a period's lines lift off the floor — it "
                f"differs by period, which is why this is not one curve. Cells with fewer than "
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
