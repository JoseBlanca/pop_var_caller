# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas"]
# ///
"""Flank-anchored oracle — adjudicate the STR delimiter divergence (step 2).

Given the reads where the two tract delimiters disagree (the dump from
`examples/ng_ssr_divergent_reads.rs`), decide *who is right* by a method independent of
either aligner: anchor the read against the locus's **unique reference flanks** and read
off the tract that sits *between* the anchors. Anchoring on the flanks — not on the motif
— pins the tract boundary exactly where the reference defines it, which is the one thing
neither aligner does unambiguously and the whole source of the boundary-half-unit
disagreement.

For each read the oracle reports the tract length between the anchors (the truth for that
read, given the reference boundary), the unit count + interruptions, and which aligner's
measurement it matches. It is deliberately slow and cautious, and research-only — it does
not belong in ng's codebase.

  uv run benchmarks/ssr_hg002/scripts/ng_ssr_delimiter_oracle.py \
      [benchmarks/ssr_hg002/results/ng_aligner_bakeoff/divergent_reads.tsv]

Caveats it is honest about: it can only adjudicate reads where BOTH flanks anchor cleanly
(reported), and the "truth" is the read's content against the *reference-defined* boundary
— it does not know the sample's true allele, only what this molecule shows.
"""
import sys
from pathlib import Path

import pandas as pd

MIN_OVERLAP = 8  # bp of flank that must lie inside the read for a confident anchor
MAX_MISMATCH_RATE = 0.25  # tolerate SNPs / sequencing error in the anchored flank


def _place_left(read: str, flank: str):
    """Place the whole left flank so its 3' end falls at read position `t` (the tract start),
    allowing the flank's 5' end to hang off the read's left edge. Returns the `t` minimising the
    mismatch rate over the overlap (ties → more overlap), or None. Using the whole 30 bp flank —
    not just its innermost bases — is what pins the boundary: the outer, non-motif bases dominate
    the match, so the anchor cannot slide into a motif-like tract edge."""
    L = len(flank)
    best = None  # (mismatch_rate, -overlap, t)
    for t in range(MIN_OVERLAP, len(read) + 1):
        overlap = min(t, L)
        if overlap < MIN_OVERLAP:
            continue
        fseg = flank[L - overlap :]
        rseg = read[t - overlap : t]
        mism = sum(a != b for a, b in zip(fseg, rseg))
        rate = mism / overlap
        if rate <= MAX_MISMATCH_RATE:
            key = (rate, -overlap, t)
            if best is None or key < best:
                best = key
    return None if best is None else best[2]


def _place_right(read: str, flank: str, start: int):
    """Place the whole right flank so its 5' end falls at read position `e` ≥ `start` (the tract
    end), allowing the flank's 3' end to hang off the read's right edge. Returns the best `e`."""
    L = len(flank)
    best = None
    for e in range(start, len(read) - MIN_OVERLAP + 1):
        overlap = min(len(read) - e, L)
        if overlap < MIN_OVERLAP:
            continue
        fseg = flank[:overlap]
        rseg = read[e : e + overlap]
        mism = sum(a != b for a, b in zip(fseg, rseg))
        rate = mism / overlap
        if rate <= MAX_MISMATCH_RATE:
            key = (rate, -overlap, e)
            if best is None or key < best:
                best = key
    return None if best is None else best[2]


def oracle_tract(read: str, left_flank: str, right_flank: str):
    """The tract the read shows between its anchored flanks. Returns a dict with the tract
    bases and length, or {'anchored': False} if either flank could not be placed."""
    if len(left_flank) < MIN_OVERLAP or len(right_flank) < MIN_OVERLAP:
        return {"anchored": False, "reason": "flank_shorter_than_min_overlap"}
    t = _place_left(read, left_flank)
    if t is None:
        return {"anchored": False, "reason": "left_flank_not_found"}
    e = _place_right(read, right_flank, t)
    if e is None:
        return {"anchored": False, "reason": "right_flank_not_found"}
    return {"anchored": True, "bases": read[t:e]}


def purity(tract: str, motif: str):
    """Tile `motif` over `tract` from position 0; return (n_mismatch, n_full_units, remainder_bp)."""
    if not motif:
        return (0, 0, len(tract))
    mism = sum(tract[i] != motif[i % len(motif)] for i in range(len(tract)))
    return (mism, len(tract) // len(motif), len(tract) % len(motif))


def measured_len(cell: str):
    """`class:bases` -> (class, len(bases)); `none` -> (None, None)."""
    if cell == "none":
        return (None, None)
    cls, _, bases = cell.partition(":")
    return (cls, len(bases))


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).resolve().parent.parent
        / "results" / "ng_aligner_bakeoff" / "divergent_reads.tsv"
    )
    if not path.is_file():
        print(f"missing {path} — run examples/ng_ssr_divergent_reads.rs first", file=sys.stderr)
        return 2
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    for c in ("period", "ref_len"):
        df[c] = df[c].astype(int)

    rows = []
    for r in df.itertuples(index=False):
        o = oracle_tract(r.read_seq, r.left_flank, r.right_flank)
        fc, fl = measured_len(r.flat_gap)
        uc, ul = measured_len(r.unit_slip)
        rec = dict(
            contig=r.contig, start=r.start, period=r.period, motif=r.motif,
            ref_len=r.ref_len, category=r.category,
            flat_class=fc, flat_len=fl, unit_class=uc, unit_len=ul,
            anchored=o["anchored"], anchor_reason=o.get("reason", "ok"),
        )
        if o["anchored"]:
            tract = o["bases"]
            mism, units, rem = purity(tract, r.motif)
            rec.update(oracle_len=len(tract), oracle_units=units, oracle_rem=rem,
                       oracle_mism=mism, oracle_tract=tract)
        rows.append(rec)
    res = pd.DataFrame(rows)

    # ---- coverage of the oracle, and its intrinsic bias ----
    n = len(res)
    anc = int(res["anchored"].sum())
    print(f"divergent reads: {n}   oracle anchored both flanks: {anc} ({100*anc/n:.0f}%)")
    print("why the rest can't be adjudicated (a read that doesn't span both flanks has no")
    print("independent full-length measurement — the long-tract / partial regime):")
    print(res.loc[~res.anchored, "anchor_reason"].value_counts().to_string())
    bc = res[(res.category == "length_differ") & (res.flat_class == "complete") & (res.unit_class == "complete")]
    print(f"\nanchoring rate within both-complete length_differ: "
          f"{int(bc.anchored.sum())}/{len(bc)} ({100*bc.anchored.mean():.0f}%)")
    band = pd.cut(bc.ref_len, [0, 15, 25, 40, 10**9], labels=["<15", "15-24", "25-39", "40+"])
    print("anchoring rate by reference tract length (the short-tract bias):")
    print(bc.groupby(band, observed=True).anchored.agg(["mean", "size"]).assign(
        mean=lambda d: (100 * d["mean"]).round().astype(int)).rename(columns={"mean": "anchored_%", "size": "n"}).to_string())
    print()

    # ---- core adjudication: both-complete length_differ, oracle-anchored ----
    core = res[
        (res.category == "length_differ")
        & (res.flat_class == "complete")
        & (res.unit_class == "complete")
        & (res.anchored)
    ].copy()
    core["flat_err"] = (core.flat_len - core.oracle_len).abs()
    core["unit_err"] = (core.unit_len - core.oracle_len).abs()

    def winner(x):
        if x.flat_err < x.unit_err:
            return "flat_closer"
        if x.unit_err < x.flat_err:
            return "unit_closer"
        return "tie"

    core["winner"] = core.apply(winner, axis=1)
    print(f"=== ADJUDICATION: both-complete length_differ, oracle-anchored ({len(core)} reads) ===")
    print("truth = tract between the reference-flank anchors in the read")
    print("who matches the oracle:", core.winner.value_counts().to_dict())
    print(f"exact-match rate  flat_gap={100*(core.flat_len==core.oracle_len).mean():.0f}%  "
          f"unit_slip={100*(core.unit_len==core.oracle_len).mean():.0f}%")
    print(f"mean abs error    flat_gap={core.flat_err.mean():.2f} bp   unit_slip={core.unit_err.mean():.2f} bp")
    print(f"flat_gap == reference length (collapse): {int((core.flat_len==core.ref_len).sum())}/{len(core)}"
          f"  (of those, oracle ≠ ref: {int(((core.flat_len==core.ref_len)&(core.oracle_len!=core.ref_len)).sum())})")
    print()
    print("winner by period:")
    print(core.groupby('period').winner.value_counts().unstack().fillna(0).astype(int).to_string())
    print()

    # ---- separate genuine error from the boundary half-unit convention ----
    # A read is a "boundary-only" disagreement when the oracle sits exactly at one aligner's
    # length and the other differs by < one period (a trailing partial unit), i.e. nobody is
    # grossly wrong — it is a definitional edge.
    def kind(x):
        d = abs(x.flat_len - x.unit_len)
        near = d < x.period or (x.period == 1 and d <= 1)
        matches_one = (x.flat_len == x.oracle_len) or (x.unit_len == x.oracle_len)
        if near and matches_one:
            return "boundary_convention (<1 unit, one matches oracle)"
        if matches_one:
            return "one_correct_other_wrong"
        return "both_off_from_oracle"
    core["kind"] = core.apply(kind, axis=1)
    print("nature of the disagreement:")
    print(core.kind.value_counts().to_string())
    print()
    real = core[core.kind == "one_correct_other_wrong"]
    if len(real):
        who = pd.Series(["unit_slip" if r.unit_len == r.oracle_len else "flat_gap"
                         for r in real.itertuples()]).value_counts().to_dict()
        print(f"among the {len(real)} genuine one-right-one-wrong reads, the CORRECT aligner was: {who}")
    print()

    # ---- interruptions: are divergent tracts impure? ----
    ancd = res[res.anchored]
    print(f"tract purity over all {len(ancd)} anchored divergent reads: "
          f"{int((ancd.oracle_mism==0).sum())} perfect, "
          f"{int((ancd.oracle_mism>0).sum())} with ≥1 interruption "
          f"(mean {ancd.oracle_mism.mean():.2f} mismatches/tract)")
    print()

    # ---- a few worked examples ----
    print("=== worked examples (oracle vs the two aligners) ===")
    show = core[core.kind == "one_correct_other_wrong"].head(5)
    for r in show.itertuples():
        tag_f = "✓" if r.flat_len == r.oracle_len else "✗"
        tag_u = "✓" if r.unit_len == r.oracle_len else "✗"
        print(f"{r.contig}:{r.start} p{r.period}({r.motif}) ref={r.ref_len}bp | "
              f"ORACLE={r.oracle_len}bp ({r.oracle_units}u+{r.oracle_rem}) | "
              f"flat={r.flat_len}{tag_f} unit={r.unit_len}{tag_u} | tract={r.oracle_tract}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
