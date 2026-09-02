#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# ///
"""Does ng's repeat-tract candidate selection offer the sequences HG002 really carries?

`spec/candidate_alleles_ssr.md` §4.1 and §5 measured that offline in August 2026, through the
existing caller's Stage-1 pileup and a scorer that no longer exists. This scores the same thing
through **the shipped `select_ssr`**: it reads the per-tract candidate dump
`examples/ng_candidate_selection_probe` writes with `NG_TRACT_DUMP`, reconstructs both of HG002's
haplotype sequences at each tract from the reference plus the phased truth VCF's edits inside it,
and reports how often every true sequence is among the candidates.

Three things it reports that the offline tables did not, because each of them decides how a
candidates-per-tract figure should be read:

  * **the denominator.** Selection only ever sees a tract the cohort merge built a locus for —
    a set already enriched for tracts that vary. A refused tract is left at its reference
    sequence and carries exactly one candidate, so it is counted, and the dump carries both.
  * **the truth floor** — the reference plus this sample's two true tract sequences,
    deduplicated. That is what a perfect selector would offer, so the excess over it is the
    only part of the figure a rule can be blamed for.
  * **where a missing true sequence was lost** — four places, and only the last two are
    selection's: the merge refused the tract, no read carried the sequence so the merge's table
    never held it, the support bar refused it, or it cleared the bar and the per-sample
    top-`ploidy` rung cut dropped it anyway.

Pass a production `.cat` catalog instead of a dump to get the truth floor and the class sizes of
*that* tract set, which is how two runs over different catalogs are made comparable.

Usage (from ssr_hg002/, or with SSR_HG002_ROOT pointing at it):
    ng_tract_candidate_recall.py <dump.tsv | catalog.cat> [...]
"""

import gzip
import os
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("SSR_HG002_ROOT", Path(__file__).resolve().parent.parent))
REFERENCE = (
    ROOT.parent
    / "giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
)
TRUTH = ROOT / "truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz"

CLASSES = (
    "homozygous",
    "heterozygous, different repeat length",
    "heterozygous, same repeat length, different spelling",
)


class IndexedFasta:
    """Random access to a FASTA through its `.fai`, 1-based and inclusive at both ends."""

    def __init__(self, path):
        self.handle = open(path, "rb")
        self.index = {}
        self.names = []
        with open(str(path) + ".fai") as fh:
            for line in fh:
                name, length, offset, line_bases, line_width = line.split("\t")[:5]
                self.index[name] = (
                    int(length),
                    int(offset),
                    int(line_bases),
                    int(line_width),
                )
                self.names.append(name)

    def fetch(self, chrom, start, end):
        length, offset, line_bases, line_width = self.index[chrom]
        start, end = max(1, start), min(length, end)
        if end < start:
            return ""
        begin = offset + (start - 1) // line_bases * line_width + (start - 1) % line_bases
        stop = offset + (end - 1) // line_bases * line_width + (end - 1) % line_bases
        self.handle.seek(begin)
        return self.handle.read(stop - begin + 1).replace(b"\n", b"").decode().upper()


def read_truth(path):
    """The phased truth records per contig, ascending by position."""
    rows = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            genotype = f[9].split(":")[0].replace("/", "|")
            if "|" not in genotype:
                continue
            left, right = genotype.split("|")[:2]
            if left == "." or right == ".":
                continue
            rows[f[0]].append(
                (int(f[1]), f[3].upper(), f[4].upper().split(","), int(left), int(right))
            )
    for value in rows.values():
        value.sort()
    return rows, {chrom: [r[0] for r in v] for chrom, v in rows.items()}


def records_over(truth, starts, chrom, start, end):
    """Truth records whose REF span touches the 1-based inclusive window."""
    rows = truth.get(chrom)
    if not rows:
        return []
    at = bisect_right(starts[chrom], end)
    return [
        rows[i]
        for i in range(max(0, at - 64), at)
        if rows[i][0] + len(rows[i][1]) - 1 >= start
    ]


def haplotype(fasta, chrom, start, end, which, rows):
    """One haplotype's tract sequence, or None where a truth record crosses the tract's end.

    The window is `[start - 1, end]` so that a left-anchored indel sitting just before the tract
    — which is where left alignment puts every repeat-length change — is applied; the anchor
    base is dropped afterwards.
    """
    window_start = start - 1
    window = fasta.fetch(chrom, window_start, end)
    pieces, at = [], window_start
    for pos, ref, alts, left, right in rows:
        allele = left if which == 0 else right
        if pos < at:
            continue  # two records overlapping on one haplotype: keep the first
        if pos < window_start or pos + len(ref) - 1 > end:
            return None
        pieces.append(window[at - window_start : pos - window_start])
        pieces.append(ref if allele == 0 else alts[allele - 1])
        at = pos + len(ref)
    pieces.append(window[at - window_start :])
    return "".join(pieces)[1:]


def read_tracts(path, fasta):
    """One record per tract: its candidates, and every allele the merge's table held."""
    tracts = {}
    if str(path).endswith(".cat"):
        contig_of = {name: at for at, name in enumerate(fasta.names)}
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] in contig_of:
                    tracts[(contig_of[f[0]], int(f[1]) + 1, int(f[2]))] = {
                        "built": True,
                        "period": len(f[3]),
                        "candidates": set(),
                        "table": {},
                    }
        return tracts
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (int(f[0]), int(f[1]), int(f[2]))
            tract = tracts.setdefault(
                key,
                {
                    "built": f[5] == "1",
                    "period": int(f[4]),
                    "candidates": set(),
                    "table": {},
                },
            )
            if f[13]:
                tract["table"][f[13]] = (f[10] == "1", float(f[12]))
                if f[9] == "1":
                    tract["candidates"].add(f[13])
    return tracts


def score(path, fasta, truth, starts):
    tracts = read_tracts(path, fasta)
    per_class = defaultdict(lambda: [0, 0])
    # Spec §4.1 scores *repeat counts*, not spellings: a rule that offers the right two
    # lengths has done its nomination job even if the second spelling at one of them is
    # missing. The two measures answer different questions and both are reported.
    lengths_offered = [0, 0]
    adjacent = [0, 0]
    lost = defaultdict(int)
    refused_shares = []
    scored = skipped = candidates = floor = beyond_truth = 0
    for (contig, start, end), tract in tracts.items():
        chrom = fasta.names[contig]
        period = tract["period"]
        reference = fasta.fetch(chrom, start, end)
        offered = tract["candidates"] or {reference}
        rows = records_over(truth, starts, chrom, start, end)
        first = haplotype(fasta, chrom, start, end, 0, rows)
        second = haplotype(fasta, chrom, start, end, 1, rows)
        if first is None or second is None:
            skipped += 1
            continue
        scored += 1
        candidates += len(offered)
        true_set = {first, second, reference}
        floor += len(true_set)
        beyond_truth += len(offered - true_set)
        if first == second:
            name = CLASSES[0]
        elif len(first) != len(second):
            name = CLASSES[1]
        else:
            name = CLASSES[2]
        per_class[name][0] += 1
        if name == CLASSES[1]:
            offered_lengths = {len(seq) for seq in offered}
            lengths_offered[0] += 1
            both_lengths = len(first) in offered_lengths and len(second) in offered_lengths
            lengths_offered[1] += int(both_lengths)
            if period and abs(len(first) - len(second)) == period:
                adjacent[0] += 1
                adjacent[1] += int(both_lengths)
        if first in offered and second in offered:
            per_class[name][1] += 1
            continue
        for missing in {first, second} - offered:
            if not tract["built"]:
                lost["the merge refused the tract, so no locus was built"] += 1
            elif missing not in tract["table"]:
                lost["no read carried it, so the merge's table never held it"] += 1
            elif tract["table"][missing][0]:
                lost["it cleared the support bar; the per-sample top-ploidy cut dropped it"] += 1
                refused_shares.append(tract["table"][missing][1])
            else:
                lost["the merge's table held it, and the support bar refused it"] += 1
                refused_shares.append(tract["table"][missing][1])

    print(f"## {path}")
    print(f"# tracts: {len(tracts)}   scored: {scored}   "
          f"skipped (a truth record crosses the tract's end): {skipped}")
    denominator = max(1, scored)
    print(f"candidates_per_tract\t{candidates / denominator:.3f}")
    print(f"truth_floor_per_tract\t{floor / denominator:.3f}"
          "\t(the reference and this sample's two true sequences, deduplicated)")
    print(f"beyond_truth_per_tract\t{beyond_truth / denominator:.3f}"
          "\t(neither a true sequence nor the reference)")
    print(f"\n{'class':<55}{'tracts':>8}{'both offered':>14}")
    for name in CLASSES:
        total, hit = per_class[name]
        print(f"{name:<55}{total:>8}{f'{100.0 * hit / total:.1f}%' if total else '-':>14}")
    print(f"\n{'both true repeat LENGTHS offered (spec §4.1s measure)':<55}"
          f"{'tracts':>8}{'both offered':>14}")
    for label, (total, hit) in (
        ("heterozygous, different repeat length", lengths_offered),
        ("  of those, one repeat apart", adjacent),
    ):
        print(f"{label:<55}{total:>8}{f'{100.0 * hit / total:.1f}%' if total else '-':>14}")
    if lost:
        print("\n# where a true tract sequence was lost, over every miss")
        for reason, count in sorted(lost.items(), key=lambda kv: -kv[1]):
            print(f"{count:>6}  {reason}")
        if refused_shares:
            refused_shares.sort()
            middle = refused_shares[len(refused_shares) // 2]
            print(f"# the ones selection dropped carried a median {100 * middle:.1f} reads "
                  "in 100 of some sample's own")
    print()


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    fasta = IndexedFasta(REFERENCE)
    truth, starts = read_truth(TRUTH)
    for path in sys.argv[1:]:
        score(path, fasta, truth, starts)


if __name__ == "__main__":
    main()
