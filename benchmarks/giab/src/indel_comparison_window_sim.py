#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["pysam"]
# ///
"""How many reads ng would credit to the reference allele at an indel if the window
it compares them over were wider than the one it uses today.

    indel_comparison_window_sim.py --reference REF.fa --sites SITES.tsv \
        --bam-template '.../bam/30x/{sample}.30x.seed42.bam' \
        [--class insertion|deletion] [--width allele|N] [--truth-gt 1/1]

WHY IT EXISTS. ng records a read as an exact reference observation when the bases it
showed across a locus's *reference* positions equal the reference there. An
insertion's reference footprint is one base — its anchor — so a read agreeing with
the reference at that single base counts as reference evidence whatever it shows
downstream, and long homozygous insertions collect enough such reads to be genotyped
heterozygous (`doc/devel/reports/reviews/ng_indel_genotypes_vs_giab_2026-09-11.md`).

This rebuilds the walker's own classification over a chosen window so the two can be
compared before any of it is built. For each read it collects what the read showed
across `anchor .. anchor + width`, the way the walker does — inserted bases anchored
inside the window included, deleted positions contributing positions but no bases —
and sorts it:

    reference   covered every position of the window, and its bases are the
                reference's
    other       covered every position and showed something else
    partial     did not reach the far side, so ng's own rule leaves it uninformative

`--width allele` makes the window the indel's own length (the change that would turn
`ReadEvent::footprint_span` for an insertion from 1 into `inserted_len + 1`);
`--width 1` reproduces today's rule and is the control to check the model against the
caller's real `AD`; `--width N` is a fixed confirmation window.

`--sites` is a TSV of `sample<TAB>chrom<TAB>pos<TAB>ref<TAB>alt<TAB>genotype`, the
shape `bcftools norm -m -any` gives after splitting to one ALT a record. It applies
none of ng's read filters (mapping quality, base masking, the depth cap), so it
models the rule rather than reproducing the caller.
"""

import argparse
import collections
from pathlib import Path

import pysam

CIGAR_MATCH = (0, 7, 8)     # M = X
CIGAR_INSERTION = 1         # I
CIGAR_SKIPS_REFERENCE = (2, 3)   # D N
CIGAR_SOFT_CLIP = 4         # S


def observation(read, anchor, width):
    """What this read showed across reference positions `[anchor, anchor + width)`,
    and whether it covered all of them — the walker's two outputs for one read."""
    reference_position = read.reference_start
    query_position = 0
    shown = []
    covered = set()
    stop = anchor + width
    for op, length in read.cigartuples or []:
        if op in CIGAR_MATCH:
            for k in range(length):
                if anchor <= reference_position + k < stop:
                    shown.append(read.query_sequence[query_position + k].upper())
                    covered.add(reference_position + k)
            reference_position += length
            query_position += length
        elif op == CIGAR_INSERTION:
            # An insertion is anchored at the reference base before it.
            if anchor <= reference_position - 1 < stop:
                shown.extend(
                    b.upper()
                    for b in read.query_sequence[query_position:query_position + length]
                )
            query_position += length
        elif op in CIGAR_SKIPS_REFERENCE:
            for k in range(length):
                if anchor <= reference_position + k < stop:
                    covered.add(reference_position + k)
            reference_position += length
        elif op == CIGAR_SOFT_CLIP:
            query_position += length
    return "".join(shown), len(covered) == width


def length_bucket(indel_length: int) -> str:
    n = abs(indel_length)
    return ("1" if n == 1 else "2-3" if n <= 3 else "4-6" if n <= 6
            else "7-12" if n <= 12 else "13+")


BUCKETS = ["1", "2-3", "4-6", "7-12", "13+"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument("--sites", required=True, type=Path,
                   help="TSV: sample chrom pos ref alt genotype")
    p.add_argument("--bam-template", required=True,
                   help="alignment path with {sample} in it")
    p.add_argument("--class", dest="indel_class", default="insertion",
                   choices=["insertion", "deletion"])
    p.add_argument("--width", default="allele",
                   help="'allele' (the indel's own length) or a fixed base count")
    p.add_argument("--truth-gt", default="1/1",
                   help="which truth genotype to restrict to")
    p.add_argument("--per-site", action="store_true")
    a = p.parse_args()

    fasta = pysam.FastaFile(str(a.reference))
    bams = {}
    totals = collections.defaultdict(lambda: [0, 0, 0, 0, 0])
    per_site = []

    for line in a.sites.read_text().splitlines():
        if not line.strip():
            continue
        sample, chrom, pos, ref, alt, genotype = line.split("\t")
        if "/".join(sorted(genotype.replace("|", "/").split("/"))) != a.truth_gt:
            continue
        indel_length = len(alt) - len(ref)
        if indel_length == 0 or (indel_length > 0) != (a.indel_class == "insertion"):
            continue
        width = (abs(indel_length) + len(ref)) if a.width == "allele" else int(a.width)
        anchor = int(pos) - 1
        reference_bases = fasta.fetch(chrom, anchor, anchor + width).upper()
        if sample not in bams:
            bams[sample] = pysam.AlignmentFile(a.bam_template.format(sample=sample))

        counts = collections.Counter()
        for read in bams[sample].fetch(chrom, anchor, anchor + 1):
            if (read.is_unmapped or read.is_duplicate
                    or read.is_secondary or read.is_supplementary):
                continue
            shown, complete = observation(read, anchor, width)
            counts["reference" if (complete and shown == reference_bases)
                   else "other" if complete else "partial"] += 1

        row = totals[length_bucket(indel_length)]
        row[0] += 1
        row[1] += sum(counts.values())
        row[2] += counts["reference"]
        row[3] += counts["other"]
        row[4] += counts["partial"]
        if a.per_site:
            per_site.append((sample, f"{chrom}:{pos}", indel_length, dict(counts)))

    print(f"{'indel len':>10} {'sites':>6} {'reads':>7} {'reference':>10} "
          f"{'other':>7} {'partial':>8} {'ref share':>10}")
    for bucket in BUCKETS:
        if bucket not in totals:
            continue
        sites, reads, reference, other, partial = totals[bucket]
        share = reference / reads if reads else 0.0
        print(f"{bucket:>10} {sites:>6} {reads:>7} {reference:>10} {other:>7} "
              f"{partial:>8} {share:>10.3f}")
    for sample, site, indel_length, counts in per_site:
        print(f"    {sample} {site} {indel_length:+d} {counts}")


main()
