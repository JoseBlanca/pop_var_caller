#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["pysam"]
# ///
"""At one indel, ask every overlapping read which of the two haplotypes its bases
fit — independently of where the aligner put it.

    read_haplotype_support.py --bam READS.bam --reference REF.fa \
        --chrom chr15 --pos 96140584 --ref-allele C --alt-allele CTGATT... [--per-read]

    read_haplotype_support.py --bam ... --reference ... --sites SITES.tsv

WHAT IT ANSWERS. A caller that says "twelve reads support the reference here" is
making a claim about the reads, and it can be checked without trusting the caller or
the aligner. This builds the two sequences the sample could have — the reference
window, and the same window with the indel applied — and slides each read along both,
counting mismatches. The indel is then *inside* the sequence being matched rather
than in a gap whose cost has to be priced, so no gap-scoring choice enters the answer.

    supports the reference     fewer mismatches against the reference haplotype
    supports the alternative   fewer against the alternative one
    uninformative              the same against both — the read cannot tell

The third bucket is the one that matters for reading a caller's `AD`. A read the
aligner placed as reference, at a repeat where the alternative haplotype is nearly
as good a match, lands there rather than in the first — so a caller crediting it to
the reference allele is stating more than its reads support.

`--sites` takes a TSV of `sample<TAB>chrom<TAB>pos<TAB>ref<TAB>alt<TAB>genotype`, the
shape `bcftools norm -m -any` output gives after splitting, and `--bam` is then a
template with `{sample}` in it.
"""

import argparse
import sys
from pathlib import Path

import pysam


def fewest_mismatches(read: str, haplotype: str) -> int:
    """How far the read is from the haplotype at the offset that suits it best.

    Ungapped: the two haplotypes already differ by the indel, so a read belonging to
    either one lines up against it with substitutions alone.
    """
    best = len(read)
    for offset in range(len(haplotype) - len(read) + 1):
        mismatches = sum(
            1 for i, base in enumerate(read) if haplotype[offset + i] != base
        )
        if mismatches < best:
            best = mismatches
            if best == 0:
                break
    return best


def verdicts_at(bam, fasta, chrom, pos, ref_allele, alt_allele, flank, per_read):
    """One site's tally. `pos` is a 1-based VCF POS."""
    start = pos - 1 - flank
    stop = pos - 1 + len(ref_allele) + flank
    window = fasta.fetch(chrom, start, stop).upper()
    anchor = pos - 1 - start
    found = window[anchor:anchor + len(ref_allele)]
    if found != ref_allele.upper():
        sys.exit(
            f"{chrom}:{pos} the reference holds {found}, not the {ref_allele} given"
        )
    reference_haplotype = window
    alternative_haplotype = (
        window[:anchor] + alt_allele.upper() + window[anchor + len(ref_allele):]
    )

    counts = {"reference": 0, "alternative": 0, "uninformative": 0}
    for read in bam.fetch(chrom, pos - 1, pos - 1 + len(ref_allele)):
        if (read.is_unmapped or read.is_duplicate
                or read.is_secondary or read.is_supplementary):
            continue
        sequence = (read.query_sequence or "").upper()
        if not sequence or len(sequence) > len(window) - 2:
            continue
        against_reference = fewest_mismatches(sequence, reference_haplotype)
        against_alternative = fewest_mismatches(sequence, alternative_haplotype)
        verdict = (
            "reference" if against_reference < against_alternative
            else "alternative" if against_alternative < against_reference
            else "uninformative"
        )
        counts[verdict] += 1
        if per_read:
            print(f"  {read.reference_start + 1}\t{read.cigarstring}\t"
                  f"vs REF {against_reference}\tvs ALT {against_alternative}\t{verdict}")
    return counts


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True,
                   help="alignments; with --sites, a template containing {sample}")
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument("--sites", type=Path,
                   help="TSV: sample chrom pos ref alt genotype")
    p.add_argument("--chrom")
    p.add_argument("--pos", type=int, help="1-based VCF POS")
    p.add_argument("--ref-allele")
    p.add_argument("--alt-allele")
    p.add_argument("--flank", type=int, default=200,
                   help="reference bases each side of the indel to match against")
    p.add_argument("--per-read", action="store_true")
    p.add_argument("--label", default="")
    a = p.parse_args()

    fasta = pysam.FastaFile(str(a.reference))
    if a.sites:
        bams = {}
        for line in a.sites.read_text().splitlines():
            if not line.strip():
                continue
            sample, chrom, pos, ref, alt, genotype = line.split("\t")
            if sample not in bams:
                bams[sample] = pysam.AlignmentFile(a.bam.format(sample=sample))
            counts = verdicts_at(bams[sample], fasta, chrom, int(pos), ref, alt,
                                 a.flank, a.per_read)
            print(f"{sample}\t{chrom}:{pos}\ttruth={genotype}\t"
                  f"ref{len(ref)}->alt{len(alt)}\treads={sum(counts.values())}\t"
                  f"reference={counts['reference']}\t"
                  f"alternative={counts['alternative']}\t"
                  f"uninformative={counts['uninformative']}")
        return

    missing = [n for n in ("chrom", "pos", "ref_allele", "alt_allele")
               if getattr(a, n) is None]
    if missing:
        p.error("without --sites, one site must be given: "
                + ", ".join("--" + n.replace("_", "-") for n in missing))
    bam = pysam.AlignmentFile(a.bam)
    counts = verdicts_at(bam, fasta, a.chrom, a.pos, a.ref_allele, a.alt_allele,
                         a.flank, a.per_read)
    print(f"{a.label or f'{a.chrom}:{a.pos}'}\t"
          f"ref{len(a.ref_allele)}->alt{len(a.alt_allele)}\t"
          f"reads={sum(counts.values())}\treference={counts['reference']}\t"
          f"alternative={counts['alternative']}\t"
          f"uninformative={counts['uninformative']}")


main()
