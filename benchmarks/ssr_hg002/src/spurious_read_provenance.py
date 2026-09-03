#!/usr/bin/env python3
"""Where a spurious second allele's reads come from — P1 of the tract-accuracy program.

At each tract the caller marks heterozygous where the truth is homozygous (the
`spurious_het` verdicts of `tract_qual_experiment.py --verdicts-out`), this pulls
the reads that carry the spurious length and asks two questions the program's
levers divide on (`doc/devel/ng/research/tract_accuracy_program.md` §3, P1):

* **Are the reads independent evidence?** Reads sharing one strand beyond what
  the tract's own strand mix allows, or sharing one exact (start, strand,
  template length) signature — the shape of a PCR duplicate family, which these
  BAMs do not flag — are one observation counted several times. That is lever
  L2's premise.
* **Does the locus keep producing the length at 300×?** The 30× alignment is a
  subsample of the 300× one, so the 30× reads alone reproduce only about a tenth
  of their 30× share there. A share that holds at 300× is the locus's own
  behaviour — its slippage, or the realigner's spelling of its reads (L3/L4).
  A share that collapses is a draw more reads dilute, which no lever fixes.

Each tract lands in exactly one bucket, in this order:

  spelling_only   the spurious sequence has the truth's length — reads binned by
                  length cannot arbitrate it; it belongs to the realigner's story
  unseen_in_raw   NO raw-aligned 30x read spells the spurious length at all —
                  the caller offered and called a length the aligner's own
                  spellings do not carry, so the read questions are vacuous and
                  the tract belongs to the realigner's story (L4); the 300x
                  count is still reported, because a length absent at 300x too
                  is the strongest form of it
  both            clustered AND persistent — printed one by one, expected rare
  clustered       a strand or duplicate-family signal at 30x, and no persistence
  locus_real      persists at 300x (at least 3 reads and at least half the 30x
                  share), no clustering signal
  sampling_noise  neither — independent reads whose share collapses at depth

The criteria live in `doc/devel/ng/research/tract_accuracy_program_report.md`
§P1 and are not knobs; nothing here is settable, so two runs cannot disagree
about what was measured. `unseen_in_raw` is the one post-registration
amendment: the pre-registered rule divided a share of zero and called the
result persistence, and the report's P1 section records the deviation.

**The probe's own control**: the set of tracts it derives must equal the
`spurious_het` rows of the verdict dump it is handed — asserted, not eyeballed.

Usage:

    spurious_read_provenance.py --reference ref.fa --truth truth.vcf.gz \\
        --query calls.vcf --confident-bed tier_sorted.bed --tract-bed tier.bed \\
        --verdicts verdicts.tsv --verdicts-arm baseline --verdicts-depth 30x \\
        --bam-low 30x.bam --bam-high 300x.bam \\
        --per-tract-out tracts.tsv --cases-out cases.txt

Needs `samtools` and `bcftools`; everything else is the standard library.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from tract_qual_experiment import (  # noqa: E402
    TractGround,
    haplotype_pairs,
    offered_sequences,
    period_class_of,
    prepared_vcf,
    read_vcf,
    records_near_each_tract,
    sample_columns_of,
    tract_reference_bases,
)

# The pre-registered persistence bar: at least this many 300x reads at the
# spurious length, and at least half the 30x share.
PERSIST_MIN_READS = 3
PERSIST_SHARE_FACTOR = 0.5
# The pre-registered strand bar: all k on one strand, and the chance of that
# under the tract's own strand mix below 1 in 20 (needs k >= 3 in practice).
STRAND_P = 0.05


@dataclass
class SpanningRead:
    """One primary mapped read spanning the tract, reduced to what P1 asks."""

    start: int
    reverse: bool
    template: int
    tract_length: int


def spanning_reads(
    bam: Path, contig: str, tract_first: int, tract_last: int
) -> list[SpanningRead]:
    """Every primary mapped read covering the tract plus one base each side.

    A read's tract length is the aligned bases it places on the tract's own
    positions, plus inserted bases anchored from the base before the tract
    through the tract's second-to-last position — a left-aligned insertion at a
    repeat's start sits on the anchor base before it. The rule is one-sided at
    the right edge on purpose, and it is the same rule at both depths, so the
    persistence comparison cancels any edge convention.
    """
    completed = subprocess.run(
        ["samtools", "view", str(bam), f"{contig}:{tract_first}-{tract_last}"],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise SystemExit(f"samtools view failed:\n{completed.stderr.strip()}")
    out: list[SpanningRead] = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        flag = int(fields[1])
        # Unmapped, secondary, supplementary, or QC-fail reads are not evidence
        # the caller admits.
        if flag & (0x4 | 0x100 | 0x200 | 0x800):
            continue
        pos = int(fields[3])
        cigar = fields[5]
        if cigar == "*":
            continue
        ref_at = pos
        in_tract = 0
        length = 0
        number = ""
        for character in cigar:
            if character.isdigit():
                number += character
                continue
            count = int(number)
            number = ""
            if character in "M=X":
                low = max(ref_at, tract_first)
                high = min(ref_at + count - 1, tract_last)
                if high >= low:
                    in_tract += high - low + 1
                ref_at += count
            elif character in "DN":
                ref_at += count
            elif character == "I":
                # Anchored between ref_at - 1 and ref_at.
                if tract_first - 1 <= ref_at - 1 <= tract_last - 1:
                    in_tract += count
            # S and H consume no reference and P nothing at all.
        aligned_last = ref_at - 1
        if pos <= tract_first - 1 and aligned_last >= tract_last + 1:
            out.append(
                SpanningRead(pos, bool(flag & 0x10), int(fields[8]), in_tract)
            )
    return out


def one_strand_probability(k_reverse: int, k: int, mix_reverse: float) -> float:
    """The chance k independent draws land all-forward or all-reverse."""
    if 0 < k_reverse < k:
        return 1.0
    return mix_reverse**k + (1.0 - mix_reverse) ** k


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--confident-bed", type=Path, required=True)
    parser.add_argument("--tract-bed", type=Path, required=True)
    parser.add_argument(
        "--verdicts", type=Path, required=True,
        help="the verdict dump whose spurious_het rows this probe must reproduce",
    )
    parser.add_argument("--verdicts-arm", required=True)
    parser.add_argument("--verdicts-depth", required=True)
    parser.add_argument("--bam-low", type=Path, required=True,
                        help="the alignment the callset was made from")
    parser.add_argument("--bam-high", type=Path, required=True,
                        help="the deep alignment for the persistence check")
    parser.add_argument("--genotype-sample", default="HG002")
    parser.add_argument("--per-tract-out", type=Path, required=True)
    parser.add_argument("--cases-out", type=Path, required=True,
                        help="a few whole cases per bucket, reads and all")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # The spurious-het tracts, derived exactly as the instrument derives them.
    # ------------------------------------------------------------------
    expected_tracts = set()
    with open(args.verdicts, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        column = {name: header.index(name) for name in
                  ("arm", "depth", "contig", "start", "end", "verdict")}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if (fields[column["arm"]] == args.verdicts_arm
                    and fields[column["depth"]] == args.verdicts_depth
                    and fields[column["verdict"]] == "spurious_het"):
                expected_tracts.add((fields[column["contig"]],
                                     int(fields[column["start"]]),
                                     int(fields[column["end"]])))
    if not expected_tracts:
        raise SystemExit(f"{args.verdicts}: no spurious_het rows for arm "
                         f"{args.verdicts_arm} at {args.verdicts_depth}")

    with tempfile.TemporaryDirectory() as scratch:
        work = Path(scratch)
        query_path = prepared_vcf(
            args.query, reference=args.reference, masks=[args.confident_bed],
            split_multiallelics=False, pass_only=False, work=work, name="query",
        )
        truth_path = prepared_vcf(
            args.truth, reference=args.reference, masks=[args.confident_bed],
            split_multiallelics=False, pass_only=True, work=work, name="truth",
        )
        query_records, _ = read_vcf(query_path, require_qual=True)
        truth_records, _ = read_vcf(truth_path, require_qual=False)
        query_column = sample_columns_of(query_path).index(args.genotype_sample)
        truth_column = sample_columns_of(truth_path).index(args.genotype_sample)
        ground = TractGround(args.tract_bed)
        bases_of = tract_reference_bases(args.reference, ground, work)
    truth_near = records_near_each_tract(truth_records, ground)
    query_near = records_near_each_tract(query_records, ground)

    found = {}
    for key, truths in truth_near.items():
        contig, start, end, period = key
        queries = query_near.get(key)
        if queries is None:
            continue
        truth_genotypes = [one.genotype_indices(truth_column) for one in truths]
        query_genotypes = [one.genotype_indices(query_column) for one in queries]
        if any(one is None for one in truth_genotypes):
            continue
        if any(one is None for one in query_genotypes):
            continue
        window = bases_of.get((contig, start, end))
        if window is None:
            continue
        first, reference = window
        tract = (start + 1, end)
        expected = haplotype_pairs(
            truths, truth_genotypes, first, reference,
            all("|" in one.samples[truth_column] for one in truths), tract,
        )
        called = haplotype_pairs(
            queries, query_genotypes, first, reference, False, tract
        )
        if expected is None or called is None or expected & called:
            continue
        offered = offered_sequences(queries, first, reference, tract)
        if any(one not in offered for pair in expected for one in pair):
            continue
        if not all(len(set(pair)) == 1 for pair in expected):
            continue
        if all(len(set(pair)) == 1 for pair in called):
            continue
        truth_sequence = sorted(expected)[0][0]
        pair = sorted(called)[0]
        wrong = [one for one in pair if one != truth_sequence]
        found[(contig, start, end)] = (period, truth_sequence, pair, wrong)

    if set(found) != expected_tracts:
        missing = expected_tracts - set(found)
        extra = set(found) - expected_tracts
        raise SystemExit(
            f"the probe's tracts do not match the verdict dump: "
            f"{len(missing)} missing, {len(extra)} extra — first of each: "
            f"{sorted(missing)[:1]} {sorted(extra)[:1]}"
        )
    print(f"control: {len(found)} tracts, exactly the dump's spurious_het rows",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # The reads.
    # ------------------------------------------------------------------
    rows = []
    buckets: Counter = Counter()
    share_bins: Counter = Counter()
    cases: dict[str, list[str]] = {}
    for (contig, start, end), (period, truth_sequence, pair, wrong) in sorted(
        found.items()
    ):
        period_class = period_class_of(period)
        tract_first, tract_last = start + 1, end
        spurious_lengths = sorted(
            {len(one) for one in wrong if len(one) != len(truth_sequence)}
        )
        if not spurious_lengths:
            bucket = "spelling_only"
            rows.append((contig, start, end, period_class, len(truth_sequence),
                         ".", 0, 0, ".", 0, 0, ".", ".", 0, 0, bucket))
            buckets[(period_class, bucket)] += 1
            continue
        # Both alleles wrong at two different lengths is two spurious lengths;
        # the tract is judged on the better-supported one and flagged in the
        # per-tract table by carrying both.
        low = spanning_reads(args.bam_low, contig, tract_first, tract_last)
        high = spanning_reads(args.bam_high, contig, tract_first, tract_last)
        support = {length: sum(1 for one in low if one.tract_length == length)
                   for length in spurious_lengths}
        spur_length = max(spurious_lengths, key=lambda one: support[one])
        carriers = [one for one in low if one.tract_length == spur_length]
        k, n = len(carriers), len(low)
        share_low = k / n if n else 0.0
        k_high = sum(1 for one in high if one.tract_length == spur_length)
        n_high = len(high)
        share_high = k_high / n_high if n_high else 0.0

        mix_reverse = (sum(1 for one in low if one.reverse) / n) if n else 0.0
        k_reverse = sum(1 for one in carriers if one.reverse)
        strand_clustered = (
            k >= 2
            and one_strand_probability(k_reverse, k, mix_reverse) < STRAND_P
        )
        family = Counter(
            (one.start, one.reverse, one.template) for one in carriers
        )
        family_size = max(family.values()) if family else 0
        family_clustered = k >= 2 and family_size >= 2 and family_size * 2 >= k
        clustered = strand_clustered or family_clustered
        persistent = (
            k_high >= PERSIST_MIN_READS
            and share_high >= PERSIST_SHARE_FACTOR * share_low
        )
        bucket = ("unseen_in_raw" if k == 0
                  else "both" if clustered and persistent
                  else "clustered" if clustered
                  else "locus_real" if persistent
                  else "sampling_noise")
        buckets[(period_class, bucket)] += 1
        share_bins[(
            period_class,
            "under 1 read in 10" if share_low < 0.10
            else "1 to 2 reads in 10" if share_low < 0.20
            else "2 to 3 reads in 10" if share_low < 0.30
            else "3 reads in 10 or more",
        )] += 1
        rows.append((
            contig, start, end, period_class, len(truth_sequence),
            "/".join(str(one) for one in spurious_lengths), n, k,
            f"{share_low:.3f}", n_high, k_high, f"{share_high:.3f}",
            f"{k_reverse}of{k}rev", family_size,
            int(clustered), bucket,
        ))
        if len(cases.setdefault(bucket, [])) < 3:
            lines = [
                f"=== {bucket}: {contig}:{tract_first}-{tract_last} "
                f"({period_class}, period {period})",
                f"    truth homozygous, {len(truth_sequence)} bases; called pair "
                f"{sorted(len(one) for one in pair)}; spurious length {spur_length}",
                f"    30x: {k} of {n} spanning reads carry it (share {share_low:.2f}); "
                f"300x: {k_high} of {n_high} (share {share_high:.2f})",
                f"    strand mix of the tract {1 - mix_reverse:.2f}F/{mix_reverse:.2f}R; "
                f"carriers {k - k_reverse}F/{k_reverse}R; "
                f"largest identical-signature family {family_size}",
                "    the carrier reads (start, strand, template length, tract length):",
            ]
            lines += [
                f"      {one.start:>12} {'R' if one.reverse else 'F'} "
                f"{one.template:>6} {one.tract_length:>4}"
                for one in carriers
            ]
            cases[bucket].append("\n".join(lines))

    with open(args.per_tract_out, "w", encoding="utf-8") as handle:
        handle.write(
            "contig\tstart\tend\tperiod_class\ttruth_len\tspurious_len\t"
            "n_low\tk_low\tshare_low\tn_high\tk_high\tshare_high\t"
            "carrier_strands\tlargest_family\tclustered\tbucket\n"
        )
        for row in rows:
            handle.write("\t".join(str(one) for one in row) + "\n")
    with open(args.cases_out, "w", encoding="utf-8") as handle:
        for bucket in ("both", "clustered", "unseen_in_raw", "locus_real",
                       "sampling_noise"):
            for case in cases.get(bucket, []):
                handle.write(case + "\n\n")

    print(f"\n=== the {len(found)} spurious heterozygotes, partitioned")
    for period_class in ("homopolymer", "period2plus"):
        for bucket in ("spelling_only", "unseen_in_raw", "clustered",
                       "locus_real", "both", "sampling_noise"):
            print(f"  {period_class:12s} {bucket:16s} "
                  f"{buckets[(period_class, bucket)]:5d}")
    print("\n=== the spurious allele's share of the tract's spanning reads at 30x")
    for period_class in ("homopolymer", "period2plus"):
        for label in ("under 1 read in 10", "1 to 2 reads in 10",
                      "2 to 3 reads in 10", "3 reads in 10 or more"):
            print(f"  {period_class:12s} {label:24s} "
                  f"{share_bins[(period_class, label)]:5d}")
    print(f"\nper-tract table: {args.per_tract_out}\ncases: {args.cases_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
