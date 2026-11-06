"""Scores one call set of the refinement benchmark (bench/refine_sim.py) by group of samples.

usage: python3 refine_score.py <calls.vcf.gz> <truth.vcf.gz> <groups.tsv> <coverage dir> <out.tsv>
                               [<loci.tsv>]

`groups.tsv` names each sample's group (`sample<TAB>group`, e.g. `sim000<TAB>parents 30x`). The
coverage directory holds `<sample>.bed`, the stretches of reference each sample's mapped reads
cover (`refine_sim`'s coverage step). Writes `group, variant_class, stage, metric, value` lines;
the benchmark turns them into result rows.

**Genotypes**, from pvc_eval's table of every sample at every locus (`compare --all-rows`), when
given: each row's outcome grouped into the suite's error types (lib/pvc_eval.py), by kind of locus
and by whether the sample had reads there — any mapped read over any base of the locus. Taken from
the alignments, not from the calls, so every caller is split alike, those that write no DP
(Beagle, STITCH) included.

**Phase**, from the calls and the truth directly, over SNP records only (every allele one base in
both), since an indel or a repeat can be spelled differently in the two
(doc/devel/ng/research/neighbour_evidence/phase_vcf_score.py, whose rule this is). A phased column
is right where the truth is heterozygous for the same two alleles and the truth's first copy carries
the allele written first. A *wrong link* is a link between two consecutive phased, truly
heterozygous records of one sample in one phase set whose rightness differs: one switch error. A
phased genotype with no PS belongs to one phase set spanning its contig, as the VCF specification
has it; that is how Beagle writes phase.
"""

import bisect
import csv
import gzip
import os
import sys
from collections import Counter, defaultdict

ERROR_TYPE = {
    "only in B": "wrong_gt_truth_ref",
    "only in A": "wrong_gt_truth_var",
    "same alleles, other genotype": "wrong_gt_truth_var",
    "different": "wrong_gt_truth_var",
    "missing in B": "no_call",
    "same genotype": "right",
    "reference in both": "right",
}


def records(path):
    """('header', samples) once, then ('record', fields) for every record."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                yield "header", fields[9:]
            else:
                yield "record", fields


def read_calls(path):
    """Per contig, the call records' positions in order, each with its alleles and every sample's
    (GT, PS)."""
    by_contig = defaultdict(lambda: ([], []))
    samples = []
    for kind, fields in records(path):
        if kind == "header":
            samples = fields
            continue
        keys = fields[8].split(":")
        ps_at = keys.index("PS") if "PS" in keys else None
        alleles = [fields[3]] + fields[4].split(",")
        phase = []
        for column in fields[9:]:
            values = column.split(":")
            ps = values[ps_at] if ps_at is not None and ps_at < len(values) else "."
            phase.append((values[0], ps))
        positions, rows = by_contig[fields[0]]
        positions.append(int(fields[1]))
        rows.append((alleles, phase))
    return samples, by_contig


def read_coverage(directory, samples):
    """Per sample and contig, the covered stretches as (starts, ends), 0-based half-open, sorted
    and disjoint."""
    coverage = {}
    for sample in samples:
        by_contig = defaultdict(lambda: ([], []))
        with open(os.path.join(directory, f"{sample}.bed")) as fh:
            for line in fh:
                chrom, start, end = line.split("\t")[:3]
                starts, ends = by_contig[chrom]
                starts.append(int(start))
                ends.append(int(end))
        coverage[sample] = by_contig
    return coverage


def has_reads(coverage, sample, chrom, start, end):
    """Whether any read of `sample` covers a base of `start..end` (1-based, inclusive)."""
    starts, ends = coverage[sample].get(chrom, ([], []))
    i = bisect.bisect_left(starts, end) - 1  # the last stretch starting before the locus ends
    return i >= 0 and ends[i] > start - 1


def genotype_counts(loci_path, coverage, group_of):
    counts = Counter()
    with open(loci_path) as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            error = ERROR_TYPE.get(rec["outcome"], "not_scored")
            sample = rec["sample_b"]
            covered = has_reads(coverage, sample, rec["chrom"], int(rec["start"]), int(rec["end"]))
            reads = "with reads" if covered else "without reads"
            counts[(group_of[sample], f"{rec['kind']}, {reads}", "error type", error)] += 1
    return counts


def phase_counts(truth_path, samples, calls, group_of):
    counts = Counter()
    truth_samples = []
    phased = defaultdict(list)  # sample -> [(contig, position, right, phase set)]
    for kind, fields in records(truth_path):
        if kind == "header":
            truth_samples = fields
            continue
        true_alleles = [fields[3]] + fields[4].split(",")
        if any(len(a) != 1 for a in true_alleles):
            continue
        positions, rows = calls.get(fields[0], ([], []))
        at = bisect.bisect_left(positions, int(fields[1]))
        call = rows[at] if at < len(positions) and positions[at] == int(fields[1]) else None
        if call is not None and any(len(a) != 1 for a in call[0]):
            call = None
        truth_column = dict(zip(truth_samples, fields[9:]))
        for i, sample in enumerate(samples):
            true_pair = [true_alleles[int(a)] for a in truth_column[sample].split("|")]
            if true_pair[0] == true_pair[1]:
                continue
            group = group_of[sample]
            counts[(group, "SNP", "phase", "true_heterozygotes")] += 1
            if call is None:
                continue
            gt, ps = call[1][i]
            if "|" not in gt or "." in gt:
                continue
            if ps == ".":
                ps = "contig"  # no PS: the contig's one phase set (VCF 4.2, the PS field)
            first, second = (call[0][int(a)] for a in gt.split("|"))
            if sorted([first, second]) != sorted(true_pair):
                counts[(group, "SNP", "phase", "phased_wrong_genotype")] += 1
                continue
            counts[(group, "SNP", "phase", "phased_right_genotype")] += 1
            phased[sample].append((fields[0], int(fields[1]), first == true_pair[0], ps))
    for sample, rows in phased.items():
        group = group_of[sample]
        for (c0, _, r0, b0), (c1, _, r1, b1) in zip(rows, rows[1:]):
            if c0 == c1 and b0 == b1:
                counts[(group, "SNP", "phase", "links_in_phase_sets")] += 1
                counts[(group, "SNP", "phase", "wrong_links")] += r0 != r1
    return counts


def main():
    calls_path, truth_path, groups_path, coverage_dir, out_path = sys.argv[1:6]
    loci_path = sys.argv[6] if len(sys.argv) > 6 else None
    group_of = {}
    with open(groups_path) as fh:
        for line in fh:
            sample, group = line.rstrip("\n").split("\t")
            group_of[sample] = group
    samples, calls = read_calls(calls_path)
    counts = phase_counts(truth_path, samples, calls, group_of)
    if loci_path:
        counts += genotype_counts(loci_path, read_coverage(coverage_dir, samples), group_of)
    # Every group also summed, as "all".
    for (group, cls, stage, metric), n in list(counts.items()):
        counts[("all", cls, stage, metric)] += n
    with open(out_path, "w") as out:
        out.write("group\tvariant_class\tstage\tmetric\tvalue\n")
        for key in sorted(counts):
            out.write("\t".join(key) + f"\t{counts[key]}\n")


if __name__ == "__main__":
    main()
