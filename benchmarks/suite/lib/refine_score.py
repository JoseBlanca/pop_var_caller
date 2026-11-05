"""Scores one call set of the refinement benchmark (bench/refine_sim.py) by group of samples.

usage: python3 refine_score.py <calls.vcf.gz> <truth.vcf.gz> <groups.tsv> <out.tsv> [<loci.tsv>]

`groups.tsv` names each sample's group (`sample<TAB>group`, e.g. `sim000<TAB>parents 30x`). Writes
`group, variant_class, stage, metric, value` lines; the benchmark turns them into result rows.

**Genotypes**, from pvc_eval's table of every sample at every locus (`compare --all-rows`), when
given: each row's outcome grouped into the suite's error types (lib/pvc_eval.py), by kind of locus
and by whether the sample had reads there — its largest DP over the call records in the locus. A
locus no call record overlaps is counted apart, so a missed variant is not taken for a sample
without reads.

**Phase**, from the calls and the truth directly, over SNP records only (every allele one base in
both), since an indel or a repeat can be spelled differently in the two
(doc/devel/ng/research/neighbour_evidence/phase_vcf_score.py, whose rule this is). A phased column
is right where the truth is heterozygous for the same two alleles and the truth's first copy carries
the allele written first. A *wrong link* is a link between two consecutive phased, truly
heterozygous records of one sample in one phase set whose rightness differs: one switch error.
"""

import bisect
import csv
import gzip
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
    """Per contig, the call records' positions in order, each with every sample's DP and, for
    SNP records, its alleles and every sample's (GT, PS)."""
    by_contig = defaultdict(lambda: ([], []))
    samples = []
    for kind, fields in records(path):
        if kind == "header":
            samples = fields
            continue
        keys = fields[8].split(":")
        dp_at = keys.index("DP") if "DP" in keys else None
        ps_at = keys.index("PS") if "PS" in keys else None
        alleles = [fields[3]] + fields[4].split(",")
        depths, phase = [], []
        for column in fields[9:]:
            values = column.split(":")
            dp = values[dp_at] if dp_at is not None and dp_at < len(values) else "."
            depths.append(int(dp) if dp.isdigit() else 0)
            ps = values[ps_at] if ps_at is not None and ps_at < len(values) else "."
            phase.append((values[0], ps))
        positions, rows = by_contig[fields[0]]
        positions.append(int(fields[1]))
        rows.append((alleles, depths, phase))
    return samples, by_contig


# How far before a locus a call record overlapping it can start: its REF's length at most. pvc
# writes no REF longer than its largest STR window, well under this.
LOOKBACK_BP = 1000


def overlapping(calls, chrom, start, end):
    """Indices of the call records whose REF overlaps `start..end` (1-based, inclusive)."""
    positions, rows = calls.get(chrom, ([], []))
    lo = bisect.bisect_left(positions, start - LOOKBACK_BP)
    hi = bisect.bisect_right(positions, end)
    return [i for i in range(lo, hi) if positions[i] + len(rows[i][0][0]) - 1 >= start]


def genotype_counts(loci_path, samples, calls, group_of):
    index = {name: i for i, name in enumerate(samples)}
    counts = Counter()
    with open(loci_path) as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            error = ERROR_TYPE.get(rec["outcome"], "not_scored")
            sample = rec["sample_b"]
            found = overlapping(calls, rec["chrom"], int(rec["start"]), int(rec["end"]))
            rows = calls[rec["chrom"]][1] if found else []
            if not found:
                reads = "no pvc record"
            elif max(rows[i][1][index[sample]] for i in found) > 0:
                reads = "with reads"
            else:
                reads = "without reads"
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
            gt, ps = call[2][i]
            if "|" not in gt or "." in gt or ps == ".":
                continue
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
    calls_path, truth_path, groups_path, out_path = sys.argv[1:5]
    loci_path = sys.argv[5] if len(sys.argv) > 5 else None
    group_of = {}
    with open(groups_path) as fh:
        for line in fh:
            sample, group = line.rstrip("\n").split("\t")
            group_of[sample] = group
    samples, calls = read_calls(calls_path)
    counts = phase_counts(truth_path, samples, calls, group_of)
    if loci_path:
        counts += genotype_counts(loci_path, samples, calls, group_of)
    # Every group also summed, as "all".
    for (group, cls, stage, metric), n in list(counts.items()):
        counts[("all", cls, stage, metric)] += n
    with open(out_path, "w") as out:
        out.write("group\tvariant_class\tstage\tmetric\tvalue\n")
        for key in sorted(counts):
            out.write("\t".join(key) + f"\t{counts[key]}\n")


if __name__ == "__main__":
    main()
