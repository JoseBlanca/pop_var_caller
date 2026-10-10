"""Right variant genotypes and false calls, per caller and run, from each run's score_cells.tsv.

    python3 scripts/ng_str_accuracy/compare.py <caller,caller,...> <run> [<run> ...]

One row per run: the variant and reference pairs (one sample at one STR segment), then for each
caller the pairs called right and the reference pairs given a false call. The 262 hexanucleotide
segments of motif AAAAAC, AAAAAG or AAAAAT are left out, as in every number of
doc/devel/reports/implementations/str_accuracy_2026-10-10.md, so runs compare with that report.

Env: STR_RUNS_DIR (default <repo>/tmp/str_accuracy/runs).
"""
import collections
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
RUNS = os.environ.get("STR_RUNS_DIR", os.path.join(ROOT, "tmp", "str_accuracy", "runs"))
LEFT_OUT = {"AAAAAC", "AAAAAG", "AAAAAT"}

callers = sys.argv[1].split(",")
print("run\tvariant_pairs\treference_pairs\t" + "\t".join(f"{c}_right\t{c}_false" for c in callers))
for run in sys.argv[2:]:
    path = os.path.join(RUNS, run, "score_cells.tsv")
    if not os.path.exists(path):
        print(f"{run}\tno score_cells.tsv", file=sys.stderr)
        continue
    counts = collections.Counter()
    variant = reference = 0
    with open(path) as handle:
        handle.readline()
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[10] in LEFT_OUT:
                continue
            counts[(f[0], f[6], f[8])] += 1
            if f[0] == "gatk":
                variant += f[6] == "variant"
                reference += f[6] == "reference"
    print(f"{run}\t{variant}\t{reference}\t" + "\t".join(
        f"{counts[(c, 'variant', 'right')]}\t{counts[(c, 'reference', 'false_call')]}" for c in callers))
