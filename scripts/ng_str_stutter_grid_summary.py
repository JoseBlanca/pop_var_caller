"""Gather a stuttering-benchmark grid into tables (doc/devel/ng/spec/str_stutter_benchmark.md §6, §9).

    python3 scripts/ng_str_stutter_grid_summary.py <grid-dir>

Reads every `<grid-dir>/<design>_<samples>s_<depth>x/score_cells.tsv` and `score_vcfeval.tsv`, and
writes `<grid-dir>/grid_cells.tsv`: one row per run, caller, segment set, motif length and outcome,
with the cell count.

Two segment sets: `all`, and `same_typing` — every segment but the AAAAAC, AAAAAG and AAAAAT
hexanucleotides, which main types as ordinary sequence (it lacks the compound-motif fix), so main
and the branch are compared on the same STR segments.

**The branch also scored with its `strBoundary` placeholders treated as absent** (callers suffixed
`_no_placeholders`). Where a generic edit crosses a segment's edge at a stutter weight of 0.95 or
more, the branch writes the STR answer and marks the flank bases that edit touched as missing for
every sample (`tract_dual_path.md` §8). The scorer counts such a cell as uncalled, which is what the
VCF says; this second count shows what the rest of the branch's answer is worth. It reuses the
scorer's own functions, so the two counts differ only in the records read.

Also writes `<grid-dir>/grid_vcfeval.tsv`: vcfeval's counts summed over each run's samples.
"""

import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ng_str_stutter_score as score  # noqa: E402

RUN = re.compile(r"^(?P<design>.+)_(?P<samples>\d+)s_(?P<depth>\d+)x$")
MAIN_TYPES_DIFFERENTLY = {"AAAAAC", "AAAAAG", "AAAAAT"}
WITH_PLACEHOLDERS = ["ng_defaults", "ng_fitted", "ng_true_stutter"]


def runs(grid):
    for name in sorted(os.listdir(grid)):
        match = RUN.match(name)
        if match and os.path.exists(os.path.join(grid, name, "score_cells.tsv")):
            yield name, match["design"], int(match["samples"]), int(match["depth"])


def without_placeholders(run, caller, samples):
    """The scorer's records for `caller` minus every `strBoundary` placeholder."""
    path = os.path.join(run, dict(score.CALLERS)[caller])
    kept = os.path.join(run, f"{caller}.no_placeholders.vcf")
    with open(path) as source, open(kept, "w") as target:
        for line in source:
            if line.startswith("#") or line.split("\t", 7)[6] != "strBoundary":
                target.write(line)
    calls = score.CallIndex(score.read_calls(kept, caller, samples))
    os.remove(kept)
    return calls


def main():
    grid = sys.argv[1]
    counts = collections.Counter()
    for name, design, samples_n, depth in runs(grid):
        run = os.path.join(grid, name)
        with open(os.path.join(run, "score_cells.tsv")) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                key = (design, samples_n, depth, row["caller"], row["period"], row["truth"], row["outcome"])
                counts[key + ("all",)] += 1
                if row["motif"] not in MAIN_TYPES_DIFFERENTLY:
                    counts[key + ("same_typing",)] += 1

        contigs = score.read_fasta(os.path.join(run, "reference.fa"))
        segments, genotypes = score.read_truth(run)
        with open(os.path.join(run, "libraries.tsv")) as handle:
            handle.readline()
            samples = [line.split("\t")[0] for line in handle]
        for caller in WITH_PLACEHOLDERS:
            if not os.path.exists(os.path.join(run, dict(score.CALLERS)[caller])):
                continue
            calls = without_placeholders(run, caller, samples)
            for segment in segments.values():
                for sample in samples:
                    truth_pair = genotypes[segment.name][sample]
                    outcome, _ = score.score_cell(caller, segment, sample, truth_pair, calls, contigs)
                    truth = "variant" if any(index != 0 for index in truth_pair) else "reference"
                    key = (design, samples_n, depth, caller + "_no_placeholders", str(segment.period),
                           truth, outcome)
                    counts[key + ("all",)] += 1
                    if segment.motif not in MAIN_TYPES_DIFFERENTLY:
                        counts[key + ("same_typing",)] += 1
        print(f"[summary] {name}", file=sys.stderr)

    with open(os.path.join(grid, "grid_cells.tsv"), "w") as out:
        out.write("design\tsamples\tdepth\tcaller\tperiod\ttruth\toutcome\tsegments\tcells\n")
        for key in sorted(counts):
            out.write("\t".join(map(str, key)) + f"\t{counts[key]}\n")

    vcfeval = collections.Counter()
    for name, design, samples_n, depth in runs(grid):
        path = os.path.join(grid, name, "score_vcfeval.tsv")
        if not os.path.exists(path):
            continue
        with open(path) as handle:
            handle.readline()
            for line in handle:
                caller, _, stretch, tp, fp, fn = line.rstrip("\n").split("\t")
                for what, value in (("true_positives", tp), ("false_positives", fp), ("false_negatives", fn)):
                    vcfeval[(design, samples_n, depth, caller, stretch, what)] += int(float(value))
    with open(os.path.join(grid, "grid_vcfeval.tsv"), "w") as out:
        out.write("design\tsamples\tdepth\tcaller\tstretch\tcount\tvalue\n")
        for key in sorted(vcfeval):
            out.write("\t".join(map(str, key)) + f"\t{vcfeval[key]}\n")


if __name__ == "__main__":
    main()
