"""Deletions across a segment's edge: could the caller's records spell them, and were carriers right?

    python3 scripts/ng_str_accuracy/edge_deletions.py <run-dir> <caller>

For every segment whose truth has an allele of kind `deletion_across_edge` (the simulator's
deletion running from the repeat into a flank), every combination of the alleles of the caller's
records over the segment's window is spelled, and the deleted sequence is looked for among them.
"Never offered" means no combination spells it — the caller had no way to be right. Then every
carrier of that allele is counted by its scored outcome. This is the measurement behind
doc/devel/reports/implementations/str_accuracy_2026-10-10.md §5's "never offered" count.
"""
import collections
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ng_str_stutter_score as s  # noqa: E402

run, caller = sys.argv[1], sys.argv[2]
contigs = s.read_fasta(os.path.join(run, "reference.fa"))
segments, genotypes = s.read_truth(run)
samples = [line.split("\t")[0] for line in list(open(os.path.join(run, "libraries.tsv")))[1:]]
calls = s.CallIndex(s.read_calls(os.path.join(run, caller + ".vcf"), caller, samples))
kind = collections.defaultdict(dict)
with open(os.path.join(run, "segment_alleles.tsv")) as handle:
    handle.readline()
    for line in handle:
        f = line.rstrip("\n").split("\t")
        kind[f[0]][int(f[9])] = f[11]
outcome = {}
with open(os.path.join(run, "score_cells.tsv")) as handle:
    handle.readline()
    for line in handle:
        f = line.rstrip("\n").split("\t")
        if f[0] == caller:
            outcome[(f[1], f[2])] = f[8]

segment_groups, pair_groups = collections.Counter(), collections.Counter()
for name, segment in segments.items():
    edge = [allele for allele, k in kind[name].items() if k == "deletion_across_edge"]
    if not edge or segment.motif in {"AAAAAC", "AAAAAG", "AAAAAT"}:
        continue
    reference = contigs[segment.chrom]
    records = calls.overlapping(segment.chrom, segment.window_start, segment.window_end)
    start = min([segment.window_start] + [r.start for r in records])
    end = max([segment.window_end] + [r.end for r in records])
    offered = set()
    for choice in itertools.product(*[range(len(r.alleles)) for r in records]):
        edits = [(r.start, r.end, r.alleles[c]) for r, c in zip(records, choice)
                 if c != 0 and r.alleles[c] != "*"]
        spelled = s.spell(reference, start, end, edits)
        if spelled:
            offered.add(spelled)
    for allele in edge:
        carriers = [smp for smp in samples if allele in genotypes[name][smp]]
        if not carriers:
            continue
        sequence = reference[start:segment.window_start] + segment.alleles[allele][0] \
            + reference[segment.window_end:end]
        key = "offered" if sequence in offered else "never offered"
        segment_groups[key] += 1
        for smp in carriers:
            pair_groups[(key, "right" if outcome[(name, smp)] == "right" else "not right")] += 1
print("segments with a carried edge-deletion allele:", dict(segment_groups))
print("carrier pairs:", dict(sorted(pair_groups.items())))
