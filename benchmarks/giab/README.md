Human benchmarks based on the Genome in a Bottle dataset (GIAB)

In these directories there are files for two different datasets:

1. per_sample: SNP calling of three samples/individual for 100 regions. 
2. mendelian: SNP calling of a trio family for 100 regions.

per_sample dataset:

- The regions were selected at random
- The regions are different for each sample.
- Directory: benchmarks/giab/per_sample
- The vcf dir includes the truth downloaded from GIAB
- In the bam dir there are BAMs with different coverages

mendelian dataset:
- The regions were selected at random
- The regions are the same for all samples.
- Directory: benchmarks/giab/mendelian
- The vcf dir includes the truth downloaded from GIAB
- HG002 is the child and HG003 and HG004 are the parents

Runners and dashboards (benchmarks/giab/src/)
---------------------------------------------

Three callers are run over the per_sample dataset, each restricted to the
sample's own confident BED, and scored against that sample's GIAB truth VCF:

- run_ours_per_sample.sh   the production caller (pileup -> .psp -> var-calling)
- run_freebayes_per_sample.sh
- run_ng_per_sample.sh     ng, the experimental caller, alignments -> VCF in
                           one process

Each takes a coverage tier and an optional sample list, e.g.

    benchmarks/giab/src/run_ng_per_sample.sh 30x HG002

Two things to hold against ng's numbers, both of which cost it recall and not
precision. It does not build loci inside tandem repeats yet, so a truth variant
inside one is unreachable for it and scores as a false negative. And nothing is
fitted: no command writes a fitted parameters file yet, so every run is
`--defaults`, with no base-quality calibration, contamination or inbreeding
coefficient.

**How much ground the first of those costs changed on 2026-09-02, and by about
seven times.** A run used to ask the repeat catalog with the floors the *file*
was stored at, so every row it held became a repeat tract of the run — 6 bases
in every 100 of these confident regions. A run now routes on the floors it is
given, defaulting to ng's measured stutter onsets, and a candidate below them is
ordinary sequence that the SNP/indel caller handles. On HG002 that is 32,577
bases against 4,930, and pooled over the three samples at 30x it took SNP recall
from 0.935 to 0.974 and indel recall from 0.818 to 0.946, for one extra false
positive. The five flags that set it are under "What counts as a repeat" in
`call-from-alignments --help`, and every run records what it used in the
parameters file it writes.

    benchmarks/giab/src/score_ng_recall.sh 30x ng

scores one results directory the way accuracy_dashboard.py does — per sample and
class, TP/FP/FN with recall and precision. Two directories under the same
coverage can be scored and compared, which is what a routing change needs.

ng_missed_sites_probe.sh separates those two kinds of miss. It takes each
caller's missed truth sites, hands ng exactly those bases, and reads from ng's
run report how many loci it built there — so a miss at a site nothing was built
for is told apart from a miss at a site the caller looked at and got wrong. It
writes results/per_sample/ng_missed_sites.tsv.

Dashboards (marimo; `uv run marimo run <file>`):

- freebayes_comparison_dashboard.py  all three callers, per coverage tier,
                                     precision/recall/F1 + QUAL and genotype
                                     concordance. Writes
                                     results/per_sample/freebayes_comparison.tsv.
- accuracy_dashboard.py              the production caller's presets alone.
- allele_balance_dashboard.py, mapq_depth_dashboard.py — single-issue studies.
