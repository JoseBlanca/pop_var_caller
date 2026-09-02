# ng's first calling benchmark — against the production caller and freebayes

**Date:** 2026-09-02.
**Callers compared:** ng (`pop_var_caller_exp call-from-alignments`, this
repository, run with `--defaults`), the production caller
(`pop_var_caller pileup` → `var-calling`, the `high-recall` preset), and
freebayes 1.3.6. GATK 4 appears in one table where a result already existed.

ng reached a working alignments-to-VCF path on 2026-09-01. This is the first
time its calls have been scored.

---

## 1. The headline

**On human data, ng recovers at least as many truth SNPs as freebayes at every
depth from 5× to 300×, and about five points fewer than the production caller
from 30× on. Its indel recall is nine to thirteen points below both. Nine in
ten of the truth variants it misses are inside tandem repeats, where it builds
no locus at all.**

ng does not
yet have a repeat-tract calling path: every tract in the ground it is asked
for is counted as ground it cannot speak for and named in its run report — 6
bases in every 100 of the GIAB confident regions, 5 in every 100 of the tomato
benchmark. A truth variant inside a tract is therefore out of reach by
construction, and the accuracy tables score it as a false negative exactly as
they would score a variant ng examined and got wrong. Section 3 separates the
two by measurement rather than by assertion.

Two other things are true of every ng number below and of none of the others':

- **Nothing is fitted.** No command writes a fitted parameters file yet, so
  every run is `--defaults` — no base-quality calibration, no contamination,
  no inbreeding coefficient. ng's own run report says so: *0 of 7 groups the
  file says were fitted*.
- **A cohort of one is ng's hardest case**, and the per-sample human benchmark
  is exactly that. Its allele frequencies come from the run itself, and at one
  sample there is nothing to pool.

---

## 2. Human accuracy across read depth — GIAB per-sample

Three samples (HG002, HG003, HG004), each on its own random 100-region
confident BED with its own GIAB truth VCF, at seven depths from 5× to 300×.
Counts are summed over the three samples. Scoring is allele concordance
(POS+REF+ALT after left-alignment and biallelic split), truth restricted to
`FILTER=PASS`. Every caller's output is gated at QUAL ≥ 30 — the production
caller by its own `--min-qual` default, freebayes and ng by a post-call filter
in their runners.

Source: `benchmarks/giab/results/per_sample/freebayes_comparison.tsv`, written
by `benchmarks/giab/src/freebayes_comparison_dashboard.py`.

### SNPs — recall

| depth | ng | production | freebayes |
|---|---:|---:|---:|
| 5× | 0.664 | 0.706 | 0.585 |
| 10× | 0.852 | 0.913 | 0.789 |
| 15× | 0.906 | 0.960 | 0.881 |
| 20× | 0.925 | 0.974 | 0.912 |
| 30× | 0.935 | 0.987 | 0.927 |
| 50× | 0.940 | 0.988 | 0.935 |
| 300× | 0.938 | 0.985 | 0.936 |

ng recovers more truth SNPs than freebayes at every depth to 30×, and the two
are within half a point of each other from 50× on. The production caller is
ahead of both throughout, by 5 points at 30× and above.

ng's recall stops climbing at about 0.94 while the production caller's goes on
to 0.99. A ceiling that does not move with depth is not a sampling limit; §3
says what it is.

### SNPs — precision, and the false positives behind it

| depth | ng | production | freebayes |
|---|---:|---:|---:|
| 5× | 0.9935 | 0.9918 | 0.9975 |
| 30× | 0.9918 | 0.9941 | 0.9984 |
| 300× | 0.9877 | 0.9951 | 0.9984 |

All three are within seven parts in a thousand of each other at every depth.
In counts, over 2 061 truth SNPs at 300×: ng emits 24 SNPs the truth set does
not carry, the production caller 10, freebayes 3.

### Indels — recall

| depth | ng | production | freebayes |
|---|---:|---:|---:|
| 5× | 0.455 | 0.518 | 0.512 |
| 10× | 0.673 | 0.767 | 0.694 |
| 15× | 0.755 | 0.876 | 0.824 |
| 20× | 0.785 | 0.921 | 0.861 |
| 30× | 0.818 | 0.930 | 0.912 |
| 50× | 0.818 | 0.946 | 0.924 |
| 300× | 0.815 | 0.939 | 0.946 |

ng is last at every depth. Its curve flattens at 0.82 from 30× on; the other
two are still at 0.93 and 0.91 there and reach 0.94 and 0.95 by 300×. Indel
precision is not what is costing it: 0.989 at 30× and 1.000 at 300×, against
the production caller's 0.997 at both.

Indels are where the missing repeat-tract path costs most, which is what a
geneticist would expect — short insertions and deletions concentrate in
homopolymers and short tandem repeats.

### Genotypes at the sites each caller found

Everything above scores whether a caller found the site. This scores whether it
got the number of copies right, over the true positives only: at each site
present in both the truth set and the caller's output, does the called genotype
match the truth's (phase ignored, so `1|0` counts as `0/1`)? Per sample, at
30×:

| class | ng | production | freebayes |
|---|---|---|---|
| SNPs | 99.9 / 100 / 100 % | 99.9 / 100 / 100 % | 99.7 / 99.3 / 99.2 % |
| indels | 78.6 / 81.1 / 78.5 % | 80.3 / 79.5 / 77.2 % | 95.3 / 97.7 / 98.9 % |

(HG002 / HG003 / HG004.)

ng's SNP genotypes are as good as the production caller's and slightly better
than freebayes'. Its indel genotypes are wrong about **one time in five**,
which is where the production caller sits too, while freebayes is wrong about
one time in fifty. ng shares this defect with the production caller rather than
introducing it, and it is separate from the repeat-tract recall gap — every
site counted here is one the caller did find.

The same panel is in the dashboard, boxed over the three samples across every
depth tier.

---

## 3. Where the misses come from

A missed truth variant can mean two different things, and the accuracy table
scores them alike:

- the caller built a locus at that position and did not call the variant — a
  genotyping miss; or
- it never built a locus there, so nothing about its genotyping is being
  measured.

To separate them, `benchmarks/giab/src/ng_missed_sites_probe.sh` takes each
caller's missed truth sites at 300×, writes them as a one-base-per-site BED,
and runs ng over exactly those bases. ng's own run report then says how many
loci it built there.

| caller | class | missed | ng builds a locus | ng builds nothing |
|---|---|---:|---:|---:|
| ng | SNPs | 128 | 11 | 117 (91%) |
| ng | indels | 53 | 1 | 52 (98%) |
| production | SNPs | 30 | 19 | 11 (37%) |
| production | indels | 20 | 4 | 16 (80%) |

**Of the 128 truth SNPs ng misses at 300×, it builds a locus at 11.** For
indels it builds one locus across 53 missed sites. So ng's recall ceiling is
the unbuilt repeat-tract path and not the genotyper: on the ground it does
build, at 300×, it recovers 1 933 of the roughly 1 944 truth SNPs within
reach.

The locus count is ng's in every row, so the production caller's rows read
differently: they say how much of *its* residual miss list also lies on ground
ng cannot reach. 16 of its 20 missed indel sites do — the two callers are
failing on the same ground, and the production caller simply reaches more of
it.

Source: `benchmarks/giab/results/per_sample/ng_missed_sites.tsv`.

**Followed up the same day**, in
[`ng_str_path_losses_2026-09-02.md`](ng_str_path_losses_2026-09-02.md), which splits every
truth variant by the kind of ground it sits on. Two of its findings change how the tables above
read:

- **On ordinary ground ng's indel recall is the best of the three** — 0.982 at 30× and 50×,
  against the production caller's 0.945 and 0.953 and freebayes' 0.931 and 0.938. Its indel
  recall is nine to thirteen points below both *overall* and above both *where it calls*.
- **The run routes about five times more ground to the unbuilt repeat path than ng's own
  calling policy says belongs there**, because it classifies with the floors the catalog file
  was stored at. Four in five of the lost variants are on ground the calling floors would leave
  on the generic path.

---

## 4. One high-coverage sample over 5 Mb — the bottle benchmark

A single HG002 CRAM restricted to 1 000 confident intervals (~5 Mb of
GRCh38), scored against the GIAB benchmark VCF for that region set. This
benchmark deliberately runs **every caller with its QUAL gate at zero**, so
each emits its full low-confidence tail and the dashboard sweeps a common
cutoff. The precision figures below are therefore the ungated ones and are not
comparable with §2's.

| caller | class | TP | FP | FN | precision | recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| production | snps | 6 687 | 1 403 | 89 | 0.827 | 0.987 | 0.900 |
| gatk | snps | 6 756 | 123 | 20 | 0.982 | 0.997 | 0.990 |
| freebayes | snps | 6 501 | 8 650 | 275 | 0.429 | 0.959 | 0.593 |
| ng | snps | 6 400 | 484 | 376 | 0.930 | 0.945 | 0.937 |
| production | indels | 1 116 | 18 | 89 | 0.984 | 0.926 | 0.954 |
| gatk | indels | 1 204 | 13 | 1 | 0.989 | 0.999 | 0.994 |
| freebayes | indels | 1 154 | 10 | 51 | 0.991 | 0.958 | 0.974 |
| ng | indels | 982 | 8 | 223 | 0.992 | 0.815 | 0.895 |

Two things this adds to §2. ng's recall gap is the same size on a region set
nine times larger — SNP recall 0.945 here against 0.938 there, indel recall
0.815 in both — so it is a property of the caller and not of the 100 regions
§2 happens to use. And with no QUAL gate at all,
ng emits 484 SNPs the truth set does not carry against the production caller's
1 403 and freebayes' 8 650: ng's low-confidence tail is much shorter than
either's.

Source: `benchmarks/human_genome_bottle/results/comparison/accuracy.tsv`,
written by `benchmarks/lib/compare_to_truth.sh`.

---

## 5. A 63-accession tomato cohort

63 *S. lycopersicum* accessions at about three reads a position, over the 80
intervals of `benchmarks/tomato1/regions.bed` (8 Mb of SL4.0). There is no
truth set here, so nothing below is an accuracy claim. What it tests is
whether ng runs a real cohort at all, and how its callset relates to the
production caller's on the same samples and the same ground.

Both callers were re-run on all 63 accessions for this comparison; the result
files that were there dated from May and held 26 samples (production) and 18
(freebayes).

**ng called the cohort.** One process, every sample's CRAM held open, CRAMs to
VCF in 295 s wall clock. Its run report accounts for the whole 8 Mb: 94.4%
called, 5.4% repeat tract it does not build, 0.1% tandem array too long to type
as callable. 249 loci were declined for being wider than
`--max-cohort-locus-span`.

The production caller's cohort step took 77 s, but that is not the comparable
number: it reads 63 pre-built `.psp` files, 62 of which were already on disk
from an earlier run and were not rebuilt or timed here. The one sample that had
to be piled up took 11 s. So the production path's CRAMs-to-VCF cost is roughly
those 77 s plus 63 pileups, and this run does not measure it.

| | production | ng |
|---|---:|---:|
| records at QUAL ≥ 30 | 189 933 | 206 873 |
| of those, carrying an indel allele | 2 838 | 4 664 |
| distinct ALT alleles | 198 974 | 230 195 |
| heterozygous genotype calls per sample per kb | 1.00 | 1.23 |
| genotypes left as no-call | 0 | 627 856 of 13.0 M (4.8%) |

**They agree on 178 464 ALT alleles** — 89.7% of the production caller's set
and 77.5% of ng's. ng carries 51 731 alleles the production caller does not
(43 235 SNPs, 8 496 indels); the production caller carries 20 510 ng does not
(17 156 SNPs, 3 354 indels).

**ng calls heterozygotes 23% more often**: 1.23 per kb per sample against 1.00,
on the same accessions and effectively the same ground. With no truth set that
is a difference and not an error, but it is the direction to be suspicious of.
Tomato accessions are largely inbred, and a `--defaults` run assumes an
inbreeding coefficient of zero — the assumption that most favours calling a
heterozygote. The production run in this comparison also fitted the coefficient
at zero, so that alone does not separate them, and neither number can be
checked against a truth set here.

**freebayes is not in this comparison, and getting it there is a job of its
own.** In one process on 63 accessions it advances at about 33 kb of reference
a minute, so the 8 Mb is roughly four hours. Splitting it into twelve
processes, one per chromosome, did not deliver: every freebayes process opens
every CRAM, and in the 16 GB dev container four of the twelve were OOM-killed,
leaving a concatenated file that looked ordinary and was missing whole
chromosomes. Both attempts are on disk under `INVALID_`/`INCOMPLETE_` names so
nothing reads them as a result. `run_freebayes.sh` now refuses a run with a
dead shard rather than concatenating what survived; the shard count has to be
chosen for the memory available. For scale, ng called the same cohort and
ground in 295 s.

Source: `benchmarks/tomato1/results/{ng,ours}/`, and
`benchmarks/tomato1/scripts/dashboard.py` for the agreement view.

---

## 6. What these benchmarks do not say

- **Nothing here scores genotypes at the cohort level.** The human scoring is
  allele concordance — a site called with the wrong number of copies still
  counts as a true positive. The per-coverage dashboard has a separate
  genotype-concordance panel for the human samples; the tomato cohort has no
  truth set at all.
- **No ng number here comes from fitted parameters**, because no command
  writes them yet. Every improvement the fit is meant to buy — base-quality
  calibration, contamination, inbreeding — is absent from all of it.
- **Repeat tracts are unscored ground for ng, not merely hard ground.** Any
  comparison of the repeat-tract fraction of a genome between species changes
  what these numbers mean: 6 bases in 100 of the GIAB confident regions, 5 in
  100 of the tomato benchmark, and both are small compared with a whole
  genome.
- **The depth range covered is 5× to 300× on human and about 3× on tomato.**
  The tomato cohort is the only low-depth multi-sample evidence, and it has no
  truth set.

---

## 7. What was added

Runners:

- `benchmarks/lib/run_ng.sh` — ng for any benchmark with a `bench.config.sh`,
  in `single` or `cohort` mode. Builds the tandem-repeat catalog if it is
  missing, and applies the benchmark's `MIN_QUAL` gate the way the freebayes
  runner does.
- `benchmarks/giab/src/run_ng_per_sample.sh` — ng over the GIAB per-sample
  coverage tiers, matching the existing per-sample runners for the other two
  callers.
- `benchmarks/giab/src/ng_missed_sites_probe.sh` — the §3 measurement.

Wiring:

- `benchmarks/lib/compare_to_truth.sh` scores `results/ng/` alongside the
  other three by default.
- `benchmarks/human_genome_bottle/bench.config.sh` and
  `benchmarks/tomato1/bench.config.sh` name a writable catalog path, since
  both references sit on a read-only mount.

Dashboards:

- `benchmarks/giab/src/freebayes_comparison_dashboard.py` — ng added as a
  fourth arm across every table and chart, plus the §3 section.
- `benchmarks/lib/comparison_dashboard.py` — ng given a stable colour; it was
  already generic over whatever `accuracy.tsv` holds.
- `benchmarks/tomato1/scripts/dashboard.py` — reworked from three fixed
  callers to a selectable set, with per-caller totals, pairwise overlap, and a
  sample-count column that makes a stale result file visible.
