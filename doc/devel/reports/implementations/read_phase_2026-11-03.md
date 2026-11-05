# Phase from reads: what was built and what it measures (2026-11-03)

*Implementation report for [`neighbour_evidence_read_phase.md`](../../ng/spec/neighbour_evidence_read_phase.md),
built on branch `neighbour-pass` in eleven commits, `850eff28` to `a5941dc1`. Not merged to main.*

## What it does, in one paragraph

A read pair comes from one molecule, so where it covers two heterozygous sites of a sample, the
alleles it shows there lie on one chromosome copy. pvc now keeps that evidence. When the cohort
merge builds a locus that becomes a VCF record, it notes which allele each of a sample's reads
showed there (the *read alleles*, at most 32 a sample). On the thread that writes records in
genome order, the reads seen at two records or more within 1,000 bases are kept as numbered
*fragments*, and the rest dropped. The fragments travel beside each record to whichever step
writes phase — in memory, or through the hidden-duplication filter's spill file — and are never
written to any file. Where the neighbour-evidence pass (*refinement*) does not run, at ploidy two,
the writer phases from reads alone: `GT` with `|` and a `PS` on every line. Where it runs, reads
break the population links they contradict, link records the population does not, and attach
heterozygotes the population cannot orient; and the log reports how often reads contradict the
population links they reach (the *read check*).

## What was built, step by step

| step | commit | what |
|---|---|---|
| 1 | `850eff28` | `pvc refine-genotypes` removed; `--neighbour-evidence` renamed `--refine-genotypes`; research scripts call calling instead |
| 2 | `9ca59984` | every minted record carries its walked region's first chain id (`chain_id_base`); psp format 1.1 stores it in the record head, and a 1.0 psp is refused with the message to regenerate it; mates pair by read group and read name |
| 3 | `ecf3b021` | the capture, inside the merge's own allele derivation; measured, then moved to loci that become records (below) |
| 4 | `31669bfd` | read alleles mapped through candidate selection and the drop of uncalled alternatives, carried on `VcfRecord`; an STR record that absorbs a generic one drops its read alleles and counts them |
| 5 | `26b33065` | the fragment keeper on the ordered thread, in both modes |
| 6 | `229628dd` | fragments beside the line to the writer; a skippable section of each spill entry, read only by pass three |
| 7 | `245a4e87`, `4c257f98` | a link's log-odds: each read allele's error floored at 0.01 and capped at ½, HapCUT2's fragment weight, the thresholds |
| 8 | `2ffaf7d5` | phase from reads alone, a stage in the writer; header lines; byte digests re-recorded |
| 9 | `6a41c7d2` | reads in the neighbour-evidence pass's link test, attachment, the read check |
| 10 | `a5941dc1` | invariance tests with the read cap binding; the simulator's `snp_rate`; scorers |

Two decisions inside the build's latitude, stated because they show in behaviour:

- **The keeper forgets a read** once its first record lies more than 1,000 bases behind the oldest
  record it holds, when none of its read alleles can still be kept. A read met again after that is
  numbered as a new fragment, so a pair with an insert over a kilobase can become two fragments. It
  can lose evidence and never join two reads.
- **The `GQ` the pass's link test reads** is computed from the posterior over every genotype, before
  the writer renormalises over the genotypes of the alleles some sample keeps. Renormalising can only
  raise it, so every link the test allows has written `GQ` 30 or more at both ends; a rare link the
  renormalised value would allow is missed.

## Where the capture runs: the measurement that decided it

Spec §2.2 asked for the capture to be measured at 200 samples before committing it to every cohort
locus, and moved to the loci that become records if it cost more than 5% of the merge's time.

On a 200-sample simulated cohort (64 founders, inbreeding 0.8, 160 samples at 3× and 40 at 30×,
2 Mb, 6 threads), capturing at every built locus cost 2.2 s of CPU: 69.5 s without and 71.7 s with,
means of three alternated runs. A profile put the merge — reading psps, closing and building loci —
at 47% of the run's busy CPU, about 34 s, so the capture was **6.4% of the merge's time**, above the
bar. It was moved: the merge hands each locus to calling with a handle on the samples' minted
records (`LocusReads`), which both drivers still hold while the locus is called, and calling asks it
only for a locus it writes. In the final build the capture is about 0.7% of the busy CPU of the same
run, roughly 3% of the merge's time.

## What it costs

Measured against the build of step 2 (psp format 1.1, no read phase), refinement running, 6 threads,
means of three alternated runs:

| simulated cohort | CPU, before → after | wall clock | peak memory | spill |
|---|---|---|---|---|
| 200 samples (40 at 30×), 2 Mb | 134.0 → 135.8 s (+1.3%) | 28.1 → 28.3 s | 411–438 MiB both ways | 38.1 → 41.8 MiB (+10%) |
| 50 samples (10 at 30×), 2 Mb | 25.2 → 25.2 s | 5.1 → 5.1 s | 258–291 MiB both ways | 8.4 → 9.0 MiB (+7%) |

In a profile of the 200-sample run, read phase is under 1% of the busy CPU: the capture 0.7%, the
keeper 0.04%, the spill's fragment section and the link sums too small to register. Peak memory did
not move beyond run-to-run noise.

**These cohorts are sparse, which is why the spill grows by a tenth and not threefold.** At one
record every 600 bases, most reads cover one record and are dropped by the keeper. At 3,000 samples
records are about 35 bases apart and nearly every read reaches a second record, so spec §4.1's
arithmetic applies: fragments roughly triple the spill at 30×. That extrapolation is the spec's, not
a measurement.

**What a posterior floor would drop** (the owner's deferred decision 7, spec §9 Q1): the share of
kept read alleles at sample columns where calling's probability of a heterozygote is below 10⁻³.

| simulated cohort | kept read alleles a 10⁻³ floor would drop |
|---|---|
| 200 samples, inbreeding 0.8, 40 at 30× | 581,253 of 716,256 (81%) |
| 50 samples, inbreeding 0.8, 10 at 30× | 85,839 of 111,095 (77%) |
| 50 samples, outbred, all at 3× | 2,904 of 57,534 (5%) |

So the floor saves most where samples are deep and confidently homozygous, and little at 3×, as the
spec predicted.

## How good the phase is

All scores count SNP records only, against the simulator's truth or, for HG002, the phase its
parents' genotypes give. A *link* is two consecutive phased records of a sample in one phase set;
it is *wrong* when the two records' written orientations disagree with the truth's. A row with fewer
than 300 links cannot show a rate under 1 in 100 and is not a decision.

### The nine cohorts of the phase output report, refinement on

Each cohort was called twice with refinement on: by this build, and by the build of step 8, whose
pass uses no reads. A link of the first run is *population* when the second run phases the same two
records together the same way, and *made by reads* otherwise; a link of the second run that the
first breaks is *broken by reads* (`read_phase_score.py`). SNPs are 1 in 1,000 bases, the
simulator's default.

| simulated cohort | wrong links, population only → with reads | links made by reads, wrong | links broken by reads, of them wrong |
|---|---|---|---|
| 50 (10 at 30×), inbreeding 0.5, 16 founders | 4 of 4,375 → 4 of 4,390 | 0 of 16 | — |
| same, inbreeding 0.8 | 2 of 2,174 → 2 of 2,179 | 0 of 9 | 0 of 2 |
| same, inbreeding 0.8, 64 founders | 6 of 1,184 → 6 of 1,187 | 0 of 4 | 0 of 1 |
| 200 (40 at 30×), inbreeding 0.8, 64 founders | 5 of 9,311 → 5 of 9,346 | 0 of 46 | 0 of 4 |
| 50 at 3×, inbreeding 0.5, 16 founders | 6 of 3,625 → 5 of 3,637 | 0 of 13 | 1 of 1 |
| 10 at 3×, inbreeding 0.95, 64 founders | none phased either way | — | — |
| 50 at 3×, outbred, 16 founders | 10 of 2,107 → 7 of 2,146 | 0 of 42 | 3 of 3 |
| 50 at 3×, outbred, 64 founders | 56 of 1,299 → 52 of 1,335 | 1 of 45 | 5 of 7 |
| 10 at 3×, outbred, 64 founders | 6 of 14 → 6 of 15 | 0 of 1 | — |

The population-only run reproduces the phase output report on the regenerated psps (for example 4 of
4,375 and 5 of 9,311 exactly). **Reads never made the population phase worse**: wrong links fell in
the three outbred-or-partly-inbred 3× cohorts and stayed within one elsewhere. Links made by reads
were wrong 1 time in 176 over all nine; no single cohort has 300 of them. The links reads broke were
mostly wrong ones in outbred cohorts (8 of 10) and right ones in inbred or partly inbred cohorts (1 of 8 wrong), where
breaking costs only block length. Each cohort phases 0.3% to 3% more true heterozygotes with reads.

### Outbred cohorts with dense SNPs, refinement on

The outbred cohorts again with SNPs 1 in 100 bases apart (the simulator's new `snp_rate`), 50
samples at 3×, 2 Mb, scored the same way:

| simulated cohort | wrong links, population only → with reads | links made by reads, wrong | links broken by reads, of them wrong | read check: contradicted of reached |
|---|---|---|---|---|
| 64 founders | 1,466 of 28,007 (5.2%) → 1,082 of 30,321 (3.6%) | 15 of 3,029 (0.5%) | 385 of 442 | 120 of 3,355 (3.6%) |
| 16 founders | 69 of 41,978 (0.16%) → 47 of 44,407 (0.11%) | 1 of 2,542 (0.04%) | 23 of 48 | 9 of 4,886 (0.18%) |

**This is the case the build was for, and both goals hold**: links made by reads are wrong under 1
time in 100, with enough links to show it (3,029 and 2,542), and the population links are better
than without reads, not merely no worse. At 64 founders the wrong population links fall by a quarter,
mostly because 87% of the links reads broke were wrong, but 3.6% is still above the 1-in-100 ceiling:
reads cannot reach most population links, as spec §2.7 foresaw (open question Q3). **The read check
tracks the true error**: 3.6% estimated against 3.9% among the population links the run kept
(1,067 of 27,292), and 0.18% against 0.11% to 0.16% at 16 founders. True heterozygotes phased rise
from 25% to 27% and from 38% to 41%. With refinement off, reads alone phase 11,851 and 8,696
heterozygotes, with 0 wrong of 6,413 and of 4,668 links.

### Phase from reads alone

The same nine cohorts with refinement off: 0 wrong of 575 links in all, every cohort under 300 links,
45 to 451 heterozygotes phased a cohort — SNPs 1 in 1,000 bases apart are mostly beyond a read
pair's reach.

Single samples, a 64-founder outbred genome with SNPs 1 in 100 bases apart, 2 Mb:

| depth | true heterozygotes phased | links | wrong |
|---|---|---|---|
| 3× | 67 of 1,165 (6%) | 35 | 0 |
| 10× | 1,287 of 2,518 (51%) | 746 | 0 |
| 30× | 1,846 of 2,680 (69%) | 1,147 | 0 |

The simulator misreads bases exactly as their quality says and has no paralogs, duplicates or
chimeras, so these measure reach and how links combine, not whether the error model holds.

GIAB HG002 in the trio benchmark's 100 random regions (452 kb, about 300 reads a position,
downsampled), called alone, the truth phase from the parents' genotypes
(`giab_trio_phase_score.py`; the GIAB v4.2.1 truth carries no phase in these regions):

| depth | heterozygous SNPs phased | links the trio scores | wrong |
|---|---|---|---|
| 300× | 196 of 289 (68%) | 107 | 0 |
| 30× | 197 of 289 (68%) | 106 | 0 |
| 10× | 106 of 289 (37%) | 52 | 0 |
| 3× | 20 of 289 (7%) | 10 | 0 |

No wrong link on real data at any depth, but no depth reaches 300 scored links, and the depths share
regions, so they cannot be pooled. 300× phases what 30× does because a sample keeps 32 read alleles a
record.

### A consequence of the ruled design, measured

**With refinement on, reads link records only next to ones the population orients.** In a
two-sample outbred cohort switched on with `--refine-genotypes on`, the consensus leans by 0.24 at
the first sample's two heterozygotes, under the 0.5 a candidate needs, so the dozen fragments that
link them phase nothing; the same cohort with refinement off writes them `0|1` in one phase set. On
the nine cohorts, reads-only phase reaches records refinement-on runs leave unphased; how many, per
cohort, is not measured here.

## Invariance

Tested end to end (`src/cli/mode_equivalence.rs`): on a fixture whose first sample carries two
heterozygous sites 40 bases apart in cis, both are written `0|1` in one phase set, and the VCF is the
same byte for byte at one thread and four, from alignments and from psps, and with the
hidden-duplication filter on (its own lines removed). On a deep version, 50 reads a position so the
read cap binds, the two routes write one VCF, and a `--regions` run over the walked region after the
repeat tract — whose reads a whole-contig run numbers from a non-zero first id — writes that region's
records exactly as the whole-contig run does. The three byte digests (`cross_platform_digests.rs`)
were re-recorded at step 8, each shown to reproduce the old checksum exactly once the new header
lines and the `PS` field are taken back out; they were recorded on macOS arm64 and still need a Linux
x86-64 run.

## Not done

- **The error floor's fit on HG002** (spec §7.2): it needs, per pair of heterozygotes, the fragments
  that cover both and which agree with the majority, which no file holds; it needs an instrumented
  run. `e_floor` stays at its starting value, 0.01.
- **The tomato read check** (§7.7): the cohort's psps would have to be regenerated in format 1.1
  from its CRAMs.
- **The share of reachable heterozygote pairs lost to STR segments** (§2.1, §7.7): unmeasured.
- **Correlated errors in segmental duplications** on HG002 with the filter on and off (§7.7).
- **3,000 samples**: the costs above are at 50 and 200; the larger figures are the spec's arithmetic.
