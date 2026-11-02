# Leaving duplicated tracts out of the stutter fit — built and measured, 2026-11-02

*Branch `stutter-fit-outliers`, from `str-inside-errors` (main e285095f plus the repaired stuttering
simulation). Design: [`../ng/spec/parameter_prepass_ssr.md`](../ng/spec/parameter_prepass_ssr.md)
§4.6.1, approved by the owner 2026-11-02. Nothing is committed.*

## 1. The answer

**The two rules fix what they were built for on HG002, and change nothing on the simulation. On
the tomato cohort the stutter-tail rule as specified removes real polymorphism.** A one-line
change to it — take each sample's genotype from its own reads alone — removes the problem and leaves
every HG002 and simulation number unchanged. **The owner chose it, and it is what the branch now
does** (§7); the specified genotype is gone from the code.

| | HG002 chr20 30× | HG002 chr20 5× | simulation, 1 sample 30× | simulation, 50 inbred samples 3× | tomato2, 59 accessions ~3× |
|---|---|---|---|---|---|
| tracts left out by depth | 278 | 218 | 5 of 3,000 | 0 | 240 of 15,778 |
| left out by the stutter tail, genotype as specified | 95 | 7 | 3 | 1 | **438** |
| left out by the stutter tail, genotype from the reads alone | 85 | 6 | 3 | 0 | 3 |
| 5-base motif slips per 100, 8 / 9 / 14 copies (fitted on every tract: 23 / 23 / 51 at 30×) | 0.8 / 1.0 / 1.7 | 2.9 / 3.6 / 3.4 | — | — | — |
| calls, missed / false (fitted on every tract) | 841 / 373 (842 / 374 prototype) | 24,913 / 4,346 (24,916 / 4,346) | 23 errors (23) | 6,233 errors (6,233) | not scored |

Both genotype choices give identical calls on HG002 at both depths and on both simulations.

## 2. What was built

**Two rules in the repeat-tract fit, in its own code** (`src/parameter_estimation/joint/ssr_fit/`):

- **The depth rule** (`tract_exclusion.rs`, `leave_out_deep_tracts`) runs once before the fit. It
  leaves out a tract whose crossing reads, summed over samples, exceed the larger of twice what one
  copy gives and the count one copy reaches by chance only 1 time in 1,000. *What one copy gives* is
  the run's typical depth times the share of reaching reads that cross a tract of that class — the
  correction for long tracts being crossed by fewer reads (spec §4.6.1, rule 1).
- **The stutter-tail rule** (`tracts_in_the_stutter_tail`, and the loop
  `fit_one_period_leaving_out_the_stutter_tail` in `penalised_fit.rs`) runs between fits of each motif
  length. Per tract, summed over samples: the reads off both alleles of the genotype each sample's
  reads favour, with no prior (§7). A tract is left out when that count is at least twice what the class's slip rate predicts
  and improbable at that rate (P < 0.001). Fit, test, remove, refit, until a round removes nothing or
  10 rounds have run. Standard errors are computed once, on the last fit.

**What the run reports.** The log, at the default level: the typical depth and the depth rule's
count, one line per motif length with the tail rule's count and rounds, and warnings for a class the
rules emptied or cut by more than 1 tract in 10. The parameters file: a new informational section,
`[repeat_tracts_left_out]` — the settings, per motif length the rounds and whether they settled, and
per class the tracts fitted and left out by the guard, by depth and by the tail. It is format 7 of
the file; a calling run reads none of it and carries it into the file it writes back.

**Tests.** Nine in `tract_exclusion.rs`. Each piece was checked by removing it and seeing its test
fail: the crossing-share correction, the chance term, the factor of two, the genotype (a heterozygote
is not stutter; with the class-prior genotype put back, the test of a tract polymorphic between inbred
samples fails), the sum over samples; the rounds' test fits a drawn class whose second tier of bad tracts only shows once
the first is gone (true rate 1.0 slips in 100: every tract 5.0, one round 2.3, settled 0.9). The
whole library suite passes in release (5,043 of 5,045); the two tests that check a debug-only panic pass in a
debug build. Three recorded checksums of a fitted cohort moved; in each, taking the new section and
the version line out of the new file reproduces the old checksum exactly, and when the section's note
was reworded for the genotype change, putting the old sentence back reproduced the intermediate ones;
so they were re-recorded
(`src/cli/cross_platform_digests.rs`).

## 3. HG002 chr20

**30×** (vcfeval against GIAB v4.2.1; calls on all of chr20 from the full psps). Slips per 100 reads:

| class | every tract | satellite removed by hand | depth, prototype | both rules, prototype | both rules, built |
|---|---:|---:|---:|---:|---:|
| homopolymer, 10 / 16 / 24 copies | 1.3 / 4.6 / 9.1 | 1.2 / 4.5 / 9.2 | 1.2 / 4.5 / 9.2 | 1.2 / 4.4 / 9.0 | 1.3 / 4.5 / 9.3 |
| 4-base motif, 22 / 24 copies | 4.6 / 38.9 | 6.8 / 32.5 | 6.6 / 32.0 | 4.2 / 10.2 | 7.8 / 27.7 |
| 5-base motif, 8 / 14 copies | 23.3 / 50.9 | 0.7 / 3.3 | 1.1 / 3.1 | 1.0 / 1.7 | 0.8 / 1.7 |
| 6-base motif, 8 / 12 copies | 3.8 / 11.0 | 3.9 / 11.6 | 3.8 / 11.0 | 4.1 / 10.6 | 4.0 / 6.2 |
| errors inside STR segments, missed / false | 118 / 123 | 110 / 119 | 109 / 119 | 108 / 118 | 108 / 118 |
| all of chr20, missed / false | | | 842 / 374 | | 841 / 373 |

The homopolymers stay inside the band measured at confident homozygous tracts (2.0, 4.6 and 9.4 in
100 below 12, at 12–17 and at 18 or more copies).

**The 4-base motif at 24 copies stays at 28 slips in 100**, where tracts known to be homozygous show
none. The prototype brought it to 10 by removing one tract at 23 copies with 9 of its 18 reads
off-length; through the curve, that pulled the neighbouring classes down. In the built version that
class's own rate is still inflated by the tract when it is tested, so the tract is not flagged. Classes
of 13+ copies hold 88 tracts on chr20; this is a thin-class effect on one chromosome, not chased here.

**5×.** Calls do not move: 24,913 missed and 4,346 false, against 24,916 and 4,346 fitted on every
tract; inside STR segments 3,145 and 804 against 3,148 and 804. The 5-base rates come down from 25,
36 and 54 to 2.9, 3.6 and 3.4 in 100 (8, 9, 14 copies), as the prototype's three-times-the-median
cut did.

## 4. The stuttering simulation

No duplications, so the rules should do nothing. They leave out 8 of 3,000 tracts at one sample 30×,
and 1 tract in the 50-sample cohort at 3×. Errors do not move (23 of 3,000 cells; 6,233 against
6,232 of 150,000). The fitted rates against the true ones (class rate times the libraries' mean
multiplier): median ratio 1.011 before and 1.008 after at 30×, 0.999 both in the cohort, and the
mean distance from the truth unchanged. A few classes move by up to 7%, in both directions; in the
cohort the one tract left out had been holding homopolymer rates up (11 copies: 2.98 → 2.77, true
2.43).

## 5. Tomato: the stutter-tail rule removes polymorphism

tomato2, 59 accessions at about 3 reads a position each, over 160 regions of 200 kb (catalog and
psps rebuilt with this branch). There is no truth, so the question is what the removed tracts look
like. *Polymorphic* below means the accessions' own most common lengths at a tract are not all the
same.

| tracts | homopolymers polymorphic | dinucleotides polymorphic | median reads, summed |
|---|---:|---:|---:|
| kept | 41 in 100 (of 12,913) | 59 in 100 (of 1,413) | 198 / 122 |
| left out by the tail rule, genotype as specified | **86 in 100 (of 324)** | **90 in 100 (of 102)** | 193 / 146 |
| left out by depth | 32 in 100 (of 129) | 56 in 100 (of 77) | 416 / 328 |

The tail rule's tracts are no deeper than the kept ones and twice as often polymorphic: they look
like tracts that differ between accessions, not like duplications. *The polymorphic share is
inflated in every row by stutter itself* — at a few reads, one slipped read can be an accession's
most common length — which is why it is 41 to 59 in 100 among kept tracts; the comparison between
rows is what carries the finding.

**The cause.** The genotype the tail rule counts against is the one the fit weighs highest, whose
prior is the class's allele shares, not the tract's own allele frequencies. At 3 reads, an
accession carrying a length rare in its class is given the class's common genotype, and its reads
count as stutter; summed over 59 accessions, a tract polymorphic between them looks improbable. In
a single sample at 30× (HG002) the reads outweigh the prior and the problem does not arise.

**What it does to the rates.** Dinucleotides of 9 to 12 copies fall 25 to 33% (9 copies: 7.31 →
5.42 slips in 100). With the genotype from the reads alone those classes stay where they were or
rise (7.81).

**The binomial check fails on tomato too.** At tracts where every accession agrees on the length and
at least 10 reads summed, off-length reads vary from tract to tract 2.1 times as much as a binomial
predicts at homopolymers (6,236 tracts), 2.7 times at dinucleotides (444) and 3.5 times at 3- to
6-base motifs (185); 64 of 6,865 lie beyond P = 0.001 where about 7 would by chance. HG002's figure
was 0.65 to 1.15. *Not separated here*: the tomato libraries' own stutter differs about fivefold, and
the fit pools them, so tracts read by different mixes of accessions can differ for that reason alone.
With the genotype from the reads alone the tail rule removes 3 tomato tracts, so this matters little
in practice; it means the 0.001 is anti-conservative on this cohort.

## 6. Tomato: what the depth rule does

240 tracts (1.5 in 100), the same with either genotype. Their summed depth is twice the kept tracts'.
It moves some rates down — a 4-base motif of 8 copies from 4.72 to 0.62 slips in 100 — and some up:
dinucleotides of 14, 16 and 20 copies rise from 17.0, 21.7 and 21.6 to 20.8, 33.6 and 32.1. **The
rises are not explained.** A collapsed duplication of two equal-length copies adds reads without
adding off-length ones, so removing it should not raise a rate; that the removed deep tracts slipped
less than their class is possible but unchecked. Tomato has no truth to settle it.

## 7. The decision: each sample's genotype from its own reads (owner, 2026-11-02)

**The tail rule counts off-length reads against the genotype each sample's reads favour, with no
prior.** The specified genotype — the fit's own weights, with a prior from the class's allele shares —
was removed from the code; the spec records the change (§4.6.1, rule 2, step 1). What it gains and
costs, measured:

- **Gain:** on tomato it removes 3 tracts where the specified genotype removed 438, almost all of
  them polymorphic, and the dinucleotide rates no longer drop by a quarter.
- **No cost measured:** identical calls on HG002 at 30× and 5× and on both simulations. HG002's rates
  move by at most 0.1 slip in 100 (5-base motif 0.9 / 1.1 / 1.7 against 0.8 / 1.0 / 1.7), except the
  4-base motif at 24 copies, 28.7 against 27.7.
- **Cost by construction:** with no prior, two lengths in one sample are always two alleles. A
  duplication whose copies put two lengths in every sample then looks heterozygous, and at a few reads
  a sample the tail rule cannot see it; only a third length shows. At those depths the depth rule
  is what catches duplications. The test of summing over samples now uses a third length.

**Considered and not built:** a prior from **the tract's own** allele frequencies, estimated across
samples, which explains polymorphism between inbred accessions without explaining a duplication's
extra lengths inside each one. It would add one case — a duplication of two lengths in an inbred
cohort at low depth — that the depth rule summed over samples already covers, at the cost of a
per-tract estimate made from the same reads, stutter included. The owner chose not to build it.

## 8. Open, and where it goes next

- **The 4-base motif at 24 copies on HG002** (§3): thin classes on one chromosome.
- **The depth rule's rises on tomato** (§6): unexplained.
- **Tomato's overdispersion** (§5) and **per-library stutter**: the owner's next item after this
  branch is fitting stutter per library (`str_library_stutter.md`); it would also say how much of
  tomato's 2–3.5-fold spread is library mixing.
- **The cheap test the owner approved**: call HG002 with the 3- to 6-base rates set near their true
  values, to size what the remaining overestimate costs.

## 9. Data and scripts

All under the worktree's `tmp/`: `run_hg002.sh` (refit, call and score HG002 chr20 from the
`str-inside-errors` worktree's psps), `run_measure.sh` and `run_sim_rest.sh` (simulation reruns,
`REDO="ng_fit ng_fitted score"`), `run_tomato.sh` (catalog, psps, both fits), `run_reads_alone.sh`
(the variant everywhere), `by_distance.py`, `rates.py`, `rate_truth.py`, `rate3.py`,
`tomato_tracts.py`. The frozen binaries are `tmp/bin/` (as specified) and `tmp/bin_reads_alone/` (the latter is what the branch now builds).
The tomato tract tables come from `examples/ng_str_tract_spread.rs`, which loads a cohort's tract
evidence as `estimate-parameters` does, runs the rules, and labels every tract.
