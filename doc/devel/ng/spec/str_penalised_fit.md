# The repeat-tract fit by penalised EM

**Status:** approved by the owner, 2026-10-14. The method is chosen
(`doc/devel/ng/research/str_fit/recommendation_2026-10-14.md`, adopted by the owner the same day); this
settles what is built. A research prototype exists and is measured; nothing is wired into production.
Companion spec: [`str_fit_to_calling.md`](str_fit_to_calling.md), what the fit hands to calling. It
replaces the unbuilt draft `str_counting_fit.md`, whose one settled requirement (§2.1, standard errors
that allow for samples sharing a locus) is §2.6 here. Architecture doc and implementation plan follow
on approval.

---

## 1. What this is

### 1.1 What we are building

**`estimate-parameters` fits, for each motif length, how repeat tracts behave along the reference
repeat count.** A *class* is every locus with one motif length (period) and one reference repeat
count, for example every dinucleotide locus with 12 repeats in the reference. For each motif length
the fit estimates, as a curve along the reference repeat count:

- **stutter**, the measurement error: the slip rate (how often a read reports a length other than its
  allele's), the shorter share (of slipped reads, how many came out shorter) and the fall-off (how
  fast slips of two, three or more repeats become rarer than slips of one);
- **the allele length variation:** over every chromosome of the cohort at every locus of the class,
  the fraction carrying each length, from 9 repeats shorter than the reference to 9 longer. The two
  outermost lengths mean "9 or further" (`str_allele_window.md` §1.2);
- **the between-locus share F:** how much of that variation lies between loci rather than among
  samples at one locus. It is the Balding–Nichols F_ST with the class's loci in place of populations:
  at 0 every locus carries the class's whole variation, at 1 each locus is fixed at one length.

Each sample's homozygote excess is an input, from the SNP/indel half of the same fit.

**How it fits, in one paragraph.** The fit repeats three steps until no number moves. (1) For each
sample at each locus, weigh every genotype it could carry by how well it explains the sample's reads
under the current curves. (2) Count, with those weights, at each reference repeat count: reads read
exactly and reads that slipped, which way and how far, chromosomes at each length, and how often two
samples at one locus carry the same length. (3) Draw each curve through those counts, close to a
repeat count backed by many reads and smoothly past one backed by few, charging it for bending; how
stiff each curve is, is estimated from the data. The method is penalised expectation maximisation
(EM): step 1 is the expectation, steps 2–3 the maximisation.

**It replaces the fit in main**, which fits each class on its own by a climb (an optimiser taking
trial steps), averaging each locus over its possible allele frequencies at 256 points, and then draws
curves through the per-class answers and blends each class with its curve
(`parameter_estimation/joint/ssr_fit.rs`; `spec/str_slippage_level_curve.md`).

### 1.2 Why — the measurements that decided it

On simulated cohorts with known truth, 1 to 2,000 samples, 3 to 300 reads a sample
(`research/str_fit/compare_D_prior_2026-10-14.md`, `cost_2026-10-14.md`):

- **It predicts loci it was not fitted on better than every alternative, in all 45 cells tested.**
  Against the counting fit (EM per class, curves afterwards): +0.10 to +0.15 log units a held-out
  locus at one sample, +0.01 to +0.05 at 500–2,000 samples. The fit in main predicted worse than both
  where it could be run.
- **It is thousands of times faster.** tomato2 (59 accessions, 32 Mb): 5.8 s at ±9 on the Mac, against
  6.7 hours for main's fit at ±6. Estimated at kimura's scale (2,169 samples, 66 million (locus,
  sample) pairs with reads): 7–23 minutes at 3 reads a sample, up to 28 at 30 reads.

### 1.3 Goals

1. Every number of §1.1 for every class with reads, at every cohort size from 1 sample to thousands
   and every depth from 3 reads to hundreds.
2. A standard error for every number it reports, honest where the data can support one (§2.6).
3. **A number the data cannot place is absent, not a fitted bound**: F at one sample (§2.5).
4. At kimura's scale: within one hour on one 64 GiB machine, the fit alone.
5. The same output bit for bit on every run and thread count (§3, first trap).

### 1.4 Non-goals, and what this does not do

- **One stutter for all read groups.** It does not fit stutter per slippage group or library
  multipliers. Production declares one slippage group today (`every_read_group_pooled`,
  `run/census_fit.rs:282-288`), so nothing is lost now; the climb's support for several groups goes
  with it. **Deferred to the library-multiplier work** (`spec/str_library_stutter.md`), which must
  bring them back on top of this fit.
- **No parametric allele-length shape.** Each of the 19 lengths keeps its own curve along the repeat
  count. A few-number shape across lengths (research avenue C) was measured and is not needed for this
  decision (`research/str_fit/trackC_tomato2_span9_2026-10-12.md`, `kimura_shapes_span9_2026-10-14.md`).
- **No sample subsets.** Main's fit grows a subset of samples per class to bound its cost
  (`SampleSubsets`, `ssr_fit.rs:1149-1157`); this fit uses every sample, within budget.
- **It does not change** the census, how loci are selected or guarded, the substitution rate inside
  tracts, or the SNP/indel half of the fit.
- **It does not correct** the bias against long alleles from reads that must cross the whole tract
  (`str_allele_window.md` §5.3), nor widen the ±8 recording.

---

## 2. The design

### 2.1 Where it sits

`fit_a_cohort` (`run/census_fit.rs:206`) calls today `ssr_fit::fit_strata_on_sample_subsets`
(`:262`) with every stratum's evidence, then `draw_and_apply_curves` runs inside it
(`ssr_fit.rs:3127-3130`). **The new fit takes that call's place**: same evidence
(`the_tract_strata_of_a_cohort`, `census_fit.rs:180-196`), same homozygote excess vector, one call per
motif length, its output turned into the `StratumOutcome`s that `StratumFits::over`
(`stratum_fits.rs:504-607`) already reads. Nothing downstream of `StratumFits` changes except what
[`str_fit_to_calling.md`](str_fit_to_calling.md) settles.

The prototype is `fit_period_by_penalised_counting` (`ssr_fit/penalised_fit.rs:328`), with its
E-step from the counting fit (`composite_fit.rs`: `expected_counts`, `SamplePairs`) and its errors
(`counting_fit_errors.rs`, `penalised_fit.rs`). It is built to be promoted, not rewritten.

### 2.2 One loop per motif length; classes are fitted together

**All classes of one motif length are fitted at once**, because every curve runs across them. Motif
lengths are independent and may run in parallel. Within one, each round's step 1 runs over classes
and their loci in parallel; steps 2–3 are cheap (1–9% of the time on the bench).

**Every class with reads at even one locus is fitted.** There is no refusal floor: a thin class's
numbers come mostly from the curve, which is the point. A class with no reads at any locus has no
numbers of its own; calling reads its values off the curves inside their fitted range and holds them
flat beyond, as today (`str_slippage_level_curve.md` §6). *Main's floor of 8 loci
(`DEFAULT_REFUSAL_FLOOR`, `ssr_fit.rs:1199`) existed to stop a per-class climb on too little; it has
no job here.*

### 2.3 The curves and how stiff they are

- **One value per reference repeat count**, from the lowest class with reads to the highest, with a
  penalty on the bend between consecutive values (a Whittaker smoother; the limiting case of
  P-splines with a knot at every repeat count).
- **Scales:** the three stutter numbers and F on the logit scale; each allele length as the log of its
  share over the reference length's share, so shares stay positive and sum to one. The two edge
  lengths get curves like any other and are never forced onto a falling tail.
- **Stiffness is estimated, separately for each of the five kinds of number** (one strength shared by
  all 18 allele-length curves): the strength that maximises the Laplace approximation to the marginal
  likelihood of that number's counts (as REML does for generalised additive models), chosen afresh
  while the fit moves and then held, so the loop converges (`penalised_fit.rs:96-118`).
- **A real kink survives.** A bend costs its square up to 2.5 standard deviations and only
  `2·ln(bend/2.5)` beyond, so a class whose many reads insist on a step keeps it (measured on a drawn
  kink: `research/str_fit/track2_penalised_em_2026-10-12.md` §3.3). The 2.5 is the existing blend's
  knee (`str_slippage_level_curve.md` §7.2), inherited; the 2 was measured against 6.25.

### 2.4 The floor under the allele lengths — decided

**Every allele length of every class gets a prior weight of one locus's worth of allele draws,
spread evenly over the 19 lengths**, on top of one draw a class (`allele_prior_tracts = 1`,
`allele_pseudo_draws = 1/19`, `penalised_fit.rs:57-75`). Without a floor, a length that reaches zero
can never return and a held-out locus carrying it costs hundreds of log units.

**The cost, decided by the owner's lead on 2026-10-14:** at 3 reads a sample in cohorts of 500 or
more, the slip rate comes out 0.3–0.7% low (largest in long, thin classes), and its errors do not
cover that. Eight reshaped floors were measured; each removed the bias only as fast as it lost
prediction (`research/str_fit/shares_prior_2026-10-14.md`). Whether 0.5% changes any call is not
measured (§5).

### 2.5 F, and when it is absent

- **F is read from pairs of samples at a locus**, alternating with the loop at most 10 times, because
  in an inbred cohort a sample's own two copies say almost nothing about it
  (`pair_alternations = 10`).
- **With no locus holding two samples with reads, F is absent** (`None`), and so is its error. The
  loop still needs a working value and uses one read off each sample's own copies
  (`working_concentration`, `penalised_fit.rs:214`); **it is never written out as an estimate**. What
  calling does with an absent F is [`str_fit_to_calling.md`](str_fit_to_calling.md)'s.

### 2.6 Standard errors

**Every reported number carries one, computed as `H⁻¹(B + P)H⁻¹`**: H the penalised curvature of the
objective, B the sum over loci of each locus's score contribution times itself, P the bend penalties.
This is the Bayesian covariance of a penalised fit with B in place of the curvature, **because the
objective treats samples at one locus as independent and they are not**; errors that ignore this held
the truth within one error only 0–2 times in 100 at 500 samples, where an honest error holds it about
68 times (`track1_counting_fit_errors_2026-10-12.md`).
Measured honest for stutter, F and allele lengths of 5% or more. **Not honest** for lengths under 1%
from 60 samples up (the truth inside one error 27–42 times in 100) nor for the homopolymer shorter
share at 3 reads (29–55); rare lengths get no error (§7).

### 2.7 The settings, as named constants

What the bench measured as the lead candidate (`bench_fixtures/fit_bench.rs:497-515`,
`examples/ng_str_fit_bench.rs:186-190`). The prototype's own defaults differ; **production uses
these**:

| setting | value | source |
|---|---|---|
| allele window | ±9 | `str_allele_window.md` |
| start | slip rate 0.05, shorter share 0.5, fall-off 0.3 | bench; soft |
| rounds at most / stop when no number moves more than | 2,000 / 10⁻⁷ | bench; soft |
| floor | one locus's worth + one draw a class | §2.4, measured |
| F from pairs, at most | 10 alternations | bench; soft |
| kink knee / tail | 2.5 / 2 | §2.3 |
| stiffness held once it moves less than | 0.1 natural-log units | measured, track 2 §10.1 |
| a class whose numbers stopped moving reuses its counts | on | measured, every number unchanged to the third decimal |
| clustering weight on allele draws | off | measured harmful at large cohorts, track 2 §2.5 |

### 2.8 Two exact speed-ups, kept

- **Each distinct row of reads is scored once.** A sample's contribution depends only on its row of
  read counts and its homozygote excess; at 3 reads tomato2 has 86 rows for each distinct one.
- **Each pair of samples is four numbers**, the coefficients of a cubic in the concentration, in place
  of a 2.9 kB table.

Both change no number beyond about 10⁻⁶ (`cost_2026-10-14.md`). **Their gain falls with depth**: at
30 reads a row repeats 3.6 times, at 300 reads 1.2 times.

### 2.9 What the run logs

Per motif length: rounds, the five stiffness strengths and where bends were eased, the time in each
step, peak memory. Per class: rows with reads, distinct rows, and whether F was absent. These replace
the climb's summaries (`ssr_fit.rs:3163-3330`).

---

## 3. What will bite

- **Output must be bit-identical across runs and thread counts.** The oracle
  (`scripts/promote_ng_oracle.sh`) pins checksums. The prototype differs by up to 2×10⁻¹⁰ between two
  runs, most likely from sums taken in thread order inside the penalised fit
  (`cost_2026-10-14.md`, "uncertain"). Sum in a fixed order, as `ssr_fit.rs` already does elsewhere.
- **The branch is 56 commits behind main** (main at 0e905bbb on 2026-10-14) and carries 61 of its own:
  the ±9 window (not on main), half of the renaming, and the research code. Merge main first.
- **The renaming is half done** on this branch (`aa42d629`): the parameters file's keys were renamed
  without a format version or reading of older files, so this branch refuses every parameters file
  written before it ("unknown field `length_spectrum_by_period`" on kimura). Fixed by
  `str_fit_to_calling.md`; build the two together.
- **Calling panics on a non-positive concentration** (`stratum_fits.rs:1051-1055`) and stores a class's
  shares and concentration as one pair (`:379-385`). A class with shares but an absent F has nowhere
  to go until the companion spec's change lands.
- **Memory at deep reads.** The prototype holds every class's distinct rows at once during a round:
  4.8 GB at 2,000 samples × 30 reads against 1.1 GB for the counting fit, which would not fit 64 GiB at
  kimura's scale. **Hold one class's rows at a time** (or bound the classes in flight).
- **The evidence is already the largest memory**: `gather_strata` holds every stratum at once
  (`ssr_fit.rs:4861-4885`), 19.4 GiB on kimura for 66.2 million pairs with reads, about 315 bytes a
  pair against tomato2's 139. Not changed here (§6).
- **The provenance written beside each stutter number assumes main's curves** (`LevelSmoothing`,
  `parameters_file/mod.rs:1190`, carries a fitted `SlippageCurve`; `LevelSource` is cell / curve /
  blend, `slippage_curve.rs:519-531`). A penalised curve is neither; the companion spec settles what
  is written.
- **`--str-param-estimates-at-once`** (`cli/estimate_parameters.rs:144-145`) sizes the climb's
  likelihood tables; it means nothing to this fit. Decided in §7.
- **A stale comment** in `fit_strata_with` says thin strata get "two shares copied from a neighbour"
  (`ssr_fit.rs:3008-3013`); nothing copies a neighbour any more (`:2955`). Do not build on it.

---

## 4. How we know it works

1. **Unit tests carried from the prototype:** with no penalty the fit equals the counting fit class by
   class; a drawn kink on many reads survives while a thin class borrows; the banded solve equals a
   dense one; the distinct-row and cubic-pair shortcuts equal the direct computation; F is absent at
   one sample.
2. **The bench reproduces the recommendation's numbers** for the production settings at the cells of
   `compare_D_prior_2026-10-14.md` (same seeds), within run-to-run noise.
3. **Determinism:** two runs, and runs at 1 and 18 threads, give the same bits.
4. **Real data, held-out:** on tomato2, held-out prediction by whole locus, the new fit against the
   counting fit, with the scorer of the bench (`fit_bench.rs`). This is the first test on real reads.
5. **Cost on kimura, measured:** the fit alone within one hour and the run's peak memory reported, at
   2,169 samples.
6. **The oracle moves, and the move is explained** before its checksums are re-recorded
   (`scripts/promote_ng_oracle.sh`; the fitted file and the calls made with it).
7. **Calling does not get worse** on HG002 against its repeat-tract truth and on the GIAB trio.

---

## 5. Open

1. **Does a slip rate 0.5% low change any call?** Measure at 500 and 2,000 simulated samples × 3 reads
   with calling. If it does, revisit the floor (§2.4). *Leaning:* it does not; settled by that
   measurement.

## 6. Deferred, with a recommended home

- **Stutter per slippage group and library multipliers** → `spec/str_library_stutter.md`, rebuilt on
  this fit.
- **Gathering evidence one motif length at a time**, to bound the evidence's memory → the census-read
  step (`gather_strata`), once kimura's 19.4 GiB matters; the reason for reading the whole band
  (neighbour borrowing, `census.rs:2194-2200`) no longer holds.
- **The read guard dropping 45,432 of kimura's 86,688 loci** (`kimura_shapes_span9_2026-10-14.md`) →
  its own question; not this fit's.
- **A parametric allele-length shape** → research avenue C, if thin classes need more than the
  curves give.

## 7. Decisions

- **2026-10-14, owner:** adopt penalised EM with the one-locus floor.
- **2026-10-14, owner:** the CLI flag `--str-param-estimates-at-once` is removed, with an error
  naming why, rather than accepted and ignored.
- **2026-10-14, owner: errors on rare allele lengths are reported as absent.** A share under 3 loci's
  worth (share × loci with reads < 3) gets no error, and the log says how many were withheld; such
  errors were measured too small (§2.6). The 3 rests on two settings measured by track 1; soft.
- **2026-10-14, owner: the overshoot at the ends of the repeat range at 1 sample × 3 reads is
  accepted** (+49% to +127% at the extreme class, `track2_penalised_em_2026-10-12.md` §10.3), and
  stated as a known limit of the fit.
- **2026-10-14, owner: the climb's code is deleted** once §4.4 and §4.7 pass, with its tests and
  examples, in one commit.
- **2026-10-14, owner:** stutter per slippage group and library multipliers go to the
  library-multiplier work (§6).

## 8. Reuse map

| what | existing code | reuse |
|---|---|---|
| the loop, curves, stiffness, errors | `ssr_fit/penalised_fit.rs` | promote; settings per §2.7 |
| step 1 (weigh genotypes) and the pair likelihood | `ssr_fit/composite_fit.rs` (`expected_counts`, `SamplePairs`) | promote |
| read probabilities per allele | `Slippage::read_probabilities`, `ssr_fit.rs` | as is |
| evidence | `the_tract_strata_of_a_cohort`, `gather_strata` | as is |
| into calling's tables | `StratumOutcome` → `StratumFits::over` | as is, plus the companion spec |
| the held-out scorer and the simulator | `bench_fixtures/{fit_bench,cohort_draw}.rs` | the acceptance tests of §4 |
| the parity oracle | the bench's `D_tract_prior` method on `str-fit-compare` | §4.2 |
