# The repeat-tract fit by penalised EM

**Status:** approved by the owner, 2026-10-14. The method is chosen
(`doc/devel/ng/research/str_fit/recommendation_2026-10-14.md`, adopted by the owner the same day); this
settles what is built. A research prototype exists and is measured; nothing is wired into production.
Companion spec: [`str_fit_to_calling.md`](str_fit_to_calling.md), what the fit hands to calling. It
replaces the unbuilt draft `str_counting_fit.md`, whose one settled requirement (§2.1, standard errors
that allow for samples sharing a locus) is §2.6 here. Architecture:
[`../arch/str_penalised_fit.md`](../arch/str_penalised_fit.md); plan:
[`../impl_plan/str_penalised_fit.md`](../impl_plan/str_penalised_fit.md).

**Amended 2026-10-14 after an adversarial review** of this spec, its companion, the architecture and
the plan. Corrections of fact are made in place. **Changes that are the lead's decisions, not the
owner's, are marked *(lead, 2026-10-14)* and listed in §7**; the owner confirmed the withheld errors of §2.6 on 2026-10-14.

**Amended 2026-10-15: a weak prior under the three stutter curves (§2.10), approved by the owner
2026-10-15.** The milestone C review found slip rates down to 4×10⁻¹⁸ that no read placed. The owner
chose the direction the same day: a weak prior on the stutter curves, like the floor §2.4 gives the
allele lengths. The lead proposed its form, weight and centre; the owner answered each open question
the same day (§2.10.9, §7). It is built in the plan's step C8. **The same day the owner also changed
which classes the stated F's median counts** (companion [`str_fit_to_calling.md`](str_fit_to_calling.md)
§3.5; plan step C9).

**Amended 2026-10-15, after the review of step C8, with four answers from the owner the same day:** a
slip rate with fewer expected slipped reads than the prior's weight gets no error (§2.10.4); at one
sample the prior's pull on two-base direction numbers is a known dishonesty of their errors (§2.6);
the prior's weight stays one read, on a three-seed sweep (§2.10.2); and the fit is made independent
of the order the psps are named in before milestone D (§3, plan step C10).

**Amended 2026-10-16, owner (§7):** where a thin class's errors start (§2.6) — the shorter share and
the fall-off from 3 loci with reads, the slip rate, F and the allele-length shares from 8 — and a rare
share's line lowered from 3 loci's worth to 1. Built in plan step D8.

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

**How it fits.** The fit repeats three steps until no number moves. (1) For each sample at each
locus, weigh every genotype it could carry by how well it explains the sample's reads under the
current curves. (2) Count, with those weights, at each reference repeat count: reads read exactly and
reads that slipped, which way and how far, chromosomes at each length, and how often two samples at
one locus carry the same length. (3) Draw each curve through those counts, close to a repeat count
backed by many reads and smoothly past one backed by few, charging it for bending; how stiff each
curve is, is estimated from the data. The method is penalised expectation maximisation (EM): step 1
is the expectation, steps 2–3 the maximisation.

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
- **Which cells those are.** No held-out score was run at 500 samples × 30 and 300 reads, nor at
  2,000 samples × 10, 30 and 300 reads, so large cohorts at deep reads are untested for accuracy; the
  2,000-sample cells of the PCR-free truth and of F = 0.4 rest on one seed each
  (`compare_D_prior_2026-10-14.md` §5). The plan's first check adds 500 × 30, 2,000 × 10 and three
  seeds at every 2,000-sample cell.
- **It is thousands of times faster.** tomato2 (59 accessions, 32 Mb): 5.8 s at ±9 on the Mac, against
  6.7 hours for main's fit at ±6. Estimated at kimura's scale (2,169 samples, 66 million (locus,
  sample) pairs with reads): 7–23 minutes at 3 reads a sample, up to 28 at 30 reads.

### 1.3 Goals

1. Every number of §1.1 for every class with reads, at every cohort size from 1 sample to thousands
   and every depth from 3 reads to hundreds. §1.2 says which of those cells were measured.
2. A standard error for every number it reports, honest where the data can support one (§2.6).
3. **A number the data cannot place is absent, not a fitted bound**: F at one sample, and F stopped at
   an end of its search (§2.5). For stutter the same goal is met by a prior and an error rule
   (§2.10, approved by the owner 2026-10-15): no stutter number is set by a numerical guard, and one that rests
   more on the prior than on the class's reads gets no error.
4. At kimura's scale: within one hour on one 64 GiB machine, the fit alone.
5. The same output bit for bit on every run and thread count (§3, first trap).

### 1.4 Non-goals, and what this does not do

- **One stutter for all read groups.** It does not fit stutter per slippage group or library
  multipliers. Production declares one slippage group today (`every_read_group_pooled`,
  `run/census_fit.rs:288`), so nothing is lost now; the climb's support for several groups goes
  with it. **Deferred to the library-multiplier work** (`spec/str_library_stutter.md`), which must
  bring them back on top of this fit. The multiplier prototype (`ssr_fit/library_stutter.rs`) uses
  the climb's quadrature (`:71`, `:107`) and is deleted with the climb *(lead, 2026-10-14)*.
- **No parametric allele-length shape.** Each of the 19 lengths keeps its own curve along the repeat
  count. A few-number shape across lengths (research avenue C) was measured and is not needed for this
  decision (`research/str_fit/trackC_tomato2_span9_2026-10-12.md`, `kimura_shapes_span9_2026-10-14.md`).
- **No sample subsets.** Main's fit grows a subset of samples per class to bound its cost
  (`SampleSubsets`, `ssr_fit.rs:1139`); this fit uses every sample, within budget.
- **It does not change** the census, how loci are selected or guarded, the substitution rate inside
  tracts, or the SNP/indel half of the fit.
- **It does not correct** the bias against long alleles from reads that must cross the whole tract
  (`str_allele_window.md` §5.3), nor widen the ±8 recording.

---

## 2. The design

### 2.1 Where it sits

`fit_a_cohort` (`run/census_fit.rs:220`) calls today `ssr_fit::fit_strata_on_sample_subsets`
(`:262`) with every stratum's evidence, then `draw_and_apply_curves` runs inside it
(`ssr_fit.rs:3124`). **The new fit takes that call's place**: same evidence
(`the_tract_strata_of_a_cohort`, `census_fit.rs:180-196`), same homozygote excess vector, one call per
motif length, its output turned into the `StratumOutcome`s that `StratumFits::over`
(`stratum_fits.rs:504-607`) already reads. Downstream of `StratumFits`, what changes is what
[`str_fit_to_calling.md`](str_fit_to_calling.md) settles, including the two places that decide
whether a stutter number's error is written (its §2.3).

The prototype is `fit_period_by_penalised_counting` (`ssr_fit/penalised_fit.rs:328`), with its
E-step from the counting fit (`composite_fit.rs`: `expected_counts`, `SamplePairs`) and its errors
(`counting_fit_errors.rs`, `penalised_fit.rs`). It is built to be promoted, not rewritten.

### 2.2 One loop per motif length; classes are fitted together

**All classes of one motif length are fitted in one loop**, because every curve runs across them:
each round updates every class's numbers before the curves are drawn again.

**Within that loop the classes are visited one at a time, and the parallelism is inside a class**,
over its distinct rows and loci; motif lengths also run one at a time *(lead, 2026-10-14)*. So at most
one class's table of distinct rows is alive at once, whatever the thread count (the memory trap of
§3). The cost: a class small enough to leave the machine idle is fitted on few threads, which at one
sample costs seconds. Steps 2–3 are cheap (1–9% of the time on the bench).

**Every class with reads at even one locus is fitted.** There is no refusal floor: a thin class's
numbers come mostly from the curve, which is the point. *Main's floor of 8 loci
(`DEFAULT_REFUSAL_FLOOR`, `ssr_fit.rs:1201`) existed to stop a per-class climb on too little; it has
no job here.*

**A class with no reads at any locus is refused, as today** (`Refused { NoSpanningReads }`,
`ssr_fit.rs:2691-2697`) *(lead, 2026-10-14: keep the climb's behaviour)*. It gets no values off the
curves. Calling scores a candidate in such a class under the shipped default stutter,
`StutterModel::hipstr_shipped()`, with provenance `Defaulted`
(`calling/inference/repeat_tract_parameters.rs:386-400`; test
`a_candidate_the_fit_never_reached_takes_the_stated_constant_and_says_so`, `:1774`). *(Corrected:
this section said calling reads such a class's values off the curves and holds them flat beyond; it
does not.)*

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
- **A known approximation in that choice:** the strength is chosen on step 1's expected counts as
  though they were observed counts (`CurveData::Binomial`, `penalised_fit.rs:1207-1210`), which
  ignores the information lost to not knowing the genotypes, largest at low depth. It is checked only
  indirectly, by held-out prediction on the bench.
- **A real kink survives.** A bend costs its square up to 2.5 standard deviations and only
  `2·ln(bend/2.5)` beyond, so a class whose many reads insist on a step keeps it (measured on a drawn
  kink: `research/str_fit/track2_penalised_em_2026-10-12.md` §3.3). The 2.5 is the existing blend's
  knee (`str_slippage_level_curve.md` §7.2), inherited; the 2 was measured against 6.25.
- **A motif length with one or two classes has no bend to penalise** (`bends()` is the classes minus
  two, saturating at zero, `penalised_fit.rs:1429-1431`), so its curves are the classes' own counts
  and the strength criterion is flat. §4.1 tests it. *(Owner, 2026-10-15:)* for the three stutter
  numbers, "own counts" means the counts plus the prior of §2.10; without it such a class with no
  slipped read gets the rate the 10⁻⁶ ridge sets (§2.10.1).

### 2.4 The floor under the allele lengths — decided

**Every allele length of every class gets a prior weight of one locus's worth of allele draws,
spread evenly over the 19 lengths**, on top of one draw a class (`allele_prior_tracts = 1`,
`allele_pseudo_draws = 1/19`, `penalised_fit.rs:57-75`). Without a floor, a length that reaches zero
can never return and a held-out locus carrying it costs hundreds of log units.

**The cost, decided by the owner's lead on 2026-10-14:** at 3 reads a sample in cohorts of 500 or
more, the slip rate comes out 0.3–0.7% of its value low (median over classes; largest in long, thin
classes), and its errors do not cover that. Eight reshaped floors were measured; each removed the
bias only as fast as it lost prediction (`research/str_fit/shares_prior_2026-10-14.md`). Whether a
slip rate 0.5% of its value low changes any call is not measured (§5).

### 2.5 F, and when it is absent

- **F is read from pairs of samples at a locus**, alternating with the loop at most 10 times, because
  in an inbred cohort a sample's own two copies say almost nothing about it
  (`pair_alternations = 10`).
- **With no locus holding two samples with reads, F is absent** (`None`), and so is its error. The
  loop still needs a working value and uses one read off each sample's own copies
  (`working_concentration`, `penalised_fit.rs:214`); **it is never written out as an estimate**. What
  calling does with an absent F is [`str_fit_to_calling.md`](str_fit_to_calling.md)'s.
- **F is also absent when its estimate stops within 1% of an end of its search** *(lead,
  2026-10-14)*. The prototype clamps the concentration to [10⁻³, 10⁴] (`penalised_fit.rs:379`,
  `:720-721`), which is F in [10⁻⁴, 0.999], and today reports a value at the clamp as fitted; a motif
  length with one shared locus, or a cohort whose loci are almost all fixed, can land there. The rule
  is the counting fit's (`concentration_is_placed`, `composite_fit.rs:83-98`;
  `track1_counting_fit_errors_2026-10-12.md` §7, "The rule now").

### 2.6 Standard errors

**Every reported number carries one, computed as `H⁻¹(B + P)H⁻¹`**: H the penalised curvature of the
objective, B the sum over loci of each locus's score contribution times itself, P the fit's prior —
the bend penalties, the stutter prior's curvature (§2.10.4) and the allele-length prior's (below).
This is the Bayesian covariance of a penalised fit with B in place of the curvature, **because the
objective treats samples at one locus as independent and they are not**; errors that ignore this held
the truth within one error only 0–2 times in 100 at 500 samples, where an honest error holds it about
68 times (`track1_counting_fit_errors_2026-10-12.md`).

**The allele-length prior's curvature is in H and in the middle term** *(owner, 2026-10-16; plan step
D7)*. The floor of §2.4 adds `q` draws to each of a class's lengths (both floors together), so the
objective the fit climbs carries `Σ q·ln s` over the lengths, `s` the shares. On the share logits,
the reference length's held at 0, its curvature is `Q·(diag(s) − s sᵀ)`, with `Q` the prior's draws
over every length and `s` the shares of the other lengths (`allele_prior_curvature`,
`penalised_fit.rs`). It enters P, so H and the middle term both get it, as the stutter prior's does
(§2.10.4); putting it in H alone gave the same coverage to within 0.01. Without it, H was not the
curvature of the objective the fit maximises, and at thin classes it was not positive: the reads'
likelihood alone curves upward along a trade between two rare allele lengths, which only the prior
holds. Every error of the motif length was then withheld, in 89 of D1's 480 coverage runs
(`reports/implementations/str_penalised_fit_thin_errors_diagnosis_2026-10-15.md`); with the
curvature, all 89 invert. The fitted numbers do not change. **The rule that withholds every error of
a motif length when H cannot be inverted stays, as a fallback no measured case reaches**
(`UndeterminedErrors`).

**What was measured** (`compare_D_prior_2026-10-14.md` §3, truth within one error, honest about 68 in
100):

- stutter and F: about honest, **except F at 10 samples × 30 reads on two-base runs: 41 in 100**;
- allele lengths of 5% or more: honest;
- **not honest:** lengths of 1–5% from 60 samples up (42–55 in 100); lengths under 1% from 60 samples
  up (27–42); the homopolymer shorter share at 3 reads (29–55);
- at the ends of the repeat range at 1 sample × 3 reads, the two-base slip rate's errors held the
  truth within two errors in 2 draws of 6 (`track2_penalised_em_2026-10-12.md` §10.3);
- **at one sample, the two-base direction numbers with the stutter prior on** *(accepted by the owner,
  2026-10-15, as the known cost of §2.10; no code change)*. The prior pulls them toward its centre
  (0.5 for the shorter share, 0.05 for the fall-off), and the errors do not widen to cover the pull: a
  penalised fit's error measures the spread of the estimate, not how far the penalty moved it. So the
  file can state a shorter share of 0.55 ± 0.08 where the truth is 0.8. Measured by step C8's review on
  the bench, counting classes whose fitted number lies more than two of its own errors from the truth,
  with the prior | without (an honest error misses about 1 time in 22):

  | cell (1 sample, PCR truth unless named) | number | with the prior | without |
  |---|---|---|---|
  | two-base × 3 reads, PCR-free (3 seeds) | shorter share | 31/63 | 18/62 |
  | two-base × 3 reads (6 seeds) | fall-off | 28/126 | 18/126 |
  | two-base × 10 reads (6 seeds) | shorter share | 21/126 | 12/126 |

  The corner was already weak without the prior (18 of 62 is about 6 times what honest errors give),
  so this is not a new kind of defect. In the PCR-free cell the prior also removes a failure: one seed's
  shorter share swung between 0 and 1 along the repeat count without it. The effect on calls was not
  measured; in the example above, at a slip rate of 0.01, the share of reads one repeat short goes
  from 0.008 to 0.0055. Plan step D1 counts these classes again, for the owner.

Rare lengths get no error: a share expected at less than one locus, share × loci with reads < 1
(§7; it was 3 until 2026-10-16).

**Where a thin class's errors start** *(owner, 2026-10-16, after plan step D1 measured them; §7)*:

| number | error reported from | constant |
|---|---|---|
| shorter share, fall-off | 3 loci with reads | `LOCI_FOR_A_DIRECTION_ERROR` |
| slip rate, F, the allele-length shares | 8 loci with reads | `LOCI_FOR_AN_ERROR` |

Measured on drawn cohorts of 2,000 samples, at 3 to 7 loci the truth lay within two of the shorter
share's and the fall-off's errors 89 to 93 times in 100, against 95 for an honest error; within two
of the slip rate's only 72 to 79 times, because the rate itself comes out low there (below the truth
in 81 to 86 classes of 100, by a median of 1.5–2.8%; §2.4's floor), not because its error is narrow.
Below 3 loci no number gets an error. The other rules of the error table (arch §1.4) — too few
slipped reads, an error that is absent, zero or wider than `NOT_PLACED` — apply after these.

*History:* from 2026-10-14 every number of a class with fewer than 8 loci with reads got no error
*(lead; confirmed by the owner the same day)*, because thin classes' errors had never been measured:
the bench drew 8 to 4,947 loci a class (`phase0_bench_2026-10-12.md:55`) and track 1 started at 20, and
the sandwich's B has rank at most the number of loci. On kimura, 107 of 141 classes had fewer than 8
loci with reads (`kimura_shapes_span9_2026-10-14.md`).

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
| prior under the three stutter curves: weight (`STUTTER_PRIOR_READS`) | one read's worth a class, each curve | §2.10.2; owner 2026-10-15, kept after C8's three-seed sweep of 0, 1 and 5 on held-out loci |
| prior under the three stutter curves: centre (`STUTTER_PRIOR_CENTRE`) | slip rate 0.10, shorter share 0.5, fall-off 0.05: calling's default stutter | §2.10; owner 2026-10-15; not a measurement |
| F from pairs, at most | 10 alternations | bench; soft |
| kink knee / tail | 2.5 / 2 | §2.3 |
| stiffness held once it moves less than | 0.1 natural-log units | measured, track 2 §10.1 |
| a class whose numbers stopped moving reuses its counts | on | measured, every number unchanged to the third decimal |
| clustering weight on allele draws | off | measured harmful at large cohorts, track 2 §2.5 |
| errors withheld below | 3 loci with reads a class for the shorter share and the fall-off, 8 for the slip rate, F and the allele-length shares; a share under one locus's worth; a slip rate, shorter share or fall-off under one expected slipped read (the prior's weight) | §2.6, §2.10.4, §7 (owner 2026-10-16); soft |

### 2.8 Two exact speed-ups, kept

- **Each distinct row of reads is scored once.** A sample's contribution depends only on its row of
  read counts and its homozygote excess; at 3 reads tomato2 has 86 rows for each distinct one.
- **Each pair of samples is four numbers**, the coefficients of a cubic in the concentration, in place
  of a 2.9 kB table.

Both change no number beyond about 10⁻⁶ (`cost_2026-10-14.md`). **Their gain falls with depth**: at
30 reads a row repeats 3.6 times, at 300 reads 1.2 times.

### 2.9 What the run logs

Per motif length: rounds, the five stiffness strengths and where bends were eased, the time in each
step, peak memory. Per class: rows with reads, distinct rows, and whether F was absent. The count of
withheld errors (§7). These replace the climb's summaries (`ssr_fit.rs:3165-3332`). *(Owner,
2026-10-15:)* the stutter prior and the classes resting on it (§2.10.5).

### 2.10 A weak prior under the three stutter curves — approved by the owner 2026-10-15

**Status.** The direction is the owner's (2026-10-15): a weak prior on the stutter curves, as §2.4
gives the allele lengths. **The form, the weight and the centre below were proposed by the lead and
approved by the owner on 2026-10-15**; the owner's answer to each open question is recorded in
§2.10.9. Line numbers in this section are at commit `47e21b8c` (step C7), except §2.10.4's, which
point at the code as built (re-pointed after the review of commits `0b12941d..53a3e4f9`).

#### 2.10.1 The problem it fixes

**Nothing holds a slip rate, a shorter share or a fall-off away from 0 or 1 except two numerical
guards.** Each round draws the three stutter curves through each class's expected counts as binomials
on the logit scale (`penalised_fit.rs:1071-1102`). Only two things stop a logit from running off:

- a ridge of 10⁻⁶ on every logit (`RIDGE`, `:291`);
- a clamp at ±40 on every Newton step (`LARGEST_UNBOUNDED`, `:295`, applied at `:2226`).

**Zero is an absorbing boundary of EM, for stutter exactly as §2.4 found for allele lengths.** At a
thin class, a read one repeat off its sample's allele can be explained as a slipped read or as a
second allele. As the slip rate falls, step 1 hands fewer reads to "slipped", and the next drawing of
the curve lowers the rate again, until a guard stops it. The direction numbers do the same: a class
whose few slipped reads were all shorter has its shorter share driven to 1.

**What the milestone C review measured** (2026-10-15, on the oracle: 4 tomato2 accessions, 20
regions, about 5 reads a sample; the review is `doc/devel/reports/reviews/str_penalised_fit_milestone_C_2026-10-15.md` on this branch, not
committed):

- 12 classes report slip rates from 1.7×10⁻³ down to 4.2×10⁻¹⁸, none of them resting on a slipped
  read.
- **Motif lengths with one or two classes** (five bases at 5 repeats; six bases at 4 and 5) have no
  bend to penalise (§2.3). With no slipped read, each rate is where the ridge balances the reads,
  `n·p = 10⁻⁶·|logit p|` for `n` reads: 9.3×10⁻⁷, 2.7×10⁻⁷ and 1.9×10⁻⁶. **Each is a function of
  the read count and of the constant 10⁻⁶, not an estimate.**
- **The three-base slip-rate curve became a straight logit line** (strength 1.9×10⁹), anchored by the
  one class with a slipped read (6 repeats, 1.0 expected slipped read) and falling 1.465 logit units
  a repeat to −40.00 at 30 repeats, the clamp.
- **Direction numbers at one slipped read**: shorter shares of 0.99998 and 0.99999 and fall-offs of
  3.5×10⁻⁵ and 1.1×10⁻⁵ (three bases at 6 repeats; four bases at 7).

**What calling does with those numbers.** Calling scores a read as `(1 − λ)·[the genotype's emission]
+ λ·(lengths the read allows)/R`, with λ = 0.2 (`DEFAULT_OUTLIER_WEIGHT`,
`calling/likelihood/ssr.rs:158`; the term at `:699`) and R the reachable lengths (36 for one
three-base candidate at 6 repeats). So every slip rate below about 10⁻⁴ scores as no stutter. Take a
homozygote with one read one repeat shorter than its allele. **On likelihood alone, before the
genotype prior, a heterozygote is preferred until the homozygote has more than about 6.2 reads at its
own length under no stutter, and about 3.4 under calling's default stutter**
(`StutterModel::hipstr_shipped`, `alignment/stutter.rs:310-319`). Tomato sits at 3 to 5 reads a
sample. And a class with no reads at all gets that default, 10% slipped (§2.2), so in the oracle a
three-base class with one locus and one read scored at 10⁻¹⁰ while its neighbour with none scored at
10⁻¹.

#### 2.10.2 What the prior is

**Each class's three binomials get one extra trial, split as a stated centre says, before the curves
are drawn.** In words: every class is treated as if it had one more read, which slipped with
probability 0.10; one more slipped read, which was shorter with probability 0.5; and one more
repeat slipped, which went one repeat further with probability 0.05. In counts, at every repeat count
with a class (`a` the weight, the centre `(p₀, h₀, f₀)`):

| curve | successes | trials | the prior adds |
|---|---|---|---|
| slip rate | slipped reads | slipped + exact reads | `a·p₀` slipped, `a·(1 − p₀)` exact |
| shorter share | shorter slips | shorter + longer slips | `a·h₀` shorter, `a·(1 − h₀)` longer |
| fall-off | repeats slipped beyond the first (`steps − stepped`) | repeats slipped (`steps`) | `a·f₀` further, `a·(1 − f₀)` stopped |

This is a Beta prior on each class's probability, with `a·p₀` and `a·(1 − p₀)` as its two
parameters. On the logit scale the fit already works on, it is one more binomial term of the same
shape as the reads' (`CurveData::Binomial`, `penalised_fit.rs:1712-1727`), so it needs no new solver.

**The design** *(lead's proposal; approved by the owner 2026-10-15, §2.10.9)*:

- **Weight `a`: one read's worth** (`STUTTER_PRIOR_READS = 1.0`), the same for all three curves and
  every class, whatever the cohort size or depth.
- **Centre: calling's default stutter, written in the fit's three numbers** (`STUTTER_PRIOR_CENTRE`):
  slip rate 0.10, shorter share 0.5, fall-off 0.05. `hipstr_shipped` puts 0.05 of reads one repeat
  short and 0.05 one long, with a one-step share of 0.95 (`alignment/stutter.rs:310-319`); calling
  turns the fit's numbers into its own by `fall-off = 1 − one-step share`
  (`calling/likelihood/stutter_rates.rs:86-104`).

**What one read's worth does, class by class.** Alone (with no curve), a class with `y` slipped reads
among `n` reads gets the slip rate `(y + 0.1)/(n + 1)`. So:

- **the prior moves any class by at most 1/(n + 1) of the distance to 0.10**: a class's own reads
  outweigh it from its second read on;
- **with no slipped read the rate is 0.1/(n + 1), never 0**, and it falls as the reads accumulate.
  Ten thousand reads with no slip give 10⁻⁵: genuinely low stutter, measured on many reads, stays low;
- **with no read at all the class is refused and calling uses 0.10 (§2.2); with one read and no
  slip it gets 0.05**, with three reads 0.025. The rate now changes smoothly with the reads behind it,
  where today it drops from 0.10 to 10⁻⁶ at the first read.

**Why the weight is one read, and why it is soft.** It follows §2.4's first floor (one allele draw a
class, measured), and one read is the smallest weight that changes calling at the thinnest classes
(case A below). A heavier prior biases moderate classes with genuinely low stutter: at five reads'
worth, a class with 1 slipped read among 1,000 reads would come out at 1.5×10⁻³ where the reads say
10⁻³, half again its value. The owner chose one read (2026-10-15) before any weight was measured
against held-out loci; step C8 of the plan then measured it.

**What the sweep measured, and the owner's answer** *(owner, 2026-10-15, after C8's review)*. Held-out
log-likelihood on the bench (`fit_bench.rs`'s scorer, about 900 held-out loci a cell, 5 folds),
weight 0 or 5 minus weight 1, at seeds 1 / 2 / 3, PCR truth, 3 reads a sample:

| cell | weight 0 − weight 1 | weight 5 − weight 1 |
|---|---|---|
| one-base, 1 sample | −0.62 / −0.98 / −0.70 | −0.16 / +1.06 / +0.33 |
| one-base, 10 samples | −0.82 / +0.13 / +0.03 | −1.44 / −0.12 / −0.74 |
| one-base, 60 samples | −1.25 / −0.26 / +0.08 | +0.23 / −0.44 / +0.43 |
| two-base, 1 sample | −0.33 / −1.39 / −0.26 | −6.13 / −3.50 / −5.49 |
| two-base, 10 samples | +0.49 / +0.49 / +0.36 | −4.43 / −3.24 / −4.22 |
| two-base, 60 samples | +0.25 / +0.21 / +0.50 | −0.76 / −1.11 / +0.11 |

Three results hold at every seed; the rest change sign from seed to seed:

- weight 5 loses 3.2–6.1 log units on two-base runs at 1 and 10 samples;
- weight 1 beats weight 0 at one sample, in both motif lengths, by 0.3–1.4 log units;
- weight 0 beats weight 1 on two-base runs at 10 and 60 samples, by 0.2–0.5 log units over about 900
  held-out loci, under one part in 1,700 a locus. This fits the pull on two-base direction numbers
  that §2.6 records.

(The 500 samples × 30 reads cell was run at seed 1 only; its differences were 0.00–0.03.) **The owner
kept the weight at one read**: weight 0 is the weight that lets thin classes collapse to near-zero
slip rates (§2.10.1), and its gain where data are plentiful is small. The weight remains a soft
setting; the sweep compared three values, not a range.

**Why that centre.** It is what calling already uses for a class with no reads, so the fit and calling
agree on the limit of no evidence **for the whole-repeat shares**: the centre gives 0.05 of reads one
repeat short and 0.05 one long, as `hipstr_shipped` does. The part-repeat shares do not agree. Calling
derives them from the fit's slip rate as `PART_REPEAT_SHARE_OF_WHOLE` (0.05) of it
(`stutter_rates_for`, `calling/likelihood/stutter_rates.rs`), which at the centre gives 0.0025 each
way, where `hipstr_shipped` has 0.01. The calling-side test pins the whole-repeat shares only. The
centre is not a measurement: `hipstr_shipped` is "HipSTR's shipped
constants and not a fit … on no organism in particular" (`calling/parameters_file/defaults.rs:114-121`).
At one read's worth that matters only where a class has a handful of reads.

#### 2.10.3 How it works with the rest of the fit

- **The bend penalty is unchanged and acts on top.** The prior is part of each class's own data; the
  curve then pulls thin classes toward their neighbours as before. In a run of thin classes the prior
  counts once per class, so a stiff curve through ten thin classes carries ten reads' worth of prior
  against those classes' pooled reads. Because every class's term now has a finite best value, **the
  best line through a tail of classes with no slipped read has a finite slope**, and the line to the
  clamp of §2.10.1 cannot form.
- **Only repeat counts with a class get the prior.** Repeat counts between classes have no counts and
  take the curve's value, as today. Adding the prior there would make them look like classes with
  reads (next bullet).
- **The prior changes how many points "have counts", and with it how a curve is fitted** (trap).
  `refit` counts the points whose trials are above zero (`has_counts`, `penalised_fit.rs:1807-1813`,
  counted at `:1921`) and, at two or fewer, fits a straight line at the top strength (`:1931`). Today
  a class with no slipped read has no trials for the shorter share and the fall-off, so a motif length
  whose slips sit in two classes is fitted as a line. With the prior every class has a trial. **This
  is intended**: it is exactly the line that ran to the clamp. Count classes, not points with trials,
  for the "one or two classes" rule of §2.3.
- **The stiffness search sees the prior as data.** The strength is chosen on the expected counts as if
  observed (§2.3); it is now chosen on the counts plus the prior. Expect the three-base slip-rate
  strength of 1.9×10⁹ to come down, since no class pulls toward −∞ any more.
- **Direction numbers when no read slipped:** the shorter share is the centre's 0.5 and the fall-off
  its 0.05, the values calling's default holds. With one slipped read, one repeat shorter, they are
  (1 + 0.5)/2 = 0.75 and 0.05/2 = 0.025, where today they are 0.99998 and 3.5×10⁻⁵. The fall-off's
  clamp at 0.99 (`:1160`) stays.
- **The ridge and the clamp stay**, as guards for F and the allele lengths. On the stutter numbers
  they no longer set anything: the prior's pull on a logit is about `a·p₀` = 0.1 a class, the ridge's
  at most 10⁻⁶ × 40.
- **The class's expected slipped reads are written without the prior** (trap). The file writes each
  rate's `expected_slipped_reads` (`calling/parameters_file/to_toml.rs:1242-1245`), and the run log
  and §2.6's error rules read it. It stays the reads' own count: a class with none says 0, not 0.1.
- **The start does not change.** The loop still starts at 0.05, 0.5, 0.3 (§2.7); the centre is not a
  starting point, and changing the start would move the fixture for no reason.
- **Allele lengths and F** are not given a new prior; they move only through step 1 seeing different
  stutter.
- **Cost:** three additions a class a round. Nothing measurable.

#### 2.10.4 Standard errors (§2.6)

- **The prior enters the errors as the bend penalty does: in H and in the middle term.** Its
  curvature at a class is `a·p(1 − p)` for the slip rate (and the same form for the other two),
  added at that class's position in P (`penalised_fit.rs:1643-1656`), which `penalised_sandwich`
  adds to both H and B (`:2699-2713`). Leaving it out of the middle term would give a number with no
  reads behind it an error near zero from B alone, which the review's first finding shows is the
  dangerous direction.
- **A shorter share or fall-off whose class has fewer than one expected slipped read (the prior's
  weight) gets no error** (`NotPlaced`) *(owner, 2026-10-15)*: its value leans on the centre and on
  its neighbours on the curve as much as on its reads, and its error would be close to the prior's
  own width (2 logit units for the shorter share), not a measurement. The shorter share is judged
  on its trials, the slipped reads with a direction, which are never more than the slipped reads
  and equal them unless a slipped read could not be given a direction. **The fall-off is judged on
  the slipped reads, not on its trials, the repeats slipped**, which a read that slipped two
  repeats or more brings to one while the slipped reads are still below it. *Corrected 2026-10-15
  after the review of commits `0b12941d..53a3e4f9`: this bullet and the code first judged the
  fall-off on its repeats slipped, against the owner's words; a lone class of 10 loci × 3 reads
  with one read two repeats short (0.99997 expected slipped reads, 3.1 repeats slipped) was written
  with a fall-off error.*
- **So does a slip rate whose class has fewer expected slipped reads than the prior's weight**
  (`NotPlaced`) *(owner, 2026-10-15, after C8's review)*. The rule counts the slip rate's successes,
  not its trials. *Corrected: this section first said that the slip rate needed no rule, because
  alone it would get a logit error of `1/√(0.1·(1 − p))` ≈ 3.2, which the `NOT_PLACED` test (3 on the
  logit scale) withholds. That figure is the curvature's alone, `H^(-1/2)`. The reported error is the
  sandwich, and with no slipped read its middle term collapses:* B sums each locus's squared score,
  and with zero events every locus's score is near zero. Take a class of 10 loci × 3 reads with no
  slip, at weight 1: the rate is 0.1/31 = 0.0032, H ≈ 31·p(1 − p) ≈ 0.1, B ≈ 10 × (3p)² ≈ 9.4×10⁻⁴,
  P = p(1 − p) ≈ 3.2×10⁻³, so the variance is (B + P)/H² ≈ 0.42 and the logit error 0.65. The class
  has 8 loci or more and 0.65 is under 3, so without this rule the file would say 0.0032 `Estimated`,
  about 0.0017–0.006 within one error, when no slip in 30 reads is consistent with any rate up to
  about 0.1. With the rule, case C2 of §2.10.6 (10,000 reads, none slipped) carries its 10⁻⁵ with no
  error: zero events place no logit error.
- Classes with fewer than 3 loci withhold every error, and classes with fewer than 8 the slip
  rate's (§2.6), before these rules apply; so the rules decide the slip rate's error from 8 loci
  and the shorter share's and fall-off's from 3. They matter most for classes with genuinely no
  slips, such as a long-motif class at a large cohort.
- **The prior adds a positive diagonal at every class**, which removes the exactly flat directions
  the review blamed for failed factorisations (`newton`, `penalised_fit.rs:2514`, `:2575-2583`). It
  is too small to fix the review's first finding (about 0.01 a class against strengths of 10¹⁰),
  which needs its own fix. *Note, 2026-10-15: the joint inversions measured to fail since then
  (`doc/devel/reports/implementations/str_penalised_fit_thin_errors_diagnosis_2026-10-15.md`) fail
  not along a stiff stutter curve but along a trade between two rare allele lengths at a thin
  class, because H leaves out the curvature of the allele-length prior of §2.4; the stutter numbers
  alone had positive curvature in every such run.* *Done 2026-10-16 (plan step D7): H and the middle
  term now hold that curvature too (§2.6), and every one of those runs inverts.*

#### 2.10.5 What the run log and the parameters file say

- **The run log states the prior once**, in the stage's opening line: "stutter: a prior of one read's
  worth a class, centred on calling's default (slip rate 0.10, shorter share 0.5, fall-off 0.05)".
  Per motif length, it counts the classes with fewer expected slipped reads than the prior's weight:
  "in N of the M classes fitted, fewer than one read slipped, too few to place their shorter share
  and fall-off, which carry no error and lean on the stutter prior and on any neighbouring repeat
  counts on the curve" (`how_the_fit_ran`, `penalised_fit/run_log.rs`). *Corrected 2026-10-15 after
  C8's review: this section said such numbers "rest on the stutter prior", which is true of a class
  fitted alone; a thin class on a curve with well-read neighbours takes most of its value from them.
  Corrected again after the review of commits `0b12941d..53a3e4f9`: the line said those numbers come
  "from their neighbours on the curve and the stutter prior, not from their own reads", which a
  class with no neighbours contradicts; one read two repeats short moved such a class's fall-off
  from the prior's 0.05 to 0.525.*
- **The parameters file gets no new key** *(owner, 2026-10-15; §2.10.9, question 3)*. §2.4's
  floor is not written either. A reader can tell a prior-led number from its row: it carries
  `expected_slipped_reads` below one and, by §2.10.4, no error on any of its three stutter numbers.
  The file's own note on stutter provenance (the `repeat_tracts` note in `to_toml.rs`) gains one
  sentence: a stratum on a penalised curve was fitted with one extra read's worth of evidence that
  behaves as calling's default stutter does, so a row with `expected_slipped_reads` below one
  carries no error on any of its three numbers, and its shorter share and fall-off lean on that
  default and on its neighbours on the curve, where it has any, as much as on its own reads. The
  weight and the centre in that sentence are built from `STUTTER_PRIOR_READS` and
  `STUTTER_PRIOR_CENTRE`, so the note cannot state other numbers than the fit uses. *That sentence
  is the lead's reading of "no key" (2026-10-15); the owner may still drop it.*

#### 2.10.6 Worked cases

Each case is one class fitted alone, as at a motif length with one or two classes; a curve through
neighbours moves them further. The calling column is how many reads at the allele's own length a
homozygote needs, on likelihood alone, before it beats the heterozygote that one read one repeat short
suggests, with calling's terms as in §2.10.1 (`tmp/stutter_prior/arith.py` on this branch). Fewer
reads needed means the short read is explained as stutter sooner.

| case | reads (slipped) | slip rate today → with the prior | shorter share; fall-off | calling: reads needed |
|---|---|---|---|---|
| **A.** 1 sample × 3 reads, a motif length with one class | 3 (0) | 4.1×10⁻⁶ (the ridge) → **0.025** | 0.5; 0.05 (today set by no read) | 6.2 → **4.9** (default stutter: 3.4) |
| **B.** the oracle's three bases at 6 repeats, alone | about 135 (1) | 7.4×10⁻³ → **8.1×10⁻³** | 0.99998 → **0.75**; 3.5×10⁻⁵ → **0.025** | 5.2 → **5.4** |
| **C1.** no slip among 10 reads | 10 (0) | 1.4×10⁻⁶ → **9.1×10⁻³** | 0.5; 0.05 | 6.2 → **5.6** |
| **C2.** no slip among 10,000 reads | 10,000 (0) | 2.0×10⁻⁹ → **1.0×10⁻⁵** | 0.5; 0.05 | 6.2 → 6.2 (no stutter, as the reads say) |
| **D.** genuinely low stutter, true rate 10⁻³ | 1,000 (1) / 10,000 (10) / 30 million (30,000) | the reads' own → +10% / +1% / +3 parts in a million | the reads' own, moved by at most 1/(slipped + 1) toward 0.5 and 0.05 | unchanged |
| **E.** 2,000 samples × 300 reads, one locus, no slip | 600,000 (0) | 4.0×10⁻¹¹ (the ridge) → **1.7×10⁻⁷** | 0.5; 0.05 | unchanged |

Read the table this way:

- **Case D is the owner's condition**: a class with genuinely low stutter and many reads keeps its low
  rate. The +10% at 1,000 reads is a tenth of the uncertainty one slipped read leaves (its own logit
  error is about 1), and the shift halves with every doubling of reads.
- **Cases A and C1 are the thin corner**: the rate is now what the prior and the few reads give
  together, and calling needs about one read fewer to call the homozygote. It does not reach the
  default stutter's 3.4, because three reads with no slip are weak evidence that stutter is low, and
  the prior is weak by design.
- **Case B moves calling the other way, slightly**: a shorter share of 0.75 rather than 1 means a read
  one repeat short is a less likely slip, so the homozygote needs 5.4 reads, not 5.2. One slipped read
  cannot say that all slips are contractions, and the number now says so.
- **The oracle's other near-zero rows**, alone: five bases at 5 repeats (about 15 reads, no slip),
  9.3×10⁻⁷ → about 6×10⁻³; six bases at 4 (about 56 reads) 2.7×10⁻⁷ → about 1.8×10⁻³; six bases at 5
  (about 7 reads) 1.9×10⁻⁶ → about 0.013. The read counts are inferred from the ridge's balance, so
  they are approximate.
- **At one sample** nothing new is needed: the prior is per class and does not use pairs of samples.
  **At thousands of samples** every class with reads has thousands of reads, and the prior's share is
  below one part in a thousand except at classes with no slip, where it sets the rate's order of
  magnitude (case E) and the error rules withhold the errors of all three stutter numbers (§2.10.4).

#### 2.10.7 What it changes in calling, at the cases the review measured

- **No rate in the file is below `0.1/(n + 1)` for its class's reads alone, and no logit sits at ±40.**
  Calling's "below about 10⁻⁴ is no stutter" now applies only to classes whose reads, in thousands,
  showed no slip.
- **The three-base patchwork becomes a gradient.** Today the classes from 7 to 30 repeats score at
  10⁻³ to 10⁻¹⁸, and those with no row at 0.10. With the prior each scores between its own reads'
  answer and 0.10, closer to 0.10 the fewer reads it has; a class with no reads still scores at the
  default. A class with one read and no slip scores at most 0.05 against the empty class's 0.10, a
  factor of two where today it is a factor of 10¹⁷.
- **Genotypes change only where a homozygote has 4 to 6 reads with one read one repeat off, at a thin
  class.** The oracle's records in three- and six-base tracts (four of them; QUAL fell in all four
  under C7) are where to look first.
- **Every class's direction numbers stop at the boundary no longer**, so calling stops treating a
  one-read class as having only contractions or only one-step slips.

#### 2.10.8 Considered and not chosen

- **A normal prior on the logit, centred on the default.** Its curvature does not vanish at low rates,
  so it pulls a genuinely low, well-measured rate up far more: at a true rate of 10⁻⁴ measured on
  100,000 reads (10 slipped), a standard deviation of 2 logit units moves it up 17%, where one
  read's worth of pseudo-counts moves it 1%. That fails the owner's condition at large cohorts.
- **One locus's worth of reads**, as §2.4's second floor sizes the allele draws. The unit that
  repeats independently for allele lengths is the locus; for stutter a read slips on its own. A class
  with one locus and 10,000 reads would weigh the prior as heavily as all of its reads.
- **A centre taken from the data** (the motif length's pooled stutter, or the cohort's). At a motif
  length with one class it is that class, so it holds nothing away from zero; pooled across motif
  lengths it is dominated by homopolymers, whose stutter is ten times a trinucleotide's. The curve
  already does the borrowing within a motif length.
- **The review's other options, a floor or "absent at a boundary"**, treat the output, not the cause.
  §2.10.4's error rule keeps the useful part of "absent": a number the reads did not place carries no
  error.

#### 2.10.9 The questions put to the owner, and the answers (2026-10-15)

The lead put five questions with a recommendation each; the owner answered all five on 2026-10-15.

1. **Weight: one read's worth, or more?** *Recommended: one.* **Owner: one read's worth a class.**
   Step C8's held-out comparison on the bench at weights 0, 1 and 5 checks it; if a weight other
   than one predicts better, that result goes back to the owner rather than changing the constant.
   *It did, in part: weight 0 predicts slightly better on two-base cohorts of 10 and 60 samples. The
   owner kept one read (2026-10-15); the sweep and the reason are in §2.10.2.*
2. **Centre: calling's default (0.10, 0.5, 0.05), or the fit's start (0.05, 0.5, 0.3)?**
   *Recommended: calling's default.* **Owner: calling's default stutter** (`hipstr_shipped`), so that
   a class with one read and a class with none agree in the limit; the start is a starting point,
   not a belief.
3. **Write the prior into the parameters file?** *Recommended: no key; the file's note says it.*
   **Owner: no new key.** The lead keeps the one explanatory sentence in the file's stutter note
   (§2.10.5); that is the lead's reading, and the owner may still drop it.
4. **The error rule of §2.10.4** (no error on a shorter share or fall-off whose class has fewer than
   one expected slipped read)? *Recommended: adopt.* **Owner: yes.** Without it a shorter share of
   0.5 that no read placed would be written with an error of 2 logit units, `Estimated`. *Extended
   by the owner after C8's review (2026-10-15) to the slip rate, counted by its slipped reads
   (§2.10.4).*
5. **Keep `RIDGE` at 10⁻⁶ for F and the allele lengths?** *Recommended: keep.* **Owner: keep.** With
   the prior the ridge no longer sets any stutter number, and changing it would move F and the
   allele lengths for a reason unrelated to this fix.

---

## 3. What will bite

- **Output must be bit-identical across runs and thread counts.** The oracle
  (`scripts/promote_ng_oracle.sh`) pins checksums. The prototype differs by up to 2×10⁻¹⁰ between two
  runs, most likely from sums taken in thread order (`cost_2026-10-14.md`, "uncertain"). One such sum
  is outside `penalised_fit.rs`: `SamplePairs::log_likelihood` sums with rayon's `.sum()`
  (`composite_fit.rs:1583-1590`), called at every grid point of F's table (`penalised_fit.rs:1081`).
  Sum in a fixed order, as `ssr_fit.rs` already does elsewhere, in every module the fit calls. **A
  2×10⁻¹⁰ difference can flip a discrete decision** (a strength held below 0.1 log units, a class's
  counts reused below the tolerance, the pair loop stopping below 0.01; `penalised_fit.rs:387`,
  `:575-578`), so a comparison before and after the fix can differ for that reason too.
- **The output must also not depend on the order the psps are named in, and today it does** *(found
  2026-10-15; the owner chose to fix it on this branch before milestone D, plan step C10)*.
  `estimate-parameters` numbers the cohort's samples and read groups in the order its psps were
  named, and the fit's sums over them run in that numbering; floating-point addition is not
  associative, so a re-ordered cohort sums differently in the last places. The same four tomato2
  psps fitted in the oracle's order (the VCF's) and in alphabetical order give two different files
  (md5 `4c4efa9…` against `5c53b9dc…`; `tmp/thread_determinism_report.md` on this branch, not
  committed). The first difference is in the SNP/indel fit, which stops at its pass limit rather than
  converged and so carries it into its returned values: the inbreeding coefficient of SRS1839219 is
  0.96645 against 0.96617. Through the repeat-tract fit, whose curve strengths near the top of their
  search sit on an almost flat marginal likelihood, it becomes a two-base fall-off strength of 5.7×10⁹
  against 1.2. Thread count is not the cause: the oracle's order gives `4c4efa9…` at 4 threads and
  at 18. The specs already require order independence: `parameter_prepass_cohort.md` test 7 ("The
  gather is order-independent — same summaries, same parameters, whatever order they are read in")
  and `parameters_file.md` §6 (the sample list: "Order is immaterial and must be"). **The defect
  predates this branch and is on main too**; the fix (put the cohort in one fixed order before
  anything is summed; architecture §2.4a) reaches main when this branch merges.
- **The renaming is half done** on this branch (`aa42d629`): the parameters file's keys were renamed
  without a format version or reading of older files, so this branch refuses every parameters file
  written before it ("unknown field `length_spectrum_by_period`" on kimura). Fixed by
  `str_fit_to_calling.md`; build the two together.
- **Calling panics on a non-positive concentration** (`stratum_fits.rs:1051-1055`) and stores a class's
  shares and concentration as one pair (`:379-385`). A class with shares but an absent F has nowhere
  to go until the companion spec's change lands.
- **Memory at deep reads.** The prototype holds every class's distinct rows at once during a round:
  4.8 GB at 2,000 samples × 30 reads against 1.1 GB for the counting fit, which fits a class's rounds
  one at a time; scaled up it would not fit 64 GiB at kimura's scale (`cost_2026-10-14.md` §9.1).
  **One class at a time (§2.2) is the bound; "one class per worker" is not**: on 18 threads that is
  every class of a motif length. And a parallel loop over classes with a parallel loop over rows
  inside it lets a thread waiting on its inner work pick up another class and build a second table;
  the climb measured exactly this, "six tables at four strata at once"
  (`every_stratum_from_every_start`, `ssr_fit.rs:2736-2740`).
- **The evidence is already the largest memory**: `gather_strata` holds every stratum at once
  (`ssr_fit.rs:4863-4887`), 19.4 GiB on kimura for 66.2 million pairs with reads, about 315 bytes a
  pair against tomato2's 139. Not changed here (§6).
- **The provenance written beside each stutter number assumes main's curves** (`LevelSmoothing`,
  `parameters_file/mod.rs:1199`, carries a fitted `SlippageCurve`; `LevelSource` is cell / curve /
  blend, `slippage_curve.rs:519-531`). A penalised curve is neither; the companion spec settles what
  is written.
- **Deleting the climb takes more than the climb.** The only producer of the period-pooled rung
  (`fit_period_allele_length_variations`, `ssr_fit.rs:3392`) calls the climb's walk (`:3429`); the
  multiplier prototype uses its quadrature (§1.4); the bench's methods `A` and `B_curves` call it
  (`fit_bench.rs:284`, `:660`); and the prototype reads its allele window from `SsrFitConfig`
  (`penalised_fit.rs:331`). What goes and what stays is settled in the architecture (§4) *(lead,
  2026-10-14)*.
- **`--str-param-estimates-at-once`** (`cli/estimate_parameters.rs:144-145`) sizes the climb's
  likelihood tables; it means nothing to this fit. Decided in §7.
- **Three cross-platform digest tests fail on this branch already** (`cli::cross_platform_digests`):
  they pin the parameters file's checksum, which the ±9 window and the rename moved without
  re-recording. Re-record them with the oracle (§4.6), not before.
- **A stale comment** in `fit_strata_with` says thin strata get "two shares copied from a neighbour"
  (`ssr_fit.rs:3006-3010`); nothing copies a neighbour any more (`:2955`). Do not build on it.

---

## 4. How we know it works

1. **Unit tests carried from the prototype:** with no penalty the fit equals the counting fit class by
   class; a drawn kink on many reads survives while a thin class borrows; the banded solve equals a
   dense one; the distinct-row and cubic-pair shortcuts equal the direct computation; F is absent at
   one sample. **Added:** F is absent when two samples share one locus and its estimate stops at an
   end of its search (§2.5); a motif length with one class, and one with two, fits and reports
   (§2.3). **Added 2026-10-15 (approved by the owner, §2.10):** a one-class motif length with no slipped read
   returns a slip rate of 0.1/(n + 1) for its n reads, a shorter share of 0.5 and a fall-off of 0.05,
   with all three errors `NotPlaced`; no stutter logit of any test ends within one unit of ±40; a
   drawn class with genuinely low stutter and many reads keeps its rate within its error.
2. **The bench reproduces the recommendation's numbers** for the production settings at the cells of
   `compare_D_prior_2026-10-14.md` (same seeds), within run-to-run noise, plus the cells §1.2 says
   were not run (500 × 30, 2,000 × 10, three seeds at every 2,000-sample cell), and the coverage of
   errors on classes with fewer than 8 loci (§2.6).
3. **Determinism:** two runs, and runs at 1 and 18 threads, give the same bits; so do two orderings
   of the same psps (§3; plan step C10).
4. **Real data, held-out:** on tomato2, held-out prediction by whole locus, the new fit against the
   counting fit, with the scorer of the bench (`fit_bench.rs`). This is the first test on real reads.
5. **Cost on kimura, measured:** the fit alone within one hour and the run's peak memory reported, at
   2,169 samples.
6. **The oracle moves, and the move is explained** before its checksums are re-recorded
   (`scripts/promote_ng_oracle.sh`; the fitted file and the calls made with it).
7. **Calling does not get worse** on HG002 against its repeat-tract truth and on the GIAB trio.

---

## 5. Open

1. **Does a slip rate 0.5% of its value low change any call?** Measure at 500 and 2,000 simulated
   samples × 3 reads. If it does, revisit the floor (§2.4). *Leaning:* it does not; settled by that
   measurement. Nothing in the project simulates reads for 500 or 2,000 samples, so the measurement is
   made on the bench's read counts: each sample's genotype weights under the fitted stutter against
   under the true stutter. That is a proxy for calling, not calling *(lead, 2026-10-14)*.

## 6. Deferred, with a recommended home

- **Stutter per slippage group and library multipliers** → `spec/str_library_stutter.md`, rebuilt on
  this fit; the climb-based prototype is deleted with the climb (§1.4).
- **Gathering evidence one motif length at a time**, to bound the evidence's memory → the census-read
  step (`gather_strata`), once kimura's 19.4 GiB matters; the reason for reading the whole band
  (neighbour borrowing, `census.rs:2239-2245`) no longer holds.
- **The read guard dropping 45,432 of kimura's 86,688 loci** (`kimura_shapes_span9_2026-10-14.md`) →
  its own question; not this fit's.
- **A parametric allele-length shape** → research avenue C, if thin classes need more than the
  curves give.
- **Fitting the period-pooled rung** → nothing fits it once the climb goes; a file that carries it
  still reads (companion §3.4). Whether to refit it with this fit is the library-multiplier work's or
  a later spec's question.

## 7. Decisions

- **2026-10-14, owner:** adopt penalised EM with the one-locus floor.
- **2026-10-14, owner:** the CLI flag `--str-param-estimates-at-once` is removed, with an error
  naming why, rather than accepted and ignored.
- **2026-10-14, owner: errors on rare allele lengths are reported as absent.** A share under 3 loci's
  worth (share × loci with reads < 3) gets no error, and the log says how many were withheld; such
  errors were measured too small (§2.6). The 3 rests on two settings measured by track 1; soft.
  [*Lowered to 1 by the owner, 2026-10-16, below.*]
- **2026-10-14, owner: the overshoot at the ends of the repeat range at 1 sample × 3 reads is
  accepted** (+49% to +127% at the extreme class, `track2_penalised_em_2026-10-12.md` §10.3), and
  stated as a known limit of the fit.
- **2026-10-14, owner: the climb's code is deleted** once §4.4 and §4.7 pass, with its tests and
  examples, in one commit.
- **2026-10-14, owner:** stutter per slippage group and library multipliers go to the
  library-multiplier work (§6).

### Decisions by the lead in response to the review, 2026-10-14 — for the owner to confirm

- **Errors withheld for every number of a class with fewer than 8 loci with reads** (§2.6), until
  measured; the plan's first check measures their coverage, and the owner may lift the rule.
  **Confirmed by the owner, 2026-10-14.** [*Narrowed by the owner, 2026-10-16, below: the shorter
  share and the fall-off get an error from 3 loci.*]
- **Classes visited one at a time within a motif length, motif lengths one at a time**, parallel
  inside a class (§2.2, §3). Thread scaling and wall time are measured again at full scale.
- **A class with no reads stays refused**, as today; the new fit produces no `Derived` class (§2.2).
- **F absent also at an end of its search** (§2.5), the counting fit's rule.
- **What goes with the climb** (§3; architecture §4): the period-pooled rung stays readable from a
  file and loses its only producer; the multiplier prototype and the bench's climb methods are
  deleted, their results kept in the research reports; `SsrFitConfig` is folded into the new fit's
  settings.
- **The 0.5% question is measured on the bench's read counts**, as a proxy for calling (§5).

### 2026-10-15: the stutter numbers' lower bound

- **Owner:** the slip rate and its two direction numbers get a weak prior on their curves, as §2.4
  gives the allele lengths, removing the absorbing boundary at zero rather than flooring the output;
  it must not mistake genuinely low stutter on many reads for missing data (§2.10.1).
- **Proposed by the lead and approved by the owner, 2026-10-15** (§2.10; the questions and answers
  in §2.10.9):
  - the form: one extra binomial trial a class on each of the three curves, split as the centre says;
  - **owner: the weight is one read's worth**, whatever the cohort or depth; step C8's sweep of 0, 1
    and 5 on held-out loci checks it, and a different winner goes back to the owner;
  - **owner: the centre is calling's default stutter** (`hipstr_shipped`): slip rate 0.10, shorter
    share 0.5, fall-off 0.05;
  - the prior in the errors' H and middle term; **owner: no error (`NotPlaced`) on a shorter share or
    fall-off whose class has fewer than one expected slipped read**;
  - **owner: no new key in the parameters file.** The run log states the prior, and the lead keeps one
    sentence in the file's stutter note (the lead's reading; the owner may drop it);
  - **owner: `RIDGE` stays at 10⁻⁶** for F and the allele lengths.
- **Owner, 2026-10-15: the stated F's median counts only classes with at least `LOCI_FOR_AN_ERROR`
  (8) loci with reads**, falling back to every class with an F when none has 8. Recorded in the
  companion, [`str_fit_to_calling.md`](str_fit_to_calling.md) §3.5; built in plan step C9.

### 2026-10-15: the owner's answers to step C8's review

Each as the lead recommended.

- **Owner: a slip rate whose class has fewer expected slipped reads than the prior's weight gets no
  error (`NotPlaced`)**, like its shorter share and fall-off (§2.10.4). It corrects §2.10.4's
  assumption that the `NOT_PLACED` threshold would withhold it: the sandwich's middle term collapses
  with zero slipped reads, so a class of 30 reads with no slip got a logit error of 0.65.
- **Owner: at one sample, the prior's pull on two-base direction numbers is accepted** as the known
  cost of §2.10 and recorded in §2.6's list of errors that are not honest; plan step D1 measures it
  again. No code change.
- **Owner: the prior's weight stays one read.** The three-seed sweep is in §2.10.2. Weight 0 is the
  weight that lets thin classes collapse to near-zero slip rates; its gain where data are plentiful
  (0.2–0.5 log units over about 900 held-out loci, two-base, 10 and 60 samples) is small.
- **Owner: the fit's dependence on the order the psps are named in is fixed on this branch before
  milestone D** (§3), as plan step C10: the cohort is put in one fixed order, samples by name and read
  groups within a sample by ID, before anything is summed. The defect is on main too.

### 2026-10-16: the owner's answers at checkpoint D

Each as the lead recommended (`reports/implementations/str_penalised_fit_checkpoint_D_2026-10-16.md`).

- **Owner: D2 and D5 pass; the climb is deleted** (plan milestone E).
- **Owner: the errors' `H` and middle term gain the allele-length prior's curvature** (§2.4's floor:
  one locus's worth of draws spread over the lengths, plus the half draw a class [*correction, review of D7: one draw a class, `ALLELE_PSEUDO_DRAWS_A_CLASS`*]), as §2.10.4 already
  does for the stutter prior. Without it `H` is not invertible at thin classes, where the reads'
  likelihood alone curves upward along a trade between two rare lengths, and every error of the motif
  length was withheld (`reports/implementations/str_penalised_fit_thin_errors_diagnosis_2026-10-15.md`:
  89 of D1's 480 coverage runs). Fitted values do not change. The rule that withholds every error
  when `H` still cannot be inverted stays as the fallback. Plan step D7.
- **Owner: where errors start, as D1 recommends** (`reports/implementations/str_penalised_fit_D1_2026-10-15.md`
  §1): the three stutter numbers (slip rate, shorter share, fall-off) get an error from **3** loci
  with reads; F and the allele-length shares keep **8**; a share's rule (`LOCI_FOR_A_SHARE_ERROR`,
  share × loci with reads) drops from 3 to **1**. The stated F's median (C9, companion §3.5) keeps
  counting classes with at least 8 loci. Plan step D8.
- **Owner, later the same day: the slip rate keeps the line at 8 loci**; only the shorter share and
  the fall-off get an error from 3. D7 recovered the errors of 33 runs at 2,000 samples that D1's
  coverage had not counted; with them, at 3 to 7 loci and 2,000 samples, the slip rate's error holds
  the truth within two errors 72 to 79 times in 100 (the shorter share's and fall-off's 89 to 93),
  because the slip rate itself is biased low there by a median of 1.5–2.8% (§2.4's floor), not
  because its error is narrow (`tmp/d7/coverage_by_number.md`, quoted in plan step D8).
  *Done 2026-10-16 (plan step D8): `LOCI_FOR_A_DIRECTION_ERROR` = 3, `LOCI_FOR_AN_ERROR` = 8,
  `LOCI_FOR_A_SHARE_ERROR` = 1; §2.6, §2.7 and arch §1.4 state the rules.*
- **Owner: the SNP/indel fit's convergence is a separate job** after this branch merges (the
  performance review's finding L5); not on this branch.

## 8. Reuse map

| what | existing code | reuse |
|---|---|---|
| the loop, curves, stiffness, errors | `ssr_fit/penalised_fit.rs` | promote; settings per §2.7 |
| step 1 (weigh genotypes) and the pair likelihood | `ssr_fit/composite_fit.rs` (`expected_counts`, `SamplePairs`) | promote |
| F's placement rule | `concentration_is_placed`, `composite_fit.rs:83-98` | as is (§2.5) |
| read probabilities per allele | `Slippage::read_probabilities`, `ssr_fit.rs` | as is |
| evidence | `the_tract_strata_of_a_cohort`, `gather_strata` | as is |
| into calling's tables | `StratumOutcome` → `StratumFits::over` | as is, plus the companion spec |
| the held-out scorer and the simulator | `bench_fixtures/{fit_bench,cohort_draw}.rs` | the acceptance tests of §4 |
| the parity oracle | the bench's `D_tract_prior` method on `str-fit-compare` | §4.2 |
