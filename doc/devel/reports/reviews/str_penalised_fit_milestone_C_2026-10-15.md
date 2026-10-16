# Milestone C review: the penalised fit made production

**Scope.** `git diff a540b074..47e21b8c`, the eight steps C0 to C7, read across steps. The per-step
reviews (`tmp/C0_review.md` to `tmp/C7_review.md`) are not repeated here. I ran no cargo, as
instructed. Every number below comes from reading the code or from the two oracle runs:
`tmp/oracle_C6/` (before C7, the climb) and `tmp/oracle_C7/` (after C7, the penalised fit). Both
oracle runs use 4 tomato2 accessions and 20 regions, at about 5 reads a sample. The small scripts
that produced the tables are in `tmp/mC/`.

**Verdict.** The wiring is correct:

- `estimate-parameters` passes `PenalisedFitConfig::production()` (`cli/estimate_parameters.rs:382`)
  into `fit_a_cohort`, and from there into `fit_repeat_tracts`.
- Every read group is pooled into one slippage group (`estimate_parameters.rs:436`).
- The run log's counts agree with the outcomes.
- I found no path where the fit panics on real input.

Three findings come from the oracle's fitted file:

- **The standard errors written for two- and three-base classes are far too small** (Blocker). At
  two-base runs they are 3 to 90 times smaller than the reads could support. At the one three-base
  class with errors they are 160 to 2 million times smaller. They are written as estimated.
- **The near-zero slip rates are not set by the reads** (Major, Ask). They are set by two numerical
  guards in `penalised_fit.rs`: a ridge of 10⁻⁶ and a clamp at ±40 on the logit scale. Calling
  reads every one of them as "no stutter".
- **A failed factorisation inside the curve fit is swallowed** (Major). When it happens, the choice
  of how stiff a curve is can drift to the top of its range.

The large moves in F and in the curves' tails are mostly what the design expects (§2 below).

---

## 1. The near-zero slip rates

### 1.1 Why the slip rate can reach 10⁻¹⁸

**The slip rate has no lower bound except numerical guards.** Each round fits the slip-rate curve
on the logit scale (`penalised_fit.rs:1086-1090`), and only two things hold the values back:

- `RIDGE = 1e-6` (`:291`) adds `½·10⁻⁶·η²` to the objective, where η is the logit of the slip rate.
  Its documentation says what it is for: "so that a number with no counts behind it … stays
  finite".
- `LARGEST_UNBOUNDED = 40` (`:295`) clamps each Newton step to ±40 (`:2226`). The documentation
  calls this "beyond which every probability here is zero or one to working precision".

The fall-off has a clamp at 0.99 and F has a search range (`:1160`, `:1162-1168`). **The slip rate
has neither.**

**The fit drifts toward zero because, for the slip rate, zero is an absorbing boundary of EM.** At a
thin class, a read of a different length can be explained either as a slipped read or as a second
allele. As the slip rate falls, the E-step hands fewer reads to "slipped". The next M-step then
lowers the rate again, and the rate keeps shrinking until a guard stops it. This is the same
mechanism spec §2.4 describes for allele lengths ("EM has an absorbing boundary at zero"). §2.4
puts a floor under allele lengths. Nothing does that for stutter.

**Two different guards set the values in the oracle file.**

- **A motif length with one class, or two:** spec §2.3 says such a motif length "has no bend to
  penalise", so its curve is the class's own counts. With zero slipped reads, the binomial term and
  the ridge balance where `n·p = 10⁻⁶·|logit p|`, with `n` the reads. That gives p ≈ 2.6×10⁻⁶ at
  5 reads and 1.6×10⁻⁷ at 100 reads. All four oracle values in this group fit that formula:
  - the five-base class at 5 repeats, 9.34×10⁻⁷, implies about 15 reads, from 7 (locus, sample)
    pairs;
  - the six-base classes at 4 and 5 repeats, 2.7×10⁻⁷ and 1.9×10⁻⁶;
  - C7's two-sample, one-locus test fixture, 5.8×10⁻⁷, implies about 25 reads.

  **None of these values is an estimate. Each is a function of the read count and of the constant
  10⁻⁶.**
- **A motif length with three or more classes, at its far end:** the three-base slip-rate curve came
  out with strength 1.9×10⁹, which makes it a straight line on the logit scale. Its logits fall by
  exactly 1.465 per repeat: −4.89 at 6 repeats, −6.36 at 7, −7.82 at 8, and so on, down to
  **−40.00 at 30 repeats, which is the clamp** (logistic(−40) = 4.248×10⁻¹⁸, the value in the
  file). Only the class at 6 repeats has a slipped read (1.0 expected). Every class above it has
  close to none, so the best line has slope −∞. **The slope written to the file is the one that
  makes the line hit −40 at the class furthest away.** Without the clamp, every three-base value
  from 7 to 30 repeats would be lower still.

### 1.2 Every near-zero row in the oracle

All loci counts come from `tmp/oracle_C7/fit.log`. Every row below has fewer than 8 loci with
reads, so **all of its errors are withheld for too few loci** (`TooFewLoci`, spec §2.6). The one
exception is the three-base class at 6 repeats, which is the anchor of the line.

| motif length, repeats | loci with reads | (locus, sample) pairs | slip rate | expected slipped reads | what sets the value |
|---|---|---|---|---|---|
| 3, 6 | 11 | 40 | 7.4×10⁻³ | 1.0 | the reads (errors estimated; see B1) |
| 3, 7 | 5 | 19 | 1.7×10⁻³ | 5.5×10⁻⁴ | the line to the clamp |
| 3, 8 | 4 | 14 | 4.0×10⁻⁴ | 0.021 | the line to the clamp |
| 3, 9 | 1 | 4 | 9.3×10⁻⁵ | 4×10⁻⁶ | the line to the clamp |
| 3, 10 | 1 | 3 | 2.2×10⁻⁵ | 2×10⁻⁶ | the line to the clamp |
| 3, 11 | 1 | 3 | 5.0×10⁻⁶ | 8×10⁻⁷ | the line to the clamp |
| 3, 17 | 1 | 1 | 7.7×10⁻¹⁰ | 2×10⁻¹³ | the line to the clamp |
| 3, 28 | 1 | 3 | 7.9×10⁻¹⁷ | 3×10⁻¹⁹ | the line to the clamp |
| 3, 30 | 1 | 1 | 4.2×10⁻¹⁸ | 5×10⁻²² | **the clamp itself** (logit −40.00) |
| 5, 5 | 2 | 7 | 9.3×10⁻⁷ | 1.2×10⁻⁷ | the ridge (one class) |
| 6, 4 | 5 | 19 | 2.7×10⁻⁷ | 2.8×10⁻⁸ | the ridge (two classes) |
| 6, 5 | 1 | 3 | 1.9×10⁻⁶ | 5×10⁻¹² | the ridge (two classes) |
| 2, 31 | 1 | 1 | 1.7×10⁻³ | 2×10⁻³ | the curve's slope from thin classes |
| 2, 36 | 1 | 1 | 2.4×10⁻⁴ | 3×10⁻⁵ | the curve's slope from thin classes |

**The direction numbers reach their boundaries too.** At three bases and 6 repeats, and at four
bases and 7 repeats, one slipped read each gives a shorter share of 0.99998 and 0.99999 and a
fall-off of 3.5×10⁻⁵ and 1.1×10⁻⁵. The fall-off is clamped to a one-step share of 0.99 in calling
(`GEOM_MAX`, `alignment/stutter.rs:108`). A shorter share of 1 means calling never expects a read
one repeat longer.

### 1.3 What calling does with 10⁻¹⁷

**Calling treats every slip rate below about 10⁻³ as no stutter, so 10⁻¹⁷, 10⁻⁷ and 10⁻⁴ all give
the same calls.** The reason is how calling scores a read against a genotype
(`calling/likelihood/ssr.rs:697-740`). The probability of a read is

`(1 − λ)·[the genotype's alleles' emission] + λ·(lengths the read allows)/R`,

with:

- λ = 0.20, the outlier weight (`DEFAULT_OUTLIER_WEIGHT`, `:158`): the share of reads assumed to
  come from somewhere else and to show any reachable length;
- R, the number of tract lengths the candidates can reach. For one three-base candidate at 6
  repeats, R is 36 (`reachable_length_changes_of`, `alignment/stutter.rs:631-654`), so the outlier
  term is about 0.0056 per read.

There is no floor on the slip rate in calling (`stutter_rates_for`, `stutter_rates.rs:86-104`).
Take a homozygote and one read one repeat shorter than its allele. That read scores:

| slip rate scored | stutter term (0.8·s·shorter·one-step) | + outlier term | ln-cost of the read against s = 0 |
|---|---|---|---|
| any s below 10⁻⁴ (all the near-zero rows) | under 8×10⁻⁵ | 0.0056 | under 0.02 nats |
| 7.4×10⁻³ (three bases at 6 repeats) | 0.0059 | 0.0115 | 0.72 nats |
| `hipstr_shipped`, 0.05 shorter × 0.95 | 0.038 | 0.0436 | 2.05 nats |

**What that costs at low depth.** Under "no stutter", one read a repeat short makes a
heterozygote preferred, on likelihood alone and before the prior, until the homozygote has about 6
reads at its own length. Under the shipped default the threshold is about 3 reads. Tomato sits at
3 to 5 reads a sample, so this decision falls exactly in that range. At 300 reads it matters less
per read, but the count of stutter reads grows in proportion.

**Calling's stutter jumps between "none" and HipSTR's 10% as a candidate allele grows by one
repeat.** A candidate is scored at the stutter of its own allele's class. A class with no row in
the file gets `StutterModel::hipstr_shipped()`, marked `Defaulted` (spec §2.2;
`repeat_tract_parameters.rs:27-31`). In this oracle the three-base classes therefore score as
follows:

- 7 to 11 repeats: the fitted near-zero values;
- 12 to 16: 10%, because they have no row (14 has no reads, so it is refused);
- 17: 7.7×10⁻¹⁰;
- 18 to 27: 10%;
- 28 and 30: near zero.

**So a class with one locus and one read gets less stutter than a class with no reads at all.**

**What the oracle saw.** Four STR records lie in three- or six-base tracts. QUAL fell in all four,
by 3.0, 5.4, 14.9 and 12.0, and no genotype changed. Over all STR records, QUAL moved in 426 records
and 2 genotypes changed, but at 4 samples and 20 regions none of that can be attributed to a
particular row.

### 1.4 The options

The decision is the owner's. These are the options, what each changes, and which section it
amends.

1. **Leave it as designed.** Nothing in the spec forbids it: §2.2 fits every class with reads, and
   §2.3 lets a one- or two-class motif length keep its own counts. What it costs:
   - the thin corner writes rates that no read placed;
   - every one of them is "no stutter" in calling;
   - they move if `RIDGE` or `LARGEST_UNBOUNDED` changes, which are numerical constants;
   - calling gets the patchwork described in §1.3.

   If kept, §3 ("What will bite") should say so, and spec §1.3 goal 3 ("a number the data cannot
   place is absent, not a fitted bound") would be knowingly broken for stutter.
2. **Mark a stutter number absent when it stops at a boundary**, as spec §2.5 does for F. Calling
   then scores the class with `hipstr_shipped`, marked `Defaulted`. The rule has to be chosen:
   - "within 1% of the clamp" catches only the class at 30 repeats, 1 of these 13 rows;
   - "the class's own expected slipped reads below 1" catches 11 of them; it misses the anchor at
     6 repeats (1.0) and the four-base class (1.0, rate 0.023).

   It amends §2.5 (a new bullet for stutter), goal 3 in §1.3, and `str_fit_to_calling.md`, which
   has to let a `Fitted` class carry `slippage: vec![None]`. Today that slot means "this group put
   no read here", so calling's `NoSlippage` reasons need a new case. The fixture does not move.
3. **A weak prior on the three stutter curves, replacing the 10⁻⁶ ridge**, for example a fraction of
   a read at the starting values of §2.7 (0.05, 0.5, 0.3), as §2.4 does for allele lengths. This
   removes the absorbing boundary at its cause:
   - zero-slip thin classes settle near the prior or their neighbours;
   - the clamp-driven line goes away;
   - rich classes move by about one part in (their reads ÷ the prior's weight).

   It amends §2.3 or §2.4 and the table in §2.7. The C0 fixture moves and must be re-recorded, and
   D1's comparison with the bench report has to be re-measured, because §2.4's measured bias was
   taken without it.
4. **A floor**, for example a slip rate of at least 10⁻³. This is the smallest change. It
   contradicts goal 3 unless the floored value is marked as floored. It amends §2.7 and §1.3.

Options 2 and 3 can be combined: absent where no read places the rate, and a prior so that the
curve's slope is not set by a guard. Whichever is chosen, the shorter share and fall-off at one
slipped read (§1.2) are the same case.

---

## 2. The other moves in the fitted file

**F on classes both fits wrote: mostly expected.** The climb's 0.015 at two bases and 10 repeats,
and its 0.0022 at 12 repeats, were its extreme values on 8-locus classes. An F near 0 means each
locus carries the class's whole spread of lengths. That is implausible in 4 inbred accessions, and
the climb's own neighbouring classes say otherwise: 0.57 at 8 repeats and 0.83 at 9. The penalised
fit reads F off pairs of samples at a locus (spec §2.5) and smooths it along the repeat count. That
gives 0.663 and 0.693, with errors of 0.08 and 0.10.

The one-base class at 15 repeats moved from 0.805 to 0.568. It has 13 loci and an error of 0.11, and
its neighbours at 12 to 14 repeats sit at 0.54 to 0.56. A move of two errors toward the neighbours
is what a curve of strength 20 does to a 13-locus class. **None of this points to a defect.**

The two-base F curve's strength came out at exactly 1.0×10¹⁰, the top of the range searched, with
19 classes. That may be what the data say, or it may be M2 below. It cannot be told apart without
the check that M2 asks for.

**The fallback F, 0.526 → 0.723: about half of the move comes from which classes are counted.** The
fallback is the median F over every class that has one (`median_concentration`,
`stratum_fits.rs:1250`). The climb wrote F for 15 classes. The penalised fit writes it for 46, and
31 of those have fewer than 8 loci. For a class like three bases at 30 repeats (1 locus, 1 sample),
F is a point on the curve beyond the last class with a pair, yet it gets one vote, the same as the
one-base class at 8 repeats with 1,214 loci. Over the same 15 classes the climb had, the new median
is 0.615.

So 0.09 of the move comes from the F values themselves, and 0.11 from the thin classes joining the
median. The fallback is now applied to the classes without F (four bases at 7 repeats, five at 5,
six at 5) and to tracts with no row. It is a design consequence, filed below as m1.

**The climb's flat borrowed tails became curves: expected in kind, but the tails' shapes come from
the same mechanism as §1.** The one-base slip-rate curve falls from 0.12 at 14 repeats to 0.041 at
21 and rises to 0.117 at 25. The two-base curve peaks at 0.23 at 16 repeats and falls to 2.4×10⁻⁴
at 36. In both, the falling stretch rests on classes whose own expected slipped reads are between
0.002 and 1.5, at 1 to 4 loci. Those classes pull the curve down for the reason in §1.1, against
the known rise of stutter with tract length. This is not a code defect. It is a property of the
model in the thin corner, and option 3 in §1.4 is the one that addresses it.

---

## 3. Findings

### Blocker

- `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:2301-2334` (called at `:1409`) —
  **[Blocker]** Standard errors at motif lengths with stiff curves are far too small, and they are
  written as estimated.
- **Confidence:** High that the errors are wrong. Medium-high on the mechanism.
- **Problem:** the comparison below is between what the file says and the smallest error the reads
  allow. A point on a smooth curve cannot be more precise than all of its motif length's reads
  pooled into one number. On the logit scale that pooled error is `1/√(slipped reads × h(1−h))`,
  where h is the shorter share (taken as ½ for the level). From the oracle file (`tmp/mC/se.py`):

  | class | reported logit error: level / shorter share / fall-off | smallest error the reads allow (pooled over the motif length) |
  |---|---|---|
  | one base, 8 to 15 repeats | 0.09–0.27 / 0.11–0.32 / 0.32–0.46 | level 0.045, shorter 0.091 (and each matches its own class's binomial, e.g. 0.108 against 0.104) |
  | two bases, 6 repeats (3 slipped reads) | 0.043 / 0.0042 / 0.0044 | level 0.127, shorter 0.374 |
  | two bases, 12 repeats (6 slipped reads) | 0.056 / 0.0031 / 0.0033 | level 0.127, shorter 0.366 |
  | three bases, 6 repeats (1 slipped read in the whole motif length) | 0.0061 / 0.00011 / 0.00012 | level 0.99, shorter about 200 |

  One-base errors are honest. Two-base errors are 3 times too small for the level and 90 times too
  small for the direction numbers. At three bases they are 160 times and 2 million times too small:
  the file states the shorter share of the three-base class at 6 repeats, measured from one slipped
  read, as 0.99998 ± 2.5×10⁻⁹.

  **The likely mechanism.** `penalised_sandwich` first tries to invert `H = A + P + 10⁻⁶·I`, where A
  is the curvature of the reads, P is the bend penalties and I is the identity matrix. If that
  fails, it retries with ridges of 10⁻⁹ to 10⁻³ times `scale`, the mean diagonal of A + P. At two and
  three bases, four of the five curves have strengths of 10⁹ to 10¹⁰, so `scale` is about 10¹⁰, and
  even the first fallback ridge (about 25 to 60 a number) is larger than the reads' own curvature (1
  to 10). Where the first inversion fails, the ridge, not the reads, sets H⁻¹, and the error shrinks
  in proportion: a ridge of 57 a point gives a level error near 0.03 at two bases, and 0.043 is
  reported.

  Why the first inversion fails: some direction has almost no information. Examples are the slope of
  a line through one class's data, or a share at logit ±10. That puts H's condition number near
  10¹⁷, beyond double precision. One-base strengths are lower (`scale` about 3×10⁸), so a fallback
  there would be invisible, which matches the one-base rows being honest.
- **Why it matters:** spec §1.3 goal 2 promises an honest error wherever the data can support one.
  These errors pass the `NOT_PLACED` test (3 on the logit scale), so they are written as
  `Estimated`, and a reader takes them as measurements. The bench (C0, D1) cannot see this, because
  its cells draw dense classes whose strengths stay low.
- **Suggested fix:**
  1. Verify first: have `penalised_sandwich` return which ridge it used, and print it per motif
     length on the oracle data. Expected: 10⁻⁶ at one base, a larger one at two and three bases.
  2. Then pick one of these (Ask; it is a design choice):
     - (a) scale the fallback ridge to A's mean diagonal, not A + P's;
     - (b) invert in two parts, the null space of P (the lines) and its range, so the 10¹⁰ never
       meets the reads' 1;
     - (c) when any ridge above 10⁻⁶ is needed, return `None` so that every number is `NotPlaced`
       rather than shrunk.

     (c) is a one-line change that is safe now. (a) or (b) gives errors back.
  3. Add a test: a drawn motif length whose curve strength reaches the top of the range, checking
     each reported logit error against the pooled bound above.

### Major

- `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:287-295, 1086-1090, 2226` — **[Major,
  Ask]** Slip rates near zero are set by `RIDGE` and `LARGEST_UNBOUNDED`, not by reads, and calling
  scores them as "no stutter".
- **Confidence:** High.
- **Problem:** see §1. In this 4-sample oracle, 12 classes write slip rates from 1.7×10⁻³ down to
  4.2×10⁻¹⁸, and none of them rests on a slipped read. Each value is either the ridge equilibrium
  (`n·p = 10⁻⁶·|logit p|`) or a straight line driven to the clamp at logit −40.
- **Why it matters:** this is the confident-and-wrong case CLAUDE.md names for a single sample at low
  depth. It also gives calling a stutter that jumps between none and 10% from one repeat count to the
  next (§1.3).
- **Suggested fix:** the owner chooses from §1.4. Whatever is chosen, add a test that pins the
  behaviour of a one-class motif length with zero slipped reads, which today returns about 10⁻⁶ by
  the ridge.

- `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:2180, 2211-2213` — **[Major, Apply]**
  A failed Cholesky factorisation in `newton` is silent: the fit returns its starting values with a
  log-determinant of 0, and the choice of strength can then favour the top of the range.
- **Confidence:** Medium. The code path is certain. That it fires on real data is not shown.
- **Problem:** `let mut ln_determinant = 0.0;` and then `let Some(factor) = matrix.cholesky() else {
  break; };`. If the first factorisation fails, `newton` returns `from` unchanged and
  `ln_determinant = 0`. In `fit_robustly_at` (`:2126-2127`) the marginal likelihood is
  `… + ½·rank·ln λ − ½·ln|H|`. The missing `ln|H|` (about rank·ln λ at a large λ) leaves a term that
  grows with λ, so that strength wins the search. The band matrix at λ = 10¹⁰ with a ridge of 10⁻⁶ is
  at the limit of double precision: rounding in the last pivots is about λ·10⁻¹⁶ = 10⁻⁶, the size of
  the ridge. So any curve with a direction no reads inform, such as a line through one class's data,
  is at risk. Three of the five two-base strengths in the oracle are exactly 1.0×10¹⁰.
- **Why it matters:** a swallowed failure that changes a fitted number with nothing in the log. It
  also feeds B1, because strengths at the top of the range make `scale` 10¹⁰.
- **Suggested fix:** on a failed factorisation, return a marginal of `f64::NEG_INFINITY` from
  `fit_robustly_at`, so that strength is never chosen. Inside `newton`, retry with the ridge raised
  tenfold, up to a limit, before giving up. Count the failures per curve and print the count in the
  run log's stiffness line. To check whether it fires on tomato2, run the oracle fit with a debug
  counter.

### Minor

- `src/parameter_estimation/joint/stratum_fits.rs:1250` and
  `penalised_fit/run_log.rs` (`class_line`, the `"F fitted"` branch) — **[Minor, Ask]** The fallback
  F and the log treat a class's F that is pure curve extrapolation as fitted.
  - **Problem:** the fallback is the median over every class with an F. 31 of the 46 votes come from
    classes with fewer than 8 loci, some with no pair of samples at all, such as three bases at 30
    repeats. Their median is 0.723. The 15 classes with 8 or more loci give 0.615. The log prints
    "F fitted" for a class whose only locus has one sample.
  - **Fix:** either weight the median by loci with reads, or take it over classes with at least
    `LOCI_FOR_AN_ERROR` loci (`fitted_here` either way). Have the log say "F from its motif
    length's curve; no pair of samples here" when the class has no locus with two samples. This
    amends `str_fit_to_calling.md`'s definition of the fallback.
- `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:1368-1369, 2301-2334` — **[Minor,
  Apply]** The standard errors build dense matrices of `(points × width)²` floats across every
  repeat count between the lowest and highest class, empty ones included. The memory note at
  `:443-450` mentions only the per-class blocks of about 900 floats.
  - **Problem:** at two bases in the oracle there are 42 points × 21 numbers = 882, or 0.78 million
    floats per matrix and about nine such matrices. That is fine. But a motif length whose classes
    span 4 to 500 repeats (long reads crossing long tracts) gives size 10,000: 800 MB per matrix and
    about 10¹² operations in `dense_inverse` and `dense_product`, single-threaded, up to five times
    over in the ridge loop.
  - **Fix:** state the bound in the doc comment. Better, build the matrices over only the points with
    a class (empty points are pinned by P alone and carry no error), or use the band structure of
    A + P.

### Nits

- `penalised_fit.rs:1-5`: the module doc still opens with "a research prototype … being promoted".
  After C7 it is the caller's fit.

---

## 4. Across steps: panics, the log, the settings

**No panic on real input found.**

| input | what happens |
|---|---|
| empty cohort | `fit_repeat_tracts` logs one line and returns an empty vector (`:605-613`, `:654-657`) |
| one sample | `between_locus_share_measurable` is false and `pair_tables` returns `None`, so the loop breaks; F is `None`, the errors never index F (`:1228`, `:1472`); tested end to end (`from_run_parameters.rs`, the one-sample trip) |
| a motif length whose classes all lack reads | `Period::new` returns `None`; every class is `Refused { NoSpanningReads }`; `curves.expect` is reached only for fitted classes (`:742`) |
| every sample's reads identical | one distinct row; values may reach ±40, where `logistic` gives finite values; a failed dense inverse gives `NotPlaced`; a shorter share of exactly 1.0 (logistic(40) rounds to 1) passes the file's `[0, 1]` check (`validate.rs:1242`) and `an_error_of` returns `NotPlaced` on a zero slope |
| huge depth | counts are `f64`; nothing on the fit's path indexes by depth |

The remaining asserts hold by construction:

- more than one slippage group (`:659`): the CLI pools every read group;
- two classes at one repeat count (`:936`): strata are keyed by class;
- more than one motif length in a fit (`:478`): grouped by motif length at `:668`.

**The run log matches the outcomes.**

- 34 classes have fewer than 8 loci: 9 at one base, 13 at two, 8 at three, and one each at four and
  five, plus two at six.
- 34 × 23 numbers (3 stutter, 19 lengths, F) = 782. Less the 3 F's that are absent and so
  `NotPlaced` (four bases at 7 repeats, five at 5, six at 5), that is 779, as logged.
- The remaining 257 withheld errors are shares in the 15 thicker classes.
- "settled" is printed only when the fit stopped before 2,000 rounds (`run_log.rs:167`).

**`production()` is what C7 passes:** `repeat_tract_config()` returns it, and no option alters it
(`estimate_parameters.rs:382-383`). The census-fit tests use it too. `--str-param-estimates-at-once`
is refused before anything is read (`:387-396`, `:403`).
