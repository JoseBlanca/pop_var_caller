# ng step 4, the STR path — D4: the search, and what it must report beside its answer

*Implementation report, 2026-08-12. Step D4 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md) — **own commit, do not
bundle**, and Milestone D's last step. One review agent, which ran the sharp control itself and
found four things needing a fix. Design authority:
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §3,
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.2. Ported from
`maximise_slip`, `starting_points` and `fit_from_starts` in
[`examples/shared/stutter_model.rs`](../../../../examples/shared/stutter_model.rs).*

## What the step is

`fit_by_multistart` in a new shared module `fitting/multistart.rs`: maximise the noise parameters
from several starting points, climbing the genotype frequencies at every trial, and return the
best-scoring **with every start's outcome beside it**. Each start is followed by a coordinate-wise
golden-section search on the scale each axis lives on — a rate on a log scale, two shares on logit
scales, so the search cannot walk out of a share's range and resolves a rate to a fixed fraction of
itself.

A ladder is the wrong instrument here for two reasons the spec gives: a flat scan over three
parameters is 4.2 million scores per (read group × stratum) with several hundred strata a sample,
and a quarter-Phred spacing is the wrong ruler for two parameters that are shares in `(0, 1)`.

## Recorded deviations from the architecture

1. **`fitting/` gained a third thing, not the two the plan anticipated.** Alongside the widened
   `fit_mixture_weights` (D3) and the new driver, it needed a **trait**: `SearchableNoise`, which
   says how a path's noise parameters are read, written and bounded on the scale each axis is
   searched on. The driver is generic over `NoiseModel`, and there is no way to move a point in
   `M::NoiseParams` without one. **This is the answer to the question Checkpoint D asks** — the
   seam was cut in nearly the right place, and the piece it was missing is the one that says what
   the parameter space *is*, as against how a cell is scored.
2. **`score_one_candidate` was extracted out of `ladder_scan.rs`** so the ladder and the search
   score a candidate through the same code rather than two copies of it. Behaviour-identical: the
   reviewer checked it assertion by assertion, and the generic path's 476 tests are unchanged.
3. **The driver reports the spread; it does not judge it.** The architecture has it flag a fit whose
   starts spanned more than `START_AGREEMENT_LIMIT`. That constant is the STR path's, and how far
   apart two answers may sit before a fit is disowned is a path-specific judgement, so the shared
   driver returns the ratio and the caller compares.
4. **A stalled inner climb is reported rather than asserted.** The sibling scan treats it as a
   debug-time bug, and that is right for three genotypes; a stratum here has 66 to 91, most of them
   heading to a frequency of zero, and reaching stillness on those takes far longer than the shared
   cap allows. `MultistartResult::every_climb_settled` carries it.

## What the review changed

**Blocker — the warm start was measured to fire zero times, and its prose claimed otherwise.** I
had added it, from the reference, to stop the inner climb exhausting its pass cap. The reviewer
instrumented it: **668 climbs, 668 rejections, not one warm start taken.** The reason is
mechanical rather than fixture-specific — 30 of the stratum's 36 genotypes have a true frequency of
zero, their responsibilities underflow to exactly zero after a few hundred passes, and the previous
answer is then no longer the strictly-interior point the climb requires. A production stratum has
more genotypes and fewer occupied lengths, so it would reject there too. **Removed**, with what was
learned written where the next reader will find it: the idea needs a floor *inside* the climb, not
a filter outside it. Its filter's comment was also wrong about the mechanism — an unfiltered warm
start does not silently pin a genotype at zero, it panics in `check_start`.

**Blocker — the spread is near one by construction, and nothing said so.** Each axis is
line-searched over its *whole* range at every sweep, so a start's own value on the headline axis is
overwritten before it is ever read. Measured on the control, three of the four starts reach a
bit-identical point and the spread is exactly 1.00. That is the search working — and it is exactly
the reading this module's own preamble warns against. `headline_spread`'s doc now says it is a
floor on the disagreement rather than four independent searches agreeing.

**Major — three of the four starting points collapsed onto one level.** I floored the starts at the
search's own lower bound, 1e-5, where the reference floors at 1e-4. At 1e-5 a stratum whose reads
all sit at the reference length gives an estimate of zero and three of the four multipliers clamp
to the same point — so the set would still be four starts and would no longer be a spread. Floored
at 1e-4 as the reference does. The same function now refuses a level estimate that is not a number,
which is what dividing by a stratum with no whole-repeat reads gives.

**Major — `spread_across` returned `NaN` where every start reached zero**, because `−∞ − (−∞)` is
`NaN` and the reference's `if low > 0.0` guard was dropped in the move to log coordinates. A `NaN`
spread compares `false` against any limit, so the fit would have been reported as *agreeing*. Fixed
by asking the trait for the headline parameter **on its natural scale** and taking the ratio
directly, which also removes a latent coupling: `exp(Δcoordinate)` is the ratio only if axis 0 is
logarithmic, and on a logit axis it is the odds ratio — 12.17 where the ratio is 8.64, a 41% error
in the one number a caller gates on.

**Major — the search returned its last iterate where the architecture asks for its best.** The
bracket's midpoint is where a line search stopped, not the best point it saw, so on an axis that is
not unimodal a sweep can move downhill and the search would report that. Now tracked and returned.

**⛦ The reviewer also found a trap in this repo's own tooling**, worth carrying: `scripts/dev.sh`
always passes `-w $PROJECT_DIR`, so `cd tmp/wt && ../../scripts/dev.sh cargo test` silently builds
**the author's checkout** rather than the worktree.

## ⚠ Six wrong claims of mine

Every figure quoted from the design documents held — 4.2 million scores, 66 to 91 genotypes, the
22-fold span, `START_AGREEMENT_LIMIT` as one quarter-Phred. Every wrong one was mine.

1. **"495-shape entry space"** — 495 shapes exist at depth 4 over nine buckets, but the fitted table
   is **222 cells** after the rounding drops the rest, and the sentence was about what a climb
   costs.
2. **"forty golden-section steps, eight sweeps"** as a cost — forty is the cap; the bracket reaches
   the tolerance in about 29 steps an axis.
3. **"nine minutes in a debug build"** — measured at **twelve** (719.6 s), and 775.7 s on the
   re-run after these fixes.
4. **"under 1e-9 of the total"** for the rounding — that is what each *dropped shape* carries; 273
   of them drop and together they hold **1.1e-8**, with the rounding over all 495 coming to
   **6.4e-8**. The follow-on "four orders below the sharpest claim" was wrong with it.
5. **"30 loci in 100 carry something other than the reference"** — 30 in 100 is the *allele*
   frequency; under Hardy-Weinberg it is **51 loci in 100**.
6. **"without the warm start a 36-genotype stratum exhausts the pass cap on the first candidate"** —
   the first candidate always climbs from uniform anyway, so the warm start could not have helped
   it; and the cap is still reached with the warm start "on".

## The oracle

**The control the plan names, and it passes.** Fill a stratum's whole entry space with each shape's
exact probability under a known truth, fit, and get the truth back: no sampling noise, so anything
but recovery is this code's fault. It recovers the level to 0.03%, the direction split to 0.00006
and the fall-off to 0.0002, no start beats the score at the truth, and the four starts agree.

**It is `#[ignore]`d and run by hand, which is a departure worth the owner's eye.** It takes 775 s
in a debug build, and the cost cannot be cut without cutting what it says: at a coarser resolution
a genuine 1% bias and the search's own step are the same size, so it would be reporting its own
grid. A second test, `every_start_is_reported_with_its_score_and_the_spread_across_them`, runs in
the suite in about 25 s and pins what the search *reports* — every start recorded, best-scoring
first, the spread measured on the level, and how the search ended.

## ⛦ The follow-up commit, and why it was needed

**The reviewer's mutation table arrived after the step was committed, and it carried the finding
the report file did not: nine of twenty-one mutations of the search passed the whole suite.** Among
them a golden section with its comparison reversed — so the search walks *downhill* — a search that
returns its starting point untouched, and an axis write that puts the direction split into the
fall-off. The sharp control catches all three; nothing that runs caught any of them, so the step
had shipped with its only real test switched off.

`the_search_moves_toward_the_truth_and_says_it_settled` closes it: every locus at the reference
length, two reads deep, starts at 0.09, 0.03, 0.01 and 0.009 against a truth of 0.0201, asserting
the level to 10% and both termination flags. I reproduced the downhill mutant against it — the
level rails at 1e-5 and the test fails.

**What it does not assert is the direction split and the fall-off**, and that is deliberate rather
than a gap: at this fixture every locus is the reference length, so how *often* a read slips is
sharply identified while how *far* it slips rests on a few reads in ten thousand. Asserting those
here would be asserting noise. They are measured by the sharp control, over a wider allele spectrum
and a deeper locus.

**The cost is real and worth stating**: this file's tests went from 27 s to 171 s. Every candidate
the search tries is a whole climb over the stratum's 36 genotypes, and no cheaper fixture both
identifies the level and converges — one at two reads with alleles either side neither settles in
five sweeps nor runs faster.

## Validation

`cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` and
`cargo test --lib --bins --tests --all-features` in the container: **3,439 lib tests, 11 ignored**.
Suite 3,518 → 3,520.
