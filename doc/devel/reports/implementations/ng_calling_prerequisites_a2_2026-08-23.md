# ng — calling prerequisites, A2: the fitted inbreeding coefficient is capped at 0.99

**2026-08-23**, branch `ng-calling-prerequisites`. Step A2 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §7 and
[`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §2.1, third bullet.

**One constant, one helper, five tests.** [A1](ng_calling_prerequisites_a1_2026-08-23.md) made
`InbreedingF` half-open `[0, 1)`; this is the change that keeps a legitimate fit from becoming a
panic on the way in — and, as the review showed, the change that matters on far more samples than
the ceiling ever would.

---

## 1. What changed and why

`fit_inbreeding` ends by handing its answer to `InbreedingF::try_new(…).expect(…)`. The answer is
the **coverage-weighted posterior occupancy** — each 100 kb window's forward–backward probability
of lying inside a run of homozygosity, weighted by the reference positions it covered, summed and
normalised (`spec/parameter_prepass_generic.md` §6.5). **Two things follow, and the second is the
larger.** A fit that reaches exactly `1.0` now panics, because A1 made that value unconstructible
and a caller cannot recover from a panic. And a fit anywhere above `0.99` — far commoner than the
ceiling, and reached by any long-selfed line — hands the prior a heterozygote branch no read
evidence can move. Neither is an error on the sample's part.

**The cap is `0.99`, and the number is not arbitrary.** The genotype prior multiplies its
heterozygote branch by `1 − F`. At `0.99` that leaves 20 on the Phred scale against every
heterozygote, which two clean alternative bases at Q30 — 60 Phred — can overturn. At the greatest
value `InbreedingF` accepts, `1 − F` is `2⁻⁵³`, about 160 Phred, and no read evidence could. So
the newtype's range removes the mathematical limit and nothing more; **capping the estimate is
this estimator's own job**, which is exactly what spec §7 says and why A1 alone was not enough.
`0.99` is production's own ceiling, `MAX_INBREEDING_COEFFICIENT`
([`inbreeding.rs:25`](../../../../src/paralog/inbreeding.rs)) — copied, with its reasoning, rather
than depended on: ng does not import from the frozen tree.

## 2. Changes made

**[`src/ng/parameter_estimation/generic/runs.rs`](../../../../src/ng/parameter_estimation/generic/runs.rs)**

- `MAX_FITTED_INBREEDING = 0.99`, private to the module, documented against both the Phred cost of
  the cap and the Phred cost of not capping.
- `capped_inbreeding(occupancy: f64) -> InbreedingF` — clamps, then constructs. The one call site
  in `fit_inbreeding` goes through it.
- `NaN` is deliberately not capped. `f64::clamp` propagates it, the constructor rejects it, and the
  `expect` fires — which is the case that convention exists for, our own arithmetic being broken.
  **`f64::min` would have been the careless choice and is silently wrong here:** it ignores `NaN`
  by definition, so a broken fit would come back as `0.99` — a confident, plausible number nothing
  downstream could question. Measured, not recalled: `NaN.min(0.99)` is `0.99` and
  `NaN.clamp(0.0, 0.99)` is `NaN`.
- Five tests. Four on the helper: the cap holds at and above the newtype's ceiling and at the
  lower bound, with the two Phred figures asserted rather than stated and `0.99` asserted against
  the literal; every value below the cap passes through bit for bit, `0.0` included; and two on
  the `NaN` route, one pinning both halves of the `min`-against-`clamp` trap and one pinning that
  the constructor is where a broken fit stops.
- **The fifth goes through `fit_inbreeding`, and it is the one that matters.** All three reviews
  found the same hole independently: with the call site restored to constructing directly — the
  whole point of the step deleted from the production path — every one of the library's 4,144
  tests stayed green.

## 3. What the reviews changed, and the one thing they changed about the design's own story

Three agents, each in its own worktree: what a fitted result now silently differs by; every number
and mechanism re-measured; and the strength of the tests under mutation. **31 mutations between
them, of which 8 survived the tests as first written.** Every survivor is closed.

**The step had no test of the thing it does.** Reverting the call site left the whole library
green, in all three reviews independently. Closed by
`the_emitted_coefficient_is_capped_where_the_fit_itself_went_higher`, which builds its genome
window by window rather than drawing one: the drawn-genome helper walks a chain whose entry rate
is `LEAVE_RUN·F/(1 − F)`, so past about `F` = 0.97 that rate exceeds one and the draw saturates —
measured, a nominal 0.995 and a nominal 0.9995 both return the same fit of 0.966.

**And the mechanism this report first gave for why the cap is needed was wrong.** A1's report and
this one both said a sample whose every window lies inside a run returns an occupancy of exactly
`1.0`. Two separate reviews built that sample and showed it is *refused* — the fit needs some
window's posterior below one half to see a second state at all, so it never reaches the
constructor. **The case that matters needs no ceiling:** a genome 3,599 windows of 3,600 inside a
run is accepted, fits above `0.99`, and lands where read evidence can no longer move the prior.
That is ordinary rather than exotic — under self-fertilisation `F = 1 − 2⁻ᵗ` passes `0.99` at the
seventh generation, so a maintained inbred line is past it and the tomato benchmark is 63 of them.
The doc comments and A1's report now say this; the old wording sent a reader hunting a symptom
that does not occur.

**Two numbers were wrong and are corrected.** "About 37 nats" is the textbook threshold at which
`1 + exp(−d)` stops rounding to one; here the normaliser is added to a chain log-likelihood in the
thousands or billions, which brings it to about 27 nats on this module's own fixture and about 16
at the totals this file quotes for a real genome. And "leaves that branch at `ln 0.01`" named the
wrong quantity: `ln 0.01` is −4.6, while the 20 Phred is `0.01` itself.

**Three smaller corrections**: the doc re-imported the guarantee spec §7 had withdrawn the same
day (production's clamp covers a *fitted* `F` only); "imported with its reasoning" read as a code
dependency that does not exist, ng not importing from the frozen tree; and the link to production's
constant used a filesystem path, which rustdoc emits verbatim into a page where it points at
nothing.

**One thing deliberately not changed.** The fit's own record, `RunsModelFit::starts_tried`, stays
uncapped. It is the only surviving trace that the cap fired, and `MAX_IDENTIFIED_START_SPREAD` is
measured across those values — nine starts between 0.94 and 1.00 spread by 0.06 and are refused,
where capping them to 0.94 and 0.99 spreads them by 0.05 and they are not. Said in the doc comment
now, since a reader comparing the two numbers otherwise finds a mismatch with nothing explaining it.

## 4. Validation

<!-- filled from the container run -->

## 5. Follow-ups

- **A capped sample and a sample genuinely fitted at `0.99` are the same number to every
  consumer, and nothing records which it was.** `Estimate::provenance` is still `FittedHere` and
  `RunsModelFit` has no flag. On a crop panel of inbred lines this is not the rare case — it is
  most of the panel — so a reader cannot tell a cohort whose coefficients were fitted from one
  whose coefficients were flattened onto the ceiling. **Recommended: a `bool` on `RunsModelFit`
  beside `undecided_windows`**, not a `Provenance` variant, since the value really was fitted
  here. Raised at Checkpoint A rather than taken, because it adds a field to a type the plan did
  not put in scope.

- **The joint route's second estimator will hit this too.** `HomozygoteExcess`
  ([`joint/fit.rs:157`](../../../../src/ng/parameter_estimation/joint/fit.rs)) accepts the closed
  `[0, 1]`, so exactly `1.0` constructs. Nothing converts it to an `InbreedingF` today; whatever
  does owes a cap of its own, which is the reason this one lives in the estimator rather than in
  the newtype.

- **A supplied coefficient would go round this ceiling**, and the newtype alone would not stop it
  — `0.999` constructs. When ng gets a command line, its parser owes this cap and not merely the
  type's range. Recorded in `spec/calling_priors.md` §7.

- **The clamp's lower bound is unreachable defence.** `StateFit.inbreeding` is already clamped to
  `[0, 1]` where it is built, so nothing negative reaches the helper. It is asserted anyway, since
  the cheapest wrong edit — `min` in place of `clamp` — changes exactly that end and the `NaN`
  behaviour with it.
