# ng — calling prerequisites, F1: the slippage lookup the caller borrows

**2026-08-24**, branch `ng-calling-prerequisites`. Step F1 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`arch/read_likelihoods.md`](../../ng/arch/read_likelihoods.md) §4.2 and
[`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2.

**The pre-pass has fitted how each library's polymerase slips; there was no way to ask for one
answer.** This is that way, and it computes nothing.

---

## 1. What changed and why

Scoring a candidate at a short tandem repeat needs three numbers for the library that produced the
read: how often a read reports a tract length other than its allele's, which way the slips go, and
how fast multi-repeat slips fall off. They are fitted per read group within a **stratum** — a motif
length and the reference's repeat count. Every piece existed;
[`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2 records the gap in so many words:
*"`StratumFits` has no wrapper type today — the pieces exist … but nothing gathers them, so whoever
writes the parameter-prepass arch doc names it."*

## 2. Changes made

**[`src/ng/parameter_estimation/joint/stratum_fits.rs`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs)**, new.

- `StratumFits::over(outcomes, slippage_group_of)` — built once per run from what `fit_strata`
  returned. It refuses two outcomes naming one stratum, and an outcome whose three per-group
  vectors are not one length.
- `StratumFits::at(read_group, period, candidate_repeats)` →
  `Result<FittedSlippage, NoSlippage>`. **The two halves of the key are taken by name**, so a
  caller cannot pass the reference tract's repeat count where the candidate's belongs — see §4.
- `FittedSlippage` carries the three numbers, where the level came from, and where the two shares
  came from.
- `NoSlippage` keeps four absences apart, two of them ordinary traffic and two of them saying the
  group map and the fit came from different runs.

**[`doc/devel/ng/arch/parameter_prepass_joint_fit.md`](../../ng/arch/parameter_prepass_joint_fit.md)**
gains §1.7 naming the type, which is what the em-loop architecture asked this document to do, and
**[`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md)** has its open item closed and its
`FrozenParameters` note corrected — it described the field as "named by its contents because no
wrapper type exists yet".

Six tests.

## 3. What this step does *not* do, and why that is the whole of its risk

**The plan says the lookup returns the stratum's `Slippage` "with `level` replaced by
`blend_level`'s value for that cell". It replaces nothing, because the replacement has already
happened** — by one of three routes, which the level's provenance names:

- a stratum fitted on its own tracts has had `blend_level` applied in place, after its period's
  curves were drawn from every stratum's own answer
  ([`ssr_fit.rs`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs)'s
  `smooth_levels_across_repeat_count`);
- a stratum too thin to fit takes its period's curve whole and never goes through `blend_level`;
- a run with curves switched off keeps its cell's own answer.

**A gather that blended again would weigh the curve against a number the curve is already inside.**
`str_slippage_level_curve.md` §5.1 does not name that act — what it forbids in so many words is a
*curve* fitted from blended values, "otherwise each round of smoothing fits a curve to the previous
round's curve, and the cells stop being evidence" — and this is the same circularity one step
downstream. So the step's requirement, "the level read off the fitted curve", is met by *not*
re-deriving, and the module's doc says so where a later editor will read it. §4 has what the
re-derivation would cost, measured.

## 4. What the reviews changed

Four agents, each in its own worktree: does the gather match a real fit, every claim re-measured,
test strength by mutation, and whether this is the shape the consumer needs. **65 mutations between
them.** Three findings changed the code, one changed its signature, and eleven changed its prose.

**The gather does match a real fit, and that was worth checking because every fixture is
hand-built.** One agent ran `fit_strata` over six drawn strata and probed all 21
`(read group, stratum)` pairs: 12 with numbers, 9 absences, **no mismatch in the slippage, the
level provenance or the shares provenance, nothing dropped and nothing duplicated**. It also
confirmed the invariant the lookup rests on — that a slippage number always has a level and a
shares provenance beside it, at the same index, in vectors of one length — on both of the fit's
paths and by measurement over every outcome.

**The signature changed, and this is the finding that mattered most.**
[`read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §4.4 is explicit: *"A read's chance of
slipping is a property of the tract it was copied from, and that is the candidate allele, not the
reference"*, and *"the stutter parameters cannot be hoisted out of the candidate loop"*. The lookup
took a `Stratum`, whose own field is `reference_repeats` — the right word on the fit's side of the
seam and the wrong number on the caller's — and **my doc comment told the caller to fill it from
the tract**, while citing the section that forbids exactly that. A caller that believed it would
score every candidate at a locus against one polymerase model; at tomato dinucleotides that is
about 6 reads in 100 slipping against about 15. `at` now takes `period` and `candidate_repeats` by
name, so the mistake is one somebody has to type on purpose, and a new test pins §4.4's own example
— a 6-repeat and a 12-repeat candidate at one tract are two strata.

**The test named as the guard for this module's one real risk could not fail.** The module exists
to *look up* the emitted slippage level rather than re-derive it, and
`a_level_the_curve_supplied_is_the_curves_own_value_at_that_repeat_count` stored
`curve.level_at(10)` as the stratum's level and then asserted the answer equalled
`curve.level_at(10)`. Blending a number with itself returns that number at any weight, so the
fixture was a fixed point of the very operation it was meant to catch. Two agents built the
re-blending defect and ran it: one saw the test pass outright, the other saw it fail by
`0.11000000000000001` against `0.11` — one unit in the last place, and `exp(ln x) == x` holds for
**69,956 of 100,000** levels drawn between 0.0001 and 0.5, so the kill was a coin flip on the
fixture's constants. The fixture is now a stratum whose level sits *between* its own fit and its
curve and equals neither, with an `assert_ne!` pinning that so it cannot degenerate again. What the
defect is worth where the level is not already the curve's: on a small real fit, re-blending the
five blended strata moved them by **0.6 % to 4.1 %**, and left the `curve_weight` in their
provenance unchanged — the number moves and the provenance does not say so.

**Four mutations survived because every fixture was degenerate**, and all four now die:

- **`shares` was never observed.** The fixture helper hard-coded `shares_provenance: vec![None; …]`
  — a shape neither of the fit's paths can produce — so returning `None` for every answer passed.
- **Per-group level provenance was never checked.** Both groups in the only two-group fixture
  carried identical provenance, so reading group 0's for every group passed. They are now told
  apart by a count rather than a variant, because a variant is what a wrong index most easily
  preserves.
- **Two outcomes naming one stratum silently lost one**, and the doc said this could not arise from
  `fit_strata`. It can: `fit_strata` returns one outcome per *evidence* it was handed, and what
  makes strata distinct is `gather_strata`, which keys them off a map. `fit_strata` is public and
  three examples and a benchmark build evidence by hand. Measured on two outcomes for one stratum:
  levels 0.05 and 0.99, one kept, nothing said. Now a release-level assertion naming the stratum.
- **`slippage_group_of` was asserted only where the answer was group 0**, which is what the
  always-zero mutant returns.

**Two behaviours were merged that should not have been.** A slippage group past the end of a
stratum's rows was reported as a quiet library. It is not: it means the group map names more groups
than the fit was run over, so the map and the fit came from different runs — the same class of fact
as an unknown read group. It has its own answer now, and `NoSlippage`'s four variants are ordered
and labelled by which two are ordinary traffic and which two say the run is not what it claims.

**Eleven prose corrections**, of which four mattered:

- **"one group per read group is the specified default"** — it is the specified *grain*. The only
  builder of that map in the tree pools every read group into one set unless told otherwise, so the
  sentence described the opposite of what every run here actually does.
- **"the circularity §5.1 forbids in so many words"** — §5.1 forbids a *curve* fitted from blended
  values. Re-blending an emitted level is the same circularity one step downstream and the spec
  does not name it. Corrected, with §5.1 quoted so a reader can see the difference.
- **"`fit_strata` applies `blend_level` in place"** — true on the fitted path only. A stratum too
  thin to fit takes its curve whole and never goes through `blend_level`; a run with curves off
  keeps its cell's own answer. The conclusion holds on all three, the mechanism on one.
- **The `ContaminationEstimate` analogy** — that type is an enum to separate *absence from a
  value*, which a `Result` already does. The structural twin is `NotIdentifiedReason`: several
  named reasons behind one absence.

Also corrected: "answers every accessor with an empty slice" (three of five), `strata()`'s claim to
distinguish a cohort with no tracts from a gather nobody built (it cannot, and both answers mean
the same thing to a caller), and a `PartialEq` comparison of two `Result`s in the pooled-groups test
that passed under a mutation making every lookup fail — two absences compare equal too.

**One check moved from the lookup to the build.** The three per-group vectors are the same length by
construction, but every field of the outcome types is public, so a hand-built outcome could be
short — and the failure would then have been an index panic inside a lookup rather than a sentence
naming the stratum. It is asserted in `over` now.

### What the reviews found for whoever writes the consumer

- **A candidate whose repeat count no kept reference tract occupies gets `NoSuchStratum`, and that
  is ordinary traffic rather than an error.** §4.4's own example is a 12-repeat candidate at a
  6-repeat tract; if no kept tract of that period has 12 repeats, there is no row. The curves that
  could answer are thrown away — `over` keeps the per-group cells, though every stratum that used a
  curve carries a copy of it in its provenance, and `str_slippage_level_curve.md` §6 designs a
  defined answer beyond the fitted range. Whether a repeat count with no tract behind it should be
  furnished at all is the spec owner's call, and it should be taken before a caller is built on
  this.
- **The fit and the lookup disagree about an unknown read group.** `gather_strata` folds one into
  slippage group 0 without a word (`unwrap_or(&0)`); `at` refuses it. So a library's reads can be
  fitted into group 0's numbers and then denied them. The map is total on one run through the walk,
  so this is dead today — which is the moment to make the two agree.
- **The STR substitution rate sits at the same `(read group, stratum)` grain** (spec §6.1) and is
  not carried. It lives on `StratumEvidence`, which is gone by the time `over` is called, so adding
  it later changes this signature.
- **Cost, measured**: two `BTreeMap` probes cost 22 ns at 300 strata and 63 read groups, against
  1.8 ns for flat vectors — about 15 s a run at the tomato cohort's ceiling against the 155 s the
  tract fit already takes. **Leave the map.** Hoisting the read-group probe out of the candidate
  loop, where it does not belong, is free and buys 13 of the 22 ns; that is the caller's to do.

## 5. Validation

All in the dev container, on the tree as committed.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::parameter_estimation::joint::stratum_fits` | **6 passed, 0 failed** |
| `cargo test --lib ng::parameter_estimation` | **733 passed, 0 failed, 6 ignored** |
| `cargo test --lib` | **4,174 passed, 0 failed, 14 ignored**, 581.84 s |
| `cargo doc --no-deps` | 24 unresolved-link errors, 12 redundant-explicit-link-target warnings — the recorded baseline, unchanged |

**The reversion test says yes and here it says nothing**, which all four agents reported the same
way: removing the module's one line from the module list takes the type and its only tests together,
because nothing else in the tree names it. That is the shape of a type built for a seam that does
not exist yet, not a defect signal. What the mutation runs measured instead is that **every one of
the six tests has a mutation only it catches**.

## 6. Follow-ups

- **A candidate repeat count no kept reference tract occupies has no answer**, and the per-period
  curves that could give one are not gathered (§4). Owed to the read-likelihoods plan, or to the
  spec's owner if the answer is that such a candidate should get nothing.
- **`gather_strata` folds an unnamed read group into slippage group 0 and `at` refuses it** (§4).
  Dead today; the two should be made to agree before it is not.
- **The STR substitution rate is at this grain and is not here** (§4). Adding it changes `over`'s
  signature, because the evidence it lives on is gone by then.
- **Nothing consumes any of this yet.** The consumer is `FrozenParameters`
  ([`calling_loop.md`](../../ng/impl_plan/calling_loop.md)) and the STR read likelihood
  ([`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md)).
