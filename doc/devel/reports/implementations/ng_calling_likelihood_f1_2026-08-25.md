# ng read likelihoods — F1: the emission seam

*Implementation report, 2026-08-25. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step F1 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone F, on
top of `7f834dfd`.*

## 1. What it is

Three types in the new [`calling/likelihood/ssr_emission.rs`](../../../../src/ng/calling/likelihood/ssr_emission.rs),
and no arithmetic at all: **the seam**. What an STR emission is handed, what it returns, and
nothing about how it decides.

| type | what it is |
|---|---|
| `SsrCandidate` | one candidate allele: its bases, and the repeat count that keys the stratum lookup |
| `SsrScoringContext` | everything a model is handed for one `(read group, candidate)` — and the only channel it has |
| `SsrEmissionModel` | the trait, with `emission` and `censored_emission` |

## 2. The three decisions the types make binding

**Every number arrives per call.** Nothing is read from global state, so the EM loop can
re-estimate the slippage numbers between iterations with no change on this side — spec §6.1's one
binding constraint on this module.

**A context is per `(read group, candidate)`, not per locus.** A read's chance of slipping is a
property of the tract it was copied from, and that is the *candidate*: 6 repeats and 12 at the same
locus are different strata slipping at measurably different rates, about 1.3-fold per repeat count
(spec §4.4). `two_candidates_at_one_locus_get_different_contexts` pins it — the four-repeat
candidate leaves more mass unplaced than the thirty-repeat one, from the same locus and the same
model.

**A repeat count is `NonZeroU32` here too**, matching the contract E2's review put on
`StutterModel::unreachable_mass`: a candidate whose tract holds no repeats is not a candidate.

## 3. Two things the seam takes rather than trusts

**The unreachable mass is read off the distribution, not passed in.** The architecture sketches it
as a plain field; taking it from a caller would be two ways in and two chances to disagree, on the
one number that keeps candidates comparable. `SsrScoringContext::new` computes it from the model it
was handed.

**The weakest warrant is folded, not asserted.** Spec §4.4 requires the weakest provenance of any
parameter entering a locus to travel onto its output. That needed an ordering, and none existed —
`Provenance` had no way to combine two. `Provenance::weaker_of` now supplies it, ordered by **the
ladder this repository already states** in `ParameterEstimationError`'s own doc: *fitted here,
borrowed from the sample's other read groups, supplied, defaulted*. So a supplied value ranks below
a borrowed one — a number the run was handed says nothing about this data, where a borrowed fit is
at least a measurement of a neighbouring grain. The placement is the ladder's, not this step's, and
the doc says where to change it if the ladder is wrong.

The fold's identity is `FittedHere`, so a context nothing weakened stays fitted;
`a_context_with_no_weakening_warrant_is_fitted` pins that, because the alternative would mark every
context of a fully-fitted run as a guess.

## 4. Deviations

- **`Provenance::weaker_of` is new, in `parameter_estimation/mod.rs`** — a shared module the plan
  does not name for this step. It is additive, and F1 is the first consumer that needs to combine
  two warrants at all.
- **A private `period_of` converts a `Motif` to the `NonZeroU8` the distribution wants.** ng's
  `Motif::period()` returns `usize`; changing it would ripple into frozen production, which shares
  the method name on its own type. One conversion, in one place, citing the invariant that makes it
  infallible.
- **The trait's methods take `&[u8]` for the observation**, as the architecture sketches. No
  evidence view is threaded yet — that is H1's, which routes complete and partial observations to
  the two methods by witness.

## 5. Validation

| command | result |
|---|---|
| `./scripts/dev.sh cargo test` — library target | **4,383 passed, 0 failed, 14 ignored** (4,379 at Milestone E's close) |
| `clippy --lib --all-features --tests -- -D warnings` | exit 0 |
| `cargo fmt --check` | exit 0 |

Four tests, all on the context's construction — there is no arithmetic here to test, which is the
point of a seam.

## 6. Deferred for the owner

Carried forward, plus one new:

1. **Spec §4.2's prose and its figures disagree** about whether a tract can lose its last repeat
   (E2's report §3). This branch follows the figures, and Milestone E's review re-derived both
   readings independently and confirmed the argument.
2. **Spec §4.2 calls the part-repeat cutoff "10 base pairs"** where it is applied to a compressed
   rank; ten of those admit about 13 base pairs at period 4.
3. **`arch/read_likelihoods.md` §4.2 sketches `stutter_rates_for` beside the distribution**, the
   placement E3 rejected on module-edge grounds and the review endorsed rejecting. It also names
   the context's field `truncated_mass_lost`, which is now `unreachable_mass`.
4. **Spec §12's fourth test asks periods 1 to 6; the tripwire runs 2 to 6**, because at period 1
   the part-repeat branch is unreachable by construction.
5. **`Provenance`'s ladder puts *supplied* below *borrowed***, and `weaker_of` follows it. If a
   run's supplied values should outrank a borrowed fit, the ladder is the place to say so.
