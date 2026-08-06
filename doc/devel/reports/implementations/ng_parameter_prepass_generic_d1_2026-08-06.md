# ng step 4, D1 — `fit_mixture_weights`, the concave climb

**Date:** 2026-08-06. **Branch:** `ng-parameter-estimation`. **Plan:**
`doc/devel/ng/impl_plan/parameter_prepass_generic.md`, Milestone D step D1. **Design:**
`doc/devel/ng/arch/parameter_prepass_generic.md` §4.1, `doc/devel/ng/spec/parameter_prepass.md`
§3.1.

## What landed

`src/ng/parameter_estimation/fitting/mixture_weights.rs`, which was documentation before this
step. One public function and one crate-internal one:

- `fit_mixture_weights(ln_likelihood_by_cell_and_genotype, cell_weights) -> SmallVec<[f64; 3]>` —
  the architecture's signature. Given how likely each genotype makes each cell, and how much each
  cell counts for, return the genotype frequencies that best explain the table.
- `climb_mixture_weights(table, cell_weights, start) -> MixtureWeightsFit` — the same climb with
  the start named and `{ weights, log_likelihood, passes, converged }` reported.

Expectation-maximization, from the uniform point by default: an expectation step forming each
cell's per-genotype responsibility in logs, a maximization step averaging those responsibilities
under the cell weights. It stops when no weight moved by more than `CLIMB_STILLNESS = 1e-13` in a
pass, or after `MAX_CLIMB_PASSES = 1_000`.

## The three decisions the architecture left open

**The likelihoods are natural logarithms.** The architecture writes the parameter as
`component_likelihoods` without saying which space it is in; §5.1 writes the scoring rule as
`ln L(cell | θ)` and the research harness's `climb_frequencies` takes logs, so logs is what the
design means. It is also the only workable choice — the corners of a deep table are outside what
`f64` holds in linear space: a site at the depth cap of 124 whose every read showed the alternative
allele is `124 · ln ε` under the homozygous-reference genotype, which is `−857` at an error rate of
0.001 and `−1428` at the ladder's floor, both zero in linear space, while the ordinary 124-read cell
beside them (no alternative read at all) is 0.883. At the ladder's floor a depth of 65 already
underflows. The parameter is named `ln_likelihood_by_cell_and_genotype` rather than
`component_likelihoods` so that a caller cannot hand it linear probabilities and get a plausible
wrong answer; `−∞` is accepted and means this genotype cannot have produced this cell, and it is
the only non-finite value accepted.

**Convergence is reported to the crate, not to the consumer.** The plan asks for convergence to
be asserted in tests rather than propagated as a flag no consumer reads. `MixtureWeightsFit` is
`pub(crate)` and carries `converged` and `passes`; the public function returns the weights alone.
`log_likelihood` is on the same struct because the profile scan of D3 needs a score to compare
rungs on and would otherwise walk the table a second time — it is computed in a final pass after
the climb stopped, so it belongs to the weights returned beside it and not to the pass before them.

**Concavity rules out a false summit, not a slow one, and exhausting the cap is a data condition
rather than a bug.** `spec/parameter_prepass.md` §3.1 proves the surface has no local maximum that
is not also global, and in the same paragraph says the rate of approach is only linear and slowest
where the components overlap — which is low coverage, this estimator's own regime. So a climb that
exhausts the cap has not found a second summit; it has run out of time on the only one, and a cap
that fired on a slow-but-correct climb would turn shallow data into a crash.

**The cap is measured, and D1's first value was too small.** Reaching `CLIMB_STILLNESS` on the
four-cell fixture takes 257 passes at a truth of `[0.60, 0.35, 0.05]`, 449 at `[0.90, 0.07, 0.03]`
and **1,234** at `[0.80, 0.02, 0.18]`. D1 shipped at 1,000, which cut the third one off — it
returned `converged = false` 3.6 × 10⁻¹⁰ short, inside the 10⁻⁹ its test asserts, and the test could
not tell because it went through `fit_mixture_weights`, which discards the flag. **This was found by
review, not by the author**, and the sentence the author wrote about it — "the same thing cannot
pass quietly again" — was false one test away from where it was written. The cap is now 10,000,
eight times the slowest measured; every recovery test goes through `climb_mixture_weights` and
asserts `converged`; and `climb_with_cap` takes the cap as an argument, so "does the answer move
when the cap is raised?" is a test rather than a recompile.

## Deviations from the architecture, recorded

1. **Parameter renamed and retyped** — `component_likelihoods: &[&[f64]]` becomes
   `ln_likelihood_by_cell_and_genotype: GenotypeLikelihoodTable<'_>`, one flat row-major buffer with
   the width named once. The rename is so a caller cannot hand it linear probabilities; the retype
   is because the row-of-slices shape **cannot be handed a reusable buffer**. A `Vec<&[f64]>`
   collected before the ladder borrows the buffer, so refilling it per rung is
   `error[E0502]`, and D3 would have to rebuild ~583 fat pointers inside the rung loop, 161 times
   per read group. The new type takes one slice and one `usize`, so there is no pair of same-typed
   values to transpose, and a raggedly-built table becomes a length the width does not divide —
   unrepresentable rather than checked on every call.
2. **A second entry point**, `climb_mixture_weights` (`pub(super)`), plus a private
   `climb_with_cap`. The architecture sketches one function; the scan needs the score, and the tests
   need the start and the cap.
3. **Panics on a malformed table.** The architecture's signature has no error type and none is
   added: a ragged table, a `NaN` likelihood, a negative cell weight or a cell no genotype could
   have produced are all faults in the code that built the table, and each of them otherwise
   leaves a `NaN` to travel through the fit as a plausible number. The messages name the cell.

**One thing to carry into D3:** `fit_mixture_weights` is `pub` only because that is what keeps the
module from being dead code until a consumer exists. When D3's scan lands it should become
`pub(crate)`, which is the honest visibility for a function whose every consumer is in this crate.

## What the tests are, and what makes each able to fail

Twenty-eight tests after review, taking `ng::parameter_estimation` from 130 to 158. The recovery
fixture is the same device the research harnesses use
(`doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md` §1): each genotype's column
is a distribution over the cells — every column verified to sum to 1.0 exactly — so weighting each
cell by `Σ_j π_j · L(cell | j)` makes `π` the maximiser **exactly**, with no sampling noise in it.
That is the infinite-genome table. The diploid recoveries assert `1e-9` and the tetraploid one
`1e-8`.

- `the_climb_recovers_the_frequencies_that_generated_the_table` — the claim the step exists for.
  The fixture is **4 cells × 3 genotypes**, deliberately unequal, so a fixture built with the two
  dimensions swapped is a buffer the width does not divide rather than a wrong number.
- `every_interior_start_reaches_the_same_summit` — five starts: one in each of the three corners'
  neighbourhoods, the uniform point, and one ordinary interior point. This is the concavity of
  `spec/parameter_prepass.md` §3.1 asserted rather than assumed, and it also asserts `converged` at
  each. The heterozygous corner is in the list because it is the coordinate the estimator exists to
  measure; the committed version of this test missed it.
- `two_truths_over_one_likelihood_table_give_two_answers` — the same likelihood table under two
  truths. What this makes reachable: a climb that ignored `cell_weights` entirely would still
  return a point on the simplex and would still pass a single-truth recovery test. The weights span
  11.5 to 1, which is what makes ignoring them fail five tests rather than drift. **This is the
  fixture that measured the cap**, and it now goes through `climb_mixture_weights` and asserts
  `converged`.
- `a_climb_stopped_early_says_so_and_scores_the_weights_it_returns` — capped at five passes, which
  is the only way to reach the two states every other fixture is built to avoid: `converged =
  false`, and a score that belongs to the weights returned rather than to the pass before them. On
  a settled climb those two scores agree to 2.3 × 10⁻¹³, so nothing else can tell them apart.
- `raising_the_pass_cap_does_not_move_the_answer` — the same climb at the cap and at ten times it,
  bit for bit. The cap is a stopping rule, not a shaping one, and this is the experiment that says
  so.
- `a_genotype_impossible_everywhere_does_not_end_the_climb_early` — a genotype that no cell can
  have produced reaches weight zero on the first pass and never moves again. It is the only table
  here where the largest move over all genotypes and the last genotype's move are different
  numbers, and the only one where a weight reaches exactly zero *during* the climb, which is the
  state `genotype_weight.ln()` has to survive.
- `no_move_away_from_the_fitted_weights_scores_higher` — eighteen steps along the simplex's edges
  from the fitted point, each scored, the smallest landing 9.13 × 10⁻³ below the summit against a
  floating-point noise floor of 1.16 × 10⁻¹⁰. It checks the answer against the objective rather
  than the fixture's construction, so it survives a change of fixture — but not the scorer, since
  both points go through the same one.
- `the_reported_score_belongs_to_the_weights_returned_with_it` — the score field D3 will consume,
  checked against `Σ_cells w · ln Σ_genotypes π·L` written out longhand in linear space. This is
  the half `no_move_away_…` cannot hold.
- `a_tetraploid_table_is_fitted_over_its_own_five_genotypes` — five genotypes over six cells.
  Nothing here is diploid, and five is past the three a `SmallVec<[f64; 3]>` holds inline, so this
  exercises the spill as well as the loop bound.
- `the_fitted_weights_are_a_point_on_the_simplex` — a maximization step that forgot to divide by
  the total weight would return counts, and every other assertion above is scale-free.
- `a_rung_loop_refills_one_buffer_and_allocates_nothing_per_rung` — the shape D3 will use,
  compiled: one buffer allocated before the ladder, refilled at every rung, borrowed by the table
  inside the loop. The committed `&[&[f64]]` cannot do this.
- `the_log_sum_is_exact_where_the_linear_sum_still_works_and_survives_where_it_does_not` — the
  shift in `ln_sum_exp` checked against shift-free arithmetic where that arithmetic works, then
  against three terms at −800, −810 and −820, which a linear sum flushes to zero.
- Fourteen refusal tests: ten on the table and its weights — length mismatch, a length the width
  does not divide, `NaN`, `+∞`, negative weight, `+∞` weight, every weight zero, empty table, no
  genotypes at all, and a weighted cell no genotype can produce — and four on the start (on a face,
  not summing to one, the wrong width, `+∞`). Two tests pin the pair that must *not* be refused: a
  `−∞` entry beside a finite one, and an all-`−∞` cell that carries no weight, whose score has to
  come back finite because `0 · −∞` is `NaN`.

A `debug_assert!` inside the loop holds the pass to monotone ascent, which is the one thing
expectation-maximization is actually proved to give.

## Validation

All in the container. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → the library binary at **3,060
passed**, 0 failed, 5 ignored, up from 3,032, with the nine other binaries unchanged at 69 passing
(the command prints no grand total); `cargo test --doc ng::parameter_estimation` → 1 passed, which
the feature's usual validation set does not cover and which the new `# Examples` block needs;
`cargo doc --no-deps --lib` at the 12-unresolved-link pre-existing baseline, none of the twelve in
`parameter_estimation`.

Both research harnesses were run green around this step, as Milestone D's and E's preconditions
require: `examples/ng_multilib_key_harness.rs` before D1 started and
`examples/ng_inbreeding_harness.rs` alongside it, both to completion.

## Open, carried forward

- **`MAX_CLIMB_PASSES` is untested against a real table**, because no real table exists until F2.
  The generic path's is 583 cells at three reads a site, which is exactly the regime
  `spec/parameter_prepass.md` §3.1 names as slowest, and the four-cell fixtures here already span
  257 to 1,234 passes. `climb_with_cap` is what makes that a measurement rather than a recompile:
  fit at 10,000 and at 100,000 and see whether the answer moves.
- **`arch/parameter_prepass_generic.md` §5.2 says the inner climb "needs none" — no cap.** It has
  one, it is measured, and one of this file's own fixtures needs 1,234 passes. The code and
  `fitting/mod.rs` now say so; the architecture does not. Left for the owner.
- **`fit_mixture_weights` should become `pub(crate)` in D3**, when a consumer exists and `pub` is no
  longer what keeps the module from being dead code.
