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
design means. It is also the only workable choice — a cell at depth 124 has a likelihood no `f64`
holds in linear space. The parameter is named `ln_likelihood_by_cell_and_genotype` rather than
`component_likelihoods` so that a caller cannot hand it linear probabilities and get a plausible
wrong answer; `−∞` is accepted and means this genotype cannot have produced this cell, and it is
the only non-finite value accepted.

**Convergence is reported to the crate, not to the consumer.** The plan asks for convergence to
be asserted in tests rather than propagated as a flag no consumer reads. `MixtureWeightsFit` is
`pub(crate)` and carries `converged` and `passes`; the public function returns the weights alone.
`log_likelihood` is on the same struct because the profile scan of D3 needs a score to compare
rungs on and would otherwise walk the table a second time — it is computed in a final pass after
the climb stopped, so it belongs to the weights returned beside it and not to the pass before them.

**The cap is generous and exhausting it is not an error.** Expectation-maximization converges
linearly at a rate set by how much the components overlap, and they overlap most at low coverage —
this estimator's own regime. A cap that fired on a slow-but-correct climb would turn shallow data
into a crash. This is not hypothetical: the tetraploid fixture below was first written with five
dosage columns that overlapped heavily, and it came back 0.6977 against a truth of 0.7000 after
1,000 passes without converging. The fixture was replaced with columns separated the way a real
dosage table's are, and that test now asserts `converged` so the same thing cannot pass quietly.

## Deviations from the architecture, recorded

1. **Parameter renamed** — `component_likelihoods` → `ln_likelihood_by_cell_and_genotype`, as
   above. Same shape, `&[&[f64]]`.
2. **A second entry point**, `climb_mixture_weights`, crate-internal. The architecture sketches one
   function; the scan needs the score and the tests need the start.
3. **Panics on a malformed table.** The architecture's signature has no error type and none is
   added: a ragged table, a `NaN` likelihood, a negative cell weight or a cell no genotype could
   have produced are all faults in the code that built the table, and each of them otherwise
   leaves a `NaN` to travel through the fit as a plausible number. The messages name the cell.

## What the tests are, and what makes each able to fail

Twenty new tests, taking `ng::parameter_estimation` from 130 to 150. The recovery fixture is the same device the research harnesses use: each
genotype's column is a distribution over the cells, so weighting each cell by
`Σ_j π_j · L(cell | j)` makes `π` the maximiser **exactly**, with no sampling noise in it. That is
the infinite-genome table, and the assertion is `1e-9`.

- `the_climb_recovers_the_frequencies_that_generated_the_table` — the claim the step exists for.
  The fixture is **4 cells × 3 genotypes**, deliberately unequal, so an index transposed between
  cell and genotype is a length mismatch rather than a wrong number.
- `every_interior_start_reaches_the_same_summit` — four starts, three of them near a corner. This
  is the concavity of `spec/parameter_prepass.md` §3.1 asserted rather than assumed, and it also
  asserts `converged` and `passes < MAX_CLIMB_PASSES` at each.
- `two_truths_over_one_likelihood_table_give_two_answers` — the same likelihood table under two
  truths. What this makes reachable: a climb that ignored `cell_weights` entirely would still
  return a point on the simplex and would still pass a single-truth recovery test.
- `no_move_away_from_the_fitted_weights_scores_higher` — eighteen steps along the simplex's edges
  from the fitted point, each scored. This is the only test that checks the answer against the
  objective instead of against the construction that built the fixture, so it survives a change of
  fixture.
- `the_reported_score_belongs_to_the_weights_returned_with_it` — the score field D3 will consume,
  checked against `Σ_cells w · ln Σ_genotypes π·L` written out longhand in linear space, not
  against the climb's own running total.
- `a_tetraploid_table_is_fitted_over_its_own_five_genotypes` — five genotypes over six cells.
  Nothing here is diploid, and five is past the three a `SmallVec<[f64; 3]>` holds inline, so this
  exercises the spill as well as the loop bound.
- `the_fitted_weights_are_a_point_on_the_simplex` — a maximization step that forgot to divide by
  the total weight would return counts, and every other assertion above is scale-free.
- `the_log_sum_is_exact_where_the_linear_sum_still_works_and_survives_where_it_does_not` — the
  shift in `ln_sum_exp` checked against shift-free arithmetic where that arithmetic works, then
  against three terms at −800, −810 and −820, which a linear sum flushes to zero.
- Ten refusal tests: eight on the table — length mismatch, ragged, `NaN`, `+∞`, negative weight,
  every weight zero, empty table, a weighted cell no genotype can produce — and two on the start
  (on a face, and not summing to one). Two tests pin the pair that must *not* be refused — a `−∞` entry beside a
  finite one, and an all-`−∞` cell that carries no weight.

A `debug_assert!` inside the loop holds the pass to monotone ascent, which is the one thing
expectation-maximization is actually proved to give.

## Validation

All in the container. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → **3,052 passed**, 0 failed, 5
ignored, up from 3,032; `cargo doc --no-deps --lib` at the 12-unresolved-link pre-existing
baseline, unchanged.

`examples/ng_multilib_key_harness.rs` was run green before the step, as Milestone D's precondition
requires.

## Open, carried forward

- `MAX_CLIMB_PASSES` is untested against a real table, because no real table exists until F3. If
  the generic path's 583-cell table at three reads a site does not settle inside 1,000 passes, that
  will show up as a fit that moves when the cap is raised — worth checking once F2 fills a
  histogram from known parameters.
