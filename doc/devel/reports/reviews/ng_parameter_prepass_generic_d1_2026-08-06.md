# ng step 4, D1 — review of `fit_mixture_weights`

**Date:** 2026-08-06. **Reviewed:** commit `b2c1a918`. **Fixes:** the commit after it.
**Agents:** three, each in its own worktree detached at `b2c1a918`, covering ten categories.

| agent | categories | method |
|---|---|---|
| reliability | `reliability`, `errors`, `refactor_safety` | 17 mutations, one at a time, each re-run |
| structure | `module_structure`, `naming`, `idiomatic`, `smells` | every finding applied and rebuilt before it was filed |
| numbers | `defaults`, `tooling`, `extras` | every quantitative claim checked against the research note, the design docs, and the code |

## Verdict

**Three Blockers, five Majors, and every one of the Blockers was a missing test.** That is the
sixth round in eight on this plan where the Blocker was an absent assertion rather than wrong code.
Seven of the reliability agent's seventeen mutations survived the committed suite; all seven now
fail, and each was re-run against the fix to prove it.

The two findings that would have travelled furthest are the ones this section is for.

**⛦ A committed fixture was riding an unconverged fit, and the report claimed the opposite.**
`two_truths_over_one_likelihood_table_give_two_answers` fits the truth `[0.80, 0.02, 0.18]` over the
four-cell table. It runs every one of the 1,000 allowed passes, returns `converged = false`, and
lands 3.6 × 10⁻¹⁰ from the answer against a tolerance of 10⁻⁹ — a margin of 2.8×. It passes because
it went through `fit_mixture_weights`, which discards the flag. The commit message and the
implementation report both closed the pass-cap argument with *"it now asserts `converged`, so the
same thing cannot pass quietly again"*, and that sentence was false one test away from where it was
written. Measured directly afterwards: reaching `CLIMB_STILLNESS` on that fixture takes **1,234
passes**, against 449 for `[0.90, 0.07, 0.03]` and 257 for `[0.60, 0.35, 0.05]`. The cap is now
10,000 — eight times the slowest measured — every recovery test goes through the reporting entry
point and asserts `converged`, and `climb_with_cap` takes the cap as an argument so *"does the
answer move when the cap is raised?"* is a test rather than a recompile.

**⛦ The stillness test can be narrowed to one genotype and nothing notices.** Replacing
`largest_move.max(…)` with the last genotype's move alone left the whole suite green. Every fixture
in the file settles all three coordinates together, so the maximum and the last one are the same
number — but a genotype that no cell can have produced reaches weight zero on the first pass and
never moves again, which is an input the doc comment explicitly calls legal. On such a table the
mutated climb reports **converged after 2 passes** with the first weight at **0.7698** against a
truth of **0.95**, while the correct one takes 241 passes and lands on 0.94999999999. That is wrong
genotype frequencies handed to the consumer as a settled answer, through the entry point that
returns the weights and nothing else. `a_genotype_impossible_everywhere_does_not_end_the_climb_early`
closes it, and it is also the only table in the file where a weight reaches exactly zero *during*
the climb — the state `genotype_weight.ln()` has to survive.

**⛦ The seam could not be handed a reusable buffer, and the compiler is the evidence.** The
architecture writes the argument as `component_likelihoods: &[&[f64]]` (§4.1). The structure agent
compiled the rung loop D3 will write against that shape and got

```
error[E0502]: cannot borrow `flat` as mutable because it is also borrowed as immutable
```

— a `Vec<&[f64]>` collected before the ladder borrows the buffer, so refilling it per rung is a
borrow error and the row index has to be rebuilt *inside* the loop: ~583 fat pointers, 161 times per
read group per sample, carrying nothing the width does not already carry. The replacement is
`GenotypeLikelihoodTable<'a>` over one flat row-major buffer with the width named once. It takes one
slice and one `usize`, so there is no pair of same-typed values to transpose; raggedness becomes
unrepresentable rather than checked on all 161 calls; and the second implementor makes the case
stronger, not weaker — an STR locus with 20 allele lengths has 210 diploid genotypes, and the row
index is pure overhead on a buffer that is already the dominant cost.

**⚠ A wrong number in the author's own prose again, and again about the author's own test's
reach.** Four rounds of eight on this plan have had one. This time it was the scratch-buffer
comment: *"the inner loop runs over every cell of the table, 161 times per fit"* swaps the fit/scan
nesting and understates the cell axis 3.6-fold. 161 is the ladder's **rung** count — one whole climb
per rung — and the table is **583** cells, walked once per pass. Three further wrong claims: the
start-independence test was described as three starts when it ran four and as covering each corner
when it missed the heterozygous one (the coordinate the estimator exists to measure); *"a cell at
depth 124 has a likelihood no `f64` can hold"* is true of the corners and not of a cell (an ordinary
124-read cell with no alternative read is 0.883); and *"3,052 passed"* is the library binary's count
attributed to a command that runs ten binaries.

**⚠ Four statements about convergence, in two files, that do not agree** — and one of them,
`FitTermination`'s *"the inner climb over the genotype frequencies is provably concave and needs no
cap"*, is now false of the code beneath it. The distinction the author meant is real and was stated
nowhere: **concavity rules out a false summit, not a slow one.** `spec/parameter_prepass.md` §3.1
proves there is no local maximum that is not global, and in the same paragraph says the rate of
approach is only linear and slowest where the components overlap — this estimator's own regime. It
is stated once now, on `MAX_CLIMB_PASSES`, and the others defer to it. **The same phrasing is
upstream in `arch/parameter_prepass_generic.md` §5.2** ("Why the cap, given the inner climb needs
none") and is left for the owner rather than edited here.

**⚠ The one thing an unverifiable anecdote cost.** The report's evidence for the pass cap was a
retired tetraploid fixture — "0.6977 against a truth of 0.7000" — that was not in the commit, not in
the file's history, and not described numerically, so the numbers agent could not reproduce it
(their reconstruction gave 0.666456). It has been replaced by the inbred fixture's 1,234 passes,
which is in the tree and re-measurable.

## Findings and what was done

### Blockers (3, all missing tests)

1. **Neither claim about the returned score is tested.** Two mutations survived: returning the
   climb's last-pass running score instead of the final `weighted_log_likelihood`, and
   `converged: converged || true`. Both survive because every fixture converged, and on a settled
   climb the two scores agree to 2.3 × 10⁻¹³ against the test's 10⁻⁶ tolerance. **Fixed** by
   `a_climb_stopped_early_says_so_and_scores_the_weights_it_returns`, which caps the climb at five
   passes — the two states are then far apart and the longhand check has teeth.
2. **A live `0 · −∞ = NaN`.** Deleting the zero-weight guard in `weighted_log_likelihood` left the
   suite green and `log_likelihood` came back `NaN`. The module already built that exact input —
   `a_cell_carrying_no_weight_is_skipped_whatever_it_says` has a zero-weight all-`−∞` row — but it
   went through `fit_mixture_weights` and never looked at the score. A `NaN` does not lose loudly in
   a profile scan; it simply never wins. **Fixed** by routing that test through the climb and
   asserting the score is finite; the guard now carries a comment saying it is not an optimisation.
3. **The stillness measure's maximum is untested** — above.

### Majors (5)

4. **The `&[&[f64]]` seam** — above. Applied, with the deviation recorded.
5. **Four validator refusals had no test** and could each be deleted or weakened silently: a `+∞`
   cell weight (the `-1.0` test is caught by `>= 0.0`, leaving `is_finite()` unpinned), a
   zero-width row, a wrong-length `start`, and a `+∞` start weight. **Fixed**: four tests, each
   verified against its mutation.
6. **The scratch-buffer comment's numbers** — above.
7. **The four disagreeing convergence statements** — above.
8. **Two behavioural defaults invisible from the public surface.** `MAX_CLIMB_PASSES` and
   `CLIMB_STILLNESS` were private, unnamed in the public function's docs, and cited no source.
   **Fixed**: both are `pub`, both are named from `fit_mixture_weights`'s docs, both cite
   `examples/ng_multilib_key_harness.rs`'s `climb_frequencies` (which caps at 400 and breaks out
   silently), and the cap's doc carries the three measured pass counts.

### Minors and nits applied

- `weights` → `genotype_weights` throughout, including the `MixtureWeightsFit` field.
  `weighted_log_likelihood` took `cell_weights` and `weights` adjacent, both `&[f64]` and one
  qualifier apart; transposed, both `zip`s truncate and it returns a score over three of the cells
  with the cell counts used as mixing weights, silently.
- `check_table` did three jobs; the table's own invariants moved into
  `GenotypeLikelihoodTable::from_natural_logs`, where they run once per table instead of 161 times
  per scan, and what is left is `check_cell_weights`.
- `climb_mixture_weights` is `pub(super)`, which is the boundary the design describes.
  `fit_mixture_weights` stays `pub` only because it is what keeps the module from being dead code
  until D3 lands; **it should become `pub(crate)` in D3**, when every consumer is in this crate.
- Bare adjectives named: `next` → `next_genotype_weights`, `moved` → `largest_move`, `largest` →
  `largest_term`.
- The monotone-ascent `debug_assert!` now skips the first pass explicitly rather than relying on
  `−∞ − ∞` to make it vacuous.
- `converged`'s doc says which direction the implication runs.
- A `# Examples` doc test on the public function, which is the only executable statement of the
  cell/genotype convention on the public surface. It runs under `cargo test --doc`, which the
  feature's normal validation set does not cover — run once here and green.
- `no_move_away_from_the_fitted_weights_scores_higher` says what it does **not** check: both points
  go through the same scorer, so a uniformly wrong `weighted_log_likelihood` passes it untouched.
  That half is `the_reported_score_belongs_to_the_weights_returned_with_it`.

### Checked and correct

The fixtures' construction holds: every column of `CELL_GIVEN_GENOTYPE` and of the tetraploid
`cell_given_dosage` sums to 1.0 with zero deviation, which is what makes the chosen truth the exact
maximiser. The cell weights span 11.5 to 1, so a climb ignoring them fails five tests. The eighteen
simplex steps all execute — none is skipped by the guard — and the smallest scores 9.13 × 10⁻³ below
the summit against a floating-point noise floor of 1.16 × 10⁻¹⁰. `cargo doc` is at its
12-unresolved-link baseline with none of the twelve in `parameter_estimation`. The Redner & Walker
and Dempster–Laird–Rubin attributions match `spec/parameter_prepass.md` §3.1 verbatim.

The `debug_assert!` question raised in the review prompt is answered: it is correct as written, and
a weight reaching exactly zero does **not** produce a `NaN`, because `−∞ + finite` is `−∞`. The one
reachable `0 · −∞` was Blocker 2.

## Verification

Container throughout. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → **3,060** in the library binary
(from 3,052) and 69 across the nine other binaries, unchanged; `cargo test --doc
ng::parameter_estimation` → 1 passed; `cargo doc --no-deps --lib` at 12 unresolved links, the
pre-existing baseline. `ng::parameter_estimation` 150 → **158** tests, `mixture_weights` 20 → **28**.

Seven mutations re-run against the fix, each failing exactly one test: the narrowed stillness
measure, the deleted zero-weight guard, the last-pass score, `converged || true`, and the three
weakened validators.

## Left for the owner

- **`arch/parameter_prepass_generic.md` §5.2 says the inner climb "needs none"** — no cap. It does
  have one, it is measured, and one of this file's own fixtures needs 1,234 passes. The code and
  `fitting/mod.rs` now say so; the architecture does not.
- **Whether 10,000 passes is enough for a real table** is still open and cannot be settled before
  F2. The generic path's table is 583 cells at three reads a site, which is exactly the regime
  §3.1 names as slowest. `climb_with_cap` is what makes that a measurement when the time comes.
