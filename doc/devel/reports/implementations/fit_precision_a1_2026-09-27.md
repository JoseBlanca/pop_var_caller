# Fit precision, step A1 — each parameter's slope at one position

**Date:** 2026-09-27. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A1.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §3.2. **Branch:** `fit-precision`.

## 1. What was built

The SNP/indel fit (`src/parameter_estimation/joint/fit.rs`) can now say, at any position, how
steeply that position's log-likelihood rises or falls as each parameter is nudged — the parameter's
**score** there. Summed over positions, the pairwise products of these scores are the information
the standard errors of steps A2 and A3 will be computed from.

The scorer is a new child module,
[fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs), so that `fit.rs`
grows only by what the scorer must read from it: four fields of the pass's per-position working
memory and five small functions beside the model they differentiate (4,880 lines before, 4,944
after). The scorer reads what `one_position` leaves in that working memory and writes nothing back,
so it cannot move a fitted number; nothing outside its own tests calls it yet — the expectation pass
will, at step A2.

A score row holds **eleven numbers for a sample-and-cohort block**: the cohort's eight (the mismapped
share, the density's two point masses and two Beta shapes, the duplicated class's share and its
carrier Beta's two shapes) and three per sample (its first read group's clean and mismapped error
rates, and its homozygote excess).

## 2. Assumptions and deviations

1. **A sample's second and later read groups get no score, because the likelihood never reads their
   rates.** `one_position` scores a sample's reads — pooled over all its read groups — under the rate
   of its first read group alone, while the maximisation fits every read group's rate from the same
   pooled tallies, so every later group's fitted rate is a copy of the first's. The spec (§3.2) sizes
   a sample's block as `1 + 2 × read groups`; under the likelihood as written the extra rows are
   identically zero, so the block here is three wide whatever a sample's read groups. The test
   `a_second_read_groups_rates_have_no_slope` pins the premise: moving a second group's rate leaves
   the log-likelihood bit-for-bit unchanged.

   **Consequence for the owner at checkpoint A.** Those read groups' errors will be *absent*. Under
   the owner's ruling that a parameter with no information is written as `defaulted` with the
   default value (spec question 1, built at step E2), each such group's fitted rate — a real number
   the fit produced, equal to its sample's first group's — would be replaced in the file by the
   default. On kimura (2,169 samples, 2,651 read groups) that is up to 482 read groups: the difference
   of the two counts, since the run log does not count later read groups directly. **Recommended:**
   write a later read group's rate with its first group's value and standard error, marked as shared
   with that group, since it is the same fitted number; the alternative, `defaulted`, discards a value
   the data did determine. This is a pre-existing property of the model, not something this step
   introduces, and the spec text about `k` is left for the owner to amend with the ruling.
2. **The two Beta shapes' slopes are the derivative of the quadrature the likelihood is computed on,
   not the digamma form spec §3.2's table names.** The pass integrates over a position's allele
   frequency on a Gauss–Jacobi rule — a fixed set of nodes and weights placed by the Beta itself.
   There are two ways to take the slope in a shape:
   - **the digamma form**: integrate `ln f − ψ(a) + ψ(a + b)` on the rule, where `ψ`, the digamma
     function, is the slope of `ln Γ`. This is what the spec's table names;
   - **the rule's own derivative**: as the shape moves, each node moves and each weight changes, and
     the likelihood is read at the moved nodes. This is exactly the slope of the log-likelihood the
     fit computes.

   Both are slopes of the same log-likelihood; they differ in how the integral over the frequency is
   taken, and only the second passes the spec's own oracle (§3.6 item 1). Measured on this step's
   drawn cohort at `a = 0.9`, against a finite difference of the log-likelihood that reads −20.9407282
   at every node count: the digamma form is **2.9% off at the shipped 16 nodes**, 4.7% at 12, 1.4% at
   24, 0.42% at 48, 0.12% at 96 and 0.049% at 160. `ln f`, unbounded at zero where the Beta puts its
   mass, is what the rule integrates badly. The rule's derivative is 9.6 × 10⁻⁹ off at 16 nodes and
   9.0 × 10⁻⁹ at 12. How the nodes and weights move depends on the shapes alone, so it is computed
   once a pass (`RuleSlopes`, by central differences of the rule); how the likelihood changes along a
   node's move is the genotype prior's derivative in the frequency, one more sum a sample. The test
   `the_digamma_form_integrated_on_the_rule_misses_the_slope_in_a_shape` pins both at 12 and 16 nodes.
3. **A finding this leaves for the owner: the maximisation's own update for the Beta shapes uses the
   digamma form** (`fit_beta_shapes`, over the pass's summed `ln f`), so the fit settles where that
   form, not the log-likelihood's slope, is zero. On this fixture the two slopes in `a` differ by 2.9%
   at the shipped 16 nodes. Fixing it would move fitted numbers and is outside this plan. **The
   measurement to bring to checkpoint A:** how far the fitted shapes sit from where the
   log-likelihood's slope in them is zero, in units of their standard errors — which needs step A3's
   errors.
4. **The error-rate slope under a depth *range* is exercised under a cap of 140**, because under the
   census's shipped cap of 124 every stored depth is one depth, so the range path never runs in a
   shipped configuration. The census ladder's first rung above 124 runs to 159 — 35 depths, more than
   the 32 the fit reserves room for — so a census whose cap leaves that rung wider than 32 depths
   panics in `next_position`; at a cap of 140 the rung is cut to a range of 16. Pre-existing, for the
   owner.
5. **No posterior is dropped for being small.** When the pass credits each branch and each node its
   share of a position's counts, it skips any whose posterior is below 10⁻¹², which costs a count
   nothing; a share's slope divides by the share, which can itself be 10⁻¹², so the scorer keeps every
   posterior and forms each as one exponential of a difference of logarithms. **No test can see this
   rule on the drawn cohort**: every share there is at least 0.006, and adopting the pass's skip was
   measured (in review) to move the summed slopes by about 10⁻¹¹ relative or not at all. The shares
   test's documentation records the limitation.

## 3. Changes

- [fit.rs](../../../../src/parameter_estimation/joint/fit.rs): `Scratch` keeps, per noise class, the
  three branch log-likelihoods before their shares multiply them (`fixed_alt_ln`, `segregating_ln`,
  `duplicated_ln`) and the position's log-likelihood (`position_ln`); `one_position` fills them from
  values it already computed. **The duplicated branch's sum keeps its original left-to-right order**
  (`ln d + summed − ln 3`), so the stored value is computed beside it rather than reused in it — a
  reassociated sum would move the branch's last bits. Beside the model they differentiate:
  `candidate_read_probability` (now also what `ReadLogs::of` calls — the same expression, so the same
  bits), its slope `candidate_read_probability_slope`, `reference_read_probability_slope`,
  `genotype_frequencies_slope_in_frequency` and `genotype_frequencies_slope_in_excess`. `mod
  information;` declares the child module.
- [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs): `PositionScores`
  (the reused row), the `cohort` and `sample` index modules, `RuleSlopes` and `SlopedRule` (how a
  Gauss–Jacobi rule's nodes and weights move with its shapes), `score_position`. The module carries a
  `#![cfg_attr(not(test), expect(dead_code))]` until A2 calls the scorer.

## 4. Tests

Each slope test compares the slope summed over a drawn four-sample, 600-position cohort against a
central difference of the pass's own log-likelihood (step 10⁻⁵ of the parameter), at parameters
chosen away from the maximum so the slopes are far from zero. Every sample's rates and excess differ
from every other's. **The disagreement is relative to the slope, or to 10⁻³ where the slope is
smaller**, and every test asserts it below 10⁻⁶.

| test | coordinates | largest disagreement (measured) |
|---|---|---|
| `the_shares_slopes_are_the_derivative_of_the_log_likelihood` | mismapped, invariant, fixed and duplicated shares | 1.41 × 10⁻⁸ |
| `the_beta_shapes_slopes_are_the_derivative_of_the_log_likelihood` (16 nodes) | density `a`, `b`, carrier `a`, `b` | 1.70 × 10⁻⁷ (carrier `a`, slope −0.203) |
| `the_error_rates_slopes_are_the_derivative_of_the_log_likelihood` | 4 samples × clean, mismapped | 1.00 × 10⁻⁸ |
| `the_error_rates_slopes_are_right_where_a_depth_is_a_range` (cap 140, 132 reads; 1,729 of 2,400 sample-positions ranged) | 4 samples × clean, mismapped | 8.63 × 10⁻⁹ |
| `the_homozygote_excess_slopes_are_the_derivative_of_the_log_likelihood` | 4 samples | 1.54 × 10⁻⁷ |
| `the_slopes_hold_with_the_duplicated_class_off` | 5 cohort, 8 rates, 4 excesses; the class's three slots exactly zero | 1.27 × 10⁻⁷ |
| `the_slopes_hold_with_held_out_reads_and_coverage_odds` (583 sample-positions with a held-out read; odds 0.25, 3, 1, 8) | carrier `a`, `b`, duplicated share, 8 rates | 1.38 × 10⁻⁷ |
| `a_second_read_groups_rates_have_no_slope` | a second group's two rates: log-likelihood unchanged bit for bit; the first group's clean rate | within 10⁻⁶ |
| `the_digamma_form_integrated_on_the_rule_misses_the_slope_in_a_shape` (12 and 16 nodes) | density `a` | digamma form 4.71 × 10⁻² and 2.88 × 10⁻² off (asserted above 3 × 10⁻² and 2 × 10⁻²); this module's 9.00 × 10⁻⁹ and 9.56 × 10⁻⁹ |
| `the_models_slopes_are_the_derivatives_of_its_functions` | the five helper slopes in `fit.rs`, over a grid | within 10⁻⁷ absolute |

Numbers from `cargo test --release --lib parameter_estimation::joint::fit::information -- --nocapture`
in the container (`tmp/fit_precision/module_a1c.log`): 10 passed.

## 5. Validation

- `cargo fmt --check`: clean.
- `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` on the committed tree
  (`tmp/fit_precision/suite_a1c.log`): 4,890 passed, 3 failed, 4 ignored. The three failures are
  pre-existing and outside this work — `examples/ng_generic_loci_dump.rs`
  (`a_deletion_across_a_region_boundary_keeps_the_support_a_single_region_walk_gives`,
  `only_the_rows_that_departed_from_the_reference_carry_a_chain_id`) and
  `examples/ng_ssr_loci_dump.rs` (`a_cap_above_the_depth_is_invisible_and_a_cap_below_it_bites`).
  **Both cross-platform checksum tests (`cli::cross_platform_digests`) pass unchanged**: no fitted
  number moved.
- The identity oracle (`scripts/promote_ng_oracle.sh`) was run on this branch's starting point
  (§6), not on this step; the checksum tests above are this step's evidence.

## 6. The identity oracle's baseline is stale on main

`scripts/promote_ng_oracle.sh` over tomato2's first 20 regions and four accessions, run on this
branch's starting point, reproduces the five default-parameter checksums in
`scripts/promote_ng_oracle.baseline` and **differs on the two fit lines**: `fitted.parameters.toml`
`c590689d40eb66359b2973c2bd1cbfcd` and `from_fit.comparable` `5cdeafd0aaa805cfd57b79e3bbb4e3d2`
(`tmp/fit_precision/oracle_main/digests.txt`). The baseline was last changed in `8eeac56d`, before the
change of 2026-09-25 (`33527ebd`) that accelerated the SNP/indel fit by extrapolating along two
successive steps (SQUAREM, Varadhan & Roland 2008); that commit re-recorded the cross-platform
checksums in `src/cli/cross_platform_digests.rs` but not this file. That it is the commit which moved
the two lines is inferred from the history, not measured. Milestone A's "nothing moves" is judged
against these two values.

## 7. Follow-ups

- Owner, checkpoint A: the later read groups (§2 item 1, with a recommendation), the maximisation's
  digamma-form update for the Beta shapes (§2 item 3), and the depth-range panic above a cap of 124
  (§2 item 4).
- Step A2 wires the scorer into the expectation pass and removes the module's `expect(dead_code)`.
  Deferred to it from this step's review: bundling the inputs `score_position` shares with
  `one_position` into one value; computing each read's slope once per (class, candidate, sample,
  genotype) rather than inside the node loop; splitting `score_position` by branch; and names for the
  parameter slots, which A3, A5 and B3 need to print.
