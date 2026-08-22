# ng genotype prior — D1: what spectrum a candidate pair predicts, in closed form

*Implementation report, 2026-08-22. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step D1 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone D. Includes the review and the
fixes applied from it.*

## 1. What it is

Given a candidate concentration pair — `α_ref` and `α_alt`, chromosomes' worth of prior belief —
`fill_expected_spectrum` says what a panel of `N` diploid individuals would look like: over many
sites, what fraction carry the alternative allele on exactly `j` of the panel's `2N` chromosomes.
`2N + 1` classes, summing to one.

Step D2 fits the pair by searching for the one whose prediction matches the spectrum the pre-pass
measured. So this function is both **D2's objective** and the way **D2's test targets** are built.

Design authority: [`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §3.2 (the two-branch
sampling), §4.1 (the projection), §12 tests 5–7; [`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §4.

## 2. The shape, and why it is summed this way

Two draws. The locus's alternative-allele frequency comes from a Beta with the candidate pair; then
each individual is drawn at that frequency under the **same two-branch model the genotype prior
itself uses** — with probability `F` its two copies are one ancestral copy counted twice, otherwise
two independent draws. Averaging over the frequency is what this computes.

**The decomposition.** Split the panel by how many individuals took the identical-by-descent
branch. With `M` of them inbred, the panel holds `2N − M` *distinct* chromosomes — the other `M`
are copies of one of those — and each distinct chromosome carries the alternative allele
independently. The class is then the number of alternative distinct chromosomes plus however many
of the duplicated ones are among them.

**Every term is a product of non-negative factors, and that is the point.** Written instead as a
polynomial in the frequency, the coefficients alternate in sign and grow about 8.3-fold per
individual — 9.5e9 by twelve individuals — and would have to cancel against a total that must come
out at one. These do not cancel at all.

## 3. The mathematics survived three independent oracles

The review's mathematics agent built three checks sharing no code with the function:

| oracle | worst relative gap |
|---|---|
| exhaustive enumeration of all `5^N` (branch, genotype) panels with exact Beta moments, `N ≤ 5` | **9.4e-15** |
| Simpson quadrature of an `N`-fold Wright convolution | **7.6e-14** |
| a Wright generating function expanded as a bivariate polynomial | 1.5e-7, and that residue is the *oracle's* own cancellation |

The loop bounds are exactly the legitimate range at all 35,301 `(N, M, a)` triples with `N ≤ 40`,
and the inner multiplicity sums to `C(D, a)` by Vandermonde.

## 4. What the review found

Three agents in isolated worktrees — mathematics and numerics, reliability, and naming with errors
and smells. The formula was right; everything below is about what the step left unguarded, what it
cost, and what it claimed.

### The cost, and why it reaches into D2

**This report's first draft said the sum is "paid once per run".** It is not. D1 is D2's objective,
so a multistart fit evaluates it on the order of a hundred times. That difference is the finding
with the longest reach, and it was worth two fixes:

- **Tabulate the log-factorials once per call** rather than recomputing `lgamma` at every term.
- **Skip branch splits too rare to matter.** For a fixed `M` the classes are themselves a
  distribution summing to one, so such a split can move no class by more than its own probability;
  the skip holds that below `1e-300`. Without it, an inbreeding coefficient of one in a million was
  five orders of magnitude slower than one of exactly zero, for the same answer — the `None` path
  fired only at exactly 0 and exactly 1.

Measured in release at `F = 0.8`, one prediction, before and after:

| individuals | before | after |
|---|---|---|
| 63 | 3.7 ms | **241 µs** |
| 200 | — | 5.8 ms |
| 400 | — | 43 ms |
| 800 | 6.55 s | **354 ms** |

About 15-fold, and still cubic — a factor of 8.2 across the last doubling, an exponent of 3.03.
**So a fit costs about a minute at 800 individuals and hours by several thousand.** Whether the
projection can run at the top of the committed cohort range is **step D2's question**, and it has
room to answer it without this function changing: bin the classes, or cap the panel it projects at.
Recorded here and raised at Checkpoint D.

### The Blocker: two release-held checks that nothing could catch

`α_ref > 0` and `α_alt ≥ 0` are the only guard in release — `lgamma`'s own check is a
`debug_assert!`. Deleting both left all eight tests green in **both** profiles. With them gone,
`α_ref = −0.5` returns nine finite, non-negative class weights totalling 1.0097: a spectrum that
looks like a spectrum. Two `#[should_panic]` tests now hold them.

### A concentration ceiling, because D2 searches this axis

The precondition admitted any finite positive concentration, and past about `1e9` the sum stops
agreeing with itself. Measured at 63 individuals, `F = 0.9`: the classes total 1.0000016 at `1e9`,
1.0022 at `1e12`, and **1,107 at `1e15`** — every entry finite and non-negative, so nothing
downstream could tell. `MAX_PROJECTION_CONCENTRATION` is now `1e6` chromosomes, which is half a
million diploid individuals and two orders past the committed range, with the total within 2e-11 of
one there.

### No oracle strictly inside the inbreeding range — proved, not argued

Both exact tests sat at `F = 0` and `F = 1`, where the triple sum collapses to one term per class.
The mathematics agent demonstrated the gap: replacing `class = draws + doubled` with `class =
draws` — deleting the doubling, which is a **different model** — still passed the sum-to-one test,
the beta-binomial test *and* the neutral-limit test, because the inner sum is a Vandermonde
identity either way.

Two oracles now close it, and they are complementary:

- **The first two moments**, which the two-branch model fixes exactly:
  `E[j] = 2N·E[p]`, free of `F`, and `E[j²] = 4N²·E[p²] + 2N(1 + F)·(E[p] − E[p²])`, which carries
  `F` linearly. Cheap at any panel size, and it works at **one individual**, where the
  doubleton test has no two classes to compare. Derived independently before adoption.
- **Exhaustive enumeration** of every `(branch, genotype)` assignment at `N ≤ 5` — `5^N` panels,
  each a monomial in the frequency averaged exactly. It shares the Beta-moment formula and nothing
  else: no split by inbred count, no hypergeometric counting, no binomial coefficients.

### One of my own new tests was asserting something untrue

The test for the branch skip first asserted that `F = 1e-6` agrees with `F = 0` to floating point.
It does not, and should not — the two genuinely differ by order `F`, measured 3.0e-10 at three
individuals. It now asserts what actually holds: that no mass is dropped (the classes still sum to
one) and that the two are within order `F`, with the real oracle for the skip being the moment
identity, which now runs at `F = 1e-6` too.

### A gate I had not been running

`cargo clippy --lib` does not lint test code, but `scripts/precommit-check.sh` and CI both run
`--all-targets`. My tests carried a `needless_range_loop`. Fixed, and **`--lib --tests` is now part
of this step's gates** rather than `--lib` alone.

### Doc claims corrected, each because it was measured false

| claim | measured |
|---|---|
| the polynomial's coefficients "reach `4^N`" | they overshoot it 564-fold by twelve individuals; growth is about 8.3× per individual |
| "Measured in this module's own tests" | no timing test existed; there is one now, `#[ignore]`d |
| "9 to 14% at tomato's fitted `F`" | spec §4.1 puts it at 12–14% at `F` 0.8–0.9; 8.6% is the `F = 0.6` figure |
| "0.28%" neutral departure | 0.272% |
| binomial coefficients "reach `10^600` at the largest panel here" | `10^36.8` at 63 individuals; `10^600` is the thousand-individual figure |
| "rather than measuring the fit against itself" | backwards — D1 supplies both D2's targets and D2's objective, so D2's tests check its *search*, not its mathematics |

**The 12–14% error is also in `arch/calling_priors.md` §4 and in this plan's step D2 line**, both of
which say 9–14%. Those are design documents and are the owner's to correct; raised at Checkpoint D.

### Deleted rather than fixed

The `individuals == 0` early return was redundant — the general sum already returns exactly 1.0 in
class 0, bit pattern and all, in every combination probed. Both agents said so; it is gone.

### Not taken

**A checked type for the concentration pair**, as step C1 got for its two copy-count arrays. The
naming agent recommended against and I agree: a swap is not silent here — `α_ref ≈ 1` against
`α_alt ≈ 6e-4` puts all the mass in the top class and every oracle fails — there is one intended
caller, and `SpectrumSeed` already holds this pair, so D2 can pass that if it wants the axis closed
for free.

## 5. Tests

Fifteen, plus one `#[ignore]`d timing test. The module goes from 62 passed debug / 54 release to
**69 / 61**, with 1 ignored.

| test | what it pins |
|---|---|
| `the_classes_of_a_panel_sum_to_one` | normalisation over six panel sizes, five diversities, three reference concentrations and five inbreeding coefficients. Worst departure 1.7e-13 against a `1e-9` budget, asserted |
| `at_no_inbreeding_the_spectrum_is_the_beta_binomial` | the exact oracle at one end of the `F` range |
| `at_full_inbreeding_only_the_even_classes_can_hold_anything` | the exact oracle at the other, and that the doubling is a doubling — 26 selfers carry the information of 26 chromosomes, not 52 |
| `the_first_two_moments_match_the_two_branch_model_at_every_inbreeding_coefficient` | **the oracle inside the range**, from `F = 1e-6` to 0.999, down to one individual |
| `every_panel_enumerated_one_by_one_gives_the_same_spectrum` | the whole model written a second way, at `N ≤ 5` |
| `doubletons_beat_singletons_only_once_the_panel_is_inbred` | the signature no independent-chromosome spectrum can produce — the reason the two-branch model is here at all |
| `the_neutral_shape_appears_in_the_small_diversity_limit` | that `θ/k` is reached in the limit, and by how much it is missed outside it |
| `a_vanishing_inbreeding_coefficient_drops_no_mass` | that the branch skip is a saving, not an approximation |
| `a_reference_concentration_at_zero_is_refused` / `a_negative_alternative_concentration_is_refused` | the two release-held value checks |
| `a_concentration_past_the_computable_range_is_refused` / `the_spectrum_still_sums_to_one_at_the_largest_concentration_accepted` | the ceiling, from both sides |
| `a_cohort_with_no_alternative_concentration_puts_every_site_in_the_monomorphic_class` | a fully invariant cohort is an answer, not a division by zero |
| `a_single_individual_still_has_a_spectrum` | the low end of the committed range |
| `a_mis_sized_class_buffer_is_refused` | the class count |
| `cost_of_one_prediction_by_panel_size` | the timings this report quotes, `#[ignore]`d |

## 6. Mutations re-run after the fixes

Every one applied against the fixed tree, and the file verified **byte-identical** to a pristine
copy afterwards — by content, not by a summary count:

| mutation | before | after |
|---|---|---|
| the doubling deleted (`class = draws`) | passed 3 of the exact tests | **killed** — 5 tests |
| `F` swapped with `1 − F` | passed sum-to-one | **killed** — 4 tests |
| the inner loop's lower bound forced to zero | survived in release | **killed** — 8 tests |
| the inner loop's upper bound lowered by one | killed | **killed** — 8 tests |
| the reference-concentration check removed | survived both profiles | **killed** — 2 tests |
| the negligible-branch skip widened to `1e-3` | — (added by the fix) | **killed** — 5 tests |

## 7. The projection was then made fast, at no measurable cost

*Added 2026-08-22, at the owner's direction to solve the scaling before building D2.*

Three ways of computing the same sum were built and compared —
[the report](../../ng/reports/spectrum_projection_cost_2026-08-22.md),
[the harness](../../../../examples/ng_spectrum_projection_cost.rs). The shipped function now writes
each term as a beta-binomial weight times a hypergeometric one, steps the second by an exact ratio
rather than exponentiating it, and drops branch splits below `1e-18` of the likeliest one.

| samples | before | after | worst class error |
|---|---|---|---|
| 400 | 43.8 ms | **5.8 ms** | 2e-13 |
| 800 | 339.6 ms | **29.9 ms** | 4e-13 |
| 1,600 | 2.1 s | **179.4 ms** | 4e-13 |
| 3,200 | 12.1 s | **960.3 ms** | 1e-12 |

`N^2.95` becomes `N^2.45`, so a fit at 3,200 samples is 2.6 minutes rather than 32 — **the whole
committed cohort range is now a once-per-run cost, and D2 can be built as spec §4.1 writes it**,
over every class including the monomorphic one. The accuracy figures are the same disagreement the
*untrimmed* version has with the term-by-term sum: floating-point accumulation, not the trim.

The term-by-term sum is kept in the test module as the oracle the fast one is checked against, and
two defects the change introduced are recorded in the report — a multiplicative walk that must
start at its mode or lose whole rows above about a thousand samples, and a `usize` subtraction that
release was wrapping while debug refused it.

## 8. Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `69 passed; 0 failed; 1 ignored` |
| `cargo test --release --lib ng::calling::genotype_prior` | 0 | `61 passed; 0 failed; 1 ignored` |
| `cargo test --lib` | 0 | see the commit message |
