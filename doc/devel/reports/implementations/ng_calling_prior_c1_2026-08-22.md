# ng genotype prior — C1: the concentration one sample gets, with its own reads taken back out

*Implementation report, 2026-08-22. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step C1 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone C. Includes the review and the
fixes applied from it.*

## 1. Plan

`fill_sample_concentration`: the run's starting concentration plus what the **other** samples showed
at this locus.

```text
α'_s(a) = seed(a) + max(0, cohort expected copies of a − this sample's own)
```

A port of `leave_one_out_alpha` ([`ssr/cohort/em.rs`](../../../../src/ssr/cohort/em.rs)) and its SNP
twin inside [`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs). Design
authority: [`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §6 and §12 tests 8–9;
[`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §3.1.

## 2. What the step is for, in one paragraph

The cohort's expected allele copies are estimated from every sample **including this one**. Leave
this sample's own contribution in, and its reads arrive twice — once through the read likelihood,
and once through the allele frequency they helped estimate. Taking them back out is what makes the
prior a prior. The spec is explicit that this is *not* the mechanism behind the GIAB 214-site
failure — the starting concentration is — and that it belongs here anyway, because using a sample's
reads twice is wrong and needs no measurement to justify (§6).

**Both ends of the committed cohort range are one formula with no branch.** At one sample the cohort
total and the sample's own copies are the same number, so the term is exactly zero and the output is
the seed bit for bit. At several thousand samples the term swamps the seed and the prior converges
on the panel's own frequencies.

## 3. The port is exact, and that was checked rather than assumed

The review re-spelled both frozen production expressions verbatim and compared bit patterns over
**20,000 random cases at 1 to 8 alleles**: ng, the STR spelling and the SNP spelling agree **bit for
bit**. Unlike step B1, there is no ulp gap between production's two copies here.

Where the two production spellings *do* differ is in what they check, and ng carries the union:

| | length checks | negative-difference check |
|---|---|---|
| `leave_one_out_alpha` (STR) | debug | none |
| the SNP twin | none | debug |
| **ng** | **release** | debug |

The length checks are promoted because `out` is the caller's buffer and is reused across loci: a
short one leaves the previous locus's entries standing in this locus's prior. Neither production
spelling can suffer that — one allocates its output, the other sizes its scratch once per locus.

**One-sample bit-exactness is real across the range, not a fixture artefact.** Checked over 20,014
magnitudes from the smallest subnormal to `f64::MAX`: `x − x` gave `+0.0` every time and `−0.0`
never, and `seed + 0.0` was bit-identical to `seed` for all 10,375 seeds at or above the floor.

## 4. What the review found, and what changed

Three agents in isolated worktrees — reliability, port fidelity and numerics, naming with errors and
defaults. **Mutation totals: 34 run, 5 survived, 1 changed no behaviour.**

### The two that could produce a silently wrong prior

**Swapping the cohort's copies and the sample's own was silent in release.** The flat-slice
signature had four `&[f64]` parameters, so the compiler could not tell two of them apart. Swapped,
the function returns the bare seed at every allele — the cohort's evidence gone, nothing raised.
The debug assertion caught it only where the gap exceeded `1e-6`, which misses one sample entirely
(the two arrays are equal by construction there) and misses small gaps everywhere.

**Fixed with two borrowing newtypes**, `CohortAlleleCopies` and `SampleAlleleCopies`, in the same
`checked` module as `Concentration` so nothing in this folder can build one without the check.
Verified: passing them in the wrong order is now `error[E0308]`, with rustc naming both types.
Neither owns anything, so nothing allocates.

**A `NaN` copy count was swallowed.** `f64::max` returns the other operand on a `NaN`, so a `NaN`
difference became `0.0` and the allele came back carrying nothing but its seed — a plausible-looking
number with the evidence quietly gone. An infinity passed in **both** profiles and left an infinite
concentration, contradicting this function's own claim to return a valid `Concentration`. The two
new types check finite-and-non-negative when they are built, in debug, which is where this module
puts every check on a value; the loop's own `ExpectedAlleleCopies` makes the same check in release,
once per locus rather than once per sample per pass.

### The Blocker: a release check with no test

Of the three release-held length checks, `own_expected_copies` was the only one untested.
Downgrading it to `debug_assert_eq!` left the module green in both profiles while the trailing
alleles silently carried the previous locus's entries. Now `a_short_own_count_array_is_refused`
kills it, and only in the release run — which is what the plan's two-profile verification earns.

### Two tests that could not fail

**`the_result_is_accepted_as_a_concentration` ran where its own mechanism was inert.** Every
difference in its fixture was non-negative, so deleting the `max(0, ·)` left it green. Its fixture
now gives one allele a noise-negative difference, which is exactly where a missing clamp pushes an
entry under the floor.

**The desync threshold could only drift in the direction that hides defects.** Tightening `-1e-6`
to `>= 0.0` was caught; widening it a thousandfold to `-1e-3` was invisible to every test in the
module. Two bracketing tests now hold both sides — **written as literals, not as fractions of the
constant**, because a test that derives its input from the constant it pins moves with it and can
never fail. (That defect was in the first draft of these tests and was caught here, not by the
review.)

### Three claims I wrote about production that were wrong

Two agents caught the same one independently.

- *"Both production spellings return a `Vec` per sample per pass."* **False for the SNP spelling**,
  which writes a reused scratch buffer sized once per locus — already the shape this function
  claims as an improvement. It is the STR spelling that collects a fresh `Vec`.
- *"Production checks the same thing in debug only, because it allocates its output."* **False**:
  the SNP spelling holds no length check at any level.
- *"Caught in debug exactly as production catches it (`em.rs`, `posterior_engine.rs`)."* **Only the
  second.** `leave_one_out_alpha` has no negative-difference check at all.

### One number the review confirmed rather than corrected

**`-1e-6` is well calibrated for this caller's range.** Modelling the exact pair production names —
a biallelic fast path forming a total by complement against a per-row accumulator — the two-path gap
reaches 3.3e-10 at 5,000 samples, about 3,000 times inside the threshold, and would not reach the
threshold until roughly two million samples. It is now a named constant,
`COUNT_PATH_DESYNC_THRESHOLD`, carrying that measurement.

### Deliberate departures, both recorded for ratification

1. **`sample_concentration` → `fill_sample_concentration`.** [`arch §3.1`](../../ng/arch/calling_priors.md)
   pins the first name. It reads as a noun for a function that returns nothing and fills a buffer,
   and *sample* reads as a verb in a statistics context — "sample the concentration" is the wrong
   sentence entirely. The folder's other three buffer-fillers are all `fill_*`. **Owed: a one-word
   edit to arch §3.1.**
2. **The two copy-count arrays are newtypes, not slices.** Arch §3.1 sketches four bare `&[f64]`.
   This is the same move review forced at step A2, where eight flat parameters became one checked
   bundle — and for the same reason, that the compiler cannot otherwise refuse a caller bug.
   **Owed: the signature in arch §3.1.**

## 5. Tests

Fourteen, of which four are debug-only. The module goes from 40 passed debug / 36 release to
**54 / 46**.

| test | what it pins |
|---|---|
| `at_one_sample_the_concentration_is_the_seed_bit_for_bit` | spec §12 test 8. Its doc now says plainly **what it does not pin** — with the two arrays equal, an implementation ignoring, halving or swapping them passes here too (measured, all three do) |
| `the_leave_one_out_term_is_the_cohorts_evidence_less_the_samples_own` | the identity itself, entry by entry — added because the monotonicity below is satisfied by "constant", "halved" and "swapped" alike |
| `raising_the_cohorts_evidence_never_lowers_an_alleles_weight` | spec §12 test 9; the accumulator is seeded below every reachable value rather than at `−∞`, so the first of five rises carries an assertion |
| `float_noise_below_zero_leaves_the_seed_untouched` | the `max(0, ·)` on noise |
| `a_difference_just_inside_the_desync_threshold_is_absorbed` | the noise side of `-1e-6` |
| `a_difference_just_outside_the_desync_threshold_is_refused_in_debug` | the defect side |
| `a_cohort_total_below_the_samples_own_is_refused_in_debug` | a real desync |
| `a_short_output_buffer_is_refused` / `a_short_cohort_count_array_is_refused` / `a_short_own_count_array_is_refused` | the three release-held length checks, one test each |
| `an_empty_seed_is_refused` | and that the message names the *seed*, not just "a concentration" |
| `a_nan_cohort_copy_count_is_refused_in_debug` / `an_infinite_own_copy_count_is_refused_in_debug` | the value checks on the two new types |
| `the_result_is_accepted_as_a_concentration` | that the output is always wrappable, with a fixture where the clamp is live |

## 6. Mutations re-run after the fixes

Every one applied against the fixed tree, in the container:

| mutation | before | after |
|---|---|---|
| the two arrays swapped at the call site | survived in release, silently returning the seed | **compile error** (`E0308`) |
| `own` length check downgraded to debug | survived both profiles | **killed** — release only |
| the desync threshold widened 1,000× | survived | **killed** — debug |
| `max(0, ·)` removed | survived | **killed** — 3 tests |
| the cohort's evidence used without subtracting own | killed | **killed** — 8 tests |
| the leave-one-out term halved | survived | **killed** — the identity test |
| the copy-count value check neutered | — (added by the fix) | **killed** — 2 tests |

## 7. Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `54 passed; 0 failed` |
| `cargo test --release --lib ng::calling::genotype_prior` | 0 | `46 passed; 0 failed` |
| `cargo test --lib` | 0 | see the commit message |

## 8. One thing no test can hold

**Spec §6's "no branch on the cohort size" is not test-enforceable.** An implementation with an
explicit `n == 1` branch returns bit-identical values, so no fixture can separate it from this one —
measured. It is a shape requirement held by review, and the test's own doc says so rather than
implying coverage it does not have.
