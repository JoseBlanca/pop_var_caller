# ng calling loop — C2: convergence, the cap, and the emitted flag

**Date:** 2026-08-25
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step C2
**Design authority:** [spec/calling_em_loop.md](../../ng/spec/calling_em_loop.md) §6, §7, §13
tests 1 and 4; [arch/calling_em_loop.md](../../ng/arch/calling_em_loop.md) §4
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`

> **Read this against the review that followed it.** Five agents returned **one Blocker and four
> Majors** on code whose tests were all green. Two changed the shipped signature: `ploidy` was a
> second source of truth beside the genotype table's own, and a mismatch was silent in both
> profiles (a `ploidy` of 64 against a diploid table reported `passes: 2, converged: true`
> where the truth is `passes: 4`); and the `.abs()` in the stopping rule **could not fail at any
> fixture in the file**, because at two alleles the two movements are equal and opposite by
> construction. The other three were untested properties, not wrong code. **Two sentences of
> this report's own were wrong**, both explaining a mechanism rather than stating a number.
> [The review](../reviews/ng_calling_loop_c2_2026-08-25.md) carries all of it; this report has
> been corrected to describe what shipped.

---

## 1. Plan

C1 built the first pass and B1/B2 built the two halves of a pass. **C2 is where the passes
become a loop**: it decides when to stop, what to report when it does not, and how many passes
it took.

Three things land, all in
[`summarise_condition.rs`](../../../../src/ng/calling/inference/summarise_condition.rs):

- **`cohort_expected_copies_have_settled`** — the stopping rule. The largest change in the
  cohort's expected allele copies between two passes, **divided by the cohort's chromosome
  count**, against the configured threshold. It takes a `Ploidy` and a sample count rather than
  the product, so the two cannot be transposed and an infinite chromosome count is not
  expressible.
- **`FrequencyLoopOutcome`** — `passes` and `converged`, the two fields that travel into
  `LocusInference` unchanged. `#[must_use]`, because with no `Result` in this module it is the
  only carrier of *did not settle*.
- **`run_frequency_loop`** — the prior-free initialisation, then seeded passes until the copies
  stop moving or the cap is reached. It reads the ploidy off the genotype table it is already
  given rather than taking one.

## 2. Changes made

### 2.1 The stopping rule, and why it is spelled `all(… < threshold)`

**The division is the load-bearing half and it is the easy one to drop.** Expected copies are a
count; the threshold is a fraction. At one diploid sample the cohort carries 2 chromosomes and
at a thousand it carries 2,000, so a criterion written on raw counts tightens by the cohort size
across exactly the range this caller commits to. Production makes the same division and its
comment gives the reason ([`posterior_engine.rs:2702`](../../../../src/var_calling/posterior_engine.rs)).

**The spelling is not a style choice.** `prepare_for_locus` fills *both* cohort rows with the
`NaN` sentinel, and `advance_cohort_expected_copies` leaves the previous row holding it until a
pass has actually advanced. `f64::max` is documented to return the *other* argument when one
side is `NaN`, and a `>` comparison against one is false — so `fold(0.0, f64::max)` and a
hand-written maximum both hand back `0.0`, and `fold(-inf, f64::max)` hands back `−∞`. All three
are below any threshold. **Three of the four natural spellings therefore report a locus settled
after one pass, having compared it against nothing**, and the failure is a genotype flagged
`converged` that was never checked. `all(|d| d < threshold)` is the one that behaves as
`previous_cohort_expected_copies`' own doc comment promises. All four are computed in
`the_fold_spellings_of_the_delta_settle_where_this_one_does_not`, so the reason is a fact in the
suite rather than a remark in a doc comment.

**`run_frequency_loop` itself never hands over such a row**, and the doc comment now says so
*(corrected after review — it claimed the opposite)*. The prior-free initialisation's M-step
writes finite copies before the first swap, so by the time the rule is first called the previous
row is a real estimate. The guarantee is for C3's final pass and D1's outer rounds, which
restart the initialisation and could reorder the swap.

**And `.abs()` is load-bearing only from three alleles on**, which is why nothing in the first
draft could have caught its removal: expected copies sum to the cohort's chromosome total, so at
two alleles the two movements are exactly equal and opposite and signed and absolute comparisons
agree by construction. At three, the reference allele can fall by more than the threshold while
every alternative rises by less.

### 2.2 The loop, and where the swap sits

One pass is: E-step every sample against the cohort row as it stands → swap that row into
`previous` and take a `NaN`-filled buffer → M-step fills it → compare. **Swapping before the
E-step is the mutation this ordering exists to refuse**, and it is silent in release (§4).

**The prior-free pass is not counted as a pass.** It mints the estimate pass 1 is compared
against, which is what makes `LocusInference::passes` documented "at least one" true and a
one-pass locus a real outcome rather than a counter that never incremented.

**The settled test comes before the cap test**, and the order is a claim rather than a
formality: a locus whose last allowed pass is the one that settles is *converged*, not capped.
§6 makes the flag a statement about the locus, so reporting the cap there would understate
every genotype at the site.

**No line branches on cohort size**, per spec §7. At one sample the loop settles in two passes
by arithmetic, not by a test.

## 3. Deviations

**One, and it is a name.** The plan's step title says "the emitted flag", and the flag itself
already existed — `LocusInference::converged` and `::passes` were built in A1. What C2 adds is
the thing that *produces* them, so the two values arrive as a small named type
(`FrequencyLoopOutcome`) rather than as a bare `(u32, bool)`. C3 reads its two fields straight
into the `LocusInference` constructor.

**And one thing deliberately not done: the pass cap has no ceiling.** `max_passes` ships at 50
with a floor of one (`NonZeroU32`) where production caps its analogue at 500. Left as it was —
the consequence of a large cap is a slow locus rather than a wrong one, and spec §12's question
4 is what sets it from a real pass-count distribution.

## 4. Tests added

**Seventeen** — ten as written, seven the review added.

| test | what it pins |
|---|---|
| `at_one_sample_the_loop_settles_on_the_second_pass` | spec §13 test 1 — the initialisation and pass 1 differ; pass 2 equals pass 1 bit for bit; `passes = 2` |
| `a_locus_that_runs_out_of_passes_is_called_and_says_so` | spec §13 test 4 — capped at 2 the same locus reports `converged = false`, at 50 it settles in 4, and the flag survives into a `LocusInference` |
| `run_frequency_loop_reports_converged_when_the_last_allowed_pass_settles` | ⚑ the settled test wins over the cap test — a locus that settles on its last allowed pass is converged, not capped |
| `the_same_frequency_scale_movement_settles_at_one_sample_and_at_a_thousand` | the division — the same fraction of each cohort's chromosomes gives the same verdict at n = 1 and n = 1,000 |
| `a_fall_larger_than_every_rise_has_not_settled` | ⚑ the `.abs()`, which only three alleles can show |
| `a_movement_exactly_at_the_threshold_has_not_settled` | ⚑ the strict `<`; `2e-3 / 2.0` is exactly `1e-3` in `f64` |
| `the_fold_spellings_of_the_delta_settle_where_this_one_does_not` | the three rejected spellings settle against the sentinel row; the shipped one does not |
| `the_loop_reproduces_the_passes_driven_by_hand` | the loop is the hand-written sequence, both cohort rows compared bitwise after four passes |
| `each_sample_is_scored_against_its_own_inbreeding_coefficient` | ⚑ the pairing the length check protects — trade two samples' coefficients without moving the samples |
| `the_loop_settles_at_three_alleles_and_the_copies_still_sum_to_the_chromosomes` | ⚑ the loop past two alleles at all, and the sum invariant that must hold at any allele count |
| `a_convergence_test_over_no_alleles_is_refused` | `all` over nothing is `true` — an empty row would settle every locus |
| `cohort_rows_of_different_lengths_are_refused` | `zip` stops at the shorter row |
| `a_convergence_test_over_a_cohort_of_no_samples_is_refused` | no samples is no chromosomes, and nothing would ever settle |
| `a_threshold_that_is_not_a_fraction_is_refused` | the check at the point of use, for a threshold built by hand |
| `an_infinite_threshold_is_refused` | ⚑ the `is_finite()` half, which `> 0.0` alone leaves unreached |
| `fewer_inbreeding_coefficients_than_samples_are_refused` | a short slice scores fewer samples and the M-step sums the rest as current |
| `more_inbreeding_coefficients_than_samples_are_refused` | ⚑ the mirror; without it the panic is a scratch-indexing message two modules away |

**⚑ marks the seven the review added**, each written against a mutation that had survived.

### The numbers behind the fixtures, measured rather than assumed

- **One sample, the initialisation against pass 1:** `[1.7279, 0.2721]` against
  `[1.9998, 0.0002]`. That is 0.272 copies of the alternative allele, against the 0.002 raw
  copies the `1e-3` threshold allows over one diploid sample's two chromosomes.
- **Three samples pulling apart** (one sample's reads favouring each of the three diploid
  genotypes by 2 nats, about 8.7 Phred): the loop settles in **4** passes at the shipped cap.
  The cohort's expected copies over those four passes are `[3.2030, 2.7970]`, `[3.2466, 2.7534]`,
  `[3.2606, 2.7394]`, `[3.2641, 2.7359]`.
- **The same three samples over three alleles**, each favouring one homozygote by the same 2
  nats: 2 passes, copies `[2.4868, 1.7566, 1.7566]`, summing to the six chromosomes three
  diploid samples carry.
- **Deleting the `/ cohort_chromosomes`** fails two tests. The division test stops at its first
  cell — one sample, 0.0018 raw copies against a `1e-3` threshold, so a locus that *has* settled
  reads as still moving; the thousand-sample cell at 1.8 raw copies would fail the same way. The
  cap test's three-sample locus goes from 4 passes to 6.
- **Moving the swap ahead of the E-step** scores every sample against a row just refilled with
  the sentinel. In debug it panics inside the prior (*"the cohort's expected allele copies …
  got `[NaN, NaN]`"*) — but that is a `debug_assert!`, so **under `--release` nothing panics**:
  the leave-one-out `max(0, ·)` absorbs the `NaN`, every sample is scored against the bare seed,
  and the loop reports `passes: 2, converged: true` where the shipped order gives
  `passes: 4, converged: false`. A converged flag, two passes early, with the cohort's evidence
  silently absent. `the_loop_reproduces_the_passes_driven_by_hand` catches it in both profiles.
- **Hoisting the inbreeding coefficient** to `inbreeding_by_sample[0]` moves the cohort's copies
  of the alternative allele by 0.38 out of six chromosomes — an allele-frequency shift of
  0.064 — and the pass count from 3 to 4.
- **A `ploidy` argument of 64 against a diploid table** reported `passes: 2, converged: true`
  where the truth is `passes: 4`, identically in debug and release, with nothing asserting.
  That argument no longer exists.

### The release-held checks

C2 has five, and the module's no-`Result` design rests on them holding in release. Downgraded
all five to `debug_assert` in one run and re-run under `--release`: `567 passed; 7 failed`, with
every check reached:

| check | tests that fail without it |
|---|---|
| the previous row names at least one allele | `a_convergence_test_over_no_alleles_is_refused` |
| the two rows are the same length | `cohort_rows_of_different_lengths_are_refused` |
| the cohort holds at least one sample | `a_convergence_test_over_a_cohort_of_no_samples_is_refused` |
| the threshold is a finite positive fraction | `a_threshold_that_is_not_a_fraction_is_refused`, `an_infinite_threshold_is_refused` |
| one inbreeding coefficient per sample | `fewer_…`, `more_inbreeding_coefficients_than_samples_are_refused` |

**And one check was deleted rather than tested**, which is the better outcome: the chromosome
count is no longer an `f64` argument, so an infinite or `NaN` one is not expressible and the
guard that refused it is gone. What is left is a `Ploidy` — at least one by construction — and a
sample count.

## 5. Validation

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4620 passed; 0 failed; 14 ignored`. Before C2: **4,603**, so C2 adds
  seventeen.
- `cargo test --release --lib ng::calling --all-features` — `574 passed; 0 failed; 3 ignored`.
  Before C2: **557**.

## 6. Trade-offs and follow-ups

- **`run_frequency_loop` has no caller outside tests.** C3 is the caller: it adds the final pass
  that reads each sample's genotype and quality off the posterior row, and turns
  `FrequencyLoopOutcome` into a `LocusInference`. D1 then assembles the two inert outer rounds
  around it.
- **A sample the candidate step set aside is still not handled**, and C2 does not change that.
  `sum_cohort_expected_copies`' doc comment already records the choice D1 owns: skip such rows
  in the sum, or never give them one. What C2 adds to the question is that the *denominator* is
  affected too — the stopping rule's chromosome count is `ploidy × the prepared sample count`,
  so a set-aside sample makes the convergence criterion marginally looser rather than being
  absent from it.
- **The threshold's floor is `1e-300` in the tests, not zero** — and it does not mean *never
  settles*, which is what the test helper using it was first called. Expectation-maximization
  reaches bitwise-identical passes: the three-sample fixture does so at pass 29 and reports
  `converged: true` however small the threshold. The helper is honest below that cap and a trap
  above it, and its name and doc comment now say so.
- **`FrequencyLoopOutcome::passes` stays a `u32`, and the review argued for `NonZeroU32`.** Its
  consumer is `LocusInference::passes`, a `u32` whose constructor asserts `passes > 0` for
  callers other than this loop. Retyping only this end would leave two adjacent types in the
  same plan disagreeing and would remove no check. **C3 is where the two meet and where the
  question belongs.**
- **The two open questions carried into C2 are unchanged by it**: whether a silent sample's
  flat-pass vote needs changing, and whether any locus shape traps the loop permanently. The
  second is now cheaper to answer — `run_frequency_loop` is the harness it needs, and the
  three-sample fixture's fixed point at pass 29 is the first datum: expectation-maximization on
  this model does reach a bitwise fixed point, so a permanent trap would have to be a cycle
  rather than a slow crawl.
