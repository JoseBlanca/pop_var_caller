# ng — calling prerequisites, A1: `InbreedingF` becomes half-open `[0, 1)`

**2026-08-23**, branch `ng-calling-prerequisites`. Step A1 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §7 and
[`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §2.1.

**One type, one constructor, one new error variant.** Nothing calls the genotype prior yet; this
step only makes the value it will read unable to carry the one number that would silence
heterozygotes.

---

## 1. What changed and why

`InbreedingF::try_new` accepted `1.0`. The genotype prior multiplies its heterozygote branch by
`1 − F`, so `F = 1` makes every heterozygote impossible: no read evidence, however clean, could
then produce a heterozygous call. The model itself survives that limit — production pins it in
`dirichlet_prior_full_inbreeding_concentrates_on_homozygotes`
([`posterior_engine.rs:4341`](../../../../src/var_calling/posterior_engine.rs)) — and what is
capped is the *estimate*: production's estimator clamps at `0.99`
([`inbreeding.rs:25`](../../../../src/paralog/inbreeding.rs)), while its engine config still
accepts `1.0` ([`posterior_engine.rs:4348`](../../../../src/var_calling/posterior_engine.rs)).
Spec §7 asks ng to move that ceiling into the type so a second estimator cannot forget it.

## 2. Changes made

**[`src/ng/types.rs`](../../../../src/ng/types.rs)**

- `InbreedingF::try_new` reuses `checked_probability` for the `[0, 1]` part and rejects exactly
  `1.0` on top of it. `checked_probability`'s behaviour is unchanged — it is shared with
  `ErrorRate` and `GenotypeFrequency`, for which `1.0` is a real answer (arch §2.1's first
  bullet); only its doc comment gained a paragraph saying why the composition is not the drift it
  warns about.
- New `DomainError::InbreedingFAtCeiling(f64)`, immediately after `DomainError::InbreedingF`.
- The acceptance assertion `InbreedingF::try_new(1.0).unwrap().get() == 1.0` moved out of
  `the_constrained_rates_accept_both_endpoints` — renamed
  `each_constrained_rate_accepts_the_endpoints_of_its_own_range`, since one of the three no longer
  accepts both and a test's name is what a failure prints — and into the rejection list beside
  `1.5`, as the
  plan directs. In its place that test asserts that the very next `f64` down still constructs, and
  pins that it *is* the next one — so exactly one value is excluded and no more.
- Both rejection messages now name `[0, 1)`. The pre-existing one said `[0, 1]`, which was true
  of it in isolation and useless in practice: someone who mistypes `1.5`, reads "not a fraction in
  [0, 1]" and retries `1.0` is refused again. New test
  `neither_inbreeding_rejection_names_the_closed_range` pins the property rather than the wording —
  that neither message names the closed range, and that only the ceiling message explains what
  `F = 1` costs. It fetches both through `try_new`, so it also pins which variant the constructor
  picks.
- `expect("a coefficient in [0, 1]")` at
  [`accumulators.rs:1029`](../../../../src/ng/parameter_estimation/generic/accumulators.rs)
  named the old range.
- `the_constrained_rates_accept_exactly_the_probabilities_and_round_trip` (a proptest over the
  whole `f64` line) asserted one acceptance range for all three types. It now carries two: closed
  for the other two, half-open for `InbreedingF`. **Left alone it would have stayed green, not
  gone flaky** — measured, a million draws from its generator produced `1.0` zero times and came
  no closer than 2.6 in a million, so the half of the property that had become false was
  unreachable. That cuts both ways, and the review caught it: the corrected expectation was
  unreachable for exactly the same reason. The generator gained two boundary arms — the literal
  `1.0` and the `f64` immediately below it — each weighted 1 against the two broad arms' 20.

**[`examples/ng_inbreeding_resolution.rs`](../../../../examples/ng_inbreeding_resolution.rs)** —
its printed sanity line said `InbreedingF` accepts `[0, 1]`, and claimed every fitted value the
example reports was inside it. The example constructs the type once, from the literal `0.5`, and
discards it; no fitted value ever meets it. The line now says both of those things and no more.

**Four documents corrected.** Each was accurate before this step and is not after it:

- [`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §7 — the guarantee that no sample
  reaches the caller at `F = 1` holds for a fitted `F` only (see §6 below), and §7 now also states
  what the half-open range does and does not buy.
- [`arch/calling_priors.md`](../../ng/arch/calling_priors.md) — its reconciliation table said the
  tightening was still owed.
- [`arch/parameter_prepass_generic.md`](../../ng/arch/parameter_prepass_generic.md) — two places:
  the "three constrained rates" comment, whose argument (three types rather than one shared
  `Probability`) rested on the ranges being indistinguishable and now says so directly, and the
  `try_new` comment thirteen lines below it, which still named `[0, 1]`.
- [`impl_plan/parameter_prepass_generic.md`](../../ng/impl_plan/parameter_prepass_generic.md) step
  A2 — a completed step, so its text stands as the record of what was built then, with a
  parenthesis pointing at the later change.

**Left alone deliberately:** the pre-change line numbers in `arch/calling_priors.md` §2.1 and in
this plan's own step A1. Both are instructions describing the tree as it was; re-pointing them at
the post-change lines would attach "today it admits `1.0`" to code that does not.

**Two edits fall outside the regions the plan promised**, and the parallel branch should know:
`checked_probability`'s doc comment and the body of the shared three-rate proptest, both in
`types.rs`. The plan's Worktree section names only the `InbreedingF` block, its boundary test and
the new variant. Neither is likely to collide — `ng-calling-foundations` adds `AlleleId` (an
unconstrained `u16`) and `Phred` (a non-negative `f32`), and neither is a `[0, 1]` fraction, so
neither belongs in that proptest or in that doc comment.

## 3. The one choice the plan left open, and how it was taken

The plan and the architecture both say "a new `DomainError` variant that says so" without saying
whether the existing `DomainError::InbreedingF` survives. **It does, and the two split the work:**
a value outside the range keeps `InbreedingF`; exactly `1` gets `InbreedingFAtCeiling`.

**The reason is what the two rejections need said to them, not which message would be true.** One
variant naming `[0, 1)` would be perfectly true of every value — that was this report's first
argument for the split, and the review showed it was a false dichotomy. The real argument is that
`1.5` is a typo, for which naming the range is the whole answer, while `1.0` is a coherent request
the model refuses, for which the answer is what refusing it means. Putting that explanation on the
shared variant would show it to someone who typed `-0.5`, where it is noise.

Reasoning from truth rather than from the reader is what left the shared message naming `[0, 1]`
until the review caught it — a user who mistyped `1.5` would have been sent back to the shell with
`1.0`, which is also refused. Both messages now name `[0, 1)`.

This is an implementation choice inside the coder's latitude, recorded here rather than escalated.

## 4. Blast radius — every constructor site in the crate

`InbreedingF::try_new` has eight call sites outside `types.rs`'s own tests. Six pass literals of
`0.0`–`0.5`, one is a test helper whose callers all pass `0.0`–`0.4`, and none is affected:

| site | value |
|---|---|
| [`generic/fallback.rs:412`](../../../../src/ng/parameter_estimation/generic/fallback.rs) | `0.4` |
| [`generic/accumulators.rs:962`](../../../../src/ng/parameter_estimation/generic/accumulators.rs) | `0.1` |
| [`generic/accumulators.rs:1029`](../../../../src/ng/parameter_estimation/generic/accumulators.rs) | `0.3` |
| [`generic/estimate.rs:384`](../../../../src/ng/parameter_estimation/generic/estimate.rs) | test helper; every caller passes `0.0`–`0.4` |
| [`generic/real_alignments.rs:926`](../../../../src/ng/parameter_estimation/generic/real_alignments.rs) | `0.0` |
| [`generic/truth_anchors.rs:338`](../../../../src/ng/parameter_estimation/generic/truth_anchors.rs) | `0.0` |
| [`examples/ng_inbreeding_resolution.rs:701`](../../../../examples/ng_inbreeding_resolution.rs) | `0.5` |

The eighth is the fitted path,
[`generic/runs.rs:634`](../../../../src/ng/parameter_estimation/generic/runs.rs), which constructs
from a coverage-weighted posterior occupancy with `.expect(…)`. That occupancy can reach exactly
`1.0` on a fully homozygous sample, so after this step the `expect` is a panic on a legitimate
fit. **A2 is the clamp that closes it**, and the plan makes it a separate commit.

**It is closer than "in principle", and nothing in the suite would see it.** Each window's
posterior is `exp(a + b − logsumexp(a + b, …))`, which returns exactly `1.0` in `f64` as soon as
the other state is about 37 nats behind — routine for a window carrying hundreds of sites, not a
rounding accident — and the sum is clamped to `[0, 1]` at
[`runs.rs:1178`](../../../../src/ng/parameter_estimation/generic/runs.rs), so an overshoot is
*converted into* exactly `1.0` rather than staying above it. What stands between that and the
panic is the used-both-states check, which needs some window below `0.5`; no reviewer could
construct an input satisfying both at once, so this is not a reachability claim. But the whole
library suite is green, so no test reaches the boundary either way: the four fits over real
alignments are `#[ignore]`d for want of a BAM, and the synthetic fits never approach it.

No production code uses this type: `InbreedingF` appears only under `src/ng/` and in one example.
Production's own `F` (`src/paralog/inbreeding.rs`, `src/ssr/cohort/`) is a bare `f64` and is
untouched.

## 5. Validation

All in the dev container, on the committed tree.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo clippy --example ng_inbreeding_resolution -- -D warnings` | clean — the stated gates build neither examples nor rustdoc, and this step changes an example |
| `cargo test --lib` | **3,939 passed, 0 failed, 11 ignored**, 636.61 s |
| `cargo doc --no-deps` | 16 unresolved-link errors, 11 "redundant explicit link target" warnings |

**The red baseline on `main` (`ee62a518`), measured before any of this was written**, so a
pre-existing failure can be told from a new one: `fmt` clean, `clippy` clean, `cargo test --lib`
**3,938 passed, 0 failed, 11 ignored** in 597.33 s, and `cargo doc --no-deps` already failing with
**the same 16 errors and 11 warnings** — none of them in a file this step touches. So the test
count moved by exactly the one test added, and the new intra-doc link to production's
`MAX_INBREEDING_COEFFICIENT` resolves rather than becoming a seventeenth broken one.

### The review

Five agents, one per category, each in its own worktree and each handed the gate output above so
it spent its time reviewing instead of re-running: what breaks when a legal value becomes illegal;
the range predicate at every float a program can construct; every claim in the author's own prose,
re-measured; fidelity to the plan and the design; and the strength of the tests under mutation.
Between them they applied **58 mutations** to this step's code.

**Nothing in the repository broke.** `cargo test --lib ng::parameter_estimation` passes 722 tests
on the changed tree, including every path that fits `F`, and the range check behaves as documented
at every float value a program can construct — both zeros, both smallest subnormals, the values
either side of one, three NaN bit patterns, both infinities, `f64::MAX` and `f64::MIN`.

**What they found were claims, and every one is fixed where it was made:** the shared rejection
message naming the old range (§2), the report's parallel-branch count (§6), the argument in §3,
the overclaim about what `[0, 1)` buys (§6), and four documents left saying `[0, 1]` (§2).

**Three mutations survived the tests as first written, and all three are now closed:**

- reverting the property test to its single shared expectation — green even at 200,000 cases,
  because the generator never reaches `1.0`. Closed by the two boundary arms.
- loosening "the nearest `f64` below one" to any value below one. Closed by asserting the value's
  successor is `1.0`.
- rewriting the ceiling message to keep both of the strings the test looked for while inverting
  its advice — it told the user to *raise* `F`. Closed by dropping the assertion on how the range
  is worded and keeping the one that matters, which is that no rejection names the closed range.

## 6. Follow-ups

- **A2** — the fitted path must clamp before constructing.

- **Every estimator still owes its own cap, and the type does not supply one.** Excluding the
  endpoint removes the mathematical limit and nothing else: the largest value the type accepts
  leaves `1 − F = 2⁻⁵³`, about 160 on the Phred scale against every heterozygote, where two clean
  alternative bases at Q30 supply 60. Production's `0.99` is 20 Phred, which evidence can
  overcome. Recorded in `spec/calling_priors.md` §7 alongside the correction below.

- **A correction to `spec/calling_priors.md` §7, made in this commit.** It said production's
  `0.99` clamp means "no sample ever reaches the caller carrying a prior that has ruled
  heterozygotes out". True of a *fitted* `F` only: `--inbreeding-coefficient` admits the closed
  `[0, 1]` ([`parsers.rs:166`](../../../../src/pop_var_caller/cli/parsers.rs), pinned admitting
  `1.0` by the test at [`:392`](../../../../src/pop_var_caller/cli/parsers.rs)) and the value
  reaches the engine as typed ([`pipeline.rs:343`](../../../../src/var_calling/pipeline.rs)).

- **Done, in the merge commit that follows this one.** `main` moved while this step was being
  reviewed: the genotype-prior branch landed, so the sites below stopped being someone else's and
  became this branch's to fix. **The list here was measured before that merge and undercounts it
  in two ways** — the `exact_spectrum` helper is driven with `1.0` by six tests, not four, and a
  fifth file, `calling/genotype_prior/hardy_weinberg.rs`, has a limit test of its own that this
  survey missed. Ten tests failed on the merge, across five constructor sites; what each became is
  in the merge commit's message.

  The list as measured on the branch, kept because it is what the survey found:

  | site on `ng-calling-prior` | tests reaching it with `1.0` |
  |---|---|
  | its own copy of the `types.rs` boundary test | 1 — the same assertion this step moved |
  | `calling/genotype_prior/dirichlet_multinomial.rs:1063` | 1 — `the_seam_rules_out_heterozygotes_at_the_greatest_coefficient_the_newtype_accepts` |
  | `calling/genotype_prior/seed_generic.rs:2249` | 1 — `a_spectrum_no_pair_can_produce_is_marked_rather_than_answered` |
  | `calling/genotype_prior/seed_generic.rs:2434` | 1 — `a_fully_inbred_panel_whose_spectrum_holds_heterozygotes_still_returns_a_pair` |
  | `calling/genotype_prior/seed_generic.rs:1212`, the `exact_spectrum` test helper | 4 — the loops at `:1245`, `:1451`, `:1560`, `:1804` each include `1.0` |

  Three of the five say in their own doc comments that they pin today's behaviour and name the
  tightening; `seed_generic.rs:2249` does not. **A further site is not affected and is easy to
  miscount:** `dirichlet_multinomial.rs`'s `mixed_row_for` helper drives
  `fill_inbreeding_mixture_log_priors` with a bare `f64`, so the loop at `:1161` reaches `F = 1`
  without constructing the newtype at all — that is the mathematical-limit test spec §12 test 3
  asks for, and it keeps working unchanged.
