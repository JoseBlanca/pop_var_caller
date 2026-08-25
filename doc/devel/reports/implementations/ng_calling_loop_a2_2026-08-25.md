# ng calling loop — A2: the calling seam, and every switch as a value

**Date:** 2026-08-25
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step A2
**Design authority:** [spec/calling_em_loop.md](../../ng/spec/calling_em_loop.md) §4.1, §5.1, §6;
[arch/calling_em_loop.md](../../ng/arch/calling_em_loop.md) §2.1, §3.1, §6.1, §6.2
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`

> **Read this against the review that followed it.** Five category reviews raised three Majors
> and thirteen Minors, and the fixes changed some of what §2 and §3 describe: three counts
> became `NonZeroU32`, the discovery bar became the merge's own `MinAltReads`, the seam now
> takes a validated token rather than a bare configuration, and the error enum shrank from
> eight variants to five because three of the values it refused stopped being expressible.
> [The review](../reviews/ng_calling_loop_a2_2026-08-25.md), and
> [what was done about it](../reviews/fixes_applied_2026-08-25_v2.md).

---

## 1. Plan

A new module, [src/ng/calling/inference/mod.rs](../../../../src/ng/calling/inference/mod.rs),
holding two things and no arithmetic:

- **`LocusGenotyper`** — the one boundary every way of handling a cohort crosses: the
  evidence, the frozen parameters, the candidate alleles, the configuration and the worker's
  scratch in; a called locus out.
- **`CallingLoopConfig`** — every switch of the design as a **value**, including the two
  outer loops that ship switched off, with the constants they inherit named and marked soft.

Plus the rule the plan is explicit about: **a setting whose body is not built is refused,
never silently ignored.**

## 2. Assumptions and deviations

Five. Four are consequences of what A1 built or of what the sibling modules shipped; the
fifth is a number this step had to choose.

### 2.1 The seam carries a type parameter the architecture's sketch does not

Architecture §3.1 writes `scratch: &mut CallingScratch`. A1's `CallingScratch` is generic
over the repeat-tract emission model's own working memory, because the `SsrRowScratch` it
holds is — so the trait is `LocusGenotyper<SsrEmissionScratch>` and an implementation is
generic over it. The loop hands that scratch straight to the row builder and never looks
inside it, so no implementation is specialised by it.

**Object safety survives**, which is the property that matters: `Box<dyn
LocusGenotyper<StutterSubstitutionScratch>>` is what a run holds when it selects an arm, and
a test pins it.

**An earlier version of this section said an associated type "would have worked equally". It
does not compile** (measured by the review's `idiomatic` agent): an implementation that works
for every model constrains the type nowhere, which is `error[E0207]`, and every implementor
would need its own type parameter plus a `PhantomData`. The plain parameter is forced rather
than preferred.

### 2.2 The trait gains `name()`

Not in the architecture. Taken by analogy with
[`GenotypePriorModel::name`](../../../../src/ng/calling/genotype_prior/mod.rs), whose own
doc gives the reason and it applies here word for word: **this seam exists to compare three
answers, and a result that cannot say which one produced it is not auditable.** Spec §12's
first open question is a three-arm comparison; an arm without a label is a number nobody can
act on.

### 2.3 `validate()` returns a `Result`, in a module that otherwise has none

Everything else in `calling/` panics, because everything else it refuses is a wiring mistake
(spec §8). **A configuration is not a wiring mistake — it is a run's request**, and the
honest answer to a request this caller cannot serve is to say so and stop, not to assert.
So there is a typed `CallingLoopConfigError`, and it is the only `Result` in the folder.

### 2.4 `SummariseConditionLoop` is not declared here

Architecture §3.1 declares it beside the trait. Its body is the plan's step D1, in
`summarise_condition.rs`, and a type shipped now would carry either an `unimplemented!()` or
a loop the plan has not reached. The seam is exercised instead by a **stand-in in the test
module** that calls every sample homozygous reference and looks at no evidence at all —
which is the pattern `genotype_prior` already uses for the same purpose, and whose own
comment says why: it is deliberately not a caller, so nothing can be mistaken for a check of
one.

### 2.5 One constant is this step's own, and nothing inherited it

`DEFAULT_DISCOVERY_MAX_ROUNDS = 4`. The spec stops a discovery loop on a round that adds
nothing, or on the allele cap, and names **no round cap**; the architecture gives the field
and no value. Four is twice the expected one or two rounds — a runaway guard, not a
measurement — and it is marked soft in place, inert while discovery ships off, and left for
the plan that switches discovery on to set.

**Its placement is the part that is not arbitrary.** `DiscoveryConfig::default()` is written
out rather than derived, so the guard is 4 and not 0. Derived, every field but the mode would
default to nothing, and a run switching the mode on would get a discovery loop that runs no
rounds and reports finding nothing — a setting silently doing nothing, which is exactly what
`validate` exists to prevent one level up.

## 3. Changes made

| what | shape |
|---|---|
| `LocusGenotyper<SsrEmissionScratch>` | `call_locus(evidence, parameters, candidates, config, scratch) -> LocusInference`, plus `name()`. `candidates` by value — a discovery round appends to the table and the prune shrinks it, so the loop owns it and hands it back inside the result |
| `CallingLoopConfig` | `convergence_threshold`, `max_passes`, `slippage_refit`, `discovery`, with `Default` and `validate()` |
| `SlippageRefitConfig` | `max_rounds` (0 = frozen), the two pull-backs, the round threshold; `is_frozen()` |
| `DiscoveryConfig` / `DiscoveryMode` / `DiscoveryBar` | `Off` / `AgainstFrozenFrequencies` / `AgainstFullConvergence`; the bar's two halves; `is_off()` |
| `CallingLoopConfigError` | `#[non_exhaustive]`, `thiserror` — eight variants as written, five after the review, because three of the values it refused stopped being expressible |
| nine `DEFAULT_*` constants | each with what it is, where it came from, and whether anything has measured it |

**The two not-yet-built settings are refused with a message that says what accepting them
would cost.** Not "unsupported": *"it ships frozen at zero rounds, and accepting this setting
would run the frozen loop and report it as the re-fitted arm"*. That is the failure being
guarded — a measurement harness that sets the re-fit on, gets the frozen loop's answers, and
reports the two arms as agreeing exactly.

**Range checks come before not-built checks**, so a configuration that is both is told about
the range first: that half is the caller's to fix today, where the other is this caller's to
build.

**The allele cap is absent, deliberately.** It belongs to candidate selection, which holds it
as a `MaxCandidateAlleles` — a type refusing anything below two. A `u16` field of the same
name here would be two spellings of one rule and would drop the check (arch §2.1). The
module's own doc says so, so that the absence reads as a decision rather than an omission.

## 4. Tests added

Six.

| test | what a wrong implementation would do |
|---|---|
| `the_shipped_configuration_is_one_this_caller_will_run` | a default and its validation drift apart — the shipped configuration is the one nobody passes explicitly, so it is the one a range check is least likely to be tried against. It also pins that `Off` is the *mode* and not a zeroed guard |
| `a_setting_whose_body_is_not_built_is_refused_rather_than_ignored` | the frozen loop's answers reported as the re-fitted arm's, so the two arms agree exactly and that reads as a finding |
| `a_value_outside_its_range_is_refused_by_the_check_it_belongs_to` | a `NaN` threshold passing a `<= 0.0` comparison, after which the loop's stopping test can never be satisfied and every locus runs to its cap |
| `zero_pull_back_is_the_free_setting_and_passes_the_range_check` | a range check that refuses zero would make one of the three settings the design compares unreachable |
| `an_out_of_range_value_outranks_a_setting_that_is_not_built` | the caller told to wait for a feature when what is wrong is a value they can fix now |
| `the_seam_is_object_safe_and_every_arm_can_name_itself` | a trait that cannot be held behind a `Box`, which would make choosing an arm a compile-time decision and the three-arm comparison need three binaries |

## 5. Validation

In the container, from this worktree's own `scripts/dev.sh`, at main's 1.98 compiler pin.

| command | result |
|---|---|
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 — **the wider scope**, which this branch can now run since taking main's pin |
| `cargo test --lib` | **4,523 passed, 0 failed, 14 ignored** |
| `cargo test --release --lib ng::calling::inference` | 6 passed, 0 failed |

A1 left the library at 4,517, so the six tests above are the whole of the difference.

## 6. Trade-offs and follow-ups

- **The `Result` is a one-off in this folder and that is a seam worth watching.** If a second
  configuration type appears, the question is whether it shares this error type or gets its
  own; today one enum with eight variants is smaller than two with four.
- **Two of the eight variants are temporary.** `SlippageRefitNotBuilt` and `DiscoveryNotBuilt`
  come out when the bodies land, and `#[non_exhaustive]` is what lets them go without breaking
  a match.
- **A discovery mode that is on with a zero round cap is not checked**, because any non-`Off`
  mode is refused before that could be reached. The check belongs with the plan that removes
  the refusal.
- **`SummariseConditionLoop` and the two assignment arms are unwritten.** The first is step
  D1; the other two, with their joint priors and the exhaustive scorer, are the bake-offs
  plan's.
