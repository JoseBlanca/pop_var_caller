# ng step 4, generic path — A1+A2+A3: the module tree and the vocabulary

**Date:** 2026-08-06
**Plan:** [parameter_prepass_generic.md](../../ng/impl_plan/parameter_prepass_generic.md), Milestone A steps A1, A2, A3
**Design authority:** [spec](../../ng/spec/parameter_prepass_generic.md) · [arch](../../ng/arch/parameter_prepass_generic.md) §2.1 · [shared framing](../../ng/spec/parameter_prepass.md) §3

---

## 1. Plan

Three of Milestone A's six steps, run as one loop of the plan-driven implementation
skill. They are bundled deliberately and the bundling is named here and in the commit:
A1 is a pure scaffold — seven files carrying module documentation and nothing a review
could bite on — and A2 and A3 are adjacent declarations of the same vocabulary, one
half shared and one half step-local. A4 is **not** in this bundle; the plan marks it
*own commit, do not bundle*, and it stays that way.

- **A1** — create `src/ng/parameter_estimation/` with `mod.rs`,
  `fitting/{mod.rs, mixture_weights.rs}` and
  `generic/{mod.rs, depth_and_alt_reads.rs, histogram.rs, runs.rs}`; wire the module
  into `ng/mod.rs`.
- **A2** — four constrained newtypes in `src/ng/types.rs`: `ErrorRate`,
  `GenotypeFrequency`, `InbreedingF`, `Ploidy`, each with a private field, a checked
  `try_new` returning `DomainError`, and `.get()`.
- **A3** — the step-local scalars: `WindowIndex`, `INBREEDING_WINDOW_BP`, the three
  `ERROR_RATE_LADDER_*_PHRED` constants and `error_rate_ladder()`.

## 2. Assumptions

Three choices the plan and architecture left open. None changes a design decision;
each is recorded because a later reader would otherwise have to re-derive it.

1. **Where the step-local scalars live.** Arch §2.1 says they "stay in
   `parameter_estimation/`" without naming a file. They are in `mod.rs`, the module's
   own surface, because `WindowIndex` is read by `generic/` and the ladder by
   `fitting/`, and a scalar read from both sub-units belongs above both.
2. **`error_rate_ladder()` builds rather than tabulates.** The three Phred constants
   are the single statement of the ladder's shape, so a stored table would be a second
   one that could disagree with them.
3. **The three `[0, 1]` rates reject non-finite values explicitly**, not only through
   the range check. The range check alone already rejects `NaN` and `+∞`; `is_finite`
   is what makes `-∞` reject by the same rule as the rest, so a rate arriving from a
   division by zero cannot enter a likelihood. This copies `FlatEmission::try_new`
   ([emission.rs:287](../../../../src/ng/alignment/emission.rs#L287)) rather than
   `MismatchFraction::try_new`, which relies on the range check alone.

## 3. Changes made

**New — [src/ng/parameter_estimation/](../../../../src/ng/parameter_estimation/)**

| file | what it holds now |
|---|---|
| `mod.rs` | the module's contract; `WindowIndex`, `INBREEDING_WINDOW_BP`, the ladder constants and `error_rate_ladder()` |
| `fitting/mod.rs` | why the mathematics is a folder: it is step 4's one swappable seam |
| `fitting/mixture_weights.rs` | the concave climb — documentation only, built in Milestone D |
| `generic/mod.rs` | why there are two accumulators and neither is derivable from the other |
| `generic/depth_and_alt_reads.rs` | the one place that decides what an alternative read is — Milestone C |
| `generic/histogram.rs` | the cell table — Milestone B |
| `generic/runs.rs` | the inbreeding coefficient's two-state chain — Milestone E |

The five files with no code yet carry their module documentation and no
`#[cfg(test)] mod tests` block. **A deviation from A1 as written**, which asks for one
per file: an empty test module is three lines that assert nothing, and the block lands
with the code it tests in the milestone that writes it. `mod.rs` has its block, because
A3 gave it something to test.

**Changed — [src/ng/types.rs](../../../../src/ng/types.rs)**

Four constrained newtypes, and three new `DomainError` variants to carry their
rejections (`GenotypeFrequency`, `InbreedingF`, `Ploidy`). `ErrorRate` reuses the
existing `DomainError::ErrorRate`, whose doc is widened to name both producers.

`Ploidy` derives `Ord` where the three rates derive only `PartialOrd`: it keys the
histogram and the emitted rate maps, and its `u8` has a total order where an `f64`
does not.

**Changed — [src/ng/mod.rs](../../../../src/ng/mod.rs)** — `pub mod
parameter_estimation;`, the four newtypes added to the re-export list, and one clause
in the module doc's roll of what has landed.

## 4. Tests added

In `types.rs` — five, all against the constructors, since that is the whole of what
these types do:

- `the_constrained_rates_accept_both_endpoints` — zero and one are real answers for
  all three rates, so a half-open check would reject valid data.
- `each_constrained_rate_rejects_out_of_range_with_its_own_variant` — the point of
  four types rather than one shared `Probability`: a message naming the wrong quantity
  would send a reader to the wrong fit.
- `the_constrained_rates_reject_nan_and_the_infinities` — assumption 3 above, pinned.
- `ploidy_rejects_zero_and_accepts_every_real_copy_number` — the one value the type
  exists to make unrepresentable.
- `ploidy_orders_by_copy_number` — a haploid region's cells must sort before a
  diploid's, since ploidy is part of both histogram keys.

In `parameter_estimation/mod.rs` — three:

- `the_error_rate_ladder_spans_phred_10_to_50_in_161_rungs` — the count and both
  endpoints. The plan's A3 oracle.
- `the_error_rate_ladder_rungs_are_a_constant_ratio_apart` — every adjacent pair
  differs by `10^0.025`, and the rates descend as the Phred ascends. This is what makes
  "one rung" a unit of distance the coupled fit can report movement in.
- `the_inbreeding_window_is_a_hundred_kilobases` — the runs model's noise floor is a
  function of how many windows a genome has, so the constant is load-bearing.

## 5. Validation results

All in the container (`./scripts/dev.sh`).

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib ng::types::` | 26 passed, 0 failed |
| `cargo test --lib ng::parameter_estimation` | 3 passed, 0 failed |
| `cargo test --all-targets --all-features` | 2,904 passed, **1 failed — pre-existing** |

**The one failure is not this step's**, and it was proved so rather than assumed:
`ng::locus_generation::pileup::parity::every_divergence_from_production_is_one_of_the_six_named_classes`
fails at seed `0x5eed0001` case 18 with `record_widen_events` 4 against production's 3.
Stashing every uncommitted change and re-running reproduces it identically at `HEAD`.
It belongs to locus generation and is reported to the owner rather than absorbed here.

Both gates were red on this branch before any of this landed — two whitespace spots and
seven lints in `examples/` — and were greened in their own commit (`ce3f0b4`) so that
this one carries only the plan's code.

## 6. Tradeoffs and follow-ups

- **Nothing computes yet.** Milestone A is types by design; the first thing that reads a
  number is Milestone D. The ladder is the exception and is tested against its own
  arithmetic.
- **`Phred` is not a type**, and will not become one. `types.rs` already has `LogProb`
  for the logarithm of a probability; a second log-scaled type in a different base would
  make a base mix-up a plausible wrong number instead of a compile error.
- **Deferred to A4:** the depth ladder, which the plan isolates in its own commit
  because its edges are a correctness parameter — sixteen bins at the same cap biases
  the fitted error rate by 0.55 rungs against 0.05 for twenty, and nothing downstream
  would show it.
