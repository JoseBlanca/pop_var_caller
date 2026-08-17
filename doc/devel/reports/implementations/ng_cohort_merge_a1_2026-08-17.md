# ng cohort merge — A1: the three parameters

*Implementation report, 2026-08-17. Step A1 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §3, §4.3, §6.1, §6.4 and
[arch](../../ng/arch/cohort_merge.md) §1.*

## 1. Plan

Land the module home `src/ng/run/cohort_merge/` and the three parameters a calling run
sets, each as a `NonZeroU32` newtype with its default: `MaxCohortLocusSpan` (50),
`MinAltObs` (2), `CohortLocusBuilderRegionsLen` (20). No behaviour — the walk that
reads them is A3 and A4.

## 2. Assumptions and recorded deviations

- **The module is declared `pub`, not `pub(crate)`.** The arch calls everything in it
  "crate-private machinery inside the two caller objects", and it will be, but those
  objects do not exist yet: `pub(crate)` items with no consumer are `dead_code`, which
  is denied on this project's clippy gate. Every other ng module is `pub mod` for the
  same reason. The narrowing belongs with the caller objects (run arch §3).
- **`DEFAULT_MIN_ALT_OBS` carries production's value, 2, by copy rather than by name.**
  The arch cites `var_calling/mod.rs:72`; the constant there is
  `DEFAULT_MIN_ALT_OBS_PER_SAMPLE`, the default of `--min-alt-obs-per-sample`, and
  production feeds that one number to two different tests — the cohort keep rule
  ng inherits (`derive_is_kept`, `cohort_integration.rs:166`) and a per-sample pre-EM
  filter (`variant_caller.rs:380`). ng's threshold is only the first, so reaching the
  constant by name would let a retune of production's per-sample filter move ng's
  cohort keep silently. The provenance is in the doc comment instead.
- **The arch declares `DEFAULT_MAX_COHORT_LOCUS_SPAN` twice** (§1, once at each of two
  adjacent doc comments). One declaration, no design content lost.
- **`ObservationReachCeiling`, the fourth type in arch §1, is not here.** A1's contract
  names three parameters; the ceiling is the observation cache's business (spec §6.4)
  and lands with it in milestone D.

## 2a. What the review changed

The step went through `rust-code-review` (six categories) and `apply-code-review-fixes`
before it was committed — [review](../reviews/ng_cohort_merge_a1_2026-08-17.md),
[fixes](../reviews/fixes_applied_ng_cohort_merge_a1_2026-08-17.md). Three things below
are the review's rather than the first draft's, and they matter enough to name here:

- the **`pub const DEFAULT: Self`** on each newtype. The first draft made a zeroed
  default a build error with three `const _: () = assert!` lines paired by hand with the
  three constants — a fourth parameter that forgot its line would have compiled and
  panicked at run time. The associated constant makes the guarantee at the call, so
  nothing has to be remembered;
- the test **`get_returns_the_wrapped_value_not_the_default`**. Every `get()` call went
  through `Default`, so an accessor ignoring its argument passed the whole suite while
  silently running every cohort at 50/2/20 whatever the operator set;
- the **`DEFAULT_MIN_ALT_OBS` doc**, which claimed production's measured justification
  for a rule that is not production's (below).

## 3. Changes made

- **`src/ng/run/mod.rs`** (new) — the calling run's module home, one line of tree.
- **`src/ng/run/cohort_merge/mod.rs`** (new) — the three newtypes, each with a typed
  `pub const DEFAULT: Self`, a `Default` returning it, a `const fn get()` giving the
  inner `u32` (the shape `ContigId` and `Position` already use in `ng/types.rs`), and
  the three readable `DEFAULT_*` `u32` constants.
- **`src/ng/mod.rs`** — `pub mod run;` and one clause in the module doc's landed list.

The two spans are separate types on purpose (spec §11): `MaxCohortLocusSpan` judges a
locus, `CohortLocusBuilderRegionsLen` sizes a builder's work, and handing one where the
other belongs is a compile error rather than a wrong number.

**A zeroed default is a build error, and that is measured, not asserted.** Each
`DEFAULT` is a `const` item calling a `const fn` whose `None` arm panics, so the arm is
evaluated when the crate is compiled. Setting `DEFAULT_MIN_ALT_OBS` to 0 and running
`./scripts/dev.sh cargo build --lib` gives:

```
error[E0080]: evaluation panicked: a cohort-merge default must be non-zero
   --> src/ng/run/cohort_merge/mod.rs:100:36
    |
100 |     pub const DEFAULT: Self = Self(non_zero_default(DEFAULT_MIN_ALT_OBS));
    |                                    ^^^^^ evaluation of `…::MinAltObs::DEFAULT` failed inside this call
```

The constants stay `u32` at their declaration, where an operator reading the source
looks for them.

## 4. Tests added

Both in `src/ng/run/cohort_merge/mod.rs`:

- `the_defaults_are_the_documented_values` — the three `pub` constants *and* the three
  typed defaults, each pinned to 50 / 2 / 20. Both spellings, because a command line
  advertises one and a run uses the other.
- `get_returns_the_wrapped_value_not_the_default` — each newtype built from a value that
  is not its default (200, 7, 100) and read back, plus `NonZeroU32::MIN` and
  `NonZeroU32::MAX`. Without it, an accessor answering with its own default passed the
  whole suite while ignoring what the operator set.

Net zero tests: the review deleted one that could not fail
(`the_region_width_is_not_the_locus_bound` — see the fixes report) and added one that
can.

## 5. Validation

All in the container (`./scripts/dev.sh`), on the tree as committed:

| command | result |
|---|---|
| `cargo fmt --check` | clean, exit 0 |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib` | see the commit message for the count; the two tests of this step run as `ng::run::cohort_merge::tests::*` |
| `cargo test --tests` | 80 integration tests across 9 targets, 0 failed |
| `cargo build --lib` with `DEFAULT_MIN_ALT_OBS = 0` | `error[E0080]` as quoted in §3 — the guarantee, exercised |

**The `--all-targets --all-features` clippy gate is red and was red before this work**:
49 errors across 20 files in `examples/`, `benches/`, and the test code of
`census_file.rs`, `ssr_fit.rs`, `open_bam.rs` and `src/ssr/cohort/sim.rs` — the last of
which ng may not edit. Confirmed pre-existing by stashing this step's changes and
re-running. None are in this step's files, and greening them is its own commit.
`cargo fmt` did green the formatting gate, which needed three unrelated files
(`benches/ng_joint_fit_perf.rs`, `examples/dhat_ng_joint_fit.rs`,
`src/ng/locus_generation/ssr.rs`) and is included here.

## 6. Tradeoffs and follow-ups

- Nothing reads these parameters yet. A3 reads the first, A4 the second; the third is
  read by the organiser in milestone E.
- No CLI wiring: the parameters belong to a calling run's command line, and the run
  objects that own a command line are out of this plan's scope
  ([plan](../../ng/impl_plan/cohort_merge.md), *Out*).
- **Owed by the design, not by this step:** `MaxCohortLocusSpan`'s effective value has
  to reach the run's output beside the failed-locus count, or two runs over the same
  records under different bounds are indistinguishable (arch §1; spec §3.1, §3.3). The
  doc comment carries the obligation; the emission step owns the surface.
- **Five defects in the design documents** were surfaced by the review and are the
  owner's to rule on — two of them reach A4. They are listed as open questions in
  [the review](../reviews/ng_cohort_merge_a1_2026-08-17.md) §4 and raised at
  Checkpoint A.
