# The hidden-duplication filter — A2: the score, and the differential that proves it is production's

**Date:** 2026-09-06
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone A, step A2 — *own commit, do not bundle*
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.2, §7, §10
**Branch:** `ng-paralog-filter`

## The answer

**ng's copied scorer and production's return the same numbers on 800 randomised loci, bit for
bit, at cohort sizes 1, 2, 10 and 63.** The score is the filter's core question: at one locus,
how much better does *two reference-collapsed gene copies piling their reads onto one position*
explain every sample's depth and allele counts than *a real variant some of them carry* does. It
comes out as a log likelihood ratio — positive favours the collapsed pair.

The comparison is on every number the scorer returns, not the ratio alone: the ratio is a
difference of two log-likelihoods, so two errors that shifted both the same way would cancel in
it. Both log-likelihoods and both sample counts are asserted too, and **equality is by bit
pattern, not by tolerance** — a port that agrees to twelve decimal places is a port that has
diverged.

**One sample is in that sweep deliberately** (spec §4). The folded site-frequency-spectrum prior
runs over `[1/2N, 1 − 1/2N]`, which at `N = 1` is the single point `[½, ½ ]`, so the grid
degenerates; production has never been run there, because its filter is a cohort filter.
**144 of the 200 one-sample loci score, and every one of them agrees with production's bit for
bit.** The other 56 are loci whose only sample was drawn absent or with a degenerate σ₀ — those
are neutral by design, and both trees agree they are.

## What was added

| file | what it is |
|---|---|
| [`src/ng/paralog/locus_score.rs`](../../../src/ng/paralog/locus_score.rs) | production's `src/paralog/locus_score.rs`, 797 lines past the header, differing on **three declared lines** — see below |
| [`src/ng/paralog/production_parity.rs`](../../../src/ng/paralog/production_parity.rs) | ng's own: the differential, its generator, and the checks on the generator itself |

`mod.rs` declares and re-exports both; `copy_fidelity.rs` guards the copy.

**Three lines differ, all of them paths into production, all declared to the guard.**
Production's transcribed test module has `use crate::paralog::{GridSpec, SfsPriorSpec};`, which
left alone would build ng's fixture from *production's* types — and would not compile, since ng's
are distinct types. Two more are rustdoc link definitions naming
`crate::paralog::SingleCopyCoverageModel`; ng re-exports its own, so left alone they would send a
reader of ng's own API documentation into the frozen tree for the type ng's scorer consumes. **A
`use crate::` sweep finds neither of those two, which is part of why the first draft of this step
missed them** — the review did not.

Each is declared to the guard as a `Repoint`, and the sanction is exact: production's line must
occur exactly once, ng's must be exactly the declared replacement, the indentation and the line
ending must match, and a repoint may only turn a `crate::paralog` path into a `crate::ng::paralog`
one — so the mechanism cannot be used to sanction a changed constant. **One of the three is a
`//!` line inside production's module header**, which the guard compares as a byte-verbatim
prefix, so substitutions are applied to the whole of production's file rather than only to the
content past the header. **Measured: `diff` reports 6 lines, 3 pairs, and all three are
declared.** A1's guard claimed these files needed no repoints; that was true of the two files it
guarded and is now false, and the module header says so.

## Assumptions and recorded deviations

**1. The differential lives in ng, as a `#[cfg(test)]` module, and reads production.** This is
the pattern `src/ng/mod.rs` already names for its two other parity oracles (`scanner_parity`,
`calling::genotype_table_parity`): nothing shipped depends on production through it, and the
direction that matters — production depending on ng — stays absent.

**2. The generator is a fixed-seed splitmix64, not a crate.** `rand` is not a dependency, and the
codebase's idiom for a reproducible test stream is a hand-rolled generator with the reason stated
([`calling/parameters_file/to_toml.rs`](../../../src/ng/calling/parameters_file/to_toml.rs),
[`parameter_estimation/subsample.rs`](../../../src/ng/parameter_estimation/subsample.rs)). A
failure reproduces from its seed with nothing installed.

**3. The plan asks for "zero-read samples included"; the stream draws five shapes, not one.**
Each is a branch of the scorer's admission rules, and a differential that missed one would be
comparing two implementations on a narrower input class than it claims. Measured on the
63-sample stream's **12,600 drawn observations**:

| shape | drawn | drawn at |
|---|---|---|
| absent from the locus | 1,598 | 1 in 8 |
| zero reads — the coverage-only case spec §3.2 rests on | 2,210 | 1 in 6 |
| `alt_reads > total_reads`, which the scorer clamps | 1,256 | 1 in 10 |
| degenerate σ₀ (zero, negative, `NaN`, infinite), which drops the sample | 1,841 | 1 in 7 |
| copy number past the winsor cap of four, which it clips | 6,990 | 5 in 9 |

**4. The σ₀ rate was wrong when first written, and measuring it is what found that.** The first
draft drew σ₀ from four cases in seven rather than four in twenty-eight. Measured on that draft's
63-sample stream: **7,230 of 12,600 samples were degenerate — 57 in 100** — against the 1 in 7 its
own comment claimed. That matters beyond the comment. Combined with the 1-in-8 absent rate it
leaves `1 − (7/8)(3/7) = 62 in 100` one-sample loci with nothing usable, and a locus with nothing
usable scores the *neutral* verdict — which two implementations agree on while computing nothing.
The differential would have been passing largely on agreement about zero.

**5. So the sweep asserts what it reached, not just that it agreed.** `samples_used > 0` per
locus is counted and pinned exactly — **144, 189, 200, 200 of 200** at cohort sizes 1, 2, 10 and
63. Pinned rather than thresholded because the stream is deterministic: a later change to the
generator that narrows the differential's reach then has to be looked at instead of absorbed. At
one sample the expected share of unscored loci is `1 − (7/8)(6/7) = 25 in 100`, against the
**28 in 100** measured.

**Counting finite ratios would not have worked**, and that was the first attempt: the neutral
verdict is `0.0`, which *is* finite, so a stream that scored nothing at all would have shown 800
finite ratios out of 800.

## Tests

**59 under `ng::paralog::`**, up from 31 at A1: 23 in `coverage_model`, **18 in `locus_score`** and
4 in `model_params`, all production's transcribed; and 14 of ng's own — 5 in `copy_fidelity`, 1 in
`mod`, and **8 in `production_parity`**:

- `the_copied_scorer_agrees_with_productions_bit_for_bit` — the 800-locus sweep, with the loci that
  actually admitted a sample pinned per cohort size at 144, 189, 200 and 200 of 200.
- `one_sample_scores_finitely_and_identically` — 200 more loci at `N = 1`, asserting that every
  locus whose single sample was usable was actually scored with it, not merely that the ratio came
  back finite.
- `a_negative_copy_number_is_winsorised_at_zero_on_both_sides` — the winsor clamp's *lower* arm.
- `a_mismatched_sigma_slice_is_neutral_on_both_sides` — spec §6 trap 3, a slice both shorter and
  **longer** than the cohort.
- `tables_built_for_another_cohort_are_neutral_on_both_sides` — the per-pass tables built for a
  different `N` than the locus carries, which is the mismatch step B1's run wiring is most likely
  to produce.
- `an_empty_carrier_set_is_neutral_on_both_sides` — the one input that makes a *scored* locus
  return `NaN` if the scorer stops refusing it, which spec §6 trap 4 forbids.
- `a_locus_with_no_usable_sample_is_neutral_on_both_sides`.
- `the_drawn_stream_contains_every_shape_the_differential_claims` — the generator's own check,
  with all six shape counts pinned exactly.

**The last four of those eight came out of the review, and each kills a mutation that the
differential as first written did not.** With `production_parity` run alone, so that the textual
guard cannot mask the result: narrowing ng's `clamp(0.0, cmax)` to `min(cmax)`; deleting the
`precompute.cohort_size != cohort_size` clause; weakening the σ₀ length check from `!=` to `<`;
deleting the `configs.is_empty()` clause. All four passed the first version of the differential and
all four fail this one.

**The differential's power was measured, not assumed.** Perturbing ng's winsor clamp from
`clamp(0.0, cmax)` to `clamp(0.0, cmax * 1.0000001)` — a change of 1 part in 10 million to one
input — makes it fail at `cohort size 1, locus 0` with `ng 85.12238217765987 (0x405547d51c0eaa65)`
against `src/paralog/ 85.12233325249034 (0x405547d44ed9aa32)`: the two differ by **4.9 × 10⁻⁵ in a
value of 85.1, about 1 part in 1.7 million**, which no tolerance-based comparison would flag. The
mutation was reverted from a backup and the restoration confirmed by `diff`.

**The guard's repoint machinery was mutation-tested too.** Changing ng's repointed import to an
unsanctioned but compiling variant (`{SfsPriorSpec, GridSpec}`, reordered) fails at
`locus_score.rs, line 743 of the copied content`; so does deleting the `Repoint` declaration while
ng's line stays repointed. Only the exact declared substitution passes, and every rejection is now
driven on synthetic strings rather than only on the real, currently-agreeing pair.

## Validation

Run in the dev container, on this tree.

    cargo test --all-features --lib "ng::paralog::"
    → test result: ok. 59 passed; 0 failed; 0 ignored; 6290 filtered out

    cargo test --all-features --lib --bins --tests
    → 6334 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 3 errors, all needless_lifetimes in cohort_merge/build.rs and serial.rs

    cargo fmt --check
    → dirty on 9 files, none of them src/ng/paralog/

The integration failure, the three lints and the nine files are `main`'s at `a33ada0f`, proved so
before A1 began by removing `pub mod paralog;` and re-running (6,275 lib tests, same single
failure). The gate — no failure `main` does not already have — holds, and the lib count moves
6,306 → 6,334.

## What this closes and what it does not

**Spec §9's second OPEN — `N = 1` — is closed at the level this step can close it.** The copied
precompute runs at one sample, produces finite ratios, and agrees with production's bit for bit
on 144 scored loci out of 200 drawn. What it does *not* say is anything about real data: plan
step D1 runs one tomato accession and reports what the filter actually does there.

**Spec §9's first OPEN — every record scored, non-SNPs on coverage alone — is untouched.** It is
step C1's, and the plan requires it confirmed with the owner before C1 is coded.

## Follow-ups

- **`GuardedCopy` does not fit step A3's copy, and that is A3's first task.** The plan puts two
  originals into one ng file — `src/paralog/prior.rs` and `calibrate.rs:64-118`, a line range in
  `src/var_calling/` — where the type holds one source, `CopyBegins` has no end marker, and three
  failure messages name `src/paralog/` as a literal. Raised at Checkpoint A rather than pre-built
  here.
- **Nothing constrains what a `Repoint` may say in general.** A repoint must now turn a
  `crate::paralog` path into a `crate::ng::paralog` one, which closes the case at hand; the
  differential would catch a repoint that changed a constant in `locus_score.rs`, and A3 brings a
  file the differential does not cover.
- **The differential runs at one parameter setting and at no cohort size near spec §4's three
  thousand.** Both are size questions rather than correctness ones and belong to the D milestone.
- **Five equivalent mutants point at dead code inside production's frozen scorer** — the
  `total_reads > 0` clause in the hom-alt guard never discriminates, `ms.sort_unstable()` never
  reorders, the `usable.is_empty()` early return is a performance guard only, and `LogSumExp`'s
  `>` against `>=` is exactly equivalent. Production's backlog, not this plan's.
