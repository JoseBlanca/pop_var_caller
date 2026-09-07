# The hidden-duplication filter — A3: how common duplications are in this run, and where to cut

**Date:** 2026-09-06
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone A, step A3
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.3, §7
**Branch:** `ng-paralog-filter`

## The answer

**Milestone A is complete: all four of the filter's statistics are in ng, and each agrees with
production's bit for bit.** A2 brought the per-locus score — how much better a collapsed pair of
gene copies explains a locus than a real variant does. This brings the two that turn that score
into a verdict:

- **the prior** — how common hidden duplications are *in this run*, fitted from the run's own
  scores by an EM rather than assumed;
- **the false-discovery curve and its cut** — with the prior in hand, each score becomes a
  probability, and the operator's target false-discovery rate resolves to a likelihood-ratio
  threshold. The default target is about 1 record in 100 of those removed being a real variant.

## What was added

| file | what it is |
|---|---|
| [`src/ng/paralog/prior.rs`](../../../src/ng/paralog/prior.rs) | production's `src/paralog/prior.rs`, **524 lines past the header, differing on none** |
| [`src/ng/paralog/calibration.rs`](../../../src/ng/paralog/calibration.rs) | **a span** of `src/var_calling/paralog_filter/calibrate.rs` — 83 lines, from its first item down to but not including the cohort inbreeding coefficient — differing on **seven declared lines** |

**One copy takes a span rather than a whole file, and that is new.** Production keeps
`ParalogCalibration` in the same file as the spill-streaming driver ng does not want, so ng takes
the file from its first item to the line before `cohort_inbreeding` — and **that helper is
deliberately left behind**: production replaced its per-individual inbreeding estimate with one
cohort number because the proxy it had was divergence-contaminated, while ng's parameters file
already carries a fitted coefficient per sample (spec §3.2). It is the one piece of the filter ng
replaces rather than ports.

### The seven lines, and the three kinds of change now sanctioned

Two are what A2 established: a path into production, and a documentation link to a production item
ng does not carry (`[`calibrate`]`, the driver left behind, which would otherwise be a link
resolving to nothing).

**Five are a third kind, new here: `pub(crate)` widened to `pub`.** Production keeps these items
inside a private module; ng's `paralog` is a public module and the run wiring that will use them
arrives at step B1, outside it. Left crate-private they are dead code — `cargo clippy -D warnings`
says so in four places — and the module would export a calibration nothing can name.

Each kind is a variant of `WhyRepointed`, and each carries its own check, so the mechanism cannot
be used to wave through a change of behaviour. The visibility check is the tightest: **the two
lines must be identical once `pub(crate)` is replaced by `pub`**, and nothing else on the line may
move. Verified by mutation: declaring a repoint that widens the visibility *and* changes the
fallback prior from `0.03` to `0.05` is refused with the reason named.

## The guard, extended for a span

`GuardedCopy` gains two fields: `production_path`, so a failure message names the file that must
not be edited rather than a directory — not every original is under `src/paralog/` any more — and
`ends_before`, the first line *past* the copy.

**A declared end marker that production no longer has is a guard failure, not a pass.** Without
that, a renamed marker would let the extraction run to the end of production's file and the
comparison would fail somewhere unrelated, with a message blaming the copy.

**One boundary is normalised, and it is worth saying which.** Production separates its items with
blank lines, so a span ends with one; ng's copy cannot keep it, because `cargo fmt` strips a
trailing blank line and the copy would then fail on every run. So a span's end is normalised on
both sides to no trailing blank lines and one final newline. It is the only place either side is
normalised at all, and it is a boundary rather than content. **Measured**: with the blank line
kept, the guard passed and `cargo fmt` then broke it; with the normalisation, both hold.

## Assumptions and recorded deviations

**1. The differential covers the prior and the curve too, which the plan does not ask for.**
A3's stated proof is "production's tests, transcribed". A2's review made the case against relying
on that alone — transcribed tests pass on both trees whatever either computes, and *"once a file is
released the differential is the only guard left"* — so `production_parity.rs` grows two tests
covering the other three copied quantities. It folds one randomised stream of likelihood ratios
into both trees' histograms and compares π, whether the EM converged, the curve at fifteen probe
ratios, and the resolved cut at five targets, all by bit pattern.

**The streams are four named shapes, not one uniform draw**, because the EM's answer depends on
the *shape* of the distribution rather than its spread: a run where nothing is duplicated, one
where about a tenth are, one where everything is, and one with `NaN` and infinities mixed in.

**2. `ParalogCalibration`'s own tests are ng's, because production's are not portable.** They live
in `calibrate.rs`'s test module, outside the copied span, and depend on production's spill and
pre-pass machinery. So `mod.rs` gains two: that a locus is flagged exactly when its tail
false-discovery rate meets the target and never on a ratio that is not a number, and that the
posterior is the logistic of the ratio shifted by the prior's log-odds, withheld where that shift
is undefined. **Both were mutation-tested**: dropping the finiteness screen from `flags`, and the
degenerate-prior screen from `posterior`, each fails its test.

**3. `cohort_inbreeding` is not copied, and that is a decision rather than an omission.** It is
the helper spec §3.2 says ng replaces with the parameters file's per-sample coefficient — it fills
a slice with one cohort-wide value repeated per sample. The span's end marker is its doc comment's
first line, and `the_item_the_span_stops_before_is_still_productions` asserts production still
declares it, so a rename cannot leave the guard green over a copy that silently grew.

## Tests

**83 under `ng::paralog::`**, up from 59: 58 are production's transcribed (23 `coverage_model`,
18 `locus_score`, **13 `prior`**, 4 `model_params`) and 25 are ng's own — 8 in `copy_fidelity`,
4 in `mod`, 13 in `production_parity`.

The twenty-four added here are production's thirteen for the prior, the curve and the threshold;
and eleven of ng's, of which **eight came out of this step's review** — see the
[fixes applied](../reviews/fixes_applied_ng_paralog_filter_a3_2026-09-06.md). The review found one
**Blocker**: two of the guard's three sanction checks tested only that the declared lines
*contained* certain substrings, so a declaration could carry an arbitrary edit into a frozen copy
with the suite green. Every kind is now exact and every one has a test in its rejecting direction.

**Both new copies are guarded, and both were mutation-tested.** Changing `DEFAULT_EM_TOL` from
`1e-9` to `1e-8` in ng's `prior.rs` fails at `prior.rs, line 6 of the copied content`; loosening
`flags`'s comparison from `<=` to `<` inside the span fails at `calibration.rs, line 63 of the
copied content — ng's copy and src/var_calling/paralog_filter/calibrate.rs have diverged`. Each
mutation was reverted from a backup and the restoration confirmed by `diff`.

## Validation

Run in the dev container, on this tree.

    cargo test --all-features --lib "ng::paralog::"
    → test result: ok. 83 passed; 0 failed; 0 ignored; 6290 filtered out

    cargo test --all-features --lib --bins --tests
    → 6358 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 3 errors, all needless_lifetimes in cohort_merge/build.rs and serial.rs

    cargo fmt --check
    → dirty on 9 files, none of them src/ng/paralog/

The integration failure, the three lints and the nine files are `main`'s at `a33ada0f`, proved so
before A1 began by removing `pub mod paralog;` and re-running. The gate — no failure `main` does
not already have — holds, and the lib count moves 6,334 → 6,352.

**Two clippy errors of ng's own appeared and were fixed rather than silenced**: the four dead-code
errors that led to the visibility widening above, and `enum_variant_names` on `WhyRepointed`,
whose variants had all begun with `A`.

## What this closes

**Milestone A is done.** Production's filter statistics are in ng, and the port's claim is
asserted two independent ways: textually, that the five copied files are production's byte for
byte but for ten declared lines — three in the score, seven in the calibration; and numerically, that the score, the prior, the curve and the
cut return the same values as production's on randomised inputs, compared by bit pattern.

## What the review changed, beyond the Blocker

- **The differential now reaches the verdict itself** — `flags` and `posterior` compared with
  production's, which it had stopped short of.
- **A bin-index error at either end of the shipped histogram is invisible through π and the
  curve**, and the reason is worth recording: the end bins sit at `±99.95`, where the logistic has
  saturated to exactly `0` and `1` in `f64`. Adding ratios past `±100` — 1,251 and 1,250 of 5,000,
  measured — did not help. A histogram over `[-6, 6]` does, and both clamp mutations die there.
- **`flags`'s target is swept rather than fixed**, which kills an implementation that reads a
  literal `0.01` in place of the operator's knob.
- **`lr_threshold` is what the run writes into the VCF header, and it is not what `flags` decides
  on.** The cut is the crossing bin's *centre*; the flag is the bin. Every ratio in the lower half
  of that bin — `0.05` wide at the shipped resolution — is flagged while sitting below the
  recorded cut. `calibration.rs` documents the two as equivalent; they are equivalent up to the
  bin, and that is now what is pinned. **Step C4 writes this number into the header.**

## Follow-ups

- **The differential still runs at one model-parameter setting** and at no cohort size near spec
  §4's three thousand. Both are size questions rather than correctness ones and belong to the plan's
  D milestone.
- **`CalibrationConfig` is copied but nothing constructs it yet.** Pass two does, at step C3.
- **Five equivalent mutants point at dead code inside production's frozen scorer**, carried from
  A2's review — production's backlog, not this plan's.
