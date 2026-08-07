# ng step 4, D3 — `fit_by_profile_scan` and the rail flag

**Date:** 2026-08-06. **Branch:** `ng-parameter-estimation`. **Plan:**
`doc/devel/ng/impl_plan/parameter_prepass_generic.md`, Milestone D step D3. **Design:**
`doc/devel/ng/arch/parameter_prepass_generic.md` §4.2,
`doc/devel/ng/spec/parameter_prepass.md` §3.1.

## What landed

**`src/ng/parameter_estimation/fitting/profile_scan.rs`** — a new file:

```rust
pub fn fit_by_profile_scan<M>(
    model: &M,
    cells: &[M::Cell],
    ladder: &[M::NoiseParams],
) -> ScanResult<M::NoiseParams>
where M: NoiseModel, M::Cell: WeightedCell, M::NoiseParams: Clone;
```

Step through the ladder; at each rung ask the model how likely each genotype makes each cell, climb
to the genotype frequencies that best explain them (D1), and keep the rung that scores highest. A
profile likelihood over the noise parameters.

**`src/ng/parameter_estimation/fitting/mod.rs`** — the `WeightedCell` trait, and `ScanResult`'s
frequencies re-keyed by ploidy.

**`src/ng/parameter_estimation/generic/noise_model.rs`** — `impl WeightedCell for Cell`.

## The four contract points, and how each is held

**Every rung is scored; no early exit.** Nobody has shown the profile curve has a single hump, and
that missing proof is the *reason* for a scan rather than a one-dimensional optimiser
(`spec/parameter_prepass.md` §3.1, §9.3). Two tests hold it: a fixture model records which rungs it
was asked about, and a deliberately two-humped curve — rung 0 a local summit, rung 3 the global one
— which the test first proves is two-humped by scanning the first two rungs alone and watching rung
0 win. Adding an early exit on the first decline fails the second.

**Ploidy travels with each cell, and the frequencies are climbed once per ploidy.** One noise
parameter is fitted across every ploidy its reads covered — chemistry does not know about
chromosomes — while a haploid region has two genotype classes and a diploid three, so they cannot
share a weight vector. A rung's score is the sum over ploidies. The test mixes haploid and diploid
cells under a model whose genotype count really does depend on ploidy, and asserts two frequency
sets of lengths 2 and 3, each summing to one.

**Ties resolve to the lower error rate**, which the plan asks be stated so two implementations
cannot differ. The scan cannot see which way its ladder runs, so what it can state is positional:
**a tie keeps the last of the tied rungs**. The generic ladder ascends in Phred and therefore
*descends* in error rate, so the last of a tied run is the lowest rate. Those are one fact split
across two files, so one test asserts both halves together — that `error_rate_ladder()` descends,
and that a wholly tied ladder returns its last rung. Reversing either without the other fails it.

**The rail flag.** The winning rung being first or last. Both ends are tested, on a ladder of five
rungs marching towards a truth that is on none of them, run forwards and reversed; the contrast that
gives it teeth is the recovery test, where the winner is interior and the flag is **false**. A
one-rung ladder is all edge and says so, which is the failure `ERROR_RATE_LADDER_RUNGS` guards
against from the other side.

## Deviations from the architecture, recorded

1. **`cells: &[M::Cell]` with an `M::Cell: WeightedCell` bound, not `&[(M::Cell, Ploidy, u64)]`.**
   The architecture sketches the tuple. Both paths' cell types already carry their ploidy and their
   site count, so the tuple would restate beside a cell what the cell knows — one more place for the
   two to disagree — and two of its three members are integers, which is the transposition this plan
   has hit four times. `WeightedCell` is the same move B2 made when it replaced the identical tuple
   with a struct.
2. **`ScanResult::frequencies` is a `BTreeMap<Ploidy, SmallVec<[f64; 3]>>`, not one vector.** The
   architecture's §4.2 says the frequencies are climbed once per ploidy and its §5.2 sketch gives
   `ScanResult` a single vector; the two cannot both hold. A single vector would mean picking one
   ploidy's answer to report and dropping the rest silently. Nothing consumed the field yet, so the
   change is free now and would not have been after E1.
3. **The width of a cell's row comes from `NoiseModel::genotypes`**, added in D2's review, and the
   scan asserts the model appended what it declared, per cell.

## Mutation record

Ten mutations, ten killed. The scan's own arithmetic is short; what it mostly does is decide, so
the mutations are decisions:

| mutation | killed by |
|---|---|
| exit the ladder at the first rung that scores worse | the two-humped curve (1) |
| a tie keeps the **first** of the tied rungs | the tie test (1) |
| the rail flag watches only the ladder's last rung | the reversed arm of the rail test (1) |
| the rail flag watches only its first | 2 |
| the rail flag is always `false` | 2 |
| the score is the last ploidy's, not the sum | the two-ploidy sum test (1) |
| the declared-versus-appended width check deleted | its refusal test (1) |
| the scratch buffer is not cleared between ploidies | 7 |
| every cell weighs one site rather than its own count | 4 |
| every cell is filed under the first cell's ploidy | 1 |

**One test was rewritten during this step because it could not fail as first written.** The rail
test originally used a two-rung ladder, where *both* rungs are edges and the flag is true whatever
the scan does. It now marches a five-rung ladder towards a truth on none of them, in both
directions. The fixture that replaced it was also wrong on its first attempt, and instructively:
blending each rung's columns towards the truth left them unnormalised, so the flattest rung won on
total column mass rather than on fit, and the closest rung lost. Each column is a genotype's
distribution over the cells and has to sum to one for the comparison to mean anything.

## Validation

All in the container. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → the library binary at **3,100
passed**, 0 failed, 5 ignored, up from 3,089; `cargo doc --no-deps --lib` at the 12-unresolved-link
pre-existing baseline. `ng::parameter_estimation` 187 → **198** tests, of which `profile_scan` is 11.

## Open, carried forward

- **Nothing has read a locus.** The scan's input is a slice of cells a caller materialises; E1 is
  the first step that builds one from an accumulator.
- **`fit_mixture_weights` can now become `pub(crate)`** — carried from D1, and D3 is the consumer
  that makes it possible. Left for the fix commit after this step's review, so that the change lands
  with a reviewer's eyes on it rather than inside the step that created the consumer.
- **The prefactor cost D2's review measured is still on the table.** Three `lgamma` per cell per
  rung, ~98% of the inner loop's arithmetic, identical at every rung. The scan is where a per-cell
  cache would live, and it does not have one.
