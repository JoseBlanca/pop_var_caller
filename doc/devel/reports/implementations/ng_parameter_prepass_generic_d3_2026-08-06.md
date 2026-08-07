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
where M: NoiseModel, M::NoiseParams: Clone;   // M::Cell: WeightedCell is on the trait
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
(`spec/parameter_prepass.md` §3.1, §9.3). **One test holds it against the natural mutation** — a
deliberately two-humped curve, rung 0 a local summit and rung 3 the global one, which the test first
proves is two-humped by scanning the first two rungs alone and watching rung 0 win. A second test
has a fixture model record which rungs it was asked about, but that one **survives** an exit taken
after the model call, which is where an exit would naturally go; it can only catch an exit taken
before it.

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
4. **`ScanResult` lives in `profile_scan.rs`, and the scan is its own file.** The architecture's
   module table put both the seam and the scan in `fitting/mod.rs`. `ScanResult` is constructed at
   exactly one place and its rail-flag doc repeated, almost word for word, the paragraph on
   `fit_by_profile_scan` — a contract stated in two files drifts. `fitting/mod.rs` now holds the
   seam and nothing else. **The architecture's module table is stale and is left for the owner**,
   as A4's `depth_bins.rs` row was.
5. **`WeightedCell` is a bound on `NoiseModel::Cell`, not on the scan.** Applied at review: with it
   on the function, a model whose cell knew neither its ploidy nor its site count compiled and
   failed only at the one call site that scanned it. On the associated type, "what must a path
   supply to be fitted?" is answered in one place — which matters, because the STR path will read
   this trait as a specification of what to build.

## Mutation record

**Twenty-eight mutations, twenty-eight killed** — ten before review, eighteen during it, every
survivor closed. The scan's own arithmetic is short; what it mostly does is decide, so the mutations
are decisions. Ten of the eighteen came from a review agent; three of those survived the committed
tests and are marked.

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
| **the zero-sites guard sums every cell rather than the ploidy's** — *survived* | now the two-ploidy refusal test (1) |
| **`genotypes > 0` deleted** — *survived* | now its refusal test (1) |
| **`REFERENCE_READ_TOLERANCE`-style: the model writes `NaN` or an all-`−∞` row at one rung** — *reached an unlocalisable message* | now the two rung-naming tests (1 each) |
| the impossible-cell check ignores the cell's weight | the empty-cell companion test (1) |
| `HashMap` in place of the `BTreeMap` the determinism claim rests on | the order-stability test (1) |
| the rail flag off by one (`== ladder.len()`) | 3 |
| the uniform start one entry too wide, or not a distribution | 13 each |
| the ladder iterated in reverse; the winner's frequencies taken from the last rung; frequencies filed under a constant ploidy; `cell_weights` built in reverse | 1–2 each |

**One test was rewritten during this step because it could not fail as first written.** The rail
test originally used a two-rung ladder, where *both* rungs are edges and the flag is true whatever
the scan does — a flag hard-coded to `true` passed it. It now marches a five-rung ladder towards a
truth on none of them, in both directions.

**The fixture that replaced it was wrong too, and the correction is worth more than the fault.** Its
`flat` baseline was one over the *genotypes* where a column is a distribution over the *cells*, so
the blend towards the truth did not preserve each column's sum. What that cost is **not** what this
report first said. Re-measured at review, with the five rung scores printed: the flattest rung —
the one carrying the most column mass — scored **worst of the five**, about 160,000 nats below the
winner. The extra mass moved the argmax **one rung inward**, from 4 to 3, and never came close to
overwhelming fit. The first telling of this said the flattest rung won, which would send the next
reader looking for a symptom that does not occur; the real one is an argmax one short of the end,
which is far easier to mistake for rounding. With `flat` at one over the cells the blend preserves
the sum by itself and the renormalising loop is a no-op — kept, so that it is a checked property of
the fixture rather than an assumption.

## Validation

All in the container. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → the library binary at **3,108
passed**, 0 failed, 5 ignored, up from 3,089; `cargo test --doc ng::parameter_estimation` → 1
passed; `cargo doc --no-deps --lib` at the 12-unresolved-link pre-existing baseline, none in
`parameter_estimation`. `ng::parameter_estimation` 187 → **206** tests, of which `profile_scan` is
20.

## Open, carried forward

- **Nothing has read a locus.** The scan's input is a slice of cells a caller materialises; E1 is
  the first step that builds one from an accumulator.
- **`fit_mixture_weights` cannot become `pub(crate)` yet, and D1's carried-forward item saying it
  can was wrong.** D3 did not become its consumer: the scan calls `climb_mixture_weights`, because it
  needs the score. Applied at review, `pub(crate)` produces `function fit_mixture_weights is never
  used` — a hard error under `-D warnings` — and breaks the `# Examples` doctest, the module's only
  public worked example, because a doctest is an external consumer. **The consumer is E2**, whose
  step 2 climbs the whole-sample frequencies once per ploidy. When E2 calls it, `pub(crate)` costs
  only the doctest, which should then move into the test module as a named test rather than be
  deleted.
- **The prefactor cost D2's review measured is still on the table.** Three `lgamma` per cell per
  rung, ~98% of the inner loop's arithmetic, identical at every rung. The scan is where a per-cell
  cache would live, and it does not have one.
