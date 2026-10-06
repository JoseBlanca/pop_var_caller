# Fit precision, step A2 — the information, summed in blocks by a pass

**Date:** 2026-09-27. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A2.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §3.2 (the block approximation), §3.6
item 5 (determinism). **Branch:** `fit-precision`.

## 1. What was built

A pass over the data can now be asked to sum the **information** the standard errors are computed
from: each position's scores (step A1) multiplied pairwise and added up over the positions. Only the
products the standard errors read are kept:

| block | size | what it pairs |
|---|---|---|
| the cohort's | 8 × 8, once | the eight cohort-level parameters with each other |
| a sample's own | 3 × 3, per sample | its first read group's two error rates and its homozygote excess, with each other |
| a sample's with the cohort's | 3 × 8, per sample | its three with the cohort's eight |

A product of two different samples' scores is never formed. That is the spec's block approximation:
two samples' parameters are related only through the cohort's, and that indirect relation is what
the error formula of step A3 accounts for. The blocks cost 264 bytes a sample (33 numbers), against
about 450 MB for the full matrix over kimura's 7,478 parameters (spec §3.1).

The sums follow the pass's existing chunking: each chunk sums its own positions, and the chunks are
joined in the same fixed order as every other count of the pass. So the result does not depend on the
number of threads.

**Nothing asks for them yet.** No pass of the fit asks (`PassKeeps::information` is false everywhere); step
A5 turns the final pass on, and step B1 the passes that refresh the errors while the fit runs. No
fitted number moves.

## 2. Assumptions and deviations

1. **The plan names a `SampleInformationBlock` of (k + 8)² per sample; what is built keeps the cohort's
   8 × 8 once, and per sample only the 3 × 3 and the 3 × 8.** A (k + 8)² block per sample would repeat
   the cohort's 8 × 8 in every sample, while the cohort's own block sums every position once. It is
   the same information; the type is `InformationSums`. k is 3 because of A1's finding: the likelihood
   reads only a sample's first read group's rates (A1 report §2 item 1).
2. **Three items deferred from the A1 review are done here, where the scorer gained its caller:**
   - **the scorer's inputs cannot disagree with the pass's**. `PassModel` bundles
     what a pass holds fixed — the parameters, the read-group map, the ploidy, the read-log tables and
     the two quadrature rules — and both `one_position` and `score_position` take it. `one_position`
     unpacks it into the same local names at its top, so its body and arithmetic are unchanged.
     `RuleSlopes` records the shapes it was built at, and `score_position` checks, in debug builds,
     that the slopes it is handed are those of the pass's own rules;
   - **each read's slope is computed once per position**: a table
     over (class, candidate, sample, genotype), filled before the branches, instead of recomputed
     inside the loop over quadrature nodes;
   - **the scorer is split by branch**: the shares, the two fixed-genotype branches,
     the segregating branch and the duplicated branch are each their own function.
3. **Names for the parameter slots move to A5**, where they are first printed; A2
   indexes the blocks by the existing slot constants.

## 3. Changes

- [fit.rs](../../../../src/parameter_estimation/joint/fit.rs): `PassModel` (built once a pass by
  `PassModel::new`); `one_position` takes it; `Statistics` gains `information: Option<InformationSums>`,
  merged in `absorb` (which refuses two chunks that disagree on keeping it); `expectation_pass` takes a
  `PassKeeps` naming what it keeps beyond its sums — the per-position posteriors it already kept on
  the final pass, and now the information — builds the scorer's tables once a pass when asked, and
  scores and adds each position after `one_position`.
- [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs): `InformationSums`
  (`new`, `add_position`, `absorb`); `ScoringTables` (the two rules' slopes and the log of each branch's
  share, once a pass); `score_position` takes the bundle and is split into five functions; `PositionScores` holds the per-position table of read slopes. The module's
  `expect(dead_code)` is gone: the pass calls the scorer.

## 4. Tests

Measured with `cargo test --release --lib parameter_estimation::joint::fit -- --nocapture` in the
container (`tmp/fit_precision/module_a2c.log`): 36 passed. The debug-build test below runs in the
full suite.

| test | what it shows |
|---|---|
| `a_pass_sums_the_scores_multiplied_pairwise` | Every entry of the pass's blocks agrees with the products formed position by position, on both ways a pass joins its chunks (in position order, as the final pass does; in the fixed halving tree, as the iterating passes do). On the A1 fixture — four samples, 600 positions, so three chunks — to **1.43 × 10⁻¹⁵ and 1.37 × 10⁻¹⁵** of the bound its row and column put on it; on one sample alone, 600 positions, to 6.40 × 10⁻¹⁵ and 6.54 × 10⁻¹⁵ (asserted below 10⁻¹²). Asking for the information leaves the log-likelihood and three other counts bit-for-bit unchanged. |
| `the_information_is_the_same_bits_at_any_pool_width` | Three samples, 3,172 positions (thirteen chunks): the blocks are present and the same bits at one, four and eight threads, on both ways of joining. |
| `a_sample_without_reads_carries_no_information` | A sample whose evidence is emptied at every position has its two blocks exactly zero, so its errors will be absent, not computed from rounding. |
| `tables_built_for_another_rule_are_refused` (debug builds) | The scorer panics when handed tables built at other shapes than the pass's own. |

Every A1 slope test still passes with the same printed disagreements as at A1 (the largest, carrier
shape `a`, 1.70 × 10⁻⁷). The review compared A1's and A2's scorers bit for bit: 0 of 600 positions
differed in any bit, in 42 configurations (depths 2, 8, 40; the duplicated class on and off; coverage
odds on and off; 8 and 16 nodes; homozygote excess at exactly 0 and 1).

## 5. Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` (`tmp/fit_precision/suite_a2b.log`): counts
  in the commit message; the three pre-existing failures only. **Both cross-platform checksum tests
  pass unchanged**, which is the check that bundling the pass's inputs moved no fitted bit; the
  review also found whole-fit dumps byte-identical before and after on three drawn cohorts, one of
  them a single sample.

## 6. What a collecting pass costs, and notes for A3 and A5

Measured by the review, one thread, depth 3, median of five interleaved runs in the container:

| samples | positions | duplicated class | plain pass | collecting pass | ratio |
|---|---|---|---|---|---|
| 4 | 200,000 | on | 0.439 s | 0.620 s | 1.41 |
| 64 | 20,000 | on | 0.401 s | 0.690 s | 1.72 |
| 64 | 20,000 | off | 0.250 s | 0.436 s | 1.75 |

The ratio grows with the samples because the scorer repeats, per sample, the loop over classes,
candidates and nodes that `one_position` runs. Spec §3.2 estimated "of the order of one pass"; at
kimura's 2,169 samples a collecting pass should cost about 1.75 plain passes, unmeasured there.

- **Memory on the final pass:** that pass holds every chunk's sums until it joins them in order —
  264 bytes × 2,169 samples × 128 chunks, about 73 MB at kimura, beside the per-position lists it
  already holds.
- **The two ways of joining agree to 1.66 × 10⁻¹⁵ relative but are not bit-equal**, so A5 must not
  expect the final pass's information to match an iterating pass's at the same parameters bit for
  bit.
- **At one sample the homozygote excess has non-zero information** (90.0 on a drawn 3,000-position
  sample), though the fit never moves it there; A3 reports it absent rather than inverting it.
- **With the duplicated class off the cohort block's last three rows and columns are zero**; A3
  restricts to the parameters that have information.
- **At a homozygote excess of exactly 1** the largest per-position score was 8.3 × 10⁵ and the largest
  block entry 6.9 × 10¹¹, all finite; nothing overflows at kimura's two million positions.
- Step A3 turns the blocks into standard errors.
