# Code Review: str_penalised_fit_D6
**Date:** 2026-10-16
**Reviewer:** rust-code-review skill (orchestrator; categories covered by one reviewer, no sub-agents, at the caller's request)
**Scope:** the uncommitted change on `str-penalised-fit` at HEAD 2c0bdf99 — plan step D6, the proxy for calling (`--calls-out`)
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the working-tree diff (`git diff`) of four files, 439 lines added, 16 removed.
- **Reviewed against:** HEAD 2c0bdf99 (a docs-only commit on top of 830049c1, the code the implementer's report names), plan step D6 of `doc/devel/ng/impl_plan/str_penalised_fit.md`, spec `doc/devel/ng/spec/str_penalised_fit.md` §5.1, and the implementer's report `tmp/d6/impl_report.md`.
- **In-scope files:**
  - [examples/ng_str_fit_bench.rs](../../../../examples/ng_str_fit_bench.rs) — `--calls-out`, `call_lines`, `call_line`, `write_call_lines`
  - [src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs) — `ClassCaller`, `drawn_genotypes_in_window`, `CallChanges`, `call_changes`, one test
  - [src/parameter_estimation/joint/ssr_fit/bench_fixtures/cohort_draw.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/cohort_draw.rs) — `DrawnPeriod::genotypes`, `draw_locus_and_genotypes`
  - [src/parameter_estimation/joint/ssr_fit/expected_counts.rs](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs) — `ScoredDistinctRows::genotypes_weighted_highest`
- **Read as direct callees, not reviewed as changes:** `DistinctRows::of`/`in_blocks_of`, `GenotypePrior::{of, totals, add_sample}`, `ReadModel::{of, score, counts_of}`, `ClassNumbers::of_truth`, `ClassTruth::shares_in_window`, `numbers_for_scoring`, `variants_of`, `library_level`, `genotype_pairs`.
- **Deliberately out of scope:** the 2,000-sample runs (not done by the implementer; the caller forbade anything above 500 samples here).
- **Categories covered:** reliability (correctness of pairing, the call, clipping, the arms; challenge tests), refactor_safety (draw unchanged), float_portability (thread-count independence of the call), errors and defaults (CLI refusals), naming/idiomatic (light). `unsafe_concurrency`, `module_structure`, `tooling` had nothing to act on: no `unsafe`, no new locks, no new modules or dependencies.

## 2. Verdict

**Approve-with-changes.** The four things the caller asked about all hold, each checked by running code, not by reading the report:

1. **The comparison is the plan's.** Calls and drawn genotypes are paired by the same (locus, sample), in the same order; the call is the argmax of exactly the expectation step's genotype weights; the ±9 window clips truth and calls alike; the three stutter arms hold the shares and F at the truth's.
2. **The draw is unchanged.** Two cells, run with a HEAD build and a build of the change, give identical rows (timings aside), truth, locus and reference-fixture files — one of them with two libraries at multipliers 1 and 1.6 and the allele reach moved to 12.
3. **Production is untouched.** Every changed line a release `pvc` build could compile sits under `#[cfg(any(test, feature = "bench-fixtures"))]`.
4. **The 500 × 3 numbers reproduce.** Two seed-1 cells (one-base and two-base) re-run from the exported tree give calls files byte-identical to the implementer's, and the two-base one is byte-identical again at one thread.

One change is wanted before merge: the test that claims to check the pairing does not (M1). It would pass with every sample's genotype shifted by one within each locus.

## 3. Execution status

All builds from exported copies under `tmp/review_d6/` (`src_copy` = the working tree, `head_copy` = `git archive HEAD`), native `cargo`, `-j 4`, `RAYON_NUM_THREADS=4` at most.

| command | exit | result |
|---|---|---|
| `cargo build --release --features bench-fixtures --example ng_str_fit_bench` (both trees) | 0, 0 | built in 2 m 36 s and 3 m 03 s |
| `cargo fmt --check` (src_copy) | 0 | clean |
| `cargo clippy -j 4 --all-targets --all-features` (src_copy) | 0 | warnings only at `expected_counts.rs:2179`, `penalised_fit.rs:2915`, `:2916`, `fit_bench.rs:465`, `edge_deletion_gathering.rs:201` (the five known, none in the diff) and in `vendor/noodles-cram` |
| `cargo check --bin pvc` (src_copy, no feature) | 0 | 27 warnings, every one in `vendor/noodles-cram`; none in `src/` |
| `cargo test --release --lib --features bench-fixtures -- fit_bench:: cohort_draw:: expected_counts::` (probe copy) | 0 | `30 passed; 0 failed` — the implementer's 29 plus the reviewer's probe |
| draw-unchanged runs, cells A and B, HEAD vs change | 0 | see §6 "Checks that passed" |
| 500 × 3 reproduction, `p2` and `p1` seed 1 | 0 | calls files byte-identical to `tmp/d6/n500/` |

Not run: `cargo doc`, `cargo audit` (no dependency change); `--all-features -D warnings` (the five known clippy findings would fail it, as on HEAD).

Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

1. **Assumed:** `draw_period` always draws at a Poisson depth, so every sample listed at a locus has at least one read and "pairs with reads" equals "samples listed". `Depth::Fixed` lists samples with no reads; it is used by `draw_stratum`, not by `draw_period`, so it does not reach `--calls-out` today. Affects no finding.

## 5. Top 3 priorities

1. **M1** — make the test's pairing check able to fail: a one-sample shift within each locus still passes it (1,743 of 2,293 calls match the shifted genotypes).
2. **Mi1** — say in the `--calls-out` documentation that an `all` row sums only the classes that arm was computed for, so two arms' `all` rows can have different `pairs`.
3. **Nits** — a CLI misuse that panics instead of a clap error, and two small prose corrections in the implementer's report.

## 6. Findings

### Checks that passed (the caller's four questions)

**Pairing.** `DistinctRows::in_blocks_of` ([expected_counts.rs](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L475-L537)) walks `evidence.tracts` in order and each tract's `samples` in order; the blocks are `par_chunks` collected in order and re-flattened with `flat_map_iter`, so `rows` is (tract, sample) order. `genotypes_weighted_highest` maps over `distinct.rows` with an order-preserving `par_iter().collect()`. On the draw side, `draw_locus_and_genotypes` pushes the genotype immediately after pushing the sample's reads, both after the `depth == 0` skip, so the two lists have the same members in the same order ([cohort_draw.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/cohort_draw.rs#L695-L730)); `draw_period` collects the loci in order with `unzip`, and `drawn_genotypes_in_window` flattens loci then samples. `drawn` is bound immutably in `main` and every method borrows `&drawn.strata`, so nothing reorders or filters the tracts between draw and call. Measured: with the true numbers on the test's two classes, 2,151 of 2,293 calls are the drawn genotype; shifting each locus's genotypes by one sample drops that to 1,743, reversing the class to 876.

**The call is the fit's weights.** In the expectation step a sample's weight on genotype g is `scaled[g] · (one_draw[g] · c₁ + two_draws[g] · c₂) / T`, split as `one + two` in `ReadModel::counts_of` ([expected_counts.rs](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L683-L706)), with `as_one = c₁/T` and `as_two = c₂/T` from `GenotypePrior::add_sample`. The new code computes exactly `scaled · (one_draw · as_one + two_draws · as_two)` from one sample's `RowWeights` ([expected_counts.rs](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L1012-L1022)). `T` is one positive number per sample, so the argmax is the posterior mode. The returned genotype is `model.genotypes[g]`, built by `genotype_pairs` with the first class not above the second, which is the ordering `drawn_genotypes_in_window` puts the drawn pair in.

**Window clipping, identical on both sides.** The calls live on 19 classes, `-9 ..= +9`, with `allele_span = len(shares) / 2 = 9`. The truth's shares are folded into the end classes by `ClassTruth::shares_in_window` with `clamp(-allele_span, allele_span)` ([cohort_draw.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/cohort_draw.rs#L627-L635)); the drawn copies are folded by the same clamp in `drawn_genotypes_in_window` ([fit_bench.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs#L1873-L1888)). Every arm is called on the same `DistinctRows` and the same window.

**The arms.** `with_stutter` builds `ClassNumbers { slippage, ..truth.clone() }`, so `level_0.995`, `fitted_level` and `fitted_stutter` keep the truth's `allele_length_variation` and `concentration`; `fitted_level` also keeps the truth's `shorter_share` and `fall_off` ([ng_str_fit_bench.rs](../../../../examples/ng_str_fit_bench.rs#L1199-L1218)). `Slippage` has exactly three fields, so `fitted_stutter` changes all three and nothing else. The reference reads were drawn through the truth's stutter only at multiplier one: `library_level` returns the class level unchanged at 1.0, and `main` refuses any other multiplier with `--calls-out`.

**Draw unchanged.** Cell A: `--period 1 --samples 200 --depth 3 --seed 7 --loci-at-most 20 --heldout-loci 5 --methods B,D_tract_prior` with `--truth-out` and `--locus-out`. Truth and locus files `cmp`-identical; rows file identical once `*_seconds` lines are removed (1,246 lines). Cell B: `--period 2 --samples 100 --seed 3 --loci-at-most 30 --library-multipliers 1,1.6 --allele-reach 12 --methods D_tract_prior --fixture-out` — the JSON fixture is equal key for key once `seconds` is removed. A fit is a deterministic function of the reads, so identical fitted numbers and held-out scores mean identical draws. Files: `tmp/review_d6/draw_check/`.

**500 × 3 reproduced.** `p2_n500_d3_pcr_s1` and `p1_n500_d3_pcr_s1` re-run with the change's binary: both calls files byte-identical to `tmp/d6/n500/`, and the `p2` rows identical apart from timings and peak memory. `p2` re-run at `RAYON_NUM_THREADS=1`: calls file byte-identical to the four-thread one. Files: `tmp/review_d6/repro/`.

**The report's figures.** Every 500 × 3 number in `tmp/d6/impl_report.md` was recomputed from the calls files (`tmp/review_d6/check_numbers.py`): the table's sums; 52 and 55 changed calls; 3.76 and 4.40 changed in 100,000; 6.27% and 4.75% wrong under the truth; 1,668 and 1,080 for "1,100 to 1,700 times smaller"; 2.16 and 1.23 in 10,000 for `fitted_stutter`; 112.9 and 161.0 in 10,000 for `fitted_all`; 42 and 63 class fits; fitted slip rate median −0.67% (−3.42% to +5.32%) and −0.63% (−5.04% to +4.10%); 43 of 52 changes in the 18-repeat class (24, 10, 9 by seed); 1,897 wrong in 10,468 pairs there at seed 1; fits 4.0 to 10.2 s, comparisons 0.4 to 0.8 s. All correct. Two prose slips are under Nits.

### Blocker

None.

### Major

- `src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs:2701` — **[Major]** M1: the test's pairing check passes with the pairing broken
- **Confidence:** High
- **Problem:** The test's documentation says "most calls under the truth are the drawn genotype, which checks that calls and drawn genotypes are paired in the same order" ([fit_bench.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs#L2658-L2662)), and the assertion is `2 * same.wrong_under_reference < same.pairs` ([fit_bench.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs#L2701-L2704)). Most samples at these two classes are homozygous for the reference (1,413 of 2,293), so a mispairing still matches most calls. Measured on the test's own draw (probe test in `tmp/review_d6/probe_copy`): correct pairing 2,151 of 2,293 calls right (142 wrong, 6.2%); every locus's genotypes rotated by one sample, 1,743 right (550 wrong, 24%) — **passes** the assertion; the whole class reversed, 876 right — fails it.
- **Why it matters:** The pairing is the one thing that, if it broke (a reordered `DistinctRows`, a filter on tracts, a change to which samples a draw lists), would turn every D6 number into noise without a panic. The code is correct today; the guard against its regression is not.
- **Suggested fix:** assert against the broken pairing directly, so the threshold needs no tuning:
  ```rust
  // Calls and genotypes paired as drawn match far more often than with each locus's
  // genotypes shifted by one sample (measured: 142 against 550 wrong of 2,293).
  let shifted: Vec<Vec<(i8, i8)>> = genotypes_by_locus
      .iter()
      .map(|locus| { let mut l = locus.clone(); if !l.is_empty() { l.rotate_left(1); } l })
      .collect();
  // ... accumulate `wrong_shifted` with call_changes(&under_truth, &under_truth,
  //     &drawn_genotypes_in_window(&shifted, ALLELE_SPAN)).wrong_under_reference
  assert!(2 * same.wrong_under_reference < wrong_shifted, "{same:?} against {wrong_shifted}");
  ```
  and correct the documentation sentence to say what the assertion checks.

### Minor

- `examples/ng_str_fit_bench.rs:1194` — **[Minor]** Mi1: an arm's `all` row can cover fewer classes than another arm's, undocumented
- **Confidence:** High
- **Problem:** `over_all[arm]` adds a class only when that arm exists for it ([ng_str_fit_bench.rs](../../../../examples/ng_str_fit_bench.rs#L1225-L1239)). `fitted_level` and `fitted_stutter` are absent for a class the method gave no stutter, and `fitted_all` for a class `numbers_for_scoring` gives no numbers, while `level_0.995` is always present. The module documentation says such a class "has no row for the arms that need it" ([ng_str_fit_bench.rs](../../../../examples/ng_str_fit_bench.rs#L108-L109)) but not that the `all` rows then sum over different classes. In the six 500 × 3 cells every arm's `pairs` agree, so no reported number is affected.
- **Why it matters:** a reader comparing `changed` between arms on the `all` rows (as the report's table does) would compare different denominators without being told, at small cohorts or with method `B`, where classes go unfitted.
- **Suggested fix:** one sentence in the module documentation's "The proxy for calling": "An `all` row sums the classes its arm was computed for; compare `pairs` before comparing two arms' `changed`." Optionally `eprintln!` when an arm's `pairs` differ from `level_0.995`'s.

### Nits

- [examples/ng_str_fit_bench.rs:1267-1270](../../../../examples/ng_str_fit_bench.rs#L1267-L1270): `--calls-out` with `--library-multipliers` other than one is refused by `assert!`, a panic with a backtrace hint, where the bench's other refusals are clap errors. A `clap::Error` (or `eprintln!` + `exit(2)`) would read as a usage error.
- [cohort_draw.rs:827](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/cohort_draw.rs#L827): the genotype offsets are kept as `i8` in every drawn run, `--calls-out` or not, so `--allele-reach 128` now panics where it used to draw. No plausible run uses a reach above 127 (the default is 14, the read span 8); a sentence on `--allele-reach` would document the bound. The kept genotypes also add 2 bytes a (locus, sample) pair to the draw-only memory runs, under a few percent of a pair's reads.
- [ng_str_fit_bench.rs:1179](../../../../examples/ng_str_fit_bench.rs#L1179): "one over every class a arm" → "an arm".
- `tmp/d6/impl_report.md`: (a) "The one-base classes of 17 to 19 repeats, whose true slip rate is the highest of the grid (0.27)" — 0.27 is the truth at every class from 16 to 21; the `level_0.995` changes sit in 18 (43), 17 (2), and one or two each in 9, 10, 11 and 15, none in 19. (b) "on HEAD 830049c1" — HEAD is now 2c0bdf99, a docs-only commit on top of it; the byte-identity claim holds against either (re-checked here against 2c0bdf99).

## 7. Out of scope observations

None.

## 8. Missing tests to add now

**`genotypes_weighted_highest`**
- `genotypes_weighted_highest_matches_the_argmax_of_a_brute_force_posterior` — input: a handful of drawn samples at one class. Computes each sample's posterior over the 190 genotypes directly (prior from `shares`, `concentration` and the sample's excess as in `expected_counts_by_tract`'s brute-force path around [expected_counts.rs:1105-1160](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L1105-L1160), times the read likelihood) and asserts the argmax agrees with `genotypes_weighted_highest`. Catches a weight that drifts from the expectation step's (a dropped `one_draw` term, `as_one`/`as_two` swapped), which the existing test would not see because it only compares the function with itself.

**`drawn_genotypes_in_window`**
- `drawn_genotypes_in_window_folds_far_copies_to_the_end_classes_and_orders_the_pair` — input `[[(-12, 3), (5, -1), (10, 10)]]` at span 9; expected `[(0, 12), (8, 14), (18, 18)]`. Catches a fold that drops instead of clamping, or an unordered pair (which would make every heterozygote drawn "larger first" a wrong call).

**Pairing (see M1)**
- `the_same_stutter_changes_no_call_and_a_wrong_one_changes_some` strengthened as M1 suggests.

## 9. What's good

- `draw_locus` is now a one-line wrapper over `draw_locus_and_genotypes`, so the random stream is provably one code path rather than two kept in step ([cohort_draw.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/cohort_draw.rs#L686-L692)).
- The call reuses `GenotypePrior::add_sample` and the scored rows instead of re-deriving the posterior, so "the fit's own weights" is true by construction ([expected_counts.rs](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L988-L1027)).
- The stutter arms are built by struct update from the truth (`..truth.clone()`), which makes "only the stutter differs" visible in one line ([ng_str_fit_bench.rs](../../../../examples/ng_str_fit_bench.rs#L1201-L1204)).
- `call_changes` asserts the three lists have one length before zipping, so a mispairing by count fails loudly ([fit_bench.rs](../../../../src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs#L1929-L1939)).

## 10. Commands to re-verify

```
# exported trees: tmp/review_d6/src_copy (working tree), tmp/review_d6/head_copy (git archive HEAD)
cd tmp/review_d6/src_copy && CARGO_TARGET_DIR=../target_new cargo build -j 4 --release --features bench-fixtures --example ng_str_fit_bench
cd tmp/review_d6/head_copy && CARGO_TARGET_DIR=../target_head cargo build -j 4 --release --features bench-fixtures --example ng_str_fit_bench

# draw unchanged (run with target_head and target_new; compare all but *_seconds / "seconds")
RAYON_NUM_THREADS=4 <bench> --period 1 --samples 200 --depth 3 --seed 7 --truth pcr --loci-at-most 20 \
  --heldout-loci 5 --methods B,D_tract_prior --out A.rows.tsv --truth-out A.truth.tsv --locus-out A.loci.tsv
RAYON_NUM_THREADS=4 <bench> --period 2 --samples 100 --depth 3 --seed 3 --truth pcr --loci-at-most 30 \
  --library-multipliers 1,1.6 --allele-reach 12 --methods D_tract_prior --fixture-out B.fixture.json

# 500 x 3 reproduction (compare with tmp/d6/n500/)
RAYON_NUM_THREADS=4 target_new/release/examples/ng_str_fit_bench --period 2 --samples 500 --depth 3 --seed 1 \
  --truth pcr --loci-at-most 100 --folds 0 --methods D_tract_prior --out r.rows.tsv --calls-out r.calls.tsv

# report figures
uv run --no-project python tmp/review_d6/check_numbers.py tmp/d6/n500

# pairing probe (M1)
cd tmp/review_d6/probe_copy && CARGO_TARGET_DIR=../target_probe cargo test -j 4 --release --lib \
  --features bench-fixtures -- review_probe --nocapture
```

`PROJECT_STATUS.md` was not updated: the caller set the worktree read-only apart from this report.
