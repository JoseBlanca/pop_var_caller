# Review: the five performance commits of the penalised STR fit (2026-10-16)

**Scope.** Branch `str-penalised-fit`, five commits:

- `ae239e26`: each round takes the log chance of a read in each length bucket once (`ln_read_chances`), instead of once per distinct row of reads.
- `ef307745`: the research bench (`examples/ng_str_fit_bench.rs`) installs mimalloc, the allocator `pvc` uses.
- `e0dc3cf2`: the standard errors score a class's rows once for all the steps in F and in the allele-length shares. `RoundModel` is split into `ReadModel` (depends on the stutter), `GenotypePrior` (depends on the shares and the spread) and `ScoredDistinctRows` (a class's rows scored at one stutter).
- `8a554b08`: the table of distinct rows is built through one scratch row and an `AHashMap`, without a heap vector per (locus, sample) pair.
- `0084426a`: each chunk's partial weights hold only the distinct rows its own samples hold.

Each commit says every fitted number is bit-identical. The review's main question was whether any input could still make the bits differ that the 153-run reference fixture would not catch. The categories of `ai/skills/rust-code-review` were covered by one reviewer, without sub-agents, as the dispatch asked. Context: `doc/devel/reports/reviews/perf_estimate-parameters_2026-10-15.md` (findings H1, H2, L1), `doc/devel/ng/arch/str_penalised_fit.md`, `doc/devel/ng/spec/str_penalised_fit.md`.

**Verdict.** No correctness defect found. The bits claim holds by construction, and two probes confirm it on inputs the fixture does not contain: before and after the five commits, at 1, 3, 7 and 8 threads. There are four small findings: one comment states a behaviour that does not happen (Minor), one field name collides with another meaning of "rows" (Minor), and two are wording or wrong numbers in comments and commit messages (Nits).

## 1. What was checked, and how

### 1.1 Why the bits cannot move, change by change

| Change | Why every sum keeps its bits |
|---|---|
| Log table (`ae239e26`) | `ln_reads` (`expected_counts.rs:578`) adds `reads × ln(½(p_a + p_b))` for the non-empty buckets in bucket order, starting from `0.0`. That is the same products in the same order as `add_ln_reads` (`ssr_fit.rs:4155`). Each log is computed from the same argument, `max(½(p_a + p_b), 10⁻³⁰⁰)`. The only change is *when* each log is taken. |
| No prior in a scored row (`e0dc3cf2`) | `counts_of` now computes `(scaled × prior) × weight` where the old code stored `scaled × prior` and multiplied it by the weight later. The products are the same, with the same grouping. `GenotypePrior::totals` sums the same products with `Iterator::sum`, as `score` did. |
| Curvature columns reordered (`penalised_fit.rs:1563-1586`) | Each entry of the curvature receives exactly two additions to `0.0`. Entry (r, c) gets one from column c and one from column r; a diagonal entry gets both from its own column. The sum (0 + a) + b equals (0 + b) + a for every pair of doubles, `-0.0` included. Steps in columns 3 and above never move the stutter, because `numbers_of` reads the stutter from `theta[0..3]` alone. So the scored table they reuse is exactly the one a fresh scoring would give, and the `debug_assert_eq!` at line 1575 guards a property that holds by construction. Release builds do not depend on it. |
| Distinct-row table (`8a554b08`) | Rows are still numbered in the order first met: `index_of.len()` before an insertion equals `reads.len()` before a push. The map is iterated only to put each row at its own number, so the hasher's random seed never reaches the output. |
| Sparse partial weights (`0084426a`) | Each row's chunk totals are still added in chunk order. What is dropped is the `+0.0` an untouched chunk used to add. Adding `+0.0` changes only `-0.0`. A running total that starts at `+0.0` can never become `-0.0` under round-to-nearest, because a sum is `-0.0` only when both terms are. The chunk size is `max(rows / 256, 4096)`, which depends on the data alone, never on the thread count. |

### 1.2 Probe 1: the old code beside the new, on edge cases

`tmp/review_perf/src_copy/src/parameter_estimation/joint/ssr_fit/expected_counts/review_probe.rs` holds the code as it stood at `53a3e4f9`, before all five commits: `RoundModel`, the dense partial weights, the `HashMap` table builder and `ordered_row` over `add_ln_reads`. It runs that code beside `0084426a`'s and compares every number through `{:?}`, which writes a double so that it reads back to the same bits (`-0.0` and `NaN` included).

The inputs are drawn classes with these cases added:

- samples with no reads (empty `by_library`);
- a sample whose only library holds an all-zero row;
- a sample with two libraries;
- homozygote excesses of `0.0`, `-0.0`, `1.0`, `0.999999999`, `10⁻³⁰⁰` and `0.37`;
- shares with exact zeros, and with one allele at share 1;
- spreads of `0`, `10⁻⁹`, `0.7`, `10³⁰⁰` and `+∞` (the last makes the two-draw chance `NaN`);
- classes from 4 rows (3 distinct) up to 14,029 rows: four chunks, where some distinct rows appear in only some chunks.

The new code ran on pools of 1, 3 and 7 threads; the old code ran at 8.

| Probe test | Result |
|---|---|
| `probe_distinct_rows_old_and_new_agree` (old builder at 1, 7 and 64 tracts a block) | same rows, same numbers, same excess bits |
| `probe_expected_counts_old_and_new_agree_on_edge_cases`: whole class and tract by tract, 180 combinations | all bit-identical |
| `probe_scored_table_reuse_matches_fresh`: one scored table reused across shares and spreads, as the standard errors use it | bit-identical |
| `probe_ordered_row_old_and_new_agree`: 10,000 random rows at five stutter levels, counts up to 2²⁶ a bucket | bit-identical |

**The probe can fail.** Merging the chunks in reverse order (`partials.iter().rev()`) failed two of the four tests. So the probe sees a change in the order of summation.

### 1.3 Probe 2: the whole fit, before and after, at two thread counts

The research bench was built at `53a3e4f9` and at `0084426a`. It ran six cells that are not in the fixture's grid, each with `--loci-at-most 100 --methods D_tract_prior --fixture-out`. The old build ran at `RAYON_NUM_THREADS=8`, the new one at 1 and at 7. Every number was compared as written text, so `-0.0` differs from `0.0` (`tmp/review_perf/cmp_exact.py`). The fixture's own comparison script, `cmpjson.py`, compares through Python floats and treats `-0.0` and `0.0` as equal.

| Cell | Fields | Differ (new at 1 thread / new at 7) |
|---|---|---|
| one-base, 200 samples × 30 reads, PCR, seed 7 | 1,135 | 0 / 0 |
| two-base, 200 × 30, PCR, F = 0.4, seed 7 | 1,695 | 0 / 0 |
| one-base, 3 × 300, PCR-free, seed 9 | 1,128 | 0 / 0 |
| two-base, 1 × 30, PCR, seed 11 | 1,658 | 0 / 0 |
| two-base, 300 × 3, PCR-free, seed 5 | 1,696 | 0 / 0 |
| one-base, 120 × 10, PCR, F = 0.4, seed 3 | 1,132 | 0 / 0 |

These cells include standard errors with F read off pairs and with F in the counting, one sample, and 300 reads a position. Outputs are in `tmp/review_perf/bench_runs/`. Nothing ran at 2,000 samples.

### 1.4 Memory: one class's table at a time

`ScoredDistinctRows<'a>` borrows one `DistinctRows` and owns only that class's scored rows. So it cannot hold two classes' tables. In the standard errors, the central table (`central`, `penalised_fit.rs:1570`) is dropped at line 1580, before each stutter step builds its own. A scored row now holds one vector of 190 doubles where it held two, so the scored table is half its former size.

### 1.5 Dead code

`add_ln_reads` (`ssr_fit.rs:4155`) is still called in production by `TractLikelihoods::of`, the climb's scoring (lines 4194 and 4203). `pooled_reads` is `#[cfg(test)]` (`expected_counts.rs:1211`), and its only caller, `expected_counts_by_tract_row_by_row`, is test-only too. Nothing dead was left in production. `RoundModel`, `one_total` and `two_total` have no references left in `src/` or `doc/devel/ng/`.

### 1.6 Verification commands

Run on an exported copy of `0084426a` (`tmp/review_perf/src_copy`), natively, `CARGO_TARGET_DIR=tmp/review_perf/target_probe`:

- `rustfmt --check --edition 2024` on `expected_counts.rs`, `penalised_fit.rs` and `examples/ng_str_fit_bench.rs`: clean.
- `cargo clippy --release --all-targets --all-features -- -D warnings` (clippy of rustc 1.98.0): 2 errors in the library (`fit_bench.rs:460` useless conversion, `run/edge_deletion_gathering.rs:201` single-pattern match) and 3 more in its tests (`expected_counts.rs:2123` disallowed `f64::powi`, `penalised_fit.rs:2915-2916` needless range loops). `git blame` puts all five in commits older than `53a3e4f9`, outside this scope. They may come from a newer clippy than CI uses; that was not checked.
- `cargo test --release --features bench-fixtures --lib ssr_fit`: `162 passed; 0 failed; 5 ignored`, which matches the commits' "162 passed".

## 2. Findings

### Minor

- `examples/ng_str_fit_bench.rs:48-53` — **[Minor]** The comment says "dhat takes the slot when both features are on", but this bench declares no dhat allocator
- **Confidence:** High
- **Problem:** The mimalloc declaration is gated on `not(feature = "dhat-heap")`, and the comment says dhat then takes the slot. No `dhat::Alloc` is declared anywhere in the bench (`grep global_allocator`: `src/main.rs`, `ng_generic_walk_probe`, `ng_call_from_psps_cost`, `dhat_ng_joint_fit`, `dhat_ng_merge`, not this bench). With `--all-features`, which CI uses, or with `--features dhat-heap`, the bench therefore runs the **system** allocator. That is the situation commit `ef307745` set out to end. `src/main.rs:19` gates mimalloc on `alloc-mimalloc` alone.
- **Why it matters:** A timing or peak-memory figure taken from a build with both features names mimalloc in the comment but measures the system allocator. `Cargo.toml:216-224` names exactly this as the case to refuse at startup.
- **Suggested fix:** Either gate on `alloc-mimalloc` alone, as `main.rs` does (the bench has no dhat allocator to collide with), or declare `dhat::Alloc` under `dhat-heap` and refuse to start with both features, as `ng_generic_walk_probe` does. The first is one line:
  ```rust
  // The allocator `pvc` ships (`src/main.rs`), so that the bench's times and peak memory describe
  // what `estimate-parameters` runs. A `#[global_allocator]` is per binary: without this line the
  // bench ran the system allocator.
  #[cfg(feature = "alloc-mimalloc")]
  #[global_allocator]
  static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;
  ```

- `src/parameter_estimation/joint/ssr_fit/expected_counts.rs:832-836, 876-897` — **[Minor]** `ScoredDistinctRows::rows` and `DistinctRows::rows` are different things, and `expected_counts` uses both on adjacent lines
- **Confidence:** High
- **Problem:** `DistinctRows::rows` is one `(distinct-row number, excess)` entry per (locus, sample) pair. `ScoredDistinctRows::rows` is one `ScoredRow` per *distinct* row. In `expected_counts`, `distinct.rows.par_chunks(chunk)` walks the first, and `self.rows[index].largest` indexes the second, inside the same closure. The loop variable `rows` (line 880) is a chunk of the first kind.
- **Why it matters:** Mixing the two is the kind of mistake that indexes the wrong table. Today it would panic or give wrong counts only under a new edit, but the reader has to keep track of which "rows" each line means.
- **Suggested fix:** Rename the field to what it holds, for example `scored` (the name the old local had), and update its doc line: "`scored`: each distinct row's [`ScoredRow`], in `distinct.reads` order".

### Nits

- `expected_counts.rs:869-870`: "with up to 256 chunks". The chunk is `max(rows / 256, 4096)` rows, so a class of at least 256 × 4,096 rows whose row count is not a multiple of 256 gets 257 chunks (for example 1,052,671 rows: chunk 4,111, 257 chunks). Say "about 256", or use `rows.div_ceil(256)` if exactly 256 matters. The same wording is in the 2026-10-15 review's L1.
- `expected_counts.rs:872` is 169 characters long. Text meant for the next line was joined onto this one, so the paragraph breaks mid-thought. Lines 555 and 559 (`ln_read_chances`' doc, from `ae239e26`) also run past 100 characters. Rewrap.
- `expected_counts.rs:746-748`, `GenotypePrior`'s doc: "Each genotype's prior … given a sample's homozygote excess". The struct holds no excess; the excess comes in through `add_sample`. Suggest: "…split into the part where the sample's two copies are one draw and the part where they are two; [`Self::add_sample`] mixes the two by the sample's homozygote excess and the spread."
- `expected_counts.rs:578-579`: `ln_reads` zips the table with the row, and only a `debug_assert_eq!` checks that their lengths match. `add_ln_reads` panicked on a row longer than the table; in a release build `ln_reads` drops the extra buckets. The lengths match by construction (both are `2 × read_span + 1` from one `evidence`). The only way to break it is to pass `ScoredDistinctRows::of` a `distinct` built from another evidence, which its doc forbids and no caller does. A single `assert_eq!(evidence.read_span …)` in `ScoredDistinctRows::of`, checked once per table, would turn that misuse into a panic again, at no cost per row.
- Commit message `e0dc3cf2`: "36 of a class's 42 steps" holds when F is read off pairs: 3 stutter numbers plus 18 shares make 21 columns and 42 steps, and the 18 shares account for 36 of them. When F is in the counting, there are 22 columns and 44 steps, and 38 of them reuse the table. The message says "steps in F or in an allele-length share", which counts F, so the figure it gives should be 38 of 44. (Commit messages are not changed here; noted for the implementation report.)
- `penalised_fit.rs:1569-1580`: when a class has no shares and no F column, so `width == 3`, the central table is scored and dropped unused. That happens only at `allele_span = 0` (one allele class); production uses 9. Building `central` only when `width > 3` would avoid it.

## 3. Quantitative claims checked

| Claim | Where | Verdict |
|---|---|---|
| Bit-identical results, thread count included | all five commits | **Checked, correct**: probes 1 and 2 above |
| "each entry receives exactly two additions to zero" | `penalised_fit.rs:1566`, `e0dc3cf2` | **Checked, correct** (§1.1) |
| "a sum begun at `+0.0` is never `-0.0`" | `expected_counts.rs:873-875` | **Checked, correct** under round-to-nearest, the only rounding mode Rust uses |
| "up to 256 chunks" | `expected_counts.rs:870`, `0084426a` | **Wrong**: up to 257 (Nits) |
| "36 of a class's 42 steps" | `e0dc3cf2` | **Correct only with F off pairs**: 38 of 44 with F in the counting (Nits) |
| 6.3 GB → 2.9 GB peak at 1,000 loci × 2,000 samples × 30 reads | `expected_counts.rs:871-872`, `0084426a` | **Not re-measured**: 2,000-sample runs were out of this review's bounds |
| "`add_ln_reads` stays for the row-by-row test oracle and the climb" | `ae239e26` | **Checked, correct** (§1.5) |

## 4. Scratch

Everything is under `tmp/review_perf/`:

- `src_copy/`: `0084426a`, with the probe module `src/parameter_estimation/joint/ssr_fit/expected_counts/review_probe.rs`. It is unhooked now. To run it again, append `#[cfg(test)] mod review_probe;` to `expected_counts.rs` and run `cargo test --release --lib review_probe`.
- `before_src/`: `53a3e4f9`.
- `run_bench.sh`, `bench_runs/`, `bench_runs.txt`, `cmp_exact.py`: probe 2.
- `target_before/`, `target_after/`, `target_probe/`: build directories, about 5 GB in all, which can be deleted.
