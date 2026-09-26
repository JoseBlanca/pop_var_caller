# Issue 3: contamination's markers in two passes

*Implementation report, 2026-09-26. Plan:
[estimation_memory.md](../../implementation_plans/estimation_memory.md) §6. Branch
`contamination-two-pass`.*

## Plan

Contamination is measured at *markers*: positions where the cohort's common alternative allele is
at a usable frequency. To choose them, and to estimate each sample's genotype at each, `markers`
(`contamination.rs`) filled two arrays of samples × positions — each sample's reads on the common
alternative allele and its depth, pooled over its read groups. At 8 bytes a sample a position that
is 1.6 GB at 100 samples and 35 GB at 2,169, over 2 million positions (plan §2).

Neither question needs them:

- **choosing a marker** reads only sums across samples — how many samples have depth there, total
  alternative reads, total depth — which can be added up one sample at a time;
- **a sample's genotype** is needed only at the markers (52,525 of 2 million positions on the
  63-sample tomato cohort).

So `markers` now:

1. finds the common alternative allele at each position, as before;
2. **first pass**, one sample at a time: pools that sample's reads into one pair of scratch arrays
   the length of the position list, and adds them into three totals (`covered`,
   `total_alternative`, `total_depth`, the sums as `u64` in sample order);
3. chooses the markers from the totals, with the same tests in the same order (mismapping
   posterior, `MIN_SAMPLES_WITH_DATA`, zero depth, `MIN_FREQUENCY`);
4. **second pass**, one sample at a time: pools the sample again and writes its dosage at each
   marker, with the same formula over the same two counts and the marker's `pooled`;
5. gathers each library's reads at the markers, unchanged.

## Assumptions and deviations

- **The steps became named functions**: `major_alleles`, `pool_one_sample`, `dosage_of`,
  `fill_the_reads_of_each_library`, and `positions_in`; `DROPPED` moved to module level. The plan did not
  prescribe the structure; each pass now calls the same pooling code, which is what makes the two
  passes agree by construction.
- **Each marker's genotype prior is computed once, when the marker is chosen**, and kept beside it
  for the second pass, as the dense layout computed it once a marker. The same call on the same
  `pooled` gives the same bits.
- **The dense `markers` is kept verbatim as `markers_dense` inside the test module**, rather than as
  a `#[cfg(test)]` item beside the new one, so it cannot be reached from the library.

## Memory

The two arrays of samples × positions (35 GB at 2,169 samples) become two scratch arrays of the
position count (u32, 8 bytes a position together) and three totals (24 bytes a position):
**64 MB at 2 million positions, whatever the number of samples** — above the plan's 40 MB
because `covered` is held as a `usize` and the plan counted it smaller; still constant in the
cohort. The common-allele table (20 bytes a position) was already there.

What remains that grows with the cohort is what the fit itself reads, and **it is twice what plan
§6 counted** (corrected after review): each marker holds a dosage a sample (8 bytes) and four
per-library read counts (`alternative`, `depth`, `depth_low`, `depth_high`, 16 bytes a library).
With one library a sample that is 24 bytes a sample a marker:

| markers | dosages | per-library counts | together, at 2,169 samples |
|---|---|---|---|
| 52,525 (the tomato count plan §6 quotes, from before mismapped positions were refused) | 0.91 GB | 1.82 GB | 2.7 GB |
| 31,758 (the tomato count with the current default refusal) | 0.55 GB | 1.10 GB | 1.7 GB |

Plan §7's expected per-sample cost for contamination counts only the dosages (0.4 MB a sample at
50,000 markers); with the counts it is 1.2 MB a sample at 50,000. This is arithmetic, not a
measurement; the simulation (plan §7) measures it.

## Cost in time

Each sample's reads are pooled twice instead of once — a walk over its depth codes and its
non-reference observations. Not measured here.

## Tests added

- `two_passes_choose_the_markers_the_dense_layout_chose` — on a 30-sample panel at three reads a
  position (one library a plant), and on one of two libraries a plant at six reads a position, the
  two-pass markers equal the dense layout's in every field, floating-point fields by their bits, at
  both grains (per read group and per sample). Each panel yields more than 100 markers (asserted).
- `two_passes_agree_at_the_coverage_floor_and_at_mismapped_positions` — a panel at a third of a read
  a position, which the test asserts holds positions covered by exactly `MIN_SAMPLES_WITH_DATA` (8)
  samples and by 7; and a panel where every tenth position has a mismapping posterior of 0.9 (above
  the refusal at 0.5) and every seventh 0.3 (kept, and weighted by it). Fewer markers come back with
  the posteriors than without (asserted).

After review ([review](../reviews/estimation_memory_issue3_2026-09-26.md),
[fixes](../reviews/fixes_applied_2026-09-26_v3.md)):

- `two_passes_agree_on_deep_reads_and_a_second_alternative_allele` — 30 samples at 150 reads a
  position, with every third sample's alternative reads at every fourth position recorded as `G`.
  The test asserts the panel reaches what the shallow ones do not — some position's common
  alternative allele is `G`, some `G` reads are left out as not the common allele's, and some depth
  code stands for a range the cap of 124 clamps — then compares.
- `two_passes_agree_on_a_panel_below_the_coverage_floor` — 1 and 5 samples: no markers either way.
- `a_dosage_whose_weights_all_underflow_is_twice_the_panel_frequency` — 1,000 alternative reads in
  2,000 give 0.6 at a panel frequency of 0.3.
- The comparison destructures `Marker`, so a field added later cannot be skipped.

Checked against defects put in by hand, the file restored byte for byte after each:

- the floor loosened from `<` to `<=`: the floor test fails ("how many markers"); the other passes,
  which is why the floor test exists;
- the second pass not refilling each sample's counts: both tests fail on the dosages;
- (after review) the common-allele filter made always true, and the cap's clamp dropped: each
  fails the deep test only — before it, both passed every test.

The existing contamination tests and **the cross-platform checksum test pass unchanged**.

## Validation

In the dev container, after review: `cargo fmt --check` and `cargo clippy --all-targets
--all-features -- -D warnings` clean; `cargo test --lib -- contamination cross_platform` 40 passed,
0 failed; `cargo test --all-targets --all-features --no-fail-fast` 4,876 passed, 3 failed — the pre-existing failures in `examples/ng_generic_loci_dump.rs` and `examples/ng_ssr_loci_dump.rs`.

## Follow-ups

- Review found an identical-output saving: compute each sample's counts at the markers inside the
  per-library pass instead of pooling every position a second time. Plan §6 fixes the per-library
  pass as unchanged, so it is left for later.
