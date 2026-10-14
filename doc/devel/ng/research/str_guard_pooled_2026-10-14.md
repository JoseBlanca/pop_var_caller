# The repeat-tract guard, tested once a tract (2026-10-14)

Branch `str-guard-pooled`, from `main` at `0e905bbb`. A port of the same fix onto `str-fit-c`,
`str-guard-pooled-fit-c`, exists so that kimura can measure it on psps `main` no longer reads.

## The defect

**The guard** leaves a repeat tract out of the slippage fit when too many of its reads differ from
the reference tract length by a number of bases that is not a whole number of repeat units. Such
reads are a sign that the tract is not slipping the way the noise model describes. The spec states
the threshold as a property of the tract
(`doc/devel/ng/spec/parameter_prepass_joint_records.md` §3.3; also `parameter_prepass_ssr.md` §4.1
and §5):

> one non-whole-repeat read in ten of the reads that differ from the reference tract length. Above
> that, the locus is not something this noise model describes

**The code tested it once for every sample's every read group**, and left the tract out if any one
failed. In `gather_strata` (`src/parameter_estimation/joint/ssr_fit.rs`, the loop over samples at
about line 4580 on `main`):

```rust
for (read_group, records) in sample.get(stratum).into_iter().flatten() {
    if records.guard_is_over_threshold(locus) {
        over_guard = true;
    }
```

and `SsrEvidence::guard_is_over_threshold` (`census.rs` lines 781–793 on `main`) counted only that
one section's reads. At about three reads a sample, one non-whole read beside no other differing
read is one in one, so a single stray read in one sample removed the tract for the whole cohort.

Nothing in the specs or the code comments says the test was meant per sample; every statement of it
is per locus or per stratum. **No other code applies the guard**: the census writer records the
guard's reads, the census file stores them, and calling never reads them.

**Before the fix, the guard kept only tracts with no non-whole read at all.** On tomato2, every one
of the 15,077 tracts the old test kept had zero non-whole reads, out of 296,956 reads that differed
from the reference length. At this depth "more than one in ten in some read group" and "at least one
anywhere" are nearly the same rule.

## The fix

Commit `49fa5304`. The guard's two counts — non-whole reads, and reads that differ by whole repeat
units — are summed over every sample and read group at the tract, and the one-in-ten test is made
once, on the sums.

- `SsrEvidence::guard_is_over_threshold` is gone. `SsrEvidence::section_guard_tally(locus)` returns
  one section's two counts as a `GuardTally`, and only `GuardTally::is_over_threshold` carries the
  test, so a per-section helper cannot be read as the tract's test.
- The guard records themselves (which reads, their length differences) are unchanged, and so is the
  psp format.
- `StratumEvidence::guard_reads` counts the same reads as before: every non-whole read in the
  stratum, at tracts left out as well as kept. Its doc said "at tracts that stayed in", which it
  never was; the doc now says what it counts.
- Tests: `one_stray_read_in_one_sample_does_not_remove_the_tract_for_the_cohort` (40 samples, one
  with a single stray read, 1 in 41 over the tract) failed before the fix and passes after;
  `a_tract_over_one_in_ten_across_the_cohort_is_still_removed`;
  `the_guard_counts_add_across_read_groups_before_the_threshold`.

`examples/ng_str_guard_census.rs` reads a cohort's tract evidence once and counts what the guard
leaves out under both tests, then gathers the evidence as `estimate-parameters` does and checks that
the run's own "left out by the guard" count equals its count under the new test.

## tomato2: before and after

59 accessions over 160 regions of 200 kb, about 3 reads a position. **The psps had to be written
again**: those in `pop_var_caller-stutter-multiplier/tmp/route2/psps_tomato2_32mb` are psp format
1.1, which `main` refuses. They were regenerated into `tmp/tomato2/psps` with this branch's
`pvc generate-psps`, the same reference, catalog, regions and CRAMs. The old psps were also counted,
with the port branch, which reads format 1.

| | old test (each read group) | new test (the tract) |
|---|---|---|
| **regenerated psps** — 15,248 tracts in 103 classes | | |
| tracts left out by the guard | 171 | **110** |
| tracts kept | 15,077 | 15,138 |
| classes with at least 8 tracts with reads (enough to fit) | 46 | 47 |
| (tract, sample) pairs with reads | 760,115 | 762,277 |
| **old psps (format 1.1)** — 15,252 tracts | | |
| tracts left out by the guard | 256 | **141** |
| tracts kept | 14,996 | 15,111 |
| classes with at least 8 tracts with reads | 45 | 46 |
| (tract, sample) pairs with reads | 757,740 | 762,081 |

- The run's own evidence line, with the fix, on the regenerated psps: "110 tract(s) left out by the
  guard; 762277 (tract, sample) pairs with reads of the 893142 the kept tracts make (85.3 in every
  100)".
- The old psps reproduce the 256 tracts the earlier tomato2 run reported. The regenerated psps hold
  4 fewer tracts and lose 171 to the old test rather than 256; they were written by a pvc six days
  newer, and I did not trace which change moves that count.
- No tract kept by the old test is left out by the new one.

**On 59 samples the fix is small: 61 tracts come back, 4 in every 1,000.** It is large where the
cohort is large, because the chance that some sample has a stray read grows with the number of
samples. On kimura's 2,169-sample tomato25 the old test left out 45,432 of 86,688 tracts (52 in
100); that cohort has not been measured with the fix (see the last section).

## What the returning tracts look like

Share of non-whole reads among the reads that differ from the reference length, over the whole
tract, on the regenerated psps:

| | tracts | none differs | 0 | up to 0.02 | 0.02–0.05 | 0.05–0.08 | 0.08–0.10 | 0.10–0.20 | 0.20–0.50 | over 0.50 | over all their reads |
|---|---|---|---|---|---|---|---|---|---|---|---|
| kept by both tests | 15,077 | 3,086 | 11,991 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 of 296,956 |
| brought back by the fix | 61 | 0 | 0 | 12 | 26 | 10 | 13 | 0 | 0 | 0 | 70 of 2,252 (3.1 in 100) |
| still left out | 110 | 0 | 0 | 0 | 0 | 0 | 0 | 19 | 36 | 55 | 1,241 of 1,789 (69 in 100) |

- **Most returning tracts are well under the threshold**: 38 of 61 are at or below 5 in 100, and
  pooled they sit at 3.1 in 100 — close to the 1.97 in 100 the SSR spec (§4.1) quotes for tomato as
  the background rate of a tract the model describes.
- **13 of 61 sit between 8 and 10 in 100**, just under the threshold. Those are messy tracts that the
  threshold, as specified, keeps.
- **What is still left out is not borderline**: 91 of 110 tracts are above 2 in 10, and 55 above
  half.
- The old psps give the same picture: 115 tracts come back, 73 of them at or below 5 in 100, 18
  between 8 and 10 in 100, pooled 144 of 4,527 (3.2 in 100).

**The returning tracts carry more differing reads.** They hold 2,252 reads that differ from the
reference length across 61 tracts, 37 a tract, against 297,000 across 15,077 tracts, 20 a tract,
among those both tests keep. A tract with more slipped reads has more chances to show one non-whole
read, so the old test removed slippery tracts preferentially. (Differing reads also include real
non-reference alleles, which this count does not separate from slippage.)

### Do slippage numbers move?

The six classes the fix brings the most tracts back to, each fitted on the tracts kept before and on
those kept now (`FIT_CLASSES=6`, homozygote excess 0.98 for every sample, otherwise the fit's
defaults). Slippage level is the share of reads that slip; ± is the fit's standard error.

| class (period : reference repeats) | tracts with reads, before → now | slippage level before | now | shorter share before → now | fall-off before → now |
|---|---|---|---|---|---|
| 2 : 8 | 154 → 160 | 0.0437 ± 0.0013 | 0.0442 ± 0.0013 | 0.879 → 0.876 | 0.031 → 0.029 |
| 2 : 10 | 77 → 83 | 0.1201 ± 0.0033 | 0.1239 ± 0.0033 | 0.837 → 0.833 | 0.099 → 0.097 |
| 2 : 11 | 63 → 67 | 0.1375 ± 0.0050 | 0.1419 ± 0.0049 | 0.730 → 0.755 | 0.145 → 0.166 |
| 2 : 12 | 44 → 49 | 0.1250 ± 0.0061 | **0.1462 ± 0.0062** | 0.807 → 0.792 | 0.196 → 0.222 |
| 2 : 15 | 26 → 30 | 0.2815 ± 0.0143 | 0.2931 ± 0.0135 | 0.840 → 0.772 | 0.262 → 0.325 |
| 2 : 18 | 26 → 30 | 0.3143 ± 0.0219 | **0.3951 ± 0.0197** | 0.811 → 0.829 | 0.343 → 0.515 |

- **The level rises in all six classes.** In the other four it moves by about one standard error or
  less. In 2:12 it rises by 0.021, about 3.5 standard errors, and in 2:18 by 0.081,
  about 4.
- Adding 4 to 6 tracts to a class of 26 to 154 moves its level that far only if the added tracts
  slip much more than the rest, which matches the previous paragraph: **the old test biased the
  slippage levels of thin classes downwards**, by removing the slippiest tracts.
- The "before" fit is on a subset of the "now" fit's tracts, so the two errors are not independent;
  "standard errors apart" here is a rough guide, not a test.

## The identity oracle

`scripts/promote_ng_oracle.sh` over its usual slice (20 regions, 4 tomato2 samples), run on the
macOS host twice: with `pvc` built from `main` at `0e905bbb`, and with this branch.

- **The fix moves nothing.** All seven checksums are identical between the two runs. Both fits leave
  out the same 7 tracts; with four samples none of them is borderline — pooled, they hold 21
  non-whole reads of 28 that differ, and 6 of the 7 are above one in two.
- **The baseline was not re-recorded.** All seven checksums already differ from
  `scripts/promote_ng_oracle.baseline` on `main` itself, before this change. The baseline was last
  re-recorded at `b4599247` (2026-10-10); the moves come from what `main` merged since — 7,448
  records then, 7,508 now — and are not explained here. They are main's to explain and re-record.
  The digests from both runs:

  ```
  slice: 20 regions, 4 alignments
  67bfe9a17baff53552a410f4d30fe77d  from_alignments.comparable
  67bfe9a17baff53552a410f4d30fe77d  from_psps.comparable
  915d9f3ff2f60881949132d6d8c3efbc  from_alignments.parameters.toml
  4acdda6810237c942d79f14384953fbe  from_alignments.windows.sorted
  0d61963f33f82e691aa8fdef322384bb  from_alignments.windows.histograms
  129492d37836c67cf4b52fc59fa2e2e7  fitted.parameters.toml
  4c06e2f5ff367fbb6926e5e3a63f9471  from_fit.comparable
  ```

## Kimura

`doc/devel/ng/research/str_guard_kimura_prompt.md` asks a session on kimura to build
`ng_str_guard_census` from `str-guard-pooled-fit-c` and run it, without the fit, on the tomato25
psps the 2026-10-13 run used. Those psps are format 1, so the port branch is the one that can read
them. On tomato2's format-1 psps the port reproduces the 256 tracts the earlier run reported, so its
before count is the same test the kimura run made. **Both branches need pushing** before the prompt
can be used.

## Files

- Fix and tests: `src/parameter_estimation/joint/census.rs`, `src/parameter_estimation/joint/ssr_fit.rs`
  (`census_file.rs` only widens a test helper's visibility).
- Counting tool: `examples/ng_str_guard_census.rs`.
- Outputs (not committed): `tmp/tomato2/guard_census.txt`, `tmp/tomato2/guard_census_fit.txt`,
  `tmp/tomato2/guard_census_oldpsps.txt`, `tmp/oracle_before/`, `tmp/oracle_after/`.
