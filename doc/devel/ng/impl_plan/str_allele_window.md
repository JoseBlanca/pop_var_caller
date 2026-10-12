# Plan: the repeat-tract fit's allele window reaches past the recorded reads

**Spec:** [`str_allele_window.md`](../spec/str_allele_window.md), approved by the owner on
2026-10-12 with both open decisions settled:

- **§5.1:** the window widens in both fits, the usual climb and the counting fit, before deciding
  whether the counting fit replaces the climb. That decision waits on the climb's cost at ±9.
- **§5.2:** calling splits an edge class's weight equally among the candidate lengths at or beyond
  the edge.

Four steps. The first two change code; the third measures; the fourth records.

## 1. The fit's window is the recorded range plus one

- `ALLELE_SPAN` (`parameter_estimation/joint/ssr_fit.rs`) becomes `RECORDED_OFFSET_RANGE + 1`, and
  its doc comment says why: the outermost class means "this far or further", because the reads
  cannot tell alleles past the recorded range apart. Both fits read it through
  `SsrFitConfig::default()`, so both widen.
- **A guard** that fails the build if the window is not wider than the recorded range.
- Doc comments and tests that quote 13 classes or 91 genotypes as the default are brought up to 19
  and 190.
- **A drawn-class test** (spec §6.1): alleles placed up to 8 repeats short, reads recorded to ±8.
  Fitted with a window of ±6 the slip rate comes out high and the −6 class takes the far alleles;
  at ±9 the slip rate is recovered within its error. Run on the counting fit, which takes
  seconds; on the climb too if it fits in a test's time, otherwise as an ignored measurement.

**Checks:** the touched tests, filtered, in release; clippy.

## 2. Calling gives an edge class's weight to every candidate at or beyond it

- `fill_candidate_weights` (`calling/genotype_prior/seed_ssr.rs`): a candidate at or past `−span`
  takes the `−span` class's weight divided by how many candidates sit at or past `−span`, and the
  same on the long side. Every other candidate is unchanged, and the floor still applies.
- This holds for files written at ±6 too: their end classes also collected every allele that far
  or further (spec §1.1), so the new reading is the truer one for them as well.
- The test that pins the old behaviour,
  `a_candidate_past_the_fits_reach_takes_the_floor_and_not_the_end_class`, is replaced by one
  that pins the split (spec §6.5), including two candidates beyond the edge sharing it.
- `LengthSpectrum::allele_span`'s doc comment, which says the far candidate takes the floor.

**Checks:** the seed tests, filtered, in release; clippy.

## 3. Measure

- **tomato2** (59 accessions, 32 Mb, psp files in `tmp/route2/psps_tomato2_32mb`): the counting fit
  at ±6 and ±9, every class's slip rate side by side. The long classes should fall as they did
  at ±10; the short ones should stay within their errors.
- **The climb's cost at ±9**, on the oracle's four accessions and on tomato2. A whole tomato2 climb
  took 6.9 hours at ±6, so tomato2 is timed on a sample of strata at both widths rather than run
  whole.
- **The oracle** (`scripts/promote_ng_oracle.sh`, the four tomato2 accessions over 160 regions):
  the record count and checksums before and after, and the explanation of what moved.
- The full test suite and a code review, once, since fitted numbers move.

## 4. Record

- The spec's status becomes "built"; its §4 records the two decisions above and §5 keeps only what
  is still open (the crossing bias, a wider recording for long tracts).
- A short report under `doc/devel/reports/implementations/` with the measurements of step 3.
- The oracle's checksums re-recorded, once the move is explained.
