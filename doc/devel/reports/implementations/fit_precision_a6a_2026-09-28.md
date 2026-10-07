# Fit precision, step A6 (first part) — each library's reads under its own error rates

**Date:** 2026-09-28. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A6,
first part. **Source:** checkpoint A's decision 1; spec [§3.2](../../ng/spec/fit_precision.md)'s
per-library block. **Branch:** `fit-precision`.

## 1. What was built

The SNP/indel fit now scores each library's reads under that library's own two error rates. A
sample's genotype at a position is still one, shared by its libraries; everything above the
genotype — the allele-frequency density, the classes, the homozygote excess — is unchanged.

Before this step a sample sequenced in two libraries had its reads **added together** at every
position and scored under its **first** library's rates, and the maximisation then credited
**every** library of the sample with those pooled reads, so both libraries came back at one rate
between them. Now:

- **the per-position reader keeps each library apart** (`EvidenceCursor::next_position`): each
  library's counts, depth range and depth weights have their own slot
  (`PositionEvidence::libraries`, indexed by read group), where they were summed into one slot a
  sample;
- **the likelihood reads each library's own rates** (`one_position`): a sample's log-likelihood
  under a genotype is the sum, over its libraries, of that library's reads under that genotype and
  that library's rates;
- **each library's tallies are credited with its own reads only**, so the maximisation fits each
  library's rates from its own evidence (the four branches of `one_position`'s attribution);
- **each library's depth range is weighted around that library's own mean depth**
  (`EvidenceCursor::mean_depth_of_each_library`), where it was the sample's: two libraries at a
  position are two Poisson draws, and their sum's centre is neither's.

**The plan's split, taken.** The step reaches the scores, the information blocks and the errors as
well, which all hold three parameters a sample. Scoping put the whole step above one reviewed
commit, so it is split where the plan says: this part changes the likelihood and moves fitted
numbers; the second part will score each library's rates, give a sample of `k` libraries its
`1 + 2k` parameters in the information and the errors, and move no fitted number. **Until then a
cohort holding a sample read from more than one library gets no standard error at all**, said as a
reason of its own, `AwaitingEachLibrarysScores` ("not computed yet where a sample has several
libraries"). Only the rates themselves lack a slope, but leaving them out of the information makes
every other error come out as if they were known: measured in review, the mismapped share's error
falls to 0.26 of its value at 4 samples and 3 reads, 0.70 at 8 samples, 0.90 at 20. (The first
version of this step marked only those rates; the review caught it.) A cohort of one library a
sample keeps every error.

## 2. What moves, and where

**Only samples read from more than one library.** For a sample of one library every line of new
arithmetic reduces to the old one: its reads are its library's reads, the first library's terms are
taken as they are and later libraries' added (so a sample of one library is scored with exactly
the old operations, in the old order), and its library's mean depth is its sample's mean depth to
the bit. The checks below say so on real data.

- **The identity oracle** (`scripts/promote_ng_oracle.sh`, four tomato accessions at about 3 reads
  a position, one library each, over tomato2's first 20 regions; `tmp/fit_precision/oracle_a6a/`):
  **all seven checksums equal step A5's and main's**, the two fit lines included
  (`fitted.parameters.toml` `c590689d40eb66359b2973c2bd1cbfcd`, `from_fit.comparable`
  `5cdeafd0aaa805cfd57b79e3bbb4e3d2`). The script's exit 1 is the stale baseline's two fit lines
  (A1 report §6), as before.
- **The cross-platform checksums did not move either, and the plan expected them to.** Their
  fixture (`cli::test_fixtures`) has one sample with two read groups, but its synthetic reads carry
  no sequencing error, so every library's two rates sit on their lower bounds whether its reads are
  pooled or not. Measured with the per-pass trace (`tmp/fit_precision/a6a_fixture_trace.tsv`): at
  the last pass of each of the three starts, all three read groups' clean rate is
  1.00043701084885278 × 10⁻⁶ and mismapped rate 1.00000375495820165 × 10⁻⁴, the bounds being 10⁻⁶
  and 10⁻⁴. **So nothing is re-recorded in this commit**, and neither checksum can see this change.
  The only evidence that the change does what it says is the drawn-cohort tests below.

So neither checksum check moves. Numbers move only in a sample read from several libraries, where
each library now gets its own rates instead of one pair between them: on kimura, 482 read groups are a
sample's second or later (the difference of 2,651 read groups and 2,169 samples), so between 483
and 964 read groups belong to samples of several libraries; the log does not count them directly.

## 3. Deviations and assumptions

1. **The step is split** (§1), as the plan allows. The interim reason `AwaitingEachLibrarysScores`
   exists for one commit and goes in the second part.
2. **Each library's depth range is weighted around its own mean depth**, which the plan does not
   name. It is the only consistent reading once a library's reads are scored apart. It moves
   numbers only where a stored depth is a range in a sample of several libraries, which needs a run
   whose depth cap is above 124: the shipped cap is 124 (`src/run/gatherer.rs`), so today it moves
   nothing.
3. **A sample's count of positions with one and with two reads** (`positions_with_reads`,
   `positions_with_two_reads`) now adds its libraries' depths, each the middle of its own range.
   For one library this is the old number; for several it differs from the middle of the summed
   range only where a range was cut at zero.
4. **The cross-platform checksums are unchanged** where the plan expected them to move (§2).
5. **The log line's "N read group(s) after a sample's first" clause is replaced**, for this commit,
   by one saying how many samples hold several libraries and why no error is computed; the second
   part removes it.
7. **A sample holding no ordinary-position section** (a census of repeat tracts only) no longer
   panics at its first read group: it contributes no reads, and its homozygote excess is fitted from
   its genotype prior alone. Reachable, untested, recorded (review Mi3).
8. **The cost, measured in review**: for one library a sample a pass is 2.3 to 3.6% slower than
   before (8 samples, 40,000 positions, 10 passes, 4 threads, median of 5, three interleaved rounds);
   two libraries at 1.5 reads a position cost 8.8% more than one at 3 reads.
6. **The drawing generator gained samples of several libraries** (`draw_cohort_of_libraries`,
   `DrawnLibrary`). A library's read count is drawn before the sample's genotype, as the sample's
   was, so a cohort of one library a sample draws exactly the numbers it drew before: every existing
   fixture is the same cohort (the full suite's measured assertions, the coverage test's seeds).

## 4. Changes

- [fit.rs](../../../../src/parameter_estimation/joint/fit.rs): `SampleAtPosition` renamed
  `LibraryAtPosition`; `PositionEvidence::libraries` (per read group) in place of `samples`;
  `EvidenceCursor` takes `group_index` and each library's coverage
  (`mean_depth_of_each_library`), and `next_position` fills one slot a library; `Scratch::new`
  takes the library count; `one_position` sums each library's reads under its own rates, and its
  four attributions credit each library's tally with its own reads at its own rate; `fit_jointly`
  computes each library's coverage (the per-sample mean stays for the log line) and marks the
  errors of samples of several libraries; the generator (`draw_cohort_of_libraries`,
  `DrawnLibrary`, `DrawnCohort::libraries_of_each_sample`); the helpers the likelihood calls name
  their argument `library`.
- [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs):
  `fill_read_slopes` scores a sample of one library under that library's rates and leaves a sample
  of several at zero; the module doc's "Which parameters a sample carries".
- [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs):
  `StandardError::AwaitingEachLibrarysScores`, `StandardErrors::with_errors_awaiting_each_librarys_scores`;
  `named` gives each read group its sample's row; `described` says why no error is computed.

## 5. Tests

| test | what it shows |
|---|---|
| `a_samples_libraries_are_scored_each_under_its_own_rates` (information.rs), replacing `a_second_read_groups_rates_have_no_slope` | Sample 0's reads split into two libraries: at one rate the log-likelihood equals the pooled one (−1.3207167459 × 10³ both, to ten figures; asserted within 10⁻¹⁰ relative); at different rates the second library's clean rate has a slope of −2,086 and its mismapped rate −72.5 (asserted above 10 in size; the old likelihood left them bit-for-bit unchanged); the sample's rate slopes are zero; its homozygote-excess slope matches a central difference to 7.8 × 10⁻⁸. |
| `a_samples_two_libraries_come_back_at_their_own_rates` (fit.rs) | Six samples of two libraries each, 12,000 positions, six reads a position a library, clean rates drawn two to four times apart, mismapped rates 0.04 and 0.10: every library's clean rate within 25% of the drawn one (worst: 0.00168 against 0.00200), each pair's ratio within a factor 1.5 of the drawn ratio (the pooled fit gives 1), each sample's mismapped rates in the drawn order. |
| `a_librarys_reads_are_scored_under_its_own_rates_and_no_other` (information.rs, from review) | A second library holding no read leaves the likelihood bit-for-bit what the first library alone gives, at two different rate pairs; exchanging two libraries' reads together with their rates leaves it bit-for-bit unchanged (−1.3378684282 × 10³ both). Scoring every library under its sample's last library's rates fails it. |
| `each_librarys_tallies_hold_its_own_reads_in_every_branch` (information.rs, from review) | Over each of 600 positions, with sample 0 split over two libraries and the duplicated class fitted (branch mass: invariant 531.8, fixed 5.4, segregating 59.3, duplicated 3.5), each library's tallies hold its own non-reference and reference reads, to 8.5 × 10⁻¹² relative (asserted below 10⁻⁹). Crediting the segregating or the duplicated branch to a sample's first library fails it. |
| `each_library_is_weighted_around_its_own_depth` (information.rs, from review) | Two samples of a 40-read and a 132-read library under a cap of 140: each library's mean depth within 5% of its draw (40.00, 128.85, 39.75, 129.00), and at all 597 ranged library-positions the weights are `fill_depth_weights` at the library's own mean. Centring on the sample's first library fails it. |
| `a_samples_two_libraries_come_back_at_their_own_rates` (fit.rs), extended in review | Also: each sample's positions with one and with two reads equal the drawn counts (its libraries added), and every error says `AwaitingEachLibrarysScores`. Counting the first library only, or dropping the marking in `fit_jointly`, fails it. |
| `the_logged_line_summarises_each_kind`, `the_trace_names_each_error_after_its_value` (standard_errors.rs) | Updated: a cohort with a sample of two libraries has every error marked, the line says how many samples and why, and every trace error is NaN; a cohort of one library a sample keeps its errors; each read group's trace error is its sample's row, and a read group no sample holds NaN. |

Numbers from `cargo test --release --lib -- a_samples_libraries_are_scored_each_under_its_own_rates
a_samples_two_libraries_come_back_at_their_own_rates --nocapture` (`tmp/fit_precision/a6a_newtests.log`).
`cargo test --release --lib parameter_estimation::joint::fit`: 63 passed, 1 ignored (59 at A5; one
test replaced, four added). The mutations the review found surviving and each test above names as
failing it were re-run against these tests: six mutations, all six fail (`tmp/fit_precision/a6f_mutations.log`).

## 6. Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` (`tmp/fit_precision/suite_a6a.log`):
  4,918 passed, 3 failed, 5 ignored; the three failures are the pre-existing ones
  (`examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1). Both
  `cli::cross_platform_digests` tests pass unchanged (§2).
- The oracle: §2. Its per-pass times, 728.5, 687.9 and 689.3 ms for the three starts, against step
  A5's 776, 710 and 674 ms on the same cohort — both runs beside a full test suite, so rough.

## 7. Follow-ups

- The second part of A6: each library's rate slopes, the `1 + 2k` sample block in the information
  and the errors, the finite-difference test per library, and each library's rate recovered within
  its errors; `AwaitingEachLibrarysScores` and its log clause go.
- For the owner at checkpoint A′: neither checksum check sees this change (§2).
