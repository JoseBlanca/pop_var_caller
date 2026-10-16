# The oracle re-recorded after D7 and D8 (2026-10-16)

Plan [`str_penalised_fit.md`](../../ng/impl_plan/str_penalised_fit.md) says D4's re-recording is
repeated after D7 and D8. This report does that, the way
[D4's report](str_penalised_fit_D4_2026-10-16.md) did. Branch `str-penalised-fit` at `909266ee`
(main, up to `98230b01`, merged in). Measured on macOS, arm64.

The identity oracle (`scripts/promote_ng_oracle.sh`) calls four tomato2 accessions over 20 regions
with default parameters, fits the parameters on them, and calls again with the fitted file; its
checksums sit in `scripts/promote_ng_oracle.baseline`.

## 1. The answer

**Only standard errors and notes moved, and each move is pinned to one commit.** No fitted number
moved, and no VCF moved — neither the two called with default parameters nor the one called with the
fitted file. Two of the seven checksums change, both parameters files; the baseline is re-recorded.

| checksum | D4 baseline | now | moved by |
|---|---|---|---|
| default-parameter VCFs (both routes), windows, histograms | — | unchanged | — |
| `from_alignments.parameters.toml` (the defaults run's file) | `4f272d7e…` | `bd7b9bea…` | D8: two note paragraphs |
| `fitted.parameters.toml` | `5ed794a0…` | `c84ce763…` | D7 and D8: error keys; D8: the same two notes |
| `from_fit.comparable` (calls made with the fitted file) | `61b09b60…` | unchanged | — |

## 2. How each commit since D4 was isolated

Five things changed the code between D4 (`932fc9b0`) and now: D7 (`680435b2`), E1, the climb deleted
(`3198d55d`), D8 (`af6c71fd`), the curve fitting deleted (`2276c93c`), and main's eight commits on
STR bundle-caller performance (`0ffd77ae..98230b01`). Each was separated by a fit on the oracle's own
psps on either side of it:

- **D7 against E1.** `680435b2` was built from a `git archive` copy (`tmp/oracle_E/d7_src`, its own
  target directory) and refitted on today's psps. Its fitted file is byte-identical to the one D8's
  work recorded from `c50d2ce6`, which includes E1 (`tmp/d8/fitted_head.parameters.toml`). So E1
  moved nothing, and every move between D4 and `c50d2ce6` is D7's.
- **D8.** `tmp/d8/fitted_d8.parameters.toml`, D8's tree, against `c50d2ce6`'s.
- **The curve deletion and main.** Today's fitted file is byte-identical to D8's. And every VCF,
  window and histogram checksum is D4's, measured before main's performance commits were merged, so
  those commits are output-identical here, as they claimed. No build of main was needed.

## 3. What moved

### 3.1 The fitted file

391 lines that are not comments in both files; 22 differ, all slippage rows (a slippage row holds one
repeat class's three stutter numbers: the share of reads that slip, the share of slips that are
shorter, and how fast longer slips fall off). Every one of the 22 is identical once its
standard-error keys are removed from both sides (`tmp/oracle_E/fit_moves.py`). No allele-length row
moved, F's errors included.

**D7: 15 rows.** D7 adds the allele-length prior's curvature to the matrix the errors are read from.

- The six two-base classes at 6, 7, 8, 9, 10 and 12 repeats **gain all three errors** (18 keys; e.g.
  6 repeats: slip rate 0.0020, shorter share 0.055, fall-off 0.041). At D4 that motif length's matrix
  could not be inverted, and the D4 fit's log says so: "motif length 2, standard errors: none
  determined for any of its numbers except F" (`tmp/oracle_D4/fit.log`, line 46). This is the
  failure D7 was written to fix.
- The 27 errors of nine rows that had them already (one-base classes at 8 to 15 repeats, the
  three-base class at 6) move by a factor between 0.9946 and 1.0019, median 0.9999.
- F's errors do not move because this fit reads F off pairs of samples with a curvature of its own,
  not off the matrix D7 changed (the same log line, and `penalised_fit.rs` around line 1950).

**D8: 7 rows.** D8 lets the shorter share and the fall-off carry an error from 3 loci with reads,
where all three numbers needed 8. Seven classes gain those two errors (14 keys): one-base at 18
repeats, two-base at 11, 13, 14, 15 and 16, three-base at 8. That is the "14 errors, 7 classes with
3 or 4 loci" D8's commit measured on the same psps. Their slip rates still carry none.

**Notes (D8), in both parameters files.** The paragraph on own-fit errors gains "Nor does a number
fitted on too few of the stratum's loci with reads … the first number below 8 loci, the other two
below 3", and the fallback-F paragraph now says the 8-loci count withholds the error on "its
between-locus share or its `share_of_reads_that_slip`" where it said "own-fit standard errors". These
are the only changes in the defaults run's file.

### 3.2 Why the calls made with the fitted file do not move

Calling reads the fitted numbers, not their errors. The errors are loaded with the file and carried
beside the numbers, but nothing in a genotype likelihood, QUAL or GQ is computed from them; the
checksum of the VCF called with the new file is D4's, `61b09b60…`, which confirms it on this data.

## 4. The new baseline lines

```
bd7b9bea209a3c754e19deea7dcfa87d  from_alignments.parameters.toml
c84ce763024442ef439033899a57f2fd  fitted.parameters.toml
```

The other five lines are unchanged. The baseline carries a paragraph naming these moves, and every
line of `tmp/oracle_E/digests.txt` matches it.

## 5. The suite

- `cargo test --release --features bench-fixtures --no-fail-fast`: 5,243 passed, 3 failed, 7 ignored.
  The three failures are the tests of a panic only debug builds raise
  (`lgamma_panics_on_non_positive_in_debug` and the two `an_offset_past_the_reference_panics`), as
  expected in release. The four integration tests in `tests/` and the doc tests ran in the same
  command and passed. (D4 counted 5,362; E1 and the curve deletion removed the climb's tests since.)
- `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo fmt --check`: clean.

## 6. Where the evidence is

All under `tmp/oracle_E/` unless named: the oracle run (log `tmp/oracle_E.log`); `fit_moves.py` and
its outputs `fit_moves_D4_vs_E.txt`, `fit_moves_D4_vs_D7.txt`, `fit_moves_D7_vs_D8.txt`,
`fit_moves_D8_vs_E.txt`; the D7 build (`d7_src/`, `d7_target/`, `d7_build.log`) and its fit
(`fitted_d7.parameters.toml`, `fit_d7.log`); `suite.log`, `clippy.log`, `fmt.log`.
