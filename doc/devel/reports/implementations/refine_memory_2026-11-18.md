# Refinement's window memory cut (2026-11-18)

Branch `refine-memory`, from main 08959627.

## Why

On 2026-11-17 the tomato25 run (2,176 samples) stopped in the writing pass: refinement's smallest
window at the start of chromosome 1, 87,050 records, needed 88.3 GB against the 32 GB budget, 466
bytes a sample a record. Measured on 100 of those samples over the same chromosome's first window
(261,277 records, 87 a kb), the budget counted 363 bytes a sample a record:

| Part | Bytes a sample a record |
|---|---|
| The core's sums, each sample's result, and the records' results | 157 (236.6 a core record) |
| Read alleles | 73 |
| Each sample's state over the rounds | 68 |
| The VCF line, counted as it arrived (about 10 held) | 39 |
| Likelihoods and shares | 25 |

## Plan

Cut what the window holds without changing what is computed, except where the owner accepted a
change: the sums and the written posteriors in single precision (owner, 2026-11-18). Positions stay
64-bit (owner: genomes larger than tomato's).

## Changes made

- [`refine.rs`](../../../../src/neighbour_evidence/refine.rs):
  - **The rounds' haplotype store is let go** once the rounds end, before any result is made.
  - **Results are made [`SOLUTIONS_AT_ONCE`] = 256 samples at a time**, each sample's state dropped
    as its result is made and each result as it is copied into the records. The window no longer
    holds every sample's state, every sample's result and the records' results at once; the
    per-sample results are no longer counted (a chunk's, like the walks' arrays).
  - **Sums over the averaged rounds and `RefinedRecord`'s posteriors and background shares are
    `f32`.** A distribution held so is read back with its largest entry set to one less the others
    (`posterior_from_single`), so `GQ` keeps the small entries' seven significant digits rather
    than the rounding of a number near one, and the entries total one as `score_best_genotype`
    asserts.
  - **`RefinedRecord`'s phase is a phase set (0 for none: a phase set is a 1-based position) and
    two `u16` alleles**, 12 bytes a sample, against 32 for `Option<WrittenPhase>`.
  - **A sample's informative positions are listed when its result is made**, not held over the
    rounds (8 bytes a sample a record); round 1 keeps only their count. Positions stay `u64`.
- [`fragments.rs`](../../../../src/read_phase/fragments.rs): **read alleles in four bytes where a
  record's fit** — the fragment as its distance from the sample's first at the record (`u16`), a
  `u8` allele, the Phred — and seven where some sample's span more than 65,535 fragments or an
  allele passes 255. `of_sample` gives a `SampleReadAlleles` view; the link sum reads it through
  `ReadAlleleList` ([`link.rs`](../../../../src/read_phase/link.rs)), the same entries in the same
  order.
- [`stage.rs`](../../../../src/neighbour_evidence/stage.rs): **lines counted at their compressed
  length** over the windows that span them, and at their whole length in a core (handed back
  rewritten and whole). No window is planned on a line not yet compressed: when one could be,
  the uncompressed lines are compressed and the window planned again. Counting a plain line at its
  whole length only overstates, so a window ready on the compressed count is ready on the plain
  one, and the windows are those of counting every line compressed whenever its batch was.
- **`--refine-window-floor <BASES>`** on both calling commands (default 500,000, at most the
  1,000,000-base margin), stated in the `##neighbourEvidence` header as `margin_floor`. The
  over-budget error names it.
- [`neighbour_evidence_pass.md`](../../ng/spec/neighbour_evidence_pass.md) §6 updated.

## Accounting, biallelic, bytes a sample a record

| | Before | After |
|---|---|---|
| Held at every record: haplotypes, followed templates, phase, positions | 2·(8+2) + 24 + 8 = 52 | 44 |
| Core: sums | 7 × 8 = 56 | 28 |
| Core: each sample's result | 24 + 32 + 32 = 88 | not held at once |
| Core: the records' result | 24 + 1 + 8 + 32 = 65 | 12 + 1 + 4 + 12 = 29 |
| Core: the rewritten line | counted in the line | the line's length |
| Read alleles | 8 a sample with any + 7 each | 12 + 4 each where they fit |
| Line, every record | as it arrived | compressed |

## Tests added or changed

- `fragments::tests::read_alleles_are_held_in_four_bytes_where_they_fit`: the edges (65,535 and
  255 fit, one more does not), every entry given back, the bytes counted.
- `refine::tests::a_phased_windows_posteriors_are_the_same_bits_everywhere`: re-recorded for `f32`
  posteriors. The phase digest test did not move.
- `stage::tests` header string: `margin_floor`.

## Validation

- `cargo test --release --lib`: 5,265 passed, 2 failed (the two tests that expect a debug-only
  panic, as on main), 9 ignored. `cargo test --release --test '*'`: 33 passed.

## Measured (100 tomato25 samples, chromosomes 0 and the first window of 1, filter off)

The same psps, before (main) and after, each window's budgeted bytes and solve time logged
(`tmp/refine_scale/cmp_*`; the logging was a measurement-only patch):

| | Before | After |
|---|---|---|
| Budgeted bytes, chromosome 1's first window (261,277 records, 2-Mb core) | 9.48 GB | 4.91 GB |
| Budgeted bytes a sample a record there | 363 | 188 |
| Budgeted bytes, chromosome 0's five windows | 2.8 to 5.5 GB | 1.4 to 2.9 GB |
| Real peak, 48 threads | 34.9 GB | 34.9 GB |
| Real peak, 8 threads | 12.5 GB | 10.9 GB |
| Solve time of the six windows, 48 threads | 147.4 s | 145.8 s |
| Solve time of the six windows, 8 threads | 374 s | 361 s |

At 100 samples the walks' per-thread arrays (about 2 kB a window record a thread, uncounted) make
most of the real peak: at 48 threads it is the same 35 GB either way. The budget falls more than
the real peak because it adds what is never held at once (the rounds' store, then the results):
it is an upper bound, and what moved is that bound, which is what stops a run. At 2,176 samples,
where the per-sample part dominates, the failed window's 466 bytes a sample a record scale to
about 240: about 46 GB for that window, and about 68 GB for an interior window at the floor (about
130,000 records), plus about 14 GB of the 48 threads' arrays.

**Output.** Over the 527,247 records written before the run was stopped, 100 samples each: 6
genotypes of 52.7 million differ, 2 `GQ` (by one), and 157 `GP` in the third decimal.

## Tradeoffs and follow-ups

- The phase's path bits stay 16 bytes a pair: the default 36 paths need more than 32 bits each.
  Packing two 36-bit masks in 9 bytes would save 7 bytes a sample a record.
- The budget still adds the rounds' store and the records' results, though the store is gone
  before the results exist; counting the larger of the two would save about 20 more.
- The walks' per-thread arrays (about 2 kB a window record a thread) are still not counted.
