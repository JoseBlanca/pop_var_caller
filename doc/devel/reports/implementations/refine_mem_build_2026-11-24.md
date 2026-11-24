# Refinement's window memory: five cuts built (2026-11-24)

Branch `refine-mem`, from main 01eb7714. It builds five items of the memory review
`doc/devel/reports/reviews/refine_memory_opportunities_2026-11-24.md`, which is on branch
`mem-hunt` (commit 45f7bc46), not on this one: items 1, 2, 3, 4 and 8 in the review's numbering,
one commit each. The
review's prototype carried hooks for its probe; the changes were written again here without them.
The probe itself is kept as a measuring tool.

**What a reader should take away.**

- Three of the cuts hold less memory, and the output stays byte for byte the same: the walks'
  per-thread arrays (item 1), the packed copy kept beside an unpacked record (item 2), and the
  prediction's sums at every core record (item 3).
- Two change only what the window budget counts, to match what a window holds: the rounds'
  haplotype store and the core's results are never held together (item 4), and a written line no
  longer holds its text (item 8).
- Items 2 and 3 lower the count too. So **every item but 1 changes the count, and where the
  budget binds, windows grow and the output moves**, as for any change of window. Where it does
  not bind, every one of 8 drawn cohorts gave the same bytes as main after each commit.
- On a drawn cohort of 20 samples, dhat's peak of requested bytes fell 26% at one thread (112.2 to
  82.9 MB) and 32% at six (219.8 to 149.6 MB). The window budget's count fell 31% (69.1 to 47.9
  MB). Wall time did not move.

Measurements are on a laptop (macOS, arm64) that other jobs shared, at most six threads.

## What changed

| Review item | Commit | What |
|---|---|---|
| — | d5546eb2 | The probe `examples/refine_mem_probe.rs`, and `NeighbourEvidenceStage::last_window`, the plan of the last window solved, which the probe reads the count from. |
| 1 | bdc4a3cc | The walks' per-thread arrays at what they use. |
| 2 | e848066e | A skipped record held whole drops its packed copy; the count follows. |
| 3 | cb6ef4bc | A phased core keeps one distribution's sums a sample. |
| 4 | 2694761c | The budget counts the rounds' store or the core's results, whichever is larger. |
| 8 | 35796d0d | A written line is counted at nothing. |
| — | (this commit) | Spec §6 and this report. |

### Item 1 — the walks' per-thread arrays

Each thread walking one sample held, at every walked record, the nine drawn paths' haplotypes,
nine state paths, nine paths' uniforms drawn in advance, and the round's input arrays grown by
doubling. Four changes (`src/neighbour_evidence/model.rs`, `refine.rs`, `random.rs`):

- `walk_and_draw_paths` hands every path's haplotypes to a reader the caller gives
  (`on_path`), the round's own first and then the extra paths in the order drawn, and returns only
  the round's own. The phase records each path as it comes, in the order it recorded them before,
  so each path takes the same bit. The eight extra paths share one array.
- A drawn path holds each state as a two-byte code (`State::code`), in the walk's own order of
  states (autozygous `0..K`, then the pairs row by row), not the 24-byte enum. A walk that draws
  asserts at most 255 templates, so every code fits; the pass uses 17.
- Each path reads its uniforms as it draws them, in the order the stream gives them; an extra
  path's stretch of its shared stream is reached with `Uniforms::skip`, which moves a SplitMix64
  state on by `n` steps in one addition.
- The round's likelihoods, shares and drawn alleles, and each path's haplotypes
  (`PathHaplotypes::clear_to_hold`), are allocated at their exact length; the state paths keep the
  capacity that `vec![Vec::with_capacity(n); m]` used to clone away.

The budget does not count these arrays (spec §6, so that windows do not depend on the thread
count), so no window moves.

### Item 2 — one form for a skipped record

With thinned margins a record the far margins skip is packed (zstd) when its line is compressed.
A window that walked it unpacked it and kept the packed bytes, so every core and near-margin record
held both forms, and the count charged both. Now (`src/neighbour_evidence/stage.rs`):

- `HeldSites::hold_whole` drops the packed bytes when it unpacks a site, and packs a site again
  when a window does not walk it. `HeldSites::packed_length` keeps each site's packed size from its
  first packing; zero marks a site never packed.
- **The count follows.** A line's held bytes no longer include its packed size. A running total,
  `packed_bytes` (and `packed_besides_joined`, less the records some sample's last written block
  joins at), holds the skipped records' packed sizes, and `WindowBytes::packed` hands it to the
  planner, which charges it beyond the near margins — where it does not charge the record's whole
  numbers and reads. A record held outside a window is charged in the form it is held in.
- The tests' check of the kept totals, run before every window is planned in every stage test,
  now also checks that every site is in exactly one form and that a packed site holds exactly the
  bytes it is counted at — that is, that packing a record again gives the bytes its first packing
  gave.

### Item 3 — one distribution's sums in a phased core

A core record held, for every sample, the sums over the averaged rounds of its posterior and its
prediction (what the neighbours alone say), and the background's share: `2G + 1` numbers at a
record of `G` genotypes. With phasing, which every run uses, the prediction is read only at the
sample's hidden records, where the fill check reads it and the written posterior is made from it;
elsewhere only the posterior is read. `KeptSums::OneDistribution` keeps the prediction's sums at a
hidden record and the posterior's elsewhere, `G + 1` numbers; a window solved without phasing,
which only tests and research probes do, keeps both (`KeptSums::Both`). The same `f32` additions
in the same order, so the same bits. `core_solve_bytes` counts a phased window's sums.

### Item 4 — the store or the results

The solve lets the rounds' haplotype store go (two rounds' haplotypes and the templates they
followed) before it makes the core's written results. The count added both. Now a window is
counted at the larger: its count is what it was, less the smaller of the store its walks hold and
the results its core writes (`windows::StoreAndResults`). The stage keeps three more running
totals — the store at every record, the store of each skipped record (besides the joined ones, as
the skipped totals are), and the results at every record. `round_store_bytes` and
`written_result_bytes` split the existing per-sample sizes into those two parts. Both are fixed
functions of the records, so windows still do not depend on the thread count; and a record adds at
least its store to the held total and at least its result to the core total, so a window's count
still only grows as it reads further, which the size ladder relies on.

This lowers the count, not the memory. The review recommended it after item 1, so that the
headroom between the count and the real peak, which the per-thread arrays fill, does not shrink;
item 1 halved those arrays.

### Item 8 — written lines

A core's lines are rewritten and written a chunk at a time, their text let go as each is written,
but they stay held for their records in the next windows' left margin, and the count kept each at
its compressed length. Once a core is written, its lines' compressed lengths now come out of the
held running totals. The check of the kept totals now also checks, in an unthinned stage, that
every line is counted at the text it holds plus its record's numbers, reads and solve state; with
the subtraction removed, the stage tests fail that check
(`a_shrunken_window_counts_the_lines_held_outside_it` among them); with it, they pass.

## Verification

**Tests.** After each commit, `cargo test --release --lib neighbour_evidence` passed: 141 tests
after item 1 (one new), 141 after item 2, 142 after item 3 (one new), 143 after items 4 and 8
(one new). They include the recorded digests of posteriors and phase, which hold. New tests:

| Test | Checks |
|---|---|
| `random::tests::skipping_lands_where_drawing_does` | `Uniforms::skip(n)` lands where `n` draws do, at 0, 1, 7 and 1,000. |
| `model::tests::a_state_code_gives_the_state_back` | Every state's code is its place in the walk's order and decodes to it, at 1, 17 and 255 templates. |
| `refine::tests::a_phased_core_keeps_one_distributions_sums` | A phased layout holds `G + 1` numbers a core site, an unphased one `2G + 1`; `core_solve_bytes(2)` counts the phased sums. |
| `windows::tests::a_window_pays_the_larger_of_store_and_results` | With a store of 30 bytes a record over 3,500 records and results of 40 over 2,000 core records, the results go uncounted; with a store of 10, the store does. |

`windows::tests::a_thinned_window_pays_skipped_solve_state_only_near_its_core` gained a case with
packed bytes, charged beyond the near margins only. `extra_paths_leave_the_round_unchanged` and
`reads_that_link_nothing_change_no_bit` read the paths through the new reader. One stage test
changed its budget: `a_shrunken_window_counts_the_lines_held_outside_it` needs windows over its
dense stretch to shrink; they now count 0.97 MB, not about 1.2 MB, so its budget went from 1 MB to
0.8 MB.

At the end, `cargo test --release --no-fail-fast` over every target: 5,374 passed, 2 failed, 9
ignored. The two failures are the tests that expect a panic only debug builds raise
(`alignment::left_align_structured::tests::an_offset_past_the_reference_panics`,
`genetics::tests::lgamma_panics_on_non_positive_in_debug`), which fail the same way on main under
`--release`.

`cargo clippy --release --all-targets`: no warning in any file this branch touches. The warnings
it prints (in `parameter_estimation::joint::ssr_fit`, `run::census_fit`,
`examples/vcf_genotype_diff.rs` and the vendored `noodles-cram`) are in files the branch does not
change.

**Output, byte for byte.** The probe draws a multi-allelic cohort (2 to 6 alleles, 43 records in 100
with more than two), reads as 300-base fragments, and hashes every line the stage writes. It was
run with the budget unbounded, so the windows are the same in every build, on main (with only the
probe's accessor added) and after each commit, at six threads:

| Samples | Contig | Records a kb | Reads a position | Core / margins | Thinning | Hash, every build |
|---|---|---|---|---|---|---|
| 1 | 300 kb | 160 | 3 | 100 / 100 kb | 5 a kb | 035d441e1f74fb4c |
| 2 | 300 kb | 160 | 1 | 100 / 100 kb | off | 865ced59f5b42e2f |
| 50 | 300 kb | 160 | 3 | 100 / 100 kb | 5 a kb | b7e1865d5a6e9d1f |
| 50 | 300 kb | 160 | 30 | 100 / 100 kb | off | 7d658716d08e3166 |
| 50 | 300 kb | 160 | 30 | 100 / 100 kb | 5 a kb | 5654ab5f6402e3a7 |
| 50 | 600 kb | 20 | 1 | 150 / 150 kb | 5 a kb | 6c89effa10fcd690 |
| 300 | 100 kb | 160 | 3 | 30 / 30 kb | 5 a kb | 609791f937c27c00 |
| 300 | 150 kb | 40 | 30 | 50 / 50 kb | off | cc9682d450f9191b |

The pass does not run on one sample in a real run (`decide` turns it off); the probe forces the
stage on, so the one-sample row exercises the code, not a configuration a user meets. The dhat and
resident-memory runs below also matched main's hash at every thread count.

**Where the budget binds** the output moves, as expected. On the 50-sample, 3-read, thinned row
with a 200-MB budget and a 30-kb floor, main planned seven windows (cores of 34 to 100 kb) and the
branch four (55 to 100 kb); the hashes differ, and 49,567 against 49,563 lines were written (a
record is dropped when no sample keeps an alternative allele).

## What the budget counts

The largest window's count on each drawn row, unbounded budget, at each commit (MB):

| Row | main | item 1 | item 2 | item 3 | item 4 | item 8 |
|---|---|---|---|---|---|---|
| 1 sample, 3 reads, thinned | 10.7 | 10.7 | 8.7 | 8.3 | 7.6 | 7.6 |
| 2 samples, 1 read | 14.8 | 14.8 | 14.8 | 14.0 | 12.6 | 12.6 |
| 50, 3 reads, thinned | 280.8 | 280.8 | 257.3 | 237.1 | 202.8 | 202.8 |
| 50, 30 reads | 433.8 | 433.8 | 433.8 | 413.5 | 379.1 | 378.1 |
| 50, 30 reads, thinned | 522.2 | 522.2 | 425.9 | 405.6 | 371.2 | 363.0 |
| 50, 20 a kb, 1 read, thinned | 56.0 | 56.0 | 53.6 | 50.0 | 43.7 | 42.7 |
| 300, 3 reads, thinned | 776.3 | 776.3 | 711.3 | 675.0 | 613.2 | 606.1 |
| 300, 30 reads | 312.9 | 312.9 | 312.9 | 298.3 | 273.4 | 273.4 |

Item 2 moves nothing without thinning, as nothing is packed; it moves most at 30 reads, where read
alleles are most of a packed record (522 to 426 MB). Item 8 moves a window only once some core has
been written before it, so the largest window, often the first, does not always show it.

**At 2,000 samples, computed** for the review's case (a window at the 500-kb floor, 160 records a
kb, tomato's mix of allele counts, margins thinned to 5 a kb: 80,000 core records, 100,500 walked,
93,000 core and near-margin records held whole), from the per-record sizes in the review's §2:

| Item | What the count no longer charges | At 2,000 samples |
|---|---|---|
| 2 | 27 bytes a sample at each record held whole | 5.0 GB |
| 3 | 24.2 bytes a sample at each core record (`G` floats, 6.06 genotypes on average) | 3.9 GB |
| 4 | the smaller of the store (33.1 bytes a sample a walked record, 6.65 GB) and the results (41.2 a sample a core record, 6.59 GB) | 6.6 GB |
| 8 | 8.1 bytes a sample at each written left-margin line | 1.3 GB |

About 17 GB of the review's 61-GB count; items 2 and 3 are held memory too, items 4 and 8 count
only. The review's real peak during the rounds was about 1.09 times the count before these
changes; that ratio was not measured again, on drawn or real data.

## Memory held

**Requested bytes (dhat).** 20 samples, 100 kb at 160 records a kb, 3 reads, a 50-kb core with
50-kb margins, thinned to 5 a kb: 16,076 lines written, about 16,700 records walked in the largest
window. Built with `--no-default-features --features dhat-heap`.

| Threads | main | branch | Change |
|---|---|---|---|
| 1 | 112.2 MB | 82.9 MB | −26% |
| 2 | 145.7 MB | 100.8 MB | −31% |
| 6 | 219.8 MB | 149.6 MB | −32% |
| bytes asked for in all | 13.4 GB | 6.2 GB | −54% |

Each thread beyond the first added 21.5 MB on main and 13.3 MB on the branch, averaged from one to
six threads; from one thread to two, 33.5 against 17.9 MB, about 2.0 against 1.07 kB a walked
record over the review's count of 16,700 walked records for this cohort — the review's 2.15 and
1.10 kB. At the review's 48 threads over its 100,500 walked
records that is about 10 GB less 5 GB, a thread's arrays not counted by the budget.

**Resident memory and wall time (mimalloc, `/usr/bin/time -l`).** 48 samples, 300 kb, 3 reads,
100-kb core and margins, thinned; and 48 samples, 200 kb, 30 reads, the same windows. Runs
alternated main and branch.

| Cohort | Threads | Peak resident, main | Peak resident, branch | Wall, main | Wall, branch |
|---|---|---|---|---|---|
| 48 samples, 3 reads | 2 | 631, 659 MB | 545, 519 MB | 38.4, 38.8 s | 38.3, 38.8 s |
| 48 samples, 3 reads | 6 | 1,064, 1,043 MB | 790, 772 MB | 17.5, 17.7 s | 17.7, 17.2 s |
| 48 samples, 30 reads | 6 | 1,150, 1,139 MB | 842, 859 MB | 14.0, 14.0 s | 13.9, 13.7 s |

The peak fell by 14–19% at two threads and 25–27% at six; the hashes matched main's in every
pair. Wall time moved by less than the spread between repeats, so item 2's extra packing of each
core record, once, costs nothing measurable here.

## Assumptions and choices

- **`walk_and_draw_paths` hands the round's own path to the reader too**, first. The phase records
  the own path's bit before the extra paths', and handing it on keeps that order without the
  caller recording it after the walk returns.
- **Packing again is assumed to give the same bytes** — zstd at one level over the same plain
  bytes, which `RecordNumbers` and `PackedRecordFragments` write back exactly as they read them.
  The count uses the first packing's size whatever happens, so the windows would not move even if
  it did not; the stage tests check that it does.
- **The research probes that plan windows (`ng_neighbour_pass_probe`, `refine_replay`) do not use
  `StoreAndResults` or the packed totals**: they count as before, both the store and the results.
  `ng_neighbour_pass_probe` solves without phasing, which keeps both distributions' sums, so the
  lowered `core_solve_bytes` undercounts it by `G` floats a sample a core record; it plans with an
  unbounded or test budget.
- **Item 4 is applied to window totals, not per record.** Counting the larger part at each core
  record would be simpler but saves less (about 5.3 against 6.6 GB in the review's case), since
  the walked records outside the core hold store and no results.

## For the owner

- **`--refine-window-bytes` was set from the old ratio of real memory to count** (16–22% above the
  count on tomato, 2026-11-23). Items 4 and 8 lower the count by about 8 GB at 2,000 samples
  without lowering what is held, so a run at the same budget now plans larger windows that hold
  more. The review's estimate is that the rounds' peak sits about at the new count plus the
  per-thread arrays and allocator slack; that should be measured on the tomato replay before the
  default budget is trusted at the same value.
- **Items 2 and 3 are lossless, items 4 and 8 change windows wherever the budget binds** — for
  tomato25's dense stretches, which is where it binds, genotypes will move as they do for any
  window change.
