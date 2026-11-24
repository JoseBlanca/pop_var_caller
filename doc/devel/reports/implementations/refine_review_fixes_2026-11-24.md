# Fixes from the refinement-memory review (2026-11-24)

Branch `refine-review-fixes`, from main 822432bb. The review is
[`refine_memory_2026-11-24.md`](../reviews/refine_memory_2026-11-24.md); it covers the three commits
that cut refinement's memory (790f4a69, 8c9ce22e, d6c543a3). Each finding was checked against the
code before it was fixed; all seven that were asked for were real. Each fix is its own commit, with
a test that fails before it and passes after, except the doc-comment nit.

Measurements are on a laptop (macOS, arm64) that other jobs shared, so times vary by up to half
between repeats; each is given as the runs it came from.

## What changed

| Finding | Commit | What |
|---|---|---|
| — | 87ca6830 | The review added; six files the fixes touch, not formatted on main, run through `rustfmt` alone, so each fix's diff shows only the fix. |
| B1 | b953c408 | A batch of held lines with no site no longer stops a thinned run. |
| M3 | 5fa56d4d | Three-byte read alleles held at their own length. |
| Nit | b4ac94f4 | `write_one_line`'s and `place_of`'s doc comments each on its own function. |
| Mi2 | 58111230 | On resume, a core floor never set on its own follows a new `--refine-window-floor`. |
| Mi1 | c56d96ec | The rule for writing a sample without reads reads every record, walked or not. |
| M2 | baa5bcb9 | The joins and skipped totals are kept as stage state, not made for each line. |
| M1 | e21c3a02 | A window's count includes what the stage holds outside it. |

### B1 — a later contig of filtered lines with thinning on

Compressing a batch of held lines also packs the skipped sites among them. With no site among the
lines, the sites came back as the empty range `0..0` and the code subtracted the first held site's
number from it — above zero on every contig after the first — so a scaffold of a few `LowQual`
records stopped the run. An empty span now packs nothing. Test
`a_later_contig_of_filtered_lines_only_passes_through_with_thinning_on`: 20 sites on one contig,
one `LowQual` line on the next, default density. Before: panic (`range start index
18446744073709551596 out of range`). After: every kept line once, in order.

### M3 — slack in the three-byte read-allele vector

Collected through `Option<Vec<_>>`, the vector grew by doubling: 1,025 entries in room for 2,048,
which the budget does not count. The entries are now checked first and collected from a slice map,
which allocates once at the exact length. The four- and seven-byte forms were already exact; test
`read_alleles_are_held_at_their_length_in_every_form` checks all three (before: tight form
`(1025, 2048)`).

### Mi2 — the core floor on resume

`NeighbourEvidenceDecision` now records `core_floor_given`. `with_window_floor` moves a core
floor that was not given; `with_core_floor` sets it and marks it given. A fresh run calls
`with_core_floor` only when `--refine-core-floor` is given, which plans the same windows as before.
The resume manifest records the flag (`core_floor_given`, optional); a manifest without it reads a
core floor different from the margin floor as given and an equal or absent one as not given, which
is how those runs had it. The resume options are applied by one function, `the_decision_resumed`.

Tests: `a_resumed_run_plans_the_windows_a_fresh_run_with_the_same_options_plans` compares a
resumed decision's sizes with a fresh one's for four option sets (before: core floor 500,000
against 800,000); `a_manifest_says_whether_the_core_floor_was_set_on_its_own` round-trips the
flag and reads both older manifest forms.

### Mi1 — evidence nearby, over every record

`informative_positions` listed only the walked sites, so under thinning a sample whose nearest
reads on one side lay among skipped far records was written `./.`. The rows it reads are kept
unpacked for skipped records (`take_numbers` leaves them), so the list is again made over every
site, at no unpacking cost. Test `a_sample_whose_far_evidence_was_skipped_is_still_written`
(30 samples at four reads a record; sample 0 without reads from record 1,000 to 1,299; core
1,000–1,099; every record beyond 50 kb of the core skipped): before, 0 of sample 0's 100 core
records called thinned against 98 unthinned; after, 75. The 23 still lost come from thinning's
other effect, which templates are chosen. At one read a record the walks fit a 30-kb segment
length, inside the 50-kb near margin, and the rule cannot reach skipped records; it matters only
for samples fitted at 300 kb.

### M2 — per-line rebuilds on the writing thread

`joined` and `skipped_besides_joined` are stage fields: extended on push (a new line is in no last
block), made again once after each window. Under `cfg(test)` the stage checks before every window
plan that the kept values equal a fresh rebuild; all stage tests pass with the check.

Byte-identical output, and time, on a drawn inbred cohort through the stage
(`tmp/probe/`, not committed; 300,000 records 30 bases apart, 10 samples, all cores):

| Density (records a kb) | Before | After | Output hash before = after |
|---|---|---|---|
| 0 (unthinned) | 12.4, 20.0 s | 14.0, 16.8 s | yes |
| 5 (default) | 54.9, 43.2 s | 7.7, 11.9 s | yes |
| 1,000,000 (every record walked) | 54.9, 52.5 s | 21.6, 18.8 s | yes |

At 1,000,000 a kb the output is also the unthinned output, byte for byte, before and after. At 100
samples over 60,000 records (one window) the time is the solve's: 32.0 s before, 35.4 s after,
same hash.

### M1 — what the stage holds outside a window

`WindowBytes` gains `outside`, running totals of what each record holds while no window reads it,
and `WindowPlan` gains `outside_bytes`, included in `bytes`. A window pays it for every record
held before or after it. The stage keeps the totals (`resident_bytes`): a line's held count less
its solve state, and for a site held packed less its numbers and reads too — so a far record a
previous window unpacked and nobody packed again counts as unpacked, which it is. They are
extended on push, recounted from the batch on compression, and recounted after each window, and
checked under `cfg(test)` against a recount before every plan.

The extra term only grows as lines are read, so the planner's rule "over budget on what has been
read stays over" holds.

How much the count moves:

- `a_shrunken_window_counts_the_lines_held_outside_it` (join cohort, 20 samples, no reads, a dense
  stretch at 10 records a kb, 60-kb margins): at a 1-MB budget 3 of 15 windows shrink, and the
  lines outside add at most 3.8% to a shrunken window's count; at full margins at most 0.1% (the
  one line read past the window). With margins allowed down to 30 kb, the smallest budget that fits
  is 1.04 MB, and its smallest window (33-kb margins) counts 3.0% more. With 20-kb floors an 800-kB
  budget stops the stage at a window counted at 840,090 bytes, the lines outside it included.
- `a_window_pays_for_the_records_held_outside_it` (windows unit test): a window that fit exactly
  at full size shrinks to 800-kb core and margins once 10 bytes a record outside it are counted.
- The review's estimate for 2,176 tomato samples at the 500-kb floor was about 16%; that cohort
  holds 46 bytes of reads a sample a record, which the drawn cohorts here do not, so the share here
  is smaller. Not measured on tomato.

Output where the budget does not bind is unchanged: the same drawn cohort at the default 32-GB
budget gives the same hash before and after at densities 0, 5 and 1,000,000 a kb. Where it binds,
the windows move: at 80 MB and 50 MB, unthinned, the hash changes (216,330 → 216,332 and 216,387
→ 216,386 lines written); at 5 a kb, 80 MB is unchanged (the far lines outside are held packed and
cost little) and 50 MB changes.

## Not done, as asked

- **Mi3** (window edges depend on zstd's compressed length): left for the owner; `Cargo.lock`
  pins zstd.
- **The review's open question 2** (1-a-kb packed and not-packed rows): not a code change.

## Tests

- Per finding, the touched modules: `neighbour_evidence` 137 passed after Mi1 (139, with
  the two M1 tests, in the full run below); `read_phase::fragments` 5 passed; the Mi2 tests in `cli::resume_calling`
  and `run::paralog_filter::resume` passed.
- `cargo test --release --no-fail-fast` at e21c3a02: library 5,308 passed, 9 ignored, 2 failed —
  `alignment::left_align_structured::tests::an_offset_past_the_reference_panics` and
  `genetics::tests::lgamma_panics_on_non_positive_in_debug`, which expect a panic that only a
  debug build raises and fail under `--release` on main too; integration tests 62 passed, none
  failed.
- `cargo clippy --release --all-targets`: no warning in a file these commits touch; the 11 it
  reports are in files they do not (`ssr_fit`, `census_fit.rs`, `examples/vcf_genotype_diff.rs`,
  vendored noodles).

## Observations for the owner

- **`a_phase_set_carries_through_short_windows` no longer shrinks windows by its budget.** Its
  2.5-MB budget was set when lines were counted plain; with lines counted compressed, its windows
  count at most about 1.2 MB (measured while building the M1 test), so the test exercises short
  cores from its sizes, not from the budget. It still passes; whether it should be retuned is a
  call for whoever owns it.
- **Resume does not re-apply the hidden `--refine-window-max`.** A fresh run given
  `--refine-window-max` below `--refine-window-floor` lowers the margin floor to it; a resumed run
  given a new `--refine-window-floor` sets the margin floor to that, above the stored margins.
  Hidden options, so not fixed here.
- **Written lines still count their text.** Lines before the core are written and their text
  freed, but `held_bytes` keeps their compressed length, inside the window and now outside it too;
  the count overstates by that much. Unchanged here.
