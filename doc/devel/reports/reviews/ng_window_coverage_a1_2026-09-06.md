# Review — window coverage A1: production's sliding window, transcribed

**Date:** 2026-09-06
**Reviewed:** commit `09bad846` on `ng-window-coverage` (plan step A1 of
[window_coverage.md](../../ng/impl_plan/window_coverage.md))
**Impl report:** [ng_window_coverage_a1_2026-09-06.md](../implementations/ng_window_coverage_a1_2026-09-06.md)
**Per-category files (audit trail):** `tmp/review_2026-09-06_window-coverage-a1/`

## 1. Scope

A PR-shaped diff: a new module and one `mod` line.

- **In scope:** `src/ng/window_coverage/mod.rs`, `accumulator.rs`, `production_parity.rs`, and
  the added line at `src/ng/mod.rs:76`.
- **Out of scope:** `src/sample_summary/coverage.rs`, read by every reviewer as the *oracle* the
  transcription is measured against, never edited (the freeze rule); the commit's doc, plan and
  `PROJECT_STATUS.md` changes, except that the implementation report's own numbers were checked;
  four checks that are red on `main` before this branch — nine files failing `cargo fmt --check`,
  three `needless_lifetimes` clippy errors in `src/ng/run/cohort_merge/`, the failing integration
  test `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`, and the example
  `ng_candidate_selection_probe`, which does not compile.

## 2. Verdict

**Approve with changes.** The transcription is faithful — an independent normalised diff against
production found no undeclared behavioural difference — and the module is where it belongs. What
the review found is almost entirely **untested surface**: nineteen mutations were run against the
twelve committed tests and **seven survived**, none of them because the mutation changed nothing.
Every survivor now has a test.

## 3. Top three

1. **The histogram's four bin-scheme fields were echoed with no test** (B1). They are how a
   consumer turns a cell index back into a depth; three separate wrong echoes left all twelve
   tests green, and a wrong `depth_bin_width` makes every fitted single-copy depth wrong by a
   constant factor with nothing to signal it.
2. **Four behaviours were pinned by the differential against production and by nothing else**
   (M5) — the overflow column's boundary, the GC clamp, lowercase bases, and the closing
   frontier — and the module's own doc schedules the differential's histogram half to stop
   applying at plan step A3, which is two steps away.
3. **Two panicking contracts had no test** (M1, M2). Production has both tests; they were written
   against the tiled accumulator ng does not transcribe, so transcribing only the `sliding_*`
   tests left them behind. A silent `depth_bin_width` of `NaN` sends every window into depth
   bin 0 and produces a plausible-looking histogram from a configuration nobody chose.

## 4. What's good

- **The differential is a real oracle, not a tautology**, and two reviewers established that
  independently: production's accumulator is a separate body of code, and one reviewer
  re-implemented the test's generator in Python and reproduced its window count exactly.
- **The transcription is faithful, checked rather than asserted.** Both bodies were stripped of
  comments, their identifiers renamed to a common alphabet and diffed; `assert_valid` and
  `cell_index` come out character-identical, and every surviving difference maps to one of the
  six deviations the implementation report declares.
- **The `# Invariants` argument holds**, including under repeated positions — which the
  differential's stream produces 1,955 times — because the front is popped on a strict `<`.
- **The freeze holds:** `git diff` over `src/sample_summary/`, `src/paralog/` and
  `src/var_calling/` is empty.

## 5. Findings

Severity codes assigned at synthesis. `Categories:` marks a finding more than one reviewer
raised.

### Blocker

**B1 — `accumulator.rs:175-182`: the histogram's bin scheme is echoed with no test.**
*Categories: reliability, refactor_safety (convergent).* `finish` copies `window_bp`, `gc_bins`,
`depth_bin_width` and `depth_bins` out of the configuration, and nothing reads any of the four —
not the eleven unit tests, not the differential, which compares only cell counts and the folded
count. Three separate mutations to that block (`gc_bins + 1`, `window_bp + 1`,
`depth_bin_width` forced to `1.0`) left all twelve tests green. **Applied:**
`finish_echoes_the_configured_bin_scheme`, whose fixture gives the four fields four different
values and deliberately avoids `depth_bin_width: 1.0` — the reviewer's own first attempt used the
shared test configuration, whose width *is* 1.0, and the mutation survived it.

### Major

**M1 — `mod.rs:64-73`: `assert_valid`, a panicking contract, had no test.** Deleting any of the
four asserts left all twelve tests green. Production has `zero_gc_bins_scheme_panics`
(`coverage.rs:769`), written against the tiled accumulator, so it did not come across.
**Applied:** four `new_panics_on_*` tests, one per assert.

**M2 — `accumulator.rs:123-127`: the coordinate-order guard had no test.** Replacing the
condition with `true` left all twelve green. **Applied:** two `observe_panics_in_debug_on_*`
tests. Whether the guard should also hold in release is a separate question, deferred below.

**M3 — `accumulator.rs`: `sliding_n_position_excluded_but_advances_frontier` could not check the
second half of its own name.** Every covered depth in its fixture is 10, so deferring a
finalisation changes neither a value nor an order. Replaying the differential's 200 streams under
both frontier rules shows they diverge on 120 of them, and the first divergence is a **repeat** at
exactly `centre + half` — a repeated position is the only input that makes the boundary visible,
and none of the eleven has one. **Applied:** the test is renamed to what it does check
(`sliding_n_position_is_excluded_from_every_window`, with the gap named in its doc), and
`a_centre_closes_when_the_frontier_reaches_its_right_edge` pins the frontier.

**M4 — `mod.rs:109-118`: the hand-written `PartialEq` and the absent-value convention it exists
for had no test and no producer.** Replacing it with a derive left all twelve green.
**Applied:** `an_absent_window_equals_itself_and_signed_zeroes_differ`.

**M5 — `production_parity.rs`: four behaviours were caught by the differential alone.** Measured:
the overflow column's boundary, case-insensitive GC, case-insensitive `N`, and the `<=` frontier
each left all eleven unit tests green. The module's own note says the histogram half of the
comparison stops applying at A3. **Applied:** unit tests for all four
(`mean_depth_on_the_top_bin_edge_lands_in_the_overflow_column`,
`gc_fraction_of_one_saturates_into_the_last_gc_bin`,
`lowercase_reference_bases_are_read_the_same_as_uppercase`, and M3's frontier test).

### Minor

**Mi1 — two items were `pub` that nothing outside the module reaches.** *Categories:
module_structure, idiomatic (convergent).* `pub mod accumulator;` gave the one exported type two
public paths, and `assert_valid` was widened from production's `pub(crate)`. **Applied:** both
narrowed.

**Mi2 — `accumulator.rs`: the histogram cell counter wraps in release.** *Categories: errors,
reliability (convergent).* `counts[…] += 1` on a `Vec<u32>`, one count per covered position: a
cell can only reach `u32::MAX` on a reference above about 4.3 Gbp — larger than tomato (0.9) or
human (3.1), not larger than the plant genomes this caller's range commitment names. A wrapped
cell reports a near-empty cell where the single-copy peak is, and the filter's fit anchors on
that mode. Inherited from production, not introduced. **Applied:** `saturating_add`, with the
bound stated on the field's doc. Widening to `u64` would double the histogram from 80.2 kB to
160 kB a sample against spec §3.4's 500 kB budget; deferred to the D3 memory measurement.

**Mi3 — `accumulator.rs`: the order guard is `debug_assert!` while the configuration guard
asserts in release.** **Deferred to plan step C2**, which supplies the caller — see §7.

**Mi4 — `accumulator.rs`: `finish` copied four configuration fields by name**, so a fifth field
would compile clean here while the histogram's doc claim to carry "the bin scheme needed to read
a cell" quietly stopped being true. A2 and A3 both add a field. **Applied:** exhaustive
destructure.

**Mi5 — the two settled values, 500 bases and 50 GC bins, existed only in prose.** The module's
first doc line says "the 500-base window" and no code carried 500; the only configuration
literals were test fixtures, so C2 would have typed four bare numbers. It matters because spec
§3.3's look-ahead is `window_bp / 2`, implemented separately in `cover` at C3, and the spec calls
that failure silent. **Applied:** `WINDOW_BP` and `GC_BINS` with their spec citation; the depth
axis deliberately has no constant, because D1 and D2 measure it.

**Mi6 — three long-lived fields kept production's shorthand.** `count` is the divisor of both
numbers the module produces and sits four lines from `windows_folded`, a different count of a
different thing. **Applied:** `half` → `half_window_bp`, `buf` → `positions`, `count` →
`summed_positions`.

**Mi7 — `production_parity.rs`: a modulus restated an array's length**, so a base added to
broaden coverage would silently never be drawn. **Applied:** the array is named and its length
read from it.

**Mi8 — `mod.rs`: "absent" is spelled as `NaN` in two independently settable public fields**, so a
half-absent pair is representable. Nothing builds one yet; A2 is where two sites start writing
`f32::NAN`. **Deferred to A2** — see §7.

**Mi9 — `accumulator.rs`: thirteen fields, six of which `reset_contig` zeroes by hand.** A future
per-contig quantity added and forgotten there would average across a contig boundary silently.
The reviewer applied a `SlidingWindow` sub-struct extraction and confirmed the differential stays
bit-identical. **Deferred to A3** — see §7.

**Mi10 — `sliding_is_deterministic` cannot fail for any defect this module can have.** No map
iteration order, no clock, no randomness. **Applied:** kept for parity with production, with a
doc line saying what it can and cannot fail on.

**Mi11 — the ramp test's GC assertion cannot fail**, because every base in its stream is `G`.
Its mean-depth half is genuinely discriminating. **Applied:** the doc says so, and M5's lowercase
test carries the mixed-GC case.

**Mi12 — `mod.rs`: a prospective mechanism written in the present tense.** The doc said both
fields *are* `NaN` when absent and that a derived `PartialEq` "would break every test that
compares two runs". Measured: nothing at A1 emits `NaN`, and swapping in `#[derive(PartialEq)]`
left all twelve tests green. **Applied:** the doc now says the convention starts at A2 and why
the impl is written ahead of it.

**Mi13 — the mutation was described as *widening* the window; it narrows it.**
`centre - (half - 1)` moves the left edge one base to the right. The same sentence was in the
implementation report and in `PROJECT_STATUS.md`. **Applied:** corrected in all three, with the
rider that at `window_bp = 1` the literal mutation underflows instead of producing a wrong
window.

**Mi14 — the named reader of `callable_positions` does not read it.** *Categories: extras, naming
(convergent).* `src/paralog/inbreeding.rs` takes the number as a *parameter* and, inside `src/`,
its `obs_het` has no caller at all. The field's actual in-crate reader is
`src/var_calling/diversity.rs:223-235`, which exhaustively destructures the histogram for a
cohort diversity denominator. The substantive half of the claim was checked and holds:
`src/paralog/coverage_model.rs` reads only `window_bp`, `gc_bins`, `depth_bins`,
`depth_bin_width` and `counts`. **Applied:** corrected in `mod.rs` and in the implementation
report.

**Mi15 — the differential never exercised a window anywhere near the configured width.** Widths
were drawn from 1 to 12; the run configures 500. One reviewer widened the sweep and ran it:
widths to 600 alone change nothing (the asserted count is the number of non-`N` covered positions
and is width-independent), and lengthening the contigs so a 500-base window genuinely slides
gives 447,581 windows, all bit-identical to production. **Applied**, and the asserted count
re-measured on this tree.

**Mi16 — `finish` lacked `#[must_use]`** although its doc's whole argument is that the tail
cannot be dropped. **Applied.**

**Mi17 — a doc link labelled a module path with a type name.** **Applied:** the link now names
the type.

**Mi18 — the ready deque's entry is 24 bytes where spec §5 budgets 12.** Not a code defect: spec
§3.6 specifies `pop_ready` returning `(GenomePosition, WindowCoverage)` and ng's `Position` is a
`u64`, so §5's own arithmetic is what does not add up. **Raised with the owner** at Checkpoint A;
plan step D3 measures it.

### Nits

Applied: `depth_cols` → `depth_columns`; the generator's four magic constants named and sourced;
`DeterministicStream::next`/`below`/`seeded` → `next_value`/`next_below`/`from_seed`, since
`next` on a non-`Iterator` reads as `Iterator::next`; `u32::try_from` rather than `as u32` where
the differential feeds production, so a wider fixture fails instead of silently comparing against
a wrapped stream; production's scheme built from ng's configuration rather than typed twice;
`assert_valid`'s doc names C2 as the boundary that supplies a configuration; the memory claim
qualified for streams with one record per position.

Considered and **not** taken, with reasons: `windows_folded` → `windows_binned` (the spec's own
vocabulary is "folded into the histogram", §3.4); an `Eq` impl beside the hand-written
`PartialEq` (production deliberately withholds `Eq`/`Hash` on a float pair, and bitwise equality
is right for comparing runs and wrong for a lookup key); rewriting the shrink-left loop as
`while … is_some_and(…)` (it is production's exact control flow, and line-for-line diffability
against the oracle is worth more at A1).

## 6. Verification

Run in the container on the fixed tree:

- `cargo test --lib --all-features` — **6,301 passed, 0 failed, 15 ignored** (6,275 on `main`
  without this module; 26 tests in it).
- `cargo clippy --lib --all-features` — no diagnostic in `src/ng/window_coverage/`; the three
  pre-existing `needless_lifetimes` warnings elsewhere are unchanged.
- `rustfmt --check` — clean on all three module files.
- **Seven mutations re-run by the orchestrator on the fixed tree, each against the accumulator's
  own unit tests with the differential excluded** — the histogram's `depth_bin_width` echo forced
  to a constant, the overflow column merged into the last regular bin, the frontier's `<=`
  narrowed to `<`, case-sensitive GC, case-sensitive `N`, a derived `PartialEq`, and the
  `window_bp` assert removed. **Each fails exactly one test** (24 passed, 1 failed), so every
  survivor the review found is now carried by a unit test rather than by the differential alone.

**The calling oracle was not re-run for this step and does not need to be**: nothing calls the
module, so no VCF byte can move. It is run either side of A3, which the plan marks as its own
commit.

## 7. Deferred, with a home

- **The order guard in release (Mi3) → plan step C2.** The reviewer's argument rests on C2's
  caller being the parallel cover; each sample's records in fact reach that sample's own
  accumulator in that sample's own coordinate order, so the premise is weaker than stated and the
  decision is better made with the caller in hand. A1 keeps production's `debug_assert!` and says
  so in the doc.
- **A named absent value (Mi8) → plan step A2**, the step that introduces it. A `pub fn absent()`
  with no caller and no producer at A1 would be API written ahead of its first use.
- **The `SlidingWindow` sub-struct (Mi9) → plan step A3.** The hazard is a per-contig field added
  and forgotten, and neither A2's floor (configuration only) nor A3's depth-scale buffer (which
  spans contigs and must *not* be reset per contig) is one. A3 reshapes the struct anyway, and the
  differential's window half guards the extraction whenever it happens — so doing it then is one
  restructure instead of two, and A1 keeps the property it exists for: a file a reader can diff
  line for line against production.
- **`counts` widened to `u64` (Mi2) → plan step D3**, priced against the measured per-sample
  memory rather than guessed.
- **Spec §5's 12-byte deque entry (Mi18) → the owner**, at Checkpoint A.
