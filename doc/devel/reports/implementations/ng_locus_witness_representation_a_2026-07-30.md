# ng — the locus witness representation, Milestone A (the vocabulary)

*Implementation report, 2026-07-30. Plan:
[locus_witness_representation.md](../../ng/impl_plan/locus_witness_representation.md) Milestone A.
Design: [spec](../../ng/spec/locus_witness_representation.md) §1, §3.1, §4, §6;
[arch](../../ng/arch/locus_witness_representation.md) §3. Branch `ng-pileup-generator`,
worktree `pop_var_caller-ng-pileup`.*

**Status: Milestone A complete, at Checkpoint A.** Six commits, one per step, no behaviour
change and — apart from one renamed column header — no output byte moved.

---

## 1. Plan

Milestone A is the plan's largest diff and its least interesting one: the crate's observation
and witness vocabulary is renamed before any type or fold changes, so that every later diff
reads as a behaviour change rather than as a rename with a behaviour change hidden inside it.
The six steps were run in plan order, each through implement → review → validate → commit:

| step | rename | commit |
|---|---|---|
| A1 | `ObservedSequence` → `SequenceObservation`; the field `observed_sequences` → `observations` | `e6571c3` |
| A2 | `ReadCoverage` → `ReadWitness`; the field `read_coverage` → `read_witness` | `2a2d361` |
| A3 | `coverage_of` → `witness_of`; `coverage_order` → `witness_order` | `b7c78b5` |
| A4 | the variant `ReadWitness::Observed` → `Partial` | `9dde99f` |
| A5 | `ObservationRow` → `KeyedObservation` (`ObservationKey` keeps its name) | `7da5f49` |
| A6 | the code's own doc comments: "row" and "cell" → "observation" | `4ad1bc1` |

## 2. Assumptions and recorded departures

Nothing here needed a decision the spec or arch had not made. Five choices were the
implementer's and are recorded rather than escalated:

- **A2 renamed three more things of the same concept than the plan enumerated**, all
  behaviour-free: the STR path's `Classified::Observed { coverage }` → `witness` (16 code
  sites in `ssr.rs`, the field holds a `ReadWitness`), `coverage_label` → `witness_label` in
  the four ng example dumps, and `parity.rs`'s `our_coverage` binding. Leaving a field named
  `coverage` holding a `ReadWitness` is the drift this milestone exists to remove.
- **The two research dumps keep their own `coverage` output column.**
  `ng_ssr_cohort_stutter` and `ng_ssr_aligner_bakeoff` print a column literally named
  `coverage`; the plan renames the column in the two *dump tools*, and moving a research
  tool's output is not this step's business.
- **A4 left the generic dump rendering `observed:<offset>+<positions>`.** That string is
  output, and Milestone A moves no output bytes; D4 rewrites the label when the dump starts
  printing the set.
- **A5 touched six sites, not the plan's 26.** The 26 counts `ObservationKey` (which does not
  move) and the two dump tools' own `ObservationRow` structs — TSV rows of a printed table,
  which A6 is explicitly told to leave alone.
- **A6 renamed comments, not identifiers.** `observation_rows`, the local `rows`, and
  `reference_row()` still say "row". Renaming the fold's locals is diff without value while
  C and D are about to rewrite that function.

## 3. Changes made

12 files, plus the four ng example dumps. No source file outside `src/ng/` and `examples/`
was touched, and no test expectation was edited — the only test-file changes are the same
renames and the reflows `cargo fmt` takes when a shorter name lets a call chain fit one line.

The one output change in the whole milestone is the STR dump's column header,
`read_coverage` → `read_witness` (A2), which Checkpoint A predicts by name.

Two hazards were live and both were caught before commit:

- **A blanket `Observed` → `Partial` would have destroyed an unrelated enum.** `ssr.rs` has
  its own `Classified::Observed`, a read-classification variant with nine construction sites.
  A4 ran against `ReadWitness::Observed` / `Self::Partial` and the prose naming the witness
  variant, never a bare `Observed`.
- **A blanket `cell` → `observation` rewrote nine `std::cell::Cell` comments in
  `generator.rs`** into sentences like "nothing else touches the observation while a read is
  being prepared". That file was reverted whole and given two targeted edits instead. The
  same pass had to keep the idiom "two `u16`s in a row", "`begin_segment` twice in a row",
  and `parity.rs`'s "the debug row", which is a row of debug output.

Five sentences in A6 were rewritten rather than substituted, because substitution produced
English no reader wants — "one `SequenceObservation` per observation" became "per accumulated
key", "the observations sitting on `Partial` observations" became "whose witness is
`Partial`", and so on.

## 4. Tests added or updated

**None, by design.** Milestone A's contract is that no test's expectations change; a test
that needed editing would mean the step did more than rename. The suite moved from 2,806
passing tests to 2,806 passing tests.

## 5. Validation results

Run in the project container, from **this worktree's** `scripts/dev.sh` (see §6 — the wrapper
is per-tree). After every one of the six steps:

- `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --lib --bins --tests --examples --all-features`: **2,806 passed, 0 failed.**
- The STR dump on a tomato CRAM (`SRR7279503.p1.bench.cram`, `SL4.0ch01`, 11,318 lines):
  byte-identical, except for A2's single header line.

At the checkpoint:

- `cargo test --release --lib ng::locus_generation`: **275 passed, 0 failed, 1 ignored** —
  including `ng_agrees_with_production_where_production_fabricated_nothing`,
  `ng_emits_the_same_bytes_in_a_second_process` and
  `every_divergence_from_production_is_one_of_the_six_named_classes`.
- The STR dump's 10 fixture tests pass.
- `cargo doc --no-deps --all-features` produces exactly the errors and warnings it produced
  before the milestone (checked against `HEAD~5`); none of them are in the renamed doc links.

**One failure exists and is pre-existing:** `cargo test --all-targets` panics in
`benches/psp_writer_perf.rs:386` with "index out of bounds: the len is 3300000 but the index
is 3300000". It was verified failing on the unmodified `HEAD` before A1, and the bench is
untouched by this work. The per-step gate is therefore
`--lib --bins --tests --examples`.

## 6. Tradeoffs and follow-ups

- **⚠ The container wrapper is per-worktree, and using the wrong one silently tests the wrong
  code.** `scripts/dev.sh` computes its project directory from its own path, so
  `/Users/jose/devel/pop_var_caller/scripts/dev.sh` builds and mounts the **main** worktree
  regardless of where it is invoked from. Run from that wrapper,
  `cargo test --release --lib ng::locus_generation` reports **202** tests and no generator
  tests at all; run from this tree's wrapper it reports **275**. Precondition 2 was only
  confirmable once the difference was found.
- **`arch/locus_generation.md` still uses the old vocabulary** — 8 occurrences of
  `ObservedSequence` / `ReadCoverage`. It is the arch doc the witness arch doc *amends*, and
  the plan's precondition 5 names only the four specs and the witness arch doc, so it is out
  of Milestone A's scope. It now reads stale against the code and should be reconciled.
- **Identifiers still saying "row"**: `observation_rows`, `reference_row()`, the local `rows`
  in `finalise`, and the STR test names `..._merges_the_two_sides_into_one_row` and
  `observed_is_sorted_by_bases_then_coverage`. Cheapest to fold into C/D, which rewrite those
  functions anyway.
- **The variant's deferred note is unchanged.** `ReadWitness::Partial`'s doc still asks to
  revisit the constructor set "when the generic path mints its first run"; D3 mints
  `from_run` and is where that note is discharged.
