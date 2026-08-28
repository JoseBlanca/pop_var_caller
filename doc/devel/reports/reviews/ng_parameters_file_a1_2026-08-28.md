# Code Review: ng parameters file — A1, the file's Rust shape

**Date:** 2026-08-28
**Reviewer:** rust-code-review skill (orchestrator) over five category agents in isolated worktrees
**Scope:** the uncommitted working-tree change of step A1 of [parameters_file.md](../../ng/impl_plan/parameters_file.md)
**Status:** Request-changes (all applied — see [the fixes report](fixes_applied_ng_parameters_file_a1_2026-08-28.md))

---

## 1. Scope

- **What was reviewed:** the working-tree diff of step A1 — a new module holding the parameters
  file's Rust shape, with `serde` derives and no reading or writing.
- **Reviewed against:** branch `ng-parameters-file`, uncommitted, on top of `main` at `a6e8472b`.
  Each agent received the diff as a patch and applied it to its own worktree.
- **In-scope files:** [src/ng/calling/parameters_file/mod.rs](../../../../src/ng/calling/parameters_file/mod.rs)
  (new, 960 lines at review time), [src/ng/calling/mod.rs](../../../../src/ng/calling/mod.rs)
  (one `pub mod` line).
- **Deliberately out of scope:** every reuse target the step reads and does not touch —
  `run_parameters.rs`, `stratum_fits.rs`, `ssr_fit.rs`, `slippage_curve.rs`, `share_curve.rs`,
  `census.rs`. Also the work the plan assigns to A2 (the shared value+warrant+count shape), A3
  (spec §5's five states as `Option`), B (the writer) and C (the reader and the real-fit round
  trip); agents were told to file any of it as a **scheduling** finding rather than as a defect,
  and three did.
- **Categories dispatched:** `reliability` (always), `naming` (**the step's central deliverable** —
  the TOML key names are the coder's proposal by the owner's decision of 2026-08-28),
  `idiomatic` (a serde data model that must survive a TOML round trip exactly), `smells`
  (twelve types mirror upstream ones; illegal states), `module_structure` (a new module in a
  folder whose charter does not obviously cover it). Not dispatched, with reasons: `errors` and
  `defaults` (no error path and no default-acting value in a step that is types only),
  `unsafe_concurrency` (none of its triggers), `tooling` (no `Cargo.toml` change), `extras` (no
  parser yet).

## 2. Verdict

**Request-changes.** One Blocker and thirteen Majors, and the Blocker is that the step's own
tests could not fail on most of what they claimed to hold.

## 3. Execution status

Run by the orchestrator in the container, on the reviewed tree, and passed to every agent so none
re-ran them:

| command | result |
|---|---|
| `cargo fmt -- --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib ng::calling::parameters_file` | 4 passed, 0 failed |
| `cargo doc --no-deps` | fails on the tree, on 25 unresolved intra-doc links in other modules; none in this one |

**Three of the five agents lost their isolation worktree before they could write their findings
file** — `smells` partway through its fifth mutation, `idiomatic` and `reliability` during probe
builds. All three returned their findings in full and the orchestrator wrote them to
`tmp/review_2026-08-28_parameters-file-a1/`. Every claim in them was made before the worktree
vanished, and every write those agents made was inside their own worktree. The orchestrator
re-checked the three claims that drove the largest changes against the source directly: the
`u32`/`u64` split at one of three `reference_repeats`, the two `usize` fields, and the fixture
that wrote a reach with no curve.

Findings labelled "Needs verification": one — the `idiomatic` agent could not run its
inline-form parse probe, so the Major about that claim is filed at Medium confidence. It is now
covered by a test.

## 4. Open questions and assumptions

1. **No step of this plan owns range checking.** Seven doc comments state a constraint the type
   permits violating — a dense `0..n` read-group axis, at least one sample, an inbreeding
   coefficient in `[0, 1)`, a curve weight in `[0, 1]`, a length spectrum that is odd and sums to
   one, a substitution rate that is a probability, a contamination table that is all-or-nothing.
   A3 owns `Option`, C owns reading; **neither owns refusing a value outside its stated range**,
   and spec §9 promises "a malformed file fails at read with a line number". Affects the
   `reliability` and `idiomatic` Majors on invariants. **Raised at Checkpoint A rather than
   settled here.**
2. **Whether `LevelProvenance::slipped_reads` and `SharesProvenance::slipped_reads` are one
   number.** Upstream their docs are near-identical and their absence conditions are worded
   differently. If they are one, the file writes it twice and a hand-edited file can make the two
   disagree. Affects the `smells` and `idiomatic` findings on the duplicated count. **The
   measurement that settles it is step C4's round trip on a real fit**, and the field now says so.

## 5. Top 3 priorities

1. **The suite killed 1 of 14 mutations** (`reliability`, Blocker). The only whole-shape test
   writes and reads through the same `serde` derive, so every renaming moves both sides at once
   and stays invisible. Dropping `rename_all` from five of the seven enums, and renaming any of
   about sixty struct fields, each changed the file on disk and left every test green.
2. **The module's own fixture wrote an illegal state** (`smells`, Major): a `curve: None` beside a
   `reach: Some(BelowFitted)` — a claim about a curve's fitted range in a row that records no
   curve. Three of four combinations of those two `Option`s are meaningless and all four were
   writable.
3. **A typo in a hand-edited file read back as a fitted fact** (`idiomatic`, Major): with no
   `deny_unknown_fields`, misspelling the `curve` table by one letter parses and yields
   `curve: None`, which the module's own doc defines as *this stratum's period had no curve*.
   `deny_unknown_fields` is a type attribute with no call-site knob, so this could not be
   deferred to the reader in step C.

## 6. Findings

**1 Blocker, 13 Majors and 32 Minors as filed** across the five agents, plus nits. The
per-finding text lives in the five audit-trail files under
`tmp/review_2026-08-28_parameters-file-a1/`, which are gitignored; this section indexes them by
severity and names the convergences, and
[the fixes report](fixes_applied_ng_parameters_file_a1_2026-08-28.md) accounts for every one.

### Blocker

- **B1** `parameters_file/mod.rs:895` — the sole whole-shape test cannot fail on two of the three
  properties its doc comment claims. *(reliability; convergent with `smells` and `idiomatic`
  cross-category notes and with `module_structure`'s first Minor.)*

### Major

- **M1** the two batching axes are interchangeable, and the fixture — both `vec![0, 1]` — makes
  exchanging them invisible even in the emitted file. It also **falsifies the fixture's own doc
  comment** that "every numeric value is distinct". *(reliability)*
- **M2** five of the seven file-facing enums have no spelling test. *(reliability, convergent:
  `module_structure`'s first Minor, `smells` and `idiomatic` cross-category)*
- **M3** `stated_concentration` carries no warrant, so the run's own fitted median and the stated
  constant — both able to be exactly 1.0 — are the same line in the file. Spec §8 names it among
  the three numbers that must be marked `Defaulted`. *(reliability)*
- **M4** seven documented invariants the types permit violating, on values that arrive from a
  hand-edited file. *(reliability, convergent: `idiomatic` Minor)* — see open question 1.
- **M5** no property-based test on a serializer/deserializer pair whose contract is a round-trip
  law. *(reliability)*
- **M6** nothing checks the claim steps B and C meet on: that the derived reader parses the
  inline-table form the module doc advertises and the derived *writer* never produces.
  *(reliability, Medium confidence)*
- **M7** `an_absent_curve_writes_no_key` passes for an implementation that never writes the key at
  all, and its substring assertion would fire spuriously on a variant spelled `its_periods_curve`.
  *(reliability, convergent: `idiomatic` cross-category)*
- **M8** twelve mirrored types with nothing that makes upstream drift a compile error. *(smells)*
- **M9** `curve` and `reach` are two independent `Option`s where one is meant, and the fixture
  writes the illegal pair. *(smells)*
- **M10** the same slipped-read count is written twice per stratum and the two can disagree.
  *(smells, Medium confidence, convergent: `idiomatic` Minor)* — see open question 2.
- **M11** no `deny_unknown_fields`: an unknown key is discarded in silence, and on an `Option`
  field a one-letter typo reads back as absence. *(idiomatic, convergent: `smells` and
  `reliability`)*
- **M12** `fitted_from` names three unrelated things — a section, whose reads a contamination
  fraction came from, and a repeat count. *(naming, convergent: `idiomatic` Minor)*
- **M13** `warrant` names both spec §1.4's four-state vocabulary word and three structures that
  hold none of the four. *(naming)*

### Minor

Thirty-two, indexed by theme. The ones that changed the code are listed in the fixes report; the
rest are recorded there as deferred with an owner.

- **Widths and types:** `reference_repeats` is `u64` in two rows and `u32` in a third *(four
  categories, convergent)*; `usize` in a persisted format *(three, convergent)*; bare primitives
  where a file-local newtype would cost nothing in the emitted bytes *(smells, idiomatic)*;
  `Copy` on a 208-byte struct and not on the 352-byte one *(idiomatic)*.
- **Names:** `source` meaning two different things in a parent block and its child; `scale` saying
  neither what it scales nor which way; five spellings of "one value per axis"; two positional
  batch arrays contradicting the module doc's own rule; the `InFile` suffix on three of seven
  twins *(three categories, convergent)*; `*Row` on two types that are not rows *(two,
  convergent)*; `regime`, `weights`, `observations`, `declared`, `coefficient`, `term`,
  `CensusBinding`, `FittedFromInputs` *(naming)*.
- **Documentation accuracy:** the round-trip test's doc claims a property `toml` 1.1 does not have
  *(idiomatic, reliability)*; `serde` emits a struct's table-valued fields last, so the `Blend`
  row's `source` is written after its siblings and a hand-written writer matching serde must do
  the same *(idiomatic)*; nothing records why the module lives under `calling/`
  *(module_structure)*.
- **Structure:** the 32 `pub` types *(module_structure — narrowing tested and **not available at
  A1**: per-type `pub(crate)` emits `private_interfaces` warnings that `-D warnings` makes fatal,
  and `pub(crate) mod` fails clippy with 32 dead-code errors, because nothing outside the tests
  constructs these yet)*.

### Nits

Grouped and not enumerated: about fifteen, mostly single names that omit a subject or a unit.
`smells` recorded that the module has no `TODO`, no `#[allow]`, no commented-out code, no dead
code and no boolean parameter, and that every `pub` item carries a doc comment that says *why*
rather than *what*.

## 7. Out of scope observations

- `cargo doc --no-deps` fails on this tree, on 25 unresolved intra-doc links in modules this step
  does not touch. Pre-existing; worth a sweep of its own.
- `module_structure` checked and did **not** file: the module's location (one non-test `use`
  today, and `calling → parameter_estimation` already exists in 13 files while the reverse is
  zero, so a top-level peer would acquire an edge into each); the directory form; and a
  `mod.rs`/`shape.rs` split, which it built and found clippy-clean but against
  `arch/module_layout.md` principle 1 and against the tree, where 11 of ng's 23 `mod.rs` files
  exceed 960 lines.

## 8. Missing tests to add now

Ten, from `reliability` and `idiomatic`. Six are in the fixes report as added; four are recorded
there as deferred with their owner named — full-precision float round-tripping (step C3's own
question, though the fixture now carries one such value), non-finite values, counts past
`i64::MAX`, and the property-based round trip (step C).

## 9. What's good

- **The weight-on-the-variant shape for a blend** (`BlendedSource::Blend { curve_weight }`) was
  called out by two agents as the right instinct; the fix for M9 is that same device extended to
  the two fields that sat outside it.
- **Every `pub` item carries a doc comment that says *why*** — `smells` named this the strongest
  quality in the change.
- **The module has no `use crate::` line at all**, which is what makes its placement cheap and its
  types independent of the calling loop's.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt -- --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file`
- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file -- --ignored regenerate` —
  rewrites the golden file after an intended change to the shape, and is the only thing that may.
