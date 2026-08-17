# Code Review: ng cohort merge — B1 (projection onto the locus span)
**Date:** 2026-08-17
**Reviewer:** rust-code-review skill (orchestrator, five category sub-agents)
**Scope:** the working-tree diff of plan step B1, captured as `edccb48a` on branch `ng-cohort-merge`
**Status:** Approve-with-changes

---

### 1. Scope

- **What was reviewed:** one step's diff — a new file plus two one-line changes.
- **Reviewed against:** `edccb48a4c271f0a0e39a22bdea2b1862d097ecc` (a `git stash create` object of the
  uncommitted step), branch `ng-cohort-merge`.
- **In-scope files:**
  - [build.rs](../../../../src/ng/run/cohort_merge/build.rs) — new: `LocusReference`, `project_into`,
    `offset_within`, `NOT_COVERED`, 12 tests;
  - [close.rs](../../../../src/ng/run/cohort_merge/close.rs) — the `span_of` visibility change and its doc;
  - [mod.rs](../../../../src/ng/run/cohort_merge/mod.rs) — `pub mod build;` and one module-doc line;
  - [ng_cohort_merge_b1_2026-08-17.md](../implementations/ng_cohort_merge_b1_2026-08-17.md) — its claims
    about the code.
- **Deliberately out of scope:** the rest of `close.rs` (committed at A4, reviewed then); the psp path;
  milestones C, D and E; the pre-existing `--all-targets` clippy failures.
- **Categories dispatched:** `reliability` (always; the step's failures are silent-wrong-answer shaped),
  `errors` (always; the step chose assertions over errors and that choice is the step),
  `naming` (always; the module's vocabulary is fixed by spec §1.3),
  `idiomatic` (always; the projection runs once per observation per locus),
  `smells` (always, carrying the claim-verification pass of skill step 8a).

**One process failure, recorded because it changes how to read the reviews.** The five agents were
dispatched **without** `isolation: "worktree"`, so they started in the main checkout rather than in
worktrees of their own. All five detected it, none detached the main checkout's `HEAD`, and each verified
its in-scope files byte-identical to `edccb48a` before reviewing; three then made their own worktrees for
mutation work. The main checkout was left as it was. Two review worktrees remain to be pruned.

### 2. Verdict

**Approve-with-changes.** No finding says the code computes a wrong projection — a reviewer re-derived
every asserted byte string by hand and all four are right. Every Major is about what the tests and the
prose fail to hold down.

### 3. Execution status

Run by the orchestrator in the container, in the main checkout:

| command | result |
|---|---|
| `cargo fmt --check` | clean, exit 0 |
| `cargo clippy --lib --all-features -- -D warnings` | clean, exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | **red, pre-existing** — 49 errors in `examples/`, `benches/` and other modules' test code, none in a file this step touches. Standing item under the block's `Open:` |
| `cargo test --lib ng::run::cohort_merge` | `ok. 50 passed; 0 failed` |
| `cargo test --lib` | `ok. 3671 passed; 0 failed; 11 ignored` (575.81 s) |

Two agents re-ran the lib clippy and the module tests in their own worktrees and got the same results.
Findings labelled "Needs verification": **0**.

### 4. Open questions and assumptions

1. **Does `LocusReference::over` accept a locus of any verdict?** The doc bounds its cost by
   `max_cohort_locus_span`, which bounds *assembly*, not closing — a `Failed` locus is uncapped.
   Raised independently by `errors` and `smells`. Affects M4.
2. **When observations arrive from a psp file, does the reference-width assertion become a `RunError`?**
   Today the only producer is this crate's generator, so a panic is right; the class changes with the
   producer. Affects M5.
3. **Is a partial observation refused, or skipped?** B1 panicked on one; `errors` argues partials are
   routine data rather than a caller's mistake. Affects M2.

### 5. Top 3 priorities

1. **M1 — a sequence is not tied to the member it came from** (`errors`, `idiomatic`, `smells`, `reliability`
   all surfaced it): mispairing two arguments of one call yields a well-formed allele and no panic.
2. **M3 — the ceiling test pins half of what it claims**, and the report says it pins both: replacing
   `span_of` with `GenomeRegion::len()` inside the projection left all 50 tests green.
3. **M6 — no fixture gives a sample two observations in one locus**, so the gather's inner loop never
   iterates twice, and that is the case spec §4.2 and plan step B3 are written around.

### 6. Findings

#### Major

**M1: build.rs:155 — nothing ties `sequence` to `member` in `project_into`, and mispairing them is silent.**
**Categories:** idiomatic (filed), errors, smells, reliability (cross-category). Confidence: High.
The two borrowed arguments are meaningful only as a pair — the member supplies the offset and the width of
reference replaced, the sequence supplies the bases — but any member of the locus may be passed with any
other member's sequence. The three assertions do not catch it: a member from the *same* locus passes contig
and reach, so the result is a sequence padded at the wrong offset — a well-formed byte string B2 would
accept as an allele. The same shape re-does two assertions, a `checked_sub` and a `usize::try_from` for
every sequence a member carries. Fix: a placement handle returned by `placing(member)`, with `project_into`
hanging off it.

**M2: build.rs:161 — the `Partial` panic fires on ordinary data, not on a caller's mistake.**
Category: errors. Confidence: High.
Every other release-level check in the file guards a state the walk cannot produce. `ReadWitness::Partial`
is a routine product of real reads — the witness type's own doc records 6,704 partials on one tomato
chromosome — and it arrives in `observations`, a public `Vec`. The instruction to reach sequences through
`complete_observations()` is a convention, not a type. `locus_generation/mod.rs:78-96` records this project
rejecting a debug-only guard in the same situation and choosing a *total* derivation rather than promoting
the guard.

**M3: build.rs:430-451 — the ceiling test pins the gather against `GenomeRegion::len()` but not the
projection, though its doc and the report both say it pins both.**
Category: reliability. Confidence: High.
The fixture projects the SNP at `u64::MAX − 2`, where `span_of` and `len()` both answer 1. The member that
separates them — the one ending *at* the ceiling — is gathered and never projected. Mutation M3 (replace
`span_of` with `len()` inside `project_into`) survived all 50 tests. In the release profile `covered` would
be 0, the suffix would become the whole reference, and that member's allele would come back `AACGT`
instead of `A`.

**M4: build.rs:30-37 — nothing pins the one property `NOT_COVERED` must have.**
Category: reliability. Confidence: High.
Mutation M7 (`0` → `b'N'`) survived. Unlike the coordinate ceiling this is reachable on real data: ng's
fetch folds every non-ACGT byte to `N`, so any locus over an assembly gap carries them, and a sentinel
spelled `b'N'` would refuse a locus its members do cover.

**M5: build.rs:169-173 — `project_into`'s contig guard has no test, and removing it produces a wrong allele
with no panic.**
Category: reliability. Confidence: High.
All three contig/containment tests go through `over`; the projection carries its own copies and only the
witness one is tested. Mutation M12 was the only survivor whose mutated run produced **no panic at all**: a
member at `contig 1:12-12` projected onto a `contig 0:10-14` locus came back padded from contig 0's
reference.

**M6: build.rs:77-78 — no fixture gives a sample more than one observation, so the gather's inner loop is
never iterated twice.**
Category: reliability. Confidence: High.
Mutation M15 (`take(1)`) survived. Spec §4.2's "two of its own observations" and plan B3 are written around
exactly this case; the existing spread test spreads across *samples*, which is the other loop.

**M7: build.rs:84 — the reference-width check is about the observation's own well-formedness, and the doc's
justification does not cover it.**
Category: errors. Confidence: Medium.
`over`'s doc says every check is against a caller's mistake. That holds for contig, reach and coverage —
properties of how the walk paired members with a region — but not for this one: the walk never reads
`reference_bases`, so this asserts a *producer's* guarantee. Today the only producer is in-process; once
observations are decoded from a psp file it is corrupt input, which is arch §5's `RunError` class.

#### Minor

**Mi1 (errors):** the uncovered-locus message names the locus but not the gap position — the catch-all
assertion with the least to say.
**Mi2 (errors):** the messages identify the sample by a bare `usize` while reading as a name; spec §12 and
arch §5 both want the sample named.
**Mi3 (errors, smells):** `over` allocates the whole span before any check, and enforces neither the width
bound nor the `Verdict::Build` its doc rests on; an absurd region aborts with no message.
**Mi4 (errors):** three `expect`s without `// PANIC-FREE:`, and one states an invariant nothing enforces
("a locus span is at most the length of a contig") — the opposite of what `offset_within`'s doc argues
about the same type twenty lines below.
**Mi5 (errors, idiomatic, smells):** overlapping members that disagree on a reference base are resolved
last-write-wins, silently — the same plausible-allele failure the four assertions exist to catch.
**Mi6 (reliability):** `project_into`'s reach guard is untested; removed, its named message becomes a raw
out-of-range slice panic.
**Mi7 (reliability):** `bases()`'s "this is also the reference allele" is a round-trip law with no test —
every projection in the file uses a non-reference sequence, so the identity is never exercised.
**Mi8 (reliability):** `the_reference_is_gathered_across_members…`'s prose says "neither member alone
carries more than three of those bases" while line 385 of the same test asserts one carries four.
**Mi9 (reliability):** `a_member_covering_the_whole_locus_projects_to_its_own_bases` is satisfied by an
implementation that ignores the reference entirely — both padding slices are empty, and the gathered
reference is built and discarded.
**Mi10 (reliability):** `close.rs`'s revised `span_of` doc still recommends `GenomeRegion::len()` for "an
observation's own span" — the written invitation to make mutation M3, which the suite did not catch.
**Mi11 (naming, errors, idiomatic, smells — convergent):** `member` names two types in one file: a
`SampleMembers` in the gather's outer loop, a `SampleLocusObservations` in `project_into`'s signature. One
`assert_eq!` uses both senses in the same call, and spec §4.2 sides with the messages.
**Mi12 (naming):** `LocusReference` reads as a Rust borrow beside `&ClosedLocus`, and is a third crate
spelling for reference bases over a stretch (`*RefSeq` on types, `reference_bases` on fields).

#### Nits

`bases` stays `mut` after the gather; a test import already supplied by `use super::*`; deref coercion in
the fixture's struct-field position; an avoidable clone in the partial test; `Clone`/`PartialEq`/`Eq`
derived with no consumer; `assert_eq!(allele.len(), 8)` implied by the line above it; two assertions about
the fixture's own bytes that no implementation can fail; `sequence()`'s doc says fields are "zeroed" while
`num_obs` is 3; the report's "Open-coded, an inner region…" missing its conditional; a rename pass over
test locals (`whole`, `short`, `early`, `overhanging`, `elsewhere`) and helpers (`closed`, `projected`).

### 6a. Mutation testing

`reliability`, in its own worktree: **17 mutations run, 5 survived, 0 changed-no-behaviour.** Every survivor
was proven with a purpose-written test that passes on the real code and fails under the mutation — M3
(`span_of`→`len()` in the projection), M7 (sentinel `b'N'`), M12 (drop the projection's contig assert), M13
(drop its reach assert), M15 (gather only each sample's first observation). The twelve killed mutants
included both slice swaps, the dropped `clear()`, `offset + sequence.bases.len()`, all four gather
assertions, and `checked_sub` → `saturating_sub`.

### 7. Out of scope observations

- `SampleLocusObservations::locus_len()`, `ClosedLocus::span()` and `span_of` are three spellings of one
  quantity, two of them outside this step (naming, cross-category).
- The module's blanket `pub` is `mod.rs`'s recorded decision, to be narrowed when the caller objects land
  (idiomatic, not filed).

### 8. Missing tests to add now

From `reliability`, beyond the fixes above: `where_members_overlap_the_last_reference_written_wins` (pins
which rule is in force, since "skip what is written" is the opposite rule) and
`project_into_reuses_a_buffer_whose_capacity_already_exceeds_the_projection` (the scratch buffer at its
widest swing).

### 9. What's good

- The `NOT_COVERED` sentinel's justification checks out against `RefSeq::fetch_into`'s canonicalisation,
  and a silent-but-total `N`-fill would have been the wrong instrument (errors).
- Every asserted byte string in the tests re-derived correctly by hand, and the gapless-coverage claim on
  `over` proved by induction over `close.rs`'s walk (smells).
- `offset_within`'s `checked_sub` is the right instrument and its doc gives an accurate reason (errors).
- The three `should_panic` tests on the gather each reach their own assertion — the preceding checks all
  pass on those fixtures, so control genuinely arrives where the name says (reliability).

### 10. Commands to re-verify

`./scripts/dev.sh cargo fmt --check`; `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings`;
`./scripts/dev.sh cargo test --lib`. The mutation driver and log are in the review worktree
`/Users/jose/devel/pop_var_caller-review-b1`; the per-category findings are in
`tmp/review_2026-08-17_cohort_merge_b1/`.
