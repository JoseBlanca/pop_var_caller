# Fix Application Report: ng_read_filtering_stages_a2_2026-08-03.md

**Date:** 2026-08-03
**Source review:** `doc/devel/reports/reviews/ng_read_filtering_stages_a2_2026-08-03.md`
**Source state reviewed against:** branch `ng-generic-perf`, base `5438927` (A1), A2 working-tree diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 1
- Minors: 7
- Nits: 1 group (6 items)

### Outcome totals
- Applied: 6
- Applied with adaptation: 1
- Already fixed: 0
- Deferred: 2
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean (5.53s)
- `cargo test --lib` → 0, **2,839 passed / 0 failed / 5 ignored** — unchanged from A1
- `cargo test --lib ng::` → 0, **1,540 passed / 0 failed / 2 ignored** — unchanged from A1
- `cargo test --examples` → 0, 52 passed / 0 failed
- `cargo test --all-targets --all-features` → not run; pre-existing panic in
  `benches/psp_writer_perf.rs:386`
- `cargo doc --no-deps` → 1, 12 unresolved links, all pre-existing, none in a touched file
- `cargo audit` → not run (not part of this project's gate)
- Performance check → **not applicable**; every `Apply` was a doc comment or a doc-link path

**The four acceptance dumps remain byte-identical** to the `8cf6f03` baseline by `cmp` after
the fixes — 251,792 / 4,406 / 1,718,914 / 11,945 lines — and the walk probe prints
`loci=236081 observations=251786 reads_admitted=54709` at `seconds=1.880`.

**The suite count did not move at A2**, in either direction. Every fix below is prose.

### Unresolved high-priority findings
None. Both `Deferred` items are Minor, both are questions for Checkpoint A, and both are
recorded in `PROJECT_STATUS.md`.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | contract list holds a universal the CRAM arm breaks | Apply | Applied | No | `aligned_reads_reader/mod.rs` | Pass | No |
| Mi1 | Minor | module doc's inventory is false | Apply | Applied | No | `aligned_reads_reader/mod.rs` | Pass | No |
| Mi2 | Minor | ten "record reader" sites the grep cannot see | Apply | Applied | No | 5 files | Pass | No |
| Mi3 | Minor | the naming justification contradicts itself | Apply w/ adaptation | Applied with adaptation | No | `aligned_reads_reader/mod.rs` | Pass | No |
| Mi4 | Minor | impl report arithmetic wrong in two places | Apply | Applied | No | the impl report | N/A | No |
| Mi5 | Minor | four `pub(crate) mod` wider than needed | Defer | Deferred | No | None | N/A | Yes — Checkpoint A |
| Mi6 | Minor | five live design docs name `record_reader/` | Defer | Deferred | No | None | N/A | Yes — Checkpoint A |
| Mi7 | Minor | two links in this plan's own docs broken by the `git mv` | Apply | Applied | No | spec + arch | Pass | No |
| Nits | Nit | 6 items | Apply (4) / Defer (2) | Applied | No | 3 files | Pass | Partly |

## 3. Questions asked and answers

None. Mi5 and Mi6 are both owner decisions and were **deferred to Checkpoint A**, where the
plan already pauses, rather than asked as separate blocking questions.

## 4. Per-finding log

### M1 — the contract list states as a universal something the CRAM arm does not do
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** High confidence, verified by the agent against `cram.rs:189`. The bullet is
  pre-existing, but A2's new "what every arm yields" paragraph is what makes the module doc
  read as a complete statement — so this step is where the incompleteness becomes misleading.
- **Implementation summary:** the bullet now states the clearing rule, then names the CRAM arm
  as **the one exception in this contract** and why it exists (a CRAM stores the read group as a
  container-level number, so it is decided at decode and travels with the record). A closing
  sentence says why the exception is written *here*: this list is where a new arm's author
  checks their arm, and a universal with a live counterexample is how an arm ends up clearing a
  field it should have carried.
- **Review suggestion used verbatim?:** Substantially — the agent supplied a diff; the closing
  rationale sentence is added.
- **Verification performed:** re-read `cram.rs:189` and the arm's own doc to confirm the
  exception is real and described consistently in both places.
- **Files changed:** `src/ng/read/input/aligned_reads_reader/mod.rs` · **Validation:** clippy +
  `cargo doc` clean · **Residual risk:** None.

### Mi1 — the module doc's inventory is false
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** "What is here so far … The CRAM arm lands in Milestone E" becomes
  "What is here", listing all three arms plus `container` — which the old text did not mention
  at all — and saying explicitly that `container` holds the CRAM decode's unit of work and **is
  not an arm** (which also answers the `module_structure` agent's fourth Minor about
  `container`'s home being undocumented). The build-order rationale is kept, in the past tense.
- **Files changed:** `src/ng/read/input/aligned_reads_reader/mod.rs` · **Residual risk:** None.

### Mi2 — ten sites spell the old type name with a space
- **Severity:** Minor · **Final status:** Applied (three-way convergent)
- **Implementation summary:** all ten replaced across `open_bam.rs` (4), `cursor.rs` (2),
  `input/mod.rs` (2), `in_memory.rs` (1) and `bam.rs` (1). The last is the one that mattered
  most — an `unreachable!` message an operator would read, naming a type that no longer exists.
  `src/psp/kind.rs:77`'s "typed-record reader" is production's and unrelated; left alone.
- **Verification performed:** `grep -rniE "record[ -]readers?" src/ng` → no matches.
- **Files changed:** five · **Follow-up:** the widened grep is written into the impl report as
  the check **A3 must use**, since the narrow one would certify it the same incomplete way.

### Mi3 — the naming justification contradicts the module title and the impl report
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied with adaptation
- **Reasoning:** two agents converged, from opposite directions: `naming` said the enum doc's
  justification ("the name says reads because that is what the file holds") is the exact premise
  the impl report used to defend keeping "records"; `extras` said A2's reword of the module
  opener to "reads" is what broke the agreement with the contract list below.
- **Implementation summary:** **adapted** — the reviewers offered either a one-sentence
  re-justification or reverting the rewordings and doing a whole-module vocabulary pass at
  Checkpoint A. Neither alone is right: reverting leaves the module silent about a distinction
  it genuinely makes, and a blanket pass would erase a real one. So the module doc now **names
  the distinction**: a caller wants *reads*, a record is the *encoding* of one, the contract is
  about finding and unpacking those, and "where this module says record it means the thing on
  disk; where it says read it means what a caller gets back." The enum's own justification is
  re-pointed from "what the file holds" to "what a caller wants", which is both true and no
  longer in competition with the other word.
- **Files changed:** `src/ng/read/input/aligned_reads_reader/mod.rs`
- **Follow-up:** the impl report's §6 note was rewritten to match; it previously defended the
  choice on the grounds the review showed were self-defeating.

### Mi4 — the impl report's arithmetic
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** "four places" → "five"; "67 sites across eleven files" → "67 type
  names plus 23 module-path sites, 90 substitutions across eleven files", with the nine/two/zero
  file split spelled out. Also added the clause the `refactor_safety` agent asked for: two of
  the 90 are inside `debug_struct` string literals — the only non-comment, non-identifier bytes
  in the diff — and nothing `{:?}`-formats these readers, so the "no behaviour change" claim is
  airtight rather than nearly so.
- **Files changed:** the impl report · **Residual risk:** None.

### Mi5 — the four `pub(crate) mod` declarations
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** the finding is correct and the agent verified plain `mod` compiles and tests
  clean. But it is a visibility change, not a rename, and Milestone A's whole discipline is that
  the diff *is* the rename. It also pairs naturally with A1's Mi2 (the `pub` on
  `NoodlesRawAlignedRead`) — **one visibility decision at Checkpoint A, not two scattered ones.**
- **Files changed:** None · **Follow-up:** Checkpoint A.

### Mi6 — five live ng design docs still name `record_reader/`
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** real, and convergent across three agents. Deferred because **A1's review found
  the same class** (`README.md`, `spec/`+`arch/read_filtering.md`, `arch/alignment_file.md`,
  `arch/read_groups.md` still saying `RawRecord`), and A3 will add a third batch. One sweep at
  Checkpoint A, over `doc/devel/ng/{spec,arch}/*.md` and `README.md`, covering all three renames,
  is cheaper and less error-prone than three partial passes — and it can be verified in one go
  with the widened grep. Historical reports and the `PROJECT_STATUS.md` narrative stay frozen.
- **Files changed:** None · **Follow-up:** Checkpoint A.

### Mi7 — two links in this plan's own spec and arch, broken by the `git mv`
- **Severity:** Minor · **Final status:** Applied
- **Reasoning:** distinct from Mi6, and not deferred with it: these are *this step's* governing
  documents and *this step* broke them. Repairing a citation path is not editing a design.
- **Implementation summary:** `spec/read_filtering_stages.md:147`'s evidence pointer now resolves
  to `aligned_reads_reader/container.rs`; `arch/read_filtering_stages.md:231`'s reconciliation
  row keeps the old type name in its "existing code" column (which is what it describes) but
  points at the live path and records that the rename landed at A2. The three *rename tables*
  that name `record_reader/` as the "before" side were correctly left alone.
- **Files changed:** the spec and the arch · **Residual risk:** None.

### Nits — 6 items
- **Final status:** Applied (4 applied, 2 deferred)
- **Applied:** four "a `AlignedReadsReader`" → "an" (the substitution carried the old article
  through, while `mod.rs` had been hand-corrected — which is what told the reviewers the others
  were unexamined); the ASCII layer diagram in `region_records.rs` re-padded (the new name is
  six characters longer, so its description column had drifted out of line with the three below
  it); the one comment line past 100 columns re-wrapped; `region_records.rs`'s opening paragraph
  re-wrapped with it.
- **Deferred:** the two hand-written `Debug` impls that pick fields without destructuring
  `Self` (pre-existing, `finish_non_exhaustive()` declares the omission, and A2 changed only the
  label string); `cursor.rs:25`'s stale milestone sentence (same class as Mi1 but in a file
  whose own doc A3 will rewrite).

## 5. Deferred findings to carry forward
- **Mi5** — narrow the four `pub(crate) mod` inside `aligned_reads_reader/`. Owner decision at
  Checkpoint A, together with A1's Mi2.
- **Mi6** — the design-doc sweep over five live ng docs. Checkpoint A, once for all three
  renames.
- **Nits (2)** — the two manual `Debug` impls; `cursor.rs`'s stale milestone sentence.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — every `Apply` was a doc comment, a prose substitution or a doc-link path.
  No executable line changed.
- **Outcome:** skipped. The walk probe was run anyway, as this plan does at every step:
  `seconds=1.880`, against `1.871` pre-fix and `1.846` at the base — run-to-run noise.

## 10. Commands run
- `cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib`, `cargo test --lib ng::`, `cargo test --examples`
- `cargo build --release --examples`
- `grep -rniE "record[ -]readers?" src/ng` (the widened check, before and after)
- the four acceptance dumps + `ng_generic_walk_probe`, `cmp`'d against the `8cf6f03` baseline

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 2,839 passed / 0 failed / 5 ignored
- `cargo test --lib ng::` → 0, 1,540 passed / 0 failed / 2 ignored
- `cargo test --examples` → 0, 52 passed / 0 failed
- `grep -rniE "record[ -]readers?" src/ng` → 1 (no matches) after the fix
- four dumps → **byte-identical**; walk probe → anchor exact

## 12. Notes
- **The most useful thing this review did was not find a bug — it found the limit of the
  step's own verification.** `grep -rn "RecordReader\|record_reader" src` returning empty was
  the evidence A2 rested on, and it could not see the type's name written with a space. Ten
  sites, one of them a panic message. The widened check is now written into the impl report so
  A3 inherits it rather than repeating the gap.
- **Mi3's fix is an adaptation, not the reviewers' proposal**, because both proposals would have
  lost something real. Recorded above and in the impl report's §6.
- Four review agents ran in isolated worktrees; per-category findings are kept as an audit trail
  in the gitignored `tmp/review_2026-08-03_ng-read-filtering-stages-a2/`.
