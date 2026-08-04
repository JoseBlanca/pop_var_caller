# Fix Application Report: ng_read_filtering_stages_a3_2026-08-03.md

**Date:** 2026-08-03
**Source review:** `doc/devel/reports/reviews/ng_read_filtering_stages_a3_2026-08-03.md`
**Source state reviewed against:** branch `ng-generic-perf`, base `db3057a` (A2), A3 working-tree diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 5 · Nits: 1 group (2 items)

### Outcome totals
- Applied: 6 · Applied with adaptation: 0 · Already fixed: 0 · **Deferred: 2**
- Disputed: 0 · Failed validation: 0 · Blocked by context mismatch: 0 · Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean (5.69s)
- `cargo test --lib` → 0, **2,839 passed / 0 failed / 5 ignored** — unchanged
- `cargo test --lib ng::` → 0, **1,540 passed / 0 failed / 2 ignored** — unchanged
- `cargo test --examples` → 0, 52 passed / 0 failed
- `cargo test --all-targets --all-features` → not run; pre-existing panic in
  `benches/psp_writer_perf.rs:386`
- `cargo doc --no-deps` → 1, 12 unresolved links, all pre-existing, none in a touched file
- `cargo audit` → not run (not part of this project's gate)
- Performance check → **not applicable**; the only executable change is a four-line wrapper
  deletion whose body moved to its callee, and it is not on a bench path

**Four dumps byte-identical** to the `8cf6f03` baseline by `cmp` after the fixes; walk probe
`loci=236081 observations=251786 reads_admitted=54709` at `seconds=1.887`.
**The suite count did not move at A3**, in either direction.

### Unresolved high-priority findings
None. Both `Deferred` items are naming objections against **architecture-prescribed names**, and
both are Checkpoint A questions for the owner.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| Mi1 | Minor | old name spelled out in the type's summary | Apply | Applied | No | `region_raw_aligned_reads.rs` | Pass | No |
| Mi2 | Minor | `fill_raw_read` a 1:1 wrapper the rename made incoherent | Apply | Applied | No | `container.rs` | Pass + mutation | No |
| Mi3 | Minor | three dangling links from the `git mv` | Apply | Applied | No | spec + arch | Pass | No |
| Mi4 | Minor | `PROJECT_STATUS.md` contradicted itself on the suite count | Apply | Applied | No | `PROJECT_STATUS.md` | N/A | No |
| Mi5 | Minor | three impl-report counts + two undisclosed edits | Apply | Applied | No | the impl report | N/A | No |
| Q1 | — | `RawReadIndex` disputed (arch §2 prescribes it) | Defer | Deferred | No | None | N/A | Yes — Checkpoint A |
| Q2 | — | `fill_raw_read`'s signature disputed (arch §2) | Defer | Deferred | No | None | N/A | Yes — Checkpoint A |
| Nits | Nit | 2 items | Apply | Applied | No | 2 files | Pass | No |

## 3. Questions asked and answers

None asked mid-run. Q1 and Q2 are **deferred to Checkpoint A**, where the plan already pauses —
see §5.

## 4. Per-finding log

### Mi1 — the renamed type's summary still spelled the old name out
- **Severity:** Minor · **Final status:** Applied
- **Reasoning:** *"The records of one region of one file"* is `RegionRecords` written out and
  reordered, so **neither** grep could see it — not the exact-identifier one, not the widened
  spelled-out one A2's review made this step inherit. And it is the first sentence on the type
  the step is named after.
- **Implementation summary:** → "The raw aligned reads of one region of one file".
- **Files changed:** `src/ng/read/input/region_raw_aligned_reads.rs` · **Residual risk:** None.

### Mi2 — `fill_raw_read` was a 1:1 wrapper whose ends the rename made disagree
- **Severity:** Minor · **Final status:** Applied
- **Reasoning:** **this reverses a deferral the impl report had already written down.** The
  agent verified `fill` had exactly one caller (the wrapper) and `fill_raw_read` exactly one
  (`cram.rs:188`), and argued that collapsing a name-only 1:1 indirection *is* a rename plus a
  four-line deletion — so deferring it would push a naming edit into a later behaviour diff,
  which is the one thing this milestone exists to prevent. That is right.
- **Implementation summary:** the shim deleted; the private `fill` promoted to
  `pub(crate) fn fill_raw_read`, keeping its doc and gaining a note that `out` is only the
  **record half** of the caller's `NoodlesRawAlignedRead` — the read-group half being this arm's
  to set at the call site, because on CRAM the group is decided at decode.
- **Verification performed:** the agent's own mutation re-run **after** the collapse:
  `container.fill_raw_read(i, …)` → `(0, …)`, marker `grep`-confirmed present, →
  `FAILED. 1536 passed; 4 failed`. Reverted, marker count back to 0, suite green.
- **Files changed:** `src/ng/read/input/aligned_reads_reader/container.rs` · **Residual risk:**
  None. Q2 below is a *separate* objection to the same function's signature, and is deferred.

### Mi3 — three dangling links created by the `git mv`
- **Severity:** Minor · **Final status:** Applied (convergent, two agents)
- **Reasoning:** repairing a citation path is not editing a design, and A2 set the precedent for
  its own move. The agent also caught the real risk: `PROJECT_STATUS.md`'s deferred-sweep note
  **predated A3** and named only A1's and A2's identifiers, so these three would have fallen
  between the deferral and the fix.
- **Implementation summary:** `spec/read_filtering_stages.md:154` (evidence pointer, stale `#L222`
  anchor dropped) and `:411`, and `arch/read_filtering_stages.md:230`, now resolve to
  `region_raw_aligned_reads.rs`; the two reconciliation rows keep the old type name in their
  "existing code" column, which is what that column describes, and record that the rename landed
  at A3. The rename *tables* naming `region_records.rs` as the "before" side were left alone.
- **Files changed:** the spec and the arch · **Residual risk:** None.

### Mi4 — `PROJECT_STATUS.md` contradicted itself on the suite count
- **Severity:** Minor · **Final status:** Applied
- **Reasoning:** my own error, from the block written at A1 and not updated when A1's fixes
  landed. One line said "suite 2,837 / `ng::` 1,538, both unchanged"; another in the same block
  correctly recorded the +2.
- **Implementation summary:** the line now states 2,839 / 1,540, says the dumps are compared with
  `cmp` rather than by line count, and attributes the +2 to A1's review explicitly.
- **Files changed:** `PROJECT_STATUS.md` · **Residual risk:** None.

### Mi5 — three impl-report counts, and two undisclosed edits
- **Severity:** Minor · **Final status:** Applied (three-way convergent)
- **Implementation summary:** "nine ng files" → **eight**; "six comment lines" → **seven blocks,
  sixteen lines**; the module-path/method-name conflation corrected. Both undisclosed edits now
  disclosed (rustfmt's `use` reflow in `cursor.rs`, and the `held` field's re-wrap, which no
  rename explains). The report's *understatement* was also corrected — `region_records` appears
  twice more in `benches/psp_reader_perf.rs` than it claimed.
- **Files changed:** the impl report · **Residual risk:** None.

### Q1 — `RawReadIndex` (deferred)
- **Final status:** Deferred · **Reasoning:** **arch §2 prescribes this name.** The `naming`
  agent's objection is substantive — "Index" names the container's *field*, not one entry, and
  the surrounding code already calls an entry `entry`; and `RawRead` abbreviates
  `RawAlignedRead`, which appears in full 85 times in `src/ng` against `RawRead` at exactly the
  six sites A3 created. But changing it means amending the architecture's rename table, which
  this skill does not do. Landed as specified. Proposed alternative: `RawAlignedReadEntry`.
- **Follow-up:** Checkpoint A.

### Q2 — `fill_raw_read`'s signature (deferred)
- **Final status:** Deferred · **Reasoning:** the name is arch §2's, and the agent's fix is to
  change the *signature* to `&mut NoodlesRawAlignedRead` so the function sets both halves. That
  is a behaviour-shaped edit, not a rename, and belongs to Milestone C at the earliest — the
  CRAM arm's read-group stamping is exactly what the contract list calls out as its documented
  exception. The wrapper collapse (Mi2) is the part that *was* in scope, and it was done.
- **Follow-up:** Checkpoint A.

### Nits
- **Applied:** `open_bam.rs:1440`, still 103 **characters** after a re-wrap whose purpose was to
  get it under 100 (the agent's character-accurate audit matters: these lines carry multibyte
  `—`, `→` and `§`); and the layer diagram's gloss, which a first attempt had cut to "this
  region's only" — a dangling possessive — now spelled in full by widening the column instead,
  matching spec §6's own diagram.

## 5. Deferred findings to carry forward
- **Q1** — `RawReadIndex` → `RawAlignedReadEntry`? Needs arch §2 amended. **Checkpoint A.**
- **Q2** — should `fill_raw_read` take `&mut NoodlesRawAlignedRead` and set both halves?
  **Checkpoint A**, and probably Milestone C if taken.

## 6. Disputed findings to return to reviewer
None. Q1 and Q2 are deferred rather than disputed: the reviewer's reasoning is sound, and the
obstacle is authority over the architecture, not correctness.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no. One four-line wrapper deleted with its body already in the callee; every
  other change is prose or a doc-link path.
- **Outcome:** skipped. The walk probe was run anyway: `seconds=1.887` against `1.922` pre-fix
  and `1.846` at the base — run-to-run noise across the milestone.

## 10. Commands run
- `cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib`, `cargo test --lib ng::`, `cargo test --examples`
- `cargo build --release --examples`
- one mutation (`fill_raw_read(i, …)` → `(0, …)`), `grep -c`-confirmed present before its run
  and confirmed absent after the revert
- the four acceptance dumps + `ng_generic_walk_probe`, `cmp`'d against the `8cf6f03` baseline

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 2,839 passed / 0 failed / 5 ignored
- `cargo test --lib ng::` → 0, 1,540 passed / 0 failed / 2 ignored
- mutation `MUT-A3` → `FAILED. 1536 passed; 4 failed`; after revert, `grep -c MUT-A3` → 0 and
  `1540 passed; 0 failed`
- four dumps → **byte-identical**; walk probe → anchor exact

## 12. Notes
- **Mi2 reverses a deferral this step's own impl report had written down**, on the reviewer's
  argument. Recorded in both documents rather than quietly changed.
- **Q1 and Q2 are the first findings in this plan that the architecture, not the code, is the
  obstacle to.** They are raised at the checkpoint rather than absorbed, because amending arch
  §2's rename table is the owner's call.
- Three review agents ran in isolated worktrees; per-category findings are kept as an audit
  trail in the gitignored `tmp/review_2026-08-03_ng-read-filtering-stages-a3/`.
