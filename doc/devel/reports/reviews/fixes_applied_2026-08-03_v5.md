# Fix Application Report: ng_read_filtering_stages_b2_2026-08-03.md

**Date:** 2026-08-03
**Source review:** `doc/devel/reports/reviews/ng_read_filtering_stages_b2_2026-08-03.md`
**Source state reviewed against:** branch `ng-generic-perf`, base `d5fb526`, B2 working-tree diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
Blockers: 0 · Majors: 2 · Minors: 8 (4 on B2, 4 from the milestone close-out) · Nits: 1

### Outcome totals
Applied: 8 · Applied with adaptation: 0 · Already fixed: 0 · **Deferred: 2** · Disputed: 0
Failed validation: 0 · Blocked by context mismatch: 0 · Superseded: 0 · Awaiting user answer: 0
(The Nit is folded into the deferred pair.)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **2,842 passed / 0 failed / 5 ignored**
- `cargo test --lib ng::` → 0, **1,543 passed / 0 failed / 2 ignored**
- `cargo test --examples` → 0, 52 passed / 0 failed
- `cargo doc --no-deps --lib` → 12 unresolved links, all pre-existing
- Performance check → not applicable; the only executable change is one field assignment
  relocated and a destructure

**Four dumps byte-identical** to the `8cf6f03` baseline by `cmp`; walk probe anchor exact at
`seconds=1.848`.

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | no multi-container multi-read-group CRAM fixture | Apply | Applied | `test_fixtures.rs`, `open_bam.rs` | Pass + mutation |
| M2 | Major | read-group *value* rested on one test | Apply | Applied | `mod.rs` | Pass + mutation |
| Mi1 | Minor | "both halves" reached fields by name | Apply | Applied | `container.rs` | Pass |
| Mi2 | Minor | `container.rs` has no test module | Defer | Deferred | None | N/A |
| Mi3 | Minor | aux-tag clear invariant untested | Defer | Deferred | None | N/A |
| Mi4 | Minor | impl report §4 over-generalised | Apply | Applied | the B2 report | N/A |
| Mi5 | Minor | **B1's deviation recorded inaccurately** | Apply | Applied | 3 documents | N/A |
| Mi6 | Minor | spec says "no code yet"; stale §1 bullet | Apply | Applied | the spec | N/A |
| Mi7 | Minor | §11 reuse map citations drifted | Apply | Applied | the spec | Links checked |
| Mi8 | Minor | "proves strictly more" still shipped | Apply | Applied | `filtering.rs` | Pass |

## 3. Questions asked and answers

None. Mi5 changes what the owner is being asked at Checkpoint B, but it is a correction of the
record rather than a new question.

## 4. Per-finding log

### M1 — the container-boundary gap
- **Final status:** Applied
- **Reasoning:** B2's whole subject is the read-group stamp, and this was the one shape of stamp
  defect the suite could not see. The reviewer proved it was not an equivalent mutant by writing
  the missing test and running it both ways.
- **Implementation:** `multi_container_cram_two_read_groups` in `test_fixtures.rs` — two `@RG`s
  of one sample, alternating record by record, so the group varies *within* a container and
  *across* the boundary — and
  `a_read_past_the_first_container_carries_its_own_read_group` in `open_bam.rs`, opened with
  `ReadGroupResolution::PerRecord` rather than `Sole` (the point is to reach the arm `Sole`
  short-circuits). Each fixture doc records why the gap existed.
- **Verification:** re-ran the reviewer's mutation myself — a hard-coded `ReadGroupId(0)` from
  the second container onwards, `grep -c`-confirmed present: **`1542 passed; 1 failed`, the new
  test alone**. Reverted; marker count 0.
- **Residual risk:** None.

### M2 — the value rested on one test
- **Final status:** Applied
- **Implementation:** `a_shared_cram_serves_each_open_only_its_own_reads` now asserts
  `collect_reads(&mut cursor)` against `vec![("a", 0), ("c", 0)]` instead of mapping the group
  away. The comment says why: this is the arm that decides the group itself, so a
  wrong-but-valid group is exactly the defect it should see.
- **Residual risk:** None.

### Mi1 — "both halves" by name
- **Final status:** Applied · exhaustive destructure of `NoodlesRawAlignedRead`, so a third
  field is a compile error here rather than a silently broken doc promise.

### Mi2, Mi3 — `container.rs` has no test module
- **Final status:** Deferred, recorded in the B2 impl report's §6.
- **Reasoning:** both are **pre-existing coverage of pre-existing behaviour** — the unnamed-record
  arm, empty spans, the clear-and-refill claim, and the aux-tag clear. B2 changed the function's
  signature, not those paths. A new test module for `container.rs` is its own piece of work, and
  bundling it into a two-file signature change would make B2's diff stop being about B2.
  Mi3 is additionally **latent, not live**: every caller passes a fresh buffer today.
- **Follow-up:** worth taking **before C2**, which moves the filtering loop and is exactly the
  change that could hand this function a buffer with a history.

### Mi4 — §4 over-generalised
- **Final status:** Applied · §4 now carries the seven-row mutation log and states plainly that
  the original conclusion was drawn from the one mutation that could not survive.

### Mi5 — B1's deviation was recorded inaccurately
- **Final status:** Applied
- **Reasoning:** **this is the finding that matters most outside the code.** Three documents said
  B1's new error variant went "against spec §1's *adds no new error*". I verified: spec §1 says
  "change the meaning of any error" — adding a variant does not — and the only "No new error
  type" sentence is arch §4, scoped to `ReadFilterError`, a different enum. The design authority
  is **silent** on adding a variant to `AlignmentFileError`, not against it. Since the owner is
  ruling on exactly this at Checkpoint B, the difference between "a rule was broken" and "a gap
  was filled" is the difference between two decisions.
- **Files changed:** `PROJECT_STATUS.md`, the B1 impl report, the B1 review and fix reports.

### Mi6, Mi7, Mi8 — design-doc and shipped-doc staleness
- **Final status:** Applied
- **Mi6:** the spec's header said "no code yet" after five code commits; it now says which
  milestones are built and warns that its present tense describes what C and D must still reach.
  Its §1 bullet described the probe loop B1 deleted, citing a line that is now an unrelated
  function — now annotated with what B1 actually did.
- **Mi7:** §11's reuse map is **the table Milestone C executes against**, and its citations had
  drifted by up to 120 lines. Refreshed, with a dated note saying so. Every `src/` link in the
  file was then checked to resolve.
- **Mi8:** `with_validated_contigs`' shipped doc still claimed the comparison "proves strictly
  more" — the claim B1's own review overturned, and the one that produced B1's Blocker. It now
  states that **two** checks replaced the loop and why the second is not redundant.

## 5. Deferred findings to carry forward
- **Mi2 / Mi3** — a test module for `container.rs`, covering the unnamed-record arm, empty
  spans, the buffer-shrink claim, and the aux-tag clear. **Before C2.**

  > **Closed at C1b (2026-08-03)**, and one word of this line was wrong. §4 above lists the third
  > item as *"the clear-and-refill claim"*; this line calls it *"the buffer-shrink claim"*. They
  > are different properties, and C1b's review caught the divergence — so C1b covers **both**:
  > `a_shorter_record_keeps_no_tail_of_the_longer_one_before_it` and
  > `shrinking_gives_back_the_slack_the_buffers_grew_by`. Under either reading the deferral is now
  > fully discharged.
  >
  > **The deferral's stated reason was also wrong**, and that matters more than the wording.
  > It said the findings were "latent, not live: every caller passes a fresh buffer today". They
  > are not latent — `ReadFilter` refills **one** buffer for a whole pass, so every read after the
  > first arrives with a history, and has since B2. Measured at C1b by instrumenting
  > `fill_raw_read` and running the CRAM cursor walk. A regression in the clears would corrupt
  > production reads today; nothing was waiting on C2.

## 6. Disputed findings
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no. One field assignment relocated, one destructure, plus tests and docs.
- **Outcome:** skipped. The walk probe was run anyway: `seconds=1.848`, in line with B1's
  post-change mean of 1.834 and well inside the session spread.

## 10. Commands run
`cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --lib`,
`cargo test --lib ng::`, `cargo test --examples`, `cargo doc --no-deps --lib`,
`cargo build --release --examples`, the four acceptance dumps and the walk probe against the
`8cf6f03` baseline, and two mutations `grep -c`-confirmed present before their runs and absent
after the reverts.

## 11. Command results
- fmt clean; clippy clean; `--lib` 2,842 / 0; `ng::` 1,543 / 0; examples 52 / 0
- four dumps **byte-identical**; walk probe anchor exact
- stale-group mutation → **killed by the new test alone** (`1542 passed; 1 failed`)
- read-group-stamp-deleted mutation → killed by two tests
- final tree: no `MUTATION` markers remain

## 12. Notes
- **A perl one-liner corrupted the spec mid-run and was caught immediately.** I used `|` as the
  substitution delimiter on patterns containing markdown table pipes; perl parsed the first
  table pipe as the delimiter and prepended four mangled rows to line 1. Restored from git and
  redone with the `Edit` tool, one row-block at a time. Recorded because the failure is silent
  in the sense that matters — it produced a file that still looked plausible at a glance.
- **The review caught an error in its own brief.** I told the reviewer that B1 deviated from spec
  §1; it checked and found no such sentence (Mi5). Worth noting that the fan-out is only worth
  its cost if the agents are willing to contradict the orchestrator.
- Two review agents ran in isolated worktrees; per-category findings are kept in the gitignored
  `tmp/review_2026-08-03_ng-read-filtering-stages-b2/`.
