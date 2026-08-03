# Fix Application Report: ng_read_filtering_stages_a1_2026-08-03.md

**Date:** 2026-08-03
**Source review:** `doc/devel/reports/reviews/ng_read_filtering_stages_a1_2026-08-03.md`
**Source state reviewed against:** branch `ng-generic-perf`, base `8cf6f03`, A1 working-tree diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 2
- Majors: 1
- Minors: 8
- Nits: 1 group (11 items)

### Outcome totals
- Applied: 10
- Applied with adaptation: 1
- Already fixed: 0
- Deferred: 4
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean (5.24s)
- `cargo test --lib` → 0, **2,839 passed / 0 failed / 5 ignored** (was 2,837)
- `cargo test --lib ng::` → 0, **1,540 passed / 0 failed / 2 ignored** (was 1,538)
- `cargo test --examples` → 0, 52 passed / 0 failed
- `cargo test --all-targets --all-features` → **not run**; aborts on a pre-existing panic in
  `benches/psp_writer_perf.rs:386`, unrelated to this work
- `cargo doc --no-deps` → 1, 12 unresolved links, all pre-existing, none in a touched file
- `cargo audit` → not run (not part of this project's gate)
- Performance check → **not applicable**; no `Apply` touched a hot path covered by `benches/`
  (the changes are doc comments, test code, and one import).

**The evidence that matters for a renames milestone:** all four acceptance dumps remain
**byte-identical** to the `8cf6f03` baseline by `cmp` after the fixes — 251,792 / 4,406 /
1,718,914 / 11,945 lines — and `ng_generic_walk_probe` still prints
`loci=236081 observations=251786 reads_admitted=54709` (`seconds=1.882`).

### The suite count moved, deliberately, and here is the accounting

The plan's Checkpoint A bar is "the suite is unchanged in count — a rename that changes a number
is not a rename". **The rename did not change it: 1,538 → 1,538 before fixes.** The review then
found two Blocker-class coverage gaps whose fix is tests, so the count is now **1,540**:

| test | why |
|---|---|
| `noodles_raw_aligned_read_decode_refuses_a_buffer_with_no_read_group` | B1 |
| `noodles_raw_aligned_read_default_reports_no_read_group` | B2 (kills both its mutations) |

and one test was **repurposed, not added** (`..._decode_errors_on_a_record_with_no_position` →
`..._decode_errors_on_a_stamped_record_with_no_position`), so M1 is net zero.

The rule's intent is that no number moves *unexplained*. These two are explained, named, and
each was verified to fail under the mutation it exists to catch. **Flagged for the owner at
Checkpoint A** — if the preference is a strictly count-neutral Milestone A, both tests lift out
cleanly and become a rider on C2, which is where this buffer changes owner anyway.

### Unresolved high-priority findings
None. The four `Deferred` items are all Minor and all have a named later home.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | `decode`'s read-group guard unpinned | Apply | Applied | No | `src/ng/read/aligned_read.rs` | Pass + mutation | No |
| B2 | Blocker | `Default` and `read_group()` unpinned | Apply | Applied | No | `src/ng/read/aligned_read.rs` | Pass + 2 mutations | No |
| M1 | Major | mis-named test never reaches its path | Apply | Applied | No | `src/ng/read/aligned_read.rs` | Pass + mutation | No |
| Mi1 | Minor | `input/mod.rs` doc misattributes the type | Apply | Applied | No | `src/ng/read/input/mod.rs` | Pass | No |
| Mi2 | Minor | `pub` carried across unexamined | Defer | Deferred | No | None | N/A | Yes — Checkpoint A |
| Mi3 | Minor | doc back-reference + "production `RecordSource`" | Apply | Applied | No | `src/ng/read/aligned_read.rs` | Pass | No |
| Mi4 | Minor | `..Default::default()` in two struct literals | Apply | Applied | No | `in_memory.rs`, `bam.rs` | Pass | No |
| Mi5 | Minor | `record_with`'s name and quals inert | Apply | Applied | No | `src/ng/read/aligned_read.rs` | Pass + mutation | No |
| Mi6 | Minor | impl report wrong in four places | Apply | Applied | No | the impl report | N/A | No |
| Mi7 | Minor | `mapq` is `MapQual` vs `u8` | Defer | Deferred | No | None | N/A | Yes — B or C |
| Mi8 | Minor | `CigarOp` back-reference into `pileup` | Defer | Deferred | No | None | N/A | Yes — layout doc |
| Nits | Nit | 11 items | Apply (8) / Defer (3) | Applied w/ adaptation | No | 4 files | Pass | Partly |

## 3. Questions asked and answers

None. Mi2 is a public-API decision and was **deferred to Checkpoint A** — where the owner is
already pausing — rather than asked as a separate blocking question. The plan-driven run
centralizes pausing at checkpoints.

## 4. Per-finding log

### B1 — `decode`'s read-group guard has no test that fails when the refusal is deleted
- **Severity:** Blocker · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** High confidence, mutation-proven by the `reliability` agent, one clearly
  correct implementation path (a test), no policy invention, no API change.
- **Implementation summary:** added
  `noodles_raw_aligned_read_decode_refuses_a_buffer_with_no_read_group`. Its fixture is a
  *decodable* record, and the test **asserts that first** — `decode_record(&record,
  ReadGroupId(0)).is_ok()` — so the missing stamp is provably the only thing that can fail it.
  That assertion is the whole difference from the test that failed to catch this.
- **Review suggestion used verbatim?:** No — adapted to use the renamed fixture helper.
- **Verification performed:** replaced the guard with
  `let read_group = self.read_group.unwrap_or(ReadGroupId(0));`, confirmed the replacement
  landed with `grep -c MUTATION-A` → `1`, ran the module: **FAILED, 8 passed / 1 failed**, and
  the failure is this test alone. Reverted; re-verified green.
- **Files changed:** `src/ng/read/aligned_read.rs`
- **Tests added:** `noodles_raw_aligned_read_decode_refuses_a_buffer_with_no_read_group`
- **Validation:** `cargo test --lib ng::read::aligned_read` → 0, 9 passed
- **Follow-up / residual risk:** None.

### B2 — the hand-written `Default` and the `read_group()` accessor are both unpinned
- **Severity:** Blocker · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** convergent (reliability + refactor_safety), both agents independently wrote the
  trap value in and got 2,837 green.
- **Implementation summary:** added `noodles_raw_aligned_read_default_reports_no_read_group`,
  reading a defaulted buffer **through the trait accessor** (`RawAlignedRead::read_group(&…)`)
  rather than through the field. That is what makes one test kill both mutations: the existing
  `in_memory.rs` / `region_records.rs` tests assert on the field, which is precisely why neither
  noticed.
- **Review suggestion used verbatim?:** No — the agent offered two tests; one through the
  accessor covers both.
- **Verification performed:** two mutations, each grep-confirmed before running.
  `Default`'s `read_group: None` → `Some(ReadGroupId(0))` (`MUTATION-B`): **FAILED, this test
  alone**. The accessor → `Some(self.read_group.unwrap_or(ReadGroupId(0)))` (`MUTATION-C`):
  **FAILED, this test alone**. Both reverted; re-verified green.
- **Files changed:** `src/ng/read/aligned_read.rs`
- **Tests added:** `noodles_raw_aligned_read_default_reports_no_read_group`
- **Validation:** `cargo test --lib ng::read::aligned_read` → 0, 9 passed
- **Follow-up / residual risk:** None.

### M1 — `..._decode_errors_on_a_record_with_no_position` tests neither thing its name claims
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** convergent, and it is the reason B1 and B2 survived to be found here.
- **Implementation summary:** **repurposed rather than deleted** — the test now stamps the
  buffer and blanks `alignment_start`, so the read-group guard cannot fire first and the
  decoder's own error is the one that surfaces. Renamed
  `..._decode_errors_on_a_stamped_record_with_no_position`, and it now asserts the **message**
  ("alignment start"), not just the kind, so it can only pass for the right reason. Its doc
  comment records what the old version actually did, so the trap is not re-set later.
- **Review suggestion used verbatim?:** No — adapted; the agent proposed a replacement, this
  keeps the test identity and fixes its subject.
- **Verification performed:** mutated `decode_record`'s `alignment_start` `ok_or_else(…)?` to
  `.map_or(1u64, …)` (`MUTATION-D`, grep-confirmed): **FAILED, 7 passed / 2 failed** — this
  test *and* `a_record_with_no_position_fails_naming_the_read`. The old version could not have
  failed under this mutation, because it never reached `decode_record`. Reverted.
- **Files changed:** `src/ng/read/aligned_read.rs`
- **Tests modified:** the renamed test
- **Validation:** `cargo test --lib ng::read::aligned_read` → 0, 9 passed
- **Follow-up / residual risk:** None.

### Mi1 — `input/mod.rs` doc says the type still lives in `filtering.rs`
- **Severity:** Minor · **Final status:** Applied (three-way convergent finding)
- **Implementation summary:** the sentence now names both homes — `RecordSource` from
  `filtering.rs`, and the buffer it fills, `RawAlignedRead`, from `aligned_read.rs` — and both
  links are crate-absolute defining paths rather than `super::` re-export paths, matching the
  sibling `record_reader/mod.rs`. The `super::` form is what let the prose drift while
  `cargo doc` stayed quiet.
- **Files changed:** `src/ng/read/input/mod.rs` · **Validation:** clippy + `cargo doc` clean
- **Residual risk:** None.

### Mi2 — `NoodlesRawAlignedRead`'s `pub` and its dead re-export
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** the agent verified the narrowing is clippy-clean, and the finding is correct.
  But `src/ng` is `pub mod ng`, so this is a **crate public-API change**, and this skill's rule
  is that public surface is not changed casually — `Ask` or `Defer`, never silent inclusion.
  A renames-only milestone is the wrong commit for it. Recorded for Checkpoint A, where the
  owner is already pausing. Note the agent also established the *other* half: the trait
  `RawAlignedRead` **cannot** be narrowed until Milestone C deletes `pub trait RecordSource`,
  so doing the struct alone now and the trait later would touch these lines twice.
- **Files changed:** None · **Follow-up:** Checkpoint A decision.

### Mi3 — the doc back-reference into `filtering.rs`, and "the **production** `RecordSource`"
- **Severity:** Minor · **Final status:** Applied
- **Reasoning:** the mislabel is a plain factual error — `RecordSource` is ng's trait, and the
  next paragraph exists to pin down exactly which side of the ng/production boundary things sit
  on. The intra-doc link is a genuine direction violation and a scheduled `cargo doc` breakage
  at Milestone C. **This reverses a call recorded in the impl report's §6**, which argued for
  keeping the link so the compiler would flag it at C; the reviewer's case — that
  `aligned_read.rs` should have zero references into `filtering.rs` so C's deletion is a
  one-file change — is stronger, and naming *what actually stamps the field* is more accurate
  prose besides.
- **Implementation summary:** three edits — the struct doc now says "whatever reads the file" /
  "whatever fills the buffer" and drops the link; `RawAlignedRead::decode`'s doc says "the read
  off the file"; the inline comment in `decode` says "whatever filled it". `aligned_read.rs` now
  has **no reference, code or doc, into `filtering.rs`**.
- **Files changed:** `src/ng/read/aligned_read.rs` · **Validation:** `cargo doc` clean
- **Residual risk:** None.

### Mi4 — `..Default::default()` in the two struct literals
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** both fields spelled out at `in_memory.rs` and `bam.rs`, each with
  a comment saying why. The irony the agent noted is real: both literals exist to plant a stale
  `read_group` and prove the reader clears it, and the spread is the one construct that would
  silently absorb a third field — on exactly the literal whose subject is the field being spread
  past. Milestone C moves this buffer onto `AlignmentCursor`.
- **Files changed:** `src/ng/read/input/record_reader/in_memory.rs`,
  `src/ng/read/input/record_reader/bam.rs` (which already imports `RecordBuf` in its test module)
- **Validation:** clippy clean, `cargo test --lib` green · **Residual risk:** None.

### Mi5 — `record_with`'s name and quality scores are inert
- **Severity:** Minor · **Final status:** Applied
- **Implementation summary:** two assertions added to the consuming test
  (`read.qname == b"r1"`, `read.qual == vec![30u8; 4]`), which also closes the gap the agent
  named: the adapter's pass-through of `qname` and `qual` was asserted nowhere (the parity
  oracle covers the free function `decode_record`, not `NoodlesRawAlignedRead::decode`). The
  helper's doc now says every field it sets is asserted somewhere, and why that matters.
- **Verification performed:** mutated the helper's quality scores to `vec![7u8; 4]`
  (`MUTATION-E`, grep-confirmed): **FAILED, this test alone**. Reverted.
- **Files changed:** `src/ng/read/aligned_read.rs` · **Residual risk:** None.

### Mi6 — four wrong claims in the impl report
- **Severity:** Minor · **Final status:** Applied (three-way convergent)
- **Implementation summary:** all four corrected in
  `ng_read_filtering_stages_a1_2026-08-03.md` — "Five call-site files" → "Six"; the `bam_record`
  claim rewritten to the verifiable one (one call site, and what clippy actually reported);
  the trait's rewritten opening sentence disclosed; §4's third table row now states what that
  test actually validated and points at this report for the repurpose.
- **Files changed:** the impl report · **Residual risk:** None.

### Mi7 — mapping quality is `MapQual` on the trait and `u8` on the struct
- **Severity:** Minor · **Final status:** Deferred
- **Reasoning:** the reviewer filed it as explicitly not-A1. `AlignedRead.mapq` is a public
  field; changing it touches every reader across read preparation and the STR path, and the
  acceptance dumps have to stay byte-identical. For B or C, whichever first edits the struct.
- **Follow-up:** carried to the plan's open items.

### Mi8 — `use crate::pileup::walker::CigarOp;`
- **Severity:** Minor · **Final status:** Deferred
- **Reasoning:** pre-existing, unchanged by this diff, and fixing it means touching frozen
  production code (~19 ng files inherit the path). The reviewer filed it as a note for
  `module_layout.md`'s open items, not a request against A1. Agreed.
- **Follow-up:** record against `arch/module_layout.md`.

### Nits — 11 items
- **Final status:** Applied with adaptation (8 applied, 3 deferred)
- **Applied:** merged the split `use` in `region_records.rs`; deleted the now-redundant
  `use noodles_sam::alignment::RecordBuf;` from `aligned_read.rs`'s test module; removed the
  four function-local flag constants in favour of production's canonical
  `FLAG_PAIRED` / `FLAG_MATE_REVERSE_STRAND` / `FLAG_REVERSE_STRAND` (the shadow was **created
  by this diff's new import**, so it is the diff's to clear); renamed `record_with` →
  `record_with_mapq_and_flags` so the call site reads; fixed the ungrammatical sentence in the
  new module doc; corrected the now half-false standing comment in `filtering.rs` ("this module
  knows what a *record* is" — A1 is precisely the commit that took that out); re-wrapped the
  over-long comment line in `read/mod.rs`; re-wrapped `input/mod.rs`'s paragraph as part of Mi1.
- **Deferred, with reasons:** `err` → `error` (subsumed — the two tests that used `err` were
  rewritten during B1/M1 and the new code uses `error` throughout); `FakeRecord` →
  `FakeRawAlignedRead` (Milestone C3 deletes it with `RecordSource` and its doubles, so the
  rename would be churn); bare `raw` bindings → `raw_read` (moved verbatim; not worth touching
  lines the milestone is trying to keep provably unchanged).

## 5. Deferred findings to carry forward
- **Mi2** — narrow `NoodlesRawAlignedRead` to `pub(crate)` and drop the dead re-export. Public
  API change; **owner decision at Checkpoint A**. Best done together with the trait, which
  cannot be narrowed until Milestone C.
- **Mi7** — `AlignedRead.mapq: u8` → `MapQual`, with one shared helper for the `0xFF` rule.
  For B or C.
- **Mi8** — lift `CigarOp` out of `pileup::walker`. For `module_layout.md`'s open items.
- **Nits (3)** — `FakeRecord`'s name (moot at C3) and the bare `raw` bindings.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — no `Apply` touched code reachable from any harness in `benches/`. The
  production-code changes are doc comments and one inline comment; everything else is test code,
  an import merge and an import deletion.
- **Baseline saved:** not applicable.
- **Outcome:** skipped — no `Apply` touched perf-sensitive code.
- **Notes:** the walk probe was run anyway, since this plan measures it at every step:
  `seconds=1.882` against `1.880` pre-fix and `1.846` at the base, all within run-to-run noise.

## 10. Commands run
- `cargo fmt`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib`
- `cargo test --lib ng::`
- `cargo test --lib ng::read::aligned_read`
- `cargo test --examples`
- `cargo build --release --examples`
- the four acceptance dumps + `ng_generic_walk_probe`, compared with `cmp` against the
  `8cf6f03` baseline
- five mutations, each `grep -c`-confirmed present before its test run and reverted after

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 2,839 passed / 0 failed / 5 ignored
- `cargo test --lib ng::` → 0, 1,540 passed / 0 failed / 2 ignored
- `cargo test --lib ng::read::aligned_read` → 0, 9 passed (was 7)
- `cargo test --examples` → 0, 52 passed / 0 failed
- four dumps → **byte-identical**; walk probe → anchor exact at `seconds=1.882`
- `MUTATION-A` → killed by B1's test only
- `MUTATION-B` → killed by B2's test only
- `MUTATION-C` → killed by B2's test only
- `MUTATION-D` → killed by M1's repurposed test *and* `a_record_with_no_position_fails_naming_the_read`
- `MUTATION-E` → killed by Mi5's strengthened test only
- final tree: `grep -c MUTATION src/ng/read/aligned_read.rs` → **0**

## 12. Notes
- **The suite count moved by +2 and that is the one thing a reviewer of Milestone A should look
  at twice.** The accounting is in §1. The rename itself was count-neutral; the review's two
  Blockers are what moved it, and both gaps are **pre-existing** — the moved impl is
  byte-identical to `8cf6f03`, as the `refactor_safety` agent proved by `sed`-and-diff. A1 is
  simply the commit that re-homes this code next to the tests meant to guard it.
- **Mi3 reverses a judgement recorded in the impl report.** Noted there and here rather than
  quietly changed.
- Five review agents ran in isolated worktrees; per-category findings are kept as an audit trail
  in the gitignored `tmp/review_2026-08-03_ng-read-filtering-stages-a1/`.
