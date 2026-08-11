# Fix Application Report: ng_parameter_prepass_generic_a4_2026-08-06.md

**Date:** 2026-08-06
**Source review:** `doc/devel/reports/reviews/ng_parameter_prepass_generic_a4_2026-08-06.md`
**Source state reviewed against:** `54378b9` (ng-parameter-estimation)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 3
- Minors: 10
- Nits: 4

### Outcome totals
- Applied: 14
- Applied with adaptation: 2
- Already fixed: 0
- Deferred: 1
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib ng::parameter_estimation` → 0, **20 passed** (was 13)
- `cargo test --all-targets --all-features` → 101, **2,920 passed, 1 failed, 5 ignored**. The failure is the pre-existing `ng::locus_generation::pileup::parity` divergence.
- `cargo doc --no-deps --lib` → 101, 12 unresolved links, all pre-existing, none in this file.
- Performance check → not applicable; nothing under `parameter_estimation/` is reachable from `benches/`.

### Unresolved high-priority findings
None. The one deferral is a design-document edit outside this skill's remit.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Validation |
|---|---|---|---|---|---|
| M1 | Major | Degeneracy guard dead; clamps ordered backwards | Apply | Applied with adaptation | Mutation-verified |
| M2 | Major | `row_start` answers where `depth_range` refuses | Apply | Applied | Two `should_panic` tests |
| M3 | Major | Report overstates the oracle (2/3/3 → 2/2/2) | Apply | Applied | Re-measured |
| Mi1 | Minor | Cap's cost quoted at the wrong bin count | Apply | Applied | Checked against research §4.3 |
| Mi2 | Minor | `EXACT_DEPTH_LIMIT` doc leads with the wrong quantity | Apply | Applied | Read back |
| Mi3 | Minor | `max_depth()`'s uniqueness claim contradicted | Apply | Applied | Constants privatised |
| Mi4 | Minor | No build-time check on the shape constants | Apply | Applied | `const _: () = assert!` |
| Mi5 | Minor | Vacuously-true `PartialEq` where `Arc::ptr_eq` is mandated | Apply | Applied | Derive dropped |
| Mi6 | Minor | "330 MB" arithmetic; "costs nothing" stated as settled | Apply | Applied | Checked against spec §9, research §4.6 |
| Mi7 | Minor | Design docs still name four files under `generic/` | Defer | Deferred | — |
| Mi8 | Minor | `pub` constants; no `#[must_use]` | Apply | Applied | clippy |
| Mi9 | Minor | Two `mut` loops where `scan` says it | Apply | Applied | 20 tests pass |
| Mi10 | Minor | `Default` and `bin_for` totality undertested | Apply | Applied | Two new tests |
| N1 | Nit | `bin_for`'s fast path unexplained | Apply with adaptation | Applied with adaptation | Kept + documented |
| N2 | Nit | The `+1` offset re-derived at three sites | Apply | Applied | `FIRST_WIDENING_BIN` |
| N3 | Nit | `.expect()` without `// PANIC-FREE:` | Apply | Applied | Read back |
| N4 | Nit | `as u16` truncation above 65,535 bins | Apply | Applied | Covered by the const assertion |

## 3. Questions asked and answers

None asked. Mi7 is recorded as an owner item — see §5.

## 4. Per-finding log

### M1 — the degeneracy guard (Applied with adaptation)

- **Implementation.** The generator moved into a private `ladder_tops(exact_limit, cap, bin_count)`, so the clamp is reachable from a test with shapes the adopted constants never produce. The clamps are reordered to `top.min(cap).max(previous + 1)`.
- **Adaptation, and it changed the shape of the fix.** The review proposed the reorder plus an assertion that the tops are strictly increasing. Writing it revealed that the reorder makes strict increase **provable by construction** — `previous + 1` always exceeds `previous` — so that assertion is unreachable, which is precisely the dead-guard defect being fixed. What is left to check is the other failure: the `+1` walk running past the cap. So `ladder_tops` asserts one thing, that the ladder ends exactly at its cap, and the doc states strict increase as a consequence of the clamp order rather than as something checked.

  This was caught by the first two `should_panic` tests failing for the wrong reason, not by reasoning ahead.
- **Verification.** Reverting the clamp order to the shipped one turns two green tests red: `ladder_tops(8, 12, 20)` then builds a ladder that *does* end at its cap but with repeated tops, so the old guard passed it. That is the defect, reproduced and now caught.
- **Tests:** `a_cap_too_low_for_the_bin_count_is_rejected`, `more_bins_than_the_cap_has_room_for_is_rejected`, `every_ladder_shape_is_strictly_increasing` (three shapes, including one where the geometric step rounds below a whole depth so the `+1` clamp does all the work).

### M2 — `row_start` (Applied)

A range assertion naming the bin and the ladder's size, matching what `depth_range` already did by accident of indexing. `depth_range` gained the same explicit assertion so the two agree and both messages name the fault. Two `should_panic` tests, one per accessor.

### M3 — the report's oracle claim (Applied)

Corrected to 2 / 2 / 2, with the narrower true statement: the two literal-bearing tests are the oracle, and the other five check internal consistency. The same statement is now in the doc comment of `the_widening_bins_top_out_at_the_measured_depths`, where the next reader meets it. `the_bin_ranges_partition_the_depths_from_zero_to_the_cap` also gained a literal `assert_eq!(edges.max_depth(), 124)` beside its constant-relative one, so a third test now catches a moved cap.

**The commit message of `54378b9` carries the wrong 2/3/3 figures and cannot be corrected without rewriting history.** Recorded here and in the impl report instead.

### Mi1–Mi6, N1–N4 — documentation and shape (Applied)

The cap's cost restated at twenty bins (0.054→0.190 rungs, 0.30%→0.88%, 3.5-fold and 2.9-fold), with the sixteen-bin figures kept and labelled as the architecture's comparison. `EXACT_DEPTH_LIMIT` re-led as "the deepest depth that keeps a bin to itself … which is nine bins". The three constants made private, which resolves Mi3 — `max_depth()` is now genuinely the only statement of the cap a consumer can reach — and its doc says so. A `const _: () = assert!` guards the shape at build time, mirroring the sibling ladder in `generic/mod.rs`, and covers the `as u16` truncation (N4) as a side effect. `PartialEq`/`Eq` dropped from `DepthBinEdges` with a doc paragraph saying `Arc::ptr_eq` is the sameness check and a value comparison cannot stand in for it. The 330 MB arithmetic now includes the ×8,000-windows step; `MAX_BINNED_DEPTH` says plainly that what losing the reads above costs is **not measured** and cites research §4.6. `scan` replaces both `mut` accumulator loops, removing one `.expect()`; the two that remain carry `// PANIC-FREE:` comments. `FIRST_WIDENING_BIN` names the offset once. `#[must_use]` on every accessor.

**N1, adapted.** The review showed a single `partition_point` over the whole ladder gives identical answers and deletes the repeated offset arithmetic. The fast path is **kept**: it is where 97 sites in 100 of a three-read cohort land, and `bin_for` runs once per covered position over hundreds of millions of them, so one comparison against four or five is worth a branch. What was missing was the reason — the code now says it, which was the reviewer's own stated alternative.

### Mi7 — the design documents (Deferred)

`arch/parameter_prepass_generic.md`'s module table names four files under `generic/` and maps §2.2 to `histogram.rs`; the plan's A1 names the same four. All three reviewers agree the fifth file is right, so the documents are behind the code. **Editing a spec, architecture or plan is outside this skill's remit** — it changes design, not code. Raised as an owner item in `PROJECT_STATUS.md`, with the suggested module-table diff in the review report.

## 5. Deferred findings to carry forward
Mi7 — add `depth_bins.rs` to the architecture's module table and the plan's A1 file list. Owner decision.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive code. Nothing under `src/ng/parameter_estimation/` is reachable from any harness in `benches/`; the module has no consumers yet. `bin_for`'s fast path was kept on a hot-path argument rather than a measurement, and is flagged as such: it costs one comparison against four or five, and will be measurable once Milestone C feeds it real loci.

## 10. Commands run
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::parameter_estimation`
- `./scripts/dev.sh cargo test --all-targets --all-features`
- `./scripts/dev.sh cargo doc --no-deps --lib`
- The clamp-order revert, to confirm the new tests catch it

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, no output
- `cargo test --lib ng::parameter_estimation` → 0, 20 passed
- `cargo test --all-targets --all-features` → 101, 2,920 passed / 1 failed (pre-existing) / 5 ignored
- `cargo doc --no-deps --lib` → 101, 12 unresolved links, all pre-existing
- clamp order reverted → 2 tests fail, confirming the guard is live

## 12. Notes

- **The defect this step was isolated for was real, and only mutation found it.** Three reviewers read the same twenty lines; the one that instrumented the generator and counted clamp firings found that both clamps fire zero times, and the one that drove other shapes found that when they do fire they produce the empty bins the comment says they prevent. Neither is visible by reading.
- **The fix's own first attempt was wrong**, and the tests caught it: reordering the clamps makes the strictly-increasing assertion unreachable, so the first version of the guard reproduced the dead-branch defect it was fixing. Recorded because it is the same failure mode one step later.
- **The oracle is two tests.** That is now stated where a reader meets it rather than only in a report — the other five check that whatever ladder is configured is internally consistent, which is a different question from whether it is the ladder that was measured.
