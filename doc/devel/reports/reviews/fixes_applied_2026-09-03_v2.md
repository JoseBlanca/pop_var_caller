# Fix Application Report: ng_psp_mode_a2_a3_a4_2026-09-03.md

**Date:** 2026-09-03
**Source review:** `doc/devel/reports/reviews/ng_psp_mode_a2_a3_a4_2026-09-03.md`
**Source state reviewed against:** branch `ng-psp-mode`, uncommitted diff over 114efe24
**Execution mode:** non-interactive (plan-driven bundled step A2+A3+A4 loop)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1 · Majors: 3 · Minors: 9 · Nits: 1 grouped set

### Outcome totals
- Applied: 12 (B1, M1, M2, M3, Mi1, Mi2, Mi4, Mi5, Mi6, Mi7, Mi8, Mi9)
- Applied with adaptation: 1 (Mi3)
- Deferred: 0 · Disputed: 0 · Failed validation: 0
- Nits: the `_bp` suffix, the redundant `.to_string()` disappearance, the `"off"` const, and
  the doc line on unconditional drops all landed inside M2's move; test locals and parameter
  naming left as-is (cosmetic, would churn the diff); checkboxes flipped after commit.

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-targets --all-features` → 0, 16 binaries; lib **6,055 passed; 0 failed**
  (was 6,052 pre-fixes; `ng::psp` 417, `ng::read::filtering` 21)
- `cargo doc` / `cargo audit` → not run (no dependency change)
- Performance check → skipped: changed code is once-per-file header encode/decode plus
  test/fixture code; the benches time per-record loops no fix touches. No baseline saved,
  consistent with the skip.

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| B1 | Blocker | provenance test blind on four of six values | Applied | filtering.rs (new pinning tests), header.rs (old test superseded) |
| M1 | Major | non-exhaustive config field access | Applied | filtering.rs (exhaustive destructure inside the moved function) |
| M2 | Major | psp→read dependency edge | Applied | filtering.rs (`provenance_parameters` + `READ_FILTER_PROVENANCE_KEYS`), read/mod.rs (export), header.rs (generic `record_parameters`; `wire_float_of` moved home — Mi8) |
| M3 | Major | duplicated-@RG-ID acceptance unpinned | Applied | header.rs (`two_read_groups_sharing_an_rg_id_round_trip`) |
| Mi1 | Minor | `Bp(65_535)` literal ×3 + mislabeling comment | Applied | both examples, the bench (named `MAX_RECORD_SPAN_CEILING`, comment corrected) |
| Mi2 | Minor | `identifier` vs wire `walk-local-id` | Applied | field renamed `walk_local_id` everywhere |
| Mi3 | Minor | stale floor key on re-record | Applied with adaptation | the key family is published (`READ_FILTER_PROVENANCE_KEYS`) and `record_parameters` documents the insert-only contract; no clearing caller exists yet (B1 of the plan calls once) |
| Mi4 | Minor | order fixture far from boundary | Applied | header.rs (identifier-repeating-zero rule row added, the 7 kept) |
| Mi5 | Minor | control-char rule pinned only at newline | Applied | header.rs (tab rule row) |
| Mi6 | Minor | no ceiling fixture at `MAX_TOML_INTEGER` | Applied | header.rs (widest-number test extended) |
| Mi7 | Minor | digits-as-string fallback untested | Applied | filtering.rs (u64::MAX assertion) |
| Mi8 | Minor | `wire_float_of` in the wrong file | Applied | moved beside `hex_of`/`digest_of` in header.rs |
| Mi9 | Minor | impl-report prose overstatements | Applied | the A2–A4 impl report corrected (fixture reach, test claims, renamed fields, re-measured counts) |

## 3. Questions asked and answers
None; open question 1 of the review (which ceiling value the real walk records) is B1-of-the-plan's
to answer and is noted in the report for that step's review.

## 4. Per-finding log — the decisions that were not verbatim

- **M2/M1 (one change):** the enumeration became `ReadFilterConfig::provenance_parameters()
  -> Vec<(String, ParameterValue)>` in `ng::read` — a stage importing psp (infrastructure) is
  the sanctioned direction, as with `ref_seq` — destructuring the config exhaustively so a new
  filter fails to compile at the recording site. psp's `WriterProvenance` keeps only the
  generic `record_parameters(entries)`. The B1-strength value tests moved with the logic to
  `filtering.rs`, where the values live.
- **Mi3 (adaptation):** rather than a `remove` call inside psp (which no longer knows the
  keys), the key family is published as `READ_FILTER_PROVENANCE_KEYS` with both sides'
  docs stating the contract: insert-only, and a re-recorder whose key set may have shrunk
  clears the family first. No second-record caller exists; the plan's B1 records once.
- **Mi2 with the idiomatic alternative declined:** the field stays (renamed `walk_local_id`)
  rather than being derived from position — a row that travels alone into E2's merge stays
  self-describing, and the position-equality redundancy remains checked on both sides.
- **Nit applied opportunistically:** `observation_reach_ceiling` → `observation_reach_ceiling_bp`,
  matching its sibling `genomic_block_size_bp` (mechanical; sed in the container, verified by
  the full suite).

## 5–8. Deferred / Disputed / Failed / Blocked
None.

## 9. Performance check
Skipped — no Apply touched per-record code (reason in §1).

## 10–11. Commands run and results
- `./scripts/dev.sh cargo fmt --check` → 0
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` → 0
- `./scripts/dev.sh cargo test --lib 'ng::psp'` → 0, "417 passed; 0 failed"
- `./scripts/dev.sh cargo test --lib 'ng::read::filtering'` → 0, "21 passed; 0 failed"
- `./scripts/dev.sh cargo test --all-targets --all-features` → 0, 16 binaries, lib
  "6055 passed; 0 failed; 14 ignored"

## 12. Notes
The nine reviewers were all interrupted mid-run by a session rate limit and resumed in place;
their worktrees and applied patches survived, and every findings file was completed after the
resume.
