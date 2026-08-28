# Fix Application Report: ng_parameters_file_b2_2026-08-28.md

**Date:** 2026-08-28
**Source review:** [ng_parameters_file_b2_2026-08-28.md](ng_parameters_file_b2_2026-08-28.md)
**Source state reviewed against:** the uncommitted B2 diff over `021f7d8a`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 2
- Majors: 5
- Minors: 8
- Nits: 3

### Outcome totals
- Applied: 12
- Applied with adaptation: 1
- Deferred: 4 (three to step B3, one to the owner at Checkpoint B)
- Awaiting the owner: 2 (the key names, and the slippage table's layout)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib ng::calling::parameters_file` → 0, **58 passed, 0 failed, 2 ignored**
- `cargo test --lib` → 0, **4,981 passed, 0 failed, 13 ignored**
- Performance check → skipped; nothing in the diff is reachable from `benches/`

### The five mutations the review proved survived were re-run against the fixed tree, and **all
five now fail a test**

Each was applied, the module suite run, the tree restored from a pristine copy and verified
byte-identical by `diff`:

| mutation | test that now fails |
|---|---|
| `was_declared_by_the_run` hard-coded `true` | `a_run_that_declared_no_batching_writes_the_flag_as_false` |
| `SharesOrigin::slipped_reads` written as `unwrap_or(0.0)` | `a_shares_origin_that_fitted_nothing_writes_no_slipped_reads_key` |
| `one_a_line` omits its `key = []` | `a_file_with_every_table_empty_writes_and_reads_back` |
| `a_float_array` returns `""` for an empty slice | the same test |
| the inbreeding rows emitted twice | `every_row_of_every_table_is_one_line` (and four others) |

The last is the point of M2: before the fix, only the unanchored lookup stood between that mutation
and a green suite, and it looked at the wrong table.

## 2. Findings table

| ID | Severity | Title | Decision | Final status |
|---|---|---|---|---|
| B1 | Blocker | two §5 absences written by untested code | Apply | Applied |
| B2 | Blocker | both empty-array paths untested | Apply | Applied |
| M1 | Major | a `u64` above `i64::MAX` has no TOML spelling | Apply | Applied with adaptation |
| M2 | Major | the layout test inspects the wrong `by_sample` | Apply | Applied |
| M3 | Major | the float doc's stated reason is false | Apply | Applied |
| M4 | Major | `shares_by_repeat_offset` does not say where it starts | Defer | **Deferred to B3** |
| M5 | Major | goal 3's own worked example does not work | Defer | **Deferred to B3 + Checkpoint B** |
| Mi1 | Minor | twenty-two keys the first reader had to ask about | Ask | **Awaiting the owner** |
| Mi2 | Minor | `reach = "inside"` — inside what | Ask | **Awaiting the owner** |
| Mi3 | Minor | "stratum" names three tables and no row | Defer | Deferred to B3 |
| Mi4 | Minor | the article convention is unstated | Apply | Applied |
| Mi5 | Minor | the `*_word` family are dropped possessives | Apply | Applied |
| Mi6 | Minor | two test-coverage gaps in the TOML pass | Apply | Applied |
| Mi7 | Minor | `[inbreeding]` called the only cohort-sized axis | Apply | Applied |
| Mi8 | Minor | the substitution row is 162 chars against §9's 146 | Note | Raised at Checkpoint B |
| — | — | four wrong numbers in the diff's prose | Apply | Applied |
| Nits | Nit | three | Apply (1) / Note (2) | Applied / noted |

## 3. Questions asked and answers

None asked mid-run. **Two are carried to Checkpoint B**: the key names the file's first reader had
to ask about (Mi1, Mi2), and whether the slippage table alone should be emitted as
`[[array-of-tables]]` headers instead of one long line a row.

## 4. Per-finding log — the ones that changed a decision

### M1 — a count no TOML integer can hold
*Applied with adaptation.* One `a_toml_integer` helper, and all eight `u64` sites routed through
it. **The reviewer proposed a `debug_assert!` plus a saturating release path; the assertion was
dropped.** It made the one path the helper exists for unreachable from the test suite — the new
test panicked on it — which is the same defect this review round is about. The value therefore
saturates in every build, and it announces itself: what a reader sees is 9,223,372,036,854,775,807,
which no evidence count can be. Refusing a count sitting at exactly `i64::MAX` on the way back in
is recorded as step C2's.

### M3 — the float doc's reason
*Applied.* The rule stands and the reason is replaced with the true one. `1` **does** deserialise
into an `f64` through this crate's derived reader — measured — so the crate's own round trip cannot
see the difference. What a bare `1` changes is the file's type for every *other* reader:
`toml::Value` types it as an integer, and so does every parser in every other language. A file
whose concentration is a float or an integer depending on who opened it is not the artefact goal 3
describes. **This matters beyond the doc**: C2's reader will be designed against the claim, and the
claim was wrong.

### M4 and M5 — deferred to B3, which is the next step
Both are things the file does not *say*, and B3 is the step that makes the file say things.
- **M4** — a comment above `length_spectrum_by_stratum` giving the offset convention. This is a
  widening of B3's contract, which reads "each defaulted value gets a comment beside it saying
  where the default came from". Recorded as a deliberate widening: it is the same machinery, one
  comment, and it closes the only edit in the file that is wrong without being invalid.
- **M5** — a comment saying that editing a value means changing its warrant to `supplied` and
  dropping its `observations`. Same widening, and it closes the hazard on the very edit spec §1.2
  goal 3 names. **What a comment cannot close is the two-hop join** — the calibration table is
  keyed by `read_group` and carries no library name, where the contamination table beside it
  carries one. That is a field on the shape, so it goes to Checkpoint B.

### Mi1 and Mi2 — the key names
*Awaiting the owner, by the plan's own rule.* The plan says the names are revised "the first time a
person reads a file this writer produced and has to ask what a key means", and until then the
provisional names stand. That reader now exists and its questions are in the review. **One of its
suggestions is disputed rather than deferred:** renaming `level_origin` to `level_warrant` would
break a settled decision — the module reserves *warrant* for the four-state ladder and says so, and
where a repeat tract's numbers record how much of a period's curve went into them, that is
*smoothing* and not a warrant.

## 5. Deferred findings to carry forward
- **M4, M5, Mi3** → step B3's comments, with M4 and M5 recorded as a widening of its contract.
- **Mi1, Mi2** → the owner's revision of the key names, at Checkpoint B.

## 6–8. Disputed / failed validation / blocked
- `level_origin` → `level_warrant`: **disputed**, reason above.
- None failed validation; none blocked.

## 9. Performance check
Skipped — nothing in the diff is reachable from any harness in `benches/`.

## 10–11. Commands run

| command | exit | result |
|---|---|---|
| `./scripts/dev.sh cargo fmt --check` | 0 | clean |
| `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `./scripts/dev.sh cargo test --lib ng::calling::parameters_file` | 0 | 58 passed, 2 ignored |
| `./scripts/dev.sh cargo test --lib` | 0 | 4,981 passed, 13 ignored |

## 12. Notes

- **Test counts:** the module went 54 → **58** passing; the library suite 4,977 → **4,981**.
- **The golden file did not change**, which is the point of the four new tests: every one of them
  varies the fixture rather than the writer, so they reach states the golden file's fixture does
  not hold.
