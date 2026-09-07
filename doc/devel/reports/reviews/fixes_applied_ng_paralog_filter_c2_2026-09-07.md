# Fix Application Report: ng_paralog_filter_c2_2026-09-07.md

**Date:** 2026-09-07
**Source review:** [ng_paralog_filter_c2_2026-09-07.md](ng_paralog_filter_c2_2026-09-07.md)
**Source state reviewed against:** `3b55b94c` on `ng-paralog-filter`; **fixes folded into that
commit**, which had not been reported as finished — the review ran while it was still the head of
the branch and nothing had landed on top.
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
Blockers 1 · Majors 7 · Minors 12 · Nits 2

### Outcome totals
Applied 20 · Applied with adaptation 1 · Deferred 3 · Disputed 0 · Failed validation 0

### Validation summary

- `cargo test --all-features --lib "ng::run::paralog_filter"` → `ok. 110 passed; 0 failed`
- `cargo test --all-features --lib --bins --tests` → **6,582 lib tests passed, 0 failed, 15
  ignored**; one integration test failed, which is `main`'s
- `cargo clippy … -- -D warnings` → **9 errors, the merge base's exactly** — the two this step had
  added are gone. Counted with a command that groups *every* error kind, not the three the baseline
  had
- `cargo fmt --check` → the same 9 files as the merge base
- Performance check → not applicable; the review measured the one hot-path question directly

**The lib count moves 6,574 → 6,582**: eight tests added.

### Unresolved high-priority findings

None. The Blocker and all seven Majors are Applied.

## 2. Findings table

| ID | Severity | Title | Final status |
|---|---|---|---|
| B1 | Blocker | the alternative count untested against unexplained reads | **Applied** |
| M1 | Major | direct mode's guard has no test | **Applied** |
| M2 | Major | two clippy errors added, gate reported clean | **Applied** |
| M3 | Major | the contig is untested | **Applied** |
| M4 | Major | the ploidy is untested | **Applied** |
| M5 | Major | the zip's unequal-length rule is untested | **Applied with adaptation** |
| M6 | Major | `finish_parking` and both error variants untested | **Applied** |
| M7 | Major | the VCF is created before the sink is chosen; `unreachable!` recovery | **Applied** |
| Mi1 | Minor | `--paralog-fdr` unvalidated | **Applied** |
| Mi2 | Minor | the default is an unnamed literal, twice | **Applied** |
| Mi3 | Minor | the guard sits below the other input checks | **Applied** |
| Mi4 | Minor | `mod.rs` says nothing fills a spill | **Applied** |
| Mi5 | Minor | `spill.rs` documents the span rule | **Applied** |
| Mi6 | Minor | `PassOneError::Vcf` duplicates its parent's message | **Deferred** |
| Mi7 | Minor | `is_repeat_tract()` read twice | **Applied** |
| Mi8 | Minor | a third variant would not be caught in `entry_for` | **Deferred** |
| Mi9 | Minor | `pub` with no caller outside the module | **Deferred** |
| Mi10 | Minor | ~40 lines duplicated across subcommands | **Deferred** |
| Mi11 | Minor | `--paralog-filter-tag` help text misleads | **Applied** |
| Mi12 | Minor | "6580 filtered out" does not reconcile | **Applied** |
| Nits | Nit | two overstated sentences in the report | **Applied** |

## 3. The mutation ledger, before and after

The review ran sixteen mutations and **five survived**. All five are now killed, each verified on
the tree being committed and restored from a byte-compared backup:

| mutation | before | after |
|---|---|---|
| the contig hardcoded to `ContigId(0)` | survived all 8 | **1 of 12 fails** |
| the ploidy hardcoded to 2 | survived all 8 | **1 of 12 fails** |
| `alt_reads` as `DP − AD[0]` | survived all 8 | **1 of 12 fails** |
| the zip padding instead of truncating | survived all 8 | unrepresentable — the lengths are asserted equal |
| direct mode's guard deleted | survived all 171 | **its own test fails** |

## 4. Adaptations

- **M5** — the review proposed an `assert_eq!` at the boundary plus a test, and that is what
  landed; the adaptation is that the assertion makes the *mutation* meaningless rather than merely
  detected. Padding versus truncating is no longer a behaviour the function can have, so the
  finding is closed by construction and the test pins the refusal instead.

## 5. Deferred, with reasons

- **Mi6 — `PassOneError::Vcf`'s message duplicating `RunError::RecordNotWritten`'s.** Real, and
  cosmetic in a path that is currently unreachable from the run (the off arm's error travels as
  `RunError`, not `PassOneError`). C4 is the step that gives `PassOneError` a live caller and can
  see which message an operator actually gets.
- **Mi8 — `entry_for` branching on a bool where a third row shape would need a `match`.** The
  producer and the writer now disagree in strictness, which is the finding. Deferred because the
  fix belongs with whoever adds the third variant — spec §8's tract-aware allele term — and doing
  it now means inventing the shape of a rule that has not been designed. The spill's decoder
  carries the same constraint and says so in a marked comment; this is the same note in a second
  place.
- **Mi9 — `pub` items with no caller outside the module.** They have a caller two steps from now,
  and `pub(crate)` would have to be widened again at C4. Recorded so the milestone's end can
  narrow the whole subtree at once, which is where the same finding on B3 was also sent.
- **Mi10 — ~40 duplicated lines across the two subcommands.** Both subcommands are ~40 lines of
  near-identical setup either side of this step's change, and the duplication predates it. Sharing
  them is a change to two commands' shape for a reason this step did not create.

## 6. What the numbers said

| claim | verdict |
|---|---|
| 171 lines; 8 tests; 7 construction sites; 6,565 → 6,574; `ref 2, alt 7`; the guard mutation | CHECKED-CORRECT |
| **"9 clippy errors, unchanged"** | **WRONG — 11.** Corrected, and the two are now removed rather than merely reported |
| **"6580 filtered out"** | **WRONG — 6,581.** Did not reconcile with 6,574 + 15 |
| **"three lines of the writing path changed"** | **WRONG** — three *edits*, six lines |
| "no branch beyond the one that chose the sink" | **WRONG** — the match runs per record |
| "before anything is opened" | **WRONG then, true now** — the guard has moved above the reference read and the psp listing |

**The wrong clippy figure is the one worth keeping.** Its cause was the orchestrator's own gate
command, which grepped for the three lint kinds the baseline already had — so a *new* kind was
excluded by construction while the count kept matching. The command in §10 of the review counts
every kind. The lesson is in `reporting-in-chat`'s failure log: a regression check that can only
see failures it has seen before is not a check.
