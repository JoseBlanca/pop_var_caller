# Fix Application Report: ng_psp_mode_c2_c3_2026-09-03.md

**Date:** 2026-09-03
**Source review:** `doc/devel/reports/reviews/ng_psp_mode_c2_c3_2026-09-03.md`
**Source state reviewed against:** 7a164243 + the uncommitted C2/C3 diff, branch `ng-psp-mode`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 7
- Minors: 10

### Outcome totals
- Applied: 7 Majors, 7 Minors
- Deferred: 3 Minors — §5
- Applied with adaptation: 1 (M7 — see the log)

### Validation summary
- `cargo fmt --check` → 0; `cargo clippy --all-targets --all-features -- -D warnings` → 0
- `cargo test --lib 'pop_var_caller_exp'` → 0, **125 passed** (was 119; the command's own tests
  are now 37, from 31)
- **Run for real** after the fixes, stdout only:
  ```
  walked 1 sample over 2 intervals, SL4.0ch01:3406887-3506886 … SL4.0ch01:13806670-13906669 — 200000 bases asked for
    SRS3394712: 193603 loci stored, 914723 bytes at tmp/c2_psps/SRS3394712.psp; spoke for 311 of 318 typed regions (199672 of 200000 bases walked, 99.8%); not stored — clusters of repeats too close together to have clean flanks: 328 bases (0.2%)
  1 psp, 914723 bytes in total
  ```
- Performance check → skipped: nothing on a `benches/` path changed.

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Title | Final status |
|---|---|---|
| M1 | "observations" names loci | Applied |
| M2 | share taken over the wrong whole | Applied |
| M3 | the C2 split undone by a bare `println!` | Applied |
| M4 | C3's guard defanged a C1 test | Applied |
| M5 | C3's headline property unpinned | Applied |
| M6 | advisory guard + colliding scratch name | Applied |
| M7 | the report's content barely pinned | Applied with adaptation |
| Mi1–Mi7 | see §4 | Applied |
| Mi8–Mi10 | see §5 | Deferred |

## 3. Questions asked and answers
None — non-interactive. The review's two open questions were resolved as recorded there: the
progress line stays but moves to stderr and shares one formatter with the report; the report
keeps five of `LocusCounts`' eight fields.

## 4. Per-finding log

- **M1 — Applied.** "observations" → "loci stored" everywhere the record count is named. The
  word now means in the report what it means in the crate.
- **M2 — Applied.** `SampleWalkOutcome::bases_walked()` is the three parts' own sum, and every
  share is of it; where it differs from the BED's ask the report says so on its own line, which
  is the shape `run/report.rs` arrived at after the 200.0% measurement. Its `describe`,
  `plural` and `share_of` are now `pub(crate)` and reused rather than reimplemented, so both
  commands render the same arithmetic in the same words.
- **M3 — Applied.** `SampleWalkOutcome::line()` is the single formatter; the progress note
  prints it to **stderr** as each sample finishes and the report prints it to stdout at the
  end. They cannot drift, and a shell capturing stdout gets the report alone.
- **M4 — Applied.** The defanged test gets `force: true` and now asserts the error is
  `Walk { .. }`, so it fails if it ever stops at the door again. Its comment records that it
  did.
- **M5 — Applied.** The stopped-walk test plants a stale `.partial` for this process before
  running, so the cleanup has something observable to remove; the assertion now fails if the
  cleanup is deleted.
- **M6 — Applied.** The scratch file is `<sample>.psp.<pid>.partial`, so two invocations of one
  sample no longer interleave into one file — each writes its own and the rename is atomic, so
  the loser replaces a whole psp with a whole psp. The comment says plainly that this is not a
  lock and racing on one sample is still a thing not to do. `exists()` → `try_exists()`, so
  "cannot tell" is an error rather than "no file".
- **M7 — Applied with adaptation.** The review asked for the eight numbers pinned in the
  command-level test; that test's fixture covers **all** of its ground, so a swap of the two
  region counts is invisible there and a guard asserting otherwise fails. The formatter is
  pinned instead on a hand-built outcome where no two numbers are equal, plus a second test
  that the two uncovered-ground clauses are absent when their counters are zero. The guard
  that caught this is kept in the report as a note.
- **Mi (applied):** the refusal now names the blocked sample, pinned by a test where only the
  *second* sample's psp is in the way; `--force`'s spelling tied to clap by a test; the two
  unhandled kinds rendered in the sibling's plain English ("clusters of repeats too close
  together to have clean flanks", "tandem arrays longer than this run types as callable"); the
  report names the ground it walked, pinned by a test; the `--force` test no longer builds a
  throwaway cohort whose `TempDir` guards die on the same line; `write!` replaces
  `push_str(&format!(…))` and the comment's wrong cure is gone.

## 5. Deferred findings to carry forward
- **Nothing says a psp was *replaced* when `--force` acts** — the run's output is identical
  either way. Worth a line, but it needs the existence check to keep its answer rather than
  short-circuit, and that is better done alongside whatever C-milestone follow-up wants a
  fuller run report.
- **`WalkReport`/`SampleWalkOutcome` are `pub` with a private constructor and no `Debug`** —
  they become reachable when something outside the command consumes a walk report; until then
  the visibility is anticipatory, as `SampleObservationGatherer`'s was.
- **`SampleWalkOutcome`'s field names diverge from `SampleWalkTallies`'** (`sample` vs
  `sample_name`, `counts` vs `regions`) — a rename worth doing across both at once.

## 6–8. Disputed / Failed validation / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched code reachable from `benches/`.

## 10–11. Commands run and results
In the container: `cargo fmt` / `--check` → 0; `cargo clippy --all-targets --all-features --
-D warnings` → 0; `cargo test --lib 'pop_var_caller_exp'` → 0 (125 passed, 1 ignored); the real
run quoted in §1.

## 12. Notes
One process error worth recording: an edit meant for `src/ng/run/report.rs` in this worktree
ran in the **main checkout** instead, because the shell's working directory had reset between
calls. It was three visibility keywords, caught on the next compile, and reverted with
`git checkout --` — the main checkout is back to its original state. Every later command in
this run uses an absolute path or an explicit `cd` to the worktree.
