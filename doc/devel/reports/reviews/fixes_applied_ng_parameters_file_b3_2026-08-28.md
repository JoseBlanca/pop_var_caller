# Fix Application Report: ng_parameters_file_b3_2026-08-28.md

**Date:** 2026-08-28
**Source review:** [ng_parameters_file_b3_2026-08-28.md](ng_parameters_file_b3_2026-08-28.md)
**Source state reviewed against:** the uncommitted B3 diff over `83934abe`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 6
- Minors: 6
- Wrong numbers in the diff's prose: 3

*(The reliability agent's final report arrived after its first draft and raised the two origins no
test reaches from a Minor to a Blocker, with a second half — the layout test's stated invariant is
false for any run that defaults a row in one of the tables it checks.)*

### Outcome totals
- Applied: 17
- Deferred: 0
- Disputed: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib ng::calling::parameters_file` → 0, **63 passed, 0 failed, 2 ignored**
- `cargo test --lib` → 0, **4,986 passed, 0 failed, 13 ignored**
- Performance check → skipped; nothing in the diff is reachable from `benches/`

### The five surviving mutations were re-run against the fixed tree, and **all five now fail a test**

Applied, module suite run, tree restored from a pristine copy and verified byte-identical by `diff`:

| mutation | test that now fails |
|---|---|
| the flat-concentration note fires whatever the warrant | `a_run_that_defaulted_nothing_writes_no_per_row_notes` |
| a per-row note counted as a row | `every_row_of_every_table_is_one_line` |
| the substitution rate's origin crossed with the inbreeding one | `every_defaulted_number_says_where_its_default_came_from` |
| the inbreeding origin crossed with the calibration one | the same |
| the two defaulted scalars' origins swapped | `every_defaulted_number_says_where_its_default_came_from` (was the golden file alone) |
| the wrapper's boundary moved by one, and its unit changed to bytes | `a_note_wraps_by_characters_at_the_width_it_states` |

## 2. Findings table

| ID | Severity | Title | Decision | Final status |
|---|---|---|---|---|
| B1 | Blocker | two origins reach no document any test writes, and the layout test's invariant is false where a row is defaulted | Apply | Applied |
| M1 | Major | contamination note: wrong substance, wrong grain | Apply | Applied |
| M2 | Major | "each value carries a warrant" is false | Apply | Applied |
| M3 | Major | `rung` defined once, means two things | Apply | Applied |
| M4 | Major | the batching note refuted by its own rows | Apply | Applied |
| M5 | Major | the flat-concentration note fires unconditionally | Apply | Applied |
| Mi1 | Minor | two origins never exercised | Apply | Applied |
| Mi2 | Minor | the wrapper's boundary and unit unpinned | Apply | Applied |
| Mi3 | Minor | each note asserted somewhere, not above its key | Apply | Applied |
| Mi4 | Minor | `[repeat_tracts]` never says what its numbers are | Apply | Applied |
| Mi5 | Minor | "one chromosome's worth" beside a value of 1.25 | Apply | Applied |
| Mi6 | Minor | seven lines to cut, five to spend | Apply | Applied |
| — | — | three wrong numbers in the report | Apply | Applied |

## 3. Questions asked and answers

None. Every finding was the coder's to fix.

## 4. Per-finding log — what changed in the file

- **M1.** The note now reads "how much of each **read group's reads** came from somebody else — one
  row a **lane**, because two lanes of one library can differ: index hopping happens on a flowcell,
  not in a tube", and the editing instruction says "to stop correcting one lane … a library
  sequenced over several lanes has a row for each". The state bullet that said "this library could
  not be measured" says "this lane".
- **M2.** "**A number that could be fitted** carries a `warrant`", followed by the rule, followed by
  a new paragraph naming the numbers that carry none — the slippage numbers, the prior's two
  concentrations, the length spectra — and saying there is nowhere in them to record an edit, so
  note it elsewhere.
- **M3.** One clause added where the second `rung` lives: "a `rung` inside one of those curves is
  not the `rung` in `[ordinary_site_prior]`: here it says what the curve itself was fitted on".
- **M4.** "A declared batching that happens to have one batch writes identical rows, and this flag
  is the only thing that tells those two apart" — the narrower claim, which the rows do not refute.
- **M5, Mi1, Mi3.** The two tests were strengthened rather than added to: the absence test now lists
  all five origin texts, and the presence test now (a) locates each note **above its own key** by
  walking back from the key through the comment lines, and (b) defaults a substitution rate and an
  inbreeding coefficient so the two origins nothing reached are exercised.
- **Mi2.** `a_note_wraps_by_characters_at_the_width_it_states` — the room a note is given is one
  line, one more character is two, and a note of that many accented characters (117 bytes) is still
  one line. Also that a row's note has six characters less room than a section's, a single word
  longer than the width, an empty note, and one of only spaces.
- **Mi4.** Five lines added defining `level`, `shorter_share` and `fall_off`, and saying why
  `slipped_reads` is fractional: it is how many reads the fitted level says slipped, not a count
  anybody labelled.
- **Mi5.** "**this many** chromosomes' worth of belief", which is true whatever the value on the
  line.
- **Mi6.** Seven lines cut — the inbreeding section's design rationale about sample order, and the
  stated-constants framing — against the five spent on Mi4.

## 5–8. Deferred / disputed / failed / blocked

None.

## 9. Performance check

Skipped — nothing in the diff is reachable from any harness in `benches/`.

## 10–11. Commands run

| command | exit | result |
|---|---|---|
| `./scripts/dev.sh cargo fmt --check` | 0 | clean |
| `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `./scripts/dev.sh cargo test --lib ng::calling::parameters_file` | 0 | 63 passed, 2 ignored |
| `./scripts/dev.sh cargo test --lib` | 0 | 4,986 passed, 13 ignored |

## 12. Notes — and one thing the fixes changed beyond the findings

**The wrapper now knows what a line costs before the words start.** The review found that
`wrapped` budgeted 78 columns of prose while a per-row note adds six characters of prefix, and the
width test capped lines at 82 — so a correct writer whose origin wrapped to 77 columns would have
failed the suite. There is now one stated width, `COMMENT_WIDTH = 80`, and two constants derived
from it for the two places a note can sit. **Every note in the file, hand-written or generated, now
goes through the wrapper**: I had been hand-tuning section notes to fit a column, which is what a
wrapper is for, and three of them drifted past it twice during this step. Lines that are laid out
on purpose — the bullet list where spec §5's three contamination states are set against each other
— are emitted as they stand.

**And the anchoring defect the review of B2 found reappeared while fixing this one.** A new
assertion searched for `by_sample = [` from the top of the file and found the batching table rather
than the inbreeding one. Same trap, same file, twenty minutes apart; it failed loudly because the
assertion under it was specific enough to notice.

- **The produced file went 148 → 162 lines** across the fixes: 81 comment lines against the 67 the
  review read, because M1, M2, M3 and Mi4 each added prose and only Mi6 removed any.
- **Every one of the four wrong comments was wrong in the same way** — a claim made confidently
  about something the file itself, or the spec, states otherwise. Three were checkable against
  lines within twenty of the comment.
