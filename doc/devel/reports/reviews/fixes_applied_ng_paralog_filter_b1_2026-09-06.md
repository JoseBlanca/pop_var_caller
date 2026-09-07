# Fix Application Report: ng_paralog_filter_b1_2026-09-06.md

**Date:** 2026-09-06
**Source review:** [ng_paralog_filter_b1_2026-09-06.md](ng_paralog_filter_b1_2026-09-06.md)
**Source state reviewed against:** `2ad2a92b` on `ng-paralog-filter`; fixes folded into that commit
**Execution mode:** non-interactive
**Overall status:** Completed

---

## The answer

**The round-trip comparator can now be shown to fail, and a field the decoder drops no longer
compiles.** Those are the Blocker and the first Major, and they were the same nine tests'
problem: `same_bits` was the only assertion nine round trips made, it compared its fields by
access rather than by destructure, and nothing anywhere showed it could return `false`. Hard-
wiring it to `return true` used to leave all seventeen tests green; it now fails one, by name.

**Two mutations that used to survive are dead**, each re-run in the container and each reverted
from a byte-compared backup:

| mutation | before | after |
|---|---|---|
| the comparator hard-wired to `return true` | all 17 pass | `holds_the_same_bits_separates_entries_that_differ_in_any_one_field` fails |
| a field added to `SpilledSample`, encoded correctly, decoded as `0`, with the fixtures and the pinned byte list updated as a coder resolving the errors would | all 17 pass | **does not compile** — two `E0027` inside the comparator |
| the line's length read as one byte instead of a varint | all 17 pass | three fail, including the new property |
| the reader's failure latch removed, so it resynchronises | *(new behaviour)* | `the_reader_stops_after_a_decode_error` fails |

**The decoder no longer pulls the rest of the file into memory on a corrupt line length.** The
extras agent measured 16,777,230 bytes read into one `Vec` before the short-read check could
fire; a 16 MiB ceiling refuses the length before a byte of it is read. Two other agents probed
the same field over a file holding *no* payload and saw no spike — correctly, and that is the
mechanism: the bound was whatever the file still held, which on a genome's spill is the rest of
it.

**Three of my own numbers were wrong** and are corrected in the implementation report and the
commit message: the report said "seven" where its own sentence lists eight; both the report and
the commit message said "fifteen round trips are blind to a layout wrong on both sides" where
there are nine; and the commit message said the swapped-floats mutation "kills only the byte-list
test" where the report's own table says two. Every other figure — thirty of thirty-three checked
— was right, and every one of the three wrong was my count of my own tests.

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 11
- Minors: 21
- Nits: 13

### Outcome totals
- Applied: 36
- Applied with adaptation: 3
- Already fixed: 0
- Deferred: 4
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo test --all-features --lib "ng::run::paralog_filter::"` → 0, `ok. 29 passed; 0 failed`
- `cargo test --all-features --lib --bins --tests` → non-zero, 6,387 lib tests pass; the one
  integration failure is `main`'s
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` → non-zero, 3 errors, all
  `needless_lifetimes` in `cohort_merge/`, `main`'s
- `cargo fmt --check` → non-zero, dirty on 9 files, none of them this module's, `main`'s
- `cargo doc --no-deps` → not run; red on `main` for 11 unresolved intra-doc links elsewhere
- `cargo audit` → not run; no offline advisory database in the container
- `cargo test --all-targets --all-features` → not run; `examples/ng_candidate_selection_probe.rs`
  has not compiled on `main` since `a33ada0f`, so the command measures nothing
- Performance check → **skipped**: nothing under `src/ng/run/paralog_filter/` is reachable from
  any harness in `benches/`, and nothing writes or reads a spill yet.

### Unresolved high-priority findings
- **M10 (part)** — whether `SpilledSample` keeps the spec's flat shape or becomes a sum type is
  the owner's call; the impossible flag pair is refused meanwhile.
- **M9 (part)** — nothing yet carries the writer's entry count to the reader, so a spill that
  lost its tail on an entry boundary still reads back as a complete, shorter spill. B2 owns the
  file's lifecycle and is where the number can be threaded.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | the comparator cannot be shown to fail | Apply | **Applied** | No | `spill/tests.rs` | Mutation killed | No |
| M1 | Major | a corrupt line length pulls in the rest of the file | Apply | **Applied** | No | `spill.rs`, `spill/tests.rs` | Pass | No |
| M2 | Major | the comparator reads its fields by access | Apply | **Applied** | No | `spill/tests.rs` | Mutation now a compile error | No |
| M3 | Major | a decode error is not terminal | Apply | **Applied** | No | `spill.rs`, `spill/tests.rs` | Mutation killed | No |
| M4 | Major | no fixture crosses the 128-byte varint boundary | Apply | **Applied** | No | `spill/tests.rs` | Mutation killed | No |
| M5 | Major | 7 of 11 error labels unpinned | Apply | **Applied** | No | `spill.rs`, `spill/tests.rs` | Pass | No |
| M6 | Major | `Io` stands for six operations | Apply | **Applied** | No | `spill.rs`, `spill/tests.rs` | Pass | No |
| M7 | Major | `Io` unreached by any test | Apply | **Applied** | No | `spill/tests.rs` | Pass | No |
| M8 | Major | a per-sample failure does not name the sample | Apply | **Applied** | No | `spill.rs`, `spill/tests.rs` | Pass | No |
| M9 | Major | no completeness witness | Apply (doc) + Defer (mechanism) | **Applied with adaptation** | No | `spill.rs` | Pass | **Yes — B2** |
| M10 | Major | the two flags admit a state the spec excludes | Apply (refusal) + Defer (sum type) | **Applied with adaptation** | No | `spill.rs`, `spill/tests.rs` | Pass | **Yes — owner** |
| M11 | Major | no property test | Apply | **Applied** | No | `spill/tests.rs` | Pass | No |
| Mi1 | Minor | "fixed-width or carries its own length" is false of varints | Apply | **Applied** | No | `spill.rs` | Pass | No |
| Mi2 | Minor | `finish`'s doc claims a guarantee it lacks | Apply | **Applied** | No | `spill.rs` | Pass | No |
| Mi3 | Minor | "spills once" is true of both sides | Apply | **Applied** | No | `mod.rs` | Pass | No |
| Mi4 | Minor | a doc cites a file not in the repository | Apply | **Applied** | No | `mod.rs` | Pass | No |
| Mi5 | Minor | the stand-in type has two public paths | Apply | **Applied** | No | `mod.rs`, `spill.rs`, `spill/tests.rs` | Pass | No |
| Mi6 | Minor | `SAMPLE_CAPACITY_HINT` names a ceiling as a hint | Apply | **Applied with adaptation** | No | `spill.rs` | Pass | No |
| Mi7 | Minor | the writer states no buffering requirement | Apply | **Applied** | No | `spill.rs` | Pass | No |
| Mi8 | Minor | a ten-byte varint past `u64` is absorbed | Apply | **Applied** | No | `spill.rs`, `spill/tests.rs` | Pass | No |
| Mi9 | Minor | the error labels are spelled three ways | Apply | **Applied** | No | `spill.rs`, `spill/tests.rs` | Pass | No |
| Mi10 | Minor | the reader allocates two `Vec`s per entry | Defer | **Deferred** | No | None | N/A | **Yes — D2** |
| Mi11 | Minor | `OutOfRange` on the read counts untested | Apply | **Applied** | No | `spill/tests.rs` | Pass | No |
| Mi12 | Minor | the boundary test misses the floats' signed edges | Apply | **Applied** | No | `spill/tests.rs` | Pass | No |
| Mi13 | Minor | two round trips cannot take a path the others do not | Apply | **Applied** | No | `spill/tests.rs` | Pass | No |
| Mi14 | Minor | `is_biallelic_snp`'s doc says the decision happens here | Apply | **Applied** | No | `spill.rs` | Pass | No |
| Mi15 | Minor | *verdict* and *cut* undefined; the two columns unnamed | Apply | **Applied** | No | `spill.rs` | Pass | No |
| Mi16 | Minor | two clauses assert importance instead of a fact | Apply | **Applied** | No | `spill.rs` | Pass | No |
| Mi17 | Minor | milestone labels in the module's status line | Apply | **Applied** | No | `mod.rs` | Pass | No |
| Mi18 | Minor | `same_bits` is not a predicate; `heads_match` reuses "head" | Apply | **Applied** | No | `spill/tests.rs` | Pass | No |
| Mi19 | Minor | `read` names a byte count where *read* means a sequencing read | Apply | **Applied** | No | `spill.rs` | Pass | No |
| Mi20 | Minor | `next_entry` duplicates `Iterator::next`, and has no `# Errors` | Apply (docs) + keep both | **Applied with adaptation** | No | `spill.rs` | Pass | No |
| Mi21 | Minor | `SpillWriter` carries an unmarked `finish` obligation | Apply | **Applied** | No | `spill.rs` | Pass | No |
| Nits | Nit | thirteen wording and naming nits | Apply | **11 Applied, 2 Deferred** | No | all three files | Pass | No |

## 3. Questions asked and answers

None — the two decisions the review raises are recorded for the milestone checkpoint rather than
asked mid-run, because neither blocks B2 or B3.

## 4. Per-finding log

Findings whose fix is a single doc-comment or name change are grouped at the end; each of them
was applied verbatim as the review proposed, with the module's own prose adjusted around it.

### B1 — the comparator cannot be shown to fail
- **Severity:** Blocker · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** `holds_the_same_bits_separates_entries_that_differ_in_any_one_field` asserts
  that an entry matches itself and then that the comparator separates eleven pairs differing in
  exactly one field each — contig, position, both flags, the line, the sample count, both floats,
  both read counts, and **an absent window against a zeroed one**, which is spec §6 trap 4 stated
  as a property of the comparator rather than of the codec.
- **Review suggestion used verbatim?** Adapted — the same eleven pairs, renamed for the renamed
  comparator.
- **Verification:** hard-wiring `holds_the_same_bits` to `return true` gives
  `FAILED. 28 passed; 1 failed`, the one being this test. Before the fix the same mutation left
  all 17 green. The mutation was reverted from a backup and the restoration confirmed by `diff`.
- **Files:** `spill/tests.rs` · **Follow-up:** None

### M1 — a corrupt line length pulls in the rest of the file
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** `MAX_LINE_BYTES = 16 MiB`, checked before the bytes are read;
  `read_line` is now its own function so the ceiling and the short-read check sit together.
- **Adaptation:** the review proposed a 1 MiB ceiling and a `vec![0u8; line_length]` +
  `read_exact`. **16 MiB, and `take(..).read_to_end(..)` kept.** A ceiling that is too tight turns
  a valid run into an error, which is worse than the exhaustion it prevents: a line is the eight
  fixed columns plus one genotype column per sample, so 16 MiB leaves room for five times spec
  §4's three thousand samples at a hundred bytes each. And `read_to_end` over a `take` grows as
  the bytes arrive, where `vec![0u8; line_length]` would reserve the full length up front — so
  keeping it means a corrupt length under the ceiling still costs only what the file holds.
- **Verification:** `a_line_longer_than_the_ceiling_is_refused_before_it_is_read` asserts
  `OutOfRange { field: "line_length", value: 4294967295 }` on a header claiming a
  four-billion-byte line followed by 4 KiB of payload.
- **Files:** `spill.rs`, `spill/tests.rs` · **Follow-up:** None

### M2 — the comparator reads its fields by access
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** both sides of `holds_the_same_bits` destructured with no `..`, and the
  per-sample half split into `sample_holds_the_same_bits`, which destructures `SpilledSample` and
  `WindowCoverage` too.
- **Verification:** the review's experiment re-run — `scoring_weight: u32` added to
  `SpilledSample`, encoded correctly, decoded as `0`, with the three fixtures given `7` and the
  pinned byte list extended, exactly as a coder resolving the compile errors would. Before the
  fix: `17 passed, 0 failed`. After: **two `E0027` inside the comparator**, so the tree does not
  build. Both files were restored from byte-compared backups.
- **Files:** `spill/tests.rs` · **Follow-up:** None

### M3 — a decode error is not terminal
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** `SpillReader` gains a `failed` latch set by the first failure, after which
  `next_entry` returns `None`; `impl FusedIterator for SpillReader`.
- **Verification:** `the_reader_stops_after_a_decode_error` feeds a three-byte corrupt prefix
  followed by a whole valid entry and asserts the second call is `None`. Removing the latch makes
  it fail; before the fix the second call returned a well-formed `SpillEntry` decoded from the
  middle of the file.
- **Files:** `spill.rs`, `spill/tests.rs` · **Follow-up:** None

### M4 — no fixture crosses the 128-byte varint boundary
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** `a_line_longer_than_one_varint_byte_round_trips` (a 300-byte line) and
  `a_cohort_of_a_thousand_samples_round_trips`.
- **Verification:** reading the line length as a single byte instead of a varint gives
  `FAILED. 26 passed; 3 failed` — both new tests and the property. Before the fix the same
  mutation left all 17 green.
- **Files:** `spill/tests.rs` · **Follow-up:** None

### M5 — seven of eleven error labels unpinned
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** `a_file_cut_anywhere_inside_a_record_names_the_field_the_bytes_ran_out_in`
  replaces the old walk, asserting the expected label at each of the nineteen offsets inside
  `a_tiny_record`'s twenty-byte encoding — so all eleven labels are pinned by one test, not
  eleven. `a_biallelic_flag_byte_that_is_neither_zero_nor_one_is_refused` covers the one flag no
  test corrupted, and `BIALLELIC_FLAG_OFFSET_IN_TINY_RECORD` sits beside its sibling so the gap
  cannot reopen. Beside it, a private `field` module holds the eleven labels, so a rename has one
  site.
- **Verification:** the nineteen offsets each produce the label the table names; the test asserts
  the encoding is twenty bytes first, so a layout change fails loudly rather than shifting the
  table under it.
- **Files:** `spill.rs`, `spill/tests.rs` · **Follow-up:** None

### M6 — `Io` stands for six operations
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** `#[from] io::Error` dropped; `Write`, `Flush` and `Read { field, source }`
  replace it, matching `VcfWriteError`'s shape next door.
- **Verification:** `append_reports_a_sink_that_refused_the_bytes_and_does_not_count_the_entry`
  asserts `Write`, `finish_reports_a_flush_that_failed` asserts `Flush`.
- **Files:** `spill.rs`, `spill/tests.rs` · **Follow-up:** None

### M7 — `Io` unreached by any test
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** a `RefusesEverything` sink whose `write` and `flush` both fail.
- **Verification:** the two tests above. Deleting `self.sink.flush()?` from `finish` used to pass
  all 17; it now fails `finish_reports_a_flush_that_failed`.
- **Files:** `spill/tests.rs` · **Follow-up:** None

### M8 — a per-sample failure does not name the sample
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** `SpillError::InSample { index, source: Box<SpillError> }` wraps whatever
  `decode_sample` returns. The record's ordinal stays the caller's, as spec §5 puts it.
- **Adaptation:** nested rather than adding an index to every variant — the alternative touches
  eleven construction sites to carry a number that is meaningless on nine of them.
- **Verification:** `a_failure_inside_the_samples_names_which_sample` cuts a three-sample entry two
  bytes into its **second** sample and asserts `index == 1`, with the offsets computed from the
  encoding rather than written down.
- **Files:** `spill.rs`, `spill/tests.rs` · **Follow-up:** None

### M9 — no completeness witness
- **Severity:** Major · **Initial decision:** Apply the doc, defer the mechanism ·
  **Final status:** Applied with adaptation
- **Implementation:** `finish`'s doc no longer claims that consuming the writer makes a forgotten
  flush visible. It now says what is true: the writer has no `Drop`, so a dropped writer leaves
  the sink to do whatever it does; **and a spill that lost its tail on an entry boundary reads
  back as a complete, shorter spill**, because nothing in the file says how many entries it should
  hold.
- **Deferred half:** threading `entries_written` to the reader so the shortfall is refused. B2
  builds the file's lifecycle and is the only place that can own the number; doing it here would
  give `SpillReader::new` a parameter no caller can supply yet.
- **Files:** `spill.rs` · **Follow-up:** **Yes — B2.**

### M10 — the two flags admit a state the spec excludes
- **Severity:** Major · **Initial decision:** Apply the refusal, defer the reshape ·
  **Final status:** Applied with adaptation
- **Implementation:** `SpillError::TractMarkedAsABiallelicSnp { contig, position }`, refused by
  `append` before anything is encoded and by `decode_entry` as soon as both flags are read. Two
  tests, one per side.
- **Deferred half:** the three-variant sum type two agents proposed. It changes spec §3.7's type
  block, which this loop may not do — raised as open question 1 for the checkpoint. The refusal is
  faithful to the spec as written (§3.2 puts repeat tracts inside "every other record") and needs
  no ruling.
- **Files:** `spill.rs`, `spill/tests.rs` · **Follow-up:** **Yes — owner.**

### M11 — no property test
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation:** `any_stream_of_entries_comes_back_bit_for_bit` over one to eight generated
  entries, with lines of 0–400 bytes and cohorts of 0–300 samples — both crossing the varint
  boundary — and floats drawn from `NaN`, a `NaN` with a payload, `±0.0`, both infinities, and any
  bit pattern at all. The generator never marks a tract as a biallelic SNP, since the writer now
  refuses that.
- **Verification:** it is one of the three tests the one-byte-line-length mutation kills.
  `proptest` was already a dev-dependency; no manifest change.
- **Note:** the mutation runs left a `proptest-regressions` seed behind — a failure case generated
  against a **deliberately broken** tree, not against shipped code. It was deleted rather than
  checked in: its "shrinks to" text describes a defect that never existed here, and the boundary it
  covers is already pinned by name in `a_line_longer_than_one_varint_byte_round_trips`.
- **Files:** `spill/tests.rs` · **Follow-up:** None

### Mi6 — `SAMPLE_CAPACITY_HINT` names a ceiling as a hint
- **Final status:** Applied with adaptation. Renamed `MAX_SAMPLES_RESERVED_UP_FRONT`, moved above
  its use site, and its doc now carries the numbers instead of pointing at a document: 16 bytes a
  sample, so 128 KiB at the cap; spec §4's largest cohort is three thousand; a corrupt count of
  4,026,531,840 would otherwise reserve about 64 GB. **Not made public** — the review suggested it
  so the policy could be reported at runtime, and nothing reports anything yet; a `pub` constant
  with no reader is surface without a caller.

### Mi10 — the reader allocates two `Vec`s per entry
- **Final status:** Deferred. The fix is a `read_into(&mut SpillEntry)` beside `next_entry`, and
  the smells agent wrote and compiled one. It is an optimisation with no measurement behind it and
  no caller to measure: plan step D2 reports the pass shares on the tomato slice, and that is the
  number that says whether two allocations a record matter. Adding the entry point now doubles the
  reader's surface for a cost nobody has weighed.

### Mi20 — `next_entry` duplicates `Iterator::next`
- **Final status:** Applied with adaptation. The `# Errors` section is added and the doc now says
  what the latch means. **Both names are kept**: `next_entry` is the named, documented operation
  with its own error contract, and the `Iterator` impl is what passes two and three will loop over.
  Removing either would cost the module the thing the other provides.

### The remaining Minors and Nits
Applied as proposed, all of them prose or naming: the layout sentence now says *self-delimiting*
(Mi1); the module header defines *verdict* and *cut* before using them and names `FILTER` and
`INFO` (Mi15); the two "this is important" clauses are gone (Mi16); `mod.rs` no longer claims a
contrast that is true of both sides (Mi3), cites no file that is not in the repository (Mi4), and
carries no milestone labels (Mi17); `WindowCoverage` is declared in `mod.rs` and no longer
re-exported from `spill` (Mi5); `SpillWriter`'s doc asks for a `BufWriter` (Mi7) and the type
carries `#[must_use]` (Mi21); `is_biallelic_snp`'s doc no longer says the decision happens here
(Mi14); a ten-byte varint above `u64::MAX` is refused (Mi8) with its own test; the eleven labels
are spelled once (Mi9); the read counts' `OutOfRange` has a test (Mi11); the boundary fixture
carries `-0.0` and both infinities (Mi12); the two line-content round trips became one loop that
says what they prove (Mi13); `read` became `bytes_read` (Mi19); `same_bits` became
`holds_the_same_bits` and its six-field local no longer calls itself *heads* (Mi18); and the
fixture helpers, the encoder helper, the offset constants and the byte-list test are renamed.

**Two nits were not applied.** The `NotABoolean` variant keeps its name while the tests call the
byte a *flag byte* — the variant names what the codec expected and the tests name what the field
is, and collapsing them would make one of the two worse. And the `as u64` on the encode side stays
where the decode side uses `try_from`: `usize` to `u64` is a lossless widening on every target
this builds for, and `try_from` there would add a `Result` that cannot be `Err`.

## 5. Deferred findings to carry forward
- **M9 (part)** — carry `entries_written` to the reader so a lost tail is refused. **B2.**
- **M10 (part)** — whether `SpilledSample` becomes a sum type. **The owner's, at the checkpoint.**
- **Mi10** — a reader-side scratch buffer. **After D2's pass-share measurement.**
- **Mi6 (part)** — making the reservation policy inspectable at runtime, when something reports.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — nothing under `src/ng/run/paralog_filter/` is reachable from any harness in
  `benches/`, and no run writes or reads a spill yet.
- **Outcome:** skipped.

## 10. Commands run
- `cargo test --all-features --lib "ng::run::paralog_filter::"` (repeatedly, including under four
  mutations)
- `cargo test --all-features --lib --bins --tests`
- `cargo clippy --lib --bins --tests --all-features -- -D warnings`
- `cargo fmt` then `git checkout --` on the nine files `main` owns
- `cargo fmt --check`

## 11. Command results
- `cargo test --all-features --lib "ng::run::paralog_filter::"` → 0, `ok. 29 passed; 0 failed;
  0 ignored; 6373 filtered out`
- `cargo test --all-features --lib --bins --tests` → 101, **6,387 lib tests passed, 0 failed, 15
  ignored**; one integration test failed,
  `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`, which is `main`'s at
  `a33ada0f`
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` → 101, 3 errors, all
  `needless_lifetimes` in `cohort_merge/build.rs` and `serial.rs`, `main`'s
- `cargo fmt --check` → 1, dirty on 9 files, none of them `src/ng/run/paralog_filter/`

The gate this plan uses — no failure `main` does not already have, and the `fmt` and `clippy`
failure sets no larger than `main`'s — holds. The lib count moves 6,375 → 6,387, and the module's
own tests 17 → 29.

## 12. Notes
- **`cargo fmt` reformats `main`'s nine files as a side effect**; they were reverted with
  `git checkout --` before staging, as at every step of this plan.
- **Three numbers in the implementation report and the commit message were wrong and are
  corrected**, all three the author's own count of the author's own tests. The corrections are
  folded into `2ad2a92b` along with these fixes.
- **The eight review worktrees were removed** and their branches deleted.
