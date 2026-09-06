# Code Review: ng_paralog_filter_b2

**Date:** 2026-09-06
**Reviewer:** rust-code-review skill (orchestrator), four sub-agents in isolated worktrees
**Scope:** step B2 of the hidden-duplication filter plan — the spill file and its guard
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** commit `f7cd840b` on branch `ng-paralog-filter` — two new files, plus
  44 net lines in `spill.rs` and 4 in `mod.rs`.
- **Reviewed against:** `f7cd840b`. Every sub-agent detached its own worktree to that commit and
  confirmed two branch-only files before starting.
- **In-scope files:**
  - [spill_file.rs](../../../../src/ng/run/paralog_filter/spill_file.rs)
  - [tests.rs](../../../../src/ng/run/paralog_filter/spill_file/tests.rs)
  - the changed parts of [spill.rs](../../../../src/ng/run/paralog_filter/spill.rs) — the reader's
    completeness check
  - [mod.rs](../../../../src/ng/run/paralog_filter/mod.rs)
  - [ng_paralog_filter_b2_2026-09-06.md](../implementations/ng_paralog_filter_b2_2026-09-06.md) — its
    quantitative claims
- **Deliberately out of scope:** the rest of `spill.rs` (reviewed at B1); `src/ng/paralog/`;
  `src/paralog/`, `src/var_calling/`, `src/sample_summary/` (production, frozen).
- **Categories dispatched:** reliability; errors + defaults + the diff's own numbers; naming +
  idiomatic + module_structure; extras (a filesystem resource guard on a whole-genome run's path)
  + refactor_safety. `tooling` skipped — `Cargo.toml` untouched; `unsafe_concurrency` skipped —
  no `unsafe`, no threads, no shared state.

### 2. Verdict

**Request-changes.** No Blocker. Six Major findings, of which one is convergent across three of
the four agents and is the one that costs results: **an append after pass one has ended silently
empties the whole spill.**

### 3. Execution status

| command | result |
|---|---|
| `cargo test --all-features --lib "ng::run::paralog_filter::"` | `ok. 41 passed; 0 failed; 6373 filtered out` |
| `cargo test --all-features --lib --bins --tests` | 6,399 lib tests passed; one integration test failed, pre-existing on `main` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 3 errors, all `needless_lifetimes` in `cohort_merge/`, pre-existing |
| `cargo fmt --check` | dirty on 9 files, none of them this commit's, pre-existing |

`--all-targets`, `cargo doc --no-deps` and `cargo audit` were not run, for B1's reasons: an
example has not compiled on `main` since `a33ada0f`, the doc build is red on other modules' links,
and the container has no offline advisory database.

**Findings labelled "Needs verification": 0.** Every finding below was produced by mutation or by
a probe run in the reviewer's own worktree, and the outputs are quoted in the per-category files
under [tmp/review_2026-09-06_ng_paralog_filter_b2/](../../../../tmp/review_2026-09-06_ng_paralog_filter_b2/).

### 4. Open questions and assumptions

1. **What should happen when something is already at the spill's path?** Spec §3.4 fixes the name
   `<output>.paralog-spill.tmp`, so two runs writing one output collide on it, and a run killed
   from outside leaves one behind. Production sidesteps this with `tempfile`, which the fixed name
   rules out. Affects **M1**, **M2**.
2. **How much larger than the output is the spill?** The module's doc claims a filesystem with
   room for the output has room for it; the spill is uncompressed and the output may be bgzf. The
   direction is certain from the code, the ratio is one measurement, and step D2 is where a real
   run exists to measure it on. Affects **M6**.

### 5. Top 3 priorities

1. **M1** — an append after pass one reopens with `truncate(true)` and destroys the spill, returning `Ok`.
2. **M3** — the completeness check is opt-in, so forgetting one call restores exactly the failure this step was built to close.
3. **M4** — the panic test replaces the process-wide panic hook, so any other test failing in that window reports no message at all.

### 6. Findings

#### Major

- **M1: [spill_file.rs:159-174](../../../../src/ng/run/paralog_filter/spill_file.rs#L159-L174) — an append after pass one has ended empties the whole spill, and returns `Ok`**
- **Categories:** idiomatic, errors, refactor_safety, reliability — convergent across all four
- **Confidence:** High — measured three times independently.
- **Problem:** `open_for_writing` uses `.truncate(true)` and is reachable again after
  `finish_writing`, because nothing records that pass one has ended. Measured: a 228-byte,
  three-record file became 77 bytes and one record; `entries_written` reported 4; the read then
  said "ends after 1 records, where 4 were written" — the *lost tail* message for a file that lost
  its head. The same mechanism makes two `SpillFile`s for one output destroy each other's file,
  against the type's own doc ("One of these exists for one run"), and lets a `read` cursor be
  invalidated under a caller's feet.
- **Why it matters:** the completeness check this step added catches the result, so no wrong VCF
  is written — but the records are already gone, and the message sends the reader hunting a
  truncation that did not happen. C3 and C4 hold the two longest-lived cursors this module will
  have.
- **Fix:** record which stage the file is in, refuse the append, and create with `create_new`
  rather than `truncate` so a collision is named instead of overwritten.

- **M2: [spill_file.rs:159-174](../../../../src/ng/run/paralog_filter/spill_file.rs#L159-L174) — a file already at the path is truncated without a word**
- **Categories:** extras, errors
- **Confidence:** High — measured.
- **Problem:** whatever is at `<output>.paralog-spill.tmp` when the run starts is destroyed by the
  first record. It can be two things and both matter: a spill left by a run that was killed, which
  is evidence about that failure, or a spill another run for the same output is writing now, whose
  records are then half gone. Measured: the second `SpillFile`'s first `append` destroyed 151 of
  228 bytes and returned `Ok(())`; dropping the first then unlinked the file the second was using,
  so the second's next `read` failed with `NotFound`.
- **Fix:** `create_new(true)`, with the `AlreadyExists` case named as its own error telling the
  operator what to remove. Open question 1.

- **M3: [spill.rs:378](../../../../src/ng/run/paralog_filter/spill.rs#L378) — the completeness check is opt-in, and forgetting it is silent**
- **Categories:** defaults, refactor_safety, extras — convergent
- **Confidence:** High — measured.
- **Problem:** `SpillReader::new` leaves `entries_expected` at `None` and `expecting` is a
  separate call. Measured: a bare `SpillReader::new` over a truncated file yielded 2 entries and
  no error, compiling clean; the same bytes through `SpillFile::read` gave the refusal. `#[must_use]`
  guards a discarded return, not an absent call, and `SpillReader` is re-exported at the module's
  front door — so the door whose default leaves the check off is the one a later step reaches for.
- **Why it matters:** this is the failure the step exists to close, one forgotten call away from
  coming back. C3 and C4 are the next callers.
- **Fix:** make the count an argument of `new`. Forgetting an argument does not compile.

- **M4: [tests.rs:156-163](../../../../src/ng/run/paralog_filter/spill_file/tests.rs#L156-L163) — the panic test replaces the process-wide panic hook, silencing every other test that fails in that window**
- **Categories:** reliability
- **Confidence:** High — measured.
- **Problem:** the test swaps the global hook for a discarding one across its `catch_unwind`, and
  the lib binary runs 6,399 tests on many threads. Two probe tests failing on the identical
  assertion — one inside the window, one outside — produced a report where the outside one carries
  its assertion text, values and `file:line`, and **the inside one has no stdout section at all**,
  just its name in the failures list. `catch_unwind` alone stops the unwind; the hook only buys a
  quieter log.
- **Fix:** leave the hook alone and accept one backtrace in the log.

- **M5: [spill_file.rs:177-190](../../../../src/ng/run/paralog_filter/spill_file.rs#L177-L190) — a failed removal is silent, where the crate warns**
- **Categories:** errors
- **Confidence:** High.
- **Problem:** `Drop` discards `fs::remove_file`'s error. The thing leaked is a file the spec
  sizes at about the VCF's, and the crate already answers this the other way:
  `reference_info.rs:1075-1081` warns to stderr from a `Drop`, guarded by
  `!std::thread::panicking()`.
- **Fix:** the same shape, guarded the same way.

- **M6: [spill_file.rs:44-47](../../../../src/ng/run/paralog_filter/spill_file.rs#L44-L47) — "a filesystem with room for the output has room for the spill" is not true for a compressed output**
- **Categories:** defaults
- **Confidence:** High on the mechanism, Medium on the size.
- **Problem:** ng compresses `.vcf.gz` / `.vcf.bgz` output; the spill is the same lines
  uncompressed plus about ten bytes a sample a record. The claim is the justification given for
  putting the file beside the output rather than anywhere else, and it is the one an operator
  would act on when sizing a filesystem.
- **Fix:** say what the spill actually holds and drop the claim. The ratio is open question 2.

#### Minor

- **Mi1: [tests.rs:297-320](../../../../src/ng/run/paralog_filter/spill_file/tests.rs#L297-L320)** — the
  lost-tail test's `set_len` shortens the file by **0 bytes**. It measures the two-record boundary
  by building a second `SpillFile` at the same path, and that one's first `append` has already
  truncated the file being measured. The assertion is sound and the intended mutation still dies,
  but by an aliasing writer rather than by the lost tail the comment describes. (naming,
  reliability, extras — convergent)
- **Mi2: [spill_file.rs:107](../../../../src/ng/run/paralog_filter/spill_file.rs#L107)** —
  `finish_writing` is documented "idempotent" and nothing tests it; deleting the guard that makes
  it so leaves all 41 green, and the mutant is not harmless: the second call re-enters
  `open_for_writing`, whose `truncate(true)` empties the spill. (reliability)
- **Mi3: [spill_file.rs:199-252](../../../../src/ng/run/paralog_filter/spill_file.rs#L199-L252)** — four
  of the five `SpillFileError` variants have no test. `Create` is reachable with no scaffolding;
  `Write` and `Flush` cannot be reached through `SpillFile` at all, because the sink is a
  hard-wired `File`. (reliability)
- **Mi4: [spill_file.rs:143-148](../../../../src/ng/run/paralog_filter/spill_file.rs#L143-L148)** — the
  `!self.exists` half of `read`'s guard has no test, and dropping it survives — so the "pass one
  has not finished" against "the file is missing" distinction that `StillWriting`'s own doc gives
  as its reason to exist is unpinned. (extras, reliability)
- **Mi5: [tests.rs:187-198](../../../../src/ng/run/paralog_filter/spill_file/tests.rs#L187-L198)** —
  `a_spill_nothing_was_written_to_leaves_nothing_behind` cannot fail for the property its message
  states: nothing ever creates a file for it to leave. (reliability)
- **Mi6: [spill.rs:407](../../../../src/ng/run/paralog_filter/spill.rs#L407)** — rewriting the
  completeness comparison from `!=` to `>` survives all 41 tests: a file holding *more* records
  than the writer counted reads back as complete. The variant is also named and worded for one
  direction only. (extras)
- **Mi7: [spill_file.rs:87-103](../../../../src/ng/run/paralog_filter/spill_file.rs#L87-L103)** — an
  entry the writer refused is counted: no spill-file test makes `append` fail, and the mutation
  that counts a refused record survives. (extras)
- **Mi8: [spill_file.rs:95-97](../../../../src/ng/run/paralog_filter/spill_file.rs#L95-L97)** — the
  `expect` in `append` cannot fire today, but carries no `// PANIC-FREE:` comment. (errors)
- **Mi9: [spill_file.rs](../../../../src/ng/run/paralog_filter/spill_file.rs)** — a wrapped
  `SpillError` renders the parent's sentence twice: measured, *"…could not be flushed to disk: the
  paralog spill could not be flushed to disk: No space left on device"*. And `SpillFileError::Write`
  carries a record the codec refused as though the filesystem had refused it. (errors)
- **Mi10:** every *read* failure escapes without the path — `read` hands back a `SpillReader` whose
  errors are bare `SpillError`, which carries no path by design, where spec §5 asks the file to be
  named on reads as well as writes. (errors)
- **Mi11:** neither error enum is `#[non_exhaustive]`, where `RunError` is. (errors,
  refactor_safety)
- **Mi12:** `read` names an operation that opens a cursor, mirroring the private
  `open_for_writing`; `StillWriting` names a situation where nothing is being written; `exists`'s
  doc says "whether the file exists on disk" where the field is really how far through its life
  the file is. (naming)
- **Mi13:** `BUFFER_BYTES`'s doc says "a syscall per field"; it is one write per *record* on the
  writing side and one read per *byte* on the reading side, which `spill.rs:343-345` already says
  correctly. (naming)
- **Mi14:** the module doc says a run that calls nothing "leaves no file at all", which
  `finish_writing` contradicts by creating one. (naming)
- **Mi15:** `spill_path_beside`'s doc names three input classes and one is tested. (reliability)
- **Mi16:** `Drop` names two fields by hand rather than destructuring `Self` — the one impl obliged
  to account for every resource-holding field. (refactor_safety)
- **Mi17:** the spill's disk cost is never reported: `entries_written` is a count, nothing exposes
  bytes, and spec §3.5's run-report list has no line for it. (extras)
- **Mi18:** `mod.rs` re-exports `SpillReader` and `SpillWriter` at the front door — the door whose
  default leaves the completeness check off. (module_structure)

#### Nits

`read` could be `entries`; `ShortFile` has no singular form in its message; `finish_writing` on an
empty run opens a file only to close it; the per-test scratch directories under `tmp/` are never
removed; and a `finish_writing` whose create failed is later reported as one that was never
called.

### 7. Out of scope observations

- **From B1's extras agent, which finished after B1's review was written:** a corrupt sample count
  builds `Vec`s bounded by the file rather than by one entry — a corrupt count over a 16 MiB file
  built 26.8 MB — because the reservation cap bounds what is *reserved* and not what the loop
  builds. And no test reaches `read_f32_bits`' non-EOF arm, so a device failure mid-float is
  reported as a file that ended. Both are B1 code found during B2's window and are fixed forward
  with B2's fixes.
- **The `cargo mutants` baseline fails out of the box** because `examples/ng_cohort_merge_real_cost.rs`
  does not compile on `main` at `a33ada0f`; the workaround is
  `additional_cargo_args = ["--all-features", "--lib"]`.

### 8. Missing tests to add now

**The file's life:** `appending_after_pass_one_has_ended_is_refused` (M1);
`something_already_at_the_path_is_refused_rather_than_overwritten` and
`a_second_spill_for_one_output_is_refused_rather_than_clobbering_the_first` (M2);
`finishing_pass_one_twice_leaves_the_records_where_they_are` (Mi2);
`a_record_the_codec_refused_is_not_counted` (Mi7);
`reading_before_pass_one_has_ended_is_refused` extended to the never-opened half (Mi4);
`a_spill_nothing_was_written_to_creates_nothing_and_removes_nothing` rewritten so it can fail (Mi5);
`a_spill_that_cannot_be_created_names_the_path_and_says_what_the_filesystem_said` (Mi3);
`a_failure_off_the_reader_can_be_given_the_file_it_happened_on` (Mi10);
`the_path_is_the_whole_output_path_with_a_suffix_added` over three input classes (Mi15).

**The codec:** `a_file_holding_more_records_than_were_written_is_refused` (Mi6);
`a_sample_count_past_the_ceiling_is_refused_before_the_samples_are_built` and
`a_source_that_fails_mid_float_says_the_read_failed_and_not_that_the_file_ended` (§7).

### 9. What's good

- **The three exits are tested separately and the separation earns its keep**: a guard that fires
  only where pass one was properly finished passes the normal-exit test and fails the other three.
- **`exists` as a field rather than a `stat` is the right call**, and for a reason the reviewers
  had to work out: a leftover from a killed run must not be unlinked, and `append` runs once per
  record. Only its name and doc were wrong.
- **The absence of an `fsync` is right**: both reads are `read(2)` in the same process served from
  the page cache, and the file is unlinked minutes later. The VCF writer's `fsync` guards a rename
  that publishes something meant to survive.
- **The `spill.rs` / `spill_file.rs` seam is at the right place** — `spill.rs` imports no
  `std::fs`, `spill_file.rs` names no field of the layout — and the `WindowCoverage` stand-in is
  placed so the eventual swap is three `use` lines.
- **Every filesystem failure already names the path and the OS reason** through the crate's error
  chain rendering: *"could not be created: Permission denied (os error 13)"*.

### 10. Commands to re-verify

- `cargo test --all-features --lib "ng::run::paralog_filter::"`
- `cargo test --all-features --lib --bins --tests`
- `cargo clippy --lib --bins --tests --all-features -- -D warnings`
- `cargo fmt --check` — and `git checkout --` the nine main-owned files afterwards
