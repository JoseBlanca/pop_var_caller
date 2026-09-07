# Fix Application Report: ng_paralog_filter_b2_2026-09-06.md

**Date:** 2026-09-06
**Source review:** [ng_paralog_filter_b2_2026-09-06.md](ng_paralog_filter_b2_2026-09-06.md)
**Source state reviewed against:** `f7cd840b` on `ng-paralog-filter`; fixes folded into that commit
**Execution mode:** non-interactive
**Overall status:** Completed

---

## The answer

**Nothing can write over the spill while its run holds it, and the run cannot write over it
itself.** A record appended after pass one has ended used to reopen the file with
`truncate(true)`: measured on the reviewed commit, a 228-byte three-record spill became 77 bytes
and one record, `entries_written` said 4, and the read then reported *"ends after 1 records, where
4 were written"* — the lost-tail message for a file that had lost its head. Three of the four
review agents found it independently. It is now refused, and the file is created with `create_new`
rather than `truncate`, so a leftover from a killed run or another run's live spill is **named**
instead of destroyed.

**The completeness check is no longer something a caller can forget.** It was a builder call on
top of `SpillReader::new`, whose default left it off; the count is now an argument, so the reader
that made a short file undetectable does not compile.

**Three of my own claims were wrong and are corrected here and in the commit message:**

| claim | right value |
|---|---|
| the guard-on-the-tidy-path mutation fails "4 of 12" | **3**, re-measured. The fourth failure was a stale file left by an earlier failing run, from before `an_output_path` learned to clear one — so the number came from the very defect that run had just exposed. |
| `spill.rs` "gains 54 lines" | **44 net** (49 added, 5 removed); 54 was the diffstat's churn |
| the lost-tail test cuts the file with `set_len` | **it cut nothing.** The test measured the two-record boundary by building a second `SpillFile` at the same path, and that one's first append had already truncated the file to exactly the length being measured. The assertion held; the step presented as its cause was inert. |

The last is the one worth keeping: it is a wrong *mechanism*, not a wrong number, and it was
found by two agents independently after the fix to the aliasing behaviour made the aliasing
visible.

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 6
- Minors: 18
- Nits: 5 (grouped)

### Outcome totals
- Applied: 24
- Applied with adaptation: 3
- Already fixed: 0
- Deferred: 2
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo test --all-features --lib "ng::run::paralog_filter::"` → 0, `ok. 51 passed; 0 failed`
- `cargo test --all-features --lib --bins --tests` → non-zero, **6,409 lib tests pass**; the one
  integration failure is `main`'s
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` → non-zero, 3 errors, all
  `needless_lifetimes` in `cohort_merge/`, `main`'s
- `cargo fmt --check` → non-zero, dirty on 9 files, none of them this module's, `main`'s
- `cargo doc --no-deps`, `cargo audit`, `--all-targets` → not run, for the reasons the review's §3
  gives
- Performance check → **skipped**: nothing under `src/ng/run/paralog_filter/` is reachable from any
  harness in `benches/`, and no run writes a spill yet.

### Unresolved high-priority findings
- **M6 (part)** — how much larger than a compressed output the spill is, measured. The claim it
  contradicted is gone; the number belongs to step D2, where a real run exists.
- **Mi17** — the spill's disk cost is not reported anywhere. Spec §3.5's run-report list has no
  line for it, so adding one is a spec question rather than a code fix; raised at Checkpoint B.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | an append after pass one empties the spill | Apply | **Applied** | `spill_file.rs`, `spill_file/tests.rs` | Mutation killed |
| M2 | Major | a file already at the path is truncated | Apply | **Applied** | `spill_file.rs`, `spill_file/tests.rs` | Mutation killed |
| M3 | Major | the completeness check is opt-in | Apply | **Applied** | `spill.rs`, both test files | Compiles only with the count |
| M4 | Major | the panic test silences other tests | Apply | **Applied** | `spill_file/tests.rs` | Pass |
| M5 | Major | a failed removal is silent | Apply | **Applied** | `spill_file.rs` | Pass |
| M6 | Major | the "room for the output" claim is false | Apply (prose) + Defer (the ratio) | **Applied with adaptation** | `spill_file.rs` | Pass |
| Mi1 | Minor | the lost-tail test's cut is inert | Apply | **Applied** | `spill_file/tests.rs` | Pass |
| Mi2 | Minor | `finish_writing`'s idempotence untested | Apply | **Applied** | `spill_file/tests.rs` | Pass |
| Mi3 | Minor | four error variants untested | Apply | **Applied with adaptation** | `spill_file/tests.rs` | Pass |
| Mi4 | Minor | `read`'s never-opened half untested | Apply | **Applied** | `spill_file/tests.rs` | Pass |
| Mi5 | Minor | "leaves nothing behind" cannot fail | Apply | **Applied** | `spill_file/tests.rs` | Pass |
| Mi6 | Minor | the completeness direction untested | Apply | **Applied** | `spill.rs`, `spill/tests.rs` | Pass |
| Mi7 | Minor | a refused record is counted | Apply | **Applied** | `spill_file/tests.rs` | Pass |
| Mi8 | Minor | the `expect` has no PANIC-FREE comment | Apply | **Applied** | `spill_file.rs` | Pass |
| Mi9 | Minor | the wrapped error's sentence is doubled | Apply | **Applied** | `spill_file.rs` | Pass |
| Mi10 | Minor | read failures carry no path | Apply | **Applied** | `spill_file.rs`, `spill_file/tests.rs` | Pass |
| Mi11 | Minor | neither error enum is `#[non_exhaustive]` | Apply | **Applied** | `spill.rs`, `spill_file.rs` | Pass |
| Mi12 | Minor | `StillWriting`, `exists`, `read` are misnamed | Apply | **Applied with adaptation** | `spill_file.rs` | Pass |
| Mi13 | Minor | `BUFFER_BYTES`'s "syscall per field" | Apply | **Applied** | `spill_file.rs` | Pass |
| Mi14 | Minor | "leaves no file at all" contradicts the code | Apply | **Applied** | `spill_file.rs` | Pass |
| Mi15 | Minor | one of three path classes tested | Apply | **Applied** | `spill_file/tests.rs` | Pass |
| Mi16 | Minor | `Drop` does not destructure `Self` | Apply | **Applied** | `spill_file.rs` | Pass |
| Mi17 | Minor | the disk cost is never reported | Defer | **Deferred** | None | N/A |
| Mi18 | Minor | the front door re-exports the reader | Apply | **Applied** | closed by M3 | Pass |
| Nits | Nit | five wording and shape nits | Apply | **4 Applied, 1 Deferred** | `spill_file.rs`, tests | Pass |
| §7a | Major (B1) | a corrupt sample count builds `Vec`s bounded by the file | Apply | **Applied** | `spill.rs`, `spill/tests.rs` | Pass |
| §7b | Minor (B1) | `read_f32_bits`' non-EOF arm untested | Apply | **Applied** | `spill/tests.rs` | Pass |

## 3. Questions asked and answers

None. Both open questions are measurements rather than decisions, and both are recorded for the
checkpoint.

## 4. Per-finding log

### M1 — an append after pass one empties the spill
- **Final status:** Applied. `SpillFile` gains a three-state `Stage` — named, being written,
  finished — so the file's position in its one life is a value rather than an inference from
  `writer.is_some()`. `append` in the finished stage returns `AppendedAfterPassOneEnded`;
  `finish_writing` moves to `Finished` *before* the flush can fail, so a failed flush cannot leave
  the spill looking as though pass one never opened it.
- **Verification:** restoring the reopen-and-truncate behaviour fails **3 of 51** —
  `appending_after_pass_one_has_ended_is_refused`,
  `something_already_at_the_path_is_refused_rather_than_overwritten` and
  `a_second_spill_for_one_output_is_refused_rather_than_clobbering_the_first`. The mutation was
  reverted from a byte-compared backup.

### M2 — a file already at the path is truncated
- **Final status:** Applied. `create_new(true)` replaces `create(true).truncate(true)`, and
  `AlreadyExists` becomes `SomethingIsAlreadyThere`, whose message tells the operator to remove the
  file or give the run its own output path.
- **Adaptation:** the review noted `tempfile` is production's answer and is ruled out by spec
  §3.4's fixed name. Naming the collision is the version that needs no spec change. The cost is
  that a spill left by a killed run blocks the next run of the same output until it is removed —
  which is the right way round: the alternative destroys it.
- **Verification:** two tests, one for a planted leftover and one for a second live `SpillFile`,
  each asserting the first file's bytes survive.

### M3 — the completeness check is opt-in
- **Final status:** Applied. `SpillReader::new(source, entries_expected)`; `expecting` is gone.
  Fourteen call sites in `spill/tests.rs` now pass the count they wrote.
- **Verification:** the omission the review measured no longer compiles.

### M4 — the panic test silences other tests
- **Final status:** Applied. The hook is left alone, with a comment saying why the quieter log is
  not worth it: the lib binary runs several thousand tests on many threads, and a test failing
  inside the window would report its name and no message.

### M5 — a failed removal is silent
- **Final status:** Applied. `Drop` warns to stderr naming the path and the reason, guarded by
  `!std::thread::panicking()`, following `reference_info.rs`'s `VerificationHandle`. It also
  destructures `Self` now (Mi16), so a field added to the type has to be accounted for in the one
  impl obliged to release every resource.

### M6 — the "room for the output" claim
- **Final status:** Applied with adaptation. The claim is replaced by what is true — the spill
  holds each line uncompressed plus about ten bytes a sample, so against a `.vcf.gz` output it is
  several times larger. **The ratio is not stated**, because it has not been measured; step D2 is
  the first run that could.

### Mi1 — the lost-tail test's cut is inert
- **Final status:** Applied. The two-record boundary is now computed in memory with a plain
  `SpillWriter` over a `Vec`, and the test asserts the boundary is shorter than the file before
  cutting — the assertion whose absence let the inert version look sound.
- **Note:** M2's fix would have turned the old version into a hard failure rather than a silent
  no-op, since the second `SpillFile` is now refused. Both were fixed together.

### Mi3 — four error variants untested
- **Final status:** Applied with adaptation. `Create` and `Read` now have tests. `Write` is
  reachable through a record the codec refuses, which is what `a_record_the_codec_refused_is_not_counted`
  drives. **`Flush` remains unreachable through `SpillFile`**, because the sink is a hard-wired
  `File`; B1's `RefusesEverything` covers the same code path one layer down, and making the sink
  injectable to reach it here would add a type parameter for a test.

### §7a — a corrupt sample count builds `Vec`s bounded by the file
- **Final status:** Applied. This is B1 code, found by B1's extras agent after B1's review had
  been written, and **fixed forward here** rather than by amending a settled step. The reservation
  cap bounded what was reserved and not what the loop built, so a corrupt count kept decoding
  until the file ran out. A `MAX_SAMPLES` ceiling of one million — three hundred times spec §4's
  largest cohort, 16 MB of `SpilledSample` at the ceiling — refuses the count itself, which is
  what `MAX_LINE_BYTES` already does for the other variable-length field.
- **Adaptation:** the existing test that drove a four-billion count now exceeds the ceiling, so it
  was moved to nine hundred thousand — under the ceiling, so it still exercises the reservation cap
  and the truncation — and the ceiling has its own test.

### The remaining Minors and Nits
Applied as proposed. `StillWriting` became `PassOneHasNotEnded`, which is what it means in both
the situations it covers; `exists` became `stage`, a three-state enum, which is what the field was
really tracking; `BUFFER_BYTES`'s doc now says one write per record and one read per byte, with
the size at which a record stops fitting the buffer; the module doc no longer claims a run that
calls nothing leaves no file; both error enums are `#[non_exhaustive]`; the doubled sentence in
`Flush`'s message is gone; and the `expect` in `append` carries its `// PANIC-FREE:` reason.

**Not applied.** `read` keeps its name rather than becoming `entries`: it is the counterpart of
`finish_writing` in the three-pass vocabulary the module's header sets up, and `entries` reads as
an accessor for something the type holds. Recorded rather than argued.

## 5. Deferred findings to carry forward
- **M6 (part)** — the spill's size against a compressed output, measured. **Step D2.**
- **Mi17** — reporting the spill's disk cost. Spec §3.5's run-report list has no line for it, so
  it is a spec question. **Checkpoint B.**

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — nothing in this module is reachable from any harness in `benches/`.
- **Outcome:** skipped.

## 10. Commands run
- `cargo test --all-features --lib "ng::run::paralog_filter::"` and
  `… "ng::run::paralog_filter::spill_file"`, repeatedly and under three mutations
- `cargo test --all-features --lib --bins --tests`
- `cargo clippy --lib --bins --tests --all-features -- -D warnings`
- `cargo fmt`, then `git checkout --` on the nine files `main` owns
- `cargo fmt --check`

## 11. Command results
- `cargo test --all-features --lib "ng::run::paralog_filter::"` → 0, `ok. 51 passed; 0 failed;
  0 ignored; 6373 filtered out`
- `cargo test --all-features --lib --bins --tests` → 101, **6,409 lib tests passed, 0 failed, 15
  ignored**; one integration test failed,
  `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`, `main`'s at `a33ada0f`
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` → 101, 3 errors, all
  `needless_lifetimes` in `cohort_merge/`, `main`'s
- `cargo fmt --check` → 1, dirty on 9 files, none of them `src/ng/run/paralog_filter/`

The gate holds. The lib count moves 6,399 → 6,409 and the module's own tests 41 → 51.

**Three mutations, each killed, each reverted from a byte-compared backup:**

| mutation | outcome |
|---|---|
| the `Drop` deleted | 4 of 19 spill-file tests fail — all four exits |
| the guard fired only where pass one was properly finished | **3** fail, and `the_file_is_gone_when_the_run_ends_normally` passes |
| an append after pass one reopening with `truncate` | 3 of 51 fail |

## 12. Notes
- **`cargo fmt` reformats `main`'s nine files as a side effect**; they were reverted with
  `git checkout --` before staging, as at every step of this plan.
- **Two findings are B1's, found late and fixed forward here** rather than by amending a settled
  step: the sample-count memory bound and the untested non-EOF read arm.
- **The four review worktrees were removed** and their branches deleted.
