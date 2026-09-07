# Code Review: ng_paralog_filter_c2

**Date:** 2026-09-07
**Reviewer:** rust-code-review skill (orchestrator), two sub-agents in isolated worktrees
**Scope:** step C2 of the hidden-duplication filter plan — pass one's sink and the flag that chooses it
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** commit `3b55b94c` on branch `ng-paralog-filter` — one new module and its
  tests, plus two CLI flags, an error variant and a guard in each of the two subcommands.
- **In-scope files:** [pass_one.rs](../../../../src/ng/run/paralog_filter/pass_one.rs),
  [tests.rs](../../../../src/ng/run/paralog_filter/pass_one/tests.rs),
  [mod.rs](../../../../src/ng/run/paralog_filter/mod.rs), the changed parts of
  [call_from_psps.rs](../../../../src/pop_var_caller_exp/call_from_psps.rs) and
  [call_from_alignments.rs](../../../../src/pop_var_caller_exp/call_from_alignments.rs), and the
  step's implementation report.
- **Out of scope:** steps A1–C1 (reviewed); `src/ng/window_coverage/` (reviewed on its own branch);
  `src/ng/paralog/*.rs` (guarded copies — findings raised, never edited); production.
- **Categories dispatched:** reliability + extras; errors + naming + idiomatic + defaults +
  refactor_safety + module_structure + smells + the diff's own numbers. `tooling` skipped —
  `Cargo.toml` untouched. `unsafe_concurrency` skipped — no `unsafe`, no threads.

### 2. Verdict

**Request-changes.** One Blocker, seven Major, twelve Minor.

**The theme is a suite that could not see its own fixtures' uniformity.** Sixteen mutations were
run against `entry_for`'s eight tests and **five survived** — the contig, the ploidy, the
alternative-count rule, the zip's length rule, and one whole CLI guard. Every survivor exploits the
same thing: every fixture in the file sat on contig 0, carried no no-call, and had no reads that no
written allele explains. The code was right in each case; the fixtures could not tell.

**And the step's own gate report was wrong.** It added two `clippy -D warnings` errors on lines it
authored and reported the gate as unchanged.

### 3. Execution status

| command | at `3b55b94c` | at the merge base `d533cc23` |
|---|---|---|
| `cargo test --all-features --lib --bins --tests` | lib `ok. 6574 passed` | `ok. 6565 passed` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | **11 errors** | 9 errors |
| `cargo fmt --check` | dirty on 9 files, none this step's | same 9 |

One integration test fails at both (`a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`),
which is `main`'s at `a33ada0f`. `--all-targets`, `cargo doc` and `cargo audit` not run, for B1's
reasons.

**Findings labelled "Needs verification": 0.** Every finding came from a mutation or a probe run in
a reviewer's own worktree.

**Mutation totals: 16 run, 5 survived, 0 changed no behaviour.**

### 4. Open questions and assumptions

1. **Is the `0.0` default acceptable where spec §3.6 says `0.01`?** Both agents were told it is a
   deliberate deviation and asked to judge the reasoning rather than the departure. Neither
   objected to the reasoning; one filed the *literal* rather than the value. Resolved: keep the
   value, name the constant.

### 5. Top 3 priorities

1. **B1** — the alternative-count rule is untested against the reads no written allele explains,
   and the wrong rule survives every fixture.
2. **M1** — direct mode's guard has no test at all; deleting it leaves the run writing an
   unfiltered VCF while ignoring what it was asked for.
3. **M2** — this step adds two clippy errors and reports the gate as clean.

### 6. Findings

#### Blocker

- **B1: [pass_one.rs:90](../../../../src/ng/run/paralog_filter/pass_one.rs#L90) — the alternative count is untested against reads no written allele explains, and the wrong rule passes everything**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** `alt_reads` is `Σ AD[1..]`. Deriving it as `DP − AD[0]` instead **survives all eight
  tests**, because every fixture builds its samples with `SampleReadCounts::new(counts, 0)` — no
  reads that no written allele explains. On real reads that count is routinely non-zero: reads that
  reached the locus and supported nothing the record lists.
- **Why it matters:** the wrong rule folds those reads into the non-reference share the scorer
  sees, at every locus of every sample, and nothing panics. The filter's whole allele signal is
  "what share of the reads is non-reference"; inflating it systematically biases every score.
- **Suggested fix:** a fixture with a non-zero unexplained count — 4 reference, 5 alternative, 6
  explained by neither, and the answer must be 5.

#### Major

- **M1: [call_from_alignments.rs](../../../../src/pop_var_caller_exp/call_from_alignments.rs) — direct mode's guard has no test, and deleting it leaves the run silently ignoring the flag.** psp mode got a refusal test; direct mode did not. With the guard removed, all 171 `pop_var_caller_exp` tests stay green and a run asking for `--paralog-fdr 0.01` writes an ordinary unfiltered VCF. Two modes, one rule, and only one of them was pinned.
- **M2: the step adds two `clippy -D warnings` errors and reports the gate as unchanged.** `variable does not need to be mutable`, on the `let mut writer` lines this step changed — the writer is now moved into the sink rather than mutated in the closure. 11 against the merge base's 9. **The orchestrator's own gate command was the cause**: it grepped for the three lint kinds the baseline already had, so a new kind was invisible by construction.
- **M3: [pass_one.rs](../../../../src/ng/run/paralog_filter/pass_one.rs) — the contig is untested.** Hardcoding `ContigId(0)` passes all eight tests; every fixture is on contig 0. A spill whose entries all claim one contig orders the whole file against the wrong sequence at pass three.
- **M4: the ploidy the run was given is untested.** Hardcoding `2` passes; ploidy shows only in how a *no-call* is spelled, and no fixture holds one.
- **M5: the zip's unequal-length rule is untested.** Padding instead of truncating passes. Judged in full under §7 — it is *short*, never mis-paired.
- **M6: `finish_parking` and both `PassOneError` variants have no test and no caller.** A sink that never flushed would leave passes two and three reading a file shorter than the run wrote, reported by the spill's completeness check as a truncation rather than as a forgotten flush.
- **M7: both subcommands create the VCF before choosing the sink**, so `pass_one.rs`'s central promise — with the filter on the VCF is never opened — is a property of the type the wiring cannot honour. Related: the writer was recovered from the sink by an `unreachable!` whose soundness rested on a guard ninety lines earlier, which becomes a post-calling-pass panic the moment C4 moves that guard.

#### Minor

- **Mi1: `--paralog-fdr` is an unvalidated `f64`.** `7`, `inf` and `nan` all parse; `-0.0` slips the `!= 0.0` guard.
- **Mi2: the `0.0` default is an unnamed literal written twice**, where every neighbouring flag uses a named constant.
- **Mi3: the guard sits ~90 lines below the function's own "everything a person typed is judged before a byte is read" block** — and by the time it runs, the reference is read and verified and every psp is open. The report's "before anything is opened" was true only of the output.
- **Mi4: [mod.rs](../../../../src/ng/run/paralog_filter/mod.rs) still says "Nothing fills a spill yet: the sink … [is a] later step"**, in the commit that adds the sink.
- **Mi5: [spill.rs](../../../../src/ng/run/paralog_filter/spill.rs) still documents the row shape as the record's *span***, which is the rule the owner replaced with the tract flag on 2026-09-07.
- **Mi6: `PassOneError::Vcf`'s message duplicates its parent `RunError::RecordNotWritten`'s**, so a failed write prints the sentence twice.
- **Mi7: `is_repeat_tract()` is read twice**, so the flag the entry carries and the rows it carries come from two answers to one question.
- **Mi8: a third `SpilledSamples` variant would not be caught in `entry_for`**, which branches on a bool — while the spill *writer* deliberately guards against exactly that.
- **Mi9: `entry_for` and `SpillingSink` are `pub` with no caller outside their module**, which is what keeps `dead_code` quiet.
- **Mi10: ~40 lines are copy-pasted across the two subcommands.**
- **Mi11: `--paralog-filter-tag` is parsed and unread**, with help text implying a non-zero target would activate it.
- **Mi12: the report's "6580 filtered out" does not reconcile** with 6,574 passing and 15 ignored; the real figure is 6,581.

#### Nits

- "no branch beyond the one that chose the sink" — the match runs per record.
- "three lines of the writing path changed in each subcommand" — three *edits*, six lines.

### 7. The zip, judged

**Short, not wrong.** Both slices are prefixes of the same sample order, so the rows the zip emits
pair the right window with the right column; a mismatch loses a suffix and never mis-pairs. At the
one production call site the lengths cannot differ — `window_coverage` and the record's columns
both come from `evidence.samples`, under an `assert_eq!` in `assemble_record`.

**C3's cohort-size check catches it only on one branch.** Measured, with three windows and one
sample column: a generic locus parks 1 row (caught) and a repeat tract parks 3 (not caught, and
correct, since tract rows never consult the columns). The cost of leaving it to C3 is that a
wiring regression is named after the whole calling pass has finished, by an error pointing at pass
two rather than at the sink.

### 8. The cost, measured

At 63 samples over 20,000 records in the container: `record_line` alone 209.8/219.1 ms, `entry_for`
261.5/274.2 ms — **13.1–13.7 µs a record against 10.5–11.0 µs for the encode it contains**, about a
quarter on top of work the *off* path already pays inside `write_record`. One allocation per record
(the rows `Vec`, 1,008 bytes at 63 samples, from a size-hint-exact iterator); `into_bytes()` reuses
the encoder's buffer and the spill writer reuses its scratch. **No finding** — the code already
takes the two cheap options.

### 9. What's good

- **The byte-identity claim is pinned by the suite, not only by reading the diff.** Making
  `accept`'s off arm a no-op fails three existing tests, including both mode-equivalence
  VCF-equality tests.
- **The written position is right and singly sourced**: `written_position` feeds the `POS` column,
  `place_of` and the entry, so the three cannot drift.
- **`into_bytes()` reuses the encoder's buffer** rather than copying the line.
- **The spill writer's shape check is a `match`**, which is what makes M8's gap in `entry_for`
  visible as an inconsistency rather than invisible.

### 10. Commands to re-verify

    ./scripts/dev.sh cargo test --all-features --lib "ng::run::paralog_filter"
    ./scripts/dev.sh cargo test --all-features --lib --bins --tests
    ./scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings 2>&1 | grep -E "^error: " | grep -v "could not compile" | sort | uniq -c
    ./scripts/dev.sh cargo fmt --check
