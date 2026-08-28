# Code Review: ng_parameters_file_b3
**Date:** 2026-08-28
**Reviewer:** rust-code-review skill (orchestrator), two agents in isolated worktrees
**Scope:** step B3 of the parameters-file plan — the provenance comments
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the uncommitted diff of step B3 over `83934abe`, applied as a patch in
  each agent's own worktree.
- **In-scope files:** [to_toml.rs](../../../../src/ng/calling/parameters_file/to_toml.rs) — the
  comment machinery, the `origins` module, a note on every section, and the new tests; and
  [testdata/every_shape_as_written.toml](../../../../src/ng/calling/parameters_file/testdata/every_shape_as_written.toml),
  regenerated.
- **Out of scope:** the writer's layout and escaping (B2, committed); the projection (B1); the
  reader (C).
- **Two agents, and only one of them reviews Rust.** A **reliability** pass over the comment
  machinery, and **the file's reader** — a geneticist who has not seen the code, re-running the
  exercise B2's review ran on the uncommented file, to ask whether the comments answer the
  questions it had. For a step whose whole output is prose addressed to that reader, that is the
  review that matters.

## 2. Verdict

**Request-changes.** No Blockers. Four Majors, all of them the same kind and none of them a Rust
defect: **the comments say things that are not true.** One of them answers the previous reader's
exact question with the wrong option.

## 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::calling::parameters_file` | 0 | 62 passed, 0 failed, 2 ignored |
| `cargo test --lib` | 0 | 4,985 passed, 0 failed, 13 ignored |

**Mutation numbers:** 34 run, **5 survived**, 4 of them changing no behaviour any fixture reaches,
1 a true survivor.

**One measurement about the suite's shape**, because it governs how the table reads:
`the_whole_shape_writes_the_documented_toml` was the only oracle for five mutations. It is a
golden-file comparison whose failure says the output *moved*, not that it is *wrong*, and its own
message invites regeneration — so for those five the last check is a human reading a diff.

## 4. Open questions

None new. **B3 confirmed one that was already open:** spec §8's headline comment — the slippage
defaults' "which alignments, at what depth, on which date" — has nowhere to attach, because a
slippage number carries a smoothing origin and no warrant. Same gap as B1's design Major and the
item standing since Checkpoint A; one ruling closes all three.

## 5. Top 3 priorities

1. **M1** — the contamination note names the wrong substance and the wrong grain.
2. **M2** — "each value carries a `warrant`" is false for most of the file's numbers, and the
   file's one editing rule rests on it.
3. **M5** — the flat-concentration note fires whatever the warrant says, and no test saw it.

## 6. Findings

### Major

**M1: the contamination note names the wrong substance and the wrong grain** — *the reader's test*.
It opened "how much of each library's **DNA** came from somebody else" and closed "to stop
correcting one **library**". The quantity is the share of that **read group's reads**
(`mod.rs:514`), and the grain is the lane — spec §3.4 makes it the whole point, with the reason
that "index hopping happens on a flowcell and not in a tube". Index hopping mislabels *reads*; it
does not put another individual's DNA in the tube. **So "DNA" is not a loose gloss: it names the one
mechanism the spec explicitly excludes.** A reader following the note thinks in libraries, edits one
row, and believes a four-lane library is uncorrected when three of its lanes still are — and the
fixture gives no signal, because every library in it has exactly one read group.

**M2: "each value carries a `warrant`" is false** — *the reader's test*. Eight numbers in the
produced file carry one. Every slippage number, both of the prior's concentrations, both length
spectra and the contamination fraction do not. The file's **one editing rule** — change the warrant
to `supplied` and delete the observations — was built on top of it, so a reader editing `level`
looks for the promised warrant, finds none, and either invents a key (which `deny_unknown_fields`
refuses) or concludes the rule did not apply.

**M3: `rung` is defined once and means two things** — *the reader's test*. The prior's ladder gets
an authoritative four-way gloss; a `ShareCurveRung` twenty lines below means what the *curve* was
fitted on. The crate already knows they are different questions (`mod.rs:940`). Giving one of the
two a confident definition and leaving the other bare makes the second **more** likely to be
misread, not less.

**M4: the batching note's justification is refuted by the rows beneath it** — *the reader's test*.
"The rows below look the same either way" — the rows below are batches 0, 1, 1. The true and
narrower claim is that a *declared* batching that happens to have one batch is indistinguishable
from an assumed one. A reader who checks a comment against the data three lines down and finds it
false stops trusting the comments that are right.

**M5: the flat-concentration note fires whatever the warrant says** — *reliability, the one true
survivor of 34 mutations*. `a_run_that_defaulted_nothing_writes_no_per_row_notes` listed two of the
five origin texts, and this was not one of them.

### Minor

- **Mi1** *reliability.* **Two origins are never exercised**, because no fixture defaults a
  substitution rate or an inbreeding coefficient — so a text put beside the wrong quantity is
  invisible there. The five origins are five interchangeable `&str`.
- **Mi2** *reliability.* **The wrapper's boundary and its unit are both unpinned**: moving `>` to
  `>=`, and counting bytes rather than characters, each survived everything. No note happens to sit
  at 78 characters or to carry a byte that is not a character.
- **Mi3** *reliability.* Swapping the two defaulted scalars' origins at their call sites was caught
  by the golden file alone — the targeted test asserted each text was *somewhere* in the document
  rather than above its own key.
- **Mi4** *the reader's test.* `[repeat_tracts]` spends its lines on where the slippage numbers
  came from and never says what `level`, `shorter_share`, `fall_off` or `slipped_reads` **are** —
  the four things a reader would change, and where both of the previous reader's unanswered
  guesses live. `slipped_reads = 8000.5` reads as a typo until somebody says it is an expected
  count.
- **Mi5** *the reader's test.* `origins::FLAT_CONCENTRATION` said "one chromosome's worth" beside a
  value of 1.25. The phrase glosses the constant 1.0; the upstream field doc gets it right with
  "**this many** chromosomes' worth".
- **Mi6** *the reader's test.* 67 comment lines was the right total, misallocated: seven to cut
  (design rationale a reader cannot act on), five to spend on the definitions above.

### What the reviews confirmed rather than faulted

- **The step's stated hazard is genuinely guarded.** A newline inside a note — a comment landing
  inside a row — failed **seven** tests. The `Option`/`Vec` absorption the brief worried about
  cannot happen: `WarrantedValue` refuses unknown fields and requires three, so a truncated inline
  table fails to parse rather than deserialising short.
- **A `"`, a `#` or a newline inside an origin string cannot truncate anything**, because `wrapped`
  splits on whitespace and both characters are ordinary inside a TOML comment.
- **Two of the three worked edits from spec §1.2 goal 3 now work from the file alone**, verified
  empirically by editing the written text and parsing it back.

### Wrong numbers in the diff's own prose (review step 8a)

Thirteen claims re-derived and correct; **two wrong and one misleading**, all the author's own:

| claim | truth |
|---|---|
| "the contamination section's note is the longest" | it is 10 lines; `[repeat_tracts]`' is 30 |
| "the **four** defaults ... in one `origins` module" | **five** constants — and the same sentence correctly says "five call sites" |
| "only **two** of those are per-row" | two *lines*, but three notes fire totalling seven lines across the file; the argument it supports — fixed cost, not per-row — is unaffected |

## 7. Missing tests to add now

`a_note_wraps_by_characters_at_the_width_it_states`; the flat-concentration text added to
`a_run_that_defaulted_nothing_writes_no_per_row_notes`; the two unreached origins exercised by
defaulting them; and each note asserted **above its own key** rather than somewhere in the file.

## 8. What's good

- **The reader's test is the right review for this step**, and it is repeatable: the same exercise
  on B2's uncommented file produced the list of questions this step was measured against.
- **Three of the previous reader's six wrongly-guessed keys are now answered**, and two of the
  three worked edits succeed from the file alone.
- **`the_comments_change_what_a_reader_learns_and_not_what_it_reads`** — stripping every comment
  must leave a file that reads back equal, and so must leaving them in. It is what makes the
  truncation hazard testable.
- **The origins are one module rather than five call sites**, so a number's origin and the sentence
  shown to a reader cannot drift apart, and step E1 has one list to reconcile.

## 9. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file`
- `./scripts/dev.sh cargo test --lib`

Audit trail: the two per-category files in the gitignored
`tmp/review_2026-08-28_parameters-file-b3/`.
