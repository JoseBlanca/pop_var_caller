# Code Review: ng_psp_mode_b1
**Date:** 2026-09-03
**Reviewer:** rust-code-review skill (orchestrator; nine category sub-agents, one isolated worktree each)
**Scope:** step B1's uncommitted diff — `src/ng/run/gatherer.rs` (new) + `src/ng/run/mod.rs` additions, on 1d19975c
**Status:** Approve-with-changes

---

### 1. Scope

- What was reviewed: the working-tree diff of plan step B1 (`SampleObservationGatherer`), exported as `tmp/review_2026-09-03_b1_gatherer/b1.patch`.
- Reviewed against: commit `1d19975c` + the patch, branch `ng-psp-mode`.
- In-scope files: [src/ng/run/gatherer.rs](../../../src/ng/run/gatherer.rs) (new), [src/ng/run/mod.rs](../../../src/ng/run/mod.rs) (module wiring + `NotOneSamplesFiles` + `PspNotWritten`).
- Deliberately out of scope: the already-committed `SegmentationInputs` lift and bench fix (reviewed implicitly as context); `src/ng/psp/` and `src/ng/read/` internals (read as callers/callees only).
- Categories dispatched (9): reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, smells, extras — every "always" category plus defaults (public inputs struct), module_structure (multi-file), extras (stable-output + diff-matches-intent). `unsafe_concurrency` skipped: no unsafe, no new locking; the only `Arc` is a pass-through handle. `tooling` skipped: no `Cargo.toml` change.

### 2. Verdict

**Approve-with-changes.** The construction, header content, and clean-path round trip are solid and mutation-pinned (the header test alone killed 9 of the 12 killed mutants across harnesses). What must change before commit: the failure paths and the applied-side settings are unpinned — four independent harnesses produced 10 surviving mutations, all in territory the six tests cannot see.

### 3. Execution status

- Each sub-agent detached its worktree at `1d19975c`, applied the patch, and verified branch-only files before reviewing; four ran scoped container builds (`DEV_CPUS=2 DEV_MEM=12g … cargo test --lib 'ng::run::gatherer'`), baseline `6 passed; 0 failed` in every one.
- Mutation totals across the four mutation harnesses (errors, reliability, defaults, refactor_safety): **23 mutants run, 10 survived, 2 changed no behaviour, 11 killed.** Every survivor was proven behaviour-changing by a probe (or is filed as the changed-no-behaviour coverage gap it is).
- Orchestrator-side verification (quoted to agents): container `cargo fmt --check` clean, `clippy --all-targets --all-features -- -D warnings` exit 0, scoped suite 6/6.
- Findings labeled "Needs verification": 0.

### 4. Open questions and assumptions

1. **Will anything read a psp-mode walk's `LocusCounts`?** Mi4 assumes plan step C2's per-sample report wants regions handled/refused. The fix (return the tally beside `WriteStats`) is cheap either way; if C2 ends up not reading it, the accessor's "final once the walk is spent" claim gets trimmed instead. Affects Mi4.
2. **The tract path through the gatherer is unexercised** (no fixture has a repeat-tract segment). The bundle-threshold refactor (M3) removes the drift channel this made dangerous; the remaining coverage gap is assumed to be B2's to close, since B2's real-CRAM oracle runs through the real catalog. Affects M3's residual risk.

### 5. Top 3 priorities

1. **B1** — a swallowed walk error seals a psp every reader accepts (mutation-proven); add the failure-propagation test.
2. **M2** — the suite pins what the header *records*, not what the walk *applies*: a gatherer walking with default filters/settings while recording the configured ones passes 6/6. Make the fixture discriminate.
3. **M1/M3** — the two hand-copied rules (progress-and-error wrap; bundle threshold) each survived mutation in every harness that tried them; pin the first with tests, collapse the second's three copies into the shared function.

### 6. Findings

#### Blocker

- `src/ng/run/gatherer.rs:292-293` — **B1: `write_psp`'s failure propagation is untested; a swallow regression seals a psp every reader accepts**
- **Categories:** errors
- **Confidence:** High
- **Problem:** No test drives `write_psp` (or the Iterator impl) through a failing walk. A mutant that discards walk errors returned `Ok(WriteStats { records: 0, blocks: 0, bytes: 4269 })` from a failed walk and sealed a footer-complete file `PspReader::open` accepts — all six tests green. The doc promise ("a reader will refuse it as interrupted") and the writer's `spent` latch both hold only if the gatherer propagates.
- **Fix:** the errors agent's verified test (`write_psp_propagates_a_walk_failure`): delete the fixture FASTA after `open` (it is opened lazily per contig), assert `Err(SourceFailed)` and that the file does not read back whole.

#### Major

- `src/ng/run/gatherer.rs:337,342` — **M1: the progress contract is untested — mutants that never advance `reached`, or stamp `NothingYet` into `SourceFailed`, survive every harness that tried them**
- **Categories:** errors, reliability, refactor_safety (convergent — each proved it independently)
- **Confidence:** High
- **Problem:** No test reads `reached()` after a draw and none produces a mid-walk failure. The walker's twin logic has exactly these tests; the gatherer copied the code without them. Under either mutant, every psp-mode walk failure reports "it had produced no observations yet" wherever it died.
- **Fix:** two tests (verified by the agents): `reached()` advances to the first observation's `reach_position`; a mid-walk failure (reads on both contigs, FASTA deleted after the first draw) names the sample and an `After(_)` progress.

- `src/ng/run/gatherer.rs:158,220` (test rule at 434-437) — **M2: the tests pin recorded settings, not applied ones — walking with `ReadFilterConfig::default()` or `PileupGeneratorConfig::default()` while recording the configured values passes 6/6**
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Fixture insensitivity: every read is MAPQ 60 (above both thresholds), and 30 bp reads cannot tell span 4,321 from 5,000. The parity test builds its oracle with the *same* settings, so a substitution shows up only if the two configs produce different observations — which this fixture cannot make them do. A regression writes a psp built from differently-filtered evidence whose header says otherwise: wrong results, no red test.
- **Fix:** make the fixture discriminate (a MAPQ-30 read the unusual filter drops; stacked reads a lowered depth cap truncates), assert the gatherer's output excludes/caps them, and guard the guard with an `assert_ne!` between a default-settings walk and the unusual one.

- `src/ng/run/gatherer.rs:224` — **M3: the bundle-threshold rule now lives in three uncoupled copies, and the suite cannot tell this one from a hard-coded constant**
- **Categories:** refactor_safety, smells, reliability (convergent)
- **Confidence:** High
- **Problem:** `Bp(segmentation.inputs().repeat_tract_criteria.classification.bundle_threshold)` is copied from `callers.rs:448` and again in the gatherer's own test; mutation-proven that hard-coding it (or reading `StrRepeatCriteria::default()`) stays green, because the test segmentation is built with default criteria and no fixture has a tract segment. `callers.rs`'s own comment predicted the drift at the second copy; this is the third.
- **Fix:** make `generic_path_generators` take `&SegmentationInputs` and derive the threshold itself, deleting the expression at both call sites and the test copy — no caller can then pass a constant. Build the gatherer's test segmentation with non-default criteria.

- `src/ng/run/gatherer.rs:129-133,286-308` — **M4: the remaining error mappings are unexercised — `NotOneSamplesFiles` via unreadable read groups, `OpeningSample`, `RecordNotWritten`, both `PspNotWritten` arms**
- **Categories:** errors, reliability (convergent)
- **Confidence:** High (create-arm fix verified; the mid-stream arms are structural)
- **Problem:** The step's contract names seven refusal shapes; tests pin three. A future edit can reroute a mapping (create-failure into the wrong variant, a dropped `map_err`) with the suite green.
- **Fix:** three cheap verified tests — nonexistent path → `NotOneSamplesFiles` naming the path; unindexed BAM with `build_index_if_missing: false` → `OpeningSample` naming the sample; `write_psp` into a missing parent directory → `PspNotWritten` carrying the exact path. `RecordNotWritten` and the seal arm have no cheap trigger; a comment beside each records that they are exercised structurally only.

#### Minor

- `src/ng/run/gatherer.rs:129-153` — **Mi1: the one-sample three-way match is the crate's third copy of the `SampleNameMismatch` listing**
- **Categories:** idiomatic, smells, refactor_safety, module_structure, errors (convergent, five of nine)
- **Fix:** `ReadGroups::only_sample()` beside the table; `open_only_sample` and the gatherer both delegate. This also retires the gatherer's unreachable `[]` arm whose `NoAlignmentFiles` message would misstate if it ever fired (three categories flagged it independently).

- `src/ng/run/gatherer.rs:88` — **Mi2: the `sample` field duplicates `header.sample`** (idiomatic). Drop the field; read the header's copy. The walker's justification (its name would otherwise vanish into the iterator) does not transfer.

- `src/ng/run/gatherer.rs:171-182` — **Mi3: silent empty-value fallbacks in the header's identity** (idiomatic, defaults, errors — convergent). `unwrap_or_default()` can record `""` as the reference name; `filter_map` can silently shorten `input_alignments`. **Fix:** a `WalkReference::fasta_path()` accessor so the `Option` never reopens, and `map`+`expect` naming the invariant for the basenames.

- `src/ng/run/gatherer.rs:285` — **Mi4: `write_psp(self)` discards the walk tally `counts()` calls "final once the walk is spent"** (idiomatic). **Fix:** return the `LocusCounts` beside `WriteStats`; consuming `self` stays (a `&mut` version called twice would seal a valid-looking empty psp).

- `src/ng/run/gatherer.rs:109-240` — **Mi5: `open` is a 132-line function** (smells). **Fix:** extract the one-sample resolution (falls out of Mi1) and the header assembly along the existing comment boundaries.

- `src/ng/run/mod.rs:1-22` — **Mi6: the module front-door doc now misstates the module's scope** (module_structure, errors). The "Landed so far" list and "Two ways out" predate the gatherer. **Fix:** extend the overview (and close the pre-existing `report` omission).

- `src/ng/run/gatherer.rs:186-216` — **Mi7: the `ContigInfo`/`ReferenceInfo` → identity mappings read fields without destructuring** (refactor_safety). A field added to the *source* types silently bypasses the only production header builder. **Fix:** exhaustive source destructures naming the deliberate discards.

- `src/ng/run/gatherer.rs:314-322` — **Mi8: the manual `Debug` picks fields without a `Self` destructure** (refactor_safety); the announced Milestone-G census field would be absorbed silently. **Fix:** destructure with named `_` discards.

- `src/ng/run/gatherer.rs:95-108,215` — **Mi9: the manifest default and the empty trailer are invisible in rustdoc** (defaults). **Fix:** a `# Defaults` paragraph on `open` naming `Manifest::as_this_build_writes_it`, and one sentence on `write_psp` for the empty trailer; plus the missing `# Errors` section on `open` (extras nit — 32 files under `src/ng` follow that convention).

- `src/ng/run/gatherer.rs:704-720` — **Mi10: the two-samples refusal test cannot see the pairing the variant exists for** (errors, reliability — convergent). A mutant reversing the `names` vector survives. **Fix:** assert the rendered `'{file}' names '{sample}'` pairs (verified to kill the mutant).

- `src/ng/run/mod.rs:381` — **Mi11: `NotOneSamplesFiles` is a dropped-apostrophe possessive that reads as two stacked plurals** (naming). **Fix:** rename to `FilesNotFromOneSample`; display text unchanged; three mechanical call sites.

- `src/ng/run/gatherer.rs:261-265` — **Mi12: `counts()` is untested** (reliability). **Fix:** assert `loci_emitted` and `regions_in` after draining, inside the existing parity test.

#### Nits

Fully-qualified `crate::ng::locus_generation::pileup::PileupGeneratorConfig` at the field and test helper where a `use` exists (four categories); redundant `use std::path::PathBuf` in tests; test bindings `a`/`b` outliving the abbreviation scope; the repeated basename chain in one test; `PspNotWritten` printing the path twice in a rendered chain (cosmetic, matches `RecordNotWritten`'s shape); truncation-on-create and exhaustion-is-final inherited contracts pinned nowhere at this level.

### 7. Out of scope observations

- `callers.rs:453` — the "this is the second place it has to reach" comment goes stale (or moot) with M3's fix; adjust when applying.
- `walker.rs:175-183` — the walker's own `Debug` has the same no-destructure shape as Mi8 (pre-existing).
- The run-level test-fixture duplication (the `[7; 16]` catalog header at its 4th copy, `index()`, `unusual_read_filters()` unreachable behind `callers`' private tests module) wants a `src/ng/run/`-level `#[cfg(test)]` fixtures module — a follow-up, natural to land with B2's shared-fixture work.
- No gatherer fixture carries a repeat-tract segment (reliability, Medium confidence): after M3's fix removes the drift channel, the residual is a coverage gap B2's real-CRAM oracle should close — B2's fixture must include tract ground.

### 8. Missing tests to add now

All named by the agents with verified or sketched bodies; the apply stage implements them: `write_psp_propagates_a_walk_failure` (B1), `reached_advances_to_the_last_observations_reach` + `a_midwalk_failure_names_the_sample_and_its_progress` (M1), the discriminating-fixture assertions (M2), `open_refuses_files_whose_headers_cannot_be_read`, `open_names_the_sample_when_a_file_will_not_open`, `write_psp_names_the_path_when_the_file_cannot_be_created` (M4), the pairing assertions (Mi10), the `counts()` assertions (Mi12).

### 9. What's good

- The header test's non-default-fixture discipline killed 9 of the 11 killed mutants across four independent harnesses — the `assert_ne!`-guarded "not the default" rule works ([gatherer.rs:559-564](../../../src/ng/run/gatherer.rs)).
- The byte-determinism claim was checked empirically, not argued: the same sample gathered twice produced byte-identical files (extras' probe), and no production line in the file reads a clock.
- The dependency direction ruled at Checkpoint A holds measurably: zero `use` of `ng::run` under `src/ng/psp/` after the patch (module_structure's grep).
- Every mutation harness restored its tree and proved it (`git apply --reverse --check`), and survivors were probe-proven behaviour-changing rather than assumed.

### 10. Commands to re-verify

- `scripts/dev.sh cargo fmt --check`
- `scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `scripts/dev.sh cargo test --lib 'ng::run::gatherer'`
- After fixes: the same scoped run, expecting the new tests listed in §8.

### Author response convention
Address findings by identifier (B1, M1–M4, Mi1–Mi12) in the fix-application report.
