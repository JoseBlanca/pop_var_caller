# Code Review: ng_psp_mode_d1
**Date:** 2026-09-04
**Reviewer:** rust-code-review skill (orchestrator; two sub-agents in isolated worktrees, seven category checklists between them)
**Scope:** step D1's uncommitted diff — the psp-backed observation source, on `638cfe38`
**Status:** Request-changes (all applied — see the fix report)

---

### 1. Scope

- Reviewed: the working-tree diff of plan step D1, exported as `tmp/review_2026-09-04_d1/d1.patch`.
- Against: commit `638cfe38` + the patch, branch `ng-psp-mode`.
- In-scope: [psp_source.rs](../../../src/ng/run/psp_source.rs) (new), [run/mod.rs](../../../src/ng/run/mod.rs), and the [implementation report](../implementations/ng_psp_mode_d1_2026-09-04.md).
- Categories (7): reliability, errors, refactor_safety, naming, idiomatic, smells, module_structure. Skipped: defaults (no configuration or default-acting value), unsafe_concurrency (no primitive; the one concurrency claim is a `Send` assertion), tooling (no manifest change), extras (applied piecemeal inside the two agents).

### 2. Verdict

**Request-changes**, on one defect both agents found independently and reached by different routes: **the doc's account of what happens after a failure was wrong for two of the three failure paths, and one of those two lost a record in silence.** Measured by the reliability agent: after an out-of-order refusal the next draw returned `Some(Ok(contig 0:201-201))` — the refused record at position 1 gone, the stream carrying on as if sound. That is the exact failure class the module's own doc says it exists to prevent.

Everything else was Minor or a nit. Both agents ran the diff's quantitative claims rather than re-reading them, and **all six checkable claims in the implementation report and the doc comments were correct** — including the fixture's 8 blocks and 2 contigs, the 459→470 test count with both sides measured, and three mechanism claims re-derived by mutation.

### 3. Execution status

- Both agents detached at `638cfe38`, applied the patch, verified `src/ng/run/psp_source.rs` existed, and restored their trees (verified by checksum).
- Agent mutations: **6 run, 2 survived** (`reached`'s `reach_position`, and the sample name in the walk-start failure).
- Orchestrator, before the review: 5 mutations run, 5 killed. After the fixes: **8 run, 7 killed** — the survivor is the one arm no fixture can reach through a `PspReader`, now marked as uncovered in the code.

### 4. Open questions and assumptions

1. **Should a refusal end the source, when the trait says a failure leaves it live?** Ruled here: it ends it, and every later draw says so rather than answering `None`. The trait's clause exists so a cover can be made again; nothing in the merge does that, and the two shapes of silence — the next record, or exhaustion — are both a sample quietly short at a locus. Recorded on the new variant.
2. **Is `over` (was `over(sample, walk)`, now `new`) part of the public surface?** No. The name it carries is whatever a caller passes, and three documents asserted the name comes from the file. `pub(crate)` keeps every caller E1/E2 will need and closes the gap.

### 5. Top 3 priorities

1. **The post-failure paragraph, and the silent record loss underneath it** — both agents, Major.
2. **`reached` unpinned against a multi-base observation**: every fixture asserting it was one base wide, so `reach_position` → `start_position` survived all eleven tests.
3. **A test named for a path it does not exercise**: the walk-start failure's error literal was untested and duplicated.

### 6. Findings

#### Major

- `psp_source.rs` — **The failure/liveness note was untrue for two of three paths, and one of them dropped the refused record.** The walk's own failure fuses the walk, so that path is the walker's documented deviation. The two refusals this file mints touched no latch: the source stayed live, the refused record vanished, and the walk carried on succeeding. **Categories:** reliability, errors, naming (both agents).

#### Minor

- **`reached` is never compared against a multi-base observation**, so the `reach_position` the comment argues for is unpinned — `a_sample()`'s records are all one base wide, where reach and start are the same number.
- **`a_walk_that_will_not_start_…` never called the constructor it was named for**, and that constructor's `RunError::SourceFailed` literal was a second copy of the one `refuse` builds; a fabricated sample name in it left the suite green.
- **The enum's doc and one variant's contradicted each other** about whose mistake a head-only record is — and the rendered message blamed a sound file, sending a reader to rebuild a psp that is not broken.
- **`StreamedRecord` was read by field access rather than destructured**, unlike `psp::reader`'s own consumer of the same shape; step E2's read-group remap lands in that exact function.
- **"Exhaustion is final" was asserted and unpinned** — no test drew again after `None` on a non-empty file.
- **`offered`'s doc named the head; the code filled it from the decoded record**, eight lines under a comment arguing the opposite choice for the neighbouring value.
- **The out-of-order message printed its two coordinates in two shapes** (`contig 0:1-1` against `contig 0 position 101`), where the store's own two out-of-order errors render symmetrically.
- **`a_psp_of_under` was a parameterised helper with one caller always passing the same argument.**

#### Nits

An unguarded `block_index()[2]` in the damage fixture; a doubled "and" in `run/mod.rs`'s landed list; a `#[derive(Debug)]` where the sibling hand-writes one and documents why (the derive also bounds the impl on `W: Debug`, which the type's whole generic-over-the-walk argument invites callers to break); "the three things the merge's trait asks for" whose third is a parameter the type discards; "and this is that day" (clear-writing Rule 6); `over`/`reading` naming the opposite roles from the walker's `new`/`over`, `reading` being the crate's only gerund constructor; the report's "the tail D2 lifts" leaving three plan labels doing the work of nouns.

### 7. Out of scope observations

- `observation_cache.rs`'s release assertion on a backwards source still stands for every non-psp source. Arch §8's item is half discharged, which the implementation report already says.
- `psp/walk.rs`'s seek arm has no fixture and says so — the reason the walk-start finding is Minor.
- Nothing asserts `Sync` for this source, and arch §2 names `S: Sync + Send` as what `merge_cohort_in_parallel` requires. Direct mode's walker is not `Sync` either and that driver is off by default, so this is a note for whoever re-opens the parallel merge, not a gap here.

### 8. Missing tests, added

Four: the refusal latch, `reached` against a six-base observation, exhaustion holding across repeated draws, and a body declined part-way through a walk (where `reached` can be wrong in a way the first-record fixture cannot see). The module goes from 11 tests to 15.

### 9. What's good

- **The coordinate-order boundary is genuinely pinned, not decorated**: one fixture carries both an overlap and a tie, and `<` → `<=` fails it.
- **The contig comparison is pinned by construction** — the round-trip fixture spans two contigs, so a check on position alone would refuse contig 1's first record.
- **One source per sample is enforced by the borrow checker**: the walk holds `&mut PspReader`, so a second source over one open psp does not compile.
- **The damaged-file test damages a real file** and asserts both that some records arrived and that what arrived equals what was stored, so a silently-wrong decode fails it too.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib 'ng::run'                 # 474 passed
./scripts/dev.sh cargo test --lib 'ng::run::psp_source'     # 15 passed
bash tmp/d1_mutations/run2.sh                                # 8 mutations, 7 killed
```
