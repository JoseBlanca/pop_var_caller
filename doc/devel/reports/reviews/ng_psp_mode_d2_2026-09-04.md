# Code Review: ng_psp_mode_d2
**Date:** 2026-09-04
**Reviewer:** rust-code-review skill (orchestrator; two sub-agents in isolated worktrees, seven category checklists between them)
**Scope:** step D2's uncommitted diff — the calling loop lifted out of `AlignedFilesVariantCaller`, on `652fea99`
**Status:** Request-changes (all applied — see the fix report)

---

### 1. Scope

- Reviewed: the working-tree diff of plan step D2, exported as `tmp/review_2026-09-04_d2/d2.patch`.
- Against: commit `652fea99` + the patch, branch `ng-psp-mode`.
- In-scope: [callers.rs](../../../src/ng/run/callers.rs) and [observation_cache.rs](../../../src/ng/run/cohort_merge/observation_cache.rs).
- Categories (7): refactor_safety, reliability, errors, naming, module_structure, idiomatic, smells. Skipped: defaults, unsafe_concurrency, tooling, extras.

### 2. Verdict

**Request-changes**, and on nothing about the lift itself. **Both agents confirmed the body did not change**, by mechanical diff rather than by reading: the lifted 100 lines against the original have exactly four hunks and all four are the intended ones. What the review found is one missing test that D2 is what finally makes writable, and a set of documentation defects — including two the step introduced and one it inherited and kept.

The second agent settled the question this step exists to answer, and settled it by compiling rather than by arguing: it wrote Milestone E's method as `PspVariantCaller` will have to write it — `Vec<PspReader>` out of `self`, one `PspObservationSource::over` each, into `ObservationCache::over`, through the lifted function, sources read back afterwards — and `cargo check --lib --tests` came back clean in 36 s. **Nothing in the signature will force Milestone E to change it.**

### 3. Execution status

- Both agents detached at `652fea99`, applied the patch, verified the marker symbol, and restored their trees (one by blob hash against the patch, one by content).
- Agent mutations: **5 run, 1 survived** — the error-precedence swap, which is the Major below.
- Orchestrator, after the fixes: that mutation re-run and now killed by the new test alone; the byte-identity oracle re-run over the reshaped `WrittenCohort`.

### 4. Open questions and assumptions

1. **Should the five tallies be duplicated between `CohortCallingOutcome` and `WrittenCohort`, or factored?** Factored, on the second agent's recommendation and its argument: making `WrittenCohort` wrap the generic outcome would make a run-report type generic over the observation source, which is the leak D2 exists to prevent, so the shared five became a non-generic value both hold. Its evidence that the duplication was already costing something: two facts lived in only one of the two copies on the day they were written.
2. **`pub` or `pub(crate)` for the three new items?** `pub(crate)`. Thirteen references, all inside `callers.rs`, one real call site — and Milestone E's caller lands in that same file.

### 5. Top 3 priorities

1. **The error-precedence rule had no test**, and the two failures can genuinely be live at once.
2. **A new method landed inside another method's doc comment**, so rustdoc documented each under the other's first paragraph.
3. **A doc claim that a test does not exist**, when it has existed since 2026-09-01.

### 6. Findings

#### Major

- `callers.rs` — **A refused record outranks a source that fails afterwards, and nothing pinned it.** Swapping the two statements left all 474 `ng::run` tests green; on a fixture built for it, the clean code returns `RecordNotWritten` at chr1:15 and the mutant returns `SourceFailed` for the same input. The logic is lifted verbatim and is correct — what was missing is the test, and **D2 is what makes it writable**, because the loop now takes any source and a `Vec`'s iterator is one. **Categories:** reliability, refactor_safety.

#### Minor

- `observation_cache.rs` — **`sample_count` was inserted into the middle of `into_sources`'s doc comment**, so rustdoc documented a `usize`-returning accessor as "The sources back, in the run's sample order" and left `into_sources` with a headless fragment. Found independently by both agents; one quoted the generated HTML.
- `callers.rs` — **the doc credited direct mode's concurrency invariance to a test in the wrong module and said the end-to-end fixture "is not built yet".** `the_parallel_cover_gives_the_serial_drivers_answer` is at the merge; the end-to-end oracle is `the_record_path_is_byte_identical_at_every_thread_count`, landed 2026-09-01, which compares VCF bytes at pools of 1/2/4/8/16. Pre-existing, but this step trimmed that paragraph to three sentences and kept the wrong claim in them.
- `callers.rs` — **`CohortCallingOutcome::sources` was documented for a failed run**, a state the type cannot represent: it is built only after both error returns. The sentence was adapted from `ObservationCache::into_sources`, where it is true.
- `callers.rs` — **nothing exercised the loop's source-agnostic claim**: every instantiation in the tree was `Source = RunWalker`.
- `callers.rs` — **the new name did not carry the file's own accumulate-versus-hand-over distinction**, so it read like a sibling of `call_cohort` rather than of `call_cohort_handing_each_record_over`.

#### Nits

`CohortCallingOutcome` derived no `Debug` where `WrittenCohort` does — measured, that blocks `expect_err` in a test; "`LocusGeneratorSettings` and its neighbours" is one neighbour, `TractGeneratorSettings`; the five shared field docs duplicated between two structs.

### 7. What the "did anything change?" pass established

- The mechanical diff of the lifted body has **four hunks, all intended**: the cache moved out, two arguments became already-references, and rustfmt reflowed one `map_err` closure. The prologue left behind is identical apart from the one deleted `run_sample_count` line.
- `cache.sample_count() == walkers.len() == self.samples.len()` is proved through `walkers()` (one push per sample, `Err` before the vector is used) and `ObservationCache::over` (one window per source). A `+ 1` mutation kills ten tests, so the substitution is load-bearing rather than cosmetic.
- The padding accessor is still minted once, at the same point, and the compiler enforces that it is minted before `walkers()` consumes the caller. Its drop order relative to the cache inverts on the error path and that is inert: `WindowedRefSeq` has no `Drop`, and the load-bearing order lives inside `SampleLocusObservationsIterator`.
- Dropping the sources on the error path was already true at `652fea99`.

### 8. Out of scope observations

- Arch §3.4's `PspVariantCaller::open` sketch still lists a `callers_in_flight: CallersInFlight` argument. Arch §8 struck that question on 2026-09-01 and no such type exists in `src/`; the sketch is the stale copy, not a change the signature forces.

### 9. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib                      # 6,124 passed
./scripts/dev.sh cargo test --lib 'ng::run'            # 475 passed
./scripts/dev.sh cargo test --tests                    # 21 binaries, all ok
./scripts/dev.sh tmp/d2_oracle.sh before | after       # the VCF byte comparison
```
