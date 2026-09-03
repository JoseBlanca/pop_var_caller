# Code Review: ng_psp_mode_b2_b3
**Date:** 2026-09-03
**Reviewer:** rust-code-review skill (orchestrator; six category sub-agents, one isolated worktree each)
**Scope:** step B2+B3's uncommitted diff — three tests in `src/ng/run/gatherer.rs` + `examples/ng_psp_gather_oracle.rs`, on 4fcdcd5b
**Status:** Approve-with-changes

---

### 1. Scope

- What was reviewed: the working-tree diff of plan steps B2 and B3, exported as `tmp/review_2026-09-03_b2_b3_oracle/b2b3.patch`.
- Reviewed against: commit `4fcdcd5b` + the patch, branch `ng-psp-mode`.
- In-scope: the three new tests (`analysed_but_empty_ground_round_trips`, `tract_bearing_ground_round_trips_the_walk`, `gathering_the_same_sample_twice_gives_identical_bytes`) in [src/ng/run/gatherer.rs](../../../src/ng/run/gatherer.rs), and the whole of [examples/ng_psp_gather_oracle.rs](../../../examples/ng_psp_gather_oracle.rs).
- Out of scope: the rest of `gatherer.rs` (reviewed by the B1 nine-agent pass), `src/ng/psp/` internals (read as callees).
- Categories dispatched (6): reliability, errors, naming, idiomatic, smells, extras. Skipped: `defaults` (the diff adds no default-acting value — the settings it uses are B1's, already pinned), `module_structure` (no module change; the example is a leaf), `unsafe_concurrency` (none), `tooling` (no manifest change).

### 2. Verdict

**Approve-with-changes.** Both oracles do what the plan asks and the mutation evidence is unusually clean: three of the new tests each killed *exactly* the defect class they were written for, and nothing else. What must change is one hole the byte-identity test leaves open (a self-stamping clock passes it vacuously) and two robustness defects in the harness itself, where a bad environment variable panics and a failing directory entry can drop a sample out of the oracle without a word.

### 3. Execution status

- Every agent detached at `4fcdcd5b`, applied the patch, and verified branch-only files first.
- Mutation totals (reliability): **6 run, 4 killed, 2 survived, 0 changed-no-behaviour**, each survivor probe-proven behaviour-changing; each revert verified by reading the diff content.
- Orchestrator-side, quoted to agents: `cargo fmt --check` clean; `clippy --all-targets --all-features -- -D warnings` exit 0; `cargo test --lib 'ng::run::gatherer'` → 16 passed. The example ran for real in the container on one tomato slice: 183,807 records equal, header equal, second gather byte-identical, exit 0.
- Findings labeled "Needs verification": 0.

### 4. Open questions and assumptions

1. **Is a tract-free slice a legal input to the oracle?** Mi2 assumes not — the harness's stated coverage of the tract path should be enforced, not printed. Applied as a hard failure; if a later run wants generic-only ground, the flag to add is an explicit one.
2. **Contig lengths narrowing to `u32`** (errors, Low confidence) is the `regions` layer's repo-wide shape, copied from the sibling example. Left alone here; it is an upstream question about the caller's committed input range, not this diff's.

### 5. Top 3 priorities

1. **M1** — assert the file carries the *caller's* timestamp; without it a gatherer that stamps its own clock passes B3 vacuously (mutation-proven).
2. **M2/M3** — the harness's own failure paths: a malformed `NG_SAMPLES` panics outside the documented exit contract, and a failing directory entry silently removes a sample from an oracle whose whole claim is "exit 0 means every sample was checked".
3. **Mi1/Mi2** — two assertions that make the tests say what their doc comments claim: `regions_handled` in the empty-ground test, and a tract-records-present check in the example.

### 6. Findings

#### Major

- `src/ng/run/gatherer.rs:906` — **M1: nothing asserts the file carries the caller's timestamp, so a clock-reading gatherer passes every test**
- **Categories:** reliability (mutation-proven survivor A)
- **Confidence:** High
- **Problem:** The byte-identity test's doc claims "nothing in the gatherer reads a clock … identity holds over the whole file, timestamp included". A mutant where `open` discards `provenance.created` and stamps the wall clock at second precision passed all 16 tests: the two gathers land in the same second, so the files are still identical. No test asserts `writer.created`.
- **Why it matters:** §12.1's byte-identity rests on the timestamp being caller-supplied; the regression would pass CI deterministically and surface only as irreproducible psp files across runs that straddle a second.
- **Fix:** one assertion, the probe that killed the mutant — read the first file's header and compare `writer.created` against the fixture's.

- `examples/ng_psp_gather_oracle.rs:114` — **M2: an unparseable `NG_SAMPLES`/`NG_REGIONS` panics instead of erroring, outside the harness's documented exits**
- **Categories:** errors, idiomatic, smells (convergent)
- **Confidence:** High
- **Problem:** `how_many` parses an environment variable with `.expect("a count")`. `NG_SAMPLES=all` aborts at exit 101 with a message naming neither the variable nor the value, where the harness documents 0/1/2 and renders every other failure with its cause chain.
- **Fix:** make `how_many` fallible and propagate through `run`, naming the variable and the offending value; refuse zero there too (Mi8).

- `examples/ng_psp_gather_oracle.rs:171` — **M3: a failing directory entry is dropped silently, so a sample can escape the oracle**
- **Categories:** errors, idiomatic (convergent)
- **Confidence:** Medium (the discard is certain; entry errors are rare)
- **Problem:** `read_dir(crams)?.filter_map(Result::ok)` swallows per-entry errors. One failing entry means that sample is never gathered or compared and the run still exits 0 — in an oracle whose contract is "exit 0 means every sample's file is its walk".
- **Fix:** collect the listing fallibly so any entry error propagates.

#### Minor

- `src/ng/run/gatherer.rs:1128` — **Mi1: the empty-ground test cannot tell "walked and found empty" from "skipped for having no reads"** (reliability, extras — convergent). It asserts `regions_in` (dispatched) but not `regions_handled` (routed to a filled generator); a walk that skipped read-less ground writes a byte-identical file, so the tally is the only observable and its discriminating half is unasserted. **Fix:** assert `regions_handled == 1`.

- `examples/ng_psp_gather_oracle.rs:240` — **Mi2: the example prints the tract-record count and never checks it** (reliability). A tract-free run prints "0 … repeat tracts" and still exits 0, so the doc's claim that the tract path is on this oracle holds only while someone reads the printout. **Fix:** refuse a walk with no tract records, naming what to widen.

- `src/ng/run/gatherer.rs:1151` — **Mi3: the tract fixture is a homopolymer, so a left/right-symmetric store defect round-trips undetectably** (reliability, mutation-proven: a flank swap survived all 16). Covered one layer down by `record.rs`'s asymmetric round-trips (3 tests fail under that mutant) and by the real-data example. **Fix:** state the limitation on the test rather than leave the next reader to rediscover it.

- `examples/ng_psp_gather_oracle.rs:27` — **Mi4: the chain-equality pin the doc cites covers only generic ground** (extras). The example's in-memory side is a second gatherer, warranted by `the_gatherer_yields_what_the_direct_walk_yields` — whose fixture has no tract segment, while the measurement it warrants counts 1,217 tract records. **Fix:** generalise the test helper `direct_walk` to take a segmentation and compare the tract walk against the bare chain too — one assert, no new fixture.

- `examples/ng_psp_gather_oracle.rs:85` — **Mi5: `NG_TWICE` — the switch that runs B3's real-read oracle — is missing from the usage message** (extras); a rerun without it silently skips B3 and prints success.

- `examples/ng_psp_gather_oracle.rs:236` — **Mi6: the middle `open()` exists only to clone a header the first gatherer already held** (idiomatic). Each open re-reads the CRAM headers — real work on real data — and the clone is billed to the comparison walk's timing.

- `examples/ng_psp_gather_oracle.rs:217` — **Mi7: `std::slice::from_ref(&path.to_path_buf())` borrows a statement-scoped temporary** (idiomatic); it compiles by temporary-lifetime extension and a small refactor breaks it. **Fix:** hoist a `let alignments = [path.to_path_buf()];`.

- `examples/ng_psp_gather_oracle.rs:183` — **Mi8: `NG_SAMPLES=0` dies blaming the CRAM directory, `NG_REGIONS=0` the BED** (errors) — context that misleads. Folded into M2's fallible `how_many`.

- `src/ng/run/gatherer.rs:1092` — **Mi9: the segmentation-building scaffold now exists three times in one test module** (smells) — the patch takes it from one copy to three. **Fix:** a file-local `build_segmentation(segments, analysed)`; the deferred run-level fixture hoist then replaces one function rather than three literals.

- `src/ng/run/gatherer.rs:1209` — **Mi10: three shapes of "open a gatherer at the unusual settings"** (smells). The "every setting differs from its default" discipline now lives in three literals that can drift. **Fix:** `open_gatherer_over(alignments, reference, segmentation)`, with `open_gatherer` delegating.

- `examples/ng_psp_gather_oracle.rs:208` — **Mi11: `gather_and_compare` is ~106 lines carrying three separable oracles** (smells) — the header check, B3's byte-identity block, and the record loop it interrupts. **Fix:** extract the latter two as named functions the module doc can cite.

- `examples/ng_psp_gather_oracle.rs:257,214,259` and `src/ng/run/gatherer.rs:1083` — **Mi12: four naming repairs** (naming): `header_the_walk_fixed` credits a different actor than the message three lines below it; the closure `open` is half a name for what the crate calls `open_gatherer`; `chr2_only` is a modifier without its noun; and `"zstd-compression-level"` is a magic key three sites must spell identically while only its *value* is exported. **Fix:** rename the three, and export `ZSTD_COMPRESSION_LEVEL_KEY` beside the value.

#### Nits

`NG_TWICE` firing on any value including `0` while the doc says `=1` (three categories); `report_first_difference` covering the record's fields without an exhaustive destructure, so a new field would be misreported; a header mismatch printing no field-level lead where a record mismatch gets one; the doc's "exits non-zero on the first sample" describing a stop-at-first loop that actually finishes every sample; `work` → `work_dir`; `read_back` naming a count here and records elsewhere; bare `first`/`second` for the two psp paths; "re-walk" in output against "the comparison walk" in the doc; `first_regions_of`'s parameter shadowing the `how_many` function; two `mut` bindings mutated only during init; a needless `Arc::clone` at the binding's last use.

### 7. Out of scope observations

- `entry.length as u32` (both this example and its sibling) narrows a `.fai` length silently; the whole `ContigBounds`/`regions` layer is `u32`. Pre-existing and repo-wide — an upstream question about the committed input range.
- `how_many`, `first_regions_of` and `main`'s cause-printer are verbatim copies from `ng_call_cohort_end_to_end`; the crate's examples are self-contained by convention, so this is a standing choice rather than this diff's defect.
- The run-level shared test-fixture module stays deferred (recorded at B1); Mi9's local helper is the interim it will absorb.

### 8. Missing tests to add now

`write_psp_records_the_callers_timestamp` (M1's assertion, folded into the byte-identity test); the `regions_handled` assertion (Mi1); the bare-chain comparison over tract ground (Mi4); the tract-records-present refusal in the example (Mi2). A tract round trip over *asymmetric* reference sequence is **not** written here — the fixture reference is a homopolymer by construction; the class is pinned by `record.rs` and by the real-data example, and Mi3 records that on the test.

### 9. What's good

- Three of the six mutations were killed by exactly one test each — the empty-file read path by `analysed_but_empty_ground_round_trips`, tract routing by the tract guard, a record-invisible trailer byte by the whole-file comparison. Each new test is individually discriminating for the path it claims.
- The whole-file byte comparison caught a defect no record-level comparison could see (mutation D: every record compared equal).
- Every number in the example's "What it measured" block was checked against the run output and all matched, including the derived ones ("200 kb" from 200000 analysed bases; "the rest are bundles and satellites" from `RegionKind`'s four variants).
- The recorded adaptation (second gatherer instead of the bare chain) is stated in the module doc at the point it matters, and the pin it cites genuinely exists.

### 10. Commands to re-verify

- `scripts/dev.sh cargo fmt --check`; `scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `scripts/dev.sh cargo test --lib 'ng::run::gatherer'`
- `NG_TWICE=1 scripts/dev.sh cargo run --release --example ng_psp_gather_oracle -- <ref> <catalog> <bed> <cram>`

### Author response convention
Address findings by identifier (M1–M3, Mi1–Mi12) in the fix-application report.
