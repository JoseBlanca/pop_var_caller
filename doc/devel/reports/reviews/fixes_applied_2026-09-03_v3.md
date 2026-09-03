# Fix Application Report: ng_psp_mode_b1_2026-09-03.md

**Date:** 2026-09-03
**Source review:** `doc/devel/reports/reviews/ng_psp_mode_b1_2026-09-03.md`
**Source state reviewed against:** 1d19975c + the uncommitted B1 diff, branch `ng-psp-mode`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 4
- Minors: 12
- Nits: 6 (grouped)

### Outcome totals
- Applied: 16 (B1, M1–M4, Mi1–Mi12 except as below; most nits)
- Applied with adaptation: 1 (M1 — see log)
- Deferred: 3 items inside findings (see §5)
- Disputed / Already fixed / Failed validation / Blocked / Superseded / Awaiting: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib 'ng::run::gatherer'` → 0, **13 passed** (was 6)
- `cargo test --lib 'ng::run'` → 0, 456 passed; `cargo test --lib 'ng::read::input'` → 0, 228 passed
- `cargo test --all-targets --all-features` → 101, red **only** on the three documented pre-existing locus-dump behaviour failures (2 in `ng_generic_loci_dump`, 1 in `ng_ssr_loci_dump` — verified by name against PROJECT_STATUS's standing entry)
- Performance check → skipped: nothing on a `benches/` path changed

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files changed |
|---|---|---|---|---|---|
| B1 | Blocker | write_psp failure propagation untested | Apply | Applied | gatherer.rs (test `write_psp_propagates_a_walk_failure`) |
| M1 | Major | progress contract untested | Apply | Applied with adaptation | gatherer.rs (2 tests) |
| M2 | Major | applied settings not pinned, only recorded | Apply | Applied | gatherer.rs (fixture + guard test) |
| M3 | Major | bundle threshold in three uncoupled copies | Apply | Applied | walker.rs, callers.rs, gatherer.rs |
| M4 | Major | remaining error mappings unexercised | Apply | Applied | gatherer.rs (3 tests) |
| Mi1 | Minor | three-way match triplicated | Apply | Applied | read_groups.rs, input/mod.rs, gatherer.rs |
| Mi2 | Minor | `sample` field duplicates `header.sample` | Apply | Applied | gatherer.rs |
| Mi3 | Minor | silent empty-identity fallbacks | Apply | Applied | walker.rs (`fasta_path`), gatherer.rs |
| Mi4 | Minor | write_psp discards the walk tally | Apply | Applied | gatherer.rs (returns `(WriteStats, LocusCounts)`) |
| Mi5 | Minor | `open` at 132 lines | Apply | Applied | gatherer.rs (`header_for` extracted; Mi1 removed the match) |
| Mi6 | Minor | run/mod.rs front-door doc stale | Apply | Applied | run/mod.rs |
| Mi7 | Minor | source types not destructured exhaustively | Apply | Applied | gatherer.rs |
| Mi8 | Minor | Debug without Self destructure | Apply | Applied | gatherer.rs |
| Mi9 | Minor | manifest/trailer defaults invisible in rustdoc | Apply | Applied | gatherer.rs (`# Defaults`, `# Errors`, trailer sentence) |
| Mi10 | Minor | pairing not asserted in two-samples test | Apply | Applied | gatherer.rs |
| Mi11 | Minor | `NotOneSamplesFiles` grammatical mush | Apply | Applied | run/mod.rs → `FilesNotFromOneSample` |
| Mi12 | Minor | `counts()` untested | Apply | Applied | gatherer.rs (asserted in parity + round-trip tests) |
| Nits | Nit | qualified paths, redundant import, a/b bindings, repeated chain, `# Errors` | Apply | Applied | gatherer.rs |
| Nits | Nit | PspNotWritten double path print; truncation/exhaustion pins | Won't fix / Deferred | see §5 | — |

## 3. Questions asked and answers
None — non-interactive; the review's two open questions were resolved by the recorded assumptions (Q1: `write_psp` now returns the tally, cheap either way; Q2: the tract-fixture gap defers to B2, below).

## 4. Per-finding log (what changed, compressed)

- **B1 — Applied.** `write_psp_propagates_a_walk_failure`: FASTA deleted after open (it is opened lazily per contig), `write_psp` must return `SourceFailed` and the file must not read back whole. The review's verified body, adapted to the renamed fixtures.
- **M1 — Applied with adaptation.** `reached_advances_to_the_last_observations_reach` as written. The mid-walk test's first draft asserted the error's progress equals the *first* observation's reach and failed against correct code — the resident chr1 window keeps yielding after the FASTA is gone, so the walk legitimately advances further before chr2's fetch fails. The landed test tracks the last clean observation's reach and asserts the error carries exactly that; it still kills both proven mutants (never-advance ⇒ `NothingYet`; hard-coded `NothingYet` in the error arm).
- **M2 — Applied.** The shared BAM fixture now carries a MAPQ-30 read (between the default floor 20 and the unusual 37) and three reads stacked at one position; `unusual_locus_generator_settings` gains `max_snp_column_depth: 2`. New guard test `the_fixture_tells_applied_settings_from_the_defaults` varies each dimension alone and requires the walks to differ — so the parity test now fails under either applied-side substitution the defaults harness proved survivable.
- **M3 — Applied.** `generic_path_generators` takes `&SegmentationInputs` and derives the radius itself; the expression is deleted at `callers.rs`, the gatherer, and the test copy (the two `DEFAULT_BUNDLE_THRESHOLD` test sites in callers.rs now pass a default-criteria segmentation's inputs — same value, 15, by construction). The stale "second place" comment is superseded by the new doc on the function.
- **M4 — Applied.** `open_refuses_files_whose_headers_cannot_be_read` (nonexistent path → `FilesNotFromOneSample`, chain names the file), `open_names_the_sample_when_a_file_will_not_open` (unindexed BAM, `build_index_if_missing: false` → `OpeningSample { sample }`), `write_psp_names_the_path_when_the_file_cannot_be_created` (missing parent dir → `PspNotWritten` carrying the exact path). `RecordNotWritten` and the seal arm have no cheap trigger; noted in the review as structural.
- **Mi1 — Applied.** `ReadGroups::only_sample()` in `read_groups.rs` owns the classification and the per-read-group listing; `SampleReads::open_only_sample` and the gatherer both delegate. The gatherer's unreachable `[]` arm is gone (the helper's `NoFiles` maps to `NoAlignmentFiles`).
- **Mi2–Mi12, nits — Applied** as the table says; the review's suggested shapes were used with only mechanical adaptation. `write_psp` now returns `(WriteStats, LocusCounts)` and both consumers of the walk's tally are asserted (`loci_emitted == records written`; `regions_in == 2`).

## 5. Deferred findings to carry forward
- **The tract path through the gatherer is unexercised** (review §7): B2's oracle fixtures must include a repeat-tract segment; the drift channel itself is closed by M3.
- **Run-level shared test-fixture module** (smells): land with B2's shared-fixture work — the gatherer is the fifth consumer of the duplicated fixtures.
- **`LocatedWalk` extraction** (the walker/gatherer mirrored `next` bodies): the missing tests now pin both copies; extract if a third copy ever appears.
- Nit-level: exhaustion/truncation pins; `PspNotWritten`'s cosmetic double path print (won't fix — consistent with `RecordNotWritten`).

## 6–8. Disputed / Failed validation / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched code reachable from `benches/` (the psp writer itself is unchanged; the gatherer is new and unbenched).

## 10–11. Commands run and results
All in the container via `scripts/dev.sh`, on the exact tree committed: `cargo fmt` / `cargo fmt --check` → 0; `cargo clippy --all-targets --all-features -- -D warnings` → 0; `cargo test --lib 'ng::run::gatherer'` → 0 (13 passed); `cargo test --lib 'ng::run'` → 0 (456); `cargo test --lib 'ng::read::input'` → 0 (228); `cargo test --all-targets --all-features` → 101 with exactly the three documented pre-existing failures, listed by name in §1.

## 12. Notes
The M1 adaptation is the one place the review's own suggested test was wrong against correct code; the report body above records why, so the next reader does not "fix" the test back.
