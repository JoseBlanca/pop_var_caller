# Fix Application Report: str_bundles_simulator_2026-10-11.md

**Date:** 2026-10-11
**Source review:** `doc/devel/reports/reviews/str_bundles_simulator_2026-10-11.md`
**Source state reviewed against:** branch `str-bundles`, HEAD ccef7dca plus the uncommitted change under review
**Execution mode:** non-interactive (no owner questions could be asked; decisions needing the owner are deferred)
**Overall status:** Partial — every correctness finding applied; four refactors and one library API change deferred

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 3
- Minors: 13
- Nits: 5 groups

### Outcome totals
- Applied: 8 (M1, M2, M3, Mi1, Mi2, Mi3, Mi6, Mi11)
- Applied with adaptation: 4 (Mi4, Mi5, Mi10, Mi13)
- Already fixed: 0
- Deferred: 4 (Mi7, Mi8, Mi9, Mi12)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0
- Nits: two of the five groups applied (the check script's error and member order; `score_callers.py`'s docstring), the rest left.

**No fix changes what the simulator writes.** Same settings with bundles laid (`samples=2 depth=1
segments_per_period=60 bundles=300 genome_bp=800000 seed=5`): the reviewed build and the fixed build
write identical `reference.fa`, `truth.vcf`, `bundles.bed`, `bundle_alleles.tsv`,
`segment_alleles.tsv` and FASTQs.

### Validation summary
- `cargo fmt --check` on the two examples → 0, clean.
- `cargo clippy` on the two examples → no warning in either file (the vendored noodles crates warn, unchanged). Not run with `-D warnings` across all targets: the vendored crates' warnings fail that gate before this change.
- `cargo test --release --example ng_str_stutter_benchmark` → 0, `17 passed; 0 failed` (14 before; 3 added).
- `cargo test --all-targets`, `cargo doc`, `cargo audit` → not run: the owner's light dev loop runs touched tests only, and no library code changed.
- Performance check → not applicable (examples and scripts only; nothing under `benches/` reaches them).

### Unresolved high-priority findings
- None. The deferred findings are all Minor.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | vcfeval cache reuses a stale `outside` | Apply | Applied | No | `scripts/ng_str_stutter_score.py` | Pass (faked-tool check) | No |
| M2 | Major | bundle member runs and stutter stretches untested | Apply | Applied | No | `examples/ng_str_stutter_benchmark.rs` | Pass; 2 mutations caught | No |
| M3 | Major | scanner copy follows library scores; no differential test | Apply | Applied | No | `examples/ng_str_stutter_benchmark.rs` | Pass; 2 mutations caught | No |
| Mi1 | Minor | `outside` overlaps nested windows | Apply | Applied | No | `scripts/ng_str_stutter_score.py` | Pass (self-check) | No |
| Mi2 | Minor | truth-VCF test's segment counter counts bundles | Apply | Applied | No | `examples/ng_str_stutter_benchmark.rs` | Pass | No |
| Mi3 | Minor | "byte for byte" reads as holding at every setting | Apply | Applied | No | `examples/ng_str_stutter_benchmark.rs` | N/A (doc) | No |
| Mi4 | Minor | hull median's cause misstated | Apply | Applied with adaptation | No | `doc/devel/ng/spec/str_stutter_benchmark.md` | Figures re-measured | Owner decision |
| Mi5 | Minor | new default changes what old commands produce | Apply | Applied with adaptation | No | `scripts/ng_str_stutter_first_pass.sh` | N/A (doc) | No |
| Mi6 | Minor | `longest_period` misreads | Apply | Applied | No | `scripts/ng_str_stutter_score.py` | Pass | No |
| Mi7 | Minor | duplicated allele-draw arithmetic | Defer | Deferred | No | None | N/A | Yes |
| Mi8 | Minor | per-member parallel vectors | Defer | Deferred | No | None | N/A | Yes |
| Mi9 | Minor | anonymous tuples for members and scanner hits | Defer | Deferred | No | None | N/A | Yes |
| Mi10 | Minor | second `members_json`, wildcard arm | Apply (part) | Applied with adaptation | No | `examples/ng_typed_region_dump.rs` | Builds | Yes (shared fn) |
| Mi11 | Minor | shape redraw unbounded | Apply | Applied | No | `examples/ng_str_stutter_benchmark.rs` | Pass; outputs unchanged | No |
| Mi12 | Minor | each call file parsed twice | Defer | Deferred | No | None | N/A | Yes |
| Mi13 | Minor | `outside` changed meaning under its old name | Apply | Applied with adaptation | No | `scripts/ng_str_stutter_score.py` | N/A (doc) | No |

## 3. Questions asked and answers

None asked (non-interactive). One decision is left to the owner under Mi4.

## 4. Per-finding log

### M1 — vcfeval cache reuses a stale `outside`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** confirmed in the code: `run_vcfeval` skipped a stretch whenever `summary.txt` existed.
- **Implementation summary:** each stretch's result now carries `inputs.sha256`, the SHA-256 of the
  truth VCF, the call file and the stretch's BED. A missing or different digest reruns vcfeval.
  `truth.vcf.gz` is rebuilt when `truth.vcf`'s digest changes, and a packed call file is indexed
  every time it is rewritten.
- **Review suggestion used verbatim?:** No — same idea, with block-wise hashing and the helpers named `file_digest`, `read_stamp`, `write_stamp`.
- **Adaptation:** None of substance.
- **Verification performed:** a faked-tool harness (`tmp/review_2026-10-11_str-bundles-sim/cache_check/check.py`, since deleted with the other scratch runs) calling `run_vcfeval` five times. New scorer, reruns of vcfeval: first call 1, unchanged inputs 0, BED changed 1, truth changed 1, an old cache with no digest 1. HEAD's scorer on the same sequence: 0 reruns for the changed BED and the changed truth — the stale reuse the finding describes.
- **Files changed:** `scripts/ng_str_stutter_score.py`
- **Tests added or modified:** none in the repository (the scorer has no test suite); the harness above.
- **Validation:** `uv run --no-project python .../cache_check/check.py` → 0, `OK`.
- **User input:** None
- **Follow-up:** None
- **Residual risk:** a 50-sample call file is hashed once per caller per scoring, seconds against vcfeval's minutes.

### M2 — bundle member runs and stutter stretches untested
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** the code was correct but no test pinned the coordinates that move.
- **Implementation summary:** new test `every_bundle_allele_keeps_its_member_runs_on_its_members`
  over 150 bundles, every allele kind: the gap bases sit exactly between consecutive runs; an untouched
  run is its motif `count` times; a member cut by an edge deletion is periodic, holds `count` whole
  copies and ends (left cut) or starts (right cut) on a whole motif; a substitution changes exactly one
  base; every pure span lies in its run, is at least two units, and repeats at the member's motif
  length. It requires at least one deletion at each edge.
- **Review suggestion used verbatim?:** N/A
- **Adaptation:** None
- **Verification performed:** passes on the code; **two mutations, each caught by this test and by
  no other**: shifting later runs by `in_repeat` instead of all the deleted bases (left edge), and not
  shortening the last run (right edge). Mutations run 2, survived 0, changed-no-behaviour 0.
- **Files changed:** `examples/ng_str_stutter_benchmark.rs`
- **Tests added or modified:** `every_bundle_allele_keeps_its_member_runs_on_its_members`
- **Validation:** `cargo test --release --example ng_str_stutter_benchmark` → 0, 17 passed.
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M3 — scanner copy follows library scores; no differential test
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** the import contradicted the copy's stated purpose, and drift in either direction
  was invisible until a benchmark run stopped.
- **Implementation summary:** the scores are now local constants `SCANNER_MATCH_REWARD` (2) and
  `SCANNER_MISMATCH_PENALTY` (7). Two tests: `the_scanner_copy_uses_the_library_scores` asserts both
  and `LONGEST_BUNDLE_GAP` against the library (`DEFAULT_BUNDLE_THRESHOLD − 1`);
  `the_scanner_copy_agrees_with_the_library_scanner` runs the copy and `find_tandem_repeats` plus the
  copy floor on 3,000 bundle-like sequences whose gap bases are drawn with no rule, and requires
  identical repeats on all of them.
- **Review suggestion used verbatim?:** Yes for the constants test; the differential test is the orchestrator's.
- **Adaptation:** None
- **Verification performed:** the differential test passes, which also answers the Rust reviewer's
  question about tie handling: on these 3,000 sequences the recursive copy and the library's
  Ruzzo–Tompa agree. Two mutations of the copy's tie handling (`>` to `>=` when picking the best
  stretch; resetting on `< 0` instead of `<= 0`) are each caught by the differential test and by no
  other test. Mutations run 2, survived 0, changed-no-behaviour 0.
- **Files changed:** `examples/ng_str_stutter_benchmark.rs`
- **Tests added or modified:** the two above.
- **Validation:** as M2.
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### Mi1 — `outside` overlaps nested windows
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** one-token fix with one correct form.
- **Implementation summary:** `cursor = max(cursor, end)`.
- **Review suggestion used verbatim?:** Yes
- **Adaptation:** None
- **Verification performed:** scorer self-check still `right`/`ref_right` everywhere on three seeds (the window BEDs are written only when vcfeval runs, so the nested case itself was reproduced by the scripts reviewer, not rerun).
- **Files changed:** `scripts/ng_str_stutter_score.py`
- **Tests added or modified:** None
- **Validation:** see §11.
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### Mi2 — truth-VCF test's segment counter counts bundles
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** `segment_records` counts IDs starting `seg`.
- **Reasoning / Verification:** test passes; the counter now measures what its name says.
- **Review suggestion used verbatim?:** No — the prefix constants were not introduced (Nit-level).
- **Adaptation:** None
- **Files changed:** `examples/ng_str_stutter_benchmark.rs`
- **Tests added or modified:** `truth_vcf_reproduces_the_truth_sequences`
- **Validation:** 17 passed.
- **User input / Follow-up / Residual risk:** None

### Mi3 — "byte for byte" reads as holding at every setting
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** module doc now says the segments' files keep their columns, are byte
  for byte as before at `bundles=0`, and are a different cohort with bundles laid.
- **Reasoning / Verification:** doc only; matches the byte-identity runs.
- **Review suggestion used verbatim?:** Close to it.
- **Adaptation:** None
- **Files changed:** `examples/ng_str_stutter_benchmark.rs`
- **Tests added or modified:** None
- **Validation:** fmt clean.
- **User input / Follow-up / Residual risk:** None

### Mi4 — hull median's cause misstated
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** the sentence was fixed; whether the draw should change is not a review's decision.
- **Implementation summary:** §11.1 now says the median (77) is more than twice HG002's (35)
  because every member's count is uniform from the floor to the segments' highest (30 or 15 copies),
  so a member averages 31 bases (1,081 members at the default seed, measured), that long hulls leave
  fewer crossing reads at low depth, and that drawing counts nearer the floor is open.
- **Adaptation:** the draw itself is unchanged.
- **Verification performed:** the 30.9-base average and all §11.1 figures re-measured at the default seed.
- **Files changed:** `doc/devel/ng/spec/str_stutter_benchmark.md`
- **Tests added or modified:** None
- **Validation:** N/A
- **User input:** None
- **Follow-up:** **owner decision** — keep the uniform draw, or draw member counts nearer the floor so the hull median approaches HG002's.
- **Residual risk:** None

### Mi5 — new default changes what old commands produce
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Implementation summary:** the first-pass grid script's header says that a simulator built since
  2026-10-11 lays 500 bundles, that the 2026-10-10 cohorts come back only with `bundles=0` (which the
  script does not pass), and that `FIT_FROM` refuses a fit from a run made before bundles.
- **Adaptation:** documented rather than changing the default or the `FIT_FROM` comparison; both
  behaviours are safe once stated, and changing them is the owner's call.
- **Files changed:** `scripts/ng_str_stutter_first_pass.sh`
- **Tests added or modified:** None
- **Validation:** N/A
- **User input / Follow-up / Residual risk:** None

### Mi6 — `longest_period` misreads
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** the files are new, so no reader depends on the old column names.
- **Implementation summary:** columns renamed `longest_member_period` and
  `longest_member_repeat_band` in `score_bundle_cells.tsv` and `score_bundle_summary.tsv`; docstring
  updated.
- **Verification performed:** header of a freshly scored `score_bundle_summary.tsv` read back.
- **Files changed:** `scripts/ng_str_stutter_score.py`
- **Tests added or modified:** None
- **Validation:** self-check on three seeds.
- **User input / Follow-up / Residual risk:** None

### Mi7 — duplicated allele-draw arithmetic
- **Severity:** Minor
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** a refactor touching the segment draw path, which the byte-identity test pins; safe in
  principle but not needed for correctness, and outside a fix run's minimal-diff rule.
- **Implementation summary / Files changed / Tests:** None
- **Follow-up:** extract shared helpers and named constants, keeping the draw order; the hash test checks it.
- **Residual risk:** a change made to one copy only would give segments and bundles different populations.

### Mi8 — per-member parallel vectors
- **Severity:** Minor
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** restructures `BundleAllele` and every use; M2's test now catches a misplaced run,
  and `a_bundle_has_one_stutter_site_per_member` asserts one span per member.
- **Follow-up:** one `Vec<AlleleMember>` when the file is next restructured.
- **Residual risk:** low, now covered by tests.

### Mi9 — anonymous tuples
- **Severity:** Minor
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** naming refactor across the file; no behaviour at stake.
- **Follow-up:** `BundleMember`, `ScannedRepeat`.
- **Residual risk:** None

### Mi10 — second `members_json`, wildcard arm
- **Severity:** Minor
- **Initial decision:** Apply (wildcard); Defer (shared function)
- **Final status:** Applied with adaptation
- **Implementation summary:** the wildcard is replaced by the three named variants, so a new
  `RegionKind` fails to compile here.
- **Adaptation:** the shared public formatter is not made: it changes the library's public API, which
  this run does not do without the owner.
- **Verification performed:** the dump builds; its output is unchanged (the bundle check reads it on three seeds).
- **Files changed:** `examples/ng_typed_region_dump.rs`
- **Follow-up:** move `members_json` into the library and call it from both, if the owner wants the API.
- **Residual risk:** the two format strings could still drift; the bundle check would fail loudly if they did.

### Mi11 — shape redraw unbounded
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** `for _ in 0..BUNDLE_SHAPE_ATTEMPTS` (100,000) with a panic naming the
  rules after it, and a `# Panics` section.
- **Verification performed:** outputs with bundles laid identical between the reviewed and fixed builds (§1).
- **Files changed:** `examples/ng_str_stutter_benchmark.rs`
- **Tests added or modified:** None
- **Validation:** 17 passed.
- **User input / Follow-up / Residual risk:** None

### Mi12 — each call file parsed twice
- **Severity:** Minor
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** runtime only; folding bundle scoring into the segment loop restructures `main` and
  must keep `score_cells.tsv` byte-identical, which is worth doing with a timing on a 50-sample run.
- **Follow-up:** score bundles inside the per-caller loop.
- **Residual risk:** longer scoring on large runs; no wrong numbers.

### Mi13 — `outside` changed meaning under its old name
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Implementation summary:** the scorer's docstring says `outside` differs between runs with and
  without bundles, and when a cached stretch is reused.
- **Adaptation:** `grid_summary` not changed to record which kind of run a row came from.
- **Files changed:** `scripts/ng_str_stutter_score.py`
- **User input / Follow-up / Residual risk:** None

### Nits
- Applied: the bundle check stops with a message naming the missing `members` column when the dump
  predates it, and compares members sorted by start; `score_callers.py`'s docstring names the bundle files.
- Left: the 1,000-base look-back, stale bundle score files on an empty bundle table, the Rust naming
  nits, the probe's spacer, the realised second-member share.

## 5. Deferred findings to carry forward
- Mi7 — shared allele-draw helpers (refactor).
- Mi8 — one vector of per-member records (refactor).
- Mi9 — named member and scanner-hit types (refactor).
- Mi12 — one parse of each call file (runtime).
- Mi10's second half — a public `members_json` shared by the CLI and the dump (library API).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — examples and scripts only; nothing under `benches/` reaches them.
- **Baseline saved:** no (not triggered).
- **Outcome:** skipped.

## 10. Commands run
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-sim/test.sh` (rustfmt, build, tests)
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-sim/mutate.sh` (four mutations, source restored and compared)
- `uv run --no-project python tmp/review_2026-10-11_str-bundles-sim/cache_check/check.py` and the same against HEAD's scorer
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-sim/final.sh` (fmt, build, clippy on the two examples, tests, `verify.sh` on seeds 1, 424242 and 77)
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-sim/same.sh` (reviewed against fixed outputs; deletes the generated runs)

## 11. Command results
- `cargo fmt --check` (two examples) → 0, clean.
- clippy, two examples → no warning in them.
- `cargo test --release --example ng_str_stutter_benchmark` → 0, 17 passed.
- Mutations → left-edge shift and right-edge end caught by M2's test only; the two tie mutations caught by M3's differential test only.
- Byte identity at `bundles=0` (`design=heavy samples=3 depth=2 segments_per_period=60 genome_bp=500000 contigs=3 inbreeding=0.1 seed=7`) against HEAD's simulator → all 14 files identical (FASTQs decompressed).
- Typing check, 600 segments and 500 bundles on 1.5 Mb → seed 1, 424242 and 77: 500 of 500 bundles typed as drawn, offset 0, no bundle invented; 600 of 600 segments as planted.
- Truth VCF scored as a caller → every bundle and segment cell `right` or `ref_right` on all three seeds (bundles: 632/868, 617/883, 596/904 right/ref_right).
- Reviewed build against fixed build, bundles laid → identical reference, truth, bundle and segment tables, FASTQs.
- Cache harness → reruns 1, 0, 1, 1, 1 (new scorer); HEAD's scorer reruns 0 for a changed BED and a changed truth.

## 12. Notes
- **`PROJECT_STATUS.md` was not updated.** It has no block for the STR-bundles work and its current
  focus stops at 2026-10-10, before the edge-deletion merge; adding a block on this branch alone, in a
  750 KB file other branches edit, is better done by the owner when the branch merges.
- The generated runs and binaries under `tmp/review_2026-10-11_str-bundles-sim/` were deleted; the
  scripts, the reviewers' findings files and the two verification transcripts are kept.
