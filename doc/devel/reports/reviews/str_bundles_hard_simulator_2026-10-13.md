# Code Review: str_bundles_hard_simulator
**Date:** 2026-10-13
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** uncommitted diff on branch `str-bundles-hard` (from main b2714614): hard-shaped STR bundles in the stuttering simulator, the new `bundle_runs.tsv`, and the redefined typing check
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the working-tree diff against b2714614 — the stuttering simulator lays a share of its STR bundles (`hard_bundle_share`, default 0.6) in hard shapes: repeats that abut, and repeats cut into pure pieces by one to three interruptions. Stutter in a hard shape is placed once per run of each allele's own sequence, found by the simulator's copy of the caller's run rule. The typing check now asks for the laid hull and for the caller's reference runs to equal the simulator's stutter sites.
- **Reviewed against:** b2714614 plus the uncommitted diff, in the worktree `pop_var_caller-str-bundles-hard`.
- **In-scope files:** `examples/ng_str_stutter_benchmark.rs`, `scripts/ng_str_stutter_check_bundles.py`, `scripts/ng_str_stutter_benchmark.sh`, `doc/devel/ng/spec/str_stutter_benchmark.md` (§11, §11.4), `doc/devel/ng/spec/str_bundle_loci.md` (§9, §10 step 5).
- **Deliberately out of scope:** `src/region_typing/reference_runs.rs` and `examples/ng_str_bundle_runs.rs` (unchanged; read only to check the copy and the output format); the experiment builds the comparison ran.
- **Categories dispatched** (three read-only sub-agents in this worktree, no git worktrees and no builds — the owner's deviation from the skill for this review):
  - reliability — always; the diff's new tests and the run-rule copy are its substance;
  - naming, idiomatic, smells (one agent) — always;
  - errors, defaults, refactor_safety, float_portability, extras (one agent) — a new setting with a default, a driver change, a claimed byte identity, and spec numbers to check against the diff's intent.
  - module_structure, unsafe_concurrency, tooling: not dispatched — one example file, no concurrency added, no `Cargo.toml` change.

The orchestrator also ran the checks the request named: the run-rule copy against the library, stutter-site bookkeeping, truth scored as a caller, the calibration table, the byte identity, the typing check, and the spec's numbers against `tmp/bench/summary.tsv`.

## 2. Verdict

**Approve-with-changes.** The simulator does what the spec says, and every one of the spec's comparison numbers matches the scored runs. The changes worth making before merge: tests that check the copied run rule on the sequences it is actually used on; one test that checks the sites against the caller's rule instead of only their properties; a driver step that fails silently; and the under-representation of bundles with no reference run, which the request asked to fix if cheap. It was cheap.

## 3. Execution status

Commands run, in the dev container of this worktree (`cargo -j 2`):

- `cargo fmt -- --check examples/ng_str_stutter_benchmark.rs` → exit 0.
- `cargo test --release -j 2 --example ng_str_stutter_benchmark` → `test result: ok. 22 passed; 0 failed`.
- `cargo clippy --release -j 2 --example ng_str_stutter_benchmark -- -D warnings` → exit 101, from a pre-existing library lint outside the diff (`src/run/edge_deletion_gathering.rs:201`, `single_match`); without `-D warnings`, zero warnings in `examples/ng_str_stutter_benchmark.rs`.
- Truth scored as a caller (`tmp/review_2026-10-13_str_bundles_hard/truth_as_caller.sh`, `SKIP_VCFEVAL=1`), on the bench runs `pcr_1s_30x` and `pcr_50s_30x`: every segment and bundle cell right — 1 sample: 1,240 variant and 1,760 reference segment cells, 200 and 300 bundle cells; 50 samples: 61,433 and 88,567 segment cells, 10,036 and 14,964 bundle cells.
- Calibration (`tmp/hg002/compare.py` on `tmp/calibrate/hard6`): reproduces §11.4's simulated column exactly — 2,923 of 3,000 typed hard; 79, 37, 17 and 19 in 100; hull 46 and 76; 6 bundles with no run.
- Not run: `cargo test --all-targets` (the request limited builds to the simulator; the library is untouched), `cargo doc`, `cargo audit`.
- Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

1. Should the simulator place stutter by anything other than the caller's own run rule? As built it cannot, so the benchmark cannot detect the rule putting stutter in the wrong stretch (affects Mi13). The review assumes this is intended and asks only that the spec say so.

## 5. Top 3 priorities

1. **M4** — bundles with no reference run: 6 of 3,000 hard shapes against HG002's 123 of 1,447. A case the caller handles differently (no stutter modelled) went untested.
2. **M1 / M2** — the tests guarding the run-rule copy never fed it a flanked allele window, and the hard-bundle stutter test passed on an empty answer.
3. **M3** — `ng_str_bundle_runs` failing stops the driver with nothing on the terminal.

## 6. Findings

### Major

- **M1: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L4092) — The run-rule agreement test never sees the inputs `stutter_sites_of` is called on**
- **Confidence:** High
- **Problem:** `the_run_rule_copy_agrees_with_the_library` compares `reference_runs` with `find_reference_runs` on reference hulls and on random bundle-like sequences. Production calls it (through `stutter_sites_of`) on every allele's whole window: 30 bases of flank either side, length changes, substitutions, a SNP beside the hull and deletions across its edge.
- **Why it matters:** the copy can drift from the library on exactly those inputs with the test green.
- **Suggested fix:** add every allele window of a cohort's bundles to the comparison.
- **Categories:** reliability; orchestrator (the request's first question).

- **M2: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L4137) — `a_hard_bundle_stutters_in_pure_stretches_at_the_floor` checks properties, not the rule**
- **Confidence:** High
- **Problem:** it asserts each site is pure, at its floor and not overlapping. A `stutter_sites_of` returning nothing passes (`hard > 50` counts bundles, not sites), as would one returning the members' pure pieces.
- **Why it matters:** the one in-file test of the hard sites cannot catch a wrong placement.
- **Suggested fix:** assert the reference allele's sites equal `find_reference_runs` on its hull, shifted by the flank; require sites to exist, and at least one bundle without any.
- **Categories:** reliability.

- **M3: [scripts/ng_str_stutter_benchmark.sh](../../../../scripts/ng_str_stutter_benchmark.sh#L131) — A failing `ng_str_bundle_runs` stops the run with no message**
- **Confidence:** High
- **Problem:** its stderr goes to `bundle_runs_check.log`; under `set -e` the driver exits, and the terminal shows only the previous stage's output. A `BIN_DIR` without the tool fails the same way.
- **Suggested fix:** on failure, print the log's tail to stderr and exit 1.
- **Categories:** errors.

- **M4: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L990) — Bundles with no reference run are almost never laid**
- **Confidence:** High (measured: 6 of 3,000 hard shapes, against 123 of HG002's 1,447 hard bundles)
- **Problem:** the interrupted repeats keep at least one piece at its floor, so nearly every hard bundle has a run. The caller calls a no-run bundle with no stutter modelled — a different path — and the benchmark said nothing about it.
- **Suggested fix:** with probability 123/1,447 interrupt every repeat and redraw until `reference_runs(hull)` is empty. The decision is drawn inside the hard draw, so a zero share stays byte-identical.
- **Categories:** orchestrator.

- **M5: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L1250) — The scanner copy cannot see the flank laid after the spacer** *(found in verification)*
- **Confidence:** High (one case reproduced)
- **Problem:** `scanner_probe` lays the hull with a spacer either side, one base per base of the first or last motif, and nothing beyond. On 3,000 hard-only shapes after the M4 fix, one, `(AAT)6 A8`, was typed with a hull five bases longer than laid: the typing's scanner read the `AAT` repeat on through the one-base spacer of the `A` tail into the flank. The simulator accepted the shape because its copy of the scanner never saw those flank bases. The separated shape uses the same probe.
- **Why it matters:** the driver's typing check catches it and stops the run, so no wrong number reaches a score, but a seed that draws such a bundle cannot be run. At this rate (1 in 3,000 hard shapes) a default run, with about 300 hard bundles, draws one about 1 time in 10. Both seeds checked here and the four 1-sample runs passed.
- **Suggested fix:** lay the spacer for every motif that reaches the hull's end, or check the accepted shape's typing again once the flank is laid and redraw the flank. Either changes `build_reference` and needs the zero-share identity re-proved.
- **Categories:** orchestrator.

### Minor

- **Mi1: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L1000) — The junction after a homopolymer tail abuts at 0.4, not 0.04.** The abutting draws were all made before the tail rule turned a repeat into a homopolymer. Fix: draw the junctions in order, after any tail rewrite. (reliability, extras)
- **Mi2: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L3969) — The zero-share test's doc comment names the allele table, which it does not hash**, and calls the hash "measured on" the old build, when it was taken from the new one after a file-by-file comparison. (reliability, extras)
- **Mi3: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L4045) — `hard_shapes_keep_their_rules` cannot tell `AbuttingAndInterrupted` from `Abutting`**; its `_` arm also covers `Separated`. (reliability, extras, naming)
- **Mi4: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L3025) — `site.pure_start - WINDOW_FLANK_BP` is unchecked** and would wrap in release if a reference site ever lay in the flank. (reliability)
- **Mi5: [scripts/ng_str_stutter_check_bundles.py](../../../../scripts/ng_str_stutter_check_bundles.py#L130) — A bundle the runs tool did not report shows as `other_runs … caller's runs None`**, which reads as a disagreement about runs. (errors)
- **Mi6: [scripts/ng_str_stutter_check_bundles.py](../../../../scripts/ng_str_stutter_check_bundles.py#L155) — The summary says the runs agree "in every matched one"**, true by construction, so it reads as a pass when bundles fail on runs. (extras)
- **Mi7: [scripts/ng_str_stutter_check_bundles.py](../../../../scripts/ng_str_stutter_check_bundles.py#L124) — A bundle missing from `bundle_runs.tsv` raises a bare `KeyError`.** (errors)
- **Mi8: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L3453) — `..Settings::default()` silently moved six existing tests onto hard shapes.** Say so where `small_cohort` builds its settings. (refactor_safety)
- **Mi9: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L281) — The default's source is cited as a gitignored script.** (defaults)
- **Mi10: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L1117) — Unnamed `30` and `3 + below(floor - 2)`.** (naming, defaults)
- **Mi11: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L1073) — The `_ => Interrupted` arm also absorbs `(false, false)`**, excluded only by an earlier `continue`. (idiomatic, refactor_safety)
- **Mi12: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L2082) — Two `(usize, usize, usize)` triples with different middle fields.** `reference_runs` returns start, length, period; `scanner_repeats` returns start, end, period. The naming agent filed this as Major; filed here as Minor, because no call site mixes them and the run-rule test would catch a swap inside `stutter_sites_of`. (smells, naming)
- **Mi13: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L2082) — The run-rule copy is one function, the library three, under another name.** A sync review re-maps it by hand. (smells)
- **Mi14: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L1250) — The spacer–hull–spacer probe is built in three places**; only one uses `scanner_probe`. (smells)
- **Mi15: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L1135) — Interruption kinds encoded as `draw.below(4)` with arms `2`, `3`, `_`.** (naming, idiomatic)
- **Mi16: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L990) — `repeats` names the list and then a copy count in `draw_hard_bundle_shape`, a ~100-line function with five rejection points.** (naming, smells)
- **Mi17: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L1108) — `draw_interrupted_member` returns an anonymous pair** whose length invariant lives in the doc comment. (smells)
- **Mi18: [scripts/ng_str_stutter_check_bundles.py](../../../../scripts/ng_str_stutter_check_bundles.py#L82) — `read_runs` returns whole rows indexed by bare column numbers**, and its docstring says it returns the runs cell. (naming)
- **Mi19: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L3212) — The share used is not recorded in the run's log.** (defaults)
- **Mi20: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L1169) — `is_rotation_of`, `majority_base` and `draw_interrupted_member` have no direct tests.** (reliability)
- **Mi21: [doc/devel/ng/spec/str_bundle_loci.md](../../../ng/spec/str_bundle_loci.md) §10 step 5 — "in hulls of 30 bases or less the two callers are level (13 and 13, 15 and 18)" overstates.** 15 against 18 of 29 cells is a gap of 3, not level. (orchestrator)
- **Mi22: [examples/ng_str_stutter_benchmark.rs](../../../../examples/ng_str_stutter_benchmark.rs#L2082) and §11.4 — The doc claims the flanks hold no run, so the sites are the hull's runs.** True of the reference allele only: in a derived allele a SNP beside the hull or an edge deletion can lengthen a run into the flank. The molecule then slipping on the longer stretch is honest PCR; the claim needs narrowing. The spec also leaves out that stutter is placed by the caller's own rule, so the benchmark cannot test that rule (open question 1). (reliability, orchestrator)

### Nits

Functions named as nouns (`stutter_sites_of`, `reference_runs`, `scanner_probe`, `majority_base`); `BundleShape::kind` lacks a doc comment; single-letter closure parameter `w` in `is_rotation_of`; `scanner_reads_it_as_separated` re-scans the shape the previous predicate scanned; `".".to_string()` allocation; the §11.4 "every file byte for byte" sentence omits that the printed summary changes; the driver needs `ng_str_bundle_runs` built even for stages that never use it.

## 7. Out of scope observations

- `src/run/edge_deletion_gathering.rs:201` — clippy `single_match` under `-D warnings` on the library, independent of this diff; fails `cargo clippy -- -D warnings` for any target. Separate fix.

## 8. Missing tests to add now

- `the_run_rule_copy_agrees_with_the_library` extended with every allele window of a cohort's bundles (M1).
- `a_hard_bundle_stutters_in_pure_stretches_at_the_floor`: reference sites equal `find_reference_runs` on the hull; some sites exist; some bundles have none (M2, M4).
- `hard_shapes_keep_their_rules`: the share of shapes with no run lies between 1 in 20 and 1 in 8 of 2,000 (M4).
- `is_rotation_of_returns_true_only_for_rotations_of_equal_length`, `majority_base_returns_earliest_on_tie`, `draw_interrupted_member_returns_pieces_summing_to_total_on_every_count` (Mi20).

## 9. What's good

- Byte identity at a zero share is designed in: no draw is spent choosing a shape when the share is 0 (`draw_cohort`), and a file-by-file comparison backs it.
- Stutter is applied last site first in `cut_molecule`, and the per-allele sites are sorted and never overlap, so a length change at one site cannot move another's coordinates.
- The typing check compares the caller's own reference runs (`ng_str_bundle_runs`) with the simulator's sites, so simulator and caller cannot silently disagree about where the reference slips.
- The spec's comparison table was re-derived from `tmp/bench/summary.tsv` and every cell (96 numbers plus the cell counts) matched.

## 10. Commands to re-verify

- `scripts/dev.sh bash tmp/review_2026-10-13_str_bundles_hard/test_sim.sh` — fmt, the simulator's tests, clippy.
- `scripts/dev.sh bash tmp/identity.sh` — byte identity at share 0 and at no bundles.
- `scripts/dev.sh bash tmp/review_2026-10-13_str_bundles_hard/verify.sh` — the typing check on two seeds and on 3,000 hard-only shapes.
- `bash tmp/review_2026-10-13_str_bundles_hard/truth_as_caller.sh <run>` — the truth scored as a caller.
