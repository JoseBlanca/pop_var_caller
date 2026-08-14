# Code Review: census_rename_a2
**Date:** 2026-08-14
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** commit `7f806fdf` — step A2 of the census plan: `records.rs` → `census.rs` plus ten type renames
**Status:** Approve-with-changes

---

### 1. Scope

- **What was reviewed:** one commit's diff — `git diff ng-joint-fit..7f806fdf`.
- **Reviewed against:** `7f806fdf` on branch `ng-census-encoding`.
- **In-scope files:** [census.rs](../../../../src/ng/parameter_estimation/joint/census.rs) (was `records.rs`), [loci.rs](../../../../src/ng/parameter_estimation/joint/loci.rs), [fit.rs](../../../../src/ng/parameter_estimation/joint/fit.rs), [ssr_fit.rs](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs), [contamination.rs](../../../../src/ng/parameter_estimation/joint/contamination.rs), [coverage.rs](../../../../src/ng/parameter_estimation/joint/coverage.rs), [mod.rs](../../../../src/ng/parameter_estimation/joint/mod.rs), and the four `examples/ng_joint_*` harnesses the rename touches.
- **Deliberately out of scope:** the coverage-by-window summary's existence (step A3 deletes it); the encoding changes (milestone B); the census file and its `Sections` shape (plan 2).
- **Categories dispatched, and why only four.** `naming` — the whole change is names. `refactor_safety` — the milestone's oracle is that no fitted number moves, so *is this behaviour-preserving* is the review's central question. `module_structure` — a file moved. `reliability` — always, and here it judges the author's position that a substitution needs no new test. **`errors`, `defaults`, `idiomatic`, `smells`, `tooling`, `extras` and `unsafe_concurrency` were not dispatched:** the diff adds no error type (`CensusError` does not exist yet), changes no default, no `Cargo.toml`, no `unsafe`, no concurrency primitive, and no expression — there is no changed surface for those checklists to bite on. Recorded as a deliberate triage, not an omission.

### 2. Verdict

**Approve-with-changes.** The rename is behaviour-preserving and complete; the changes asked for are the prose the `sed` pass damaged, four hardenings the reviewers converged on, and two tests for branches that mutated freely.

### 3. Execution status

Run by the orchestrator in the dev container:

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo check --all-targets` | `Finished dev profile … in 22.08s`, no errors |
| `cargo test --lib` | `3581 passed; 0 failed; 11 ignored`, 469.33 s |
| `cargo test --all-targets` | every test target ok (3,581 library + 61 integration), then `benches/psp_writer_perf.rs:386` panics — `index out of bounds: the len is 3300000 but the index is 3300000` |
| `cargo clippy --all-targets --all-features -- -D warnings` | red |

**Two gates are red and both are red identically on the parent branch `ng-joint-fit` before this commit** — checked, not assumed. Clippy: `this function has too many arguments (9/7)` and `useless use of vec!` twice in `lib test`, plus errors in `examples/ng_duplicated_class_harness.rs`. The bench: the same panic at the same line. Neither is in a line this commit changed. `cargo doc --no-deps` and `cargo audit` were not run.

Findings labelled "Needs verification": **0**. Every finding below was demonstrated by a probe that fails one way and passes the other.

### 4. Open questions and assumptions

1. **Are flank bases in scope for the STR difference list?** `TractDifference::offset`'s doc promises negative offsets in the left flank and offsets past the tract in the right; no code path can produce either, because `add_ssr` reads offsets from a `zip` over the tract alone. Either the writer must compare against the flanks, or the doc must narrow to `0..len`. Affects **B2**. This is a design question and is raised with the owner rather than answered here.
2. **Is the difference list read by anything today?** If it is not, B1's fix moves no fitted number and can land in any step; if it is, it must not land inside milestone A. Affects **B1**.

### 5. Top 3 priorities

1. **B1** — two distinct reads at one repeat tract come back as read 0 twice, which is exactly the confusion the read index exists to prevent. Pre-existing; demonstrated by a test that fails against unmodified code.
2. **M1/M2** — the manual `PartialEq` and `first_disagreement` field lists. Milestone B adds fields to these very structs, and a field that drops out of the comparison pools incomparable evidence without a word. **Fixed in this commit.**
3. **M3/M4** — the `kept_loci` refusal and `thin_to_cap` both mutated freely with the whole suite green. `thin_to_cap` is the one function whose stated contract is the byte-identity this milestone is measured by. **Fixed in this commit.**

### 6. Findings

#### Blocker

**B1: src/ng/parameter_estimation/joint/census.rs:1137 — a repeat tract's difference reads are numbered per observation, so two reads become one (pre-existing)**
**Categories:** reliability. **Confidence:** High — a test written by the reviewer fails against unmodified code.
`add_ssr` numbers each difference's read with `for copy in 0..reads`, where `reads` is *one observation's* count. The counter restarts at every observation, so two distinct read sequences each carrying one interruption both come back as read 0. `TractDifference::read`'s own doc promises "which of this locus's reads carried it, in the locus's own read order". A consumer grouping by read then sees one read with two interruptions — evidence of a bad read — where the truth is one interruption on two reads, which is evidence of an allele. The existing fixture hands the writer a single observation with a count of two, the one regime where the two numberings coincide.
**Not caused by this commit** — `add_ssr`'s body is byte-identical across the rename. **Not fixed here:** it changes what the census records, which milestone A may not do. Raised at Checkpoint A.

**B2: src/ng/parameter_estimation/joint/census.rs:1419 — `a_flank_difference_and_an_interior_one_come_back_apart` asserts only over its own literals, and the property has no producer (pre-existing)**
**Categories:** reliability. **Confidence:** High.
The test builds two `TractDifference` values and asserts `-2 < 0` and `4 ∈ 0..12`. It calls nothing. Worse, the writer cannot emit a negative offset at all, so the doc-comment's flank range describes a state no code path produces. Spec §7.3 names flank-against-interior as the assertion an interrupted-repeat model rests on, and the suite reports it green. **Not fixed here** — the repair needs open question 1 answered first.

#### Major

**M1: src/ng/parameter_estimation/joint/loci.rs:457,513 — manual `PartialEq` impls list fields by hand — FIXED**
**Categories:** refactor_safety. **Confidence:** High — mutation: dropping `self.scan == other.scan` left `80 passed; 0 failed`.
Both impls spelled out a field list, so a field added to `CatalogBuildSettings` or `SelectionTerms` compiles clean and is simply not compared. These are the values that let the fit **refuse** to pool two runs walked on different machines. **Milestone B adds encoding parameters to exactly these terms.** Fixed by destructuring `Self` without `..`, so the next added field stops the function compiling.

**M2: src/ng/parameter_estimation/joint/census.rs:599 and loci.rs:531 — `first_disagreement` walks a hand-written field list — FIXED**
**Categories:** refactor_safety. **Confidence:** High.
Same failure on the path that produces the operator-facing message. Fixed the same way. The doc comments' hard-coded counts ("thirteen values", "the seven") were replaced with phrasing that cannot go stale.

**M3: src/ng/parameter_estimation/joint/census.rs:603 — the `kept_loci` refusal had no test and the mutation survived — FIXED**
**Categories:** refactor_safety, reliability (convergent). **Confidence:** High — disabling the branch left `80 passed; 0 failed`.
`CensusLociDigest` is the witness that the selection rule really produced the same list of positions on both machines — the check that makes *"the same questions put to every sample"* true. Every fixture in the module built its terms with the same empty digest, so no test ever held two that differed. Fixed by `a_different_set_of_kept_loci_is_refused_and_named`.

**M4: src/ng/parameter_estimation/joint/census.rs:871 — `thin_to_cap` had no test, and its contract is this milestone's oracle — FIXED**
**Categories:** refactor_safety, reliability (convergent). **Confidence:** High — `.round()` → `.floor()` left all 25 census tests green while `thin_to_cap(3, 8, 5)` moved from 2 to 1.
Fixed by `a_thinned_share_rounds_to_nearest_and_never_loses_the_last_read`, a table over every branch: under the cap, zero depth, a share that rounds up, one that rounds down, the cap clamp, and the floor that keeps one stray read at 400 reads a position. **One number in the reviewer's proposed test was wrong and was corrected before use:** it asserted `thin_to_cap(20, 40, 5) == 5`, where 20×5/40 is 2.5 and rounds to 3.

**M5: src/ng/parameter_estimation/joint/census.rs:1339 — `read_groups_fold_by_addition` queried only the index carrying data — FIXED**
**Categories:** reliability. **Confidence:** High.
An implementation that ignored the index entirely returned the same answer. Fixed by asserting the empty position too.

**M6: src/ng/parameter_estimation/joint/census.rs:1359 — the four STR states are asserted against hand-set private fields, so the writer is never shown to produce any of them (pre-existing, deferred)**
**Categories:** reliability. **Confidence:** High.
Two writer branches are consequently dead to the suite: the covering-not-crossing path (no fixture builds a `ReadWitness::Partial`) and the reads-without-observation path (both locus helpers set that field to zero). Following the branches by hand shows the gap is not benign: a locus whose only read trips the guard reads back as `NoRead` — "walked, and no read reached the locus" — for a locus reads did reach and did cross. Deferred: writing it needs new fixture machinery, and B4 rewrites this record's shape.

**M7: src/ng/parameter_estimation/joint/census.rs:270 — `GenericEvidence::from_parts` documents two panics and no test provokes either (pre-existing, deferred)**
**Categories:** reliability. **Confidence:** High.

**M8: src/ng/parameter_estimation/joint/census.rs:973 — `add_generic`'s multi-base branch and its `ReadWitness::Partial` handling are never executed (pre-existing, deferred)**
**Categories:** reliability. **Confidence:** High. Wide loci are the ordinary case for an indel-bearing generic locus, and the `Partial` arm's clamps have no fixture.

**M9: src/ng/parameter_estimation/joint/census.rs:1376 — `the_difference_list_tells_one_interruption_on_two_reads_from_two_errors` also asserts only over its own literals (pre-existing, deferred)**
**Categories:** reliability. **Confidence:** High. Filed below Blocker only because the property is genuinely covered by a writer-driven test beside it. Tied to B1's and B2's repair.

#### Minor

**Mi1: five prose sites the `sed` pass damaged or left behind — FIXED**
**Categories:** naming, module_structure, refactor_safety (convergent).
`fit.rs:4` displayed the retired module name in a live intra-doc link; `loci.rs:465` had been rewritten into a claim untrue of the function it documents; `loci.rs:526,529` still said "identities" on the very method the rename exists to fix; `loci.rs:835` gained a doubled "the"; and `census.rs`'s seven section banners still divided the file by "record" while its own module doc had been rewritten to say "evidence". Also fixed: the test name `an_identity_equals_itself…`, the helper `generic_records()`, `loci.rs:638`'s "the record writer", `ssr_fit.rs:1198`'s hyphenation, and `doc/devel/ng/README.md:150`'s stale `SelectionIdentity`.

**Mi2: src/ng/parameter_estimation/joint/census.rs:1307 — a comment claimed a corner the assertion did not check — FIXED**
**Categories:** reliability. The comment said "reads, none non-reference" above an assertion on the never-walked sentinel. The fixture now carries a fifth position and asserts the real corner.

**Mi3: src/ng/parameter_estimation/joint/census.rs:1046 — length check then literal index — FIXED**
**Categories:** refactor_safety. Replaced by a slice pattern, which is what the compiler can check.

**Mi4: src/ng/parameter_estimation/joint/census.rs:1438 — the guard-threshold fixture leaves three of the function's decisions unexercised (pre-existing, deferred)**
**Categories:** reliability. The locus filter, the offset-0 exclusion and the exact one-in-ten boundary all mutate freely.

**Mi5: src/ng/parameter_estimation/joint/census.rs:1617 — the writer is only ever exercised on one contig (pre-existing, deferred)**
**Categories:** reliability. A comparison that dropped the contig would record another chromosome's reads at a position and never panic.

**Mi6: src/ng/parameter_estimation/joint/census.rs:1563 — `summing_windows_weights_each_by_its_own_position_count` passes under an unweighted mean (pre-existing)**
**Categories:** reliability. **Not fixed, and it will not need fixing: step A3 deletes `CoverageByWindow` entirely.** Recorded so the defect is not carried forward if any part survives.

#### Nits

`fit.rs:3165` reaches a sibling through `super::super::` inside its test module (pre-existing). The sibling import of `census` is written as `super::census` in two files and crate-absolute in four (pre-existing; both legal). `census.rs:1220`'s `terms: RecordingTerms { selection: self.terms, … }` puts two different types under one word. `JointFitError::IdentityMismatch` is now the only name an operator reads that still says "identity" — left because `arch/parameter_prepass_joint_loci.md` specifies that variant name, so the design has to move first. `contamination.rs:813` and `ssr_fit.rs:814` hold local `terms` arrays of log-likelihood summands in files where `terms` now also names the recording terms.

### 7. Out of scope observations

Three documentation-accuracy items, and **the first two are in files this run may not edit** — they are raised with the owner at Checkpoint A instead.

- `doc/devel/ng/arch/parameter_prepass_joint_records.md:568-585` — the rename's own authority still reads *"The code has not been renamed yet, and this document leads it"*, its table is headed *"the code today"*, and *"It waits rather than colliding … Do it in one pass when that work is committed"* describes work that is now done. Whoever reads it next for milestone B is told the code carries names it no longer has.
- `doc/devel/ng/arch/module_layout.md:73-74` — step 4's sub-units are listed as `fitting/`, `generic/`, `ssr/`. The tree holds a fourth, `joint/`, with 9,875 lines across six files. Pre-existing and general: the same doc has no `repeat_catalog/` entry either.
- The four `examples/ng_joint_*records*` files keep the old vocabulary in their filenames while their contents say `census`. Left alone because the impl plan's own preconditions name `examples/ng_joint_records_walk.rs`.
- About forty local bindings in `census.rs`, `fit.rs` and `ssr_fit.rs` are still named `records` where they hold a `GenericEvidence` or an `SsrEvidence`. Both reviewers raised it and both routed it away from this commit; it belongs with whoever schedules milestone B's sweep.

### 8. Missing tests to add now

Added in this commit:

- `a_different_set_of_kept_loci_is_refused_and_named` — input class: two sets of recording terms differing only in the digest of the loci actually kept. Catches the `kept_loci` branch being disabled, which pools samples indexed against different position lists.
- `a_thinned_share_rounds_to_nearest_and_never_loses_the_last_read` — input class: `thin_to_cap` across every branch it has. Catches rounding downwards, the cap clamp being dropped, and a single stray read at a very deep position being rounded away.

Named, specified and **deferred**, each to the step that changes the code it covers:

- `two_reads_carrying_different_interruptions_are_two_reads_and_not_one` — B1's repair.
- `a_position_above_the_cap_records_the_capped_depth_and_thinned_counts` — **milestone B, steps B2 and B3**, which change exactly this behaviour and must assert it.
- `the_four_str_states_come_out_of_the_writer` — **step B4**, which rewrites this record.
- `a_read_that_slipped_by_part_of_a_unit_lands_in_the_guard_with_its_size`, `the_guard_threshold_counts_only_this_locus_and_only_the_differing_reads`, `a_kept_position_on_the_edge_of_a_walked_stretch_is_marked`, `a_locus_with_both_crossing_and_covering_reads_reads_as_crossed`, `from_parts_refuses_unsorted_entries`, `from_parts_refuses_an_entry_past_the_end`, `a_wide_locus_gives_every_position_its_depth_and_no_allele`, `a_partial_read_adds_depth_only_where_it_covered`, `an_offset_bucket_saturates_rather_than_wrapping_at_u16_max` — bodies are in [tmp/review_2026-08-14_census-a2/reliability.md](../../../../tmp/review_2026-08-14_census-a2/reliability.md).

### 9. What's good

- The rename was verified behaviour-preserving **mechanically rather than by reading**: both sides normalised, the substitutions applied to the old side, and the residual diff enumerated — 27 hunks, 14 prose-only, 13 non-comment, none of them an expression. That method is worth reusing on any future rename.
- `git mv` was used, so `git show --find-renames` reports `R095` and `git log --follow` keeps the file's history.
- `same_criteria` in [loci.rs](../../../../src/ng/parameter_estimation/joint/loci.rs) already solved the "a new field silently drops out" hazard the safe way, with a comment saying so — which is what made M1 obvious once someone looked at its neighbours.
- The five-bit packing round trip is asserted at **every bit offset within a byte**, by repeating the ladder eight times, rather than at one convenient index.

### 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo check --all-targets`
- `./scripts/dev.sh cargo test --lib` — expect 3,583 after this commit's two new tests
- `./scripts/dev.sh cargo test --lib ng::parameter_estimation::joint::census` — 27 tests
