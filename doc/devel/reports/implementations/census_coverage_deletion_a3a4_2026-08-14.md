# A3 and A4 — the coverage-by-window summary is gone, and nothing moved

**Plan:** [census_rename_and_encoding.md](../../ng/impl_plan/census_rename_and_encoding.md), milestone A,
steps A3 and A4.
**Design authority:** [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§4 (why it went and what replaces it) and [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§1.6.
**Date:** 2026-08-14.

---

## 1. What was deleted

`coverage.rs` — 566 lines that walked the reference and the loci in parallel to build one mean
read depth per 500 bp window, plus a depth-against-GC curve — and with it:

| gone | where |
|---|---|
| `CoverageAccumulator`, `CoverageGrid` | `joint/coverage.rs`, whole file |
| `CoverageByWindow`, `WINDOW_DEPTH_SCALE`, `MIN_ALIGNED_BASES_PER_WINDOW` | `joint/census.rs` |
| `SampleCensusEvidence::coverage` | `joint/census.rs` |
| `RecordingTerms::coverage_window`, and its branch of `first_disagreement` | `joint/census.rs` — the terms are now **twelve**, not thirteen |
| `CensusWriter::finish`'s parameter | `joint/census.rs` — `finish()` takes nothing |
| the `COVERAGE_DISCRIMINATOR` path and `coverage_odds()` | `examples/ng_joint_records_walk.rs` |

Measured over `src/` and `examples/`: **+43, −996, so 953 lines net**, of which 566 are the module
itself (`git diff --numstat HEAD -- src examples`).

**What was deliberately kept.** `JointFitConfig::coverage_odds` stays. It is the fit's own input
for the duplicated-copy class, not the summary; the harness now always hands it an empty vector,
because the thing that used to fill it is gone. Spec §4 says what replaces it — the position's own
depth, which the census already carries — and that is milestone B's ladder extension, not this step.

## 2. A4 — the oracle, and the one thing it caught

**The plan asks for byte-identical output on two real cohorts. It is byte-identical on every
fitted number, and the only lines that differ are the deletion showing through.** Comparing the
pre-rename binary against the post-deletion one on the small tomato oracle (eight accessions from
2.4 to 30.6 reads a position, six of the bench BED's spans, a 60,000-position census), every
differing line is one of five shapes:

| differing line | why |
|---|---|
| `coverage grid  N windows of 500 bp over N generic bases` | deleted |
| `coverage  N MB over N windows of 500 bp`, once a sample | deleted |
| the `; median window depth …` clause on each walk line | deleted |
| `TOTAL  N MB held for this sample` | no longer adds the summary's bytes |
| `the identity check` → `the recording-terms check`, `thirteen values` → `twelve values` | A2's vocabulary, and A3's count |

**Not one fitted number moved**: no read group's error rate, no sample's heterozygosity or
homozygote excess, no frequency density, no contamination estimate, and no repeat-tract stratum row.

### 2.1 The oracle is not the 63-accession cohort, and that is a recorded deviation

The plan's A4 says *"re-run both cohorts"*. The full tomato cohort takes **over two and a half
hours** — the walk is 625 s of it and the two fits are the rest — which does not belong in a loop
whose steps are minutes of editing. **Owner's call, 2026-08-14: build a smaller oracle.**

What replaces it runs in **88 s** and is in `tmp/oracle_tomato.sh`. What it gives up is stated
rather than hidden: at eight samples the duplicated-copy class is below the twenty-five samples it
needs to be identified at all, and the contamination estimator draws its allele frequencies from
the panel it is given. **It is a valid identity check and not a valid reading of those two
numbers.** For a rename and a deletion, which are structural, that is the whole of what A4 needs.

### 2.2 The fit's trajectory depends on the container's CPU count — MEASURED

The comparison first came back with one line differing that should not have: the fit's
largest-move column at pass 1, `615262.063054` against `615262.050187`. It was not the edit.

- **Same binary, run twice at the same CPU count: byte-identical output.** So the harness is
  deterministic.
- **Same binary at 4 container CPUs and at 6: those two values exactly.** The parallel reduction
  sums the same terms in a different order, and the eighth significant digit moves.

The converged parameters agreed at both counts, so this is a property of the trajectory rather
than of the answer. But it means **a comparison run at a different CPU count is comparing two
things**, and that is a trap the next person would fall into. `tmp/run_oracle.sh` now pins
`DEV_CPUS=4` and says why.

## 3. Validation

Run in the dev container.

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | clean; the only warnings are pre-existing dead code in `examples/dhat_ng_merge.rs` |
| `cargo test --lib ng::parameter_estimation::joint::census` | `24 passed; 0 failed` — 27 before, and the three that went are the coverage summary's own |
| the small tomato oracle, before and after | §2's table |

## 4. Follow-ups

- **`each_of_the_four_unit_values_refuses_on_its_own`** is the renamed test: the coverage window
  was one of the five and is now one of four.
- The reviewed-and-deferred defects from A2 stand unchanged; they are listed in
  [the A2 fix report](census_rename_a2_fixes_2026-08-14.md).
