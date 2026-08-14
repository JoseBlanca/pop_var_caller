# A2 — fixes applied from the review

**Review:** [census_rename_a2_2026-08-14.md](../reviews/census_rename_a2_2026-08-14.md).
**Commit:** amended into `refactor(ng): A2 …`, so the step stays one commit.
**Date:** 2026-08-14.

---

## What was applied, and why each one was in scope

Milestone A's oracle is that two real cohorts return the same fitted numbers, so the bar for
applying a finding inside this step is **it must not change what the walk records**. Everything below
clears that bar: it is prose, an exhaustive destructure that compares the same fields, a slice pattern
in place of a length test, or a new test.

| finding | what was done |
|---|---|
| **M1** — `impl PartialEq` for `CatalogBuildSettings` and `SelectionTerms` listed fields by hand | both destructure `Self` without `..`, so a field added in milestone B stops them compiling instead of dropping out of the comparison |
| **M2** — both `first_disagreement` functions did the same | same fix; the doc comments' hard-coded counts ("thirteen values", "the seven") replaced with phrasing that cannot go stale |
| **M3** — the `kept_loci` refusal had no test | `a_different_set_of_kept_loci_is_refused_and_named` |
| **M4** — `thin_to_cap` had no test | `a_thinned_share_rounds_to_nearest_and_never_loses_the_last_read`, a table over all six branches |
| **M5** — `read_groups_fold_by_addition` queried only the index carrying data | asserts the empty position too |
| **Mi1** — prose the `sed` pass damaged or left behind | eleven sites, listed in the review |
| **Mi2** — a comment claiming a corner its assertion did not check | the fixture grew a fifth position and the assertion now checks it |
| **Mi3** — length check then literal index | slice pattern |

**One number in a reviewer's proposed test was wrong and was corrected before it was written down.**
The suggested table asserted `thin_to_cap(20, 40, 5) == 5`; 20 × 5 / 40 is 2.5, which rounds to 3, and
the clamp it was meant to demonstrate needs `thin_to_cap(400, 400, 20) == 20` instead. Every assertion
in the committed test was run before this sentence was written.

## What was not applied, and why

**Two Blockers, both pre-existing, both left standing deliberately.**

- **A repeat tract's difference reads are numbered per observation** (`add_ssr`), so two distinct reads
  each carrying one interruption come back as read 0 twice. Fixing it changes what the census records,
  which milestone A may not do, and it is not what any step of this plan is for. **Raised at
  Checkpoint A.**
- **The flank-against-interior test asserts only its own literals, and the writer cannot emit a flank
  offset at all.** The repair needs a design answer — are flank bases in scope for the difference
  list? — so it is a question for the owner, not a patch. **Raised at Checkpoint A.**

**Four Majors and two Minors deferred to the step that changes the code they cover**: the four STR
states through the writer (B4 rewrites that record), the depth cap at the writer (B2 and B3 change
exactly that behaviour and must assert it), `from_parts`'s two documented panics, `add_generic`'s
wide-locus and partial-read branches, the guard threshold's three unexercised decisions, and the
writer's single contig. Named with bodies in the review's §8.

**One Minor needs no fix at all:** the coverage-summary fixture that passes under an unweighted mean —
step A3 deletes the type.

**Two documentation items are in files this run may not edit** — `arch/parameter_prepass_joint_records.md`
§4, which still says the rename has not happened, and `arch/module_layout.md`, which lists three
sub-units under step 4 where the tree holds four. Raised at Checkpoint A.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib ng::parameter_estimation::joint::census` | `27 passed; 0 failed` (25 before) |
| `cargo test --lib ng::parameter_estimation::joint::` | `82 passed; 0 failed` (80 before) |
| `cargo test --lib` | see the amended commit message |
