# Code Review: str_partial_flank
**Date:** 2026-10-21
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** commit 2cc8ea86 on branch `str-partial-flank` — a partial read compared with an edge-deletion allele followed by the flank it leaves (`src/calling/likelihood/ssr_emission.rs`)
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the diff `b0162d47..2cc8ea86 -- src/`, one file. It replaces the prototype `censored_with_the_flank_beyond` (commit 707d725f, behind `PVC_EXP_PARTIAL_FLANK`) with the shipped rule: `remaining_flank_ahead`, `censored_into_the_remaining_flank`, `censored_from_the_tract_start_into_the_flank`, a new `then` parameter on `letters_over`, and `letters_of_stretching_then` in place of `prefix_of_stretching` + `substitution_probability` at three sites. Five new tests.
- **Reviewed against:** 2cc8ea86 (HEAD of `str-partial-flank`), compared with b0162d47.
- **In-scope files:** [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs); its caller [src/calling/likelihood/ssr.rs](../../../../src/calling/likelihood/ssr.rs) `fill_emissions` (unchanged, read for the half-and-half case); the contract, [doc/devel/ng/spec/str_edge_deletion_alleles.md](../../../ng/spec/str_edge_deletion_alleles.md) §5.6, from "Reads that end inside the segment" to "How we know it works".
- **Deliberately out of scope:** STR bundles (`bundle_censored_emission`, untouched); complete-read scoring in `fill_emissions` (untouched).
- **Categories applied:** reliability (correctness of the sum, edge cases, test challenge), errors (the too-short guard), naming, idiomatic, refactor_safety (the bit-for-bit claim), float_portability (log-space sums reach genotype likelihoods), smells, extras (diff matches stated intent). Skipped: unsafe_concurrency (none in scope), module_structure (one file), tooling (no Cargo change), defaults (no new default-acting value).
- **Deviation from the skill's procedure:** the categories were run by one reviewer, not fanned out to sub-agents in separate worktrees, because the caller restricted the work to this one worktree. Mutation testing was done on a copy of the crate under `tmp/review/crate/` (and the pre-change code under `tmp/review/crate_old/`), so no source file was edited. `PROJECT_STATUS.md` has no block for this feature and was not edited (the caller asked for the report only); add a block when the branch is merged.

## 2. Verdict

**Approve-with-changes.** The rule matches the spec, and every candidate it does not reach scores bit for bit as before — checked by running the pre-change and post-change code on the same 5,280 cases. The findings are documentation that the change made stale, one test whose name claims more than it checks, and three missing tests. None changes a result.

## 3. Execution status

Commands run (natively, `CARGO_TARGET_DIR` under `tmp/review/`):

| Command | Exit | Result |
|---|---|---|
| `cargo test --lib calling::likelihood` | 0 | `test result: ok. 286 passed; 0 failed; 0 ignored; 0 measured; 4981 filtered out` |
| `cargo test --lib` | 101 | `test result: FAILED. 5259 passed; 1 failed; 7 ignored` — the one failure, `parameter_estimation::joint::fit::tests::a_frequency_that_is_not_a_probability_does_not_reach_the_caller` (`panicked at src/genetics.rs:83:5: Beta shapes must be above zero, got 0 and 0`), fails identically at b0162d47: pre-existing, unrelated |
| `rustfmt --edition 2024 --check src/calling/likelihood/ssr_emission.rs` | 0 | no diff |
| `cargo fmt --check` | 1 | diffs only in `examples/dup_filter_allele_balance_grid.rs`, `src/paralog/locus_score.rs`, `src/paralog/mod.rs` — files this commit does not touch |
| `cargo clippy --lib --tests -- -D warnings` | 0 | clean |
| `cargo clippy --all-targets -- -D warnings` | 101 | `error: could not compile pop_var_caller (example "dhat_ng_psp") due to 14 previous errors` — an example this commit does not touch |

Probes and mutations (all in the scratch copies):

- **Old against new, bit for bit.** A grid of 4,800 partial reads (8 tracts — pure, interrupted, ending part-way through a unit, period 1 to 4 — 60 reads each, both held sides, five flank situations: no deletion with and without flanks, a deletion on the held side, and a deletion that consumes the whole 15-base stored flank or more) and 480 complete reads, scored by b0162d47 and by 2cc8ea86: **all 5,280 values have identical bits.**
- **Mutations against the five new tests:** 5 run, 4 killed, 1 survived, 0 that changed no behaviour. The survivor (dropping `!flank.is_empty()` from `remaining_flank_ahead`) changes 588 of the grid's 4,800 values by a few units in the last place: see Mi3.
- **A rounding mutation in `letters_of_stretching_then`** (`exp(total_ln · (1 + 4ε))`) is caught by two older tests (`a_censored_read_of_one_admissible_length_scores_what_the_complete_read_does`, `a_one_run_bundle_scores_as_the_segment_it_is`), not by the new test that claims to check "as before": see Mi2.

Not run: `cargo doc --no-deps`, `cargo audit` (no dependency or public-API change).

Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

None that block. One assumption, stated: a partial read's `bases` never exceed the reference locus length, so under a deletion candidate `ℓ ≤ |a| + deleted tract bases`; the spec's argument that "a stretching whose spelling plus F is shorter than the read does not arise" rests on it, and the code handles the case anyway (scores zero).

## 5. Top 3 priorities

1. **Mi1** — `letters_over`'s "One function, two questions" section now says something false about the path this commit adds.
2. **Mi2** — the test named "scores as before" does not compare with the code before; pin a few recorded values so the bit-for-bit promise has a guard of its own.
3. **Mi3** — the empty-remaining-flank fallback and a both-sides deletion have no test; the fallback's removal survives every test.

## 6. Findings

### Blocker

None.

### Major

None.

### Minor

#### Mi1: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L1117-L1125) — `letters_over`'s doc says a censored read is shorter than the stretching; on the new path it need not be
- **Confidence:** High
- **Problem:** The section reads: "A censored read is shorter, and that truncation is the whole difference between the two callers." Under `censored_into_the_remaining_flank` the floor admits stretchings *shorter* than the read (Δ with `|a| + Δ + |F| ≥ ℓ`), and the comparison runs past the stretching into `then`. The `# What follows the stretching` section added below it is correct; the paragraph above it now contradicts it.
- **Why it matters:** a reader who trusts the first section will believe `observed.len() ≤ stretching.len()` is an invariant of `letters_over` and may add an index or a slice that panics on the new path.
- **Suggested fix:**
  ```rust
  /// read is exactly as long as the stretched candidate — `bp_diff` was derived from it — so
  /// truncating to its length is a no-op there. A censored read compared with the tract alone is
  /// shorter, and that truncation is the whole difference between the two callers; one compared
  /// with the tract and the flank its candidate leaves may be longer than the stretching, and the
  /// comparison continues into `then` (below).
  ```

#### Mi2: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L6402) — `a_candidate_without_a_flank_deletion_ahead_of_the_read_scores_as_before` does not compare with "before"
- **Confidence:** High
- **Problem:** Both arms of the comparison run the post-commit code: "before" is the same candidate with no deletion and no flanks. The test therefore proves that `remaining_flank_ahead` returns `None` for these candidates (mutation M1, scoring every candidate with its flank, is killed by it and only by it), not that the refactor of `letters_over` left their scores unchanged. A rounding change in `letters_of_stretching_then` passes it: `exp(total_ln · (1 + 4ε))` leaves this test green and is caught only by two older cross-path tests.
- **Why it matters:** the commit's main safety claim — "all other candidates score bit for bit as before", which the spec's byte-identical oracle runs rest on — has no test of its own. It holds today (5,280 of 5,280 values identical, §3), but nothing pins it.
- **Suggested fix:** either rename to what it checks (`a_candidate_without_a_flank_deletion_ahead_of_the_read_is_not_reached_by_the_flank_rule`) and say so in the doc, or keep the name and add recorded bits from b0162d47 for three or four of its cases — one on the pure-tract factorised path (`CHR20_TRACT`, count 8), one on the exact sum (`CHR20_ALLELE`, count 6), one right-held — asserted with `to_bits()`. The values can be produced by running the case list at b0162d47 (`tmp/review/crate_old/` holds that tree with the probe).

#### Mi3: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L851-L858) — three behaviours of the rule have no test; one of them survives its own removal
- **Confidence:** High
- **Problem:**
  1. **The empty-remaining-flank fallback.** `remaining_flank_ahead` returns `None` when the deletion takes the whole stored flank or more, keeping the shipped scoring. Removing `&& !flank.is_empty()` leaves all five new tests green, and it does change output: 588 of the 4,800 grid values move by a few units in the last place (example: `GAGAGAGAGAGAGAGA`, left-held read `GAGACAG`, deletion 15 of 15 — `0x3f35b5866fcf4574` becomes `…4578`), because the full sum and the factorised form round differently on a pure tract.
  2. **A candidate deleting flank bases on both sides.** Each held side must read only its own side's deletion. A probe confirms it (`{left: 2, right: 5}` scores bit for bit as `{0, 5}` from the left and as `{2, 0}` from the right, read lengths 1 to 16), but no test states it.
  3. **The half-and-half case** in `fill_emissions` (`likelihood/ssr.rs:927`): a whole-locus witness with no recorded side, scored against a deletion candidate, takes the rule from the left half only when the deletion is on the right. Unchanged code, but it is the one place both of the new function's branches meet in one emission.
- **Why it matters:** item 1 is the only branch that keeps the bit-for-bit promise for deletions that run past the stored flank; a later edit that "simplifies" the condition would move bytes on the oracle with no test failing.
- **Suggested fix:** add the three tests listed under §8.

#### Mi4: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L6252-L6266) — the first new test's doc attributes the fixture's number to the real site, and understates which expansions read as the reference
- **Confidence:** High
- **Problem:** Three things in one doc comment.
  - "1.2 × 10⁻³ against the allele's tract alone — the score that wrote the homozygous site `0/1`." The 1.2 × 10⁻³ is this fixture's, under `a_model()`; the spec and `censored_into_the_remaining_flank`'s own doc give the real site's score as 0.003 under the fitted model. A reader comparing the two docs sees two different "scores that wrote the site".
  - "the unchanged allele and a one-unit expansion both read as the reference does". Every whole-unit expansion does: `GA` × (6 + k) + `G`, then `AGACAG…`, spells `GA` × 8 over the read's 16 bases for every k from 0 to 10. Only the part-unit changes and the contractions differ.
  - A broken wrap: "The allele's score is the full" ends a line and "sum, written out by hand" starts the next one.
- **Why it matters:** the doc is the explanation of why the allele scores as the reference does; it names the wrong set of terms.
- **Suggested fix:**
  ```rust
  /// Measured: 0.96686 under the allele, 0.96697 under the reference, 1.2 × 10⁻³ against the
  /// allele's tract alone — this fixture's version of the slippage-only score that wrote the
  /// homozygous site `0/1`. The allele's score is the full sum, written out by hand: the
  /// unchanged allele and every whole-unit expansion read as the reference does over the read's
  /// sixteen bases, and no factorised shortcut stands in for the sum.
  ```

## Nits

- `letters_over(…, then: &[u8], …)` and `letters_of_stretching_then`: `then` is a conjunction, not a name for the value. `followed_by` (or `flank_after`) says what it is; the function could be `letters_against_stretching_followed_by`.
- [ssr_emission.rs:1224](../../../../src/calling/likelihood/ssr_emission.rs#L1224): the `debug_assert!` message still says "stretchings at least as long as what the read showed"; with `then` it is "stretchings that, with what follows, are at least as long …", as the doc above it already says.
- [ssr_emission.rs:6454](../../../../src/calling/likelihood/ssr_emission.rs#L6454): "as the prototype did" names code that no longer exists in the tree. Say "as a rule that stretches the tract at its end would".
- [ssr_emission.rs:1229](../../../../src/calling/likelihood/ssr_emission.rs#L1229): `letters_of_stretching_then` looks up `substitution_scores(context)` (two logarithms, per its own doc) on every call — once per placement per length change on the new path, which on a pure tract now walks the whole support instead of one factorised comparison. Inherited from `substitution_probability`, so not a regression; hoisting the lookup into the two censored sums would remove it.
- [ssr_emission.rs:949](../../../../src/calling/likelihood/ssr_emission.rs#L949): 104 characters; the rest of the file wraps doc comments at 100.
- `censored_into_the_remaining_flank`'s doc uses "stretching", "rendering" and "placement" for the candidate after a slip, within three paragraphs. One term would read more easily.

## 7. Out of scope observations

- `parameter_estimation::joint::fit::tests::a_frequency_that_is_not_a_probability_does_not_reach_the_caller` fails at b0162d47 and at HEAD (`Beta shapes must be above zero, got 0 and 0`, `src/genetics.rs:83`). Separate fix.
- `cargo fmt --check` reports diffs in `src/paralog/locus_score.rs`, `src/paralog/mod.rs` and `examples/dup_filter_allele_balance_grid.rs`; `cargo clippy --all-targets` fails on the example `dhat_ng_psp` (14 errors, e.g. `let_unit_value`). Neither is touched by this commit.

## 8. Missing tests to add now

For `remaining_flank_ahead` / `censored_emission`:

- **`a_deletion_that_takes_the_whole_stored_flank_scores_as_the_tract_alone`** — input class: right deletion ≥ the stored right flank (15 of 15, and 20 of 15), left-held reads of several lengths on a pure tract. Assert `to_bits()` equal to the same candidate with no flanks. Catches: dropping `!flank.is_empty()` (verified: 588 of 4,800 grid values move).
- **`a_candidate_deleting_both_flanks_reads_each_side_from_its_own_deletion`** — `{left: 2, right: 5}` against `{0, 5}` held left and `{2, 0}` held right, `to_bits()` equal. Catches: reading the deletion of the wrong side, or the sum of both.
- **`a_whole_locus_witness_with_no_side_takes_the_flank_rule_from_one_half_only`** (in `likelihood/ssr.rs`) — a pre-2.1 whole-locus partial witness against a right-deletion candidate; expected `0.5 · censored(Left) + 0.5 · censored(Right)`, where only the left half reaches the flank. Catches: a future change to `fill_emissions` that applies one side's rule to both halves.

For the bit-for-bit promise: see Mi2 (recorded bits from b0162d47).

## 9. What's good

- `remaining_flank_ahead` reuses `flank_beyond_the_hull`, so a segment's partial reads and a bundle's see the same flank under one definition, as the spec requires ([ssr_emission.rs:851](../../../../src/calling/likelihood/ssr_emission.rs#L851)).
- The new path routes through the existing `letters_over`, so slips are placed per run exactly as for every other partial read; test 4 proves the difference from the prototype's stretch-at-the-end with a hand sum (1.1 × 10⁻³ against 8 × 10⁻¹⁴, both reproduced).
- The right-held case reuses the `MirroredCensoredRead` scratch, extended by one buffer, so the new path allocates nothing per call.
- The mirror test asserts `to_bits()` equality between a fixture and its mirror image, which is what catches an unreversed flank (mutation M4 is killed by it and only by it).
- Every number quoted in the new doc comments was reproduced: 0.96686 / 0.96697 / 1.2 × 10⁻³ (measured 0.9668622, 0.9669661, 1.229560 × 10⁻³); 1.1 × 10⁻³ / 8 × 10⁻¹⁴ (1.082026 × 10⁻³, 8.039494 × 10⁻¹⁴); "about half the one-unit expansion's probability" (0.5 × 2.21 × 10⁻³ × 0.9792 = 1.082 × 10⁻³); "`GCA`, since the tract is fourteen bases long" (`CAGCAGTTCAGCAG` + 3 → `…CAGGCA`); "leaves two flank bases" (`TTGCAC`[4..] = `AC`).

## 10. Commands to re-verify

- `cargo test --lib calling::likelihood`
- `cargo clippy --lib --tests -- -D warnings`
- `rustfmt --edition 2024 --check src/calling/likelihood/ssr_emission.rs`
- Old-against-new bits (review probe, in the scratch copies):
  `cd tmp/review/crate && REVIEW_PROBE_OUT=$PWD/../bits_new.txt CARGO_TARGET_DIR=$PWD/../target cargo test --lib review_probe_bits_grid`,
  the same in `tmp/review/crate_old` with `CARGO_TARGET_DIR=$PWD/../target_old` (a separate target directory is required: sharing one let cargo reuse the other tree's test binary), then `diff tmp/review/bits_old.txt tmp/review/bits_new.txt`.
- Mutations: `tmp/review/mutate.sh <name> <perl substitution>` (runs the five new tests against one mutant).

### Author response convention

Address each finding by its identifier (Mi1–Mi4) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`.

## 11. Fixes applied (2026-10-21)

All four findings and four of the six nits were applied in the commit after 2cc8ea86. None changes a
score: every change is a name, a doc comment or a test.

- **Mi1.** `letters_over`'s "One function, two questions" now says a read compared with the tract and
  the remaining flank may be longer than the stretching, and the comparison continues into what
  follows it.
- **Mi2.** The test is renamed `a_candidate_without_a_flank_deletion_ahead_of_the_read_is_not_reached_by_the_rule`,
  and a new test, `a_partial_read_off_the_rule_scores_the_bits_recorded_before_it`, pins four cases
  (pure tract, interrupted tract, a tract ending part-way through a unit, and that one from the
  right) to bit patterns recorded by running the same cases on b0162d47.
- **Mi3.** Three tests added: `a_deletion_that_takes_the_whole_stored_flank_scores_as_the_tract_alone`
  (on a pure tract, so the tract-alone route is the factorised one; dropping `!flank.is_empty()` now
  fails it), `a_candidate_deleting_both_flanks_reads_each_side_from_its_own_deletion` (reading the
  near side's deletion instead fails it and five others), and, in `likelihood/ssr.rs`,
  `a_whole_tract_read_with_no_side_takes_each_sides_own_flank_rule` (the row's cell is half of each
  side's score, to the bit, and the left half is more than a hundred times the right).
- **Mi4.** The doc now says the 1.2 × 10⁻³ is under the test's model (the site's own was 0.003 a
  read), that every whole-unit expansion reads as the reference does over the read's sixteen bases,
  and the broken wrap is gone.
- **Nits applied:** `then` → `followed_by` and `letters_of_stretching_then` →
  `letters_against_stretching_followed_by`; the `debug_assert` message names what follows the
  stretching; "as the prototype did" → "a rule that stretched the tract at its end"; the over-long
  doc line is rewrapped.
- **Nits not applied:** hoisting `substitution_scores` out of the per-comparison call (inherited by
  every partial read, not introduced here; a performance change belongs with a measurement), and the
  "stretching / rendering / placement" wording, which follows the three words' existing use in this
  file — a placement is where a slip lands, a rendering the bases that gives, a stretching the
  candidate after any length change.

`cargo test --lib calling::likelihood`: 290 passed. `cargo clippy --all-targets --all-features -- -D
warnings`: clean.
