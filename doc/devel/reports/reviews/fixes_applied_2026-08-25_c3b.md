# Fixes applied — ng calling loop C3b review

**Review:** [`ng_calling_loop_c3b_2026-08-25.md`](ng_calling_loop_c3b_2026-08-25.md).
**Date:** 2026-08-25. **Everything the review raised was applied; nothing was deferred and nothing
was escalated as a design question.**

## What changed

**Six tests added, one fixture changed, one assertion pair added** — the whole of the reliability
half:

| finding | fix |
|---|---|
| B1 — no multi-sample genotype quality pinned | `two_samples_take_their_confidence_from_their_own_posterior_rows`: two samples 4 nats apart in opposite directions, both qualities pinned to measured values (11.553 and 17.697 Phred) |
| B2 — the ninth count not tied to the primary alternative | two assertions added to the existing three-allele fixture: both calls are `2/2`, and the expectation is 17.0 where copies of allele 1 would give 0 |
| M1 — the reference row's two counts never told apart | the hand-computed fixture's reference row becomes 3 reads, 2 forward, 1 placed left; five assertions and the doc comment's four numbers updated |
| M2 — the set-aside sample's exclusion from the *choice* untested | `the_primary_alternative_ignores_a_set_aside_samples_reads`, at three alleles |
| M3 — the tract path's site quality unasserted | compared against an independently prepared fold, as the SNP/indel test does |
| M4 — four pass-through values unpinned | `the_pass_carries_the_region_the_outcome_and_the_provenance_onto_the_record` |
| Minor — no boundary test on the site-quality ceiling | `a_site_quality_at_the_ceiling_is_accepted` |
| reliability's own suggestions | `at_ploidy_one_a_call_expects_all_or_none_of_its_reads` and a proptest on `mint_genotype` over the whole `(ploidy, allele count)` grid |

**Five craft fixes**, all of which the craft agent had already compiled in its own worktree:

- `pooled_primary_alternative` → `pool_reads_and_pick_primary_alternative`;
- three `// PANIC-FREE:` markers, naming `MAX_ALLELE_COUNT` and `prepare_for_locus` rather than
  restating the invariant;
- `add_called_sample` takes `Ploidy` instead of `f64`, and the value stays a `Ploidy` from the
  genotype table to the division it is the divisor of;
- the two co-dependent `Option`s become one `artifact_pool` value bound by `and_then`;
- `LocusEvidence::Generic { region: _, per_sample }` names the field it discards;
- the second `use` block inside the test module folded into the first, and `repeat_n` imported.

**Ten corrected claims**, in doc comments and in the implementation report — three mechanisms and
seven numbers. The three mechanisms are listed as M5–M7 in the review; the most consequential is
that **the plan's blocker note about the repeat-tract read-likelihood row is stale**, and this
commit corrects the plan as well as the comment that repeated it.

## What was not changed, and why

- **Three cross-category observations were left alone**: release-profile wrapping on `u64` read
  counters (needs ~4×10⁹ observations at one locus), the site quality's quadratic fold (a
  bakeoffs-plan measurement, `calling_quality.md` §13's Q3), and `ArtifactTestCounts` deriving
  `Copy` at ~72 bytes (the previous commit's).
- **The three judgement calls the review was asked to examine were all upheld** — the `Option`, the
  private quality field with its `pub(crate)` reader, and the nine arguments. Each is argued in the
  review; none needed a change.

## Validation after the fixes

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4672 passed; 0 failed; 14 ignored` (4,644 before C3b).
- `cargo test --release --lib ng::calling --all-features` — `626 passed; 0 failed; 3 ignored`.
- The seven release-held checks downgraded to `debug_assert` in one run: `618 passed; 8 failed`,
  every check reached. Files restored from byte-identical backups and the release suite re-run green
  before the commit.
