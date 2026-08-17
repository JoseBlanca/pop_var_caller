# Fixes applied — ng cohort merge, B3

*2026-08-17, against [the review](ng_cohort_merge_b3_2026-08-17.md). All in
`src/ng/run/cohort_merge/build.rs`.*

## Behaviour changed

- **A sequence no read is behind contributes no allele, on either branch.** It used to be
  interned by the one-record branch — an allele nothing supports, carrying a quality nobody
  measured — while the composed branch panicked on it. Both now skip it, and
  `share_of_one_read`'s check is a documented backstop carrying the "must become a `RunError`
  on the psp path" obligation its three siblings carry (arch §5).
- **Two silent fallbacks removed.** The quality no longer falls back to `0.0` when a read has
  no sightings — which cannot happen, and `ln P(error) = 0` is the *worst* quality
  expressible, not a neutral one — and the division is no longer guarded by a condition the
  assertion above already rules out.
- **The records travel inside `AlleleBacking::OneRead`.** A sighting is a pair of indices and
  means nothing without the slice it indexes; the caller used to supply its own, right only
  by coincidence of expression.

## Claims corrected

- **`placed_left` is counted against each record's own position**, not the locus's anchor,
  and pooling it across records mixes as many questions as the sample has records — an
  approximation bounded by the locus width. The field said the locus's anchor and the
  division's justification covered only the two quantities that *are* read-invariant.
- **`round_to_u32`'s doc was wrong about the language.** Rust's float-to-integer `as` has
  saturated since 1.45, so these helpers are not a repair; what they buy is a visible
  boundary and a stated answer for `NaN`. Pinned by
  `the_rounding_of_a_divided_count_saturates_at_both_ends`.
- **The zero-read guard's comment told the wrong story.** It said the division "would come
  back as an infinity and poison every score, silently". Measured by removing the guard: the
  `-inf` is *discarded* by the `max` and the allele reports a plausible quality over reads
  that were never counted — worse than the story, and now the story.
- **`num_reads` was labelled "Exact." while built with `saturating_add`**;
  `reads_composed_across_records` did not say it saturates; `reads_removed_as_evidence` named
  one of the two removal cases.
- **One row per allele pools the read groups**, ending a cross `SequenceObservation` keeps
  deliberately. Recorded at the type with what wants it and the shape that would restore it.
- **The dense support row's cost is recorded at the field with its measurement** — 614 MB for
  one observation at 4,000 samples each showing a distinct allele — and raised at Checkpoint
  B rather than decided.

## Tests

97 in the module, from 94 at review time; three added.

- `one_observations_reads_split_across_two_alleles_and_each_takes_its_own_share` closes both
  Blockers with one fixture: a record with two sequences, three reads taking two paths, the
  weakest sighting **first** for two of them and **last** for the third, and a second sample
  whose own removal count differs from the cohort's. Both mutations were re-run against it:
  keeping the last sighting's quality **fails** it, and reading the record's first sequence
  instead of the sighted one **fails** it. Before, each left all 94 tests green.
- `the_rounding_of_a_divided_count_saturates_at_both_ends` — `NaN`, both infinities, negative,
  both boundaries, and the half-way case.
- `a_sequence_with_no_reads_behind_it_contributes_no_allele_on_either_branch`.

## Shape, not behaviour

`SupportSums` replaces three copies of the five-field list, with `add` and `divide_counts_by`
destructured so a sixth sum is a compile error rather than a silent zero; `ShownBy` became
`AlleleBacking::{OneSequence, OneRead}`; the per-allele tally is hoisted and refilled per
sample; a dead `resize` is gone; `CohortObservation::over` destructures the table.

## Not applied

- **The dense-row memory shape** (sparse row, allele cap, or price it in spec §8) — the
  owner's, at Checkpoint B.
- **Renaming the `AlleleSupport`/`SampleSupport` pair** — `SampleSupport` is arch §4's name.
- **A private field to protect `CohortObservation`'s parallel-vector invariant** — it changes
  the shape arch §4 declares.
- **Deleting `round_to_u32`/`round_to_u64`** — they now say what they do and are pinned.

## Validation

In the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — 97 passed, 0 failed.
- `cargo test --lib` — 3,720 passed, 0 failed, 11 ignored (3,717 at review time).
</content>
