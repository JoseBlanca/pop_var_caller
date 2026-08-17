# Fixes applied — ng cohort merge, C1

*2026-08-17, against [the review](ng_cohort_merge_c1_2026-08-17.md). All in
`src/ng/run/cohort_merge/build.rs`.*

## Tests, and the mutation each one kills

Three added (104 → 111 in the module, with the review's naming pass and the earlier B-step
tests included in that total). Each was run against the mutation it exists to kill, and each
mutation now fails:

- `a_locus_opening_on_the_regions_first_base_belongs_to_that_region` — widening
  `start < builder_region.start` to `<=` **passed all 104 tests** before, and loses roughly
  one locus in twenty from a run with nothing to say so.
- `a_locus_on_another_contig_is_not_this_regions`, rebuilt: the other contig's locus now sits
  at position 80 against a region ending at 50, so a break that compared positions alone
  truncates the whole region's output. Before, both loci sat at position 12 and the contig
  terms were unexercised.
- `the_keep_threshold_a_builder_is_given_is_the_one_the_walk_uses` — three non-reference reads
  built at the default of two and dropped at four. Before, `min_alt_obs` was never passed at a
  value of its own, so hardcoding the default inside `build_region` survived.

## Claims corrected

- **The failed spans' second job does not exist under this function's input contract.** Every
  builder is handed everything overlapping its ground, so two loci owned by different regions
  cannot overlap: what keeps the spans here is spec §3.3's count, not displacement. The doc
  says so, and hands the displacement question to the organiser (E2) rather than repeating
  the architecture's justification as though it were established here.
- **The cost of a long prefix is recorded with its measurement** — 3.3 µs per base at 63
  samples, 40 µs at 250, and what that means for a megabase — beside the note that the
  observation cache is what removes it.
- Three loose words: an 8-base record called a *chain*, a sample said to cover a locus's
  "tail" where it covers the eighth base of thirteen, and the organiser credited with summing
  a count `RegionOutcome` does not hold.

## Names

`RegionOutcome::observations` → `cohort_observations`, `failed` → `failed_locus_spans`,
`build_region`'s `region` parameter → `builder_region`. The type keeps arch §4's name; the two
field renames are a departure from its spelling, recorded here — `outcome.failed` reads as a
boolean, and `observations` collides with the meaning the module's own doc fixes for the word
(one sample's record over one stretch).

## Not applied

- **Renaming the type to `BuiltRegion`** — arch §4 fixes `RegionOutcome`, and "outcome" is not
  the word that misleads.
- **A test that the walk stops early rather than reading on** — the reviewer showed it is
  unpinnable through `RegionOutcome`, since `break` and `continue` return byte-identical
  outcomes here; the doc says what is observable instead.

## Validation

In the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — 111 passed, 0 failed.
</content>
