# Fixes applied — ng cohort merge D2 (the merge read through the cache)

*2026-08-17. Against [the D2 review](ng_cohort_merge_d2_2026-08-17.md); step D2 of
[the plan](../../ng/impl_plan/cohort_merge.md).*

**Everything actionable was applied.** Seven Majors, twelve Minors and the nits; nothing
deferred to the owner, nothing disputed. The module's tests go from 154 to 164.

## The Majors

- **M1, the inverted analysed region.** The guard both drivers share now refuses a region whose
  ends are the wrong way round, and is named for the two rules it enforces —
  `refuse_malformed_analysed_regions`. **This changes the oracle too**, which is the point: the
  two drivers must be given the same input or the byte-identity claim means nothing. No legal
  input changes behaviour; an illegal one now panics where it used to be read two ways.
  `an_analysed_region_with_inverted_ends_is_refused` pins it, and the division's own reading of
  such a region is pinned in `organise.rs` as the belt-and-braces it now is.
- **M2, the division pinned at the divider rather than at the driver.**
  `the_width_the_caller_asks_for_is_the_width_the_driver_builds_in` merges one fixture at two
  widths and reads the cache: **60 observations held at one region for the stretch, 2 at
  twenty-base regions.** It fails for a driver that ignores the width whatever eviction does.
- **M3, covering the whole analysed region.** `the_window_stays_short_up_to_a_failure` reads the
  window mid-stretch, at a source failure — the only place from outside where the drawing pace
  is visible.
- **M4, the overlap guard's boundary.** `analysed_regions_sharing_one_base_are_refused`.
- **M5, a failed merge leaves the cache advanced.** Documented on the driver, in the terms the
  review measured: the same cache cannot retry the same ground, because a second merge over it
  comes back short and says `Ok`. A run that means to retry builds a new cache. Making it
  unrepresentable is the organiser's, at milestone E.
- **M6, the held count over one sample.** Two tests in `organise.rs` —
  `the_held_count_sums_every_samples_window` (zero before anything is drawn, two across two
  samples afterwards) and `the_held_count_is_zero_over_no_samples`.
- **M7, the ordering panic against a `Result` signature.** Named on the driver: not every
  source-side failure arrives as `Err`, why that is right today, and what changes when
  observations come from a psp file.

## The Minors and nits

- **Mi1** — the eviction test now says what actually produces its number: the records at 581 and
  591, both inside the last building region, with the source spent after 591; and that at
  five-base regions the same fixture ends holding none, so the number moves with the width rather
  than being a bound. The report carries the same correction.
- **Mi2 + Mi3 together** — `both_drivers_agree` becomes
  `the_outcome_both_drivers_agree_on`: it asserts inside, renders locus by locus so a
  disagreement names the first entry that differs, and returns the outcome the two agree on. Five
  call sites lose their own comparison.
- **Mi4** — `building_regions_of` moves to `organise.rs`, beside the cache, because handing
  regions out is the organiser's job at milestone E; its two tests move with it.
- **Mi5** — the file header now describes both drivers, says which is the reference, and carries
  the vocabulary.
- **Mi6** — *analysed region* and *building region* are defined at the head of the file that
  turns on them.
- **Mi7** — the driver's doc gives the cost beside the saving: one cover per building region,
  616 µs at 1,000 samples and 2.87 ms at 3,000 on 20-base regions.
- **Mi8** — `held_observations` → `held_observations_len`.
- **Mi9** — the coordinate-ceiling test takes three regions instead of collecting an unbounded
  iterator, and says why.
- **Mi10** — a `#[cfg(test)] mod fixtures` in `mod.rs` now holds `region`, `region_on`,
  `position_on` and `SourceFailed`; `organise.rs` and `serial.rs` read them from there. `build.rs`
  and `close.rs` still carry their own copies and are named as the next two to fold in.
- **Mi11** — the eviction test states that the evict-then-cover *order* is not what it pins, and
  what the order costs.
- **Mi12** — the report's counts and table are corrected, and the four wrong claims are named in
  a section of their own.
- **Nits** — the width has one name per level (`cohort_locus_builder_regions_len` at the
  driver, `building_region_width` at the divider, `bases_per_region` for the raw count); the
  `width - 1` carries the reason it cannot underflow; the redundant second merge in one test is
  gone.

## Not done, and why

- **A method on `RegionOutcome` for gathering one region's outcome**, which would make the two
  drivers' accumulation structural rather than duplicated. `RegionOutcome` lives in `build.rs`,
  committed at C1 and outside this step.
- **Bounding `E` so the cache can mint its own ordering error.** It changes every `impl` block
  and every `E` in `organise.rs`, and the `RunError` it should become does not exist yet. Named
  in the doc instead.
- **Folding `build.rs` and `close.rs` onto the shared fixtures.** Recorded as an open item.

## Mutation re-run

Thirteen mutations, rewritten for the new code: **all thirteen fail at least one test.** Five of
them are the survivors the three reviewers found — the overlap guard's shared base, the inverted
region, the width ignored at the call site, the cover over the whole analysed region, and the
held count over one sample.

## Validation

`cargo fmt --check` clean; `cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::run::cohort_merge` → `164 passed; 0 failed`; the whole library suite green.
