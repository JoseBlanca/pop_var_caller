# Fixes applied — ng cohort merge, C2

*2026-08-17, against [the review](ng_cohort_merge_c2_2026-08-17.md). In
`src/ng/run/cohort_merge/serial.rs`, with one doc correction in `build.rs`.*

## Behaviour changed

- **Analysed regions that overlap, repeat or descend are refused**, at release level, with the
  message naming both regions. They used to produce duplicated observations silently — the
  locus in the shared ground built once per region — which no consumer can tell from a cohort
  that varied at two places.

## Tests, and the mutation each kills

115 in the module, from 111 at review time.

- **The end-to-end fixture now discriminates the quality rule.** One of the substitution
  sample's three reads carries a bad base (quality 5) at position 110, so its six records
  inside the locus no longer share one quality sum. Flipping the composition from the weakest
  sighting to the best now **fails** the test; before, every base was quality 30 and the two
  rules were indistinguishable there.
- `the_same_loci_come_out_however_the_analysed_ground_is_divided` — sixty loci over six
  hundred bases, built as 1, 6, 60 and 120 regions. This is spec §15's regression anchor, and
  the file claiming to be the oracle had no test of it.
- `analysed_regions_that_overlap_are_refused` and `analysed_regions_out_of_order_are_refused`.
- `the_keep_threshold_reaches_the_builders` — three reads built at the default of two and
  dropped at four; every other test here ran at the default, so dropping the argument shipped
  green.
- The failed-spans test now puts **one failed locus in each of the two regions**. With both in
  the first, stopping after one region or walking them backwards passed it. Stopping after the
  first region now fails four tests.

## Claims corrected

- **Both fixture positions in the end-to-end test's own doc were five bases low** — 107 and
  105–109 where the reads mint 112 and 110–114 — consistent with each other and contradicted
  by the assertion below them.
- **The driver's cost with many regions is stated with its measurement**: the same 20,000
  observations cost 5.4 ms in one region and 184 ms in a thousand, and `build_region`'s claim
  that the serial driver's prefix is "empty by construction" was true only at one region.
- **`alleles_of_sample`'s gap explanation is corrected in `build.rs`.** The generic mint writes
  a record at every covered position — thirty for a thirty-base read, as the new fixture shows
  — so a gap between a sample's records is ground no read covered, not ground where its reads
  agreed with the reference. The filling stays, because this walk does not depend on which
  generator minted its input.

## Not applied

- **Slicing the observations per region inside the driver.** That is the observation cache's
  job (milestone D); doing it here would build a second, untested windowing and make the oracle
  the thing most likely to be wrong.

## Validation

In the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — 115 passed, 0 failed.
- `cargo test --lib` — 3,738 passed, 0 failed, 11 ignored.
</content>
