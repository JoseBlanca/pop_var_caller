# STR bundles called as one STR locus — step 4 of `str_bundle_loci.md` (2026-10-11)

*Branch `str-bundles`, on 8fee2920, uncommitted. The design is
[`str_bundle_loci.md`](../../ng/spec/str_bundle_loci.md) §7; every choice it left open is written
there as a "Changed during implementation (2026-10-11, step 4)" note, and the measurements are in
its §10, step 4. This report says where the code is and what to look at first.*

## 1. What was built

An STR bundle — repeats within 15 bases of each other, typed as one region whose hull has clean
flanks — is now called. Before this, calling set every bundle aside.

| Piece (spec) | Where | Reuse of the segment's code |
|---|---|---|
| Run parse and run counts (§7.1) | [`region_typing/run_parse.rs`](../../../../src/region_typing/run_parse.rs) | new: an edit-distance alignment against the template `fixed (motif)* part fixed …` |
| Candidates (§7.2) | [`allele_candidates/ssr.rs`](../../../../src/calling/allele_candidates/ssr.rs), `select_ssr_bundle` | the segment's passes, keyed differently: `build_ladder_keyed`, `promote_rungs_for_cohort_by` with `RungNeighbours`, `admit_promoted_sequences_on` with `SelectionGround`, `locus_is_periodic_on_a_grid_of` |
| Read likelihood (§7.3) | [`likelihood/ssr_emission.rs`](../../../../src/calling/likelihood/ssr_emission.rs), `bundle_emission` | the segment's `segment_tract_into`, `place_slip_in_each_run`, `resize_at_the_end` and substitution factor, per run; the row (`likelihood/ssr.rs`) unchanged but for shifting a respelled allele's runs |
| Per-run stutter (§7.3) | [`inference/repeat_tract_parameters.rs`](../../../../src/calling/inference/repeat_tract_parameters.rs), `gather_for_bundle`, `bundle_candidates` | the segment's lookups, per run |
| Prior (§7.4) | [`genotype_prior/seed_ssr.rs`](../../../../src/calling/genotype_prior/seed_ssr.rs), `fill_ssr_bundle_seed` | the segment's weight reading, one factor per run |
| Evidence and loop | `LocusEvidence::Ssr { shape: SsrLocusShape }` ([`calling/mod.rs`](../../../../src/calling/mod.rs)); [`summarise_condition.rs`](../../../../src/calling/inference/summarise_condition.rs) | the segment's loop; no slippage re-fit at bundles |
| Driver | [`run/callers.rs`](../../../../src/run/callers.rs), `call_one_ssr_bundle_locus` | `call_one_ssr_locus`'s shape; complete observations only |
| Record (§7.5) | [`vcf/mod.rs`](../../../../src/vcf/mod.rs) `TractAnnotation::bundle`, [`vcf/encode.rs`](../../../../src/vcf/encode.rs), [`vcf/header.rs`](../../../../src/vcf/header.rs), [`vcf/assemble.rs`](../../../../src/vcf/assemble.rs), [`run/records.rs`](../../../../src/run/records.rs) | the segment's widening and gathering, over the bundle's flanks |

Two decisions of the session lead: the psp reader refuses a 1.x file by its version
(`OLDEST_READABLE_MAJOR = 2`, [`psp/header.rs`](../../../../src/psp/header.rs)); and a read ending
inside a bundle's hull is kept in the psp as a partial, as at a segment, so it counts toward the
hidden-duplication filter's coverage windows, and calling leaves it out
([`locus_generation/ssr.rs`](../../../../src/locus_generation/ssr.rs)).

## 2. The regressions that pin it

- **A one-run bundle spanning its hull scores every read as the segment does, bit for bit**
  (`a_one_run_bundle_scores_as_the_segment_it_is`: four hulls, over 1,000 read–allele pairs), is
  seeded as the segment is (`a_one_run_bundles_seed_is_its_runs_segment_seed`), and is called as the
  segment is through the whole loop — genotypes, GQs, site quality, expected copies, warrant
  (`a_one_run_bundle_is_called_as_the_segment_over_its_bases`).
- Two runs are checked against the slip-set sum written out by hand, the part-repeat regime by
  hand, and a bundle with no run.
- Every STR segment cell of the three benchmark runs is unchanged against main.

## 3. Validation

`cargo test --release --lib`: 5,060 passed, 3 failed — the three debug-only `should_panic` tests
known to fail in release. Integration tests: all pass. `cargo clippy --all-targets`: two warnings,
both in code this work did not touch. `cargo fmt`: applied.

The cross-platform VCF checksums moved; with the four changed header lines put back, each file
hashes to its recorded checksum, and they were re-recorded
([`cli/cross_platform_digests.rs`](../../../../src/cli/cross_platform_digests.rs)).

## 4. For the owner

1. **The periodicity gate costs bundles 31 of HG002's misses** and moves nothing on the simulator;
   switched off at bundles, HG002 goes from 824 / 158 / 149 to 842 / 169 / 120 (right / wrong /
   missed) with false calls unchanged. Built as the spec says; the switch is a decision.
2. **Deletions across a hull's edge are spelled inside the hull** by the flat-gap delimiter, and
   the generic path calls the same event beside it: 24 of the 51 wrong cells at 50 samples. Step 3's
   generator, not this step's caller.
3. **Few crossing reads** decide the rest of HG002's misses and the 5× gap — the reads ending inside
   are now in the psp, ready for the censored model §11's fifth decision deferred.
