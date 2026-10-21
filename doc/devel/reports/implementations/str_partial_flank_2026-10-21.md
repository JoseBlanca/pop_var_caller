# A partial read compared with an edge-deletion allele and the flank it leaves — implementation

**Date:** 2026-10-21. **Branch:** `str-partial-flank`, from `generic-indel-error` (b0162d47). Not merged.
**Contract:** [`str_edge_deletion_alleles.md`](../../ng/spec/str_edge_deletion_alleles.md) §5.6, "Reads that
end inside the segment: compared with the allele and the flank it leaves" (amended 2026-10-21).
**Measurements:** [`ng_generic_indel_error_2026-10-21.md`](../ng_generic_indel_error_2026-10-21.md) §5.2,
"Shipped". **Review:** [`reviews/str_partial_flank_2026-10-21.md`](../reviews/str_partial_flank_2026-10-21.md).

## What it does

At an STR segment, a read that stopped inside the repeat had held one flank and run towards the
other. When the candidate allele it is compared with deletes flank bases on that side (an edge
deletion), the read is now compared with the allele's repeat followed by the reference flank the
deletion leaves there, instead of with the repeat alone. Before, the bases such a read showed past the
allele's shorter repeat could only be explained as slippage, and at a homozygous edge deletion a dozen
such reads at 60× outweighed the complete reads and made the call `0/1`. Every other candidate — the
reference, plain repeat alleles, an edge deletion on the side the read held — is scored bit for bit as
before.

## Changes

All in [`src/calling/likelihood/ssr_emission.rs`](../../../../src/calling/likelihood/ssr_emission.rs):

- `remaining_flank_ahead` decides whether the rule reaches a candidate: a flank deletion on the side
  the read runs towards and something of the stored flank left past it.
- `censored_into_the_remaining_flank` (both sides; a right-held read is reversed, with its tract, the
  motif in the phase the tract ends on, and the flank) and
  `censored_from_the_tract_start_into_the_flank` (the sum over every reachable length change, with no
  factorised shortcut, admitting every stretching that with the flank covers the read).
- `letters_over` takes what follows the stretching (`followed_by`), so a slip is placed exactly as for
  every other partial read — in each run of an interrupted tract with equal weight — before the flank
  is appended. `letters_against_stretching_followed_by` replaces the prefix-then-compare pair at its
  three sites; with nothing following it adds the same scores in the same order as before.
- The mirrored buffers live in the model's existing scratch (`MirroredCensoredRead.flank`): nothing is
  allocated per read.
- Removed: the prototype `censored_with_the_flank_beyond`, its `PVC_EXP_PARTIAL_FLANK` switch and the
  deferred `all` scope.

## Assumptions

- **A deletion that leaves nothing of the stored flank** takes the tract-alone comparison. With an
  empty flank the rule's sum is the same quantity added in another order; keeping the shipped one
  keeps the last bits. The spec does not say; it is pinned by a test.
- **A candidate deleting on both sides** is read from each held flank by its far side's deletion alone,
  as the spec's per-side rule implies.

## Tests

Nine, each shown to fail when the rule it pins is reversed (mutations run on copies, recorded in the
review): the chr20:53879889 fixture against a hand-written sum and the reference; its right-held
mirror, to the bit; candidates the rule does not reach; four partial reads off the rule against bit
patterns recorded on b0162d47; slips placed per run on an interrupted tract; a stretching too short
even with the flank; a deletion taking the whole stored flank; deletions on both sides; and, in
`likelihood/ssr.rs`, a whole-tract read with no recorded side, half from each side's rule.

## Validation

- `cargo test --lib calling::likelihood`: 290 passed. The full library suite: 5,264 pass and one fails,
  `parameter_estimation::joint::fit::tests::a_frequency_that_is_not_a_probability_does_not_reach_the_caller`,
  which fails identically at b0162d47. `examples/ng_generic_loci_dump.rs` does not compile as a test
  target at b0162d47 either.
- `cargo clippy --all-targets --all-features -- -D warnings`: clean. `rustfmt --edition 2024` on the
  two touched files.
- The tomato identity oracle: b0162d47 reproduced all seven recorded checksums in the dev container;
  this build moves the three VCFs (no genotype, record, parameters file or window line), re-recorded
  with a note in `scripts/promote_ng_oracle.baseline`.
- Calls and scores on chr20, the tandem-repeat benchmark, the GIAB panel, the 5 Mb set and the
  stuttering benchmark: report §5.2.

## Trade-offs and follow-ups

- **One SNP lost** at a compound heterozygote on HG002's tandem-repeat benchmark at 30×
  (chr11:26874157); the prototype loses it too, and why the deletion wins as homozygous there is not
  traced.
- Comparing every partial read with every candidate followed by its flank stays deferred to
  `read_likelihoods.md` §5.2 (spec §9).
- `substitution_scores` is still looked up once per comparison on every partial read (inherited).
- `PROJECT_STATUS.md` has no block for this work; add one at merge if the owner wants it.
