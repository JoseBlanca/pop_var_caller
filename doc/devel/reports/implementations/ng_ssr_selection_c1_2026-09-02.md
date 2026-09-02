# C1 — which repeat counts one sample puts forward

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone C step C1. **Design:**
[`arch/candidate_alleles_ssr.md`](../../ng/arch/candidate_alleles_ssr.md) §3.1;
[`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md) §4, §4.1.
**Module:** [`src/ng/calling/allele_candidates/ssr.rs`](../../../../src/ng/calling/allele_candidates/ssr.rs),
`promote_rungs_for_sample`, with one buffer added to the shared scratch.

---

## What landed

**`promote_rungs_for_sample`** takes one sample's length histogram and its spanning-read total and
returns the rung indices that sample puts forward, ascending, at most `ploidy` of them.

**The rule is the shared one, asked of a rung.** A repeat count is nominated when this sample's
reads at it reach `max(2 reads, ceil(share × this sample's spanning reads))` —
[`MinAltReads::reached_by`], the same predicate the merge asks of a sample's non-reference reads
and the ordinary path asks of one sequence. Neither the numerator nor the denominator is computed
here: the numerator is B2's histogram and the denominator is `compared_reads_of`, so there is one
producer of each.

**Then the best `ploidy` of them, ties to the shorter repeat count.** The buffer holds two orders in
turn, exactly as the ordinary path's cap does — ranked while the cut chooses, then sorted back into
ascending rung order, because that is the order everything downstream reads.

**Production's clear-peak test is not called and not ported.** It nominates a length only if its
reads exceed *both* neighbours by more than three, so a heterozygote whose two copies differ by one
repeat resolves nothing — each length has the other beside it. Nothing in this function reads a
neighbour.

## Two properties the range demands, both asserted

**A cohort of one and a cohort of a thousand give a sample the same answer.** The denominator is
that sample's own spanning reads and no term of the bar reads the cohort, so nothing about who else
is in the run can change what this sample nominates.
`a_samples_nomination_is_the_same_alone_and_in_a_cohort` runs the same sample alone and beside a
neighbour carrying 400 reads at a third length, and gets the same two rungs.

**A sample with no spanning reads nominates nothing, without dividing by its zero.** The floor is
at least one read and `reached_by` tests the floor first by integer comparison, so a sample whose
reads all stopped inside the tract falls out before the share is computed —
`a_sample_with_only_partials_nominates_nothing`.

## Tests — 9 new

| test | what it pins |
|---|---|
| `a_sample_with_equal_reads_at_adjacent_lengths_nominates_both` | **spec §13's test that production cannot pass**: 150 reads at ten repeats and 150 at eleven promote both |
| `a_diploid_sample_promotes_its_two_best_supported_rungs` | the cut keeps 30 and 40 reads over 4 and 6, and returns them in **ascending rung order** where the ranking had them reversed |
| `a_triploid_sample_promotes_three` | the ploidy is the run's, and is the only thing that changes |
| `a_tie_in_support_goes_to_the_shorter_repeat_count` | twenty reads each, one copy, the shorter wins |
| `a_rung_below_the_bar_is_not_nominated_even_with_copies_to_spare` | the bar and the cut are two questions |
| `the_bar_is_a_share_of_the_samples_own_spanning_reads` | two reads in ten clear a fifth; the same two in a hundred do not |
| `a_samples_nomination_is_the_same_alone_and_in_a_cohort` | the cohort-size property above |
| `a_sample_with_only_partials_nominates_nothing` | the zero-denominator case |
| `nominating_for_a_second_sample_leaves_none_of_the_firsts_rungs` | the buffer is one per worker |
| `a_histogram_that_is_not_this_ladders_is_refused` | one entry per rung is what makes a rung index mean one thing |

## What the mutations found

Three deliberate defects, applied, run, and copied back — all caught:

| mutation | outcome |
|---|---|
| the tie goes to the longer repeat count | caught — `a_tie_in_support_goes_to_the_shorter_repeat_count` |
| the promoted rungs are left in ranked order | caught — 2 tests |
| the promoted buffer is not emptied between samples | caught — `nominating_for_a_second_sample_leaves_none_of_the_firsts_rungs` |

**The second mutation is why the diploid fixture's read counts are what they are.** It was first
written with the better-supported rung at the lower index, so the ranked order and the ascending
order were the same list and dropping the sort back changed nothing. The numbers were swapped —
30 reads at four repeats, 40 at five — so the cut ranks them 2 then 1 and the sort back has a
failing state. The tenth test, on the per-sample buffer, was added for the same reason: every other
fixture used a fresh buffer, so an append in place of a refill was invisible.

## Validation

All in the container (`./scripts/dev.sh`):

- `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **5,950 passed, 0 failed, 14 ignored**;
  `ng::calling::allele_candidates` at **128**, from 118 at B2.
