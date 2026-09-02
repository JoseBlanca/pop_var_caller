# C2 — the ±1 rescue, and the cohort's union

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone C step C2, and **Checkpoint C**. **Design:**
[`arch/candidate_alleles_ssr.md`](../../ng/arch/candidate_alleles_ssr.md) §3.1;
[`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md) §4, §4.1.
**Module:** [`src/ng/calling/allele_candidates/ssr.rs`](../../../../src/ng/calling/allele_candidates/ssr.rs),
`rescue_occupied_neighbours` and `promote_rungs_for_cohort`.

---

## What landed

**`rescue_occupied_neighbours`.** A diploid sample that resolved only one repeat length has a copy
unaccounted for, and what it most likely carries there is a neighbour of what it did resolve — a
second allele one repeat away, hidden under the first by stutter. So each resolved length's `±1`
neighbours are put forward, **but only where the cohort's reads actually reached that length**.
This is production's rescue and its `occupied` test
([`candidate_set.rs:239-258`](../../../../src/ssr/cohort/candidate_set.rs)), ported unchanged, and
it is the one part of production's nomination ng keeps.

**`promote_rungs_for_cohort`.** Every covering sample is asked its own question against its own
reads, and the cohort's promoted set is the **union**. A union and not a vote, for the same reason
one sample reaching the bar admits a sequence for the whole cohort on the ordinary path: an allele
one accession of sixty-three carries is still an allele, and a rule needing two would delete
exactly the rare variation a cohort is sequenced to find. The cost of a union is candidates, and
that is what the cap is for.

## The two silent failures this step is isolated for

The plan calls both out, and each is now a test with a failing state.

**Dropping the occupancy test invents a length.** Without it the rescue offers a repeat count
nothing in the run has ever seen, at every under-resolved sample of every tract. The test asks *the
ladder has a rung at that count **and** its cohort reads are non-zero* — production's
`cohort_support(length) > 0`. Those two conditions come apart at exactly one rung, which B1 already
pinned: the merge interns the reference tract whether or not a read landed on it, so the
reference's length always has a rung and may have no reads.

**Firing the rescue on a resolved sample widens every locus.** Up to two extra rungs a sample, and
every sequence they carry is a column in every genotype table for the life of the locus. The guard
is `promoted.len() < ploidy`, which is production's `peaks.len() < ploidy` — equivalent because a
promoted list shorter than `ploidy` is one the top-`ploidy` cut never truncated.

**The rescue is not itself capped at `ploidy`**, which is production's behaviour ported unchanged:
a sample that resolved one length of three occupied ones can come back with three. The cap that
does bind is the shared one over *sequences*, applied at admission in D1.

## Tests — 8 new

| test | what it pins |
|---|---|
| `a_sample_resolving_one_length_gains_its_occupied_neighbour_only` | three repeats is occupied by another sample and is offered; five holds only the reference and no reads, and is refused |
| `a_sample_resolving_two_lengths_gains_no_neighbour` | a live neighbour stays cut when the sample resolved its ploidy |
| `the_neighbour_is_one_repeat_away_and_not_one_rung_away` | rungs at 3, 4 and 6 — the rung *adjacent in the ladder* is two repeats away and is not reached |
| `a_zero_repeat_length_has_no_lower_neighbour` | a tract a deletion removed entirely, where the lower neighbour would underflow |
| `the_cohorts_promoted_set_is_the_union_across_samples` | three samples, one length each, all three promoted |
| `a_cohort_of_one_promotes_exactly_what_that_sample_does` | the union of one sample is that sample |
| `a_silent_sample_neither_contributes_nor_blocks` | a partial-only sample does not stop the samples after it |
| `a_second_locus_union_leaves_no_flag_of_the_first` | the flags are a per-worker buffer |

## What the mutations found

Three deliberate defects, applied, run, and copied back — all caught, and the first two are the
plan's own named silent failures:

| mutation | outcome |
|---|---|
| the rescue drops the occupancy test | caught — `a_sample_resolving_one_length_gains_its_occupied_neighbour_only` |
| the rescue fires whatever the sample resolved | caught — 2 tests |
| the union keeps only the last sample | caught — `the_cohorts_promoted_set_is_the_union_across_samples` |

Two fixtures were wrong on the first run and the failures were the tests' own, not the code's: one
gave a sample its support rows out of ascending allele order, which `one_run_per_allele` refuses by
design, and one expected a rescue at a sample that had in fact resolved its single copy. Both are
recorded here rather than quietly fixed, because a fixture that has to be corrected is a fixture
whose expected value was written before it was computed.

## Checkpoint C

Nomination is complete: the per-sample bar over rungs, the top-`ploidy` cut, the `±1` rescue and
the cohort's union. **The adjacent-length heterozygote test is green** —
`a_sample_with_equal_reads_at_adjacent_lengths_nominates_both`, from C1: a sample with 150 reads at
ten repeats and 150 at eleven puts both forward, where production's clear-peak rule resolves
neither because each length has the other beside it.

## Validation

All in the container (`./scripts/dev.sh`):

- `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **5,958 passed, 0 failed, 14 ignored**;
  `ng::calling::allele_candidates` at **136**, from 128 at C1 and 93 before this work.
- `cargo doc --no-deps` — 26 unresolved-link errors, unchanged from the pre-change tree, none in
  these files.
