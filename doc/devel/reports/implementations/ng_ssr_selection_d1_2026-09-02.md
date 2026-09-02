# D1 — which sequences on a promoted rung are admitted

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone D step D1. **Design:**
[`arch/candidate_alleles_ssr.md`](../../ng/arch/candidate_alleles_ssr.md) §3.2;
[`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md) §5, §6, §8.
**Module:** [`src/ng/calling/allele_candidates/ssr.rs`](../../../../src/ng/calling/allele_candidates/ssr.rs),
`admit_promoted_sequences`.

---

## What landed

A promoted rung says *this length is worth calling over*; it does not say which sequences at that
length are. **Every sequence on a promoted rung faces the shared support rule, asked of the
sequence** — the same bar the ordinary path asks, already folded into the per-allele summaries by
`summarise_alleles`. No representative is privileged and no recurrence term applies. Then the
shared cap and truncation, the shared leftover, and the reference tract admitted first and exempt
from both the bar and the promotion test.

**What that replaces.** Production promotes a rung's best-supported sequence unconditionally and
makes any sibling clear three further gates — 8 reads, **3 distinct samples**, and a tenth of the
rung's reads. The three-sample term has no cohort-size clamp, so **below three samples no second
spelling can ever be promoted**: the mechanism is absent rather than strict. Spec §5's measurement
of the class it matters for, at 300× on HG002 — a heterozygote whose two copies are the same length
spelled differently, 296 of that sample's 695 heterozygous tracts — is 35.1% for production against
86.1% for this rule, with a ceiling of 93.6% set by what some read actually carried. **Those are the
spec's offline numbers, not this code's**; reproducing them through the shipped module is
Milestone E.

**The leftover is filled although this path's likelihood never reads its pool.** Spec §8 is
explicit that a read no candidate explains is already carried by the junk term, spread over the
tract lengths the stutter model can reach — so `q_sum` here is computed and unread. What *is* read
is the other count: a sample whose own earned sequence the cap cut carries something the locus is
no longer called over, so every genotype the caller could form for it is wrong and it must be
emitted as missing. That rule is identical on both paths, so the leftover is built by the same
`leftover_of` rather than by a second walk.

## Tests — 8 new

| test | what it pins |
|---|---|
| `both_spellings_of_a_promoted_length_are_admitted` | one sample, two five-repeat spellings, both in — the case production cannot reach below three samples |
| `a_sequence_on_an_unpromoted_rung_is_not_admitted` | nine reads at six repeats, on a length two copies had already accounted for |
| `a_sequence_below_the_bar_on_a_promoted_rung_is_not_admitted` | a rung is a length, not a licence for every spelling at it |
| `the_reference_tract_is_admitted_although_no_read_reached_it` | exempt from the bar *and* from promotion |
| `above_the_cap_the_worst_ranked_sequences_are_cut_and_the_locus_is_still_called` | truncation, not refusal, with the count reported |
| `a_sample_whose_earned_sequence_the_cap_cut_is_marked_uncallable` | the leftover count the tract likelihood does read |
| `the_candidate_table_carries_the_tracts_kind` | `Ssr` with the motif, not `Generic` |
| `the_admitted_order_is_the_merge_tables_and_not_the_ladders` | the ALT column's order is the merge's |

## What the mutations found, and one fixture that was wrong

Three deliberate defects, all caught: admission ignoring the promotion flag; admission ignoring the
support bar on a promoted rung; and the candidate table minted `Generic` instead of the tract's own
kind.

**One fixture failed on the first run and the failure was the test's, not the code's** — worth
recording because the mechanism is one a reader would get wrong the same way. The cap fixture had
three samples at 60, 40 and 10 reads and expected the ten-read sample to lose its sequence. It did
not: the sample at 60 reads was *heterozygous* in that fixture, four of its 64 reads on the
reference, so its within-sample share was 0.94 against the two shallower samples' 1.00 — and the
cap's **first** ranking key is the largest share of one sample's reads, not the read total. The
deep sample lost its allele to a sample carrying ten reads. That is the shared ranking behaving
exactly as its own documentation describes ("a sample sequenced at 3 reads outranks every
heterozygous sample sequenced at 300"), so the fixture was made homozygous throughout: all three
shares are then 1.0, the cohort read total decides, and which sample is cut follows from the
numbers on the page.

## Validation

All in the container (`./scripts/dev.sh`):

- `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **5,966 passed, 0 failed, 14 ignored**;
  `ng::calling::allele_candidates` at **144**, from 136 at C2.
