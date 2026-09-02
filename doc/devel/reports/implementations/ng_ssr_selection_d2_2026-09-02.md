# D2 — the periodicity verdict, and where its grid is anchored

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone D step D2. **Design:**
[`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md) §7 and
[`arch/candidate_alleles_ssr.md`](../../ng/arch/candidate_alleles_ssr.md) §3.3, **both amended by
this step** with a dated note recording the owner's decision below.
**Module:** [`src/ng/calling/allele_candidates/ssr.rs`](../../../../src/ng/calling/allele_candidates/ssr.rs),
`locus_is_periodic` and the entry point `select_ssr`.

---

## What the verdict decides

A stretch the catalog called a repeat tract, whose reads sit at lengths the motif cannot explain,
is not something this caller's model describes: the stutter distribution is written on
whole-repeat and part-repeat regimes and the prior's ladder on repeat counts. So it is refused —
the reference tract alone comes back, under `SelectionVerdict::NotPeriodic`, and no other length is
ever called there. A refusal that still yields a usable table, so what the run does with such a
locus stays emission's decision.

A read is **off the grid** when the difference between its tract length and the anchor is not a
whole number of motif units. A sample is non-periodic when more than one read in ten of its
spanning reads is off the grid, and the locus is refused only when **no** sample is periodic — the
same "one sample suffices" shape as the support bar.

## The anchor — an owner's decision, and it is not what the documents said

**The design documents anchor the grid at zero, and that refuses a real class of tract.** Arch §3.3
and spec §7 say a read is off-grid when its length is *"not a whole number of motif units from the
ladder's mode"*, and spec §3 adds that ng measures this in units where production measures it in
bases. Taken literally the mode cancels: it is a repeat count, so its length in bases is a whole
number of units, and the test collapses to *the read's length is not a whole multiple of the
period*.

**Why that loses tracts, measured rather than argued.** The catalog trims every tract back to whole
motif copies **at both ends**, but a length-changing interruption inside puts the two ends out of
phase with each other, so the tract's own reference length is then not a multiple of the period.
Such a tract is admitted whenever the interruption is late enough to clear the catalog's purity
floor of 0.8. Run through the catalog's own `minimal_trim` and `recompute_purity` on 2026-09-02:

| tract | trimmed length | length mod period | purity |
|---|---|---|---|
| `AT` × 6, clean | 12 | 0 | 1.000 |
| `AT` × 4, one extra base, `AT` × 4 | 17 | 1 | 0.471 |
| `AT` × 8, one extra base, `AT` × 2 | 21 | 1 | 0.762 |
| **`AT` × 20, one extra base, `AT` × 4** | **49** | **1** | **0.816** |

The last row clears the 0.8 floor, so the catalog admits it. Under a zero-anchored grid every read
at that tract's reference length is off the grid, every sample is non-periodic, and **the locus is
refused and never called**.

**Three anchors were put to the owner and the third was taken (2026-09-02).** Zero, as the
documents say; the commonest observed length in bases, as production does
([`candidate_set.rs:114-145`](../../../../src/ssr/cohort/candidate_set.rs), whose own comment says
the anchor is the modal length "not zero, so an interrupted repeat sitting at an odd reference
length stays periodic"); or **the reference tract's own length**. The third keeps those tracts *and*
is a property of the locus rather than of the reads — the commonest observed length moves with
depth, so a shallow sample could shift which lengths count as on-grid — and it is the same quantity
the genotype prior was re-indexed onto on 2026-08-27, offset from the reference tract length. Spec
§7 and arch §3.3 now carry a dated note saying so.

## Two counting rules the design left open, decided here

**A sample with no spanning reads is not asked, rather than counted as periodic.** Counting it as
periodic would make one silent sample enough to save every locus in a run — and a tract too long
for a read to span produces exactly such a sample, so the verdict would become unreachable wherever
coverage is thin, which is where it is most needed.

**A locus no sample spanned is periodic by default.** There is nothing to judge it on, and refusing
it would be a verdict about coverage rather than about the tract.

**A homopolymer can never be off the grid**, since every difference is a whole number of one-base
units. Production short-circuits period 1 explicitly; here it falls out of the arithmetic, and
`a_homopolymer_is_always_periodic` is what says the arithmetic really does it.

## `select_ssr` — the entry point

Four passes in a fixed order: the periodicity verdict first of all, then the fold, the ladder,
nomination and admission. **The tract's detail now comes off `CohortObservation::kind`**, which the
observations plan landed on `main` earlier the same day — arch §1 wrote the separate `detail`
argument as the shape to use "until it lands", and it has. A locus of the wrong kind panics rather
than being scored against the other path's model.

## Tests — 9 new

| test | what it pins |
|---|---|
| `a_tract_at_an_odd_reference_length_is_periodic_about_its_own_reference` | the 49-base tract above, periodic here and refused under a zero anchor |
| `a_tract_whose_reads_are_off_the_motif_grid_is_refused_but_still_yields_the_reference` | the refusal is a table, not an error |
| `one_periodic_sample_saves_a_locus_every_other_sample_fails` | the cohort verdict takes one sample |
| `the_off_grid_share_decides_at_one_read_in_ten` | 4 of 40 passes, 5 of 40 does not |
| `a_homopolymer_is_always_periodic` | period 1 |
| `a_sample_with_no_spanning_reads_does_not_vote` | a silent sample does not save the locus |
| `a_locus_no_sample_spanned_is_periodic_by_default` | nothing to judge |
| `the_entry_point_refuses_a_locus_that_is_not_a_tract` | a routing bug is loud |
| `the_entry_point_offers_both_spellings_of_one_length` | end to end through `select_ssr` |

## What the mutations found

Four deliberate defects, all caught, each by exactly one test:

| mutation | outcome |
|---|---|
| **the grid anchored at zero — what the documents literally say** | caught — `a_tract_at_an_odd_reference_length_is_periodic_about_its_own_reference` |
| the share boundary excludes exactly one read in ten | caught — `the_off_grid_share_decides_at_one_read_in_ten` |
| a sample with no spanning reads votes periodic | caught — `a_sample_with_no_spanning_reads_does_not_vote` |
| a locus nobody spanned is refused | caught — `a_locus_no_sample_spanned_is_periodic_by_default` |

The first is the one worth naming: **reverting to the design documents' own rule now fails a test**,
which is what makes the amendment to spec §7 and arch §3.3 enforceable rather than a note somebody
has to remember.

## Validation

All in the container (`./scripts/dev.sh`):

- `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **5,975 passed, 0 failed, 14 ignored**;
  `ng::calling::allele_candidates` at **153**, from 144 at D1.
- `cargo doc --no-deps` — 26 unresolved-link errors, unchanged from the pre-change tree, none in
  these files.
