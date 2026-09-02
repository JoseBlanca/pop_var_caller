# E1 — the existing caller's rules switched back in, on 269 tomato tracts

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone E step E1 (executed from [`calling_loop_ssr.md`](../../ng/impl_plan/calling_loop_ssr.md)
Milestone B). **Design:** [`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md)
§10; [`arch/candidate_alleles_ssr.md`](../../ng/arch/candidate_alleles_ssr.md) *Test & bench shape*.
**Module:**
[`src/ng/calling/allele_candidates/ssr_production_differential.rs`](../../../../src/ng/calling/allele_candidates/ssr_production_differential.rs)
(test-only), with its fixture under
[`testdata/`](../../../../src/ng/calling/allele_candidates/testdata/).

---

## What landed

**One end of the differential spec §10 asks for.** Three of the existing caller's rules are
replaced on purpose — its clear-peak nomination, its cohort-summed depth gate, and its
same-length sibling bar — so a byte-identical parity test is impossible by construction. What
exists now instead is a test-only re-implementation of those three, driving **ng's own fold and
ng's own ladder**, required to reproduce the existing caller's candidate set on real reads.

**The comparison runs against the existing caller's own code, not against frozen expectations.**
`assemble_candidates` and `build_rungs` are `pub(crate)`, so the test builds each fixture tract
twice — once as a `CohortObservation` and once as a `CohortLocus` — and compares the two answers
live. Nothing in the fixture can go stale, and a change on either side turns it red.

**The rescue is called, not copied.** Spec §10's reuse map records the `±1` occupied-neighbour
rescue as ported unchanged, so the arm calls the shipped `rescue_occupied_neighbours` rather than
writing production's version a second time. The other three rules are written out, with their
constants retyped rather than read from `CandidateCfg`, so that a change to the existing caller's
development defaults turns this red instead of quietly moving both sides together.

## The fixture

`testdata/tomato_tract_alleles.csv` and `testdata/tomato_tract_rows.csv` — the merge's allele
table and every covering accession's reads at **269 repeat tracts of the 51-accession tomato
panel**, from the first 400 intervals of `benchmarks/ssr_tomato1/ssr_regions.bed`, written by
`examples/ng_candidate_selection_probe`'s new `NG_TRACT_DUMP` and `NG_TRACT_ROWS` outputs. 366 KB
for 1,376 alleles and 12,223 evidence rows; a sequence is named by its index into its tract's
allele table rather than by its bases, which is what keeps the rows file checkable-in.

**Two things it drops, and neither is free in general.** Read groups are pooled to one count an
accession, which is the shape the existing caller's own Stage-1 evidence has; and the merge's
partial observations are not carried, so an accession's compared reads here are the sum of its
rows. None of the three rules under test reads either — nomination counts reads at a length, the
depth gate sums the cohort, the sibling bar divides by a rung's own total.

## What the two rules actually do differently, measured

On the same 269 tracts at ploidy 2, **the shipped rules narrow 184 of them differently** from the
existing caller's, to **602 candidate sequences against 668**. The three numbers are pinned in the
test rather than bounded, because they are the measurement the fixture exists to make.

The direction matters: spec §4.1 and §5 measure the replacement as *cheaper* in candidates as well
as better in recall, on human data at 300×. On a 51-accession panel at about three reads a
position it stays cheaper — 602 against 668 — which is the other corner of the range the caller
has to work across.

## Three of the five rules are invisible on a 51-accession panel, so three tests exist for them

**Predicted before running, and the predictions were what found the gaps.** Six mutations were
run against the arm. Three failed immediately: the prominence constant 3 → 0, the sibling
read floor 8 → 2, and dropping the rung representative's unconditional promotion. **Three
survived**, and each survivor is a fact about the rule rather than a thin fixture:

- **the cohort depth gate** (10 reads summed across the panel) — fifty-one accessions at three
  reads clear it at every tract, so raising or lowering it changes nothing;
- **the three-distinct-accessions recurrence term** — at three reads an accession, eight cohort
  reads on a sequence already implies three accessions showed it, so the term never decides;
- **the periodicity gate** — no tract of the fixture is refused by the existing caller's own
  mode-anchored measure.

Three tests were added and all six mutations now fail:

| test | what it reaches |
|---|---|
| `productions_rules_are_reproduced_with_one_accession_too_where_two_of_them_first_bite` | the same 269 tracts with one accession: the depth gate then refuses **every one** of them down to the reference tract, where the whole panel refuses none |
| `productions_recurrence_term_is_reproduced_on_a_sibling_two_accessions_carry` | hand-built: a sibling with 20 reads from two accessions, clear of the eight-read floor and of a tenth of its rung's 120, refused for the third accession alone |
| `productions_periodicity_gate_is_reproduced_on_an_off_grid_tract` | hand-built: 12 of 22 reads one base off a dinucleotide grid |

**The one-accession test is the range check `CLAUDE.md` asks for**, and it lands on the exact
sentence spec §6 makes: the same tract is refused alone and admitted in company.

## Tests — 5 new

Two run over the real fixture (the reproduction, and the "the two rules must actually differ"
counter-test), three over the cases the panel cannot reach. All five in
`ssr_production_differential.rs`.

## Deviations from the plan

- **The arm lives in its own `#[cfg(test)]` file** rather than inside `ssr.rs`'s `mod tests`. It
  is 636 lines with a fixture reader in it, and `ssr.rs` is already the largest file in the
  module. It is a sibling module under `allele_candidates`, so it still reaches the `pub(super)`
  ladder and fold the arm has to drive.
- **The candidate-selection probe's whole tract-arm diff lands here**, although only two of its
  new outputs — the per-tract candidate dump and the per-accession evidence rows — are what this
  step's fixture is cut from. The rest (a support-share override, and the count of catalog tracts
  the merge refused) is Milestone E's measurement machinery and is used by E2; separating them
  inside one function would have been a contrivance.
- **The comparison is on the candidate *set*, not on its order.** The existing caller walks its
  nominated lengths in support order and ng's rescue returns its rungs ascending, so the two push
  the same sequences in different orders. The shipped path's own order is the merge table's and is
  pinned by the module's unit tests.

## Validation

`cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib` 5,999 passed / 0 failed / 14 ignored, all in the container.
