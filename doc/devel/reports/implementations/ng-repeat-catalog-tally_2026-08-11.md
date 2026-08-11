# ng repeat catalog — the tally: what a consumer keeps when it stops walking the reference

*Implementation report, 2026-08-11. Wiring work towards
[`spec/typed_regions.md`](../../ng/spec/typed_regions.md)'s step 3 reading the catalog instead of
scanning ([`spec/repeat_catalog.md`](../../ng/spec/repeat_catalog.md) §8, first deferral). Design:
spec §5.1 and §3.1.*

## Why this had to land before any consumer switches

A walk over the reference counts as it goes — how many STR loci, how many bundles, how many bases of
repeat that yielded no locus, and what admission turned down and why — and its consumers print that
at the end of a run (`type-regions`, `ng_str_stutter_by_library`, `ng_ssr_loci_dump`). The catalog
had no equivalent, so a consumer that swapped a scan for the file lost its report. This is the
counterpart, and the question it turned on was **which of those counters the file can honestly
answer**.

## The answer, and the one omission

**Every counter reproduces exactly except one: repeats turned down for touching a contig's very
first or last base.** The file holds no repeat with less than 15 bases of sequence beside it (spec
§4.1), so a repeat sitting on base 1 was never written down and no reader can see it. A counter for
it would report `0` for every genome — zero by construction, not zero in the genome.

**So there is no such field.** [`CatalogRegionCounts`](../../../../src/ng/repeat_catalog/tally.rs)
carries the walk's eight region counters and four of its five rejection reasons; the fifth has no
slot. A field that cannot exist cannot be misread as a measurement, and a field holding `0` will be.

The cost is that the type is not the walk's, so a switching consumer edits its printing line. Three
consumers read these counts today.

Two differences follow at a contig's ends, and both are stated in the tests rather than assumed:

- a repeat 1 to 14 bases from a contig end is a **locus** to a live scan and absent from the file
  altogether, so the loci and generic counts differ there too. Measured over whole genomes: 7 such
  tracts in tomato SL4.00 and 857 in GRCh38;
- everything else agrees exactly, because both sides run the same pre-screen and the same admission
  gates — the file merely supplies the whole-motif cut, the motif and the purity instead of the
  bases.

## A defect this found in the walk, and it is fixed here

**The walk was charging a base covered by two overlapping repeats twice.** Its accumulator summed
the clipped coverage spans; the pre-screen removes period-multiple re-detections of *one* tract, not
two different repeats that intersect, so the sum charges a shared base once per tract over it.

It is a defect rather than a difference of opinion: the walk's own test for this number computes it
as `coverage_runs(&cleaned)` — merged — and the satellite cap in the same function merges before it
measures. The accumulator was the one place that did not, and its fixture had no overlap, so nothing
disagreed until the catalog's tally was compared against it and came out a base short over 200 kb of
sequence.

Fixed by unioning before summing (`covered_bp`), with a unit test that a shared base counts once.
Window-invariance is untouched: cores tile and do not overlap, so the union of the per-core unions is
the contig's union whatever `window_bp` is.

**How much it was off by: 5,704 bases on tomato chromosome 1, 371,085 against a true 365,381** — one
and a half per cent of the number, and every walk that has ever printed it was that much high.

## What the numbers look like on a real chromosome

Tomato chromosome 1 (90 Mb), walked at the calling floors, which is the measurement that decided how
much any of this matters:

| emitted | | turned down, in bases | |
|---|---|---|---|
| STR loci | 57,540 | no whole-motif cut | 207,854 |
| bundles | 2,910 (105,891 bp) | too impure | 33,253 |
| satellites | 26 (38,230 bp) | too few copies | 15,743 |
| generic stretches | 60,477 | motif is itself a repeat | 1,862 |
| repeat bases with no locus | 365,381 | **touches a contig end** | **0** |

**All four of the reasons the catalog can answer fire on real sequence.** A doc comment in
`segment_criteria.rs` still says four of the five are structurally zero and only the contig-end one
is reachable; that was true before the copy floor moved to the cut tract (2026-07-20) and it is now
the reverse. Left for the review pass, since it is prose rather than behaviour.

## Assumptions and deviations

1. **Only `genome_segments` carries a tally.** `str_loci` stops before the generic stretches are
   built, so a tally there would report a generic count of zero — the same silent-zero failure the
   whole design is avoiding. A consumer of `str_loci` counts loci itself.
2. **`finish_from_row` returns the walk's own `RejectionReason`**, not a second vocabulary. The two
   tallies are compared against each other, so one name per gate is what makes the comparison mean
   anything.
3. **The tally is the whole contig's even for a region read.** That is not a choice: the walk widens
   every requested region to its whole contig before scanning (`scan_set`), so both sides type the
   same sequence and the numbers agree without a special case.
4. **Two of admission's gates cannot be driven by a fixture reference.** An impure tract and a
   compound motif do not survive to admission from a scan — the detector emits the purest
   sub-segment and only primitive motifs — and 2 Mb of random sequence produced neither. They are
   pinned by handing admission the rows directly, the way `segment_criteria`'s own compound-motif
   test does.

## Changes made

- [`src/ng/repeat_catalog/tally.rs`](../../../../src/ng/repeat_catalog/tally.rs) — new:
  `CatalogRegionCounts`, `CatalogRejectionCounts`, `ContigTally`, and the module doc that says which
  question the file cannot answer and why there is no field for it.
- [`src/ng/repeat_catalog/segments.rs`](../../../../src/ng/repeat_catalog/segments.rs) — admission
  records *why* it turned a repeat down instead of dropping the reason; the segment functions hand
  back a `ContigTally` beside their regions.
- [`src/ng/repeat_catalog/reader.rs`](../../../../src/ng/repeat_catalog/reader.rs) —
  `GenomeSegments::counts()`, and the per-region tally the walk keeps.
- [`src/ng/region_typing/mod.rs`](../../../../src/ng/region_typing/mod.rs) — `covered_bp`, and the
  accumulator uses it (the defect above).

## Tests added

Eight, and the fixture grew one contig.

- **The two tallies agree, counter for counter**, over six contigs of the differential's fixture plus
  200 kb of sequence with no designed structure. That last contig is there because the hand-built
  ones drive no rejection at all, which would have made the comparison assert four zeroes against
  four zeroes; over 200 kb the detector produces tracts that lose a copy to the cut and tracts with
  no whole-motif boundaries at all, and the test asserts both counters are non-zero before comparing
  them.
- **The differences are named and measured** on the two contigs with structure at their ends: 60
  bases charged to the contig-end rejection that the file cannot see, one locus the walk finds and
  the file does not, and the generic stretch that locus splits in two.
- **Each gate charges its own counter, and charges the detected length** — four rows handed to
  admission directly, including the impure and compound cases a scan never reaches.
- **A row that clears every gate charges nothing**, so the counters are not simply charging
  everything they see.
- **Two tracts sharing a base cover it once**, at both levels: `covered_bp` directly, and through
  the derived segmentation.
- **A repeat out of the period range is not a rejection** — it is out of the question being asked,
  which is how `classify` treats it too.
- **The whole reference is one span per non-empty contig**, the rule that makes the `spans` counter
  the same number on both sides.

### What mutation testing found

Five mutations, four caught immediately:

- summing the coverage instead of unioning it → the tally comparison fails (this is how the walk's
  defect was found in the first place);
- charging the cut length instead of the detected one → the gate test fails;
- dropping the compound charge → the gate test fails;
- forgetting to cancel a locus's bases from the no-locus gap → the tally comparison fails.

The fifth survived and is why there is an eighth test: **deleting the whole-reference `spans` count
left everything green**, because every tally test asked for regions rather than for the whole
reference. `the_whole_reference_is_one_span_per_non_empty_contig` now fails on that mutation.

## Validation

In the dev container: `cargo fmt` clean; `cargo clippy --lib --tests --all-features -- -D warnings`
clean; `cargo test --lib` **3,304 passed**; every integration test green, the differential now **9
passed**.

**Pre-existing and untouched:** `cargo clippy --all-targets` still warns in
`examples/shared/stutter_model.rs`, `examples/shared/stutter_table.rs` and
`examples/ng_str_stutter_rate.rs`. None is in the library or the tests.

## Follow-ups

- **The consumer swap itself.** `GenomeSegments::counts()` exists and nothing calls it yet; the
  fourteen call sites that build a walk are the next job.
- **`segment_criteria.rs`'s stale gate archaeology** — the comment saying four of five rejection
  reasons are structurally zero. The measurement above says otherwise, and the fix is prose.
- **`segments_of_rows` and `locus_span` have no callers.** Noticed while threading the tally through;
  left for the review pass rather than deleted in passing.
