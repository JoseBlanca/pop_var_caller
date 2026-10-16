# The repeat-tract fit's allele window must reach past the recorded reads

**Status:** draft, 2026-10-12, for the owner's review. Nothing here is built. It amends
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1 and §8.1 (the ±6 decision) and
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §3.2 (the two widths), and
it touches calling's length prior ([`population_diversity.md`](population_diversity.md) §4.2,
`calling/genotype_prior/seed_ssr.rs`).

---

## 1. What this is

### 1.1 The problem

**The repeat-tract fit lets an allele sit at most 6 repeats from the reference length
(`ALLELE_SPAN = 6`, `parameter_estimation/joint/ssr_fit.rs:934`), but the census records how far a
read sits out to 8 (`RECORDED_OFFSET_RANGE = 8`, `parameter_estimation/joint/census.rs:428`).** The
design requires the opposite order. Both specs that set the two widths say the allele window "must
reach beyond" the recorded range: that is what lets the fit attribute a read in an end bucket to a
distant allele instead of to a long slip ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)
§4.1, [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §3.2). The rule held
when the recording was ±4. On 2026-08-20 the recording was widened to ±8, to stop it
under-measuring the slip rate of long homopolymers
([`str_slippage_level_curve.md`](str_slippage_level_curve.md) §11). The allele window stayed at ±6.

**What that does.** A chromosome 7 or more repeats from the reference has no allele class to sit in.
The fit puts it in the ±6 class instead. Its reads, which the census records at offsets up to ±8, are
then explained the only way left: as slips of one or two repeats from ±6.
- The ±6 class collects every allele that far or further.
- The fitted stutter rises to explain the reads beyond it.

**Measured.**
- **kimura, 2,169 tomato accessions** (`tmp/str_shapes_kimura_2026-10-11.md` on kimura). In
  homopolymer classes from 15 repeats up, the share at −6 is larger than the share at −5: 0.070
  against 0.013 at 22 repeats, and 0.114 against 0.021 at 23. Both classes have 100 or more tracts.
- **tomato2, 59 accessions over 32 Mb.** Refitted with the window at ±10 instead of ±6 (the counting
  fit, `ssr_fit/composite_fit.rs`), the slip rate:

| class | tracts | slip rate, window ±6 | slip rate, window ±10 |
|---|---|---|---|
| homopolymer, 8–15 repeats | 88–4,947 | — | within 3% of ±6's |
| dinucleotide, 6–15 repeats | 20–426 | — | within 3% of ±6's |
| dinucleotide, 18 repeats | 25 | 0.299 | 0.223 |
| dinucleotide, 21 repeats | 13 | 0.376 | 0.167 |
| dinucleotide, 24 repeats | 11 | 0.235 | 0.071 |
| dinucleotide, 25 repeats | 21 | 0.528 | 0.228 |
| homopolymer, 21 repeats | 16 | 0.107 | 0.092 |

  How far a slip reaches also moves where tracts are short. At homopolymer 8 repeats (4,947 tracts)
  the fall-off is 0.19 at ±6 and 0.11 at ±10: reads 7 or 8 repeats short had been forced to count as
  long slips.

**The ±6 decision itself was sound for the data it was measured on.** It cleared a threshold: the cost
of loci outside the window is nothing up to about 8% of loci outside, and wrecks the fit beyond about
19% ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §8.1). That was a fraction over a whole
cohort, where most tracts are short. A long homopolymer class on its own sits in that 8–19% band:
kimura's 23-repeat class has 0.114 of its alleles at −6 or further.

### 1.2 What changes

1. **The allele window becomes ±9: the recorded range plus one.**
2. **The two outermost classes mean "this far or further".** The class at −9 holds every allele 9 or
   more repeats shorter, and +9 every allele 9 or more longer.
3. **The window is tied to the recorded range in code**, so widening one cannot again leave the other
   behind.

**Why ±9, and not wider.** The census saturates at ±8: a read from an allele 9 repeats short and one
from 15 repeats short both land in the −8 bucket. Alleles at ±9 and beyond therefore produce
identical reads, so the data cannot tell them apart. One class per side holds all of them, and every
further class would be a number the reads cannot place. ±9 is the widest window the recording
supports.

**Why not narrow the recording instead.** The recording was widened from ±4 for a measured reason:
±4 under-measured the slip rate of 30-repeat homopolymers 2.26-fold
([`str_slippage_level_curve.md`](str_slippage_level_curve.md) §11).

---

## 2. What it touches

- **The repeat-tract fit** (`ssr_fit.rs`): `ALLELE_SPAN` becomes `RECORDED_OFFSET_RANGE + 1`. Every
  class carries 19 allele lengths where it carried 13, so 190 genotypes a sample where it carried
  91. Each class's likelihood table and the climb's time per round grow by about 2.1 times.
- **The counting fit** (`ssr_fit/composite_fit.rs`, a prototype): takes the window as a setting
  already, and costs little more at ±9.
- **The parameters file:** each class's `shares_by_repeat_offset` has 19 entries where it had 13. The
  format already allows any odd count of at least three (`calling/parameters_file/validate.rs:1312`),
  and calling reads the window from the count (`stratum_fits.rs:328`), so no version change. Files
  written at ±6 still read.
- **Calling's length prior** (`calling/genotype_prior/seed_ssr.rs:288-300`): today a candidate length
  beyond the window takes a floor weight, because the end class means "exactly ±6". Under §1.2 the
  end class means "this far or further". Its weight therefore belongs to every candidate at or
  beyond ±9, shared among them (open question §5.2).
- **Numbers that move:** every repeat-tract number of a cohort whose alleles reach 7 or more repeats
  from the reference, and the calls made with them. The oracle's checksums move; the move is
  measured and explained before they are re-recorded.

---

## 3. What will bite

- **The usual fit's average over tract frequencies gets worse with more allele classes.** At 256
  points it is already several log-likelihood units a tract too low at 13 classes, and the error grows
  with the class count (`reports/implementations/fit_precision_quadrature_average_2026-10-01.md` §3:
  at 20 samples, at most 0.03 a tract too low at 3 classes and up to 8.9 at 13). At 19 classes it has not been measured,
  and it will be worse. This change and the counting fit (§5.1) are therefore not independent.
- **The edge classes are catch-alls, not ordinary lengths.** Anything that treats the shares as a
  distribution over single lengths must handle them apart. That covers a shape fitted across lengths
  (the route-2 research) and calling's prior (§2). A falling tail cannot follow an edge that holds
  0.40 of a class, as tomato2's dinucleotide 25 does at ±10.
- **Some of the edge mass may not be alleles.** A read counts towards a tract only if it crosses the
  whole tract, and a longer allele is crossed by fewer reads. A long tract's reads therefore
  over-represent its shorter alleles. That bias is not measured, and this change does not address
  it (§5.3).
- **The guard against stale constants.** `ALLELE_SPAN`'s doc comment still says the recorded range is
  ±4 (`ssr_fit.rs:931`). Replacing the number without tying it to `RECORDED_OFFSET_RANGE` repeats the
  defect at the next widening.

---

## 4. Decided by the owner

- **2026-10-12:** fix the window first, before the counting fit replaces the climb and before the
  library multipliers resume.
- **2026-10-12, §5.1 settled:** widen the window in both fits now. Whether the counting fit replaces
  the climb is decided afterwards, on the climb's measured cost at ±9.
- **2026-10-12, §5.2 settled:** calling splits an edge class's weight equally among the candidate
  lengths at or beyond the edge.

---

## 5. Open

1. **Does the counting fit replace the climb first?** Settled 2026-10-12 (§4): both fits widen
   now; the order is decided on the climb's measured cost at ±9.
2. **How calling shares an edge class's weight.** Settled 2026-10-12 (§4): equally among the
   candidate lengths at or beyond the edge. With no candidate that far, the edge's weight goes
   unused, as a share's does today for a length with no candidate.
3. **The crossing bias for long alleles** (§3). Measure the chance a read crosses a tract of each length
   at the cohort's read length. Then decide whether the census should correct for it, or the fit
   should weight by it. **Not in this change.**
4. **Should the recording widen with tract length?** For the longest tracts, ±8 may itself be narrow:
   tomato2's dinucleotide 25 puts 0.40 at the edge even at ±10. **Not in this change**; it needs the
   crossing bias first.

---

## 6. How we know it works

1. **A drawn class with alleles 8 repeats from the reference** (unit test). At a window of ±6 it shows
   the defect: the slip rate is fitted high and the −6 share collects the far alleles. At ±9 the slip
   rate is recovered within its error, and the −9 class holds the far alleles' mass.
2. **A guard test:** the allele window is wider than the recorded range.
3. **tomato2 and kimura:** the long classes' slip rates at ±9 against ±6, reported as in §1.1; the
   short classes within their errors.
4. **The oracle cohort:** its checksums move. The move is measured with
   `scripts/promote_ng_oracle.sh` and explained before re-recording. The four accessions' classes are
   short, so the move is expected to be small.
5. **Calling:** a tract with a candidate 10 repeats short gets the edge's weight, not the floor (§5.2).
