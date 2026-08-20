# ng — what splits when the contamination fraction moves to the read group

*Steps A1, A2 and the archive half of B2 of
[`../impl_plan/contamination_read_group_grain.md`](../impl_plan/contamination_read_group_grain.md).
Reading and counting only — **no estimator code was changed and nothing was fitted**. 2026-08-20.*

---

## 1. The answer in one paragraph

**Of the eleven quantities the contamination estimator computes, three change when the fraction
moves from the sample to the read group, and all three are read counts.** Everything that decides
*which positions are markers*, *what allele frequency each sample is judged against*, and *which
samples can be judged at all* is computed from a sample's reads pooled across its libraries and
stays exactly where it is. §6 of the plan is confirmed, with one correction in its favour: a fourth
quantity — the per-base error rate the fraction is scored with — is **already fitted per read group**
and is currently thrown away, so the split repairs an approximation rather than introducing one.

**The panel is not partitioned by this change, so the +0.015 measured for panel partitioning does not
apply.** §3 gives the mechanism. What the split does cost is stated in §4: the same number of markers,
divided reads.

---

## 2. Parameter by parameter

Line numbers are
[`joint/contamination.rs`](../../../../src/ng/parameter_estimation/joint/contamination.rs) unless
another file is named. "Sample-pooled" means *computed from all of a sample's reads, whichever library
they came from* — which is what the code does today by summing over a sample's sections.

| what it is | computed from | grain today | grain after | why |
|---|---|---|---|---|
| **which non-reference base the cohort carries at a position** (`major`, [L463–483](../../../../src/ng/parameter_estimation/joint/contamination.rs#L463)) | every read of every sample | position × cohort | **unchanged** | a property of the population at that position; no sample's reads are separated from another's in it |
| **whether a position is a marker at all** (`covered ≥ MIN_SAMPLES_WITH_DATA`, `MIN_FREQUENCY`, [L530–549](../../../../src/ng/parameter_estimation/joint/contamination.rs#L530)) | sample-pooled depth | position × cohort | **unchanged** | keeping this sample-pooled is what makes the two grains give the identical marker set, which is the regression guarantee of plan step D5 |
| **the panel-wide allele frequency** (`pooled`, [L546](../../../../src/ng/parameter_estimation/joint/contamination.rs#L546)) | every read of every sample | position × cohort | **unchanged** | it is also the frequency the *contaminant* is drawn against (spec §3.4.3: whoever else was on the plate), which is not this sample's property at all |
| **how likely a position is mismapped** (`cleanliness`, [L521–529](../../../../src/ng/parameter_estimation/joint/contamination.rs#L521)) | the joint fit, over the whole cohort | position × cohort | **unchanged** | but see §5 — this is the one place where read length can still bite |
| **each sample's expected number of non-reference copies** (`dosage`, [L551–569](../../../../src/ng/parameter_estimation/joint/contamination.rs#L551)) | sample-pooled counts | sample × position | **unchanged** | a genotype is the individual's, and both libraries of one plant carry it; pooling is not an approximation here, it is the right estimator |
| **each sample's coordinates on the panel's axes of variation** (`ancestry_coordinates`, [L402](../../../../src/ng/parameter_estimation/joint/contamination.rs#L402)) | the dosages above | sample | **unchanged** | ancestry is the individual's. **This is the whole frequency half of the model and it never sees a read group** |
| **how much of its own fitted frequency a sample supplies** (`coordinate_leverage`, [L403](../../../../src/ng/parameter_estimation/joint/contamination.rs#L403)) | the coordinates only | sample | **unchanged** | see §6 — this is A2 |
| **the straight line each position's frequency is in those axes** (`fitted_lines`, [L404](../../../../src/ng/parameter_estimation/joint/contamination.rs#L404)) | every sample's dosages | position × cohort | **unchanged** | fitted from all 50 samples before and all 50 after; **nothing here shrinks** |
| **the panel's spread along each axis** (`spread`, [L407–417](../../../../src/ng/parameter_estimation/joint/contamination.rs#L407)) | the coordinates | cohort | **unchanged** | |
| **homozygote excess** (`hom_excess[sample]`, [L433](../../../../src/ng/parameter_estimation/joint/contamination.rs#L433)) | the joint fit | sample | **unchanged** | inbreeding is the individual's; two libraries of one plant share it |
| **the reads the fraction is fitted from** (`alternative`, `depth_low`, `depth_high`, [L485–514](../../../../src/ng/parameter_estimation/joint/contamination.rs#L485), consumed at [L790–817](../../../../src/ng/parameter_estimation/joint/contamination.rs#L790)) | one sample's reads, **summed over its libraries** | sample | **→ read group** | this is the change |
| **the per-base error rate the reads are scored with** (`error_rate[sample]`, [L434](../../../../src/ng/parameter_estimation/joint/contamination.rs#L434)) | the joint fit | **read group already**, then discarded | **→ read group** | see §2.1 |

### 2.1 The error rate is already at the new grain and is currently thrown away

The joint fit fits the clean error rate **per read group** — `Parameters.clean` is documented *"per
read group, in the order `groups` lists them"*
([`fit.rs:992`](../../../../src/ng/parameter_estimation/joint/fit.rs#L992)), which is
[`parameter_prepass.md`](../spec/parameter_prepass.md) §1's table obeyed. The call into the
contamination estimator then does this
([`fit.rs:1278–1280`](../../../../src/ng/parameter_estimation/joint/fit.rs#L1278)):

```rust
let per_sample_error: Vec<f64> = (0..lent.len())
    .map(|s| parameters.clean[group_index[s][0]])
    .collect();
```

**A sample's first library's error rate is used for all of its libraries.** On the benchmark cohorts
that is exact, because every sample has one library. On a sample with two libraries at different
error rates it is wrong, and it is wrong in the direction that matters: the estimator reads a fraction
out of the gap between how many disagreeing reads a genotype predicts and how many are there, and the
error rate is what sets the prediction.

*(The same `group_index[s][0]` shorthand appears in the main fit's own scoring —
[`fit.rs:2050`](../../../../src/ng/parameter_estimation/joint/fit.rs#L2050),
[`2079`](../../../../src/ng/parameter_estimation/joint/fit.rs#L2079),
[`2145`](../../../../src/ng/parameter_estimation/joint/fit.rs#L2145),
[`2222`](../../../../src/ng/parameter_estimation/joint/fit.rs#L2222) — where a sample's reads are
scored at its first read group's rate while the *tallies* are spread over all of them. That is a
separate approximation in a different estimator and this plan does not touch it, but it is the same
shape and somebody should know it is there.)*

### 2.2 The evidence is already stored at read-group grain — nothing has to be plumbed

The census keeps one section per read group:
`SectionKey::Generic(ReadGroupId)`, and a sample's lent evidence is
`SampleGenericSections<'a> = Vec<(ReadGroupId, &GenericEvidence)>`
([`census.rs:1075–1131`](../../../../src/ng/parameter_estimation/joint/census.rs#L1075)). The
estimator discards the identifier in exactly two loops — `for (_, group) in sections` at
[L465](../../../../src/ng/parameter_estimation/joint/contamination.rs#L465) and
[L490](../../../../src/ng/parameter_estimation/joint/contamination.rs#L490) — and those two `_`
patterns are the entire extent of the module's read-group awareness. **No file format changes, no new
field is recorded, no walk is re-run.**

---

## 3. Why the +0.015 does not transfer

[`joint_contamination_2026-08-12.md`](joint_contamination_2026-08-12.md) §2 states what its three
arms varied: they differ *"only in the frequency each sample's genotype is scored against"*. The arm
that cost +0.015 on every sample's estimate, and put 41 to 47 of 50 clean samples over a 1% threshold
(§3's table there), estimated each subpopulation's allele frequency **from that subpopulation's twelve
members** instead of from all fifty. What shrank was the number of *samples* behind a *frequency*.

Under this change the frequency at a position is still fitted from every sample in the panel
(`fitted_lines`, §2's row eight), each sample still gets one set of coordinates from all of its own
reads (`ancestry_coordinates`, row six), and the panel is not divided into anything. **The number of
samples behind every frequency is identical before and after.** The +0.015 is a measurement of a
partitioned panel and this change does not partition the panel.

**The objection is written down in a current spec and should be struck.**
`read_likelihoods.md` §3.6 argues against the read-group grain in these words: *"Splitting one sample's reads by read group and fitting each separately would hand every fit
less data, and that document's own measurement says what too little data does here: partitioning a
panel manufactures contamination, adding about 0.015 to every sample's estimate."* The first clause is
true and is §4 below; the second names a different measurement. Plan step D3 already amends that
section.

---

## 4. What the split does cost, stated as what it is

**The marker set does not shrink and the reads behind one fraction do.** Marker selection is
sample-pooled (§2, row two), so a two-library sample's fraction is fitted over the same positions as
before — but at each of them, from that library's reads only. A sample whose reads are spread evenly
over `L` libraries gives each fraction about `1/L` of the depth per position that the per-sample
estimate had.

**That is the depth axis, which is already measured for this estimator, not a new unknown.** From the
module header and [`census_depth_resolution_2026-08-16.md`](census_depth_resolution_2026-08-16.md): a
genuinely 3%-contaminated sample reads 0.0263 at 30 reads a position with 39 clean samples at exactly
0.0000, and reads 0.0102 at three reads a position with the worst clean sample at 0.0003. So the
estimator still separates the contaminated sample from the panel at three reads a position — thirty
times the worst clean one — while understating the value by a factor of three. **Plan step C2 is
therefore a sweep along an axis whose two ends are already known**, and the question it answers is
where between them a library falls.

Two facts sharpen how far down that axis a real library goes. In the tomato archive a library is one
file (§5), so a two-library sample at tomato's three reads a position leaves each fraction about 1.5
reads a position — below anything this estimator has been measured at. And the one sample with 42
libraries would leave each fraction a fortieth of the sample's depth, which no threshold can survive;
that is what plan step C3 exists to decide.

**One small gain in the other direction.** Above eight reads a position the census stores a depth
*range* rather than a number, and the estimator sums over it — the width of that range was the whole
cause of a floor of 0.025 on a panel with nobody contaminated
([`census_depth_resolution_2026-08-16.md`](census_depth_resolution_2026-08-16.md)). Today a
two-library sample has its two ranges **added** ([L500–503](../../../../src/ng/parameter_estimation/joint/contamination.rs#L500)),
so the uncertainty widens; under the split each fraction sees its own library's narrower range. The
size of that gain is not measured here.

---

## 5. Read length: the position set is already common, and the residue is elsewhere

Plan §6 asks whether the set of positions each fraction is fitted over has to be held common between
two libraries of different read length. **It already is, by construction.** Each read group's section
is a depth code for every kept position in the same order
([`census.rs`](../../../../src/ng/parameter_estimation/joint/census.rs), `PackedDepthCodes` over the
kept-position list), and `fit_alpha` drops the positions a read group has no read at, once, before
the search ([L784–787](../../../../src/ng/parameter_estimation/joint/contamination.rs#L784)). A short
library covers fewer positions; it does not get a different list of them. **Nothing has to be built to
hold the position set common, and holding it common costs nothing.**

**The residue is the mismapping exclusion, not the position set.** A position is dropped as mismapped
on a probability the joint fit computes **for the position across the whole cohort**
(`cleanliness`, §2 row four). Short reads anchor worse than long ones, so a position that is mismapped
*only in the 101 bp library* keeps a low cohort-wide probability, survives as a marker, and puts
unexpected reads into that one library — which is the contamination signature, in exactly one read
group of one sample. **The split gives the estimator the resolution to see that and no way to tell it
from real contamination of that library.** This is what plan step C4 must measure, and the check it
describes — two simulated libraries differing only in read length, with equal planted fractions — is
the right one.

**How often it can happen, measured (§6.1 below): 15 of the 157 multi-library samples in the tomato
archive have libraries at different modal read lengths** — about 9 samples in 1,000 of the archive.

---

## 6. A2 — the refusal criterion transfers unchanged, and one new refusal is needed

**`MAX_LEVERAGE` refuses a sample, and it should keep refusing the sample.** How much of its own
fitted allele frequency a sample supplies is computed from the coordinates alone
(`coordinate_leverage`, [L618–644](../../../../src/ng/parameter_estimation/joint/contamination.rs#L618)) —
one number per sample for the whole run, before a single position is fitted. The coordinates are
sample-pooled and stay so, so the number is the same before and after, and a sample whose fitted
frequency is mostly its own echo is unfittable **in every one of its libraries at once**. Refusing all
of them together is the right answer and it needs no new code: the check at
[L422](../../../../src/ng/parameter_estimation/joint/contamination.rs#L422) simply runs per sample and
its verdict is copied to that sample's read groups.

**What has no equivalent today is a per-read-group refusal for want of data.** The current refusals
are all panel-wide or sample-wide: fewer than two samples (`NoPanel`), fewer than 100 markers in the
whole panel (`TooFewMarkers`, [L392–399](../../../../src/ng/parameter_estimation/joint/contamination.rs#L392)),
and the leverage refusal. **Nothing anywhere checks whether one sample has enough reads**, because
until now a sample either had the panel's markers or the panel had none. A read group covering four
markers would today be given a number from those four. That is plan step C3, and it is not optional.

### 6.1 The archive half of B2 — how many libraries, how far apart, how much data each

Counted from `tmp/tomato_manifest.tsv`, the 2,085-file manifest produced by
[`rick_sample_manifest.sh`](../../../../benchmarks/ssr_tomato1/scripts/rick_sample_manifest.sh) over
the 68-project tomato archive. This reproduces
[`read_groups.md`](../spec/read_groups.md) §1's counts exactly, which is why the rest of the table can
be trusted.

| | count |
|---|---:|
| files | 2,085 |
| samples | 1,707 |
| samples with one library | 1,550 |
| samples with two | 133 |
| samples with three | 20 |
| samples with 7, 16, 16 and 42 | 4 |
| **multi-library samples** | **157** |
| …whose libraries differ in modal read length | **15** |
| …whose libraries span more than one project | 0 |
| …whose libraries differ in declared platform | 0 |
| files per library, over the multi-library samples | median 1, max 2 |

The differing read-length pairs are 150/151 nine times, and then 101/301, 101/301, 54/100, 101/251,
76/101/114 and 100/251. **The 150-against-151 pairs are not a read-length difference in any sense
that matters**; the six others are, and one of them is the 42-library sample.

**A "library" in this archive is a sequencing run.** `LB` values are `<project>_<run>` strings written
by a re-headering script — 94% of them, per [`read_groups.md`](../spec/read_groups.md) §1 — and both
local benchmark files confirm the shape: `@RG ID:SRR5079859 SM:SRS1839214 LB:PRJNA353161_SRR5079859`.
**This makes the grain more right rather than less**: index hopping, which is the mechanism that makes
two libraries of one plant differ, happens within a sequencing run, and run is the unit these
identifiers carry.

**Median one file per library means the split divides a sample's reads by its library count**, with no
library carrying most of them. That is what §4's arithmetic assumed and it is now checked rather than
assumed.

### 6.2 What the archive half of B2 could **not** answer here, and what it would take

**Whether libraries of one sample come out with different fractions is unanswered.** It needs the
estimator run over each library separately, which needs the archive's alignment files; those are on
the Linux box the survey was run from and are not on this machine. What is here is 51 tomato accessions
(`benchmarks/ssr_tomato1/crams`) and 63 (`benchmarks/tomato1/crams`), **every one of them a single
library**, so nothing local can be asked this question.

Running it needs a panel, not a sample: the estimator has no answer for one sample
(`NotIdentified { NoPanel }`), and the fraction is read out of how a sample differs from the panel's
fitted frequencies. So the run is *one cohort of the archive's samples, with the multi-library ones
in it*, not 157 separate runs.

---

## 7. What this changes about the plan

- **The instrument for plan step B1 is the wrong file in the plan.** §7 names
  [`ng_joint_contamination_harness.rs`](../../../../examples/ng_joint_contamination_harness.rs), which
  is a **standalone re-implementation** of the estimator — its own `fit_alpha`, `pc_lines`,
  `coordinate_leverage` — and produced the frequency-arm numbers of
  [`joint_contamination_2026-08-12.md`](joint_contamination_2026-08-12.md) §3 and §4. Extending it
  would measure a model of the estimator rather than the estimator.
  [`ng_joint_contamination_control.rs`](../../../../examples/ng_joint_contamination_control.rs) draws
  a cohort with contamination planted at a chosen fraction, writes it as the census records a walk
  would have written, and runs the shipped `fit_jointly` and `fit_contamination` on it — its own
  header says *"what is compared is the shipped code under different settings, not a
  re-implementation of it"*. **That is the one to extend**, and the extension is one read group per
  library where it now writes `SectionKey::Generic(ReadGroupId(0))`
  ([L409](../../../../examples/ng_joint_contamination_control.rs#L409)).

- **Plan step D5's "no change at all" is a guarantee of the design, not a hope.** Because marker
  selection, the frequencies, the coordinates and the leverage are all sample-pooled and unchanged, a
  sample with one read group is scored on byte-identical inputs to today's. The benchmark result is
  not merely expected to be unchanged; a difference would mean the split leaked into one of §2's
  unchanged rows.
