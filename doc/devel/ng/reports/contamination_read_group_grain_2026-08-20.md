# ng — the contamination fraction at the read group, and what splitting it costs

*Steps B1, C1 and C2 of
[`../impl_plan/contamination_read_group_grain.md`](../impl_plan/contamination_read_group_grain.md),
with the decomposition it rests on in
[`contamination_grain_decomposition_2026-08-20.md`](contamination_grain_decomposition_2026-08-20.md).
Everything below is the shipped estimator under two settings, on panels drawn as the records a walk
would have written. 2026-08-20.*

---

## 1. What was measured, and the one sentence

**A plant's DNA is prepared into libraries, and a library is what a second plant's DNA gets into.**
Two libraries made from one plant can therefore carry different amounts of it, which one number for
the whole plant cannot say. The estimator now fits one fraction per read group; both grains are
selectable in one build ([`ContaminationGrain`](../../../../src/ng/parameter_estimation/joint/contamination.rs)).

**The finding: splitting the fraction to the library costs almost nothing, because what the split
takes away is not what limits the estimate.** A library holding three reads a position returns 0.026
when it is the plant's only library, 0.046 when it is half of a six-read plant, and 0.057 when it is
a quarter of a twelve-read plant — the same reads in the library being measured, against a planted
0.060. What sets the accuracy is the depth of everything the split leaves alone: the plant's
genotype and the panel's fitted allele frequencies, which keep every read the plant produced.

---

## 2. How the panels were drawn

[`ng_joint_contamination_control.rs`](../../../../examples/ng_joint_contamination_control.rs) draws a
cohort and hands it to `fit_jointly` and `fit_contamination` — **the shipped code, not a model of
it.** 400,000 positions, 40 plants in four subpopulations at `F_st` 0.20, 3.3% of positions mismapped
at an error rate of 0.024 against the ordinary 0.003, and the markers those positions produce dropped
as the estimator ships. About 8,800 of the 400,000 positions survive as markers in every cell.

**Each plant's reads are divided equally between its libraries and its total depth is held
constant**, so what a row varies is how one plant's reads are divided, not how many it has. One
plant — always sample 0 — has **6% of its first library's reads** coming from another plant of the
panel, drawn afresh per read. Every other library of every plant is clean.

`LIBRARIES`, `LIBRARY_ALPHAS`, `LIBRARY_DEPTHS`, `SEED` and `SWEEP` set all of this; with `LIBRARIES`
unset the program draws what it always drew, generator call for generator call, and reproduces the
table already published from it.

---

## 3. The sweep

Two seeds where both finished. **The contaminated library always carries 0.060 of its own reads**; the
column after it is what the clean libraries of the same plant returned.

| plant's depth | libraries | reads a position in the contaminated library | **by library** | by plant | worst clean library | worst clean plant |
|---:|---:|---:|---|---|---:|---:|
| 3 | 1 | 3 | 0.0261, 0.0263 | *the same number* | 0.0024, 0.0011 | 0.0024, 0.0011 |
| 3 | 2 | 1.5 | 0.0267, 0.0235 → 0.0000 | 0.0119, 0.0164 | 0.0043, 0.0048 | 0.0012, 0.0028 |
| 3 | 4 | 0.75 | 0.0186, 0.0213 → 0.0000 | 0.0069, 0.0065 | 0.0098, 0.0045 | 0.0007, 0.0020 |
| 6 | 1 | 6 | 0.0438, 0.0446 | *the same number* | 0.0000, 0.0000 | 0.0000, 0.0000 |
| 6 | 2 | 3 | 0.0455, 0.0381 → 0.0000 | 0.0216, 0.0223 | 0.0042, 0.0040 | 0.0007, 0.0003 |
| 6 | 4 | 1.5 | 0.0394, 0.0502 → 0.0000 | 0.0101, 0.0111 | 0.0101, 0.0087 | 0.0000, 0.0000 |
| 12 | 1 | 12 | 0.0534, 0.0529 | *the same number* | 0.0000 | 0.0000 |
| 12 | 2 | 6 | 0.0512 → 0.0000 | 0.0243 | 0.0004 | 0.0000 |
| 12 | 4 | 3 | 0.0573 → 0.0000 | 0.0114 | 0.0052 | 0.0000 |
| 30 | 1 | 30 | 0.0521 | *the same number* | 0.0000 | 0.0000 |
| 30 | 2 | 15 | 0.0545 → 0.0000 | 0.0251 | 0.0000 | 0.0000 |
| 30 | 4 | 7.5 | 0.0607 → 0.0000 | 0.0114 | 0.0000 | 0.0000 |

*"→ 0.0000" is what every clean library of the contaminated plant returned — exactly zero in all
twelve multi-library cells.*

### 3.1 The value follows the plant's depth, not the library's

Read down the third column and find the rows where **the contaminated library holds the same three
reads a position**:

| | plant's depth | libraries | fitted |
|---|---:|---:|---:|
| the library is the whole plant | 3 | 1 | **0.0261** |
| …half of it | 6 | 2 | **0.0455** |
| …a quarter of it | 12 | 4 | **0.0573** |

The reads the fraction is fitted from are identical in all three; the answer more than doubles, and
the third is within a twentieth of the planted 0.060. **A library's own read count is not what limits
this estimate.** The genotype the stray reads are surprising, and the allele frequency that genotype
is drawn against, are both fitted from every read the plant produced and from the whole panel — and
neither of those changes grain. At three reads a position a plant's genotype is barely known, and it
is that uncertainty, not the read share, that attenuates the fraction.

**This is [`contamination_grain_decomposition_2026-08-20.md`](contamination_grain_decomposition_2026-08-20.md)
§3's argument turned into a number.** The measured cost of partitioning a *panel* to fit frequencies
was about +0.015 on every sample; the cost of partitioning one plant's *reads* to fit a fraction is,
over most of this table, smaller than the seed-to-seed spread. *(Plant depth and panel depth move
together in this sweep — every plant is drawn at the same depth — so it does not separate "this
plant's genotype is better known" from "the panel's frequencies are better fitted". Both are things
the split leaves alone, which is the point being made.)*

### 3.2 One number per plant is the average, and it is wrong by the library count

At two libraries the plant grain returns 0.012 to 0.025 where the dirty library carries 0.060; at
four, 0.0065 to 0.0114. **It divides by how many libraries the plant has**, which is what an average
over reads must do. A plant with one contaminated library in four therefore looks about four times
cleaner than it is, and a threshold set on the plant will miss it.

### 3.3 What it costs: a higher noise floor, and only at low depth

The worst *uncontaminated library* anywhere in the panel sits above the worst *uncontaminated plant*,
because a library holds fewer reads and there are more of them to be the worst:

| plant's depth | worst clean library | worst clean plant |
|---:|---:|---:|
| 3 | 0.0043 – 0.0098 | 0.0007 – 0.0028 |
| 6 | 0.0040 – 0.0101 | 0.0000 – 0.0007 |
| 12 | 0.0004 – 0.0052 | 0.0000 |
| 30 | 0.0000 | 0.0000 |

**Above about twelve reads a plant the floor is back to zero.** Below it the floor rises, and the
weakest corner of the whole sweep is a plant at three reads a position divided into four libraries —
**0.75 reads a position each** — where the contaminated library reads 0.0186 to 0.0213 against a floor
of 0.0045 to 0.0098. That is a separation of two- to fourfold, against eleven-fold at twelve reads and
four libraries. It is still the right answer and it is the thinnest evidence in the table.

**Tomato's three reads a position with a single library is the common case and is unaffected** — one
library a plant gives the two grains the identical number, asserted as equality rather than closeness
by `one_library_a_plant_gives_the_two_grains_the_same_number`.

---

## 4. A library nobody could measure comes back clean, and says so

**The default assumption is that nothing is contaminated, and too little evidence lands near it**
(owner, 2026-08-20). With few markers covered the likelihood barely moves with the fraction, and the
check against zero the search ends with then keeps zero. A library drawn at a fortieth of its plant's
reads — half a read a position — **read 0.0080**, below a 1% threshold, having seen 2,875 of the
panel's markers against 7,669 for the full-depth libraries
(`a_library_with_almost_no_reads_reads_clean_and_says_how_little_it_saw`).

**So there is no separate refusal for a thin library, and there does not need to be — but the
evidence travels with the number.** Every fitted fraction now carries `markers_with_reads` and
`reads_on_markers`, its own counts, and `source`, saying whether it was fitted from that library's
reads or from every read of the plant. **A library measured clean and a library nobody could measure
both read near zero and are told apart by those counts**, which is what plan §9 asks for and what a
bare zero cannot give.

---

## 5. What changed in the code

- **[`joint/contamination.rs`](../../../../src/ng/parameter_estimation/joint/contamination.rs)** fits
  one fraction per read group. `Units` says which of a plant's sections pour into which fraction, so
  **one gather serves both grains** and the comparison in §3 is exact rather than two runs that might
  differ for another reason. The gather is now two passes: the first pools each plant's libraries over
  every kept position, as before, to choose the markers and fit the frequencies; the second fills the
  per-library counts at the surviving markers only, so the memory cost is one number per library per
  *marker* rather than per *position*.
- **The error rate is taken from the read group.** The joint fit has always produced it at that grain
  and this module used the plant's *first* read group's rate for all of them
  ([`fit.rs:1278`](../../../../src/ng/parameter_estimation/joint/fit.rs#L1278)) — exact only because
  every benchmark plant has one library.
- **`JointFit::contamination` is now a list per sample**, one entry per read group.
- **The `estimate-contamination` report has a row per library**, carrying the `@RG` id and library
  name, the two evidence counts, and which reads the number was fitted from.

**What did not change**: which positions become markers, the allele frequency at each, each plant's
place on the panel's axes of variation, the line fitted through those places, the leverage refusal,
and the homozygote excess. All are computed from a plant's reads pooled over its libraries, and
pooling is the right estimator for each of them.

---

## 6. What this does not answer

- **A library that mismaps more than its sibling.** Plan step C4 asks for a read-length check, and
  **the plan has the mechanism wrong** (owner, 2026-08-20). The read-length trap
  [`../spec/parameter_prepass.md`](../spec/parameter_prepass.md) §5 records belongs to the STR path:
  a library of 100 bp reads spans fewer long repeat tracts than one of 250 bp, so its slippage rate is
  measured over a different mix of loci. **Contamination is fitted from ordinary positions and has no
  such mix**: every position is a position, the census gives every read group a depth code for the
  same list of them, and a library with shorter reads simply has no read at more of them — which
  `fit_alpha` already drops, per read group, before it searches.

  **What is real is mappability, of which read length is one cause among several** (single-ended
  reads, insert size, a different aligner). A library whose reads anchor worse carries more reads from
  elsewhere in the genome, and a read from elsewhere carries an allele the plant should not have —
  the contamination signature. The defence is §3's mismapping exclusion, and **it is a per-position
  probability fitted across the whole cohort**, so a position that only misbehaves for one library
  keeps a low cohort-wide probability, survives as a marker, and puts stray-looking reads into that
  library alone. **So the experiment worth running plants extra unexpected reads in one library of a
  clean plant** — which the drawing program can express, where read lengths are not something it
  models at all. Until it is run, a fraction that stands out in exactly one library of a plant should
  be read as *either* contamination *or* worse mapping in that library.

  For scale on how often two libraries of one plant could differ this way at all: 15 of the 157
  multi-library tomato samples have libraries at different modal read lengths, and nine of those
  fifteen differ only as 150 against 151.
- **Real multi-library data.** No cohort under `benchmarks/` has any: every sample of tomato1, tomato2
  and the GIAB trio was sequenced from one library. The owner ruled the archive run out of scope on
  2026-08-20 — the grain follows from how contamination happens, not from a measurement.
- **What the caller does with it.** The number exists to change genotypes, and `read_likelihoods.md`
  §3.6 consumes a fraction per read group already; this removes the approximation recorded there, and
  nothing yet measures the effect on calls.

  **That document is not amended yet, and the reason is branch topology rather than oversight.** It
  was an untracked draft in the main checkout when this work began and has since been committed on
  `ng-str-slippage-curve`, so it is not on this branch to edit. Bringing a copy across would make two
  branches add the same 154 kB file independently, which merges as a whole-file conflict rather than
  as hunks. **Three things in its §3.6 have to change once that branch lands**, recorded here so they
  are not lost:

  1. *"It gives one number per sample"* — it now gives one per read group, and so do
     [`../spec/parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §3.1 and §3.4.
  2. **The approximation it records** — *a per-sample estimate applied to every read group of that
     sample* — is gone, which is the thing that section asked to be told.
  3. **The objection it raises is wrong and should be struck rather than softened.** It argues that
     splitting a sample's reads by read group *"would hand every fit less data"* and cites the +0.015
     that partitioning a panel costs. That measurement is about splitting the *panel* to estimate
     *frequencies* from a twelfth of it; nothing here splits the panel, and §5 above lists what does
     and does not change grain.
