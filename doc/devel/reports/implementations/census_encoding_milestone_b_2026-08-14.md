# B5 — milestone B asserted: what the encoding moved, parameter by parameter

**Plan:** [census_rename_and_encoding.md](../../ng/impl_plan/census_rename_and_encoding.md) step B5.
**Compares:** the branch after B4 (`fb944a65`) against the state before B1 (`24eb4409`), both
built and run in this session.
**Date:** 2026-08-14.

---

## 1. The answer

**On the eight tomato accessions — 2.4 to 30.6 reads a position — not one fitted number moved.**
Not the error rates, not a sample's heterozygosity, not the frequency density, not a
contamination estimate, not a repeat-tract stratum row, and not the log-likelihood.

**On the GIAB trio — a few hundred reads a position — every fitted parameter is unchanged to the
precision the harness prints, and the log-likelihood moved by 759 nats.**

That is the split the specification predicts, and §3 says why it is even narrower than
predicted.

## 2. The two cohorts, parameter by parameter

Both runs at four container CPUs through `tmp/run_oracle.sh`, which pins the count because the
parallel reduction sums in a different order at a different one.

### 2.1 Tomato, eight accessions over twelve spans, 59,900 kept positions

| parameter | before B1 | after B4 |
|---|---|---|
| read group 0, error rate at an ordinary position | 0.00444 | **0.00444** |
| …at a mismapped one | 0.0965 | **0.0965** |
| positions mismapped | 0.0104 | **0.0104** |
| carrying only the reference | 0.9833 | **0.9833** |
| only a non-reference base | 0.00012 | **0.00012** |
| the segregating density | Beta(2.631, 17.222) | **Beta(2.631, 17.222)** |
| positions a sample carries an extra copy of | 0.00651 | **0.00651** |
| the share of the panel carrying one | Beta(0.869, 2.512) | **Beta(0.869, 2.512)** |
| expected heterozygosity | 3.621 /kb | **3.621 /kb** |
| every sample's heterozygosity, all eight | 0.000 to 0.192 /kb | **identical, all eight** |
| every sample's homozygote excess, all eight | 0.922 to 1.000 | **identical, all eight** |
| contamination, the four estimated | median 0.0127, highest 0.0564 | **identical** |
| positions more likely mismapped than not | 565 of 59,900 | **565 of 59,900** |
| every repeat-tract stratum row | 26 strata | **identical, all 26** |
| **log-likelihood** | 72,042 | **72,042** |

**The whole diff of the two runs is three things, none of them a fitted number**: the sparse
list's size in memory, the tract record's size in memory, and one table header.

### 2.2 The GIAB trio, three samples over 100 spans, 59,737 kept positions

| parameter | before B1 | after B4 |
|---|---|---|
| read group 0, error rate at an ordinary position | 0.00491 | **0.00491** |
| …at a mismapped one | 0.0481 | **0.0481** |
| positions mismapped | 0.0083 | **0.0083** |
| carrying only the reference | 0.9983 | **0.9983** |
| only a non-reference base | 0.00020 | **0.00020** |
| the segregating density | Beta(3.459, 6.224) | **Beta(3.459, 6.224)** |
| positions a sample carries an extra copy of | 0.00006 | **0.00006** |
| expected heterozygosity | 0.638 /kb | **0.638 /kb** |
| HG002 / HG003 / HG004 heterozygosity | 0.786 / 0.601 / 0.766 /kb | **identical** |
| homozygote excess, every sample | 0.000 | **0.000** |
| contamination | not identified, all three | **not identified, all three** |
| positions more likely mismapped than not | 442 of 59,737 | **442 of 59,737** |
| every repeat-tract stratum row | 32 strata | **identical, all 32** |
| reads that reached a tract without crossing it | 102,308 | **102,308** |
| **log-likelihood** | −859,732 | **−858,973** |
| positions above the cap | 99.9% | **99.8%** |

**The fit's path to those numbers differs early and not late.** At pass 1 the first sample's
heterozygosity reads 1.1807 /kb where it read 1.1802; by pass 25 the two runs print the same
0.7858 and stay together to pass 60. The only column still differing at the end is the
log-likelihood.

## 3. Why the movement is narrower than the specification predicts

The specification predicts change "confined to positions carrying more than 124 reads a
position". **Measured, the change is confined to positions carrying between 98 and 124**, and
above the cap the new encoding and the old agree exactly:

- **above 124 reads** the walk used to clip the depth to 124, whose bin read back as `124..=124`;
  now it keeps the true depth, whose bin is clamped to `124..=124` when a count is divided by
  it. Identical.
- **between 98 and 124** nothing was clipped, but the old rule read every code in the ladder's
  top bin back as exactly 124. Now it reads back as 98 to 124. **That is the correction**, and
  on the trio it is about 1 position in 1,000 — the occupancy table's 99.9% above the cap
  becoming 99.8%.

759 nats over 59,737 positions × 3 samples is **0.004 nats a position-sample**. The
log-likelihood is the only number that moves because it is the only one summed over every
position rather than fitted from them.

## 4. What the milestone bought, measured

Two of the four steps were free by construction (the ladder's ten new rungs cost no bits; the
true depth costs nothing to record). The two encoding steps shrink what a cohort holds:

| | before B1 | after B4 |
|---|---:|---:|
| HG004's sparse list, 47,552 entries | 0.571 MB | **0.380 MB** |
| HG004's whole census | 0.622 MB | **0.432 MB** |
| tomato's deepest accession, whole census | 0.139 MB | **0.108 MB** |
| a tract, on a census with 13 tracts a stratum | 25.0 bytes | **21.0** |
| a tract, on a census with 6.75 tracts a stratum (the trio) | 25.0 bytes | **25.2** |

**The tract record grows on a census too small to amortise a per-stratum count** — B4's report
gives the break-even, about seven tracts a stratum — and shrinks by a quarter at the scale the
census runs at. §5 reads that off the full cohort.

## 5. The confirming run on the full cohort

*(To be completed when the run lands: the 63 accessions over 8 Mb at a two-million-position
census, against the A1 baseline in `tmp/baseline_a1/tomato.txt`. It is the run that carries a
census of 1,999,404 positions and 4,164 tracts rather than the oracle's 59,900 and 349, so it
is where the per-tract and per-entry sizes are read at scale — and where the duplicated-copy
class and the contamination estimator have the panel they need.)*

**A caveat on that comparison, stated before it is made.** The A1 baseline predates the pinning
of the container's CPU count, so its fit trajectory may differ from this run's in the last
digits for that reason alone — measured at 615262.063054 against 615262.050187 at four CPUs
against six. **The converged parameters agreed at both counts**, so the comparison is on the
fitted numbers and not on the trajectory column.
