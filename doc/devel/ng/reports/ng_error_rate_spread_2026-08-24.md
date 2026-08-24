# ng — how far apart are 63 libraries' error rates, and what borrowing one costs

**2026-08-24**, branch `ng-error-rate-routes`. Measured with
[`examples/ng_histogram_error_rates.rs`](../../../../examples/ng_histogram_error_rates.rs) on the
63-accession tomato cohort (`benchmarks/tomato1`, 8.0 Mb of BED in 80 spans, 2.5× to 28.6× a
position).

---

## 1. The answer

**A library that cannot fit its own error rate should keep borrowing its siblings' — and the reason
is not that borrowing is accurate, it is that the alternatives are worse.**

What the rate governs is the calibration scale the read likelihood divides into each read's own
error ([`read_likelihoods.md`](../spec/read_likelihoods.md) §3.2), so the question that decides it
is *how far does the average charged error end up from the library's own measured rate*:

| what a thin library gets | how far the charged average lands from its own rate |
|---|---|
| **the unweighted mean of the other libraries** — what the code does | **1.51×**, median over the 63 |
| the default of 0.001 | 2.99×, against the cohort's median fitted rate |
| no rescaling at all — a scale of one | ≈ **5×** |

## 2. Libraries differ by 15.6-fold, which is why the question looked open

Across 63 preparations of one crop on one instrument — the friendliest case borrowing will ever
get — the fitted marginal error rates run:

| | rate | Phred |
|---|---|---|
| lowest (`SRS3394642`) | 5.957 × 10⁻⁴ | 32.3 |
| median | 2.985 × 10⁻³ | 25.2 |
| highest (`SRS3394695`) | 9.271 × 10⁻³ | 20.3 |

**15.6-fold end to end**, and it is data rather than grid: the fit steps a ladder from Phred 10 to
50 in quarter-Phred rungs, so adjacent rungs differ by 6%.

Hand each library the unweighted mean of the other 62 and compare against its own fitted rate:
**median factor 1.51, worst 6.06** (`SRS3394642`, own 5.957 × 10⁻⁴ against a borrowed
3.609 × 10⁻³). **46 of the 63 would land within a factor of two of their own rate; 17 would not.**

That number on its own reads as an argument against borrowing, and it was taken that way for one
draft of this report. It is the wrong comparison, because it prices borrowing against a rate a thin
library by definition does not have.

## 3. Why the alternatives are worse, and by how much

**Read quality scores systematically overstate quality, and that is the whole reason the scale
exists.** The cohort's median fitted error rate is 3.0 in a thousand; the mean minted per-read error
over the same cohort is 6.0 in ten thousand
([`ng_prereq_closeout_two_averages_2026-08-24.md`](../../reports/implementations/ng_prereq_closeout_two_averages_2026-08-24.md)).
So a typical library's reads claim to be about **five times** better than the pre-pass finds them.

- **No rescaling** hands that factor of five straight through, in the direction that makes reads
  look cleaner than they are.
- **The default of 0.001** sits 2.99× from the cohort's median fitted rate.
- **Borrowing** is 1.51×.

Borrowing is better than the next best answer by a factor of two, and better than the one that
sounds most cautious by a factor of three.

**One library's scale, end to end, because a ratio of two medians is not a measurement of anything
in particular.** `SRS3394712` (`SRR7279481`): fitted marginal rate 4.591 × 10⁻³ against a mean
minted per-read error of 6.859 × 10⁻⁴ — **a scale of 6.7**. That is the first time both halves of a
calibration scale have existed on real data.

## 4. A larger finding from the same run: the second class of site is refused 6 times in 10

The pre-pass models a sample as a mixture of ordinary sites and a smaller class of mismapped,
noisier ones. **In 39 of the 63 accessions that class was refused outright** — the sample asked for
a noise level outside the range the model can represent, so it was fitted as though every site were
ordinary (`site_noise_off_the_ladder`). Where the class did fit, it covers 2.1 to 22.2 sites in a
thousand, median 6.5.

**Four of the 63 ran the coupled error-rate/genotype-frequency fit out of iterations** without
converging: `SRS3394687`, `SRS3394598`, `SRS3394641_SRR7279529`, `SRS3394560`. All four *did* get a
second class, so this is not the same failure.

The previously recorded evidence for this was two of three tomato accessions railing at the model's
ceiling
([`ng_noise_model_extension_n5_2026-08-10.md`](../../reports/implementations/ng_noise_model_extension_n5_2026-08-10.md)).
It is now measured over a whole cohort, and **it is the thing the census route's cohort-wide noisy
share is meant to fix**, which makes it the sharpest question the route comparison can answer.

**Nothing here is an artefact of a fit that failed quietly.** All 63 rates were fitted from the
library's own sites (`Provenance::FittedHere`, none borrowed or defaulted in this cohort), and none
railed at either end of the error-rate ladder.

## 5. What the run was

One walk per accession over the same 8.0 Mb of BED, typed against the reference's repeat catalog so
that repeat tracts go to the repeat path and only ordinary sequence reaches the accumulator. The fit
is the pre-pass's own public entry point,
[`estimate_generic_parameters`](../../../../src/ng/parameter_estimation/generic/estimate.rs);
nothing in the driver reimplements it.

**Three properties of the run a reader of these numbers needs.**

- **The inbreeding coefficient is supplied at zero, not fitted.** The runs model needs 3,000
  separate 100 kb windows each holding a site and the tomato BED touches about 80. The error rate is
  fitted jointly with the sample's genotype frequencies, so those frequencies were fitted under a
  stated `F`.
- **Every site is scored at ploidy 2.** `ConstantPloidy` is what the pre-pass builds; the BED is
  autosomal.
- **The rate is the marginal one** — the rate a read disagrees at a site drawn at random, which is
  what a sample emits. It is the only one of the three numbers in the mixture that both error-rate
  routes can produce, which is why the route comparison runs on it.

**The tomato repeat catalog did not exist before this run.** The typed-region stream reads a catalog
built beside the reference, and `$HOME/genomes` is mounted read-only. `repeat-catalog --output`
writes one anywhere, so it lives beside the CRAMs at
`benchmarks/tomato1/crams/S_lycopersicum_chromosomes.4.00.repeats.parquet` — 6,375,702 repeats over
13 contigs, 61 MB, and the reader checks it against the reference's contig digests wherever it sits.

## 6. What this does not settle

- **It is one cohort of one crop.** The spread across libraries of *different* preparations, or of a
  species with a worse reference, is unmeasured and could be wider.
- **The `Supplied` rung is untested here**, because no run supplied a rate.
- **The route comparison is not in this report**, and it is not settled either:
  [`ng_error_rate_routes_2026-08-24.md`](ng_error_rate_routes_2026-08-24.md) records what the cohort
  route gave on the same 63 accessions, and why the answer arrived as one number instead of 63.
