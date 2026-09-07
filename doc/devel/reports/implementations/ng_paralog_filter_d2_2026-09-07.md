# The hidden-duplication filter — D2: six samples, and where the time goes

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone D, step D2
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §4, §5, §8
**Branch:** `ng-paralog-filter`

## The answer

**The filter removes 36 records of 2,311, and they are not scattered: 34 of the 36 sit inside four
stretches of DNA totalling 1,227 bases, out of the 199,672 the run called.** A reference-collapsed
duplication is a region, so a run of adjacent flagged positions is its signature. Drawing 36 of the
run's 2,311 written records at random 200,000 times, **not one draw was this clustered**: a random
36 puts 12.6 of its members into multi-member clusters on average, against 34 here.

**Scoring is 5% of the filtered part of the run**, so spec §8's parallel scoring is not worth a plan
at this size. **The constraint at the top of the range is disk, not time**: the spill would need
about 350 GB at 3,000 samples and five million records, against about 14 hours of one core.

| | |
|---|---|
| records without the filter | 2,311 |
| records scored | 2,311 — none unscorable, and **no record had no alternative allele** |
| records removed at `--paralog-fdr 0.01` | **36** (1.6 in 100) |
| fitted duplication rate `π` | 0.020634, **fitted from this run**; the EM converged |
| the cut | likelihood ratio 6.2500 |
| samples with a coverage model | **6 of 6** — none rejected, for any reason |
| ratios past the histogram's ±100 ends | **4, all at the negative end** |

```text
##paralogFilter=target_fdr=0.0100;pi=0.020634;lr_cut=6.2500;em_converged=true;samples_with_coverage_model=6/6
```

The off run still reproduces the pre-filter baseline: 2,311 records, sha256
`84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d`.

## ⚠ The six accessions are five plants

**`SRS3394712` and `SRS3394712_SRR7279484` are two sequencing runs of one biosample.** The tomato1
cohort contains biosamples sequenced more than once, and every run's CRAM carries the same sample
name; [`rename_dup_samples.sh`](../../../../benchmarks/tomato1/scripts/rename_dup_samples.sh)
rewrites the second and later runs to `<sample>_<run>` so a cohort VCF does not see duplicate
columns — its own header says so. Six of the 63 accessions in the full cohort are such pairs.

**This was not known when D2 was designed and it changes two things.** The plan's standing
"six-accession slice" is **five plants, one of them sequenced twice**, so the effective cohort at
every locus is five; and where the two runs of that plant agree, that is a consistency check on the
pipeline rather than two samples corroborating each other. It is stated here rather than worked
around: the ground is what the plan's baseline was recorded on, and moving it would move the sha256
every step in this plan checks against.

## The cohort's depth range, and what the fit made of it

The run now prints what each sample's fit came to. Beside it, that sample's **own median window
depth over the 2,311 records the run wrote** — a model-free figure, taken from the run's window dump:

| sample | fit: one copy, reads a window | fit: σ₀ | its own median window depth |
|---|---|---|---|
| SRS3394712 | 5.22 | 0.300 | 7.23 |
| SRS3394606 | **29.24** | 0.174 | 26.08 |
| SRS3394713 | 15.36 | 0.228 | 11.36 |
| SRS3394712_SRR7279484 | **3.02** | 0.302 | 3.96 |
| SRS3394714 | 20.66 | 0.158 | 20.32 |
| SRS3394711 | 13.09 | 0.183 | 14.72 |

**The measured depth range across the cohort is 6.6×** (26.08 against 3.96). The fitted one-copy
values span 9.7×, but that is a ratio of two fitted parameters and not a statement about the reads;
the two differ because the fit anchors on the mode of a sample's whole depth histogram while the
median above is over written variant sites only, which are a biased subset and biased differently
per sample. **All six fits were accepted**, which is the result: `CLAUDE.md` asks that a method work
from a few reads a position to several hundred, and one cohort of six already spans a factor of six
and a half.

**σ₀ falls as depth rises** — 0.300 and 0.302 at the two shallowest against 0.158 and 0.174 at the
two deepest. Over six points that is Spearman −0.886, two-tailed p = 0.033 on all 720 permutations:
real at the 5% level with no margin, and **partly arithmetic rather than discovery**, since σ₀ is a
scatter of *relative* depth and so carries a 1/√depth counting term by construction.

Two per-sample windows in a hundred are absent (277 of 13,866 sample-record pairs, 2.00%).

## What it removed, and why it looks real

**36 records, all biallelic SNPs.** The run wrote 2,217 biallelic SNPs, 27 multiallelic sites, 22
deletions, 21 repeat tracts, 15 insertions and 9 equal-length multi-base substitutions; **nothing
but a biallelic SNP was flagged.** Counting the flagged by kind and asking what that says about
spec §8's tract-aware allele term is step D4's job, and these counts are too small to carry it
alone.

**Where they are:**

| flagged records | over | at |
|---|---|---|
| 14 | 547 bp | 13,808,178–13,808,724 |
| 8 | 90 bp | 13,889,118–13,889,207 |
| 7 | 149 bp | 3,411,867–3,412,015 |
| 5 | 441 bp | 3,471,149–3,471,589 |
| 1 | — | 3,459,371 |
| 1 | — | 13,887,809 |

**34 of the 36 fall in four stretches totalling 1,227 bases of the 199,672 the run called.** The
median gap between adjacent flagged records is 19 bp. **Against the right null** — 36 drawn at
random from the 2,311 records the run actually wrote, which are themselves not uniformly spread —
200,000 draws never reached 34 in multi-member clusters; the mean is 12.6 and the highest 28.

**The flagged records are nine times more heterozygous than the rest.** 67 of 216 genotype calls at
flagged records are heterozygous, **31.0%**, against 474 of 13,361 — **3.55%** — at kept records, so
**8.7 times**. Restricted like for like, to the 2,018 kept records where all six samples are called
and all six windows are present (the flagged set has no no-calls and no absent windows), kept falls
to **2.70%** and the ratio is **11.5 times**. Heterozygotes are counted by the two alleles
differing, not by the string `0/1`, so the 21 multiallelic heterozygotes are included.

**Their coverage sits where the model says.** Measured against each sample's own median window depth
— model-free, and not the model's copy number — the sample-record pairs at flagged records run to a
median of **1.22** and a ninetieth percentile of **2.71**, against **1.00** and **1.50** at kept
records. So flagged records are not uniformly over-covered: most of their samples sit at ordinary
depth and the ones carrying the duplication sit near three times theirs.

**The largest cluster, sample by sample.** The 14 records at 13,808,178–13,808,724, each sample's
median window depth there against its own usual depth:

| | SRS3394712 | SRS3394606 | SRS3394713 | …_SRR7279484 | SRS3394714 | SRS3394711 |
|---|---|---|---|---|---|---|
| × its own usual depth | **2.43** | 1.22 | **1.93** | **2.75** | 0.83 | 0.99 |
| genotype | 0/0 | 0/0 | **0/1** | 0/0 | 0/0 | 0/0 |

**Three samples are locally over-covered, and the heterozygous one is among them** at 1.93 times its
own usual depth. That is the coincidence the whole model rests on — the heterozygous sample being
the over-covered one — visible in the output. The other two over-covered samples are the two runs of
the duplicated plant, so they are one plant, not two.

**An earlier draft of this report said the opposite**, that the over-covered samples were homozygous
reference and the heterozygous sample was at ordinary depth, and called it the one thing that did not
fit. That was an artefact of dividing by the fitted one-copy value: `SRS3394713`'s fit sits 35% above
its own median depth at these positions, so a window at 1.93 times its usual depth reads as only 1.4
copies. **Dividing by a fitted parameter and calling the result "copies" is what produced a wrong
mechanism**, and it is why every over-coverage figure in this report is now against each sample's own
measured median instead, with what that is and is not stated.

**What the run still cannot say** is the model's own copy number, which divides by the fitted level
*times a GC multiplier the run does not print*. Recovering it needs the run to record the relative
copy number it used, per record per sample — carried forward for D4, which needs it to explain
ratios by record kind.

## The ratio distribution, and the histogram's ends

Percentiles are the value at sorted index `int(p·n)`:

| | ratio |
|---|---|
| lowest | −104.9288 |
| 25th percentile | −14.5168 |
| median | −10.8659 |
| 75th percentile | −7.9443 |
| 95th percentile | −2.0141 |
| 99th percentile | 10.7679 |
| highest | 34.7583 |

**The four records past the `[-100, 100]` range are all at the negative end; none at the positive
one.** C3's review left this for D2 and D3: the ratio grows with the cohort while the range is
fixed, and the hazard is that real variants and duplications saturate the *same* edge, share one bin
and leave the target-FDR knob nothing to move. At six samples they do not. The negative edge is the
*confidently a real variant* end, where being folded into an end bin changes no verdict. **The
positive end is 65 nats clear**: the largest ratio in the run is 34.8 against an edge at 100.

The growth is real: the largest ratio was **21.6 at one sample** and is **34.8 at six**. C3's review
measured 1,663 on one duplication-shaped record at 63 samples, so a 63-accession run is where the
positive edge comes into play. **D2 settles only that six samples do not reach it.**

## Where the time goes — and what actually binds at scale

```text
time from the start of the calling pass, 988ms in all: calling 934ms (95%), fitting the coverage
models 203µs (0%), scoring 47ms (5%), writing 7ms (1%) — the run's startup before the calling
pass, which reads the reference and the catalogue and opens every input, is not counted here
```

**The shares are of the filter's own window, not of the process.** The clock starts when the spill
is named, just before the calling loop — so reading the reference, building the segmentation,
opening every input and fitting the run's parameters all fall outside both the numerator and the
denominator. What is inside `calling` is the calling loop and the spill's own flush. The whole
process takes **3.50 s** (median of five), so the four passes are 0.99 s of it: **scoring is 5% of
the filtered part and 1.3% of the process**, and the share reproduced at 5–6% in every one of ten
runs made in review.

**So parallel scoring is not worth a plan at this size.** Where it becomes worth one: 2,311 records
× 6 samples is 13,866 sample-records in **47 ms**, or **3.4 µs each**.

| run | sample-records | scoring, one thread |
|---|---|---|
| this one | 13,866 | 47 ms |
| the same slice at 63 accessions | 145,593 | ~0.5 s |
| spec §4's largest — 3,000 samples, 5 million records | 1.5 × 10¹⁰ | **about 14 hours** |

Pass two really is O(records × samples): both grids are fixed constants, so the per-record work is
linear in the sample count. **But that hour count is a lower bound, not an estimate.** At six
samples the precompute's Wright table is 29 kB and sits in cache; at 3,000 samples it is 14.4 MB
streamed per record, and the per-record carrier table is allocated and zeroed each time. The
per-sample constant grows with the cohort, so the real figure is above this one — the opposite of
how "order of magnitude" is usually read.

**And the disk binds first.** The spill held **647 kiB** for this run. Decomposed by measurement: a
record costs **134 bytes of fixed line** plus about **25 bytes a sample** — 14.2 in the VCF line's
own sample column and about 11 in the spill's per-sample row. At 3,000 samples that is **75 kB a
record**, and five million records is **about 350 GB of scratch beside the output**. A run needing
half a day of one core and a third of a terabyte of temporary disk is constrained by the second, and
**spec §5 prices memory and time and does not mention it.**

For scale here: the spill is 6.6 times the *filtered* run's own gzipped output and 7.8 times the
pre-filter run's, the two differing because the filtered output carries `PARALOG_LR` and
`PARALOG_POST` on every record. It is deleted when the run ends, whatever way it ends.

## Memory: a difference smaller than the noise, and I got its sign wrong

D2 asks for peak resident against the filter-off run. **One run of each cannot answer it.** Five
replicates of each arm here, and eight of each in review, every run measured in its own process so
no high-water mark can mask another's:

| | this report, n = 5 | review, n = 8 |
|---|---|---|
| filter off | 418.6–446.1 MB (median 436.7) | 409.9–424.6 (median 417.4) |
| filter on, dropping | 408.6–433.2 MB (median 422.1) | 409.5–443.9 (median 426.1) |

**The two sets disagree on the sign.** Mine had the filter-on runs 15 MB lower; the review's had them
10 MB higher, at U = 11 against a two-tailed 5% critical U of 13 for eight against eight. **A draft
of this report offered a mechanism for the direction it happened to measure** — that with the filter
on the VCF writer is not open during the calling pass — and that was a story fitted to noise. Both
differences are 2–3% of a 420 MB footprint.

**The honest statement: the filter's peak-memory cost at six samples over 2,311 records is below a
run-to-run spread of about 30 MB, and neither its size nor its sign is established here.** Wall is
likewise unchanged, 3.50 s against 3.54 s at the median on runs that vary by 2% among themselves.

## What this step changed in the caller

Two lines were added to the run report:

- **where the run's time went**, per pass with each pass's share, saying explicitly that the shares
  are of the filter's window rather than of the process — the difference is a factor of three and a
  half here;
- **how much disk the parked records took**, which B2's review filed as unreported.

`WhereTheTimeWent` hangs off `FilteredRun`; the calling pass is timed **from the spill's own birth**,
the one place both subcommands go through on the way into pass one. Timing it in each subcommand
would have put the same clock in the two copies of the wiring that C4's review already filed as
duplicated. `SpillFile` gained the two fields behind it: when parking began, and the size read
**after** the flush, so it is the whole file. A failed stat leaves `None` rather than ending the run.

**Two tests.** The timing one pins what a fixture can: that the four durations are four separate
clocks and that the four printed shares sum to about a hundred, which a line dividing by the wrong
total would fail. The spill-size one asserts the size is above zero — the value a read taken *before*
the flush would give — and that the printed figure is the one the field holds.

## Deviations from the plan, recorded

1. **"Calling" is the calling pass plus what surrounds it**, from the spill's creation to the fit
   beginning, and the run's startup is outside every figure. The plan asks for wall per pass; this
   gives the last three exactly and the first as a window that also holds opening the alignments and
   the merge. Nothing in the filter can see the process's own start.
2. **Peak resident is a range over replicates**, not the single figure the harness prints.
3. **Over-coverage is reported against each sample's own median window depth**, not against the
   fitted one-copy level, because the model's own copy number needs a GC multiplier the run does not
   print. The fitted values are still reported — they are the fit's outcome — but nothing divides by
   them and calls the result copies.
4. **The record-kind counts are here; the by-kind analysis is D4's**, which wants HG002 beside them.

## What D2 settles, and what it leaves

- **Spec §8's parallel-scoring question**: 5% of the filtered run at six samples, so no plan needed
  at this size; about 14 hours as a lower bound at 3,000 samples and five million records, which is
  where one belongs.
- **⚠ A finding the plan did not ask for: the spill is the constraint at scale, not the scoring.**
  About 350 GB of scratch at the top of spec §4's range, from a measured 134 bytes a record plus 25
  a sample. Spec §5 prices memory and time and not this.
- **The ±100 histogram range does not bite at six samples**, and the hazard named was about the
  positive edge, which is 65 nats clear. A 63-accession run is where to look again.
- **The "records with no alternative allele still count towards π" question is moot here**: there are
  none, and this cohort's ground produces none.
- **Whether the flagged clusters are real duplications is still inference** — clustered beyond any of
  200,000 random draws, 8.7 to 11.5 times more heterozygous, and over-covered in the samples that
  carry the heterozygote. There is no tomato truth set. **D4's comparison against production's filter
  over these same accessions is the nearest external check this plan has.**
- **The effective cohort is five plants**, and D4 should say so wherever it counts samples.
