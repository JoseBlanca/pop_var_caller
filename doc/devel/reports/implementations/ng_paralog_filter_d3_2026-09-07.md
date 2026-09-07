# The hidden-duplication filter — D3: HG002, three hundred reads a position

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone D, step D3
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §4, and **§3.2, which
this step's measurement contradicts**
**Branch:** `ng-paralog-filter`

## The answer

**The fit accepts at three hundred reads a window. The filter's target is met on substitutions and
missed by seven times on the file as a whole, and the whole of the miss is indels.**

Against GIAB's own HG002 v4.2.1 benchmark, both sides normalised against the same reference:

| what the filter removed | records | of which GIAB says are real |
|---|---|---|
| **substitutions** | 265 | **2 — 0.8 in 100** |
| **indels** | **19** | **17 — 89 in 100** |
| everything removed | 284 | 19 — **6.7 in 100**, against a target of **1 in 100** |

**The removed indels are true heterozygous indels, and the reason is mechanical.** Their
alternative-read fractions cluster on **0.32 to 0.37** — a third, not a half — because reads
carrying an indel are harder to place, so at 300× a real heterozygous indel reads near a third. The
model's carrier configurations include `(T=3, m=1)`, whose expected fraction is exactly one third.
**The filter is reading alignment bias as a collapsed duplication.**

**It bears on spec §3.2**, which rules that an insertion or a deletion carries the allele signal
because "a read's allele is read off directly" there. At this depth it is not: the fraction is
biased by how the read aligns.

**Not acted on, by the owner's ruling (2026-09-07): "we'll study the giab problem and we'll fix it
if necessary. The main objective right now is to have a working filter so we can carry out the real
experiments and decide how to tweak it."** So this is recorded as a measurement to come back to, and
neither §3.2 nor the code changes here. What it costs to leave alone is one number: at high depth
the filter removes about nine real indels for every one it should, and indels are 3% of what it
removes.

The rest of the run:

| | |
|---|---|
| records without the filter | 8,242 |
| records scored | 8,242 — none unscorable |
| records removed at `--paralog-fdr 0.01` | 284 (3.4 in 100) |
| **the fit** | **accepts: one copy at 302.45 reads a window, σ₀ 0.092** |
| fitted duplication rate `π` | 0.038073, fitted from this run; the EM converged |
| the cut | likelihood ratio 4.6500 |
| ratios past the histogram's ±100 ends | 2,768 — **all at the negative end, none at the positive** |

```text
##paralogFilter=target_fdr=0.0100;pi=0.038073;lr_cut=4.6500;em_converged=true;samples_with_coverage_model=1/1
```

## What was run

**HG002, one sample, through the stored-file route the plan names.**
`benchmarks/human_genome_bottle/crams/HG002_reads_selected_1000_rg.cram`, pre-restricted to the
benchmark's 1,000 high-confidence intervals — about 5 Mb of GRCh38 — against the GRCh38 no-alt
analysis set and its repeat catalogue. `generate-psps` walked it once (167 s, peak 776 MB) and
`call-from-psps` read the stored file three times: filter off, dropping, tagging.

**Psp mode is the point, not a detail.** The two subcommands carry their own copy of the filter's
wiring, and C4's review filed that duplication as untested. D1 and D2 ran `call-from-alignments`;
this is the first real-data run of the other copy, and the only earlier exercise of it is a
command-line fixture over a cohort that writes no records.

**Spec §4's high-depth question is answered.** It says the histogram's depth axis scales with the
sample "so the fit is not rejected for range as production's was on a 100× human sample". The fit
accepts and reports what it came to: **one copy at 302.45 reads a window, σ₀ 0.092**, against
tomato's 0.158 to 0.302 across six accessions. A one-copy window's relative depth scatters about a
third as much here, which is what three hundred reads buys.

## The external check, and the failure it found

The benchmark's reference callset is **GIAB's own HG002 benchmark**, not a derived one — its header
records how it was made:

```text
##bcftools_viewCommand=view -R HG002_bench_azar_sorted_1000.bed -Oz -o … HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
```

so it is `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` restricted to the bench intervals and nothing
else. **An earlier draft of this report warned the reader it was "one step removed from the
consortium's release"; that was wrong, and it told a reader to discount the only external check in
this plan for a reason that does not exist.**

**Both sides are normalised before comparing** — `bcftools norm -f <the analysis set> -m -any` on
the truth set and on the filter's output — so an indel written one way matches the same indel
written another, and multiallelic records are split. Without that, exact key matching undercounts
agreement at indels, which is where the whole finding lives.

**The result, again:** substitutions 2 of 265, indels 17 of 19, everything 19 of 284 — **6.7 in
100 against a target of 1 in 100**. On the records it keeps the callset agrees at 97.6% of
substitutions and 97.3% of indels.

**The twenty-two removed indel records, with the fraction of reads carrying the alternative:**

```text
chr11:32724581  AT->A          AD 207,119   0.365      chr4:106325519 TA->T        AD 214,110  0.340
chr21:41354221  CAT->C         AD 203,116   0.364      chr6:160524933 TCT->TC      AD 197,100  0.337
chr22:45432773  CTGG->C        AD 187,109   0.368      chr4:113516633 AT->A        AD 223,114  0.338
chr17:30786697  A->AATCT       AD 159,89    0.359      chr4:129160033 CCATGGTATT->C AD 200,101 0.336
chr13:23222437  TGAGA->T       AD 204,108   0.346      chr15:50131750 AATATT->A    AD 228,113  0.331
chr10:42337684  CATT->C        AD 197,93    0.321      chr6:24080605  AG->A        AD 211,103  0.328
chr10:23553951  AT->A          AD 199,96    0.325      chr3:135003417 GA->G        AD 199,91   0.314
chr1:101051327  GTA->G         AD 233,103   0.307      chr14:39509446 T->TATAC     AD 174,74   0.298
chr12:68514352  CG->GG         AD 191,68    0.263      chr2:172341451 TTC->TTT     AD 164,54   0.248
chr17:37058143  G->GTATTTATTT  AD 182,54    0.229      chr8:67662433  T->TAC       AD 193,57   0.228
chr4:136395577  TC->CC         AD 168,47    0.219      chr9:131452079 AATAAGATTG->A AD 232,46  0.165
```

Sixteen of the twenty-two sit between 0.29 and 0.37. **A real heterozygote gives half.** The gap is
alignment bias, and the model has a carrier configuration sitting exactly in it.

## What it removed, and what the coverage half contributed

**284 records: 262 biallelic SNPs, 15 deletions, 4 insertions, 3 equal-length substitutions — and
no repeat tracts**, though the run wrote 580. **Every one of the 284 is a heterozygote**; not one
homozygous call was flagged.

**The verdicts rest on allele balance, and the coverage half contributes no positive evidence.**

| | alternative-read fraction | window depth ÷ the fitted 302.45 | window GC |
|---|---|---|---|
| the 284 removed | median **0.249** (p10 0.183, p90 0.311) | median 0.92 | median **0.473** |
| the 7,958 kept | median **0.483** among heterozygotes (p10 0.430, p90 0.532) | median 0.99 | median 0.406 |

The fraction is `AD` summed over the alternatives divided by the whole of `AD`, which is what the
scorer reads; a heterozygote here is a kept record whose two called alleles differ, 4,952 of them.
(A slightly different selection moves the ninetieth percentile by about 0.004, so the definition
matters more than the digit.)

**The depth column is not the model's copy number, and this report will not pretend otherwise.**
The scorer divides a window's depth by `single_copy_scale × gc_multiplier(gc)`, and **the run prints
the scale but not the GC curve**, so a record's copy number cannot be recovered from the output.
What can be said from measurement: the removed records' raw depth is 8% below the fitted scale
*and* their GC is 0.067 higher than the kept records', where Illumina depth falls with GC — so the
two corrections push opposite ways and the removed records sit near one copy rather than clearly
above or below it. **An earlier draft divided by the scale alone, called the result copies, and
concluded the removed records were under-covered.** They are not; they are ordinary.

**Either way the load-bearing point holds: at 1.0 copies there is no excess coverage, and a
duplication carrier is expected at 1.5.** The coverage half is voting *against* every one of these
verdicts and losing.

**Why it loses, with the arithmetic rather than a gesture.** For the top-scoring record — 72
alternative reads of 470 — the best-scoring carrier configuration is `(T=4, m=1)`, and against a
real heterozygote's half it gains **111.5 nats** on the alleles while paying **34.4 nats** on the
coverage: a margin of about **3.2 to 1**, not the order of magnitude an earlier draft claimed. The
conclusion is the same and the size is not.

**What that means, and why it is the owner's call.** Production's spec §1 keeps coverage in the
score because it is the signal that separates a collapsed duplication from an introgression: an
introgressed haplotype also gives odd allele balance, and only the duplication raises depth. **On a
single sample at this depth that safeguard is outvoted three to one.** The removed substitutions
are overwhelmingly not GIAB variants, so removing them looks right; what is not established is that
they are duplications rather than mapping artefacts or contamination, which give the same signature
— and the indels show the same machinery removing things that are simply real.

**One thing the earlier draft claimed and the data does not support**: that the removed records'
fractions sit on the model's `m/T` grid. Across the range they occupy, the grid values
`1/8, 1/6, 1/4, 1/3, 1/2` are never more than 0.042 apart, so *every* value lands near one; and the
284 have a median of 0.249, which is a quarter, not the sixth the ten highest happen to show. What
the data does support is the separation from a half, which does not by itself name a cause.

## Repeat tracts: none flagged, and none close

**580 tracts were written and none was removed.** The highest tract ratio in the run is **1.4843**
against a cut of 4.6500, so none came within three nats. **331 of the 580 share one ratio exactly,
−0.6572** — a tract is scored on coverage alone (spec §3.2), its spilled rows carry no read counts,
and where the coverage says one copy the ratio is a constant. The deepest tract window is 428 reads,
1.42 times the fitted scale.

**So "no tract was flagged" is a fact about this run's coverage and not about tracts.** A tract here
cannot respond to allele balance at all, which is the half that removed every one of the 284. D4
counts this beside D2's 21 tracts.

## The histogram's ends — the question C3's review left for D3

**2,768 of the 8,242 ratios ran past the `[-100, 100]` range, and every one is at the negative
end.** None reached the positive one; the largest ratio in the run is **84.16**, still 16 nats short
of the edge.

C3's review named the hazard: the ratio grows with the cohort while the range is fixed, and if real
variants *and* duplications saturate the same edge they share one bin and the target-FDR knob has
nothing to move. **They do not, and the reason is structural.** The negative end is where a record
is confidently a real variant, and a third of this file lands there because at 300 reads an ordinary
heterozygote at half is overwhelming evidence against every carrier configuration. A record already
far past the cut on the safe side does not change its verdict by being folded into an end bin.

| | ratio |
|---|---|
| lowest | −348.37 |
| 25th percentile | −193.84 |
| median | −28.63 |
| 75th percentile | −19.92 |
| 95th percentile | −0.66 |
| 99th percentile | 21.12 |
| highest | 84.16 |

**The positive edge at a large cohort is still untested.** The largest ratio was 21.6 at one tomato
accession, 34.8 at six, 84.2 here at one deep sample; C3's review measured 1,663 on one
duplication-shaped record at 63 samples. **Depth pushes the negative tail out and cohort size pushes
the positive one**, and no run in this plan has been large enough to test the second.

## The recorded cut and the flag disagree, on real data

A3's review recorded that `lr_threshold` is the crossing bin's *centre* while the flag is decided on
the bin, so a record in the lower half of that bin is removed while sitting below the number the
header quotes. **It happened here**: the header says `lr_cut=4.6500` and the lowest-scoring removed
record is at **4.6098**. 283 of the 284 are at or above the quoted cut, no kept record is, and the
bins are 0.1 wide over `[-100, 100]` at 2,000 of them. One record in 8,242.

## Time, memory and disk

```text
time from the start of the calling pass, 6.87s in all: calling 6.83s (99%), fitting the coverage
models 60µs (0%), scoring 33ms (0%), writing 7ms (0%)
```

**Scoring is 33 ms of 6.87 s — half a percent**, against 5% on the six-accession tomato slice. The
direction is what D2's rate predicts: scoring grows with records × samples, and this run has 3.6
times the records but a sixth of the samples. It reinforces D2's answer to spec §8 rather than
adding to it.

**The spill held 1.3 MiB** for 8,242 records at one sample — 165 bytes a record, against tomato's
287 at six, consistent with D2's fixed-line-plus-per-sample decomposition.

**Peak resident was 13.4 GB on all three runs alike**, against 776 MB for the walk. **That is the
calling pass and not the filter** — the filter-off arm peaks the same, 13,428 MB against 13,428 —
so it is outside this plan, but it is a large number for one sample over 5 Mb and nothing else in
the D milestone would have surfaced it.

## The three relations, on real reads and on the other subcommand

All three hold: the filter changed no record it did not flag (10,565 lines compared), dropping is
tagging minus the tagged (7,958), and each of the 284 flagged records is its off-run record with the
verdict spliced in.

The off run's own baseline: 8,242 records, sha256
`0c991c8f43c30d3318049c974a0be3a8ed174463f3ae200aa3d89ea04876b7e6` on everything but
`##commandline`.

## Deviations from the plan, recorded

1. **The harness gained a psp mode** (`NG_MODE=psps`), which walks the alignments once and calls the
   stored files three times. The plan names the HG002 psp as its ground.
2. **The truth comparison is beyond what D3 asks for**, and it is where the step's finding came
   from. D3 wants the fit's verdict, the records dropped and the ten highest ratios; a truth set
   exists for this benchmark and nowhere else in the plan.
3. **Peak resident, the spill's size and the pass shares are reported** though D3 asks for none of
   them; the harness measures them anyway, and the 13.4 GB is worth someone's attention.

## What D3 settles, and what it leaves

- **Spec §4's high-depth question**: the fit accepts at 302 reads a window, σ₀ 0.092.
- **Indels: 17 of the 19 removed are GIAB variants**, because alignment bias puts a true
  heterozygous indel's alternative fraction on the model's `(3,1)` carrier value. **Recorded, not
  acted on** — the owner's ruling is that a working filter comes first and this is studied later.
  The number to come back to is 89 in 100, at 300 reads a position, on one sample.
- **The target is met on substitutions** — 0.8 removed in 100 were real, against 1 in 100 — and
  missed on the file as a whole at 6.7 in 100.
- **⚠ Open: on one high-depth sample the coverage half is outvoted about three to one** and does not
  discriminate. That is the signal spec §1 relies on to tell a duplication from an introgression.
- **The `[-100, 100]` range does not bite from depth**; the positive edge at a large cohort is still
  untested.
- **A tract cannot be flagged on this data**, and cannot respond to allele balance at all. D4.
