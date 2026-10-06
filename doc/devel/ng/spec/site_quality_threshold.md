# ng step 11b — which called sites are written: the site-quality threshold

*Design spec, 2026-10-06, with the code (`--min-site-quality`). Settles the one question
[`calling_quality.md`](calling_quality.md) §12 left to this step: below what site quality a called
locus is not written. Reads on: `calling_quality.md` §3.5 (one quality field, and the defect behind
that rule) and §6 (the artifact correction whose output this reads).*

---

## 1. The decision

**A called locus whose site quality is below `--min-site-quality` is not written. The default is 1;
zero writes every called locus.** The locus is dropped, not written with a filter tag, and the run
report counts how many were left out (`loci left out for a site quality below --min-site-quality`).
`--explain-loci` names the outcome `below_min_site_quality` and prints the quality.

Both calling commands carry it, `call-from-alignments` and `call-from-psps`. It applies to every
locus kind, ordinary loci and repeat tracts alike.

## 2. Which number the threshold reads

**The site quality the record would carry**: the baseline the calling loop computes, after the
artifact correction subtracts its two penalties (`calling_quality.md` §6). It is the `QUAL` the
file would show. `calling_quality.md` §3.5 records why this is a rule and not a convenience: the
production caller's gate once read the baseline while `QUAL` showed the corrected number, and
wrote sites `PASS` with a `QUAL` of 0 — 40 false positives at 30× on GIAB HG002.

The gate sits where the record is built, after the correction and before the window coverage is
taken (`run/callers.rs`, both drivers). A locus below it never reaches the file or the
hidden-duplication filter, which runs on written records.

## 3. Why 1, and not GATK's 30

**The site quality already grows with depth, so the threshold does not need to.** The calls a
threshold should remove are the ones whose artifact penalties outweigh their evidence; the
correction floors those at zero whatever the depth. A threshold that also has to stand in for
depth, as 30 does, removes true calls at low depth where every quality is small.

Measured on the three GIAB samples (HG002, HG003, HG004), each over its own 100 regions, at 5×,
10×, 30× and 300×, against GIAB v4.2.1 with vcfeval. These samples share no ground, so every locus
is called from one sample. False and true calls removed by each threshold:

| depth | below 1 | below 10 | below 30 |
|---|---|---|---|
| 5× | 141 / 15 | 192 / 33 | 212 / 117 |
| 10× | 215 / 17 | 247 / 39 | 275 / 131 |
| 30× | 142 / 2 | 147 / 6 | 148 / 18 |
| 300× | 32 / 3 | 33 / 3 | 33 / 3 |

Raising the cut from 1 to 10 removes about as many true calls as false ones at 5× and 10×; 30
costs more than a hundred true calls at each of those depths, about 6 in 100 of them. The true
calls below 1 are thin: repeat tracts at quality 0.6 and 0.8, two SNPs with 8 reads at a 300×
site, and a heterozygote with 18 variant reads of 151 whose allele-balance penalty is 221.

**At the other end of the range**, 50 tomato accessions at about 3× over tomato1's 80 regions,
there is no truth set, so the calls were compared sample by sample with GATK's joint calls on the
same accessions (vcfeval, an allele in either zygosity counting as agreement). Below 1 removes
26,696 sample genotypes GATK does not call and 382 it does; below 30, 33,236 and 951. In a cohort
the site quality pools every sample, so the threshold costs little there at any of these values.
Part of the tomato sites at quality 0 are indels in repeats whose base qualities the old CRAMs
carry rewritten by BAQ (`samtools calmd -Ar`); GATK does not call them either.

With the threshold at its default, the GIAB runs score (missed / false calls, SNPs and indels
together): 5× 746 / 315 → 761 / 173; 10× 200 / 308 → 217 / 93; 30× 49 / 164 → 51 / 22; 300× 40 /
42 → 43 / 10.

## 4. Dropped, not tagged

**Dropped**, by the owner's decision (2026-10-06): a user who wants every called site sets
`--min-site-quality 0`. The alternative, keeping the record with a `LowQual` filter, scores the
same in a benchmark (hap.py and vcfeval treat a filtered record as not called) and keeps the call
visible; it was not taken. GATK drops below its threshold too.

## 5. Open

- **One threshold for every kind of locus.** Repeat tracts' qualities come from the same site
  quality, and the two true tract calls below 1 on GIAB are genuine low-confidence calls; whether
  tracts want their own value has not been measured.
- **Large cohorts at low depth beyond tomato.** The cohort measurement is one panel of 50
  accessions at about 3×. A cohort of thousands pools more evidence per site, which only lowers
  what a threshold of 1 costs, but it has not been run.
