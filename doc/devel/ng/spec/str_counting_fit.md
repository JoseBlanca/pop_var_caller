# The counting fit as the repeat-tract fit

**Status:** draft started 2026-10-12. Nothing here is built beyond the prototype
(`parameter_estimation/joint/ssr_fit/composite_fit.rs`). This spec is written with the owner, one
requirement at a time; the sections below are the ones settled so far. It follows
[`str_allele_window.md`](str_allele_window.md), which widened both fits to ±9.

---

## 1. What this is

**The repeat-tract fit estimates, for each class of repeat tracts, how reads misreport an allele
and which allele lengths the class's loci carry.** A class is every locus with one motif length and
one reference repeat count. Today it climbs to the best-scoring numbers by an average over each
locus's possible allele frequencies, which costs hours on a cohort of 59 accessions and is several
log-likelihood units a locus wrong at 13 allele lengths
(`reports/implementations/fit_precision_quadrature_average_2026-10-01.md`).

**The counting fit estimates the same numbers by expectation maximisation:** it weighs every
genotype each sample could carry at each locus by how well it explains the sample's reads, counts
reads and alleles with those weights, and repeats until no number moves. On tomato2 (59
accessions, 32 Mb) it takes 70 s at the ±9 window where the climb would take about 22 hours
(`reports/implementations/str_allele_window_2026-10-12.md`).

What it estimates for each class:

- **Stutter**, three numbers: the slip rate (the chance a read reports a length other than its
  allele's), the shorter share (of slipped reads, the fraction that came out shorter) and the
  fall-off (how fast slips of two repeats become rarer than slips of one).
- **The allele length variation:** over every chromosome at every locus of the class, the fraction
  carrying each repeat count from −9 to +9 of the reference length.
- **The between-locus share:** how much of that variation lies between loci rather than among
  samples at one locus. It is the Balding–Nichols F_ST with the class's loci in place of
  populations, `F = 1 / (1 + concentration)`, the concentration being the Dirichlet total the
  maths uses. tomato2's classes sit at F ≈ 0.67–0.71.

Each sample's homozygote excess is an input, from the SNP/indel half of the fit.

---

## 2. Requirements settled with the owner

### 2.1 It reports a standard error for every number it estimates

**Settled 2026-10-12.** The climb computes one for every number, so that a reader can tell a
well-measured number from a poorly measured one. The counting fit must do the same before it
replaces the climb. Today its error columns are empty.

What the errors are for, and what they are not for:

- **The three stutter numbers' errors are written into the parameters file**, as each class's
  own-fit errors (format version 2, `fit_precision.md` §5.2).
- **The allele length variation's and the between-locus share's errors are not in the file today.**
  The climb computes them and reports them only in its log. Whether the file should carry them is
  open (§3).
- **The curves across repeat count do not read them.** A class's own slip rate is blended with its
  motif length's curve by an error computed from how many of its reads the fit says slipped
  (`str_slippage_level_curve.md` §7, `slippage_curve.rs` `FittedCell`). The counting fit already
  has that count: the slip rate times the reads crossing the class's loci.

**The errors must allow for samples at one locus sharing that locus's alleles.** The counting fit's
objective treats each sample at a locus as if it drew its alleles from the class's variation
independently of the other samples there. That costs precision, not correctness of the estimates,
but an error computed as if the samples were independent would come out too small. The natural
way to get errors that hold is to resample whole loci, not samples. Which way exactly — a
bootstrap over loci, or the sandwich formula a composite likelihood takes — is open (§3).

---

## 3. Open

1. **How the errors are computed** (§2.1): resampling loci, or the composite likelihood's sandwich
   formula. Either must be checked against drawn classes, where the true spread of the estimates
   is known.
2. **Whether the parameters file carries errors for the allele length variation and the
   between-locus share**, as it does for the stutter numbers.
