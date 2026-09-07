# ng — the hidden-duplication filter: dropping calls that a collapsed paralog explains better than a variant

*Status: draft, 2026-09-06 — the design settled with the owner in the conversation preparing this
port. No code yet. Consumes [`window_coverage.md`](window_coverage.md), which produces the numbers
this scores. The statistics are production's, ported unchanged: the model is
[`doc/devel/specs/hidden_paralog_filter.md`](../../specs/hidden_paralog_filter.md), reformulated
for one sample up by
[`hidden_paralog_single_sample_scoring.md`](../../architecture/hidden_paralog_single_sample_scoring.md);
this document does not re-derive either and cites them where the model is meant. Implementation
plan: [`../impl_plan/hidden_paralog_filter.md`](../impl_plan/hidden_paralog_filter.md). No
architecture document: §3's type blocks are the code shape the plan cites. This is
[`ng_proposal.md`](ng_proposal.md)'s step 11a.*

---

## 1. What this is

**A filter that removes calls better explained by two reference-collapsed gene copies piling their
reads onto one position than by a real variant.** A hidden duplication shows two footprints — a
sample's read depth around the locus is about twice its one-copy level, and the fixed difference
between the two copies reads out as a heterozygote — and the strongest evidence is the two together:
*the heterozygous samples are the over-covered ones.* Production's spec §1–§3 says why coverage is
the one signal that does not also flag an introgression, and this document takes that as settled.

Per written record, the filter computes a likelihood ratio between two stories over every sample's
coverage and allele counts, estimates from the whole run how common paralogs are, turns each ratio
into a probability, and drops — or, on request, tags — the records above a false-discovery cut the
operator sets. Default target: about 1 record in 100 of those removed was a real variant.

**What ng changes about production's filter is where its numbers come from and when they are
known, not what it computes.** Production fits each sample's coverage model at reader-open from a
histogram the pileup stored, scores every record inline, and spills once. ng has no stored
histogram: [`window_coverage.md`](window_coverage.md) accumulates it during the calling pass, so the
model exists only when that pass ends, and the filter becomes three passes over the called records
instead of two (§2). In exchange it runs identically in both modes.

### 1.1 Goals

1. **Production's statistic, bit for bit** — the same per-locus likelihood ratio on the same
   inputs, so the tomato2 validation and the single-sample reformulation carry over.
2. **Both modes, one path.** Direct mode and psp mode run the same three passes on the same
   numbers.
3. **One sample to several thousand, three reads to three hundred**, each end with a stated
   behaviour (§4).
4. **Drop by default; tag on request.** A dropped record leaves its parameters in the header; a
   tagged one carries its score.
5. **Every written record is scored** — everything on coverage and allele balance except a
   repeat tract, which is scored on coverage alone (§3.2). Repeat tracts are first-class.
6. **Memory flat in the record count and in the sample count**: one record in hand at a time,
   plus per-sample constants and one number per record.
7. **Off means off**: with the filter disabled the run writes byte for byte what it writes today.

### 1.2 Non-goals, and what this does not do

- **No new statistics.** The likelihood ratio, the prior's estimation, the FDR curve are
  copied from `src/paralog/`. The one thing that is not a straight copy — scoring a repeat tract
  on coverage alone — is a *configuration* of the copied scorer, shown in §3.2 to be the same
  likelihood with the allele counts marginalised out; nothing new is derived.
- **No mapping-quality term in the score.** `MQDIFF` is already an INFO field
  ([`vcf_output.md`](vcf_output.md) §6) and stays a downstream signal, for production's reason:
  it is shared with introgression.
- **No per-sample inbreeding derivation.** ng's parameters file already carries a fitted
  coefficient per sample ([`parameters_file.md`](parameters_file.md) §3.5); the filter reads it.
- **No histogram in the psp** — [`window_coverage.md`](window_coverage.md) §8.
- **It does not decide which loci exist.** The keep rule and everything before a record is
  built are the merge's ([`cohort_merge.md`](cohort_merge.md) §4.3).

### 1.3 Vocabulary

The model's terms — relative copy number, σ₀, the SFS prior, `F`, the empirical-Bayes prior π,
FDR and q-value — are defined in production's spec §0 and are used here with those meanings.
Three that this document adds:

- **The spill** — an ephemeral file, written once and read twice, holding every record the run
  would have written together with what the filter needs to score it. It exists because a
  record's verdict depends on every other record's score (§2).
- **Pass one, two, three** — the calling pass that writes the spill; the scoring pass that reads
  it, computes every likelihood ratio and calibrates; the writing pass that reads it again and
  applies the verdict.
- **The coverage-only score** — the likelihood ratio for a record whose allele counts are not
  handed to the scorer (§3.2).
- **A generic locus** and **a repeat tract** — the project's existing pair, and the one that
  decides which signals a record is scored on. A generic locus is anything that is not a tract: a
  SNP of any allele count, an insertion, a deletion. Its reads report their allele directly, so it
  is scored on coverage *and* allele balance. A tract's reads are moved between length alleles by
  slippage, so it is scored on coverage alone (§3.2).

---

## 2. Where it sits — three passes, and why not two

**The verdict on the first record depends on the last record's score.** The prior π is fitted by
an EM over every record's likelihood ratio, and the FDR cut is a quantile of the resulting
posteriors, so no record can be dropped or tagged until all have been scored. Production's module
doc says it in one line: the filter "cannot fit in the current single streaming pass"
([`paralog_filter/mod.rs`](../../../../src/var_calling/paralog_filter/mod.rs)). Something must hold the
records between scoring and verdict; production holds them in a spill file, and so does ng (§3.4).

**Production scores inline and needs one spill read; ng cannot, and needs two.** Production's
per-sample coverage model is fitted at reader-open from the psp's stored histogram, so a worker
scores each record as it is called and the spill carries the ratio
([`write_pass.rs:1-19`](../../../../src/var_calling/paralog_filter/write_pass.rs)). In ng the
histogram is complete only when the calling pass ends ([`window_coverage.md`](window_coverage.md)
§1), so:

| pass | reads | does | holds |
|---|---|---|---|
| **one** — calling | the sources, as today | writes each record and its per-sample window pairs to the spill; the cache accumulates the histograms | the cache's per-sample state |
| *between* | the histograms | fits each sample's coverage model | one model per sample |
| **two** — scoring | the spill | scores every record; folds each ratio into the LR histogram; estimates π; builds the FDR curve; resolves the cut | one `f64` per record |
| **three** — writing | the spill | applies the verdict to each record and writes the survivors | one record |

The extra pass is one sequential read of a file about the size of the VCF, which the owner priced
as small beside the psps. The trailer histogram would buy it back and is deferred for the reason
[`window_coverage.md`](window_coverage.md) §8 gives.

**Both modes reach the filter through one function.** Direct mode and psp mode already share
[`call_cohort_from_sources_handing_each_record_over`, `callers.rs:799`](../../../../src/ng/run/callers.rs),
which hands each assembled record to a sink in genome order. With the filter on, that sink is the
spill writer; with it off, the VCF writer as today. The histograms come out of the cache with the
sources when that function returns.

---

## 3. The design

### 3.1 The per-sample coverage model — production's fit, copied

**What one copy's depth looks like in this sample, as a function of GC** — the one-copy depth
level, a GC-bias multiplier, and σ₀, the scatter of a one-copy window's relative depth. Fitted per
sample from its histogram by [`SingleCopyCoverageModel::fit`, `coverage_model.rs:227`](../../../../src/paralog/coverage_model.rs),
copied into ng with its configuration and its tests. Its defaults — the single-copy band
`0.4..1.6`, the minimum bin count 50, the smoothing window, the mode/median guard `[0.5, 1.5]`,
the overflow guard at a fifth — are **inherited from the tomato2 prototype and production, not
re-measured**.

**A rejected fit is not an error.** The fit refuses a sample it cannot anchor — no windows, a
peak in the bottom bin, a mode/median ratio outside the guard, too many windows past the
histogram's range — and production carries such a sample as *absent*: it contributes no coverage
evidence to any record ([`prepass.rs`](../../../../src/var_calling/paralog_filter/prepass.rs)). ng
does the same, and **the run report names each rejected sample and the reason**, so a run over a
cohort whose fits all failed says so rather than silently filtering nothing.

**When**: once, between pass one and pass two, from the histograms
[`window_coverage.md`](window_coverage.md) §3.5 hands over.

### 3.2 What is scored, and with what — decision: every record, repeat tracts on coverage alone

Per record, per sample, the scorer takes four numbers
([`SampleObservation`, `locus_score.rs:43`](../../../../src/paralog/locus_score.rs)): the relative
copy number, the ALT reads, the total reads, and the sample's inbreeding coefficient.

- **Relative copy number** = the sample's window mean depth at the locus, divided by the model's
  expected one-copy depth at the window's GC
  ([`relative_copy_number`, `coverage_model.rs:302`](../../../../src/paralog/coverage_model.rs)).
  The pair comes from the spill; **a sample with an absent pair, or no model, is absent** from
  the record's score.

  **One depth for the whole locus, however many positions it occupies — and the stored number
  already is one.** What each record carries per sample is the *mean* read depth of the window
  centred on its first base ([`window_coverage.md`](window_coverage.md) §3.5), so it is a
  per-position figure already and nothing has to be divided by the span. The window is 500 bases
  wide, so for a SNP, an indel or an ordinary tract the choice of base inside the span moves the
  window by a fraction of a percent of its width. It is a real approximation for a long deletion,
  where the choice of end moves the window appreciably and where the deletion itself depresses
  depth across its own span in exactly the samples carrying it — which is not a duplication
  signal. That is why D4 counts deletions apart from tracts.
- **Inbreeding coefficient** — ng's own, per sample: the parameters file's `[inbreeding]` row,
  fitted per sample from runs of homozygosity by the pre-pass and exposed as
  [`RunParameters::inbreeding_coefficient_by_sample`, `run_parameters.rs:695`](../../../../src/ng/calling/run_parameters.rs).
  Production replaced its per-individual estimate with one cohort number because the
  heterozygosity proxy it had was divergence-contaminated (single-sample doc §5); ng has the
  quantity production wanted, and the copied scorer already takes it as a per-sample slice.
- **ALT and total reads — this is where the record's kind enters, and it is the one decision
  here that goes beyond production.** Production scores only **biallelic SNPs** — one reference
  base, one single-base alternative — and never flags anything else
  ([`is_biallelic_snp`, `calibrate.rs:135`](../../../../src/var_calling/paralog_filter/calibrate.rs)),
  because its allele term models a SNP: a real heterozygote sits near half, a collapsed copy's
  difference at `m/T`. The owner's ruling is that repeat tracts are first-class and the filter
  must apply to them.

  **The line is drawn by slippage, not by the locus's span and not by its alleles.** The
  scorer's allele half compares two stories about *which allele a read carries*: under a real
  heterozygote the reads split near half, under a collapsed duplication at some whole `m/T`.
  Both assume a read's allele is observed correctly. **At a repeat tract it is not.** Slippage
  during amplification makes a read report a tract one or two units from the molecule it came
  from, so reads belonging to one length allele appear at its neighbours, and the observed split
  is smeared by an amount that has nothing to do with copy number. Feeding that to the scorer
  would have it read the smear as evidence about copies. Pooling to *reference length against
  everything else* does not help, because slippage crosses that boundary too. Nowhere else does
  this happen: at a SNP, an insertion, a deletion or a multiallelic site a read's allele is read
  off directly. So:

  - **Every record that is not a repeat tract** — a biallelic SNP, a multiallelic SNP, an
    insertion, a deletion: `alt_reads =` the sum of every alternative's reads,
    `total_reads = AD[0] + alt_reads`, from the record's sample column
    ([`SampleReadCounts::allele_reads`, `vcf/mod.rs:597`](../../../../src/ng/vcf/mod.rs)).
    A sample with `total_reads == 0` is absent from the score, as in production.
    **Alternatives are pooled**, which the H2 story permits unchanged: it asks whether the
    non-reference share sits at some whole `m/T`, and two copies of three carrying two *different*
    non-reference bases give the same two-thirds as two copies carrying one.
  - **A repeat tract**: `alt_reads = 0`, `total_reads = 0`, **and the sample is not skipped for
    it.** The scorer's allele term is
    `alt·ln(vaf) + (total − alt)·ln(1 − vaf)` in both hypotheses
    ([`locus_score.rs:333-345`](../../../../src/paralog/locus_score.rs) for H1,
    [`:386-399`](../../../../src/paralog/locus_score.rs) for H2); at zero reads it is zero under
    every genotype and every carrier configuration, so the genotype and carrier sums collapse to
    one and the ratio rests on the coverage term alone. That is the same likelihood with the
    allele counts marginalised out — no new model, less power. The homozygous-alternative veto
    needs five reads and is inert here, which is right.

  **⚠ The zeros are an abstention, not a measurement.** A sample at a repeat tract has reads and
  its call rests on them; what it does not have is an allele split this scorer's model fits. The
  `0, 0` is the filter declining to use that evidence, and it must never be confused with a sample
  that observed nothing — which is exactly why §3.7 gives the two kinds of record different types
  rather than one type with fields that sometimes mean nothing (trap 1).

  **This is a deferral, not a judgement that a tract's alleles say nothing.** Two copies of a
  tract drifted to different lengths give a bimodal length distribution, which is real evidence of
  a collapsed duplication and arguably better than a SNP's. What is missing is a model of it under
  slippage — §8's tract-aware allele term — and until that exists, feeding tract counts to the SNP
  model would be applying a model known to be wrong there rather than using more evidence
  (the owner, 2026-09-07).

  **The record already carries the flag that decides this.** Every record says whether it is a
  repeat tract, because the writer needs that to order the file (a tract may share a position with
  the generic locus owning its anchor base). Nothing new is computed or stored, and no test reads
  the record's alleles — which is what strikes trap 2.

  **Why not production's rule.** A collapsed duplication leaves the same coverage footprint at a
  tract as at a SNP, and a filter that never looked at tracts would leave every such tract in the
  file while dropping its SNP neighbours. **What is unmeasured**: how many wide-locus records the
  coverage-only score flags on real data, and whether their depth is biased — and the two wide
  kinds are biased *differently*, so they are counted separately (step D4). A tract's observation
  depth runs below its read depth ([`window_coverage.md`](window_coverage.md) §3.1); a deletion
  depresses depth across its own span in the samples that carry it, by construction. **Confirmed
  by the owner, 2026-09-07** — including this span-based line, which replaces an earlier
  biallelic-SNP-versus-everything-else draft.

- **A record with no alternative allele is not scored.** A refused tract written `ALT .` with
  every sample no-called ([`vcf_output.md`](vcf_output.md) §8) establishes nothing to flag.

**The score itself** is [`score_locus_for_paralogy`, `locus_score.rs:248`](../../../../src/paralog/locus_score.rs)
with its per-pass [`ParalogScorePrecompute`](../../../../src/paralog/locus_score.rs), copied. Its
constants — the carrier copy numbers `{3, 4, 6, 8}`, the winsor cap at four copies, the grids,
the pseudocount — are production's ([`paralog/mod.rs`](../../../../src/paralog/mod.rs)),
**prototype-tuned on tomato2 and inherited**. A record with no usable sample is unscored, which is
*kept, never flagged*: its ratio is `NaN`, which the histogram refuses to fold and the verdict
never flags ([`ParalogCalibration::flags`, `calibrate.rs:97`](../../../../src/var_calling/paralog_filter/calibrate.rs)).

### 3.3 From scores to a cut — production's prior and FDR, copied

Pass two folds every finite ratio into a fixed-size histogram of ratios
([`ParalogLrHistogram`, `prior.rs:108`](../../../../src/paralog/prior.rs)), estimates π from it by
EM ([`ParalogPrior::estimate`, `prior.rs:181`](../../../../src/paralog/prior.rs)), builds the
tail-FDR curve over the same bins ([`ParalogFdrCurve`, `prior.rs:230`](../../../../src/paralog/prior.rs)),
and resolves the least-stringent ratio at which the FDR meets the target
([`lr_threshold_for_fdr`, `prior.rs:274`](../../../../src/paralog/prior.rs)). All copied. When the
EM does not converge the prior falls back to `0.03`
([`DEFAULT_FALLBACK_PARALOG_PRIOR`, `calibrate.rs:30`](../../../../src/var_calling/paralog_filter/calibrate.rs)),
the run warns, and the header says so.

**Pass two also keeps every record's ratio in memory** — one `f64` per record, in spill order —
so pass three applies the cut to the value the histogram was built from rather than re-scoring.
That is production's own argument for storing the ratio: the cut then matches the histogram *by
construction*. Forty megabytes at five million records, flat in the sample count.

### 3.4 The spill — a framed binary file of finished lines

**Pass one writes, for each record the run would have written, a framed entry holding the
record's VCF line as the encoder produced it and what the scorer needs beside it.**

```text
entry :=
  contig          varint     -- the three fields the writer's ordering check reads
  position        varint     -- (`place_of`, writer.rs)
  is_repeat_tract u8         -- and it also selects the per-sample shape below (§3.2)
  line            bytes      -- length-prefixed; the record's line, no newline
  samples         varint     -- the run's sample count, dense
  per sample, when is_repeat_tract is 0:
    gc_fraction   u32        -- an f32's bits; NaN = absent
    mean_depth    u32        -- an f32's bits
    ref_reads     varint     -- AD[0]
    alt_reads     varint     -- every alternative's reads, summed (§3.2)
  per sample, when it is 1:
    gc_fraction   u32        -- an f32's bits; NaN = absent
    mean_depth    u32        -- an f32's bits
```

**One flag, doing both jobs, because they are the same question.** The tract flag is already in
the entry for the writer's ordering rule — a tract may share a position with the generic locus
owning its anchor base — and §3.2's rule is *is this a repeat tract*. Storing a second flag beside
it would be one more byte a record and one more way for two fields to disagree, so the row shape
is read off this one.

**A tract's row carries no read counts at all, and that is the point.** Writing `0, 0` there would
put the filter's abstention in the same two fields that elsewhere hold a measurement, where nothing
downstream can tell it from a sample that observed nothing — and production's scorer drops a sample
at `total == 0`, so the confusion silently empties every tract's score (§6 trap 1). Two shapes make
the mistake unspellable, and they make the file smaller: a tract costs eight bytes a sample instead
of ten or more.

**Why the line and not the record.** The record type has ten fields behind a constructor that
asserts every parity between them ([`VcfRecord::new`, `vcf/mod.rs:187`](../../../../src/ng/vcf/mod.rs));
a codec over it is ten field codecs and a reconstruction that must satisfy every assertion, and a
field the record gains later is a field the codec silently drops. The line is the record's own
bytes: pass three writes them back unchanged except for the two columns the verdict touches, so
byte-identity for every unflagged record holds by construction, and nothing about the record's
shape has to be mirrored. The per-sample side fields are the four numbers the scorer reads,
written so the scorer never parses text — and at a repeat tract there are only two of them.

**The codec is hand-rolled over the psp's varint primitives**
([`encode_u64_leb128`, `psp/varint.rs:46`](../../../../src/psp/varint.rs)), with **the encoder
destructuring its input exhaustively** — the codebase's idiom for a struct that must not gain a
field silently ([`cohort_merge/mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs)'s `render`;
[`types.rs:233`](../../../../src/var_calling/types.rs)) — and **floats written by bit pattern**, so
an absent pair survives. **No manifest, no layout check, no version**: the file is created, read
twice and deleted inside one run, and its only reader is the process that wrote it. That is the
owner's ruling, and the distinction from the psp — a user-facing artifact that has to interpret
itself — is the reason.

**Where it lives**: `<output>.paralog-spill.tmp`, beside the output on the same filesystem, in the
convention the writer already uses for `<output>.tmp`. **Deleted when the run ends, success or
failure**, by a guard that runs on every exit path the process reaches; a run killed from outside
leaves it, as it leaves `<output>.tmp`. Streaming compression is not built: the file is about the
VCF's size, and the owner ruled it is not worth building unless it is easy — it is a later change
to one writer and one reader if it becomes so.

**Alternatives that lost.** *Serde with a binary codec crate*: it would work, and what it buys —
a field cannot be forgotten — the destructure already gives without a new dependency in a manifest
where every entry carries a justification. *The VCF as its own spill*: it rewrites the file anyway,
so the I/O is the same, and it makes an uncalibrated VCF a thing that exists on disk and can be
consumed by mistake. *Holding records in memory*: the entry is per sample per record, so it is the
cohort size times the record count.

### 3.5 Pass three — the verdict, and what the file says

Per entry, in spill order, with its ratio from pass two:

- **flagged** ⇔ the ratio is finite and its tail FDR is at or below the target
  ([`flags`, `calibrate.rs:97`](../../../../src/var_calling/paralog_filter/calibrate.rs));
- **dropped** if flagged and `--paralog-filter-tag` is not given: the line is not written and
  the drop is counted;
- **tagged** if flagged and the flag is given: the `FILTER` column becomes `hiddenParalog`, or
  `<existing>;hiddenParalog` where the record already carried a filter, and the line is written;
- **written unchanged** otherwise — except that every scored record, flagged or not, gains two
  INFO fields.

**The two INFO fields**, on every record with a finite ratio, absent on an unscored one:

```text
##INFO=<ID=PARALOG_LR,Number=1,Type=Float,Description="Log likelihood ratio of a hidden paralog over a real variant, from coverage and allele balance across samples">
##INFO=<ID=PARALOG_POST,Number=1,Type=Float,Description="Posterior probability that this locus is a hidden paralog, under the run's fitted paralog rate">
##FILTER=<ID=hiddenParalog,Description="Better explained by a reference-collapsed duplication than by a variant; tail FDR at or below the run's target">
```

`PARALOG_POST` is the id [`vcf_output.md`](vcf_output.md) §6 reserved for this document. The
ratio is written so that a run that drops can be audited against one that tags.

**The header carries the run's calibration**, production's line with one addition:

```text
##paralogFilter=target_fdr=0.0100;pi=0.031250;lr_cut=4.2117;em_converged=true;samples_with_coverage_model=61/63
```

Production writes the first four ([`paralog_provenance`, `write_pass.rs:120`](../../../../src/var_calling/paralog_filter/write_pass.rs))
because a dropped record leaves no per-record trace; the fifth says how much of the cohort the
coverage evidence rested on. ng's header writes its meta lines as `##key=value`
([`header_text`, `header.rs:431`](../../../../src/ng/vcf/header.rs)); this is one more.

**The writer gains one entry point**: write a line whose contig, position and tract-ness are
given, running the same ordering check `write_record` runs
([`writer.rs:114`](../../../../src/ng/vcf/writer.rs)). Patching the line touches columns seven and
eight only — `FILTER` and `INFO`, tab-separated, both always present in ng's records — and a
record whose `INFO` is `.` has it replaced rather than appended to.

**The run report** gains: records dropped as hidden paralogs, or tagged; samples whose coverage
model was rejected, each with its reason; π, the cut, and whether the EM converged.

### 3.6 The command line

| flag | default | meaning |
|---|---|---|
| `--paralog-fdr <f>` | `0.01` | the target false-discovery rate among the records removed; `0` turns the filter off |
| `--paralog-filter-tag` | off | keep flagged records, on the `hiddenParalog` filter, instead of dropping them |

Both on `call-from-psps` and `call-from-alignments`
([`CallFromPspsArgs`, `call_from_psps.rs:75`](../../../../src/pop_var_caller_exp/call_from_psps.rs)),
under an `Advanced` heading beside the merge's knobs. The window and histogram parameters
([`window_coverage.md`](window_coverage.md) §3.6) are not exposed in this document's scope; they are
compiled defaults until a measurement says one needs to move.

**Off means the spill is never created and the sink is the VCF writer**, so a run with
`--paralog-fdr 0` is the run before this work, byte for byte.

### 3.7 The types

```rust
/// One sample's four numbers at one record, or absent.
/// Production's `SampleObservation`, copied — the scorer's input.
pub struct SampleObservation {
    pub relative_copy_number: f64,
    pub alt_reads: u32,
    pub total_reads: u32,
    pub inbreeding_coefficient: f64,
}

/// What pass two needs from pass one, per sample, for the whole run.
pub struct ParalogScoringContext {
    /// One fitted model per sample, `None` where the fit was rejected (§3.1).
    coverage_models: Vec<Option<SingleCopyCoverageModel>>,
    /// The parameters file's coefficient per sample (§3.2).
    inbreeding: Vec<f64>,
    /// σ₀ per sample, `NaN` where absent — the scorer's parallel slice.
    single_copy_depth_sd: Vec<f64>,
    /// Production's per-pass tables, built once.
    precompute: ParalogScorePrecompute,
}

/// One spill entry, as §3.4 lays it out.
pub struct SpillEntry {
    pub contig: ContigId,
    pub position: Position,
    pub is_repeat_tract: bool,
    pub line: Vec<u8>,
    /// **Which signals this record's samples carry** — the variant agrees with
    /// `is_repeat_tract` above, and the writer refuses an entry where it does not (§3.2).
    /// A scorer cannot ask a tract for an allele split slippage has smeared.
    pub samples: SpilledSamples,
}

pub enum SpilledSamples {
    /// A generic locus: a SNP of any allele count, an insertion, a deletion. Both signals;
    /// a sample at zero total reads is absent from the score, as in production.
    GenericLocus(Vec<GenericLocusSample>),
    /// A repeat tract. Coverage alone; **no sample is skipped**, because there is no read
    /// count that could be zero.
    RepeatTract(Vec<RepeatTractSample>),
}

pub struct GenericLocusSample {
    pub window: WindowCoverage,   // the pair, NaN = absent
    pub ref_reads: u32,
    pub alt_reads: u32,           // every alternative's reads, summed
}

pub struct RepeatTractSample {
    pub window: WindowCoverage,   // the pair, NaN = absent
}

/// What pass two settles for pass three. Production's `ParalogCalibration`, copied,
/// plus the per-record ratios in spill order.
pub struct ParalogVerdicts {
    pub calibration: ParalogCalibration,
    pub ratios: Vec<f64>,
}
```

---

## 4. One sample and three thousand, three reads and three hundred

- **One sample.** No gate and no branch: the single-sample reformulation removed the
  minimum-sample count, and the score self-gates — when coverage cannot tell one copy from two,
  both stories fit the data equally and the ratio sits near zero, which never crosses a cut
  (single-sample doc §3). The SFS prior's grid degenerates to one point at `N = 1`
  (`[1/2N, 1 − 1/2N]` is `[½, ½]`); **the plan runs one sample explicitly** (step D1) because the
  copied precompute has not been exercised there. Calibration is soft — the ratio distribution is
  less bimodal, so π is noisier and the fallback fires more often — and the header says which
  happened.
- **Three thousand samples.** Per record: 3,000 spill sample entries (about 20 bytes each) and a
  score over 3,000 samples, on one thread in pass two. The precompute tables are
  `grid points × samples`. Nothing is held across records but the ratios. **Pass two's wall at
  that size is unmeasured**; the plan measures it at 63 and reports (step D2), and parallelising
  the scoring is deferred until that number exists (§7).
- **Three reads a position.** σ₀ is large, the score self-gates, and the filter removes little;
  the tomato slice at this depth is the plan's first real run. A heterozygote's allele balance
  from three reads is weak evidence and the model already weighs it as such.
- **Three hundred reads a position.** The histogram's depth axis scales with the sample
  ([`window_coverage.md`](window_coverage.md) §3.4), so the fit is not rejected for range as
  production's was on a 100× human sample. HG002 is the plan's check that the fit accepts and the
  verdicts are sane (step D3).

## 5. Cross-cutting concerns

- **Memory** — one spill entry in hand at a time in passes one and three; the ratios (8 bytes a
  record) and the per-sample context through pass two. The histograms are
  [`window_coverage.md`](window_coverage.md)'s and are dropped once the models are fitted.
- **Time** — two sequential reads of a VCF-sized file, and one score per record. Unmeasured; the
  plan reports pass shares on the tomato slice.
- **Errors** — a spill write or read failure is a `RunError` naming the file; a rejected fit is
  a reported fact, not an error; EM non-convergence is a warning with a fallback. A run that
  fails after pass one deletes its spill on the way out.
- **Concurrency** — none in passes two and three; they are sequential by design and the output
  is a function of the spill's order, which is genome order. Scoring in parallel is deferred (§7).
- **Determinism** — the standing oracle: byte-identical output at any thread count, with the
  filter on and off.

## 6. Traps — what will bite the coder

1. **Do not skip a zero-read sample at a repeat tract** — and §3.7's types are what stop you.
   Production's `build_observation` returns `None` at `total == 0`
   ([`calibrate.rs:173`](../../../../src/var_calling/paralog_filter/calibrate.rs)), which is right
   when the sample observed nothing and wrong when the *filter* chose not to read its alleles. A
   tract's samples have reads and their calls rest on them; what they lack is a split the scorer's
   model fits. Storing that abstention as `0, 0` in the fields a measurement lives in is what makes
   the two indistinguishable, so `SpilledSamples::RepeatTract` carries no read counts and the skip
   cannot reach it. **The failure this prevents is silent**: every tract scores on nothing and is
   never flagged, and the file looks correct.
2. **~~The SNP test reads the record's alleles before padding.~~** *Struck, 2026-09-07.* Which
   signals a record carries is decided by `is_repeat_tract`, which the entry already holds for the
   writer's ordering rule (§3.2), so no test looks at an allele at all and there is nothing to get
   wrong about `VcfRecord::alleles()` being the span's own sequences while `REF`/`ALT` gain an
   anchor base at encode time.
3. **`F` per sample, and the slice's length is the cohort's.** The scorer returns a neutral
   score on a length mismatch rather than failing ([`locus_score.rs:248`](../../../../src/paralog/locus_score.rs)).
   Build the slice from the parameters file in the run's sample order and assert the length.
4. **`NaN` is *unscored*, and it must stay `NaN` through the spill, the ratio vector and the
   verdict.** Any step that turns it into zero makes a record with no evidence look neutral and
   fold it into π.
5. **Byte-identity is by construction only for the columns the verdict does not touch.** The
   oracle is: filter on, at a target no record reaches, then strip the two INFO keys — the result
   must equal the filter-off file exactly (§9).
6. **The writer's ordering check must still run in pass three**, from the three head fields the
   entry carries; bypassing `write_record` must not bypass its check.
7. **A pre-existing `FILTER` value is joined, not replaced.** A record on `EMNoConv` that is
   also flagged carries both, `;`-separated, which the VCF grammar allows and ng's header
   declares.

## 7. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the model constants and grids | [`src/paralog/mod.rs`](../../../../src/paralog/mod.rs) | copied into `src/ng/paralog/`, production untouched |
| the coverage model fit | [`coverage_model.rs`](../../../../src/paralog/coverage_model.rs) with its tests | copied; reads ng's histogram type, which is production's shape |
| the per-locus score and its precompute | [`locus_score.rs`](../../../../src/paralog/locus_score.rs) with its tests | copied; **the parity oracle is bit-identical ratios against production's on identical synthetic inputs** |
| the prior, the EM, the FDR curve | [`prior.rs`](../../../../src/paralog/prior.rs) with its tests | copied |
| the calibration and the verdict | [`ParalogCalibration`, `flags`, `posterior`, `calibrate.rs:64-118`](../../../../src/var_calling/paralog_filter/calibrate.rs) | copied |
| the header line | [`paralog_provenance`, `write_pass.rs:120`](../../../../src/var_calling/paralog_filter/write_pass.rs) | copied, one field added |
| the spill's framing idea | [`spill.rs`](../../../../src/var_calling/paralog_filter/spill.rs) | the shape (length-framed, written once, read twice, deleted) is kept; the payload is §3.4's, not production's record codec |
| the varint codec | [`psp/varint.rs`](../../../../src/psp/varint.rs) | called as-is |
| the sink both modes share | [`call_cohort_from_sources_handing_each_record_over`, `callers.rs:799`](../../../../src/ng/run/callers.rs) | the filter is a sink; the function does not change |
| the inbreeding coefficients | [`RunParameters::inbreeding_coefficient_by_sample`, `run_parameters.rs:695`](../../../../src/ng/calling/run_parameters.rs) | read as-is |
| the reference draft of the maths | `benchmarks/tomato2/src/build_paralog_lr.py`, `build_paralog_eb.py`, `build_gc_normalization.py` | what production was validated against; not an oracle for ng |

**What is not an oracle**: production's filter on the same slice. Its depth is raw depth where
ng's is observation depth, its `F` is one cohort number where ng's is per sample, and it scores
SNPs only. The plan runs both on the tomato slice as a *sanity comparison* — the flagged sets
should mostly agree — and reports the difference rather than chasing it.

## 8. Deferred, with a recommended home

- **The coverage-only score's power at tracts and indels**, and whether the tract depth caveat
  biases it. Home: the plan's D4 measurement first; if it flags nothing or flags wildly, the
  decision in §3.2 comes back here.
- **A tract-aware allele term** — an H2 story for a collapsed copy whose difference is a repeat
  length. Home: a later document under this one, once the coverage-only measurement says whether
  it is needed.
- **Scoring pass two in parallel.** Home: the plan's D2 measurement decides whether it is worth
  a step.
- **A mapping-quality rescue of the diverging-paralog blind spot** — production's spec §7 names
  it and keeps it out of the score. Home: a separate filter, if ever.
- **Streaming compression of the spill** — one writer and one reader, if it proves easy.

## 9. Resolved decisions & open questions

- **Three passes, not two — resolved** (§2). The trailer alternative is
  [`window_coverage.md`](window_coverage.md) §8's.
- **A spill, not memory and not the VCF — resolved** (§3.4).
- **Lines in the spill, not records — resolved** (§3.4).
- **Drop by default, tag on a flag — resolved** (the owner, 2026-09-06). Production drops with
  no escape hatch.
- **`F` per sample from the parameters file — resolved** (§3.2).
- **Every record scored — resolved** (the owner, 2026-09-07). It is the same likelihood with the
  allele counts marginalised out where they are not read, and the owner's ruling is that tracts
  are first-class. Step D4 still counts what it flags by kind and can send it back.
- **Where the line falls — resolved, after two wrong drafts** (the owner, 2026-09-07). The first
  draft split biallelic SNPs from everything else, which was the shape of production's test rather
  than a reason. The second split one genomic position from more than one, and the owner refuted
  it: **once the depth of a multi-position locus is taken as one number for the whole span, span
  distinguishes nothing** — and a deletion's reads report their allele as directly as a SNP's, so
  it never distinguished the allele signal either. The rule is **repeat tract against everything
  else** (§3.2), because slippage is what actually stops a read reporting its allele, and it
  happens at tracts and nowhere else. Multiallelic SNPs, insertions and deletions all carry the
  allele signal. Two things fall out: the decision needs no look at the record's alleles, which
  strikes trap 2; and it needs no new field, because `is_repeat_tract` is already in the entry for
  the writer's ordering rule.
- **Scoring a tract on coverage alone is a deferral, not a verdict on its evidence** (the owner,
  2026-09-07: *"at least as a first approximation, because the STR loci with stuttering are very
  noisy"*). §8's tract-aware allele term is the piece that would use it, and D4's counts say
  whether it is worth building.
- **OPEN — `N = 1`.** Leaning: the copied code handles it; step D1 is the check.

## 10. How we know it works

- **Bit-identical ratios**: ng's copied scorer against production's on identical synthetic
  inputs, over the copied tests and a randomised differential.
- **Off is byte-identical**: `--paralog-fdr 0` reproduces the pre-filter VCF on the run fixtures
  and the tomato slice.
- **Unflagged is byte-identical**: filter on at an unreachable target, the two INFO keys
  stripped, equals the filter-off file.
- **Tag equals drop plus the tagged**: the drop-mode file is the tag-mode file with every
  `hiddenParalog` line removed.
- **Mode equivalence with the filter on**: direct mode and psp mode produce the same file
  ([`scripts/ng_mode_equivalence_oracle.sh`](../../../../scripts/ng_mode_equivalence_oracle.sh)).
- **The spill round-trips**: encode → decode equality on entries with absent pairs, `NaN`
  preserved by bits.
- **The ends of the range run**: one tomato accession; six over the two 100 kb intervals
  ([`scripts/ng_fit_stage_end_to_end.sh`](../../../../scripts/ng_fit_stage_end_to_end.sh)'s
  ground); HG002. Each reports records dropped, samples rejected, and π.
