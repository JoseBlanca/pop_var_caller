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
5. **Every written record is scored, on coverage and allele balance**, repeat tracts included
   (§3.2). Repeat tracts are first-class.
6. **Memory flat in the record count and in the sample count**: one record in hand at a time,
   plus per-sample constants and one number per record.
7. **Off means off**: with the filter disabled the run writes byte for byte what it writes today.

### 1.2 Non-goals, and what this does not do

- **No new statistics, with one exception.** The likelihood ratio, the prior's estimation, the
  FDR curve are copied from `src/paralog/`, and every record goes through them as production's
  biallelic SNPs do. The exception is the allele term's dispersion (§3.2, adopted 2026-10-21).
  Production's binomial was tuned on tomato2 at about six reads a sample; at 300x it was checked
  only for overall accuracy on the 100-region GIAB panel, which holds too few duplicated sites to
  show what the split does there.
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
- **A generic locus** and **a repeat tract** — the project's existing pair. A generic locus is
  anything that is not a tract: a SNP of any allele count, an insertion, a deletion. Both are
  scored on coverage *and* allele balance (§3.2); the pair matters to the spill only for the
  writer's ordering rule (§3.4).

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

### 3.2 What is scored, and with what — decision: every record, on coverage and its read split

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

  **Every record is scored on both signals, repeat tracts included** (owner, 2026-10-10; first
  decided 2026-10-08 on branch `tract-ceiling`, which called each STR segment both ways and was not
  merged — this rule was ported from it alone):
  `alt_reads =` the sum of every alternative's reads, `total_reads = AD[0] + alt_reads`, from the
  record's sample column ([`SampleReadCounts::allele_reads`, `vcf/mod.rs:597`](../../../../src/ng/vcf/mod.rs)),
  however the record was called. A sample with `total_reads == 0` is absent from the score, as in
  production. **Alternatives are pooled**, which the H2 story permits unchanged: it asks whether
  the non-reference share sits at some whole `m/T`, and two copies of three carrying two
  *different* non-reference bases give the same two-thirds as two copies carrying one.

  **Why a repeat tract's split is read, though slippage smears it.** The scorer's allele half
  compares two stories about *which allele a read carries*: under a real heterozygote the reads
  split near half, under a collapsed duplication at some whole `m/T`. Both assume a read's allele
  is observed correctly, and at a repeat tract it is not always: slippage during amplification
  makes a read report a tract one or two units from the molecule it came from. Until 2026-10-08 a
  tract was therefore scored on coverage alone — the same likelihood with the allele counts
  marginalised out (owner, 2026-09-07). But the STR path spells a segment whole, so a SNP inside
  a segment is written as part of an STR record's alleles, and scoring that record on coverage alone
  takes away the allele evidence the same SNP would have as a generic record — above all the
  homozygous-alternative veto, since a sample confidently carrying only the alternative is close to
  impossible under a collapsed duplication. Measured on branch `tract-ceiling`, whose combined STR
  records carried those SNPs the same way: on 63 tomato accessions, 135 SNP genotypes GATK also
  calls were lost to STR records the filter dropped (at SL4.0ch02:14821092 the STR record's
  likelihood ratio was 7.8; the generic record for the same SNP passed); scored from the split, 10
  were. What it does on main's own STR records is measured in the port's report. The run calling every STR segment as generic ground,
  which scored every segment from its split, had been the best of the routings measured on tomato
  and on GIAB's STR segments.

  **What is given up.** An STR record's `AD` counts only reads spanning the whole segment, and a
  read that stuttered onto another listed length is counted to that length. It moves reads between
  alleles, and between reference and alternative where a listed length is the reference's
  neighbour, by an amount that has nothing to do with copy number. How much that costs is
  unmeasured until data with known stutter exists (§8).

  **Why not production's rule.** A collapsed duplication leaves the same coverage footprint at a
  tract as at a SNP, and a filter that never looked at tracts would leave every such tract in the
  file while dropping its SNP neighbours. **What is unmeasured**: whether the depth of the wide
  kinds is biased — and the two wide kinds are biased *differently*, so they are counted separately
  (step D4). A tract's observation depth runs below its read depth
  ([`window_coverage.md`](window_coverage.md) §3.1); a deletion depresses depth across its own span
  in the samples that carry it, by construction.

- **How the split is scored.** Under both stories, a sample's split is scored as a
  beta-binomial with the story's share as its mean and a fixed dispersion ρ = 0.0015: the
  correlation between two reads at one site. A real heterozygote's reads do not split like
  independent coin tosses at exactly ½; library amplification, mapping and the local sequence
  make some sites lean one way. On clean GIAB ground at 300x the lean widens the binomial's
  variance by 1.5 (ρ 0.0012 to 0.0017 in three samples), and no heterozygote of 1,096 fell
  below an ALT share of 0.35. The same ρ applies to the homozygotes' error share ε and to every
  collapsed configuration's `m/T`, because read sampling does not depend on which copy a read
  came from. ρ is fixed: below 30 reads it widens the binomial by under 5% and cannot be
  measured, and the three samples agree. It becomes fitted, per sample from records at ALT
  shares 0.25 to 0.75 with 50 reads or more, if a second high-depth dataset with a different
  library preparation scatters more than twice as widely.

  **What the split may say on its own.** At high depth a split that a single copy cannot
  produce is evidence of another copy's reads even at normal window depth, and the filter
  removes such records: on GIAB at 300x every false SNP it removed sat at 0.70 to 1.28 copies.
  The dispersion bounds that evidence to about 1/ρ reads' worth (670) and no more. It is not set
  larger to make depth the gate again, for the reason production's §7 kept the depth likelihood
  Normal: the tail is the signal. Introgressions stay safe because they are single copy, so
  their reads split like any heterozygote's, which is what the one-copy story now describes.

  **What it costs, stated.** A real heterozygote whose short reads include another copy's
  reads (a sample-specific duplication, a long tandem repeat) looks like a false call to one
  sample's depth and split, and may be removed. GIAB's own short-read callers are filtered at
  such sites. At one sample nothing in this score separates the two.

  *Adopted by the owner 2026-10-21 (§11); until then the split was a strict binomial, as in
  production.*

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
  is_repeat_tract u8         -- 0 or 1
  line            bytes      -- length-prefixed; the record's line, no newline
  samples         varint     -- the run's sample count, dense
  per sample:
    gc_fraction   u32        -- an f32's bits; NaN = absent
    mean_depth    u32        -- an f32's bits
    ref_reads     varint     -- AD[0]
    alt_reads     varint     -- every alternative's reads, summed (§3.2)
```

**The tract flag is there for the writer's ordering rule only** — a tract may share a position
with the generic locus owning its anchor base. Every record's samples have one shape (§3.2). Until
2026-10-09 a tract could be spilled with its window alone, and the flag also chose that shape.

**Why the line and not the record.** The record type has ten fields behind a constructor that
asserts every parity between them ([`VcfRecord::new`, `vcf/mod.rs:187`](../../../../src/ng/vcf/mod.rs));
a codec over it is ten field codecs and a reconstruction that must satisfy every assertion, and a
field the record gains later is a field the codec silently drops. The line is the record's own
bytes: pass three writes them back unchanged except for the two columns the verdict touches, so
byte-identity for every unflagged record holds by construction, and nothing about the record's
shape has to be mirrored. The per-sample side fields are the four numbers the scorer reads,
written so the scorer never parses text.

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
    /// One per sample, in the run's order: both signals; a sample at zero total reads is
    /// absent from the score, as in production (§3.2).
    pub samples: Vec<SpilledSample>,
}

pub struct SpilledSample {
    pub window: WindowCoverage,   // the pair, NaN = absent
    pub ref_reads: u32,
    pub alt_reads: u32,           // every alternative's reads, summed
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
  from three reads is weak evidence and the model already weighs it as such. The allele term's
  dispersion (§3.2) barely acts here (the variance factor is 1.003 at three reads); over a cohort
  the small per-sample changes add up: on 59 tomato accessions at about five reads a sample it
  un-flags 37 of 2,914 records, all within 0.75 of the cut.
- **Three hundred reads a position.** The histogram's depth axis scales with the sample
  ([`window_coverage.md`](window_coverage.md) §3.4), so the fit is not rejected for range as
  production's was on a 100× human sample. HG002 is the plan's check that the fit accepts and the
  verdicts are sane (step D3). **Until 2026-10-17 that scale came from the run's first 10,000
  windows**, and a GIAB run over chromosome 22 alone at 201–223 reads a position had every
  sample's fit refused for range; the scale is fitted from every window since
  ([report](../../reports/ng_coverage_depth_scale_2026-10-17.md)). **At this depth the split
  alone can decide a record at normal window depth, and that is intended** (§3.2): it is how the
  filter removes false calls on duplicated ground whose extra reads did not raise the window's
  depth. The dispersion caps how far.

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

1. **~~Do not skip a zero-read sample at a repeat tract.~~** *Struck, 2026-10-09.* It guarded the
   coverage-alone score of a tract, whose samples carried no read counts; every record now carries
   them (§3.2), so a sample at zero total reads is one that observed nothing, everywhere.
2. **~~The SNP test reads the record's alleles before padding.~~** *Struck, 2026-09-07.* Which
   signals a record carries was decided by `is_repeat_tract`, and since 2026-10-09 every record
   carries both (§3.2), so no test looks at an allele at all and there is nothing to get
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

- **The score's power at tracts and indels**, and whether the tract depth caveat biases it. Home: the plan's D4 measurement first; if it flags nothing or flags wildly, the
  decision in §3.2 comes back here.
- **A tract-aware allele term** — an H2 story for a collapsed copy whose difference is a repeat
  length. Home: a later document under this one, once the stuttering benchmark says whether it is
  needed. **Until then the generic term stands in** (owner, 2026-10-08 and 2026-10-10; §3.2): a
  run writes SNPs inside STR records, and their allele evidence was worth more there than the smear
  stutter adds. Its measure is how many STR records the filter drops on stuttering data, which the
  stuttering benchmark ([`str_stutter_benchmark.md`](str_stutter_benchmark.md)) can give.
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
- **~~Scoring a tract on coverage alone~~** *(superseded 2026-10-08/09, below)* **is a deferral, not a verdict on its evidence** (the owner,
  2026-09-07: *"at least as a first approximation, because the STR loci with stuttering are very
  noisy"*). §8's tract-aware allele term is the piece that would use it, and D4's counts say
  whether it is worth building.
- **An STR record is scored from its read split — resolved** (the owner, 2026-10-10, ported from
  branch `tract-ceiling` where it was decided 2026-10-08; §3.2). It replaces the coverage-alone
  rule of the entry above.
- **OPEN — `N = 1`.** Leaning: the copied code handles it; step D1 is the check.

## 10. How we know it works

- **Bit-identical ratios**: ng's copied scorer against production's on identical synthetic
  inputs, over the copied tests and a randomised differential.
- **Off is byte-identical**: `--paralog-fdr 0` reproduces the pre-filter VCF on the run fixtures
  and the tomato slice.
- **Unflagged is byte-identical**: filter on at an unreachable target, the two INFO keys **and
  the filter's four header lines** stripped, equals the filter-off file. The four are the
  `##paralogFilter=` line and the three declarations §3.5 gives.

  **Amended 2026-09-07 (the owner), because as first written this could not pass.** §3.5's
  declarations are emitted only on a run where the filter ran — they have to be, or the *first*
  oracle above breaks on the header alone, and that one is what the whole plan rests on. So the
  on-run's header carries four lines the off-run's does not, and stripping `INFO` keys does not
  remove them. The two oracles as first written were incompatible; this is the one that gives
  way, because what it is really asserting is that **the filter changes no record it does not
  flag**, and the header is not a record.
- **Tag equals drop plus the tagged**: the drop-mode file is the tag-mode file with every
  `hiddenParalog` line removed.
- **Mode equivalence with the filter on**: direct mode and psp mode produce the same file
  ([`scripts/ng_mode_equivalence_oracle.sh`](../../../../scripts/ng_mode_equivalence_oracle.sh)).
- **The spill round-trips**: encode → decode equality on entries with absent pairs, `NaN`
  preserved by bits.
- **The ends of the range run**: one tomato accession; six over the two 100 kb intervals
  ([`scripts/ng_fit_stage_end_to_end.sh`](../../../../scripts/ng_fit_stage_end_to_end.sh)'s
  ground); HG002. Each reports records dropped, samples rejected, and π.

## 11. The allele term's dispersion — adopted 2026-10-21

*The owner adopted ρ = 0.0015 on 2026-10-21. Evidence:
[`ng_dup_filter_allele_balance_2026-10-21.md`](../../reports/ng_dup_filter_allele_balance_2026-10-21.md).
The text now lives in §1.2, §3.2 and §4.*

**What prompted it.** At 300x the filter removed real heterozygotes at a window depth of about one
copy, on their read split alone (report §2.2). Production's spec §7 says "a lone heterozygote at
*normal* coverage is kept", and the single-sample architecture §3 says the score self-gates because
"the reads are few". Neither holds at high depth, and neither sentence describes ng's filter. With a
strict binomial the highest score a split alone could reach at one copy's depth was 0.7 at 30 reads,
43 at 300 and 209 at 1,000; with the dispersion it is 0.7, 23 and 53.

**What was not changed.** No mapping-quality term (§1.2 stands). No change to the depth
likelihood, the configurations, the prior, or the cut. A larger dispersion that would make depth the
gate again was measured and declined (report §5.2).
