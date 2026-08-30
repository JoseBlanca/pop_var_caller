# ng — the VCF a run writes: one file for SNPs, indels and repeat tracts

*Design spec, 2026-08-30. **No code yet — this settles the format.** One document, one module:
`src/ng/vcf/` (the writer's stream wiring — which stage feeds it — stays in the run's own
documents; see §1.2).*

*Reads on: [`calling_quality.md`](calling_quality.md) §3.5 (the one quality field this file's QUAL
column is), [`cohort_merge.md`](cohort_merge.md) §4.2 (the allele table and per-sample support every
record is written from), [`typed_regions.md`](typed_regions.md) §2.3 (the reference partition that
keeps two records off the same base), and [`parameters_file.md`](parameters_file.md) §7 (the file
written beside this one). Read by: whoever builds emission — steps 11a and 11b of
[`ng_proposal.md`](ng_proposal.md) — which own what is *dropped*; this document owns what an
emitted record *looks like*.*

*Production's equivalents are two separate writers producing two separate files:
[`src/vcf/header.rs`](../../../../src/vcf/header.rs) +
[`record_encode.rs`](../../../../src/vcf/record_encode.rs) (the SNP/indel file) and
[`src/ssr/cohort/vcf_out.rs`](../../../../src/ssr/cohort/vcf_out.rs) (the repeat-tract file).
Everything said about them is a record of what they do, not a proposal to change them — `src/vcf/`,
`src/ssr/` and `src/var_calling/` are frozen production.*

---

## 1. What this is

**Production writes two VCFs per cohort: one from the SNP/indel caller, one from the repeat-tract
caller. ng writes one.** The two production files were never designed against each other — they
disagree on which depth `DP` counts, on whether a filtered locus is written or vanishes, on whether
a thin sample is no-called or force-called, and on what the header says about where the file came
from. A user who wants "the variants" must merge them with outside tooling, and the merge has no
rules because nobody wrote any.

ng removes the reason for the split. Both of its locus kinds already share one internal shape — a
reference span, a table of sequence alleles, and per-sample support against that table
([`cohort_merge.md`](cohort_merge.md) §4.2) — and the reference partition guarantees a repeat
tract's bases and a generic locus's bases are never the same bases
([`typed_regions.md`](typed_regions.md) §2.3). So one record shape covers both, records interleave
in genome order, and a repeat-tract record differs from a SNP record only by carrying a few extra
fields.

### 1.1 Goals

1. **One file, one record shape.** A repeat-tract record and a SNP record differ in annotations,
   never in structure. Any VCF-consuming tool that understands one understands the other.
2. **Every disagreement between the two production formats is resolved by a stated rule**, not by
   whichever writer happened to run. §2 lists the disagreements; §5–§8 resolve each.
3. **A record carries what a downstream filter needs**, in the annotate-and-defer style of
   freebayes, GATK and GangSTR — the two quality penalties, the mapping-quality contrasts, the
   repeat annotations. What ng drops inline is step 11's decision; what survives must not need a
   re-run to be re-filtered.
4. **The file identifies its run.** Reference, command line, caller version, and the parameters
   file written beside it ([`parameters_file.md`](parameters_file.md) §7) — production's tract
   file carries none of these.
5. **Degrade across the committed range** — one sample to several thousand, three reads a position
   to several hundred (`CLAUDE.md`). §10.
6. **Byte-identical output at any worker count**, the run-wide rule
   ([`run_streaming.md`](run_streaming.md) §1 goal 4), which for a text format means the encoding
   of every number is part of the format (§11).

### 1.2 Non-goals, and what this document does not do

- **It does not decide what is dropped.** The emission threshold, the inline artifact drops, and
  the hard-drop-versus-annotate choice are steps 11a and 11b
  ([`ng_proposal.md`](ng_proposal.md)). This document defines the FILTER vocabulary those steps
  write into (§8) and the invariant they must keep (§5, QUAL), nothing more.
- **It does not compute any value.** QUAL, GQ and the penalties are
  [`calling_quality.md`](calling_quality.md)'s; AF/AC/AN come off the calling loop; AD/DP come off
  the merge. This document says how each is *spelled*, and which record carries which.
- **It does not define the tract quality.** `calling_quality_ssr.md` is unwritten; when it exists,
  its numbers land in the same columns with the same meanings (§5).
- **No gVCF.** Reference-confidence blocks, `<NON_REF>`, per-position invariant output — a
  different output mode with its own memory shape. Deferred with a home (§13).
- **No phasing.** ng computes none ([`ng_proposal.md`](ng_proposal.md) step 10 — a genuine gap,
  owned there). Every genotype is written unphased, `/` never `|`.

### 1.3 Vocabulary

- **generic record** — a record from the SNP/indel path: a cohort locus the generic caller scored.
  Covers SNPs, indels, and multi-variant loci chained by overlap
  ([`cohort_merge.md`](cohort_merge.md) §4.1).
- **tract record** — a record from the repeat-tract path: one catalogued repeat tract, its alleles
  full tract sequences ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)).
- **the run's sample order** — the fixed order the run holds its samples in, the same order every
  per-sample column and every per-sample loop uses
  ([`run_streaming.md`](run_streaming.md) §6.2). The sample columns of this file are that order.
- **no-call** — a sample the record carries no genotype for, written `./.` with every per-sample
  field it cannot honestly fill written `.`. Distinct from a confident `0/0`, and never
  substituted for one.

---

## 2. The two production formats, and where they disagree

What follows is the record of what exists, condensed to the points a unified format must rule on.
Full field lists are in the two writers' own code (§ preamble links).

**The SNP/indel file** (`src/vcf/`): VCFv4.4 via noodles. INFO `AF,AC,AN,DP,CA,MQRef,MQAlt,MQDiff,
MQDiffT,PARALOG_POST`; FORMAT `GT:GQ:DP:AD` (+`GP` opt-in); FILTER only `PASS`/`EMNoConv` — every
other rejection is a silent drop. Header carries `##source`, `##commandline`, `##paralogFilter`,
contigs with md5. Indels anchor-based and left-aligned. Every sample always gets a concrete
genotype — there is no no-call. `DP` (per sample) is the sum of `AD`, by construction.

**The repeat-tract file** (`src/ssr/cohort/vcf_out.rs`): VCFv4.4 hand-formatted. INFO `PERIOD`
only; FORMAT `GT:GQ:REPCN:DP:AD`; FILTERs `notPeriodic`, `tooManyAlleles`, `lowDepth` — a filtered
locus **is written**, with every sample no-called. Header carries contigs (no md5) and a warning
line; no source, no command line, no reference. Alleles are full tract sequences; a full-tract
deletion takes the base to the left as anchor, or the letter `N` when the tract starts the contig.
Thin or imbalanced samples are no-called (`./.`). `DP` is all spanning reads, so `DP − ΣAD` is the
reads no called allele explains.

The disagreements, each resolved at the section named:

| the two files disagree on | SNP/indel file | tract file | ng — resolved at |
|---|---|---|---|
| can a sample be no-called? | never | yes, `./.` | yes, both kinds — §7 |
| per-sample `DP` counts what? | reads on the record's alleles (`= ΣAD`) | every spanning read (`≥ ΣAD`) | every read observed — §7 |
| a locus that fails a locus-level check | vanishes | written, on its FILTER | format supports both; policy is step 11's — §8 |
| header provenance | source + command line | none | full, plus the parameters file — §4 |
| anchor for an empty allele | always present (anchor-indel form) | left base, or an invented `N` | left base; right base at contig start; never `N` — §5 |
| QUAL when in doubt | recomputed at write time (the shipped defect) | computed at write time | written once upstream, never recomputed — §5 |

---

## 3. One record shape — the decision

**Every record is: a reference span, sequence alleles over that span, one site quality, and
per-sample support.** A tract record is that shape plus repeat annotations; nothing else marks it.

**Alleles are literal sequences everywhere.** Both production paths already do this (the tract file
follows HipSTR's convention), and the alternative — GangSTR's symbolic `<STR12>` ALTs — was
rejected because it cannot spell two same-length alleles that differ in sequence, which are 43% of
HG002's heterozygous tracts and a class ng's candidate selection was specifically built to keep
([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §5). A consumer that only wants lengths
reads `REPCN` (§7) and never parses a sequence.

**A tract record is marked by the `STR` flag, and its annotations travel with it.** One INFO flag
(`STR`), always accompanied by `RU` and `PERIOD` (§6) — GATK's own convention for repeat records
(its `STR` flag beside `RU`), so tools that know that output already understand this one, and
selecting either kind of record is one expression (`bcftools view -i 'STR=1'` for the tracts,
`-e` for the rest). The flag and the annotations are written from the same branch of one writer,
so they cannot disagree. Decided with the owner, 2026-08-30 (§14 Q1). *Rejected:* a
`INFO/LOCUS_TYPE` string enum — an invented spelling where a field convention already exists.

**Records never claim the same base.** Within a path, loci are disjoint by construction
([`cohort_merge.md`](cohort_merge.md) §4.1; tracts are catalogue entries, disjoint by the
partition); across paths, the reference partition is exclusive and complete
([`typed_regions.md`](typed_regions.md) §2.3). The one place two records can share a *POS* — while
still not sharing a base — is the anchor shift, and §5 defines the order there.

---

## 4. The header

In emission order:

```
##fileformat=VCFv4.4
##source=ng {CARGO_PKG_VERSION}
##commandline={argv, joined by single spaces}
##reference={the reference path the run was given}
##parametersFile={the parameters TOML written beside this file, by file name}
##contig=<ID={name},length={length}[,md5={md5}]>        one per contig, reference order
##INFO=...   ##FILTER=...   ##FORMAT=...                the declarations of §6–§8
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{samples, in the run's sample order}
```

- **`##parametersFile` is new — neither production file has it, and it is the line that makes a
  run reproducible from its output directory.** [`parameters_file.md`](parameters_file.md) §7 has
  every run write the file; this line names which one. File name, not path — the pair travels as a
  directory, and an absolute path would be stale the first time it moved.
- **Contig lines come from the run's `ReferenceInfo`** ([`reference_info.md`](reference_info.md)):
  name and length always, `md5` when the reference FASTA was read (a run driven from a `.fai`
  alone has none — the attribute is omitted, not invented).
- **Every FILTER, INFO and FORMAT id the file can carry is declared, including `PASS`.** The tract
  writer's choice to leave `PASS` undeclared is legal VCF but gratuitous; declared costs one line.
- **Sample names are the run's sample order** and must match the merge's — production's header
  builder rejects a mismatch and duplicates ([`header.rs:41-48`](../../../../src/vcf/header.rs)),
  and ng keeps both refusals.
- **No warning lines.** Production's tract file writes `##ssrCallWarning=` (the apparent-F_IS
  notice). In ng that finding is a property of the *fit*, and the parameters file is where a
  fitted number's warrant lives ([`parameters_file.md`](parameters_file.md)); a VCF header is
  where nobody looks for it. The warning still reaches stderr at fit time.

---

## 5. The fixed columns

- **CHROM** — the contig's name, from the run's contig table. An id outside the table is a bug and
  panics (production's rule, both writers).
- **POS** — 1-based start of the record's reference span, after the anchor rule below.
- **ID** — `.` always. Neither production writer fills it; nothing consumes it. (The repeat
  catalog could stamp a stable tract id here one day — deferred, §13.)
- **REF** — the reference bases of the whole span, uppercase `{A,C,G,T,N}`. For a tract record
  that is the whole tract; for a generic record the whole chained locus.
- **ALT** — every non-reference allele in the record's allele-table order, comma-joined, each a
  full-span literal sequence. A multi-allelic site is **one record** — never split, never
  symbolic. `.` when a written record has no alternative (only reachable for a filtered tract
  locus, §8).
- **QUAL** — **the one site quality the record carries, already corrected — written, never
  recomputed.** [`calling_quality.md`](calling_quality.md) §3.5 is the design rule and records the
  shipped production defect (gate on one number, write another: 40–64 false positives on GIAB)
  that makes it a rule. The writer formats the field; it must not derive it. Numeric always, `0.0`
  floor, capped at 9999 upstream; for tract records the number is `calling_quality_ssr.md`'s to
  define and this column's meaning does not move.
- **FILTER** — §8.

**The anchor rule — one rule for both kinds.** VCF cannot write an empty allele. Whenever any
allele of the record would be empty (an insertion's or deletion's nature on the generic path; a
full-tract deletion on the tract path), every allele is prefixed with the reference base
immediately to the **left** of the span and POS moves one left. When the span starts at position 1
of the contig, the padding base is instead the reference base to the **right**, appended — the
VCF 4.4 rule for events at position 1. **Production's tract writer instead invents the letter `N`
there** ([`vcf_out.rs:405-435`](../../../../src/ssr/cohort/vcf_out.rs)) — a base the reference
does not contain, at an unshifted POS; not ported. Generic alleles arrive already left-aligned
([`cohort_merge.md`](cohort_merge.md) §4.2 — unification depends on it), so the anchor here is
presentation, not normalisation.

**Ordering, and the one legal POS tie.** Records are written in (contig order, POS) order,
**non-decreasing**, and two records may share a POS only through the anchor shift: a tract whose
record moved one left can land on the POS of the generic locus that owns the anchor base. At most
one generic and one tract record can meet this way (loci are disjoint within each path), and the
**generic record is written first** — its span genuinely starts there; the tract's starts one
base later. Production's generic writer enforces *strictly* increasing POS
([`writer.rs:216-227`](../../../../src/vcf/writer.rs)) and could, because it never shared a file
with the tract records; ng's writer relaxes exactly this far and no further — a third record at
the same POS is a bug.

---

## 6. INFO — what a site carries

Declarations, with what the value is and which records carry it:

```
##INFO=<ID=AF,Number=A,Type=Float,Description="Fitted frequency of each ALT allele, from the calling loop's converged pass">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Copies of each ALT allele in the called genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total called allele copies (no-call samples excluded)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Sum of the samples' DP">
##INFO=<ID=ABPEN,Number=1,Type=Float,Description="Phred subtracted from QUAL by the allele-balance artifact test">
##INFO=<ID=SPPEN,Number=1,Type=Float,Description="Phred subtracted from QUAL by the strand and read-position artifact test">
##INFO=<ID=MQREF,Number=1,Type=Float,Description="Cohort-pooled mean mapping quality of reads supporting REF">
##INFO=<ID=MQALT,Number=A,Type=Float,Description="Cohort-pooled mean mapping quality of reads supporting each ALT">
##INFO=<ID=MQDIFF,Number=A,Type=Float,Description="MQALT minus MQREF per ALT; negative means ALT reads map worse (multi-mapper fingerprint)">
##INFO=<ID=STR,Number=0,Type=Flag,Description="This record is a repeat-tract locus">
##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit of the tract, reference strand">
##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat unit length in bases">
```

- **`AF`, `AC`, `AN`, `DP`** — the standard quartet, semantics as production's generic file: `AF`
  is the loop's fitted frequency (not derived from genotypes), `AC` is counted from the called
  genotypes, and `AC`-consistency with `AN` is asserted, not hoped. **One change: `AN` counts
  called samples only.** Production's generic `AN` is always `samples × ploidy` because it cannot
  no-call; ng can (§7), and an `AN` that counted absent genotypes would be wrong at every thin
  locus.
- **`ABPEN`, `SPPEN` — the two artifact penalties become annotations, and this is the decision
  [`calling_quality.md`](calling_quality.md) §12 left to this document.** Emitted on every generic
  record the correction ran on (0.0 when a test charged nothing). Two reasons: the uncorrected
  quality stays recoverable as `QUAL + ABPEN + SPPEN` wherever the floor did not bite — the §3.5
  rule keeps one *authoritative* number without hiding what was subtracted; and it is the
  annotate-and-defer convention every other caller follows
  ([`ng_proposal.md`](ng_proposal.md) step 11a), at a cost of two floats a record. Absent on tract
  records until the tract quality document defines its own corrections.
- **`MQREF`, `MQALT`, `MQDIFF`** — ported from the generic file; the inputs (per-observation
  `mapq_sum` over `num_obs`) are already on every observation, both paths, so tract records carry
  them too — production's tract file simply never had them. A key whose pooled denominator is zero
  is omitted (`MQREF`) or written `.` per entry (`MQALT`, `MQDIFF`), production's rule.
  **`MQDiffT` (the Welch's t companion) is not ported** — two spellings of one contrast, and the
  production filter built on this family is the one that costs SNP recall at the knee
  ([`mapq_diff_filter`](../../reports/) history); one honest contrast plus a downstream threshold
  is enough surface.
- **`STR`, `RU`, `PERIOD` — the tract marker and its annotations, always together (§3).** The
  flag is GATK's convention and exists so that pulling either kind of record out of the file is
  one filter expression, not a presence test on a string field. `PERIOD` is production's one
  tract INFO field, kept (HipSTR's tag). **`RU` is new: production never writes the motif itself
  anywhere in its file** — a consumer wanting it must re-scan the reference — and the motif is on
  every catalogue entry already. `PERIOD = len(RU)` is redundant by construction and kept anyway:
  it is the field dumpSTR-style tools key on, and one integer per tract record is the whole cost.

**Not carried over, by name:** `CA` (production's chain-anchor flag — the chain machinery has no
ng counterpart); `PARALOG_POST` (the hidden-paralog score — ng has no paralog filter yet; the id
is reserved for step 11a's document, §13); `MQDiffT` (above).

---

## 7. FORMAT — what a sample carries

```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, unphased">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred probability the called genotype is wrong, capped at 99">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Reads this sample observed at the locus, whether or not a written allele explains them">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Reads whose sequence matched each allele exactly, REF first">
##FORMAT=<ID=REPCN,Number=.,Type=Integer,Description="Repeat copy number of each called allele, GT order">
```

The FORMAT string is `GT:GQ:DP:AD` on a generic record and `GT:GQ:DP:AD:REPCN` on a tract record.

- **`GT`** — allele indices in canonical non-decreasing order, `/`-joined, ploidy from the run's
  parameters. Decoding through the shared genotype-order table
  ([`genotype_table.rs`](../../../../src/ng/calling/genotype_table.rs)) — the same enumeration
  `GP` would use if it is ever added, and the same at any ploidy.
- **`GQ`** — the number [`calling_quality.md`](calling_quality.md) §4 defines, integer-rounded,
  cap 99.
- **`DP` — every read observation the sample had at the locus.** The tract file's semantics,
  extended to both kinds, and a real change to the generic one: production's generic `DP` sums
  only the record's alleles, so a sample whose reads half-support an allele that candidate
  selection truncated shows a clean `DP == ΣAD` and the disagreement is invisible. Under this
  rule `DP − ΣAD` is, on every record, *how many of the sample's reads no written allele
  explains* — stutter at a tract, a dropped candidate or noise at a generic locus — which is
  exactly the per-sample artifact signal a downstream filter can use. The number is in hand: the
  merge holds the sample's full observation list before selection narrows it.
- **`AD` — reads whose observed sequence matched the allele, exactly.** Counts, not model
  responsibilities. **This is a deliberate break from production's tract file**, whose `AD` is
  the EM's per-allele responsibility split, rounded — a model output dressed as a count, which
  changes when the genotype does and fills only the called alleles' slots. ng's merge carries
  honest per-allele read counts for every allele of the table
  ([`cohort_merge.md`](cohort_merge.md) §4.2); the file reports those. A read showing a sequence
  the record does not carry is in `DP` and in no `AD` slot.
- **`REPCN`** — tract records only: each *called* allele's tract length in whole repeat units,
  in `GT`'s (sorted) index order — production computes it identically but leaves it in candidate
  order while sorting `GT`, so the two fields' entries need not correspond; ng aligns them. What
  a length-space consumer reads instead of parsing sequences (GangSTR's tag, and the field the
  benchmark concordance already maps HipSTR's `GB` onto).
- **A no-call sample writes `./.` with `.` in every other field it cannot fill** — `GQ` and
  `REPCN` always `.` (no genotype, no quality of it), `DP` and `AD` still written when the sample
  had reads (the evidence exists; the call does not). Three ways a sample gets here, one
  spelling: no coverage at the locus; ruled uncallable by candidate selection
  ([`calling_em_loop.md`](calling_em_loop.md) §5.0); or — when step 11 adopts it — a per-sample
  quality floor, the tract file's `no_call_gq` mechanism. **A sample is never force-called to
  `0/0` for lack of evidence** — the run-streaming leaning
  ([`run_streaming.md`](run_streaming.md) §11 Q5), made a format rule here.

**Not carried over:** `GP` (opt-in posterior rows — the loop's posteriors live in reused scratch
and retaining them is a real memory decision; deferred, §13).

---

## 8. FILTER — the vocabulary, and whose the policy is

```
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=EMNoConv,Description="The calling loop did not converge within its pass cap">
##FILTER=<ID=notPeriodic,Description="Tract allele-length distribution inconsistent with the motif period">
##FILTER=<ID=tooManyAlleles,Description="More candidate alleles segregate than the caller admits">
##FILTER=<ID=lowDepth,Description="Insufficient cohort depth to call the tract">
```

Five ids: production's four (one generic, three tract — the tract three are
[`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §7's verdicts) in one namespace, plus
`PASS`. A record carries exactly one value; there is no `.` FILTER.

**The format's rule: a locus-level refusal that is written at all is written on its FILTER, with
every sample no-called, `ALT .` where no alternative was established, and `QUAL 0`.** That is the
tract file's convention, and it is the right one for a locus the caller *looked at and could not
call* — a vanished locus is indistinguishable from an uncovered one, a filtered record says "here,
and don't trust it". **Whether each refusal is written or dropped is step 11's decision, not
this document's** — production answers it both ways today (the generic file drops everything, the
tract file writes its three), and [`ng_proposal.md`](ng_proposal.md) step 11a flags exactly this
annotate-versus-drop split as a choice to make deliberately. The vocabulary above is what that
decision writes into; new ids (an artifact drop that becomes a tag, a paralog verdict) are added
by the step that creates them.

---

## 9. Which loci appear at all

- **Invariant loci are never written.** Both production files agree, ng keeps it: no gVCF, no
  reference blocks (§1.2, §13). A locus every sample matched the reference at was already
  discarded by the merge's variability filter ([`cohort_merge.md`](cohort_merge.md) §4.3).
- **A locus called all-hom-ref by the loop is not written** (it established no variant; its
  absence means "nothing here", exactly as production's generic path counts it).
- **Below-threshold and artifact-dropped loci** — step 11b's threshold reads the QUAL this file
  would have written ([`calling_quality.md`](calling_quality.md) §12), and what it drops never
  reaches the writer. Not this document's rule; named so the reader knows where it lives.
- **Filtered-but-written tract loci** — §8's shape, if step 11 keeps production's behaviour.

---

## 10. One sample and three thousand, three reads and three hundred

**One sample.** Nothing branches. `AN = 2` (or `0` with the sample no-called — a legal, honest
record if a filtered tract locus is written there). The tract FILTERs still fire — their depth and
periodicity tests are cohort-denominated but defined at `N = 1`
([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §9).

**Three thousand samples.** A record's line is dominated by its sample columns — order of 20
bytes a sample, so ~60 KB a record, linear in `N` and bgzf-compressed on disk (§11). Nothing in
the format is quadratic and no header line grows with the cohort beyond the `#CHROM` line itself.
The real pressure at that scale is *which loci are written*, and that is step 11b's threshold, not
the format.

**Three reads a position.** No-calls are common and honest (§7) — this is where "never force-call
`0/0`" earns its place, and where `AN` counting only called samples keeps `AF`'s denominator
meaning something. `ABPEN`/`SPPEN` sit at or near zero by the ramp
([`calling_quality.md`](calling_quality.md) §7) and the file says so rather than omitting them.

**Three hundred reads.** The penalties carry the weight, and `DP − ΣAD` is large exactly at the
artifact shapes worth flagging. Wide tract records stay one line; nothing in the format degrades
with depth.

---

## 11. Cross-cutting: determinism, encoding, file mechanics

- **Byte-determinism is part of the format.** Same run inputs, same bytes, at any worker count.
  Consequence: every float's rendering is fixed — QUAL to one decimal (both production writers'
  choice), `AF`/`MQ*`/penalties to a stated precision the arch doc pins — and every list order in
  this document (allele-table order, GT-sorted order, the run's sample order) is load-bearing,
  not advisory.
- **Missing is `.`, absent is absent.** A key whose value is undefined for this record is omitted
  (`MQREF` with no reference reads); a defined key with an undefined *entry* writes `.` in that
  slot (`MQALT` for an ALT nobody's reads reached). Production's generic rule, kept, because the
  two mean different things to a parser.
- **Writing.** The writer consumes the run's genome-ordered record stream
  ([`run_streaming.md`](run_streaming.md) §10) and writes as it reads, enforcing §5's ordering
  rule. Bytes go to `{output}.tmp`, then flush, bgzf EOF block where applicable, fsync file and
  parent directory, atomic rename — production's sink, ported
  ([`sink.rs`](../../../../src/vcf/sink.rs)). bgzf when the name ends `.vcf.gz`/`.vcf.bgz`,
  plain text otherwise.
- **Record encoding may be parallelised** (production formats records off-thread and writes in
  order, [`writer.rs:401-425`](../../../../src/vcf/writer.rs)); an implementation detail with one
  format consequence — formatting must be a pure function of the record, which the determinism
  rule already requires.

---

## 12. Reuse map

| what | production code | how ng reuses it |
|---|---|---|
| header building over noodles, metadata validation (dup samples, dup contigs, contig overflow) | [`header.rs`](../../../../src/vcf/header.rs) | ported; declarations replaced by §4/§6–§8's set |
| GT decode through the genotype-order table; AC/AN tally with the range assertion | [`record_encode.rs:485-514,614-625`](../../../../src/vcf/record_encode.rs) | ported unchanged over ng's table |
| INFO omission rules for the MQ family | [`record_encode.rs:429-474`](../../../../src/vcf/record_encode.rs) | ported; `MQDiffT` left behind (§6) |
| tract record shape: full-sequence alleles, anchor on empty allele, `REPCN` | [`vcf_out.rs:372-505`](../../../../src/ssr/cohort/vcf_out.rs) | shape ported into the one record shape; the `N` anchor fallback replaced by the position-1 right-anchor rule (§5); `REPCN` re-ordered to match GT (§7) |
| dense no-call row (`./.` placeholder per absent sample) | [`vcf_out.rs:447-450`](../../../../src/ssr/cohort/vcf_out.rs) | ported as §7's no-call spelling, minus the tract file's "DP also `.`" — ng writes evidence it has |
| ordered writer, atomic tmp-rename sink, suffix-switched bgzf | [`writer.rs`](../../../../src/vcf/writer.rs), [`sink.rs`](../../../../src/vcf/sink.rs) | ported; strict-increasing POS relaxed to §5's one legal tie |
| one function for the gated and written QUAL | [`record_encode.rs:260-278`](../../../../src/vcf/record_encode.rs) | **the property, not the mechanism**: ng writes the stored corrected quality and recomputes nothing ([`calling_quality.md`](calling_quality.md) §3.5) |

**Not reused, by name:** the tract file's EM-deconvolved `AD` (§7); its `N` anchor (§5); its
`##ssrCallWarning` line (§4); the generic file's silent-drop-only FILTER philosophy (§8 — the
vocabulary is format, the policy is step 11's); `PVC_*` environment overrides — anything that
changes this file's bytes is typed configuration, the repository rule.

---

## 13. Deferred, with a recommended home

- **gVCF / reference-confidence output.** A second output mode, not a variant of this one — it
  changes which loci exist, not how they are spelled. **Home:** step 11b's document, which owns
  emission; this format gains `<NON_REF>` machinery only if that document asks.
- **`PARALOG_POST`.** Reserved id, semantics as production's. **Home:** step 11a's document, with
  ng's paralog filter when one exists.
- **`GP` (per-sample genotype posteriors).** Wants the loop to retain per-sample posterior rows
  it currently overwrites ([`calling_quality.md`](calling_quality.md) §3.1) — a memory decision,
  not a formatting one. **Home:** a measurement, then this document gains the field.
- **A stable tract id in the ID column.** The repeat catalog has one per tract; nothing consumes
  it yet. **Home:** here, the day a consumer names itself.
- **Tract-side artifact penalties.** `ABPEN`/`SPPEN` slots exist; what fills them at a tract is
  `calling_quality_ssr.md`'s, unwritten (§1.2).
- **Per-sample analysed-region no-calls** — a sample called over ground it never analysed
  ([`run_streaming.md`](run_streaming.md) §11 Q5). The spelling is §7's no-call; the plumbing
  (per-sample analysed regions reaching emission) is that question's.

---

## 14. Open questions

**Q1 — is one interleaved file what downstream tract tooling can read?** **RESOLVED — one
interleaved file (owner, 2026-08-30), with the `STR` flag added so either kind of record is one
filter expression away** (`bcftools view -i 'STR=1'` / `-e 'STR=1'`, §3). What remains is a
verification, not a decision: run the benchmark concordance tooling
(`benchmarks/lib/ssr_concordance.py`) and one external dumpSTR-class tool over a merged fixture
file — a tool that proves rigid gets a one-line pre-split, not a format change.

**Q2 — does anything downstream miss `MQDiffT`?** OPEN, low stakes. It is dropped on the argument
that one contrast suffices (§6); the June MAPQ-filter investigation used the plain diff, but if a
future filter wants the variance-aware form the sums it needs (`mapq_sum_sq`) are on the
observations and the field can return. **Settled by:** whoever builds step 11a's MAPQ filter,
against measurements, not preference.

**Q3 — `DP` for a no-call sample: evidence or `.`?** Decided as *written when reads exist* (§7),
recorded because production's tract file answers the other way (all five fields `.`). The
argument: a no-called sample with 3 reads and one with 0 are different facts a filter can use, and
blanking both loses it. The cost is a reader that assumed `./.` ⟹ all-`.`; none is known.
**Reopens** only if Q1's survey finds such a reader.

---

## 15. How we know it works

1. **Round-trip under a strict parser.** Every fixture file parses under noodles *and* under
   `bcftools view` with no warning — both kinds of record, the anchor cases, the POS tie, the
   filtered tract locus, the no-call spellings. The tract file's conventions were never pushed
   through an external parser in anger; this is where interleaving would break first.
2. **The partition invariant on output:** over a whole-genome fixture run, no two records share a
   reference base, and the only POS ties are §5's generic-then-tract pairs. Checkable by a linear
   scan of the file alone. The same scan asserts the marker's consistency: every record carrying
   `STR` carries `RU` and `PERIOD`, and no record without it carries either.
3. **Differential against production, where semantics were kept.** On the same cohort, ng's
   generic records against production's file: `GT`, `AC`, `AN`, `AF` agree; `DP`/`AD` differ only
   at loci where selection dropped an allele (the §7 change, visible and countable). Tract
   records against the tract file: alleles and `GT` agree; `AD` differs by the §7 rule; the `N`
   anchor cases differ by the §5 rule; each divergence class counted, none absorbed.
4. **Byte-identity across worker counts** on a cohort fixture — the §11 rule, asserted the way
   [`run_streaming.md`](run_streaming.md) §12 asserts it for the run.
5. **`AN` honesty:** on a fixture with no-calls, `AN + 2×(no-call samples) = 2×samples` (diploid),
   and `AC ≤ AN` per allele — the §6 change, which a force-calling regression would break first.
6. **`DP − ΣAD` is the dropped-allele mass:** construct a locus where selection truncates a
   supported allele and assert the difference lands in every carrying sample's `DP` and no
   sample's `AD`.
