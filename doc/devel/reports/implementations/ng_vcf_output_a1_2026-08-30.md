# ng VCF module, step A1 — what a record carries

**Date:** 2026-08-30. **Branch:** `ng-vcf-output`. **Plan:**
[`../../ng/impl_plan/vcf_output.md`](../../ng/impl_plan/vcf_output.md) step A1.
**Spec:** [`../../ng/spec/vcf_output.md`](../../ng/spec/vcf_output.md) §3, §5–§8.
**Module:** `src/ng/vcf/` (new).

---

## What landed

The carrier type for one output record, and nothing that writes bytes. `VcfRecord` holds
everything the file needs once the locus and its reads are gone: the reference span, the allele
table, the cohort's expected copies per allele, one column per sample of the run, the pooled
mapping qualities per allele, the padding base an empty allele needs, the corrected site
quality, the two artifact penalties, one filter verdict, and — only at a repeat tract — the
motif.

**31 tests, all passing.** `cargo fmt --check`, `cargo clippy --all-targets --all-features -D
warnings` and `cargo test --lib ng::vcf` are green in the container (the test filter also
matches three of production's `var_calling::vcf_writer` tests by substring, which is why the
runner reports 34).

## Why there is a carrier type at all

Three of the file's columns are read counts — `DP`, `AD`, and the pooled mapping qualities
behind `MQREF`/`MQALT`/`MQDIFF` — and **no type downstream of the merge holds them**.
`LocusInference` carries calls and qualities; the reads live in the merge's `SampleSupport` and
are released with the locus. So the counts are summed in the worker while the evidence is in
hand, and travel on the record. Same shape, same reason, as the quality module's nine-number
artifact summary (`calling_quality.md` §3.3).

The review found that **the same argument applies to two things the plan had not listed**, and
both were Blockers. They are the substantive content of this step and are recorded in §"What
the review changed" below.

## Departures from the plan, and why

**1. `converged` is carried as a filter verdict, not as a boolean field.** The plan's A1 lists
`converged` among the record's contents and says dropping an item is a stop-and-ask. The fact
still reaches the output — as `FilterVerdict::EmDidNotConverge` — but not as its own field.
Reason: spec §8 gives a record exactly one `FILTER` value, so a boolean beside it would be a
second spelling of one fact, which is the shape that lets two fields disagree. **Recorded in
the type's doc comment, and raised at Checkpoint A rather than settled here.**

**2. `PERIOD` and the `STR` flag are not fields.** `PERIOD` is the motif's length and the flag
is the presence of the annotation, so Checkpoint A's "a tract annotation without its sibling" is
unrepresentable rather than merely refused — the stronger form of what the checkpoint asks for.

**3. The names are not the plan's.** The plan speaks of the record and its per-sample entries
without fixing identifiers. What shipped: `VcfRecord`, `SampleColumn` (it is a column of the
file, not a sample), `SampleReadCounts` (says what the value is), `FilterVerdict` (the verdict
written into `FILTER`, not a thing that filters records), `SampleCall`, `MapqPool`,
`PaddingBase`, `TractAnnotation`.

## What the review changed

Three reviewers ran in isolated worktrees over the step's diff. Two reported before this
report was written; the third (mutation testing) was still running and its findings will be
folded into the next step if they land after the commit.

### Two Blockers, both "a required column has no home"

**`AF` had nowhere to live.** Spec §6 requires `INFO/AF` on every record and fixes its
provenance: the calling loop's fitted frequency, *not* a number derived from the called
genotypes. `LocusInference::cohort_expected_copies` is that quantity, and its own doc comment
forbids the fallback — a call has already discarded the uncertainty the fit carries. The record
had no such vector and the locus is released after assembly, so step B2 would have had to omit
a required field or write a quietly wrong one. **Fixed:** `expected_copies`, one entry per
allele, checked parallel to the table.

**The anchor rule had no base to anchor with.** Writing a full-tract deletion means padding
every allele with a flanking reference base, and that base lies *outside* the record's span.
Nothing downstream holds reference sequence to fetch it: `ReferenceInfo` carries contig
geometry and digests rather than bases, and a run driven from a `.fai` alone has no sequence at
all — while spec §11 requires the encoder to be a pure function of the record. So step B1 could
not have been written from a `VcfRecord`, and the available fallback is production's invented
`N`, which spec §5 explicitly does not port. **Fixed:** `PaddingBase::{Left, Right}`, resolved
in the worker, carried on the record, and refused unless it is present exactly when some allele
is empty and on the side the span's position allows.

### Three Majors

- **The reference allele's length was never tied to its span.** A `REF` of three bases over a
  one-base region describes different ground than the record's own `POS` claims, and the line
  still parses. Now asserted.
- **The mapping-quality pools could count different reads than `AD`.** In production these are
  literally one field (`record_encode.rs:391-401` sums the same `num_obs` that becomes `AD`), so
  two totals mean an `MQDIFF` computed over reads the `AD` beside it does not describe. Now
  asserted, and two of this file's own fixtures had to be corrected to satisfy it.
- **A doc comment contradicted production.** It said the pools cover "every sample the locus was
  called on"; production pools unconditionally over the cohort, and spec §7 applies the same
  principle in writing a no-called sample's evidence beside its `./.`. Corrected — this one
  matters because step D1 fills the pools from that sentence, and the two readings differ
  wherever no-calls are common, which is the three-reads-a-position corner.

### One Major that is a design question, not a coding slip

**The no-call I reused meant something narrower than the file needs.**
`SampleGenotypeCall::Missing` is defined as *candidate selection cut an allele this sample's own
reads had earned*, on the SNP/indel path only — and a sample with **no reads at all** is
deliberately `Called`, not `Missing`, because the prior alone decides it
(`calling/mod.rs:2624-2630`). But spec §7 says a sample is never force-called for lack of
evidence and names three routes to `./.`, and spec §8 writes a refused **tract** locus with
every sample no-called — a shape `LocusInference::new` actively refuses.

**Resolved for A1 by giving the module its own `SampleCall`**, whose `NoCall` carries the
file's meaning. That is contained: no shared type changed. **What is not resolved is the
conversion** — which called samples become `./.` — and in particular what happens to a sample
the loop called from the prior alone. That is step D1's, and it is raised at Checkpoint A.

### Minors applied

`region` made private (an invariant checked once and then left mutable is checked nowhere);
`penalties` → `artifact_penalties`; `TractAnnotation` now holds the crate's `Motif` rather than
its own `Box<[u8]>`, which removes an allocation per tract record, removes a duplicated
emptiness check, and makes a zero divisor unrepresentable; `MapqPool::mean` takes `self`;
`unexplained_reads` no longer has an `expect` (it stays in `u32` and saturates, with the
invariant named in a `PANIC-FREE` comment); the module doc's list of what the encoder derives
was wrong by four entries and now states the closed list.

### One clippy failure, found and fixed twice

`WrittenRecord::new`'s eight arguments tripped `clippy::too_many_arguments`, which `-D warnings`
makes a build failure. Fixed with an `#[allow]` carrying a stated reason, matching
`LocusInference::new` and two other sites in `calling/mod.rs`; the parameters have distinct
types, so no two can be transposed silently. After the Blocker fixes it takes ten.

A second clippy failure came from the fix itself: `#[must_use]` on `SampleReadCounts::new` made
the `should_panic` test that calls it for its panic alone a lint error. The module's
constructors now carry no `#[must_use]`, matching `LocusInference::new`.

## The factual claims in the doc comments

The third reviewer checked every claim the doc comments make about production and about the
design documents. **Nine of ten hold**: that production recomputes its correction at encode
time; that its emission gate read the baseline while the corrected number was written; that the
window was sixteen days; all four false-positive counts (40 at 30×, 64 at 50×, 14 and 14 after
the repair); that its SNP/indel writer cannot no-call at all; the four production filter names
and which writer each comes from; that its generic writer drops rejections silently while its
tract writer writes them with every sample no-called; that it writes `PERIOD` and never the
motif; and that its `REPCN` truncates. **The tenth was wrong** and is the mapping-quality
pooling scope corrected above.

## How A1 is verified

Checkpoint A asks for two things. Both are met, and by different means:

**Every shape the format must express is built:** both kinds of locus, a sample with reads and
no call beside one with neither, a refused tract locus (reference allele alone, `ALT .`, quality
zero, every sample `./.`), a full-tract deletion away from a contig's start and another at its
first base.

**Every shape it forbids is refused or unrepresentable.** Thirteen refusals, each with a
`should_panic` test naming that assertion's own message — and since a `should_panic` test whose
expected message does not match *fails*, the passing suite is itself the evidence that each
refusal fires for its stated reason rather than incidentally. Two states are stronger than
refused: a motif cannot be empty, and the `STR` flag cannot disagree with `RU`/`PERIOD`, because
neither is a separate field.

## What A1 does not do

No bytes. No header. No encoding, no anchor application, no ordering, no I/O. The padding base
is *carried*, not applied — applying it is step B1, which is now possible.
