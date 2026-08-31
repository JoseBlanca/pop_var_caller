# ng VCF module, step A3 — what the header states, and what it may not

**Date:** 2026-08-30. **Branch:** `ng-vcf-output`. **Plan:**
[`../../ng/impl_plan/vcf_output.md`](../../ng/impl_plan/vcf_output.md) step A3.
**Spec:** [`../../ng/spec/vcf_output.md`](../../ng/spec/vcf_output.md) §4.
**Completes Milestone A**, with [A1](ng_vcf_output_a1_2026-08-30.md) and
[A2](ng_vcf_output_a2_2026-08-30.md).

---

## What landed

`src/ng/vcf/header.rs`: the header's *content* — contigs, the sample names in the run's sample
order, the command line, the reference path, and the parameters file written beside the VCF —
with a checked constructor. No text: rendering is step C1.

**45 tests across the module, all passing** (12 here, 33 on the record), with `cargo fmt
--check`, `cargo clippy --all-targets --all-features -D warnings` and `cargo test --lib ng::vcf`
green in the container. The runner reports 48 because the filter `ng::vcf` also matches three of
production's `var_calling::vcf_writer` tests by substring.

## Errors here, assertions there — the same module, deliberately split

The record type panics on a bad record. This one returns a `Result`. The line between them is
**where the bad value came from**:

- every invariant on a record binds two things the same worker built moments apart — the allele
  table against the counts pooled over it, a genotype against the table it indexes — so a
  violation is a wiring defect in this crate and nobody downstream could act on a `Result`;
- header metadata is the run's inputs: sample names come from the alignment files, contigs from
  the reference, the command line from whoever typed it. **Two files naming the same sample is a
  run someone can fix**, so it is an error with the offending name in it.

Production draws the line in the same place, and its header builder returns
`Result<(), VcfWriteError>` for these same checks.

## The refusals, and why each one is not merely tidiness

- **No samples.** A VCF names its samples in the `#CHROM` line and a cohort has at least one.
- **Two samples of one name.** The sample columns are *positional*, so nothing reading the file
  could tell the two apart. That makes the output ambiguous rather than wrong, which is worse:
  a wrong file can be spotted.
- **An empty sample or contig name** — a column heading that names nothing, and a contig no
  record's `CHROM` could name.
- **Two contigs of one name.** Every `CHROM` would be ambiguous.
- **A contig longer than `i32::MAX`.** A VCF states a contig length as a 32-bit signed integer,
  so a longer one cannot be written honestly and is refused rather than truncated into a
  plausible smaller number. **This catches a corrupt index, not a large genome:** human
  chromosome 1 is about 249 million bases against a ceiling of about 2.15 billion, more than
  eight times the largest real chromosome — which is asserted in a test rather than claimed
  here, so the sentence cannot go stale.

**An empty contig list is accepted**, matching production: a run whose reference states no
contigs has nothing to say about them, which is a strange run rather than an unwritable header.

## Two things the type decides rather than stores

**`##source` is derived**, from the crate version at compile time, so a file cannot claim to
have been written by something other than the binary that wrote it.

**`##parametersFile` is a file name, not a path.** The VCF and its parameters file travel as a
directory; an absolute path would be stale the first time the pair moved. It is the line that
makes a run reproducible from its own output, and **neither production writer has anything like
it** — the SNP/indel one carries `##source` and `##commandline`, the repeat-tract one carries no
provenance at all.

## What the contig digest says, and when it says nothing

`HeaderContig.md5` is `Option`, and `None` is honest rather than missing: a run driven from a
`.fai` alone never read the bases, so it has no digest to state and the attribute is left off
the line rather than invented. That is a third behaviour, not a copy of either production
writer — its SNP/indel writer always states the digest, its repeat-tract writer never does.

## One thing this type cannot check, stated where it is asked for

The sample names must be **in the run's sample order** — the same order every record's sample
columns are in. Nothing here can verify that: a permutation is still a list of distinct names.
It is the run's to get right, and the requirement is written into the parameter's own
documentation rather than left implied, because the failure is silent — every genotype in the
file would land under the wrong sample name, and the file would parse.

## Milestone A is complete

The checkpoint asks that the types express both record kinds, a no-call and a filtered tract
locus, and that a record the spec forbids is unrepresentable or refused. All three steps'
reports record how. **One design question is open and is the checkpoint's to settle** — which
called samples become `./.` — carried from A1's review and stated in the plan's Checkpoint A
note.
