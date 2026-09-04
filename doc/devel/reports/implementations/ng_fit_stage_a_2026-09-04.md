# The fit stage — Milestone A: a census built from a stored psp, and it is the same census

**Date:** 2026-09-04
**Plan steps:** [parameter_prepass_runs.md](../../ng/impl_plan/parameter_prepass_runs.md) Milestone A, steps A1 and A2
**Spec:** `parameter_prepass_joint_records.md` §6.1, §7.12; `parameter_prepass_census_sites.md` §2, §3
**Branch:** `ng-psp-mode`

## The answer

**One sample's census built while its reads are walked, and built again afterwards from the psp
that walk wrote, are the same file byte for byte.** That holds for both samples of the varying
fixture cohort — a 600-base non-repetitive reference carrying a deliberate ten-copy `GT` tract,
two samples varying in different places, three read groups between them.

That is the question the milestone exists to answer, and it is a question about the **psp**
rather than about the new code: a census needs a read's read group, its per-position witness,
its observation count and its length at a repeat tract, and a record format that dropped any of
them would still read back, still call, and still produce a census — one that differs only here.

## What the two producers share, and why that matters

Both build their writer through **one function**, `CensusPlan::writer_for`, which this milestone
factored out of `SampleObservationGatherer::open`. Eight arguments settle what a census records
— the kept loci, the read groups, the contig lookup, the recording terms, the depth ladder, the
two caps — and a second hand-written call site would make the comparison a test of whether two
lists were kept in step rather than a test of the format.

Beyond that the two differ in one thing only: where the loci come from. Both feed
`CensusWriter::add_locus`, and both feed it a `SampleLocusObservations` — what the walk yields,
and what a psp decodes back to.

## What the comparison would catch, measured rather than assumed

Four deliberate defects were written into the psp-driven producer, each run against the
milestone's own tests, each then removed
([`scripts/ng_census_agreement_mutations.sh`](../../../scripts/ng_census_agreement_mutations.sh)):

| the defect | what the tests did |
|---|---|
| the producer skips repeat-tract loci | all 3 fail |
| one read is lost at every locus | 2 of 3 fail |
| every read is credited to read group 0 | 1 of 3 fails |
| one read's minted error arrives one step off | **all 3 pass** |

**The third fails on only one of the two samples, and that is right**: the fixture's second
sample declares a single read group, so re-labelling its reads to group 0 changes nothing about
what it recorded. The first sample carries more than one and separates them.

**The fourth is not a hole in the test — it is a fact about the census**, and it confirms from
the code what the plan reasoned to. A census section holds a depth code for every kept position
and read group plus the non-reference allele observations, and **no per-read quality at all**.
So a change to a read's Σ `ln ε` cannot move a census byte. That is exactly why
`RunParameters::assemble`'s per-read-group minted-error totals cannot come out of a census file
as the format stands, which is the open design question Milestone E carries
([plan](../../ng/impl_plan/parameter_prepass_runs.md) §3.4, step E2).

## The digest has to come off the file, and a test says the two routes agree

A census names the psp it was built from by a digest of that psp's header and its record count.
The walk-time producer takes the digest from what the writer hands back, because `PspWriter::create`
records the compression level into the header **before** encoding it — so a digest of the header
a walk *holds* names a file that does not exist, and every freshness check would answer *rebuild*
for ever, silently. That defect was found and fixed during the walk stage.

The new producer has no writer to ask, so `psp::header_digest` reads the header's bytes straight
off the file. A test asserts that this equals the digest of the same header re-encoded, so a
format change that made a header round-trip to different bytes cannot break the pairing without
saying so.

**The record count is counted while walking**, not read from the file's own index, so a psp whose
index disagrees with its blocks cannot produce a census claiming a count its records do not
support.

## What is not shown here

**Nothing yet reads a census back at command level, and no command builds one from a psp** —
those are Milestone B. The producer is a library function with tests; `generate-psps` still
writes the census it always wrote, unchanged, and both routes stay
([plan](../../ng/impl_plan/parameter_prepass_runs.md) §3.1, the owner's ruling of 2026-09-04).

**The two routes have not been timed against each other.** That is step B3, and it needs the
command.

## Validation

`cargo test --lib` in the container: **6,236 passed, 0 failed, 15 ignored** — the 6,229 the
branch stood at, plus this milestone's seven. `cargo fmt` and
`cargo clippy --all-targets --all-features -D warnings` are clean.
