# The fit stage — Milestone C: a parameters file produced from data

**Date:** 2026-09-05
**Plan steps:** [parameter_prepass_runs.md](../../ng/impl_plan/parameter_prepass_runs.md) Milestone C, steps C1–C6, and Milestone E brought forward
**Spec:** `parameter_prepass_joint_records.md` §6.1, §6.2; `parameters_file.md`; `ordinary_site_prior_moments.md` §6
**Branch:** `ng-psp-mode`

## The answer

**`estimate-parameters` fits a cohort from its census files and writes the parameters file a
calling run scores with.** Before this, a run had two sources for its numbers and neither was a
fit: the constants compiled into the binary, or a file somebody handed it.

The file names the reference, the samples, the read groups and the census its numbers came from,
and one cohort fitted twice writes the same bytes.

## The two things that had to be built first, and neither was in the plan

**A cohort of censuses could not be assembled at all.** Every census numbers its read groups from
zero, because a walk sees one sample — so on a two-sample cohort both claim read group 0 and
`CohortCensusEvidence::new` refuses them, correctly, as libraries that would be fitted as one.
That is psp mode's normal state and not a corner: the advertised way to walk a cohort is one
invocation a sample. It was measured on the fixture before anything was changed.

The owner's ruling was to put the read groups' names in the census — the `@RG ID` and the
library — so a cohort merges on those and the fit's input stays the censuses (C1). The merge then
renumbers onto run-wide identifiers in (sample order, group order), the rule
`ReadGroups::of_merged_tables` already uses for alignment files, and refuses instead the thing no
run can produce: two samples declaring one `@RG ID` (C2).

**`RunParameters::assemble` refuses a fitted error rate with no minted read-error total**, where
plan §3.4 had assumed the pair would fall back to a defaulted calibration. It does not — it
panics, and C5's own tests found it on their first run. The owner's answer was to bring Milestone
E forward: the census now accumulates the totals as its loci go past, and the base-quality
calibration is **fitted** rather than defaulted.

## Where the accumulation went, and why it is better than the step asked for

Milestone E's step said to accumulate the totals while `generate-census` reads the psps. It went
into `CensusWriter` instead. **Both producers feed that writer**, so the walk-time census and the
psp-built one accumulate the same totals — and §7.12's byte-for-byte agreement, which still
holds, now checks them too.

It reuses `minted_error_by_read_group` unchanged, so the number means what it means on the
per-sample route: a read at a position counted once for every position it is seen at, at generic
loci, over complete witnesses, before the per-position depth cap.

**Every generic locus the walk hands over, not only the kept ones.** Restricting it to the
census's kept positions would be a second definition of a per-read-group total, and the one number
it feeds — how far a library's own base qualities may be trusted — is a property of the library
rather than of which positions were kept.

## Two round trips that were not symmetric, and how they showed

**The totals live in a table of their own** rather than as a field beside the read-group names.
Writing one entry per declared group turns *no entry* into *an entry that saw no reads* on every
round trip, and those are different claims: a group nothing was accumulated for, against a group
whose reads all carried a quality of zero.
`every_corner_state_survives_a_round_trip` caught exactly that — a hand-built value whose map was
empty came back full of zeros.

**A census that names one read group twice is refused**, not read last-one-wins: two entries under
one identifier would name a section's read group two ways, and only one could reach a cohort's
merge.

## What the fit needs besides the censuses

**The psps, but only their headers.** A census names the psp it was built from, and a census whose
psp is not beside it is refused rather than trusted — nothing can then check it and nothing can
rebuild it. What is taken from each is one short read: the header's digest, and the ground the
walk covered.

**The check is on the header alone, and that is stated rather than implied.** A census names its
psp by the header digest *and* the record count, and a psp carries no record count anywhere a
reader can reach cheaply — its footer holds block and byte counts, and the per-block counts sit
inside each block's compressed stream. `freshness_by_header` says so in its own documentation, and
names what it leaves out: a psp whose header is unchanged and whose records are not, which only
`PspWriter::append` can produce and nothing in the shipped commands calls.

**The reference and the catalog**, because a census stores a repeat tract by its index within its
stratum and nothing else. The selection has to be rebuilt, and **the rebuild is checked against a
digest every census carries** before a tract is read. That guard caught the first version of C4's
own test, which rebuilt at the fixtures' 300-position budget while the walk had used the shipped
two million.

## What the file still declares rather than fits

**The inbreeding coefficient.** It is fitted from a sample's own windowed genome histogram, which
is the other pre-pass route. `--inbreeding` states one for the cohort and the file records it as
`supplied`, so a reader can tell it from a measurement.

## A fixture fact worth carrying

**The plain on-disk cohort has no repeat tract at all.** Its reference is all `A`, so every base is
a homopolymer, the whole genome routes to the repeat path, and the selection keeps nothing —
0 strata over 0 tracts, measured. Every test of the fit and of `estimate-parameters` uses the
varying cohort, which carries a deliberate ten-copy `GT` tract.

## Validation

`cargo test --lib` in the container: **6,271 passed, 0 failed, 15 ignored** — the 6,229 the branch
stood at when this plan began, plus this plan's forty-two. `cargo fmt` and
`cargo clippy --all-targets --all-features -D warnings` are clean.
