# psp mode — F1: `call-from-psps`, and what a run over stored files says about each sample

**Date:** 2026-09-04
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone F, step F1
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §2, §5.3, §6.1, §6.2, §12
**Branch:** `ng-psp-mode`

## The answer

**psp mode calls from the command line, and what it writes is direct mode's VCF.** Six tomato
accessions (SRS3394712, SRS3394687, SRS3394549, SRS3394640, SRS3394615, SRS3394559) over the
first two 100 kb intervals of `benchmarks/tomato1/regions.bed`, walked into psps by
`generate-psps` and then called by `call-from-psps`: **599 records, and every byte of the file
except one header line is `call-from-alignments`' own** — sha256 `fd677c91…` on both sides after
that line is removed. The parameters file written beside it is byte-identical too.

The one line that differs is `##commandline`, which records the command a person typed, and the
two commands are different commands. That is F2's oracle passing early; F2 owns pinning it as a
test.

**What the run says about each sample is what it measured, and the owner ruled what that should
be** (2026-09-04). A psp carries no count of what its walk kept or dropped — those tallies live
in the read cursor of a walk that ran in another process — so the report states how many stored
loci this run drew out of each file and how many reads went into the comparison at one of them,
and it names any file whose walk applied read filters the rest of the cohort did not. On the six
accessions:

```
samples: 6 — 6 whose stored file gave this run loci, 0 whose file held none over this ground
  SRS3394549: 167648 loci read, 3.4 reads a locus compared with the reference
  SRS3394559: 179258 loci read, 2.7 reads a locus compared with the reference
  SRS3394615: 185152 loci read, 3.1 reads a locus compared with the reference
  SRS3394640: 174514 loci read, 2.8 reads a locus compared with the reference
  SRS3394687: 182842 loci read, 2.9 reads a locus compared with the reference
  SRS3394712: 193603 loci read, 8.5 reads a locus compared with the reference
```

**That average is not depth and the line does not call it depth.** The number summed is the
record head's `reads-compared-with-reference`, the keep rule's own denominator, and
[its field](../../../../src/ng/psp/record.rs) lists what it leaves out: reads a filter turned
away, reads the per-position cap discarded, reads that covered the locus and produced no
observation, and reads whose witness stopped inside it. At a repeat tract 40 reads cover but only
22 anchor both borders of, this says 22. The first draft called it *depth* in three places; both
reviews caught it, and the wording now says what it is.

## What was built

**`call-from-psps`** ([`src/pop_var_caller_exp/call_from_psps.rs`](../../../../src/pop_var_caller_exp/call_from_psps.rs)):
one `--psp` per sample or a directory of them, `--parameters`/`--defaults` exactly as direct
mode, the same VCF, parameters file and run report. **No `--regions`**, and that is spec §5.3
rather than an omission: a psp records the ground its walk covered, the cohort is refused unless
every file agrees about it, and a flag here could only let a person contradict the files. A
directory contributes the `.psp` files directly inside it in name order, so two runs naming one
directory open the same cohort in the same order.

**The per-sample tallies.** [`StoredSampleTallies`](../../../../src/ng/run/psp_source.rs) counts,
in the psp source, the loci handed over and the compared reads summed across them —
**incremented where a record is handed over, below every refusal**, so a record the source turns
down is not a locus the run read. [`StoredCohortTallies`](../../../../src/ng/run/psp_caller.rs)
pairs each spent source's counts with its open file, with the pairing asserted on the sample name
the way the source-building loop already asserts its own.

**One report, two modes.** [`RunReport`](../../../../src/ng/run/report.rs) now holds
`WhatEachSampleDid` — the walk's tallies, or the stored files' — and everything that comes out of
the calling itself is stated identically either way. A run over stored files prints the ground it
called over and **does not partition it**: the called / repeat-cluster / long-array base split is
a walk's own region tally, no psp records one, and a zero there would read as *measured and none*.

## Deviations from the plan, recorded

**Three lifts the step's wording did not ask for, all made because a second copy would be a place
for the two modes to drift apart while both kept running** — which is exactly what F2 compares:

1. **The numbers a calling run scores with** moved to
   [`src/pop_var_caller_exp/calling_run.rs`](../../../../src/pop_var_caller_exp/calling_run.rs):
   the parameters and their refusals, the calling-loop settings read from the `NG_*` measurement
   switches, the round width chosen from the cohort's size, the ploidy, the two output refusals,
   the VCF header metadata and the report printer. It is
   [`run_ground`](../../../../src/pop_var_caller_exp/run_ground.rs)'s sibling — that module owns
   the ground a run speaks for, this one the numbers it speaks with — and both keep their
   refusals `#[error(transparent)]` at the commands, so one mistake reads the same however it was
   reached.
2. **The on-disk cohort fixture** moved to
   [`src/pop_var_caller_exp/test_fixtures.rs`](../../../../src/pop_var_caller_exp/test_fixtures.rs),
   which closes a Milestone C carry-forward. It was worth doing rather than merely tidy: the two
   copies were **not the same cohort** — `call-from-alignments`' gave both samples no reads at
   all, so the only test that drove that whole command drove it over a cohort with nothing in it.
3. **The read filters' key prefix** (`READ_FILTER_PROVENANCE_PREFIX`) joined their key list, so a
   reader can pick the family out of a header another build of ng wrote.

**Direct mode is byte-identical across all of it**, which is the point of recording the lifts
rather than escalating them: the six-accession VCF, its parameters file and its whole run report
are unchanged (`tmp/d2_oracle.sh`, 599 records, sha256 `869b8058…`).

## What is knowingly left

- **A psp written against a reference with a different contig count gets the catalog's refusal,
  not its own.** The segmentation is built over the psps' ground before `PspVariantCaller::open`
  compares each file's contig table with the reference, so `check_scope` speaks first and says
  *a region names contig 4, and the catalog holds 2* instead of naming the psp. Still a refusal
  before a block is decoded; only the wording and a wasted catalog scan differ.
- **The per-sample section has no cap.** At the thousand-sample end of ng's range it is a
  thousand lines, and the list of files that held nothing is one long line. Direct mode's
  per-sample section has the same shape, so capping one alone would make the two disagree; it
  belongs to both or to neither.
- **The disagreement line can only fire across ng builds today.** `generate-psps` exposes no
  read-filter flag and walks with the compiled-in policy, so two psps differ here only when the
  builds that wrote them differ — a cohort walked over months, which is the cohort psp mode
  exists to make possible. Whether those settings become flags is Milestone C's open question and
  the owner's call.

## How it is verified

| what | how |
|---|---|
| psp mode produces direct mode's VCF | `tmp/f2_oracle.sh` — 599 records, byte-identical apart from `##commandline`, parameters file identical |
| direct mode unchanged by the lifts | `tmp/d2_oracle.sh` — VCF, parameters file and run report byte-identical |
| the command's refusals fire | six tests driving `run_call_from_psps` itself: output a directory, output in a missing directory, parameters file overwritten, ploidy against the file, another cohort's samples, routing that is not the walk's |
| the supplied-parameters path | `the_parameters_a_previous_run_wrote_bind_to_this_cohort`, and §12.5's sample-list refusal |
| the report's psp branch | nine tests in `ng::run::report`, including both orders of the read-filter comparison |
| the counters sit below every refusal | two tests in `ng::run::psp_source` — the out-of-order refusal and the unbuilt body |
| each sample's counts under its own name | `psp_caller`'s cohort fixture now stores two loci for one sample and one for the other |

`cargo test --lib` in the container: **6,214 passed, 0 failed, 15 ignored**. `cargo test --tests`:
one failure, `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`, which is main's
and pre-dates this branch. `cargo fmt --check` and
`cargo clippy --all-targets --all-features -- -D warnings` clean.
