# psp mode — F2: the oracle that justifies the design

**Date:** 2026-09-04
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone F, step F2
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §12.3, §1.1 goal 1
**Branch:** `ng-psp-mode`

## The answer

**One cohort called two ways gives one VCF, on real reads and on a fixture, and the comparison
can now see four kinds of defect it could not see when it was written.**

On real reads: six tomato accessions over the first two 100 kb intervals of
`benchmarks/tomato1/regions.bed`, walked into psps and called — **599 records, byte-identical to
`call-from-alignments`' own apart from `##commandline`** (sha256 `fd677c91…` on both sides once
that line is removed), with the parameters file beside each identical too. The command is
`scripts/ng_mode_equivalence_oracle.sh <reference> <catalog> <regions> <out-dir> <cram>…`, which
lives in the repository so the run reproduces from a fresh checkout.

In the suite: `pop_var_caller_exp::mode_equivalence`, three tests, comparing the two routes'
VCFs **whole, with nothing filtered out** — inside one test process both routes record the same
`##commandline`, so the exemption the script takes is not needed there.

## What the exemption is, and what it is not

`##commandline` records the command a person typed. `call-from-alignments` and `call-from-psps`
are different commands, so at a shell the line differs by construction rather than by anything
the calls decided. Dropping it is the same shape as spec §12.1's exemption of the psp header's
timestamp: a field a comparison deliberately skips, named, rather than a comparison quietly
weakened. Every other header line — the contigs with their checksums, the `INFO` and `FORMAT`
declarations, the sample columns in their order — and every record is compared byte for byte.

## The fixture, and why the existing one could not be used

**`a_cohort_on_disk`'s reference is all `A`s**, so the catalog routes every base to the
repeat-tract path and a run over it writes no VCF record at all. Two empty files are equal for
the wrong reason. The new `a_varying_cohort_on_disk` is built to discriminate, and **every one of
its four discriminating properties was added because a defect was measured surviving without
it**:

| the fixture carries | because without it, this passed |
|---|---|
| two samples varying at different positions | one sample's observations given to the other — every record's two columns were the same string |
| a deliberate repeat tract with a length variant | every stored locus written as `Generic`, discarding every motif and flank |
| the alternative reads leaning to one strand | the stored forward-strand read count zeroed on write |
| three read groups, and a run under parameters that score them differently | the walk-local-to-run-wide read-group renumbering deleted outright |

**The last one is not a fixture property alone.** Under `--defaults` every read group is scored
with the same numbers, so which group a read belongs to reaches no genotype and dropping the
renumbering changes nothing — measured, with three read groups already in place. The second test
gives the three groups base-quality multipliers of 0.25, 2.5 and 4.0, which puts read-group
identity into every likelihood; that test is also the only place psp mode's supplied-parameters
path is compared against direct mode's.

**Two numbers in the fixture were measured rather than chosen.** At a read stride of forty it
gives 1.4 reads a locus and the run writes no record at all — the cohort is analysed, the loci
are built, and two reads are not evidence of a variant. At five it is 18.2 reads a locus and
three records, at exactly the three places the fixture varies. And the reference generator's
freedom from accidental repeats holds by one base: its longest tandem stretch is 5.5 copies of a
period-2 motif against a six-copy floor, which is why `the_fixtures_ground_is_typed_as_its_doc_says`
asserts the segmentation — three regions, only the middle one a tract, at the coordinates the
fixture declares — rather than trusting the note.

## Where this oracle stops

- **It compares VCFs, so a stored locus that produces no record is not compared.** 578 of the
  fixture's 581 stored loci a sample carry no variant, and a psp route that dropped one of them
  passes. What holds a psp equal to the walk *field for field* is Milestone B's oracle
  (`examples/ng_psp_gather_oracle.rs`), not this one.
- **Neither route fits parameters**, so what only a fit reads — the stored sum of squared
  mapping qualities, and the count of reads that covered a locus without producing an
  observation — can be destroyed on write with both comparisons green. A property of a run that
  does not fit rather than of the fixture; the real-data script shares it.
- **Quantisation drift in the stored base qualities is invisible.** Zeroing the quality sum is
  caught; adding one quantisation step to it is not, because the step is below what a genotype
  quality resolves.

## Deviation from the plan, recorded

The plan asks for the comparison "on the run fixtures **and** on six tomato accessions over
400 kb". The tomato half is a script rather than a test, because the data it needs is 123 MB and
is not in the repository; the script is, and the report quotes what it printed. The fixture half
is what keeps the claim from silently stopping being true between such runs.

## How it is verified

`cargo test --lib` in the container: **6,217 passed, 0 failed, 15 ignored**. `cargo fmt --check`
and `cargo clippy --all-targets --all-features -- -D warnings` clean. The real-data oracle as
quoted above.
