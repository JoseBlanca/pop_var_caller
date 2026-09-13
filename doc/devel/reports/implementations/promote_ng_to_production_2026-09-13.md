# Promoting ng to the project's caller — implementation report

**Date:** 2026-09-13. **Branch:** `promote-ng`, from `0ba40e79` (the baseline) on. **Plan:**
[promote_ng_to_production.md](../../implementation_plans/promote_ng_to_production.md), whose §5
records every step, its deviation from the plan if any, and what was given up.

## What changed

The project had two callers: the older one — "production" in the code's comments — and ng, written
beside it. **Production is deleted, and ng is the only caller.** ng's modules moved from `src/ng/`
to the crate root, and its binary, `pop_var_caller_exp`, took the name `pop_var_caller`. The
subcommands keep their names.

It happened in five milestones, each ending at a checkpoint the owner reviewed:

| milestone | what it did | what proved it |
|---|---|---|
| A — baseline | recorded what ng computes before anything moved (below) | two runs of one binary gave identical digests |
| B — sever | ng stopped importing anything production owned; what it used was copied in | digests unchanged; tests unchanged but for the copies' own |
| C — freeze | every test that ran production as an oracle compares against production's recorded answers, runs a test-only copy of the function, or was deleted with its reason | digests unchanged; a sweep for production paths over ng's code found nothing |
| D — delete | production's library, binary, tests, benches, examples, scripts and benchmark drivers, its lint exemptions and unused dependencies | digests unchanged; every cargo gate green |
| E — move | modules to the crate root, the binary's rename, CI's test filter | digests unchanged but for `##commandline`, which the comparison excludes |

## The identity check

At every step, one script
([promote_ng_oracle.sh](../../../../scripts/promote_ng_oracle.sh)) called four tomato accessions
(`SRR5079906`, `SRR5079864`, `SRR5079878`, `SRR5079876`) over 20 regions of SL4.00 both ways — from
their CRAMs in one process, and through per-sample psp files — and hashed five outputs. All five
were identical to the baseline after every step, including the last (Milestone E's run, with the
renamed binary):

| output | md5 |
|---|---|
| VCF from alignments, `##commandline` removed (7,787 records) | `b74f3edf9026fa756621e2492baf019a` |
| VCF from psp files, `##commandline` removed | `b74f3edf9026fa756621e2492baf019a` |
| parameters file the run fitted | `5e6d874da147f898eda134ce342424b8` |
| window-coverage rows, sorted | `25d388c83819321a023deb57663c2208` |
| per-sample window histograms | `259baf040610b33f7ead595da7b26a3b` |

The psp files themselves are not hashed, because a psp header carries the time it was written;
everything a psp carries that reaches a call shows in the VCF from psp files. What this check covers
is one corner — four samples at tomato's low depth — and it is a check that nothing *moved*, not
that the calls are right.

## Sizes

| | before (`0ba40e79`'s parent) | after |
|---|---:|---:|
| Rust lines under `src/` | 409,119 | 314,430 |
| `cargo test --lib --tests` | 6,370: 6,358 passed, 1 failed, 11 ignored | 4,791: 4,787 passed, 0 failed, 4 ignored |
| binaries | 2 | 1 |
| examples | 86 | 62 |
| criterion benches | 12 | 6 |

The tests, reconciled from each step's commit message: Milestone B added 52, 43 of them the copied
code's own tests (B2's 9, B4's 10, B6's 24); Milestone C removed 63 net (oracle harnesses with no ng half left, and the copy-fidelity guards
that compared source text); Milestone D deleted 1,568 with production (50 integration tests, 1,517
library tests of which 4 were ignored, and one test of dead `bam` code). The one failing test at the
baseline was main's, fixed at B8.

Across the branch, 563 files changed, 74,573 lines inserted and 136,341 deleted; most of the
insertions are the recorded fixtures and the moved files' new paths.

## Deviations the plan records

- **C grew from 13 steps to 26** when a sweep for production paths found oracles outside the files
  named "parity"; three further findings in it were shipped code, not tests.
- **Where recording was impossible or wasteful, a copy replaced a fixture** (C12, C22, C24): a
  property test draws inputs a recording cannot answer, and one function was already copied.
- **D was reordered.** Production's library modules import each other in a loop, and each partial
  deletion left dead code that `clippy -D warnings` rejects, so the library went in one commit, after
  its consumers outside the library.
- **One piece of shared code was pruned**, a per-worker read source only production's STR caller
  used, because the compiler refused it; two more `bam` modules are left callerless (below).
- **E3 landed before E2**, since CI's filter matched no test once E1 landed. E3 also broke the CI
  file's YAML for one commit, fixed forward.

## Handed on

From the plan's §8, each with its reason there:

- a check of the aligner's band margin that needs no other caller — C10 lost the long randomised
  run against production that validated it;
- why a cohort of confident homozygous-variant samples has its site quality corrected from 900 to 0
  by a step documented as skipping it (found at C7; ng reproduces it alone);
- three `#[should_panic]` tests that expect a debug-only panic and fail under `--release`, which keep
  CI's release-mode step scoped to `calling::`;
- pruning `bam`, `fasta` and `regions` of items nothing calls — the largest are `bam::segment_reader`
  and `bam::segment_merge`, hidden from the dead-code lint by an allow in their own files;
- moving `doc/devel/ng/` up and repairing its paths, and renaming the `ng_` prefix off examples,
  tests, benches and scripts.
