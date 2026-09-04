# psp mode — F3: the remaining run-level invariances

**Date:** 2026-09-04
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone F, step F3
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §12.2, §12.6, §12.7, §12.9
**Branch:** `ng-psp-mode`

## The answer

**All four hold.** Three are tests in `pop_var_caller_exp::mode_equivalence`; the fourth is a
script, and the reason it cannot be a test is worth knowing.

| what | how it is held |
|---|---|
| §12.6 the order the psps are named in does not change the calls | `the_order_the_psps_are_named_in_does_not_change_the_calls` |
| §12.7 a cohort walked one sample at a time calls what one invocation calls | `a_cohort_walked_one_sample_at_a_time_calls_what_one_invocation_calls` |
| §12.9 ground analysed and found empty is not ground never looked at | `a_sample_that_analysed_ground_and_found_nothing_is_not_a_sample_that_never_looked` |
| §12.2 the VCF does not depend on the thread count | `scripts/ng_psp_concurrency_invariance.sh` — **599 records at 1, 2, 4 and 8 threads, identical apart from `##commandline`** |

## §12.6 — and it is compared sample for sample, not line for line

The sample columns follow the order the `--psp` flags were given, so the two files' *text*
differs by construction. What must match is each sample's own fields at each locus, which is what
the test compares — and it asserts the two files differ first, so a comparison that had somehow
been given one order twice would fail rather than pass.

**What it proves is that nothing joins by position.** Every per-sample number a run carries — the
inbreeding coefficient, the read-group calibration, the column a genotype is written into — is
keyed by name, and a run that joined any of them by argument order would give a reordered cohort
another sample's numbers with nothing failing. In the file that failure is invisible: every field
would be in range and every column would carry a genotype.

## §12.7 — compared against the one-invocation cohort, not merely run

Each sample is walked in its own `generate-psps` invocation, into a directory of its own — the
way this command's own help tells a person to spread a cohort — and the resulting VCF must equal
the one the same cohort gives when one invocation walks both.

**The fixture makes the failure reachable.** A gatherer sees one sample's files, so it numbers
that sample's read groups from zero; separately-walked samples therefore collide on `0`. This
cohort's first sample declares two read groups and its second one, so the two tables collide on
`0` *and* disagree about what `1` means — a run that failed to merge them into one run-wide
numbering could not produce this VCF.

## §12.9 — the assertion is on the file, because the VCF cannot carry it

A sample with no reads over a stretch and a sample that never looked at it produce the same
absence of records, which is the whole reason a psp records its analysed regions. So the test
asserts what is outside the VCF: the psp of a sample with no reads **holds no record and its
header still claims the whole ground**, and the run calls that cohort rather than refusing it or
dropping the sample — both samples keep their column.

## §12.2 — a script, and the flag is why

`--threads` builds rayon's **global** pool, which a process may build once. A unit test sweeping
thread counts would build the pool on its first call and silently run every later count at the
first one's width, reporting a sweep it did not do. One invocation a thread count is the only
honest measurement, so it is
`scripts/ng_psp_concurrency_invariance.sh <reference> <catalog> <out-dir> <psp>…`.

Measured on the six tomato accessions' psps: **599 records at every thread count, and the files
are byte-identical apart from `##commandline`** — which carries `--threads N` and therefore
differs by construction. That last point was found by running the script without the exemption
and reading the diff: the four files differed on that line alone, and the script's own comment
had claimed the thread count did not appear in the header. Both the comment and the comparison
are fixed.

## How it is verified

`cargo test --lib` in the container: **6,220 passed, 0 failed, 15 ignored**. `cargo test --tests`:
one failure, `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`, which is main's
and pre-dates this branch. `cargo fmt --check` and
`cargo clippy --all-targets --all-features -- -D warnings` clean.
