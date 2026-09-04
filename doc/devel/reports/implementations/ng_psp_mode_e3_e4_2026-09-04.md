# ng psp mode — E3+E4: stored files call, and the reach ceiling has a reader

**Date:** 2026-09-04
**Plan steps:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone E, steps E3 and E4 — **bundled deliberately**: E4 is an accessor and two deferral notes discharged, and its whole point is that the calling stage E3 builds is the thing that now holds the number.
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §3.1, §5.3; [cohort_merge.md](../../ng/spec/cohort_merge.md) §13
**Branch:** `ng-psp-mode`

## Plan

E3: drive the lifted calling loop from a cohort of stored files. E4: the reach ceiling the psp
header has carried since A3 gets the reader `cohort_merge.md` §13 routed here.

## Changes made

- **`PspVariantCaller::call_cohort_handing_each_record_over`** — one source per open psp, each
  carrying its own map into the run's read-group numbering, into `ObservationCache::over`, into
  the body direct mode drives. Everything after the sources cannot tell it is reading files.
- **`PspVariantCaller::observation_reach_ceiling`**, and the two deferral notes discharged: the
  one in `observation_cache.rs`'s `cover` and `cohort_merge.md` §13's first bullet.

### What psp mode does *not* return, and why it is not an omission

Direct mode's `WrittenCohort` carries per-sample walk tallies beside the calling tallies. A
walker knows what its walk saw — the regions it handled, its generators' counts, its read
filters' drops; a psp source knows none of it, because the walk happened in another process and
what it counted is in that file's provenance rather than in its records. So this returns the
calling tallies alone, and what a run over stored files can say about each sample is the run
report's question at F1.

### The reach ceiling has a reader and still no consumer

§13 said a reader "takes the maximum over the cohort's files and sizes its cache accordingly".
The first half landed at E1 and is now exposed; **the second half is a no-op and the code now
says so** — the observation cache has no capacity to size, it grows on demand. The number is
there for the psp-mode performance work rather than for the merge, and both the spec bullet and
the cache's own note record that rather than leaving a reader to look for the sizing.

## Tests added

One, and `cargo test --lib 'ng::run'` goes from 498 to 499: **a cohort of stored samples calls
and hands its records over** — two samples, two variants each, records leaving in genome order,
the count and the records handed over pinned as one fact.

**A pairing assertion rather than a test**, at the one place two per-sample lists are walked
together: the psp's own sample name against the read-group entry's. A mis-pairing would hand one
sample's records another sample's read-group numbers — every number in range, scored against
another individual's calibration, nothing failing. It cannot fire on data; it fires on an edit.

### Mutation pass

Three, all killed (`tmp/e1_mutations/e3.sh`):

| mutation | killed by |
|---|---|
| the sources built with no read-group map | the calling test |
| the psps paired with the read-group tables in reverse | the pairing assertion |
| only the first sample becomes a source | the calling test |

## Validation results

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib 'ng::run'` — **499 passed**, 0 failed.

## Tradeoffs and follow-ups

- **The genotypes are not compared against direct mode's here.** That is spec §12.3's
  mode-equivalence oracle and it is F2's — this pins the join, not the arithmetic.
