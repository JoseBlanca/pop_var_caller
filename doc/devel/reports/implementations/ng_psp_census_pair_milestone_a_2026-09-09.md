# The census lives inside the psp — Milestone A: written into the trailer, and read back

**Date:** 2026-09-09
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestone A
**Spec:** [psp_census_pair.md](../../ng/spec/psp_census_pair.md) §3, §3.1, §5, §9, §11;
[psp_file_format.md](../../ng/spec/psp_file_format.md) §3.4
**Branch:** `census-vs-psp-perf`, on top of `9fdc427f`

---

## The tree this milestone is built on is not green, and each step is judged against it

Recorded once, on the untouched tree at `9fdc427f`, before the first line was written. Every
step below is compared against this and not against green.

| gate | on the untouched tree |
|---|---|
| `cargo fmt --check` | red in **4 files**: `examples/ng_call_cohort_end_to_end.rs`, `examples/ng_call_from_psps_cost.rs`, `examples/ng_census_locus_spans.rs`, `src/ng/psp/block.rs` |
| `cargo check --all-targets --keep-going` | red: **4 examples** do not compile — `ng_call_cohort_end_to_end` (1 error), `ng_candidate_selection_probe` (4), `ng_cohort_merge_parallel_cost` (8), `ng_cohort_merge_real_cost` (11) |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | red: **14 errors** in 6 files, none in a file this milestone touches except `psp_writer_line.rs:360`'s pre-existing `unused_must_use` |
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | red in exactly **1 test**: `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` (`tests/ng_calling_loop_calls_genotypes.rs:1241`). 6,672 lib tests pass; every other target passes |

**The four broken examples were not repaired, and three of them are why.** `ObservationCache::over`
has grown a `Box<dyn MergeReference>` argument, `MemberRange` has lost its `observations` field, and
`CohortObservation::over` now takes a `ClosedLocus` and a `WindowedCohort` where the probe holds a
`ClosedLocusRanges` — those are decisions about what the probes should measure, not signatures to
adapt. `ng_call_cohort_end_to_end`'s single closure-arity error is mechanical, but fixing one of
four does not make `--all-targets` green, so it was left as well.

**So the gate for each step is `--lib --bins --tests`, plus a `check --all-targets --keep-going`
that must fail on exactly those four examples and no more.** That second half is what catches a
change breaking a fifth.

---

## A1 — the walk hands the census to `finish`

**Committed:** see `git log` for `feat(ng): A1`.

### What it does

`PspWriterLine::finish` now takes the file's closing payload and passes it to
`PspWriter::finish`. `SampleObservationGatherer::write_psp` finishes its `CensusWriter` and
encodes the census **before** it seals the file, and hands the bytes across the writing thread's
queue. A walk opened without a census plan closes with nothing, which is what every psp in this
tree carried until now.

The census is written with **no pileup identity**. Spec §3: that identity exists to catch a census
paired with a psp it was not built from, and a census that *is* its psp's trailer cannot be paired
with anything else. The field was already an `Option`, so nothing about the format changed.

### The one deviation from the plan, and why

**The plan has A1 delete `write_psp`'s census-path argument; it is deleted in A2 instead, and for
this one step the census is written twice — into the trailer and into the file beside the psp.**

The plan's A1 cannot be committed green. Deleting the argument stops `generate-psps` writing
`<sample>.census`, and six tests in `generate_psps/tests.rs` plus three in `gatherer.rs` assert
that file exists — rewriting them is A2's own first sentence. So A1 keeps the sidecar untouched and
adds the trailer beside it; A2 removes the sidecar, the argument and those tests together. The
intent of both steps is unchanged and each is independently green.

**What the deviation costs while it stands:** the census is encoded twice per sample (the same
`SampleCensusEvidence`, two destinations), and a failure writing the sidecar now throws away a psp
that is already whole — `generate-psps` deletes the part-written psp on any error out of
`write_psp`, which was right while a psp without a census was useless. Both are recorded in the
code at the site, and both go with the sidecar in A2.

### What was measured

- **`cargo test --lib --bins --tests --all-features --no-fail-fast`: 6,677 lib tests pass against
  the baseline's 6,672** — the five this step adds — and the one pre-existing integration failure
  is the only failure, unchanged.
- **`cargo clippy --lib --bins --tests --all-features -- -D warnings`: the same 11 error locations
  in the same 6 files as the baseline.** No new one.
- **`cargo check --all-targets --keep-going`: the same 4 examples fail, with the same error
  counts** (1, 4, 8, 11).
- **`cargo fmt --check`: the same 4 files as the baseline**; the three this step touches are clean.

### The tests it adds, and what each would catch

Five, and two of them exist because the review asked what a wrong implementation would still pass.

- **`the_file_carries_the_trailer_the_line_was_asked_to_seal_with`** (`psp_writer_line.rs`) — the
  round trip across the thread seam: the payload handed to `finish` is the payload
  `PspReader::trailer` gives back. Nothing else in the pipeline notices an empty trailer, because
  a psp carrying one is a well-formed psp.
- **`a_trailer_of_megabytes_round_trips_through_the_line`** — **31 bytes was the largest trailer
  anything in this tree had ever written**; every other `finish` under test passes an empty slice
  or a short literal. A whole-genome census is a few megabytes of positions plus tens of megabytes
  of tracts (spec §3.1), so the size class this change moves the trailer into was the one class
  with no test at all. The payload is 4,000,000 varying bytes rather than a constant, so a defect
  that wrote one page twice would show.
- **`the_psps_trailer_is_the_census_the_walk_built`** (`gatherer.rs`) — the trailer decodes to the
  census the walk accumulated, with the identity absent.
- **`a_walk_with_no_census_plan_leaves_the_trailer_empty`** — the other half: no plan, no payload.
- **`a_walk_that_fails_with_a_census_plan_leaves_no_readable_psp`** — the existing failure test
  opens its gatherer without a plan, so the encoding step this change adds was not on any failing
  path. What must not happen is a psp sealed from a half-fed census: it would read back whole and
  a fit would score a cohort against it without complaint.

**One assertion in the third test was wrong on its first draft, and the mutation that found it was
run rather than reasoned about.** The draft asserted that some kept position was not marked
*never walked* — which sounds like "the walk fed the census" and is not. Building the writer marks
every kept position of the analysed ground as walked ([`gatherer.rs:294`](../../../../src/ng/run/gatherer.rs)),
so a census closed **before** the walk fed it comes back walked everywhere with a depth of zero.
Measured: moving `self.census.take().map(CensusWriter::finish)` above the walk loop leaves
`gatherer::census_tests` at **6 passed, 0 failed** with that assertion, while turning
`census_from_psp`'s two parity tests red. Rewritten as *some kept position has a depth above zero*,
the same mutation fails the test at `gatherer.rs`'s assertion. The mutation was reverted and the
module's tests re-run on the restored tree before anything was staged.

### What the review found, and what was done

Two reviewers over the working-tree diff, read-only, one on reliability and errors and one on
concurrency, smells and refactor safety. Four Major, five Minor and seven Nits between them; the
report is [psp_census_pair_milestone_a_2026-09-09.md](../reviews/psp_census_pair_milestone_a_2026-09-09.md).
Everything actionable was applied in the same step. The three carried forward — the sidecar's
paired-`Option` signature, `write_loci`'s `unreachable!`, and `RunError`'s boxed `dyn Error` — are
pre-existing shapes in code A2 deletes or does not touch.

Two of the fixes changed the design of the code rather than its prose:

- **`PspTrailer` replaces the bare `Vec<u8>`.** `line.finish(Vec::new())` and *the caller forgot
  the census* were the same value to the compiler, and the product of that mistake is a psp that
  opens, indexes and reports the right record count with no census in it — found at the parameters
  fit, a stage and a re-walk from its cause. Closing with nothing is now a sentence the author has
  to write: `PspTrailer::Nothing`.
- **`RunError::CensusNotEncoded` replaces `CensusNotWritten` for the in-memory encode.** The old
  variant rendered as *the census at …/zeta.psp could not be written*: the wrong operation, a path
  belonging to a different file, and at that moment the file there is a half-written psp. The new
  one names the sample and no path, because no file was touched.
