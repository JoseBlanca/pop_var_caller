# C1 — the driver calls a repeat tract

**Date:** 2026-09-02. **Plan:** [`calling_loop_ssr.md`](../../ng/impl_plan/calling_loop_ssr.md)
Milestone C step C1. **Design:** [`spec/calling_loop_ssr.md`](../../ng/spec/calling_loop_ssr.md)
§3.2. **Modules:** [`src/ng/run/callers.rs`](../../../../src/ng/run/callers.rs) (the dispatch)
and [`src/ng/calling/evidence_shaping.rs`](../../../../src/ng/calling/evidence_shaping.rs)
(`SsrEvidenceScratch`).

---

## What landed

**A repeat tract is now called through its own model instead of being set aside.** Both drivers
— the one that collects loci and the one that hands records over as it finishes them — go
through `call_one_cohort_locus`, which branches on the observation's kind: the SNP/indel arm is
what it was, and the tract arm runs `select_ssr` → `shape_ssr_locus` → the same
`genotyper.call_locus`. The guard that set every tract aside is gone; what the set-aside count
now holds is bundles alone.

**Four ends, not two.** The drivers used to read an `Option`: a call, or nobody to call. A tract
adds two more outcomes that are not the same fact — a repeat cluster nothing builds a caller for,
and a tract the model cannot describe — so the call site returns a four-way `LocusOutcome` and
each driver counts the three that produce no record separately. The report's own partition is
C3's.

## The one thing the design did not anticipate: the merge does not hand over what the tract model reads

`shape_ssr_locus` takes the STR generator's own observation rows. **The merge does not carry
them**: it interns each distinct sequence into the locus's allele table and folds every sample's
reads onto `(allele, read group)` pairs, and the driver has the merged object and nothing else.
Spec §2 lists evidence shaping as *built*, and it is — what was missing is the join in front of
it.

`SsrEvidenceScratch::rebuild` is that join, and **it is exact because a repeat tract is one
record per sample.** The merge folds two of a sample's rows together only where two of its
*records* fall inside one cohort locus, which cannot happen here, so each merge row is exactly
one generator observation. The bases come from the allele table, the read group is the row's own
key, the read count is its `num_reads`, and the witness is which of the merge's two lists the row
was in — `supported` is a read that spanned the locus, `partials` one whose reads ran out inside
it, carrying the positions it did witness.

**Four counters do not survive on a partial and none of them is read at a tract.** The merge
keeps a partial's reads and its error mass and drops its forward-strand count, its two MAPQ
moments and its left-placed count. The repeat-tract model reads an observation's bases, its read
group and its read count and nothing else — the site-quality artifact correction, which is what
reads the other four, skips a tract by design. They are filled with zeros rather than invented,
and the type's own documentation says so, so a later consumer that starts reading them sees zero
rather than a plausible wrong number.

## A candidate carrying no whole repeat stops the locus, and its frequency is measured

Selection counts a candidate's whole repeats by flooring its length by the motif's period, so a
sequence shorter than one copy of the unit comes back as zero. The evidence's counts are
`NonZeroU32`, because the stutter ladder is written in whole repeats and such a candidate sits
below its bottom rung. `repeat_counts_the_tract_model_can_take` is the conversion, and returning
nothing is the refusal.

**How often that fires, measured on HG002's 50,000-region Tier set through ng's own catalog: 1 of
17,315 kept candidates at 30×, and none at all at 50× or 300×** (E2's runs). Refusing is right at
one in ten thousand; if it stops being one in ten thousand the change to make is to the evidence's
count type, not to this refusal, and both the helper's documentation and the outcome's say so.

## Tests — 4 new, 2 reversed

| test | what it pins |
|---|---|
| `a_merge_row_is_rebuilt_as_the_observation_the_generator_minted` | the join: each row's own allele bases, its own read group, its own reads, one list per **run** sample |
| `a_partial_row_is_rebuilt_as_a_partial_with_its_witnessed_positions` | the witness survives — a partial scored as complete is scored as a *short* allele |
| `one_scratch_over_two_loci_holds_only_the_seconds_rows` | the buffers are emptied between loci, so no tract is scored against the last one's reads |
| `a_candidate_with_no_whole_repeat_is_not_convertible_for_the_tract_model` | the refusal, and that an empty list is not one |

**Two tests of the parallel plan's assert the opposite of what they did, and the reversal is this
step.** `a_tract_a_sample_varies_at_is_built_and_set_aside_uncalled` is now
`..._and_called_through_the_tract_path`, and the thread-invariance sweep's non-vacuity anchor
moved from the set-aside count — zero at a tract now — to a record over the tract's own ground.

**The reversed driver test says *which model* called the tract, not merely that something did**,
and that assertion exists because a mutation predicted to survive would have. The SNP/indel path
would also produce a record over that ground — it is what the replaced guard existed to prevent —
so the test reads the kind selection stamped on the candidate table.

## Mutation testing — five run, five killed, one hole recorded

Predicted before running, which is what found the missing assertion above.

| mutation | outcome |
|---|---|
| the rebuild takes the reference's bases instead of the row's allele | killed |
| the buffers are not cleared between loci | killed |
| a partial is rebuilt as complete | killed |
| a tract is dispatched to the SNP/indel arm | killed — **by the assertion added after predicting it would survive** |
| a zero repeat count becomes one instead of refusing | killed |

**The hole: nothing exercises the bundle arm**, because nothing in the run builds a
`LocusKind::SsrBundle` — the bundle generator is deferred and has no spec. Sending bundles to
any other arm cannot be caught by a fixture today.

## Validation

`cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib` 6,008 passed / 0 failed / 14 ignored, all in the container.
