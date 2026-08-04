# Review — ng read filtering in stages, D1: the first filter runs with no reference

**Date:** 2026-08-04 · **Branch:** `ng-generic-perf` · **Base:** `f5630f8`
**Scope:** the working-tree diff for plan step **D1** — as submitted, 53 added lines in
`src/ng/read/filtering.rs` (one `#[cfg(test)] mod the_first_filter_needs_no_reference`).
**Design authority:** [spec](../../ng/spec/read_filtering_stages.md) §5, §8 ·
[arch](../../ng/arch/read_filtering_stages.md) §3.2 ·
[plan](../../ng/impl_plan/read_filtering_stages.md) D1, Checkpoint D.
**Fixes:** [`fixes_applied_2026-08-04_v1.md`](fixes_applied_2026-08-04_v1.md).

Three `general-purpose` agents, each in its own git worktree, each detached at `f5630f8` with the
staged diff applied by `git apply`. Raw findings kept under
`tmp/review_2026-08-04_d1-first-filter-no-reference/`:
`reliability.md`, `naming-smells.md`, `intent-refactor.md`.

Every agent was given the mutation mandate explicitly and told to contradict the author's
justification prose. All three did, and all three were right.

---

## 1. Verdict

**4 Major, 8 Minor, 8 nits. No Blocker.** The diff was tests-only and could not change compiled
behaviour, which is what kept it out of Blocker territory. Everything else about it was wrong.

The three reviews converge on one sentence, which the reliability agent put plainly:

> **Is D1 as built a test that can fail? Yes — but never alone, and never first.** I could not
> construct a single mutation that D1 catches and the tree at `f5630f8` does not. Its unique
> detection power, measured, is zero.

The step was **rebuilt rather than patched** — see the fixes report.

## 2. The four Majors

**M1 — the checkpoint mutation was already caught, twice, one of them by production.**
Adding `_reference: &impl RawRefSeq` to `verdict_on_raw_read` gives four `E0061` sites: the two in
the new module, `mod tests`'s `pre` helper, and `cursor.rs:706`. That last one is not a test, so
`cargo build` fails on the named mutation with D1 deleted. Reported independently by all three
agents, with identical output. *(Reliability also ran the generic-bound variant —
`fn verdict_on_raw_read<R: RawRefSeq + Default>` — and got the same picture in `E0283`.)*

**M2 — "its import list is the assertion" is false, and the property escapes through a name in
that very list.** `ReadFilterConfig` is imported by the module and is the route: the reliability
agent put a working `&'static dyn RawRefSeq` on the config behind `Default`, had the first filter
fetch bases from it on every call, and asserted inside the filter that it really read them. Result:
spec §5's property false, **2,867 passed, 0 failed**, D1 green. `ReadFilterConfig::default()`
launders exactly the requirement the module exists to forbid.

**M3 — the property is not pinned by scope, and a routine refactor un-pins it silently.** The
refactor-safety agent replaced the module's import list with `use super::*;` (which pulls in
`RawRefSeq`) and everything stayed green; moving the body verbatim into `mod tests` also stayed
green. A signature assertion, by contrast, cannot be repaired by an import.

**M4 — the alarm is not new, and the capability spec §5 claims is still unavailable.** At `7e8cfce`
(before C2) the function was already `verdict_pre_decode(flag, mapq, &config)` — reference-free, and
already called that way by `mod tests`. D1 transliterated would have passed before the change it
protects. And the reference bound did not disappear at C2; it moved onto the struct:
`AlignmentCursor<R: RawRefSeq>`, `AlignmentFile::cursor<R: RawRefSeq + ContigTable>`,
`SampleReads::cursor`. **Spec §5's three named callers still cannot filter a file's reads without
producing a reference.** What C2 bought is one layer down — `RegionRawAlignedReads` and
`AlignedReadsReader` carry no bound — so the capability is a reference-free *pass*, which the
submitted test did not exercise.

## 3. The Minors and nits, grouped

- **Three false factual claims in the prose**, each caught by two agents independently: the
  header-wide contig fetch went at **B1**, not C2 (spec §5's own amendment says so); "this module is
  the only thing that says so" is contradicted by `verdict_on_raw_read`'s doc at `filtering.rs:210`
  and by arch §3.2; and "one of the two things … (spec §8)" cites §8 for a count only the plan
  gives — §8 says three.
- **The test name described the signature, not the assertions** — a claim no `assert!` in the body
  could fail on.
- **The module name was an assertion sentence** where every other named test module in `src/ng`
  (`copy_fidelity`, `leftmost_property`, `scanner_parity`, …) is a noun phrase.
- **Both assertions were verbatim duplicates** of `pre_decode_keeps_a_clean_primary_read` and the
  first arm of `each_flag_bit_drops_to_its_own_bucket`, and the doc's defence of the duplication
  rested on M2/M3, which do not hold.
- **30 lines of `///` on a `cfg(test)` item**, which `cargo doc` never renders.
- **Mixed `super::` / `crate::` import paths** inside the list the comment called "the assertion".

## 4. What this review cost the author, and what it bought

The submitted step took the plan's D1 sentence literally — *"construct it and drive it without a
`RawRefSeq` in scope"* — and produced a module whose entire mechanism was its import list. Three
independent agents measured that mechanism at zero and one of them broke the property through it.

What replaced it has one measured, unique alarm: with a reference field added to
`ReadFilterConfig`, `Default` supplying it and `filtering.rs`'s `post_config` repaired the way its
author would repair it, the rebuilt file is the **sole remaining compile failure** in the tree. The
fixes report has the output.

**This is the fifth consecutive step in this plan whose review found a test that could not fail,
and the fourth where the author's own justification prose was overturned.** The pattern is now
well enough established to be worth naming at Checkpoint D.

## 5. Open, carried to Checkpoint D

- **Spec §5 overclaims.** "The reference stops being a precondition for filtering at all" is true of
  the filter, the reader and the narrowing, and false of `AlignmentCursor`. Either the spec says so,
  or the cursor gains a reference-free construction. Recorded in the new module's doc; no code or
  doc change made, because this skill does not edit design docs.
- **Plan D1's wording, "construct it"**, predates C2 making the first filter a free function.
- **The `const` witness is `#[cfg(test)]`.** Beside `verdict_on_raw_read` it would make `cargo build`
  carry the signature check; that is a line in a shipping file, which Milestone D's tests-only scope
  does not cover.
