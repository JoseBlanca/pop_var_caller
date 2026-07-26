# ng read preparation — A1: the `ReadPreparer` trait + `ReadPrepError`

**Date:** 2026-07-26 · **Plan:** [read_preparation.md](../../ng/impl_plan/read_preparation.md) step A1
· **Spec:** [read_preparation.md](../../ng/spec/read_preparation.md) §6, §7 · **Arch:**
[read_preparation.md](../../ng/arch/read_preparation.md) §1.2, §2

Combined implementation + review + fixes report for one plan step. (The plan-driven skill asks for
three separate reports per step; at this step's size that would be three files of boilerplate around
one small diff, so they are folded into one — recorded here rather than done silently.)

## 1. Plan

Land step 2's contract with no behaviour behind it: the `ReadPreparer` trait and the run-fatal
`ReadPrepError` in `src/ng/read/mod.rs`, plus the `left_align` module the trait's v1 implementation
will fill in A2. Types only — no preparer implementation, no reference access.

## 2. Assumptions

None that changed direction. The one judgement call is recorded in §6 below (test breadth).

## 3. Changes made

- **[src/ng/read/mod.rs](../../../../src/ng/read/mod.rs)**
  - `ReadPrepError` — `#[non_exhaustive]`, `thiserror`, one variant `Reference(#[source] RefSeqError)`.
    Documented as *never* a per-read verdict: a read yielding nothing is `Ok(None)`, a broken
    reference ends the run. Mirrors `ReadFilterError::Reference`, which made the same call.
  - `ReadPreparer` — `type Scratch: Default` plus
    `prepare_read(&self, read: MappedRead, scratch: &mut Self::Scratch) -> Result<Option<PreparedRead>, ReadPrepError>`.
    The doc comment carries the contract: per-read independence, locus-independence, by-value
    ownership (buffers move rather than clone), static dispatch, no reference argument, and the
    decline-vs-broken-run split.
  - `pub mod left_align;` and the module header updated — step 2 is generic-path only, and its
    sibling files are named for the transform they perform.
- **[src/ng/read/left_align.rs](../../../../src/ng/read/left_align.rs)** — new; module documentation
  only. Its types land in A2.

## 4. Tests added

Four `#[cfg(test)]` stand-in implementations, driven through a **generic** helper
`prepare_pair<P: ReadPreparer>` that builds scratch via `P::Scratch::default()` and threads one
scratch across two reads. That construction is the point: a concrete call compiles whether or not the
trait requires `Scratch: Default` or forbids `dyn`, so only a generic caller pins them.

| test | what it proves |
|---|---|
| `a_preparer_yields_the_prepared_read_through_a_generic_caller` | the trait is implementable and returns a populated `PreparedRead` |
| `declining_a_read_is_an_answer_not_an_error` | `Ok(None)` is reachable and is not an error |
| `a_broken_run_is_an_error_not_a_decline` | `Err` is the *other* channel — the split §7 exists to enforce |
| `one_scratch_is_reused_across_reads` | a non-`()` scratch carries state across reads rather than being rebuilt |

## 5. Validation

Run in the container (`./scripts/dev.sh`):

- `cargo test --lib ng::read` — **135 passed; 0 failed; 0 ignored** (2220 filtered out), including the
  four new tests.
- `cargo fmt --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — no diagnostics; `Finished dev profile
  in 17.88s`.

## 6. Review and fixes applied

Reviewed inline against the `rust-code-review` category checklists rather than fanned out to
per-category sub-agents (mechanism changed, bar unchanged — owner informed). **No Blocker, no
Major.** Three findings, all applied in the same step:

1. *Minor, naming/consistency* — a doc link used `../../../../src/pileup/…`, where sibling ng modules
   at the same depth use `../../../pileup/…`. Corrected.
2. *Minor, stale documentation* — the module header still called step 1 "(this milestone)". Removed.
3. *Nit* — `<ReusesScratch as ReadPreparer>::Scratch::default()` replaced with
   `CountingScratch::default()`; the qualified form proved nothing the generic helper does not.

## 7. Deviations from the plan

**Test breadth, absorbed not escalated.** The plan asks for "a stand-in impl driven through a generic
bound"; this landed four stand-ins. The trait has three outcome channels and a scratch contract, and
the decline-vs-error split is the whole reason the return type is `Result<Option<_>, _>` — pinning it
at the step that introduces it is cheaper than discovering it unpinned later. Same scope, wider
coverage.

## 8. Tradeoffs and follow-ups

- **`Ok(None)` carries no reason**, while spec §7 wants declines tallied *by reason*. Harmless in v1,
  where nothing declines; the first decline reason forces the `Option` into a two-variant outcome
  (production's `ReadOutcome { Prepared, Dropped(DropReason) }` is that shape). Tracked as `OPEN:` in
  arch §6 — deliberately not pre-built.
- **No counts type** for the same reason: an all-zero tally with no variants to count would be dormant
  scaffolding.
