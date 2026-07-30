# Code Review: ng generic locus generator — the generator, Milestone A

**Date:** 2026-07-29
**Reviewer:** rust-code-review skill (orchestrator) — 9 category sub-agents
**Scope:** the Milestone A diff, `58b3db4..942d59f` (six commits, A0–A5)
**Status:** Request-changes → **fixes applied**, see §11

---

## 1. Scope

- **What was reviewed:** a branch diff — six commits on `ng-generic`, `da778ab` (A0) … `942d59f` (A5).
- **Reviewed against:** `942d59f`, worktree `/Users/jose/devel/pop_var_caller-ng-generic`.
- **In-scope files:** [copy_fidelity.rs](../../../../src/ng/locus_generation/pileup/copy_fidelity.rs),
  [errors.rs](../../../../src/ng/locus_generation/pileup/errors.rs),
  [genome_walk.rs](../../../../src/ng/locus_generation/pileup/genome_walk.rs),
  [mock_reference.rs](../../../../src/ng/locus_generation/pileup/mock_reference.rs),
  [mod.rs](../../../../src/ng/locus_generation/pileup/mod.rs),
  [open_record.rs](../../../../src/ng/locus_generation/pileup/open_record.rs),
  [parity.rs](../../../../src/ng/locus_generation/pileup/parity.rs), and the plan doc.
- **Deliberately out of scope:** `src/pileup/**` (production, frozen — the oracle, not the
  subject); `cigar_cursor.rs`, `decompose.rs`, `active_read_set.rs`, `chain_id_allocator.rs`,
  `tests.rs` (still verbatim copies — reviewed for *fidelity*, not content); Milestones B–D,
  which are not built and whose absence is not a finding.
- **Categories dispatched:** reliability, errors, naming, idiomatic, refactor_safety,
  module_structure, smells, defaults, extras. `unsafe_concurrency` skipped (no `unsafe`,
  threads, atomics or `async` in the diff); `tooling` skipped (no `Cargo.toml` change).

**A method flaw in this review, recorded because it affects how its evidence should be
read.** All nine agents were dispatched in parallel against **one shared worktree**, and
several of them mutate files as their checklist requires. They collided: agents report their
edits being overwritten mid-run, `src/` reverted under them three times, one build hitting a
truncated file, and one agent's baseline failing on another's marker. Five agents detected
this and re-ran their experiments in private detached worktrees; that is why their numbers
are trustworthy. **Both Blockers below were therefore re-verified serially by the
orchestrator**, at a clean tree, before being accepted — and both reproduced exactly. Future
milestone reviews should give each mutating agent its own worktree.

## 2. Verdict

**Request-changes** — two Blockers, both in the milestone's *own verification apparatus*
rather than in the walk. All fixes have since been applied and are verified in §11.

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib` | `2671 passed; 0 failed; 5 ignored` |
| `cargo test --all-targets --all-features` | green but for the pre-existing, unrelated panic at `benches/psp_writer_perf.rs:386` |
| `cargo doc --no-deps` | 13 unresolved intra-doc links — 12 pre-existing and out of scope, **1 in scope** (Mi1) |
| `cargo audit` | **not run** — `error: no such command: audit`; cargo-audit is not installed in the container. No dependency changes in this diff |

Findings labeled "Needs verification": **0**. Every finding below was either reproduced by
its agent in a clean worktree or re-verified by the orchestrator.

## 4. Open questions and assumptions

1. **Is `--max-record-span` allowed to exceed `u16::MAX`?** (B/M1.) The run encoding is
   `u16`; the flag is an unbounded `u32`. Either the flag gains a ceiling at config
   validation (C1's step) or `ReadCoverage` needs a wider run. Until then a footprint above
   65,535 silently reports a truncated `positions_covered`.
2. **Should the contributor skip in `refold_live_reads` be pinned, or stay recorded?** It is
   unpinned — mutating it away leaves the suite green — and given the carry correction its
   only residue is bucket-*creation* order. Affects M8 and the note now in the code.
3. **Is bucket-creation order part of ng's contract before B2's sort lands?** Determines
   whether (2) is a test or a comment.

## 5. Top 3 priorities

1. **B1** — the permanent anchor is blind to `placed_start`, the one field this milestone
   made fragile. It absorbed 2,542 records carrying a real wrong number while printing "same
   support totals".
2. **B2** — nothing asserts an *emitted* record has no unsupported allele bucket, and the
   differential is blind by construction, so A3's eviction could be moved before the fold
   loop with the whole suite green.
3. **M1** — `coverage_of` saturates silently at `u16::MAX` while the cap that bounds it is an
   unbounded CLI flag.

## 6. Findings

### Blocker

**B1: [parity.rs:1138](../../../../src/ng/locus_generation/pileup/parity.rs#L1138) — the
permanent anchor omits `placed_start` from its evidence check, and absorbs real wrong
numbers.**
**Categories:** refactor_safety, reliability (convergent).
`total_support` summed six of production's seven support scalars, read field-by-field rather
than destructured. `classify_record`'s `EvidenceIntact` verdict — the class the anchor
*tolerates* — therefore never compared `placed_start`, nor chain ids at all.

Injecting a genuine defect (`finalise` stops counting `placed_start` for the REF bucket) left
`ng_walks_identically_to_production_on_complete_reads` **green**, moving the tolerated class
from 264 records (0.15 %) to 2,806 (1.58 %) — 2,542 wrong records absorbed — while its own
census line went on printing "same support totals" and "Every other record is identical,
field for field". Only `placed_left_and_placed_start_are_per_record`, an **inherited** test
from production's suite, noticed. Verbatim, at a clean tree:

```
complete-reads differential: 177782 records compared over 1600 cases, 16001 of them
multi-base. 2806 (1.58%) hold the same evidence — same reference bytes, same support totals
— with some rows' bases differing …
test … ng_walks_identically_to_production_on_complete_reads ... ok
```

This is the test that **replaced** the retired stage-1 differential. `placed_start` is
precisely the field A1 made fragile: ng dropped it from its own stats and `finalise`
reconstructs it from a per-read flag.

*Fix:* make the sum a struct built by an **exhaustive destructure** of production's
`AlleleSupportStats`, so a field added to that type stops the file compiling; and check chain
ids too — by **subset**, not equality (see B1a).

**B1a — a discovery made while fixing B1, worth recording as its own result.** Requiring
chain-id *equality* fails, and correctly: at `seed 0x5eed0001 case 11, record 30` production
folds a read into its REF bucket (`num_obs: 2`, id dropped by the `allele_index == 0` rule)
having missed an insertion it never re-folded, while ng emits the nine bases the read
actually witnessed as its own row carrying chain id `6`. That is the defect being fixed
showing up in the ids. The invariant that *does* hold is directional — ng's REF bucket holds
the full reference bytes, so any partial or unseen-event witness lands outside it and keeps
its id; **ng's id set is a superset of production's**, and an id production has that ng lacks
would mean ng lost a read's identity.

**B2: [open_record.rs:1479](../../../../src/ng/locus_generation/pileup/open_record.rs#L1479)
— nothing asserts an emitted record has no unsupported allele bucket.**
**Categories:** reliability.
Moving `evict_unsupported_alleles` to *before* the contributor fold loop leaves the whole
suite green (`150 passed; 0 failed`). Buckets emptied by the widen are still caught, because
the widen runs first — but buckets emptied by the **fold loop itself**, when a contributor
re-folds into a different bucket at that position, survive to `finalise` and are emitted.

The gap is structural, not accidental, and that is what makes it a Blocker rather than a
missing test: A3's own eviction fixtures reach into `OpenPileupRecordTable` while the record
is still *open*, and the differential is blind **by construction**, because `comparable`
drops unsupported non-REF buckets on both sides — precisely so ng evicting them and
production keeping them is not read as a divergence. The one projection that makes the two
walkers comparable is also the one that hides whether ng evicted anything at all.

That accumulation is what spec §7 says to *design for rather than discover*, against a
`find_allele_index` that is a linear scan with a full byte compare run once per (record,
contributor) per position.

*Fix:* assert it on the records that leave the walker, with a vacuity guard that the fixture
actually produces records having somewhere to move a read to.

### Major

**M1: [open_record.rs:155](../../../../src/ng/locus_generation/pileup/open_record.rs#L155) —
`coverage_of` saturates silently at `u16::MAX`, and its stated reason for safety is the wrong
bound.** **Categories:** defaults, errors, reliability (convergent, three ways).
The doc argued saturation was unreachable "because production's `max_record_span` is 5000".
The bound is not that constant — it is `--max-record-span`, an unbounded `u32` CLI flag.
Probe output: a 70,000-position witness returns `positions_covered: 65535`; a run starting at
offset 69,999 returns `offset_in_locus: 65535`. A wrong number, no error.
*Fix applied:* the false claim corrected, a `debug_assert` stating the envelope, and the
enforcement assigned to **C1** — the step that turns these constants into
`PileupGeneratorConfig` — since a ceiling belongs with the knob, not with the reader.

**M2:
[open_record.rs:1220](../../../../src/ng/locus_generation/pileup/open_record.rs#L1220) —
`refold_live_reads` re-places a read by assigning two fields, so a new `FoldedReadState`
field is silently carried stale.** **Categories:** refactor_safety, idiomatic (convergent).
A new field errors at the fold literal in `fold_read_into_record` and **nowhere** at the
re-place. That asymmetry is the A3 defect re-armed: `witnessed` going stale across a widen is
exactly how this went wrong the first time, and the field's own doc records the consequence
as "a wrong depth with no error".
*Fix applied:* the state is rebuilt from an exhaustive destructure. Verified — a probe field
now produces **three** errors (fold literal, re-place destructure, re-place rebuild) where it
produced one.

**M3: [open_record.rs:635](../../../../src/ng/locus_generation/pileup/open_record.rs#L635) —
`widen`'s doc comment is still production's, describing behaviour A3 inverted.**
**Categories:** smells.
"Existing alleles are rewritten by appending the new reference bases" is exactly what A3
stopped doing. The body says the opposite in bold seventy lines below, a test pins the
behaviour the doc denies, and the doc never mentions the live-read re-fold that gave the
function two new parameters — on the one function whose changed semantics *are* Milestone A.

**M4: [parity.rs:1513](../../../../src/ng/locus_generation/pileup/parity.rs#L1513) — the
fabrication census asserts a floor but no ceiling.** **Categories:** reliability.
Its headline could be driven from 1.9 % to 91.8 % by a genuine defect and the test still
passed, reporting the defect *as the measurement*. A census whose only guard is "> 0" cannot
fail upwards, and this one's number is a deliverable (D3).

**M5: [open_record.rs:384](../../../../src/ng/locus_generation/pileup/open_record.rs#L384) —
`..RecordWitness::default()` in `finalise`.** **Categories:** defaults, refactor_safety
(convergent).
The one construction site in the diff that opted out of the rule stated twice elsewhere in
the same file. B2 explicitly reshapes `RecordWitness`; a field added then compiles here and
arrives as a silent `0`.

**M6: [open_record.rs:381](../../../../src/ng/locus_generation/pileup/open_record.rs#L381) —
`coverage_of` is never exercised on the real walk.** **Categories:** reliability.
`RecordWitness` is computed for every record and dropped by both walker call sites until B2
consumes it, so the function's only coverage was its own two hand fixtures.
*Fix applied:* a `debug_assert` in `finalise` that every folded read resolves to exactly one
coverage class — which runs on every record of every parity run, a few hundred thousand of
them.

**M7: [mod.rs:137](../../../../src/ng/locus_generation/pileup/mod.rs#L137) — ng's *fork* of
`DEFAULT_MAX_ACTIVE_READS` is indistinguishable by name from production's.**
**Categories:** defaults. **Recorded, not fixed** — C1's mandate is "production's `pub
const`s **by name**", and the module header already anticipates the two diverging. Harmless
today (pinned equal by a test); it becomes a trap the moment they differ.

**M8: [open_record.rs:1163](../../../../src/ng/locus_generation/pileup/open_record.rs#L1163)
— `refold_live_reads`' `ids.sort_unstable()` is unpinned, and unpinnable by the oracle.**
**Categories:** reliability. **Recorded, not fixed.** Determinism *was* verified empirically
by the extras agent — five processes with five different `AHash` seeds, with a canary printed
to prove the seed varied, gave byte-identical record digests on every fixture — but no test
fails if the sort is removed, because bucket-creation order is not compared. Related to open
question 3.

### Minor

- **Mi1** — `mod.rs:29` linked `[`parity`]`, a `#[cfg(test)]` module: the one in-scope
  `cargo doc` failure. *(module_structure, and independently the orchestrator.)* **Fixed.**
- **Mi2** — `coverage_of`'s clamp was one-sided: the left edge was pulled up to `record_pos`
  but never pushed down to `record_end_exclusive`, so an extent right of the footprint
  yielded an out-of-range `offset_in_locus` with `positions_covered: 0` — a run of no
  positions, which both `RefSpan` and `ReadCoverage` document as not existing. Unreachable
  from the fold; the *stated* reason for safety was not the real one. *(defaults, errors.)*
  **Fixed** — clamped at both ends, with a `debug_assert` on intersection.
- **Mi3** — the `usize::MAX` sentinel in `evict_unsupported_alleles` makes "evicted" a
  representable index, so a missed remap is a plausible subscript rather than a type error;
  `locus_generation/mod.rs:79` records this project hitting that trap in release twice.
  *(idiomatic.)* **Fixed** — `Vec<Option<usize>>` + `expect`, and the third counter dropped.
- **Mi4** — `mock_reference.rs:306` claimed "`PileupRecord` has no `PartialEq`"; it has a
  hand-written one, and `parity.rs` relies on it. Two new files in one milestone
  contradicting each other about one type. *(refactor_safety.)* **Fixed.**
- **Mi5** — `mod.rs:82`'s claim that "ng's non-test code imports neither
  `MultiChromRefFetcher` nor `ChromRefFetchError`" is false as scoped:
  `src/ng/raw_chrom_reader.rs:50` does, in public signatures. True of *this module* only.
  *(extras.)* **Fixed** — rescoped.
- **Mi6** — `copy_fidelity.rs:6-7`'s "no file is both [a copy and ng's own]" became false at
  A0, three lines above the release table that documents it. *(module_structure.)* **Fixed**
  — the third state ("released") is now named.
- **Mi7** — `src/ng/mod.rs:21-25` still called the whole module "a verbatim copy … kept
  provably identical". *(module_structure.)* **Fixed.**
- **Mi8** — two parity tests named for an identity neither asserts:
  `ng_walks_identically_to_production_on_complete_reads` counts `EvidenceIntact` instead of
  failing it *and* requires that class to be non-empty, so a truly identical run would fail
  the test named for identity; `..._on_real_reads` says in its own comment that the walkers
  differ on purpose. *(naming.)* **Fixed** — renamed to
  `ng_holds_the_same_evidence_as_production_on_complete_reads` and
  `ng_diverges_from_production_on_real_reads_only_where_a_read_did_not_witness`.
- **Mi9** — stale identifier names in doc comments: `apply_events_to_ref(_into)` for ng's
  `apply_events_into` at four ng-side sites (the three naming *production's* function are
  correct), `[ProductionView]` for `SharedReference`, and a `cargo test` example naming a
  deleted test. *(naming.)* **Fixed.**
- **Mi10** — `process_position` computed `rec_end` with a plain `+` where every other reader
  goes through `footprint_end_exclusive()`'s documented `saturating_add` (Mi8); a wrapped end
  yields an empty event window and a silently unfolded read. *(idiomatic.)* **Fixed.**
- **Mi11** — the eviction predicate is spelled twice, in `evict_unsupported_alleles` and in
  `parity::comparable_exact_q_sum`. If they drift, the differential misattributes ng's own
  eviction as a walker divergence. *(smells.)* **Recorded** — the two operate on different
  types (`OpenAllele` vs `AlleleObservation`), so sharing needs a trait or a macro; the risk
  is now stated in both places.
- **Mi12** — `open_record.rs:29`'s `use super::super::{…}` is the only file-level two-level
  path walk in `src/ng/`, in a block that is crate-absolute for everything else.
  *(module_structure.)* **Recorded.**
- **Mi13** — `coverage_of` takes two positional `u32` record coordinates; transposing them
  compiles. This is the hazard that made plan 1 introduce `LocusLen`. *(naming.)*
  **Partially addressed** — a `debug_assert` catches the transposition in debug builds; a
  footprint newtype is the real fix and is left as a decision (open question, §4).
- **Mi14** — the arch doc's *Module home* inventory is stale: `mock_reference.rs` missing,
  `mod.rs` described as holding the deleted shim, "eight of those are copies" now five.
  *(module_structure.)* **Recorded** — design docs are not this skill's to edit.
- **Mi15** — `find_or_create_allele_index` is `#[cfg(test)]` code kept alive by a test of
  itself while both hot-path fold sites hand-inline a borrowed copy of it. *(smells.)*
  **Recorded.**

### Nits

`RefSpan` is the third meaning of "ref span" in the crate (`var_calling::types::RefSpan`, and
`OpenPileupRecord::ref_span()` returning a *length*) — recorded rather than renamed, since
arch §1.2 names the type verbatim; `Q_SUM_GRAIN` holds the reciprocal of the grain;
`RecordFoldState` rhymes with the unrelated `FoldedReadState`; `affected: Vec::new()` is the
one per-position allocation on a struct built entirely of hoisted buffers; `process_position`
is 155 lines and four levels deep; `WalkerError::Fasta` concatenates its `#[source]` into
`Display`; `copy_fidelity`'s sanctioned-addition strip is unbounded.

## 7. Out of scope observations

- `benches/` cannot benchmark ng's walker at all — `pileup_walker_scaling.rs` drives
  production's `run` only — so this milestone's cost was unmeasurable with the repo's own
  tooling and the extras agent had to write a throwaway probe. Worth committing a bench
  before Milestone B adds more work here.
- `cargo doc --no-deps` has 12 pre-existing unresolved intra-doc links in the `ssr`, `em` and
  `sfs` modules.
- `benches/psp_writer_perf.rs:386` panics; pre-existing and unrelated.

## 8. Missing tests to add now

Added in §11: `ng_emits_no_allele_bucket_without_support` (B2), the census ceiling (M4), the
`finalise` coverage-class invariant (M6), and the chain-id subset check (B1a).

Still owed, and deliberately deferred:

- `every_folded_read_is_a_complete_witness_on_the_complete_reads_fixture` — the **D1**
  sharpening the code already asks for: replace `classify_record`'s frontier with the spec's
  own definition of the anchor class, which A4's `coverage_of` finally makes computable. Needs
  the witness tally to escape through the walker's public surface, which `tests.rs`'s verbatim
  guard still holds shut.
- A test pinning `refold_live_reads`' read-id sort (M8), if open question 3 says
  bucket-creation order is part of the contract.

## 9. What's good

- **The mutation discipline in the commit messages is real.** The reliability agent re-ran
  seven rows spanning A2–A5 and every one reproduced, including A5's `left: 4 / right: 1` to
  the digit; the extras agent independently corroborated four more claims. No table row was
  found false. The one row that says "caught by nothing" is honest and was verified honest.
- **`copy_fidelity.rs` demonstrably works, and proved it under adversarial conditions** — it
  caught a *different* agent's one-comment edit to a guarded file, from outside, with a
  message naming the original.
- **A0's error work is clean**: both lossy spots in the deleted adaptor are gone with nothing
  replacing them, and nothing was widened — `UnknownContig` now survives where it used to be
  folded into `Io { NotFound }`.
- **Spec §7's cost prediction was confirmed and contained by measurement rather than
  assertion**: longest allele list per record 8→13, 13→26, 18→41, but total alleles only
  +0.8–2.3 % and wall +2.3–6.4 %, not depth-driven.
- **Determinism is verified, not asserted** — five processes, five different `AHash` seeds,
  byte-identical digests, with a canary proving the seed varied.

## 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib
./scripts/dev.sh cargo doc --no-deps          # 12 pre-existing failures, none in scope
PVC_PARITY_CASES=5000 ./scripts/dev.sh cargo test --release --lib ng::locus_generation::pileup::parity
```

## 11. Fixes applied

All Blockers and Majors M1–M6 applied, plus Minors Mi1–Mi10. Each fix that claims to close a
hole was verified by re-running the injection that opened it:

| finding | verification |
|---|---|
| B1 | the same `placed_start` injection now **fails** the anchor; the tolerated class stays at 264 (0.15 %) unmutated |
| B1a | chain-id equality fails on a real record; the subset check passes, and the record is quoted in the code |
| B2 | moving the eviction before the fold loop now **fails** `ng_emits_no_allele_bucket_without_support`, naming the emitted bucket `("CCATAG")` with `num_obs: 0` |
| M2 | a probe field on `FoldedReadState` now produces **3** compile errors, up from 1 |
| M6 | the `#[expect(dead_code)]` on `read_group` **fired** when the exhaustive destructure began reading it — removed rather than downgraded, with B1's remaining debt stated |
| M1, Mi2 | `coverage_of`'s envelope and intersection are now `debug_assert`ed; enforcement assigned to C1 |

`cargo fmt --check` clean, `cargo clippy --all-targets --all-features -- -D warnings` clean,
`cargo test --lib` **2672 passed; 0 failed**, `cargo doc --no-deps` down to 12 failures with
none in scope.

M7, M8 and Minors Mi11–Mi15 are recorded rather than applied, each with its reason above.
