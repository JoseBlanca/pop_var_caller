# ng generic locus generator — the generator, Milestone A: the fabrication stops

**Date:** 2026-07-29 · **Plan:**
[locus_generation_pileup_generator.md](../../ng/impl_plan/locus_generation_pileup_generator.md)
steps A0–A5 · **Spec:** [locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md)
§4, §6, §7, §8 · **Arch:** [locus_generation_pileup.md](../../ng/arch/locus_generation_pileup.md) §1.2

Implementation report for Milestone A of plan 3 of 3. Six commits, one per step; A0 is a pure
refactor and A2, A3 each got their own commit because a mis-derived extent is a wrong depth with
no error, and `git bisect` has to be able to find which one.

## 1. Plan

Take ownership of the copied walker outright, then make it stop fabricating evidence.

Production's haplotype builder emits a reference base for every offset in a record's footprint that
no read event covered. At a six-base deletion locus a read that witnessed two bases is folded as a
full witness of a six-base haplotype it never saw, and `widen` extends that fabrication
**retroactively**, appending new reference bases to every allele bucket — including the buckets of
reads that have already expired. ng's rule replaces both mechanisms with one sentence: **nothing is
ever written into an observation that its read did not witness.**

## 2. Assumptions

One, and it was wrong: that spec §4's *"a live read re-folds against the wider window"* described
what the code would do. It does not — see §6.1. Everything else the plan asserted held.

## 3. Changes made

### A0 — the reference adaptor deleted, `copy_fidelity` narrowed (`da778ab`)

ng's walker fetches through ng's own `RefSeq`. `RefSeqFetcher`, `to_chrom_ref_fetch_error` and both
of that translation's lossy spots — a contig *name* rendered as an id, a `u64 → u32` narrowing —
are gone; `WalkerError::Fasta`'s source becomes `RefSeqError`; both fetch sites use `fetch_into`, so
`widen` writes into a buffer the open-record table owns rather than allocating a `Vec<u8>` per call.

**It went first because it is the last step the full differential can prove free**, and it did:
1,010,515 records over 20,000 release cases, zero divergences.

### A1 — the state the rule needs (`b47bd3b`)

`RefSpan { start, end }` (1-based inclusive, matching `PreparedRead`'s own coordinates);
`FoldedReadState` gains `witnessed` and `read_group`; ng's own `AlleleSupportStats`, production's
minus `placed_start`. Types only — `witnessed` is filled with the record's whole footprint, which is
exactly what production's fill assumes every folded read saw, so A1 asserts nothing new.

### A2 — the builder stops filling (`a760ab5`)

`apply_events_to_ref_into` becomes `apply_events_into`: it emits **only what the events cover** and
returns the extent, or `None` when the witnessed positions are non-contiguous.

The three traps the spec names were handled rather than discovered: `ref_seq` is still needed for an
indel's **anchor base** when no `Match` emitted it (the one recorded residual); `events_overlapping`
does not clip a deletion, so the extent is intersected with `[record_pos, record_end)`; and
`bases.len()` is **not** `positions_covered`.

### A3 — `widen` grows the REF bucket only (`36b2626`)

`alleles[0]` grows; no other bucket does. Every **live** folded read is re-placed against the wider
window, carrying its existing contribution; contributors are skipped; buckets no read is folded into
are evicted at the end of the fold.

### A4 — coverage resolved at `finalise` (`592f25d`)

`coverage_of(witnessed, record_pos, record_end_exclusive) → ReadCoverage`, called from `finalise` for
every folded read. A read's `witnessed` extent is absolute; what it *means* is relative to a
footprint that grows until the record closes, so a read that was a complete witness becomes
`Observed` when the record widens under it — **with nothing about the read having changed**, and with
no re-fold that would notice, since the read may have expired in between.

This is also the moment the shared type's own note anticipated: `ReadCoverage`'s
`from_left`/`from_right` cannot express a run flush with neither border, and the generic path mints
those ordinarily.

### A5 — the no-observation path (`942d59f`)

A read whose witnessed positions are non-contiguous yields no observation, and **which** reads those
were is recorded as a per-record set of ids. Both arrivals record: the fold loop, and
`refold_live_reads` for a read that folded contiguously and became non-contiguous when the record
widened across its interior gap.

## 4. The oracle, rebuilt before the old one was retired

The full stage-1 differential cannot survive a deliberate divergence, so plan 2's `parity.rs`
differential is **replaced, not deleted**:

- **`ng_walks_identically_to_production_on_complete_reads`** — the permanent anchor, over a generator
  whose reads span their contig end to end and are silent nowhere inside it. It asserts the
  multi-base record count, because a fixture of one-base records would hold green against any
  implementation at all.
- **`ng_diverges_from_production_only_where_a_read_did_not_witness`** — the census on the general
  fixture, asserting that A2–A5 change no record's existence, anchor, REF bytes, error item or
  `RunSummary`, and that **something** must differ.

Both compare through `comparable`, a projection applied to **both** sides: non-REF buckets with no
support dropped, the rest sorted by bases, `q_sum` rounded to 1e-9. That rounding is a **sixth**
divergence class where the spec names five, and it is named rather than absorbed — production keeps
an emptied bucket alive and accumulates `+q −q +q` into it where ng evicts and recreates it, so the
same read's sum starts from `0.0`. `float_only_divergences` counts how often it fires.

## 5. Validation

| | |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib` | **2672 passed**, 0 failed (2648 at the start of the plan) |
| `cargo test --lib ng::locus_generation::pileup` | 151 passed |
| `cargo test --all-targets --all-features` | green but for the pre-existing, unrelated panic at `benches/psp_writer_perf.rs:386` |
| `cargo doc --no-deps` | 12 unresolved intra-doc links, all pre-existing and out of scope. A thirteenth, in scope, was introduced at A0 and fixed by the review |
| `cargo audit` | **not run** — not installed in the container; no dependency changes in this milestone |
| host-native `cargo test --lib` | **2672 passed**, 0 failed — the same result outside the container |

**At soak scale, host-native** (`PVC_PARITY_CASES=5000 cargo test --release --lib
ng::locus_generation::pileup::parity`), all eight parity tests green:

| | |
|---|---|
| complete-reads anchor | **2,253,903** records over 20,000 cases, 197,380 of them multi-base. **3,073 (0.14 %)** in the one tolerated class; every other record identical field for field. 521 agree only after `q_sum` is rounded to 1e-9 |
| the eviction census (new, from the review) | **1,008,679** emitted records, 30,747 with more than one non-REF bucket; **none** carried an unsupported one |
| the fabrication census | **19,703 of 1,010,515 records (1.9 %)** carried bases production credited to a read that had not witnessed them — **D3's headline on synthetic data** |

> **`PVC_PARITY_CASES` does not reach the container.** `scripts/dev.sh` forwards only
> `CARGO_TARGET_DIR` and `HOME`, so `PVC_PARITY_CASES=5000 ./scripts/dev.sh cargo test …`
> silently runs the *default* 1,600 cases and finishes in under a second. Soak runs must go
> host-native, or the wrapper needs to forward the variable. Worth knowing before someone
> reports a soak they did not run.

**Mutation testing is the real validation here**, because every failure mode in this milestone is a
wrong number rather than a crash. Twenty mutations across the six steps, each applied, run, and
reverted. Three were found to be caught by **nothing** and had fixtures written for them:

| step | mutation | before | after |
|---|---|---|---|
| A2 | the `None` path skipping the subtract | nothing | a record widened across a read's interior `N` |
| A3 | `evict_unsupported_alleles` not called | nothing | both eviction tests |
| A3 | the eviction's remap neutered | nothing | the eviction fixture + 3 more |
| A3 | the `index == 0` REF guard removed | nothing | `eviction_keeps_the_ref_bucket_even_with_no_observations` |

The review then found three more that nothing caught, two of them in the oracle itself (§7):

| where | mutation | before | after |
|---|---|---|---|
| `parity.rs` | `finalise` stops counting `placed_start` for the REF bucket | the anchor **passed**, absorbing 2,542 wrong records | the anchor fails |
| `open_record.rs` | `evict_unsupported_alleles` moved before the fold loop | nothing | `ng_emits_no_allele_bucket_without_support` |
| `open_record.rs` | a field added to `FoldedReadState` | 1 compile error | 3 |

## 6. Deviations from the plan

### 6.1 `widen` re-places every live folded read — option (b), owner, 2026-07-29

Spec §4 asserts *"a live read re-folds against the wider window"*. **That is false** for a live read
with no event anchored at the widening position — the ordinary state of a read sitting inside its
own deletion, which is the ordinary state at exactly the long-deletion loci this port exists to fix.
Production's `process_position` re-folds only the *contributors* at the current walker position, so
such a read is live, is in the record, and never re-folds. Production hides it by appending the
reference bases to every bucket; with REF-only widening the read's `witnessed` extent would stay
pinned to the pre-widen footprint, and A4 resolves coverage from exactly that — **a wrong depth,
with no error**.

So `widen` takes the active set and re-places every live folded read. Two corrections landed on top,
both found by the differential:

1. **Contributors are skipped** — the fold loop re-folds them a moment later, and is the only place
   the mate-overlap decision is replayed.
2. **The re-placement carries the read's existing contribution** rather than recomputing it. `q_sum`
   encodes decisions the *walk* took at earlier positions and replayed into the fold; recomputing
   would silently undo every reconciliation made earlier in the record's footprint.

### 6.2 `placed_start` is reconstructed at the `PileupRecord` boundary, not dropped (A1)

The plan has A1 remove the field from ng's stats, which it does — but `finalise` still returns
production's `PileupRecord`, which *has* it, and writing a zero would have blinded the A2/A3 oracle
on every field of every record and broken an inherited test inside the still-verbatim `tests.rs`.
A per-read flag reconstructs the per-bucket count **exactly**. Both die at Milestone B.

### 6.3 A2 took A5's subtract half

`apply_events_into` returning `None` also takes the prior contribution off its bucket, which the plan
assigns to A5. In the safe direction, and not a corner: a read can fold contiguously and *become*
non-contiguous when the record widens across its interior gap, so a bare `continue` would strand a
live contribution for a read with no row. A5 delivered the other half — *which* reads.

### 6.4 `genome_walk.rs` released from `copy_fidelity` at A0

It carried the `F: MultiChromRefFetcher` bound, so it could not stay verbatim. The plan's
parenthetical names two files, but its own "leaving … five still guarded" enumeration already
excludes it, so the five that remain are exactly as planned.

### 6.5 `RecordFoldState` grouped three parameters (A5)

`fold_read_into_record` reached eight parameters, over clippy's threshold. The three record fields it
mutates always travel together and always come from the same destructure, so they were grouped rather
than the lint silenced.

### 6.6 The parity anchor's frontier is weaker than the spec's definition, deliberately

`classify_record` accepts "same reference bytes, same support totals, some rows' bases differing"
where the spec's own definition of the anchor class is "every folded read is `Complete`". A4's
`coverage_of` now makes that computable, but the predicate needs the witness tally to escape through
the walker's public surface, which `tests.rs`'s verbatim guard still holds shut. Tagged **D1** in the
code, which is where the plan puts the permanent anchor.

## 7. Review, and what it changed

Nine category agents over the milestone diff
([report](../reviews/ng_locus_generation_pileup_generator_a_2026-07-29.md)): **2 Blockers, 8
Majors, 15 Minors.** Both Blockers were re-verified serially by the orchestrator at a clean
tree before being accepted, because the fan-out had a method flaw (below).

**Both Blockers were in the milestone's own verification apparatus, not in the walk — the
sixth consecutive milestone on this branch where the seam is weaker than its documentation.**

- **The permanent anchor was blind to `placed_start`.** `total_support` summed six of
  production's seven support scalars, read field by field rather than destructured, so
  `classify_record`'s tolerated class never compared it — nor chain ids at all. Injecting a
  real defect (`finalise` stops counting `placed_start` for the REF bucket) left the anchor
  **green** while moving that class from 264 records (0.15 %) to 2,806 (1.58 %): **2,542
  wrong records absorbed**, with the census line still printing "same support totals" and
  "Every other record is identical, field for field". Only an *inherited* test from
  production's own suite noticed. And `placed_start` is exactly the field §6.2 above made
  fragile. It is now a struct built by an exhaustive destructure of production's type.
- **Nothing asserted that an *emitted* record has no unsupported allele bucket**, so A3's
  eviction could be moved to before the fold loop with the whole suite green. The gap was
  structural: A3's own fixtures reach into the table while the record is still open, and the
  differential is blind **by construction**, since `comparable` drops unsupported non-REF
  buckets on both sides — the one projection that makes the walkers comparable is the one
  that hides whether ng evicted anything.

**A result fell out of fixing the first.** Requiring chain-id *equality* fails, and
correctly: at `seed 0x5eed0001 case 11, record 30`, production folds a read into its REF
bucket — id dropped by the `allele_index == 0` rule — having missed an insertion it never
re-folded, while ng emits the nine bases the read actually witnessed, carrying chain id `6`.
That is the defect being fixed showing up in the ids. The invariant that holds is
directional: **ng's chain-id set is a superset of production's**, and that is what is now
asserted.

Six Majors applied: the `u16` saturation whose stated safety argument cited the wrong bound
(`--max-record-span` is an unbounded `u32` flag, not the 5000 constant); `refold_live_reads`
re-placing a read by assigning two fields, which is the A3 defect re-armed and now produces
three compile errors instead of one when a field is added; `widen`'s doc header still
describing production's append-to-every-bucket that A3 inverted; the fabrication census
asserting a floor but no ceiling (its headline could be driven to 91.8 % and the test still
passed); `..RecordWitness::default()`; and `coverage_of` never being exercised on the real
walk, now covered by an invariant that runs on every record of every parity run.

Ten Minors applied, most of them **documentation that had stopped being true** — the
recurring failure mode: the `cargo doc`-breaking link, the false "`PileupRecord` has no
`PartialEq`", the overstated import claim, `copy_fidelity`'s "no file is both", `src/ng/mod.rs`
still calling the directory a verbatim copy, and two test names promising an identity neither
asserts. M7, M8 and Mi11–Mi15 are recorded with reasons rather than applied.

**What the review confirmed rather than found**, and which is worth as much:

- Every mutation-table row re-run — seven spanning A2–A5 — **reproduced**, including A5's
  `left: 4 / right: 1` to the digit. No claim in the four tables was false.
- `copy_fidelity.rs` works under adversarial conditions: it caught a *different* agent's
  one-comment edit to a guarded file, from outside, naming the original.
- **Spec §7's cost prediction, measured for the first time.** Longest allele list per record
  8→13, 13→26, 18→41 — but total alleles only +0.8–2.3 % and wall time +2.3–6.4 %, *not*
  depth-driven (+2.7 % at 200×). Confirmed and contained. The largest single component is
  §6.2's transitional `vec![0; alleles.len()]`: deleting only that `Vec` moves the overhead
  to +0.7–3.4 % across four fixtures. Recorded rather than optimised, because B2 deletes the
  code outright.
- **Determinism verified, not asserted**: five processes, five different `AHash` seeds, a
  canary proving the seed varied, byte-identical digests on every fixture.
- A 24,000-case debug soak (3.9 M records) of the randomised differential: no `debug_assert`
  in A1–A5 is reachable on legal input.

**A method flaw in the review itself, recorded so it is not repeated.** All nine agents ran
in parallel against **one shared worktree**, and several mutate files by design. They
collided — edits overwritten mid-run, `src/` reverted under an agent three times, one build
hitting a truncated file. Five agents detected it and moved to private detached worktrees,
which is why their numbers stand; both Blockers were re-verified serially regardless. **Future
milestone reviews must give each mutating agent its own worktree.**

## 8. Checkpoint A — what is open

None of these blocks Milestone B.

1. ~~**Is `--max-record-span` allowed to exceed `u16::MAX`?**~~ **Resolved (owner,
   2026-07-29): cap the knob.** `PileupGeneratorConfig` will reject `max_record_span >
   u16::MAX`, written into the plan's **C1** step. The reasoning is that the cap constrains
   nothing real — a locus is at most ~100 bp of reference, and a 5,000 bp record is already
   unreachable with Illumina reads, so production's default is generous fifty-fold and this
   ceiling six hundred-fold; widening the run to `u32` would touch the shared locus type and
   the STR generator to buy a range no data can occupy. **This makes `max_record_span` the one
   knob where ng's constant is not simply production's**, against C1's "by name" default.
   A4's `debug_assert` stays as the invariant's statement; C1 is its enforcement.
2. **Is bucket-creation order part of ng's contract before B2's sort lands?** Two findings
   hang on it: `refold_live_reads`' read-id sort is unpinned by any test (determinism *was*
   verified empirically across five hash seeds), and the contributor skip in the same function
   is unpinned for the same reason — given the carry (§6.1), its only residue is that order.
   Both are recorded in the code. If the answer is "no, B2's sort is the contract", both stay
   comments.
3. **`coverage_of` takes two positional `u32` record coordinates** and transposing them
   compiles — the hazard that made plan 1 introduce `LocusLen`. A `debug_assert` catches it in
   debug builds; a footprint newtype is the real fix and is a small type addition, not an
   implementation choice.
4. **ng's fork of `DEFAULT_MAX_ACTIVE_READS` is indistinguishable by name from production's**,
   while C1's mandate is "production's `pub const`s **by name**". Harmless today, pinned equal
   by a test.
5. **Nothing in `benches/` can drive ng's walker** — this milestone's cost had to be measured
   with a throwaway probe. Worth committing a bench before Milestone B adds more work here.
6. **The arch doc's *Module home* inventory is stale** (`mock_reference.rs` missing, `mod.rs`
   described as holding the deleted shim, "eight of those are copies" now five). Not this
   skill's to edit.
7. **The four Checkpoint-A questions from plan 2 are still open** and still not blocking.
