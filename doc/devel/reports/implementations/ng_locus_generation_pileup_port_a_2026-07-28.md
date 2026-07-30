# ng generic locus generator — the port, Milestone A: the copy

**Date:** 2026-07-28 · **Plan:**
[locus_generation_pileup_port.md](../../ng/impl_plan/locus_generation_pileup_port.md) steps A1–A4 ·
**Spec:** [locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md) §3, §6, §12 ·
**Arch:** [locus_generation_pileup.md](../../ng/arch/locus_generation_pileup.md) *Module home*, §1.3

Implementation report for Milestone A of plan 2 of 3. One report for the milestone, four commits.

## 1. Plan

Copy production's pileup walker into `src/ng/`, **changing nothing about what it computes**, so that
every later change is a measured delta against a baseline that is provably production. The copy is
the mechanism, not a fallback: production is frozen, and ng needs its own read type anyway (§6), and
every module that looked reusable names `PreparedRead` in its signatures — so the read type reaches
all of them whatever else is decided.

The milestone's gate is spec §12: **production's own test suite, unmodified, green against the
copy.**

## 2. Assumptions

None that changed direction. Four judgement calls are recorded in §6; the largest is what "the
copied suite" in A4 means, which §5 answers with numbers rather than with the plan's wording.

## 3. Changes made

### A1 — ng's `PreparedRead` (`9bfd483`)

- **[src/ng/read/prepared_read.rs](../../../../src/ng/read/prepared_read.rs)** (new) —
  `PreparedRead`, `MateRole` and `ReadLengthError`, transcribed from
  [pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs) and extended with
  `read_group: ReadGroupId`. Deliberately **not** `#[non_exhaustive]`, for production's own stated
  reason, which this step then immediately collected on (§6.2).
- **`PreparedRead::from_production`** — production's prepared read plus the group it had nowhere to
  put. A destructure, so a field added to production's type stops it compiling rather than being
  dropped in silence; the mirror of `AlignedRead::into_mapped_read` on the way out.
- **[read/mod.rs](../../../../src/ng/read/mod.rs)** — `ReadPreparer::prepare_read` returns ng's type.
- **[read/left_align.rs](../../../../src/ng/read/left_align.rs)** — `into_prepared` reads the group
  off the `AlignedRead` before the conversion drops it, and re-attaches it. **The per-field wiring
  stays production's**: the read is still built by `prepare_passthrough`, its `--no-baq` arm, which
  is what keeps the read-preparation parity fixture exact on every field this step does not compute.

### A2 — the seven files, verbatim (`8b6307d`)

5,495 source lines into [src/ng/locus_generation/pileup/](../../../../src/ng/locus_generation/pileup/)
— 5,503 landed, the difference being the seam lines below — with `driver.rs` landing as
`genome_walk.rs`. Still emitting production's `PileupRecord`.

The copies reach their shared vocabulary through `super::`, and **leaving those paths alone is what
keeps the transcription verbatim**. ng's new `pileup/mod.rs` therefore stands in for production's
`walker/mod.rs`, drawing each name from wherever it now lives: `PreparedRead` / `MateRole` /
`ReadLengthError` from ng's read module, and `CigarOp` / `WalkerConfig` / the `DEFAULT_*` constants
from production, **by name rather than by literal**, so there is one source of truth until ng
deliberately diverges.

**The complete seam-edit record is `pileup/copy_fidelity.rs`**, which compares all eight copies
against their originals textually, in the tree, on every `cargo test`. That is the authority; the
prose below is a summary of it, and an earlier draft of this section claimed to be a closed list and
was not — it omitted the four `use crate::ng::types::ReadGroupId;` imports those fixture lines
require, `open_record.rs`'s `MockFasta` repoint, `tests.rs`'s doc block, and the import reorder
described in §5. The edits fall in three classes, and the test sanctions exactly those three:

- **two module-path repoints** — `cigar_cursor.rs`'s test module reaches ng's `decompose` oracle
  rather than production's (which is in a private module and does not resolve from here), and
  `open_record.rs` takes `MockFasta` from ng's copied suite. Matched *by presence, not by position*,
  because repointing a path changes where rustfmt sorts it.
- **28 `read_group: ReadGroupId(0),` lines** in `#[cfg(test)]` fixture literals, plus the five
  imports they require.
- **one module doc note**, at the head of `tests.rs`. ng may only *append* to production's header;
  the fidelity test asserts production's own header survives verbatim.

### A3 — the reference shim (`6e44051`)

`RefSeqFetcher<R: RefSeq>` in `pileup/mod.rs`. Semantically empty by design — the walker asks for
"uppercase ASCII over `{A,C,G,T,N}`, canonicalised by the fetcher implementation" and
`RefSeq::fetch_into` promises exactly that — and the claim is **checked, not asserted**:
`the_shim_canonicalises_what_it_serves` drives a contig written `acgtRYacgt` through it and requires
`ACGTNNACGT` back. Read preparation's uppercase divergence (`read_preparation.md` §6) does not recur
here.

A newtype rather than a blanket `impl<R: RefSeq> MultiChromRefFetcher for R`: most `RefSeq`
implementations will never play this role.

The error translation is an exhaustive match, so a variant added to either enum stops it compiling
instead of falling into a catch-all that reports the wrong failure mode. `UnknownContig` maps to
`Io`, which is where production's own trait documentation puts it. Two things are lost, both message
text only: `ChromRefFetchError` names contigs where `RefSeqError` numbers them — and `RefSeq` alone
cannot resolve one to the other, the contig table being a separate capability by design — and the
`u64 → u32` narrowing in the out-of-bounds arm, on a window that has already been rejected. Nothing
computes on either; the walker's `WalkerError::Fasta { chrom_id, .. }` carries the id anyway.

### A4 — the copied suite is green (this commit)

[pileup/tests.rs](../../../../src/ng/locus_generation/pileup/tests.rs) — production's
`walker/tests.rs`, copied. `pub(crate)`, exactly as production's is, so `MockFasta` / `snp_read` /
`paired_snp_reads` are reachable from B1's parity harness. `open_record.rs`'s import of `MockFasta`
is repointed from production's test module to this one — the fixtures are identical today and will
not be after plan 3.

**Strip the 24 `read_group:` lines, the import they require and the module doc note, and this file
is identical to production's** — checked by `pileup/copy_fidelity.rs` on every `cargo test`, not by
a one-shot `diff` into a scratch directory. (The hand-run `diff` quoted in an earlier draft of this
section left *two* hunks, not one: the doc block and the `ReadGroupId` import.)

## 4. Tests added

The milestone's own tests are few by design — the point of the step is the *inherited* suite. What
was added covers only the two pieces that are ng's rather than production's.

| test | what it proves |
|---|---|
| `the_conversion_from_productions_read_moves_every_field` | the destructure makes `from_production` exhaustive; only this makes it *correct*. Every field gets a distinct value, because `chrom_id` and `alignment_start` are both `u32` and swapping them would compile. |
| `every_production_mate_role_maps_to_its_counterpart` | a conversion that collapsed two roles would silently disable the mate-overlap tie-break rather than fail. |
| `the_read_group_rides_through_both_paths` (left_align.rs) | the whole point of ng owning the read type. Driven on **both** arms — the indel arm rewrites the CIGAR first — with a non-zero id, since `ReadGroupId(0)` is also what a defaulted field reads as. |
| the six `length()` / `MateRole` tests | the transcription of the one function in `prepared_read.rs` that computes anything, including that the `seq`/`bq_baq` check runs *before* the CIGAR check (the walker's message differs between them). |
| `the_shim_canonicalises_what_it_serves` + four siblings | §3's A3 — the shim's claim, its 1-based contract, and each error arm landing in its own variant. |

## 5. Validation — and what "the copied suite" actually is

**The plan's "46 tests in `walker/tests.rs`" is wrong, and the number that matters is three
numbers.** Established by counting rather than by reading:

| set | count | where it runs |
|---|---|---|
| `walker/tests.rs`, the end-to-end suite | **44** | copied to `pileup/tests.rs`; **name for name identical** to production's, checked by diffing `cargo test -- --list` output |
| inline `#[cfg(test)]` modules in the seven copied files | **70 `#[test]` markers, 69 in any one profile** | `subtract_contribution`'s debug and release tests are mutually exclusive by `cfg(debug_assertions)`; **both sides were run** — `cargo test --lib` and `cargo test --release --lib` are each 118 |
| ng's own additions (A1's 8, A3's 5, left_align's 1) | 14 | — |

So Milestone A's gate reads: **113 inherited tests (44 + 69) pass unmodified against the copy in a
debug build, 114 counting the release-only sibling**, and the 44-test end-to-end suite is identical
to production's name for name.

Run in the container (`./scripts/dev.sh`):

- `cargo fmt --all -- --check` — exit 0 **at `8fdba95` and at the milestone's head**. It was *not*
  green at `8b6307d` (A2), which landed `cigar_cursor.rs` with its test-module imports out of
  rustfmt order; A3's `cargo fmt` corrected it, and A3's commit message does not mention touching a
  copied file. Found by the review, verified by checking A2's version of that file back into the
  tree and re-running. The reorder is now a recorded, tested divergence (§3), and it is what
  `copy_fidelity.rs` caught on its very first run.
- `cargo clippy --all-targets --all-features -- -D warnings` — no diagnostics. *(Which is why
  `driver.rs:647`'s `#[allow(clippy::ptr_arg)]` had to survive the copy; it did, being inside the
  verbatim text.)*
- `cargo test --all-targets --all-features` — **2641 passed / 0 failed / 4 ignored** after the
  review fixes (2504 → 2513 → 2582 → 2587 → 2631 → 2641).
- `cargo test --release --lib ng::locus_generation::pileup` — 118 passed at A4, covering the
  release-only half of the `cfg` pair.

Two standard commands are excepted by hand, red independently of this work and tracked under
PROJECT_STATUS *Standing project-wide items*: the `--all-targets` run panics in
`benches/psp_writer_perf.rs:386` (the test totals quoted are from that same run, whose every other
target is green), and `cargo doc --no-deps` — **not re-run here**, known-red on 11 unresolved
intra-doc links.

## 6. Deviations from the plan

Four, all minor, none reaching past the step.

1. **`cigar_cursor.rs`'s test module reaches ng's `decompose`, not production's.** The plan permits
   module-path edits; this one is also forced — `pileup::walker::decompose` is a private module and
   does not resolve from `src/ng/`. It is the right target anyway: `decompose` is the oracle the
   cursor is parity-tested against, and testing ng's cursor against production's oracle would leave
   ng's own copy of the oracle untested.
2. **Four `#[cfg(test)]` fixture builders in the copies, and 24 in `tests.rs`, name `read_group`.**
   The plan's "only the module paths, and `PreparedRead` resolving to ng's" does not literally cover
   a struct literal gaining a field. It is forced by the type rather than chosen — and it is
   *exactly* what not making `PreparedRead` `#[non_exhaustive]` was for (A1's own instruction: "a new
   field should break every construction site"). Every one is `ReadGroupId(0)`; no fixture's
   behaviour moves.
3. **`walker/tests.rs` is copied, which the arch doc's file inventory does not list.** A4 says "all
   46 inherited tests pass unmodified" and the plan's verification table says "green **against the
   copy**", which is only meaningful if the suite runs against ng's walker; production's own tests
   exercise production's code and would prove nothing about the transcription. The arch inventory
   names the seven copies, `mod.rs` and `parity.rs`, and is silent on the suite. Copying it is
   additive and is what makes the gate a gate.
4. **`mod ref_seq_fetcher_tests`, not `mod tests`, for A3's tests.** The copied suite takes the
   `tests` name, so the two suites stay comparable module for module. *(Landed as `shim_tests` and
   renamed to name its subject when the review pointed at the crate's `mod baq_tests` precedent.)*
5. **`copy_fidelity.rs` is a ninth file the arch inventory does not list either.** Added after the
   review: it is ng's own, not a copy, and it is what turns "verbatim" from a claim into a test.
   It cannot live inside `tests.rs`, because `tests.rs` is one of the files it checks.

## 7. Review, and what it changed

Reviewed the same day over the whole milestone diff:
[ng_locus_generation_pileup_port_a_2026-07-28.md](../reviews/ng_locus_generation_pileup_port_a_2026-07-28.md)
— 9 categories, **0 Blockers, 3 Majors, 19 Minors**, verdict Approve-with-changes. The copied files
were scoped out for *content* and checked for *fidelity*; four categories confirmed independently
that no unsanctioned edit had reached them.

**All three Majors were the seams being weaker than their own documentation**, which is this
branch's recurring pattern:

- **M1 — `assert_same_prepared_read` claimed a field addition would force it to be updated, and did
  not.** Twelve independent `assert_eq!`s, no destructure, no `PartialEq` behind them — while
  `from_production` *does* destructure, so a new field would have been carried into ng's type and
  then silently escaped the thing the module calls "the port anchor". Both sides now destructure
  exhaustively with no `..`, and `read_group` is ignored **by name**.
- **M2 — the shim was never driven through the walker it exists to feed.** Every test called
  `RefSeqFetcher::fetch` directly, so the 113 inherited tests proved the walk only over `MockFasta`.
  `a_walk_over_the_shim_matches_the_same_walk_over_the_mock_reference` now runs the same reads
  through both fetchers, with a deletion so `widen`'s `fetch(chrom_id, old_end, extra_len)` — an
  *exclusive* end passed as a *1-based* start — is on the path. **Mutation-verified:** a `+ 1` on
  the shim's start fails 6 of 10 tests, serving one base short fails 5.
- **M3 — the verbatim property had no check that survived the branch.** The only evidence was a
  `diff` into a gitignored scratch directory. `pileup/copy_fidelity.rs` now compares all eight pairs
  textually on every `cargo test`, sanctioning exactly three classes of divergence and failing with
  a message that names the production original and says the file must not be edited.
  **It found a real divergence on its first run** — A2's `cigar_cursor.rs` import order, silently
  corrected by A3's `cargo fmt` and recorded nowhere (§5). **Mutation-verified:** five mutations,
  five failures — a reworded comment, a reworded doc comment, a rewritten repoint, two reordered
  production lines, and a deleted blank line. A sixth mutation *passed* and correctly so: it moved
  an ng-*added* line, whose position is not part of the copy.

Sixteen more findings applied, the substantive ones being: the `super::` vocabulary demoted from
`pub` to `pub(crate)` (it was minting public ng-flavoured aliases for frozen production types);
`RefSeqFetcher` narrowed to `pub(crate)` with a `pub(crate)` field, matching arch §1.3, its inert
`Copy` dropped and its `R: RefSeq` bound moved to the impl; **two overstated
compile-time-enforcement claims corrected** (`From<MateRole>` and `to_chrom_ref_fetch_error` each
check only the *matched* side — the `MateRole` one is now true in both directions via a
`#[cfg(test)]` reverse conversion, the other is a wording fix because an into-mapping owes no
coverage); `MateRole`'s predicates made exhaustive `match`es, since this is the enum ng owns *in
order to extend*; `length()` pinned against production's own method over every op class; the
`RefSeqError::Io` arm and the `narrow` saturation given the tests they lacked; and
`PLACEHOLDER_READ_GROUP` named, because `ReadGroupId(0)` is the run's *first* read group and reads
identically to a defaulted field. Suite 2631 → **2641**.

Three findings are **not** applied and are Checkpoint A questions instead (§8): the per-file copy
banners, the "46 tests" in the plan and the spec, and the arch doc's file inventory.

## 8. Checkpoint A — four questions for the owner

Each reaches past the step, which is why none was decided in code.

1. **The copy banners — RESOLVED: no banner** (owner, 2026-07-29). It would widen the plan.s
   sanctioned edit set and cost `genome_walk.rs` and `errors.rs` their byte-identity, and
   `copy_fidelity.rs` answers the same worry by failing the build with a message that names the
   production original. The reason is recorded in that file, so the decision is met where the
   question arises.
2. **The "46 tests" — RESOLVED: corrected** (owner, 2026-07-29). Plan `:50`, `:75`, `:114` and spec
   `:1050`, `:1054` now read 113 (44 end-to-end + 69 inline), each with a dated note saying what the
   number used to be and when it was counted.
3. **The arch inventory — RESOLVED: fixed** (owner, 2026-07-29). `tests.rs` and `copy_fidelity.rs`
   added, the count corrected to eight copies, with a dated note on why neither is optional. The
   private-`RefSeqFetcher` mismatch was closed in the review pass by narrowing the code to match.
4. **RESOLVED: neither — it is deleted** (owner, 2026-07-29, executed by plan 3's A0 in `da778ab`).
   ng's `open_record.rs` takes a `RefSeq` directly, and the newtype, `to_chrom_ref_fetch_error` and
   both of that translation's lossy spots (a contig *name* rendered as an id, a `u64 → u32`
   narrowing) went with it; no file in `pileup/` imports `MultiChromRefFetcher` or
   `ChromRefFetchError` outside `#[cfg(test)]`. **A rename would have been the wrong fix to the
   right observation:** the adapter existed because the copies' signatures were production's, which
   was a consequence of transcribing verbatim, and the moment the two walkers began to diverge the
   adapter had nothing left to adapt. **With this, all four of Milestone A's questions are closed**
   and plan 2 carries nothing.

   The question as it was asked, kept because the reasoning is what decided it:
   **Should `RefSeqFetcher` be renamed, and should it move to its own file?** The review found the
   name was **deliberately retired** in this codebase for a different concept
   (`fasta/fetcher.rs:20-23`, a 2026-05-23 review), and that the name does not say which way the
   adapter runs. Separately, `pileup/mod.rs` is now the one file that both stands in for
   production's `walker/mod.rs` and holds ng-original code. Arch §1.3 names the type and puts it in
   `mod.rs`, so neither was changed.

## 9. What Checkpoint A hands to Milestone B

- **`walker::tests` really is reachable**, and B1's fixture plan rests on it: ng's copy is
  `pub(crate)` under `#[cfg(test)]`, and `open_record.rs` already imports `MockFasta` across the
  module boundary. Either copy's fixtures can feed the differential; the harness will use one
  `Vec<PreparedRead>` for both walkers, since preparing separately would inject step 2's uppercase
  divergence into a comparison that is about the walk.
- **The two walkers are now genuinely two.** Nothing prevents them drifting except the differential
  B1 builds, and the stage-1 differential **dies in plan 3** by design — B3 is the last moment the
  baseline can be banked.
- **A differential that goes green on the first run is the warning, not the result.** This branch has
  shipped three tests-that-cannot-fail across two milestones. B2 is where that is settled, with its
  own commit.
