# Code Review: ng_locus_generation_pileup_port_a
**Date:** 2026-07-28
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Milestone A of the ng generic locus generator port — production's pileup walker copied verbatim into `src/ng/locus_generation/pileup/`, plus ng's own `PreparedRead` and the `RefSeq` → `MultiChromRefFetcher` shim
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the diff `6d3fe18..HEAD` on branch `ng-pileup-generator` — four commits,
  `9bfd483` (A1), `8b6307d` (A2), `6e44051` (A3), `8fdba95` (A4). 8,124 insertions over 16 files.
- **Reviewed against:** `8fdba95`, branch `ng-pileup-generator`.
- **In-scope files** — ng's own code and the seams, the only places a finding can be actionable:
  - [prepared_read.rs](../../../../src/ng/read/prepared_read.rs) (new)
  - [read/mod.rs](../../../../src/ng/read/mod.rs), [left_align.rs](../../../../src/ng/read/left_align.rs), [left_align_parity.rs](../../../../src/ng/read/left_align_parity.rs) (changed)
  - [pileup/mod.rs](../../../../src/ng/locus_generation/pileup/mod.rs) (new)
  - [locus_generation/mod.rs](../../../../src/ng/locus_generation/mod.rs) (one added line)
- **Deliberately out of scope:** the *contents* of the seven copied walker files and `tests.rs` under
  `src/ng/locus_generation/pileup/`. They are a deliberate verbatim copy of frozen production code;
  a finding that proposed improving them would destroy the property the milestone exists to
  establish. Their **fidelity** was checked instead, mechanically, by every category. Also out of
  scope: `src/pileup/`, `src/psp/`, `src/var_calling/`, `src/vcf/` (frozen, so findings there are not
  actionable), and `benches/psp_writer_perf.rs` (pre-existing failure).
- **Categories dispatched (9):** `reliability` (always), `errors` (always — the shim's error
  translation is the densest new surface), `naming` (always), `defaults` (a re-exported constant
  family and a placeholder id), `idiomatic` (always), `refactor_safety` (always — and the highest
  value here: a 5,495-line duplicate plus four cross-type conversions), `module_structure` (the diff
  spans a new folder and changes ng's public surface), `smells` (always), `extras` (a PR-shaped diff,
  so "diff matches stated intent" is live and is the milestone's whole claim).
  `unsafe_concurrency` was **not** dispatched: no `unsafe`, no threading, no channels, and the only
  `Arc` is `PreparedRead::qname`, transcribed unchanged and never shared across threads by this diff.

## 2. Verdict

**Approve-with-changes.** The milestone did what it claimed. The freeze holds absolutely, the copy is
verbatim, and the inherited suite is green against it. The three Majors are all about the *seams*
being weaker than their own documentation says — which is this branch's recurring failure mode, and
is why they are Majors rather than Minors.

## 3. Execution status

Run in the container per `CLAUDE.md`; verbatim output in
`tmp/review_2026-07-28_ng-pileup-port-a/verification.txt`.

| command | exit | result |
|---|---|---|
| `cargo fmt --all -- --check` | 0 | no output |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | no diagnostics |
| `cargo test --all-targets --all-features` | 101 | lib suite **2631 passed / 0 failed / 4 ignored**, every integration target green; the sole failure is the **pre-existing** panic at `benches/psp_writer_perf.rs:386` |
| `cargo test --release --lib ng::locus_generation::pileup` | 0 | 118 passed — the release half of the `cfg(debug_assertions)` test pair |
| `cargo doc --no-deps` | — | **not run**: known-red on 11 unresolved intra-doc links |
| `cargo audit` | — | **not run**: not installed in the dev container |

Both skipped commands are tracked under PROJECT_STATUS *Standing project-wide items* and are red
independently of this work.

Findings labelled "Needs verification": **0**. Every finding cites a location a sub-agent read, and
the two language-level claims behind Mi5 were checked with a standalone `rustc` compile rather than
asserted.

## 4. Open questions and assumptions

1. **May the copied files carry a "this is a copy, do not edit" banner?** (affects **M3**.) Adding
   one to each of the seven widens the plan's sanctioned edit set — A2 permits "the module paths, and
   `PreparedRead` resolving to ng's" and nothing else — and it costs the byte-identity of
   `genome_walk.rs` and `errors.rs`, currently the cheapest possible fidelity statement. An owner
   decision, not the implementer's.
2. **The plan and the spec still say "46 tests".** (affects **Mi14**.) The real numbers are 44
   end-to-end + 69 inline; three lines of the plan and two of the spec carry the wrong one. Design
   documents are not this skill's to edit.
3. **The arch doc's *Module home* inventory omits `tests.rs` and specifies a private
   `RefSeqFetcher`.** (affects **Mi15**, **Mi9**.) Both are arch-doc edits.
4. **Should `RefSeqFetcher` be renamed and/or moved to its own file?** (affects **Mi16**, **Mi17**.)
   The arch doc names the type and puts it in `mod.rs`; the review found the name was *deliberately
   retired* in this codebase for a different concept, and that a separate file would make the
   copied/ours boundary a file boundary.

## 5. Top 3 priorities

1. **M1** — `assert_same_prepared_read` documents compile-time enforcement it does not have, and it is
   now the only place a production field can vanish from the parity comparison silently.
2. **M2** — `RefSeqFetcher` is the milestone's only authored *behaviour* and is never driven through
   the walker it exists to feed.
3. **M3** — the copies carry no in-file marker and no mechanical fidelity check; the only evidence
   they are verbatim is a `diff` in a gitignored scratch directory.

## 6. Findings

### Blocker

None.

### Major

**M1: src/ng/read/left_align_parity.rs:190 — `assert_same_prepared_read` claims a field addition forces it to be updated; nothing enforces that**
**Categories:** refactor_safety, reliability (convergent) · **Confidence:** High

The doc says "Listing the fields means a field added to production forces this to be updated rather
than passing vacuously." It does not. The body is twelve independent `assert_eq!` calls; a thirteenth
field on either type leaves this compiling and passing, uncompared. This is the one place the
compiler-enforcement story breaks — `from_production` *does* destructure exhaustively, so a new field
is forced through the conversion into ng's type and then silently escapes the thing the module calls
"the port anchor". **Fix:** exhaustive destructures of both sides, no `..`, `read_group: _` named.

**M2: src/ng/locus_generation/pileup/mod.rs:114 — `RefSeqFetcher` is never driven through the walker it exists to feed**
**Categories:** reliability · **Confidence:** High

`grep` returns only the definition, the impl, and five isolated `shim_tests` construction sites.
Nothing composes it with `run`/`PileupWalker`, so the 113 inherited tests prove the walk only over
`MockFasta`. Untested in particular is [open_record.rs:378](../../../../src/ng/locus_generation/pileup/open_record.rs#L378)'s
`widen`, which passes a record's **exclusive** footprint end as a **1-based** start — a disagreement
there yields wrong REF bases with no error. The reviewer hand-checked that `MockFasta` and
`validate_window` agree on the bounds predicate today, which is why this is Major and not Blocker.
**Fix:** a shim-vs-`MockFasta` differential walk over the same reads, with a deletion so `widen` is
on the path, using `PileupRecord`'s `Debug` as the oracle.

**M3: src/ng/locus_generation/pileup/*.rs — the copies carry no in-file marker, and the fidelity evidence does not survive the branch**
**Categories:** smells, refactor_safety, reliability (convergent) · **Confidence:** High

Six of the seven copies say nothing about being copies; `genome_walk.rs` is the worst case, since its
header is production's verbatim and nothing in it records either that it is a copy *or* that it was
renamed — and the rename is what a reader needs to find the file to diff against. Meanwhile the only
artefact establishing the verbatim property is `tmp/…/copy_fidelity.txt`, and `.gitignore:3` is
`/tmp/`. Nothing in `precommit-check.sh` or CI compares the two trees; an edit on either side
compiles and both suites stay green, because they are independent copies.
**Fix, in two halves:** a committed, mechanical fidelity check (a test or a script) — and, separately,
the per-file banner, which is Open question 1 because it widens the plan's sanctioned edit set.

### Minor

**Mi1: src/ng/locus_generation/pileup/mod.rs:171 — the `RefSeqError::Io` arm has no test.** `shim_tests`
covers `OutOfBounds`, `InvalidStart` and `UnknownContig`; `Io` is unreachable from `InMemoryRefSeq`,
the only reference those tests use. It is the arm that fires on the one failure a user actually hits —
a broken reference file. *Fix:* a `BrokenRefSeq` stub. (reliability)

**Mi2: src/ng/locus_generation/pileup/mod.rs:149 — the `narrow` saturation invariant is documented and unpinned.**
Reachable through the shim's own signature: `fetch(0, u32::MAX, 10)` makes `end` exceed `u32`, and
only `narrow` keeps it from wrapping to a value *below* its own `start`. (reliability)

**Mi3: src/ng/read/prepared_read.rs:211 — nothing checks the transcribed `length()` against the original.**
This is the one file that was hand-transcribed rather than byte-copied, and `length()` is the one
function in it that computes anything. Production's `walker/mod.rs` has no inline `mod tests`, so no
test came across with the type; the five new ones pin ng against ng. *Fix:* a table-driven agreement
test across every op class and both failure modes. (reliability)

**Mi4: src/ng/read/left_align_parity.rs:217 — the `mate_role` assertion routes both sides through the same conversion.**
`ours.mate_role` already *is* `production.mate_role.into()`, so a `From` that collapsed two roles
would make both sides equally wrong. Independently covered by
`every_production_mate_role_maps_to_its_counterpart`, so nothing is unproven — but a parity test that
reads as self-contained and silently depends on a test in another file is the shape that decays.
(reliability, smells — convergent)

**Mi5: src/ng/locus_generation/pileup/mod.rs:66 — the `super::` vocabulary block is `pub`.**
It exists only so `super::Foo` resolves inside private child modules, and a `pub(crate)` (or private)
binding does that — verified with a standalone `rustc` compile. As written it mints public paths for
frozen production items under an ng namespace, plus a *third* public path to ng's own `PreparedRead`.
Nothing outside the module consumes any of them. One wrinkle constrains the fix, also verified:
`DEFAULT_MAX_SNP_COLUMN_DEPTH` and `DEFAULT_MAX_INDEL_COLUMN_DEPTH` are reached by no copy, so
demoting them trips `unused_imports` under `-D warnings` — drop them or keep them `pub` with a
stated reason. (module_structure, idiomatic — convergent)

**Mi6: src/ng/locus_generation/pileup/mod.rs:44-47 — the "one source of truth" claim holds for four of five constants.**
`DEFAULT_MAX_ACTIVE_READS` is re-exported from ng's *copy* of `chain_id_allocator.rs`, so two
definitions of `= 4096` now exist. The duplication is forced by the verbatim rule; the doc describing
it as if it were not is the defect. Visible in one line of the copy,
`chain_id_allocator.rs:132`, whose first operand is ng's and second production's.
*Fix:* correct the doc and add a test pinning the two equal. (defaults, module_structure, naming — convergent)

**Mi7: src/ng/read/prepared_read.rs:96-103 — the `#[non_exhaustive]`/no-`Default` decision is a `//` block, which rustdoc drops.**
The two decisions governing how a caller must construct the type are visible only in source. The style
is production's, but ng already rewrote the text, so promoting it costs nothing. (defaults)

**Mi8: four sites — `ReadGroupId(0)` reads as a real group, not a placeholder.**
`types.rs:171-176` makes 0 the run's *first* `@RG`; the crate already recorded this exact trap at
`filtering.rs:494-501` ("zero is not a sentinel"). Only one of the four in-scope sites says
"placeholder". The failure mode is a vacuously-passing assertion — precisely the check the read group
was added to make possible. *Fix:* a named `PLACEHOLDER_READ_GROUP`. **Not** in the copied fixtures,
where that literal *is* the transcription seam. (defaults)

**Mi9: src/ng/locus_generation/pileup/mod.rs:114 — the shim is `pub` with a `pub` field where the arch doc specifies a private newtype, and the deviation is unrecorded.**
`arch:217` says `struct RefSeqFetcher<R: RefSeq>(R);`. Nothing in this milestone needs either `pub`;
B1's `parity.rs` is a sibling inside `pileup/` and reaches a private item through `super::`. A `pub`
tuple field also lets a caller reach past the shim to the underlying `RefSeq`, which is the opposite
of "moves bytes and decides nothing". (extras)

**Mi10: src/ng/read/prepared_read.rs:72 and pileup/mod.rs:135 — two conversions claim enforcement "on either side"; only the source side is checked.**
A variant added to *ng's* `MateRole`, or to `ChromRefFetchError`, compiles silently. For `MateRole` a
`#[cfg(test)]` reverse conversion makes the claim true; for `to_fetch_error` no construct can — an
into-mapping owes no coverage — so that one is a wording fix. This branch has already had to correct
one overstated invariant (`6d3fe18`); these are two more of the same kind.
(refactor_safety, smells — convergent)

**Mi11: src/ng/read/prepared_read.rs:58 — ng's `MateRole` predicates use `matches!`.**
A fourth variant is silently classified as paired-and-not-first, and those two predicates drive
mate-lookup registration and the equal-BQ tie-break. Verbatim from production — but this is the enum
ng owns *in order to extend it*, and production's is frozen where ng's is not. (refactor_safety)

**Mi12: src/ng/locus_generation/pileup/mod.rs:64-65 — the facade claims to mirror production's `walker/mod.rs` exactly; it omits `indel_norm` with no note.**
The omission is correct — nothing in the seven reaches it — but this is the file whose job is to make
the two module trees comparable side by side. (smells)

**Mi13: src/ng/locus_generation/pileup/mod.rs:116 — `iter_bases` is left on its default, which materialises a whole contig.**
Dead today. But the shim is the designated adapter between ng's references and anything wanting a
`MultiChromRefFetcher`, and the next consumer through a `WindowedRefSeq` would get exactly the
whole-contig residency that reference exists to avoid, with no error and no test failing.
(smells, errors — convergent)

**Mi14: the implementation report's "per commit: `cargo fmt --check` — exit 0" is false at A2.**
A2 landed `cigar_cursor.rs` with its test-module imports out of rustfmt order; A3 fixed it, and A3's
message describes only the shim. The reviewer checked A2's version of that file back into the tree
and ran the command: `FMT_EXIT=1`, the diff being exactly that reorder. In a milestone whose
deliverable is a closed edit list, an unlisted edit — even a cosmetic reorder — is the class of slip
the list exists to catch. Related: **"Edits inside the copies, and this is the whole list" is not the
whole list** — it omits the four `use crate::ng::types::ReadGroupId;` imports, the A3 reorder, the
`open_record.rs` `MockFasta` repoint, and `tests.rs`'s doc block. (extras)

**Mi15: the plan and the spec still say "46 tests" in five places.** The diff's only plan edits are
four checkbox flips, so A4 is ticked ✅ against a criterion this milestone disproved, and the
correction lives only in a dated report. Plan `:50`, `:75`, `:114`; spec `:1050`, `:1054`. See Open
question 2. (extras)

**Mi16: src/ng/locus_generation/pileup/mod.rs:114 — `RefSeqFetcher` does not say which way it adapts, and re-uses a retired name.**
Both traits are fetch-shaped, so the name reads equally as "fetches `RefSeq`s". The crate's three
other fetcher types keep `RefFetcher`/`ChromRefFetcher` as the head noun. Worse,
[fasta/fetcher.rs:20-23](../../../../src/fasta/fetcher.rs#L20) records that a `RefSeqFetcher` trait
was **deliberately retired** in a 2026-05-23 review; a grep now returns both the retirement note and
a live type with an unrelated meaning. The arch doc names the type, so see Open question 4.
(naming, smells — convergent)

**Mi17: src/ng/locus_generation/pileup/mod.rs:82-268 — ng's own code shares the file whose job is to mirror production's `walker/mod.rs`.**
268 lines, ~180 of them ng-original. The one file a reader would diff against `walker/mod.rs` to
check the mirroring is the one file where the mirroring deliberately does not hold. A separate
`ref_seq_fetcher.rs` makes the copied/ours boundary a file boundary — the cheapest kind to respect.
See Open question 4. (smells)

**Mi18: src/ng/mod.rs — ng's own statement of the freeze is stale, and the largest dependency is invisible from the entry point.**
It states the policy as "does not edit `src/ssr/` or `src/regions.rs`", naming neither `src/pileup/`
nor the other three modules this milestone is frozen against, and its "Landed so far" list omits
`locus_generation` entirely — so a reader does not learn that a 5,495-line duplication exists at all.
(module_structure)

**Mi19: src/ng/read/prepared_read.rs — `ReadLengthError` implements neither `Display` nor `std::error::Error` and is not `#[non_exhaustive]`,** unlike every other ng error type. Faithful to
production, so filed as additive-and-deferrable rather than as something to fix now: filing it against
the transcription would attack the property the milestone exists to establish. (errors)

### Nits

Collected, not enumerated. `iter_bases`' inherited default deserves one doc line
(reliability, smells). `the_shim_canonicalises_what_it_serves` should say what distinguishes it from
`ref_seq.rs`'s own canonicalisation test, or a later reader deletes it as a duplicate.
`the_read_group_rides_through_both_paths` uses `ReadGroupId(7)` twice, where two different ids would
additionally rule out a preparer that latched the first group it saw. A `matches!` pattern uses a bare
`..` where `chrom_name: _` would flag a rename. `to_fetch_error` names half its target type where the
sibling `to_aligned_read` names its in full. `passthrough` is a noun-named function wrapping a
verb-named one. `shim_tests` names the pattern rather than the subject; the crate's precedent
(`mod baq_tests`) is subject-named. The `Copy` derive on `RefSeqFetcher` is inert — no `RefSeq` impl
in the tree is `Copy` — while advertising bitwise-copy semantics for a reference accessor of unbounded
size; and its `R: RefSeq` bound belongs on the impl, not the struct definition. The import block is
split by the re-export block. `crate::pileup::walker::{MateRole, PreparedRead}` is spelled
fully-qualified at a dozen sites where two aliased imports would make "ours vs production's" read at a
glance. Two sibling files import `PreparedRead` through the inner module rather than the curated
`crate::ng::read::` surface. The new `passthrough` helper hardcodes `chrom_id` `0`, ignoring
`read.ref_id`, where the shipping path derives it — the literal predates the diff, but hoisting it
into a helper is what hid it. Three transcribed doc references point at `ia/reviews/`, a directory
that no longer exists (the files are under `doc/devel/reports/reviews/`). `#[cfg(test)] pub(crate) mod
tests;` is justified as "reachable from the parity harness", but the harness is a *descendant* and
would reach a private `mod tests` — "mirrors production's declaration" is the true reason.
`copy_fidelity.txt` names no generating command, so it cannot be re-derived from the artefact alone.
Five small prose slips in the impl report: the strip-and-diff claim leaves two hunks not one; "the
five `length()`/`MateRole` tests" is six; `cargo doc` is quoted as an observed failure where the
record says not-run; and "5,495 lines" is the *source* count, 5,503 landed.

## 7. Out of scope observations

- **[locus_generation/mod.rs:681](../../../../src/ng/locus_generation/mod.rs#L681)** — the standing
  `FIXME(pileup-generator)` on `SampleLocusObservationsIterator`'s field-drop order names *this
  plan's* generator as the trigger that turns a latent silent under-report of step-1 drop tallies into
  a live one. Untouched by this diff — no generator holds a read stream yet — but it becomes
  actionable the moment plan 3 wires the walk to `SampleReads`. Raised by three categories
  independently. **Follow-up:** whoever sequences plan 3.
- **[left_align.rs](../../../../src/ng/read/left_align.rs)** — a pre-existing
  `expect("ref_id fits u32")` outside the changed lines; noted, not filed.

## 8. Missing tests to add now

Grouped by function under test. Full code for each is in the per-category files under
`tmp/review_2026-07-28_ng-pileup-port-a/`.

**`RefSeqFetcher::fetch` / the shim as a whole**
- `a_walk_over_the_shim_matches_the_same_walk_over_the_mock_reference` — same reads, two fetchers,
  `PileupRecord`'s `Debug` as the oracle, a deletion so `widen` is on the path. Catches any
  coordinate-convention disagreement between ng's reference and the walker (M2).
- `a_broken_reference_read_is_reported_as_io_with_its_source_intact` — the `RefSeqError::Io` arm, via
  a stub `RefSeq`. Catches a dropped `source` or a mis-routed variant (Mi1).
- `an_out_of_bounds_window_wider_than_u32_saturates_rather_than_wrapping` — `fetch(0, u32::MAX, 10)`.
  Catches a truncating cast reporting an `end` below its own `start` (Mi2).
- `a_zero_length_window_is_served_rather_than_refused` — at, and one past, the last base. Catches a
  divergence from `MockFasta` on the empty window the walker's exclusive ends sit at.

**`PreparedRead::length`**
- `length_agrees_with_productions_on_every_op_mix` — table-driven against production's own method,
  every op class and both failure modes. The only check on the one hand-transcribed file (Mi3).
- `a_read_with_no_bases_reports_a_length_of_zero` — the degenerate read; `0`/empty are the boundary
  classes `length()` currently misses.

**`LeftAlignPreparer::prepare_read`**
- `two_reads_from_different_groups_keep_their_own` — two distinct non-zero groups through one preparer
  and one scratch. Catches a preparer that latches the first group it sees, which is exactly the bug
  `LeftAlignScratch` invites and which the current single-id test cannot see.

**The copy itself**
- a fidelity test over all eight pairs (M3), whose failure message names the production original and
  says the file must not be edited.

## 9. What's good

- **`PreparedRead::from_production` destructures rather than field-copies**
  ([prepared_read.rs:258](../../../../src/ng/read/prepared_read.rs#L258)) — a field added to
  production's type is a compile error, not a silent drop. It is the pattern M1 asks the parity
  comparison to adopt.
- **`MateRole` is textbook co-dependent-bool elimination** — production's paired/first-of-pair bits
  collapsed into three variants so "solo but first-of-pair" is unrepresentable.
- **The fixtures give same-typed fields distinct values on purpose**
  (`the_conversion_from_productions_read_moves_every_field` uses `2`/`101`/`140` for three `u32`s;
  `only_the_read_consuming_ops_count_toward_the_length` gives every CIGAR op a distinct count) — the
  reliability pass put all 14 new tests through the single-character-mutation question and found none
  unfalsifiable.
- **The shim's semantic-emptiness claim is checked, not asserted** — a contig written `acgtRYacgt`
  must come back `ACGTNNACGT`, which fails against the one mutation that matters (repointing at
  `RawRefSeq`).
- **The `DEFAULT_*` constants are reached by name, not by literal**, so production and ng cannot drift
  silently — for four of the five, which is what makes Mi6 worth stating rather than a quibble.

## 10. Commands to re-verify

Reviewer ran (container, `./scripts/dev.sh`): `cargo fmt --all -- --check`,
`cargo clippy --all-targets --all-features -- -D warnings`,
`cargo test --all-targets --all-features`, `cargo test --release --lib ng::locus_generation::pileup`,
`cargo test --lib -- --list`.

Introduced by this review: the eight-pair fidelity diff (M3), and the tests in §8.

### Author response convention

Address each finding by identifier (`M1`, `Mi7`, …) with `fixed in <commit>` / `disputed because …` /
`deferred to <issue>` / `won't fix because …`. Answer the open questions in §4 first.
