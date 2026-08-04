# Code Review: ng_read_filtering_stages_a3 (and Milestone A close-out)
**Date:** 2026-08-03
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A3 — `RegionRecords` → `RegionRawAlignedReads` (file follows the type),
`fill_record` → `fill_raw_read`, `RecordIndex` → `RawReadIndex` — plus an audit of Milestone A
as a whole, this being its checkpoint
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the A3 working-tree diff, exported as
  `tmp/review_2026-08-03_ng-read-filtering-stages-a3/a3.patch` and re-applied by each agent onto
  a detached `db3057a` (A2). One agent additionally audited `8cf6f03..HEAD` + the working tree —
  the whole of Milestone A.
- **In-scope files:** `src/ng/read/input/region_raw_aligned_reads.rs`,
  `src/ng/read/input/aligned_reads_reader/{container,cram,mod}.rs`,
  `src/ng/read/input/{cursor,mod,open_bam}.rs`, `src/ng/read/filtering.rs`, the impl report, and
  for the close-out the plan, spec, arch and `PROJECT_STATUS.md`.
- **Deliberately out of scope:** Milestones B–D; `PspReader::region_records`, a public method of
  the frozen production `.psp` reader and a different concept entirely.
- **Categories dispatched:** `refactor_safety`, `naming`, `extras` (+ close-out).
  `module_structure` was not re-dispatched: A3 moves one file within a directory A2's
  `module_structure` agent had just audited, and its findings there (the `pub(crate) mod`
  question, `container`'s home) are already deferred to this checkpoint.

## 2. Verdict

**Approve-with-changes.** A3 is provably pure. What the review adds is one leftover neither grep
could see, one wrapper the rename made incoherent, three dangling links the `git mv` created,
and — from the close-out — a correction to `PROJECT_STATUS.md`'s own arithmetic.

**Two of its findings are naming objections against names the architecture prescribes.** Those
are not this skill's to overrule; they are raised for the owner in §4.

## 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib` | 0 | 2,839 / 0 failed / 5 ignored — **unchanged from A2** |
| `cargo test --lib ng::` | 0 | 1,540 / 0 failed / 2 ignored — **unchanged from A2** |
| `cargo test --examples` | 0 | 52 passed / 0 failed |

Four dumps byte-identical to the `8cf6f03` baseline by `cmp`; walk probe anchor exact.
Findings labeled "Needs verification": **0**.

**Purity established by inversion, as at A2.** Reverse-substituting all four identifiers across
the eight changed files and diffing against `db3057a`: **four files come back byte-identical**;
`filtering.rs` and `open_bam.rs` differ only by comment re-wrap; `cursor.rs` by re-wrap plus one
rustfmt-split `use` that gained a trailing comma; the moved file by re-wrap plus its title and
diagram. No executable statement, signature, visibility, struct-literal field or string literal
changed anywhere.

**Production is untouched** — `git diff db3057a -- src/psp src/pop_var_caller src/var_calling
src/regions.rs benches examples tests` is empty, and `PspReader::region_records` still stands at
14 occurrences in `src/psp/reader.rs`, 19 across five files. The scoping held.

**Three mutations, each `grep`-asserted present before its run, all killed:**
`fill_raw_read(i, …)` → `(0, …)` (4 failures); `continue_into` also repositioning the reader (25
failures, including `continuing_into_a_region_does_not_move_the_reader`);
`RawReadIndex.mapping_quality` ← `None` (4 failures). The first was re-run by the orchestrator
**after** the wrapper collapse below and still kills.

## 4. Open questions and assumptions

**Both are for the owner at Checkpoint A. Both are objections to names the architecture
specifies, so neither was acted on.**

1. **`RawReadIndex` — arch §2 prescribes it; the `naming` agent says it fails twice.** "Index"
   is not what the value *is* — it is one per-record **entry**, which the surrounding code
   already calls `entry` (field doc, and the binding at `container.rs:210`), living inside a
   field that is itself named `index`. And `RawRead` is a shorthand for A1's `RawAlignedRead`:
   85 occurrences of the full name in `src/ng`, against `RawRead` appearing at exactly the six
   sites A3 created. Proposed: `RawAlignedReadEntry` — which would need arch §2's table amended
   too. Landed as specified; raised here.
2. **`fill_raw_read` takes `&mut RecordBuf`, not a raw aligned read.** Its one caller must fill
   the other half on the very next line (`buf.read_group = …`). `fill_record` was true of that
   signature; `fill_raw_read` is not, quite. The `naming` agent's fix is to change the signature
   to `&mut NoodlesRawAlignedRead` and set both fields — which is beyond a rename, so it was not
   done. The wrapper *collapse* was (see Mi2). Arch §2 prescribes the name.

## 5. Top 3 priorities

1. **Mi1** — the renamed type's own summary sentence still spelled the old name out in a
   different word order, where no grep could see it.
2. **Mi2** — `fill_raw_read` was a one-line wrapper over a private `fill`, and the rename made
   the two ends contradict each other.
3. **Mi3** — the `git mv` left three dangling links, in this step's own governing documents.

## 6. Findings

### Minor

**Mi1: src/ng/read/input/region_raw_aligned_reads.rs:52 — the renamed type's summary still says the old name, in a word order no grep matches**
**Categories:** naming · **Confidence:** High
*"The records of one region of one file"* — `RegionRecords` spelled out and reordered, so
neither `region[ _-]?records?` nor the exact-identifier grep can see it. It is the **first
sentence on the type the step is named after**. Fixed to "The raw aligned reads of one region of
one file". This is the second time in this milestone that a grep-shaped check missed prose:
A2's review found ten "record reader" sites the same way.

**Mi2: src/ng/read/input/aligned_reads_reader/container.rs:151-154 — `fill_raw_read` was a 1:1 wrapper whose two ends the rename made disagree**
**Categories:** refactor_safety · **Confidence:** High
`fill` had exactly one caller — the wrapper — and `fill_raw_read` exactly one — `cram.rs:188`.
A reader following `fill_raw_read` → `fill` was told the operation is about a raw read and then
that it is about nothing in particular. The impl report had deferred this; the agent's argument
for doing it here is right and was accepted: **collapsing a 1:1 name-only indirection is a
rename plus a four-line deletion**, which is precisely what this milestone is for — deferring it
would put a naming edit into a later behaviour diff.

**Mi3: three dangling markdown links, created by this step, in this step's own spec and arch**
**Categories:** refactor_safety, extras (convergent) · **Confidence:** High
`spec/read_filtering_stages.md:154` and `:411`, `arch/read_filtering_stages.md:230` all target
`src/ng/read/input/region_records.rs`, which no longer exists. Distinct from the deferred
design-doc sweep, and the agent was right that `PROJECT_STATUS.md`'s deferral note **predated
A3** and named only A1's and A2's identifiers — so these would have fallen between the two.
A dangling link breaks on click rather than reading oddly.

**Mi4: PROJECT_STATUS.md — the block contradicted itself on the suite count**
**Categories:** extras (close-out) · **Confidence:** High
One line said "suite **2,837** / `ng::` **1,538**, both unchanged" while another in the same
block correctly recorded the +2 from A1's review. My own error, introduced when the block was
written at A1 and not updated when the fixes landed.

**Mi5: the impl report's counts, three of them**
**Categories:** refactor_safety, extras, naming (three-way convergent) · **Confidence:** High
"nine ng files" is **eight** (the moved file counted under both names); "six comment lines
re-wrapped" is **seven blocks, sixteen lines**; and "the module path (4 sites in ng — the other
19 in the tree are production's)" conflates two quantities — the ng four are module paths, the
production nineteen are occurrences of a *method name* and include no module path. Two further
edits went undisclosed: rustfmt reflowed a `use` in `cursor.rs` and added a trailing comma, and
the `held` field's doc was re-wrapped although it contains no renamed identifier. Everything
else in the report reproduced exactly — including, as the close-out agent noted, an
*understatement*: `region_records` appears twice more in `benches/psp_reader_perf.rs` than the
report claimed.

### Nits

`open_bam.rs:1440` was still 103 **characters** after the re-wrap that existed to get it under
100 — a character-accurate audit matters here because these lines carry multibyte `—`, `→` and
`§`, so byte counts mislead. And the layer diagram's gloss had been shortened to "this region's
only", a dangling possessive borrowing its noun from the row above, where the spec's own diagram
spells it in full — fixed by widening the column instead of cutting the words.

## 7. Milestone A close-out

**All eight of arch §2's rename rows landed, under the names the design specified.**
`grep -rnw` for the eight old names across `src/ng` is empty. Spec §6's `RecordSource` →
**deleted** row and arch §1's trait-impl-to-inherent are C3's and correctly untouched.

**The milestone's claim — nothing behaves differently — holds.** The close-out agent read the
full 1,742-line A1+A2 source diff: nothing changes production behaviour. **But the two exceptions
named in its brief were not the only non-rename changes** — there are five more, all test-only
and all documented in `fixes_applied_2026-08-03.md`: a repurposed test, two added assertions,
four function-local flag constants swapped for production's canonical ones (verified
value-identical), and two `..Default::default()` spreads spelled out. Worth stating plainly at
the checkpoint rather than leaving "two exceptions" on the record.

**The evidence, across all three steps:** the four dumps byte-identical by `cmp` at every step —
251,792 / 4,406 / 1,718,914 / 11,945 — and the walk probe's anchor exact every time
(`seconds`: baseline 1.846, A1 1.880, A2 1.871, A3 1.887, same machine and session — noise).
Suite 2,839 / `ng::` 1,540, moved only by A1's two review-added tests.

**Deferred to the checkpoint, deliberately batched:** the two visibility questions
(`NoodlesRawAlignedRead`'s `pub`; `aligned_reads_reader/`'s four `pub(crate) mod`); one
design-doc sweep covering all three renames; and the two arch-prescribed names §4 disputes.
Also noted: `arch/read_filtering_stages.md` §7's open item *"whether `region_records.rs` is
renamed on disk or the type simply moves"* is the question A3 has now answered — but §7 is a
design document, so it was left for the owner rather than edited.

## 8. Missing tests to add now

**None.** A3 adds no code path and changes no invariant; all three renamed items were shown
pinned by live tests, by mutation. The wrapper collapse was re-verified by the orchestrator
after the fact with the same mutation.

## 9. What's good

- **Two reviewers independently converged on inversion as the proof technique** for a rename
  step, and it is strictly better than re-running the dumps: it shows *which* bytes changed
  rather than that the output did not.
- **The scoping was the right call and was reasoned about in advance.** `region_records` really
  is a live production API; a repo-wide substitution inside a renames-only milestone would have
  renamed it, and the impl report named the hazard before the reviewers found it.
- **A character-accurate line-length audit**, rather than a byte-accurate one, on comments full
  of em-dashes and section signs.
- **Every mutation was `grep`-asserted present before its run** — the `cargo fmt` hazard this
  branch has been bitten by, handled without being reminded.

## 10. Commands to re-verify

```
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --lib
cargo test --lib ng::
grep -rn   "RegionRecords\|region_records\|fill_record\|RecordIndex" src/ng     # exact
grep -rniE "region[ -]records|record[ -]index|fill[ -]record|record[ -]readers?" src/ng
grep -rc   "region_records" src/psp/reader.rs                                  # must stay 14
```
Plus the four acceptance dumps and the walk probe against the `8cf6f03` baseline.
