# Code Review: ng_read_filtering_stages_b1
**Date:** 2026-08-03
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step B1 — the per-contig fetch loop becomes a contig-table comparison in
`AlignmentFile::cursor`
**Status:** Request-changes (all applied; see the fix report)

---

## 1. Scope

- **What was reviewed:** the B1 working-tree diff, exported as
  `tmp/review_2026-08-03_ng-read-filtering-stages-b1/b1.patch` and re-applied by each agent onto
  a detached `bfb54dd`.
- **In-scope files:** `src/ng/read/input/{open_bam,cursor,mod,sample_cursor,test_fixtures,
  region_raw_aligned_reads}.rs`, `src/ng/read/filtering.rs`,
  `src/ng/locus_generation/pileup/generator.rs`, `src/fasta/mod.rs`, and the impl report.
- **Categories dispatched:** `reliability` (**the first behaviour-changing step in this plan** —
  what the check does and does not prove), `errors` + the public-API items of `defaults` (a new
  failure path on a `pub` API), `refactor_safety` + `extras` (scope creep in both directions,
  and the impl report's accuracy). `naming` and `module_structure` were not dispatched: B1 moves
  no module and introduces one identifier.

## 2. Verdict

**Request-changes.** The change is right in direction and its arithmetic claim is right, but as
first written it **removed a guarantee it claimed to strengthen**, and its central property —
order-sensitivity — was pinned by nothing.

All findings have been applied; the shipped code is not what was reviewed. The fix report has
the accounting.

## 3. Execution status

Each agent reproduced the orchestrator's numbers independently: `cargo test --lib` 2,837,
`ng::` 1,538, fmt and clippy clean. Findings labeled "Needs verification": **0** — every finding
below was produced by running a mutation or measuring a behaviour, not by reading.

## 4. Open questions

1. **`AlignmentFileError::CursorAccessorContigTable` is a new variant on a `pub` enum, and spec
   §1 says this change adds no new error.** Applied on the `errors` reviewer's argument (§6, M2)
   and recorded as a deviation. **Owner's call at Checkpoint B.**
2. **Spec §9 Q2's cost arithmetic is wrong** in a way two agents confirmed independently. Left
   for the owner, as a design document.

## 5. Top 3

1. **B1** — the order-sensitivity of the comparison, which is the entire reason it is an
   equality test, had no test; an order-insensitive rewrite passed all 1,538.
2. **M1** — a matching table does not mean the accessor can serve the bases, so the fail-fast
   the loop provided was *gone*, not moved. Measured.
3. **M2** — two distinct checks shared one error variant, and the composed message was false.

## 6. Findings

### Blocker

**B1: src/ng/read/input/open_bam.rs — the check's order-sensitivity has no test at `cursor`**
**Categories:** reliability · **Confidence:** High (mutation-verified)

`over_records`' own doc gives the reason the check is an equality and not a resolvability test:
*"A permuted list would resolve on every fetch and make filter #8 compare each read against the
wrong contig's bases, silently."* `with_validated_contigs`' doc repeats it. The new test
asserted only wrong-names and right-names-wrong-lengths — **both order-blind**.

The agent replaced `first_disagreement` with an order-insensitive `HashMap<&str, u64>` lookup
emitting the same message strings. It **survived all 1,538 `ng::` tests**, the new test
included.

**Failure scenario:** this is a performance branch, and an ordered walk of 2,580 entries per
cursor is exactly the shape someone optimises into a hash lookup. That refactor produces wrong
variant calls — filter #8 scoring every read against another chromosome's bases — with no error,
no panic, and a green suite.

**Fix:** the permutation case. Verified to fail under the mutation and pass on the shipped code.

### Major

**M1: src/ng/read/input/open_bam.rs — table equality does not imply the accessor can serve the bases; the loop's fail-fast was removed, not relocated**
**Categories:** reliability · **Confidence:** High (measured)

`InMemoryRefSeq::from_named_contigs` derives its table from the bytes it stores, so for it
equality implies resolvability. `ResidentRefSeq::new` and `WindowedRefSeq::new` /
`with_shared_index` take the `ContigList` as an **independent constructor argument** — nothing
ties it to the bytes on disk. The agent built a `WindowedRefSeq` over a FASTA holding `chr1`
only, behind a table naming `chr1`/`chr2` exactly as the fixture BAM does, and measured:

```
PROBE: old per-contig probe on chr2 would have succeeded? false
PROBE: cursor accepted the accessor? true
PROBE: reads=0 err=Some("… reference read failed on ContigId(1) … contig chr2 not in FASTA index")
```

**Failure scenario:** a stale `.fai`, or a `.fai` whose FASTA was replaced. The run now fails
after arbitrary work, per chromosome, under a top-level message naming the **BAM** — and on a
`--regions` run whose reads are all dropped before filter #8, never at all.

This also makes the impl report's "proves strictly more" false, and leaves
`AlignmentFileError::Reference` with no producer while its doc still promises the fail-fast.

**Fix (applied):** keep the comparison, and add **one** zero-length fetch on this cursor's own
contig — one `open(2)` instead of 2,580. Restores the guarantee, revives the variant, makes the
doc true again.

**M2: src/ng/read/input/open_bam.rs — two distinct checks shared one error variant, and the composed message was false**
**Categories:** errors · **Confidence:** High (triggered, message captured verbatim)

The first implementation raised the open gate's `ContigReconcile` and encoded which check fired
by prefixing the free-form `detail`. The message an operator reads:

> alignment file '…/sample.bam' does not match the reference contig table: the accessor passed
> to cursor() is over a different table: name disagreement at index 0 ('chr1' vs 'not_chr1')

**The headline is false at that point in the run** — `open` has already proved the file *does*
match. Two clauses that both say "table" contradict each other before the true one arrives, and
the two failures are discriminable only by substring, which the new test then did.

`errors.md`'s split test is explicit: *"if two operations have different failure modes **or**
should display different messages, split them"* — and these already displayed different
messages, through a string prefix rather than the type.

**Fix (applied):** `AlignmentFileError::CursorAccessorContigTable`. See open question 1.

**M3: src/ng/read/input/open_bam.rs — neither the operand order nor the disambiguating prefix was asserted**
**Categories:** reliability · **Confidence:** High (mutation-verified)

Two recorded design decisions with zero regression protection: swapping the operands *and*
deleting the prefix survived all 1,538 tests, because the assertions matched only
`contains("name disagreement")`. **Fix (applied):** the assertions now pin the rendered shape.

### Minor

**Mi1** — `AlignmentFileError::Reference` had no producer, and its doc asserted a fail-fast the
code no longer performed ("surfaced here rather than at the first read"). Resolved by M1's fix,
which makes it true again; doc updated to scope it to one contig.
**Mi2** — three test call sites construct cursors via `over_records`, bypassing the check, with
nothing asserting their fixtures satisfy its precondition — the exact state the whole tree was
in before B1. **Fix (applied):** a standing guard test.
**Mi3** — the `@SQ`-count mismatch, `first_disagreement`'s third branch, was unreached at
`cursor` by any test, and is the input class the *deleted* test actually drove. **Applied.**
**Mi4** — `ReadFilterError::Reference`'s doc explained itself through the deleted
`ReadFilter::new` and read backwards. **Applied.**
**Mi5** — two broken intra-doc links to the deleted `ReadFilter::new`, neither among the 12
pre-existing. **Applied.**
**Mi6** — the new module comment in `filtering.rs` ended "the alias lives in the test module"
while the same diff removed that alias, and a comment 1,500 lines later said no test needs a
header. **Applied.**
**Mi7** — `SampleReads::cursor` gained a failure mode its doc did not mention. **Applied.**
**Mi8** — the impl report's own defects: "four imports" (three); the citation `:652` (the
comment spans `:622-655`); the "~18 µs → 46 ms" replacement estimate inherits the source
comment's `open(2)` attribution, which cannot hold for `with_shared_index`; and the deleted-test
accounting was half wrong. **All applied.**

### What the reviewers verified with no finding

`over_records`' infallibility is genuine — four call sites audited, no information swallowed.
The `+ ContigTable` bound reached everything it should, with no stub impls (both spies
delegate). The fixture change is a **fix, not a weakening**, verified three ways. The two
constructor-comparison tests genuinely have no possible successor, and the input class they
carried (a MAPQ-0 record through the full iterator) is still covered — mutating `record_drop`'s
`LowMapq` arm is caught by `a_walk_charges_every_drop_reason_by_hand_count`. The deletions of
`ReadFilter::new` and `RecordSource::header` are **compiler-forced**, not scope creep: both
produce `dead_code` under `-D warnings` once the loop goes.

## 7. Out of scope observations

- `ResidentRefSeq::new` and `WindowedRefSeq::new` accept a `ContigList` unrelated to the bytes
  they will serve — a constructor that can build a lying accessor. M1's check contains the
  damage at the cursor; the constructors are an API-design question of their own.
- 21 pre-existing `..Default::default()` test literals in the in-scope files, none introduced
  here.

## 8. Missing tests

All five were written, run, and mutation-verified by the reviewers, and all five are applied.
See §5 of the impl report for the table and the four mutations.

## 9. What's good

- **Three reviewers, three different instruments**: one measured a real accessor failing, one
  captured verbatim error messages by triggering both paths, one inverted the diff. None
  reviewed by reading.
- **The estimate correction was independently confirmed** by two agents reading the cited
  comment in place.
- **The fixture defect was found by the change, not by the tests** — and the reviewer proved the
  repair was a repair by reverting it and counting the failures back to exactly 23.

## 10. Commands to re-verify

```
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --lib
cargo test --lib ng::read::input::open_bam
```
Plus the four acceptance dumps and the walk probe against the `8cf6f03` baseline, and the
timing A/B against `bfb54dd`.
