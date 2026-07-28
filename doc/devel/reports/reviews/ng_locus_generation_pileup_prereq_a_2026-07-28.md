# Code Review: ng generic locus generator — prerequisites, Milestone A

**Date:** 2026-07-28
**Reviewer:** rust-code-review skill (orchestrator), 6 category sub-agents
**Scope:** the Milestone A diff — `src/ng/read/input/` returns an owned region stream
**Status:** Approve-with-changes — **all findings applied or recorded; see §11**

---

### 1. Scope

- **What was reviewed:** a diff — `e1718d3..270b9d6`, three commits (`d1db8dd` A1, `60acd15` A2,
  `270b9d6` A3).
- **Reviewed against:** branch `ng-pileup-generator` at `270b9d6`.
- **In-scope files:** [open_bam.rs](../../../../src/ng/read/input/open_bam.rs),
  [region_query.rs](../../../../src/ng/read/input/region_query.rs),
  [merge.rs](../../../../src/ng/read/input/merge.rs),
  [mod.rs](../../../../src/ng/read/input/mod.rs).
- **Deliberately out of scope:** `src/pileup/`, `src/psp/`, `src/var_calling/`, `src/vcf/`
  (production, frozen — this work edits none of it); `src/ng/locus_generation/` (plans 2 and 3);
  the doc and `PROJECT_STATUS` edits inside the range.
- **Categories dispatched:** `unsafe_concurrency` (the change introduces three `Arc`s into a
  `Mutex`-pooled, `Sync` type), `reliability` (the milestone's whole claim is "nothing moved"),
  `refactor_safety` (a representation change to a shipped module), `idiomatic` (`&Arc<Self>` is a
  receiver choice worth challenging), `naming` + `module_structure` (many doc comments rewritten),
  `smells` + `extras` (hot path; and "diff matches stated intent" against the plan, arch and impl
  report). Not dispatched: `defaults` (no default-acting value changed), `errors` (no error path
  touched), `tooling` (no `Cargo.toml` change).

**One process deviation, recorded not hidden.** The plan-driven skill prescribes a review per
step; this is one review over the milestone's three commits. A1 and A2 are one change taken in two
commits, and the reviewable question — "is the ownership change correct?" — is only answerable over
both. Three fan-outs would have produced three copies of the same findings on a ~200-line
representational diff.

### 2. Verdict

**Approve-with-changes.** No Blockers. One Major, which was real: the milestone's one
newly-reachable *run-time* behaviour was asserted nowhere, and the test that claimed to cover it
could not. It is now covered by a test that was **mutation-verified** to fail (§3).

The change itself is sound. Three categories independently confirmed the core properties: no `Arc`
cycle, lock discipline unchanged, every production `Arc::new` once per open (never per query), no
`Arc::clone` added inside any read loop, `Drop` ordering unchanged, and the `Sync` requirement not
merely preserved but incidentally strengthened — `RegionReads: Send` now proves
`AlignmentFile: Send + Sync` through the `Arc` where it previously proved `Sync` alone.

### 3. Execution status

Run in the container (`./scripts/dev.sh`), quoted verbatim:

- `cargo fmt --all --check` — exit 0, no output.
- `cargo clippy --all-targets --all-features -- -D warnings` — no diagnostics.
- `cargo test --all-features` — `test result: ok. 2489 passed; 0 failed; 4 ignored; 0 measured;
  0 filtered out; finished in 42.29s` at review time; **2493 after the fixes**.
- Also verified natively on the host: `cargo test --lib` — 2489 passed at review time.

**Not run, with reason:** `cargo test --all-targets --all-features` (pre-existing panic in
`benches/psp_writer_perf.rs`) and `cargo doc --no-deps` (11 pre-existing unresolved intra-doc
links) — both tracked under PROJECT_STATUS *Standing project-wide items*. `cargo audit` not run;
no dependency changed.

**Mutation verification of the fix.** `RegionReads::drop`'s body was replaced with a discard and
the suite re-run. `a_stream_outliving_every_other_handle_still_banks_its_reader_and_tally`
**FAILED**, along with 11 pre-existing tests — while
`a_region_stream_outlives_the_sample_reads_it_was_made_from` and
`a_merged_region_stream_outlives_…` both stayed **green**, confirming the review's claim precisely:
the detached tests are compile-time anchors and cannot see the drop path. `test result: FAILED. 118
passed; 12 failed`. The probe was reverted and the suite is green again (130 passed in
`ng::read::input`).

Findings labeled "Needs verification": 0.

### 4. Open questions and assumptions

1. **Should `AlignmentFile::open` return `Arc<Self>`?** (affects M4) Ten call sites now chain
   `.map(Arc::new)`, because the `Arc`-ness is an invariant of *using* the type but is enforced by
   convention at each site. **Not applied** — it is a change to a `pub` constructor's signature,
   beyond the plan's ask, which makes it a checkpoint decision rather than an implementer's.
2. **Should the share move inside `ReadGroupResolution`?** (affects M5) `PerRecord(Arc<[…]>)`
   instead of `Arc<ReadGroupResolution>` would drop the wrapper from three fields and both source
   constructors, since `Sole(ReadGroupId)` is `Copy`-sized. **Not applied** — it changes a shipped
   type's shape in `read_groups.rs`, outside A2's representational remit.
3. **Who folds in the two specs the change invalidated?** (affects M3) The plan-driven skill
   forbids this run from editing design docs. Recorded as a `SPEC-FOLLOWUP` for the owner, in the
   project's established idiom.

### 5. Top 3 priorities

1. **M1 (Major)** — the orphaned tally: the one new run-time behaviour was untested and
   mis-documented. *Applied: a test that observes it through a retained handle, mutation-verified.*
2. **M2** — three doc comments the change made false, surfaced by five categories independently.
   *Applied.*
3. **M6/M7** — the `Merged` arm (the k-file production shape) and the interleaved-queries property
   were both uncovered. *Applied as two tests.*

### 6. Findings

#### Major

**M1: [src/ng/read/input/mod.rs:1226](../../../../src/ng/read/input/mod.rs#L1226) — a stream
outliving its `SampleReads` silently discards that query's step-1 tally, and the new test
institutionalises the ordering without asserting anything about it**
**Categories:** reliability, unsafe_concurrency, refactor_safety, extras, smells (convergent, 5)
**Confidence:** High.

After `drop(sample)` the `Arc<AlignmentFile>` is reachable only through the stream —
`SampleReads::files` is private and `counts()` needs `&self` — so the tally `RegionReads::drop`
folds in lands in an object nobody can read, and is freed with it. The test's own doc-comment
claimed it reached that path; it did not, and could not. **Proven by mutation** (§3).

*Applied.* Three parts: `a_stream_outliving_every_other_handle_still_banks_its_reader_and_tally` in
`open_bam`'s tests, which keeps a second `Arc` on purpose and asserts `Arc::strong_count`, the
returned reader and the folded tally; the overclaiming comment on the `mod.rs` test replaced with
what it actually proves (a compile-time bit) plus a pointer to the test that proves the rest; and
the consequence documented on both `AlignmentFile::counts` and `SampleReads::reads_in_region`.

**Verification pointer worth carrying forward, out of scope here:**
`SampleLocusObservationsIterator` ([locus_generation/mod.rs:485](../../../../src/ng/locus_generation/mod.rs#L485))
declares `reads: SampleReads` **before** `generators`, and Rust drops fields in declaration order —
so once a generator holds a stream, the sample dies first and every still-held stream's tally
becomes unobservable. That turns this from latent to live. Recorded in PROJECT_STATUS for plan 3.

#### Minor

**M2: [open_bam.rs:72](../../../../src/ng/read/input/open_bam.rs#L72),
[region_query.rs:417](../../../../src/ng/read/input/region_query.rs#L417),
[mod.rs:327](../../../../src/ng/read/input/mod.rs#L327) — three doc statements the change
invalidated**
**Categories:** naming, smells, refactor_safety, idiomatic, reliability (convergent, 5)
**Confidence:** High. *Applied.*

`AlignmentFile`'s type doc said *"sharing is by reference"* — contradicted 370 lines later by its
own method doc; `CramRegionSource.resolution` still said *"Borrowed from the `AlignmentFile`"* while
its type is now an `Arc` and its BAM sibling had been updated in the same commit; `SampleReads`'s
*"Not `Clone` — it owns k files"* no longer follows from a `Vec<Arc<_>>`. All three are the first
lines a reader of each type sees. The last also records a real loss: non-`Clone` used to be a
compile error and is now a decision.

**M3: [open_bam.rs:457](../../../../src/ng/read/input/open_bam.rs#L457) — the module's own specs
still state the signatures this diff replaced**
**Category:** extras. **Confidence:** High. *Recorded, not applied (see §4.3).*

`spec/alignment_file.md:379` still gives `reads_in_region(&self) -> Result<RegionReads<'_>, …>` and
`:389` re-argues the `&self` receiver as load-bearing; `spec/sample_reads.md:180` still returns
`SampleRegionReads<'_>`. `arch/alignment_file.md:51-66` is stale too but independently, from
earlier work.

**M4: [open_bam.rs:294](../../../../src/ng/read/input/open_bam.rs#L294) — `open` returns a bare
`Self` that can no longer do the type's main job**
**Category:** smells. **Confidence:** High. *Recorded, not applied (see §4.1).*

**M5: [open_bam.rs:101](../../../../src/ng/read/input/open_bam.rs#L101) — one `Arc` is avoidable**
**Category:** idiomatic. **Confidence:** High. *Recorded, not applied (see §4.2).*

**M6: [mod.rs:1226](../../../../src/ng/read/input/mod.rs#L1226) — the outlives property covered the
`Single` arm only**
**Categories:** reliability, smells, extras (convergent, 3). **Confidence:** High. *Applied* as
`a_merged_region_stream_outlives_the_sample_reads_it_was_made_from`. A cohort sample is normally
several experiment files, so `Merged` — holding k readers and k tallies past the sample's death —
is the arm that matters in production and was the untested one.

**M7: [mod.rs:1278](../../../../src/ng/read/input/mod.rs#L1278) —
`a_held_region_stream_can_be_resumed_across_separate_borrows` did not exercise the half it was named
for**
**Category:** reliability. **Confidence:** High. *Applied.* `next_qname`'s `&SampleReads` was
unused, so nothing in the body could behave differently if the sample were not lent at all. Renamed
to `a_region_stream_can_be_stored_in_a_struct_without_a_lifetime` — what it does prove — with the
unused parameter's purpose spelled out so nobody "cleans it up"; and the falsifiable property it
was reaching for added beside it as `a_second_query_does_not_disturb_a_stream_already_held`, which
catches a pool handing out a reader on loan or two live queries sharing a cursor.

**M8: [open_bam.rs:1457](../../../../src/ng/read/input/open_bam.rs#L1457) — no `Send` assertion for
the types the change made holdable**
**Categories:** unsafe_concurrency, reliability, naming, idiomatic (convergent, 4). **Confidence:**
High. *Applied* as `a_sample_region_stream_is_send_in_both_arms`.

**M9: [open_bam.rs:686](../../../../src/ng/read/input/open_bam.rs#L686) — the manual `Debug`
hand-picks fields without destructuring**
**Category:** refactor_safety. **Confidence:** High. *Applied.* `resolution`, one of the two changed
fields, flowed through with no compile-time check. Now an exhaustive `let Self { … }` destructure,
so a field added, removed or retyped is a compile error there rather than a silent omission.

**M10: [mod.rs:44](../../../../src/ng/read/input/mod.rs#L44) — a peer→pipeline-stage back-reference**
**Category:** module_structure. **Confidence:** High. **Pre-existing, not introduced here.**
`crate::pop_var_caller::common::format_md5_hex` is imported by ng read input and by
`ng/reference_info.rs`. Suggested home: `src/fasta/`. *Recorded as out-of-scope.*

**M11: [open_bam.rs:1250](../../../../src/ng/read/input/open_bam.rs#L1250) — the new test
re-implemented `drain` verbatim**
**Category:** smells. **Confidence:** High. *Applied* — extracted `collect_reads`, so "the same
reads" now means literally the same code path on both sides of the comparison.

#### Nits

Applied: the unused-parameter comment (M7). Not applied, recorded here: the same §2.2 rationale
paragraph is now written seven times and one canonical statement with cross-references would keep
the copies from drifting the way M2's three already had; `Ordering::Relaxed` on `readers_opened`
lacks the justification comment the checklist wants (the ordering is correct — an independent
counter read after a scope join); the "two mutexes, never held together" invariant is true but
unstated; `region_query.rs`'s test helper `fixture()` is a bare noun; three `Arc::clone`s per query
where a source reaching through the file would need one — **accepted deliberately**, because
`region_query.rs` currently knows nothing about `AlignmentFile` and its unit tests build both
sources from a bare header and a bare resolution; taking the file would invert that dependency to
save two atomic increments on a path that then decodes thousands of records.

### 7. Out of scope observations

- **M10** above (`format_md5_hex`), pre-existing.
- `arch/alignment_file.md:51-66` describes fields that no longer exist (`path: PathBuf`,
  `header: sam::Header`, a `sample_name`) — drift from earlier work, not this diff's.
- `SampleLocusObservationsIterator`'s field drop order — see M1's verification pointer. Plan 3's.

### 8. Missing tests to add now

All four the review proposed were applied, plus the mutation probe that validated the first:

| test | input class | bug it catches |
|---|---|---|
| `a_stream_outliving_every_other_handle_still_banks_its_reader_and_tally` | stream is the last owner, file observable through a retained `Arc` | any weakening of `RegionReads::drop` on the detached path; a future non-owning handle (`Weak`) — `Arc::strong_count` pins the ownership claim |
| `a_merged_region_stream_outlives_the_sample_reads_it_was_made_from` | k = 2, the `Merged` arm, drained after `drop(sample)` | a merge keeping per-call state tied to the sample; an interleave regression visible only once the streams outlive their owner |
| `a_second_query_does_not_disturb_a_stream_already_held` | two live streams over one `SampleReads`, interleaved | a pool handing out a reader on loan; a cursor or scratch shared between live queries |
| `a_sample_region_stream_is_send_in_both_arms` | compile-time | an `Rc` or non-`Sync` field entering the merge, which would otherwise surface at the first parallel call site a milestone away |

Not added: a CRAM-flavoured detached test (reliability's #4). The CRAM path shares the same
`Arc<AlignmentFile>` and `Arc<sam::Header>` machinery, and the `DecodedContainer` travels with the
pooled reader rather than with the opener, so the BAM detached test plus `t8` cover the mechanism;
recorded here rather than written.

### 9. What's good

- **The oracle was chosen before the change and honoured**: 2487 → 2489 with the parity tests
  unmodified is a stronger statement about a representational change than any new test could make,
  and `t5`'s linear-scan oracle was deliberately *not* pulled onto the `Arc` path, so subject and
  oracle stayed independent.
- **The `&Arc<Self>` receiver survived three independent challenges** and is documented with its
  cost (the narrowing) rather than only its benefit.
- **`BorrowedReader` was left borrowed** — a smaller change than the plan asked for, correctly
  reasoned and self-reported.
- **Every `Arc::clone` is spelled `Arc::clone(&x)`, never `.clone()`**, so the sharing is visible at
  every site.

### 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --all --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --all-features` (expect 2493 lib tests)
- Mutation probe (manual): discard `RegionReads::drop`'s body and confirm
  `a_stream_outliving_every_other_handle_still_banks_its_reader_and_tally` fails.

### 11. Author response

- **M1** — fixed in the fixes commit; mutation-verified.
- **M2, M6, M7, M8, M9, M11** — fixed in the fixes commit.
- **M3, M4, M5** — deferred to Checkpoint A; they are owner decisions (§4), recorded in
  PROJECT_STATUS *Open*.
- **M10** and the §7 items — pre-existing or later plans'; recorded, not touched.

Per-category audit trail: `tmp/review_2026-07-28_ng-owned-region-stream/` (gitignored).
