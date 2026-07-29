# ng generic locus generator — prerequisites, Milestone A: an owned region stream

**Date:** 2026-07-28 · **Plan:**
[locus_generation_pileup_prerequisites.md](../../ng/impl_plan/locus_generation_pileup_prerequisites.md)
steps A1–A3 · **Spec:** [locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md) §2, §7
· **Arch:** [locus_generation_pileup.md](../../ng/arch/locus_generation_pileup.md) §2.2, §4 P1

Implementation report for Milestone A of plan 1 of 3. One report for the milestone rather than one
per step: the three steps are one change taken in three commits, and three reports would be three
copies of the same argument.

## 1. Plan

`src/ng/read/input/` must hand back a region stream that does **not** borrow the `SampleReads` it
was made from. `LocusGenerator::next_locus` lends `&SampleReads` per call and carries no lifetime
parameter, so a generator cannot hold a borrowed stream between calls — and the generic locus
generator, unlike the STR one, yields many loci per segment and therefore must. *Resumable +
borrowed stream + no lifetime on `Self`* is not expressible; the borrow is the cheapest of the three
to give up (arch §2.2).

Three borrows had to go. The arch doc's first draft named two; the third arrived with the
`ng-read-groups` merge.

## 2. Assumptions

None that changed direction. Two judgement calls are recorded in §6.

## 3. Changes made

### A1 — `Arc<sam::Header>` (`d1db8dd`)

- **[open_bam.rs](../../../../src/ng/read/input/open_bam.rs)** — `AlignmentFile.header` becomes
  `Arc<sam::Header>`, built once at `open` and cloned out per query.
- **[region_query.rs](../../../../src/ng/read/input/region_query.rs)** — `BamRegionSource` and
  `CramRegionSource` hold `Arc<sam::Header>` in place of `&'a sam::Header`.

An **independent** `Arc`, not a reference into an `Arc`'d file — that is what keeps the sources from
becoming self-referential once the file itself is shared (arch §2.2's own warning). Precedent in the
same struct: `CramRegionSource.entries: Arc<[crai::Record]>`.

### A2 — the file and the read-group resolution (`60acd15`)

- **`AlignmentFile.resolution: Arc<ReadGroupResolution>`**, and both region sources hold a clone.
  This is the third borrow. `Arc` rather than a clone of the value: `ReadGroupResolution::PerRecord`
  owns a `Box<[(Box<str>, RecordOwner)]>`, so cloning it would allocate a table per query on a path
  that runs ~10⁶ times.
- **`RegionReads` holds `Arc<AlignmentFile>`** and loses its lifetime parameter. It needs the file at
  `Drop` — to return the pooled reader and fold the query's tally — and that is precisely the point
  at which the file must outlive the borrow that made the stream.
- **`AlignmentFile::reads_in_region` takes `self: &Arc<Self>`.** One atomic increment per query,
  against a query that then decodes thousands of records.
- **`SampleReads.files: Vec<Arc<AlignmentFile>>`**; `RegionSource`, `MergedRegionReads` and
  `SampleRegionReads` all lose their lifetime parameter.

`SampleReads::reads_in_region` still takes `&self`. Only the *returned type* stops borrowing.

### A3 — the property, asserted (this commit)

Two tests in [mod.rs](../../../../src/ng/read/input/mod.rs), plus the unchanged suite as the oracle.

## 4. Tests added

| test | what it proves |
|---|---|
| `a_region_stream_outlives_the_sample_reads_it_was_made_from` | the `SampleReads` is **dropped before the first read is pulled** and the stream still yields exactly the control drain. Does not compile at all if the stream borrows (`cannot move out of sample because it is borrowed`); and the file must still exist at `Drop`, where the pooled reader goes back and the tally is folded in. |
| `a_held_region_stream_can_be_resumed_across_separate_borrows` | the shape that needs it: a **lifetime-free struct** holding the stream, with `begin`/`next_qname` each lent `&SampleReads` separately — `LocusGenerator::next_locus` in miniature. A compile-time anchor: a borrowed stream forces a lifetime parameter onto the struct and the test does not build. |

The first compares against a control drain rather than counting reads, because "the stream survives"
and "the stream survives and still yields the right reads" are different claims.

**What the milestone's oracle really is: nothing moved.** This is a representational change, so the
existing tests are the assertion. `cargo test --all-features` went 2487 → 2487 (A1, A2) → 2489 (A3,
the two new tests), with the BAM/CRAM parity oracle
(`t8_a_cram_yields_the_same_ordered_reads_as_the_same_bam`), the region-query-vs-linear-scan oracle
(`t5_…`) and the whole read-input suite **unmodified**. A moved read would be a behaviour change, not
a representation change.

## 5. Validation

Run in the container (`./scripts/dev.sh`), per commit:

- `cargo fmt --all --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — no diagnostics.
- `cargo test --all-features` — **2489 passed / 0 failed / 4 ignored** at A3 (2487 at A1 and A2).

Two standard commands are excepted by hand, red independently of this work and tracked under
PROJECT_STATUS *Standing project-wide items*: `cargo test --all-targets --all-features` (panics in
`benches/psp_writer_perf.rs`) and `cargo doc --no-deps` (11 unresolved intra-doc links).

## 6. Deviations from the plan

All three minor, all recorded rather than escalated. *(The third was added after the review, which
found it: the first draft of this section listed two.)*

1. **`BorrowedReader` keeps its borrow.** The plan lists it alongside `RegionReads` as becoming
   `Arc<AlignmentFile>`. It does not need to: it is created and `take()`n inside `reads_in_region`
   with nothing fallible between, so its lifetime never reaches the returned type. Leaving it
   borrowed avoids an `Arc` clone per query for a value that is consumed immediately. The intent —
   no borrow escapes — is unchanged.
2. **The receiver is `&Arc<Self>`, not `&self`.** The plan does not say how `RegionReads` comes by
   its `Arc`; this is the only way that does not force `AlignmentFile::open` to return an `Arc`. The
   cost is a real API narrowing — a bare `&AlignmentFile` can no longer query — which is documented
   on the method and absorbed by the test helpers (`opened_over` and `OpenedFixture.file` now hand
   out `Arc<AlignmentFile>`).
3. **`RegionSource` does not hold an `Arc<AlignmentFile>` either.** Plan A2 names it beside
   `BorrowedReader` and `RegionReads`; it carries the header and resolution `Arc`s instead, which
   serves the same end — no lifetime — while keeping `region_query.rs` ignorant of `AlignmentFile`,
   as its unit tests require (they build both sources from a bare header and a bare resolution).
   The cost is two atomic increments per query that reaching through the file would have saved.

## 7. Review, and what it changed

Reviewed the same day over the whole milestone diff:
[ng_locus_generation_pileup_prereq_a_2026-07-28.md](../reviews/ng_locus_generation_pileup_prereq_a_2026-07-28.md)
— 6 categories, **0 Blockers, 1 Major, 11 Minor**, verdict Approve-with-changes.

**The Major was mine and it was a real one.** §4's claim above — that the outlives test "reaches"
the drop path — was false: after `drop(sample)` nothing can observe the tally, so gutting
`RegionReads::drop` left the test green. Fixed by asserting the property where it is observable, at
the `AlignmentFile` level through a deliberately retained second handle, and **mutation-verified**:
with `RegionReads::drop`'s body discarded, `a_stream_outliving_every_other_handle_still_banks_its_reader_and_tally`
fails while both detached tests stay green — exactly the review's claim. Six more findings applied
(three stale doc comments, the `Merged` arm, the interleaved-queries property, the `Send` anchor,
an exhaustive destructure in the manual `Debug`, a de-duplicated test helper). Suite 2489 → 2493.

## 8. Checkpoint A — the three decisions, and what the owner chose

All three were raised at the checkpoint rather than decided in code, because each reached past the
step's remit. **The owner took all three (2026-07-28)**, and each landed as its own commit.

1. **`AlignmentFile::open` returns `Arc<Self>`** — *"you can change to pub, no problem"*. The
   `Arc`-ness was an invariant of *using* the type that ten call sites each had to remember by
   chaining `.map(Arc::new)`; stating it once in the constructor makes it un-skippable and deletes
   all ten. The four call sites that only inspect the error are unaffected, and no `Arc` is
   allocated on the failure path.
2. **The share moved inside `ReadGroupResolution`, and the wrapper is gone** — *"I don't like the
   idea of an `Arc` in a struct that will have millions of objects, you might consider modifying the
   type."* Right, and it was the sharper reading of the finding: a region source is built **per
   query**, ~10⁶ times a run, so `Arc<ReadGroupResolution>` charged a pointer chase and an atomic
   pair to every one of them — including the overwhelmingly common `Sole` case, where there is
   nothing to share at all. `PerRecord(Box<[…]>)` became `PerRecord(Arc<[…]>)` (settled at open,
   read but never written), which makes the enum cheap to clone by value; `AlignmentFile.resolution`
   and both region sources now hold it **directly, with no wrapper**.
   `Arc<sam::Header>` stays, and the contrast is the justification: it is noodles' type, not ours to
   reshape, and it has no cheap-clone form — which is the case an `Arc` is actually for.
3. **The two spec fold-ins landed** — *"you can edit the docs, in this case."*
   [`spec/alignment_file.md`](../../ng/spec/alignment_file.md) §3.4 now gives the `&Arc<Self>`
   receiver and the `Arc<AlignmentFile>` return, with a dated fold-in note preserving the original
   *shared, not `&mut`* argument (which is unchanged) and adding why the receiver is now shared
   **ownership**; [`spec/sample_reads.md`](../../ng/spec/sample_reads.md) §3.4 records that only the
   *returned type* stopped borrowing, plus the tally caveat a caller has to know.

**And the hazard for plan 3 is now a marker in the code**, at the owner's instruction:
`SampleLocusObservationsIterator` ([locus_generation/mod.rs](../../../../src/ng/locus_generation/mod.rs))
declares `reads: SampleReads` **before** `generators`, and Rust drops fields in declaration order —
so once a generator holds a region stream, the sample dies first and that stream's step-1 tally
becomes unobservable at drop. A `FIXME(pileup-generator)` on the field states the mechanism, the two
cheap fixes, the test shape it needs, and **that the comment is to be deleted once it is fixed**.
Latent today; live the moment Milestone A's capability is used.

Milestone B is next, and B1 — the `ReadCoverage` reshape — is a silent-failure step that gets its
own commit with the STR dump green before and after.
