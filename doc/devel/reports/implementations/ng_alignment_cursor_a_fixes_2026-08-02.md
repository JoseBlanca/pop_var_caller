# ng — the alignment cursor, Milestone A: review fixes applied

*Companion to [the review](../reviews/ng_alignment_cursor_a_2026-08-02.md) and
[the implementation report](ng_alignment_cursor_a_2026-08-02.md). A section per step.*

---

## A1

**All findings applied**, with **two disputed on measurement** and one recorded as a stated
gap rather than a test. Nothing was deferred.

The suite went from 9 tests to 19, and four of the new ones failed when first written —
each failure a real property of the code or the fixture that had been assumed rather than
checked.

### The Blocker

**B1 — the typed-region walk could drop three quarters of its regions silently.**

Added `the_typed_walk_visits_every_region_the_stream_produced`, which counts the same typed
stream a second time without the walk in the way and compares three things: the region
count, the base pairs those regions cover, and the loci (one per base, since the fixture's
reads are pure Match and tile the contig).

**Verified by re-running the review's own mutation.** With
`walks.into_iter().flatten().take(1)`:

```
the walk saw 1 Generic regions; the stream produced 4
test result: FAILED. 18 passed; 1 failed
```

### The Majors

| finding | what was done |
|---|---|
| **M1** the discarded `stream.counts()` | `regions_in`, `regions_handled`, `loci_emitted` and **both unhandled counters** now live in `ProbeReport` and print. `the_walk_accounts_for_every_region_and_every_locus` asserts the spec §13.2 partition across three run shapes. |
| **M2** 3 of 16 printed keys pinned | `the_report_prints_exactly_the_documented_key_set` pins all twenty-one, in order. `each_printed_key_carries_its_own_counter` gives every field a distinct value so a transposition inside `render` cannot pass. The copy out of `PileupGeneratorCounts` is now one destructure, so a counter added upstream stops the build here. |
| **M3** untested contig filter | Hoisted to `contig_is_selected`, one body instead of two, with a unit test and an end-to-end test through **both** region sources. |
| **M4** the record-span knob | **Disputed on measurement — see below.** |
| **M5** untested fixture, false tiling claim | Three `assert!`s in `SyntheticSample::build`, a `stride()` method the tiling rule is stated in, and five tests including three `should_panic`. |
| **M6** the wrong `Cargo.toml` claim | Rewritten to what was measured. |
| **M7** the allocator guard's justification | Rewritten to say what it catches (a contradiction the operator typed) **and what it does not** (the silent case, because nothing in the output names the allocator). |

### Two findings disputed on measurement

**M4 — the record-span knob is unobservable on this fixture, and the proposed test asserted
an effect that does not occur.** The review suggested asserting that a 200 bp halo and the
default 5,000 bp one admit different reads. Measured at 500 bp regions, the two produce
**identical reports** — every counter, including `reads_admitted=425` and
`records_outside_region=1288`. That is not the knob failing to reach the generator: the walk
stops at the region's end and the read stream is lazy, so halo records are never pulled
unless something reaches back into the region, and nothing here does — these reads are 150 bp
of pure Match with no long deletion. Exercising the halo needs a fixture containing the thing
the halo exists for. **Recorded as a comment stating the gap**, because an assertion that
happens to hold would claim coverage this file does not have.

**M2's walk-side transposition cannot be caught here, and the code now says so.** Swapping
two assignments out of `PileupGeneratorCounts` stays green because seven of those eight
counters are zero on this fixture. The destructure stops a counter being *forgotten*; it
cannot stop two being *crossed*. Written down rather than left as implied coverage.

### Four assertions that failed when first written

Each one taught something the review and the implementation had both assumed:

1. **`regions_handled == regions_in` is false by design.** On the typed path, `4 of 7`
   regions reach a generator — the rest are repeat regions with no generator filled. The
   correct invariant is the three-way partition, which is why both unhandled counters were
   added rather than just the two the review named.
2. **An observation is a folded row, not a read.** At 10× the fixture yields **5,332**
   observations over 5,000 loci, not ~50,000: ten reads agreeing with the reference fold into
   one row. The assertion is now a band — surplus between 1 and the number of reads written —
   with the fold explained. The review's suggested `loci * 5` would have been wrong in the
   same way the original `> loci` was weak.
3. **The fixture does not tile at coverage 1.** 33 reads at a stride of 152 over 150 bp reads
   covers **4,950 of 5,000** positions. `build` now refuses any geometry whose stride exceeds
   its read length.
4. **The exhaustive `ProbeReport` literal earned its keep immediately** — adding the two
   unhandled counters failed to compile at exactly the test that pins the printed surface.

### One finding the orchestrator's own mutation testing added

After applying the fixes, the four mutations were re-run. Two still passed. One was M2's
transposition (above); the other was **loosening `contig_is_selected` to a prefix match** —
which my test cases missed, because none of them had a contig name *extending* the filter.
`chr21_KI270872v1_alt` is exactly that, and exactly what the anchor's reference is full of.
Added, and the mutation now fails.

### Validation

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test` (plain) — reaches the probe's **19 passed**.
- `cargo test --lib --tests --examples --all-features` — `2770 passed; 0 failed; 5 ignored`.
- `cargo bench --bench ng_generic_pileup_perf -- --test` — six cases, all `Success`.
- **The chr21 anchor, after every fix:** `loci=236081 observations=251786
  reads_admitted=54709` — unchanged. The new lines are additive:
  `regions_in=205875 regions_handled=102938 unhandled_not_implemented=102659
  unhandled_out_of_scope=278 loci_emitted=236081`, which partitions exactly and records for
  the first time that **half of chromosome 21's typed regions have no generator filled**.

### Not applied

Nothing. Two items are recorded in the review as out of scope for A1:

- `benches/psp_writer_perf.rs:386` panics under `cargo test --all-targets` (pre-existing,
  untouched by this branch, and the cause of a wrong measurement in this review).
- `examples/ng_generic_loci_dump.rs` holds a third copy of `write_fasta` / `write_bam`;
  converting it belongs to Milestone F1.

---

## A2

**All findings applied.** One was reached independently by all three agents; three mutations
that survived the first draft now fail.

### The Major, found three ways

**The docs said `contigs` is "the reference's list" in three places. It is the file's.** The
value is `contig_list(&header)`, reconciled against the reference by a gate that treats an
absent `M5` as a wildcard — so a file declaring digests against a `.fai`-only reference
passes, and what is stored is the file's claim.

**And the test named after the property could not fail.** Two agents mutated the field to be
literally what the docs claimed; `the_contig_list_is_the_reference_own_list_in_its_own_order`
passed under both the true and the false version, because its whole-list `assert_eq!` runs
through `ContigEntry`'s digest-wildcarding `PartialEq`.

Applied: the field doc, the accessor doc and all three tests now say what the value is —
the file's list, whose **names, lengths and order** the gate proves are the reference's, with
the digests explicitly the file's. The tests compare those three fields directly rather than
through `PartialEq`, and use a **distinct digest per contig** so an order or a projection
error cannot pass.

**The same sentence is in arch §2.1 and spec §8.** Not edited — raised at Checkpoint A.

### Three mutations that survived, and now do not

| mutation | caught by |
|---|---|
| populate `contigs` from the reference instead of the file | `the_contig_list_and_the_digest_list_index_alike`, `…_carries_the_digests_the_file_declared_not_the_reference_s` |
| build `sq_md5s` in reverse order | `the_contig_list_and_the_digest_list_index_alike` |
| `Debug`'s contig count `+ 999` | `the_debug_line_counts_the_contigs_the_file_actually_has` (new) |

### The error type

- `WrongChromosome` gained the **file path**. Two bare `ContigId`s are indices, and an index
  means nothing without naming the table it indexes; a run holds up to 320 cursors. The
  message is now `cursor on '/data/sample.bam' covers contig 20 but the region is on contig
  7`.
- `Io` → **`ReadRecord`**, and its message joins the module's house shape
  (`reading alignment file '…' failed`). The rule — name the failed operation, not the
  mechanism — is one every other error enum in the crate follows.
- `cursor` / `requested` → **`cursor_contig` / `requested_contig`**: two same-typed fields
  whose names did not tell them apart, in the variant whose whole job is telling them apart.

### The doc-comment claims

- `sq_md5s_by_file` **does not exist**; the real consumer is `SampleReads::assembly_inputs`.
- `(spec §4)` cited "Error model"; corrected to `alignment_file.md` §3.1, check 2.
- 35,228 now names its mode, because spec §11.5 says both figures are in circulation and B3
  has to choose between 34,633 and 35,228.
- The module's **rustdoc summary line** now says what the file contains. The index showed
  `cursor` as *"A reader that stays where it is."* while the file held an error enum.
- Two `unresolved link` rustdoc **errors** — intra-doc-link brackets on filesystem paths,
  under the crate's `broken_intra_doc_links = "deny"` — replaced with plain backticked paths.
  Confirmed gone: `cargo doc --no-deps --lib --all-features` reports no `cursor.rs` error.

### Applied with a stated reason not to go further

**`contigs` and `sq_md5s` are duplicated state, and my justification for keeping both was
false.** An agent deleted the field, derived the accessor, and got a green suite ~18 lines
shorter. But that route changes `check_assembly` — a `pub` function — and eight call sites,
which is more than the step that *introduced* the duplication should carry. The field doc now
states the real reason (the accessor lends a slice; `assembly_inputs` lends one from every
open file at once), records the measured alternative, and names it a follow-up.

### Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ng::` — **1477 passed; 0 failed; 2 ignored**.
- `cargo doc --no-deps --lib --all-features` — no `cursor.rs` errors; the 12 that remain are
  pre-existing in files this branch does not touch.

---

## A3

**All findings applied. Thirteen mutations were run against the first draft and six survived;
all six now fail.**

| mutation | now caught by |
|---|---|
| `begin_region` rewinds only when the previous region was drained | `a_reader_abandoned_part_way_through_still_rewinds` |
| the enum's `begin_region` swallows the call | `the_enum_forwards_every_contract_method_to_its_arm`, rewritten to reposition **after** a drain |
| the clone loses `alignment_start` | `the_whole_record_survives_the_clone` |
| the clone keeps only the name | `the_whole_record_survives_the_clone` |
| `read_group` cleared only for the first record | `every_record_of_a_pass_comes_out_with_no_read_group` |
| `header()` returns a different header of the same shape | `the_header_is_the_one_the_reader_was_built_with` |

The root cause of three of them was one helper: `drain` compared **name lists**, so anything
that preserved names passed. It now delegates to `drain_records`, and its doc says why names
alone are not enough.

### The contract's false claims

- **"reusing the buffer's allocations"** — true of a file arm, false of this one, and the
  difference is now stated with its mechanism (`RecordBuf` derives `Clone`, so `clone_from`
  is the default `*self = source.clone()`). The claim that a real arm has "the same cost
  shape" is deleted: it is the opposite.
- **"A reader holds only its position"** — now says "and (from Milestone C) the single record
  the sorted early stop consumes without yielding", which is what arch §1.3 requires and what
  C2 will implement.
- **"stated once so two arms cannot drift"** — the heading now says what the list is: a place
  to check an arm against, not a mechanism, with a note that what will actually hold the arms
  together is the oracle at spec §11.3 and a shared harness when the second arm lands.
- **"a record replayed from memory"** — spec §5's sentence is about **reads**, and a replayed
  read skips decode and filtering entirely. Corrected, with the distinction spelled out.

### Removed rather than fixed

**`other_sample_records`.** It returned a constant `0` justified by a rule that does not
apply: `RecordReader` is an enum with inherent methods and does not implement `RecordSource`.
Nothing calls it until `RegionRecords` at C1, which resolves read groups itself and will
answer for its own skipping.

### Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ng::` — **1488 passed; 0 failed; 2 ignored**.

---

## A4

**All findings applied. Five mutations were run against the first draft and none survived** —
unusual for this milestone, and it is because the step's own tests were built from a finding
rather than from an assumption. What the review corrected was an assertion that could not
fail and a doc comment aimed at the wrong half of the problem.

| finding | what was done |
|---|---|
| **M1** the fuse's clean-EOF half | `source_mut`'s doc now names the clean end of input as *the* case, carries the measured "2 reads then 0", and says plainly that the accessor is not enough for `move_to_region`. The remedy is left open as a Milestone B decision. |
| **Mi1** two byte-identical records | `named_fake`, and the closing assertion now compares read **names** — verified to fail under a no-op `rewind` on that assertion alone. |
| **Mi2** accumulation untested | four lines: after a rewind and a replay, `duplicate` and `kept` both reach 2. |
| **Mi3** `position` collides with the genomic meaning | renamed `records_consumed`, with a doc saying what it is not. |
| **Mi4** no fatal-error reposition test | `a_filter_stopped_by_a_fatal_error_stays_finished_after_a_reposition`. |

### Re-run after the fixes

| mutation | result |
|---|---|
| `rewind` is a no-op | **4 tests fail** |
| `records_consumed` always reports 0 | 1 fails |
| `next()` no longer fuses on clean EOF | 1 fails |

### Not applied — a stop-and-ask

The review's second half of M1 asks for edits to `spec/alignment_cursor.md` §3 and
`arch/alignment_cursor.md` §2.3, which claim the filter seam is *solved* by this accessor.
They are design documents; this loop does not edit them. **Raised at Checkpoint A.**

### Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ng::` — **1492 passed; 0 failed; 2 ignored**.
