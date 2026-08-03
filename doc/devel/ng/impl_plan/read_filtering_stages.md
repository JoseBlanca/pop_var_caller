# ng — read filtering in stages: implementation plan

**Status:** draft, 2026-08-03. Build order for dividing step 1 into two filters and a
conversion, and giving the loop to the cursor. The design is settled in
[`../spec/read_filtering_stages.md`](../spec/read_filtering_stages.md) (what and why) and
[`../arch/read_filtering_stages.md`](../arch/read_filtering_stages.md) (types and interfaces);
both of the spec's open questions were resolved by the owner on 2026-08-03. **This turns that
design into order and decides nothing.** If a step turns out to need a decision, stop and take it
back to the spec.

Follows [`alignment_cursor.md`](alignment_cursor.md), whose Milestone F left the read path with
one way to read a file and the two vestiges this plan removes.

---

## Scope

**In.** The renames that make the raw read say what it is; moving it and its conversion into
`read/aligned_read.rs`; replacing the per-contig fetch loop with a contig-table comparison;
moving the filtering loop out of `ReadFilter` and into `AlignmentCursor`; deleting `ReadFilter`,
`FilterState`, the `RecordSource` trait and its two test doubles, `with_validated_contigs` and
`ReadFilterBuffers`; `AlignmentCursor::reset_counts`; and the three tests the design owes.

**Out (later plans).**

- **The file owning its reference and building accessors itself** — spec §10. It would close the
  last gap in the contig check and drop a type parameter from both generators; it is a change to
  `AlignmentFile`, `SampleReads::open` and every tool, with its own reasons.
- **A whole-file read path** — spec §10. Filtering without a reference is not a way to *read* a
  whole file, which ng does not have.
- **`read/filtering/` as a folder** — spec §10. This change makes the file smaller; revisit only
  if that stops being true.
- **Any change to which filters run, their thresholds or their order** — that is
  [`read_filtering.md`](../spec/read_filtering.md)'s, and it is unchanged.

## Principles (how the order was chosen)

- **A rename never travels with a behaviour change.** Milestone A is renames only, byte-identity
  proven, so that when B and C change behaviour the diff is the behaviour.
- **The cheap independent win first.** The contig check (B) does not depend on the loop moving
  and pays for itself immediately, so it lands before the large step rather than behind it.
- **Isolate the step whose failure is silent.** C2 moves the tally. A wrong tally changes no
  output and no dump — the four acceptance dumps cannot see it — so it lands as its own commit
  with its own oracle, never bundled.
- **Build the replacement before deleting what it replaces.** C1 gives the in-memory reader a
  scripted error *before* C3 deletes the test doubles that raise fatal errors today. Otherwise
  three error-path tests vanish quietly, which is exactly how this branch has lost coverage
  before.
- **Verify against ground truth.** Every step is proven by output identity on real data — the
  four dumps and the walk probe — not by the unit suite, which cannot see a conversion hoisted
  above the flag checks.
- **Incremental, with pauses.** One milestone, then stop.

## Preconditions (already in place)

- **The design is settled.** Spec §9 records both questions resolved: the cursor holds the
  reference bases and the buffer, both filters stay plain functions; the contig check becomes a
  table comparison.
- **The suite is green at `72c6089`:** `cargo test --lib ng::` gives **1,538 passed**;
  `cargo test --lib` gives **2,837**.
- **The four acceptance dumps and the walk probe are reproducible**, and every step below is
  checked against them: `ng_generic_loci_dump` / `ng_ssr_loci_dump` on HG002 chromosome 21
  (251,792 and 4,406 lines) and on tomato `SL4.0ch01` (1,718,914 and 11,945), and
  `ng_generic_walk_probe` on chromosome 21 printing
  `loci=236081 observations=251786 reads_admitted=54709`.
- **The open gate reconciles the file against the reference**, which is what B leans on:
  `AlignmentFile::open` compares the `@SQ` list to the reference's contig table with
  `ContigList::first_disagreement` ([`open_bam.rs:206`](../../../../src/ng/read/input/open_bam.rs#L206)).
- **Every accessor implements `ContigTable`** — `InMemoryRefSeq`, `ResidentRefSeq`,
  `WindowedRefSeq` and the three test spies ([`ref_seq.rs:257`](../../../../src/ng/ref_seq.rs#L257)).
- **The validation gate is not the default one.** `cargo test --release` is red on a clean tree
  (four tests assert on `debug_assert!` messages); `cargo test --all-targets` aborts on a
  pre-existing panic in `benches/psp_writer_perf.rs:386`; `cargo doc` has 12 pre-existing
  unresolved links. Use `cargo test --lib`, `cargo test --examples`, and
  `cargo clippy --all-targets --all-features -- -D warnings`, in debug.

---

## The steps

### Milestone A — the renames (no behaviour change)

- ✅ **A1.** `RawRecord` → `RawAlignedRead` and `NoodlesRawRecord` → `NoodlesRawAlignedRead`,
  **moved to `read/aligned_read.rs`** beside `AlignedRead` and the conversion. Its doc gains the
  fact that an unmapped read is one of these. *Depends:* —. *Source:* spec §6, arch §2, §3.1.
- ✅ **A2.** `RecordReader` → `AlignedReadsReader`, its three arms to
  `BamAlignedReadsReader` / `CramAlignedReadsReader` / `InMemoryAlignedReadsReader`, and the
  module `record_reader/` → `aligned_reads_reader/`. Each reader's doc must now state that what
  it yields is undecoded, because the name no longer says so. *Depends:* A1. *Source:* spec §6.
- ✅ **A3.** `RegionRecords` → `RegionRawAlignedReads` (file follows the type),
  `DecodedContainer::fill_record` → `fill_raw_read`, `RecordIndex` → **`PackedReadEntry`**.
  *Depends:* A2. *Source:* spec §6, arch §2. — *`RecordIndex`'s new name was `RawReadIndex` when
  this step landed; the owner revised it at Checkpoint A on the review's argument (arch §2 has
  the reasoning).*

> **Checkpoint A:** nothing behaves differently. The four dumps are byte-identical, the probe
> prints the anchor, the suite is unchanged in count. Pause for review.

### Milestone B — the contig check becomes a comparison

- ✅ **B1.** `AlignmentFile::cursor` compares the accessor's contig table against the file's with
  `ContigList::first_disagreement`, and the per-contig fetch loop stops running on the cursor
  path. Adds `+ ContigTable` to `cursor`'s `R`, which propagates to `SampleReads::cursor` and
  both generators' signatures. **Its own commit** — a check that never fires looks exactly like a
  check that works, so it ships with a test that a mismatched accessor is refused, and that test
  is mutation-verified (disable the comparison; it must fail). *Depends:* A3. *Source:* spec §9
  Q2, arch §6.
- ✅ **B2.** `DecodedContainer::fill_raw_read` takes `&mut NoodlesRawAlignedRead` instead of
  `&mut RecordBuf`, and sets **both** its fields — the record and the read group — so the CRAM
  arm stops stamping the group on the line after the call. **Added by the owner at Checkpoint A**
  (2026-08-03), from A3's review: the name says it fills a raw aligned read and it fills half of
  one, which is the record-versus-read confusion Milestone A removed everywhere else. It is a
  signature change rather than a rename, so it could not travel with A3. Small and behaviour-free
  — the container already holds both halves — but it moves *where* CRAM's read-group exception
  lives, so it lands as its own commit and is proven by the four dumps like everything else.
  *Depends:* A3, and nothing in B1. *Source:* arch §2.

> **Checkpoint B:** the dumps and the anchor unchanged, and the walk probe's `seconds` measured
> before and after. The estimate says ~130 ms per cursor comes off; it has never been measured,
> so record what actually happens rather than the estimate. Pause for review.

### Milestone C — the loop moves into the cursor

- ☐ **C1.** `InMemoryAlignedReadsReader` gains a scripted error — a read position at which it
  returns `Err` instead of a record. Re-point the three fatal-error tests
  (`read_filter_source_read_error_is_fatal`, `read_filter_decode_error_is_fatal`,
  `read_filter_reference_error_mid_stream_is_fatal`) through it, so they run the real chain
  rather than a double that bypasses two layers. **Before C3, which deletes the doubles.**
  *Depends:* A2. *Source:* spec §6, §8.
- ☐ **C2. The cursor takes over the loop — its own commit, do not bundle.** `AlignmentCursor`
  gains the record buffer, the reference, the fetch scratch buffer, the config, the tally and a
  `failed` flag, and calls the two filters and the conversion itself. `ReadFilter`, `FilterState`,
  `restart_after_end_of_input`, `has_failed`, `source_mut`, `with_validated_contigs` and
  `ReadFilterBuffers` are all deleted. **The tally is the silent surface** — a wrong fold changes
  no output and no dump — so its oracle is `a_walk_charges_every_drop_reason_by_hand_count`, green
  before *and* after, plus the `other_sample` rider still landing on the first entry.
  *Depends:* B1, C1. *Source:* spec §5, §7; arch §3.4.
- ☐ **C3.** Delete the `RecordSource` trait and its two doubles; `RegionRawAlignedReads`'s trait
  implementation becomes inherent methods. Nothing generic consumes it once C2 lands.
  *Depends:* C2, C1. *Source:* spec §6, arch §3.3.
- ☐ **C4.** `AlignmentCursor::reset_counts`, with a test that a fresh window starts empty and that
  nothing else on the cursor moves. Small and additive; kept out of C2 so that step stays about
  the loop. *Depends:* C2. *Source:* spec §7, arch §3.4.

> **Checkpoint C:** the dumps and the anchor unchanged, `filtering.rs` holding only the
> keep-or-drop rules and their thresholds, and `cargo clippy --all-targets --all-features -- -D
> warnings` clean. Pause for review.

### Milestone D — the two tests output identity cannot see

Both are **work** properties: hoisting the conversion above the flag checks changes no output, so
nothing already in the tree would notice.

- ☐ **D1.** The first filter runs **with no reference at all** — construct it and drive it
  without a `RawRefSeq` in scope. Untested, spec §5's capability quietly stops being true.
  *Depends:* C2. *Source:* spec §8.
- ☐ **D2.** The conversion is asked for nothing when a read fails the first filter. **Use a read
  that would fail to convert:** unmapped, with no alignment start. Filter #5 drops it before the
  conversion; if the conversion were hoisted it would raise a fatal decode error instead. So the
  assertion is a clean drop charged to `unmapped` and no error — and no test double is needed.
  *Depends:* C2. *Source:* spec §8.

> **Checkpoint D:** both new tests mutation-verified — hoist the conversion above the first
> filter and D2 must fail; give the first filter a reference requirement and D1 must not compile.
> Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the four dumps byte-identical, the probe's anchor exact, the suite count unchanged — a rename that changes a number is not a rename |
| B | the same, plus a mutation-verified refusal of a mismatched accessor, plus the probe's `seconds` recorded before and after |
| C | the same, plus `a_walk_charges_every_drop_reason_by_hand_count` green before and after C2, plus clippy clean once the deletions land |
| D | each new test failing under the mutation it names |

**The oracle is this code before the change.** There is no second implementation to differ from,
so every step is measured against the four dumps and the walk probe — `ng_generic_loci_dump` and
`ng_ssr_loci_dump` on HG002 chromosome 21 and tomato `SL4.0ch01`, and
`ng_generic_walk_probe`'s `loci=236081 observations=251786 reads_admitted=54709`.

**Two things no unit test can prove, so they are checked by hand at every checkpoint:** the dumps'
byte-identity on real data, and — at B — what the contig check actually costs, because the only
number anyone has for it is arithmetic on a micro-measurement.

## Out of scope (next plans)

- **The file owning its reference**, so a caller cannot hand a cursor an accessor over a
  different FASTA at all. Spec §10 states the shape; it drops the accessor factory from
  `SampleReads::cursor` and a type parameter from both generators, and it is the natural
  companion to the parallel fan-out rather than to this change.
- **The parallel fan-out** — one cursor per worker. Unchanged by this plan, and easier after it:
  the cursor's state is all in one type instead of split across it and a filter.
