# The hidden-duplication filter — B1: the spill entry, and a codec that keeps an absence absent

**Date:** 2026-09-06
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone B, step B1
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.4, §3.7
**Branch:** `ng-paralog-filter`

## The answer

**A called record's finished VCF line, and the four numbers per sample the scorer reads, now
survive a round trip through bytes and back with every bit intact — including the `NaN` that
means *this sample has no usable window here*.**

The filter's cut cannot be resolved until every record has been called, so ng cannot write the
VCF as it calls. Pass one parks each finished record instead: the line the encoder produced,
the three fields the writer's ordering check reads, and per sample the window's GC fraction and
mean depth plus the reads behind the two alleles. This step is that parking format — the entry
type and the code that turns it into bytes and back. Nothing writes one yet; the file's
lifecycle is B2 and the sink that fills it is C2.

**The absence is the point.** A sample with no window at a locus is carried as `NaN` in both
floats, and any step that turns it into a zero makes a sample with no evidence look like one
with average coverage, which then folds into the run's estimate of how common duplications are
(spec §6 trap 4). So both floats are written as bit patterns and read back as bit patterns, and
every test compares by bit pattern rather than by `==` — which is not fastidiousness: `NaN ==
NaN` is false, so an equality-based round-trip test would fail on the one input the file exists
to carry, and a numeric comparison would pass on a codec that zeroed every absence.

## What was added

| file | what it is |
|---|---|
| [`src/ng/run/paralog_filter/mod.rs`](../../../src/ng/run/paralog_filter/mod.rs) | 49 lines: the module's declaration, what the three passes are, and the stand-in `WindowCoverage` |
| [`src/ng/run/paralog_filter/spill.rs`](../../../src/ng/run/paralog_filter/spill.rs) | 581 lines: `SpillEntry`, `SpilledSample`, `SpillError`, `SpillWriter`, `SpillReader`, and the codec |
| [`src/ng/run/paralog_filter/spill/tests.rs`](../../../src/ng/run/paralog_filter/spill/tests.rs) | 920 lines: 29 tests |

`src/ng/run/mod.rs` gains `pub mod paralog_filter;`. Nothing else in the tree is touched.

The layout is spec §3.4's, field for field, and
`the_encoding_matches_the_layout_the_spec_fixes` writes the whole of a small entry's encoding out
as a byte list, so the field order and each field's width are pinned by a test rather than by
the module header alone. The byte list is also where the offsets used by the tests that corrupt
one byte are read off.

## Assumptions and recorded deviations

**1. `WindowCoverage` is declared here, and it is scheduled for deletion.** The type belongs to
the window coverage design, which puts it in `src/ng/window_coverage/` on branch
`ng-window-coverage` — neither that module nor its specification is on `main`, and neither is
this branch's to create. So its two fields are declared in `mod.rs` under a doc comment saying
so, and the codec `use`s it from there. The codec reads those two fields and nothing else, so
the swap at the rebase is an import and a deletion. **The declaration is in `mod.rs` rather than
in the codec file** because a type scheduled for deletion should not also have two public paths,
and the crate already carries an unrelated `WindowCoverage` in `sample_summary::coverage`.

**2. There is no per-entry length frame.** Production's spill puts a `u32` length in front of
every record ([`spill.rs`](../../../src/var_calling/paralog_filter/spill.rs), `write_frame`), and
§7's reuse map says the *shape* is kept — length-framed, written once, read twice, deleted. But
§3.4's layout block, which is what the plan cites as the code shape, shows a length prefix on
the **line** and on nothing else, and every other field is fixed-width or self-delimiting. So an
entry ends where the next begins, and the reader detects a clean end of stream by finding no
byte where the next entry's first would be. The cost of the choice is that a corrupt entry
cannot be skipped over to reach the next one; the file is written and read inside one process,
so there is nothing to recover to.

**3. The two floats are four fixed little-endian bytes, not varints.** §3.4 writes `u32` for the
two float fields where every integer field beside them says `varint`, and an `f32`'s bit pattern
is not small: `NaN` is `0x7FC00000`, which LEB128 spells in five bytes. Fixed is both shorter and
what the contrast in the spec's own table says.

**4. The reader takes a `BufRead`, not a `Read`.** Variable-length integers are read one byte at
a time, because a varint's length is only known from its bytes. Buffered that is a bounds check
per byte; unbuffered it would be a read syscall per byte. The writer stays on plain `Write`.

**5. Four failure modes the spec does not name are refused rather than absorbed**, all of them
corruption a well-formed spill cannot contain: a flag byte that is neither `0` nor `1`
(`NotABoolean`); a decoded count or index too large for the field it belongs to (`OutOfRange`);
a ten-byte variable-length integer holding more than a `u64`, whose high bits the psp primitive
drops without complaint (`OverlongVarint`); and a record marked as **both** a repeat tract and a
biallelic SNP (`TractMarkedAsABiallelicSnp`), which spec §3.2 excludes by scoring every
non-SNP record — repeat tracts among them — on coverage alone. The alternative in each case —
`byte != 0`, a wrapping cast, six bits quietly discarded, a tract handed the SNP allele term —
turns a corrupt file into a plausible one. The last is refused by the writer as well as the
reader, so a producer cannot create the pair in the first place.

**6. A decoded sample count reserves at most 8,192 entries up front, and a decoded line length
is refused past 16 MiB.** Both are bounds against a corrupt number, and they are bounded
differently because the two fields fail differently. A count reserves before anything is read, so
four billion would become about a 64 GB allocation: `MAX_SAMPLES_RESERVED_UP_FRONT` caps the
reservation at 8,192 — 128 KiB, against spec §4's largest cohort of three thousand — and a larger
run decodes correctly and pays `Vec` growth. A length is *read*, and `read_to_end` appends every
byte the file still holds before the short-read check can fire, so on a spill the size of a
genome's VCF the bound was the rest of the file: `MAX_LINE_BYTES` refuses a length past 16 MiB
before a byte of it is read. A line is the eight fixed columns plus one genotype column per
sample, so 16 MiB leaves room for five times spec §4's cohort at a hundred bytes each.

**7. The spec's `SpillEntry` is copied field for field, with no `PartialEq`.** Deriving one over
a struct holding `NaN` would give a type whose equality is false against itself; the tests use
`holds_the_same_bits` instead, which destructures both sides and compares the floats by bit
pattern.

**8. The spec's §3.7 type block is kept as it stands, and the review argued for changing it.**
Two review agents proposed replacing `SpilledSample`'s `alt_reads` and `SpillEntry`'s
`is_biallelic_snp` with a sum type, so that a scorer reading an alternative count on a record
that has no single alternative would not compile. That changes the spec's type sketch, which this
loop may not do; **it is raised at Checkpoint B.** What landed instead is the half that needs no
ruling — the impossible combination is refused by both sides.

## Tests

**29 under `ng::run::paralog_filter::`**, all of them ng's own — nothing here is a copy. Twelve
of the twenty-nine came out of this step's review; see the
[fixes applied](../reviews/fixes_applied_ng_paralog_filter_b1_2026-09-06.md).

**The comparator has its own negative test, and that is the review's Blocker.** Every round trip
asserts one thing — that the entry came back holding the same bits — and hard-wiring
that comparison to `true` used to leave all seventeen tests green, so four of the nine were
asserting nothing but *the decode did not panic*.
`holds_the_same_bits_separates_entries_that_differ_in_any_one_field` now separates eleven pairs
differing in one field each, including **an absent window against a zeroed one**, which is spec
§6 trap 4 stated as a property of the comparator rather than of the codec. It is also
destructured on both sides, so a field either struct gains has to be answered for there or it
drops out of every round trip at once.

**Ten round trips**, each asserting the entry came back bit for bit and left nothing behind
it: a biallelic SNP with one covered sample and one with no window; a `NaN` with a payload of its
own; a record with **no samples at all**; a record whose `FILTER` is already `EMNoConv` and one
whose `INFO` is `.`; a repeat tract; an empty line; every field at the edge of its range, `-0.0`
and both infinities among them; **a 300-byte line** and **a thousand-sample cohort**, which are
the two places a length prefix crosses 128 and every fixture before them sat below it; and a
three-record stream, which comes back in the order it went in.

**Seventeen more on the file's shape and its refusals**: the encoding written out byte by byte;
an empty file yielding no entries; a file cut at every one of the nineteen offsets inside a
record **naming the field the bytes ran out in**, which pins all eleven of the decoder's labels
in one test; a failure inside the samples naming *which* sample; both flag bytes set to `2`; a
record marked as both a repeat tract and a biallelic SNP, refused by the writer and by the
reader; a contig of 2³²; a read count of 2³²; an eleven-byte varint and a ten-byte one holding
more than a `u64`; a line length past the ceiling; a sample count of four billion; a reader that
stops after a decode error; and a sink that refuses every byte and every flush.

**One property**, `any_stream_of_entries_comes_back_bit_for_bit`, over generated entries with
lines of 0–400 bytes and cohorts of 0–300 samples — both crossing the varint boundary — and
floats drawn from `NaN`, a `NaN` with a payload, `±0.0`, both infinities and any bit pattern at
all.

**Six mutations were run against the codec, each killed, each reverted from a byte-compared
backup.** Run as `cargo test --all-features --lib "ng::run::paralog_filter::"` in the container:

| mutation | what happened |
|---|---|
| the two floats written **and** read in the opposite order — a codec that agrees with itself and not with §3.4 | 15 of 17 passed; the byte-list test and `a_sample_count_larger_than_the_file…` failed |
| a decoded `NaN` turned into `0.0` — spec §6 trap 4, exactly | 7 of 17 failed, `a_snp_with_a_sample_that_has_no_window…` among them |
| the line's length check deleted, so a short read is accepted | 1 failed: the test that names the field |
| the comparator hard-wired to `return true` | 1 of 29 fails — before the review's fix, **all 17 passed** |
| a field added to `SpilledSample`, encoded correctly and decoded as `0`, with the fixtures and the pinned byte list updated as a coder resolving the errors would | **does not compile** — two `E0027` inside the comparator; before the fix, all 17 passed |
| the line's length read as one byte instead of a varint | 3 of 29 fail — before the review's fix, **all 17 passed** |

**The first is why the byte-list test is there**: the nine round trips are blind to a layout that
is wrong in the same way on both sides, and so are six of the other tests — fifteen in all.

**The third is why the truncated-line test names its field.** With the length check gone, the
decoder runs on into what should have been the sample count, finds nothing and reports a
truncated *sample count* — so a test asserting only that the file failed still passes; only one
that names the field catches it. That is now every cut, not one.

**The fifth is the review's point about the compiler.** The encoder's exhaustive destructure
forces a decoder to *mention* a new field; it does not force it to decode it correctly, and
before the comparator was destructured nothing did.

## Validation

Run in the dev container, on this tree.

    cargo test --all-features --lib "ng::run::paralog_filter::"
    → test result: ok. 29 passed; 0 failed; 0 ignored; 6373 filtered out

    cargo test --all-features --lib --bins --tests
    → 6387 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 3 errors, all needless_lifetimes in cohort_merge/build.rs and serial.rs

    cargo fmt --check
    → dirty on 9 files, none of them src/ng/run/paralog_filter/

The integration failure, the three lints and the nine files are `main`'s at `a33ada0f` and are
this branch's inherited gate (A1's report says how they were proved pre-existing). The gate — no
failure `main` does not already have, and the `fmt` and `clippy` failure sets no larger than
`main`'s — holds, and the lib count moves 6,358 → 6,387.

**`cargo fmt` reformats `main`'s nine files as a side effect and they were reverted with
`git checkout --` before staging**, so this step's diff is three new files and one line of
`src/ng/run/mod.rs`.

## Tradeoffs and follow-ups

- **`MAX_SAMPLES_RESERVED_UP_FRONT` is a guess at where a cohort stops**, not a limit: a run of
  more than 8,192 samples decodes correctly and pays one `Vec` growth per record. It exists only
  to keep a corrupt count from becoming an allocation.
- **Nothing yet tells the reader how many entries the file should hold**, so a spill that lost
  its tail on an entry boundary reads back as a complete, shorter spill — and pass three would
  then write a VCF short by its last records, with no panic and no message. The writer counts
  them; B2 owns the file's lifecycle and is where the number can be carried across.
- **The reader allocates two `Vec`s per entry** where the writer reuses one scratch buffer: at
  three thousand samples that is 48 KB a record, twice over. Whether it matters is what step D2's
  pass shares say, so the fix is not built blind.
- **The codec is not compressed.** Spec §3.4 defers that to "if it proves easy" — one writer and
  one reader — and this shape leaves it that way.
- **`WindowCoverage`'s deletion is the rebase's first task**, together with the import that
  replaces it.
- **Nothing writes or reads a spill file yet.** B2 gives it a path and a guard that unlinks it on
  every exit path; C2 gives it a producer.
