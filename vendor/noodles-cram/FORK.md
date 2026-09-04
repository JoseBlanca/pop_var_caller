# What this copy has that crates.io does not

This is `noodles-cram` 0.93.0 — the same version `Cargo.toml` names for every other noodles
crate — extracted from the cargo registry and then changed. **This file is the whole list of
changes.** Anything not named here is upstream's, byte for byte, and the claim is checkable:

```
diff -r ~/.cargo/registry/src/index.crates.io-*/noodles-cram-0.93.0 vendor/noodles-cram
```

`[patch.crates-io]` in the workspace manifest is what points the build at this copy rather than
at the registry.

## Why a copy at all

ng reads a CRAM by driving `Slice::decode_blocks` and `Slice::records` directly, and the CRAM
read path is about a third of a single-sample calling run. Two of the costs found there are
inside noodles and reachable in no other way: neither the block decompressors nor the record's
field accessors take a caller-supplied buffer, and the record's fields are `pub(crate)`. See
`doc/devel/ng/research/cram_read_path_2026-09-04.md`.

## The changes

### 1. Per-codec counters, behind the `perf-counters` feature — measurement only

`src/perf.rs` (new), declared in `src/lib.rs`; `Block::decode` in
`src/io/reader/container/block.rs` split into a timed wrapper and a `decode_inner` carrying
the body it always had; the feature declared in this crate's `Cargo.toml`.

**With the feature off, `decode` is the upstream function under a different name and nothing
else changes** — no counter, no clock read, no atomic. With it on, each block decode records
its compression method, its compressed and inflated sizes, and its elapsed nanoseconds, so a
run can say which codec its decompression time went into. That question has no other answer:
a sampling profile shows `Block::decode` as one frame whatever method the block used.

This is a measuring instrument, not a change to what noodles does, and it is not the kind of
thing to send upstream.

### 2. rANS 4x8 order-1 skips the symbol contexts a block never uses — a decode-speed fix

`src/codecs/rans_4x8/decode/order_1.rs`:
`build_cumulative_frequencies_symbols_table` now takes the frequency table as well as the
cumulative one, and leaves a context's 4,096-entry lookup at the zeros the allocation already
holds when every frequency in that context is zero. `order_0.rs` and `order_1.rs` also carry
the feature-gated phase counters that show the cost.

Measured on 60 containers of a whole-genome tomato CRAM: block inflation 0.347 s → 0.268 s.
**This one is a plain defect and belongs upstream** — it needs no new API and changes no
behaviour a caller can observe.

### 3. A decoded record's fields, written into buffers the caller owns — a new API

`src/record/direct.rs` (new), declared in `src/record.rs`; `src/record/sequence.rs` gains
`iter_concrete`, and its `iter` module and that module's `Iter` are widened to `pub(crate)` so
the bases can be iterated without a `Box`.

`Record`'s only published accessors are `sam::alignment::Record`'s, and each returns a
`Box<dyn …>`: eight heap allocations per record before a byte is copied, whatever the caller
wants. The new methods — `write_bases_into`, `write_quality_scores_into`, `write_cigar_into`,
`name_bytes`, `read_group_index` — write into a `Vec` the caller supplies and reuses.

Measured: building ng's read from each record 0.355 s → 0.140 s, 2.5×, over 600,000 reads
whose every field hashes the same both ways.

**This one does not belong upstream as it stands.** It is an interface noodles has no other
caller for, and proposing it is a conversation about a `write_*_into` family across the
alignment-record traits, not a bug report.

### 4. Records whose auxiliary tags are never read — a new API, and a check that says when it is safe

`src/io/reader/container/slice.rs` gains `Slice::records_discarding_tags` beside `records`
(both now call one private `read_records`); `src/io/reader/container/slice/tag_streams.rs`
(new) decides whether the tags may be left unread; `slice/records.rs` gains a `TagPolicy` and
acts on it.

A caller that reads no tags — ng reads none, since a CRAM stores the read group as a number
rather than a tag — was paying for every one of them. Skipping them is not simply a matter of
not asking: CRAM decodes tags and data series out of the same streams, and a stream is a
cursor, so a skipped read that something else depended on leaves every value after it wrong
and says nothing. `tags_can_be_left_unread` decides the question from the compression header
before any record is read, and answers *no* unless every tag reads only external blocks that
no data series reads. When it answers no, the tags are decoded and their values discarded,
which is smaller but always safe.

Measured, whole-genome tomato: the read path 0.570 s → 0.454 s. On a GIAB human CRAM written
by a different aligner, 0.551 s → 0.505 s. Records identical in both.

**Not upstream as it stands**, for the same reason as change 3: it is an API with one caller.
The `tag_streams` check would be the interesting part of any upstream proposal, since it is
what makes the option safe to offer at all.
