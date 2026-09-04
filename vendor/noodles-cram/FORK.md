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

### 2. rANS 4x8 builds one packed lookup per used context, not three arrays over all 256

`src/codecs/rans_4x8/decode/table.rs` (new), and the decode loops in `decode/order_0.rs` and
`decode/order_1.rs` rewritten around it. `order_0.rs` also loses `build_cumulative_frequencies`
and `build_cumulative_frequencies_symbols_table`, which nothing needs any more.

Decoding a byte needs three numbers keyed on the same slot of the same context — which symbol
owns the slot, that symbol's frequency, and where its range starts — and they lived in three
arrays, all built over all 256 contexts whether the block used them or not. They are now one
32-bit word per slot (range start, frequency less one, symbol), and only the contexts whose
frequencies are non-zero get a table at all. A block of DNA bases uses about five of the 256
and a block of quality scores a few dozen; the rest can never be reached, because the decoder
indexes context *i* only after emitting symbol *i*.

The frequency is stored less one because it does not otherwise fit: frequencies are normalised
to sum to 4,096, so a single-symbol context has a frequency of 4,096 and needs thirteen bits.
Its range then starts at zero and every other frequency is at most 4,095, so subtracting one
makes symbol, frequency and range start fit in exactly 32 bits with nothing spare.

Measured over 60 containers of a whole-genome tomato CRAM, 1,351 rANS blocks inflating 74.5 MB:

| | building the tables | decoding the symbols | block inflation, all codecs |
|---|---:|---:|---:|
| before | 0.161 s | 0.196 s | 0.347 s |
| after | **0.015 s** | 0.183 s | **0.227 s** |

**Almost all of it is the table building, and that is the finding.** Three further attempts at
the decode loop each moved it by 2 % or less and were dropped: the packed lookup itself
(0.196 → 0.186 s), removing the inner bounds check by giving each context a fixed-size row
(0.186 s, unmoved), and renormalising without an `io::Result` per byte (0.186 → 0.183 s). At
352 MB/s the loop is bound by the dependency chain — each byte's table row is chosen by the
byte before it — and the format's four interleaved states are all the parallelism there is.
The packed table and the fixed-size rows are kept anyway, because they are what makes building
only the used contexts natural; the renormalise change was reverted for earning nothing.

**Upstream: yes.** It changes no API and no observable behaviour, and the 223 tests include two
new ones for the packing, one of them for the single-symbol context that does not fit in twelve
bits.

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

### 5. Decoding against a window of the reference instead of a whole contig — a new API

`src/io/reader/container/slice.rs`: `Slice::reference_span` says which bases a slice needs
before any are fetched; `Slice::records_over_window` takes those bases and decodes against
them; `ReferenceSequence` gains a `Window` variant carrying the window's first position.
`src/record.rs` and `src/record/direct.rs` index it exactly as they already index a CRAM's
own embedded reference, which is the same problem — bases whose first byte is not position 1.

A CRAM decodes a mapped read as differences from the reference, and noodles fetched the whole
contig from a `fasta::Repository` to do it. **`Repository::get` also clones the sequence into
its cache, so the peak is two copies of the contig**, and the cache never evicts. A caller
reading the genome in coordinate order does not need any of that: it knows one slice ahead
which bases it will want, because the slice header states the span before a block is decoded.

Measured, decoding the same containers into the same records, `/usr/bin/time -l`:

| | whole contig | per-slice window |
|---|---:|---:|
| tomato chromosome 1 (90.9 Mb) | 198.2 MB | **31.0 MB** |
| human chromosome 1 (249.0 Mb) | 492.8 MB | **76.5 MB** |

Wall time is unchanged on tomato (0.453 s against 0.454 s) and 6 % worse on the human file,
whose slices span 2.5 Mb each so the windows are large and re-read per slice.

**A window that does not cover the slice is refused, not truncated** — the message names both
spans — because the failure it prevents is silent: every read decoded against the wrong bases.

**Not upstream as it stands**, for the same reason as changes 3 and 4. The shape may be worth
proposing, since it is what a coordinate-ordered reader wants and noodles already does the
same coordinate rebasing for embedded references.
