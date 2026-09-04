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
