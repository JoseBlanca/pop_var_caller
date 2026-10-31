# What this copy has that crates.io does not

This is `noodles-bam` 0.89.0 — the version `Cargo.toml` names — extracted from the cargo registry
and then changed. **This file is the whole list of changes.** Anything not named here is
upstream's, byte for byte, and the claim is checkable:

```
diff -r ~/.cargo/registry/src/index.crates.io-*/noodles-bam-0.89.0 vendor/noodles-bam
```

`[patch.crates-io]` in the workspace manifest, and in `fuzz/Cargo.toml`, is what points the build
at this copy rather than at the registry. The registry's `.cargo-ok` and `.cargo_vcs_info.json`
are the extraction's own bookkeeping and are left out.

## The whole diff, file by file

| file | | change |
|---|---|---|
| `FORK.md` | this file | — |
| `src/io/reader/record.rs` | `read_record`, and a test | 1 |
| `src/io/reader/header/reference_sequences.rs` | `read_reference_sequences`, and a test | 2 |
| `src/io/reader/header/reference_sequences/reference_sequence.rs` | `read_name`, and a test | 3 |

## Why a copy at all

**A BAM says how long each of its parts is, and upstream believes it before reading them.** A
record's block size, the header's reference-sequence count and each reference name's length were
each turned into a buffer or a map of that size before a byte of it was read. A damaged or badly
written BAM therefore made pvc ask for up to tens of gigabytes. Without a memory limit that ends
in a clean error, because the operating system lends memory nobody touches; **under one —
`ulimit -v`, or a cluster's `h_vmem` — pvc aborted with no word about which file was at fault.**
Measured on 2026-10-31 with pvc's soak build: a header declaring 2^31 reference sequences aborted
under an 8 GB cap, and a record declaring 4.26 GB under a 4 GB cap. Both were found by fuzzing
(`fuzz/README.md`). Upstream 0.98.0 still reads a record and a reference name this way; it no
longer reserves for the reference-sequence count.

Every change keeps upstream's results on a sound file: the same bytes, the same errors' kinds on
a short one — `UnexpectedEof` — and only the memory taken before the error differs.

## The changes

**1. A record is read into a buffer that grows as the bytes arrive.** `buf.resize(block_size,
0)` then `read_exact` became `buf.clear()` and `read_to_end` through `take(block_size)`, with a
short read refused as `UnexpectedEof`. The buffer is the one every record of a walk reuses, so
after the first few records it already holds the capacity a record needs and nothing is
allocated; what changes is that a damaged size can no longer reserve 4 GiB before the read fails.

**2. The reference-sequence map reserves at most 4,096 entries up front.** It was sized by the
declared count. 4,096 covers GRCh38 with its decoys and alternates (about 3,400); a header with
more grows the map as entries are read.

**3. A reference name is read the way change 1 reads a record**, for the same reason.

## Upgrading noodles

When the workspace moves to a newer noodles, re-extract that version, re-apply these three
changes, and check whether upstream has made any of them unnecessary — 0.98.0 already made
change 2 so.
