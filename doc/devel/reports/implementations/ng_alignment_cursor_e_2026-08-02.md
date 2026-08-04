# ng — the alignment cursor, Milestone E: the CRAM arm

*Plan: [impl_plan/alignment_cursor.md](../../ng/impl_plan/alignment_cursor.md).
Design: [spec](../../ng/spec/alignment_cursor.md), [arch](../../ng/arch/alignment_cursor.md).
Branch `ng-generic-perf`.*

A cursor can now read a CRAM. Both formats go through one cursor, one filter and one set of
rules; the only thing that differs is finding and unpacking, which is what the design said all
along.

---

## What a CRAM makes different

**A BAM is read record by record; a CRAM is read container by container.** A container holds
around ten thousand records and has to be decompressed and decoded whole before any one of them
can be looked at. So this reader's unit of work is a container: decode one, serve its records in
order, move to the next.

Everything else follows from that.

### The rule that had to be got right

A record reader **positions; it never bounds**. After being pointed at a region, reading on must
yield every record to the end of the chromosome — because the cursor above serves a forward
region by not repositioning at all. A reader that stopped at the previous region's end would
lose every record past it, silently, for every region after the first.

This is where the new arm differs most from the per-region CRAM source it will replace, and it
is not a re-layering: that source stops its index walk at the first container beginning past the
region, and filters each container's records against it. Both belong to a reader that answers
one region and is then thrown away.

### Finding where to start

A binary search on the position each container starts at, then a walk back over earlier
containers whose span reaches into the region.

Both halves are needed. The index is sorted by where a container *starts*, so the search finds
the first one starting at or after the region — but a container beginning earlier can still hold
records reaching into it, and those records are ours to serve.

Walking from the first entry instead would also be correct, and is what the per-region source
does. It cannot be done here: that source is built once per region, while this one is
repositioned on every jump of a walk that may make millions.

### One record, one read group

Every other arm hands a record up with nothing attached, and the layer above resolves the read
group from the record's `RG` tag. **This arm attaches it, because a CRAM has no `RG` tag to
resolve** — the read group is a number the container carries, an index into the header's `@RG`
list. Deciding it while the container is decoded is what lets every auxiliary tag be dropped,
and re-inflating a number into a string so the layer above could parse it back would be
perverse.

So the answer travels with the record, and the layer above uses it when it is there and resolves
when it is not. That is the one exception in a contract that otherwise reads as absolute, and it
is written down where the contract is.

---

## The oracle

**A run of regions through one CRAM cursor, against a linear scan of the whole file** — the same
shape as the BAM one, sharing nothing with what it checks: the scan never seeks and never opens
the index.

It runs on a fixture of 25,000 reads over 400,000 bases, which is **three containers**. A
fixture inside one container exercises the decode and none of the walk — the CRAM shape of the
trap the BAM tests already carry, where every region resolved to one chunk and an oracle passed
with the defect it existed to catch.

**Its regions include a probe just inside each container's first base, taken from the index
itself** rather than hard-coded, so it stays a boundary whatever the writer's container size
turns out to be. That detail is load-bearing: without it every region sits comfortably inside a
container, and the walk-back can be deleted with the whole suite green. It was checked, and it
could. With them, two mutations die — bounding the walk at the region's end, and dropping the
walk-back.

### ⚠ One guard is untested, and the fixture is why

A container may hold several slices, each its own index entry sharing the container's offset.
Decoding it once per entry would serve its records twice. noodles' writer puts one slice in each
container, so no fixture this project can build reaches that branch — deleting it leaves the
suite green. samtools writes multi-slice containers, so the guard is needed for real input and
cannot be exercised by input we can produce. Testing it needs a fixture builder that can write
several slices per container, or a committed CRAM from another writer.

---

## The first CRAM measurement

CRAM was unmeasured in the performance review. These are its first numbers, on a tomato sample
from the `ssr_tomato1` cohort — a 37 MB targeted slice, aligned to SL4.0.

### The read path

`ng_cursor_vs_query`: the same typed regions, the same reads, once through each path, with every
read compared over its whole content in an untimed second pass.

| file | chromosome | regions | before | after | |
|---|---|---:|---:|---:|---:|
| SRR5079860 | ch01 | 263,800 | 12.78 s | 0.54 s | **23.7×** |
| SRR5079860 | ch01 (repeat) | 263,800 | 12.81 s | 0.54 s | 23.7× |
| SRR5080000 | ch01 | 263,800 | 14.32 s | 0.59 s | 24.5× |
| SRR5079860 | ch02 | 160,298 | 9.12 s | 0.35 s | 26.4× |

Every run reports `agreement=exact` — the reads are identical, compared element by element.

The cursor decodes **80,510** reads to serve **1,444,307**, reusing on 263,799 of 263,800
regions and jumping once.

### End to end

The whole generic walk, which is what a user actually runs:

| | seconds | peak RSS | loci |
|---|---:|---:|---:|
| a query per region | 9.51 | 238 MB | 1,711,775 |
| one cursor per chromosome | 6.76 | 228 MB | 1,711,775 |

**1.41×, with identical output** — 1,711,775 loci, 1,718,908 observations, 105,894 reads
admitted, on both.

### Why 1.41× end to end when the read path is 23.7×

Two reasons, and the first is the interesting one.

**The old CRAM path already kept its last container.** Spec §1 says so: *"the CRAM path already
retains across queries through the pool"*. So on CRAM the retention this whole design is about
was partly there — the per-region query re-consulted the index and re-narrowed, but often did
not re-decode. There was less to win than on BAM, where nothing was kept and the same records
were decoded thirty times over.

**And the walk is a larger share here.** This chromosome yields 1.7 M loci against chromosome
21's 236 k, so the per-region read cost is diluted by far more walking. The BAM figure at
Milestone D was 2.41× on the same shape of change.

Both numbers are real and they measure different things. Quote the read-path figure only with
the region shape it was taken on.

### ⚠ The memory number is about the reference, not the reads

228 MB looks alarming beside the 21 MB a BAM walk of chromosome 21 peaks at, and almost all of
it is one thing: **a CRAM decodes against the reference, so the chromosome's bases are
resident.** SL4.0 chromosome 1 is 90.9 Mb, and the repository holds it for as long as a cursor
on that chromosome lives.

That is the cost spec §10 describes and the reason a cursor covers one chromosome: on CRAM a
chromosome change means re-reading hundreds of megabytes. It is also what the deferred
per-chromosome reference registry (spec §12) is for — its trigger was named as *the first
parallel run over CRAM*, and this is the number that would multiply.

The cursor does not make it worse: 238 → 228 MB.

---

## Verification

`cargo test --lib ng::` **1,572 passed**; clippy with warnings as errors, `cargo fmt --check`
and `cargo test --examples` clean.

**Nothing moved on BAM.** Both dumps stay byte-identical to binaries built from `ee0c94b`, and
`ng_generic_walk_probe` on chromosome 21 prints `loci=236081 observations=251786
reads_admitted=54709`.
