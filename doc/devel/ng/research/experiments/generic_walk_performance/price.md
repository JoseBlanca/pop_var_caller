# What one contributor visit costs the ng generic walk, re-taken on the shipped binary

Measured on `1e5ffa8`, host-native builds, `instructions retired` from `/usr/bin/time -l`,
floor-subtracted, min of 3 runs. Wall clock is never quoted: on this host two runs of an
*identical* binary have reported six-point swings.

Vocabulary, defined once and used unchanged from the 2026-08-05 census so the two sets of
numbers are comparable:

- a **column** is one covered reference base the walk turns over — one call of
  `WalkerState::process_position`, including the 33 columns in 10,000 where every active
  read is silent (a deletion interior, an `N`-skip) and no visit happens;
- a **contributor visit** is one (column × read) pair that survives the mate-overlap
  collapse and the per-position depth cap — the unit the fold iterates over;
- **depth** means contributor visits per column, which is what the walk pays for, not the
  sample's nominal coverage;
- the **ordinary-column path** is `fast_column::try_ordinary_column`, the scalar path that
  landed today; the **general path** is everything else.

Every figure below was measured in this session unless it says *cited*, in which case it
comes from the brief or from the 2026-08-05 census beside it.

---

## The answer

**The factor of three to six is gone.** On a BAM at the 30× target a contributor visit now
costs **659 instructions** against the 300–600 that a hand-count of the work it performs
gives (cited) — a factor of 1.1 to 2.2, where before today it was 1,312 against the same
300–600, a factor of 2.2 to 4.4. The walk is close to the cost of its own arithmetic on the
per-visit axis, and no structural change reaches a large win there. What is left is smaller
and sits on a different axis: the **depth-independent per-column term, 4,533 instructions**,
which three rounds of work have moved by 14 % while the per-visit term fell by half, and
which is **a quarter of what a 30× covered base costs** (4,533 of 18,116) against a
twenty-fifth of what a 130× one costs (4,533 of 112,360). It is a breadth term, and 30× WGS
is the workload where breadth dominates.

---

## Q2 — the price, before and beside

| | measured at `6fbbd09` (cited) | measured now at `1e5ffa8` | change |
|---|---:|---:|---:|
| per contributor visit, BAM | 1,312 | **659** | −50 % |
| per contributor visit, CRAM | 1,781 | **1,103** | −38 % |
| per column, depth-independent | 5,260 | **4,533** | −14 % |
| hand-count of the work a visit performs (cited) | 300–600 | 300–600 | — |

**Which term moved is the informative part, and the answer is the per-visit one.** A fast
path that answers 7,378 columns in 10,000 at ~130× and 6,987 in 10,000 at 30× cut the
per-visit price by half on a BAM and by 38 % on a CRAM, and cut the per-column fixed cost by
14 %. That is the expected shape and it is worth stating plainly: the ordinary-column path
replaced the *per-read* work — a cursor re-query, a hash-map insert keyed by read id, an
element in each of two depth-sized sorts, a re-derived allele string — with six scalars per
read. It did **not** replace what a column costs regardless of its reads: fetching the
reference base, building the emitted locus, and allocating the boxed slices the locus
carries. The fast lane still does all of that; it just does it in its own module.

### The marginals it was fitted from

Differencing two runs of one fixture at two loci counts cancels start-up and reference
loading. Instructions are the min of 3; columns and visits come from a default-off
`walk-reprice` build whose instruction counts are never quoted.

| marginal | Δ instructions | Δ visits | depth | instr/visit | instr/column | instr/visit at `6fbbd09` (cited) |
|---|---:|---:|---:|---:|---:|---:|
| tomato ~130× CRAM, region stream, 1 M → 2 M loci | 129.257 G | 118,793,400 | 82.5 | **1,088** | 89,758 | 1,772 |
| tomato ~130× CRAM, whole contig, 1 M → 2 M | 112.807 G | 98,178,421 | 97.8 | **1,149** | 112,360 | 1,835 |
| HG002 30× BAM, whole contig, 770 k → 1.54 M | 14.023 G | 15,966,198 | 20.6 | **878** | 18,116 | 1,567 |
| HG002 30× BAM, region stream, 770 k → 1.54 M | 42.594 G | 21,562,917 | 17.3 | 1,975 | 34,124 | 2,686 |
| tomato sparse CRAM, whole contig, 400 k → 800 k | 3.967 G | 1,949,612 | 4.9 | 2,035 | 9,896 | — |

Two of those Δ-visit counts — 118,793,400 and 98,178,421 — are **identical to the digit**
with the ones the 2026-08-05 census measured on the same marginals. The walk folds exactly
the same visits; only what each costs has changed. That is as close to a controlled
before-and-after as this instrument offers.

The last row is the contaminated one and the row above it is the fit's low-depth point.

### The split

Two depths on the same reference, the same contig and the same file format — the sparse
tomato bench CRAM and the big one over `SL4.0ch01`, both in whole-contig mode:

```
sparse:   9,896 instructions per column at depth  4.86
big:    112,360 instructions per column at depth 97.79
        → per column 4,533 + per visit 1,103 × depth        (was 5,260 + 1,781 × depth)
```

Feeding the HG002 BAM's whole-contig marginal — 18,116 instructions per column at depth
20.63 — through the same 4,533-instruction column term gives **659 instructions per visit on
a BAM**, against 1,103 on a CRAM. The 444-instruction gap is the input format: a CRAM pays
rANS entropy decoding, a BAM pays BGZF inflation, and Q1 below prices both.

So at the 30× target one covered reference base costs

```
4,533 + 20.6 × 659 = 18,116 instructions        (before: 5,260 + 20.6 × 1,312 = 32,287)
```

— **44 % less than this morning, on the same seventeen-to-twenty-one reads.**

### The caveat that has to travel with every 30× number

**The human fixture is tandem-repeat-targeted, not whole-genome.** Walking `chr1` through the
shipped region stream costs 34,124 instructions per column at depth 17.3; walking the same
BAM in whole-contig mode costs 18,116 per column at depth 20.6. Correcting the second for the
depth difference leaves **18,215 instructions per column that are region setup rather than
walking — 53 % of that fixture's instructions**, against 18,480 and 40 % at `6fbbd09`
(cited). The absolute overhead is unchanged to within 1.4 %; its *share* rose because the
walk beside it got cheaper. Every 30× figure above says which mode it came from, and the
region-stream row is quoted only to size this tax.

---

## Q1 — where the time goes now

`sample` at 1 ms against a `[profile.profiling]` build. **That profile sets `lto = false` and
`codegen-units = 16`, so its self-times say *where* work is and their sizes do not transfer to
the release binary the instruction counts above were taken on.** This caveat applies to every
percentage in this section.

### tomato ~130× CRAM, `SL4.0ch01`, shipped region-stream mode, 3 M loci, 14,700 self samples

Self time, top of stack:

| site | share | what it is |
|---|---:|---|
| `noodles_cram rans_4x8 order_1::decode` | 11.46 % | CRAM entropy decoding |
| `<deduplicated_symbol>` | 8.27 % | three things — resolved below |
| `open_record::process_position` | 6.38 % | the general fold |
| `_xzm_free` | 5.50 % | the allocator |
| `fast_column::try_ordinary_column` | 5.22 % | the ordinary-column path |
| noodles CRAM `sequence::iter::Iter::next` | 3.88 % | CRAM record building |
| `PileupWalker::next` | 3.52 % | the walk loop itself |
| `_xzm_xzone_malloc_tiny` | 3.29 % | the allocator |
| `_platform_memmove` | 3.18 % | — |
| `CigarCursor::match_at` | 3.01 % | the fast lane's per-read query |
| `OpenPileupRecord::finalise_recycling` | 2.82 % | closing a general-path record |
| `RecordBuf::try_from_alignment_record` | 2.80 % | noodles record conversion |
| `CigarCursor::events_at` | 2.63 % | the general path's per-read query |
| `open_record::apply_events_into` | 2.49 % | the general fold |
| `ng::tandem_repeat::find_tandem_repeats` | 0.29 % | region typing |

Inclusive, from the call graph:

| subtree | share |
|---|---:|
| getting reads out of the file (`PreparedSampleReads::next`) | **52.5 %** |
|  └ CRAM record decode + conversion (`CramAlignedReadsReader::read_next`) | 45.4 % |
| the general fold (`process_position`) | 12.1 % |
| the ordinary-column path (`try_ordinary_column`) | 10.3 % |
| chain-id allocation (`allocate_for_read`) | 5.8 % |

**On a CRAM at ~130×, more of the run is now spent reading the file than walking it.** That
is a consequence of the walk getting 38 % cheaper around a decoder that did not change.

### The `<deduplicated_symbol>` line, resolved

The linker folds identically-compiled functions onto one address, so `sample` prints one
ambiguous name. Resolving each frame's address against `nm` (load slide `0x4124000`,
script `resolve_dedup.py` beside this file) splits the 8.27 % three ways:

| what it actually is | share of the line | share of the run | reached from |
|---|---:|---:|---|
| `hashbrown::HashMap::retain` | 57.0 % | 4.71 % | `ChainIdAllocator::evict_stale_pending` |
| `core::hash::BuildHasher::hash_one` | 21.1 % | 1.74 % | noodles CRAM `Byte::decode_take` |
| `core::slice::sort::smallsort::small_sort_general` | 11.1 % | 0.92 % | the two remaining sorts |
| `small_sort_network`, `HashMap::remove`, `RawTable::clear`, `HashMap::insert` | 10.8 % | 0.89 % | — |

`pending_mates.retain` is the only `retain` over a map anywhere in ng's walk
(`chain_id_allocator.rs:431`), so that first row is `evict_stale_pending` unambiguously.
**Half the line is one function, and it is lever 1.**

### HG002 30× BAM, `chr1`, whole-contig mode, 1,055 self samples

Small sample; a second run over every contig of the same BAM (16,206 self samples) agrees
within two points on each row and is the corroboration, not the headline — it carries an
artefact of its own, 12.6 % spent re-reading the BAM's 3,366-contig header once per contig,
which no real run pays.

| subtree, inclusive | chr1 | all contigs, header cost removed |
|---|---:|---:|
| getting reads out of the file | 26.3 % | 24.8 % |
|  └ BGZF inflate | 10.3 % | 12.1 % |
| the ordinary-column path | 21.3 % | 20.3 % |
| the general fold | 15.7 % | 19.9 % |
| fetching reference bases (`RawChromReader::fetch`) | 12.8 % | 11.6 % |
| `evict_stale_pending` (from the resolved dedup line) | — | 1.6 % |

**The shipped region-stream run on this fixture is not this profile.** Run through the region
stream, `chr1` spends **64.8 % of its self time in `find_tandem_repeats`** — the scan that
types regions. That is not a walk cost and it is not a finding about the walk: the BAM covers
only the GIAB tandem-repeat tiers, so typing runs over all 248 Mbp of `chr1` to deliver
1.5 Mbp of covered bases, a ratio of 165 to 1. On the actual target — one human sample at 30×
covering the whole genome — that ratio is 1 to 1, and the tomato fixture is the one that
shows it: `find_tandem_repeats` is **0.29 % of the run** there, against a cited 0.17 %.

---

## Q3 — the census, re-taken against a walk with two paths

The census (`PVC_COLUMN_CENSUS=1`) is incremented **after** the ordinary-column path has had
its turn and returned, so since that path landed its counters describe the columns the fast
lane handed back, not the walk. The visit-side denominator was missing entirely; the
default-off `walk-reprice` feature supplies it.

### How much the ordinary-column path takes

| | tomato ~130× CRAM | HG002 30× BAM |
|---|---:|---:|
| columns walked | 1,449,893 | 2,499,261 |
| answered by the ordinary-column path | **7,378 in 10,000** | **6,987 in 10,000** |
| contributor visits | 123,613,443 | 43,084,509 |
| visits absorbed by the ordinary-column path | **7,272 in 10,000** | **6,390 in 10,000** |

Both column counts are identical to the 2026-08-05 census's, and both visit counts are too.

### How the columns that still reach the general path differ

**They are deeper, and at 30× substantially so.**

| mean depth of… | tomato ~130× | HG002 30× |
|---|---:|---:|
| every column | 85.3 | 17.2 |
| a column the ordinary path answered | 84.0 | 15.8 |
| a column that reached the general path | **89.8** | **22.7** |

At 30× the general path holds 27 columns in 100 but **36 visits in 100**, because its columns
carry 44 % more reads than the fast lane's. At ~130× it holds 26 columns in 100 and 27 visits
in 100, so there the gap is 7 % and barely tilts the visit share. The reason is the
fast lane's second condition: it requires *every* read the active set holds to have a CIGAR
free of `I` and `D`, so the more reads a column holds the likelier one of them disqualifies
it. That is a per-column all-or-nothing test over a per-read property, and depth is what
makes it bind.

**They do not carry wider records.** Measured at the close, over the general path's own
records only:

| | tomato ~130× | HG002 30× |
|---|---:|---:|
| records the general path closed | 370,980 | 675,946 |
| reference bases they span, summed | 373,859 | 681,783 |
| **mean footprint** | **1.008 bases** | **1.009 bases** |

So the general path is running machinery built for a record that widens, on a record that is
one reference base wide 992 times in 1,000. The excess is 2,879 bases on tomato and 5,837 on
HG002 — if every widened record is exactly two bases, that is 78 records in 10,000 and 86 in
10,000 respectively.

**Why each column was handed back.** Shares of the general path's columns; a column can fail
more than one test, so they sum past 100 %.

| the column was handed back because… | tomato ~130× | HG002 30× | of *all* columns |
|---|---:|---:|---|
| some active read's CIGAR carries an insertion or deletion | 65.9 % | 70.0 % | 1,706 / 1,915 in 10,000 |
| some active read carries a deletion specifically | 44.6 % | 45.5 % | 1,155 / 1,244 in 10,000 |
| two contributors share a chain id (mates overlap) | 40.5 % | 38.4 % | 1,048 / 1,052 in 10,000 |
| a contributor's own event here is an indel | 0.9 % | 0.6 % | 24 / 16 in 10,000 |
| a record is already open over this base | 0.7 % | 0.8 % | 19 / 21 in 10,000 |
| the per-position depth cap bound | 0 % | 0 % | 0 |
| the column mixes read groups | 0 % | 0 % | 0 |

Two reads of this. First, **the indel test is what the general path is a minority of**: two
thirds of its columns are there because some read *held at that column* carries an indel
somewhere in its CIGAR, not because anything indel-shaped happens at that base — a
contributor's own event is an indel at only 24 columns in 10,000. Second, the depth cap and
the read-group test never fired on either fixture, so two of the seven conditions are
currently free.

---

## The three levers, re-measured

### Lever 1 — the chain-id allocator's sweep over its pending-mate table

`ChainIdAllocator::evict_stale_pending` runs a `retain` over the whole map of first mates
awaiting a partner, once per admitted read. **It did not get worse on ordinary data, and it
did get worse in the one place the raised ceiling lets the map grow.**

Priced in release instructions by running the sweep a second time over the same map and
differencing (min of 3). The second pass touches each entry through `std::hint::black_box`:
without that, a `retain` whose closure is a constant `true` has no side effect and LLVM
deletes the loop — the first attempt measured the sweep at 0.04 % of a run the profiler puts
at 4.7 %, which is the shape of a deleted loop, not a cheap one.

| | HG002 30× `chr1` | tomato ~130×, first 1 M loci | tomato ~130×, the deep region |
|---|---:|---:|---:|
| mean entries the sweep scans, per admitted read | **43** | **160** | **306** |
| the same, through the shipped region stream | 21 | — | — |
| peak entries held (`pending_mates_high_water`) | 132 on chr21 | 363 | **11,384** |
| instructions the sweep costs | 0.215 G | 4.674 G | 23.559 G |
| **share of the walk** | **0.71 %** | **3.72 %** | **9.90 %** |
| for comparison, `6fbbd09` (cited) | — | 3.6 % | — |

The deep region is `SL4.0ch01` from 32.5 Mbp, 1.5 M loci, reached with `PVC_PROBE_FROM_BP`; its
share is against the walk alone, with the 7.534 G cost of reading the reads in front of it
subtracted. Its peak of 11,384 pending mates is the figure the brief cites, reproduced.

The per-entry cost is not constant: 19.5 instructions where the map averages 160 entries,
41.6 where it averages 306 and peaks at 11,384. A map of 11,384 entries is about 546 kB and
does not fit this host's 128 kB L1 data cache, which is the likeliest reason. The profiler
agrees with the doubling to within a factor of 1.3: `HashMap::retain` is 4.71 % of the
3 M-locus tomato run against 3.72 % measured on the 1 M prefix.

### Lever 2 — the columns the ordinary-column path hands back for mate overlap

Two contributors from one fragment need their qualities reconciled against each other, which
the scalar path has no term for, so it hands the column back.

| | tomato ~130× | HG002 30× |
|---|---:|---:|
| columns handed back **only** for mate overlap | **864 in 10,000** | **807 in 10,000** |
| column coverage today | 73.8 % | 69.9 % |
| column coverage if those were handled | 82.4 % | 77.9 % |
| **visit coverage today** | **72.7 %** | **63.9 %** |
| **visit coverage if those were handled** | **83.3 %** | **77.0 %** |
| gain, in points of visit coverage | **+10.5** | **+13.1** |
| for comparison, sized at 300× (cited) | +28 | — |

So the lever is worth 10 to 13 points of visit coverage at the two depths that matter, against
the 28 points it was sized at on 300× data — it shrinks as depth falls, and the target is
30×. Reported, not proposed; it has been declined twice on complexity.

### Lever 3 — input decoding is format, not code

Inclusive share of the profiling build's main thread, so the reader can separate walk cost
from decode cost. Same `lto = false` caveat.

| | tomato ~130× CRAM | HG002 30× BAM |
|---|---:|---:|
| getting reads out of the file, all of it | **52.5 %** | **26.3 %** |
| the record decoder itself | 45.4 % (CRAM) | 13.5 % (BAM) |
| entropy/compression codec inside it | 12.6 % (rANS 4×8) | 10.3 % (BGZF inflate) |
| cited at `6fbbd09` | ~10 % of the walk | — |

The cited ~10 % and the 45.4 % measured here are not the same measurement — the 2026-08-05
census put `decode_container_at` at 26.0 % of its own profiling run — but the direction is
unambiguous and was predicted: **a decoder that did not change is a larger share of a walk
that got 38 % cheaper.** On a CRAM the walk now has less than half the run left to act on.

---

## What is no longer true, and what still is

**Re-taken: the working set.** The active set can now hold 32,768 reads instead of 4,096, so
the 2026-08-05 finding that a column fits in this host's 128 kB L1 data cache needed
re-checking. `ActiveRead` is still **184 bytes** — none of the four reverted narrowing
attempts changed it — so the ~590 bytes per read the earlier table totals still stands, and
the L1 boundary is 222 reads.

| | HG002 30× `chr1` | tomato ~130×, 1 M loci | tomato ~130×, 2 M loci | tomato ~130×, deep region |
|---|---:|---:|---:|---:|
| mean reads held per column | 17.3 | 86.0 | 84.7 | 86.9 |
| working set at the mean | 10.0 kB | 49.6 kB | 48.8 kB | 50.1 kB |
| **columns whose working set left L1** | **0 of 2,499,261** | **0 of 1,449,893** | 72 of 2,889,945 | **2,271 of 1,860,533** |
| as a frequency | 0 in 10,000 | 0 in 10,000 | 0.2 in 10,000 | **12 in 10,000** |

**The conclusion holds at both target depths and stops holding in the deep region.** At 30×
and over the first million loci at ~130×, not one column's working set left L1. In the deep
region 12 columns in 10,000 do, and the deepest column there holds 10,744 reads — a working
set of about 6.0 MB, forty-eight times L1 and inside the 16 MB shared L2. So cache behaviour still
does not explain the per-visit price anywhere the walk normally is, and the struct-narrowing
experiments stay refuted; the deep region is the one place where that argument would have to
be re-made, and it is 12 columns in 10,000 of one contig.

**Not re-run, per the brief:** the four struct-narrowing attempts, the per-read CIGAR hint,
the sorted-`Vec` fold table against a permuted arrival order, the three fold-table size
estimators, `mimalloc`, noodles' lazy record type, caller-side region coalescing, run-length
compression of the repeat scanner, and candidate-site-only emission.

---

## Q4 — the verdict

**There is no longer headroom of the kind a structural change reaches on the per-visit axis,
and the remaining opportunity is a quarter of a 30× covered base, sitting in the
depth-independent per-column term.**

The per-visit case first, because it is the one that closes. A contributor visit costs 659
instructions on a BAM and 1,103 on a CRAM, against 300–600 for the work it performs (cited):
pull at most one match from a cursor, compare one base against one allele bucket, add six
scalars, and — on the general path only — insert one hash-map entry. On a BAM that is a
factor of 1.1 to 2.2. The 2026-08-05 census called the same gap "a factor of three to six"
and argued a specialised path would close it; the path was built, and it did. The remaining
444 instructions between BAM and CRAM are the decoder, not the walk.

What did not move is the fixed cost of turning over one column: 5,260 instructions this
morning, **4,533 now**, a 14 % cut against 50 % on the visit. At ~130× that term is a
twenty-fifth of a column and irrelevant. **At 30× it is a quarter of one** — 4,533 of the
18,116 instructions a covered base costs — and 30× whole-genome is the stated target. Its
size in the shipped walk is 4,533 × the number of covered bases: on a 3 Gbp human genome at
30×, about **1.4 × 10¹³ instructions of column overhead** before a single read is looked at.

Where it sits, named without designing anything: every covered reference base produces one
`SampleLocusObservations`, and that object owns a boxed slice of reference bases, a `Vec` of
observations, and per observation a boxed slice of allele bytes and a `Vec` of chain ids. At
the cited 340 instructions per allocation, 4,533 instructions is the price of roughly
thirteen of them. The ordinary-column path allocates the same object as the general path;
that is why specialising the column moved the per-visit term and left this one alone. The
profile is consistent: the allocator's `free` and `malloc_tiny` together are 8.8 % of the
tomato run's self time, fourth and eighth in the ranking.

Two smaller things are worth their one line each and no more. The general path answers 27
columns in 100 at 30× and 26 in 100 at ~130×, and the records it closes are **1.008 reference
bases wide** — general machinery on a one-base record, and two thirds of those columns are
there because a read held at that column carries an indel somewhere else in its CIGAR. And
`evict_stale_pending` is 0.71 % of the walk at 30× but **9.90 % in the deep region** where
the raised hold ceiling lets its map reach 11,384 entries.

On a CRAM, none of this reaches more than half the run: getting reads out of the file is
52.5 % of the profiled time, and that is the decoder's cost, not the walk's.

---

## Method, and the gates

**Instrument.** `instructions retired` from `/usr/bin/time -l`, floor-subtracted, min of 3.
Three consecutive runs of one binary agreed to better than one part in a thousand throughout
and usually to one in ten thousand. Wall-clock A/Bs and criterion's `change:` line are not
admissible on this host.

**Floors, measured per binary per fixture with `PVC_PROBE_MAX_LOCI=1`:**

| fixture | floor | cited |
|---|---:|---|
| tomato big CRAM, `SL4.0ch01` | 1.3089 G | 1.306 G ✓ |
| HG002 30× BAM, `chr1` | **0.3500 G** | 0.349 G ✓ — **not** the 1.900 G that circulated |
| tomato sparse bench CRAM, `SL4.0ch01` | 1.1549 G | — |
| tomato big CRAM from 32.5 Mbp (reads in front are still read) | 7.5339 G | — |

Every marginal above differences two runs of one fixture, so the floor cancels and is quoted
only for reference. `PVC_TRUST_REFERENCE_INDEX=1` throughout — the mode is echoed as
`reference_check=trusted_unverified` — so no run here is compared against one that verified
the FASTA.

**Instrumentation, all default-off, and the release binary provably unchanged.** Two cargo
features were added, `walk-reprice` (the counters the fast lane put out of the census's reach)
and `walk-evict-doubling` (the sweep run a second time, to price it). The default release
build's `__TEXT,__text` section is **byte-identical** to the one built from a clean tree
before any of this was written, for all three examples the gates use — `ng_generic_walk_probe`,
`ng_generic_loci_dump`, `ng_ssr_loci_dump`. That is the standard the earlier census set and it
is met.

**Gates, all on the default build:**

- `ng_generic_loci_dump` chr21 — 251,792 lines, `diff` against the stored baseline shows
  **exactly one line**, the `record_widen_events=423 → 425` header, which the owner accepted;
- `ng_ssr_loci_dump` chr21 — 4,406 lines, byte-identical;
- `ng_generic_loci_dump` `SL4.0ch01` — 1,718,914 lines, one line, `record_widen_events=622 → 621`;
- `ng_ssr_loci_dump` `SL4.0ch01` — 11,945 lines, byte-identical;
- probe counters on chr21 exact: `loci=236081 observations=251786 reads_admitted=54709
  fast_columns=262498 mate_overlap_positions=39312`;
- `cargo test --lib` — 2,893 passed, 1 failed, and the failure is
  `parity::every_divergence_from_production_is_one_of_the_six_named_classes`, the accepted
  clean-tree divergence;
- `cargo test --examples` — 33 targets, all ok;
- `cargo clippy --all-targets --all-features -- -D warnings` — clean, which covers both new
  features since `--all-features` enables them;
- `cargo doc --no-deps` — 12 unresolved links, the baseline, all in `ssr` and `em` modules.

**One incidental finding, not acted on:** `cargo fmt --check` is **not** clean on `1e5ffa8`.
Seven hunks in `open_record.rs` and one in `cigar_cursor.rs` differ, in code this session
never touched. The two hunks this session's own edits introduced were fixed; the eight
pre-existing ones were left, because reformatting files this review only reads would put
noise in a diff whose whole point is that the shipped binary did not change.

**Raw output beside this file.**

- `instructions_raw.txt` — every `instructions retired` reading, three per configuration, in
  run order, including the failed first attempt at lever 1 and why it failed.
- `sample_tomato_130x.txt`, `sample_hg002_30x_region.txt`, `sample_hg002_30x_whole.txt`,
  `sample_hg002_30x_allcontigs.txt` — the four profiles; `nm_prof.txt` is the symbol table
  the deduplicated-symbol split resolves against.
- `counts/` — the probe's own output and census counters at every endpoint;
  `columns_walked.txt` collects the column counts the price table divides by.
- `gates/` — the four dump diffs (`diff_ssr_*.txt` are empty, `diff_generic_*.txt` are the
  one accepted header line each) and the chr21 probe counters.
- `text_section_sha256.txt` — the `__TEXT,__text` hashes before and after the
  instrumentation, for all three example binaries. All three pairs match.
- `walk_reprice_instrumentation.patch` — **the instrumentation itself, 176 added lines and
  no line removed.** It was written in a throwaway worktree, so this patch is the only copy;
  `git apply` it on `1e5ffa8` to reproduce anything here.
- `price.py` (the price table's arithmetic), `resolve_dedup.py` (the deduplicated-symbol
  split), `rank_sample.py` and `tree_sample.py` (self and inclusive profile rankings),
  `prof_*.sh` and `columns.sh` (the run drivers).
