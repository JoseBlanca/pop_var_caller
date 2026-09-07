# A4 and A5 — a CRAM block holding reads from several chromosomes

**Date:** 2026-09-07
**Branch:** `ng-cram-window`
**Plan:** [cram_reference_window.md](../../ng/impl_plan/cram_reference_window.md), Milestone A′
**Design:** the owner's ruling at Checkpoint A, 2026-09-07 — support such a block by decoding it
twice, rather than refusing the file

## Verdict

**A CRAM whose block holds reads from several chromosomes now decodes, and holds no chromosome
while doing it.** Before this it aborted the process — `expect("invalid reference sequence name")`
inside noodles — which is what A2 and A3 cost when they stopped handing whole chromosomes to the
decoder. 6,283 library tests pass, up 3; the fork's own 298 still pass.

The two steps went in one commit: A4's fork API has no caller until A5 and A5 cannot compile
without it, and a commit that builds is worth more here than the split.

## What a block like this is, and why it had no answer

A CRAM stores a read as its differences from the reference, so rebuilding one needs the bases
under it. Reads are stored in blocks, and a block's header normally names the one chromosome and
the one stretch its reads cover — which is what A2 made the decode fetch, about 9 kb on a real
file, in place of a whole chromosome.

**A block whose reads sit on several chromosomes names none of them**: no id, no start, no span.
noodles then resolves each mapped read's *whole* chromosome from a repository, and `expect`s a
hit. Against the empty repository a windowed caller passes, that is an abort.

## How the two passes work, and what they cost

The block's header cannot say what to fetch, so the records are asked instead.

1. `Slice::record_extents` decodes the records with **no reference attached at all** and reports
   one extent per chromosome touched — first to last position. This works because a record's
   chromosome, start and CIGAR come out of the file; only its *bases* are rebuilt against a
   reference, and nothing here asks for bases.
2. ng fetches those stretches, one buffer each, reusing the buffers across the container's blocks.
3. `Slice::records_over_windows` decodes again, each record resolving to its own chromosome's
   window.

**The blocks are decompressed once, not twice.** `decode_blocks` already runs before any record is
read, and both passes read its output — so what repeats is the record decode, the cheaper half.
The owner's ruling budgeted for two decompressions; it costs less than that.

**Not taken:** the `.crai` carries each chromosome's stretch for these blocks — the fixture's index
reads `chrom=0 start=300 span=101` and `chrom=1 start=1 span=110` against one container offset — so
the windows could be fetched with no second pass at all. Rejected because it would make the bases a
read is rebuilt from depend on the index being right, where today the index is trusted only for
where to seek.

## Three API changes in the fork, and one of them is a bug fix

`FORK.md` gains change 6.

- **`Slice::reference_span` → `Slice::reference_extent`.** The old `Option` returned `None` for
  both "unmapped, needs no bases" and "several chromosomes, needs one window each", and **that
  fold is what shipped the abort**: A2 read `None` as the first and took a path that is only safe
  for it. Three named states — `OneSequence`, `SeveralSequences`, `Unmapped` — cannot be folded
  that way, and the ng call site is now a `match` with no default arm.
- **`Slice::record_extents`** — the first pass. Returns the extents rather than the records,
  deliberately: nothing was resolved for those records, so asking one for its bases would panic,
  and a function that hands one back invites exactly that.
- **`Slice::records_over_windows`** — the second pass, one window per chromosome.

## The refusal, and why it carries more weight here

**A record whose chromosome has no window, or whose span its window falls short of, is refused
with both spans named.** On the single-chromosome path that check has company: the block header
stores a reference MD5 over its own stretch, and noodles validates the window against it. **A
several-chromosome block has no MD5** — there is no single stretch to compute one over — so the
bound is the only guard between a short window and reads rebuilt against whatever is in the
buffer.

Measured: deleting the bound turns a one-base-short window into an index-out-of-range panic inside
`record/sequence/iter.rs`, not a wrong read — but a panic reached through arbitrary earlier decoding
is not the failure to ship, and a longer window would not panic at all.

## The fixture, and why it is committed rather than built

`src/ng/read/input/testdata/multi_contig_slice.{fa,fa.fai,sam,cram,cram.crai}`, with a README
beside them. The CRAM writer this project vendors puts one chromosome in every block, so no fixture
built the way every other CRAM fixture here is built can reach this code.

**Written by samtools 1.16.1 with no writer options set.** That is the load-bearing part: the
merging is samtools' own default, so the fixture is evidence about files this will actually meet.
24 chromosomes of 400 bases with three 60-base reads each gives **4 blocks, 2 of them spanning
several chromosomes, the largest holding 20**. Two chromosomes alone do *not* trigger it — samtools
needs several under-full blocks in a row (`cram_encode.c`, `c->curr_rec < c->max_rec/4+10`) — which
is why the fixture has 24 and why the realistic case is a draft or scaffold-level assembly rather
than a chromosome boundary.

**The oracle is the SAM samtools was given**, read by the test and compared read for read: name,
position and bases. It is independent of anything ng or noodles does. The reference is pseudo-random
from a fixed seed, because against an all-`A` reference a read rebuilt at the wrong offset comes
back identical and no assertion can fail.

## Tests, and what each fails for

| test | fails when |
|---|---|
| `a_cram_block_spanning_several_contigs_rebuilds_the_reads_its_sam_holds` | any of the 72 reads comes back on the wrong chromosome, at the wrong position, or with bases rebuilt against the wrong window |
| `a_short_window_on_a_block_spanning_several_contigs_is_refused` | a window one base short of what a record needs is decoded against instead of refused |

### Mutation-verified

Each mutation applied alone, reverted from a copy taken beforehand, and the revert confirmed by
`grep` before the next run.

| mutation | result |
|---|---|
| the `SeveralSequences` arm takes the pre-A′ path — no windows, empty repository | the decode aborts: `panicked … invalid reference sequence name`, which is the defect reproduced exactly |
| every window fetched one base late | caught, but by the *single*-chromosome blocks' stored digest first (`reference sequence checksum mismatch`) — so this mutation does not test what it looks like it tests |
| the window fetched one base late **on multi-chromosome blocks only**, which carry no digest | `r02_0's bases were rebuilt against the wrong reference window` |
| the fork's per-record bound check always says the window covers the record | index-out-of-range panic in `record/sequence/iter.rs` |

The second row is worth keeping: a mutation that fails for the wrong reason reads as proof and
is not.

## Validation

Run in the container on this tree.

- `cargo test --lib` — **6,283 pass**, 0 failed, 15 ignored, in 61.13 s.
- `cargo test --lib ng::read` — **308 pass**, re-run as the last thing before `git add`.
- `cargo test` in `vendor/noodles-cram` — **298 pass** (223 + 75), unchanged by the new API.
- `cargo clippy --lib --tests --all-features -- -D warnings` — the three lifetime-elision errors
  in `run/cohort_merge/`, all `main`'s.
- `cargo fmt --all -- --check` — nine files drift, all `main`'s.

## What is not proven here

**No file with such a block has been read at scale.** The fixture is 24 contigs and 72 reads. What
a draft assembly of 50,000 scaffolds costs — how often the second pass runs, and what it adds to a
whole run — is not measured, and Milestone B measures the tomato whole-genome CRAM, which has no
such block. If that number matters, it wants a fragmented-reference CRAM and its own run.

**The extents pass assumes a record's CIGAR needs no reference.** That is how noodles is written —
`read_record` fills the features and `alignment_end()` is computed from them — and the fixture's
reads come back at the right positions, which would fail if it were false. It is not separately
asserted.
