# The 268 true tract sequences no read reported

`ng_tract_candidate_recall.py` scores 20,006 of GIAB's 20,204 Tier tracts at 30x and finds 434
places where a sequence HG002 really carries was not among ng's candidate alleles. 268 of those
are the bucket nobody had opened: ng's cohort-merge allele table never held the sequence at all,
so no support bar and no ploidy cut could have dropped it. This is what those 268 are.

Everything below is measured on `tmp/attrib/tier_30x_candidates.tsv` (the 30x candidate dump),
`bam/30x` and `bam/300x` of the GIAB Tier BAMs, and the phased truth VCF
`HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz`. It is a fact about this one high-coverage human
sample against this one benchmark; it says nothing directly about a 3x cohort.

## The answer in one line

**Half of the 268 are not ng's to fix, a quarter are, and the quarter that is is worth about 46
tracts in 20,006 — around 2 wrong genotypes in every 1,000 tracts.** The single most useful number
is a different one: **at 195 of the 268 ng had already offered the correct repeat LENGTH**, and
**174 of the 268 differ from the reference by no length at all**. What is missing at three quarters
of them is a spelling, not a length.

## The split

| bucket | count | share of 268 |
|---|---|---|
| (a) unobservable — the allele cannot be spanned by a 150-base read, or nothing spans the tract | 14 | 1 in 19 |
| (b) alignment loss — reads carry the sequence, ng's table does not hold it | 67 | 1 in 4 |
| (c) not in the reads at all, at 300x — truth-side | 121 | 45 in 100 |
| ambiguous — every base the truth needs is in the reads, but no read spells the tract that way | 66 | 1 in 4 |

Of the 67 in (b), **46 are unambiguous at 30x**: two or more spanning reads at 30x reproduce the
missing tract sequence exactly. 8 more have exactly one such read. The remaining 13 have none at
30x and several at 300x — those are the coverage coin-flips, and they are the *only* part of the
268 that low depth alone explains.

Of the 46, ng had already offered the right tract length at 32. **At 14 it offered neither the
sequence nor its length** — those 14 are outright wrong genotypes, recoverable at 30x, and they
are listed in full below.

### How each bucket was decided

Two tests, both proof against where an aligner chose to place an indel.

* Where the missing sequence is the reference length (173 of 268, so the truth records inside the
  tract are substitutions): find every reference position where the missing sequence differs from
  the reference and count the reads' aligned bases there. A base inside an `M` block cannot be
  moved by indel placement, so this is exact. "Supported" means every differing position is
  carried by at least 3 reads.
* Where the length differs (94 of 268): a read's implied tract length is (read bases aligned
  inside `[start-12, end+12]`) − 24. Wherever inside that padded window the aligner put the indel,
  the count is the same. "Supported" means at least 3 reads at the missing length.

A case lands in (b) only if it passes the test *and* at least 2 reads reproduce the whole tract
string exactly. Cases that pass the positional test but whose exact string no read spells are the
66 "ambiguous" — see below.

## (a) Too long to span: 13 tracts, plus 1 with no spanning read

Reads are 150 bases; naming a tract's length needs the tract plus flank on both sides. With 20
bases of flank the ceiling is 110 bases. 13 of the 268 missing alleles are longer than that:

    chr6:35794661-35794694    12,256 bp     0 spanning reads even at 300x
    chr7:105908612-105908623     564 bp   171 spanning reads at 300x, 0 carry it
    chr4:926876-926883           267 bp    47
    chr4:147095874-147095889     180 bp   235
    chr14:103145125-103145140    159 bp   216
    chr3:20018677-20018702       158 bp   191
    chr11:1543024-1543051        154 bp   120
    chr19:55931710-55931725      151 bp    65
    chr12:577417-577442          116 bp    86
    chr4:135007835-135007846     116 bp   167
    chr4:183625314-183625341     112 bp   139
    chr12:7957526-7957585        112 bp    73
    chr1:8276856-8276879         757 bp   105 spanning reads at 300x, 0 carry it

At every one of these, and at 300x where the median tract has 115 spanning reads, **not one read
carries the sequence**. That is exactly what "too long to span" predicts, and it is why this
bucket needs no further work. (chr1:8276856 appears twice in the 268 because both of HG002's
haplotypes were missing there: the 757-base one, listed here, and a 24-base one that 20 reads at
300x do carry — that second one is counted in bucket (b), not here.)

Only 4 of the 13 are homopolymers or dinucleotides shorter than 30 reference bases — these are
tracts where HG002 carries a large expansion relative to a short reference tract, which is the
hardest case for any short-read caller and not a defect in ng.

## Low coverage: a contributing factor, not a deciding one

Depth at the tracts that lost a sequence, measured as the reads in ng's own merge table:

| | median | p25 | p75 | under 10 reads |
|---|---|---|---|---|
| the 268 | 16 | 11 | 22 | 53 of 268 (1 in 5) |
| the 8,280 built tracts where both true sequences were offered | 21 | 16 | 27 | 537 of 8,280 (1 in 15) |

So the 268 sit about 5 reads shallower and are three times as likely to be under 10 reads. But
that difference does not carry the bucket: of the 121 in (c), the *300x* BAM gives a median 115
spanning reads and a minimum of 16, and the sequence is still absent. Coverage explains 13 of the
268 cleanly — the ones invisible at 30x and present at 300x.

## The shape of the 268 against the tracts that lost nothing

Comparison set: the 8,280 tracts the merge built where both true sequences were offered.

**Period.** Dinucleotides are over-represented and homopolymers under-represented:

| period | the 268 | comparison |
|---|---|---|
| 1 | 107 (40 in 100) | 5,301 (64 in 100) |
| 2 | 124 (46 in 100) | 2,060 (25 in 100) |
| 3 | 11 | 273 (3 in 100) |
| 4 | 20 (7 in 100) | 548 (7 in 100) |
| 5 | 2 | 72 |
| 6 | 4 | 26 |

**Length.** The 268 are longer tracts: median 26 reference bases against 18, p90 60 against 40.
In repeat units the gap is smaller — median 16 units against 14 — so most of the excess is the
period-2 enrichment, not more repeats.

**Distance from the reference, in repeat units.** This is the surprise:

    same length as the reference   174 of 268
    one repeat away                 16
    two to five repeats away        33
    more than five repeats away     45

**174 of 268 — nearly two in three — are the reference length.** The truth record inside the tract
is a substitution, not a repeat-count change. And **195 of 268 have the missing allele's length
already in ng's merge table**; **180 have it already among the offered candidates.** Whatever else
these are, they are mostly not repeat-length recall failures.

## (b) The finding that matters: 46 tracts where the reads say it plainly

Three worked examples, all verifiable with `samtools view` on the coordinates given.

### chr3:33877690-33877700 — an 11-base poly-A, het, and ng shortens one haplotype

    reference tract   AAAAAAAAAAA        (11 bp)
    truth hap0        AAAAAAAAAAA        (11 bp)
    truth hap1        CAAAAAAAAAA        (11 bp; truth VCF chr3:33877690 A>C, 0|1)

    reads spanning the tract at 30x (n=23)
        13  AAAAAAAAAAA     <- hap0
        10  CAAAAAAAAAA     <- hap1
    at 300x (n=211)
       106  AAAAAAAAAAA     <- hap0
       103  CAAAAAAAAAA     <- hap1

    ng's merge table at this tract
        AAAAAAAAAAA   11 bp   15 reads   kept
        AAAAAAAAAA    10 bp   16 reads   kept

The 10 reads that carry `CAAAAAAAAAA` were counted by ng as a 10-base poly-A: the leading C is
dropped rather than kept as an interruption, and the tract is reported one repeat short. HG002's
second haplotype is the same length as the reference and ng calls it one repeat shorter.

### chr3:37126860-37126871 — a homozygous non-reference tract ng cannot name

    reference tract   AAAAAAAAAAAA         (12 bp)
    truth (1|1)       AAAAAAAAAAAAAGAAA    (17 bp; truth VCF chr3:37126868 A>AAAAAG)

    reads spanning the tract at 30x (n=11)
        10  AAAAAAAAAAAAAGAAA    <- the truth, 17 bp
         1  AAAAAAAAAAAA         <- reference

    at 300x (n=138)
       120  AAAAAAAAAAAAAGAAA
        10  AAAAAAAAAAAA

    ng's merge table
        AAAAAAAAAAAA    12 bp    0 reads   kept (reference)
        AAAAAAAAAAAAA   13 bp   14 reads   kept

10 of 11 reads at 30x, and 120 of 138 at 300x, spell the tract as a 13-base A-run followed by
`GAAA`. ng's table holds a bare 13-base A-run: the run is counted and the `GAAA` tail inside the
tract's span is discarded. The tract is homozygous non-reference and no candidate ng offers can
be right.

### chr11:37147255-37147268 — the same mechanism, and the error variant kept instead

    reference tract   AAAAAAAAAAAAAA       (14 bp)
    truth (1|1)       AAAAAAAAAAAAAGAGA    (17 bp; truth VCF chr11:37147267 A>AGAG)

    reads at 30x (n=14)
        12  AAAAAAAAAAAAAGAGA    <- the truth
         2  AAAAAAAAAAAAAGAGG    <- one base of sequencing error off it

    ng's merge table
        AAAAAAAAAAAAAA     14 bp    0 reads   kept (reference)
        AAAAAAAAAAAAA      13 bp   19 reads   kept, 0.90 share
        AAAAAAAAAAAAAAGA   16 bp    1 read
        AAAAAAAAAAAAAGAGG  17 bp    1 read

ng collapsed 19 of the 21 reads onto a 13-base A-run, and the only 17-base allele it does hold is
the *error* spelling, carried by one read. The true 17-mer that 12 of 14 reads carry is absent.

### How general is this mechanism?

Not universal. Across 267 of the 268 (the 12 kb allele excluded), ng's table does hold an
interrupted allele — one with a base off the motif — at 189 tracts, so ng is not blindly
truncating everywhere. The pattern "the missing allele is interrupted and every read-supported ng
allele at that tract is a clean motif run" holds at 51 of 267, and 17 of those sit in bucket (b).
So interruption handling is one mechanism among several, worth perhaps 17 to 30 of the 268, not
all 67.

### The 14 where ng offered neither the sequence nor its length, at 30x

Each row: tract, period, reference length, missing-allele length, reads in ng's table, and the
reads at 30x that carry the missing sequence exactly out of those spanning the tract.

    chr4:64382192-64382235   per2  44 -> 18 bp   ng 5 reads    17 of 18 reads carry it
    chr3:37126860-37126871   per1  12 -> 17 bp   ng 14 reads   10 of 11
    chr11:37147255-37147268  per1  14 -> 17 bp   ng 21 reads   12 of 14
    chr12:51201535-51201544  per1  10 -> 13 bp   ng 33 reads   11 of 27
    chr8:110988781-110988795 per1  15 -> 19 bp   ng 23 reads    7 of 15
    chr1:160764897-160764918 per1  22 -> 23 bp   ng 19 reads    5 of 13
    chr6:138955441-138955457 per1  17 -> 26 bp   ng 16 reads    4 of 12
    chr6:160587752-160587799 per2  48 -> 56 bp   ng 12 reads    4 of 10
    chr12:64561988-64562023  per4  36 -> 48 bp   ng 15 reads    4 of 10
    chr1:156758593-156758603 per1  11 -> 25 bp   ng 17 reads    3 of 10
    chr3:193742235-193742246 per1  12 -> 22 bp   ng 19 reads    3 of 14
    chr3:197288720-197288744 per1  25 -> 35 bp   ng 10 reads    3 of 5
    chr7:117768566-117768584 per1  19 -> 27 bp   ng 10 reads    3 of 6
    chr5:55443178-55443193   per1  16 -> 30 bp   ng 30 reads    2 of 16

Ten of the fourteen are homopolymers, and in every one the reads report a tract *longer* than the
reference while ng's table does not hold that length. `chr4:64382192` is the loudest: 17 of 18
spanning reads at 30x agree on an 18-base allele where the reference tract is 44 bases, and ng's
table has 5 reads on it.

## (c) 121 where the sequence is not in the reads at all

At these, 300x gives a median 115 spanning reads (minimum 16) and not one carries the truth
sequence. Two shapes:

**27 are a single substitution the reads simply do not show.** The clearest is
`chr1:147531181-147531195`, a poly-A tract. The truth VCF calls `chr1:147531194 A>G` as `1|1` —
homozygous. `samtools mpileup -Q 0 -d 400` over the 300x BAM at that position returns 220 reads,
of which **2** carry G. Whatever HG002's assembly says, these Illumina reads say the position is A.
Other members of the same shape: `chr4:117740548 T>C` 0 of 229 reads; `chr4:147095874 C>T`
0 of 228; `chr7:155206174 G>C` 0 of 143; `chr5:84242812 A>C` 0 of 128.

**The rest are substitution walls written alongside an indel that sits outside the tract.** The
GIAB tandem-repeat truth set represents one assembly-level event as a cluster: an indel plus a run
of neighbouring SNVs. `chr2:183286191-183286210` carries an insertion at 183286190 and nine
alternating `T>A`/`A>T` records inside the tract, all on the same haplotype. Applying only the
records that fall inside `[start, end]` builds a rotated 20-base string that no read carries; the
reads at 300x split cleanly 99 / 40 / 29 between three sequences, none of which is it. 40 of the
121 are walls of 3 or more substitutions with zero read support at at least one of them.

To test whether the *truth* or our tract-local reconstruction of it is at fault, both haplotypes
were rebuilt over `[start-25, end+25]` from every truth record in that wider window and matched
against reads spanning it. At 22 of the 121, both wider haplotypes are carried by 3 or more reads
— so at those the truth is right and the tract-local reconstruction is what failed. At 48 neither
wider haplotype is carried either, so the disagreement is with the truth set itself.

**A concrete defect in the recall script, worth fixing whatever else is done.**
`records_over(truth, starts, chrom, start, end)` keeps a truth record only when its REF span
reaches `start`. An insertion anchored on the base at `start - 1` — which is where left alignment
puts every repeat-length gain at a tract's left edge — has REF span `[start-1, start-1]` and is
dropped. But `haplotype()`'s own docstring says its window begins at `start - 1` precisely so that
such a record *is* applied. The two disagree, and the substitutions the truth set wrote alongside
that insertion get applied without it. 35 of the 268 have such a dropped length-changing record,
and so do 2,737 of the 19,613 tracts where nothing was missing — so it is not the main cause of
the 268, but it does manufacture some of them. Widening the filter to `>= start - 1` changes the
totals (434 misses becomes 564), so it is not a one-line fix; it needs the window and the
overlap rule reworked together.

## The 66 ambiguous ones

Every base the truth needs is present in the reads at 300x, but no read spells the whole tract that
way. Looking at them individually, they are spelling and rotation disagreements at the same length:

    chr6:41399585-41399641   82 of 82 reads   CCTCCT...CTTCTTCTTCCTCCTCCTCCT   (51 bp)
                             truth            CCTCCT...CTTCTTCCTCCTCCTCCTCCT   (51 bp)

    chr14:48430082-48430107  96 of 170 reads  TATATATATATATATATTTT             (20 bp)
                             truth            TATATATATATATATATATT             (20 bp)

    chr10:4503054-4503105    72 of 130 reads at 48 bp; the truth's 48-mer is a rotation of it

Reads and truth agree on the length at 61 of these 66, and ng had already offered that length. If
a tract genotype is scored on repeat length, these are not errors.

## What is worth doing

1. **The 46 unambiguous (b) cases, and inside them the 14 where the length is wrong too.** That is
   about 2 wrong genotypes in every 1,000 tracts at 30x — small, but the reads are unambiguous and
   the failure is ng's. The `chr3:37126860` and `chr11:37147255` examples point at one specific
   thing: how a non-motif base inside the tract span is handled when the tract length is measured.
2. **Nothing in (a) or (c).** 135 of the 268 — half — are an allele no short read can span or a
   sequence the reads do not contain at 300x. Work spent there buys nothing.
3. **Before any of it, settle whether tract genotypes are scored on sequence or on length.**
   174 of the 268 are the reference length and 195 have their length already in ng's table. If the
   score is length-based, three quarters of this bucket is not a bucket.

## Scripts

All under `tmp/agent_unseen/`, all run with `uv run --no-project python`:

* `emit_missing.py` — the recall script's reconstruction, extended to write every miss (and every
  clean tract, as the comparison set) rather than only counting them.
  `emit_missing.py <dump.tsv> <outdir>` writes `missing.tsv` and `clean.tsv`.
* `shape.py`, `shape2.py` — the distributions and the first-cut split.
* `bam_check.py` — reads' tract strings over the tract's reference span.
* `bam_check2.py` — the same over a padded window (stricter; retained for comparison).
* `bam_check3.py` — the two placement-proof tests that decide the buckets.
* `wide_window.py` — truth rebuilt over `[start-25, end+25]` and matched against reads.
* `anchor_check.py` — the dropped left-edge truth record count.
* `purity.py` — motif-interruption test on the missing allele and on ng's table.
* `inspect.py` — everything known about one tract: `inspect.py chr3 33877690 33877700`.

Intermediate tables under `tmp/agent_unseen/out/`; `classified.tsv` is the per-case bucket
assignment with the read evidence at 30x and 300x.
