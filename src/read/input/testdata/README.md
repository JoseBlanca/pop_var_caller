# `multi_contig_slice.*` — a CRAM block holding reads from several chromosomes

**Why this file is committed rather than built by a test.** The CRAM writer this project vendors
puts one chromosome in every block, so nothing in the tree can produce this shape. `samtools` can,
and does so on its own.

**What it is.** 24 short chromosomes (`mr00`…`mr23`, 400 bases each, pseudo-random from a fixed
seed) with three 60-base reads on each — one read per chromosome carrying a planted substitution,
the rest exact copies of the reference under them. Coordinate-sorted.

**Written with samtools 1.16.1 and no writer options at all**, which matters: the merging is
samtools' own default behaviour, not something the fixture asked for. It merges blocks across
chromosomes once several in a row would be under-full (`htslib/cram/cram_encode.c`, the
`c->curr_rec < c->max_rec/4+10` test), so a reference of many short contigs — a draft assembly, a
scaffold-level reference — produces these as a matter of course. Two chromosomes alone are not
enough to trigger it; 24 are.

The result, read back from the `.crai` by grouping its lines on (container offset, landmark):
**4 blocks, 2 of them spanning several chromosomes, the largest holding 20.**

```
samtools faidx multi_contig_slice.fa
samtools view -C -T multi_contig_slice.fa -o multi_contig_slice.cram multi_contig_slice.sam
samtools index multi_contig_slice.cram
```

**The reference is pseudo-random rather than a run of one base**, for the same reason A2's fixture
is: against an all-`A` reference a read rebuilt at the wrong offset comes out identical, and no
assertion can fail.

**`multi_contig_slice.sam` is the ground truth**, not a by-product: it is the input samtools was
given, so comparing decoded reads against it is independent of anything ng or noodles does. The
tests that use it are `a_cram_block_spanning_several_contigs_rebuilds_the_reads_its_sam_holds` and
`a_short_window_on_a_block_spanning_several_contigs_is_refused`, both in
`../open_bam.rs`.

**How a `.crai` reveals such a block**, since the obvious check does not: the reference id `-2`
that the block's own header carries never reaches the index. htslib writes one index line per
chromosome instead, all sharing the container offset and the landmark
(`htslib/cram/cram_index.c`, `cram_index_build_multiref`). So group the lines and count distinct
ids; searching for `-2` finds nothing however many there are.
