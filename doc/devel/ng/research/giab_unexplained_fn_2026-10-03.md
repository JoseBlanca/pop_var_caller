# Missed GIAB indels nobody had explained: where each true allele is lost

**Date:** 2026-10-03
**Branch:** `giab-unexplained-fn`. This report is the only file the branch adds.
**Caller:** `call-from-psps` at commit `25bfe313`, the commit the user benchmarked.
**Data:** GIAB HG002, HG003 and HG004 at 300 reads a position, the `per_sample` benchmark's 100
regions per sample, scored against GIAB v4.2.1.

## The answers

| case | where the true allele is lost | what kind of problem | recommendation |
|---|---|---|---|
| **1.** 31-base deletion, chr10:4845453 (HG004) | The merge refuses to build the locus. The deletion is one copy of a two-copy, 31-base tandem duplication, so it can be placed anywhere along 62 bases. The locus is therefore 62 bases wide, past the 50-base limit (`--max-cohort-locus-span`). | A deliberate limit, set at 50 bases and never measured. | Raise the default to 100 once a large low-depth cohort shows what it costs there. On this run, 100 recovers the deletion and changes nothing else. |
| **2.** CCTTC repeat tract, chr1:206838736 (HG002) | The periodicity check refuses the tract: 8 of its 55 spanning reads (15%) have lengths off the 5-base grid, against an allowance of 10%. Five of the eight end only 2 to 7 bases past the tract, inside reference sequence that is itself CCTTC-like, so the tract length measured from them is wrong. | A defect in our code: reads with no unique anchor at one end are counted as spanning. The 10% allowance it trips was inherited and never measured. | Count a read as spanning only if it extends past the flank by enough bases to place the tract's end uniquely. Name refused tracts by position in the run report. |
| **3.** The three sites the paralog filter removed | **Not the filter.** The walk that writes the psps caps reads at 250 per read group at positions where some read has an insertion or deletion, but at 1,000 elsewhere. At these sites the cap thins the deletion reads to about 56%, while the reference reads are counted at the next positions, where the 1,000 cap applies. A 47–50% heterozygote is stored as 31–37%, and at 1.2–1.5 times the usual depth the filter reads that as a duplication. | A defect in our code. It is also the cause of the "known defect" the filter's D3 report put down to alignment bias. The truth set is right: no GIAB difficult-region list contains the sites, and every read has mapping quality 60. | Apply one cap to every position of an indel's record. Dropping the separate 250 cap does it: on this run it recovers the three sites and one more true deletion, and costs nothing in file size or time. |
| **4.** "Same variant written another way" | Not lost. At five truth indels, GIAB writes an indel and a nearby SNP as two records. We write one repeat-tract record whose alleles carry both. vcfeval matches all five. Exact position-and-allele matching misses all five. | A scoring artefact, not a calling error. | Score with vcfeval, or with hap.py's `--engine=vcfeval`. No caller change. |

**What happened next is in §5.** Case 3's fix is done and measured (branch `one-read-cap`). For
case 2, every rule tried cost more right genotypes on the HG002 tandem-repeat benchmark than the
one tract it recovers, so none is committed. Refused tracts are now listed by position in the run
report (branch `refused-tracts-listed`). In case 2, the +6 lengths turned out to be the mapper's
spelling, which the caller trusts; its own re-alignment measured the truth.

Two corrections to the brief, both from the data:

- **GIAB calls the three paralog-filter sites heterozygous (0/1), not homozygous.** Its own read
  counts across all its datasets put the variant in 47–50% of reads.
- **The run report does list case 1's locus as too wide.** Reproduced at `25bfe313`, the report
  reads `loci the merge declined to assemble for being too wide: 1 — chr10:4845453-4845514 (62
  bases)`, both over the full regions and over a 6-kb window. If the user's report said 0, their
  run differed from this one in some way not identified here.

## What was run

**The reads at every case site come from one sample only.** Each sample's alignment file holds
reads only over its own 100 regions, and the regions do not overlap between samples. So a joint
call over the three is, at every site in this report, a one-sample call with two samples that
have no reads.

`call-from-psps` refuses a cohort whose psps were walked over different ground. So every sample
was walked over the union of the three BED files: 299 intervals, 1,425,377 bases. The commands
were `generate-psps` for each sample's alignment file, then `call-from-psps --defaults` on the
three psps. A parameters file fitted on 1.4 Mb would mean little. The filter's coverage fit is
made inside the calling run either way, and it fitted one copy at 301, 306 and 311 reads per
window for HG002, HG003 and HG004.

The dev loop used one 6-kb window per case (four windows, 20,528 bases). Every case-1 to case-3
number below was checked in both the windowed run and the full one. The exceptions are stated
where they occur.

**Scoring.** hap.py could not be run on this machine: its container image is built for x86-64
only, and this Mac has no x86 emulation. vcfeval (RTG Tools 3.13) was run per sample against
GIAB v4.2.1, restricted to that sample's own regions. Beside it, an exact match was taken: both
files normalised with `bcftools norm`, then intersected on position, reference and alternative
allele. This is the method of `benchmarks/giab/src/score_ng_recall.sh`. Over the full run,
vcfeval counts 29, 5 and 22 missed variants (SNPs and indels together) in HG002, HG003 and HG004.

**The user's read counts reproduce exactly:**

| site | user | this run |
|---|---|---|
| chr10:4845453 deletion, once the locus is built (case 1) | 61 deletion, 83 reference | 61, 83 |
| chr1:206838736 tract, reads at the two truth lengths (case 2) | 24 and 20 | 24 and 20 |
| chr1:101051327 (case 3) | 233 / 103 | 233 / 103 |
| chr3:125575780 (case 3) | 199 / 118 | 199 / 118 |
| chr3:125576513 (case 3) | 253 / 114 | 253 / 114 |

To see inside the caller, temporary prints were added at three places. The first is the
periodicity check, `locus_is_periodic` in `src/calling/allele_candidates/ssr.rs`. The second is
the merge's per-sample allele derivation, `alleles_of_sample` in `src/run/cohort_merge/build.rs`.
The third is the walk's per-position read cap, in `src/locus_generation/pileup/genome_walk.rs`.
One experiment changed a constant, `DEFAULT_MAX_INDEL_COLUMN_DEPTH`. None of this is committed.
The patch is at `tmp/probes.patch` in the branch's worktree.

## 1. The 31-base deletion at chr10:4845453

### Where it is lost

The deletion removes `GGGCTGCCAGCATGCTGTCACCCCTCAATGG`. In the reference that sequence occurs
twice in a row, at 4845454–4845484 and 4845485–4845515:

```text
4845440 AGTGAGGGTTGTGAGGGCTGCCAGCATGCTGTCACCCCTC
4845480 AATGGGGGCTGCCAGCATGCTGTCACCCCTCAATGGGGGC
```

Deleting either copy gives the same haplotype. So the reads place the deletion anywhere from its
leftmost position to its rightmost, and the locus that holds every placement runs from 4845453
to 4845514: **62 bases**. The merge refuses any ordinary (non-repeat-tract) locus wider than
`--max-cohort-locus-span`, 50 by default (`too_wide` in `src/run/cohort_merge/close.rs`). The
refusal is counted and listed in the run report. Nothing is written.

The same arithmetic holds for any deletion of one copy of a two-copy duplication. The locus is
twice the deletion's length, so at the default limit **no such deletion longer than 25 bases can
be called.**

### The evidence

Re-calling the same psps with `--max-cohort-locus-span 70`:

```text
chr10  4845453  AGGGCTGCC…CCCCTCAATGGGGGCTGCC…CCCCTCAATG  AGGGCTGCC…CCCCTCAATG
       QUAL 741.5  PASS  HG004 0/1  AD 83,61
```

The genotype is GIAB's (0/1), and vcfeval matches it.

Raising the limit on the full run cost nothing measurable:

| `--max-cohort-locus-span` | loci refused as too wide | records written | vcfeval missed / false, all three samples |
|---|---|---|---|
| 50 (default) | 1 (this one) | 2,451 | 56 / 53 |
| 100 | 0 | 2,452 | 55 / 53 |
| 200 | 0 | 2,452 | 55 / 53 |

### The four leads in the brief

None of them is involved, because the locus is refused before any of their code runs:

- **reads removed by the merge** act on a locus being built, and this one is not built;
- **a deletion cut at a region's edge** would give a shorter deletion. The record arrives whole:
  62 bases, and called correctly once allowed;
- **a called locus with no variant genotype written**: the locus is never called;
- **the depth ceiling** (1,000 reads per read group): HG004 has 277 reads at the site.

### Classification and recommendation

**A deliberate limit, and the measurement it was waiting for now exists for one corner.** The
code calls the 50-base default "the owner's number, unmeasured and soft", and says the
measurement that would settle it is how much real signal sits just above 50. On three samples at
300×, one true deletion sits above 50, and nothing false does.

That is not the whole range. A wide locus is mainly a risk at low depth in a large cohort: there,
overlapping observations from many samples chain into one long locus. **Recommendation: raise
the default to 100 once a run on the 63-sample tomato cohort at about 3× shows how many loci the
change admits, and that calling them does not add false calls.** Re-calling stored psps under a
new limit needs no new walk, so this is cheap.

## 2. The CCTTC repeat tract at chr1:206838736

### Where it is lost

The catalog's tract is chr1:206838725–206838799: 75 bases, 15 copies of `CCTTC`. GIAB's two
alleles, at 206838736, are deletions of 5 and of 10 bases, i.e. one and two whole units. **So
both of the brief's suspicions are refuted.** The truth alleles are whole numbers of units, and
the catalog's tract is in phase: 75 bases is exactly 15 units.

The periodicity check refuses the tract. It asks each sample what share of its spanning reads
differ from the reference tract by something other than whole 5-base units. If every sample is
over 10%, the tract is refused. HG002, the only sample with reads, has 55 spanning reads:

| reads | length change | what they are |
|---|---|---|
| 22 + 2 | −5 | truth allele (two of the reads carry one substitution) |
| 18 + 2 | −10 | truth allele (two with one substitution) |
| 1 each | −20, −25 | on the grid |
| 1 | −35 | on the grid |
| 4 + 1 | **+6** | **off the grid**: tract read past its end (below) |
| 1 each | −6, −8, −27 | **off the grid**: one read each with a ragged end |
| 0 | 0 | the reference length (GIAB: 0 reference reads) |

**Off the grid: 8 of 55 reads, 15%, against an allowance of 10% (5.5 reads).** The tract is
refused, called as the reference allele alone, and so nothing is written. The run report says
only `not called — the reads do not vary in whole motif units (notPeriodic): 1`.

### Why five reads are 6 bases too long

The reference continues past the tract's end with sequence close to the motif, impure enough
that the catalog did not include it:

```text
206838770 CCTTCCCTTCCCTTCCCTTCCCTTCCCTTC | CTCTCTCCTTCCCTTCCTTTCTC…
                       tract ends at 799 ^   the 15-base right anchor is CTCTCTCCTTCCCTT
```

The five +6 reads are the five reads that start before the tract and end only 2 to 7 bases past
it, by the mapper's alignment. The mapper gave each a 6-base insertion just before its end
(`130M6I12M`, `128M6I12M2S` twice, `125M6I17M` twice). Each read's last 12 to 17 bases stand
where the tract's 15-base right anchor should be, and that anchor sequence is itself CCTTC-like,
so the mapper could not place the end of the tract in these reads.

**The +6 is the mapper's spelling, not the caller's measurement.** The caller re-aligns every
spanning read against the tract and 15 bases of flank each side, and that re-alignment measured
all five reads at −10, the truth (a temporary print in the walk, 2026-10-03). But a complete read
whose mapper alignment carries an insertion or deletion is *spelled* by that alignment, which
needs to reach only one base past each end of the tract; the re-alignment then decides only
whether the read is complete. That rule was adopted by the tract-accuracy program because the
mapper's account keeps junction variants the re-alignment destroys. Here it kept the mapper's
misplacement instead. The catalog only guarantees that no other *catalogued* tract lies within 15
bases of this one, not that the 15 bases are unlike the motif.

Without those five reads, the off-grid share is 3 of 50, 6%, and the tract passes.

### What calling it would have produced

With the check bypassed, the same psps give:

```text
chr1  206838725  <75-base tract>  <tract −10>,<tract −5>  QUAL 573.9  PASS
      HG002 1/2  AD 0,18,22
```

That is GIAB's genotype: one haplotype missing two units, the other missing one. It does not need
an ordinary indel locus. The repeat-tract path calls it correctly once the check lets it through.

### Classification and recommendation

**A defect in our code, compounded by an unmeasured threshold.** The defect is that a read whose
end lies inside an anchor that looks like the motif is spelled by its mapper, which could not place
it. The threshold, one read in ten, is documented as inherited from the deleted caller and "never
measured, by them or by us".

*The recommendations below were written before the repairs were measured. §5 has what happened.*

Recommendations, most useful first:

1. **Count a read as spanning only when it reaches far enough past the tract to place the
   tract's end uniquely.** For example: a minimum number of matched bases beyond the anchor, or
   a margin between the best and second-best placement of the tract's end. At this tract, the
   five +6 reads end 12 to 17 bases into a 15-base anchor and would not count.
2. **List refused tracts by position in the run report**, as too-wide loci are listed. Today a
   refused tract is a count of one in a report covering 1.4 Mb.
3. Measure the 10% allowance. Fix 1 removes this tract's reason to need it.

## 3. The three sites the paralog filter removed

### Where the allele is lost

**In the walk that writes the psp, not in the filter.** The walk caps how many reads of one read
group it uses at one position. The cap is 1,000, and **250 at any position where some read has an
insertion or deletion**: `DEFAULT_MAX_INDEL_COLUMN_DEPTH`, taken from samtools' indel depth limit
(`column_depth_cap` in `src/locus_generation/pileup/genome_walk.rs`). At 300× a deletion's anchor
position holds about 430 reads, so the 250 cap applies there. It keeps 250 reads chosen by a hash
of the read name, so it thins both alleles alike: 183 deletion reads become 103.

The record does not take its reference count from that position alone. The reference allele is
also seen at the positions the deletion spans, where no read has an indel and the 1,000 cap keeps
everything. **So the deletion keeps 56% of its reads, while the reference keeps nearly all of
its own.** The column-by-column counts, from the temporary print in the walk:

| site | position | reads there | cap | reads kept | deletion reads before → after |
|---|---|---|---|---|---|
| chr1:101051327 | anchor | 428 | 250 | 250 | 183 → **103** |
| | next base | 236 | 1,000 | 236 | — |
| chr3:125575780 | anchor | 378 | 250 | 250 | 174 → **118** |
| | next base | 203 | 1,000 | 203 | — |
| chr3:125576513 | anchor | 478 | 250 | 250 | 220 → **114** |
| | next base | 256 | 1,000 | 256 | — |

The psp record at chr1:101051327 then holds 103 complete reads with the deletion and 233 without.
The merge removes no read at any of the three sites: each is a single record for its sample, and
the merge consults individual reads only when a sample has several records in one locus.

How the three sites look at each stage:

| site | fragments in the alignment file, deletion / none | stored with the 250 cap | stored with one 1,000 cap | filter's probability of duplication: 250 → 1,000 |
|---|---|---|---|---|
| chr1:101051327 (HG002) | 185 / 210 — **47%** | 103 / 233 — **31%** | 183 / 233 — 44% | 1.000 → passes |
| chr3:125575780 (HG004) | 175 / 181 — **49%** | 118 / 199 — **37%** | 174 / 199 — 47% | 0.993 → passes |
| chr3:125576513 (HG004) | 217 / 220 — **50%** | 114 / 253 — **31%** | 220 / 253 — 47% | 1.000 → passes |

Fragments are counted by read name among reads that cross the deletion with at least 5 bases on
the left and 8 on the right; no fragment's two mates disagree. The aligner placed the deletion
consistently: 184 of 186, 174 of 178 and 220 of 221 deletion reads carry it at the same position.

### Why the filter then removes them

The filter compares two explanations of a site. Either it is an ordinary variant in a sample with
two copies of the region, or the region is duplicated in the sample and collapsed onto one place
in the reference. A heterozygote in three copies out of four puts the variant in a third of the
reads, and puts 1.33 times the usual depth on the site. These sites have 1.43, 1.22 and 1.54 times
their sample's fitted one-copy depth (429, 380 and 478 reads, against 301 and 311). With the
stored fraction at 31–37%, the duplication explanation wins. With the true fraction, 44–47%, it
does not. The depth excess alone does not carry the verdict.

### Is it the truth set?

No. GIAB calls all three 0/1, from 47–50% variant reads across its datasets. None of the three
sites lies in any of these GIAB GRCh38 stratifications (v3.0):

- segmental duplications;
- the union of all difficult regions;
- the union of low-mappability regions and segmental duplications;
- the regions where collapsed duplications cause false positives;
- HG002's and HG004's own CNV and structural-variant regions.

Every read at all three sites has mapping quality 60.

### This is the cause of the recorded "known defect"

`doc/devel/reports/implementations/ng_paralog_filter_d3_2026-09-07.md` found that at 300× the
filter removes about nine real indels for every one it should. It found the removed indels'
variant fractions clustered at 0.32–0.37. It explained that as alignment bias: "reads carrying an
indel are harder to place". **These three sites contradict that explanation.** The alignment
file has the deletion in 47–50% of fragments, and the third is made by our cap. D3 measured on a
larger HG002 file at the same depth, so its 19 removed indels were very likely thinned the same
way. That is an inference, not a re-measurement of D3's file.

### What one cap does to the whole run

The same full run, re-walked with `DEFAULT_MAX_INDEL_COLUMN_DEPTH` raised from 250 to 1,000. That
makes the cap 1,000 at every position:

| | 250 at indel positions (shipped) | 1,000 everywhere |
|---|---|---|
| indels the filter removes | 6 | 3 |
| SNPs the filter removes | 86 | 86 |
| vcfeval missed, HG002 / HG003 / HG004 | 29 / 5 / 22 | 27 / 5 / 20 |
| vcfeval false, HG002 / HG003 / HG004 | 25 / 14 / 14 | 24 / 14 / 14 |
| psp sizes | 18.15 / 13.67 / 17.32 MB | 18.11 / 13.63 / 17.27 MB |
| calling time | 4 s | 4 s |

Of 335 indel calls (one per sample per site), 24 gain more than 5 percentage points of variant
reads when the cap is lifted, and none lose. Two genotypes change. One is a fourth recovered
truth variant, below. The other is chr1:243535155, a false call either way.

**The fourth variant is one of the homozygous deletions the brief puts down to partial reads.**
chr1:21684679, a 6-base deletion GIAB calls 1/1, is written 0/1 with read counts 0 reference and
232 deletion. With one cap it is written 1/1 with 0 and 360. The 250 cap thins the reads that
carry the whole deletion. Reads that end inside the deletion are recorded at positions the cap
does not reach, so they keep their full number, and the partial-read scoring then counts them as
reference. The cap is therefore a second cause of that pattern, beside the partial-read scoring
another session is fixing. The two fixes are independent, and either one alone moves this site.

### Classification and recommendation

**A defect in our code, in the walk that writes the psps.** Neither the filter nor the truth set
is at fault.

**Recommendation: one cap for every position an indel's record covers.** The simplest form drops
the separate 250 cap and uses `--max-reads-per-position` everywhere. Its stated reason, that indel
evidence in a homopolymer saturates early, is about work saved. On this run it saved none: file
sizes and calling time are unchanged. If a tighter cap at indels is wanted, it has to be applied
across all of the record's positions, so that every allele is thinned alike.

After the fix, D3's measurement should be re-run on its 1,000-region HG002 file, and its "alignment
bias" explanation corrected. Three indels are still removed at one cap (chr7:85066911, chr16:88779362
and chr22:41687198). The first is a GIAB variant with 53 of 219 reads carrying the insertion; it
was not traced here.

## 4. The "same variant written another way" cases

### What they are

Five truth indels are missed by exact matching and matched by vcfeval. Each is the same shape.
GIAB writes an indel and a SNP a few bases away as two records. We write one repeat-tract record,
and its alleles carry both changes:

| sample | GIAB | ours |
|---|---|---|
| HG003 | chr11:123865154 `CGTGTGTGT`→`CGTGTGT`,`C` (2/1), and 123865171 G→A (0/1) | chr11:123865155, 38-base `GT` tract, two alleles: one with the A substitution, one shorter (1/2) |
| HG004 | chr8:24760294 `TAC`→`T` (0/1), and 24760296 `CACAT`→`C` (0/1) | chr8:24760292 `TATACAC`→`TATAC` (1/1), and a 2-base deletion in the `AT` tract at 24760299 (0/1) |
| HG004 | chr14:90749402 `CTT`→`CT`,`C` (2/1), and 90749422 T→G (0/1) | chr14:90749403, 20-base `T` tract: one allele 2 bases shorter and ending in G, one 1 base shorter (1/2) |
| HG004 | chr10:4853442 `CA`→`C` (0/1), and 4853458 A→C (1/1) | chr10:4853443, 16-base `A` tract: one allele with the C, one also a base shorter (1/2) |
| HG004 | chr20:48207407 `GCACA`→`GCACACA`,`G` (2/1), and 48207435 A→G (0/1) | chr20:48207408, 28-base `CA` tract: one allele 2 bases longer, one 4 shorter and ending in G (1/2) |

In every case the two haplotypes are the same sequence, and vcfeval, which compares haplotypes,
counts them as correct. Exact matching of normalised records misses all five: normalising does
not split a multi-base allele into an indel and a separate substitution.

**Which three of these five the user's hap.py run counted could not be checked here.** hap.py
does not run on this machine (see *What was run*). Its default comparison engine matches some
re-written variants and not others, so its count can be below five. The case-4 genotypes are all
correct by vcfeval.

### Classification and recommendation

**A scoring artefact. The calls are right.** Writing a repeat tract as one record is what the
VCF spec (`doc/devel/ng/spec/vcf_output.md`) decides: the merged locus is written whole and shared
bases are not trimmed. A repeat tract that has both changed length and taken a substitution is
one event in our model.

**Recommendation: score GIAB runs with vcfeval, or with `hap.py --engine=vcfeval`.** That makes
these five correct, and no caller change is needed. Splitting tract records into an indel plus
substitutions on output would only help a comparison tool that does not compare haplotypes. It
would also make the record no longer say that the substitution and the length change are on the
same allele.

## 5. What was done afterwards (2026-10-03)

The owner asked for two of the recommendations: one read cap at every position (case 3), and, for
case 2, counting a read as spanning only when it reaches far enough past the tract, with refused
tracts listed by position.

### Case 3 — one read cap: done, branch `one-read-cap`

The walk's separate 250-read cap at positions with an indel is removed; `--max-reads-per-position`
applies everywhere. Measured with the branch's own binary:

| | before | after |
|---|---|---|
| this report's run, vcfeval missed / false (three samples) | 56 / 53 | 52 / 52 |
| D3's run (HG002, 1,000 regions): true GIAB variants among what the filter removes | 19 of 280 (6.8 in 100) | 4 of 265 (1.5 in 100) |
| … indels the filter removes, of which GIAB variants | 18, 17 | 3, 2 |
| tomato oracle (4 accessions, about 3×): calls with default parameters | | byte-identical |
| … with fitted parameters | | 1 fitted number moves by 2 in 100,000; no genotype or filter changes; largest QUAL change 0.1 |

The D3 report carries a correction: its "alignment bias" was this cap. One case stays open and is
recorded in spec `locus_generation_pileup.md` §4: the same imbalance returns where an anchor
position is deeper than the cap and the positions the indel spans are not. At the defaults that
needs more than 1,000 reads in a read group, and both calling commands drop such loci
(`--max-read-group-depth`, 1,000).

### Case 2 — no rule shipped; refused tracts listed, branch `refused-tracts-listed`

Four rules were tried on GIAB's HG002 tandem-repeat benchmark (36,497 truth records), at 30× and
50×. Each recovers this tract. Each costs more elsewhere:

| rule | right genotypes gained or lost, over the four cells (30×/50×, homopolymer/period 2+) |
|---|---|
| use the mapper's spelling only with at least 8 bases of overhang past the tract | −12 (−9, −1, −1, −1) |
| … only with the whole 15-base flank of overhang | −48 (−21, −6, −10, −11) |
| count a read with less than 8 bases of overhang as partial, not spanning (the rule as recommended above) | −142, and 62 homopolymer tracts no longer called at 30× |
| use the re-alignment's spelling where the mapper's is off the motif grid and the re-alignment's is on it | +14 (0, +6, 0, +8), with 25 and 21 more wrong calls at 30× and 50×: tracts whose truth is a sequence change, which the periodicity check had been refusing |

8 bases is the least that covers all five reads here (their overhangs are 1 to 7). A rule of 4
bases changed almost nothing (−1, +1, 0, +2) and does not recover the tract. None of the four is
committed; the experimental code is at `tmp/sweep_rules.patch` in the branch's worktree. **The owner accepted case 2 as a known miss** (2026-10-03): the spelling rules stay, and the run
report's list of refused tracts will show on a larger run how often it happens.

The second half of the recommendation is done. The run report now names refused tracts by
position, the first five in genome order, under the count it always printed. Sorted at printing,
so it is the same at any thread count. On this report's run it lists `chr1:206838725-206838799`
and `chr10:68636982-68637019`. The second is a 38-base TG tract in HG004 where 201 of 254 reads
are 3 bases short; GIAB has a 5-base deletion beside it, `chr10:68636979 CGTTGT→C` (1/1), which
vcfeval also counts as missed. It was not traced here.

## Reproducing this

All paths are in the branch's worktree, `../pop_var_caller-giab-unexplained`, under `tmp/`,
which git ignores:

- `tmp/union.bed` — the union of the three samples' region BEDs; `tmp/windows.bed` — the four
  6-kb case windows;
- `tmp/run.sh OUTDIR [call-from-psps flags]` — `generate-psps` for the three samples if
  `OUTDIR/psps` is missing (`BED=` selects the ground), then `call-from-psps --defaults`;
- `tmp/score.sh CALLS.vcf OUTDIR` — vcfeval per sample; `tmp/exact.sh OUTDIR` — the exact-match
  comparison beside vcfeval's missed indels;
- `tmp/probes.patch` — the temporary prints and the cap experiment. With it applied, the prints
  are turned on by `PVC_DEBUG_PERIODIC`, `PVC_DEBUG_SKIP_PERIODIC`, `PVC_DEBUG_LOCUS=<position
  substrings>` and `PVC_DEBUG_COLUMN=<positions>`.

The stratification BEDs were downloaded from GIAB's
`release/genome-stratifications/v3.0/GRCh38/` into `tmp/strat/`.
