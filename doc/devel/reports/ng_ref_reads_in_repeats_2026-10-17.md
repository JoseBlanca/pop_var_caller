# Reads that stop inside a repeat are counted as the reference

**Date:** 2026-10-17. Branch `ref-reads-in-repeats`, from main 354754b8. Prompted by a colleague's
GIAB report on pvc 0e905bbb (`tmp/giab_findings/colleague_report_pvc0e905bbb.md` in the main
checkout, sections 3 and 4a). Nothing that matters here changed between 0e905bbb and 354754b8.

**The defect is real, on both calling paths.** At an indel inside a repeat, pvc counts as a
*reference* read some reads that cannot show the indel at all: their sequence is the same whichever
haplotype they came from. At a homozygous indel each such read is evidence for a reference allele
that is not there, and enough of them turn `1/1` into `0/1` at GQ 99. That is what the colleague
saw: at 13 homozygous indels at 300×, pvc counted between 1.7 and 28 times as many reference reads
as cross the indel.

A prototype fix on this branch removes most of those reads: of the four ways they arise (§1), it
answers three and not the fourth. On simulated reads from
four of the colleague's sites, every homozygote it called `0/1` is now `1/1` at 30× and 200×, and
every simulated heterozygote stays `0/1`. On the local GIAB and tomato benchmarks it changes almost
nothing else (§4). It is not merged; §7 lists what needs a decision.

**Words this report uses.**

- **Slide extent.** An indel inside a repeat can be written at several places with the same
  result: inserting `AAAC` into `…AAACAAAC…` gives the same sequence wherever it goes. The slide
  extent is the stretch of reference over which the indel can be placed, from the first base after
  its left-most anchor to the last base it can reach. A read tells the indel's alleles apart only if
  it shows at least one base past the slide extent. The colleague's report calls this the *sliding
  region*.
- **Complete and partial observations.** A read that covers every position of a locus is a
  *complete* observation: its bases there are taken as one allele's spelling. A read that covers
  only some positions is a *partial* one, and calling counts it for every allele it could have come
  from (`read_likelihoods.md` §5.2–5.3).
- **The two paths.** The *STR path* builds a locus over each short tandem repeat (STR) tract the
  repeat catalog lists (`repeat_tract` in `--explain-loci`). The *generic path* builds loci over the
  rest of the reference (`generic`). Region typing cuts the reference into STR tracts and the
  generic regions between them.

---

## 1. How the reads get counted as the reference

There are four ways, two on each path. All four have one root: **the locus ends before the indel's
slide extent does**, and a read that crosses the locus is taken to have seen the indel or its
absence.

| | Path | What happens | Reported site it explains |
|---|---|---|---|
| A | STR | A read stops inside a short array that continues the tract's periodicity past its edge. The array is too short for the catalog to list, so the tract's flank lies inside it. The read matches its flank and is a complete reference observation. | chr22:18094763, 38187308, 42616310 |
| B | STR | The mapper wrote the indel inside that array, not at its left-most place in the tract. The STR path spells the tract from the read's own alignment, which shows the tract as the reference. | the same three |
| C | generic | The record carrying an indel is widened to the indel's slide extent, then cut at the generic region's last base, because a generic region always ends where an STR tract begins. A read that stops (or is soft-clipped) between the cut and the slide end covers the record. | chr22:20516001, 33464883; HG004 chr10:68636979 and chr5:80361671 on local reads |
| D | generic | An insertion anchored on the last base of an STR tract is moved onto the next generic region's first base (`claimed_junction_insertion`), and that record was not widened to the insertion's slide extent. | chr22:36082092 |

At chr22:18094763, for example, the STR locus is the nine-`T` tract at 18094759-767. The reference
then reads `CTTTTCTTTTCTTTTCTTC`, and an inserted `TTTTC` can be placed anywhere up to 18094785, 18
bases past the tract. A read ending inside that stretch shows exactly the reference's bases, whether
or not it carries the insertion.

**36082092 is mechanism D, and the colleague's doubt about it is settled.** Its generic locus
(36082093-101) does start exactly at a region border: a `T10` tract ends at 36082092. The insertion
`GAGA` is anchored on the tract's last base and so is re-homed onto the generic region's first base,
where it kept its own five positions. It was never widened to the `GA` run it slides through. The
locus the colleague saw covering the slide extent was widened by other reads' placements of the same
insertion, which is why its width varied with depth: in the simulation it was 36082093-097 at 30×
and 36082093-103 at 200×.

**The low-depth all-`ALT` cases (colleague's §4a) are not this, as far as can be checked here.** The
two the report names (HG004 chr22:43084962 `CTTTTTTTTT>C` and 45059992 `CAAAAAAAAAAAAA>C`) are
deletions wholly inside a long homopolymer tract, and neither can slide past the tract's edge.
Simulated heterozygotes at both sites at 10× are called `0/1` before and after the change. The other
27 cases are not named in the report; settling them needs the colleague's per-site table (§7).

## 2. How it was confirmed

**None of the colleague's sites fall inside the local GIAB data.** The local per-sample alignments
cover 100 random regions per sample, and none of those reaches the reported chr22 sites.

**So the sites were reproduced from simulated reads.** For each site, 2 kb of GRCh38 either side was
cut out, a haplotype carrying the truth indel was written, 150-base read pairs were drawn from it
with `wgsim` (0.2% substitution errors, no other variants), and the pairs were mapped back to the
cut reference with `minimap2 -ax sr`. pvc then called them with `--defaults --explain-loci`. A
simulated homozygote is the cleanest test there is: **every read carries the indel, so every
reference read pvc counts is a wrong one.**

At 200× the simulated homozygotes reproduce the colleague's symptom at all three STR sites. Each is
called `0/1` at GQ 99, with about a quarter of reads counted as the reference: 47 of 196 at
18094763, 41 of 184 at 38187308 and 42 of 169 at 42616310.

At 18094763 the reads were sorted by where they stop and where the mapper put the insertion, and
each group was called on its own. 18 reference reads come from reads ending inside the array
(mechanism A) and 17 from carriers whose mapper put the insertion inside it (mechanism B). The other
12 appear only when the groups are called together, and were not traced further.

**On real reads, the generic path shows it too.** Among the 119 truth homozygous-alternative indels
in the local GIAB regions, two sit where a generic region ends at an STR tract and the indel slides
into the tract:

- **HG004 chr10:68636979, `CGTTGT>C` sliding into a `(TG)n` tract.** The reference and the deleted
  haplotype agree for two bases after the anchor and first differ at 68636982. The mapper soft-clips
  each carrier read exactly there, so its alignment ends at 68636981, the cut end of the generic
  record. At 300× pvc counts 39 such reference reads where no read crossing the slide extent is
  reference (the generic record is `0/1` at GQ 99), and 3 at 30× (`0/1`, GQ 35).
- **HG004 chr5:80361671, `AT>A` sliding into a `(TTTA)n` tract.** 13 reference reads counted at 300×
  against none crossing. Still called `1/1`.

No local site has an STR tract followed by an unlisted array that a homozygous indel slides into, so
mechanisms A and B could not be seen on local real reads.

**How common the STR geometry is.** For each of chr22's 26,650 STR tracts, the reference alone tells
how far an indel of one to six bases starting inside the tract can slide past either edge. 16,720
tracts (63%) have no such continuation. 1,738 (6.5%) have one of 10 bases or more, and 417 (1.6%)
one of 15 or more, which runs past the 15-base flank the STR path fetches. The four STR sites in the
colleague's table reach 16 to 24 bases past the tract's right edge.

## 3. The fix prototyped, path by path

### 3.1 STR path

**A reference read that stops inside the slide becomes a partial observation** (mechanism A). For
each tract, the reference says how far an indel starting inside it could slide past each edge
(`SlideReach`, in `src/locus_generation/ssr/slide_reach.rs`). An indel of `p` bases slides one base
further for each following base that equals the base `p` before it. A continuation counts only when
it holds for at least one more whole copy, so a single base that happens to match does not count.
A complete observation that spells the tract as the reference, from a read whose alignment stops
before one base past that reach, is kept as a partial observation, anchored on the flank it did
reach. Calling then counts it as "at least this long", which a carrier's read would also show.

**Mechanism B is not fixed.** Moving a read's indel to its left-most place before the tract is
spelled — the left-alignment the generic path applies to every read — was built and measured, and
removed. It rescued none of the 17 mechanism-B reads in the simulation, because each carries a
sequencing error inside the inserted copy that pins the mapper's placement where no exact move
reaches. And it changed no byte of the calls on the GIAB per-sample benchmark at 30× and 300× or
on the 63 tomato accessions. §5 says what is left.

**What it costs, at both ends of the range.**

- **It is the same for every sample of a cohort**, from one sample to thousands, because the reach
  is measured on the reference. Every sample's reads at a locus face the same test.
- **At low depth, few reads change.** On tomato at about three reads a position (8 accessions,
  4 Mb), 390 of 121,011 complete STR observations became partial: 3 in 1,000.
- **What a demoted read loses** is its evidence against longer alleles of the tract itself: as a
  partial it says "at least this many bases", so it no longer counts against an expansion. It still
  counts against contractions. Reads that spell the tract some other way, one copy short say, are
  never demoted: they saw the tract end where a carrier's read would not.

### 3.2 Generic path

**The claimed junction insertion is widened to its slide extent** like every other indel
(mechanism D). A read starting on the region's first base, one base after the insertion's anchor,
is not counted at that record, since it may have started inside the inserted copies.

**A record cut at the region's end keeps the uncut end as an "evidence end"** (mechanism C). The
footprint keeps the cut: the bases past it belong to the next region, as the owner ruled on
2026-09-02. A read that covers the footprint but whose alignment stops before the evidence end is
not counted at that record. It is set aside, not made partial: the cut record spells the indel by
the bases inside the region only, so a carrier's sequence past the footprint is not "the allele
followed by the reference", and the partial rule would find the read compatible with the reference
alone.

Two guards decide which slides set an evidence end. Both came from a false call the first prototype
made (§6):

- **only an indel whose first changed base is inside the region.** A deletion anchored on the
  region's last base removes only the next tract's bases. Left-aligned stutter in that tract is
  anchored there, and must not count.
- **only a slide end that at least two reads in the sample ask for**, so that one read with an odd
  insertion cannot set reads aside for everyone.

**What it costs, at both ends of the range.**

- **It is decided within each sample**, because the generic walk is per sample.
- **At about three reads a position a heterozygote often brings one carrier read, and then nothing
  is set aside.** The defect then stays as it was in that sample; this is the cost of the two-read
  guard.
- **At low depth, few reads move.** On tomato, 16 reads were set aside over 8 accessions and 4 Mb.

## 4. Measured before and after

**Simulated homozygotes, every read a carrier.**

Reference reads counted, and the genotype. At 18094763 the 19 left are the 17 mechanism-B reads and
2 others; mechanism A's 18 are gone.

| Site | Path | Reference reads at 200×, before → after | 200× genotype | 30× genotype, before → after |
|---|---|---|---|---|
| 18094763 `T>TTTTTC` | STR | 47 → 19 | `0/1` → `1/1` | `0/1` → `1/1` (GQ 28) |
| 38187308 `A>AAAAG` | STR | 41 → 15 | `0/1` → `1/1` | `0/1` → `1/1` (GQ 48) |
| 42616310 `A>AAAAG` | STR | 42 → 8 | `0/1` → `1/1` | `0/1` → `1/1` (GQ 42) |
| 36082092 `T>TGAGA` | generic | 11 → 0 | `1/1` | `0/1` → `1/1` (GQ 78) |
| 20516001 `A>AAAAC` | generic | 18 → 10 | `1/1` | `1/1`, GQ 15 → 63 |
| 33464883 `GCA>G` | generic | 14 → 2 | `1/1` | `1/1`, GQ 30 → 99 |

**Simulated heterozygotes** at 10× and 200× at the same sites are `0/1` before and after. At 200× the
counts move toward the even split a heterozygote should show: at 38187308, 117 reference against 85
alternative reads becomes 86 against 85.

**The local GIAB benchmark** (HG002, HG003 and HG004, each over its own 100 regions, at 30× and 300×;
genotypes scored against the truth after normalising both sides). Exactly one genotype changes across
the six runs: HG004 at 30×, chr10:68636979, `0/1` → `1/1`, which is right. Every other count of right,
wrong, missed and false calls is identical. About 40 records change their QUAL or GQ only.

**Tomato** (all 63 accessions at about three reads a position, the first 40 benchmark regions, 4 Mb).
37 of 6.2 million sample genotypes change, in 12 records; 13 records disappear and 12 appear. There
is no truth set, so the changes were compared with GATK's joint calls over 40 of the accessions: 63
changed sample genotypes move toward GATK's call and 71 away. That says the change does little on
tomato, not which side is right.

## 5. What the prototype does not fix

- **Mechanism B: a carrier read whose mapper put the indel inside the array.** In the simulation
  every such read had a sequencing error inside the inserted copy, which pins the mapper's placement
  past the tract. They are most of the reference reads left in §4's table: 19 at 18094763, 15 at
  38187308 and 8 at 42616310 at 200×. With mechanism A's reads gone they no longer change the
  genotype there.
- **The same thing on the generic path at chr10:68636979 at 300×.** 19 carrier reads ending about ten
  bases into the `(TG)n` tract were written as `138M1I9M` and the like, with the insertion just past
  the cut. That generic record stays `0/1`. It does not reach the VCF: the record pvc writes there is
  the STR tract's edge deletion, which was `1/1` at GQ 99 before and after.
- **Mechanisms A and B on the colleague's real reads are untested.** Only simulated reads reproduce
  them here.

## 6. What went wrong on the way, so it is not tried again

The first prototype had two more rules. Both made false calls and both are removed.

- **On the STR path: "an indel past the tract that cannot be moved into it means the read's own
  spelling is untrustworthy, so let the repeat delimiter spell it, and set the read aside if that
  says reference."** On tomato SL4.0ch06:14526115 the reads of one haplotype carry real insertions 1
  and 7 bases past a `T9` tract. Their own alignment spells the tract correctly as the reference, but
  the rule set every one of them aside. The tract was left with no reference read, and a false
  `T9 → T9GT` was called at allele frequency 0.98, at a site where GATK has no variant. Across tomato
  this version moved 423 sample genotypes away from GATK's call and 221 toward it.
- **On the generic path: "a read reaches only as far as its first indel past the region's end."** At
  GIAB HG004 chr10:4848077, the base before a `T17` tract, every left-aligned stutter deletion of the
  tract is anchored on that base. The rule set aside 244 of 252 reference reads, and a false `TC>T`
  was written at QUAL 6.

## 7. Decisions for the owner

**1. Whether to merge this as the fix.** I recommend merging the three parts on this branch — the STR
demotion (A), the generic evidence end (C) and the junction insertion's slide (D) — once two things
that are mine to do are done: unit tests for each, and one existing test's fixture.
`a_locus_on_a_new_chromosome_mints_a_cursor_and_keeps_the_old_tallies` builds its tract on a contig
that is all `A`, so the tract slides to both contig ends and its one read is now set aside. The test
needs plain sequence around its tract, not a change to the rule. The trade-off of merging is the
two-read guard on the generic path: at about three reads a position it often leaves the defect in
place for a heterozygous sample, in exchange for never setting reads aside on one odd read's say-so.

**2. Whether to go further for the residual reads of §5.** The only way I see to recover a carrier
whose error pins the mapper's placement is to compare the read against each candidate haplotype:
realignment. The STR path's own delimiter does that over the tract and its 15-base flank, but a
placement past the flank is out of its reach, and letting it re-spell such reads is the rule that
made the false tomato call of §6. The other way is to make region
typing extend an STR tract over the short array that continues it, so the tract's locus covers the
whole slide extent. That changes where every affected region border falls, on every sample and every
reference. **I recommend neither for now.** The residual reads no longer change a genotype in the
simulation, and the colleague's real reads (decision 3) would show whether they do in practice.

**3. The colleague's files, to test this on the reported sites.** psp files are not enough, because
the change acts while psps are written. What is needed:

- their chr22 alignments for HG002, HG003 and HG004 at full depth (the 300× tier) and the 10× draw
  (seed 42), with indexes;
- their `het_for_hom_reads.tsv` (§3 of their report) and `hom_for_het.tsv` (§4), which list the 38
  and 29 sites with the read counts to compare against.

## 8. Reproducing

All scratch is under the worktree's `tmp/`, not committed:

- `tmp/repro/`: `make_sites.py` cuts the sites and writes the alternative haplotypes;
  `simulate.sh` draws and maps reads; `run_pvc.sh` calls them; `summarise.py` and `compare.py` read
  the `--explain-loci` output; `split_by_bucket.sh` calls one group of reads at a time.
- `tmp/real/`: `run_real.sh` calls the local GIAB samples; `informative.py` counts reads crossing
  each truth indel's slide extent beside pvc's counts; `score_gt.sh` scores genotypes against the
  truth; `slide_past_edges.py` is the chr22 survey of §2.
- `tmp/tomato/`: `run.sh` calls the cohort; `psp_census.sh` writes and dumps psps for the read
  census; `gatk_cmp.py` compares changed genotypes with GATK's.
