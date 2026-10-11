# ng — a simulated stuttering benchmark: STR segments of known genotype, read through libraries that stutter as real ones do

*Design spec, 2026-10-09. For the owner's review before any code. It is milestone 2 of
`tract_dual_path.md` (branch `tract-ceiling`) §13, built here because the colleague's simulation is
weeks away (owner, 2026-10-09). Reads on: `tract_dual_path.md` (branch `tract-ceiling`) §2 (what the
real benchmarks cannot show), [`read_likelihoods.md`](read_likelihoods.md) (the stutter model ng
calls with), and the existing simulator
[`examples/ng_tract_simulator.rs`](../../../../examples/ng_tract_simulator.rs), which this extends.*

---

*Note, 2026-10-10: `tract_dual_path.md`, the design this benchmark was built to test (calling each STR
segment both ways), lives on branch `tract-ceiling`, which the owner decided not to merge after this
benchmark's first pass
([report](../../reports/implementations/str_stutter_benchmark_first_pass_2026-10-10.md)). The
benchmark itself is kept on main; references to that spec below point to that branch.*

---

## 1. What this is

A program that makes a cohort whose STR genotypes we chose, sequences it through libraries that
stutter at rates taken from the literature and from this project's own fits, and hands the reads to
three callers — ng, GATK HaplotypeCaller and HipSTR — so that every call can be scored against the
truth.

**Stutter** is the copying error at a repeat: while a library is amplified, the polymerase
sometimes adds or drops one whole copy of the motif, so a read reports an STR one unit longer or
shorter than the DNA it came from. It is why ng calls each STR segment twice and weighs the two
answers by how much the reads stutter.

**Why it is needed.** Neither benchmark the caller is measured on can say whether that weighing
works where stutter is real. GIAB's human libraries are PCR-free and stutter little; the tomato
libraries do stutter, but there is no truth for their STR lengths (GATK's indel calls in repeats are
not admissible as truth there). So the STR model's side of the design has never been tested
(`tract_dual_path.md` §2), and three numbers in it are judgements waiting on this benchmark:
the 0.95 weight that hands a boundary-crossing segment to the STR answer (§8), the cost of reading
an STR record's allele counts in the duplication filter, and whether one strongly stuttering
library can be outvoted by many clean ones in a large cohort.

**What it answers, per caller, at every combination of stutter level, period, repeat count, depth
and cohort size:**

1. How often each sample's genotype at an STR segment is right — the whole sequence, and the length
   alone.
2. How often a SNP inside or beside an STR segment is called right.
3. How many false calls each makes at STR segments, and how confident they are.

## 2. What is simulated, in words

- **A reference genome made for the purpose**: contigs of ordinary sequence with STR segments laid
  into it — every motif length from one to six bases, at repeat counts from the calling floor up to
  where stutter stops rising — each with at least 150 bases of non-repetitive sequence either side,
  so the region typing calls each one an STR segment with clean flanks.
- **A cohort drawn from a population**: at every STR segment a set of alleles and their
  frequencies, and each sample's two copies drawn from those frequencies. Alleles differ in length
  by whole repeat units, and some carry a substitution inside the segment. Beside the segments,
  SNPs and short indels, some of them deletions running from a segment into its flank.
- **Libraries that stutter**: each sample has a library, and each library a stutter level. A read
  copies one of its sample's two alleles; it may then gain or lose whole units at each STR segment it
  crosses, at the library's rate for that motif length and repeat count; then sequencing errors.
- **Paired-end reads over the whole genome**, at the depth asked for, so reads start anywhere: some
  span a segment whole and many end inside it, as on a real run.

## 3. Stutter: the model and its levels

### 3.1 The model

For each read crossing an STR segment of motif length `p` and `n` repeats, in a library at level `L`:

1. **Whether it stutters**: with probability `s(p, n, L) × m`, where `s` is the level's rate for
   that motif length and repeat count (§3.2) and `m` is the library's own multiplier (§3.3).
2. **Which way**: shorter with probability `c(L)`, longer otherwise.
3. **By how much**: one unit with probability 0.85, two with 0.85 × 0.15, three with 0.15² × 0.85,
   and so on, capped at four units.

A read crossing several segments stutters at each independently. **This is the shape of ng's own
read model**, which is the point of a test: a caller whose model matches the data should do well,
and §3.3 adds what the model does not assume.

### 3.2 The rate, by motif length and repeat count

**Stutter rises with the number of repeats and falls with motif length**, in every dataset. The rate
at the calling floor is about one thirtieth of the rate on long segments; it rises evenly on a log
scale to a plateau, reached at 16 repeats for mono- and dinucleotides and at 10 for longer motifs,
and stays there.

The plateau `s_max(p, L)`, the share of reads that stutter on a long segment:

| motif length | none | PCR-free | PCR | heavy PCR |
|---|---|---|---|---|
| 1 | 0 | 0.065 | 0.19 | 0.30 |
| 2 | 0 | 0.015 | 0.15 | 0.25 |
| 3 | 0 | 0.005 | 0.05 | 0.08 |
| 4 | 0 | 0.003 | 0.015 | 0.03 |
| 5 | 0 | 0.003 | 0.015 | 0.03 |
| 6 | 0 | 0.002 | 0.015 | 0.03 |
| share shorter, `c(L)` | — | 0.70 | 0.78 | 0.80 |

**Where the numbers come from:**

- **Mono- and dinucleotides, PCR-free and PCR: this project's own fits**, which estimate the rate
  per motif length and repeat count on real reads. GIAB HG002–HG004 at 30× (PCR-free): homopolymers
  2 reads in 1,000 at 8 repeats rising to 66 in 1,000 at 17 and beyond; dinucleotides about 1 in
  1,000 at 6 and 7 repeats, the only lengths it saw. Four tomato accessions (PCR libraries):
  homopolymers 5 in 1,000 at 8 repeats to 190 in 1,000 at 16 and beyond; dinucleotides 7 in 1,000
  at 6 repeats to 147 in 1,000 at 15 and beyond. Shorter share: 0.68 on GIAB, 0.73 and 0.82 on
  tomato. And **51 tomato accessions measured directly** (each read's length against the
  sample's most common length at the locus, 46,800 loci a sample; project memory, 2026-07-26):
  homopolymers 0.3% below 10 bases, 6.5% at 10–14, 22.7% at 15–19; dinucleotides 0.9% at 10–14
  bases, 6.6% at 15–19, 12.6% at 20–29; trinucleotides 1.3% at 20–29 bases and 5.1% at 30 and
  beyond; tetra- and pentanucleotides under 0.5% everywhere; hexanucleotides about 1%.
- **Longer motifs, and the PCR-free dinucleotide plateau: Gymrek lab (2014)**, the share of
  stuttered reads over all loci of each motif length in one human sample sequenced with PCR
  (HGDP00542, 26×) and one without (38×): with PCR 0.17, 0.038, 0.011, 0.0069, 0.0074 and 0.014
  for motif lengths one to six; PCR-free 0.027, 0.011, 0.0029, 0.0016, 0.0015 and 0.0010. Those are
  averages over loci of every length, so the plateaus above are set at about twice them for motif
  lengths four to six, and for trinucleotides at the tomato measurement's 5%.
- **The fold between PCR and PCR-free**: 4.6 averaged over motif lengths (Gymrek lab), and up to 9
  (Fungtammasan et al. 2015, plasmid data, dinucleotides). The table's PCR-free to PCR fold is 3 to
  10 by motif length.
- **Direction and size**: losing one unit is the most common error at every motif length, then
  gaining one, then larger changes (Fungtammasan et al. 2015); the PCR-free data contract slightly
  less (Gymrek lab). The 0.85 one-unit share is a judgement consistent with that ordering; neither
  source gives the number.
- **Heavy PCR** is not measured anywhere above: it stands for a low-input or many-cycle library, set
  at 1.6 to 2 times PCR, to see where each caller breaks.

### 3.3 What varies between libraries, and what ng's model does not assume

- **Libraries differ.** Each library's multiplier `m` is drawn around one (log-normal, a standard
  deviation of 0.4 on the log scale, so 95 libraries in 100 sit between 0.45 and 2.2 times their
  level's rate). On the 51 tomato accessions the share of stuttered reads ran from 0.37% to 1.98%
  between samples, about fivefold, far beyond what the read counts allow by chance; a standard
  deviation of 0.4 gives that spread over 51 draws. The spread partly followed each library's
  duplicate rate, a sign of PCR (rank correlation 0.48).
- **One outlier library**: a cohort design (§4) where every library is PCR-free except one at heavy
  PCR. This is the case `tract_dual_path.md` flags as a known limit — the stuttering library's
  evidence outweighed by many clean ones.
- **Not modelled** (§7): stutter that depends on the flanks or on interruptions inside the repeat,
  partial-unit slips, PCR chimeras, and stutter that differs between the two alleles of one sample.

## 4. The cohort, the depths and the loci

### 4.1 Designs

| design | libraries | why |
|---|---|---|
| clean | level none | a control: every STR error is the callers' own |
| PCR-free | PCR-free, multipliers drawn | GIAB-like |
| PCR | PCR, multipliers drawn | tomato-like |
| heavy | heavy PCR, multipliers drawn | where stutter dominates |
| one outlier | PCR-free, one sample heavy PCR | the large-cohort limit |

**Cohort sizes: 1, 10 and 50 samples. Depths: 5, 15, 30 and 60 reads a position.** That spans the
range the caller is built for at both ends that can be run in hours: one sample is the hardest case
and must be in it; 50 is the size of the tomato comparison; 5× and 60× bracket the depths both
benchmarks have. Thousands of samples are out of reach for GATK on this schedule and are left out.

**The full grid is 5 designs × 3 cohort sizes × 4 depths = 60 runs per caller.** I propose a first
pass of 16 — PCR-free and PCR, at 1 and 50 samples, at 5, 15, 30 and 60× — and the rest only where
the first pass shows a difference worth placing.

### 4.2 The loci

Per run, on about 4 Mb of simulated genome:

- **3,000 STR segments**: 500 per motif length, their reference repeat counts spread evenly from the
  calling floor (8, 6, 6, 6, 5, 4 copies for motif lengths one to six) to 30 copies for mono- and
  dinucleotides and 15 for the rest, motifs drawn from the common ones of each length.
- **Alleles at each segment**: two to six, by a draw; each differs from the reference by a whole
  number of units, one or two units most often (a geometric step around the reference); a minor
  allele frequency drawn from a Beta distribution so most variation is rare and some common, as in
  a population. **One segment in ten** has an allele carrying a substitution inside the repeat;
  **one in ten** a SNP within 30 bases of its edge; **one in thirty** a deletion that begins in the
  segment and ends in the flank, or the reverse.
- **Genotypes**: each sample's two alleles drawn from the frequencies, with an inbreeding
  coefficient — 0 by default, and 0.9 as an option for a selfing crop like tomato.
- **Away from the segments**: SNPs at one per thousand bases and short indels at one per ten
  thousand, drawn the same way. They are controls, and they give ng's parameter fit the ordinary
  ground its census needs.

### 4.3 The reads

Paired-end, 2 × 150 bases, fragment lengths drawn around 350 (standard deviation 50), positions
uniform over the genome, one read group a sample. Sequencing errors: substitutions at the rate the
base quality states, qualities drawn around Q32 with a falling tail along the read. The simulator
writes FASTQ, and **minimap2** maps it (§8).

## 5. The callers, and how each is run

- **ng**, both ways it can be run: `call-from-alignments --defaults`, and `generate-psps` +
  `estimate-parameters` + `call-from-psps` with the fitted parameters. The fitted run is the one
  that counts; `--defaults` shows what the fit buys.
- **GATK 4.6 HaplotypeCaller** in its standard joint pipeline, already used in this project's
  tomato benchmarks (`benchmarks/tomato1/scripts/perf_gatk_joint.py`): per sample
  `HaplotypeCaller -ERC GVCF`, then `CombineGVCFs` and `GenotypeGVCFs` over the cohort, with default
  settings and no filtering.

- **HipSTR v0.6.2**, the STR specialist whose read model ng's STR path follows, already built in
  this project (`HipSTR/HipSTR` in the main checkout) and run on tomato by
  `benchmarks/lib/run_hipstr.sh`: all samples jointly, over the simulator's STR segments given as
  its region BED, default settings. It calls STR segments only, so it is scored on them alone
  (§6); its calls carry SNPs inside a segment as part of the STR allele, which the scoring reads.

## 6. Scoring

**The unit is one sample at one STR segment.** For each, the truth is the sample's two sequences
over the segment and the base before it. A caller's answer is built from its records overlapping
that stretch, the same way `tract_dual_path.md` §6 compares two answers: each allele's difference
from the reference written as edits slid to their left-most placement, counted per copy. So an STR
record and two generic records saying the same thing score alike.

- **Genotype right**: the called edits equal the truth's. **Length right**: the two called lengths
  equal the truth's, whatever the sequence — the measure an STR-only tool is usually judged by.
- **No call**: a sample uncalled, or no record where the truth has a variant.
- **False call**: a variant called where the sample carries none — counted with its genotype
  quality, so each caller's confidence in its errors can be read off.
- **SNPs and indels away from segments, and the SNPs inside them**: RTG vcfeval against the truth
  VCF, as on GIAB, restricted to the segments and to everything else separately.

Every count is broken down by design, motif length, repeat count (in bands), depth and cohort size,
and given as counts, not only shares. HipSTR is scored on the STR segments only; ng and GATK on
everything. ng's `--explain-loci` output is kept, so every stretch where
ng's two answers disagree can be checked against the truth: whether the weight favoured the right
model.

## 7. What this does not test

- **Alignment beyond what §8 decides**: with an aligner, its own misplacements are in; without one,
  every read is placed where it came from.
- **Stutter that depends on anything but motif length, repeat count and library** (§3.3), interrupted
  repeats and satellites. STR bundles — repeats a few bases apart — were added on 2026-10-11 (§11).
- **Real genomes' sequence context**: GC content, mappability, paralogs. The duplication filter will
  see a genome with no duplications, so its effect here is only its false drops.
- **Cohorts of thousands** (§4.1).

A caller that does well here has done well on the problem as these numbers state it; one that does
badly has a problem that real data will not hide.

## 8. Alignment: minimap2 (owner, 2026-10-09)

The simulator writes FASTQ and **minimap2** maps it (`-ax sr`), because minimap2 is the aligner this
caller's production data are mapped with (owner). Every aligner habit at repeats is then in the
data — indels shifted to one end of a repeat, reads clipped rather than gapped, reads ending inside
a segment placed with no anchor — which matters because ng's generic path takes alignments as given
while its STR path re-aligns. Cost: `minimap2` added to the dev container (one package, and the
image rebuild of about ten minutes), and minutes of mapping a run. Writing each read's true
alignment instead was considered and set aside: it is the easiest version of the problem.

## 9. How it is built

1. **The simulator**: a new program, `ng_str_stutter_benchmark`, beside `ng_tract_simulator`,
   taking its deterministic draw and its BAM writing; settings as `key=value` (design, samples,
   depth, seed). It writes the reference, the reads as FASTQ, the truth VCF, each
   sample's truth sequences at every segment, the segments' BED and the stutter each library was
   given.
2. **One driver script** that runs a configuration end to end: simulate, map with minimap2, sort and
   index, run ng both ways, GATK and HipSTR, score — and one that runs the grid. `minimap2` goes
   into the Containerfile first.
3. **The scorer** (§6), as a program in this repository.
4. **The report**: per configuration, the tables of §6 for both callers, then a short account of
   what the grid shows.

**Unit tests** on the draw: the stutter rate a library produces, counted over many simulated reads,
within its sampling error of the table; truth sequences reproduced from the truth VCF.

**Cost**: a run at 50 samples and 30× is about 6 Gb of reads; GATK's per-sample GVCF step is the
slowest part, roughly minutes a sample on a few Mb. The first pass of 16 runs is a few hours.

## 10. Sources

- Gymrek lab, *Illumina PCR-free protocol greatly reduces stutter noise at STRs* (2014):
  stutter share by motif length, with and without PCR, on HGDP00542 —
  <https://gymreklab.com/science/2014/04/04/stutternoise.html>.
- Fungtammasan A. et al., *Accurate typing of short tandem repeats from genome-wide sequencing data
  and its applications*, Genome Research 25:736 (2015): PCR-containing against PCR-free error rates,
  error direction and size — <https://pmc.ncbi.nlm.nih.gov/articles/PMC4417121>.
- Willems T. et al., *Genome-wide profiling of heritable and de novo STR variations* (HipSTR),
  Nature Methods 14:590 (2017): the per-locus stutter model with a geometric step size —
  <https://www.nature.com/articles/nmeth.4267>.
- This project's fitted parameters: GIAB HG002–HG004 at 30× (`tract_dual_path` milestone C fits)
  and four tomato accessions (the portable-float measurements of 2026-09-14), the
  `slippage_by_stratum_and_group` rows.

## 11. STR bundles (added 2026-10-11)

*Built for [`str_bundle_loci.md`](str_bundle_loci.md) §9–§10 step 1, so that calling a bundle can be
measured under PCR stutter and in a cohort.*

**An STR bundle is two or more repeats within 15 bases of each other**, so that none of them has the
15 repeat-free bases either side an STR segment needs; the typing calls the stretch from the first
repeat's start to the last one's end — the **hull** — one bundle, and each repeat in it a
**member**. The simulator now lays bundles beside the segments.

### 11.1 What is laid

- **How many**: the setting `bundles` (default **500**, one for every six segments at the default
  3,000; `0` lays none). With `bundles=0` every file the simulator wrote before is written byte for
  byte as before, for the same seed and settings: everything about the bundles — their shapes,
  alleles and genotypes — is drawn from a second random stream, seeded from the run's seed, so the
  first stream makes exactly the draws it made before. Checked on 2026-10-11 by running the build
  before bundles and this one on one configuration (2 samples, 600 segments, 0.6 Mb) and comparing
  the checksums of all twelve files, reads decompressed; the unit test
  `no_bundles_draws_the_cohort_drawn_before_bundles` pins the truth VCF and the reference.
- **Members**: two, or three in 15 bundles in 100 (HG002's tandem-repeat regions have about one in
  eleven; set higher so three-member bundles are numerous enough to count). Each member's motif
  length is uniform over one to six, its motif drawn from the segments' list, its reference repeat
  count uniform from the calling floor (8, 6, 6, 6, 5, 4) to the segments' highest (30 for motif
  lengths one and two, 15 for the rest) — at or above the floor, so the typing counts each member
  as a repeat.
- **Gaps**: 1 to 14 bases between consecutive members, uniform. Each gap base is drawn so it
  continues neither neighbour's motif, as the spacer beside a segment does: the first bases differ
  from the next bases a copy of the left motif would have, the last from the bases a copy of the
  right motif would have before it; and no three copies of any unit of one to six bases end on a
  gap base.
- **Hull**: at most 100 bases — the default `--max-str-bundle-len` of `str_bundle_loci.md` §5.2 —
  with 150 bases of clean flank either side, spacer included, as a segment has.
- **The typing must see exactly the members.** A shape is kept only if the simulator's own copy of
  the scanner's arithmetic (each base scored +2 when it equals the base one motif length before it,
  −7 otherwise; maximal scoring stretches; the copy floors) finds the members and nothing else. This
  is what refuses two members of one motif length — whatever their motifs, since the scanner
  compares bases a motif length apart — whose gap is too short to keep them apart: the scanner would
  read them as one repeat. The copy lives in the simulator so the fixture does not change when the
  caller's code does.

The rules are applied by redrawing the members and gaps (never their number), so they shape the
distribution: with the default seed (the shapes do not depend on the genome's size), of 500 bundles 81 have three members; members by motif length
one to six are 248, 185, 235, 167, 143 and 103 (the hull cap thins long motifs); 93 of 581 pairs of
neighbouring members share a motif length; the hull's median is 77 bases, its 90th percentile 96.
**That is more than twice HG002's median of 35 bases**, because each member's repeat count is drawn
uniform from the calling floor to the segments' highest — 30 copies for motif lengths one and two,
15 for the rest, so a hexanucleotide member can be 90 bases — and a member averages 31 bases (1,081
members at the default seed). A long hull leaves fewer reads that cross it whole, so at low depth
these bundles are harder than HG002's; whether to draw counts nearer the floor is open.

**Alleles**, as for a segment (§4.2): two to six, the non-reference ones sharing a Beta-drawn
frequency. A length change moves one member by whole units of its own motif (one most often), and in
one allele in four a second member too; no member falls below one copy. Derived alleles, at the
segments' shares: one bundle in ten gains a substitution inside a member of at least three copies,
one in ten a SNP within 30 bases of the hull's edge, one in thirty a deletion of two to ten bases
across the hull's left or right edge, cutting into the first or the last member.

**Stutter**: one stutter site per member in each copy of the genome, each at its member's own rate
— its motif length and its repeat count in that allele (§3.1). A molecule crossing a bundle
stutters in each member independently.

### 11.2 The files

The bundles have **files of their own**, beside the segments' three, so that the segments' files are
unchanged (a kind column on them would change them at `bundles=0`) and every reader of them — the
driver, `--explain-loci`, HipSTR's region file, the report scripts — still reads segments only.
Per-member values are listed in member order, comma-separated; coordinates are zero-based half-open.

- `bundles.bed` — a `#` header line, then `chrom start end name periods motifs reference_repeats
  member_starts member_ends`: the hull, the bundle's name (`bun00000`), then each member's motif
  length, motif, reference repeat count and span.
- `bundle_alleles.tsv` — one row per allele, as `segment_alleles.tsv`: `bundle chrom window_start
  window_end hull_start hull_end members periods motifs reference_repeats allele repeats change
  frequency sequence`. The window is the hull and 30 bases either side; `repeats` is the allele's
  count in each member; `change` takes the segments' values.
- `bundle_genotypes.tsv` — `bundle sample first_copy second_copy`, indices into the allele table.
- `truth.vcf` — one record per bundle window any sample varies at, spelled as a segment's (trimmed to
  the stretch its alleles differ over), with ID the bundle's name and INFO `BUNDLE`, `MOTIFS`
  (`Number=.`), `CHANGE` and `MEMBER_REPEATS` (`Number=R`, each allele's member counts joined by
  `_`, e.g. `8_10,9_10`). The three new keys are declared in the header only when the file has a
  bundle.

### 11.3 The checks and the scoring

- **Typing** (`scripts/ng_str_stutter_check_bundles.py`, run by the driver's typing check): every
  simulated bundle must be typed as exactly one `ssr_bundle` overlapping its hull, with as many
  members, in order, each of the drawn motif length and with a start and an end within a tolerance
  of the drawn member's — **0 bases by default**. The typed members are the repeat detector's
  intervals, which differ from a drawn run only when a gap base continues a motif, and the gap rule
  stops that; on one run (the default seed, 600 segments and 500 bundles on 1.5 Mb) all 500 matched with no offset, and the
  typing made no bundle the simulator did not lay. The segments are checked as before. The typing
  dump (`examples/ng_typed_region_dump.rs`) gained a `members` column for this, a bundle's members
  as a JSON array of `{"start","end","period"}` as `pvc type-regions` writes them.
- **HipSTR is not given the bundles.** A HipSTR region carries one period and HipSTR fits one stutter
  model to it, so a bundle given as its hull would have every member's stutter priced at one
  member's period; its calls would measure that, not HipSTR. GATK and ng are scored on the bundles.
- **Scoring**: a cell is one sample at one bundle, scored exactly as a segment's cell (§6) — the
  scorer compares sequences over the window and needs no member structure. Bundle cells go to
  `score_bundle_cells.tsv` and `score_bundle_summary.tsv`, so `score_cells.tsv` and
  `score_summary.tsv` hold segments only, as before. One motif length cannot describe a bundle, so it
  is summarised by its **number of members** and its **longest member** — the one with the most
  reference bases — with that member's motif length and repeat band (§6). vcfeval gains a third
  stretch, `bundles`, and `outside` now leaves out the bundles' windows as well as the segments'.
