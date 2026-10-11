# ng — an STR bundle called as one STR locus, with stutter in each of its runs

*Design, 2026-10-11. The direction (an STR bundle path with a size cap) and the points in §11 are
the owner's decisions of 2026-10-11; the rest is this design's proposal. Reads on: [`typed_regions.md`](typed_regions.md) §2.4 (what a bundle is),
[`locus_generation_ssr.md`](locus_generation_ssr.md) §3–§4 (how a read is measured across an STR
segment), [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md), [`read_likelihoods.md`](read_likelihoods.md)
§4, [`str_edge_deletion_alleles.md`](str_edge_deletion_alleles.md), [`vcf_output.md`](vcf_output.md),
[`str_stutter_benchmark.md`](str_stutter_benchmark.md); the measurement that decided it is
[`../impl_plan/str_bundles_measurement.md`](../impl_plan/str_bundles_measurement.md) §4.*

---

## 1. The problem

**An STR bundle is two or more repeats within 15 bases of each other**, so none of them has the 15
repeat-free bases on each side that an STR segment needs (`typed_regions.md` §2.4). The bundle's
region is the **hull** of its repeats, so the bases between repeats lie inside it, and the hull as a
whole does have clean flanks. Each repeat the typing finds in a bundle is a **member**.

**Nothing calls a bundle today.** The typing emits them, the walk has no generator for them, and
calling sets them aside (`locus_generation.md` §2, the owner's ruling of 2026-09-12: "a wrong call
in a repeat cluster is worse than a missing one"). On GIAB HG002's tandem-repeat benchmark regions
that leaves 2,422 truth variant records in 2,473 bundles uncalled, against 8,534 in 20,183 STR
segments: about one variant per bundle against 0.4 per segment.

**Size.** 2,207 of the 2,473 bundles have two members, 229 three, 37 more. The median hull is 35
bases, 90% are 71 or shorter, 2,410 are 100 or shorter, 10 are longer than 150.

## 2. The decision

**Each bundle of at most 100 bases is one STR locus.** A read that crosses the hull is spelled over
it, anchored on the clean flanks outside it, exactly as a read is spelled over an STR segment. Each
allele is that spelling. What a bundle locus adds to an STR segment is that **its sequence is split
into runs, each of one motif** (§4, §6.1), so that:

- stutter is priced in every run, in that run's own motif and at that run's own rate;
- each allele has a repeat count per run, for the prior and the VCF.

**A bundle longer than 100 bases is typed as a satellite** and not called, for the reason a single
repeat longer than 100 bases is one (`--max-str-len`): with 150-base reads, a read crosses a repeat
and a usable anchor each side only up to about 100 bases. **The limit is the user's to change**, as
the repeat criteria are, with its own option, `--max-str-bundle-len` (§5.2).

Every other part of the STR path is reused: the per-sample tallies, the edge-deletion rule at the
hull's two outer edges, the candidate bar, the EM loop, the record and its gathering.

## 3. Why not the alternatives

**The generic path over the bundle** was built and measured (the plan's §4). On HG002 (one sample,
PCR-free) it looked good: 901 of 1,125 bundles with a truth variant right at 30×, 16 false calls in
1,342 without one. It fails where the library stutters. Calling the stuttering benchmark's STR
segments with the generic path instead of the STR path, 50 samples, PCR, 30×:

| | variant cells right | false calls, of ~87,800 cells with no variant |
|---|---|---|
| STR path | 61,370 of 62,127 | 79 |
| generic path | 58,477 of 62,290 | 9,097 |

8,279 of those 9,097 false calls are on mono- and dinucleotides of 16 copies or more, and 293 of
HG002's bundles hold such a run. The generic path takes stutter reads for an allele, as GATK does on
the same reads (9,448). On the simulator's bundles (§9; 50 samples, PCR, 30×) it makes a false call
in 756 of 14,543 bundle cells with no variant, 690 of them in bundles with such a run — against 61
of 88,330 for the STR path on segments in the same runs.

**Each member as its own STR segment** is ruled out by the anchor. The STR path aligns a read
against the repeat and 15 repeat-free bases either side (`locus_generation_ssr.md` §4,
`flank_bp ≤ bundle_threshold`); inside a bundle that flank is partly the neighbouring member. A
length change in the neighbour moves the anchor; a read one `CA` short at `(CA)8 GT (CA)10` is
evidence for both loci at once; and a deletion across the gap is the edge-deletion problem between
two STR loci.

**HipSTR's design** — the whole cluster genotyped jointly as one haplotype — is what this is, within
the STR path (`typed_regions.md` §10 named it as the only option that recovers the loci).

## 4. Vocabulary

- **Member.** One repeat the typing found in a bundle: a span and a motif (§5.1). Members are what
  make a bundle a bundle; they are **not** where stutter is placed, because they overlap and nest
  (§5.1, the review's finding).
- **Hull.** The bundle's region: the earliest member start to the furthest member end.
- **Run.** A stretch of pure, whole copies of one of the members' motifs, at least the typing's copy
  floor for that motif's length (`--min-copies`) — where a slip of the polymerase can happen. The
  hull's **reference runs** are found once, from the reference sequence (§6.1).
- **Fixed sequence.** Every hull base outside the runs. No stutter is modelled in it.
- **Run counts.** For an allele, how many copies of each reference run's motif it carries at that
  run, in order — `(20, 8, 8)` for `(CT)20 (AT)8 T8`. The bundle's counterpart of an STR segment's
  repeat count.
- **Run parse.** How an allele's sequence is split into its runs and fixed sequence (§7.1).

## 5. Typing

**5.1 A member carries its motif.** Today a bundle carries its members as bare detector intervals —
start, end, period — split off before the steps that give a segment its motif and its trimmed span
(`segment_criteria.rs`, `classify` step 3). The catalog rows already hold each member's trimmed span
and motif and drop them on the way (`repeat_catalog/segments.rs`, `admit`). The change: each member
goes through the same end-trim and motif extraction an STR segment does, and keeps the result —
`RegionKind::SsrBundle { members: [BundleMember { start, end, motif }] }`. Members are **not** put
through the copy floor, the purity gate or the flank test: membership stays what `typed_regions.md`
§2.4 defines, and only the description of each member is richer. Both routes to the typing — the
scan and the catalog — produce the same members, and the existing differential between them checks
it.

*Changed during implementation (2026-10-11).* Two things this section left open. **A member with no
whole-motif cut** — an imperfect repeat with no two clean copies near one end — keeps its detected
span; it is still a member, and its motif is the same either way. **The hull stays the detected
spans' hull** — from the earliest detected start to the furthest detected end — so a member whose
ragged end was cut lies strictly inside it. On HG002's regions, 385 of the 2,410 bundles end past
every member's cut end, and 7 start before every member's cut start. Keeping the hull is what leaves
every other region of the typing unchanged.

*Found in review (2026-10-11); settled by the owner the same day.* **Members are not ordered,
disjoint repeats.** The repeat finder tolerates impurities, so one member can run through sequence
that is really another repeat. On HG002's 2,410 bundles, **551 hold one member wholly inside
another** and **183 more have consecutive members that overlap**. Two of them:

- chr1:8,007,443, `(CT)20` then `ATATATATATATATAT` then `TTTTTTTTT`: typed as a `CT` member over all 64 bases and a `T` member inside
  it; the `AT` run is no member at all.
- chr1:28,487,365, `AAAAAAAAATAAAAATAAAATAAAAATAAA`: an `A` member over bases 1–9 and an `AAAAAT`
  member over bases 5–30.

**Decision (owner, 2026-10-11): stutter is placed in runs found from the reference sequence, and the
members serve only as the list of motifs to look for** (§6.1). The typing is left as built.

**5.2 The cap.** A hull longer than `--max-str-bundle-len` becomes a `Satellite`. No re-absorption
follows: anything within 15 bases of a member would have joined the bundle, so nothing is within a
satellite's reach that was not already inside the hull.

- **Its own option, default 100** (owner, 2026-10-11: the user can change the limit, as they can
  change what counts as an STR). Not `--max-str-len` itself: the two answer the same read-length
  question today, but a user who widens what a single repeat may be — for long reads, say — need
  not also widen bundles, whose model is newer and less measured.
- **It is a typing setting, so it travels where the other repeat criteria travel**: the three
  subcommands that type (`generate-psps`, `call-from-alignments`, `call-from-psps`), the typing a
  psp records, and the typing block of the parameters file (`parameters_file.md`, beside
  `max_str_len` and `bundle_threshold`). Psps typed with different limits are different typings,
  and are refused together as any two typings are.
- **0 turns bundle calling off**: every bundle becomes a satellite and none is called, which is
  today's result for calls.

*Changed during implementation (2026-10-11).* Where the cap is recorded, in detail:

- **The repeat catalog's own file does not record it.** Its header predates the cap and the builder
  never applies it, so every catalog reads back with the default 100 as its build setting; two
  catalogs never differ on it, and no existing catalog has to be rebuilt.
- **A psp records it** (`max-str-bundle-len-bp` in the header's repeat criteria). A psp written before
  the cap has no such key and is **refused, naming the key**, rather than read as 100: its bundles
  were typed with no cap, a typing no value of the key describes. Such psps must be written again.
- **The parameters file's key is optional on reading**: a file written before the cap says nothing
  about it and is compared on the other axes only. The census digests of the routing criteria and of
  the catalog's build settings are taken over their printed form, which now carries the cap, so they
  move (the cross-platform checksums of the parameters file moved for this and nothing else). **So a
  parameters file fitted before the cap is demoted to `supplied` when a run with a census uses it**
  (`parameters_file.md` §6), although its numbers were fitted from the same loci: bundles never
  enter the census (§6.3).
- **A windowed read of the catalog reaches as far as the cap** (`reach_into_window`). A bundle is a
  chain of repeats, longer than any of them; a read padded only by the longest tract cut a chain
  short and could type a satellite as a bundle (found in review, 2026-10-11).

**5.3 A contig's ends.** A bundle whose hull has no flank on one side — it starts on the contig's
first base or ends on its last — is dropped, as an STR segment is (`FlankClamped`), and counted as
that rejection, as a segment is. Its bases become generic, as today.

*Changed during implementation (2026-10-11).* Before this, such a bundle was kept. The cap is asked
first, so a hull over the cap at a contig end is a satellite, which needs no flank. Both routes
apply the rule. It rarely fires on the catalog's, which keeps a row only with `min_flank_bp` bases
beside its cut span, but a hull is drawn over detected spans, which can reach further.

**5.4 The tallies** (`TypedRegionCounts`) count capped bundles as satellites; `ssr_bundles` counts
what remains.

## 6. Reads and the psp file

### 6.1 The hull's reference runs, and measuring a read

**The reference runs are found once per bundle, when its locus is built**, from the hull's reference
bases and the members' motifs:

1. For each distinct motif among the members, find every maximal stretch of the hull that repeats
   with that motif's period and whose unit is that motif or one of its rotations (`CT` also matches
   `TCTC…`). Its copies are its length divided by the period, rounded down; a part copy at its end
   stays in the run, as a reference segment keeps its part copy.
2. Keep the stretches with at least the typing's copy floor for that period (`--min-copies`,
   default 6, 4, 4, 3, 3, 3 copies for periods 1 to 6).
3. Choose runs **longest first** (in bases; a tie goes to the earlier start). A run that overlaps
   one already chosen gives up the shared bases, cut back to whole copies, and is kept if it still
   reaches the copy floor. The chosen runs never overlap.

The two examples of §5.1 come out as `(CT)20`, `(AT)8`, `T8` — the `T` where the `AT` run and the
`T` run meet goes to the longer `AT` run — and as `A9` with the rest fixed sequence (its `AAAAAT`
units reach two clean copies, below the floor of three for period 6). A bundle can come out with no run at all
— a cluster of impure repeats, none with a stretch at the copy floor. It is still called, with no
stutter modelled anywhere in it; how many there are is counted on HG002 at step 3.

The bundle generator is then the STR segment's generator run over the hull. **Three things differ:**

- **The delimiter** — the aligner that finds where in the read the repeat starts and ends — is the
  **flat-gap delimiter** (`alignment/ssr_best_path_flat_gap.rs`, `SsrFlatGapAligner`), which needs
  only the hull's boundaries, not a motif (`alignment/mod.rs`, `RepeatGeometry`). The default
  delimiter for segments prices whole-unit slips of one motif and cannot describe a hull of
  several. Whether a per-run version of it earns its keep is measured, not assumed (§10, step 3).
- **The edge-escape check** (`reference_run_escapes_the_tract`) uses, at each edge of the hull, the
  motif of the run that reaches that edge, and none where fixed sequence reaches it.
- **Only reads that cross the whole hull are used.** An STR segment also uses a read that ends
  inside the repeat, as a lower bound on its length. In a bundle, such a read gives a lower bound on
  one run's count and says nothing about the runs it did not reach; using it needs a censored
  emission over run counts, which this design leaves out (§11).

Everything else is unchanged and motif-free: the spelling of the read over the hull, the flank
deletion at the two outer edges in its smallest form (`str_edge_deletion_alleles.md` §4), the
per-sample tally keyed on spelling, flank deletion and read group.

### 6.2 The psp file

**`LocusKind::SsrBundle` gets a payload**: the two flanks, as `SsrDetail` has, and the reference runs
as `(offset in the hull, length, motif)`. The members are not stored: once the runs are found, nothing
downstream reads them. Its tag (`KIND_SSR_BUNDLE = 2`) exists and no file has ever carried a bundle
record, so no reader has an expectation to break — but a reader built before this change would read
the tag and no payload, so the format version must make such a reader refuse a file with bundle
records (`psp_file_format.md` §4.3). The typing a psp records changes too, so new psps are needed for
any run that is compared.

The typing dump (`typed_regions_cli.md` §3.1) writes each member's motif in its `members` cell,
beside its start, end and period. *(As built, 2026-10-11: the start and end are the member's cut
span, 0-based half-open as `pvc type-regions` writes its rows, where before they were the detected interval.)*

### 6.3 The parameter census

**Bundles stay out of the stutter census**, as they are today (`census.rs`, the `SsrBundle` arm).
A run is priced with the stutter fitted for STR segments of its period and repeat count (§7.3). This
keeps the fit's ground as it is, and needs no new stratum. (The fit is not run in this work.)

## 7. Calling

### 7.1 The run parse, and run counts

**Every allele is parsed into the reference runs and the fixed sequence once, when it becomes a
candidate.** The parse is the split that spells the allele as: each reference run as whole copies of
its motif (plus the part copy the reference run has at its end, if any), in the reference's order;
the fixed sequence between them as the reference has it; and the flanks as the allele's flank
deletion leaves them — with the fewest mismatching bases. Ties go to the run counts nearest the
reference's. It is a small dynamic program: a hull of at most 100 bases, and a few runs.

The parse gives the allele's **run counts** (§4) and each run's place inside the allele, which is
where a stutter slip is placed (§7.3). An allele with a substitution or a part-unit change still
parses — the mismatches are its substitutions, and a part-unit change shows as one; it simply costs
more mismatches than a clean allele.

This replaces, for bundles, the segment's `repeat_count_of_bases` (allele length ÷ period), which
has no meaning with two periods.

### 7.2 Candidates

The STR segment's selection (`candidate_alleles_ssr.md`), with each single-motif rule generalised:

- **The bar is unchanged**: a sequence is a candidate when at least two crossing reads of one sample
  show it, and at least a tenth of that sample's crossing reads.
- **A rung of the ladder is a run-count vector**: two alleles with the same run counts are two
  sequences on one rung, as two spellings of one repeat count are today.
- **The periodicity check** asks whether an allele's length change from the reference is a sum of
  whole units of the runs' periods — `Σ kᵣ·periodᵣ` for some whole numbers `kᵣ`, one per run `r` —
  instead of a multiple of one period. With no run, every length change is off the grid.
- **The ±1-unit rescue** offers one unit more and one unit less in each run separately.
- **The tie-break toward the purest allele** counts the parse's mismatches instead of bases off the
  one motif.

### 7.3 Read likelihood: stutter in every run

**A read is scored against an allele as the allele with any combination of slips, one per run at
most a few units each, then read with the substitution rate.** Each run slips independently, by `kᵣ`
units of its own motif, with the probability its own stutter model gives for its period and its
count in this allele. So the probability that a read shows the allele changed by `Δ` bases is the
sum, over every set of slips `(k₁, …, k_R)` whose lengths add up to `Δ` (`Σ kᵣ·periodᵣ = Δ`), of the
product of each run's probability of its `kᵣ` — including `kᵣ = 0`, the probability that run did not
slip.

- Only the slip sets whose total equals the read's length difference can spell it, so the sum is
  over few terms: for two runs of periods 2 and 3 and `Δ = 6`, the sets are `(3, 0)` and `(0, 2)`
  within ±4 units.
- Within one run, a slip is placed as in an interrupted segment today: at each distinct place in
  that run, weighted equally (`ssr_emission.rs`, the placements over runs). That machinery already
  splits one motif's tract into runs and fixed interruptions; the change is a motif per run.
- A read whose length difference no slip set reaches is scored by the part-repeat regime, as a read
  off the grid is today (`read_likelihoods.md` §4.2).
- The outlier term's reachable set is the union of the slip sets' totals.

**Stutter parameters are looked up per run**: read group, the run's period, its count in this allele
— the same key an STR segment uses, so no new stratum and no new parameter. `TractScoringFits` grows
from `read groups × candidates` to `read groups × candidates × runs`.

### 7.4 The prior

**The seed of an allele is the product of its runs' seeds**, each from the length spectrum of its own
stratum (period, reference count), at its own offset from the reference count
(`population_diversity.md` §4). For one run this is today's seed; with no run, every allele takes the
flat seed. Treating the runs as independent is the plain assumption with nothing measured against
it; a bundle is one haplotype and its runs' lengths may well be correlated.

### 7.5 The record

One multi-allelic record over the hull, widened over the flanks by edge deletions exactly as an STR
segment's record (`str_edge_deletion_alleles.md` §5.8), and gathered with the generic records the
same way (§6 there). The repeat annotation becomes per run:

- `RU` and `PERIOD` hold one value per reference run, in order (`RU=CT,AT,T;PERIOD=2,2,1`). Their
  header `Number` goes from `1` to `.`, which changes the header of every VCF, segments included. A
  bundle with no run writes neither.
- **Run counts per called allele go in a field of their own, `REPCNB`**, because `REPCN` is one
  integer per called allele: a per-sample string, one value per called allele, the run counts joined
  by `_` — `20_8_8,21_8_8` for a sample carrying `(20, 8, 8)` and `(21, 8, 8)`. Written on bundle
  records only; `REPCN` is left off them, and stays an integer on segment records.
- `--explain-loci` writes the bundle's locus type, `repeat_bundle`, which exists, and its outcome.

## 8. One sample and thousands, three reads and three hundred

- **One sample.** The candidate bar is the sample's own reads, as for segments; with no cohort the
  prior is the fitted spectrum, or the stated flat one where nothing was fitted. Nothing here needs
  a cohort.
- **Thousands of samples.** The cost a locus adds is the slip sets per read and candidate: a few
  terms per pair of runs, and a hull of at most 100 bases holds few runs. Memory per locus is the
  per-sample tallies, as for a segment.
- **Three reads a position.** Only crossing reads count (§6.1), and a hull of 100 bases leaves a
  150-base read about 20 starting positions from which it crosses the hull and 15 bases of flank
  each side, against about 100 for a 20-base segment. At 3× a long bundle will often have no
  crossing read in a sample, and is uncalled there rather than called from partial reads. Short
  bundles — the median is 35 bases — are close to segments.
- **Three hundred reads.** The candidate bar's tenth of the sample's crossing reads is what keeps
  stutter out of the candidates as depth grows, as for segments; the generic path's failure was
  exactly a bar with no stutter model behind it.

## 9. The benchmark: bundles in the stuttering simulator

**The stuttering benchmark placed no bundles** (`str_stutter_benchmark.md` §7), and without it
nothing here could be measured under PCR stutter or in a cohort. Built first (§10, step 1):

- **The reference.** A share of the shapes the simulator lays are bundles: two members (most) or
  three, gaps of 1 to 14 bases — never 0 between two members of one motif, which the typing would
  read as one repeat — and hulls of at most 100 bases. Gap bases are drawn so they extend neither
  neighbour's motif, as the spacer does today. Members are drawn to resemble HG002's bundles
  (`str_stutter_benchmark.md` §11.1): motif lengths in HG002's shares (half mononucleotides, three
  in ten dinucleotides), most members at or just above the calling floor, and 7 in 100 mono- and
  dinucleotide members long (16 to 30 copies). The median hull is 36 bases against HG002's 34, and
  about one bundle in eight holds a mono- or dinucleotide run of 16 copies or more, as in HG002
  (284 of 2,410).
- **Alleles.** Length changes drawn per member; derived alleles (a substitution inside a member, a
  SNP beside the hull, a deletion across an outer edge) as for segments.
- **Stutter.** One stutter site per member in each molecule, each at its member's own rate. The
  molecule cutter already applies several sites a molecule overlaps.
- **Truth.** One record per bundle window. The bundles' table and allele tables carry per-member
  motifs and counts — in files of their own (`bundles.bed`, `bundle_alleles.tsv`,
  `bundle_genotypes.tsv`) rather than a kind column on the segments', so that the segments' files
  stay byte-identical when no bundle is laid (`str_stutter_benchmark.md` §11.2).
- **Scoring.** A cell is one sample at one bundle, scored exactly as a segment's cell — the
  scorer compares sequences over the window and needs no member structure — and reported
  separately from segments.
- **Checks.** The driver's typing check also expects the bundles, as `ssr_bundle` regions. HipSTR
  is left off bundles.

**What it does not lay yet: bundles whose repeats abut, overlap or carry impurities** — the 734 of
HG002's 2,410 bundles of §5.1. Its bundles are pure members separated by gaps, the easy case for the
run finder. Extending it to those shapes, with stutter in the runs the caller would find, is step 5
(§10; owner, 2026-10-11).

## 10. Building and measuring it

**Built as the only behaviour, with no switch**, as the edge-deletion work was. Each step is
measured against main and against the step before; the STR segments' numbers must not move at any
step.

1. **The simulator's bundles** (§9). Checked by typing: every simulated bundle is typed as one
   bundle with the drawn members. Then two baselines on the new runs, at 1 and 50 samples, PCR and
   PCR-free, 5× and 30×: main (bundles uncalled) and the generic-path patch of the plan's §4.
   *Done 2026-10-11*: on 50 samples, PCR, 30×, the generic path gets 10,071 of 10,457 bundle cells
   with a variant right and calls a variant in 756 of 14,543 without one.
2. **Typing** (§5). On HG002's regions: 2,410 bundles, 63 new satellites, every segment unchanged;
   the VCF unchanged, as nothing calls a bundle yet. *Done 2026-10-11.*
3. **Reference runs, reads and psp** (§6). Counted on HG002: runs per bundle, and bundles with no
   run. Counted on the benchmark psps: crossing reads per bundle against the truth's carriers. Here
   the flat-gap delimiter is compared with a per-run one on the simulated bundles, by how often the
   spelling equals the read's true sequence.
4. **Calling and the record** (§7). Measured on both benchmarks: HG002's bundles at 30× and 50×
   against the generic path's 901 right, 178 wrong, 46 missed, 16 false calls at 30×; the
   stuttering benchmark's bundles against both baselines of step 1, false calls first. The
   cross-platform VCF checksums move here (the header's `Number`, §7.5), and are re-recorded only
   after showing this change alone moves them.
5. **The simulator's harder bundles** (§9): abutting, overlapping and impure repeats, measured as in
   step 4.

## 11. Decided (owner, 2026-10-11)

1. **The size limit is the user's**: `--max-str-bundle-len`, default 100 (§5.2).
2. **Stutter is modelled in every run** (§7.3), and **every allele has a count per run** (§7.1).
3. **Runs, not members, carry the stutter**: runs of pure whole copies of the members' motifs, found
   in the reference hull, non-overlapping, longest first, each at least the typing's copy floor
   (§6.1). The members only list the motifs to look for.
4. **Run counts are written as `REPCNB`** on bundle records (§7.5); `REPCN` keeps its type.
5. **Reads that end inside a bundle are left out at first** (§6.1), and measured at 3–5× once the
   rest works. They are where segments gain most at low depth, so if bundles lag segments at 5×
   this is the first thing to add.
6. **One name, *STR bundle*,** in the glossary, the specs and the code; the glossary's *cluster of
   repeats* is retired.
7. **The simulator is extended to abutting, overlapping and impure bundles** once the caller works
   (§9, §10 step 5).
