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
- **Run.** A stretch of pure, whole copies of one motif of period 1 to 6, at least the typing's copy
  floor for that motif's length (`--min-copies`) — where a slip of the polymerase can happen.
  *(Changed in review, 2026-10-11: any motif, not only the members'; §6.1.)* The
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
*Changed in review, the same day: the members are not consulted at all; every motif of period 1 to
6 is looked for (§6.1).*

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

*Changed during implementation (2026-10-11, step 3).* **The count is a scan's only.** A scan charges
a dropped bundle's whole hull to `FlankClamped`, where the hull starts, when the stretch asked about
holds that base — as it charges a segment. The catalog's tally has no contig-end counter, by design
(`repeat_catalog/tally.rs`: a row at a contig end is never stored, so such a counter would read zero
for every genome), and the rare bundle the catalog's route drops there goes uncounted, as it did
before.

**5.4 The tallies** (`TypedRegionCounts`) count capped bundles as satellites; `ssr_bundles` counts
what remains.

## 6. Reads and the psp file

### 6.1 The hull's reference runs, and measuring a read

**The reference runs are found once per bundle, when its locus is built**, from the hull's reference
bases alone:

1. For every period from 1 to 6, find every maximal stretch of the hull that repeats with that
   period, counted under its smallest period (`ATATAT` is a period-2 stretch, not a period-4 or
   period-6 one) and made of `A`, `C`, `G` and `T` only. A motif and its rotations are one motif
   (`CT` also matches `TCTC…`). Its copies are its length divided by the period, rounded down; a part
   copy at its end stays in the run, as a reference segment keeps its part copy.
2. Keep the stretches with at least the typing's copy floor for that period (`--min-copies`,
   default 6, 4, 4, 3, 3, 3 copies for periods 1 to 6). *(As built: the calling default has been
   8, 6, 6, 6, 5, 4 since 2026-08-10, `MinCopies::default`; the run takes whatever the ground was
   typed with.)*
3. Choose runs **longest first** (in bases; a tie goes to the earlier start). A run that overlaps
   one already chosen gives up the shared bases, cut back to whole copies, and is kept if it still
   reaches the copy floor. The chosen runs never overlap.

The two examples of §5.1 come out as `(CT)20`, `(AT)8`, `T8` — the `T` where the `AT` run and the
`T` run meet goes to the longer `AT` run — and as `A9` with the rest fixed sequence (its `AAAAAT`
units reach two clean copies, below the floor of three for period 6). A bundle can come out with no run at all
— a cluster of impure repeats, none with a stretch at the copy floor. It is still called, with no
stutter modelled anywhere in it; how many there are is counted on HG002 at step 3.

*Changed during implementation (2026-10-11): every motif, not the members' motifs.* Step 1 first
looked only for the members' motifs, and **that missed the `AT` run of this section's own first
example**. At chr1:8,007,443 the typing lists two members, `CT` over all 64 bases and `T` inside it,
and no `AT` member (§5.1). So the hull came out as `(CT)20` at 0 and `T9` at 55, and the 16-base
`AT` run between them was fixed sequence, with no stutter priced in it. Looking for every motif
of period 1 to 6 (owner, in review) gives `(CT)20`, `(AT)8`, `T8` as written, and the second example
still comes out as `A9` and fixed sequence. Both are pinned by tests
(`region_typing/reference_runs.rs`). On HG002's 2,410 bundles the change finds 3,892 runs instead
of 3,767, and bundles with no run fall from 141 to 130 (§10, step 3).

What the rule as written left open, and how it was settled in the code
(`find_reference_runs`):

- **A stretch is at least one whole copy and one more base**: a single isolated copy, with nothing
  repeating, is not a stretch. It matters only under a copy floor of 1.
- **A motif and its rotations** (`CA` and `AC`) are one motif: a stretch is found once, in the
  phase its first base gives it.
- **Ties after length and start** go to the shorter period. Two stretches alike in length, start and
  period have the same first bases, so the same motif, and nothing is left to break.
- **"Cut back to whole copies"**: a stretch that gives up bases at its end is cut back to whole
  copies counted from its start; one that gives them up at its start is re-read from its new first
  base — still whole copies of a rotation of its motif — and keeps its natural end's part copy.
  `(AT)8` in the first example is the latter: the `AT` stretch starts on the last `T` of `(CT)20`,
  gives it up, and is read as `AT` from the next base.
- **A stretch split in two** by a run chosen before it (possible only when that run was itself cut
  back) keeps each piece that reaches the floor.
- **A run's motif is phase-faithful**: the run's own first `period` bases, so it can be any rotation
  of the unit (`TC` for a run `TCTC…`).

The bundle generator is then the STR segment's generator run over the hull. **Three things differ:**

- **The delimiter** — the aligner that finds where in the read the repeat starts and ends — is the
  **flat-gap delimiter** (`alignment/ssr_best_path_flat_gap.rs`, `SsrFlatGapAligner`), which needs
  only the hull's boundaries, not a motif (`alignment/mod.rs`, `RepeatGeometry`). The default
  delimiter for segments prices whole-unit slips of one motif and cannot describe a hull of
  several. Whether a per-run version of it earns its keep is measured, not assumed (§10, step 3).
  *(Changed in review, 2026-10-11: the flat-gap delimiter is built for bundles with the segments'
  firm-anchor test, `SsrFlatGapAligner::with_anchor_rule`; §10, step 3.)*
- **The edge-escape check** (`reference_run_escapes_the_tract`) uses, at each edge of the hull, the
  motif of the run that reaches that edge, and none where fixed sequence reaches it.
- **Only reads that cross the whole hull are used.** An STR segment also uses a read that ends
  inside the repeat, as a lower bound on its length. In a bundle, such a read gives a lower bound on
  one run's count and says nothing about the runs it did not reach; using it needs a censored
  emission over run counts, which this design leaves out (§11).
  *(Changed 2026-10-12, owner: such reads are used, scored from the flank they held — §7.6.)*

Everything else is unchanged and motif-free: the spelling of the read over the hull, the flank
deletion at the two outer edges in its smallest form (`str_edge_deletion_alleles.md` §4), the
per-sample tally keyed on spelling, flank deletion and read group.

*Changed during implementation (2026-10-11, step 3).* As built (`locus_generation/ssr_bundle.rs`),
with these details the section left open:

- **The flanks are the segments'**: 15 bases each side by default (`SsrGeneratorConfig`'s
  `flank_bp`, held at or under the bundle radius as for a segment), clamped at a contig end.
- **A read that ends inside the hull is dropped and counted**, as `ended_inside_a_bundle`, a fifth
  no-observation reason; a read that never entered the hull is still `OutsideTract`. On the
  stuttering benchmark (one sample, PCR, 30×) a bundle cell has 20 crossing reads and 23 reads that
  yield no observation, against 20 and 6 at a segment.
- **The geometry's motif is the first member's.** The flat-gap delimiter reads none, but the
  geometry type carries one.
- **The runs are found after the hull's bases are fetched and before any read is measured**, so the
  edge-escape check has them; a bundle with no run at an edge asks nothing there.

*Changed during implementation (2026-10-11, step 4; session lead's decision).* **A read that ends
inside the hull is kept in the psp as a partial observation, as at a segment, and calling leaves it
out** (`call_one_ssr_bundle_locus` took each sample's complete observations only; since 2026-10-12 it
takes them all, §7.6). Dropping it
had kept it out of the hull's depth too, which the hidden-duplication filter's coverage windows
read (`window_coverage/depth.rs`): a bundle's depth was its crossing reads alone, about half the
reads there (step 3's note above). Kept as a partial, it reports depth at the positions it
witnessed, as a segment's partial does. The `ended_inside_a_bundle` no-observation count is gone
with the drop; such reads are counted as partial observations. A bundle record's `DP` now counts
them as a segment record's does. What it moved is in §10, step 4. *(Since 2026-10-12 calling
scores them too, §7.6.)*

*Found in review (2026-10-11, step 4).* **A deletion across a hull's outer edge is recorded as a
flank deletion even when the mapper wrote it some other way.** The simulator's edge deletions
at bundles sit on a length allele, so the read is often the reference's length or longer, and
the mapper writes the event as substitutions (`CTCCGA|AAAC…` read as `CTCAAC|AAAC…`) or as a
deletion inside the flank and an insertion at the hull's first base (`10M 4D 2M 6I …`). No
deletion then crosses the hull's edge, so the read's own alignment records nothing, and the
flat-gap delimiter — which prices nothing inside the hull and a stiff gap in the flanks — either
deleted the flank bases and spelled the read's extra bases into the hull, or could not anchor
that flank at all. The first described the event twice, as hull bases in the bundle record and
as the flank change the generic path calls beside it; the second made every carrier read a
partial, which calling leaves out. Together: 24 of the 51 wrong bundle cells at 50 samples, PCR,
30×, in three bundles. A segment's own delimiter prices an insertion inside the repeat that is
not whole units, so it keeps such a flank in place and leaves the change to the generic path,
once.

So at a bundle, **where the delimiter itself says a flank lost bases** — it anchored a complete
read with fewer read bases on that side than the flank holds, or could not anchor that side of a
read whose own alignment deletes bases in that flank — the read's bases over the window are
read as the flank's longest stretch matching from its far end, the deleted bases between it and
the hull, and the spelling, put in smallest form (`locus_generation/ssr.rs`,
`flank_deletion_the_delimiter_implies`). Both flanks must keep at least five matching bases
(`MIN_ALIGNED_FLANK_BASES`), the read's own alignment must cover the whole window and reach five
aligned bases beyond the hull on each side, and a side the delimiter anchored in full keeps its
junction. A flank the delimiter left intact never yields a deletion, so a substituted flank base
reads as it does at a segment. Segments are untouched. What it moved is in §10, step 4.

### 6.2 The psp file

**`LocusKind::SsrBundle` gets a payload**: the two flanks, as `SsrDetail` has, and the reference runs
as `(offset in the hull, length, motif)`. The members are not stored: once the runs are found, nothing
downstream reads them. Its tag (`KIND_SSR_BUNDLE = 2`) exists and no file has ever carried a bundle
record, so no reader has an expectation to break — but a reader built before this change would read
the tag and no payload, so the format version must make such a reader refuse a file with bundle
records (`psp_file_format.md` §4.3). The typing a psp records changes too, so new psps are needed for
any run that is compared.

*Changed during implementation (2026-10-11, step 3).* **Where the payload sits, and the version.**

- **The payload is one byte string at the very end of every record body**, `str-bundle-payload`,
  empty on every record that is not a bundle — not fields after the kind tag. That follows format
  1.2's flank deletions: a field every record carries is one a reader can step over without knowing
  it, which a field written only for one kind is not (`psp_file_format.md` §4.3, `RecordLayout`).
  It costs one byte a record before compression: on the stuttering benchmark's one-sample, 30×
  psp, the file grew from 35,116,701 to 35,309,906 bytes (0.55%), the 500 bundle records' own
  observations included.
- **Format 2.0, not 1.3.** A reader accepts every later minor of its own major, so only a major
  bump makes a reader built before this change refuse these files. Two things make refusal right.
  The payload is the body's last field, so such a reader would skip it as an unknown final field
  and read every bundle as one with no runs. And it ignores the header's `max-str-bundle-len-bp`
  (§5.2), so it would accept a typing it cannot reproduce.
- **A 1.x header passes a 2.0 reader's version check, but in practice every 1.x psp must be written
  again.** The typing section refuses a psp without `max-str-bundle-len-bp` (§5.2), and no
  released build wrote that key: only this branch's builds from 270da29f to format 2.0 did. A 1.x
  file that does pass holds no bundle record, and one found in it is refused as damage.
  *(Corrected in review, 2026-10-11: this said a 2.0 reader "still reads 1.x files".)*
- `pvc inspect-psp` writes a bundle's runs in its `motif` cell as `motif:offset:length`, joined by
  `,`, and the hull's flanks in the flank cells.

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

*Changed during implementation (2026-10-11, step 4).* As built (`region_typing/run_parse.rs`,
`parse_into_runs`), with what the section left open:

- **An edit is a substituted, an inserted or a deleted base, one each**; "fewest mismatching
  bases" is the fewest edits. The parse is an edit-distance alignment of the allele against the
  template `fixed₀ (motif₁)* part₁ fixed₁ … fixed_R`, each run allowed any whole number of copies,
  zero included. Ties go to the smallest `Σ |countᵣ − reference countᵣ|`; what ties remain go to
  the reading the alignment meets first, which is fixed.
- **A run's place** is the stretch of the allele its copies and part copy were read from, **and any
  base inserted where the run meets fixed sequence** — so a run spanning the whole hull has the
  whole allele as its place, which is what lets a one-run bundle score as a segment (§7.3). Where
  two runs meet with no fixed base between, such a base goes to the earlier run.
- **A bundle with no run** parses every allele as fixed sequence: no counts, and its edits are the
  allele's edit distance from the reference hull.
- The flanks are not parsed: an allele's bases are the hull's spelling, and its flank deletion
  travels beside them (§6.1).

### 7.2 Candidates

The STR segment's selection (`candidate_alleles_ssr.md`), with each single-motif rule generalised:

- **The bar is unchanged**: a sequence is a candidate when at least two crossing reads of one sample
  show it, and at least a tenth of that sample's crossing reads.
- **A rung of the ladder is a run-count vector**: two alleles with the same run counts are two
  sequences on one rung, as two spellings of one repeat count are today.
- **The periodicity check** asks whether an allele's length change from the reference is a sum of
  whole units of the runs' periods — `Σ kᵣ·periodᵣ` for some whole numbers `kᵣ`, one per run `r` —
  instead of a multiple of one period. With no run, every length change is off the grid.
  *(Removed 2026-10-12, owner: a bundle has no periodicity check — see the end of this section.)*
- **The ±1-unit rescue** offers one unit more and one unit less in each run separately.
- **The tie-break toward the purest allele** counts the parse's mismatches instead of bases off the
  one motif.

*Changed during implementation (2026-10-11, step 4).* As built (`select_ssr_bundle`, sharing the
segment selector's passes through `build_ladder_keyed`, `promote_rungs_for_cohort_by` and
`admit_promoted_sequences_on`):

- **The periodicity check is the greatest common divisor of the runs' periods**: a length change
  is `Σ kᵣ·periodᵣ` for some whole numbers `kᵣ` exactly when that divisor divides it. With no run
  the divisor is zero and only no change is on the grid, so a no-run bundle whose reads change
  length in more than a tenth of a sample's reads (and two reads) is called over the reference
  alone, as `notPeriodic`. A deletion across the hull's edge never counts as off the grid, as at a
  segment. *(Removed 2026-10-12, below.)*
- **Rungs** are the distinct run-count vectors of the merge's alleles, numbered in ascending order;
  the rescue's neighbours are the vectors one copy either side in one run, among those.
- **Every bundle allele can be scored** (a run it holds no copy of cannot slip), so a bundle has no
  counterpart of a segment's "candidate with no whole repeat" and nothing stops its locus.
- **The discovery pre-pass is not run on bundles.** It is off by default for segments, and its bar
  is written in a segment's repeat counts.
- `notPeriodic` and `tooManyAlleles` reach a bundle record's `FILTER` as a segment's. *(Since
  2026-10-12 only `tooManyAlleles` can, below.)*

*Changed 2026-10-12 (owner): **a bundle has no periodicity check.** The fixed sequence between and
around a bundle's runs carries real indels of any length, so a length change that is not a sum of
whole run units is evidence of an allele, not a sign that the locus is not a repeat. At step 4 the
check refused 31 of HG002's 149 missed bundles, among them a 9-base insertion in a bundle with no
run and a one-base insertion in fixed sequence, and refused nothing on the simulator. Selection
at a bundle now goes straight to the rungs; `notPeriodic` is never a bundle's verdict, and
`tooManyAlleles` still is. Segments keep the check unchanged. What it moved is in §10, step 4b.*

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

*Changed during implementation (2026-10-11, step 4).* As built (`ssr_emission.rs`,
`bundle_emission`; `repeat_tract_parameters.rs`, `gather_for_bundle`):

- **Each run's slip reaches the model's own cutoff**, `MAX_WHOLE_REPEAT_SLIP` (10) units either way,
  and must leave the run a whole copy; a run the allele holds no copy of does not slip. The
  example above, `(3, 0)` and `(0, 2)`, is incomplete even within ±4: `(−3, 4)` also totals 6.
- **Where several runs slip at once**, every combination of their places is one rendering of the
  allele, weighted by the product of each run's equal weights.
- **The part-repeat regime** prices a change no slip set reaches as a sequencing indel in one run:
  that run's model for a change that is not whole units of its period, resized at the run's end,
  times every other run's chance of not slipping, averaged over the runs for which the change is
  not whole units. A run for which it is whole units, past the cutoff, gives nothing, as at a
  segment. With no run, any length change scores zero and goes to the outlier term.
- **The outlier term's reachable set** is the slip sets' totals and every run's part-repeat changes
  — for one run, exactly the segment's set.
- **From three runs on, a branch of slip sets whose probability cannot reach 10⁻¹² is not walked**
  (`NEGLIGIBLE_SLIP_SET_PROBABILITY`); four runs would otherwise walk thousands of sets a read and
  five a hundred thousand. Two runs are walked whole.
- **One substitution rate per cell, the mean of the runs' fitted rates** (each looked up as the
  run's stutter is), because the letters are compared over the whole allele at once; a bundle with
  no run takes the default rate and its cells say so. A cell's fallback counts — slippage,
  substitution rate — count a cell once where any run of it fell back.
- **The per-locus slippage re-fit does not run at a bundle** (it is off by default for segments).
- **No censored read**: calling hands a bundle only reads that crossed the hull (§6.1).
  *(Changed 2026-10-12: it hands them all, §7.6.)*
- *Found in review (2026-10-11).* **A run the allele holds no copy of is a factor of one** in
  every slip set and in the part-repeat regime, not its model's no-slip share: that share (0.88
  under the shipped model, which such a run is given) scaled down every read of an allele that
  lost a whole run, 0.02 over 30 reads. **A run with no copy, or one whose stretch is shorter
  than the contraction, is not a part-repeat place**: it rendered an allele of another length
  than the read, a panic in a debug build and in release a place scoring zero that still
  counted in the average. **Where every slip set reaching a change was pruned**, the read scores
  zero rather than falling to the part-repeat regime, which is for changes no slip set reaches.
  None of the three moved a benchmark cell (§10, step 4).
- **The regression this had to pass**: a bundle of one run spanning its hull scores every read
  exactly as the segment over the same bases, bit for bit — checked over four hulls, every allele a
  few copies either side, interrupted and part-unit alleles, and reads from twelve copies short to
  twelve long (`a_one_run_bundle_scores_as_the_segment_it_is`).

*Changed 2026-10-12 (owner; the review's finding M3).* **A length change that no single run's slip
makes is priced as a sequencing indel as well as by the slip sets that reach it.** Before, the
part-repeat regime priced only a change no slip set reached at all. At a bundle whose runs' periods
share no divisor above one — periods 2 and 3, say — every change is some slip set's total, so a
one-base indel was priced only as two runs slipping opposite ways: about 0.0024 under the shipped
constants, against 0.0095 for the same indel at a dinucleotide segment. Now, for a change `Δ` that
is not whole units of any run within that run's own reach, the emission is the slip sets' sum
**plus** the part-repeat regime's price, the one described above. A change one run makes alone is
that run's slip and takes no such term, as at a segment: one base more at a bundle with a
homopolymer run is the homopolymer slipping. With one run nothing changes, since every change that
run's slip does not make is one no slip set makes; the one-run regression above still holds bit for
bit. Where every slip set reaching a change was pruned (three runs or more), a change one run makes
alone scores zero as before, and one no single run makes takes the part-repeat price. Pinned by
`a_change_no_single_run_makes_is_also_priced_as_a_sequencing_indel` and
`a_change_one_run_makes_alone_takes_no_indel_term`.

### 7.4 The prior

**The seed of an allele is the product of its runs' seeds**, each from the length spectrum of its own
stratum (period, reference count), at its own offset from the reference count
(`population_diversity.md` §4). For one run this is today's seed; with no run, every allele takes the
flat seed. Treating the runs as independent is the plain assumption with nothing measured against
it; a bundle is one haplotype and its runs' lengths may well be correlated.

*Changed during implementation (2026-10-11, step 4)* (`seed_ssr.rs`, `fill_ssr_bundle_seed`). **The
product is of the runs' shapes; the concentration is the smallest of the runs'** — multiplying the
seeds outright would hold a two-run bundle with the square of a run's conviction. **A run on the
stated flat rung contributes no factor**; where no run has a fitted spectrum, or there is no run,
the shape is flat over the candidates at the stated concentration, as a segment's flat rung is. The
record's length-spectrum rung is the weakest of the runs'. For one run this is today's seed, bit for
bit (`a_one_run_bundles_seed_is_its_runs_segment_seed`). Each run's offsets are measured from its
own count in the reference allele's parse.

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

*Changed during implementation (2026-10-11, step 4).* **A bundle with no run writes no `STR`, `RU`,
`PERIOD` or `REPCNB`** (`vcf_output.md` §6: the first three travel together); its record is still
the repeat path's for the gathering and the writer's ordering. Each allele's run counts are read off
the alleles the record keeps, before any widening. The explain outcome `bundle_set_aside` is gone;
the run report counts bundles among the tract outcomes and says how many were bundles. **The
cross-platform VCF checksums moved with the header alone**: with the four changed header lines put
back (`STR`'s description, `RU` and `PERIOD`'s `Number`, the new `REPCNB`), each new file hashes to
its old checksum; they were re-recorded (`cli/cross_platform_digests.rs`).

### 7.6 Reads that end inside a bundle

*Added 2026-10-12. The owner's decision of 2026-10-11 (§11, fifth item) was to add these reads
first if bundles lagged at low depth, and they did: at one sample, PCR, 5×, step 4 got 98 of 208
bundle cells with a variant right against the generic path's 115.*

**What such a read says.** A read that holds one of the hull's flanks and ends inside the hull
shows the start of what its carrier has from that flank on. Coming from the left flank, it gets
past some runs and fixed bases whole and then stops — inside a run, inside fixed sequence, or a few
bases past the hull's far edge, too few for the delimiter to anchor the far flank. So:

- **about every run and fixed base it got past, it says what a complete read says** — their
  spelling, with the stutter of those runs priced exactly as at a complete read;
- **about the run it ended inside, it gives a lower bound**, as a segment's censored read does on
  its one run: every length of that run long enough to hold what the read showed agrees with it;
- **about the runs it never reached, it says nothing.**

A read from the right flank is the same read from the other end.

*Changed 2026-10-12 (session lead, the review's findings M2 and M3; the rule as first built is
kept below it).* **The rule, as it stands.** For a read that held the left flank and showed bases
`w`, the carrier is read piece by piece from that flank: the fixed piece before the first run,
the first run, the fixed piece after it, and so on, then the right flank. **Each run makes one
change, independently of the others**: a slip of `k` whole units at its stutter model's price
`P_r(k)`, placed at each place in the run as at a complete read, or a sequencing indel of `Δ`
bases that are not whole units, at that model's price for `Δ` (`|Δ| ≤ 2 · MAX_PART_REPEAT_SLIP`,
the stretch long enough to lose it), resized at the run's end — the places and prices the
part-repeat regime of §7.3 uses in that run.

```text
P(w | allele) = Σ over every way the runs the read reaches can change   Π over those runs   price of its change
                                                                       × letters(w | the first |w| bases of that carrier)
```

**The sum stops where the read's bases run out.** On every path through the runs' changes, the
pieces after the point where `|w|` bases have been laid down are not compared, and contribute
exactly one: a run the read never reached says nothing, whatever count the allele holds there.
So the three readings above are now exact:

- about the runs and fixed bases it got past, it is evidence as at a complete read;
- about the run it ended inside, a lower bound: every change of that run long enough to hold
  what it showed gives the same letters, and only that run's own total price over those changes
  separates its counts;
- about the runs it never reached, nothing: **two alleles that differ only beyond everything the
  read can reach score it bit for bit alike** (`a_partial_read_scores_alike_against_alleles_that_differ_beyond_its_reach`).

A read that held the right flank is walked from its last base backwards, the left flank last.

**How this differs from a complete read, and why.** For one run it is the segment's censored sum,
term for term (`a_one_run_bundles_partial_read_scores_as_the_segments`, to twelve digits). With
several runs, a complete read prices a sequencing indel only for a length change no single run's
slip makes, averaged over the runs it could sit in, with the other runs not slipping (§7.3). That
rule reads every run of the allele — which changes count as made by one run, how many places
share the average, each other run's chance of not slipping — so a read that stopped short of a
run cannot follow it and stay neutral about that run. Here each run the read reaches offers its
own sequencing indels beside its slips, at its own model's price. **The two rules do not price
the same renderings** *(corrected 2026-10-12, follow-up review)*: a partial read lets every run it
reaches take a sequencing indel, whatever the other runs do, and adds them up; a complete read
takes one only where no single run's slip reaches the length, averages it over the runs and has
the other runs not slip. So some renderings the complete rule gives nothing get a price here — at
`(CA)8 GT (TTG)5`, the `CA` run three bases longer, a length the `TTG` run's slip reaches alone.
Each run's changes are still a distribution (they sum to at most one, and its
indels skip every whole-unit length, so no length is counted twice), and the partial read's sum
holds every term of the complete read's: a partial read is at least as likely as the complete read
of the same bases, as at a segment
(`a_two_run_partial_read_is_at_least_as_likely_as_the_complete_read`). What the change moved on
the benchmark is in §10, step 4c.

**The rate the letters are compared at** *(2026-10-12, follow-up review)*. The substitution rate is
fitted per repeat count, and a cell's one rate — what a complete read is compared at — is the mean
of its runs' (§7.3). Compared at that mean, a partial read would hear about a run it never reached
through the far run's count. So a partial read compares **each run's bases at that run's own
rate**, and **every other base** — the fixed pieces and the flank beyond the hull — at one rate per
read group that no candidate's counts move: the mean of the runs' rates at their reference counts
(`a_partial_read_is_compared_at_the_rates_of_what_it_reached`). Where no rate is fitted, every one
of these is the defaulted rate, and the read is compared exactly as before.

**What made the read depend on the runs it never reached (finding M2).** Two things, both in the
rule as first built. Each run past the read's end contributed its total probability over the
slips it can make, and that total depends on its count: a run of one copy cannot contract, and a
run with no copy is a factor of one. And the part-repeat regime was taken whole from the complete
read: which length changes it priced, how many runs shared its average and each other run's
chance of not slipping all read the runs the read never reached. Under the shipped constants the
read `GG (CA)8 G` stopping before the `TTG` run scored 0.863, 0.811, 0.852 and 0.854 against
alleles with 0, 1, 2 and 5 `TTG` copies; it now scores 0.960567 against those four and against
nine, equal to thirteen digits
(`a_read_stopping_before_a_run_scores_alleles_differing_there_alike`; the carriers whose `CA` run
contracted reach the `TTG` run, and they mismatch the read at two bases or more).

**The rule as first built (2026-10-12, before the review)**, kept for the record. For a read that
held the left flank and showed bases `w`:

```text
P(w | allele) = Σ over every slip set (k₁, …, k_R)   Π_r P_r(kᵣ) × letters(w | the first |w| bases of: the allele with those slips, then the right flank)
              + Σ over every length change Δ no single run's slip makes   its part-repeat price (§7.3) × letters(w | the first |w| bases of: that rendering, then the right flank)
```

— over every slip set, not only those of one total, because the read does not say how long its
carrier's hull is. This is the complete read's emission of §7.3 summed over every complete read
whose first bases are `w`, which is what a segment's censored term is (`read_likelihoods.md`
§5.2). The three readings above are properties of the sum, not rules added to it:

- a run past the read's end changes no base compared, so its slips sum to its total probability
  over them, whatever count the allele holds there;
- every rendering of the run the read ended inside that is long enough to hold what it showed
  gives the same letters;
- a rendering too short for what the read showed is compared, past its end, with the bases that
  follow it in the carrier — the next fixed piece, or the flank — so a read showing more of a run
  than an allele can reach scores low against it.

A read that held the right flank is compared from its last base backwards: against the last `|w|`
bases of the left flank followed by the rendering.

As built (`ssr_emission.rs`, `bundle_censored_emission`), with what the rule leaves open:

- **Which flank the read held is read off its witness** (`ReadWitness`): one that starts at the
  hull's first position and stops short of its last held the left flank, and the mirror image the
  right. A witness over the whole hull is a read that reached the far edge without anchoring
  there, from either flank; nothing in the observation says which, so it is scored both ways and
  the two averaged. *(Changed 2026-10-12, finding M3: the classifier records the flank, below.)*
- **The far flank is the reference's, less the bases the allele's own edge deletion takes on that
  side** (`str_edge_deletion_alleles.md`). The held flank is not compared: the delimiter anchored
  it. A read that showed more bases than the rendering and the flank hold together scores zero for
  that rendering.
- **One run over the hull, with no flank bases to read on into, scores as the segment over the same
  bases**, to twelve digits (`a_one_run_bundles_partial_read_scores_as_the_segments`; not bit for
  bit, since the segment takes the letters out of a pure tract's sum). Two runs are checked against
  the sum written out by hand, from each flank and from both
  (`a_partial_read_sums_every_rendering_from_the_flank_it_held`).
- **It is computed run by run from the held flank**, not set by set: what follows a run depends
  only on how many of the read's bases the pieces before it took, so that is remembered, and a read
  costs a few slips per run per such offset rather than every combination of slips.
  *(Follow-up review, 2026-10-12.)* The walk now goes forward: each state — a piece reached with
  so many of the read's bases spent — carries the probability of reaching it and hands it on. **At
  a bundle of two runs or more, a state no more probable than what is left of a budget of 10⁻¹² is
  dropped** (`NEGLIGIBLE_PARTIAL_READ_MASS`): everything after a state is at most one, so the read
  comes out at most 10⁻¹² low and never high — against an outlier floor near `0.2 / 40`, the bound
  a complete read's slip sets are pruned to (§7.3). What it drops are mostly the states a
  sequencing indel shifts the read into, a few mismatches deep long before the read's bases are
  spent. A bundle of one run is walked whole and stays the segment's sum to twelve digits.
  *(Review, 2026-10-12.)* What does not depend on the read is built once and kept: an allele's
  layout from one flank — its pieces, its runs' slips and renderings, its part-repeat renderings —
  is reused for every read the row scores against that allele from that flank, found by comparing
  everything it was built from; the two per-base scores are looked up once a read; a fixed piece
  reached at one offset by two paths is compared once, as a run is; and a comparison's sum of match
  scores up to its first mismatch is taken from a table. Each gives the same numbers bit for bit;
  the calls are byte-identical (§10, step 4b).
- **A bundle with no run** compares the read's letters with the allele and the flank alone.

*Found in review (2026-10-12)* (`doc/devel/reports/reviews/str_bundles_partial_reads_2026-10-12.md`,
M2 and M3). *Both fixed the same day (session lead): the first by the rule at the top of this
section, the second by the paragraph after this list.*

- **A read says a little about the runs it never reached.** Each such run contributes its total
  probability over the slips it can make, and that total depends on its count: a run with one copy
  cannot contract, and a run with no copy is a factor of one. Under the shipped constants a read
  `GG (CA)8 G` that stops before the `TTG` run scores 0.863 against an allele with no `TTG` copy,
  0.811 with one, 0.852 with two and 0.854 with five or more — about 6 in 100 between one copy and
  none, in the same direction for every such read. A complete read's emission loses the same mass;
  the claim above that the read "says nothing" about those runs is true only up to it. The choice
  is a factor of one for every unreached run, or this paragraph as the record of the limitation.
- **A read whose witness covers the whole hull is scored at half weight.** The read classifier
  knows which flank it anchored, but the witness keeps only the positions, so §7.6 averages both
  readings; the wrong one is about 4 parts in 10¹³ of the right one, so the read scores half of
  what the same bases score as a left-only read, up to ln 2 = 0.69 nats against the outlier floor.
  Such reads were 185,000 of 1.6 million partial emissions at the 50-sample cost benchmark.

*Changed 2026-10-12 (session lead, finding M3).* **The flank such a read held is carried to
scoring.** The classifier lays a partial read's reach down from the border it anchored and clamps
it to the locus, so a read that ran past the far edge covers every position and its witness no
longer says its side; the classifier still knows it. At a bundle, and only where the witness
covers the whole hull, it now writes the side down on the observation
(`SequenceObservation::held_flank`, part of the observation's identity, so the same bases read
from the two sides stay two observations). The cohort merge carries it on the partial row, and
the emission reads it: a whole-hull read is scored from the flank it held. The average of both
readings is kept only where the side is genuinely unknown — a psp written before format 2.1, or a
witness touching neither border, which no generator mints. Segments record nothing, since they
read every partial from the tract's start (below), and their observations are unchanged.

In the psp the side is format 2.1's new last body field, `observation-held-flanks`: a list of
`(observation, flank)` pairs, empty — one byte — on nearly every record. A minor bump, as format
1.2's flank deletions were (`psp_file_format.md` §4.3): a 2.0 reader walks past it as an unknown
final field and scores such a read as 2.0 did, from both flanks one half each — a proper mixture
over a side it does not know, not a misreading. `pvc inspect-psp` writes it in a new last column,
`held_flank` (columns version 3). On the benchmark's new psps every partial over a whole hull
carries its flank, and none falls back to the average: 184 such reads at one sample 5×, 1,030 at
30×, 50,776 at 50 samples 30× — about one in eight of the bundles' partial reads, split evenly
between the two flanks (§10, step 4c).

**Candidate selection counts these reads nowhere, as at a segment.** A segment's partial reads
nominate no rung and count toward no denominator (`candidate_alleles_ssr.md` §8: "a read that ran
out inside a tract gives a lower bound on the length, not a length"); the merge holds them off the
support rows that selection reads, at bundles exactly as at segments. So the candidates are still
what two crossing reads of one sample, and a tenth of them, show, and a partial read only weighs
between those candidates. At low depth that leaves a bundle whose other allele no two crossing
reads show with no candidate for it, whatever its partial reads show.

**One difference from a segment, left as it is.** A segment reads every partial's bases as the
start of its tract, whichever flank the read held (`censored_emission`). For a read that held the
right flank of a pure tract that is the same comparison when the read showed whole copies and a
phase-shifted one otherwise; segments are not changed here. *(Measured in review, 2026-10-12, on
one simulated sample, PCR, 30×: 14,175 of the segments' 22,661 right-flank partial reads are
compared out of phase — none at period 1, half at period 2, three in four at periods 4 to 6 — 27 in
100 of all 51,818 partial observations. Such a read mismatches nearly every base against every
candidate, so it falls to the outlier floor and its evidence is lost. A separate change.)*

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
  bundles — the median is 35 bases — are close to segments. *(Since 2026-10-12 partial reads are
  scored (§7.6), but they still nominate no candidate, so the bundle is still uncalled where no two
  crossing reads show its other allele.)*
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
   *Done 2026-10-11.* **HG002** (the 2,410 bundles, calling floors; `examples/ng_str_bundle_runs.rs`),
   runs of every motif (§6.1): 130 bundles have no run, 869 one, 1,234 two, 155 three, 20 four, 2
   five; the 3,892 runs are 2,384 of period 1, 1,268 of period 2, 83, 98, 35 and 24 of periods 3 to
   6. With the members' motifs only, as first built: 141, 922, 1,210, 124, 12 and 1 bundles, 3,767
   runs.
   **The stuttering benchmark** (PCR, 30×; each crossing read's spelling put back into the truth's
   window and compared with the sample's two true alleles; cells with a SNP beside the repeat left
   out). Per 1,000 crossing reads:

   | | spelled as a true allele | one whole-unit indel from one (stutter's shape) | same length, another base | a part of the hull only | another length |
   |---|---|---|---|---|---|
   | bundles, 1 sample (10,859 reads) | 921 | 23 | 45 | 7.4 | 3.2 |
   | bundles, 1 sample, with the anchor test (9,910) | 939 | 16 | 44 | 0.2 | 0.7 |
   | segments, 1 sample (59,653) | 932 | 28 | 39 | 0.0 | 0.8 |
   | bundles, 50 samples (542,911) | 906 | 38 | 45 | 7.7 | 3.6 |
   | bundles, 50 samples, with the anchor test (497,333) | 922 | 32 | 44 | 0.4 | 1.7 |
   | segments, 50 samples (2,970,633) | 902 | 59 | 38 | 0.2 | 1.6 |

   The flat-gap delimiter spells bundles as often right as the default delimiter spells segments,
   and stutter-shaped reads are no more frequent, so **no per-run delimiter was built**. One
   failure was its own: about 1 crossing read in 100 at a bundle was spelled over **part of the
   hull only** (or part plus an error base) — a read that ends inside the hull, reported complete
   because the flat-gap delimiter calls a flank anchored from the alignment's shape, not from
   matched flank bases. A read at chr1:144,818 shows how: it carries `A11` against the hull's
   `A10 CTACTC A8` and ends there; the path takes ten `A`s, deletes the rest of the hull, puts
   the eleventh `A` against the right flank's first base and runs off the frame's end, and the
   flank counts as anchored on one mismatched base. **Fixed in review (owner, 2026-10-11)**: the
   bundle generator builds the flat-gap delimiter with the segments' firm-anchor test
   (`SsrFlatGapAligner::with_anchor_rule`, the rule of `alignment/ssr_unit_robust.rs`), so such a
   read is one-sided and dropped as ending inside the hull. The rows "with the anchor test" above
   are that build: part-of-hull reads fall to 0.2 and 0.4 in 1,000, as low as at segments, and the
   segments' rows do not move. A bundle cell now has 20 crossing reads and 23 that yield no
   observation (a segment cell 20 and 6).
   **Nothing else moved** on the one-sample 30× run: every non-bundle record of the psp is
   byte-identical with main's build (4,120,085 dump rows), and the VCF's genotypes, QUALs and
   filters are identical; the hidden-duplication filter's `PARALOG_LR` moved in 67 of 2,396 records
   (66 with the anchor test) and its cut from 31.75 to 31.85, because its coverage now sees reads
   at bundle hulls, which `window_coverage/depth.rs` already counted for any bundle record. That is
   correct counting, not double counting: a hull is one typed region, so no other record reports
   its positions, which before had no record at all and were absent from the windows. Its depth
   there is the crossing reads only, about half the reads present (the review's Mi10). The
   cross-platform checksums did not move.
4. **Calling and the record** (§7). Measured on both benchmarks: HG002's bundles at 30× and 50×
   against the generic path's 901 right, 178 wrong, 46 missed, 16 false calls at 30×; the
   stuttering benchmark's bundles against both baselines of step 1, false calls first. The
   cross-platform VCF checksums move here (the header's `Number`, §7.5), and are re-recorded only
   after showing this change alone moves them.
   *Built 2026-10-11; first measurements, before review* (true stutter, no fit; HG002 by
   `tmp/region_confusion.py`'s per-region verdicts over main's typing of 2,473 bundles):

   | | bundle cells with a variant: right / wrong / missed | false calls, of cells with no variant |
   |---|---|---|
   | 1 sample, PCR, 30×: this build | 205 / 3 / 0 of 208 | 0 of 292 |
   | same: generic path / GATK | 206 / 1 / 1 — 207 / 1 / 0 | 0 — 0 |
   | 1 sample, PCR, 5×: this build | 98 / 30 / 80 of 208 | 0 of 292 |
   | same: generic path / GATK | 115 / 34 / 59 — 121 / 32 / 55 | 0 — 0 |
   | 50 samples, PCR, 30×: this build | 10,349 / 51 / 35 of 10,435 | 3 of 14,565 |
   | same: generic path / GATK | 10,071 / 229 / 157 — 10,097 / 290 / 49 | 756 — 734 |
   | HG002, 30×: this build | 824 / 158 / 149 of 1,131 | 32 of 1,342 |
   | same: generic path | 901 / 184 / 46 | 16 |

   **Every STR segment cell is unchanged** against main on all three benchmark runs (3,000, 3,000
   and 150,000 cells, outcome and GQ). At 50 samples the false calls fall from the generic path's
   756 to 3. At 5× bundles lag, as §11's fifth decision foresaw: 80 cells have no call, against 59
   for the generic path, because a cell with no crossing read has none. **On HG002 this build is
   below the generic path**: of its 149 misses, 42 are hulls over 100 bases, now satellites; 31 are
   refused as `notPeriodic` — real alleles such as a 9-base insertion in a bundle with no run and a
   one-base insertion in fixed sequence; the rest have too few crossing reads to call. With the
   periodicity gate switched off at bundles (an experiment, not built), HG002 moves to 842 / 169 /
   120 with the same 32 false calls, and the three benchmark runs do not move. Of the 32 false
   calls, most are low-quality heterozygotes over an interior substitution of an impure
   homopolymer (`TTTTTTTTCTTT…` against two reads of pure `T`, QUAL 3). On the benchmark, 24 of
   the 51 wrong cells at 50 samples carry a deletion across the hull's edge, which the flat-gap
   delimiter spells as bases inside the hull and the generic path calls beside it, so two records
   describe one event.
   **Keeping the reads that end inside the hull** (§6.1) moved no genotype, filter or `AD` on the
   one-sample 30× run (step 3's psps against step 4's, both called by step 4); it raised `DP` at
   204 bundle records and `PARALOG_LR` at 249 records (204 bundles and 45 records whose windows
   overlap one), every one upward, median 0.63, at most 3.77; the filter's cut moved 32.25 →
   32.15, and nothing was tagged or dropped either way.
   *After review (2026-10-11)*: the review's fixes — the edge deletion recorded as a flank
   deletion (§6.1) and the three emission corrections (§7.3) — measured on new psps with the same
   tools. **Every STR segment cell is again unchanged against main** (3,000, 3,000 and 150,000
   cells).

   | | bundle cells with a variant: right / wrong / missed | false calls |
   |---|---|---|
   | 1 sample, PCR, 30×: before review → after | 205 / 3 / 0 → 207 / 1 / 0 of 208 | 0 → 0 of 292 |
   | 1 sample, PCR, 5× | 98 / 30 / 80 → 98 / 30 / 80 | 0 → 0 |
   | 50 samples, PCR, 30× | 10,349 / 51 / 35 → 10,374 / 26 / 35 of 10,435 | 3 → 3 of 14,565 |
   | HG002, 30× | 824 / 158 / 149 → 822 / 160 / 149 of 1,131 | 32 → 32 of 1,342 |

   All of the movement is the edge-deletion fix; the emission corrections moved no cell. At 50
   samples the bundle cells carrying a deletion across the edge went from 27 not right (24 wrong,
   3 with the right lengths only) to 2 wrong, of 374. On HG002 two regions went from right to
   wrong. At chr4:28,575,007 the bundle record now spells the truth's flank SNP `G>C` beside the
   hull as one allele with the deletion inside it, widened over that base, where before the
   generic path called the SNP; the sequence is the truth's, and the region scorer, which reads
   the record against the region's bounds, does not credit it. At chr3:81,721,205 the haplotype
   with a one-base deletion just past the hull had two crossing reads; one now carries the flank
   deletion and the other does not, so each key has one read and neither clears the two-read
   bar — the split the smallest form exists to prevent, here between a read the rule reaches and
   one it does not.

   **Cost.** With one thread, 50 samples at 30×, the 500 bundle loci of the benchmark add 1.45 s
   of calling against the build that sets them aside (1.98 s against 0.53 s over the same
   psps), about 2.9 ms a bundle locus; 2,000 segment regions take 4.0 s, generic loci in their
   flanks included, so at most 2 ms a segment locus.

   **Step 4b (2026-10-12): the three owner's changes** — no periodicity check at bundles (§7.2),
   reads that end inside the hull scored (§7.6), and the indel term for a change no single run
   makes (§7.3). Measured as step 4 was (true stutter, no fit; new psps; HG002 per region over
   main's typing), each change added in turn on one build with switches, which were then removed;
   the final build reproduces the last row's VCFs byte for byte. Bundle cells with a variant:
   right / wrong / missed; then false calls.

   | | 1 sample, PCR, 5× (208; 292) | 1 sample, PCR, 30× (208; 292) | 50 samples, PCR, 30× (10,435; 14,565) | HG002 30× (1,131; 1,342) |
   |---|---|---|---|---|
   | step 4, reviewed | 98 / 30 / 80; 0 | 207 / 1 / 0; 0 | 10,374 / 26 / 35; 3 | 822 / 160 / 149; 32 |
   | no periodicity check | 98 / 30 / 80; 0 | 207 / 1 / 0; 0 | 10,374 / 26 / 35; 3 | 840 / 171 / 120; 32 |
   | + partial reads | 99 / 29 / 80; 0 | 207 / 1 / 0; 0 | 10,400 / 13 / 22; 1 | 862 / 150 / 119; 20 |
   | + indel term | 99 / 29 / 80; 0 | 207 / 1 / 0; 0 | 10,400 / 13 / 22; 1 | 862 / 150 / 119; 20 |
   | generic path | 115 / 34 / 59; 0 | 206 / 1 / 1; 0 | 10,071 / 229 / 157; 756 | 901 / 184 / 46; 16 |

   **Every STR segment cell is unchanged** against main on all three benchmark runs. The indel
   term moved no cell and no region; it changed the qualities of 7 records at 50 samples and 6 on
   HG002. Partial reads leave the one-sample 5× misses where they were — 80 cells with no call —
   because they nominate no candidate (§7.6): a cell whose other allele no two crossing reads show
   still has no candidate for it.

   **Cost of the partial reads.** With one thread, 50 samples at 30×, calling the 500 bundle loci
   took 1.9 s without them and 14–17 s with them (two runs each), about 25 ms more a bundle locus:
   each partial observation is laid out again for every candidate and walked over every slip of
   the runs it reaches; nothing is shared between observations of one candidate. Where that time
   goes was not profiled. *(Review, 2026-10-12: profiled — two logarithms per letter comparison
   took 17 in 100 of the time, and the layout rebuilt per read most of another 15. With the
   layouts kept, the scores looked up once a read and fixed pieces memoised (§7.6), the same
   calling takes 5.7–9.4 s
   of CPU against 13.9–23.8 s before and 1.9–3.0 s for step 4, the ranges spanning the host's
   background load: about 8–9 ms more a bundle locus, from about 24. The calls are byte-identical
   on all three benchmark runs and HG002 30×. What is left is about 200 letter comparisons per
   partial read per candidate; `str_bundles_partial_reads_2026-10-12.md` §6.1.)*

   **Step 4c (2026-10-12): the review's findings M2 and M3** — a partial read neutral about the
   runs it never reached, and the flank a whole-hull partial held carried to scoring (§7.6). New
   psps (format 2.1) for every run; the step 4b build reads them as 2.0 would, and calls them as
   before, so both rows are over the same psps. Bundle cells with a variant: right / wrong /
   missed; then false calls.

   | | 1 sample, PCR, 5× | 1 sample, PCR, 30× | 50 samples, PCR, 30× | HG002 30× |
   |---|---|---|---|---|
   | step 4b, reviewed | 99 / 29 / 80; 0 | 207 / 1 / 0; 0 | 10,400 / 13 / 22; 1 | 862 / 150 / 119; 20 |
   | step 4c | 99 / 29 / 80; 0 | 207 / 1 / 0; 0 | 10,407 / 11 / 17; 1 | 862 / 151 / 118; 20 |

   **Every STR segment cell is unchanged** against main on all three runs (3,000, 3,000 and
   150,000 cells). Records whose line changed, mostly in quality: 59 of 1,898, 151 of 2,591, 452 of
   6,005, and 1,336 of 20,348 on HG002. No whole-hull partial fell back to the average of both
   flanks (§7.6).

   **Cost.** One thread, 50 samples at 30×, the 500 bundle loci, the two builds run alternately on
   the same psps three times each on a host under heavy other load (load average 55): 12.8–13.3 s of
   CPU for step 4b, 17.6–18.0 s for step 4c — about 1.37 times, some 9–10 ms more a bundle locus.
   The extra work is the sequencing indels each run now offers inside the run-by-run walk: each
   lands the next piece at a new offset, which the walk then follows through every later run.
   *(The absolute times are about twice the review's session for the same step 4b build, from the
   load; the ratio is the measurement.)*

   *(Follow-up review, 2026-10-12: `str_bundles_partial_reads_followup_2026-10-12.md`.)* The HG002
   column is not over shared psps: each build calls HG002 from the BAM, and its step 4b cell is
   from an earlier run; the step 4c VCF also holds one record more (20,348 against 20,347). The
   review's two changes — each run's bases compared at that run's substitution rate (§7.6) and the
   forward walk that drops states within 10⁻¹² — **change no call**: the VCFs before and after are
   byte for byte the same on all three runs and on HG002 30×, so the table above stands. (Every rate
   on these runs is the defaulted constant, so the rate change cannot show here; it shows only on a
   run with fitted rates.) **Cost**, the same 500 bundle loci at one thread, three alternating runs
   of each build: 6.8–6.9 s of CPU with both changes, against 16.4–17.1 s for step 4c as first built
   and 12.0–12.3 s for step 4b — 0.41 times step 4c's and 0.57 times step 4b's, about 20 ms less a
   bundle locus. The forward walk alone, dropping nothing, measured 14.1–14.4 s: the dropped states
   are what buys most of it.

5. **The simulator's harder bundles** (§9): abutting, overlapping and impure repeats, measured as in
   step 4.

## 11. Decided (owner, 2026-10-11)

1. **The size limit is the user's**: `--max-str-bundle-len`, default 100 (§5.2).
2. **Stutter is modelled in every run** (§7.3), and **every allele has a count per run** (§7.1).
3. **Runs, not members, carry the stutter**: runs of pure whole copies of one motif, found in the
   reference hull, non-overlapping, longest first, each at least the typing's copy floor (§6.1).
   *(In review, the same day: every motif of period 1 to 6 is looked for, not only the members',
   which missed the `AT` run of §5.1's first example.)*
4. **Run counts are written as `REPCNB`** on bundle records (§7.5); `REPCN` keeps its type.
5. *Changed 2026-10-12 (owner): reads that end inside a bundle are scored, as a lower bound from
   the flank they held (§7.6).* The first decision left them out until the rest worked and
   bundles were measured at low depth; at 5× bundles lagged the generic path (§10, step 4), which
   was the condition for adding them.
6. **One name, *STR bundle*,** in the glossary, the specs and the code; the glossary's *cluster of
   repeats* is retired.
7. **The simulator is extended to abutting, overlapping and impure bundles** once the caller works
   (§9, §10 step 5).
