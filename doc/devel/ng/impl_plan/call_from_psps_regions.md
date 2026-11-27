# `call-from-psps --regions` — design and implementation plan

**Status:** built 2026-11-27, steps 1–7 (report: `doc/devel/reports/implementations/call_from_psps_regions_2026-11-27.md`). Branch `psp-call-regions`, cut from `main` at 95569d47.

## 1. What this is for

The owner's cluster workflow:

1. each sample writes **one whole-genome psp** with `generate-psps`;
2. `estimate-parameters` fits **once over all of them**, so the rarer STR classes have as much data as
   the cohort holds;
3. **one `call-from-psps` job per chromosome**, each reading the same whole-genome psps with
   `--regions chrN.bed`;
4. the per-chromosome VCFs are joined with `bcftools concat`.

Today step 3 is impossible: `call-from-psps` has no `--regions`, by design. The module comment
(`src/cli/call_from_psps.rs:11-17`) and spec `run_streaming.md` §5.3 (line 683, *"the analysed
regions come from the headers, not from a flag"*) give the reason: a flag could contradict the files.
**A set of regions that lies inside the ground the psps cover contradicts nothing**, and spec §11
question 5 (line 1180) already names this case as legitimate (*"calling one chromosome from
whole-genome psps"*). This plan adds the flag on that condition.

**Vocabulary used below.**

- **The psps' ground**, *G*: the regions every psp of the cohort records as walked. The cohort's files
  must already agree on it (`RunError::AnalysedRegionsDiffer`, `src/run/mod.rs:595`).
- **The requested regions**, *R*: the BED given to `--regions`, after the BED reader has sorted and
  merged it (`src/regions.rs:113-145`).
- **A cut**: an end of a span of *R* that does not fall at a contig's end — a place where *R* stops
  and *G* goes on.
- **A locus's first base**: the first reference position of a cohort locus, the merge's unit of
  calling. The merge already gives every locus to the building region holding its first base
  (`src/run/cohort_merge/frontier.rs:17-45`).

## 2. The owner's hard requirement: the subset check

**Every base of *R* must lie in *G*.** If any requested span, or part of one, lies outside, the run
stops before any calling, naming each offending piece and what the psps cover. *R* equal to *G* is
allowed and calls exactly what a run without the flag calls.

### 2.1 Where it sits

Today the order is: command line → log opened → **reference read and verified**
(`call_from_psps.rs:426-446`) → `--explain-loci` BED → **psps opened** (`:467`) → catalog and
segments → parameters → calling.

Opening the psps reads only each file's footer, index and header and decodes no block
(`src/run/psp_caller.rs:31-34`). It needs nothing from the reference: the per-file comparison with
the reference's checksums happens later, in `PspVariantCaller::open` (`psp_caller.rs:515-526`). So
**the cohort can be opened before the reference is read**, and the BED resolved against the psps'
own contig table (`Header::contigs`, `src/psp/header.rs:48`, through
`segmentation_section::contig_bounds_of`, `src/psp/segmentation_section.rs:417`). That table is the
reference's — every file is checked against the reference before a block is decoded — so a contig
name means the same thing either way, and the refusal arrives in seconds instead of after a
whole-genome FASTA has been read and checksummed.

New order: command line → log → **psps opened → `--regions` read and checked** → reference →
`--explain-loci` → catalog and segments → parameters → calling.

### 2.2 The check and its message

A set difference, *R* minus *G*, computed on the merged spans. Each piece left over is reported
with what *G* holds on the same contig.

```
--regions chr_extra.bed asks for ground the psps were not written over, so nothing was called:
  chr3:1-2000000            the psps hold nothing on chr3
  chr2:900001-1000000       on chr2 the psps hold chr2:1-900000
The psps were written over 2 intervals, chr1:1-1000000 … chr2:1-900000 (1.9 Mb).
Remove these from the BED, or write the psps again over ground that includes them.
```

- Coordinates are 1-based and inclusive, as every other message and the run report write them
  (`src/run/report.rs:833-848`, `name_of`). The psps' ground is summarised with the report's own
  `describe` (`report.rs:819`).
- At most ten pieces are listed, then `… and N more`.
- A contig name the psps' table does not have is the BED reader's existing refusal
  (`BedError::UnknownContig`, `src/regions.rs:363-365`), raised before the subset check.
- A span overhanging its contig's end is clamped by the BED reader (`regions.rs:110-112, `:277``), so a BED
  written as `chr1 0 999999999` asks for all of `chr1` and passes against whole-contig psps.

**Home:** a new `GroundError::RegionsOutsideThePsps` in `src/cli/run_ground.rs`, rendered
transparently like the module's other refusals (`run_ground.rs:15-17`), and a function there that
reads the BED, checks it and returns *R*. The set difference itself is a method on `GenomeRegions`
(`src/region_typing/mod.rs:83`), `parts_outside(&self, covering: &GenomeRegions) -> Vec<GenomeRegion>`.

## 3. Where the run's ground flows, and what changes

The ground is read from the headers at `call_from_psps.rs:469` and goes to six places.

| Consumer | Today | With `--regions` |
|---|---|---|
| Segmentation (`run_ground::segments_cut_with`, `call_from_psps.rs:475-494`) | cut over *G* | cut over *R*. In psp mode the segments themselves are read by nothing but the run report's satellite count (`segments()` is used only by the walker and the gatherer, `src/run/walker.rs:357`, `src/run/gatherer.rs:332`); what the merge walks is `Segmentation::analysed_regions` (`src/run/segments.rs:176`). Cutting over *R* also reads only *R*'s part of the catalog. |
| The per-file check in `PspVariantCaller::open` (`psp_caller.rs:515-526`) | `first_difference` compares the setting, catalog, criteria **and the analysed regions** (`src/segmentation_inputs.rs:190-208`) | compare the first three as before; replace the regions comparison by "the run's regions lie inside every file's". The ground half gets its own `RunError` (defensive: the command has already refused). |
| The merge (`merge_and_call_streaming` and the round driver, `src/run/cohort_merge/stream.rs:139`) | walks *G* | walks *R*: it already takes the analysed regions from the segmentation. |
| The psp sources (`PspSummarySource::over`, `src/run/psp_source.rs:408`) | every file read from its first block | start at *R*, and skip what lies between its spans — §4. |
| Run report (`report.rs:303-311`) and first log lines (`call_from_psps.rs:670-711`) | *"analysed ground: … as every file's header records them"* | name both: the ground called (*R*, with the BED's path) and the ground the psps cover (*G*). |
| VCF header (`calling_run::header_for`, `src/cli/calling_run.rs:1862`) | every reference contig, the command line | **unchanged.** Every job then writes the same `##contig` lines, which `bcftools concat` wants; the `##commandline` line records `--regions`. |

**The hidden-duplication filter and refinement need no code change**, but they work per run: the
filter fits each sample's depth model from the windows this run read, and estimates how common
hidden duplications are from this run's calls (`calling_run.rs:2097-2117`). A per-chromosome job
therefore fits them per chromosome. §5.3 says what that changes.

## 4. Reading only the requested ground

### 4.1 What happens if nothing is changed

A psp source walks its file from the first block (`psp_source.rs:408-416`, `reader.rs:312-320`).
With the merge walking chromosome 5 only, the observation cache
(`src/run/cohort_merge/observation_cache.rs`) would draw every record of chromosomes 1 to 4 at the
first cover, hold them, and measure window coverage over them — fetching those chromosomes' reference
bases to do it (`observation_cache.rs:1285-1296`, `:1714-1790`). That is slower than the
whole-genome run and puts other chromosomes' windows into the hidden-duplication filter's depth model.
So the sources must skip.

### 4.2 Seek to the first span

The block index exists for this: one entry per block, decoded at open, *"so a reader can start at
any block without reading what comes before it"* (`src/psp/index.rs:1-2`). `PspReader::records_from`
(`src/psp/reader.rs:396-418`) already turns a position into the earliest block that can hold a record
starting there, with one binary search, and is tested (`src/psp/walk.rs:716-1891`). No run calls it
today.

**Its one trap is documented at `reader.rs:383-394`: it selects records by where they start.** A
record that begins in an earlier block and reaches into the span would be missed. So the source
seeks to the first span's start **minus the cohort's observation reach ceiling** — the widest span
any record in these files can have (`OpenPspCohort::observation_reach_ceiling`,
`psp_caller.rs:265`) — which reaches back over every record that can overlap it.

**Cost:** one binary search over at most ~14,000 index entries, and at most one block read that
holds nothing the run keeps. A chromosome job then decompresses that chromosome's blocks and no
others. The prefetch thread needs no change: it starts at the block holding the merge's first
position (`src/run/psp_prefetch.rs:121-142`, `:277-305`). It will read up to one share past the
last span's end; capping it there is optional.

### 4.3 Drop what the merge has already passed

The merge drivers evict, before every round's cover, each held record ending before the round's
first base (`stream.rs:461-478`; `rounds.rs:193-195`; `serial.rs:320`, `:1103`). **Add a hint to the
source interface** (`ObservationSource`, `observation_cache.rs:190`):

```rust
/// Records ending before `position` will not be read. A source may drop them unread.
fn nothing_wanted_before(&mut self, position: GenomePosition) {}
```

`ObservationCache::evict_before` and `evict_before_in_parallel` (`observation_cache.rs:1127`,
`:1161`) pass their position to every source. The walker keeps the default and ignores it. The psp
source records it as a floor and, in `next_summary` (`psp_source.rs:518-575`), drops each record
ending before the floor. It drops it the way the depth ceiling already drops a record — bytes out of
the arena, `HeldLiveSets::skip` (`psp_source.rs:552-562`) — so the live read sets stay consistent.

**Why this changes nothing in a run without `--regions`.** A cover draws each sample until its next
record starts past the chain's reach, which is at least half a window (250 bases) past the round's
last base (`observation_cache.rs:1298-1361`). So when the next round evicts at its own first base,
the source's next undrawn record already starts after that base and the floor drops nothing. It
drops records only where the analysed ground has a gap: between two spans of *R*, and before the
first.

**Gaps between spans are read through, not skipped by a second seek.** Re-seeking a source mid-run
means restarting its walk under an arena and live sets the cache still points into; the per-chromosome
workflow has one span per job and never needs it. A BED of many intervals spread over the genome
costs at most one read of the file from its first span on — never more than the whole-genome run
reads. A mid-run re-seek is a follow-up for when a measured BED needs it.

## 5. Edges: what a cut does, and the rule

### 5.1 How the walking commands treat a BED edge

`generate-psps --regions` and `call-from-alignments --regions` cut the reference into segments over
the BED (spec `run_streaming.md` §4.2). A repeat tract, bundle or satellite that overlaps the BED is
**kept whole, even past the edge**; ordinary sequence is **clipped** to the BED
(`src/region_typing/mod.rs:693-731`). The pileup walker then reads a halo past each segment's end so
that a record anchored inside is completed (`src/locus_generation/pileup/generator.rs:577-580`,
`:683-706`). So a walk over *R* holds records that start inside *R*'s segments, and none that start
in the clipped-off ordinary sequence.

### 5.2 Recommended rule: a run over *R* is the slice of the whole-ground run that *R* owns

**A locus is called by the run whose regions hold its first base. It is called from the same
evidence the run over all of *G* would use.** Concretely: the merge walks *R*; the sources hand it
every stored record that ends at or after the current span's first base, and do not stop at a
span's end, so a locus starting just inside the end is still built from the records that reach past
it.

**This is not a new rule. The merge already does it at every analysed region's edge.** A builder
discards a locus that opened before its own ground (`src/run/cohort_merge/build.rs:971-979`), the
last building region of an analysed region claims no locus that starts past its end
(`serial.rs:1383`, `a_locus_past_the_analysed_ground_is_not_claimed_by_the_last_building_region`),
and dividing the analysed ground differently changes no locus
(`serial.rs:765`, `the_same_loci_come_out_however_the_analysed_ground_is_divided`). What is new is
only that the psps hold records outside the analysed ground, and §4 keeps those from costing time.

**What it gives, measured against the two references a run could be compared with:**

| | Whole contigs in *R* | A cut inside a contig |
|---|---|---|
| **Against `call-from-psps` on psps walked over *R*** | **identical VCF** | differs within about one locus of the cut (below) |
| **Against one run over all of *G*, the per-*R* VCFs concatenated** | identical records with the hidden-duplication filter off; with it on, see §5.3 | identical calls, except as listed below |

**Why whole contigs are exact.** Nothing crosses a contig: a segment never does (`segments.rs:161`),
a cohort locus never does, and neither does a read, a phase set, or the gathering of an STR record
widened over its flank. The one thing the two walks number differently — the read identifiers, where
a whole-genome walk starts chromosome 5's reads after chromosome 4's — is already compared by rank
from each record's first identifier, and a test holds a `--regions` walk and a whole walk to the
same reads (`src/cli/mode_equivalence.rs:863-920`, `the_cap_keeps_the_same_ranks_whole_contig_and_regions`).

**Where a cut inside a contig differs from a walk over *R*.** Each difference is confined to one
locus span of the cut — at most `--max-locus-span` (default 100 bases) — or, for the last item, to
the run as a whole:

1. **A repeat tract starting before the BED's first base and reaching into it** is in a walk-over-*R*
   VCF, which keeps the tract whole. Under this rule it belongs to the run whose regions hold its
   first base. **This is the one place the rule disagrees with `call-from-alignments --regions`**,
   and it is what keeps `bcftools concat` free of duplicate records: with overlap-selection, two
   adjacent jobs would both write that tract.
2. **A deletion read across the BED's first base.** The whole-ground walk stored one record starting
   before the cut; this run gives the locus it opens to the previous region. A walk over *R* starts
   afresh at the cut and can call a variant in the deletion's first bases.
3. **A locus that starts just before the BED's last base** is built here from records reaching past
   the cut, as the whole-ground run builds it. A walk over *R* has no records past the cut to chain
   with.

**Where the concatenation differs from one run over *G* at a cut inside a contig:**

- **A phase set spanning the cut is split in two**, because each run phases only its own records.
- **An STR record widened over its flank is not merged with a generic record on the other side of the
  cut** (`src/run/edge_deletion_gathering.rs:1-30`, which holds records 256 bases back,
  `HOLD_REACH` at `:44`). The concatenation can then hold two records over one base, which one run
  would have merged.
- **The window depth at loci in the first 250 bases after the cut** is measured without the records
  before it. It reaches no VCF field (`src/vcf/assemble.rs:73-79`) and only the hidden-duplication
  filter reads it.

**Why not restrict `--regions` to whole contigs.** Arbitrary intervals are safe: every difference
above is near a cut, none produces a genotype the evidence does not support, and the same effects
already follow from `--regions` on the two walking commands. They cost nothing extra either: §4's
seek is by position, not by contig. So intervals are allowed, and the documentation recommends whole
contigs for splitting a genome into jobs because that is where the concatenation is exact.

### 5.3 What per-run filtering changes, even for whole contigs

The hidden-duplication filter is fitted per run (§3). A per-chromosome job fits each sample's depth
model on that chromosome's windows, takes the chromosome's rate of hidden duplications as the rate
every local estimate is pulled toward (`--hidden-dup-window`, `calling_run.rs:424-435`), and controls
the false-discovery rate within the chromosome. Controlling it at 1% in each job keeps it at about 1%
over their union. Two cases deserve a line in the user documentation: **a sex chromosome or an
organelle**, whose one-copy depth differs from the autosomes', gets a depth model fitted to itself;
and **a short contig** gives the fit few windows and few calls, so its estimates lean on the
defaults. Neither is new — `call-from-alignments --regions` has both today.

## 6. The parameters file needs no change

A parameters file binds to the reference and to the cohort's samples and read groups
(`refuse_if_not_this_runs_inputs`, `src/calling/parameters_file/bindings.rs:584`). It also records
the census it was fitted from, whose terms include a digest of the analysed regions
(`SelectionTerms::analysed_regions`, `src/parameter_estimation/joint/loci.rs:524-525`). **A calling
run never compares that census**: `calling_run::run_parameters` passes `None` for it
(`calling_run.rs:1720-1721`), so the analysed-region digest is not checked at calling. A
whole-genome fit therefore already binds to a one-chromosome run, as it already does for
`call-from-alignments --regions`. Step 5's test pins this, so that a later change binding the census
at calling fails a test rather than this workflow.

## 7. Implementation steps

Each step builds, passes its own tests and the touched modules' existing tests, and changes no VCF
of a run without `--regions`.

1. **The set difference.** `GenomeRegions::parts_outside` in `src/region_typing/mod.rs`. Unit tests:
   *R* inside *G* gives nothing; *R* equal to *G* gives nothing; overhang at the left end, at the right
   end, and at both; a span bridging the gap between two spans of *G* gives the gap; a contig *G* does
   not cover gives the whole span; edges one base apart (*R* ending on *G*'s last base, *R* starting
   one past it).
2. **The refusal and the flag.** `--regions` on `CallFromPspsArgs` (`call_from_psps.rs:77`);
   `GroundError::RegionsOutsideThePsps` and the reading-and-checking function in `run_ground.rs`;
   the cohort opened before the reference (§2.1). Tests in `src/cli/call_from_psps/tests.rs`:
   replace `the_ground_cannot_be_narrowed_by_a_flag` (`tests.rs:78-91`) with a test that the flag
   parses; psps walked over a BED and a `--regions` partly outside it are refused, the message names
   the piece outside and the psps' ground, **and no VCF is written**; a BED naming a contig the
   reference lacks is refused by the BED reader's message.
3. **The caller accepts ground inside the files'.** Split `SegmentationInputs::first_difference` so
   `PspVariantCaller::open` compares setting, catalog and criteria, then checks the run's regions lie
   inside every file's (new `RunError` variant). `call_from_psps.rs` cuts the segmentation over *R*.
   Tests in `psp_caller.rs` beside `a_checked_cohort_opens_as_a_caller` (`:2338`): a caller over a
   sub-ground opens; one whose ground reaches outside is refused naming it; the existing
   different-catalog and different-criteria refusals still fire.
4. **The sources skip.** `ObservationSource::nothing_wanted_before` with a no-op default; the two
   `evict_before` functions call it; `PspSummarySource` seeks at construction to the first span's
   start minus the reach ceiling (`psp_caller.rs:706-736` builds the sources) and drops records
   below the floor. Tests: in `psp_source.rs`, on the multi-block fixture (`a_psp_of`, `:1334`), a
   source started at a position and floored there yields exactly the records ending at or after it,
   and begins at the expected block (`blocks_begun`); in `serial.rs`, a source holding records in a
   gap between two analysed regions gives the same loci as one holding none; the existing driver
   comparisons stay green, which is the "drops nothing without `--regions`" claim.
5. **End-to-end equivalence**, in `src/cli/mode_equivalence.rs`, which already walks psps over a BED
   (`walked_into_psps_over`, `:193`):
   - **Whole contigs, exact.** A new two-contig varying cohort in `src/cli/test_fixtures.rs` (the
     present one has one contig, `test_fixtures.rs:127`). `--regions` naming the second contig on
     whole-genome psps writes the **same VCF, every line**, as psps walked over that contig alone —
     with the hidden-duplication filter on and refinement on, so the per-run fits are compared too.
   - **The partition.** On the one-contig cohort, cut at a base inside ordinary sequence and at a base
     inside the repeat tract (`TRACT`, `test_fixtures.rs:137`): the two runs' records, concatenated,
     equal the whole run's records, with the filter off. Assert the record count first, as the
     module's other oracles do (`mode_equivalence.rs:260-268`).
   - ***R* equal to *G*** writes the same VCF as no `--regions`.
   - **A whole-genome parameters file calls a one-contig run** (§6).
6. **Report and log.** First log lines and the run report name *R* with its BED and *G*
   (`call_from_psps.rs:670-711`, `report.rs:303-311`). Test: extend `the_report_a_person_sees`
   (`tests.rs:280`).
7. **Docs** (§8).

## 8. Documentation to change

Code and developer docs, in this branch:

- `src/cli/call_from_psps.rs:11-17` — the module section "The ground is the files', not a flag's":
  the files say what ground exists, `--regions` chooses a part of it, and the subset rule is why the
  original argument still holds. Also the order in `run_call_from_psps`'s `# Errors` (`:367-373`).
- `src/run/psp_caller.rs:11-14` — "a run over psps is told nothing".
- `CallFromPspsArgs::regions` help text: *"A BED file of the regions to call, inside the ground the
  psps were written over; without it, all of that ground. A locus belongs to the run whose regions
  hold its first base, so runs over the pieces of a genome, joined with `bcftools concat`, write each
  locus once. Splitting by whole chromosomes gives exactly the calls of one run."*
- `doc/devel/ng/spec/run_streaming.md` §5.3 (line 683) and §6.2 (the cohort check stays; add the
  subset rule against the run); §11 question 5 (line 1180): the one-chromosome case is settled by
  this; a cohort walked over different ground stays open.
- `doc/devel/ng/spec/command_line.md` §4.2 table (line 156): `--regions` used by
  `call-from-psps` too; §5.2 (line 343), which lists `--regions` as alignments-only.

User docs — **not edited here**, because the uncommitted `user-docs` branch edits the same files:

- `doc/user/calling.md:77-78` — replace *"`call-from-psps` has no `--regions`"* with the flag, the
  subset rule, and the per-chromosome cluster recipe (generate whole-genome psps, fit once, one call
  job per chromosome, `bcftools concat`). Say that whole chromosomes are exact, what a cut inside a
  chromosome changes (§5.2), and that the hidden-duplication filter fits per job (§5.3).
- `doc/user/options.md` — `--regions` under `call-from-psps`.
- `doc/user/recipes.md` (on `user-docs`) — the cluster recipe, if that file carries recipes.
- `doc/user/limits.md:36` — add the refusal for `--regions` outside the psps.

## 9. Decisions for the owner

1. **The edge rule: a locus belongs to the run holding its first base** (recommended), against
   the walking commands' rule that a repeat tract overlapping the BED is reported whole. The
   recommended rule writes every locus exactly once across adjacent jobs; the other writes a
   straddling tract in both jobs and diverges from the whole run. The cost is that
   `call-from-psps --regions` and `call-from-alignments --regions` disagree about a tract that
   starts before the BED's first base.
2. **Open the psps before reading the reference** (recommended), so the subset refusal comes in
   seconds. It moves a psp-open failure ahead of a reference-read failure in the order of refusals;
   nothing else changes.
3. **The hidden-duplication filter fits per job** (§5.3). Recommended: accept and document. Fitting
   it once over the joined VCF would be a separate feature.
