# ng — the cohort merge: types and interfaces

*Status: architecture draft (2026-08-17), companion to the spec
[`../spec/cohort_merge.md`](../spec/cohort_merge.md) (the design and every *why*) and to the run
architecture [`run_streaming.md`](run_streaming.md) (the objects this machinery lives inside),
with the shared vocabulary in [`ng_step_interfaces.md`](ng_step_interfaces.md). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): nouns for types,
verbs for functions, **STR** in prose ↔ `ssr` in code. Signatures are illustrative; the
**contract** is the deliverable. This document does not re-argue a decision — the spec section is
cited instead.*

Everything here is crate-private machinery inside the two caller objects
([`run_streaming.md`](run_streaming.md) §3); nothing in this module is public API of the crate.

## Module home

`src/ng/run/cohort_merge/`, a folder inside the run module
([`run_streaming.md`](run_streaming.md)'s tree):

```
src/ng/run/cohort_merge/
├── mod.rs               – the five run parameters, the shared analysed-region guard,
│                          the module's test fixtures
├── close.rs             – LocusCloser: the reach walk that closes loci and judges them,
│                          the two per-locus verdicts
├── build.rs             – build_region: assembling survivors into CohortObservations;
│                          CohortObservation, SampleSupport, RegionOutcome live here
├── observation_cache.rs – ObservationCache: one forward reader per sample and the window
│                          the builders read; building_regions_of, the geometry
├── organise.rs          – Organiser: ordered release and overlap resolution
├── serial.rs            – the oracle, and the same merge read through the cache
└── parallel.rs          – the builders working a round of regions at a time
```

**Amended 2026-08-18, after milestone E**, in three places: the cache is its own file, not part
of `organise.rs` (the argument for keeping them together assumed the organiser would be the
cache's only writer, and the cached serial driver is a second one); the observation types live
with the code that builds them rather than in `mod.rs`; and the two drivers of milestone C and
the parallel arrangement of milestone E are files this tree did not have.

`close.rs` holds everything decided before anything is assembled; `build.rs`
holds the only code that touches heavy columns. **This revises
[`run_streaming.md`](run_streaming.md) (arch) §3.2**, whose `calling.rs` sketched a `k_way_merge`
free function: that function is superseded by this module — `call_vars_in_segment` calls
`build_region` (or consumes its output) instead. Recorded rather than quietly changed, because
that document still shows the old shape.

---

## The algorithm, end to end

One pass through, before any type appears. The *why* for each step is the spec section named;
nothing here re-argues one.

1. **The organiser divides the genome into building regions** — `cohort_locus_builder_regions_len`
   bases each, 20 by default — and never puts a boundary inside an STR or bundle segment, where a
   builder would have no locus to start (spec §6.1).
2. **It fills the observation cache** to cover the regions in play, drawing each sample's one
   forward reader along. This is the module's memory, which is why the regions are short and the
   builders stay close together (spec §6.4).
3. **A builder takes one region and walks the samples' observations merged by position** (spec
   §4.1). It keeps where the open locus starts, how far it reaches, and — **per sample** — that
   sample's non-reference reads and the reads they were compared against. Each observation
   extends the reach if it goes further; when the next one starts beyond the reach, the locus
   closes.
4. **It judges each closed locus, width first.** Wider than `max_cohort_locus_span` → failed:
   not assembled, counted, and its ground still displaces overlapping loci (spec §3.2). Otherwise
   no single sample reaching `min_alt_reads` → dropped, not counted, ground judged empty (spec
   §4.3). Otherwise built.
5. **It assembles the survivors** — every sample's observations projected onto the locus span and
   identical projections unified into one allele table, each sample's support expressed against it
   (spec §4.2).
6. **It starts loci only inside its own region but follows one past the end** if a deletion carries
   it there, so builders overlap by design (spec §6.1).
7. **The organiser resolves the overlaps**: of two overlapping loci the one starting earlier
   stands, the other is dropped — a failed locus wins the same way an emitted one does (spec §6.1).
8. **It releases in genome order**, holding a region until its predecessor has arrived, sums the
   failed counts, and evicts the cache behind what it has released (spec §6.3).

**Nothing coordinates the builders.** They share no mutable state, exchange no messages, and never
learn whether their start position was inside someone else's locus. Correctness comes from the
ownership rule plus step 7, not from where the work was cut (spec §5).

## What to read in production first

This module is largely production's algorithm with its columns removed. Read these before writing
it; each row says what carries over and what does not.

| production | where | how it relates |
|---|---|---|
| the reach walk that groups positions into variant groups | `derive_is_kept` [`cohort_integration.rs:166-187`](../../../../src/var_calling/cohort_integration.rs), `reach` [`:46-48`](../../../../src/var_calling/cohort_integration.rs) | **the closest ancestor** — step 3 is this walk, consuming merged observations instead of decoded columns |
| projecting each sample's alleles onto the group span and unifying them | `PerGroupMerger` [`per_group_merger.rs:585`](../../../../src/var_calling/per_group_merger.rs), module doc [`:1-20`](../../../../src/var_calling/per_group_merger.rs) | step 5 is this, minus the likelihoods it also computes — here the builder emits evidence only |
| the k-way merge across samples | `PerPositionMerger` [`per_position_merger.rs:145`](../../../../src/var_calling/per_position_merger.rs) | **upstream of this module**, and what makes step 3 a single pass ([`run_streaming.md`](run_streaming.md) §3.2) |
| emitting out-of-order results in genome order | `VcfWriter`'s reorder map [`vcf_writer.rs:162-176`](../../../../src/var_calling/vcf_writer.rs), gap guard [`:152-158`](../../../../src/var_calling/vcf_writer.rs) | step 8's structure, carried whole |
| the keep threshold | `max_nonref_obs` per sample in `derive_is_kept`, default 2 [`var_calling/mod.rs:72`](../../../../src/var_calling/mod.rs) | step 4's second verdict, per sample as production's is, with a share added for depth |
| cutting work at gaps no group can span | `merge_block_ranges` [`cohort_integration.rs:403-430`](../../../../src/var_calling/cohort_integration.rs) | **not carried** — this design overlaps and resolves instead (spec §5) |
| the staged producer/caller pipeline | [`pipeline.rs:1-30`](../../../../src/var_calling/pipeline.rs) | **not carried** — production's fold is serial on one thread; here it is what parallelises |
| deferring the expensive columns until the keep decision | `TwoPhaseSegment` [`sample_reader.rs:698`](../../../../src/var_calling/sample_reader.rs), `set_variable_rows` [`:789`](../../../../src/var_calling/sample_reader.rs) | the psp path only; the direct path defers nothing (§2) |

## 1. The constants (`mod.rs`)

Two of them are spans with different jobs and are separate types so they cannot be swapped
(spec §11). The third and fourth are the two halves of the keep rule, bundled into one value so
that no call site can pass a floor and a share that were not chosen together. **The fifth was
added by milestone E** and is described at the end of this section.

```rust
/// The policy bound: the widest cohort locus the caller undertakes to build, in reference
/// bases. A wider locus fails — counted, never built (spec §3.2, the owner's 2026-08-17
/// ruling).
///
/// **A command-line parameter of a calling run**, default 50. Never recorded in a psp — the
/// files hold what the generator minted — so re-calling under a new value needs no re-walk
/// (spec §3.1). It *is* recorded in the run's output beside the failed-locus count, because
/// it decides which ground was refused and two runs over the same psps under different
/// values are otherwise indistinguishable (spec §3.1, §3.3).
pub struct MaxCohortLocusSpan(pub NonZeroU32);

/// **How much non-reference evidence one sample must show for the cohort to build a locus**
/// (spec §4.3). A locus no single sample reaches is dropped: not assembled, not emitted, and
/// **not counted** as a failure — a failure is ground the caller refused, and this is ground it
/// judged empty.
///
/// Two numbers, for the two ends of the depth axis: the floor decides at low coverage, where a
/// share of three reads rounds to nothing, and the share decides at high coverage, where two
/// reads out of three hundred is the error rate rather than an allele. Asked of the sample's own
/// reads, never the cohort's — a share of the cohort's would raise its bar as samples are added
/// while a rare allele's evidence stays where it is (spec §4.3).
pub struct MinAltReads {
    pub floor: MinAltObs,
    pub share: MinAltReadShare,
}

impl MinAltReads {
    /// The floor, or the share of `reads_compared_with_reference` rounded up, whichever is more.
    pub fn required_of(self, reads_compared_with_reference: u32) -> u32;
    /// Whether one sample's counts reach it — `required_of` with the comparison folded in.
    pub fn reached_by(self, non_reference_reads: u32, reads_compared_with_reference: u32) -> bool;
}

/// The floor half. **A command-line parameter**, and the reason it exists is measured rather
/// than aesthetic: in production it removes a large number of very low-quality variants and
/// improves performance substantially. Its cost is that a variant no sample showed twice is
/// unrecoverable.
pub struct MinAltObs(pub NonZeroU32);

/// The share half, as a fraction of one. **A command-line parameter.** Constructed through a
/// checked `new`, which refuses anything that is not a fraction of one rather than clamping.
pub struct MinAltReadShare(f64);

/// Production's default, carried over ([`var_calling/mod.rs:72`](../../../../src/var_calling/mod.rs)).
pub const DEFAULT_MIN_ALT_OBS: u32 = 2;

/// Two reads in a hundred — the owner's number (spec §4.3), chosen against loci-kept counts on
/// two benchmarks a hundredfold apart in depth and not yet against a truth set.
pub const DEFAULT_MIN_ALT_READ_SHARE: f64 = 0.02;

/// How many reference bases one builder's region covers (spec §6.1). **A command-line
/// parameter, default 200** (raised from 20 on 2026-08-18, spec §14 question 1). It is not
/// derived from `max_cohort_locus_span`: what it really costs is the observation cache, which
/// must cover every region in play at once, so `builders × this` is the ground held resident
/// (spec §6.4, §8).
pub struct CohortLocusBuilderRegionsLen(pub NonZeroU32);
pub const DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN: u32 = 200;

/// The default, and it is soft — the owner's number, unmeasured (spec §14 question 3).
pub const DEFAULT_MAX_COHORT_LOCUS_SPAN: u32 = 50;

/// Its default. The owner's number, unmeasured — soft (spec §14 Q2).
pub const DEFAULT_MAX_COHORT_LOCUS_SPAN: u32 = 50;

/// The widest reference span any observation in the input can have. **No rule in this module
/// depends on it** (spec §5): it bounds only how far past a builder's ground the observation
/// cache may have to reach, which is a memory fact (spec §6.4, §8). Direct mode reads the
/// generator's `max_record_span` (`pileup/generator.rs:93`).
pub struct ObservationReachCeiling(pub NonZeroU32);
```

The psp header therefore gains one recorded value — the writing run's ceiling. **This revises
[`run_streaming.md`](run_streaming.md) (arch) §4's `PspHeader`, flagged there and in spec §13;
the run spec's §6.1 does not carry the field yet.**

---

**The fourth, added 2026-08-18 with the parallel arrangement.**

```rust
/// How many building regions the merge works at once (spec §6.2).
///
/// **A count of regions in flight, not of threads.** The threads come from rayon's pool; what
/// this number sets is the ground the observation cache must hold, which is
/// `cohort_locus_builder_regions_in_flight × cohort_locus_builder_regions_len` bases plus the
/// tail of the observations reaching past it (spec §6.4, §8).
pub struct CohortLocusBuilderRegionsInFlight(pub NonZeroUsize);
```

**It is the only one of the four with no constant default**, and deliberately: what it should
be depends on the machine's cores and on how much memory the cohort's width leaves, neither of
which is knowable where the other three defaults are written. A run given no value takes
`one_per_worker_thread()` — one region in flight per thread in rayon's pool. Like the other
three, its resolved value belongs in the run's output, because two runs differing only in it
differ in memory and in nothing a reader of the output can see.

## 2. What a builder reads

**A builder never touches a reader.** It reads observations from the cache the organiser owns (§4),
which holds one forward reader per sample and is the only thing in this module that pulls. So this
document declares no source trait: `ObservationSource` belongs to
[`run_streaming.md`](run_streaming.md) arch §2, which owns it because it is what lets one calling
stage serve both a walk and a psp reader. Here it appears only as what the cache was handed.

**It needs two numbers per position, and both belong on the observation type rather than in this
module.** Neither is a fact about merging; both are facts about an observation, and one of them
already has a second caller.

```rust
impl SequenceObservation {
    /// Whether this observation's sequence is the reference's over the locus it belongs to.
    ///
    /// **The one definition of "non-reference" in the codebase.** `CensusWriter::add_generic`
    /// makes this comparison inline today
    /// ([`census.rs:2084`](../../../../src/ng/parameter_estimation/joint/census.rs)) and should
    /// call this instead — two spellings of one test are two things that can disagree.
    pub fn matches_reference(&self, reference_bases: &[u8]) -> bool;
}

impl SampleLocusObservations {
    /// The last reference base this observation covers.
    ///
    /// `region.end` today, but named: production's grouping arithmetic is
    /// `pos + max(span, 1) − 1` (`reach`,
    /// [`cohort_integration.rs:46-48`](../../../../src/var_calling/cohort_integration.rs)), and
    /// whoever compares the two should find one place where they agree rather than an
    /// open-coded expression at each use.
    pub fn reach(&self) -> Position;

    /// Reads here whose sequence is not the reference's — `matches_reference` over the members,
    /// summing `num_obs` on those that differ.
    pub fn non_reference_reads(&self) -> u32;
}
```

**The predicate is the shared thing, not the sum.** The census needs it per read group, this
module needs a flat total, and the two must agree on what "non-reference" means — so
`matches_reference` is where the definition lives and each caller keeps its own sum.

**Both are computed as the builder walks, not gathered into columns first.** A builder advances
through the samples' observations merged by position and keeps three running values: where the
current locus starts, how far it reaches, and its non-reference total. Each observation it passes
extends the reach if it goes further and adds to the total; when the next observation starts beyond
the reach, the locus is closed and judged. That is production's group walk
(`derive_is_kept`, [`cohort_integration.rs:166-187`](../../../../src/var_calling/cohort_integration.rs))
with the columns removed — extending the reach as it scans is what makes a second pass over the
locus unnecessary, however far a late observation widens it.

**Production materialises the columns; here they would be waste.** `CohortSpanFold` holds three
parallel vectors over a whole block
([`cohort_integration.rs:64-78`](../../../../src/var_calling/cohort_integration.rs)) because a psp
reader has just decoded that block and the columns are what it decoded. A builder in the direct
path is handed observations, so building vectors from them only to walk the vectors would allocate
per region for nothing. **There is no `PositionSummaries` type.**

**When the psp path arrives it may hand these over precomputed**, since a file can store what a
walk must derive, and the deferred half is fetched only for the loci that survive (run spec §3.3;
production's shape, [`sample_reader.rs:698-712,789`](../../../../src/var_calling/sample_reader.rs)).
The walk above is unchanged by that: it consumes a position, a reach and a count per observation,
and does not care which of them were read and which were computed.

---

## 3. Closing loci and judging them (`close.rs`)

One walk over the merged observations closes each locus and judges it. No cohort summary is
materialised (§2).

```rust
/// Walks the samples' observations for one building region, merged by position, and closes
/// each locus as the reach stops growing. Yields the loci in genome order.
///
/// **A locus is closed when the next observation starts beyond the current reach.** Until then
/// each observation extends the reach if it goes further and adds its non-reference reads to
/// the running total — production's group walk with the columns removed
/// ([`derive_is_kept`, `cohort_integration.rs:166-187`](../../../../src/var_calling/cohort_integration.rs)).
pub struct LocusCloser<'a> { /* the per-sample cursors, the open locus's start, reach and total */ }

impl<'a> Iterator for LocusCloser<'a> {
    /// Every closed locus, judged: what it covers, its verdict, and the observations in it.
    type Item = ClosedLocus<'a>;
}

/// A locus as it comes off the walk, before anything is assembled.
pub struct ClosedLocus<'a> {
    pub region: GenomeRegion,
    /// The members, borrowed from the cache — never copied out of it (spec §6.4).
    pub members: &'a [SampleObservationRef],
    pub verdict: Verdict,
}

/// The two verdicts, in the order they are decided (spec §4.3).
pub enum Verdict {
    /// Wider than `max_cohort_locus_span`: not assembled, counted, and its span still
    /// displaces overlapping loci (spec §3.2, §6.1).
    Failed,
    /// No single sample reached `min_alt_reads`: not assembled, not counted — ground judged
    /// empty, not refused.
    TooQuiet,
    /// Assemble it (§4).
    Build,
}
```

**Contract.** The walk is a pure function of the observations it is given and the two parameters:
the same input yields the same loci with the same verdicts, whatever building region it was
started from and whatever the reach did along the way. Reach arithmetic saturates, as production's
does ([`:46-48`](../../../../src/var_calling/cohort_integration.rs)). A locus's members are
borrowed from the cache, so closing costs no copy.

**The merge across samples happens here, over the cache's windows.** The cache holds one window
per sample (§4); the closer keeps a cursor into each and takes the lowest position among their
heads, refilling only the one it took — the argmin k-way merge the read layer already uses over
per-file streams ([`sample_reads.md`](sample_reads.md) §4). It is not a separate stage and nothing
materialises a merged sequence: the walk *is* the merge.

**One pass is enough only because the observations arrive merged by position.** A deletion widens
the locus as it is read, and the observations the widening pulls in — other samples' records at
positions now inside the locus — have not been consumed yet, so the walk simply keeps taking them
while the next position is within the reach, each possibly pushing the reach further again.

**Write it as a loop over samples at a position and that property is gone.** Then a widening sends
you back over the samples you have already passed, because they may hold records in the newly
covered bases, and each pass may widen it again. It is the merge order that buys the single pass,
not the walk itself — this is why production's version is one `while` with a growing `group_end`
and not a fixpoint, its columns being position-ordered across the cohort.

---

## 4. The builder and the organiser (`build.rs`, `organise.rs`)

A builder owns one region and works alone. The organiser holds every builder's output, resolves
the overlaps between them, and emits in genome order (spec §6.1, §6.3).

```rust
/// One builder's whole job over one region (spec §6.2): derive the position summaries from
/// the cache, fold them, close loci, apply both verdicts, and assemble the survivors.
///
/// **Takes the cache by shared reference and mutates nothing.** It holds no reader, seeks
/// nothing, and cannot disturb another builder (spec §6.4, goal 1).
///
/// **Starts a locus only inside `region`, and may finish one outside it.** The first half is
/// what assigns ownership — a locus belongs to the builder whose region holds its first
/// position — and the second is what keeps a locus whole when a deletion carries it past the
/// end (spec §6.1). So the returned loci may reach beyond `region`, and a builder reads past
/// its own end to follow them.
pub fn build_region(
    region: &GenomeRegion,
    cache: &ObservationCache,
    bound: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
) -> RegionOutcome;

/// What one region delivers to the organiser — exactly one per region, even when every field
/// is empty: the drain's index sequence must stay gapless and the counts must all arrive to be
/// summed (spec §6.3, §3.3).
pub struct RegionOutcome {
    /// The survivors, in genome order. Loci are disjoint within one region, so position is a
    /// total order (spec §9).
    pub observations: Vec<CohortObservation>,
    /// **The ground of the loci that failed `max_cohort_locus_span`, in genome order.** A failed locus is an
    /// ordinary locus everywhere but emission (spec §3.2): it owns its ground and displaces
    /// what overlaps it, so the organiser needs its span even though nothing is built for it.
    /// Without these spans the loci a neighbouring builder built inside that ground — from a
    /// view that never saw what opened there — would survive with nothing to displace them.
    pub failed: Vec<GenomeRegion>,
}

/// Every sample's evidence over one cohort locus, against one allele table — the answer to
/// the run arch's "OPEN: CohortObservation's shape" (run arch §8).
pub struct CohortObservation {
    /// The locus span: first position to furthest reach, ≤ `max_cohort_locus_span` (spec §4.1).
    pub region: GenomeRegion,
    /// The distinct alleles the samples showed, each a sequence over the **whole** span —
    /// every member observation projected to the span's width and identical projections
    /// unified (spec §4.2). The reference allele is one of them. Exact-match unification is
    /// sound only because indels were left-aligned upstream
    /// ([`read/left_align.rs:92`](../../../../src/ng/read/left_align.rs)); without that the
    /// same deletion becomes two alleles and one variant's evidence splits in half.
    pub alleles: Vec<Box<[u8]>>,
    /// Indexed by the run's sample order. A sample with no coverage over the span has no
    /// support at all, which is a different fact from reference-only support and stays one.
    pub per_sample: Vec<SampleSupport>,
}

/// One sample's evidence at one locus, expressed against the locus's allele table.
pub struct SampleSupport {
    /// Parallel to `CohortObservation::alleles`. Where two of this sample's own observations
    /// projected onto the same allele, their counts and per-read moments are summed; support
    /// is never merged *across* alleles, because a genotype likelihood needs them apart
    /// (spec §4.2).
    pub per_allele: Vec<SequenceObservation>,
    /// Reads that covered the locus and produced no observation — carried through from the
    /// members, not re-derived.
    pub reads_without_observation: u32,
}
```

**Contract, `build_region`.** The outcome is a pure function of the region's records and `max_cohort_locus_span`
(spec §9): same input, same observations, same failed spans, at any builder count and any division
into regions. Assembly happens only for loci that passed both verdicts — a failed locus costs no
heavy decode (spec §3.2). Cursor requests stay monotonic (spec §11). Errors end the region and name
the sample and span (run spec §9).

```rust
/// Every sample's observations over the ground currently assigned to builders (spec §6.4).
///
/// **One reader per sample for the whole run, forward only.** Builders read this and never a
/// reader: they hold no observations, seek nothing, and cannot mutate it. The organiser is the
/// only writer — it draws the readers forward when it hands out a region, and drops ground once
/// the loci over it are resolved and released.
///
/// **This is the module's dominant memory** (spec §8), which is why regions are short: what it
/// spans is `builders × CohortLocusBuilderRegionsLen`, plus the tail of observations reaching
/// past that span.
pub struct ObservationCache { /* per sample: a forward reader and a window of observations */ }

impl ObservationCache {
    /// Draw every sample forward until `region` is covered, and far enough past it to hold
    /// what a locus starting inside it can reach. Called by the organiser before the region
    /// is handed out; the readers it draws from are the run's, one per sample
    /// ([`run_streaming.md`](run_streaming.md) arch §2).
    fn cover(&mut self, region: GenomeRegion) -> Result<(), RunError>;

    /// Every sample's observations overlapping `span`, for the length of the call — a builder's
    /// only way in, and read-only by construction.
    fn with_observations<R>(&self, span: GenomeRegion, f: impl FnOnce(&[&[SampleLocusObservations]]) -> R) -> R;

    /// Drop everything before `position`. Called once the organiser has released every locus
    /// that could have started there.
    fn evict_before(&mut self, position: GenomePosition);
}

/// Holds the builders' outcomes and the observation cache, resolves the overlaps between
/// neighbouring regions, and releases loci in genome order (spec §6.1, §6.3, §6.4).
pub struct Organiser { /* next_expected, held: BTreeMap<RegionIndex, RegionOutcome>, … */ }

impl Organiser {
    pub fn submit(&mut self, index: RegionIndex, outcome: RegionOutcome);

    /// Everything now releasable, in genome order. **A region is releasable only once its
    /// predecessor has arrived**, because that is what says whether a locus owned earlier
    /// covers this region's first loci (spec §6.3). Failed spans participate in that
    /// resolution and are never released.
    pub fn drain_ready(&mut self) -> impl Iterator<Item = CohortObservation> + '_;

    /// Summed across every region — the number spec §3.3 requires to reach the run summary.
    pub fn failed_locus_count(&self) -> u64;

    /// How many loci — emitted and failed alike — were dropped because an earlier locus already
    /// owned the ground they started on (spec §6.1). **Expected to be zero**, and that is the
    /// point: under `build_region`'s input contract two loci owned by different regions cannot
    /// overlap, so this is what a run would say if that argument failed on real data.
    pub fn displaced_locus_count(&self) -> u64;

    /// Nothing outstanding: every region the run handed out released, every released locus
    /// taken. `regions_handed_out` is how many building regions the run dealt out, their
    /// indexes being `0..regions_handed_out`.
    pub fn is_finished(&self, regions_handed_out: u64) -> bool;

    /// End the run: what it did not emit, or a refusal naming what would have been lost.
    pub fn finish(self, regions_handed_out: u64) -> Result<MergeTally, RunEndedShort>;
}

/// What a finished run has to say about the loci it did not emit. **Both counts leave with the
/// organiser**, because `finish` consumes it: a caller that read them afterwards could not, and
/// one that had to read them first would lose them by taking the obvious order.
pub struct MergeTally {
    pub failed_loci: u64,
    pub displaced_loci: u64,
}
```

**Amended 2026-08-18, after milestone E**, in three places, and each was forced by something
this sketch was written before knowing:

- **`is_finished` and `finish` take how many regions the run handed out.** Without it the
  organiser cannot see a gap at the *tail* of a run — indexes that never submitted, with no
  later index behind them to hold — so a run that lost its last regions was indistinguishable
  from one that finished. All three of E1's review agents found that independently.
- **`finish` returns a `MergeTally`** rather than nothing. `displaced_locus_count` exists to be
  noticed, and a number readable only *before* the consuming call is one the natural call order
  loses.
- **`cache()` is gone, and the organiser does not hold the cache.** The cache is generic over
  its source's error type, because the run's `ObservationSource` and `RunError` do not exist
  yet; making the organiser generic over the same parameters would push that genericity into a
  type with no other use for it. Drawing the readers forward is the driver's
  (`parallel::merge_cohort_in_parallel`), which is also what evicts. **Still owed** when the
  run's own types land.

**Contract, `Organiser`.** Where two loci overlap, the one whose first position is earlier stands
and the other is dropped, whether the earlier one is emitted or failed — one rule, no special case
(spec §6.1). The drain holds at most one region per builder, since a builder cannot start another
until its outcome is taken. Submission order does not affect the released sequence; only the
regions' own order does, which is what spec §9's determinism argument rests on.

---

## 5. Errors (`mod.rs`, extending `RunError`)

A failed locus is **not** an error — it is a counted outcome of a healthy run (spec §3.2, §12).
Two variants join the run's `RunError` (run arch §5):

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RunError {
    // … run arch §5's variants …

    /// A stored observation is wider than its file's recorded reach ceiling: the file is
    /// inconsistent with its own header, and cut safety rested on that ceiling — the stream
    /// can no longer be trusted (spec §3.2).
    #[error("sample {sample}: observation at {region} exceeds the file's recorded reach ceiling")]
    ObservationExceedsReachCeiling { sample: String, region: GenomeRegion },

    /// The drain finished with region indexes never delivered — a builder dropped its
    /// result instead of shipping exactly one `RegionOutcome`, which would silently
    /// truncate the output and the failed-locus total. Production's `MissingChunks`, kept at
    /// release level (`vcf_writer.rs:152-158`; spec §6.3).
    #[error("{count} region result(s) never emitted — a gap stalled the ordered drain")]
    MissingRegionResults { count: usize },
}
```

**Amended 2026-08-18: `MissingRegionResults` became `RunEndedShort`, a three-variant enum**, and
lives in `organise.rs` until `RunError` exists to fold it into.

```rust
#[non_exhaustive]
pub enum RunEndedShort {
    /// Regions the run handed out whose loci never reached the caller: the first index that
    /// stalled the drain, and how many regions from it onwards were never released.
    RegionsNeverReleased { first_stalled: RegionIndex, regions: u64 },
    /// Loci released in order and never taken. Nothing stalled; the caller stopped draining.
    LociNeverDrained { loci: u64 },
    /// Both at once, which one count could not express.
    RegionsNeverReleasedAndLociNeverDrained { first_stalled: RegionIndex, regions: u64, loci: u64 },
}
```

Three things the single struct could not carry. **There are two ways to end short, not one**:
production emits from inside its own submit and has no second step to forget, whereas here
taking the released loci is the caller's own call, so a run can also end with loci released and
never taken. **They can happen together**, which is why the third variant exists. And **one
count told the wrong story**: against a run that had merely stopped draining, the single
message said "a gap stalled the ordered drain" when no gap had. Each variant now names only
what happened, and the two with a stall carry `first_stalled` — the index an operator can map
to a building region and to a builder, which a bare count cannot.

**What stays a panic rather than joining this enum:** a region submitted twice, and one
submitted after it was released. Both are bugs in whoever hands the regions out rather than
facts about the data, and both are caught mid-flight, where the release order is already wrong
and nothing coherent can follow — which is what separates them from a run that ends short,
caught at the end where what was lost can still be named and reported.

---

## 6. Design decisions — decided

- **A locus wider than `max_cohort_locus_span` fails: its span reported in `RegionOutcome`, never built, nothing
  downstream.** The owner's 2026-08-17 ruling — spec §3.2.
- **Two span constants, two types.** `MaxCohortLocusSpan` decides verdicts;
  `ObservationReachCeiling` bounds the cache's reach and nothing else — spec §11.
- **The ceiling is read, never set here**: the generator's `max_record_span` in direct mode, the
  header's recorded value in psp mode, maximum across files, no refusal — spec §5.
- **Regions are short and fixed, about twice `max_cohort_locus_span`; a locus is owned by the region holding
  its first position, and the organiser drops what overlaps an earlier owner.** A segment cannot
  be the unit — it can be a whole chromosome, and its loci would all wait on it — spec §6.1.
- **The observation cache is the organiser's, and it is why builders need no reader of their own.**
  Upstream produces in one forward pass; builders want several places at once. One reader per
  sample fills a shared window, builders read it, and nothing seeks — spec §6.4.
- **`cohort_locus_builder_regions_len` is its own parameter, not twice `max_cohort_locus_span`.**
  Its cost is the cache's span, which has nothing to do with how wide a locus may be — spec §6.4.
- **Builders unify alleles; they do not choose them.** Projection onto the locus span and
  exact-match unification produce the allele table, because that is evidence-shaping and it
  parallelises with the locus. Which alleles are worth calling, and how they are written out, stay
  with the calling steps — spec §4.2, §13.
- **A failed locus resolves overlaps exactly like an emitted one**, so `RegionOutcome` carries
  failed spans and the organiser has one rule rather than two — spec §3.2, §6.1.
- **Closure is uncapped; `max_cohort_locus_span` is a per-locus verdict after it.** Loci stay disjoint, so
  members are moved, not cloned — spec §4.1, §4.2.
- **Verdict order is bound, then variability**, so refused reference-only ground is counted —
  spec §4.3.
- **The keep rule is exact "any non-reference observation in any sample"**, not production's
  `max`-approximated threshold — spec §4.3.
- **Exactly one `RegionOutcome` per region, empty included; gap guarded at release level** —
  spec §6.3.
- **The run's source gains a second way to be asked, when the psp path lands** — revision of run arch §2, recorded above —
  spec §1.3, §4.4.
- **`PspHeader` gains the writing run's reach ceiling** — revision of run arch §4, flagged —
  spec §5, §13.
- **Zero samples refused at caller construction** — spec §7.2.

---

## 7. Reconciliation with existing code

Every row read at the cited line.

| this doc | existing code | how they meet |
|---|---|---|
| saturating reach | `reach` [`cohort_integration.rs:46-48`](../../../../src/var_calling/cohort_integration.rs) | copied, saturating form and all |
| `loci_under`'s closure | `derive_is_kept`'s group walk [`cohort_integration.rs:166-187`](../../../../src/var_calling/cohort_integration.rs) | chaining rule unchanged; ng adds `max_cohort_locus_span` verdict per closed locus |
| the ceiling's source | `max_record_span` [`pileup/generator.rs:93,141`](../../../../src/ng/locus_generation/pileup/generator.rs), ceiling [`:45`](../../../../src/ng/locus_generation/pileup/generator.rs), production default [`walker/mod.rs:67`](../../../../src/pileup/walker/mod.rs) | read and recorded; never re-decided here |
| the physical assumption's prior home | fixture comment "MGS … ≥ ref spans" [`cohort_integration.rs:1664`](../../../../src/var_calling/cohort_integration.rs) | promoted from a test comment to `ObservationReachCeiling` |
| assembling survivors (psp path) | `set_variable_rows` [`sample_reader.rs:789`](../../../../src/var_calling/sample_reader.rs) | the psp reader's `observations_at` follows its one-column-at-a-time shape; bytes are the encoding spec's |
| the ordered drain | `VcfWriter`'s reorder map [`vcf_writer.rs:168-176`](../../../../src/var_calling/vcf_writer.rs), guard [`:152-158`](../../../../src/var_calling/vcf_writer.rs) | carried whole in the callers' skeleton, keyed by region index; `MissingRegionResults` is `MissingChunks` renamed |
| members | `SampleLocusObservations` [`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs) | moved into `CohortObservation::per_sample` unchanged; `reference_bases` ([`:46`](../../../../src/ng/locus_generation/mod.rs)) is what the caller's projection will pad with |
| the ground a region covers | `TypedRegion` [`region_typing/mod.rs:144`](../../../../src/ng/region_typing/mod.rs), `RegionKind` [`:168`](../../../../src/ng/region_typing/mod.rs) | regions are cut inside these; a segment boundary is an unconditional cut and needs no check (spec §5) |
| coordinates | `Position` [`ng/types.rs:34`](../../../../src/ng/types.rs), `GenomeRegion` [`:79`](../../../../src/ng/types.rs) | used as-is; 1-based inclusive, the `+ 1` lives in `GenomeRegion::len` |
| k-stream walking during the fold | `MergedRegionReads` ([`sample_reads.md`](sample_reads.md) §4) | the argmin-over-heads, keys-beside-heads layout, item swapped to light-column rows |

---

## 8. Open items

Genuinely open design questions (each argued in the spec):

- **OPEN: `region_width`** — spec §14 Q1; settled by the
  segment-length measurement.
- **OPEN: how failed loci surface beyond the count** — spec §14 Q4; leaning warn-log capped, a
  BED sidecar if operators need to intersect refused ground.
- **OPEN: calling inside the builder vs after the drain** — spec §14 Q5; decides whether the
  drain's item is `RegionOutcome` as above or a called equivalent carrying the same counts.
- **OPEN: a cohort-scaled keep threshold at the far end** — spec §14 Q3; leaning no.

Implementation-time confirmations (nothing to decide, pin when coding):

- **`CensusWriter::add_generic` should be moved onto `matches_reference`** (§2) rather than keeping
  its inline comparison. It is a change to shipped, tested code that this module does not need in
  order to work, so it is a tidy-up to do alongside the first builder rather than a prerequisite —
  and the census's own tests are what say it is safe.

- The off-by-one edges of the `ceiling − 1` window and of reach at contig ends — write the §15
  boundary tests first, arithmetic second.
  behind the cursor; the fold consumes either.
- `CohortObservation::per_sample`'s index ↔ the run's sample order: pin to the callers' one
  sample table at construction.
- Fold scratch reuse across loci and regions (production's double-buffer pattern) — measure
  before adding anything cleverer.

**Owed after milestone E (2026-08-18):**

- **The organiser should hold the cache**, as §4 originally sketched. It cannot while the cache
  is generic over its source's error type; when `run_streaming.md` arch §2's `ObservationSource`
  and §5's `RunError` land, that genericity goes and `cache()` can come back.
- **`RunEndedShort` should fold into `RunError`** with `ObservationExceedsReachCeiling`, at the
  same moment.
- **`ObservationCache::cover` and `evict_before` are `pub(super)`** and should become private to
  `observation_cache.rs` once the organiser is their only caller — which is the same moment
  again.
- **The round's tail is unmeasured**, and spec §6.2 gates the two alternatives to it (an
  `RwLock` over the cache, or windows handed out as owned copies) on measuring it first. What
  *is* measured, on a fabricated 2,000-base stretch with a record every four bases and 8
  threads: at 1,000 samples the whole cached merge takes 111 ms of which the builders are 10 ms,
  so the parallel arrangement can only reach what is left — 82 ms at 16 regions in flight, and
  26 of the 29 ms saved come from evicting once per round rather than once per region, not from
  the builders. **The next performance question this module has is `cover`, not the round.**
- Where the summed `FailedLocusCount` rides at the run's finish — beside the read-filter
  tallies (run spec §8); the reporting surface is the emission step's (spec §13).

---

## 9. Test shape

Unit tests beside each file; spec §15's table maps each production ancestor to its ng test, and
its "new tests the ruling requires" list is the failed-locus suite (`close.rs`: verdicts, counts,
scheduling-invariance; `build.rs`: suppression is whole-locus, neighbours untouched). The
regression anchor is partition-invariance — one builder over whole segments equals many builders
over any division into regions, observations and failed spans both — with the k = 1 multi-region run
first among the integration tests (spec §7.2), and the eager whole-region build kept as the
two-phase path's byte-identity oracle (spec §15). Run-level oracles are inherited unchanged from
[`run_streaming.md`](run_streaming.md) spec §12.
