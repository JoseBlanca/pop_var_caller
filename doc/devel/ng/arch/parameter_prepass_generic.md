# ng step 4 — the generic prior estimates: types & interfaces

*Status: architecture draft (2026-08-05), companion to the spec
[`../spec/parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md) (the design and
its rationale) and to its shared framing
[`../spec/parameter_prepass.md`](../spec/parameter_prepass.md). Shared arch docs:
[`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary + step traits) and
[`module_layout.md`](module_layout.md) (the `src/ng/` tree). `ng_step_interfaces.md` line 343
gives `SampleSummarizer` a signature the spec replaces; §2.4 and §3 below are the replacement. Naming
follows [`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): domain nouns
for types, verbs for functions, newtypes for domain scalars. Signatures are illustrative; the
**contract** is the deliverable. Every "why" is the spec's — this doc does not re-argue them.*

**Scope: the SNP/indel path only** — the two generic histograms and the four numbers fitted from
them. The STR histogram, the two censuses and the cohort gather are separate units; the fitting
machinery in §4 is shared with the STR path and is specified here because this path is its first
consumer.

## Module home

`src/ng/parameter_estimation/` — the module [`module_layout.md`](module_layout.md) reserves for
step 4, named for what it owns rather than for when it runs. Two sub-units, split so that the
shaping of data and the mathematics on it never live in one file:

```
src/ng/parameter_estimation/
├── mod.rs                    – step 4's surface; routes one locus by LocusKind (§3)
├── fitting/                  – the mathematics. Knows nothing about markers, loci or windows
│   ├── mod.rs                – NoiseModel, the profile scan (§4.2)
│   └── mixture_weights.rs    – the concave climb (§4.1), shared by two fits
└── generic/
    ├── mod.rs                – the two accumulators and what each fits (§3, §5)
    ├── depth_and_alt_reads.rs – data shaping: one locus → one cell key (§2.3)
    ├── depth_bins.rs         – the binning rule: which depths share a bin (§2.2)
    ├── histogram.rs          – the cell table, and the fold over it (§2.2)
    └── runs.rs               – the inbreeding coefficient: a two-state HMM over windows (§5.3)
```

**Why the binning rule is its own file** (added 2026-08-06, when Milestone A4 built it; an
earlier draft of this table put it in `histogram.rs`). It fits neither side of `generic/`'s
data-shaping-versus-mathematics split: it is a fixed rule, not a shaping of this sample's data
and not a fit. §2.2 names three consumers for it inside `generic/` — the cell table, the
subsampling cap of §2.3, and the memory arithmetic of spec §9 — and it has no consumer outside
step 4. Keeping it separate is also what lets `DepthBinEdges` refuse `PartialEq` and `Clone`
without those refusals reading as arbitrary restrictions on the table (§2.2).

`fitting/` is a folder rather than a file because it is the one place with a genuine swappable
seam — one trait, one implementation here and a second on the STR path (spec
[`parameter_prepass.md`](../spec/parameter_prepass.md) §3.2). Everything under `generic/` has no
competing implementation and carries no trait ceremony.

## 1. What this unit consumes and produces

**Input:** the locus stream, `Iterator<Item = Result<SampleLocusObservations, LocusGenerationError>>`
([`locus_generation/mod.rs:701`](../../../../src/ng/locus_generation/mod.rs)). This unit **borrows
each item and passes it on untouched** — it never takes ownership, so it composes into an existing
consumer loop rather than replacing one.

**Output:** four numbers per sample, each carrying where it came from and how much data stood
behind it (spec [`parameter_prepass.md`](../spec/parameter_prepass.md) §6).

## 1.1 The entry point

**Two ways in, because the step has two callers.** A run that walks the loci for step 4 alone wants
one function; a run that folds step 4 into an existing consumer loop wants the accumulator, so the
loci can be passed on untouched.

```rust
/// Everything the fits need that is not in the loci themselves.
pub struct GenericEstimationConfig {
    /// The sample's name as the alignment files declared it — only ever used to name the
    /// sample in an error message and in the emitted summary. **A `String`, not a newtype
    /// and not `SampleIdentity`**: ng has no `SampleId`, and the one identity type it does
    /// have (`read/input/mod.rs:424`) compares files by pointer, which answers "are these
    /// the same open sample?" rather than "what is this sample called".
    pub sample_name: String,
    pub read_groups: Vec<ReadGroupId>,
    pub ploidy: Arc<dyn PloidyMap>,
    pub inbreeding: InbreedingMode,
    pub edges: Arc<DepthBinEdges>,
    /// The read-admission policy the loci were produced under, recorded so that every
    /// emitted `ε` says what population of reads it describes (spec §2).
    pub read_admission: ReadFilterConfig,
}

/// Walk a sample's loci and return its generic parameters. The whole step in one call,
/// for a caller that has nothing else to do with the stream.
pub fn estimate_generic_parameters(
    loci: impl Iterator<Item = Result<SampleLocusObservations, LocusGenerationError>>,
    config: &GenericEstimationConfig,
) -> Result<GenericSampleParameters, ParameterEstimationError>;

/// The same reduction, for a caller that drove the accumulator itself — one per region
/// shard, merged (§3). This is the half that does no I/O.
impl GenericAccumulators {
    pub fn estimate(&self, config: &GenericEstimationConfig)
        -> Result<GenericSampleParameters, ParameterEstimationError>;
}
```

**Contract.** `estimate_generic_parameters` is `estimate` over one accumulator fed by the stream, so
the two cannot diverge. A `LocusGenerationError` in the stream is fatal and propagates — the loci a
walk failed to produce are missing evidence, not zero evidence, and a rate fitted over a truncated
genome is wrong in a way nothing downstream would announce.

**Where the output goes is not this module's decision, and must be said rather than assumed.**
`GenericSampleParameters` is one of three things a sample contributes; the STR path and the two
censuses produce the others, and step 4's own surface (`parameter_estimation/mod.rs`) assembles them
into the `SampleSummary` that `CohortEstimator` consumes unchanged
([`ng_step_interfaces.md`:351](ng_step_interfaces.md), spec
[`parameter_prepass.md`](../spec/parameter_prepass.md) §7). This document specifies the generic
third and stops there.

## 2. Types

### 2.1 Scalar newtypes

Two groups, and the split matters: the first four are **shared vocabulary** that step 4 happens to
be the first to need, and they go in `src/ng/types.rs` beside `GenomeRegion`, `ReadGroupId` and
`Bp`. The rest are step-4's own and stay in this module.

**Shared — extend `types.rs`.** Every one of them has a consumer outside this step: the likelihood
(step 7) reads the error rate, the genotype prior (step 8) reads the genotype frequencies and the
inbreeding coefficient, and ploidy reaches all of them. They are also the `genotype`/`params`
cluster [`module_layout.md`](module_layout.md) anticipates splitting `types.rs` into — worth
grouping them together now so that split has its seed.

```rust
/// Three constrained rates, one type each — **not one shared `Probability`**. They are
/// all fractions in `[0, 1]`, so a single type would let an inbreeding coefficient be
/// handed to something expecting an error rate and compile. `types.rs` already plans this
/// split: `DomainError`'s doc names `AlleleFreq`, `InbreedingF` and `Theta` as the
/// variants to come, and `DomainError::ErrorRate(f64)` is already there.
///
/// Each follows `MismatchFraction`'s shape (types.rs:243): private field, checked
/// constructor, `.get()`.
pub struct ErrorRate(f64);            // per-base, per read group — DomainError::ErrorRate exists
pub struct GenotypeFrequency(f64);    // how common one genotype is: π_hom_ref/π_het/π_hom_alt
pub struct InbreedingF(f64);          // the name types.rs already anticipates

/// `try_new` is the **boundary** constructor: it rejects a value outside `[0, 1]` rather
/// than coercing it, for values arriving from outside the program — the inbreeding
/// coefficient a user may supply on the command line is the one such path here (spec §6.4).
///
/// The fits construct through the same door and `.expect()`, because a frequency off the
/// simplex or a rate outside `[0, 1]` means our own arithmetic is broken, and there is
/// nothing a caller could do about it. Each such site states the invariant it is relying
/// on, in the `// PANIC-FREE:` style the codebase already uses — which is what
/// `FlatEmission::try_new(EPS).expect("eps in [0, 1]")` already does today.
///
/// **The name is ng's house convention, not std's**, and it is settled: `MismatchFraction`,
/// `ReadBases` and `FlatEmission` all spell their checked constructor `try_new`. std would
/// more often call this `new` and return the `Result` (`CString::new`, `NonZeroU32::new`),
/// reserving `try_` for the fallible twin of an infallible constructor. Consistency inside
/// ng wins; a fourth constrained newtype spelling it differently would only make a reader
/// look for a difference that is not there.
impl InbreedingF {
    pub fn try_new(x: f64) -> Result<Self, DomainError>;
    pub fn get(self) -> f64;
}

/// How many copies of the genome an individual has at this region. An input to every
/// fit, never a global constant — it varies by region within one genome (spec
/// parameter_prepass §3).
///
/// Constrained, unlike the unchecked newtypes elsewhere in `types.rs`: the likelihood
/// divides by `P`, so a zero is a division by zero rather than an odd answer.
pub struct Ploidy(u8);
impl Ploidy {
    pub fn try_new(copies: u8) -> Result<Self, DomainError>;   // rejects 0
    pub fn get(self) -> u8;
}
```

**Step-4's own — stay in `parameter_estimation/`.** Nothing outside this step has a use for them,
and pushing them into the shared vocabulary would put a step-4 model's parameter where every module
sees it — the same argument that keeps `window()` off `SampleLocusObservations` (§6).

```rust
/// Which fixed-width window of the reference a locus falls in:
/// `region.start / INBREEDING_WINDOW_BP` within a contig. Windows never span contigs.
/// Local to this step: the window exists to serve the runs model (§5.3) and nothing else.
pub struct WindowIndex(pub u32);

/// The width of an inbreeding window. Fixed, not a knob: a window size is not a quantity
/// a user is in a position to choose (spec §4).
pub const INBREEDING_WINDOW_BP: Bp = Bp(100_000);

/// The error rates the scan steps through: Phred 10 down to Phred 50 in quarter-Phred
/// steps — 161 rungs, each about 6% apart in probability, spanning 0.1 to 0.00001
/// (spec parameter_prepass §3).
///
/// **Phred appears here and nowhere else.** The rungs themselves are probabilities; the
/// Phred scale is how the ladder is *spaced*, because these rates span orders of magnitude
/// and the distance that matters is a ratio. There is deliberately no `Phred` newtype:
/// `types.rs` already has `LogProb` for the logarithm of a probability, and a second
/// log-scaled type in a different base would make a base mix-up a plausible wrong number
/// instead of a compile error — the very hazard `LogProb` exists to prevent.
pub const ERROR_RATE_LADDER_MIN_PHRED: f32 = 10.0;
pub const ERROR_RATE_LADDER_MAX_PHRED: f32 = 50.0;
pub const ERROR_RATE_LADDER_STEP_PHRED: f32 = 0.25;

pub fn error_rate_ladder() -> Vec<ErrorRate>;   // 161 rungs, ascending in Phred
```

### 2.2 The cell table

**What it is: a tally of what the sites looked like.** Each covered position reduces to two
numbers — how many reads covered it, and how many of those showed something other than the
reference — and the table counts how many positions showed each combination. Ten positions:

| positions | depth | alt reads |
|---|---:|---:|
| 1, 2, 6, 8, 10 | 12 | 0 |
| 3 | 11 | 0 |
| 4 | 12 | 1 |
| 5 | 13 | 0 |
| 7 | 30 | 15 |
| 9 | 11 | 1 |

become six entries: `(12,0)→5`, `(11,0)→1`, `(12,1)→1`, `(13,0)→1`, `(30,15)→1`, `(11,1)→1`. Ten
positions became six numbers; 800 million become a few hundred. **It works because a site enters
the likelihood only through those two numbers** (spec §4) — five sites that all looked like
`(12, 0)` are scored once and multiplied by five, and which sites they were is never asked again.
Fitting the error rate is then 161 candidate rates × 583 cells ≈ **94,000 evaluations**, against
161 × 800 million without the tally.

**Why depths are grouped.** Keeping every exact depth to 100 needs 5,151 entries — `101 × 102 / 2`,
since alt count runs 0 to depth — which at eight bytes a cell is the 330 MB per tomato sample the
spec rejects (spec §9; an earlier figure here said 165 MB, from before the per-cell depth sums the
same section now prices). Depth 100 and 105 say
almost the same thing while depth 2 and 3 say very different things at 3 reads per plant, so bins
are one-per-depth at the bottom and widen going up (spec §4, and the binning rule in
[`../spec/parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md) §4).
**Twenty bins: exact integers to 8, then eleven widening ones to a cap of 124** — the ladder §2.2's
`DepthBinEdges::new` measures and adopts:

```text
bin 0: depth 0     bin  8: depth 8         bin 14: depth 29–36
bin 1: depth 1     bin  9: depth 9–10      bin 15: depth 37–46
bin 2: depth 2     bin 10: depth 11–13     bin 16: depth 47–59
…                  bin 11: depth 14–17     bin 17: depth 60–75
bin 7: depth 7     bin 12: depth 18–22     bin 18: depth 76–97
                   bin 13: depth 23–28     bin 19: depth 98–124
```

That is `Σ (bin's top depth + 1)` = **583 cells** — 45 in the exact bins and 538 above them —
because a bin's row must be as wide as its deepest site's alternative count.

**One more thing has to be in the key, and it buys identifiability rather than precision.** With
the library forgotten, a read shows a non-reference base at the share-weighted rate
`Σ_g w_g·p_j(ε_g)`, and because `p_j` is a straight line in `ε` that equals `p_j(ε̄)` at the single
mean `ε̄ = Σ_g w_g·ε_g`. So a key of total depth and total alternative count sees `ε̄` and **nothing
else** about the individual rates: the likelihood is exactly flat along every combination holding
`ε̄` fixed, and no amount of genome separates the libraries. Keeping which library each alternative
read came from is what breaks that flatness (spec §1).

So the key carries the **alternative reads' attribution** for `k ≤ 4`, and there is no second arm.
Spec §1 carries the measurement: scored as a likelihood the key is exactly unbiased in every world
the harness covers; scored by the average-share plug-in an earlier draft used it is 68% out in
heterozygosity at three reads, on two libraries with the *same* error rate.

```rust
/// The deepest a site is entered at — **the ladder's own top**, not an independent number.
/// `DepthBinEdges::max_depth()` returns it, so the two cannot disagree; an earlier draft
/// set it to 300 while the ladder ended at 124, which left depths 125–300 with no bin and
/// made every cell count in both documents wrong.
///
/// A deeper site is **subsampled** down to it. Losing the extra reads costs nothing: at
/// 124 reads a heterozygote shows about 62 alternative reads and a homozygous-reference
/// site about 0.12, so the genotype is already certain and more depth cannot make it so.
///
/// **Draw the kept alternative count, do not rescale it.** The two are not the same thing,
/// and the difference lands on the `k = 0, 1, 2` cells that carry every bit of the
/// error-rate evidence:
///
/// - *Rescaling and rounding to nearest* is not a subsample at all. A 500-read site with
///   one alternative read becomes `(124, 1)` — an alternative fraction of 1/124 against a
///   true 1/500 — and the bias reverses sign at the depth where a lone alternative read
///   stops surviving the round.
/// - *Rescaling with a stochastic round* — take the floor, add one with probability equal
///   to the fraction — fixes the **mean** and breaks the **spread**. Thinning `k` by a
///   factor `r` this way gives a variance of about `r²·Var[k]`, where a real subsample of
///   the same reads gives `r·Var[k]`. At the 124-from-500 ratio that is four times too
///   narrow, and the fit reads an under-dispersed alternative-count distribution as a
///   cleaner read group than the data supports.
/// - **Subsampling the reads is exact**: keep 124 of the site's reads and count the
///   alternative ones among them, which is a hypergeometric draw and marginally the
///   binomial the likelihood assumes. No rescaling, no correction, nothing to argue about.
///
/// Seed the draw from the locus position rather than a shared stream, so a region-sharded
/// walk and a single-threaded one keep the same reads and `merge` stays exact.
///
/// **Where this happens is `count_whole_site` and `count_by_read_group` (§2.3), not
/// `add_site`.** The histogram records the depth it was handed, so the cap must fire
/// before the pair is built — otherwise the depth in `depth_sums` and the depth the
/// alternative count belongs to are different numbers.
///
/// Subsampled sites are counted in `AccumulationCounts` (§3), because a run where most
/// sites are subsampled is one where the ladder, not the data, is setting the depth.
pub fn max_site_depth(edges: &DepthBinEdges) -> u32;   // = edges.max_depth()

/// How many alternative reads keep their library attribution. Above this the site pools
/// them.
///
/// **A precision choice, not a correctness one**, which is what changed when the scoring
/// rule became a likelihood: the fit is exactly unbiased at any value of this bound, so
/// what it trades is cells against how sharply the two libraries can be told apart. Four
/// is the measured default — at three reads a bound of two is equally unbiased on 28%
/// fewer cells, and neither loses measurable precision against scoring every read against
/// its own library (spec §1, §12.8).
///
/// An earlier draft called four "measured, not argued" and left an `OPEN:` about raising
/// it near Phred 20 and depth 300. Both belonged to the plug-in scoring rule that draft
/// used: the bound was carrying the weight of a broken score, and neither the value nor
/// the open question survives the fix.
pub const MAX_ATTRIBUTED_ALT_READS: u32 = 4;

/// Which depth bin a site's depth falls in.
pub struct DepthBin(pub u16);

/// The tally's key. Both arms are **one entry per site** — this is a tally, and sites with
/// identical keys are counted together. At one read group both collapse to the second, so
/// 1,550 of the 1,707 samples in the archive survey are keyed exactly as they are today.
///
/// **The attributed arm keeps a depth bin and not an exact depth, and that is now
/// measured rather than assumed.** §5.1's closed form takes `n` as the site's depth, so an
/// implementation hands it the bin mean above; the question was whether a rule certified
/// unbinned survives it. It does: under the adopted ladder the fit's asymptotic bias stays
/// at 0.054 rungs and 0.3% across twenty worlds, against exactly zero unbinned
/// ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §4.3). So
/// this arm keeps a `DepthBin`, and the alternative designs the open question named — an
/// exact depth on this arm, an exact ladder much further up, or a score taking the bin's
/// whole depth distribution — are all unnecessary.
///
/// **Two arms, not three.** An earlier draft carried a `WholeBreakdown` arm holding every
/// library's depth at sites of total depth four or less. It existed to keep the
/// average-share plug-in from misbehaving at shallow depth; scoring the key as a
/// likelihood removes the reason for it (spec §1), and with it a `SmallVec` of per-library
/// depths on the hottest path in the accumulator.
pub enum SiteKey {
    /// At most `MAX_ATTRIBUTED_ALT_READS` alternative reads: the depth pools, the
    /// attribution of the alternative reads survives. In read-group order, so the key is
    /// canonical.
    Attributed { depth_bin: DepthBin, alt_by_group: SmallVec<[(ReadGroupId, u8); 2]> },
    /// Too many alternative reads to attribute: everything pools. These are the cells
    /// where the genotype is certain and the error rate is not being estimated from them
    /// anyway.
    Pooled { depth_bin: DepthBin, alt_reads: u32 },
}

/// The binning rule: where one bin ends and the next begins. Built once per run and
/// shared by every histogram in step 4, so two accumulators cannot drift apart and their
/// cells stay comparable — which is why histograms hold it by `Arc` rather than a copy.
///
/// Named for depths and not generic over "any binned quantity" on purpose: a repeat count
/// and a depth are both `u32`, and edges that accept either would let the two be
/// transposed silently (§2.1's argument, applied to the binning rule).
/// Also carries the **row offsets** — where each bin's row starts in a histogram's flat
/// `counts` vector. They are a pure function of the edges, so they are computed once per
/// run and shared rather than stored per window. (The saving is modest — one offset per
/// bin is ~15 entries, so 8,000 copies would be about 1 MB against 37 MB. The reason to
/// share is that one binning rule cannot then disagree with another, which `merge` proves
/// with `Arc::ptr_eq`.)
pub struct DepthBinEdges { /* the edges, the depth → DepthBin lookup, the row offsets */ }

impl DepthBinEdges {
    /// **The default ladder: exact integers to 8, then geometrically widening bins to a cap
    /// of 124, twenty bins in all.** Measured rather than argued, and the measurement
    /// changed the shape: the edges are a **correctness** parameter, not the memory-only
    /// knob an earlier draft assumed
    /// ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §4.3).
    ///
    /// Across twenty worlds — error-rate ratios of 1 and 4, mean depths 3 to 60, even and
    /// 90/10 read splits — this ladder's worst asymptotic bias is 0.054 rungs of the
    /// error-rate ladder and 0.3% in each genotype frequency, on 495 cells at 20 reads
    /// against an exact table's 2,542. Three alternatives are worse by an order of
    /// magnitude and one of them is the obvious economy:
    ///
    /// | ladder | worst `ε̄` | worst `π_hom_alt` |
    /// |---|---:|---:|
    /// | **exact≤8, 20 bins, cap 124 — adopted** | **0.054 rungs** | **0.30%** |
    /// | exact≤8, 16 bins, cap 124 | 0.545 | 1.83% |
    /// | exact≤16, 20 bins, cap 124 | 0.981 | 5.53% |
    /// | exact≤8, 16 bins, cap 300 | 1.038 | 7.96% |
    ///
    /// **The cap competes for bins.** Raising it from 124 to 300 at a fixed bin count
    /// doubles the error-rate bias and quadruples the one in `π_hom_alt` — measured on data
    /// where no site is deeper than 125, so the reach is spent on depths nothing occupies
    /// and paid for out of the depths everything occupies. **The bin count and the cap are
    /// one decision, not two.**
    ///
    /// **Where a ladder can hurt is 10 to 30 reads a site**, not the extremes: at 3 reads
    /// 97 sites in 100 are never binned at all, and at 60 the genotype is certain whatever
    /// the exact depth. Any replacement must be checked in that band, which is where an
    /// ordinary whole-genome run sits — checking it at tomato's 3 reads would pass anything.
    pub fn new() -> Self;
    pub fn bin_for(&self, depth: u32) -> DepthBin;
    pub fn row_start(&self, bin: DepthBin) -> usize;
    pub fn cell_count(&self) -> usize;               // the flat vector's length
    /// The depths this bin holds. A range rather than one endpoint, and `RangeInclusive`
    /// rather than a `(u32, u32)` pair, because the sole consumer is a row width where an
    /// off-by-one silently mis-sizes the table — inclusivity belongs in the type.
    pub fn depth_range(&self, bin: DepthBin) -> RangeInclusive<u32>;
    pub fn bin_count(&self) -> usize;
}

/// A dense, ragged table of site counts: for each depth bin, one counter per alt-read
/// count from 0 to that bin's top depth (`*depth_range(bin).end() + 1` wide). Ragged
/// because alt count cannot exceed depth, and a rectangular table would waste half of
/// itself.
///
/// **Dense and index-ordered on purpose** (§6): iteration order is fixed, so every sum
/// over cells is deterministic without sorting; merging is element-wise addition on
/// integers, so a region-sharded walk is exact whatever order its shards finish in.
/// **Holds the edges by `Arc`**, because it needs them to index its own `counts`: the
/// row offsets live there. A sample has ~8,000 of these — one per 100 kb window — so
/// that is 8,000 pointers to one shared object, 64 kB against 37 MB. Handles rather than
/// copies, so two histograms cannot be binned differently and `merge` can prove it with
/// `Arc::ptr_eq` instead of comparing lengths and hoping.
/// **Generic over its counter width, and the two widths are not interchangeable.** A
/// window accumulates in `u32` and the whole-sample fold produces `u64`, because the two
/// hold quantities that differ by four orders of magnitude:
///
/// | | site count | Σ of exact depths |
/// |---|---:|---:|
/// | one 100 kb window | ≤ 100,000 | ≤ 1.2 × 10⁷ |
/// | folded human genome | 3.1 × 10⁹ | **3.1 × 10¹¹** |
/// | `u32` ceiling | 4.29 × 10⁹ | 4.29 × 10⁹ |
///
/// The depth sum is the one that matters: folded over a human genome it exceeds `u32` by
/// **seventy-two fold**, so a fold that widened only the site counts would silently wrap
/// the very quantity `mean_depth_in_cell` exists to hold — and a wrapped mean depth is a
/// cell scored at the wrong depth, which §2.2's own argument says rails a fit. An earlier
/// draft widened the counts and left the depth sums alone.
pub struct DepthAltHistogram<C: CellCounter = u32> {
    counts: Vec<C>,               // pooled arm: flat, rows located via edges.row_start
    depth_sums: Vec<C>,           // pooled arm: Σ of the exact depths, **per cell**
    fine: BTreeMap<SiteKey, (C, C)>,  // the attributed arm: count and depth sum, as above.
                                      // Sparse, and empty at one read group.
    edges: Arc<DepthBinEdges>,
}

/// What a cell may be counted in. Two implementors, `u32` and `u64`, and no third — the
/// trait exists to make the widening at the fold explicit in the type rather than to
/// invite arithmetic over "any integer".
pub trait CellCounter: Copy + Into<u64> { /* add, zero */ }

impl<C: CellCounter> DepthAltHistogram<C> {
    pub fn new(edges: Arc<DepthBinEdges>) -> Self;

    /// Add one site to the tally: find its depth bin, bump that cell's counter, and add
    /// the depth to the bin's running sum. The bin is derived here, from the **exact**
    /// depth `DepthAndAltReads` carries — which is also what `mean_depth_in_cell` needs.
    ///
    /// Takes the pair whole rather than two `u32`s, so that a depth and an alt-read count
    /// cannot be handed over transposed.
    ///
    /// **The depth recorded is the one the site was entered at**, which is the capped one
    /// where the cap fired (§2.3). Recording the pre-cap depth would put a mean depth in
    /// the cell that no site in it actually had, and reintroduce the `n < k` divergence
    /// below by a different route.
    pub fn add_site(&mut self, counts: DepthAndAltReads);
    pub fn merge(&mut self, other: &Self);           // element-wise; panics unless the edges are the same object

    /// Every cell, materialised once. **A `Vec` rather than an iterator, because the
    /// profile scan re-walks these 161 times** — once per rung of the error-rate ladder —
    /// and an iterator would either re-derive the attributed arm's keys on every pass or
    /// need a `Clone` bound that its `BTreeMap` walk cannot cheaply satisfy. A few hundred
    /// to a few thousand cells is kilobytes; building it once and lending a slice is
    /// simpler than either alternative (§4.2).
    ///
    /// The ploidy travels with each cell because one error rate is fitted per read group
    /// **across** ploidies (§6), so a single scan sees cells of more than one.
    pub fn cells(&self, ploidy: Ploidy) -> Vec<(SiteKey, Ploidy, u64)>;

    /// How many **loci** entered — not how many reference positions they covered. A
    /// generic locus widened to an indel's reference span
    /// ([`open_record.rs:573`](../../../../src/ng/locus_generation/pileup/open_record.rs))
    /// is one entry here and several covered positions.
    pub fn total_loci(&self) -> u64;

    /// How many reference **positions** those loci covered — `Σ region.len()`, accumulated
    /// alongside. This is what §5.3 weights window posteriors by, so that `F` is a
    /// fraction of the analysable genome rather than of the locus list, and so that
    /// windows dense in indels are not under-weighted.
    pub fn total_covered_positions(&self) -> u64;

    /// The depth a **cell's** sites are scored at: the mean of the exact depths that
    /// landed in that cell. On the histogram, not on the edges — it is a property of this
    /// sample's data, not of the binning rule.
    ///
    /// **Per cell, not per bin, and that is a correctness requirement rather than a
    /// refinement.** Binning the depth while keeping the alternative count exact means a
    /// bin covering depths 100–124 holds cells up to `alt = 124`. Scored at the *bin's*
    /// mean — necessarily below its top — a homozygous-non-reference site at depth 124
    /// gets `n − k = −12`, and its term `(ε/3)^(n−k)` grows without limit as `ε` falls.
    /// Per cell the problem cannot arise: every site in cell `(bin, k)` has `n ≥ k`, so
    /// the cell mean does too, and the measurement confirms it — **zero truth mass on a
    /// cell scored below its own alternative count, on every ladder and at every depth**
    /// ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §4.2).
    ///
    /// **What the bin mean actually costs is worse to live with than what an earlier draft
    /// of this paragraph claimed.** That draft said the objective becomes unbounded and
    /// every fit rails to the ladder's floor. It does not: the 0.3% of sites whose term
    /// grows as `ε` falls are outweighed by the sites showing one or two alternative reads,
    /// whose terms fall faster, so the objective stays bounded. The fit lands **5.2 rungs
    /// below the true error rate and 29% below the true `π_hom_alt`** and reports nothing
    /// (research note §4.5). A railed fit announces itself — `argmax_at_ladder_end` (§4.2)
    /// exists for that. A rate a quarter low announces nothing, so the rail flag is not the
    /// protection here; the per-cell mean is.
    ///
    /// The cost is a second counter per cell rather than per bin — eight bytes a cell in
    /// a window, four for the site count and four for the depth sum. That is ~4.7 kB per
    /// window, ~37 MB per tomato sample and ~145 MB per human one, which is what spec
    /// §9's table prices.
    ///
    /// **Defined for the attributed arm too, and by the same rule.** Those cells carry a
    /// `DepthBin` as well, so they need a depth to be scored at; the `fine` map holds the
    /// depth sum alongside the count for exactly this. The divergence the paragraph above
    /// describes cannot reach them — an attributed cell has at most four alternative reads
    /// and a depth above the pooling threshold, so `n ≥ k` with room to spare — but a cell
    /// scored at a made-up depth is wrong whether or not it diverges.
    pub fn mean_depth_in_cell(&self, cell: &SiteKey) -> f64;
}
```

Which is used like this — the whole loop the two types exist for:

```rust
let edges = Arc::new(DepthBinEdges::default());            // once per run, shared by every worker
let mut histogram = DepthAltHistogram::new(edges.clone()); // one atomic increment, here only

for locus in locus_stream {                                // one per covered position
    let counts = depth_and_alt_reads::count_whole_site(&locus);  // DepthAndAltReads { 12, 1 }
    histogram.add_site(counts);                     // a deref to bin it; no atomic
}

// Scoring one candidate error rate: a pass over cells, not over the genome. The cells are
// materialised once and re-walked at every rung of the ladder (§4.2).
let cells = histogram.cells(ploidy);
let mut score = 0.0;
for (cell, cell_ploidy, sites) in &cells {
    let depth = histogram.mean_depth_in_cell(cell);   // ≥ this cell's alt count by construction
    score += *sites as f64 * ln_likelihood(depth, cell, *cell_ploidy, error_rate, frequencies);
}
```

**Contract.** `total_loci()` counts entries; `total_covered_positions()` counts reference bases, and
the two differ because a generic locus can be widened to an indel's reference span. Everything
weighted by "how much genome is here" — the inbreeding denominator above all — takes the second.
`merge` is associative and exact for both arms: the pooled table is element-wise integer addition
and the attributed map is a key-wise sum, so two histograms built from disjoint region shards merge
to the histogram of the union.

### 2.3 One locus → one cell

**What it does:** reduces a locus to the pair the histogram is keyed on. This is the only place
that decides what counts as an alternative read, which is why it is its own file rather than a
method on the locus type — the locus cannot know a model's answer to that question.

```rust
/// One site's evidence, reduced: how many reads covered it, and how many of those showed
/// something other than the reference. `alt_reads ≤ depth` always.
///
/// Named for the two numbers it holds rather than for where they came from, so that a
/// call site reading `-> DepthAndAltReads` needs nothing else to know what it got.
/// **Defined in `histogram.rs`**, because it is the histogram's key; this module produces
/// it, so the dependency points from the shaping to the storage and not back.
pub struct DepthAndAltReads { pub depth: u32, pub alt_reads: u32 }

/// Count one generic locus's reads, at the two grains the two histograms need: split
/// between the read groups that covered the site, or the site taken **whole** — the
/// spec's word for the un-split site (§4). The module name carries the nouns, so these
/// say only which grain they count at.
///
/// **Complete witnesses only.** A read that spanned only part of the locus witnessed
/// neither the reference allele nor an alternative one at the positions it missed, and
/// scoring it as either is the mis-scoring `complete_observations` exists to guard
/// (locus_generation/mod.rs:134). Loci one base wide — the overwhelming majority on this
/// path — lose nothing, since a read either covers the base or does not.
///
/// **`reads_without_observation` does not enter the depth** (locus_generation/mod.rs:53).
/// Those reads covered the locus and witnessed nothing, so counting them would assert
/// they showed the reference — and they did not show anything. The depth here is reads
/// whose evidence is present, which is what every likelihood in this document conditions
/// on. They are counted and emitted instead: a locus where many reads say nothing is one
/// where the depth is not what it appears, and that is a fact about the mapping rather
/// than about the genotype. *This also settles what was an open item, and it had to go
/// this way:* the field is a bare scalar with no read-group attribution, so including it
/// would make §9's cell-for-cell equality between the two histograms fail by construction.
///
/// **`reads_discarded_by_cap` does not skip the locus** (locus_generation/mod.rs:56). An
/// earlier draft skipped any locus whose reads an upstream cap had subsampled. That was
/// depth-dependent selection of exactly the kind spec `parameter_prepass.md` §2.1 exists to
/// remove: ng's generic pileup caps indel-bearing columns at 250 reads against 8,000 for
/// the rest ([walker/mod.rs:89](../../../../src/pileup/walker/mod.rs)), so at 300× every
/// indel-bearing column would be dropped and at 30× none would — and the coverage-invariance
/// anchor in §9 would be measuring the skip rule rather than the estimator.
///
/// **Entering it at the observed depth would be correct if the cap took a random
/// subsample. It does not.** Every likelihood here conditions on depth, so a site whose
/// reads were subsampled to 250 is a site observed at depth 250 and its alternative count
/// is a binomial draw at that depth — *provided which reads survived is independent of
/// what they show*. The walk truncates instead:
///
/// ```text
/// contributors.truncate(cap);            // genome_walk.rs:870
/// ```
///
/// and `contributors` is built by walking the active-read set, whose own doc states
/// "Iteration order is admission order" ([active_read_set.rs:66](../../../../src/ng/locus_generation/pileup/active_read_set.rs)).
/// Reads are admitted as the walk reaches their alignment start, so **the reads kept are
/// the earliest-starting ones** — which means the capped column sits further into each
/// retained read, nearer its 3' end, where Illumina error rates are highest. The cap does
/// not just lower the depth; it enriches for the noisier observations of that position and
/// biases `ε` upward there.
///
/// **Where it bites, and where it does not.** The cap is 250 on indel-bearing columns and
/// 8,000 on the rest ([walker/mod.rs:83,89](../../../../src/pileup/walker/mod.rs)), so it
/// fires on the indel-bearing columns of any deep sample — essentially all of them at
/// HG002's 300× — and never at tomato's 3 reads a site.
///
/// **Decision: enter the locus anyway, and count it.** Skipping it would be exactly the
/// depth-dependent selection spec `parameter_prepass.md` §2.1 exists to remove, and it
/// would be a *worse* selection than the one it avoided: at 300× every indel-bearing
/// column would go, and at 30× none would. So the compromise is bounded and reported
/// rather than hidden — `AccumulationCounts::loci_with_upstream_subsample` (§3), and a run
/// where that count is a large share of the loci is one whose `ε` describes the ends of
/// reads more than it describes the chemistry.
///
/// **The real fix is one line upstream, and it belongs to locus generation**: make the
/// truncation a subsample seeded from the column's position rather than a prefix of the
/// admission order. That keeps the walk deterministic and reproducible — the property the
/// prefix was presumably chosen for — while making the retained reads independent of where
/// the column falls within them. Until it lands, this counter is the only thing that says
/// how much of a fit rests on it.
pub fn count_by_read_group(
    locus: &SampleLocusObservations,
    out: &mut Vec<(ReadGroupId, DepthAndAltReads)>,   // scratch, cleared then filled
);
pub fn count_whole_site(locus: &SampleLocusObservations) -> DepthAndAltReads;
```

### 2.4 The fitted output

```rust
/// Where a parameter came from. Not an error condition — a parameter fitted on 80,000
/// reads and one borrowed from a neighbour are both usable and the consumer must be able
/// to tell them apart (spec parameter_prepass §6).
pub enum Provenance { FittedHere, Borrowed, Defaulted, Supplied }

/// A fitted number with its warrant. Generic over the rate, so an `Estimate<ErrorRate>`
/// and an `Estimate<InbreedingF>` stay unmixable — the warrant travels without erasing
/// which quantity it is a warrant for. No uncertainty interval: these are priors, and §6
/// of the shared spec says why one is not emitted.
pub struct Estimate<T> {
    pub value: T,
    pub provenance: Provenance,
    pub observations: u64,      // reads for a noise rate, sites for a per-site rate
}

/// Everything the generic path estimates for one sample. Named for what it holds, not
/// "priors": half of it is noise-model terms, and what a caller builds into a prior is
/// the caller's design (spec parameter_prepass §1.2, non-goals).
pub struct GenericSampleParameters {
    /// One per read group — chemistry belongs to the library, not the individual, and
    /// does not vary with ploidy, so there is one of these however many ploidies the
    /// genome holds.
    pub error_rate: BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
    /// The genotype frequencies, **one set per ploidy present**. A genome with a haploid
    /// sex chromosome has two entries; today's runs have one. Heterozygosity and the
    /// homozygous-non-reference rate are read off each set (§2.4).
    pub rates: BTreeMap<Ploidy, Estimate<SampleRates>>,
    /// The fraction of the analysable genome lying in runs of homozygosity (spec §6).
    /// Diploid-only, and `None` when no diploid region exists: `F` above diploidy needs
    /// several identity-by-descent coefficients and is deferred (spec §7).
    pub inbreeding: Option<Estimate<InbreedingF>>,
    /// What the runs model fitted alongside `F`, when it ran (§5.3).
    pub runs_model: Option<RunsModelFit>,
    /// How the coupled error-rate/frequency fit ended (§5.2). One iteration and converged
    /// on a single-library sample, where the alternation is exact.
    pub coupled_fit: FitTermination,
}
```

## 3. The accumulators

**Two objects, differing only in how a site is keyed.** The read-group histogram enters a site
**once per read group that covered it**: a site with 20 reads from one library and 10 from another
becomes two entries, one at depth 20 and one at depth 10. The windowed histogram enters that same
site **once, at depth 30**.

**Which table a parameter comes from follows from what the parameter is a property of.** The error
rate describes the **chemistry** — how the DNA was prepared and sequenced — so two libraries of one
sample can genuinely have different error rates, and it must be estimated per read group.
Heterozygosity describes the **individual**: one genome has one heterozygosity however many
libraries were used to read it, so it must be estimated per sample and must not be split by library
at all.

**Neither table can stand in for the other.** The read-group table is keyed at the grain the error
rate needs, and splitting a site's reads across two entries costs it nothing, because what it counts
is reads. That same split is fatal to the other rate: a site covered by two libraries has become two
shallower entries and its total depth is gone, so heterozygosity cannot be counted there at all.
Both objects are built, and neither is derived from the other (spec §1, §5.1).

**Both are `BTreeMap` rather than `HashMap`.** A `HashMap` iterates in an
arbitrary order that varies between runs, which would break two things: the runs model reads windows
**in order** along the chromosome (§5.3), and every fit is a floating-point sum over cells, which is
not associative — so an unstable order would let a fitted rate wobble between runs of identical
data. Sorting on `(ContigId, WindowIndex)` gives genome order for free. The O(log n) lookup is ~13
comparisons at 8,000 windows, against a value that is a 4.7 kB table.

```rust
/// The read-group histogram: one cell table per read group, a site entering once per
/// group that covered it. Kilobytes.
pub struct ReadGroupHistograms { by_group: BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram> }

/// The windowed histogram: one cell table per 100 kb window, a site entering once at its
/// total depth. 37 MB per tomato sample, 145 MB per human one (spec §9).
pub struct WindowedHistograms { by_window: BTreeMap<(ContigId, WindowIndex, Ploidy), DepthAltHistogram> }

/// Which set of genotypes a region's sites are drawn from. **Part of both histograms'
/// keys**, because a haploid region has two genotype classes and a diploid three: pooling
/// them would score a site against the wrong set. One `ε` is still fitted per read group
/// across all of them — chemistry does not know about ploidy — but each cell is scored at
/// its own (§5.1), and the genotype frequencies come out one set per ploidy (§2.4).
///
/// **What this module needs is only a lookup**, `region → Ploidy`. Where that comes from —
/// a command-line flag, a BED, a per-contig default — is not settled and is not this
/// module's business; a run with no such input hands over a constant and every key takes
/// the same value, which is today's behaviour exactly.
pub trait PloidyMap { fn ploidy_at(&self, region: GenomeRegion) -> Ploidy; }

/// Whether this run fits the inbreeding coefficient or was handed one.
///
/// **It changes what is accumulated, which is why it is a constructor argument and not a
/// flag read at the end.** `Fitted` keys the windowed histogram by window, as §5.3's runs
/// model needs. `Supplied` drops the window key — the object itself stays, because §5's
/// per-site rates still need whole sites — collapsing it from ~37 MB to a few kB per
/// sample (spec §6.4). The supplied value is emitted with `Provenance::Supplied`.
pub enum InbreedingMode { Fitted, Supplied(InbreedingF) }

/// Step 4's generic front door: holds both, and adds one locus to both.
///
/// **Holds the run's binning rule and hands a handle to every histogram it creates**, so
/// that all of a worker's tables — and every other worker's — are binned by one shared
/// object. Cloning the handle is an atomic increment, paid once per histogram at
/// creation; binning a site is a dereference and costs no atomic at all.
pub struct GenericAccumulators {
    edges: Arc<DepthBinEdges>,
    /* the two histogram collections above */
}

impl GenericAccumulators {
    /// One per **region shard**. Every shard must be handed the same edges, or their
    /// cells mean different things and `merge` is meaningless — which is the whole reason
    /// the rule is shared rather than rebuilt per worker.
    /// `read_groups` is the sample's own set: at one read group no site ever takes the
    /// attributed arm, since a lone library's attribution carries nothing the pooled key
    /// does not — so 1,550 of 1,707 samples pay nothing for the multi-library machinery.
    pub fn new(
        edges: Arc<DepthBinEdges>,
        read_groups: &[ReadGroupId],
        ploidy: Arc<dyn PloidyMap>,
        inbreeding: InbreedingMode,
    ) -> Self;

    /// Add one locus to both histograms — once per read group that covered it in the
    /// first, once at its total depth in the second. **Borrows** — the caller keeps the
    /// locus and passes it on unchanged. A locus whose `kind` is not generic is ignored.
    /// The locus's ploidy is resolved here, from its region, and joins both keys.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations);

    /// Combine a shard's accumulators into these. Associative and exact (§2.2).
    pub fn merge(&mut self, other: Self);

    /// Sum the windows into one whole-sample table, for one ploidy. Free and exact, which
    /// is why no third object is accumulated (spec §1).
    ///
    /// **Both counters widen to `u64` here, and here only.** Per window neither can
    /// overflow `u32` — 100,000 sites and a depth sum under 1.2 × 10⁷ — but the fold over
    /// a human genome reaches 3.1 × 10⁹ sites against a 4.29 × 10⁹ ceiling and
    /// **3.1 × 10¹¹ in the depth sum**, which is seventy-two times over. §2.2's table has
    /// the arithmetic; an earlier draft widened the site counts and left the depth sums at
    /// `u32`, which would have wrapped the mean depth every cell is scored at.
    pub fn whole_sample_histogram(&self, ploidy: Ploidy) -> DepthAltHistogram<u64>;

    /// Everything the accumulator did to a locus other than enter it as it arrived, plus
    /// the overlap counter that must read zero (§2.2, §2.3, and the contract below).
    pub fn adjustments(&self) -> &AccumulationCounts;
}

/// **Everything the accumulator did to a locus other than enter it as it arrived.** Every
/// field is a plain sum, so it merges with the rest of the accumulator; every field is
/// reported, because the alternative is a fit that quietly describes a different
/// population of reads than the caller will see (spec §2).
#[derive(Default)]
pub struct AccumulationCounts {
    /// Loci whose reads a cap upstream had already subsampled — `reads_discarded_by_cap`
    /// non-zero. **Entered anyway, at the depth observed** (§2.3): a site subsampled to
    /// 250 reads is a site observed at depth 250, and skipping it would be exactly the
    /// depth-dependent selection this step exists to remove.
    pub loci_with_upstream_subsample: u64,
    /// Reads that covered a locus and witnessed nothing, summed over loci. Not part of any
    /// depth (§2.3), and a locus with many of them is one whose depth is not what it looks.
    pub reads_without_observation: u64,
    /// Sites this step subsampled down to the ladder's top depth (§2.2). A run where most
    /// sites appear here is one where the ladder, not the data, is setting the depth.
    pub sites_subsampled_to_cap: u64,
    /// Loci that began before the previous locus on their contig ended. **Must be zero** —
    /// see the partition invariant below.
    pub loci_overlapping_previous: u64,
}
```

**Contract.** `add_locus` never allocates after the first locus of a given key. It reads the locus
only through `complete_observations()` and `region`. The two histograms are independent: neither
is derivable from the other once a sample has two read groups, and at one read group they are
equal cell for cell, which is the assertion §9 pins.

**The invariant `add_locus` depends on: generic loci partition the covered positions.** A locus can
span several reference bases — `open_record.rs` widens a record to an indel's reference span
([open_record.rs:573](../../../../src/ng/locus_generation/pileup/open_record.rs),
`footprint_end_exclusive`) — but no two loci may cover the same position. Overlap would enter a
site into the windowed histogram twice, breaking "a site enters once, whole" (spec §4), and
`num_obs_along_locus`'s own doc warns that overlapping generic loci double-count if summed.

**So `add_locus` does not de-duplicate; it counts.** It keeps the previous locus's end per contig —
one comparison per locus — and tallies any locus that starts before it. Counting rather than
asserting, because a debug-only guard compiles out of the release build this repo actually runs, a
trap it has recorded hitting twice
([locus_generation/mod.rs:87-90](../../../../src/ng/locus_generation/mod.rs), the
`num_obs_along_locus` clamp). Counting rather than repairing, because a de-duplication rule would
hide an upstream bug behind a plausible number; a counter that must read zero cannot.

**And it keeps the spans it saw, so that the counter still merges — as a list, not as a pair.** A
"previous locus end" is sequential state, and everything else in this accumulator is a sum of parts:
two shards that split a contig each begin with no previous locus, so neither ever checks the seam
between them, and adding their two counters is not the counter of the union.

The fix is to record, per contig, the **span each shard covered** — its first start and its last
end — and to keep those spans as a list that `merge` concatenates, checking for overlap once at the
end by sorting. An earlier draft collapsed them to a single first-start and last-end pair per
contig and checked the seam during `merge`, which is wrong in a way that only shows up under
scheduling: merging three contiguous shards `[0,100)`, `[200,300)`, `[100,200)` in that order gives
`[0,300)` after the first two, and the third then reads as an overlap that does not exist. **A
counter that must be zero cannot be allowed a false positive**, and the list is bounded by shards ×
contigs — tens to hundreds of pairs — so nothing is saved by collapsing it. With that, merging
shards gives exactly what a single-shard walk would have given, this counter included, whatever
order the shards finish in.

**Concurrency.** One `GenericAccumulators` per region shard, merged at the end — the shape rayon's
`fold`/`reduce` consumes. No trait, no adapter: the accumulator is a plain owned struct and the
"pass through untouched" property is the `&` on `add_locus`, not an iterator wrapper (§6).

## 4. The fitting machinery

Shared with the STR path. It knows nothing about markers, loci or windows — it is given a table of
numbers and returns the values that best explain them.

### 4.1 The piece the two fits share

**What it does:** given, for each cell, how likely each genotype makes that cell, and a weight per
cell, find the genotype frequencies that best explain the whole table. This is the climb the spec
proves cannot get stuck on a false summit (spec [`parameter_prepass.md`](../spec/parameter_prepass.md)
§3.1).

```rust
/// Fit the mixing weights of a finite mixture whose component likelihoods are already
/// known and fixed. Returns one weight per genotype, summing to 1.
///
/// Two consumers, both fitting a free point on the simplex:
///   1. the error-rate scan — called at every rung (§5.1);
///   2. the sample's rates — called once on the whole-sample table (§5.2).
///
/// **Not the runs model.** An earlier draft listed it as a third consumer. Its two states
/// are a *constrained* parameterisation — inside `(1−f−h, h, f)`, outside
/// `((1−f)², 2f(1−f), f²)`, sharing one `f` (spec §6.1) — which is a two-dimensional
/// surface inside the simplex, not a free point on it. A free-simplex maximiser cannot
/// impose the tie, and the concavity guarantee below is a statement about the simplex that
/// does not transfer to a curve inside it. The runs model maximises over `(f, h)` itself
/// (§5.3).
///
/// **Convergence failure is a bug, not a data condition** — the surface is concave, so a
/// climb has no legitimate reason to stall. Asserted in tests, not propagated as a flag
/// no consumer would read (spec parameter_prepass §3.1).
pub fn fit_mixture_weights(
    component_likelihoods: &[&[f64]],   // per cell, per genotype — fixed during the climb
    cell_weights: &[f64],               // site count, or site count × a posterior
) -> SmallVec<[f64; 3]>;
```

### 4.2 The noise model, and the scan over it

```rust
/// What a path assumes can go wrong with a read. The one seam between the two paths:
/// same procedure, two models (spec parameter_prepass §3.2).
pub trait NoiseModel {
    /// The cell type this model's histogram is keyed on — `(depth, alt-count)` here,
    /// a table of repeat-length offsets on the STR path.
    type Cell;
    /// The noise parameters being scanned — one error rate here, three on the STR path.
    type NoiseParams;

    /// How likely each genotype makes this cell, at these noise parameters.
    fn genotype_likelihoods(
        &self, cell: &Self::Cell, noise: &Self::NoiseParams, ploidy: Ploidy, out: &mut Vec<f64>,
    );
}

/// Step through a ladder of noise parameters; at each rung climb to the best genotype
/// frequencies (§4.1) and score; return the best rung with the frequencies found there.
/// A profile likelihood, and a single flat pass — there is no refinement stage, because
/// the spacing is already finer than a caller can feel (spec parameter_prepass §3.1).
///
/// **Takes a slice of cells, not a histogram and not an iterator.** Not a histogram,
/// because an earlier draft took `&DepthAltHistogram` — the generic path's concrete
/// table — which made the "shared with the STR path" claim false: the STR accumulator is
/// repeat-length offsets plus two composition counters, nothing that fits that type. Not
/// an iterator, because the scan re-walks the cells once per rung, and a slice is
/// re-walkable without asking the caller for `Clone` on a `BTreeMap` traversal or
/// re-deriving the attributed arm's keys 161 times. The caller materialises once
/// (`DepthAltHistogram::cells`, §2.2) and lends.
///
/// **Ploidy travels with each cell rather than being one argument**, because one error
/// rate is fitted per read group across every ploidy that group covered (§6): a haploid
/// sex chromosome and the diploid autosomes were prepared by the same chemistry. Each
/// cell is scored against its own genotype set, and the frequencies are climbed **once
/// per ploidy** on that ploidy's cells — a haploid cell has two genotypes to mix and a
/// diploid three, so they cannot share a weight vector.
pub fn fit_by_profile_scan<M: NoiseModel>(
    model: &M,
    cells: &[(M::Cell, Ploidy, u64)],
    ladder: &[M::NoiseParams],
) -> ScanResult<M::NoiseParams>;
```

**Contract.** `fit_by_profile_scan` scores every rung — none is skipped and no early exit is
taken, because nobody has shown the curve has a single hump (spec
[`parameter_prepass.md`](../spec/parameter_prepass.md) §9.3). The generic ladder is 161 rungs at
quarter-Phred spacing ([`parameter_prepass.md`](../spec/parameter_prepass.md) §3).

**It reports whether its answer sat on the ladder's edge, and that is not decoration.** A read group
whose true error rate lies outside Phred 10–50 — a bad run, heavy contamination, or any of the
arithmetic failures a scan can suffer — has its answer silently clamped to an endpoint and emitted
as `FittedHere` with a large observation count. `ScanResult` therefore carries the winning rung, its
frequencies, its score, **and a flag for an endpoint argmax**, which the summary surfaces. One bit,
and it is the only thing standing between a railed fit and a plausible-looking number.

**Static dispatch, deliberately: `<M: NoiseModel>` and not `&dyn NoiseModel`.** The compiler emits
one specialised copy of the scan per noise model, with `M::Cell` substituted, so sharing the
procedure across the two paths costs nothing at run time — no indirect call and no barrier to
inlining in a loop that runs ~75,000 times per fit. The price is a second copy of the code in the
binary, which at two implementations is not worth weighing.

## 5. The four fits

### 5.1 The per-base error rate

`fit_by_profile_scan` over `ReadGroupHistograms`, once per read group, with the generic
`NoiseModel` — a read over a reference copy shows another base with probability `ε`, one over an
alternative copy reverts with `ε/3` (spec §2). Only `ε` is kept from this table; the genotype
frequencies it climbs to at each rung are a means, not an output (spec §3).

**How a cell is scored, and it is a likelihood rather than an approximation.** `p_j(ε)` is the
per-read probability above at `j` alternative copies, `w_g` is library `g`'s share of the sample's
reads, `n` the cell's depth and `k_g` the alternative reads attributed to library `g`:

```text
                                     n!                                                       n−k
ln L(cell | θ)  =  ln  Σ  π_j  ────────────── · Π (w_g·p_j(ε_g))^{k_g} · ( Σ w_g·(1 − p_j(ε_g)) )
                       j        Π k_g! (n−k)!   g                          g
```

A `Pooled` cell is the same expression with the `G` alternative categories collapsed into one, which
leaves a binomial at the share-weighted rate `Σ_g w_g·p_j(ε_g)`.

**Never invent a per-library depth.** The cell has forgotten how the depth split between libraries,
and the tempting repair — give each library `n̂_g = w_g·n` — makes the score stop being a
probability: it does not sum to one over the cell space, and it charges a library for reference
reads it did not have whenever `k_g` exceeds its average share. The expression above sums over the
split instead of guessing it. Spec §1 prices the difference at 68% in heterozygosity at three reads,
on two libraries with the *same* error rate, and it does not shrink as data accumulates. The three
checks in spec §12.8 — sums to one, no negative reference count, exact when the error rates are
equal — are what a reviewer should run against any replacement.

### 5.2 The sample's two rates, and the coupled fit

The same scan, on `whole_sample_histogram()`, keeping the frequencies rather than the noise
parameter. But the two fits are coupled — a higher error rate explains the same alternative reads as
less real variation — and they read two different tables (spec §5.1). So neither can be fitted
without the other, and the loop below is how that is resolved.

```rust
/// Fit the error rates and the sample's genotype frequencies together, alternating
/// between the two tables until neither moves.
///
/// **One call for the whole sample, not one per ploidy.** The error rates come out per
/// read group and span every ploidy that group covered — chemistry does not know about
/// chromosomes — while the genotype frequencies come out one set per ploidy present,
/// because a haploid region has two genotype classes and a diploid three. An earlier
/// signature took a `Ploidy` and returned one rate map, which would have produced a
/// different error rate per ploidy and contradicted §6's decision record.
///
/// **Fallible**, because a sample too thin to fit its genotype frequencies at some ploidy
/// is a real condition and `GenotypeFrequenciesNotFittable` exists for it (§5.4).
pub fn fit_coupled(
    accumulators: &GenericAccumulators,
) -> Result<CoupledFit, ParameterEstimationError>;

pub struct CoupledFit {
    pub error_rate: BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
    pub rates: BTreeMap<Ploidy, Estimate<SampleRates>>,
    pub termination: FitTermination,
}

/// The sample's genotype frequencies at one ploidy: one entry per number of alternative
/// copies, `0..=P`. A vector rather than two named fields, because at `P = 4` there are
/// five and the intermediate dosages have no diploid name (spec §7); the diploid readouts
/// are the accessors below.
pub struct SampleRates {
    pub ploidy: Ploidy,
    pub by_alt_copies: SmallVec<[GenotypeFrequency; 3]>,   // sums to 1
}

impl SampleRates {
    /// Diploid only — how often the individual's two copies differ. `None` above `P = 2`,
    /// where the replacement is gene diversity (spec §7, deferred).
    pub fn observed_heterozygosity(&self) -> Option<GenotypeFrequency>;
    /// How often *every* copy is non-reference: the last entry, at any ploidy.
    pub fn homozygous_non_reference_rate(&self) -> GenotypeFrequency;
}

/// What one scan returned. The score is what makes "best-scoring iterate" a defined
/// comparison in the loop below.
pub struct ScanResult<P> {
    pub noise: P,                          // the winning rung
    pub frequencies: SmallVec<[f64; 3]>,   // the frequencies climbed to at it
    pub log_likelihood: f64,
    pub argmax_at_ladder_end: bool,        // §4.2 — the rail flag
}

/// How the alternation ended. Emitted rather than discarded, because a fit that ran out
/// of iterations is still a number a caller would otherwise consume as if it had settled.
pub struct FitTermination { pub iterations: u32, pub converged: bool }

/// The cap. Generous: a single-library sample converges in one iteration, and only the
/// 157-in-1,707 multi-library samples iterate at all.
pub const MAX_COUPLED_FIT_ITERATIONS: u32 = 20;
```

**The loop, and what "hold fixed" means in each half.** One iteration is:

1. **Each read group's error rate, from the read-group table**, at the genotype frequencies the
   previous iteration produced. This is `fit_by_profile_scan`, so it *does* re-climb the frequencies
   at every rung of the ladder — that is what a profile likelihood is — but the frequencies it
   arrives at are **discarded**. Only the winning rung is kept. What is "held fixed" between
   iterations is therefore the frequencies used in step 2, not anything inside the scan.
2. **The genotype frequencies, from the whole-sample table**, at the rates step 1 just produced.
   One climb per ploidy present (§4.1), because a haploid region has two genotype classes and a
   diploid three.

An earlier version of this paragraph said the frequencies were held fixed while the rates were
fitted, which contradicts the scan it calls — and the contradiction mattered, because the claim
below about single-library samples depends on the scan re-climbing.

**Stop when every read group's winning rung is the same as last iteration's.** The scan returns a
rung index, so "moves by less than one rung" and "does not move" are the same condition, and only
the second is testable — an earlier version stated the first, which a fit oscillating between two
adjacent rungs would never satisfy. Rung stability is the right resolution anyway:
[spec §3](../spec/parameter_prepass.md) argues a difference finer than one rung is smaller than a
caller can feel, and the measurement agrees — worlds that ran 200 iterations without meeting a
movement tolerance of 10⁻¹² were already at the truth to better than a thousandth of a rung
(spec §5.1). A loop that oscillates between two rungs hits the cap and the best-scoring iterate is
kept, which is the behaviour already specified below.

**The alternation converges on the truth, and that is measured rather than argued** (spec §5.1): from
a start at three times the true rates and half the true frequencies, the fixed point is exact in all
25 worlds tried.

**At one read group it ends after one iteration**, which is the reason the decision is low-risk. The
two tables are then the same table (§1), so step 1's scan — which climbs to the best frequencies at
every rung — already returns the joint maximum, and step 2 re-derives what it just found. That is
1,550 of the 1,707 samples in the archive survey.

**Why alternate rather than scan both jointly.** GATK scans jointly, on a 2-D grid of variant prior
crossed with error rate
([`DragstrParametersEstimator.java:196-219`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/tools/dragstr/DragstrParametersEstimator.java)),
and can afford to because its genotype frequencies are **not free**: it scans one variant prior and
derives het and hom-var from it by a fixed het:hom ratio (lines 195-199). That is exactly the tie
spec §11.4 rejects, since it holds only where the inbreeding coefficient is zero. GATK also has one
table per stratum and no read-group split, so it never meets this problem at all. HipSTR, which does
face it, alternates: its `train()` loop runs an expectation step and then updates the genotype
frequencies and the noise model in turn, looping until either the likelihood or **the parameters**
stop moving
([`em_stutter_genotyper.cpp:170-226`](../../../../HipSTR/src/em_stutter_genotyper.cpp),
`max_param_diff = 0.0001`). The parameter-movement half of that test is the rule adopted above.

**Why the cap, given the inner climb needs none.** [Spec §3.1](../spec/parameter_prepass.md) says a
climb that fails to converge is a bug rather than a data condition — but that is the **inner** climb
over the frequencies, which is provably concave. The **outer** alternation has no such proof: in the
standard vocabulary it is *block coordinate ascent*, which converges to a stationary point rather
than provably to the maximum. The literature's specific warning is about **more than two** blocks;
we have exactly two, one of them concave, which is the well-behaved case rather than the cautioned
one. HipSTR caps its loop anyway and returns `false`, and even notes its likelihood is occasionally
non-monotone because of pseudocounts. So: cap, keep the **best-scoring** iterate rather than the
last, and report `FitTermination` — a non-converged fit that arrived silently would be consumed as
though it had settled. This is the same "no silent caps" convention the codebase uses for dropped
reads.

### 5.3 Inbreeding

```rust
/// The fraction of the analysable genome lying in runs of homozygosity, from a two-state
/// hidden Markov model over windows. Sequential, over ~8,000 tomato windows; irrelevant
/// against the walk's cost.
///
/// Takes the **per-read-group** error rates, not one pooled rate: a site with few
/// alternative reads keeps which library each came from (§2.2) and each must be weighed
/// against its own library's rate, by §5.1's expression. Where the site pooled, the
/// share-weighted mean rate is what that expression reduces to, so no separate rule is
/// needed — and those are the cells the error rate is not estimated from anyway.
/// **Fallible**, and the failure is not "too little data": it is a search that never found
/// the second state. Spec §6.5 — a fit whose starting points all guessed the two states
/// too far apart empties the inside state and returns `F` = 0 with every appearance of
/// having converged, so the error path is the only thing between that and a caller.
pub fn fit_inbreeding(
    windows: &WindowedHistograms,
    error_rate: &BTreeMap<ReadGroupId, ErrorRate>,
    ploidy: Ploidy,
) -> Result<(Estimate<InbreedingF>, RunsModelFit), ParameterEstimationError>;

/// The starting points the fit climbs from. **They must disagree about the state
/// separation, not only about `F`** — that is the whole content of this type, and spec §6.5
/// carries the measurement: starts sharing one separation guess return `F` = 0.0000 on a
/// genome 26% covered by runs, converged and silent.
pub struct RunsModelStarts {
    /// The inside state's heterozygote rate as a fraction of the outside one. Default
    /// `[0.05, 1.0/3.0, 0.75]`.
    pub separations: SmallVec<[f64; 3]>,
    /// The implied inside fraction each start begins at. Default `[0.05, 0.5, 0.75]`.
    pub implied_f: SmallVec<[f64; 3]>,
}

/// What the runs model fitted alongside `F`, and how it terminated. Emitted because every
/// one of these is a number someone will want when `F` looks wrong.
pub struct RunsModelFit {
    pub outside_het: f64,             // `Hout`
    pub outside_hom_alt: f64,         // `Aout`
    pub inside_het_floor: f64,        // `h` — absorbs false heterozygotes inside a run
    pub inside_hom_alt: f64,          // `Ain`
    pub enter_run_per_base: f64,      // `tAZ`
    pub leave_run_per_base: f64,      // `tHW`
    pub termination: FitTermination,

    /// **How the search went, and it is the only thing that separates a real `F` = 0 from
    /// a failed search.** Both leave the inside state empty and its frequencies at their
    /// starting values; only the scores say whether a better answer was looked for and
    /// rejected. Sorted by score, best first.
    pub starts_tried: SmallVec<[StartOutcome; 9]>,

    /// The noise floor at this run's window count — what `F` comes back as on a genome
    /// with no runs at all. About 0.01 at 8,000 windows and 0.003 at 31,000 (spec §6.1).
    /// An `F` below it is *nothing detected*, and a consumer that cannot see it cannot
    /// know that.
    pub resolution: f64,

    /// Windows whose posterior landed between 0.01 and 0.99. **Zero at 100 kb**, which is
    /// the measurement saying the chain's transitions changed no window's classification
    /// (spec §6.1). Non-zero is not a fault; it is the chain earning its keep.
    pub undecided_windows: u32,
}

/// One starting point's result.
pub struct StartOutcome {
    pub separation: f64,
    pub implied_f: f64,
    pub inbreeding: f64,
    pub log_likelihood: f64,
}
```

**Contract, and the five places it is easy to get wrong** (argued in spec §6.1 and §6.5):

- **Climb from every start in `RunsModelStarts`, keep the best-scoring, and report them all.** This
  is the first item because it is the one that produces a confident wrong number: nine starts
  spanning the separation return `F` = 0.2634 where five sharing one separation guess return
  `F` = 0.0000 on the same data, converged (spec §6.5). Reject the fit —
  `InbreedingStatesNotSeparated` — when no start left posterior mass on both states, since then the
  answer is *not found* rather than *zero*.
- **Each state carries its own three genotype frequencies, fitted freely** — outside
  `(1 − Hout − Aout, Hout, Aout)`, inside `(1 − h − Ain, h, Ain)`. **Not tied through one allele
  frequency**: that is `bcftools roh`'s form, correct there because its `f` is per site, and with one
  genome-wide `f` it forces `F = 0.57` on an outbred genome (spec §6.1). The identifying constraint
  is the ordering `h << Hout`, applied by relabelling after the fit — nothing in Baum–Welch imposes
  it. `h` is fitted rather than zero: at zero one collapsed paralog inside a run costs the whole
  heterozygote-against-homozygous-reference ratio at that site, about 125 for one alternative read
  of three and past 10³⁵ for fifteen of thirty.
- **The emission is a sum over the window's cells**, never a per-window heterozygote count.
  Dividing by a site count would make a thinly-covered window look autozygous.
- **Both transition rates are fitted, per base**, and converted to a per-window probability. Fixing
  them would set the chain's stationary inside-probability and so assume the shape of the answer —
  but **that stationary probability is not `F`** (spec §6.1); `F` is the coverage-weighted posterior
  occupancy below, and the two differ by 3.5–11% on a finite genome. Fitting the rates makes this
  Baum–Welch: it climbs to a stationary point on a surface that is not concave, so it takes the same
  treatment as §5.2's loop — capped, best-scoring iterate kept, termination reported.
- **The chain covers every window of a contig**, absent ones included as empty, and restarts at each
  contig boundary. The accumulator holds only windows that received a site, so iterating it directly
  would step across an unmappable megabase in one transition and run the end of one chromosome into
  the start of the next.

`F` is the state posteriors weighted by each window's `total_covered_positions()` — reference
positions, not loci, so that a window dense in widened indel loci is not under-weighted (§2.2). That
is the definition, and it is the one that recovers a drawn genome's realised autozygous fraction to
four decimal places (spec §6.5).

**`MIN_WINDOWS_TO_FIT_INBREEDING` is not the same kind of floor as `MIN_SITES_TO_FIT`.** Below a few
thousand windows the estimator's own noise swamps the signal: a genome generated with **no runs at
all** returned `F` averaging 0.23 at 1,200 windows, and 0.84 on one seed of eight (spec §6.1). A
tomato genome is 8,004 windows and a human 31,000, so no real run is near the limit — but a
development fixture or a region-restricted run is, and the number it would produce looks like any
other. Fail rather than emit.

```rust
/// Fewest windows the runs model will accept. Below this its noise floor is of the same
/// size as the answer (spec §6.1, §6.5).
pub const MIN_WINDOWS_TO_FIT_INBREEDING: usize = 3_000;
```

### 5.4 When there is not enough data

**The three parameters fail differently, and only one of them may be guessed.**

```rust
/// Fewest sites a fit will accept before it borrows or fails. *Soft* — §4.1 of the spec
/// shows 6 million read observations pin `ε` to one part in eighty, so this is a floor
/// against noise rather than a precision target.
pub const MIN_SITES_TO_FIT: u64 = 10_000;

/// The error rate used when none can be fitted and none was supplied. *Soft*, and the
/// only defaulted parameter on this path: chemistry varies far less between runs than
/// biology does between samples, so a stated constant is defensible here and is not for
/// inbreeding (below).
pub const DEFAULT_ERROR_RATE: f64 = 0.001;

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ParameterEstimationError {
    /// The sample has too few sites to fit its genotype frequencies at this ploidy.
    /// Not recoverable here: there is no sibling to borrow from — a sample has one
    /// heterozygosity — and no constant worth inventing, since it is the biology.
    #[error("sample {sample}: {sites} sites at ploidy {ploidy} is too few to fit genotype \
             frequencies (need {MIN_SITES_TO_FIT}); supply them or drop the sample")]
    GenotypeFrequenciesNotFittable { sample: String, ploidy: Ploidy, sites: u64 },

    /// `F` was to be fitted and the runs model had too little to run on. **Deliberately
    /// has no default**: inbreeding is the parameter that differs most between an
    /// outcrosser and a selfing landrace, so any constant would be wrong for half the
    /// runs — and the cohort's diversity divides by `1 − F`, so a wrong one is amplified
    /// (spec §6.3). Supply it (`InbreedingMode::Supplied`) or accept the failure.
    #[error("sample {sample}: {windows} usable windows is too few to fit the inbreeding \
             coefficient (need {MIN_WINDOWS_TO_FIT_INBREEDING}); supply one instead")]
    InbreedingNotFittable { sample: String, windows: usize },

    /// Every starting point emptied one of the two states, so no separation was found.
    /// **This is not `F` = 0** — an outcrosser and a failed search leave the same fitted
    /// values, and only the scores across starts tell them apart (spec §6.5). Returning
    /// zero here is the one way this estimator produces a confident wrong number.
    #[error("sample {sample}: the runs model found no second state from any of {starts} \
             starting points — this is a search that failed, not an inbreeding \
             coefficient of zero; widen the state separations or supply F")]
    InbreedingStatesNotSeparated { sample: String, starts: usize },

    /// A domain invariant was violated on the way — a rate off `[0, 1]`, a ploidy of zero.
    #[error(transparent)]
    Domain(#[from] DomainError),
}
```

**The error rate has a ladder of fallbacks, and every rung is visible in the output.**

1. **Fitted here** from that read group's own sites.
2. **Borrowed** — a read group below `MIN_SITES_TO_FIT` takes the mean of the sample's
   other groups' fitted rates. Chemistry differs between libraries, which is the whole reason for
   the read-group grain, so this is a compromise and is marked as one: `Provenance::Borrowed`.
3. **Supplied**, if the run was given one.
4. **Defaulted** to `DEFAULT_ERROR_RATE`, when *every* group is thin and nothing was supplied.

That is what gives `Provenance::Borrowed` and `Defaulted` producers on this path; before this section
they were declared and unreachable. **A consumer that treats all four alike is the failure the
provenance exists to prevent** — a defaulted error rate is a guess, and a caller that cannot tell it
from a rate measured on 80,000 reads will trust both equally.

**Ties on the ladder resolve to the lower error rate**, stated so two implementations cannot differ.
Given §2.2's fixed iteration order the scores are bit-reproducible, so this is about agreement
between implementations rather than between runs.

## 6. Design decisions — decided

- **The accumulator is a plain owned struct, not an iterator adapter — decided.** `add_locus(&…)`
  borrows and the caller keeps the item, which is the pass-through property. An adapter owning the
  iterator would fight the per-shard `create/fill/merge` shape that sharded accumulation needs.
- **Dense ragged cell table, not a hash map — decided.** The key space is small and dense, and a
  fixed iteration order is the whole of the determinism requirement (spec
  [`parameter_prepass.md`](../spec/parameter_prepass.md) §6); a hash map would cost the ordering
  and buy nothing.
- **No trait over the accumulators — decided.** Nothing generic drives them: the walk knows which
  objects it is filling. The five step-4 accumulators also split on lifetime — three reduce when
  their sample ends, the two censuses stay raw to the gather — so one trait would carry a
  `reduce()` two implementors must not have. Shared method *names* are a convention, not a
  contract.
- **The runs model climbs from several starts spanning the state separation — decided.** Starts that
  disagree only about `F` are not a spread: they miss a genome whose real states sit close together
  identically, and "keep the best-scoring" then has nothing better to pick (spec §6.5). The default
  set is three separations × three implied `F`, and it costs seconds on 8,000 windows.
- **`F` is the coverage-weighted posterior occupancy, not the transition rates' ratio — decided.**
  The ratio is a property of the fitted model; `F` is a property of this genome, and only the second
  recovers the realised autozygous fraction (spec §6.1, §6.5).
- **The alt-count decision lives here, not on `SampleLocusObservations` — decided.** Whether a
  partial witness counts as evidence is a modelling choice; the locus type cannot answer it. The
  unambiguous reductions stay on the locus type, which already has them.
- **Windowing is this unit's, not the locus type's — decided.** `INBREEDING_WINDOW_BP` is a step-4
  parameter serving the runs model; no other consumer of locus generation has a use for it, so a
  `window()` method there would push a step-4 concept into shared vocabulary.
- **Complete witnesses only — decided.** See `depth_and_alt_reads`'s contract (§2.3).
- **The coupled fit alternates between the two tables, capped at 20 iterations — decided.** The two
  parameters trade off but are read off different tables, so neither can be fitted first. Rejected:
  GATK's joint 2-D scan, which is only affordable because it ties the genotype frequencies to one
  variant prior by a fixed ratio — the tie spec §11.4 rejects. HipSTR faces the same coupling and
  alternates, stopping on parameter movement, which is the rule adopted (§5.2).
- **A coarsened cell is scored by summing over what the key forgot, never by plugging in an
  average — decided.** The attributed cell's likelihood has a closed form (§5.1), so this costs the
  same as the plug-in an earlier draft used and is exactly unbiased where that one is 68% out in
  heterozygosity at three reads (spec §1). Rejected: `n̂_g = w_g·n`, which is not a probability over
  the cell space and manufactures a two-fold error-rate difference between identical libraries.
- **The whole-breakdown arm is gone, and `MAX_ATTRIBUTED_ALT_READS` is now a precision knob —
  decided.** Both existed to keep that plug-in in bounds; with the score fixed the fit is unbiased
  at any bound, so what four buys is sharpness, not correctness (§2.2, spec §1).
- **The depth ladder is exact integers to 8, twenty bins, cap 124 — decided, and it is a
  correctness choice.** Measured across twenty worlds: this ladder biases the error rate by 0.054
  rungs and each genotype frequency by 0.3%, where sixteen bins at the same cap costs 0.55 rungs and
  1.8% and a cap of 300 at sixteen bins costs 1.04 rungs and 8.0% (§2.2, research note §4.3).
  Rejected: treating the bin count and the cap as independent knobs — a cap buys reach out of the
  same twenty bins that buy resolution, and on data no deeper than 125 the reach is worth nothing.
- **Ploidy is part of both histogram keys — decided.** The likelihood sums over `P + 1` genotypes,
  so cells from a haploid region and a diploid one cannot share a table. Keying it now costs a
  tuple; retro-fitting it would change the accumulator, the merge, the fold and every fit signature
  (spec parameter_prepass §3, which pays a loop bound to avoid exactly that). One `ε` is still
  fitted per read group across all ploidies — chemistry does not know about ploidy — with each cell
  scored at its own.

## 7. Reconciliation with existing code

Every row read before it was written.

| this doc | existing code | action |
|---|---|---|
| the input stream | `SampleLocusObservationsIterator` ([locus_generation/mod.rs:701](../../../../src/ng/locus_generation/mod.rs)) | consume by reference; unchanged |
| the locus | `SampleLocusObservations` ([locus_generation/mod.rs:40](../../../../src/ng/locus_generation/mod.rs)) | reuse as-is |
| per-allele support | `SequenceObservation` ([locus_generation/mod.rs:157](../../../../src/ng/locus_generation/mod.rs)) | reuse; `read_group` is already part of its identity |
| the scorable subset | `complete_observations()` ([locus_generation/mod.rs:134](../../../../src/ng/locus_generation/mod.rs)) | call directly — §2.3's guard |
| marker routing | `LocusKind` ([locus_generation/mod.rs:213](../../../../src/ng/locus_generation/mod.rs)) | match on it in `mod.rs` |
| `WindowIndex` source | `GenomeRegion` ([types.rs:79](../../../../src/ng/types.rs)) | derive from `region.start` |
| `ReadGroupId` | [types.rs:199](../../../../src/ng/types.rs) (`pub struct ReadGroupId(pub u32)`) | reuse as-is |
| read-group identity | `ReadGroup` ([read/input/read_groups.rs:46](../../../../src/ng/read/input/read_groups.rs)) | reuse — resolves a group to sample/library/experiment |
| `ErrorRate`, `GenotypeFrequency`, `InbreedingF`, `Ploidy` | not present as types; `DomainError::ErrorRate` already is ([types.rs:275](../../../../src/ng/types.rs)) | **new in `types.rs`** — shared vocabulary, consumed by steps 7 and 8. Each copies `MismatchFraction`'s checked-constructor shape ([types.rs:242](../../../../src/ng/types.rs)); `DomainError`'s doc already names `InbreedingF` as expected |
| `WindowIndex`, `DepthBin`, `DepthBinEdges`, `DepthAltHistogram`, `CellCounter`, `SiteKey`, `AccumulationCounts`, `CoupledFit`, `RunsModelStarts`, `StartOutcome` | not present | **new, and step-4-local** — no consumer outside this step (§2.1) |
| the sample's name in `GenericEstimationConfig` | `SampleIdentity` ([read/input/mod.rs:424](../../../../src/ng/read/input/mod.rs)) | **not reused** — it compares open files by pointer, answering "the same sample?" rather than "called what?". A `String` from the alignment header, used only in messages and the summary (§1.1) |
| the ladder's log scale | `LogProb` ([types.rs:227](../../../../src/ng/types.rs)) | **not reused, and no `Phred` twin added** — a second log-scaled probability type would make a base mix-up a plausible wrong number rather than a compile error, which is what `LogProb`'s own doc says it exists to prevent (§2.1) |
| `fit_mixture_weights` | `run_em_loop` / `GenotypeEmModel` ([var_calling/posterior_engine.rs:2635 (GenotypeEmModel), 2733 (run_em_loop)](../../../../src/var_calling/posterior_engine.rs)) | **shape to copy, not code to call** — production's EM is tied to its own genotype model, and ng does not depend on production |
| the profile scan | GATK `DragstrParametersEstimator.java` (vendored under `gatk/`) | algorithm copied, not code |
| the runs model | `bcftools roh` ([vcfroh.c:476-499](../../../../bcftools/vcfroh.c), emissions; [452-472](../../../../bcftools/vcfroh.c), distance-scaled transitions) | reference form for §5.3; not a dependency |
| `SampleSummarizer` | [`ng_step_interfaces.md`:343](ng_step_interfaces.md) — `fn summarize(&self, confident: &[ConfidentGenotype]) -> SampleSummary` | **re-specified.** No rough caller exists, so there are no confident genotypes; the input is `GenericAccumulators` (spec parameter_prepass §2.3) |
| site reduction | `SiteCounts::from_record` ([sample_summary/het.rs:146](../../../../src/sample_summary/het.rs)) | **not reusable** — it takes a `PileupRecord`, production's type, and returns `None` for a pure-reference column, which is the evidence this step exists to keep. Shape only |
| the three-genotype score | `observe_site` ([sample_summary/het.rs:266](../../../../src/sample_summary/het.rs)) | shape only — it gates on `min_depth` and classifies, which is the bias the spec removes |

## 8. Open items

- **CLOSED: `DepthBinEdges`' exact edges, and they were never only a memory choice.** The
  adopted ladder is **exact integers to 8, then geometrically widening bins to a cap of 124,
  twenty bins in all** (§2.2). The measurement is
  [research note](../research/parameter_estimator_experiments_2026-08-06.md) §4.3: sixteen
  bins at the same cap costs 0.55 rungs of the error-rate ladder and 1.8% of `π_hom_alt`
  where twenty costs 0.05 and 0.3%, so an earlier draft's premise here — that the edges buy
  memory and not accuracy — was wrong. What remains for implementation time is the memory
  measurement of [`parameter_prepass.md`](../spec/parameter_prepass.md) §10.6, which now
  prices 583 cells a window rather than the 465 an earlier draft of §9's table assumed.
- **CLOSED: what depth binning does to the multi-library score.** Nothing measurable under
  the adopted ladder: the fit's asymptotic bias is 0.054 rungs against exactly zero unbinned
  (research note §4.3), so `SiteKey::Attributed` keeps a `DepthBin` (§2.2) and no arm needs
  an exact depth. Two adjacent facts came out of the same measurement and belong here: the
  *cap* competes for bins, so a cap of 300 at sixteen bins is four times worse in
  `π_hom_alt` than a cap of 124 (§2.2); and the band a ladder can hurt is 10 to 30 reads a
  site, so any replacement checked only at tomato's 3 reads would pass whatever it did.
- `OPEN:` **what a site deeper than the cap costs.** §2.2 specifies subsampling it down by a
  hypergeometric draw, and no harness implements one — the worlds measured above all sit
  below the cap, so the subsampling rule is the one depth mechanism with no measurement
  behind it. It fires only on samples above ~124×, which is HG002 and not tomato.
  **Settled by:** adding the draw to
  [`ng_multilib_key_harness.rs`](../../../../examples/ng_multilib_key_harness.rs) and running
  a mean depth of 300 against the exact answer.

## 9. Test & bench shape

Unit tests beside each file; the acceptance tests are the spec's §12 and
[`parameter_prepass.md`](../spec/parameter_prepass.md) §10.

**Every recovery test in both documents generates its data from the model it then fits**, so a
shared misspecification cancels and passes. Those tests catch gross bugs and cannot catch bias,
which is the failure this step exists to remove. The anchors below are the only evidence in the
design that does not come from the model itself.

**The truth set is the whole-genome one**, `benchmarks/giab/all_bench_regions/` — the HG002 v4.2.1
small-variant VCF and its confident BED. *Not* `benchmarks/ssr_hg002/`, which an earlier draft named:
that is the **tandem-repeat** benchmark, 36,497 records all inside repeat tracts over 6 Mb of
scattered intervals, and region typing routes those tracts to the STR path. It can supply neither a
genome-wide error rate nor a single 100 kb window for the runs model. **Reads are the gap**: the
repo's whole-genome truth has no matching whole-genome alignment, so the coverage sweep below needs
one fetched or the anchor restricted to the regions we hold. That is a data question, not a design
one, and it must be settled before the anchors are called done.

**`F` is the one parameter the restriction cannot rescue.** The HG002 alignments here —
`benchmarks/giab/per_sample/bam/`, 100 randomly selected regions — exist to make development
testable, not to estimate anything, and `F` is the parameter where that distinction bites: at
~1,200 windows the runs model's own noise returns `F` averaging 0.23 on a genome with no runs at all
(spec §6.1). Restricting the anchor to those regions therefore does not weaken it, it voids it. The
`F ≈ 0` check below **needs a whole-genome alignment and has no substitute**; the other three
parameters are fine on a region subset because they are per-read or per-site rates rather than a
fraction of the genome.

- **Model-free values for all four parameters, from GIAB truth.** The error rate is non-reference
  bases over total bases at truth homozygous-reference positions — a count, no model and no fit.
  Heterozygosity is the truth het count over the confident regions' length; the
  homozygous-non-reference rate is the 1/1 count over the same. `F ≈ 0`, because HG002 is not
  consanguineous. The known caveat is that confident regions are the *easy* regions, so the truth
  error rate is a lower bound and the truth heterozygosity is depleted of hard sequence: this bounds
  the fitted values rather than pinning them. A fitted error rate that came out *below* the
  model-free one on easy regions would be an unambiguous bug.
- **Invariance to coverage, which needs no truth at all.** Fit all four parameters on the same HG002
  alignment downsampled to 300×, 30×, 10× and 3×. **Restricted to the confident BED at every rung**,
  or the arms compare different site sets rather than different depths. Same genome, so every one of them must come out flat —
  an error rate is per read and the other three are properties of a genome, and none of them may
  depend on how deeply it was sequenced. Any slope is bias, and its sign names the mechanism. The
  tomato equivalent is free: plot each sample's fitted heterozygosity against its mean depth across
  the cohort, where a biological quantity has no business correlating with library yield.

**Two harnesses carry evidence that does not come from the model either**, and both live in
`examples/` with their measurements written up in
[`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md):

- [`ng_multilib_key_harness.rs`](../../../../examples/ng_multilib_key_harness.rs) — the
  multi-library score's bias, computed **exactly** by weighting each cell with its probability under
  a known truth, so bias is separated from sampling noise by construction. Three algebraic
  assertions gate it, none needing a fit: the rule sums to one over the cell space, no cell is
  charged a negative count of reference reads, and with the libraries' error rates equal the rule
  reproduces the exact per-library likelihood. **Any change to §5.1's expression re-runs these
  first.** It also carries the depth-ladder sweep behind §2.2 (`--only=binning`) and the cost of
  assuming a heterozygote is a half (`--only=balance`).
- [`ng_inbreeding_harness.rs`](../../../../examples/ng_inbreeding_harness.rs) — `F` against a
  simulated genome's realised autozygous fraction, plus the three things §5.3's contract turns on:
  the noise floor by window count, the false-heterozygote robustness, and what happens when the
  starting points do not span the state separation. Its `--only=binning` section refits one genome
  through every candidate ladder and both mean-depth grains, which is what says `F` does not move.

**Each has a control whose answer must be exactly zero, and both earned it.** The multi-library
harness runs every world with the libraries' error rates equal; the binning sweeps run the exact
ladder, where a bin is one depth and the binned code must reproduce the unbinned answer. Two of the
four wrong findings these harnesses have produced were caught by such a control and by nothing else
(research note §6) — including one where the *generator*, not the estimator, put every site in a
window into a single cell above 40 reads a site and the runs model correctly reported one long run.

Neither is a check that the *model* is right — nothing generated from a model can be. They check
that the estimator on top of it throws nothing away.

Five further anchors, all cheap:

- **The two histograms agree on a single-library sample** — fold the windowed one over its windows
  and compare cell for cell against the read-group one. Exact equality, no simulated truth, and it
  runs on the tomato and HG002 cohorts as they stand (spec §12.6).
- **Sharded accumulation is exact** — one sample walked in one region and in many must give
  identical histograms, which the integer merge makes a real assertion rather than a tolerance.
- **`adjustments().loci_overlapping_previous` is zero** on every fixture and on both real cohorts (§3). This
  is the cheapest place the partition invariant is checked against real data, and a non-zero count
  is a bug report against locus generation rather than something for this unit to absorb.
- **The fit recovers known parameters**, at ploidy 2 and 4, from an accumulator filled directly —
  no reads, no reference (spec parameter_prepass §10.1). **Run it at high depth too**: at tomato's
  3 reads every site sits in a one-per-depth bin, so the binning bug `mean_depth_in_cell` fixes
  (§2.2) is invisible below ~100×.
- **No fit rails.** Assert `cell.alt_reads() ≤ mean_depth_in_cell(cell)` at every scored cell, and
  that no scan's argmax is an endpoint of its ladder (§4.2). These are the two ways this estimator
  produces a confident wrong number rather than a visible failure.

No `bench/`: this unit has no competing implementations. Its cost claims are memory, and
[`parameter_prepass.md`](../spec/parameter_prepass.md) §10.6 measures those on real runs — and must be re-run against §2.2's per-cell depth sums, which §9's table does not yet include.
