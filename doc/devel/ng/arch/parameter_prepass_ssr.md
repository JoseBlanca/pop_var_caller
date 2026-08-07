# ng step 4 — the STR stutter parameters: types & interfaces

*Status: architecture draft (2026-08-06), companion to the spec
[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) (the design and its
rationale) and to its shared framing [`../spec/parameter_prepass.md`](../spec/parameter_prepass.md).
The sibling path's architecture, [`parameter_prepass_generic.md`](parameter_prepass_generic.md),
owns the fitting machinery this unit consumes and the scalar newtypes it shares; read its §4 before
this document's §3. Shared arch docs: [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary +
step traits) and [`module_layout.md`](module_layout.md) (the `src/ng/` tree). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): domain nouns for types,
verbs for functions, newtypes for domain scalars. Signatures are illustrative; the **contract** is
the deliverable. Every "why" is the spec's — this doc does not re-argue them.*

**Scope: the STR path only** — one accumulator, keyed by `(read group, motif period, repeat count)`,
and the four numbers fitted from it. The generic histograms, the two censuses and the cohort gather
are separate units.

**The types in §2 have been through the measurements**, which is why this document exists after them
rather than before: the entry key, the offset origin and the end-bucket scoring were all moved by
[`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§6, and writing them earlier would have meant unwinding them.

## Module home

`src/ng/parameter_estimation/ssr/` — beside `generic/` under the module
[`module_layout.md`](module_layout.md) reserves for step 4. **STR-ness is a property of an
implementation, not a top-level division** (`module_layout.md` principle 2), which is why this is a
sibling folder rather than a separate tree.

```
src/ng/parameter_estimation/
├── mod.rs                    – step 4's surface; routes one locus by LocusKind
├── fitting/                  – the mathematics, shared. NoiseModel, the concave climb,
│                               and the two search drivers (§3)
├── generic/                  – the SNP/indel path
└── ssr/
    ├── mod.rs                – the accumulator, and what is fitted from it (§4)
    ├── locus_offsets.rs      – data shaping: one STR locus → one entry key (§2.3)
    ├── stratum_table.rs      – the sparse table of locus shapes, per stratum (§2.2)
    └── slippage.rs           – the NoiseModel implementation and its search (§3, §4.1)
```

A folder rather than a file for the reason `generic/` is one: the shaping of data and the
mathematics on it never live together. There is no trait over the accumulator — nothing generic
drives it.

## 1. What this unit consumes and produces

**Input:** the same locus stream the generic path reads,
`Iterator<Item = Result<SampleLocusObservations, LocusGenerationError>>`
([`locus_generation/mod.rs:701`](../../../../src/ng/locus_generation/mod.rs)). This unit **borrows
each item and passes it on untouched**, and **ignores every locus whose `kind` is not
`LocusKind::Ssr`** ([`locus_generation/mod.rs:217`](../../../../src/ng/locus_generation/mod.rs)).

**Output:** four numbers per `(read group, stratum)` — a slippage rate, a direction split, a
distance decay and a substitution rate — each carrying where it came from and how many loci stood
behind it, plus **one summary across strata** that a person reads (§4.3).

```rust
/// Everything the STR fits need that is not in the loci themselves.
pub struct SsrEstimationConfig {
    /// As the alignment files declared it — only ever used to name the sample in an
    /// error message and in the emitted summary. A `String` for the reason
    /// `GenericEstimationConfig` gives (generic arch §1.1).
    pub sample_name: String,
    pub read_groups: Vec<ReadGroupId>,
    pub ploidy: Arc<dyn PloidyMap>,
    /// The read-admission policy the loci were produced under, recorded so that every
    /// emitted rate says what population of reads it describes (spec §4.1).
    pub read_admission: ReadFilterConfig,
}

/// Walk a sample's loci and return its STR parameters. The whole unit in one call, for a
/// caller with nothing else to do with the stream.
pub fn estimate_ssr_parameters(
    loci: impl Iterator<Item = Result<SampleLocusObservations, LocusGenerationError>>,
    config: &SsrEstimationConfig,
) -> Result<SsrSampleParameters, SsrEstimationError>;

/// The same reduction, for a caller that drove the accumulator itself — one per region
/// shard, merged (§4). The half that does no I/O.
impl SsrAccumulators {
    pub fn estimate(&self, config: &SsrEstimationConfig)
        -> Result<SsrSampleParameters, ParameterEstimationError>;
}
```

**Contract.** `estimate_ssr_parameters` is `estimate` over one accumulator fed by the stream, so the
two cannot diverge. A `LocusGenerationError` is fatal and propagates: loci a walk failed to produce
are missing evidence, not zero evidence, and a stratum fitted over a truncated set of loci is wrong
in a way nothing downstream would announce. `SsrSampleParameters` is one of three things a sample
contributes; step 4's own surface assembles them into the `SampleSummary` that `CohortEstimator`
consumes unchanged ([`ng_step_interfaces.md`:351](ng_step_interfaces.md)).

## 2. Types

### 2.1 Scalar newtypes

**Shared — extend `types.rs`.** `SsrPeriod` is STR domain vocabulary with consumers in steps 6 and 7,
and [`module_layout.md`](module_layout.md) already names it as expected there beside `Bp` and
`LogProb`. The rest are step-4's own.

```rust
/// A repeat unit's length in bases, 1..=MAX_MOTIF_LEN. **The stratification axis's first
/// half** (spec §4). Constrained, because a period of zero would divide by zero when a
/// tract length is converted to a repeat count.
///
/// `Motif::period()` (types.rs:369) returns a `usize` today; this is the checked domain
/// type it should have returned, and `Motif` gains `ssr_period()` beside it rather than
/// changing the existing accessor's type under its callers.
pub struct SsrPeriod(u8);
impl SsrPeriod {
    pub fn try_new(bases: u8) -> Result<Self, DomainError>;   // rejects 0 and > MAX_MOTIF_LEN
    pub fn get(self) -> u8;
}

/// How many whole motif copies a tract holds — **the reference tract's**, not the
/// sample's. The stratification axis's second half, and a pure function of the reference
/// (spec §4.1), so every sample strata-fies identically and a cohort can compare strata.
pub struct RepeatCount(pub u32);

/// One group of loci that gets its own fitted parameters (spec §4). Ordered by
/// `(period, repeats)` so that §4.2's monotonicity walk visits neighbours in order.
pub struct Stratum { pub period: SsrPeriod, pub repeats: RepeatCount }
```

**Step-4's own — stay in `parameter_estimation/ssr/`.**

```rust
/// How far a read's tract sits from the **reference** tract length, in whole motif
/// copies. Negative is shorter. The origin is the reference and not each locus's modal
/// observed length: measured, and the modal origin returns a slippage rate 50% to 408%
/// high with the direction asymmetry destroyed (spec §4.1).
pub struct WholeRepeatOffset(pub i8);

/// The offsets an entry records, either side of the origin — `±OFFSET_HALF_RANGE`, with
/// the end buckets absorbing everything beyond.
///
/// **Measured to matter far less than it looks.** Against the reference origin the end
/// buckets absorb whole alleles and not only far slips, so a narrow range looks dangerous;
/// with the ends scored by their marginal (§3), a range of **±1** against alleles reaching
/// ±3 still returns the slippage rate to within 0.05% and both shares to within 0.002
/// (spec §4.1). What a narrow range costs is the heterozygosity that falls out of the
/// fitted genotype frequencies — 1.5% at ±1 — which this path does not emit.
///
/// **Four is comfortable on real data**: the saturating end buckets take **0.89% of reads**
/// across GRCh38's typed tracts and **0.14%** across tomato's 138 million (research note
/// §6.8). Nothing is piling up against the ends.
///
/// **The width that *is* load-bearing is `ALLELE_OFFSET_LIMIT`, and it is not this one.**
pub const OFFSET_HALF_RANGE: i8 = 4;

/// How far from the reference tract length the fit may place an allele — the support of
/// the genotype frequencies of §2.4, and **a separate width from the recorded range,
/// wider than it**.
///
/// It is what lets the marginal rule attribute an end bucket to a distant allele instead
/// of to a far slip, which is why it and not `OFFSET_HALF_RANGE` decides the answer. Too
/// narrow and the fit charges a distant allele to slippage; too wide and a thin stratum
/// fits frequencies for lengths no locus carries — and the count grows as
/// `A(A+1)/2` (§2.4).
///
/// **It is clipped at the low end, so `A` is a per-stratum quantity and not a constant.**
/// An allele cannot be shorter than nothing: a stratum at 4 repeats reaches only −4, so it
/// has 11 lengths and 66 genotypes where a stratum at 6 or more has the full 13 and 91.
/// The clip is not a special case to remember at the fit — it is what the support *is* —
/// but it is the reason `fit_mixture_weights` is handed a length rather than assuming one.
///
/// **Six comes from the measured distribution, and it is a threshold to clear rather than a
/// number to tune.** A locus's modal observed length against the reference, on HG002 at
/// 300×: 88.9% sit exactly at the reference, ±4 holds 99%, ±12 holds 99.9%, ±19 holds
/// 99.99%; tomato is tighter, 95.7% at zero and ±1 holding 99% (research note §6.8). Six
/// covers all but about one human locus in 200.
///
/// **What a locus outside the support costs is nothing, and then everything** (research
/// note §6.4.1): leaving 2.5% of loci out costs 0.1% of the slippage rate, 7.9% costs 2.5%,
/// and **19.3% costs +499%, with the direction asymmetry destroyed** — 0.17 becoming 0.47,
/// the same collapse the modal origin produced. Six sits a fifth of the way to the row that
/// is already free, so its exact value does not matter; what matters is never choosing one
/// narrow enough to leave a tenth of the loci out.
pub const ALLELE_OFFSET_LIMIT: i8 = 6;

/// Which bucket an offset falls in: `offset.clamp(±OFFSET_HALF_RANGE) + OFFSET_HALF_RANGE`,
/// so `0 ..= 2·OFFSET_HALF_RANGE`. The two end buckets are saturating.
pub struct OffsetBucket(pub u8);
pub const OFFSET_BUCKETS: usize = (2 * OFFSET_HALF_RANGE + 1) as usize;

/// Reads entered from one locus. A deeper locus is entered from a **random subsample** of
/// its reads down to this, seeded from the locus's position so a region-sharded walk and a
/// single-threaded one keep the same reads and `merge` stays exact.
///
/// A subsample is exact rather than approximate: thinning a locus's reads uniformly leaves
/// the bucket counts distributed exactly as they would be at the lower depth, which is the
/// same argument `max_site_depth` rests on (generic arch §2.2). What it costs is precision.
///
/// **It is not the memory knob an earlier draft of this comment called it.** Measured on
/// HG002 at 300× over the GIAB tandem-repeat set: uncapped, the table is **0.43 entries a
/// locus** — 12,727 for 29,811 loci, 0.36 MB — because most loci at a clean tract are
/// "every read at the reference length" and what separates two entries is mostly their
/// depth (spec §4.1). Deep data deduplicates; it does not give one entry per locus.
///
/// **Nor is it a correctness limit.** It looked like one — the exact-bias measurements
/// stopped at 12 reads a locus, so above that the scoring rule was unmeasured. Measured
/// since, on a narrower recorded range that keeps the cell space affordable: the rule is
/// **exactly unbiased at every depth to 45 reads a locus**, 0.00% on the rate and 0.0000 on
/// both shares (research note §6.8). Raising this constant does not step outside the
/// evidence.
///
/// **What it does decide is the counter width, and the two are one decision.**
/// `LocusShape` holds its counts in `u8`; a cap above 255 wraps them silently. Twelve is a
/// low starting value whose only remaining cost is the precision of the reads it drops —
/// a trade, not a limit. Spec §8.8.
pub const MAX_LOCUS_READS: u32 = 12;
```

### 2.2 The entry, and the table of entries

**What an entry is: one locus's shape.** How many of that locus's reads showed each whole-repeat
offset from the reference tract length, plus how many showed a length that is not a whole number of
copies. Loci whose shapes are identical are counted together, exactly as the generic path counts
sites that looked alike — the table holds *how many loci in this stratum looked like this*, and
which loci they were is never asked again.

**Why the locus and not the read, in one line, because it is the revision that moved this type:** a
read carries no genotype, so a table of reads holds the allele spectrum convolved with the slippage
kernel, and the fitted slippage rate then moves **333-fold with the starting point** (spec §4.1).

```rust
/// One locus's reads, laid out across the offset buckets. `counts.iter().sum() +
/// not_whole_repeat == depth`, always — the invariant the "no bucket is charged a negative
/// number of reads" gate checks (spec §10.1).
///
/// **The depth is exact, not binned, and that follows from the cap rather than from a
/// separate decision.** `MAX_LOCUS_READS` bounds it at a dozen values, so there is nothing
/// for a ladder to save. If the cap ever rises far, the generic path's measurement applies
/// and is not free: its depth ladder turned out to be a **correctness** parameter, sixteen
/// bins costing 0.55 rungs of the error-rate ladder where twenty cost 0.05
/// (research note §4.3).
///
/// Ordered and hashable so it can key a table and so iteration order is fixed, which is
/// the whole of the determinism requirement (shared spec §6).
pub struct LocusShape {
    /// Reads at each bucket, in bucket order. Index `OFFSET_HALF_RANGE` is the origin.
    pub counts: [u8; OFFSET_BUCKETS],
    /// Reads differing from the reference by something that is not a whole number of
    /// copies. **A guard, not a parameter**: it factorises out of the likelihood exactly,
    /// so nothing is estimated from it and it disturbs nothing (spec §4.1).
    pub not_whole_repeat: u8,
}

/// One stratum's evidence, for one read group.
///
/// **A locus covered by two read groups makes two entries, and that is sound.** Each
/// entry's own distribution is correctly specified — the genotype is drawn once for the
/// locus and enters both through the same mixture — so the product over a locus's entries
/// is a composite likelihood and the split costs precision, not correctness (spec §4.1,
/// and the generic path's measurement of the same split in research note §2.6). What must
/// not be split is a locus's reads *within* one read group, which is what the entry key
/// exists to prevent.
///
/// **Sparse, not dense, and that is forced rather than chosen.** The generic path's cell
/// space is `(depth bin × alt count)` — a few hundred cells, dense and ragged. Here an
/// entry is a whole locus's split across ten buckets (nine offsets and the guard), so the
/// possible space is `C(depth + 9, 9)` at one depth — 220 at three reads, 293,930 at
/// twelve — and only a small, data-dependent corner of it is ever occupied. A `BTreeMap`
/// rather than a `HashMap`, for the ordering reason the generic accumulators give.
///
/// **A stratum's entries are bounded twice, and which bound binds depends on the cap.**
/// By its loci, obviously; and by the shapes a capped locus can take at all, which is
/// `C(cap + 10, 10) − 1` and **does not depend on the genome**:
///
/// | `MAX_LOCUS_READS` | shapes a stratum can ever hold |
/// |---:|---:|
/// | 4 | 1,000 |
/// | 8 | 43,757 |
/// | 12 | 646,645 |
/// | 20 | 30,045,014 |
///
/// So a small cap gives a hard per-stratum ceiling — at four reads a stratum cannot exceed
/// a thousand entries however many million loci it holds — while from about twelve upward
/// the ceiling stops binding and the locus count is what limits the table. The measured
/// numbers sit far below both (spec §4.1), so neither bound is the operative constraint on
/// real data; they are here because a reader sizing the structure will ask.
pub struct StratumTable {
    loci: BTreeMap<LocusShape, u32>,
    /// The composition channel: two running counts, pooled over the stratum's reads.
    /// `ε` is fitted from these alone, in closed form (spec §4.1).
    bases_compared: u64,
    bases_mismatched: u64,
}

impl StratumTable {
    /// Add one locus. Never allocates after the first locus of a given shape.
    pub fn add_locus(&mut self, shape: LocusShape, bases_compared: u32, bases_mismatched: u32);
    /// Element-wise; associative and exact, so shards merge to the table of the union.
    pub fn merge(&mut self, other: &Self);
    /// Every distinct shape and how many loci had it, materialised once. **A `Vec` rather
    /// than an iterator**, because the search re-walks these once per candidate and a
    /// `BTreeMap` walk is not cheaply re-walkable (generic arch §2.2's argument).
    pub fn shapes(&self) -> Vec<(LocusShape, u64)>;
    pub fn loci(&self) -> u64;
    /// Mismatched over compared — the maximum-likelihood substitution rate, and a division
    /// rather than a search (spec §4.1).
    pub fn substitution_rate(&self) -> Option<ErrorRate>;
    /// The share of reads differing from the origin — **the reference tract length, not the
    /// allele** — that did so by something other than a whole number of copies. §5's
    /// diagnostic, and the number `GUARD_SHARE_LIMIT` is compared against.
    ///
    /// The allele is what the *model* means by a slip and is not available here, so a real
    /// non-reference allele enters this denominator and not its numerator: the reported
    /// share is diluted relative to the model's, never inflated (spec §4.1). A stratum that
    /// crosses the limit on this number has crossed it on the model's too.
    pub fn not_whole_repeat_share(&self) -> f64;
}

/// Above this share of the reads that differ from the origin, the stratum is one this
/// noise model does not describe, and its fitted slippage is mostly mis-modelled ordinary
/// indel however many loci stood behind it (spec §5).
///
/// **One in ten, and the bands it separates are 0.9% against 33.8% and 58.5%** — a factor
/// of three either way, so a stratum crossing it is unambiguous. *Soft*: three bands of one
/// dataset, and the number moves if the per-stratum distribution turns out continuous.
pub const GUARD_SHARE_LIMIT: f64 = 0.10;
```

### 2.3 One locus → one entry

**What it does:** reduces one `LocusKind::Ssr` locus to the shape the table is keyed on, and to the
two composition counts. Its own file because it is the only place that decides what a read's repeat
count *is*, which the locus type cannot answer.

```rust
/// The reference tract's repeat count and period — the stratum this locus belongs to.
/// Derived from the **reference**: `reference_bases.len() / motif.period()`, both of which
/// the locus carries (`SampleLocusObservations::reference_bases`,
/// `SsrDetail::motif`). A tract whose reference length is not a whole number of copies is
/// counted and skipped, not rounded.
pub fn stratum_of(locus: &SampleLocusObservations) -> Option<Stratum>;

/// Reduce one locus to its shape.
///
/// **Complete witnesses only** (`complete_observations()`,
/// locus_generation/mod.rs:134). A partial witness saw part of the tract, so its length is
/// a lower bound; scoring it as a length is the mis-scoring that iterator exists to guard,
/// and on this path it would read as a read that lost repeats — a direct bias in the
/// direction split, which is the parameter §3 of the spec exists to protect.
///
/// **`reads_without_observation` does not enter the depth** and is counted instead
/// (locus_generation/mod.rs:53); those reads covered the tract and witnessed nothing.
///
/// **`reads_discarded_by_cap` does not skip the locus** (locus_generation/mod.rs:56): the
/// generator's own reservoir cap is a random subsample (`DEFAULT_SSR_MAX_READS_PER_LOCUS`,
/// locus_generation/ssr.rs:125), so a locus it fired on is a locus observed at the lower
/// depth. It is counted, because a run where it fires everywhere is a run whose depths are
/// the cap's and not the data's.
///
/// A read's offset is `(observation.bases.len() − reference_bases.len()) / period`, whole
/// only when that difference divides by the period; otherwise the read goes to
/// `not_whole_repeat`.
pub fn shape_of(
    locus: &SampleLocusObservations,
    read_group: ReadGroupId,
    cap: u32,
) -> LocusShape;

/// The composition channel for one locus: bases compared and bases mismatched, over the
/// complete witnesses of one read group.
///
/// `OPEN:` **where the mismatch count comes from.** Nothing on the locus carries one —
/// `SequenceObservation` holds `q_sum`, a quality sum, not a count
/// (locus_generation/mod.rs:186). Two ways to close it, and the choice is impl-time:
/// compare `observation.bases` against the motif tiled to that length here, which is exact
/// for a perfect tract and charges an impure tract's interruptions to `ε`; or have the STR
/// aligner emit the count it already computes while scoring
/// (`alignment/emission.rs:43`, the match/mismatch pick). **The second is preferable and
/// larger**: it touches step 3's output. The first is what this unit can do alone, and its
/// caveat — that `ε` then absorbs tract impurity, which `SsrSegment::purity_fraction()`
/// (region_typing/segment_criteria.rs:209) makes measurable per stratum — belongs in the
/// emitted provenance either way.
pub fn composition_of(
    locus: &SampleLocusObservations,
    read_group: ReadGroupId,
) -> (u32, u32);
```

### 2.4 The fitted output

```rust
/// The three numbers that describe how a read's repeat count moves away from its allele.
/// Named for what each **is**, not for its place in a formula.
pub struct SlippageModel {
    /// How often a read shows a length other than its allele's.
    pub slip_rate: SlipRate,
    /// Of the reads that slipped, the share that **gained** repeats. 0.17 at tomato
    /// dinucleotides — a read is 4.9 times as likely to lose a repeat (spec §3).
    pub gain_share: SlipGainShare,
    /// Of the reads that slipped, the chance of a second step given a first. **One number
    /// for both directions** (spec §3), and 0.065 to 0.12 on real data.
    pub step_decay: SlipStepDecay,
}

/// Three constrained rates, one type each — not one shared `Probability`, for the reason
/// generic arch §2.1 gives: they are all fractions in `[0, 1]`, so one type would let a
/// direction split be handed to something expecting a slippage rate and compile. Each
/// follows `MismatchFraction`'s shape (types.rs:243): private field, `try_new`, `.get()`.
pub struct SlipRate(f64);
pub struct SlipGainShare(f64);
pub struct SlipStepDecay(f64);

/// What one stratum's fit returned.
pub struct StratumFit {
    pub stratum: Stratum,
    pub slippage: Estimate<SlippageModel>,
    pub substitution: Estimate<ErrorRate>,
    /// The allele-length distribution the fit weighed the genotype against — one entry per
    /// unordered pair of allele offsets (§3). Emitted because it is what a reader needs
    /// when a slippage rate looks wrong, and because the gather has a use for it (spec §6).
    pub genotypes: Vec<(WholeRepeatOffset, WholeRepeatOffset, f64)>,
    /// §5's diagnostic. Above `GUARD_SHARE_LIMIT` the fitted slippage describes ordinary
    /// indel rather than repeat slippage.
    pub not_whole_repeat_share: f64,
    /// How the search went, and it is what separates a fitted number from a stopped
    /// search: the slippage each starting point reached and what it scored, best first.
    pub starts_tried: SmallVec<[SlippageStart; 4]>,
    /// Which strata this fit's loci actually came from — itself where it was fitted in
    /// place, its neighbours where it borrowed, both where a merge fired (§4.2).
    pub fitted_over: SmallVec<[Stratum; 2]>,
}

pub struct SlippageStart { pub start: SlippageModel, pub reached: SlippageModel, pub log_likelihood: f64 }

/// Everything the STR path estimates for one sample.
pub struct SsrSampleParameters {
    pub by_stratum: BTreeMap<(ReadGroupId, Stratum), StratumFit>,
    /// **The thing a person reads.** Several hundred fits per sample make a per-stratum
    /// record a file nobody opens (spec §4.4).
    pub summary: BTreeMap<ReadGroupId, StratumFitSummary>,
}
```

## 3. The fitting machinery, and what the STR path needs that the generic path does not

The generic architecture's §4 specifies `fitting/` as path-agnostic and this unit is its second
consumer. **Two of the three pieces transfer unchanged; the third does not, and that is a finding
about the shared design rather than a licence to restate it here.**

**Transfers unchanged — `NoiseModel`.** Its associated `Cell` and `NoiseParams` are exactly the seam
this path needs:

```rust
impl NoiseModel for SsrNoiseModel {
    type Cell = LocusShape;
    type NoiseParams = SlippageModel;

    /// How likely each genotype makes this locus's shape, at these slippage parameters.
    ///
    /// A genotype is an unordered pair of allele offsets. Each read picks one of the two
    /// copies and then slips, so the per-bucket probability is the average of the two
    /// copies' slip kernels, and a shape's probability is one multinomial over the buckets.
    ///
    /// **An end bucket's probability is the sum over every offset it absorbs**, never the
    /// probability of sitting exactly on the edge. Measured: the marginal is exactly
    /// unbiased at every range tried; plugging in the edge fails the sums-to-one gate
    /// (0.9488 at ±1) and, rescaled, costs +33% of the slippage rate where 30 in 100
    /// slipped reads take a second step (spec §4.1).
    fn genotype_likelihoods(
        &self, cell: &LocusShape, noise: &SlippageModel, ploidy: Ploidy, out: &mut Vec<f64>,
    );
}
```

**Transfers with one signature change — `fit_mixture_weights`.** The climb over the genotype
frequencies is the same concave problem and the same code. Its declared return type is not:

```rust
// generic arch §4.1, today
pub fn fit_mixture_weights(..) -> SmallVec<[f64; 3]>;
```

Three is the diploid generic path's genotype count. **At `ALLELE_OFFSET_LIMIT = 6` the support runs
±6 around the reference length, so a stratum has up to 13 allele lengths and up to `13·14/2 = 91`
genotypes** — *up to*, because an allele cannot be shorter than nothing: a stratum at 4 repeats
reaches only −4, giving 11 lengths and 66 genotypes, and only strata at 6 repeats and above get the
full 13. So the count is a per-stratum quantity bounded by 91, and **the return type has to widen to
a `Vec<f64>` or be generic in its inline capacity**. One line, and it belongs in the shared module
rather than being copied here. *(The spec works its examples at nine lengths and 45 genotypes, which
is the same arithmetic at a narrower support.)*

**Does not transfer — `fit_by_profile_scan`.** It steps a ladder end to end because nobody had shown
the profile curve has a single hump. Two things stop it here, both in spec §4.2: a flat scan over
three slippage parameters is 4.2 million scores **per (read group × stratum)** with several hundred
strata per sample; and a quarter-Phred ladder is the wrong spacing for two parameters that are
shares in `(0, 1)` rather than rates spanning orders of magnitude. So `fitting/` gains a second
driver beside the scan:

```rust
/// Maximise the noise parameters from several starting points, climbing the genotype
/// frequencies at every trial (`fit_mixture_weights`), and return the best-scoring with
/// **every start's outcome beside it**.
///
/// The starts must disagree about **every** axis the fit can stick on, not only the
/// headline one — the trap the runs model fell into, where five starts sharing one guess
/// at a nuisance axis returned a confident zero (generic spec §6.5). `SLIPPAGE_STARTS`
/// below is that set for this path.
///
/// **Convergence failure is a data condition here, unlike the inner climb.** The climb
/// over the frequencies is concave and a failure there is a bug; the outer search over the
/// slippage parameters has no such proof, so it is capped, the best-scoring iterate is
/// kept, and the termination is reported — the same treatment as the generic path's
/// coupled fit (generic arch §5.2).
pub fn fit_by_multistart<M: NoiseModel>(
    model: &M,
    cells: &[(M::Cell, Ploidy, u64)],
    starts: &[M::NoiseParams],
) -> MultistartResult<M::NoiseParams>;

/// What a multi-start search returned: the winner, and every start's outcome beside it.
pub struct MultistartResult<P> {
    pub best: P,
    pub frequencies: Vec<f64>,
    pub log_likelihood: f64,
    /// Sorted by score, best first. Never empty.
    pub starts: SmallVec<[StartOutcome<P>; 4]>,
    /// The highest-to-lowest ratio across `starts`, in the headline parameter, against
    /// `START_AGREEMENT_LIMIT`.
    pub starts_disagreed: bool,
}
pub struct StartOutcome<P> { pub from: P, pub reached: P, pub log_likelihood: f64 }

/// Four starts, each disagreeing about the rate, the direction and the decay at once —
/// the four the harness reports every number from (spec §4.1, §4.2).
///
/// **The rate is a multiplier on a moment estimate, not an absolute value**, because a
/// stratum's rate spans twenty-two-fold across repeat counts (spec §4) and a fixed ladder
/// of absolute rates would start every stratum in the wrong place. The moment estimate is
/// the share of reads sitting off the origin — an over-estimate, since it counts real
/// alleles too, which is why the spread runs below it as well as above.
pub const SLIPPAGE_STARTS: [(f64, f64, f64); 4] = [
    // (multiplier on the moment estimate, gain share, step decay)
    (3.0, 0.20, 0.03),
    (1.0 / 3.0, 0.80, 0.40),
    (1.0, 0.50, 0.15),
    (0.3, 0.35, 0.08),
];

/// How far two starting points may land apart before the fit is reported as not having
/// found an answer. **One quarter-Phred — 6% — borrowed from the generic path's ladder
/// spacing**, which is the argued size of the finest difference a caller can feel
/// (shared spec §3). Two starts closer than that agreed; two starts further apart did not.
pub const START_AGREEMENT_LIMIT: f64 = 1.06;
```

**Contract.** `fit_by_multistart` returns the highest-scoring outcome **and the full set**, and it
flags a fit whose starts spanned more than `START_AGREEMENT_LIMIT` in the rate — the flag §4.3's
summary counts and §4.2's `SlippageNotIdentified` is raised from. It never returns a single answer
without the set behind it: an answer with no spread beside it is indistinguishable from a search
that never looked, which is the failure the generic path's runs model produced as a confident zero.

## 4. The accumulator and the four fits

```rust
/// Step 4's STR front door: one table per (read group, stratum).
pub struct SsrAccumulators {
    by_stratum: BTreeMap<(ReadGroupId, Stratum), StratumTable>,
    ploidy: Arc<dyn PloidyMap>,
    counts: SsrAccumulationCounts,
}

impl SsrAccumulators {
    /// One per **region shard**.
    pub fn new(read_groups: &[ReadGroupId], ploidy: Arc<dyn PloidyMap>) -> Self;

    /// Add one locus. **Borrows** — the caller keeps it and passes it on unchanged. A
    /// locus whose `kind` is not `LocusKind::Ssr` is ignored; one whose reference tract is
    /// not a whole number of motif copies is counted and skipped (§2.3).
    pub fn add_locus(&mut self, locus: &SampleLocusObservations);

    /// Combine a shard's tables into these. Associative and exact.
    pub fn merge(&mut self, other: Self);

    /// Everything the accumulator did to a locus other than enter it as it arrived.
    pub fn adjustments(&self) -> &SsrAccumulationCounts;
}

/// Every field a plain sum, so it merges; every field reported, because the alternative is
/// a fit that quietly describes a different population of reads than the caller will see.
#[derive(Default)]
pub struct SsrAccumulationCounts {
    /// Loci this unit subsampled down to `MAX_LOCUS_READS`. A run where most loci appear
    /// here is one where the cap, not the data, is setting the depth.
    pub loci_subsampled_to_cap: u64,
    /// Loci the generator's own reservoir had already subsampled
    /// (`reads_discarded_by_cap`). Entered anyway, at the depth observed (§2.3).
    pub loci_with_upstream_subsample: u64,
    /// Reads that covered a tract and witnessed nothing. Not part of any depth.
    pub reads_without_observation: u64,
    /// Reads whose witness was partial, so their length is a lower bound. Excluded (§2.3),
    /// and a large share here means the read length is short against these tracts.
    pub reads_with_partial_witness: u64,
    /// Loci whose reference tract is not a whole number of motif copies, so no stratum
    /// holds them. **Should be near zero**; region typing delimits on whole copies.
    pub loci_without_whole_repeat_reference: u64,
}
```

**Contract.** `add_locus` reads the locus only through `region`, `reference_bases`, `kind` and
`complete_observations()`. **Unlike the generic path, overlapping loci are not a hazard here**: an
STR locus is one delimited tract and two tracts do not overlap, so there is no partition invariant
to count against. What replaces it is `loci_without_whole_repeat_reference`, which must read near
zero and is a bug report against region typing if it does not.

### 4.1 The four fits, in order

1. **The substitution rate.** `bases_mismatched / bases_compared` per stratum. A division (spec
   §4.1), not a search, and it needs none of the other three.
2. **The three slippage parameters.** `fit_by_multistart` over that stratum's shapes, with
   `fit_mixture_weights` climbing the genotype frequencies at each trial.
3. **Borrowing**, for a stratum below `MIN_LOCI_TO_FIT`: take the neighbouring repeat counts at the
   same period, marked `Provenance::Borrowed` with `fitted_over` naming them.
4. **The monotonicity walk**, last, because it reads the fitted sequence: visit each period's strata
   in repeat-count order and where a fitted `slip_rate` falls below its predecessor's, merge the two
   strata's tables and refit, repeating until the sequence rises.

```rust
/// Fewest loci a stratum needs before it is fitted rather than borrowing. *Soft.*
pub const MIN_LOCI_TO_FIT: u64 = 1_000;

/// What a merge costs, so the constant above can be read against it: two strata pooled
/// return close to the loci-weighted mean of their rates, so each carries its own distance
/// from it — about a quarter of the rate at a 1.5-fold difference and half at two-fold
/// (spec §4.3). On real strata slippage rises about 1.3-fold per repeat count, so
/// borrowing one neighbour costs on the order of 15 to 25%.
pub fn merge_and_refit_monotone(fits: &mut BTreeMap<Stratum, StratumFit>, model: &SsrNoiseModel);
```

### 4.2 When there is not enough data

```rust
/// **This path's own error enum, not a set of variants added to the generic path's.** The
/// two units fail differently — a stratum too thin to fit has neighbours to borrow from
/// and a sample's heterozygosity does not — and step 4's own surface (`mod.rs`) wraps both
/// for a caller that drove the whole step.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum SsrEstimationError {
    /// Every stratum at this period was below `MIN_LOCI_TO_FIT`, so there is no
    /// neighbour to borrow from. **Deliberately has no default**: a slippage rate spans
    /// twenty-two-fold across repeat counts within one dataset (spec §4), so any constant
    /// would be wrong for most strata.
    #[error("sample {sample}: period {period} has no stratum with {MIN_LOCI_TO_FIT} loci \
             to fit or borrow from; supply the parameters or drop the period")]
    NoFittableStratumAtPeriod { sample: String, period: u8 },

    /// The search reached materially different answers from different starting points, so
    /// what it returned is where it stopped rather than what the data says. **Not "too
    /// little data"** — that is `Borrowed`. Spec §4.2.
    #[error("sample {sample}: stratum (period {period}, {repeats} repeats) reached slippage \
             rates spanning {spread:.1}x across {starts} starting points")]
    SlippageNotIdentified { sample: String, period: u8, repeats: u32, spread: f64, starts: usize },

    #[error(transparent)]
    Domain(#[from] DomainError),
}
```

### 4.3 The summary, which is the part a person reads

**Several hundred fits per sample against the generic path's four**, so the diagnostics aggregate
rather than accumulate (spec §4.4). The per-stratum `StratumFit` is still written — a fit that looks
wrong has to be traceable — but nothing downstream is expected to read it.

```rust
/// One read group's fits, summarised. Every field answers a question a person would
/// otherwise have to grep several hundred records for.
pub struct StratumFitSummary {
    pub strata_fitted_here: u32,
    pub strata_borrowed: u32,
    /// The merged sets, named — a merge is a claim about two strata at once.
    pub strata_merged: Vec<SmallVec<[Stratum; 2]>>,
    /// Fits whose starting points disagreed by more than `START_AGREEMENT_LIMIT` in the
    /// slippage rate, and the worst of them. This is the diagnostic the four starts exist
    /// to produce. **Not "one step of the rate's resolution"**, which an earlier draft said
    /// and which names nothing on this path: it searches rather than scanning, so it has no
    /// rungs — the limit is the generic path's ladder spacing, borrowed (§3).
    pub strata_with_disagreeing_starts: u32,
    pub worst_start_disagreement: Option<(Stratum, f64)>,
    /// Strata above `GUARD_SHARE_LIMIT` — the ones this noise model does not describe.
    pub strata_above_guard_limit: u32,
    pub worst_guard_share: Option<(Stratum, f64)>,
    /// How many loci stood behind the thinnest and the thickest fit.
    pub loci_behind_fits: (u64, u64),
}
```

## 5. Design decisions — decided

- **An entry is a locus, not a read — decided.** A read carries no genotype, so a per-read tally
  holds the allele spectrum convolved with the slippage kernel and the fitted rate moves 333-fold
  with the starting point; keyed by locus the same fit is exactly unbiased (spec §4.1). This is the
  STR twin of the generic path's rejection of a windowed histogram keyed per read group.
- **Offsets are measured from the reference tract length — decided**, closing what the spec carried
  as `OPEN` before its 2026-08-06 revision, and against the leaning it recorded. The modal origin returns
  a slippage rate 50% to 408% high and a direction split of 0.48 where the truth is 0.17 — a
  1.1-fold asymmetry where the truth is 4.9-fold (spec §4.1). It also makes this table and the STR
  census share an origin, which they did not before.
- **The stratum key is the reference tract's repeat count — decided**, and it follows from the
  origin: both are pure functions of the reference, so every sample strata-fies identically and a
  cohort can compare strata (§2.1).
- **An end bucket is scored by summing over what it absorbs, never by plugging in its edge —
  decided.** Rejected: plugging in the edge, which sums to 0.9488 rather than one at ±1 and so is
  not a likelihood; and the rescaled plug-in, which is proper and costs +33% of the rate where
  slippage is heaviest (spec §4.1).
- **The substitution rate is a division, not an axis of the search — decided.** The length channel
  and the composition channel factorise exactly, so the two pooled counters are a sufficient
  statistic (spec §4.1). This is also what corrects the shared spec's "the STR path's three, scanned
  together" arithmetic.
- **The search is several starts rather than a flat scan — decided**, against
  [`../spec/parameter_prepass.md`](../spec/parameter_prepass.md) §3.1's rule for this path only. The
  scan's justification was that nobody had shown the profile curve has one hump; measured, it does,
  on two worlds and one axis (spec §4.2). The scan would also be 4.2 million scores per stratum.
- **Genotype frequencies are fitted freely over unordered allele pairs — decided**, matching the
  generic path's §11.4 rather than tying them through one allele frequency, which would presume the
  inbreeding coefficient is zero. `OPEN:` whether a tied form using the sample's fitted `F` should
  be used where a stratum is too thin for the free one (spec §8.6).
- **A locus's reads are capped by subsampling, not by dropping the locus — decided.** Dropping deep
  loci would be depth-dependent selection, which is the bias step 4 exists to remove; a uniform
  subsample is exact and costs precision only (§2.1).
- **No trait over the accumulator — decided.** Nothing generic drives it, and the walk knows which
  object it is filling. `fitting/` is the only place with a genuine swappable seam.

## 6. Reconciliation with existing code

Every row read before it was written.

| this doc | existing code | action |
|---|---|---|
| the input stream | `SampleLocusObservationsIterator` ([locus_generation/mod.rs:701](../../../../src/ng/locus_generation/mod.rs)) | consume by reference; unchanged |
| the locus | `SampleLocusObservations` ([locus_generation/mod.rs:40](../../../../src/ng/locus_generation/mod.rs)) | reuse as-is |
| routing to this path | `LocusKind::Ssr(SsrDetail)` ([locus_generation/mod.rs:217, 229](../../../../src/ng/locus_generation/mod.rs)) | match on it; `SsrDetail::motif` is the period source |
| the reference tract | `SampleLocusObservations::reference_bases` ([locus_generation/mod.rs:46](../../../../src/ng/locus_generation/mod.rs)) | the origin and the stratum's repeat count both derive from its length |
| a read's observed tract | `SequenceObservation::bases` ([locus_generation/mod.rs:159](../../../../src/ng/locus_generation/mod.rs)) — "allele content, in **read** coordinates" | the offset is `(bases.len() − reference_bases.len()) / period` |
| the scorable subset | `complete_observations()` ([locus_generation/mod.rs:134](../../../../src/ng/locus_generation/mod.rs)) | call directly — §2.3's guard, and on this path a partial reads as a lost repeat |
| per-read-group support | `SequenceObservation::read_group` ([locus_generation/mod.rs:178](../../../../src/ng/locus_generation/mod.rs)) — already part of the identity, and its doc names this path as the near-term consumer | reuse as-is |
| read support | `SequenceObservation::num_obs` ([locus_generation/mod.rs:181](../../../../src/ng/locus_generation/mod.rs)) — "the whole support on the STR path" | the bucket increment |
| the upstream cap | `reads_discarded_by_cap` ([locus_generation/mod.rs:56](../../../../src/ng/locus_generation/mod.rs)), set by `DEFAULT_SSR_MAX_READS_PER_LOCUS` ([locus_generation/ssr.rs:125](../../../../src/ng/locus_generation/ssr.rs)) | count it; do not skip the locus |
| the motif | `Motif` ([types.rs:338](../../../../src/ng/types.rs)), `period()` ([types.rs:369](../../../../src/ng/types.rs)) returning `usize` | reuse; add `SsrPeriod` beside it rather than changing the accessor's type under its callers |
| the tract's coordinates and purity | `SsrSegment` ([region_typing/segment_criteria.rs:141](../../../../src/ng/region_typing/segment_criteria.rs)), `period()` `:215`, `tract_len()` `:223`, `purity_fraction()` `:209` | **not consumed by this unit** — the locus carries what is needed. `purity_fraction` is what would make §2.3's `ε`-absorbs-impurity caveat measurable |
| `ReadGroupId` | [types.rs:199](../../../../src/ng/types.rs) | reuse as-is |
| the read-admission policy | `ReadFilterConfig` ([read/filtering.rs:63](../../../../src/ng/read/filtering.rs)) | carried into `SsrEstimationConfig` and emitted with the parameters, for the reason generic spec §2 gives — an `ε` describes the reads that survived admission |
| `DomainError` | [types.rs:268](../../../../src/ng/types.rs) — `#[non_exhaustive]`, doc already says later constrained types add their own variants | extend with the three slippage newtypes' variants |
| the checked-newtype shape | `MismatchFraction` ([types.rs:243](../../../../src/ng/types.rs)) | copy: private field, `try_new`, `.get()` |
| `ErrorRate`, `Estimate<T>`, `Provenance`, `Ploidy`, `PloidyMap` | [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1, §2.4, §3 — new there, not yet in code | consume; this unit adds no second copy |
| `NoiseModel` | [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §4.2 | implement — this is its second implementation and the reason it is a trait |
| `fit_mixture_weights` | [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §4.1 | reuse, **after widening its return type past three genotypes** (§3) |
| `fit_by_profile_scan` | [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §4.2 | **not used**; `fitting/` gains `fit_by_multistart` beside it (§3) |
| the mismatch score | `alignment/emission.rs:43` (the match/mismatch pick) | **not reused today** — it scores, it does not count. `OPEN:` §2.3 |
| the per-sample summary | [`ng_step_interfaces.md`:343](ng_step_interfaces.md) | re-specified with the generic path's: the input is accumulated statistics, not `&[ConfidentGenotype]` (shared spec §2.3) |
| production's stutter pre-pass | [`src/ssr/cohort/prepass.rs`](../../../../src/ssr/cohort/prepass.rs) | **not reused** — it pools reads from loci that passed a confident-genotype gate, which is both biases this design removes. Frozen production |

## 7. Open items

- **Impl-time confirmation, not an open design question: `MAX_LOCUS_READS`.** Neither of the two
  things that would have made it one survives — the table's memory is measured and small, and the
  scoring rule is exactly unbiased to 45 reads a locus (§2.1). What is left is the precision of the
  reads it drops, against the `u8` ceiling it shares with `LocusShape`.
- `OPEN:` **where the mismatch count comes from** — derived here under a tiling assumption, or
  emitted by the aligner that already computes it (§2.3). The second is better and touches step 3.
- `OPEN:` **free genotype frequencies or an allele spectrum plus the sample's `F`** for a thin
  stratum (§5, spec §8.6). **Settled by:** the exact-bias harness, which fits both against one truth.
- **Impl-time confirmation.** `MIN_LOCI_TO_FIT`, and the resolution against which "starting points
  disagreed" is judged.
- **A change this unit forces on the shared module**, and it is not optional: `fit_mixture_weights`
  returns `SmallVec<[f64; 3]>`, which cannot hold a stratum's genotype count (§3).

## 8. Test & bench shape

Unit tests beside each file; the acceptance tests are the spec's §10 and
[`../spec/parameter_prepass.md`](../spec/parameter_prepass.md) §10.

**The harness is the evidence that does not come from the model**, and it already exists:
[`ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs), written up in
[`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§6. **Any change to `genotype_likelihoods` or to the entry key re-runs its three algebraic gates
first** — the rule sums to one over the entry space, no bucket is charged a negative number of
reads, and a silent kernel puts every locus's reads on its own alleles. Each is one line and each
rejects a broken rule without fitting anything.

Five anchors:

- **The control, and it must read exactly zero.** Key by locus with the reference origin, generate
  and fit under the same key: 0.000% on the rate, 0.0000 on both shares, four starts agreeing to
  1.000×. Anything else is the harness's fault, not the estimator's.
- **Agreement with truth where truth exists.** On HG002, the fitted parameters against those
  measured directly on known-homozygous loci — 2.0% at six or more repeats and a 3.4× direction
  split (spec §10.3). This is the only check in the design that does not generate its data from the
  model it then fits, and it is the test production's estimator fails.
- **Sharded accumulation is exact** — one sample walked in one region and in many gives identical
  tables. Integer entry counts make it an equality rather than a tolerance, which is why the read
  cap is seeded from the locus's position (§2.1).
- **`adjustments().loci_without_whole_repeat_reference` is near zero** on both real cohorts. A large
  count is a bug report against region typing.
- **The summary is what fails, not a per-stratum record.** Feed a deliberately unfittable stratum —
  loci generated with the slippage rate at zero and the alleles spread — and assert it appears in
  `StratumFitSummary`. A record only a debugger would open does not satisfy this (spec §10.6).

No `bench/`: this unit has no competing implementations, and its cost claims are measured rather
than open — [`ng_str_table_memory.rs`](../../../../examples/ng_str_table_memory.rs) drives the real
region-typing walk and STR locus generator over real alignments and prices the table on them
(research note §6.8). **Re-run it when the entry key changes**, which is the one change that could
move the numbers.
