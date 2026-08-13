# ng — the joint fit, which loci are kept: types & interfaces

*Status: **partly built** (2026-08-12) — `loci.rs` exists and holds the generic rule, the
unambiguous-base mask and the threshold arithmetic, with the properties §6 asks for asserted in its
own tests. `KeptLoci`, `SelectionIdentity` and `KeptLociDigest` are still the contract below and not
yet code. Companion to the spec
[`../spec/parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md) (the design and
its rationale) and to the shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md)
(vocabulary) and [`module_layout.md`](module_layout.md) (the `src/ng/` tree). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): nouns for types, verbs
for functions, **STR** in prose ↔ `ssr` in code. Signatures are illustrative; the **contract** is the
deliverable. Every "why" lives in the spec — this doc does not re-argue one.*

*One of three covering the per-site route: this one, [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md)
(what is stored at a kept locus) and [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md)
(the estimator).*

## Module home

`src/ng/parameter_estimation/joint/` — a sibling of `generic/` and `ssr/` under step 4, because it is
a second route to the same parameters rather than a different step (`module_layout.md`,
*parameter_estimation/*). This document owns one file of it.

```
src/ng/parameter_estimation/joint/
├── mod.rs      – the route's surface; re-exports the three units below
├── loci.rs     – THIS DOCUMENT: which loci are kept, and the identity that travels with them
├── records.rs  – what is recorded at each (its own arch doc)
└── fit.rs      – the estimator (its own arch doc)
```

**One file, not a folder, and no trait.** The selection has exactly one implementation and no
competing candidate: the generic rule is a hash comparison and the STR rule is two calls on the
repeat catalog. A trait here would be ceremony over a pure function
([`module_structure.md`](../../../../ai/skills/rust-code-review/code_review/module_structure.md);
spec §1).

**Nothing is added to `repeat_catalog/`.** The STR sampler this consumer asked for already ships —
`RepeatCatalog::sample_loci_per_stratum` and `count_loci_per_stratum`
([`src/ng/repeat_catalog/reader.rs:242,260`](../../../../src/ng/repeat_catalog/reader.rs)) — so this
unit states a policy and calls them (spec §3.2).

---

## 1. Types

### 1.1 `KeptLoci` — the selection, as one value

The two sets, and the counts a consumer is obliged to reweight by (spec §3.5). Held together because
they are produced by one pass over one set of inputs and checked for agreement as one thing.

```rust
/// The loci every sample keeps raw evidence at, for one run.
///
/// **Reproducible, not transported**: two samples that build this from the same
/// [`SelectionIdentity`] get byte-identical contents, which is what lets samples be
/// walked on different machines (spec §1). Nothing here is written to a sample's
/// records — [`SelectionIdentity`] is (spec §5.1).
pub struct KeptLoci {
    /// The generic set, ascending in genome order. The fit indexes records by
    /// position in this vector, so **the order is part of the contract**, not an
    /// artifact of how it was built.
    generic: Vec<GenomePosition>,
    /// The STR set, keyed by stratum, as the catalog's sampler returns it.
    ssr: StratumSample,
    /// How many loci each stratum **holds**, against how many were kept. Anything
    /// pooled across strata is biased without this and silently so (spec §3.5).
    ssr_stratum_counts: StratumCounts,
}
```

**Contract.** `generic` is sorted, duplicate-free, and every entry lies inside the analysed regions.
Building it twice from one `SelectionIdentity` yields equal values, whatever the region sharding or
thread count — the property spec §7.1 tests. `ssr_stratum_counts` is over the same `scope` as the
sample, never over the whole reference (spec §3.3).

### 1.1a `SelectableRegions` and `UnambiguousRuns` — the domain, and its length, as one value

**Built.** The two failures this pair exists to make unrepresentable are the same failure at two
levels: computing the threshold from one notion of "how much genome is in play" and sweeping over
another. Under a `--regions` BED the count then comes out of the genome's length; with an
`N`-masked domain it comes out of the contig's.

```rust
/// The stretches a run may keep positions from, **and their total length**. One value
/// carries both, so `threshold_for`'s denominator cannot be stated separately from what
/// the sweep reads (spec §2, §7.1, §7.3).
pub struct SelectableRegions { /* sorted, disjoint, non-inverted */ }

impl SelectableRegions {
    pub fn new(spans: Vec<GenomeRegion>) -> Result<Self, SelectionError>;
    pub fn spans(&self) -> &[GenomeRegion];
    pub fn total_length(&self) -> Bp;
}

/// The runs of `A`/`C`/`G`/`T`, collected as the reference streams past — a
/// `ReferenceBasesObserver`, so the mask rides on the pass that computes the digests and
/// costs no second read (spec §2).
pub struct UnambiguousRuns { /* … */ }

impl UnambiguousRuns {
    pub fn ambiguous_bases(&self) -> u64;
    pub fn into_selectable(self) -> Result<SelectableRegions, SelectionError>;
}
```

**Contract.** `new` sorts, then refuses an inverted or an overlapping span — a position offered twice
is a silent double weight in every rate fitted from the set, and no later check would see it.
`total_length` is the sum of the spans it holds and is never recomputed from a contig table.

### 1.2 `SelectionIdentity` — the seven values that travel

What a sample carries so the fit can refuse to pool mismatched runs (spec §5.1). **This is the type
that crosses the machine boundary**, not `KeptLoci`.

```rust
/// Everything that decides which loci a run keeps. Two samples that disagree on any
/// field selected different loci, so the fit must refuse rather than average
/// (spec §5.1).
///
/// **`PartialEq` and not `Eq`, and that is forced rather than chosen**: `StrRepeatCriteria`
/// carries a purity floor and a score floor, which are floats. **Compare them by bit
/// pattern, not by `==`** — this is an identity check, so "the same value" means the same
/// bits, and the `f64` `==` that `derive(PartialEq)` generates answers `false` for two
/// `NaN` floors that came from the same configuration file.
#[derive(Clone, PartialEq)]
pub struct SelectionIdentity {
    /// The run's selection seed. A different seed selects a disjoint set and both
    /// look well-formed.
    pub seed: u64,
    /// Content digest of the reference, as `reference_info` computes it.
    pub reference: ReferenceDigest,
    /// Content digest of the analysed region set. The likeliest to differ by
    /// accident, because a BED feels like a runtime convenience and is not.
    pub analysed_regions: RegionSetDigest,
    /// What the catalog file was built at — floors, period range, flank, and the two
    /// Ruzzo–Tompa weights. A different weighting is a different set of tracts, not
    /// a subset of one.
    pub catalog_built_under: CatalogBuildSettings,
    /// What this run asked the catalog for: copy floors, purity floor, satellite cap,
    /// bundle radius.
    pub ssr_criteria: StrRepeatCriteria,
    /// The generic target position count. Two samples must agree even though a smaller
    /// target's positions *are* a subset of a larger one's: the records array holds no
    /// coordinates, so the same offset means a different locus under a different target
    /// (spec §4.3.4).
    pub generic_target: u64,
    /// The STR per-stratum cap.
    pub ssr_cap: usize,
}
```

**Contract.** `PartialEq` **is** the compatibility check — there is no tolerance and no partial
match, because a set difference is meaningless rather than noisy (spec §5). The fit compares every
sample's against the first and fails on the first disagreement, naming the field.

**OPEN (impl-time):** `ReferenceDigest` and `RegionSetDigest` are named here as the types
`reference_info` already produces; pin them when coding rather than minting new ones.

### 1.3 `KeptLociDigest` — the eighth value, which checks the answer

The seven above say what a run was *asked*. This says what it *produced* (spec §5.1, census sites
§5.2). Its input is fed by the record writer, not by re-running the rule — a digest derived from the
rule a second time proves only that the rule is deterministic.

**Built, 2026-08-12**, in two types rather than one: the filling half is mutable and the finished
half is not, which is what lets the finished one be compared and stops a consumer feeding a digest it
was handed.

```rust
/// A witness that two samples really did keep the same loci.
pub struct KeptLociDigest {
    whole: [u8; 16],
    /// One entry per **occupied** megabase, in genome order — 800 on tomato at a
    /// two-million budget.
    per_block: Vec<BlockDigest>,
}

/// One megabase, named — so a mismatch says where and not only that.
pub struct BlockDigest { pub contig: ContigId, pub megabase: u32, pub digest: u64 }

/// Fills one as the record writer walks the kept loci.
pub struct KeptLociDigester { /* … */ }

impl KeptLociDigester {
    /// Absorb the `index`-th kept locus. Called by the record writer, in index order.
    pub fn observe(&mut self, index: usize, locus: GenomePosition);
    /// How many have been absorbed — checked against `KeptLoci::generic().len()`.
    pub fn observed(&self) -> usize;
    pub fn finish(self) -> KeptLociDigest;
}

impl KeptLociDigest {
    /// The first megabase two digests disagree on — the whole reason for blocking.
    pub fn first_disagreement(&self, other: &Self) -> Option<(ContigId, u32)>;
}
```

**Contract.** `observe` must be called exactly once per kept locus, in ascending index order; it
**panics** when handed any other index, so a writer that skips or reorders fails loudly rather than
producing a plausible digest.

**Sixteen bytes and MD5, not thirty-two.** This is the digest the same crate already computes for a
reference's content, and the failure it guards is two cooperating runs of one pipeline disagreeing
about a list — not an adversary constructing a collision. A second hash family here would be a second
thing to keep in step for no property anyone can name.

**The per-block entry carries its coordinates.** A bare `u64` per block leaves a reader counting
entries to find out which megabase entry 412 was, and the blocks are only the *occupied* megabases,
so that count is not a coordinate. Three fields is 6.4 kB more on tomato and it is the difference
between "the digests differ" and "they differ on chromosome 3, megabase 41".

---

## 2. Interfaces

### 2.1 Building the selection

```rust
/// Derive the kept loci from the run's inputs alone — no read is opened.
///
/// # Errors
///
/// [`SelectionError::CatalogTooRestrictive`] when the run asks for a copy floor or
/// flank the catalog was not built at; the catalog refuses rather than serving a short
/// list, and the refusal is passed through with both values (spec §7.6).
pub fn select_kept_loci(
    identity: &SelectionIdentity,
    catalog: &RepeatCatalog,
    analysed: &SelectableRegions,
    unambiguous: &SelectableRegions,
) -> Result<KeptLoci, SelectionError>;
```

**Contract.** Pure over its inputs; no interior mutability, no ordering dependence, no I/O beyond the
catalog handle it is given. Runs in one forward pass over the catalog and holds `cap` values per
stratum plus the generic positions.

**Two region sets, and the intersection happens inside.** `analysed` is the run's own territory and
`unambiguous` is `UnambiguousRuns::into_selectable`'s output for the same reference. Handing the
caller the job of intersecting them is the failure this module is shaped around: the threshold's
denominator has to be the same object the sweep walks, and `SelectableRegions::intersect` is what
keeps the narrowed length travelling with the narrowed spans.

**The generic domain excludes repeat tracts, and so does its denominator.** Spec §2's rule says
tracts are excluded and separately that the threshold comes from the masked length; those are the
same length or the realised count falls short by the tract fraction. So the domain is the analysed
regions, masked, cut by the catalog, and reduced to the `Generic` pieces — and `threshold_for` is
computed from *that*.

**The STR scope is the analysed regions unmasked.** A tract in the catalog is a stretch of real
sequence by construction, and cutting the scope at every `N` run would split scope spans without
changing which tracts they hold.

### 2.2 The generic rule, exposed on its own

The hash rule is separated from the sweep over regions so a test can ask about one position without
building a genome's worth (spec §7.1).

**Built**, in `loci.rs`.

```rust
/// One contig's name absorbed once, so the per-position work is a clone and one `u64`.
/// **The value is identical to `hash_position`'s**, which its own test pins — the whole
/// selection rests on those two agreeing.
pub struct ContigHasher { /* … */ }
impl ContigHasher {
    pub fn new(contig_name: &str, seed: u64) -> Self;
    pub fn hash_position(&self, position: u64) -> u64;
}
pub fn hash_position(contig_name: &str, position: u64, seed: u64) -> u64;

/// Keep this position? `hash(contig, position, seed) < threshold` (spec §2).
pub fn keeps_position(hasher: &ContigHasher, position: Position, threshold: u64) -> bool;

/// The threshold that yields about `target` positions out of `selectable_length`.
///
/// **Computed in 128 bits**: `2^64 · target` overflows a `u64`. Saturates to "keep
/// everything" when `target >= selectable_length`, which is the right answer rather than
/// a case to guard (spec §2).
pub fn threshold_for(target: u64, selectable_length: Bp) -> u64;

/// Every position the rule keeps, ascending. `contig_names` is indexed by `ContigId`,
/// as the reference's contig table is; a span naming a contig the table does not hold
/// is an error and never a skip.
pub fn select_generic_positions(
    regions: &SelectableRegions,
    contig_names: &[String],
    seed: u64,
    threshold: u64,
) -> Result<Vec<GenomePosition>, SelectionError>;
```

**Decided — the same hash construction the catalog's sampler uses**, keyed by the contig's **name**
as `hash_locus` is ([`src/ng/repeat_catalog/strata.rs:158`](../../../../src/ng/repeat_catalog/strata.rs)),
so both halves of the selection rest on one uniformity assumption rather than two (spec §2). *An
earlier draft of this section proposed lifting `hash_locus` to `pub(crate)` and rekeying it by
`ContigId`. Neither happened, and neither should: an index is a property of the file's contig order,
where a name is a property of the genome, and the per-position cost that motivated the index is
removed by `ContigHasher` instead — the name is hashed once per contig and the hasher cloned per
position.*

**Measured, on real references rather than a fixture** (`examples/ng_joint_loci_probe.rs`): tomato
returns 2,002,505 positions for a target of 2,000,000 and GRCh38 returns 1,999,981; the gaps are
geometric to three figures at both; the per-contig counts give chi-square 19.8 on 12 degrees of
freedom and 127.9 on 127; and a sweep over shards cut at coordinates the rule knows nothing about
returns an identical vector.

### 2.3 Errors

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum SelectionError {
    /// A selectable span ends before it starts, or two of them overlap. **Built** — a
    /// position offered twice is a silent double weight and nothing downstream sees it.
    InvertedRegion { region: GenomeRegion },
    OverlappingRegions { first: GenomeRegion, second: GenomeRegion },
    /// A span names a contig the reference table does not hold. **Built** — skipping it
    /// silently drops a whole chromosome from the selection, identically in every
    /// sample, so no cross-sample check could notice.
    UnknownContig { contig: ContigId },
    /// The run asked the catalog for tracts it was not built to hold. Carries the axis
    /// and both values; the run stops rather than proceeding on a short list, which
    /// would be a wrong per-stratum total nothing downstream could notice (spec §3.3).
    #[error("the repeat catalog cannot serve this run's criteria")]
    CatalogTooRestrictive(#[source] RepeatCatalogError),
    /// Two samples' [`SelectionIdentity`] disagree. Names the field.
    #[error("samples were selected under different {field}; they hold different loci")]
    IdentityMismatch { field: &'static str },
    /// The digests of two samples' kept loci differ although their identities agree —
    /// a hash function or a threshold's arithmetic changed underneath them (§1.3).
    #[error("two samples' kept loci differ in the megabase block starting at {block_start}")]
    DigestMismatch { block_start: GenomePosition },
}
```

**Fatal, never a warning.** Both mismatch variants stop the fit; averaging over mismatched sets
produces a number with no meaning (spec §5).

---

## 3. Design decisions — decided

- **`KeptLoci` holds the loci; `SelectionIdentity` travels.** The loci are reproducible and
  coordinates cost 8× the records they index — spec §5 (census sites §5).
- **The STR half is the catalog's, not this unit's.** `sample_loci_per_stratum` ships with the
  ordering, exactness and refusal properties this consumer needs, asserted in its own spec §10 —
  spec §3.2. *No fork: the alternative rules were rejected in the spec before the sampler was built.*
- **The per-stratum counts are a stored field of the selection, not a note.** Reweighting is one line
  of arithmetic and silent when omitted — spec §3.5.
- **The digest is fed by the writer, not derived from the rule.** A re-derived digest proves
  determinism, not agreement — spec §5.1.
- **`PartialEq` is the compatibility check.** No partial matches, because the failure is
  meaninglessness rather than noise — spec §5.
- **OPEN:** the generic budget's default. Spec §6 question 3 has the experiment; the type takes a
  `u64` either way, so nothing here waits on it. **What 2026-08-13 settles is which way to be wrong**:
  a run writing records to a file writes the largest census it will plausibly want, because a smaller
  target's positions are a subset of a larger one's under `keeps_position` and can be taken without
  reading anything again, while a larger one needs the records rebuilt — spec §4.3.4.

---

## 4. Reconciliation with existing code

| this doc | existing code | how they meet |
|---|---|---|
| the STR sampler | `RepeatCatalog::sample_loci_per_stratum` ([`repeat_catalog/reader.rs:260`](../../../../src/ng/repeat_catalog/reader.rs)) | called, not reimplemented; returns `(StratumCounts, StratumSample)` from one pass, which is exactly `KeptLoci`'s two STR fields |
| the per-stratum totals | `RepeatCatalog::count_loci_per_stratum` ([`repeat_catalog/reader.rs:242`](../../../../src/ng/repeat_catalog/reader.rs)) | it delegates to the sampler at `cap = 0`, so asking for both costs one pass — take the pair from `sample_loci_per_stratum` and never call this one separately |
| `StratumCounts`, `StratumSample` | [`repeat_catalog/strata.rs:15,49`](../../../../src/ng/repeat_catalog/strata.rs) | used as-is; `KeptLoci` stores them rather than converting |
| the hash | `hash_locus` ([`repeat_catalog/strata.rs:158`](../../../../src/ng/repeat_catalog/strata.rs), called from `StratumSampler::offer:110`) | the generic rule uses the same function and seed convention. **It is private today and keyed by contig *name*** — lifting it to `pub(crate)` and taking a `ContigId` is this unit's one edit outside its own file |
| `StrRepeatCriteria` | [`repeat_catalog/criteria.rs`](../../../../src/ng/repeat_catalog/criteria.rs), taken by `str_loci` ([`reader.rs:220`](../../../../src/ng/repeat_catalog/reader.rs)) | stored in `SelectionIdentity` by value; it is already the run's whole STR policy |
| the catalog-vs-reference check | `RepeatCatalog::open_checking_against_reference` ([`repeat_catalog/reader.rs:53`](../../../../src/ng/repeat_catalog/reader.rs)) | happens before this unit runs, which is why the catalog's own identity is **not** in `SelectionIdentity` (spec §5.1) |
| `GenomePosition`, `ContigId`, `Position`, `Bp` | [`src/ng/types.rs:60,13,34,174`](../../../../src/ng/types.rs) | used as-is; this unit seeds no new scalar |
| the scope of a run | `ReadScope` (taken by `str_loci`/`count_loci_per_stratum`, [`reader.rs:220,242`](../../../../src/ng/repeat_catalog/reader.rs)) | `AnalysedRegions` must convert into one, so a `--regions` run counts and samples inside the BED — spec §3.3 |

---

## 5. Open items

- **`OPEN:` how the analysed region set is represented** (`AnalysedRegions` above). Step 4 has no
  region type of its own today and the catalog takes `ReadScope`; pin it against whatever the driver
  passes when the driver exists.
- **Impl-time:** the digest's hash function and width. 32 bytes and per-megabase blocks are the
  contract; which construction fills them is a code-review question.
- **Impl-time:** whether `KeptLoci::generic` is a `Vec<GenomePosition>` or a per-contig run-length
  form. The contract is *ascending, indexable*; the memory shape is the implementer's.

---

## 6. Test and bench shape

Tests live beside the code in `joint/loci.rs`'s `#[cfg(test)] mod tests`. **The regression anchor is
identity rather than a fixture**: build the selection twice, once whole and once from region shards
at several thread counts, and compare — spec §7.1, which no fixture file is needed for. The
`--regions` case is the one most likely to be written against the contig table by reflex, so it is
asserted by name. The digest's own test plants a wrong position in an otherwise correct walk and
requires the block digest to name the megabase, then requires a rule-derived digest to **pass** that
same test — the check that the check is the right one (spec §7.8).
