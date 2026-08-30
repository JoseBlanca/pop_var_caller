//! Step 4's front door on the SNP/indel path: the two tables a sample's loci are
//! poured into, and everything the pour had to do to a locus other than enter it as it
//! arrived.
//!
//! **Two objects, differing only in how a site is keyed.** The read-group table enters
//! a site once per read group that covered it — a site with 20 reads from one library
//! and 10 from another becomes two entries, at depth 20 and at depth 10 — because an
//! error rate describes the **chemistry**, and two libraries of one sample can
//! genuinely differ. The windowed table enters that same site **once, at depth 30**,
//! because heterozygosity describes the **individual**, and one genome has one
//! heterozygosity however many libraries were used to read it.
//!
//! **Neither can stand in for the other.** Splitting a site's reads across two entries
//! costs the error rate nothing, since what it counts is reads. That same split is
//! fatal to the other rate: a site covered by two libraries has become two shallower
//! entries and its total depth is gone, so heterozygosity cannot be counted there at
//! all. Both are built, and neither is derived from the other
//! (`spec/parameter_prepass_generic.md` §1, §5.1).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §3.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use crate::ng::locus_generation::{LocusKind, SampleLocusObservations};
use crate::ng::types::{Bp, ContigId, GenomeRegion, InbreedingF, Ploidy, ReadGroupId};

use super::calibration::{MintedReadErrors, fold_into, minted_error_by_read_group};
use super::depth_and_alt_reads::{
    CountedSite, count_by_read_group, count_whole_site, count_whole_site_by_library,
};
use super::depth_bins::DepthBinEdges;
use super::histogram::{DepthAltHistogram, fold_windows_of_one_ploidy};
use super::{INBREEDING_WINDOW_BP, WindowIndex};

/// Which set of genotypes a region's sites are drawn from.
///
/// **Part of both tables' keys**, because a haploid region has two genotype classes and
/// a diploid three: pooling them would score a site against the wrong set. One error
/// rate is still fitted per read group across all of them — chemistry does not know
/// about chromosomes — but each cell is scored at its own ploidy, and the genotype
/// frequencies come out one set per ploidy.
///
/// **What this module needs is only a lookup**, `region -> Ploidy`. Where that comes
/// from — a flag, a BED, a per-contig default — is not settled and is not this module's
/// business (`arch/parameter_prepass_generic.md` §3).
pub trait PloidyMap: Send + Sync {
    fn ploidy_at(&self, region: GenomeRegion) -> Ploidy;
}

/// One ploidy everywhere — **today's behaviour exactly**, and the only implementation
/// this plan builds.
#[derive(Copy, Clone, Debug)]
pub struct ConstantPloidy(pub Ploidy);

impl PloidyMap for ConstantPloidy {
    fn ploidy_at(&self, _region: GenomeRegion) -> Ploidy {
        self.0
    }
}

/// Whether this run fits the inbreeding coefficient or was handed one.
///
/// **It changes what is accumulated, which is why it is a constructor argument and not
/// a flag read at the end.** `Fitted` keys the windowed table by window, as the runs
/// model needs. `Supplied` drops the window key — the object itself stays, because the
/// per-site rates still need whole sites — collapsing it from about 37 MB to a few kB
/// per sample (`spec/parameter_prepass_generic.md` §6.4).
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum InbreedingMode {
    Fitted,
    Supplied(InbreedingF),
}

/// Which window a site's evidence is filed under.
///
/// **Two arms, and which one a run uses is [`InbreedingMode`]'s to decide.** A run
/// fitting `F` needs the genome cut into 100 kb windows, because the runs model reads
/// them in order along the chromosome; a run handed `F` needs none of that and keeps
/// one table for the whole sample. Ordering is by contig and then by window, so a
/// `BTreeMap` keyed on this walks the genome in order for free.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum WindowKey {
    /// One 100 kb window of one contig.
    InWindow {
        contig: ContigId,
        window: WindowIndex,
    },
    /// The whole sample in one table — `F` was supplied, so no window key is kept.
    WholeSample,
}

impl WindowKey {
    /// Which window a locus starts in.
    ///
    /// **ng's [`Position`] is 1-based, and this is where that is dealt with.** A naive
    /// `start / INBREEDING_WINDOW_BP` would put positions 1 to 99,999 in window 0 and
    /// give every later window a full 100,000 — a first window one base short, and a
    /// resolution that differs between the first window of a contig and all the rest.
    /// Subtracting one first gives every window exactly 100,000 positions, window 0
    /// holding 1 to 100,000.
    ///
    /// A locus is filed by where it **starts**, so a locus spanning a window boundary —
    /// an indel widened across one — belongs wholly to the window it opens in. That is
    /// a rounding of at most a few bases against a 100,000-base window, and the
    /// alternative, splitting a locus's evidence across two windows, would break "a
    /// site enters once, whole".
    fn of_locus(region: GenomeRegion) -> Self {
        let from_contig_start = region.start.get().saturating_sub(1);
        let window = from_contig_start / INBREEDING_WINDOW_BP.get();
        Self::InWindow {
            contig: region.contig,
            window: WindowIndex(u32::try_from(window).unwrap_or(u32::MAX)),
        }
    }
}

/// **Everything the accumulator did to a locus other than enter it as it arrived.**
///
/// Every field is a plain sum, so it merges with the rest of the accumulator; every
/// field is reported, because the alternative is a fit that quietly describes a
/// different population of reads than the caller will see
/// (`spec/parameter_prepass.md` §2).
#[derive(Clone, Default, PartialEq, Eq, Debug)]
pub struct AccumulationCounts {
    /// Loci whose reads a cap **upstream** had already subsampled. Entered anyway, at
    /// the depth observed: a site subsampled to 250 reads is a site observed at depth
    /// 250, and skipping it would be the depth-dependent selection this step exists to
    /// remove. A run where this is a large share of the loci is one whose error rate
    /// describes the ends of reads more than it describes the chemistry.
    pub loci_with_upstream_subsample: u64,
    /// Reads that covered a locus and witnessed nothing, summed over loci. Not part of
    /// any depth, and a locus with many of them is one whose depth is not what it looks.
    pub reads_without_observation: u64,
    /// Sites **this** step subsampled down to the ladder's top depth. A run where most
    /// sites appear here is one where the ladder, not the data, is setting the depth.
    pub sites_subsampled_to_cap: u64,
    /// Loci that began before the previous locus on their contig ended. **Must be
    /// zero** — see [`GenericAccumulators::add_locus`].
    pub loci_overlapping_previous: u64,
    /// Shard spans that overlap the span before them on a contig. **Must be zero**, and
    /// **a lower bound rather than a count** when it is not: the check walks the sorted
    /// spans in pairs, so three shards all covering one stretch report 1 rather than 3.
    /// Zero-detection is exact, which is the property that matters — and while it is
    /// zero,
    /// `loci_overlapping_previous` is the exact count a single unsharded walk would have
    /// produced: no locus of one shard can overlap a locus of another if the shards'
    /// spans are disjoint.
    ///
    /// It is a separate counter rather than folded into the one above because it counts
    /// a different thing — a seam, not a locus. Counting the *loci* that overlap across
    /// a seam would need the loci, and a shard keeps only its span.
    pub shard_spans_overlapping: u64,
}

/// The stretch of one contig one shard's loci covered.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct ContigSpan {
    contig: ContigId,
    first_start: u64,
    last_end: u64,
}

/// Step 4's generic front door: holds both tables, and adds one locus to both.
///
/// **Holds the run's binning rule and hands a handle to every table it creates**, so
/// that all of a worker's tables — and every other worker's — are binned by one shared
/// object. Cloning the handle is an atomic increment paid once per table; binning a
/// site is a dereference and costs no atomic at all.
///
/// **One of these per region shard, merged at the end** — the shape rayon's
/// `fold`/`reduce` consumes. No trait and no adapter: the accumulator is a plain owned
/// struct, and the "passes the locus on untouched" property is the `&` on
/// [`GenericAccumulators::add_locus`], not an iterator wrapper.
pub struct GenericAccumulators {
    edges: Arc<DepthBinEdges>,
    /// Whether this sample has more than one library, which is what decides whether a
    /// site's alternative reads keep the library they came from. **A property of the
    /// sample, not of the site** — at one read group an attributed key carries nothing
    /// a pooled one does not, and entering such a sample through the attributed door
    /// would build a sparse map that `spec/parameter_prepass_generic.md` §9 does not
    /// price and would break §12.6's cell-for-cell equality between the two tables.
    multi_library: bool,
    ploidy: Arc<dyn PloidyMap>,
    inbreeding: InbreedingMode,
    /// The read-group table. **`u64`, because it accumulates over the whole sample and
    /// no fold widens it**: on a human sample its busiest cell's depth sum passes
    /// 4.29 × 10⁹ about a third of the way through the run (owner's call, 2026-08-06).
    by_group: BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
    /// The windowed table, used when `F` is being fitted. `u32`: 100 kb of sites can reach
    /// neither ceiling, and a sample holds about 8,000 of these, which is where the narrower
    /// counter buys something — 37 MB rather than 74.
    ///
    /// **Empty when `F` was supplied**, because the collapsed arm is not a window and the
    /// `u32` justification above does not survive the collapse. See [`Self::whole_sample`].
    by_window: BTreeMap<(WindowKey, Ploidy), DepthAltHistogram<u32>>,
    /// The collapsed table, used when `F` was **supplied** — one per ploidy, no window key.
    ///
    /// **`u64`, and the width is the whole reason this is a second field rather than another
    /// key in `by_window`.** That map's `u32` is justified by a window holding at most
    /// 100,000 sites; dropping the window key pours the entire genome into one table, and
    /// tomato's ~800 Mb at 30× puts well over 300 million sites into the modal-depth
    /// zero-alternative cell. Measured before this split: a panic at **143,165,576 sites at
    /// depth 30 in one cell** — loud rather than wrapping, which is the counter's overflow
    /// guard working, but every supplied-`F` run on a real genome died there.
    ///
    /// **Widening costs nothing here and would have cost 37 MB there.** This map holds one
    /// table per ploidy — a few kB — where `by_window` holds about 8,000 of them, which is
    /// the memory `InbreedingMode::Supplied` exists to save in the first place
    /// (`spec/parameter_prepass_generic.md` §6.4, §9).
    ///
    /// **The two are never both populated**: [`Self::add_locus`] routes by
    /// [`InbreedingMode`], and `merge` refuses two accumulators whose modes differ.
    whole_sample: BTreeMap<Ploidy, DepthAltHistogram<u64>>,
    counts: AccumulationCounts,
    /// One entry per contig this shard saw. Bounded by shards × contigs — tens to
    /// hundreds of pairs — so nothing is saved by collapsing it.
    covered_spans: Vec<ContigSpan>,
    /// Where the previous locus on each contig ended. **Sequential state, and
    /// deliberately not merged**: two shards that split a contig each begin with no
    /// previous locus, so neither ever checks the seam between them. The seam is the
    /// span list's business.
    previous_end: BTreeMap<ContigId, u64>,
    by_group_scratch: Vec<(ReadGroupId, CountedSite)>,
    alt_by_group_scratch: Vec<(ReadGroupId, u32)>,
    /// **The read likelihood's error-rate scale, denominator half** — Σ `q_sum` and the read
    /// count per read group, over exactly the sites [`Self::read_group_histograms`] counted
    /// (`spec/read_likelihoods.md` §3.2, and [`calibration`](super::calibration) for what the
    /// average is and why).
    ///
    /// **Beside the histograms rather than inside them**, because nothing here is fitted from it:
    /// the fit reads the tables and this is divided into the fit's answer afterwards. One entry per
    /// read group the sample covered — bytes.
    minted_errors: BTreeMap<ReadGroupId, MintedReadErrors>,
    minted_errors_scratch: Vec<(ReadGroupId, MintedReadErrors)>,
}

impl GenericAccumulators {
    /// One per region shard.
    ///
    /// **Every shard must be handed the same `edges`**, or their cells mean different
    /// things and `merge` is meaningless — which is the whole reason the rule is shared
    /// rather than rebuilt per worker, and why `merge` proves it by pointer identity.
    ///
    /// `read_groups` is the sample's own set. At one read group no site ever takes the
    /// attributed arm, so 1,550 of the 1,707 samples in the tomato archive survey pay
    /// nothing for the multi-library machinery.
    #[must_use]
    pub fn new(
        edges: Arc<DepthBinEdges>,
        read_groups: &[ReadGroupId],
        ploidy: Arc<dyn PloidyMap>,
        inbreeding: InbreedingMode,
    ) -> Self {
        Self {
            edges,
            multi_library: read_groups.len() > 1,
            ploidy,
            inbreeding,
            by_group: BTreeMap::new(),
            by_window: BTreeMap::new(),
            whole_sample: BTreeMap::new(),
            counts: AccumulationCounts::default(),
            covered_spans: Vec::new(),
            previous_end: BTreeMap::new(),
            by_group_scratch: Vec::new(),
            alt_by_group_scratch: Vec::new(),
            minted_errors: BTreeMap::new(),
            minted_errors_scratch: Vec::new(),
        }
    }

    /// Add one locus to both tables — once per read group that covered it in the first,
    /// once at its total depth in the second.
    ///
    /// **Borrows**: the caller keeps the locus and passes it on unchanged. A locus whose
    /// kind is not [`LocusKind::Generic`] is ignored, because the STR path has its own
    /// accumulator and its own stutter model.
    ///
    /// **It tallies rather than repairs.** The invariant it depends on is that generic
    /// loci *partition* the covered positions: a locus can span several reference bases
    /// — the pileup widens a record to an indel's reference span — but no two loci may
    /// cover the same position, or a site enters the windowed table twice and "a site
    /// enters once, whole" is broken. So this keeps the previous locus's end per contig,
    /// one comparison per locus, and counts any locus that starts before it.
    ///
    /// Counting rather than asserting, because a `debug_assert!` compiles out of the
    /// release build this repo actually runs — a trap it has recorded hitting twice.
    /// Counting rather than de-duplicating, because a de-duplication rule would hide an
    /// upstream bug behind a plausible number, and a counter that must read zero cannot.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations) {
        if !matches!(locus.kind, LocusKind::Generic) {
            return;
        }

        let region = locus.region;
        let ploidy = self.ploidy.ploidy_at(region);
        let covered = Bp(region.len());
        self.note_position(region);

        if locus.reads_discarded_by_cap > 0 {
            self.counts.loci_with_upstream_subsample += 1;
        }
        self.counts.reads_without_observation += u64::from(locus.reads_without_observation);

        // The scratch buffers are lifted out so that the counting functions can borrow
        // `self.edges` while filling them, and put back so the next locus reuses the
        // allocation.
        let edges = Arc::clone(&self.edges);
        let mut by_group_scratch = std::mem::take(&mut self.by_group_scratch);
        count_by_read_group(locus, &edges, &mut by_group_scratch);
        for &(group, counted) in &by_group_scratch {
            self.by_group
                .entry((group, ploidy))
                .or_insert_with(|| DepthAltHistogram::new(Arc::clone(&edges)))
                .add_site(counted.counts(), covered);
        }
        self.by_group_scratch = by_group_scratch;

        // **The scale's denominator, over the sites the table above just counted.** Folded here
        // rather than beside the fit because this is the only place that sees a read's minted
        // error at all — by the time the fit runs, the observations are cells.
        let mut minted = std::mem::take(&mut self.minted_errors_scratch);
        minted_error_by_read_group(locus, &mut minted);
        fold_into(&mut self.minted_errors, &minted);
        self.minted_errors_scratch = minted;

        // **The site is counted once and filed once, into whichever of the two tables this
        // run's mode uses.** The counting is identical; only the width of what receives it
        // differs, which is why the two arms below duplicate a call rather than a rule.
        let mut alt_by_group = std::mem::take(&mut self.alt_by_group_scratch);
        let whole = if self.multi_library {
            count_whole_site_by_library(locus, &edges, &mut alt_by_group)
        } else {
            alt_by_group.clear();
            count_whole_site(locus, &edges)
        };
        let multi_library = self.multi_library;
        match self.inbreeding {
            InbreedingMode::Fitted => {
                let table = self
                    .by_window
                    .entry((WindowKey::of_locus(region), ploidy))
                    .or_insert_with(|| DepthAltHistogram::new(Arc::clone(&edges)));
                if multi_library {
                    table.add_attributed_site(whole.counts(), &alt_by_group, covered);
                } else {
                    table.add_site(whole.counts(), covered);
                }
            }
            InbreedingMode::Supplied(_) => {
                let table = self
                    .whole_sample
                    .entry(ploidy)
                    .or_insert_with(|| DepthAltHistogram::new(Arc::clone(&edges)));
                if multi_library {
                    table.add_attributed_site(whole.counts(), &alt_by_group, covered);
                } else {
                    table.add_site(whole.counts(), covered);
                }
            }
        }
        self.alt_by_group_scratch = alt_by_group;
        if whole.subsampled_from().is_some() {
            self.counts.sites_subsampled_to_cap += 1;
        }
    }

    /// Note where this locus sits, for the two overlap checks.
    fn note_position(&mut self, region: GenomeRegion) {
        let (start, end) = (region.start.get(), region.end.get());
        match self.previous_end.get_mut(&region.contig) {
            Some(previous) => {
                if start <= *previous {
                    self.counts.loci_overlapping_previous += 1;
                }
                *previous = (*previous).max(end);
            }
            None => {
                self.previous_end.insert(region.contig, end);
            }
        }

        match self
            .covered_spans
            .iter_mut()
            .find(|span| span.contig == region.contig)
        {
            Some(span) => {
                span.first_start = span.first_start.min(start);
                span.last_end = span.last_end.max(end);
            }
            None => self.covered_spans.push(ContigSpan {
                contig: region.contig,
                first_start: start,
                last_end: end,
            }),
        }
    }

    /// Combine a shard's accumulators into these. Associative and exact.
    ///
    /// # Panics
    ///
    /// **If the two shards were built with anything different that changes what their
    /// cells mean — all four of `new`'s arguments.** Different edges objects or ploidy
    /// maps, checked by pointer identity; different library counts or inbreeding modes,
    /// checked by value. Every one of those mismatches is silent if it merges: two
    /// shards handed different ploidy maps produce cells scored against different sets
    /// of genotype classes and nothing downstream can tell, which is this module's
    /// declared failure mode on the argument that was, until this check, the one left
    /// unguarded.
    pub fn merge(&mut self, other: Self) {
        assert!(
            Arc::ptr_eq(&self.edges, &other.edges),
            "these two shards were binned by different ladder objects, so their cells \
             are not the same cells"
        );
        assert_eq!(
            self.multi_library, other.multi_library,
            "these two shards disagree about whether the sample has more than one \
             library, so one keyed its sites by library and the other did not"
        );
        assert_eq!(
            self.inbreeding, other.inbreeding,
            "these two shards disagree about whether F is fitted, so one keyed its \
             sites by window and the other did not"
        );
        assert!(
            Arc::ptr_eq(&self.ploidy, &other.ploidy),
            "these two shards were handed different ploidy maps, so their cells were \
             scored against different sets of genotype classes — a haploid region has \
             two and a diploid three, and pooling them scores a site against the wrong \
             set"
        );

        for (key, table) in other.by_group {
            match self.by_group.get_mut(&key) {
                Some(mine) => mine.merge(&table),
                None => {
                    self.by_group.insert(key, table);
                }
            }
        }
        for (ploidy, table) in other.whole_sample {
            match self.whole_sample.get_mut(&ploidy) {
                Some(existing) => existing.merge(&table),
                None => {
                    self.whole_sample.insert(ploidy, table);
                }
            }
        }
        // **The scale's denominator merges like the tables and for the same reason**: a shard saw
        // some of the sample's sites and the average has to be over all of them. Both of its
        // numbers are integers, so this merge is exact and
        // `three_shards_merged_in_every_order_give_the_same_counters` (below, in this file's test
        // module) covers it rather than passing over it — see [`MintedReadErrors`] for why an
        // `f64` sum here would have broken that contract silently.
        for (group, totals) in other.minted_errors {
            self.minted_errors.entry(group).or_default().add(totals);
        }
        for (key, table) in other.by_window {
            match self.by_window.get_mut(&key) {
                Some(mine) => mine.merge(&table),
                None => {
                    self.by_window.insert(key, table);
                }
            }
        }

        // Destructured exhaustively on purpose: a counter added later stops this
        // compiling rather than being silently left out of every merge. All four of
        // these sums could be deleted with the suite green until
        // `merge_carries_every_counter_and_both_tables` existed.
        let AccumulationCounts {
            loci_with_upstream_subsample,
            reads_without_observation,
            sites_subsampled_to_cap,
            loci_overlapping_previous,
            // Never stored — `adjustments` derives it from the spans, which merge below.
            shard_spans_overlapping: _,
        } = other.counts;
        let counts = &mut self.counts;
        counts.loci_with_upstream_subsample += loci_with_upstream_subsample;
        counts.reads_without_observation += reads_without_observation;
        counts.sites_subsampled_to_cap += sites_subsampled_to_cap;
        counts.loci_overlapping_previous += loci_overlapping_previous;
        self.covered_spans.extend(other.covered_spans);
    }

    /// Everything the accumulator did to a locus other than enter it as it arrived, plus
    /// the two counters that must read zero.
    ///
    /// **Owned rather than borrowed, because the seam count is derived here.** The
    /// shard spans are checked for overlap once, at read time, by sorting — not during
    /// `merge`. Collapsing them to a first-start and last-end pair per contig and
    /// checking the seam as shards arrive reports an overlap that does not exist:
    /// merging three contiguous shards `[0,100)`, `[200,300)`, `[100,200)` in that order
    /// gives `[0,300)` after the first two, and the third then reads as overlapping. **A
    /// counter that must be zero may not have false positives**, and the list is bounded
    /// by shards × contigs, so nothing is saved by collapsing it.
    #[must_use]
    pub fn adjustments(&self) -> AccumulationCounts {
        let mut spans = self.covered_spans.clone();
        spans.sort_unstable();
        let overlapping = spans
            .windows(2)
            .filter(|pair| {
                pair[0].contig == pair[1].contig && pair[1].first_start <= pair[0].last_end
            })
            .count();

        AccumulationCounts {
            shard_spans_overlapping: overlapping as u64,
            ..self.counts.clone()
        }
    }

    /// The read-group table: one cell table per `(read group, ploidy)`, a site entering
    /// once per group that covered it. Kilobytes.
    #[must_use]
    pub fn read_group_histograms(
        &self,
    ) -> &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>> {
        &self.by_group
    }

    /// The windowed table: one cell table per `(window, ploidy)`, a site entering once
    /// at its total depth. In genome order, which is what the runs model reads it in.
    #[must_use]
    pub fn windowed_histograms(&self) -> &BTreeMap<(WindowKey, Ploidy), DepthAltHistogram<u32>> {
        &self.by_window
    }

    /// **The read likelihood's error-rate scale, denominator half** — per read group, over
    /// exactly the sites [`Self::read_group_histograms`] counted.
    ///
    /// The numerator is the error rate
    /// [`fit_read_group_error_rates`](super::read_group_error_rate::fit_read_group_error_rates)
    /// fits from those same tables, and the scale is the one divided by the other
    /// (`spec/read_likelihoods.md` §3.2). **Nothing here is fitted from this**, which is why it
    /// sits beside the tables rather than inside them and why the fit does not take it.
    ///
    /// A read group the sample covered with no complete observation anywhere has no entry, which
    /// is the difference between a library that was measured and one that was not — a scale
    /// wanting a denominator gets no number rather than a one.
    #[must_use]
    pub fn minted_errors(&self) -> &BTreeMap<ReadGroupId, MintedReadErrors> {
        &self.minted_errors
    }

    /// Sum this sample's windows into one whole-sample table, for one ploidy — free and
    /// exact, which is why no third accumulator is kept.
    #[must_use]
    pub fn whole_sample_histogram(&self, ploidy: Ploidy) -> DepthAltHistogram<u64> {
        match self.inbreeding {
            InbreedingMode::Fitted => fold_windows_of_one_ploidy(
                &self.edges,
                self.by_window
                    .iter()
                    .filter(|((_, at), _)| *at == ploidy)
                    .map(|(_, table)| table),
            ),
            // **A fold over nought or one table, and it is still a fold.** The collapsed
            // arm is already whole-sample and already `u64`, so this copies rather than
            // widens — but going through the same function is what keeps one definition of
            // what "the whole sample at this ploidy" means. `Option` is an iterator of at
            // most one item, so a ploidy this sample never saw folds to an empty table,
            // which is the same answer the arm above gives.
            InbreedingMode::Supplied(_) => {
                fold_windows_of_one_ploidy(&self.edges, self.whole_sample.get(&ploidy))
            }
        }
    }

    /// Every ploidy this sample's loci landed at.
    ///
    /// **Its own accessor because the two modes keep their sites in different maps**, and a
    /// caller deriving the set from `windowed_histograms().keys()` — which is what the
    /// coupled fit did — sees nothing at all when `F` was supplied.
    #[must_use]
    pub fn ploidies(&self) -> BTreeSet<Ploidy> {
        match self.inbreeding {
            InbreedingMode::Fitted => self.by_window.keys().map(|&(_, ploidy)| ploidy).collect(),
            InbreedingMode::Supplied(_) => self.whole_sample.keys().copied().collect(),
        }
    }

    /// What this run was told about inbreeding, which decides whether the windowed table
    /// is keyed by window at all.
    #[must_use]
    pub fn inbreeding_mode(&self) -> InbreedingMode {
        self.inbreeding
    }

    /// **Reference positions this sample's loci covered, over every ploidy** — the warrant a
    /// supplied inbreeding coefficient carries, since it says how much genome the number is
    /// being applied to.
    ///
    /// **From the windowed table and not the read-group one**, because a site covered by two
    /// libraries enters the read-group table twice and its positions would be counted twice
    /// with it. The windowed table enters each site once, whole, which is the whole reason
    /// both tables exist.
    #[must_use]
    pub fn covered_positions(&self) -> u64 {
        self.by_window
            .values()
            .map(DepthAltHistogram::total_covered_positions)
            .sum::<u64>()
            + self
                .whole_sample
                .values()
                .map(DepthAltHistogram::total_covered_positions)
                .sum::<u64>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{ReadWitness, SequenceObservation};
    use crate::ng::types::Position;

    fn ladder() -> Arc<DepthBinEdges> {
        Arc::new(DepthBinEdges::new())
    }

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("a positive copy number")
    }

    fn group(id: u32) -> ReadGroupId {
        ReadGroupId(id)
    }

    fn observation(bases: &[u8], read_group: ReadGroupId, num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.into(),
            read_witness: ReadWitness::Complete,
            read_group,
            num_obs,
            num_fwd: 0,
            // **Awkward on purpose, and it used to be `0.0`.** The calibration totals
            // (`Self::minted_errors`) are summed from this field, and with every fixture at zero
            // there was nothing for a wrong sum to lose: a fold routed through `f64` passed the
            // six-order merge test, which is the one test whose whole subject is that merging is
            // order-independent. These decimals do not add associatively as floats.
            q_sum: crate::ng::types::SummedLogError::from_nats(
                -0.1 * f64::from(num_obs) - 0.07 * f64::from(read_group.get() + 1),
            ),
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn site(
        contig: u32,
        start: u64,
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(contig),
                start: Position(start),
                end: Position(start),
            },
            reference_bases: b"A".as_slice().into(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// The run's one ploidy map, shared by every shard exactly as a real run shares it.
    ///
    /// **A `thread_local` rather than a fresh `Arc` per shard**, because `merge` proves
    /// two shards were handed the same map by pointer identity: minting one per shard
    /// here would have made the merge tests exercise a configuration no run can produce,
    /// and they failed the moment that check was added — which is the check working.
    fn shared_ploidy() -> Arc<dyn PloidyMap> {
        thread_local! {
            static PLOIDY: Arc<dyn PloidyMap> = Arc::new(ConstantPloidy(
                Ploidy::try_new(2).expect("a positive copy number"),
            ));
        }
        PLOIDY.with(Arc::clone)
    }

    fn accumulators(edges: &Arc<DepthBinEdges>, groups: &[ReadGroupId]) -> GenericAccumulators {
        GenericAccumulators::new(
            Arc::clone(edges),
            groups,
            shared_ploidy(),
            InbreedingMode::Fitted,
        )
    }

    /// A sample's sites in genome order, over several windows and both alleles.
    ///
    /// **One contig, deliberately.** An earlier version alternated two contigs, and
    /// that made the sharding test below unable to fail: shards chunked out of it never
    /// had two spans of the *same* contig adjacent in the merged list, so dropping the
    /// sort in `adjustments` — the exact fault the architecture warns about — left the
    /// suite green. A shard is a stretch of one contig; the fixture should look like
    /// one.
    fn a_walk() -> Vec<SampleLocusObservations> {
        (0..12u32)
            .map(|at| {
                let start = 1 + u64::from(at) * 40_000;
                site(
                    0,
                    start,
                    vec![
                        observation(b"A", group(0), 20 - at),
                        observation(b"T", group(0), at),
                    ],
                )
            })
            .collect()
    }

    /// **A non-generic locus changes nothing.** The STR path has its own accumulator and
    /// its own stutter model; a tract entering this one would be scored against a
    /// per-base error rate that does not describe it.
    #[test]
    fn a_locus_that_is_not_generic_changes_nothing() {
        let edges = ladder();
        let mut accumulators = accumulators(&edges, &[group(0)]);

        let mut tract = site(0, 500, vec![observation(b"A", group(0), 10)]);
        tract.kind = LocusKind::SsrBundle;
        accumulators.add_locus(&tract);

        assert!(accumulators.read_group_histograms().is_empty());
        assert!(accumulators.windowed_histograms().is_empty());
        assert_eq!(accumulators.adjustments(), AccumulationCounts::default());
    }

    /// **The two tables hold the same sites, keyed two ways.** On a single-library
    /// sample they hold the same *numbers* too, which is the identity
    /// `spec/parameter_prepass_generic.md` §12.6 pins: the read-group table equals the
    /// windowed one folded over its windows, cell for cell.
    #[test]
    fn on_a_single_library_sample_the_two_tables_agree_cell_for_cell() {
        let edges = ladder();
        let mut accumulators = accumulators(&edges, &[group(0)]);
        for locus in a_walk() {
            accumulators.add_locus(&locus);
        }

        let by_group = accumulators.read_group_histograms();
        assert_eq!(by_group.len(), 1, "one library, one ploidy");
        let read_group_table = &by_group[&(group(0), diploid())];
        let folded = accumulators.whole_sample_histogram(diploid());

        assert_eq!(read_group_table.cells(diploid()), folded.cells(diploid()));
        assert_eq!(read_group_table.total_loci(), folded.total_loci());
        assert_eq!(
            read_group_table.total_covered_positions(),
            folded.total_covered_positions()
        );
        assert_eq!(folded.total_loci(), 12);
    }

    /// **§12.6 above the cap, which is where it was never checked.** The shipped test
    /// for "the read-group table equals the windowed one folded, on a single-library
    /// sample" tops out at depth 20 — below the cap, where neither grain draws. Above
    /// it both draw, and they agree only because the draw is a pure function of
    /// `(seed, population, successes, draws)` and a lone read group presents the same
    /// four. That is an argument; this is the test.
    #[test]
    fn the_two_tables_agree_cell_for_cell_above_the_cap_too() {
        let edges = ladder();
        let mut accumulators = accumulators(&edges, &[group(0)]);

        for at in 0..40u32 {
            let depth = 300 + at * 17;
            accumulators.add_locus(&site(
                0,
                1 + u64::from(at) * 5_000,
                vec![
                    observation(b"A", group(0), depth - depth / 3),
                    observation(b"T", group(0), depth / 3),
                ],
            ));
        }

        let counts = accumulators.adjustments();
        assert_eq!(counts.sites_subsampled_to_cap, 40, "every site was capped");

        let read_group_table = &accumulators.read_group_histograms()[&(group(0), diploid())];
        let folded = accumulators.whole_sample_histogram(diploid());
        assert_eq!(read_group_table.cells(diploid()), folded.cells(diploid()));
        assert_eq!(
            read_group_table.total_covered_positions(),
            folded.total_covered_positions()
        );
    }

    /// **Two loci that overlap are counted, not repaired.** Overlap would enter a site
    /// into the windowed table twice, breaking "a site enters once, whole" — and a
    /// de-duplication rule would hide the upstream bug behind a plausible number.
    #[test]
    fn two_overlapping_loci_are_counted() {
        let edges = ladder();
        let mut accumulators = accumulators(&edges, &[group(0)]);

        let mut first = site(0, 100, vec![observation(b"A", group(0), 10)]);
        first.region.end = Position(120);
        let second = site(0, 110, vec![observation(b"A", group(0), 10)]);
        let third = site(0, 200, vec![observation(b"A", group(0), 10)]);

        for locus in [&first, &second, &third] {
            accumulators.add_locus(locus);
        }

        let counts = accumulators.adjustments();
        assert_eq!(counts.loci_overlapping_previous, 1);
        assert_eq!(counts.shard_spans_overlapping, 0, "one shard, no seam");
    }

    /// **Three shards merged in every order give the same counters and the same cells**,
    /// including the counter that must read zero. That is what makes a region-sharded
    /// walk exact: the shards are added as integers, so the order they finish in cannot
    /// move a fitted rate.
    ///
    /// The shards are deliberately handed over **out of genome order** in some
    /// permutations, which is the case that broke the collapsed first-start/last-end
    /// check: merging `[0,100)`, `[200,300)` and then `[100,200)` gives `[0,300)` after
    /// the first two, and the third then reads as an overlap that does not exist.
    #[test]
    fn three_shards_merged_in_every_order_give_the_same_counters() {
        let edges = ladder();
        let walk = a_walk();
        let shards: Vec<Vec<&SampleLocusObservations>> =
            walk.chunks(4).map(|chunk| chunk.iter().collect()).collect();

        let mut whole = accumulators(&edges, &[group(0)]);
        for locus in &walk {
            whole.add_locus(locus);
        }

        for order in [
            [0, 1, 2],
            [2, 1, 0],
            [1, 2, 0],
            [0, 2, 1],
            [1, 0, 2],
            [2, 0, 1],
        ] {
            let mut merged = accumulators(&edges, &[group(0)]);
            for at in order {
                let mut shard = accumulators(&edges, &[group(0)]);
                for locus in &shards[at] {
                    shard.add_locus(locus);
                }
                merged.merge(shard);
            }

            assert_eq!(merged.adjustments(), whole.adjustments(), "order {order:?}");
            assert_eq!(
                merged.adjustments().loci_overlapping_previous,
                0,
                "the counter that must read zero, order {order:?}"
            );
            assert_eq!(merged.adjustments().shard_spans_overlapping, 0);
            assert_eq!(
                merged.whole_sample_histogram(diploid()).cells(diploid()),
                whole.whole_sample_histogram(diploid()).cells(diploid()),
                "order {order:?}"
            );
            // **The scale's denominator too, and it is the one field here that could have failed
            // this.** Both of its numbers are integers, so the equality is exact rather than
            // approximate; a running `f64` sum would have agreed with the unsharded walk in some
            // orders and not others, and no other assertion in this test would have noticed.
            assert_eq!(
                merged.minted_errors(),
                whole.minted_errors(),
                "the minted-error totals, order {order:?}"
            );
        }
    }

    /// **Every counter and both tables survive a merge**, which nothing checked: the
    /// sharding test above compares five zeros against five zeros, because its walk is
    /// twelve clean loci, and the two tests that do exercise counters never merge. All
    /// four counter sums could be deleted, and the whole read-group table dropped, with
    /// the suite green — the second of those means a region-sharded walk producing an
    /// **empty** read-group table, which is every error rate this step exists to fit.
    #[test]
    fn merge_carries_every_counter_and_both_tables() {
        let edges = ladder();

        let shard_of = |contig: u32, start: u64| {
            let mut shard = accumulators(&edges, &[group(0)]);
            let mut with_upstream_cap = site(contig, start, vec![observation(b"A", group(0), 20)]);
            with_upstream_cap.reads_discarded_by_cap = 7;
            with_upstream_cap.reads_without_observation = 3;
            shard.add_locus(&with_upstream_cap);
            // Deeper than the cap, so this step subsamples it.
            shard.add_locus(&site(
                contig,
                start + 10,
                vec![
                    observation(b"A", group(0), 400),
                    observation(b"T", group(0), 100),
                ],
            ));
            // Overlaps the one before it.
            let mut overlapping = site(contig, start + 10, vec![observation(b"A", group(0), 9)]);
            overlapping.region.end = Position(start + 12);
            shard.add_locus(&overlapping);
            shard
        };

        let mut merged = shard_of(0, 1_000);
        merged.merge(shard_of(1, 1_000));

        let counts = merged.adjustments();
        assert_eq!(counts.loci_with_upstream_subsample, 2, "one per shard");
        assert_eq!(counts.reads_without_observation, 6);
        assert_eq!(counts.sites_subsampled_to_cap, 2);
        assert_eq!(counts.loci_overlapping_previous, 2);
        assert_eq!(counts.shard_spans_overlapping, 0, "different contigs");

        let by_group = merged.read_group_histograms();
        assert_eq!(by_group.len(), 1, "one library, one ploidy");
        assert_eq!(
            by_group[&(group(0), diploid())].total_loci(),
            6,
            "three loci from each shard reached the read-group table"
        );
        assert_eq!(merged.whole_sample_histogram(diploid()).total_loci(), 6);

        // **The calibration totals, by value and not by agreement.** The assertion in the
        // six-order test above compares one merge against another, so it catches an
        // order-dependent fold and nothing else — a build that dropped the sum entirely would
        // satisfy it, because both sides would be zero.
        //
        // Each shard contributes 20 + 400 + 100 + 9 = 529 reads, and `observation` gives each
        // one `−0.1 × num_obs − 0.07`, so the four sum to −53.18 a shard and −106.36 over both.
        let minted = merged.minted_errors();
        assert_eq!(minted.len(), 1, "one library");
        let group = minted[&group(0)];
        assert_eq!(
            group.reads(),
            1_058,
            "every read of both shards, and note what that means about the depth cap: the \
             500-read locus entered the histogram thinned to 124, and all 500 of its reads \
             count here. The cap draws on counts and never looks at a quality, so it cannot \
             bias a mean — only a sum, which this is not",
        );
        // **The tolerance is the bound `MintedReadErrors` documents**, half a unit of 2⁻²⁰, so
        // this assertion is what checks that claim rather than a number picked to pass. The miss
        // measured here is 1.226 × 10⁻⁹ — eight fixed-point conversions' worth of rounding (four
        // observations in each of two shards) spread over 1,058 reads — which is 389 times
        // inside it, and any wrong sum would be orders of magnitude outside.
        let mean = group.mean_log_error().expect("reads were seen");
        assert!(
            (mean - (-106.36 / 1_058.0)).abs() < 0.5 / (1_u64 << 20) as f64,
            "the mean log error came back {mean}",
        );
    }

    /// **Both boundaries of the overlap comparison.** A locus starting exactly where the
    /// previous one ended overlaps by one position and is counted; one starting the
    /// position after does not. And a locus wholly *contained* in its predecessor is
    /// counted too, and does not shorten the reach the next locus is checked against —
    /// which a naive `previous = end` would have done.
    #[test]
    fn touching_and_contained_loci_are_both_counted() {
        let edges = ladder();

        let stretch = |start: u64, end: u64| {
            let mut locus = site(0, start, vec![observation(b"A", group(0), 8)]);
            locus.region.end = Position(end);
            locus
        };

        let mut touching = accumulators(&edges, &[group(0)]);
        touching.add_locus(&stretch(100, 200));
        touching.add_locus(&stretch(200, 300));
        assert_eq!(
            touching.adjustments().loci_overlapping_previous,
            1,
            "position 200 is in both"
        );

        let mut abutting = accumulators(&edges, &[group(0)]);
        abutting.add_locus(&stretch(100, 200));
        abutting.add_locus(&stretch(201, 300));
        assert_eq!(
            abutting.adjustments().loci_overlapping_previous,
            0,
            "abutting is not overlapping"
        );

        let mut contained = accumulators(&edges, &[group(0)]);
        contained.add_locus(&stretch(100, 300));
        contained.add_locus(&stretch(150, 200));
        contained.add_locus(&stretch(250, 400));
        assert_eq!(
            contained.adjustments().loci_overlapping_previous,
            2,
            "the contained locus, and the one still inside the first's reach"
        );
    }

    /// `merge`'s four guards each refuse a shard that would make the merged cells mean
    /// something other than what they say. The ploidy map is the one that was unguarded
    /// until this review: two shards handed different maps merged silently, and every
    /// cell one of them filled was scored against the wrong set of genotype classes.
    #[test]
    #[should_panic(expected = "different ploidy maps")]
    fn merging_shards_handed_different_ploidy_maps_is_refused() {
        let edges = ladder();
        let mut mine = accumulators(&edges, &[group(0)]);
        mine.merge(GenericAccumulators::new(
            Arc::clone(&edges),
            &[group(0)],
            Arc::new(ConstantPloidy(diploid())),
            InbreedingMode::Fitted,
        ));
    }

    #[test]
    #[should_panic(expected = "different ladder objects")]
    fn merging_shards_binned_by_different_ladders_is_refused() {
        let edges = ladder();
        let mut mine = accumulators(&edges, &[group(0)]);
        mine.merge(accumulators(&ladder(), &[group(0)]));
    }

    #[test]
    #[should_panic(expected = "whether the sample has more than one library")]
    fn merging_shards_that_disagree_about_the_library_count_is_refused() {
        let edges = ladder();
        let mut mine = accumulators(&edges, &[group(0)]);
        mine.merge(accumulators(&edges, &[group(0), group(1)]));
    }

    #[test]
    #[should_panic(expected = "whether F is fitted")]
    fn merging_shards_that_disagree_about_inbreeding_is_refused() {
        let edges = ladder();
        let mut mine = accumulators(&edges, &[group(0)]);
        mine.merge(GenericAccumulators::new(
            Arc::clone(&edges),
            &[group(0)],
            shared_ploidy(),
            InbreedingMode::Supplied(InbreedingF::try_new(0.1).expect("a coefficient")),
        ));
    }

    /// **Shards whose spans overlap are reported**, because while that counter is
    /// non-zero the locus counter is not the count an unsharded walk would have given —
    /// a locus of one shard may overlap a locus of another and neither shard could see
    /// it.
    #[test]
    fn shards_covering_the_same_stretch_are_reported_as_a_seam() {
        let edges = ladder();
        let mut merged = accumulators(&edges, &[group(0)]);

        for start in [100u64, 300] {
            let mut shard = accumulators(&edges, &[group(0)]);
            for at in 0..3u64 {
                shard.add_locus(&site(
                    0,
                    start + at * 100,
                    vec![observation(b"A", group(0), 8)],
                ));
            }
            merged.merge(shard);
        }

        let counts = merged.adjustments();
        assert_eq!(counts.loci_overlapping_previous, 0, "neither shard saw one");
        assert_eq!(counts.shard_spans_overlapping, 1, "but their spans do");
    }

    /// The windows are 100 kb and a locus is filed by where it **starts**. ng's
    /// `Position` is 1-based, so position 1 opens window 0 and position 100,001 opens
    /// window 1 — every window holding exactly 100,000 positions, including the first.
    #[test]
    fn a_locus_is_filed_by_the_window_its_start_falls_in() {
        let at = |start: u64| match WindowKey::of_locus(GenomeRegion {
            contig: ContigId(3),
            start: Position(start),
            end: Position(start),
        }) {
            WindowKey::InWindow { window, contig } => {
                assert_eq!(contig, ContigId(3));
                window.get()
            }
            WindowKey::WholeSample => unreachable!("the mode is Fitted here"),
        };

        assert_eq!(at(1), 0, "the first position of a contig");
        assert_eq!(at(100_000), 0, "the last position of window 0");
        assert_eq!(at(100_001), 1);
        assert_eq!(at(200_000), 1);
        assert_eq!(at(200_001), 2);
    }

    /// **A supplied `F` drops the window key**, collapsing the windowed object from one
    /// table per 100 kb to one for the sample. The object itself stays, because the
    /// per-site rates still need whole sites.
    ///
    /// **And the collapsed arm is a different table with a wider counter, which is the
    /// property this asserts rather than a detail of where it lives.** The windowed map's
    /// `u32` is justified by a window holding at most 100,000 sites; a table with no window
    /// key holds the genome, and before the split that overflowed at 143 million sites at
    /// depth 30 in one cell. So the assertion is that the `u32` map receives **nothing** in
    /// this mode — a site landing there is a site on the path that overflowed.
    #[test]
    fn a_supplied_inbreeding_coefficient_collapses_the_windowed_table() {
        let edges = ladder();
        let supplied = InbreedingF::try_new(0.3).expect("a coefficient in [0, 1)");
        let mut with_supplied_f = GenericAccumulators::new(
            Arc::clone(&edges),
            &[group(0)],
            shared_ploidy(),
            InbreedingMode::Supplied(supplied),
        );
        for locus in a_walk() {
            with_supplied_f.add_locus(&locus);
        }

        assert!(
            with_supplied_f.windowed_histograms().is_empty(),
            "the u32-counted map must receive nothing when there is no window key to bound \
             what one of its tables can hold"
        );
        assert_eq!(
            with_supplied_f.ploidies(),
            BTreeSet::from([diploid()]),
            "and the ploidies are still findable, which the coupled fit needs"
        );
        assert_eq!(
            with_supplied_f
                .whole_sample_histogram(diploid())
                .total_loci(),
            12,
            "one table, and it still holds every site"
        );

        // Against the same walk with F fitted, which is where the 37 MB goes.
        let mut fitted = accumulators(&edges, &[group(0)]);
        for locus in a_walk() {
            fitted.add_locus(&locus);
        }
        assert!(
            fitted.windowed_histograms().len() > 1,
            "fitted keeps one table per window"
        );
        assert_eq!(
            fitted.whole_sample_histogram(diploid()).cells(diploid()),
            with_supplied_f
                .whole_sample_histogram(diploid())
                .cells(diploid()),
            "and the fold over them is the same table either way"
        );

        // **The scale's denominator does not know about `F`, and this is what says so.** The
        // calibration fold sits before `add_locus`'s mode branch, so both modes must produce the
        // same totals — and until this assertion existed, moving the fold inside the `Fitted`
        // arm left all 274 tests in this module's parent green while a supplied-`F` run silently
        // shipped an empty denominator. Every other accumulator test builds through a helper
        // that hard-codes `Fitted`.
        assert_eq!(
            fitted.minted_errors(),
            with_supplied_f.minted_errors(),
            "the minted-error totals are the same whether F was fitted or supplied"
        );
        assert!(
            !with_supplied_f.minted_errors().is_empty(),
            "and they are not empty, or the assertion above would hold by both sides being so"
        );
    }

    /// A multi-library sample keys its windowed sites by which library each alternative
    /// read came from, and the read-group table still enters each group at its own
    /// depth.
    #[test]
    fn a_multi_library_sample_keeps_the_attribution_in_the_windowed_table() {
        let edges = ladder();
        let mut accumulators = accumulators(&edges, &[group(0), group(1)]);

        accumulators.add_locus(&site(
            0,
            1_000,
            vec![
                observation(b"A", group(0), 12),
                observation(b"T", group(0), 1),
                observation(b"A", group(1), 7),
                observation(b"T", group(1), 1),
            ],
        ));

        assert_eq!(accumulators.read_group_histograms().len(), 2);
        let windowed = accumulators.whole_sample_histogram(diploid());
        let cells = windowed.cells(diploid());
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0].key.alt_reads(), 2);
        assert!(
            cells[0].key.attribution() != super::super::histogram::Attribution::Pooled,
            "two libraries showed one alternative read each"
        );
    }

    /// The tallies that describe what the pour had to live with: an upstream cap, reads
    /// that witnessed nothing, and this step's own subsampling.
    #[test]
    fn the_adjustments_report_what_the_pour_had_to_live_with() {
        let edges = ladder();
        let mut accumulators = accumulators(&edges, &[group(0)]);

        let mut upstream_capped = site(0, 100, vec![observation(b"A", group(0), 20)]);
        upstream_capped.reads_discarded_by_cap = 30;
        upstream_capped.reads_without_observation = 4;
        accumulators.add_locus(&upstream_capped);

        accumulators.add_locus(&site(
            0,
            200,
            vec![
                observation(b"A", group(0), 400),
                observation(b"T", group(0), 100),
            ],
        ));

        let counts = accumulators.adjustments();
        assert_eq!(counts.loci_with_upstream_subsample, 1);
        assert_eq!(counts.reads_without_observation, 4);
        assert_eq!(counts.sites_subsampled_to_cap, 1, "the 500-read site");
        assert_eq!(counts.loci_overlapping_previous, 0);
    }
}
