//! Filling the census during the walk.
//!
//! **The writer is the one part of the census that knows how a walk produces loci** — the
//! locus stream, the reference, the base-quality and misread accumulators. The data types it
//! fills live in the parent module, which the fit, the file codec and contamination read
//! without depending on any of that.

use std::collections::BTreeMap;

use crate::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use crate::parameter_estimation::calibration::MintedReadErrors;
use crate::parameter_estimation::depth_bins::DepthBinEdges;
use crate::parameter_estimation::joint::loci::{CensusLoci, CensusLociDigester, SelectionTerms};
use crate::parameter_estimation::misread_after_runs::MisreadCounter;
use crate::ref_seq::RefSeq;
use crate::repeat_catalog::StratumCounts;
use crate::types::{ContigId, GenomePosition, GenomeRegion, Position, ReadGroupId};

use super::{
    AlleleObservation, DepthCap, DepthCode, DepthLadderDigest, GenericEvidence, GuardObservation,
    NamedReadGroup, ObservedAllele, ReadCap, RecordingTerms, SampleCensusEvidence, Section,
    SectionKey, SelectionTermsDigest, SsrEvidence, Stratum, TractDifference,
};

/// Fills a sample's records as the walk hands it loci.
///
/// It borrows each locus rather than taking it, so the walk can hand the same locus on to the
/// rest of the run. It knows which loci are kept and ignores the rest.
pub struct CensusWriter {
    /// The kept generic positions, in selection order — the index a record entry carries.
    generic_loci: Vec<GenomePosition>,
    /// The kept STR loci as regions, in genome order.
    ssr_loci: Vec<GenomeRegion>,
    /// Which stratum each of those loci is in, in the same order. **The writer carries it
    /// because the tract half is stored one stratum at a time**, and the locus stream says which
    /// tract a read reached, never which stratum it belongs to.
    ssr_stratum_of: Vec<Stratum>,
    /// Where each of those loci sits **within its own stratum**, in the same order — the index
    /// every entry of an [`SsrEvidence`] carries, since a section holds one stratum's tracts and
    /// a genome-wide index would name a tract it does not have.
    ssr_index_in_stratum: Vec<u32>,
    /// How many tracts each stratum holds, which is the length a section of it is built at
    /// whether or not a read ever reaches one of them.
    ssr_stratum_size: BTreeMap<Stratum, usize>,
    edges: DepthBinEdges,
    generic: BTreeMap<ReadGroupId, GenericEvidence>,
    /// Per read group, then per stratum. **Two levels rather than one key**, because the
    /// walk asks twice for "every read group seen so far" and the outer map answers it.
    ssr: BTreeMap<ReadGroupId, BTreeMap<Stratum, SsrEvidence>>,
    /// Non-reference observations before they are sorted — the walk arrives in position
    /// order, so this is already sorted in practice and the sort at `finish` is a guard.
    pending: BTreeMap<ReadGroupId, Vec<AlleleObservation>>,
    digester: CensusLociDigester,
    /// How far the digester has been fed, so a locus arriving out of order is caught.
    digested: usize,
    terms: SelectionTerms,
    read_cap: ReadCap,
    depth_cap: DepthCap,
    stratum_counts: StratumCounts,
    /// Every read group the sample declares. **Held, rather than discovered from the
    /// observations**, so a group that put no read at a position still gets its zero there:
    /// discovering them would make a group's record start at its first read, and every
    /// position before that indistinguishable from never walked.
    read_groups: Vec<ReadGroupId>,
    /// **Who each of those groups is** — the `@RG ID` and the library, which travel into the
    /// census because a cohort of censuses can only be merged on them.
    declared: BTreeMap<ReadGroupId, NamedReadGroup>,
    /// **Σ `ln P(this read is wrong)` and the read count, per read group**, accumulated as the
    /// loci go past. What a base-quality calibration is fitted from, and the half of it a census
    /// carries — see [`SampleCensusEvidence::minted_read_errors`].
    minted: BTreeMap<ReadGroupId, MintedReadErrors>,
    /// Scratch for the accumulator, cleared and refilled once a locus rather than allocated.
    minted_scratch: Vec<(ReadGroupId, MintedReadErrors)>,
    /// **The count of misreads after runs**, where the walk handed the writer its reference
    /// ([`Self::count_misreads_after_runs`]); `None` writes an empty table.
    misread_counter: Option<MisreadCounter>,
    sample: String,
    /// Ground whose loci this census does not count, in genome order: the satellites a
    /// `--call-satellites` walk handed to the generic path ([`Self::leave_out`]). Empty otherwise.
    left_out: Vec<GenomeRegion>,
}

/// How many of `reads` survive when a position's `depth` reads are thinned down to `cap`.
///
/// **Proportional and rounded to nearest, never a draw.** A region-sharded walk and a
/// single-threaded one must produce byte-identical records, and a share needs no seed to do
/// that. Rounding to nearest keeps a single stray read at a very deep position rather than
/// discarding it, which matters because those reads are what the error rate is fitted from.
///
/// **It returns the byte the record stores**, which is safe for the reason [`DepthCap`]
/// exists: a count cannot exceed its position's depth, the surviving share cannot exceed the
/// cap, and the cap cannot exceed 255. The conversion says so rather than truncating — a
/// truncated count is a wrong allele fraction and has no symptom.
///
/// # Panics
///
/// When the surviving count will not fit in a byte, which means an allele carried more reads
/// than the position it sits at held.
pub(super) fn thin_to_cap(reads: u32, depth: u32, cap: DepthCap) -> u8 {
    let cap = cap.get();
    let surviving = if depth <= cap || depth == 0 {
        reads
    } else {
        let scaled = f64::from(reads) * f64::from(cap) / f64::from(depth);
        (scaled.round() as u32).min(cap).max(u32::from(reads > 0))
    };
    u8::try_from(surviving).unwrap_or_else(|_| {
        panic!(
            "{surviving} reads on one allele where the position held {depth} and the cap is \
             {cap}: a count cannot exceed the depth it was thinned to, and the record holds \
             one byte"
        )
    })
}

impl CensusWriter {
    #[allow(
        clippy::too_many_arguments,
        reason = "every one is a distinct input the records cannot be built without; \
                  bundling them into a config struct would only move the list"
    )]
    pub fn new(
        sample: String,
        loci: &CensusLoci,
        // **The sample's read groups and who each one is, in one argument** — two lists that
        // had to agree would be a way for a census to record a library under another's
        // identifier, with no symptom.
        declared: BTreeMap<ReadGroupId, NamedReadGroup>,
        contig_of: &dyn Fn(&str) -> Option<ContigId>,
        terms: SelectionTerms,
        edges: DepthBinEdges,
        read_cap: ReadCap,
        depth_cap: DepthCap,
    ) -> Self {
        // **The stratum travels beside the locus from here on**, because the locus stream the
        // walk hands over names a tract by where it sits and never by which stratum it is in,
        // and two of a tract's counts are kept per stratum.
        let mut ssr_loci: Vec<(GenomeRegion, Stratum)> = loci
            .ssr()
            .iter_sorted()
            .into_iter()
            .flat_map(|((period, reference_repeats), segments)| {
                segments.iter().map(move |segment| {
                    (
                        segment,
                        Stratum {
                            period,
                            reference_repeats,
                        },
                    )
                })
            })
            .filter_map(|(segment, stratum)| {
                contig_of(segment.chrom()).map(|contig| {
                    (
                        GenomeRegion {
                            contig,
                            start: Position(segment.start()),
                            end: Position(segment.end()),
                        },
                        stratum,
                    )
                })
            })
            .collect();
        ssr_loci.sort_unstable_by_key(|(r, _)| (r.contig.get(), r.start.get(), r.end.get()));
        let ssr_stratum_of: Vec<Stratum> = ssr_loci.iter().map(|(_, stratum)| *stratum).collect();
        // **Numbered within its stratum in genome order**, so that a section's tracts run in the
        // same order the whole selection does. The fit sums over a stratum's tracts, and a sum of
        // floats is not associative — an order that depended on how the strata were interleaved
        // would move the answer's last digits for no reason anybody could name.
        let mut ssr_stratum_size: BTreeMap<Stratum, usize> = BTreeMap::new();
        let mut ssr_index_in_stratum: Vec<u32> = Vec::with_capacity(ssr_stratum_of.len());
        for stratum in &ssr_stratum_of {
            let so_far = ssr_stratum_size.entry(*stratum).or_insert(0);
            ssr_index_in_stratum.push(*so_far as u32);
            *so_far += 1;
        }
        let ssr_loci = ssr_loci.into_iter().map(|(region, _)| region).collect();
        Self {
            generic_loci: loci.generic().to_vec(),
            ssr: BTreeMap::new(),
            ssr_loci,
            ssr_stratum_of,
            ssr_index_in_stratum,
            ssr_stratum_size,
            edges,
            generic: BTreeMap::new(),
            pending: BTreeMap::new(),
            digester: CensusLociDigester::new(),
            digested: 0,
            terms,
            read_cap,
            depth_cap,
            stratum_counts: loci.ssr_stratum_counts().clone(),
            read_groups: declared.keys().copied().collect(),
            declared,
            minted: BTreeMap::new(),
            minted_scratch: Vec::new(),
            misread_counter: None,
            sample,
            left_out: Vec::new(),
        }
    }

    /// Record this locus if it is a kept one.
    ///
    /// **Every kept locus gets an entry whether or not a read reached it** — the entry is the
    /// denominator, and a locus in a region never walked keeps [`DepthCode::NeverWalked`], so
    /// the three states survive.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations) {
        if self.is_left_out(locus.region) {
            return;
        }
        // **Every generic locus the walk hands over, not only the kept ones.** This is the
        // accumulator the per-sample calibration pre-pass uses, with its own unit unchanged: a
        // read at a position, counted once for every position it is seen at, at generic loci,
        // over complete witnesses, before the per-position depth cap. Restricting it to the
        // census's kept positions would be a second definition of a per-read-group total, and
        // the one number it feeds — how far a library's own base qualities may be trusted — is
        // a property of the library rather than of which positions were kept.
        crate::parameter_estimation::calibration::minted_error_by_read_group(
            locus,
            &mut self.minted_scratch,
        );
        for (group, totals) in self.minted_scratch.drain(..) {
            self.minted.entry(group).or_default().add(totals);
        }

        if let Some(counter) = &mut self.misread_counter {
            counter.add_locus(locus);
        }

        match &locus.kind {
            LocusKind::Generic => self.add_generic(locus),
            LocusKind::Ssr(_) => self.add_ssr(locus),
            LocusKind::SsrBundle(_) => {}
        }
    }

    /// **Count nothing from loci that start on `spans`**, which must be in genome order and not
    /// overlap.
    ///
    /// For the satellites a `--call-satellites` walk hands to the generic path
    /// (`command_line.md` §4.4.1): the option changes nothing in the fit, so their loci must not
    /// reach it. The positions the census keeps are drawn from the catalog's ordinary sequence
    /// and never lie there, but [`add_locus`](Self::add_locus) also tallies every generic
    /// locus's read errors for the base-quality calibration, and satellite loci in that tally
    /// moved each library's error multiplier in its seventh figure on the test cohort.
    pub fn leave_out(&mut self, spans: &[GenomeRegion]) {
        self.leave_out_spans(spans);
    }

    /// **Count misreads after runs from here on, reading the walk's reference through
    /// `reference`** (`str_context_error_fit.md` §3). Without it the census carries an empty
    /// table, which every fit reads as nothing counted.
    pub fn count_misreads_after_runs(&mut self, reference: Box<dyn RefSeq + Send>) {
        self.misread_counter = Some(MisreadCounter::new(reference));
    }

    fn leave_out_spans(&mut self, spans: &[GenomeRegion]) {
        self.left_out = spans.to_vec();
    }

    /// Whether `region` starts on ground [`leave_out`](Self::leave_out) named.
    fn is_left_out(&self, region: GenomeRegion) -> bool {
        let at = (region.contig, region.start);
        let first = self
            .left_out
            .partition_point(|span| (span.contig, span.end) < at);
        self.left_out
            .get(first)
            .is_some_and(|span| (span.contig, span.start) <= at)
    }

    /// This stretch of genome was walked, whatever the walk found in it.
    ///
    /// **Without this the three states of §1.1 collapse to two on real data.** The generic
    /// locus generator emits a locus only where a read reached, so a position no read reached
    /// produces nothing and would keep [`DepthCode::NeverWalked`] — which is the code for a
    /// bug, a region the run never opened. Measured on tomato SRR7279482 at 25×, that is
    /// 93,150 of 1,999,404 kept positions, 1 in 21, every one of them data being reported as
    /// a defect.
    ///
    /// Call it for each region handed to the walk, before or after the loci from it: a real
    /// depth overwrites a zero and a zero never overwrites a real depth.
    pub fn mark_walked(&mut self, region: GenomeRegion) {
        let zero = DepthCode::Binned(self.edges.bin_for(0));
        let first = self.generic_loci.partition_point(|kept| {
            (kept.contig.get(), kept.position.get()) < (region.contig.get(), region.start.get())
        });
        for index in first..self.generic_loci.len() {
            let kept = self.generic_loci[index];
            if kept.contig != region.contig || kept.position.get() > region.end.get() {
                break;
            }
            for group in &self.read_groups {
                let records = self
                    .generic
                    .entry(*group)
                    .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
                if records.depth.get(index) == DepthCode::NeverWalked {
                    records.depth.set(index, zero);
                }
            }
        }
    }

    fn add_generic(&mut self, locus: &SampleLocusObservations) {
        // **A locus the selection kept nothing in writes nothing, so it is left before the
        // depth array is built.** Every write below sits inside `for offset in 0..span`
        // behind `generic_index`, which answers `None` at every position of such a locus —
        // so the depth array is built, filled from every observation, and thrown away. On
        // this 2 Mb fixture that is 4 loci in 5; on a whole-genome walk it is 399 in 400,
        // because the census budget is the same two million positions over four hundred
        // times the ground.
        if !self.keeps_a_position_in(locus.region) {
            return;
        }
        // **Depth is per read group, not pooled.** `num_obs_along_locus` sums the sample's
        // observations, which is the wrong grain for a record keyed by read group — and it
        // is the grain the error rate is fitted at, so pooling here would score every read
        // against a rate fitted from a depth its own library never had.
        let span = locus.region.len() as usize;
        let mut depths: BTreeMap<ReadGroupId, Vec<u32>> = BTreeMap::new();
        for group in &self.read_groups {
            depths.insert(*group, vec![0; span]);
        }
        for observation in &locus.observations {
            let Some(depth) = depths.get_mut(&observation.read_group) else {
                continue;
            };
            match &observation.read_witness {
                ReadWitness::Complete => {
                    for slot in depth.iter_mut() {
                        *slot = slot.saturating_add(observation.num_obs);
                    }
                }
                ReadWitness::Partial { positions } => {
                    for (start, end) in positions.runs() {
                        let from = (start as usize).min(span);
                        let to = (end as usize).min(span).max(from);
                        for slot in &mut depth[from..to] {
                            *slot = slot.saturating_add(observation.num_obs);
                        }
                    }
                }
            }
        }

        // The whole locus's observations describe the whole locus's span; only a
        // single-base locus lets "which base did this read show here" be answered from them.
        // A wider one contributes its depth, which is what the fit reads it for.
        let single_base = span == 1;
        for offset in 0..span {
            let position = GenomePosition {
                contig: locus.region.contig,
                position: Position(locus.region.start.get() + offset as u64),
            };
            let Some(index) = self.generic_index(position) else {
                continue;
            };
            self.feed_digest(index, position);

            // **Every read group gets an entry at every walked position**, whether or not it
            // put a read there: the entry is the denominator, and a group with no read here
            // saw zero, which is data and not the absence of it.
            //
            // **The depth recorded is the one the position actually held** (§2.2, changed
            // 2026-08-14). The cap below still thins the allele counts, which is what keeps
            // them inside their own field and what bounds the sparse list at high depth; it
            // no longer clips this. Nothing is lost by that and one thing is gained: a
            // consumer asking *how many copies of this stretch does the sample carry* reads
            // a depth that can answer, where a clipped one went blind at twice the cap.
            // A consumer needing the counts' own denominator computes it with
            // [`DepthCap::denominator_for`], the cap travelling in the recording terms.
            for (group, depth) in &depths {
                let records = self
                    .generic
                    .entry(*group)
                    .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
                records
                    .depth
                    .set(index, DepthCode::Binned(self.edges.bin_for(depth[offset])));
            }
            if !single_base {
                continue;
            }
            for observation in locus.complete_observations() {
                // The one definition of "non-reference", shared with the cohort merge
                // rather than spelled here — two spellings of one test are two things
                // that can disagree. Complete observations only, which is what lets the
                // whole locus's reference bases be the right stretch to compare against.
                if observation.matches_reference(&locus.reference_bases) {
                    continue;
                }
                let allele = match &*observation.bases {
                    [base] => ObservedAllele::of_base(*base),
                    // An insertion or a deletion at a one-base locus: not one of the four,
                    // and the fit scores it as the fifth rather than guessing a base.
                    _ => ObservedAllele::Other,
                };
                // **Thinned by the same ratio as the depth**, so the fractions the reads showed
                // survive the cap exactly. Deterministic rather than a draw: a region-sharded
                // walk and a single-threaded one must produce byte-identical records, and a
                // proportional share needs no seed to do that.
                let group_depth = depths
                    .get(&observation.read_group)
                    .map_or(0, |depth| depth[offset]);
                let reads = thin_to_cap(observation.num_obs, group_depth, self.depth_cap);
                if reads == 0 {
                    continue;
                }
                self.pending
                    .entry(observation.read_group)
                    .or_default()
                    .push(AlleleObservation {
                        index: index as u32,
                        allele,
                        reads,
                    });
            }
        }
    }

    fn add_ssr(&mut self, locus: &SampleLocusObservations) {
        let Some(index) = self.ssr_index(locus.region) else {
            return;
        };
        let LocusKind::Ssr(detail) = &locus.kind else {
            return;
        };
        let period = detail.motif.period().max(1) as i64;
        let reference_length = locus.region.len() as i64;
        let stratum = self.ssr_stratum_of[index];
        // The index every entry below carries: this tract's place **in its own stratum**, since
        // that is the whole of what a section holds.
        let local = self.ssr_index_in_stratum[index] as usize;
        let tracts_here = self.ssr_stratum_size[&stratum];
        // **How many of this read group's reads at this locus have been numbered already.**
        // The read a difference sits on is numbered within the *locus*, not within the
        // observation it arrived in: the walk folds reads carrying identical bases into one
        // observation, so a counter that restarted at each observation gave two reads that
        // differ from the reference in different places the same number, and a consumer
        // grouping by it saw one bad read where the truth is an allele on two reads.
        //
        // It advances over the reads that were **compared** — the unslipped ones that reach
        // the difference loop below — because those are the only reads that can carry a
        // difference, and every number spent elsewhere is headroom lost from a byte.
        let mut numbered: BTreeMap<ReadGroupId, u32> = BTreeMap::new();

        for observation in &locus.observations {
            let records = self
                .ssr
                .entry(observation.read_group)
                .or_default()
                .entry(stratum)
                .or_insert_with(|| SsrEvidence::never_walked(tracts_here));
            records.walked.set(local);
            let reads = u16::try_from(observation.num_obs).unwrap_or(u16::MAX);
            if observation.read_witness != crate::locus_generation::ReadWitness::Complete {
                // A read that covered the tract without crossing it reports no length; it is a
                // lower bound. **Counted into this stratum's total rather than at the locus**
                // (spec §3): the fit summed the per-locus version the moment it read it, and
                // the loss runs along repeat count, which is what a stratum is.
                records.covering_not_crossing += u64::from(reads);
                continue;
            }
            // **A read missing flank bases is a different sequence, not a slip** (owner,
            // 2026-10-10; `str_edge_deletion_alleles.md` §5.3): its spelling is the repeat a
            // deletion into the flank left behind, and counting its length here would charge
            // that deletion to the polymerase.
            if !observation.flank_deletion.is_none() {
                continue;
            }
            let difference = observation.bases.len() as i64 - reference_length;
            if difference % period == 0 {
                records.offsets[local].add((difference / period) as i32, reads);
            } else {
                records.guard.push(GuardObservation {
                    locus: local as u32,
                    length_difference: difference as i32,
                    reads,
                });
            }
            // **The difference list, and the one case it can be built from.** A read whose
            // tract is the reference's length lines up with it base for base, so a mismatch is
            // read off directly. A read that slipped does not: which of its bases sits over
            // which reference base is the aligner's answer and not this writer's, and inventing
            // one would put a whole tract's worth of manufactured mismatches into the error
            // rate. Such a read contributes its length to the offsets and **nothing to the
            // denominator**, so the rate stays a ratio of two quantities counted over the same
            // reads.
            if observation.bases.len() as i64 != reference_length {
                continue;
            }
            // **Per stratum, for the same reason**: it is the denominator of a rate fitted per
            // stratum, and per locus it was 1.85 MB a read group on tomato that nothing read.
            records.bases_compared += u64::from(reads) * observation.bases.len() as u64;
            // **One entry per read, not one per distinct sequence.** The walk has already
            // folded reads carrying the same bases into a single observation with a count, and
            // collapsing them here too would lose the fact this channel exists for: that two
            // reads carried the *same* interruption, which is what makes an interruption an
            // allele rather than an error.
            let numbered_here = numbered.entry(observation.read_group).or_insert(0);
            let first_read = *numbered_here;
            *numbered_here += u32::from(reads);
            for (offset, (read_base, reference_base)) in observation
                .bases
                .iter()
                .zip(locus.reference_bases.iter())
                .enumerate()
            {
                if read_base.eq_ignore_ascii_case(reference_base) {
                    continue;
                }
                for copy in 0..u32::from(reads) {
                    records.differences.push(TractDifference {
                        locus: local as u32,
                        read: u16::try_from(first_read + copy).unwrap_or(u16::MAX),
                        offset: i16::try_from(offset).unwrap_or(i16::MAX),
                        base: ObservedAllele::of_base(*read_base),
                    });
                }
            }
        }
        // Reads that reached the locus and produced no observation at all. **Every read group
        // seen so far learns the walk reached this tract**, and the count itself is charged
        // once — to the lowest-numbered group, since the reads cannot be attributed to one and
        // the fit adds every group's total into the stratum's anyway.
        if locus.reads_without_observation > 0 {
            for strata in self.ssr.values_mut() {
                strata
                    .entry(stratum)
                    .or_insert_with(|| SsrEvidence::never_walked(tracts_here))
                    .walked
                    .set(local);
            }
            if let Some(strata) = self.ssr.values_mut().next() {
                strata
                    .entry(stratum)
                    .or_insert_with(|| SsrEvidence::never_walked(tracts_here))
                    .covering_not_crossing += u64::from(locus.reads_without_observation);
            }
        }
    }

    fn feed_digest(&mut self, index: usize, position: GenomePosition) {
        // The walk arrives in genome order, so the digest is fed in index order; a locus
        // reached twice is fed once.
        while self.digested <= index {
            let at = self.digested;
            self.digester.observe(at, self.generic_loci[at]);
            self.digested += 1;
        }
        debug_assert_eq!(self.generic_loci[index], position);
    }

    /// Whether the selection kept any generic position inside `region`.
    ///
    /// **Exact rather than conservative**, which is what lets [`add_generic`](Self::add_generic)
    /// leave on it: the kept positions are sorted by contig and position, so the first one at
    /// or after the region's start settles the question in one binary search. Same shape as
    /// [`mark_walked`](Self::mark_walked)'s scan, which walks forward from the same point.
    fn keeps_a_position_in(&self, region: GenomeRegion) -> bool {
        let start = (region.contig.get(), region.start.get());
        let at = self
            .generic_loci
            .partition_point(|kept| (kept.contig.get(), kept.position.get()) < start);
        self.generic_loci.get(at).is_some_and(|kept| {
            kept.contig == region.contig && kept.position.get() <= region.end.get()
        })
    }

    fn generic_index(&self, position: GenomePosition) -> Option<usize> {
        self.generic_loci
            .binary_search_by_key(&(position.contig.get(), position.position.get()), |kept| {
                (kept.contig.get(), kept.position.get())
            })
            .ok()
    }

    fn ssr_index(&self, region: GenomeRegion) -> Option<usize> {
        self.ssr_loci
            .binary_search_by_key(&(region.contig.get(), region.start.get()), |kept| {
                (kept.contig.get(), kept.start.get())
            })
            .ok()
    }

    /// The finished evidence, with the digest of what was actually written.
    pub fn finish(mut self) -> SampleCensusEvidence {
        // Every kept locus is digested, including any the walk never reached — the digest
        // witnesses the selection, and a run that stopped early must not produce a short
        // digest that happens to match another short one.
        while self.digested < self.generic_loci.len() {
            let at = self.digested;
            self.digester.observe(at, self.generic_loci[at]);
            self.digested += 1;
        }
        // **Every declared read group ends with a record set, even one that recorded
        // nothing.** A missing map entry and a set of never-walked entries are the same
        // situation, and only the second says so — a fit handed the first would have to
        // guess whether the group exists.
        for group in &self.read_groups {
            self.generic
                .entry(*group)
                .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
        }
        for (group, mut entries) in std::mem::take(&mut self.pending) {
            entries.sort_unstable_by_key(|entry| (entry.index, entry.allele));
            let records = self
                .generic
                .entry(group)
                .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
            records.non_reference = entries;
        }

        // **The same sections in every sample, and that is what a census is.** Every declared
        // read group gets an ordinary-position section and one section per stratum the selection
        // holds, whether or not a read ever reached one of them — a stratum with no tracts read
        // is a section holding zero, and an absent section would say instead that this sample was
        // never asked. The fit enumerates sections; two samples enumerating different ones would
        // be two samples answering different questions.
        let mut sections: BTreeMap<SectionKey, Section> = BTreeMap::new();
        for (group, records) in std::mem::take(&mut self.generic) {
            sections.insert(SectionKey::Generic(group), Section::Generic(records));
        }
        let mut tracts = std::mem::take(&mut self.ssr);
        for group in &self.read_groups {
            let mut theirs = tracts.remove(group).unwrap_or_default();
            for (stratum, size) in &self.ssr_stratum_size {
                let records = theirs
                    .remove(stratum)
                    .unwrap_or_else(|| SsrEvidence::never_walked(*size));
                sections.insert(SectionKey::Ssr(*group, *stratum), Section::Ssr(records));
            }
        }
        // **The declared groups are all the groups there are.** The generic half already drops
        // an observation whose read group the sample did not declare, so tract evidence from one
        // would be half a read group — recorded at the tracts and absent at the positions its
        // own error rate is fitted from.
        assert!(
            tracts.is_empty(),
            "read group {} put reads on a repeat tract without being declared by the sample",
            tracts.keys().next().expect("a group").get()
        );

        SampleCensusEvidence::resident(
            self.sample,
            RecordingTerms {
                selection: SelectionTermsDigest::of(&self.terms),
                kept_loci: self.digester.finish(),
                ssr_stratum_counts: self.stratum_counts,
                read_cap: self.read_cap,
                depth_ladder: DepthLadderDigest::of(&self.edges),
                depth_cap: self.depth_cap,
            },
            self.declared,
            self.minted,
            sections,
        )
        .with_misread_after_runs(
            self.misread_counter
                .map(MisreadCounter::finish)
                .unwrap_or_default(),
        )
    }
}
