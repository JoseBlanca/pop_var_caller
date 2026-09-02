//! **The input edge: one merged cohort locus, shaped into what the calling loop reads.**
//!
//! The merge hands over a [`CohortObservation`] — the locus's unified allele table and, for
//! each **covering** sample, its reads folded onto those alleles. Candidate selection then
//! narrows the allele table and says, per covering sample, what it dropped. What the loop
//! reads is neither of those: it is [`LocusEvidence`], one entry per sample **of the run**,
//! over the *candidate* allele ids.
//!
//! This module is the join, and it is data-shaping only — no arithmetic, and nothing here
//! decides anything (`doc/devel/ng/arch/calling_em_loop.md` §2, §7).
//!
//! # Three joins, and none of them is enforced by the types
//!
//! **The two per-sample lists are in different orders.** `CohortObservation::per_sample`
//! holds only the samples that covered the locus, each naming its own index in the run's
//! sample order; `LocusEvidence`'s list is one entry per run sample, so a sample that covered
//! nothing gets [`GenericSampleEvidence::empty`] — an empty sum is zero, so every genotype
//! scores alike and the prior decides alone, which is the right answer rather than a special
//! case (`doc/devel/ng/spec/read_likelihoods.md` §3.3). `LocusSelection::unmatched` is parallel
//! to the **merge's** covering samples, not to the run's
//! (`doc/devel/ng/arch/candidate_alleles.md` §5.1).
//!
//! **The row builder needs ascending `(allele, read group)`, and today's selection already
//! gives it.** `select_generic` ranks the alternatives only to decide *which* survive a binding
//! cap, and then puts the survivor list back into the merge table's own index order before
//! admitting — its own
//! `the_survivors_of_a_binding_cap_are_admitted_in_the_merge_tables_order` pins that, and was
//! written after a reviewer deleted the restoring sort and every test still passed. So the
//! remapping is order-preserving, the merge's rows arrive ascending on
//! `(merge allele, read group)`, and the narrowing hands them back ascending on the candidate
//! key with nothing to do.
//!
//! **This module sorts anyway, and the sort is insurance rather than description.**
//! [`AlleleRemap::admit`](super::allele_candidates::AlleleRemap::admit) constrains only the
//! *candidate* ids — dense and ascending — and says nothing about which table index each is
//! admitted for, so a selection that admitted by rank would be type-legal and would hand this
//! module a permutation. The read likelihood's sum must run in a fixed order
//! (`doc/devel/ng/spec/read_likelihoods.md` §8), and [`GenericSampleEvidence::new`] checks that
//! in **debug builds only** — a release run would not say so, which is the second reason to
//! sort here rather than to leave it to an assertion.
//!
//! **The pooled leftover is selection's number, and adding the narrowing's own would double
//! it.** Both are Σ `ln P(error)` over exactly the rows whose allele has no candidate —
//! selection sums them at selection time
//! ([`UnmatchedSupport::q_sum`](super::allele_candidates::UnmatchedSupport::q_sum), *"the sum
//! of the dropped rows' own `q_sum`, to the bit"*) and
//! [`GenericObservation::fill_from_supported_alleles`] sums the same rows again as it walks
//! them. This module takes selection's and **checks the narrowing's against it**, which is a
//! free cross-check between two independently written walks over one set of rows.
//!
//! # What it allocates
//!
//! **Per worker, not per locus, for everything it fills** — the narrowed rows, the two
//! per-sample maps and the remapping as a slice all live in [`GenericEvidenceScratch`] and are
//! cleared and refilled at each locus.
//!
//! **The one per-locus allocation is the caller's list of views**, and two things make it
//! unavoidable rather than one. A [`GenericLocusSample`] borrows the rows the scratch holds, so
//! a buffer holding both would have to name its own lifetime — which is why the list cannot
//! live *inside* [`GenericEvidenceScratch`]. And a `Vec` is **invariant in its element type**,
//! so a caller cannot hold one across two loci either: the element type names the lifetime of
//! the borrow, and reusing the list holds the first locus's borrow open into the second, which
//! is an `E0499` on the one-call spelling and an `E0502` on the two-step one. The list is a few
//! dozen bytes per sample, against a called locus's own output — one `Genotype` per sample —
//! which is the same order.

use std::num::NonZeroU32;

use super::SsrSampleEvidence;
use super::allele_candidates::LocusSelection;
use super::{GenericLocusSample, GenericObservation, GenericSampleEvidence, LocusEvidence};
use crate::ng::locus_generation::{SequenceObservation, SsrDetail};
use crate::ng::run::cohort_merge::build::CohortObservation;
use crate::ng::types::{AlleleId, GenomeRegion};

/// How far selection's pooled leftover and this narrowing's own sum may differ before the two
/// are treated as disagreeing, **relative to the larger of them**.
///
/// The two walk the same rows and add them up differently — selection sums each allele's rows
/// and then adds the per-allele totals, this one adds row by row — and floating-point addition
/// is not associative, so they agree only to rounding.
///
/// **Relative and not absolute, because `q_sum` is a sum of logarithms with no bound.** The two
/// walks were measured against each other over their exact association orders: 200 dropped rows
/// pooling to −3,913 differ by 1.4e-12, 600 rows at −12,314 by 1.1e-11, and 2,000 rows at
/// −4.0e6 by 3.7e-9 — which an absolute `1e-9` would report as a selection defect that did not
/// happen. **For a sample with one read group the two are bit-identical**, because selection's
/// per-allele group is then a single row, and that is most samples of most runs.
const LEFTOVER_AGREEMENT_TOLERANCE: f64 = 1e-12;

/// What one sample's dropped reads left behind: the pooled error mass, and whether the cap
/// cut a sequence this sample's own reads had earned.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct SampleLeftover {
    /// Σ `ln P(error)` over this sample's reads that showed no candidate allele.
    unmatched_q_sum: f64,
    /// Whether the locus is no longer called over something this sample carries
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0).
    genotype_must_be_missing: bool,
}

/// **Every buffer shaping one locus's evidence fills** — allocated once per worker and reused
/// at each locus, so the shaping itself costs no allocation.
///
/// **Amortised, not allocation-free**: a buffer grows when a locus is wider, or a cohort
/// larger, than any this worker has met, and capacity never comes back down — the same shape
/// and the same reason as [`CallingScratch`](super::CallingScratch)'s.
#[derive(Default, Debug)]
pub struct GenericEvidenceScratch {
    /// One entry per allele of the **merge's** table: the candidate id it now has, or nothing
    /// where selection dropped it. A flattened [`AlleleRemap`], because the narrowing takes a
    /// slice and the remapping answers one index at a time.
    ///
    /// [`AlleleRemap`]: super::allele_candidates::AlleleRemap
    candidate_of_merge_allele: Vec<Option<AlleleId>>,
    /// One entry per **run** sample, in the run's own sample order.
    each_run_sample: Vec<NarrowedRunSample>,
    /// Which locus these buffers were last narrowed for, or nothing before the first
    /// [`narrow`](GenericEvidenceScratch::narrow) — what
    /// [`fill_views`](GenericEvidenceScratch::fill_views) checks the observation it is handed
    /// against.
    ///
    /// **The region and not the covering count**, which is what this was first checked on: two
    /// *different* loci with the same number of covering samples are the common case at a
    /// cohort where most samples cover most loci, not the corner, and filling one locus's views
    /// against another's observation pairs narrowed reads with the wrong partials — both legal
    /// evidence, and neither the caller's.
    region_narrowed_for: Option<GenomeRegion>,
}

/// **One run sample's share of the narrowing** — its rows, where the merge kept it, and what
/// its dropped reads left behind.
///
/// **The three travel in one entry rather than in three lists indexed alike**, for the reason
/// [`GenericLocusSample`] holds a sample's evidence and its ruling together: a positional join
/// between per-sample lists is what this module exists to remove, and three of them here would
/// be three more of exactly that.
#[derive(Debug, Default)]
struct NarrowedRunSample {
    /// This sample's reads, narrowed onto the candidates — the storage its view borrows.
    ///
    /// **A buffer per sample rather than one flat buffer with spans**, because the narrowing
    /// this calls clears what it fills; spans would need a second, appending spelling of it,
    /// and two functions that must stay identical is the worse trade.
    ///
    /// **What it costs at a large cohort, stated rather than implied**: no buffer is ever
    /// shrunk, so the rows held at high-water are Σ over samples of *that sample's widest
    /// locus*, where a flat buffer would hold the largest single locus's total. The first is
    /// never smaller, and the two coincide only if every sample's widest locus is the same one.
    /// At a few dozen accessions the gap is bounded; at a thousand samples it is the number to
    /// come back for.
    rows: Vec<GenericObservation>,
    /// Which entry of the merge's covering list this sample is, or nothing where it covered
    /// nothing.
    merge_entry: Option<usize>,
    /// A sample that covered nothing has the default: no dropped mass and nothing set aside,
    /// which is what having shown no read means.
    leftover: SampleLeftover,
}

impl GenericEvidenceScratch {
    /// **Narrow one merged locus onto selection's candidates**, filling this worker's buffers.
    ///
    /// `run_sample_count` is the run's own sample count, which the merge does not carry: its
    /// list holds covering samples only.
    ///
    /// # Panics
    ///
    /// Held in release, because each is a caller bug whose symptom is a wrong genotype rather
    /// than a crash (`doc/devel/ng/spec/calling_em_loop.md` §8):
    ///
    /// - **the selection must be of this locus's own merged table.** Its remapping is indexed
    ///   by the merge's allele indices, so one built against a different locus maps a row onto
    ///   whatever allele happens to sit at that index.
    /// - **selection's per-sample leftovers must be one per covering sample.** That list is
    ///   parallel to the merge's covering samples and to nothing else
    ///   (`arch/candidate_alleles.md` §5.1); a length mismatch means the two were assembled
    ///   from different loci, and a positional join between them would pair one sample's
    ///   dropped reads with another's.
    /// - **every covering sample must name a sample of the run**, in strictly ascending order.
    ///   The merge builds them that way; a positional join here would otherwise write one
    ///   sample's evidence into another's slot.
    /// - **the two independently summed leftovers must agree.** Selection's pool and the
    ///   narrowing's own dropped mass are Σ over the same rows, so a disagreement means one of
    ///   the two walks changed its mind about which rows have no candidate — and the loop
    ///   would then score every sample's reads against a leftover that is too large or too
    ///   small by exactly the rows they disagree about.
    pub fn narrow(
        &mut self,
        observation: &CohortObservation,
        selection: &LocusSelection,
        run_sample_count: usize,
    ) {
        assert!(
            run_sample_count > 0,
            "a run has at least one sample, so a locus shaped for none is a run whose sample \
             order went missing rather than a locus nobody covered — which is one empty view \
             per sample"
        );
        let remap = selection.remap();
        assert_eq!(
            remap.table_len(),
            observation.alleles.len(),
            "the selection's remapping covers {} merge alleles and this locus holds {}, so \
             the two belong to different loci",
            remap.table_len(),
            observation.alleles.len()
        );
        assert_eq!(
            selection.unmatched().len(),
            observation.per_sample.len(),
            "selection's per-sample leftovers run parallel to the merge's covering samples: \
             {} leftovers arrived and {} samples covered {}",
            selection.unmatched().len(),
            observation.per_sample.len(),
            observation.region
        );

        // **Destructured rather than named through `self`**, so that a buffer added later fails
        // to compile here instead of silently carrying the previous locus's contents into the
        // next — the same device `SelectionScratch::reset_for` and the merge's own assembly use,
        // and for the same reason.
        let Self {
            candidate_of_merge_allele,
            each_run_sample,
            region_narrowed_for,
        } = self;

        *region_narrowed_for = Some(observation.region);
        candidate_of_merge_allele.clear();
        candidate_of_merge_allele
            .extend((0..remap.table_len()).map(|allele| remap.candidate_for(allele)));

        // `resize_with` truncates as well as extends, so nothing follows it.
        each_run_sample.resize_with(run_sample_count, NarrowedRunSample::default);
        // **Every entry reset, not only the covering samples'.** A sample that covered nothing
        // at this locus may have covered the last one, and its buffer would still hold that
        // locus's rows. **They would not be *scored***, and the first version of this comment
        // said they would: `fill_views` reads a sample's rows only in the covering arm, so a
        // stale row is unread rather than wrong. What the reset buys is that the buffer says
        // what the locus says — and `rows_left_for` is how a test sees it, because no other
        // route can.
        for sample in &mut *each_run_sample {
            let NarrowedRunSample {
                rows,
                merge_entry,
                leftover,
            } = sample;
            rows.clear();
            *merge_entry = None;
            *leftover = SampleLeftover::default();
        }

        let mut previous_run_sample: Option<usize> = None;
        for (merge_entry, support) in observation.per_sample.iter().enumerate() {
            let run_sample = support.sample;
            assert!(
                run_sample < run_sample_count,
                "the merge's entry {merge_entry} at {} names sample {run_sample} of a run of \
                 {run_sample_count}",
                observation.region
            );
            assert!(
                previous_run_sample.is_none_or(|previous| previous < run_sample),
                "the merge's covering samples are in ascending sample order and entry \
                 {merge_entry} at {} names sample {run_sample} after {:?}",
                observation.region,
                previous_run_sample
            );
            previous_run_sample = Some(run_sample);

            let narrowed = &mut each_run_sample[run_sample];
            let dropped_by_the_narrowing = GenericObservation::fill_from_supported_alleles(
                &support.supported,
                candidate_of_merge_allele,
                &mut narrowed.rows,
            );
            // **Ascending on the candidate key — and today this sorts nothing.** Selection puts
            // its survivors back into the merge table's index order before admitting, so the
            // remapping is order-preserving and the rows arrive sorted. It is here as insurance
            // against that contract changing: `AlleleRemap::admit` constrains only the candidate
            // ids, so a rank-order selection would be type-legal, and
            // `GenericSampleEvidence::new`'s order check is a `debug_assert` that a release run
            // would not raise (`spec/read_likelihoods.md` §8).
            narrowed
                .rows
                .sort_unstable_by_key(|row| (row.allele, row.read_group));

            let leftover = selection.unmatched()[merge_entry];
            let disagreement = (leftover.q_sum - dropped_by_the_narrowing).abs();
            let scale = leftover
                .q_sum
                .abs()
                .max(dropped_by_the_narrowing.abs())
                .max(1.0);
            assert!(
                disagreement <= LEFTOVER_AGREEMENT_TOLERANCE * scale,
                "sample {run_sample}'s dropped reads pool to {} by selection's own sum and to \
                 {dropped_by_the_narrowing} by this narrowing's, a disagreement of \
                 {disagreement} over what should be the same rows — too large to be the two \
                 walks' addition order, so they disagree about which of the merge's alleles \
                 have no candidate",
                leftover.q_sum
            );
            narrowed.merge_entry = Some(merge_entry);
            narrowed.leftover = SampleLeftover {
                // **Selection's, not the sum of the two.** They are Σ over the same rows; see
                // the module's own note.
                unmatched_q_sum: leftover.q_sum,
                genotype_must_be_missing: leftover.genotype_must_be_missing(),
            };
        }
    }

    /// The rows this narrowing left for one run sample.
    ///
    /// **Test-only, and it exists so that the buffer-reuse test can fail.**
    /// [`fill_views`](Self::fill_views) reads a sample's rows only when that sample covered the
    /// locus, so a stale row left in a non-covering sample's buffer is invisible by every other
    /// route — and a test asserting on the views alone passes whether the reset happened or not.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn rows_left_for(&self, run_sample: usize) -> &[GenericObservation] {
        &self.each_run_sample[run_sample].rows
    }

    /// **The per-sample views the loop reads**, one per sample of the run, filled into `out`.
    ///
    /// `out` is the caller's and is cleared here. It is the shaping's one per-locus
    /// allocation, and the module's own note says why it cannot live in this type.
    ///
    /// # Panics
    ///
    /// If `observation` is not the one [`narrow`](Self::narrow) was last called with, **checked
    /// by the locus's own region** rather than by a count of anything: the rows are per run
    /// sample and the partials are borrowed from the *merge's* entries, so filling against a
    /// different locus pairs one locus's narrowed reads with another's partials — both legal
    /// evidence, and neither the caller's. It also refuses a scratch that was never narrowed at
    /// all, which a count cannot.
    pub fn fill_views<'a>(
        &'a self,
        observation: &'a CohortObservation,
        out: &mut Vec<GenericLocusSample<'a>>,
    ) {
        assert_eq!(
            self.region_narrowed_for,
            Some(observation.region),
            "these buffers were narrowed for {:?} and the views are being filled for {}",
            self.region_narrowed_for,
            observation.region
        );
        out.clear();
        out.reserve(self.each_run_sample.len());
        for (run_sample, narrowed) in self.each_run_sample.iter().enumerate() {
            let leftover = narrowed.leftover;
            let evidence = match narrowed.merge_entry {
                Some(merge_entry) => {
                    // **The join itself, at every sample.** The rows are this run sample's and
                    // the partials come from the merge's entry, so the two are one sample's
                    // evidence only if that entry still names this sample. The region check
                    // above rules out a different locus; this rules out an observation of the
                    // same locus whose covering list was rebuilt.
                    assert_eq!(
                        observation.per_sample[merge_entry].sample, run_sample,
                        "the buffers were narrowed with covering entry {merge_entry} as run \
                         sample {run_sample}, and the observation at {} names sample {} there",
                        observation.region, observation.per_sample[merge_entry].sample
                    );
                    GenericSampleEvidence::new(
                        &narrowed.rows,
                        leftover.unmatched_q_sum,
                        &observation.per_sample[merge_entry].partials,
                    )
                }
                // **A sample that covered nothing, and this is not a special case**: an empty
                // sum is zero, so every genotype scores alike and the prior decides alone
                // (`doc/devel/ng/spec/read_likelihoods.md` §3.3).
                None => GenericSampleEvidence::empty(),
            };
            out.push(GenericLocusSample {
                evidence,
                genotype_must_be_missing: leftover.genotype_must_be_missing,
            });
        }
    }
}

/// **One merged locus as the loop reads it**, over selection's candidates and the run's sample
/// order — the whole of the input edge in one call.
///
/// The two halves are separate because the views borrow the narrowed rows, and a caller that
/// wants to hold the evidence has to hold the buffers too. This is the ordinary spelling.
///
/// **`views` must be a fresh `Vec` at each locus, declared inside the loop body**, and that is
/// not this function's signature imposing it: the element type names the lifetime of the
/// buffers the views borrow, and a `Vec` is invariant in its element type — so a list held
/// across two loci holds the first locus's borrow of `shaping` open into the second, which is
/// an `E0499` here and an `E0502` on the two-step form. Calling
/// [`narrow`](GenericEvidenceScratch::narrow) and
/// [`fill_views`](GenericEvidenceScratch::fill_views) separately does not lift it. That fresh
/// `Vec` is the per-locus allocation the module note names.
#[must_use]
pub fn shape_generic_locus<'a>(
    shaping: &'a mut GenericEvidenceScratch,
    observation: &'a CohortObservation,
    selection: &LocusSelection,
    run_sample_count: usize,
    views: &'a mut Vec<GenericLocusSample<'a>>,
) -> LocusEvidence<'a> {
    shaping.narrow(observation, selection, run_sample_count);
    shaping.fill_views(observation, views);
    LocusEvidence::generic(observation.region, views)
}

/// **One repeat tract's evidence, per sample of the run** — the repeat-tract half of the input
/// edge.
///
/// **Almost nothing to do, and that is the design rather than an accident.** The STR generator
/// already keys its observations on `(bases, witness, read group)` and carries a read count, so
/// the aggregation contract the row needs holds by construction and there is nothing to reshape
/// (`doc/devel/ng/arch/calling_em_loop.md` §2.2). What this adds is the run's sample order: the
/// generator produces one list per sample it saw, and the loop reads one entry per sample of the
/// run.
///
/// `observations_of_each_run_sample` is that list, **one entry per run sample**, empty where the
/// sample showed nothing. Unlike the SNP/indel path there is no covering-samples list to join
/// against — the caller holds the run's samples already — and **no sample is set aside at a
/// tract**: a discovery round there can put back a length the cap cut, so nobody is locked out of
/// the locus for the rest of its calling (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1).
///
/// `candidate_repeat_counts` is the one per-locus input a tract needs beyond its candidates, and
/// it travels because **it is not derivable from the candidates' bases**: an interrupted tract
/// holds fewer whole repeats than its length suggests, and the slippage fit is keyed by the count
/// (`doc/devel/ng/spec/read_likelihoods.md` §4.4). Its producer is repeat-tract candidate
/// selection, which is unwritten, so today every caller supplies it.
///
/// # Panics
///
/// If either list is empty — by [`LocusEvidence::ssr`], which owns both refusals. A run has at
/// least one sample and every sample gets an entry, so an empty list is evidence that went
/// missing rather than a tract nobody covered, which is one empty entry per sample.
#[must_use]
pub fn shape_ssr_locus<'a>(
    region: GenomeRegion,
    observations_of_each_run_sample: &'a [&'a [SequenceObservation]],
    detail: &'a SsrDetail,
    candidate_repeat_counts: &'a [NonZeroU32],
    views: &'a mut Vec<SsrSampleEvidence<'a>>,
) -> LocusEvidence<'a> {
    // **The emptiness check is [`LocusEvidence::ssr`]'s**, not restated here: a run has at
    // least one sample and every sample gets an entry, so an empty list is evidence that went
    // missing rather than a tract nobody covered — and that constructor refuses it by name.
    views.clear();
    views.reserve(observations_of_each_run_sample.len());
    views.extend(
        observations_of_each_run_sample
            .iter()
            .map(|observations| SsrSampleEvidence::new(observations, detail)),
    );
    LocusEvidence::ssr(region, views, detail, candidate_repeat_counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::ng::calling::CandidateAlleles;
    use crate::ng::calling::allele_candidates::{AlleleRemap, SelectionVerdict, UnmatchedSupport};
    use crate::ng::locus_generation::LocusKind;
    use crate::ng::run::cohort_merge::build::{AlleleSupport, SampleSupport, SupportedAllele};
    use crate::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError};

    /// Two candidates' repeat counts, different from each other — the shape a tract's evidence
    /// carries beside its samples until repeat-tract selection produces them.
    ///
    /// **Supplied, not selected**: the repeat-tract half of candidate selection is unwritten, so
    /// a fixture states its candidates' repeat counts and a reader must not take them for a
    /// step's output.
    fn repeat_counts() -> Vec<NonZeroU32> {
        vec![
            NonZeroU32::new(6).expect("six repeats"),
            NonZeroU32::new(7).expect("seven repeats"),
        ]
    }

    fn region() -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(3),
            start: Position(940),
            end: Position(940),
        }
    }

    /// One merge row: which of the merge's alleles, from which read group, and its reads.
    fn merge_row(allele: usize, read_group: u32, num_reads: u32, q_sum: f64) -> SupportedAllele {
        SupportedAllele {
            allele,
            read_group: ReadGroupId(read_group),
            support: AlleleSupport {
                num_reads,
                num_fwd: num_reads / 2,
                q_sum,
                mapq_sum: 60 * num_reads,
                mapq_sum_sq: u64::from(num_reads) * 3_600,
                placed_left: num_reads / 2,
            },
        }
    }

    /// One covering sample's merge entry.
    fn covering(sample: usize, supported: Vec<SupportedAllele>) -> SampleSupport {
        SampleSupport {
            sample,
            supported,
            partials: Vec::new(),
            reads_composed_across_records: 0,
            reads_removed_as_evidence: 0,
            reads_without_observation: 0,
        }
    }

    /// A merged locus over `alleles` many sequences, covered by `per_sample`.
    fn merged_locus(alleles: usize, per_sample: Vec<SampleSupport>) -> CohortObservation {
        let sequences = [b"A".as_slice(), b"T", b"C", b"G"];
        CohortObservation {
            region: region(),
            alleles: sequences[..alleles].iter().map(|s| Box::from(*s)).collect(),
            per_sample,
        }
    }

    /// What selection would hand back: `kept` maps each merge allele to its candidate id, and
    /// `unmatched` is **parallel to the merge's covering samples**.
    fn selection_of(kept: &[Option<u16>], unmatched: Vec<UnmatchedSupport>) -> LocusSelection {
        let mut alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
        let mut remap = AlleleRemap::with_all_dropped(kept.len());
        let mut admitted = vec![None; kept.len()];
        // The reference is never dropped and is always candidate 0.
        remap.admit(0, AlleleId::REFERENCE);
        admitted[0] = Some(0_u16);
        // Admit the alternatives in *candidate id* order, which is selection's ranking order
        // and need not be the merge's.
        let mut by_candidate: Vec<(u16, usize)> = kept
            .iter()
            .enumerate()
            .skip(1)
            .filter_map(|(table_index, candidate)| candidate.map(|id| (id, table_index)))
            .collect();
        by_candidate.sort_unstable();
        for (candidate, table_index) in by_candidate {
            let sequences = [b"A".as_slice(), b"T", b"C", b"G"];
            let minted = alleles.admit(Box::from(sequences[table_index]));
            assert_eq!(minted.get(), candidate, "the fixture's ids are consistent");
            remap.admit(table_index, minted);
            admitted[table_index] = Some(candidate);
        }
        let covering_samples = unmatched.len();
        LocusSelection::new(
            alleles,
            SelectionVerdict::Selected,
            unmatched,
            remap,
            covering_samples,
        )
    }

    fn leftover(q_sum: f64, earned_reads_cut_by_the_cap: u32) -> UnmatchedSupport {
        UnmatchedSupport {
            num_reads: 0,
            q_sum,
            earned_reads_cut_by_the_cap,
        }
    }

    /// **A sample that covered nothing gets an empty view, and it keeps its place in the run's
    /// order.**
    ///
    /// The merge's list holds covering samples only, each naming its own run index; the loop's
    /// holds one entry per sample of the run. A join by position rather than by that index
    /// would give sample 1 sample 2's reads — and nothing downstream would object, because both
    /// are legal evidence.
    #[test]
    fn a_sample_that_covered_nothing_keeps_its_place_and_shows_nothing() {
        let observation = merged_locus(
            2,
            vec![
                covering(0, vec![merge_row(0, 0, 4, -4.0)]),
                covering(2, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0); 2]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 3);
        let mut views = Vec::new();
        shaping.fill_views(&observation, &mut views);

        assert_eq!(views.len(), 3, "one entry per sample of the run");
        assert_eq!(views[0].evidence.supported.len(), 1);
        assert_eq!(views[0].evidence.supported[0].allele, AlleleId::REFERENCE);
        assert!(
            views[1].evidence.supported.is_empty() && views[1].evidence.partials.is_empty(),
            "sample 1 covered nothing, so it shows nothing — and it is still sample 1"
        );
        assert_eq!(views[2].evidence.supported.len(), 1);
        assert_eq!(views[2].evidence.supported[0].allele, AlleleId(1));
    }

    /// **The narrowed rows come back ascending on the candidate key even when the remapping
    /// permutes** — the insurance case, and **a shape today's selection cannot produce.**
    ///
    /// `select_generic` puts its survivors back into the merge table's index order before
    /// admitting, so on every input the shipped producer can build this sort has nothing to do.
    /// The fixture hand-builds what a rank-order selection would hand over — merge allele 1 as
    /// candidate 2, merge allele 2 as candidate 1 — because `AlleleRemap::admit` constrains only
    /// the candidate ids and would accept it. **So this test pins the insurance, not today's
    /// behaviour**, and saying which it is matters: the read likelihood's sum must run in a
    /// fixed order (`spec/read_likelihoods.md` §8), and `GenericSampleEvidence::new` checks that
    /// only in debug.
    #[test]
    fn the_narrowed_rows_are_sorted_on_the_candidate_key_and_not_the_merges() {
        let observation = merged_locus(
            3,
            vec![covering(
                0,
                vec![
                    merge_row(0, 0, 4, -4.0),
                    merge_row(1, 0, 5, -5.0),
                    merge_row(2, 0, 6, -6.0),
                ],
            )],
        );
        // Merge allele 1 → candidate 2, merge allele 2 → candidate 1: the ranking reversed
        // the two alternatives.
        let selection = selection_of(&[Some(0), Some(2), Some(1)], vec![leftover(0.0, 0)]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 1);
        let mut views = Vec::new();
        shaping.fill_views(&observation, &mut views);

        let rows = views[0].evidence.supported;
        assert_eq!(
            rows.iter().map(|row| row.allele.get()).collect::<Vec<_>>(),
            vec![0, 1, 2],
            "ascending on the candidate id"
        );
        // And the sort really moved something: candidate 1 is the merge's allele 2, whose six
        // reads now sit before candidate 2's five.
        assert_eq!(
            rows.iter().map(|row| row.num_reads).collect::<Vec<_>>(),
            vec![4, 6, 5],
            "candidate 1 is the merge's third sequence, which showed six reads"
        );
    }

    /// **The sort runs for every covering sample, not only the first** — and it takes a sample
    /// at a *shifted* run index with *two* alternatives to say so.
    ///
    /// A permutation of the allele mapping is distinguishable from a relabelling only where a
    /// sample carries at least two alternatives, and a per-sample loop body is pinned only where
    /// at least two samples run through it at a shifted index. Measured on the module's other
    /// fixtures: restricting the sort to the first covering entry left every one of them green.
    #[test]
    fn a_shifted_covering_sample_with_two_alternatives_is_sorted_like_the_first() {
        let rows = vec![
            merge_row(0, 0, 4, -4.0),
            merge_row(1, 0, 5, -5.0),
            merge_row(2, 0, 6, -6.0),
        ];
        let observation = merged_locus(3, vec![covering(1, rows.clone()), covering(3, rows)]);
        // Merge allele 1 → candidate 2, merge allele 2 → candidate 1.
        let selection = selection_of(&[Some(0), Some(2), Some(1)], vec![leftover(0.0, 0); 2]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 4);
        let mut views = Vec::new();
        shaping.fill_views(&observation, &mut views);

        for run_sample in [1, 3] {
            let ids: Vec<u16> = views[run_sample]
                .evidence
                .supported
                .iter()
                .map(|row| row.allele.get())
                .collect();
            assert_eq!(
                ids,
                vec![0, 1, 2],
                "sample {run_sample}'s rows are ascending on the candidate id"
            );
            assert_eq!(
                views[run_sample].evidence.supported[1].num_reads, 6,
                "candidate 1 is the merge's third sequence, which showed six reads"
            );
        }
    }

    /// **Two read groups of one allele stay apart and stay ordered**, because the likelihood
    /// pools an observation's reads into one term only if every one of them would get the same
    /// number — and two lanes have different error rates.
    #[test]
    fn two_read_groups_of_one_allele_stay_apart_and_in_group_order() {
        let observation = merged_locus(
            2,
            vec![covering(
                0,
                vec![
                    merge_row(0, 0, 4, -4.0),
                    merge_row(0, 1, 3, -3.0),
                    merge_row(1, 1, 6, -6.0),
                    merge_row(1, 0, 5, -5.0),
                ],
            )],
        );
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0)]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 1);
        let mut views = Vec::new();
        shaping.fill_views(&observation, &mut views);

        let key: Vec<(u16, u32)> = views[0]
            .evidence
            .supported
            .iter()
            .map(|row| (row.allele.get(), row.read_group.get()))
            .collect();
        assert_eq!(key, vec![(0, 0), (0, 1), (1, 0), (1, 1)]);
    }

    /// **The pooled leftover is selection's own number, and the narrowing's is a check on it
    /// rather than a second helping.**
    ///
    /// Both are Σ `ln P(error)` over exactly the rows whose allele has no candidate. Adding
    /// them would double the mass every genotype of that sample is charged — the same number
    /// under every genotype, so it cancels in the *genotype* comparison and does not cancel in
    /// the data likelihood, which emission and the site quality read.
    #[test]
    fn the_leftover_is_selections_number_and_not_twice_it() {
        let observation = merged_locus(
            3,
            vec![covering(
                0,
                vec![
                    merge_row(0, 0, 4, -4.0),
                    merge_row(1, 0, 5, -5.0),
                    // Merge allele 2 is dropped: its reads are the leftover.
                    merge_row(2, 0, 2, -7.5),
                ],
            )],
        );
        let selection = selection_of(&[Some(0), Some(1), None], vec![leftover(-7.5, 0)]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 1);
        let mut views = Vec::new();
        shaping.fill_views(&observation, &mut views);

        assert_eq!(views[0].evidence.supported.len(), 2, "one dropped row");
        assert_eq!(
            views[0].evidence.unmatched_q_sum, -7.5,
            "selection's pool, not -15.0"
        );
    }

    /// Two independently written walks over one set of rows must reach the same total. A
    /// disagreement means one of them changed its mind about which of the merge's alleles have
    /// no candidate, and the loop would score every sample against a leftover wrong by exactly
    /// those rows.
    #[test]
    #[should_panic(expected = "they disagree about which of the merge's alleles")]
    fn a_leftover_that_disagrees_with_the_narrowings_own_sum_is_refused() {
        let observation = merged_locus(
            3,
            vec![covering(
                0,
                vec![merge_row(0, 0, 4, -4.0), merge_row(2, 0, 2, -7.5)],
            )],
        );
        // Selection says it dropped 3 nats' worth; the rows say 7.5.
        let selection = selection_of(&[Some(0), Some(1), None], vec![leftover(-3.0, 0)]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 1);
    }

    /// **The uncallable ruling travels by the merge's covering order and lands in the run's.**
    ///
    /// Selection's leftovers are parallel to the merge's covering samples; the flag has to
    /// reach the *run* sample that entry names. Here the second covering entry is run sample 2,
    /// and it is the one whose earned allele the cap cut.
    #[test]
    fn the_uncallable_ruling_lands_on_the_run_sample_the_merge_entry_names() {
        let observation = merged_locus(
            2,
            vec![
                covering(0, vec![merge_row(0, 0, 4, -4.0)]),
                covering(2, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        let selection = selection_of(
            &[Some(0), Some(1)],
            vec![leftover(0.0, 0), leftover(0.0, 9)],
        );
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 3);
        let mut views = Vec::new();
        shaping.fill_views(&observation, &mut views);

        assert!(!views[0].genotype_must_be_missing);
        assert!(
            !views[1].genotype_must_be_missing,
            "a sample that covered nothing lost nothing"
        );
        assert!(
            views[2].genotype_must_be_missing,
            "the second covering entry is run sample 2"
        );
    }

    /// A sample's partial reads travel to the view untouched — they are borrowed from the
    /// merge's own rows rather than copied, because a partial's bases and witness are what the
    /// censored term reads and nothing here reshapes them.
    #[test]
    fn the_partial_observations_are_borrowed_from_the_merges_own_rows() {
        use crate::ng::locus_generation::WitnessedLocusPositions;
        use crate::ng::run::cohort_merge::build::PartialObservation;

        let mut support = covering(0, vec![merge_row(0, 0, 4, -4.0)]);
        support.partials.push(PartialObservation {
            witnessed_in_locus: WitnessedLocusPositions::from_half_open_runs([(0_u16, 1_u16)])
                .expect("one witnessed position"),
            read_group: ReadGroupId(0),
            bases: Box::from(b"A".as_slice()),
            num_reads: 3,
            q_sum: -9.0,
        });
        let observation = merged_locus(2, vec![support]);
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0)]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 1);
        let mut views = Vec::new();
        shaping.fill_views(&observation, &mut views);

        assert_eq!(views[0].evidence.partials.len(), 1);
        assert_eq!(views[0].evidence.partials[0].num_reads, 3);
        assert!(
            std::ptr::eq(
                views[0].evidence.partials.as_ptr(),
                observation.per_sample[0].partials.as_ptr()
            ),
            "borrowed, not copied"
        );
    }

    /// **A worker's buffers are reused, and a sample that covered the last locus but not this
    /// one must show nothing.** A stale row would be the previous locus's reads scored at this
    /// one — finite, plausible and wrong.
    #[test]
    fn a_sample_that_covered_the_last_locus_but_not_this_one_shows_nothing() {
        let mut shaping = GenericEvidenceScratch::default();

        let first = merged_locus(2, vec![covering(1, vec![merge_row(1, 0, 6, -6.0)])]);
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0)]);
        shaping.narrow(&first, &selection, 2);
        let mut views = Vec::new();
        shaping.fill_views(&first, &mut views);
        assert_eq!(views[1].evidence.supported.len(), 1);

        let mut second = merged_locus(2, vec![covering(0, vec![merge_row(0, 0, 4, -4.0)])]);
        second.region.start = Position(1_500);
        second.region.end = Position(1_500);
        shaping.narrow(&second, &selection, 2);
        // **Asserted on the buffer, not on the view.** `fill_views` reads a sample's rows only
        // when that sample covered the locus, so a view alone comes back empty whether the
        // buffer was reset or not — measured: deleting the whole reset loop leaves every other
        // test in this module green.
        assert!(
            shaping.rows_left_for(1).is_empty(),
            "sample 1's buffer still holds the first locus's row"
        );
        let mut views = Vec::new();
        shaping.fill_views(&second, &mut views);
        assert!(
            views[1].evidence.supported.is_empty(),
            "sample 1 covered the first locus and not the second"
        );
    }

    /// **The views must be filled for the locus the buffers were narrowed for**, and the check
    /// is on the locus's own region rather than on a count of anything: the row buffers are per
    /// run sample and the partials are borrowed from the *merge's* entries, so filling against
    /// a different locus pairs one locus's narrowed reads with another's partials — both legal
    /// evidence, and neither the caller's.
    ///
    /// **A count would let the likely half through.** Two *different* loci with the same number
    /// of covering samples are the common case at a cohort where most samples cover most loci,
    /// not the corner — so this fixture's second locus differs only in where it is.
    #[test]
    #[should_panic(expected = "the views are being filled for")]
    fn views_filled_for_a_different_locus_than_the_buffers_hold_are_refused() {
        let narrowed = merged_locus(
            2,
            vec![
                covering(0, vec![merge_row(0, 0, 4, -4.0)]),
                covering(1, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0); 2]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&narrowed, &selection, 2);

        // The same shape, the same covering count, a different place.
        let mut elsewhere = merged_locus(
            2,
            vec![
                covering(0, vec![merge_row(0, 0, 4, -4.0)]),
                covering(1, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        elsewhere.region.start = Position(1_500);
        elsewhere.region.end = Position(1_500);
        let mut views = Vec::new();
        shaping.fill_views(&elsewhere, &mut views);
    }

    /// Filling views from a scratch that was never narrowed at all is refused — the failure a
    /// count of covering samples could not see, because a fresh scratch and a locus nobody
    /// covered agree at zero.
    #[test]
    #[should_panic(expected = "the views are being filled for")]
    fn views_filled_from_a_scratch_that_was_never_narrowed_are_refused() {
        let observation = merged_locus(2, vec![]);
        let shaping = GenericEvidenceScratch::default();
        let mut views = Vec::new();
        shaping.fill_views(&observation, &mut views);
    }

    /// **A worker's loop over loci compiles, with the view list declared inside the body** —
    /// which is what "callable" has to mean, and what a single call outside a loop does not
    /// say.
    ///
    /// The list's element type names the lifetime of the buffers the views borrow, and a `Vec`
    /// is invariant in its element type, so a list hoisted out of the loop holds the first
    /// locus's borrow of the scratch open into the second. Compiled and rejected while this
    /// test was written: hoisting `views` is an `E0499` here and an `E0502` on the two-step
    /// form. **That caller now exists** — `tests/ng_calling_loop_calls_genotypes.rs`, which runs
    /// this shape end to end over a cohort locus and asserts the genotypes that come out.
    #[test]
    fn a_worker_shapes_one_locus_after_another_on_one_scratch() {
        let mut first = merged_locus(
            2,
            vec![
                covering(0, vec![merge_row(0, 0, 4, -4.0)]),
                covering(2, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        first.region.start = Position(100);
        first.region.end = Position(100);
        let second = merged_locus(2, vec![covering(1, vec![merge_row(1, 0, 5, -5.0)])]);
        let selection_of_first = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0); 2]);
        let selection_of_second = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0)]);

        let mut shaping = GenericEvidenceScratch::default();
        let mut regions = Vec::new();
        for (observation, selection) in [
            (&first, &selection_of_first),
            (&second, &selection_of_second),
        ] {
            // **Inside the body, and it has to be.**
            let mut views = Vec::new();
            let evidence = shape_generic_locus(&mut shaping, observation, selection, 3, &mut views);
            assert_eq!(evidence.sample_count(), 3);
            regions.push(evidence.region());
        }
        assert_eq!(regions, vec![first.region, second.region]);
    }

    /// A run has at least one sample, so a locus shaped for none is a run whose sample order
    /// went missing — refused where the count arrives rather than two calls later, where the
    /// message would be about something else.
    #[test]
    #[should_panic(expected = "a run has at least one sample")]
    fn a_locus_shaped_for_a_run_of_no_samples_is_refused() {
        let observation = merged_locus(2, vec![]);
        let selection = selection_of(&[Some(0), Some(1)], Vec::new());
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 0);
    }

    /// **A repeat tract's evidence is the generator's own observations, in the run's order** —
    /// and **no sample is set aside there**, however little it showed
    /// (`spec/calling_em_loop.md` §5.0.1).
    #[test]
    fn a_tract_takes_the_generators_observations_and_sets_no_sample_aside() {
        use crate::ng::locus_generation::{ReadWitness, SequenceObservation, SsrDetail};
        use crate::ng::types::Motif;

        let detail = SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        };
        let seen = [SequenceObservation {
            bases: b"ATAT".to_vec(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: 5,
            num_fwd: 3,
            q_sum: SummedLogError::from_nats(-5.0),
            mapq_sum: 300,
            mapq_sum_sq: 18_000,
            placed_left: 2,
            chain_ids: Vec::new(),
        }];
        let per_run_sample: [&[SequenceObservation]; 2] = [&seen, &[]];
        let mut views = Vec::new();
        let counts = repeat_counts();
        let evidence = shape_ssr_locus(region(), &per_run_sample, &detail, &counts, &mut views);

        assert_eq!(
            evidence.sample_count(),
            2,
            "one entry per sample of the run"
        );
        match evidence {
            LocusEvidence::Ssr { per_sample, .. } => {
                assert_eq!(per_sample[0].observations.len(), 1);
                assert!(
                    per_sample[1].observations.is_empty(),
                    "a sample that showed nothing at the tract is still a sample of the run"
                );
            }
            LocusEvidence::Generic { .. } => panic!("a tract is not a SNP/indel locus"),
        }
    }

    /// A run has at least one sample, so an empty list is evidence that went missing rather
    /// than a tract nobody covered — which is one empty entry per sample.
    #[test]
    #[should_panic(expected = "names no sample")]
    fn a_tract_whose_evidence_names_no_sample_is_refused() {
        use crate::ng::types::Motif;
        let detail = SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        };
        let mut views = Vec::new();
        let _ = shape_ssr_locus(region(), &[], &detail, &repeat_counts(), &mut views);
    }

    /// **A sample ruled uncallable at one locus is callable again at the next** — the reset the
    /// whole per-sample entry gets, seen from the one place a non-covering sample's leftover is
    /// still visible.
    ///
    /// `fill_views` writes `genotype_must_be_missing` for **every** sample, covering or not, so
    /// a leftover that survived into the next locus emits a missing genotype at a locus where
    /// the sample simply showed nothing. Measured on the version that resized instead of
    /// resetting: `Vec::resize` fills only the slots it adds, so the flag stayed and the whole
    /// suite was green.
    #[test]
    fn a_sample_ruled_uncallable_at_one_locus_is_callable_at_the_next() {
        let mut shaping = GenericEvidenceScratch::default();
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 7)]);

        let first = merged_locus(2, vec![covering(1, vec![merge_row(1, 0, 6, -6.0)])]);
        shaping.narrow(&first, &selection, 2);
        let mut views = Vec::new();
        shaping.fill_views(&first, &mut views);
        assert!(
            views[1].genotype_must_be_missing,
            "the cap cut a sequence sample 1 had earned at the first locus"
        );

        // The next locus, which sample 1 does not cover at all.
        let mut second = merged_locus(2, vec![covering(0, vec![merge_row(0, 0, 4, -4.0)])]);
        second.region.start = Position(1_500);
        second.region.end = Position(1_500);
        let one_sample = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0)]);
        shaping.narrow(&second, &one_sample, 2);
        let mut views = Vec::new();
        shaping.fill_views(&second, &mut views);
        assert!(
            !views[1].genotype_must_be_missing,
            "sample 1 showed nothing here, which is not the same as being ruled uncallable"
        );
    }

    /// **The allele mapping is the locus's own, and the previous locus's is not left under it.**
    ///
    /// The remapping is rebuilt per locus into one buffer; extending without clearing would
    /// leave the first locus's entries at indices the second locus's rows also name, so a row
    /// would be narrowed onto whatever allele the *last* locus put there. Here the two loci map
    /// merge allele 1 to different candidates, so the stale mapping is visible in the ids that
    /// come back.
    #[test]
    fn the_allele_mapping_is_this_locus_own_and_not_the_last_ones() {
        let mut shaping = GenericEvidenceScratch::default();

        let first = merged_locus(3, vec![covering(0, vec![merge_row(1, 0, 5, -5.0)])]);
        // Merge allele 1 is candidate 2 here.
        let first_selection = selection_of(&[Some(0), Some(2), Some(1)], vec![leftover(0.0, 0)]);
        shaping.narrow(&first, &first_selection, 1);
        let mut views = Vec::new();
        shaping.fill_views(&first, &mut views);
        assert_eq!(views[0].evidence.supported[0].allele, AlleleId(2));

        // The same allele at the next locus is candidate 1.
        let mut second = merged_locus(3, vec![covering(0, vec![merge_row(1, 0, 5, -5.0)])]);
        second.region.start = Position(1_500);
        second.region.end = Position(1_500);
        let second_selection = selection_of(&[Some(0), Some(1), Some(2)], vec![leftover(0.0, 0)]);
        shaping.narrow(&second, &second_selection, 1);
        let mut views = Vec::new();
        shaping.fill_views(&second, &mut views);
        assert_eq!(
            views[0].evidence.supported[0].allele,
            AlleleId(1),
            "the second locus's own mapping, not the first's"
        );
    }

    /// **The covering entry a row was narrowed under must still name the same run sample.**
    ///
    /// The region check rules out a *different* locus. This rules out an observation of the
    /// same locus whose covering list was rebuilt — where the rows are one sample's and the
    /// partials another's, both legal evidence and neither the caller's.
    #[test]
    #[should_panic(expected = "and the observation at")]
    fn views_filled_against_a_rebuilt_covering_list_are_refused() {
        let narrowed = merged_locus(
            2,
            vec![
                covering(0, vec![merge_row(0, 0, 4, -4.0)]),
                covering(1, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0); 2]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&narrowed, &selection, 3);

        // The same region and the same covering count, and the second entry is a different
        // sample.
        let rebuilt = merged_locus(
            2,
            vec![
                covering(0, vec![merge_row(0, 0, 4, -4.0)]),
                covering(2, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        let mut views = Vec::new();
        shaping.fill_views(&rebuilt, &mut views);
    }

    /// A remapping built against a different locus's allele table is refused: it is indexed by
    /// the merge's allele indices, so it would map a row onto whatever allele sits there.
    #[test]
    #[should_panic(expected = "belong to different loci")]
    fn a_remapping_of_another_locus_is_refused() {
        let observation = merged_locus(3, vec![covering(0, vec![merge_row(0, 0, 4, -4.0)])]);
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0)]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 1);
    }

    /// Selection's leftovers run parallel to the merge's covering samples and to nothing else.
    #[test]
    #[should_panic(expected = "run parallel to the merge's covering samples")]
    fn leftovers_that_are_not_one_per_covering_sample_are_refused() {
        let observation = merged_locus(
            2,
            vec![
                covering(0, vec![merge_row(0, 0, 4, -4.0)]),
                covering(1, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0)]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 2);
    }

    /// A covering sample naming an index past the run is refused rather than dropped.
    #[test]
    #[should_panic(expected = "of a run of")]
    fn a_covering_sample_past_the_end_of_the_run_is_refused() {
        let observation = merged_locus(2, vec![covering(5, vec![merge_row(0, 0, 4, -4.0)])]);
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0)]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 2);
    }

    /// The merge's covering samples are in ascending sample order, and the join to selection's
    /// leftovers is positional — so an unsorted list would pair one sample's dropped reads with
    /// another's.
    #[test]
    #[should_panic(expected = "ascending sample order")]
    fn covering_samples_out_of_order_are_refused() {
        let observation = merged_locus(
            2,
            vec![
                covering(1, vec![merge_row(0, 0, 4, -4.0)]),
                covering(0, vec![merge_row(1, 0, 6, -6.0)]),
            ],
        );
        let selection = selection_of(&[Some(0), Some(1)], vec![leftover(0.0, 0); 2]);
        let mut shaping = GenericEvidenceScratch::default();
        shaping.narrow(&observation, &selection, 2);
    }
}
