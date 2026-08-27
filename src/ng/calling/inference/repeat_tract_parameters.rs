//! What one repeat tract is scored under — the assembly between the run's fitted parameters
//! and the STR row.
//!
//! The repeat-tract row takes a **scoring context per `(read group, candidate)`**
//! ([`SsrScoringContext`]), a length support for its junk term, and one outlier weight
//! (`doc/devel/ng/spec/read_likelihoods.md` §4). Every one of those numbers already exists
//! somewhere — the slippage fit, the substitution-rate map, the stutter distribution's own
//! reachability rule — and until this module nothing outside `likelihood/ssr.rs`'s own tests
//! had ever put them together. That gap is why a tract could not be genotyped while its row
//! was shipped and merged.
//!
//! # What is looked up per cell, and why it cannot be looked up per locus
//!
//! **A read's chance of slipping is a property of the tract it was copied from, and that is
//! the candidate allele** (spec §4.4). A candidate of 6 repeats and one of 12 at the same
//! tract are drawn from different strata and slip at measurably different rates, so the
//! stutter parameters cannot be hoisted out of the candidate loop; and slippage is a property
//! of the chemistry, so they cannot be hoisted out of the read-group loop either. The table is
//! `read groups × candidates` and that shape is the model rather than a convenience.
//!
//! # The three answers this module owes, and what it answers
//!
//! **Two of the three fitted numbers a scoring context carries can be missing on perfectly good
//! data, and the outlier weight — which is not on a context at all, but beside them on
//! [`SsrLocusParameters`] — is not fitted anywhere.** Each is answered with a **stated constant
//! and a warrant that says so**, never with a silent zero:
//!
//! - **A candidate whose stratum the slippage fit has no numbers for.**
//!   [`NoSlippage`]'s own documentation says a caller has to have an answer, because *"a
//!   candidate several repeats from its reference tract's length can land here on perfectly
//!   good data"*. It gets [`StutterModel::hipstr_shipped`] and [`Provenance::Defaulted`]. The
//!   four absences it names are not alike, and the two that mean *the run is not what it
//!   claims* are counted apart
//!   ([`cells_whose_read_group_the_fit_does_not_describe`](TractScoringFits::cells_whose_read_group_the_fit_does_not_describe)).
//! - **A `(read group, candidate stratum)` pair the substitution-rate fit has no entry for** —
//!   ordinarily because the candidate's stratum is not in the fit at all, the same reason the
//!   line above gives. The pre-pass emits this rate as [`Provenance::FittedHere`] or not at
//!   all, so there is *"no rung below it for this parameter"*; this module is that rung, and a
//!   cell that reaches it takes [`DEFAULT_SSR_SUBSTITUTION_RATE`] and
//!   [`Provenance::Defaulted`].
//! - **The outlier weight**, which no fit produces anywhere: [`DEFAULT_OUTLIER_WEIGHT`],
//!   inherited from production at 0.01 and declared inherited (spec §4.5).
//!
//! # Why the outlier weight is not in any context's warrant
//!
//! **It is one run-wide constant and a warrant is per `(read group, candidate)`.** Folding a
//! constant that is defaulted everywhere into the per-cell warrant would make *every* repeat
//! tract's call [`Provenance::Defaulted`], and the distinction spec §4.4 says the warrant
//! exists to carry — a genotype resting on a fitted direction split against one resting on a
//! borrowed one — would be gone at every tract in every run. The same line is already drawn
//! one level down:
//! [`PART_REPEAT_SHARE_OF_WHOLE`](crate::ng::calling::likelihood::stutter_rates::PART_REPEAT_SHARE_OF_WHOLE)
//! is a placeholder inside every stutter model this module builds **from a fit** — the
//! defaulted ones come from [`StutterModel::hipstr_shipped`], whose part-repeat shares are its
//! own literals — and no provenance mentions it either. What the constant is owed instead is a line in the run's output, which is where
//! `spec/read_likelihoods.md` §3.6 already puts the contamination fraction.
//!
//! # The contaminant seed, and why it is frozen where the SNP/indel path's is not
//!
//! **Spec §4.5.1's third term is built here**, on every run whose fit found a fraction: how
//! common each of this tract's reachable lengths is in the contaminating population. It is the
//! genotype prior's own belief about which lengths this tract can be — the stratum's fitted
//! **length spectrum** — converted from the per-candidate shape the prior builds into a
//! distribution over lengths, which is what `c · seed(o)` asks for.
//!
//! **The cohort's own fitted frequencies at the locus would be the natural source and are
//! refused.** They are specific to the locus, which is the first thing the term needs; they fail
//! the second, because contamination is frozen before calling and they are what the caller
//! rewrites at every pass. The fitted spectrum meets both — it is indexed from *this tract's*
//! reference length, and it does not move while the loop iterates. **So a repeat tract's whole
//! row is frequency-free even under contamination**, which is the opposite of the SNP/indel
//! path, where `q(o)` is the loop's own estimate and moves with it (spec §3.6).
//!
//! # What this module does not build
//!
//! **Nor does anything here score a tract.** What this module assembles is what a tract's row
//! is scored *from*; the walk that hands each sample's reads to that row belongs to the calling
//! loop's driver (`inference::summarise_condition`).

use crate::ng::alignment::StutterModel;
use crate::ng::calling::genotype_prior::fill_seed_share_per_candidate;
use crate::ng::calling::likelihood::ssr::{
    DEFAULT_OUTLIER_WEIGHT, SsrContaminationMixture, SsrLocusParameters, SsrScoringContextTable,
};
use crate::ng::calling::likelihood::ssr_emission::{
    SsrCandidate, SsrScoringContext, fill_reachable_lengths,
};
use crate::ng::calling::likelihood::stutter_rates::stutter_model_for;
use crate::ng::calling::{CandidateAlleles, ContaminationView, FrozenParameters};
use crate::ng::locus_generation::LocusKind;
use crate::ng::parameter_estimation::Provenance;
use crate::ng::parameter_estimation::joint::stratum_fits::{
    FittedSlippage, LengthSpectrum, NoSlippage,
};
use crate::ng::parameter_estimation::ssr::RepeatCount;
use crate::ng::types::{ErrorRate, Motif, ReadGroupId};

use std::num::NonZeroU32;

/// The per-base substitution rate a tract is scored under where the fit has no entry for it —
/// **0.001**.
///
/// **Soft, and the only rung below the fit for this parameter.** The pre-pass emits the rate as
/// [`Provenance::FittedHere`] or not at all
/// (`parameter_estimation::ssr::substitution_rate_of`), so a `(read group, candidate stratum)`
/// the fit never accumulated — ordinarily a candidate several repeats from its tract's length,
/// whose stratum is not in the fit — has nowhere to fall. This is that floor, and a cell that
/// reaches it is marked [`Provenance::Defaulted`].
///
/// **It is *defined as* the SNP/indel path's default
/// ([`DEFAULT_ERROR_RATE`](crate::ng::parameter_estimation::generic::DEFAULT_ERROR_RATE)), so
/// editing that constant moves this one**, and it is still not the same parameter.
/// `doc/devel/ng/spec/read_likelihoods.md` §4.3 forbids tying the two *fitted* rates — each
/// absorbs what its own model cannot otherwise explain — and nothing here ties them: wherever
/// either is measured, its own measurement is used. One definition is deliberate at the point
/// where **neither** was measured, so that a run cannot end up defaulting its two error
/// parameters to two different guesses.
///
/// **The argument for the number is thinner than it looks, and that is worth saying.**
/// `parameter_prepass_ssr.md` §4.5 requires the two rates to agree to within a quarter-Phred
/// *where a stratum barely slips* — which is a statement about low-slippage strata that were
/// measured, not about strata that were not measured at all, which is the condition this is
/// reached under. And base quality inside tracts is systematically worse than outside them
/// (§4.1), so 0.001 is very likely optimistic at a tract. Nothing here measures by how much.
pub const DEFAULT_SSR_SUBSTITUTION_RATE: f64 =
    crate::ng::parameter_estimation::generic::DEFAULT_ERROR_RATE;

/// **The fitted numbers one repeat tract's scoring contexts are built from**, one pair per
/// `(read group, candidate)`, together with the length support the tract's junk term is spread
/// over.
///
/// # Why this is a type and not a `Vec<SsrScoringContext>`
///
/// **A context borrows its stutter model**, so the models and the contexts cannot live in one
/// struct — a struct holding both would have to refer to itself. This owns the models; the
/// contexts are built from a borrow of it by [`Self::scoring_contexts`] and live as long as
/// that borrow. It is the shape the row's own signature forces, not a choice made here.
///
/// # What is reused across loci and what is not
///
/// Everything this owns is a buffer that [`Self::gather_for_locus`] clears and refills, so one
/// of these per worker allocates on the first few tracts and then stops. **The contexts are the
/// exception**: they borrow this, so they cannot outlive one locus and
/// [`Self::scoring_contexts`] allocates a vector per tract. `#![forbid(unsafe_code)]` is what
/// closes the usual escape — a buffer re-lent under a shorter lifetime — and the cost is one
/// allocation of `read groups × candidates` contexts per repeat tract.
///
/// # Its cost is on the run's read-group axis, and that is worth knowing before a large run
///
/// The table covers **every read group of the run**, not the ones whose reads reached this
/// tract, because [`SsrScoringContextTable::of`] indexes by [`ReadGroupId`] directly and
/// because the contamination half of the mixture is checked against exactly this width
/// (`likelihood::ssr`'s row). So a tract costs `read groups × candidates` stratum lookups and
/// that many stutter models, whatever covered it: 6 lookups at one library and six candidates,
/// 6,000 at a thousand libraries. Nothing in this repository's benchmark cohorts reaches the
/// second — the tomato panel is 63 accessions of one library each — and no measurement here
/// says what it costs against the row's own work at that size.
#[derive(Clone, Debug, Default)]
pub struct TractScoringFits {
    /// The tract's repeat unit — **`None` until a tract has been gathered**, which is what
    /// every accessor here refuses on.
    ///
    /// **Stored rather than taken again**, and that is a correction rather than a
    /// convenience: while [`Self::scoring_contexts`] took a second motif, a caller could gather
    /// a mononucleotide tract and score it as a dinucleotide, and the context would then report
    /// an unreachable mass computed under one period beside a stutter model looked up under
    /// another — measured at 1.2 in a hundred million million against 2.0 in a hundred, a
    /// factor of 1.7 × 10¹², with no panic. The same argument
    /// [`SsrScoringContext::new`] makes for taking the mass from the distribution rather than
    /// from its caller: a fact taken twice is a fact that can disagree.
    motif: Option<Motif>,
    /// Whether the run's parameter fit found contamination — **what decides whether this tract
    /// gets the three-term form or the two**, and what
    /// [`Self::contaminant_length_frequencies`] was built, or left empty, from. Its one other
    /// reader checks that the fractions the row is handed came from the same run.
    run_fitted_contamination: bool,
    /// `read groups × candidates`, read-group-major — the same order
    /// [`SsrScoringContextTable`] indexes.
    stutter: Vec<StutterModel>,
    /// Parallel to [`Self::stutter`]: the per-base substitution rate for that cell.
    substitution_rate: Vec<ErrorRate>,
    /// Parallel to [`Self::stutter`]: the weaker of that cell's two warrants, which is what the
    /// context carries and what the locus's own warrant is folded from.
    warrant: Vec<Provenance>,
    /// The tract lengths the junk term is spread over — ascending, without repeats, and a
    /// property of the candidate set and the two slip cutoffs with no cohort in it (spec §4.5).
    reachable_lengths: Vec<u32>,
    /// **How common each of those lengths is in the contaminating population** — parallel to
    /// [`Self::reachable_lengths`], entry for entry, summing to one. **Empty on a run whose fit
    /// found no contamination**, where there is no third term to spread.
    ///
    /// It is the prior's own belief about which lengths this tract can be, converted from the
    /// per-candidate shape the prior builds into a distribution over lengths
    /// (`doc/devel/ng/spec/read_likelihoods.md` §4.5.1). Most entries are zero: the reachable
    /// support is every length the slip cutoffs admit from any candidate, and only the
    /// candidates themselves carry mass.
    contaminant_length_frequencies: Vec<f64>,
    /// How many candidates each read group's row covers — the stride, held so that the one
    /// spelling of it is this type's.
    candidates: usize,
    /// How many read groups the table covers — **held rather than divided out of the cell
    /// count**, so that an ungathered value asks for zero rows instead of dividing by zero.
    read_groups: usize,
    /// How many cells took [`StutterModel::hipstr_shipped`] because the slippage fit had no
    /// numbers for them.
    slippage_defaulted: usize,
    /// Of those, how many were defaulted because **the fit does not describe this run's read
    /// groups** — a library the pre-pass never saw, or a slippage map naming more groups than
    /// the fit was run over.
    ///
    /// **Separated because the other two absences are ordinary and these two are not.**
    /// [`NoSlippage`](crate::ng::parameter_estimation::joint::stratum_fits::NoSlippage) calls
    /// them *"the run is not what it claims"*: a candidate several repeats from its tract's
    /// length lands in no stratum on perfectly good data, but a read group the fit never named
    /// means the parameters and the reads came from different runs. Folding both into one count
    /// would let the second arrive as routine.
    slippage_defaulted_by_an_unknown_read_group: usize,
    /// How many cells took [`DEFAULT_SSR_SUBSTITUTION_RATE`] because the fit has no entry for
    /// them.
    substitution_defaulted: usize,
    /// The candidates' repeat counts as plain numbers, for the prior's seed builder — a type
    /// conversion held here rather than allocated per tract.
    ///
    /// **Empty on an uncontaminated run**, where the prior's shape is still read — every tract
    /// is seeded from it — but not through this type: what is not needed there is the
    /// *normalised per-candidate share* below, which only the contaminant term takes.
    candidate_repeat_counts: Vec<u32>,
    /// The prior's share for each candidate, summing to one, before it is scattered onto the
    /// length support. Empty on an uncontaminated run, which has no third term to spread.
    share_of_each_candidate: Vec<f64>,
}

impl TractScoringFits {
    /// Read the run's fitted parameters for one tract's candidates, clearing whatever the last
    /// tract left.
    ///
    /// `candidates` is in **genotype-table allele order**, which is the order the row indexes
    /// them in; [`tract_candidates`] is the one builder that pairs them with the bases the
    /// candidate table holds.
    ///
    /// # Panics
    ///
    /// If `candidates` is empty. A locus is called over at least its reference allele, and an
    /// empty candidate set would make the table's stride zero — which
    /// [`SsrScoringContextTable::new`] refuses one step later, naming the table rather than the
    /// locus.
    ///
    /// If the run has no read groups, which [`FrozenParameters`] already refuses at
    /// construction; restated here because the table's shape is the product of the two and a
    /// zero on either axis is an empty table that would fail at whichever tract first held a
    /// read.
    pub fn gather_for_locus(
        &mut self,
        motif: &Motif,
        candidates: &[SsrCandidate<'_>],
        tract_prior: TractPrior<'_>,
        parameters: &FrozenParameters<'_>,
    ) {
        let reference_repeats = tract_prior.reference_repeats;
        assert!(
            !candidates.is_empty(),
            "a repeat tract is called over at least its reference allele, so a locus with no \
             candidates is a candidate set that went missing on the way in"
        );
        // **The reference tract's repeat count is a second spelling of `candidates[0]`, and
        // this is what stops the two coming apart.** The reference allele is id 0 of every
        // candidate table and these candidates are in that order, so the caller has the number
        // already; it is taken as an argument because *which* repeat count keys the prior's
        // length spectrum is the one thing at a tract that is easy to get wrong, and naming it
        // makes the wrong one something somebody has to type. A candidate's count passed here
        // re-centres the spectrum on that candidate and flattens the prior, and the seed still
        // sums to one — which is why `fill_seed_share_per_candidate`'s own documentation lists
        // it among the mistakes nothing inside it can catch. Here it is catchable.
        assert_eq!(
            reference_repeats.get(),
            candidates[0].repeat_count.get(),
            "this tract's reference allele holds {} whole repeats and the prior was asked to \
             centre on {}: a candidate's count passed as the reference tract's re-centres the \
             fitted length spectrum on that candidate",
            candidates[0].repeat_count.get(),
            reference_repeats.get()
        );
        let read_groups = parameters.read_group_count();
        // **Held in debug only, and no test can reach it**: `FrozenParameters` refuses an empty
        // calibration list at construction, and that list is the axis this counts, so a run
        // that reached here has at least one read group. What it guards is that refusal being
        // relaxed later, which would otherwise surface as an empty context table at whichever
        // tract first held a read.
        debug_assert!(
            read_groups > 0,
            "every read of the run belongs to a read group and a run has at least one, so a \
             tract scored against no read group is a run whose read-group axis went missing"
        );

        self.stutter.clear();
        self.substitution_rate.clear();
        self.warrant.clear();
        self.motif = Some(*motif);
        self.run_fitted_contamination = !parameters.contamination_is_absent();
        self.candidates = candidates.len();
        self.read_groups = read_groups;
        self.slippage_defaulted = 0;
        self.slippage_defaulted_by_an_unknown_read_group = 0;
        self.substitution_defaulted = 0;

        let period = motif.ssr_period();
        for group in 0..read_groups {
            // `ReadGroupId` is the run's own index over the same axis `read_group_count`
            // counts, so this is the identity rather than a mapping. Checked rather than cast:
            // `group as u32` would silently make read group 2³² into read group 0 and score a
            // library against another's polymerase.
            let read_group = ReadGroupId(
                u32::try_from(group).expect("a run has fewer read groups than a u32 can name"),
            );
            for candidate in candidates {
                let repeats = candidate.repeat_count.get();
                let slippage =
                    parameters
                        .ssr_slippage_fits()
                        .at(read_group, period.get(), u64::from(repeats));
                let (stutter, slippage_warrant) = match &slippage {
                    Ok(fitted) => (stutter_model_for(&fitted.slippage), warrant_of(fitted)),
                    Err(absence) => {
                        self.slippage_defaulted += 1;
                        if matches!(
                            absence,
                            NoSlippage::UnknownReadGroup | NoSlippage::GroupNotInTheFit { .. }
                        ) {
                            self.slippage_defaulted_by_an_unknown_read_group += 1;
                        }
                        (StutterModel::hipstr_shipped(), Provenance::Defaulted)
                    }
                };
                let rate =
                    parameters.ssr_substitution_rate_at(read_group, period, RepeatCount(repeats));
                let (substitution_rate, substitution_warrant) = match rate {
                    Some(fitted) => (fitted.value, fitted.provenance),
                    None => {
                        self.substitution_defaulted += 1;
                        (default_substitution_rate(), Provenance::Defaulted)
                    }
                };
                self.stutter.push(stutter);
                self.substitution_rate.push(substitution_rate);
                self.warrant
                    .push(slippage_warrant.weaker_of(substitution_warrant));
            }
        }

        fill_reachable_lengths(candidates, motif, &mut self.reachable_lengths);
        self.fill_contaminant_length_frequencies(candidates, tract_prior, parameters);
    }

    /// **Turn the prior's belief about this tract's lengths into the third term of the row's
    /// mixture** — how common each reachable length is in the contaminating population.
    ///
    /// Left empty on a run whose fit found no contamination, where there is no third term.
    ///
    /// # The two halves speak different units, and this is where they meet
    ///
    /// **The prior speaks in whole repeats.** Its length spectrum is indexed by offset from the
    /// reference tract's repeat count, so what it hands back is one share per *candidate*, each
    /// candidate placed by the repeat count the locus generator measured for it.
    ///
    /// **The reads speak in bases.** An observation shows a byte length, and the reachable
    /// support the outlier term is spread over is a list of byte lengths — so this term has to
    /// be one too, or the row's three terms would not be probabilities of the same event.
    ///
    /// So each candidate's share is added at the support entry its **bases** land on. **Two
    /// candidates of one byte length therefore share one entry**, and they sum into it rather
    /// than each taking the full share — which is what
    /// [`SsrContaminationMixture::contaminant_length_frequencies`]' own documentation asks for.
    ///
    /// **An interrupted tract is exactly that case, and what it costs is worth stating
    /// precisely.** Such a candidate can spell as many bases as a clean one while holding fewer
    /// whole repeats, so the prior places the two at different offsets and gives them different
    /// shares. The *read likelihood* separates them easily — a read carrying the interruption
    /// scores higher against the interrupted allele by about 28 Phred per distinguishing base at
    /// an error rate of 1 in 200 (spec §4.6). It is this term that cannot: a contaminating read
    /// shows a length, and two spellings of one length are one length. **How the prior should
    /// divide a length class between two such candidates is a separate question and is not
    /// answered here** — it is `spec/calling_priors.md`'s open question 3, stated in its §5.2,
    /// and keying this to lengths is what keeps it in one place.
    ///
    /// # The shares are renormalised over the candidates, and that raises the term
    ///
    /// The prior's shape is divided by its total **over this locus's candidates only**, so the
    /// spectrum's mass at lengths no candidate carries is spread onto the ones that do rather
    /// than dropped. The row's own contract forces it — a distribution that did not sum to one
    /// would make the three terms incomparable — and it is what
    /// [`fill_seed_share_per_candidate`] is documented to return.
    ///
    /// **What it costs, with its size.** On this module's own two-candidate fixture — reference
    /// 6 repeats, candidates at 6 and 7, the fit's spectrum putting 0.44 of a stratum's
    /// chromosomes at the reference length — the seed says 0.8 rather than 0.44. At a fitted
    /// fraction of 5 in 100 the mixture then credits 0.040 of a read at that length to the
    /// contaminant where the spectrum alone would credit 0.022: **1.8 times the fitted weight**.
    /// The effective fraction at a candidate length is `c` divided by the share of the
    /// spectrum's mass the candidate set covers, and a locus whose candidates cover little of
    /// what the stratum spreads over inflates it most.
    ///
    /// **What that is not.** Spec §4.5.1 says a contaminating read at a length no candidate
    /// carries "falls to the outlier floor instead — which is where they go today, so nothing
    /// is lost", and that stays true: such a length gets no entry here. What moves is the
    /// weight *between* the candidates' own lengths, and it moves toward them.
    ///
    /// # Why the shares are the prior's rather than the loop's own frequencies
    ///
    /// The cohort's fitted frequencies at this locus are the natural answer to *how common is
    /// each length here*, and they are refused: contamination is frozen before calling, and
    /// those are what the loop rewrites at every pass (spec §4.5.1). The fitted length spectrum
    /// is specific to this locus — it is indexed from this tract's own reference length — and it
    /// does not move while the loop iterates. **So a repeat tract's whole row is frequency-free
    /// even under contamination**, which is the opposite of the SNP/indel path, where `q(o)` is
    /// the loop's own estimate and moves with it (§3.6).
    ///
    /// # Panics
    ///
    /// If a candidate's own byte length is not in the reachable support. It always is — the
    /// support is built from each candidate's length plus every slip the cutoffs admit, and a
    /// slip of nothing is among them — so this is a check on that construction rather than on
    /// the caller.
    fn fill_contaminant_length_frequencies(
        &mut self,
        candidates: &[SsrCandidate<'_>],
        tract_prior: TractPrior<'_>,
        parameters: &FrozenParameters<'_>,
    ) {
        // **Cleared before the early return, so an uncontaminated tract cannot carry a
        // contaminated one's seed.** Only the first of the three is load-bearing in a run — a
        // stale `contaminant_length_frequencies` of the same width would survive the `resize`
        // below and the `+=` would accumulate onto it — because a run is contaminated or not
        // as a whole and cannot change between two of its loci. The other two are cleared
        // here because this type is public and default-constructible, and because the two
        // fields say they are empty on an uncontaminated run.
        self.contaminant_length_frequencies.clear();
        self.candidate_repeat_counts.clear();
        self.share_of_each_candidate.clear();
        if parameters.contamination_is_absent() {
            return;
        }

        // The prior's shape, one share per candidate, summing to one.
        self.candidate_repeat_counts.extend(
            candidates
                .iter()
                .map(|candidate| candidate.repeat_count.get()),
        );
        self.share_of_each_candidate
            .resize(candidates.len(), f64::NAN);
        fill_seed_share_per_candidate(
            &self.candidate_repeat_counts,
            tract_prior.reference_repeats.get(),
            tract_prior.length_spectrum,
            &mut self.share_of_each_candidate,
        );

        self.contaminant_length_frequencies
            .resize(self.reachable_lengths.len(), 0.0);
        for (candidate, share) in candidates.iter().zip(&self.share_of_each_candidate) {
            let spelled = candidate.bases.len() as u32;
            let at = self
                .reachable_lengths
                .binary_search(&spelled)
                .unwrap_or_else(|_| {
                    panic!(
                        "a candidate spelling {spelled} bases is not among the {} lengths this \
                         tract can reach, and every candidate's own length is reachable by \
                         construction — the support and the candidates were built from \
                         different sets",
                        self.reachable_lengths.len()
                    )
                });
            self.contaminant_length_frequencies[at] += share;
        }
    }

    /// The scoring contexts, in the order [`SsrScoringContextTable`] indexes them.
    ///
    /// **The candidates must be the ones [`Self::gather_for_locus`] read**, because the
    /// unreachable mass each context carries is computed from the candidate's own repeat count
    /// against the model looked up for it. Handing a different set would pair one candidate's
    /// mass with another's model — a plausible number, and no panic. **The motif is not asked
    /// for at all**, for the same reason one field over: it is the one this was gathered under.
    ///
    /// # Panics
    ///
    /// If nothing has been gathered. A default-constructed value has no motif and no cells, and
    /// the contexts it would return are an empty table that the row refuses several frames
    /// later.
    ///
    /// If `candidates` is not the length this was gathered for. That is the mispairing above,
    /// caught where the two are still distinguishable rather than at whichever genotype first
    /// read the wrong cell.
    #[must_use]
    pub fn scoring_contexts<'a>(
        &'a self,
        candidates: &[SsrCandidate<'_>],
    ) -> Vec<SsrScoringContext<'a>> {
        let motif = self.gathered_motif();
        assert_eq!(
            candidates.len(),
            self.candidates,
            "these fits were gathered for {} candidates and {} were handed to the contexts, so \
             one of the two belongs to a different locus",
            self.candidates,
            candidates.len()
        );
        let mut contexts = Vec::with_capacity(self.stutter.len());
        // **Read-group-major, candidate within it** — the order
        // [`SsrScoringContextTable::of`] indexes, and the order `gather_for_locus` filled.
        // Written as the same two nested loops rather than as one walk over the cells, so that
        // the two orders are the same statement in both places.
        for group in 0..self.read_groups {
            for (candidate, allele) in candidates.iter().enumerate() {
                let cell = group * self.candidates + candidate;
                contexts.push(SsrScoringContext::new(
                    motif,
                    &self.stutter[cell],
                    allele,
                    self.substitution_rate[cell],
                    [self.warrant[cell]],
                ));
            }
        }
        contexts
    }

    /// **Everything the row takes at this tract** — the candidates, their scoring contexts, the
    /// outlier weight, the lengths it is spread over, and, on a run whose fit found
    /// contamination, the third term of the mixture.
    ///
    /// **The third term is present exactly when the run's fit found a fraction**, which is spec
    /// §4.5.1's rule and the same one the SNP/indel path follows: contamination is a property of
    /// the sample rather than of the marker, so a caller that corrects for it at one kind of
    /// locus and not the other is treating one number as two.
    ///
    /// `contamination_of_each_read_group` is the run's own list, in read-group order — the
    /// fractions the pre-pass fitted, each carrying whose reads it was fitted from. It must be
    /// empty exactly when the fit found nothing, which is what
    /// [`FrozenParameters::contamination_is_absent`] answers.
    ///
    /// # Panics
    ///
    /// If `contamination_of_each_read_group` disagrees with what this was gathered under about
    /// whether the run is contaminated. The seed was built, or not built, from that same
    /// predicate one call earlier, so a disagreement means the two came from different runs —
    /// and the row would then be handed a fraction against an empty seed, or a seed nothing
    /// scales.
    ///
    /// If `candidates` or `contexts` is not the shape this was gathered for — the same
    /// mispairing [`Self::scoring_contexts`] refuses, restated because a caller may hold
    /// contexts from an earlier tract. **The second check catches a truncated slice rather than
    /// another tract's contexts**: two tracts of one run share the read-group count, so their
    /// cell counts differ only where their candidate counts do, which the first check already
    /// covers.
    #[must_use]
    pub fn locus_parameters<'a>(
        &'a self,
        candidates: &'a [SsrCandidate<'a>],
        contexts: &'a [SsrScoringContext<'a>],
        contamination_of_each_read_group: &'a [ContaminationView],
    ) -> SsrLocusParameters<'a> {
        assert_eq!(
            contamination_of_each_read_group.is_empty(),
            !self.run_fitted_contamination,
            "these fits were gathered for a run whose parameter fit {} contamination, and a \
             contamination list of {} read groups reached the row — so the seed and the \
             fractions came from different runs",
            if self.run_fitted_contamination {
                "found"
            } else {
                "found no"
            },
            contamination_of_each_read_group.len()
        );
        // **And it must be one fraction per read group of the run**, which is the axis the
        // context table below is indexed on. A shorter list is not caught by the emptiness
        // check above and would surface only if a read from a high-numbered library happened to
        // arrive at this tract — at which point the row refuses it several frames away, naming
        // the table rather than the locus.
        assert!(
            contamination_of_each_read_group.is_empty()
                || contamination_of_each_read_group.len() == self.read_groups,
            "these fits cover {} read groups and {} contamination fractions reached the row, so \
             one of the two describes a different run",
            self.read_groups,
            contamination_of_each_read_group.len()
        );
        assert_eq!(
            candidates.len(),
            self.candidates,
            "these fits were gathered for {} candidates and {} candidates reached the row",
            self.candidates,
            candidates.len()
        );
        assert_eq!(
            contexts.len(),
            self.stutter.len(),
            "these fits cover {} (read group, candidate) cells and {} contexts reached the row, \
             so the contexts are a truncated slice or were built from a different tract",
            self.stutter.len(),
            contexts.len()
        );
        SsrLocusParameters {
            candidates,
            contexts: SsrScoringContextTable::new(contexts, self.candidates),
            outlier_weight: DEFAULT_OUTLIER_WEIGHT,
            reachable_lengths: &self.reachable_lengths,
            contamination: self
                .run_fitted_contamination
                .then(|| SsrContaminationMixture {
                    fraction_of_each_read_group: contamination_of_each_read_group,
                    contaminant_length_frequencies: &self.contaminant_length_frequencies,
                }),
        }
    }

    /// The tract lengths the junk term is spread over.
    #[inline]
    #[must_use]
    pub fn reachable_lengths(&self) -> &[u32] {
        &self.reachable_lengths
    }

    /// How many `(read group, candidate)` cells this tract was gathered over.
    #[inline]
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.stutter.len()
    }

    /// How many cells found no numbers in the slippage fit and took the stated constant.
    ///
    /// **Ordinary rather than an error**, and worth reporting for the same reason the
    /// contamination fraction is: a tract scored mostly on defaulted slippage and one scored on
    /// a fit are different claims about the same genotype.
    #[inline]
    #[must_use]
    pub fn cells_with_no_fitted_slippage(&self) -> usize {
        self.slippage_defaulted
    }

    /// Of those, how many were defaulted because the fit does not describe this run's read
    /// groups — **the half of the absence that is not ordinary**, and the number a run should
    /// act on rather than record.
    #[inline]
    #[must_use]
    pub fn cells_whose_read_group_the_fit_does_not_describe(&self) -> usize {
        self.slippage_defaulted_by_an_unknown_read_group
    }

    /// How many cells found no fitted substitution rate and took
    /// [`DEFAULT_SSR_SUBSTITUTION_RATE`].
    #[inline]
    #[must_use]
    pub fn cells_with_no_fitted_substitution_rate(&self) -> usize {
        self.substitution_defaulted
    }

    /// **The weakest warrant behind any parameter that reached this tract** — what the locus's
    /// record is entitled to claim (spec §4.4).
    ///
    /// **Combined, never branched on**: a call resting on one fitted number and one borrowed
    /// one is a borrowed call.
    ///
    /// # Panics
    ///
    /// If nothing has been gathered. **The fold's identity is the strongest rung on the
    /// ladder**, so an ungathered value would answer [`Provenance::FittedHere`] over zero cells
    /// — the best warrant a call can have, from having read nothing. A default-constructed
    /// value is reachable from outside this module, so this is held rather than assumed.
    #[must_use]
    pub fn weakest_warrant(&self) -> Provenance {
        assert!(
            !self.warrant.is_empty(),
            "no tract has been gathered, so there is no warrant to report — and the fold's \
             identity would answer `FittedHere` over no parameters at all"
        );
        self.warrant
            .iter()
            .copied()
            .fold(Provenance::FittedHere, Provenance::weaker_of)
    }

    /// The motif this was gathered under.
    ///
    /// # Panics
    ///
    /// If nothing has been gathered.
    fn gathered_motif(&self) -> &Motif {
        self.motif.as_ref().unwrap_or_else(|| {
            panic!(
                "no tract has been gathered, so there is no motif to score one under — \
                 `gather_for_locus` is what fills these fits"
            )
        })
    }
}

/// **What a repeat tract's genotype prior believes about its lengths** — the two things every
/// consumer of that belief needs, looked up **once** per tract and handed to both.
///
/// # Why one value rather than two arguments
///
/// The prior's belief is used twice at each tract: the genotype prior's own seed, and the third
/// term of the read-likelihood mixture, which asks how common each length is in the
/// contaminating population. Both read the same fitted **length spectrum**, and both must read
/// the *same* one — the run reports which rung of the tract ladder answered, and a run that
/// reported one rung while scoring against another would be saying something false about its own
/// calls.
///
/// **Two lookups keyed identically is a coincidence somebody can break**; one lookup passed to
/// both is not. So the caller looks it up and this carries the pair.
#[derive(Debug, Clone, Copy)]
pub struct TractPrior<'a> {
    /// **The tract's own reference repeat count**, which is what the spectrum's offsets are
    /// measured from — not any candidate's. It is entry 0 of the candidate table, and
    /// [`TractScoringFits::gather_for_locus`] checks that it is.
    pub reference_repeats: RepeatCount,
    /// How that stratum's chromosomes are spread over tract lengths, and how strongly the fit
    /// holds them — with the rung of the tract ladder it came from on it.
    pub length_spectrum: LengthSpectrum<'a>,
}

/// **What one cell's slippage numbers are entitled to claim**, from where the fit says they
/// came.
///
/// The fit reports three sources per number — the cell's own fit, its period's curve, or a
/// blend of the two — and each carries the share the curve took
/// ([`LevelSource::curve_weight`](crate::ng::parameter_estimation::joint::slippage_curve::LevelSource::curve_weight)).
/// **A number the curve entered at all is [`Provenance::Borrowed`]**: the curve is drawn
/// through other strata, so a cell reading off it is reading a neighbouring grain's
/// measurement, which is exactly what that rung means. Only a cell whose level and both shares
/// are wholly its own is [`Provenance::FittedHere`].
///
/// **The threshold is zero rather than a share**, deliberately: any cut-off between "mostly its
/// own" and "mostly the curve's" would be a number invented here, and the share itself is not
/// lost — it stays on the `LevelProvenance` the fit emitted, for a consumer that wants to weigh
/// it (`str_slippage_level_curve.md` §8).
///
/// **The two shares are asked separately**, because the fit reports them separately: a cell can
/// fit its own contraction bias and read its fall-off off a curve, and a warrant that looked at
/// one of them would call such a cell fitted.
///
/// **A cell with slippage numbers and no shares provenance is [`Provenance::Borrowed`].** Both
/// of the fit's own paths set the slippage numbers and the shares provenance from one mask, so
/// this is a `StratumFits` assembled by hand rather than anything a run produces; the numbers
/// are still the fit's, so `Defaulted` would be wrong, and the cell cannot claim its shares
/// were fitted here.
fn warrant_of(fitted: &FittedSlippage) -> Provenance {
    let level_is_its_own = fitted.level.source.curve_weight() == 0.0;
    let shares_are_its_own = fitted.shares.is_some_and(|shares| {
        shares.shorter_share.source.curve_weight() == 0.0
            && shares.fall_off.source.curve_weight() == 0.0
    });
    if level_is_its_own && shares_are_its_own {
        Provenance::FittedHere
    } else {
        Provenance::Borrowed
    }
}

/// [`DEFAULT_SSR_SUBSTITUTION_RATE`] as the checked type.
///
/// The constant is a probability by inspection, so the conversion cannot fail; done here once
/// rather than at the one call site so that a future edit to the constant fails in one place
/// with a message naming it.
fn default_substitution_rate() -> ErrorRate {
    ErrorRate::try_new(DEFAULT_SSR_SUBSTITUTION_RATE)
        .expect("the default STR substitution rate is a probability")
}

/// **The candidate table as the emission sees it** — the bases the locus is called over, paired
/// with how many repeats each candidate's tract holds.
///
/// # Why the repeat counts arrive rather than being counted here
///
/// **A tract's repeat count is not its byte length divided by the period.** An interrupted
/// tract holds fewer whole repeats than its bases would suggest, and the count is what keys the
/// stratum lookup — so counting it here would both re-measure something the locus generator has
/// already measured and get the interrupted case wrong (spec §7, and
/// [`SsrCandidate::repeat_count`]'s own contract).
///
/// **Today they come from a fixture**, because the repeat-tract half of candidate selection is
/// unwritten (`candidate_alleles_ssr.md`). When it lands, they come from it, and this signature
/// does not change.
///
/// # Panics
///
/// If the two lists are different lengths, or if `alleles` is not a single repeat tract's
/// table. **The two ways the second can fail are different failures and say so.** A SNP/indel
/// table reaching here is a locus routed to the wrong read model, which would otherwise surface
/// as a stutter score over bases that are not a tract. A repeat **bundle** is not misrouted —
/// the calling seam sends a bundle down the repeat path deliberately, which is how every other
/// consumer of [`LocusKind`] groups the two — but nothing here scores one, and no producer
/// emits a bundle into calling today.
#[must_use]
pub fn tract_candidates<'a>(
    alleles: &'a CandidateAlleles,
    repeat_counts: &[NonZeroU32],
) -> Vec<SsrCandidate<'a>> {
    assert!(
        !matches!(alleles.kind(), LocusKind::SsrBundle),
        "the locus at this table is a bundle of repeat tracts, and nothing scores a bundle \
         yet: this path is written for one tract, whose candidates are its own lengths. The \
         calling seam sends a bundle down the repeat path on purpose, so this is a gap in \
         what the repeat path covers rather than a locus routed to the wrong read model"
    );
    assert!(
        matches!(alleles.kind(), LocusKind::Ssr(_)),
        "these candidates were generated at a {:?} locus and are being scored as a repeat \
         tract, so the locus was routed to the wrong read model",
        alleles.kind()
    );
    assert_eq!(
        alleles.len(),
        repeat_counts.len(),
        "this tract is called over {} candidates and {} repeat counts were supplied, so one of \
         the two belongs to a different locus",
        alleles.len(),
        repeat_counts.len()
    );
    alleles
        .iter()
        .zip(repeat_counts)
        .map(|(bases, repeat_count)| SsrCandidate {
            bases,
            repeat_count: *repeat_count,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::ng::calling::ContaminationView;
    use crate::ng::calling::GenotypeTable;
    use crate::ng::calling::ReadGroupCalibration;
    use crate::ng::calling::SsrSampleEvidence;
    use crate::ng::calling::genotype_prior::{SeedRegime, SpectrumSeed};
    use crate::ng::calling::likelihood::SsrRowScratch;
    use crate::ng::calling::likelihood::ssr::genotype_log_likelihood_row;
    use crate::ng::calling::likelihood::ssr_emission::{
        StutterSubstitutionEmission, StutterSubstitutionScratch,
    };
    use crate::ng::locus_generation::{ReadWitness, SequenceObservation, SsrDetail};
    use crate::ng::parameter_estimation::Estimate;
    use crate::ng::parameter_estimation::joint::census::Stratum as FitStratum;
    use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;
    use crate::ng::parameter_estimation::joint::share_curve::ShareSource;
    use crate::ng::parameter_estimation::joint::slippage_curve::{
        LevelSource, RiseShape, SlippageCurve,
    };
    use crate::ng::parameter_estimation::joint::ssr_fit::{
        DerivedStratum, LevelProvenance, ShareProvenance, SharesProvenance, Slippage,
        StratumOutcome,
    };
    use crate::ng::parameter_estimation::joint::stratum_fits::StratumFits;
    use crate::ng::parameter_estimation::ssr::{Stratum as SsrStratum, StratumKey};
    use crate::ng::types::{InbreedingF, LogProb, Ploidy, SsrPeriod};

    use std::collections::BTreeMap;

    /// **A dinucleotide, so the reachable support is wider than the candidate set.** Period 1
    /// would make every one-base step a whole repeat and hide the part-repeat half of the
    /// distribution.
    const MOTIF: &[u8] = b"AC";

    /// The two candidates every fixture here is called over, in genotype-table allele order:
    /// **6 repeats and 11**, five apart, so they land in different strata and a lookup that
    /// hoisted itself out of the candidate loop returns a different number rather than the same
    /// one.
    const CANDIDATE_REPEATS: [u32; 2] = [6, 11];

    /// **Three read groups against two candidates**, deliberately unequal. At an equal shape a
    /// table filled read-group-major and one filled candidate-major are the same length and the
    /// same set of cells, so a transposition passes every shape check; at 3 × 2 the two orders
    /// disagree cell by cell.
    const READ_GROUPS: usize = 3;

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("a diploid")
    }

    /// **Borrowed rather than returned by value**, because a context borrows the motif it was
    /// built with and a temporary would be freed at the end of the statement that built it.
    fn motif() -> &'static Motif {
        static PARSED: std::sync::OnceLock<Motif> = std::sync::OnceLock::new();
        PARSED.get_or_init(|| Motif::new(MOTIF).expect("a valid dinucleotide motif"))
    }

    fn period() -> SsrPeriod {
        motif().ssr_period()
    }

    /// The reference tract's own repeat count — entry 0 of [`CANDIDATE_REPEATS`], because the
    /// reference allele is id 0 of every candidate table. It is what the prior's length
    /// spectrum measures its offsets from, and it is not any other candidate's.
    fn reference_repeats() -> RepeatCount {
        RepeatCount(CANDIDATE_REPEATS[0])
    }

    /// What the genotype prior believes about a tract of this reference length, looked up the
    /// way the driver looks it up.
    fn tract_prior<'a>(
        reference_repeats: RepeatCount,
        parameters: &FrozenParameters<'a>,
    ) -> TractPrior<'a> {
        TractPrior {
            reference_repeats,
            length_spectrum: parameters
                .ssr_length_spectrum_at(motif().ssr_period(), reference_repeats),
        }
    }

    fn repeats(count: u32) -> NonZeroU32 {
        NonZeroU32::new(count).expect("a candidate always holds a repeat")
    }

    /// The bases of a tract holding `count` whole copies of the motif.
    fn tract(count: u32) -> Vec<u8> {
        MOTIF.repeat(count as usize)
    }

    /// The two candidates' bases, owned, so the borrows below stay simple.
    fn candidate_bases() -> Vec<Vec<u8>> {
        CANDIDATE_REPEATS.iter().copied().map(tract).collect()
    }

    fn candidates(bases: &[Vec<u8>]) -> Vec<SsrCandidate<'_>> {
        bases
            .iter()
            .zip(CANDIDATE_REPEATS)
            .map(|(bases, count)| SsrCandidate {
                bases,
                repeat_count: repeats(count),
            })
            .collect()
    }

    /// **A slippage level that is different in every cell of the table**, and different for a
    /// reason on each axis: longer tracts slip more, and one library slips more than another.
    ///
    /// **Keyed by the candidate's repeat count and not by its position**, which is what a
    /// stratum lookup does — a fixture keyed by position would make permuting the candidates
    /// change the answer for a reason that has nothing to do with the lookup.
    fn level_at(slippage_group: u32, repeat_count: u32) -> f64 {
        0.02 + 0.004 * f64::from(repeat_count) + 0.003 * f64::from(slippage_group)
    }

    fn slippage_at(slippage_group: u32, repeat_count: u32) -> Slippage {
        Slippage {
            level: level_at(slippage_group, repeat_count),
            // Contraction-biased, as every real fit is, and not a half — so a dropped or
            // swapped direction share changes a number.
            shorter_share: 0.83,
            fall_off: 0.35,
        }
    }

    /// A level the cell fitted itself, carrying a distinguishable read count.
    fn fitted_level(slipped_reads: f64) -> LevelProvenance {
        LevelProvenance {
            source: LevelSource::Cell,
            curve: None,
            reach: None,
            slipped_reads: Some(slipped_reads),
        }
    }

    /// A level the cell read whole off its period's curve — the fit's numbers, and not the
    /// cell's own measurement.
    fn level_off_the_curve() -> LevelProvenance {
        LevelProvenance {
            source: LevelSource::Curve,
            curve: Some(SlippageCurve {
                rise_shape: RiseShape::new(1.0).expect("a rise shape of one is on the grid"),
                intercept: 0.01,
                slope: 0.01,
                fitted_from: 8,
                fitted_to: 12,
                held_out_error: 0.077,
                cells: 5,
            }),
            reach: None,
            slipped_reads: None,
        }
    }

    /// **The two shares are named separately and never defaulted to one value**, because the
    /// fit reports them separately: a cell can fit its own contraction bias and read its
    /// fall-off off a curve. A helper taking one source made those two conjuncts of
    /// `warrant_of` indistinguishable, and dropping either of them left all nineteen tests
    /// green.
    fn shares_from(shorter: ShareSource, fall_off: ShareSource) -> SharesProvenance {
        let of = |source| ShareProvenance {
            source,
            curve: None,
            reach: None,
        };
        SharesProvenance {
            slipped_reads: Some(400.0),
            shorter_share: of(shorter),
            fall_off: of(fall_off),
        }
    }

    /// One stratum's row at the fixture's own period, with a slippage number for every
    /// slippage group named, and both shares fitted by the cell.
    fn stratum_row(
        repeat_count: u32,
        groups: usize,
        level: impl Fn(u32) -> LevelProvenance,
    ) -> StratumOutcome {
        stratum_row_at(
            MOTIF.len() as u8,
            repeat_count,
            groups,
            level,
            ShareSource::Stratum,
            ShareSource::Stratum,
        )
    }

    /// The same, with the period and each share's source named — **the period because it is 2
    /// in every other fixture and 2 is also the candidate count**, so a lookup keyed by a
    /// hard-coded 2 was indistinguishable from one keyed by the motif until this existed.
    fn stratum_row_at(
        period: u8,
        repeat_count: u32,
        groups: usize,
        level: impl Fn(u32) -> LevelProvenance,
        shorter: ShareSource,
        fall_off: ShareSource,
    ) -> StratumOutcome {
        let slippage: Vec<Option<Slippage>> = (0..groups as u32)
            .map(|group| Some(slippage_at(group, repeat_count)))
            .collect();
        let level_provenance: Vec<Option<LevelProvenance>> =
            (0..groups as u32).map(|group| Some(level(group))).collect();
        let shares_provenance = vec![Some(shares_from(shorter, fall_off)); groups];
        StratumOutcome::Derived(Box::new(DerivedStratum {
            stratum: FitStratum {
                period,
                reference_repeats: u64::from(repeat_count),
            },
            slippage,
            level_provenance,
            shares_provenance,
            tracts_of_its_own: 4,
            reads_crossing: 40,
        }))
    }

    /// **Every read group in its own slippage group**, which is the specified grain — so the
    /// read-group axis of the lookup is live and a cell that read another group's numbers
    /// comes back with a different level.
    fn one_group_each(read_groups: usize) -> BTreeMap<ReadGroupId, u32> {
        (0..read_groups as u32)
            .map(|group| (ReadGroupId(group), group))
            .collect()
    }

    /// The whole fit: both candidates' strata, every group fitted on its own reads.
    fn fits_for_both_candidates() -> StratumFits {
        StratumFits::over(
            &[
                stratum_row(CANDIDATE_REPEATS[0], READ_GROUPS, |group| {
                    fitted_level(400.0 + f64::from(group))
                }),
                stratum_row(CANDIDATE_REPEATS[1], READ_GROUPS, |group| {
                    fitted_level(800.0 + f64::from(group))
                }),
            ],
            one_group_each(READ_GROUPS),
        )
    }

    /// The substitution rate for one cell — different on both axes, for the same reason the
    /// slippage level is.
    fn substitution_rate_at(read_group: u32, repeat_count: u32) -> f64 {
        1e-3 * (1.0 + 0.1 * f64::from(repeat_count) + f64::from(read_group))
    }

    /// A rate for every `(read group, candidate)` this fixture reaches.
    fn all_substitution_rates() -> BTreeMap<StratumKey, Estimate<ErrorRate>> {
        let mut rates = BTreeMap::new();
        for group in 0..READ_GROUPS as u32 {
            for repeat_count in CANDIDATE_REPEATS {
                rates.insert(
                    StratumKey {
                        read_group: ReadGroupId(group),
                        stratum: SsrStratum::new(period(), RepeatCount(repeat_count)),
                        ploidy: diploid(),
                    },
                    Estimate {
                        value: ErrorRate::try_new(substitution_rate_at(group, repeat_count))
                            .expect("a valid rate"),
                        provenance: Provenance::FittedHere,
                        observations: 40_000,
                    },
                );
            }
        }
        rates
    }

    fn a_contamination_view(fraction: f64) -> ContaminationView {
        ContaminationView {
            fraction,
            markers_with_reads: 1_000,
            reads_on_markers: 30_000,
            source: ContaminationSource::ThisReadGroupsReads,
        }
    }

    fn calibrations() -> Vec<ReadGroupCalibration> {
        vec![ReadGroupCalibration::defaulted(); READ_GROUPS]
    }

    fn outbred(samples: usize) -> Vec<InbreedingF> {
        vec![InbreedingF::try_new(0.0).expect("an outbred sample"); samples]
    }

    fn run<'a>(
        calibration: &'a [ReadGroupCalibration],
        inbreeding: &'a [InbreedingF],
        fits: &'a StratumFits,
        rates: &'a BTreeMap<StratumKey, Estimate<ErrorRate>>,
    ) -> FrozenParameters<'a> {
        FrozenParameters::uncontaminated(
            calibration,
            inbreeding,
            SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
            fits,
            rates,
            diploid(),
        )
    }

    /// **Every cell of the table carries its own read group's and its own candidate's fitted
    /// numbers** — six cells, six different slippage levels and six different substitution
    /// rates, none of them equal to another.
    ///
    /// **This is the test the table's shape exists for**, and both ways of getting the shape
    /// wrong were run rather than reasoned about. Reading the cells candidate-major
    /// (`candidate * read_groups + group`) while `gather_for_locus` filled them
    /// read-group-major leaves all six models present and fails at read group 0 / candidate 1,
    /// which then carries **read group 1's slippage level, 0.067, where its own is 0.064**.
    /// Hoisting the lookup out of the candidate loop — one lookup per read group, at candidate
    /// 0's repeat count — fails at the same cell with **0.044 against 0.064**, candidate 0's
    /// stratum standing in for candidate 1's.
    #[test]
    fn every_cell_carries_its_own_read_groups_and_its_own_candidates_numbers() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        for group in 0..READ_GROUPS as u32 {
            for (candidate, repeat_count) in CANDIDATE_REPEATS.iter().copied().enumerate() {
                let cell = table.of(ReadGroupId(group), candidate);
                assert_eq!(
                    *cell.stutter,
                    stutter_model_for(&slippage_at(group, repeat_count)),
                    "read group {group}, candidate {candidate} ({repeat_count} repeats) is \
                     scored under another cell's slippage"
                );
                assert_eq!(
                    cell.substitution_rate.get(),
                    substitution_rate_at(group, repeat_count),
                    "read group {group}, candidate {candidate} is scored under another cell's \
                     substitution rate"
                );
                assert_eq!(cell.weakest_provenance, Provenance::FittedHere);
            }
        }
        assert_eq!(gathered.cell_count(), READ_GROUPS * CANDIDATE_REPEATS.len());
        assert_eq!(gathered.cells_with_no_fitted_slippage(), 0);
        assert_eq!(gathered.cells_with_no_fitted_substitution_rate(), 0);
        assert_eq!(gathered.weakest_warrant(), Provenance::FittedHere);
    }

    /// **The unreachable mass is the candidate's own**, so a context pairing one candidate's
    /// mass with another's model is visible without scoring anything.
    ///
    /// The two candidates are 6 and 11 repeats of one dinucleotide, and the mass the stutter
    /// distribution cannot place differs between them: **2.02 in 10,000 at 6 repeats against
    /// 1.85 in a million at 11**, about 109-fold, because a shorter tract runs out of repeats
    /// to lose sooner.
    #[test]
    fn each_context_carries_its_own_candidates_unreachable_mass() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        let short = table.of(ReadGroupId(0), 0).unreachable_mass;
        let long = table.of(ReadGroupId(0), 1).unreachable_mass;
        assert!(
            short > long * 50.0,
            "the shorter candidate loses about 109 times as much mass as the longer one, \
             measured: {short} against {long}"
        );
    }

    /// **A candidate whose stratum the fit never reached is scored under the stated constant,
    /// and the cell says so.**
    ///
    /// The fit here holds only the 6-repeat stratum, so all three read groups' cells for the
    /// 11-repeat candidate fall through — which is the ordinary case `NoSlippage::NoSuchStratum`
    /// documents, not a broken run.
    #[test]
    fn a_candidate_the_fit_never_reached_takes_the_stated_constant_and_says_so() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = StratumFits::over(
            &[stratum_row(CANDIDATE_REPEATS[0], READ_GROUPS, |group| {
                fitted_level(400.0 + f64::from(group))
            })],
            one_group_each(READ_GROUPS),
        );
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        for group in 0..READ_GROUPS as u32 {
            assert_eq!(
                *table.of(ReadGroupId(group), 0).stutter,
                stutter_model_for(&slippage_at(group, CANDIDATE_REPEATS[0])),
                "the fitted candidate lost its own numbers"
            );
            assert_eq!(
                *table.of(ReadGroupId(group), 1).stutter,
                StutterModel::hipstr_shipped(),
                "the candidate with no stratum should take the stated constant"
            );
            assert_eq!(
                table.of(ReadGroupId(group), 1).weakest_provenance,
                Provenance::Defaulted
            );
        }
        assert_eq!(gathered.cells_with_no_fitted_slippage(), READ_GROUPS);
        assert_eq!(gathered.weakest_warrant(), Provenance::Defaulted);
    }

    /// **A read group the run has and the fit never named takes the stated constant too**, and
    /// it is a different absence from the one above: the stratum is there, the library is not.
    ///
    /// Both come back `Defaulted` because there is one rung below the fit, and both are counted
    /// together — what tells them apart is `NoSlippage`, which the run's own report reads.
    #[test]
    fn a_read_group_the_fit_never_named_takes_the_stated_constant() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        // Two of the run's three read groups are in the fit's map; the third is a library the
        // pre-pass did not know existed.
        let fits = StratumFits::over(
            &[
                stratum_row(CANDIDATE_REPEATS[0], 2, |group| {
                    fitted_level(400.0 + f64::from(group))
                }),
                stratum_row(CANDIDATE_REPEATS[1], 2, |group| {
                    fitted_level(800.0 + f64::from(group))
                }),
            ],
            one_group_each(2),
        );
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        for candidate in 0..CANDIDATE_REPEATS.len() {
            assert_eq!(
                *table.of(ReadGroupId(2), candidate).stutter,
                StutterModel::hipstr_shipped(),
                "the unknown library's cells should take the stated constant"
            );
        }
        assert_eq!(
            *table.of(ReadGroupId(1), 1).stutter,
            stutter_model_for(&slippage_at(1, CANDIDATE_REPEATS[1])),
            "a library the fit does name keeps its own numbers"
        );
        assert_eq!(
            gathered.cells_with_no_fitted_slippage(),
            CANDIDATE_REPEATS.len()
        );
    }

    /// **A cell with no fitted substitution rate takes
    /// [`DEFAULT_SSR_SUBSTITUTION_RATE`] and is marked defaulted**, while its slippage stays
    /// the fit's.
    ///
    /// The two lookups are independent, and this is what says so: the slippage numbers below
    /// are all the fit's own and only the rate falls through.
    #[test]
    fn a_cell_with_no_fitted_substitution_rate_takes_the_stated_constant_and_says_so() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        // Only read group 0's two cells were measured.
        let mut rates = all_substitution_rates();
        rates.retain(|key, _| key.read_group == ReadGroupId(0));
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        assert_eq!(
            table.of(ReadGroupId(0), 1).substitution_rate.get(),
            substitution_rate_at(0, CANDIDATE_REPEATS[1])
        );
        assert_eq!(
            table.of(ReadGroupId(0), 1).weakest_provenance,
            Provenance::FittedHere
        );
        for group in 1..READ_GROUPS as u32 {
            for (candidate, repeat_count) in CANDIDATE_REPEATS.iter().copied().enumerate() {
                let cell = table.of(ReadGroupId(group), candidate);
                assert_eq!(cell.substitution_rate.get(), DEFAULT_SSR_SUBSTITUTION_RATE);
                assert_eq!(cell.weakest_provenance, Provenance::Defaulted);
                assert_eq!(
                    *cell.stutter,
                    stutter_model_for(&slippage_at(group, repeat_count)),
                    "the slippage lookup is independent of the rate lookup and must not have \
                     fallen through with it"
                );
            }
        }
        assert_eq!(
            gathered.cells_with_no_fitted_substitution_rate(),
            (READ_GROUPS - 1) * CANDIDATE_REPEATS.len()
        );
        assert_eq!(gathered.cells_with_no_fitted_slippage(), 0);
        assert_eq!(gathered.weakest_warrant(), Provenance::Defaulted);
    }

    /// **A level read off its period's curve is borrowed, not fitted here** — the numbers are
    /// the fit's, and they are a neighbouring grain's measurement rather than this cell's.
    ///
    /// The distinction is the one spec §4.4 says the run's output must be able to make: *"a
    /// genotype resting on a direction split borrowed from two repeat counts away is
    /// distinguishable in the run's output from one resting on a fit"*.
    #[test]
    fn a_level_read_off_the_curve_is_borrowed_rather_than_fitted_here() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = StratumFits::over(
            &[
                stratum_row(CANDIDATE_REPEATS[0], READ_GROUPS, |group| {
                    fitted_level(400.0 + f64::from(group))
                }),
                stratum_row_at(
                    MOTIF.len() as u8,
                    CANDIDATE_REPEATS[1],
                    READ_GROUPS,
                    |_| level_off_the_curve(),
                    ShareSource::Curve,
                    ShareSource::Curve,
                ),
            ],
            one_group_each(READ_GROUPS),
        );
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        assert_eq!(
            table.of(ReadGroupId(0), 0).weakest_provenance,
            Provenance::FittedHere
        );
        assert_eq!(
            table.of(ReadGroupId(0), 1).weakest_provenance,
            Provenance::Borrowed
        );
        // **The numbers are still the curve's own**, and are used exactly as a fitted cell's
        // — spec §4.4 forbids down-weighting a borrowed value.
        assert_eq!(
            *table.of(ReadGroupId(0), 1).stutter,
            stutter_model_for(&slippage_at(0, CANDIDATE_REPEATS[1]))
        );
        assert_eq!(gathered.cells_with_no_fitted_slippage(), 0);
        assert_eq!(gathered.weakest_warrant(), Provenance::Borrowed);
    }

    /// **A cell whose level is its own and whose shares came off a curve is borrowed too.**
    ///
    /// The two are smoothed on separate curves and the fit reports them separately, so a
    /// warrant that read only the level would call this cell fitted — and the direction split
    /// it is scoring reads under would be two repeat counts away with nothing saying so.
    #[test]
    fn a_cell_that_borrowed_only_its_shares_is_borrowed() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = StratumFits::over(
            &[
                stratum_row_at(
                    MOTIF.len() as u8,
                    CANDIDATE_REPEATS[0],
                    READ_GROUPS,
                    |group| fitted_level(400.0 + f64::from(group)),
                    ShareSource::Curve,
                    ShareSource::Curve,
                ),
                stratum_row(CANDIDATE_REPEATS[1], READ_GROUPS, |group| {
                    fitted_level(800.0 + f64::from(group))
                }),
            ],
            one_group_each(READ_GROUPS),
        );
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        assert_eq!(
            table.of(ReadGroupId(0), 0).weakest_provenance,
            Provenance::Borrowed,
            "a cell whose direction split came off a curve has not fitted its own shares"
        );
        assert_eq!(
            table.of(ReadGroupId(0), 1).weakest_provenance,
            Provenance::FittedHere
        );
    }

    /// **The lengths the junk term is spread over come from the candidates and the two slip
    /// cutoffs, and from nothing any sample showed.**
    ///
    /// Measured on this fixture — a dinucleotide with candidates of 6 and 11 repeats: the
    /// support holds **41 lengths, every whole number of bases from 2 to 42**, against the two
    /// the candidates themselves spell (12 bases and 22). That width is the shape rather than a
    /// defect: it is every length either candidate can be stretched or trimmed to under the two
    /// slip cutoffs.
    #[test]
    fn the_junk_terms_lengths_come_from_the_candidates_and_the_cutoffs() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );

        let lengths = gathered.reachable_lengths();
        assert_eq!(lengths.len(), 41);
        assert_eq!(*lengths.first().expect("a non-empty support"), 2);
        assert_eq!(*lengths.last().expect("a non-empty support"), 42);
        assert!(lengths.windows(2).all(|pair| pair[0] < pair[1]));
        for count in CANDIDATE_REPEATS {
            let spelled = count * MOTIF.len() as u32;
            assert!(
                lengths.contains(&spelled),
                "a candidate's own length {spelled} must be reachable"
            );
        }
    }

    /// **Gathering a second tract leaves nothing of the first** — the buffers are cleared and
    /// refilled, which is what makes one of these per worker rather than one per locus.
    ///
    /// **The order is the trap, and the first draft had it backwards.** With a fully fitted
    /// tract first and a defaulting one second, a dropped `clear()` on the rates or the
    /// warrants is invisible — the stale entries sit past the ones the second tract reads — and
    /// a dropped counter reset adds nothing to zero. Reviewed by mutation: deleting
    /// `substitution_rate.clear()`, `warrant.clear()`, or either counter reset left all
    /// nineteen tests green. **So the defaulting tract goes first.** Each deletion now fails
    /// here: a stale rate reads 0.001 where the fit says 0.0016, a stale warrant reads
    /// `Defaulted` where the cell is fitted, and a stale counter reads 3 where the second tract
    /// defaults nothing.
    ///
    /// **The reachable support needs the same care**: both tracts reach 41 lengths, so a length
    /// buffer appended to rather than cleared would still be a legal width. What separates them
    /// is *which* lengths — 20 to 60 for the 20-repeat tract and 2 to 42 for the other.
    #[test]
    fn gathering_a_second_tract_leaves_nothing_of_the_first() {
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        // A tract of 20 repeats: no stratum in the fit, no rate in the map, so every cell of it
        // is defaulted on both lookups.
        let lone_bases = tract(20);
        let lone = [SsrCandidate {
            bases: &lone_bases,
            repeat_count: repeats(20),
        }];
        let mut gathered = TractScoringFits::default();
        // Its own reference count, because it is its own tract — the gather refuses a
        // reference count that is not its first candidate's.
        gathered.gather_for_locus(
            motif(),
            &lone,
            tract_prior(RepeatCount(20), &parameters),
            &parameters,
        );
        assert_eq!(gathered.cell_count(), READ_GROUPS);
        assert_eq!(gathered.cells_with_no_fitted_slippage(), READ_GROUPS);
        assert_eq!(
            gathered.cells_with_no_fitted_substitution_rate(),
            READ_GROUPS
        );
        assert_eq!(gathered.reachable_lengths().len(), 41);
        assert!(gathered.reachable_lengths().contains(&60));
        assert!(!gathered.reachable_lengths().contains(&12));

        let bases = candidate_bases();
        let alleles = candidates(&bases);
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );

        assert_eq!(gathered.cell_count(), READ_GROUPS * CANDIDATE_REPEATS.len());
        assert_eq!(gathered.cells_with_no_fitted_slippage(), 0);
        assert_eq!(gathered.cells_with_no_fitted_substitution_rate(), 0);
        assert_eq!(gathered.reachable_lengths().len(), 41);
        assert!(gathered.reachable_lengths().contains(&12));
        assert!(
            !gathered.reachable_lengths().contains(&60),
            "the first tract's own lengths must not survive into the second's support"
        );

        // **Read back per cell, not only in aggregate.** The counters and the widths above are
        // what a stale rate or a stale warrant slips past.
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());
        for group in 0..READ_GROUPS as u32 {
            for (candidate, repeat_count) in CANDIDATE_REPEATS.iter().copied().enumerate() {
                let cell = table.of(ReadGroupId(group), candidate);
                assert_eq!(
                    cell.substitution_rate.get(),
                    substitution_rate_at(group, repeat_count),
                    "read group {group}, candidate {candidate} carries a rate from the first \
                     tract"
                );
                assert_eq!(
                    cell.weakest_provenance,
                    Provenance::FittedHere,
                    "read group {group}, candidate {candidate} carries a warrant from the \
                     first tract"
                );
                assert_eq!(
                    *cell.stutter,
                    stutter_model_for(&slippage_at(group, repeat_count))
                );
            }
        }
    }

    /// One observation of a tract, seen by `reads` reads that spanned the whole of it.
    fn spanning(bases: &[u8], reads: u32, read_group: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.to_vec().into_boxed_slice(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(read_group),
            num_obs: reads,
            num_fwd: reads,
            q_sum: -10.0 * f64::from(reads),
            mapq_sum: 60 * reads,
            mapq_sum_sq: u64::from(reads) * 3_600,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn detail() -> SsrDetail {
        SsrDetail {
            motif: *motif(),
            left_flank: Box::from(&b"GGGGG"[..]),
            right_flank: Box::from(&b"TTTTT"[..]),
        }
    }

    /// Score one sample's reads at the tract under what this module assembled.
    fn score_row(
        gathered: &TractScoringFits,
        alleles: &[SsrCandidate<'_>],
        observations: &[SequenceObservation],
    ) -> Vec<LogProb> {
        let detail = detail();
        let evidence = SsrSampleEvidence::new(observations, &detail);
        let contexts = gathered.scoring_contexts(alleles);
        let table = GenotypeTable::build(diploid(), alleles.len());
        let view = table.view();
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
        genotype_log_likelihood_row(
            &StutterSubstitutionEmission,
            &evidence,
            gathered.locus_parameters(alleles, &contexts, &[]),
            &view,
            &mut out,
            &mut scratch,
        );
        out
    }

    /// **The assembled parameters score a real row, and the reads pick the genotype** — the
    /// whole point of the step, in one test.
    ///
    /// Twenty reads all showing the 11-repeat candidate's tract, at a diploid biallelic locus
    /// whose genotypes are `0/0`, `0/1`, `1/1` in table order. The homozygous 11-repeat call
    /// wins: measured **−2.51 nats against −16.37 for the heterozygote and −161.35 for the
    /// homozygous 6-repeat call**, so the reads separate the three by 14 and 159 nats and
    /// nothing in the assembly is a placeholder that would leave them alike.
    #[test]
    fn the_assembled_parameters_score_a_row_and_the_reads_pick_the_genotype() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );

        let observations = [spanning(&tract(CANDIDATE_REPEATS[1]), 20, 0)];
        let row = score_row(&gathered, &alleles, &observations);

        assert_eq!(row.len(), 3);
        assert!(
            row.iter().all(|score| score.get().is_finite()),
            "every genotype gets a finite score: {row:?}"
        );
        let best = row
            .iter()
            .enumerate()
            .max_by(|left, right| left.1.get().total_cmp(&right.1.get()))
            .expect("a row of three")
            .0;
        assert_eq!(
            best, 2,
            "twenty reads of the 11-repeat tract should pick the homozygous 11-repeat \
             genotype, and the row scored {row:?}"
        );
    }

    /// **The candidate axis of the lookup reaches the row's own numbers.**
    ///
    /// The same reads scored twice: once under the fit above, where the two candidates sit in
    /// different strata and slip at different rates, and once under a fit that gives both
    /// candidates the 6-repeat stratum's numbers. The heterozygous genotype's score moves by
    /// **0.221 nats** — 0.96 Phred, on two slippage levels 0.044 and 0.064 apart — so the
    /// candidate axis reaches the row's own arithmetic and not only its parameter table.
    #[test]
    fn scoring_both_candidates_under_one_stratum_changes_the_row() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let rates = all_substitution_rates();

        let fits = fits_for_both_candidates();
        let parameters = run(&calibration, &inbreeding, &fits, &rates);
        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );

        // The same table with the longer candidate's stratum carrying the shorter one's
        // numbers — what a lookup keyed by the reference tract rather than by the candidate
        // would produce.
        let flattened = StratumFits::over(
            &[
                stratum_row(CANDIDATE_REPEATS[0], READ_GROUPS, |group| {
                    fitted_level(400.0 + f64::from(group))
                }),
                {
                    let mut row = stratum_row(CANDIDATE_REPEATS[1], READ_GROUPS, |group| {
                        fitted_level(800.0 + f64::from(group))
                    });
                    if let StratumOutcome::Derived(derived) = &mut row {
                        for (group, slot) in derived.slippage.iter_mut().enumerate() {
                            *slot = Some(slippage_at(group as u32, CANDIDATE_REPEATS[0]));
                        }
                    }
                    row
                },
            ],
            one_group_each(READ_GROUPS),
        );
        let flat_parameters = run(&calibration, &inbreeding, &flattened, &rates);
        let mut flat_gathered = TractScoringFits::default();
        flat_gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &flat_parameters),
            &flat_parameters,
        );

        let observations = [
            spanning(&tract(CANDIDATE_REPEATS[0]), 10, 0),
            spanning(&tract(CANDIDATE_REPEATS[1]), 10, 0),
        ];
        let row = score_row(&gathered, &alleles, &observations);
        let flat_row = score_row(&flat_gathered, &alleles, &observations);

        let moved = (row[1].get() - flat_row[1].get()).abs();
        assert!(
            moved > 0.1,
            "the candidate's own stratum must reach the row, and dropping the axis moves the \
             heterozygote by exactly nothing: it moved {moved} nats, against 0.221 measured"
        );
    }

    /// **The outlier weight is the inherited constant, and it is in no cell's warrant.**
    ///
    /// Every cell of this fixture is fitted, and the locus's warrant comes back `FittedHere`
    /// even though the weight the row is spread over is defaulted everywhere — which is the
    /// decision this module's own documentation argues for, and the thing that would break
    /// silently if the constant were folded into the fold.
    #[test]
    fn the_outlier_weight_is_the_inherited_constant_and_is_in_no_warrant() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let locus = gathered.locus_parameters(&alleles, &contexts, &[]);

        assert_eq!(locus.outlier_weight, DEFAULT_OUTLIER_WEIGHT);
        assert_eq!(locus.outlier_weight, 0.01);
        assert!(locus.contamination.is_none());
        assert_eq!(gathered.weakest_warrant(), Provenance::FittedHere);
    }

    /// **A tract of another period reaches both lookups with its own period** — and this test
    /// exists because it is the one axis every other fixture holds fixed.
    ///
    /// Found by mutation: the module's motif is a dinucleotide, `MOTIF.len()` is 2, and 2 is
    /// also the candidate count — so replacing `motif.ssr_period()` with a hard-coded period of
    /// 2 left all nineteen other tests green. On a trinucleotide it defaults every cell of the
    /// table: **six cells of six take the stated constants where the fit has numbers for all
    /// six.**
    ///
    /// It also runs the module against a period it has otherwise never seen, which `CLAUDE.md`
    /// asks for on its own account — mononucleotide and trinucleotide tracts are not the corner
    /// case a dinucleotide fixture makes them look like.
    #[test]
    fn a_tract_of_another_period_is_scored_under_its_own_period() {
        let trinucleotide = Motif::new(b"ACG").expect("a valid trinucleotide motif");
        let bases: Vec<Vec<u8>> = CANDIDATE_REPEATS
            .iter()
            .map(|count| b"ACG".repeat(*count as usize))
            .collect();
        let alleles: Vec<SsrCandidate<'_>> = bases
            .iter()
            .zip(CANDIDATE_REPEATS)
            .map(|(bases, count)| SsrCandidate {
                bases,
                repeat_count: repeats(count),
            })
            .collect();

        let fits = StratumFits::over(
            &[
                stratum_row_at(
                    3,
                    CANDIDATE_REPEATS[0],
                    READ_GROUPS,
                    |group| fitted_level(400.0 + f64::from(group)),
                    ShareSource::Stratum,
                    ShareSource::Stratum,
                ),
                stratum_row_at(
                    3,
                    CANDIDATE_REPEATS[1],
                    READ_GROUPS,
                    |group| fitted_level(800.0 + f64::from(group)),
                    ShareSource::Stratum,
                    ShareSource::Stratum,
                ),
            ],
            one_group_each(READ_GROUPS),
        );
        let mut rates = BTreeMap::new();
        for group in 0..READ_GROUPS as u32 {
            for repeat_count in CANDIDATE_REPEATS {
                rates.insert(
                    StratumKey {
                        read_group: ReadGroupId(group),
                        stratum: SsrStratum::new(
                            trinucleotide.ssr_period(),
                            RepeatCount(repeat_count),
                        ),
                        ploidy: diploid(),
                    },
                    Estimate {
                        value: ErrorRate::try_new(substitution_rate_at(group, repeat_count))
                            .expect("a valid rate"),
                        provenance: Provenance::FittedHere,
                        observations: 40_000,
                    },
                );
            }
        }
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            &trinucleotide,
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );

        assert_eq!(
            gathered.cells_with_no_fitted_slippage(),
            0,
            "a trinucleotide tract's own period must reach the slippage lookup"
        );
        assert_eq!(
            gathered.cells_with_no_fitted_substitution_rate(),
            0,
            "a trinucleotide tract's own period must reach the substitution-rate lookup"
        );
        assert_eq!(gathered.weakest_warrant(), Provenance::FittedHere);
    }

    /// **A cell that fitted one direction share and read the other off a curve is borrowed.**
    ///
    /// Found by mutation: the fixture helper used to set both shares from one source, so
    /// dropping either half of `warrant_of`'s test left all nineteen tests green. The fit
    /// reports the two separately — a stratum can fit its contraction bias on its own reads and
    /// take its fall-off from its period's curve, which needs about three times as much
    /// evidence — so a cell in that state is real, and calling it fitted would say the run
    /// measured something it read off a neighbour.
    #[test]
    fn a_cell_that_fitted_one_share_and_borrowed_the_other_is_borrowed() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = StratumFits::over(
            &[
                stratum_row_at(
                    MOTIF.len() as u8,
                    CANDIDATE_REPEATS[0],
                    READ_GROUPS,
                    |group| fitted_level(400.0 + f64::from(group)),
                    ShareSource::Stratum,
                    ShareSource::Curve,
                ),
                stratum_row_at(
                    MOTIF.len() as u8,
                    CANDIDATE_REPEATS[1],
                    READ_GROUPS,
                    |group| fitted_level(800.0 + f64::from(group)),
                    ShareSource::Curve,
                    ShareSource::Stratum,
                ),
            ],
            one_group_each(READ_GROUPS),
        );
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        assert_eq!(
            table.of(ReadGroupId(0), 0).weakest_provenance,
            Provenance::Borrowed,
            "a fall-off read off a curve is not a cell that fitted its own shares"
        );
        assert_eq!(
            table.of(ReadGroupId(0), 1).weakest_provenance,
            Provenance::Borrowed,
            "a contraction bias read off a curve is not a cell that fitted its own shares"
        );
    }

    /// **A cell with slippage numbers and no shares provenance at all is borrowed** — the
    /// branch `warrant_of` documents and that no `StratumFits` the fit itself builds can
    /// produce, since both of its paths set the numbers and the shares from one mask.
    ///
    /// Tested against the function directly, because the state is not reachable through
    /// `gather_for_locus`: the numbers are the fit's, so `Defaulted` would be wrong, and the
    /// cell cannot claim its shares were fitted here.
    #[test]
    fn a_cell_whose_shares_have_no_provenance_is_borrowed() {
        let fitted_shares = FittedSlippage {
            slippage: slippage_at(0, CANDIDATE_REPEATS[0]),
            level: fitted_level(400.0),
            shares: Some(shares_from(ShareSource::Stratum, ShareSource::Stratum)),
        };
        assert_eq!(warrant_of(&fitted_shares), Provenance::FittedHere);
        let no_shares = FittedSlippage {
            shares: None,
            ..fitted_shares
        };
        assert_eq!(warrant_of(&no_shares), Provenance::Borrowed);
    }

    /// **A substitution rate the fit did not measure here carries its own warrant through.**
    ///
    /// The pre-pass emits this rate as `FittedHere` or not at all today, so no run can reach
    /// this — but the slippage side already has a borrowed rung and this one is a candidate for
    /// the same. Found by mutation: writing `Provenance::FittedHere` in place of the estimate's
    /// own warrant left every other test green, so the day the rate gains a rung the assembly
    /// would have laundered it.
    #[test]
    fn a_rate_the_fit_borrowed_is_carried_as_borrowed() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let mut rates = all_substitution_rates();
        let borrowed = StratumKey {
            read_group: ReadGroupId(1),
            stratum: SsrStratum::new(period(), RepeatCount(CANDIDATE_REPEATS[0])),
            ploidy: diploid(),
        };
        rates
            .get_mut(&borrowed)
            .expect("the fixture measures every cell")
            .provenance = Provenance::Borrowed;
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());

        assert_eq!(
            table.of(ReadGroupId(1), 0).weakest_provenance,
            Provenance::Borrowed
        );
        assert_eq!(
            table.of(ReadGroupId(1), 1).weakest_provenance,
            Provenance::FittedHere,
            "the cell beside it measured its own rate"
        );
        assert_eq!(
            gathered.cells_with_no_fitted_substitution_rate(),
            0,
            "a borrowed rate is a measurement, not an absence"
        );
        assert_eq!(gathered.weakest_warrant(), Provenance::Borrowed);
    }

    /// **A run of one library is gathered like any other** — and it is every sample of both
    /// benchmark cohorts here, where every other fixture in this module has three.
    ///
    /// At one read group the table is one row, so a cell index and a candidate index are the
    /// same number; nothing here can tell a transposition apart, which is why the three-group
    /// fixtures exist beside it. What this pins is that the one-library shape works at all.
    #[test]
    fn a_run_of_one_library_is_gathered_like_any_other() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = StratumFits::over(
            &[
                stratum_row(CANDIDATE_REPEATS[0], 1, |group| {
                    fitted_level(400.0 + f64::from(group))
                }),
                stratum_row(CANDIDATE_REPEATS[1], 1, |group| {
                    fitted_level(800.0 + f64::from(group))
                }),
            ],
            one_group_each(1),
        );
        let rates = all_substitution_rates();
        let calibration = vec![ReadGroupCalibration::defaulted()];
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        assert_eq!(gathered.cell_count(), CANDIDATE_REPEATS.len());
        assert_eq!(gathered.cells_with_no_fitted_slippage(), 0);
        assert_eq!(gathered.weakest_warrant(), Provenance::FittedHere);

        let contexts = gathered.scoring_contexts(&alleles);
        let table = SsrScoringContextTable::new(&contexts, alleles.len());
        assert_eq!(table.read_group_count(), 1);
        assert_eq!(
            *table.of(ReadGroupId(0), 1).stutter,
            stutter_model_for(&slippage_at(0, CANDIDATE_REPEATS[1]))
        );
    }

    /// **A library the fit's own map does not name is counted apart from a stratum that has no
    /// numbers**, because the two absences mean different things: the first says the parameters
    /// and the reads came from different runs, the second is ordinary.
    #[test]
    fn a_library_the_fit_does_not_describe_is_counted_apart() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        // **The two counts must come out different, and the first draft of this fixture made
        // them equal — two of each — so narrowing the matcher to the ordinary absence still
        // reported the right number.** The fit here names one of the run's three libraries and
        // holds one of the two candidates' strata, which gives **four** cells lost to an
        // unknown library against **one** lost to a stratum with no numbers.
        let fits = StratumFits::over(
            &[stratum_row(CANDIDATE_REPEATS[0], 1, |group| {
                fitted_level(400.0 + f64::from(group))
            })],
            one_group_each(1),
        );
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );

        // Read groups 1 and 2 are unknown to the fit, so all four of their cells are lost that
        // way — an unknown library is answered before the stratum is even looked up. Read group
        // 0 is named, and loses only the 11-repeat candidate, whose stratum the fit lacks.
        assert_eq!(gathered.cells_with_no_fitted_slippage(), 5);
        assert_eq!(
            gathered.cells_whose_read_group_the_fit_does_not_describe(),
            4,
            "only the unknown libraries' cells are the kind of absence a run should act on"
        );
        assert_eq!(gathered.weakest_warrant(), Provenance::Defaulted);
    }

    /// **A run whose fit found contamination gets the three-term form, and the third term is a
    /// distribution over this tract's reachable lengths.**
    ///
    /// Spec §4.5.1 puts it on wherever the pre-pass emits a fraction. What it must be is a
    /// probability over the *same* support the outlier term is spread over — otherwise the row's
    /// three terms are not probabilities of one event — so this checks the width, the total, and
    /// where the mass actually sits.
    ///
    /// **The mass sits at the candidates' byte lengths and nowhere else.** This fixture's two
    /// candidates hold 6 and 11 whole `AC` repeats, so they spell **12 and 22 bases**, and the
    /// reachable support is far wider than that — every length the slip cutoffs admit from
    /// either. Every other entry is zero, which is the shape rather than a defect: a
    /// contaminating read at a length no candidate carries falls to the outlier floor, which is
    /// where §4.5.1 puts it.
    ///
    /// The fit behind this fixture carries no length spectrum, so the prior answers from the
    /// ladder's bottom rung — a flat shape — and the two candidates take **half each**. A fit
    /// with a spectrum would divide it unevenly; what would not change is the width or the
    /// total.
    #[test]
    fn a_contaminated_run_gets_a_seed_over_the_lengths_its_candidates_reach() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(READ_GROUPS);
        let views = vec![a_contamination_view(0.03); READ_GROUPS];
        let batching = crate::ng::calling::tests::one_batch(READ_GROUPS, READ_GROUPS);
        let parameters = FrozenParameters::new(
            &calibration,
            &views,
            &batching,
            &inbreeding,
            SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
            &fits,
            &rates,
            diploid(),
        );

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let locus = gathered.locus_parameters(&alleles, &contexts, &views);

        let contamination = locus
            .contamination
            .expect("a run whose fit found a fraction is scored on the three-term form");
        let seed = contamination.contaminant_length_frequencies;
        assert_eq!(
            seed.len(),
            gathered.reachable_lengths().len(),
            "the seed is keyed to the same support the outlier term is spread over"
        );
        assert!(
            (seed.iter().sum::<f64>() - 1.0).abs() < 1e-12,
            "the seed is how common each length is, so it sums to one: {seed:?}"
        );
        assert!(
            gathered.reachable_lengths().len() > 4,
            "the support must be wider than the candidate set for this test to say anything"
        );

        for (length, share) in gathered.reachable_lengths().iter().zip(seed) {
            let a_candidates_own_length = *length == 12 || *length == 22;
            assert_eq!(
                *share > 0.0,
                a_candidates_own_length,
                "length {length} carries {share} of the seed"
            );
            if a_candidates_own_length {
                assert!(
                    (*share - 0.5).abs() < 1e-12,
                    "a flat prior over two candidates gives each half: {share}"
                );
            }
        }
    }

    /// **A reference repeat count that is not the reference allele's is refused**, and this is
    /// the one mistake at a tract that would otherwise come back as a plausible number.
    ///
    /// The prior's length spectrum is indexed by offset from the tract's *reference* length, so
    /// a candidate's count passed in its place re-centres the whole shape on that candidate: the
    /// seed still sums to one, every genotype still gets a prior, and the locus is called under
    /// a belief about which lengths are common that the fit never expressed.
    /// `fill_seed_share_per_candidate`'s own documentation lists it among the mistakes nothing
    /// inside it can catch. Here the candidate table is in hand and it is catchable.
    #[test]
    #[should_panic(expected = "re-centres the fitted length spectrum")]
    fn a_reference_count_that_is_not_the_reference_alleles_is_refused() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(READ_GROUPS);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        // The *second* candidate's count, where the reference allele's belongs.
        let wrong = RepeatCount(CANDIDATE_REPEATS[1]);
        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(wrong, &parameters),
            &parameters,
        );
    }

    /// **The other direction of the same check, and it is the one that loses a measurement.**
    /// An uncontaminated gather handed a run's fitted fractions would return the two-term form
    /// and drop them — a genotype computed as though the library were clean, with nothing in
    /// the output saying so, which is the failure spec §3.6 exists to prevent.
    #[test]
    #[should_panic(expected = "came from different runs")]
    fn an_uncontaminated_gather_handed_fractions_is_refused() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(READ_GROUPS);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let views = vec![a_contamination_view(0.03); READ_GROUPS];
        let _ = gathered.locus_parameters(&alleles, &contexts, &views);
    }

    /// **A fraction list that is not one per read group of the run is refused**, which the
    /// emptiness check above cannot catch.
    ///
    /// A short list is the dangerous shape: the row indexes it by [`ReadGroupId`], so it would
    /// score every library the list does cover and fail only if a read from one past its end
    /// happened to reach this tract.
    #[test]
    #[should_panic(expected = "describes a different run")]
    fn a_fraction_list_of_the_wrong_width_is_refused() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(READ_GROUPS);
        let views = vec![a_contamination_view(0.03); READ_GROUPS];
        let batching = crate::ng::calling::tests::one_batch(READ_GROUPS, READ_GROUPS);
        let parameters = FrozenParameters::new(
            &calibration,
            &views,
            &batching,
            &inbreeding,
            SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
            &fits,
            &rates,
            diploid(),
        );

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let _ = gathered.locus_parameters(&alleles, &contexts, &views[..1]);
    }

    /// **The seed and the fractions must come from one run.** The seed is built, or left empty,
    /// from the same predicate that says whether the fractions exist, so a disagreement means
    /// the two were assembled against different runs — and the row would then scale an empty
    /// distribution by a fraction, or hold a distribution nothing scales.
    #[test]
    #[should_panic(expected = "came from different runs")]
    fn a_contaminated_gather_handed_no_fractions_is_refused() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(READ_GROUPS);
        let views = vec![a_contamination_view(0.03); READ_GROUPS];
        let batching = crate::ng::calling::tests::one_batch(READ_GROUPS, READ_GROUPS);
        let parameters = FrozenParameters::new(
            &calibration,
            &views,
            &batching,
            &inbreeding,
            SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
            &fits,
            &rates,
            diploid(),
        );

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let _ = gathered.locus_parameters(&alleles, &contexts, &[]);
    }

    /// **Nothing gathered means no warrant to report**, rather than the strongest rung on the
    /// ladder from having read nothing.
    ///
    /// `TractScoringFits` derives `Default` and is public, so this state is reachable from
    /// outside; the fold's identity is `FittedHere`, which is the answer that would otherwise
    /// come back for a tract nobody read a parameter for.
    #[test]
    #[should_panic(expected = "no tract has been gathered")]
    fn an_ungathered_value_reports_no_warrant() {
        let _ = TractScoringFits::default().weakest_warrant();
    }

    /// **Nor any contexts** — the same state, at the other public entry point, where it
    /// otherwise divided the cell count by a stride of zero.
    #[test]
    #[should_panic(expected = "no tract has been gathered")]
    fn an_ungathered_value_builds_no_contexts() {
        let _ = TractScoringFits::default().scoring_contexts(&[]);
    }

    /// **The defaulted substitution rate is 0.001**, asserted as a literal.
    ///
    /// Every other test compares a cell against the constant, which passes whatever the
    /// constant is set to — so this is the one place that would fail if the number moved
    /// without anyone deciding to move it.
    #[test]
    fn the_defaulted_substitution_rate_is_one_error_in_a_thousand_bases() {
        assert_eq!(DEFAULT_SSR_SUBSTITUTION_RATE, 0.001);
    }

    /// **A tract with no candidates is refused**, rather than gathered into a table with a
    /// stride of zero.
    #[test]
    #[should_panic(expected = "at least its reference allele")]
    fn a_tract_with_no_candidates_is_refused() {
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);
        TractScoringFits::default().gather_for_locus(
            motif(),
            &[],
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
    }

    /// **Contexts asked for a different candidate set than the fits were gathered over are
    /// refused** — the mispairing that would otherwise put one candidate's unreachable mass
    /// beside another's model, as a plausible number.
    #[test]
    #[should_panic(expected = "belongs to a different locus")]
    fn contexts_for_another_candidate_set_are_refused() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let _ = gathered.scoring_contexts(&alleles[..1]);
    }

    /// **A candidate set narrower than the one these fits were gathered over is refused by the
    /// row's parameters**, which is the mistake `scoring_contexts` refuses one call earlier —
    /// restated here because a caller holding contexts from this tract and candidates from
    /// another would otherwise reach the row with the two disagreeing.
    #[test]
    #[should_panic(expected = "candidates reached the row")]
    fn a_narrower_candidate_set_is_refused_by_the_rows_parameters() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let _ = gathered.locus_parameters(&alleles[..1], &contexts, &[]);
    }

    /// **Contexts built at an earlier tract are refused by the row's parameters**, which is a
    /// different mistake from the one above: the candidates match and the contexts do not.
    #[test]
    #[should_panic(expected = "contexts reached the row")]
    fn contexts_from_another_tract_are_refused() {
        let bases = candidate_bases();
        let alleles = candidates(&bases);
        let fits = fits_for_both_candidates();
        let rates = all_substitution_rates();
        let calibration = calibrations();
        let inbreeding = outbred(1);
        let parameters = run(&calibration, &inbreeding, &fits, &rates);

        let mut gathered = TractScoringFits::default();
        gathered.gather_for_locus(
            motif(),
            &alleles,
            tract_prior(reference_repeats(), &parameters),
            &parameters,
        );
        let contexts = gathered.scoring_contexts(&alleles);
        let _ = gathered.locus_parameters(&alleles, &contexts[..alleles.len()], &[]);
    }

    /// **A SNP/indel candidate table is refused as a repeat tract** — a locus routed to the
    /// wrong read model, which would otherwise be scored as stutter over bases that are not a
    /// tract.
    #[test]
    #[should_panic(expected = "routed to the wrong read model")]
    fn a_snp_locus_is_refused_as_a_repeat_tract() {
        let alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
        let _ = tract_candidates(&alleles, &[repeats(1)]);
    }

    /// **A bundle of repeat tracts is refused, and it is not the same refusal as a misrouted
    /// SNP.** The calling seam sends a bundle down the repeat path deliberately — that is how
    /// every consumer of `LocusKind` groups the two — so a bundle reaching here is a gap in
    /// what this path covers, not a locus handed to the wrong read model. The message says
    /// which, because the two ask a reader to do different things.
    ///
    /// Nothing constructs a bundle into calling today; this is what keeps the two messages from
    /// merging back into one when something does.
    #[test]
    #[should_panic(expected = "nothing scores a bundle")]
    fn a_bundle_of_repeat_tracts_is_refused_as_one_tract() {
        let alleles = CandidateAlleles::new(Box::from(tract(6).as_slice()), LocusKind::SsrBundle);
        let _ = tract_candidates(&alleles, &[repeats(6)]);
    }

    /// **A repeat-count list that does not match the candidate table is refused**, because a
    /// short list would silently drop the last candidates and a long one belongs to another
    /// locus.
    #[test]
    #[should_panic(expected = "belongs to a different locus")]
    fn a_repeat_count_list_of_the_wrong_length_is_refused() {
        let alleles =
            CandidateAlleles::new(Box::from(tract(6).as_slice()), LocusKind::Ssr(detail()));
        let _ = tract_candidates(&alleles, &[repeats(6), repeats(11)]);
    }

    /// **The candidates come back in the table's own allele order**, reference first, each
    /// paired with the repeat count supplied for it.
    #[test]
    fn the_candidates_keep_the_tables_allele_order() {
        let mut alleles =
            CandidateAlleles::new(Box::from(tract(6).as_slice()), LocusKind::Ssr(detail()));
        alleles.admit(Box::from(tract(11).as_slice()));
        let built = tract_candidates(&alleles, &[repeats(6), repeats(11)]);
        assert_eq!(built.len(), 2);
        assert_eq!(built[0].repeat_count.get(), 6);
        assert_eq!(built[0].bases, tract(6).as_slice());
        assert_eq!(built[1].repeat_count.get(), 11);
        assert_eq!(built[1].bases, tract(11).as_slice());
    }
}
