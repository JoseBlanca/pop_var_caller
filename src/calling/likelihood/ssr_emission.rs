//! The STR emission seam — the one surface of this step that swaps.
//!
//! *How probable is one observed sequence, given one candidate allele?* At a repeat tract
//! that question has more than one defensible answer, and the caller must be able to run a
//! second one against the default without touching anything around it. This module holds
//! the **seam**: what an emission is handed, what it returns, and nothing about how it
//! decides (`doc/devel/ng/spec/read_likelihoods.md` §2.4, §4.1).
//!
//! Everything the model reads arrives **per call**, in [`SsrScoringContext`]. Nothing is read
//! from global state and nothing is fitted here — the whole point of the seam is that the
//! EM loop can re-estimate the slippage numbers between iterations with no change on this
//! side (spec §6.1).
//!
//! # The two questions a model answers, and why they are separate methods
//!
//! A read that spanned the whole tract pins what the sample carries there. A read that
//! entered the tract and **ran off its own end** proves only that the tract is *at least* as
//! long as what it showed — the statistician's word is **censored**. Those are different
//! questions about the same candidate, so they are [`SsrEmissionModel::emission`] and
//! [`SsrEmissionModel::censored_emission`]. Routing between them is the row's job, from the
//! observation's witness. The censored method is handed the witness's positions too, which only
//! an STR bundle reads: they say which flank the read held.
//!
//! # What a context is built per, and what must not be hoisted
//!
//! **Per `(read group, candidate)`** — not per locus. A read's chance of slipping is a
//! property of the tract it was copied from, and that is the **candidate** allele: a
//! candidate of 6 repeats and one of 12 at the same locus are drawn from different strata and
//! slip at measurably different rates, about 1.3-fold per repeat count over the measured
//! range (spec §4.4). So the stutter parameters cannot be hoisted out of the candidate loop.
//! **The lookup can** — it is a small table indexed by period and repeat count — and that is
//! the distinction a coder has to keep.

use std::num::{NonZeroU8, NonZeroU32};

use crate::alignment::StutterModel;
use crate::alignment::emission::{BaseScores, Emission, FlatEmission};
use crate::alignment::stutter::{MAX_PART_REPEAT_SLIP, MAX_WHOLE_REPEAT_SLIP};
use crate::float;
use crate::locus_generation::{LocusLen, WitnessedLocusPositions};
use crate::parameter_estimation::Provenance;
use crate::region_typing::reference_runs::ReferenceRun;
use crate::region_typing::run_parse::RunParse;
use crate::types::{BaseQual, ErrorRate, Motif};

/// One candidate allele, as the emission sees it.
///
/// Two fields, and the second is not derivable from the first: **the repeat count keys the
/// stratum lookup**, and it is the candidate's own count rather than the reference's
/// (spec §4.4). Counting it from `bases` would mean re-measuring a tract the locus generator
/// has already measured, which is the duplication spec §7 puts on the alignment module's side
/// of the boundary — this type *consumes* a measurement and never makes one. **Nor is it
/// derivable**: an interrupted tract's byte length divided by the period is not its repeat
/// count.
#[derive(Debug, Clone, Copy)]
pub struct SsrCandidate<'a> {
    /// **The tract, and only the tract** — the repeat run as a carrier of this allele has
    /// it, without the flanks.
    ///
    /// **This differs from the generic path deliberately, and the difference is the locus's
    /// own shape.** There an allele is the whole locus as a carrier has it, because the
    /// varying regions sit inside surrounding context that alleles disagree about. Here the
    /// locus *is* the tract: the STR generator slices exactly the tract into an observation's
    /// bases ([`locus_generation/ssr.rs`](../../../../src/locus_generation/ssr.rs)), the
    /// flanks are shared by every candidate by construction, and comparing them would be
    /// comparing two copies of the same context.
    ///
    /// So an observation and a candidate are on the same footing here, which is what lets the
    /// substitution term compare them once the stutter factor has equalised their lengths.
    pub bases: &'a [u8],
    /// How many repeats this candidate's tract holds. **Non-zero**: a candidate whose tract
    /// holds no repeats is not a candidate, which is the same contract
    /// [`StutterModel::unreachable_mass`] states from the distribution's side.
    pub repeat_count: NonZeroU32,
    /// **The flank bases this allele's deletion takes beside the tract** — none but for an
    /// allele that is a deletion running from the repeat into a flank (`str_edge_deletion_alleles.md`
    /// §5.6). A read whose own flank deletion differs is not this allele with stutter: the
    /// emission is zero, and the read is left to the outlier term.
    pub flank_deletion: crate::locus_generation::FlankDeletion,
    /// The locus's reference flanks, left and right — what a read and an allele whose flank
    /// extents differ are compared over (`str_edge_deletion_alleles.md` §5.6). Empty where the
    /// locus carries none, as in the tests that predate flank extents.
    pub reference_flanks: (&'a [u8], &'a [u8]),
    /// **At an STR bundle, the allele's runs** — the bundle's reference runs and this allele's run
    /// parse (`str_bundle_loci.md` §7.1) — and `None` at a segment. Where it is set, the emission
    /// prices stutter in every run rather than in one motif (§7.3), and
    /// [`repeat_count`](Self::repeat_count) is not read.
    pub bundle: Option<BundleCandidate<'a>>,
}

/// **One STR bundle allele's runs, as the emission places stutter in them** — the bundle's
/// reference runs, in hull order, and the allele's run parse: each run's count and its place among
/// the allele's bases (`str_bundle_loci.md` §7.1).
#[derive(Debug, Clone, Copy)]
pub struct BundleCandidate<'a> {
    /// The bundle's reference runs, in hull order — each run's motif.
    pub runs: &'a [ReferenceRun],
    /// The allele read against them: one place and count per run.
    pub parse: &'a RunParse,
    /// **Bases the scored sequence carries before the allele's own first base** — zero, except
    /// where a read and the allele were respelled over a wider flank extent and the allele was
    /// given reference flank bases in front (`fill_emissions` in `likelihood/ssr.rs`). Each run's
    /// place is shifted by it.
    pub leading_flank_bases: u32,
    /// **The hull's length in reference positions** — the reference allele's length, and the axis
    /// a partial read's witness is measured on. It is what says which flank a read that ended
    /// inside the hull held (`str_bundle_loci.md` §7.6).
    pub reference_hull_len: u32,
}

/// **`bases`, read as the sequence over `extent`, respelled over the wider `target` extent**:
/// the reference flank bases between the two are added on each side. `None` when `target` does
/// not contain `extent` or reaches past the flanks the locus holds.
pub(crate) fn respelled_over(
    bases: &[u8],
    extent: crate::locus_generation::FlankDeletion,
    target: crate::locus_generation::FlankDeletion,
    (left_flank, right_flank): (&[u8], &[u8]),
    out: &mut Vec<u8>,
) -> Option<()> {
    let (left, right) = (usize::from(extent.left), usize::from(extent.right));
    let (to_left, to_right) = (usize::from(target.left), usize::from(target.right));
    if to_left < left
        || to_right < right
        || to_left > left_flank.len()
        || to_right > right_flank.len()
    {
        return None;
    }
    out.clear();
    out.extend_from_slice(&left_flank[left_flank.len() - to_left..left_flank.len() - left]);
    out.extend_from_slice(bases);
    out.extend_from_slice(&right_flank[right..to_right]);
    Some(())
}

/// Everything a model is handed for one `(read group, candidate)`, and the only channel it
/// has.
///
/// **The tier-two seam** (spec §6.1): every number here may be re-estimated between the
/// caller's iterations, and the emission never asks where any of them came from. That is what
/// makes the EM loop's re-fitting a change of nobody's code but its own.
#[derive(Debug, Clone, Copy)]
pub struct SsrScoringContext<'a> {
    /// The repeating unit, whose length is the period the stutter distribution is indexed by.
    pub motif: &'a Motif,
    /// How likely each length change is, for **this candidate's** stratum — built per
    /// `(read group, candidate)` and never hoisted out of the candidate loop (spec §4.4).
    pub stutter: &'a StutterModel,
    /// The per-base substitution rate for this read group and stratum.
    ///
    /// **Never the SNP/indel path's ε, and never a read's own summed quality** — spec §4.3's
    /// closed question Q6, and the reason is a unit mismatch rather than a preference. A
    /// read's error probability is a per-*read* number, the chance it is wrong somewhere; the
    /// substitution term needs a per-*base* rate, applied once for each of the tract's twenty
    /// or forty bases. Using the first as the second overcharges by the tract's length.
    ///
    /// The two rates are separate fitted parameters that are never tied: each absorbs what
    /// its own model cannot otherwise explain, and forcing one number to carry both would
    /// make each model wrong in a way neither could report.
    pub substitution_rate: ErrorRate,
    /// The mass the stutter distribution cannot place for this candidate
    /// ([`StutterModel::unreachable_mass`]) — computed and carried, never assumed negligible.
    ///
    /// **It travels because the row compares candidates.** A model that loses mass on some
    /// candidates and not others is comparing them on different scales, and at period 1 the
    /// loss is 2 in 100 rather than the 1-in-10¹³ a cutoff tail costs (spec §4.2).
    pub unreachable_mass: f64,
    /// **At an STR bundle, one stutter model per reference run**, in hull order, each looked up by
    /// this read group, the run's period and the run's count in this candidate
    /// (`str_bundle_loci.md` §7.3); empty at a segment, whose one model is
    /// [`stutter`](Self::stutter). A bundle's `stutter` is its first run's model — what prices a
    /// read whose flank extent differs from the allele's — or the shipped constants where the
    /// bundle has no run.
    pub run_stutter: &'a [StutterModel],
    /// The weakest warrant behind any parameter in this context.
    ///
    /// **The model never branches on it; it propagates** (spec §4.4). A stratum whose numbers
    /// were borrowed is used exactly as a fitted one, with no down-weighting — a borrowed
    /// value is the best estimate available and discounting it would mean inventing a
    /// penalty. But the fact travels, so a call resting on borrowed parameters is
    /// distinguishable in the run's output from one resting on a fit, without re-running
    /// anything.
    pub weakest_provenance: Provenance,
}

impl<'a> SsrScoringContext<'a> {
    /// Build a context for one candidate, taking the unreachable mass from the distribution
    /// rather than from the caller.
    ///
    /// **The mass and the model must agree**, and letting a caller pass both separately is
    /// two chances to disagree. The only inputs that are genuinely the caller's are the
    /// parameters' warrants, which come from the fits the numbers were read out of —
    /// combined here with [`Provenance::weaker_of`], because a context resting on one fitted
    /// number and one borrowed number is a borrowed context.
    #[must_use]
    pub fn new(
        motif: &'a Motif,
        stutter: &'a StutterModel,
        candidate: &SsrCandidate<'_>,
        substitution_rate: ErrorRate,
        parameter_warrants: impl IntoIterator<Item = Provenance>,
    ) -> Self {
        let weakest_provenance = parameter_warrants
            .into_iter()
            .fold(Provenance::FittedHere, Provenance::weaker_of);
        Self {
            motif,
            stutter,
            substitution_rate,
            unreachable_mass: stutter.unreachable_mass(period_of(motif), candidate.repeat_count),
            run_stutter: &[],
            weakest_provenance,
        }
    }

    /// **Build a context for one STR bundle candidate** (`str_bundle_loci.md` §7.3): one stutter
    /// model per reference run, `run_stutter`, in hull order, and `stutter` the model a flank
    /// length change is priced under (see [`Self::run_stutter`]). `motif` is not read by the
    /// bundle emission; it is carried because the field is shared with segments.
    ///
    /// **The unreachable mass is the chance that some run slips where it cannot**, from each run's
    /// own: `1 − Π (1 − loss of run r)` over the runs that hold a copy. It is reported, never read
    /// by the scoring, as for a segment.
    ///
    /// # Panics
    ///
    /// If `candidate` carries no bundle, or `run_stutter` is not one model per run.
    #[must_use]
    pub fn for_bundle(
        motif: &'a Motif,
        stutter: &'a StutterModel,
        run_stutter: &'a [StutterModel],
        candidate: &SsrCandidate<'_>,
        substitution_rate: ErrorRate,
        parameter_warrants: impl IntoIterator<Item = Provenance>,
    ) -> Self {
        let bundle = candidate
            .bundle
            .expect("a bundle's scoring context is built for a bundle candidate");
        assert_eq!(
            run_stutter.len(),
            bundle.runs.len(),
            "a bundle candidate is scored with one stutter model per reference run"
        );
        let weakest_provenance = parameter_warrants
            .into_iter()
            .fold(Provenance::FittedHere, Provenance::weaker_of);
        let kept = bundle
            .runs
            .iter()
            .zip(bundle.parse.runs())
            .zip(run_stutter)
            .filter_map(|((run, placed), model)| {
                NonZeroU32::new(placed.copies)
                    .map(|copies| 1.0 - model.unreachable_mass(period_of(&run.motif()), copies))
            })
            .product::<f64>();
        Self {
            motif,
            stutter,
            substitution_rate,
            unreachable_mass: (1.0 - kept).max(0.0),
            run_stutter,
            weakest_provenance,
        }
    }
}

/// A motif's period as the stutter distribution wants it.
///
/// [`Motif::new`] rejects an empty motif and one longer than six bases, so the period is in
/// `1..=6` by construction and neither conversion can fail — which is what
/// [`StutterModel::probability`]'s doc means when it says "the conversion at a real call site
/// cannot fail". Done once here rather than at every consumer.
#[inline]
fn period_of(motif: &Motif) -> NonZeroU8 {
    let period = u8::try_from(motif.period()).expect("a motif is at most six bases");
    NonZeroU8::new(period).expect("a motif is at least one base")
}

/// **The tract lengths this locus's candidates can reach**, ascending and without repeats,
/// written into a buffer the caller reuses (cleared first).
///
/// # What it is for, and why it is not a count of what anybody showed
///
/// Two of the row's terms are spread over this support. The **outlier weight** is uniform over
/// it — a read from a paralogous tract or a chimera could have shown any of these lengths and
/// nothing about the genotype says which. And the **contamination seed** is a distribution over
/// it, one entry per length.
///
/// **Production spreads its outlier weight over `D`, the number of distinct sequences the whole
/// cohort showed at the locus, and that is the defect this repairs** (spec §4.5). A single
/// sample showing two sequences got a junk floor of 0.005 and a 63-accession panel showing
/// twenty got 0.0005 — ten times lower — so a sample's own genotype likelihood moved when an
/// unrelated sample was added to the run. That is not a property a per-sample likelihood may
/// have.
///
/// **So this asks the candidates and the cutoffs and nothing else.** No observation reaches it,
/// no sample count, and — because it goes through
/// [`StutterModel::reachable_length_changes_of`], which takes no model — none of the fitted
/// rates either, so it does not move between the caller's iterations.
///
/// **What that does and does not buy, measured.** What it removes is the dependence on *what
/// samples showed*: a sample's own row no longer moves when an unrelated sample's reads join
/// the locus. What it does not remove is every trace of the cohort, because the candidate set
/// is itself cohort-derived — a locus that admits one more candidate reaches more lengths.
/// Measured over 1 to 20 candidates, the junk floor swings by about **2.2 to 2.6 fold**, against
/// production's **10 fold** for the same growth in what the cohort showed, and one extra
/// candidate at a five-candidate dinucleotide locus moves it by 5 parts in 100. That is the
/// honest size of the repair: much smaller and no longer a function of the reads, not zero.
///
/// A candidate contributes the lengths its own tract can be stretched or trimmed to. **The
/// guard is against a contraction, not a stretching** — only trimming can take a tract to no
/// bases at all, and the distribution's own reachability rule already stops it one repeat
/// earlier on the whole-repeat branch.
///
/// **An STR bundle candidate reaches the totals of its slip sets and its runs' part-repeat
/// changes** (`str_bundle_loci.md` §7.3); `motif` is not read for it.
pub fn fill_reachable_lengths(candidates: &[SsrCandidate<'_>], motif: &Motif, out: &mut Vec<u32>) {
    out.clear();
    let period = period_of(motif);
    let mut bundle_changes = Vec::new();
    for candidate in candidates {
        let held = candidate.bases.len() as i64;
        let mut reach = |bp_diff: i64| {
            let reached = held + bp_diff;
            if reached > 0 {
                out.push(reached as u32);
            }
        };
        match candidate.bundle {
            Some(bundle) => {
                bundle_changes.clear();
                bundle_reachable_length_changes(bundle, &mut bundle_changes);
                bundle_changes.iter().copied().for_each(&mut reach);
                // A bundle allele with no base still has itself to reach — a read of it is empty.
                if held == 0 {
                    out.push(0);
                }
            }
            None => StutterModel::reachable_length_changes_of(period, candidate.repeat_count)
                .for_each(reach),
        }
    }
    out.sort_unstable();
    out.dedup();
}

/// `Lr(observation | one candidate allele)` — **the only part that differs between models**.
///
/// Everything around it is shared: the copy-weighted mixture over a genotype's alleles, the
/// outlier term, the caching, the logarithm. A second model swaps this trait and nothing
/// else, which is what makes the comparison behind spec §4.1 a fair one.
///
/// # Probability space, not log space
///
/// Both methods return a **linear probability**, floored. The row takes one logarithm per
/// observation per genotype, after the mixture — putting the log inside would mean taking it
/// per allele instead, and spec §2.1's junk term needs a logarithm around a *sum* over
/// alleles, which a per-allele logarithm cannot express.
///
/// # The scratch is the model's own
///
/// An implementation that needs working memory declares its shape; the row owns one and
/// hands it back on every call, so nothing allocates per observation per candidate. A model
/// that needs none says `type Scratch = ()`.
pub trait SsrEmissionModel {
    /// Working memory this model reuses across calls. `()` for a model that needs none.
    type Scratch: Default;

    /// The probability that one copy of `candidate` produced this whole observed sequence.
    ///
    /// `observation` is the read's own bases over the **tract**, which is what the STR
    /// generator slices — the same footing as [`SsrCandidate::bases`], so the two are
    /// comparable without either being re-measured.
    fn emission(
        &self,
        observation: &[u8],
        candidate: &SsrCandidate<'_>,
        context: &SsrScoringContext<'_>,
        scratch: &mut Self::Scratch,
    ) -> f64;

    /// The probability that one copy of `candidate` produced a read that showed **at least**
    /// this much and then ran out — `P(length ≥ what was witnessed | candidate)` times the
    /// letter match over what was witnessed (spec §5.2).
    ///
    /// **Not a shorter complete observation.** A read that ran out inside the tract has not
    /// shown a shorter allele; it has shown a lower bound, and scoring it as though it pinned
    /// a length would score it as evidence *for* a short allele — the trap spec §5.1 names.
    ///
    /// **What a lower bound is safe about, and what it is not.** It never scores *below* the
    /// complete read of the same bases, so it cannot be mistaken for evidence of a short
    /// allele; that is the property §5.1 needs, and it is tested. It is **not** always the
    /// *less discriminating* of the two: where one candidate is shorter than the stretch the
    /// read got through and the other is longer, the censored read separates them **further**
    /// — it needs no stutter at all under the longer candidate, so it collects the same-length
    /// share, while the complete read needs a one-repeat change under each and the two differ
    /// only by the fitted direction ratio. Measured at 5.661 nats against 1.586 by
    /// `a_censored_read_out_discriminates_a_complete_one_where_the_candidates_straddle_it`.
    /// That is real information — a lower bound rules out everything below it — and it is the
    /// evidence §5.1 turned these reads on to collect.
    ///
    /// *(Spec §5.2 and §12's thirteenth test claimed the stronger property until 2026-08-25,
    /// when this measurement was put to the owner and both were corrected. §5.2's correction
    /// box is the place to read the argument.)*
    ///
    /// **`witnessed_positions` says which flank the read held**, and only an STR bundle reads it
    /// (`str_bundle_loci.md` §7.6): there the read's bases are the start of the carrier's hull
    /// for a read that held the left flank and its end for one that held the right, and the two
    /// are compared differently. A segment reads every partial's bases as the start of the
    /// tract, whichever flank it held, as it always has.
    fn censored_emission(
        &self,
        witnessed_prefix: &[u8],
        witnessed_positions: &WitnessedLocusPositions,
        candidate: &SsrCandidate<'_>,
        context: &SsrScoringContext<'_>,
        scratch: &mut Self::Scratch,
    ) -> f64;
}

/// **Model A** — the default STR emission: how likely a length change is, times how likely
/// the letters are.
///
/// # Two factors, each unable to answer the other's question
///
/// The **stutter factor** decides *how long* the read is: it stretches or trims the candidate
/// to the read's own length, using the distribution
/// ([`StutterModel`]) the pre-pass fitted for this candidate's stratum. The **substitution
/// factor** then decides *which letters*, over two sequences that are the same length by
/// construction — so there is nothing left for it to insert or delete.
///
/// **Keeping them apart is what makes their two parameters separately estimable.** If the
/// substitution factor could also delete bases inside the tract, a read one repeat short
/// would have two explanations — a polymerase slip, or a sequencing deletion — and nothing in
/// the data would choose between them: raising the slippage rate and lowering the indel rate
/// would describe the reads exactly as well as the reverse (spec §4.3).
///
/// # The substitution factor is composed, never re-implemented
///
/// It is [`FlatEmission`] under the context's fitted per-stratum rate — `1 − ε` per matching
/// base and `ε/3` per mismatching one, the alignment module's own comparison
/// (spec §4.3, §7, §9). This model resizes and calls; it does not score letters itself.
///
/// # Where a slip can land
///
/// In a pure tract, adding a repeat anywhere gives the same sequence and there is nothing to
/// choose. In an **interrupted** tract the placements give genuinely different sequences, so
/// the model enumerates them and averages with equal weight.
///
/// **ng enumerates for whole-repeat slips only, and resizes a part-repeat change at the
/// tract's end in a single placement.** That is production's split, stated here rather than
/// left for a coder to guess, and it is a simplification of the same class as the fixed
/// part-repeat share (spec §4.2).
///
/// **A slip that cannot be reached at all scores zero** — contracting away more repeats than
/// a run holds. The mass that costs is not assumed away: it is what
/// [`SsrScoringContext::unreachable_mass`] reports.
#[derive(Debug, Default, Clone, Copy)]
pub struct StutterSubstitutionEmission;

/// Working memory Model A reuses across calls: the candidate's run structure, the placements a
/// whole-repeat slip could land in, and the buffer a part-repeat resize is rendered into.
///
/// Held by the row and handed back on every call, so the buffers themselves are not
/// reallocated per observation per candidate.
///
/// **Two allocations do survive that, and are named here because the sentence above used to
/// claim more than it delivers.** `segment_tract_into` builds a fresh `Vec<u8>` for every
/// interruption it finds, and clearing a `Vec<Vec<u8>>` drops each inner buffer, so
/// `placements` keeps only its spine between calls and `render_tract` allocates one buffer per
/// placement. Both are **per call** rather than per length change — the censored term's exact
/// sum segments once and then walks the support — and closing them needs `render_tract` to
/// write into a reused buffer, which is a change of its own.
#[derive(Debug, Default)]
pub struct StutterSubstitutionScratch {
    /// The candidate split into its runs and interruptions — **a property of the candidate,
    /// so it survives a whole walk over the length changes** rather than being rebuilt for
    /// each one.
    segments: Vec<TractSegment>,
    /// One entry per placement-distinct realisation of the slip.
    placements: Vec<Vec<u8>>,
    /// The candidate resized to the read's length, for the single-placement branch.
    resized: Vec<u8>,
    /// What an STR bundle candidate's scoring reuses ([`bundle_emission`]).
    bundle: BundleScoringScratch,
}

/// Working memory for scoring one read against one STR bundle candidate: each run split into its
/// whole copies and interruptions, each run's placements of the slip being priced, the slip set
/// being visited, and the candidate rendered with it.
#[derive(Debug, Default)]
struct BundleScoringScratch {
    /// Per run, its stretch of the candidate split as [`segment_tract_into`] splits a tract.
    run_segments: Vec<Vec<TractSegment>>,
    /// Per run, the distinct renderings of its stretch under the slip being priced.
    run_placements: Vec<Vec<Vec<u8>>>,
    /// The candidate rendered with one placement chosen in every run.
    rendered: Vec<u8>,
    /// One run's stretch resized by a part-repeat change.
    resized_stretch: Vec<u8>,
    /// Per run, which of its placements the rendering being scored took.
    chosen_placement: Vec<usize>,
    /// The slip sets reaching the read's length change, `runs` entries each, one after another.
    slip_sets: Vec<i64>,
    /// The slip set being built while the sets are enumerated.
    slip_set: Vec<i64>,
    /// What scoring a read that ended inside the hull reuses ([`bundle_censored_emission`]).
    censored: CensoredScratch,
}

/// **How many laid-out alleles a row keeps between reads** ([`CensoredScratch::layouts`]): enough
/// for every candidate of a typical bundle from both flanks, and for a sample whose reads come
/// from a few libraries. A row that needs more lays some out again; nothing else changes.
const KEPT_CENSORED_LAYOUTS: usize = 16;

/// Working memory for the reads that ended inside an STR bundle's hull
/// ([`bundle_censored_emission`]).
///
/// **The layout of an allele does not depend on the read**, only on the allele, its runs, the
/// flank the read held and the run models it is priced under — so it is built once and kept for
/// every read the row scores against that allele from that flank, rather than once a read
/// (str_bundles_partial_reads review, 2026-10-12: rebuilding it was a quarter of the partial
/// reads' cost). A kept layout is found by comparing those inputs themselves, never an address,
/// so a layout can only be reused for an allele it was built from.
#[derive(Debug, Default)]
struct CensoredScratch {
    /// The layouts kept, at most [`KEPT_CENSORED_LAYOUTS`].
    layouts: Vec<CensoredLayout>,
    /// The kept layout replaced next when a new one is needed and every slot is taken.
    next_replaced: usize,
    /// The read's bases, from the flank it held.
    read: Vec<u8>,
    /// `memo[run × read length + offset]`: the read's probability from that run on, given the
    /// read's bases before `offset` are spent; then `memo[(runs + piece) × read length + offset]`,
    /// the same from a fixed piece on. `NaN` until computed.
    memo: Vec<f64>,
    /// A run's stretch split into whole copies and interruptions, while a layout is built.
    segments: Vec<TractSegment>,
    /// A run's renderings under one slip, while a layout is built.
    placements: Vec<Vec<u8>>,
    /// The read's comparisons, with the sums of agreeing bases kept.
    matched: MatchedPrefixLetters,
}

/// **What a layout was built from** — compared in full before a kept layout is reused.
#[derive(Debug, Default)]
struct CensoredLayoutInputs {
    held: Option<HeldFlank>,
    allele: Vec<u8>,
    flank_beyond_the_hull: Vec<u8>,
    runs: Vec<ReferenceRun>,
    runs_in_allele: Vec<crate::region_typing::run_parse::RunInAllele>,
    run_models: Vec<StutterModel>,
}

impl CensoredLayoutInputs {
    fn are(
        &self,
        held: HeldFlank,
        allele: &[u8],
        flank_beyond_the_hull: &[u8],
        bundle: BundleCandidate<'_>,
        context: &SsrScoringContext<'_>,
    ) -> bool {
        self.held == Some(held)
            && self.allele == allele
            && self.flank_beyond_the_hull == flank_beyond_the_hull
            && self.runs == bundle.runs
            && self.runs_in_allele == bundle.parse.runs()
            && self.run_models == context.run_stutter
    }

    fn set(
        &mut self,
        held: HeldFlank,
        allele: &[u8],
        flank_beyond_the_hull: &[u8],
        bundle: BundleCandidate<'_>,
        context: &SsrScoringContext<'_>,
    ) {
        self.held = Some(held);
        self.allele.clear();
        self.allele.extend_from_slice(allele);
        self.flank_beyond_the_hull.clear();
        self.flank_beyond_the_hull
            .extend_from_slice(flank_beyond_the_hull);
        self.runs.clear();
        self.runs.extend_from_slice(bundle.runs);
        self.runs_in_allele.clear();
        self.runs_in_allele.extend_from_slice(bundle.parse.runs());
        self.run_models.clear();
        self.run_models.extend_from_slice(context.run_stutter);
    }
}

/// **An STR bundle allele laid out from the flank a partial read held** — the pieces a carrier
/// shows in the order the read meets them, so the read is compared from its first base whichever
/// flank it held. From the right flank every piece is reversed, and so is the read. It also holds
/// the part-repeat renderings, each followed (or, from the right, preceded) by the flank.
///
/// Built by [`CensoredLayout::lay_out`], which clears and refills every buffer.
#[derive(Debug, Default)]
struct CensoredLayout {
    /// What it was built from.
    inputs: CensoredLayoutInputs,
    /// The fixed pieces' bases, one after another; the last piece carries the flank beyond the
    /// hull on the side the read ran towards.
    fixed_bases: Vec<u8>,
    /// Each fixed piece's bases in `fixed_bases`, in reading order — one more than the runs.
    fixed: Vec<std::ops::Range<usize>>,
    /// Every run's renderings' bases, one after another.
    rendering_bases: Vec<u8>,
    /// Each rendering's bases in `rendering_bases`.
    renderings: Vec<std::ops::Range<usize>>,
    /// Each slip a run can make: its probability and its renderings.
    slips: Vec<CensoredSlip>,
    /// Each run's slips in `slips`, in reading order.
    slips_of_run: Vec<std::ops::Range<usize>>,
    /// `unreached_mass[i]`: the product, over the runs from reading position `i` on, of each
    /// run's total probability over its slips — what those runs contribute when the read ended
    /// before reaching them. One more entry than the runs; the last is one.
    unreached_mass: Vec<f64>,
    /// Every length change no single run's slip makes, in increasing order, with its places.
    part_repeat_changes: Vec<PartRepeatChange>,
    /// Every place of those changes with a non-zero price, in run order within a change.
    part_repeat_places: Vec<PartRepeatPlace>,
    /// The carriers of those places, one after another, in the read's own orientation: the
    /// rendering then the flank from the left, the flank then the rendering from the right.
    part_repeat_carrier_bases: Vec<u8>,
}

/// One slip a run can make, as [`CensoredLayout`] holds it.
#[derive(Debug, Clone)]
struct CensoredSlip {
    /// The run's probability of this slip.
    probability: f64,
    /// Its renderings in [`CensoredLayout::renderings`], each weighted equally.
    renderings: std::ops::Range<usize>,
    /// The length of every rendering.
    length: usize,
}

/// One length change no single run's slip makes, as [`CensoredLayout`] holds it.
#[derive(Debug, Clone)]
struct PartRepeatChange {
    /// How many runs it can be a sequencing indel in — the average's denominator, counting a
    /// place whose price is zero.
    places: usize,
    /// Its places with a non-zero price, in [`CensoredLayout::part_repeat_places`].
    priced: std::ops::Range<usize>,
}

/// One run a part-repeat change can be a sequencing indel in, as [`CensoredLayout`] holds it.
#[derive(Debug, Clone)]
struct PartRepeatPlace {
    /// That run's price for the change, times the other runs' chance of not slipping.
    probability: f64,
    /// The carrier's bases in [`CensoredLayout::part_repeat_carrier_bases`].
    carrier: std::ops::Range<usize>,
}

impl SsrEmissionModel for StutterSubstitutionEmission {
    type Scratch = StutterSubstitutionScratch;

    fn emission(
        &self,
        observation: &[u8],
        candidate: &SsrCandidate<'_>,
        context: &SsrScoringContext<'_>,
        scratch: &mut Self::Scratch,
    ) -> f64 {
        if let Some(bundle) = candidate.bundle {
            return bundle_emission(observation, candidate.bases, bundle, context, scratch);
        }
        let period = period_of(context.motif);
        let bp_diff = observation.len() as i64 - candidate.bases.len() as i64;

        // How likely a length change of this size is at all. Zero past the cutoffs, and a
        // read that lands there falls to the row's outlier term rather than being explained
        // away as an implausibly large slip.
        let length_probability = context.stutter.probability(bp_diff, period);
        if length_probability <= 0.0 {
            return 0.0;
        }

        let motif = context.motif.as_bytes();
        segment_tract_into(candidate.bases, motif, &mut scratch.segments);
        let Some(letters) = letters_over(
            observation,
            candidate.bases,
            motif,
            bp_diff,
            context,
            scratch,
        ) else {
            // Unreachable from this candidate — more repeats contracted than any run holds.
            // `SsrScoringContext::unreachable_mass` is what accounts for it.
            return 0.0;
        };

        length_probability * letters
    }

    fn censored_emission(
        &self,
        witnessed_prefix: &[u8],
        witnessed_positions: &WitnessedLocusPositions,
        candidate: &SsrCandidate<'_>,
        context: &SsrScoringContext<'_>,
        scratch: &mut Self::Scratch,
    ) -> f64 {
        if let Some(bundle) = candidate.bundle {
            return bundle_censored_emission(
                witnessed_prefix,
                witnessed_positions,
                candidate,
                bundle,
                context,
                scratch,
            );
        }
        let period = period_of(context.motif);
        let motif = context.motif.as_bytes();
        let witnessed_bases = witnessed_prefix.len();

        // The read got through `witnessed_bases` bases of tract, so the tract is at least that
        // long — which for this candidate means a length change of at least this much. It is
        // negative whenever the candidate is already long enough on its own, and that is the
        // ordinary case rather than an edge one.
        let smallest_bp_diff = witnessed_bases as i64 - candidate.bases.len() as i64;

        if is_a_pure_tract(candidate.bases, motif) {
            // **The factorised form, and on a pure tract it is exact.** Every stretching of a
            // pure tract agrees on its first `witnessed_bases` bases — stretching appends or
            // trims at the end — so the letters come out of the sum and what is left is the
            // length factor (spec §5.2).
            let length_probability_at_least =
                context.stutter.probability_at_least_this_much_longer(
                    smallest_bp_diff,
                    period,
                    candidate.repeat_count,
                );
            // An early exit rather than a correction: the letter factor is finite, so falling
            // through would multiply out to this same zero. It saves a resize.
            if length_probability_at_least <= 0.0 {
                return 0.0;
            }
            resize_at_the_end(
                candidate.bases,
                motif,
                smallest_bp_diff,
                &mut scratch.resized,
            );
            return length_probability_at_least
                * substitution_probability(witnessed_prefix, &scratch.resized, context);
        }

        // **The exact sum, on an interrupted tract.** Here the stretchings give genuinely
        // different first-`witnessed_bases` bases, so the letters cannot come out of the sum:
        // an interruption sits at a different offset in each one. Spec §5.2 bounds what the
        // factorisation would cost here at `log(3(1−ε)/ε)` per distinguishing base — 6.4 nats
        // at an error rate of 1 in 200 — and says to pay for the sum instead.
        //
        // **What that costs is one letter comparison per placement per length change**, not
        // one per length change: the support is at most 41 changes wide and an interrupted
        // tract offers one placement per run, so a two-run tract pays at most 82 comparisons
        // of the tract's length. It reaches only candidates that are interrupted at all — a
        // pure one took the branch above.
        //
        // **The runs are found once**, before the walk: they are a property of the candidate,
        // and the candidate does not move inside the loop.
        segment_tract_into(candidate.bases, motif, &mut scratch.segments);
        let mut total = 0.0;
        for (bp_diff, length_probability) in context
            .stutter
            .reachable_length_changes(period, candidate.repeat_count)
        {
            // Only the stretchings that leave the tract at least as long as what the read
            // showed. Everything else is a tract the read could not have run out inside — and
            // this filter is also what makes every prefix below well defined, because a
            // stretching admitted here is by construction at least `witnessed_bases` long.
            if bp_diff < smallest_bp_diff || length_probability <= 0.0 {
                continue;
            }
            let Some(letters) = letters_over(
                witnessed_prefix,
                candidate.bases,
                motif,
                bp_diff,
                context,
                scratch,
            ) else {
                // Unreachable from this candidate's runs. A slip lands in *one* run, and the
                // distribution's own reachability rule sees only the total, so the support can
                // offer a contraction no single run can absorb — the recorded open question.
                // Skipping is an early exit rather than a correction: an empty placement list
                // contributes zero to this sum either way.
                continue;
            };
            total += length_probability * letters;
        }
        total
    }
}

/// **How probable the letters are, given this candidate stretched by `bp_diff`** — compared
/// over `observed`'s own length, and `None` when no run of the candidate can absorb the
/// change.
///
/// # One function, two questions
///
/// **The comparison stops at `observed.len()`, and that is what lets one function serve both
/// [`SsrEmissionModel::emission`] and [`SsrEmissionModel::censored_emission`].** A complete
/// read is exactly as long as the stretched candidate — `bp_diff` was derived from it — so
/// truncating to its length is a no-op there. A censored read is shorter, and that truncation
/// is the whole difference between the two callers. Written once because the alternative is
/// two copies of the whole-versus-part dispatch, the equal-weight average and the unreachable
/// rule, and the next change to where a slip may land would have to reach both.
///
/// # The caller owes the segmentation
///
/// `scratch.segments` must already hold this candidate's runs ([`segment_tract_into`]). That
/// is what lets a caller walk many length changes against one candidate without re-splitting a
/// tract that does not move: the censored term's exact sum pays for it once and then walks the
/// support.
///
/// # Where a slip can land
///
/// A whole-repeat change enumerates one placement per run it could have landed in and averages
/// them with equal weight; a part-repeat change is resized at the tract's end in a single
/// placement (spec §4.2).
fn letters_over(
    observed: &[u8],
    candidate: &[u8],
    motif: &[u8],
    bp_diff: i64,
    context: &SsrScoringContext<'_>,
    scratch: &mut StutterSubstitutionScratch,
) -> Option<f64> {
    let period_bases = motif.len() as i64;

    if bp_diff == 0 {
        // The one placement is the candidate itself, so there is nothing to render.
        return Some(substitution_probability(
            observed,
            prefix_of_stretching(candidate, observed.len())?,
            context,
        ));
    }

    if bp_diff % period_bases == 0 {
        // A whole-repeat change: every run the slip could have landed in gives its own
        // sequence, and they are averaged with equal weight.
        place_slip_in_each_run(
            &scratch.segments,
            motif,
            bp_diff / period_bases,
            &mut scratch.placements,
        );
        if scratch.placements.is_empty() {
            return None;
        }
        let each = 1.0 / scratch.placements.len() as f64;
        Some(
            scratch
                .placements
                .iter()
                .map(
                    |placement| match prefix_of_stretching(placement, observed.len()) {
                        Some(prefix) => each * substitution_probability(observed, prefix, context),
                        None => 0.0,
                    },
                )
                .sum(),
        )
    } else {
        // A part-repeat change: one placement, resized at the tract's end.
        resize_at_the_end(candidate, motif, bp_diff, &mut scratch.resized);
        Some(
            match prefix_of_stretching(&scratch.resized, observed.len()) {
                Some(prefix) => substitution_probability(observed, prefix, context),
                None => 0.0,
            },
        )
    }
}

/// The first `length` bases of a stretching, or `None` if the stretching is shorter than that.
///
/// **It cannot be `None` at either call site, and it is fallible anyway.** The guarantee spans
/// three statements and two functions — `censored_emission`'s floor filter admits only
/// stretchings at least as long as what the read showed, and a complete read is exactly as
/// long as its own stretching — which is more derivation than a bare index should rest on in a
/// path that runs once per observation per candidate. Indexing would abort a calling run;
/// this scores nothing, which is what [`substitution_probability`] already does when handed
/// two lengths that disagree, so the two branches of the sum now fail the same way.
#[inline]
fn prefix_of_stretching(stretching: &[u8], length: usize) -> Option<&[u8]> {
    let prefix = stretching.get(..length);
    debug_assert!(
        prefix.is_some(),
        "the floor filter admits only stretchings at least as long as what the read showed"
    );
    prefix
}

/// Whether this candidate's tract is whole copies of the motif and nothing else.
///
/// **It decides which of spec §5.2's two forms scores a censored read**, and that is the only
/// thing it is for: on a pure tract every stretching agrees on the witnessed prefix, so the
/// letters factor out exactly; on an interrupted one they do not, because an interruption
/// sits at a different offset in each stretching.
///
/// **It must agree with [`segment_tract`], which is what actually places a slip.** That is
/// asserted rather than argued — `a_tract_is_pure_exactly_when_the_segmenter_finds_one_run`
/// checks the two against each other over every tract of eight bases or fewer on a two-letter
/// alphabet — because `segment_tract`'s own documentation calls its greedy left-to-right phase
/// *a choice*, and a later change to that choice must not be able to leave the two classifying
/// different tracts.
///
/// **The length check is not redundant with the chunk walk**, and it is the conjunct a reader
/// is most likely to think is. `chunks_exact` drops the remainder on the floor, so without it
/// `CAGCAGTT` under motif `CAG` comes back pure — and it is not, because its stretchings are
/// built from its runs rather than by continuing the tiling off the end. Measured on that
/// tract, a censored read of `CAGCAGCAG` scores 3.369e-3 down the interrupted route against
/// 3.939e-10 down the pure one.
fn is_a_pure_tract(candidate: &[u8], motif: &[u8]) -> bool {
    candidate.len().is_multiple_of(motif.len())
        && candidate
            .chunks_exact(motif.len())
            .all(|chunk| chunk == motif)
}

/// `P(the letters differ this way)` for two sequences the stutter factor has already made the
/// same length — [`FlatEmission`] under the context's fitted rate, in probability space.
///
/// **Unequal lengths cannot arrive here**: every branch above resizes first, which is the
/// separation spec §4.3 requires. A mismatch in length would mean the resize was wrong, so it
/// is an assertion rather than a silent zero.
fn substitution_probability(
    observation: &[u8],
    resized_candidate: &[u8],
    context: &SsrScoringContext<'_>,
) -> f64 {
    letters_under(substitution_scores(context), observation, resized_candidate)
}

/// **The two per-base scores of [`substitution_probability`]** under the context's rate — two
/// logarithms, so a caller comparing one read many times looks them up once and calls
/// [`letters_under`] (a read that ended inside an STR bundle, [`bundle_censored_emission`]).
fn substitution_scores(context: &SsrScoringContext<'_>) -> BaseScores {
    let flat = FlatEmission::try_new(context.substitution_rate.get())
        .expect("an ErrorRate is a probability by construction");
    // A flat model's two scores do not vary with quality, so they are resolved once and the
    // loop is one comparison per base.
    // A flat model ignores the quality entirely, so any value resolves the same two
    // scores; zero is the one that says "this is not read from the base".
    flat.scores_for(BaseQual(0))
}

/// [`substitution_probability`] under scores already looked up ([`substitution_scores`]).
fn letters_under(scores: BaseScores, observation: &[u8], resized_candidate: &[u8]) -> f64 {
    debug_assert_eq!(
        observation.len(),
        resized_candidate.len(),
        "the stutter factor equalises the lengths before the letters are compared"
    );
    if observation.len() != resized_candidate.len() {
        return 0.0;
    }
    let total_ln: f64 = observation
        .iter()
        .zip(resized_candidate)
        .map(|(&read_base, &candidate_base)| scores.pick(read_base, candidate_base))
        .sum();
    float::exp(total_ln)
}

/// **[`letters_under`] for many comparisons of one read under one pair of scores, without
/// re-adding the bases that agree.** A comparison's log-probability is a sum taken base by base
/// from the first; up to its first mismatch every term is the match score, so that much of the sum
/// is the sum of `i` match scores — the same number for every comparison that agrees on its first
/// `i` bases. Those sums are kept, with their exponentials, and a comparison adds only from its first
/// mismatch on. **The result is bit for bit [`letters_under`]'s**: the same additions in the same
/// order, from the same starting value (`a_matched_prefix_gives_letters_under_bit_for_bit`).
///
/// It is what a read that ended inside an STR bundle is compared with ([`bundle_censored_emission`]):
/// there most comparisons are a run's rendering against the same run in the read, which agree to
/// their end (str_bundles_partial_reads review, 2026-10-12: the base-by-base sum was most of what was
/// left of those reads' cost once their layouts were kept).
#[derive(Debug, Default)]
struct MatchedPrefixLetters {
    /// The scores the sums are of; `None` before the first use.
    scores: Option<BaseScores>,
    /// `sums[i]`: the sum of `i` match scores, taken as [`letters_under`] takes it.
    sums: Vec<f64>,
    /// `letters[i]`: `exp(sums[i])`.
    letters: Vec<f64>,
}

impl MatchedPrefixLetters {
    /// Make ready for comparisons of up to `longest` bases under `scores`.
    fn prepare(&mut self, scores: BaseScores, longest: usize) {
        let same_scores = self.scores.is_some_and(|kept| {
            kept.match_ln.to_bits() == scores.match_ln.to_bits()
                && kept.mismatch_ln.to_bits() == scores.mismatch_ln.to_bits()
        });
        if !same_scores {
            self.scores = Some(scores);
            self.sums.clear();
            self.letters.clear();
        }
        if self.sums.is_empty() {
            // The empty sum, whatever starting value `Sum` folds from.
            let empty: f64 = std::iter::empty::<f64>().sum();
            self.sums.push(empty);
            self.letters.push(float::exp(empty));
        }
        while self.sums.len() <= longest {
            let next = self.sums[self.sums.len() - 1] + scores.match_ln;
            self.sums.push(next);
            self.letters.push(float::exp(next));
        }
    }

    /// [`letters_under`]`(scores, observation, resized_candidate)`, for comparisons no longer than
    /// [`Self::prepare`] was given.
    fn letters(&self, observation: &[u8], resized_candidate: &[u8]) -> f64 {
        debug_assert_eq!(
            observation.len(),
            resized_candidate.len(),
            "the stutter factor equalises the lengths before the letters are compared"
        );
        if observation.len() != resized_candidate.len() {
            return 0.0;
        }
        let scores = self.scores.expect("prepared before the first comparison");
        let agreed = agreeing_prefix_len(observation, resized_candidate);
        if agreed == observation.len() {
            return self.letters[agreed];
        }
        let total_ln = observation[agreed..]
            .iter()
            .zip(&resized_candidate[agreed..])
            .map(|(&read_base, &candidate_base)| scores.pick(read_base, candidate_base))
            .fold(self.sums[agreed], |sum, score| sum + score);
        float::exp(total_ln)
    }
}

/// How many bases two equally long sequences agree on before their first difference.
fn agreeing_prefix_len(one: &[u8], other: &[u8]) -> usize {
    let mut agreed = 0;
    for (one, other) in one.as_chunks::<8>().0.iter().zip(other.as_chunks::<8>().0) {
        let differ = u64::from_le_bytes(*one) ^ u64::from_le_bytes(*other);
        if differ != 0 {
            return agreed + (differ.trailing_zeros() / 8) as usize;
        }
        agreed += 8;
    }
    agreed
        + one[agreed..]
            .iter()
            .zip(&other[agreed..])
            .take_while(|(one, other)| one == other)
            .count()
}

/// Every placement-distinct realisation of `candidate` with `repeats` whole repeats added to
/// one of its runs, written into `out` (cleared first).
///
/// A pure tract is one run, so it yields one sequence. An interrupted tract yields one per run
/// the slip could land in, **deduplicated** — two runs of the same length in a symmetric tract
/// give the same bytes, and counting that twice would weight it double.
///
/// A run that cannot absorb the change contributes nothing, and a candidate where **no** run
/// can yields an empty list, which the caller reads as unreachable. Two limits make a run
/// unable to: it cannot go below zero repeats, and **the tract as a whole must keep at least
/// one** — the same rule [`StutterModel::unreachable_mass`] applies from the distribution's
/// side, so the mass the report calls unplaced is the mass the scoring actually declines to
/// place.
///
/// **On a pure tract the two agree exactly. On an interrupted one they can differ**, and in a
/// stated direction: a slip lands in *one* run, so a tract of two runs of two repeats cannot
/// lose three even though it holds four — while `unreachable_mass` sees only the total and
/// counts that contraction as reachable. **The report therefore understates the loss on
/// interrupted tracts.** Closing that needs the run structure to reach the distribution,
/// which takes only a period and a repeat count today; it is recorded as an open question
/// rather than papered over.
/// **Test-only since the scoring paths began reusing a segmentation.** `place_slip_in_each_run`
/// is what `emission` and `censored_emission` call; this wrapper exists for the tests and the
/// comparison model, which score one length change at a time and have no buffer to keep.
#[cfg(test)]
fn enumerate_placements(candidate: &[u8], motif: &[u8], repeats: i64, out: &mut Vec<Vec<u8>>) {
    if repeats == 0 {
        out.clear();
        out.push(candidate.to_vec());
        return;
    }
    place_slip_in_each_run(&segment_tract(candidate, motif), motif, repeats, out);
}

/// The same enumeration over a tract whose runs have already been found — what a caller uses
/// when it walks many length changes against **one** candidate, so the split is paid for once
/// rather than once a change.
///
/// `repeats` is non-zero: the identity has no slip to place, and [`enumerate_placements`]
/// answers it without splitting anything.
fn place_slip_in_each_run(
    segments: &[TractSegment],
    motif: &[u8],
    repeats: i64,
    out: &mut Vec<Vec<u8>>,
) {
    debug_assert_ne!(repeats, 0, "the identity has no slip to place");
    out.clear();
    let run_count = segments
        .iter()
        .filter(|segment| matches!(segment, TractSegment::Run(_)))
        .count();

    let held_in_total: usize = segments
        .iter()
        .filter_map(|segment| match segment {
            TractSegment::Run(repeats) => Some(*repeats),
            TractSegment::Interruption(_) => None,
        })
        .sum();

    for target in 0..run_count {
        let held = segments
            .iter()
            .filter_map(|segment| match segment {
                TractSegment::Run(repeats) => Some(*repeats),
                TractSegment::Interruption(_) => None,
            })
            .nth(target)
            .expect("the run count bounds the index");
        let Some(resized) = (held as i64).checked_add(repeats).filter(|held| *held >= 0) else {
            continue;
        };
        // **The tract must keep a repeat**, which is the same rule
        // [`StutterModel::unreachable_mass`] applies from the distribution's side: a tract
        // holding none is not this locus any more. Without it the two halves disagree by one
        // step on a pure tract — the report would call a total contraction unreachable while
        // the scoring placed mass on it.
        if held_in_total + resized as usize - held == 0 {
            continue;
        }
        let rendered = render_tract(segments, motif, target, resized as usize);
        if !out.contains(&rendered) {
            out.push(rendered);
        }
    }
}

/// One stretch of a tract: a run of whole motif copies, or the bases between two runs.
#[derive(Debug, Clone, PartialEq, Eq)]
enum TractSegment {
    /// How many whole copies of the motif this run holds.
    Run(usize),
    /// Bases that are not a whole copy — an interruption, or a partial copy at either end.
    Interruption(Vec<u8>),
}

/// Split a tract into runs of the motif and the interruptions between them, matching the
/// motif greedily from the left.
///
/// **Left-to-right and greedy, which is a choice about interrupted tracts.** A different
/// phase would split some tracts differently; this is production's rule, and the two must
/// agree while the models are meant to.
#[cfg(test)]
fn segment_tract(candidate: &[u8], motif: &[u8]) -> Vec<TractSegment> {
    let mut segments = Vec::new();
    segment_tract_into(candidate, motif, &mut segments);
    segments
}

/// The same split, written into a buffer the caller owns and reuses (cleared first).
///
/// The scoring paths take this one. A candidate's runs do not change while the model walks the
/// length changes that candidate could show, so splitting once a call rather than once a
/// change is the difference between one split and up to forty-one of them.
fn segment_tract_into(candidate: &[u8], motif: &[u8], segments: &mut Vec<TractSegment>) {
    segments.clear();
    let period = motif.len();
    let mut run = 0usize;
    let mut interruption: Vec<u8> = Vec::new();
    let mut at = 0usize;

    while at < candidate.len() {
        if at + period <= candidate.len() && &candidate[at..at + period] == motif {
            if !interruption.is_empty() {
                segments.push(TractSegment::Interruption(std::mem::take(
                    &mut interruption,
                )));
            }
            run += 1;
            at += period;
        } else {
            if run > 0 {
                segments.push(TractSegment::Run(run));
                run = 0;
            }
            interruption.push(candidate[at]);
            at += 1;
        }
    }
    if run > 0 {
        segments.push(TractSegment::Run(run));
    }
    if !interruption.is_empty() {
        segments.push(TractSegment::Interruption(interruption));
    }
}

/// Render the segments back to bytes, with the run at `target` holding `resized` repeats
/// instead of its own.
fn render_tract(segments: &[TractSegment], motif: &[u8], target: usize, resized: usize) -> Vec<u8> {
    let mut out = Vec::new();
    let mut seen = 0usize;
    for segment in segments {
        match segment {
            TractSegment::Run(held) => {
                let repeats = if seen == target { resized } else { *held };
                for _ in 0..repeats {
                    out.extend_from_slice(motif);
                }
                seen += 1;
            }
            TractSegment::Interruption(bases) => out.extend_from_slice(bases),
        }
    }
    out
}

/// Resize a candidate by `bp_diff` bases **at the tract's end**, in one placement — the
/// part-repeat branch.
///
/// Lengthening continues the motif's tiling from the phase the tract ended on, so the added
/// bases are the ones a slipped polymerase would have written. Shortening drops bases from the
/// end. Production's rule, and a simplification: a part-repeat change could have landed
/// anywhere, and enumerating those placements is the follow-up spec §4.2 records.
fn resize_at_the_end(candidate: &[u8], motif: &[u8], bp_diff: i64, out: &mut Vec<u8>) {
    out.clear();
    let period = motif.len();
    if bp_diff >= 0 {
        out.extend_from_slice(candidate);
        for step in 0..bp_diff as usize {
            out.push(motif[(candidate.len() + step) % period]);
        }
    } else {
        let keep = candidate
            .len()
            .saturating_sub(bp_diff.unsigned_abs() as usize);
        out.extend_from_slice(&candidate[..keep]);
    }
}

// ── STR bundles: stutter in every run (`str_bundle_loci.md` §7.3) ──

/// **The probability that one copy of an STR bundle allele produced this whole read** — Model A
/// with a slip allowed in every run of the bundle at once (`str_bundle_loci.md` §7.3).
///
/// # The sum over slip sets
///
/// Each run slips independently, by `kᵣ` whole units of its own motif, with the probability its
/// own stutter model gives — the model of its period at its count in this allele, `kᵣ = 0` being
/// the chance it did not slip. A read whose length differs from the allele's by `Δ` bases can come
/// from any **slip set** `(k₁, …, k_R)` whose lengths add up to it, `Σ kᵣ·periodᵣ = Δ`, so the
/// emission is the sum over those sets of the product of the runs' probabilities, times the
/// letters of the read against the allele with those slips in place:
///
/// ```text
/// P(read | allele) = Σ over slip sets   Π_r P_r(kᵣ)  ×  letters(read | allele with the slips)
/// ```
///
/// - **Each run's slip is bounded by the model's own reach**, [`MAX_WHOLE_REPEAT_SLIP`] units
///   either way, and a run must keep a whole copy; a run the allele holds no copy of does not
///   slip.
/// - **Within a run, a slip lands as in an interrupted segment**: at each distinct place in that
///   run's stretch of the allele, weighted equally ([`place_slip_in_each_run`]). Where several
///   runs slip at once, every combination of their places is one rendering, weighted by the
///   product of the runs' equal weights.
/// - **A length change no single run's slip makes is also a sequencing indel** — the part-repeat
///   regime, as a read off the grid is at a segment: an indel in one run, priced by that run's
///   model for a change that is not whole units of its motif and resized at the run's end, while
///   the other runs did not slip. Every run for which `Δ` is not whole units is such a place,
///   weighted equally. It is added to whatever slip sets of several runs reach `Δ` (owner,
///   2026-10-12): at runs of periods 2 and 3 every change is some set's total, and a one-base
///   indel priced only as two runs slipping opposite ways was a quarter of a segment's price for
///   it. A change one run makes alone is that run's slip, as at a segment. A bundle with no run
///   has no such place, and a read of another length scores nothing — the outlier term takes it.
///
/// **A bundle of one run that spans the whole hull scores exactly as the segment over the same
/// bases**: one slip set per whole-unit `Δ`, with the run's own placements, and the part-repeat
/// regime resized at the allele's end (`a_one_run_bundle_scores_as_the_segment_it_is`).
///
/// # Cost
///
/// A slip set is one run's change per run, so two runs give at most `2·MAX_WHOLE_REPEAT_SLIP + 1`
/// sets for a `Δ` and three a few hundred. **From three runs on, a branch whose probability
/// cannot reach [`NEGLIGIBLE_SLIP_SET_PROBABILITY`] is not walked**: four runs would otherwise walk
/// thousands of sets a read, and five (2 of HG002's 2,410 bundles) a hundred thousand, each worth
/// less than a millionth of the row's outlier floor.
fn bundle_emission(
    observation: &[u8],
    candidate: &[u8],
    bundle: BundleCandidate<'_>,
    context: &SsrScoringContext<'_>,
    scratch: &mut StutterSubstitutionScratch,
) -> f64 {
    let runs = bundle.runs;
    let placed = bundle.parse.runs();
    debug_assert_eq!(runs.len(), placed.len(), "one place per reference run");
    debug_assert_eq!(context.run_stutter.len(), runs.len(), "one model per run");
    let bp_diff = observation.len() as i64 - candidate.len() as i64;
    let shift = bundle.leading_flank_bases as usize;
    let stretch_of = |run: usize| -> std::ops::Range<usize> {
        shift + placed[run].start as usize..shift + placed[run].end as usize
    };
    let scratch = &mut scratch.bundle;

    // Each run's stretch, split into whole copies and interruptions once for the whole walk.
    scratch.run_segments.resize_with(runs.len(), Vec::new);
    for (run, reference_run) in runs.iter().enumerate() {
        segment_tract_into(
            &candidate[stretch_of(run)],
            reference_run.motif().as_bytes(),
            &mut scratch.run_segments[run],
        );
    }

    let ranges = run_slip_ranges(bundle, context);
    let prune = runs.len() >= 3;
    scratch.slip_sets.clear();
    scratch.slip_set.clear();
    collect_slip_sets(
        &ranges,
        bp_diff,
        prune,
        1.0,
        &mut scratch.slip_set,
        &mut scratch.slip_sets,
    );

    // **A change no single run's slip makes is also a sequencing indel** (owner, 2026-10-12): at
    // `(CA)n (AAT)m` a one-base change is reached only by the two runs slipping opposite ways,
    // and priced by that alone it was a quarter of what a segment pays for a one-base indel.
    let reached_by_one_run = a_single_run_slip_reaches(&ranges, bp_diff);
    if scratch.slip_sets.is_empty() {
        // Every slip set reaching the change may have been pruned as negligible; a change one
        // run alone makes is then worth less than the threshold, and is no sequencing indel.
        if reached_by_one_run {
            return 0.0;
        }
        return part_repeat_emission(
            candidate,
            bundle,
            context,
            bp_diff,
            (&mut scratch.resized_stretch, &mut scratch.rendered),
            |rendered| substitution_probability(observation, rendered, context),
        );
    }

    let mut total = 0.0;
    let set_len = runs.len().max(1);
    let set_count = if runs.is_empty() {
        1
    } else {
        scratch.slip_sets.len() / set_len
    };
    for set in 0..set_count {
        let mut weight = 1.0;
        for (run, range) in ranges.iter().enumerate() {
            let slip = scratch.slip_sets[set * set_len + run];
            weight *= range.probability[(slip - range.lowest) as usize];
        }
        if weight <= 0.0 {
            continue;
        }
        // Each run's renderings under its slip of this set.
        scratch.run_placements.resize_with(runs.len(), Vec::new);
        let mut reachable = true;
        for (run, reference_run) in runs.iter().enumerate() {
            let slip = scratch.slip_sets[set * set_len + run];
            let placements = &mut scratch.run_placements[run];
            if slip == 0 {
                placements.clear();
                placements.push(candidate[stretch_of(run)].to_vec());
            } else {
                place_slip_in_each_run(
                    &scratch.run_segments[run],
                    reference_run.motif().as_bytes(),
                    slip,
                    placements,
                );
            }
            if placements.is_empty() {
                reachable = false;
                break;
            }
        }
        if !reachable {
            continue;
        }
        let letters = letters_over_every_combination(
            observation,
            candidate,
            &stretch_of,
            runs.len(),
            context,
            scratch,
        );
        total += weight * letters;
    }
    if !reached_by_one_run {
        total += part_repeat_emission(
            candidate,
            bundle,
            context,
            bp_diff,
            (&mut scratch.resized_stretch, &mut scratch.rendered),
            |rendered| substitution_probability(observation, rendered, context),
        );
    }
    total
}

/// **Each run's range of slips and the probability of each**, in hull order — what both of a
/// bundle's emissions walk. A run keeps a whole copy and slips at most
/// [`MAX_WHOLE_REPEAT_SLIP`] units either way; a run the allele holds no copy of cannot slip.
fn run_slip_ranges(
    bundle: BundleCandidate<'_>,
    context: &SsrScoringContext<'_>,
) -> Vec<RunSlipRange> {
    let reach = i64::from(MAX_WHOLE_REPEAT_SLIP);
    bundle
        .runs
        .iter()
        .zip(bundle.parse.runs())
        .zip(context.run_stutter)
        .map(|((reference_run, in_allele), model)| {
            let period = reference_run.period() as i64;
            let copies = i64::from(in_allele.copies);
            let (lowest, highest) = if copies == 0 {
                (0, 0)
            } else {
                ((1 - copies).max(-reach), reach)
            };
            // A run the allele holds no copy of cannot slip, so it keeps its length for
            // certain: a factor of one, not its model's no-slip share (found in review,
            // 2026-10-11 — that share, 0.88 under the shipped model, scaled down every read
            // of an allele that lost a whole run).
            let probability = if copies == 0 {
                vec![1.0]
            } else {
                (lowest..=highest)
                    .map(|slip| model.probability(slip * period, period_of(&reference_run.motif())))
                    .collect()
            };
            RunSlipRange {
                period,
                lowest,
                highest,
                probability,
            }
        })
        .collect()
}

/// **Whether one run, slipping alone, changes the allele's length by `bp_diff`** — whole units of
/// its period, within its range. The identity is reached by every bundle, with no run too.
fn a_single_run_slip_reaches(ranges: &[RunSlipRange], bp_diff: i64) -> bool {
    bp_diff == 0
        || ranges.iter().any(|range| {
            bp_diff % range.period == 0
                && (range.lowest..=range.highest).contains(&(bp_diff / range.period))
        })
}

/// The letters of `observation` against the candidate rendered with every combination of the
/// runs' placements in [`BundleScoringScratch::run_placements`], each combination weighted by the
/// product of the runs' equal weights — the bundle's spelling of [`letters_over`]'s average.
fn letters_over_every_combination(
    observation: &[u8],
    candidate: &[u8],
    stretch_of: &impl Fn(usize) -> std::ops::Range<usize>,
    run_count: usize,
    context: &SsrScoringContext<'_>,
    scratch: &mut BundleScoringScratch,
) -> f64 {
    scratch.chosen_placement.clear();
    scratch.chosen_placement.resize(run_count, 0);
    let each: f64 = scratch.run_placements[..run_count]
        .iter()
        .map(|placements| 1.0 / placements.len() as f64)
        .product();
    let mut letters = 0.0;
    loop {
        scratch.rendered.clear();
        let mut from = 0;
        for run in 0..run_count {
            let stretch = stretch_of(run);
            scratch
                .rendered
                .extend_from_slice(&candidate[from..stretch.start]);
            scratch
                .rendered
                .extend_from_slice(&scratch.run_placements[run][scratch.chosen_placement[run]]);
            from = stretch.end;
        }
        scratch.rendered.extend_from_slice(&candidate[from..]);
        letters += each * substitution_probability(observation, &scratch.rendered, context);

        // The next combination, the last run's placement turning fastest.
        let mut run = run_count;
        loop {
            if run == 0 {
                return letters;
            }
            run -= 1;
            scratch.chosen_placement[run] += 1;
            if scratch.chosen_placement[run] < scratch.run_placements[run].len() {
                break;
            }
            scratch.chosen_placement[run] = 0;
        }
    }
}

/// **The part-repeat regime of a bundle**: a length change no single run's slip makes, priced as
/// a sequencing indel in one run — that run's model for a change that is not whole units of its
/// motif, resized at the run's end, times the other runs' chance of not slipping — averaged with
/// equal weight over the runs for which the change is not whole units. Zero where there is no
/// such run.
///
/// `letters` prices the read against each rendering. A read that ended inside the hull takes the
/// renderings from [`part_repeat_places`] instead, once per allele ([`CensoredLayout::lay_out`]).
fn part_repeat_emission(
    candidate: &[u8],
    bundle: BundleCandidate<'_>,
    context: &SsrScoringContext<'_>,
    bp_diff: i64,
    buffers: (&mut Vec<u8>, &mut Vec<u8>),
    mut letters: impl FnMut(&[u8]) -> f64,
) -> f64 {
    let mut total = 0.0;
    let places = part_repeat_places(
        candidate,
        bundle,
        context,
        bp_diff,
        buffers,
        |probability, rendered| total += probability * letters(rendered),
    );
    if places == 0 {
        return 0.0;
    }
    total / places as f64
}

/// **Every place a part-repeat change of `bp_diff` can be a sequencing indel in** — each run for
/// which the change is not whole units, holding a copy and long enough to take it. Hands `place`
/// each one's price (the run's model for the change times the other runs' chance of not slipping)
/// and the allele rendered with it, in run order, skipping a place whose price is zero; returns
/// how many places there are, those included.
fn part_repeat_places(
    candidate: &[u8],
    bundle: BundleCandidate<'_>,
    context: &SsrScoringContext<'_>,
    bp_diff: i64,
    (resized_stretch, rendered): (&mut Vec<u8>, &mut Vec<u8>),
    mut place: impl FnMut(f64, &[u8]),
) -> usize {
    let runs = bundle.runs;
    let placed = bundle.parse.runs();
    let shift = bundle.leading_flank_bases as usize;
    let mut places = 0usize;
    for (run, reference_run) in runs.iter().enumerate() {
        let period = period_of(&reference_run.motif());
        if bp_diff % i64::from(period.get()) == 0 {
            continue;
        }
        // A run with no copy cannot take a sequencing indel, and a stretch shorter than the
        // contraction cannot be shortened by it: neither is a place (found in review,
        // 2026-10-11 — such a place rendered an allele of another length than the read, and in
        // a release build scored zero while still counting in the average). At a segment the
        // stretch is the whole allele, so neither case arises there.
        let stretch = shift + placed[run].start as usize..shift + placed[run].end as usize;
        if placed[run].copies == 0 || (stretch.len() as i64) + bp_diff < 0 {
            continue;
        }
        places += 1;
        let model = &context.run_stutter[run];
        let mut probability = model.probability(bp_diff, period);
        for (other, (other_model, other_placed)) in
            context.run_stutter.iter().zip(placed).enumerate()
        {
            // A run with no copy cannot slip: certain to keep its length.
            if other != run && other_placed.copies > 0 {
                probability *= other_model.same_length_share();
            }
        }
        if probability <= 0.0 {
            continue;
        }
        resize_at_the_end(
            &candidate[stretch.clone()],
            reference_run.motif().as_bytes(),
            bp_diff,
            resized_stretch,
        );
        rendered.clear();
        rendered.extend_from_slice(&candidate[..stretch.start]);
        rendered.extend_from_slice(resized_stretch);
        rendered.extend_from_slice(&candidate[stretch.end..]);
        place(probability, rendered);
    }
    places
}

/// The range of slips one run can make, and the probability of each, `lowest` first.
struct RunSlipRange {
    /// The run's motif length in bases.
    period: i64,
    /// The largest contraction, as a negative number of units: the run keeps a whole copy.
    lowest: i64,
    /// The largest expansion, in units.
    highest: i64,
    /// `probability[slip − lowest]`: the run's model at that slip, no slip included.
    probability: Vec<f64>,
}

/// **The smallest probability a slip set must be able to reach to be walked**, from three runs on
/// ([`bundle_emission`]). A slip set's emission is at most its probability, and the row adds an
/// outlier floor near `0.2 / 40` to every read (`likelihood/ssr.rs`), so a set below a millionth of
/// a millionth changes no genotype's likelihood in its tenth significant digit.
pub const NEGLIGIBLE_SLIP_SET_PROBABILITY: f64 = 1e-12;

/// Write every slip set of `ranges` whose lengths add up to `bp_diff` into `out`, one after
/// another, `ranges.len()` entries each — runs in hull order. With `prune`, a branch whose
/// probability cannot reach [`NEGLIGIBLE_SLIP_SET_PROBABILITY`] is not walked. With no run, the
/// one empty set reaches a change of zero and nothing else, and is written as one entry, `0`.
fn collect_slip_sets(
    ranges: &[RunSlipRange],
    bp_diff: i64,
    prune: bool,
    probability_so_far: f64,
    building: &mut Vec<i64>,
    out: &mut Vec<i64>,
) {
    let run = building.len();
    if ranges.is_empty() {
        if bp_diff == 0 {
            out.push(0);
        }
        return;
    }
    if run == ranges.len() {
        if bp_diff == 0 {
            out.extend_from_slice(building);
        }
        return;
    }
    let rest = &ranges[run + 1..];
    let rest_lowest: i64 = rest.iter().map(|range| range.lowest * range.period).sum();
    let rest_highest: i64 = rest.iter().map(|range| range.highest * range.period).sum();
    let rest_most_probable: f64 = rest
        .iter()
        .map(|range| range.probability.iter().copied().fold(0.0, f64::max))
        .product();
    let range = &ranges[run];
    for slip in range.lowest..=range.highest {
        let left = bp_diff - slip * range.period;
        if left < rest_lowest || left > rest_highest {
            continue;
        }
        let probability = probability_so_far * range.probability[(slip - range.lowest) as usize];
        if prune && probability * rest_most_probable < NEGLIGIBLE_SLIP_SET_PROBABILITY {
            continue;
        }
        building.push(slip);
        collect_slip_sets(ranges, left, prune, probability, building, out);
        building.pop();
    }
}

/// **The length changes a read of this bundle candidate can show**, as bases added to its length:
/// the totals of every slip set (`str_bundle_loci.md` §7.3, *the union of the slip sets' totals*),
/// and each run's part-repeat changes — the bundle's [`StutterModel::reachable_length_changes_of`].
/// For a bundle of one run these are exactly that run's own.
fn bundle_reachable_length_changes(bundle: BundleCandidate<'_>, out: &mut Vec<i64>) {
    let reach = i64::from(MAX_WHOLE_REPEAT_SLIP);
    let mut totals: std::collections::BTreeSet<i64> = std::iter::once(0).collect();
    for (reference_run, in_allele) in bundle.runs.iter().zip(bundle.parse.runs()) {
        let Some(copies) = NonZeroU32::new(in_allele.copies) else {
            continue;
        };
        let period = reference_run.period() as i64;
        let lowest = (1 - i64::from(copies.get())).max(-reach);
        totals = totals
            .iter()
            .flat_map(|&total| (lowest..=reach).map(move |slip| total + slip * period))
            .collect();
    }
    out.extend(totals);
    for (reference_run, in_allele) in bundle.runs.iter().zip(bundle.parse.runs()) {
        let Some(copies) = NonZeroU32::new(in_allele.copies) else {
            continue;
        };
        let period = period_of(&reference_run.motif());
        out.extend(
            StutterModel::reachable_length_changes_of(period, copies)
                .filter(|change| change % i64::from(period.get()) != 0),
        );
    }
}

// ── STR bundles: reads that end inside the hull (`str_bundle_loci.md` §7.6) ──

/// **Which flank of an STR bundle's hull a read that ended inside the hull held**, read off its
/// witness. A witness from the hull's first position that stops short of its last is a read that
/// held the left flank, and the mirror image one that held the right. A witness over the whole
/// hull is a read that reached the far edge without anchoring there, from either flank; nothing
/// in it says which.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HeldFlank {
    Left,
    Right,
    Either,
}

impl HeldFlank {
    fn of(positions: &WitnessedLocusPositions, hull_len: u32) -> Self {
        let hull_len = LocusLen::from_positions(u64::from(hull_len));
        match (
            positions.is_flush_left(),
            positions.is_flush_right(hull_len),
        ) {
            (true, false) => Self::Left,
            (false, true) => Self::Right,
            _ => Self::Either,
        }
    }
}

/// **The probability that one copy of an STR bundle allele produced a read that held one flank
/// and ended inside the hull** (`str_bundle_loci.md` §7.6) — the bundle's counterpart of a
/// segment's censored emission.
///
/// # What the read says
///
/// A read that held the left flank shows the start of what its carrier has from that flank on:
/// the hull with whatever slips its runs made, then the right flank. Its bases are compared with
/// that sequence's first bases — so it is evidence about every run and every fixed base it got
/// past, priced exactly as a complete read prices them, slips included; a lower bound on the run
/// it ended inside, since every length of that run long enough to hold what it showed agrees with
/// it; and nothing about the runs it never reached, each of which contributes its total
/// probability over its slips whatever the allele holds there. A read that held the right flank
/// is the same read from the other end. One whose witness covers the whole hull is either, each
/// weighted one half.
///
/// ```text
/// P(read | allele) = Σ over every slip set   Π_r P_r(kᵣ) × letters(read | first bases of the allele with the slips, then the flank)
///                  + Σ over every change no single run makes   its part-repeat price × letters(read | first bases of that rendering, then the flank)
/// ```
///
/// That is a complete read's emission summed over every complete read whose first bases are these
/// ones, which is what a segment's censored term is, and for one run over the whole hull with no
/// flank bases it is the segment's own sum (`a_one_run_bundles_partial_read_scores_as_the_segments`).
///
/// # Cost
///
/// A slip set's runs past the read's end do not change the bases compared, so the sum is walked
/// run by run from the held flank, and the read's probability from a run on depends only on how
/// many of the read's bases the pieces before it took. That is remembered, so a read costs a few
/// slips per run per such offset — never the product of every run's slips.
fn bundle_censored_emission(
    witnessed: &[u8],
    witnessed_positions: &WitnessedLocusPositions,
    candidate: &SsrCandidate<'_>,
    bundle: BundleCandidate<'_>,
    context: &SsrScoringContext<'_>,
    scratch: &mut StutterSubstitutionScratch,
) -> f64 {
    let scratch = &mut scratch.bundle;
    match HeldFlank::of(witnessed_positions, bundle.reference_hull_len) {
        held @ (HeldFlank::Left | HeldFlank::Right) => {
            censored_from_one_flank(witnessed, candidate, bundle, context, held, scratch)
        }
        HeldFlank::Either => {
            let left = censored_from_one_flank(
                witnessed,
                candidate,
                bundle,
                context,
                HeldFlank::Left,
                scratch,
            );
            let right = censored_from_one_flank(
                witnessed,
                candidate,
                bundle,
                context,
                HeldFlank::Right,
                scratch,
            );
            0.5 * (left + right)
        }
    }
}

/// [`bundle_censored_emission`] for a read that held `held`, which is `Left` or `Right`.
fn censored_from_one_flank(
    witnessed: &[u8],
    candidate: &SsrCandidate<'_>,
    bundle: BundleCandidate<'_>,
    context: &SsrScoringContext<'_>,
    held: HeldFlank,
    scratch: &mut BundleScoringScratch,
) -> f64 {
    assert_eq!(
        bundle.leading_flank_bases, 0,
        "a read that ended inside the hull is compared with the allele as the table holds it"
    );
    let beyond = flank_beyond_the_hull(candidate, held);
    let censored = &mut scratch.censored;
    let kept = censored.layouts.iter().position(|layout| {
        layout
            .inputs
            .are(held, candidate.bases, beyond, bundle, context)
    });
    let at = match kept {
        Some(at) => at,
        None => {
            let at = if censored.layouts.len() < KEPT_CENSORED_LAYOUTS {
                censored.layouts.push(CensoredLayout::default());
                censored.layouts.len() - 1
            } else {
                let at = censored.next_replaced;
                censored.next_replaced = (at + 1) % KEPT_CENSORED_LAYOUTS;
                at
            };
            censored.layouts[at].lay_out(
                candidate.bases,
                bundle,
                context,
                beyond,
                held,
                LayoutBuffers {
                    segments: &mut censored.segments,
                    placements: &mut censored.placements,
                    resized_stretch: &mut scratch.resized_stretch,
                    rendered: &mut scratch.rendered,
                },
            );
            at
        }
    };
    let layout = &censored.layouts[at];

    // The two per-base scores, looked up once for every comparison below, and the sums of
    // agreeing bases under them.
    censored
        .matched
        .prepare(substitution_scores(context), witnessed.len());
    censored.read.clear();
    if held == HeldFlank::Right {
        censored.read.extend(witnessed.iter().rev());
    } else {
        censored.read.extend_from_slice(witnessed);
    }
    censored.memo.clear();
    censored.memo.resize(
        (2 * layout.slips_of_run.len() + 1) * censored.read.len(),
        f64::NAN,
    );
    let mut total = probability_on_from_fixed_piece(
        layout,
        &censored.read,
        &mut censored.memo,
        0,
        0,
        &censored.matched,
    );

    // The part-repeat renderings: one run's sequencing indel, the others not slipping, for every
    // change no single run's slip makes — the same set and price a complete read takes, each
    // compared from the flank the read held. Consecutive places that show the read the same bases
    // (an indel past where the read ended) are compared once.
    let mut last_compared: Option<(&[u8], f64)> = None;
    for change in &layout.part_repeat_changes {
        let mut priced = 0.0;
        for place in &layout.part_repeat_places[change.priced.clone()] {
            let carrier = &layout.part_repeat_carrier_bases[place.carrier.clone()];
            let letters = match carrier.len().checked_sub(witnessed.len()) {
                // The carrier is shorter than what the read showed.
                None => 0.0,
                Some(from) => {
                    let shown = if held == HeldFlank::Left {
                        &carrier[..witnessed.len()]
                    } else {
                        &carrier[from..]
                    };
                    match last_compared {
                        Some((bases, letters)) if bases == shown => letters,
                        _ => {
                            let letters = censored.matched.letters(witnessed, shown);
                            last_compared = Some((shown, letters));
                            letters
                        }
                    }
                }
            };
            priced += place.probability * letters;
        }
        total += if change.places == 0 {
            0.0
        } else {
            priced / change.places as f64
        };
    }
    total
}

/// **The reference flank a carrier of this allele has past the hull, on the side a read that held
/// `held` ran towards** — less any flank bases the allele's own deletion takes.
fn flank_beyond_the_hull<'a>(candidate: &SsrCandidate<'a>, held: HeldFlank) -> &'a [u8] {
    let (left, right) = candidate.reference_flanks;
    let deletion = candidate.flank_deletion;
    match held {
        HeldFlank::Left => right.get(usize::from(deletion.right)..).unwrap_or_default(),
        _ => left
            .get(..left.len().saturating_sub(usize::from(deletion.left)))
            .unwrap_or_default(),
    }
}

/// The working buffers [`CensoredLayout::lay_out`] builds a layout with.
struct LayoutBuffers<'a> {
    /// A run's stretch split into whole copies and interruptions.
    segments: &'a mut Vec<TractSegment>,
    /// A run's renderings under one slip.
    placements: &'a mut Vec<Vec<u8>>,
    /// One run's stretch resized by a part-repeat change.
    resized_stretch: &'a mut Vec<u8>,
    /// The allele rendered with that change.
    rendered: &'a mut Vec<u8>,
}

impl CensoredLayout {
    /// Lay `allele` out from the flank `held` (`Left` or `Right`): its fixed pieces, the flank
    /// beyond the hull after the last of them, every slip each run can make with its renderings,
    /// and every part-repeat rendering with the flank — all in the order a read from that flank
    /// meets them, reversed from the right except the part-repeat carriers, which are compared
    /// with the read as it is.
    fn lay_out(
        &mut self,
        allele: &[u8],
        bundle: BundleCandidate<'_>,
        context: &SsrScoringContext<'_>,
        beyond: &[u8],
        held: HeldFlank,
        buffers: LayoutBuffers<'_>,
    ) {
        let LayoutBuffers {
            segments,
            placements,
            resized_stretch,
            rendered,
        } = buffers;
        self.inputs.set(held, allele, beyond, bundle, context);
        let ranges = run_slip_ranges(bundle, context);
        let from_right = held == HeldFlank::Right;
        let placed = bundle.parse.runs();
        let run_count = placed.len();

        // The fixed pieces: before the first run, between runs, after the last.
        self.fixed_bases.clear();
        self.fixed.clear();
        for at in 0..=run_count {
            let piece = if from_right { run_count - at } else { at };
            let start = if piece == 0 {
                0
            } else {
                placed[piece - 1].end as usize
            };
            let end = if piece == run_count {
                allele.len()
            } else {
                placed[piece].start as usize
            };
            let first = self.fixed_bases.len();
            if from_right {
                self.fixed_bases.extend(allele[start..end].iter().rev());
                if piece == 0 {
                    self.fixed_bases.extend(beyond.iter().rev());
                }
            } else {
                self.fixed_bases.extend_from_slice(&allele[start..end]);
                if piece == run_count {
                    self.fixed_bases.extend_from_slice(beyond);
                }
            }
            self.fixed.push(first..self.fixed_bases.len());
        }

        // Every slip each run can make, with its renderings.
        self.rendering_bases.clear();
        self.renderings.clear();
        self.slips.clear();
        self.slips_of_run.clear();
        for at in 0..run_count {
            let run = if from_right { run_count - 1 - at } else { at };
            let stretch = &allele[placed[run].start as usize..placed[run].end as usize];
            let motif = bundle.runs[run].motif();
            let range = &ranges[run];
            if placed[run].copies > 0 {
                segment_tract_into(stretch, motif.as_bytes(), segments);
            }
            let first_slip = self.slips.len();
            for slip in range.lowest..=range.highest {
                let probability = range.probability[(slip - range.lowest) as usize];
                if probability <= 0.0 {
                    continue;
                }
                let first_rendering = self.renderings.len();
                if slip == 0 {
                    push_rendering(
                        &mut self.rendering_bases,
                        &mut self.renderings,
                        stretch,
                        from_right,
                    );
                } else {
                    place_slip_in_each_run(segments, motif.as_bytes(), slip, placements);
                    for placement in placements.iter() {
                        push_rendering(
                            &mut self.rendering_bases,
                            &mut self.renderings,
                            placement,
                            from_right,
                        );
                    }
                }
                // A slip no place in the run can absorb is not made, as at a complete read.
                if self.renderings.len() == first_rendering {
                    continue;
                }
                let length = self.renderings[first_rendering].len();
                self.slips.push(CensoredSlip {
                    probability,
                    renderings: first_rendering..self.renderings.len(),
                    length,
                });
            }
            self.slips_of_run.push(first_slip..self.slips.len());
        }

        self.unreached_mass.clear();
        self.unreached_mass.resize(run_count + 1, 1.0);
        for at in (0..run_count).rev() {
            let run_mass: f64 = self.slips[self.slips_of_run[at].clone()]
                .iter()
                .map(|slip| slip.probability)
                .sum();
            self.unreached_mass[at] = run_mass * self.unreached_mass[at + 1];
        }

        // The part-repeat renderings, each with the flank on the side the read ran towards. A
        // part-repeat change is at most `2 · MAX_PART_REPEAT_SLIP` bases at any period of two or
        // more.
        self.part_repeat_changes.clear();
        self.part_repeat_places.clear();
        self.part_repeat_carrier_bases.clear();
        let widest = 2 * i64::from(MAX_PART_REPEAT_SLIP);
        for bp_diff in -widest..=widest {
            if a_single_run_slip_reaches(&ranges, bp_diff) {
                continue;
            }
            let first_place = self.part_repeat_places.len();
            let places = part_repeat_places(
                allele,
                bundle,
                context,
                bp_diff,
                (&mut *resized_stretch, &mut *rendered),
                |probability, rendered| {
                    let first = self.part_repeat_carrier_bases.len();
                    if from_right {
                        self.part_repeat_carrier_bases.extend_from_slice(beyond);
                        self.part_repeat_carrier_bases.extend_from_slice(rendered);
                    } else {
                        self.part_repeat_carrier_bases.extend_from_slice(rendered);
                        self.part_repeat_carrier_bases.extend_from_slice(beyond);
                    }
                    self.part_repeat_places.push(PartRepeatPlace {
                        probability,
                        carrier: first..self.part_repeat_carrier_bases.len(),
                    });
                },
            );
            self.part_repeat_changes.push(PartRepeatChange {
                places,
                priced: first_place..self.part_repeat_places.len(),
            });
        }
    }
}

/// **The read's probability from fixed piece `at` of `layout` on**, its first `offset` bases
/// already accounted for by the pieces before. `read` is in the layout's orientation.
fn probability_on_from_fixed_piece(
    layout: &CensoredLayout,
    read: &[u8],
    memo: &mut [f64],
    at: usize,
    offset: usize,
    matched: &MatchedPrefixLetters,
) -> f64 {
    let read_len = read.len();
    if offset >= read_len {
        return layout.unreached_mass[at];
    }
    // A fixed piece is reached at one offset by several paths as often as a run is.
    let key = (layout.slips_of_run.len() + at) * read_len + offset;
    if !memo[key].is_nan() {
        return memo[key];
    }
    let piece = layout.fixed[at].clone();
    let compared = piece.len().min(read_len - offset);
    let letters_here = matched.letters(
        &read[offset..offset + compared],
        &layout.fixed_bases[piece.start..piece.start + compared],
    );
    let next = offset + piece.len();
    let total = if next >= read_len {
        letters_here * layout.unreached_mass[at]
    } else if at == layout.slips_of_run.len() {
        // The carrier's sequence ran out before the read did.
        0.0
    } else {
        letters_here * probability_on_from_run(layout, read, memo, at, next, matched)
    };
    memo[key] = total;
    total
}

/// **The read's probability from run `at` of `layout` on**, its first `offset` bases already
/// accounted for: each slip of the run, its renderings' letters averaged, times what follows it.
///
/// Renderings that show the read the same bases — every slip of a pure run long enough to hold
/// the rest of the read — are compared once: the letters are a function of those bases alone.
fn probability_on_from_run(
    layout: &CensoredLayout,
    read: &[u8],
    memo: &mut [f64],
    at: usize,
    offset: usize,
    matched: &MatchedPrefixLetters,
) -> f64 {
    let read_len = read.len();
    let key = at * read_len + offset;
    if !memo[key].is_nan() {
        return memo[key];
    }
    let mut last_compared: Option<(&[u8], f64)> = None;
    let mut total = 0.0;
    for slip in &layout.slips[layout.slips_of_run[at].clone()] {
        let compared = slip.length.min(read_len - offset);
        let each = 1.0 / slip.renderings.len() as f64;
        let letters: f64 = layout.renderings[slip.renderings.clone()]
            .iter()
            .map(|bases| {
                let shown = &layout.rendering_bases[bases.start..bases.start + compared];
                let letters = match last_compared {
                    Some((before, letters)) if before == shown => letters,
                    _ => {
                        let letters = matched.letters(&read[offset..offset + compared], shown);
                        last_compared = Some((shown, letters));
                        letters
                    }
                };
                each * letters
            })
            .sum();
        let next = offset + slip.length;
        let rest = if next >= read_len {
            layout.unreached_mass[at + 1]
        } else {
            probability_on_from_fixed_piece(layout, read, memo, at + 1, next, matched)
        };
        total += slip.probability * letters * rest;
    }
    memo[key] = total;
    total
}

/// Append one rendering's bases, reversed for a read from the right flank, and its range.
fn push_rendering(
    rendering_bases: &mut Vec<u8>,
    renderings: &mut Vec<std::ops::Range<usize>>,
    bases: &[u8],
    reversed: bool,
) {
    let first = rendering_bases.len();
    if reversed {
        rendering_bases.extend(bases.iter().rev());
    } else {
        rendering_bases.extend_from_slice(bases);
    }
    renderings.push(first..rendering_bases.len());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::MarginalAligner;
    use crate::alignment::ssr_marginal_sequence::{SequenceMarginalScratch, SsrSequenceMarginal};
    use crate::alignment::stutter::MAX_WHOLE_REPEAT_SLIP;
    use crate::calling::likelihood::stutter_rates::stutter_model_for;
    use crate::parameter_estimation::joint::ssr_fit::Slippage;

    /// **Model B — the comparator, and an oracle rather than an alternative.**
    ///
    /// `Σ_n S(n) · avg_v align(observation | candidate ⊕ n repeats)`: the whole-repeat
    /// stutter mass marginalised **outside** a sequence-versus-sequence align that absorbs
    /// whatever length is left over as slop at the tract's ends.
    ///
    /// # How it differs from Model A, which is the whole point
    ///
    /// **Model A picks one length change** — the one the observation actually shows — scores
    /// its probability, and compares letters over sequences it has already made equal. **Model
    /// B sums over every whole-repeat change**, and lets the aligner explain the residual.
    /// Where A has an explicit part-repeat branch with its own two parameters, B has none: a
    /// change that is not a whole number of repeats reaches B only as end-gap slop on a
    /// neighbouring term.
    ///
    /// So the two explain a read's length by genuinely different routes, and that is what
    /// makes agreement between them evidence. **What they share is the placement enumeration**
    /// — production shares it too — so the independence is in how length is explained, not in
    /// where a slip may land.
    ///
    /// **Test-only, deliberately.** It is worth more as a second opinion than as a model
    /// anyone runs; production keeps its own Model B the same way (spec §9).
    #[derive(Debug, Default, Clone, Copy)]
    struct ClassicEmissionOracle;

    #[derive(Debug, Default)]
    struct ClassicOracleScratch {
        placements: Vec<Vec<u8>>,
        aligner: SequenceMarginalScratch,
    }

    impl SsrEmissionModel for ClassicEmissionOracle {
        type Scratch = ClassicOracleScratch;

        fn emission(
            &self,
            observation: &[u8],
            candidate: &SsrCandidate<'_>,
            context: &SsrScoringContext<'_>,
            scratch: &mut Self::Scratch,
        ) -> f64 {
            let period = period_of(context.motif);
            let period_bases = i64::from(period.get());
            let motif = context.motif.as_bytes();
            let aligner = SsrSequenceMarginal::try_new(context.substitution_rate.get())
                .expect("an ErrorRate is a probability by construction");

            let cutoff = i64::from(MAX_WHOLE_REPEAT_SLIP);
            let mut total = 0.0;
            for repeats in -cutoff..=cutoff {
                let length_probability =
                    context.stutter.probability(repeats * period_bases, period);
                if length_probability <= 0.0 {
                    continue;
                }
                enumerate_placements(candidate.bases, motif, repeats, &mut scratch.placements);
                if scratch.placements.is_empty() {
                    continue;
                }
                let each = 1.0 / scratch.placements.len() as f64;
                let letters: f64 = scratch
                    .placements
                    .iter()
                    .map(|placement| {
                        each * float::exp(
                            aligner
                                .marginal_probability(
                                    observation,
                                    placement,
                                    (),
                                    &mut scratch.aligner,
                                )
                                .get(),
                        )
                    })
                    .sum();
                total += length_probability * letters;
            }
            total
        }

        fn censored_emission(
            &self,
            _witnessed_prefix: &[u8],
            _witnessed_positions: &WitnessedLocusPositions,
            _candidate: &SsrCandidate<'_>,
            _context: &SsrScoringContext<'_>,
            _scratch: &mut Self::Scratch,
        ) -> f64 {
            unimplemented!("the oracle scores complete observations only")
        }
    }

    fn a_motif(bases: &[u8]) -> Motif {
        Motif::new(bases).expect("a valid test motif")
    }

    fn repeats(count: u32) -> NonZeroU32 {
        NonZeroU32::new(count).expect("a test candidate always holds a repeat")
    }

    fn a_model() -> StutterModel {
        stutter_model_for(&Slippage {
            level: 0.02,
            shorter_share: 0.83,
            fall_off: 0.35,
        })
    }

    /// **The context takes its unreachable mass from the distribution, not from the caller.**
    /// Two ways in would be two chances to disagree, and the row compares candidates on the
    /// strength of this number.
    #[test]
    fn the_context_reads_the_unreachable_mass_off_the_distribution() {
        let motif = a_motif(b"CAG");
        let model = a_model();
        let candidate = SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: b"AAACAGCAGCAGCAGTTT",
            repeat_count: repeats(4),
        };
        let context = SsrScoringContext::new(
            &motif,
            &model,
            &candidate,
            ErrorRate::try_new(0.001).expect("a valid rate"),
            [Provenance::FittedHere],
        );
        assert_eq!(
            context.unreachable_mass,
            model.unreachable_mass(period_of(&motif), repeats(4))
        );
    }

    /// **The mass differs between candidates at one locus**, which is the whole reason a
    /// context is per candidate rather than per locus. A four-repeat tract cannot lose as
    /// many repeats as a thirty-repeat one, so it leaves more unplaced.
    #[test]
    fn two_candidates_at_one_locus_get_different_contexts() {
        let motif = a_motif(b"CAG");
        let model = a_model();
        let rate = ErrorRate::try_new(0.001).expect("a valid rate");

        let short = SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: b"AAACAGCAGCAGCAGTTT",
            repeat_count: repeats(4),
        };
        let long = SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: b"AAACAGCAGCAGCAGTTT",
            repeat_count: repeats(30),
        };
        let short_context =
            SsrScoringContext::new(&motif, &model, &short, rate, [Provenance::FittedHere]);
        let long_context =
            SsrScoringContext::new(&motif, &model, &long, rate, [Provenance::FittedHere]);

        assert!(
            short_context.unreachable_mass > long_context.unreachable_mass,
            "four repeats left {} unplaced, thirty left {}",
            short_context.unreachable_mass,
            long_context.unreachable_mass
        );
    }

    /// **The weakest warrant wins, and one borrowed parameter is enough.** A context resting
    /// on a fitted rate and a borrowed slippage row is a borrowed context; stamping it
    /// `FittedHere` would launder the weaker of the two.
    #[test]
    fn the_context_carries_the_weakest_warrant_that_entered_it() {
        let motif = a_motif(b"CA");
        let model = a_model();
        let candidate = SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: b"AACACACATT",
            repeat_count: repeats(3),
        };
        let rate = ErrorRate::try_new(0.001).expect("a valid rate");

        for (warrants, expected) in [
            (
                vec![Provenance::FittedHere, Provenance::FittedHere],
                Provenance::FittedHere,
            ),
            (
                vec![Provenance::FittedHere, Provenance::Borrowed],
                Provenance::Borrowed,
            ),
            (
                vec![Provenance::Borrowed, Provenance::Supplied],
                Provenance::Supplied,
            ),
            (
                vec![Provenance::Supplied, Provenance::Defaulted],
                Provenance::Defaulted,
            ),
            (
                vec![Provenance::Defaulted, Provenance::FittedHere],
                Provenance::Defaulted,
            ),
        ] {
            let context =
                SsrScoringContext::new(&motif, &model, &candidate, rate, warrants.clone());
            assert_eq!(
                context.weakest_provenance, expected,
                "warrants {warrants:?} gave {:?}",
                context.weakest_provenance
            );
        }
    }

    /// **No warrants at all is `FittedHere`**, which is the identity of the fold rather than
    /// an opinion: a context nothing weakened is as well founded as its inputs. Pinned
    /// because the alternative — defaulting to `Defaulted` — would mark every context of a
    /// fully-fitted run as a guess.
    #[test]
    fn a_context_with_no_weakening_warrant_is_fitted() {
        let motif = a_motif(b"CA");
        let model = a_model();
        let candidate = SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: b"AACACACATT",
            repeat_count: repeats(3),
        };
        let context = SsrScoringContext::new(
            &motif,
            &model,
            &candidate,
            ErrorRate::try_new(0.001).expect("a valid rate"),
            [],
        );
        assert_eq!(context.weakest_provenance, Provenance::FittedHere);
    }
    /// Score one observation against one candidate under Model A, at a stated substitution
    /// rate — the shape every test below uses.
    fn score(
        observation: &[u8],
        candidate_bases: &[u8],
        repeat_count: u32,
        motif: &[u8],
        model: &StutterModel,
        substitution_rate: f64,
    ) -> f64 {
        let motif = a_motif(motif);
        let candidate = SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: candidate_bases,
            repeat_count: repeats(repeat_count),
        };
        let context = SsrScoringContext::new(
            &motif,
            model,
            &candidate,
            ErrorRate::try_new(substitution_rate).expect("a valid rate"),
            [Provenance::FittedHere],
        );
        let mut scratch = StutterSubstitutionScratch::default();
        StutterSubstitutionEmission.emission(observation, &candidate, &context, &mut scratch)
    }

    /// Rank the candidates by one model's emission, best first.
    fn ranking<M: SsrEmissionModel>(
        model: &M,
        observation: &[u8],
        candidates: &[(Vec<u8>, u32)],
        motif: &Motif,
        stutter: &StutterModel,
        epsilon: f64,
        scratch: &mut M::Scratch,
    ) -> (Vec<usize>, Vec<f64>) {
        let scores: Vec<f64> = candidates
            .iter()
            .map(|(bases, count)| {
                let candidate = SsrCandidate {
                    bundle: None,
                    reference_flanks: (&[], &[]),
                    flank_deletion: crate::locus_generation::FlankDeletion::NONE,
                    bases,
                    repeat_count: repeats(*count),
                };
                let context = SsrScoringContext::new(
                    motif,
                    stutter,
                    &candidate,
                    ErrorRate::try_new(epsilon).expect("a valid rate"),
                    [Provenance::FittedHere],
                );
                model.emission(observation, &candidate, &context, scratch)
            })
            .collect();
        let mut order: Vec<usize> = (0..scores.len()).collect();
        order.sort_by(|a, b| scores[*b].partial_cmp(&scores[*a]).expect("no NaN scores"));
        (order, scores)
    }

    /// **The independent-implementation check: the two models rank candidates alike wherever
    /// the observation's length differs from the candidate's by a whole number of repeats.**
    ///
    /// That is the case Model B is built for, and the two reach it by genuinely different
    /// routes — A picks the single length change the read shows and compares letters over
    /// equal-length sequences; B sums over *every* whole-repeat change and lets the aligner
    /// absorb the residual. Agreement across that grid is evidence about both.
    ///
    /// Measured here: the full ranking is identical for every observation, and the winning
    /// score agrees to better than one part in ten thousand.
    #[test]
    fn the_two_models_rank_candidates_alike_on_whole_repeat_differences() {
        let stutter = a_model();
        let motif = a_motif(b"CAG");
        let epsilon = 0.001;
        let candidates: Vec<(Vec<u8>, u32)> = (3..=6)
            .map(|count| (b"CAG".repeat(count), count as u32))
            .collect();

        // Whole numbers of repeats, plus one carrying a substitution — every one a whole-repeat
        // difference from every candidate.
        let observations: Vec<Vec<u8>> = vec![
            b"CAGCAGCAG".to_vec(),
            b"CAGCAGCAGCAG".to_vec(),
            b"CAGCAGCAGCAGCAG".to_vec(),
            b"CAGCTGCAGCAG".to_vec(),
        ];

        for observation in &observations {
            let mut a_scratch = StutterSubstitutionScratch::default();
            let mut b_scratch = ClassicOracleScratch::default();
            let (a_order, a_scores) = ranking(
                &StutterSubstitutionEmission,
                observation,
                &candidates,
                &motif,
                &stutter,
                epsilon,
                &mut a_scratch,
            );
            let (b_order, b_scores) = ranking(
                &ClassicEmissionOracle,
                observation,
                &candidates,
                &motif,
                &stutter,
                epsilon,
                &mut b_scratch,
            );

            assert_eq!(
                a_order,
                b_order,
                "the two models ranked {} differently: {a_scores:?} against {b_scores:?}",
                String::from_utf8_lossy(observation)
            );

            let winner = a_order[0];
            let relative = (a_scores[winner] - b_scores[winner]).abs() / a_scores[winner];
            assert!(
                relative < 1e-4,
                "the winning scores differ by {relative} on {}",
                String::from_utf8_lossy(observation)
            );
        }
    }

    /// **Where the two part company, and why that is the choice rather than a defect.**
    ///
    /// A read whose length is *not* a whole number of repeats from the candidate is the one
    /// case the two models explain differently by construction. **Model A charges the fitted
    /// part-repeat share**; **Model B has no such branch at all** and absorbs the odd base as
    /// slop at the tract's end, at the flat per-base rate.
    ///
    /// At ε = 0.001 the slop route is about nine times cheaper than the part-repeat event, so
    /// **B keeps the shorter candidate and A prefers the one needing the smaller stutter
    /// event**. Measured on a ten-base read against a three-repeat and a four-repeat
    /// candidate: A scores 1.094e-4 and 1.869e-4 — the four-repeat candidate wins — while B
    /// scores 9.706e-4 and 1.167e-5, and the three-repeat candidate wins.
    ///
    /// **Model A is the one to trust here**, and this divergence is the reason: an explicit
    /// part-repeat branch with its own two fitted parameters is what the comparison behind
    /// spec §4.1 chose it for. The test exists so the difference stays a measured, deliberate
    /// thing rather than a surprise the first time a ranking disagrees.
    #[test]
    fn the_two_models_part_company_on_a_part_repeat_length() {
        let stutter = a_model();
        let motif = a_motif(b"CAG");
        let epsilon = 0.001;
        let candidates: Vec<(Vec<u8>, u32)> =
            vec![(b"CAGCAGCAG".to_vec(), 3), (b"CAGCAGCAGCAG".to_vec(), 4)];
        // Three repeats and one base over — a part-repeat difference from both candidates.
        let observation = b"CAGCAGCAGC";

        let mut a_scratch = StutterSubstitutionScratch::default();
        let mut b_scratch = ClassicOracleScratch::default();
        let (a_order, a_scores) = ranking(
            &StutterSubstitutionEmission,
            observation,
            &candidates,
            &motif,
            &stutter,
            epsilon,
            &mut a_scratch,
        );
        let (b_order, b_scores) = ranking(
            &ClassicEmissionOracle,
            observation,
            &candidates,
            &motif,
            &stutter,
            epsilon,
            &mut b_scratch,
        );

        assert_eq!(
            a_order[0], 1,
            "Model A should prefer the four-repeat candidate"
        );
        assert_eq!(
            b_order[0], 0,
            "Model B should prefer the three-repeat candidate"
        );

        // The sizes, so the mechanism is checkable and not merely asserted.
        assert!((a_scores[0] - 1.094e-4).abs() < 1e-7, "{:?}", a_scores);
        assert!((a_scores[1] - 1.869e-4).abs() < 1e-7, "{:?}", a_scores);
        assert!((b_scores[0] - 9.706e-4).abs() < 1e-6, "{:?}", b_scores);
        assert!((b_scores[1] - 1.167e-5).abs() < 1e-7, "{:?}", b_scores);

        // B's route is the flat per-base rate; A's is the fitted part-repeat product.
        let slop_route = b_scores[0] / stutter.same_length_share();
        let part_repeat_route =
            stutter.part_repeat_longer_share() * stutter.part_repeat_one_step_share();
        assert!(
            slop_route > part_repeat_route * 8.0,
            "slop {slop_route} against a part-repeat event {part_repeat_route}"
        );
    }

    /// **Spec §12's first test: a read matching its allele is dominated by the same-length
    /// term.** A read identical to the candidate scores at least
    /// `same_length_share × (1 − ε)^length` — to within rounding, since the score reaches that
    /// product by a sum of logarithms — and within 5% of it: the two factors and nothing else, with
    /// no placement multiplicity to dilute it.
    #[test]
    fn a_read_identical_to_its_candidate_is_dominated_by_the_same_length_term() {
        let model = a_model();
        let tract = b"CAGCAGCAGCAGCAGCAG";
        let epsilon = 0.001;
        let scored = score(tract, tract, 6, b"CAG", &model, epsilon);

        let floor = model.same_length_share() * float::powi(1.0 - epsilon, tract.len() as i32);
        // The score and the floor reach the same product by different routes (a sum of logarithms
        // against repeated multiplication), so they may round a few units in the last place apart.
        assert!(
            scored >= floor * (1.0 - 1e-12),
            "{scored} fell below {floor}"
        );
        assert!(
            scored <= floor * 1.05,
            "{scored} exceeded {floor} by more than a twentieth"
        );
    }

    /// **Spec §12's second test: direction and size are ordered as fitted.** With a
    /// contraction-biased split, a read one repeat short outscores one a repeat long; and one
    /// repeat outscores two. Both are properties of the fitted numbers, so a model that
    /// dropped the direction split or the size decay would fail here rather than silently
    /// score a symmetric distribution.
    #[test]
    fn a_shorter_read_outscores_a_longer_one_and_one_repeat_outscores_two() {
        let model = a_model();
        let candidate = b"CAGCAGCAGCAGCAGCAG";
        let epsilon = 0.001;

        let one_short = score(b"CAGCAGCAGCAGCAG", candidate, 6, b"CAG", &model, epsilon);
        let one_long = score(
            b"CAGCAGCAGCAGCAGCAGCAG",
            candidate,
            6,
            b"CAG",
            &model,
            epsilon,
        );
        assert!(
            one_short > one_long,
            "a repeat short scored {one_short}, a repeat long {one_long}"
        );

        let two_short = score(b"CAGCAGCAGCAG", candidate, 6, b"CAG", &model, epsilon);
        assert!(
            one_short > two_short,
            "one repeat short scored {one_short}, two short {two_short}"
        );
    }

    /// **Spec §12's third test: a whole repeat beats a single stray base** — and under the
    /// **corrected** condition, which is a comparison of *products* rather than of the
    /// direction shares alone.
    ///
    /// "Any part-repeat share below the whole-repeat one" suffices only while the two one-step
    /// shares are tied, and spec §10 schedules untying them, so the test states the condition
    /// it actually depends on and asserts it holds for this fixture before asserting the
    /// consequence. The specification's own counter-example: a whole-repeat share of 0.02 at a
    /// one-step share of 0.1 gives 0.002, against a part-repeat share of 0.019 at 0.95 giving
    /// 0.018 — nine times higher, both inside the clamps.
    #[test]
    fn a_whole_repeat_beats_a_stray_base_when_the_products_say_so() {
        let model = a_model();
        let candidate = b"CAGCAGCAGCAGCAGCAG";
        let epsilon = 0.001;

        let whole_product = model.whole_repeat_longer_share() * model.whole_repeat_one_step_share();
        let part_product = model.part_repeat_longer_share() * model.part_repeat_one_step_share();
        assert!(
            part_product < whole_product,
            "this fixture does not meet the condition: {part_product} against {whole_product}"
        );

        let whole_repeat_longer = score(
            b"CAGCAGCAGCAGCAGCAGCAG",
            candidate,
            6,
            b"CAG",
            &model,
            epsilon,
        );
        let one_base_longer = score(
            b"CAGCAGCAGCAGCAGCAGC",
            candidate,
            6,
            b"CAG",
            &model,
            epsilon,
        );
        assert!(
            whole_repeat_longer > one_base_longer,
            "a whole repeat scored {whole_repeat_longer}, a stray base {one_base_longer}"
        );
    }

    /// **An interrupted tract has more than one place a slip could land, and they are averaged
    /// with equal weight.** Two runs either side of an interruption give two distinct
    /// sequences for a one-repeat expansion; a read matching one of them is scored at half the
    /// letters' probability, because the model does not know which run slipped.
    #[test]
    fn an_interrupted_tract_averages_over_the_placements_a_slip_could_land_in() {
        let mut placements = Vec::new();
        enumerate_placements(b"CACACATTCACA", b"CA", 1, &mut placements);
        assert_eq!(placements.len(), 2, "{placements:?}");
        assert!(placements.contains(&b"CACACACATTCACA".to_vec()));
        assert!(placements.contains(&b"CACACATTCACACA".to_vec()));

        // A pure tract has one placement, so nothing is diluted there.
        let mut pure = Vec::new();
        enumerate_placements(b"CACACACACA", b"CA", 1, &mut pure);
        assert_eq!(pure.len(), 1);
        assert_eq!(pure[0], b"CACACACACACA".to_vec());
    }

    /// **A contraction no single run can absorb is unreachable, and scores exactly zero** —
    /// not a small number. The mass that costs is what `SsrScoringContext::unreachable_mass`
    /// reports, which is why it is reported rather than assumed away.
    ///
    /// **It takes an interrupted tract to reach this branch through `emission` at all**, and
    /// that is worth stating: a slip lands in *one* run, so a tract of two runs of two
    /// repeats cannot lose three even though it holds four. On a **pure** tract the branch is
    /// unreachable from `emission`, because the deepest contraction an observation can ask
    /// for is the one that leaves it empty — and a single run can always give up everything
    /// it holds.
    #[test]
    fn a_contraction_no_single_run_can_absorb_scores_zero() {
        let model = a_model();
        // Two runs of two repeats, four in total: no run can give up three.
        let interrupted = b"CACATTCACA";
        let scored = score(b"CACA", interrupted, 4, b"CA", &model, 0.001);
        assert_eq!(scored, 0.0);

        let mut placements = Vec::new();
        enumerate_placements(interrupted, b"CA", -3, &mut placements);
        assert!(placements.is_empty(), "{placements:?}");

        // Two repeats it can give up, from either run — and emptying one run is not the same
        // sequence as emptying the other, so both placements stand: `TTCACA` and `CACATT`.
        enumerate_placements(interrupted, b"CA", -2, &mut placements);
        assert_eq!(placements.len(), 2, "{placements:?}");
        assert!(placements.contains(&b"TTCACA".to_vec()));
        assert!(placements.contains(&b"CACATT".to_vec()));

        // **A pure tract may not vanish entirely**, which is the rule `unreachable_mass`
        // states from the distribution's side. Four repeats can give up three, not four.
        enumerate_placements(b"CACACACA", b"CA", -3, &mut placements);
        assert_eq!(placements.len(), 1, "{placements:?}");
        enumerate_placements(b"CACACACA", b"CA", -4, &mut placements);
        assert!(placements.is_empty(), "{placements:?}");
    }

    /// **A part-repeat change is resized at the tract's end, in one placement** — production's
    /// split, and a simplification the specification records rather than hides. Lengthening
    /// continues the motif's tiling from the phase the tract ended on.
    #[test]
    fn a_part_repeat_change_is_resized_at_the_tract_end() {
        let mut resized = Vec::new();

        resize_at_the_end(b"CAGCAGCAG", b"CAG", 1, &mut resized);
        assert_eq!(resized, b"CAGCAGCAGC".to_vec());

        resize_at_the_end(b"CAGCAGCAG", b"CAG", 2, &mut resized);
        assert_eq!(resized, b"CAGCAGCAGCA".to_vec());

        resize_at_the_end(b"CAGCAGCAG", b"CAG", -1, &mut resized);
        assert_eq!(resized, b"CAGCAGCA".to_vec());
    }

    /// **The two factors are separable, and the substitution one is the alignment module's.**
    /// At a candidate the read matches exactly, the score is the same-length share times the
    /// flat model's own product over the tract — so a change to either factor moves the score
    /// by exactly that factor and nothing else.
    #[test]
    fn the_score_is_the_length_factor_times_the_letter_factor() {
        let model = a_model();
        let tract = b"CAGCAGCAGCAGCAGCAG";

        for epsilon in [1e-4, 1e-3, 1e-2] {
            let scored = score(tract, tract, 6, b"CAG", &model, epsilon);
            let letters = float::powi(1.0 - epsilon, tract.len() as i32);
            let expected = model.same_length_share() * letters;
            assert!(
                (scored - expected).abs() <= expected * 1e-12,
                "at ε={epsilon}: {scored} against {expected}"
            );
        }

        // One mismatching base costs exactly one `ε/3` in place of one `1 − ε`.
        let epsilon = 1e-3;
        let mismatched = b"TAGCAGCAGCAGCAGCAG";
        let scored = score(mismatched, tract, 6, b"CAG", &model, epsilon);
        let expected = model.same_length_share()
            * float::powi(1.0 - epsilon, tract.len() as i32 - 1)
            * (epsilon / 3.0);
        assert!(
            (scored - expected).abs() <= expected * 1e-12,
            "{scored} against {expected}"
        );
    }

    /// Score a read that ran out inside the tract, at a stated substitution rate — the
    /// counterpart of [`score`], and the shape every censored test below uses.
    fn score_censored(
        witnessed_prefix: &[u8],
        candidate_bases: &[u8],
        repeat_count: u32,
        motif: &[u8],
        model: &StutterModel,
        substitution_rate: f64,
    ) -> f64 {
        let motif = a_motif(motif);
        let candidate = SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: candidate_bases,
            repeat_count: repeats(repeat_count),
        };
        let context = SsrScoringContext::new(
            &motif,
            model,
            &candidate,
            ErrorRate::try_new(substitution_rate).expect("a valid rate"),
            [Provenance::FittedHere],
        );
        let mut scratch = StutterSubstitutionScratch::default();
        StutterSubstitutionEmission.censored_emission(
            witnessed_prefix,
            &WitnessedLocusPositions::from_half_open_runs([(0, 1)]).expect("a witness"),
            &candidate,
            &context,
            &mut scratch,
        )
    }

    /// The cell of a sweep where a censored read and a complete one disagreed most about how
    /// far apart two candidates are — carried as three named numbers rather than a tuple,
    /// because the cell is read sixty lines from where it is built.
    #[derive(Clone, Copy, Default)]
    struct WidestDisagreement {
        /// How far apart the two separations were, in nats.
        gap: f64,
        /// What the complete read made of the two candidates, in nats.
        complete_separation: f64,
        /// What the censored read made of them.
        censored_separation: f64,
    }

    /// A tract of `repeat_count` whole copies of `motif`.
    fn a_tract(motif: &[u8], repeat_count: usize) -> Vec<u8> {
        motif.repeat(repeat_count)
    }

    /// **Where the read's own length admits exactly one stretching, a censored read scores
    /// bit for bit what the complete read of the same bases scores** — spec §12's twelfth
    /// test, second half.
    ///
    /// **This is a test of the tail's arithmetic and not of a tolerance.** The widest
    /// stretching either branch reaches is [`MAX_WHOLE_REPEAT_SLIP`] whole repeats, so a read
    /// that got through exactly that much tract leaves the tail one term — and a tail summed
    /// from the distribution's own terms is then that term itself, to the bit. A tail written
    /// as a telescoped geometric difference would agree only to within rounding, which is why
    /// it is not written that way.
    ///
    /// Both branches are exercised: the pure tract takes the factorised form, the interrupted
    /// one the exact sum over stretchings.
    #[test]
    fn a_censored_read_of_one_admissible_length_scores_what_the_complete_read_does() {
        let model = a_model();
        let epsilon = 1e-3;

        // Pure: the factorised branch.
        for (motif, repeat_count) in [
            (b"A".as_slice(), 6u32),
            (b"CA".as_slice(), 5),
            (b"CAG".as_slice(), 4),
        ] {
            let candidate = a_tract(motif, repeat_count as usize);
            let widest = motif.len() * MAX_WHOLE_REPEAT_SLIP as usize;
            let observation = a_tract(
                motif,
                repeat_count as usize + MAX_WHOLE_REPEAT_SLIP as usize,
            );
            assert_eq!(observation.len(), candidate.len() + widest);

            let complete = score(
                &observation,
                &candidate,
                repeat_count,
                motif,
                &model,
                epsilon,
            );
            let censored = score_censored(
                &observation,
                &candidate,
                repeat_count,
                motif,
                &model,
                epsilon,
            );
            assert!(complete > 0.0, "the fixture stopped scoring anything");
            assert_eq!(
                censored.to_bits(),
                complete.to_bits(),
                "motif {}: censored {censored}, complete {complete}",
                String::from_utf8_lossy(motif)
            );
        }

        // Interrupted: the exact-sum branch. Two runs of two `CAG` around a `TT`, so a
        // whole-repeat slip has two placements and the letters cannot factor out.
        let candidate = b"CAGCAGTTCAGCAG";
        let mut observation = a_tract(b"CAG", 12);
        observation.extend_from_slice(b"TT");
        observation.extend_from_slice(&a_tract(b"CAG", 2));
        assert_eq!(observation.len(), candidate.len() + 30);

        let complete = score(&observation, candidate, 4, b"CAG", &model, epsilon);
        let censored = score_censored(&observation, candidate, 4, b"CAG", &model, epsilon);
        assert!(complete > 0.0, "the interrupted fixture scores nothing");
        assert_eq!(
            censored.to_bits(),
            complete.to_bits(),
            "interrupted: censored {censored}, complete {complete}"
        );
    }

    /// **A lower bound is never less likely than the exact length it bounds.** The censored
    /// score sums the complete read's own term together with every longer stretching, so it
    /// can only be larger — and where the read's length admits one stretching alone the two
    /// meet, which the test above pins to the bit.
    ///
    /// Pinned because it is the cheapest check that the tail sums **upward** from its floor: a
    /// tail that summed the wrong direction, or that made its floor strict, breaks it at once.
    ///
    /// **It does not catch a floor dropped altogether**, and an earlier version of this comment
    /// claimed it did. Dropping the floor makes the tail *larger*, which an inequality in this
    /// direction accepts by construction. Measured: replacing the floor with `i64::MIN` leaves
    /// this test green and fails four others, of which
    /// `a_censored_read_past_every_stretching_scores_nothing` is the one to read.
    #[test]
    fn a_censored_read_is_at_least_as_likely_as_the_complete_read_of_the_same_bases() {
        let model = a_model();
        let epsilon = 1e-3;
        let motif = b"CA";

        for repeat_count in [2u32, 4, 8] {
            let candidate = a_tract(motif, repeat_count as usize);
            for witnessed_repeats in 1..=(repeat_count as usize + 4) {
                let observation = a_tract(motif, witnessed_repeats);
                let complete = score(
                    &observation,
                    &candidate,
                    repeat_count,
                    motif,
                    &model,
                    epsilon,
                );
                let censored = score_censored(
                    &observation,
                    &candidate,
                    repeat_count,
                    motif,
                    &model,
                    epsilon,
                );
                assert!(
                    censored >= complete,
                    "{repeat_count} repeats, read of {witnessed_repeats}: censored {censored} \
                     below complete {complete}"
                );
            }
        }
    }

    /// **Spec §12's thirteenth test, measured — and what the measurement changed about it.**
    ///
    /// For two candidates the read **outgrew** — both shorter than the stretch it got through
    /// — a partial observation and the complete observation of the same bases separate them by
    /// **very nearly the same amount**, not by an amount ordered one way.
    ///
    /// The reason is the geometric's memorylessness: above the read's own length the tail is
    /// proportional to the point probability, so the two log-likelihood ratios cancel down to
    /// the same number. What breaks the exact equality is the **part-repeat branch**, whose
    /// re-indexing means its terms do not line up with the whole-repeat ones — and it breaks
    /// it in *both* directions, so "no larger for the partial" is false here too, by a
    /// whisker.
    ///
    /// **Measured over this grid, the two separations differ by at most 0.043 nats — 1.4 parts
    /// in a hundred of a separation of 3.15 nats — and the partial is the larger of the two at
    /// that cell.** That is the honest form of the property, and the assertion below is a
    /// two-sided bound for that reason rather than the one-sided claim the specification
    /// states.
    ///
    /// The unrestricted claim fails far more loudly, and the test after this one is the
    /// counterexample.
    #[test]
    fn a_censored_and_a_complete_read_separate_outgrown_candidates_alike() {
        let epsilon = 1e-3;
        let motif = b"CA";
        let mut widest_gap = WidestDisagreement::default();

        for model in [a_model(), StutterModel::hipstr_shipped()] {
            for witnessed_repeats in [6usize, 8, 10] {
                let observation = a_tract(motif, witnessed_repeats);
                for shorter in 2u32..=4 {
                    for longer in (shorter + 1)..=5 {
                        let short_tract = a_tract(motif, shorter as usize);
                        let long_tract = a_tract(motif, longer as usize);

                        let complete_shorter =
                            score(&observation, &short_tract, shorter, motif, &model, epsilon);
                        let complete_longer =
                            score(&observation, &long_tract, longer, motif, &model, epsilon);
                        let censored_shorter = score_censored(
                            &observation,
                            &short_tract,
                            shorter,
                            motif,
                            &model,
                            epsilon,
                        );
                        let censored_longer = score_censored(
                            &observation,
                            &long_tract,
                            longer,
                            motif,
                            &model,
                            epsilon,
                        );
                        if complete_shorter <= 0.0 || complete_longer <= 0.0 {
                            continue;
                        }

                        let complete_separation =
                            (float::ln(complete_shorter) - float::ln(complete_longer)).abs();
                        let censored_separation =
                            (float::ln(censored_shorter) - float::ln(censored_longer)).abs();
                        let gap = (censored_separation - complete_separation).abs();
                        assert!(
                            gap < 5e-2,
                            "read of {witnessed_repeats} repeats, candidates {shorter} and \
                             {longer}: the censored read separated them by \
                             {censored_separation} nats against the complete read's \
                             {complete_separation}"
                        );
                        if gap > widest_gap.gap {
                            widest_gap = WidestDisagreement {
                                gap,
                                complete_separation,
                                censored_separation,
                            };
                        }
                    }
                }
            }
        }

        // **The grid has to reach a cell where the two genuinely differ**, or the property is
        // being read off arithmetic that happened to cancel — and the cell it reaches is
        // asserted rather than described, so it cannot go stale. At the widest cell the
        // complete read separates the two candidates by 3.149 nats and the partial by 3.193:
        // the partial is the **larger** of the two, by 1.4 parts in a hundred.
        let WidestDisagreement {
            gap,
            complete_separation,
            censored_separation,
        } = widest_gap;
        assert!(
            (gap - 0.043_354).abs() < 1e-5,
            "the two separations differed by at most {gap} nats"
        );
        assert!(
            (complete_separation - 3.149_47).abs() < 1e-4
                && (censored_separation - 3.192_82).abs() < 1e-4,
            "the widest cell moved: complete {complete_separation}, censored \
             {censored_separation}"
        );
    }

    /// **A censored read can separate two candidates further than a complete one, measured
    /// rather than argued** — and it does so whenever the two candidates **straddle** the
    /// stretch the read got through. The reason is not a parameter choice but the shape of the
    /// question a censored read asks.
    ///
    /// A read that got through ten bases of tract, against a candidate that already holds
    /// twelve, needs nothing to have happened at all — the tract is at least ten because it is
    /// twelve — so it collects the same-length share, which is most of the distribution.
    /// Against a candidate holding eight it needs an expansion, which is rare. The complete
    /// read of those same ten bases needs a one-repeat change either way, and those two shares
    /// differ only by the fitted direction ratio. So the censored read separates the two
    /// candidates **further** than the complete read does.
    ///
    /// That is real information rather than a defect: a lower bound rules out everything below
    /// it, and it is exactly the evidence spec §5.1 turned these reads on to collect.
    ///
    /// **This test exists because the specification once claimed the opposite.** §5.2 said a
    /// partial is *always* less discriminating and §12's thirteenth test asked for that without
    /// restriction; both were corrected on 2026-08-25 after this fixture measured them. The
    /// numbers stay asserted here so the correction cannot quietly come undone.
    #[test]
    fn a_censored_read_out_discriminates_a_complete_one_where_the_candidates_straddle_it() {
        // The contraction-biased fitted row, not HipSTR's shipped one: that ships equal
        // direction shares, so the complete read separates the two candidates by exactly
        // zero and the counterexample would be true for an uninteresting reason.
        let model = a_model();
        let epsilon = 1e-3;
        let motif = b"CA";

        let observation = a_tract(motif, 5); // ten bases of tract
        let shorter = a_tract(motif, 4); // eight bases — the read outgrew it
        let longer = a_tract(motif, 6); // twelve bases — the read did not reach its end

        let complete_shorter = score(&observation, &shorter, 4, motif, &model, epsilon);
        let complete_longer = score(&observation, &longer, 6, motif, &model, epsilon);
        let censored_shorter = score_censored(&observation, &shorter, 4, motif, &model, epsilon);
        let censored_longer = score_censored(&observation, &longer, 6, motif, &model, epsilon);

        let complete_separation = (float::ln(complete_shorter) - float::ln(complete_longer)).abs();
        let censored_separation = (float::ln(censored_shorter) - float::ln(censored_longer)).abs();

        assert!(
            censored_separation > complete_separation,
            "the counterexample stopped being one: censored {censored_separation} nats against \
             complete {complete_separation}"
        );
        // **The sizes, asserted rather than described.** The complete read separates the two
        // candidates by 1.586 nats — the direction ratio alone, since one candidate needs a
        // one-repeat expansion and the other a one-repeat contraction. The partial separates
        // them by 5.661, because against the longer candidate it needs nothing to have
        // happened at all.
        assert!(
            (complete_separation - 1.5856).abs() < 1e-3,
            "complete separation {complete_separation}"
        );
        assert!(
            (censored_separation - 5.6605).abs() < 1e-3,
            "censored separation {censored_separation}"
        );
    }

    /// **A read that got through more tract than any stretching of this candidate could
    /// produce scores zero**, and falls to the row's outlier term — the same fate a complete
    /// read of an impossible length meets.
    #[test]
    fn a_censored_read_past_every_stretching_scores_nothing() {
        let model = a_model();
        let motif = b"CA";
        let candidate = a_tract(motif, 4);
        // Eight bases of tract plus the widest whole-repeat expansion, and one repeat more.
        let observation = a_tract(motif, 4 + MAX_WHOLE_REPEAT_SLIP as usize + 1);
        assert_eq!(
            score_censored(&observation, &candidate, 4, motif, &model, 1e-3),
            0.0
        );
    }

    /// **A tract is pure exactly when the segmenter finds one run and nothing else**, which is
    /// the agreement the whole routing rests on: the factorised form is only exact if every
    /// stretching really does share the witnessed prefix, and what builds a stretching is
    /// `segment_tract`.
    ///
    /// Asserted rather than argued, because `segment_tract`'s own documentation calls its
    /// greedy left-to-right phase *a choice* — so a later change to that choice must fail here
    /// rather than silently send interrupted tracts down the exact branch, or worse, pure ones
    /// down neither.
    ///
    /// **The named fixture is the one the length check alone catches.** `CAGCAGTT` under motif
    /// `CAG` passes the chunk walk — `chunks_exact` drops the trailing `TT` on the floor — and
    /// is plainly not a whole number of copies.
    #[test]
    fn a_tract_is_pure_exactly_when_the_segmenter_finds_one_run() {
        let motif = b"CAG";
        let candidate = b"CAGCAGTT";
        assert!(
            candidate
                .chunks_exact(motif.len())
                .all(|chunk| chunk == motif),
            "the fixture must be one the chunk walk alone calls pure, or this tests nothing"
        );
        assert!(
            !is_a_pure_tract(candidate, motif),
            "the length check is what has to catch this"
        );

        // Exhaustive over a two-letter alphabet: motifs to three bases, tracts to eight.
        let mut checked = 0usize;
        for motif_len in 1..=3usize {
            for motif_code in 0..(1usize << motif_len) {
                let motif: Vec<u8> = (0..motif_len)
                    .map(|at| {
                        if (motif_code >> at) & 1 == 0 {
                            b'A'
                        } else {
                            b'C'
                        }
                    })
                    .collect();
                for tract_len in 1..=8usize {
                    for tract_code in 0..(1usize << tract_len) {
                        let tract: Vec<u8> = (0..tract_len)
                            .map(|at| {
                                if (tract_code >> at) & 1 == 0 {
                                    b'A'
                                } else {
                                    b'C'
                                }
                            })
                            .collect();
                        let segments = segment_tract(&tract, &motif);
                        let one_run_only = matches!(segments.as_slice(), [TractSegment::Run(_)]);
                        assert_eq!(
                            is_a_pure_tract(&tract, &motif),
                            one_run_only,
                            "motif {} tract {}: {segments:?}",
                            String::from_utf8_lossy(&motif),
                            String::from_utf8_lossy(&tract)
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert_eq!(checked, 7_140, "the sweep changed size");
    }

    /// **A pure candidate takes the factorised route, and the assertion is bitwise for a
    /// reason.** Both routes are correct to about one part in 10¹⁵, so a tolerance cannot tell
    /// them apart; the product `length factor × letters` is reproduced to the bit only by the
    /// route that computes it that way.
    ///
    /// Without this, forcing every candidate down the exact sum passes the whole suite — the
    /// branch is pinned in one direction only, and the untested direction is the one that
    /// costs a placement enumeration per length change on the majority of candidates.
    #[test]
    fn a_pure_candidate_is_scored_by_the_factorised_route() {
        let model = a_model();
        let epsilon = 1e-3;
        for (motif_bytes, repeat_count) in [(b"A".as_slice(), 1u32), (b"CA".as_slice(), 4)] {
            let motif = a_motif(motif_bytes);
            let candidate = a_tract(motif_bytes, repeat_count as usize);
            assert!(
                is_a_pure_tract(&candidate, motif_bytes),
                "the fixture must be pure"
            );

            for witnessed_len in 0..=(candidate.len() + 8) {
                let witnessed: Vec<u8> = a_tract(motif_bytes, 20)
                    .into_iter()
                    .take(witnessed_len)
                    .collect();
                let scored = score_censored(
                    &witnessed,
                    &candidate,
                    repeat_count,
                    motif_bytes,
                    &model,
                    epsilon,
                );

                let smallest_bp_diff = witnessed.len() as i64 - candidate.len() as i64;
                let length_factor = model.probability_at_least_this_much_longer(
                    smallest_bp_diff,
                    period_of(&motif),
                    repeats(repeat_count),
                );
                let mut stretched = Vec::new();
                resize_at_the_end(&candidate, motif_bytes, smallest_bp_diff, &mut stretched);
                let candidate_view = SsrCandidate {
                    bundle: None,
                    reference_flanks: (&[], &[]),
                    flank_deletion: crate::locus_generation::FlankDeletion::NONE,
                    bases: &candidate,
                    repeat_count: repeats(repeat_count),
                };
                let context = SsrScoringContext::new(
                    &motif,
                    &model,
                    &candidate_view,
                    ErrorRate::try_new(epsilon).expect("a valid rate"),
                    [Provenance::FittedHere],
                );
                let factorised =
                    length_factor * substitution_probability(&witnessed, &stretched, &context);
                assert_eq!(
                    scored.to_bits(),
                    factorised.to_bits(),
                    "motif {} at {witnessed_len} bases witnessed: {scored} against {factorised}",
                    String::from_utf8_lossy(motif_bytes)
                );
            }
        }
    }

    /// **Every stretching of an interrupted tract is cut back to what the read witnessed**, and
    /// that cut is the only thing standing between the exact sum and two sequences of
    /// different lengths — which score nothing at all.
    ///
    /// **The fixture is chosen so that the cut does real work.** Eleven bases witnessed of a
    /// fourteen-base candidate: the deepest admitted stretching lands exactly on the prefix —
    /// the floor is `−3` and the tract gives up exactly three bases — and every other one
    /// overshoots it and must be cut. The bit-for-bit identity test cannot reach this, because
    /// there the read is exactly as long as the single stretching admitted and both cuts are
    /// full length.
    ///
    /// It is also the only test that puts the placement average and the part-repeat arm of the
    /// sum under a real proper prefix.
    #[test]
    fn the_exact_sum_cuts_every_stretching_back_to_the_witnessed_prefix() {
        let model = a_model();
        let candidate = b"CAGCAGTTCAGCAG"; // two runs of two around a TT
        assert!(!is_a_pure_tract(candidate, b"CAG"));
        let witnessed = b"CAGCAGTTCAG";

        let motif = a_motif(b"CAG");
        let smallest_bp_diff = witnessed.len() as i64 - candidate.len() as i64;
        let mut admitted = 0usize;
        let mut part_repeat_terms = 0usize;
        let mut cut_terms = 0usize;
        for (bp_diff, probability) in model.reachable_length_changes(period_of(&motif), repeats(4))
        {
            if bp_diff < smallest_bp_diff || probability <= 0.0 {
                continue;
            }
            admitted += 1;
            if bp_diff % 3 != 0 {
                part_repeat_terms += 1;
            }
            let stretched_len = candidate.len() as i64 + bp_diff;
            assert!(
                stretched_len >= witnessed.len() as i64,
                "stretching by {bp_diff} is shorter than the prefix and was admitted anyway"
            );
            if stretched_len > witnessed.len() as i64 {
                cut_terms += 1;
            }
        }
        assert!(admitted > 1, "the sum has only {admitted} term(s)");
        assert!(
            cut_terms >= admitted - 1,
            "only {cut_terms} of {admitted} terms are actually cut, so the slicing is barely \
             exercised"
        );
        assert!(
            part_repeat_terms > 0,
            "the fixture never enters the part-repeat arm of the sum"
        );

        let scored = score_censored(witnessed, candidate, 4, b"CAG", &model, 1e-3);
        assert!(
            scored > 0.0 && scored < 1.0,
            "the exact sum over cut stretchings came back {scored}"
        );
    }

    /// **A read that witnessed nothing rules nothing out**, so its score is the whole mass the
    /// model can place on this candidate — and on a pure tract that is `1 − unreachable_mass`
    /// to the bit.
    ///
    /// **On an interrupted tract it falls short, and the shortfall is the recorded open
    /// question rather than a defect here.** A slip lands in one run, so two runs of two
    /// repeats cannot lose three even though the tract holds four; the exact sum drops that
    /// term while `unreachable_mass`, which sees only the total, counts it as reachable. The
    /// size is pinned so it cannot grow unnoticed.
    #[test]
    fn a_read_that_witnessed_nothing_scores_the_whole_reachable_mass() {
        let model = a_model();

        let pure = a_tract(b"CA", 4);
        let scored = score_censored(b"", &pure, 4, b"CA", &model, 1e-3);
        let reachable = 1.0 - model.unreachable_mass(period_of(&a_motif(b"CA")), repeats(4));
        assert_eq!(
            scored.to_bits(),
            reachable.to_bits(),
            "{scored} against {reachable}"
        );

        let interrupted = b"CAGCAGTTCAGCAG";
        let scored = score_censored(b"", interrupted, 4, b"CAG", &model, 1e-3);
        let reachable = 1.0 - model.unreachable_mass(period_of(&a_motif(b"CAG")), repeats(4));
        let unplaced_by_the_runs = reachable - scored;
        assert!(
            (unplaced_by_the_runs - 1.3218e-3).abs() < 1e-6,
            "the per-run shortfall moved: {scored} against {reachable}, short by \
             {unplaced_by_the_runs}"
        );
    }

    /// **A candidate holding one repeat can only be stretched.** Both contraction branches
    /// close at once — no whole repeat can go, and none of the non-multiples inside it can
    /// either — so the read is scored over the same-length term and expansions and nothing
    /// else.
    ///
    /// One repeat is the boundary the contraction rule turns on, and no other censored test
    /// uses it.
    #[test]
    fn a_one_repeat_candidate_has_nothing_to_contract() {
        let model = a_model();
        let motif = a_motif(b"CAG");
        let contractions = model
            .reachable_length_changes(period_of(&motif), repeats(1))
            .filter(|(bp_diff, probability)| *bp_diff < 0 && *probability > 0.0)
            .count();
        assert_eq!(contractions, 0, "a one-repeat tract has nothing to give up");

        let scored = score_censored(b"CAG", b"CAG", 1, b"CAG", &model, 1e-3);
        let reachable = 1.0 - model.unreachable_mass(period_of(&motif), repeats(1));
        let letters = float::powi(1.0 - 1e-3, 3);
        assert!(
            (scored - reachable * letters).abs() <= scored * 1e-12,
            "{scored} against {}",
            reachable * letters
        );
    }

    /// **One wrong base costs one `ε/3` in place of one `1 − ε`, and the length factor is
    /// untouched** — the two halves separate on the censored side exactly as they do on the
    /// complete one.
    ///
    /// Every other censored fixture hands the letter half a perfect match, so it is at its
    /// maximum in all of them and this is the only place it is asked to do anything. Without
    /// it, an implementation that compared only the first base, or compared against the wrong
    /// placement where all placements happen to agree, passes.
    #[test]
    fn a_mismatching_base_inside_the_witnessed_prefix_is_charged_once() {
        let model = a_model();
        let epsilon = 1e-3;
        let candidate = a_tract(b"CAG", 6);

        let clean = score_censored(b"CAGCAGCAG", &candidate, 6, b"CAG", &model, epsilon);
        let dirty = score_censored(b"CAGCTGCAG", &candidate, 6, b"CAG", &model, epsilon);
        let expected = clean * (epsilon / 3.0) / (1.0 - epsilon);
        assert!(
            (dirty - expected).abs() <= expected * 1e-12,
            "{dirty} against {expected}"
        );
    }

    /// **The support the distribution hands out and the placements the scoring can build are
    /// the same rule on a pure tract, and the support is the looser of the two otherwise** —
    /// the direction `place_slip_in_each_run`'s own documentation states.
    ///
    /// **This is what stops the two from being changed apart.** The reachability rule is
    /// written on the distribution's side, where it sees only a total repeat count, and here,
    /// where it sees the runs. Before this test, moving either one alone left the other green:
    /// the distribution's tests never call the placements and the placement test pins its
    /// boundary with hand-written literals.
    #[test]
    fn the_placements_and_the_support_agree_on_a_pure_tract() {
        let model = a_model();
        let motif = a_motif(b"CA");
        let period = period_of(&motif);
        let mut placements = Vec::new();

        for (bases, repeat_count, pure) in [
            (b"CACACACA".as_slice(), 4u32, true),
            (b"CACATTCACA".as_slice(), 4u32, false),
        ] {
            assert_eq!(is_a_pure_tract(bases, motif.as_bytes()), pure);
            for (bp_diff, probability) in
                model.reachable_length_changes(period, repeats(repeat_count))
            {
                if probability <= 0.0 || bp_diff % 2 != 0 {
                    continue;
                }
                enumerate_placements(bases, motif.as_bytes(), bp_diff / 2, &mut placements);
                if pure {
                    assert!(
                        !placements.is_empty(),
                        "pure tract: the support offers {bp_diff} and the placements refuse it"
                    );
                }
            }

            // And nothing outside the support is placeable, on either kind of tract: one
            // repeat deeper than the tract can give up, and one repeat past the cutoff.
            for repeats_moved in [
                -(i64::from(repeat_count)),
                i64::from(MAX_WHOLE_REPEAT_SLIP) + 1,
            ] {
                enumerate_placements(bases, motif.as_bytes(), repeats_moved, &mut placements);
                let in_support = model
                    .reachable_length_changes(period, repeats(repeat_count))
                    .any(|(bp_diff, probability)| {
                        bp_diff == repeats_moved * 2 && probability > 0.0
                    });
                assert!(
                    !in_support,
                    "{repeats_moved} repeats is in the support after all"
                );
                if repeats_moved < 0 {
                    assert!(
                        placements.is_empty(),
                        "{repeats_moved} repeats is outside the support but placeable"
                    );
                }
            }
        }
    }

    /// **Where the scoring and the support disagree about what a tract can reach, measured —
    /// because they do, and the disagreement is not this step's to settle.**
    ///
    /// `censored_emission` takes its length changes from
    /// [`StutterModel::reachable_length_changes`], which applies the contraction rule: a read
    /// of this candidate must still show a repeat. `emission` does not — on the part-repeat
    /// branch it asks only whether the distribution gives the change a non-zero probability,
    /// and that question does not know the tract's length. So a two-base `CA` tract scores a
    /// one-base read at 5.390e-4 while `unreachable_mass` counts that same mass as unplaced.
    ///
    /// **All 22 disagreements over this grid run one way**: `emission` places mass the support
    /// refuses, never the reverse, and every one is a part-repeat contraction that would leave
    /// the tract below one repeat. The number is asserted so it cannot grow unnoticed and so
    /// that whoever settles the contraction question in spec §4.2 finds this rather than
    /// rediscovering it.
    ///
    /// **This is a fact about `emission`, which milestone F shipped, not about the censored
    /// term.** Repairing it changes scores at one- and two-repeat tracts, which is a decision
    /// about the model rather than a coding slip.
    #[test]
    fn the_scoring_and_the_support_disagree_only_where_the_open_question_says_they_do() {
        let model = a_model();
        let mut disagreements = Vec::new();

        for motif_bytes in [
            b"A".as_slice(),
            b"CA".as_slice(),
            b"CAG".as_slice(),
            b"CAGT".as_slice(),
        ] {
            let motif = a_motif(motif_bytes);
            for repeat_count in 1u32..=4 {
                let candidate = a_tract(motif_bytes, repeat_count as usize);
                let admitted: Vec<i64> = model
                    .reachable_length_changes(period_of(&motif), repeats(repeat_count))
                    .filter(|(_, probability)| *probability > 0.0)
                    .map(|(bp_diff, _)| bp_diff)
                    .collect();

                for bp_diff in -(candidate.len() as i64)..=40 {
                    let observed_len = candidate.len() as i64 + bp_diff;
                    if observed_len < 0 {
                        continue;
                    }
                    let observation: Vec<u8> = a_tract(motif_bytes, 60)
                        .into_iter()
                        .take(observed_len as usize)
                        .collect();
                    let scored = score(
                        &observation,
                        &candidate,
                        repeat_count,
                        motif_bytes,
                        &model,
                        1e-3,
                    );
                    if (scored > 0.0) != admitted.contains(&bp_diff) {
                        assert!(
                            scored > 0.0,
                            "the support admits {bp_diff} at {repeat_count} repeats of {} and \
                             the scoring refuses it — that is the other direction, and it \
                             would mean the censored term is charging mass nothing can place",
                            String::from_utf8_lossy(motif_bytes)
                        );
                        assert!(
                            bp_diff < 0 && bp_diff % (motif_bytes.len() as i64) != 0,
                            "a disagreement that is not a part-repeat contraction: {bp_diff} at \
                             {repeat_count} repeats of {}",
                            String::from_utf8_lossy(motif_bytes)
                        );
                        disagreements.push(bp_diff);
                    }
                }
            }
        }

        assert_eq!(
            disagreements.len(),
            22,
            "the disagreement between the scoring and the support changed size"
        );
    }

    // ── STR bundles (`str_bundle_loci.md` §7.3) ──

    use crate::region_typing::run_parse::parse_into_runs;

    fn reference_run(offset: u64, length: u64, motif: &[u8]) -> ReferenceRun {
        ReferenceRun::new(offset, length, a_motif(motif)).expect("a run of whole copies")
    }

    /// A second, different model, so a run scored under the wrong one shows.
    fn another_model() -> StutterModel {
        stutter_model_for(&Slippage {
            level: 0.07,
            shorter_share: 0.7,
            fall_off: 0.5,
        })
    }

    /// `P(read | allele)` for an STR bundle allele, with one model per run.
    fn bundle_scored(
        observation: &[u8],
        allele: &[u8],
        hull: &[u8],
        runs: &[ReferenceRun],
        models: &[StutterModel],
        rate: f64,
    ) -> f64 {
        let parse = parse_into_runs(allele, hull, runs);
        let candidate = SsrCandidate {
            bases: allele,
            repeat_count: repeats(1),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            reference_flanks: (&[], &[]),
            bundle: Some(BundleCandidate {
                runs,
                parse: &parse,
                leading_flank_bases: 0,
                reference_hull_len: hull.len() as u32,
            }),
        };
        let motif = a_motif(b"A");
        let fallback = StutterModel::hipstr_shipped();
        let context = SsrScoringContext::for_bundle(
            &motif,
            models.first().unwrap_or(&fallback),
            models,
            &candidate,
            ErrorRate::try_new(rate).expect("a valid rate"),
            [Provenance::FittedHere],
        );
        StutterSubstitutionEmission.emission(
            observation,
            &candidate,
            &context,
            &mut StutterSubstitutionScratch::default(),
        )
    }

    /// `P(read | allele)` for the same bases scored as an STR segment of `motif`.
    fn segment_scored(
        observation: &[u8],
        allele: &[u8],
        motif: &[u8],
        model: &StutterModel,
    ) -> f64 {
        let motif = a_motif(motif);
        let candidate = SsrCandidate {
            bases: allele,
            repeat_count: repeats((allele.len() / motif.period()).max(1) as u32),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            reference_flanks: (&[], &[]),
            bundle: None,
        };
        let context = SsrScoringContext::new(
            &motif,
            model,
            &candidate,
            ErrorRate::try_new(1e-3).expect("a valid rate"),
            [Provenance::FittedHere],
        );
        StutterSubstitutionEmission.emission(
            observation,
            &candidate,
            &context,
            &mut StutterSubstitutionScratch::default(),
        )
    }

    /// **A bundle of one run that spans its whole hull scores every read exactly as the segment
    /// over the same bases does**, bit for bit (`str_bundle_loci.md` §7.3): the regression the
    /// generalisation has to pass. Over three hulls — a pure dinucleotide, one ending on a part
    /// copy, a mononucleotide — every allele a few copies either side, an interrupted one and a
    /// part-unit one, against every read from twelve copies short to twelve long, part-unit
    /// lengths and substitutions included.
    #[test]
    fn a_one_run_bundle_scores_as_the_segment_it_is() {
        let model = a_model();
        let mut compared = 0;
        let mut nonzero = 0;
        for (hull, motif) in [
            (b"CA".repeat(8), &b"CA"[..]),
            ([b"CA".repeat(8), b"C".to_vec()].concat(), &b"CA"[..]),
            (b"A".repeat(12), &b"A"[..]),
            (b"AAG".repeat(6), &b"AAG"[..]),
        ] {
            let runs = [reference_run(0, hull.len() as u64, motif)];
            let period = motif.len();
            let mut alleles: Vec<Vec<u8>> = (0..5)
                .map(|extra| {
                    let copies = hull.len() / period - 2 + extra;
                    [
                        motif.repeat(copies),
                        hull[hull.len() / period * period..].to_vec(),
                    ]
                    .concat()
                })
                .collect();
            let mut interrupted = hull.clone();
            interrupted[period * 3] = b'T';
            alleles.push(interrupted);
            alleles.push([&hull[..], b"G"].concat());
            for allele in &alleles {
                let mut reads: Vec<Vec<u8>> = Vec::new();
                for copies in 1..=allele.len() / period + 12 {
                    let read = motif.repeat(copies);
                    reads.push([&read[..], &hull[hull.len() / period * period..]].concat());
                    reads.push([&read[..], b"T"].concat());
                    reads.push(read);
                }
                let mut substituted = allele.clone();
                substituted[1] = b'G';
                reads.push(substituted);
                reads.push(allele.clone());
                for read in &reads {
                    let as_segment = segment_scored(read, allele, motif, &model);
                    let as_bundle = bundle_scored(
                        read,
                        allele,
                        &hull,
                        &runs,
                        std::slice::from_ref(&model),
                        1e-3,
                    );
                    assert_eq!(
                        as_bundle.to_bits(),
                        as_segment.to_bits(),
                        "read {} against allele {} of hull {}: {as_bundle} as a bundle, \
                         {as_segment} as a segment",
                        String::from_utf8_lossy(read),
                        String::from_utf8_lossy(allele),
                        String::from_utf8_lossy(&hull),
                    );
                    compared += 1;
                    nonzero += usize::from(as_segment > 0.0);
                }
            }
        }
        assert!(
            compared > 1_000 && nonzero > 300,
            "{compared} compared, {nonzero} non-zero"
        );
    }

    /// `(1 − ε)` per matching base and `ε/3` per mismatching one — the substitution factor,
    /// written out independently of `FlatEmission`.
    fn letters_by_hand(read: &[u8], rendered: &[u8], rate: f64) -> f64 {
        assert_eq!(read.len(), rendered.len());
        read.iter()
            .zip(rendered)
            .map(|(a, b)| if a == b { 1.0 - rate } else { rate / 3.0 })
            .product()
    }

    /// **Two runs: the emission is the sum over every slip set reaching the read's length, of the
    /// two runs' probabilities times the letters** — checked against the sum written out by hand.
    /// `GG (CA)8 GT (TTG)5 C` read with three `CA` copies more (`Δ = 6`): the sets are
    /// `(−6, 6), (−3, 4), (0, 2), (3, 0), (6, −2), (9, −4)` — `(−9, 8)` would empty the `CA` run.
    /// Each run is priced under its own model.
    #[test]
    fn a_two_run_bundle_sums_every_slip_set_under_each_runs_own_model() {
        let hull_of = |ca: usize, ttg: usize| {
            [
                &b"GG"[..],
                &b"CA".repeat(ca),
                b"GT",
                &b"TTG".repeat(ttg),
                b"C",
            ]
            .concat()
        };
        let hull = hull_of(8, 5);
        let runs = [reference_run(2, 16, b"CA"), reference_run(20, 15, b"TTG")];
        let models = [a_model(), another_model()];
        let rate = 1e-3;
        let read = hull_of(11, 5);
        let scored = bundle_scored(&read, &hull, &hull, &runs, &models, rate);

        let two = NonZeroU8::new(2).unwrap();
        let three = NonZeroU8::new(3).unwrap();
        let mut expected = 0.0;
        let mut sets = 0;
        for ca_slip in -7_i64..=10 {
            for ttg_slip in -4_i64..=10 {
                if ca_slip * 2 + ttg_slip * 3 != 6 {
                    continue;
                }
                sets += 1;
                let rendered = hull_of((8 + ca_slip) as usize, (5 + ttg_slip) as usize);
                expected += models[0].probability(ca_slip * 2, two)
                    * models[1].probability(ttg_slip * 3, three)
                    * letters_by_hand(&read, &rendered, rate);
            }
        }
        assert_eq!(sets, 6);
        assert!(
            (scored - expected).abs() <= 1e-12 * expected,
            "{scored} against {expected} summed by hand"
        );
        // The read is the `(3, 0)` set exactly; the others add next to nothing to it.
        let main = models[0].probability(6, two)
            * models[1].probability(0, three)
            * letters_by_hand(&read, &read, rate);
        assert!(scored >= main && scored - main <= 1e-6 * main);
    }

    /// **A length change no slip set reaches takes the part-repeat regime**: with two runs of
    /// period 2, one extra base is a sequencing indel in either run, priced by that run's model
    /// and resized at its end while the other run did not slip, the two places weighted equally.
    #[test]
    fn a_change_no_slip_set_reaches_is_a_part_repeat_change_in_one_run() {
        let hull = [&b"CA".repeat(8)[..], b"GG", &b"AT".repeat(6)].concat();
        let runs = [reference_run(0, 16, b"CA"), reference_run(18, 12, b"AT")];
        let models = [a_model(), another_model()];
        let rate = 1e-3;
        let read = [&b"CA".repeat(8)[..], b"C", b"GG", &b"AT".repeat(6)].concat();
        let scored = bundle_scored(&read, &hull, &hull, &runs, &models, rate);

        let two = NonZeroU8::new(2).unwrap();
        let in_the_ca_run = models[0].probability(1, two)
            * models[1].same_length_share()
            * letters_by_hand(&read, &read, rate);
        let at_rendering = [&b"CA".repeat(8)[..], b"GG", &b"AT".repeat(6), b"A"].concat();
        let in_the_at_run = models[1].probability(1, two)
            * models[0].same_length_share()
            * letters_by_hand(&read, &at_rendering, rate);
        let expected = (in_the_ca_run + in_the_at_run) / 2.0;
        assert!(
            (scored - expected).abs() <= 1e-12 * expected,
            "{scored} against {expected}"
        );
    }

    /// **A bundle with no run models no stutter**: a read of the allele's own length is scored on
    /// its letters alone, and a read of any other length scores nothing — the row's outlier term
    /// takes it.
    #[test]
    fn a_bundle_with_no_run_scores_letters_alone_and_no_length_change() {
        let hull = b"ACGTTGCAAGCTAGCATG".to_vec();
        let rate = 1e-3;
        let exact = bundle_scored(&hull, &hull, &hull, &[], &[], rate);
        assert_eq!(exact, letters_by_hand(&hull, &hull, rate));
        let mut substituted = hull.clone();
        substituted[4] = b'A';
        let one_off = bundle_scored(&substituted, &hull, &hull, &[], &[], rate);
        assert!((one_off - letters_by_hand(&substituted, &hull, rate)).abs() < 1e-15);
        assert_eq!(bundle_scored(&hull[1..], &hull, &hull, &[], &[], rate), 0.0);
    }

    /// **A run the allele holds no copy of does not slip.** The `CA` run is gone from this allele,
    /// so a read with one `CA` copy is not stutter of it.
    #[test]
    fn a_run_with_no_copy_in_the_allele_cannot_slip() {
        let hull = [&b"GG"[..], &b"CA".repeat(8), b"GT", &b"TTG".repeat(5), b"C"].concat();
        let runs = [reference_run(2, 16, b"CA"), reference_run(20, 15, b"TTG")];
        let models = [a_model(), another_model()];
        let without_ca = [&b"GG"[..], b"GT", &b"TTG".repeat(5), b"C"].concat();
        let one_ca = [&b"GG"[..], b"CA", b"GT", &b"TTG".repeat(5), b"C"].concat();
        let scored = bundle_scored(&one_ca, &without_ca, &hull, &runs, &models, 1e-3);
        // Not the `CA` run gaining a copy; only the part-repeat regime of the `TTG` run, which
        // a two-base change reaches as a sequencing indel there.
        let three = NonZeroU8::new(3).unwrap();
        // Exactly the `TTG` run's part-repeat price, its stretch resized at its end, with the
        // copyless `CA` run a factor of one and no place of its own.
        let ttg = 4..4 + 15;
        let mut resized = Vec::new();
        resize_at_the_end(&without_ca[ttg.clone()], b"TTG", 2, &mut resized);
        let rendered = [&without_ca[..ttg.start], &resized, &without_ca[ttg.end..]].concat();
        let expected = models[1].probability(2, three) * letters_by_hand(&one_ca, &rendered, 1e-3);
        assert!(
            (scored - expected).abs() <= 1e-12 * expected,
            "{scored} against {expected}"
        );
    }

    /// **A run the allele holds no copy of is a factor of one**, not its model's no-slip share:
    /// a read of the allele's own bases scores the other run's no-slip share times the letters,
    /// as it would were the copyless run never in the hull (found in review, 2026-10-11: it was
    /// scaled by the copyless run's 0.88 on every read).
    #[test]
    fn a_run_with_no_copy_contributes_a_factor_of_one() {
        let hull = [&b"GG"[..], &b"CA".repeat(8), b"GT", &b"TTG".repeat(5), b"C"].concat();
        let runs = [reference_run(2, 16, b"CA"), reference_run(20, 15, b"TTG")];
        let models = [a_model(), another_model()];
        let without_ca = [&b"GG"[..], b"GT", &b"TTG".repeat(5), b"C"].concat();
        let scored = bundle_scored(&without_ca, &without_ca, &hull, &runs, &models, 1e-3);
        let three = NonZeroU8::new(3).unwrap();
        let expected =
            models[1].probability(0, three) * letters_by_hand(&without_ca, &without_ca, 1e-3);
        assert!(
            (scored - expected).abs() <= 1e-12 * expected,
            "{scored} against {expected}"
        );
    }

    /// **A run too short for a part-repeat change, or with no copy, is not a place for it.**
    /// `GG GT (AAAG)4 C` has lost its `CA` run; a read one base short (`Δ = −1`, reached by no
    /// slip set of periods 2 and 4) is priced in the `AAAG` run alone. The `CA` run was a place
    /// too, rendering an allele one base longer than the read: a panic in a debug build, and in
    /// release half the read's value (found in review, 2026-10-11).
    #[test]
    fn a_run_with_no_copy_is_not_a_place_for_a_part_repeat_change() {
        let hull = [
            &b"GG"[..],
            &b"CA".repeat(8),
            b"GT",
            &b"AAAG".repeat(4),
            b"C",
        ]
        .concat();
        let runs = [reference_run(2, 16, b"CA"), reference_run(20, 16, b"AAAG")];
        let models = [a_model(), another_model()];
        let allele = [&b"GG"[..], b"GT", &b"AAAG".repeat(4), b"C"].concat();
        let read = [&b"GG"[..], b"GT", &b"AAAG".repeat(3), b"AAA", b"C"].concat();
        let scored = bundle_scored(&read, &allele, &hull, &runs, &models, 1e-3);
        let four = NonZeroU8::new(4).unwrap();
        let expected = models[1].probability(-1, four) * letters_by_hand(&read, &read, 1e-3);
        assert!(
            (scored - expected).abs() <= 1e-12 * expected,
            "{scored} against {expected}"
        );
    }

    /// **The lengths a bundle candidate reaches are its slip sets' totals and its runs'
    /// part-repeat changes**; for one run spanning the hull they are the segment's own.
    #[test]
    fn a_bundles_reachable_lengths_are_its_slip_set_totals_and_part_changes() {
        let hull = b"CA".repeat(8);
        let runs = [reference_run(0, 16, b"CA")];
        let parse = parse_into_runs(&hull, &hull, &runs);
        let as_bundle = SsrCandidate {
            bases: &hull,
            repeat_count: repeats(8),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            reference_flanks: (&[], &[]),
            bundle: Some(BundleCandidate {
                runs: &runs,
                parse: &parse,
                leading_flank_bases: 0,
                reference_hull_len: hull.len() as u32,
            }),
        };
        let as_segment = SsrCandidate {
            bundle: None,
            ..as_bundle
        };
        let (mut from_bundle, mut from_segment) = (Vec::new(), Vec::new());
        fill_reachable_lengths(&[as_bundle], &a_motif(b"CA"), &mut from_bundle);
        fill_reachable_lengths(&[as_segment], &a_motif(b"CA"), &mut from_segment);
        assert_eq!(from_bundle, from_segment);

        let hull = [&b"CA".repeat(8)[..], b"GT", &b"TTG".repeat(5)].concat();
        let runs = [reference_run(0, 16, b"CA"), reference_run(18, 15, b"TTG")];
        let parse = parse_into_runs(&hull, &hull, &runs);
        let candidate = SsrCandidate {
            bases: &hull,
            bundle: Some(BundleCandidate {
                runs: &runs,
                parse: &parse,
                leading_flank_bases: 0,
                reference_hull_len: hull.len() as u32,
            }),
            ..as_segment
        };
        let mut reached = Vec::new();
        fill_reachable_lengths(&[candidate], &a_motif(b"A"), &mut reached);
        let held = hull.len() as u32;
        // One base longer is no slip set of whole units of 2 and 3 alone — it is (2, −1).
        for change in [1_u32, 5, 6, 30] {
            assert!(reached.contains(&(held + change)), "{change} longer");
        }
        // The `CA` run can lose at most 7 copies and the `TTG` run 4: 26 bases shorter, no more.
        assert!(reached.contains(&(held - 26)));
        assert!(!reached.contains(&(held - 27)) && !reached.contains(&(held - 28)));
    }

    /// **Bases the scored sequence carries in front of the allele move the runs with them**: a
    /// read and an allele both given two flank bases score as the two without them, times the two
    /// matching bases' letters.
    #[test]
    fn leading_flank_bases_shift_every_run_with_the_allele() {
        let hull = [&b"GG"[..], &b"CA".repeat(8), b"GT", &b"TTG".repeat(5), b"C"].concat();
        let runs = [reference_run(2, 16, b"CA"), reference_run(20, 15, b"TTG")];
        let models = [a_model(), another_model()];
        let rate = 1e-3;
        let read = [&b"GG"[..], &b"CA".repeat(9), b"GT", &b"TTG".repeat(5), b"C"].concat();
        let plain = bundle_scored(&read, &hull, &hull, &runs, &models, rate);

        let parse = parse_into_runs(&hull, &hull, &runs);
        let given = [&b"TA"[..], &hull].concat();
        let read_given = [&b"TA"[..], &read].concat();
        let candidate = SsrCandidate {
            bases: &given,
            repeat_count: repeats(1),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            reference_flanks: (&[], &[]),
            bundle: Some(BundleCandidate {
                runs: &runs,
                parse: &parse,
                leading_flank_bases: 2,
                reference_hull_len: hull.len() as u32,
            }),
        };
        let motif = a_motif(b"A");
        let context = SsrScoringContext::for_bundle(
            &motif,
            &models[0],
            &models,
            &candidate,
            ErrorRate::try_new(rate).unwrap(),
            [Provenance::FittedHere],
        );
        let shifted = StutterSubstitutionEmission.emission(
            &read_given,
            &candidate,
            &context,
            &mut StutterSubstitutionScratch::default(),
        );
        let expected = plain * (1.0 - rate) * (1.0 - rate);
        assert!(
            (shifted - expected).abs() <= 1e-12 * expected,
            "{shifted} against {expected}"
        );
    }

    /// **Comparing from the first mismatch on gives `letters_under`'s result bit for bit**: every
    /// length from no base to 40, every place of the first mismatch, mismatches after it at
    /// pseudo-random places, under three rates — and again after the scores change, which must
    /// rebuild the kept sums rather than reuse them.
    #[test]
    fn a_matched_prefix_gives_letters_under_bit_for_bit() {
        let mut matched = MatchedPrefixLetters::default();
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 33
        };
        let mut compared = 0;
        for rate in [1e-3, 0.02, 0.3, 1e-3] {
            let scores = FlatEmission::try_new(rate)
                .expect("a valid rate")
                .scores_for(BaseQual(0));
            matched.prepare(scores, 40);
            for length in 0..=40_usize {
                let read: Vec<u8> = (0..length).map(|_| b"ACGT"[next() as usize % 4]).collect();
                for first_mismatch in 0..=length {
                    let mut other = read.clone();
                    if first_mismatch < length {
                        other[first_mismatch] = if read[first_mismatch] == b'A' {
                            b'C'
                        } else {
                            b'A'
                        };
                        for base in &mut other[first_mismatch + 1..] {
                            if next() % 3 == 0 {
                                *base = b"ACGT"[next() as usize % 4];
                            }
                        }
                    }
                    assert_eq!(
                        matched.letters(&read, &other).to_bits(),
                        letters_under(scores, &read, &other).to_bits(),
                        "rate {rate}, {length} bases, first mismatch at {first_mismatch}"
                    );
                    compared += 1;
                }
            }
        }
        assert_eq!(compared, 4 * (41 * 42 / 2));
        assert_eq!(agreeing_prefix_len(b"ACGTACGTACGTA", b"ACGTACGTACGTA"), 13);
        assert_eq!(agreeing_prefix_len(b"ACGTACGTACGTA", b"ACGTACGTAAGTA"), 9);
        assert_eq!(agreeing_prefix_len(b"ACGTACGTACGTA", b"ACGTACCTACGTA"), 6);
    }

    // ── STR bundles: a length change no single run makes (owner, 2026-10-12) ──

    /// **A one-base change no run can make alone is priced as a sequencing indel as well as by
    /// the runs slipping opposite ways** (the review's M3, settled by the owner 2026-10-12). At
    /// `GG (CA)8 GT (TTG)5 C` one base more is `(−1, 1)`, `(2, −1)` and the other sets of periods 2
    /// and 3 — and also a one-base indel in either run, which is what it usually is. Before, the
    /// slip sets alone priced it.
    ///
    /// **Two reads, so that both halves of the sum are pinned** (str_bundles_partial_reads review,
    /// 2026-10-12: with the first read alone the slip sets were 1 part in 10¹⁷ of the total, and an
    /// emission that priced only the indel passed). The first shows the `CA` run's own one-base
    /// indel, which the indel term explains; the second shows the `(−1, +1)` slip set, which only
    /// the slip sets explain.
    #[test]
    fn a_change_no_single_run_makes_is_also_priced_as_a_sequencing_indel() {
        let hull = two_run_hull(8, 5);
        let runs = two_runs();
        let models = [a_model(), another_model()];
        let rate = 1e-3;
        let (two, three) = (NonZeroU8::new(2).unwrap(), NonZeroU8::new(3).unwrap());
        let ca_resized = [
            &b"GG"[..],
            &b"CA".repeat(8),
            b"C",
            b"GT",
            &b"TTG".repeat(5),
            b"C",
        ]
        .concat();
        let ttg_resized = [
            &b"GG"[..],
            &b"CA".repeat(8),
            b"GT",
            &b"TTG".repeat(5),
            b"T",
            b"C",
        ]
        .concat();
        let one_slip_set = two_run_hull(7, 6);
        for (read, the_indel_explains_it) in [(&ca_resized, true), (&one_slip_set, false)] {
            let scored = bundle_scored(read, &hull, &hull, &runs, &models, rate);
            let mut slip_sets = 0.0;
            for ca_slip in -7_i64..=10 {
                for ttg_slip in -4_i64..=10 {
                    if ca_slip * 2 + ttg_slip * 3 == 1 {
                        slip_sets += models[0].probability(ca_slip * 2, two)
                            * models[1].probability(ttg_slip * 3, three)
                            * letters_by_hand(
                                read,
                                &two_run_hull_counts(8 + ca_slip, 5 + ttg_slip),
                                rate,
                            );
                    }
                }
            }
            let as_an_indel = (models[0].probability(1, two)
                * models[1].same_length_share()
                * letters_by_hand(read, &ca_resized, rate)
                + models[1].probability(1, three)
                    * models[0].same_length_share()
                    * letters_by_hand(read, &ttg_resized, rate))
                / 2.0;
            let expected = slip_sets + as_an_indel;
            assert!(
                (scored - expected).abs() <= 1e-12 * expected,
                "{scored} against {expected}"
            );
            // Each read is explained by one half of the sum alone, so leaving that half out would
            // move the score by far more than the tolerance.
            let (explaining, other) = if the_indel_explains_it {
                (as_an_indel, slip_sets)
            } else {
                (slip_sets, as_an_indel)
            };
            assert!(
                explaining > 1e6 * other,
                "{explaining} explains the read, {other} does not"
            );
        }
    }

    /// **A change one run makes alone is not also an indel**: one base more at a bundle with a
    /// homopolymer run is that run slipping, as at a homopolymer segment.
    ///
    /// **The read is the rendering a wrongly added indel term would explain** — the `TTG` run one
    /// base longer at its end — so that such a term would dominate the score (str_bundles_partial_reads
    /// review, 2026-10-12: against a read of the `A` run's slip it was 5 parts in 10⁴³ of the total,
    /// and the test could not fail).
    #[test]
    fn a_change_one_run_makes_alone_takes_no_indel_term() {
        let hull = [&b"GG"[..], &b"A".repeat(10), b"GT", &b"TTG".repeat(5), b"C"].concat();
        let runs = [reference_run(2, 10, b"A"), reference_run(14, 15, b"TTG")];
        let models = [a_model(), another_model()];
        let rate = 1e-3;
        let read = [
            &b"GG"[..],
            &b"A".repeat(10),
            b"GT",
            &b"TTG".repeat(5),
            b"T",
            b"C",
        ]
        .concat();
        let wrong_indel_term = models[1].probability(1, NonZeroU8::new(3).unwrap())
            * models[0].same_length_share()
            * letters_by_hand(&read, &read, rate);
        let scored = bundle_scored(&read, &hull, &hull, &runs, &models, rate);
        let (one, three) = (NonZeroU8::new(1).unwrap(), NonZeroU8::new(3).unwrap());
        let mut expected = 0.0;
        for a_slip in -9_i64..=10 {
            for ttg_slip in -4_i64..=10 {
                if a_slip + ttg_slip * 3 == 1 {
                    let rendered = [
                        &b"GG"[..],
                        &b"A".repeat((10 + a_slip) as usize),
                        b"GT",
                        &b"TTG".repeat((5 + ttg_slip) as usize),
                        b"C",
                    ]
                    .concat();
                    expected += models[0].probability(a_slip, one)
                        * models[1].probability(ttg_slip * 3, three)
                        * letters_by_hand(&read, &rendered, rate);
                }
            }
        }
        assert!(
            (scored - expected).abs() <= 1e-12 * expected,
            "{scored} against {expected}"
        );
        assert!(
            wrong_indel_term > 1e6 * expected,
            "a wrongly added indel term ({wrong_indel_term}) would swamp the price ({expected})"
        );
    }

    // ── STR bundles: reads that end inside the hull (`str_bundle_loci.md` §7.6) ──

    /// The reference flanks the partial-read tests carry.
    const LEFT_FLANK: &[u8] = b"ACGTACGTAC";
    const RIGHT_FLANK: &[u8] = b"GTCAGTCAGT";

    /// `GG (CA)ca GT (TTG)ttg C`: two runs, of periods 2 and 3, fixed bases around and between.
    fn two_run_hull(ca: usize, ttg: usize) -> Vec<u8> {
        [
            &b"GG"[..],
            &b"CA".repeat(ca),
            b"GT",
            &b"TTG".repeat(ttg),
            b"C",
        ]
        .concat()
    }

    fn two_run_hull_counts(ca: i64, ttg: i64) -> Vec<u8> {
        two_run_hull(ca as usize, ttg as usize)
    }

    /// The runs of `two_run_hull(8, 5)`.
    fn two_runs() -> [ReferenceRun; 2] {
        [reference_run(2, 16, b"CA"), reference_run(20, 15, b"TTG")]
    }

    fn witness(positions: (u16, u16)) -> WitnessedLocusPositions {
        WitnessedLocusPositions::from_half_open_runs([positions]).expect("a witness")
    }

    /// `P(read that ended inside the hull | allele)` for an STR bundle allele, at a substitution
    /// rate of 1 in 1,000.
    fn bundle_censored_scored(
        witnessed: &[u8],
        positions: &WitnessedLocusPositions,
        allele: &[u8],
        hull: &[u8],
        runs: &[ReferenceRun],
        models: &[StutterModel],
        flanks: (&[u8], &[u8]),
    ) -> f64 {
        let parse = parse_into_runs(allele, hull, runs);
        let candidate = SsrCandidate {
            bases: allele,
            repeat_count: repeats(1),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            reference_flanks: flanks,
            bundle: Some(BundleCandidate {
                runs,
                parse: &parse,
                leading_flank_bases: 0,
                reference_hull_len: hull.len() as u32,
            }),
        };
        let motif = a_motif(b"A");
        let fallback = StutterModel::hipstr_shipped();
        let context = SsrScoringContext::for_bundle(
            &motif,
            models.first().unwrap_or(&fallback),
            models,
            &candidate,
            ErrorRate::try_new(1e-3).expect("a valid rate"),
            [Provenance::FittedHere],
        );
        StutterSubstitutionEmission.censored_emission(
            witnessed,
            positions,
            &candidate,
            &context,
            &mut StutterSubstitutionScratch::default(),
        )
    }

    /// **A one-run bundle over its hull, with no flank bases to read on into, scores a read that
    /// ended inside it as the segment over the same bases does** — the segment's censored sum,
    /// over pure, interrupted and mononucleotide tracts, every read length from one copy (and two
    /// bases at least) to past the allele's end. (Equal to twelve digits rather than bit for bit:
    /// the segment factors the letters out of a pure tract's sum.)
    #[test]
    fn a_one_run_bundles_partial_read_scores_as_the_segments() {
        let model = a_model();
        let mut compared = 0;
        let mut nonzero = 0;
        for (hull, motif) in [
            (b"CA".repeat(8), &b"CA"[..]),
            (b"A".repeat(12), &b"A"[..]),
            (b"AAG".repeat(6), &b"AAG"[..]),
        ] {
            let runs = [reference_run(0, hull.len() as u64, motif)];
            let period = motif.len();
            let mut alleles: Vec<Vec<u8>> = (0..5)
                .map(|extra| motif.repeat(hull.len() / period - 2 + extra))
                .collect();
            let mut interrupted = hull.clone();
            interrupted[period * 3] = b'T';
            alleles.push(interrupted);
            for allele in &alleles {
                let copies = (allele.len() / period) as u32;
                // From one copy on: a rendering shorter than a copy is where the segment, which stops
                // its part-repeat contractions earlier than its complete reads do, sums less.
                for shown in period.max(2)..allele.len() + 3 * period {
                    let mut reads: Vec<Vec<u8>> =
                        vec![motif.iter().copied().cycle().take(shown).collect()];
                    reads.push(
                        allele
                            .iter()
                            .chain(motif)
                            .copied()
                            .cycle()
                            .take(shown)
                            .collect(),
                    );
                    for read in &reads {
                        let as_segment = score_censored(read, allele, copies, motif, &model, 1e-3);
                        let as_bundle = bundle_censored_scored(
                            read,
                            &witness((0, 1)),
                            allele,
                            &hull,
                            &runs,
                            std::slice::from_ref(&model),
                            (&[], &[]),
                        );
                        assert!(
                            (as_bundle - as_segment).abs() <= 1e-12 * as_segment.max(as_bundle),
                            "read {} against allele {}: {as_bundle} as a bundle, {as_segment} as a \
                             segment",
                            String::from_utf8_lossy(read),
                            String::from_utf8_lossy(allele),
                        );
                        compared += 1;
                        nonzero += usize::from(as_segment > 1e-6);
                    }
                }
            }
        }
        assert!(
            compared > 300 && nonzero > 100,
            "{compared} compared, {nonzero} non-zero"
        );
    }

    /// What a read that ended inside `two_run_hull(ca, ttg)` sums, written out by hand: every slip
    /// set, each run's whole range, and every change no single run makes as a sequencing indel in
    /// each run it can be one in — each rendering priced by `letters`.
    fn two_run_partial_by_hand(
        ca: i64,
        ttg: i64,
        models: &[StutterModel; 2],
        letters: impl Fn(&[u8]) -> f64,
    ) -> f64 {
        let (two, three) = (NonZeroU8::new(2).unwrap(), NonZeroU8::new(3).unwrap());
        let (ca_lowest, ttg_lowest) = ((1 - ca).max(-10), (1 - ttg).max(-10));
        let mut total = 0.0;
        for ca_slip in ca_lowest..=10 {
            for ttg_slip in ttg_lowest..=10 {
                total += models[0].probability(ca_slip * 2, two)
                    * models[1].probability(ttg_slip * 3, three)
                    * letters(&two_run_hull_counts(ca + ca_slip, ttg + ttg_slip));
            }
        }
        let allele = two_run_hull_counts(ca, ttg);
        let ca_stretch = 2..2 + 2 * ca as usize;
        let ttg_stretch = ca_stretch.end + 2..ca_stretch.end + 2 + 3 * ttg as usize;
        for bp_diff in -20_i64..=20 {
            let one_run_makes_it = bp_diff == 0
                || (bp_diff % 2 == 0 && (ca_lowest..=10).contains(&(bp_diff / 2)))
                || (bp_diff % 3 == 0 && (ttg_lowest..=10).contains(&(bp_diff / 3)));
            if one_run_makes_it {
                continue;
            }
            let (mut places, mut sum) = (0.0, 0.0);
            for (stretch, motif, own, other, period) in [
                (ca_stretch.clone(), &b"CA"[..], &models[0], &models[1], two),
                (
                    ttg_stretch.clone(),
                    &b"TTG"[..],
                    &models[1],
                    &models[0],
                    three,
                ),
            ] {
                if bp_diff % i64::from(period.get()) == 0 || stretch.len() as i64 + bp_diff < 0 {
                    continue;
                }
                places += 1.0;
                let mut resized = Vec::new();
                resize_at_the_end(&allele[stretch.clone()], motif, bp_diff, &mut resized);
                let rendered =
                    [&allele[..stretch.start], &resized, &allele[stretch.end..]].concat();
                sum += own.probability(bp_diff, period)
                    * other.same_length_share()
                    * letters(&rendered);
            }
            if places > 0.0 {
                total += sum / places;
            }
        }
        total
    }

    /// `letters` for a read that held the left flank: its bases against the first ones of the
    /// rendering followed by the right flank.
    fn from_the_left<'a>(read: &'a [u8]) -> impl Fn(&[u8]) -> f64 + 'a {
        move |rendered| {
            let carrier = [rendered, RIGHT_FLANK].concat();
            if carrier.len() < read.len() {
                return 0.0;
            }
            letters_by_hand(read, &carrier[..read.len()], 1e-3)
        }
    }

    /// `letters` for a read that held the right flank: its bases against the last ones of the left
    /// flank followed by the rendering.
    fn from_the_right<'a>(read: &'a [u8]) -> impl Fn(&[u8]) -> f64 + 'a {
        move |rendered| {
            let carrier = [LEFT_FLANK, rendered].concat();
            if carrier.len() < read.len() {
                return 0.0;
            }
            letters_by_hand(read, &carrier[carrier.len() - read.len()..], 1e-3)
        }
    }

    /// **A read that ended inside a two-run hull sums every slip set and every sequencing indel,
    /// compared from the flank it held** — checked against that sum written out by hand, for a read
    /// that held the left flank, one that held the right, and one whose witness covers the whole
    /// hull (either, one half each).
    #[test]
    fn a_partial_read_sums_every_rendering_from_the_flank_it_held() {
        let hull = two_run_hull(8, 5);
        let runs = two_runs();
        let models = [a_model(), another_model()];
        let flanks = (LEFT_FLANK, RIGHT_FLANK);
        let close = |scored: f64, expected: f64| {
            assert!(
                (scored - expected).abs() <= 1e-12 * expected,
                "{scored} against {expected} by hand"
            );
        };
        // Held the left flank; a `CA` copy more than the reference, then into the `TTG` run.
        let from_left = two_run_hull(9, 5)[..27].to_vec();
        for (ca, ttg) in [(8, 5), (9, 5), (9, 6), (9, 1)] {
            let allele = two_run_hull_counts(ca, ttg);
            close(
                bundle_censored_scored(
                    &from_left,
                    &witness((0, 27)),
                    &allele,
                    &hull,
                    &runs,
                    &models,
                    flanks,
                ),
                two_run_partial_by_hand(ca, ttg, &models, from_the_left(&from_left)),
            );
        }
        // Held the right flank; began inside the `CA` run of an allele with a `TTG` copy more.
        let from_right = two_run_hull(8, 6)[10..].to_vec();
        for (ca, ttg) in [(8, 5), (8, 6), (3, 6)] {
            let allele = two_run_hull_counts(ca, ttg);
            close(
                bundle_censored_scored(
                    &from_right,
                    &witness((7, 36)),
                    &allele,
                    &hull,
                    &runs,
                    &models,
                    flanks,
                ),
                two_run_partial_by_hand(ca, ttg, &models, from_the_right(&from_right)),
            );
        }
        // Crossed the whole hull and three bases of the right flank without anchoring there.
        let either = [&hull[..], &RIGHT_FLANK[..3]].concat();
        let allele = two_run_hull(8, 5);
        close(
            bundle_censored_scored(
                &either,
                &witness((0, 36)),
                &allele,
                &hull,
                &runs,
                &models,
                flanks,
            ),
            0.5 * (two_run_partial_by_hand(8, 5, &models, from_the_left(&either))
                + two_run_partial_by_hand(8, 5, &models, from_the_right(&either))),
        );
    }

    /// **What a read that ended inside the hull says, in plain terms**: it separates alleles that
    /// differ in a run it got past, as a complete read does; it is a lower bound on the run it
    /// ended inside, so every allele long enough there agrees with it; and it rules out an allele
    /// shorter there than what it showed.
    #[test]
    fn a_partial_read_is_a_lower_bound_on_the_run_it_ended_in() {
        let hull = two_run_hull(8, 5);
        let runs = two_runs();
        let models = [a_model(), another_model()];
        let flanks = (LEFT_FLANK, RIGHT_FLANK);
        // `GG (CA)9 GT TTGTT` — got past the `CA` run, ended five bases into the `TTG` run.
        let read = two_run_hull(9, 5)[..27].to_vec();
        let scored = |ca: usize, ttg: usize| {
            bundle_censored_scored(
                &read,
                &witness((0, 27)),
                &two_run_hull(ca, ttg),
                &hull,
                &runs,
                &models,
                flanks,
            )
        };
        let carrier = scored(9, 5);
        // Every `TTG` count from two copies up holds the five bases it showed; what separates them
        // is only how likely each is to have contracted below two.
        for ttg in [6, 9] {
            let ratio = scored(9, ttg) / carrier;
            assert!(
                (ratio - 1.0).abs() < 1e-2,
                "(9, {ttg}): {ratio} of the carrier's"
            );
        }
        // One `TTG` copy is shorter than what it showed; one `CA` copy fewer is a run it got past.
        assert!(
            scored(9, 1) < 0.05 * carrier,
            "{} against {carrier}",
            scored(9, 1)
        );
        assert!(
            scored(8, 5) < 0.01 * carrier,
            "{} against {carrier}",
            scored(8, 5)
        );
    }

    /// What a read that ended inside a bundle of **pure** runs sums, written out by hand for any
    /// number of runs: every slip set (each run's whole range, by an odometer over the runs), and
    /// every change no single run makes as a sequencing indel in each run it can be one in — each
    /// rendering priced by `letters`. `fixed` holds the fixed pieces around and between the runs,
    /// one more than `runs`; each run is its motif and its copies in the allele (at least one).
    fn partial_by_hand(
        fixed: &[&[u8]],
        runs: &[(&[u8], i64)],
        models: &[StutterModel],
        letters: impl Fn(&[u8]) -> f64,
    ) -> f64 {
        assert_eq!(fixed.len(), runs.len() + 1);
        let period = |motif: &[u8]| NonZeroU8::new(motif.len() as u8).unwrap();
        let render = |counts: &[i64]| -> Vec<u8> {
            let mut out = fixed[0].to_vec();
            for (run, &(motif, _)) in runs.iter().enumerate() {
                out.extend(motif.repeat(counts[run] as usize));
                out.extend_from_slice(fixed[run + 1]);
            }
            out
        };
        let lowest: Vec<i64> = runs
            .iter()
            .map(|&(_, copies)| (1 - copies).max(-10))
            .collect();
        let mut total = 0.0;
        let mut slips = lowest.clone();
        'sets: loop {
            let mut probability = 1.0;
            for (run, &(motif, _)) in runs.iter().enumerate() {
                probability *=
                    models[run].probability(slips[run] * motif.len() as i64, period(motif));
            }
            let counts: Vec<i64> = runs
                .iter()
                .zip(&slips)
                .map(|(&(_, copies), &slip)| copies + slip)
                .collect();
            total += probability * letters(&render(&counts));
            for run in 0..runs.len() {
                if slips[run] < 10 {
                    slips[run] += 1;
                    continue 'sets;
                }
                slips[run] = lowest[run];
            }
            break;
        }
        let copies: Vec<i64> = runs.iter().map(|&(_, copies)| copies).collect();
        let allele = render(&copies);
        for bp_diff in -20_i64..=20 {
            let one_run_makes_it = bp_diff == 0
                || runs.iter().zip(&lowest).any(|(&(motif, _), &lowest)| {
                    let period = motif.len() as i64;
                    bp_diff % period == 0 && (lowest..=10).contains(&(bp_diff / period))
                });
            if one_run_makes_it {
                continue;
            }
            let (mut places, mut sum) = (0.0, 0.0);
            let mut start = fixed[0].len();
            for (run, &(motif, copies)) in runs.iter().enumerate() {
                let stretch = start..start + motif.len() * copies as usize;
                start = stretch.end + fixed[run + 1].len();
                if bp_diff % motif.len() as i64 == 0 || stretch.len() as i64 + bp_diff < 0 {
                    continue;
                }
                places += 1.0;
                let mut price = models[run].probability(bp_diff, period(motif));
                for (other, model) in models.iter().enumerate() {
                    if other != run {
                        price *= model.same_length_share();
                    }
                }
                let mut resized = Vec::new();
                resize_at_the_end(&allele[stretch.clone()], motif, bp_diff, &mut resized);
                let rendered =
                    [&allele[..stretch.start], &resized, &allele[stretch.end..]].concat();
                sum += price * letters(&rendered);
            }
            if places > 0.0 {
                total += sum / places;
            }
        }
        total
    }

    /// `letters` for a read that held the left flank, against `rendered` then `right_flank`.
    fn held_left<'a>(read: &'a [u8], right_flank: &'a [u8]) -> impl Fn(&[u8]) -> f64 + 'a {
        move |rendered| {
            let carrier = [rendered, right_flank].concat();
            if carrier.len() < read.len() {
                return 0.0;
            }
            letters_by_hand(read, &carrier[..read.len()], 1e-3)
        }
    }

    /// `letters` for a read that held the right flank, against `left_flank` then `rendered`.
    fn held_right<'a>(read: &'a [u8], left_flank: &'a [u8]) -> impl Fn(&[u8]) -> f64 + 'a {
        move |rendered| {
            let carrier = [left_flank, rendered].concat();
            if carrier.len() < read.len() {
                return 0.0;
            }
            letters_by_hand(read, &carrier[carrier.len() - read.len()..], 1e-3)
        }
    }

    /// **A read that ends inside fixed sequence, exactly where a run or a fixed piece ends, or
    /// that is longer than any carrier sums as the rule says** — from each flank, against alleles
    /// longer, shorter and equal in each run (str_bundles_partial_reads review, 2026-10-12: the
    /// first two-run test had reads ending inside a run and in the flank only).
    #[test]
    fn a_partial_read_ending_in_fixed_sequence_or_at_a_boundary_sums_every_rendering() {
        let hull = two_run_hull(8, 5);
        let runs = two_runs();
        let models = [a_model(), another_model()];
        let flanks = (LEFT_FLANK, RIGHT_FLANK);
        // `two_run_hull(9, 5)`: `CA` run 2..20, `GT` 20..22, `TTG` run 22..37, `C` 37.
        let carrier = two_run_hull(9, 5);
        let too_long = [&carrier[..], RIGHT_FLANK, &b"A".repeat(70)].concat();
        let too_long_right = [&b"A".repeat(70)[..], LEFT_FLANK, &carrier[..]].concat();
        let mut compared = 0;
        for (ca, ttg) in [(8, 5), (9, 1), (3, 6)] {
            let allele = two_run_hull_counts(ca, ttg);
            let fixed: [&[u8]; 3] = [b"GG", b"GT", b"C"];
            let pure_runs = [(&b"CA"[..], ca), (&b"TTG"[..], ttg)];
            for end in [21, 20, 22] {
                let read = &carrier[..end];
                let scored = bundle_censored_scored(
                    read,
                    &witness((0, 30)),
                    &allele,
                    &hull,
                    &runs,
                    &models,
                    flanks,
                );
                let expected =
                    partial_by_hand(&fixed, &pure_runs, &models, held_left(read, RIGHT_FLANK));
                assert!(
                    (scored - expected).abs() <= 1e-12 * expected && expected > 0.0,
                    "({ca}, {ttg}) from the left to {end}: {scored} against {expected}"
                );
                let read = &carrier[end..];
                let scored = bundle_censored_scored(
                    read,
                    &witness((10, 36)),
                    &allele,
                    &hull,
                    &runs,
                    &models,
                    flanks,
                );
                let expected =
                    partial_by_hand(&fixed, &pure_runs, &models, held_right(read, LEFT_FLANK));
                assert!(
                    (scored - expected).abs() <= 1e-12 * expected && expected > 0.0,
                    "({ca}, {ttg}) from the right from {end}: {scored} against {expected}"
                );
                compared += 2;
            }
            for (read, positions) in [(&too_long, (0, 30)), (&too_long_right, (10, 36))] {
                let scored = bundle_censored_scored(
                    read,
                    &witness(positions),
                    &allele,
                    &hull,
                    &runs,
                    &models,
                    flanks,
                );
                assert_eq!(scored, 0.0, "({ca}, {ttg}): longer than any carrier");
            }
        }
        assert_eq!(compared, 18);
    }

    /// **Over three runs a read sums every slip set, and the memo is read back**: at
    /// `G (CA)6 TT (AAG)4 CC (TG)5 A`, the `CA` run three copies longer with the `AAG` run two
    /// shorter reaches the `TG` run at the same offset as no slip at all, so the run-by-run sum
    /// answers that place from what it remembered (str_bundles_partial_reads review, 2026-10-12:
    /// with two runs no remembered value is ever used).
    #[test]
    fn a_partial_read_over_three_runs_sums_every_slip_set() {
        let fixed: [&[u8]; 4] = [b"G", b"TT", b"CC", b"A"];
        let render = |ca: usize, aag: usize, tg: usize| {
            [
                fixed[0],
                &b"CA".repeat(ca),
                fixed[1],
                &b"AAG".repeat(aag),
                fixed[2],
                &b"TG".repeat(tg),
                fixed[3],
            ]
            .concat()
        };
        let hull = render(6, 4, 5);
        let runs = [
            reference_run(1, 12, b"CA"),
            reference_run(15, 12, b"AAG"),
            reference_run(29, 10, b"TG"),
        ];
        let models = [a_model(), another_model(), a_model()];
        let flanks = (LEFT_FLANK, RIGHT_FLANK);
        for (ca, aag, tg) in [(6, 4, 5), (7, 4, 5), (6, 3, 6)] {
            let allele = render(ca, aag, tg);
            let pure_runs = [
                (&b"CA"[..], ca as i64),
                (&b"AAG"[..], aag as i64),
                (&b"TG"[..], tg as i64),
            ];
            let from_left = &hull[..34];
            let from_right = &hull[6..];
            for (read, positions, expected) in [
                (
                    from_left,
                    (0, 34),
                    partial_by_hand(
                        &fixed,
                        &pure_runs,
                        &models,
                        held_left(from_left, RIGHT_FLANK),
                    ),
                ),
                (
                    from_right,
                    (6, 40),
                    partial_by_hand(
                        &fixed,
                        &pure_runs,
                        &models,
                        held_right(from_right, LEFT_FLANK),
                    ),
                ),
            ] {
                let scored = bundle_censored_scored(
                    read,
                    &witness(positions),
                    &allele,
                    &hull,
                    &runs,
                    &models,
                    flanks,
                );
                assert!(
                    (scored - expected).abs() <= 1e-12 * expected && expected > 0.0,
                    "({ca}, {aag}, {tg}) witness {positions:?}: {scored} against {expected}"
                );
            }
        }
    }

    /// **A read whose witness covers the whole hull is the average of both readings**, checked
    /// where both explain it: a pure `(CA)10` hull between flanks that continue the run, so the
    /// read `(CA)11` fits from either end. (The two-run test's right-flank reading of such a read
    /// was 4 parts in 10¹³ of the left one, too small for the tolerance to see.)
    #[test]
    fn an_either_read_averages_both_readings_where_both_explain_it() {
        let hull = b"CA".repeat(10);
        let runs = [reference_run(0, 20, b"CA")];
        let models = [a_model()];
        let flank = b"CACACACACA";
        let read = b"CA".repeat(11);
        let fixed: [&[u8]; 2] = [b"", b""];
        let left = partial_by_hand(
            &fixed,
            &[(&b"CA"[..], 10)],
            &models,
            held_left(&read, flank),
        );
        let right = partial_by_hand(
            &fixed,
            &[(&b"CA"[..], 10)],
            &models,
            held_right(&read, flank),
        );
        let scored = bundle_censored_scored(
            &read,
            &witness((0, 20)),
            &hull,
            &hull,
            &runs,
            &models,
            (flank, flank),
        );
        let expected = 0.5 * (left + right);
        assert!(
            (scored - expected).abs() <= 1e-12 * expected,
            "{scored} against {expected}"
        );
        assert!(
            right > 1e-3 * left && left > 1e-3 * right,
            "both readings explain the read: {left} and {right}"
        );
    }
}
