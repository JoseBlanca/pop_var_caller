//! The repeat-tract half of the joint fit: how often a polymerase slips, per stratum.
//!
//! While a polymerase copies a repeat tract it sometimes **slips**, adding or dropping a whole
//! repeat unit, so a read reports a tract one unit longer or shorter than the DNA it came
//! from. This module estimates how often that happens, from every sample's reads at the same
//! tracts at once. It is not the variant caller: it runs once over the cohort, before anything
//! is called, and hands the caller the numbers it will assume.
//!
//! Design: [`spec/parameter_prepass_joint_fit.md`] §4, §4.1 and §4.2. The
//! [`generic`](super::fit) half runs first and this one depends on it in **one direction
//! only** — it takes each sample's homozygote excess and gives nothing back — so a run may
//! drop the ordinary-position records before reading a single tract.
//!
//! # Vocabulary
//!
//! - **Tract** — one repeat region of the reference; what the records call an STR locus.
//! - **Stratum** — every tract sharing a motif length and a reference repeat count. Slippage
//!   depends on repeat count more than on anything else, so a stratum is the unit that is
//!   fitted. Tomato holds 462,701 kept tracts in 141 strata.
//! - **Offset** — a read's tract length minus the reference tract's, in whole repeat units.
//!   The records store `-4 … +4` with the ends saturating
//!   ([`RECORDED_OFFSET_RANGE`](super::census::RECORDED_OFFSET_RANGE)).
//! - **Length spectrum** — how the stratum's chromosomes are spread over the tract lengths.
//! - **Concentration** — how monomorphic the stratum's tracts are. Small means most tracts are
//!   fixed at one length while the stratum as a whole spans many.
//!
//! # What is fitted
//!
//! Per (read group × stratum), three numbers describing slippage: **how often** a read slips,
//! **which way** — at tomato dinucleotides a slipped read shows a shorter tract 4.9 times as
//! often as a longer one — and **how fast** two-unit slips fall off against one-unit slips.
//! Per stratum, the length spectrum and the concentration. **None of them is per tract.**
//!
//! # Two things the design settles, and this module obeys
//!
//! - **A tract's own length frequencies are a latent vector**, drawn from a fitted per-stratum
//!   Dirichlet and integrated away — never a parameter of the tract. Fitting them per tract
//!   directly moves the slippage level 333-fold depending only on where the search starts
//!   (spec §4.1).
//! - **The integral is a fixed 256-point numerical integration** over that Dirichlet. The
//!   earlier design enumerated the cases where a tract is fixed at one length or segregates
//!   exactly two; it cannot represent a tract carrying three, and the fitted slippage absorbs
//!   the difference — +23.7% where 18% of tracts carry three or more, +722% where nearly all
//!   do. It is withdrawn (spec §4.2).
//!
//! # Where the estimator came from
//!
//! It is lifted from `examples/ng_joint_str_harness.rs`, the program the design was measured
//! with, rather than written again — a second implementation is two things to keep agreeing.
//! The harness's `library` mode fits the same draw both ways and prints the two side by side.
//!
//! Three things are generalised in the lift, and each is the specification catching up with
//! the harness rather than a new idea:
//!
//! 1. **Alleles reach further than the read buckets.** The records store `±4`; the lengths the
//!    fit may place allele mass on reach `±6` (`parameter_prepass_joint_records.md` §3.2),
//!    which is what lets an end bucket be attributed to a distant allele rather than to a far
//!    slip. With the two spans equal the arithmetic is the harness's exactly.
//! 2. **The homozygote excess is per sample**, as it arrives from the ordinary-position half,
//!    where the harness supplied one number for the whole panel.
//! 3. **Slippage is per read group**, as spec §4 has it, where the harness had one set of
//!    slippage numbers. Read groups are named in **slippage groups** so a run may pool them;
//!    one group per read group is the specified default and the widest one.
//!
//! [`spec/parameter_prepass_joint_fit.md`]: ../../../../doc/devel/ng/spec/parameter_prepass_joint_fit.md

use std::collections::BTreeMap;

use rayon::prelude::*;

use crate::ng::parameter_estimation::joint::census::{
    CensusError, CohortCensusEvidence, RECORDED_OFFSET_RANGE, SsrEvidence, SsrLocusState,
};
use crate::ng::parameter_estimation::joint::loci::CensusLoci;
use crate::ng::types::{ContigId, ReadGroupId};

// ---------------------------------------------------------------------
// What a read does
// ---------------------------------------------------------------------

/// The three numbers describing how a polymerase slips, for one read group in one stratum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Slippage {
    /// How often a read reports a tract length other than its allele's.
    pub level: f64,
    /// Of the reads that slip, the share showing a **shorter** tract. Tomato dinucleotides sit
    /// near 0.83 — 2,438 shorter against 501 longer.
    pub shorter_share: f64,
    /// How fast two-repeat slips fall off against one-repeat slips.
    pub fall_off: f64,
}

impl Slippage {
    /// `P(a read reports each bucket | the allele sits at `allele_offset`)`.
    ///
    /// Both the allele and the buckets are counted in whole repeat units from the reference
    /// tract length. The buckets run `-read_span … +read_span`; the allele may sit outside
    /// them, which is the point of the two spans being different.
    ///
    /// **The end buckets get their marginal**, never the probability of sitting exactly on the
    /// edge: the outermost bucket is *at least this many repeats short*, so every step that
    /// would land at or beyond it — including the tail of the geometric fall-off, and the
    /// unslipped read of an allele that is itself outside the recorded range — is summed into
    /// it. Measured at a recorded range of `±1` on a stratum whose alleles reach three repeats
    /// either side, the marginal rule returns the slippage level to within 0.05% where
    /// plugging in the edge costs 33% of it (`parameter_prepass_joint_records.md` §3.2).
    pub fn read_probabilities(&self, allele_offset: i32, read_span: i32) -> Vec<f64> {
        let buckets = (2 * read_span + 1) as usize;
        let mut out = vec![0.0; buckets];
        let slot = |offset: i32| (offset.clamp(-read_span, read_span) + read_span) as usize;

        out[slot(allele_offset)] += 1.0 - self.level;

        // Shorter, then longer. Each direction gives its own weight to every step that lands
        // inside the buckets, and the whole remaining tail to the end bucket it saturates
        // into. `steps` is how many steps are still inside; a step past that cannot be told
        // from any later one.
        for (direction, share) in [(-1_i32, self.shorter_share), (1, 1.0 - self.shorter_share)] {
            let inside_steps = if direction < 0 {
                (allele_offset + read_span).max(0)
            } else {
                (read_span - allele_offset).max(0)
            };
            for step in 1..=inside_steps {
                let weight = (1.0 - self.fall_off) * self.fall_off.powi(step - 1);
                out[slot(allele_offset + direction * step)] += self.level * share * weight;
            }
            // Everything at least `inside_steps + 1` steps away, whose weight telescopes to
            // `fall_off^inside_steps`, lands in the end bucket.
            let tail = self.fall_off.powi(inside_steps);
            out[slot(allele_offset + direction * (inside_steps + 1))] += self.level * share * tail;
        }
        out
    }
}

// ---------------------------------------------------------------------
// The evidence one stratum brings
// ---------------------------------------------------------------------

/// Which stratum a tract is in: its motif length and the reference's repeat count.
///
/// **It lives in the census** ([`census::Stratum`](super::census::Stratum)), because the
/// census is keyed by it — two of a tract record's counts are held per stratum rather than per
/// locus — and re-exported here, where the fit reads it, so a consumer of either module names
/// one type.
pub use super::census::Stratum;
pub use super::share_curve::{
    DEFAULT_FALL_OFF, DEFAULT_SHORTER_SHARE, FittedShare, ShareCurve, ShareCurveConfig,
    ShareCurveSource, ShareShape, ShareSource, blend_share, share_curve_for_a_period,
};
pub use super::slippage_curve::{
    CurveReach, FittedCell, LevelSource, PeriodCurves, SlippageCurve, SlippageCurveConfig,
    blend_level, choose_rise_shape,
};

/// One sample's spanning reads at one tract, split by the slippage group that produced them.
///
/// **Only the samples that put a read on the tract are here.** A sample with no read
/// contributes a likelihood of exactly one whatever the parameters are, so leaving it out is
/// not an approximation — and at three reads a site it is most of the panel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleTractReads {
    pub sample: u32,
    /// Per slippage group that put a read here, its reads in each bucket.
    pub by_group: Vec<(u32, Vec<u32>)>,
}

/// Every sample's reads at one tract.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TractReads {
    pub samples: Vec<SampleTractReads>,
}

/// One stratum's tracts, and what the fit needs to know about the shape of them.
#[derive(Debug, Clone, PartialEq)]
pub struct StratumEvidence {
    pub stratum: Stratum,
    pub tracts: Vec<TractReads>,
    /// How many buckets a read's offset was recorded in, `2 × span + 1`.
    pub read_span: i32,
    /// How many slippage groups the cohort declares. Every group is fitted whether or not it
    /// put a read in this stratum; one that did not is returned as not fitted.
    pub groups: usize,
    /// Tracts left out because more than one read in ten of those differing from the reference
    /// length differed by a **non-whole** number of repeat units — the guard's threshold. Such
    /// a tract is not something this noise model describes, and the fit says so rather than
    /// fitting it (`parameter_prepass_joint_records.md` §3.3).
    pub tracts_over_guard_threshold: u64,
    /// Reads that reached a tract and crossed no whole copy of it, so reported no length.
    /// **Not evidence about slippage and not dropped silently**: a tract longer than a read is
    /// never crossed in any sample at any depth, so this count runs along the repeat-count
    /// axis and a stratum unreadable at this read length must not look like one that was
    /// merely unlucky with coverage.
    pub reads_reaching_not_crossing: u64,
    /// Reads whose tract differed from the reference by a non-whole number of repeat units, at
    /// tracts that stayed in. A diagnostic; nothing about slippage is estimated from them.
    pub guard_reads: u64,
    /// Bases of tract sequence a read was compared against, over every tract and sample — the
    /// denominator of this stratum's substitution rate.
    pub bases_compared: u64,
    /// Of those bases, how many the read disagreed with — the numerator of this stratum's
    /// substitution rate. **Counted only on reads whose tract was the reference's length**,
    /// which is where a mismatch can be read off at all (`census::TractDifference`).
    pub mismatching_bases: u64,
}

impl StratumEvidence {
    /// Tracts carrying at least one spanning read — the count the per-stratum floor is
    /// measured in.
    pub fn tracts_with_reads(&self) -> usize {
        self.tracts.iter().filter(|t| !t.samples.is_empty()).count()
    }

    /// Spanning reads, over every tract and sample.
    pub fn spanning_reads(&self) -> u64 {
        self.tracts
            .iter()
            .flat_map(|tract| tract.samples.iter())
            .flat_map(|sample| sample.by_group.iter())
            .flat_map(|(_, counts)| counts.iter())
            .map(|reads| u64::from(*reads))
            .sum()
    }

    /// Spanning reads whose tract was **not** the reference's length, over every tract and
    /// sample.
    ///
    /// **This is not the count of slipped reads and must not be read as one.** A read sits off
    /// the reference length because the polymerase slipped *or* because the chromosome it came
    /// from genuinely carries another length, and this count cannot tell the two apart — at a
    /// polymorphic tract most of it is the second. It is here because it is the observable that
    /// bounds how sharply a stratum can determine its slippage numbers: a stratum with none of
    /// these determines nothing.
    pub fn reads_off_reference_length(&self) -> u64 {
        self.tracts
            .iter()
            .flat_map(|tract| tract.samples.iter())
            .flat_map(|sample| sample.by_group.iter())
            .flat_map(|(_, counts)| {
                counts
                    .iter()
                    .enumerate()
                    .filter(|(bucket, _)| *bucket as i32 != self.read_span)
            })
            .map(|(_, reads)| u64::from(*reads))
            .sum()
    }

    /// Mismatching bases over bases compared — the stratum's substitution rate, which is one
    /// division and needs none of the other numbers (spec §4.2).
    ///
    /// `None` where no read was compared against a tract at all.
    pub fn substitution_rate(&self) -> Option<f64> {
        (self.bases_compared > 0)
            .then(|| self.mismatching_bases as f64 / self.bases_compared as f64)
    }

    /// Which slippage groups put a read in this stratum.
    fn groups_with_reads(&self) -> Vec<bool> {
        let mut seen = vec![false; self.groups];
        for tract in &self.tracts {
            for sample in &tract.samples {
                for (group, counts) in &sample.by_group {
                    if counts.iter().any(|reads| *reads > 0) {
                        seen[*group as usize] = true;
                    }
                }
            }
        }
        seen
    }
}

// ---------------------------------------------------------------------
// What comes back
// ---------------------------------------------------------------------

/// Everything one stratum's fit produces.
#[derive(Debug, Clone, PartialEq)]
pub struct StratumFit {
    pub stratum: Stratum,
    /// Per slippage group, its three slippage numbers — `None` where that group put no read in
    /// this stratum. **An absent group is not a fitted zero**: a group with no reads here has
    /// no slippage estimate, and saying so is the difference between missing and quiet.
    pub slippage: Vec<Option<Slippage>>,
    /// How the stratum's chromosomes are spread over the allele lengths, indexed from
    /// `-allele_span` to `+allele_span` in whole repeat units.
    pub length_spectrum: Vec<f64>,
    /// How monomorphic the stratum's tracts are. Small means most tracts carry one length.
    pub concentration: f64,
    /// The mean log-likelihood a tract, at the returned parameters.
    pub log_likelihood_a_tract: f64,
    /// Tracts the fit actually read — its own if it stood alone, its own plus its neighbours'
    /// if it borrowed.
    pub tracts_fitted: usize,
    /// The neighbouring repeat counts this stratum borrowed tracts from, empty when it stood
    /// on its own. Same period throughout: slippage is not comparable across motif lengths.
    pub borrowed: Vec<u64>,
    /// Whether the climb settled or ran out of rounds. **Running out is never reported as
    /// convergence.**
    pub converged: bool,
    /// Tracts **this stratum itself** holds with at least one spanning read, whatever it read to
    /// produce its answer.
    ///
    /// **Distinct from [`StratumFit::tracts_fitted`], and the difference is the whole point.** A
    /// stratum with eight tracts of its own that borrowed its way to a thousand has an answer
    /// resting on its neighbours, and a consumer told only the second number cannot see that
    /// (`str_slippage_level_curve.md` §8).
    pub tracts_of_its_own: usize,
    /// Reads that crossed a whole tract of **this stratum itself**, over every sample and group.
    pub reads_crossing: u64,
    /// Per slippage group, where that group's emitted slippage *level* came from — `None` where
    /// the group put no read in this stratum, matching [`StratumFit::slippage`] index for index.
    ///
    /// **The level is the only one of the four numbers a curve supplies.** The direction split
    /// and the fall-off are still the cell's own or its neighbours', so this says nothing about
    /// them.
    pub level_provenance: Vec<Option<LevelProvenance>>,
    /// Per slippage group, where that group's direction split and fall-off came from.
    pub shares_provenance: Vec<Option<SharesProvenance>>,
}

/// Where one slippage group's level at one stratum came from, and what stood behind it.
///
/// **This replaces what `Provenance::Borrowed` used to mark.** After the level becomes a curve, a
/// value fitted from 8,000 slipped reads and one interpolated across a gap look identical in the
/// number alone, and the mechanism that used to distinguish them no longer sets the level
/// (`str_slippage_level_curve.md` §8).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LevelProvenance {
    /// The cell's own fit, the curve, or a blend — and for a blend, the share the curve carried.
    pub source: LevelSource,
    /// The curve that supplied it, absent when this stratum's period had none. It carries the
    /// curve's own held-out error and how many cells stood behind it, so a consumer can tell a
    /// curve through twenty-three cells from one through four.
    pub curve: Option<SlippageCurve>,
    /// Whether this stratum's repeat count sat inside the curve's fitted range. `None` where
    /// there is no curve. **A level held at a fitted end is under-stated in a known direction**
    /// (`str_slippage_level_curve.md` §6).
    pub reach: Option<CurveReach>,
    /// How many of this stratum's own reads **its own fitted level** said slipped, and `None`
    /// where the stratum has no level of its own because it borrowed.
    ///
    /// **The stratum's own level, not the emitted one**, because this is the evidence that stood
    /// behind the cell — it is what set how precisely the stratum could determine its own answer,
    /// and it is the weight the blend gave that answer. Computed from the emitted level it would
    /// be partly a property of the curve, which is the thing it exists to be weighed against.
    ///
    /// **Absent is not zero.** A stratum that borrowed has reads of its own — they are in
    /// [`StratumFit::reads_crossing`] — but no level of its own to say how many of them slipped.
    pub slipped_reads: Option<f64>,
}

/// Where a stratum's direction split and fall-off came from, and what stood behind them.
///
/// **The two shares are smoothed exactly as the level is** — each gets its period's curve, and
/// each stratum departs from that curve by how much evidence it has
/// (`str_slippage_level_curve.md` §5.1). What differs between the three numbers is the shape
/// their curve may take and how their own precision is computed, and nothing else.
///
/// **This replaced a gate with a cliff.** A stratum used to keep its own two shares above 4,000
/// slipped reads and take one named neighbour's whole below it. Across both cohorts one motif
/// period of twelve ever cleared that floor, so 69 of HG002's strata and every one of tomato's
/// got nothing at all.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SharesProvenance {
    /// How many of this stratum's own reads **its own fitted level** said slipped, and `None`
    /// where nothing was fitted here.
    ///
    /// **Both shares are proportions over the reads that slipped**, so this one count sets how
    /// precisely the stratum holds either of them, and it is what the blend weighed its own
    /// answer by. It is the stratum's own level rather than the emitted one: computed from the
    /// emitted level it would be partly a property of the curve, which is the thing it exists to
    /// be weighed against.
    pub slipped_reads: Option<f64>,
    /// Where the share of slipped reads showing a *shorter* tract came from.
    pub shorter_share: ShareProvenance,
    /// Where the fall-off — how much rarer a two-unit slip is than a one-unit slip — came from.
    pub fall_off: ShareProvenance,
}

/// Where one of the two shares came from: the stratum's own fit, its period's curve, or a blend.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShareProvenance {
    /// Which of the three, and for a blend the share the curve carried.
    pub source: ShareSource,
    /// The curve that supplied it, absent where this stratum's period had none. It carries its
    /// own held-out error, how many strata stood behind it, and which rung of the fallback
    /// ladder produced it.
    pub curve: Option<ShareCurve>,
    /// Whether this stratum's repeat count sat inside the curve's fitted range; `None` where
    /// there is no curve. **A share held at a fitted end is the end stratum's answer**, not this
    /// stratum's.
    pub reach: Option<CurveReach>,
}

impl SharesProvenance {
    /// The provenance a stratum's own fit starts with, before any curve is drawn.
    fn own(slipped_reads: f64) -> Self {
        let own = ShareProvenance {
            source: ShareSource::Stratum,
            curve: None,
            reach: None,
        };
        Self {
            slipped_reads: Some(slipped_reads),
            shorter_share: own,
            fall_off: own,
        }
    }
}

/// A stratum whose numbers were all supplied from elsewhere — nothing about it was fitted.
///
/// **All three of its slippage numbers are its period's curves**, which is a complete parameter
/// set for the read likelihood and is what lets a stratum below the refusal floor be emitted at
/// all (`str_slippage_level_curve.md` §1.1, §5.1).
///
/// **It carries no length spectrum, no concentration and no log-likelihood**, because there was
/// no fit to produce them. A separate shape rather than a [`StratumFit`] with those fields left
/// empty, so that a consumer cannot read a spectrum that was never estimated.
#[derive(Debug, Clone, PartialEq)]
pub struct DerivedStratum {
    pub stratum: Stratum,
    /// Per slippage group, its three slippage numbers — `None` where that group put no read here.
    pub slippage: Vec<Option<Slippage>>,
    /// Where each group's level came from; always the curve, since there was no fit to blend.
    pub level_provenance: Vec<Option<LevelProvenance>>,
    /// Where each group's two shares came from; always their period's curve, since there was no
    /// fit here to blend with it.
    pub shares_provenance: Vec<Option<SharesProvenance>>,
    /// Tracts this stratum holds with at least one spanning read — the evidence that was too
    /// thin to fit, and which a consumer still has to be able to see.
    pub tracts_of_its_own: usize,
    /// Reads that crossed a whole tract of it.
    pub reads_crossing: u64,
}

/// Why a stratum produced no fit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StratumRefusal {
    /// Not one tract in the stratum carried a spanning read.
    NoSpanningReads,
    /// The stratum, with everything it could borrow, still holds fewer tracts than the floor
    /// says can carry an answer.
    BelowTheFloor { tracts: usize, floor: usize },
}

/// One stratum's answer: what was fitted, or why nothing was.
#[derive(Debug, Clone, PartialEq)]
pub enum StratumOutcome {
    /// Something was fitted from this stratum's own tracts.
    Fitted(Box<StratumFit>),
    /// Nothing was fitted here; every number came from the curve and from a neighbour.
    Derived(Box<DerivedStratum>),
    /// There is no answer, and saying so is the point.
    Refused {
        stratum: Stratum,
        tracts: usize,
        reason: StratumRefusal,
    },
}

impl StratumOutcome {
    /// Which stratum this is about, whichever of the three it is.
    pub fn stratum(&self) -> Stratum {
        match self {
            Self::Fitted(fit) => fit.stratum,
            Self::Derived(derived) => derived.stratum,
            Self::Refused { stratum, .. } => *stratum,
        }
    }

    /// The slippage numbers a consumer would use, per slippage group — empty when there are
    /// none.
    ///
    /// **A fitted stratum and a derived one are the same thing to the read likelihood**, which
    /// looks up three numbers per candidate and does not ask where they came from
    /// (`read_likelihoods.md` §4.4). What differs is the provenance beside them, and that is why
    /// the two are separate variants rather than one with empty fields.
    pub fn slippage(&self) -> &[Option<Slippage>] {
        match self {
            Self::Fitted(fit) => &fit.slippage,
            Self::Derived(derived) => &derived.slippage,
            Self::Refused { .. } => &[],
        }
    }

    /// Where each group's level came from — empty for a refusal.
    pub fn level_provenance(&self) -> &[Option<LevelProvenance>] {
        match self {
            Self::Fitted(fit) => &fit.level_provenance,
            Self::Derived(derived) => &derived.level_provenance,
            Self::Refused { .. } => &[],
        }
    }

    /// Where each group's two shares came from — empty for a refusal.
    pub fn shares_provenance(&self) -> &[Option<SharesProvenance>] {
        match self {
            Self::Fitted(fit) => &fit.shares_provenance,
            Self::Derived(derived) => &derived.shares_provenance,
            Self::Refused { .. } => &[],
        }
    }

    /// Tracts this stratum holds with at least one spanning read — its own, never a pooled set's.
    pub fn tracts_of_its_own(&self) -> usize {
        match self {
            Self::Fitted(fit) => fit.tracts_of_its_own,
            Self::Derived(derived) => derived.tracts_of_its_own,
            Self::Refused { tracts, .. } => *tracts,
        }
    }
}

// ---------------------------------------------------------------------
// What the run was asked for
// ---------------------------------------------------------------------

/// How many points the integral over a tract's length frequencies is taken on.
///
/// **Fixed, whatever the number of length classes.** A tensor grid needs
/// `nodes^(classes − 1)` points — 576 at three classes and 1.1 × 10¹¹ at thirteen — where this
/// needs the same 256 at every class count. Measured against the grid at three classes it
/// returns the same answers to within 0.3 percentage points on the concentration and 0.2 on
/// every slippage number, in half the time (spec §4.2).
pub const QUADRATURE_POINTS: usize = 256;

/// How far either side of the reference length the fit may place allele mass.
///
/// **Wider than the recorded offsets' `±4`**, and that is what lets an end bucket be
/// attributed to a distant allele rather than to a far slip
/// (`parameter_prepass_joint_records.md` §3.2).
pub const ALLELE_SPAN: i32 = 6;

/// Where one run of the climb starts.
///
/// **The starting points are part of the estimator, not a tuning detail.** They are spread
/// over the slippage level and over how monomorphic tracts are assumed to be, which are the
/// two the objective is least likely to have a single peak in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StartingPoint {
    pub slippage_level: f64,
    pub concentration: f64,
}

impl StartingPoint {
    /// The three the harness climbs from, and the ones every measurement on this path was made
    /// with.
    pub fn spanning_the_monomorphic_range() -> Vec<Self> {
        vec![
            Self {
                slippage_level: 0.02,
                concentration: 0.3,
            },
            Self {
                slippage_level: 0.10,
                concentration: 3.0,
            },
            Self {
                slippage_level: 0.30,
                concentration: 30.0,
            },
        ]
    }
}

#[derive(Debug, Clone)]
pub struct SsrFitConfig {
    /// How far either side of the reference length allele mass may sit.
    pub allele_span: i32,
    pub quadrature_points: usize,
    pub starting_points: Vec<StartingPoint>,
    /// How many rounds of coordinate ascent one start gets.
    pub max_rounds: u32,
    /// A round stops the climb when it improved the mean log-likelihood a tract by less than
    /// this.
    pub stillness: f64,
    /// How the two **shares** are smoothed across repeat count once every stratum has its own
    /// answer; see [`share_curve`](super::share_curve).
    pub share_curve: ShareCurveConfig,
    /// How many tracts a stratum must reach, its borrowings included, before anything is
    /// fitted at all.
    pub refusal_floor: usize,
    /// How the slippage *level* is smoothed across repeat count once every stratum has its own
    /// answer; see [`slippage_curve`](super::slippage_curve).
    ///
    /// **`draw_curves: false` is the arm the parity oracle runs.** Nothing about how a stratum is
    /// fitted depends on this, so a stratum's own level moving between the two arms is a defect
    /// in the plumbing rather than a consequence of the design.
    pub curve: SlippageCurveConfig,
}

impl Default for SsrFitConfig {
    fn default() -> Self {
        Self {
            allele_span: ALLELE_SPAN,
            quadrature_points: QUADRATURE_POINTS,
            starting_points: StartingPoint::spanning_the_monomorphic_range(),
            max_rounds: 5,
            stillness: 1e-6,
            share_curve: ShareCurveConfig::default(),
            refusal_floor: DEFAULT_REFUSAL_FLOOR,
            curve: SlippageCurveConfig::default(),
        }
    }
}

/// How few tracts leave a stratum with nothing worth fitting at all.
///
/// **8, lowered from 50 on 2026-08-20 because a stratum's answer no longer stands alone.** Under
/// the curves every stratum's three numbers are its own answer blended with its motif period's
/// curve, weighted by how precisely it holds them, so a thin stratum's noisy answer costs a
/// consumer nothing — and refusing to fit it costs the period a contributing stratum it could
/// have drawn its curve through.
///
/// **8 is where a thin fit stops breaking and starts merely being noisy**, measured on drawn
/// strata at both ends of the range this caller works over
/// (`examples/ng_ssr_thin_stratum_gate.rs`, 30 draws a row, a truth of 2 reads slipping in 100):
///
/// | tracts | the level came back below a tenth of the truth | median level | 1 in 10 came back above |
/// |---:|---|---:|---:|
/// | 3 | **27%** of one deep sample's fits, 3% of a 63-sample cohort's | 0.0167, 0.0204 | 2.17×, 2.03× |
/// | 5 | **20%**, 0% | 0.0229, 0.0212 | 2.34×, 1.81× |
/// | **8** | **3%**, 0% | 0.0208, 0.0172 | 1.67×, 1.24× |
/// | 12 | 0%, 0% | 0.0204, 0.0196 | 1.59×, 1.41× |
/// | 50 | 0%, 0% | 0.0207, 0.0203 | 1.15×, 1.25× |
///
/// **A collapsed fit excludes itself, which is why the failure is survivable rather than
/// dangerous.** A level that comes back at zero puts zero slipped reads behind the stratum, so it
/// carries no weight in either share's curve and is dropped outright from the level's. What a
/// lower floor risks is the other tail — a level fitted high is weighted high — and at 8 tracts
/// one fit in ten comes back 1.67 times the truth against 1.15 at 50.
///
/// **What it reaches, on the two cohorts' real tables.** Tomato goes from 15 of its 49 populated
/// strata carrying a full parameter set to **38**, because its dinucleotides reach four
/// contributing strata for the first time. HG002's count does not move — its thin periods are
/// thinner than this floor — but its trinucleotide curve is drawn through 10 strata rather than
/// 4, and 4 is where its held-out error was 32%.
///
/// **Why not lower still.** At 5 tracts one fit in five collapses on a single deep sample, which
/// is the shape HG002 has; and going to 5 or 3 buys 10 more strata on HG002 and 8 more on tomato,
/// from periods whose curves would then rest on strata a fifth of which said nothing.
///
/// *Whether the climb "converged" is not the alarm the specification expected it to be*: at 400
/// tracts on a 63-sample cohort only 83% of fits settle within their rounds, and that row's
/// median level is 1.5% from the truth. Convergence counts rounds, not quality.
pub const DEFAULT_REFUSAL_FLOOR: usize = 8;

// ---------------------------------------------------------------------
// Fitting one stratum
// ---------------------------------------------------------------------

/// What the climb carries between evaluations.
#[derive(Debug, Clone, PartialEq)]
struct Parameters {
    /// Per slippage group.
    slippage: Vec<Slippage>,
    /// Over the allele classes, `-allele_span … +allele_span`.
    length_spectrum: Vec<f64>,
    concentration: f64,
}

impl Parameters {
    fn start(start: StartingPoint, groups: usize, classes: usize) -> Self {
        Self {
            slippage: vec![
                Slippage {
                    level: start.slippage_level,
                    shorter_share: 0.5,
                    fall_off: 0.3,
                };
                groups
            ],
            length_spectrum: vec![1.0 / classes as f64; classes],
            concentration: start.concentration,
        }
    }
}

/// Fit one stratum, on its own tracts alone.
///
/// `homozygote_excess` is one number a sample, indexed the way [`SampleTractReads::sample`] is:
/// how far short of the heterozygote proportions the population's allele frequencies predict
/// that sample falls. It arrives from the ordinary-position half and is **held, never fitted
/// here**.
///
/// # Panics
///
/// When a tract names a sample the homozygote-excess list does not reach — the two are one
/// cohort described twice, and a mismatch is a wiring error rather than data.
pub fn fit_stratum(
    evidence: &StratumEvidence,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Option<StratumFit> {
    fit_pooled(evidence, &[], homozygote_excess, config)
}

/// Fit `evidence`, whose tracts may already include borrowed ones, recording where from.
fn fit_pooled(
    evidence: &StratumEvidence,
    borrowed: &[u64],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Option<StratumFit> {
    let classes = (2 * config.allele_span + 1) as usize;
    let genotypes = genotype_pairs(classes);
    let live_groups = evidence.groups_with_reads();
    if !live_groups.iter().any(|live| *live) {
        return None;
    }

    let mut best: Option<(Parameters, f64, bool)> = None;
    for start in &config.starting_points {
        let mut parameters = Parameters::start(*start, evidence.groups, classes);
        let mut scorer = Scorer::new(evidence, homozygote_excess, &genotypes, config);
        let mut score = scorer.score(&parameters);
        let mut converged = false;
        for _ in 0..config.max_rounds {
            let before = score;
            climb_one_round(&mut parameters, &mut scorer, &live_groups, classes);
            score = scorer.score(&parameters);
            if score - before < config.stillness {
                converged = true;
                break;
            }
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_score, _)| score > *best_score)
        {
            best = Some((parameters, score, converged));
        }
    }

    let (parameters, score, converged) = best.expect("at least one starting point");
    Some(StratumFit {
        stratum: evidence.stratum,
        slippage: live_groups
            .iter()
            .enumerate()
            .map(|(group, live)| live.then_some(parameters.slippage[group]))
            .collect(),
        length_spectrum: parameters.length_spectrum,
        concentration: parameters.concentration,
        log_likelihood_a_tract: score,
        tracts_fitted: evidence.tracts_with_reads(),
        borrowed: borrowed.to_vec(),
        converged,
        // **What the stratum holds on its own is not knowable here.** `evidence` may already be
        // the pooled set, so these two are placeholders that `fit_strata` replaces with the
        // receiving stratum's own counts, exactly as it replaces `stratum` and `borrowed`.
        tracts_of_its_own: evidence.tracts_with_reads(),
        reads_crossing: evidence.spanning_reads(),
        // Every stratum starts owning its own shares, with the slipped-read count they rest on.
        // `fit_strata` re-emits them through their period's curves.
        shares_provenance: live_groups
            .iter()
            .enumerate()
            .map(|(group, live)| {
                live.then(|| {
                    SharesProvenance::own(
                        parameters.slippage[group].level * evidence.spanning_reads() as f64,
                    )
                })
            })
            .collect(),
        // Every level starts as the stratum's own fit. Drawing curves is step B3's; until then
        // this records the truth, which is that no curve touched it.
        level_provenance: live_groups
            .iter()
            .enumerate()
            .map(|(group, live)| {
                live.then(|| LevelProvenance {
                    source: LevelSource::Cell,
                    curve: None,
                    reach: None,
                    slipped_reads: Some(
                        parameters.slippage[group].level * evidence.spanning_reads() as f64,
                    ),
                })
            })
            .collect(),
    })
}

/// One pass of coordinate ascent over everything the stratum fits.
///
/// **The slippage numbers move first and the frequencies after**, because moving slippage is
/// what invalidates the cached read likelihoods; doing it the other way round would rebuild
/// them on every spectrum coordinate.
fn climb_one_round(
    parameters: &mut Parameters,
    scorer: &mut Scorer<'_>,
    live_groups: &[bool],
    classes: usize,
) {
    for (group, live) in live_groups.iter().enumerate() {
        if !live {
            continue;
        }
        for which in 0..3 {
            let current = read_slippage(&parameters.slippage[group], which);
            let moved = climb_scalar(
                |x| {
                    let mut trial = parameters.clone();
                    write_slippage(&mut trial.slippage[group], which, expit(x));
                    scorer.score(&trial)
                },
                logit(current),
                3.0,
            );
            write_slippage(&mut parameters.slippage[group], which, expit(moved));
        }
    }

    // The spectrum, one class at a time on a log scale, renormalised each time.
    for class in 0..classes {
        let current = parameters.length_spectrum[class].max(1e-9).ln();
        let moved = climb_scalar(
            |x| {
                let mut trial = parameters.clone();
                trial.length_spectrum[class] = x.exp();
                normalise(&mut trial.length_spectrum);
                scorer.score(&trial)
            },
            current,
            2.0,
        );
        parameters.length_spectrum[class] = moved.exp();
        normalise(&mut parameters.length_spectrum);
    }

    let moved = climb_scalar(
        |x| {
            let mut trial = parameters.clone();
            trial.concentration = x.exp();
            scorer.score(&trial)
        },
        parameters.concentration.ln(),
        2.5,
    );
    parameters.concentration = moved.exp();
}

fn read_slippage(slippage: &Slippage, which: usize) -> f64 {
    match which {
        0 => slippage.level,
        1 => slippage.shorter_share,
        _ => slippage.fall_off,
    }
}

fn write_slippage(slippage: &mut Slippage, which: usize, value: f64) {
    match which {
        0 => slippage.level = value,
        1 => slippage.shorter_share = value,
        _ => slippage.fall_off = value,
    }
}

fn normalise(weights: &mut [f64]) {
    let total: f64 = weights.iter().sum();
    for weight in weights.iter_mut() {
        *weight /= total;
    }
}

// ---------------------------------------------------------------------
// Thin strata: every slippage number from its period's curve
// ---------------------------------------------------------------------

/// Fit every stratum on its own tracts, then draw a curve a motif period through what they
/// measured and re-emit every number through it.
///
/// **Three things happen, in this order.** Every stratum holding at least
/// [`SsrFitConfig::refusal_floor`] tracts of its own is fitted from them and no others'. Then one
/// curve a motif period is drawn for each of the three slippage numbers, through the strata that
/// stood on their own tracts and weighted by how precisely each holds its own answer. Then every
/// stratum's numbers are re-emitted as a blend of its own answer and its period's curve, and a
/// stratum too thin to have been fitted at all takes the curves whole.
///
/// **All three curves are drawn before either blend runs**, so no curve is ever fitted to a
/// number another curve emitted — the circularity `str_slippage_level_curve.md` §5.1 forbids.
///
/// **Nothing pools tracts and nothing copies a neighbour's shares any more.** Pooling was how a
/// thin stratum got an answer at all, and copying was how it got its two shares; both are
/// replaced by a curve through every stratum of the period, weighted, which reaches strata
/// neither could. Removing pooling also removed the run's expensive arm — 1,036.8 s against
/// 155.5 s on the same cohort (`str_slippage_level_curve.md` §5.1).
pub fn fit_strata(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Vec<StratumOutcome> {
    // **Every stratum is fitted on its own tracts and no other's.** Pooling a thin stratum's
    // tracts with its neighbours' and refitting used to be how it got an answer at all; it now
    // gets its level from its period's curve and its two shares copied from a neighbour that
    // measured them well, so there is nothing left for a pooled fit to supply
    // (`str_slippage_level_curve.md` §5.1). What that removes is the run's expensive arm: the
    // same cohort took 1,036.8 s with pooled borrowing against 155.5 s without.
    let mut outcomes: Vec<StratumOutcome> = strata
        .iter()
        .map(|evidence| {
            if evidence.tracts_with_reads() == 0 {
                return StratumOutcome::Refused {
                    stratum: evidence.stratum,
                    tracts: 0,
                    reason: StratumRefusal::NoSpanningReads,
                };
            }
            if evidence.tracts_with_reads() < config.refusal_floor {
                // Too thin to fit anything of its own. It is not refused yet — the curve and a
                // neighbour's shares may still furnish it, which is what D3 below decides.
                return StratumOutcome::Refused {
                    stratum: evidence.stratum,
                    tracts: evidence.tracts_with_reads(),
                    reason: StratumRefusal::BelowTheFloor {
                        tracts: evidence.tracts_with_reads(),
                        floor: config.refusal_floor,
                    },
                };
            }
            match fit_pooled(evidence, &[], homozygote_excess, config) {
                Some(fit) => StratumOutcome::Fitted(Box::new(fit)),
                None => StratumOutcome::Refused {
                    stratum: evidence.stratum,
                    tracts: evidence.tracts_with_reads(),
                    reason: StratumRefusal::NoSpanningReads,
                },
            }
        })
        .collect();

    // **The curves are drawn after every stratum has its own answer, never during.** Stage one
    // is untouched by this — a stratum's own fitted numbers are the same whether curves are drawn
    // or not, which is the property the parity oracle checks.
    if config.curve.draw_curves {
        // **All three curves are drawn here, before either blend runs**, and both read fits that
        // nothing has touched. Drawing a curve after the levels had been blended would fit a
        // curve to the previous curve's output, which is the circularity
        // `str_slippage_level_curve.md` §5.1 forbids in so many words.
        let levels = draw_a_curve_a_period(&outcomes, config);
        let shares = draw_share_curves_a_period(&outcomes, config);

        smooth_levels_across_repeat_count(&mut outcomes, &levels, config);
        smooth_shares_across_repeat_count(&mut outcomes, &shares, config);
        // **Last, because it needs all three**: a stratum too thin to fit anything of its own is
        // furnished from its period's curves rather than refused.
        derive_thin_strata(&mut outcomes, strata, &levels, &shares);
    }
    outcomes
}

/// Turn a stratum too thin to fit anything of its own into one furnished from elsewhere.
///
/// **This is what spec §1.1's first goal asks for: every populated stratum carries a level.** All
/// three of its numbers are its period's curves. Since nothing about it was estimated it comes
/// back as a [`DerivedStratum`] — no length spectrum, no concentration, no log-likelihood —
/// rather than as a [`StratumFit`] with those left empty.
///
/// **Two refusals survive**: a stratum no read crossed, and one whose period has no *level*
/// curve. The second is the only floor left standing — a period needs
/// [`SlippageCurveConfig::min_cells_for_a_curve`] strata fitted on their own tracts before a
/// level curve is drawn at all, where the two shares always have a curve to give.
fn derive_thin_strata(
    outcomes: &mut [StratumOutcome],
    strata: &[StratumEvidence],
    curves: &BTreeMap<u8, PeriodCurves>,
    share_curves: &BTreeMap<(u8, usize), SharesCurves>,
) {
    let evidence_of: BTreeMap<Stratum, &StratumEvidence> = strata
        .iter()
        .map(|evidence| (evidence.stratum, evidence))
        .collect();

    for outcome in outcomes.iter_mut() {
        let StratumOutcome::Refused {
            stratum, reason, ..
        } = outcome
        else {
            continue;
        };
        // A stratum no read crossed has nothing to furnish, and says so.
        if matches!(reason, StratumRefusal::NoSpanningReads) {
            continue;
        }
        let stratum = *stratum;
        let Some(evidence) = evidence_of.get(&stratum) else {
            continue;
        };
        let with_reads = evidence.groups_with_reads();
        let repeats = stratum.reference_repeats;

        let mut slippage: Vec<Option<Slippage>> = vec![None; with_reads.len()];
        let mut level_provenance: Vec<Option<LevelProvenance>> = vec![None; with_reads.len()];
        let mut shares_provenance: Vec<Option<SharesProvenance>> = vec![None; with_reads.len()];
        let mut furnished_any = false;

        for (group, live) in with_reads.iter().enumerate() {
            if !live {
                continue;
            }
            let curve = curves
                .get(&stratum.period)
                .and_then(|period| period.by_group.get(group))
                .and_then(Option::as_ref);
            let (Some(curve), Some(shares)) = (curve, share_curves.get(&(stratum.period, group)))
            else {
                continue;
            };
            slippage[group] = Some(Slippage {
                level: curve.level_at(repeats),
                shorter_share: shares.shorter_share.share_at(repeats),
                fall_off: shares.fall_off.share_at(repeats),
            });
            level_provenance[group] = Some(LevelProvenance {
                source: LevelSource::Curve,
                curve: Some(*curve),
                reach: Some(curve.reach(repeats)),
                // Nothing was fitted here, so there is no level of its own to count against.
                slipped_reads: None,
            });
            shares_provenance[group] = Some(SharesProvenance {
                slipped_reads: None,
                shorter_share: ShareProvenance {
                    source: ShareSource::Curve,
                    curve: Some(shares.shorter_share),
                    reach: Some(shares.shorter_share.reach(repeats)),
                },
                fall_off: ShareProvenance {
                    source: ShareSource::Curve,
                    curve: Some(shares.fall_off),
                    reach: Some(shares.fall_off.reach(repeats)),
                },
            });
            furnished_any = true;
        }

        if furnished_any {
            *outcome = StratumOutcome::Derived(Box::new(DerivedStratum {
                stratum,
                slippage,
                level_provenance,
                shares_provenance,
                tracts_of_its_own: evidence.tracts_with_reads(),
                reads_crossing: evidence.spanning_reads(),
            }));
        }
    }
}

/// One curve a motif period, drawn through the strata fitted from their own tracts.
///
/// **Only a stratum fitted from its own tracts feeds a curve.** A stratum furnished from
/// elsewhere carries a curve's answer already, and fitting a curve through it would be fitting a
/// curve to its own output — the circularity `str_slippage_level_curve.md` §4 exists to prevent.
fn draw_a_curve_a_period(
    outcomes: &[StratumOutcome],
    config: &SsrFitConfig,
) -> BTreeMap<u8, PeriodCurves> {
    let groups = outcomes
        .iter()
        .map(|outcome| outcome.slippage().len())
        .max()
        .unwrap_or(0);
    if groups == 0 {
        return BTreeMap::new();
    }

    // One list of contributing cells per slippage group, per motif period.
    let mut by_period: BTreeMap<u8, Vec<Vec<FittedCell>>> = BTreeMap::new();
    for outcome in outcomes.iter() {
        let StratumOutcome::Fitted(fit) = outcome else {
            continue;
        };
        if !fit.borrowed.is_empty() {
            continue;
        }
        let cells = by_period
            .entry(fit.stratum.period)
            .or_insert_with(|| vec![Vec::new(); groups]);
        for (group, slippage) in fit.slippage.iter().enumerate() {
            let Some(slippage) = slippage else { continue };
            cells[group].push(FittedCell {
                repeats: fit.stratum.reference_repeats,
                level: slippage.level,
                slipped_reads: slippage.level * fit.reads_crossing as f64,
            });
        }
    }

    by_period
        .iter()
        .filter_map(|(period, cells)| {
            choose_rise_shape(cells, &config.curve)
                .ok()
                .map(|curves| (*period, curves))
        })
        .collect()
}

/// Draw one curve per motif period and re-emit every stratum's level through it.
///
/// **Only the level moves.** The direction split, the fall-off and the length spectrum are the
/// stratum's own and are not touched (`str_slippage_level_curve.md` §1.2).
///
/// **Only a stratum fitted from its own tracts feeds a curve.** A stratum that borrowed carries
/// its neighbours' slippage already, and fitting a curve through it would be fitting a curve to
/// its own output — the circularity `str_slippage_level_curve.md` §4 exists to prevent. So a run
/// that draws curves is meant to fit stage one with borrowing off; with borrowing on, few strata
/// contribute and the curves are correspondingly thin.
fn smooth_levels_across_repeat_count(
    outcomes: &mut [StratumOutcome],
    curves: &BTreeMap<u8, PeriodCurves>,
    config: &SsrFitConfig,
) {
    for outcome in outcomes.iter_mut() {
        let StratumOutcome::Fitted(fit) = outcome else {
            continue;
        };
        let Some(period_curves) = curves.get(&fit.stratum.period) else {
            continue;
        };
        let repeats = fit.stratum.reference_repeats;
        let reads_crossing = fit.reads_crossing as f64;
        // **A stratum that borrowed brings no level of its own to be weighed.** Its pooled level
        // is its neighbours' — the thing the curve replaces — so the curve supplies the level
        // outright and the pooled fit keeps only the two shares. That is the whole of
        // `str_slippage_level_curve.md` §5: the level stops borrowing and nothing else does.
        let borrowed_from_neighbours = !fit.borrowed.is_empty();
        for group in 0..fit.slippage.len() {
            let Some(own) = fit.slippage[group] else {
                continue;
            };
            let curve = period_curves.by_group.get(group).and_then(Option::as_ref);
            let cell = (!borrowed_from_neighbours).then_some(FittedCell {
                repeats,
                level: own.level,
                slipped_reads: own.level * reads_crossing,
            });
            let Some(blended) = blend_level(cell, curve, repeats, &config.curve) else {
                continue;
            };
            if let Some(slippage) = fit.slippage[group].as_mut() {
                slippage.level = blended.level;
            }
            if let Some(provenance) = fit.level_provenance[group].as_mut() {
                provenance.source = blended.source;
                provenance.curve = curve.copied();
                provenance.reach = blended.reach;
                provenance.slipped_reads = cell.map(|cell| cell.slipped_reads);
            }
        }
    }
}

/// One motif period and slippage group's two share curves.
#[derive(Debug, Clone, Copy, PartialEq)]
struct SharesCurves {
    shorter_share: ShareCurve,
    fall_off: ShareCurve,
}

/// One curve for each share, per motif period and slippage group.
///
/// **Only a stratum fitted from its own tracts feeds a curve**, and it feeds it with its own
/// share and its own slipped-read count — never a blended one, or each round of smoothing would
/// fit a curve to the previous round's curve.
///
/// **A curve always comes back for a period that has a populated stratum**, even one where
/// nothing was fitted: `share_curve_for_a_period` falls back to the run's other periods and then
/// to a built-in constant, recording which in the curve's own provenance. These numbers are a
/// prior the read likelihood consults, so answering coarsely beats refusing.
fn draw_share_curves_a_period(
    outcomes: &[StratumOutcome],
    config: &SsrFitConfig,
) -> BTreeMap<(u8, usize), SharesCurves> {
    let groups = outcomes
        .iter()
        .map(|outcome| outcome.slippage().len())
        .max()
        .unwrap_or(0);
    if groups == 0 {
        return BTreeMap::new();
    }

    // Every stratum's own two shares, keyed by motif period and slippage group.
    let mut fitted: BTreeMap<(u8, usize), (Vec<FittedShare>, Vec<FittedShare>)> = BTreeMap::new();
    for outcome in outcomes.iter() {
        let StratumOutcome::Fitted(fit) = outcome else {
            continue;
        };
        // **A stratum that read another's tracts carries another's shares**, so it consumes a
        // curve and does not feed one. Nothing pools tracts today; the guard keeps the rule true
        // by construction rather than by the absence of pooling.
        if !fit.borrowed.is_empty() {
            continue;
        }
        for (group, slippage) in fit.slippage.iter().enumerate() {
            let Some(slippage) = slippage else { continue };
            let slipped_reads = slippage.level * fit.reads_crossing as f64;
            let repeats = fit.stratum.reference_repeats;
            let here = fitted.entry((fit.stratum.period, group)).or_default();
            here.0.push(FittedShare {
                repeats,
                share: slippage.shorter_share,
                slipped_reads,
            });
            here.1.push(FittedShare {
                repeats,
                share: slippage.fall_off,
                slipped_reads,
            });
        }
    }

    // Every period a populated stratum sits at, whether or not anything there was fitted.
    let mut wanted: Vec<(u8, usize)> = Vec::new();
    for outcome in outcomes.iter() {
        if matches!(
            outcome,
            StratumOutcome::Refused {
                reason: StratumRefusal::NoSpanningReads,
                ..
            }
        ) {
            continue;
        }
        for group in 0..groups {
            wanted.push((outcome.stratum().period, group));
        }
    }
    wanted.sort_unstable();
    wanted.dedup();

    wanted
        .into_iter()
        .map(|(period, group)| {
            let empty = (Vec::new(), Vec::new());
            let here = fitted.get(&(period, group)).unwrap_or(&empty);
            // The same slippage group at every *other* motif period — the rung below this
            // period's own strata, and the one the curve records as crossing periods.
            let mut elsewhere: (Vec<FittedShare>, Vec<FittedShare>) = (Vec::new(), Vec::new());
            for ((other_period, other_group), shares) in &fitted {
                if *other_group == group && *other_period != period {
                    elsewhere.0.extend_from_slice(&shares.0);
                    elsewhere.1.extend_from_slice(&shares.1);
                }
            }
            (
                (period, group),
                SharesCurves {
                    shorter_share: share_curve_for_a_period(
                        &here.0,
                        &elsewhere.0,
                        DEFAULT_SHORTER_SHARE,
                        &config.share_curve,
                    ),
                    fall_off: share_curve_for_a_period(
                        &here.1,
                        &elsewhere.1,
                        DEFAULT_FALL_OFF,
                        &config.share_curve,
                    ),
                },
            )
        })
        .collect()
}

/// Re-emit every fitted stratum's two shares through its period's curves.
///
/// **Only the two shares move.** The level has already had its own curve, and the length spectrum
/// and concentration are the stratum's own and are not touched.
///
/// **The weight each stratum's own answer carries is its own slipped-read count**, read from the
/// provenance rather than recomputed from the emitted level — by this point the level has been
/// blended, and a weight computed from it would be partly a property of the level's curve.
fn smooth_shares_across_repeat_count(
    outcomes: &mut [StratumOutcome],
    curves: &BTreeMap<(u8, usize), SharesCurves>,
    config: &SsrFitConfig,
) {
    for outcome in outcomes.iter_mut() {
        let StratumOutcome::Fitted(fit) = outcome else {
            continue;
        };
        let repeats = fit.stratum.reference_repeats;
        let period = fit.stratum.period;
        for group in 0..fit.slippage.len() {
            let Some(own) = fit.slippage[group] else {
                continue;
            };
            let Some(shares) = curves.get(&(period, group)) else {
                continue;
            };
            // **A stratum that read another's tracts brings no shares of its own to weigh**, so
            // the curve supplies them outright — the same rule the level follows.
            let own_slipped_reads = (fit.borrowed.is_empty())
                .then(|| {
                    fit.shares_provenance[group].and_then(|provenance| provenance.slipped_reads)
                })
                .flatten();

            let blend = |share: f64, curve: &ShareCurve| {
                blend_share(
                    own_slipped_reads.map(|slipped_reads| FittedShare {
                        repeats,
                        share,
                        slipped_reads,
                    }),
                    Some(curve),
                    repeats,
                    &config.share_curve,
                )
            };
            let (Some(shorter), Some(fall_off)) = (
                blend(own.shorter_share, &shares.shorter_share),
                blend(own.fall_off, &shares.fall_off),
            ) else {
                continue;
            };

            if let Some(slippage) = fit.slippage[group].as_mut() {
                slippage.shorter_share = shorter.share;
                slippage.fall_off = fall_off.share;
            }
            if let Some(provenance) = fit.shares_provenance[group].as_mut() {
                if own_slipped_reads.is_none() {
                    provenance.slipped_reads = None;
                }
                provenance.shorter_share = ShareProvenance {
                    source: shorter.source,
                    curve: Some(shares.shorter_share),
                    reach: shorter.reach,
                };
                provenance.fall_off = ShareProvenance {
                    source: fall_off.source,
                    curve: Some(shares.fall_off),
                    reach: fall_off.reach,
                };
            }
        }
    }
}

// ---------------------------------------------------------------------
// The likelihood
// ---------------------------------------------------------------------

/// The genotypes of a diploid over `classes` allele lengths, unordered pairs.
fn genotype_pairs(classes: usize) -> Vec<(usize, usize)> {
    (0..classes)
        .flat_map(|first| (first..classes).map(move |second| (first, second)))
        .collect()
}

/// One tract's read likelihoods, in the form the integral sweeps over cheaply.
///
/// **Rescaled out of log space, once.** The inner loop is *for each point of the integral, for
/// each sample, sum over genotypes*, so doing it in logs would cost one `exp` per genotype per
/// sample **per point** — which at 256 points is where the whole run's time would go.
/// Subtracting each sample's largest log-likelihood makes the inner sum a plain dot product,
/// and the offsets are added back at the end.
struct TractLikelihoods {
    /// One row a sample-with-reads, over the genotypes, **laid end to end in one buffer**: row
    /// `r` is `scaled[r * width..][..width]`.
    ///
    /// **Flat and not a vector of vectors**, because the innermost loop sweeps a row against the
    /// genotype prior in vector lanes: separate heap blocks a sample cost a dependent load before
    /// each row and leave the rows scattered, where one buffer keeps them contiguous and lets the
    /// loop walk the whole tract without leaving cache.
    scaled: Vec<f64>,
    /// How wide one row is — the number of genotype pairs.
    width: usize,
    /// That sample's homozygote excess, in the same order.
    homozygote_excess: Vec<f64>,
    /// Σ over those samples of the log-likelihood each row was divided by.
    ln_offset: f64,
}

impl TractLikelihoods {
    fn of(
        tract: &TractReads,
        per_group_allele: &[Vec<Vec<f64>>],
        genotypes: &[(usize, usize)],
        homozygote_excess: &[f64],
    ) -> Self {
        let width = genotypes.len();
        let mut scaled: Vec<f64> = Vec::with_capacity(tract.samples.len() * width);
        let mut excess = Vec::with_capacity(tract.samples.len());
        let mut ln_offset = 0.0;
        for sample in &tract.samples {
            let start = scaled.len();
            scaled.extend(genotypes.iter().map(|(first, second)| {
                let mut total = 0.0;
                for (group, counts) in &sample.by_group {
                    let per_allele = &per_group_allele[*group as usize];
                    for (bucket, reads) in counts.iter().enumerate() {
                        if *reads == 0 {
                            continue;
                        }
                        let probability =
                            0.5 * (per_allele[*first][bucket] + per_allele[*second][bucket]);
                        total += f64::from(*reads) * probability.max(1e-300).ln();
                    }
                }
                total
            }));
            let row = &mut scaled[start..];
            let largest = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            ln_offset += largest;
            for value in row {
                *value = (*value - largest).exp();
            }
            excess.push(
                *homozygote_excess
                    .get(sample.sample as usize)
                    .expect("every sample in a tract has a homozygote excess"),
            );
        }
        Self {
            scaled,
            width,
            homozygote_excess: excess,
            ln_offset,
        }
    }
}

/// One point of the simplex a tract's length frequencies are integrated over, with the weight
/// it stands for.
///
/// **The genotype prior is carried here rather than rebuilt per tract**, because it is a function
/// of the point alone: a tract enters it only through the read likelihoods it is multiplied by.
/// Building it per tract cost `points × genotypes` fills for every tract of every objective
/// evaluation, where `points × genotypes` once a quadrature is the whole of it.
struct Quadrature {
    /// `point × classes` — the allele-length frequencies at each point, **laid out flat**, in
    /// the same shape [`independent`](Self::independent) already uses.
    ///
    /// **One buffer and not one `Vec` a point.** A vector of vectors cost one heap block per
    /// point per rebuild, and the climb rebuilds the quadrature whenever it moves the spectrum
    /// or the concentration: measured with dhat over one profiling run, that was 1,942,272 of
    /// the repeat-tract fit's 4,078,263 blocks — 48% of every allocation the half makes — for
    /// 7,587 rebuilds of 256 points each. Flat, a rebuild is one block.
    frequencies: Vec<f64>,
    /// How many allele classes one point holds — the stride of `frequencies`.
    classes: usize,
    ln_weight: f64,
    /// `point × genotypes` — the chance of drawing this genotype from two independent draws at
    /// this point, laid out flat.
    independent: Vec<f64>,
    /// The genotype slot of each homozygous pair, in allele order.
    ///
    /// **The by-descent half of the prior is zero everywhere else**, so it is a sum over the
    /// thirteen allele classes wearing a ninety-one-slot loop until this list splits it out. The
    /// slot's own value at a point is `frequencies[point][class]`, so nothing else need be stored.
    diagonal: Vec<usize>,
}

/// `ln P(this tract's panel | parameters)`.
fn ln_tract(
    likelihoods: &TractLikelihoods,
    quadrature: &Quadrature,
    genotypes: &[(usize, usize)],
) -> f64 {
    let width = genotypes.len();
    let mut terms = Vec::with_capacity(quadrature.frequencies.len() / quadrature.classes);
    // The genotype prior splits into the part that comes from the two copies being identical
    // by descent and the part that comes from two independent draws, so a sample's own
    // homozygote excess weights two dot products rather than rebuilding the prior per sample.
    // **Both halves are built once with the quadrature**, not once a tract.
    for (index, point) in quadrature
        .frequencies
        .chunks_exact(quadrature.classes)
        .enumerate()
    {
        let independent = &quadrature.independent[index * width..][..width];
        // **One logarithm a point, not one a sample.** The samples' likelihoods multiply, so the
        // product is carried directly and rescaled when it is about to underflow — the same trick
        // the ordinary-position half uses (`fit.rs`'s `RESCALE`), and worth the cohort size: eight
        // logarithms a point become one here, and three thousand become one at the top of the
        // range this caller is for.
        let mut product = 1.0_f64;
        let mut scaled_by = 0.0_f64;
        let mut vanished = false;
        for (row, excess) in likelihoods
            .scaled
            .chunks_exact(likelihoods.width)
            .zip(&likelihoods.homozygote_excess)
        {
            // **Four running sums rather than one, and that is the whole trick.** A single
            // accumulator makes ninety-one additions that each wait for the one before, so the
            // loop runs at the latency of an addition and the machine's vector lanes sit idle —
            // and the compiler may not split it itself, because reassociating floating-point
            // addition changes the answer and Rust does not allow it uninvited. Splitting it here
            // says which association we want, in the source, where it is reproducible.
            let mut lanes = wide::f64x4::ZERO;
            let mut weights = independent.chunks_exact(4);
            let mut values = row[..width].chunks_exact(4);
            for (weight, value) in weights.by_ref().zip(values.by_ref()) {
                let weight = wide::f64x4::new([weight[0], weight[1], weight[2], weight[3]]);
                let value = wide::f64x4::new([value[0], value[1], value[2], value[3]]);
                lanes += weight * value;
            }
            let parts = lanes.to_array();
            let mut at_random = (parts[0] + parts[1]) + (parts[2] + parts[3]);
            for (weight, value) in weights.remainder().iter().zip(values.remainder()) {
                at_random += weight * value;
            }
            // **Only the homozygous slots**: the by-descent prior is zero at every heterozygous
            // pair, so seventy-eight of ninety-one products were a multiply by zero.
            let mut by_descent = 0.0;
            for (class, slot) in quadrature.diagonal.iter().enumerate() {
                by_descent += point[class] * row[*slot];
            }
            let sum = excess * by_descent + (1.0 - excess) * at_random;
            if sum <= 0.0 {
                vanished = true;
                break;
            }
            product *= sum;
            if product < 1.0 / RESCALE {
                product *= RESCALE;
                scaled_by -= LN_RESCALE;
            }
        }
        terms.push(if vanished {
            f64::NEG_INFINITY
        } else {
            quadrature.ln_weight + product.ln() + scaled_by
        });
    }
    ln_sum_exp(&terms) + likelihoods.ln_offset
}

/// The scale the running product above is multiplied back up by when it is about to underflow,
/// and its logarithm. Both are `fit.rs`'s values, and deliberately so: the two halves of the fit
/// carry the same rescale and should one day share it.
const RESCALE: f64 = 1e150;
const LN_RESCALE: f64 = 345.398_899_014_487; // ln(1e150)

/// The evidence with the read likelihoods already computed for one set of slippage numbers.
///
/// **Held between evaluations**, because the climb moves the spectrum and the concentration
/// far more often than it moves slippage, and neither changes a read likelihood. Recomputing
/// them inside the search is what once made this program too slow to run.
struct Prepared {
    slippage: Vec<Slippage>,
    tracts: Vec<TractLikelihoods>,
}

/// The evidence plus whichever slippage the last question was about.
struct Scorer<'a> {
    evidence: &'a StratumEvidence,
    homozygote_excess: &'a [f64],
    genotypes: &'a [(usize, usize)],
    allele_span: i32,
    quadrature_points: usize,
    prepared: Option<Prepared>,
    /// The last integral built, with the spectrum and concentration it was built for.
    ///
    /// **The slippage climb asks about dozens of parameter sets that leave the tract's length
    /// frequencies alone**, and building the integral costs 19 ms at thirteen allele classes —
    /// which on a thin stratum is more than the likelihood it feeds.
    held_quadrature: Option<(Vec<f64>, f64, Quadrature)>,
}

impl<'a> Scorer<'a> {
    fn new(
        evidence: &'a StratumEvidence,
        homozygote_excess: &'a [f64],
        genotypes: &'a [(usize, usize)],
        config: &SsrFitConfig,
    ) -> Self {
        Self {
            evidence,
            homozygote_excess,
            genotypes,
            allele_span: config.allele_span,
            quadrature_points: config.quadrature_points,
            prepared: None,
            held_quadrature: None,
        }
    }

    /// The mean log-likelihood a tract at these parameters.
    fn score(&mut self, parameters: &Parameters) -> f64 {
        self.refresh(&parameters.slippage);
        let stale = !self
            .held_quadrature
            .as_ref()
            .is_some_and(|(spectrum, held, _)| {
                *held == parameters.concentration && *spectrum == parameters.length_spectrum
            });
        if stale {
            self.held_quadrature = Some((
                parameters.length_spectrum.clone(),
                parameters.concentration,
                dirichlet_points(
                    &parameters.length_spectrum,
                    parameters.concentration,
                    self.quadrature_points,
                    self.genotypes,
                ),
            ));
        }
        let (_, _, quadrature) = self.held_quadrature.as_ref().expect("built above");
        let prepared = self.prepared.as_ref().expect("refreshed above");
        // Across tracts in parallel, but summed back in tract order: a parallel float sum
        // reorders the additions run to run, and a fit's whole output is a difference between
        // one set of parameters and another.
        let per_tract: Vec<f64> = prepared
            .tracts
            .par_iter()
            .map(|tract| ln_tract(tract, quadrature, self.genotypes))
            .collect();
        per_tract.iter().sum::<f64>() / per_tract.len().max(1) as f64
    }

    fn refresh(&mut self, slippage: &[Slippage]) {
        if self
            .prepared
            .as_ref()
            .is_some_and(|held| held.slippage == slippage)
        {
            return;
        }
        let per_group_allele: Vec<Vec<Vec<f64>>> = slippage
            .iter()
            .map(|group| {
                (-self.allele_span..=self.allele_span)
                    .map(|allele| group.read_probabilities(allele, self.evidence.read_span))
                    .collect()
            })
            .collect();
        let tracts = self
            .evidence
            .tracts
            .par_iter()
            .map(|tract| {
                TractLikelihoods::of(
                    tract,
                    &per_group_allele,
                    self.genotypes,
                    self.homozygote_excess,
                )
            })
            .collect();
        self.prepared = Some(Prepared {
            slippage: slippage.to_vec(),
            tracts,
        });
    }
}

// ---------------------------------------------------------------------
// The integral over a tract's length frequencies
// ---------------------------------------------------------------------

/// The first primes, one a stick-breaking dimension. Thirteen allele classes need twelve.
const HALTON_BASES: [usize; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];

/// A Dirichlet with this mean and concentration, over a **fixed** low-discrepancy point set
/// pushed through the stick-breaking Beta quantiles.
///
/// **The points are fixed and the quantile map is continuous in the concentration**, so the
/// objective is a smooth function of it rather than a jittery one — the same uniforms are
/// reused at every value the climb tries, which is what makes this quadrature rather than
/// Monte Carlo and what stops the search chasing sampling noise.
fn dirichlet_points(
    spectrum: &[f64],
    concentration: f64,
    points: usize,
    genotypes: &[(usize, usize)],
) -> Quadrature {
    let classes = spectrum.len();
    let alpha: Vec<f64> = spectrum
        .iter()
        .map(|weight| (concentration * weight).max(1e-3))
        .collect();
    // **The two shapes are recomputed per point, and that measured as free.** Hoisting them and
    // their log-Beta to one per stick was tried on 2026-08-15 and moved the repeat-tract fit from
    // 19.7 s to 19.8 s — nothing, because the log-Beta is already computed once a bisection rather
    // than once a step, and a suffix sum over thirteen classes is twelve additions.
    //
    // **One buffer, filled in place, rather than one `Vec` a point.** Each point's stick-breaking
    // is untouched — same values, same order, same `beta_quantile` calls — and the points remain
    // independent, so this is exactly the same arithmetic written into a row of a flat buffer
    // instead of into a fresh allocation. What it removes is `points` heap blocks a rebuild.
    let mut frequencies = vec![0.0_f64; points * classes];
    frequencies
        .par_chunks_mut(classes)
        .enumerate()
        .for_each(|(point, piece)| {
            let mut remaining = 1.0;
            for stick in 0..classes - 1 {
                let a = alpha[stick];
                let b: f64 = alpha[stick + 1..].iter().sum();
                let uniform =
                    van_der_corput(point + 1, HALTON_BASES[stick.min(HALTON_BASES.len() - 1)]);
                let share = beta_quantile(uniform, a, b);
                piece[stick] = remaining * share;
                remaining *= 1.0 - share;
            }
            piece[classes - 1] = remaining;
        });
    // The genotype prior at every point, and where the homozygous pairs sit. Both depend on the
    // points alone, so this is the one place they are built.
    let mut independent = vec![0.0_f64; points * genotypes.len()];
    for (index, point) in frequencies.chunks_exact(classes).enumerate() {
        let row = &mut independent[index * genotypes.len()..][..genotypes.len()];
        for (slot, (first, second)) in genotypes.iter().enumerate() {
            row[slot] = if first == second {
                point[*first] * point[*first]
            } else {
                2.0 * point[*first] * point[*second]
            };
        }
    }
    // **Sized before it is filled**: thirteen classes give thirteen homozygous slots, and a
    // `collect` from a filter reaches that by three reallocations rather than one.
    let mut diagonal: Vec<usize> = Vec::with_capacity(classes);
    diagonal.extend(
        genotypes
            .iter()
            .enumerate()
            .filter(|(_, (first, second))| first == second)
            .map(|(slot, _)| slot),
    );

    Quadrature {
        frequencies,
        classes,
        ln_weight: -(points as f64).ln(),
        independent,
        diagonal,
    }
}

/// The `index`-th term of the van der Corput sequence in `base` — one coordinate of a Halton
/// point.
fn van_der_corput(mut index: usize, base: usize) -> f64 {
    let (mut fraction, mut result) = (1.0 / base as f64, 0.0);
    while index > 0 {
        result += (index % base) as f64 * fraction;
        index /= base;
        fraction /= base as f64;
    }
    result
}

// ---------------------------------------------------------------------
// Small numerics
// ---------------------------------------------------------------------

fn ln_sum_exp(values: &[f64]) -> f64 {
    let largest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !largest.is_finite() {
        return largest;
    }
    largest + values.iter().map(|v| (v - largest).exp()).sum::<f64>().ln()
}

fn ln_gamma(x: f64) -> f64 {
    // Lanczos, g = 7, n = 9 — the same coefficients the harnesses carry.
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut series = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        series += coefficient / (x + index as f64);
    }
    let t = x + 7.5;
    0.5 * std::f64::consts::TAU.ln() + (x + 0.5) * t.ln() - t + series.ln()
}

/// `ln B(a, b)`, the constant in front of the incomplete Beta.
///
/// **Symmetric in its two shapes**, up to the order the two subtractions happen in, which is why
/// one value serves both sides of the argument swap in
/// [`regularised_incomplete_beta_with`].
fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b)
}

/// `I_x(a, b)`, by its continued fraction — enough for a Beta quantile by bisection, with
/// `ln B(a, b)` handed in.
///
/// **The bisection above moves only `x`.** Recomputing `ln_beta` per step cost three `ln_gamma`
/// calls — each nine divisions, two logarithms, and below `x < 0.5` a sine and a recursion — on
/// every one of sixty steps, at every one of 256 quadrature points, for each of twelve
/// stick-breaking dimensions: 552,960 calls a quadrature build where twelve do.
fn regularised_incomplete_beta_with(x: f64, a: f64, b: f64, ln_beta: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    // **The swap is tested before the front factor is built, not after.** On this branch `front`
    // is dead, and building it costs two logarithms and an exponential for nothing.
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularised_incomplete_beta_with(1.0 - x, b, a, ln_beta);
    }
    let front = (a * x.ln() + b * (1.0 - x).ln() + ln_beta).exp() / a;
    let (mut f, mut c, mut d) = (1.0_f64, 1.0_f64, 0.0_f64);
    for index in 0..=200 {
        let m = index / 2;
        let numerator = if index == 0 {
            1.0
        } else if index % 2 == 0 {
            let m = m as f64;
            (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
        } else {
            let m = m as f64;
            -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
        };
        d = 1.0 + numerator * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        d = 1.0 / d;
        c = 1.0 + numerator / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        let step = c * d;
        f *= step;
        if (1.0 - step).abs() < 1e-12 {
            break;
        }
    }
    front * (f - 1.0)
}

/// The `p`-th quantile of `Beta(a, b)`, by bisection on the cumulative distribution.
fn beta_quantile(p: f64, a: f64, b: f64) -> f64 {
    beta_quantile_with(p, a, b, ln_beta(a, b))
}

/// How narrow the bracket has to get before the bisection stops.
///
/// **Sixty halvings of the unit interval reach 2⁻⁶⁰ ≈ 9 × 10⁻¹⁹, which is below what an `f64`
/// holds near 1**, so the last twenty of them moved the answer by nothing while each cost a whole
/// continued fraction of up to 201 terms. This stops at a width the answer can still carry — a
/// stick-breaking share is reported to four decimals and the quantities fitted from it to three.
const QUANTILE_TOLERANCE: f64 = 1e-12;

/// The same, with `ln B(a, b)` handed in — for a caller inverting the same distribution at many
/// probabilities, which is what a quadrature build does 256 times a stick.
fn beta_quantile_with(p: f64, a: f64, b: f64, ln_beta: f64) -> f64 {
    let (mut low, mut high) = (0.0_f64, 1.0_f64);
    for _ in 0..60 {
        if high - low < QUANTILE_TOLERANCE {
            break;
        }
        let middle = 0.5 * (low + high);
        if regularised_incomplete_beta_with(middle, a, b, ln_beta) < p {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}

fn logit(p: f64) -> f64 {
    let p = p.clamp(1e-9, 1.0 - 1e-9);
    (p / (1.0 - p)).ln()
}

fn expit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Golden-section on one coordinate, over a bracket of `span` either side of `start`.
fn climb_scalar(mut score: impl FnMut(f64) -> f64, start: f64, span: f64) -> f64 {
    const GOLDEN: f64 = 0.618_033_988_749_895;
    let (mut low, mut high) = (start - span, start + span);
    let (mut left, mut right) = (high - GOLDEN * (high - low), low + GOLDEN * (high - low));
    let (mut at_left, mut at_right) = (score(left), score(right));
    for _ in 0..16 {
        if at_left > at_right {
            high = right;
            right = left;
            at_right = at_left;
            left = high - GOLDEN * (high - low);
            at_left = score(left);
        } else {
            low = left;
            left = right;
            at_left = at_right;
            right = low + GOLDEN * (high - low);
            at_right = score(right);
        }
    }
    if at_left > at_right { left } else { right }
}

// ---------------------------------------------------------------------
// Reading the records
// ---------------------------------------------------------------------

/// Which stratum each kept STR locus is in, in the order the records index them.
///
/// **The record entry carries an index and nothing else** — no coordinates, no stratum — so
/// the order has to be rebuilt from the same kept-loci object the writer was given. That is
/// the loci of every stratum flattened together and sorted by contig, start and end.
///
/// A locus whose contig `contig_of` does not resolve is dropped, exactly as the writer drops
/// it, so the two lists stay the same length.
pub fn strata_of_kept_loci(
    loci: &CensusLoci,
    contig_of: &dyn Fn(&str) -> Option<ContigId>,
) -> Vec<Stratum> {
    let mut with_position: Vec<((u32, u64, u64), Stratum)> = loci
        .ssr()
        .iter_sorted()
        .into_iter()
        .flat_map(|((period, reference_repeats), segments)| {
            segments.iter().filter_map(move |segment| {
                contig_of(segment.chrom()).map(|contig| {
                    (
                        (contig.get(), segment.start(), segment.end()),
                        Stratum {
                            period,
                            reference_repeats,
                        },
                    )
                })
            })
        })
        .collect();
    with_position.sort_unstable_by_key(|(position, _)| *position);
    with_position
        .into_iter()
        .map(|(_, stratum)| stratum)
        .collect()
}

/// Gather one stratum's evidence a locus at a time, from a whole cohort's records.
///
/// `slippage_group_of` names, for each read group, which set of slippage numbers its reads are
/// drawn under. **One group per read group is the specified grain**; a run that knows several
/// read groups ran on one machine may pool them, and one that pools everything is saying it
/// cannot tell them apart.
///
/// Every sample must hold evidence for the same STR loci in the same order, which the
/// recording-terms check the cohort makes at its door has already refused to let fail silently.
///
/// **`strata` is one entry per kept tract, in genome order** — the stratum each tract is in.
/// The census stores a tract under an index within its own stratum, so this list is also what
/// says how many tracts each stratum holds, and a section of a different length means the loci
/// and the evidence were built from different selections.
///
/// # Panics
///
/// When a sample's section for a stratum is not as long as that stratum's share of `strata`.
pub fn gather_strata(
    cohort: &mut CohortCensusEvidence,
    strata: &[Stratum],
    slippage_group_of: &BTreeMap<ReadGroupId, u32>,
) -> Result<Vec<StratumEvidence>, CensusError> {
    let groups = slippage_group_of
        .values()
        .map(|group| *group as usize + 1)
        .max()
        .unwrap_or(1);
    let buckets = (2 * RECORDED_OFFSET_RANGE + 1) as usize;

    // How many tracts each stratum holds, which is the length its sections are built at.
    let mut tracts_in: BTreeMap<Stratum, usize> = BTreeMap::new();
    for stratum in strata {
        *tracts_in.entry(*stratum).or_insert(0) += 1;
    }
    let band: Vec<Stratum> = tracts_in.keys().copied().collect();
    let names: Vec<String> = cohort.sample_names().map(str::to_string).collect();

    // **The whole band at once, and that is this step's shape rather than the design's.** The
    // slippage fit borrows a thin stratum from its neighbours across the whole set
    // (`fit_strata`), so it is handed every stratum together; how many may be resident at once
    // is a measurement the fit specification owns (§11, questions 8 and 10) and not something
    // this call can decide.
    cohort.with_strata(&band, |lent| {
        // Each sample's tract sections, gathered by stratum. **One row a sample and not one
        // value**, because a stratum is fitted from every sample with reads in it at once.
        let by_sample: Vec<BTreeMap<Stratum, Vec<(ReadGroupId, &SsrEvidence)>>> = lent
            .iter()
            .zip(&names)
            .map(|(sections, name)| {
                let mut by_stratum: BTreeMap<Stratum, Vec<(ReadGroupId, &SsrEvidence)>> =
                    BTreeMap::new();
                for (group, stratum, records) in sections {
                    assert_eq!(
                        records.len(),
                        tracts_in.get(stratum).copied().unwrap_or(0),
                        "sample {} holds {} tracts at period {} and {} repeats where the \
                         selection has {}",
                        name,
                        records.len(),
                        stratum.period,
                        stratum.reference_repeats,
                        tracts_in.get(stratum).copied().unwrap_or(0)
                    );
                    by_stratum
                        .entry(*stratum)
                        .or_default()
                        .push((*group, records));
                }
                by_stratum
            })
            .collect();

        tracts_in
            .iter()
            .map(|(stratum, tracts)| {
                let mut evidence = StratumEvidence {
                    stratum: *stratum,
                    tracts: Vec::new(),
                    read_span: RECORDED_OFFSET_RANGE,
                    groups,
                    tracts_over_guard_threshold: 0,
                    reads_reaching_not_crossing: 0,
                    guard_reads: 0,
                    bases_compared: 0,
                    mismatching_bases: 0,
                };
                // **The reads that reached a tract and crossed none of it are counted per
                // stratum by the writer**, so they are read once here rather than accumulated a
                // locus at a time (spec §3). Every sample's every read group contributes its own
                // total.
                for sample in &by_sample {
                    for (_, records) in sample.get(stratum).into_iter().flatten() {
                        evidence.reads_reaching_not_crossing += records.covering_not_crossing();
                        // **The substitution rate's two counts are per section**, which is one
                        // read group's tracts for one stratum — the grain the rate is fitted at
                        // (`census::SsrEvidence::bases_compared`). A tract dropped by the guard
                        // below still contributed to them; the rate is a property of the
                        // sequence read, not of which tracts the slippage fit kept.
                        evidence.bases_compared += records.bases_compared();
                        evidence.mismatching_bases += records.differences().len() as u64;
                    }
                }

                for locus in 0..*tracts {
                    let mut reads = TractReads::default();
                    let mut over_guard = false;
                    for (sample_index, sample) in by_sample.iter().enumerate() {
                        let mut by_group: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
                        for (read_group, records) in sample.get(stratum).into_iter().flatten() {
                            if records.guard_is_over_threshold(locus) {
                                over_guard = true;
                            }
                            evidence.guard_reads += records
                                .guard()
                                .iter()
                                .filter(|entry| entry.locus as usize == locus)
                                .map(|entry| u64::from(entry.reads))
                                .sum::<u64>();
                            if records.state(locus) != SsrLocusState::Crossed {
                                continue;
                            }
                            let group = *slippage_group_of.get(read_group).unwrap_or(&0);
                            let counts = by_group.entry(group).or_insert_with(|| vec![0; buckets]);
                            let offsets = records.offsets(locus);
                            for (bucket, count) in counts.iter_mut().enumerate() {
                                *count +=
                                    u32::from(offsets.at(bucket as i32 - RECORDED_OFFSET_RANGE));
                            }
                        }
                        if !by_group.is_empty() {
                            reads.samples.push(SampleTractReads {
                                sample: sample_index as u32,
                                by_group: by_group.into_iter().collect(),
                            });
                        }
                    }
                    if over_guard {
                        evidence.tracts_over_guard_threshold += 1;
                        continue;
                    }
                    evidence.tracts.push(reads);
                }
                evidence
            })
            .collect()
    })
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Drawn strata: one generator, shared by the positive control and the bench
// ---------------------------------------------------------------------

/// A stratum drawn at a known truth, for anything that has to fit evidence it already knows
/// the answer to.
///
/// **Compiled under `cfg(test)` and under the `bench-fixtures` feature, and nowhere else.** Two
/// callers need a drawn stratum and they need the *same* one: this module's positive control
/// ([`fit_stratum`] must return the numbers a draw was made at) and `benches/ng_joint_fit_perf.rs`
/// (the fit must be timed on evidence with no CRAM behind it). A second generator would be a
/// second thing to keep agreeing, and a benchmark drawn differently from the oracle would be
/// timing a workload no test has ever checked.
///
/// Nothing here is production code: a release build without the feature compiles none of it.
#[cfg(any(test, feature = "bench-fixtures"))]
pub mod bench_fixtures {
    use super::*;

    /// A reproducible stream, the same one the harnesses use.
    struct Rng(u64);

    impl Rng {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        fn uniform(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1_u64 << 53) as f64
        }

        fn gamma(&mut self, shape: f64) -> f64 {
            if shape < 1.0 {
                let u = self.uniform().max(1e-300);
                return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
            }
            let d = shape - 1.0 / 3.0;
            let c = 1.0 / (9.0 * d).sqrt();
            loop {
                let x = self.normal();
                let v = (1.0 + c * x).powi(3);
                if v <= 0.0 {
                    continue;
                }
                let u = self.uniform().max(1e-300);
                if u.ln() < 0.5 * x * x + d - d * v + d * (v.ln()) {
                    return d * v;
                }
            }
        }

        fn normal(&mut self) -> f64 {
            let u1 = self.uniform().max(1e-300);
            let u2 = self.uniform();
            (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
        }

        fn dirichlet(&mut self, alpha: &[f64]) -> Vec<f64> {
            let draws: Vec<f64> = alpha.iter().map(|a| self.gamma(*a).max(1e-300)).collect();
            let total: f64 = draws.iter().sum();
            draws.into_iter().map(|d| d / total).collect()
        }

        fn categorical(&mut self, weights: &[f64]) -> usize {
            let mut u = self.uniform();
            for (index, weight) in weights.iter().enumerate() {
                u -= weight;
                if u <= 0.0 {
                    return index;
                }
            }
            weights.len() - 1
        }
    }

    /// Draw one stratum: `tracts` tracts, `samples` samples, `depth` reads a sample a tract.
    ///
    /// `span` is both the read span and the number of allele classes the spectrum must carry
    /// (`2 × span + 1`), so a caller fitting at [`SsrFitConfig::allele_span`] draws at the same
    /// span and hands the fit a spectrum of that length.
    ///
    /// **Every sample gets `depth` reads at every tract**, where a real cohort at three reads a
    /// position puts a read at a tract in a minority of its samples. So a drawn stratum of `n`
    /// samples costs more a tract than a recorded one of `n` samples — which is the right way
    /// round for a benchmark, and the wrong way round for extrapolating a wall time to a cohort.
    #[allow(
        clippy::too_many_arguments,
        reason = "the drawn stratum's own parameters, and the two axes CLAUDE.md §0 commits to \
                  (tracts, samples) are two of them"
    )]
    pub fn draw_stratum(
        slippage: Slippage,
        spectrum: &[f64],
        concentration: f64,
        homozygote_excess: f64,
        tracts: usize,
        samples: usize,
        depth: u32,
        span: i32,
        seed: u64,
    ) -> StratumEvidence {
        let classes = spectrum.len();
        let buckets = (2 * span + 1) as usize;
        let per_allele: Vec<Vec<f64>> = (0..classes)
            .map(|class| slippage.read_probabilities(class as i32 - span, span))
            .collect();
        let mut rng = Rng(seed);
        let mut drawn = Vec::with_capacity(tracts);
        for _ in 0..tracts {
            let alpha: Vec<f64> = spectrum.iter().map(|q| concentration * q).collect();
            let frequencies = rng.dirichlet(&alpha);
            let mut reads = TractReads::default();
            for sample in 0..samples {
                let first = rng.categorical(&frequencies);
                let second = if rng.uniform() < homozygote_excess {
                    first
                } else {
                    rng.categorical(&frequencies)
                };
                let mut counts = vec![0_u32; buckets];
                for _ in 0..depth {
                    let allele = if rng.uniform() < 0.5 { first } else { second };
                    counts[rng.categorical(&per_allele[allele])] += 1;
                }
                reads.samples.push(SampleTractReads {
                    sample: sample as u32,
                    by_group: vec![(0, counts)],
                });
            }
            drawn.push(reads);
        }
        StratumEvidence {
            stratum: Stratum {
                period: 2,
                reference_repeats: 10,
            },
            tracts: drawn,
            read_span: span,
            groups: 1,
            tracts_over_guard_threshold: 0,
            reads_reaching_not_crossing: 0,
            guard_reads: 0,
            // A drawn stratum has no sequence behind it, so it has no substitution rate. Zero
            // bases compared is what `substitution_rate()` returns `None` for, which is the
            // honest answer here rather than a fitted zero.
            bases_compared: 0,
            mismatching_bases: 0,
        }
    }

    /// The three-class spectrum every measurement on this path was made against: most
    /// chromosomes at the reference length, one repeat either side carrying most of the rest.
    pub fn spectrum_of(classes: usize) -> Vec<f64> {
        let middle = classes / 2;
        let mut spectrum: Vec<f64> = (0..classes)
            .map(|class| 0.55_f64.powi((class as i32 - middle as i32).abs()))
            .collect();
        normalise(&mut spectrum);
        spectrum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bench_fixtures::{draw_stratum, spectrum_of};

    /// Every read distribution is a distribution: it sums to one, whatever the allele's
    /// distance from the recorded range.
    #[test]
    fn a_read_distribution_sums_to_one_from_every_allele() {
        let slippage = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        for allele in -6..=6 {
            let total: f64 = slippage.read_probabilities(allele, 4).iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-12,
                "allele {allele} gave {total}, not one"
            );
        }
    }

    /// **The end bucket carries the marginal**: an allele sitting three repeats outside the
    /// recorded range puts nearly everything in the end bucket rather than nothing there.
    #[test]
    fn an_allele_outside_the_recorded_range_lands_in_the_end_bucket() {
        let slippage = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let reads = slippage.read_probabilities(-6, 4);
        assert!(
            reads[0] > 0.9,
            "the shortest bucket took {}, not nearly everything",
            reads[0]
        );
        assert!(reads[8] < 0.01, "the far end took {}", reads[8]);
    }

    /// **The positive control.** A stratum drawn with a known truth comes back with it: this
    /// is the run that says the estimator has power, and without it a clean-looking answer
    /// cannot be told from one with no information in it.
    #[test]
    fn a_drawn_stratum_returns_the_numbers_it_was_drawn_with() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let spectrum = spectrum_of(3);
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 1_500, 20, 6, 1, 11);
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let fitted = fit_stratum(&evidence, &vec![0.4; 20], &config).expect("reads were drawn");
        let slippage = fitted.slippage[0].expect("the one group has reads");

        assert!(
            (slippage.level - truth.level).abs() / truth.level < 0.10,
            "slippage level {} against a truth of {}",
            slippage.level,
            truth.level
        );
        assert!(
            (slippage.shorter_share - truth.shorter_share).abs() < 0.05,
            "shorter-share {} against a truth of {}",
            slippage.shorter_share,
            truth.shorter_share
        );
        assert!(
            (fitted.concentration - 0.5).abs() / 0.5 < 0.25,
            "concentration {} against a truth of 0.5",
            fitted.concentration
        );
    }

    /// A stratum fitted on its own tracts, built directly so the smoothing can be exercised
    /// without paying for the climb.
    fn fitted_at(period: u8, repeats: u64, level: f64, reads: u64) -> StratumOutcome {
        StratumOutcome::Fitted(Box::new(StratumFit {
            stratum: Stratum {
                period,
                reference_repeats: repeats,
            },
            slippage: vec![Some(Slippage {
                level,
                shorter_share: 0.7,
                fall_off: 0.3,
            })],
            length_spectrum: vec![1.0],
            concentration: 0.6,
            log_likelihood_a_tract: -1.0,
            tracts_fitted: 500,
            borrowed: Vec::new(),
            converged: true,
            tracts_of_its_own: 500,
            reads_crossing: reads,
            level_provenance: vec![Some(LevelProvenance {
                source: LevelSource::Cell,
                curve: None,
                reach: None,
                slipped_reads: Some(level * reads as f64),
            })],
            shares_provenance: vec![Some(SharesProvenance::own(level * reads as f64))],
        }))
    }

    /// Draw the level's curves and re-emit every level through them, as `fit_strata` does.
    fn smooth_levels(outcomes: &mut [StratumOutcome], config: &SsrFitConfig) {
        let curves = draw_a_curve_a_period(outcomes, config);
        smooth_levels_across_repeat_count(outcomes, &curves, config);
    }

    fn level_of(outcome: &StratumOutcome) -> f64 {
        outcome.slippage()[0].expect("a group with numbers").level
    }

    fn provenance_of(outcome: &StratumOutcome) -> LevelProvenance {
        outcome.level_provenance()[0].expect("a group with numbers")
    }

    /// **The property the whole change rests on: smoothing moves the level and nothing else.**
    /// A stratum's direction split, fall-off, spectrum and concentration are its own answer and
    /// must be the same number whether curves are drawn or not.
    #[test]
    fn smoothing_moves_the_level_and_leaves_every_other_number_alone() {
        let cells: Vec<StratumOutcome> = (8..=20)
            .map(|repeats| {
                // A straight line with one cell knocked 40% off it, so smoothing has work to do.
                let on_the_line = 0.005 * repeats as f64 - 0.035;
                let level = if repeats == 14 {
                    on_the_line * 0.6
                } else {
                    on_the_line
                };
                fitted_at(1, repeats, level, 200_000)
            })
            .collect();

        let mut unsmoothed = cells.clone();
        let mut smoothed = cells.clone();
        let config = SsrFitConfig::default();
        smooth_levels(&mut smoothed, &config);

        // The "off" arm is the same call with the switch down; it must change nothing at all.
        let off = SsrFitConfig {
            curve: SlippageCurveConfig {
                draw_curves: false,
                ..SlippageCurveConfig::default()
            },
            ..SsrFitConfig::default()
        };
        if off.curve.draw_curves {
            smooth_levels(&mut unsmoothed, &off);
        }
        assert_eq!(unsmoothed, cells, "the switch down must move nothing");

        for (before, after) in cells.iter().zip(&smoothed) {
            let (StratumOutcome::Fitted(before), StratumOutcome::Fitted(after)) = (before, after)
            else {
                panic!("both arms are fitted");
            };
            let (was, now) = (
                before.slippage[0].expect("fitted"),
                after.slippage[0].expect("fitted"),
            );
            assert_eq!(was.shorter_share, now.shorter_share);
            assert_eq!(was.fall_off, now.fall_off);
            assert_eq!(before.concentration, after.concentration);
            assert_eq!(before.length_spectrum, after.length_spectrum);
            assert_eq!(before.tracts_of_its_own, after.tracts_of_its_own);
        }

        // The knocked-down cell is pulled back toward its neighbours, and every cell records
        // that a curve had a say.
        let knocked = smoothed
            .iter()
            .find(|outcome| matches!(outcome, StratumOutcome::Fitted(fit) if fit.stratum.reference_repeats == 14))
            .expect("the cell at fourteen repeats");
        let on_the_line = 0.005 * 14.0 - 0.035;
        assert!(
            level_of(knocked) > on_the_line * 0.6,
            "the knocked-down cell should be pulled up, not left at {}",
            level_of(knocked)
        );
        let provenance = provenance_of(knocked);
        assert!(matches!(provenance.source, LevelSource::Blend { .. }));
        assert_eq!(provenance.reach, Some(CurveReach::Inside));
        let curve = provenance.curve.expect("a curve stood behind it");
        assert_eq!(curve.cells, 13);
        assert_eq!((curve.fitted_from, curve.fitted_to), (8, 20));
    }

    /// A period with too few strata to draw a curve keeps every level exactly as fitted.
    #[test]
    fn a_period_below_the_cell_floor_is_left_entirely_alone() {
        let cells: Vec<StratumOutcome> = (8..=10)
            .map(|repeats| fitted_at(1, repeats, 0.002 * repeats as f64, 100_000))
            .collect();
        let mut smoothed = cells.clone();
        smooth_levels(&mut smoothed, &SsrFitConfig::default());
        assert_eq!(smoothed, cells);
    }

    /// **A stratum that borrowed must not feed the curve**, or the curve is fitted to its own
    /// output. It still reads the curve — it is the borrowing this replaces.
    #[test]
    fn a_borrowed_stratum_reads_the_curve_but_does_not_feed_it() {
        let mut cells: Vec<StratumOutcome> = (8..=20)
            .map(|repeats| fitted_at(1, repeats, 0.005 * repeats as f64 - 0.035, 200_000))
            .collect();
        // One more stratum, wildly off the line, marked as having borrowed.
        let mut borrower = fitted_at(1, 21, 0.9, 200_000);
        if let StratumOutcome::Fitted(fit) = &mut borrower {
            fit.borrowed = vec![20];
        }
        cells.push(borrower);

        let mut smoothed = cells.clone();
        smooth_levels(&mut smoothed, &SsrFitConfig::default());

        let last = smoothed.last().expect("the borrower");
        let curve = provenance_of(last).curve.expect("a curve");
        assert_eq!(
            curve.cells, 13,
            "the borrower's own level must not be one of the cells behind the curve"
        );
        assert_eq!((curve.fitted_from, curve.fitted_to), (8, 20));
        // What the borrower's own level then becomes is
        // `a_borrowed_stratum_takes_the_curves_level_and_keeps_the_pooled_shares`; what this
        // test owns is that its level never reached the cells the curve was drawn through.
        assert!(
            !cells.iter().take(13).any(|outcome| level_of(outcome) > 0.5),
            "no contributing cell carries the borrower's level"
        );
    }

    /// **The narrowing B4 is for: a stratum that borrowed takes the curve's level outright and
    /// keeps the pooled fit's shares.** Its pooled level is its neighbours' — the very thing the
    /// curve replaces — so weighing it against the curve would be weighing the curve against a
    /// blurred copy of itself.
    #[test]
    fn a_borrowed_stratum_takes_the_curves_level_and_keeps_the_pooled_shares() {
        let mut cells: Vec<StratumOutcome> = (8..=20)
            .map(|repeats| fitted_at(1, repeats, 0.005 * repeats as f64 - 0.035, 200_000))
            .collect();
        // A borrower at 21 repeats, carrying a pooled level far off the line and pooled shares.
        let mut borrower = fitted_at(1, 21, 0.9, 200_000);
        if let StratumOutcome::Fitted(fit) = &mut borrower {
            fit.borrowed = vec![20];
            let slippage = fit.slippage[0].as_mut().expect("a fitted group");
            slippage.shorter_share = 0.81;
            slippage.fall_off = 0.42;
            fit.level_provenance[0]
                .as_mut()
                .expect("a fitted group")
                .slipped_reads = None;
        }
        cells.push(borrower);

        let mut smoothed = cells.clone();
        smooth_levels(&mut smoothed, &SsrFitConfig::default());

        let last = smoothed.last().expect("the borrower");
        let StratumOutcome::Fitted(fit) = last else {
            panic!("the borrower is fitted");
        };
        let curve = provenance_of(last).curve.expect("a curve stood behind it");

        // The level is the curve's, whole — not a blend with the pooled 0.9.
        assert_eq!(provenance_of(last).source, LevelSource::Curve);
        assert_eq!(level_of(last), curve.level_at(21));
        assert_eq!(provenance_of(last).reach, Some(CurveReach::AboveFitted));
        assert_eq!(provenance_of(last).slipped_reads, None);

        // The two shares are untouched: borrowing still supplies them.
        let slippage = fit.slippage[0].expect("a fitted group");
        assert_eq!(slippage.shorter_share, 0.81);
        assert_eq!(slippage.fall_off, 0.42);
    }

    /// Re-emitting the levels leaves a stratum with no fit of its own untouched. Turning it into
    /// a stratum furnished from its period's curves is a separate step — `derive_thin_strata` —
    /// and it needs all three curves, not the level's alone.
    #[test]
    fn smoothing_the_levels_leaves_a_stratum_with_no_fit_alone() {
        let mut cells: Vec<StratumOutcome> = (8..=20)
            .map(|repeats| fitted_at(1, repeats, 0.005 * repeats as f64 - 0.035, 200_000))
            .collect();
        cells.push(StratumOutcome::Refused {
            stratum: Stratum {
                period: 1,
                reference_repeats: 21,
            },
            tracts: 3,
            reason: StratumRefusal::BelowTheFloor {
                tracts: 3,
                floor: 50,
            },
        });
        let mut smoothed = cells.clone();
        smooth_levels(&mut smoothed, &SsrFitConfig::default());
        assert_eq!(smoothed.last(), cells.last());
    }

    /// **The whole point of the change, end to end: a stratum too thin to be fitted alone comes
    /// back with all three of its slippage numbers drawn from its period's curves.** Standing
    /// alone it is refused and gets nothing at all.
    #[test]
    fn a_stratum_too_thin_to_stand_alone_ends_up_with_all_three_numbers_from_curves() {
        let spectrum = spectrum_of(3);
        // Five fat strata whose levels rise with repeat count, and one far too thin to be fitted
        // on its own tracts.
        let mut strata = Vec::new();
        for (repeats, tracts, level, seed) in [
            (8_u64, 120_usize, 0.06, 3_u64),
            (9, 120, 0.08, 5),
            (10, 120, 0.10, 7),
            (11, 120, 0.12, 9),
            (12, 120, 0.14, 11),
            (13, 8, 0.16, 13),
        ] {
            let truth = Slippage {
                level,
                shorter_share: 0.83,
                fall_off: 0.25,
            };
            let mut evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, tracts, 8, 6, 1, seed);
            evidence.stratum = Stratum {
                period: 2,
                reference_repeats: repeats,
            };
            strata.push(evidence);
        }
        let base = SsrFitConfig {
            allele_span: 1,
            max_rounds: 1,
            refusal_floor: 50,
            starting_points: vec![StartingPoint {
                slippage_level: 0.10,
                concentration: 3.0,
            }],
            ..SsrFitConfig::default()
        };

        let outcomes = fit_strata(&strata, &vec![0.4; 8], &base);

        // The five that clear the refusal floor are fitted from their own tracts and no others'.
        for outcome in outcomes.iter().take(5) {
            let StratumOutcome::Fitted(fit) = outcome else {
                panic!("a hundred and twenty tracts clears the floor");
            };
            assert!(fit.borrowed.is_empty());
            assert_eq!(fit.tracts_fitted, fit.tracts_of_its_own);
        }

        // **Eight tracts is far below the refusal floor, so nothing was fitted from them — and
        // the stratum still comes back with a complete set of numbers.**
        let StratumOutcome::Derived(thin) = &outcomes[5] else {
            panic!(
                "a stratum with reads is furnished, not refused: {:?}",
                outcomes[5]
            );
        };
        assert_eq!(thin.stratum.reference_repeats, 13);
        assert_eq!(thin.tracts_of_its_own, 8);
        assert!(thin.reads_crossing > 0);

        // Its level is the curve's, held at the top of the range the curve was drawn over.
        let provenance = thin.level_provenance[0].expect("the only group put reads here");
        assert_eq!(provenance.source, LevelSource::Curve);
        assert_eq!(provenance.slipped_reads, None, "nothing was fitted here");
        let curve = provenance.curve.expect("its period has a curve");
        assert_eq!(
            curve.cells, 5,
            "only the five strata fitted on their own tracts feed it"
        );
        assert_eq!((curve.fitted_from, curve.fitted_to), (8, 12));
        assert_eq!(provenance.reach, Some(CurveReach::AboveFitted));

        let numbers = thin.slippage[0].expect("a furnished group");
        assert_eq!(numbers.level, curve.level_at(13));
        assert!(numbers.level > 0.0 && numbers.level < 1.0);

        // **Its two shares are its period's curves, the same treatment the level gets.** Nothing
        // was fitted here, so there is no own answer to blend with and the curve is taken whole.
        let shares = thin.shares_provenance[0].expect("a furnished group");
        assert_eq!(shares.slipped_reads, None, "nothing was fitted here");
        assert_eq!(shares.shorter_share.source, ShareSource::Curve);
        assert_eq!(shares.fall_off.source, ShareSource::Curve);

        let split_curve = shares.shorter_share.curve.expect("its period has a curve");
        let fall_off_curve = shares.fall_off.curve.expect("its period has a curve");
        assert_eq!(numbers.shorter_share, split_curve.share_at(13));
        assert_eq!(numbers.fall_off, fall_off_curve.share_at(13));
        assert_eq!(split_curve.source, ShareCurveSource::ThisPeriod);
        assert_eq!(
            split_curve.strata, 5,
            "only the five strata fitted on their own tracts feed it"
        );

        // Both shares stay proportions, and near the truth the strata were drawn from.
        assert!(numbers.shorter_share > 0.0 && numbers.shorter_share < 1.0);
        assert!(numbers.fall_off > 0.0 && numbers.fall_off < 1.0);
    }

    /// **The refusal floor is the measured one**, and the test carries the number so that moving
    /// it is a deliberate act rather than a typo. Below it a stratum's own fit is refused; it can
    /// still be furnished from its period's curves.
    #[test]
    fn nothing_is_fitted_below_eight_tracts_by_default() {
        assert_eq!(SsrFitConfig::default().refusal_floor, DEFAULT_REFUSAL_FLOOR);
        assert_eq!(DEFAULT_REFUSAL_FLOOR, 8);
    }

    /// A stratum with no reads at all is refused rather than fitted, and it is refused by name.
    #[test]
    fn a_stratum_with_no_reads_is_refused_and_not_fitted() {
        let evidence = StratumEvidence {
            stratum: Stratum {
                period: 2,
                reference_repeats: 9,
            },
            tracts: vec![TractReads::default(); 20],
            read_span: 4,
            groups: 1,
            tracts_over_guard_threshold: 0,
            reads_reaching_not_crossing: 40,
            guard_reads: 0,
            bases_compared: 0,
            mismatching_bases: 0,
        };
        let outcomes = fit_strata(&[evidence], &[0.4], &SsrFitConfig::default());
        assert!(matches!(
            &outcomes[0],
            StratumOutcome::Refused {
                reason: StratumRefusal::NoSpanningReads,
                ..
            }
        ));
    }

    /// **Every stratum's two shares are re-emitted through its period's curves, and how far each
    /// moves is set by its own evidence.** The rule this replaced was a gate: a stratum either
    /// measured its own shares on 4,000 slipped reads or took one neighbour's whole.
    #[test]
    fn a_stratum_departs_from_its_periods_share_curve_by_how_much_evidence_it_has() {
        let spectrum = spectrum_of(3);
        let truth = Slippage {
            level: 0.10,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        // Five strata at one period, the last of them holding a twentieth of the others' tracts.
        let mut strata = Vec::new();
        for (repeats, tracts, seed) in [
            (8_u64, 400_usize, 3_u64),
            (9, 400, 5),
            (10, 400, 7),
            (11, 400, 9),
            (12, 60, 11),
        ] {
            let mut evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, tracts, 8, 6, 1, seed);
            evidence.stratum = Stratum {
                period: 2,
                reference_repeats: repeats,
            };
            strata.push(evidence);
        }
        let config = SsrFitConfig {
            allele_span: 1,
            max_rounds: 1,
            refusal_floor: 50,
            starting_points: vec![StartingPoint {
                slippage_level: 0.10,
                concentration: 3.0,
            }],
            ..SsrFitConfig::default()
        };
        let outcomes = fit_strata(&strata, &vec![0.4; 8], &config);

        let shares_of = |index: usize| match &outcomes[index] {
            StratumOutcome::Fitted(fit) => fit.shares_provenance[0].expect("a fitted group"),
            other => panic!("stratum {index} is fitted: {other:?}"),
        };

        // Every one of them is a blend of its own answer and its period's curve, and each says
        // how much of the curve it took.
        for index in 0..5 {
            let shares = shares_of(index);
            assert!(
                matches!(shares.shorter_share.source, ShareSource::Blend { .. }),
                "stratum {index} came back as {:?}",
                shares.shorter_share.source
            );
            assert!(shares.slipped_reads.expect("its own fit") > 0.0);
        }

        // **The thin one takes more of the curve than the fat ones**, because it holds its own
        // answer less precisely — and that is the whole of the rule that replaced the gate.
        let fat = shares_of(0).shorter_share.source.curve_weight();
        let thin = shares_of(4).shorter_share.source.curve_weight();
        assert!(
            thin > fat,
            "the thin stratum took {thin:.3} of its curve and the fat one {fat:.3}"
        );
    }

    /// **A stratum that read another's tracts consumes the share curves and does not feed
    /// them**, exactly as it does for the level: its shares are the pooled set's, so letting them
    /// into the fit would be fitting a curve to a neighbour's answer.
    #[test]
    fn a_stratum_that_borrowed_reads_the_share_curves_but_does_not_feed_them() {
        let mut cells: Vec<StratumOutcome> = (8..=12)
            .map(|repeats| fitted_at(2, repeats, 0.02 * (repeats - 6) as f64, 40_000))
            .collect();
        let without_it = draw_share_curves_a_period(&cells, &SsrFitConfig::default());

        let mut pooled = fitted_at(2, 13, 0.16, 40_000);
        if let StratumOutcome::Fitted(fit) = &mut pooled {
            fit.borrowed = vec![12];
            if let Some(slippage) = fit.slippage[0].as_mut() {
                slippage.shorter_share = 0.05;
                slippage.fall_off = 0.95;
            }
        }
        cells.push(pooled);
        let with_it = draw_share_curves_a_period(&cells, &SsrFitConfig::default());
        assert_eq!(
            without_it.get(&(2, 0)),
            with_it.get(&(2, 0)),
            "a stratum that borrowed moved its period's share curves"
        );

        // And it takes the curves whole, with no own answer weighed against them.
        let config = SsrFitConfig::default();
        smooth_shares_across_repeat_count(&mut cells, &with_it, &config);
        let StratumOutcome::Fitted(fit) = &cells[5] else {
            panic!("the pooled stratum is still a fit");
        };
        let shares = fit.shares_provenance[0].expect("a fitted group");
        assert_eq!(shares.shorter_share.source, ShareSource::Curve);
        assert_eq!(shares.fall_off.source, ShareSource::Curve);
        assert_eq!(shares.slipped_reads, None);
    }

    /// **A curve is fitted to what the strata measured, never to what a curve emitted.** Drawing
    /// the level's curve after the levels had been blended would fit the second curve to the
    /// first one's output; `fit_strata` draws all three before either blend runs, so the curve a
    /// thin stratum is furnished from is the one the fitted strata's own levels give.
    #[test]
    fn the_curve_a_thin_stratum_reads_is_fitted_to_unblended_levels() {
        let mut cells: Vec<StratumOutcome> = (8..=12)
            .map(|repeats| fitted_at(2, repeats, 0.02 * (repeats - 6) as f64, 40_000))
            .collect();
        let before = draw_a_curve_a_period(&cells, &SsrFitConfig::default());

        let config = SsrFitConfig::default();
        let shares = draw_share_curves_a_period(&cells, &config);
        smooth_levels_across_repeat_count(&mut cells, &before, &config);
        smooth_shares_across_repeat_count(&mut cells, &shares, &config);

        // Refitting now — which is what the code did before the three curves were hoisted out —
        // gives a different curve, because every level it reads has already been smoothed once.
        let after = draw_a_curve_a_period(&cells, &config);
        let line = |curves: &BTreeMap<u8, PeriodCurves>| {
            let period = curves.get(&2).expect("period 2 has a curve");
            let curve = period.by_group[0].expect("the only group has a line");
            (curve.intercept, curve.slope)
        };
        assert_ne!(
            line(&before),
            line(&after),
            "a second round of smoothing should move the curve, or this test proves nothing"
        );
    }
}
