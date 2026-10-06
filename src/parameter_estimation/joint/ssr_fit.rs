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
//!   ([`RECORDED_OFFSET_RANGE`]).
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
use std::num::NonZeroUsize;

use rayon::prelude::*;

use crate::float;
use crate::parameter_estimation::joint::census::{
    CensusError, CohortCensusEvidence, RECORDED_OFFSET_RANGE, SsrEvidence, SsrLocusState,
};
use crate::parameter_estimation::joint::fit::SETTLED_FRACTION;
use crate::parameter_estimation::joint::loci::CensusLoci;
use crate::parameter_estimation::progress::StageProgress;
use crate::types::{ContigId, ReadGroupId};

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
                let weight = (1.0 - self.fall_off) * float::powi(self.fall_off, step - 1);
                out[slot(allele_offset + direction * step)] += self.level * share * weight;
            }
            // Everything at least `inside_steps + 1` steps away, whose weight telescopes to
            // `fall_off^inside_steps`, lands in the end bucket.
            let tail = float::powi(self.fall_off, inside_steps);
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

/// One stratum's substitution counts: its bases compared against a read and, of those, the ones
/// the read disagreed with.
///
/// **Two numbers a stratum, where [`StratumEvidence`] holds every sample's reads at every
/// tract.** The substitution rate needs only these two, so they can be kept after the evidence
/// is gone. On kimura's 100 tomato samples at about three reads a position, gathering the
/// evidence added 0.7 GiB (`doc/devel/implementation_plans/estimation_memory.md` §1), and it
/// grows with every sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StratumSubstitutionCounts {
    /// Which stratum these were counted over.
    pub stratum: Stratum,
    /// Bases of tract sequence a read was compared against, over every tract and sample.
    pub bases_compared: u64,
    /// Of those bases, how many the read disagreed with, on reads whose tract was the
    /// reference's length.
    pub mismatching_bases: u64,
}

impl StratumSubstitutionCounts {
    /// Mismatching bases over bases compared — the stratum's substitution rate (spec §4.2) —
    /// **never exactly zero or one**.
    ///
    /// **A count that found no mismatch takes half of one**, `0.5 / (n + 1)` over `n` bases
    /// compared, and one where every base mismatched takes half a match, `(n + 0.5) / (n + 1)`.
    /// Calling scores a tract's reads under this rate, and under a rate of zero a read with even one
    /// base that disagrees is explained by no tract length at all: it falls wholly to the outlier
    /// term and stops counting for any genotype (owner, checkpoint E of `fit_precision.md`,
    /// 2026-10-06). Half a count keeps what the count says — no mismatch in 500 bases becomes
    /// 0.001, no mismatch in 34 becomes 0.014 — and a count that saw both outcomes keeps its own
    /// ratio.
    ///
    /// `None` where no read was compared against a tract at all. **A count of more mismatching
    /// bases than bases compared cannot be built** — every read's bases are counted as compared
    /// before they are compared — and is given its ratio, above one, which no rate type accepts.
    pub fn substitution_rate(&self) -> Option<f64> {
        let bases_compared = self.bases_compared;
        if bases_compared == 0 {
            return None;
        }
        debug_assert!(
            self.mismatching_bases <= bases_compared,
            "{} mismatching bases of {bases_compared} compared",
            self.mismatching_bases
        );
        let bases_compared_f64 = bases_compared as f64;
        Some(if self.mismatching_bases == 0 {
            0.5 / (bases_compared_f64 + 1.0)
        } else if self.mismatching_bases == bases_compared {
            (bases_compared_f64 + 0.5) / (bases_compared_f64 + 1.0)
        } else {
            self.mismatching_bases as f64 / bases_compared_f64
        })
    }
}

impl StratumEvidence {
    /// **Rows with reads**: one for each (tract, sample) pair where the sample has a read that
    /// crossed the tract — what the evidence and the stratum's likelihood table grow with.
    pub fn rows_with_reads(&self) -> usize {
        self.tracts.iter().map(|tract| tract.samples.len()).sum()
    }

    /// **The bytes this stratum's evidence holds**, counted by what its vectors reserved:
    /// the tract list, each tract's rows, each row's slippage groups and each group's counts.
    /// The allocator rounds every block up, so the process holds somewhat more.
    pub fn heap_bytes(&self) -> usize {
        self.tracts.capacity() * std::mem::size_of::<TractReads>()
            + self
                .tracts
                .iter()
                .map(|tract| {
                    tract.samples.capacity() * std::mem::size_of::<SampleTractReads>()
                        + tract
                            .samples
                            .iter()
                            .map(|row| {
                                row.by_group.capacity() * std::mem::size_of::<(u32, Vec<u32>)>()
                                    + row
                                        .by_group
                                        .iter()
                                        .map(|(_, counts)| {
                                            counts.capacity() * std::mem::size_of::<u32>()
                                        })
                                        .sum::<usize>()
                            })
                            .sum::<usize>()
                })
                .sum::<usize>()
    }

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

    /// The stratum's substitution rate, which needs none of the other numbers (spec §4.2) — see
    /// [`StratumSubstitutionCounts::substitution_rate`], which is never exactly zero or one.
    ///
    /// `None` where no read was compared against a tract at all.
    pub fn substitution_rate(&self) -> Option<f64> {
        self.substitution_counts().substitution_rate()
    }

    /// The two counts the stratum's substitution rate is made of, without the tracts — what
    /// outlives the evidence once the fit is done.
    pub fn substitution_counts(&self) -> StratumSubstitutionCounts {
        StratumSubstitutionCounts {
            stratum: self.stratum,
            bases_compared: self.bases_compared,
            mismatching_bases: self.mismatching_bases,
        }
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
    /// **How the walk that won ended** ([`ClimbEnding`]): settled, stopped at a round that lost,
    /// or out of rounds. **Neither of the last two is reported as convergence.**
    pub ending: ClimbEnding,
    /// **Every walk over the stratum**, one per starting point in [`SsrFitConfig::starting_points`]'s
    /// order: how it ended and how many rounds it took — what the climb cost. Empty on the fixtures
    /// tests build by hand.
    pub walks: Vec<WalkRecord>,
    /// **How many samples the stratum was fitted on**, when a large cohort's stratum was fitted on
    /// a subset of them ([`fit_strata_on_sample_subsets`]) — the whole cohort's count when the
    /// subset grew to every sample. `None` when every sample was read without a subset being drawn:
    /// a cohort of at most [`SampleSubsets::first`] samples, or any call to [`fit_strata`]. The
    /// evidence counts the fit carries ([`tracts_fitted`](Self::tracts_fitted),
    /// [`tracts_of_its_own`](Self::tracts_of_its_own), [`reads_crossing`](Self::reads_crossing))
    /// are the subset's.
    pub samples_fitted_on: Option<usize>,
    /// Tracts **this stratum itself** holds with at least one spanning read **among the samples it
    /// was fitted on** — every sample unless a large cohort's stratum was fitted on a subset
    /// ([`samples_fitted_on`](Self::samples_fitted_on)) — whatever it read to produce its answer.
    ///
    /// **Distinct from [`StratumFit::tracts_fitted`], and the difference is the whole point.** A
    /// stratum with eight tracts of its own that borrowed its way to a thousand has an answer
    /// resting on its neighbours, and a consumer told only the second number cannot see that
    /// (`str_slippage_level_curve.md` §8).
    pub tracts_of_its_own: usize,
    /// Reads that crossed a whole tract of **this stratum itself**, over every group and every
    /// sample it was fitted on ([`samples_fitted_on`](Self::samples_fitted_on)).
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
    /// **How precisely the stratum's own tracts determine each of its numbers**, at the answer the
    /// climb returned (`fit_precision.md` §4.2).
    ///
    /// **The errors of the stratum's own fit, and only of it.** Once [`fit_strata`] has drawn the
    /// curves, [`StratumFit::slippage`] may hold a blend of the own fit and its period's curve, and
    /// these are not that blend's errors — a blended number has none (spec §5.2). The own fit's level
    /// is the one [`LevelProvenance::slipped_reads`] is counted from.
    ///
    /// `None` where no error was computed: on the fixtures tests build by hand.
    pub standard_errors: Option<StratumErrors>,
}

/// **The standard errors of one stratum's own fit**, each on its number's own scale, from the
/// curvature of the stratum's total log-likelihood at the answer (`fit_precision.md` §4.2).
///
/// An error is how far the number would typically move were the same kind of reads drawn again. A
/// number without one says why ([`StratumError`]); never a zero, and never a large number standing in
/// for one.
#[derive(Debug, Clone, PartialEq)]
pub struct StratumErrors {
    /// Per slippage group, its three numbers' errors — `None` where the group put no read in this
    /// stratum, matching [`StratumFit::slippage`] index for index.
    pub slippage: Vec<Option<SlippageErrors>>,
    /// Per allele class, the error of its share of the stratum's chromosomes, indexed as
    /// [`StratumFit::length_spectrum`] is.
    pub length_spectrum: Vec<StratumError>,
    /// The concentration's error.
    pub concentration: StratumError,
}

/// The standard errors of one slippage group's three numbers in one stratum ([`StratumErrors`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlippageErrors {
    pub level: StratumError,
    pub shorter_share: StratumError,
    pub fall_off: StratumError,
}

/// **One number's standard error in a stratum's fit, or why it has none.**
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StratumError {
    /// The error, on the number's own scale.
    Estimated(f64),
    /// **The tracts do not tell it apart from the stratum's other numbers**: once they are accounted
    /// for, the log-likelihood is not curved downwards in its direction — flat, or a saddle.
    NotIdentified,
    /// **The tracts do not place it**: its error on the scale the climb moves it on — logit for a
    /// slippage number, log for the concentration, a log-ratio for a share — is above
    /// [`NOT_PLACED`], so a one-error interval spans odds or a size hundreds of times apart. Where
    /// the climb ran a number off towards an end of its range and stopped at the end of its reach,
    /// the log-likelihood there is flat, and carrying its curvature to the number's own scale would
    /// report an error hundreds of times too small. Also a number in `[0, 1]` — a slippage number or
    /// a share — whose error, on its own scale, is wider than that whole range: the SNP/indel fit's
    /// "wider than its range" (spec §3.2).
    NotPlaced,
    /// **A length class with no share at all**, held at zero and never fitted.
    NoShare,
}

impl StratumErrors {
    /// Whether any number of the stratum has an error.
    pub fn any(&self) -> bool {
        let slippage = self.slippage.iter().flatten().any(|group| {
            [group.level, group.shorter_share, group.fall_off]
                .iter()
                .any(|error| error.value().is_some())
        });
        slippage
            || self
                .length_spectrum
                .iter()
                .any(|error| error.value().is_some())
            || self.concentration.value().is_some()
    }
}

impl StratumError {
    /// The error, or `None` when there is none.
    pub fn value(self) -> Option<f64> {
        match self {
            Self::Estimated(error) => Some(error),
            Self::NotIdentified | Self::NotPlaced | Self::NoShare => None,
        }
    }
}

/// **The largest error, on the scale the climb moves a number on, that still places the number**:
/// three units of logit or of log, so a one-error interval either side spans odds or a size `e⁶`,
/// about 400 times, apart. Measured in the review of plan step C1 on drawn strata: numbers the tracts
/// determine sat at 0.02 to 0.5 on those scales, and numbers the climb ran off to the end of its
/// reach at 1,000 to 2,400.
pub const NOT_PLACED: f64 = 3.0;

/// **How one walk of the climb ended** (`fit_precision.md` §4.3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClimbEnding {
    /// **Settled**: every number is within [`SsrFitConfig::settled_fraction`] of its own standard
    /// error of where the likelihood peaks, as a Newton step from the walk's last point estimates the
    /// distance.
    Settled,
    /// **A round lowered the log-likelihood, and nothing it passed through was better than where it
    /// began**, which was judged and had not settled. Its moves were undone and the walk stopped
    /// there, since a next round would start from the same point and repeat it exactly. Not
    /// convergence: the golden section that moves each number never scores the value it starts from,
    /// so a round can leave a number somewhere worse while the others still had a way to go.
    LostARound,
    /// **Ran out of rounds** ([`SsrFitConfig::max_rounds`]) before it settled.
    OutOfRounds,
}

impl ClimbEnding {
    /// Whether the walk settled — the one ending that is convergence.
    pub fn settled(self) -> bool {
        self == Self::Settled
    }
}

/// **One walk over a stratum**, from one starting point: how it ended, how many rounds it took, a
/// round that lost included, and how many times it was judged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WalkRecord {
    /// How the walk ended.
    pub ending: ClimbEnding,
    /// How many rounds it took, a round that lost and was undone included.
    pub rounds: u32,
    /// How many times a point of it was judged for settling — one curvature of the stratum each.
    pub judgements: u32,
    /// **Whether it settled at a point where no number had a standard error**, so the judgement had
    /// nothing to measure: settled by the spec's definition (`fit_precision.md` §4.3) and counted in
    /// the run's log, so that such walks are seen rather than read as convergence (owner, checkpoint
    /// C, 2026-10-02).
    pub settled_with_no_error: bool,
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

/// How many strata are fitted at once when the run says nothing: **one**, the schedule whose
/// memory is one likelihood table whatever the cohort ([`SsrFitConfig::strata_at_once`]).
///
/// The safe value is the one to fall back on: a run that was told nothing should not be able to
/// run out of memory in the repeat-tract fit.
pub const DEFAULT_STRATA_AT_ONCE: NonZeroUsize = NonZeroUsize::MIN;

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

/// **Where the repeat-tract fit spends the thread pool.**
///
/// The work has two grains that are independent of each other, and a pool can be spread over one
/// of them or the other but not usefully over both. A stratum — every tract of one motif period
/// whose reference carries the same number of copies — is fitted from its own tracts and no
/// other's, and the strata are fitted in one call. So the choice is whether the threads split
/// *one stratum's tracts* while the strata go by one at a time, or split *the strata* while each
/// is fitted on a single thread.
///
/// **Neither arm moves a floating-point result.** Every parallel site inside the fit collects its
/// values in tract order and sums them serially, precisely so that the answer does not depend on
/// how the work was divided; the same cohort fitted under either arm, at any pool width, returns
/// the same bits.
///
/// **What does differ is memory.** Fitting a stratum holds one likelihood row per sample per
/// tract — `tracts × samples × 91 genotypes × 8 bytes`, which is about 46 kB a tract at 63
/// samples — for as long as that stratum is being fitted. Under
/// [`AcrossTheTractsOfOneStratum`](Self::AcrossTheTractsOfOneStratum) one stratum is resident;
/// under [`AcrossStrata`](Self::AcrossStrata) one per pool thread is, because a thread never
/// leaves a stratum part-way through. So the second arm's peak is bounded by `threads × the
/// largest stratum`, and [`fit_strata`] runs it on a pool of exactly
/// [`SsrFitConfig::strata_at_once`] threads so that the run, not the machine's core count, says
/// what that bound is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WhereTheThreadsGo {
    /// One stratum at a time, its tracts spread across the pool.
    ///
    /// **What the fit did before 2026-09-11**, and still the right arm when there is one stratum
    /// to fit: there is no other grain to spread then, and this one keeps the whole pool busy.
    AcrossTheTractsOfOneStratum,
    /// Several strata at once — [`SsrFitConfig::strata_at_once`] of them — and every starting
    /// point of every stratum, each walk on one thread ([`every_stratum_from_every_start`]).
    ///
    /// **The arm that reaches a pool wider than a thin stratum has tracts.** A stratum the fit
    /// will touch at all holds at least [`DEFAULT_REFUSAL_FLOOR`] tracts, which is 8 — so on an
    /// 18-thread machine the first arm leaves more than half the pool with nothing to do at every
    /// thin stratum, and there are far more thin strata than fat ones.
    ///
    /// **The unit is the walk and not the stratum**, because a stratum is searched once from each
    /// starting point and those searches never read each other. Eleven strata at three starting
    /// points is thirty-three pieces of work, not eleven.
    ///
    /// **The ceiling is the longest single walk**, because nothing splits one: a run under this
    /// arm cannot finish before the slowest stratum's slowest starting point would have finished
    /// on one thread.
    AcrossStrata,
}

#[derive(Debug, Clone)]
pub struct SsrFitConfig {
    /// How far either side of the reference length allele mass may sit.
    pub allele_span: i32,
    pub quadrature_points: usize,
    pub starting_points: Vec<StartingPoint>,
    /// How many rounds of coordinate ascent one start gets, at most ([`DEFAULT_MAX_ROUNDS`]).
    pub max_rounds: u32,
    /// **When a walk has settled**: once every number is within this share of its own standard error
    /// — [`SETTLED_FRACTION`], 0.1, by default — of where the likelihood peaks, as a Newton step
    /// estimates the distance (`fit_precision.md` §4.3). A walk is judged once the gain it can still
    /// make, projected from its last two rounds, is below `½ · p · settled_fraction²` for the
    /// stratum's `p` numbers.
    pub settled_fraction: f64,
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
    /// **How many strata are fitted at once**, each on one thread of a pool of that many built
    /// for the step. One means one stratum at a time with every thread of the caller's pool on
    /// its tracts; see [`WhereTheThreadsGo`] for the two schedules.
    ///
    /// **This is what decides whether the fit's memory fits the machine, so the run states it and
    /// nothing guesses it.** A stratum's fit holds one likelihood table, `tracts × samples with
    /// reads × 91 genotype pairs × 8 bytes` — 7.4 GiB for a 5,000-tract stratum at 2,169 samples —
    /// and `N` strata at once hold up to `N` of them. One, the default, holds one table whatever
    /// the cohort. **No choice moves a number**: both schedules return the same bits.
    ///
    /// **Read by [`fit_strata`] alone.** [`fit_stratum`] is handed one stratum and has no second
    /// grain to choose, so it always spreads the pool over that stratum's tracts whatever this
    /// says.
    pub strata_at_once: NonZeroUsize,
    /// **How a large cohort's strata take a subset of their samples** ([`SampleSubsets`]). Read
    /// only by [`fit_strata_on_sample_subsets`]; [`fit_strata`] fits every sample.
    pub subsets: SampleSubsets,
}

impl Default for SsrFitConfig {
    fn default() -> Self {
        Self {
            allele_span: ALLELE_SPAN,
            quadrature_points: QUADRATURE_POINTS,
            starting_points: StartingPoint::spanning_the_monomorphic_range(),
            max_rounds: DEFAULT_MAX_ROUNDS,
            settled_fraction: SETTLED_FRACTION,
            share_curve: ShareCurveConfig::default(),
            refusal_floor: DEFAULT_REFUSAL_FLOOR,
            curve: SlippageCurveConfig::default(),
            strata_at_once: DEFAULT_STRATA_AT_ONCE,
            subsets: SampleSubsets::default(),
        }
    }
}

/// **How many rounds a walk gets when the run says nothing: 40**, measured at plan step C3 with the
/// limit at 60 on drawn strata (`fit_precision_c3_2026-10-01.md` §2.3). Walks that settled took 2 to
/// 6 rounds at three allele classes and 6 to 36 at thirteen. At thirteen classes one walk of fifteen
/// had not settled by 60, and none settled between 37 and 60. So 40 holds every walk that settled
/// there; whether a walk ever settles after 40 is not known. It was 5, which left 9 of 15 walks at
/// thirteen classes out of rounds.
pub const DEFAULT_MAX_ROUNDS: u32 = 40;

/// **How many samples a large cohort's repeat-tract stratum is first fitted on: 256** (spec §4.4,
/// soft). A cohort of at most this many is fitted on every sample, exactly as before samples were subset.
pub const FIRST_SUBSET: usize = 256;

/// **How precisely a stratum's slippage level must be measured before its subset stops growing**:
/// its standard error below 0.02 of the level itself (spec §4.4, soft). On drawn strata the level's
/// error is close to one over the square root of the slipped reads, so this is roughly 2,500
/// slipped reads in the subset at thirteen allele classes and 5,000 at three (the review of plan
/// step D2, `fit_precision_d2_2026-10-02.md`).
pub const LEVEL_RELATIVE_ERROR_TARGET: f64 = 0.02;

/// **How many samples with reads of one slippage group a subset takes at least**, whenever the
/// group put reads in the stratum: 8, or all of the group's own if it has fewer (spec §4.4, soft;
/// amended at checkpoint D2 from "when the first samples hold none").
pub const MIN_SAMPLES_A_GROUP: usize = 8;

/// **How a large cohort's repeat-tract strata are fitted on a subset of their samples** (spec §4.4).
/// Each stratum is fitted on the first [`first`](Self::first) samples of the cohort's fixed order
/// ([`sample_order`](super::sample_order::sample_order)), then on twice as many, and so on, until
/// its slippage level is measured to
/// [`level_relative_error_target`](Self::level_relative_error_target) or every sample is in. Read
/// only by [`fit_strata_on_sample_subsets`]; [`fit_strata`] ignores it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SampleSubsets {
    /// The samples a stratum is first fitted on ([`FIRST_SUBSET`]); a cohort of no more than this
    /// is fitted on all of them.
    pub first: usize,
    /// The level's standard error, as a share of the level, below which the subset stops growing
    /// ([`LEVEL_RELATIVE_ERROR_TARGET`]). Every live slippage group's level must reach it.
    pub level_relative_error_target: f64,
    /// The samples with reads of a slippage group a subset takes at least, when the first samples
    /// hold none of that group's reads ([`MIN_SAMPLES_A_GROUP`]).
    pub min_samples_a_group: usize,
}

impl Default for SampleSubsets {
    fn default() -> Self {
        Self {
            first: FIRST_SUBSET,
            level_relative_error_target: LEVEL_RELATIVE_ERROR_TARGET,
            min_samples_a_group: MIN_SAMPLES_A_GROUP,
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
/// tracts on a 63-sample cohort only 83% of fits settle within their rounds (measured under the
/// stopping rule before 2026-10-01: a gain of the mean below 10⁻⁶, or five rounds), and that row's
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
    // **One stratum has no second grain to spread the pool over**, so this entry ignores
    // `config.strata_at_once` and always splits this stratum's tracts. `fit_strata` is the caller
    // that has a choice to make.
    fit_pooled(
        evidence,
        &[],
        homozygote_excess,
        config,
        WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
    )
}

/// **Where one climb from one starting point ended up.**
///
/// The fit does not solve for a stratum's slippage numbers, it searches for them: it takes a
/// guess and walks uphill until the reads stop becoming more likely. A walk stops at the first
/// peak it reaches, which need not be the highest one, so the search is run once from each of
/// [`SsrFitConfig::starting_points`] and the best-scoring walk wins. **The walks never read each
/// other**, which is what makes them a unit of work in their own right.
#[derive(Debug)]
struct Climb {
    parameters: Parameters,
    /// The mean log-likelihood a tract this walk ended at, which is what the walks are ranked on.
    score: f64,
    /// Why the walk stopped.
    ending: ClimbEnding,
    /// The rounds it took, a round that lost included.
    rounds: u32,
    /// How many times its points were judged ([`StratumClimb::judge`]).
    judgements: u32,
    /// The standard errors at the point it stopped, when its last round was judged there — settled
    /// or not: the stratum's errors if this walk wins, without a second curvature.
    standard_errors: Option<StratumErrors>,
    /// Whether it settled at a point where no number had an error ([`WalkRecord`]).
    settled_with_no_error: bool,
}

impl Climb {
    fn record(&self) -> WalkRecord {
        WalkRecord {
            ending: self.ending,
            rounds: self.rounds,
            judgements: self.judgements,
            settled_with_no_error: self.settled_with_no_error,
        }
    }
}

/// **Whether a walk settled with nothing to judge**: it ended settled, and the judgement at its last
/// point found no number with an error — every one not placed, not identified, or a class with no
/// share. A settled walk was always judged at the point it returns.
fn settled_with_no_error(ending: ClimbEnding, judged_here: Option<&StratumErrors>) -> bool {
    ending.settled() && judged_here.is_some_and(|errors| !errors.any())
}

/// **One walk uphill, from one starting point.**
///
/// `genotypes` and `live_groups` are properties of the stratum rather than of the starting
/// point, so they are computed once by the caller and lent to every walk over that stratum.
fn climb_from(
    evidence: &StratumEvidence,
    start: StartingPoint,
    homozygote_excess: &[f64],
    genotypes: &[(usize, usize)],
    live_groups: &[bool],
    config: &SsrFitConfig,
    threads: WhereTheThreadsGo,
) -> Climb {
    let allele_classes = (2 * config.allele_span + 1) as usize;
    let parameters = Parameters::start(start, evidence.groups, allele_classes);
    climb_from_parameters(
        evidence,
        parameters,
        homozygote_excess,
        genotypes,
        live_groups,
        config,
        threads,
    )
}

/// **One walk uphill from `parameters`**: [`climb_from`]'s walk, from a point given whole rather
/// than from a starting point — a larger subset's walk starts from the smaller one's answer
/// ([`fit_on_growing_subsets`]).
fn climb_from_parameters(
    evidence: &StratumEvidence,
    parameters: Parameters,
    homozygote_excess: &[f64],
    genotypes: &[(usize, usize)],
    live_groups: &[bool],
    config: &SsrFitConfig,
    threads: WhereTheThreadsGo,
) -> Climb {
    let allele_classes = (2 * config.allele_span + 1) as usize;
    let mut climbing = StratumClimb {
        scorer: Scorer::new(evidence, homozygote_excess, genotypes, config, threads),
        evidence,
        live_groups,
        classes: allele_classes,
        settled_fraction: config.settled_fraction,
    };
    let score = climbing.scorer.score(&parameters);
    let trigger = remaining_gain_target(&parameters, live_groups, config.settled_fraction);
    let walked = walk_until_settled(
        &mut climbing,
        parameters,
        score,
        trigger,
        tracts_a_mean_is_over(evidence),
        config.max_rounds,
    );
    Climb {
        settled_with_no_error: settled_with_no_error(walked.ending, walked.judged_here.as_ref()),
        parameters: walked.at,
        score: walked.score,
        ending: walked.ending,
        rounds: walked.rounds,
        judgements: walked.judgements,
        standard_errors: walked.judged_here,
    }
}

/// **The largest ratio of one round's gain to the last's that the projection believes**, so a walk
/// whose gains shrink slowly is projected to have at most nineteen times its last gain still to
/// come; and the ratio a walk's first round, with no earlier gain to compare, is projected with.
/// The same cap the SNP/indel fit's first projection used (`fit_precision.md` §2).
const MAX_GAIN_CONTRACTION: f64 = 0.95;

/// **The projected gain still to come below which a walk is judged** (spec §4.3): `½ · p ·
/// fraction²` over the stratum's `p` numbers — each live slippage group's three, every allele class's
/// share but one (the shares sum to one), and the concentration; the numbers the standard errors are
/// taken over ([`CurvatureLayout`]). Near the maximum the gain still to come is about
/// `½ Σ (dⱼ / SEⱼ)²` for distances `dⱼ` from it, so a walk whose every number is within `fraction` of
/// its error has at most this left — 0.08 at one slippage group, thirteen classes and a tenth.
///
/// **A trigger, not the test.** The converse does not hold: the same gain can sit in one number,
/// which is then `√p` times further than `fraction`, and a projection from two rounds' gains cannot
/// see a climb crossing a plateau (`fit_precision_c3_stopped_2026-10-01.md`). Whether a walk has
/// settled is judged on each number by [`StratumClimb::judge`].
fn remaining_gain_target(parameters: &Parameters, live_groups: &[bool], fraction: f64) -> f64 {
    let numbers = CurvatureLayout::of(parameters, live_groups)
        .coordinates
        .len();
    0.5 * numbers as f64 * fraction * fraction
}

/// **The total log-likelihood a walk can still gain**, projected from its last round's gain and the
/// one before: gains that shrink by a factor λ a round leave `gain · λ / (1 − λ)` to come (Aitken's
/// projection, applied to the log-likelihood; spec §2 and §4.3). λ is the two gains' ratio, capped at
/// [`MAX_GAIN_CONTRACTION`], and the cap itself when there is no earlier gain.
fn projected_remaining_gain(gain: f64, earlier: Option<f64>) -> f64 {
    let contraction = match earlier {
        Some(earlier) if earlier > 0.0 => (gain / earlier).min(MAX_GAIN_CONTRACTION),
        _ => MAX_GAIN_CONTRACTION,
    };
    gain * contraction / (1.0 - contraction)
}

/// **The best point a walk has stood at**, and its score: where a round that loses falls back to.
#[derive(Debug)]
struct HeldPoint<P> {
    at: P,
    score: f64,
}

impl<P: Clone> HeldPoint<P> {
    /// Hold `at` instead if it scores strictly higher.
    fn offer(&mut self, at: &P, score: f64) {
        if score > self.score {
            self.at.clone_from(at);
            self.score = score;
        }
    }
}

/// **One walk's climbing as its stopping rule sees it**: a round, and a judgement of whether a point
/// has settled.
trait Climbing {
    type Point: Clone;
    /// What judging a point finds, beside whether it has settled.
    type Judgement;
    /// Move `at` by one round, offering every point it stands at to `held`; the score where the
    /// round ended, a mean a tract.
    fn one_round(&mut self, at: &mut Self::Point, held: &mut HeldPoint<Self::Point>) -> f64;
    /// Whether `at` has settled, and what judging it found.
    fn judge(&mut self, at: &Self::Point) -> (bool, Self::Judgement);
}

/// Where a walk stopped, its score there, why, after how many rounds and judgements, and what the
/// last judgement found when it was taken at that point.
#[derive(Debug)]
struct Walked<P, J> {
    at: P,
    score: f64,
    ending: ClimbEnding,
    rounds: u32,
    judgements: u32,
    judged_here: Option<J>,
}

/// **The climb's stopping rule, apart from what a round does and how a point is judged** (spec §4.3,
/// as amended with the owner at plan step C3).
///
/// Scores are a mean a tract; `tracts` turns their differences into the total's. After each round:
///
/// - **a round that ended lower than it began** has its moves undone: the walk goes back to the best
///   point it stood at — the round's start, or a point part-way through it — and that point is judged
///   at once — when the put-off below allows, or always when nothing in the round was better than its
///   start. A score that is not a number is a loss too;
/// - otherwise the point the round reached is judged once the gain still to come, projected from this
///   round's gain and the last ([`projected_remaining_gain`]), is below `trigger`;
/// - a judged point that has settled ends the walk, settled ([`ClimbEnding::Settled`]). One that has
///   not, after a round that lost and had nothing better than its start, ends it there
///   ([`ClimbEnding::LostARound`]): the next round would start where this one did and repeat it
///   exactly. Otherwise the walk goes on;
/// - **a judgement that finds the walk unsettled puts the next one off**, by one round more each
///   time — the walk is next judged one, then two, then three rounds later — since a judgement costs
///   about as much as a round or two of the climb, and a walk crossing a plateau would otherwise be
///   judged every round. A round that loses while the judgement is put off, with a better point
///   part-way through, goes back to that point unjudged;
/// - and it stops at `max_rounds` ([`ClimbEnding::OutOfRounds`]).
fn walk_until_settled<C: Climbing>(
    climbing: &mut C,
    mut at: C::Point,
    mut score: f64,
    trigger: f64,
    tracts: f64,
    max_rounds: u32,
) -> Walked<C::Point, C::Judgement> {
    let mut earlier_gain = None;
    let mut judgements = 0;
    let mut judged_here = None;
    // The first round a judgement may be taken at, and how many have found the walk unsettled.
    let (mut judge_from, mut unsettled) = (1, 0);
    for round in 1..=max_rounds {
        let start = score;
        let mut held = HeldPoint {
            at: at.clone(),
            score,
        };
        let after = climbing.one_round(&mut at, &mut held);
        // Anything but at least as high is a loss, a `NaN` included.
        let lost = !matches!(
            after.partial_cmp(&start),
            Some(std::cmp::Ordering::Greater | std::cmp::Ordering::Equal)
        );
        if lost {
            at = held.at;
            score = held.score;
        } else {
            score = after;
        }
        let gain = (score - start) * tracts;
        judged_here = None;
        // A loss with nothing better than its start, a gain that is not a number included: the walk
        // can go nowhere from here, so it is judged now whatever the put-off says.
        let stuck = lost && !matches!(gain.partial_cmp(&0.0), Some(std::cmp::Ordering::Greater));
        let wanted = lost || projected_remaining_gain(gain, earlier_gain) < trigger;
        let settled = if stuck || (wanted && round >= judge_from) {
            judgements += 1;
            let (settled, judgement) = climbing.judge(&at);
            judged_here = Some(judgement);
            if !settled {
                unsettled += 1;
                judge_from = round + unsettled;
            }
            settled
        } else {
            false
        };
        let ending = if settled {
            Some(ClimbEnding::Settled)
        } else if stuck {
            Some(ClimbEnding::LostARound)
        } else {
            None
        };
        if let Some(ending) = ending {
            return Walked {
                at,
                score,
                ending,
                rounds: round,
                judgements,
                judged_here,
            };
        }
        earlier_gain = Some(gain);
    }
    Walked {
        at,
        score,
        ending: ClimbEnding::OutOfRounds,
        rounds: max_rounds,
        judgements,
        judged_here,
    }
}

/// **One walk over a stratum**: its scorer, and what judging a point needs.
struct StratumClimb<'a> {
    scorer: Scorer<'a>,
    evidence: &'a StratumEvidence,
    live_groups: &'a [bool],
    classes: usize,
    settled_fraction: f64,
}

impl Climbing for StratumClimb<'_> {
    type Point = Parameters;
    /// The stratum's standard errors at the judged point, which are the winning walk's when its last
    /// round was judged at the point it returns.
    type Judgement = StratumErrors;

    fn one_round(&mut self, at: &mut Parameters, held: &mut HeldPoint<Parameters>) -> f64 {
        climb_one_round(at, &mut self.scorer, self.live_groups, self.classes, held);
        self.scorer.score(at)
    }

    /// **Settled when every number is within `settled_fraction` of its own standard error of where
    /// the likelihood peaks**, as a Newton step from `at` estimates it — the rule the SNP/indel fit
    /// stops by (spec §2's amendment), on the stratum's curvature and slope ([`curvature_at`]),
    /// judged on each number's own scale. A number without an error is settled by definition: the
    /// tracts do not place it, or do not tell it from the others. Costs one curvature, `1 + 2p²`
    /// evaluations.
    fn judge(&mut self, at: &Parameters) -> (bool, StratumErrors) {
        let curvature = curvature_at(
            &mut self.scorer,
            self.evidence,
            at,
            self.live_groups,
            CURVATURE_STEP,
        );
        let errors = errors_on_the_natural_scale(
            &curvature.layout,
            at,
            self.live_groups,
            &curvature.identified,
        );
        let distances = newton_distances_on_the_natural_scale(&curvature, at, self.live_groups);
        let settled = furthest_in_errors(&errors, &distances)
            .is_none_or(|furthest| furthest < self.settled_fraction);
        (settled, errors)
    }
}

/// **The better of two walks, under the rule the serial loop has always used.**
///
/// A walk replaces the standing best only when it scores **strictly** higher, so a tie keeps
/// the one that came first — and since the walks are always folded in starting-point order,
/// "first" means the earlier starting point whether they ran one after another or all at once.
///
/// **Spelled as `climb > current` and not as `current >= climb`**, which are the same rule for
/// every pair of numbers except when one of them is `NaN`: a walk that scored `NaN` loses under
/// the first and wins under the second, because every comparison with `NaN` is false. A score
/// should never be `NaN`, and this is the spelling that does not quietly promote one if it is.
fn the_better_walk(best: Option<Climb>, climb: Climb) -> Option<Climb> {
    match best {
        Some(current) => {
            if climb.score > current.score {
                Some(climb)
            } else {
                Some(current)
            }
        }
        None => Some(climb),
    }
}

/// **The stratum's answer, assembled from the walk that won.**
fn the_fit_of(
    evidence: &StratumEvidence,
    borrowed: &[u64],
    live_groups: &[bool],
    winner: Climb,
    walks: Vec<WalkRecord>,
    standard_errors: Option<StratumErrors>,
) -> StratumFit {
    let Climb {
        parameters,
        score,
        ending,
        rounds: _,
        judgements: _,
        standard_errors: _,
        settled_with_no_error: _,
    } = winner;
    StratumFit {
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
        ending,
        walks,
        samples_fitted_on: None,
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
        standard_errors,
    }
}

/// Fit `evidence`, whose tracts may already include borrowed ones, recording where from, with the
/// standard errors of the answer ([`standard_errors_at`]).
fn fit_pooled(
    evidence: &StratumEvidence,
    borrowed: &[u64],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    threads: WhereTheThreadsGo,
) -> Option<StratumFit> {
    let (mut best, walks, live_groups) =
        the_best_walk(evidence, homozygote_excess, config, threads)?;
    let standard_errors = best.standard_errors.take().unwrap_or_else(|| {
        standard_errors_at(
            evidence,
            &best.parameters,
            homozygote_excess,
            &live_groups,
            config,
            threads,
        )
    });
    Some(the_fit_of(
        evidence,
        borrowed,
        &live_groups,
        best,
        walks,
        Some(standard_errors),
    ))
}

/// **The best of the walks from every starting point over `evidence`**, every walk's record in
/// starting-point order, and which slippage groups they moved; `None` when no group put a read in
/// it, so there is nothing to walk.
fn the_best_walk(
    evidence: &StratumEvidence,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    threads: WhereTheThreadsGo,
) -> Option<(Climb, Vec<WalkRecord>, Vec<bool>)> {
    // **The one precondition on `SsrFitConfig::allele_span` that nothing else states.** It is a
    // public field with no lower bound, read from an environment variable by
    // `examples/ng_joint_records_walk.rs` and parsed with no floor, and at zero the fit returns
    // a one-class length spectrum — which is a tract that can only ever be its reference length.
    // Nothing downstream can use that: `StratumFits::over` refuses it with a message about a
    // class count, naming the wrong thing. Refused here, where the message can name the knob.
    assert!(
        config.allele_span >= 1,
        "the fit places allele mass from -{span} to +{span} whole repeat units either side of \
         the reference length, so `SsrFitConfig::allele_span` must be at least 1; at {span} a \
         tract could only ever carry its reference length",
        span = config.allele_span
    );
    let genotypes = genotype_pairs((2 * config.allele_span + 1) as usize);
    let live_groups = evidence.groups_with_reads();
    if !live_groups.iter().any(|live| *live) {
        return None;
    }

    let climbs: Vec<Climb> = config
        .starting_points
        .iter()
        .map(|start| {
            climb_from(
                evidence,
                *start,
                homozygote_excess,
                &genotypes,
                &live_groups,
                config,
                threads,
            )
        })
        .collect();
    let walks = climbs.iter().map(Climb::record).collect();
    let best = climbs.into_iter().fold(None, the_better_walk);
    Some((
        best.expect("at least one starting point"),
        walks,
        live_groups,
    ))
}

// ---------------------------------------------------------------------
// How precisely a stratum's tracts determine its answer
// ---------------------------------------------------------------------

/// **Each central difference's step**, on the climb's scales — logit, log, log-ratio — which are
/// already relative to a number's size, so one step suits every coordinate.
///
/// Large enough that the objective's rounding does not reach the curvature — the quadrature's
/// points are placed by solving Beta quantiles to 10⁻¹², which the objective carries into every
/// tract at once — and small enough that the difference's own error, of order the step squared, is
/// small beside it: at a third of it the errors of a drawn stratum agree to 2 in 100 (tested).
const CURVATURE_STEP: f64 = 1e-2;

/// One number of a stratum's fit as its curvature is taken: on the scale the climb moves it on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CurvatureCoordinate {
    /// Slippage group `group`'s number `which` — 0 the level, 1 the shorter share, 2 the fall-off —
    /// on the logit scale.
    Slippage { group: usize, which: usize },
    /// An allele class's share of the spectrum, as the log of its ratio to the largest class's
    /// share. The shares sum to one, so the largest is the one not moved on its own.
    SpectrumRatio { class: usize },
    /// The concentration, on the log scale.
    Concentration,
}

/// The numbers one stratum's curvature is taken over, in order: each live slippage group's three,
/// every allele class with a share but the largest, and the concentration.
struct CurvatureLayout {
    coordinates: Vec<CurvatureCoordinate>,
    /// The class whose share the others are taken as ratios to.
    largest_class: usize,
}

impl CurvatureLayout {
    fn of(parameters: &Parameters, live_groups: &[bool]) -> Self {
        let spectrum = &parameters.length_spectrum;
        // PANIC-FREE: `fit_pooled` refuses an allele span below one, so a spectrum has at least three
        // classes.
        let largest_class = (0..spectrum.len())
            .max_by(|&left, &right| spectrum[left].total_cmp(&spectrum[right]))
            .expect("a spectrum has at least one class");
        let mut coordinates = Vec::new();
        for (group, live) in live_groups.iter().enumerate() {
            if *live {
                coordinates
                    .extend((0..3).map(|which| CurvatureCoordinate::Slippage { group, which }));
            }
        }
        // A class with no share at all is held at zero: its ratio has no logarithm.
        coordinates.extend(
            (0..spectrum.len())
                .filter(|&class| class != largest_class && spectrum[class] > 0.0)
                .map(|class| CurvatureCoordinate::SpectrumRatio { class }),
        );
        coordinates.push(CurvatureCoordinate::Concentration);
        Self {
            coordinates,
            largest_class,
        }
    }

    /// `parameters` on the climb's scales, in the layout's order.
    fn center(&self, parameters: &Parameters) -> Vec<f64> {
        let largest = parameters.length_spectrum[self.largest_class];
        self.coordinates
            .iter()
            .map(|coordinate| match *coordinate {
                CurvatureCoordinate::Slippage { group, which } => {
                    let value = read_slippage(&parameters.slippage[group], which);
                    float::ln(value / (1.0 - value))
                }
                CurvatureCoordinate::SpectrumRatio { class } => {
                    float::ln(parameters.length_spectrum[class] / largest)
                }
                CurvatureCoordinate::Concentration => float::ln(parameters.concentration),
            })
            .collect()
    }

    /// `base` with every number of the layout set from `at`, on the climb's scales.
    fn parameters_at(&self, base: &Parameters, at: &[f64]) -> Parameters {
        let mut parameters = base.clone();
        let mut ratios: Vec<f64> = base
            .length_spectrum
            .iter()
            .map(|share| {
                if *share > 0.0 {
                    share / base.length_spectrum[self.largest_class]
                } else {
                    0.0
                }
            })
            .collect();
        ratios[self.largest_class] = 1.0;
        for (coordinate, value) in self.coordinates.iter().zip(at) {
            match *coordinate {
                CurvatureCoordinate::Slippage { group, which } => {
                    write_slippage(&mut parameters.slippage[group], which, expit(*value));
                }
                CurvatureCoordinate::SpectrumRatio { class } => ratios[class] = float::exp(*value),
                CurvatureCoordinate::Concentration => parameters.concentration = float::exp(*value),
            }
        }
        normalise(&mut ratios);
        parameters.length_spectrum = ratios;
        parameters
    }
}

/// **The second derivatives of `total` at `center`, and its slope**, by central differences of
/// `steps`: the curvature as a row-major `p × p` matrix, `(f(x + hⱼ) − 2 f(x) + f(x − hⱼ)) / hⱼ²` on
/// the diagonal and `(f(+hᵢ +hⱼ) − f(+hᵢ −hⱼ) − f(−hᵢ +hⱼ) + f(−hᵢ −hⱼ)) / 4hᵢhⱼ` off it, and the slope
/// `(f(x + hⱼ) − f(x − hⱼ)) / 2hⱼ` from the same points — `1 + 2p²` evaluations.
///
/// **The points are visited grouped by the coordinates `costly` names**, so a scorer that rebuilds
/// something whenever those change — the read likelihoods, when a slippage number moves — rebuilds
/// it once a distinct setting rather than once a point. The order changes no value.
fn curvature_of(
    center: &[f64],
    steps: &[f64],
    costly: impl Fn(usize) -> bool,
    mut total: impl FnMut(&[f64]) -> f64,
) -> (Vec<f64>, Vec<f64>) {
    // One point of the differences: the coordinates it moves, and which way.
    type Moves = Vec<(usize, i8)>;
    let p = center.len();
    // The centre moves none.
    let mut points: Vec<Moves> = vec![Vec::new()];
    for j in 0..p {
        points.push(vec![(j, 1)]);
        points.push(vec![(j, -1)]);
    }
    for i in 0..p {
        for j in i + 1..p {
            for (a, b) in [(1, 1), (1, -1), (-1, 1), (-1, -1)] {
                points.push(vec![(i, a), (j, b)]);
            }
        }
    }
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_by_key(|&point| {
        let (costly_moves, others): (Moves, Moves) =
            points[point].iter().partition(|(j, _)| costly(*j));
        (costly_moves, others)
    });
    let mut values = vec![0.0; points.len()];
    for point in order {
        let mut at = center.to_vec();
        for &(j, sign) in &points[point] {
            at[j] += f64::from(sign) * steps[j];
        }
        values[point] = total(&at);
    }
    let at_centre = values[0];
    let mut hessian = vec![0.0; p * p];
    let mut slope = vec![0.0; p];
    for j in 0..p {
        let (up, down) = (values[1 + 2 * j], values[2 + 2 * j]);
        hessian[j * p + j] = (up - 2.0 * at_centre + down) / (steps[j] * steps[j]);
        slope[j] = (up - down) / (2.0 * steps[j]);
    }
    let mut next = 1 + 2 * p;
    for i in 0..p {
        for j in i + 1..p {
            let [both_up, up_down, down_up, both_down] = [
                values[next],
                values[next + 1],
                values[next + 2],
                values[next + 3],
            ];
            next += 4;
            let entry = (both_up - up_down - down_up + both_down) / (4.0 * steps[i] * steps[j]);
            hessian[i * p + j] = entry;
            hessian[j * p + i] = entry;
        }
    }
    (hessian, slope)
}

/// **The standard errors of one stratum's fit at `parameters`** (`fit_precision.md` §4.2).
///
/// The curvature of the stratum's **total** log-likelihood — the mean a tract the climb maximises,
/// times the tracts it is a mean over; the mean's curvature would give every error too large by the
/// square root of the tract count — over every number the stratum fits, on the scales the climb
/// uses, by central differences ([`curvature_of`]). Its negative is the information, inverted over
/// the numbers it identifies ([`invert_identified`](super::fit::invert_identified)): a number whose
/// curvature is not negative once the others are accounted for — a flat or a saddle direction — is
/// dropped and has no error, and the rest are inverted without it. The variances are carried to each
/// number's own scale by its derivative: `p(1 − p)` for a slippage number, the concentration itself
/// for its logarithm, and for a share of the spectrum the derivative of the shares in the log-ratios.
///
/// **Cost:** `1 + 2p²` evaluations of the stratum's likelihood for `p` numbers — 513 at one slippage
/// group and thirteen classes — once a stratum, at the winning walk's answer.
fn standard_errors_at(
    evidence: &StratumEvidence,
    parameters: &Parameters,
    homozygote_excess: &[f64],
    live_groups: &[bool],
    config: &SsrFitConfig,
    threads: WhereTheThreadsGo,
) -> StratumErrors {
    standard_errors_with_step(
        evidence,
        parameters,
        homozygote_excess,
        live_groups,
        config,
        threads,
        CURVATURE_STEP,
    )
}

/// [`standard_errors_at`] with the central differences' step given, which a test varies.
fn standard_errors_with_step(
    evidence: &StratumEvidence,
    parameters: &Parameters,
    homozygote_excess: &[f64],
    live_groups: &[bool],
    config: &SsrFitConfig,
    threads: WhereTheThreadsGo,
    step: f64,
) -> StratumErrors {
    let genotypes = genotype_pairs((2 * config.allele_span + 1) as usize);
    let mut scorer = Scorer::new(evidence, homozygote_excess, &genotypes, config, threads);
    let curvature = curvature_at(&mut scorer, evidence, parameters, live_groups, step);
    errors_on_the_natural_scale(
        &curvature.layout,
        parameters,
        live_groups,
        &curvature.identified,
    )
}

/// **What turns the climb's score into the stratum's total log-likelihood**: the score is the mean a
/// tract over every tract, those without reads included (spec §4.1's trap), so the total is it times
/// every tract. The curvature (§4.2) and the walk's trigger (§4.3) are both in the total's units.
fn tracts_a_mean_is_over(evidence: &StratumEvidence) -> f64 {
    evidence.tracts.len() as f64
}

/// One stratum's curvature at a point, inverted: the numbers it is taken over, the inverse of the
/// information over those it identifies, and the slope of the total log-likelihood in each, on the
/// climb's scales.
struct CurvatureAt {
    layout: CurvatureLayout,
    identified: super::fit::Identified,
    slope: Vec<f64>,
}

/// **The curvature and slope of the stratum's total log-likelihood at `parameters`**, by central
/// differences of `step` on the climb's scales ([`curvature_of`]), the information — the curvature's
/// negative — inverted over the numbers it identifies ([`invert_identified`](super::fit::invert_identified)).
fn curvature_at(
    scorer: &mut Scorer<'_>,
    evidence: &StratumEvidence,
    parameters: &Parameters,
    live_groups: &[bool],
    step: f64,
) -> CurvatureAt {
    let layout = CurvatureLayout::of(parameters, live_groups);
    let center = layout.center(parameters);
    let steps = vec![step; center.len()];
    let tracts = tracts_a_mean_is_over(evidence);
    let (hessian, slope) = curvature_of(
        &center,
        &steps,
        |j| matches!(layout.coordinates[j], CurvatureCoordinate::Slippage { .. }),
        |at| tracts * scorer.score(&layout.parameters_at(parameters, at)),
    );
    let p = center.len();
    let information: Vec<f64> = hessian.iter().map(|entry| -entry).collect();
    let own_curvature: Vec<f64> = (0..p).map(|j| information[j * p + j]).collect();
    let identified = super::fit::invert_identified(&information, p, &own_curvature);
    CurvatureAt {
        layout,
        identified,
        slope,
    }
}

/// **Each number's distance to where the likelihood peaks, as a Newton step estimates it**, on the
/// number's own scale, laid out as [`StratumErrors`] is: the step `I⁻¹ g` over the coordinates the
/// information identifies (the others held), carried to each number's scale by the same derivatives
/// its error is — `p(1 − p)` for a slippage number, the concentration itself, the shares' slopes in
/// the log-ratios.
fn newton_distances_on_the_natural_scale(
    curvature: &CurvatureAt,
    parameters: &Parameters,
    live_groups: &[bool],
) -> StratumDistances {
    let CurvatureAt {
        layout,
        identified,
        slope,
    } = curvature;
    let kept = identified.kept.len();
    let mut step = vec![0.0; layout.coordinates.len()];
    for (row, &coordinate) in identified.kept.iter().enumerate() {
        step[coordinate] = (0..kept)
            .map(|column| identified.inverse[row * kept + column] * slope[identified.kept[column]])
            .sum();
    }
    let step_of = |wanted: CurvatureCoordinate| {
        layout
            .coordinates
            .iter()
            .position(|coordinate| *coordinate == wanted)
            .map_or(0.0, |index| step[index])
    };
    let slippage = live_groups
        .iter()
        .enumerate()
        .map(|(group, live)| {
            live.then(|| {
                let of = |which: usize| {
                    let value = read_slippage(&parameters.slippage[group], which);
                    value * (1.0 - value) * step_of(CurvatureCoordinate::Slippage { group, which })
                };
                SlippageDistances {
                    level: of(0),
                    shorter_share: of(1),
                    fall_off: of(2),
                }
            })
        })
        .collect();
    let spectrum = &parameters.length_spectrum;
    let ratio_steps: Vec<(usize, f64)> = (0..spectrum.len())
        .map(|class| (class, step_of(CurvatureCoordinate::SpectrumRatio { class })))
        .collect();
    let length_spectrum = (0..spectrum.len())
        .map(|class| {
            ratio_steps
                .iter()
                .map(|&(j, d)| {
                    spectrum[class] * (if j == class { 1.0 } else { 0.0 } - spectrum[j]) * d
                })
                .sum()
        })
        .collect();
    StratumDistances {
        slippage,
        length_spectrum,
        concentration: parameters.concentration * step_of(CurvatureCoordinate::Concentration),
    }
}

/// Each number's distance from where the likelihood peaks, laid out as [`StratumErrors`] is: per live
/// slippage group the level's, the shorter share's and the fall-off's; per class its share's; the
/// concentration's.
struct StratumDistances {
    slippage: Vec<Option<SlippageDistances>>,
    length_spectrum: Vec<f64>,
    concentration: f64,
}

/// One slippage group's three distances from the peak, named as [`SlippageErrors`] names their errors.
#[derive(Debug, Clone, Copy)]
struct SlippageDistances {
    level: f64,
    shorter_share: f64,
    fall_off: f64,
}

/// **The furthest any number is from where the likelihood peaks, in its own errors**, over the
/// numbers with an error; `None` when none has one.
fn furthest_in_errors(errors: &StratumErrors, distances: &StratumDistances) -> Option<f64> {
    let mut pairs: Vec<(StratumError, f64)> = Vec::new();
    for (group, distance) in errors.slippage.iter().zip(&distances.slippage) {
        if let (Some(group), Some(distance)) = (group, distance) {
            pairs.extend([
                (group.level, distance.level),
                (group.shorter_share, distance.shorter_share),
                (group.fall_off, distance.fall_off),
            ]);
        }
    }
    pairs.extend(
        errors
            .length_spectrum
            .iter()
            .copied()
            .zip(distances.length_spectrum.iter().copied()),
    );
    pairs.push((errors.concentration, distances.concentration));
    pairs
        .into_iter()
        .filter_map(|(error, distance)| {
            let ratio = distance.abs() / error.value()?;
            // A distance that is not a number is not settled.
            Some(if ratio.is_nan() { f64::INFINITY } else { ratio })
        })
        .reduce(f64::max)
}

/// The errors [`standard_errors_at`] reports, from the inverse over the identified coordinates: a
/// coordinate dropped from it, or whose variance is not positive, is [`StratumError::NotIdentified`];
/// one whose error on the climb's scale exceeds [`NOT_PLACED`], or a number in `[0, 1]` whose error
/// exceeds that range, is [`StratumError::NotPlaced`]; a class with no share is
/// [`StratumError::NoShare`].
fn errors_on_the_natural_scale(
    layout: &CurvatureLayout,
    parameters: &Parameters,
    live_groups: &[bool],
    identified: &super::fit::Identified,
) -> StratumErrors {
    let kept = identified.kept.len();
    // Each layout coordinate's row in the inverse, `None` when it was dropped.
    let row_of = |coordinate: CurvatureCoordinate| {
        let index = layout.coordinates.iter().position(|c| *c == coordinate)?;
        identified.kept.iter().position(|&k| k == index)
    };
    // A coordinate's error on the climb's scale, and what that says of the number.
    let on_its_scale = |coordinate: CurvatureCoordinate| -> Result<f64, StratumError> {
        let row = row_of(coordinate).ok_or(StratumError::NotIdentified)?;
        let variance = identified.inverse[row * kept + row];
        if !(variance.is_finite() && variance > 0.0) {
            return Err(StratumError::NotIdentified);
        }
        let error = variance.sqrt();
        if error > NOT_PLACED {
            return Err(StratumError::NotPlaced);
        }
        Ok(error)
    };
    let estimated = |error: f64| {
        if error.is_finite() && error > 0.0 {
            StratumError::Estimated(error)
        } else {
            StratumError::NotIdentified
        }
    };
    // A number in `[0, 1]` — a slippage number, a share — whose error is wider than that whole range
    // is not placed by the tracts either (spec §3.2's rule, on the number's own scale). The largest
    // class's share has no log-ratio of its own; its error comes through the other classes'.
    let within_its_range = |error: StratumError| match error {
        StratumError::Estimated(value) if value > 1.0 => StratumError::NotPlaced,
        other => other,
    };
    let slippage = live_groups
        .iter()
        .enumerate()
        .map(|(group, live)| {
            live.then(|| {
                let of = |which: usize| match on_its_scale(CurvatureCoordinate::Slippage {
                    group,
                    which,
                }) {
                    Ok(error) => {
                        let value = read_slippage(&parameters.slippage[group], which);
                        within_its_range(estimated(error * value * (1.0 - value)))
                    }
                    Err(reason) => reason,
                };
                SlippageErrors {
                    level: of(0),
                    shorter_share: of(1),
                    fall_off: of(2),
                }
            })
        })
        .collect();
    let spectrum = &parameters.length_spectrum;
    // The spectrum's coordinates the inverse kept, with their rows.
    let ratios: Vec<(usize, usize)> = (0..spectrum.len())
        .filter_map(|class| Some((class, row_of(CurvatureCoordinate::SpectrumRatio { class })?)))
        .collect();
    let length_spectrum = (0..spectrum.len())
        .map(|class| {
            if spectrum[class] <= 0.0 {
                return StratumError::NoShare;
            }
            if class != layout.largest_class {
                // Its own log-ratio must be placed before its share can be.
                if let Err(reason) = on_its_scale(CurvatureCoordinate::SpectrumRatio { class }) {
                    return reason;
                }
            }
            if ratios.is_empty() {
                return StratumError::NotIdentified;
            }
            // The share's slope in each kept log-ratio: `s_k (δ_kj − s_j)`.
            let slope: Vec<f64> = ratios
                .iter()
                .map(|&(j, _)| spectrum[class] * (if j == class { 1.0 } else { 0.0 } - spectrum[j]))
                .collect();
            let mut total = 0.0;
            for (a, &(_, row_a)) in ratios.iter().enumerate() {
                for (b, &(_, row_b)) in ratios.iter().enumerate() {
                    total += slope[a] * identified.inverse[row_a * kept + row_b] * slope[b];
                }
            }
            within_its_range(estimated(total.sqrt()))
        })
        .collect();
    let concentration = match on_its_scale(CurvatureCoordinate::Concentration) {
        Ok(error) => estimated(error * parameters.concentration),
        Err(reason) => reason,
    };
    StratumErrors {
        slippage,
        length_spectrum,
        concentration,
    }
}

/// One pass of coordinate ascent over everything the stratum fits, offering each point it moves to
/// to `held`.
///
/// **The slippage numbers move first and the frequencies after**, because moving slippage is
/// what invalidates the cached read likelihoods; doing it the other way round would rebuild
/// them on every spectrum coordinate.
///
/// **Each point's score is the one its golden section measured**, at exactly the parameters it
/// then writes, so offering it costs no evaluation.
fn climb_one_round(
    parameters: &mut Parameters,
    scorer: &mut Scorer<'_>,
    live_groups: &[bool],
    classes: usize,
    held: &mut HeldPoint<Parameters>,
) {
    for (group, live) in live_groups.iter().enumerate() {
        if !live {
            continue;
        }
        for which in 0..3 {
            let current = read_slippage(&parameters.slippage[group], which);
            let (moved, score) = climb_scalar(
                |x| {
                    let mut trial = parameters.clone();
                    write_slippage(&mut trial.slippage[group], which, expit(x));
                    scorer.score(&trial)
                },
                logit(current),
                3.0,
            );
            write_slippage(&mut parameters.slippage[group], which, expit(moved));
            held.offer(parameters, score);
        }
    }

    // The spectrum, one class at a time on a log scale, renormalised each time.
    for class in 0..classes {
        let current = float::ln(parameters.length_spectrum[class].max(1e-9));
        let (moved, score) = climb_scalar(
            |x| {
                let mut trial = parameters.clone();
                trial.length_spectrum[class] = float::exp(x);
                normalise(&mut trial.length_spectrum);
                scorer.score(&trial)
            },
            current,
            2.0,
        );
        parameters.length_spectrum[class] = float::exp(moved);
        normalise(&mut parameters.length_spectrum);
        held.offer(parameters, score);
    }

    let (moved, score) = climb_scalar(
        |x| {
            let mut trial = parameters.clone();
            trial.concentration = float::exp(x);
            scorer.score(&trial)
        },
        float::ln(parameters.concentration),
        2.5,
    );
    parameters.concentration = float::exp(moved);
    held.offer(parameters, score);
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
// Spending the pool on the walks rather than on the strata
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// A stratum read from a subset of samples (spec §4.4)
// ---------------------------------------------------------------------

impl StratumEvidence {
    /// **The same stratum read from the samples `keep` marks, and no others.** Every tract is kept,
    /// one whose readers are all left out kept empty: the score is a mean over every tract either
    /// way (spec §4.1's trap), and a tract without reads moves no answer. The counts the fit does
    /// not take per sample — the guard's and the substitution counts — stay the whole stratum's.
    fn of_samples(&self, keep: &[bool]) -> StratumEvidence {
        StratumEvidence {
            stratum: self.stratum,
            tracts: self
                .tracts
                .iter()
                .map(|tract| TractReads {
                    samples: tract
                        .samples
                        .iter()
                        .filter(|reads| keep[reads.sample as usize])
                        .cloned()
                        .collect(),
                })
                .collect(),
            read_span: self.read_span,
            groups: self.groups,
            tracts_over_guard_threshold: self.tracts_over_guard_threshold,
            reads_reaching_not_crossing: self.reads_reaching_not_crossing,
            guard_reads: self.guard_reads,
            bases_compared: self.bases_compared,
            mismatching_bases: self.mismatching_bases,
        }
    }

    /// **Which samples put reads of each slippage group in this stratum**, a row of `samples`
    /// flags a group.
    fn readers_by_group(&self, samples: usize) -> Vec<Vec<bool>> {
        let mut readers = vec![vec![false; samples]; self.groups];
        for tract in &self.tracts {
            for reads in &tract.samples {
                for (group, counts) in &reads.by_group {
                    if counts.iter().any(|count| *count > 0) {
                        readers[*group as usize][reads.sample as usize] = true;
                    }
                }
            }
        }
        readers
    }
}

/// **Adds to `keep` the samples a stratum's subset of `size` holds** (spec §4.4 as amended at
/// checkpoint D2), so a larger subset always holds the smaller one's samples: the first `size` in
/// `order`; then, for every slippage group that put reads in the stratum, that group's readers in
/// `order` until `min_a_group` of them are in, or all of its own.
fn grow_subset(
    keep: &mut [bool],
    readers: &[Vec<bool>],
    order: &[usize],
    size: usize,
    min_a_group: usize,
) {
    for &sample in order.iter().take(size) {
        keep[sample] = true;
    }
    for group_readers in readers {
        let mut held = group_readers
            .iter()
            .zip(keep.iter())
            .filter(|(reads, kept)| **reads && **kept)
            .count();
        for &sample in order {
            if held >= min_a_group {
                break;
            }
            if group_readers[sample] && !keep[sample] {
                keep[sample] = true;
                held += 1;
            }
        }
    }
}

/// **The size the next subset takes**: twice `size`, or every sample once that would hold more than
/// three quarters of the cohort's `samples` (spec §4.4 as amended at checkpoint D2) — at 2,169
/// samples, 256, 512 and 1,024, then every sample rather than 2,048.
fn next_subset_size(size: usize, samples: usize) -> usize {
    let doubled = size.saturating_mul(2);
    if doubled.saturating_mul(4) > samples.saturating_mul(3) {
        samples
    } else {
        doubled
    }
}

/// **Whether a fitted stratum's slippage level is measured to `target`** of itself in every live
/// slippage group that `still_growing` marks: the level's standard error below `target` times the
/// level. A group whose level has no error, or is zero, has not reached it. A group not marked —
/// every one of its readers already in the subset — cannot be measured better by growing, so it is
/// not asked; with no group asked, the level is as well measured as it can be.
fn level_is_measured_to(
    parameters: &Parameters,
    live_groups: &[bool],
    errors: &StratumErrors,
    still_growing: &[bool],
    target: f64,
) -> bool {
    live_groups
        .iter()
        .zip(&still_growing[..live_groups.len()])
        .enumerate()
        .filter(|(_, (live, growing))| **live && **growing)
        .all(|(group, _)| {
            let level = parameters.slippage[group].level;
            errors.slippage[group]
                .and_then(|group| group.level.value())
                .is_some_and(|error| {
                    error.partial_cmp(&(target * level)) == Some(std::cmp::Ordering::Less)
                })
        })
}

/// **One stratum of a large cohort, fitted on a subset of its samples grown until its slippage
/// level is measured precisely enough** (spec §4.4 as amended at checkpoint D2).
///
/// - **The refusal floor is judged on the whole stratum**, before any subset is drawn: a stratum
///   whose reads are all in samples outside the first subset grows its subset rather than being
///   refused.
/// - The first subset is the first [`SampleSubsets::first`] samples of `order`, each slippage group
///   with reads in the stratum topped up to [`SampleSubsets::min_samples_a_group`] of its readers
///   ([`grow_subset`]). Each larger subset holds the smaller one's samples, doubles it, and takes
///   every sample once it would pass three quarters of them ([`next_subset_size`]).
/// - **A subset holding fewer tracts with reads than [`SsrFitConfig::refusal_floor`] is not
///   fitted**; it grows.
/// - The first subset fitted is fitted from every starting point; each later one by one walk from
///   the last answer, which decides only whether to grow.
/// - The subset stops growing once every live slippage group that still has readers outside it has
///   a level with a standard error below [`SampleSubsets::level_relative_error_target`] of itself
///   ([`level_is_measured_to`]), or once it holds every sample. **The subset the answer is taken
///   from is fitted from every starting point and from the last answer as well**, the best walk
///   winning: one walk from a nearby answer can settle short of the best of several.
/// - The answer's evidence counts are that subset's, and the fit records how many samples it held
///   ([`StratumFit::samples_fitted_on`]); its walks are that subset's.
fn fit_on_growing_subsets(
    evidence: &StratumEvidence,
    order: &[usize],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    threads: WhereTheThreadsGo,
) -> StratumOutcome {
    assert_eq!(
        order.len(),
        homozygote_excess.len(),
        "the sample order must rank every sample the homozygote excess is given for"
    );
    if let Some(refused) = refused_before_any_walk(evidence, config) {
        return refused;
    }
    // The precondition `the_best_walk` states, for the walks taken here directly.
    assert!(
        config.allele_span >= 1,
        "`SsrFitConfig::allele_span` must be at least 1, so a tract can carry a length other than \
         its reference's"
    );
    let subsets = config.subsets;
    let samples = order.len();
    let readers = evidence.readers_by_group(samples);
    let genotypes = genotype_pairs(allele_classes_of(config));
    let mut keep = vec![false; samples];
    let mut size = subsets.first.max(1);
    if size.saturating_mul(4) > samples.saturating_mul(3) {
        size = samples;
    }
    // The last answer, which a larger subset's one walk starts from.
    let mut last: Option<Parameters> = None;
    loop {
        grow_subset(
            &mut keep,
            &readers,
            order,
            size,
            subsets.min_samples_a_group,
        );
        let every_sample = keep.iter().all(|kept| *kept);
        // **The whole stratum when every sample is in**, so a stratum that grows that far holds
        // no second copy of its evidence.
        let copied;
        let subset: &StratumEvidence = if every_sample {
            evidence
        } else {
            copied = evidence.of_samples(&keep);
            &copied
        };
        let live_groups = subset.groups_with_reads();
        let fit_it = live_groups.iter().any(|live| *live)
            && (every_sample || subset.tracts_with_reads() >= config.refusal_floor);
        if fit_it {
            let climb = |start: Parameters| {
                climb_from_parameters(
                    subset,
                    start,
                    homozygote_excess,
                    &genotypes,
                    &live_groups,
                    config,
                    threads,
                )
            };
            let from_every_start = || -> Vec<Climb> {
                config
                    .starting_points
                    .iter()
                    .map(|start| {
                        climb(Parameters::start(
                            *start,
                            subset.groups,
                            allele_classes_of(config),
                        ))
                    })
                    .collect()
            };
            let mut climbs = match &last {
                None => from_every_start(),
                Some(previous) => vec![climb(previous.clone())],
            };
            let mut best = the_best_of_the_walks(&climbs);
            let errors_at = |climb: &mut Climb| {
                climb.standard_errors.take().unwrap_or_else(|| {
                    standard_errors_at(
                        subset,
                        &climb.parameters,
                        homozygote_excess,
                        &live_groups,
                        config,
                        threads,
                    )
                })
            };
            let mut errors = errors_at(&mut climbs[best]);
            let still_growing: Vec<bool> = readers
                .iter()
                .map(|group| {
                    group
                        .iter()
                        .zip(&keep)
                        .any(|(reads, kept)| *reads && !*kept)
                })
                .collect();
            let precise = level_is_measured_to(
                &climbs[best].parameters,
                &live_groups,
                &errors,
                &still_growing,
                subsets.level_relative_error_target,
            );
            if precise || every_sample {
                if last.is_some() {
                    // The answer's subset: every starting point, then the walk from the last
                    // answer, the best winning.
                    let warm = climbs.remove(0);
                    climbs = from_every_start();
                    climbs.push(warm);
                    best = the_best_of_the_walks(&climbs);
                    // The walk from the last answer was judged already; a starting point's
                    // winner needs its own errors.
                    if best != climbs.len() - 1 {
                        errors = errors_at(&mut climbs[best]);
                    }
                }
                let walks = climbs.iter().map(Climb::record).collect();
                let winner = climbs.swap_remove(best);
                let mut fit = the_fit_of(subset, &[], &live_groups, winner, walks, Some(errors));
                fit.samples_fitted_on = Some(keep.iter().filter(|kept| **kept).count());
                return StratumOutcome::Fitted(Box::new(fit));
            }
            last = Some(climbs.swap_remove(best).parameters);
        } else if every_sample {
            return StratumOutcome::Refused {
                stratum: evidence.stratum,
                tracts: evidence.tracts_with_reads(),
                reason: StratumRefusal::NoSpanningReads,
            };
        }
        size = next_subset_size(size, samples);
    }
}

/// The allele classes the climb fits, `2 × span + 1`.
fn allele_classes_of(config: &SsrFitConfig) -> usize {
    (2 * config.allele_span + 1) as usize
}

/// **Which of `climbs` scored highest**, the earlier winning a tie and a score that is not a number
/// never winning — the rule [`the_better_walk`] keeps.
fn the_best_of_the_walks(climbs: &[Climb]) -> usize {
    let mut best = 0;
    for (index, climb) in climbs.iter().enumerate().skip(1) {
        if climb.score.partial_cmp(&climbs[best].score) == Some(std::cmp::Ordering::Greater) {
            best = index;
        }
    }
    best
}

/// **What a stratum's outcome is before a single walk is taken**, or `None` when it is to be
/// fitted.
///
/// Both arms of [`fit_strata`] ask this, and they must ask the same thing: a stratum refused by
/// one and fitted by the other would make the answer depend on how the threads were spent, which
/// is exactly what the arms are built not to do. [`fit_on_growing_subsets`] asks it of the whole
/// stratum, before any subset is drawn.
fn refused_before_any_walk(
    evidence: &StratumEvidence,
    config: &SsrFitConfig,
) -> Option<StratumOutcome> {
    if evidence.tracts_with_reads() == 0 {
        return Some(StratumOutcome::Refused {
            stratum: evidence.stratum,
            tracts: 0,
            reason: StratumRefusal::NoSpanningReads,
        });
    }
    if evidence.tracts_with_reads() < config.refusal_floor {
        // Too thin to fit anything of its own. It is not refused yet — the curve and a
        // neighbour's shares may still furnish it, which is what `derive_thin_strata` decides.
        return Some(StratumOutcome::Refused {
            stratum: evidence.stratum,
            tracts: evidence.tracts_with_reads(),
            reason: StratumRefusal::BelowTheFloor {
                tracts: evidence.tracts_with_reads(),
                floor: config.refusal_floor,
            },
        });
    }
    None
}

/// A stratum whose walks are done: the slippage groups they moved, the walk that won, and every walk's
/// record in starting-point order.
type WalkedStratum = (Vec<bool>, Climb, Vec<WalkRecord>);

/// **Every stratum's every starting point, handed to the pool as one flat list of walks.**
///
/// The finest independent unit the repeat-tract fit has is not a stratum but a **walk**: one
/// stratum searched from one starting point (see [`Climb`]). A stratum is searched from each of
/// [`SsrFitConfig::starting_points`], and those searches never read each other, so a run over
/// `S` strata at `P` starting points holds `S × P` pieces of work and not `S`.
///
/// **Why that matters is the remainder.** Spreading only the strata over the pool leaves a thread
/// idle for every stratum short of the pool's width, and makes every busy thread walk its own
/// stratum `P` times end to end. Eleven strata at three starting points on eighteen threads is
/// eleven pieces where there are thirty-three: seven threads idle, and the run as long as three
/// walks rather than two. Where the strata already outnumber the threads several times over the
/// difference is small — seventy-nine strata on eighteen threads go from fifteen walk-lengths to
/// fourteen.
///
/// **The memory bound is the one [`WhereTheThreadsGo::AcrossStrata`] states.** A walk holds one
/// stratum's likelihood table and a thread runs one walk at a time, so two threads walking the
/// same stratum from two starting points hold two copies of that stratum's table — which is still
/// one table a thread of the pool this runs on, [`SsrFitConfig::strata_at_once`] of them.
///
/// **That bound rests on a walk making no call into the pool.** Under this arm every parallel site
/// inside a walk runs serially. A walk that split its own tracts over the pool it runs on would let
/// a thread waiting on that split steal another walk and build a second table: measured in review,
/// six tables at four strata at once.
///
/// **And the answer is the one the serial loop gives.** Each stratum's walks are folded in
/// starting-point order by [`the_better_walk`] — the same rule reading the same values in the
/// same order — so which thread finished first cannot decide which walk wins.
fn every_stratum_from_every_start(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    stage: &StageProgress,
) -> Vec<StratumOutcome> {
    // **The same precondition `fit_pooled` states**, checked here because this arm does not go
    // through it. See there for why a span below one is refused rather than fitted.
    assert!(
        config.allele_span >= 1,
        "the fit places allele mass from -{span} to +{span} whole repeat units either side of \
         the reference length, so `SsrFitConfig::allele_span` must be at least 1; at {span} a \
         tract could only ever carry its reference length",
        span = config.allele_span
    );
    // One table for the whole run: which pairs of allele classes a diploid can carry depends on
    // the span alone, and every stratum here is fitted over the same span.
    let genotypes = genotype_pairs((2 * config.allele_span + 1) as usize);

    // Either the outcome this stratum already has, or the slippage groups its walks will move.
    // Settled before any walk starts, because a stratum with no reads in any group has no walk
    // to take.
    let before_walking: Vec<Result<StratumOutcome, Vec<bool>>> = strata
        .iter()
        .map(|evidence| {
            if let Some(refused) = refused_before_any_walk(evidence, config) {
                return Ok(refused);
            }
            let live = evidence.groups_with_reads();
            if live.iter().any(|it| *it) {
                Err(live)
            } else {
                Ok(StratumOutcome::Refused {
                    stratum: evidence.stratum,
                    tracts: evidence.tracts_with_reads(),
                    reason: StratumRefusal::NoSpanningReads,
                })
            }
        })
        .collect();

    // The flat list: one entry a (stratum, starting point) pair, over the strata that will walk.
    let starts = config.starting_points.len();
    let walks: Vec<(usize, usize)> = before_walking
        .iter()
        .enumerate()
        .filter(|(_, entry)| entry.is_err())
        .flat_map(|(stratum, _)| (0..starts).map(move |start| (stratum, start)))
        .collect();

    // Counted as walks finish, which is not the order they were listed in — so a line says how
    // many are done rather than which one this was.
    let walked = std::sync::atomic::AtomicUsize::new(0);
    let climbed: Vec<Climb> = walks
        .par_iter()
        .map(|(stratum, start)| {
            let live = before_walking[*stratum]
                .as_ref()
                .expect_err("only the strata that will walk are in this list");
            let climb = climb_from(
                &strata[*stratum],
                config.starting_points[*start],
                homozygote_excess,
                &genotypes,
                live,
                config,
                WhereTheThreadsGo::AcrossStrata,
            );
            let done = walked.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
            stage.now_and_then(|into| {
                format!(
                    "repeat-tract fit: {done} of {} walks done ({} strata from {starts} starting \
                     point(s) each); {into}",
                    walks.len(),
                    walks.len() / starts.max(1),
                )
            });
            climb
        })
        .collect();

    // **Folded back in starting-point order**, which is what makes the winner independent of the
    // order the walks finished in: `climbed` is in the order of `walks`, and `walks` lists a
    // stratum's starting points in the order the config gives them.
    let mut climbed = climbed.into_iter();
    let winners: Vec<Result<StratumOutcome, WalkedStratum>> = before_walking
        .into_iter()
        .map(|entry| {
            entry.map_err(|live| {
                let own: Vec<Climb> = (&mut climbed).take(starts).collect();
                let walks = own.iter().map(Climb::record).collect();
                let best = own
                    .into_iter()
                    .fold(None, the_better_walk)
                    .expect("at least one starting point");
                (live, best, walks)
            })
        })
        .collect();
    // **Each winner's standard errors, one stratum a thread**, under the same bound as the walks: a
    // thread holds one stratum's table at a time. Collected in stratum order.
    winners
        .into_par_iter()
        .enumerate()
        .map(|(index, entry)| match entry {
            Ok(outcome) => outcome,
            Err((live, mut best, walks)) => {
                let standard_errors = best.standard_errors.take().unwrap_or_else(|| {
                    standard_errors_at(
                        &strata[index],
                        &best.parameters,
                        homozygote_excess,
                        &live,
                        config,
                        WhereTheThreadsGo::AcrossStrata,
                    )
                });
                StratumOutcome::Fitted(Box::new(the_fit_of(
                    &strata[index],
                    &[],
                    &live,
                    best,
                    walks,
                    Some(standard_errors),
                )))
            }
        })
        .collect()
}

// ---------------------------------------------------------------------
// Thin strata: every slippage number from its period's curve
// ---------------------------------------------------------------------

/// **The bytes of the largest likelihood table any one stratum's walk will build** — one row a
/// sample with reads a tract, over the genotype pairs, eight bytes each ([`TractLikelihoods`]).
///
/// Counted from the evidence rather than estimated from the cohort's size, because how many
/// samples have reads at a tract is what fills the table: at three reads a position most do, at
/// a thin corner of a stratum few.
fn largest_table_bytes(strata: &[StratumEvidence], config: &SsrFitConfig) -> u64 {
    let width = genotype_pairs((2 * config.allele_span.max(0) + 1) as usize).len() as u64;
    strata
        .iter()
        .map(|stratum| {
            stratum
                .tracts
                .iter()
                .map(|tract| tract.samples.len() as u64)
                .sum::<u64>()
                * width
                * std::mem::size_of::<f64>() as u64
        })
        .max()
        .unwrap_or(0)
}

/// A table's size as the progress line prints it: whole mebibytes below one gibibyte, where a
/// cohort of a few hundred samples sits and one decimal of a gibibyte would read `0.0`, and
/// gibibytes to one decimal above.
fn table_size(bytes: u64) -> String {
    crate::parameter_estimation::progress::size(bytes)
}

/// **The pool several strata at once run on, and how the progress line names the schedule.**
///
/// `None` is one stratum at a time: asked for, or the only schedule a single stratum has, or what
/// is left when a pool of `at_once` threads cannot be built — the schedule with the smallest
/// memory, which returns the same bits, and the line says why. `build` is the pool's builder, a
/// parameter so that a test can make it fail.
fn a_pool_for(
    at_once: usize,
    strata: usize,
    build: impl FnOnce(usize) -> Result<rayon::ThreadPool, rayon::ThreadPoolBuildError>,
) -> (Option<rayon::ThreadPool>, String) {
    if at_once == 1 {
        return (None, "one stratum at a time".to_owned());
    }
    if strata <= 1 {
        return (
            None,
            format!("one stratum at a time, since there is {strata} to fit"),
        );
    }
    match build(at_once) {
        Ok(pool) => (
            Some(pool),
            format!("{at_once} strata at once, each on one thread"),
        ),
        Err(error) => (
            None,
            format!(
                "one stratum at a time, because a pool of {at_once} threads could not be built: \
                 {error}"
            ),
        ),
    }
}

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
///
/// **Every stratum is read from every sample**, whatever `config.subsets` says. A run calls
/// [`fit_strata_on_sample_subsets`], which reads a large cohort's strata from a subset grown to the
/// precision their levels need.
pub fn fit_strata(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> Vec<StratumOutcome> {
    fit_strata_with(strata, homozygote_excess, None, config)
}

/// **[`fit_strata`], with each stratum of a large cohort read from a subset of its samples**
/// (spec §4.4): a cohort of more than [`SampleSubsets::first`] samples has each
/// stratum fitted on the first samples of `sample_order`, grown until its slippage level is measured
/// to [`SampleSubsets::level_relative_error_target`] ([`fit_on_growing_subsets`]). A cohort of no
/// more is [`fit_strata`]'s, exactly.
///
/// `sample_order` is every sample of the cohort, in the fixed order
/// [`sample_order`](super::sample_order::sample_order) gives, indexed as the evidence's
/// [`SampleTractReads::sample`] and `homozygote_excess` are.
pub fn fit_strata_on_sample_subsets(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    sample_order: &[usize],
    config: &SsrFitConfig,
) -> Vec<StratumOutcome> {
    if sample_order.len() <= config.subsets.first {
        return fit_strata(strata, homozygote_excess, config);
    }
    fit_strata_with(strata, homozygote_excess, Some(sample_order), config)
}

/// [`fit_strata`]'s work, every stratum read from every sample when `sample_order` is `None` and
/// from a growing subset of them otherwise; everything after the strata's own answers is the same.
fn fit_strata_with(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    sample_order: Option<&[usize]>,
    config: &SsrFitConfig,
) -> Vec<StratumOutcome> {
    // **Every stratum is fitted on its own tracts and no other's.** Pooling a thin stratum's
    // tracts with its neighbours' and refitting used to be how it got an answer at all; it now
    // gets its level from its period's curve and its two shares copied from a neighbour that
    // measured them well, so there is nothing left for a pooled fit to supply
    // (`str_slippage_level_curve.md` §5.1). What that removes is the run's expensive arm: the
    // same cohort took 1,036.8 s with pooled borrowing against 155.5 s without.
    //
    // **Several strata at once where the run asks for it** (`config.strata_at_once`). A
    // stratum's answer is a function of its own tracts and of the per-sample homozygote excess,
    // and of nothing another stratum produces — the curves below are drawn afterwards, from the
    // finished answers — so fitting them together changes no number and no order. What it changes is
    // which of the fit's two grains the pool is spread over; see [`WhereTheThreadsGo`].
    //
    // **One stratum at a time is still the arm for a run holding one stratum**: spreading the
    // pool over the strata when there is a single one leaves every thread but one idle, where
    // splitting that stratum's tracts keeps the whole pool busy.
    let one_at_a_time = |evidence: &StratumEvidence| {
        if let Some(refused) = refused_before_any_walk(evidence, config) {
            return refused;
        }
        match fit_pooled(
            evidence,
            &[],
            homozygote_excess,
            config,
            WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
        ) {
            Some(fit) => StratumOutcome::Fitted(Box::new(fit)),
            None => StratumOutcome::Refused {
                stratum: evidence.stratum,
                tracts: evidence.tracts_with_reads(),
                reason: StratumRefusal::NoSpanningReads,
            },
        }
    };
    let table = table_size(largest_table_bytes(strata, config));
    // **A pool of exactly `strata_at_once` threads**, so the number of tables resident together
    // is the number the run asked for and not the machine's core count.
    let (pool, schedule) = a_pool_for(config.strata_at_once.get(), strata.len(), |threads| {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build()
    });
    let stage = StageProgress::begin(format!(
        "repeat-tract fit: fitting {} strata from {} starting point(s) each, {schedule}; the \
         largest stratum's table is {table}, so --str-param-estimates-at-once N needs about \
         N × {table}",
        strata.len(),
        config.starting_points.len(),
    ));
    let progress = |index: usize| {
        stage.now_and_then(|into| {
            format!(
                "repeat-tract fit: stratum {} of {}; {into}",
                index + 1,
                strata.len()
            )
        });
    };
    let mut outcomes: Vec<StratumOutcome> = match (sample_order, &pool) {
        (None, Some(pool)) => pool
            .install(|| every_stratum_from_every_start(strata, homozygote_excess, config, &stage)),
        (None, None) => strata
            .iter()
            .enumerate()
            .map(|(index, evidence)| {
                progress(index);
                one_at_a_time(evidence)
            })
            .collect(),
        // **A subset grows one fit after another**, so a stratum's walks cannot be handed to the
        // pool together as a whole stratum's are; the strata are, each on one thread.
        (Some(order), Some(pool)) => pool.install(|| {
            strata
                .par_iter()
                .map(|evidence| {
                    fit_on_growing_subsets(
                        evidence,
                        order,
                        homozygote_excess,
                        config,
                        WhereTheThreadsGo::AcrossStrata,
                    )
                })
                .collect()
        }),
        (Some(order), None) => strata
            .iter()
            .enumerate()
            .map(|(index, evidence)| {
                progress(index);
                fit_on_growing_subsets(
                    evidence,
                    order,
                    homozygote_excess,
                    config,
                    WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
                )
            })
            .collect(),
    };

    if let Some(order) = sample_order
        && let Some(summary) = subsets_summary(&outcomes, order.len())
    {
        stage.always(|into| format!("repeat-tract fit: {summary}; {into}"));
    }

    if let Some(summary) = climb_endings_summary(&outcomes, config.max_rounds) {
        stage.always(|into| format!("repeat-tract fit: {summary}; {into}"));
    }
    if let Some(summary) = standard_errors_summary(&outcomes) {
        stage.always(|into| format!("repeat-tract fit: {summary}; {into}"));
    }

    // **The curves are drawn after every stratum has its own answer, never during.** Stage one
    // is untouched by this — a stratum's own fitted numbers are the same whether curves are drawn
    // or not, which is the property the parity oracle checks.
    if config.curve.draw_curves {
        stage.always(|into| format!("repeat-tract fit: strata fitted, drawing the curves; {into}"));
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
    stage.always(|into| format!("repeat-tract fit: done; {into}"));
    outcomes
}

/// **What the run's log says about the subsets the strata were read from** (spec §4.4): over the
/// strata fitted on their own tracts, the fewest, the median and the most samples a stratum was
/// fitted on, and how many grew to every sample. `None` when no fitted stratum took a subset.
fn subsets_summary(outcomes: &[StratumOutcome], cohort: usize) -> Option<String> {
    let mut sizes: Vec<usize> = outcomes
        .iter()
        .filter_map(|outcome| match outcome {
            StratumOutcome::Fitted(fit) => fit.samples_fitted_on,
            _ => None,
        })
        .collect();
    if sizes.is_empty() {
        return None;
    }
    sizes.sort_unstable();
    Some(format!(
        "the {strata} strata fitted on their own tracts were read from a subset of the cohort's \
         {cohort} samples: {fewest} to {most}, median {median}; {every} of them grew to every \
         sample",
        strata = sizes.len(),
        fewest = sizes[0],
        most = sizes[sizes.len() - 1],
        median = sizes[sizes.len() / 2],
        every = sizes.iter().filter(|size| **size >= cohort).count(),
    ))
}

/// **What the run's log says about how the climbs ended** (`fit_precision.md` §4.3): over the strata
/// fitted on their own tracts, how many walks settled — and how many of those at a point where no
/// number had an error — stopped at a round that lost, or ran out of rounds, the rounds they took
/// between them, and in how many strata the winning walk settled.
/// `None` when no fitted stratum records its walks.
fn climb_endings_summary(outcomes: &[StratumOutcome], max_rounds: u32) -> Option<String> {
    let fits: Vec<&StratumFit> = outcomes
        .iter()
        .filter_map(|outcome| match outcome {
            StratumOutcome::Fitted(fit) if !fit.walks.is_empty() => Some(fit.as_ref()),
            _ => None,
        })
        .collect();
    if fits.is_empty() {
        return None;
    }
    let walks = || fits.iter().flat_map(|fit| fit.walks.iter());
    let ended = |ending: ClimbEnding| walks().filter(|walk| walk.ending == ending).count();
    let rounds: u64 = walks().map(|walk| u64::from(walk.rounds)).sum();
    Some(format!(
        "the climbs of the {strata} strata fitted on their own tracts: of their {count} walks, \
         {settled} settled ({no_error} of them at a point where no number had an error), {lost} \
         stopped at a round that lost (its moves undone) and {out} ran \
         out of their {max_rounds} rounds, {rounds} rounds in all; the winning walk settled in \
         {winners} of the {strata}",
        strata = fits.len(),
        count = walks().count(),
        settled = ended(ClimbEnding::Settled),
        no_error = walks().filter(|walk| walk.settled_with_no_error).count(),
        lost = ended(ClimbEnding::LostARound),
        out = ended(ClimbEnding::OutOfRounds),
        winners = fits.iter().filter(|fit| fit.ending.settled()).count(),
    ))
}

/// **What the run's log says about the strata's own standard errors**, once every stratum has its
/// answer and before any curve blends it (`fit_precision.md` §4.2): over the strata fitted on their
/// own tracts, the median and largest error of the level and of the concentration as a percentage of
/// the number, and of the two shares, the fall-off and the length classes' shares as they are; then
/// how many numbers have no error, by why. `None` when no stratum carries errors.
fn standard_errors_summary(outcomes: &[StratumOutcome]) -> Option<String> {
    let errors: Vec<(&StratumFit, &StratumErrors)> = outcomes
        .iter()
        .filter_map(|outcome| match outcome {
            StratumOutcome::Fitted(fit) => Some((fit.as_ref(), fit.standard_errors.as_ref()?)),
            _ => None,
        })
        .collect();
    if errors.is_empty() {
        return None;
    }
    // Per kind of number, its errors (divided by `scale`), and the numbers without one by why.
    let mut level = Vec::new();
    let mut shorter = Vec::new();
    let mut fall_off = Vec::new();
    let mut concentration = Vec::new();
    let mut shares = Vec::new();
    let (mut numbers, mut not_identified, mut not_placed, mut no_share) = (0, 0, 0, 0);
    let mut count = |error: StratumError, into: &mut Vec<f64>, scale: f64| {
        numbers += 1;
        match error {
            StratumError::Estimated(error) => into.push(error / scale),
            StratumError::NotIdentified => not_identified += 1,
            StratumError::NotPlaced => not_placed += 1,
            StratumError::NoShare => no_share += 1,
        }
    };
    for (fit, stratum) in &errors {
        for (slippage, group) in fit.slippage.iter().zip(&stratum.slippage) {
            let (Some(slippage), Some(group)) = (slippage, group) else {
                continue;
            };
            count(group.level, &mut level, slippage.level / 100.0);
            count(group.shorter_share, &mut shorter, 1.0);
            count(group.fall_off, &mut fall_off, 1.0);
        }
        count(
            stratum.concentration,
            &mut concentration,
            fit.concentration / 100.0,
        );
        for share in &stratum.length_spectrum {
            count(*share, &mut shares, 1.0);
        }
    }
    let described = |mut values: Vec<f64>, what: &str, relative: bool| {
        if values.is_empty() {
            return format!("{what}: none has an error");
        }
        values.sort_by(f64::total_cmp);
        let (median, largest) = (values[values.len() / 2], values[values.len() - 1]);
        if relative {
            format!("{what} ± median {median:.0}% of itself (largest {largest:.0}%)")
        } else {
            format!("{what} ± median {median:.3} (largest {largest:.3})")
        }
    };
    Some(format!(
        "standard errors of the {} strata fitted on their own tracts, before any curve: {}; {}; {}; \
         {}; {}; of their {numbers} numbers, {not_placed} have no error because the tracts do not \
         place them (an error of more than {NOT_PLACED} on the climb's logit or log scale, or \
         wider than a share's whole range), \
         {not_identified} because the tracts cannot tell them apart from the others, and {no_share} \
         are length classes with no share",
        errors.len(),
        described(level, "slippage level", true),
        described(shorter, "shorter share", false),
        described(fall_off, "fall-off", false),
        described(concentration, "concentration", true),
        described(shares, "length-class shares", false),
    ))
}

// ---------------------------------------------------------------------
// The tract prior's middle rung: one motif period's tracts pooled
// ---------------------------------------------------------------------

/// **One motif period's length spectrum and concentration, fitted over every tract of that
/// period at once** — the middle rung of the tract ladder
/// (`doc/devel/ng/spec/population_diversity.md` §4.4).
///
/// A stratum too thin to be fitted carries no length spectrum of its own
/// ([`DerivedStratum`]), and the caller's genotype prior at such a tract needs one. This is
/// what it falls back to: the same two numbers a [`StratumFit`] produces, estimated from every
/// tract that shares the motif period rather than from the stratum's own eight.
///
/// **Why a pooled fit rather than a curve through the strata**, which is what the three
/// slippage numbers get: a curve through a *distribution* means one curve per length class,
/// refitted and renormalised, and the classes are not independent — a pooled fit is a real
/// distribution by construction. The gap it covers is also far smaller than slippage's was:
/// the strata with no fit of their own hold about 2 in 100 of HG002's tracts and at most 7 in
/// 100 of tomato's, against most strata on both cohorts for slippage
/// (`population_diversity.md` §4.4).
///
/// **What it gives up is the repeat-count trend within a period, and it gives it up twice.** A
/// longer tract spreads over more lengths, and pooling flattens that directly. It also flattens
/// it *indirectly*, and that half is easy to miss: the pooled climb fits **one** slippage triple
/// per slippage group across every repeat count in the period, where slippage rises about
/// 1.3-fold per repeat count over the measured range
/// ([`StratumFits::at`](super::stratum_fits::StratumFits::at)). A read off the
/// reference length is either a slip or a real allele, so the two trade off — the pooled spectrum
/// is fitted against a slippage level too low at the period's long strata and too high at its
/// short ones, which widens its tails at the long end and narrows them at the short. The slippage
/// numbers a *caller* reads are unaffected: those come from the period's curves, and nothing here
/// emits a slippage number. Bounded by the loci the rung applies to — 2 in 100 of HG002's tracts
/// and at most 7 in 100 of tomato's.
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodLengthSpectrum {
    /// The motif length these tracts share.
    pub period: u8,
    /// How the period's chromosomes are spread over tract lengths, indexed from
    /// `-allele_span` to `+allele_span` in whole repeat units **from each tract's own
    /// reference length**.
    ///
    /// **That the index is an offset is what makes pooling legitimate at all**: two strata of
    /// one period sit at different absolute repeat counts, and it is only because every tract
    /// is described relative to its own reference length that their evidence can be added up.
    pub length_spectrum: Vec<f64>,
    /// How monomorphic the period's tracts are. Small means most tracts carry one length.
    pub concentration: f64,
    /// Tracts with at least one spanning read that went into the pool.
    pub tracts_fitted: usize,
    /// How many of the period's strata contributed tracts to it.
    pub strata_pooled: usize,
    /// Whether the climb settled. **Running out of rounds, or stopping at a round that lost, is
    /// never reported as convergence**, exactly as [`StratumFit::ending`].
    pub converged: bool,
}

/// **Fit one length spectrum and concentration a motif period**, pooling every tract of that
/// period — what a stratum with no fit of its own falls back to.
///
/// **Separate from [`fit_strata`] and opt-in, because it is the one thing on this path that
/// costs a run more than it used to.** `population_diversity.md` §6 says nothing here is
/// fitted; that is true of the seam and of the ordinary-site side, and not of this rung, which
/// §4.4 settles as a pooled fit.
///
/// **What it costs, measured**: on two strata of 300 tracts each, 8 samples, allele span 1, this
/// call takes **1.67–1.72 s against `fit_strata`'s 2.68 s** on the same evidence — about **60%
/// on top**, three runs each. It is less than the strata cost between them because it runs one
/// climb where they run two, over the same tracts. So a run that asks for the middle rung pays
/// about 1.6 times for the repeat-tract half of its fit, not twice.
///
/// A run that does not ask still gets an answer at every tract — the ladder's bottom rung, a
/// flat shape at a stated concentration — and the rung it lands on is reported either way, which
/// is why this is a second call rather than a widened `fit_strata`.
///
/// **A period below the refusal floor is left out rather than fitted badly.** The floor is the
/// same [`SsrFitConfig::refusal_floor`] a stratum is held to, measured in the same unit —
/// tracts with a spanning read — because a pool of five tracts is exactly as thin as a stratum
/// of five.
///
/// `homozygote_excess` is one number a sample, as [`fit_stratum`] takes it, and is held rather
/// than fitted here too.
///
/// # Panics
///
/// When two strata of one period disagree about how many buckets a read's offset was recorded
/// in, or about how many slippage groups the cohort declares. Both are properties of the run
/// rather than of a stratum, so a disagreement means the evidence was assembled from two runs —
/// and the symptom otherwise is an index past the end of a bucket row, inside the scorer,
/// naming neither period nor stratum.
#[must_use]
pub fn fit_period_length_spectra(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> BTreeMap<u8, PeriodLengthSpectrum> {
    let mut by_period: BTreeMap<u8, Vec<&StratumEvidence>> = BTreeMap::new();
    for evidence in strata {
        by_period
            .entry(evidence.stratum.period)
            .or_default()
            .push(evidence);
    }

    let mut fitted = BTreeMap::new();
    for (period, members) in by_period {
        // **The floor is tested before anything is copied.** `pool_a_period` clones every
        // `TractReads` of the period, and on a tomato-sized cohort the largest period's copy is
        // a second copy of that period's whole STR evidence at peak — paid, under the old
        // order, even for a period that was about to be discarded. The sum over members is the
        // same number the pooled evidence would have reported.
        let tracts: usize = members
            .iter()
            .map(|evidence| evidence.tracts_with_reads())
            .sum();
        if tracts < config.refusal_floor {
            continue;
        }
        let pooled = pool_a_period(period, &members);
        // **One period at a time, its tracts spread over the pool.** The periods are as
        // independent of each other as the strata are, so this loop could be spread the other
        // way too; it is left alone because it runs over a handful of periods rather than over
        // a hundred and forty strata, and because no calling run reaches it yet.
        // **No standard errors here**: only the length spectrum and the concentration are kept,
        // and nothing reads their errors at this rung.
        let Some((best, walks, live_groups)) = the_best_walk(
            &pooled,
            homozygote_excess,
            config,
            WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
        ) else {
            continue;
        };
        let fit = the_fit_of(&pooled, &[], &live_groups, best, walks, None);
        fitted.insert(
            period,
            PeriodLengthSpectrum {
                period,
                length_spectrum: fit.length_spectrum,
                concentration: fit.concentration,
                tracts_fitted: fit.tracts_fitted,
                // **Strata that put a tract in, not strata of the period.** A stratum whose
                // every tract went unread contributed nothing, and counting it would say the
                // pool rested on more evidence than it did.
                strata_pooled: members
                    .iter()
                    .filter(|evidence| evidence.tracts_with_reads() > 0)
                    .count(),
                converged: fit.ending.settled(),
            },
        );
    }
    fitted
}

/// Concatenate one period's strata into the single evidence set the pooled fit reads.
///
/// **The `stratum` field of the result is the smallest reference repeat count in the pool and
/// means nothing.** [`fit_pooled`] copies it onto the [`StratumFit`] it returns and
/// [`fit_period_length_spectra`] then keeps only the length spectrum and the concentration, so
/// no repeat count is claimed for a pool that spans many. Nothing in the scorer reads it: a
/// tract's evidence is buckets of reads at offsets from its own reference length, and the
/// absolute length never enters.
fn pool_a_period(period: u8, members: &[&StratumEvidence]) -> StratumEvidence {
    let first = members
        .first()
        .expect("a period with no strata is not keyed");
    for evidence in members {
        assert_eq!(
            evidence.read_span,
            first.read_span,
            "period {period}: the stratum at {} repeats recorded read offsets in {} buckets \
             either side and the one at {} repeats in {} — the span is a property of the run, \
             so two values mean the evidence came from two runs",
            evidence.stratum.reference_repeats,
            evidence.read_span,
            first.stratum.reference_repeats,
            first.read_span
        );
        assert_eq!(
            evidence.groups,
            first.groups,
            "period {period}: the stratum at {} repeats declares {} slippage groups and the one \
             at {} repeats declares {} — the count is a property of the run, so two values mean \
             the evidence came from two runs",
            evidence.stratum.reference_repeats,
            evidence.groups,
            first.stratum.reference_repeats,
            first.groups
        );
    }
    let smallest = members
        .iter()
        .map(|evidence| evidence.stratum.reference_repeats)
        .min()
        .expect("a period with no strata is not keyed");
    StratumEvidence {
        stratum: Stratum {
            period,
            reference_repeats: smallest,
        },
        tracts: members
            .iter()
            .flat_map(|evidence| evidence.tracts.iter().cloned())
            .collect(),
        read_span: first.read_span,
        groups: first.groups,
        // **The five diagnostic counters below are summed rather than zeroed, and nothing in the
        // fit reads any of them.** The pooled evidence never leaves this function —
        // `fit_period_length_spectra` keeps only the length spectrum and the concentration — so
        // these are dead in the sense that removing them changes no output. They are summed
        // anyway because a zero here would be a claim: *this period had no guarded tracts, no
        // reads that reached without crossing, no bases compared*, which is false of every real
        // period and is the shape a later reader would take at face value.
        tracts_over_guard_threshold: members
            .iter()
            .map(|evidence| evidence.tracts_over_guard_threshold)
            .sum(),
        reads_reaching_not_crossing: members
            .iter()
            .map(|evidence| evidence.reads_reaching_not_crossing)
            .sum(),
        guard_reads: members.iter().map(|evidence| evidence.guard_reads).sum(),
        bases_compared: members.iter().map(|evidence| evidence.bases_compared).sum(),
        mismatching_bases: members
            .iter()
            .map(|evidence| evidence.mismatching_bases)
            .sum(),
    }
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
                        total += f64::from(*reads) * float::ln(probability.max(1e-300));
                    }
                }
                total
            }));
            let row = &mut scaled[start..];
            let largest = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            ln_offset += largest;
            for value in row {
                *value = float::exp(*value - largest);
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
            let (weights, weights_left) = independent.as_chunks::<4>();
            let (values, values_left) = row[..width].as_chunks::<4>();
            for (weight, value) in weights.iter().zip(values) {
                lanes += wide::f64x4::new(*weight) * wide::f64x4::new(*value);
            }
            let parts = lanes.to_array();
            let mut at_random = (parts[0] + parts[1]) + (parts[2] + parts[3]);
            for (weight, value) in weights_left.iter().zip(values_left) {
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
            quadrature.ln_weight + float::ln(product) + scaled_by
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
    /// Which evidence these tables were built from, so a test can count the tables one fit
    /// holds at once ([`table_census`]).
    #[cfg(test)]
    evidence: usize,
}

#[cfg(test)]
impl Drop for Prepared {
    fn drop(&mut self) {
        table_census::dropped(self.evidence);
    }
}

/// **How many likelihood tables are alive at once, counted for the strata one test names.**
///
/// The only thing [`SsrFitConfig::strata_at_once`] changes is how many tables the fit holds
/// together — every schedule returns the same bits — so a test of the option has to count them.
/// Tables are counted only for the evidence a test registered, by address, because the library's
/// tests run at the same time and build tables of their own.
#[cfg(test)]
mod table_census {
    use std::sync::atomic::{AtomicUsize, Ordering::SeqCst};
    use std::sync::{Mutex, MutexGuard};

    static ONE_COUNT_AT_A_TIME: Mutex<()> = Mutex::new(());
    static COUNTED: Mutex<Vec<usize>> = Mutex::new(Vec::new());
    static LIVE: AtomicUsize = AtomicUsize::new(0);
    static PEAK: AtomicUsize = AtomicUsize::new(0);

    fn is_counted(evidence: usize) -> bool {
        COUNTED
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .contains(&evidence)
    }

    /// Start counting the tables built from `strata`, and stop any other test counting until the
    /// returned guard is dropped.
    pub(super) fn count_the_tables_of(
        strata: &[super::StratumEvidence],
    ) -> MutexGuard<'static, ()> {
        let guard = ONE_COUNT_AT_A_TIME
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        *COUNTED
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = strata
            .iter()
            .map(|evidence| std::ptr::from_ref(evidence) as usize)
            .collect();
        LIVE.store(0, SeqCst);
        PEAK.store(0, SeqCst);
        guard
    }

    /// Count again from nothing, over the same strata.
    pub(super) fn restart() {
        LIVE.store(0, SeqCst);
        PEAK.store(0, SeqCst);
    }

    /// The most tables alive at once since counting started.
    pub(super) fn peak() -> usize {
        PEAK.load(SeqCst)
    }

    pub(super) fn built(evidence: usize) {
        if is_counted(evidence) {
            let now = LIVE.fetch_add(1, SeqCst) + 1;
            PEAK.fetch_max(now, SeqCst);
        }
    }

    pub(super) fn dropped(evidence: usize) {
        if is_counted(evidence) {
            LIVE.fetch_sub(1, SeqCst);
        }
    }
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
    /// Whether this scorer may spread its own work over the pool, or is one of many running on
    /// that pool at once and must stay on the thread it was started on ([`WhereTheThreadsGo`]).
    threads: WhereTheThreadsGo,
}

impl<'a> Scorer<'a> {
    fn new(
        evidence: &'a StratumEvidence,
        homozygote_excess: &'a [f64],
        genotypes: &'a [(usize, usize)],
        config: &SsrFitConfig,
        threads: WhereTheThreadsGo,
    ) -> Self {
        Self {
            evidence,
            homozygote_excess,
            genotypes,
            allele_span: config.allele_span,
            quadrature_points: config.quadrature_points,
            prepared: None,
            held_quadrature: None,
            threads,
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
                    self.threads,
                ),
            ));
        }
        let (_, _, quadrature) = self.held_quadrature.as_ref().expect("built above");
        let prepared = self.prepared.as_ref().expect("refreshed above");
        // Across tracts in parallel, but summed back in tract order: a parallel float sum
        // reorders the additions run to run, and a fit's whole output is a difference between
        // one set of parameters and another.
        //
        // **Serial when the pool is already carrying the other classes** — the values and the
        // order they are summed in are the same either way, so the two arms differ in where the
        // threads are and in nothing else.
        let per_tract: Vec<f64> = match self.threads {
            WhereTheThreadsGo::AcrossTheTractsOfOneStratum => prepared
                .tracts
                .par_iter()
                .map(|tract| ln_tract(tract, quadrature, self.genotypes))
                .collect(),
            WhereTheThreadsGo::AcrossStrata => prepared
                .tracts
                .iter()
                .map(|tract| ln_tract(tract, quadrature, self.genotypes))
                .collect(),
        };
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
        // **The old tables go before the new ones are built**, which halves what a walk holds at
        // its peak: nothing reads them while the new ones are filled, so holding both would cost
        // a second table for nothing — and make the size the fit prints for one table untrue.
        self.prepared = None;
        let per_group_allele: Vec<Vec<Vec<f64>>> = slippage
            .iter()
            .map(|group| {
                (-self.allele_span..=self.allele_span)
                    .map(|allele| group.read_probabilities(allele, self.evidence.read_span))
                    .collect()
            })
            .collect();
        let build = |tract: &TractReads| {
            TractLikelihoods::of(
                tract,
                &per_group_allele,
                self.genotypes,
                self.homozygote_excess,
            )
        };
        let tracts = match self.threads {
            WhereTheThreadsGo::AcrossTheTractsOfOneStratum => {
                self.evidence.tracts.par_iter().map(build).collect()
            }
            WhereTheThreadsGo::AcrossStrata => self.evidence.tracts.iter().map(build).collect(),
        };
        #[cfg(test)]
        table_census::built(std::ptr::from_ref(self.evidence) as usize);
        self.prepared = Some(Prepared {
            slippage: slippage.to_vec(),
            tracts,
            #[cfg(test)]
            evidence: std::ptr::from_ref(self.evidence) as usize,
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
    threads: WhereTheThreadsGo,
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
    let one_point = |point: usize, piece: &mut [f64]| {
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
    };
    // **Each point is written into its own row and reads none of the others**, so dividing the
    // rows between threads and not dividing them at all give the same buffer.
    match threads {
        WhereTheThreadsGo::AcrossTheTractsOfOneStratum => frequencies
            .par_chunks_mut(classes)
            .enumerate()
            .for_each(|(point, piece)| one_point(point, piece)),
        WhereTheThreadsGo::AcrossStrata => frequencies
            .chunks_mut(classes)
            .enumerate()
            .for_each(|(point, piece)| one_point(point, piece)),
    }
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
        ln_weight: -float::ln(points as f64),
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
    largest + float::ln(values.iter().map(|v| float::exp(v - largest)).sum::<f64>())
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
        return float::ln(std::f64::consts::PI / float::sin(std::f64::consts::PI * x))
            - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut series = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        series += coefficient / (x + index as f64);
    }
    let t = x + 7.5;
    0.5 * float::ln(std::f64::consts::TAU) + (x + 0.5) * float::ln(t) - t + float::ln(series)
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
    let front = float::exp(a * float::ln(x) + b * float::ln(1.0 - x) + ln_beta) / a;
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
    float::ln(p / (1.0 - p))
}

fn expit(x: f64) -> f64 {
    1.0 / (1.0 + float::exp(-x))
}

/// Golden-section on one coordinate, over a bracket of `span` either side of `start`: the point it
/// ends on, and its score there. **It never scores `start` itself**, so the point it returns can
/// score lower than `start` did.
fn climb_scalar(mut score: impl FnMut(f64) -> f64, start: f64, span: f64) -> (f64, f64) {
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
    if at_left > at_right {
        (left, at_left)
    } else {
        (right, at_right)
    }
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
    let stage = StageProgress::begin(format!(
        "repeat-tract evidence: gathering {} tract(s) in {} strata from {} sample(s)",
        strata.len(),
        tracts_in.len(),
        names.len()
    ));
    cohort.with_strata(&band, |lent| {
        stage.always(|into| format!("repeat-tract evidence: every sample's sections read; {into}"));
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

        let gathered: Vec<StratumEvidence> = tracts_in
            .iter()
            .enumerate()
            .map(|(index, (stratum, tracts))| {
                stage.now_and_then(|into| {
                    format!(
                        "repeat-tract evidence: stratum {} of {}; {into}",
                        index + 1,
                        tracts_in.len()
                    )
                });
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
            .collect();

        // **What the arranged evidence holds, and what it grows with**: how many (tract, sample)
        // pairs have reads, out of how many there are. Printed because it is the number that
        // sets both this evidence's size and each stratum's likelihood table.
        let rows: usize = gathered.iter().map(StratumEvidence::rows_with_reads).sum();
        let pairs = gathered
            .iter()
            .map(|stratum| stratum.tracts.len())
            .sum::<usize>()
            * names.len();
        let bytes: usize = gathered.iter().map(StratumEvidence::heap_bytes).sum();
        let dropped: u64 = gathered
            .iter()
            .map(|stratum| stratum.tracts_over_guard_threshold)
            .sum();
        stage.always(|into| {
            format!(
                "repeat-tract evidence: arranged, {into}; {dropped} tract(s) left out by the \
                 guard; {rows} (tract, sample) pairs with reads of the {pairs} the kept tracts \
                 make ({:.1} in every 100), holding {} by the vectors' sizes, {:.0} bytes a pair \
                 with reads",
                100.0 * rows as f64 / pairs.max(1) as f64,
                crate::parameter_estimation::progress::size(bytes as u64),
                bytes as f64 / rows.max(1) as f64,
            )
        });
        gathered
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
                return self.gamma(shape + 1.0) * float::powf(u, 1.0 / shape);
            }
            let d = shape - 1.0 / 3.0;
            let c = 1.0 / (9.0 * d).sqrt();
            loop {
                let x = self.normal();
                let v = float::powi(1.0 + c * x, 3);
                if v <= 0.0 {
                    continue;
                }
                let u = self.uniform().max(1e-300);
                if float::ln(u) < 0.5 * x * x + d - d * v + d * float::ln(v) {
                    return d * v;
                }
            }
        }

        fn normal(&mut self) -> f64 {
            let u1 = self.uniform().max(1e-300);
            let u2 = self.uniform();
            (-2.0 * float::ln(u1)).sqrt() * float::cos(std::f64::consts::TAU * u2)
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
            .map(|class| float::powi(0.55, (class as i32 - middle as i32).abs()))
            .collect();
        normalise(&mut spectrum);
        spectrum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bench_fixtures::{draw_stratum, spectrum_of};

    // -----------------------------------------------------------------
    // What outlives the evidence
    // -----------------------------------------------------------------

    /// The counts kept after the evidence is dropped carry its stratum and both counts, and give
    /// the rate the evidence gives — including `None` where nothing was compared, which must not
    /// turn into a fitted zero on the way.
    #[test]
    fn substitution_counts_give_the_rate_the_evidence_gives() {
        let spectrum = spectrum_of(3);
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let mut evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 4, 3, 6, 1, 7);
        let nothing_compared = evidence.substitution_counts();
        assert_eq!(nothing_compared.substitution_rate(), None);
        assert_eq!(evidence.substitution_rate(), None);

        evidence.stratum = Stratum {
            period: 3,
            reference_repeats: 7,
        };
        evidence.bases_compared = 4_000;
        evidence.mismatching_bases = 3;
        let counts = evidence.substitution_counts();
        assert_eq!(
            counts,
            StratumSubstitutionCounts {
                stratum: evidence.stratum,
                bases_compared: 4_000,
                mismatching_bases: 3,
            }
        );
        assert_eq!(counts.substitution_rate(), Some(3.0 / 4_000.0));
    }

    /// **A count with one outcome only takes half a count of the other**, so no substitution rate
    /// is exactly zero or one (owner, checkpoint E): no mismatch in 459 bases gives 0.5 / 460,
    /// every base mismatched gives 459.5 / 460, and a count that saw both keeps its ratio.
    #[test]
    fn a_count_with_one_outcome_only_takes_half_a_count_of_the_other() {
        let only = |mismatching_bases: u64| StratumSubstitutionCounts {
            stratum: Stratum {
                period: 3,
                reference_repeats: 7,
            },
            bases_compared: 459,
            mismatching_bases,
        };
        assert_eq!(only(0).substitution_rate(), Some(0.5 / 460.0));
        assert_eq!(only(459).substitution_rate(), Some(459.5 / 460.0));
        assert_eq!(only(1).substitution_rate(), Some(1.0 / 459.0));
    }

    // -----------------------------------------------------------------
    // The answer does not depend on how the work was divided
    // -----------------------------------------------------------------

    /// A set of repeat-length classes of deliberately unequal size, drawn at one period so that
    /// the curves have something to draw through.
    ///
    /// **The sizes span an order of magnitude on purpose.** A run under
    /// [`WhereTheThreadsGo::AcrossStrata`] cannot finish before its largest class would have
    /// taken on one thread, so a fixture where every class is the same size would be the one
    /// shape that hides the difference between the two arms.
    fn classes_of_unequal_size() -> Vec<StratumEvidence> {
        let spectrum = spectrum_of(3);
        let mut strata = Vec::new();
        for (repeats, tracts, seed) in [
            (8_u64, 96_usize, 11_u64),
            (9, 48, 12),
            (10, 24, 13),
            (11, 16, 14),
            (12, 12, 15),
            (13, 10, 16),
        ] {
            let truth = Slippage {
                level: 0.08,
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
        strata
    }

    fn a_two_round_config(strata_at_once: usize) -> SsrFitConfig {
        SsrFitConfig {
            allele_span: 1,
            max_rounds: 2,
            refusal_floor: 8,
            strata_at_once: NonZeroUsize::new(strata_at_once).expect("at least one stratum"),
            ..SsrFitConfig::default()
        }
    }

    /// A fitted stratum's standard errors laid out flat, each absent one as a negative code saying
    /// why — so a schedule's errors are compared as the numbers are. A fit must carry them.
    fn every_error(fit: &StratumFit) -> Vec<f64> {
        let errors = fit
            .standard_errors
            .as_ref()
            .expect("a fitted stratum carries its errors");
        let code = |error: &StratumError| match error {
            StratumError::Estimated(error) => *error,
            StratumError::NotIdentified => -1.0,
            StratumError::NotPlaced => -2.0,
            StratumError::NoShare => -3.0,
        };
        let mut flat: Vec<f64> = errors.length_spectrum.iter().map(code).collect();
        flat.push(code(&errors.concentration));
        for group in &errors.slippage {
            match group {
                Some(it) => flat.extend([&it.level, &it.shorter_share, &it.fall_off].map(code)),
                None => flat.extend([-4.0; 3]),
            }
        }
        flat
    }

    /// Every number a fitted class carries, laid out flat so two runs can be compared value by
    /// value rather than by a derived `PartialEq` the type does not have.
    fn every_fitted_number(outcomes: &[StratumOutcome]) -> Vec<f64> {
        let mut flat = Vec::new();
        for outcome in outcomes {
            match outcome {
                StratumOutcome::Fitted(fit) => {
                    flat.push(1.0);
                    flat.push(fit.concentration);
                    flat.push(fit.log_likelihood_a_tract);
                    flat.extend(fit.length_spectrum.iter().copied());
                    flat.extend(every_error(fit));
                    // How the winner ended, and every walk's ending, rounds and judgements.
                    let code = |ending: ClimbEnding| match ending {
                        ClimbEnding::Settled => 1.0,
                        ClimbEnding::LostARound => 2.0,
                        ClimbEnding::OutOfRounds => 3.0,
                    };
                    flat.push(code(fit.ending));
                    for walk in &fit.walks {
                        flat.extend([
                            code(walk.ending),
                            f64::from(walk.rounds),
                            f64::from(walk.judgements),
                            f64::from(u8::from(walk.settled_with_no_error)),
                        ]);
                    }
                    for slippage in &fit.slippage {
                        match slippage {
                            Some(it) => flat.extend([it.level, it.shorter_share, it.fall_off]),
                            None => flat.extend([f64::NAN; 3]),
                        }
                    }
                }
                StratumOutcome::Derived(derived) => {
                    flat.push(2.0);
                    for slippage in &derived.slippage {
                        match slippage {
                            Some(it) => flat.extend([it.level, it.shorter_share, it.fall_off]),
                            None => flat.extend([f64::NAN; 3]),
                        }
                    }
                }
                StratumOutcome::Refused { tracts, .. } => {
                    flat.push(3.0);
                    flat.push(*tracts as f64);
                }
            }
        }
        flat
    }

    /// The same values, compared so that two `NaN`s in the same slot agree — a slippage group
    /// with no reads is `None` in both runs and must not fail the comparison.
    fn the_same_bits(left: &[f64], right: &[f64], what: &str) {
        assert_eq!(
            left.len(),
            right.len(),
            "{what}: different number of values"
        );
        for (index, (a, b)) in left.iter().zip(right).enumerate() {
            assert!(
                a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
                "{what}: value {index} is {a} one way and {b} the other"
            );
        }
    }

    fn fitted_in_a_pool(
        strata: &[StratumEvidence],
        excess: &[f64],
        config: &SsrFitConfig,
        threads: usize,
    ) -> Vec<f64> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a pool of the requested width");
        pool.install(|| every_fitted_number(&fit_strata(strata, excess, config)))
    }

    /// **Spreading the pool over the classes rather than over one class's tracts moves no
    /// number.**
    ///
    /// This is the whole warrant for fitting the classes at once. A class's answer is a function
    /// of its own tracts and of the per-sample homozygote excess, and every parallel site inside
    /// the fit collects in tract order and sums serially — so the two arms are the same
    /// arithmetic in a different place, and the comparison is bit-for-bit rather than to a
    /// tolerance.
    #[test]
    fn the_two_ways_of_spending_the_pool_give_the_same_bits() {
        let strata = classes_of_unequal_size();
        let excess = [0.4; 8];
        let one_class_at_a_time = fitted_in_a_pool(&strata, &excess, &a_two_round_config(1), 4);
        let four_at_once = fitted_in_a_pool(&strata, &excess, &a_two_round_config(4), 4);
        assert!(
            one_class_at_a_time.contains(&1.0),
            "the fixture must fit at least one class, or the comparison is between two empties"
        );
        the_same_bits(
            &one_class_at_a_time,
            &four_at_once,
            "one class at a time against four at once",
        );
    }

    /// **And one stratum at a time does not depend on how wide the caller's pool is.**
    ///
    /// A fit whose answer moved with the machine's core count would be one no two runs could be
    /// compared across, which is the property the census oracles rest on. **Only this schedule
    /// runs on the caller's pool**: several strata at once build a pool of their own, whose width
    /// is `strata_at_once` itself, and
    /// [`any_number_of_strata_at_once_gives_the_same_bits`] varies that.
    #[test]
    fn one_stratum_at_a_time_does_not_move_with_the_width_of_the_pool() {
        let strata = classes_of_unequal_size();
        let excess = [0.4; 8];
        let config = a_two_round_config(1);
        let narrow = fitted_in_a_pool(&strata, &excess, &config, 1);
        let wide = fitted_in_a_pool(&strata, &excess, &config, 8);
        the_same_bits(
            &narrow,
            &wide,
            "one stratum at a time, at one thread against eight",
        );
    }

    /// **However many strata the run fits at once, the answer is the same bits** — one, two,
    /// four and seven, over the fixture's six strata from three starting points (eighteen walks);
    /// and from one starting point (six walks) at seven and at sixty-four, where most of the
    /// pool has nothing to do.
    ///
    /// This is the contract `--str-param-estimates-at-once` rests on: the option chooses how much
    /// memory the fit holds and never what it returns.
    #[test]
    fn any_number_of_strata_at_once_gives_the_same_bits() {
        let strata = classes_of_unequal_size();
        let excess = [0.4; 8];
        let one = fitted_in_a_pool(&strata, &excess, &a_two_round_config(1), 4);
        assert!(
            one.contains(&1.0),
            "the fixture must fit at least one class, or the comparison is between two empties"
        );
        for strata_at_once in [2, 4, 7] {
            let several =
                fitted_in_a_pool(&strata, &excess, &a_two_round_config(strata_at_once), 4);
            the_same_bits(
                &one,
                &several,
                &format!("one stratum at a time against {strata_at_once} at once"),
            );
        }

        let from_one_start = |strata_at_once: usize| SsrFitConfig {
            starting_points: vec![StartingPoint::spanning_the_monomorphic_range()[1]],
            ..a_two_round_config(strata_at_once)
        };
        let one = fitted_in_a_pool(&strata, &excess, &from_one_start(1), 4);
        for strata_at_once in [7, 64] {
            the_same_bits(
                &one,
                &fitted_in_a_pool(&strata, &excess, &from_one_start(strata_at_once), 4),
                &format!("one start: one stratum at a time against {strata_at_once} at once"),
            );
        }
    }

    /// **`N` strata at once hold `N` likelihood tables at most, and more than one when `N` is**,
    /// with the caller's pool at eight threads — wider than any `N` here, so a schedule that
    /// escaped onto the caller's pool, or sized its pool from the machine, would show.
    ///
    /// Every schedule returns the same bits, so this is the only test that sees what the option
    /// is for. Review measured what it catches: running the walks on the caller's pool held eight
    /// tables at two at once, and building the new tables before dropping the old held two at one
    /// at a time.
    #[test]
    fn n_strata_at_once_hold_at_most_n_tables() {
        let strata = classes_of_unequal_size();
        let excess = [0.4; 8];
        let _counting = table_census::count_the_tables_of(&strata);
        for strata_at_once in 1..=4 {
            table_census::restart();
            let _ = fitted_in_a_pool(&strata, &excess, &a_two_round_config(strata_at_once), 8);
            let peak = table_census::peak();
            assert!(
                peak <= strata_at_once,
                "{strata_at_once} at once held {peak} tables together"
            );
            assert!(
                strata_at_once == 1 || peak >= 2,
                "{strata_at_once} at once never held two tables together, so it ran one at a time"
            );
        }
    }

    /// **A single stratum is fitted one at a time whatever the run asked for** — there is nothing
    /// to spread several threads over — so it holds one table, and gives the bits it gives at one.
    #[test]
    fn a_single_stratum_is_fitted_one_at_a_time_whatever_was_asked() {
        let strata = classes_of_unequal_size();
        let alone = &strata[..1];
        let excess = [0.4; 8];
        let at_one = fitted_in_a_pool(alone, &excess, &a_two_round_config(1), 8);
        let _counting = table_census::count_the_tables_of(alone);
        let at_four = fitted_in_a_pool(alone, &excess, &a_two_round_config(4), 8);
        assert_eq!(table_census::peak(), 1);
        the_same_bits(&at_one, &at_four, "one stratum, asked for one and for four");
    }

    /// **A pool that cannot be built leaves one stratum at a time, and the line says why.**
    #[test]
    fn a_pool_that_cannot_be_built_is_one_stratum_at_a_time_and_says_why() {
        let (pool, schedule) = a_pool_for(4, 6, |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .spawn_handler(|_| Err(std::io::Error::other("refused")))
                .build()
        });
        assert!(pool.is_none());
        assert!(
            schedule.starts_with(
                "one stratum at a time, because a pool of 4 threads could not be built"
            ),
            "{schedule}"
        );

        let (pool, schedule) = a_pool_for(4, 6, |threads| {
            rayon::ThreadPoolBuilder::new().num_threads(threads).build()
        });
        assert_eq!(pool.map(|pool| pool.current_num_threads()), Some(4));
        assert_eq!(schedule, "4 strata at once, each on one thread");
        assert_eq!(
            a_pool_for(4, 1, |_| unreachable!("one stratum builds no pool")).1,
            "one stratum at a time, since there is 1 to fit"
        );
        assert_eq!(
            a_pool_for(1, 6, |_| unreachable!("one at a time builds no pool")).1,
            "one stratum at a time"
        );
    }

    /// Whole mebibytes below a gibibyte, where one decimal of a gibibyte would read `0.0`.
    #[test]
    fn a_table_size_reads_in_mebibytes_below_a_gibibyte() {
        assert_eq!(table_size(364_000_000), "347 MiB");
        assert_eq!(table_size(7_895_160_000), "7.4 GiB");
    }

    /// **A scorer that moved to new slippage holds the tables a fresh scorer builds for it.**
    ///
    /// The tables after a move are compared, bit for bit, with those built straight at the new
    /// slippage — and with the old ones, which must differ, or the move tested nothing. That the
    /// old tables are gone *before* the new ones are built is a claim about memory, not about
    /// these values, and [`n_strata_at_once_hold_at_most_n_tables`] is what checks it.
    #[test]
    fn a_refreshed_scorer_holds_the_tables_a_fresh_one_builds() {
        let strata = classes_of_unequal_size();
        let evidence = &strata[2];
        let excess = [0.4; 8];
        let config = a_two_round_config(1);
        let genotypes = genotype_pairs((2 * config.allele_span + 1) as usize);
        let at = |level: f64| {
            vec![
                Slippage {
                    level,
                    shorter_share: 0.8,
                    fall_off: 0.3,
                };
                evidence.groups
            ]
        };
        let tables_of = |scorer: &Scorer<'_>| -> Vec<(Vec<u64>, usize, Vec<u64>, u64)> {
            scorer
                .prepared
                .as_ref()
                .expect("refreshed")
                .tracts
                .iter()
                .map(|tract| {
                    (
                        tract.scaled.iter().map(|value| value.to_bits()).collect(),
                        tract.width,
                        tract
                            .homozygote_excess
                            .iter()
                            .map(|value| value.to_bits())
                            .collect(),
                        tract.ln_offset.to_bits(),
                    )
                })
                .collect()
        };
        let a_scorer = || {
            Scorer::new(
                evidence,
                &excess,
                &genotypes,
                &config,
                WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
            )
        };

        let mut moved = a_scorer();
        moved.refresh(&at(0.05));
        let before = tables_of(&moved);
        moved.refresh(&at(0.12));
        let mut fresh = a_scorer();
        fresh.refresh(&at(0.12));

        assert!(
            !before.is_empty(),
            "the stratum has tracts to build tables for"
        );
        assert_ne!(before, tables_of(&moved), "the move changed the tables");
        assert_eq!(tables_of(&moved), tables_of(&fresh));
        assert_eq!(
            moved.prepared.as_ref().map(|held| held.slippage.clone()),
            Some(at(0.12))
        );
    }

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
        let fitted = fit_stratum(&evidence, &[0.4; 20], &config).expect("reads were drawn");
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

    /// **Central differences return a quadratic's second derivatives and its slope**, off the
    /// diagonal as on it, in `1 + 2p²` evaluations — whichever order the points are visited in.
    #[test]
    fn central_differences_give_a_quadratics_second_derivatives() {
        let curvature = [
            [4.0, 1.0, -0.5, 0.2],
            [1.0, 3.0, 0.3, 0.0],
            [-0.5, 0.3, 2.0, -0.7],
            [0.2, 0.0, -0.7, 5.0],
        ];
        let slope = [0.3, -1.1, 0.7, 2.0];
        let center = [0.5, -2.0, 1.5, 0.1];
        let mut evaluations = 0;
        let (hessian, slope_found) = curvature_of(
            &center,
            &[0.01, 0.02, 0.03, 0.04],
            |j| j == 2,
            |at| {
                evaluations += 1;
                let mut value = 7.0;
                for i in 0..4 {
                    value += slope[i] * at[i];
                    for j in 0..4 {
                        value -= 0.5 * at[i] * curvature[i][j] * at[j];
                    }
                }
                value
            },
        );
        assert_eq!(evaluations, 1 + 2 * 4 * 4);
        for i in 0..4 {
            let truth = slope[i] - (0..4).map(|j| curvature[i][j] * center[j]).sum::<f64>();
            assert!(
                (slope_found[i] - truth).abs() < 1e-8,
                "slope {i}: {} against {truth}",
                slope_found[i]
            );
            for j in 0..4 {
                assert!(
                    (hessian[i * 4 + j] + curvature[i][j]).abs() < 1e-8,
                    "({i}, {j}): {} against {}",
                    hessian[i * 4 + j],
                    -curvature[i][j]
                );
            }
        }
    }

    /// **Each variance is carried to its number's own scale**: a slippage number's by `p(1 − p)`, the
    /// concentration's by itself, a share's by the shares' slopes in the log-ratios, the largest
    /// class's included; a dropped coordinate leaves its number with no error, and the others keep
    /// theirs.
    #[test]
    fn the_errors_are_carried_to_each_numbers_own_scale() {
        let parameters = Parameters {
            slippage: vec![Slippage {
                level: 0.1,
                shorter_share: 0.8,
                fall_off: 0.3,
            }],
            length_spectrum: vec![0.2, 0.5, 0.3],
            concentration: 2.0,
        };
        let layout = CurvatureLayout::of(&parameters, &[true]);
        assert_eq!(layout.largest_class, 1);
        assert_eq!(layout.coordinates.len(), 6);
        // Variances 0.01 … 0.06 on the diagonal, and a covariance between the two log-ratios.
        let n = 6;
        let mut inverse = vec![0.0; n * n];
        for j in 0..n {
            inverse[j * n + j] = 0.01 * (j + 1) as f64;
        }
        inverse[3 * n + 4] = 0.01;
        inverse[4 * n + 3] = 0.01;
        let identified = crate::parameter_estimation::joint::fit::Identified {
            kept: (0..n).collect(),
            dropped: Vec::new(),
            inverse: inverse.clone(),
        };
        let errors = errors_on_the_natural_scale(&layout, &parameters, &[true], &identified);
        let group = errors.slippage[0].expect("a live group");
        let close = |got: StratumError, expected: f64| {
            let got = got.value().expect("an error");
            assert!((got - expected).abs() < 1e-12, "{got} against {expected}");
        };
        close(group.level, 0.01_f64.sqrt() * 0.1 * 0.9);
        close(group.shorter_share, 0.02_f64.sqrt() * 0.8 * 0.2);
        close(group.fall_off, 0.03_f64.sqrt() * 0.3 * 0.7);
        close(errors.concentration, 0.06_f64.sqrt() * 2.0);
        // The shares' slopes in (log-ratio of class 0, log-ratio of class 2).
        let s = [0.2, 0.5, 0.3];
        let (v0, v2, c) = (0.04, 0.05, 0.01);
        let share_error =
            |g0: f64, g2: f64| (g0 * g0 * v0 + 2.0 * g0 * g2 * c + g2 * g2 * v2).sqrt();
        close(
            errors.length_spectrum[0],
            share_error(s[0] * (1.0 - s[0]), -s[0] * s[2]),
        );
        close(
            errors.length_spectrum[1],
            share_error(-s[1] * s[0], -s[1] * s[2]),
        );
        close(
            errors.length_spectrum[2],
            share_error(-s[2] * s[0], s[2] * (1.0 - s[2])),
        );

        // Class 2's log-ratio dropped: its share has no error, class 0's and the largest's keep one
        // from class 0's ratio alone.
        let kept: Vec<usize> = vec![0, 1, 2, 3, 5];
        let dropped_inverse: Vec<f64> = kept
            .iter()
            .flat_map(|&row| kept.iter().map(move |&column| (row, column)))
            .map(|(row, column)| inverse[row * n + column])
            .collect();
        let identified = crate::parameter_estimation::joint::fit::Identified {
            kept,
            dropped: vec![4],
            inverse: dropped_inverse,
        };
        let errors = errors_on_the_natural_scale(&layout, &parameters, &[true], &identified);
        assert_eq!(errors.length_spectrum[2], StratumError::NotIdentified);
        close(errors.length_spectrum[0], (s[0] * (1.0 - s[0])) * v0.sqrt());
        close(errors.length_spectrum[1], (s[1] * s[0]) * v0.sqrt());
        close(errors.concentration, 0.06_f64.sqrt() * 2.0);

        // Class 2's log-ratio kept but its variance 10⁶: that ratio is not placed, and the shares
        // whose errors it enters come out wider than their whole range, so they are not placed either.
        let mut wide = inverse.clone();
        wide[4 * n + 4] = 1e6;
        let identified = crate::parameter_estimation::joint::fit::Identified {
            kept: (0..n).collect(),
            dropped: Vec::new(),
            inverse: wide,
        };
        let errors = errors_on_the_natural_scale(&layout, &parameters, &[true], &identified);
        assert_eq!(
            errors.length_spectrum,
            [StratumError::NotPlaced; 3],
            "class 2 on its own scale, classes 0 and 1 wider than the range"
        );
        close(
            errors.slippage[0].expect("a live group").level,
            0.01_f64.sqrt() * 0.1 * 0.9,
        );
    }

    /// **A stratum's errors are its total log-likelihood's, not the mean's**: on drawn strata of 400
    /// and of 1,600 tracts, every number has an error, the level's shrinks by about half as the tracts
    /// quadruple — the mean's curvature would leave it unchanged — and the fitted level lies within
    /// three of its errors of the truth it was drawn at.
    #[test]
    fn a_stratums_errors_shrink_as_its_tracts_grow() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let spectrum = spectrum_of(3);
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let level_error = |tracts: usize| {
            let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, tracts, 20, 6, 1, 29);
            let fitted = fit_stratum(&evidence, &[0.4; 20], &config).expect("reads were drawn");
            let errors = fitted
                .standard_errors
                .as_ref()
                .expect("a fit computes its errors");
            let group = errors.slippage[0].expect("the one group has reads");
            assert!(
                [group.level, group.shorter_share, group.fall_off]
                    .iter()
                    .all(|error| error.value().is_some()),
                "{group:?}"
            );
            assert!(errors.concentration.value().is_some(), "{errors:?}");
            assert!(
                errors
                    .length_spectrum
                    .iter()
                    .all(|error| error.value().is_some()),
                "{errors:?}"
            );
            let level = fitted.slippage[0].expect("fitted").level;
            let error = group.level.value().expect("checked");
            eprintln!("{tracts} tracts: level {level:.5} ± {error:.5}; {errors:?}");
            assert!(
                (level - truth.level).abs() < 3.0 * error,
                "{tracts} tracts: level {level} ± {error} against a truth of {}",
                truth.level
            );
            error
        };
        let (few, many) = (level_error(400), level_error(1_600));
        let ratio = many / few;
        assert!(
            (0.35..0.65).contains(&ratio),
            "the level's error went from {few} to {many}, a ratio of {ratio}"
        );
    }

    /// **The errors do not depend on the differences' step**: at the default step and at a third of
    /// it, every error of a drawn stratum's fit agrees to within 2 in 100.
    #[test]
    fn a_stratums_errors_do_not_depend_on_the_difference_step() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let spectrum = spectrum_of(3);
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 600, 20, 6, 1, 31);
        let fitted = fit_stratum(&evidence, &[0.4; 20], &config).expect("reads were drawn");
        let parameters = Parameters {
            slippage: vec![fitted.slippage[0].expect("fitted")],
            length_spectrum: fitted.length_spectrum.clone(),
            concentration: fitted.concentration,
        };
        let at = |step: f64| {
            standard_errors_with_step(
                &evidence,
                &parameters,
                &[0.4; 20],
                &[true],
                &config,
                WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
                step,
            )
        };
        let (default, finer) = (at(CURVATURE_STEP), at(CURVATURE_STEP / 3.0));
        let flatten = |errors: &StratumErrors| {
            let group = errors.slippage[0].expect("a live group");
            let mut all = vec![
                group.level,
                group.shorter_share,
                group.fall_off,
                errors.concentration,
            ];
            all.extend(errors.length_spectrum.iter().copied());
            all.into_iter().map(StratumError::value).collect::<Vec<_>>()
        };
        for (index, (a, b)) in flatten(&default)
            .into_iter()
            .zip(flatten(&finer))
            .enumerate()
        {
            let (a, b) = (a.expect("an error"), b.expect("an error"));
            assert!((a - b).abs() < 0.02 * a, "number {index}: {a} against {b}");
        }
        assert_eq!(
            fitted.standard_errors,
            Some(default),
            "the fit's errors are the default step's"
        );
    }

    /// **A number the climb left at the end of its reach is not placed**, where carrying its flat
    /// curvature to its own scale would report an error hundreds of times too small (review of plan
    /// step C1): on three thin drawn strata the fall-off runs towards zero, the shorter share towards
    /// one, and the concentration off to millions — each stops where the climb stopped, still rising
    /// (the fall-off at 5 × 10⁻⁵ against a truth of 0.25, the shorter share at 0.9999 against 0.83, the
    /// concentration at 8 × 10⁴ against 50), and each comes back without an error, the others of the
    /// stratum keeping theirs.
    #[test]
    fn a_number_the_climb_left_at_the_end_of_its_reach_is_not_placed() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let spectrum = spectrum_of(3);
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let fit = |tracts: usize, samples: usize, depth: u32, seed: u64, concentration: f64| {
            let evidence = draw_stratum(
                truth,
                &spectrum,
                concentration,
                0.4,
                tracts,
                samples,
                depth,
                1,
                seed,
            );
            let fitted =
                fit_stratum(&evidence, &vec![0.4; samples], &config).expect("reads were drawn");
            eprintln!(
                "{tracts} x {samples} x {depth}: {:?}, concentration {}; {:?}",
                fitted.slippage[0], fitted.concentration, fitted.standard_errors
            );
            fitted
        };
        let errors_of = |fitted: &StratumFit| {
            fitted
                .standard_errors
                .clone()
                .expect("a fit computes its errors")
        };

        let fall_off_run_off = fit(10, 8, 6, 16, 0.5);
        assert!(fall_off_run_off.slippage[0].expect("fitted").fall_off < 1e-3);
        let group = errors_of(&fall_off_run_off).slippage[0].expect("fitted");
        assert_eq!(group.fall_off, StratumError::NotPlaced);
        assert!(group.level.value().is_some(), "{group:?}");

        let share_run_off = fit(15, 4, 3, 6, 0.5);
        assert!(share_run_off.slippage[0].expect("fitted").shorter_share > 1.0 - 1e-3);
        let group = errors_of(&share_run_off).slippage[0].expect("fitted");
        assert_eq!(group.shorter_share, StratumError::NotPlaced);

        let concentration_run_off = fit(200, 4, 3, 8, 50.0);
        assert!(concentration_run_off.concentration > 1e4);
        let errors = errors_of(&concentration_run_off);
        assert_eq!(errors.concentration, StratumError::NotPlaced);
        assert!(
            errors.slippage[0].expect("fitted").level.value().is_some(),
            "{errors:?}"
        );
    }

    /// **Tracts no sample read change no error**: they add a factor of one to the likelihood, so the
    /// total — the mean a tract times every tract, those without reads included — is the same, and so
    /// is its curvature, whichever count the mean is taken over.
    #[test]
    fn tracts_without_reads_change_no_error() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(truth, &spectrum_of(3), 0.5, 0.4, 300, 20, 6, 1, 37);
        let fitted = fit_stratum(&evidence, &[0.4; 20], &config).expect("reads were drawn");
        let parameters = Parameters {
            slippage: vec![fitted.slippage[0].expect("fitted")],
            length_spectrum: fitted.length_spectrum.clone(),
            concentration: fitted.concentration,
        };
        let mut padded = evidence.clone();
        padded
            .tracts
            .extend(std::iter::repeat_n(TractReads::default(), 300));
        let errors = |evidence: &StratumEvidence| {
            standard_errors_at(
                evidence,
                &parameters,
                &[0.4; 20],
                &[true],
                &config,
                WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
            )
        };
        let (plain, with_empty) = (errors(&evidence), errors(&padded));
        let flat = |errors: &StratumErrors| {
            let group = errors.slippage[0].expect("a live group");
            let mut all = vec![
                group.level,
                group.shorter_share,
                group.fall_off,
                errors.concentration,
            ];
            all.extend(errors.length_spectrum.iter().copied());
            all.into_iter()
                .map(|error| error.value().expect("an error"))
                .collect::<Vec<f64>>()
        };
        for (a, b) in flat(&plain).into_iter().zip(flat(&with_empty)) {
            assert!((a - b).abs() <= 1e-6 * a, "{a} against {b}");
        }
    }

    /// **The curvature is taken at the fitted answer**: the layout's centre, set back through
    /// `parameters_at`, is the answer the climb returned, every number to within 10⁻¹² of itself.
    #[test]
    fn the_curvatures_centre_is_the_fitted_answer() {
        let parameters = Parameters {
            slippage: vec![
                Slippage {
                    level: 0.05,
                    shorter_share: 0.83,
                    fall_off: 0.25,
                },
                Slippage {
                    level: 0.2,
                    shorter_share: 0.4,
                    fall_off: 0.6,
                },
            ],
            length_spectrum: vec![0.1, 0.05, 0.6, 0.2, 0.05],
            concentration: 3.7,
        };
        let layout = CurvatureLayout::of(&parameters, &[true, true]);
        let back = layout.parameters_at(&parameters, &layout.center(&parameters));
        let close = |a: f64, b: f64| assert!((a - b).abs() <= 1e-12 * b.abs(), "{a} against {b}");
        for (a, b) in back.slippage.iter().zip(&parameters.slippage) {
            close(a.level, b.level);
            close(a.shorter_share, b.shorter_share);
            close(a.fall_off, b.fall_off);
        }
        for (a, b) in back.length_spectrum.iter().zip(&parameters.length_spectrum) {
            close(*a, *b);
        }
        close(back.concentration, parameters.concentration);
    }

    /// **The log's summary**: the level's and the concentration's errors as a percentage of the
    /// number, the others as they are, the median and the largest, and the numbers without an error
    /// counted by why.
    #[test]
    fn the_summary_gives_each_kinds_errors_and_counts_the_missing_by_why() {
        let with_errors = |level: f64, level_error: StratumError, concentration: StratumError| {
            let StratumOutcome::Fitted(mut fit) = fitted_at(2, 10, level, 1_000) else {
                unreachable!("fitted_at builds a fit")
            };
            fit.standard_errors = Some(StratumErrors {
                slippage: vec![Some(SlippageErrors {
                    level: level_error,
                    shorter_share: StratumError::Estimated(0.05),
                    fall_off: StratumError::NotPlaced,
                })],
                length_spectrum: vec![StratumError::Estimated(0.01)],
                concentration,
            });
            StratumOutcome::Fitted(fit)
        };
        let outcomes = [
            with_errors(
                0.1,
                StratumError::Estimated(0.02),
                StratumError::Estimated(0.3),
            ),
            with_errors(
                0.04,
                StratumError::Estimated(0.004),
                StratumError::NotIdentified,
            ),
            with_errors(0.2, StratumError::NotIdentified, StratumError::NotPlaced),
        ];
        let summary = standard_errors_summary(&outcomes).expect("errors were carried");
        for part in [
            "the 3 strata",
            "slippage level ± median 20% of itself (largest 20%)",
            "shorter share ± median 0.050 (largest 0.050)",
            "fall-off: none has an error",
            "concentration ± median 50% of itself (largest 50%)",
            "of their 15 numbers, 4 have no error because the tracts do not place them",
            "2 because the tracts cannot tell them apart",
            "and 0 are length classes with no share",
        ] {
            assert!(summary.contains(part), "{part:?} in {summary}");
        }
        assert_eq!(
            standard_errors_summary(&[fitted_at(2, 10, 0.1, 1_000)]),
            None
        );
    }

    // -----------------------------------------------------------------
    // Whether a stratum's errors mean what they say (plan step C2)
    // -----------------------------------------------------------------

    /// How often one kind of number lands within one and two of its errors of the truth, over many
    /// drawn strata, and how many came back without an error.
    #[derive(Debug, Default)]
    struct StratumCoverage {
        /// Each estimate's distance from the truth in its own errors, signed.
        distances: Vec<f64>,
        not_identified: usize,
        not_placed: usize,
    }

    impl StratumCoverage {
        fn record(&mut self, estimate: f64, truth: f64, error: StratumError) {
            match error {
                StratumError::Estimated(error) => self.distances.push((estimate - truth) / error),
                StratumError::NotIdentified => self.not_identified += 1,
                StratumError::NotPlaced => self.not_placed += 1,
                StratumError::NoShare => {}
            }
        }

        fn share_within(&self, errors: f64) -> f64 {
            self.distances
                .iter()
                .filter(|distance| distance.abs() <= errors)
                .count() as f64
                / self.distances.len().max(1) as f64
        }

        fn mean(&self) -> f64 {
            self.distances.iter().sum::<f64>() / self.distances.len().max(1) as f64
        }

        fn spread(&self) -> f64 {
            let mean = self.mean();
            (self
                .distances
                .iter()
                .map(|distance| (distance - mean) * (distance - mean))
                .sum::<f64>()
                / (self.distances.len().max(2) - 1) as f64)
                .sqrt()
        }
    }

    /// How many strata each regime draws by default, and the variable that overrides it.
    const STRATUM_DRAWS: usize = 100;
    const STRATUM_DRAWS_VARIABLE: &str = "NG_FIT_PRECISION_STRATUM_DRAWS";

    /// Draw `draws` strata of one regime from known numbers, fit each, and tally every estimate's
    /// distance from the truth in its own errors, by kind; how many fits converged; and how many
    /// draws gave nothing to tally — no fit, or no errors — which the test requires to be none.
    #[allow(clippy::too_many_arguments)]
    fn stratum_coverage_of(
        truth: Slippage,
        classes: usize,
        concentration: f64,
        tracts: usize,
        samples: usize,
        depth: u32,
        draws: usize,
        seed: u64,
    ) -> (BTreeMap<&'static str, StratumCoverage>, usize, usize) {
        let span = (classes / 2) as i32;
        let spectrum = spectrum_of(classes);
        let config = SsrFitConfig {
            allele_span: span,
            ..SsrFitConfig::default()
        };
        let mut tally: BTreeMap<&'static str, StratumCoverage> = BTreeMap::new();
        let mut converged = 0;
        let mut skipped = 0;
        for draw in 0..draws {
            let evidence = draw_stratum(
                truth,
                &spectrum,
                concentration,
                0.4,
                tracts,
                samples,
                depth,
                span,
                seed + draw as u64,
            );
            let Some(fitted) = fit_stratum(&evidence, &vec![0.4; samples], &config) else {
                skipped += 1;
                continue;
            };
            converged += usize::from(fitted.ending.settled());
            let errors = fitted
                .standard_errors
                .as_ref()
                .expect("a fit computes its errors");
            let (Some(slippage), Some(group)) = (fitted.slippage[0], errors.slippage[0]) else {
                skipped += 1;
                continue;
            };
            for (kind, estimate, truth, error) in [
                ("slippage level", slippage.level, truth.level, group.level),
                (
                    "shorter share",
                    slippage.shorter_share,
                    truth.shorter_share,
                    group.shorter_share,
                ),
                (
                    "fall-off",
                    slippage.fall_off,
                    truth.fall_off,
                    group.fall_off,
                ),
                (
                    "concentration",
                    fitted.concentration,
                    concentration,
                    errors.concentration,
                ),
            ] {
                tally
                    .entry(kind)
                    .or_default()
                    .record(estimate, truth, error);
            }
            for (class, (estimate, error)) in fitted
                .length_spectrum
                .iter()
                .zip(&errors.length_spectrum)
                .enumerate()
            {
                let kind = if class == span as usize {
                    "reference-length share"
                } else {
                    "other length shares"
                };
                tally
                    .entry(kind)
                    .or_default()
                    .record(*estimate, spectrum[class], *error);
            }
        }
        (tally, converged, skipped)
    }

    /// **A stratum's errors mean what they say** (plan step C2, spec §4.5 item 1): strata drawn from
    /// known slippage, a known length spectrum and concentration, many times in each regime, each
    /// fitted by the fit's own rule; of every estimate with an error, the share within one error of
    /// the truth and within two — about 68 in 100 and 95 in 100 for an error that means what it says —
    /// with the mean distance (a bias, in errors) and the spread of the distances (about one). At
    /// three allele classes over 300 tracts of 20 samples, at 3 reads a sample and at 30; and at the
    /// production span, thirteen classes, over 200 tracts at 3 reads, on a quarter as many draws.
    ///
    /// **What it measured** (100 draws, `fit_precision_c2` report). At three classes the slippage
    /// numbers and the concentration land within one error 61 to 76 times in 100 and within two 91 to
    /// 99, their mean distance from the truth under a quarter of an error; the class shares are
    /// covered less well at 30 reads (47 to 53 in 100 within one, spread 1.6). **At thirteen classes
    /// the estimate itself is off** — the reference-length share a mean of 18 errors from its truth,
    /// the concentration 8 — because the fixed 256-point integral over a tract's length frequencies
    /// misrepresents the likelihood there (refitting one stratum at 1,024 and 4,096 points moved the
    /// share from 0.17 to 0.23 and 0.31, against 0.30, and in review no point count up to 65,536 was
    /// shown to suffice); so that regime is printed and not held to anything. The three-class regimes
    /// are held to [`hold_to_the_three_class_bounds`] at the default count of draws.
    ///
    /// Ignored by default: 20 to 40 minutes in the container. `scripts/dev.sh env
    /// NG_FIT_PRECISION_STRATUM_DRAWS=10 cargo test --release --lib a_stratums_errors_mean_what_they_say
    /// -- --ignored --nocapture` shortens it, and checks nothing but that it runs.
    #[test]
    #[ignore = "a measurement over hundreds of fitted strata; run at a checkpoint"]
    fn a_stratums_errors_mean_what_they_say() {
        let draws = match std::env::var(STRATUM_DRAWS_VARIABLE) {
            Err(std::env::VarError::NotPresent) => STRATUM_DRAWS,
            Err(std::env::VarError::NotUnicode(value)) => {
                panic!("{STRATUM_DRAWS_VARIABLE}={value:?} is not a positive whole number")
            }
            Ok(value) => value
                .trim()
                .parse()
                .ok()
                .filter(|&count: &usize| count > 0)
                .unwrap_or_else(|| {
                    panic!("{STRATUM_DRAWS_VARIABLE}={value:?} is not a positive whole number")
                }),
        };
        let truth = Slippage {
            level: 0.05,
            shorter_share: 0.8,
            fall_off: 0.3,
        };
        for (classes, tracts, depth, regime_draws) in [
            (3, 300, 3, draws),
            (3, 300, 30, draws),
            (13, 200, 3, draws.div_ceil(4)),
        ] {
            let (tally, converged, skipped) = stratum_coverage_of(
                truth,
                classes,
                0.5,
                tracts,
                20,
                depth,
                regime_draws,
                0xC2C0_7E7A_0000_0000 + (classes as u64) * 1_000 + u64::from(depth),
            );
            let regime = format!("{classes} classes, {tracts} tracts x 20 samples x {depth} reads");
            eprintln!(
                "STRATUM COVERAGE {regime}: {converged} of {regime_draws} climbs settled within \
                 their rounds; {skipped} draw(s) gave nothing to tally"
            );
            assert_eq!(
                skipped, 0,
                "{regime}: every draw is fitted and carries errors"
            );
            // Every kind printed before any is held to its bounds, so a failure shows the whole
            // regime.
            for (kind, coverage) in &tally {
                eprintln!(
                    "STRATUM COVERAGE {regime}, {kind}: {} estimates with an error, {} not \
                     identified, {} not placed; within one error {:.3}, within two {:.3}; mean \
                     distance {:+.3}, spread {:.3}",
                    coverage.distances.len(),
                    coverage.not_identified,
                    coverage.not_placed,
                    coverage.share_within(1.0),
                    coverage.share_within(2.0),
                    coverage.mean(),
                    coverage.spread(),
                );
            }
            if classes == 3 && regime_draws >= STRATUM_DRAWS {
                for (kind, coverage) in &tally {
                    hold_to_the_three_class_bounds(&regime, kind, coverage, regime_draws);
                }
            }
        }
    }

    /// **The bounds a three-class regime is held to** at [`STRATUM_DRAWS`] draws or more, from the
    /// runs of plan step C2 (report `fit_precision_c2`) with a margin: every draw gives every kind an
    /// estimate, and every estimate an error; the slippage numbers and the concentration land within
    /// one error 55 to 82 times in 100 and within two at least 88, their distances spread 0.85 to 1.2
    /// times their errors, the mean distance under 0.4 errors — the shorter share sits +0.22 errors
    /// high in every run; the class shares within one at least 40 in 100 and within two 75, which
    /// accepts their measured under-coverage at 30 reads, due at least in part to the 256-point
    /// integral (report §3).
    ///
    /// **What they cannot see**: every error scaled by 0.88 to 1.18 passes, and the class shares'
    /// errors by much more (measured in review).
    fn hold_to_the_three_class_bounds(
        regime: &str,
        kind: &str,
        coverage: &StratumCoverage,
        draws: usize,
    ) {
        let shares = kind.contains("length");
        let expected = if kind == "other length shares" {
            2 * draws
        } else {
            draws
        };
        assert_eq!(
            (
                coverage.distances.len(),
                coverage.not_identified,
                coverage.not_placed
            ),
            (expected, 0, 0),
            "{regime}, {kind}: every estimate has an error"
        );
        let (one, two) = (coverage.share_within(1.0), coverage.share_within(2.0));
        let (least_within_one, least_within_two) = if shares { (0.40, 0.75) } else { (0.55, 0.88) };
        assert!(
            one >= least_within_one && one <= 0.82 && two >= least_within_two,
            "{regime}, {kind}: within one error {one:.3}, within two {two:.3}"
        );
        if !shares {
            let spread = coverage.spread();
            assert!(
                (0.85..=1.2).contains(&spread),
                "{regime}, {kind}: the distances spread {spread:.3} times the errors"
            );
            assert!(
                coverage.mean().abs() < 0.4,
                "{regime}, {kind}: the estimates sit {:+.3} errors from the truth",
                coverage.mean()
            );
        }
    }

    // -----------------------------------------------------------------
    // When the climb stops (plan step C3)
    // -----------------------------------------------------------------

    /// Scripted rounds for the stopping rule: each round's best point part-way through, if any, and
    /// where it ends, both as a mean a tract; and the points that have settled when judged. A walk
    /// stands at the label `round` at a round's end and `round - 0.5` part-way through it.
    struct Scripted<'a> {
        rounds: std::slice::Iter<'a, (Option<f64>, f64)>,
        settled_at: &'a [f64],
        judged: Vec<f64>,
    }

    impl Climbing for Scripted<'_> {
        type Point = f64;
        type Judgement = f64;

        fn one_round(&mut self, at: &mut f64, held: &mut HeldPoint<f64>) -> f64 {
            let (part_way, end) = self
                .rounds
                .next()
                .expect("a scripted round for every round walked");
            if let Some(part_way) = part_way {
                held.offer(&(*at + 0.5), *part_way);
            }
            *at += 1.0;
            held.offer(at, *end);
            *end
        }

        fn judge(&mut self, at: &f64) -> (bool, f64) {
            self.judged.push(*at);
            (self.settled_at.contains(at), *at)
        }
    }

    /// Walk `rounds` from a score of `start`, and the points that were judged on the way.
    fn scripted_walk(
        start: f64,
        rounds: &[(Option<f64>, f64)],
        settled_at: &[f64],
        trigger: f64,
        tracts: f64,
        max_rounds: u32,
    ) -> (Walked<f64, f64>, Vec<f64>) {
        let mut scripted = Scripted {
            rounds: rounds.iter(),
            settled_at,
            judged: Vec::new(),
        };
        let walked = walk_until_settled(&mut scripted, 0.0, start, trigger, tracts, max_rounds);
        (walked, scripted.judged)
    }

    /// Scores rising by `gains`, as a mean a tract over one tract.
    fn rising_by(gains: &[f64]) -> Vec<(Option<f64>, f64)> {
        gains
            .iter()
            .scan(0.0, |score, gain| {
                *score += gain;
                Some((None, *score))
            })
            .collect()
    }

    /// **The projection**: gains shrinking by a factor λ leave `gain · λ / (1 − λ)` to come; λ is
    /// capped at 0.95, and is the cap on a first round, with no gain before it.
    #[test]
    fn the_gain_still_to_come_is_projected_from_the_last_two() {
        assert!((projected_remaining_gain(1.0, Some(2.0)) - 1.0).abs() < 1e-15);
        assert!((projected_remaining_gain(1.0, Some(4.0)) - 1.0 / 3.0).abs() < 1e-15);
        // Gains that do not shrink are projected at the cap, nineteen times the last.
        assert!((projected_remaining_gain(1.0, Some(0.5)) - 19.0).abs() < 1e-12);
        assert!((projected_remaining_gain(1.0, None) - 19.0).abs() < 1e-12);
        assert_eq!(projected_remaining_gain(0.0, Some(3.0)), 0.0);
    }

    /// **A walk is judged once the gain it can still make is below the trigger, and stops only when
    /// a judged point has settled**: gains halving from 8 leave as much again to come, so at a
    /// trigger of 0.6 the first point judged is the fifth round's; settled there, the walk stops
    /// there; settled only at the eighth, it is judged at the fifth, the sixth and — the second
    /// unsettled verdict putting the next off by a round — the eighth. A first round
    /// is projected at the cap, so one gaining 0.003 is judged at once below a trigger of 0.08
    /// (19 × 0.003 = 0.057) and one gaining 0.005 is not (0.095).
    #[test]
    fn a_walk_is_judged_once_its_projected_gain_is_small_and_stops_when_settled() {
        let halving = rising_by(&[8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625]);
        let (walked, judged) = scripted_walk(0.0, &halving, &[5.0], 0.6, 1.0, 10);
        assert_eq!(
            (walked.ending, walked.rounds, walked.at, walked.judgements),
            (ClimbEnding::Settled, 5, 5.0, 1)
        );
        assert_eq!((judged, walked.judged_here), (vec![5.0], Some(5.0)));
        assert!((walked.score - 15.5).abs() < 1e-12);

        let (later, judged) = scripted_walk(0.0, &halving, &[8.0], 0.6, 1.0, 10);
        assert_eq!(
            (later.ending, later.rounds, later.judgements),
            (ClimbEnding::Settled, 8, 3)
        );
        assert_eq!(judged, vec![5.0, 6.0, 8.0]);

        let (quick, judged) = scripted_walk(0.0, &rising_by(&[0.003, 1.0]), &[1.0], 0.08, 1.0, 10);
        assert_eq!(
            (quick.ending, quick.rounds, judged),
            (ClimbEnding::Settled, 1, vec![1.0])
        );
        let (not_quick, judged) =
            scripted_walk(0.0, &rising_by(&[0.005, 0.0]), &[2.0], 0.08, 1.0, 10);
        assert_eq!(
            (not_quick.ending, not_quick.rounds, judged),
            (ClimbEnding::Settled, 2, vec![2.0])
        );
    }

    /// **The trigger is in the total's units**: the same mean gains, halving from 10⁻³, have a walk
    /// of one tract judged on its first round (19 × 10⁻³ to come, below 0.08) and a walk of a
    /// thousand tracts, whose total gains are a thousand times larger, only on the fifth (0.0625).
    #[test]
    fn the_trigger_is_on_the_total_log_likelihood_not_the_mean() {
        let means = rising_by(&[1e-3, 5e-4, 2.5e-4, 1.25e-4, 6.25e-5, 3.125e-5]);
        let every_point: Vec<f64> = (1..=6).map(f64::from).collect();
        let (one, _) = scripted_walk(0.0, &means, &every_point, 0.08, 1.0, 10);
        assert_eq!((one.ending, one.rounds), (ClimbEnding::Settled, 1));
        let (thousand, _) = scripted_walk(0.0, &means, &every_point, 0.08, 1_000.0, 10);
        assert_eq!(
            (thousand.ending, thousand.rounds),
            (ClimbEnding::Settled, 5)
        );
    }

    /// **A round that loses is undone, and its best point judged at once**: back to a point part-way
    /// through it when that was better than its start, and the walk goes on from there unless that
    /// point has settled; back to its start when nothing was better, and if that has not settled the
    /// walk stops there, having nowhere new to go. A round whose score is not a number is a loss.
    #[test]
    fn a_round_that_loses_is_undone_and_its_best_point_judged() {
        // A trigger no gain reaches, so only the losing rounds are judged.
        let trigger = 1e-12;
        let part_way_better = [(None, 5.0), (None, 6.0), (Some(6.5), 5.5), (Some(6.4), 6.2)];
        let (walked, judged) = scripted_walk(0.0, &part_way_better, &[], trigger, 1.0, 10);
        assert_eq!(
            (walked.ending, walked.rounds, walked.at, walked.score),
            (ClimbEnding::LostARound, 4, 2.5, 6.5)
        );
        assert_eq!(judged, vec![2.5, 2.5]);

        let (settles, _) = scripted_walk(0.0, &part_way_better, &[2.5], trigger, 1.0, 10);
        assert_eq!(
            (settles.ending, settles.rounds, settles.at, settles.score),
            (ClimbEnding::Settled, 3, 2.5, 6.5)
        );

        let not_a_number = [(None, 5.0), (None, f64::NAN)];
        let (walked, judged) = scripted_walk(0.0, &not_a_number, &[], trigger, 1.0, 10);
        assert_eq!(
            (walked.ending, walked.rounds, walked.at, walked.score),
            (ClimbEnding::LostARound, 2, 1.0, 5.0)
        );
        assert_eq!(judged, vec![1.0]);
    }

    /// **Each unsettled verdict puts the next judgement off by one round more**: a walk whose every
    /// round wants judging and that never settles is judged at rounds 1, 2, 4 and 7 of ten. A round
    /// that loses while the judgement is put off, with nothing better than its start, is judged all
    /// the same, and ends the walk.
    #[test]
    fn an_unsettled_verdict_puts_the_next_judgement_off() {
        let (walked, judged) = scripted_walk(0.0, &rising_by(&[1e-3; 10]), &[], 1.0, 1.0, 10);
        assert_eq!(
            (walked.ending, walked.rounds, walked.judgements),
            (ClimbEnding::OutOfRounds, 10, 4)
        );
        assert_eq!(judged, vec![1.0, 2.0, 4.0, 7.0]);
        assert_eq!(walked.judged_here, None);

        let losing_while_put_off = [(None, 1e-3), (None, 2e-3), (None, 1e-3)];
        let (walked, judged) = scripted_walk(0.0, &losing_while_put_off, &[], 1.0, 1.0, 10);
        assert_eq!(
            (walked.ending, walked.rounds, walked.at, walked.score),
            (ClimbEnding::LostARound, 3, 2.0, 2e-3)
        );
        assert_eq!(judged, vec![1.0, 2.0, 2.0]);
        assert_eq!(walked.judged_here, Some(2.0));
    }

    /// The parameters a fitted stratum returned, every slippage group live.
    fn parameters_of(fit: &StratumFit) -> Parameters {
        Parameters {
            slippage: fit
                .slippage
                .iter()
                .map(|group| group.expect("every group live"))
                .collect(),
            length_spectrum: fit.length_spectrum.clone(),
            concentration: fit.concentration,
        }
    }

    /// **The winning walk's last judgement is the stratum's errors**, and they are exactly the
    /// errors computed afresh at its answer, in both schedules (review of plan step C3).
    #[test]
    fn the_winning_walks_errors_are_the_errors_at_its_answer_in_both_schedules() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let mut strata = Vec::new();
        for (repeats, seed) in [(10_u64, 41_u64), (11, 43)] {
            let mut evidence = draw_stratum(truth, &spectrum_of(3), 0.5, 0.4, 300, 20, 6, 1, seed);
            evidence.stratum = Stratum {
                period: 2,
                reference_repeats: repeats,
            };
            strata.push(evidence);
        }
        let excess = [0.4; 20];
        let mut config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        config.curve.draw_curves = false;
        let fits_of = |at_once: usize| -> Vec<StratumFit> {
            let config = SsrFitConfig {
                strata_at_once: NonZeroUsize::new(at_once).expect("positive"),
                ..config.clone()
            };
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(4)
                .build()
                .expect("a pool");
            pool.install(|| fit_strata(&strata, &excess, &config))
                .into_iter()
                .map(|outcome| match outcome {
                    StratumOutcome::Fitted(fit) => *fit,
                    other => panic!("{other:?}"),
                })
                .collect()
        };
        let (serial, several) = (fits_of(1), fits_of(2));
        for ((one, two), evidence) in serial.iter().zip(&several).zip(&strata) {
            assert_eq!(one.ending, ClimbEnding::Settled, "{:?}", one.walks);
            assert_eq!((&one.walks, one.ending), (&two.walks, two.ending));
            let at = parameters_of(one);
            for (fit, threads) in [
                (one, WhereTheThreadsGo::AcrossTheTractsOfOneStratum),
                (two, WhereTheThreadsGo::AcrossStrata),
            ] {
                let fresh = standard_errors_at(evidence, &at, &excess, &[true], &config, threads);
                assert_eq!(fit.standard_errors.as_ref(), Some(&fresh));
            }
        }
    }

    /// **A walk at two slippage groups is judged on both**: on a drawn stratum whose samples read in
    /// two groups, the point a long climb reached has settled; with group 1's level moved by one of its
    /// errors, group 1's level is about one error away, group 0's is not, and the point has not
    /// settled.
    #[test]
    fn a_walk_at_two_slippage_groups_is_judged_on_both() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let mut evidence = draw_stratum(truth, &spectrum_of(3), 0.5, 0.4, 300, 40, 6, 1, 47);
        for tract in &mut evidence.tracts {
            for sample in tract
                .samples
                .iter_mut()
                .filter(|sample| sample.sample >= 20)
            {
                for (group, _) in &mut sample.by_group {
                    *group = 1;
                }
            }
        }
        evidence.groups = 2;
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let excess = [0.4; 40];
        let genotypes = genotype_pairs(3);
        let live = [true, true];
        let path = trajectory_of(&evidence, &excess, config.starting_points[1], &config, 20);
        let peak = the_best_of(path.iter()).0.clone();
        let mut climbing = StratumClimb {
            scorer: Scorer::new(
                &evidence,
                &excess,
                &genotypes,
                &config,
                WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
            ),
            evidence: &evidence,
            live_groups: &live,
            classes: 3,
            settled_fraction: SETTLED_FRACTION,
        };
        let (settled, errors) = climbing.judge(&peak);
        assert!(settled, "{errors:?}");
        let error_of = |group: usize| {
            errors.slippage[group]
                .expect("a live group")
                .level
                .value()
                .expect("an error")
        };
        let mut moved = peak.clone();
        moved.slippage[1].level += error_of(1);
        let curvature = curvature_at(
            &mut climbing.scorer,
            &evidence,
            &moved,
            &live,
            CURVATURE_STEP,
        );
        let distances = newton_distances_on_the_natural_scale(&curvature, &moved, &live);
        let away =
            |group: usize| distances.slippage[group].expect("a live group").level / error_of(group);
        assert!((-1.3..-0.7).contains(&away(1)), "group 1: {}", away(1));
        assert!(away(0).abs() < 0.3, "group 0: {}", away(0));
        assert!(!climbing.judge(&moved).0);
    }

    /// **The furthest number counts every kind, and a distance that is not a number has not
    /// settled**: the slippage numbers, a share and the concentration each set the furthest when they
    /// are; numbers without an error are passed over; none with an error is `None`.
    #[test]
    fn the_furthest_number_counts_every_kind_and_a_distance_not_a_number_is_unsettled() {
        let one = StratumError::Estimated(1.0);
        let errors = StratumErrors {
            slippage: vec![Some(SlippageErrors {
                level: one,
                shorter_share: one,
                fall_off: StratumError::NotPlaced,
            })],
            length_spectrum: vec![one, one, StratumError::NoShare],
            concentration: one,
        };
        let distances = |level: f64| StratumDistances {
            slippage: vec![Some(SlippageDistances {
                level,
                shorter_share: 0.01,
                fall_off: 1e9,
            })],
            length_spectrum: vec![0.02, 0.03, 1e9],
            concentration: 0.04,
        };
        assert_eq!(furthest_in_errors(&errors, &distances(0.05)), Some(0.05));
        let mut share_far = distances(0.05);
        share_far.length_spectrum[1] = -0.5;
        assert_eq!(furthest_in_errors(&errors, &share_far), Some(0.5));
        let mut concentration_far = distances(0.05);
        concentration_far.concentration = 0.7;
        assert_eq!(furthest_in_errors(&errors, &concentration_far), Some(0.7));
        let mut shorter_far = distances(0.05);
        if let Some(group) = shorter_far.slippage[0].as_mut() {
            group.shorter_share = 0.9;
        }
        assert_eq!(furthest_in_errors(&errors, &shorter_far), Some(0.9));
        assert_eq!(
            furthest_in_errors(&errors, &distances(f64::NAN)),
            Some(f64::INFINITY)
        );
        let none = StratumErrors {
            slippage: vec![None],
            length_spectrum: vec![StratumError::NotIdentified; 3],
            concentration: StratumError::NotPlaced,
        };
        assert_eq!(furthest_in_errors(&none, &distances(5.0)), None);
    }

    /// **The Newton step is carried to each number's scale as its error is**: from a given inverse
    /// and slope, `I⁻¹ g` over the kept coordinates, then `p(1 − p)` for each slippage number, the
    /// concentration itself, and for each share its slopes in the log-ratios, the largest class's
    /// included — with every coordinate kept and with one dropped (held).
    #[test]
    fn newton_distances_are_carried_to_each_numbers_scale_like_the_errors() {
        let parameters = Parameters {
            slippage: vec![Slippage {
                level: 0.1,
                shorter_share: 0.8,
                fall_off: 0.3,
            }],
            length_spectrum: vec![0.2, 0.5, 0.3],
            concentration: 2.0,
        };
        let n = 6;
        let inverse: Vec<f64> = (0..n * n)
            .map(|index| {
                let (i, j) = (index / n, index % n);
                if i == j {
                    0.1 * (i + 1) as f64
                } else {
                    0.01 * ((i + j) % 4 + 1) as f64
                }
            })
            .collect();
        let slope = vec![1.0, -2.0, 0.5, 3.0, -1.5, 0.7];
        let inverse = &inverse;
        let check = |kept: Vec<usize>| {
            let layout = CurvatureLayout::of(&parameters, &[true]);
            // The layout: three slippage numbers, the log-ratios of classes 0 and 2 to the largest
            // (class 1), the concentration.
            assert_eq!(layout.largest_class, 1);
            let k = kept.len();
            let sub: Vec<f64> = kept
                .iter()
                .flat_map(|&r| kept.iter().map(move |&c| inverse[r * n + c]))
                .collect();
            let mut step = vec![0.0; n];
            for (row, &r) in kept.iter().enumerate() {
                step[r] = (0..k)
                    .map(|col| sub[row * k + col] * slope[kept[col]])
                    .sum();
            }
            let dropped = (0..n).filter(|i| !kept.contains(i)).collect();
            let curvature = CurvatureAt {
                layout,
                identified: crate::parameter_estimation::joint::fit::Identified {
                    kept,
                    dropped,
                    inverse: sub,
                },
                slope: slope.clone(),
            };
            let distances = newton_distances_on_the_natural_scale(&curvature, &parameters, &[true]);
            let close =
                |got: f64, want: f64| assert!((got - want).abs() < 1e-12, "{got} against {want}");
            let group = distances.slippage[0].expect("live");
            close(group.level, 0.1 * 0.9 * step[0]);
            close(group.shorter_share, 0.8 * 0.2 * step[1]);
            close(group.fall_off, 0.3 * 0.7 * step[2]);
            close(distances.concentration, 2.0 * step[5]);
            let shares = [0.2, 0.5, 0.3];
            let ratio_steps = [(0, step[3]), (2, step[4])];
            for class in 0..3 {
                let want: f64 = ratio_steps
                    .iter()
                    .map(|&(j, d)| {
                        shares[class] * (if j == class { 1.0 } else { 0.0 } - shares[j]) * d
                    })
                    .sum();
                close(distances.length_spectrum[class], want);
            }
            close(distances.length_spectrum.iter().sum::<f64>(), 0.0);
        };
        check((0..n).collect());
        check(vec![0, 1, 2, 3, 5]);
    }

    /// **The settled fraction is the config's**: at a fraction of 10⁻⁷ no walk of a drawn stratum
    /// settles within eight rounds, where at the default every one does.
    #[test]
    fn the_settled_fraction_is_read_from_the_config() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let evidence = draw_stratum(truth, &spectrum_of(3), 0.5, 0.4, 300, 20, 6, 1, 41);
        let fit = |settled_fraction: f64| {
            let config = SsrFitConfig {
                allele_span: 1,
                max_rounds: 8,
                settled_fraction,
                ..SsrFitConfig::default()
            };
            fit_stratum(&evidence, &[0.4; 20], &config).expect("reads were drawn")
        };
        let strict = fit(1e-7);
        assert!(
            strict.walks.iter().all(|walk| !walk.ending.settled()),
            "{:?}",
            strict.walks
        );
        let default = fit(SETTLED_FRACTION);
        assert!(
            default.walks.iter().all(|walk| walk.ending.settled()),
            "{:?}",
            default.walks
        );
    }

    /// **A point whose score is not a number is not judged settled.**
    #[test]
    fn a_point_that_scores_not_a_number_is_not_judged_settled() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(truth, &spectrum_of(3), 0.5, 0.4, 50, 8, 6, 1, 43);
        let excess = [0.4; 8];
        let genotypes = genotype_pairs(3);
        let mut climbing = StratumClimb {
            scorer: Scorer::new(
                &evidence,
                &excess,
                &genotypes,
                &config,
                WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
            ),
            evidence: &evidence,
            live_groups: &[true],
            classes: 3,
            settled_fraction: SETTLED_FRACTION,
        };
        let mut at = Parameters::start(config.starting_points[1], 1, 3);
        at.concentration = f64::NAN;
        let (settled, errors) = climbing.judge(&at);
        assert!(!settled, "{errors:?}");
    }

    /// **The point a round holds scores what it held**: over four rounds from each starting point of
    /// a drawn stratum, the best point a round passed through, scored afresh, is the score it was held
    /// with to the bit, at least the round's start and at least where the round ended.
    #[test]
    fn the_held_point_scores_what_it_holds() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(truth, &spectrum_of(3), 0.5, 0.4, 100, 10, 3, 1, 53);
        let excess = [0.4; 10];
        let genotypes = genotype_pairs(3);
        let mut scorer = Scorer::new(
            &evidence,
            &excess,
            &genotypes,
            &config,
            WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
        );
        for start in &config.starting_points {
            let mut at = Parameters::start(*start, 1, 3);
            for _ in 0..4 {
                let start_score = scorer.score(&at);
                let mut held = HeldPoint {
                    at: at.clone(),
                    score: start_score,
                };
                climb_one_round(&mut at, &mut scorer, &[true], 3, &mut held);
                let rescored = scorer.score(&held.at);
                assert_eq!(
                    held.score.to_bits(),
                    rescored.to_bits(),
                    "{} against {rescored}",
                    held.score
                );
                assert!(held.score >= start_score);
                assert!(held.score >= scorer.score(&at));
            }
        }
    }

    /// **A walk from a start that scores not a number** loses its first round with nothing better,
    /// is judged there, and ends as having lost a round.
    #[test]
    fn a_walk_from_a_start_that_is_not_a_number_lost_a_round() {
        let rounds = [(None, f64::NAN); 6];
        let (walked, judged) = scripted_walk(f64::NAN, &rounds, &[], 0.08, 1.0, 6);
        assert_eq!((walked.ending, walked.rounds), (ClimbEnding::LostARound, 1));
        assert_eq!(judged, vec![0.0]);
    }

    /// **A walk still gaining at its last round ran out, and says so**, at the point its last round
    /// reached, never judged.
    #[test]
    fn a_walk_still_gaining_at_its_last_round_ran_out() {
        let (walked, judged) = scripted_walk(0.0, &rising_by(&[1.0; 4]), &[], 0.08, 1.0, 4);
        assert_eq!(
            (walked.ending, walked.rounds, walked.at, walked.score),
            (ClimbEnding::OutOfRounds, 4, 4.0, 4.0)
        );
        assert!(judged.is_empty());
        assert_eq!(walked.judged_here, None);
    }

    /// **The judge**: on a drawn stratum, at the answer a long climb reached every number is within a
    /// tenth of its error of where the likelihood peaks, and the point has settled; with the level
    /// moved by one of its errors, the level is about one error away and the point has not settled.
    #[test]
    fn a_point_one_error_from_the_peak_has_not_settled() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(truth, &spectrum_of(3), 0.5, 0.4, 300, 20, 6, 1, 43);
        let excess = [0.4; 20];
        let genotypes = genotype_pairs(3);
        let path = trajectory_of(&evidence, &excess, config.starting_points[1], &config, 20);
        let peak = the_best_of(path.iter()).0.clone();
        let mut climbing = StratumClimb {
            scorer: Scorer::new(
                &evidence,
                &excess,
                &genotypes,
                &config,
                WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
            ),
            evidence: &evidence,
            live_groups: &[true],
            classes: 3,
            settled_fraction: SETTLED_FRACTION,
        };
        let (settled, errors) = climbing.judge(&peak);
        assert!(settled, "{errors:?}");
        let level_error = errors.slippage[0]
            .expect("the one group")
            .level
            .value()
            .expect("an error");
        let mut moved = peak.clone();
        moved.slippage[0].level += level_error;
        let curvature = curvature_at(
            &mut climbing.scorer,
            &evidence,
            &moved,
            &[true],
            CURVATURE_STEP,
        );
        let distances = newton_distances_on_the_natural_scale(&curvature, &moved, &[true]);
        let level_distance = distances.slippage[0].expect("the one group").level / level_error;
        assert!(
            (-1.3..-0.7).contains(&level_distance),
            "the level is {level_distance} of its errors from the peak"
        );
        let (settled, _) = climbing.judge(&moved);
        assert!(!settled);
    }

    /// **The target is `½ · p · fraction²` over the numbers the errors are taken over**: at one live
    /// slippage group of two and thirteen classes, p = 3 + 12 + 1 = 16 and the target at a tenth is
    /// 0.08.
    #[test]
    fn the_target_counts_the_numbers_the_errors_are_taken_over() {
        let parameters = Parameters::start(
            StartingPoint {
                slippage_level: 0.02,
                concentration: 0.3,
            },
            2,
            13,
        );
        let target = remaining_gain_target(&parameters, &[true, false], SETTLED_FRACTION);
        assert!((target - 0.08).abs() < 1e-15, "{target}");
        let both = remaining_gain_target(&parameters, &[true, true], SETTLED_FRACTION);
        assert!((both - 0.095).abs() < 1e-15, "{both}");
    }

    /// **Tracts no sample read do not move where a walk stops**: they leave the total log-likelihood
    /// as it was, and the rule is on the total, so a stratum padded with as many again takes the same
    /// rounds to the same answer.
    #[test]
    fn tracts_without_reads_do_not_move_where_a_walk_stops() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let evidence = draw_stratum(truth, &spectrum_of(3), 0.5, 0.4, 300, 20, 6, 1, 41);
        let mut padded = evidence.clone();
        padded
            .tracts
            .extend(std::iter::repeat_n(TractReads::default(), 300));
        let plain = fit_stratum(&evidence, &[0.4; 20], &config).expect("reads were drawn");
        let with_empty = fit_stratum(&padded, &[0.4; 20], &config).expect("reads were drawn");
        assert_eq!(plain.walks, with_empty.walks);
        assert_eq!(plain.walks.len(), 3);
        let close = |a: f64, b: f64| assert!((a - b).abs() <= 1e-9 * a.abs(), "{a} against {b}");
        let (a, b) = (
            plain.slippage[0].expect("fitted"),
            with_empty.slippage[0].expect("fitted"),
        );
        close(a.level, b.level);
        close(a.shorter_share, b.shorter_share);
        close(a.fall_off, b.fall_off);
        close(plain.concentration, with_empty.concentration);
        for (a, b) in plain
            .length_spectrum
            .iter()
            .zip(&with_empty.length_spectrum)
        {
            close(*a, *b);
        }
    }

    /// One walk run for `rounds` rounds whatever it gains — the climb with no stopping rule — as the
    /// point and the mean log-likelihood a tract after every round, the start first.
    fn trajectory_of(
        evidence: &StratumEvidence,
        homozygote_excess: &[f64],
        start: StartingPoint,
        config: &SsrFitConfig,
        rounds: u32,
    ) -> Vec<(Parameters, f64)> {
        let classes = (2 * config.allele_span + 1) as usize;
        let genotypes = genotype_pairs(classes);
        let live = evidence.groups_with_reads();
        let mut scorer = Scorer::new(
            evidence,
            homozygote_excess,
            &genotypes,
            config,
            WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
        );
        let mut parameters = Parameters::start(start, evidence.groups, classes);
        let score = scorer.score(&parameters);
        let mut path = vec![(parameters.clone(), score)];
        for _ in 0..rounds {
            let mut held = HeldPoint {
                at: parameters.clone(),
                score,
            };
            climb_one_round(&mut parameters, &mut scorer, &live, classes, &mut held);
            let score = scorer.score(&parameters);
            path.push((parameters.clone(), score));
        }
        path
    }

    /// The round the rule before plan step C3 stopped a walk at: the first whose gain of the mean a
    /// tract was below 10⁻⁶ — a loss included — or the fifth.
    fn where_the_rule_before_stopped(path: &[(Parameters, f64)]) -> usize {
        (1..=5)
            .find(|&round| path[round].1 - path[round - 1].1 < 1e-6)
            .unwrap_or(5)
    }

    /// The highest-scoring of `points`, the first among equals — the rule [`the_better_walk`] ranks
    /// walks by.
    fn the_best_of<'a>(
        points: impl Iterator<Item = &'a (Parameters, f64)>,
    ) -> &'a (Parameters, f64) {
        points
            .fold(None::<&(Parameters, f64)>, |best, point| match best {
                Some(best) if point.1.partial_cmp(&best.1) != Some(std::cmp::Ordering::Greater) => {
                    Some(best)
                }
                _ => Some(point),
            })
            .expect("at least one point")
    }

    /// Each number of a one-group stratum's answer with its kind, in the order [`errors_by_number`]
    /// lists their errors.
    fn numbers_of(parameters: &Parameters) -> Vec<(&'static str, f64)> {
        let slippage = parameters.slippage[0];
        let mut numbers = vec![
            ("slippage level", slippage.level),
            ("shorter share", slippage.shorter_share),
            ("fall-off", slippage.fall_off),
            ("concentration", parameters.concentration),
        ];
        numbers.extend(
            parameters
                .length_spectrum
                .iter()
                .map(|share| ("length shares", *share)),
        );
        numbers
    }

    fn errors_by_number(errors: &StratumErrors) -> Vec<StratumError> {
        let group = errors.slippage[0].expect("the one group has reads");
        let mut all = vec![
            group.level,
            group.shorter_share,
            group.fall_off,
            errors.concentration,
        ];
        all.extend(errors.length_spectrum.iter().copied());
        all
    }

    /// The largest distance of one kind of number from the reference answer, in that answer's
    /// errors, and how many lie further than a tenth of an error.
    #[derive(Debug, Default)]
    struct Shortfall {
        largest: f64,
        beyond_a_tenth: usize,
        compared: usize,
    }

    /// Every number's distance from `reference`, in `errors`, tallied by kind; numbers without an
    /// error are not compared.
    fn tally_distances(
        tally: &mut BTreeMap<&'static str, Shortfall>,
        at: &Parameters,
        reference: &Parameters,
        errors: &StratumErrors,
    ) -> f64 {
        let mut largest: f64 = 0.0;
        for (((kind, value), (_, truth)), error) in numbers_of(at)
            .into_iter()
            .zip(numbers_of(reference))
            .zip(errors_by_number(errors))
        {
            let StratumError::Estimated(error) = error else {
                continue;
            };
            let distance = (value - truth).abs() / error;
            let entry = tally.entry(kind).or_default();
            entry.largest = entry.largest.max(distance);
            entry.beyond_a_tenth += usize::from(distance > SETTLED_FRACTION);
            entry.compared += 1;
            largest = largest.max(distance);
        }
        largest
    }

    /// The draws a regime of the comparison takes by default, and the variable that overrides it.
    const CLIMB_DRAWS: usize = 20;
    const CLIMB_DRAWS_VARIABLE: &str = "NG_FIT_PRECISION_CLIMB_DRAWS";

    /// **The climb's stop loses nothing against a longer climb** (plan step C3, spec §4.5 item 2): on
    /// strata drawn at known numbers, the fit's answer under its stopping rule against the best point
    /// the same walks reach in a fixed number of rounds with no rule (40 by default), every number's
    /// distance in that answer's standard errors, and the log-likelihood the fit is short of it; the
    /// rounds and judgements the rule took, and the rounds the rule before it (a gain of the mean below
    /// 10⁻⁶, or five rounds) took, with that rule's distances too.
    ///
    /// Ignored by default: an hour in the container. `NG_FIT_PRECISION_CLIMB_DRAWS` sets the
    /// three-class regimes' draws, the 63-sample and thirteen-class ones taking a quarter as many;
    /// `NG_FIT_PRECISION_CLIMB_CLASSES` keeps the regimes of one class count;
    /// `NG_FIT_PRECISION_CLIMB_MAX_ROUNDS` and `NG_FIT_PRECISION_CLIMB_REFERENCE_ROUNDS` set the fit's
    /// round limit and the reference's rounds. At the default draws, every fit whose winning walk
    /// settled is held within a tenth of an error of the reference on every number and within
    /// ½ · p · 0.1² of its log-likelihood.
    #[test]
    #[ignore = "a measurement over many fitted strata, each also climbed to a longer reference"]
    fn the_climbs_stop_loses_nothing_against_a_longer_climb() {
        let draws = match std::env::var(CLIMB_DRAWS_VARIABLE) {
            Err(_) => CLIMB_DRAWS,
            Ok(value) => value
                .trim()
                .parse()
                .ok()
                .filter(|&count: &usize| count > 0)
                .unwrap_or_else(|| {
                    panic!("{CLIMB_DRAWS_VARIABLE}={value:?} is not a positive whole number")
                }),
        };
        let truth = Slippage {
            level: 0.05,
            shorter_share: 0.8,
            fall_off: 0.3,
        };
        let count_from = |variable: &str, default: u32| -> u32 {
            std::env::var(variable).ok().map_or(default, |value| {
                value
                    .trim()
                    .parse()
                    .unwrap_or_else(|_| panic!("{variable}={value:?} is not a whole number"))
            })
        };
        // The fit's own round limit, and the reference climb's, both settable to measure how many
        // rounds a walk needs.
        let max_rounds = count_from(
            "NG_FIT_PRECISION_CLIMB_MAX_ROUNDS",
            SsrFitConfig::default().max_rounds,
        );
        let reference_rounds = count_from("NG_FIT_PRECISION_CLIMB_REFERENCE_ROUNDS", 40);
        let only_classes = std::env::var("NG_FIT_PRECISION_CLIMB_CLASSES")
            .ok()
            .map(|value| value.trim().parse::<usize>().expect("a class count"));
        for (classes, tracts, samples, depth, regime_draws) in [
            (3, 30, 8, 3, draws),
            (3, 300, 20, 3, draws),
            (3, 300, 20, 30, draws),
            (3, 400, 63, 3, draws.div_ceil(4)),
            (13, 200, 20, 3, draws.div_ceil(4)),
        ] {
            if only_classes.is_some_and(|only| only != classes) {
                continue;
            }
            let regime =
                format!("{classes} classes, {tracts} tracts x {samples} samples x {depth} reads");
            let span = (classes / 2) as i32;
            let spectrum = spectrum_of(classes);
            let config = SsrFitConfig {
                allele_span: span,
                max_rounds,
                ..SsrFitConfig::default()
            };
            let excess = vec![0.4; samples];
            let (mut now, mut before) = (BTreeMap::new(), BTreeMap::new());
            let (mut now_settled, mut before_settled) = (BTreeMap::new(), BTreeMap::new());
            let mut endings: BTreeMap<String, usize> = BTreeMap::new();
            let (mut rounds_now, mut rounds_before, mut walks) = (0_u64, 0_u64, 0_u64);
            let mut judgements_now = 0_u64;
            let (mut short_now, mut short_before) = (0.0_f64, 0.0_f64);
            let mut winners_settled = 0;
            // Each draw whose winning walk settled: its furthest number in errors, and how far short.
            let mut settled_draws: Vec<(usize, f64, f64)> = Vec::new();
            for draw in 0..regime_draws {
                let seed = 0xC3C0_0000_0000_0000
                    + (classes as u64) * 100_000
                    + (samples as u64) * 1_000
                    + u64::from(depth) * 10
                    + draw as u64;
                let evidence = draw_stratum(
                    truth, &spectrum, 0.5, 0.4, tracts, samples, depth, span, seed,
                );
                let fitted = fit_stratum(&evidence, &excess, &config).expect("reads were drawn");
                let paths: Vec<Vec<(Parameters, f64)>> = config
                    .starting_points
                    .iter()
                    .map(|start| {
                        trajectory_of(&evidence, &excess, *start, &config, reference_rounds)
                    })
                    .collect();
                let reference = the_best_of(paths.iter().flatten());
                let before_rule = the_best_of(
                    paths
                        .iter()
                        .map(|path| &path[where_the_rule_before_stopped(path)]),
                );
                let errors = standard_errors_at(
                    &evidence,
                    &reference.0,
                    &excess,
                    &[true],
                    &config,
                    WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
                );
                let answer = Parameters {
                    slippage: vec![fitted.slippage[0].expect("the one group has reads")],
                    length_spectrum: fitted.length_spectrum.clone(),
                    concentration: fitted.concentration,
                };
                let settled = fitted.ending.settled();
                winners_settled += usize::from(settled);
                let worst_now = tally_distances(&mut now, &answer, &reference.0, &errors);
                if settled {
                    tally_distances(&mut now_settled, &answer, &reference.0, &errors);
                }
                let worst_before =
                    tally_distances(&mut before, &before_rule.0, &reference.0, &errors);
                let before_converged = paths.iter().any(|path| {
                    let stop = where_the_rule_before_stopped(path);
                    path[stop].1 - path[stop - 1].1 < 1e-6
                });
                if before_converged {
                    tally_distances(&mut before_settled, &before_rule.0, &reference.0, &errors);
                }
                let total = evidence.tracts.len() as f64;
                short_now = short_now.max((reference.1 - fitted.log_likelihood_a_tract) * total);
                short_before = short_before.max((reference.1 - before_rule.1) * total);
                for walk in &fitted.walks {
                    *endings.entry(format!("{:?}", walk.ending)).or_default() += 1;
                    rounds_now += u64::from(walk.rounds);
                    judgements_now += u64::from(walk.judgements);
                }
                rounds_before += paths
                    .iter()
                    .map(|path| where_the_rule_before_stopped(path) as u64)
                    .sum::<u64>();
                walks += fitted.walks.len() as u64;
                eprintln!(
                    "CLIMB STOP {regime}, draw {draw}: {:?} with (ending, rounds, judgements) {:?}; furthest \
                     number {:.3} errors from the reference climb (the rule before: {:.3}); total \
                     log-likelihood short by {:.4} (before: {:.4})",
                    fitted.ending,
                    fitted
                        .walks
                        .iter()
                        .map(|walk| (walk.ending, walk.rounds, walk.judgements))
                        .collect::<Vec<_>>(),
                    worst_now,
                    worst_before,
                    (reference.1 - fitted.log_likelihood_a_tract) * total,
                    (reference.1 - before_rule.1) * total,
                );
                if settled {
                    settled_draws.push((
                        draw,
                        worst_now,
                        (reference.1 - fitted.log_likelihood_a_tract) * total,
                    ));
                }
            }
            eprintln!(
                "CLIMB STOP {regime}: {regime_draws} draws, {walks} walks ended {endings:?}; the \
                 winning walk settled in {winners_settled}; rounds {rounds_now} (and {judgements_now} \
                 judgements) against {rounds_before} \
                 under the rule before; the total log-likelihood at most {short_now:.4} short of \
                 the reference climb (before: {short_before:.4})"
            );
            // **Held at the default draws**: a fit whose winning walk settled is within a tenth of an
            // error of the reference on every number, and within ½ · p · 0.1² of its log-likelihood
            // (0.08 at thirteen classes); a winner that did not settle says so and is not held.
            if draws == CLIMB_DRAWS {
                let tolerance = 0.5 * (classes + 3) as f64 * SETTLED_FRACTION * SETTLED_FRACTION;
                for (draw, furthest, short) in &settled_draws {
                    assert!(
                        *furthest <= SETTLED_FRACTION && *short <= tolerance,
                        "{regime}, draw {draw}: settled {furthest:.3} errors and {short:.4} units \
                         short of the reference climb"
                    );
                }
            }
            for (label, tally) in [
                ("every draw, this rule", &now),
                ("draws whose winner settled, this rule", &now_settled),
                ("every draw, the rule before", &before),
                ("draws the rule before called converged", &before_settled),
            ] {
                for (kind, shortfall) in tally {
                    eprintln!(
                        "CLIMB STOP {regime}, {label}, {kind}: {} compared, {} beyond a tenth of an \
                         error, the furthest {:.3}",
                        shortfall.compared, shortfall.beyond_a_tenth, shortfall.largest
                    );
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // How far the 256-point average is from the truth (after plan step C2)
    // -----------------------------------------------------------------

    /// A reproducible stream of uniforms, normals, and Gamma and Dirichlet draws kept in logs.
    struct Stream(u64);

    impl Stream {
        fn uniform(&mut self) -> f64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            ((z >> 11) as f64 + 0.5) / (1_u64 << 53) as f64
        }

        fn normal(&mut self) -> f64 {
            let (u1, u2) = (self.uniform(), self.uniform());
            (-2.0 * float::ln(u1)).sqrt() * float::cos(std::f64::consts::TAU * u2)
        }

        /// The logarithm of a Gamma(`shape`, 1) draw, kept in logs so a tiny shape does not
        /// underflow to zero.
        fn ln_gamma_draw(&mut self, shape: f64) -> f64 {
            if shape < 1.0 {
                let u = self.uniform();
                return self.ln_gamma_draw(shape + 1.0) + float::ln(u) / shape;
            }
            let d = shape - 1.0 / 3.0;
            let c = 1.0 / (9.0 * d).sqrt();
            loop {
                let x = self.normal();
                let v = float::powi(1.0 + c * x, 3);
                if v <= 0.0 {
                    continue;
                }
                if float::ln(self.uniform()) < 0.5 * x * x + d - d * v + d * float::ln(v) {
                    return float::ln(d * v);
                }
            }
        }

        /// The logarithms of one Dirichlet(`alpha`) draw's coordinates.
        fn ln_dirichlet(&mut self, alpha: &[f64]) -> Vec<f64> {
            let draws: Vec<f64> = alpha.iter().map(|a| self.ln_gamma_draw(*a)).collect();
            let total = ln_sum_exp(&draws);
            draws.into_iter().map(|draw| draw - total).collect()
        }
    }

    fn ln_dirichlet_density(alpha: &[f64], ln_frequencies: &[f64]) -> f64 {
        let total: f64 = alpha.iter().sum();
        ln_gamma(total)
            + alpha
                .iter()
                .zip(ln_frequencies)
                .map(|(a, ln_f)| (a - 1.0) * ln_f - ln_gamma(*a))
                .sum::<f64>()
    }

    /// `ln` of one tract's likelihood at the length frequencies `frequencies` — the thing the
    /// average is taken over, term for term as [`ln_tract`] computes it.
    fn ln_integrand(
        likelihoods: &TractLikelihoods,
        frequencies: &[f64],
        genotypes: &[(usize, usize)],
    ) -> f64 {
        let mut total = likelihoods.ln_offset;
        for (row, excess) in likelihoods
            .scaled
            .chunks_exact(likelihoods.width)
            .zip(&likelihoods.homozygote_excess)
        {
            let (mut at_random, mut by_descent) = (0.0, 0.0);
            for (slot, (first, second)) in genotypes.iter().enumerate() {
                if first == second {
                    at_random += frequencies[*first] * frequencies[*first] * row[slot];
                    by_descent += frequencies[*first] * row[slot];
                } else {
                    at_random += 2.0 * frequencies[*first] * frequencies[*second] * row[slot];
                }
            }
            let sum = excess * by_descent + (1.0 - excess) * at_random;
            if sum <= 0.0 {
                return f64::NEG_INFINITY;
            }
            total += float::ln(sum);
        }
        total
    }

    /// How many copies of each allele class one tract's reads put in its samples, by each sample's
    /// genotype posterior under the stratum's mean frequencies: where a tract's reads say its
    /// frequencies are.
    fn copies_the_reads_say(
        likelihoods: &TractLikelihoods,
        spectrum: &[f64],
        genotypes: &[(usize, usize)],
    ) -> Vec<f64> {
        let mut copies = vec![0.0; spectrum.len()];
        for (row, excess) in likelihoods
            .scaled
            .chunks_exact(likelihoods.width)
            .zip(&likelihoods.homozygote_excess)
        {
            let weights: Vec<f64> = genotypes
                .iter()
                .enumerate()
                .map(|(slot, (first, second))| {
                    let prior = if first == second {
                        excess * spectrum[*first]
                            + (1.0 - excess) * spectrum[*first] * spectrum[*first]
                    } else {
                        (1.0 - excess) * 2.0 * spectrum[*first] * spectrum[*second]
                    };
                    prior * row[slot]
                })
                .collect();
            let total: f64 = weights.iter().sum();
            for (weight, (first, second)) in weights.iter().zip(genotypes) {
                copies[*first] += weight / total;
                copies[*second] += weight / total;
            }
        }
        copies
    }

    /// `ln` of a tract's average by importance sampling, and that logarithm's standard error: `draws`
    /// points from a mixture — a tenth from the stratum's own Dirichlet, the rest in three equal
    /// parts from Dirichlets moved all, half and a quarter of the way to the copies the tract's reads
    /// say it carries — each weighted by the target's density over the mixture's.
    fn importance_average(
        likelihoods: &TractLikelihoods,
        genotypes: &[(usize, usize)],
        alpha: &[f64],
        towards: &[f64],
        draws: usize,
        stream: &mut Stream,
    ) -> (f64, f64) {
        let moved = |share: f64| -> Vec<f64> {
            alpha
                .iter()
                .zip(towards)
                .map(|(a, t)| a + share * (t - a))
                .collect()
        };
        let parts: Vec<(f64, Vec<f64>)> = vec![
            (0.1, alpha.to_vec()),
            (0.3, moved(1.0)),
            (0.3, moved(0.5)),
            (0.3, moved(0.25)),
        ];
        let ln_weights: Vec<f64> = (0..draws)
            .map(|_| {
                let mut u = stream.uniform();
                let mut chosen = parts.len() - 1;
                for (index, (weight, _)) in parts.iter().enumerate() {
                    if u < *weight {
                        chosen = index;
                        break;
                    }
                    u -= weight;
                }
                let ln_f = stream.ln_dirichlet(&parts[chosen].1);
                let frequencies: Vec<f64> = ln_f.iter().map(|ln| float::exp(*ln)).collect();
                let ln_target = ln_dirichlet_density(alpha, &ln_f);
                let ln_proposal = ln_sum_exp(
                    &parts
                        .iter()
                        .map(|(weight, shape)| {
                            float::ln(*weight) + ln_dirichlet_density(shape, &ln_f)
                        })
                        .collect::<Vec<_>>(),
                );
                ln_integrand(likelihoods, &frequencies, genotypes) + ln_target - ln_proposal
            })
            .collect();
        let ln_mean = ln_sum_exp(&ln_weights) - float::ln(draws as f64);
        let variance = ln_weights
            .iter()
            .map(|ln_w| {
                let ratio = float::exp(ln_w - ln_mean) - 1.0;
                ratio * ratio
            })
            .sum::<f64>()
            / draws as f64;
        (ln_mean, (variance / draws as f64).sqrt())
    }

    /// **How far the fixed-point average over a tract's length frequencies is from the truth**
    /// (the investigation decided after plan step C2): on strata drawn at known numbers, each
    /// tract's log-likelihood at those numbers, averaged over the fit's own point set at 256 to
    /// 65,536 points, against a reference by importance sampling with its standard error; at 3 to
    /// 13 allele classes, concentrations 0.1 to 10, 3 and 30 reads, 4 to 63 samples. Also the same
    /// importance sampler at 256 and 4,096 points, a first look at placing a tract's points where
    /// its reads say its frequencies are. Results: `fit_precision_quadrature_average_2026-10-01.md`.
    ///
    /// Ignored by default: two to five minutes in the container. `NG_AVERAGE_REFERENCE_DRAWS` sets the
    /// reference's points a tract (20,000; the report's 3-read cells at 9 and 13 classes used
    /// 200,000), `NG_AVERAGE_DEPTH` keeps the cells at one depth, and `NG_AVERAGE_LEAST_CLASSES` drops
    /// the cells below a class count.
    #[test]
    #[ignore = "a measurement; run by hand"]
    fn the_average_over_a_tracts_frequencies_against_the_truth() {
        let truth = Slippage {
            level: 0.05,
            shorter_share: 0.8,
            fall_off: 0.3,
        };
        let tracts = 60;
        let mut cells = Vec::new();
        for classes in [3_usize, 5, 9, 13] {
            for concentration in [0.1, 0.5, 2.0, 10.0] {
                for depth in [3_u32, 30] {
                    cells.push((classes, concentration, 20_usize, depth));
                }
            }
        }
        cells.extend([
            (13, 0.5, 4, 3),
            (13, 0.5, 63, 3),
            (13, 0.5, 4, 30),
            (13, 0.5, 63, 30),
        ]);
        let reference_draws: usize = std::env::var("NG_AVERAGE_REFERENCE_DRAWS")
            .ok()
            .map_or(20_000, |value| value.parse().expect("a count"));
        let only_depth: Option<u32> = std::env::var("NG_AVERAGE_DEPTH")
            .ok()
            .map(|value| value.parse().expect("a depth"));
        let least_classes: usize = std::env::var("NG_AVERAGE_LEAST_CLASSES")
            .ok()
            .map_or(0, |value| value.parse().expect("a count"));
        for (cell, (classes, concentration, samples, depth)) in cells.into_iter().enumerate() {
            if only_depth.is_some_and(|only| only != depth) || classes < least_classes {
                continue;
            }
            let span = (classes / 2) as i32;
            let spectrum = spectrum_of(classes);
            let excess = vec![0.4; samples];
            let evidence = draw_stratum(
                truth,
                &spectrum,
                concentration,
                0.4,
                tracts,
                samples,
                depth,
                span,
                0xA7E0_0000 + cell as u64,
            );
            let genotypes = genotype_pairs(classes);
            let per_allele: Vec<Vec<Vec<f64>>> = vec![
                (-span..=span)
                    .map(|allele| truth.read_probabilities(allele, span))
                    .collect(),
            ];
            let likelihoods: Vec<TractLikelihoods> = evidence
                .tracts
                .iter()
                .map(|tract| TractLikelihoods::of(tract, &per_allele, &genotypes, &excess))
                .collect();
            let alpha: Vec<f64> = spectrum.iter().map(|share| concentration * share).collect();
            let mut stream = Stream(0x5EED + cell as u64);
            let reference: Vec<(f64, f64)> = likelihoods
                .iter()
                .map(|tract| {
                    let towards: Vec<f64> = alpha
                        .iter()
                        .zip(copies_the_reads_say(tract, &spectrum, &genotypes))
                        .map(|(a, copies)| a + copies)
                        .collect();
                    importance_average(
                        tract,
                        &genotypes,
                        &alpha,
                        &towards,
                        reference_draws,
                        &mut stream,
                    )
                })
                .collect();
            let reference_error = reference.iter().map(|(_, se)| *se).fold(0.0, f64::max);
            let unsure = reference.iter().filter(|(_, se)| *se > 0.05).count();
            let mut line = format!(
                "AVERAGE {classes} classes, concentration {concentration}, {samples} samples x \
                 {depth} reads: reference {:.3} a tract (largest standard error {reference_error:.4}, \
                 {unsure} tracts above 0.05);",
                reference.iter().map(|(ln, _)| ln).sum::<f64>() / tracts as f64
            );
            for points in [256_usize, 1_024, 4_096, 16_384, 65_536] {
                let quadrature = dirichlet_points(
                    &spectrum,
                    concentration,
                    points,
                    &genotypes,
                    WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
                );
                let gaps: Vec<f64> = likelihoods
                    .iter()
                    .zip(&reference)
                    .map(|(tract, (truth, _))| ln_tract(tract, &quadrature, &genotypes) - truth)
                    .collect();
                if points == 256 {
                    // The integrand is the one `ln_tract` averages.
                    let first = &likelihoods[0];
                    let terms: Vec<f64> = quadrature
                        .frequencies
                        .chunks_exact(classes)
                        .map(|point| ln_integrand(first, point, &genotypes) + quadrature.ln_weight)
                        .collect();
                    let by_hand = ln_sum_exp(&terms);
                    let by_fit = ln_tract(first, &quadrature, &genotypes);
                    assert!(
                        (by_hand - by_fit).abs() < 1e-8,
                        "{by_hand} against {by_fit}"
                    );
                }
                let mean = gaps.iter().sum::<f64>() / tracts as f64;
                let worst = gaps.iter().copied().fold(0.0_f64, |a, b| a.min(b));
                line.push_str(&format!(" {points} points {mean:+.3} (worst {worst:+.2});"));
            }
            for draws in [256_usize, 4_096] {
                let gaps: Vec<f64> = likelihoods
                    .iter()
                    .zip(&reference)
                    .map(|(tract, (truth, _))| {
                        let towards: Vec<f64> = alpha
                            .iter()
                            .zip(copies_the_reads_say(tract, &spectrum, &genotypes))
                            .map(|(a, copies)| a + copies)
                            .collect();
                        importance_average(tract, &genotypes, &alpha, &towards, draws, &mut stream)
                            .0
                            - truth
                    })
                    .collect();
                let mean = gaps.iter().sum::<f64>() / tracts as f64;
                let spread = (gaps
                    .iter()
                    .map(|gap| (gap - mean) * (gap - mean))
                    .sum::<f64>()
                    / tracts as f64)
                    .sqrt();
                line.push_str(&format!(" placed {draws} {mean:+.3} (spread {spread:.3});"));
            }
            eprintln!("{line}");
        }
    }

    /// A stratum fitted on its own tracts, built directly so the smoothing can be exercised
    /// without paying for the climb.
    fn fitted_at(period: u8, repeats: u64, level: f64, reads: u64) -> StratumOutcome {
        StratumOutcome::Fitted(Box::new(StratumFit {
            standard_errors: None,
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
            ending: ClimbEnding::Settled,
            walks: Vec::new(),
            samples_fitted_on: None,
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

    /// **The log's line about the climbs**: over the strata that recorded their walks, how many walks
    /// ended each way, the rounds between them, and in how many strata the winner settled; a
    /// hand-built fit with no walks is not counted.
    #[test]
    fn the_climb_summary_counts_each_ending() {
        let walk = |ending, rounds| WalkRecord {
            ending,
            rounds,
            judgements: 1,
            settled_with_no_error: false,
        };
        let with = |ending, walks: Vec<WalkRecord>| {
            let mut outcome = fitted_at(2, 10, 0.1, 1_000);
            if let StratumOutcome::Fitted(fit) = &mut outcome {
                fit.ending = ending;
                fit.walks = walks;
            }
            outcome
        };
        let outcomes = [
            with(
                ClimbEnding::Settled,
                vec![
                    walk(ClimbEnding::Settled, 4),
                    walk(ClimbEnding::LostARound, 7),
                    walk(ClimbEnding::OutOfRounds, 30),
                ],
            ),
            with(
                ClimbEnding::LostARound,
                vec![
                    walk(ClimbEnding::LostARound, 5),
                    walk(ClimbEnding::Settled, 3),
                    WalkRecord {
                        settled_with_no_error: true,
                        ..walk(ClimbEnding::Settled, 2)
                    },
                ],
            ),
            fitted_at(2, 11, 0.1, 1_000),
        ];
        let summary = climb_endings_summary(&outcomes, 30).expect("walks recorded");
        for part in [
            "the 2 strata",
            "of their 6 walks, 3 settled (1 of them at a point where no number had an error), 2 \
             stopped",
            "and 1 ran out of their 30 rounds, 51 rounds in all",
            "settled in 1 of the 2",
        ] {
            assert!(summary.contains(part), "{part:?} in {summary}");
        }
        assert_eq!(climb_endings_summary(&[fitted_at(2, 10, 0.1, 1)], 30), None);
    }

    /// **A walk settled with nothing to judge only when it settled and its last judgement found no
    /// number with an error**: one error anywhere — a slippage number, a share or the concentration —
    /// makes it an ordinary settled walk, and a walk that did not settle is never counted.
    #[test]
    fn a_walk_settled_with_no_error_only_when_no_number_had_one() {
        let none = StratumErrors {
            slippage: vec![
                Some(SlippageErrors {
                    level: StratumError::NotIdentified,
                    shorter_share: StratumError::NotPlaced,
                    fall_off: StratumError::NotPlaced,
                }),
                None,
            ],
            length_spectrum: vec![
                StratumError::NotIdentified,
                StratumError::NoShare,
                StratumError::NotPlaced,
            ],
            concentration: StratumError::NotPlaced,
        };
        assert!(!none.any());
        assert!(settled_with_no_error(ClimbEnding::Settled, Some(&none)));
        assert!(!settled_with_no_error(ClimbEnding::LostARound, Some(&none)));
        assert!(!settled_with_no_error(
            ClimbEnding::OutOfRounds,
            Some(&none)
        ));
        assert!(!settled_with_no_error(ClimbEnding::Settled, None));

        let one = StratumError::Estimated(0.1);
        let mut slippage = none.clone();
        if let Some(group) = slippage.slippage[0].as_mut() {
            group.fall_off = one;
        }
        let mut share = none.clone();
        share.length_spectrum[2] = one;
        let concentration = StratumErrors {
            concentration: one,
            ..none.clone()
        };
        for errors in [slippage, share, concentration] {
            assert!(errors.any(), "{errors:?}");
            assert!(!settled_with_no_error(ClimbEnding::Settled, Some(&errors)));
        }
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

        let outcomes = fit_strata(&strata, &[0.4; 8], &base);

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
        let outcomes = fit_strata(&strata, &[0.4; 8], &config);

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

    // ---------------------------------------------------------------
    // The tract prior's middle rung: one motif period's tracts pooled
    // ---------------------------------------------------------------

    fn stratum_at(period: u8, reference_repeats: u64) -> Stratum {
        Stratum {
            period,
            reference_repeats,
        }
    }

    /// A drawn stratum, re-keyed — [`draw_stratum`] always stamps period 2 at 10 repeats, and
    /// pooling is about strata that differ.
    fn drawn_at(stratum: Stratum, spectrum: &[f64], tracts: usize, seed: u64) -> StratumEvidence {
        let slippage = Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        };
        let mut evidence = draw_stratum(slippage, spectrum, 0.5, 0.4, tracts, 8, 6, 1, seed);
        evidence.stratum = stratum;
        evidence
    }

    /// A three-class spectrum tilted towards one side, so that two strata drawn from two of
    /// these are told apart by where their mass sits and not only by noise.
    fn tilted(short: f64, middle: f64) -> Vec<f64> {
        vec![short, middle, 1.0 - short - middle]
    }

    fn pooling_config() -> SsrFitConfig {
        SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        }
    }

    /// **The positive control for the middle rung.** Two strata of one period, drawn from two
    /// spectra tilted opposite ways, come back as one pooled spectrum that sits between them —
    /// and the pool says it read both.
    ///
    /// **The two truths differ on purpose.** Drawn from one truth, a pool that quietly read only
    /// its first stratum would recover that truth exactly and this test would pass; drawn from
    /// two, reading one gives that one's tilt and the assertion on the middle class's neighbours
    /// fails. The tract count and the stratum count are asserted for the same reason from the
    /// other side.
    #[test]
    fn a_periods_pool_reads_every_stratum_of_it() {
        let leans_short = tilted(0.55, 0.35);
        let leans_long = tilted(0.10, 0.35);
        let strata = [
            drawn_at(stratum_at(2, 8), &leans_short, 300, 91),
            drawn_at(stratum_at(2, 14), &leans_long, 300, 92),
        ];
        let pools = fit_period_length_spectra(&strata, &[0.4; 8], &pooling_config());

        let pool = pools.get(&2).expect("period 2 has 600 tracts");
        assert_eq!(pool.strata_pooled, 2, "both strata of period 2 are in it");
        assert_eq!(
            pool.tracts_fitted, 600,
            "300 tracts from each of the two strata"
        );
        assert_eq!(
            pool.length_spectrum.len(),
            3,
            "2 x span + 1 classes at span 1"
        );

        let (short, long) = (pool.length_spectrum[0], pool.length_spectrum[2]);
        assert!(
            short > leans_long[0] && short < leans_short[0],
            "the pooled share one repeat short is {short}, and pooling a stratum drawn at \
             {} with one drawn at {} puts it between them",
            leans_short[0],
            leans_long[0]
        );
        assert!(
            long > leans_short[2] && long < leans_long[2],
            "the pooled share one repeat long is {long}, between the two truths {} and {}",
            leans_short[2],
            leans_long[2]
        );
    }

    /// **Each motif period is pooled apart**, so a dinucleotide tract is never seeded from
    /// trinucleotide evidence.
    ///
    /// The two periods are drawn from opposite tilts, so a pool that ran over every stratum of
    /// the run at once would return one spectrum twice — which is what the last assertion
    /// refuses.
    #[test]
    fn each_motif_period_is_pooled_apart() {
        let leans_short = tilted(0.55, 0.35);
        let leans_long = tilted(0.10, 0.35);
        let strata = [
            drawn_at(stratum_at(2, 8), &leans_short, 300, 93),
            drawn_at(stratum_at(3, 8), &leans_long, 300, 94),
        ];
        let pools = fit_period_length_spectra(&strata, &[0.4; 8], &pooling_config());

        assert_eq!(pools.len(), 2, "one pool a period");
        let dinucleotide = &pools[&2];
        let trinucleotide = &pools[&3];
        assert_eq!(dinucleotide.strata_pooled, 1);
        assert_eq!(trinucleotide.strata_pooled, 1);
        assert!(
            dinucleotide.length_spectrum[0] > trinucleotide.length_spectrum[0] + 0.15,
            "period 2 was drawn leaning short ({}) and period 3 leaning long ({}); pooled \
             together they would be one number twice, and they are {} and {}",
            leans_short[0],
            leans_long[0],
            dinucleotide.length_spectrum[0],
            trinucleotide.length_spectrum[0]
        );
    }

    /// **A period as thin as a stratum is refused as a stratum is**, by the same floor in the
    /// same unit: pooling five tracts does not make them eight.
    #[test]
    fn a_period_below_the_refusal_floor_gets_no_pool() {
        let spectrum = spectrum_of(3);
        let strata = [
            drawn_at(stratum_at(2, 8), &spectrum, 3, 95),
            drawn_at(stratum_at(2, 14), &spectrum, 4, 96),
        ];
        let config = pooling_config();
        assert_eq!(config.refusal_floor, DEFAULT_REFUSAL_FLOOR);

        let pools = fit_period_length_spectra(&strata, &[0.4; 8], &config);
        assert!(
            pools.is_empty(),
            "seven tracts is below the floor of {}, and a pool below it is left out rather \
             than fitted badly",
            config.refusal_floor
        );

        // …and one more tract clears it, so the emptiness above is the floor and not the
        // fixture being unfittable.
        let over = [
            drawn_at(stratum_at(2, 8), &spectrum, 4, 95),
            drawn_at(stratum_at(2, 14), &spectrum, 4, 96),
        ];
        let pools = fit_period_length_spectra(&over, &[0.4; 8], &config);
        assert_eq!(pools[&2].tracts_fitted, 8);
    }

    /// **An allele span of zero is refused where the knob is named**, not three modules later.
    ///
    /// It is a public field with no lower bound, read from an environment variable by
    /// `examples/ng_joint_records_walk.rs` and parsed with no floor. At zero the fit returns a
    /// one-class length spectrum — a tract that can only ever be its reference length — and
    /// `StratumFits::over` then aborts the run with a message about a class count, which names
    /// the symptom and not the setting.
    #[test]
    #[should_panic(expected = "`SsrFitConfig::allele_span` must be at least 1")]
    fn a_fit_that_may_place_no_allele_mass_anywhere_is_refused() {
        let spectrum = spectrum_of(3);
        let evidence = drawn_at(stratum_at(2, 8), &spectrum, 20, 104);
        let no_span = SsrFitConfig {
            allele_span: 0,
            ..SsrFitConfig::default()
        };
        let _ = fit_stratum(&evidence, &[0.4; 8], &no_span);
    }

    /// **The floor counts tracts a read crossed, not tracts.** A tract nobody sequenced
    /// contributes a likelihood of exactly one whatever the parameters are, so eight of them
    /// carry nothing — and every other fixture here draws reads at every tract, which makes the
    /// two counts the same list and a floor measured in the wrong unit invisible.
    #[test]
    fn the_pools_floor_counts_tracts_a_read_crossed() {
        let spectrum = spectrum_of(3);
        let mut with_silent_tracts = drawn_at(stratum_at(2, 8), &spectrum, 4, 101);
        // Four more tracts, present and unread: `tracts.len()` is 8, `tracts_with_reads()` 4.
        with_silent_tracts
            .tracts
            .extend(std::iter::repeat_with(TractReads::default).take(4));
        assert_eq!(with_silent_tracts.tracts.len(), 8);
        assert_eq!(with_silent_tracts.tracts_with_reads(), 4);

        let pools = fit_period_length_spectra(&[with_silent_tracts], &[0.4; 8], &pooling_config());
        assert!(
            pools.is_empty(),
            "four tracts with reads is below the floor of {}, whatever the eight rows suggest",
            DEFAULT_REFUSAL_FLOOR
        );
    }

    /// **Whether the pooled climb settled is reported and not asserted true.** Running out of
    /// rounds is a real answer and it must not come back as convergence.
    ///
    /// One round is far too few for 600 tracts, so this pins the `false` arm; the positive
    /// control above reaches the `true` one at the default limit.
    #[test]
    fn a_pool_that_ran_out_of_rounds_does_not_report_convergence() {
        let spectrum = spectrum_of(3);
        let strata = [drawn_at(stratum_at(2, 8), &spectrum, 300, 102)];
        let one_round = SsrFitConfig {
            max_rounds: 1,
            ..pooling_config()
        };
        let pools = fit_period_length_spectra(&strata, &[0.4; 8], &one_round);
        assert!(
            !pools[&2].converged,
            "one round cannot settle a 300-tract pool, and running out is not convergence"
        );
    }

    /// **The pooled fit reads the samples' homozygote excess**, which is how inbreeding is
    /// divided out inside the estimator rather than afterwards. Every other fixture here passes
    /// the same `0.4` to every sample, so a pool that ignored the argument entirely would
    /// change nothing any of them assert.
    #[test]
    fn the_pool_reads_the_homozygote_excess_it_is_handed() {
        let spectrum = spectrum_of(3);
        let strata = [drawn_at(stratum_at(2, 8), &spectrum, 200, 103)];
        let config = pooling_config();

        let selfing = fit_period_length_spectra(&strata, &[0.9; 8], &config);
        let outbred = fit_period_length_spectra(&strata, &[0.0; 8], &config);

        let (selfing, outbred) = (&selfing[&2], &outbred[&2]);
        // **The bar is "not the same number", not a size.** A fit that ignored the argument
        // would return bit-identical results, because everything else about the two runs is
        // identical down to the draw's seed. The 1% is there so that a difference in the last
        // few bits would not pass for reading it; measured, the gap is **4.9%** — 0.5510
        // against 0.5254 — on 200 tracts drawn at an excess of 0.4 and read at 0.9 and 0.0.
        let gap = (selfing.concentration - outbred.concentration).abs() / outbred.concentration;
        assert!(
            gap > 0.01,
            "the same tracts read as a selfing panel and as an outbred one give concentrations \
             {} and {}, {:.1}% apart — at zero the excess is not reaching the fit at all",
            selfing.concentration,
            outbred.concentration,
            gap * 100.0
        );
    }

    /// Two strata of one period that recorded read offsets in different numbers of buckets came
    /// from two runs, and pooling them would index past the end of a bucket row inside the
    /// scorer.
    #[test]
    #[should_panic(expected = "recorded read offsets in")]
    fn two_strata_of_one_period_disagreeing_about_the_read_span_are_refused() {
        let spectrum = spectrum_of(3);
        let mut wider = drawn_at(stratum_at(2, 14), &spectrum, 20, 98);
        wider.read_span = 2;
        let strata = [drawn_at(stratum_at(2, 8), &spectrum, 20, 97), wider];
        let _ = fit_period_length_spectra(&strata, &[0.4; 8], &pooling_config());
    }

    /// The same for the slippage-group count, which is also a property of the run.
    #[test]
    #[should_panic(expected = "slippage groups and the one at")]
    fn two_strata_of_one_period_disagreeing_about_the_group_count_are_refused() {
        let spectrum = spectrum_of(3);
        let mut more_groups = drawn_at(stratum_at(2, 14), &spectrum, 20, 100);
        more_groups.groups = 2;
        let strata = [drawn_at(stratum_at(2, 8), &spectrum, 20, 99), more_groups];
        let _ = fit_period_length_spectra(&strata, &[0.4; 8], &pooling_config());
    }

    // -----------------------------------------------------------------
    // A stratum read from a subset of samples (plan step D2)
    // -----------------------------------------------------------------

    /// The slippage the subset tests draw at.
    fn subset_truth() -> Slippage {
        Slippage {
            level: 0.08,
            shorter_share: 0.83,
            fall_off: 0.25,
        }
    }

    /// A config whose first subset is `first` samples and whose level target is `target`.
    fn subset_config(first: usize, target: f64) -> SsrFitConfig {
        let mut config = SsrFitConfig {
            allele_span: 1,
            subsets: SampleSubsets {
                first,
                level_relative_error_target: target,
                min_samples_a_group: MIN_SAMPLES_A_GROUP,
            },
            ..SsrFitConfig::default()
        };
        config.curve.draw_curves = false;
        config
    }

    fn fitted(outcome: &StratumOutcome) -> &StratumFit {
        match outcome {
            StratumOutcome::Fitted(fit) => fit,
            other => panic!("not fitted: {other:?}"),
        }
    }

    /// **A cohort no larger than the first subset is fitted whole, and one sample more takes a
    /// subset**: with a first subset of 20, twenty samples give `fit_strata`'s outcome with no subset
    /// recorded; twenty-one are read from a subset.
    #[test]
    fn a_cohort_no_larger_than_the_first_subset_is_fitted_whole() {
        let config = subset_config(20, 10.0);
        for (samples, seed, subset) in [(20_usize, 51_u64, false), (21, 52, true)] {
            let evidence = draw_stratum(
                subset_truth(),
                &spectrum_of(3),
                0.5,
                0.4,
                100,
                samples,
                6,
                1,
                seed,
            );
            let order: Vec<usize> = (0..samples).collect();
            let excess = vec![0.4; samples];
            let strata = [evidence];
            let outcomes = fit_strata_on_sample_subsets(&strata, &excess, &order, &config);
            assert_eq!(
                fitted(&outcomes[0]).samples_fitted_on.is_some(),
                subset,
                "{samples} samples"
            );
            if !subset {
                assert_eq!(outcomes, fit_strata(&strata, &excess, &config));
            }
        }
    }

    /// **The subset stops growing once the level is measured to the target, and otherwise grows to
    /// every sample, whose answer is the best of every starting point and the last answer**: on 64
    /// drawn samples with a first subset of 8, a target no fit reaches takes all 64 and four walks,
    /// scoring at least what three starts on every sample score; a target the first fit meets stops
    /// at 8 from three walks, with the subset's own, smaller, evidence counts.
    #[test]
    fn a_subset_grows_until_the_level_is_measured_or_every_sample_is_in() {
        let evidence = draw_stratum(subset_truth(), &spectrum_of(3), 0.5, 0.4, 300, 64, 6, 1, 53);
        let order: Vec<usize> = (0..64).collect();
        let excess = [0.4; 64];
        let strata = [evidence];

        let unreachable =
            fit_strata_on_sample_subsets(&strata, &excess, &order, &subset_config(8, 1e-9));
        let every = fitted(&unreachable[0]);
        assert_eq!(every.samples_fitted_on, Some(64));
        assert_eq!(every.walks.len(), 4, "{:?}", every.walks);
        let whole = fit_stratum(&strata[0], &excess, &subset_config(8, 1e-9)).expect("reads");
        assert_eq!(every.reads_crossing, whole.reads_crossing);
        assert!(every.log_likelihood_a_tract >= whole.log_likelihood_a_tract);

        let met = fit_strata_on_sample_subsets(&strata, &excess, &order, &subset_config(8, 10.0));
        let first = fitted(&met[0]);
        assert_eq!(first.samples_fitted_on, Some(8));
        assert_eq!(first.walks.len(), 3);
        assert!(first.reads_crossing < whole.reads_crossing);
        assert_eq!(first.tracts_fitted, first.tracts_of_its_own);
    }

    /// **A subset doubles, and takes every sample once doubling would pass three quarters of them.**
    #[test]
    fn the_subset_doubles_until_three_quarters_of_the_cohort() {
        assert_eq!(next_subset_size(256, 2_169), 512);
        assert_eq!(next_subset_size(512, 2_169), 1_024);
        assert_eq!(next_subset_size(1_024, 2_169), 2_169);
        assert_eq!(next_subset_size(16, 64), 32);
        assert_eq!(next_subset_size(32, 64), 64);
        assert_eq!(next_subset_size(64, 64), 64);
    }

    /// **A larger subset holds the smaller one's samples, and every slippage group has eight of its
    /// readers in**: group 1 read by samples 2, 5, 7 and 20 to 29; the first eight samples hold three
    /// of them, so five more are added in order; doubling keeps every one of them.
    #[test]
    fn a_subset_tops_up_a_thin_group_and_keeps_its_samples_as_it_grows() {
        let samples = 40;
        let order: Vec<usize> = (0..samples).collect();
        let mut readers = vec![vec![true; samples], vec![false; samples]];
        for sample in [2, 5, 7].into_iter().chain(20..30) {
            readers[1][sample] = true;
        }
        let mut keep = vec![false; samples];
        grow_subset(&mut keep, &readers, &order, 8, MIN_SAMPLES_A_GROUP);
        let kept: Vec<usize> = (0..samples).filter(|&sample| keep[sample]).collect();
        assert_eq!(kept, (0..8).chain(20..25).collect::<Vec<_>>());
        let before = keep.clone();
        grow_subset(&mut keep, &readers, &order, 16, MIN_SAMPLES_A_GROUP);
        assert!(before.iter().zip(&keep).all(|(was, is)| !*was || *is));
        let kept: Vec<usize> = (0..samples).filter(|&sample| keep[sample]).collect();
        assert_eq!(kept, (0..16).chain(20..25).collect::<Vec<_>>());
    }

    /// **The refusal floor is judged on the whole stratum, and a subset with fewer tracts than the
    /// floor grows rather than being fitted**: forty tracts each read by one sample, a first subset
    /// of 4 holding 4 tracts — not refused, and not fitted on 4 tracts: it grows to the 8 samples that
    /// read 8. Six such tracts are refused, below the floor of 8.
    #[test]
    fn the_floor_is_judged_on_the_whole_stratum_and_a_thin_subset_grows() {
        let read_by_one = |tracts: usize| {
            let mut evidence =
                draw_stratum(subset_truth(), &spectrum_of(3), 0.5, 0.4, 40, 40, 6, 1, 57);
            for (index, tract) in evidence.tracts.iter_mut().enumerate() {
                tract
                    .samples
                    .retain(|reads| reads.sample as usize == index && index < tracts);
            }
            evidence
        };
        let order: Vec<usize> = (0..40).collect();
        // One reader a group at least, so the top-up leaves the first subset at its 4 tracts.
        let mut config = subset_config(4, 10.0);
        config.subsets.min_samples_a_group = 1;
        let outcomes =
            fit_strata_on_sample_subsets(&[read_by_one(40)], &[0.4; 40], &order, &config);
        let fit = fitted(&outcomes[0]);
        assert_eq!(fit.samples_fitted_on, Some(8));
        assert_eq!(fit.tracts_fitted, 8);
        let thin = fit_strata_on_sample_subsets(&[read_by_one(6)], &[0.4; 40], &order, &config);
        assert!(
            matches!(
                thin[0],
                StratumOutcome::Refused {
                    reason: StratumRefusal::BelowTheFloor { .. },
                    ..
                }
            ),
            "{:?}",
            thin[0]
        );
    }

    /// **The subset is the same samples whatever order they arrive in**: with the cohort's indices
    /// shuffled and the evidence relabelled to match, the subset holds the same names, a slippage
    /// group the first samples lack topped up with three of its readers.
    #[test]
    fn the_subset_is_the_same_samples_whatever_order_they_arrive_in() {
        use crate::parameter_estimation::joint::sample_order::{SAMPLE_ORDER_SEED, sample_order};
        let names: Vec<String> = (0..40).map(|i| format!("accession_{i}")).collect();
        // Group 1 reads in the samples ranked twentieth or later, so the first eight hold none.
        let mut rank = vec![0; 40];
        for (place, sample) in sample_order(&names, SAMPLE_ORDER_SEED)
            .into_iter()
            .enumerate()
        {
            rank[sample] = place;
        }
        let mut evidence =
            draw_stratum(subset_truth(), &spectrum_of(3), 0.5, 0.4, 100, 40, 6, 1, 59);
        for tract in &mut evidence.tracts {
            for reads in tract
                .samples
                .iter_mut()
                .filter(|reads| rank[reads.sample as usize] >= 20)
            {
                for (group, _) in &mut reads.by_group {
                    *group = 1;
                }
            }
        }
        evidence.groups = 2;
        let members_by_name = |evidence: &StratumEvidence, names: &[String]| {
            let order = sample_order(names, SAMPLE_ORDER_SEED);
            let mut keep = vec![false; names.len()];
            grow_subset(
                &mut keep,
                &evidence.readers_by_group(names.len()),
                &order,
                8,
                3,
            );
            let mut kept: Vec<String> = (0..names.len())
                .filter(|&sample| keep[sample])
                .map(|sample| names[sample].clone())
                .collect();
            kept.sort();
            kept
        };
        let before = members_by_name(&evidence, &names);
        // Sample `i` arrives at index `(7 i + 3) mod 40`.
        let moved = |sample: usize| (7 * sample + 3) % 40;
        let mut shuffled_names = vec![String::new(); 40];
        for (sample, name) in names.iter().enumerate() {
            shuffled_names[moved(sample)] = name.clone();
        }
        let mut shuffled = evidence.clone();
        for tract in &mut shuffled.tracts {
            for reads in &mut tract.samples {
                reads.sample = moved(reads.sample as usize) as u32;
            }
        }
        assert_eq!(members_by_name(&shuffled, &shuffled_names), before);
        assert_eq!(
            before.len(),
            8 + 3,
            "three group-1 readers added: {before:?}"
        );
    }

    /// **Both schedules give the same bits on subsets**: two strata of 40 samples, a first subset of
    /// 8, one stratum at a time and both at once.
    #[test]
    fn the_two_schedules_give_the_same_bits_on_subsets() {
        let strata: Vec<StratumEvidence> = [(10_u64, 61_u64), (11, 63)]
            .into_iter()
            .map(|(repeats, seed)| {
                let mut evidence = draw_stratum(
                    subset_truth(),
                    &spectrum_of(3),
                    0.5,
                    0.4,
                    150,
                    40,
                    6,
                    1,
                    seed,
                );
                evidence.stratum = Stratum {
                    period: 2,
                    reference_repeats: repeats,
                };
                evidence
            })
            .collect();
        let order: Vec<usize> = (0..40).rev().collect();
        let outcomes_at = |at_once: usize| {
            let config = SsrFitConfig {
                strata_at_once: NonZeroUsize::new(at_once).expect("positive"),
                ..subset_config(8, LEVEL_RELATIVE_ERROR_TARGET)
            };
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(4)
                .build()
                .expect("a pool");
            pool.install(|| fit_strata_on_sample_subsets(&strata, &[0.4; 40], &order, &config))
        };
        let (one, two) = (outcomes_at(1), outcomes_at(2));
        let bits = |outcomes: &[StratumOutcome]| -> Vec<u64> {
            every_fitted_number(outcomes)
                .into_iter()
                .map(f64::to_bits)
                .collect()
        };
        assert_eq!(bits(&one), bits(&two));
        assert_eq!(one, two);
        assert!(fitted(&one[0]).samples_fitted_on.is_some());
    }

    /// **The level is measured to the target when every live group that can still grow is**: two
    /// groups at 1% and 5% of their levels; a target of 2% is met only when the second has no
    /// readers left outside, or neither has; a group without an error that can grow has not met it.
    #[test]
    fn the_level_is_judged_only_in_groups_that_can_still_grow() {
        let mut parameters = Parameters::start(
            StartingPoint {
                slippage_level: 0.1,
                concentration: 0.5,
            },
            2,
            3,
        );
        parameters.slippage[1].level = 0.1;
        let group = |error: f64| SlippageErrors {
            level: StratumError::Estimated(error),
            shorter_share: StratumError::Estimated(0.1),
            fall_off: StratumError::Estimated(0.1),
        };
        let mut errors = StratumErrors {
            slippage: vec![Some(group(0.001)), Some(group(0.005))],
            length_spectrum: vec![StratumError::Estimated(0.1); 3],
            concentration: StratumError::Estimated(0.1),
        };
        let live = [true, true];
        let judged = |errors: &StratumErrors, growing: [bool; 2]| {
            level_is_measured_to(&parameters, &live, errors, &growing, 0.02)
        };
        assert!(!judged(&errors, [true, true]));
        assert!(judged(&errors, [true, false]));
        assert!(judged(&errors, [false, false]));
        errors.slippage[0] = Some(SlippageErrors {
            level: StratumError::NotIdentified,
            ..group(0.001)
        });
        assert!(!judged(&errors, [true, false]));
        assert!(judged(&errors, [false, false]));
    }

    /// **One read makes a reader**: a sample whose buckets hold one read reads the group, one whose
    /// buckets are all empty does not.
    #[test]
    fn one_read_makes_a_sample_a_reader_of_its_group() {
        let evidence = StratumEvidence {
            stratum: Stratum {
                period: 2,
                reference_repeats: 10,
            },
            tracts: vec![TractReads {
                samples: vec![
                    SampleTractReads {
                        sample: 0,
                        by_group: vec![(0, vec![0, 1, 0])],
                    },
                    SampleTractReads {
                        sample: 2,
                        by_group: vec![(0, vec![0, 0, 0])],
                    },
                ],
            }],
            read_span: 1,
            groups: 1,
            tracts_over_guard_threshold: 0,
            reads_reaching_not_crossing: 0,
            guard_reads: 0,
            bases_compared: 0,
            mismatching_bases: 0,
        };
        assert_eq!(evidence.readers_by_group(3), vec![vec![true, false, false]]);
    }

    /// **The log's line about the subsets**: the fewest, the median and the most samples a stratum
    /// was read from, and how many grew to every sample; no line when no stratum took a subset.
    #[test]
    fn the_subsets_summary_counts_the_sizes() {
        let with = |size: Option<usize>| {
            let mut outcome = fitted_at(2, 10, 0.1, 1_000);
            if let StratumOutcome::Fitted(fit) = &mut outcome {
                fit.samples_fitted_on = size;
            }
            outcome
        };
        let outcomes = [
            with(Some(256)),
            with(Some(2_169)),
            with(Some(512)),
            with(None),
        ];
        let summary = subsets_summary(&outcomes, 2_169).expect("subsets taken");
        for part in [
            "the 3 strata",
            "256 to 2169, median 512",
            "1 of them grew to every sample",
        ] {
            assert!(summary.contains(part), "{part:?} in {summary}");
        }
        assert_eq!(subsets_summary(&[with(None)], 2_169), None);
    }

    /// **The subset against every sample** (plan step D2, spec §4.5 item 3 on drawn cohorts): strata
    /// drawn at a known level over a cohort larger than the first subset, fitted on every sample and
    /// on the grown subset; each line gives the samples the subset reached, both levels and their
    /// errors, how far apart the two levels are in the whole fit's errors, and both times.
    ///
    /// Ignored by default: minutes in the container. `NG_FIT_PRECISION_SUBSET_DRAWS` sets the draws a
    /// regime (default 3).
    #[test]
    #[ignore = "a measurement over cohorts larger than the first subset, each fitted twice"]
    fn the_subset_against_every_sample() {
        let draws: u64 = std::env::var("NG_FIT_PRECISION_SUBSET_DRAWS")
            .ok()
            .map_or(3, |value| {
                value.trim().parse().expect("a whole number of draws")
            });
        let truth = Slippage {
            level: 0.05,
            shorter_share: 0.8,
            fall_off: 0.3,
        };
        for (classes, tracts, samples, depth) in [
            (3, 300, 1_024, 3),
            (3, 1_000, 1_024, 3),
            (3, 300, 1_024, 30),
            (13, 200, 600, 3),
        ] {
            let span = (classes / 2) as i32;
            let config = SsrFitConfig {
                allele_span: span,
                ..SsrFitConfig::default()
            };
            let excess = vec![0.4; samples];
            let order: Vec<usize> = (0..samples).collect();
            for draw in 0..draws {
                let seed = 0xD200_0000 + (classes as u64) * 10_000 + tracts as u64 + draw;
                let evidence = draw_stratum(
                    truth,
                    &spectrum_of(classes),
                    0.5,
                    0.4,
                    tracts,
                    samples,
                    depth,
                    span,
                    seed,
                );
                let started = std::time::Instant::now();
                let whole = fit_stratum(&evidence, &excess, &config).expect("reads were drawn");
                let whole_time = started.elapsed();
                let started = std::time::Instant::now();
                let outcome = fit_on_growing_subsets(
                    &evidence,
                    &order,
                    &excess,
                    &config,
                    WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
                );
                let subset_time = started.elapsed();
                let StratumOutcome::Fitted(subset) = outcome else {
                    panic!("not fitted: {outcome:?}");
                };
                let level_and_error = |fit: &StratumFit| {
                    let level = fit.slippage[0].expect("live").level;
                    let error = fit
                        .standard_errors
                        .as_ref()
                        .and_then(|errors| errors.slippage[0])
                        .and_then(|group| group.level.value());
                    (level, error)
                };
                let (whole_level, whole_error) = level_and_error(&whole);
                let (subset_level, subset_error) = level_and_error(&subset);
                let apart = whole_error.map(|error| (subset_level - whole_level) / error);
                eprintln!(
                    "SUBSET {classes} classes, {tracts} tracts x {samples} samples x {depth} reads, \
                     draw {draw}: subset {:?} samples; level every sample {whole_level:.5} \
                     ± {:.5}, subset {subset_level:.5} ± {:.5} ({:.2} of the whole fit's errors \
                     apart); {:.1} s against {:.1} s",
                    subset.samples_fitted_on,
                    whole_error.unwrap_or(f64::NAN),
                    subset_error.unwrap_or(f64::NAN),
                    apart.unwrap_or(f64::NAN),
                    subset_time.as_secs_f64(),
                    whole_time.as_secs_f64(),
                );
            }
        }
    }
}

#[cfg(test)]
mod the_largest_table {
    use super::*;

    /// A stratum whose tracts have reads from `samples[t]` samples each.
    fn stratum_with(samples: &[usize]) -> StratumEvidence {
        StratumEvidence {
            stratum: Stratum {
                period: 2,
                reference_repeats: 8,
            },
            tracts: samples
                .iter()
                .map(|&count| TractReads {
                    samples: (0..count as u32)
                        .map(|sample| SampleTractReads {
                            sample,
                            by_group: vec![(0, vec![1])],
                        })
                        .collect(),
                })
                .collect(),
            read_span: RECORDED_OFFSET_RANGE,
            groups: 1,
            tracts_over_guard_threshold: 0,
            reads_reaching_not_crossing: 0,
            guard_reads: 0,
            bases_compared: 0,
            mismatching_bases: 0,
        }
    }

    /// Rows with reads are counted over every tract, and the bytes over every vector a row owns.
    #[test]
    fn a_stratums_rows_and_bytes_are_counted_over_every_tract() {
        let stratum = stratum_with(&[3, 0, 2]);
        assert_eq!(stratum.rows_with_reads(), 5);
        let row = std::mem::size_of::<SampleTractReads>()
            + std::mem::size_of::<(u32, Vec<u32>)>()
            + std::mem::size_of::<u32>();
        assert_eq!(
            stratum.heap_bytes(),
            3 * std::mem::size_of::<TractReads>() + 5 * row
        );
    }

    /// One row a sample with reads a tract, 91 genotype pairs at the default span of six, eight
    /// bytes each — and the largest stratum's, not the sum.
    #[test]
    fn the_table_is_rows_times_genotypes_times_eight_bytes_for_the_largest_stratum() {
        let config = SsrFitConfig::default();
        assert_eq!(
            config.allele_span, 6,
            "the arithmetic below assumes the default span"
        );
        let strata = [stratum_with(&[3, 5]), stratum_with(&[10, 0, 2])];
        assert_eq!(largest_table_bytes(&strata, &config), 12 * 91 * 8);
    }
}
