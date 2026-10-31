//! The repeat-tract half of the joint fit: how often a polymerase slips, per stratum.
//!
//! While a polymerase copies a repeat tract it sometimes **slips**, adding or dropping a whole
//! repeat unit, so a read reports a tract one unit longer or shorter than the DNA it came
//! from. This module estimates how often that happens, from every sample's reads at the same
//! tracts at once. It is not the variant caller: it runs once over the cohort, before anything
//! is called, and hands the caller the numbers it will assume.
//!
//! Design: [`spec/parameter_prepass_joint_fit.md`] §4 for the evidence and what is handed on,
//! `spec/str_penalised_fit.md` for the fit. The
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
//! - **Allele length variation** — how the stratum's chromosomes are spread over the tract lengths.
//! - **Concentration** — how monomorphic the stratum's tracts are. Small means most tracts are
//!   fixed at one length while the stratum as a whole spans many.
//!
//! # What is fitted
//!
//! Per (slippage group × stratum), three numbers describing slippage: **how often** a read slips,
//! **which way** — at tomato dinucleotides a slipped read shows a shorter tract 4.9 times as
//! often as a longer one — and **how fast** two-unit slips fall off against one-unit slips.
//! Per stratum, the allele length variation and the concentration. **None of them is per tract.**
//!
//! # Where the fit lives
//!
//! This file holds the vocabulary every repeat-tract fit shares: the evidence
//! ([`StratumEvidence`], gathered by [`gather_strata`]), what a fit returns ([`StratumOutcome`],
//! [`StratumFit`], [`StratumErrors`]) and where each number came from ([`LevelProvenance`],
//! [`SharesProvenance`]). The fit itself is [`penalised_fit::fit_repeat_tracts`]: every class of
//! one motif length fitted at once by penalised expectation maximisation (spec
//! `str_penalised_fit.md`). Its first step and the likelihood of pairs of samples are in
//! [`expected_counts`]; the counting fit there is the zero-penalty answer its tests compare
//! against.
//!
//! **The fit before it, the climb, was deleted on 2026-10-16** (`arch/str_penalised_fit.md`
//! §4.1): a coordinate climb per stratum over a 256-point average of each tract's likelihood,
//! with curves drawn across repeat count afterwards; the code that drew and blended those curves
//! went later the same day. What older parameters files carry from it —
//! curve-sourced and blended numbers ([`SlippageCurve`], [`ShareCurve`]), strata whose numbers
//! all came from those curves, and period-pooled allele length variations — is still read, by
//! [`StratumFits::of_gathered_rows`](super::stratum_fits::StratumFits::of_gathered_rows), which
//! takes the file's rows and never a [`StratumOutcome`].
//!
//! [`spec/parameter_prepass_joint_fit.md`]: ../../../../doc/devel/ng/spec/parameter_prepass_joint_fit.md

use std::collections::BTreeMap;
use std::sync::Arc;

use rayon::prelude::*;

pub mod expected_counts;
pub mod penalised_fit;

use crate::float;
use crate::parameter_estimation::joint::census::{
    CensusError, CohortCensusEvidence, GuardTally, RECORDED_OFFSET_RANGE, SsrEvidence,
    SsrLocusState,
};
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
pub use super::share_curve::{ShareCurve, ShareCurveSource, ShareShape, ShareSource};
pub use super::slippage_curve::{CurveReach, LevelSource, SlippageCurve};

/// **The highest slippage level a library's multiplier may lift a read to**
/// (`str_library_stutter.md` §2.1). A level of one would say no read reports its allele.
///
/// **It caps the lift, never the class's own level**: a library is scored at
/// `min(m · L, max(L, MULTIPLIED_LEVEL_CEILING))`, so a multiplier of one scores every read at
/// exactly the class rate whatever that rate is, and the fit with every multiplier one is today's
/// to the bit (spec §3.1). Named apart from the curves' `LEVEL_CEILING` (0.999), which bounds a
/// level read off a curve and is a different number.
pub const MULTIPLIED_LEVEL_CEILING: f64 = 0.95;

/// **The slippage level a library's reads are scored at**: the class rate `class_level` times the
/// library's `multiplier`, held below [`MULTIPLIED_LEVEL_CEILING`] unless the class rate itself is
/// above it.
///
/// **A multiplier of one returns the class rate itself**, whatever it is — a rate that is not a
/// number included, which the `min` and `max` below would otherwise turn into the ceiling.
///
/// # Panics
///
/// In a debug build, when `multiplier` is not a finite positive number.
#[must_use]
pub fn library_level(class_level: f64, multiplier: f64) -> f64 {
    debug_assert!(
        multiplier.is_finite() && multiplier > 0.0,
        "a library's stutter multiplier is a finite positive number, not {multiplier}"
    );
    if multiplier == 1.0 {
        return class_level;
    }
    (multiplier * class_level).min(MULTIPLIED_LEVEL_CEILING.max(class_level))
}

/// **Which library of [`StratumEvidence::libraries`] a sample's reads came through** — its position
/// in that list. A type of its own so that it cannot be mistaken for the sample or the slippage
/// group, which are numbered the same way.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LibraryIndex(pub u32);

impl LibraryIndex {
    /// The position in [`StratumEvidence::libraries`].
    #[must_use]
    pub fn get(self) -> usize {
        self.0 as usize
    }
}

/// **One library of the cohort** — a read group, the grain the census keeps apart — and how its
/// reads are scored (`str_library_stutter.md` §1.2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Library {
    /// The read group the library is.
    pub read_group: ReadGroupId,
    /// The slippage group whose class rates the library's reads are scored at.
    pub slippage_group: u32,
    /// The library's stutter relative to its slippage group's class rates: a read slips with
    /// probability [`library_level`]`(class rate, stutter_multiplier)`. **Nothing in a run sets it
    /// now**: the penalised fit holds every library at one, and the fit that estimated it was
    /// deleted with the climb on 2026-10-16. Drawn cohorts of the fitting bench set it to give
    /// libraries a known stutter of their own.
    pub stutter_multiplier: f64,
}

impl Library {
    /// A library of `slippage_group` with a multiplier of one.
    #[must_use]
    pub fn in_group(read_group: ReadGroupId, slippage_group: u32) -> Self {
        Self {
            read_group,
            slippage_group,
            stutter_multiplier: 1.0,
        }
    }
}

/// One sample's spanning reads at one tract, split by the library that produced them.
///
/// **Only the samples that put a read on the tract are here.** A sample with no read
/// contributes a likelihood of exactly one whatever the parameters are, so leaving it out is
/// not an approximation — and at three reads a site it is most of the panel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleTractReads {
    pub sample: u32,
    /// Per library that put a read here, in increasing library order, which library it is and its
    /// reads in each bucket.
    pub by_library: Vec<(LibraryIndex, Vec<u32>)>,
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
    /// **Every library of the cohort**, indexed as [`SampleTractReads::by_library`] names them:
    /// which slippage group each reads and the multiplier its reads are scored at.
    ///
    /// **One list a run, shared by every stratum** rather than copied into each: the multipliers
    /// are the run's, and a fit that updates them builds one new list and hands the same one to
    /// every stratum, so no stratum can be left scoring at an old multiplier. Every library's
    /// slippage group is below [`groups`](Self::groups) ([`check_libraries`](Self::check_libraries)).
    pub libraries: Arc<[Library]>,
    /// Tracts left out because, over every sample and read group, more than one read in ten of
    /// those differing from the reference length differed by a **non-whole** number of repeat
    /// units — the guard's threshold. Such a tract is not something this noise model describes,
    /// and the fit says so rather than fitting it (`parameter_prepass_joint_records.md` §3.3).
    pub tracts_over_guard_threshold: u64,
    /// Reads that reached a tract and crossed no whole copy of it, so reported no length.
    /// **Not evidence about slippage and not dropped silently**: a tract longer than a read is
    /// never crossed in any sample at any depth, so this count runs along the repeat-count
    /// axis and a stratum unreadable at this read length must not look like one that was
    /// merely unlucky with coverage.
    pub reads_reaching_not_crossing: u64,
    /// Reads whose tract differed from the reference by a non-whole number of repeat units, at
    /// every tract of the stratum, those the guard left out included. A diagnostic; nothing
    /// about slippage is estimated from them.
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
    /// the tract list, each tract's rows, each row's libraries and each library's counts. The
    /// allocator rounds every block up, so the process holds somewhat more. **The library list is
    /// not counted**: one list is shared by every stratum of a run.
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
                                row.by_library.capacity()
                                    * std::mem::size_of::<(LibraryIndex, Vec<u32>)>()
                                    + row
                                        .by_library
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

    /// Tracts carrying at least one spanning read — the count the penalised fit spreads its
    /// allele-length prior over, and the one the bench's
    /// `COUNTING_FIT_FLOOR` is measured in.
    pub fn tracts_with_reads(&self) -> usize {
        self.tracts.iter().filter(|t| !t.samples.is_empty()).count()
    }

    /// Spanning reads, over every tract and sample.
    pub fn spanning_reads(&self) -> u64 {
        self.tracts
            .iter()
            .flat_map(|tract| tract.samples.iter())
            .flat_map(|sample| sample.by_library.iter())
            .flat_map(|(_, counts)| counts.iter())
            .map(|reads| u64::from(*reads))
            .sum()
    }

    #[cfg(test)]
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

    /// **Refuses evidence whose libraries do not fit together**: a library in a slippage group at or
    /// past [`groups`](Self::groups), or a sample's reads naming a library the list does not hold or
    /// naming libraries out of increasing order.
    ///
    /// Called where a fit first takes a stratum, so a builder that sets one of the two fields and not
    /// the other is stopped there, by name, rather than as an index out of bounds deep in scoring.
    ///
    /// # Panics
    ///
    /// On any of the three.
    pub fn check_libraries(&self) {
        if let Some(library) = self
            .libraries
            .iter()
            .find(|library| library.slippage_group as usize >= self.groups)
        {
            panic!(
                "period {} at {} repeats: read group {:?} is in slippage group {} of the {} the \
                 evidence declares",
                self.stratum.period,
                self.stratum.reference_repeats,
                library.read_group,
                library.slippage_group,
                self.groups
            );
        }
        for row in self.tracts.iter().flat_map(|tract| &tract.samples) {
            let in_order = row.by_library.windows(2).all(|pair| pair[0].0 < pair[1].0);
            let in_range = row
                .by_library
                .last()
                .is_none_or(|(library, _)| library.get() < self.libraries.len());
            assert!(
                in_order && in_range,
                "period {} at {} repeats: sample {}'s reads name libraries {:?}, out of order or \
                 past the {} the evidence lists",
                self.stratum.period,
                self.stratum.reference_repeats,
                row.sample,
                row.by_library
                    .iter()
                    .map(|(library, _)| library.0)
                    .collect::<Vec<_>>(),
                self.libraries.len()
            );
        }
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
    pub allele_length_variation: Vec<f64>,
    /// How monomorphic the stratum's tracts are. Small means most tracts carry one length.
    ///
    /// **`None` where the fit could not measure it** (`str_penalised_fit.md` §2.5). The penalised
    /// fit reads the concentration from pairs of samples at one locus, so it is absent in two
    /// cases:
    /// - the class has no locus holding two samples with reads;
    /// - the estimate stopped within 1% of an end of its search range, so the data did not place
    ///   it.
    ///
    /// An absent concentration is not a fitted zero.
    pub concentration: Option<f64>,
    /// **How many samples the stratum was fitted on**, when it was fitted on a subset of the
    /// cohort; `None` when every sample was read. The penalised fit reads every sample, so it
    /// always writes `None`. The fit before it (the climb, deleted 2026-10-16) fitted a large
    /// cohort on a subset, and a parameters file it wrote still carries the count in its rows'
    /// `samples_fitted_on`.
    pub samples_fitted_on: Option<usize>,
    /// Tracts **this stratum itself** holds with at least one spanning read **among the samples it
    /// was fitted on** — every sample unless it was fitted on a subset
    /// ([`samples_fitted_on`](Self::samples_fitted_on)).
    pub tracts_of_its_own: usize,
    /// Reads that crossed a whole tract of **this stratum itself**, over every group and every
    /// sample it was fitted on ([`samples_fitted_on`](Self::samples_fitted_on)).
    pub reads_crossing: u64,
    /// Per slippage group, where that group's emitted slippage *level* came from — `None` where
    /// the group put no read in this stratum, matching [`StratumFit::slippage`] index for index.
    ///
    /// **The level's only**: the direction split and the fall-off have their own provenance in
    /// [`StratumFit::shares_provenance`]. From the penalised fit, a point on the class's motif
    /// length's penalised curve.
    pub level_provenance: Vec<Option<LevelProvenance>>,
    /// Per slippage group, where that group's direction split and fall-off came from.
    pub shares_provenance: Vec<Option<SharesProvenance>>,
    /// **How precisely the stratum's own reads determine each of its numbers**, at the fitted answer
    /// (`str_penalised_fit.md` §2.6).
    ///
    /// **The errors of the stratum's own fit, and only of it.** A number an older file took from a
    /// parametric curve, or blended with one, has none (`str_slippage_level_curve.md` §5.2). The
    /// own fit's level is the one [`LevelProvenance::slipped_reads`] is counted from.
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
    /// [`StratumFit::allele_length_variation`] is.
    pub allele_length_variation: Vec<StratumError>,
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
    ///
    /// **Only the counting fit's research errors produce it** (`counting_fit_errors`, which the
    /// fitting bench and two examples call). The penalised fit, the one a run uses, never does:
    /// it reports that case as [`StratumError::NotPlaced`] (the mapping in
    /// `doc/devel/ng/arch/str_penalised_fit.md` §1.4), so no parameters file or run log carries it.
    NotIdentified,
    /// **The tracts do not place it**: its error on the scale the number is fitted on — logit for a
    /// slippage number, log for the concentration, a log-ratio for a share — is above
    /// [`NOT_PLACED`], so a one-error interval spans odds or a size hundreds of times apart. Where
    /// a fit ran a number off towards an end of its range, the log-likelihood there is flat, and
    /// carrying its curvature to the number's own scale would report an error hundreds of times
    /// too small. Also a number in `[0, 1]` — a slippage number or
    /// a share — whose error, on its own scale, is wider than that whole range: the SNP/indel fit's
    /// "wider than its range" (spec §3.2).
    NotPlaced,
    /// **A length class with no share at all**, held at zero and never fitted.
    NoShare,
    /// **Withheld: too few loci for the error to be trusted.** Two rules, for two different
    /// reasons (`str_penalised_fit.md` §2.6, §7):
    ///
    /// - **A class's numbers below a count of loci with reads**: its slip rate, F and
    ///   allele-length shares below 8 loci, its shorter share and fall-off below 3. Measured on
    ///   drawn cohorts of 2,000 samples at 3 to 7 loci, the truth lay within two of the direction
    ///   numbers' errors 89 to 93 times in 100, against 95 for an honest error (89 to 90 over
    ///   every cohort size from 1 to 2,000 samples, 82 in the lowest cell; spec §2.6); within
    ///   two of the slip rate's only 72 to 79 times in 100, because the rate itself comes out low
    ///   on so few loci. Below 3 no number's error is reported.
    /// - **A share whose expected count of loci — the share times the class's loci with reads —
    ///   is below 1.** These errors were measured, and they were too small.
    TooFewLoci,
}

impl StratumErrors {
    #[cfg(test)]
    /// Whether any number of the stratum has an error.
    pub fn any(&self) -> bool {
        let slippage = self.slippage.iter().flatten().any(|group| {
            [group.level, group.shorter_share, group.fall_off]
                .iter()
                .any(|error| error.value().is_some())
        });
        slippage
            || self
                .allele_length_variation
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
            Self::NotIdentified | Self::NotPlaced | Self::NoShare | Self::TooFewLoci => None,
        }
    }
}

/// **The largest error, on the scale a number is fitted on, that still places the number**: three
/// units of logit or of log, so a one-error interval either side spans odds or a size `e⁶`, about
/// 400 times, apart. Measured on drawn strata with the fit before the penalised one (the climb,
/// deleted 2026-10-16; `fit_precision.md`, the review of its step C1): numbers the tracts
/// determined sat at 0.02 to 0.5 on those scales, and numbers the climb ran off to the end of its
/// reach at 1,000 to 2,400.
pub const NOT_PLACED: f64 = 3.0;

/// Where one slippage group's level at one stratum came from, and what stood behind it.
///
/// **This replaces what `Provenance::Borrowed` used to mark.** After the level becomes a curve, a
/// value fitted from 8,000 slipped reads and one interpolated across a gap look identical in the
/// number alone, and the mechanism that used to distinguish them no longer sets the level
/// (`str_slippage_level_curve.md` §8).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LevelProvenance {
    /// The cell's own fit, the curve, a blend — and for a blend, the share the curve carried —
    /// or a point on a penalised curve.
    pub source: LevelSource,
    /// The curve that supplied it, absent when this stratum's period had none, and absent for a
    /// point on a penalised curve, which is not a parametric curve. It carries the curve's own
    /// held-out error and how many cells stood behind it, so a consumer can tell a curve through
    /// twenty-three cells from one through four.
    pub curve: Option<SlippageCurve>,
    /// Whether this stratum's repeat count sat inside the curve's fitted range. `None` where
    /// there is no curve; always inside for a point on a penalised curve, whose range is the
    /// classes with reads. **A level held at a fitted end is under-stated in a known direction**
    /// (`str_slippage_level_curve.md` §6).
    pub reach: Option<CurveReach>,
    /// How many of this stratum's own reads slipped: from the penalised fit, the class's expected
    /// slipped reads at the fit's last E-step
    /// ([`PenalisedClass::expected_slipped_reads`](penalised_fit::PenalisedClass::expected_slipped_reads)).
    ///
    /// **The stratum's own reads, not the emitted level times them**, because this is the evidence
    /// that stood behind the cell — it is what set how precisely the stratum could determine its
    /// own answer. The emitted level is a point on a curve, so for a thin class it is mostly its
    /// neighbours'.
    ///
    /// **Absent is not zero.** In a file the climb wrote, `None` where the stratum borrowed its
    /// level: it has reads of its own — they are in [`StratumFit::reads_crossing`] — but no level
    /// of its own to say how many of them slipped.
    pub slipped_reads: Option<f64>,
}

/// Where a stratum's direction split and fall-off came from, and what stood behind them.
///
/// **From the penalised fit, each share is a point on its own penalised curve** across the
/// motif length's repeat counts, as the level is, each curve with its own stiffness. A file the
/// climb wrote may instead carry the stratum's own share, its period's parametric curve, or a
/// blend of the two (`str_slippage_level_curve.md` §5.1).
///
/// **The climb's curves replaced a gate with a cliff.** A stratum used to keep its own two shares
/// above 4,000 slipped reads and take one named neighbour's whole below it. Across both cohorts
/// one motif period of twelve ever cleared that floor, so 69 of HG002's strata and every one of
/// tomato's got nothing at all.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SharesProvenance {
    /// How many of this stratum's own reads slipped: from the penalised fit, the class's expected
    /// slipped reads at the fit's last E-step
    /// ([`PenalisedClass::expected_slipped_reads`](penalised_fit::PenalisedClass::expected_slipped_reads)),
    /// the same count as [`LevelProvenance::slipped_reads`]. In a file the climb wrote, `None`
    /// where nothing was fitted here.
    ///
    /// **Both shares are proportions over the reads that slipped**, so this one count sets how
    /// precisely the stratum holds either of them. It is the stratum's own reads rather than the
    /// emitted level times them: the emitted level is a point on a curve, partly its neighbours'.
    pub slipped_reads: Option<f64>,
    /// Where the share of slipped reads showing a *shorter* tract came from.
    pub shorter_share: ShareProvenance,
    /// Where the fall-off — how much rarer a two-unit slip is than a one-unit slip — came from.
    pub fall_off: ShareProvenance,
}

/// Where one of the two shares came from: the stratum's own fit, its period's curve, or a blend.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShareProvenance {
    /// The stratum's own fit, the curve, a blend — and for a blend the share the curve carried —
    /// or a point on a penalised curve.
    pub source: ShareSource,
    /// The curve that supplied it, absent where this stratum's period had none or the share is a
    /// point on a penalised curve. It carries its
    /// own held-out error, how many strata stood behind it, and which rung of the fallback
    /// ladder produced it.
    pub curve: Option<ShareCurve>,
    /// Whether this stratum's repeat count sat inside the curve's fitted range; `None` where
    /// there is no curve. **A share held at a fitted end is the end stratum's answer**, not this
    /// stratum's.
    pub reach: Option<CurveReach>,
}

/// Why a stratum produced no fit.
///
/// **One reason is left**: the climb also refused strata below its floor of tracts, and the
/// penalised fit has no floor, so a refused stratum now always has no tract with reads. Nothing
/// reads the reason or the count; both stay so that a refusal printed for debugging says why.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StratumRefusal {
    /// Not one tract in the stratum carried a spanning read.
    NoSpanningReads,
}

/// One stratum's answer: what was fitted, or why nothing was.
#[derive(Debug, Clone, PartialEq)]
pub enum StratumOutcome {
    /// Something was fitted from this stratum's own tracts.
    Fitted(Box<StratumFit>),
    /// There is no answer, and saying so is the point.
    Refused {
        stratum: Stratum,
        tracts: usize,
        reason: StratumRefusal,
    },
}

impl StratumOutcome {
    /// Which stratum this is about, whichever of the two it is.
    pub fn stratum(&self) -> Stratum {
        match self {
            Self::Fitted(fit) => fit.stratum,
            Self::Refused { stratum, .. } => *stratum,
        }
    }

    /// The slippage numbers a consumer would use, per slippage group — empty when there are
    /// none.
    pub fn slippage(&self) -> &[Option<Slippage>] {
        match self {
            Self::Fitted(fit) => &fit.slippage,
            Self::Refused { .. } => &[],
        }
    }

    /// Where each group's level came from — empty for a refusal.
    pub fn level_provenance(&self) -> &[Option<LevelProvenance>] {
        match self {
            Self::Fitted(fit) => &fit.level_provenance,
            Self::Refused { .. } => &[],
        }
    }

    /// Where each group's two shares came from — empty for a refusal.
    pub fn shares_provenance(&self) -> &[Option<SharesProvenance>] {
        match self {
            Self::Fitted(fit) => &fit.shares_provenance,
            Self::Refused { .. } => &[],
        }
    }
}

// ---------------------------------------------------------------------
// What the run was asked for
// ---------------------------------------------------------------------

/// How far either side of the reference length the fit may place allele mass: **the recorded
/// offsets' range plus one**, so ±9 while reads are recorded to ±8.
///
/// **Wider than the recorded range**, and that is what lets a read in an end bucket be attributed
/// to a distant allele rather than to a far slip (`parameter_prepass_joint_records.md` §3.2). A
/// window narrower than the recording forces every allele past it into its outermost class and
/// counts their reads as stutter (`str_allele_window.md` §1.1).
///
/// **No wider than the recorded range plus one**, because the reads can barely tell apart
/// alleles past it: an unslipped read from 9 repeats short and one from 15 repeats short both
/// land in the −8 bucket, and only reads that slipped two or more repeats longer differ. So **the outermost classes mean "this far or further"**: −9 holds every allele 9 or
/// more repeats shorter, +9 every one 9 or more longer, and calling shares their weight among
/// the candidates there (`str_allele_window.md` §1.2).
///
/// Tied to [`RECORDED_OFFSET_RANGE`] rather than written as a number, so that widening the
/// recording widens this too. It was a number once, and the recording grew from ±4 to ±8 while
/// this stayed at ±6.
pub const ALLELE_SPAN: i32 = RECORDED_OFFSET_RANGE + 1;

// The guard the paragraph above relies on, checked when the crate compiles.
const _: () = assert!(
    ALLELE_SPAN > RECORDED_OFFSET_RANGE,
    "the allele window must reach past the recorded offsets"
);

/// Scales `weights` to sum to one. Only the bench's cohort draws and the tests use it.
#[cfg(test)]
fn normalise(weights: &mut [f64]) {
    let total: f64 = weights.iter().sum();
    for weight in weights.iter_mut() {
        *weight /= total;
    }
}

/// **How many allele classes a window of `allele_span` repeats either side holds**: `2 × span + 1`,
/// from `-allele_span` to `+allele_span`.
fn allele_classes_in(allele_span: i32) -> usize {
    (2 * allele_span + 1) as usize
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

/// One tract's read likelihoods, in the form an integral over the tract's length frequencies
/// sweeps over cheaply. The fitting bench's held-out score
/// (`bench_fixtures::fit_bench::score_locus`) is what reads them now; the climb's 256-point
/// average, which they were built for, went with the climb.
///
/// **Rescaled out of log space, once.** The inner loop is *for each point of the integral, for
/// each sample, sum over genotypes*, so doing it in logs would cost one `exp` per genotype per
/// sample **per point**. Subtracting each sample's largest log-likelihood makes the inner sum a
/// plain dot product, and the offsets are added back at the end.
///
/// **Compiled only for the bench and the tests**, with the per-library scoring it relies on:
/// production fits collapse every library to one distribution before scoring.
#[cfg(test)]
struct TractLikelihoods {
    /// One row a sample-with-reads, over the genotypes, **laid end to end in one buffer**: row
    /// `r` is `scaled[r * width..][..width]`.
    ///
    /// **Flat and not a vector of vectors**, because the innermost loop sweeps a row against the
    /// genotype prior in vector lanes: separate heap blocks a sample cost a dependent load before
    /// each row and leave the rows scattered, where one buffer keeps them contiguous and lets the
    /// loop walk the whole tract without leaving cache.
    scaled: Vec<f64>,
    /// Σ over those samples of the log-likelihood each row was divided by.
    ln_offset: f64,
}

/// **The read-length distributions a stratum's reads are scored under**: one for each slippage group
/// and multiplier that some library is scored at, and which of them each library uses.
///
/// **One distribution for all the libraries of a slippage group that share a multiplier**, not one
/// for each library: a sample whose libraries use the same distribution has their counts summed
/// before the logarithm, as the evidence did when it summed them into the slippage group. With
/// every multiplier one, the distributions are the slippage groups' and every sum is the one the fit
/// made before libraries were kept apart, to the bit (`str_library_stutter.md` §3.1).
struct ReadDistributionsByLevel {
    /// Per distribution, per allele offset, the probability of each read bucket. **In
    /// slippage-group order, then by multiplier**, so a sample's rows are summed in the order the
    /// slippage groups were.
    read_distribution_by_allele: Vec<Vec<Vec<f64>>>,
    /// Per library, the position of its distribution in `read_distribution_by_allele`. Read only
    /// by the per-library scoring of [`TractLikelihoods`], so compiled only where that is.
    #[cfg(test)]
    distribution_of_library: Vec<usize>,
}

impl ReadDistributionsByLevel {
    fn of(slippage: &[Slippage], libraries: &[Library], allele_span: i32, read_span: i32) -> Self {
        // Each distinct (slippage group, multiplier), the multiplier by its bits so that a
        // multiplier sorts and compares the same way in both places.
        let key = |library: &Library| {
            (
                library.slippage_group as usize,
                library.stutter_multiplier.to_bits(),
            )
        };
        let mut group_and_multiplier: Vec<(usize, u64)> = libraries.iter().map(key).collect();
        group_and_multiplier.sort_unstable_by(|left, right| {
            left.0
                .cmp(&right.0)
                .then(f64::from_bits(left.1).total_cmp(&f64::from_bits(right.1)))
        });
        group_and_multiplier.dedup();
        #[cfg(test)]
        let distribution_of_library = libraries
            .iter()
            .map(|library| {
                let (group, multiplier) = key(library);
                // PANIC-FREE: every library's key was put in the list just above.
                group_and_multiplier
                    .binary_search_by(|(held_group, held_multiplier)| {
                        held_group.cmp(&group).then(
                            f64::from_bits(*held_multiplier).total_cmp(&f64::from_bits(multiplier)),
                        )
                    })
                    .expect("every library's slippage group and multiplier is listed")
            })
            .collect();
        let read_distribution_by_allele = group_and_multiplier
            .iter()
            .map(|(group, multiplier)| {
                let class = slippage[*group];
                let scored = Slippage {
                    level: library_level(class.level, f64::from_bits(*multiplier)),
                    shorter_share: class.shorter_share,
                    fall_off: class.fall_off,
                };
                (-allele_span..=allele_span)
                    .map(|allele| scored.read_probabilities(allele, read_span))
                    .collect()
            })
            .collect();
        Self {
            read_distribution_by_allele,
            #[cfg(test)]
            distribution_of_library,
        }
    }

    /// The read-length distribution from each allele that `library`'s reads are scored under.
    #[cfg(test)]
    fn of_library(&self, library: LibraryIndex) -> &[Vec<f64>] {
        &self.read_distribution_by_allele[self.distribution_of_library[library.get()]]
    }

    /// **One sample's reads, one row for each distribution its libraries use**, in distribution
    /// order, for a sample with several libraries. Libraries using one distribution have their
    /// counts summed into `summed_by_distribution`, a scratch buffer whose count vectors are reused
    /// from one sample to the next; the rows returned are the part of it in use, so they last until
    /// it is next used.
    #[cfg(test)]
    fn rows_by_distribution<'a>(
        &self,
        sample: &SampleTractReads,
        summed_by_distribution: &'a mut Vec<(usize, Vec<u32>)>,
    ) -> &'a [(usize, Vec<u32>)] {
        let mut used = 0;
        for (library, counts) in &sample.by_library {
            let distribution = self.distribution_of_library[library.get()];
            if let Some((_, sum)) = summed_by_distribution[..used]
                .iter_mut()
                .find(|(held, _)| *held == distribution)
            {
                for (into, reads) in sum.iter_mut().zip(counts) {
                    *into += reads;
                }
            } else if let Some((held, sum)) = summed_by_distribution.get_mut(used) {
                *held = distribution;
                sum.clone_from(counts);
                used += 1;
            } else {
                summed_by_distribution.push((distribution, counts.clone()));
                used += 1;
            }
        }
        let rows = &mut summed_by_distribution[..used];
        rows.sort_unstable_by_key(|(distribution, _)| *distribution);
        rows
    }
}

/// **Adds to `total` the log-likelihood of one library's `counts` under the genotype of alleles
/// `first` and `second`**, each read drawn from either allele with even odds; `per_allele` is the
/// library's read-length distribution from each allele.
#[cfg(test)]
#[inline(always)]
fn add_ln_reads(
    total: &mut f64,
    per_allele: &[Vec<f64>],
    counts: &[u32],
    first: usize,
    second: usize,
) {
    for (bucket, reads) in counts.iter().enumerate() {
        if *reads == 0 {
            continue;
        }
        let probability = 0.5 * (per_allele[first][bucket] + per_allele[second][bucket]);
        *total += f64::from(*reads) * float::ln(probability.max(1e-300));
    }
}

#[cfg(test)]
impl TractLikelihoods {
    /// `distributions` holds the read-length distributions the tract's libraries are scored under.
    fn of(
        tract: &TractReads,
        distributions: &ReadDistributionsByLevel,
        genotypes: &[(usize, usize)],
    ) -> Self {
        let width = genotypes.len();
        let mut scaled: Vec<f64> = Vec::with_capacity(tract.samples.len() * width);
        let mut ln_offset = 0.0;
        let mut summed_by_distribution: Vec<(usize, Vec<u32>)> = Vec::new();
        for sample in &tract.samples {
            let start = scaled.len();
            // **A sample with one library — every sample of a cohort of one read group a sample —
            // runs the loop the fit ran before libraries were kept apart**, with its distribution
            // fetched once: going through the general path below cost that case 1.75 times the
            // time (review of step H1, 400 tracts of 200 samples at three reads).
            if let [(library, counts)] = sample.by_library.as_slice() {
                let per_allele = distributions.of_library(*library);
                scaled.extend(genotypes.iter().map(|(first, second)| {
                    let mut total = 0.0;
                    add_ln_reads(&mut total, per_allele, counts, *first, *second);
                    total
                }));
            } else {
                let rows = distributions.rows_by_distribution(sample, &mut summed_by_distribution);
                scaled.extend(genotypes.iter().map(|(first, second)| {
                    let mut total = 0.0;
                    for (distribution, counts) in rows {
                        let per_allele = &distributions.read_distribution_by_allele[*distribution];
                        add_ln_reads(&mut total, per_allele, counts, *first, *second);
                    }
                    total
                }));
            }
            let row = &mut scaled[start..];
            let largest = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            ln_offset += largest;
            for value in row {
                *value = float::exp(*value - largest);
            }
        }
        Self { scaled, ln_offset }
    }
}

// ---------------------------------------------------------------------
// Small numerics
// ---------------------------------------------------------------------

// Only the cohort draw and the tests use this.

#[cfg(test)]
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
/// **Every read group is kept apart as a library of its own** ([`StratumEvidence::libraries`],
/// one list shared by every stratum), each at a stutter multiplier of one. A read group
/// `slippage_group_of` does not name is put in slippage group 0, and the log counts such read
/// groups unless the map is empty.
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
    // **Every read group of the cohort is a library of its own**, in the cohort's read-group order;
    // one `slippage_group_of` does not name is in group 0, as its reads always were.
    let libraries: Arc<[Library]> = cohort
        .read_groups()
        .iter()
        .map(|read_group| {
            Library::in_group(
                *read_group,
                *slippage_group_of.get(read_group).unwrap_or(&0),
            )
        })
        .collect();
    let library_of: BTreeMap<ReadGroupId, LibraryIndex> = libraries
        .iter()
        .enumerate()
        .map(|(index, library)| (library.read_group, LibraryIndex(index as u32)))
        .collect();

    // How many tracts each stratum holds, which is the length its sections are built at.
    let mut tracts_in: BTreeMap<Stratum, usize> = BTreeMap::new();
    for stratum in strata {
        *tracts_in.entry(*stratum).or_insert(0) += 1;
    }
    let band: Vec<Stratum> = tracts_in.keys().copied().collect();
    let names: Vec<String> = cohort.sample_names().map(str::to_string).collect();

    // **The whole band at once, and that is this step's shape rather than the design's.** The
    // repeat-tract fit draws each motif length's numbers as curves across its classes
    // (`penalised_fit::fit_repeat_tracts`), so it is handed every stratum together; reading the
    // census one motif length at a time is the census-read step's (`str_penalised_fit.md` §6).
    let stage = StageProgress::begin(format!(
        "repeat-tract evidence: gathering {} tract(s) in {} strata from {} sample(s)",
        strata.len(),
        tracts_in.len(),
        names.len()
    ));
    let undeclared = libraries
        .iter()
        .filter(|library| !slippage_group_of.contains_key(&library.read_group))
        .count();
    if undeclared > 0 && !slippage_group_of.is_empty() {
        stage.always(|into| {
            format!(
                "repeat-tract evidence: {undeclared} read group(s) of {} have no slippage group \
                 declared and are scored in group 0; {into}",
                libraries.len()
            )
        });
    }
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
                    libraries: Arc::clone(&libraries),
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
                    // **Summed over every sample and read group, then tested once**: the guard
                    // is a test of the tract (spec §3.3). Tested a section at a time, one stray
                    // read in one sample of thousands removed the tract for the whole cohort.
                    let mut guard = GuardTally::default();
                    for (sample_index, sample) in by_sample.iter().enumerate() {
                        let mut by_library: BTreeMap<LibraryIndex, Vec<u32>> = BTreeMap::new();
                        for (read_group, records) in sample.get(stratum).into_iter().flatten() {
                            guard.add(records.section_guard_tally(locus));
                            if records.state(locus) != SsrLocusState::Crossed {
                                continue;
                            }
                            let library = *library_of
                                .get(read_group)
                                .expect("every section's read group is one of the cohort's");
                            let counts = by_library
                                .entry(library)
                                .or_insert_with(|| vec![0; buckets]);
                            let offsets = records.offsets(locus);
                            for (bucket, count) in counts.iter_mut().enumerate() {
                                *count +=
                                    u32::from(offsets.at(bucket as i32 - RECORDED_OFFSET_RANGE));
                            }
                        }
                        if !by_library.is_empty() {
                            reads.samples.push(SampleTractReads {
                                sample: sample_index as u32,
                                by_library: by_library.into_iter().collect(),
                            });
                        }
                    }
                    evidence.guard_reads += guard.non_whole_reads;
                    if guard.is_over_threshold() {
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
// Drawn strata: one generator, shared by the tests and the bench
// ---------------------------------------------------------------------

/// A stratum drawn at a known truth, for anything that has to fit evidence it already knows
/// the answer to.
///
/// **Compiled under `cfg(test)` and under the `bench-fixtures` feature, and nowhere else.** The
/// fit's tests and the research bench (`examples/ng_str_fit_bench.rs`) need drawn strata, and they
/// need the *same* ones: a second generator would be a second thing to keep agreeing, and a bench
/// drawn differently from the tests would be measuring a workload no test has ever checked.
///
/// Nothing here is production code: a release build without the feature compiles none of it.
#[cfg(test)]
pub mod bench_fixtures {
    use super::*;

    pub mod cohort_draw;

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

        /// A Poisson number of mean `mean`: by multiplying uniforms below a mean of 10, by Hörmann's
        /// transformed rejection (PTRS, the method numpy uses) above it.
        fn poisson(&mut self, mean: f64) -> u32 {
            if mean <= 0.0 {
                return 0;
            }
            if mean < 10.0 {
                let limit = float::exp(-mean);
                let mut product = self.uniform();
                let mut count = 0;
                while product > limit {
                    product *= self.uniform();
                    count += 1;
                }
                return count;
            }
            let root = mean.sqrt();
            let ln_mean = float::ln(mean);
            let b = 0.931 + 2.53 * root;
            let a = -0.059 + 0.02483 * b;
            let inverse_alpha = 1.1239 + 1.1328 / (b - 3.4);
            let v_r = 0.9277 - 3.6224 / (b - 2.0);
            loop {
                let u = self.uniform() - 0.5;
                let v = self.uniform();
                let us = 0.5 - u.abs();
                let k = ((2.0 * a / us + b) * u + mean + 0.43).floor();
                if us >= 0.07 && v <= v_r {
                    return k as u32;
                }
                if k < 0.0 || (us < 0.013 && v > us) {
                    continue;
                }
                if float::ln(v) + float::ln(inverse_alpha) - float::ln(a / (us * us) + b)
                    <= -mean + k * ln_mean - ln_gamma(k + 1.0)
                {
                    return k as u32;
                }
            }
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
    /// (`2 × span + 1`), so a caller fitting at a window of `span` draws at the same
    /// span and hands the fit a spectrum of that length.
    ///
    /// **Every sample gets `depth` reads at every tract**, where a real cohort at three reads a
    /// position puts a read at a tract in a minority of its samples. So a drawn stratum of `n`
    /// samples costs more a tract than a recorded one of `n` samples — which is the right way
    /// round for a benchmark, and the wrong way round for extrapolating a wall time to a cohort.
    ///
    /// **Each sample is a library of its own**, library `s` read group `s`, every one in slippage
    /// group 0 at a multiplier of one ([`one_library_a_sample`]).
    ///
    /// **One class of the count-level simulator ([`cohort_draw`]) at a fixed depth**: the locus draw
    /// is the simulator's, taking its random numbers in the order this function always took them,
    /// so a seed draws what it drew before the simulator existed (tested).
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
        let per_allele: Vec<Vec<f64>> = (0..classes)
            .map(|class| slippage.read_probabilities(class as i32 - span, span))
            .collect();
        let per_library = [per_allele];
        let excess = vec![homozygote_excess; samples];
        // **One library for the reads' distribution, `samples` in the evidence**: each sample is a
        // library of its own below, every one scored at the same stutter.
        let drawing = cohort_draw::LocusDrawing {
            per_library: &per_library,
            libraries: 1,
            homozygote_excess: &excess,
            depth: cohort_draw::Depth::Fixed(depth),
            buckets: (2 * span + 1) as usize,
        };
        let mut rng = Rng(seed);
        let mut drawn = Vec::with_capacity(tracts);
        for _ in 0..tracts {
            let alpha: Vec<f64> = spectrum.iter().map(|q| concentration * q).collect();
            let frequencies = cohort_draw::draw_locus_frequencies(&mut rng, &alpha);
            let mut reads = cohort_draw::draw_locus(&mut rng, &frequencies, &drawing);
            for sample in &mut reads.samples {
                sample.by_library[0].0 = LibraryIndex(sample.sample);
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
            libraries: one_library_a_sample(samples),
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

    /// **`samples` libraries, one a sample**: library `s` is read group `s`, in slippage group 0 at a
    /// multiplier of one.
    pub fn one_library_a_sample(samples: usize) -> Arc<[Library]> {
        (0..samples as u32)
            .map(|sample| Library::in_group(ReadGroupId(sample), 0))
            .collect()
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

    /// **An error withheld for too few loci is no error**: it carries no value, and a stratum
    /// whose every number was withheld has none to report.
    #[test]
    fn an_error_withheld_for_too_few_loci_has_no_value() {
        assert_eq!(StratumError::TooFewLoci.value(), None);
        let withheld = SlippageErrors {
            level: StratumError::TooFewLoci,
            shorter_share: StratumError::TooFewLoci,
            fall_off: StratumError::TooFewLoci,
        };
        let errors = StratumErrors {
            slippage: vec![Some(withheld)],
            allele_length_variation: vec![StratumError::TooFewLoci; 3],
            concentration: StratumError::TooFewLoci,
        };
        assert!(!errors.any());
    }

    // -----------------------------------------------------------------
    // The guard: one test over the whole locus
    // -----------------------------------------------------------------

    /// A cohort of `samples` samples, one read group each, holding one dinucleotide tract.
    /// Every sample has two reads at the reference length and one a repeat short; the sample
    /// `stray` also has one read that differs by a non-whole number of repeats.
    fn cohort_with_one_stray_read(samples: u32, stray: u32) -> CohortCensusEvidence {
        use crate::parameter_estimation::joint::census::{
            GuardObservation, NamedReadGroup, OffsetCounts, SampleCensusEvidence, Section,
            SectionKey, WalkedBits,
        };
        use crate::parameter_estimation::joint::census_file::tests_support::terms;

        let stratum = Stratum {
            period: 2,
            reference_repeats: 6,
        };
        let censuses = (0..samples)
            .map(|sample| {
                let group = ReadGroupId(sample);
                let mut offsets = OffsetCounts::default();
                offsets.add(0, 2);
                offsets.add(-1, 1);
                let mut walked = WalkedBits::none_of(1);
                walked.set(0);
                let guard = if sample == stray {
                    vec![GuardObservation {
                        locus: 0,
                        length_difference: 1,
                        reads: 1,
                    }]
                } else {
                    Vec::new()
                };
                let records =
                    SsrEvidence::from_parts(vec![offsets], 0, walked, 0, guard, Vec::new());
                let name = format!("s{sample:04}");
                SampleCensusEvidence::resident(
                    name.clone(),
                    terms(),
                    NamedReadGroup::drawn_for(&name, [group]),
                    BTreeMap::new(),
                    BTreeMap::from([(SectionKey::Ssr(group, stratum), Section::Ssr(records))]),
                )
            })
            .collect();
        CohortCensusEvidence::new(censuses).expect("every sample recorded one way")
    }

    /// **The guard is a test of the locus, not of one sample's reads at it** (spec
    /// `parameter_prepass_joint_records.md` §3.3). Here one sample in forty has a single read
    /// that is not a whole number of repeats off. Pooled over the cohort that is 1 such read
    /// among 41 that differ from the reference length — under one in ten — so the tract is
    /// kept. Tested sample by sample, that one sample's reads are 1 in 2 and the whole
    /// cohort would lose the tract.
    #[test]
    fn one_stray_read_in_one_sample_does_not_remove_the_tract_for_the_cohort() {
        let stratum = Stratum {
            period: 2,
            reference_repeats: 6,
        };
        let mut cohort = cohort_with_one_stray_read(40, 17);
        let groups: BTreeMap<ReadGroupId, u32> = (0..40).map(|g| (ReadGroupId(g), 0)).collect();
        let gathered = gather_strata(&mut cohort, &[stratum], &groups).expect("resident");
        assert_eq!(gathered[0].tracts_over_guard_threshold, 0);
        assert_eq!(gathered[0].tracts.len(), 1, "the tract is kept");
        assert_eq!(gathered[0].guard_reads, 1);
    }

    /// And the pooled test still removes a tract whose cohort-wide share is over one in ten: a
    /// cohort of the one sample with the stray read is 1 non-whole read in 2 that differ.
    #[test]
    fn a_tract_over_one_in_ten_across_the_cohort_is_still_removed() {
        let stratum = Stratum {
            period: 2,
            reference_repeats: 6,
        };
        let mut cohort = cohort_with_one_stray_read(1, 0);
        let groups = BTreeMap::from([(ReadGroupId(0), 0)]);
        let gathered = gather_strata(&mut cohort, &[stratum], &groups).expect("resident");
        assert_eq!(gathered[0].tracts_over_guard_threshold, 1);
        assert!(gathered[0].tracts.is_empty());
    }

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

    /// Two slippage groups' numbers, at different levels, shares and fall-offs.
    fn two_groups_slippage() -> [Slippage; 2] {
        [
            Slippage {
                level: 0.1,
                shorter_share: 0.8,
                fall_off: 0.25,
            },
            Slippage {
                level: 0.23,
                shorter_share: 0.6,
                fall_off: 0.3,
            },
        ]
    }

    /// One tract read by one sample through two libraries, `first`'s reads through library 0 and
    /// `second`'s through library 1.
    fn one_sample_through_two_libraries(first: &[u32], second: &[u32]) -> TractReads {
        TractReads {
            samples: vec![SampleTractReads {
                sample: 0,
                by_library: vec![
                    (LibraryIndex(0), first.to_vec()),
                    (LibraryIndex(1), second.to_vec()),
                ],
            }],
        }
    }

    /// A tract's likelihood table and its offset, as bits.
    fn the_bits_of(likelihoods: &TractLikelihoods) -> Vec<u64> {
        let mut bits: Vec<u64> = likelihoods.scaled.iter().map(|x| x.to_bits()).collect();
        bits.push(likelihoods.ln_offset.to_bits());
        bits
    }

    /// **Libraries at different levels are scored under different distributions, and libraries at
    /// one level under one**: two libraries of one group at one multiplier share a distribution, a
    /// library at another multiplier has its own, and libraries of two groups never share — and the
    /// distributions run in slippage-group order, then by multiplier.
    #[test]
    fn a_read_distribution_is_one_slippage_group_at_one_multiplier() {
        let slippage = two_groups_slippage();
        let library = |group: u32, multiplier: f64| Library {
            read_group: ReadGroupId(0),
            slippage_group: group,
            stutter_multiplier: multiplier,
        };
        let distributions = ReadDistributionsByLevel::of(
            &slippage,
            &[
                library(1, 1.0),
                library(0, 2.0),
                library(0, 1.0),
                library(0, 1.0),
            ],
            1,
            1,
        );
        assert_eq!(distributions.distribution_of_library, vec![2, 1, 0, 0]);
        let doubled = Slippage {
            level: 0.2,
            ..slippage[0]
        };
        assert_eq!(
            distributions.read_distribution_by_allele[1][1],
            doubled.read_probabilities(0, 1)
        );
        assert_eq!(
            distributions.read_distribution_by_allele[2][1],
            slippage[1].read_probabilities(0, 1)
        );
    }

    /// **A sample's rows are summed in distribution order, not in the order its libraries are
    /// listed**: a sample whose library 0 is in slippage group 1 and library 1 in group 0 scores
    /// to the bit as the same reads listed the other way round.
    #[test]
    fn tract_likelihoods_do_not_depend_on_the_order_a_samples_libraries_are_listed_in() {
        let slippage = two_groups_slippage();
        let genotypes = genotype_pairs(3);
        let (one, other) = ([3, 7, 2], [5, 1, 9]);
        let inverted = [
            Library::in_group(ReadGroupId(0), 1),
            Library::in_group(ReadGroupId(1), 0),
        ];
        let in_order = [
            Library::in_group(ReadGroupId(0), 0),
            Library::in_group(ReadGroupId(1), 1),
        ];
        let listed_against_the_groups = TractLikelihoods::of(
            &one_sample_through_two_libraries(&one, &other),
            &ReadDistributionsByLevel::of(&slippage, &inverted, 1, 1),
            &genotypes,
        );
        let listed_with_the_groups = TractLikelihoods::of(
            &one_sample_through_two_libraries(&other, &one),
            &ReadDistributionsByLevel::of(&slippage, &in_order, 1, 1),
            &genotypes,
        );
        assert_eq!(
            the_bits_of(&listed_against_the_groups),
            the_bits_of(&listed_with_the_groups)
        );
    }

    /// **Two libraries of one sample and one slippage group at different multipliers are scored
    /// apart**: the sample scores to the bit as if the second library were in a group of its own
    /// at the doubled level, and not as the two libraries' reads summed at the class rate.
    #[test]
    fn a_samples_libraries_at_two_multipliers_are_two_rows() {
        let class = two_groups_slippage()[0];
        let genotypes = genotype_pairs(3);
        let (one, other) = ([3, 7, 2], [5, 1, 9]);
        let at_two_multipliers = [
            Library::in_group(ReadGroupId(0), 0),
            Library {
                stutter_multiplier: 2.0,
                ..Library::in_group(ReadGroupId(1), 0)
            },
        ];
        let scored = TractLikelihoods::of(
            &one_sample_through_two_libraries(&one, &other),
            &ReadDistributionsByLevel::of(&[class], &at_two_multipliers, 1, 1),
            &genotypes,
        );
        let doubled = Slippage {
            level: class.level * 2.0,
            ..class
        };
        let in_two_groups = TractLikelihoods::of(
            &one_sample_through_two_libraries(&one, &other),
            &ReadDistributionsByLevel::of(
                &[class, doubled],
                &[
                    Library::in_group(ReadGroupId(0), 0),
                    Library::in_group(ReadGroupId(1), 1),
                ],
                1,
                1,
            ),
            &genotypes,
        );
        assert_eq!(the_bits_of(&scored), the_bits_of(&in_two_groups));
        let summed: Vec<u32> = one.iter().zip(&other).map(|(a, b)| a + b).collect();
        let at_the_class_rate = TractLikelihoods::of(
            &TractReads {
                samples: vec![SampleTractReads {
                    sample: 0,
                    by_library: vec![(LibraryIndex(0), summed)],
                }],
            },
            &ReadDistributionsByLevel::of(&[class], &[Library::in_group(ReadGroupId(0), 0)], 1, 1),
            &genotypes,
        );
        assert_ne!(the_bits_of(&scored), the_bits_of(&at_the_class_rate));
    }

    /// **Two libraries of one slippage group at a multiplier of one score to the bit as one
    /// library holding their summed reads**, which is the claim on [`ReadDistributionsByLevel`]:
    /// libraries sharing a distribution have their counts summed before the logarithm.
    #[test]
    fn tract_likelihoods_of_two_libraries_at_multiplier_one_equal_their_summed_reads() {
        let class = two_groups_slippage()[0];
        let genotypes = genotype_pairs(3);
        let (one, other) = ([3, 7, 2], [5, 1, 9]);
        let kept_apart = TractLikelihoods::of(
            &one_sample_through_two_libraries(&one, &other),
            &ReadDistributionsByLevel::of(
                &[class],
                &[
                    Library::in_group(ReadGroupId(0), 0),
                    Library::in_group(ReadGroupId(1), 0),
                ],
                1,
                1,
            ),
            &genotypes,
        );
        let summed = TractLikelihoods::of(
            &TractReads {
                samples: vec![SampleTractReads {
                    sample: 0,
                    by_library: vec![(LibraryIndex(0), vec![8, 8, 11])],
                }],
            },
            &ReadDistributionsByLevel::of(&[class], &[Library::in_group(ReadGroupId(0), 0)], 1, 1),
            &genotypes,
        );
        assert_eq!(the_bits_of(&kept_apart), the_bits_of(&summed));
    }

    /// **A library in a slippage group the evidence does not declare is refused by name**, by the
    /// check the fit runs where it takes the stratum.
    #[test]
    #[should_panic(expected = "is in slippage group 1 of the 1 the evidence declares")]
    fn a_library_past_the_declared_groups_is_refused() {
        let mut evidence = draw_stratum(
            two_groups_slippage()[0],
            &spectrum_of(3),
            0.5,
            0.4,
            8,
            2,
            3,
            1,
            5,
        );
        let mut libraries = evidence.libraries.to_vec();
        libraries[1].slippage_group = 1;
        evidence.libraries = libraries.into();
        evidence.check_libraries();
    }

    /// **A sample's reads naming a library the evidence does not list are refused by name.**
    #[test]
    #[should_panic(expected = "out of order or past the 2 the evidence lists")]
    fn reads_naming_a_library_past_the_list_are_refused() {
        let mut evidence = draw_stratum(
            two_groups_slippage()[0],
            &spectrum_of(3),
            0.5,
            0.4,
            8,
            2,
            3,
            1,
            5,
        );
        evidence.tracts[0].samples[0].by_library[0].0 = LibraryIndex(2);
        evidence.check_libraries();
    }

    /// **A multiplier lifts a level no higher than the ceiling, and never lowers a class's own**:
    /// one leaves any level as it is, including one above the ceiling and one that is not a number.
    #[test]
    fn a_library_level_is_the_class_rate_times_the_multiplier_below_the_ceiling() {
        assert_eq!(library_level(0.1, 1.0), 0.1);
        assert_eq!(library_level(0.1, 0.5), 0.05);
        assert_eq!(library_level(0.1, 3.0), 0.1 * 3.0);
        assert_eq!(library_level(0.5, 3.0), MULTIPLIED_LEVEL_CEILING);
        assert_eq!(library_level(0.97, 1.0), 0.97);
        assert_eq!(library_level(0.97, 2.0), 0.97);
        assert_eq!(library_level(0.97, 0.5), 0.485);
        assert!(library_level(f64::NAN, 1.0).is_nan());
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
}

#[cfg(test)]
mod evidence_size {
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
                            by_library: vec![(LibraryIndex(0), vec![1])],
                        })
                        .collect(),
                })
                .collect(),
            read_span: RECORDED_OFFSET_RANGE,
            groups: 1,
            libraries: Arc::from([Library::in_group(ReadGroupId(0), 0)]),
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
            + std::mem::size_of::<(LibraryIndex, Vec<u32>)>()
            + std::mem::size_of::<u32>();
        assert_eq!(
            stratum.heap_bytes(),
            3 * std::mem::size_of::<TractReads>() + 5 * row
        );
    }
}
