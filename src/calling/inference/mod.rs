//! Step 9 — turning one locus's evidence and the run's frozen parameters into genotypes.
//!
//! **This file is the seam and its configuration; the loops themselves arrive beside it.**
//! [`LocusGenotyper`] is the one boundary everything the calling arms share crosses — the
//! evidence, the parameters, the candidate alleles and the worker's scratch — and
//! [`CallingLoopConfig`] is every switch of the design, held as a **value** rather than as a
//! code path.
//!
//! ## Two loops, one inside the other, and the outer one ships switched off
//!
//! The inner loop is the one that runs: guess the allele frequencies, work out each
//! sample's genotype probabilities under that guess, add them up to get better frequencies,
//! repeat until they stop moving. Around it sits **the slippage round**, which re-fits how
//! often a read gains or loses a whole repeat from the locus's own reads and runs the
//! frequency loop again on the new numbers ([`SlippageRefitConfig`], `max_rounds` 0, so a
//! single pass through its body).
//!
//! **Allele discovery is not a loop here.** It admits the tract lengths that recur too often
//! in one sample to be slippage ([`DiscoveryConfig`]), and it runs as a **pre-pass inside
//! candidate selection** ([`DiscoveryMode::BeforeTheLoop`], wired in
//! [`select_ssr`](crate::calling::allele_candidates::ssr::select_ssr)). Building it
//! showed the eligible set is a function of the observations and the
//! candidate table alone — no posterior enters it — so a discovery round wrapped around the
//! frequency loop would admit nothing the pre-pass had not. The two round-wrapped settings
//! the spec described were therefore examined and removed
//! (`doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md` §6.5).
//!
//! The slippage round's three settings — *frozen*, *pulled back part way* and *free* — are
//! three values of one code path rather than three implementations to keep in step
//! (`doc/devel/ng/spec/calling_em_loop.md` §5.1). Its body is built: the arithmetic in
//! [`slippage_refit`], the round itself in the loop's driver. [`LocusGenotyper::call_locus`]
//! takes a [`RunnableCallingLoopConfig`], which has one constructor and it is
//! [`CallingLoopConfig::validate`], so a value out of range cannot reach the loop at all.
//!
//! ## What is deliberately not configured here
//!
//! **The allele cap.** How many alleles a locus may be called over belongs to candidate
//! selection, which owns it together with the support bar
//! ([`CandidateSelectionConfig`](crate::calling::allele_candidates::CandidateSelectionConfig)),
//! and states it as a
//! [`MaxCandidateAlleles`](crate::calling::allele_candidates::MaxCandidateAlleles) — a
//! type that refuses anything below two, because a cap admitting no alternative is refusal
//! under another name. A second field of the same name here would be two spellings of one
//! rule, and the weaker of the two (`doc/devel/ng/arch/calling_em_loop.md` §2.1).
//!
//! **The same argument decides how discovery's evidence bar is held.** A read floor plus a
//! share of one sample's own reads, both of which must clear, is exactly the rule the merge
//! already ships as [`MinAltReads`] — so discovery reuses that type rather than restating it
//! in a `u32` and an `f64`. The *numbers* stay discovery's own, because spec §4.1's third
//! open question sweeps them independently of selection's; what is shared is the type, and
//! with it the impossibility of a negative share or a floor of zero reads.

pub mod repeat_tract_parameters;
pub(crate) mod slippage_refit;
pub mod summarise_condition;

use std::num::NonZeroU32;

use crate::calling::{
    CallingScratch, CandidateAlleles, FrozenParameters, LocusEvidence, LocusInference,
};
use crate::cohort_observation::{MinAltObs, MinAltReadShare, MinAltReads};

/// Stop the frequency loop when the largest change in the cohort's expected allele copies,
/// **divided by the number of chromosomes in the cohort**, falls below this.
///
/// **The division is what makes the number mean the same thing across the range this caller
/// commits to.** Expected copies are a count and this threshold is a fraction: `1e-3` in raw
/// copies says something different at one sample from what it says at a thousand, so a
/// criterion written on counts tightens by the cohort size
/// (`doc/devel/ng/spec/calling_em_loop.md` §6).
///
/// **Inherited from production and soft.** Production relaxed it from `1e-4` after tomato
/// records plateaued in the `5e-4 … 1e-3` band and hit the pass cap without crossing it
/// (`d9e7b076:src/var_calling/posterior_engine.rs`). Nothing has measured it on this caller's
/// range — a thousand samples, or three reads a position — and what would set it is the
/// distribution of [`LocusInference::passes`] across a real panel (spec §12, question 4).
pub const DEFAULT_CONVERGENCE_THRESHOLD: f64 = 1e-3;

/// The largest convergence threshold that is a stopping rule rather than an immediate exit.
///
/// **A threshold of 1 stops every locus after one pass and calls it settled.** The quantity
/// compared against it is a per-allele change already divided by the cohort's chromosomes, so
/// it lies within `[0, 1]` and cannot exceed 1 — which makes any threshold at or above 1
/// satisfied unconditionally, by a loop that has done nothing.
///
/// **Inherited with the value it bounds.** Production validates the same field within
/// `(0, 0.1]` for this reason, and its own constant says so in as many words; spec §6 gives
/// the division's purpose as keeping the threshold *"and its validation range"* meaningful
/// across the scale. `0.1` is loose but not degenerate — a hundred times the default, which
/// is headroom for a run that wants to stop early without reaching the degenerate case.
pub const CONVERGENCE_THRESHOLD_RANGE_MAX: f64 = 0.1;

/// Give up on the frequency loop after this many passes, emit the locus with
/// [`LocusInference::converged`] false, and carry on.
///
/// **A locus that will not settle is emitted and flagged, never dropped and never fatal** —
/// production retired its non-convergence error for exactly this reason, so that one hard
/// site does not kill a whole cohort run (spec §6).
///
/// **Inherited and soft, and it is a ceiling rather than a tuned value:** production's own
/// comment records the observed need as 3 to 5 passes on the GATK reference data, so 50 is
/// ten times what was seen.
pub const DEFAULT_MAX_PASSES: u32 = 50;

/// Ship with the locus's slippage numbers **frozen** at what the parameter fit measured for
/// its stratum: no re-fit rounds, and so no rebuild of the likelihood table.
///
/// **Frozen is this code path at zero rounds**, not a second implementation. Whether a locus
/// should be allowed to pull its own numbers away from its stratum's is unmeasured — the
/// reads a tract's numbers would be re-fitted from are the very reads being genotyped, so
/// they can end up describing this locus's noise rather than its chemistry (spec §5.1, and
/// §12's question 2 is the measurement).
pub const DEFAULT_SLIPPAGE_REFIT_ROUNDS: u32 = 0;

/// How hard a re-fit is pulled back toward what the parameter fit measured for the locus's
/// stratum, for the **two numbers that say what a slip looks like when one happens**: which
/// way it goes, and how fast longer slips get rarer. In pseudo-counts.
///
/// The crate's names for the pair are `Slippage::shorter_share` — the share of slips that
/// shorten the tract rather than lengthen it — and `Slippage::fall_off`, how quickly a slip
/// of two repeats becomes rarer than a slip of one.
///
/// Production's value. A tract with few slipped reads cannot overcome it and collapses back
/// to its stratum's numbers, which is the behaviour the pull-back exists for; zero is
/// HipSTR's setting, where the locus's own reads set the numbers outright. **Inherited and
/// never measured here.**
pub const DEFAULT_DIRECTION_AND_FALL_OFF_PULL_BACK_PSEUDOCOUNTS: f64 = 50.0;

/// The same pull-back for **how often a read slips at all** — the crate's `Slippage::level` —
/// counted in slipped reads rather than in pseudo-counts. Production's value, inherited and
/// never measured here.
///
/// **What it is pulled back *toward* moved after the spec was written.** The slippage level
/// is no longer a per-cell number but a curve in repeat count, fitted once per motif period
/// and read off at the cell (`doc/devel/ng/spec/str_slippage_level_curve.md`), so the target
/// is a point on a fitted line rather than a cell's own estimate.
pub const DEFAULT_LEVEL_PULL_BACK_SLIPPED_READS: f64 = 20.0;

/// Stop the slippage rounds when every re-fitted number moves less than this between rounds.
///
/// **The outer rounds need their own stopping rule and this is not the frequency loop's.**
/// They stop on their own numbers' movement or at the round cap — production's rule rather
/// than HipSTR's likelihood test, because a likelihood test would need the table rebuilt in
/// order to be read (spec §6).
///
/// **The value is production's, which uses it twice**: one threshold on the shape
/// coefficients and one on the level multiplier (`EmCfg::dev_default`,
/// `d9e7b076:src/ssr/cohort/em.rs`). ng holds the two at one number until something measures them
/// apart.
pub const DEFAULT_ROUND_CONVERGENCE_THRESHOLD: f64 = 1e-3;

/// How many of one sample's reads must show a tract length before discovery may admit it
/// as an allele.
///
/// **Both halves of the shipped bar must clear, so that a single stray read cannot mint an
/// allele.** Inherited from HipSTR's high-depth human setting and **soft** — the two halves
/// bind at opposite ends of the depth range, and below about 13 reads only this one does
/// anything, because 2 reads already clears 15% (spec §4.1). That is a property of the
/// shipped numbers rather than a guarantee of the type: [`MinAltReads`] admits a share of
/// zero, and spec §4.1's third open question sweeps both halves.
pub const DEFAULT_DISCOVERY_MIN_READS: u32 = 2;

/// What share of one sample's **tract-spanning** reads must show a length before discovery
/// may admit it. HipSTR's, inherited and soft — see [`DEFAULT_DISCOVERY_MIN_READS`] for
/// why the pair cannot be right at both ends of the depth range.
///
/// **Of one sample's reads, not of the cohort's.** HipSTR's other admission route is
/// cohort-wide — above 5% of samples, or 5% of reads — and ng builds neither, so the two
/// readings are live and this is the one that applies.
pub const DEFAULT_DISCOVERY_MIN_SPANNING_READ_SHARE: f64 = 0.15;

/// How much evidence one sample must show before discovery admits a tract length:
/// [`DEFAULT_DISCOVERY_MIN_READS`] reads, or
/// [`DEFAULT_DISCOVERY_MIN_SPANNING_READ_SHARE`] of that sample's tract-spanning reads,
/// whichever is more.
///
/// **The merge's own type, with discovery's own numbers.** Sharing the type is what stops a
/// negative share or a floor of zero reads being expressible at all; keeping the numbers
/// separate is what lets spec §4.1's third open question sweep them without moving the rule
/// that decides whether a locus is built.
pub const DEFAULT_DISCOVERY_BAR: MinAltReads = MinAltReads {
    floor: MinAltObs(non_zero_default(DEFAULT_DISCOVERY_MIN_READS)),
    share: MinAltReadShare::new_or_panic(DEFAULT_DISCOVERY_MIN_SPANNING_READ_SHARE),
};

/// Turn a default written as a plain `u32` into the [`NonZeroU32`] the field holds.
///
/// **Defaults are written as `u32` because that is where an operator reading the source looks
/// for them**, and a `NonZeroU32` literal is not readable at its declaration. The conversion
/// panics at compile time on a zero, so a default that broke its own field's invariant would
/// fail the build rather than at run time. The merge keeps a private twin of this for the
/// same reason.
const fn non_zero_default(default_value: u32) -> NonZeroU32 {
    match NonZeroU32::new(default_value) {
        Some(value) => value,
        None => panic!("a calling-loop default that fills a non-zero field must be non-zero"),
    }
}

/// Which of the slippage re-fit's two pull-backs a refusal is about.
///
/// **A type rather than a string**, so that a caller can act on the answer and the compiler
/// checks the match. The two are pulled back in different units and toward different things,
/// which is why they are two settings rather than one.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum PullBack {
    /// What a slip looks like when one happens — which way it goes, and how fast longer
    /// slips get rarer. Weighed in pseudo-counts.
    DirectionAndFallOff,
    /// How often a read slips at all. Weighed in slipped reads.
    Level,
}

impl std::fmt::Display for PullBack {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::DirectionAndFallOff => "direction-and-fall-off",
            Self::Level => "level",
        })
    }
}

/// Re-fitting the locus's own slippage numbers: **how many rounds, and how hard the locus is
/// allowed to pull them away from its stratum's.**
///
/// How often a read gains or loses a whole repeat is measured before calling starts, pooled
/// over every tract that shares a motif length and a repeat count. The case for letting one
/// tract depart from its class is that a tract can behave unlike it — an interruption, a
/// nearby indel, somatic instability. The case against is that the reads it would be
/// re-fitted from are the reads being genotyped. **Nobody has measured which effect is
/// larger**, so ng ships frozen, builds the machinery for the other two settings, and hands
/// the choice to a measurement (spec §5.1, §12's question 2).
///
/// The three settings are three values of this one type: **frozen** at zero rounds,
/// **pulled back part way** at production's pseudo-counts, and **free** at zero pull-back.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SlippageRefitConfig {
    /// Zero is frozen — the default, and the whole of the frozen setting: no rounds, and so
    /// no rebuild of the likelihood table. Production caps at three.
    ///
    /// Defaults to [`DEFAULT_SLIPPAGE_REFIT_ROUNDS`] (0). **The one count here that is not a
    /// [`NonZeroU32`]**, because zero is this field's shipped setting rather than a mistake.
    pub max_rounds: u32,
    /// Pseudo-counts pulling the re-fitted direction split and fall-off back toward the
    /// stratum's values. Zero is the free setting.
    ///
    /// Defaults to [`DEFAULT_DIRECTION_AND_FALL_OFF_PULL_BACK_PSEUDOCOUNTS`] (50).
    pub direction_and_fall_off_pull_back_pseudocounts: f64,
    /// Slipped reads pulling the re-fitted level back toward the fitted curve's value at the
    /// locus's cell. Zero is the free setting.
    ///
    /// Defaults to [`DEFAULT_LEVEL_PULL_BACK_SLIPPED_READS`] (20).
    pub level_pull_back_slipped_reads: f64,
    /// Stop when every re-fitted number moves less than this between rounds.
    ///
    /// Defaults to [`DEFAULT_ROUND_CONVERGENCE_THRESHOLD`] (`1e-3`).
    pub round_convergence_threshold: f64,
}

impl SlippageRefitConfig {
    /// Frozen: the numbers the parameter fit measured for the stratum, unchanged.
    pub const DEFAULT: Self = Self {
        max_rounds: DEFAULT_SLIPPAGE_REFIT_ROUNDS,
        direction_and_fall_off_pull_back_pseudocounts:
            DEFAULT_DIRECTION_AND_FALL_OFF_PULL_BACK_PSEUDOCOUNTS,
        level_pull_back_slipped_reads: DEFAULT_LEVEL_PULL_BACK_SLIPPED_READS,
        round_convergence_threshold: DEFAULT_ROUND_CONVERGENCE_THRESHOLD,
    };

    /// Whether any round would run. False is the shipped setting.
    #[inline]
    #[must_use]
    pub fn is_frozen(&self) -> bool {
        self.max_rounds == 0
    }
}

impl Default for SlippageRefitConfig {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Whether a locus may gain alleles its selection missed.
///
/// A repeat tract can hide a real allele under stutter: a short allele in a sample that also
/// carries a long one looks exactly like a contraction slip, and a diploid sample's nomination
/// puts forward at most two lengths, so a third never reaches the table however many reads
/// carry it (spec §4.1).
///
/// **Discovery runs before the frequency loop, inside candidate selection**
/// ([`BeforeTheLoop`](Self::BeforeTheLoop), wired in
/// [`select_ssr`](crate::calling::allele_candidates::ssr::select_ssr)). Spec §4.1 had it look
/// against the converged posteriors, in rounds around the loop. Building it
/// showed that ng's tract locus generator realigns every read before the caller sees it, so
/// an observation's bases *are* the tract sequence a retrace would imply: the eligible set is
/// a function of the observations and the candidate table alone, and no posterior enters it.
/// A round around the loop would therefore admit nothing the pre-pass had not, and the two
/// round-wrapped settings were removed
/// (`doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md` §6.5).
///
/// No variant is the default: [`DiscoveryConfig::DEFAULT`] is the one place that says which
/// setting ships.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum DiscoveryMode {
    /// No discovery. The alleles a locus is called over are what candidate selection settled
    /// on, and they do not change during the loop. A measurement arm (`NG_TRACT_DISCOVERY=0`).
    Off,
    /// **The shipped setting.** Discover before the loop runs at all, inside candidate
    /// selection: a tract sequence one sample showed with reads clearing the bar is admitted
    /// as a candidate, its reads counted in `AD` rather than the per-sample leftover, and the
    /// loop then runs once over the widened table.
    BeforeTheLoop,
}

/// Discovering alleles from the calling: whether it runs, and how much evidence a length
/// needs.
///
/// **On by default** since the owner adopted it on 2026-09-03. What it costs is one walk over
/// the merge's rows per tract in candidate selection, not a second convergence. How often an
/// allele is actually hidden under stutter is still unmeasured on this project's data
/// (spec §4, §12's question 3; research doc §6.5).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct DiscoveryConfig {
    /// Whether discovery runs at all. Defaults to
    /// [`DiscoveryMode::BeforeTheLoop`] — **the shipped setting since the owner's adoption
    /// of 2026-09-03** (the tract-accuracy program's L7,
    /// `doc/devel/ng/research/tract_accuracy_program_report.md`): measured at +0.11/+0.04
    /// points at 30× and +0.13/+0.03 at 50× with the spurious-heterozygote class enlarged by
    /// 0 and 2 tracts, and behaviourally inert on the 63-accession tomato cohort at ~3 reads
    /// (2 records of 228,852 changed, both coherent admissions).
    pub mode: DiscoveryMode,
    /// How much evidence one sample must show before a length is admitted. Defaults to
    /// [`DEFAULT_DISCOVERY_BAR`] — 2 reads or 15% of that sample's tract-spanning reads,
    /// whichever is more.
    pub bar: MinAltReads,
}

impl DiscoveryConfig {
    /// The pre-pass, at the shipped bar — what every run gets since the L7 adoption; a run
    /// that wants discovery off says so (`NG_TRACT_DISCOVERY=0` until the parameters-file
    /// key this switch owes exists).
    pub const DEFAULT: Self = Self {
        mode: DiscoveryMode::BeforeTheLoop,
        bar: DEFAULT_DISCOVERY_BAR,
    };

    #[cfg(test)]
    /// Whether discovery is switched off. False is the shipped setting.
    #[inline]
    #[must_use]
    pub fn is_off(&self) -> bool {
        matches!(self.mode, DiscoveryMode::Off)
    }

    /// **The bar candidate selection's discovery pre-pass applies, or `None` where discovery
    /// is off** — what the run driver hands to
    /// [`SsrSelectionConfig::discovery`](crate::calling::allele_candidates::ssr::SsrSelectionConfig::discovery).
    #[inline]
    #[must_use]
    pub fn pre_pass_bar(&self) -> Option<MinAltReads> {
        match self.mode {
            DiscoveryMode::BeforeTheLoop => Some(self.bar),
            DiscoveryMode::Off => None,
        }
    }
}

impl Default for DiscoveryConfig {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Every switch the calling loop has, as **values**.
///
/// **Nothing here is a code path**, and that is the point: the slippage re-fit's three
/// settings are configurations of one implementation, so the shipped behaviour and the
/// measured alternatives cannot drift apart (spec §5.1).
///
/// **This type is what a run builds; it is not what the loop takes.** The loop takes a
/// [`RunnableCallingLoopConfig`], which [`Self::validate`] is the only way to make — so a
/// value out of range is refused once, where it is built, rather than depending on every
/// call site to remember a check.
///
/// The two frequency-loop constants are inherited from production and neither has been
/// measured against this caller's range; both are marked soft where they are defined.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct CallingLoopConfig {
    /// Stop when the largest change in the cohort's expected allele copies **over the
    /// cohort's chromosomes** falls below this — see [`DEFAULT_CONVERGENCE_THRESHOLD`] for
    /// why the division is load-bearing, and [`CONVERGENCE_THRESHOLD_RANGE_MAX`] for why
    /// there is a ceiling.
    ///
    /// Defaults to [`DEFAULT_CONVERGENCE_THRESHOLD`] (`1e-3`).
    pub convergence_threshold: f64,
    /// Emit with [`LocusInference::converged`] false after this many passes.
    ///
    /// Defaults to [`DEFAULT_MAX_PASSES`] (50). A [`NonZeroU32`], because a cap of zero
    /// passes would call no genotype anywhere.
    pub max_passes: NonZeroU32,
    /// Whether the locus may re-fit its own slippage numbers. Frozen by default.
    pub slippage_refit: SlippageRefitConfig,
    /// Whether candidate selection runs the STR discovery pre-pass. On by default
    /// ([`DiscoveryConfig::DEFAULT`]).
    pub discovery: DiscoveryConfig,
}

impl CallingLoopConfig {
    /// What ships: the frequency loop with the slippage round inert, after candidate
    /// selection's discovery pre-pass.
    pub const DEFAULT: Self = Self {
        convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
        max_passes: non_zero_default(DEFAULT_MAX_PASSES),
        slippage_refit: SlippageRefitConfig::DEFAULT,
        discovery: DiscoveryConfig::DEFAULT,
    };

    /// Check this configuration once, and hand back the token the loop takes.
    ///
    /// **The refusals are the point of the step, so they are not skippable.**
    /// [`LocusGenotyper::call_locus`] takes a [`RunnableCallingLoopConfig`] and this is its
    /// only constructor, so a configuration this caller will not honour cannot reach a loop
    /// at all — rather than depending on a caller remembering to ask.
    ///
    /// **What it refuses is a value outside its range** — a convergence threshold that is not
    /// a fraction on the frequency scale, a negative pull-back. Fewer of these than there
    /// were: a floor of zero reads, a share outside `[0, 1]` and a pass cap of zero are no
    /// longer refused here because they are no longer *expressible*, held instead by
    /// [`MinAltReads`] and by [`NonZeroU32`].
    ///
    /// # Errors
    ///
    /// [`CallingLoopConfigError`], which names which setting and what it was.
    pub fn validate(self) -> Result<RunnableCallingLoopConfig, CallingLoopConfigError> {
        if !(self.convergence_threshold.is_finite()
            && self.convergence_threshold > 0.0
            && self.convergence_threshold <= CONVERGENCE_THRESHOLD_RANGE_MAX)
        {
            return Err(CallingLoopConfigError::ConvergenceThresholdOutOfRange {
                threshold: self.convergence_threshold,
            });
        }
        let refit = &self.slippage_refit;
        for (pull_back, which) in [
            (
                refit.direction_and_fall_off_pull_back_pseudocounts,
                PullBack::DirectionAndFallOff,
            ),
            (refit.level_pull_back_slipped_reads, PullBack::Level),
        ] {
            if !(pull_back.is_finite() && pull_back >= 0.0) {
                return Err(CallingLoopConfigError::PullBackOutOfRange { which, pull_back });
            }
        }
        if !(refit.round_convergence_threshold.is_finite()
            && refit.round_convergence_threshold > 0.0)
        {
            return Err(
                CallingLoopConfigError::RoundConvergenceThresholdOutOfRange {
                    threshold: refit.round_convergence_threshold,
                },
            );
        }
        Ok(RunnableCallingLoopConfig(self))
    }
}

impl Default for CallingLoopConfig {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// A configuration this caller has agreed to run.
///
/// **It carries no state of its own — it is the same configuration with a promise attached**,
/// and the promise is that [`CallingLoopConfig::validate`] accepted it. That constructor is
/// the only one, so a loop handed one of these cannot be running a value out of range, and no
/// implementation has to remember to check.
///
/// The wrapped settings are read straight through, so loop code writes
/// `config.convergence_threshold` as it would on the plain configuration.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct RunnableCallingLoopConfig(CallingLoopConfig);

impl RunnableCallingLoopConfig {
    #[cfg(test)]
    /// The settings, unwrapped.
    #[inline]
    #[must_use]
    pub fn settings(&self) -> CallingLoopConfig {
        self.0
    }
}

impl std::ops::Deref for RunnableCallingLoopConfig {
    type Target = CallingLoopConfig;

    #[inline]
    fn deref(&self) -> &CallingLoopConfig {
        &self.0
    }
}

impl Default for RunnableCallingLoopConfig {
    /// What ships, which needs no fallible path: that it validates is pinned by
    /// `the_shipped_configuration_is_one_this_caller_will_run`.
    fn default() -> Self {
        Self(CallingLoopConfig::DEFAULT)
    }
}

/// A calling-loop configuration this caller will not run.
///
/// **A configuration is a run's request, not a caller bug**, which is why these are returned
/// rather than asserted — everything else in `calling/` panics, because everything else it
/// refuses is a wiring mistake (`doc/devel/ng/spec/calling_em_loop.md` §8).
///
/// `PartialEq` but not `Eq`: every variant carries the offending `f64` verbatim,
/// including the not-a-number that may be why it was refused.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum CallingLoopConfigError {
    /// `convergence_threshold` is not a fraction on the allele-frequency scale.
    #[error(
        "`convergence_threshold` is a fraction on the allele-frequency scale and must be \
         finite and within (0, {CONVERGENCE_THRESHOLD_RANGE_MAX}], not {threshold}; the \
         loop's movement is already divided by the cohort's chromosomes, so a threshold of 1 \
         stops every locus after one pass and flags it as settled"
    )]
    ConvergenceThresholdOutOfRange { threshold: f64 },
    /// One of the two pull-backs in `slippage_refit` is not a weight.
    #[error(
        "the {which} pull-back is a weight — pseudo-counts, or slipped reads — so it must be \
         finite and at or above zero, not {pull_back}; zero is the free setting, where the \
         locus's own reads set the numbers outright"
    )]
    PullBackOutOfRange { which: PullBack, pull_back: f64 },
    /// `slippage_refit.round_convergence_threshold` is not a positive finite number.
    #[error(
        "`slippage_refit.round_convergence_threshold` must be finite and above zero, not \
         {threshold}"
    )]
    RoundConvergenceThresholdOutOfRange { threshold: f64 },
}

/// **One locus, the whole cohort, in and calls out** — the one boundary every way of handling
/// a cohort crosses.
///
/// Everything the alternatives share goes through here: what each sample's reads showed, what
/// the parameter pre-pass froze, which alleles the locus is called over, how the loop is
/// configured, and the worker's scratch. What differs is what happens inside, and the design
/// has three answers to compare (`doc/devel/ng/spec/calling_em_loop.md` §12, question 1):
///
/// - **summarise and condition** — estimate one cohort-wide set of allele frequencies and
///   score each sample against it, subtracting that sample's own contribution back out so its
///   reads cannot count twice. This is the default and what the rest of the design describes.
/// - **score whole assignments** — take one genotype for every sample as a single object and
///   score it under a prior over the whole cohort's allele counts, which needs no such
///   subtraction because it never estimates a frequency. Two of these, differing in the
///   prior.
///
/// **The type parameter is the repeat-tract emission model's own working memory**, carried
/// because the scratch carries it; an implementation is generic over it, since the loop hands
/// it straight to the row builder and never looks inside. It is a parameter rather than an
/// associated type because an associated type does not compile here: an implementation that
/// works for every model constrains it nowhere, which is `error[E0207]`.
///
/// **`candidates` is taken by value**, so that the loop can hand the locus's allele table back
/// inside the [`LocusInference`] without copying it. The loop does not change the table.
pub trait LocusGenotyper<SsrEmissionScratch> {
    /// Call one locus.
    ///
    /// The evidence, the alleles and the parameters must agree with each other — one run-wide
    /// sample order indexes them all, and the evidence's path must match the allele table's
    /// kind. [`LocusEvidence::assert_matches_locus_and_run`] and
    /// [`CallingScratch::prepare_for_locus`] are where those are checked, and an
    /// implementation calls both rather than trusting its caller. **The configuration needs
    /// no such call**: it arrives already checked, because
    /// [`CallingLoopConfig::validate`] is the only way to make the type this takes.
    fn call_locus(
        &self,
        evidence: &LocusEvidence<'_>,
        parameters: &FrozenParameters<'_>,
        candidates: CandidateAlleles,
        config: &RunnableCallingLoopConfig,
        scratch: &mut CallingScratch<SsrEmissionScratch>,
    ) -> LocusInference;

    #[cfg(test)]
    /// Which way of handling the cohort this is, for a run to record beside the genotypes it
    /// produced.
    ///
    /// **This seam exists to compare three answers, and a result that cannot say which one
    /// produced it is not auditable.** The genotype prior carries the same method for the same
    /// reason ([`GenotypePriorModel::name`](crate::calling::genotype_prior::GenotypePriorModel::name)):
    /// the measurement this seam is built for is a difference between arms, and an arm without
    /// a label is a number nobody can act on.
    fn name(&self) -> &'static str;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calling::{ExpectedAlleleCopies, SampleGenotypeCall};
    use crate::locus_generation::LocusKind;
    use crate::parameter_estimation::Provenance;
    use crate::types::{AlleleId, ContigId, GenomeRegion, Genotype, Phred, Position};

    /// What ships runs, and every value it ships with is the constant that names it.
    ///
    /// **This is the test that would fail if a default and its validation drifted apart**,
    /// which is the failure mode a validated configuration with a `Default` impl always has:
    /// the shipped configuration is the one nobody passes explicitly, so it is the one a
    /// range check is least likely to be tried against.
    #[test]
    fn the_shipped_configuration_is_one_this_caller_will_run() {
        let config = CallingLoopConfig::default();
        assert!(config.validate().is_ok());
        assert_eq!(config, CallingLoopConfig::DEFAULT);
        assert_eq!(*RunnableCallingLoopConfig::default(), config);

        // **Compared against the numbers, not against the constants they came from.**
        // `assert_eq!(field, DEFAULT_X)` is an identity — `Default` read `DEFAULT_X` a
        // moment earlier — so it holds for whatever value the constant is edited to.
        // Measured: written that way, changing the round threshold to 0.5 or the discovery
        // share to 0.9 each left the whole suite green.
        // These are the numbers the design's open questions will move deliberately, and a
        // deliberate move should have to edit a test that states the old value.
        assert_eq!(config.convergence_threshold, 1e-3);
        assert_eq!(config.max_passes.get(), 50);
        assert!(config.slippage_refit.is_frozen());
        assert_eq!(config.slippage_refit.max_rounds, 0);
        assert_eq!(
            config
                .slippage_refit
                .direction_and_fall_off_pull_back_pseudocounts,
            50.0
        );
        assert_eq!(config.slippage_refit.level_pull_back_slipped_reads, 20.0);
        assert_eq!(config.slippage_refit.round_convergence_threshold, 1e-3);
        assert!(!config.discovery.is_off());
        assert_eq!(config.discovery.mode, DiscoveryMode::BeforeTheLoop);

        // **The pre-pass is the shipped mode since the L7 adoption (2026-09-03).**
        assert_eq!(config.discovery.bar.floor.get(), 2); // HipSTR's read floor
        assert_eq!(config.discovery.bar.share.get(), 0.15); // HipSTR's share of one sample
    }

    /// **The shipped default is not the one from next door**, and the two must be free to
    /// move apart.
    ///
    /// Discovery's bar and candidate selection's are the same *rule* — a read floor, or a
    /// share of one sample's own reads, whichever is more — and share a type for that reason.
    /// They are not the same *numbers*: spec §4.1's third open question sweeps discovery's
    /// pair independently of the rule that decides whether a locus is built at all. This
    /// pins that they are separately stated, so a sweep of one cannot move the other.
    #[test]
    fn discovery_and_candidate_selection_share_a_rule_and_not_its_numbers() {
        use crate::calling::allele_candidates::DEFAULT_MIN_ALLELE_SUPPORT;

        assert_ne!(DEFAULT_DISCOVERY_BAR, DEFAULT_MIN_ALLELE_SUPPORT);
        assert_ne!(
            DEFAULT_DISCOVERY_BAR.share.get(),
            DEFAULT_MIN_ALLELE_SUPPORT.share.get()
        );
    }

    /// The built outer settings validate, alone and together: re-fit rounds, discovery's
    /// pre-pass, and the pair a measurement arm would set at once.
    #[test]
    fn the_built_outer_settings_validate() {
        let mut refitting = CallingLoopConfig::default();
        refitting.slippage_refit.max_rounds = 3;
        assert!(refitting.validate().is_ok());

        let mut discovering_before_the_loop = CallingLoopConfig::default();
        discovering_before_the_loop.discovery.mode = DiscoveryMode::BeforeTheLoop;
        assert!(discovering_before_the_loop.validate().is_ok());

        let mut not_discovering = CallingLoopConfig::default();
        not_discovering.discovery.mode = DiscoveryMode::Off;
        assert!(not_discovering.validate().is_ok());

        let mut both = refitting;
        both.discovery.mode = DiscoveryMode::BeforeTheLoop;
        assert!(both.validate().is_ok());
    }

    /// **What the driver reads off a validated configuration**: the pre-pass mode hands
    /// selection its bar, and off hands it nothing.
    #[test]
    fn only_the_pre_pass_mode_hands_selection_a_bar() {
        let mut config = DiscoveryConfig::DEFAULT;
        assert_eq!(
            config.pre_pass_bar(),
            Some(DEFAULT_DISCOVERY_BAR),
            "the shipped default IS the pre-pass since the L7 adoption, carrying \
             discovery's own bar rather than selection's"
        );

        config.mode = DiscoveryMode::Off;
        assert_eq!(config.pre_pass_bar(), None, "off runs no pre-pass");
    }

    /// Every range check refuses the value it is for, and each names a different variant.
    ///
    /// **The not-a-number cases are the ones a comparison-only check misses.** `NaN` fails
    /// every ordinary comparison, so a check written as `threshold <= 0.0` lets it through
    /// and the loop's stopping test can then never be satisfied — the locus runs to its pass
    /// cap and is emitted unconverged, with nothing saying why.
    ///
    /// **And the upper bound is not decoration.** The loop's movement is a per-allele change
    /// already divided by the cohort's chromosomes, so it cannot exceed 1: at a threshold of
    /// 1 every locus in the run stops after one pass and is emitted flagged as settled, which
    /// is a whole run of confident wrong claims with nothing failing.
    #[test]
    fn a_value_outside_its_range_is_refused_by_the_check_it_belongs_to() {
        for threshold in [0.0, -1e-3, f64::NAN, f64::INFINITY, 1.0, 0.5, 1e9] {
            let config = CallingLoopConfig {
                convergence_threshold: threshold,
                ..CallingLoopConfig::default()
            };
            assert!(
                matches!(
                    config.validate(),
                    Err(CallingLoopConfigError::ConvergenceThresholdOutOfRange { .. })
                ),
                "a threshold of {threshold} was accepted"
            );
        }
        // The ceiling itself is legal — it is loose, not degenerate.
        let at_the_ceiling = CallingLoopConfig {
            convergence_threshold: CONVERGENCE_THRESHOLD_RANGE_MAX,
            ..CallingLoopConfig::default()
        };
        assert!(at_the_ceiling.validate().is_ok());

        let mut negative_shape = CallingLoopConfig::default();
        negative_shape
            .slippage_refit
            .direction_and_fall_off_pull_back_pseudocounts = -1.0;
        assert!(matches!(
            negative_shape.validate(),
            Err(CallingLoopConfigError::PullBackOutOfRange {
                which: PullBack::DirectionAndFallOff,
                ..
            })
        ));

        // **Both halves of every finiteness check, and infinity is the half a fixture
        // reaches last.** `is_finite()` weakened to `!is_nan()` still refuses a
        // not-a-number and admits an infinity, so a not-a-number case alone leaves that
        // weakening alive — measured: it survived on both pull-backs and on the round
        // threshold until these cases were added.
        for level in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
            let mut bad_level = CallingLoopConfig::default();
            bad_level.slippage_refit.level_pull_back_slipped_reads = level;
            assert!(
                matches!(
                    bad_level.validate(),
                    Err(CallingLoopConfigError::PullBackOutOfRange {
                        which: PullBack::Level,
                        ..
                    })
                ),
                "a level pull-back of {level} was accepted"
            );
        }

        for threshold in [0.0, -1e-3, f64::NAN, f64::INFINITY] {
            let mut bad_round = CallingLoopConfig::default();
            bad_round.slippage_refit.round_convergence_threshold = threshold;
            assert!(
                matches!(
                    bad_round.validate(),
                    Err(CallingLoopConfigError::RoundConvergenceThresholdOutOfRange { .. })
                ),
                "a round threshold of {threshold} was accepted"
            );
        }
    }

    /// The number a refusal carries is the number that was refused, and it reaches the
    /// sentence a run logs.
    ///
    /// **Nothing checked this before.** Every test reaching an `f64`-carrying variant matched
    /// it with `{ .. }`, which reads no payload, so an implementation answering every refusal
    /// with `ConvergenceThresholdOutOfRange { threshold: 0.0 }` passed the whole suite — and
    /// a refusal naming the wrong number sends an operator to the wrong field. `thiserror`
    /// checks that the format string compiles, not that the field it names is the offending
    /// one.
    #[test]
    fn a_refused_value_travels_in_the_error_and_into_its_message() {
        let over_the_ceiling = CallingLoopConfig {
            convergence_threshold: 0.25,
            ..CallingLoopConfig::default()
        };
        let refusal = over_the_ceiling.validate().expect_err("above the ceiling");
        assert_eq!(
            refusal,
            CallingLoopConfigError::ConvergenceThresholdOutOfRange { threshold: 0.25 }
        );
        assert!(refusal.to_string().contains("0.25"), "{refusal}");

        let mut negative = CallingLoopConfig::default();
        negative
            .slippage_refit
            .direction_and_fall_off_pull_back_pseudocounts = -7.5;
        let refusal = negative.validate().expect_err("a negative weight");
        assert_eq!(
            refusal,
            CallingLoopConfigError::PullBackOutOfRange {
                which: PullBack::DirectionAndFallOff,
                pull_back: -7.5
            }
        );
        assert!(refusal.to_string().contains("-7.5"), "{refusal}");
    }

    /// **Three things `validate` used to refuse are now values nobody can write**, and that is
    /// a stronger guarantee than a check.
    ///
    /// A pass cap of zero, a discovery bar of zero reads and a share outside `[0, 1]` were
    /// range checks and error variants of their own. They are now
    /// held by [`NonZeroU32`] and by [`MinAltReads`]'s own constructors, so the illegal value
    /// never exists to be checked — which no mutation of `validate` can undo.
    #[test]
    fn the_counts_and_the_bar_refuse_their_illegal_values_at_construction() {
        assert_eq!(NonZeroU32::new(0), None);
        assert_eq!(MinAltReadShare::new(-0.01), None);
        assert_eq!(MinAltReadShare::new(1.01), None);
        assert_eq!(MinAltReadShare::new(f64::NAN), None);
        assert_eq!(MinAltReadShare::new(f64::INFINITY), None);

        // And the shipped bar answers the question discovery asks of one sample:
        // 2 reads, or 15% of that sample's tract-spanning reads, whichever is more.
        assert_eq!(DEFAULT_DISCOVERY_BAR.required_of(3), 2);
        assert_eq!(DEFAULT_DISCOVERY_BAR.required_of(100), 15);
    }

    /// **Zero pull-back is a setting, not an error** — it is the free arm, where the locus's
    /// own reads set its slippage numbers outright, and a range check that refused it would
    /// make one of the three settings the design compares unreachable.
    ///
    /// The shape it is checked in is the one a measurement uses: free **with rounds on**,
    /// which is HipSTR's setting and now a configuration this caller runs.
    #[test]
    fn zero_pull_back_is_the_free_setting_and_passes_the_range_check() {
        let mut free = CallingLoopConfig::default();
        free.slippage_refit
            .direction_and_fall_off_pull_back_pseudocounts = 0.0;
        free.slippage_refit.level_pull_back_slipped_reads = 0.0;
        assert!(free.validate().is_ok());

        // And with rounds on — the shape a measurement arm actually sets — it still runs.
        free.slippage_refit.max_rounds = 1;
        assert!(free.validate().is_ok());
    }

    /// The pull-back an operator reads about in a log line is named in words, not as a Rust
    /// identifier.
    #[test]
    fn a_refusal_reads_as_a_sentence_rather_than_as_an_identifier() {
        let mut negative = CallingLoopConfig::default();
        negative
            .slippage_refit
            .direction_and_fall_off_pull_back_pseudocounts = -1.0;
        let sentence = negative.validate().expect_err("out of range").to_string();
        assert!(
            sentence.contains("direction-and-fall-off pull-back"),
            "{sentence}"
        );
    }

    /// A site quality standing in for one the worker computed — this stand-in looks at no
    /// evidence, so it has none to compute one from.
    fn a_worker_written_site_quality() -> Phred {
        Phred::try_new(37.0).expect("a legal quality, and below the site-quality ceiling")
    }

    /// A stand-in implementation, existing only so the seam is exercised by something. It
    /// calls one genotype for every sample and looks at no evidence at all — deliberately not
    /// a caller, so nothing here can be mistaken for a check of one.
    struct EveryoneIsHomozygousReference;

    impl<S> LocusGenotyper<S> for EveryoneIsHomozygousReference {
        fn name(&self) -> &'static str {
            "everyone is homozygous reference (test stand-in)"
        }

        fn call_locus(
            &self,
            evidence: &LocusEvidence<'_>,
            parameters: &FrozenParameters<'_>,
            candidates: CandidateAlleles,
            _config: &RunnableCallingLoopConfig,
            scratch: &mut CallingScratch<S>,
        ) -> LocusInference {
            evidence.assert_matches_locus_and_run(&candidates, parameters);
            let table = crate::calling::GenotypeTable::build(parameters.ploidy(), candidates.len());
            scratch.prepare_for_locus(evidence.sample_count(), &candidates, &table.view());

            let copies = vec![
                f64::from(parameters.ploidy().get()) * evidence.sample_count() as f64,
                0.0,
            ];
            let per_sample = (0..evidence.sample_count())
                .map(|_| SampleGenotypeCall::Called {
                    genotype: Genotype::new(vec![AlleleId::REFERENCE; 2]),
                    genotype_quality: Phred::try_new(30.0).expect("a legal quality"),
                    // This stub calls every sample homozygous reference without reading a
                    // likelihood, so it claims the reads said something rather than
                    // pretending to have checked.
                    reads_were_uninformative: false,
                    genotype_rows: crate::calling::placeholder_genotype_rows(),
                })
                .collect();
            let expected = ExpectedAlleleCopies::new(copies, &candidates);
            LocusInference::new(
                evidence.region(),
                candidates,
                per_sample,
                expected,
                true,
                1,
                Provenance::FittedHere,
                None,
                a_worker_written_site_quality(),
                None,
            )
        }
    }

    /// The seam takes everything an arm needs and is object-safe, so a run can choose one at
    /// run time and record which it chose.
    ///
    /// **Object safety is the property worth pinning.** The design compares three ways of
    /// handling a cohort; if the trait could not be held behind a `Box`, choosing between them
    /// would be a compile-time decision and the comparison would need three binaries.
    ///
    /// **And the configuration it takes cannot be one this caller refused** — the seam takes
    /// the token, and validation is the only way to make one.
    #[test]
    fn the_seam_is_object_safe_and_every_arm_can_name_itself() {
        let arm: Box<dyn LocusGenotyper<()>> = Box::new(EveryoneIsHomozygousReference);
        assert!(!arm.name().is_empty());

        let mut alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
        alleles.admit(Box::from(b"T".as_slice()));
        let calibration = [crate::calling::ReadGroupCalibration::unscaled()];
        let inbreeding = [crate::types::InbreedingF::try_new(0.0).expect("legal")];
        let strata = crate::parameter_estimation::joint::stratum_fits::StratumFits::over(
            &[],
            std::collections::BTreeMap::new(),
        );
        let substitution = std::collections::BTreeMap::new();
        let parameters = FrozenParameters::uncontaminated(
            &calibration,
            &inbreeding,
            crate::calling::genotype_prior::SpectrumSeed::new(
                1.0,
                1e-3,
                crate::calling::genotype_prior::SeedRegime::NeutralShape,
            ),
            &strata,
            &substitution,
            crate::types::Ploidy::try_new(2).expect("a diploid"),
        );
        let per_sample = [crate::calling::GenericLocusSample {
            evidence: crate::calling::GenericSampleEvidence::empty(),
            genotype_must_be_missing: false,
        }];
        let region = GenomeRegion {
            contig: ContigId(0),
            start: Position(10),
            end: Position(10),
        };
        let evidence = LocusEvidence::generic(region, &per_sample);
        let mut scratch = CallingScratch::<()>::default();
        let config = CallingLoopConfig::default()
            .validate()
            .expect("the shipped configuration runs");

        let called = arm.call_locus(&evidence, &parameters, alleles, &config, &mut scratch);
        assert_eq!(called.region, region);
        assert_eq!(called.per_sample.len(), 1);
        assert_eq!(called.cohort_expected_copies().copies(), [2.0, 0.0]);

        // The token reads straight through to the settings the loop needs.
        assert_eq!(config.max_passes.get(), DEFAULT_MAX_PASSES);
        assert_eq!(config.settings(), CallingLoopConfig::DEFAULT);
    }
}
