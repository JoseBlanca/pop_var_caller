//! Step 9 — turning one locus's evidence and the run's frozen parameters into genotypes.
//!
//! **This file is the seam and its configuration; the loops themselves arrive beside it.**
//! [`LocusGenotyper`] is the one boundary everything the calling arms share crosses — the
//! evidence, the parameters, the candidate alleles and the worker's scratch — and
//! [`CallingLoopConfig`] is every switch of the design, held as a **value** rather than as a
//! code path.
//!
//! ## Three loops, one inside the next, and two of them ship switched off
//!
//! The innermost loop is the one that runs: guess the allele frequencies, work out each
//! sample's genotype probabilities under that guess, add them up to get better frequencies,
//! repeat until they stop moving. Around it sit two more, and both default to a single pass
//! through their bodies:
//!
//! - **the slippage round**, which re-fits how often a read gains or loses a whole repeat
//!   from the locus's own reads and runs the frequency loop again on the new numbers
//!   ([`SlippageRefitConfig`], `max_rounds` 0);
//! - **the discovery round**, which looks at what the converged answer is explaining as
//!   slippage and admits the tract lengths that recur too often in one sample to be slippage
//!   ([`DiscoveryConfig`], [`DiscoveryMode::Off`]).
//!
//! Both are built as configurations of one code path rather than as three implementations,
//! which is what makes *frozen*, *pulled back part way* and *free* three values instead of
//! three code paths to keep in step (`doc/devel/ng/spec/calling_em_loop.md` §5.1).
//!
//! **The slippage round's body is built** — the arithmetic in [`slippage_refit`], the round
//! itself in the loop's driver — and a non-zero `max_rounds` now validates. **Discovery's is
//! not, and a run that asks for it is still refused** — see [`CallingLoopConfig::validate`].
//! The refusal is not advisory: [`LocusGenotyper::call_locus`] takes a
//! [`RunnableCallingLoopConfig`], which has one constructor and it is the check, so a setting
//! this caller will not honour cannot reach the loop at all.
//!
//! ## What is deliberately not configured here
//!
//! **The allele cap.** How many alleles a locus may be called over belongs to candidate
//! selection, which owns it together with the support bar
//! ([`CandidateSelectionConfig`](crate::ng::calling::allele_candidates::CandidateSelectionConfig)),
//! and states it as a
//! [`MaxCandidateAlleles`](crate::ng::calling::allele_candidates::MaxCandidateAlleles) — a
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

pub mod discovery;
pub mod repeat_tract_parameters;
pub(crate) mod slippage_refit;
pub mod summarise_condition;

use std::num::NonZeroU32;

use crate::ng::calling::{
    CallingScratch, CandidateAlleles, FrozenParameters, LocusEvidence, LocusInference,
};
use crate::ng::run::cohort_merge::{MinAltObs, MinAltReadShare, MinAltReads};

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
/// (`src/var_calling/posterior_engine.rs`). Nothing has measured it on this caller's range —
/// a thousand samples, or three reads a position — and what would set it is the distribution
/// of [`LocusInference::passes`] across a real panel (spec §12, question 4).
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
/// `src/ssr/cohort/em.rs`). ng holds the two at one number until something measures them
/// apart.
pub const DEFAULT_ROUND_CONVERGENCE_THRESHOLD: f64 = 1e-3;

/// How many of one sample's reads must show a tract length before a discovery round may admit
/// it as an allele.
///
/// **Both halves of the shipped bar must clear, so that a single stray read cannot mint an
/// allele.** Inherited from HipSTR's high-depth human setting and **soft** — the two halves
/// bind at opposite ends of the depth range, and below about 13 reads only this one does
/// anything, because 2 reads already clears 15% (spec §4.1). That is a property of the
/// shipped numbers rather than a guarantee of the type: [`MinAltReads`] admits a share of
/// zero, and spec §4.1's third open question sweeps both halves.
pub const DEFAULT_DISCOVERY_MIN_READS: u32 = 2;

/// What share of one sample's **tract-spanning** reads must show a length before a discovery
/// round may admit it. HipSTR's, inherited and soft — see [`DEFAULT_DISCOVERY_MIN_READS`] for
/// why the pair cannot be right at both ends of the depth range.
///
/// **Of one sample's reads, not of the cohort's.** HipSTR's other admission route is
/// cohort-wide — above 5% of samples, or 5% of reads — and ng builds neither, so the two
/// readings are live and this is the one that applies.
pub const DEFAULT_DISCOVERY_MIN_SPANNING_READ_SHARE: f64 = 0.15;

/// How much evidence one sample must show before a discovery round admits a tract length:
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

/// How many discovery rounds a locus may run before it stops looking.
///
/// A round that admits nothing ends the loop on its own, so this is the runaway guard rather
/// than the usual stopping rule: the expected count is one or two, and a locus that keeps
/// finding one more allele is the shape that would hurt (spec §4.1).
///
/// **This number is not inherited from anywhere — it is this step's, and it is soft.** The
/// spec stops a discovery loop on a round that adds nothing or on the allele cap, and names
/// no round cap; the architecture gives the field and no value. Four is twice the expected
/// one or two, which is what a runaway guard should be and not a measurement. It is inert
/// while discovery ships off, and the plan that switches discovery on
/// (`doc/devel/ng/impl_plan/calling_bakeoffs.md`) is what should set it.
pub const DEFAULT_DISCOVERY_MAX_ROUNDS: u32 = 4;

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

/// Whether the loop may **add** alleles a locus's converged answer suggests it missed, and
/// against what.
///
/// A repeat tract can hide a real allele under stutter: a short allele in a sample that also
/// carries a long one looks exactly like a contraction slip. After the loop converges, the
/// model's own attribution can be retraced — where it says *this read slipped*, the tract
/// length that read implies is counted per sample, and a length that recurs past a bar is
/// admitted as a candidate (spec §4.1).
///
/// **The middle setting may be the answer, which is why there are three.** Discovering
/// against frequencies held fixed at what one convergence produced makes each round a single
/// scoring pass rather than a whole convergence, and a final convergence at the end puts back
/// what that gives up.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub enum DiscoveryMode {
    /// **The default.** The alleles a locus is called over are settled before the first pass
    /// and do not change during it.
    #[default]
    Off,
    /// Converge once, hold the frequencies at what that produced, and discover against them
    /// — each round is a scoring pass. Converge once more at the end on the final allele set.
    AgainstFrozenFrequencies,
    /// Converge fully on every round. The most expensive setting, and the one the spec's
    /// pseudocode writes out.
    AgainstFullConvergence,
}

impl std::fmt::Display for DiscoveryMode {
    /// What an operator should read in a log line, rather than the Rust identifier.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Off => "no allele discovery",
            Self::AgainstFrozenFrequencies => {
                "discovering alleles against frozen allele frequencies"
            }
            Self::AgainstFullConvergence => "discovering alleles against a full convergence",
        })
    }
}

/// Discovering alleles from the calling: whether, against what, and how much evidence a
/// length needs.
///
/// Off by default, and the reason is a cost-against-benefit judgement rather than a
/// mechanical objection. **The cost is certain and paid at every locus** — a round is a whole
/// extra run of the loop, because the lengths it looks for are the ones the *converged*
/// answer is explaining as slippage, so a locus where nothing is found still pays to
/// establish that. **The benefit is unmeasured**: nothing on this project's data says how
/// often an allele is actually hidden under stutter (spec §4, §12's question 3).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct DiscoveryConfig {
    /// Whether discovery runs at all, and against what. Defaults to [`DiscoveryMode::Off`].
    pub mode: DiscoveryMode,
    /// How much evidence one sample must show before a length is admitted. Defaults to
    /// [`DEFAULT_DISCOVERY_BAR`] — 2 reads or 15% of that sample's tract-spanning reads,
    /// whichever is more.
    pub bar: MinAltReads,
    /// The runaway guard — a round that admits nothing already ends the loop.
    ///
    /// Defaults to [`DEFAULT_DISCOVERY_MAX_ROUNDS`] (4), **which is this step's own number
    /// and inherited from nowhere** — see the constant for why, and for why it is soft.
    /// A [`NonZeroU32`], so a switched-on discovery loop that would run no rounds and report
    /// finding nothing is not a value anyone can write; [`DiscoveryMode::Off`] is how
    /// discovery is switched off.
    pub max_rounds: NonZeroU32,
}

impl DiscoveryConfig {
    /// Off, with the bar and the runaway guard already set to what a run switching it on
    /// would want.
    pub const DEFAULT: Self = Self {
        mode: DiscoveryMode::Off,
        bar: DEFAULT_DISCOVERY_BAR,
        max_rounds: non_zero_default(DEFAULT_DISCOVERY_MAX_ROUNDS),
    };

    /// Whether any round would run. True is the shipped setting.
    #[inline]
    #[must_use]
    pub fn is_off(&self) -> bool {
        matches!(self.mode, DiscoveryMode::Off)
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
/// settings and discovery's three are configurations of one implementation, so the shipped
/// behaviour and the measured alternatives cannot drift apart (spec §5.1, §4.1).
///
/// **This type is what a run builds; it is not what the loop takes.** The loop takes a
/// [`RunnableCallingLoopConfig`], which [`Self::validate`] is the only way to make — so a
/// setting this caller will not honour is refused once, where it is built, rather than
/// depending on every call site to remember a check.
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
    /// Whether the loop may add alleles it decides the first pass missed. Off by default.
    pub discovery: DiscoveryConfig,
}

impl CallingLoopConfig {
    /// What ships: the frequency loop alone, with both outer rounds inert.
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
    /// Two kinds of refusal, and only the second is permanent.
    ///
    /// **A value outside its range** — a convergence threshold that is not a fraction on the
    /// frequency scale, a negative pull-back. Fewer of these than there were: a floor of zero
    /// reads, a share outside `[0, 1]`, a pass cap of zero and a discovery round cap of zero
    /// are no longer refused here because they are no longer *expressible*, held instead by
    /// [`MinAltReads`] and by [`NonZeroU32`].
    ///
    /// **A setting whose body is not built yet.** Discovery ships as a value with no
    /// implementation behind it; the machinery arrives with the measurements that decide its
    /// default (`doc/devel/ng/impl_plan/calling_bakeoffs.md`). Until then a run that asks for
    /// it is stopped here, because the alternative — accepting the setting and running the
    /// loop without discovery anyway — would report a measurement of an arm that was never
    /// run. **The slippage re-fit's body is built** (the tract-accuracy program's L3,
    /// `doc/devel/ng/research/tract_accuracy_program_report.md`), so a non-zero `max_rounds`
    /// is now a setting this caller honours and validates.
    ///
    /// **Range refusals come first**, so a configuration that is both out of range and
    /// unbuilt hears about the range: that half is the caller's to fix today, where the other
    /// is this caller's to build.
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
        // The one not-yet-built setting comes last, so a configuration that is *both*
        // out of range and unbuilt is told about the range first.
        if !self.discovery.is_off() {
            return Err(CallingLoopConfigError::DiscoveryNotBuilt {
                mode: self.discovery.mode,
            });
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
/// the only one, so a loop handed one of these cannot be running a setting whose body does
/// not exist, and no implementation has to remember to check.
///
/// The wrapped settings are read straight through, so loop code writes
/// `config.convergence_threshold` as it would on the plain configuration.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct RunnableCallingLoopConfig(CallingLoopConfig);

impl RunnableCallingLoopConfig {
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
/// `PartialEq` but not `Eq`: three of the four variants carry the offending `f64` verbatim,
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
    /// Allele discovery was asked for and its body is not built.
    #[error(
        "`discovery.mode` asked for {mode}, which is not built yet; it ships off, and \
         accepting this setting would run the loop without discovery and report it as the \
         discovery arm"
    )]
    DiscoveryNotBuilt { mode: DiscoveryMode },
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
/// **`candidates` is taken by value.** The locus's allele table is not a borrow of the
/// evidence: a discovery round appends to it and the final prune shrinks it, so the loop owns
/// it and hands it back inside the [`LocusInference`].
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

    /// Which way of handling the cohort this is, for a run to record beside the genotypes it
    /// produced.
    ///
    /// **This seam exists to compare three answers, and a result that cannot say which one
    /// produced it is not auditable.** The genotype prior carries the same method for the same
    /// reason ([`GenotypePriorModel::name`](crate::ng::calling::genotype_prior::GenotypePriorModel::name)):
    /// the measurement this seam is built for is a difference between arms, and an arm without
    /// a label is a number nobody can act on.
    fn name(&self) -> &'static str;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::{ExpectedAlleleCopies, SampleGenotypeCall};
    use crate::ng::locus_generation::LocusKind;
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::types::{AlleleId, ContigId, GenomeRegion, Genotype, Phred, Position};

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
        // Measured: written that way, changing the round threshold to 0.5, the discovery
        // share to 0.9 and the discovery round cap to 1 each left the whole suite green.
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
        assert!(config.discovery.is_off());
        assert_eq!(config.discovery.mode, DiscoveryMode::Off);

        // **Off is the mode, not a zeroed guard.** Every other field is already what a run
        // switching discovery on would want, so turning the mode on cannot leave a loop that
        // runs no rounds and reports finding nothing.
        assert_eq!(config.discovery.bar.floor.get(), 2); // HipSTR's read floor
        assert_eq!(config.discovery.bar.share.get(), 0.15); // HipSTR's share of one sample
        assert_eq!(config.discovery.max_rounds.get(), 4); // twice the one or two expected
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
        use crate::ng::calling::allele_candidates::DEFAULT_MIN_ALLELE_SUPPORT;

        assert_ne!(DEFAULT_DISCOVERY_BAR, DEFAULT_MIN_ALLELE_SUPPORT);
        assert_ne!(
            DEFAULT_DISCOVERY_BAR.share.get(),
            DEFAULT_MIN_ALLELE_SUPPORT.share.get()
        );
    }

    /// The one setting whose body is not built — discovery — is refused by name, **not run as
    /// the default**; the slippage re-fit's body is built, so a non-zero round count now
    /// validates.
    ///
    /// **Silently honouring the default instead is the failure the refusal exists to stop**,
    /// and it is worse than a crash: a measurement harness would set discovery on, get the
    /// plain loop's answers back, and report them as the discovery arm's. The two arms would
    /// then agree exactly, which reads as a finding.
    #[test]
    fn a_setting_whose_body_is_not_built_is_refused_rather_than_ignored() {
        // The re-fit's rounds are the setting that used to be refused here
        // (`SlippageRefitNotBuilt`, retired when the body was built): asking for three rounds
        // is now a configuration this caller runs.
        let mut refitting = CallingLoopConfig::default();
        refitting.slippage_refit.max_rounds = 3;
        assert!(refitting.validate().is_ok());

        for mode in [
            DiscoveryMode::AgainstFrozenFrequencies,
            DiscoveryMode::AgainstFullConvergence,
        ] {
            let mut discovering = CallingLoopConfig::default();
            discovering.discovery.mode = mode;
            assert_eq!(
                discovering.validate(),
                Err(CallingLoopConfigError::DiscoveryNotBuilt { mode })
            );
        }
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

    /// With both outer loops asked for at once, discovery — the one whose body is not built —
    /// is the refusal reported; the re-fit's rounds no longer stand in its way.
    ///
    /// **This is the case a measurement harness hits first**, since a bake-off arm sets both.
    /// While both bodies were unbuilt the re-fit was the one reported; building it moved the
    /// refusal to the setting that still has none, which is what a harness needs to hear.
    #[test]
    fn both_outer_settings_at_once_report_the_unbuilt_discovery() {
        let mut both = CallingLoopConfig::default();
        both.slippage_refit.max_rounds = 3;
        both.discovery.mode = DiscoveryMode::AgainstFullConvergence;
        assert_eq!(
            both.validate(),
            Err(CallingLoopConfigError::DiscoveryNotBuilt {
                mode: DiscoveryMode::AgainstFullConvergence
            })
        );
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

    /// **Four things `validate` used to refuse are now values nobody can write**, and that is
    /// a stronger guarantee than a check.
    ///
    /// A pass cap of zero, a discovery round cap of zero, a discovery bar of zero reads and a
    /// share outside `[0, 1]` were four range checks and two error variants. They are now
    /// held by [`NonZeroU32`] and by [`MinAltReads`]'s own constructors, so the illegal value
    /// never exists to be checked — which no mutation of `validate` can undo.
    #[test]
    fn the_counts_and_the_bar_refuse_their_illegal_values_at_construction() {
        assert_eq!(NonZeroU32::new(0), None);
        assert_eq!(MinAltReadShare::new(-0.01), None);
        assert_eq!(MinAltReadShare::new(1.01), None);
        assert_eq!(MinAltReadShare::new(f64::NAN), None);
        assert_eq!(MinAltReadShare::new(f64::INFINITY), None);

        // And the shipped bar answers the question a discovery round asks of one sample:
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

    /// A configuration that is both out of range and not-yet-built is told about the range
    /// first, because that half is the caller's to fix today.
    ///
    /// **Checked across every range refusal, not at one example.** Each of the three is
    /// crossed with the unbuilt setting, so the ordering is pinned as a property rather than
    /// by a single pairing.
    #[test]
    fn an_out_of_range_value_outranks_a_setting_that_is_not_built() {
        /// One way to put a value out of range, named so a failure says which.
        type BreakOneValue = (&'static str, fn(&mut CallingLoopConfig));

        let out_of_range: [BreakOneValue; 3] = [
            ("threshold", |config| {
                config.convergence_threshold = f64::NAN
            }),
            ("pull-back", |config| {
                config
                    .slippage_refit
                    .direction_and_fall_off_pull_back_pseudocounts = -1.0;
            }),
            ("round threshold", |config| {
                config.slippage_refit.round_convergence_threshold = 0.0;
            }),
        ];

        for (which, break_it) in out_of_range {
            let mut config = CallingLoopConfig::default();
            break_it(&mut config);
            // The unbuilt setting asked for at the same time — and re-fit rounds beside it,
            // which are a built setting and must not change what is reported.
            config.slippage_refit.max_rounds = 2;
            config.discovery.mode = DiscoveryMode::AgainstFullConvergence;

            let refusal = config
                .validate()
                .expect_err("this configuration is refused");
            assert!(
                !matches!(refusal, CallingLoopConfigError::DiscoveryNotBuilt { .. }),
                "the {which} problem should outrank the unbuilt setting, got {refusal}"
            );
        }
    }

    /// The two settings an operator reads about in a log line are named in words, not as Rust
    /// identifiers.
    #[test]
    fn a_refusal_reads_as_a_sentence_rather_than_as_an_identifier() {
        let mut discovering = CallingLoopConfig::default();
        discovering.discovery.mode = DiscoveryMode::AgainstFrozenFrequencies;
        let refusal = discovering.validate().expect_err("not built");
        let sentence = refusal.to_string();
        assert!(
            sentence.contains("discovering alleles against frozen allele frequencies"),
            "{sentence}"
        );
        assert!(!sentence.contains("AgainstFrozenFrequencies"), "{sentence}");

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
            let table =
                crate::ng::calling::GenotypeTable::build(parameters.ploidy(), candidates.len());
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
        let calibration = [crate::ng::calling::ReadGroupCalibration::defaulted()];
        let inbreeding = [crate::ng::types::InbreedingF::try_new(0.0).expect("legal")];
        let strata = crate::ng::parameter_estimation::joint::stratum_fits::StratumFits::over(
            &[],
            std::collections::BTreeMap::new(),
        );
        let substitution = std::collections::BTreeMap::new();
        let parameters = FrozenParameters::uncontaminated(
            &calibration,
            &inbreeding,
            crate::ng::calling::genotype_prior::SpectrumSeed::new(
                1.0,
                1e-3,
                crate::ng::calling::genotype_prior::SeedRegime::NeutralShape,
            ),
            &strata,
            &substitution,
            crate::ng::types::Ploidy::try_new(2).expect("a diploid"),
        );
        let per_sample = [crate::ng::calling::GenericLocusSample {
            evidence: crate::ng::calling::GenericSampleEvidence::empty(),
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
