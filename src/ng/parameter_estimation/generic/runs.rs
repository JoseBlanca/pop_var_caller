//! The inbreeding coefficient: a two-state hidden Markov model over windows, deciding
//! which of them lie in a **run of homozygosity** — a stretch of genome where both
//! copies descend from one recent ancestor, and which therefore carries almost no
//! heterozygotes. The rest of the genome carries them at the sample's own rate.
//!
//! The chain walks the windows of a contig deciding which state each is in, and the
//! coefficient is the share of the analysable genome the inside state claims —
//! weighted by reference positions rather than by loci, so a window dense in widened
//! indel loci is not under-weighted (`spec/parameter_prepass_generic.md` §6).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §5.3. The fit itself lands
//! in Milestone E; the types below are Milestone A.

use std::collections::BTreeMap;

use smallvec::SmallVec;

use crate::ng::parameter_estimation::fitting::mixture_weights::{
    GenotypeLikelihoodTable, ln_sum_exp, weighted_log_likelihood,
};
use crate::ng::parameter_estimation::fitting::{FitTermination, NoiseModel};
use crate::ng::parameter_estimation::generic::SampleRates;
use crate::ng::parameter_estimation::generic::accumulators::WindowKey;
use crate::ng::parameter_estimation::generic::histogram::DepthAltHistogram;
use crate::ng::parameter_estimation::generic::noise_model::{
    SampleLibraryNoise, SubstitutionNoiseModel,
};
use crate::ng::parameter_estimation::{Estimate, ParameterEstimationError, Provenance};
use crate::ng::types::{ContigId, InbreedingF, Ploidy};

/// Fewest windows the runs model will accept.
///
/// **Not the same kind of floor as
/// [`MIN_SITES_TO_FIT`](super::MIN_SITES_TO_FIT).** That one guards against a rate
/// fitted from too little; this one guards against an answer that is entirely the
/// estimator's own noise. A genome generated with **no runs at all** returned `F`
/// averaging 0.23 at 1,200 windows, and 0.84 on one seed of eight
/// (`spec/parameter_prepass_generic.md` §6.5; §6.1 carries the resolution figures). A tomato genome is 8,004 windows and a
/// human 31,000, so no real run comes near this — but a development fixture or a
/// region-restricted run does, and the number it would produce looks like any other.
/// Fail rather than emit.
pub const MIN_WINDOWS_TO_FIT_INBREEDING: usize = 3_000;

/// The starting points the runs model climbs from.
///
/// **They must disagree about how far apart the two states are, not only about how much
/// of the genome is inside a run — and that is the whole content of this type.** Starts
/// that differ only in the second are not a spread: they all miss a genome whose real
/// states sit close together, in the same way, and "keep the best-scoring" then has
/// nothing better to pick from. Measured on a genome 26% covered by runs: nine starts
/// spanning the separation return `F` = 0.2634 against a realised 0.2629, where five
/// sharing one separation guess return `F` = 0.0000 — converged, and silent
/// (`research/parameter_estimator_experiments_2026-08-06.md` §3.4).
///
/// The default is three separations crossed with three inside fractions, and it costs
/// seconds on 8,000 windows.
///
/// **A caller who replaces the default must keep the separations distinct.** Nothing
/// here enforces it — the fit is what rejects a start set that leaves one state empty,
/// with `InbreedingStatesNotSeparated` rather than a zero.
#[derive(Clone, PartialEq, Debug)]
pub struct RunsModelStarts {
    /// The inside state's heterozygote rate **as a fraction of the outside state's**.
    ///
    /// **A smaller number means the two states are guessed further apart**, which is the
    /// opposite of what the word "separation" suggests: 0.05 guesses the inside state at
    /// a twentieth of the outside rate, and 0.75 at three quarters of it. Getting this
    /// backwards is how a start set ends up spanning nothing.
    ///
    /// Defaults to `[0.05, 1/3, 0.75]`.
    pub separations: SmallVec<[f64; 3]>,
    /// The fraction of the genome each start begins by assuming lies inside a run.
    /// Defaults to `[0.05, 0.5, 0.75]`.
    pub implied_f: SmallVec<[f64; 3]>,
}

impl Default for RunsModelStarts {
    /// Three separations crossed with three inside fractions — nine starts. The values
    /// are on the two fields; the property that matters is that the *separations*
    /// differ, not that there are nine.
    fn default() -> Self {
        Self {
            separations: SmallVec::from_slice(&[0.05, 1.0 / 3.0, 0.75]),
            implied_f: SmallVec::from_slice(&[0.05, 0.5, 0.75]),
        }
    }
}

/// One starting point's result, kept so that a search which found nothing can be told
/// apart from a genome that has nothing to find.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct StartOutcome {
    pub separation: f64,
    pub implied_f: f64,
    pub inbreeding: f64,
    pub log_likelihood: f64,
}

/// What the runs model fitted alongside `F`, and how it terminated.
///
/// Emitted because every one of these is a number someone will want when `F` looks
/// wrong.
#[derive(Clone, PartialEq, Debug)]
pub struct RunsModelFit {
    /// The heterozygote rate outside a run — the sample's ordinary heterozygosity.
    pub outside_het: f64,
    /// The homozygous-non-reference rate outside a run.
    pub outside_hom_alt: f64,
    /// The heterozygote rate **inside** a run — the small residue of apparent
    /// heterozygotes a truly autozygous stretch still shows, from collapsed paralogs and
    /// mismapping.
    ///
    /// **A fitted rate, and "floor" names its role rather than its arithmetic**: it is
    /// what stops the inside state charging an impossible cost. Fitted rather than fixed
    /// at zero, because at zero one collapsed paralog inside a run costs the whole
    /// heterozygote-against-homozygous-reference ratio at that site — about 125 for one
    /// alternative read of three, and past 10³⁵ for fifteen of thirty.
    pub inside_het_floor: f64,
    /// The homozygous-non-reference rate inside a run.
    pub inside_hom_alt: f64,
    /// The fitted per-base chance of entering a run.
    pub enter_run_per_base: f64,
    /// The fitted per-base chance of leaving one.
    pub leave_run_per_base: f64,
    pub termination: FitTermination,

    /// **How the search went, and the only thing separating a real `F` = 0 from a
    /// failed one.** Both leave the inside state empty and its frequencies at their
    /// starting values; only the scores say whether a better answer was looked for and
    /// rejected. Sorted by score, best first — read the winner through
    /// [`Self::best_start`] rather than indexing, so the ordering is stated in one place.
    pub starts_tried: SmallVec<[StartOutcome; 9]>,

    /// **A reference point, not a threshold**: what `F` came back as on genomes with no runs
    /// at this window count *in the one regime it was measured in* — 100,000 sites a window,
    /// which is about a hundred heterozygotes. About 0.01 at tomato's 8,004 windows and 0.003
    /// at a human genome's 31,000. See [`resolution_at`], which says what it can and cannot
    /// support.
    ///
    /// **It is not a line below which nothing was detected, and it was demoted from claiming
    /// to be one** (owner, 2026-08-09). Measured on a fixture at a fifth of that evidence —
    /// 3,600 windows of 400 sites — three of the five accepted fits on genomes with **no runs
    /// at all** returned `F` = 0.0322, 0.0359 and 0.0505 against a reported resolution of
    /// 0.0291. A consumer treating it as a floor would have called those detections.
    ///
    /// **What does separate a real `F` from nothing found is
    /// [`RunsModelFit::spread_across_tied_starts`]**, which is measured on both populations
    /// and is a property of *this* fit rather than an interpolation from a different regime.
    /// The failure this number was standing in for is now refused outright
    /// (`MAX_IDENTIFIED_START_SPREAD`), so what remains here is a comparison a reader can make
    /// against other runs, not a test a fit has passed.
    pub resolution: f64,

    /// Windows whose posterior landed between 0.01 and 0.99 — the ones the chain rather
    /// than their own reads decided. Non-zero is not a fault; it is the chain earning its
    /// keep.
    ///
    /// **Zero at 100 kb on the research harness's genomes**, which is the measurement saying
    /// the transitions changed no window's classification there. Do not read it as the
    /// expected value on real data: this chain includes every window of a contig, absent ones
    /// as empty, and an empty window emits nothing under either state — so it is decided by
    /// the chain by definition. A genome with unmappable stretches has a non-zero count for
    /// that reason alone.
    pub undecided_windows: u32,
}

/// How many Baum–Welch passes one starting point is allowed.
///
/// The research harness's own cap (`examples/ng_inbreeding_harness.rs`, `fit_all_starts`).
/// Generous: the recovery runs on this file's own fixtures settle in four to six passes.
pub const MAX_RUNS_MODEL_ITERATIONS: u32 = 300;

/// How long a run of homozygosity a starting point assumes, in windows.
///
/// **A start's two transition rates are derived from this and its
/// [`RunsModelStarts::implied_f`]**: `leave = 1 / MEAN_RUN_WINDOWS_AT_START` and
/// `enter = leave · F₀ / (1 − F₀)`, which is the rule the research harness generates its
/// scenarios by (`base_scenario`). Twenty windows is 2 Mb at
/// [`INBREEDING_WINDOW_BP`](super::INBREEDING_WINDOW_BP), which is the **shortest** of the
/// three run lengths the harness's scenarios use — it also uses 30 and 50.
///
/// **It is a start and not an assumption.** Both rates are fitted, so this decides how
/// many passes a start takes and not what it converges to — which is why the spread that
/// matters is over [`RunsModelStarts::separations`] and not over this.
const MEAN_RUN_WINDOWS_AT_START: f64 = 20.0;

/// A window whose posterior lands outside `[UNDECIDED_LOW, UNDECIDED_HIGH]` was settled by
/// its own evidence; one inside was settled by the chain.
const UNDECIDED_LOW: f64 = 0.01;
const UNDECIDED_HIGH: f64 = 0.99;

/// The floor a rate is clamped to, so that a state the fit emptied can still be re-entered
/// and so that `ln 0` never reaches the chain.
const RATE_FLOOR: f64 = 1e-9;

/// How close the two states' heterozygote rates may come before `F` stops being identified
/// at all: the inside rate must be **at most this fraction** of the outside one.
///
/// **The design's identifying constraint is `h << Hout`** (spec §6.1) — much less than, not
/// merely less than — and this is that constraint given a number. Research note §3.1 proves
/// why it is not a matter of precision: at coincident states every window's emission is the
/// same under either state, so the likelihood is *exactly* flat in `F` and what a fit
/// returns is its own accident.
///
/// **Where the number comes from.** Measured across this file's own fixtures, a legitimate
/// fit's inside rate sits at 0.036 to 0.045 of its outside rate with no false heterozygotes,
/// rising to **0.842** under a floor of spurious heterozygotes five times the real rate —
/// which is the most extreme case spec §6.2 asks this estimator to survive. The two
/// confident wrong answers on genomes with no runs at all sat at **0.929 and 0.968**. Nine
/// tenths is between them, nearer the failures than the extreme legitimate case.
///
/// **The asymmetry is what settles it.** Refusing a fit costs an error the caller answers by
/// supplying `F`; accepting one costs `F` = 0.9995 on an outcrossing genome, and a cohort's
/// diversity divides by `1 − F`. Spec §6.1's own instruction for this parameter is *fail
/// rather than emit*.
pub const MAX_IDENTIFIED_STATE_RATIO: f64 = 0.9;

/// How far the **tied** starting points may disagree about `F` before `F` stops being
/// identified: the largest fitted value minus the smallest, over the starts that scored
/// within [`MAX_TIED_START_LOG_LIKELIHOOD_GAP`] of the best.
///
/// **The third refusal, and it catches a failure the other two cannot see.**
/// [`MAX_IDENTIFIED_STATE_RATIO`] refuses a fit whose two states came out too close to tell
/// apart, and the used-both-states check refuses a search that collapsed. Neither sees the
/// failure measured in `research/inbreeding_resolution_2026-08-09.md`: on a genome with **no
/// runs at all**, at five heterozygotes a window, the chain manufactures a separation out of
/// sampling noise and returns `F` = 0.99 with both other checks satisfied. On the harness's
/// two worst cells the fitted state ratio averages 0.31 and 0.62 — nowhere near the 0.9 that
/// would refuse it (§2); the fixture behind this file's own test is further out still, its
/// winning start at a ratio of 0.086. It is not a floor that can be raised either: 3,000
/// windows × 5 heterozygotes fails where 6,000 × 5 does not, so what has to be detected is
/// the failure mode and not a magnitude (§1).
///
/// **The two populations do not overlap, and the reason is not margin.** Where the runs are
/// real, the tied starts do not merely agree closely — they agree **exactly**. Measured at
/// 0.0000 over the 160 fits of §3 and §4 (twelve window-count shapes plus twenty run shapes),
/// and at 3.7 × 10⁻¹¹ on this file's own fixture where the assertion is on a printed value
/// rather than a rounded table. A review swept a selfing landrace to a realised `F` of 0.9719
/// across three evidence levels and false-heterozygote floors of 0 to 5 times the real rate:
/// **every one at 0.0000**.
///
/// **Where a genome has no runs, there is no cluster to sit outside of.** Twenty-four seeds at
/// this file's fixture shape, all drawn with no runs: seven refused for never separating, and
/// the rest spread **0.0213, 0.0264, 0.0323, 0.0337, 0.0483, 0.0521, 0.0541, 0.0584, 0.0621,
/// 0.0675, 0.0722, 0.0740, 0.1027, 0.1144, 0.1368, 0.1375, 0.2749** — a continuum, with no gap
/// at a twentieth or anywhere else. **So this threshold is not separating two clusters. It
/// chooses how aggressively to refuse a population that is unidentified at every point in
/// it**, and the only bound measurement gives is that it must be above zero and below the
/// worst case. At 0.05 it refuses 41 of the 46 no-runs genomes drawn for it and accepts five
/// returning `F` between 0.0019 and 0.0505.
///
/// **A twentieth rather than a thousandth, and the reason is real data.** Nothing measured
/// distinguishes them — both accept every legitimate fit and refuse every catastrophic one —
/// so what decides it is that the legitimate population's *exactly* 0.0000 is a property of
/// genomes drawn from the model the fit assumes, and real data will not be that clean. The
/// cost of the choice runs the other way from spec §6.1's *fail rather than emit*, and it is
/// bounded: what a twentieth lets through instead of a thousandth is an outbred genome
/// answered at an `F` of a few hundredths.
pub const MAX_IDENTIFIED_START_SPREAD: f64 = 0.05;

/// How far behind the best-scoring start another start may sit and still count as **tied**
/// with it, in nats of log-likelihood.
///
/// **Without this, the spread above refuses the robustness spec §6.2 requires.** A genome
/// 30% covered by runs, with a floor of spurious heterozygotes five times the real rate,
/// fits perfectly: six of the nine starts land on `F` = 0.3157 against a realised 0.3122.
/// The other three collapse to `F` = 0.0000 — and they score **1,473 nats worse**. The raw
/// range across all nine is 0.3157, which a threshold of a twentieth refuses, and what it
/// would be refusing is three starts the fit has already rejected by a likelihood ratio of
/// e^1473.
///
/// **The failure looks nothing like that, and the score is what says so.** On a genome with
/// no runs at all the nine starts return seven distinct values — `F` = 0.0010 from three of
/// them, then 0.8497, 0.6015, 0.5722, 0.0003, 0.0051 and 0.0159 — and every one is within
/// **0.91 nats** of the best. They are not competing answers; they are the same answer,
/// because research note §3.1's proof says the likelihood is flat in `F` when the two states
/// coincide. A start beaten by 0.91 nats has been beaten by about five to two, which decides
/// nothing.
///
/// So the two populations are separated by the **score gap** and not by the spread: 0.91
/// nats where the answer is not identified, 1,473 where a start was genuinely rejected — a
/// factor of 1,600. Ten nats sits between them, eleven times above the first and a hundred
/// and forty-seven times below the second.
///
/// **Why an absolute number of nats is the right shape, and it is the opposite of the obvious
/// argument.** A log-likelihood difference between two genuinely different fits *does* grow
/// with the genome — each window contributes, so 1,473 nats on this fixture would be about
/// 2,900 on one twice as long. What does **not** grow is the gap in the failure case: there
/// the two fits are the same distribution, so their difference is sampling noise and stays of
/// order one at any genome size. An absolute cut therefore separates them at every scale, and
/// a cut proportional to the total would drift into the failure population as the genome grew.
/// (The two totals measured here differ eightfold, −225,690 against −1,802,300, and the two
/// gaps differ by 1,600 — so proportionality does not describe them either way.)
pub const MAX_TIED_START_LOG_LIKELIHOOD_GAP: f64 = 10.0;

// **The bounds the two constants above were measured between, checked at build time.** Each
// is a number a later edit could move for a plausible-sounding reason without re-running
// anything, and each has a measured side it must not cross; a compile error says so at the
// moment of the edit, where a red test says so after a build. The same device guards the
// error-rate ladder's rung count (`generic/mod.rs`).
//
// **The spread's two bounds are not alike in strength, and saying so is the point.** The lower
// one is tight: the measured legitimate population sits at 3.7 × 10⁻¹¹, so anything at or
// below that refuses fits that recovered their genome. The upper one is **weak** — it only
// says the criterion still refuses the worst case ever measured. There is no tight upper
// bound to be had, because the no-runs population is a continuum with no gap in it; what
// picks a twentieth out of that interval is the judgement recorded on the constant, not a
// measurement. An earlier version of this assert claimed 0.1147 as "the narrowest
// disagreement measured on a genome with no runs" and was wrong: 0.1147 was the narrowest
// among the fits this very threshold had *refused*, so the sample was censored by the number
// it was being used to justify. Uncensored, the same shape gives 0.0213.
const _: () = assert!(
    MAX_IDENTIFIED_START_SPREAD > 1e-6,
    "the tied starts of a fit that recovered its genome agree to 3.7e-11, and a threshold at \
     or below that refuses every fit this estimator is for"
);
const _: () = assert!(
    MAX_IDENTIFIED_START_SPREAD < 0.8494,
    "0.8494 is the narrowest spread measured on a genome with no runs that the fit would \
     otherwise have answered near one; at or above it this check stops refusing its own worst \
     case"
);
const _: () = assert!(
    MAX_TIED_START_LOG_LIKELIHOOD_GAP > 0.91,
    "starts 0.91 nats apart could not be told apart by the data and must all count as tied, \
     or the failure this exists for is invisible"
);
const _: () = assert!(
    MAX_TIED_START_LOG_LIKELIHOOD_GAP < 1473.0,
    "a start beaten by 1,473 nats has been rejected by the data, and counting it as tied \
     refuses the fit spec §6.2 requires the estimator to survive"
);

/// What `F` came back as on genomes with **no runs at all**, at a given window count, in the
/// **one regime it was measured in**: 100,000 sites a window, about a hundred heterozygotes.
///
/// **A reference point and not a threshold** (owner, 2026-08-09). An earlier version of this
/// sentence said *an `F` below it means nothing detected*, and that is not a property it has —
/// it is a function of the window count alone, and the noise floor is not. See the two
/// paragraphs at the end for what was measured against it, and
/// [`RunsModelFit::spread_across_tied_starts`] for the quantity that *does* separate a real
/// `F` from a search that found nothing.
///
/// **Interpolated from the eight-seed means of research note §3.6 rather than computed
/// from a formula**, because there is no formula: the number is what a two-state model
/// reads out of ordinary sampling wobble in window heterozygote counts, and it was
/// measured at four window counts and nowhere else.
///
/// | windows | mean `F` on a genome with no runs |
/// |---:|---:|
/// | 1,200 | 0.226 |
/// | 4,800 | 0.017 |
/// | 19,200 | 0.0040 |
/// | 76,800 | 0.0016 |
///
/// Log–log linear between those points and flat outside them. That reproduces both figures
/// the design docs quote — **0.009970 at tomato's 8,004 windows** against their "about
/// 0.01", and **0.002914 at a human genome's 31,000** against their "0.003" — which is what
/// makes this an interpolation of the measurement rather than a second, disagreeing
/// statement of it.
///
/// **⚠ It is a function of the window *count* alone, and the noise floor is not.** §3.6's
/// eight seeds were drawn at 100,000 sites a window; how much a two-state model can read out
/// of sampling wobble depends on how much evidence each window carries, and this function
/// cannot see that. **And no second argument would repair it**, which is the finding that
/// settled this: research note §1 measures the floor across three window counts crossed with
/// four evidence levels and it is not smooth in either — 3,000 windows × 5 heterozygotes fails
/// where 6,000 × 5 does not, and 3,000 × 100 behaves perfectly where 12,000 × 5 did not. What
/// varies is not a magnitude but whether a rare catastrophic mode fires.
///
/// **⚠ And it is exceeded in practice, which is why it is not a threshold.** Measured on a
/// fixture at a fifth of §3.6's evidence — 3,600 windows of 400 sites — of the five fits this
/// module accepted on genomes with **no runs at all**, three returned `F` = 0.0322, 0.0359 and
/// 0.0505 against the 0.0291 reported here. A consumer using this as a floor would have called
/// all three detections.
///
/// **What the number is still good for** is comparing one run against another: two samples at
/// the same window count are being read at the same resolution, and a fitted `F` far above it
/// means more than one just above it. What it is not good for is deciding whether a single fit
/// found anything — the fit itself now answers that, by refusing outright when the starting
/// points that scored alike disagree (`MAX_IDENTIFIED_START_SPREAD`), which is the failure
/// this number was standing in for.
///
/// **Kept rather than removed** (owner, 2026-08-09): spec §6.5 names it as one of three things
/// that go out with `F`, nothing reads it yet, and the repair its ⚠ asks for — the floor as a
/// function of evidence per window — is a measurement §1 shows cannot exist in that form.
#[must_use]
pub fn resolution_at(windows: usize) -> f64 {
    /// The measured points, ascending in window count.
    const MEASURED: [(f64, f64); 4] = [
        (1_200.0, 0.226),
        (4_800.0, 0.017),
        (19_200.0, 0.0040),
        (76_800.0, 0.0016),
    ];

    let windows = (windows as f64).max(1.0);
    if windows <= MEASURED[0].0 {
        return MEASURED[0].1;
    }
    for pair in MEASURED.windows(2) {
        let ((low_windows, low_f), (high_windows, high_f)) = (pair[0], pair[1]);
        if windows <= high_windows {
            let along = (windows.ln() - low_windows.ln()) / (high_windows.ln() - low_windows.ln());
            return (low_f.ln() + along * (high_f.ln() - low_f.ln())).exp();
        }
    }
    MEASURED[MEASURED.len() - 1].1
}

impl RunsModelFit {
    /// The starting point that scored best — the one `F` was taken from.
    ///
    /// The single reader of `starts_tried`'s "best first" ordering, so that the ordering
    /// is a property of this file rather than of every call site that indexes. `None`
    /// only if no start was tried at all, which the fit does not produce.
    #[must_use]
    pub fn best_start(&self) -> Option<&StartOutcome> {
        self.starts_tried.first()
    }

    /// **How far the tied starting points disagreed about `F`** — the criterion
    /// [`MAX_IDENTIFIED_START_SPREAD`] is applied to, reported so that a consumer can see
    /// how much margin an accepted fit had.
    ///
    /// A method rather than a field, so the number a caller reads and the number the fit was
    /// accepted on cannot drift apart.
    #[must_use]
    pub fn spread_across_tied_starts(&self) -> f64 {
        spread_across_tied_starts(&self.starts_tried)
    }

    /// The starting points the spread above was taken over: those within
    /// [`MAX_TIED_START_LOG_LIKELIHOOD_GAP`] of the best.
    ///
    /// **Exposed because a spread of zero means two different things and this is what tells
    /// them apart.** Nine starts agreeing is the margin the threshold rests on; *one* start
    /// compared with itself is a criterion that did nothing, and both report 0.0000. It
    /// happens: at a realised `F` of 0.95 with a floor of spurious heterozygotes five times
    /// the real rate, one start in nine is tied. Without this, no test and no consumer can
    /// see the difference.
    #[must_use]
    pub fn tied_starts(&self) -> &[StartOutcome] {
        tied_starts(&self.starts_tried)
    }
}

/// The starts that scored within [`MAX_TIED_START_LOG_LIKELIHOOD_GAP`] of the best.
///
/// **A prefix, because `starts_tried` is best-first** — the ordering is the type's
/// documented one, so the tied starts are its leading run and this cannot silently pick up
/// a stray start from further down.
///
/// A start beaten by more than the gap is not a competing answer, it is a rejected one, and
/// including it in the spread is what would refuse a fit that recovered its genome. A start
/// whose likelihood is `NaN` is treated as tied, so a fit that produced one is refused
/// rather than waved through — which is the direction spec §6.1 asks for.
fn tied_starts(starts: &[StartOutcome]) -> &[StartOutcome] {
    let Some(best) = starts.first() else {
        return starts;
    };
    let beaten = starts.iter().position(|start| {
        best.log_likelihood - start.log_likelihood > MAX_TIED_START_LOG_LIKELIHOOD_GAP
    });
    &starts[..beaten.unwrap_or(starts.len())]
}

/// The largest fitted `F` minus the smallest, over the tied starts.
///
/// Zero for an empty slice, which the fit does not produce — it refuses an empty start set
/// before this is reached.
///
/// **A `NaN` among the fitted values makes the whole spread `NaN`**, which the caller treats
/// as a refusal. It has to be written out rather than left to `f64::max`, which *ignores*
/// `NaN` and would have dropped such a start from the range silently — and if every tied
/// start were `NaN` the fold would return `−∞`, compare false against the threshold, and
/// **accept**. That is the exact inversion of [`tied_starts`]' rule for an unreadable score,
/// in the one place a fit is waved through rather than refused.
fn spread_across_tied_starts(starts: &[StartOutcome]) -> f64 {
    let tied = tied_starts(starts);
    if tied.is_empty() {
        return 0.0;
    }
    let mut highest = f64::NEG_INFINITY;
    let mut lowest = f64::INFINITY;
    for start in tied {
        if start.inbreeding.is_nan() {
            return f64::NAN;
        }
        highest = highest.max(start.inbreeding);
        lowest = lowest.min(start.inbreeding);
    }
    highest - lowest
}

/// Fit the fraction of the analysable genome lying in runs of homozygosity.
///
/// **Diploid only.** Above two copies `F` needs several identity-by-descent coefficients
/// and is deferred (spec §7), and below two there are no heterozygotes to be short of.
///
/// `windowed_histograms` is the accumulator's windowed table. Only its
/// [`WindowKey::InWindow`] entries at `ploidy` take part: a run handed `F` keeps one table
/// for the whole sample and has no chain to walk.
///
/// `noise` carries **every library's share and its own error rate**, not one pooled rate.
/// A site with few alternative reads keeps which library each came from, and each must be
/// weighed against its own library's rate; where a site pooled, the same expression reduces
/// to the share-weighted mean, so no second rule is needed (arch §5.3).
///
/// `outside_rates` is where each start's outside state begins — the sample's own genotype
/// frequencies, as the coupled fit left them. A start's inside state is that state's
/// heterozygote rate multiplied by one of [`RunsModelStarts::separations`].
///
/// # Errors
///
/// [`ParameterEstimationError::InbreedingNotFittable`] below
/// [`MIN_WINDOWS_TO_FIT_INBREEDING`] **windows that hold sites** — the estimator's own
/// noise is the size of the answer there, and the number it would produce looks like any
/// other.
///
/// [`ParameterEstimationError::InbreedingStatesNotSeparated`] when no start left posterior
/// mass on both states, or when every start's two states came out within
/// [`MAX_IDENTIFIED_STATE_RATIO`] of each other. **Never a zero**: a failed search and an
/// outcrossing genome leave identical fitted values, and only this error tells them apart.
///
/// [`ParameterEstimationError::InbreedingStartsDisagree`] when the starts landed further
/// than [`MAX_IDENTIFIED_START_SPREAD`] apart. **Three refusals and three different
/// failures**: a search that collapsed, two states that coincide, and a chain that found a
/// different answer from every start — which is what reading sampling noise looks like, and
/// what the first two are blind to.
///
/// # Panics
///
/// If `ploidy` is not two, if `starts` is empty on either axis, or if a start's separation
/// or implied fraction is not strictly inside `(0, 1)`.
pub fn fit_inbreeding(
    sample: &str,
    windowed_histograms: &BTreeMap<(WindowKey, Ploidy), DepthAltHistogram<u32>>,
    noise: &SampleLibraryNoise,
    ploidy: Ploidy,
    outside_rates: &SampleRates,
    starts: &RunsModelStarts,
) -> Result<(Estimate<InbreedingF>, RunsModelFit), ParameterEstimationError> {
    assert_eq!(
        ploidy.get(),
        2,
        "the inbreeding coefficient is diploid-only: above two copies it needs several \
         identity-by-descent coefficients and below two there are no heterozygotes"
    );
    assert!(
        !starts.separations.is_empty() && !starts.implied_f.is_empty(),
        "a runs model with no starting points cannot report what it looked for"
    );

    let chain = WindowChain::over(windowed_histograms, ploidy, noise);
    if chain.windows_holding_sites < MIN_WINDOWS_TO_FIT_INBREEDING {
        return Err(ParameterEstimationError::InbreedingNotFittable {
            sample: sample.to_string(),
            windows: chain.windows_holding_sites,
            floor: MIN_WINDOWS_TO_FIT_INBREEDING,
        });
    }

    let outside = [
        outside_rates.by_alt_copies()[0].get(),
        outside_rates.by_alt_copies()[1].get(),
        outside_rates.by_alt_copies()[2].get(),
    ];

    let mut outcomes: Vec<(RunsModelStart, StateFit)> = Vec::new();
    for &separation in &starts.separations {
        for &implied_f in &starts.implied_f {
            let start = RunsModelStart::new(separation, implied_f, outside);
            let fitted = chain.fit_from(&start, MAX_RUNS_MODEL_ITERATIONS);
            outcomes.push((start, fitted));
        }
    }
    // Best first, which is the ordering `starts_tried` documents and `best_start` reads.
    outcomes.sort_by(|left, right| {
        right
            .1
            .log_likelihood
            .total_cmp(&left.1.log_likelihood)
            .then(left.0.separation.total_cmp(&right.0.separation))
    });

    let starts_tried: SmallVec<[StartOutcome; 9]> = outcomes
        .iter()
        .map(|(start, fitted)| StartOutcome {
            separation: start.separation,
            implied_f: start.implied_f,
            inbreeding: fitted.inbreeding,
            log_likelihood: fitted.log_likelihood,
        })
        .collect();

    if !outcomes.iter().any(|(_, fitted)| fitted.separated_states) {
        return Err(ParameterEstimationError::InbreedingStatesNotSeparated {
            sample: sample.to_string(),
            starts: starts_tried.len(),
        });
    }

    // **The third refusal, and the one the other two are blind to.** A genome with no runs
    // at all can have its chain manufacture a separation out of sampling noise: both states
    // are used, the two of them come out far enough apart to clear
    // `MAX_IDENTIFIED_STATE_RATIO` with room to spare, and `F` comes back near one. What
    // gives it away is that the starts which cannot be told apart by score then land in
    // completely different places — 0.0000 across 160 fits wherever the runs are real, and up
    // to 0.9980 wherever they are not
    // (`research/inbreeding_resolution_2026-08-09.md` §4, §6).
    //
    // **Checked after the separation checks and not before**, so that a collapsed search —
    // whose starts also disagree — is reported as the collapsed search it is rather than as
    // a disagreement about a value none of them found.
    let tied = tied_starts(&starts_tried).len();
    let spread = spread_across_tied_starts(&starts_tried);
    // `is_nan()` first and not folded into the comparison: every comparison against `NaN` is
    // false, so `spread > threshold` alone would **accept** a fit whose starts disagreed by an
    // unreadable amount.
    if spread.is_nan() || spread > MAX_IDENTIFIED_START_SPREAD {
        return Err(ParameterEstimationError::InbreedingStartsDisagree {
            sample: sample.to_string(),
            tied_starts: tied,
            starts: starts_tried.len(),
            spread,
            threshold: MAX_IDENTIFIED_START_SPREAD,
        });
    }

    let best =
        best_separated_start(&outcomes).expect("some start separated the states, checked above");

    let inbreeding = InbreedingF::try_new(best.inbreeding)
        .expect("a coverage-weighted posterior occupancy is a fraction");
    let fit = RunsModelFit {
        outside_het: best.frequencies[OUTSIDE][1],
        outside_hom_alt: best.frequencies[OUTSIDE][2],
        inside_het_floor: best.frequencies[INSIDE][1],
        inside_hom_alt: best.frequencies[INSIDE][2],
        enter_run_per_base: per_base(best.enter_run),
        leave_run_per_base: per_base(best.leave_run),
        termination: best.termination,
        starts_tried,
        // **The windows that hold sites, not the chain's length.** The chain pads every
        // contig from window zero, so a region-restricted run's chain is far longer than
        // its evidence: 3,601 windows of sites plus one contig whose only site sits at
        // window 8,000 gives a chain of 11,601, and reading the resolution off that
        // reports 0.006768 where the truth is 0.029067 — understated 4.3-fold, on the
        // one number a consumer uses to decide whether a small `F` means anything.
        resolution: resolution_at(chain.windows_holding_sites),
        undecided_windows: best.undecided_windows,
    };

    Ok((
        Estimate {
            value: inbreeding,
            provenance: Provenance::FittedHere,
            observations: chain.covered_positions_total(),
        },
        fit,
    ))
}

/// A per-**window** transition probability as a per-**base** one.
///
/// The chain steps window to window, so what Baum–Welch fits is the chance of switching
/// state across one window of [`INBREEDING_WINDOW_BP`](super::INBREEDING_WINDOW_BP) bases.
/// This is `1 − e^{−rL} = p` inverted — the *continuous* switching rate, not the discrete
/// `1 − (1 − r)^L = p`, which an earlier version of this comment claimed. The two agree to
/// about seven significant figures at the rates involved, so the number is the same either
/// way and only the identity was misstated; it is `p / L` to first order.
///
/// **Reported and not used.** Nothing in the fit reads a per-base rate; it is on
/// [`RunsModelFit`] because a caller comparing two samples' runs wants a number that does
/// not move with the window size.
fn per_base(per_window: f64) -> f64 {
    let length = super::INBREEDING_WINDOW_BP.get() as f64;
    -(1.0 - per_window.min(1.0 - RATE_FLOOR)).ln() / length
}

/// **The ordering constraint, applied by relabelling after the fit.** The state with the
/// **lower** heterozygote rate is the one inside a run; nothing in Baum–Welch says so, and
/// a fit that landed the other way round would report the sample's ordinary heterozygosity
/// as its inside floor and `1 − F` as `F`.
///
/// Everything indexed by state moves together — the frequencies, the two transition rates
/// and every window's posterior — which is why this is one function and not three swaps at
/// the call site. Returns whether it fired.
///
/// **Its own function because otherwise it is unreachable.** Every start builds its inside
/// state at a *fraction* of the outside heterozygote rate, so the fit begins on the right
/// side of the constraint and no drawn genome in the tests crosses over; deleting the
/// relabelling left every one of them green. What it guards against is a fit that does
/// cross, and this is where that is stated and checked.
fn relabel_if_inverted(
    frequencies: &mut [[f64; GENOTYPES]; STATES],
    enter_run: &mut f64,
    leave_run: &mut f64,
    posteriors: &mut [[f64; STATES]],
) -> bool {
    if frequencies[INSIDE][1] <= frequencies[OUTSIDE][1] {
        return false;
    }
    frequencies.swap(INSIDE, OUTSIDE);
    std::mem::swap(enter_run, leave_run);
    for window in posteriors.iter_mut() {
        window.swap(INSIDE, OUTSIDE);
    }
    true
}

/// The best-scoring start **among those that left posterior mass on both states**, and not
/// simply the best-scoring start.
///
/// A collapsed fit can outscore a real one — its inside state is free to take whatever
/// frequencies fit the whole genome — so picking the best of everything would return the
/// zero [`ParameterEstimationError::InbreedingStatesNotSeparated`] exists to refuse.
///
/// `outcomes` is best-scoring first, which is the ordering [`RunsModelFit::starts_tried`]
/// documents; this takes the first that separated.
///
/// **Its own function because otherwise it is untestable.** On every drawn genome in the
/// tests all nine starts land on the same answer with the same score, so "keep the
/// best-scoring" survives being replaced by "keep the first" — the rule can be checked
/// here, over outcomes written down, and nowhere else.
fn best_separated_start(outcomes: &[(RunsModelStart, StateFit)]) -> Option<&StateFit> {
    outcomes
        .iter()
        .find(|(_, fitted)| fitted.separated_states)
        .map(|(_, fitted)| fitted)
}

/// The state a window is in. `INSIDE` is a run of homozygosity.
const INSIDE: usize = 0;
const OUTSIDE: usize = 1;
const STATES: usize = 2;
/// Three genotypes: homozygous reference, heterozygous, homozygous non-reference.
const GENOTYPES: usize = 3;

/// One starting point, expanded from a `(separation, implied_f)` pair.
#[derive(Copy, Clone, Debug)]
struct RunsModelStart {
    separation: f64,
    implied_f: f64,
    /// `[state][genotype]`.
    frequencies: [[f64; GENOTYPES]; STATES],
    enter_run: f64,
    leave_run: f64,
}

impl RunsModelStart {
    /// # Panics
    ///
    /// If `separation` or `implied_f` is not strictly inside `(0, 1)`. A separation of one
    /// makes the two states identical, at which point `F` is not identified at all — the
    /// likelihood is exactly flat in it — and a separation of zero fixes the inside
    /// heterozygote rate at the value the fit exists to estimate.
    fn new(separation: f64, implied_f: f64, outside: [f64; GENOTYPES]) -> Self {
        assert!(
            separation > 0.0 && separation < 1.0,
            "a separation of {separation} is not a fraction of the outside heterozygote \
             rate strictly between nothing and all of it"
        );
        assert!(
            implied_f > 0.0 && implied_f < 1.0,
            "a start that assumes {implied_f} of the genome is inside a run assumes the \
             answer"
        );

        let inside_het = outside[1] * separation;
        let inside = [1.0 - inside_het - outside[2], inside_het, outside[2]];
        let leave_run = 1.0 / MEAN_RUN_WINDOWS_AT_START;
        let enter_run = leave_run * implied_f / (1.0 - implied_f);

        Self {
            separation,
            implied_f,
            frequencies: [inside, outside],
            enter_run: enter_run.clamp(RATE_FLOOR, 1.0 - RATE_FLOOR),
            leave_run: leave_run.clamp(RATE_FLOOR, 1.0 - RATE_FLOOR),
        }
    }
}

/// What one starting point's Baum–Welch arrived at.
struct StateFit {
    /// `[state][genotype]`, after the ordering constraint.
    frequencies: [[f64; GENOTYPES]; STATES],
    enter_run: f64,
    leave_run: f64,
    /// The coverage-weighted posterior occupancy of the inside state — **this** is `F`.
    inbreeding: f64,
    log_likelihood: f64,
    termination: FitTermination,
    undecided_windows: u32,
    /// Whether some window's posterior favoured each state. **The one thing separating a
    /// genuine `F` near zero from a search that never found the second state**: both leave
    /// the inside state nearly empty, and only this says whether any window was ever
    /// assigned to it.
    separated_states: bool,
}

/// The windows of a sample in genome order, **including the ones that received no site**,
/// each with the likelihood table its cells make.
///
/// **Absent windows are in the chain and that is not a detail.** The accumulator holds only
/// windows a locus reached, so walking it directly would step across an unmappable megabase
/// in one transition and run the end of one chromosome into the start of the next. An empty
/// window emits nothing — its log-emission is zero under either state — but it still costs
/// a transition, which is exactly what a gap in the evidence should cost.
struct WindowChain {
    /// One entry per window of the chain, in genome order.
    windows: Vec<ChainWindow>,
    /// Where each contig's stretch of `windows` begins. The chain restarts at every one.
    contig_starts: Vec<usize>,
    /// How many windows actually received a site — what
    /// [`MIN_WINDOWS_TO_FIT_INBREEDING`] is checked against, because a chain padded with
    /// empty windows has no more evidence in it than the sites it holds.
    windows_holding_sites: usize,
}

/// One window's evidence: how likely each genotype makes each of its cells, how many sites
/// landed in each, and how many reference positions the window covered.
struct ChainWindow {
    /// Row-major, `cells × GENOTYPES`, natural logs. Empty for a window with no sites.
    ln_likelihood_row_major: Vec<f64>,
    /// One site count per cell, in the same order.
    cell_weights: Vec<f64>,
    /// **Reference positions, not loci** — what `F` is weighted by, so that a window dense
    /// in widened indel loci is not under-weighted (spec §6.5).
    covered_positions: f64,
}

impl WindowChain {
    /// Build the chain from the accumulator's windowed table.
    ///
    /// Each contig runs from **window 0** to the last window that received a site, so a
    /// leading gap is as present as an interior one.
    fn over(
        windowed_histograms: &BTreeMap<(WindowKey, Ploidy), DepthAltHistogram<u32>>,
        ploidy: Ploidy,
        noise: &SampleLibraryNoise,
    ) -> Self {
        let model = SubstitutionNoiseModel;
        let mut by_contig: BTreeMap<ContigId, BTreeMap<u32, &DepthAltHistogram<u32>>> =
            BTreeMap::new();
        for (&(key, at), table) in windowed_histograms {
            if at != ploidy {
                continue;
            }
            if let WindowKey::InWindow { contig, window } = key {
                by_contig
                    .entry(contig)
                    .or_default()
                    .insert(window.get(), table);
            }
        }

        let mut windows = Vec::new();
        let mut contig_starts = Vec::new();
        let mut windows_holding_sites = 0;

        for tables in by_contig.values() {
            let last = *tables.keys().next_back().expect("a contig with no windows");
            contig_starts.push(windows.len());
            for index in 0..=last {
                let Some(table) = tables.get(&index) else {
                    windows.push(ChainWindow {
                        ln_likelihood_row_major: Vec::new(),
                        cell_weights: Vec::new(),
                        covered_positions: 0.0,
                    });
                    continue;
                };
                let cells = table.cells(ploidy);
                let mut ln_likelihood_row_major =
                    Vec::with_capacity(cells.len() * model.genotypes(ploidy));
                let mut cell_weights = Vec::with_capacity(cells.len());
                for cell in &cells {
                    model.append_genotype_likelihoods(
                        cell,
                        noise,
                        ploidy,
                        &mut ln_likelihood_row_major,
                    );
                    cell_weights.push(cell.sites as f64);
                }
                if !cell_weights.is_empty() {
                    windows_holding_sites += 1;
                }
                windows.push(ChainWindow {
                    ln_likelihood_row_major,
                    cell_weights,
                    covered_positions: table.total_covered_positions() as f64,
                });
            }
        }

        Self {
            windows,
            contig_starts,
            windows_holding_sites,
        }
    }

    /// How long the chain is, padding included.
    ///
    /// **Test-only, and that is the finding rather than the tidy-up.** It was what the
    /// resolution was read from until the review showed the padding inflates it; the
    /// quantity the fit uses is [`WindowChain::windows_holding_sites`], and nothing outside
    /// the tests now asks how long the chain is.
    #[cfg(test)]
    fn windows(&self) -> usize {
        self.windows.len()
    }

    fn covered_positions_total(&self) -> u64 {
        self.windows
            .iter()
            .map(|window| window.covered_positions as u64)
            .sum()
    }

    /// Whether the window at `position` opens a contig — where the chain restarts.
    fn opens_a_contig(&self, position: usize) -> bool {
        self.contig_starts.binary_search(&position).is_ok()
    }

    /// `ln P(window | state)` — **a sum over the window's cells**, never a per-window
    /// heterozygote count. Dividing by a site count would make a thinly-covered window look
    /// autozygous (spec §6.1).
    ///
    /// A window with no sites scores zero under either state, which is what "no evidence"
    /// has to be for the forward recursion to leave the chain's own transitions in charge
    /// of it.
    fn ln_emission(
        &self,
        position: usize,
        frequencies: &[f64; GENOTYPES],
        scratch: &mut [f64; GENOTYPES],
    ) -> f64 {
        let window = &self.windows[position];
        if window.cell_weights.is_empty() {
            return 0.0;
        }
        weighted_log_likelihood(
            GenotypeLikelihoodTable::from_natural_logs(&window.ln_likelihood_row_major, GENOTYPES),
            &window.cell_weights,
            frequencies,
            scratch,
        )
    }

    /// One starting point's Baum–Welch: forward–backward over the chain, then re-estimate
    /// each state's three genotype frequencies and both transition rates, until the
    /// log-likelihood and the rates both stop moving.
    ///
    /// **It climbs to a stationary point on a surface that is not concave**, which is why
    /// the caller runs it from several starts and why the start is part of the answer. The
    /// ordering constraint — the state with the *lower* heterozygote rate is the one inside
    /// a run — is applied by relabelling at the end, because nothing in Baum–Welch imposes
    /// it.
    fn fit_from(&self, start: &RunsModelStart, max_iterations: u32) -> StateFit {
        let n = self.windows.len();
        let mut frequencies = start.frequencies;
        let mut enter = start.enter_run;
        let mut leave = start.leave_run;

        let mut gamma = vec![[0.0f64; STATES]; n];
        let mut alpha = vec![[f64::NEG_INFINITY; STATES]; n];
        let mut beta = vec![[0.0f64; STATES]; n];
        let mut ln_emission = vec![[0.0f64; STATES]; n];
        let mut scratch = [0.0f64; GENOTYPES];

        let mut log_likelihood = f64::NEG_INFINITY;
        let mut previous = f64::NEG_INFINITY;
        let mut iterations = 0;
        let mut converged = false;

        while iterations < max_iterations {
            iterations += 1;

            let stationary = enter / (enter + leave);
            let ln_initial = [
                stationary.max(RATE_FLOOR).ln(),
                (1.0 - stationary).max(RATE_FLOOR).ln(),
            ];
            // `ln_transition[from][to]`.
            let ln_transition = [
                [
                    (1.0 - leave).max(RATE_FLOOR).ln(),
                    leave.max(RATE_FLOOR).ln(),
                ],
                [
                    enter.max(RATE_FLOOR).ln(),
                    (1.0 - enter).max(RATE_FLOOR).ln(),
                ],
            ];

            for (position, emissions) in ln_emission.iter_mut().enumerate() {
                for (state, slot) in emissions.iter_mut().enumerate() {
                    *slot = self.ln_emission(position, &frequencies[state], &mut scratch);
                }
            }

            // Forward, in log space, restarting at every contig boundary.
            for position in 0..n {
                let opens = self.opens_a_contig(position);
                for state in 0..STATES {
                    alpha[position][state] = if opens {
                        ln_initial[state] + ln_emission[position][state]
                    } else {
                        let reached: [f64; STATES] = std::array::from_fn(|from| {
                            alpha[position - 1][from] + ln_transition[from][state]
                        });
                        ln_sum_exp(&reached) + ln_emission[position][state]
                    };
                }
            }

            // Backward.
            for position in (0..n).rev() {
                let closes = position + 1 == n || self.opens_a_contig(position + 1);
                if closes {
                    beta[position] = [0.0; STATES];
                } else {
                    for state in 0..STATES {
                        let onward: [f64; STATES] = std::array::from_fn(|to| {
                            ln_transition[state][to]
                                + ln_emission[position + 1][to]
                                + beta[position + 1][to]
                        });
                        beta[position][state] = ln_sum_exp(&onward);
                    }
                }
            }

            // Posteriors, and the total log-likelihood as the sum over contigs — one term
            // per contig, taken where the chain ends, because each contig is its own chain.
            log_likelihood = 0.0;
            for position in 0..n {
                let norm = ln_sum_exp(&[
                    alpha[position][INSIDE] + beta[position][INSIDE],
                    alpha[position][OUTSIDE] + beta[position][OUTSIDE],
                ]);
                for state in 0..STATES {
                    gamma[position][state] =
                        (alpha[position][state] + beta[position][state] - norm).exp();
                }
                if position + 1 == n || self.opens_a_contig(position + 1) {
                    log_likelihood += norm;
                }
            }

            // Expected transition counts.
            let mut transitions = [[0.0f64; STATES]; STATES];
            for position in 0..n.saturating_sub(1) {
                if self.opens_a_contig(position + 1) {
                    continue;
                }
                let norm = ln_sum_exp(&[
                    alpha[position][INSIDE] + beta[position][INSIDE],
                    alpha[position][OUTSIDE] + beta[position][OUTSIDE],
                ]);
                for from in 0..STATES {
                    for to in 0..STATES {
                        transitions[from][to] += (alpha[position][from]
                            + ln_transition[from][to]
                            + ln_emission[position + 1][to]
                            + beta[position + 1][to]
                            - norm)
                            .exp();
                    }
                }
            }

            let next_frequencies = self.reestimate_frequencies(&gamma, &frequencies);
            let leaving = transitions[INSIDE][INSIDE] + transitions[INSIDE][OUTSIDE];
            let entering = transitions[OUTSIDE][INSIDE] + transitions[OUTSIDE][OUTSIDE];
            let next_leave = if leaving > 0.0 {
                (transitions[INSIDE][OUTSIDE] / leaving).clamp(RATE_FLOOR, 1.0 - RATE_FLOOR)
            } else {
                leave
            };
            let next_enter = if entering > 0.0 {
                (transitions[OUTSIDE][INSIDE] / entering).clamp(RATE_FLOOR, 1.0 - RATE_FLOOR)
            } else {
                enter
            };

            let moved = (next_enter / enter)
                .ln()
                .abs()
                .max((next_leave / leave).ln().abs());
            frequencies = next_frequencies;
            enter = next_enter;
            leave = next_leave;

            // **A relative tolerance on this log-likelihood is a trap, and it bit the
            // research harness.** The total scales with sites × depth — it is −1.5 × 10⁹ at
            // three reads a site and −2.0 × 10¹⁰ at sixty — so a tolerance of 10⁻⁷ of the
            // total is 151 nats at three reads and 2,013 at sixty. Expectation-maximization
            // stopped after four passes with the answer still moving by hundreds of nats,
            // reported `converged`, and returned `F` = 1.0000 on a genome 32% covered by
            // runs. The tolerance is relative to the **per-window** log-likelihood, which
            // is the quantity that does not grow with the genome, and floored so that it
            // cannot ask for more than floating point can deliver.
            let per_window = log_likelihood.abs() / (n as f64).max(1.0);
            let tolerance = (1e-9 * per_window).max(1e-6) * n as f64;
            if (log_likelihood - previous).abs() < tolerance && moved < 1e-6 {
                previous = log_likelihood;
                converged = true;
                break;
            }
            previous = log_likelihood;
        }
        let log_likelihood = if converged { previous } else { log_likelihood };

        relabel_if_inverted(&mut frequencies, &mut enter, &mut leave, &mut gamma);

        let covered: f64 = self.windows.iter().map(|w| w.covered_positions).sum();
        let inbreeding = if covered > 0.0 {
            gamma
                .iter()
                .zip(&self.windows)
                .map(|(posterior, window)| posterior[INSIDE] * window.covered_positions)
                .sum::<f64>()
                / covered
        } else {
            0.0
        };

        let undecided_windows = gamma
            .iter()
            .filter(|posterior| {
                posterior[INSIDE] > UNDECIDED_LOW && posterior[INSIDE] < UNDECIDED_HIGH
            })
            .count() as u32;

        // **Two conditions, and each catches a failure the other cannot see.**
        //
        // *Both states were used.* Read as "some window's posterior favoured each" rather
        // than as "`F` is above some floor": a genuinely outbred genome still hard-assigns
        // a percent or so of its windows to the inside state — that is the resolution floor
        // of research note §3.6 — while a fit whose **search** collapsed assigns none at
        // all and returns `F` at the rate clamp, with the inside heterozygote rate left at
        // exactly its starting value (§3.4's fingerprint).
        //
        // *And the two states are actually different.* This is the other failure, and it
        // arrives from the opposite side. Research note §3.1 is a **proof** rather than a
        // measurement: if both states carry the same genotype frequencies then every
        // window's emission is the same under either, the observed sequence is independent
        // draws from one distribution, and `P(data | θ)` does not contain the transition
        // rates at all — the likelihood is exactly flat in `F`. What a fit returns there is
        // an accident. Measured on eight genomes drawn with **no runs at all**: six behave
        // (five return `F` between 0.0013 and 0.0100, below their 0.029 resolution, and one
        // fails the used-both-states check), and **two return `F` = 0.9995 and 0.8475** —
        // confident, converged, and wrong by the whole range of the parameter. Those two are
        // exactly the fits whose states coincided, at 0.97 and 0.93 of each other, where
        // every legitimate fit measured here sits at 0.84 or below.
        //
        // The used-both-states check cannot see it: on coincident states the assignment
        // happily uses both, because they are the same distribution — it is splitting noise.
        // Only the *states* say so, and the design's own identifying constraint is
        // `h << Hout` (spec §6.1), which is strictly stronger than the `h <= Hout` the
        // relabelling imposes.
        let used_both_states = gamma.iter().any(|posterior| posterior[INSIDE] > 0.5)
            && gamma.iter().any(|posterior| posterior[INSIDE] < 0.5);
        let states_differ =
            frequencies[INSIDE][1] <= MAX_IDENTIFIED_STATE_RATIO * frequencies[OUTSIDE][1];
        let separated_states = used_both_states && states_differ;

        StateFit {
            frequencies,
            enter_run: enter,
            leave_run: leave,
            inbreeding: inbreeding.clamp(0.0, 1.0),
            log_likelihood,
            termination: FitTermination {
                iterations,
                converged,
            },
            undecided_windows,
            separated_states,
        }
    }

    /// The maximization step over the genotype frequencies: each state's expected genotype
    /// counts under its own window posteriors, normalised.
    ///
    /// **Each state's three frequencies are fitted freely**, not tied through one allele
    /// frequency — that is `bcftools roh`'s form, and with one genome-wide allele frequency
    /// it forces `F` = 0.57 on an outbred genome (spec §6.1). What identifies the states is
    /// the ordering constraint, applied after the fit.
    ///
    /// **Not [`fit_mixture_weights`](super::super::fitting::mixture_weights::fit_mixture_weights).**
    /// That fits one free point on the simplex from one table; here there are two points,
    /// each weighted by a window posterior that the chain itself produced, so the concavity
    /// D1 relies on does not apply (arch §4.1).
    fn reestimate_frequencies(
        &self,
        gamma: &[[f64; STATES]],
        frequencies: &[[f64; GENOTYPES]; STATES],
    ) -> [[f64; GENOTYPES]; STATES] {
        let mut counts = [[0.0f64; GENOTYPES]; STATES];
        for (window, posterior) in self.windows.iter().zip(gamma) {
            for state in 0..STATES {
                let weight = posterior[state];
                if weight < 1e-12 {
                    continue;
                }
                for (cell, &sites) in window.cell_weights.iter().enumerate() {
                    if sites == 0.0 {
                        continue;
                    }
                    let row =
                        &window.ln_likelihood_row_major[cell * GENOTYPES..(cell + 1) * GENOTYPES];
                    // In log space, because a deep cell's likelihood underflows: a site of
                    // depth 124 showing every read alternative is `−857` under the
                    // homozygous-reference genotype at an error rate of 0.001, which is
                    // zero in linear space beside a neighbour at 0.883.
                    let ln_joint: [f64; GENOTYPES] = std::array::from_fn(|genotype| {
                        frequencies[state][genotype].max(f64::MIN_POSITIVE).ln() + row[genotype]
                    });
                    let ln_total = ln_sum_exp(&ln_joint);
                    if !ln_total.is_finite() {
                        continue;
                    }
                    for genotype in 0..GENOTYPES {
                        counts[state][genotype] +=
                            weight * sites * (ln_joint[genotype] - ln_total).exp();
                    }
                }
            }
        }

        let mut next = *frequencies;
        for state in 0..STATES {
            let total: f64 = counts[state].iter().sum();
            if total > 0.0 {
                for genotype in 0..GENOTYPES {
                    // Floored, so that a genotype no window supports cannot pin itself at
                    // zero for the rest of the fit — the expectation step multiplies by it.
                    next[state][genotype] = (counts[state][genotype] / total).max(1e-12);
                }
            }
        }
        next
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The default starts span the state separation, and that is the property the
    /// type exists for.** Three distinct separations is what makes them a spread; three
    /// distinct inside fractions sharing one separation would not be, and the
    /// difference is `F` = 0.2634 against a confident `F` = 0.0000 on the same genome.
    #[test]
    fn the_default_starts_disagree_about_the_state_separation_not_only_about_f() {
        let starts = RunsModelStarts::default();

        assert_eq!(starts.separations.len(), 3);
        assert_eq!(starts.implied_f.len(), 3);
        for pair in starts.separations.windows(2) {
            assert!(
                pair[0] < pair[1],
                "the separations must differ, and ascend: {:?}",
                starts.separations
            );
        }
        for pair in starts.implied_f.windows(2) {
            assert!(pair[0] < pair[1], "{:?}", starts.implied_f);
        }
    }

    /// Every start is a real guess: a separation is a fraction of the outside rate and
    /// an inside fraction is a fraction of the genome, so both live strictly inside
    /// `(0, 1)`. A separation of 1 would make the two states identical, and one of 0
    /// would fix the inside state's heterozygote rate at zero — the tie the fit is
    /// written to avoid.
    #[test]
    fn every_default_start_lies_strictly_inside_the_unit_interval() {
        let starts = RunsModelStarts::default();

        for guess in starts.separations.iter().chain(starts.implied_f.iter()) {
            assert!(*guess > 0.0 && *guess < 1.0, "{guess} is not a real guess");
        }
        assert_eq!(starts.separations.len() * starts.implied_f.len(), 9);
    }

    /// **The two fields hold the same type and the same shape, so nothing but their
    /// values tells them apart.** Both are three ascending fractions in `(0, 1)`, so
    /// every property asserted above survives swapping them — and a swap puts the spread
    /// on the wrong axis, which is precisely the `F` = 0.0000 this type exists to
    /// prevent. The values are therefore pinned.
    #[test]
    fn the_default_starts_are_the_values_the_design_specifies_on_each_axis() {
        let starts = RunsModelStarts::default();

        assert_eq!(starts.separations.as_slice(), &[0.05, 1.0 / 3.0, 0.75]);
        assert_eq!(starts.implied_f.as_slice(), &[0.05, 0.5, 0.75]);
    }

    /// `starts_tried` is documented best-first, and `best_start` is the one reader of
    /// that ordering — so a call site cannot index the wrong end.
    #[test]
    fn best_start_reads_the_head_of_the_best_first_ordering() {
        let outcome = |inbreeding: f64, log_likelihood: f64| StartOutcome {
            separation: 1.0 / 3.0,
            implied_f: 0.5,
            inbreeding,
            log_likelihood,
        };
        let fit = RunsModelFit {
            outside_het: 0.0105,
            outside_hom_alt: 0.0010,
            inside_het_floor: 0.0008,
            inside_hom_alt: 0.0009,
            enter_run_per_base: 1e-7,
            leave_run_per_base: 3e-7,
            termination: FitTermination {
                iterations: 12,
                converged: true,
            },
            starts_tried: SmallVec::from_slice(&[
                outcome(0.2634, -1.20e9),
                outcome(0.0000, -1.41e9),
            ]),
            resolution: 0.01,
            undecided_windows: 0,
        };

        assert_eq!(fit.best_start().map(|s| s.inbreeding), Some(0.2634));
        assert!(
            RunsModelFit {
                starts_tried: SmallVec::new(),
                ..fit
            }
            .best_start()
            .is_none()
        );
    }

    /// **The spread is the full range over the *tied* starts** — not the gap between the
    /// best two, and not the range over all of them.
    ///
    /// Five properties, each a rule the others cannot see, written on outcomes in best-score
    /// order as the fit produces them.
    ///
    /// **These four outcomes are written down rather than fitted, and the two log-likelihood
    /// gaps in them are the measured ones** — 0.02 nats for a start that cannot be told from
    /// the best, 1,473 for one the data rejected (`research/inbreeding_resolution` §6.2).
    ///
    /// - The widest pair is deliberately **not** the top two: taking the best and
    ///   second-best would report 0.02 where the answer is 0.30.
    /// - The fourth start is 1,473 nats behind, so it is excluded. Including it would report
    ///   0.32 and refuse this fit, which is the failure the tie rule exists to prevent.
    /// - The three that are kept sit 0.00, 0.02 and 0.30 behind, which brackets the tie gap
    ///   from below: a rule that kept only the exactly-equal starts would report 0.00 here.
    /// - A second slice straddles the gap at 5 and 50 nats, which is what makes the **value**
    ///   of [`MAX_TIED_START_LOG_LIKELIHOOD_GAP`] load-bearing rather than merely bounded.
    ///   Without it the constant survives anywhere from 0.92 to 1,472.9 — the whole range its
    ///   build-time asserts permit — because no other fixture has a start between 0.3 and
    ///   1,473 nats behind.
    /// - A `NaN` score counts as tied, so the disagreement it carries refuses the fit. That
    ///   is stated in [`tied_starts`]' doc and nothing else checks it: inverting the
    ///   comparison to `!(gap <= limit)`, which reads as equivalent, drops the `NaN` start
    ///   instead.
    #[test]
    fn the_spread_is_the_range_over_the_starts_that_scored_alike() {
        let outcome = |inbreeding: f64, log_likelihood: f64| StartOutcome {
            separation: 1.0 / 3.0,
            implied_f: 0.5,
            inbreeding,
            log_likelihood,
        };
        let starts = [
            outcome(0.60, -1_802_299.7),
            outcome(0.58, -1_802_299.72),
            outcome(0.30, -1_802_300.0),
            outcome(0.28, -1_803_772.9),
        ];

        assert_eq!(tied_starts(&starts).len(), 3);
        assert!((spread_across_tied_starts(&starts) - 0.30).abs() < 1e-12);

        // Five nats behind is tied — beaten by 150 to 1, which settles nothing. Fifty is
        // not: 5 × 10²¹ to 1.
        let straddling = [
            outcome(0.30, -1_000.0),
            outcome(0.31, -1_005.0),
            outcome(0.90, -1_050.0),
        ];
        assert_eq!(tied_starts(&straddling).len(), 2);
        assert!((spread_across_tied_starts(&straddling) - 0.01).abs() < 1e-12);

        // An unreadable score cannot be told apart from the best, so it ties, and the
        // disagreement it brings with it refuses — spec §6.1's direction.
        let unreadable = [outcome(0.10, -1.0), outcome(0.90, f64::NAN)];
        assert_eq!(tied_starts(&unreadable).len(), 2);
        assert!(spread_across_tied_starts(&unreadable) > MAX_IDENTIFIED_START_SPREAD);

        // An unreadable *`F`* refuses too, and by a different route: `f64::max` ignores
        // `NaN`, so a fold would have dropped this start from the range and reported 0.0.
        assert!(
            spread_across_tied_starts(&[outcome(0.10, -1.0), outcome(f64::NAN, -1.0)]).is_nan()
        );

        // One start is no spread at all, and no start is none either — the fit refuses an
        // empty start set long before this.
        assert_eq!(spread_across_tied_starts(&[outcome(0.26, -1.0)]), 0.0);
        assert_eq!(spread_across_tied_starts(&[]), 0.0);
    }

    /// The two thresholds are the measured ones, and the values themselves are pinned
    /// rather than left to a doc comment that a later edit can contradict. **The bounds
    /// they have to sit between are checked at build time instead** — see the
    /// `const _: () = assert!` block beside the constants, which makes moving either one
    /// outside the measurement an `error[E0080]` rather than a red test.
    ///
    /// **A pin is a weak test and it is worth saying which part is weak.** It catches a
    /// value that was changed without the doc being changed; it says nothing about whether
    /// the value is right. What constrains the spread from both sides is
    /// `the_spread_threshold_answers_the_measured_spread_below_it_and_refuses_the_one_above`,
    /// and what constrains the tie gap from both sides is the straddling slice above.
    #[test]
    fn the_start_agreement_thresholds_are_the_measured_ones() {
        assert!((MAX_IDENTIFIED_START_SPREAD - 0.05).abs() < 1e-12);
        assert!((MAX_TIED_START_LOG_LIKELIHOOD_GAP - 10.0).abs() < 1e-12);
    }

    /// **Both sides of the threshold, on numbers the drawn genomes cannot reach.** Every
    /// fixture in this file either agrees to 3.7 × 10⁻¹¹ or disagrees by more than 0.1, so
    /// the whole interval between is untested by them — a review moved the constant to
    /// 0.001 and to 0.1146 and the entire suite stayed green both times.
    ///
    /// The two spreads below are measured ones, from an uncensored sweep of twenty-four
    /// genomes drawn with **no runs at all** at twenty heterozygotes a window: 0.0213 was
    /// answered and 0.0521 refused, and they are adjacent in that sweep's ordering. **What
    /// they bracket is a judgement, not a boundary in the data** — the no-runs spreads run
    /// 0.0213, 0.0264, 0.0323, 0.0337, 0.0483, 0.0521, 0.0541 … in an unbroken line, so
    /// there is no gap here for a threshold to sit in. These two pin where it was put.
    #[test]
    fn the_spread_threshold_answers_the_measured_spread_below_it_and_refuses_the_one_above() {
        let outcome = |inbreeding: f64| StartOutcome {
            separation: 1.0 / 3.0,
            implied_f: 0.5,
            inbreeding,
            log_likelihood: -1_000.0,
        };

        let answered = spread_across_tied_starts(&[outcome(0.30), outcome(0.30 + 0.0213)]);
        assert!(
            answered <= MAX_IDENTIFIED_START_SPREAD,
            "a spread of {answered} has to be inside the threshold of \
             {MAX_IDENTIFIED_START_SPREAD}; tightening it below this refuses a fit whose \
             starts nearly agree"
        );

        let refused = spread_across_tied_starts(&[outcome(0.30), outcome(0.30 + 0.0521)]);
        assert!(
            refused > MAX_IDENTIFIED_START_SPREAD,
            "a spread of {refused} has to be outside the threshold of \
             {MAX_IDENTIFIED_START_SPREAD}; loosening it past this answers a genome that \
             has no runs in it at all"
        );
    }

    /// The runs model's floor is stated where the model lives, not beside the
    /// site-count floor it is unlike: that one guards against a rate fitted from too
    /// little, this one against an answer that is entirely the estimator's own noise.
    #[test]
    fn the_window_floor_is_the_measured_one() {
        assert_eq!(MIN_WINDOWS_TO_FIT_INBREEDING, 3_000);
    }

    /// **The resolution reproduces the two figures the design docs quote**, which is what
    /// makes an interpolation of research note §3.6 one statement of that measurement
    /// rather than a second, disagreeing one.
    ///
    /// Both are asserted to two significant figures, which is the precision the docs state
    /// them at: "about 0.01 at tomato's 8,004 windows and 0.003 at a human genome's
    /// 31,000".
    #[test]
    fn the_resolution_reproduces_the_window_counts_the_design_quotes() {
        let tomato = resolution_at(8_004);
        let human = resolution_at(31_000);

        assert!(
            (tomato - 0.01).abs() < 0.0005,
            "tomato's 8,004 windows give {tomato}, where the design says about 0.01"
        );
        assert!(
            (human - 0.003).abs() < 0.0005,
            "a human genome's 31,000 give {human}, where the design says about 0.003"
        );
        // The measured points come back exactly, and the curve descends between them.
        assert!((resolution_at(4_800) - 0.017).abs() < 1e-12);
        assert!((resolution_at(19_200) - 0.0040).abs() < 1e-12);
        assert!(tomato > human, "more windows, finer resolution");
        assert!(
            resolution_at(200) >= resolution_at(1_200),
            "flat below the range"
        );
    }

    // -----------------------------------------------------------------
    // The fit, against a drawn genome whose realised autozygous fraction is known.
    // -----------------------------------------------------------------

    use std::sync::Arc;

    use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
    use crate::ng::parameter_estimation::generic::expected_counts::alternative_read_probability;
    use crate::ng::parameter_estimation::generic::histogram::DepthAndAltReads;
    use crate::ng::parameter_estimation::generic::{SampleRates, WindowIndex};
    use crate::ng::types::{Bp, ContigId, ErrorRate, GenotypeFrequency, ReadGroupId};

    /// SplitMix64 — the research harness's generator, so a drawn genome here and one there
    /// are drawn the same way.
    struct SplitMix64(u64);

    impl SplitMix64 {
        fn new(seed: u64) -> Self {
            Self(seed)
        }

        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        fn unit(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }

        /// How many of `n` trials succeed at probability `p`. Linear in `n`, which is fine
        /// because `n` here is a window's site count and the loop below draws one window's
        /// cells with a handful of calls.
        fn binomial(&mut self, n: u64, p: f64) -> u64 {
            if p <= 0.0 {
                return 0;
            }
            if p >= 1.0 {
                return n;
            }
            (0..n).filter(|_| self.unit() < p).count() as u64
        }
    }

    /// What a drawn genome is, and what its realised autozygous fraction turned out to be.
    ///
    /// **Realised and never nominal.** A finite genome does not have the `F` its transition
    /// rates imply — the states are drawn, so the fraction actually covered by runs is its
    /// own random variable — and scoring against the nominal value reads sampling as bias
    /// (spec §6.5).
    struct DrawnGenome {
        windowed: BTreeMap<(WindowKey, Ploidy), DepthAltHistogram<u32>>,
        realised_f: f64,
        windows: usize,
    }

    /// Tomato-shaped but scaled so a test can run it: **twelve contigs of 300 windows**,
    /// 3,600 in all against the floor of 3,000, and 400 sites a window at ten reads a site.
    ///
    /// **The heterozygosity is raised from tomato's 1 per kb to 5%, and that is what pays
    /// for the smaller window.** The chain classifies a window on how many heterozygotes it
    /// holds, so what matters is the *count*: the research harness gets ~100 per window from
    /// 100,000 sites at 1 per kb, and this gets ~20 from 400 sites at 5%. Fewer, so the
    /// tolerances below are looser than the harness's four decimal places — which is the
    /// harness's claim at full scale and not this fixture's.
    const CONTIGS: u32 = 12;
    const WINDOWS_PER_CONTIG: u32 = 300;
    const SITES_PER_WINDOW: u64 = 400;
    const WINDOW_DEPTH: u32 = 10;
    const ERROR_RATE: f64 = 0.001;
    const OUTSIDE_HET: f64 = 0.05;
    const HOM_ALT: f64 = 0.02;
    const INSIDE_HET: f64 = 0.002;
    /// One window in thirty ends a run, so a run is 3 Mb on average.
    const LEAVE_RUN: f64 = 1.0 / 30.0;

    fn ploidy2() -> Ploidy {
        Ploidy::try_new(2).expect("a positive copy number")
    }

    /// The fixture at its own shape — 400 sites a window at 5% heterozygous, which is about
    /// twenty heterozygotes a window.
    fn draw_genome(nominal_f: f64, false_het_floor: f64, seed: u64) -> DrawnGenome {
        draw_genome_at(
            SITES_PER_WINDOW,
            OUTSIDE_HET,
            nominal_f,
            false_het_floor,
            seed,
        )
    }

    /// **The same genome at a stated evidence level**, because how many heterozygotes a
    /// window holds is what the chain classifies it on, and one of the tests below needs a
    /// window that holds few. `sites_per_window × outside_het` is that count.
    fn draw_genome_at(
        sites_per_window: u64,
        outside_het: f64,
        nominal_f: f64,
        false_het_floor: f64,
        seed: u64,
    ) -> DrawnGenome {
        let edges = Arc::new(DepthBinEdges::new());
        let mut rng = SplitMix64::new(seed);
        let enter = LEAVE_RUN * nominal_f / (1.0 - nominal_f);
        let diploid = ploidy2();

        let state_frequencies = |inside: bool| -> [f64; GENOTYPES] {
            let het = if inside { INSIDE_HET } else { outside_het } + false_het_floor;
            [1.0 - het - HOM_ALT, het, HOM_ALT]
        };
        // The cell probabilities under each state: `Σ_j π_j · Binomial(k; depth, p_j(ε))`.
        let cell_probabilities = |inside: bool| -> Vec<f64> {
            let frequencies = state_frequencies(inside);
            (0..=WINDOW_DEPTH)
                .map(|alt_reads| {
                    let mut total = 0.0;
                    for (dosage, &frequency) in frequencies.iter().enumerate() {
                        let p = alternative_read_probability(dosage as u8, diploid, ERROR_RATE);
                        let mut term = (1.0 - p).powi(WINDOW_DEPTH as i32);
                        for step in 1..=alt_reads {
                            term *= f64::from(WINDOW_DEPTH - step + 1) / f64::from(step) * p
                                / (1.0 - p);
                        }
                        total += frequency * term;
                    }
                    total
                })
                .collect()
        };
        let inside_cells = cell_probabilities(true);
        let outside_cells = cell_probabilities(false);

        let mut windowed = BTreeMap::new();
        let mut inside_positions = 0u64;
        let mut total_positions = 0u64;

        for contig in 0..CONTIGS {
            // Each contig starts in the stationary distribution of the chain.
            let mut inside = rng.unit() < enter / (enter + LEAVE_RUN);
            for window in 0..WINDOWS_PER_CONTIG {
                let cells = if inside {
                    &inside_cells
                } else {
                    &outside_cells
                };
                let mut table = DepthAltHistogram::new(Arc::clone(&edges));

                // One multinomial draw over the window's cells, by sequential conditional
                // binomials — the exact draw. It is not cheaper per site than drawing
                // each site: `binomial(n, p)` consumes `n` uniforms, so a 400-site window
                // costs about 440 draws. What it buys is that no site is ever *placed*
                // one at a time, so the cell counts are exact rather than accumulated.
                let mut left = sites_per_window;
                let mut remaining: f64 = cells.iter().sum();
                for (alt_reads, &probability) in cells.iter().enumerate() {
                    if left == 0 {
                        break;
                    }
                    let in_cell = if remaining > 0.0 {
                        rng.binomial(left, (probability / remaining).clamp(0.0, 1.0))
                    } else {
                        0
                    };
                    remaining -= probability;
                    left -= in_cell;
                    for _ in 0..in_cell {
                        table
                            .add_site(DepthAndAltReads::new(WINDOW_DEPTH, alt_reads as u32), Bp(1));
                    }
                }

                total_positions += table.total_covered_positions();
                if inside {
                    inside_positions += table.total_covered_positions();
                }
                windowed.insert(
                    (
                        WindowKey::InWindow {
                            contig: ContigId(contig),
                            window: WindowIndex(window),
                        },
                        diploid,
                    ),
                    table,
                );

                inside = if inside {
                    rng.unit() >= LEAVE_RUN
                } else {
                    rng.unit() < enter
                };
            }
        }

        DrawnGenome {
            realised_f: inside_positions as f64 / total_positions as f64,
            windows: windowed.len(),
            windowed,
        }
    }

    fn sample_rates(false_het_floor: f64) -> SampleRates {
        sample_rates_at(OUTSIDE_HET, false_het_floor)
    }

    /// The frequencies the fit starts its outside state from, at a stated heterozygosity —
    /// the partner of [`draw_genome_at`], since a genome drawn at one rate and started at
    /// another is a different experiment.
    fn sample_rates_at(outside_het: f64, false_het_floor: f64) -> SampleRates {
        let het = outside_het + false_het_floor;
        SampleRates::try_new(
            ploidy2(),
            [1.0 - het - HOM_ALT, het, HOM_ALT]
                .into_iter()
                .map(|f| GenotypeFrequency::try_new(f).expect("a frequency"))
                .collect(),
        )
        .expect("a well-formed diploid set")
    }

    fn one_library() -> SampleLibraryNoise {
        SampleLibraryNoise::single(
            ReadGroupId(1),
            ErrorRate::try_new(ERROR_RATE).expect("a probability"),
        )
    }

    /// **E3's oracle: `F` recovers the realised autozygous fraction of a drawn genome**, at
    /// four nominal levels from 0.05 to 0.60.
    ///
    /// **Scored against the realised fraction and never the nominal one.** A finite genome
    /// does not have the `F` its transition rates imply, so comparing against the nominal
    /// value reads sampling as bias — the drawn fractions here miss their nominal values by
    /// as much as 0.027 while the fit is within 0.02 of what was drawn.
    ///
    /// **The tolerance is this fixture's, not the harness's.** The research note's four
    /// decimal places come from 8,004 windows of 100,000 sites each; 3,600 windows of 400
    /// sites hold about a fifth as many heterozygotes per window, so the chain has less to
    /// classify each one on. What is being checked here is that the estimator is the one
    /// the harness measured — that it recovers a drawn fraction rather than a nominal one,
    /// across a five-fold range of it.
    #[test]
    fn f_recovers_a_drawn_genomes_realised_autozygous_fraction() {
        let starts = RunsModelStarts::default();
        let mut widest_miss: f64 = 0.0;
        for (seed, nominal) in [0.05_f64, 0.15, 0.30, 0.60].into_iter().enumerate() {
            let genome = draw_genome(nominal, 0.0, 0xE3_0000 + seed as u64);
            assert!(
                genome.windows >= MIN_WINDOWS_TO_FIT_INBREEDING,
                "{} windows",
                genome.windows
            );
            widest_miss = widest_miss.max((genome.realised_f - nominal).abs());

            let (estimate, fit) = fit_inbreeding(
                "drawn",
                &genome.windowed,
                &one_library(),
                ploidy2(),
                &sample_rates(0.0),
                &starts,
            )
            .expect("a drawn genome of 3,600 windows is fittable");

            let fitted = estimate.value.get();
            assert!(
                (fitted - genome.realised_f).abs() < 0.02,
                "nominal {nominal}: fitted {fitted} against a realised {}",
                genome.realised_f
            );
            assert!(
                fit.inside_het_floor < fit.outside_het,
                "the ordering constraint puts the lower heterozygote rate inside: \
                 inside {}, outside {}",
                fit.inside_het_floor,
                fit.outside_het
            );
            assert_eq!(fit.starts_tried.len(), 9);
            // **The margin `MAX_IDENTIFIED_START_SPREAD` is set against, on this file's own
            // genomes rather than only on the research harness's.** All nine starts landing
            // together is what makes a twentieth a threshold with three orders of magnitude
            // to spare rather than a tuned constant — and it is what a mutation that made
            // the starts disagree would break here before it broke anything downstream.
            assert!(
                fit.spread_across_tied_starts() < 0.0001,
                "nominal {nominal}: the nine starts spread by {}",
                fit.spread_across_tied_starts()
            );
            // **A spread of zero over one start is not agreement, it is an empty
            // comparison** — and both print 0.0000. Asserting the tie set's size is what
            // makes the line above say something; at a realised `F` of 0.95 under a
            // five-fold false-heterozygote floor only one start in nine ties, and there the
            // criterion contributes nothing.
            assert_eq!(
                fit.tied_starts().len(),
                9,
                "nominal {nominal}: the margin above is claimed over all nine starts"
            );
        }

        // **The premise, asserted rather than asserted about.** Scoring against the drawn
        // fraction only says something different from scoring against the nominal one if
        // the two differ by more than the tolerance above — and the doc's claim about how
        // far they differ is exactly the kind of number that goes stale silently.
        assert!(
            widest_miss > 0.02,
            "the drawn fractions sit within the tolerance of their nominal values \
             ({widest_miss}), so this test would pass against either and says nothing \
             about which one it scores"
        );
        assert!(
            widest_miss < 0.05,
            "the widest miss is {widest_miss}, and the doc above says 0.027"
        );
    }

    /// **A floor of false heterozygotes does not move `F`.** Collapsed paralogs and
    /// mismapping lift *both* states together, so what the chain reads — the gap between
    /// them — is unchanged. That is spec §6.2's second reason for choosing this estimator
    /// over the ratio one, and the research harness measures it up to five times the real
    /// rate (§3.3).
    ///
    /// The floor here runs to **five times** the outside heterozygote rate.
    ///
    /// **The four levels are not the same genome, though they share a seed**, and saying so
    /// matters because the comparison would be a different one if they were. Raising the
    /// floor changes the cell probabilities, which changes how many draws each window's
    /// multinomial consumes, which shifts the whole downstream stream — the realised
    /// fractions come out 0.2056, 0.3153, 0.3253 and 0.3122. So this is not "one genome
    /// fitted four times"; it is four genomes of the same construction, each scored against
    /// **its own** realised fraction, which is what the assertion below does.
    #[test]
    fn a_floor_of_false_heterozygotes_does_not_move_f() {
        let starts = RunsModelStarts::default();
        let mut fitted_at = Vec::new();
        for multiple in [0.0_f64, 1.0, 3.0, 5.0] {
            let floor = multiple * OUTSIDE_HET;
            let genome = draw_genome(0.30, floor, 0xE3_1000);

            let (estimate, _) = fit_inbreeding(
                "floored",
                &genome.windowed,
                &one_library(),
                ploidy2(),
                &sample_rates(floor),
                &starts,
            )
            .expect("fittable");

            fitted_at.push((multiple, estimate.value.get(), genome.realised_f));
        }

        for &(multiple, fitted, realised) in &fitted_at {
            assert!(
                (fitted - realised).abs() < 0.05,
                "at {multiple}× the real heterozygote rate the fit gave {fitted} against a \
                 realised {realised}; all levels: {fitted_at:?}"
            );
        }
    }

    /// **Starts that share one separation guess return a failed search, not a zero.** This
    /// is the one place the estimator produces a confident wrong number, and the error path
    /// is the only thing between it and a caller.
    ///
    /// The genome is 30% covered by runs with a floor of false heterozygotes three times
    /// the real rate, which is what pulls the two states close together — the shape spec
    /// §6.5 measures. Handed three starts that all guess the inside state at a twentieth of
    /// the outside rate, every one of them empties the inside state.
    #[test]
    fn starts_sharing_one_separation_report_a_failed_search_rather_than_zero() {
        let floor = 3.0 * OUTSIDE_HET;
        let genome = draw_genome(0.30, floor, 0xE3_2000);
        let narrow = RunsModelStarts {
            separations: SmallVec::from_slice(&[0.05]),
            implied_f: SmallVec::from_slice(&[0.05, 0.5, 0.75]),
        };

        let narrow_result = fit_inbreeding(
            "narrow",
            &genome.windowed,
            &one_library(),
            ploidy2(),
            &sample_rates(floor),
            &narrow,
        );
        let spread_result = fit_inbreeding(
            "spread",
            &genome.windowed,
            &one_library(),
            ploidy2(),
            &sample_rates(floor),
            &RunsModelStarts::default(),
        );

        let spread = spread_result.expect("the default starts span the separation");
        assert!(
            (spread.0.value.get() - genome.realised_f).abs() < 0.05,
            "the spread starts gave {} against a realised {}",
            spread.0.value.get(),
            genome.realised_f
        );
        let narrow_error =
            narrow_result.expect_err("the narrow starts must not return a number at all");
        assert!(
            matches!(
                narrow_error,
                ParameterEstimationError::InbreedingStatesNotSeparated { starts: 3, .. }
            ),
            "the wrong failure: {narrow_error}"
        );
        assert!(
            narrow_error
                .to_string()
                .contains("not an inbreeding coefficient"),
            "the message has to say what it is not: {narrow_error}"
        );
    }

    /// Below the window floor the fit fails rather than emitting, and the message names the
    /// count and the floor.
    #[test]
    fn a_genome_below_the_window_floor_is_refused() {
        let genome = draw_genome(0.30, 0.0, 0xE3_3000);
        let few: BTreeMap<(WindowKey, Ploidy), DepthAltHistogram<u32>> = genome
            .windowed
            .into_iter()
            .take(MIN_WINDOWS_TO_FIT_INBREEDING - 1)
            .collect();

        let error = fit_inbreeding(
            "small",
            &few,
            &one_library(),
            ploidy2(),
            &sample_rates(0.0),
            &RunsModelStarts::default(),
        )
        .expect_err("2,999 windows is below the floor");

        let message = error.to_string();
        assert!(message.contains("2999"), "{message}");
        assert!(message.contains("3000"), "{message}");
    }

    /// **`F` is weighted by reference positions, not by loci** — so a window dense in
    /// widened indel loci is not under-weighted. Two windows holding the same number of
    /// loci but different spans must not count equally.
    ///
    /// Asserted on the weighting directly, because the two coincide on every fixture above:
    /// each locus there spans one position, which is what makes a chain built from them
    /// unable to tell the two rules apart.
    #[test]
    fn the_occupancy_is_weighted_by_covered_positions_and_not_by_loci() {
        let edges = Arc::new(DepthBinEdges::new());
        let diploid = ploidy2();
        let mut narrow = DepthAltHistogram::new(Arc::clone(&edges));
        let mut wide = DepthAltHistogram::new(Arc::clone(&edges));
        for _ in 0..10 {
            narrow.add_site(DepthAndAltReads::new(WINDOW_DEPTH, 0), Bp(1));
            wide.add_site(DepthAndAltReads::new(WINDOW_DEPTH, 0), Bp(50));
        }

        let windowed = BTreeMap::from([
            (
                (
                    WindowKey::InWindow {
                        contig: ContigId(0),
                        window: WindowIndex(0),
                    },
                    diploid,
                ),
                narrow,
            ),
            (
                (
                    WindowKey::InWindow {
                        contig: ContigId(0),
                        window: WindowIndex(1),
                    },
                    diploid,
                ),
                wide,
            ),
        ]);

        let chain = WindowChain::over(&windowed, diploid, &one_library());

        assert_eq!(chain.windows(), 2);
        assert_eq!(chain.covered_positions_total(), 10 + 500);
        assert_eq!(chain.windows_holding_sites, 2);
    }

    /// **Absent windows are in the chain**, so the model cannot step across an unmappable
    /// megabase in one transition, and the chain restarts at every contig boundary rather
    /// than running the end of one chromosome into the start of the next.
    ///
    /// Two contigs, each with a hole and a leading gap: the chain holds every window from
    /// zero to the last one that received a site, and only the two that did count towards
    /// the floor.
    #[test]
    fn the_chain_covers_every_window_of_a_contig_and_restarts_at_each_boundary() {
        let edges = Arc::new(DepthBinEdges::new());
        let diploid = ploidy2();
        let mut windowed = BTreeMap::new();
        for (contig, window) in [(0u32, 3u32), (1, 5)] {
            let mut table = DepthAltHistogram::new(Arc::clone(&edges));
            table.add_site(DepthAndAltReads::new(WINDOW_DEPTH, 0), Bp(1));
            windowed.insert(
                (
                    WindowKey::InWindow {
                        contig: ContigId(contig),
                        window: WindowIndex(window),
                    },
                    diploid,
                ),
                table,
            );
        }

        let chain = WindowChain::over(&windowed, diploid, &one_library());

        // Windows 0..=3 of the first contig and 0..=5 of the second.
        assert_eq!(chain.windows(), 4 + 6);
        assert_eq!(chain.windows_holding_sites, 2);
        assert_eq!(chain.contig_starts, vec![0, 4]);
        assert!(chain.opens_a_contig(0) && chain.opens_a_contig(4));
        assert!(!chain.opens_a_contig(3) && !chain.opens_a_contig(5));
    }

    /// A window that received no site emits nothing under either state — which is what lets
    /// the chain's own transitions decide it, rather than its emptiness arguing for one
    /// state over the other.
    #[test]
    fn an_absent_window_emits_the_same_under_both_states() {
        let edges = Arc::new(DepthBinEdges::new());
        let diploid = ploidy2();
        let mut table = DepthAltHistogram::new(Arc::clone(&edges));
        table.add_site(DepthAndAltReads::new(WINDOW_DEPTH, 0), Bp(1));
        let windowed = BTreeMap::from([(
            (
                WindowKey::InWindow {
                    contig: ContigId(0),
                    window: WindowIndex(1),
                },
                diploid,
            ),
            table,
        )]);
        let chain = WindowChain::over(&windowed, diploid, &one_library());
        let mut scratch = [0.0; GENOTYPES];

        let inside = [1.0 - INSIDE_HET - HOM_ALT, INSIDE_HET, HOM_ALT];
        let outside = [1.0 - OUTSIDE_HET - HOM_ALT, OUTSIDE_HET, HOM_ALT];

        assert_eq!(chain.ln_emission(0, &inside, &mut scratch), 0.0);
        assert_eq!(chain.ln_emission(0, &outside, &mut scratch), 0.0);
        assert_ne!(
            chain.ln_emission(1, &inside, &mut scratch),
            chain.ln_emission(1, &outside, &mut scratch),
            "the window that holds a site must tell the two states apart"
        );
    }

    /// A start's separation and implied fraction have to be real guesses. A separation of
    /// one makes the two states identical, at which point the likelihood is exactly flat in
    /// `F` and what a fit returns is its own starting point.
    #[test]
    #[should_panic(expected = "is not a fraction of the outside heterozygote rate")]
    fn a_separation_of_one_is_refused() {
        let _ = RunsModelStart::new(1.0, 0.5, [0.93, 0.05, 0.02]);
    }

    #[test]
    #[should_panic(expected = "assumes the answer")]
    fn a_start_assuming_the_whole_genome_is_inside_a_run_is_refused() {
        let _ = RunsModelStart::new(0.05, 1.0, [0.93, 0.05, 0.02]);
    }

    /// A start's two transition rates come from its implied fraction and the stated mean
    /// run length, so that the stationary inside probability of the chain it starts from is
    /// the fraction it was asked for.
    #[test]
    fn a_starts_transition_rates_imply_the_fraction_it_was_given() {
        for implied_f in [0.05, 0.5, 0.75] {
            let start = RunsModelStart::new(1.0 / 3.0, implied_f, [0.93, 0.05, 0.02]);
            let stationary = start.enter_run / (start.enter_run + start.leave_run);

            assert!(
                (stationary - implied_f).abs() < 1e-12,
                "a start asked for {implied_f} begins at {stationary}"
            );
            assert!((start.leave_run - 1.0 / MEAN_RUN_WINDOWS_AT_START).abs() < 1e-12);
        }
    }

    /// **The ordering constraint fires when the fit lands the wrong way round, and moves
    /// everything indexed by state together.** Deleting it left every drawn genome green,
    /// because every start builds its inside state *below* the outside one and no fixture
    /// crosses over — so this is where the constraint is checked.
    ///
    /// The two transition rates and the window posteriors are swapped as well as the
    /// frequencies: a relabelling that moved only the frequencies would report a run's
    /// entry rate as its exit rate and `1 − F` as `F`, which is the same wrong answer by a
    /// longer route.
    #[test]
    fn the_ordering_constraint_swaps_every_thing_indexed_by_state() {
        // Inverted: the state at `INSIDE` carries the *higher* heterozygote rate.
        let mut frequencies = [[0.90, 0.08, 0.02], [0.97, 0.01, 0.02]];
        let mut enter = 0.01;
        let mut leave = 0.05;
        let mut posteriors = [[0.9, 0.1], [0.2, 0.8]];

        assert!(relabel_if_inverted(
            &mut frequencies,
            &mut enter,
            &mut leave,
            &mut posteriors
        ));

        assert_eq!(frequencies, [[0.97, 0.01, 0.02], [0.90, 0.08, 0.02]]);
        assert!((enter - 0.05).abs() < 1e-15 && (leave - 0.01).abs() < 1e-15);
        assert_eq!(posteriors, [[0.1, 0.9], [0.8, 0.2]]);
    }

    /// And it leaves a fit that already sits the right way round alone — including when
    /// the two states' heterozygote rates are *equal*, where a `<` would swap on every
    /// call and a `<=` does not.
    #[test]
    fn the_ordering_constraint_leaves_a_correctly_labelled_fit_alone() {
        let ordered = [[0.97, 0.01, 0.02], [0.90, 0.08, 0.02]];
        let mut frequencies = ordered;
        let mut enter = 0.01;
        let mut leave = 0.05;
        let mut posteriors = [[0.9, 0.1]];

        assert!(!relabel_if_inverted(
            &mut frequencies,
            &mut enter,
            &mut leave,
            &mut posteriors
        ));
        assert_eq!(frequencies, ordered);

        let mut tied = [[0.94, 0.04, 0.02], [0.94, 0.04, 0.02]];
        assert!(!relabel_if_inverted(
            &mut tied,
            &mut enter,
            &mut leave,
            &mut posteriors
        ));
    }

    /// **The start kept is the best-scoring one that separated the states, not the
    /// best-scoring one.** A collapsed fit can outscore a real one, because its inside
    /// state is free to take whatever frequencies fit the whole genome — so "best" alone
    /// would return the zero the error path exists to refuse.
    ///
    /// On every drawn genome in this file all nine starts land on the same answer with the
    /// same score, so the rule survives being replaced by "keep the first" there. Here the
    /// outcomes are written down: the top-scoring one collapsed, and the answer must be the
    /// separated one below it.
    #[test]
    fn the_start_kept_is_the_best_scoring_one_that_separated_the_states() {
        let outcome = |inbreeding: f64, log_likelihood: f64, separated: bool| StateFit {
            frequencies: [[0.97, 0.01, 0.02], [0.90, 0.08, 0.02]],
            enter_run: 0.01,
            leave_run: 0.05,
            inbreeding,
            log_likelihood,
            termination: FitTermination {
                iterations: 4,
                converged: true,
            },
            undecided_windows: 0,
            separated_states: separated,
        };
        let start = RunsModelStart::new(1.0 / 3.0, 0.5, [0.90, 0.08, 0.02]);
        // Best first, as `starts_tried` is ordered — and the best one collapsed.
        let outcomes = vec![
            (start, outcome(0.0, -1.0e9, false)),
            (start, outcome(0.2634, -1.2e9, true)),
            (start, outcome(0.2600, -1.3e9, true)),
        ];

        let kept = best_separated_start(&outcomes).expect("one of them separated");

        assert!(
            (kept.inbreeding - 0.2634).abs() < 1e-15,
            "kept {} — the top scorer collapsed, and the best separated one is 0.2634",
            kept.inbreeding
        );
        assert!(best_separated_start(&[]).is_none());
        assert!(
            best_separated_start(&[(start, outcome(0.0, -1.0e9, false))]).is_none(),
            "a set in which nothing separated has no answer, which is what the error \
             variant is for"
        );
    }

    /// The inbreeding coefficient is diploid-only: above two copies it needs several
    /// identity-by-descent coefficients, and below two there are no heterozygotes.
    #[test]
    #[should_panic(expected = "diploid-only")]
    fn a_non_diploid_ploidy_is_refused() {
        let _ = fit_inbreeding(
            "tetraploid",
            &BTreeMap::new(),
            &one_library(),
            Ploidy::try_new(4).expect("a positive copy number"),
            &sample_rates(0.0),
            &RunsModelStarts::default(),
        );
    }

    /// Build a chain of alternating blocks whose windows differ in how many reference
    /// positions they cover, so that the four ways of averaging a window posterior give
    /// four different numbers.
    ///
    /// `wide_block_is_outside` puts the fifty-positions-a-locus windows on the
    /// heterozygote-rich state; otherwise they sit on the heterozygote-free one.
    fn blocks_of_unequal_span(
        contigs: u32,
        windows_per_block: u32,
    ) -> BTreeMap<(WindowKey, Ploidy), DepthAltHistogram<u32>> {
        let edges = Arc::new(DepthBinEdges::new());
        let diploid = ploidy2();
        let mut windowed = BTreeMap::new();
        for contig in 0..contigs {
            for block in 0..3u32 {
                let autozygous = block != 1;
                let span = if autozygous { Bp(1) } else { Bp(50) };
                for offset in 0..windows_per_block {
                    let mut table = DepthAltHistogram::new(Arc::clone(&edges));
                    for site in 0..200u32 {
                        let alt = if !autozygous && site % 2 == 0 { 5 } else { 0 };
                        table.add_site(DepthAndAltReads::new(WINDOW_DEPTH, alt), span);
                    }
                    windowed.insert(
                        (
                            WindowKey::InWindow {
                                contig: ContigId(contig),
                                window: WindowIndex(block * windows_per_block + offset),
                            },
                            diploid,
                        ),
                        table,
                    );
                }
            }
        }
        windowed
    }

    /// **`F` is the coverage-weighted posterior occupancy, and this is the only test that
    /// says so.** Every other fixture in this module gives each window the same number of
    /// sites and one reference position per site, so coverage weighting, site-count
    /// weighting, a plain mean over windows and the transition rates' ratio are all the
    /// same number there — a fit cannot tell them apart, and three separate defects
    /// injected into `fit_from` left the whole module green.
    ///
    /// One contig, three blocks of twenty windows: heterozygote-free, heterozygote-rich,
    /// heterozygote-free. The rich block spans fifty reference positions a locus and the
    /// other two span one. So forty of sixty windows are inside a run, but they hold
    /// 8,000 of 208,000 covered positions:
    ///
    /// | how the occupancy is averaged | `F` |
    /// |---|---:|
    /// | by covered positions — the definition | 0.038 |
    /// | by site count | 0.667 |
    /// | by window, unweighted | 0.667 |
    /// | `enter / (enter + leave)` — spec §6.1's "their ratio" | 0.66 |
    #[test]
    fn f_is_the_coverage_weighted_occupancy_and_not_the_three_things_it_coincides_with() {
        let windowed = blocks_of_unequal_span(1, 20);
        let chain = WindowChain::over(&windowed, ploidy2(), &one_library());
        let start = RunsModelStart::new(1.0 / 3.0, 0.5, [0.70, 0.28, 0.02]);
        let fitted = chain.fit_from(&start, MAX_RUNS_MODEL_ITERATIONS);

        let by_covered_positions = 8_000.0 / 208_000.0;
        let stationary = fitted.enter_run / (fitted.enter_run + fitted.leave_run);
        println!(
            "F {:.5}; covered-position {by_covered_positions:.5}; site-count 0.66667; \
             stationary {stationary:.5}",
            fitted.inbreeding
        );

        assert!(
            (fitted.inbreeding - by_covered_positions).abs() < 0.01,
            "F came back {}, not the covered-position occupancy {by_covered_positions}; \
             the site-count and unweighted averages are 0.6667 and the transition rates' \
             ratio is {stationary}",
            fitted.inbreeding
        );
    }

    /// **The resolution belongs to the windows that hold evidence, not to the padded
    /// chain.** A region-restricted run pads every contig from window zero, so the chain
    /// can be several times longer than the evidence in it — and `resolution` is what a
    /// consumer compares `F` against to decide *nothing detected*. Reading it off the
    /// padded length understates the floor.
    ///
    /// Here 3,600 windows of evidence sit beside one contig whose only site is at window
    /// 8,000: 3,601 windows hold sites, 11,601 are in the chain, and the two floors are
    /// 0.029 against 0.0068 — a factor of four.
    #[test]
    fn the_resolution_is_read_from_the_windows_that_hold_sites() {
        let genome = draw_genome(0.30, 0.0, 0xE3_0000);
        let mut windowed = genome.windowed;
        let edges = Arc::new(DepthBinEdges::new());
        let mut lone = DepthAltHistogram::new(Arc::clone(&edges));
        lone.add_site(DepthAndAltReads::new(WINDOW_DEPTH, 0), Bp(1));
        windowed.insert(
            (
                WindowKey::InWindow {
                    contig: ContigId(CONTIGS),
                    window: WindowIndex(8_000),
                },
                ploidy2(),
            ),
            lone,
        );

        let (_, fit) = fit_inbreeding(
            "padded",
            &windowed,
            &one_library(),
            ploidy2(),
            &sample_rates(0.0),
            &RunsModelStarts::default(),
        )
        .expect("3,601 windows of evidence is fittable");

        assert!(
            (fit.resolution - resolution_at(3_601)).abs() < 1e-9,
            "resolution came back {} — the padded chain's {} rather than the {} of the \
             3,601 windows that hold sites",
            fit.resolution,
            resolution_at(11_601),
            resolution_at(3_601)
        );
    }

    /// The window floor counts windows that **hold sites**, so a region-restricted run
    /// cannot pass it on padding. Ten windows of evidence, and one at index 9,999 that
    /// pads the contig to ten thousand.
    #[test]
    fn the_window_floor_cannot_be_passed_on_padding() {
        let edges = Arc::new(DepthBinEdges::new());
        let diploid = ploidy2();
        let mut windowed = BTreeMap::new();
        for window in [0u32, 1, 2, 3, 4, 5, 6, 7, 8, 9_999] {
            let mut table = DepthAltHistogram::new(Arc::clone(&edges));
            table.add_site(DepthAndAltReads::new(WINDOW_DEPTH, 0), Bp(1));
            windowed.insert(
                (
                    WindowKey::InWindow {
                        contig: ContigId(0),
                        window: WindowIndex(window),
                    },
                    diploid,
                ),
                table,
            );
        }

        let error = fit_inbreeding(
            "padded",
            &windowed,
            &one_library(),
            diploid,
            &sample_rates(0.0),
            &RunsModelStarts::default(),
        )
        .expect_err("ten windows of evidence is below the floor, however much padding");

        assert!(
            matches!(
                error,
                ParameterEstimationError::InbreedingNotFittable { windows: 10, .. }
            ),
            "the floor counted padding: {error}"
        );
    }

    /// **A genome with no runs at all must not come back near one**, which is the failure
    /// this estimator produces from the side the collapsed-search check cannot see.
    ///
    /// Research note §3.1 is a proof rather than a measurement: when the two states carry
    /// the same genotype frequencies the likelihood is **exactly flat** in `F`, so what a
    /// fit returns is an accident of where it started. Measured over eight genomes drawn
    /// with no runs, two came back at `F` = 0.9995 and 0.8475 — and both are exactly the
    /// fits whose two states coincided, at 0.968 and 0.929 of each other, where every
    /// legitimate fit in this file sits at 0.842 or below.
    ///
    /// **Swept over twenty-four seeds after the check went in**, of which twenty-three drew
    /// no runs at all: **nine are now refused and fourteen answered, and the largest `F` any
    /// of them returns is 0.042** — against a reported resolution of 0.029, and in line with
    /// research note §3.6's own floor for this window count (mean 0.017, worst seed 0.086 at
    /// 4,800 windows). Nothing comes back near one any more.
    ///
    /// **That nine in twenty-three are refused is the check working, not over-firing.** On a
    /// genome with no runs the two states genuinely coincide, so `F` is not identified at
    /// all (§3.1's proof), and spec §6.1's instruction for this parameter is *fail rather
    /// than emit*. It does change what an outbred sample looks like from the outside: an
    /// error saying *supply `F`* rather than a small number. Milestone G3, which expects
    /// `F` ≈ 0 on HG002, should expect either.
    ///
    /// The seed below is the worse of the two original failures. **One and not both**,
    /// because a fit whose states never separate runs to the iteration cap from all nine
    /// starts, and the second seed costs another sixty seconds of the suite to say the same
    /// thing — it is refused too, measured in the sweep above.
    #[test]
    fn a_genome_with_no_runs_does_not_come_back_almost_entirely_autozygous() {
        let seed = 3u64;
        let genome = draw_genome(0.001, 0.0, 0xE3_9000 + seed * 7919);
        assert_eq!(
            genome.realised_f, 0.0,
            "the seed has to have drawn no runs at all for this to say anything"
        );

        let error = fit_inbreeding(
            "outbred",
            &genome.windowed,
            &one_library(),
            ploidy2(),
            &sample_rates(0.0),
            &RunsModelStarts::default(),
        )
        .expect_err("a genome with no runs must not come back with a number");

        assert!(
            matches!(
                error,
                ParameterEstimationError::InbreedingStatesNotSeparated { .. }
            ),
            "the wrong failure: {error}"
        );
        assert!(
            error.to_string().contains("not an inbreeding coefficient"),
            "the message has to say what it is not: {error}"
        );
    }

    /// **The failure the other two checks are blind to, and the check that sees it.**
    ///
    /// A genome with **no runs at all** at five heterozygotes a window — a quarter of what
    /// the fixtures above carry, and what a shallow or region-restricted run leaves. Both
    /// separation checks pass: some window's posterior favours each state, and the two
    /// states come out far enough apart to clear [`MAX_IDENTIFIED_STATE_RATIO`]. The chain
    /// has manufactured a separation out of sampling noise, and it looks entirely real
    /// (`research/inbreeding_resolution_2026-08-09.md` §2).
    ///
    /// **That this error and not another one is the assertion.** The disagreement check runs
    /// last, so reaching it is itself the proof that the other two passed — which is why the
    /// test asserts the variant rather than merely asserting a refusal.
    ///
    /// **Eight seeds were swept at this shape and seven are refused here**, with spreads from
    /// 0.1147 to 0.8494; this is the widest. The eighth is the test below it: the one seed
    /// that drew a run.
    #[test]
    fn starts_that_land_in_different_places_are_refused_rather_than_averaged() {
        let genome = draw_genome_at(200, 0.025, 0.001, 0.0, 0xE3_A000 + 7 * 7919);
        assert_eq!(
            genome.realised_f, 0.0,
            "the seed has to have drawn no runs at all for this to say anything"
        );

        let error = fit_inbreeding(
            "noisy",
            &genome.windowed,
            &one_library(),
            ploidy2(),
            &sample_rates_at(0.025, 0.0),
            &RunsModelStarts::default(),
        )
        .expect_err("nine starts landing 0.85 apart have not identified an F");

        let ParameterEstimationError::InbreedingStartsDisagree {
            sample,
            tied_starts,
            starts,
            spread,
            threshold,
        } = &error
        else {
            panic!("the wrong failure, so this says nothing about the spread: {error}");
        };
        assert!(
            *spread > 0.5,
            "the nine starts have to land far apart for this fixture to be the failure it \
             is named for: {spread}"
        );
        assert_eq!(sample, "noisy");
        assert_eq!(*starts, 9);
        assert_eq!(
            *tied_starts, 9,
            "every start scored within the tie gap on this fixture, which is what makes the \
             disagreement a disagreement rather than a set of rejected competitors"
        );
        assert!((*threshold - MAX_IDENTIFIED_START_SPREAD).abs() < 1e-12);

        // **The message, which is the whole point of the variant.** A caller who reads this
        // as *nothing was found* will divide the cohort's diversity by `1 − F` with a
        // defaulted `F`, which is the harm the third refusal exists to prevent. Its sibling
        // `InbreedingStatesNotSeparated` is guarded on exactly this and it must be too.
        let message = error.to_string();
        assert!(
            message.contains("not an inbreeding coefficient of zero"),
            "the message has to say what it is not: {message}"
        );
        assert!(
            message.contains("supply F"),
            "the message has to say what to do about it: {message}"
        );
        assert!(
            message.contains("noisy"),
            "a cohort run of hundreds needs the sample named: {message}"
        );
    }

    /// **And the same evidence level answers where a run is genuinely there.** The one seed
    /// in eight at this shape that drew any run at all drew a small one — 0.28% of the
    /// genome — and the fit returns it to four decimal places with every start agreeing.
    ///
    /// This is the half that stops the criterion being a way to refuse everything hard: five
    /// heterozygotes a window is where the estimator fails on a genome with nothing to find,
    /// and it is *not* where it fails on a genome with something to find. The margin is the
    /// whole threshold — the starts agree to 0.0000 against a criterion of
    /// [`MAX_IDENTIFIED_START_SPREAD`].
    #[test]
    fn a_real_run_at_the_same_evidence_level_is_still_answered_and_the_starts_agree() {
        let genome = draw_genome_at(200, 0.025, 0.001, 0.0, 0xE3_A000 + 3 * 7919);
        assert!(
            genome.realised_f > 0.0,
            "the seed has to have drawn a run for this to say anything: {}",
            genome.realised_f
        );

        let (estimate, fit) = fit_inbreeding(
            "small run",
            &genome.windowed,
            &one_library(),
            ploidy2(),
            &sample_rates_at(0.025, 0.0),
            &RunsModelStarts::default(),
        )
        .expect("a genome that has a run is identified even at five heterozygotes a window");

        assert!(
            (estimate.value.get() - genome.realised_f).abs() < 0.001,
            "fitted {} against a realised {}",
            estimate.value.get(),
            genome.realised_f
        );
        assert!(
            fit.spread_across_tied_starts() < 0.0001,
            "the starts have to agree for the margin claimed on \
             MAX_IDENTIFIED_START_SPREAD to be this fixture's: {}",
            fit.spread_across_tied_starts()
        );
        assert_eq!(
            fit.tied_starts().len(),
            9,
            "all nine tie here, so the agreement above is over all of them"
        );
    }

    /// And neither check refuses a fit that **is** identified — including the hardest one
    /// this file has, where a floor of spurious heterozygotes five times the real rate pulls
    /// the two states to 0.842 of each other. A threshold set too tight on either would turn
    /// spec §6.2's own robustness claim into an error.
    ///
    /// **This fixture is why the spread is taken over the tied starts and not over all of
    /// them.** Six of the nine agree on `F` = 0.3157 against a realised 0.3122; the other
    /// three collapse to zero. The raw range is therefore 0.3157 — six times the threshold —
    /// and a criterion that used it would refuse the very fit spec §6.2 requires to survive.
    /// The three collapsed starts score **1,473 nats worse**, so the data has already
    /// rejected them, and the tie rule is what says so. Both numbers are asserted below,
    /// because it is their *difference* that carries the argument: a mutation that widened
    /// the tie gap past 1,473 would leave every other test in this file green.
    #[test]
    fn neither_check_refuses_the_closest_legitimate_fit() {
        let floor = 5.0 * OUTSIDE_HET;
        let genome = draw_genome(0.30, floor, 0xE3_1000);

        let (estimate, fit) = fit_inbreeding(
            "floored",
            &genome.windowed,
            &one_library(),
            ploidy2(),
            &sample_rates(floor),
            &RunsModelStarts::default(),
        )
        .expect("a fit at five times the real heterozygote rate is still identified");

        let ratio = fit.inside_het_floor / fit.outside_het;
        assert!(
            ratio > 0.8 && ratio <= MAX_IDENTIFIED_STATE_RATIO,
            "the closest legitimate fit sits at {ratio}, which has to be inside the \
             threshold of {MAX_IDENTIFIED_STATE_RATIO} for this test to say anything"
        );
        assert!((estimate.value.get() - genome.realised_f).abs() < 0.05);

        assert!(
            fit.spread_across_tied_starts() < 0.0001,
            "the starts that scored alike have to agree: {}",
            fit.spread_across_tied_starts()
        );
        assert_eq!(
            fit.tied_starts().len(),
            6,
            "six of the nine tie and three are the rejected collapse — a tie set of one \
             would satisfy the line above while comparing a start with itself"
        );
        let raw_range = fit
            .starts_tried
            .iter()
            .map(|start| start.inbreeding)
            .fold(f64::NEG_INFINITY, f64::max)
            - fit
                .starts_tried
                .iter()
                .map(|start| start.inbreeding)
                .fold(f64::INFINITY, f64::min);
        assert!(
            raw_range > MAX_IDENTIFIED_START_SPREAD,
            "the range over *all* the starts is {raw_range}, which is inside the threshold \
             of {MAX_IDENTIFIED_START_SPREAD} — so this fixture no longer says that taking \
             the spread over the tied starts is what saves it"
        );
    }
}
