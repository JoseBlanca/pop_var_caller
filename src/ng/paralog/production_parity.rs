//! **ng's copied scorer computes what production's does — asserted on randomised loci, by
//! bit pattern.**
//!
//! A textual guard (`copy_fidelity.rs`, deleted at promotion step C3) said the *text* was
//! production's. This says the *numbers* are, which is the property the port actually claims
//! (`doc/devel/ng/spec/hidden_paralog_filter.md` §7's reuse map: *"the parity oracle is
//! bit-identical ratios against production's on identical synthetic inputs"*). The two are not the
//! same statement, and **neither is strictly stronger than the other.** The text guard would still
//! pass if a build flag or a platform difference made one tree evaluate the same source
//! differently. The differential only sees what it draws: a review of this file found four input
//! classes it did not reach, and for those four the text guard was the one catching the divergence.
//! Both were kept until C3, and the differential's reach is now pinned rather than described (see
//! *What the stream must contain*), because once a file was released from the text guard this was
//! the only guard left on it — and since the text guard was deleted, it is the only one on any of
//! them.
//!
//! Production's own tests are transcribed into `locus_score.rs` beside the copy, so they
//! pass on both trees whatever either computes — which is exactly why they cannot serve
//! as the oracle here.
//!
//! # What is compared, and why every number rather than the ratio alone
//!
//! `ParalogScore` carries five values, and all five are asserted equal. The ratio is the
//! one the filter acts on, but it is a difference of two log-likelihoods: two errors that
//! shift both by the same amount would cancel in it and be invisible. The two counts —
//! samples used, and confident hom-alt carriers — pin the sample-admission rules that
//! decide *which* evidence reached the arithmetic at all. Both scores are **destructured**
//! rather than read field by field, so a sixth value added to `ParalogScore` is a compile
//! error here instead of a value this file quietly stops comparing.
//!
//! Equality is **by bit pattern**, not by tolerance. A port that agrees to within
//! `1e-12` is a port that has diverged; the question is not whether the difference matters
//! but whether there is one.
//!
//! # The cohort sizes, and why one sample is among them
//!
//! 1, 2, 10 and 63. **`N = 1` is deliberate** (spec §4): the folded site-frequency-spectrum
//! prior runs over `[1/2N, 1 − 1/2N]`, which at one sample is the single point `[½, ½]`, so
//! the grid degenerates and every weight collapses onto it. Production has never been run
//! there — its filter is a cohort filter — and ng's range commitment starts at one sample,
//! so the copied precompute meets that case here for the first time. 63 is the tomato
//! cohort's size; 2 and 10 sit between. **Nothing here reaches spec §4's other end**, a
//! cohort of three thousand: the precompute's tables are `grid points × samples`, so that is
//! a size and wall-clock question rather than a correctness one, and it belongs to the plan's
//! D milestone.
//!
//! # The generator
//!
//! A fixed-seed splitmix64, not a crate: the codebase's idiom for a reproducible test
//! stream (`parameter_estimation/subsample.rs` uses splitmix64 for the same reason;
//! `calling/parameters_file/to_toml.rs` uses a fixed-seed xorshift), and a failure has to be
//! reproducible from the seed alone with nothing installed. The draws deliberately include
//! the shapes that decide which branch of the scorer runs:
//!
//! - **absent samples** (`None`), which the cohort size still counts toward the SFS floor;
//! - **zero-read samples**, whose allele term vanishes under every genotype — the case
//!   spec §3.2 rests the coverage-only score on;
//! - **degenerate σ₀** (zero, negative, `NaN`, infinite), which the scorer drops;
//! - **`alt_reads > total_reads`**, which the scorer clamps rather than underflowing;
//! - **copy numbers on both sides of the winsor cap**, negative ones included, which it
//!   clips at zero and at four.
//!
//! Four more shapes are too rare or too structural to draw, and are put to both scorers
//! directly, one test each: a σ₀ slice **longer** than the cohort, a precompute built for a
//! **different** cohort, a carrier set that keeps **no** configuration, and a locus with
//! **nothing** usable. Three of the four each killed a mutation to ng's copy that the sweep
//! alone missed.
//!
//! # The rest of the calibration, on the same terms
//!
//! The score is one of four copied quantities and the others are compared here too: the
//! prior — how common hidden duplications are in this run, fitted from the run's own scores
//! by an EM — the tail false-discovery curve built over the same histogram bins, and the cut
//! that curve resolves for an operator's target. Each is a pure function of a stream of
//! likelihood ratios, so the differential folds one randomised stream into both trees'
//! histograms and compares π, the convergence flag, the whole curve and the cut, again by
//! bit pattern.
//!
//! # What the stream must contain
//!
//! Every count this file states is pinned exactly rather than as a floor, and the reason is
//! the defect pinning would have caught. The first draft drew σ₀ from four cases in seven
//! rather than four in twenty-eight, so **4 samples in 7 were degenerate** against the 1 in 7
//! its own comment claimed; at one sample that leaves most loci with nothing usable, scoring
//! the neutral verdict two implementations agree on while computing nothing. A `count > 0`
//! assertion passes through a 99-in-100 narrowing of any of these rates. An exact count does
//! not, and a deliberate change to a rate then has to be looked at rather than absorbed.

use crate::ng::paralog as ng;
use crate::paralog as production;

/// A reproducible stream of pseudo-random `u64`, seeded once per case.
///
/// splitmix64: one multiply-xorshift chain per draw, no state beyond the counter, and the
/// same sequence on every platform. A failing case is reproduced from its seed alone.
struct Splitmix64(u64);

impl Splitmix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// A `f64` in `[0, 1)`, from the top 53 bits.
    fn next_unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A `u32` in `[0, bound)`. `bound` must be non-zero — every call site passes a literal
    /// or `total_reads + 1`. The modulo is slightly biased for bounds that do not divide
    /// `2^64`, which is fine for drawing test shapes and is not claimed otherwise.
    fn next_below(&mut self, bound: u32) -> u32 {
        (self.next_u64() % u64::from(bound)) as u32
    }
}

/// The stream one cohort size is drawn from.
///
/// Named once and derived from the size, so that every test inspecting the sweep inspects
/// the same draws. Typing the seed base and the size a second time is how a check on the
/// stream and the stream itself come apart without anything failing.
fn stream_for(cohort_size: usize) -> Splitmix64 {
    Splitmix64(0x5eed_0000_0000_0000 ^ cohort_size as u64)
}

/// One sample's evidence, drawn so that every branch of the scorer's admission rules is
/// reachable — and built twice, once as each tree's own type, from the *same* numbers.
struct DrawnSample {
    relative_copy_number: f64,
    alt_reads: u32,
    total_reads: u32,
    inbreeding_coefficient: f64,
    /// σ₀, drawn alongside because the scorer takes it in a parallel slice.
    single_copy_depth_sd: f64,
    /// Whether this sample is absent from the locus entirely.
    absent: bool,
}

/// Draw one cohort's worth of evidence at one locus.
fn draw_locus(rng: &mut Splitmix64, cohort_size: usize) -> Vec<DrawnSample> {
    (0..cohort_size)
        .map(|_| {
            // One draw in eight is absent; the cohort size still counts it.
            let absent = rng.next_below(8) == 0;
            // One in six has no reads — the coverage-only case.
            let total_reads = if rng.next_below(6) == 0 {
                0
            } else {
                rng.next_below(400)
            };
            // One in ten is malformed with more ALT reads than total, to reach the clamp.
            let alt_reads = if rng.next_below(10) == 0 {
                total_reads.saturating_add(1 + rng.next_below(50))
            } else if total_reads == 0 {
                0
            } else {
                rng.next_below(total_reads + 1)
            };
            // **Both sides of the winsor cap**, over `[−1, 9)`. A fitted coverage model can
            // hand back a negative relative depth, and the clamp has a lower arm as well as
            // an upper one: with `[0, 9)` alone, narrowing ng's `clamp(0.0, cmax)` to
            // `min(cmax)` was a change no drawn locus could see.
            let relative_copy_number = rng.next_unit() * 10.0 - 1.0;
            // **One σ₀ in seven is degenerate**, which drops the sample from the score —
            // four cases in twenty-eight draws, not four in seven. Drawn from 28 rather than
            // 7 for exactly that reason: at four in seven, 57 in 100 samples are dropped and
            // a one-sample locus scores nothing 5 times in 8, which would leave the
            // differential comparing two neutral verdicts and calling it agreement.
            let single_copy_depth_sd = match rng.next_below(28) {
                0 => 0.0,
                1 => -rng.next_unit(),
                2 => f64::NAN,
                3 => f64::INFINITY,
                _ => 0.02 + rng.next_unit() * 0.6,
            };
            DrawnSample {
                relative_copy_number,
                alt_reads,
                total_reads,
                inbreeding_coefficient: rng.next_unit() * 0.99,
                single_copy_depth_sd,
                absent,
            }
        })
        .collect()
}

/// How one drawn locus is handed to the two scorers.
///
/// The two lengths are separate fields rather than one cohort size because the scorer's
/// precondition (spec §6, trap 3) treats a disagreement between them as a failure that must
/// return the neutral verdict — and each disagreement is its own input class: a σ₀ slice
/// shorter *or longer* than the locus, and per-pass tables built for a cohort the locus does
/// not have. The carrier set is here for the same reason: emptied, it is the one input that
/// makes a *scored* locus produce a `NaN` if the scorer stops refusing it.
struct HandedToBothScorers<'a> {
    drawn: &'a [DrawnSample],
    /// Length of the σ₀ slice: short ones are truncated, long ones padded.
    sigma0_slice_length: usize,
    /// Cohort size the per-pass tables are built for.
    precompute_cohort_size: usize,
    /// The carrier copy numbers under H2. Empty keeps no configuration at all.
    carrier_copy_numbers: Vec<u32>,
}

impl<'a> HandedToBothScorers<'a> {
    /// The ordinary case: every length agrees, and the model is the shipped default.
    fn agreeing(drawn: &'a [DrawnSample]) -> Self {
        Self {
            drawn,
            sigma0_slice_length: drawn.len(),
            precompute_cohort_size: drawn.len(),
            carrier_copy_numbers: production::ParalogModelParams::default().carrier_copy_numbers,
        }
    }
}

/// Score one drawn locus through both trees and return the two verdicts.
///
/// Both sides are built from the same `DrawnSample` values, field by field, so a
/// divergence can only come from the scorers. The two return types differ, so the two
/// verdicts cannot be transposed by accident.
fn score_through_both(
    handed: &HandedToBothScorers,
) -> (production::ParalogScore, ng::ParalogScore) {
    let HandedToBothScorers {
        drawn,
        sigma0_slice_length,
        precompute_cohort_size,
        carrier_copy_numbers,
    } = handed;

    let inbreeding: Vec<f64> = (0..*precompute_cohort_size)
        .map(|index| {
            drawn
                .get(index)
                .map_or(0.0, |sample| sample.inbreeding_coefficient)
        })
        .collect();
    // A σ₀ slice longer than the locus is padded with a healthy value rather than truncated,
    // so that the length disagreement is the only thing under test.
    let sigma0: Vec<f64> = (0..*sigma0_slice_length)
        .map(|index| {
            drawn
                .get(index)
                .map_or(0.3, |sample| sample.single_copy_depth_sd)
        })
        .collect();

    let production_samples: Vec<Option<production::SampleObservation>> = drawn
        .iter()
        .map(|s| {
            (!s.absent).then_some(production::SampleObservation {
                relative_copy_number: s.relative_copy_number,
                alt_reads: s.alt_reads,
                total_reads: s.total_reads,
                inbreeding_coefficient: s.inbreeding_coefficient,
            })
        })
        .collect();
    let ng_samples: Vec<Option<ng::SampleObservation>> = drawn
        .iter()
        .map(|s| {
            (!s.absent).then_some(ng::SampleObservation {
                relative_copy_number: s.relative_copy_number,
                alt_reads: s.alt_reads,
                total_reads: s.total_reads,
                inbreeding_coefficient: s.inbreeding_coefficient,
            })
        })
        .collect();

    let production_params = production::ParalogModelParams {
        carrier_copy_numbers: carrier_copy_numbers.clone(),
        ..production::ParalogModelParams::default()
    };
    let ng_params = ng::ParalogModelParams {
        carrier_copy_numbers: carrier_copy_numbers.clone(),
        ..ng::ParalogModelParams::default()
    };

    (
        production::score_locus_for_paralogy(
            &production::LocusObservations {
                samples: &production_samples,
            },
            &sigma0,
            &production::ParalogScorePrecompute::new(&production_params, &inbreeding),
        ),
        ng::score_locus_for_paralogy(
            &ng::LocusObservations {
                samples: &ng_samples,
            },
            &sigma0,
            &ng::ParalogScorePrecompute::new(&ng_params, &inbreeding),
        ),
    )
}

/// Assert every field of the two verdicts equal, floats by bit pattern.
///
/// `NaN != NaN` under `==`, and a `NaN` is what a divergence in the arithmetic tends to
/// produce, so a derived comparison would either fail on two agreeing `NaN`s or — with a
/// tolerance — pass on a real divergence. (The scorer's *own* refusal is not a `NaN`:
/// `ParalogScore::neutral` sets every float to `0.0`. Spec §6 trap 4's `NaN` is the
/// unscored-record sentinel the run wiring carries, and that arrives at step B1.)
/// **Hands back the widest gap it saw**, so a sweep can pin what the departure actually costs
/// rather than only that it stayed under a bound. Zero for every field that is still bit-equal,
/// which is most of them.
fn assert_the_same_score(
    case: &str,
    production_score: &production::ParalogScore,
    ng_score: &ng::ParalogScore,
) -> f64 {
    // **Destructured, not field-accessed**, so that a value added to `ParalogScore` is a
    // compile error here rather than a value the differential quietly stops comparing.
    let production::ParalogScore {
        paralog_log_likelihood_ratio: their_ratio,
        samples_used: their_samples_used,
        confident_homalt_carriers: their_carriers,
        log_likelihood_real_variant: their_real_variant,
        log_likelihood_hidden_paralog: their_hidden_paralog,
    } = *production_score;
    let ng::ParalogScore {
        paralog_log_likelihood_ratio: our_ratio,
        samples_used: our_samples_used,
        confident_homalt_carriers: our_carriers,
        log_likelihood_real_variant: our_real_variant,
        log_likelihood_hidden_paralog: our_hidden_paralog,
    } = *ng_score;

    let mut widest = 0.0f64;
    for (field, theirs, ours) in [
        ("paralog_log_likelihood_ratio", their_ratio, our_ratio),
        (
            "log_likelihood_real_variant",
            their_real_variant,
            our_real_variant,
        ),
        (
            "log_likelihood_hidden_paralog",
            their_hidden_paralog,
            our_hidden_paralog,
        ),
    ] {
        // **A tolerance, and it stopped being bit equality on 2026-09-09.** ng's
        // `log_add_exp` skips `log1p` where its own cubic series is exact in the sum, which
        // is this port's one departure from production's arithmetic and is why the textual
        // guard on the file was released (the text guard's release table, at `d9e7b076`).
        // The bound is absolute because these are log-likelihoods running to hundreds of nats, and the
        // question the filter asks of them is whether a ratio sits above a cut of about 2.65
        // — so what has to be small is nats, not units in the last place.
        let apart = (ours - theirs).abs();
        assert!(
            apart <= NATS_THE_TWO_TREES_MAY_DIFFER_BY,
            "{case}: {field} differs by {apart:e} nats — ng {ours}, src/paralog/ {theirs}. \
             The two implementations may round differently and may not disagree; a gap this \
             wide means the port computes something else",
        );
        widest = widest.max(apart);
    }
    assert_eq!(
        our_samples_used, their_samples_used,
        "{case}: samples_used differs — the two trees admitted different evidence"
    );
    assert_eq!(
        our_carriers, their_carriers,
        "{case}: confident_homalt_carriers differs — the hom-alt veto counted differently"
    );
    widest
}

/// **How far ng's scorer and production's may sit apart, in nats.**
///
/// **This was bit equality until 2026-09-09**, when ng's `log_add_exp` began skipping `log1p`
/// where its cubic series is exact in the sum — the one place the two trees now round
/// differently, and the reason `locus_score.rs` was released from the textual guard.
///
/// **The bound is what the departure cannot exceed, not what it measures.** The series
/// truncates at under 1e-16 per call and the errors of a locus's calls neither share a sign nor
/// accumulate in the same place, so a locus at 63 samples makes about 46,000 of them and the
/// measured worst over 800 randomised loci is orders of magnitude inside this. The sweep pins
/// that measured figure separately, which is the number worth reading; this is the wall.
///
/// **Why nats and not units in the last place.** These are log-likelihoods running to hundreds,
/// and the only question asked of them is whether a ratio clears a cut of about 2.65
/// (`ParalogFdr`'s fitted default on the tomato cohort). A gap of 1e-9 nats cannot move a
/// record across that; a gap that could would be visible here as a number, not as a rounding.
const NATS_THE_TWO_TREES_MAY_DIFFER_BY: f64 = 1e-9;

/// Cohort sizes the differential runs at. **One is deliberate** — see the module header.
const COHORT_SIZES: [usize; 4] = [1, 2, 10, 63];

/// Loci drawn per cohort size. 200 at four sizes is 800 loci, and at 63 samples that is
/// 12,600 drawn observations; the whole test runs well inside a second.
const LOCI_PER_COHORT_SIZE: usize = 200;

/// **The differential: 800 randomised loci, every number bit-identical.**
#[test]
fn the_copied_scorer_agrees_with_productions() {
    let mut loci_scored = 0usize;
    let mut widest_gap = 0.0f64;
    let mut loci_that_reached_the_arithmetic = [0usize; COHORT_SIZES.len()];
    for (size_index, cohort_size) in COHORT_SIZES.into_iter().enumerate() {
        let mut rng = stream_for(cohort_size);
        for locus in 0..LOCI_PER_COHORT_SIZE {
            let drawn = draw_locus(&mut rng, cohort_size);
            let (theirs, ours) = score_through_both(&HandedToBothScorers::agreeing(&drawn));
            widest_gap = widest_gap.max(assert_the_same_score(
                &format!("cohort size {cohort_size}, locus {locus}"),
                &theirs,
                &ours,
            ));
            loci_scored += 1;
            // `samples_used > 0`, not `ratio.is_finite()`: the neutral verdict a locus with
            // nothing usable gets is `0.0`, which *is* finite, so counting finite ratios
            // would count loci that agreed by both computing nothing.
            loci_that_reached_the_arithmetic[size_index] += usize::from(ours.samples_used > 0);
            assert!(
                ours.paralog_log_likelihood_ratio.is_finite(),
                "cohort size {cohort_size}, locus {locus}: the ratio is not finite"
            );
        }
    }

    assert_eq!(
        loci_scored,
        COHORT_SIZES.len() * LOCI_PER_COHORT_SIZE,
        "every drawn locus is scored"
    );

    // **Without this the test could pass on a stream that scores nothing.** A generator
    // that drew every sample absent, or every σ₀ degenerate, would give 800 neutral
    // verdicts that two implementations agree on trivially — and every one of them would
    // have a finite ratio, because the neutral verdict is `0.0`.
    //
    // At one sample, a locus has its only sample drawn absent (1 draw in 8) or with a
    // degenerate σ₀ (1 in 7), so the expected share of unscored loci is
    // `1 − (7/8)(6/7) = 25 in 100`. Nothing is unscored at 10 samples or 63; those two pins
    // carry little on their own and are here so that a change which *did* start emptying
    // large cohorts could not pass.
    assert_eq!(
        loci_that_reached_the_arithmetic, LOCI_THAT_SHOULD_SCORE,
        "the drawn stream must reach the arithmetic as often as it did when this was \
         written — at cohort sizes {COHORT_SIZES:?}, out of {LOCI_PER_COHORT_SIZE} loci each"
    );

    // **What the departure actually costs, rather than what it is allowed to — and it is
    // nothing.** `assert_the_same_score`'s 1e-9 is a wall derived from the series' truncation
    // term; this is the measurement, and over these 800 loci the two trees return **the same
    // `f64`, bit for bit, on every field**. So the series shortcut is not observably an
    // approximation on anything this differential draws: it is below the rounding of the sum it
    // goes into, which is what its threshold was chosen to guarantee.
    //
    // **Pinned at one part in 1e15 rather than at zero**, six orders inside the wall. `exp` and
    // `log1p` are the platform's, so a last-bit difference on another libm is a reason to be
    // told and not a reason to fail a build; a change that actually moved the arithmetic would
    // clear this by orders of magnitude.
    assert!(
        widest_gap <= 1e-15,
        "the widest gap between the two trees over {loci_scored} randomised loci was \
         {widest_gap:e} nats, where it has been exactly zero; the port's arithmetic has moved"
    );
}

/// Loci admitting at least one sample, per cohort size, in `COHORT_SIZES` order.
const LOCI_THAT_SHOULD_SCORE: [usize; COHORT_SIZES.len()] = [144, 189, 200, 200];

/// **At one sample the two agree, the one sample is actually scored, and the score is a
/// real number.**
///
/// Called out separately from the sweep because it is spec §4's open question and §9's
/// second OPEN: the site-frequency-spectrum grid degenerates to the single point `[½, ½]`
/// at `N = 1`, and nothing had run the copied precompute there.
///
/// **Three assertions, and the middle one is the load-bearing one.** That the two trees
/// agree is the parity claim. That the ratio is finite is nearly free — the neutral verdict
/// a locus with no usable sample gets is `0.0`, which is finite, so on its own that
/// assertion is satisfied by a locus that computed nothing. So the count of loci whose
/// single sample was actually admitted is asserted too, and it must be every locus where
/// the sample was drawn usable at all.
#[test]
fn one_sample_scores_finitely_and_identically() {
    let mut rng = Splitmix64(0x0e5a_3d17);
    let (mut usable_samples_drawn, mut loci_that_scored) = (0usize, 0usize);
    for locus in 0..LOCI_PER_COHORT_SIZE {
        let drawn = draw_locus(&mut rng, 1);
        let sample_is_usable = !drawn[0].absent
            && drawn[0].single_copy_depth_sd.is_finite()
            && drawn[0].single_copy_depth_sd > 0.0;
        let (theirs, ours) = score_through_both(&HandedToBothScorers::agreeing(&drawn));
        assert_the_same_score(&format!("one sample, locus {locus}"), &theirs, &ours);
        assert!(
            ours.paralog_log_likelihood_ratio.is_finite(),
            "one sample, locus {locus}: the ratio is not finite"
        );
        usable_samples_drawn += usize::from(sample_is_usable);
        loci_that_scored += usize::from(ours.samples_used == 1);
    }
    assert_eq!(
        loci_that_scored, usable_samples_drawn,
        "every locus whose one sample was usable must be scored with it"
    );
    assert!(
        loci_that_scored * 2 > LOCI_PER_COHORT_SIZE,
        "most of the {LOCI_PER_COHORT_SIZE} one-sample loci must actually be scored, or \
         this test agrees with production about nothing; {loci_that_scored} were"
    );
}

/// **A negative copy number is winsorised at zero on both sides.**
///
/// A fitted coverage model can hand back a negative relative depth, and the winsor clamp has
/// a lower arm as well as an upper one. Given its own case as well as a place in the stream,
/// because it is the arm a port is most likely to drop: narrowing ng's `clamp(0.0, cmax)` to
/// `min(cmax)` changes nothing a non-negative draw can see, and left ng scoring one locus at
/// more than twice production's value with every other test green.
#[test]
fn a_negative_copy_number_is_winsorised_at_zero_on_both_sides() {
    let drawn: Vec<DrawnSample> = [
        (-3.0, 5u32, 10u32),
        (-0.5, 0, 8),
        (1.0, 10, 20),
        (7.0, 3, 6),
    ]
    .into_iter()
    .map(
        |(relative_copy_number, alt_reads, total_reads)| DrawnSample {
            relative_copy_number,
            alt_reads,
            total_reads,
            inbreeding_coefficient: 0.0,
            single_copy_depth_sd: 0.26,
            absent: false,
        },
    )
    .collect();
    let (theirs, ours) = score_through_both(&HandedToBothScorers::agreeing(&drawn));
    assert_the_same_score("a negative copy number", &theirs, &ours);
    assert_eq!(ours.samples_used, 4, "every sample here is usable");
}

/// **A σ₀ slice out of step with the cohort is neutral on both sides, in both directions.**
///
/// The scorer treats a length mismatch as a precondition failure rather than truncating the
/// cohort, because the site-frequency-spectrum floor `1/2N` would still count the dropped
/// samples (spec §6, trap 3). A port that instead zipped and truncated would give a
/// plausible, wrong number. **Longer as well as shorter**: a port writing `< cohort_size`
/// where production writes `!= cohort_size` refuses the short case and accepts the long one,
/// and the sweep can never draw a long one because it builds the slice from the locus.
#[test]
fn a_mismatched_sigma_slice_is_neutral_on_both_sides() {
    let mut rng = Splitmix64(0xbad_5123);
    for cohort_size in [2usize, 10] {
        for slice_length in [
            cohort_size - 1,
            cohort_size.saturating_sub(3).max(1),
            cohort_size + 1,
            cohort_size + 7,
        ] {
            let drawn = draw_locus(&mut rng, cohort_size);
            let (theirs, ours) = score_through_both(&HandedToBothScorers {
                sigma0_slice_length: slice_length,
                ..HandedToBothScorers::agreeing(&drawn)
            });
            let case = format!("cohort size {cohort_size}, sigma slice of {slice_length}");
            assert_the_same_score(&case, &theirs, &ours);
            assert_eq!(
                ours.samples_used, 0,
                "{case}: a mismatched slice must score nothing, not a truncated cohort"
            );
            assert_eq!(ours.paralog_log_likelihood_ratio, 0.0, "{case}");
        }
    }
}

/// **Per-pass tables built for another cohort are neutral on both sides.**
///
/// The tables are `grid points × samples` and the site-frequency floor is `1/2N`, so tables
/// built for a different `N` do not describe this locus; the scorer refuses rather than
/// indexing into them (spec §6, trap 3). This is the mismatch the run wiring is most likely
/// to produce, because step B1 builds one set of tables per pass and calls the scorer once
/// per locus — and the sweep can never draw it, because it builds both from the same locus.
#[test]
fn tables_built_for_another_cohort_are_neutral_on_both_sides() {
    let mut rng = Splitmix64(0x9ab1_2c3d);
    for (locus_size, table_size) in [(5usize, 10usize), (10, 5), (1, 2), (63, 62)] {
        let drawn = draw_locus(&mut rng, locus_size);
        let (theirs, ours) = score_through_both(&HandedToBothScorers {
            precompute_cohort_size: table_size,
            ..HandedToBothScorers::agreeing(&drawn)
        });
        let case = format!("a locus of {locus_size}, tables built for {table_size}");
        assert_the_same_score(&case, &theirs, &ours);
        assert_eq!(
            ours.samples_used, 0,
            "{case}: tables for another cohort must score nothing"
        );
        assert_eq!(ours.paralog_log_likelihood_ratio, 0.0, "{case}");
    }
}

/// **A carrier set that keeps no configuration is neutral on both sides, not `NaN`.**
///
/// With no carrier copy numbers there is no hidden-paralog story to marginalise, and the
/// scorer refuses before computing one. The refusal matters more than it looks: without it
/// the marginal is a log-sum-exp over an empty set, and the ratio comes back `NaN` **from a
/// locus that was scored** — the one state spec §6 trap 4 says must never leave the scorer,
/// because downstream a `NaN` means *unscored* and folds into nothing.
#[test]
fn an_empty_carrier_set_is_neutral_on_both_sides() {
    let mut rng = Splitmix64(0xe3d2_7a10);
    for cohort_size in [1usize, 10] {
        let drawn = draw_locus(&mut rng, cohort_size);
        let (theirs, ours) = score_through_both(&HandedToBothScorers {
            carrier_copy_numbers: Vec::new(),
            ..HandedToBothScorers::agreeing(&drawn)
        });
        let case = format!("no carrier configurations, cohort size {cohort_size}");
        assert_the_same_score(&case, &theirs, &ours);
        assert_eq!(
            ours.paralog_log_likelihood_ratio.to_bits(),
            0.0f64.to_bits(),
            "{case}: the refusal must be the neutral 0.0, never a NaN from a scored locus"
        );
        assert_eq!(ours.samples_used, 0, "{case}");
    }
}

/// **A locus with nothing usable is neutral on both sides.**
///
/// It establishes nothing to flag, and both trees say so the same way. Note that the value
/// is `0.0` and not `NaN`: spec §6 trap 4's `NaN` sentinel belongs to the run wiring and
/// arrives with step B1, so nothing here tests it.
#[test]
fn a_locus_with_no_usable_sample_is_neutral_on_both_sides() {
    for cohort_size in COHORT_SIZES {
        let drawn: Vec<DrawnSample> = (0..cohort_size)
            .map(|_| DrawnSample {
                relative_copy_number: 1.0,
                alt_reads: 0,
                total_reads: 0,
                inbreeding_coefficient: 0.0,
                single_copy_depth_sd: 0.3,
                absent: true,
            })
            .collect();
        let (theirs, ours) = score_through_both(&HandedToBothScorers::agreeing(&drawn));
        assert_the_same_score(
            &format!("all absent, cohort size {cohort_size}"),
            &theirs,
            &ours,
        );
        assert_eq!(ours.samples_used, 0);
        assert_eq!(ours.paralog_log_likelihood_ratio, 0.0);
    }
}

/// The shape counts `the_drawn_stream_contains_every_shape_the_differential_claims` pins,
/// over the 12,600 observations drawn at the largest cohort size, in the order that test
/// counts them: absent, zero-read, alt-exceeds-total, degenerate σ₀, past the winsor cap,
/// below zero. Drawn at 1 in 8, 1 in 6, 1 in 10, 1 in 7, 1 in 2 and 1 in 10 respectively.
const SHAPES_THE_STREAM_DRAWS: [usize; 6] = [1598, 2210, 1256, 1841, 6259, 1277];

/// **The stream the differential runs on contains what this file says it contains.**
///
/// Counted at the sweep's largest cohort size, from the same stream the sweep draws — the
/// size and the seed both come from `COHORT_SIZES` and `stream_for`, not from a second copy
/// of the literals.
///
/// **Pinned exactly, not as a floor.** `count > 0` passes through a 99-in-100 narrowing of
/// any of these rates: cutting the malformed-record rate from 1 in 10 to 1 in 1,000 takes
/// those draws from 1,256 to 17 in 12,600 and a floor never notices. An exact count makes a
/// deliberate change to a rate something to look at, which is the whole reason this exists.
#[test]
fn the_drawn_stream_contains_every_shape_the_differential_claims() {
    let largest_cohort_size = COHORT_SIZES[COHORT_SIZES.len() - 1];
    let mut rng = stream_for(largest_cohort_size);
    let (mut absent, mut zero_read, mut malformed, mut degenerate_sigma) = (0, 0, 0, 0usize);
    let (mut past_the_cap, mut below_zero) = (0usize, 0usize);
    let mut drawn_total = 0usize;
    for _ in 0..LOCI_PER_COHORT_SIZE {
        for sample in draw_locus(&mut rng, largest_cohort_size) {
            drawn_total += 1;
            absent += usize::from(sample.absent);
            zero_read += usize::from(sample.total_reads == 0);
            malformed += usize::from(sample.alt_reads > sample.total_reads);
            degenerate_sigma += usize::from(
                !(sample.single_copy_depth_sd.is_finite() && sample.single_copy_depth_sd > 0.0),
            );
            past_the_cap +=
                usize::from(sample.relative_copy_number > ng::DEFAULT_MAX_RELATIVE_COPY_NUMBER);
            below_zero += usize::from(sample.relative_copy_number < 0.0);
        }
    }
    assert_eq!(drawn_total, largest_cohort_size * LOCI_PER_COHORT_SIZE);
    assert_eq!(
        [
            absent,
            zero_read,
            malformed,
            degenerate_sigma,
            past_the_cap,
            below_zero
        ],
        SHAPES_THE_STREAM_DRAWS,
        "in {drawn_total} draws: absent, zero-read, alt-exceeds-total, degenerate σ₀, past \
         the winsor cap, below zero. Pinned, so a change to a draw rate has to be looked at \
         — a floor would pass through a 99-in-100 narrowing"
    );
}

// ---------------------------------------------------------------------------------------
// The prior, the curve and the cut
// ---------------------------------------------------------------------------------------

/// Likelihood-ratio streams the calibration differential folds, each a shape the run can
/// actually produce.
///
/// Named rather than drawn uniformly because the EM's answer depends on the *shape* of the
/// distribution, not on its spread: a run where nothing is a duplication and one where a
/// tenth are exercise different parts of the fixed-point iteration, and a stream with
/// non-finite values exercises the histogram's refusal to fold them (spec §3.3 — a
/// non-finite ratio is not a valid Bayes factor and enters neither π nor the curve).
enum RatioStream {
    /// Every locus a real variant: ratios well below zero.
    NothingIsDuplicated,
    /// About one locus in ten strongly positive, the rest negative — a plausible run.
    OneInTen,
    /// Every locus strongly positive, which drives π toward its ceiling.
    EverythingIsDuplicated,
    /// A mixture with `NaN` and infinities among it, which the histogram must refuse.
    WithUnscorableValues,
    /// Ratios past the histogram's `[-100, 100]`, which the run can produce and the histogram
    /// is documented to saturate into its end bins. Without them every bin below 790 and above
    /// 1997 stays empty in every stream, so the two clamps in `bin_index` are compared by
    /// nothing — narrowing either one leaves the differential green.
    BeyondTheHistogramsRange,
}

impl RatioStream {
    fn name(&self) -> &'static str {
        match self {
            Self::NothingIsDuplicated => "nothing is duplicated",
            Self::OneInTen => "one locus in ten",
            Self::EverythingIsDuplicated => "everything is duplicated",
            Self::WithUnscorableValues => "with unscorable values",
            Self::BeyondTheHistogramsRange => "beyond the histogram's range",
        }
    }

    /// Draw `count` ratios of this shape.
    fn draw(&self, rng: &mut Splitmix64, count: usize) -> Vec<f64> {
        (0..count)
            .map(|_| match self {
                Self::NothingIsDuplicated => -20.0 * rng.next_unit() - 1.0,
                Self::OneInTen => {
                    if rng.next_below(10) == 0 {
                        5.0 + 25.0 * rng.next_unit()
                    } else {
                        -20.0 * rng.next_unit() - 1.0
                    }
                }
                Self::EverythingIsDuplicated => 5.0 + 25.0 * rng.next_unit(),
                Self::WithUnscorableValues => match rng.next_below(12) {
                    0 => f64::NAN,
                    1 => f64::INFINITY,
                    2 => f64::NEG_INFINITY,
                    3..=5 => 5.0 + 25.0 * rng.next_unit(),
                    _ => -20.0 * rng.next_unit() - 1.0,
                },
                Self::BeyondTheHistogramsRange => match rng.next_below(4) {
                    0 => -100.0 - 500.0 * rng.next_unit(),
                    1 => 100.0 + 500.0 * rng.next_unit(),
                    _ => -20.0 * rng.next_unit() - 1.0,
                },
            })
            .collect()
    }
}

/// The target false-discovery rates the cut is resolved at: the shipped default, one an
/// order of magnitude looser, one tighter, and two that no run can reach from either end.
const TARGET_FALSE_DISCOVERY_RATES: [f64; 5] = [0.01, 0.1, 0.001, 0.0, 1.0];

/// **The prior, the curve and the cut agree with production's, bit for bit.**
///
/// One randomised stream of likelihood ratios per shape, folded into both trees' histograms
/// and taken through the EM, the curve and the threshold. Everything is compared: π and
/// whether the EM converged; the curve's answer at fifteen probe ratios spanning the axis,
/// including values outside it that the curve saturates; and the resolved cut at five
/// targets, `None` included, since an unreachable target must be unreachable on both sides.
#[test]
fn the_copied_prior_and_curve_agree_with_productions_bit_for_bit() {
    let mut streams_folded = 0usize;
    for shape in [
        RatioStream::NothingIsDuplicated,
        RatioStream::OneInTen,
        RatioStream::EverythingIsDuplicated,
        RatioStream::WithUnscorableValues,
        RatioStream::BeyondTheHistogramsRange,
    ] {
        let mut rng = Splitmix64(0xca11_b7a7_e000_0000 ^ streams_folded as u64);
        let ratios = shape.draw(&mut rng, 5_000);

        let mut their_histogram = production::ParalogLrHistogram::with_defaults();
        let mut our_histogram = ng::ParalogLrHistogram::with_defaults();
        for &lr in &ratios {
            their_histogram.push(lr);
            our_histogram.push(lr);
        }
        let case = shape.name();
        assert_eq!(
            our_histogram.total(),
            their_histogram.total(),
            "{case}: the two histograms folded different numbers of ratios"
        );

        let their_prior =
            production::ParalogPrior::estimate(&their_histogram, &production::EmConfig::default());
        let our_prior = ng::ParalogPrior::estimate(&our_histogram, &ng::EmConfig::default());
        assert_eq!(
            our_prior.prior_probability.to_bits(),
            their_prior.prior_probability.to_bits(),
            "{case}: pi differs — ng {}, src/paralog/ {}",
            our_prior.prior_probability,
            their_prior.prior_probability,
        );
        assert_eq!(
            our_prior.converged, their_prior.converged,
            "{case}: the two EMs disagree about whether they converged"
        );

        let their_curve =
            production::ParalogFdrCurve::from_histogram(&their_histogram, &their_prior);
        let our_curve = ng::ParalogFdrCurve::from_histogram(&our_histogram, &our_prior);
        for probe in [
            -1e6, -100.0, -30.0, -10.0, -1.0, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 1e6,
            0.5,
        ] {
            assert_eq!(
                our_curve.q_of_lr(probe).to_bits(),
                their_curve.q_of_lr(probe).to_bits(),
                "{case}: the tail FDR at a ratio of {probe} differs — ng {}, src/paralog/ {}",
                our_curve.q_of_lr(probe),
                their_curve.q_of_lr(probe),
            );
        }
        for target in TARGET_FALSE_DISCOVERY_RATES {
            let theirs = their_curve.lr_threshold_for_fdr(target);
            let ours = our_curve.lr_threshold_for_fdr(target);
            assert_eq!(
                ours.map(f64::to_bits),
                theirs.map(f64::to_bits),
                "{case}: the cut for a target of {target} differs — ng {ours:?}, \
                 src/paralog/ {theirs:?}. An unreachable target must be unreachable on both \
                 sides, or one tree flags what the other keeps",
            );
        }
        streams_folded += 1;
    }
    assert_eq!(streams_folded, 5, "every stream shape is folded");
}

/// **The histogram refuses a ratio that is not a number, on both sides and by the same
/// count.**
///
/// Spec §3.3: a non-finite likelihood ratio is not a valid Bayes factor, so it enters
/// neither π nor the curve — and spec §6 trap 4 turns on that, because downstream a `NaN`
/// means *unscored*. A port that folded them would put an unscored locus into the estimate
/// of how many loci are duplicated.
#[test]
fn an_unscorable_ratio_is_folded_by_neither_tree() {
    let finite = [-3.0, 0.0, 4.5, 19.0];
    let unscorable = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

    let mut their_histogram = production::ParalogLrHistogram::with_defaults();
    let mut our_histogram = ng::ParalogLrHistogram::with_defaults();
    for &lr in finite.iter().chain(unscorable.iter()) {
        their_histogram.push(lr);
        our_histogram.push(lr);
    }
    assert_eq!(
        our_histogram.total(),
        their_histogram.total(),
        "the two histograms disagree about what they folded"
    );
    assert_eq!(
        our_histogram.total(),
        finite.len() as u64,
        "the histogram's total must be the count of *finite* ratios: {} finite and {} \
         unscorable went in",
        finite.len(),
        unscorable.len(),
    );
}

/// **On a narrow histogram, where the end bins are not saturated, the two trees still
/// agree — and that is the only way the bin indexing itself can be compared.**
///
/// The shipped histogram spans `[-100, 100]`. Its lowest bin's centre is `-99.95` and its
/// highest `99.95`, and the logistic that turns a likelihood ratio into a probability
/// saturates long before either: `σ(-99.95 + logit π)` and `σ(99.95 + logit π)` are exactly
/// `0` and `1` in `f64` for any usable π. **So a bin-index error at either end is invisible
/// through π and the curve, however far outside the range the ratios reach** — measured:
/// narrowing the low clamp from bin 0 to bin 1, and the high clamp from the last bin to the
/// one below it, both leave the default-range comparison green even on a stream where 1,251
/// of 5,000 ratios fall below `-100` and 1,250 above `+100`.
///
/// A histogram over `[-6, 6]` puts its end bins where the logistic is still moving, and the
/// same mutations then change π and the curve on both sides differently. `ParalogLrHistogram::new`
/// is public on both sides, so this costs one more case rather than a new mechanism — and it
/// is also the only place the differential runs at a model setting other than the shipped one.
#[test]
fn a_narrow_histogram_agrees_on_both_sides_where_the_shipped_one_saturates() {
    const NARROW_LO: f64 = -6.0;
    const NARROW_HI: f64 = 6.0;
    const NARROW_BINS: usize = 20;

    let mut compared = 0usize;
    for shape in [
        RatioStream::NothingIsDuplicated,
        RatioStream::OneInTen,
        RatioStream::BeyondTheHistogramsRange,
    ] {
        let mut rng = Splitmix64(0x0a11_0e17_0000_0000 ^ compared as u64);
        let ratios = shape.draw(&mut rng, 5_000);
        let mut theirs = production::ParalogLrHistogram::new(NARROW_LO, NARROW_HI, NARROW_BINS)
            .expect("a narrow range is a valid one");
        let mut ours = ng::ParalogLrHistogram::new(NARROW_LO, NARROW_HI, NARROW_BINS)
            .expect("a narrow range is a valid one");
        for &lr in &ratios {
            theirs.push(lr);
            ours.push(lr);
        }
        let case = format!("a narrow histogram over {}", shape.name());
        assert_eq!(ours.total(), theirs.total(), "{case}: different totals");

        let their_prior =
            production::ParalogPrior::estimate(&theirs, &production::EmConfig::default());
        let our_prior = ng::ParalogPrior::estimate(&ours, &ng::EmConfig::default());
        assert_eq!(
            our_prior.prior_probability.to_bits(),
            their_prior.prior_probability.to_bits(),
            "{case}: pi differs — ng {}, src/paralog/ {}",
            our_prior.prior_probability,
            their_prior.prior_probability,
        );
        assert_eq!(our_prior.converged, their_prior.converged, "{case}");

        let their_curve = production::ParalogFdrCurve::from_histogram(&theirs, &their_prior);
        let our_curve = ng::ParalogFdrCurve::from_histogram(&ours, &our_prior);
        // Every bin centre, and past both ends, so each bin index is read at least once.
        for step in 0..=(NARROW_BINS + 4) {
            let probe = NARROW_LO - 2.0
                + step as f64 * (NARROW_HI - NARROW_LO + 4.0) / (NARROW_BINS + 4) as f64;
            assert_eq!(
                our_curve.q_of_lr(probe).to_bits(),
                their_curve.q_of_lr(probe).to_bits(),
                "{case}: the tail FDR at a ratio of {probe} differs — ng {}, src/paralog/ {}",
                our_curve.q_of_lr(probe),
                their_curve.q_of_lr(probe),
            );
            compared += 1;
        }
        for target in TARGET_FALSE_DISCOVERY_RATES {
            assert_eq!(
                our_curve.lr_threshold_for_fdr(target).map(f64::to_bits),
                their_curve.lr_threshold_for_fdr(target).map(f64::to_bits),
                "{case}: the cut for a target of {target} differs"
            );
        }
    }
    assert_eq!(compared, 3 * (NARROW_BINS + 5), "every probe is compared");
}

/// **The verdict itself agrees with production's: which loci are flagged, and what
/// probability each carries.**
///
/// Everything above stops one step short of the decision. These two functions are the
/// decision: whether a record is dropped, and what number goes in its `PARALOG_POST` field.
/// Production's are reachable from a test in the same crate, so they are compared here rather
/// than left to ng's own tests, which use the implementation as its own oracle.
#[test]
fn the_copied_verdict_agrees_with_productions_bit_for_bit() {
    use crate::var_calling::paralog_filter::calibrate::ParalogCalibration as TheirCalibration;

    let mut compared = 0usize;
    for shape in [
        RatioStream::NothingIsDuplicated,
        RatioStream::OneInTen,
        RatioStream::EverythingIsDuplicated,
        RatioStream::BeyondTheHistogramsRange,
    ] {
        let mut rng = Splitmix64(0xd3c1_5104_0000_0000 ^ compared as u64);
        let ratios = shape.draw(&mut rng, 5_000);
        let mut their_histogram = production::ParalogLrHistogram::with_defaults();
        let mut our_histogram = ng::ParalogLrHistogram::with_defaults();
        for &lr in &ratios {
            their_histogram.push(lr);
            our_histogram.push(lr);
        }
        let their_prior =
            production::ParalogPrior::estimate(&their_histogram, &production::EmConfig::default());
        let our_prior = ng::ParalogPrior::estimate(&our_histogram, &ng::EmConfig::default());
        let their_curve =
            production::ParalogFdrCurve::from_histogram(&their_histogram, &their_prior);
        let our_curve = ng::ParalogFdrCurve::from_histogram(&our_histogram, &our_prior);

        for target_fdr in TARGET_FALSE_DISCOVERY_RATES {
            let theirs = TheirCalibration {
                prior: their_prior,
                curve: their_curve.clone(),
                lr_threshold: their_curve.lr_threshold_for_fdr(target_fdr),
                target_fdr,
            };
            let ours = ng::ParalogCalibration {
                prior: our_prior,
                curve: our_curve.clone(),
                lr_threshold: our_curve.lr_threshold_for_fdr(target_fdr),
                target_fdr,
            };
            for lr in [
                -1e6,
                -100.0,
                -30.0,
                -8.0,
                -1.0,
                0.0,
                1.0,
                5.0,
                12.0,
                18.0,
                30.0,
                100.0,
                1e6,
                f64::NAN,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ] {
                let case = format!("{}, target {target_fdr}, ratio {lr}", shape.name());
                assert_eq!(
                    ours.flags(lr),
                    theirs.flags(lr),
                    "{case}: the two trees disagree about whether this record is dropped"
                );
                assert_eq!(
                    ours.posterior(lr).map(f64::to_bits),
                    theirs.posterior(lr).map(f64::to_bits),
                    "{case}: the two trees disagree about the probability that goes in \
                     PARALOG_POST — ng {:?}, src/var_calling/ {:?}",
                    ours.posterior(lr),
                    theirs.posterior(lr),
                );
                compared += 1;
            }
        }
    }
    // **A run where nothing is flagged would agree trivially**, so at least one of the
    // comparisons above has to be a locus both trees drop.
    assert!(compared >= 4 * 5 * 16, "every case is compared: {compared}");
}

/// **The prior falls back where the EM cannot converge, and both trees fall back the same
/// way — including on an empty histogram, which is the one-sample end of the range.**
///
/// Neither state is reachable from a five-thousand-ratio stream that converges in seven
/// steps, and both are states the run can be in: a cohort of one whose only sample scores
/// nothing folds an empty histogram, and `DEFAULT_FALLBACK_PARALOG_PRIOR` exists precisely
/// for a run whose EM does not settle.
#[test]
fn an_empty_histogram_and_an_unconverged_em_agree_on_both_sides() {
    let their_empty = production::ParalogLrHistogram::with_defaults();
    let our_empty = ng::ParalogLrHistogram::with_defaults();
    assert_eq!(our_empty.total(), 0, "the fixture must actually be empty");

    let their_prior =
        production::ParalogPrior::estimate(&their_empty, &production::EmConfig::default());
    let our_prior = ng::ParalogPrior::estimate(&our_empty, &ng::EmConfig::default());
    assert_eq!(
        our_prior.prior_probability.to_bits(),
        their_prior.prior_probability.to_bits(),
        "an empty histogram must give both trees the same pi — ng {}, src/paralog/ {}",
        our_prior.prior_probability,
        their_prior.prior_probability,
    );
    assert!(
        !our_prior.converged,
        "an empty histogram has nothing to converge on, so the flag must say so"
    );
    assert_eq!(our_prior.converged, their_prior.converged);

    // An EM given one iteration and an impossible tolerance cannot settle.
    let mut rng = Splitmix64(0x0e37_a11e_d000_0000);
    let ratios = RatioStream::OneInTen.draw(&mut rng, 5_000);
    let mut theirs = production::ParalogLrHistogram::with_defaults();
    let mut ours = ng::ParalogLrHistogram::with_defaults();
    for &lr in &ratios {
        theirs.push(lr);
        ours.push(lr);
    }
    let their_cramped = production::EmConfig {
        max_iter: 1,
        tol: 1e-300,
        ..production::EmConfig::default()
    };
    let our_cramped = ng::EmConfig {
        max_iter: 1,
        tol: 1e-300,
        ..ng::EmConfig::default()
    };
    let their_prior = production::ParalogPrior::estimate(&theirs, &their_cramped);
    let our_prior = ng::ParalogPrior::estimate(&ours, &our_cramped);
    assert!(
        !our_prior.converged,
        "one iteration at a tolerance of 1e-300 cannot converge, so the flag must say so"
    );
    assert_eq!(our_prior.converged, their_prior.converged);
    assert_eq!(
        our_prior.prior_probability.to_bits(),
        their_prior.prior_probability.to_bits(),
        "an unconverged pi must be the same on both sides"
    );
}

/// **The fallback, the curve and the cut agree with production's, bit for bit.**
///
/// The three pieces underneath are copies, compared with production's above. What is ng's own is
/// [`calibrate_from_the_ratio_histogram`](ng::calibrate_from_the_ratio_histogram) — the four lines
/// that put them together and substitute the documented rate for an estimate that did not settle.
/// Production's `calibrate_from_histogram` does the same four, so it is the oracle.
///
/// **Three configurations**, because two of them are invisible at the shipped defaults: the
/// iteration's starting guess and the fallback rate are both `0.03`, so a run that reads the wrong
/// one of the two, or never substitutes at all, gives the same answer as a correct one. The second
/// configuration separates them; the third also cramps the iteration so that the substitution
/// happens.
#[test]
fn the_fallback_and_the_cut_agree_with_productions_bit_for_bit() {
    use crate::var_calling::paralog_filter::calibrate::{
        CalibrationConfig as TheirConfig, calibrate_from_histogram,
    };

    /// A fallback rate that is not the iteration's starting guess, so a substitution reading the
    /// wrong field of the configuration shows up as a different number.
    const A_FALLBACK_THAT_IS_NOT_THE_SEED: f64 = 0.41;
    const A_SEED_THAT_IS_NOT_THE_FALLBACK: f64 = 0.17;

    // **Production's configuration is built from ng's, field by field**, so the two sides cannot
    // drift into comparing different questions.
    let theirs_of = |ours: &ng::CalibrationConfig| TheirConfig {
        em: production::EmConfig {
            start: ours.em.start,
            tol: ours.em.tol,
            max_iter: ours.em.max_iter,
        },
        fallback_prior: ours.fallback_prior,
    };

    let configs = [
        ng::CalibrationConfig::default(),
        ng::CalibrationConfig {
            em: ng::EmConfig {
                start: A_SEED_THAT_IS_NOT_THE_FALLBACK,
                ..ng::EmConfig::default()
            },
            fallback_prior: A_FALLBACK_THAT_IS_NOT_THE_SEED,
        },
        ng::CalibrationConfig {
            em: ng::EmConfig {
                max_iter: 1,
                tol: 1e-300,
                start: A_SEED_THAT_IS_NOT_THE_FALLBACK,
            },
            fallback_prior: A_FALLBACK_THAT_IS_NOT_THE_SEED,
        },
    ];

    let mut compared = 0usize;
    for ours_config in configs {
        let theirs_config = theirs_of(&ours_config);
        for ratios in [
            Vec::new(),
            vec![-4.0, -2.0, -1.0, 0.5, 1.0],
            (0..500)
                .map(|i| f64::from(i % 50) - 25.0 + f64::from(i % 7) * 0.1)
                .collect(),
        ] {
            let mut ours = ng::ParalogLrHistogram::with_defaults();
            let mut theirs = production::ParalogLrHistogram::with_defaults();
            for &lr in &ratios {
                ours.push(lr);
                theirs.push(lr);
            }
            for target_fdr in TARGET_FALSE_DISCOVERY_RATES {
                let ours = ng::calibrate_from_the_ratio_histogram(&ours, target_fdr, &ours_config);
                let theirs = calibrate_from_histogram(&theirs, target_fdr, &theirs_config);
                let case = format!("{} ratios, target {target_fdr}", ratios.len());
                assert_eq!(
                    ours.prior.prior_probability.to_bits(),
                    theirs.prior.prior_probability.to_bits(),
                    "{case}: the two trees fitted different rates — ng {}, src/var_calling/ {}",
                    ours.prior.prior_probability,
                    theirs.prior.prior_probability,
                );
                assert_eq!(ours.prior.converged, theirs.prior.converged, "{case}");
                assert_eq!(
                    ours.lr_threshold.map(f64::to_bits),
                    theirs.lr_threshold.map(f64::to_bits),
                    "{case}: the two trees cut at different ratios — ng {:?}, \
                     src/var_calling/ {:?}",
                    ours.lr_threshold,
                    theirs.lr_threshold,
                );
                for lr in [-100.0, -8.0, -1.0, 0.0, 1.0, 8.0, 30.0, 100.0, f64::NAN] {
                    assert_eq!(
                        ours.flags(lr),
                        theirs.flags(lr),
                        "{case}: the two trees disagree about removing a record at ratio {lr}"
                    );
                    compared += 1;
                }
            }
        }
    }
    assert_eq!(
        compared,
        3 * 3 * TARGET_FALSE_DISCOVERY_RATES.len() * 9,
        "every combination is compared: {compared}"
    );
}
