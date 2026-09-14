//! **ng's copied scorer computes what production's did — asserted on randomised loci, against
//! production's answers frozen in a fixture.**
//!
//! A textual guard (`copy_fidelity.rs`, deleted at promotion step C3) said the *text* was
//! production's. This says the *numbers* are, which is the property the port actually claims
//! (`doc/devel/ng/spec/hidden_paralog_filter.md` §7's reuse map: *"the parity oracle is
//! bit-identical ratios against production's on identical synthetic inputs"*). The differential
//! only sees what it draws, so its reach is pinned rather than described (see *What the stream
//! must contain*).
//!
//! Production's own tests are transcribed into `locus_score.rs` beside the copy, so they
//! pass on both trees whatever either computes — which is exactly why they cannot serve
//! as the oracle here.
//!
//! # Production's side is a fixture since promotion step C4
//!
//! Until 2026-09-13 every test here ran production's filter (`src/paralog/`,
//! `src/var_calling/paralog_filter/calibrate.rs`) beside ng's on the same inputs. Production is
//! being deleted, so each value production returned was written once, at commit `d9e7b076`, to
//! [`testdata/production_parity_answers.tsv`](testdata/production_parity_answers.tsv) — 6,494
//! values, keyed by the test section and the case — and the tests now compare ng against that file.
//! The inputs are still generated here, from the same seeds, so **a change to a generator changes
//! the inputs without changing the answers, and the comparison fails wherever the answer depends on
//! the draws**: the fixture is only valid for exactly the draws that wrote it. The three refusal
//! cases — a mismatched σ₀ slice, tables for another cohort, no carrier set — answer the neutral
//! `0.0` whatever is drawn, so they pin only that the refusal happens. Each test also asserts it
//! read every answer its section holds, so a loop that silently shrank would fail rather than
//! compare less.
//!
//! # What is compared, and why every number rather than the ratio alone
//!
//! `ParalogScore` carries five values, and all five are compared. The ratio is the one the
//! filter acts on, but it is a difference of two log-likelihoods: two errors that shift both by
//! the same amount would cancel in it and be invisible. The two counts — samples used, and
//! confident hom-alt carriers — pin the sample-admission rules that decide *which* evidence
//! reached the arithmetic at all. ng's score is **destructured** rather than read field by field,
//! so a sixth value added to `ParalogScore` is a compile error here instead of a value this file
//! quietly stops comparing.
//!
//! The three log-likelihoods are compared within [`NATS_THE_TWO_TREES_MAY_DIFFER_BY`] and
//! the prior and the verdict's probabilities within [`PROBABILITIES_MAY_DIFFER_BY_UNITS`] of each
//! other; the counts, the flags and the cut **by bit pattern**.
//!
//! # The cohort sizes, and why one sample is among them
//!
//! 1, 2, 10 and 63. **`N = 1` is deliberate** (spec §4): the folded site-frequency-spectrum
//! prior runs over `[1/2N, 1 − 1/2N]`, which at one sample is the single point `[½, ½]`, so
//! the grid degenerates and every weight collapses onto it. Production had never been run
//! there — its filter was a cohort filter — and ng's range commitment starts at one sample,
//! so the copied precompute met that case here for the first time. 63 is the tomato
//! cohort's size; 2 and 10 sit between. **Nothing here reaches spec §4's other end**, a
//! cohort of three thousand: the precompute's tables are `grid points × samples`, so that is
//! a size and wall-clock question rather than a correctness one.
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
//! Four more shapes are too rare or too structural to draw, and are put to the scorer
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
//! likelihood ratios, so the differential folds one randomised stream into ng's histogram and
//! compares π, the convergence flag, the curve and the cut against production's: the flag and the
//! cut by bit pattern, π and the curve to within rounding.
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

use std::collections::{HashMap, HashSet};
use std::sync::{Mutex, OnceLock};

use crate::paralog as ng;

// ---------------------------------------------------------------------------------------
// Production's answers
// ---------------------------------------------------------------------------------------

/// The fixture, parsed once: key → value, both borrowed from the embedded file.
fn production_answers() -> &'static HashMap<&'static str, &'static str> {
    static ANSWERS: OnceLock<HashMap<&'static str, &'static str>> = OnceLock::new();
    ANSWERS.get_or_init(|| {
        include_str!("testdata/production_parity_answers.tsv")
            .lines()
            .filter(|line| !line.starts_with('#') && !line.is_empty())
            .map(|line| {
                line.split_once('\t')
                    .unwrap_or_else(|| panic!("a fixture line without a tab: {line:?}"))
            })
            .collect()
    })
}

/// Every key some test has read, so each section can check none of its answers went unread.
fn answers_read() -> &'static Mutex<HashSet<&'static str>> {
    static READ: OnceLock<Mutex<HashSet<&'static str>>> = OnceLock::new();
    READ.get_or_init(|| Mutex::new(HashSet::new()))
}

/// Production's answer under `key`, as written. A missing key is a failure, not a skip: it
/// means the test built a case production was never asked about.
fn answer(key: &str) -> &'static str {
    let (&key, &value) = production_answers().get_key_value(key).unwrap_or_else(|| {
        panic!("production's answers have no `{key}` — the case is not the one frozen")
    });
    // Recovered from a poisoned lock rather than panicking on it: a set of keys stays valid
    // whatever another test was doing when it failed, and one failure should report once.
    answers_read()
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .insert(key);
    value
}

fn their_f64(key: &str) -> f64 {
    let bits = answer(key);
    f64::from_bits(
        u64::from_str_radix(bits, 16)
            .unwrap_or_else(|_| panic!("`{key}`: not a bit pattern: {bits}")),
    )
}

fn their_optional_f64(key: &str) -> Option<f64> {
    match answer(key) {
        "none" => None,
        _ => Some(their_f64(key)),
    }
}

fn their_bool(key: &str) -> bool {
    match answer(key) {
        "true" => true,
        "false" => false,
        other => panic!("`{key}`: not a boolean: {other}"),
    }
}

fn their_count(key: &str) -> u64 {
    let count = answer(key);
    count
        .parse()
        .unwrap_or_else(|_| panic!("`{key}`: not a count: {count}"))
}

/// Assert that every answer under `section/` was read. Called at the end of each test.
#[track_caller]
fn assert_every_answer_was_read(section: &str) {
    let prefix = format!("{section}/");
    // The guard is dropped before the assertion, so a failure here does not poison the lock.
    let unread: Vec<&str> = {
        let read = answers_read()
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        production_answers()
            .keys()
            .filter(|key| key.starts_with(&prefix) && !read.contains(*key))
            .copied()
            .collect()
    };
    assert!(
        unread.is_empty(),
        "{} of production's answers in `{section}` were never compared, e.g. {:?} — the test \
         now asks fewer questions than the fixture was written for",
        unread.len(),
        &unread[..unread.len().min(3)],
    );
}

/// Production's `ParalogScore` for one case, read back from the fixture.
struct TheirScore {
    paralog_log_likelihood_ratio: f64,
    samples_used: u64,
    confident_homalt_carriers: u64,
    log_likelihood_real_variant: f64,
    log_likelihood_hidden_paralog: f64,
}

fn their_score(key: &str) -> TheirScore {
    TheirScore {
        paralog_log_likelihood_ratio: their_f64(&format!("{key}/ratio")),
        samples_used: their_count(&format!("{key}/samples_used")),
        confident_homalt_carriers: their_count(&format!("{key}/confident_homalt_carriers")),
        log_likelihood_real_variant: their_f64(&format!("{key}/real_variant")),
        log_likelihood_hidden_paralog: their_f64(&format!("{key}/hidden_paralog")),
    }
}

// ---------------------------------------------------------------------------------------
// The generator
// ---------------------------------------------------------------------------------------

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
/// reachable.
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

/// How one drawn locus is handed to the scorer.
///
/// The two lengths are separate fields rather than one cohort size because the scorer's
/// precondition (spec §6, trap 3) treats a disagreement between them as a failure that must
/// return the neutral verdict — and each disagreement is its own input class: a σ₀ slice
/// shorter *or longer* than the locus, and per-pass tables built for a cohort the locus does
/// not have. The carrier set is here for the same reason: emptied, it is the one input that
/// makes a *scored* locus produce a `NaN` if the scorer stops refusing it.
struct HandedToTheScorer<'a> {
    drawn: &'a [DrawnSample],
    /// Length of the σ₀ slice: short ones are truncated, long ones padded.
    sigma0_slice_length: usize,
    /// Cohort size the per-pass tables are built for.
    precompute_cohort_size: usize,
    /// The carrier copy numbers under H2. Empty keeps no configuration at all.
    carrier_copy_numbers: Vec<u32>,
}

impl<'a> HandedToTheScorer<'a> {
    /// The ordinary case: every length agrees, and the model is the shipped default.
    fn agreeing(drawn: &'a [DrawnSample]) -> Self {
        Self {
            drawn,
            sigma0_slice_length: drawn.len(),
            precompute_cohort_size: drawn.len(),
            carrier_copy_numbers: ng::ParalogModelParams::default().carrier_copy_numbers,
        }
    }
}

/// Score one drawn locus through ng's copy.
///
/// Built exactly as production's side was built when its answers were written: the same
/// inbreeding vector, the same σ₀ padding, the same model parameters but for the carrier set.
fn score_with_ng(handed: &HandedToTheScorer) -> ng::ParalogScore {
    let HandedToTheScorer {
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

    let samples: Vec<Option<ng::SampleObservation>> = drawn
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
    let params = ng::ParalogModelParams {
        carrier_copy_numbers: carrier_copy_numbers.clone(),
        ..ng::ParalogModelParams::default()
    };

    ng::score_locus_for_paralogy(
        &ng::LocusObservations { samples: &samples },
        &sigma0,
        &ng::ParalogScorePrecompute::new(&params, &inbreeding),
    )
}

/// Assert ng's verdict matches production's, every field.
///
/// The three log-likelihoods within [`NATS_THE_TWO_TREES_MAY_DIFFER_BY`]; the two counts
/// exactly. **Hands back the widest gap it saw**, so a sweep can pin what the departure actually
/// costs rather than only that it stayed under a bound.
fn assert_the_same_score(case: &str, theirs: &TheirScore, ng_score: &ng::ParalogScore) -> f64 {
    // **Destructured, not field-accessed**, so that a value added to `ParalogScore` is a
    // compile error here rather than a value the differential quietly stops comparing.
    let ng::ParalogScore {
        paralog_log_likelihood_ratio: our_ratio,
        samples_used: our_samples_used,
        confident_homalt_carriers: our_carriers,
        log_likelihood_real_variant: our_real_variant,
        log_likelihood_hidden_paralog: our_hidden_paralog,
    } = *ng_score;

    let mut widest = 0.0f64;
    for (field, theirs, ours) in [
        (
            "paralog_log_likelihood_ratio",
            theirs.paralog_log_likelihood_ratio,
            our_ratio,
        ),
        (
            "log_likelihood_real_variant",
            theirs.log_likelihood_real_variant,
            our_real_variant,
        ),
        (
            "log_likelihood_hidden_paralog",
            theirs.log_likelihood_hidden_paralog,
            our_hidden_paralog,
        ),
    ] {
        // **A tolerance, and it stopped being bit equality on 2026-09-09.** ng's
        // `log_add_exp` skips `log1p` where its own cubic series is exact in the sum, and since
        // 2026-09-14 ng's `exp`, `ln` and `log1p` are libm's where production's were glibc's.
        let apart = (ours - theirs).abs();
        assert!(
            apart <= NATS_THE_TWO_TREES_MAY_DIFFER_BY,
            "{case}: {field} differs by {apart:e} nats — ng {ours}, production {theirs}. \
             The two implementations may round differently and may not disagree; a gap this \
             wide means the port computes something else",
        );
        widest = widest.max(apart);
    }
    assert_eq!(
        our_samples_used as u64, theirs.samples_used,
        "{case}: samples_used differs — ng admitted different evidence"
    );
    assert_eq!(
        our_carriers as u64, theirs.confident_homalt_carriers,
        "{case}: confident_homalt_carriers differs — the hom-alt veto counted differently"
    );
    widest
}

/// **How far ng's scorer and production's may sit apart, in nats.**
///
/// **This was bit equality until 2026-09-09**, when ng's `log_add_exp` began skipping `log1p`
/// where its cubic series is exact in the sum. Since 2026-09-14 the two trees also round
/// differently because ng's `exp`, `ln` and `log1p` go through [`crate::float`] (libm).
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

/// **The differential: 800 randomised loci, every number production's.**
#[test]
fn the_copied_scorer_agrees_with_productions() {
    let mut loci_scored = 0usize;
    let mut widest_gap = 0.0f64;
    let mut loci_that_reached_the_arithmetic = [0usize; COHORT_SIZES.len()];
    for (size_index, cohort_size) in COHORT_SIZES.into_iter().enumerate() {
        let mut rng = stream_for(cohort_size);
        for locus in 0..LOCI_PER_COHORT_SIZE {
            let drawn = draw_locus(&mut rng, cohort_size);
            let ours = score_with_ng(&HandedToTheScorer::agreeing(&drawn));
            widest_gap = widest_gap.max(assert_the_same_score(
                &format!("cohort size {cohort_size}, locus {locus}"),
                &their_score(&format!("sweep/{cohort_size}/{locus}")),
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
    assert_every_answer_was_read("sweep");

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

    // **What the departure actually costs, rather than what it is allowed to.**
    // `assert_the_same_score`'s 1e-9 is a wall derived from the series' truncation term; this is
    // the measurement. Until the maths moved to libm, ng returned **production's `f64`, bit for
    // bit, on every field** of these 800 loci, so the series shortcut is not observably an
    // approximation on anything this differential draws: it is below the rounding of the sum it
    // goes into, which is what its threshold was chosen to guarantee.
    //
    // **Pinned at 1e-10 nats, ten times inside the wall.** Until step B4 of the portable-float
    // plan this gap was exactly zero. ng's `exp`, `ln` and `log1p` now go through
    // [`crate::float`] (libm), production's answers were written with glibc's, and over these
    // loci the widest field moved by 3.6e-12 nats (measured 2026-09-14, the same on macOS and
    // Linux). A change that actually moved the arithmetic would clear 1e-10 by orders of
    // magnitude.
    assert!(
        widest_gap <= 1e-10,
        "the widest gap between ng and production over {loci_scored} randomised loci was \
         {widest_gap:e} nats, where it has been 3.6e-12 since the maths moved to libm; the \
         port's arithmetic has moved"
    );
}

/// Loci admitting at least one sample, per cohort size, in `COHORT_SIZES` order.
const LOCI_THAT_SHOULD_SCORE: [usize; COHORT_SIZES.len()] = [144, 189, 200, 200];

/// **At one sample ng agrees with production, the one sample is actually scored, and the score
/// is a real number.**
///
/// Called out separately from the sweep because it is spec §4's open question and §9's
/// second OPEN: the site-frequency-spectrum grid degenerates to the single point `[½, ½]`
/// at `N = 1`, and nothing had run the copied precompute there.
///
/// **Three assertions, and the middle one is the load-bearing one.** That the two agree is the
/// parity claim. That the ratio is finite is nearly free — the neutral verdict a locus with no
/// usable sample gets is `0.0`, which is finite, so on its own that assertion is satisfied by a
/// locus that computed nothing. So the count of loci whose single sample was actually admitted
/// is asserted too, and it must be every locus where the sample was drawn usable at all.
#[test]
fn one_sample_scores_finitely_and_identically() {
    let mut rng = Splitmix64(0x0e5a_3d17);
    let (mut usable_samples_drawn, mut loci_that_scored) = (0usize, 0usize);
    for locus in 0..LOCI_PER_COHORT_SIZE {
        let drawn = draw_locus(&mut rng, 1);
        let sample_is_usable = !drawn[0].absent
            && drawn[0].single_copy_depth_sd.is_finite()
            && drawn[0].single_copy_depth_sd > 0.0;
        let ours = score_with_ng(&HandedToTheScorer::agreeing(&drawn));
        assert_the_same_score(
            &format!("one sample, locus {locus}"),
            &their_score(&format!("one_sample/{locus}")),
            &ours,
        );
        assert!(
            ours.paralog_log_likelihood_ratio.is_finite(),
            "one sample, locus {locus}: the ratio is not finite"
        );
        usable_samples_drawn += usize::from(sample_is_usable);
        loci_that_scored += usize::from(ours.samples_used == 1);
    }
    assert_every_answer_was_read("one_sample");
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

/// **A negative copy number is winsorised at zero, as production did.**
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
    let ours = score_with_ng(&HandedToTheScorer::agreeing(&drawn));
    assert_the_same_score(
        "a negative copy number",
        &their_score("negative_copy_number/locus"),
        &ours,
    );
    assert_every_answer_was_read("negative_copy_number");
    assert_eq!(ours.samples_used, 4, "every sample here is usable");
}

/// **A σ₀ slice out of step with the cohort is neutral, in both directions, as production's
/// was.**
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
        // Keyed by position as well as length: at two samples the first two lengths are both 1.
        for (slice_index, slice_length) in [
            cohort_size - 1,
            cohort_size.saturating_sub(3).max(1),
            cohort_size + 1,
            cohort_size + 7,
        ]
        .into_iter()
        .enumerate()
        {
            let drawn = draw_locus(&mut rng, cohort_size);
            let ours = score_with_ng(&HandedToTheScorer {
                sigma0_slice_length: slice_length,
                ..HandedToTheScorer::agreeing(&drawn)
            });
            let case = format!("cohort size {cohort_size}, sigma slice of {slice_length}");
            assert_the_same_score(
                &case,
                &their_score(&format!(
                    "mismatched_sigma/{cohort_size}/{slice_index}/{slice_length}"
                )),
                &ours,
            );
            assert_eq!(
                ours.samples_used, 0,
                "{case}: a mismatched slice must score nothing, not a truncated cohort"
            );
            assert_eq!(ours.paralog_log_likelihood_ratio, 0.0, "{case}");
        }
    }
    assert_every_answer_was_read("mismatched_sigma");
}

/// **Per-pass tables built for another cohort are neutral, as production's were.**
///
/// The tables are `grid points × samples` and the site-frequency floor is `1/2N`, so tables
/// built for a different `N` do not describe this locus; the scorer refuses rather than
/// indexing into them (spec §6, trap 3). This is the mismatch the run wiring is most likely
/// to produce, because it builds one set of tables per pass and calls the scorer once per
/// locus — and the sweep can never draw it, because it builds both from the same locus.
#[test]
fn tables_built_for_another_cohort_are_neutral_on_both_sides() {
    let mut rng = Splitmix64(0x9ab1_2c3d);
    for (locus_size, table_size) in [(5usize, 10usize), (10, 5), (1, 2), (63, 62)] {
        let drawn = draw_locus(&mut rng, locus_size);
        let ours = score_with_ng(&HandedToTheScorer {
            precompute_cohort_size: table_size,
            ..HandedToTheScorer::agreeing(&drawn)
        });
        let case = format!("a locus of {locus_size}, tables built for {table_size}");
        assert_the_same_score(
            &case,
            &their_score(&format!("other_cohort_tables/{locus_size}/{table_size}")),
            &ours,
        );
        assert_eq!(
            ours.samples_used, 0,
            "{case}: tables for another cohort must score nothing"
        );
        assert_eq!(ours.paralog_log_likelihood_ratio, 0.0, "{case}");
    }
    assert_every_answer_was_read("other_cohort_tables");
}

/// **A carrier set that keeps no configuration is neutral, not `NaN`, as production's was.**
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
        let ours = score_with_ng(&HandedToTheScorer {
            carrier_copy_numbers: Vec::new(),
            ..HandedToTheScorer::agreeing(&drawn)
        });
        let case = format!("no carrier configurations, cohort size {cohort_size}");
        assert_the_same_score(
            &case,
            &their_score(&format!("empty_carrier_set/{cohort_size}")),
            &ours,
        );
        assert_eq!(
            ours.paralog_log_likelihood_ratio.to_bits(),
            0.0f64.to_bits(),
            "{case}: the refusal must be the neutral 0.0, never a NaN from a scored locus"
        );
        assert_eq!(ours.samples_used, 0, "{case}");
    }
    assert_every_answer_was_read("empty_carrier_set");
}

/// **A locus with nothing usable is neutral, as production's was.**
///
/// It establishes nothing to flag. Note that the value is `0.0` and not `NaN`: spec §6 trap 4's
/// `NaN` sentinel belongs to the run wiring, so nothing here tests it.
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
        let ours = score_with_ng(&HandedToTheScorer::agreeing(&drawn));
        assert_the_same_score(
            &format!("all absent, cohort size {cohort_size}"),
            &their_score(&format!("no_usable_sample/{cohort_size}")),
            &ours,
        );
        assert_eq!(ours.samples_used, 0);
        assert_eq!(ours.paralog_log_likelihood_ratio, 0.0);
    }
    assert_every_answer_was_read("no_usable_sample");
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

/// Fold `ratios` into a fresh histogram of ng's, at the shipped range.
fn our_histogram_of(ratios: &[f64]) -> ng::ParalogLrHistogram {
    let mut histogram = ng::ParalogLrHistogram::with_defaults();
    for &lr in ratios {
        histogram.push(lr);
    }
    histogram
}

/// **Probabilities agree to within rounding, not bit for bit, since step B4 of the portable-float
/// plan.** Production's answers were written with glibc's `exp` and `ln`; ng's go through
/// [`crate::float`], the libm crate, which rounds some last places the other way, and the EM carries
/// such a unit forward. Every *decision* — which record is dropped, whether the EM converged, how
/// many ratios were folded, the cut itself — is still compared exactly.
///
/// **How far they moved, measured 2026-09-14** by counting, in these tests, every probability whose
/// bits differ from production's: 17 of the 465 compared, identically on macOS and Linux — 10 by one
/// representable number, 5 by two, 1 by eight, and π by twenty-two (a relative 3.6e-15). The
/// allowance of 64 leaves room for that and still refuses any real change to the arithmetic, which
/// moves a posterior by a relative 1e-3 or more.
const PROBABILITIES_MAY_DIFFER_BY_UNITS: u64 = 64;

/// Both `NaN`; or the same infinity; or finite, of one sign, and at most
/// [`PROBABILITIES_MAY_DIFFER_BY_UNITS`] representable numbers apart (`+0.0` and `−0.0` are equal).
/// Counting representable numbers rather than a relative tolerance keeps the check as tight next to
/// 1 (π can be `1 − 2e-13`) as next to 1e-56, and an infinity can match only an infinity.
fn within_rounding(ours: f64, theirs: f64) -> bool {
    if ours.is_nan() || theirs.is_nan() {
        return ours.is_nan() && theirs.is_nan();
    }
    if !ours.is_finite() || !theirs.is_finite() {
        return ours == theirs;
    }
    ours == theirs
        || (ours.is_sign_negative() == theirs.is_sign_negative()
            && ours.to_bits().abs_diff(theirs.to_bits()) <= PROBABILITIES_MAY_DIFFER_BY_UNITS)
}

/// The comparison refuses what a port defect would produce: an infinity against a finite value, a
/// sign change, a number past the allowance, and a value against no value.
#[test]
fn within_rounding_refuses_infinities_sign_changes_and_distant_values() {
    let half = 0.5_f64;
    let steps_away = |value: f64, steps: u64| f64::from_bits(value.to_bits() + steps);
    assert!(within_rounding(half, half));
    assert!(within_rounding(
        half,
        steps_away(half, PROBABILITIES_MAY_DIFFER_BY_UNITS)
    ));
    assert!(!within_rounding(
        half,
        steps_away(half, PROBABILITIES_MAY_DIFFER_BY_UNITS + 1)
    ));
    assert!(within_rounding(0.0, -0.0));
    assert!(!within_rounding(1e-300, -1e-300));
    assert!(within_rounding(f64::NAN, f64::NAN));
    assert!(!within_rounding(f64::NAN, half));
    assert!(within_rounding(f64::INFINITY, f64::INFINITY));
    assert!(!within_rounding(f64::INFINITY, half));
    assert!(!within_rounding(f64::INFINITY, f64::MAX));
    assert!(!within_rounding(f64::NEG_INFINITY, f64::INFINITY));
    assert!(within_rounding(1.0 - 2e-13, steps_away(1.0 - 2e-13, 1)));
    assert!(!within_rounding(1.0 - 2e-13, 1.0));
    assert!(optional_within_rounding(None, None));
    assert!(!optional_within_rounding(Some(half), None));
}

/// [`within_rounding`] for an optional value: both absent, or both present and within rounding.
fn optional_within_rounding(ours: Option<f64>, theirs: Option<f64>) -> bool {
    match (ours, theirs) {
        (None, None) => true,
        (Some(ours), Some(theirs)) => within_rounding(ours, theirs),
        _ => false,
    }
}

/// Assert ng's fitted prior matches production's answer under `key`: π to within rounding, and
/// whether the EM converged exactly.
#[track_caller]
fn assert_the_same_prior(case: &str, key: &str, ours: &ng::ParalogPrior) {
    let theirs = their_f64(&format!("{key}/prior_probability"));
    assert!(
        within_rounding(ours.prior_probability, theirs),
        "{case}: pi differs — ng {}, production {theirs}",
        ours.prior_probability,
    );
    assert_eq!(
        ours.converged,
        their_bool(&format!("{key}/converged")),
        "{case}: ng and production disagree about whether the EM converged"
    );
}

/// **The prior, the curve and the cut agree with production's, to within rounding.**
///
/// One randomised stream of likelihood ratios per shape, folded into ng's histogram and taken
/// through the EM, the curve and the threshold. Everything is compared: the number of ratios
/// folded; π and whether the EM converged; the curve's answer at fifteen probe ratios spanning
/// the axis, including values outside it that the curve saturates; and the resolved cut at five
/// targets, `None` included, since an unreachable target must be unreachable for both.
#[test]
fn the_copied_prior_and_curve_agree_with_productions_to_within_rounding() {
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
        let our_histogram = our_histogram_of(&ratios);
        let case = shape.name();
        let key = format!("prior_and_curve/{case}");
        assert_eq!(
            our_histogram.total(),
            their_count(&format!("{key}/total")),
            "{case}: ng folded a different number of ratios"
        );

        let our_prior = ng::ParalogPrior::estimate(&our_histogram, &ng::EmConfig::default());
        assert_the_same_prior(case, &key, &our_prior);

        let our_curve = ng::ParalogFdrCurve::from_histogram(&our_histogram, &our_prior);
        for probe in [
            -1e6, -100.0, -30.0, -10.0, -1.0, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 1e6,
            0.5,
        ] {
            let theirs = their_f64(&format!("{key}/q_of_lr/{probe:?}"));
            assert!(
                within_rounding(our_curve.q_of_lr(probe), theirs),
                "{case}: the tail FDR at a ratio of {probe} differs — ng {}, production {theirs}",
                our_curve.q_of_lr(probe),
            );
        }
        for target in TARGET_FALSE_DISCOVERY_RATES {
            let theirs = their_optional_f64(&format!("{key}/lr_threshold_for_fdr/{target:?}"));
            let ours = our_curve.lr_threshold_for_fdr(target);
            assert_eq!(
                ours.map(f64::to_bits),
                theirs.map(f64::to_bits),
                "{case}: the cut for a target of {target} differs — ng {ours:?}, \
                 production {theirs:?}. An unreachable target must be unreachable for both, \
                 or ng flags what production kept",
            );
        }
        streams_folded += 1;
    }
    assert_eq!(streams_folded, 5, "every stream shape is folded");
    assert_every_answer_was_read("prior_and_curve");
}

/// **The histogram refuses a ratio that is not a number, by production's count.**
///
/// Spec §3.3: a non-finite likelihood ratio is not a valid Bayes factor, so it enters
/// neither π nor the curve — and spec §6 trap 4 turns on that, because downstream a `NaN`
/// means *unscored*. A port that folded them would put an unscored locus into the estimate
/// of how many loci are duplicated.
#[test]
fn an_unscorable_ratio_is_folded_by_neither_tree() {
    let finite = [-3.0, 0.0, 4.5, 19.0];
    let unscorable = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

    let our_histogram = our_histogram_of(&[&finite[..], &unscorable[..]].concat());
    assert_eq!(
        our_histogram.total(),
        their_count("unscorable_ratio/total"),
        "ng and production disagree about what the histogram folded"
    );
    assert_eq!(
        our_histogram.total(),
        finite.len() as u64,
        "the histogram's total must be the count of *finite* ratios: {} finite and {} \
         unscorable went in",
        finite.len(),
        unscorable.len(),
    );
    assert_every_answer_was_read("unscorable_ratio");
}

/// **On a narrow histogram, where the end bins are not saturated, ng still agrees with
/// production — and that is the only way the bin indexing itself can be compared.**
///
/// The shipped histogram spans `[-100, 100]`. Its lowest bin's centre is `-99.95` and its
/// highest `99.95`, and the logistic that turns a likelihood ratio into a probability
/// saturates long before either: `σ(-99.95 + logit π)` and `σ(99.95 + logit π)` are exactly
/// `0` and `1` in `f64` for any usable π. **So a bin-index error at either end is invisible
/// through π and the curve, however far outside the range the ratios reach** — measured:
/// narrowing the low clamp from bin 0 to bin 1, and the high clamp from the last bin to the
/// one below it, both left the default-range comparison green even on a stream where 1,251
/// of 5,000 ratios fall below `-100` and 1,250 above `+100`.
///
/// A histogram over `[-6, 6]` puts its end bins where the logistic is still moving, and the
/// same mutations then change π and the curve. `ParalogLrHistogram::new` is public, so this
/// costs one more case rather than a new mechanism — and it is also the only place the
/// differential runs at a model setting other than the shipped one.
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
        let mut ours = ng::ParalogLrHistogram::new(NARROW_LO, NARROW_HI, NARROW_BINS)
            .expect("a narrow range is a valid one");
        for &lr in &ratios {
            ours.push(lr);
        }
        let case = format!("a narrow histogram over {}", shape.name());
        let key = format!("narrow_histogram/{}", shape.name());
        assert_eq!(
            ours.total(),
            their_count(&format!("{key}/total")),
            "{case}: different totals"
        );

        let our_prior = ng::ParalogPrior::estimate(&ours, &ng::EmConfig::default());
        assert_the_same_prior(&case, &key, &our_prior);

        let our_curve = ng::ParalogFdrCurve::from_histogram(&ours, &our_prior);
        // Every bin centre, and past both ends, so each bin index is read at least once.
        for step in 0..=(NARROW_BINS + 4) {
            let probe = NARROW_LO - 2.0
                + step as f64 * (NARROW_HI - NARROW_LO + 4.0) / (NARROW_BINS + 4) as f64;
            let theirs = their_f64(&format!("{key}/q_of_lr/{probe:?}"));
            assert!(
                within_rounding(our_curve.q_of_lr(probe), theirs),
                "{case}: the tail FDR at a ratio of {probe} differs — ng {}, production {theirs}",
                our_curve.q_of_lr(probe),
            );
            compared += 1;
        }
        for target in TARGET_FALSE_DISCOVERY_RATES {
            assert_eq!(
                our_curve.lr_threshold_for_fdr(target).map(f64::to_bits),
                their_optional_f64(&format!("{key}/lr_threshold_for_fdr/{target:?}"))
                    .map(f64::to_bits),
                "{case}: the cut for a target of {target} differs"
            );
        }
    }
    assert_eq!(compared, 3 * (NARROW_BINS + 5), "every probe is compared");
    assert_every_answer_was_read("narrow_histogram");
}

/// **The verdict itself agrees with production's: which loci are flagged, and what
/// probability each carries.**
///
/// Everything above stops one step short of the decision. These two functions are the
/// decision: whether a record is dropped, and what number goes in its `PARALOG_POST` field.
#[test]
fn the_copied_verdict_agrees_with_productions_to_within_rounding() {
    let mut compared = 0usize;
    for shape in [
        RatioStream::NothingIsDuplicated,
        RatioStream::OneInTen,
        RatioStream::EverythingIsDuplicated,
        RatioStream::BeyondTheHistogramsRange,
    ] {
        let mut rng = Splitmix64(0xd3c1_5104_0000_0000 ^ compared as u64);
        let ratios = shape.draw(&mut rng, 5_000);
        let our_histogram = our_histogram_of(&ratios);
        let our_prior = ng::ParalogPrior::estimate(&our_histogram, &ng::EmConfig::default());
        let our_curve = ng::ParalogFdrCurve::from_histogram(&our_histogram, &our_prior);

        for target_fdr in TARGET_FALSE_DISCOVERY_RATES {
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
                let key = format!("verdict/{}/{target_fdr:?}/{lr:?}", shape.name());
                assert_eq!(
                    ours.flags(lr),
                    their_bool(&format!("{key}/flags")),
                    "{case}: ng and production disagree about whether this record is dropped"
                );
                let theirs = their_optional_f64(&format!("{key}/posterior"));
                assert!(
                    optional_within_rounding(ours.posterior(lr), theirs),
                    "{case}: ng and production disagree about the probability that goes in \
                     PARALOG_POST — ng {:?}, production {theirs:?}",
                    ours.posterior(lr),
                );
                compared += 1;
            }
        }
    }
    assert!(compared >= 4 * 5 * 16, "every case is compared: {compared}");
    assert_every_answer_was_read("verdict");
}

/// **The prior falls back where the EM cannot converge, as production's did — including on an
/// empty histogram, which is the one-sample end of the range.**
///
/// Neither state is reachable from a five-thousand-ratio stream that converges in seven
/// steps, and both are states the run can be in: a cohort of one whose only sample scores
/// nothing folds an empty histogram, and `DEFAULT_FALLBACK_PARALOG_PRIOR` exists precisely
/// for a run whose EM does not settle.
#[test]
fn an_empty_histogram_and_an_unconverged_em_agree_on_both_sides() {
    let our_empty = ng::ParalogLrHistogram::with_defaults();
    assert_eq!(our_empty.total(), 0, "the fixture must actually be empty");

    let our_prior = ng::ParalogPrior::estimate(&our_empty, &ng::EmConfig::default());
    assert_the_same_prior(
        "an empty histogram",
        "empty_and_unconverged/empty",
        &our_prior,
    );
    assert!(
        !our_prior.converged,
        "an empty histogram has nothing to converge on, so the flag must say so"
    );

    // An EM given one iteration and an impossible tolerance cannot settle.
    let mut rng = Splitmix64(0x0e37_a11e_d000_0000);
    let ratios = RatioStream::OneInTen.draw(&mut rng, 5_000);
    let ours = our_histogram_of(&ratios);
    let our_cramped = ng::EmConfig {
        max_iter: 1,
        tol: 1e-300,
        ..ng::EmConfig::default()
    };
    let our_prior = ng::ParalogPrior::estimate(&ours, &our_cramped);
    assert!(
        !our_prior.converged,
        "one iteration at a tolerance of 1e-300 cannot converge, so the flag must say so"
    );
    assert_the_same_prior(
        "an unconverged EM",
        "empty_and_unconverged/cramped",
        &our_prior,
    );
    assert_every_answer_was_read("empty_and_unconverged");
}

/// **The fallback, the curve and the cut agree with production's, to within rounding.**
///
/// The three pieces underneath are copies, compared with production's above. What is ng's own is
/// [`calibrate_from_the_ratio_histogram`](ng::calibrate_from_the_ratio_histogram) — the four lines
/// that put them together and substitute the documented rate for an estimate that did not settle.
/// Production's `calibrate_from_histogram` did the same four, so its answers are the oracle.
///
/// **Three configurations**, because two of them are invisible at the shipped defaults: the
/// iteration's starting guess and the fallback rate are both `0.03`, so a run that reads the wrong
/// one of the two, or never substitutes at all, gives the same answer as a correct one. The second
/// configuration separates them; the third also cramps the iteration so that the substitution
/// happens. Production was handed the same three, field for field.
#[test]
fn the_fallback_and_the_cut_agree_with_productions_to_within_rounding() {
    /// A fallback rate that is not the iteration's starting guess, so a substitution reading the
    /// wrong field of the configuration shows up as a different number.
    const A_FALLBACK_THAT_IS_NOT_THE_SEED: f64 = 0.41;
    const A_SEED_THAT_IS_NOT_THE_FALLBACK: f64 = 0.17;

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
    for (config_index, config) in configs.into_iter().enumerate() {
        for ratios in [
            Vec::new(),
            vec![-4.0, -2.0, -1.0, 0.5, 1.0],
            (0..500)
                .map(|i| f64::from(i % 50) - 25.0 + f64::from(i % 7) * 0.1)
                .collect(),
        ] {
            let histogram = our_histogram_of(&ratios);
            for target_fdr in TARGET_FALSE_DISCOVERY_RATES {
                let ours = ng::calibrate_from_the_ratio_histogram(&histogram, target_fdr, &config);
                let case = format!("{} ratios, target {target_fdr}", ratios.len());
                let key = format!(
                    "fallback_and_cut/{config_index}/{}/{target_fdr:?}",
                    ratios.len()
                );
                assert_the_same_prior(&case, &key, &ours.prior);
                let theirs = their_optional_f64(&format!("{key}/lr_threshold"));
                assert_eq!(
                    ours.lr_threshold.map(f64::to_bits),
                    theirs.map(f64::to_bits),
                    "{case}: ng and production cut at different ratios — ng {:?}, \
                     production {theirs:?}",
                    ours.lr_threshold,
                );
                for lr in [-100.0, -8.0, -1.0, 0.0, 1.0, 8.0, 30.0, 100.0, f64::NAN] {
                    assert_eq!(
                        ours.flags(lr),
                        their_bool(&format!("{key}/flags/{lr:?}")),
                        "{case}: ng and production disagree about removing a record at ratio {lr}"
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
    assert_every_answer_was_read("fallback_and_cut");
}
