//! The stutter model — how likely each length change is, for one repeat.
//!
//! **Not an alignment algorithm**, and it lives here anyway because *two* consumers share
//! it: the two-penalty best-path aligner in this module (the alignment spec §4.2), and the
//! genotyping likelihood outside it. One description, so the two cannot drift apart. It is
//! HipSTR's model.
//!
//! **The distribution is owned by `doc/devel/ng/spec/read_likelihoods.md` §4.2**, which
//! states it in full; the alignment spec's §5.2 is repointed there and no longer restates
//! it. Cite the first when this module's behaviour is in question.
//!
//! Everything here is a **linear probability**, not a logarithm — unlike HipSTR's own
//! fields, which are logs (`log_equal_`), and matching production's `stutter_pmf`, which
//! also returns a linear value. The module carries no `_ln` names for that reason; a value
//! that ever does become a logarithm must say so in its name (the crate-wide convention).
//!
//! ## The two regimes are the model's defining structure
//!
//! A read's length change is either a whole number of repeats or it is not, and the two mean
//! different things:
//!
//! - **a whole-repeat change** — a whole number of repeats. This is slippage, the common
//!   event, and its size is measured in **repeats**.
//! - **a part-repeat change** — not a whole number of repeats. This is a sequencing indel or
//!   an interruption, not slippage; it is rarer, it gets its **own** parameters, and its size
//!   is measured in **base pairs**.
//!
//! *(HipSTR calls these two **in frame** and **out of frame**. This project does not use
//! those words — `doc/devel/ng/spec/read_likelihoods.md` §1.3 bans them, because *frame* is
//! borrowed from coding sequence and was read here as meaning inside the tract against in
//! the flanks, which is a different distinction. HipSTR's field names are kept in the doc
//! comments below for whoever ports from it.)*
//!
//! Each regime splits again by direction, because stutter is asymmetric: **losing repeats is
//! more common than gaining them**.
//!
//! **Which regime applies is decided by arithmetic alone** — is the change a multiple of
//! the period? — and **never by what was actually inserted**. So an insertion that happens
//! to be period-sized is treated as slippage whether or not its bases are the motif. The
//! composition is caught downstream instead, as a base-error mismatch against the re-tiled
//! candidate. The mis-routing is worst at period 1, where it catches roughly three of every
//! four single-base insertions — which is why the algorithm-3-versus-4 comparison
//! (Milestone D) must score an indel *of the repeat's own base* separately from an indel of
//! a different base, or the two effects cancel in the average (the alignment spec §4.2).

use std::num::NonZeroU8;

/// Largest slip this model scores at all; anything past it is **zero**, so an implausibly
/// large change is not explained away as stutter — such a read falls to the genotyping's
/// outlier handling instead (`doc/devel/ng/spec/read_likelihoods.md` §4.2).
///
/// Copied from production's `MAX_SLIP` rather than imported: ng is a from-scratch caller
/// that does not depend on production (owner, 2026-07-16). **The value must stay equal to
/// production's** while the two models are meant to agree, and
/// `the_copied_cutoff_still_equals_productions` is what enforces that rather than trusting
/// this sentence. Production's is a compile-time array bound described there as a
/// provisional choice awaiting recalibration on the simulator.
///
/// # The one constant, two scales — inherited deliberately
///
/// Production applies this single number to the **repeat** count on the whole-repeat branch
/// and to the re-indexed **base-pair** count on the part-repeat branch. Those are different
/// scales, and `doc/devel/ng/spec/read_likelihoods.md` §4.2 has **decided against inheriting
/// that**: ng is to carry two cutoffs named for what they count — `max_whole_repeat_slip` in
/// repeats and `max_part_repeat_slip` in base pairs — because 10 repeats at a hexamer is 60
/// base pairs and 10 base pairs is not the same claim. **This constant is what stands until
/// they land** (plan step E2).
///
/// **Until then the behaviour is production's**, and here is what that costs. Part-repeat
/// changes are cut off about **`period − 1` times sooner** in real terms:
/// the re-indexed size is `Δ − Δ/period`, so a cutoff of 10 re-indexed steps admits roughly
/// `10 · period/(period − 1)` base pairs, against `10 · period` in whole repeats. At period
/// 4 that is about 13 bp for part-repeat changes against 40 bp for whole-repeat ones; **at
/// period 2 the effect vanishes** (about 20 against 20), and it grows with the period.
/// Changing this is a **behaviour change to a model the genotyping shares**, and this plan's
/// rule is transcribe first, change separately with its own evidence. Recorded as a
/// follow-up, not fixed here.
pub const MAX_SLIP: u32 = 10;

/// Lower bound on a one-step share, and the floor under
/// [`StutterModel::same_length_share`]. Production uses one constant for both, so this does
/// too. Named for HipSTR's `geom` fields, which is what a one-step share is called there.
///
/// **Public so the genotyping likelihood can name it rather than spell a second copy.** Its
/// contract says every probability is floored before a logarithm and asks for the floors as
/// named constants with their reasons (`doc/devel/ng/spec/read_likelihoods.md` §8); the
/// clamp on a one-step share is this one, and two spellings of one number are two things
/// that can drift apart. The clamping itself stays here, where the distribution is — a
/// consumer reads the value to document and to test against, never to apply.
pub const GEOM_MIN: f64 = 0.01;
/// Upper bound on a one-step share. See [`GEOM_MIN`], including why it is public.
pub const GEOM_MAX: f64 = 0.99;

/// The six stutter rates a [`StutterModel`] is built from, **named** rather than positional.
///
/// The names are load-bearing, not decoration. As six positional `f64` arguments, a
/// longer/shorter or whole-repeat/part-repeat transposition is invisible at every call site
/// *and* in a test suite whose fixtures happen to give the two members of a pair the same
/// value — which is exactly what the first version of this module's tests did. HipSTR keeps
/// the two one-step shares genuinely independent
/// (`doc/devel/ng/spec/read_likelihoods.md` §4.2), so the pairing matters.
///
/// The same-length share is **not** here: it is derived from these four masses, and taking
/// it as an input would let it disagree with the values that define it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StutterRates {
    /// Share of reads a whole repeat longer than the allele, at any size. HipSTR's `in_up_`.
    pub whole_repeat_longer_share: f64,
    /// Share of reads a whole repeat shorter than the allele, at any size. Usually exceeds
    /// [`Self::whole_repeat_longer_share`] — stutter is contraction-biased. HipSTR's
    /// `in_down_`.
    pub whole_repeat_shorter_share: f64,
    /// **Of the reads that slipped by whole repeats**, the share that moved by exactly one
    /// repeat. HipSTR's `in_geom_`.
    ///
    /// A share of the slips, *not* a rate of decline — see [`StutterModel::new`] for the
    /// trap.
    pub whole_repeat_one_step_share: f64,
    /// Share of reads longer by part of a repeat, at any size. HipSTR's `out_up_`.
    pub part_repeat_longer_share: f64,
    /// Share of reads shorter by part of a repeat, at any size. HipSTR's `out_down_`.
    pub part_repeat_shorter_share: f64,
    /// **Of the reads that changed by part of a repeat**, the share that moved by exactly one
    /// base. HipSTR's `out_geom_`.
    ///
    /// HipSTR keeps this independent of [`Self::whole_repeat_one_step_share`]; production
    /// ties the two to a single value, which is an undeclared placeholder rather than a
    /// fitted result (`doc/devel/ng/spec/read_likelihoods.md` §4.2).
    pub part_repeat_one_step_share: f64,
}

/// How likely each length change is, for one repeat.
///
/// # How it is built: per locus, from the *reference* allele
///
/// This model belongs to a **locus**, and is built from that locus's stutter shape and its
/// **reference** allele length — never from a candidate allele's length. The distinction is
/// not pedantic: stutter rises with allele length, so a per-*candidate* slip level is a real
/// and useful thing, but it belongs to the **genotyping**, which scores reads against
/// candidates. This module *measures*; letting the measurement's own model vary with the
/// candidate being tested would bias the ruler toward the answer under test (arch §5, §2.4).
///
/// **This type fits nothing**, which `doc/devel/ng/spec/read_likelihoods.md` §7 states as
/// the model's side of the parameter boundary: it reads its seven numbers frozen. Production
/// derives them per call from a per-locus stutter shape plus a per-read stutter level; ng's
/// come from the parameter pre-pass, per read group per stratum (§4.2 of the same document).
/// **The constructor that adapts a fit to these rates is `stutter_rates_for`, which plan
/// step E3 lands**; until then a caller builds [`StutterRates`] itself.
///
/// # Contract: the same-length share is floored, so **do not test that the five sum to one**
///
/// The same-length share is *defined* as whatever the four direction shares leave —
/// `1` minus the whole-repeat longer and shorter shares and the part-repeat longer and
/// shorter shares — **but it is floored**. When the floor binds, the five values sum to
/// slightly more than one. That is deliberate: the floor is what stops a hostile parameter
/// combination producing a negative probability. It means "the five sum to one" must **not**
/// be written as a test, and `the_five_shares_do_not_sum_to_one_when_the_floor_binds` exists
/// to make that explicit rather than merely absent.
///
/// (HipSTR instead *asserts* the shares sum below one at construction and never clamps.
/// Both disciplines are defensible; they differ, and only the floor matches the code being
/// ported — arch §2.4.)
///
/// # Why the fields are private
///
/// The clamps **are** the contract. Arch §2.4 sketches these as public fields, but public
/// fields would let a caller assemble a model whose one-step shares sit at 0 or 1, or whose
/// same-length share is negative, and every guarantee above would be a comment rather than a
/// fact. Construction goes through [`StutterModel::new`], which applies them; the seven
/// values are readable through accessors. (Arch states its signatures are illustrative and
/// the contract is the deliverable.)
///
/// Not `Copy`: at seven `f64`s this is 56 bytes, well past the size where an implicit copy
/// is free, and the per-call repeat context holds it **by reference** by design — so `Copy`
/// would buy nothing and hide the cost of the cases it did serve.
#[derive(Debug, Clone, PartialEq)]
pub struct StutterModel {
    same_length_share: f64,
    whole_repeat_longer_share: f64,
    whole_repeat_shorter_share: f64,
    whole_repeat_one_step_share: f64,
    part_repeat_longer_share: f64,
    part_repeat_shorter_share: f64,
    part_repeat_one_step_share: f64,
}

impl StutterModel {
    /// Build a model from the six rates. The same-length share is derived, as one minus the
    /// four direction shares, then floored.
    ///
    /// # The trap this type exists to make hard
    ///
    /// **A one-step share is the share of the slips that moved by exactly one step, not a
    /// per-step multiplier.** Of the reads that moved at all, this fraction moved one step;
    /// the mass at the next step is multiplied by `(1 − share)` for each step after, so a
    /// **larger** share concentrates mass on *single*-repeat slips — HipSTR ships 0.95,
    /// meaning nineteen slips in twenty are exactly one repeat. If a parameter arrives
    /// expressed as the probability of *continuing* to the next step (mean size
    /// `1/(1 − decay)`), it is the **complement**: `share = 1 − decay`. Getting this
    /// backwards inverts the size distribution — large slips become common — and **nothing
    /// crashes** (`doc/devel/ng/spec/read_likelihoods.md` §4.2).
    /// `a_single_repeat_slip_outweighs_a_larger_one` fails if it is inverted.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if any rate is not a finite value in `[0, 1]`. In release the
    /// assertion is compiled out and the values are sanitized instead (below), so the call
    /// never panics there.
    ///
    /// # What sanitizing does, and why it is not the same thing as the clamps
    ///
    /// Two different jobs, easy to conflate. The **clamps on the one-step shares and the
    /// floor under the same-length share are ported model behaviour** — production does
    /// exactly this — and they bind only at degenerate values. The **mass sanitizing is a
    /// release-mode backstop** for a violated precondition, and it exists because without it
    /// this type would break its own stated contract: an unchecked share of `2.0` makes
    /// `probability` return `1.9`, a negative share returns a negative probability, and
    /// `NaN` poisons every score while *reporting a healthy floor* — because `f64::max`
    /// absorbs `NaN` and `f64::clamp` passes it straight through, so the constructor would
    /// look fine and only the results would be wrong. A non-finite rate becomes `0` (no
    /// stutter), the no-information end of the scale.
    #[must_use]
    pub fn new(rates: StutterRates) -> Self {
        debug_assert!(
            every_rate(&rates)
                .iter()
                .all(|rate| rate.is_finite() && (0.0..=1.0).contains(rate)),
            "stutter rates must be finite probabilities in [0, 1]: {rates:?}"
        );
        Self::sanitized(rates)
    }

    /// Everything [`Self::new`] does **except** the debug assertion — that is, exactly what a
    /// release build runs.
    ///
    /// Separated so the release path is reachable from a test: the assertion fires in the
    /// test profile, so `new` cannot be used to check what happens once it is gone.
    fn sanitized(rates: StutterRates) -> Self {
        // Destructured exhaustively, with no `..`: a seventh rate added later must break
        // this line rather than be dropped here in silence. Measured before the change —
        // adding a field to `StutterRates` failed to compile at the two named-row
        // constructors and **nowhere else**, so the new rate would have been neither
        // validated nor sanitized nor stored.
        let StutterRates {
            whole_repeat_longer_share,
            whole_repeat_shorter_share,
            whole_repeat_one_step_share,
            part_repeat_longer_share,
            part_repeat_shorter_share,
            part_repeat_one_step_share,
        } = rates;

        let whole_repeat_longer_share = sanitize_direction_share(whole_repeat_longer_share);
        let whole_repeat_shorter_share = sanitize_direction_share(whole_repeat_shorter_share);
        let part_repeat_longer_share = sanitize_direction_share(part_repeat_longer_share);
        let part_repeat_shorter_share = sanitize_direction_share(part_repeat_shorter_share);

        Self {
            same_length_share: (1.0
                - whole_repeat_longer_share
                - whole_repeat_shorter_share
                - part_repeat_longer_share
                - part_repeat_shorter_share)
                .clamp(GEOM_MIN, 1.0),
            whole_repeat_longer_share,
            whole_repeat_shorter_share,
            whole_repeat_one_step_share: sanitize_one_step_share(whole_repeat_one_step_share),
            part_repeat_longer_share,
            part_repeat_shorter_share,
            part_repeat_one_step_share: sanitize_one_step_share(part_repeat_one_step_share),
        }
    }

    /// HipSTR's **shipped** default parameters, as a matched set.
    ///
    /// HipSTR has **two** parameter sets, and mixing them yields a pairing that exists
    /// nowhere — `doc/devel/ng/spec/alignment.md` §5.2 records that an earlier draft of the
    /// spec did exactly that. These constructors exist so the two rows cannot be crossed by
    /// hand: this one is whole-repeat one-step share 0.95 with 0.05/0.05, part-repeat
    /// one-step share 0.95 with 0.01/0.01.
    ///
    /// Note the shipped row makes expansion and contraction **equal**. That symmetry is a
    /// starting point, not a claim — HipSTR's *fitted* values are contraction-biased.
    #[must_use]
    pub fn hipstr_shipped() -> Self {
        Self::new(StutterRates {
            whole_repeat_longer_share: 0.05,
            whole_repeat_shorter_share: 0.05,
            whole_repeat_one_step_share: 0.95,
            part_repeat_longer_share: 0.01,
            part_repeat_shorter_share: 0.01,
            part_repeat_one_step_share: 0.95,
        })
    }

    /// HipSTR's EM **starting point** — an initialisation its fitting immediately moves away
    /// from, and **not** a default. Provided as a matched set for the same reason as
    /// [`Self::hipstr_shipped`]: whole-repeat one-step share 0.9 with 0.1/0.1, part-repeat
    /// one-step share 0.8 with 0.01/0.01.
    ///
    /// Note this row keeps the two one-step shares **different** (0.9 against 0.8), which is
    /// the shape HipSTR actually uses — production ties them to one value, an undeclared
    /// placeholder (`doc/devel/ng/spec/read_likelihoods.md` §4.2).
    #[must_use]
    pub fn hipstr_em_start() -> Self {
        Self::new(StutterRates {
            whole_repeat_longer_share: 0.1,
            whole_repeat_shorter_share: 0.1,
            whole_repeat_one_step_share: 0.9,
            part_repeat_longer_share: 0.01,
            part_repeat_shorter_share: 0.01,
            part_repeat_one_step_share: 0.8,
        })
    }

    /// `P(length change)` for a change of `bp_diff` bases on a repeat of period `period`.
    ///
    /// Both regimes are a direction share times a geometric over size, following
    /// `doc/devel/ng/spec/read_likelihoods.md` §4.2 term by term. Zero beyond [`MAX_SLIP`],
    /// so an implausibly large slip is not explained away.
    ///
    /// # Why part-repeat sizes are re-indexed
    ///
    /// The part-repeat geometric is indexed by `e = Δ − Δ/period` (truncated division), not
    /// by `Δ`. The reason is **not** double-counting — the two regimes are disjoint by
    /// construction, since a change is a part-repeat one precisely when it is *not* a
    /// multiple of the period, so no length can reach both. What the re-indexing does is
    /// **compress the ranks**: it maps the part-repeat values onto consecutive integers so
    /// the geometric's support has no gaps. At period 3 the part-repeat values 1, 2, 4, 5, 7
    /// map to e = 1, 2, 3, 4, 5. Without it the geometric would be evaluated at indices that
    /// skip the multiples, distorting the distribution.
    ///
    /// # Why `period` is a [`NonZeroU8`]
    ///
    /// The whole-repeat test divides by the period, so a zero would divide by zero — and a
    /// period-zero repeat is not a repeat in the first place. A `debug_assert!` would be the
    /// module's usual answer (arch §3), but here the illegal state is cheap to make
    /// **unrepresentable**, which is better than making it testable: there is no release
    /// path to get wrong, and no guard that compiles out. [`Motif`] already guarantees a
    /// period in `1..=6`, so the conversion at a real call site cannot fail.
    ///
    /// [`Motif`]: crate::ng::types::Motif
    /// [`NonZeroU8`]: std::num::NonZeroU8
    #[must_use]
    pub fn probability(&self, bp_diff: i64, period: NonZeroU8) -> f64 {
        let period = i64::from(period.get());

        if bp_diff % period == 0 {
            // A whole-repeat change: slippage, sized in whole repeats.
            let repeats = bp_diff / period;
            if repeats == 0 {
                self.same_length_share
            } else {
                Regime {
                    longer_share: self.whole_repeat_longer_share,
                    shorter_share: self.whole_repeat_shorter_share,
                    one_step_share: self.whole_repeat_one_step_share,
                }
                .probability(repeats)
            }
        } else {
            // A part-repeat change: a sequencing indel or an interruption, sized in base
            // pairs and rank-compressed (see above). The size is never zero here — `bp_diff`
            // is not a multiple of the period, so truncated division cannot consume all of
            // it.
            Regime {
                longer_share: self.part_repeat_longer_share,
                shorter_share: self.part_repeat_shorter_share,
                one_step_share: self.part_repeat_one_step_share,
            }
            .probability(bp_diff - bp_diff / period)
        }
    }

    /// Share of reads showing the allele's own length — neither longer nor shorter.
    /// **Floored** — see the type's contract. HipSTR's `log_equal_`, exponentiated.
    #[must_use]
    pub fn same_length_share(&self) -> f64 {
        self.same_length_share
    }
    /// Share of reads a whole repeat longer than the allele, at any size. HipSTR's `in_up_`.
    #[must_use]
    pub fn whole_repeat_longer_share(&self) -> f64 {
        self.whole_repeat_longer_share
    }
    /// Share of reads a whole repeat shorter than the allele, at any size. Usually exceeds
    /// [`Self::whole_repeat_longer_share`] — stutter is contraction-biased. HipSTR's
    /// `in_down_`.
    #[must_use]
    pub fn whole_repeat_shorter_share(&self) -> f64 {
        self.whole_repeat_shorter_share
    }
    /// **Of the reads that slipped by whole repeats**, the share that moved by exactly one
    /// repeat. HipSTR's `in_geom_`.
    #[must_use]
    pub fn whole_repeat_one_step_share(&self) -> f64 {
        self.whole_repeat_one_step_share
    }
    /// Share of reads longer by part of a repeat, at any size. HipSTR's `out_up_`.
    #[must_use]
    pub fn part_repeat_longer_share(&self) -> f64 {
        self.part_repeat_longer_share
    }
    /// Share of reads shorter by part of a repeat, at any size. HipSTR's `out_down_`.
    #[must_use]
    pub fn part_repeat_shorter_share(&self) -> f64 {
        self.part_repeat_shorter_share
    }
    /// **Of the reads that changed by part of a repeat**, the share that moved by exactly one
    /// base. HipSTR's `out_geom_`.
    #[must_use]
    pub fn part_repeat_one_step_share(&self) -> f64 {
        self.part_repeat_one_step_share
    }
}

/// The three shares one regime is scored from, **named at the call site** rather than passed
/// positionally.
///
/// Both regimes are the same shape over different parameters, so the arithmetic is written
/// once — but as three bare `f64` arguments in a row, a mispairing is invisible to the
/// compiler *and* to a reader, which is the hazard [`StutterRates`]'s own doc argues about
/// thirty lines above. Built as a literal at each of the two call sites, a crossed pair reads
/// as two disagreeing words on one line (`longer_share: self.whole_repeat_shorter_share`).
#[derive(Debug, Clone, Copy)]
struct Regime {
    /// The share of reads longer than the allele in this regime, at any size.
    longer_share: f64,
    /// The share shorter, at any size.
    shorter_share: f64,
    /// Of the reads that moved at all in this regime, the share that moved by one step.
    one_step_share: f64,
}

impl Regime {
    /// This regime's whole answer: pick the direction by the sign of `steps`, drop anything
    /// past the cutoff, and evaluate the geometric.
    ///
    /// Shared by both regimes, which also makes the "size is at least one" precondition
    /// provable in one place. `steps` is a **repeat** count on the whole-repeat branch and a
    /// re-indexed **base-pair** count on the part-repeat branch; [`MAX_SLIP`] is applied to
    /// whichever it is, which is the two-scale inheritance recorded on that constant.
    #[inline]
    fn probability(self, steps: i64) -> f64 {
        debug_assert!(
            steps != 0,
            "a zero-size change is the same-length case, not a regime"
        );
        let size = steps.unsigned_abs();
        if size > u64::from(MAX_SLIP) {
            return 0.0;
        }
        let share = if steps < 0 {
            self.shorter_share
        } else {
            self.longer_share
        };
        // `size >= 1` here, so the exponent cannot underflow. `unsigned_abs` also means
        // `i64::MIN` is safe, which a `-steps - 1` form would not be. The `MAX_SLIP` early
        // return above is what makes the cast safe: whatever replaces that cutoff must keep
        // the guard.
        share * self.one_step_share * (1.0 - self.one_step_share).powi((size - 1) as i32)
    }
}

/// Every rate a [`StutterRates`] carries, as an array — destructured exhaustively so that a
/// rate added later cannot slip past [`StutterModel::new`]'s validation unnoticed.
#[inline]
fn every_rate(rates: &StutterRates) -> [f64; 6] {
    let StutterRates {
        whole_repeat_longer_share,
        whole_repeat_shorter_share,
        whole_repeat_one_step_share,
        part_repeat_longer_share,
        part_repeat_shorter_share,
        part_repeat_one_step_share,
    } = *rates;
    [
        whole_repeat_longer_share,
        whole_repeat_shorter_share,
        whole_repeat_one_step_share,
        part_repeat_longer_share,
        part_repeat_shorter_share,
        part_repeat_one_step_share,
    ]
}

/// A direction share, made safe for release: `NaN` becomes zero (no stutter — the
/// no-information end), and anything outside `[0, 1]` is clamped into it.
#[inline]
fn sanitize_direction_share(share: f64) -> f64 {
    if share.is_nan() {
        0.0
    } else {
        share.clamp(0.0, 1.0)
    }
}

/// A one-step share, held strictly inside `(0, 1)`. `NaN` becomes [`GEOM_MIN`];
/// `f64::clamp` would otherwise pass it straight through.
#[inline]
fn sanitize_one_step_share(share: f64) -> f64 {
    if share.is_nan() {
        GEOM_MIN
    } else {
        share.clamp(GEOM_MIN, GEOM_MAX)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A repeat period, which the type system requires to be non-zero.
    fn period(bases: u8) -> NonZeroU8 {
        NonZeroU8::new(bases).expect("a test period is never zero")
    }

    /// A fixture whose six rates are **all different**, so a transposition of any pair
    /// changes an answer.
    ///
    /// This exists because the first version of these tests could not see one. Every fixture
    /// then used 0.05/0.05 for the whole-repeat shares, 0.01/0.01 for the part-repeat ones,
    /// and 0.95 for *both* one-step shares — so swapping longer with shorter, or the
    /// whole-repeat one-step share with the part-repeat one, passed all twelve tests
    /// byte-identically. Contraction-biased in both regimes, as HipSTR's fitted values are.
    fn all_distinct() -> StutterModel {
        StutterModel::new(StutterRates {
            whole_repeat_longer_share: 0.03,
            whole_repeat_shorter_share: 0.07,
            whole_repeat_one_step_share: 0.95,
            part_repeat_longer_share: 0.004,
            part_repeat_shorter_share: 0.012,
            part_repeat_one_step_share: 0.8,
        })
    }

    /// **Every accessor returns its own rate, on a fixture where all six differ.**
    ///
    /// `all_distinct` exists so that a transposition cannot hide, and until this test it
    /// could: seven tests used the fixture and every one of them went through
    /// [`StutterModel::probability`], while every test that read an **accessor** used a
    /// fixture whose longer and shorter shares are equal — `0.05/0.05`, `0.1/0.1`,
    /// `0.01/0.01`. Measured by this step's review: making `part_repeat_longer_share()`
    /// return the shorter field left the whole library green at 4,354 passing tests, while
    /// the accessor returned 0.012 where 0.004 is right. The two part-repeat accessors have
    /// no caller outside this module yet, and the genotyping likelihood is the named coming
    /// one — so the crossing would have been waiting for it rather than caught before it.
    #[test]
    fn every_accessor_returns_its_own_rate() {
        let model = all_distinct();
        assert_eq!(model.whole_repeat_longer_share(), 0.03);
        assert_eq!(model.whole_repeat_shorter_share(), 0.07);
        assert_eq!(model.whole_repeat_one_step_share(), 0.95);
        assert_eq!(model.part_repeat_longer_share(), 0.004);
        assert_eq!(model.part_repeat_shorter_share(), 0.012);
        assert_eq!(model.part_repeat_one_step_share(), 0.8);
        assert!((model.same_length_share() - (1.0 - 0.03 - 0.07 - 0.004 - 0.012)).abs() < 1e-15);
    }

    /// The published whole-repeat formula, term by term:
    /// `share · one_step · (1 − one_step)^(n − 1)`, in **both** directions with different
    /// shares.
    #[test]
    fn the_whole_repeat_branch_reproduces_the_published_formula() {
        let model = all_distinct();
        let period = period(3);
        for repeats in 1..=5i64 {
            let decay = 0.05f64.powi((repeats - 1) as i32);
            let bp = repeats * i64::from(period.get());
            assert!((model.probability(bp, period) - 0.03 * 0.95 * decay).abs() < 1e-15);
            assert!((model.probability(-bp, period) - 0.07 * 0.95 * decay).abs() < 1e-15);
        }
    }

    /// The published part-repeat formula, on the **re-indexed** size `e = Δ − Δ/period`,
    /// in **both** directions.
    ///
    /// The direction half is the point: with the two part-repeat shares equal — as every
    /// fixture in the first version of this file had — swapping them passed the entire
    /// suite, because no assertion anywhere used a negative part-repeat change.
    #[test]
    fn the_part_repeat_branch_reproduces_the_published_formula_in_both_directions() {
        let model = all_distinct();
        let period = period(3);
        for bp_diff in [1i64, 2, 4, 5, 7] {
            let effective = bp_diff - bp_diff / i64::from(period.get());
            let decay = 0.2f64.powi((effective - 1) as i32);
            assert!((model.probability(bp_diff, period) - 0.004 * 0.8 * decay).abs() < 1e-15);
            assert!((model.probability(-bp_diff, period) - 0.012 * 0.8 * decay).abs() < 1e-15);
        }
    }

    /// **Direction asymmetry must be expressible in both regimes.** Reads lose repeats more
    /// often than they gain them — at tomato dinucleotides by a factor of 4.9, 2,438 reads
    /// against 501 (`doc/devel/ng/spec/read_likelihoods.md` §4.2) — and a direction-symmetric
    /// model would be a step backwards from production's scoring, which already carries the
    /// asymmetry (the alignment spec §4.2, which requires a two-penalty aligner to inherit
    /// it).
    #[test]
    fn a_contraction_outscores_an_expansion_of_the_same_size_in_both_regimes() {
        let model = all_distinct();
        let period = period(4);
        let repeat = i64::from(period.get());

        // A whole-repeat change.
        assert!(model.probability(-repeat, period) > model.probability(repeat, period));
        // A part-repeat change — the case nothing used to cover.
        assert!(model.probability(-1, period) > model.probability(1, period));
    }

    /// **Rank compression, the reason for the re-indexing.** At period 3 the part-repeat
    /// changes 1, 2, 4, 5, 7 must map onto consecutive geometric steps 1, 2, 3, 4, 5 — so
    /// the support has no gaps. Indexing by Δ itself would skip the multiples and distort
    /// the distribution. Negatives mirror exactly, because Rust truncates toward zero.
    #[test]
    fn part_repeat_sizes_compress_onto_consecutive_ranks() {
        let period = 3i64;
        let ranks: Vec<i64> = [1i64, 2, 4, 5, 7]
            .iter()
            .map(|&bp_diff| bp_diff - bp_diff / period)
            .collect();
        assert_eq!(ranks, vec![1, 2, 3, 4, 5]);

        let negative_ranks: Vec<i64> = [-1i64, -2, -4, -5, -7]
            .iter()
            .map(|&bp_diff| bp_diff - bp_diff / period)
            .collect();
        assert_eq!(negative_ranks, vec![-1, -2, -3, -4, -5]);
    }

    /// **The inverted one-step-share trap** (`doc/devel/ng/spec/read_likelihoods.md` §4.2). A
    /// one-step share of 0.95 means nineteen slips in twenty are exactly one repeat, so a
    /// one-repeat slip must outweigh a two-repeat slip by exactly `1/(1 − share) = 20`. Read
    /// as a *decay* instead — the complement, 0.05 — the distribution inverts and large slips
    /// become common, and nothing crashes.
    ///
    /// **The quantitative assertion is the one doing the work.** Monotonicity alone would
    /// *not* catch the inversion: at a share of 0.05 the sequence still decreases, just far
    /// more slowly. The ratio is what pins it.
    #[test]
    fn a_single_repeat_slip_outweighs_a_larger_one() {
        let model = StutterModel::hipstr_shipped();
        let period = period(4);
        let repeat = i64::from(period.get());
        let one_repeat = model.probability(repeat, period);
        let two_repeats = model.probability(2 * repeat, period);
        let three_repeats = model.probability(3 * repeat, period);

        assert!(one_repeat > two_repeats && two_repeats > three_repeats);
        assert!((one_repeat / two_repeats - 20.0).abs() < 1e-9);
        assert!((two_repeats / three_repeats - 20.0).abs() < 1e-9);
    }

    /// **The same-length floor, and why "the five sum to one" must not be a test** (arch
    /// §2.4). With four direction shares summing past one, the derived same-length share would
    /// go negative; the floor stops that, and the consequence is that the five values then
    /// sum to *more* than one. Asserted here so the absence of a sums-to-one test reads as a
    /// decision rather than an oversight.
    #[test]
    fn the_five_shares_do_not_sum_to_one_when_the_floor_binds() {
        let hostile = StutterModel::new(StutterRates {
            whole_repeat_longer_share: 0.5,
            whole_repeat_shorter_share: 0.5,
            whole_repeat_one_step_share: 0.95,
            part_repeat_longer_share: 0.1,
            part_repeat_shorter_share: 0.1,
            part_repeat_one_step_share: 0.95,
        });
        assert!(
            hostile.same_length_share() > 0.0,
            "the floor must keep the same-length share positive"
        );
        assert_eq!(hostile.same_length_share(), 0.01);

        let total = hostile.same_length_share()
            + hostile.whole_repeat_longer_share()
            + hostile.whole_repeat_shorter_share()
            + hostile.part_repeat_longer_share()
            + hostile.part_repeat_shorter_share();
        assert!(total > 1.0, "expected the floor to push the total past one");
    }

    /// With well-behaved shares the floor does not bind, and the same-length share is exactly
    /// the remainder.
    #[test]
    fn the_same_length_share_is_the_remainder_when_the_floor_does_not_bind() {
        let model = StutterModel::hipstr_shipped();
        assert!((model.same_length_share() - (1.0 - 0.05 - 0.05 - 0.01 - 0.01)).abs() < 1e-15);
        assert_eq!(model.probability(0, period(4)), model.same_length_share());
    }

    /// The one-step shares are held strictly inside `(0, 1)`, so neither a certainty nor an
    /// impossibility is expressible — the *clamps carry weight* trap, which the alignment
    /// spec's §5.2 keeps because the read-likelihood spec does not state it. Asserted
    /// against the **contractual values** 0.01 and 0.99, not against the constants the
    /// implementation uses — otherwise editing a constant would move both sides of the
    /// assertion.
    #[test]
    fn the_one_step_shares_are_held_strictly_inside_zero_and_one() {
        let extreme = StutterModel::new(StutterRates {
            whole_repeat_longer_share: 0.05,
            whole_repeat_shorter_share: 0.05,
            whole_repeat_one_step_share: 0.0,
            part_repeat_longer_share: 0.01,
            part_repeat_shorter_share: 0.01,
            part_repeat_one_step_share: 1.0,
        });
        assert_eq!(extreme.whole_repeat_one_step_share(), 0.01);
        assert_eq!(extreme.part_repeat_one_step_share(), 0.99);
    }

    /// **Every ill-formed rate must still yield a probability**, because the debug assertion
    /// that rejects them is compiled out of the release build this project runs. Each of
    /// these produced a non-probability before the sanitizing: a share of 2.0 gave 1.9, a
    /// negative share gave a negative score, and `NaN` poisoned every score *while the
    /// constructor reported a healthy floor* — `f64::max` absorbs `NaN` and `f64::clamp`
    /// passes it through, so only the results were wrong.
    ///
    /// Goes through `sanitized`, which is `new` minus the debug assertion — the release path
    /// is unreachable through `new` in this profile.
    #[test]
    fn ill_formed_rates_still_yield_probabilities() {
        let well_formed = StutterRates {
            whole_repeat_longer_share: 0.05,
            whole_repeat_shorter_share: 0.05,
            whole_repeat_one_step_share: 0.95,
            part_repeat_longer_share: 0.01,
            part_repeat_shorter_share: 0.01,
            part_repeat_one_step_share: 0.95,
        };
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.5, 2.0, 1e300] {
            for slot in 0..6 {
                let mut rates = well_formed;
                match slot {
                    0 => rates.whole_repeat_longer_share = bad,
                    1 => rates.whole_repeat_shorter_share = bad,
                    2 => rates.whole_repeat_one_step_share = bad,
                    3 => rates.part_repeat_longer_share = bad,
                    4 => rates.part_repeat_shorter_share = bad,
                    _ => rates.part_repeat_one_step_share = bad,
                }
                let model = StutterModel::sanitized(rates);
                assert!(
                    model.same_length_share().is_finite(),
                    "rate {bad} in slot {slot}"
                );
                for period in (1..=6u8).map(period) {
                    for bp_diff in [-9i64, -4, -1, 0, 1, 4, 9] {
                        let probability = model.probability(bp_diff, period);
                        assert!(
                            probability.is_finite() && (0.0..=1.0).contains(&probability),
                            "rate {bad} in slot {slot} gave {probability} at Δ={bp_diff}"
                        );
                    }
                }
            }
        }
    }

    /// **Past the cutoff the answer is zero, not a small number** — an implausibly large
    /// change must not be explained away as stutter. Pinned on **both** sides of the
    /// boundary, in both regimes and both directions, so a `>` / `>=` slip cannot survive.
    #[test]
    fn slips_past_the_cutoff_score_zero() {
        let model = all_distinct();
        let period = period(4);
        let repeat = i64::from(period.get());
        let last = i64::from(MAX_SLIP) * repeat;

        // Whole-repeat: exactly at the cutoff is still scored; one repeat past is zero.
        assert!(model.probability(last, period) > 0.0);
        assert!(model.probability(-last, period) > 0.0);
        assert_eq!(model.probability(last + repeat, period), 0.0);
        assert_eq!(model.probability(-last - repeat, period), 0.0);

        // Part-repeat: the boundary is on the *re-indexed* size. Δ = 13 → e = 10 (scored),
        // Δ = 14 → e = 11 (zero).
        assert_eq!(13 - 13 / repeat, 10);
        assert_eq!(14 - 14 / repeat, 11);
        assert!(model.probability(13, period) > 0.0);
        assert!(model.probability(-13, period) > 0.0);
        assert_eq!(model.probability(14, period), 0.0);
        assert_eq!(model.probability(-14, period), 0.0);
    }

    /// **The one constant applied at two scales**, reproduced from production and pinned so
    /// the asymmetry is visible rather than surprising: on the whole-repeat branch the cutoff
    /// counts *repeats*, on the part-repeat branch it counts *re-indexed base pairs*, so
    /// part-repeat changes are cut off about `period − 1` times sooner in real terms.
    #[test]
    fn the_cutoff_counts_repeats_on_one_branch_and_base_pairs_on_the_other() {
        let model = all_distinct();

        // Period 4: 40 bp accepted as whole repeats, but 14 bp already rejected as a
        // part-repeat change.
        assert!(model.probability(40, period(4)) > 0.0);
        assert_eq!(model.probability(14, period(4)), 0.0);

        // Period 2: the effect very nearly vanishes — both scales admit about 20 bp, which
        // is why the rationale says `period − 1` and not `period`.
        assert!(model.probability(20, period(2)) > 0.0);
        assert_eq!(19 - 19 / 2, 10);
        assert!(model.probability(19, period(2)) > 0.0);
        assert_eq!(model.probability(21, period(2)), 0.0);
    }

    /// Period 1 is the case the comparison most needs, not the one to skip: **every**
    /// change is a whole-repeat one there, so the part-repeat branch is unreachable and the
    /// arithmetic regime split collapses. What does *not* collapse is direction, size decay
    /// and placement multiplicity (the alignment spec §4.2, which works through what the
    /// collapse costs).
    #[test]
    fn every_change_is_a_whole_repeat_change_at_period_one() {
        let model = all_distinct();
        for bp_diff in [-3i64, -1, 1, 2, 5] {
            let share = if bp_diff < 0 { 0.07 } else { 0.03 };
            let steps = bp_diff.unsigned_abs();
            let expected = share * 0.95 * 0.05f64.powi((steps - 1) as i32);
            assert!((model.probability(bp_diff, period(1)) - expected).abs() < 1e-15);
        }
        // Direction asymmetry survives at period 1 — the one-penalty model cannot express
        // it at all, which is what Milestone D's comparison is about.
        assert!(model.probability(-1, period(1)) > model.probability(1, period(1)));
    }

    /// Every probability the model returns is a real probability — finite and in `[0, 1]` —
    /// across both regimes, both directions, several periods, and past the cutoff.
    #[test]
    fn every_probability_is_finite_and_within_zero_and_one() {
        for model in [
            StutterModel::hipstr_shipped(),
            StutterModel::hipstr_em_start(),
            all_distinct(),
        ] {
            for period in (1..=6u8).map(period) {
                for bp_diff in -60i64..=60 {
                    let probability = model.probability(bp_diff, period);
                    assert!(
                        probability.is_finite(),
                        "non-finite at Δ={bp_diff}, period {period:?}"
                    );
                    assert!(
                        (0.0..=1.0).contains(&probability),
                        "probability {probability} out of range at Δ={bp_diff}, period {period:?}"
                    );
                }
            }
        }
    }

    /// The two HipSTR parameter rows are matched sets, and mixing them yields a pairing that
    /// exists nowhere (`doc/devel/ng/spec/alignment.md` §5.2 records an earlier draft of the
    /// spec doing exactly that). The named constructors are what stop that happening by hand;
    /// this pins their contents.
    #[test]
    fn the_two_hipstr_parameter_sets_are_kept_as_matched_rows() {
        let shipped = StutterModel::hipstr_shipped();
        assert_eq!(shipped.whole_repeat_one_step_share(), 0.95);
        assert_eq!(shipped.part_repeat_one_step_share(), 0.95);
        assert_eq!(shipped.whole_repeat_longer_share(), 0.05);
        assert_eq!(shipped.whole_repeat_shorter_share(), 0.05);
        assert_eq!(shipped.part_repeat_longer_share(), 0.01);

        let em_start = StutterModel::hipstr_em_start();
        assert_eq!(em_start.whole_repeat_one_step_share(), 0.9);
        assert_eq!(em_start.part_repeat_one_step_share(), 0.8);
        assert_eq!(em_start.whole_repeat_longer_share(), 0.1);
        assert_eq!(em_start.part_repeat_longer_share(), 0.01);

        // The rows differ in their one-step shares — which is the whole reason to keep them
        // apart, and also why HipSTR treats the two as independent.
        assert_ne!(
            em_start.whole_repeat_one_step_share(),
            em_start.part_repeat_one_step_share()
        );
    }

    /// **The copied cutoff must not drift from production's.** The doc on [`MAX_SLIP`]
    /// asserts they stay equal while the two models are meant to agree; this makes that true
    /// rather than aspirational. A **test-only** reference, so shipping ng code still depends
    /// on nothing in production.
    #[test]
    fn the_copied_cutoff_still_equals_productions() {
        assert_eq!(
            MAX_SLIP as usize,
            crate::ssr::cohort::param_estimation::MAX_SLIP,
            "ng's MAX_SLIP has drifted from production's"
        );
    }
}
