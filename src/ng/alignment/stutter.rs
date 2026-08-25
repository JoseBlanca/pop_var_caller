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

use std::num::{NonZeroU8, NonZeroU32};

/// Largest **whole-repeat** slip this model scores at all, counted in **repeats**; anything
/// past it is **zero**, so an implausibly large change is not explained away as stutter —
/// such a read falls to the genotyping's outlier handling instead
/// (`doc/devel/ng/spec/read_likelihoods.md` §4.2).
///
/// # Inherited, not measured
///
/// **The value is production's provisional 10, declared inherited rather than fitted.**
/// Neither cutoff has a source in the parameters the pre-pass fits, and production's own
/// comment calls its 10 "a provisional choice"
/// ([`param_estimation.rs`](../../../../src/ssr/cohort/param_estimation.rs)). Copied rather
/// than imported: ng is a from-scratch caller that does not depend on production (owner,
/// 2026-07-16), and `the_copied_cutoffs_still_equal_productions` keeps the two equal rather
/// than trusting this sentence.
///
/// **What should eventually set it is the mass it discards**, which
/// [`StutterModel::unreachable_mass`] now computes per candidate.
pub const MAX_WHOLE_REPEAT_SLIP: u32 = 10;

/// Largest **part-repeat** change this model scores at all, counted in **re-indexed steps**
/// — the scale `Δ − Δ/period`, which is what makes the geometric's support gapless (see
/// [`StutterModel::probability`]).
///
/// **A step is not a base pair, and saying so is this constant's whole point.** Removing the
/// multiples of the period compresses the ranks, so ten steps admit about
/// `10 · period/(period − 1)` base pairs — roughly 13 at period 4, 20 at period 2. Spec §4.2
/// calls this "`max_part_repeat_slip = 10` base pairs"; that is the unit the split exists to
/// stop conflating, and the correction is owed to that document.
///
/// # Why this is a second constant and not the same one
///
/// Production applies **one** number to the repeat count on the whole-repeat branch and to
/// the re-indexed base-pair count on the part-repeat branch. Those are different scales, and
/// `doc/devel/ng/spec/read_likelihoods.md` §4.2 decided against inheriting that: **10
/// repeats at a hexamer is 60 base pairs, and 10 base pairs is not the same claim.**
///
/// Naming them separately does not by itself change any number — both are 10, so this step
/// is behaviour-preserving — but it makes the two independently settable by whoever measures
/// them, which one shared constant did not. What the split costs today, kept here because it
/// is the thing a measurement would move: a part-repeat change is cut off about
/// `period − 1` times sooner in real terms — about 13 bp against 40 bp at period 4; **at
/// period 2 the effect vanishes** (about 20 against 20), and it grows with the period.
///
/// Inherited on the same terms as [`MAX_WHOLE_REPEAT_SLIP`] — see there.
pub const MAX_PART_REPEAT_SLIP: u32 = 10;

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
    /// look fine and only the results would be wrong.
    ///
    /// **What each ill-formed value becomes**, because "sanitized" on its own says nothing
    /// about which end of the scale it lands at. `NaN` and `−∞` become the **least**-stutter
    /// end — `0` in a direction slot, [`GEOM_MIN`] in a one-step slot — which is the
    /// no-information answer. **`+∞` and any overshoot clamp to the top of the slot's range
    /// instead**, `1.0` or [`GEOM_MAX`]: they are not missing information, they are a
    /// magnitude past the range, and the nearest legal value is the honest reading of them.
    /// `sanitizing_maps_each_ill_formed_rate_to_its_documented_value` is the table, and it
    /// exists because a sanitizer that sent `NaN` to the *most*-stutter end satisfied every
    /// assertion this module had.
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
    /// `doc/devel/ng/spec/read_likelihoods.md` §4.2 term by term. Zero beyond each branch's
    /// own cutoff — [`MAX_WHOLE_REPEAT_SLIP`] in repeats, [`MAX_PART_REPEAT_SLIP`] in
    /// re-indexed base pairs — so an implausibly large slip is not explained away.
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
                    max_steps: MAX_WHOLE_REPEAT_SLIP,
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
                max_steps: MAX_PART_REPEAT_SLIP,
            }
            .probability(bp_diff - bp_diff / period)
        }
    }

    /// **The mass this model never scores for a candidate of this shape** — one minus the
    /// total it puts on everything a read of that candidate could actually show.
    ///
    /// # Why it has to be reported rather than assumed small
    ///
    /// A model that quietly loses mass on some candidates and not others is **comparing them
    /// on different scales**, and the genotyping is a comparison between candidates. Spec
    /// §4.2 requires this computed and surfaced per candidate for exactly that reason, and
    /// spec §12's fifth test pins that it is computed and surfaced — **not that it is
    /// small**, because it is not always small.
    ///
    /// # The three things that are missing, and their sizes
    ///
    /// 1. **The cutoffs.** Each branch drops everything past its own limit, which costs that
    ///    branch's direction shares times `(1 − one_step_share)^cutoff`. At HipSTR's shipped
    ///    one-step share of 0.95 over ten steps that is a factor of `0.05^10`, around 1 part
    ///    in 10¹³ — negligible. At a one-step share of 0.5 it is `0.5^10`, about 1 in a
    ///    thousand of that branch's mass, which is not.
    /// 2. **Contractions the tract is too short to reach.** A read of this candidate must
    ///    still show a repeat, so at most `repeat_count − 1` of them can go: the whole-repeat
    ///    contraction geometric is cut there whenever that is below the cutoff, and the
    ///    part-repeat one at the `(repeat_count − 1) · (period − 1)` non-multiples inside the
    ///    part of the tract that can go. **This is the term that varies per candidate.**
    ///    Measured here on the shortest tract the copy floors admit — four repeats, at
    ///    hexamers, at a slippage level of 2 in 100 split 4:1 toward contraction — it is
    ///    **2.0 parts in a million** at a one-step share of 0.95 and **2.0 parts in a
    ///    thousand** at 0.5. A thousandfold range over one fitted parameter, which is why the
    ///    step chance is fitted and not defaulted;
    ///    `the_shortest_admissible_tract_loses_the_size_the_specification_states` pins both.
    ///
    ///    **Which boundary this is, because the specification says two things.** Spec §4.2's
    ///    prose calls the unreachable slips "contracting away *more repeats than exist*",
    ///    which would make losing all `repeat_count` reachable; its worked figures call the
    ///    unreachable tail at a four-repeat tract "a contraction of **four repeats or
    ///    more**", which makes losing all of them unreachable. **The two readings differ by
    ///    one step, and only the second reproduces the two sizes the specification states** —
    ///    exactly, to the digits quoted. This follows the figures. Flipping it is one
    ///    subtraction, and the test above would need its two numbers changed with it.
    /// 3. **At period 1, the whole part-repeat branch.** Every length change is a whole
    ///    number of repeats when the motif is one base, so the part-repeat branch is
    ///    unreachable and *all* of its mass is lost — 2 in 100 at HipSTR's shipped values.
    ///    This is the largest of the three by far, and the one most easily mistaken for a
    ///    defect. It is a property of the model, and reporting it is what keeps a
    ///    mononucleotide candidate comparable with a dinucleotide one.
    ///
    /// # Why `repeat_count` is a [`NonZeroU32`]
    ///
    /// **A candidate with no repeats is not a candidate**, and the rule above says why: a
    /// read of this candidate must still show a repeat, so a tract holding none can show
    /// nothing at all and there is no distribution to account for. As a plain `u32` the
    /// arithmetic answered for zero exactly as it does for one — `saturating_sub(1)` makes
    /// them the same input — and the two are nothing like the same answer: at HipSTR's
    /// shipped values, period 3, one repeat loses **1.68 in 100** while a zero-repeat tract
    /// under this function's own rule leaves **99.6 in 100** unreachable. That is the largest
    /// mis-scaling this function can produce, and it would have arrived silently.
    ///
    /// The same reasoning as [`Self::probability`]'s period: the illegal state is cheap to
    /// make **unrepresentable**, which is better than making it testable — there is no
    /// release path to get wrong and no guard that compiles out.
    ///
    /// [`NonZeroU32`]: std::num::NonZeroU32
    ///
    /// # The floor case
    ///
    /// [`Self::same_length_share`] is floored, so a hostile parameter row can make the five
    /// shares sum to **more** than one and the arithmetic here go negative. There is no
    /// negative mass to report, so the answer is clamped at zero — a degenerate model has no
    /// truncation to account for, it has a construction problem, and
    /// `a_floored_model_reports_no_loss_rather_than_a_negative_one` pins that.
    ///
    /// # Closed form, checked against enumeration
    ///
    /// This sums five geometric tails rather than walking the support, because it is called
    /// once per candidate per read group. `the_reported_loss_equals_one_minus_the_reachable_sum`
    /// is the guard: it enumerates the support by calling [`Self::probability`] over every
    /// length a read could show and requires the two to agree. Deriving the closed form and
    /// checking it by enumeration is the point — a version that enumerated would make that
    /// test compare a thing with itself.
    #[must_use]
    pub fn unreachable_mass(&self, period: NonZeroU8, repeat_count: NonZeroU32) -> f64 {
        let period_bases = u32::from(period.get());
        let repeat_count = repeat_count.get();

        // A geometric's first `steps` terms, `1 − (1 − one_step)^steps` — which is already
        // `0.0` at no steps at all, since `powi(0)` is one.
        let reached = |one_step: f64, steps: u32| -> f64 {
            1.0 - (1.0 - one_step).powi(i32::try_from(steps).unwrap_or(i32::MAX))
        };

        // How far a read of this candidate can contract: the tract keeps at least one
        // repeat, so `repeat_count - 1` of them can go. See the doc above for why this
        // boundary and not `repeat_count`.
        let contractable_repeats = repeat_count.saturating_sub(1);

        // Whole repeats: expansion is unbounded above, contraction stops at the tract's own
        // length.
        let whole_longer = self.whole_repeat_longer_share
            * reached(self.whole_repeat_one_step_share, MAX_WHOLE_REPEAT_SLIP);
        let whole_shorter = self.whole_repeat_shorter_share
            * reached(
                self.whole_repeat_one_step_share,
                MAX_WHOLE_REPEAT_SLIP.min(contractable_repeats),
            );

        // Part repeats: at period 1 there are none at all, because every change is a whole
        // number of repeats. Otherwise the reachable contractions are the non-multiples of
        // the period inside the contractable part of the tract — `(period - 1)` of them per
        // repeat, the same count the whole-repeat rule leaves.
        let (part_longer_steps, part_shorter_steps) = if period_bases == 1 {
            (0, 0)
        } else {
            (
                MAX_PART_REPEAT_SLIP,
                MAX_PART_REPEAT_SLIP.min(contractable_repeats.saturating_mul(period_bases - 1)),
            )
        };
        let part_longer = self.part_repeat_longer_share
            * reached(self.part_repeat_one_step_share, part_longer_steps);
        let part_shorter = self.part_repeat_shorter_share
            * reached(self.part_repeat_one_step_share, part_shorter_steps);

        let reachable =
            self.same_length_share + whole_longer + whole_shorter + part_longer + part_shorter;
        (1.0 - reachable).max(0.0)
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
#[derive(Debug)]
struct Regime {
    /// The share of reads longer than the allele in this regime, at any size.
    longer_share: f64,
    /// The share shorter, at any size.
    shorter_share: f64,
    /// Of the reads that moved at all in this regime, the share that moved by one step.
    one_step_share: f64,
    /// The largest size this regime scores, **on this regime's own scale** — repeats on the
    /// whole-repeat branch, re-indexed base pairs on the part-repeat one. Carried per regime
    /// rather than read from one shared constant, which is the whole of what E2 changes.
    max_steps: u32,
}

impl Regime {
    /// This regime's whole answer: pick the direction by the sign of `steps`, drop anything
    /// past the cutoff, and evaluate the geometric.
    ///
    /// Shared by both regimes, which also makes the "size is at least one" precondition
    /// provable in one place. `steps` is a **repeat** count on the whole-repeat branch and a
    /// re-indexed **base-pair** count on the part-repeat branch, and [`Self::max_steps`] is
    /// the cutoff for whichever scale this regime counts in.
    #[inline]
    fn probability(self, steps: i64) -> f64 {
        debug_assert!(
            steps != 0,
            "a zero-size change is the same-length case, not a regime"
        );
        let size = steps.unsigned_abs();
        if size > u64::from(self.max_steps) {
            return 0.0;
        }
        let share = if steps < 0 {
            self.shorter_share
        } else {
            self.longer_share
        };
        // `size >= 1` here, so the exponent cannot underflow. `unsigned_abs` also means
        // `i64::MIN` is safe, which a `-steps - 1` form would not be. The `max_steps` early
        // return above is what makes the cast safe: any future cutoff must keep the guard.
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

    /// A candidate's repeat count, which the type system requires to be non-zero — a tract
    /// holding no repeats is not a candidate ([`StutterModel::unreachable_mass`]).
    fn repeats(count: u32) -> NonZeroU32 {
        NonZeroU32::new(count).expect("a test candidate always holds a repeat")
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

    /// **Rank compression, the reason for the re-indexing**, asserted *through the model*.
    /// At period 3 the part-repeat changes 1, 2, 4, 5, 7 must land on consecutive geometric
    /// steps 1, 2, 3, 4, 5, so the support has no gaps; indexing by Δ itself would skip the
    /// multiples and distort the distribution.
    ///
    /// **The ratio form is the point.** This test's predecessor recomputed
    /// `bp_diff - bp_diff / period` in its own body and compared the result with
    /// `[1, 2, 3, 4, 5]` — a statement about Rust's truncating division that built no
    /// `StutterModel` and called nothing in this module, so it held whatever `probability`
    /// did. Measured by this step's review across 50 mutations, it never once killed one,
    /// including the two that corrupt the very re-indexing it was named for. Asserting the
    /// *ratio between consecutive values* instead re-spells nothing the implementation says,
    /// so a shared misconception cannot hide in it: each step down must be exactly
    /// `1 − one_step_share`, in both directions.
    #[test]
    fn part_repeat_probabilities_step_down_one_rank_at_a_time() {
        let model = all_distinct();
        let period = period(3);
        let ratio = 1.0 - model.part_repeat_one_step_share();

        for pair in [1i64, 2, 4, 5, 7].windows(2) {
            let (small, large) = (pair[0], pair[1]);
            for sign in [1i64, -1] {
                let step = model.probability(sign * large, period)
                    / model.probability(sign * small, period);
                assert!(
                    (step - ratio).abs() < 1e-12,
                    "Δ={} to Δ={} moved by {step}, not one rank ({ratio})",
                    sign * small,
                    sign * large
                );
            }
        }
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
        let last = i64::from(MAX_WHOLE_REPEAT_SLIP) * repeat;

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

    /// **Two constants at two scales**, pinned so the asymmetry is visible rather than
    /// surprising: on the whole-repeat branch the cutoff counts *repeats*, on the part-repeat
    /// branch *re-indexed base pairs*, so part-repeat changes are cut off about `period − 1`
    /// times sooner in real terms.
    ///
    /// **What no test here can check while both constants hold 10.** Which constant feeds
    /// which branch is not observable when they are equal — measured: wiring the part-repeat
    /// regime to [`MAX_WHOLE_REPEAT_SLIP`] leaves every test in this module green. That is
    /// inherent rather than a gap to fill: the split is about making the two *settable*
    /// apart, and the moment a measurement moves either one, the boundaries asserted below
    /// become the check. Named here so a later reader does not mistake the silence for
    /// coverage.
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
    ///
    /// **All twelve values, which it did not always cover.** It used to name five of
    /// `hipstr_shipped`'s six and four of `hipstr_em_start`'s, and this step's review
    /// measured what that left open: editing `hipstr_em_start`'s whole-repeat shorter share
    /// from 0.1 to 0.2, or its part-repeat shorter share from 0.01 to 0.02, changed the
    /// distribution and left the whole suite green. A test whose stated job is that a row
    /// cannot be edited by hand has to name every number in it.
    #[test]
    fn the_two_hipstr_parameter_sets_are_kept_as_matched_rows() {
        let shipped = StutterModel::hipstr_shipped();
        assert_eq!(shipped.whole_repeat_longer_share(), 0.05);
        assert_eq!(shipped.whole_repeat_shorter_share(), 0.05);
        assert_eq!(shipped.whole_repeat_one_step_share(), 0.95);
        assert_eq!(shipped.part_repeat_longer_share(), 0.01);
        assert_eq!(shipped.part_repeat_shorter_share(), 0.01);
        assert_eq!(shipped.part_repeat_one_step_share(), 0.95);

        let em_start = StutterModel::hipstr_em_start();
        assert_eq!(em_start.whole_repeat_longer_share(), 0.1);
        assert_eq!(em_start.whole_repeat_shorter_share(), 0.1);
        assert_eq!(em_start.whole_repeat_one_step_share(), 0.9);
        assert_eq!(em_start.part_repeat_longer_share(), 0.01);
        assert_eq!(em_start.part_repeat_shorter_share(), 0.01);
        assert_eq!(em_start.part_repeat_one_step_share(), 0.8);

        // The rows differ in their one-step shares — which is the whole reason to keep them
        // apart, and also why HipSTR treats the two as independent.
        assert_ne!(
            em_start.whole_repeat_one_step_share(),
            em_start.part_repeat_one_step_share()
        );
    }

    /// **What sanitizing turns each ill-formed rate into**, not merely that what comes out is
    /// still a probability. `ill_formed_rates_still_yield_probabilities` reaches this path for
    /// all 36 slot-and-value combinations and asserts finiteness and range only — so a
    /// sanitizer that mapped `NaN` to `1.0`, the *most* stutter rather than the least,
    /// satisfied it just as well. Measured by this step's review: that mutation, and the
    /// one-step equivalent sending `NaN` to [`GEOM_MAX`], both survived the whole suite.
    ///
    /// Note what the table says about `+∞`: on a direction share it clamps to **1.0**, not to
    /// zero. `StutterModel::new`'s "the least-stutter end" holds for `NaN` and `−∞`; an
    /// overshoot is a magnitude past the range rather than missing information, and clamps to
    /// the nearest legal value.
    ///
    /// The last two rows are well formed as a direction share and outside the one-step
    /// bounds, which is the only thing that separates the real clamp from one that fires
    /// at the exact endpoints alone — a share of 0.001 reaches an aligner as an
    /// almost-free slip extension, `ln(1 − 0.001) ≈ −0.001`.
    #[test]
    fn sanitizing_maps_each_ill_formed_rate_to_its_documented_value() {
        let well_formed = StutterRates {
            whole_repeat_longer_share: 0.05,
            whole_repeat_shorter_share: 0.05,
            whole_repeat_one_step_share: 0.95,
            part_repeat_longer_share: 0.01,
            part_repeat_shorter_share: 0.01,
            part_repeat_one_step_share: 0.95,
        };
        // (what goes in, what a direction slot becomes, what a one-step slot becomes)
        let table = [
            (f64::NAN, 0.0, GEOM_MIN),
            (f64::NEG_INFINITY, 0.0, GEOM_MIN),
            (-0.5, 0.0, GEOM_MIN),
            (f64::INFINITY, 1.0, GEOM_MAX),
            (2.0, 1.0, GEOM_MAX),
            (1e300, 1.0, GEOM_MAX),
            (0.001, 0.001, GEOM_MIN),
            (0.999, 0.999, GEOM_MAX),
        ];
        for (bad, as_direction, as_one_step) in table {
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
                let (got, want) = match slot {
                    0 => (model.whole_repeat_longer_share(), as_direction),
                    1 => (model.whole_repeat_shorter_share(), as_direction),
                    2 => (model.whole_repeat_one_step_share(), as_one_step),
                    3 => (model.part_repeat_longer_share(), as_direction),
                    4 => (model.part_repeat_shorter_share(), as_direction),
                    _ => (model.part_repeat_one_step_share(), as_one_step),
                };
                assert_eq!(got, want, "rate {bad} in slot {slot}");
            }
        }
    }

    /// `Regime::probability`'s comment says `unsigned_abs` keeps `i64::MIN` safe where a
    /// `-steps - 1` form would overflow — panicking in debug and wrapping in release.
    /// **Nothing exercised that**: the widest change any test passed was ±60. A doc-comment
    /// invariant with no test is a claim, and this makes it a fact.
    #[test]
    fn the_extreme_length_changes_score_zero_rather_than_overflowing() {
        let model = all_distinct();
        for period in (1..=6u8).map(period) {
            for bp_diff in [i64::MIN, i64::MIN + 1, i64::MAX - 1, i64::MAX] {
                assert_eq!(
                    model.probability(bp_diff, period),
                    0.0,
                    "Δ={bp_diff} at period {period:?} is far past the cutoff"
                );
            }
        }
    }

    proptest::proptest! {
        /// **Every row of rates yields a distribution, not only the three this file names.**
        ///
        /// The sweep above walks three hand-built models over 2,178 points — three corners of
        /// a six-dimensional rate space. That was enough while every consumer used
        /// hard-coded parameters, and it stops being enough at plan step E3, where
        /// `stutter_rates_for` starts handing this type **fitted** numbers: a row the
        /// pre-pass produced has never been through this suite.
        ///
        /// Two invariants, both of which the genotyping likelihood will rely on. Every value
        /// is a real probability — a negative or a `NaN` here poisons a log downstream while
        /// reporting nothing. And **past the cutoff the answer is exactly zero**, on whichever
        /// scale the regime counts, so an implausibly large change falls to the outlier term
        /// instead of being explained away as stutter.
        #[test]
        fn any_rates_yield_a_distribution_that_is_zero_past_the_cutoff(
            longer in 0.0f64..=1.0,
            shorter in 0.0f64..=1.0,
            whole_one_step in 0.0f64..=1.0,
            part_longer in 0.0f64..=1.0,
            part_shorter in 0.0f64..=1.0,
            part_one_step in 0.0f64..=1.0,
            bp_diff in -400i64..=400,
            period_bases in 1u8..=6,
        ) {
            let model = StutterModel::new(StutterRates {
                whole_repeat_longer_share: longer,
                whole_repeat_shorter_share: shorter,
                whole_repeat_one_step_share: whole_one_step,
                part_repeat_longer_share: part_longer,
                part_repeat_shorter_share: part_shorter,
                part_repeat_one_step_share: part_one_step,
            });
            let period = period(period_bases);
            let probability = model.probability(bp_diff, period);

            proptest::prop_assert!(
                probability.is_finite() && (0.0..=1.0).contains(&probability),
                "Δ={} period={} gave {}",
                bp_diff,
                period_bases,
                probability
            );

            // The cutoff counts repeats on one branch and re-indexed base pairs on the other,
            // so the size is whichever the regime is sized in.
            let period_bases = i64::from(period_bases);
            let size = if bp_diff % period_bases == 0 {
                (bp_diff / period_bases).unsigned_abs()
            } else {
                (bp_diff - bp_diff / period_bases).unsigned_abs()
            };
            let cutoff = if bp_diff % period_bases == 0 {
                MAX_WHOLE_REPEAT_SLIP
            } else {
                MAX_PART_REPEAT_SLIP
            };
            if size > u64::from(cutoff) {
                proptest::prop_assert_eq!(
                    probability,
                    0.0,
                    "Δ={} is {} steps past the cutoff and was not zero",
                    bp_diff,
                    size
                );
            }
        }

        /// **The reported loss equals what an enumeration measures, for every row of rates**
        /// — not only the eleven the fixed sweep names.
        ///
        /// The sweep is thorough where it looks, but it looks at two hand-built rows. The
        /// consumer that arrives at plan step E3 hands this type **fitted** numbers, and a
        /// fitted row has never been through a fixed list. This is also the shape of test
        /// that would have found the zero-repeat case without anyone thinking to write it.
        ///
        /// **The direction ranges stay under a half each so the same-length floor does not
        /// bind**, which is what makes the equality exact rather than clamped; the floored
        /// case is `a_floored_model_reports_no_loss_rather_than_a_negative_one`'s.
        #[test]
        fn any_rates_report_the_loss_an_enumeration_measures(
            longer in 0.0f64..=0.4,
            shorter in 0.0f64..=0.4,
            whole_one_step in GEOM_MIN..=GEOM_MAX,
            part_longer in 0.0f64..=0.05,
            part_shorter in 0.0f64..=0.05,
            part_one_step in GEOM_MIN..=GEOM_MAX,
            repeat_count in 1u32..=40,
            period_bases in 1u8..=6,
        ) {
            let model = StutterModel::new(StutterRates {
                whole_repeat_longer_share: longer,
                whole_repeat_shorter_share: shorter,
                whole_repeat_one_step_share: whole_one_step,
                part_repeat_longer_share: part_longer,
                part_repeat_shorter_share: part_shorter,
                part_repeat_one_step_share: part_one_step,
            });
            let period = period(period_bases);
            let repeat_count = repeats(repeat_count);
            let reported = model.unreachable_mass(period, repeat_count);
            let enumerated = (1.0 - reachable_mass(&model, period, repeat_count)).max(0.0);
            proptest::prop_assert!(
                (reported - enumerated).abs() < 1e-12,
                "period {} at {} repeats: reported {}, enumerated {}",
                period_bases,
                repeat_count,
                reported,
                enumerated
            );
        }
    }

    /// **Spec §12's fourth test: the distribution sums to one over its whole support.**
    ///
    /// Sum `probability` over every length change the model scores, for periods 2 to 6 and
    /// direction splits from symmetric to five-to-one, and the total must be one. **No
    /// production test pins this**, and it is the test that catches three silent failures at
    /// once: a one-step share read as its complement, a mis-set same-length share, and a
    /// part-repeat re-indexing off by one.
    ///
    /// **How each of the three shows up as a shortfall.** A complemented one-step share is
    /// the loud one: at a fitted 0.95 read as 0.05, the ten scored steps hold only
    /// `1 − 0.95^10`, four tenths of a branch's mass, so six tenths of every slip vanishes. A
    /// same-length share that is not the remainder moves the total by exactly its own error.
    /// And a part-repeat geometric indexed by `Δ` rather than by the compressed rank skips
    /// the multiples of the period, so its weights sum to less than one.
    ///
    /// **Period 1 is excluded and gets its own assertion below**, because there the
    /// part-repeat branch is unreachable by construction and the total is *supposed* to fall
    /// short by exactly that branch's mass — which
    /// `a_mononucleotide_candidate_loses_the_whole_part_repeat_mass` already pins.
    ///
    /// The tolerance is the cutoffs' own tail, `(1 − one_step)^10`, which is what an
    /// untruncated sum would recover. At the fitted-like shares swept here it is below
    /// `1e-9`; the second assertion states the exact form, cross-checking this sum against
    /// [`StutterModel::unreachable_mass`], which is derived a completely different way.
    #[test]
    fn the_distribution_sums_to_one_over_its_whole_support() {
        // Symmetric, then two-to-one, then five-to-one toward contraction.
        for shorter_of_the_slips in [0.5, 2.0 / 3.0, 5.0 / 6.0] {
            for level in [0.02, 0.2] {
                for one_step in [0.9, 0.95, 0.99] {
                    let model = StutterModel::new(StutterRates {
                        whole_repeat_longer_share: level * (1.0 - shorter_of_the_slips),
                        whole_repeat_shorter_share: level * shorter_of_the_slips,
                        whole_repeat_one_step_share: one_step,
                        part_repeat_longer_share: level * 0.05 * (1.0 - shorter_of_the_slips),
                        part_repeat_shorter_share: level * 0.05 * shorter_of_the_slips,
                        part_repeat_one_step_share: one_step,
                    });
                    for period_bases in 2..=6u8 {
                        let period = period(period_bases);
                        // A tract long enough that no contraction is out of reach, so the
                        // only thing missing is the cutoffs' tail.
                        let repeat_count = repeats(200);
                        let total = reachable_mass(&model, period, repeat_count);

                        assert!(
                            (total - 1.0).abs() < 1e-9,
                            "period {period_bases}, level {level}, one-step {one_step}, \
                             split {shorter_of_the_slips}: summed to {total}"
                        );

                        // Exactly, against the model's own account of what it discards.
                        let lost = model.unreachable_mass(period, repeat_count);
                        assert!(
                            (total - (1.0 - lost)).abs() < 1e-12,
                            "period {period_bases}: summed {total}, reported loss {lost}"
                        );
                    }
                }
            }
        }
    }

    /// **A complemented one-step share is what the sums-to-one test is for.** Reading a
    /// fitted fall-off as a one-step share — the trap spec §4.2 names first — leaves the ten
    /// scored steps holding four tenths of each branch's mass instead of all of it, so the
    /// distribution sums visibly short. Measured here rather than asserted: this is the
    /// shortfall the tripwire above would report.
    #[test]
    fn a_complemented_one_step_share_makes_the_distribution_sum_short() {
        let level = 0.2;
        let complemented = StutterModel::new(StutterRates {
            whole_repeat_longer_share: level * 0.17,
            whole_repeat_shorter_share: level * 0.83,
            whole_repeat_one_step_share: 0.05, // the fitted 0.95, read backwards
            part_repeat_longer_share: level * 0.05 * 0.17,
            part_repeat_shorter_share: level * 0.05 * 0.83,
            part_repeat_one_step_share: 0.05,
        });
        let total = reachable_mass(&complemented, period(3), repeats(200));
        let short_by = 1.0 - total;

        // `1 − 0.95^10` of each branch survives, so 0.95^10 of the whole slip mass is gone:
        // 0.2 × 1.05 × 0.5987 ≈ 0.1257.
        assert!(
            (short_by - 0.125_74).abs() < 1e-4,
            "a complemented share left the total short by {short_by}"
        );
        assert!(
            short_by > 1e-9,
            "the tripwire's tolerance must be far below this"
        );
    }

    /// **The copied cutoffs must not drift from production's.** Both are inherited from its
    /// single provisional 10, and their docs say so; this makes that true rather than
    /// aspirational. A **test-only** reference, so shipping ng code still depends on nothing
    /// in production.
    ///
    /// **ng carries two where production carries one**, named for the scale each counts in.
    /// Splitting the name is what makes them independently settable by whoever measures them;
    /// until someone does, both must equal the number they were copied from.
    #[test]
    fn the_copied_cutoffs_still_equal_productions() {
        assert_eq!(
            MAX_WHOLE_REPEAT_SLIP as usize,
            crate::ssr::cohort::param_estimation::MAX_SLIP,
            "ng's whole-repeat cutoff has drifted from production's"
        );
        assert_eq!(
            MAX_PART_REPEAT_SLIP as usize,
            crate::ssr::cohort::param_estimation::MAX_SLIP,
            "ng's part-repeat cutoff has drifted from production's"
        );
    }

    /// Every length a read of a candidate with `repeat_count` repeats of `period` bases could
    /// actually show: the tract keeps at least one repeat, so the deepest contraction is
    /// `repeat_count − 1` repeats, and nothing beyond either cutoff is scored. Written as an
    /// enumeration, deliberately — it is the oracle for
    /// [`StutterModel::unreachable_mass`], which uses a closed form.
    fn reachable_mass(model: &StutterModel, period: NonZeroU8, repeat_count: NonZeroU32) -> f64 {
        let period_bases = i64::from(period.get());
        let deepest_contraction = i64::from(repeat_count.get() - 1) * period_bases;
        // Far past both cutoffs on either scale, so the sum is complete.
        let widest = period_bases * i64::from(MAX_WHOLE_REPEAT_SLIP + MAX_PART_REPEAT_SLIP + 4)
            + deepest_contraction;
        (-deepest_contraction..=widest)
            .map(|bp_diff| model.probability(bp_diff, period))
            .sum()
    }

    /// **Spec §12's fifth test: truncation removes a stated mass, and the model reports it.**
    ///
    /// Sum the distribution over everything a read of this candidate could show, subtract
    /// from one, and require the *reported* loss to equal that — for every candidate length,
    /// every period, and one-step shares across the clamped range.
    ///
    /// **What this pins is that the loss is computed and surfaced, not that it is small.**
    /// An earlier version of the specification compared it against "a named bound" and named
    /// none, which is unrunnable. The size is a property of the parameters, and this sweep
    /// spans it: the assertions below record the extremes it actually reaches.
    ///
    /// The oracle enumerates by calling `probability`; the implementation sums five
    /// geometric tails. That is the whole point of the test — two routes to one number, so a
    /// closed form that drops a term has something to disagree with.
    #[test]
    fn the_reported_loss_equals_one_minus_the_reachable_sum() {
        let mut smallest = f64::INFINITY;
        let mut largest = 0.0f64;

        for one_step in [GEOM_MIN, 0.1, 0.5, 0.9, 0.95, GEOM_MAX] {
            for part_one_step in [GEOM_MIN, 0.5, GEOM_MAX] {
                let model = StutterModel::new(StutterRates {
                    whole_repeat_longer_share: 0.004,
                    whole_repeat_shorter_share: 0.016,
                    whole_repeat_one_step_share: one_step,
                    part_repeat_longer_share: 0.0002,
                    part_repeat_shorter_share: 0.0008,
                    part_repeat_one_step_share: part_one_step,
                });
                for period_bases in 1..=6u8 {
                    let period = period(period_bases);
                    for repeat_count in [1u32, 2, 3, 4, 5, 6, 9, 10, 11, 30].map(repeats) {
                        let reported = model.unreachable_mass(period, repeat_count);
                        let enumerated = 1.0 - reachable_mass(&model, period, repeat_count);
                        assert!(
                            (reported - enumerated).abs() < 1e-12,
                            "period {period_bases}, {repeat_count} repeats, one-step \
                             {one_step}/{part_one_step}: reported {reported}, \
                             enumerated {enumerated}"
                        );
                        if period_bases > 1 {
                            smallest = smallest.min(reported);
                        }
                        largest = largest.max(reported);
                    }
                }
            }
        }

        // The sweep's own extremes, asserted so the range is a measurement rather than a
        // sentence — **2 in 100 down to nothing at all**, across parameters that are all
        // inside §8's clamps. The widest is a single-repeat tract at the slowest fall-off:
        // it can contract by nothing, so its whole contraction mass is unreachable, and at
        // period 1 the part-repeat branch goes too. The narrowest is exactly zero, where a
        // one-step share of 0.99 puts every reachable step's mass within a rounding of one.
        assert!(
            (largest - 2.061_752_830_003_5e-2).abs() < 1e-14,
            "the widest loss this sweep reaches is {largest}"
        );
        assert_eq!(
            smallest, 0.0,
            "the narrowest loss this sweep reaches is {smallest}"
        );
    }

    /// **The two sizes spec §4.2 states, reproduced.** The shortest tract the copy floors
    /// admit is four repeats, at hexamers; at a slippage level of 2 in 100 split 4:1 toward
    /// contraction, the mass the model cannot place is **2.0 parts in a million** at a
    /// one-step share of 0.95 and **2.0 parts in a thousand** at 0.5.
    ///
    /// **These two numbers are what settle which boundary the specification meant.** Its
    /// prose says the unreachable slips are those "contracting away more repeats than
    /// exist", which would let a four-repeat tract lose all four; its figures call the
    /// unreachable tail "a contraction of four repeats or more", which does not. Only the
    /// second reproduces the sizes above — the first gives 1.0 in ten million and 1.0 in a
    /// thousand, each one step of the geometric away. So this test is the record of a reading
    /// as much as a guard on the arithmetic.
    #[test]
    fn the_shortest_admissible_tract_loses_the_size_the_specification_states() {
        let at_one_step_share = |one_step: f64| {
            StutterModel::new(StutterRates {
                whole_repeat_longer_share: 0.004,
                whole_repeat_shorter_share: 0.016,
                whole_repeat_one_step_share: one_step,
                part_repeat_longer_share: 0.0002,
                part_repeat_shorter_share: 0.0008,
                part_repeat_one_step_share: one_step,
            })
            .unreachable_mass(period(6), repeats(4))
        };

        let at_shipped_share = at_one_step_share(0.95);
        assert!(
            (at_shipped_share - 2.0e-6).abs() < 1e-14,
            "expected 2 parts in a million, got {at_shipped_share}"
        );

        let at_slow_share = at_one_step_share(0.5);
        assert!(
            (at_slow_share - 2.0e-3).abs() < 1e-5,
            "expected 2 parts in a thousand, got {at_slow_share}"
        );

        // A thousandfold range over one fitted parameter — the reason the step chance is
        // fitted rather than defaulted.
        assert!(at_slow_share / at_shipped_share > 900.0);
    }

    /// **At period 1 the whole part-repeat branch is unreachable**, because every length
    /// change is a whole number of repeats when the motif is one base — so all of that
    /// branch's mass is lost, and the report must say so. At HipSTR's shipped values that is
    /// 2 in 100, far larger than either cutoff's tail, and it is the loss most easily
    /// mistaken for a defect.
    #[test]
    fn a_mononucleotide_candidate_loses_the_whole_part_repeat_mass() {
        let model = StutterModel::hipstr_shipped();
        // Long enough that no contraction is out of reach, so the part-repeat branch is the
        // only thing missing.
        let lost = model.unreachable_mass(period(1), repeats(40));
        let part_repeat_mass = model.part_repeat_longer_share() + model.part_repeat_shorter_share();
        assert!(
            (lost - part_repeat_mass).abs() < 1e-12,
            "period 1 lost {lost}, and the part-repeat mass is {part_repeat_mass}"
        );
        assert!((part_repeat_mass - 0.02).abs() < 1e-12);

        // The same model at period 2 loses essentially nothing: the branch is reachable
        // there, and both cutoffs' tails are 0.05^10.
        assert!(model.unreachable_mass(period(2), repeats(40)) < 1e-11);
    }

    /// **A short tract loses the contractions it cannot reach**, and that is the term that
    /// varies per candidate. A one-repeat tract cannot contract at all — a read of it must
    /// still show a repeat — while a thirty-repeat tract reaches every step of both
    /// geometrics.
    #[test]
    fn a_short_tract_loses_more_than_a_long_one() {
        let model = StutterModel::hipstr_shipped();
        let short = model.unreachable_mass(period(3), repeats(1));
        let long = model.unreachable_mass(period(3), repeats(30));
        assert!(
            short > long,
            "a one-repeat tract lost {short}, a thirty-repeat one {long}"
        );

        // A one-repeat tract can contract by nothing at all, so **both** contraction shares
        // are unreachable in full — 0.05 whole-repeat plus 0.01 part-repeat, six parts in a
        // hundred, which is nothing like negligible. That is the case this report exists for.
        let both_contraction_shares =
            model.whole_repeat_shorter_share() + model.part_repeat_shorter_share();
        assert!(
            (short - both_contraction_shares).abs() < 1e-12,
            "a one-repeat tract lost {short}, against {both_contraction_shares} expected"
        );
        assert!((short - 0.06).abs() < 1e-12);

        // Thirty repeats reach every step of both geometrics, so all that is missing is the
        // two cutoffs' tails at 0.05^10 — eleven orders of magnitude smaller.
        assert!(long < 1e-13, "a thirty-repeat tract lost {long}");
    }

    /// **A degenerate model reports no loss rather than a negative one.** The same-length
    /// share is floored, so four direction shares summing past one make the five sum to
    /// *more* than one — and one minus that is negative. There is no negative mass to
    /// report: such a row has a construction problem, not a truncation to account for.
    #[test]
    fn a_floored_model_reports_no_loss_rather_than_a_negative_one() {
        let hostile = StutterModel::new(StutterRates {
            whole_repeat_longer_share: 0.5,
            whole_repeat_shorter_share: 0.5,
            whole_repeat_one_step_share: 0.95,
            part_repeat_longer_share: 0.1,
            part_repeat_shorter_share: 0.1,
            part_repeat_one_step_share: 0.95,
        });
        // **Against the enumeration, not against a range.** `lost >= 0.0 && lost.is_finite()`
        // is satisfied by a constant, by an absolute value, and by any of the five terms
        // being wrong; comparing with the oracle pins the clamp exactly and everything else
        // with it.
        let mut clamped_cells = 0;
        for period_bases in 1..=6u8 {
            for repeat_count in [1u32, 5, 30] {
                let period = period(period_bases);
                let repeat_count = repeats(repeat_count);
                let lost = hostile.unreachable_mass(period, repeat_count);
                let raw = 1.0 - reachable_mass(&hostile, period, repeat_count);
                assert!(
                    (lost - raw.max(0.0)).abs() < 1e-12,
                    "period {period_bases}, {repeat_count} repeats reported {lost} \
                     against {raw}"
                );
                if raw < 0.0 {
                    clamped_cells += 1;
                    assert_eq!(lost, 0.0);
                }
            }
        }

        // **The clamp must actually fire somewhere in this sweep**, or the test has stopped
        // covering what it is named for — and it must not fire everywhere, or the sweep has
        // stopped covering the honest answer. **Measured: 12 of the 18 cells clamp.** The six
        // that do not are the ones where genuinely unreachable mass exceeds this row's
        // over-allocation: at period 1 with a single repeat the model can place only 0.51 —
        // no contraction is reachable and the part-repeat branch does not exist — so 0.49 is
        // a real loss rather than a negative one to clamp away.
        assert_eq!(
            clamped_cells, 12,
            "the clamp fired at {clamped_cells} of 18 cells"
        );
    }
}
