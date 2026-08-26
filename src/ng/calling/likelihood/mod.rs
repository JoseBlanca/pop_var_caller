//! Step 7 — how probable this sample's reads are, given a genotype.
//!
//! One call answers *given that this sample's copies of this locus are such-and-such, how
//! probable are the reads we actually saw here?* — one number per candidate genotype, in
//! log space. The caller adds those numbers to the genotype prior
//! ([`super::genotype_prior`], step 8) and normalises; the sum is the posterior ng emits
//! (`doc/devel/ng/spec/read_likelihoods.md` §1).
//!
//! ## Two paths, one formula
//!
//! Both paths compute the same thing — a copy-weighted mixture over the genotype's
//! alleles, wrapped in one logarithm per observation — and differ in a single term, the
//! **emission**: the probability that one copy of one allele produced one observed
//! sequence. At an ordinary site a read can show a base the individual does not carry, and
//! one rate covers it. At a repeat tract a read can also show a whole repeat more or fewer
//! than the individual carries, because the copying enzyme slipped, and that happens
//! thousands of times more often than a base is misread (spec §2.1, §2.2).
//!
//! ## What this module never does
//!
//! **It fits nothing.** Every number it reads is fitted before calling starts, by
//! [`crate::ng::parameter_estimation`]: the per-read-group error rate and the scale that
//! calibrates it, the contamination fraction, and the four numbers that say how a repeat
//! tract slips — which are fitted per **stratum**, one *(motif length, repeat count)* cell,
//! because how much a tract slips depends on how many repeats it has more than on anything
//! else. This module does not choose the candidate alleles, does not run the caller's loop,
//! does not compute a posterior or a genotype, and does not decide whether a site is
//! emitted (spec §1.2).
//!
//! ## Vocabulary a reader needs before the types below
//!
//! An **observation** is one distinct sequence the reads showed at a locus, with the
//! number of reads that showed it — not one read, and not one allele. Its identity is
//! `(bases, read witness, read group)`, so the same sequence seen by a read that spanned
//! the whole locus and by one that ran out inside it is two observations. A **complete**
//! observation pins what the sample carries there; a **partial** one proves only that the
//! locus is *at least* what the read got through, which the statistician calls
//! **censored** (spec §1.3).
//!
//! A **read group** is one `@RG` — one lane, one chemistry. It is part of an observation's
//! identity here for a reason that is easy to lose: the formula puts a logarithm outside a
//! sum, so an observation's reads may be pooled into one term **only if every one of them
//! would have got the same number**, and two reads showing the same bases from two lanes
//! have different error rates (spec §2.3).
//!
//! An allele's **projection** is the bases it spells over one locus's span — what the
//! merge unified the samples' reads onto, so that two samples showing the same sequence
//! land on one allele. It is a whole-span thing, which is why a read that saw only part of
//! the span cannot be compared against it (spec §5.3).
//!
//! ## Two allele tables, and they are not the same numbering
//!
//! **This is the trap this module's types are shaped around.** The cohort merge unifies
//! every distinct sequence the whole cohort showed into one table, and its rows index that
//! table. Candidate *selection* — step 6, which does not exist yet — then keeps some of
//! those alleles and drops the rest, and **dropping allele *k* renumbers every allele above
//! it**. [`AlleleId`] names an index into the surviving table, the one a genotype is built
//! over, and nothing in the two numberings makes them agree.
//!
//! So there is no conversion from a merge row to an [`AlleleId`] that this module could
//! perform on its own, and it does not offer one: the mapping arrives as an argument
//! ([`GenericObservation::fill_from_supported_alleles`]), and the reads on alleles the
//! mapping drops are exactly the reads whose quality the formula pools into
//! [`GenericSampleEvidence::unmatched_q_sum`].
//!
//! ## What a row function promises
//!
//! One call computes **one sample's whole row** — one log-probability per candidate
//! genotype, parallel to the loop's genotype table. Both paths promise the same six things,
//! and each is a property some later step's test has to keep:
//!
//! - **A pure function** of the evidence, the candidates and the frozen parameters. No
//!   clock, no random numbers, nothing read from global state.
//! - **Bit-identical at any thread count**, which means the sum over observations runs in
//!   one fixed order — the order the merge sorted them into. Floating-point addition is not
//!   associative, so an order that varies is a run whose genotypes vary (spec §8).
//! - **The row itself allocates nothing**: every buffer it works in is the caller's, held
//!   once per worker and reused across every sample of every locus ([`SsrRowScratch`], and
//!   the generic path's own from Milestone D). Production lifted exactly these out of its
//!   iteration after a profile put the allocator's self-time at about a sixth of cycles.
//!   **The buffers themselves are amortised rather than allocation-free** — they grow when a
//!   sample or locus is wider than any before it and never shrink — and
//!   [`GenericEvidenceBuffer`], which is where the evidence lives rather than where the row
//!   scribbles, is the same.
//! - **An empty evidence row is all zeros, with no branch to make it so.** An empty sum is
//!   zero, so a sample that showed nothing scores every genotype alike and the prior decides
//!   — which is the right answer rather than a special case (spec §3.3).
//! - **A mis-shaped input is a caller bug and asserts, structurally in release too**: a row
//!   whose length is not the genotype count, a candidate with no stratum entry, a parameter
//!   that is not finite. Production holds the analogous assertion in release because a
//!   scratch array too short for the allele count would otherwise be indexed out of bounds
//!   silently (`per_group_merger.rs:1963`).
//! - **Every probability is floored before a logarithm**, so one impossible observation
//!   yields a finite very negative number instead of turning a sample's whole row into
//!   `NaN`. The floor is [`MIN_BASE_ERROR`] here, and the geometric clamps in the stutter
//!   distribution. **Floored and not capped**, which C1 corrected on 2026-08-24 and the owner
//!   settled the same day: production pairs that floor with a ceiling at a half, and a ceiling
//!   binds on a single read and not on the fold of that read with others — a non-linear
//!   function of a per-read quality, which is exactly what the aggregation contract forbids.
//!   The ceiling and the method that applied it are deleted;
//!   `tests::what_the_row_charges_a_poor_read_is_not_capped` is what stops them coming back.
//!
//! ## Three tiers of parameter, and only the middle one is open
//!
//! Spec §6.1 sorts every number this module reads into three tiers. The tier decides where
//! the number travels, which is why it is documented on the types rather than in prose
//! somewhere:
//!
//! | tier | which numbers | how they arrive |
//! |---|---|---|
//! | **frozen for the whole run** | the per-read-group error rate and its calibration scale, the STR substitution rate, the contamination **fraction** | fields of the run-level views — [`ReadGroupCalibration`], [`ContaminationView`] — which nothing downstream may write |
//! | **re-estimable per locus, off by default** | the slippage level, the direction split, the fall-off | per call, inside the STR scoring context, so the caller's loop can re-fit them with no change here. This is the one constraint spec §6.1 makes binding |
//! | **re-estimated every iteration** | the per-locus allele frequencies | the prior reads them — and so does one term here: the contaminating population's frequency for the allele an observation shows, which moves with the loop (spec §3.6, corrected 2026-08-24) |
//!
//! [`ContaminationMixture`] is where the first and third tiers meet: it holds the frozen
//! fraction beside the frequency that moves, so a caller cannot hand the row one without the
//! other. **What follows for the loop is that a row cannot be cached across iterations
//! wherever contamination is on** — the emission reads no frequency and is still computed
//! once per `(sample, observation, candidate)`, but the row's own arithmetic moves.
//!
//! **The two halves of the contamination mixture sit in different tiers, and that is the
//! whole of what changed on 2026-08-24.** How contaminated a library is is a property of
//! that library and of nothing else, so the fraction is frozen before the loop starts. What
//! a contaminating read *shows* is a property of the locus, and a per-locus quantity is what
//! a per-locus loop is for.

pub mod generic;
pub mod ssr;
pub mod ssr_emission;
pub mod stutter_rates;

use crate::ng::calling::genotype_prior::{COUNT_PATH_DESYNC_THRESHOLD, SampleAlleleCopies};
use crate::ng::locus_generation::{ReadWitness, SequenceObservation, SsrDetail};
use crate::ng::parameter_estimation::generic::calibration::MintedReadErrors;
use crate::ng::parameter_estimation::joint::contamination::{
    ContaminationEstimate, ContaminationSource,
};
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::run::cohort_merge::build::{PartialObservation, SupportedAllele};
use crate::ng::types::{
    AlleleId, BatchId, BatchOfEachReadGroup, BatchOfEachSample, ErrorRate, Ploidy, ReadGroupId,
};

/// The widest ploidy either row builds a copy-share table for.
///
/// **`Ploidy::try_new` rejects only zero**, so seventeen copies is constructible and would
/// index past the array. Both rows assert on it, in release as well as debug, rather than
/// panicking with `index out of bounds` — every other caller bug in these files says in a
/// sentence what went wrong, and this one used to be the exception.
///
/// **Shared rather than written once per path**, because the two rows build the same table
/// from the same ploidy and a limit that differed between them would be a limit nobody could
/// state.
pub(crate) const MAX_PLOIDY_COPIES: usize = 16;

/// `k / P` for every copy count a genotype can carry — the chance a read was copied from a
/// copy carrying a given allele, for `k` copies of it out of `P` copies of the genome.
///
/// **Built once and shared by both rows**, because they build the same table from the same
/// ploidy for the same reason, and a limit or a formula that differed between them would be one
/// nobody could state. The entry for `k = 0` is never read — a read from no copy is the error
/// side on one path and contributes nothing to the sum on the other — and is filled with a
/// value that would be visible if it ever were.
///
/// # Panics
///
/// On a ploidy past [`MAX_PLOIDY_COPIES`]. `Ploidy::try_new` rejects only zero, so seventeen
/// copies is constructible and would index past the array; the rows say in a sentence what went
/// wrong rather than panicking with `index out of bounds`.
pub(crate) fn copy_shares(ploidy: Ploidy) -> [f64; MAX_PLOIDY_COPIES + 1] {
    let copies_of_the_genome = usize::from(ploidy.get());
    assert!(
        copies_of_the_genome <= MAX_PLOIDY_COPIES,
        "a sample with {copies_of_the_genome} copies of its genome is past the \
         {MAX_PLOIDY_COPIES} this row builds a copy-share table for"
    );
    let copies_of_the_genome = f64::from(ploidy.get());
    let mut shares = [f64::NAN; MAX_PLOIDY_COPIES + 1];
    for (copies, share) in shares.iter_mut().enumerate().skip(1) {
        *share = copies as f64 / copies_of_the_genome;
    }
    shares
}

/// The clamps a geometric success probability is held inside, **re-exported and not copied**.
///
/// **Why these are re-exported where the floor below is copied, which looks
/// inconsistent and is not.** These live in `alignment/stutter.rs`, which is ng's own code, so
/// there is one implementation and one name and a consumer points at it — two spellings of one
/// number are two things that can drift, and the tree already shows what that costs: production
/// keeps a *third* private copy of these same two values in its own stutter model with nothing
/// connecting it to either. The floor below comes from `src/var_calling/`, which is frozen
/// production that ng does not depend on, so it is copied deliberately and pinned by a test
/// — the same discipline `alignment/stutter.rs` uses for the slip cutoff it inherited.
pub use crate::ng::alignment::stutter::{GEOM_MAX, GEOM_MIN};

/// Floor under any per-read error probability before a logarithm is taken of it.
///
/// **Inherited from production at `1e-12` and declared inherited**
/// (`var_calling::contamination_estimation::MIN_BASE_ERROR`), not measured here.
/// `tests::the_floors_still_equal_productions` pins the equality, so the two cannot drift
/// while they are meant to describe one model assumption.
///
/// What it buys: an observation the genotype cannot explain at all yields about −27.6 nats
/// rather than negative infinity, so one such read makes a genotype very unlikely instead of
/// making the sample's whole row `NaN` (spec §8).
pub const MIN_BASE_ERROR: f64 = 1e-12;

/// §3.2's calibration: one multiplier per read group, so that the average error the model
/// charges that group's reads is the rate the parameter fit measured.
///
/// ```text
///                    fitted error rate for this read group
/// scale  =  ──────────────────────────────────────────────────────────
///           geometric mean of the minted error over that group's reads
/// ```
///
/// **The caller has no per-read error probability to rescale**, and it is worth being exact
/// about that because the looser picture is what made the arithmetic mean look reachable. The
/// merge keeps, per `(allele, read group)` per locus, how many reads support it and the sum of
/// their *log* error probabilities; the reads are gone from there on. Exactly one average is
/// recoverable from those two numbers, and it is the geometric one — which is why the scale is
/// **one addition in log space per observation** and never a multiplication read by read.
/// Nothing about the arithmetic changes: scaling every read and scaling their geometric mean
/// are the same operation, `exp(Σ ln(s·ε) / n) = s · exp(Σ ln ε / n)` (owner, 2026-08-24).
///
/// **Why both halves and not either one alone.** The reads' own qualities carry the only
/// information that tells one read from another, and at three reads a position that is the
/// whole call — one alternative read at Phred 40 and one at Phred 13 are not the same
/// evidence, and a single fitted rate says they are. The fitted rate carries the only
/// information about what the instrument cannot see — a mismapped read, a chimera, a
/// chimeric fragment, DNA from another individual — and it was measured on this data rather
/// than asserted by the machine. The scale keeps the shape from the first and the size from
/// the second (spec §3.2).
///
/// **The specification originally asked for the arithmetic mean** (corrected 2026-08-24,
/// owner), and the sharper statement of why it could not be that is the paragraph above:
/// there is nowhere in the model a per-read `ε` survives to be averaged that way. The walk
/// sums the *logarithms* and throws the individual reads away, and `Σ ε` cannot be recovered
/// from `Σ ln ε`. Taking the geometric mean is
/// also the self-consistent choice rather than a concession, because it is what the model
/// charges an observation — `exp(q_sum / num_reads)` — so a scale built from an arithmetic
/// mean and applied to a geometric one would not make the calibrated property hold in the
/// model's own terms.
///
/// **The two are far apart on real reads, and the size is now measured rather than assumed:
/// a factor of 25.2 on the 63-accession tomato cohort and 44.1 on HG002 at 300×**
/// (`examples/ng_minted_error_means.rs`, 2026-08-24; spec §3.2 carries the table). Per read
/// group on tomato the ratio runs 22.7 to 37.0, median 24.4, with no read group anywhere near
/// one. So building the scale from the arithmetic mean would have divided every charged error
/// by 25 to 44 — **14 to 16 Phred**, every read treated as that much cleaner than the fit
/// measured it to be.
///
/// **And a second reason that does not depend on the self-consistency argument at all: the
/// arithmetic mean is mostly measuring how often mates overlap.** A read's minted error is
/// the worse of its base and mapping qualities, and the mate-overlap rule silences the losing
/// mate of an overlapping pair by giving it base quality Phred 0 — an error probability of
/// exactly one, on a read that still counts. `ln 1 = 0`, so such a read adds nothing to the
/// log sum and a whole unit to the probability sum. Measured: 9 read-positions in 1,000 on
/// HG002 carry an error of exactly one, and they are **73% of Σ ε**; on tomato, 7 in 1,000 and
/// 47%.
///
/// # The depth cap, and what it costs to leave it
///
/// **The fitted rate and this denominator are not fitted over quite the same reads.** The
/// error-rate histogram thins every position to at most 124 reads before fitting; the
/// accumulator thins nothing. Per site the cap is harmless — the draw is on counts and never
/// looks at a read's quality — but across sites it re-weights, because a 500-read position
/// casts 500 votes in the denominator and 124 in the population the numerator came from, and
/// deep positions are not a random sample of the genome.
///
/// **Measured, on the deepest real sample there is:** on HG002's benchmark regions at 300×,
/// where the cap fires nearly everywhere, thinning the denominator to the same cap moves it
/// from 2.9055 × 10⁻⁴ to 2.9862 × 10⁻⁴ — **2.7%, which is 0.12 Phred**. On the tomato cohort
/// it moves it by a factor of 1.0000: 228,468,065 of 228,492,796 read-positions on the
/// deepest accession are under the cap.
///
/// **Decided by the owner on 2026-08-24: the denominator stays unthinned, and the 2.7% is
/// carried knowingly.** The scale is applied at calling time to *every* read the caller sees,
/// not to a subsample, so the average it is built from is over every read too. **Both options
/// are 2.7% wrong about something** — thinning would make spec §3.2's one-site-set requirement
/// literally true at the price of making the property that requirement exists to serve wrong by
/// the same amount — and this one is wrong about how the *fit* weights deep sites against
/// shallow ones, which is the fit's question rather than a second place that cancels it.
///
/// Below about 124 reads a position the two answers are identical, so on tomato-like data the
/// choice is free. Revisiting it is one multiply per site and a re-run.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct ReadGroupCalibration {
    /// The multiplier, applied to each read's own error probability. **One** where nothing
    /// was fitted, which leaves the qualities exactly as the instrument reported them.
    pub scale: f64,
    /// Where the scale came from. **`Defaulted` is not an error condition and it is not a
    /// detail**: a run calibrated against a measurement and a run trusting the instrument
    /// are otherwise indistinguishable in the output, and spec §3.2 requires them not to be.
    pub provenance: Provenance,
}

impl ReadGroupCalibration {
    /// The calibration for a read group the parameter fit measured, from its fitted rate and
    /// the two sums the accumulator kept.
    ///
    /// **`None` for any of three reasons, and they are one answer: there is no scale to be
    /// had.** The accumulator saw no read, so there is no average to divide by; or that
    /// average underflowed to zero, which no real base quality reaches but a hand-built
    /// fixture can; or the fitted rate itself is zero, which would make the scale zero and
    /// charge every read of the library the floor — maximal confidence about every base,
    /// from a number that says the fit found no errors at all. Take
    /// [`defaulted`](Self::defaulted) in all three, which is the honest answer and says so in
    /// the provenance.
    ///
    /// **A zero fitted rate is refused here rather than trusted**, on the same reasoning
    /// [`ContaminationView`] uses for a fraction that could not be identified: *absent* and
    /// *zero* are different answers, and only one of them is safe to multiply by.
    ///
    /// The accumulator is `parameter_estimation::generic::calibration`'s, which sums the
    /// per-read log error in fixed point precisely so that merging shards in different
    /// orders gives the same denominator — the same determinism requirement the row itself
    /// works under (spec §8).
    ///
    /// **The calibrated property holds to the accumulator's quantum, not exactly, and the
    /// difference is the fixed point.** In real arithmetic, scaling every read by one
    /// multiplier multiplies their geometric mean by it, so the calibrated average *is* the
    /// fitted rate with nothing to bound. The denominator that reaches this function has been
    /// rounded to units of 2⁻²⁰ nats, and [`MintedReadErrors`]'s own documentation bounds the
    /// resulting miss on the **mean** at 2⁻²¹ ≈ 4.768 × 10⁻⁷ nats — a bound that is attained
    /// rather than approached, and that does not grow with the length of the run. A relative
    /// shift of that size in the mean log error is a relative shift of that size in the rate,
    /// so the average the model charges sits within **five parts in ten million** of the
    /// fitted rate. Measured at 4 in a thousand over three reads spanning Phred 40 to Phred
    /// 13: **4.7 × 10⁻¹⁰ against a rate of 4 × 10⁻³**, a quarter of the bound.
    /// `tests::the_scale_makes_the_charged_average_the_fitted_rate` asserts against the bound
    /// rather than against zero, and names where it comes from.
    #[must_use]
    pub fn from_fitted_rate(rate: &Estimate<ErrorRate>, minted: MintedReadErrors) -> Option<Self> {
        let mean_minted_error = minted.mean_error_probability()?;
        (mean_minted_error > 0.0 && rate.value.get() > 0.0).then(|| Self {
            scale: rate.value.get() / mean_minted_error,
            // **The rate's own provenance and not `FittedHere`.** A rate borrowed from a
            // sibling read group makes a *borrowed* calibration, and stamping this one
            // `FittedHere` would launder it — which is the failure `Provenance`'s own
            // documentation says it exists to prevent. The scale adds no warrant of its own:
            // it is a ratio, so it is exactly as well founded as its numerator.
            provenance: rate.provenance,
        })
    }

    /// A read group the parameter fit emitted no rate for: **the qualities are used as
    /// reported**, and the provenance says so.
    ///
    /// **The case is too little data to fit from, not a sample the fit found difficult.** Two
    /// of five real tomato alignments ask for a noisier class than the pre-pass's model
    /// covers and are refused it — and they still get an error rate, the one-rate answer, so
    /// their calibration is `FittedHere` and not this. *Spec §3.2 puts those two samples on
    /// this side of the line and it is wrong to; the correction is owed there.*
    ///
    /// What must not happen is that a defaulted calibration be invisible: a run calibrated
    /// against a measurement and a run trusting the instrument are otherwise indistinguishable
    /// in the output (spec §3.2).
    #[must_use]
    pub fn defaulted() -> Self {
        Self {
            scale: 1.0,
            provenance: Provenance::Defaulted,
        }
    }

    /// The scale in log space — what spec §3.3's log-space form adds, once per observation.
    ///
    /// §3.3 charges an unexplained observation `q_sum + n·(log scale − log m)`, so the
    /// multiplier never appears as a multiplier there: it is one addition per observation.
    ///
    /// **The row no longer takes this route** (C1, 2026-08-24). Spec §3.6's mixture is
    /// evaluated in probability space, so the row multiplies by the scale rather than adding
    /// its logarithm — through [`charged_error`](Self::charged_error). What still reads this
    /// is the independent §3.3 oracle the row's own tests check the mixture against at zero
    /// contamination, which is the one place the log-space form has to be spelled out.
    ///
    /// **Exactly zero for a defaulted calibration**, since `ln 1 = 0` — so a read group the
    /// pre-pass emitted no rate for is charged exactly what its reads were minted with, with no
    /// arithmetic in between that could round.
    #[must_use]
    pub fn log_scale(&self) -> f64 {
        self.scale.ln()
    }

    /// **The error probability the row actually charges an observation** — this read group's
    /// scale times the geometric mean of its reads' minted error, floored against zero and
    /// **not capped from above**.
    ///
    /// # Why it floors and does not cap
    ///
    /// Production pairs its floor with a ceiling at a half, and ng adopted both at A2 before
    /// anything consumed them. **The ceiling cannot be applied to what the row charges**,
    /// and the reason is a stated requirement rather than a preference: spec §2.3 asks that
    /// no term be a non-linear function of a per-read quality, because the merge hands the
    /// row a *fold* of reads and `q_sum` recovers only their geometric mean. `min(x, ½)` is
    /// exactly such a function. **Measured on the row's own aggregation fixture**, which
    /// alternates reads between Phred 93 and Phred 1: a Phred-1 read is minted at an error of
    /// 0.794, so the ceiling replaces it with 0.5 and charges that read 0.46 nats less. It
    /// binds on every one of those reads taken singly and on none of the folds — a fold's
    /// geometric mean over the alternation is `2 × 10⁻⁵` — so at the 300-read end of that
    /// fixture, 150 such reads, pooling would move the answer by 69 nats where the property
    /// being pinned is agreement to a relative 2 × 10⁻¹⁴.
    ///
    /// **The floor is different in kind, which is why it stays.** It is not reached by any
    /// quality the read preparation admits — Phred 93 is an error of `5.0 × 10⁻¹⁰`, and at the
    /// smallest scale the row's fixtures use, 0.37, that is `1.9 × 10⁻¹⁰`, 185 times the floor
    /// — so it changes no answer that can occur and exists only so a logarithm never sees a
    /// zero (spec §8).
    ///
    /// **So this may return a number above one**, where a poor read is scaled up. Spec §3.3's
    /// log-space form does the same thing: it charges `q_sum + n·log scale`, which is
    /// positive under exactly the same conditions. That is a property of the model as
    /// specified, not something introduced here, and a row that clamped it would silently
    /// disagree with §3.3 wherever the clamp bit.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on three caller mistakes, and the first two are
    /// checked here because `f64::max` cannot be trusted to carry them (spec §8 puts a
    /// parameter outside its declared range on the assertion side):
    ///
    /// - **`num_reads == 0`.** An observation with no reads behind it is a row the merge does
    ///   not build, and it has no average error to take.
    /// - **A `q_sum` that is not a finite number at or below zero.** `q_sum` is a sum of
    ///   logarithms of probabilities, so it cannot be positive; a positive one returns a
    ///   *positive* log-likelihood, measured at `+48.90` on the row's own fixture. And
    ///   `x.max(y)` returns `y` when `x` is `NaN`, so **a `NaN` would come back as
    ///   `MIN_BASE_ERROR` — the most confident error probability this module admits** rather
    ///   than as a `NaN` anything downstream could see. That is the one place this differs
    ///   from [`charged_error`](Self::charged_error), whose
    ///   `f64::clamp` passes a `NaN` through, and it is why the check is here and not there.
    /// - **A scale that is not finite and positive.** [`scale`](Self::scale) is a public
    ///   field, so [`from_fitted_rate`](Self::from_fitted_rate)'s own guard can be bypassed
    ///   by building the struct literally, which the tests in this file do.
    #[must_use]
    pub fn charged_error(&self, q_sum: f64, num_reads: u32) -> f64 {
        assert!(
            num_reads > 0,
            "an observation with no reads behind it has no average error; the merge builds \
             no such row"
        );
        assert!(
            q_sum.is_finite() && q_sum <= 0.0,
            "an observation's summed log error is {q_sum}, and a sum of logarithms of \
             probabilities is a finite number at or below zero"
        );
        assert!(
            self.scale.is_finite() && self.scale > 0.0,
            "this read group's calibration scale is {}, and a scale is a finite positive \
             multiplier",
            self.scale
        );
        let mean_log_error = q_sum / f64::from(num_reads);
        (self.scale * mean_log_error.exp()).max(MIN_BASE_ERROR)
    }
}

/// §3.6's mixture inputs for one read group, frozen before the caller's loop starts.
///
/// **What is *not* here is the other half of the mixture.** The contaminating population's
/// frequency for the allele an observation shows is a property of the locus being called,
/// re-read from the caller's own estimate every iteration, so it is looked up where the
/// frequency is rather than frozen on this view (spec §3.6, corrected 2026-08-24; an earlier
/// design had three allele-class frequencies here and they are deleted).
///
/// **The two counts are not diagnostics.** A read group with too little evidence comes back
/// with a fraction near zero, because the likelihood barely moves with the fraction and the
/// search keeps zero — which is the right default for a value this term multiplies. A read
/// group that was measured and found clean comes back near zero as well. Those are different
/// claims, and the counts are the only thing that tells them apart, so they travel beside the
/// value rather than being summarised away.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct ContaminationView {
    /// `c` — the share of this read group's reads that came from another individual.
    pub fraction: f64,
    /// How many of the panel's varying positions this read group put a read on.
    pub markers_with_reads: u64,
    /// How many reads it put there. Both names are [`ContaminationEstimate`]'s own.
    pub reads_on_markers: u64,
    /// **Whose reads the fraction was fitted from**, and the third thing spec §3.6 says must
    /// travel beside it: `c` can be four different things and a consumer cannot act on them
    /// alike. A fraction fitted from this library's own reads and one fitted from every read
    /// of the plant and copied onto it are different claims — the first can say two libraries
    /// of one sample differ, the second cannot — and nothing downstream can tell them apart
    /// from the value.
    ///
    /// A plant sequenced from one library gets the same number either way, and that is every
    /// sample of every benchmark cohort here; the distinction bites on the plant with two.
    pub source: ContaminationSource,
}

impl ContaminationView {
    /// What the parameter fit's estimate says, or `None` where it identified nothing.
    ///
    /// **`None` is *absent*, not a fitted zero**, and the distinction is the reason the
    /// estimator returns an enum rather than an `Option<f64>`. At one sample there is no
    /// panel to be surprised by, so contamination is not estimable at all and the plain
    /// formula of spec §3.3 is what runs — which is the *simple* case for this model, not
    /// the weak one.
    #[must_use]
    pub fn of_estimate(estimate: &ContaminationEstimate) -> Option<Self> {
        match estimate {
            ContaminationEstimate::Estimated {
                alpha,
                source,
                markers_with_reads,
                reads_on_markers,
                ..
            } => Some(Self {
                fraction: *alpha,
                markers_with_reads: *markers_with_reads,
                reads_on_markers: *reads_on_markers,
                source: *source,
            }),
            ContaminationEstimate::NotIdentified { .. } => None,
        }
    }

    /// Whether anything at all stood behind this fraction.
    ///
    /// **This is the whole of what the counts can be asked without inventing a threshold.**
    /// A read group that touched no marker was not measured, whatever its fraction says. How
    /// *many* markers a fraction needs before it is a measurement rather than a shrug is a
    /// number nobody here has, so this predicate does not pretend to one: a consumer wanting
    /// a stronger test reads the counts and says what it chose.
    #[must_use]
    pub fn was_measured(&self) -> bool {
        self.markers_with_reads > 0 && self.reads_on_markers > 0
    }
}

/// Floor under a contaminant allele frequency, `1e-12`.
///
/// **It is defensive, not statistical, and the difference decides its size.** A batch that
/// never shows an allele is not proof the contaminating population lacks it — with 63 diploid
/// samples an unseen allele could still sit near 1 in 40 — but **expressing that uncertainty
/// here would be the wrong direction to be wrong in.** Spec §3.6 names the one failure to
/// watch: contamination attributed too readily suppresses real heterozygotes, by explaining
/// their alternative reads as somebody else's. A floor large enough to encode "we might just
/// have missed it" does exactly that at every candidate the cohort is thin on, which is most
/// of them.
///
/// So this floor's whole job is to keep a frequency strictly positive, and it is set where it
/// **cannot compete with the route it sits beside**: at a 3% contamination fraction it
/// contributes `3 × 10⁻¹⁴` against a misread's `ε̄/m`, which is `3 × 10⁻⁵` even at a middling
/// Phred 40 — **nine orders of magnitude below**.
/// `tests::the_frequency_floor_cannot_outweigh_a_misread` pins it.
///
/// **Settled by the owner, 2026-08-24: keep it very low.** The alternative was a statistical
/// reading — a pseudocount over the batch's copies, which for a 63-sample cohort would put an
/// unseen allele near 1 in 42 rather than 1 in a trillion. That is the honest expression of what
/// a finite batch does not know, and it is the wrong thing to put here, because it would say
/// *this read might well be contamination* at every candidate the cohort is thin on. **A floor
/// this low cannot compete with the route beside it, and that is its whole job.**
pub const MIN_CONTAMINANT_FREQUENCY: f64 = 1e-12;

/// Fill one locus's expected allele copies, added up per sequencing batch — the first half of
/// spec §3.6's `q(o)`, and the half that does not depend on which sample is being scored.
///
/// **Split from the frequency because the two have different grains.** The copies a batch holds
/// are a property of the locus, so they are summed once per locus per iteration; the frequency a
/// given sample's contaminant is drawn against leaves that sample's own copies out, so it is a
/// property of `(locus, sample)` and is filled by
/// [`fill_contaminant_allele_frequencies`] below. This is the same split the genotype prior
/// makes between the cohort's copies and the leave-one-out concentration it hands each sample.
///
/// `expected_copies_by_sample` is sample-major, `allele_count` to a sample, in the run's sample
/// order — the shape the loop holds. `out` is batch-major, `allele_count` to a batch.
///
/// # Panics
///
/// **In release as well as debug**, on shapes that do not agree; on a sample in a batch the run
/// did not declare; on a batch no sample is in — **which is not a lookup reading rubbish but
/// something worse**, since the row would never be written and would leave through the
/// frequency's no-evidence fallback, indistinguishable from a batch that really was sequenced
/// and really showed nothing; and on a copy count that is not finite and at or above zero, the
/// check [`ExpectedAlleleCopies`](crate::ng::calling::ExpectedAlleleCopies) makes for the same
/// numbers, restated here because this is the other door they come through.
pub fn fill_batch_allele_copies(
    expected_copies_by_sample: &[f64],
    batch_of_each_sample: BatchOfEachSample<'_>,
    allele_count: usize,
    out: &mut [f64],
) {
    let BatchOfEachSample(batch_of_each_sample) = batch_of_each_sample;
    assert!(
        allele_count > 0,
        "a locus is called over at least its reference allele"
    );
    assert_eq!(
        expected_copies_by_sample.len(),
        batch_of_each_sample.len() * allele_count,
        "the copies are {} entries for {} samples at {allele_count} alleles each, which is not \
         a whole cohort",
        expected_copies_by_sample.len(),
        batch_of_each_sample.len()
    );
    assert_eq!(
        out.len() % allele_count,
        0,
        "the copy table is {} entries and a batch's row is {allele_count} alleles wide",
        out.len()
    );
    let batch_count = out.len() / allele_count;
    assert!(
        batch_count > 0,
        "a run has at least one sequencing batch, the default being one that holds all of it"
    );
    for batch in 0..batch_count {
        assert!(
            batch_of_each_sample
                .iter()
                .any(|declared| declared.get() as usize == batch),
            "the copy table holds {batch_count} batches and no sample ran in batch {batch}; a \
             batch nobody was sequenced in is a batching that does not describe this run"
        );
    }

    out.fill(0.0);
    for (sample, (batch, own_copies)) in batch_of_each_sample
        .iter()
        .zip(expected_copies_by_sample.chunks_exact(allele_count))
        .enumerate()
    {
        let batch = batch.get() as usize;
        assert!(
            batch < batch_count,
            "sample {sample} says it ran in batch {batch}, and the run declared {batch_count}"
        );
        let batch_row = &mut out[batch * allele_count..(batch + 1) * allele_count];
        for (slot, &copies) in batch_row.iter_mut().zip(own_copies) {
            assert!(
                copies.is_finite() && copies >= 0.0,
                "sample {sample} is expected to carry {copies} copies of an allele, and a copy \
                 count is finite and at or above zero"
            );
            *slot += copies;
        }
    }
}

/// Fill one sample's contaminant allele frequencies, per sequencing batch — spec §3.6's `q(o)`,
/// **with this sample's own copies left out of its own batch**.
///
/// **This is a lookup and not a fit.** The caller's loop already estimates how many copies of
/// each allele every sample carries; a batch's frequency is those copies added up over the
/// samples that ran together, divided by the copies they hold in total. Nothing is estimated
/// here that the loop has not estimated already, which is why it is recomputed every iteration
/// rather than frozen beside the contamination fraction (spec §3.6, corrected 2026-08-24).
///
/// # Why the sample leaves itself out
///
/// **A contaminating read is somebody else's, by definition**, so the population it is drawn
/// against must not include the individual being scored. Without the subtraction a sample's own
/// alternative reads are partly explained by its own alternative allele — it is its own
/// contaminant — and the effect scales with how much of its batch it is: one part in 63 in a
/// tomato-sized cohort, and **all of it for a sample alone in its batch**, where an alternative
/// homozygote would return `q(alt) = 1`.
///
/// The genotype prior faced this exactly and answered it the same way: `fill_sample_concentration`
/// subtracts a sample's own expected copies before using the cohort's as that sample's prior.
/// This is that subtraction, one axis over (owner, 2026-08-24; C2's review found the batch of one).
///
/// **A sample alone in its batch therefore has nothing left**, and takes the no-evidence
/// fallback below — the reference, and the floor elsewhere. That is the conservative answer and
/// the right one: a library with no neighbours has no contaminating population to speak of.
///
/// # What a batch with nothing left gets, and why the count comes back
///
/// A batch whose remaining samples hold no copies — no coverage at this locus, or only the
/// sample being scored — cannot have a frequency, and dividing by its zero total would give
/// `NaN`. It gets the reference allele and nothing else: **the honest statement of what a batch
/// with no evidence says about a contaminant is what the reference says**, and it keeps every row
/// a distribution. An allele the batch's other samples never show gets
/// [`MIN_CONTAMINANT_FREQUENCY`], never zero.
///
/// **The return is how many batches took that fallback**, because otherwise it is invisible: at
/// a locus with one candidate allele a batch full of evidence and a batch with none are
/// bit-identical, both `[1.0]`, and at any allele count the fallback row is a legal-looking
/// frequency. A caller that wants to say *this locus was scored against a contaminant nobody
/// measured* has no other way to know (C2's review). **The loop calls this and does not yet
/// report it**: what the run says about the contamination it used is the step after the loop's
/// wiring, and the call site binds this rather than dropping it so that a reader sees it is owed
/// rather than handled.
///
/// # Panics
///
/// **In release as well as debug**, on shapes that do not agree; on a batch the copy table has no
/// row for; and on a batch whose copies sum to something that is not finite. **In debug**, on a
/// sample holding materially more copies than the batch it is part of — the sample's own copies
/// are one addend of that sum, so the two count paths have gone out of step, which is
/// `fill_sample_concentration`'s own check for the same subtraction.
#[must_use]
pub fn fill_contaminant_allele_frequencies(
    batch_copies: &[f64],
    own_copies: SampleAlleleCopies<'_>,
    own_batch: BatchId,
    allele_count: usize,
    out: &mut [f64],
) -> usize {
    assert!(
        allele_count > 0,
        "a locus is called over at least its reference allele"
    );
    assert_eq!(
        own_copies.get().len(),
        allele_count,
        "this sample's copies cover {} alleles and the locus is called over {allele_count}",
        own_copies.get().len()
    );
    assert_eq!(
        batch_copies.len() % allele_count,
        0,
        "the copy table is {} entries and a batch's row is {allele_count} alleles wide",
        batch_copies.len()
    );
    assert_eq!(
        out.len(),
        batch_copies.len(),
        "the frequency table is {} entries and the copy table it is built from is {}",
        out.len(),
        batch_copies.len()
    );
    let batch_count = batch_copies.len() / allele_count;
    let own_batch = own_batch.get() as usize;
    assert!(
        own_batch < batch_count,
        "this sample ran in batch {own_batch}, and the copy table holds {batch_count} batches"
    );

    out.copy_from_slice(batch_copies);
    // **The subtraction, and only on this sample's own batch.** `max(0, ·)` for the reason the
    // prior gives: the difference is exact algebra and inexact arithmetic, so a batch holding
    // only this sample can come out a hair below zero.
    for (allele, (slot, &own)) in out[own_batch * allele_count..(own_batch + 1) * allele_count]
        .iter_mut()
        .zip(own_copies.get())
        .enumerate()
    {
        let leaving_this_sample_out = *slot - own;
        debug_assert!(
            leaving_this_sample_out > COUNT_PATH_DESYNC_THRESHOLD,
            "this sample is expected to carry {own} copies of allele {allele} and its whole \
             batch {slot}; the sample's own copies are one addend of the batch's, so the two \
             count paths have gone out of step"
        );
        *slot = leaving_this_sample_out.max(0.0);
    }

    let mut batches_with_no_evidence = 0;
    for (batch, row) in out.chunks_exact_mut(allele_count).enumerate() {
        let total: f64 = row.iter().sum();
        // **A total that is not finite is arithmetic that went wrong, not a deep batch.** Each
        // slot is checked as it accumulates, so a row of finite copies whose sum overflows
        // reaches here — and then every ratio is `finite / inf`, which is zero, and every
        // allele is lifted to the floor. The row would come back finite, plausible, and saying
        // the neighbours carry nothing (C2's review, which built the case: one sample holding
        // 1e308 copies of each of two alleles).
        assert!(
            total.is_finite(),
            "batch {batch} holds {total} copies in total, which is arithmetic that went wrong \
             upstream rather than a batch with a great many samples in it"
        );
        if total <= 0.0 {
            // Nothing left to read a frequency off. The reference is what a batch with no
            // evidence has to say — and for this sample's own batch, that is what a library
            // with no neighbours gets.
            batches_with_no_evidence += 1;
            row.fill(MIN_CONTAMINANT_FREQUENCY);
            row[usize::from(AlleleId::REFERENCE.get())] =
                1.0 - MIN_CONTAMINANT_FREQUENCY * (allele_count - 1) as f64;
            continue;
        }
        for slot in row.iter_mut() {
            *slot = (*slot / total).clamp(MIN_CONTAMINANT_FREQUENCY, 1.0);
        }
    }

    // A postcondition rather than a check on the caller: every route above either divides by a
    // total already asserted finite and positive, or fills the row with constants. It is kept
    // because it is the one statement of what this function promises, and it is `debug` only
    // because it sweeps the whole table once per sample to prove something no input can break.
    debug_assert!(
        out.chunks_exact(allele_count)
            .all(|row| row.iter().all(|f| f.is_finite())),
        "a contaminant frequency came out as an infinity or a `NaN`, which is arithmetic that \
         went wrong rather than a thin batch"
    );
    batches_with_no_evidence
}

/// **The half of spec §3.6's mixture that is frozen for the run** — a contamination fraction per
/// read group, and the sequencing batch each of those read groups ran in.
///
/// # Why it is a type and not two arguments
///
/// **Because the checks between them are worth making once and the loop would make them per
/// row.** With contamination on, the mixture's *other* half moves: `q(o)` is recomputed for
/// every sample of every pass, and each recomputation needs a fresh [`ContaminationMixture`]
/// over the refilled table. Checking the fractions against the batching there costs
/// `read groups × batches` per sample per pass.
///
/// **Both sides of that comparison, computed rather than asserted.** At a thousand libraries on
/// four plates the checks are about 5,000 operations, so once a sample once a pass over a
/// thousand samples and seven passes is about 3.5 × 10⁷ at one locus — against the row
/// assembly's own `samples × observations × genotypes × passes`, about 1.5 × 10⁶ at ten
/// observations and a six-allele diploid. **Roughly twenty times the work the assembly does.**
/// Held apart, the expensive checks are made once per assembly and
/// [`with_frequencies`](Self::with_frequencies) is what a row costs.
///
/// **What is still per row, and where it would bite.** `with_frequencies` range-checks the
/// table it is handed, which is `batches × alleles`, and the fill above it rewrites the whole
/// table per sample for the same reason. Under the shipped default that is one batch — an
/// allele's worth of work — and it stays small for a plate-sized batching. It becomes the
/// dominant cost only where a run declares roughly as many batches as it has samples, which
/// nothing produces today, and closing it means changing what
/// [`fill_contaminant_allele_frequencies`] writes rather than what this checks.
///
/// It also states in the type which half is which, which is the distinction §6.1's tier table
/// turns on: the fraction is a property of the *library* and frozen before calling; the
/// frequency is a property of the *locus* and moves with the loop.
#[derive(Copy, Clone, Debug)]
pub struct FrozenContamination<'a> {
    fractions: &'a [ContaminationView],
    batch_of_each_read_group: &'a [BatchId],
    batch_count: usize,
}

impl<'a> FrozenContamination<'a> {
    /// A fraction per read group, the batch each ran in, and how many batches the run declares.
    ///
    /// # Panics
    ///
    /// On a fraction outside `[0, 1)`; on the two slices disagreeing about how many read groups
    /// the run has; on a run of no batches or no fractions — the uncontaminated case is
    /// [`ContaminationMixture::uncontaminated`] and not an empty list; on a read group naming a
    /// batch the run did not declare; and on a declared batch no read group ran in.
    ///
    /// **The last of those is the one that fails silently.** A run that declared three batches
    /// against a batching naming only the first constructs cleanly otherwise: every read group
    /// reads row 0, `batch_count` reports three, and **a run that declared batches is scored
    /// against the cohort frequency with nothing said** (C2's review).
    #[must_use]
    pub fn new(
        fractions: &'a [ContaminationView],
        batch_of_each_read_group: BatchOfEachReadGroup<'a>,
        batch_count: usize,
    ) -> Self {
        let BatchOfEachReadGroup(batch_of_each_read_group) = batch_of_each_read_group;
        assert!(
            !fractions.is_empty(),
            "a mixture with no fractions and no frequencies is the uncontaminated case, and it \
             is spelled `ContaminationMixture::uncontaminated` — one named way to say it, so \
             that a caller reaches the decision rather than the shortest thing that compiles"
        );
        assert_eq!(
            batch_of_each_read_group.len(),
            fractions.len(),
            "every read group runs in exactly one batch, so the batching holds {} entries \
             against {} read-group fractions",
            batch_of_each_read_group.len(),
            fractions.len()
        );
        assert!(
            batch_count > 0,
            "a contaminated run has at least one sequencing batch, the default being one that \
             holds all of it"
        );
        for (read_group, batch) in batch_of_each_read_group.iter().enumerate() {
            assert!(
                (batch.get() as usize) < batch_count,
                "read group {read_group} says it ran in batch {batch}, and this mixture holds \
                 {batch_count} batches"
            );
        }
        for batch in 0..batch_count {
            assert!(
                batch_of_each_read_group
                    .iter()
                    .any(|declared| declared.get() as usize == batch),
                "the run declares {batch_count} batches and no read group ran in batch \
                 {batch}; a batching that names fewer batches than the frequency table has rows \
                 would score every library against the first"
            );
        }
        for (read_group, view) in fractions.iter().enumerate() {
            assert!(
                (0.0..1.0).contains(&view.fraction),
                "read group {read_group} has a contamination fraction of {}, and a fraction is \
                 a share of that group's reads that came from somebody else: a number at or \
                 above zero, below one, and not a `NaN`. A whole library of another \
                 individual's DNA is not a sample of this one",
                view.fraction
            );
        }
        Self {
            fractions,
            batch_of_each_read_group,
            batch_count,
        }
    }

    /// This half, plus one locus's contaminant allele frequencies — the whole mixture a row
    /// reads.
    ///
    /// **What a row costs**, and deliberately: the checks here are the ones that depend on the
    /// table just filled, and nothing that was already settled at construction is made again.
    ///
    /// # Panics
    ///
    /// On a frequency table that is not `batch_count × allele_count` entries, and on an entry
    /// outside `[0, 1]`.
    #[must_use]
    pub fn with_frequencies(
        self,
        contaminant_allele_frequencies: &'a [f64],
        allele_count: usize,
    ) -> ContaminationMixture<'a> {
        assert!(
            allele_count > 0,
            "a locus is called over at least its reference allele, so a mixture cannot be for \
             none"
        );
        assert_eq!(
            contaminant_allele_frequencies.len(),
            self.batch_count * allele_count,
            "the frequency table holds {} entries and the run declares {} batches over \
             {allele_count} alleles",
            contaminant_allele_frequencies.len(),
            self.batch_count
        );
        for (slot, &frequency) in contaminant_allele_frequencies.iter().enumerate() {
            assert!(
                (0.0..=1.0).contains(&frequency),
                "the contaminating population shows allele {} of batch {} at a frequency of \
                 {frequency}, which is not a frequency",
                slot % allele_count,
                slot / allele_count
            );
        }
        ContaminationMixture {
            fractions: self.fractions,
            batch_of_each_read_group: self.batch_of_each_read_group,
            contaminant_allele_frequencies,
            allele_count,
        }
    }

    /// How many sequencing batches the run declares.
    #[must_use]
    pub fn batch_count(&self) -> usize {
        self.batch_count
    }
}

/// **Both halves of spec §3.6's mixture, for one locus** — how much of each read group came
/// from somebody else, and what that somebody else's DNA would have shown here.
///
/// The mixture the row evaluates is
///
/// ```text
///     n_o · log[ (1 − c) · own(o | g)  +  c · q(o) ]
/// ```
///
/// and the two halves reach the row from different places, which is the reason they are
/// gathered into one type rather than passed as two loose slices. **`c` is a property of the
/// library**: fitted before calling and frozen for the whole run, one per read group. **`q(o)`
/// is a property of the locus**: the frequency, among the samples sequenced in the same batch,
/// of the allele this observation shows — the caller's own estimate, which moves every
/// iteration (spec §3.6, corrected 2026-08-24; §6.1's first tier holds the fraction only).
/// Holding them together is what lets one construction check that they describe the same run,
/// instead of the row rediscovering it per observation.
///
/// **The two halves are also indexed differently, and that is the shape to know.** The
/// fraction is one number per read group. The frequency is one *row* per sequencing batch,
/// `allele_count` wide, and which row a read group reads is what the batching says — because
/// the population a contaminant is drawn from is the samples that ran beside it, not the
/// cohort (spec §3.6). **The default batching is one batch holding the run**, so the table is
/// a single row and every read group reads the cohort frequency; a run that declares no
/// batching loses nothing it had, and nothing here branches on the batching's absence.
/// [`fill_batch_allele_copies`] sums each batch's copies once per locus per iteration, and
/// [`fill_contaminant_allele_frequencies`] turns them into this table **once per sample**,
/// leaving that sample's own copies out of its own batch — a contaminant is somebody else, so
/// the individual being scored is not in the population its contaminant is drawn from.
///
/// **What the default does cost is that it leaves no trace**, and that belongs beside the
/// claim that it loses nothing. A defaulted `[ONLY, ONLY]` and a run that genuinely declared
/// one batch holding two read groups are *the same value*, so nothing downstream can tell a
/// batching that was stated from one that was assumed — where the other half of this mixture
/// carries [`ContaminationSource`] for exactly that reason. The distinction belongs to
/// [`SequencingBatches::is_default`](crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches::is_default),
/// which answers it — built at the calling loop's E2a, and reachable from the run's parameters
/// through `RunParameters::sequencing_batches`. **Nothing reports it yet**, which is the step
/// after this one.
///
/// **[`uncontaminated`](Self::uncontaminated) is the uncontaminated case and is not a special path.** With no
/// fractions the mixture collapses to `own(o | g)` and the row computes spec §3.3 — to a few
/// units in the last place, because §8 evaluates in probability space and §3.3 in log space, so
/// there is an `exp`/`log` round trip between the two forms. That is why the row has **no
/// `c == 0` branch**: the two are the same algebra, and production keeps such a branch only
/// because its own two forms genuinely differ (spec §3.6).
///
/// **At one sample there is no mixture at all.** Contamination is a comparison between samples,
/// so a single-sample run has nothing to estimate `c` from and gets [`uncontaminated`](Self::uncontaminated) — *emit
/// it as absent*, not a fitted zero.
#[derive(Copy, Clone, Debug)]
pub struct ContaminationMixture<'a> {
    fractions: &'a [ContaminationView],
    batch_of_each_read_group: &'a [BatchId],
    /// Batch-major, one row of `allele_count` per batch — the same shape, and for the same
    /// reason, as the error-spread table the row reads beside it.
    contaminant_allele_frequencies: &'a [f64],
    allele_count: usize,
}

impl<'a> ContaminationMixture<'a> {
    /// Nothing is contaminated, so the row computes spec §3.3.
    ///
    /// This is what a run whose parameter fit emitted no fraction gets, and what a
    /// single-sample run gets.
    #[must_use]
    pub fn uncontaminated() -> Self {
        Self {
            fractions: &[],
            batch_of_each_read_group: &[],
            contaminant_allele_frequencies: &[],
            allele_count: 0,
        }
    }

    /// A fraction per read group, the batch each read group ran in, and one frequency per
    /// candidate allele **per batch** — batch-major, `allele_count` to a row.
    ///
    /// **The batch is the population a contaminant is drawn from**, which is why the
    /// frequencies have a second axis and the fractions do not. A contaminating read is far
    /// likelier to have come from a neighbour on the same run than from a random member of the
    /// species (spec §3.6), so `q(o)` is the allele's frequency among the samples sequenced
    /// beside this one. **The default batching puts every read group in [`BatchId::ALL_TOGETHER`]**, so
    /// the table is one row and every read group reads the cohort frequency — a run that
    /// declares no batching loses nothing it had, and nothing here branches on that.
    ///
    /// # Panics
    ///
    /// On a fraction outside `[0, 1)` or a frequency outside `[0, 1]`; on the three slices not
    /// agreeing about how many read groups, batches and alleles there are; and on any of them
    /// being empty while the others are not. **These are checked once here rather than per
    /// observation**, which is spec §8's rule for a parameter outside its declared range: the
    /// model has no failure mode of its own that is not a caller bug.
    ///
    /// **A fraction of exactly one is refused rather than floored**, and it is the one bound
    /// that does arithmetic rather than documentation. `own(o | g)` is positive — the copy
    /// share is at least `1/P` and the charged error is floored — so the mixture can only
    /// reach zero, and its logarithm negative infinity, when `1 − c` underflows to zero and
    /// the contaminant's frequency for that allele is zero as well. A library that is
    /// *entirely* another individual's DNA is not a sample of the individual it is labelled
    /// with, and no estimator here can return it.
    #[must_use]
    pub fn new(
        fractions: &'a [ContaminationView],
        batch_of_each_read_group: BatchOfEachReadGroup<'a>,
        contaminant_allele_frequencies: &'a [f64],
        allele_count: usize,
    ) -> Self {
        // **Both halves or neither, checked before anything else.** Each half on its own has a
        // named refusal below, and both of those say *this half is missing*; only the pair says
        // that the caller meant the uncontaminated case and half-spelled it.
        assert_eq!(
            fractions.is_empty(),
            contaminant_allele_frequencies.is_empty(),
            "a mixture needs both halves or neither: {} read-group fractions against {} allele \
             frequencies",
            fractions.len(),
            contaminant_allele_frequencies.len()
        );
        assert!(
            !fractions.is_empty(),
            "a mixture with no fractions and no frequencies is the uncontaminated case, and it \
             is spelled `ContaminationMixture::uncontaminated` — one named way to say it, so \
             that a caller reaches the decision rather than the shortest thing that compiles"
        );
        assert!(
            allele_count > 0,
            "a locus is called over at least its reference allele, so a mixture cannot be for \
             none"
        );
        assert_eq!(
            contaminant_allele_frequencies.len() % allele_count,
            0,
            "the frequency table is {} entries and a batch's row is {allele_count} alleles \
             wide, so it is not a whole number of batches",
            contaminant_allele_frequencies.len()
        );
        FrozenContamination::new(
            fractions,
            batch_of_each_read_group,
            contaminant_allele_frequencies.len() / allele_count,
        )
        .with_frequencies(contaminant_allele_frequencies, allele_count)
    }

    /// Whether this is the uncontaminated case, in which the row computes spec §3.3.
    ///
    /// **The branch that decides which formula runs deserves a name.** The row used to ask
    /// this as `allele_count() == 0`, which announces at the call site only that a table is
    /// empty (C1's review).
    #[must_use]
    pub fn is_absent(&self) -> bool {
        self.fractions.is_empty()
    }

    /// How many candidate alleles this mixture holds a contaminant frequency for — the count
    /// the row checks against the locus it is calling.
    #[must_use]
    pub fn allele_count(&self) -> usize {
        self.allele_count
    }

    /// How many sequencing batches the run declared — `1` under the default batching, and `0`
    /// where there is no mixture.
    #[must_use]
    pub fn batch_count(&self) -> usize {
        if self.is_absent() {
            return 0;
        }
        self.contaminant_allele_frequencies.len() / self.allele_count
    }

    /// How many read groups this mixture holds a fraction for — the count the row checks
    /// against the calibration the run supplied.
    ///
    /// **A mixture is dense over read groups**, because [`fraction_of`](Self::fraction_of)
    /// indexes it by [`ReadGroupId`]. Without this the row could only discover a mixture built
    /// for a different run *lazily*, when some observation happened to name a read group past
    /// the end — so a locus whose reads all came from the first few groups would pass, and the
    /// mismatch would surface at whichever locus first reached further, or never (C1's review).
    #[must_use]
    pub fn read_group_count(&self) -> usize {
        self.fractions.len()
    }

    /// `c` for one read group — zero where there is no mixture at all.
    ///
    /// # Panics
    ///
    /// On a read group the mixture has no entry for, once there is a mixture. A row scoring a
    /// read group the parameters do not cover is the same caller bug as a row scoring one the
    /// calibration does not cover, and silently reading it as clean would be a genotype quietly
    /// moved rather than a run that stopped.
    #[must_use]
    pub fn fraction_of(&self, read_group: ReadGroupId) -> f64 {
        if self.fractions.is_empty() {
            return 0.0;
        }
        let at = read_group.get() as usize;
        self.fractions
            .get(at)
            .unwrap_or_else(|| {
                panic!(
                    "read group {at} has no contamination fraction; the run supplied {}",
                    self.fractions.len()
                )
            })
            .fraction
    }

    /// `q(o)` — **how often the samples sequenced beside this read group carry the allele this
    /// observation showed**; zero where there is no mixture.
    ///
    /// **Both of its arguments do work.** The allele is what the read showed; the read group
    /// decides *whose* frequency is read, because the batching says which samples ran beside
    /// it. Under the default batching every read group lands on the same row and this is the
    /// cohort frequency.
    ///
    /// **Zero is unreachable from the producer and reachable in a test**, which is worth
    /// keeping apart: [`fill_contaminant_allele_frequencies`] floors an allele the batch's
    /// *other* samples never showed at [`MIN_CONTAMINANT_FREQUENCY`]. What this accessor
    /// promises is only that it returns what it was given.
    ///
    /// # Panics
    ///
    /// On a read group or an allele the mixture has no entry for, once there is a mixture. The
    /// row checks both counts before its observation walk, so reaching either from a run means
    /// the row and the mixture disagree about the locus.
    #[must_use]
    pub fn contaminant_frequency_of(&self, read_group: ReadGroupId, allele: AlleleId) -> f64 {
        if self.contaminant_allele_frequencies.is_empty() {
            return 0.0;
        }
        let group = read_group.get() as usize;
        let batch = self
            .batch_of_each_read_group
            .get(group)
            .unwrap_or_else(|| {
                panic!(
                    "read group {group} is in no sequencing batch; the run batched {}",
                    self.batch_of_each_read_group.len()
                )
            })
            .get() as usize;
        let allele = usize::from(allele.get());
        assert!(
            allele < self.allele_count,
            "allele {allele} has no contaminant frequency; this mixture is for \
             {} alleles",
            self.allele_count
        );
        self.contaminant_allele_frequencies[batch * self.allele_count + allele]
    }
}

/// The merge's fold of every read that showed one allele from one read group — one term of
/// the SNP/indel formula's first two sums.
///
/// **Four numbers where the merge's own row carries eight** — [`SupportedAllele`]'s allele
/// and read group, and the six of [`AlleleSupport`] — and the four that are dropped are
/// dropped on purpose: the forward-strand count, the two mapping-quality moments and the
/// read-position count belong to the site filters that run *after* genotyping, and none of
/// them enters a likelihood (spec §1.4). What is kept is what the formula reads, and the
/// type is what stops a filter statistic reaching it.
///
/// **`q_sum` is a sum of logarithms, not a probability.** It is Σ `ln P(this read is
/// wrong)` over the reads folded here, which is exactly the shape the formula needs: a read
/// the genotype cannot explain is charged its own log error, so summing the logs and
/// charging the sum are the same arithmetic as charging each read separately. That is what
/// makes the likelihood of a list of reads and the likelihood of the merge's fold of those
/// same reads agree to the last bit.
///
/// `Copy`, because it is four scalars the row builder indexes in an inner loop and nothing
/// it owns needs dropping. `tests::the_observation_stays_four_scalars_wide` pins the size, so
/// that a field which would make copying expensive — anything owning heap — has to be argued
/// for rather than arriving unnoticed.
///
/// [`AlleleSupport`]: crate::ng::run::cohort_merge::build::AlleleSupport
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct GenericObservation {
    /// Which allele these reads showed — `a(o)` in the formula, an index into the locus's
    /// [`CandidateAlleles`](super::CandidateAlleles).
    ///
    /// **A candidate index and not the merge's own**, which is why nothing here can mint
    /// one: see the module's *Two allele tables* note.
    pub allele: AlleleId,
    /// Which read group they came from — `r(o)`, and **part of the identity**, never a
    /// label to fold away. It selects the calibration scale and the contamination fraction
    /// this term is charged under (spec §2.3).
    pub read_group: ReadGroupId,
    /// How many reads showed it — `n_o`.
    pub num_reads: u32,
    /// Σ `ln P(error)` over those reads — `q_sum_o`, straight off the merge's row.
    pub q_sum: f64,
    /// Of those reads, how many were on the forward strand.
    ///
    /// **The read model does not use this, and it is carried anyway** — the same reason
    /// [`Self::q_sum`]'s own comment gives for keeping a term that cancels. The site
    /// quality's artifact correction compares the strand split of the reference allele's
    /// reads against the alternative's, and the evidence it compares is released with the
    /// locus (`doc/devel/ng/spec/calling_quality.md` §3.3). This view is where that
    /// specification says the counts travel, and the merge's `AlleleSupport` already has
    /// them; anything that gathered them a second way would be a second answer to the same
    /// question.
    pub forward_reads: u32,
    /// Of those reads, how many started strictly left of the record they were seen at —
    /// freebayes' `placedLeft`, and the read-position half of the same correction.
    ///
    /// **Counted against each record's own position, not the cohort locus's first base**, by
    /// the mint that produced it; nothing re-anchors it, so where a locus spans several of a
    /// sample's records this mixes as many questions as it has records. Carried forward from
    /// [`AlleleSupport::placed_left`](crate::ng::run::cohort_merge::AlleleSupport::placed_left),
    /// whose own comment is the one to read before drawing a conclusion from it.
    pub placed_left_reads: u32,
}

impl GenericObservation {
    /// The six numbers of one merge row, under the candidate id selection gave that row's
    /// allele.
    ///
    /// **The id is an argument because this module cannot compute it.** The row's own
    /// `allele` field indexes the merge's unification table and [`AlleleId`] indexes the
    /// candidate table, and only selection knows how one maps onto the other (see the
    /// module's *Two allele tables* note).
    ///
    /// **Four of the six are the read model's and two are the site quality's**
    /// ([`Self::forward_reads`], [`Self::placed_left_reads`]). They are copied here rather
    /// than fetched later because the merge row they come from is released with the locus.
    #[must_use]
    pub fn of_supported_allele(row: &SupportedAllele, allele: AlleleId) -> Self {
        Self {
            allele,
            read_group: row.read_group,
            num_reads: row.support.num_reads,
            q_sum: row.support.q_sum,
            forward_reads: row.support.num_fwd,
            placed_left_reads: row.support.placed_left,
        }
    }

    /// Narrow one sample's merge rows into `out` — cleared first — under selection's
    /// mapping, and return the pooled quality of the rows it dropped.
    ///
    /// `candidate_of_merge_allele` is indexed by the merge's own allele index and says which
    /// candidate that allele became, or `None` where selection dropped it.
    ///
    /// **The return value is not the whole of [`GenericSampleEvidence::unmatched_q_sum`],
    /// and calling it that would be wrong.** It is the part this call discarded — Σ `q_sum`
    /// over the rows whose allele has no candidate. Selection pools its own leftovers too,
    /// and the specification that says how is not written; a caller adds this to whatever
    /// selection hands it. What matters here is that a function that drops rows says what it
    /// dropped, because those reads are evidence and the data likelihood is compared between
    /// loci (spec §3.3).
    ///
    /// **`out` is caller scratch, cleared and refilled**, so a worker holds one buffer
    /// across every sample of every locus and the row function still allocates nothing
    /// (spec §8). There is no borrow to be had instead: a merge row is 48 bytes and this one
    /// is 32 — it was 24 before the two site-quality counters joined it — so no
    /// reinterpretation of the merge's `Vec` exists and somebody has to fill a parallel
    /// buffer.
    ///
    /// **The output order is the input order**, so `out` is ascending on
    /// `(candidate allele, read group)` exactly when selection's mapping keeps alleles in
    /// order — which the renumbering does, since it drops alleles and slides the rest down.
    /// [`GenericSampleEvidence::new`] is where that promise is checked.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on a row whose merge allele index is past the end of
    /// `candidate_of_merge_allele`. That is a mapping built against a different locus's
    /// allele table, and the alternative to panicking is scoring reads against whatever
    /// allele happens to sit at that index (spec §8; `per_group_merger.rs:1963` is the
    /// precedent for holding a structural assertion in release).
    pub fn fill_from_supported_alleles(
        rows: &[SupportedAllele],
        candidate_of_merge_allele: &[Option<AlleleId>],
        out: &mut Vec<Self>,
    ) -> f64 {
        out.clear();
        let mut dropped_q_sum = 0.0;
        for row in rows {
            let mapped = candidate_of_merge_allele
                .get(row.allele)
                .unwrap_or_else(|| {
                    panic!(
                        "merge allele {} has no entry in a mapping of {} alleles, so the \
                         mapping was built against a different locus's table",
                        row.allele,
                        candidate_of_merge_allele.len()
                    )
                });
            match mapped {
                Some(allele) => out.push(Self::of_supported_allele(row, *allele)),
                None => dropped_q_sum += row.support.q_sum,
            }
        }
        dropped_q_sum
    }
}

/// One sample's evidence at one cohort locus, as the SNP/indel row consumes it.
///
/// **A view, and the two halves borrow different things.** `partials` borrows the merge's
/// own rows, which is why there is one [`PartialObservation`] type rather than two.
/// `supported` cannot: the merge's row is twice the width of the one the formula reads, so
/// it borrows a staging buffer some caller filled with
/// [`GenericObservation::fill_from_supported_alleles`] and holds across loci. Either way the
/// row function allocates nothing per sample (spec §8).
#[derive(Copy, Clone, Debug)]
pub struct GenericSampleEvidence<'a> {
    /// One entry per `(allele, read group)` this sample's **complete** reads showed, in the
    /// ascending pair order the merge builds them in.
    ///
    /// **The order is a determinism requirement, not tidiness.** Floating-point addition is
    /// not associative, so two runs that summed a sample's observations in different orders
    /// would produce log-likelihoods differing in the last bits, and the run's genotypes
    /// with them (spec §8). The merge sorts on the pair once per sample, on the one path
    /// both its branches converge to (`AlleleTable::assemble` in
    /// `src/ng/run/cohort_merge/build.rs`), and its own
    /// `the_rows_are_ordered_by_allele_then_read_group` pins it against a fixture that
    /// arrives out of order — but that is a test in another module and this view is built
    /// from a buffer, so [`new`](Self::new) checks the promise here as well.
    pub supported: &'a [GenericObservation],
    /// Σ `ln P(error)` pooled over reads matching **no** candidate allele — the support of
    /// the merge's alleles that candidate selection dropped, folded here by that selection.
    ///
    /// **The same number for every genotype, so it cancels in genotyping, and it is kept
    /// anyway**: the data likelihood also feeds emission and QUAL, where an absolute value
    /// is compared between loci (spec §3.3's `q_sum_other`).
    ///
    /// **Nothing in this crate produces the whole of it yet.** It is created by candidate
    /// *selection*, which has no specification — "whoever specifies selection owes the
    /// pool" — so the row takes it as an input and tests supply it from fixtures.
    /// [`GenericObservation::fill_from_supported_alleles`] returns the part it can see, the
    /// quality of the rows the mapping dropped.
    pub unmatched_q_sum: f64,
    /// The partial observations, bases and witnessed positions intact.
    ///
    /// **Not folded onto alleles**, because there is no allele: a partial's bases cannot be
    /// compared against a whole-span allele. Padding one out to the span and interning it
    /// would put a sequence in the table that no molecule carried, and it would read as a
    /// *short* allele — the one direction the model must not be biased in. What a candidate
    /// is scored against is its projection *restricted to the positions the read witnessed*,
    /// so those positions have to survive to the scoring (spec §5.1, §5.3).
    ///
    /// **They are a set of runs with holes in it, not one run**, because the generic fold
    /// mints witnesses that way (spec §5.3, corrected 2026-08-24). So the restricted
    /// projection is a **gather** and cannot be a subslice of the allele's bases — it needs a
    /// buffer sized by the widest witness, which is the generic row's own scratch and lands
    /// with Milestone D.
    ///
    /// **The witnessed positions and the bases are on different axes, and their lengths are
    /// not interchangeable.** The witness counts *locus positions*; the bases are what the
    /// read showed over them, so the two differ by the net indel the read carried — a read
    /// carrying a two-base insertion and a two-base deletion inside the stretch comes back
    /// with as many bases as positions and is still not a positional match for any of them
    /// ([`PartialObservation::bases`]). Scoring indexes the *allele's* projection with the
    /// witness, never the partial's own bases.
    pub partials: &'a [PartialObservation],
}

impl<'a> GenericSampleEvidence<'a> {
    /// A sample's evidence from the three parts the merge and candidate selection hand
    /// over.
    ///
    /// A constructor rather than a struct literal at each call site so that a fourth part,
    /// when one arrives, cannot be forgotten at one of them.
    ///
    /// # Panics
    ///
    /// In debug builds, on `supported` rows that are not strictly ascending on
    /// `(allele, read group)`. The merge builds them that way and the determinism of the
    /// sum rests on it (spec §8), but this constructor takes any slice: a call site that
    /// concatenated two samples' rows, or filtered and re-spliced them, or applied a
    /// selection mapping that reordered the alleles, would break it with nothing saying so.
    /// Debug-only because the check is linear in the rows and the property is the merge's to
    /// hold, not this type's to enforce at every call.
    #[must_use]
    pub fn new(
        supported: &'a [GenericObservation],
        unmatched_q_sum: f64,
        partials: &'a [PartialObservation],
    ) -> Self {
        debug_assert!(
            supported
                .windows(2)
                .all(|pair| (pair[0].allele, pair[0].read_group)
                    < (pair[1].allele, pair[1].read_group)),
            "the merge builds one row per (allele, read group) in ascending pair order and \
             the determinism of the sum rests on it (spec §8); these are not ascending: \
             {supported:?}"
        );
        Self {
            supported,
            unmatched_q_sum,
            partials,
        }
    }

    /// A sample that showed nothing at this locus.
    ///
    /// **Not a special case in the row function**, and that is the point of naming it: an
    /// empty sum is zero, so every genotype scores zero and the prior decides alone, which
    /// is the right answer rather than a branch (spec §3.3).
    #[must_use]
    pub fn empty() -> Self {
        Self {
            supported: &[],
            unmatched_q_sum: 0.0,
            partials: &[],
        }
    }
}

/// One sample's evidence at one repeat tract, as the STR row consumes it.
///
/// **The locus generator's own type, unchanged.** [`SequenceObservation`] already keys on
/// `(bases, witness, read group)` and carries a read count, so the aggregation contract
/// holds by construction and there is nothing for this module to re-shape (arch §2.2).
#[derive(Copy, Clone, Debug)]
pub struct SsrSampleEvidence<'a> {
    /// Every observation at the tract, complete and partial alike, in the order the
    /// generator sorted them — which is the fixed order spec §8 requires the sum to run in.
    ///
    /// **Reaching this field directly is how a partial gets scored as though it were
    /// complete**, which mis-scores it as a *short* allele, because its bases are a prefix
    /// of the truth (spec §5.1). [`complete_observations`](Self::complete_observations) and
    /// [`partial_observations`](Self::partial_observations) split it by witness and carry
    /// each observation's position in this slice, which is the key an emission cache is
    /// built on — so the guarded route is also the useful one.
    ///
    /// **Calling them a guard would overstate them, and an earlier version of this comment
    /// did.** The field is public, so nothing stops a consumer walking it — and the
    /// generator's own `SampleLocusObservations::observations` is public too, so the
    /// iterator these are named after was never a guard either. What the two methods buy is
    /// that scoring a partial as complete has to be written rather than fallen into, and
    /// that the split has one spelling instead of one per caller.
    pub observations: &'a [SequenceObservation],
    /// The tract's motif and flanks — what the emission needs to know which repeat it is
    /// scoring against.
    pub detail: &'a SsrDetail,
}

impl<'a> SsrSampleEvidence<'a> {
    /// A sample's observations at one tract, against that tract's repeat detail.
    #[must_use]
    pub fn new(observations: &'a [SequenceObservation], detail: &'a SsrDetail) -> Self {
        Self {
            observations,
            detail,
        }
    }

    /// The observations whose reads spanned the whole tract — the ones an emission may
    /// score directly — each with its position in
    /// [`observations`](Self::observations).
    ///
    /// **The position travels because the STR row's emission cache is keyed by it.** That
    /// cache is what makes a row cost `observations × candidates` rather than
    /// `observations × genotypes`, which is the design and not an optimisation (spec §8);
    /// an iterator yielding only the observation would leave the row re-walking the
    /// unguarded field to recover the key.
    ///
    /// **It splits on the witness and never on the bases.** The same bases seen by a read
    /// that spanned the tract and by one that ran out are two observations, and only the
    /// witness separates them (spec §1.3) — so the plausible shortcut, *is this shorter
    /// than the tract?*, puts both on the same side and is wrong.
    ///
    /// The `use<'a>` is load-bearing: in edition 2024 an `impl Trait` return captures every
    /// lifetime in scope, so without it the iterator borrows `self` and cannot outlive the
    /// view — which the row builder needs it to.
    pub fn complete_observations(
        &self,
    ) -> impl Iterator<Item = (usize, &'a SequenceObservation)> + use<'a> {
        self.observations
            .iter()
            .enumerate()
            .filter(|(_, observation)| observation.read_witness == ReadWitness::Complete)
    }

    /// The observations whose reads ran out inside the tract — the ones the censored term
    /// scores, and nothing else may — each with its position in
    /// [`observations`](Self::observations).
    ///
    /// **This is not a small corner of the evidence.** A tract's locus is as wide as the
    /// tract, so over half the overlapping reads are partial at a 60-base tract, and an
    /// allele longer than a read can only ever be witnessed partially (spec §5.4.1).
    ///
    /// **An exhaustive match rather than `!= Complete`**, so that a third [`ReadWitness`]
    /// variant is a compile error here rather than silently joining the censored term. The
    /// censored term is the one that must never receive an observation it was not written
    /// for, and this is the single place that decides what reaches it.
    pub fn partial_observations(
        &self,
    ) -> impl Iterator<Item = (usize, &'a SequenceObservation)> + use<'a> {
        self.observations
            .iter()
            .enumerate()
            .filter(|(_, observation)| match observation.read_witness {
                ReadWitness::Partial { .. } => true,
                ReadWitness::Complete => false,
            })
    }
}

/// Where one sample's narrowed merge rows live while the SNP/indel row reads them.
///
/// **One per worker, reused across every sample of every locus** — which is the point, and
/// is why it is cleared and refilled rather than freshly made.
///
/// **It is deliberately not the row's scratch, and calling it that was a mistake this step
/// corrected.** The evidence view *borrows* this buffer, so it is still borrowed while the
/// row runs — and a row taking `&mut` the same object as the evidence borrows cannot be
/// called at all. Compiled, the mistake is `error[E0499]: cannot borrow as mutable more than
/// once`, at the first call site the next milestone writes. The two have different lifetimes
/// of use: what the row reads has to outlive the call, and what the row scribbles in does
/// not. So they are two types, and **the row's own scratch arrives with the step that first
/// needs one** — Milestone D, which adds a compatibility cache per `(partial, allele)` and a
/// gather buffer for a witness with holes in it. Inventing it empty here would be shape
/// without substance.
///
/// **Amortised, not allocation-free**, and the difference is worth stating rather than
/// implying. The buffer grows only when a sample is wider than every sample this worker has
/// met — measured over ten samples of 3, 7, 2, 11, 4, 25, 5, 9, 40 and 6 rows, five
/// reallocations, and none at all on a second pass. Capacity never comes back down, so a
/// worker ends up holding its widest sample's row count for its lifetime.
#[derive(Default, Debug)]
pub struct GenericEvidenceBuffer {
    supported: Vec<GenericObservation>,
}

impl GenericEvidenceBuffer {
    /// Narrow one sample's merge rows into the buffer under selection's mapping, and hand
    /// back the rows together with the quality of the ones the mapping dropped.
    ///
    /// The dropped quality comes back rather than being folded in, because it is only *part*
    /// of the pooled leftover the formula wants — selection pools its own, and the caller is
    /// the one holding both halves.
    ///
    /// # Panics
    ///
    /// As [`GenericObservation::fill_from_supported_alleles`], on a mapping that does not
    /// cover the merge's table.
    pub fn narrow_supported(
        &mut self,
        rows: &[SupportedAllele],
        candidate_of_merge_allele: &[Option<AlleleId>],
    ) -> (&[GenericObservation], f64) {
        let dropped = GenericObservation::fill_from_supported_alleles(
            rows,
            candidate_of_merge_allele,
            &mut self.supported,
        );
        (&self.supported, dropped)
    }
}

/// The buffers the STR row works in, held by the caller so the row allocates nothing.
///
/// **No borrow conflict here, unlike the generic path's**: the STR evidence borrows the locus
/// generator's own observations, not anything of this type's, so a caller can hold the
/// evidence and hand the row `&mut` this at the same time.
///
/// `ModelScratch` is the emission's own scratch — its placement buffer and its alignment
/// matrix (spec §8) — which is an associated type on the emission seam because each model
/// keeps its own shape. **The parameter is the model's scratch and not the model**, which is
/// worth the longer name: a stateless model that happened to derive `Default` would satisfy
/// `SsrRowScratch<ThatModel>` and compile, leaving the row with an empty model where its
/// placement buffer should be and an allocation per call.
///
/// **`emissions` is the cache that decides what a row costs.** One entry per
/// `(observation, candidate)`, filled once and read by every genotype that carries the
/// candidate, so a row costs `observations × candidates` rather than
/// `observations × genotypes` — 6 calls an observation against 21 at six candidates and a
/// diploid, a factor of 3.5. Spec §8
/// calls that the design and not an optimisation, which is why the cache is a field of the
/// scratch rather than a local inside whoever fills it.
///
/// **The observation axis is every observation, partials included.** That is the axis the two
/// filters' positions are on — both enumerate the whole slice and then filter — so a cache
/// sized over the complete observations alone is indexed past its end by the first complete
/// observation sitting above a partial. [`prepare_emissions`](Self::prepare_emissions) takes
/// the evidence rather than a count so that a caller cannot pick the wrong one, and
/// [`emission_at`](Self::emission_at) is the only spelling of the index.
///
/// **Amortised, not allocation-free**: measured over twelve loci of 24 to 840 slots, six
/// reallocations cold and none on a second pass, with capacity never coming back down.
#[derive(Default, Debug)]
pub struct SsrRowScratch<ModelScratch> {
    emissions: Vec<f64>,
    candidates: usize,
    model: ModelScratch,
}

impl<ModelScratch> SsrRowScratch<ModelScratch> {
    /// Size the emission cache for one locus and clear it to a stated value.
    ///
    /// **Cleared and not merely resized.** `Vec::resize` leaves the leading entries as they
    /// were, so a locus with as many `(observation, candidate)` pairs as the last one would
    /// silently reuse the last one's emissions — a wrong likelihood at every sample, with
    /// nothing failing.
    ///
    /// **The fill value is the caller's, and it must not be zero.** An unwritten slot has to
    /// hold something the row cannot mistake for a real score, and zero is exactly that
    /// mistake: a slip the candidate cannot reach legitimately scores zero (spec §4.2), so a
    /// cache pre-filled with zeros makes *never computed* and *computed as impossible* the
    /// same value.
    ///
    /// **The observation count comes from the evidence rather than from the caller**, because
    /// the caller has two plausible numbers to choose between and only one of them is right —
    /// see the type's own documentation.
    pub fn prepare_emissions(
        &mut self,
        evidence: &SsrSampleEvidence<'_>,
        candidates: usize,
        fill: f64,
    ) {
        self.candidates = candidates;
        self.emissions.clear();
        self.emissions
            .resize(evidence.observations.len() * candidates, fill);
    }

    /// One cached emission, by the position the evidence's filters yield and the candidate.
    ///
    /// **The one spelling of the index**, so the stride cannot be written two ways. A dense
    /// counter over the complete observations alone — the obvious thing to reach for inside a
    /// filtered loop — addresses the wrong row for every observation above a partial, and
    /// nothing about that fails.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on a position or candidate past what the cache was
    /// prepared for. Indexing past the end is what the wrong observation count produces, and
    /// silently reading a neighbouring candidate's emission is worse than stopping.
    pub fn emission_at(&self, observation: usize, candidate: usize) -> f64 {
        self.emissions[self.slot(observation, candidate)]
    }

    /// Write one cached emission. See [`emission_at`](Self::emission_at) for the index.
    pub fn set_emission(&mut self, observation: usize, candidate: usize, value: f64) {
        let slot = self.slot(observation, candidate);
        self.emissions[slot] = value;
    }

    /// How many slots the cache holds — one per `(observation, candidate)`.
    pub fn emission_count(&self) -> usize {
        self.emissions.len()
    }

    /// The whole cache, for a filler that walks it in order rather than by index.
    pub fn emissions_mut(&mut self) -> &mut [f64] {
        &mut self.emissions
    }

    /// The emission model's own scratch — its placement buffer and its alignment matrix.
    pub fn model_scratch_mut(&mut self) -> &mut ModelScratch {
        &mut self.model
    }

    fn slot(&self, observation: usize, candidate: usize) -> usize {
        assert!(
            candidate < self.candidates,
            "candidate {candidate} is past the {} this cache was prepared for",
            self.candidates
        );
        let slot = observation * self.candidates + candidate;
        assert!(
            slot < self.emissions.len(),
            "observation {observation} addresses slot {slot}, past the {} this cache was \
             prepared for — the cache was sized over the complete observations rather than \
             over all of them",
            self.emissions.len()
        );
        slot
    }
}

/// **The run's per-read-group calibration, and the one spelling of the lookup.**
///
/// A wrapper rather than a bare slice because it is what the *frequency-free* half of the row
/// is handed: [`generic::fill_generic_emissions`] takes this and not
/// [`ReadGroupParameters`], so that a fill whose results are cached across the loop's passes
/// cannot reach a quantity that moves between them (spec §6.1).
#[derive(Copy, Clone, Debug)]
pub struct ReadGroupCalibrations<'a>(&'a [ReadGroupCalibration]);

impl<'a> ReadGroupCalibrations<'a> {
    /// One calibration per read group, indexed by [`ReadGroupId`].
    #[must_use]
    pub fn over(calibration: &'a [ReadGroupCalibration]) -> Self {
        Self(calibration)
    }

    /// How many read groups the run supplied a calibration for.
    ///
    /// **No `is_empty` beside it**, which clippy would normally ask for: an empty one is a run
    /// whose read-group axis went missing, and `FrozenParameters` refuses it at construction, so
    /// the question has one answer everywhere it could be asked.
    #[must_use]
    pub fn count(&self) -> usize {
        self.0.len()
    }

    /// One read group's calibration.
    ///
    /// # Panics
    ///
    /// On a read group the run supplied no calibration for. A row scoring evidence from a group
    /// the parameters do not cover is a caller bug, and reading it as uncalibrated would be a
    /// genotype quietly moved rather than a run that stopped.
    #[must_use]
    pub fn of(&self, read_group: ReadGroupId) -> &'a ReadGroupCalibration {
        let at = read_group.get() as usize;
        self.0.get(at).unwrap_or_else(|| {
            panic!(
                "read group {at} has no calibration; the run supplied {}",
                self.0.len()
            )
        })
    }
}

/// Everything the SNP/indel row knows about a read group — its calibration and its share of
/// somebody else's DNA.
///
/// **They travel together because they are read off the same axis and have to agree about it.**
/// Both are dense over [`ReadGroupId`], and a run whose calibration covers ten read groups and
/// whose mixture covers four is a caller bug the row can only find lazily — at whichever locus
/// first holds a read from the fifth, or never. Pairing them puts that check at construction,
/// which is where [`ContaminationMixture::new`] already puts the checks between the mixture's
/// own two halves.
///
/// It also keeps the row's argument list inside what a reader can hold: the row takes what the
/// sample showed, what the locus is, what the run's chemistry is, its scratch, and its output.
#[derive(Copy, Clone, Debug)]
pub struct ReadGroupParameters<'a> {
    calibration: ReadGroupCalibrations<'a>,
    contamination: ContaminationMixture<'a>,
}

impl<'a> ReadGroupParameters<'a> {
    /// # Panics
    ///
    /// **In release as well as debug**, when the two disagree about how many read groups the
    /// run has — unless there is no mixture at all, which is the uncontaminated case and says
    /// nothing about read groups.
    #[must_use]
    pub fn new(
        calibration: &'a [ReadGroupCalibration],
        contamination: ContaminationMixture<'a>,
    ) -> Self {
        assert!(
            contamination.is_absent() || contamination.read_group_count() == calibration.len(),
            "the mixture holds a fraction for {} read groups and the run supplied {} \
             calibrations, so one of them belongs to a different run",
            contamination.read_group_count(),
            calibration.len()
        );
        Self {
            calibration: ReadGroupCalibrations::over(calibration),
            contamination,
        }
    }

    /// A run with no contamination estimate — a single sample, or a fit that emitted none.
    #[must_use]
    pub fn uncontaminated(calibration: &'a [ReadGroupCalibration]) -> Self {
        Self::new(calibration, ContaminationMixture::uncontaminated())
    }

    /// One read group's calibration.
    ///
    /// # Panics
    ///
    /// On a read group the run supplied no calibration for. A row scoring evidence from a group
    /// the parameters do not cover is a caller bug, and reading it as uncalibrated would be a
    /// genotype quietly moved rather than a run that stopped.
    #[must_use]
    pub fn calibration_of(&self, read_group: ReadGroupId) -> &'a ReadGroupCalibration {
        self.calibration.of(read_group)
    }

    /// Just the calibration half — what the row's frequency-free fill is handed.
    #[must_use]
    pub fn calibrations(&self) -> ReadGroupCalibrations<'a> {
        self.calibration
    }

    /// The mixture, for the two halves the row reads per observation.
    #[must_use]
    pub fn contamination(&self) -> ContaminationMixture<'a> {
        self.contamination
    }

    /// How many read groups the run has.
    #[must_use]
    pub fn read_group_count(&self) -> usize {
        self.calibration.count()
    }
}

/// The SNP/indel row's **emission cache** — everything about one sample's evidence at one locus
/// that reads no allele frequency, so that it can be computed once and read by every assembly of
/// the genotype-likelihood row — of which a contaminated locus makes one per pass.
///
/// Two things live here, and they are the whole of what spec §6.1 calls the emission on this
/// path:
///
/// - **one verdict per `(partial observation, candidate allele)`** — whether that allele could
///   have produced that read;
/// - **one charged error per observation**, complete and partial alike — `ε̄`, the geometric mean
///   of the reads' own error probabilities scaled by the read group's calibration.
///
/// # Why they are cached rather than recomputed
///
/// Two reasons, and the second is the one this type's shape is for.
///
/// **Within one row**, whether an allele could have produced a partial read is a property of the
/// two of them and not of the genotype being scored. Every genotype carrying that allele asks
/// the same question, so without a cache the answer is recomputed `genotypes` times instead of
/// once: at a six-allele diploid, 21 times over rather than 6. It is the same argument the STR
/// path's emission cache rests on, one axis narrower.
///
/// **Across the loop's passes**, with contamination on the row is no longer a constant: `q(o)`
/// moves with the frequencies (spec §3.6). What does not move is everything in here, so the
/// driver fills this once per locus and assembles a row from it per pass — which is what keeps
/// the expensive
/// half, `candidates × Σ_s (observations in sample s)`, independent of the pass count
/// (spec §6.1). **One of these per scratch row**, therefore, and not one reused across samples:
/// a shared one would be overwritten by the next sample before the next pass could read it.
///
/// # Why it is a separate type from [`GenericEvidenceBuffer`]
///
/// **The evidence view borrows that buffer**, so it is still borrowed while the row runs, and a
/// row taking `&mut` the same object the evidence borrows cannot be called at all. A2 recorded
/// that and deferred this type to the step that first had something to put in it, which is this
/// one.
///
/// # What is not here, and was expected to be
///
/// A2 and the plan both expected a **gather buffer** — somewhere to assemble an allele's bases
/// restricted to a witness with holes in it. There is none, because the comparison turned out
/// not to need one: an allele is the whole locus as a carrier has it, so a partial read is
/// checked against the allele's *prefix* or *suffix* and nothing is assembled
/// ([`generic::allele_is_compatible_with_partial`]). The buffer was a consequence of reading the rule as a
/// positional restriction; it is not one (owner, 2026-08-24).
#[derive(Default, Debug)]
pub struct GenericRowScratch {
    /// Partial-major, one row of `alleles` per partial observation.
    compatible: Vec<bool>,
    /// One entry per observation this sample showed: the complete observations first, in the
    /// evidence's own order, then the partials in theirs.
    ///
    /// **One buffer over both, rather than two**, because the two walks are the same walk with
    /// a different spread: an unexplained complete observation is charged `ε̄/m` and an
    /// unexplained partial `ε̄/1`, and both numbers are the same `ε̄` scaled.
    charged_error: Vec<f64>,
    supported: usize,
    alleles: usize,
}

impl GenericRowScratch {
    /// Size the cache for one sample at one locus and clear it.
    ///
    /// **The verdicts are filled `false` rather than left alone**, so one that was never written
    /// reads as *this allele cannot have produced this read* — the conservative direction, and
    /// one an incomplete fill shows up as a partial that stopped constraining anything rather
    /// than as one that silently constrained the wrong genotype.
    ///
    /// **The charged errors are filled with `NaN`, and that is the opposite choice for the
    /// opposite reason.** There is no conservative charged error: zero is a real value — a
    /// library whose reads are never wrong — so a slot the fill skipped would take every
    /// genotype that fails to explain the observation to `−∞` and look like evidence. It is the
    /// same argument [`SsrRowScratch::prepare_emissions`] makes for its own fill value.
    pub fn prepare_emissions(&mut self, supported: usize, partials: usize, alleles: usize) {
        self.alleles = alleles;
        self.supported = supported;
        self.compatible.clear();
        self.compatible.resize(partials * alleles, false);
        self.charged_error.clear();
        self.charged_error.resize(supported + partials, f64::NAN);
    }

    /// One complete observation's charged error `ε̄`, by its position in
    /// [`GenericSampleEvidence::supported`].
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on a position past what the cache was prepared for —
    /// which is what an assembly paired with another sample's fill produces, and reading a
    /// neighbouring observation's `ε̄` is worse than stopping.
    #[must_use]
    pub fn charged_error_of_supported(&self, observation: usize) -> f64 {
        assert!(
            observation < self.supported,
            "complete observation {observation} is past the {} this cache was prepared for",
            self.supported
        );
        self.charged_error[observation]
    }

    /// Write one complete observation's charged error. See
    /// [`charged_error_of_supported`](Self::charged_error_of_supported) for the index.
    pub fn set_supported_charged_error(&mut self, observation: usize, charged_error: f64) {
        assert!(
            observation < self.supported,
            "complete observation {observation} is past the {} this cache was prepared for",
            self.supported
        );
        self.charged_error[observation] = charged_error;
    }

    /// One partial read's charged error, **already divided by its spread of one**, by its
    /// position in [`GenericSampleEvidence::partials`].
    ///
    /// # Panics
    ///
    /// As [`charged_error_of_supported`](Self::charged_error_of_supported). **The two indices
    /// are separate accessors and not one axis with an offset**, because the partials sit above
    /// the complete observations in one buffer: a caller adding the offset by hand would address
    /// a real value from the wrong half whenever it got the base wrong.
    #[must_use]
    pub fn charged_error_of_partial(&self, partial: usize) -> f64 {
        self.charged_error[self.partial_slot(partial)]
    }

    /// Write one partial read's charged error. See
    /// [`charged_error_of_partial`](Self::charged_error_of_partial) for the index.
    pub fn set_partial_charged_error(&mut self, partial: usize, charged_error: f64) {
        let slot = self.partial_slot(partial);
        self.charged_error[slot] = charged_error;
    }

    /// How many complete observations the cache was prepared for.
    #[must_use]
    pub fn supported_count(&self) -> usize {
        self.supported
    }

    /// How many partial reads the cache was prepared for.
    #[must_use]
    pub fn partial_count(&self) -> usize {
        self.charged_error.len() - self.supported
    }

    /// How many candidate alleles the cache was prepared for.
    #[must_use]
    pub fn allele_count(&self) -> usize {
        self.alleles
    }

    fn partial_slot(&self, partial: usize) -> usize {
        let slot = self.supported + partial;
        assert!(
            slot < self.charged_error.len(),
            "partial {partial} is past the {} this cache was prepared for",
            self.charged_error.len() - self.supported
        );
        slot
    }

    /// Whether this allele could have produced this partial read.
    #[must_use]
    pub fn is_compatible(&self, partial: usize, allele: AlleleId) -> bool {
        self.compatible[self.slot(partial, allele)]
    }

    /// Record one verdict.
    pub fn set_compatible(&mut self, partial: usize, allele: AlleleId, compatible: bool) {
        let slot = self.slot(partial, allele);
        self.compatible[slot] = compatible;
    }

    /// How many verdicts the cache is sized for.
    #[must_use]
    pub fn verdict_count(&self) -> usize {
        self.compatible.len()
    }

    /// **The one spelling of the two-dimensional index**, asserting in release for
    /// [`SsrRowScratch`]'s reason: a cache sized for one locus and read at another's stride
    /// returns a real verdict from the wrong row rather than running off the end.
    fn slot(&self, partial: usize, allele: AlleleId) -> usize {
        let allele = usize::from(allele.get());
        assert!(
            allele < self.alleles,
            "allele {allele} is past the {} this cache was prepared for",
            self.alleles
        );
        let slot = partial * self.alleles + allele;
        assert!(
            slot < self.compatible.len(),
            "partial {partial} is past the {} this cache was prepared for",
            self.compatible.len() / self.alleles.max(1)
        );
        slot
    }
}

/// Fixtures both this module's tests and `generic.rs`'s need, in one place.
///
/// **They were copied byte-identically into both test modules and would have drifted** — the
/// batching is exactly the kind of fixture that gets adjusted on one side of a change (C2's
/// review).
#[cfg(test)]
pub(crate) mod test_batching {
    use super::{BatchId, BatchOfEachReadGroup, ContaminationMixture, ContaminationView};

    /// Enough entries that any fixture's read groups fit; a mixture takes the prefix it needs.
    /// It cannot truncate silently — slicing past 16 panics, and the widest fixture uses 4.
    pub(crate) static ALL_IN_ONE_BATCH: [BatchId; 16] = [BatchId::ALL_TOGETHER; 16];

    /// A mixture under **the default batching** — one batch holding every read group, so every
    /// observation reads the same contaminant frequencies. That is what a run which declares no
    /// batching gets, and it is the shape all but the batch-specific fixtures want.
    pub(crate) fn one_batch<'a>(
        fractions: &'a [ContaminationView],
        frequencies: &'a [f64],
    ) -> ContaminationMixture<'a> {
        ContaminationMixture::new(
            fractions,
            BatchOfEachReadGroup(&ALL_IN_ONE_BATCH[..fractions.len()]),
            frequencies,
            frequencies.len(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::test_batching::{ALL_IN_ONE_BATCH, one_batch};
    use super::*;
    use crate::ng::locus_generation::WitnessedLocusPositions;
    use crate::ng::run::cohort_merge::build::AlleleSupport;
    use crate::ng::types::Motif;

    /// A merge row with every field distinguishable, so a field crossed with another one
    /// changes an answer. `num_fwd`, the mapping moments and `placed_left` are deliberately
    /// not zero: a fixture whose uninteresting fields are all zero cannot show that a
    /// mapping of the interesting ones picked them up by accident.
    fn supported_row(
        allele: usize,
        read_group: u32,
        num_reads: u32,
        q_sum: f64,
    ) -> SupportedAllele {
        SupportedAllele {
            allele,
            read_group: ReadGroupId(read_group),
            support: AlleleSupport {
                num_reads,
                num_fwd: 7,
                q_sum,
                mapq_sum: 611,
                mapq_sum_sq: 37_442,
                placed_left: 3,
            },
        }
    }

    /// Selection kept every merge allele, in order — the mapping a locus gets before
    /// anything is pruned.
    fn kept_in_order(alleles: usize) -> Vec<Option<AlleleId>> {
        (0..alleles)
            .map(|index| Some(AlleleId(index as u16)))
            .collect()
    }

    /// A partial whose witnessed stretch and bases are **deliberately different lengths** —
    /// five locus positions against seven bases, the shape a read carrying a two-base
    /// insertion inside the stretch comes back with.
    ///
    /// A fixture where the two happen to be equal is the one a consumer indexing the bases
    /// by a locus offset would pass, and the merge's own documentation says outright that
    /// equality does not license indexing one with the other.
    fn partial_row(num_reads: u32, q_sum: f64) -> PartialObservation {
        PartialObservation {
            witnessed_in_locus: WitnessedLocusPositions::one_run_from_offset_and_length(2, 5)
                .expect("a five-position run from offset two is a witness"),
            read_group: ReadGroupId(1),
            bases: b"ACGTTAA"[..].into(),
            num_reads,
            q_sum,
        }
    }

    /// The tract's repeat detail, shared by the repeat-path fixtures.
    fn ssr_detail() -> SsrDetail {
        SsrDetail {
            motif: Motif::new(b"AC").expect("AC is a motif"),
            left_flank: b"GGTT"[..].into(),
            right_flank: b"TTGG"[..].into(),
        }
    }

    /// A read that ran out after `positions` positions of the tract.
    fn ran_out_after(positions: u16) -> ReadWitness {
        ReadWitness::Partial {
            positions: WitnessedLocusPositions::one_run_from_offset_and_length(0, positions)
                .expect("a run from offset zero over some positions is a witness"),
        }
    }

    fn ssr_observation(bases: &[u8], witness: ReadWitness, num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.into(),
            read_witness: witness,
            read_group: ReadGroupId(2),
            num_obs,
            num_fwd: 5,
            q_sum: -13.5,
            mapq_sum: 240,
            mapq_sum_sq: 14_400,
            placed_left: 1,
            chain_ids: Vec::new(),
        }
    }

    #[test]
    fn the_view_of_a_merge_row_keeps_the_four_numbers_the_formula_reads() {
        let observation =
            GenericObservation::of_supported_allele(&supported_row(3, 9, 17, -41.25), AlleleId(3));

        assert_eq!(observation.allele, AlleleId(3));
        assert_eq!(observation.read_group, ReadGroupId(9));
        assert_eq!(observation.num_reads, 17);
        assert_eq!(observation.q_sum, -41.25);
    }

    /// **And the two the formula does not read**, which is the whole reason they are here:
    /// the site quality's artifact correction compares the reference allele's strand and
    /// read-position split against the alternative's, and the merge row they come from is
    /// released with the locus (`doc/devel/ng/spec/calling_quality.md` §3.3).
    ///
    /// The fixture's `num_fwd` is 7 and its `placed_left` is 3 — deliberately different from
    /// each other and from every other field, so a copy that crossed the two, or took
    /// `num_reads` for either, changes this assertion.
    #[test]
    fn the_view_also_keeps_the_two_counts_the_site_quality_reads() {
        let observation =
            GenericObservation::of_supported_allele(&supported_row(3, 9, 17, -41.25), AlleleId(3));

        assert_eq!(observation.forward_reads, 7);
        assert_eq!(observation.placed_left_reads, 3);
    }

    /// The type's doc promises plain scalars and cheap copying. A field owning heap would
    /// make both false, and would also break the no-allocation contract the row function
    /// works under (spec §8), so the width is pinned rather than described.
    #[test]
    fn the_observation_stays_a_handful_of_scalars_wide() {
        // **32, and it was 24 before the two site-quality counters joined it.** Not 34: the
        // `f64` wants eight-byte alignment, so the two-byte allele id and the four
        // four-byte counts are padded out to twenty-four before it. The two counters cost
        // eight bytes on a struct the merge already pays 48 for.
        assert_eq!(std::mem::size_of::<GenericObservation>(), 32);
        assert_eq!(std::mem::align_of::<GenericObservation>(), 8);
        assert!(!std::mem::needs_drop::<GenericObservation>());
    }

    /// The merge's allele index is **not** the candidate id, and the mapping between them is
    /// the argument. Here selection has dropped merge allele 1, so merge allele 2 becomes
    /// candidate 1 — and a fill that had assumed the identity would have scored those reads
    /// against the wrong sequence with nothing saying so.
    #[test]
    fn a_dropped_allele_renumbers_the_ones_above_it_and_its_quality_comes_back() {
        let rows = [
            supported_row(0, 1, 11, -3.5),
            supported_row(1, 1, 4, -19.0),
            supported_row(2, 1, 6, -8.25),
        ];
        let mapping = [Some(AlleleId(0)), None, Some(AlleleId(1))];
        let mut out = Vec::new();

        let dropped = GenericObservation::fill_from_supported_alleles(&rows, &mapping, &mut out);

        assert_eq!(dropped, -19.0);
        assert_eq!(
            out.iter().map(|o| o.allele).collect::<Vec<_>>(),
            vec![AlleleId(0), AlleleId(1)]
        );
        assert_eq!(
            out.iter().map(|o| o.num_reads).collect::<Vec<_>>(),
            vec![11, 6]
        );
    }

    /// Nothing is dropped when selection keeps everything, and the pooled quality is then
    /// zero rather than the sum of what was kept — a fill that added on the wrong branch
    /// would return −30.75 here.
    #[test]
    fn keeping_every_allele_drops_no_quality() {
        let rows = [
            supported_row(0, 1, 11, -3.5),
            supported_row(1, 2, 4, -27.25),
        ];
        let mut out = Vec::new();

        let dropped =
            GenericObservation::fill_from_supported_alleles(&rows, &kept_in_order(2), &mut out);

        assert_eq!(dropped, 0.0);
        assert_eq!(out.len(), 2);
    }

    /// `out` is caller scratch: it is cleared, so a buffer carrying the previous sample's
    /// rows comes back holding this sample's and nothing else.
    #[test]
    fn the_scratch_buffer_holds_no_trace_of_the_previous_sample() {
        let previous = [supported_row(0, 4, 99, -1.0), supported_row(1, 4, 98, -2.0)];
        let current = [supported_row(0, 1, 11, -3.5)];
        let mut out = Vec::new();

        let _ =
            GenericObservation::fill_from_supported_alleles(&previous, &kept_in_order(2), &mut out);
        let _ =
            GenericObservation::fill_from_supported_alleles(&current, &kept_in_order(2), &mut out);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].num_reads, 11);
    }

    /// A mapping too short for the merge's table was built at another locus. Scoring those
    /// reads against whatever allele sat at that index is the failure this refuses, and it
    /// refuses it in every build profile — `[profile.release]` turns debug assertions off.
    #[test]
    #[should_panic(expected = "was built against a different locus's table")]
    fn a_mapping_that_does_not_cover_the_merges_table_is_a_caller_bug() {
        let rows = [supported_row(2, 1, 6, -8.25)];
        let mut out = Vec::new();

        let _ = GenericObservation::fill_from_supported_alleles(&rows, &kept_in_order(2), &mut out);
    }

    /// The last allele the mapping covers is not itself out of range — a bound written with
    /// `<` where `<=` was meant, or the other way about, would reject it.
    #[test]
    fn the_last_allele_the_mapping_covers_is_accepted() {
        let rows = [supported_row(1, 1, 6, -8.25)];
        let mut out = Vec::new();

        let dropped =
            GenericObservation::fill_from_supported_alleles(&rows, &kept_in_order(2), &mut out);

        assert_eq!(dropped, 0.0);
        assert_eq!(out[0].allele, AlleleId(1));
    }

    /// The constructor carries all three parts through, **and the slices are compared whole
    /// rather than counted**: a length hides a row dropped and another duplicated, or two
    /// reordered, and the order is a determinism requirement rather than tidiness (spec §8).
    ///
    /// The pooled leftover is the part that can go missing without any genotype moving — it
    /// is the same number in every row of the sample — so nothing but a direct comparison
    /// would notice, and it is the term the data likelihood needs to be comparable between
    /// loci.
    #[test]
    fn the_constructor_keeps_the_pooled_leftover_and_both_slices_whole() {
        let supported = [
            GenericObservation::of_supported_allele(&supported_row(0, 1, 11, -3.5), AlleleId(0)),
            GenericObservation::of_supported_allele(&supported_row(1, 2, 4, -19.0), AlleleId(1)),
        ];
        let partials = [partial_row(6, -22.5), partial_row(2, -7.25)];

        let evidence = GenericSampleEvidence::new(&supported, -7.5, &partials);

        assert_eq!(evidence.unmatched_q_sum, -7.5);
        assert_eq!(evidence.supported, &supported[..]);
        assert_eq!(evidence.partials, &partials[..]);
    }

    #[test]
    fn a_sample_that_showed_nothing_carries_no_reads_and_no_pooled_leftover() {
        let evidence = GenericSampleEvidence::empty();

        assert_eq!(evidence.unmatched_q_sum, 0.0);
        assert!(evidence.supported.is_empty());
        assert!(evidence.partials.is_empty());
    }

    /// The merge's order is the view's order, and this constructor takes any slice. A call
    /// site that concatenated two samples' rows, or applied a selection mapping that
    /// reordered the alleles, would break the order with nothing saying so — and
    /// floating-point addition is not associative, so the run's genotypes move with it.
    /// **Debug-profile only, following the convention `genotype_prior/` already uses.** The
    /// check this reaches is a `debug_assert!` by a recorded decision, not by oversight — see
    /// [`GenericSampleEvidence::new`]'s doc: the scan is linear in the rows and the order is
    /// the merge's property to hold rather than this constructor's to re-establish at every
    /// call. Without the `cfg` the test fails under `--release` with *"test did not panic as
    /// expected"*, which is what kept `cargo test --release --lib ng::calling` red and so kept
    /// **every** release-held assertion in this folder outside CI's reach.
    ///
    /// **This attribute is not a way to quiet a failing test.** It says the check under it is
    /// deliberately debug-only. A test on a check that is meant to hold in release must not
    /// carry it, or the release gate stops meaning anything.
    #[test]
    #[should_panic(expected = "not ascending")]
    #[cfg(debug_assertions)]
    fn rows_out_of_pair_order_are_refused_in_debug() {
        let out_of_order = [
            GenericObservation::of_supported_allele(&supported_row(1, 0, 3, -1.0), AlleleId(1)),
            GenericObservation::of_supported_allele(&supported_row(0, 0, 3, -1.0), AlleleId(0)),
        ];

        let _ = GenericSampleEvidence::new(&out_of_order, 0.0, &[]);
    }

    /// The read group is the second half of the key, so rows ascending on the allele alone
    /// are not ascending — a check that compared only the allele would let them through.
    /// **One `(allele, read group)` pair appearing twice is what concatenating two samples'
    /// rows produces**, and the sum would then charge that evidence twice. The order check has
    /// to be strict rather than merely non-decreasing, and both fixtures below are strictly
    /// descending, so neither puts an equal pair in front of it.
    /// **Debug-profile only** — the check it reaches is a `debug_assert!` by the recorded
    /// decision in [`GenericSampleEvidence::new`]'s doc, not by oversight. See
    /// [`rows_out_of_pair_order_are_refused_in_debug`] for why the attribute is here and what
    /// it must never be used for.
    #[test]
    #[should_panic(expected = "not ascending")]
    #[cfg(debug_assertions)]
    fn one_pair_appearing_twice_is_refused_in_debug() {
        let duplicated = [
            GenericObservation::of_supported_allele(&supported_row(1, 0, 3, -1.0), AlleleId(1)),
            GenericObservation::of_supported_allele(&supported_row(1, 0, 5, -2.0), AlleleId(1)),
        ];

        let _ = GenericSampleEvidence::new(&duplicated, 0.0, &[]);
    }

    /// **Debug-profile only** — the check it reaches is a `debug_assert!` by the recorded
    /// decision in [`GenericSampleEvidence::new`]'s doc, not by oversight. See
    /// [`rows_out_of_pair_order_are_refused_in_debug`] for why the attribute is here and what
    /// it must never be used for.
    #[test]
    #[should_panic(expected = "not ascending")]
    #[cfg(debug_assertions)]
    fn rows_out_of_read_group_order_within_one_allele_are_refused_in_debug() {
        let out_of_order = [
            GenericObservation::of_supported_allele(&supported_row(2, 9, 3, -1.0), AlleleId(2)),
            GenericObservation::of_supported_allele(&supported_row(2, 4, 3, -1.0), AlleleId(2)),
        ];

        let _ = GenericSampleEvidence::new(&out_of_order, 0.0, &[]);
    }

    #[test]
    fn the_repeat_view_tells_complete_reads_from_reads_that_ran_out() {
        let observations = [
            ssr_observation(b"ACACACAC", ReadWitness::Complete, 12),
            ssr_observation(b"ACAC", ran_out_after(4), 5),
            ssr_observation(b"ACACAC", ReadWitness::Complete, 3),
        ];
        let detail = ssr_detail();
        let evidence = SsrSampleEvidence::new(&observations, &detail);

        let complete: Vec<(usize, u32)> = evidence
            .complete_observations()
            .map(|(at, o)| (at, o.num_obs))
            .collect();
        let partial: Vec<(usize, u32)> = evidence
            .partial_observations()
            .map(|(at, o)| (at, o.num_obs))
            .collect();

        assert_eq!(complete, vec![(0, 12), (2, 3)]);
        assert_eq!(partial, vec![(1, 5)]);
        // The tract's own detail is the whole of what the emission is told about which
        // repeat it is scoring against, and nothing else here reads it back.
        assert_eq!(evidence.detail, &detail);
        assert_eq!(evidence.observations, &observations[..]);
    }

    /// The same bases seen by a read that spanned the tract and by one that ran out are
    /// **two** observations (spec §1.3), and only the witness separates them. A filter that
    /// read the bases instead — the plausible "is this shorter than the tract?" shortcut —
    /// puts both on the same side, and every other fixture here would still pass, because
    /// in each of them the partial's bases are also the shortest.
    #[test]
    fn two_observations_of_the_same_bases_split_on_the_witness_alone() {
        let observations = [
            ssr_observation(b"ACACAC", ReadWitness::Complete, 9),
            ssr_observation(b"ACACAC", ran_out_after(6), 4),
        ];
        let detail = ssr_detail();
        let evidence = SsrSampleEvidence::new(&observations, &detail);

        assert_eq!(
            evidence
                .complete_observations()
                .map(|(_, o)| o.num_obs)
                .collect::<Vec<_>>(),
            vec![9]
        );
        assert_eq!(
            evidence
                .partial_observations()
                .map(|(_, o)| o.num_obs)
                .collect::<Vec<_>>(),
            vec![4]
        );
    }

    /// Every observation is on exactly one side, and the positions the two yield are the
    /// positions in the slice. That is what lets a complete term and a censored term be
    /// summed without double-counting a read or losing one, and what lets one cache serve
    /// both — and it is a property of the *pair*: each filter's own test passes just as well
    /// if the two overlap.
    #[test]
    fn the_two_repeat_filters_partition_every_observation_and_agree_on_its_position() {
        let observations = [
            ssr_observation(b"ACACACAC", ReadWitness::Complete, 12),
            ssr_observation(b"ACAC", ran_out_after(4), 5),
            ssr_observation(b"ACACAC", ReadWitness::Complete, 3),
        ];
        let detail = ssr_detail();
        let evidence = SsrSampleEvidence::new(&observations, &detail);

        let mut seen: Vec<usize> = evidence
            .complete_observations()
            .chain(evidence.partial_observations())
            .map(|(at, _)| at)
            .collect();
        seen.sort_unstable();

        assert_eq!(seen, vec![0, 1, 2]);
        for (at, observation) in evidence
            .complete_observations()
            .chain(evidence.partial_observations())
        {
            assert_eq!(observation.num_obs, observations[at].num_obs);
        }
    }

    // ---- A2: the parameter views, the floors, the scratch ----

    /// The floor is inherited from production and the doc comment says so. Two spellings of
    /// one model assumption are two things that can drift, and nothing else would notice: a
    /// floor moved to `1e-9` changes an impossible observation's charge by 7 nats, which is a
    /// plausible-looking number rather than a failure.
    ///
    /// **Production's ceiling is deliberately not adopted and so is not pinned here** — see
    /// `what_the_row_charges_a_poor_read_is_not_capped`.
    #[test]
    fn the_floor_still_equals_productions() {
        assert_eq!(
            MIN_BASE_ERROR,
            crate::var_calling::contamination_estimation::MIN_BASE_ERROR
        );
    }

    /// The two clamps are named here by re-export rather than copied.
    ///
    /// **The equality half of this cannot fail while the re-export stands** — it compares an
    /// alias against its own source — and it is here for the day somebody replaces the
    /// re-export with a local constant, which is the drift the tree already shows: production's
    /// `ssr/cohort/read_model/hipstr.rs` holds a third private copy of these same two numbers
    /// with nothing connecting it to either. The `const` block below is the half that does
    /// independent work, and the distribution's own tests are what catch a moved value.
    #[test]
    fn the_geometric_clamps_are_the_stutter_distributions_own() {
        assert_eq!(GEOM_MIN, crate::ng::alignment::stutter::GEOM_MIN);
        assert_eq!(GEOM_MAX, crate::ng::alignment::stutter::GEOM_MAX);
        // A geometric success probability held strictly inside (0, 1): both clamps must be
        // in range and the low one below the high one, or the clamp inverts silently.
        const { assert!(GEOM_MIN > 0.0 && GEOM_MIN < GEOM_MAX && GEOM_MAX < 1.0) };
    }

    /// A rate the pre-pass fitted from this read group's own sites, with the read count it
    /// stood on. `observations` is deliberately not zero: it travels on the estimate and a
    /// consumer reading it off the wrong field would be invisible against a zero.
    fn fitted_rate(rate: f64) -> Estimate<ErrorRate> {
        Estimate {
            value: ErrorRate::try_new(rate).expect("a fixture rate is a probability"),
            provenance: Provenance::FittedHere,
            observations: 84_113,
        }
    }

    /// The same rate, borrowed from the sample's other read groups because this one had too
    /// little data — the case the calibration must not launder into a measurement.
    fn borrowed_rate(rate: f64) -> Estimate<ErrorRate> {
        Estimate {
            value: ErrorRate::try_new(rate).expect("a fixture rate is a probability"),
            provenance: Provenance::Borrowed,
            observations: 311,
        }
    }

    /// The log error of a read whose mapping quality is Phred 60, which is where BWA tops out.
    /// Every base quality this file's fixtures use is worse than that, so the base quality is
    /// what the mint returns — which is the ordinary case and the one worth exercising.
    const MAPPING_QUALITY_60_LOG_ERROR: f64 = -6.0 * std::f64::consts::LN_10;

    /// One read group's minted per-read errors, at the given base qualities, spread over two
    /// orders of magnitude on purpose — Phred 40, 30 and 13 is the spread the calibration
    /// exists to preserve. Returned as `(q_sum, num_reads)`, the shape the merge hands over.
    ///
    /// **It mints through the walk's own function rather than rolling the arithmetic**, which
    /// is spec §3.2's second requirement: the quantity the accumulator averages and the
    /// quantity the model charges must be computed by the same function, or the scale
    /// calibrates against a different definition of *how wrong is this read* than the one it
    /// is applied to. A local `10^(-q/10)` differs from the crate's table by one or two units
    /// in the last place — harmless in size, and a second definition of the quantity all the
    /// same.
    fn minted_reads_at_phred(base_qualities: &[u8]) -> (f64, u32) {
        let q_sum: f64 = base_qualities
            .iter()
            .map(|&quality| {
                crate::ng::locus_generation::pileup::minted_ln_read_error(
                    quality,
                    MAPPING_QUALITY_60_LOG_ERROR,
                )
            })
            .sum();
        (q_sum, base_qualities.len() as u32)
    }

    /// **The property the scale exists for** (spec §12 test 10): after scaling, the average
    /// error the model charges that read group's reads is the rate the fit measured.
    ///
    /// **To the accumulator's quantum, and not exactly** — and **both** bounds asserted below
    /// are needed. In real arithmetic the identity is exact: scaling every read by one
    /// multiplier multiplies their geometric mean by it, and the geometric mean is what both
    /// the denominator and the charge are. What moves the last few digits is that the
    /// denominator arrives rounded to units of `MintedReadErrors::LOG_ERROR_QUANTUM`, which
    /// shifts the mean log error by at most half a unit and the charged rate by that relative
    /// amount. Measured here: 4.7 × 10⁻¹⁰ against a rate of 4 × 10⁻³.
    ///
    /// **The derived bound alone would not have caught the accumulator getting coarser**,
    /// because it scales with the quantum it is derived from — a review mutation making the
    /// accumulator four times coarser left every test in this file and in `calibration.rs`
    /// green while pushing the real gap past the figure this module's own documentation
    /// quotes. So the absolute bound is asserted beside it, and it is the number in that
    /// documentation.
    #[test]
    fn the_scale_makes_the_charged_average_the_fitted_rate() {
        let (q_sum, num_reads) = minted_reads_at_phred(&[40, 30, 13]);
        let fitted = fitted_rate(0.004);
        let rate = fitted.value.get();
        let calibration = ReadGroupCalibration::from_fitted_rate(
            &fitted,
            MintedReadErrors::of_observation(q_sum, num_reads),
        )
        .expect("three reads are a denominator");

        // **Through `charged_error`, which is what the row charges** — not through the capped
        // reading beside it. Spec §3.2's calibrated property is *the average charged error
        // equals the measured rate*, and the cap is what breaks it, so pinning the property on
        // the capped function tests a claim only the other one makes. This fixture's reads sit
        // well inside the cap, so both return the same number here and the test passed either
        // way; that is exactly why it had to be moved rather than left to fail one day. (C1,
        // 2026-08-24 — the row stopped using the capped reading and this test did not follow.)
        let charged = calibration.charged_error(q_sum, num_reads);
        let relative_gap = (charged - rate).abs() / rate;

        assert!(
            relative_gap <= MintedReadErrors::LOG_ERROR_QUANTUM,
            "the calibrated average is {charged} and the fitted rate {}, {relative_gap} apart \
             relatively, against the accumulator's own quantum",
            rate
        );
        assert!(
            relative_gap <= 5e-7,
            "the module documents this as within five parts in ten million; it is \
             {relative_gap}"
        );
        assert_eq!(calibration.provenance, Provenance::FittedHere);
    }

    /// The identity is about **a read group's** average, not about one observation's charge,
    /// and those come apart the moment a group has more than one observation. Three
    /// observations of one read group here, at three different qualities, folded into the
    /// denominator the way the accumulator folds them — and it is the read-weighted geometric
    /// mean of the three charges that has to come back as the fitted rate, not any one of
    /// them.
    ///
    /// Nothing else in this file exercises that: the test above builds the calibration from a
    /// single observation carrying the whole group, where "the group's mean is the rate" and
    /// "this observation's charge is the rate" are the same sentence.
    #[test]
    fn the_average_is_the_groups_and_not_one_observations() {
        let observations = [
            minted_reads_at_phred(&[40, 40, 38]),
            minted_reads_at_phred(&[30]),
            minted_reads_at_phred(&[13, 15]),
        ];
        let mut group = MintedReadErrors::default();
        for &(q_sum, num_reads) in &observations {
            group.add(MintedReadErrors::of_observation(q_sum, num_reads));
        }
        let fitted = fitted_rate(0.004);
        let rate = fitted.value.get();
        let calibration = ReadGroupCalibration::from_fitted_rate(&fitted, group)
            .expect("six reads are a denominator");

        // The read-weighted geometric mean of the charges: Σ n·ln(charge) ÷ Σ n, exponentiated.
        let mut weighted_log_sum = 0.0;
        let mut reads = 0u32;
        for &(q_sum, num_reads) in &observations {
            let charge = calibration.charged_error(q_sum, num_reads);
            weighted_log_sum += f64::from(num_reads) * charge.ln();
            reads += num_reads;
        }
        let group_average = (weighted_log_sum / f64::from(reads)).exp();

        let relative_gap = (group_average - rate).abs() / rate;
        assert!(
            relative_gap <= 5e-7,
            "the group's charged average is {group_average}, the fitted rate {}",
            rate
        );
        // …and no single observation's charge is the fitted rate, which is the point.
        let single = calibration.charged_error(observations[1].0, observations[1].1);
        assert!((single - rate).abs() > 1e-4);
    }

    /// The residual above is the accumulator's rounding and nothing else, and **a log sum the
    /// fixed point can hold exactly closes it entirely**: −7 nats is a whole number of
    /// quanta, so it survives the round trip through 2⁻²⁰ units unchanged and the charged
    /// average is the fitted rate to the last bits of the double.
    ///
    /// That is what says the residual is quantisation rather than a mistake in the algebra.
    #[test]
    fn the_identity_is_exact_where_the_accumulators_rounding_does_not_bite() {
        // e^−7 is about 9 in ten thousand — a read a shade better than Phred 30.
        let q_sum = -7.0;
        let fitted = fitted_rate(0.004);
        let rate = fitted.value.get();
        let calibration = ReadGroupCalibration::from_fitted_rate(
            &fitted,
            MintedReadErrors::of_observation(q_sum, 1),
        )
        .expect("one read is a denominator");

        let charged = calibration.charged_error(q_sum, 1);

        assert!(
            (charged - rate).abs() <= f64::EPSILON * rate,
            "the calibrated average is {charged}, the fitted rate {}",
            rate
        );
    }

    /// The scale preserves the *shape* — which is the half the fitted rate cannot supply.
    /// Two observations of the same read group whose reads differ by 27 Phred must still
    /// differ after calibration, and by the same ratio they differed by before: one
    /// multiplier cannot flatten them.
    #[test]
    fn calibration_moves_the_size_and_leaves_the_ratio_between_reads_alone() {
        let (group_q_sum, group_reads) = minted_reads_at_phred(&[40, 30, 13]);
        let fitted = fitted_rate(0.004);
        let calibration = ReadGroupCalibration::from_fitted_rate(
            &fitted,
            MintedReadErrors::of_observation(group_q_sum, group_reads),
        )
        .expect("three reads are a denominator");

        let (good_q_sum, good_reads) = minted_reads_at_phred(&[40]);
        let (poor_q_sum, poor_reads) = minted_reads_at_phred(&[13]);
        let good = calibration.charged_error(good_q_sum, good_reads);
        let poor = calibration.charged_error(poor_q_sum, poor_reads);

        // Phred 13 against Phred 40 is a factor of 10^2.7 either side of the scale.
        assert!(
            (poor / good - 10f64.powf(2.7)).abs() < 1e-6,
            "ratio {}",
            poor / good
        );
        assert!(good < poor);
        assert!(calibration.scale != 1.0);
    }

    /// A read group with no read behind it has no denominator, and a scale needs one. The
    /// alternative — treating the absence as an average — divides the fitted rate by nothing.
    #[test]
    fn a_read_group_with_no_reads_yields_no_scale() {
        let fitted = fitted_rate(0.004);

        assert!(
            ReadGroupCalibration::from_fitted_rate(&fitted, MintedReadErrors::default()).is_none()
        );
    }

    /// **The other way a denominator can be useless: it underflowed to zero.** No real base
    /// quality reaches a mean log error of −1,000 nats, but the accumulator's constructor is
    /// public and a fixture can hand it one — and without the guard the scale comes back
    /// infinite, which charges every read of the library the ceiling.
    ///
    /// The `no reads` test above cannot reach this branch, because it stops at the `?` one
    /// line earlier.
    #[test]
    fn a_denominator_that_underflowed_to_zero_yields_no_scale() {
        let fitted = fitted_rate(0.004);
        let underflowed = MintedReadErrors::of_observation(-1000.0, 1);

        assert_eq!(underflowed.mean_error_probability(), Some(0.0));
        assert!(ReadGroupCalibration::from_fitted_rate(&fitted, underflowed).is_none());
    }

    /// **A fitted rate of zero is refused, and it is the numerator's turn to be guarded.** A
    /// scale of zero charges every read of that library `MIN_BASE_ERROR` — maximal confidence
    /// about every base, from a fit that says it found no errors at all. *Absent* and *zero*
    /// are different answers here for the same reason they are on the contamination fraction.
    #[test]
    fn a_fitted_rate_of_zero_yields_no_scale_rather_than_a_scale_of_zero() {
        let (q_sum, num_reads) = minted_reads_at_phred(&[30, 30]);
        let zero = fitted_rate(0.0);

        assert!(
            ReadGroupCalibration::from_fitted_rate(
                &zero,
                MintedReadErrors::of_observation(q_sum, num_reads)
            )
            .is_none()
        );
    }

    /// **The calibration's warrant is its rate's warrant, and stamping `FittedHere` on a
    /// borrowed rate would launder it.** A read group with too little data of its own takes
    /// the mean of the sample's other groups; chemistry differs between libraries, which is
    /// the whole reason for the read-group grain, so that is a compromise and the output has
    /// to be able to say so. The scale adds no warrant of its own — it is a ratio, exactly as
    /// well founded as its numerator.
    #[test]
    fn a_borrowed_rate_makes_a_borrowed_calibration() {
        let (q_sum, num_reads) = minted_reads_at_phred(&[40, 30, 13]);
        let borrowed = borrowed_rate(0.004);

        let calibration = ReadGroupCalibration::from_fitted_rate(
            &borrowed,
            MintedReadErrors::of_observation(q_sum, num_reads),
        )
        .expect("three reads are a denominator");

        assert_eq!(calibration.provenance, Provenance::Borrowed);
        // …and the scale itself is the same number either way: only the warrant differs.
        let fitted = fitted_rate(0.004);
        let from_fitted = ReadGroupCalibration::from_fitted_rate(
            &fitted,
            MintedReadErrors::of_observation(q_sum, num_reads),
        )
        .expect("three reads are a denominator");
        assert_eq!(calibration.scale, from_fitted.scale);
        assert_ne!(calibration.provenance, from_fitted.provenance);
    }

    /// **Held in every profile.** An observation with no reads behind it is a row the merge
    /// does not build, and without the check a release build returns `NaN` at a zero `q_sum`
    /// — which `f64::clamp` passes straight through the floors this module promises — or
    /// `MIN_BASE_ERROR` at a negative one, which nothing downstream could tell from a real
    /// charge.
    #[test]
    #[should_panic(expected = "has no average error")]
    fn an_observation_with_no_reads_is_a_caller_bug() {
        let _ = ReadGroupCalibration::defaulted().charged_error(-7.0, 0);
    }

    /// Where nothing was fitted the qualities are used as reported — and the run's output has
    /// to be able to say so, which is what the provenance is for. A defaulted calibration
    /// that claimed `FittedHere` would make a calibrated run and a trusting one look alike.
    #[test]
    fn a_defaulted_calibration_charges_the_reads_exactly_what_they_were_minted_with() {
        let (q_sum, num_reads) = minted_reads_at_phred(&[30, 30]);
        let defaulted = ReadGroupCalibration::defaulted();

        assert_eq!(defaulted.scale, 1.0);
        assert_eq!(defaulted.provenance, Provenance::Defaulted);
        assert!((defaulted.charged_error(q_sum, num_reads) - 0.001).abs() < 1e-15);
    }

    /// A charge that could reach zero would reach negative infinity through the logarithm and
    /// turn a sample's whole row into `NaN`. The floor is what stops it.
    ///
    /// **And nothing stops it from above, which is a decision rather than an omission**: a
    /// ceiling binds on a single read and not on the fold of that read with others, and spec
    /// §2.3 forbids any term that is a non-linear function of a per-read quality. So a scale of
    /// `1e30` charges a great deal more than one — and §3.3's log-space form, which charges
    /// `q_sum + n·log scale`, is positive under exactly the same conditions.
    #[test]
    fn a_charge_is_floored_from_below_and_deliberately_not_capped_from_above() {
        let tiny = ReadGroupCalibration {
            scale: 1e-30,
            provenance: Provenance::Supplied,
        };
        let huge = ReadGroupCalibration {
            scale: 1e30,
            provenance: Provenance::Supplied,
        };
        let (q_sum, num_reads) = minted_reads_at_phred(&[30]);

        assert_eq!(tiny.charged_error(q_sum, num_reads), MIN_BASE_ERROR);
        assert!(tiny.charged_error(q_sum, num_reads).ln().is_finite());
        assert!(
            huge.charged_error(q_sum, num_reads) > 1.0,
            "an uncapped charge at a scale of 1e30 is {}",
            huge.charged_error(q_sum, num_reads)
        );
    }

    /// *Not identified* and *zero* are different answers, and a caller told "no
    /// contamination" would act on it. At one sample there is no panel to be surprised by,
    /// so the plain formula runs — absent, never a fitted zero.
    #[test]
    fn a_fraction_that_could_not_be_identified_is_absent_rather_than_zero() {
        let refused = ContaminationEstimate::NotIdentified {
            reason:
                crate::ng::parameter_estimation::joint::contamination::NotIdentifiedReason::NoPanel,
        };

        assert!(ContaminationView::of_estimate(&refused).is_none());
    }

    /// The two libraries that must not look alike: one measured over 4,102 markers and found
    /// clean, one that touched none. **Both carry a fraction of zero** — the search keeps
    /// zero where the likelihood is flat — so the value alone cannot tell them apart and the
    /// counts are the only thing that can.
    #[test]
    fn measured_clean_and_never_measured_carry_the_same_fraction_and_different_counts() {
        let clean = ContaminationEstimate::Estimated {
            alpha: 0.0,
            source: crate::ng::parameter_estimation::joint::contamination::ContaminationSource::ThisReadGroupsReads,
            panel_markers: 5_310,
            markers_with_reads: 4_102,
            reads_on_markers: 51_884,
            leverage: 0.031,
        };
        let unmeasurable = ContaminationEstimate::Estimated {
            alpha: 0.0,
            source: crate::ng::parameter_estimation::joint::contamination::ContaminationSource::ThisReadGroupsReads,
            panel_markers: 5_310,
            markers_with_reads: 0,
            reads_on_markers: 0,
            leverage: 0.031,
        };

        let clean = ContaminationView::of_estimate(&clean).expect("an estimate is a view");
        let unmeasurable =
            ContaminationView::of_estimate(&unmeasurable).expect("an estimate is a view");

        assert_eq!(clean.fraction, unmeasurable.fraction);
        assert!(clean.was_measured());
        assert!(!unmeasurable.was_measured());
    }

    /// **Both counts have to be non-zero, and the two fixtures above cannot show it** — in
    /// each of them the counts move together, where *and* and *or* agree. A review mutation
    /// swapping one for the other survived every test in this file.
    ///
    /// Neither half of the pair is reachable from the estimator today; they are here because
    /// the fields are public and the predicate is what a consumer will trust.
    #[test]
    fn a_view_missing_either_count_was_not_measured() {
        let markers_but_no_reads = ContaminationView {
            fraction: 0.0,
            markers_with_reads: 4_102,
            reads_on_markers: 0,
            source: ContaminationSource::ThisReadGroupsReads,
        };
        let reads_but_no_markers = ContaminationView {
            fraction: 0.0,
            markers_with_reads: 0,
            reads_on_markers: 51_884,
            source: ContaminationSource::ThisReadGroupsReads,
        };

        assert!(!markers_but_no_reads.was_measured());
        assert!(!reads_but_no_markers.was_measured());
    }

    /// The fraction and the two counts are read off the estimate's own fields, and all three
    /// are integers or a fraction that a transposition would carry silently — so the fixture
    /// gives each a value nothing else here has.
    #[test]
    fn the_view_reads_the_fraction_and_both_counts_off_the_estimate() {
        let estimate = ContaminationEstimate::Estimated {
            alpha: 0.0628,
            source: crate::ng::parameter_estimation::joint::contamination::ContaminationSource::ThisReadGroupsReads,
            panel_markers: 5_310,
            markers_with_reads: 4_102,
            reads_on_markers: 51_884,
            leverage: 0.031,
        };

        let view = ContaminationView::of_estimate(&estimate).expect("an estimate is a view");

        assert_eq!(view.fraction, 0.0628);
        assert_eq!(view.markers_with_reads, 4_102);
        assert_eq!(view.reads_on_markers, 51_884);
    }

    /// The staging buffer belongs to the scratch, so a worker holds one across every sample.
    /// Two samples through one scratch: the second must see its own rows and none of the
    /// first's.
    #[test]
    fn one_buffer_narrows_two_samples_without_mixing_them() {
        let mut buffer = GenericEvidenceBuffer::default();
        let first = [supported_row(0, 1, 11, -3.5), supported_row(1, 1, 4, -19.0)];
        let second = [supported_row(1, 1, 6, -8.25)];
        let mapping = [Some(AlleleId(0)), Some(AlleleId(1))];

        let (rows, dropped) = buffer.narrow_supported(&first, &mapping);
        assert_eq!(rows.len(), 2);
        assert_eq!(dropped, 0.0);

        let (rows, dropped) = buffer.narrow_supported(&second, &mapping);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].num_reads, 6);
        assert_eq!(dropped, 0.0);
    }

    /// **The dropped quality has to come back through the wrapper too**, and the test above
    /// cannot show it: its mapping keeps every allele, so `0.0` is the only answer it could
    /// give and asserting it is a transcription. Here selection drops merge allele 1, and its
    /// −19 is the quality of reads that are still evidence — they belong in the pooled
    /// leftover the data likelihood compares between loci.
    #[test]
    fn the_wrapper_carries_the_dropped_quality_out() {
        let mut buffer = GenericEvidenceBuffer::default();
        let rows = [
            supported_row(0, 1, 11, -3.5),
            supported_row(1, 1, 4, -19.0),
            supported_row(2, 1, 6, -8.25),
        ];
        let mapping = [Some(AlleleId(0)), None, Some(AlleleId(1))];

        let (kept, dropped) = buffer.narrow_supported(&rows, &mapping);

        assert_eq!(dropped, -19.0);
        assert_eq!(kept.len(), 2);
    }

    /// **The shape the row function needs, and the reason the evidence buffer is not the
    /// row's scratch.** A caller narrows, builds the view over what came back, and then hands
    /// the row its own scratch — which it can only do if the two are different objects,
    /// because the view is still borrowing the first while the row wants `&mut` the second.
    ///
    /// This test is a compile-time check wearing a runtime test's clothes: with the buffer and
    /// the row's scratch merged into one type, the third line is
    /// `error[E0499]: cannot borrow as mutable more than once`, which is where the next
    /// milestone would have discovered it.
    #[test]
    fn one_caller_holds_the_evidence_and_a_row_scratch_at_once() {
        /// Stands in for the row function's signature: reads the evidence, writes a row,
        /// scribbles in its own scratch.
        fn stub_row(
            evidence: &GenericSampleEvidence<'_>,
            out: &mut [f64],
            scratch: &mut SsrRowScratch<()>,
        ) {
            scratch.prepare_emissions(&SsrSampleEvidence::new(&[], &ssr_detail()), 0, f64::NAN);
            out[0] = f64::from(evidence.supported.len() as u32);
        }

        let mut buffer = GenericEvidenceBuffer::default();
        let mut row_scratch: SsrRowScratch<()> = SsrRowScratch::default();
        let rows = [supported_row(0, 1, 11, -3.5)];
        let mapping = [Some(AlleleId(0))];
        let mut out = [0.0];

        let (supported, dropped) = buffer.narrow_supported(&rows, &mapping);
        let evidence = GenericSampleEvidence::new(supported, dropped, &[]);
        stub_row(&evidence, &mut out, &mut row_scratch);

        assert_eq!(out[0], 1.0);
    }

    /// **Cleared, not resized.** A locus with as many `(observation, candidate)` pairs as the
    /// last one would otherwise reuse the last one's emissions — a wrong likelihood at every
    /// sample, and nothing failing. The fixture makes the two loci the same size on purpose,
    /// because that is the only case where the defect is invisible.
    #[test]
    fn the_emission_cache_carries_nothing_over_between_two_loci_of_one_size() {
        let three = [
            ssr_observation(b"ACACAC", ReadWitness::Complete, 3),
            ssr_observation(b"ACAC", ran_out_after(4), 2),
            ssr_observation(b"ACACACAC", ReadWitness::Complete, 5),
        ];
        let two = [
            ssr_observation(b"ACACAC", ReadWitness::Complete, 3),
            ssr_observation(b"ACAC", ran_out_after(4), 2),
        ];
        let detail = ssr_detail();
        let mut scratch: SsrRowScratch<()> = SsrRowScratch::default();

        scratch.prepare_emissions(&SsrSampleEvidence::new(&three, &detail), 2, f64::NAN);
        for (at, slot) in scratch.emissions_mut().iter_mut().enumerate() {
            *slot = at as f64 + 1.0;
        }
        assert_eq!(scratch.emission_count(), 6);

        scratch.prepare_emissions(&SsrSampleEvidence::new(&two, &detail), 3, -1.0);

        assert_eq!(scratch.emission_count(), 6);
        assert!(scratch.emissions_mut().iter().all(|&slot| slot == -1.0));
    }

    /// **The cache is sized over every observation, partials included**, because that is the
    /// axis the two filters' positions are on — both enumerate the whole slice and then
    /// filter. Sized over the complete ones alone it holds 8 slots here, and the complete
    /// observation at position 2 addresses slot 8, past the end.
    ///
    /// The fixture puts a partial *between* two complete observations on purpose: with the
    /// partials at the end, the complete half would write entirely in bounds and nothing would
    /// fail until the censored half ran.
    #[test]
    fn the_cache_is_sized_by_the_positions_the_filters_yield() {
        let observations = [
            ssr_observation(b"ACACACAC", ReadWitness::Complete, 12),
            ssr_observation(b"ACAC", ran_out_after(4), 5),
            ssr_observation(b"ACACAC", ReadWitness::Complete, 3),
        ];
        let detail = ssr_detail();
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        let mut scratch: SsrRowScratch<()> = SsrRowScratch::default();

        scratch.prepare_emissions(&evidence, 4, f64::NAN);

        assert_eq!(evidence.complete_observations().count(), 2);
        assert_eq!(scratch.emission_count(), 12);
        // Every position either filter yields addresses a slot the cache holds.
        for (at, _) in evidence
            .complete_observations()
            .chain(evidence.partial_observations())
        {
            scratch.set_emission(at, 3, at as f64);
            assert_eq!(scratch.emission_at(at, 3), at as f64);
        }
    }

    /// The index has one spelling, and a candidate past what the cache was prepared for is a
    /// caller bug rather than a neighbouring candidate's emission read silently.
    #[test]
    #[should_panic(expected = "is past the")]
    fn a_candidate_past_the_prepared_width_is_a_caller_bug() {
        let observations = [ssr_observation(b"ACAC", ReadWitness::Complete, 4)];
        let detail = ssr_detail();
        let mut scratch: SsrRowScratch<()> = SsrRowScratch::default();
        scratch.prepare_emissions(&SsrSampleEvidence::new(&observations, &detail), 2, 0.0);

        let _ = scratch.emission_at(0, 2);
    }

    /// An unwritten slot holds what the caller asked for. **Zero is the one value it must not
    /// silently become**: a slip the candidate cannot reach legitimately scores zero, so a
    /// cache pre-filled with zeros makes *never computed* and *computed as impossible* the
    /// same number. A `prepare_emissions` that ignored its argument would pass every other
    /// test here.
    #[test]
    fn an_unwritten_emission_slot_holds_the_value_the_caller_asked_for() {
        let observations = [
            ssr_observation(b"ACAC", ReadWitness::Complete, 4),
            ssr_observation(b"ACACAC", ReadWitness::Complete, 2),
        ];
        let detail = ssr_detail();
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        let mut scratch: SsrRowScratch<()> = SsrRowScratch::default();

        scratch.prepare_emissions(&evidence, 2, f64::NAN);
        assert!(scratch.emissions_mut().iter().all(|slot| slot.is_nan()));

        scratch.prepare_emissions(&evidence, 2, -1.0);
        assert!(scratch.emissions_mut().iter().all(|&slot| slot == -1.0));
    }

    /// A locus with no observations or no candidates sizes the cache to nothing rather than
    /// leaving the previous locus's entries reachable.
    #[test]
    fn an_empty_locus_leaves_an_empty_emission_cache() {
        let observations = [ssr_observation(b"ACAC", ReadWitness::Complete, 4)];
        let detail = ssr_detail();
        let mut scratch: SsrRowScratch<()> = SsrRowScratch::default();

        scratch.prepare_emissions(&SsrSampleEvidence::new(&observations, &detail), 4, 1.0);
        assert_eq!(scratch.emission_count(), 4);

        scratch.prepare_emissions(&SsrSampleEvidence::new(&[], &detail), 4, 1.0);
        assert_eq!(scratch.emission_count(), 0);

        scratch.prepare_emissions(&SsrSampleEvidence::new(&observations, &detail), 0, 1.0);
        assert_eq!(scratch.emission_count(), 0);
    }

    /// The model's own scratch — its placement buffer and its alignment matrix — is the other
    /// half of what the caller holds, and it has to survive a locus so a worker allocates once.
    /// Nothing else here touches it, and deleting the accessor leaves the tests green.
    #[test]
    fn the_models_own_scratch_is_reachable_and_survives_a_locus() {
        let detail = ssr_detail();
        let mut scratch: SsrRowScratch<Vec<u8>> = SsrRowScratch::default();

        scratch.model_scratch_mut().extend_from_slice(b"ACGT");
        scratch.prepare_emissions(&SsrSampleEvidence::new(&[], &detail), 3, 0.0);

        assert_eq!(scratch.model_scratch_mut().as_slice(), b"ACGT");
    }

    // ---- C1: the contamination mixture's two halves ----

    fn a_read_group_contaminated_at(fraction: f64) -> ContaminationView {
        ContaminationView {
            fraction,
            markers_with_reads: 1_000,
            reads_on_markers: 9_000,
            source: ContaminationSource::ThisReadGroupsReads,
        }
    }

    /// **Nothing contaminated reads as nothing contaminated from both accessors**, which is
    /// what makes the row's single code path spec §3.3: `1 − 0` times the sample's own term,
    /// plus `0` times whatever the contaminant would have shown.
    #[test]
    fn no_mixture_gives_a_zero_fraction_and_a_zero_frequency() {
        let mixture = ContaminationMixture::uncontaminated();

        assert_eq!(mixture.allele_count(), 0);
        assert_eq!(mixture.fraction_of(ReadGroupId(7)), 0.0);
        assert_eq!(
            mixture.contaminant_frequency_of(ReadGroupId(7), AlleleId(3)),
            0.0
        );
    }

    /// The three counts a mixture answers about its own shape, including the batch count,
    /// which nothing else in the crate calls yet.
    #[test]
    fn a_mixture_reports_its_own_shape() {
        let fractions = [
            a_read_group_contaminated_at(0.01),
            a_read_group_contaminated_at(0.02),
            a_read_group_contaminated_at(0.03),
        ];
        let batching = [BatchId(0), BatchId(1), BatchId(1)];
        let frequencies = [0.8, 0.2, 0.4, 0.6];
        let mixture =
            ContaminationMixture::new(&fractions, BatchOfEachReadGroup(&batching), &frequencies, 2);

        assert_eq!(mixture.read_group_count(), 3);
        assert_eq!(mixture.allele_count(), 2);
        assert_eq!(mixture.batch_count(), 2);
        assert!(!mixture.is_absent());

        let absent = ContaminationMixture::uncontaminated();
        assert_eq!(absent.batch_count(), 0);
        assert!(absent.is_absent());
    }

    /// **The batch-major stride, at a shape where transposing it would still be legal.** Three
    /// batches and three alleles make the table square, so reading it as allele-major returns a
    /// real frequency from the wrong batch rather than running off the end — the failure the
    /// square case is the only one that can hide.
    #[test]
    fn the_frequency_table_is_batch_major_and_a_square_table_says_so() {
        let fractions = [
            a_read_group_contaminated_at(0.01),
            a_read_group_contaminated_at(0.01),
            a_read_group_contaminated_at(0.01),
        ];
        let batching = [BatchId(0), BatchId(1), BatchId(2)];
        // Batch 0 favours the reference, batch 1 the first alternative, batch 2 the second.
        let frequencies = [0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8];
        let mixture =
            ContaminationMixture::new(&fractions, BatchOfEachReadGroup(&batching), &frequencies, 3);

        for batch in 0..3u32 {
            for allele in 0..3u16 {
                let expected = if u32::from(allele) == batch { 0.8 } else { 0.1 };
                assert_eq!(
                    mixture.contaminant_frequency_of(ReadGroupId(batch), AlleleId(allele)),
                    expected,
                    "batch {batch}, allele {allele}"
                );
            }
        }
    }

    #[test]
    fn a_mixture_hands_back_the_fraction_and_frequency_it_was_given() {
        let fractions = [
            a_read_group_contaminated_at(0.0),
            a_read_group_contaminated_at(0.04),
        ];
        let frequencies = [0.7, 0.25, 0.05];
        let mixture = one_batch(&fractions, &frequencies);

        assert_eq!(mixture.allele_count(), 3);
        assert_eq!(mixture.fraction_of(ReadGroupId(1)), 0.04);
        assert_eq!(
            mixture.contaminant_frequency_of(ReadGroupId(1), AlleleId(2)),
            0.05
        );
    }

    /// **A fraction of exactly one is the one bound that does arithmetic**: `1 − c` underflows
    /// to zero, and an allele the contaminant never shows would take the mixture to zero and
    /// its logarithm to negative infinity — the one input that could, since every other term
    /// is floored or is a copy share.
    #[test]
    #[should_panic(expected = "is not a sample of this one")]
    fn a_library_that_is_entirely_somebody_else_is_refused() {
        let fractions = [a_read_group_contaminated_at(1.0)];
        let frequencies = [1.0];

        let _ = one_batch(&fractions, &frequencies);
    }

    #[test]
    #[should_panic(expected = "is not a sample of this one")]
    fn a_negative_fraction_is_refused() {
        let fractions = [a_read_group_contaminated_at(-0.01)];
        let frequencies = [1.0];

        let _ = one_batch(&fractions, &frequencies);
    }

    #[test]
    #[should_panic(expected = "which is not a frequency")]
    fn a_frequency_past_one_is_refused() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let frequencies = [0.4, 1.6];

        let _ = one_batch(&fractions, &frequencies);
    }

    /// **The lower half of the same bound, which had no test until C1's review widened the
    /// range to `NEG_INFINITY..=1.0` and watched the whole module stay green.**
    ///
    /// `new` is public and this assertion is the only guard. A negative frequency reaches the
    /// row as `c · q(o)` below zero, which can drive the mixture negative — and `ln` of a
    /// negative number is `NaN`, which makes every comparison in an argmax false, so a
    /// genotype is picked with nothing to say it was picked from garbage. Measured under the
    /// widened range at a frequency of −1.0: a row of `[NaN, NaN, -inf]`.
    #[test]
    #[should_panic(expected = "which is not a frequency")]
    fn a_negative_frequency_is_refused() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let frequencies = [1.4, -1.0];

        let _ = one_batch(&fractions, &frequencies);
    }

    #[test]
    #[should_panic(expected = "which is not a frequency")]
    fn a_frequency_that_is_not_a_number_is_refused() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let frequencies = [0.5, f64::NAN];

        let _ = one_batch(&fractions, &frequencies);
    }

    /// Half a mixture is not a mixture: fractions with no frequencies would read every allele
    /// as one the contaminant never shows, and frequencies with no fractions would be a table
    /// nothing multiplies. Either is a caller that built one half and forgot the other.
    #[test]
    #[should_panic(expected = "needs both halves or neither")]
    fn fractions_without_frequencies_are_refused() {
        let fractions = [a_read_group_contaminated_at(0.02)];

        let _ = ContaminationMixture::new(
            &fractions,
            BatchOfEachReadGroup(&ALL_IN_ONE_BATCH[..1]),
            &[],
            1,
        );
    }

    #[test]
    #[should_panic(expected = "needs both halves or neither")]
    fn frequencies_without_fractions_are_refused() {
        let frequencies = [0.5, 0.5];

        let _ = ContaminationMixture::new(&[], BatchOfEachReadGroup(&[]), &frequencies, 2);
    }

    /// **One named way to say uncontaminated.** Two empty slices satisfy every other check in
    /// `new`, so without this a caller building both halves from a fit that returned nothing
    /// lands on the clean case through a constructor whose name says the opposite — and never
    /// reads the paragraph where the decision is written down (C1's review).
    #[test]
    #[should_panic(expected = "is spelled `ContaminationMixture::uncontaminated`")]
    fn two_empty_halves_are_refused_because_that_case_has_a_name() {
        let _ = ContaminationMixture::new(&[], BatchOfEachReadGroup(&[]), &[], 1);
    }

    /// The two accessors' own panics, which the row's eager checks make unreachable from a
    /// run — so this is the only place they are exercised.
    #[test]
    #[should_panic(expected = "has no contamination fraction")]
    fn a_read_group_the_mixture_does_not_cover_is_a_caller_bug() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let frequencies = [0.9, 0.1];

        let _ = one_batch(&fractions, &frequencies).fraction_of(ReadGroupId(4));
    }

    #[test]
    #[should_panic(expected = "has no contaminant frequency")]
    fn an_allele_the_mixture_does_not_cover_is_a_caller_bug() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let frequencies = [0.9, 0.1];

        let _ = one_batch(&fractions, &frequencies)
            .contaminant_frequency_of(ReadGroupId(0), AlleleId(5));
    }

    // ---- D1: the row's own scratch ----

    /// **The stride assertions are the whole reason `slot` exists**, and removing both left the
    /// module green (D1's review): a cache prepared for one locus and read at another's stride
    /// returns a real verdict from the wrong row rather than running off the end.
    #[test]
    #[should_panic(expected = "allele 4 is past the 3")]
    fn the_scratch_refuses_an_allele_past_the_stride_it_was_prepared_for() {
        let mut scratch = GenericRowScratch::default();
        scratch.prepare_emissions(0, 2, 3);

        let _ = scratch.is_compatible(0, AlleleId(4));
    }

    #[test]
    #[should_panic(expected = "partial 5 is past the 2")]
    fn the_scratch_refuses_a_partial_past_the_rows_it_was_prepared_for() {
        let mut scratch = GenericRowScratch::default();
        scratch.prepare_emissions(0, 2, 3);

        let _ = scratch.is_compatible(5, AlleleId(0));
    }

    /// **One scratch serves every sample of a locus, and every locus of a worker** — which is
    /// what it is for, and which no call site exercised, since all of them stood up a fresh one.
    ///
    /// Preparing it again resizes it and clears what the previous sample wrote, so a sample with
    /// fewer partials than the one before cannot read a verdict left by it.
    #[test]
    fn one_scratch_is_reused_across_samples_and_forgets_the_last() {
        let mut scratch = GenericRowScratch::default();

        scratch.prepare_emissions(0, 3, 2);
        assert_eq!(scratch.verdict_count(), 6);
        for partial in 0..3 {
            for allele in 0..2u16 {
                scratch.set_compatible(partial, AlleleId(allele), true);
            }
        }

        // A narrower sample at the same locus.
        scratch.prepare_emissions(0, 1, 2);
        assert_eq!(scratch.verdict_count(), 2);
        assert!(
            !scratch.is_compatible(0, AlleleId(0)) && !scratch.is_compatible(0, AlleleId(1)),
            "preparing again must forget what the last sample wrote"
        );

        // And a wider locus afterwards, which grows it.
        scratch.prepare_emissions(0, 2, 5);
        assert_eq!(scratch.verdict_count(), 10);
        assert!((0..2).all(|partial| {
            (0..5u16).all(|allele| !scratch.is_compatible(partial, AlleleId(allele)))
        }));
    }

    /// A sample with no partial reads prepares an empty cache rather than keeping the last
    /// sample's — the case a locus where only some samples ran out reaches at every worker.
    #[test]
    fn a_sample_with_no_partials_prepares_an_empty_cache() {
        let mut scratch = GenericRowScratch::default();

        scratch.prepare_emissions(0, 2, 3);
        scratch.set_compatible(1, AlleleId(2), true);
        scratch.prepare_emissions(0, 0, 3);

        assert_eq!(scratch.verdict_count(), 0);
    }

    /// **The charged errors are held on the same buffer as the verdicts, and the two halves of
    /// it are addressed by separate accessors** — the partials sit above the complete
    /// observations, so a caller adding the offset by hand would read a real value from the
    /// wrong half whenever it got the base wrong.
    #[test]
    fn the_two_halves_of_the_emission_cache_are_addressed_apart() {
        let mut scratch = GenericRowScratch::default();
        scratch.prepare_emissions(2, 3, 2);

        assert_eq!(scratch.supported_count(), 2);
        assert_eq!(scratch.partial_count(), 3);
        assert_eq!(scratch.allele_count(), 2);

        for at in 0..2 {
            scratch.set_supported_charged_error(at, 0.1 * (at + 1) as f64);
        }
        for at in 0..3 {
            scratch.set_partial_charged_error(at, 10.0 + at as f64);
        }
        assert_eq!(scratch.charged_error_of_supported(0), 0.1);
        assert_eq!(scratch.charged_error_of_supported(1), 0.2);
        assert_eq!(scratch.charged_error_of_partial(0), 10.0);
        assert_eq!(scratch.charged_error_of_partial(2), 12.0);
    }

    /// An unwritten charged error is a `NaN` and not a zero, and the direction matters: **zero
    /// is a real value** — a library whose reads are never wrong — so a slot the fill skipped
    /// would take every genotype that fails to explain the observation to `−∞` and read as
    /// evidence.
    #[test]
    fn an_unfilled_charged_error_is_not_a_number() {
        let mut scratch = GenericRowScratch::default();
        scratch.prepare_emissions(2, 1, 2);
        assert!(scratch.charged_error_of_supported(1).is_nan());
        assert!(scratch.charged_error_of_partial(0).is_nan());

        scratch.set_supported_charged_error(1, 0.5);
        scratch.prepare_emissions(2, 1, 2);
        assert!(
            scratch.charged_error_of_supported(1).is_nan(),
            "preparing again must forget what the last sample wrote"
        );
    }

    /// A complete observation past the cache's own count is refused rather than read out of the
    /// partials above it — which is exactly what an assembly paired with another sample's fill
    /// would reach.
    #[test]
    #[should_panic(expected = "complete observation 2 is past the 2")]
    fn a_complete_observation_past_the_prepared_cache_is_a_caller_bug() {
        let mut scratch = GenericRowScratch::default();
        scratch.prepare_emissions(2, 3, 2);
        let _ = scratch.charged_error_of_supported(2);
    }

    /// And the same on the writing side, which is the door a fill comes through.
    #[test]
    #[should_panic(expected = "complete observation 5 is past the 1")]
    fn writing_a_complete_observation_past_the_prepared_cache_is_a_caller_bug() {
        let mut scratch = GenericRowScratch::default();
        scratch.prepare_emissions(1, 0, 2);
        scratch.set_supported_charged_error(5, 0.1);
    }

    /// A partial past the cache's own count runs off the end of the buffer rather than into a
    /// neighbour, and is refused by name.
    #[test]
    #[should_panic(expected = "partial 3 is past the 3")]
    fn a_partial_past_the_prepared_cache_is_a_caller_bug() {
        let mut scratch = GenericRowScratch::default();
        scratch.prepare_emissions(2, 3, 2);
        let _ = scratch.charged_error_of_partial(3);
    }

    // ---- C2: the frequency the contaminant is drawn against ----

    /// Four diploids at a two-allele locus: two carrying the reference twice, one
    /// heterozygous, one carrying the alternative twice.
    const FOUR_DIPLOIDS: [f64; 8] = [
        2.0, 0.0, // sample 0 — reference homozygote
        1.0, 1.0, // sample 1 — heterozygote
        0.0, 2.0, // sample 2 — alternative homozygote
        2.0, 0.0, // sample 3 — reference homozygote
    ];

    /// **Both stages, as the caller's loop will run them**: sum each batch's copies once for the
    /// locus, then leave one sample out of its own batch and normalise. Returns how many batches
    /// had nothing left to read a frequency off.
    fn frequencies_leaving_out(
        expected_copies_by_sample: &[f64],
        batch_of_each_sample: &[BatchId],
        allele_count: usize,
        sample: usize,
        out: &mut [f64],
    ) -> usize {
        let mut batch_copies = vec![f64::NAN; out.len()];
        fill_batch_allele_copies(
            expected_copies_by_sample,
            BatchOfEachSample(batch_of_each_sample),
            allele_count,
            &mut batch_copies,
        );
        fill_contaminant_allele_frequencies(
            &batch_copies,
            SampleAlleleCopies::new(
                &expected_copies_by_sample[sample * allele_count..(sample + 1) * allele_count],
            ),
            batch_of_each_sample[sample],
            allele_count,
            out,
        )
    }

    /// **The default batching gives the frequency over the rest of the cohort**, which is what
    /// makes it lose nothing: one batch, and every sample but the one being scored.
    ///
    /// Across all four samples the locus holds five reference copies and three alternative.
    /// Scoring sample 0 — a reference homozygote — takes its two reference copies back out,
    /// leaving three of each: **the contaminant is drawn against 1 in 2, not the cohort's 5 in
    /// 8**, because a contaminant is somebody else and this sample is not somebody else.
    #[test]
    fn one_batch_gives_the_frequency_over_the_rest_of_the_cohort() {
        let mut out = [f64::NAN; 2];

        let no_evidence =
            frequencies_leaving_out(&FOUR_DIPLOIDS, &[BatchId::ALL_TOGETHER; 4], 2, 0, &mut out);

        assert_eq!(out, [0.5, 0.5]);
        assert_eq!(no_evidence, 0);
        // And it is a distribution, which every row has to be for `c · q` to be a probability.
        assert!((out.iter().sum::<f64>() - 1.0).abs() < 1e-15);
    }

    /// **Leaving the sample out is not a rounding correction — it moves the frequency by how
    /// much of its batch the sample is.** Scoring sample 1, the heterozygote, against the whole
    /// cohort would give 3 in 8 for the alternative; against everyone else it is 2 in 6.
    ///
    /// At 63 samples that difference is one part in 63 and would not decide a call. **At four
    /// it is 0.042 in absolute frequency, and at one sample in a batch it is everything** —
    /// which is why this is the subtraction and not a nicety (owner, 2026-08-24).
    #[test]
    fn leaving_the_sample_out_moves_the_frequency_by_its_own_share() {
        let mut out = [f64::NAN; 2];

        frequencies_leaving_out(&FOUR_DIPLOIDS, &[BatchId::ALL_TOGETHER; 4], 2, 1, &mut out);

        let over_everyone_else = out[1];
        let over_the_whole_cohort = 3.0 / 8.0;
        assert!((over_everyone_else - 2.0 / 6.0).abs() < 1e-15);
        assert!(
            (over_everyone_else - over_the_whole_cohort).abs() > 0.04,
            "leaving one of four samples out moved the alternative's frequency from \
             {over_the_whole_cohort} to {over_everyone_else}"
        );
    }

    /// **Two batches at one locus, and the samples in each answer for their own.** Scoring
    /// sample 1 — the heterozygote, in the first batch — its own batch has only a reference
    /// homozygote left, while the second batch is untouched and half alternative.
    ///
    /// A contaminant is a neighbour on the same run, so a library that ran beside a reference
    /// homozygote must not be told the other batch's number, nor the whole cohort's.
    #[test]
    fn two_batches_at_one_locus_give_two_different_frequencies() {
        // Samples 0 and 1 ran together; samples 2 and 3 ran together.
        let batches = [BatchId(0), BatchId(0), BatchId(1), BatchId(1)];
        let mut out = [f64::NAN; 4];

        frequencies_leaving_out(&FOUR_DIPLOIDS, &batches, 2, 1, &mut out);

        // Batch 0, less the sample being scored, is one reference homozygote; batch 1 is a
        // reference homozygote and an alternative homozygote, untouched.
        assert_eq!(out[0], 1.0);
        assert_eq!(out[1], MIN_CONTAMINANT_FREQUENCY);
        assert_eq!(out[2..], [0.5, 0.5]);
    }

    /// **An allele its own batch never shows is floored, never zeroed.** Every other sample of
    /// this batch is a reference homozygote, so the alternative has no copies among the
    /// neighbours at all — and a frequency of exactly zero would be the claim that they
    /// *cannot* carry it, which two samples cannot establish.
    ///
    /// The floor is [`MIN_CONTAMINANT_FREQUENCY`] and it is deliberately tiny; the constant's
    /// own documentation says why a statistical size would be the wrong direction to be wrong
    /// in, and [`the_frequency_floor_cannot_outweigh_a_misread`] measures what that buys.
    #[test]
    fn an_allele_the_batch_never_shows_is_floored_rather_than_zeroed() {
        let three_reference_homozygotes = [2.0, 0.0, 2.0, 0.0, 2.0, 0.0];
        let mut out = [f64::NAN; 2];

        frequencies_leaving_out(
            &three_reference_homozygotes,
            &[BatchId::ALL_TOGETHER; 3],
            2,
            0,
            &mut out,
        );

        assert_eq!(out[0], 1.0);
        assert_eq!(out[1], MIN_CONTAMINANT_FREQUENCY);
        assert!(out[1] > 0.0, "a floored frequency is still a frequency");
    }

    /// **The floor's whole job is to stay out of the way**, and this is the measurement its
    /// documentation rests on. At a 3% contamination fraction a floored allele contributes
    /// `3 × 10⁻¹⁴` to the mixture, where the route it sits beside — the read simply being a
    /// misread — contributes `ε̄/m`, which is `3 × 10⁻⁵` even at a middling Phred 40. **Nine
    /// orders of magnitude below**, so a floored allele cannot make a read look like
    /// contamination.
    #[test]
    fn the_frequency_floor_cannot_outweigh_a_misread() {
        let from_the_contaminant = 0.03 * MIN_CONTAMINANT_FREQUENCY;
        let from_a_misread_at_phred_40 = 1e-4 / 3.0;

        // A billion times, measured — the bound is set two orders under that so it pins the
        // claim without pinning the arithmetic's last digits.
        assert!(
            from_a_misread_at_phred_40 / from_the_contaminant > 1e7,
            "the misread route is {} times the floored contaminant route, and the floor is \
             only defensive if that is large",
            from_a_misread_at_phred_40 / from_the_contaminant
        );
    }

    /// **A sample alone in its batch has no neighbours, and now says so.** Leaving it out of
    /// its own batch leaves nothing, so the row takes the no-evidence fallback: the reference,
    /// and the floor elsewhere.
    ///
    /// **Before the subtraction it returned its own genotype** — a lone heterozygote came back
    /// `[0.5, 0.5]`, so its own alternative reads were half-explained as a neighbour's, which
    /// is a sample acting as its own contaminant. C2's review found it; the owner settled it
    /// the same day. This is the conservative answer and the right one: a library with nobody
    /// beside it has no contaminating population to be drawn from.
    #[test]
    fn a_sample_alone_in_its_batch_has_no_neighbours_and_gets_the_reference() {
        let one_heterozygote = [1.0, 1.0];
        let mut out = [f64::NAN; 2];

        let no_evidence =
            frequencies_leaving_out(&one_heterozygote, &[BatchId::ALL_TOGETHER], 2, 0, &mut out);

        assert_eq!(no_evidence, 1);
        assert_eq!(
            out,
            [1.0 - MIN_CONTAMINANT_FREQUENCY, MIN_CONTAMINANT_FREQUENCY]
        );
    }

    /// **A batch is divided by its own copies, not by how many samples are in it** — and the
    /// fixtures that came before this one were blind to the difference, because every sample in
    /// them carried exactly two copies. Mutation testing found it: a divisor of *two copies per
    /// sample* passed all 127 module tests (C2's review).
    ///
    /// This cohort mixes ploidies, which is what separates the two. Scoring the alternative
    /// homozygote leaves a tetraploid at 3:1 beside a diploid reference homozygote — five
    /// reference copies and one alternative out of six, `[0.833, 0.167]`. Counting samples
    /// instead would divide by four and give `[1.25, 0.25]`, clamped to `[1.0, 0.25]`, **a row
    /// summing to 1.25**.
    #[test]
    fn a_batch_is_divided_by_its_copies_and_not_by_its_sample_count() {
        let mixed_ploidy = [
            3.0, 1.0, // a tetraploid carrying one alternative copy
            2.0, 0.0, // a diploid reference homozygote
            0.0, 2.0, // a diploid alternative homozygote — the one being scored
        ];
        let mut out = [f64::NAN; 2];

        frequencies_leaving_out(&mixed_ploidy, &[BatchId::ALL_TOGETHER; 3], 2, 2, &mut out);

        assert!((out[0] - 5.0 / 6.0).abs() < 1e-15 && (out[1] - 1.0 / 6.0).abs() < 1e-15);
        assert!(
            (out.iter().sum::<f64>() - 1.0).abs() < 1e-15,
            "a frequency table's row is a distribution, and this one sums to {}",
            out.iter().sum::<f64>()
        );
    }

    /// A batch whose remaining samples hold no copies at all has no frequency to read off, and
    /// dividing by its zero total would give `NaN` — which the row would take a logarithm of,
    /// and a `NaN` makes every comparison in an argmax false, so a genotype would be picked in
    /// silence. It gets the reference instead.
    #[test]
    fn a_batch_with_no_coverage_gets_the_reference_rather_than_a_nan() {
        let nothing_sequenced = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut out = [f64::NAN; 3];

        let no_evidence = frequencies_leaving_out(
            &nothing_sequenced,
            &[BatchId::ALL_TOGETHER; 2],
            3,
            0,
            &mut out,
        );

        assert_eq!(no_evidence, 1);
        assert!(out.iter().all(|f| f.is_finite()));
        assert_eq!(out[1], MIN_CONTAMINANT_FREQUENCY);
        assert_eq!(out[2], MIN_CONTAMINANT_FREQUENCY);
        assert!((out.iter().sum::<f64>() - 1.0).abs() < 1e-15);
    }

    /// The loop's copies are expectations over genotype posteriors, not counts, so they need
    /// not be whole — and the frequency is still the ratio.
    #[test]
    fn fractional_copies_are_a_frequency_like_any_other() {
        let two_uncertain_samples = [1.5, 0.5, 0.25, 1.75];
        let mut out = [f64::NAN; 2];

        frequencies_leaving_out(
            &two_uncertain_samples,
            &[BatchId::ALL_TOGETHER; 2],
            2,
            0,
            &mut out,
        );

        // Only the second sample is left: a quarter of a reference copy against one and
        // three-quarters of the alternative.
        assert!((out[0] - 0.125).abs() < 1e-15 && (out[1] - 0.875).abs() < 1e-15);
    }

    /// How many batches had nothing to read a frequency off, which is otherwise invisible: the
    /// fallback row is a legal frequency, and at one candidate allele it is bit-identical to a
    /// batch full of evidence.
    #[test]
    fn the_fill_says_how_many_batches_had_no_evidence() {
        // Samples 0 and 1 ran together; sample 2 ran on its own.
        let copies = [2.0, 0.0, 1.0, 1.0, 0.0, 2.0];
        let batches = [BatchId(0), BatchId(0), BatchId(1)];
        let mut out = [f64::NAN; 4];

        // Scoring sample 0 leaves a heterozygote in its batch and sample 2 alone in the other.
        assert_eq!(
            frequencies_leaving_out(&copies, &batches, 2, 0, &mut out),
            0
        );
        assert_eq!(out[..2], [0.5, 0.5]);

        // Scoring sample 2 empties its own batch, since nobody else was sequenced with it.
        assert_eq!(
            frequencies_leaving_out(&copies, &batches, 2, 2, &mut out),
            1
        );
        assert_eq!(out[..2], [0.75, 0.25]);
        assert_eq!(out[2], 1.0 - MIN_CONTAMINANT_FREQUENCY);
    }

    /// **The two laws this producer obeys, over cohorts nobody chose.**
    ///
    /// Every row is a distribution — the frequencies of one batch sum to one, up to whatever
    /// the never-seen floor lifted — and **the answer depends only on the ratios**, so
    /// multiplying every copy in the cohort by the same positive number changes nothing. The
    /// second is what says this is a frequency and not a count, and it is the law a divisor
    /// taken from the wrong place breaks.
    ///
    /// Hand-built fixtures pin values; this pins the shape of the function over 3 to 20
    /// samples, 1 to 5 alleles and 1 to 4 batches, which is the range around the corner the
    /// benchmarks sit in.
    #[test]
    fn every_row_is_a_distribution_and_only_the_ratios_matter() {
        use proptest::prelude::*;

        proptest!(|(
            sample_count in 3usize..20,
            allele_count in 1usize..5,
            batch_count in 1usize..4,
            scored in 0usize..3,
            copies in proptest::collection::vec(0.0f64..8.0, 3 * 5 * 20),
            scale in 0.25f64..40.0,
        )| {
            // Every declared batch must hold a sample, so deal them round-robin.
            let batches: Vec<BatchId> = (0..sample_count)
                .map(|sample| BatchId((sample % batch_count) as u32))
                .collect();
            let copies = &copies[..sample_count * allele_count];

            let mut out = vec![f64::NAN; batch_count * allele_count];
            frequencies_leaving_out(copies, &batches, allele_count, scored, &mut out);

            for row in out.chunks_exact(allele_count) {
                let total: f64 = row.iter().sum();
                prop_assert!(
                    (total - 1.0).abs() < 1e-9,
                    "a batch's frequencies sum to {total}, and a row is a distribution"
                );
                prop_assert!(row.iter().all(|&f| (0.0..=1.0).contains(&f)));
            }

            let scaled: Vec<f64> = copies.iter().map(|c| c * scale).collect();
            let mut from_scaled = vec![f64::NAN; batch_count * allele_count];
            frequencies_leaving_out(&scaled, &batches, allele_count, scored, &mut from_scaled);
            for (plain, scaled) in out.iter().zip(&from_scaled) {
                prop_assert!(
                    (plain - scaled).abs() < 1e-12,
                    "scaling every copy by {scale} moved a frequency from {plain} to {scaled}"
                );
            }
        });
    }

    // ---- C2: what each half of the producer refuses ----

    /// A negative copy count is arithmetic that went wrong upstream, and the finiteness check
    /// alone does not catch it — mutation testing dropped the `>= 0.0` half and the module
    /// stayed green, because every other fixture that reaches this guard uses a `NaN`.
    #[test]
    #[should_panic(expected = "copy count is finite and at or above zero")]
    fn a_negative_copy_count_is_a_caller_bug() {
        let mut out = [f64::NAN; 2];

        fill_batch_allele_copies(
            &[2.0, -1.0],
            BatchOfEachSample(&[BatchId::ALL_TOGETHER]),
            2,
            &mut out,
        );
    }

    #[test]
    #[should_panic(expected = "copy count is finite and at or above zero")]
    fn a_copy_count_that_is_not_a_number_is_a_caller_bug() {
        let mut out = [f64::NAN; 2];

        fill_batch_allele_copies(
            &[2.0, f64::NAN],
            BatchOfEachSample(&[BatchId::ALL_TOGETHER]),
            2,
            &mut out,
        );
    }

    /// A copy table that is not a whole number of rows would leave the last batch short, and
    /// `chunks_exact_mut` drops the remainder — so the row would keep whatever it held. The
    /// check is release-mode and had no test until mutation testing neutered it (C2's review).
    #[test]
    #[should_panic(expected = "a batch's row is 2 alleles wide")]
    fn a_copy_buffer_that_is_not_whole_rows_is_a_caller_bug() {
        let mut out = [f64::NAN; 3];

        fill_batch_allele_copies(
            &[2.0, 0.0],
            BatchOfEachSample(&[BatchId::ALL_TOGETHER]),
            2,
            &mut out,
        );
    }

    #[test]
    #[should_panic(expected = "which is not a whole cohort")]
    fn copies_that_are_not_one_row_per_sample_are_a_caller_bug() {
        let mut out = [f64::NAN; 2];

        fill_batch_allele_copies(
            &[2.0, 0.0, 1.0],
            BatchOfEachSample(&[BatchId::ALL_TOGETHER; 2]),
            2,
            &mut out,
        );
    }

    #[test]
    #[should_panic(expected = "and the run declared 1")]
    fn a_sample_in_a_batch_the_run_did_not_declare_is_a_caller_bug() {
        let mut out = [f64::NAN; 2];

        fill_batch_allele_copies(
            &[2.0, 0.0, 1.0, 1.0],
            BatchOfEachSample(&[BatchId(0), BatchId(1)]),
            2,
            &mut out,
        );
    }

    /// **A batch nobody was sequenced in is not a thin batch, and without this check the two
    /// are the same row.** The sum never writes a batch no sample names, so `out.fill(0.0)`
    /// leaves it at a zero total and it leaves through the frequency's no-evidence fallback —
    /// coming back as a legal-looking frequency that says the neighbours carry the reference.
    /// C2's review built it: a table sized for three batches against a batching naming only the
    /// first returned `[0.999999999999, 1e-12]` for the other two, with nothing raised.
    #[test]
    #[should_panic(expected = "no sample ran in batch 1")]
    fn a_declared_batch_nobody_was_sequenced_in_is_a_caller_bug() {
        let mut out = [f64::NAN; 4];

        fill_batch_allele_copies(
            &[2.0, 0.0],
            BatchOfEachSample(&[BatchId::ALL_TOGETHER]),
            2,
            &mut out,
        );
    }

    /// **Every slot is checked as it accumulates and the total is not**, so a row of finite
    /// copies whose sum overflows would reach the division — where `finite / inf` is zero, every
    /// allele is lifted to the floor, and the row comes back finite and plausible and says the
    /// neighbours carry nothing. C2's review built it with one sample holding `1e308` copies of
    /// each of two alleles.
    #[test]
    #[should_panic(expected = "arithmetic that went wrong upstream")]
    fn a_batch_whose_copies_overflow_their_sum_is_a_caller_bug() {
        let mut out = [f64::NAN; 2];

        frequencies_leaving_out(
            &[1e308, 1e308, 1e308, 1e308],
            &[BatchId::ALL_TOGETHER; 2],
            2,
            0,
            &mut out,
        );
    }

    /// A sample that claims more copies than the batch it is one addend of means the two count
    /// paths have gone out of step — the check `fill_sample_concentration` makes for the same
    /// subtraction, at the same threshold.
    /// **Debug-profile only** — `fill_contaminant_allele_frequencies`' own `# Panics` says so
    /// in as many words (*"**In debug**, on a sample holding materially more copies than the
    /// batch it is part of"*), mirroring `fill_sample_concentration`'s check for the same
    /// subtraction. See [`rows_out_of_pair_order_are_refused_in_debug`] for why the attribute
    /// is here and what it must never be used for.
    #[test]
    #[should_panic(expected = "count paths have gone out of step")]
    #[cfg(debug_assertions)]
    fn a_sample_carrying_more_than_its_whole_batch_is_refused_in_debug() {
        let batch_copies = [2.0, 0.0];
        let impossible = [4.0, 0.0];
        let mut out = [f64::NAN; 2];

        let _ = fill_contaminant_allele_frequencies(
            &batch_copies,
            SampleAlleleCopies::new(&impossible),
            BatchId::ALL_TOGETHER,
            2,
            &mut out,
        );
    }

    #[test]
    #[should_panic(expected = "the copy table holds 1 batches")]
    fn a_sample_in_a_batch_the_copy_table_lacks_is_a_caller_bug() {
        let batch_copies = [2.0, 0.0];
        let own = [1.0, 0.0];
        let mut out = [f64::NAN; 2];

        let _ = fill_contaminant_allele_frequencies(
            &batch_copies,
            SampleAlleleCopies::new(&own),
            BatchId(1),
            2,
            &mut out,
        );
    }

    // ---- C2: what the mixture refuses ----

    /// A mixture whose batching names a batch its frequency table has no row for would read
    /// another batch's frequencies, or past the end. Checked at construction, so the row's
    /// per-observation lookup cannot reach it.
    #[test]
    #[should_panic(expected = "this mixture holds 1 batches")]
    fn a_batching_past_the_frequency_table_is_a_caller_bug() {
        let fractions = [
            a_read_group_contaminated_at(0.02),
            a_read_group_contaminated_at(0.02),
        ];
        let frequencies = [0.9, 0.1];

        let _ = ContaminationMixture::new(
            &fractions,
            BatchOfEachReadGroup(&[BatchId(0), BatchId(1)]),
            &frequencies,
            2,
        );
    }

    #[test]
    #[should_panic(expected = "every read group runs in exactly one batch")]
    fn a_batching_that_misses_a_read_group_is_a_caller_bug() {
        let fractions = [
            a_read_group_contaminated_at(0.02),
            a_read_group_contaminated_at(0.02),
        ];
        let frequencies = [0.9, 0.1];

        let _ = ContaminationMixture::new(
            &fractions,
            BatchOfEachReadGroup(&[BatchId::ALL_TOGETHER]),
            &frequencies,
            2,
        );
    }

    /// **A mixture whose batching names fewer batches than its table has rows scores every
    /// library against the first**, and constructs cleanly without this check — the mirror of
    /// the producer's own, and the reason a defaulted batching cannot slip past both ends
    /// (C2's review).
    #[test]
    #[should_panic(expected = "no read group ran in batch 1")]
    fn a_mixture_whose_batching_leaves_a_row_unread_is_a_caller_bug() {
        let fractions = [
            a_read_group_contaminated_at(0.02),
            a_read_group_contaminated_at(0.02),
        ];
        let frequencies = [0.9, 0.1, 0.5, 0.5];

        let _ = ContaminationMixture::new(
            &fractions,
            BatchOfEachReadGroup(&[BatchId::ALL_TOGETHER, BatchId::ALL_TOGETHER]),
            &frequencies,
            2,
        );
    }

    /// **The frozen half of the mixture can be built on its own**, because the loop builds it
    /// once and hands it a refilled frequency table at every sample of every pass — and a run
    /// of no batches is refused there rather than at whichever row first indexed the table.
    #[test]
    #[should_panic(expected = "at least one sequencing batch")]
    fn a_frozen_contamination_over_no_batches_is_a_caller_bug() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let _ = FrozenContamination::new(
            &fractions,
            BatchOfEachReadGroup(&[BatchId::ALL_TOGETHER]),
            0,
        );
    }

    /// And a locus of no alleles, refused where the frequencies join it.
    #[test]
    #[should_panic(expected = "at least its reference allele")]
    fn a_frozen_contamination_over_a_locus_of_no_alleles_is_a_caller_bug() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let _ = FrozenContamination::new(
            &fractions,
            BatchOfEachReadGroup(&[BatchId::ALL_TOGETHER]),
            1,
        )
        .with_frequencies(&[], 0);
    }

    /// **A frequency table that is not this run's batches by this locus's alleles is refused
    /// where the two meet.**
    ///
    /// `ContaminationMixture::new` derives the batch count *from* the table's length, so it can
    /// only catch a length that is not a whole number of batches. The two-step door the loop
    /// uses knows the batch count already — it comes from the run's batching — so it catches the
    /// case the one-step door cannot: a table of the right shape for a different run.
    #[test]
    #[should_panic(expected = "the frequency table holds 4 entries and the run declares 1")]
    fn a_frequency_table_for_another_runs_batches_is_a_caller_bug() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let frozen = FrozenContamination::new(
            &fractions,
            BatchOfEachReadGroup(&[BatchId::ALL_TOGETHER]),
            1,
        );
        let _ = frozen.with_frequencies(&[0.9, 0.1, 0.5, 0.5], 2);
    }

    /// **The same frozen half, folded twice over two different frequency tables** — which is
    /// what the loop does at every pass, and what the split exists for: the checks between the
    /// fractions and the batching are made once and not once a sample.
    #[test]
    fn one_frozen_half_serves_two_frequency_tables() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let frozen = FrozenContamination::new(
            &fractions,
            BatchOfEachReadGroup(&[BatchId::ALL_TOGETHER]),
            1,
        );
        assert_eq!(frozen.batch_count(), 1);

        let first = frozen.with_frequencies(&[0.9, 0.1], 2);
        let second = frozen.with_frequencies(&[0.2, 0.8], 2);
        assert_eq!(
            first.contaminant_frequency_of(ReadGroupId(0), AlleleId(1)),
            0.1
        );
        assert_eq!(
            second.contaminant_frequency_of(ReadGroupId(0), AlleleId(1)),
            0.8
        );
        assert_eq!(first.fraction_of(ReadGroupId(0)), 0.02);
        assert_eq!(second.fraction_of(ReadGroupId(0)), 0.02);
    }

    #[test]
    #[should_panic(expected = "not a whole number of batches")]
    fn a_frequency_table_that_is_not_whole_rows_is_a_caller_bug() {
        let fractions = [a_read_group_contaminated_at(0.02)];
        let frequencies = [0.9, 0.1, 0.5];

        let _ = ContaminationMixture::new(
            &fractions,
            BatchOfEachReadGroup(&[BatchId::ALL_TOGETHER]),
            &frequencies,
            2,
        );
    }

    // ---- C1: what the row charges a wrong read ----

    /// **A poor read is charged what it was minted with, past a half and all**, which is the
    /// property a ceiling would break.
    ///
    /// A read minted at Phred 1 is an error of 0.794. Production's ceiling would return 0.5 —
    /// 0.46 nats away — and spec §2.3 forbids it: the merge hands the row *folds* of reads and
    /// `q_sum` recovers only their geometric mean, so no term may be a non-linear function of a
    /// per-read quality, and `min(x, ½)` is exactly such a function. On the row's own
    /// aggregation fixture that cap binds on every Phred-1 read taken singly and on none of the
    /// folds, whose geometric mean over the 93/1 alternation is `2 × 10⁻⁵`; at the fixture's
    /// 300-read end, 150 such reads, it would move the answer by 69 nats where the property is
    /// pinned to a relative `2 × 10⁻¹⁴`.
    ///
    /// The ceiling and the method that applied it were deleted on the owner's decision
    /// (2026-08-24) once nothing outside this module's tests read either. **This test is what
    /// is left of them**, so that reintroducing a cap fails here rather than moving genotypes.
    #[test]
    fn what_the_row_charges_a_poor_read_is_not_capped() {
        let calibration = ReadGroupCalibration::defaulted();
        let minted_at_phred_1 = -1.0 / 10.0 * std::f64::consts::LN_10;

        let charged = calibration.charged_error(minted_at_phred_1, 1);

        assert!((charged - 0.794_328_2).abs() < 1e-6, "charged {charged}");
        assert!(
            charged > 0.5,
            "a charge past a half is what a cap would have removed, and this is {charged}"
        );
        assert!(
            ((charged / 0.5_f64).ln() - 0.462_89).abs() < 1e-4,
            "production's ceiling would be worth {} nats on this read",
            (charged / 0.5_f64).ln()
        );
    }

    /// **The floor is the half of production's pair that ng adopted**, and it is not reached by any
    /// quality the read preparation admits: Phred 93 at the smallest scale the row's fixtures
    /// use is 185 times it. So it changes no answer that can occur, and exists only so that a
    /// logarithm never sees a zero (spec §8).
    #[test]
    fn the_floor_holds_and_sits_far_below_the_best_read_there_is() {
        let calibration = ReadGroupCalibration {
            scale: 0.37,
            provenance: Provenance::FittedHere,
        };
        let minted_at_phred_93 = -9.3 * std::f64::consts::LN_10;

        let best_there_is = calibration.charged_error(minted_at_phred_93, 1);
        assert!(
            (best_there_is / MIN_BASE_ERROR - 185.4).abs() < 0.1,
            "the best read there is sits {} times the floor",
            best_there_is / MIN_BASE_ERROR
        );

        // And the floor does hold, for a q_sum no read could carry.
        assert_eq!(calibration.charged_error(-1_000.0, 1), MIN_BASE_ERROR);
    }

    /// The charge scales exactly with the read group's multiplier, which is what makes the
    /// mixture's error half agree with spec §3.3's `q_sum + n·log scale` to a round trip
    /// rather than to a tolerance chosen by hand.
    #[test]
    fn the_charge_is_the_scale_times_the_geometric_mean() {
        let q_sum = -6.0;
        let reads = 2;
        let scaled = ReadGroupCalibration {
            scale: 2.5,
            provenance: Provenance::FittedHere,
        };

        let charged = scaled.charged_error(q_sum, reads);
        let by_hand = 2.5 * (-3.0_f64).exp();

        assert!((charged - by_hand).abs() <= f64::EPSILON * by_hand);
    }

    #[test]
    #[should_panic(expected = "has no average error")]
    fn a_charge_for_an_observation_with_no_reads_is_a_caller_bug() {
        let _ = ReadGroupCalibration::defaulted().charged_error(-7.0, 0);
    }

    /// **`f64::max` returns the other operand when one is `NaN`**, so without the check a
    /// `NaN` summed log error comes back as [`MIN_BASE_ERROR`] — the most confident error
    /// probability the module admits — and the row is finite, plausible and wrong.
    ///
    /// This is the one place [`ReadGroupCalibration::charged_error`] and
    /// [`ReadGroupCalibration::charged_error`] needed different guards:
    /// `f64::clamp` propagates a `NaN`, `f64::max` swallows it. Found by C1's review.
    #[test]
    #[should_panic(expected = "at or below zero")]
    fn a_summed_log_error_that_is_not_a_number_is_a_caller_bug() {
        let _ = ReadGroupCalibration::defaulted().charged_error(f64::NAN, 3);
    }

    /// A `q_sum` above zero is not a sum of logarithms of probabilities, and it yields a
    /// *positive* log-likelihood contribution — measured at `+48.90` on the row's own fixture
    /// before the check went in.
    #[test]
    #[should_panic(expected = "at or below zero")]
    fn a_positive_summed_log_error_is_a_caller_bug() {
        let _ = ReadGroupCalibration::defaulted().charged_error(3.0, 3);
    }

    /// `scale` is a public field, so `from_fitted_rate`'s own guard is bypassable by building
    /// the struct literally — which every fixture in this file does.
    #[test]
    #[should_panic(expected = "a finite positive multiplier")]
    fn a_scale_that_is_not_a_positive_number_is_a_caller_bug() {
        let broken = ReadGroupCalibration {
            scale: -1.0,
            provenance: Provenance::FittedHere,
        };

        let _ = broken.charged_error(-7.0, 3);
    }

    /// A `NaN` fraction is refused with a message that says what a fraction is, not only that
    /// a whole library of somebody else's DNA is not a sample — the three mistakes this bound
    /// catches are `NaN`, negative, and at-or-above one, and a reader hitting the first
    /// should not be told about the third.
    #[test]
    #[should_panic(expected = "not a `NaN`")]
    fn a_fraction_that_is_not_a_number_is_refused() {
        let fractions = [a_read_group_contaminated_at(f64::NAN)];
        let frequencies = [1.0];

        let _ = one_batch(&fractions, &frequencies);
    }
}
