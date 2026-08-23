//! Step 8 — how likely each genotype is **before any read is looked at**.
//!
//! At one locus, for one sample, this module produces one log-probability per candidate
//! genotype. The calling loop multiplies that by what the reads say and normalises; the
//! result is the posterior the caller emits (`doc/devel/ng/spec/calling_priors.md` §1).
//!
//! ## The belief has two sources, and both are measured before calling starts
//!
//! **How variable the population is.** In a nearly-invariant population almost every
//! sample carries two reference copies, and a single read showing something else is more
//! likely to be an error than a variant. In a diverse one it is not.
//!
//! **How homozygous this individual runs.** A selfing tomato accession is homozygous
//! nearly everywhere; an outbred human is not. That is the inbreeding coefficient, a
//! property of the sample rather than of the locus.
//!
//! Both arrive frozen from the parameter pre-pass — or, where the pre-pass could fit
//! nothing at all, the diversity is a species-range guess the run has to report as such
//! ([`ExpectedHeterozygosity::SPECIES_FALLBACK`](crate::ng::types::ExpectedHeterozygosity::SPECIES_FALLBACK)).
//! **This module fits nothing.** The one
//! thing it learns — how common each allele is *at this locus* — it learns from the other
//! samples in the cohort while the loop runs, and the sample's own contribution is
//! subtracted out so its reads cannot arrive twice (spec §6).
//!
//! ## Two functions with one contract between them
//!
//! 1. **Build a concentration** — per sample, per pass: the run's seed plus what the
//!    *other* samples showed here. One addition per allele.
//! 2. **Turn a concentration into a log-prior row** — one number per candidate genotype.
//!    Only the differences between them matter: the loop adds what the reads say and
//!    rescales the row to sum to one, so a constant shared by every entry cancels. This is
//!    the costly half — one `lgamma` per allele a genotype carries a copy of, plus one
//!    `logsumexp` per homozygous genotype (arch §3.2) — and the loop runs it once per
//!    sample per pass.
//!
//! ## Two words this folder leans on
//!
//! A **concentration** is one positive number per allele, read as *chromosomes the prior
//! behaves as though it had already seen*. Their ratio is the frequency it expects, their
//! sum is how much conviction that is (spec §1). Reading them as chromosome counts is what
//! makes the cohort term obvious: observed allele copies are added straight onto them,
//! because they are the same unit.
//!
//! The **seed** is the concentration a run starts from, before any locus is looked at: the
//! pre-pass's answer to *how variable is this population?*, written as chromosome counts.
//! Each locus's concentration is the seed plus what the other samples showed at that locus.
//! **Nothing refines the seed** — unlike the starting values the STR cohort's EM also calls
//! seeds (`src/ssr/cohort/em_init.rs`), which an iteration then moves, this one is frozen
//! for the whole run.
//!
//! ## What lands where
//!
//! - [`dirichlet_multinomial`] — the log-prior row itself. The Dirichlet-multinomial is
//!   what you get by drawing the locus's allele frequencies from a Dirichlet and then
//!   drawing this sample's few allele copies from those frequencies. Because it *averages
//!   over* the unknown frequencies rather than fixing them at one estimate —
//!   marginalizing — it is the default here (plan steps B1–B2).
//! - [`seed_generic`] — the SNP/indel starting point, read off the pre-pass's fitted frequency
//!   spectrum (plan step D). *Generic* is the crate's word for the non-STR path, as in
//!   `parameter_estimation::generic`, so it pairs with [`seed_ssr`] on the same axis.
//! - [`seed_ssr`] — the STR starting point: mass falling off geometrically from the
//!   cohort's modal repeat count, totalling what the cohort's measured repeat diversity
//!   implies (plan step E). **STR** in prose, `ssr` in module paths, as everywhere in ng.
//! - [`hardy_weinberg`] — the comparator: Hardy–Weinberg at a single estimated frequency,
//!   plugged in as though it were the truth, kept only so the change the marginalized prior
//!   makes stays measurable (plan step F). Named for the distribution, like its sibling
//!   [`dirichlet_multinomial`]; the plug-in character survives in the type it holds,
//!   `PlugInWrightPrior`.
//!
//! Build order and step contracts: `doc/devel/ng/impl_plan/calling_prior.md`. The design
//! and every why: `doc/devel/ng/spec/calling_priors.md`; the types and this seam:
//! `doc/devel/ng/arch/calling_priors.md`.

pub mod dirichlet_multinomial;
pub mod hardy_weinberg;
pub mod seed_generic;
pub mod seed_ssr;

pub use dirichlet_multinomial::MarginalizedDirichletPrior;
pub use seed_generic::{
    FittedSpectrum, VariantClass, fill_locus_concentration, project_spectrum_seed,
};

use crate::genetics::MIN_ALT_CONCENTRATION;
use crate::ng::types::InbreedingF;

/// The three types whose invariants are checked at construction, in a module of their own so
/// that **nothing else in this folder can build one without the check**.
///
/// The nesting is load-bearing rather than tidy, and it was measured. A private field is
/// visible to a module's *descendants*, and `dirichlet_multinomial`, `hardy_weinberg`,
/// `seed_generic` and `seed_ssr` are all descendants of `genotype_prior` — so with these
/// types declared directly in `genotype_prior`, a struct literal in any of those four files
/// compiles and skips the constructor entirely. Verified: a probe in
/// `dirichlet_multinomial.rs` built a `PriorRow` field by field, compiled, and ran. Those four
/// files are the ones that hold every implementation of [`GenotypePriorModel`] and every seed
/// builder, which is to say the only callers the checks exist for. One level of nesting makes
/// them siblings of this module instead of descendants, and the literal stops compiling.
mod checked {
    use super::{MIN_ALT_CONCENTRATION, SeedRegime};
    use crate::ng::types::{AlleleId, LogProb};

    /// How many chromosomes the prior behaves as though it had already seen, one strictly
    /// positive number per allele of the locus's table — in the same order as
    /// [`CandidateAlleles`](crate::ng::calling::CandidateAlleles), so entry 0 is the reference
    /// allele's.
    ///
    /// Their **ratio** is the allele frequency the prior expects, `α_a / Σα`. Their **sum** is
    /// how much conviction that is, because the variance of a Dirichlet frequency is
    /// `p(1 − p)/(Σα + 1)` — a larger sum is a tighter belief about the same expectation. So
    /// `(1, 0.005)` says *the alternative allele is expected at about one in two hundred, held
    /// with one chromosome's worth of conviction* (`doc/devel/ng/spec/calling_priors.md` §1).
    ///
    /// **Reading them as chromosome counts is what makes the cohort term obvious**: the other
    /// samples' observed allele copies are added straight onto these numbers, because they are
    /// the same unit.
    ///
    /// **It borrows; it never owns.** The builders fill a buffer the calling loop owns and this
    /// wraps it, so nothing allocates per sample per pass
    /// (`doc/devel/ng/spec/calling_priors.md` §8). A `Concentration` is therefore only as alive
    /// as the pass that built it.
    ///
    /// **The name is unhelpful and standard.** It sounds like a description of how concentrated
    /// the distribution is, which is true only of the sum; the literature and production's own
    /// code both call it this (`src/genetics.rs`), so the project keeps it.
    #[derive(Copy, Clone, PartialEq, Debug)]
    pub struct Concentration<'a>(&'a [f64]);

    impl<'a> Concentration<'a> {
        /// Wrap a filled buffer.
        ///
        /// **Empty is refused in release.** A locus always has a reference allele, so a
        /// zero-length concentration is a wiring bug and not a thin locus — and it is not merely
        /// a bad value: [`Self::allele_count`] is where [`PriorRow`] gets the allele count it
        /// measures the genotype table against, so a zero here would make three of that type's
        /// four checks degenerate into `0 == 0`. It is the checks' own premise, which is why
        /// production asserts the same thing in release (`src/genetics.rs`, `n_alleles > 0`).
        ///
        /// **The per-entry check is debug-only**, for the reason production gives for its own:
        /// a bad entry degrades a log-prior but cannot mis-shape any output. It is **tighter
        /// than production's `α > 0`** — every entry must clear [`MIN_ALT_CONCENTRATION`]
        /// (`1e-12`), so `1e-13` passes production and panics here. That is deliberate: every
        /// seed builder floors its output at that constant, so a smaller positive value means a
        /// builder was skipped rather than that the arithmetic drifted.
        #[inline]
        pub fn new(values: &'a [f64]) -> Self {
            assert!(
                !values.is_empty(),
                "a concentration needs one entry per allele, and every locus has a reference allele"
            );
            debug_assert!(
                values
                    .iter()
                    .all(|&a| a.is_finite() && a >= MIN_ALT_CONCENTRATION),
                "every concentration entry must be finite and at least MIN_ALT_CONCENTRATION \
                 ({MIN_ALT_CONCENTRATION:e}), got {values:?}"
            );
            Self(values)
        }

        /// The numbers themselves, one per allele, reference first.
        #[inline]
        pub fn get(self) -> &'a [f64] {
            self.0
        }

        /// How many alleles the locus is called over — the length of the table this is parallel
        /// to.
        #[inline]
        pub fn allele_count(self) -> usize {
            self.0.len()
        }
    }

    /// Everything one call to a [`GenotypePriorModel`] reads and writes at one locus for one
    /// sample, with every shape check already run.
    ///
    /// **The checks cannot be skipped, and that is the whole reason this type exists.**
    /// [`Self::new`] is the only way to build one and it runs them, so no implementation is
    /// reachable with mis-matched buffers. A trait cannot ask that of a method body — measured:
    /// while the checks were a helper function the trait's doc merely told implementations to
    /// call, deleting the call from this module's own stand-in left every test in this file
    /// passing. The failure that guards against is the one
    /// `doc/devel/ng/spec/calling_priors.md` §8 names as the worst available here: a short
    /// coefficient array lets the row loop **silently truncate** and corrupt every downstream
    /// genotype index without panicking, exactly as production records
    /// (`src/genetics.rs`).
    ///
    /// **Held in release, not debug**, for that same reason — a wrong genotype that does not
    /// crash costs more than these four integer comparisons. A debug test run cannot tell
    /// `assert_eq!` from `debug_assert_eq!`, so the "in release" half is pinned by running this
    /// module's tests a second time under `cargo test --release --lib ng::calling::genotype_prior`,
    /// where downgrading them fails. Keep both runs in this step's verification.
    ///
    /// **It borrows six buffers and owns none**, so constructing one per sample per pass costs
    /// no allocation — which is what lets the checks live here rather than in prose.
    pub struct PriorRow<'a> {
        concentration: Concentration<'a>,
        genotype_allele_counts: &'a [u32],
        log_multinomial_coeffs: &'a [f64],
        homozygous_allele_for: &'a [Option<AlleleId>],
        per_allele_scratch: &'a mut [f64],
        out: &'a mut [LogProb],
    }

    impl<'a> PriorRow<'a> {
        /// Check the six buffers against each other and bundle them.
        ///
        /// The two yardsticks are the concentration's allele count and `out`'s length; every
        /// other length is measured against those, and **each message names both buffers it
        /// compared**, because the caller reuses one row buffer across loci and an untrimmed one
        /// is the likeliest mis-shape here — a message that named only the array it found short
        /// would send the reader to the one that is correct.
        ///
        /// - `genotype_allele_counts` — `genotype_count × allele_count`, row-major: how many
        ///   copies of each allele each genotype carries.
        /// - `log_multinomial_coeffs` — one per genotype: `ln` of how many orderings of the
        ///   genome's copies spell it.
        /// - `homozygous_allele_for` — `Some(a)` where every copy is allele `a`. **This is the
        ///   one homozygous test in the caller**: the inbreeding mixture's second branch fires on
        ///   it and nothing else decides homozygosity, which is what gives the above-diploidy
        ///   question one place to change (spec §3.3).
        /// - `per_allele_scratch` — one `f64` per allele of working space. Its contents on entry
        ///   are ignored and on exit unspecified; it exists so an implementation can hold a
        ///   per-allele quantity without allocating. **Its length must match exactly**, not
        ///   merely suffice: an implementation that reduced over the whole slice rather than over
        ///   the first `allele_count` entries would fold stale values into the prior and be
        ///   silently wrong.
        /// - `out` — one [`LogProb`] per genotype, filled by the implementation.
        #[inline]
        pub fn new(
            concentration: Concentration<'a>,
            genotype_allele_counts: &'a [u32],
            log_multinomial_coeffs: &'a [f64],
            homozygous_allele_for: &'a [Option<AlleleId>],
            per_allele_scratch: &'a mut [f64],
            out: &'a mut [LogProb],
        ) -> Self {
            let allele_count = concentration.allele_count();
            let genotype_count = out.len();
            assert!(
                genotype_count > 0,
                "a locus always has at least the all-reference genotype, so a zero-length prior \
                 row is a wiring bug, not a thin locus"
            );
            assert_eq!(
                log_multinomial_coeffs.len(),
                genotype_count,
                "one log multinomial coefficient per genotype: `out` is sized for {genotype_count} \
                 genotypes and `log_multinomial_coeffs` holds {}",
                log_multinomial_coeffs.len()
            );
            assert_eq!(
                homozygous_allele_for.len(),
                genotype_count,
                "one homozygous lookup per genotype: `out` is sized for {genotype_count} genotypes \
                 and `homozygous_allele_for` holds {}",
                homozygous_allele_for.len()
            );
            let expected_counts = genotype_count
                .checked_mul(allele_count)
                .expect("genotype count times allele count overflows a usize");
            assert_eq!(
                genotype_allele_counts.len(),
                expected_counts,
                "the allele counts are genotype_count × allele_count, row-major: {genotype_count} \
                 genotypes from `out` × {allele_count} alleles from the concentration is \
                 {expected_counts}, and `genotype_allele_counts` holds {}",
                genotype_allele_counts.len()
            );
            assert_eq!(
                per_allele_scratch.len(),
                allele_count,
                "the scratch holds one entry per allele: the concentration covers {allele_count} \
                 alleles and `per_allele_scratch` holds {}",
                per_allele_scratch.len()
            );
            debug_assert!(
                homozygous_allele_for
                    .iter()
                    .flatten()
                    .all(|a| usize::from(a.0) < allele_count),
                "a homozygous lookup names an allele the concentration does not cover \
                 ({allele_count} alleles): {homozygous_allele_for:?}"
            );
            // The premise [`Self::ploidy`] rests on, checked rather than assumed. A genotype
            // table's rows all sum to the ploidy, but `genotype_allele_counts` arrives as a bare
            // slice with no tie to any table, and until B2 nothing read its *values*. Now the
            // inbreeding mixture adds `lgamma(Σα + m) − lgamma(Σα)` to one branch only, so a wrong
            // `m` is not a shared constant the row can carry — it re-weights the mixture. Measured
            // on a diploid biallelic table whose first row was edited to `[6, 0]`: the
            // homozygous-reference entry moved 5.89 nats, with nothing raised in either profile.
            let ploidy: u32 = genotype_allele_counts[..allele_count].iter().sum();
            // Held in release, because a zero here does not corrupt the row loudly: it makes the
            // correction `lgamma(Σα) − lgamma(Σα)`, which is zero, so the mixture silently reverts
            // to the unscaled one this step exists to fix. `Ploidy` refuses zero outright, and so
            // does this.
            assert!(
                ploidy > 0,
                "every genotype carries at least one copy of the genome, so the first genotype's \
                 counts cannot sum to zero: {:?} over {allele_count} alleles",
                &genotype_allele_counts[..allele_count]
            );
            // Debug-only, like every other check on the *values* in these buffers: it is O(rows)
            // and a disagreement mis-weights a prior rather than mis-shaping any output.
            debug_assert!(
                genotype_allele_counts
                    .chunks_exact(allele_count)
                    .all(|counts| counts.iter().sum::<u32>() == ploidy),
                "every genotype's copy counts must sum to the same ploidy ({ploidy} from the first \
                 genotype); they do not: {genotype_allele_counts:?} over {allele_count} alleles"
            );
            Self {
                concentration,
                genotype_allele_counts,
                log_multinomial_coeffs,
                homozygous_allele_for,
                per_allele_scratch,
                out,
            }
        }

        /// The chromosomes' worth of belief attached to each allele.
        #[inline]
        pub fn concentration(&self) -> Concentration<'a> {
            self.concentration
        }

        /// How many candidate genotypes the locus has — the length of the row.
        #[inline]
        pub fn genotype_count(&self) -> usize {
            self.log_multinomial_coeffs.len()
        }

        /// How many copies of the genome the sample has at this locus.
        ///
        /// Read off the first genotype's copy counts rather than carried as a field: every
        /// genotype's counts sum to the ploidy — that is what makes a genotype a genotype — so
        /// the first row is as good as a stored number.
        ///
        /// **That premise is checked, not assumed.** The counts arrive as a bare slice with no tie
        /// to any genotype table, so [`Self::new`] refuses a first row summing to zero in release
        /// and, in debug, a table whose genotypes disagree on the total. Without those, a mis-built
        /// count array would re-weight the inbreeding mixture — the correction it feeds sits on one
        /// branch only — and nothing would be raised.
        #[inline]
        pub fn ploidy(&self) -> u32 {
            self.genotype_allele_counts[..self.concentration.allele_count()]
                .iter()
                .sum()
        }

        /// Copies of each allele per genotype, `genotype_count × allele_count`, row-major.
        #[inline]
        pub fn genotype_allele_counts(&self) -> &'a [u32] {
            self.genotype_allele_counts
        }

        /// `ln` of how many orderings of the genome's copies spell each genotype.
        #[inline]
        pub fn log_multinomial_coeffs(&self) -> &'a [f64] {
            self.log_multinomial_coeffs
        }

        /// `Some(a)` where every copy is allele `a` — the caller's one homozygous test.
        #[inline]
        pub fn homozygous_allele_for(&self) -> &'a [Option<AlleleId>] {
            self.homozygous_allele_for
        }

        /// The per-allele working space and the row to fill, borrowed together because an
        /// implementation writes both in one pass and neither can be reborrowed while the other
        /// is held.
        #[inline]
        pub fn scratch_and_out(&mut self) -> (&mut [f64], &mut [LogProb]) {
            (self.per_allele_scratch, self.out)
        }
    }

    /// How many copies of each allele **the whole cohort** is expected to carry at this locus,
    /// this sample's own included — parallel to the locus's allele table, reference first.
    ///
    /// Expected, not counted: the entries are sums over the samples' current genotype posteriors,
    /// so no genotype is called to produce them, which is what lets them be used at low coverage
    /// (`doc/devel/ng/spec/calling_priors.md` §6).
    ///
    /// **It is a type of its own only so that it cannot be passed where
    /// [`SampleAlleleCopies`] belongs.** The two are the same shape and the same unit, the
    /// subtraction between them is the whole of the leave-one-out term, and swapping them at a
    /// call site would silently return the bare seed at every allele — the cohort term gone, with
    /// nothing raised. Measured on the flat-slice version this replaces: swapped, it did exactly
    /// that, and in release nothing caught it.
    #[derive(Copy, Clone, PartialEq, Debug)]
    pub struct CohortAlleleCopies<'a>(&'a [f64]);

    /// How many copies of each allele **this one sample** is expected to carry at this locus —
    /// the part of [`CohortAlleleCopies`] that came from it, and the part that has to come back
    /// off before the cohort's evidence can be used as this sample's prior.
    #[derive(Copy, Clone, PartialEq, Debug)]
    pub struct SampleAlleleCopies<'a>(&'a [f64]);

    macro_rules! allele_copies_impl {
        ($name:ident, $whose:literal) => {
            impl<'a> $name<'a> {
                /// Wrap a filled buffer.
                ///
                /// **The entries are checked in debug only**, which is where this module puts
                /// every check on a *value*. They are counts of genome copies: a negative, an
                /// infinity or a `NaN` is arithmetic that went wrong upstream rather than a
                /// low-coverage answer. `NaN` is the one that must not pass quietly — it is
                /// swallowed by the `max(0, ·)` in the leave-one-out term, which returns the
                /// other operand on a `NaN`, so the allele would silently come back carrying
                /// nothing but its seed.
                ///
                /// The loop's own
                /// [`ExpectedAlleleCopies`](crate::ng::calling::ExpectedAlleleCopies) makes the
                /// same check in **release** when it is built, once per locus rather than once
                /// per sample per pass; this is the cheaper restatement for the buffers that
                /// reach here.
                #[inline]
                pub fn new(copies: &'a [f64]) -> Self {
                    assert!(
                        !copies.is_empty(),
                        concat!(
                            "the ",
                            $whose,
                            " expected allele copies need one entry per allele, and every locus \
                             has a reference allele"
                        )
                    );
                    debug_assert!(
                        copies.iter().all(|c| c.is_finite() && *c >= 0.0),
                        concat!(
                            "the ",
                            $whose,
                            " expected allele copies are counts of genome copies: every entry \
                             must be finite and at or above zero, got {:?}"
                        ),
                        copies
                    );
                    Self(copies)
                }

                /// The copies themselves, one per allele, reference first.
                #[inline]
                pub fn get(self) -> &'a [f64] {
                    self.0
                }

                /// How many alleles these copies run parallel to.
                #[inline]
                pub fn allele_count(self) -> usize {
                    self.0.len()
                }
            }
        };
    }

    allele_copies_impl!(CohortAlleleCopies, "cohort's");
    allele_copies_impl!(SampleAlleleCopies, "sample's own");

    /// The SNP/indel seed: two numbers for the whole run, and where they came from.
    ///
    /// Two and not one per locus, because **the frequency spectrum is a property of the genome
    /// rather than of a locus** (spec §6). What varies locus by locus is what the other samples
    /// showed there, and that is added on afterwards.
    ///
    /// On a neutral panel the pair lands at `(1, θ)`. Neither number is a knob: `α_ref = 1` is
    /// the neutral `1/p` density written as a Dirichlet, and it is what holds the
    /// heterozygote-to-homozygous-alternative prior ratio near 2:1. Exactly, that ratio is
    /// `2·α_ref : (1 + α_alt_total)`, so it reaches 2:1 only as the alternative concentration
    /// goes to zero — but at a human `θ` of 1 in 1,000 it is 1.998:1 and on a panel ten times as
    /// diverse it is 1.98:1, which is why "at every realistic diversity" is close enough
    /// (spec §2.3, §4).
    ///
    /// **Private fields with a checking constructor**, like [`Concentration`] above: the two
    /// numbers are concentrations in the same unit, and the only other check between them and
    /// `lgamma` is compiled out of release builds.
    #[derive(Copy, Clone, PartialEq, Debug)]
    pub struct SpectrumSeed {
        alpha_ref: f64,
        alpha_alt_total: f64,
        regime: SeedRegime,
    }

    impl SpectrumSeed {
        /// Build a seed, refusing the two values that are not concentrations at all.
        ///
        /// **The alternative total may be exactly zero** — a fully invariant cohort is a real
        /// answer — and it is not floored here: the flooring at [`MIN_ALT_CONCENTRATION`]
        /// belongs to the per-locus expansion, which is where the total is split across however
        /// many alternative alleles a locus carries (arch §4).
        #[inline]
        pub fn new(alpha_ref: f64, alpha_alt_total: f64, regime: SeedRegime) -> Self {
            assert!(
                alpha_ref.is_finite() && alpha_ref > 0.0,
                "the reference concentration must be finite and strictly positive, got {alpha_ref}"
            );
            assert!(
                alpha_alt_total.is_finite() && alpha_alt_total >= 0.0,
                "the alternative concentration total must be finite and non-negative, got \
                 {alpha_alt_total}"
            );
            Self {
                alpha_ref,
                alpha_alt_total,
                regime,
            }
        }

        /// The reference allele's concentration.
        #[inline]
        pub fn alpha_ref(self) -> f64 {
            self.alpha_ref
        }

        /// The concentration shared out across whatever alternative alleles a locus carries.
        /// Splitting it keeps a site's total polymorphism independent of how many alleles it
        /// happens to have — a triallelic site is not twice as polymorphic as a biallelic one
        /// merely for holding a third allele (spec §4).
        #[inline]
        pub fn alpha_alt_total(self) -> f64 {
            self.alpha_alt_total
        }

        /// Which of the three kinds of information produced the pair. A seed cannot be built
        /// without one.
        #[inline]
        pub fn regime(self) -> SeedRegime {
            self.regime
        }
    }
}

pub use checked::{CohortAlleleCopies, Concentration, PriorRow, SampleAlleleCopies, SpectrumSeed};

/// Where the run's seed came from.
///
/// **It has to reach the run's output.** Two runs that used different information are
/// otherwise indistinguishable in what they emit, which is the complaint
/// `doc/devel/ng/spec/calling_priors.md` §4 makes about production's own fallback.
///
/// **Every variant is a branch on what the pre-pass had, never on how many samples there
/// are.** A single sample arrives with a fitted spectrum; a cohort of five arrives without
/// one, because the pre-pass emits the spectrum as *absent* below a panel-size floor rather
/// than as a thin estimate. Nothing downstream may test the cohort size (spec §4.1).
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum SeedRegime {
    /// Read off the pre-pass's fitted frequency spectrum by the projection of spec §4.1.
    ///
    /// `regularizer_site_weight` is how many sites' worth of pseudo-counts held the fitted
    /// spectrum at the neutral shape, and `census_sites_outweigh_regularizer` says whether
    /// the real sites won. **Report that per allele-count class, not as one panel-wide
    /// ratio**: on the panel spec §4.1 measures, the aggregate was 3,100 to 1 while the
    /// thinnest class held two sites and was outweighed only 39 to 1, and the tail is where
    /// the regularizer binds.
    ///
    /// `spectrum_match` says whether the two numbers actually reproduce what was measured or
    /// are the closest the family could reach. **A run that used a compromised starting point
    /// and one that matched must not look the same in the output**, which is the complaint this
    /// whole enum exists to answer.
    FittedSpectrum {
        regularizer_site_weight: f64,
        census_sites_outweigh_regularizer: bool,
        spectrum_match: SpectrumMatch,
    },
    /// No spectrum was emitted, so the pair is the neutral `(1, θ)` at the heterozygosity the
    /// pre-pass **did** fit. A spectrum too thin to emit and a panel with nothing to fit
    /// carry the same information about shape, which is why one variant covers both.
    NeutralShape,
    /// The same neutral `(1, θ)` pair as [`SeedRegime::NeutralShape`], but `θ` itself is a
    /// guess rather than a fit — too few sites, or no inbreeding coefficient for the sample —
    /// so the run is on
    /// [`ExpectedHeterozygosity::SPECIES_FALLBACK`](crate::ng::types::ExpectedHeterozygosity::SPECIES_FALLBACK),
    /// a species-range value taken from human data. **This is the variant that must never be
    /// silent.**
    FallbackDiversity,
}

/// Whether the two numbers the fit returned reproduce the spectrum it was given, or are only the
/// closest the two-parameter family could reach.
///
/// **The fit always returns a pair, and sometimes no pair is right.** Two ways that happens, and
/// neither is exotic:
///
/// - a panel whose alleles sit mostly at middling frequency — the shape
///   `doc/devel/ng/spec/calling_priors.md` §4.1 names as the one two parameters cannot hold;
/// - a panel at an inbreeding coefficient of exactly 1, where the model puts no weight at all on
///   an odd number of chromosomes carrying the allele, so a measured spectrum holding any
///   heterozygote cannot have come from any pair.
///
/// **Reported rather than hidden**, for the reason spec §12 test 11 gives about the repeat-tract
/// seed: what a builder must not do is return the closest it can reach as though it had met the
/// target.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum SpectrumMatch {
    /// The pair reproduces the measured spectrum, to the resolution the fit was asked for.
    Reproduced,
    /// **No pair of concentrations can produce the measured spectrum**, so what came back is the
    /// closest the family reaches. Detected by the winning pair predicting effectively nothing
    /// for a class the measurement gives real weight to.
    Unreproducible,
    /// **The best pair sits on the edge of the range the fit searches**, so a better one may lie
    /// outside it and what came back is a boundary rather than a summit. A fully invariant
    /// cohort reaches this legitimately: its answer is an alternative concentration of zero, and
    /// the search floors the ratio at `1e-9`.
    AtSearchLimit,
}

/// How far below the cohort's total a sample's own copies may sit before it stops being rounding
/// and starts being a defect.
///
/// The sample's own copies are one non-negative addend of the cohort total, so the true difference
/// cannot be negative and anything below zero is floating-point noise — **up to a point**. Past
/// this one it means the two paths that produce the counts have gone out of step; production names
/// the pair that can disagree in its own engine, a biallelic fast path against a per-row
/// accumulator (`src/var_calling/posterior_engine.rs`).
///
/// **The value is production's, and it holds across this caller's range with room to spare.**
/// Measured on the worst leave-one-out cancellation at ploidy scale: the two-path gap reaches
/// 3.3e-10 at 5,000 samples — about 3,000 times inside this threshold — and would not reach it
/// until roughly two million samples, far past the several thousand the caller commits to.
const COUNT_PATH_DESYNC_THRESHOLD: f64 = -1e-6;

/// Fill `out` with what this sample's prior at this locus is worth in chromosomes, before any of
/// its own reads are looked at: the run's starting concentration plus **what the other samples
/// showed here**.
///
/// ```text
/// α'_s(a) = seed(a) + max(0, cohort expected copies of a − this sample's own)
/// ```
///
/// ## Why the sample's own copies come off
///
/// **Not as a refinement — it is what makes the prior a prior.** The cohort's expected copies are
/// estimated from every sample including this one. Leave the sample's own contribution in and its
/// reads arrive twice: once through the read likelihood, and once through the allele frequency
/// they helped estimate. A genuinely homozygous-variant sample would push the frequency estimate
/// only to a diluted value and then be told that value made it heterozygous
/// (`doc/devel/ng/spec/calling_priors.md` §6).
///
/// **It is not what fixed the 214 sites of §2.2, and the spec is explicit that it is not.** The
/// spec's counterfactual runs both ways: with the subtraction in place *and* a reference
/// concentration of 10, the heterozygous prior comes out an order of magnitude further wrong than
/// the 22:1 that run actually met — so the subtraction does not repair it — while at a reference
/// concentration of 1 *without* the subtraction the failure disappears. The starting concentration
/// is the repair. This is here because using a sample's reads twice is wrong, which needs no
/// measurement to justify.
///
/// ## Both ends of the cohort range, with no branch on the cohort size
///
/// **At one sample the cohort term is exactly zero**, because the cohort total and the sample's
/// own copies are the same number — so `out` is the seed bit for bit, reached by arithmetic and
/// not by a test of `n`. That is the correct answer rather than a degraded one: **a single genome
/// carries no information about how common an allele is at a particular locus.** What it does
/// carry is how variable the genome is on average, which is `θ`, and which the seed already holds.
///
/// **At several thousand samples the cohort term swamps the seed** and the prior converges on the
/// panel's own frequencies. One formula covers both ends (spec §6).
///
/// ## The `max(0, ·)`, and the one thing it hides
///
/// It absorbs floating-point noise on the difference and nothing else; a materially negative
/// difference is refused instead, at [`COUNT_PATH_DESYNC_THRESHOLD`]. **It also returns the other
/// operand on a `NaN`**, which would turn a `NaN` copy count into an allele silently carrying
/// nothing but its seed — which is why [`CohortAlleleCopies`] and [`SampleAlleleCopies`] check
/// their entries when they are built rather than leaving it to this loop.
///
/// ## Shape and cost
///
/// Fills the caller's buffer and **allocates nothing**, so it costs no allocation per sample per
/// pass. Production's two spellings differ here and only one of them allocates: the STR cohort's
/// `leave_one_out_alpha` (`src/ssr/cohort/em.rs`) collects a fresh `Vec`, while the SNP engine
/// already writes a reused scratch buffer (`src/var_calling/posterior_engine.rs`), which is the
/// shape this follows — spec §8 records production lifting exactly these buffers out of its loop
/// after a profile put the allocator's own self-time at about one cycle in six.
///
/// **The three lengths are checked in release, and neither production spelling checks them at
/// all** — the STR one asserts in debug, the SNP one not at any level. Both can afford that: one
/// allocates its output, the other sizes its scratch once per locus. Here `out` is the caller's
/// and is reused across loci, so a short one would leave the previous locus's entries standing in
/// this locus's prior, which is the silent failure this module refuses everywhere.
///
/// The result is a valid [`Concentration`]: every entry is its seed entry plus a finite
/// non-negative number, so a seed at or above `MIN_ALT_CONCENTRATION` cannot produce one below it.
///
/// **The seed arrives as a [`Concentration`] and not as a bare slice, and that is the one thing
/// here a caller cannot get wrong by omission.** The only ways to obtain one are
/// [`fill_locus_concentration`] and [`Concentration::new`], so a loop that forgot to fill the
/// locus's seed and passed its scratch buffer straight in no longer compiles. It used to: a
/// buffer of zeros is the right length and its entries are legal floats, so every check on this
/// path passed and the prior's row came back `[NaN, −inf, NaN]`. `Concentration`'s own per-entry
/// check cannot close that, because it is a `debug_assert!` and release compiles it out.
pub fn fill_sample_concentration(
    seed: Concentration<'_>,
    cohort_copies: CohortAlleleCopies<'_>,
    own_copies: SampleAlleleCopies<'_>,
    out: &mut [f64],
) {
    let allele_count = seed.allele_count();
    assert_eq!(
        cohort_copies.allele_count(),
        allele_count,
        "one cohort copy count per allele: the seed covers {allele_count} alleles and the \
         cohort's copies hold {}",
        cohort_copies.allele_count()
    );
    assert_eq!(
        own_copies.allele_count(),
        allele_count,
        "one own copy count per allele: the seed covers {allele_count} alleles and the sample's \
         own copies hold {}",
        own_copies.allele_count()
    );
    assert_eq!(
        out.len(),
        allele_count,
        "one output entry per allele: the seed covers {allele_count} alleles and `out` holds {}",
        out.len()
    );

    for (allele, (((slot, &seed_a), &cohort_a), &own_a)) in out
        .iter_mut()
        .zip(seed.get())
        .zip(cohort_copies.get())
        .zip(own_copies.get())
        .enumerate()
    {
        let leave_one_out = cohort_a - own_a;
        debug_assert!(
            leave_one_out > COUNT_PATH_DESYNC_THRESHOLD,
            "the cohort's expected copies of allele {allele} came out materially below this \
             sample's own ({cohort_a} against {own_a}, difference {leave_one_out}); the sample's \
             own copies are one addend of the total, so the two count paths have gone out of step"
        );
        *slot = seed_a + leave_one_out.max(0.0);
    }
}

/// One sample's log-prior over every candidate genotype at one locus — the seam step 8's
/// two competing answers sit behind.
///
/// **The two answers, and why the losing one is kept.** The default integrates over the
/// unknown allele frequency; the comparator estimates the frequency and substitutes the
/// estimate as though it were the truth. The second undercounts homozygotes by exactly the
/// variance of that frequency, and on GIAB's three samples, each called on its own at 5×,
/// the difference between them was 11 points of genotype accuracy at true variants —
/// 83.6% against 94.6% (spec §2.2). The comparator ships only so that measurement stays
/// re-runnable.
///
/// ## The contract
///
/// **Values are log-priors up to one additive constant shared by every genotype**, not
/// normalised log-probabilities. The genotype-independent term of the exact
/// Dirichlet-multinomial is dropped because it cancels when the loop rescales the row, so an
/// entry may be positive and the row does not sum to anything in particular (spec §3.1).
///
/// **Every entry is finite. A genotype the prior rules out carries a very negative number, never
/// `−∞`** — the probability is floored at [`PROBABILITY_FLOOR`](crate::genetics::PROBABILITY_FLOOR)
/// before the logarithm, so it lands near `−691` and the row always has a maximum the loop can
/// subtract (spec §8, arch §1.1).
///
/// **This is a rule for every implementation, which is why it is here and not in one
/// implementation's tests.** Production's own mixture writes `−∞` at this point, so the default
/// implementation departs from what it ports; the comparator arriving at plan step F1 is instead
/// ported from `wright_genotype_log_priors`, which already floors. Two priors are compared behind
/// this seam to attribute a difference in genotypes to the *model*, and that only works if they
/// agree on this. The floor changes no call either way — moving a genotype off it would take read
/// evidence worth about 3,000 Phred (owner, 2026-08-22).
///
/// **Same inputs, bit-identical rows, at any thread count.** No RNG, no clock, no
/// thread-dependent iteration order.
///
/// **Nothing allocates.** Every buffer is the caller's, which is why the row and the
/// per-allele working space arrive inside [`PriorRow`] rather than coming back as return
/// values. Production lifted exactly these buffers out of its own loop for a measured
/// reason: a profile put the allocator's own self-time at about one cycle in six before the
/// lift (spec §8).
///
/// **A mis-shaped input is a caller bug, so there is no `Result` here or anywhere in this
/// module.** [`PriorRow::new`] has already refused it before an implementation is entered.
pub trait GenotypePriorModel {
    /// Fill the row with one log-prior per candidate genotype, in the table's genotype order.
    ///
    /// The flat views come from the locus's
    /// [`GenotypeTableView`](crate::ng::calling::GenotypeTableView); taking them flat rather
    /// than taking the view itself is what keeps this module free of a back-reference into
    /// its caller (arch §0, §7).
    fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, inbreeding: InbreedingF);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::GenotypeTable;
    use crate::ng::types::Ploidy;
    use crate::ng::types::{AlleleId, LogProb};
    use std::sync::Arc;

    /// A stand-in implementation, existing only so the seam is exercised by something. It
    /// copies the multinomial coefficients and ignores everything else — deliberately not a
    /// prior, so nothing here can be mistaken for a check of the real one.
    struct CoefficientOnlyPrior;

    impl GenotypePriorModel for CoefficientOnlyPrior {
        fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, _inbreeding: InbreedingF) {
            let coefficients = row.log_multinomial_coeffs();
            let (_, out) = row.scratch_and_out();
            for (slot, &coeff) in out.iter_mut().zip(coefficients) {
                *slot = LogProb(coeff);
            }
        }
    }

    fn table_for(copies: u8, alleles: usize) -> Arc<GenotypeTable> {
        GenotypeTable::build(Ploidy::try_new(copies).unwrap(), alleles)
    }

    /// A flat concentration of the right length for a table: reference at 1, the rest sharing
    /// a small total, which is the shape `seed_for_locus` will produce.
    fn concentration_for(alleles: usize) -> Vec<f64> {
        let mut values = vec![1e-3 / (alleles.max(2) - 1) as f64; alleles];
        values[0] = 1.0;
        values
    }

    /// Run the stand-in over a table and hand back the row it wrote.
    fn row_from_stand_in(copies: u8, alleles: usize) -> Vec<LogProb> {
        let table = table_for(copies, alleles);
        let view = table.view();
        let concentration = concentration_for(alleles);
        let mut per_allele_scratch = vec![0.0; alleles];
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut row = PriorRow::new(
            Concentration::new(&concentration),
            view.genotype_allele_counts(),
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &mut per_allele_scratch,
            &mut out,
        );
        CoefficientOnlyPrior.fill_genotype_log_priors(&mut row, InbreedingF::try_new(0.0).unwrap());
        out
    }

    /// **The one thing this step can get wrong that nothing else would catch until the loop
    /// is written**: the bundle names four flat views by type, and a real [`GenotypeTable`]'s
    /// accessors are passed straight in with nothing adapting between them. If
    /// `homozygous_alleles()` yielded bare ids, or the copy counts were `u16`, this stops
    /// compiling — which is the point.
    ///
    /// **Triallelic, not biallelic, and that is load-bearing.** A biallelic diploid row is
    /// `[0, ln 2, 0]`, a palindrome, so an implementation walking the coefficients backwards
    /// would produce it unchanged. The triallelic row is not its own reverse, so a reversed
    /// or permuted walk cannot pass. The two views the row's values never touch are pinned
    /// against the layout the bundle documents — row-major copy counts, reference first, and
    /// the homozygous lookup naming the allele every copy is.
    #[test]
    fn the_prior_trait_takes_the_genotype_tables_views_unadapted() {
        let ln2 = 2.0_f64.ln();
        // Triallelic diploid in VCF genotype order: 0/0, 0/1, 1/1, 0/2, 1/2, 2/2. Only the
        // three heterozygotes have two orderings, and they sit at 1, 3 and 4.
        assert_eq!(
            row_from_stand_in(2, 3),
            vec![
                LogProb(0.0),
                LogProb(ln2),
                LogProb(0.0),
                LogProb(ln2),
                LogProb(ln2),
                LogProb(0.0),
            ]
        );

        let table = table_for(2, 3);
        let view = table.view();
        assert_eq!(
            view.genotype_allele_counts(),
            &[2, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 2]
        );
        assert_eq!(
            view.homozygous_alleles(),
            &[
                Some(AlleleId(0)),
                None,
                Some(AlleleId(1)),
                None,
                None,
                Some(AlleleId(2)),
            ]
        );
    }

    /// The seam is usable through a trait object, which is what lets a run select between the
    /// two implementations without the calling loop being generic over the choice.
    #[test]
    fn the_prior_trait_is_reachable_through_a_trait_object() {
        let table = table_for(2, 2);
        let view = table.view();
        let concentration = concentration_for(2);
        let mut per_allele_scratch = [0.0; 2];
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut row = PriorRow::new(
            Concentration::new(&concentration),
            view.genotype_allele_counts(),
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &mut per_allele_scratch,
            &mut out,
        );
        let model: &dyn GenotypePriorModel = &CoefficientOnlyPrior;
        model.fill_genotype_log_priors(&mut row, InbreedingF::try_new(0.5).unwrap());

        assert_eq!(out, vec![LogProb(0.0), LogProb(2.0_f64.ln()), LogProb(0.0)]);
    }

    /// **The shapes this caller commits to, not only the one the tests are convenient at.**
    /// `CLAUDE.md` puts polyploids in scope and the homozygous lookup is the one thing spec
    /// §3.3 defers above diploidy, so the bundle has to accept every shape the table builds:
    /// a haploid sample, a monomorphic locus with nothing but the reference, a tetraploid,
    /// and a locus at the shipping candidate cap.
    #[test]
    fn the_bundle_accepts_every_ploidy_and_allele_count_the_table_builds() {
        // (copies, alleles, genotypes) — the genotype counts are C(alleles + copies − 1, copies).
        for (copies, alleles, genotypes) in [
            (1_u8, 1_usize, 1_usize),
            (1, 2, 2),
            (2, 1, 1),
            (2, 6, 21),
            (3, 3, 10),
            (4, 3, 15),
            (8, 2, 9),
        ] {
            let row = row_from_stand_in(copies, alleles);
            assert_eq!(
                row.len(),
                genotypes,
                "ploidy {copies} over {alleles} alleles should give {genotypes} genotypes"
            );
            assert!(
                row.iter().all(|p| p.get().is_finite()),
                "ploidy {copies} over {alleles} alleles left a non-finite entry: {row:?}"
            );
        }
    }

    /// A concentration reports the numbers it wraps and how many alleles they cover, and
    /// wrapping copies nothing — the loop's buffer is the only storage. The pointer
    /// comparison is what catches an accessor that copies: comparing values alone passes on
    /// one that returns a fresh allocation with the same contents.
    #[test]
    fn a_concentration_borrows_the_buffer_it_is_given() {
        let buffer = [1.0, 2e-3, 4e-3];
        let concentration = Concentration::new(&buffer);
        assert_eq!(concentration.allele_count(), 3);
        assert_eq!(concentration.get(), &buffer);
        assert_eq!(concentration.get().as_ptr(), buffer.as_ptr());
    }

    /// A locus always has a reference allele, so an empty concentration is a wiring bug and
    /// not a thin locus. Refused in **release** as well as debug, because the allele count it
    /// reports is what every shape check below measures against.
    #[test]
    #[should_panic(expected = "every locus has a reference allele")]
    fn an_empty_concentration_is_refused() {
        let _ = Concentration::new(&[]);
    }

    /// An entry at or below zero is refused in debug: `lgamma` is defined only for a positive
    /// argument, and the alternative concentrations are floored precisely so it stays finite
    /// when the fitted diversity is exactly zero — a fully invariant cohort.
    #[test]
    #[should_panic(expected = "MIN_ALT_CONCENTRATION")]
    #[cfg(debug_assertions)]
    fn a_non_positive_concentration_entry_is_refused_in_debug() {
        let _ = Concentration::new(&[1.0, 0.0]);
    }

    /// The floor is [`MIN_ALT_CONCENTRATION`] itself, not merely "greater than zero": a value
    /// below it would have escaped the flooring every seed builder applies, so it says the
    /// caller skipped a step rather than that the arithmetic drifted.
    #[test]
    #[should_panic(expected = "MIN_ALT_CONCENTRATION")]
    #[cfg(debug_assertions)]
    fn a_concentration_entry_under_the_floor_is_refused_in_debug() {
        let _ = Concentration::new(&[1.0, MIN_ALT_CONCENTRATION / 2.0]);
    }

    /// `+∞` is the one bad entry the floor comparison alone does not catch — it is greater
    /// than the floor — so the `is_finite` half of the check has its own case. It matters
    /// because `lgamma(+∞)` is `+∞`, and the mixture's `logsumexp` over a `+∞` returns `NaN`,
    /// losing a whole sample's row rather than one genotype's.
    #[test]
    #[should_panic(expected = "must be finite")]
    #[cfg(debug_assertions)]
    fn an_infinite_concentration_entry_is_refused_in_debug() {
        let _ = Concentration::new(&[1.0, f64::INFINITY]);
    }

    /// A row of no genotypes is the mis-shape the concentration's own emptiness check has an
    /// answer for and the row's did not: with `out` empty, three of the four length checks
    /// degenerate to `0 == 0` and an implementation writes nothing at all. Every locus has at
    /// least the all-reference genotype at any ploidy, so this is a wiring bug too.
    #[test]
    #[should_panic(expected = "at least the all-reference genotype")]
    fn an_empty_row_is_refused() {
        let concentration = [1.0];
        let mut per_allele_scratch = [0.0; 1];
        let mut out: [LogProb; 0] = [];
        let _ = PriorRow::new(
            Concentration::new(&concentration),
            &[],
            &[],
            &[],
            &mut per_allele_scratch,
            &mut out,
        );
    }

    /// Build the bundle from a well-formed diploid biallelic call with one slice replaced, so
    /// each case below differs from a passing call by exactly one length.
    fn bundle_with(
        coefficients: &[f64],
        homozygous: &[Option<AlleleId>],
        counts: &[u32],
        scratch: &mut [f64],
        out: &mut [LogProb],
    ) {
        let concentration = [1.0, 1e-3];
        let _ = PriorRow::new(
            Concentration::new(&concentration),
            counts,
            coefficients,
            homozygous,
            scratch,
            out,
        );
    }

    /// The four length checks pass together on a well-formed call — the control the four
    /// refusals below are measured against, on the same fixture rather than a wider one.
    #[test]
    fn the_row_length_checks_pass_on_a_well_formed_call() {
        let table = table_for(2, 2);
        let view = table.view();
        let mut scratch = [0.0; 2];
        let mut out = [LogProb(0.0); 3];
        bundle_with(
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            view.genotype_allele_counts(),
            &mut scratch,
            &mut out,
        );
    }

    /// A short coefficient array is the mis-shape production names: it would let a row loop
    /// truncate silently and corrupt every downstream genotype index.
    #[test]
    #[should_panic(expected = "one log multinomial coefficient per genotype")]
    fn a_short_coefficient_array_is_refused() {
        let table = table_for(2, 2);
        let view = table.view();
        let mut scratch = [0.0; 2];
        let mut out = [LogProb(0.0); 3];
        bundle_with(
            &view.log_multinomial_coeffs()[..2],
            view.homozygous_alleles(),
            view.genotype_allele_counts(),
            &mut scratch,
            &mut out,
        );
    }

    /// The homozygous lookup is the inbreeding mixture's only branch, so a short one would
    /// silently give the last genotypes the random-mating term alone.
    #[test]
    #[should_panic(expected = "one homozygous lookup per genotype")]
    fn a_short_homozygous_lookup_is_refused() {
        let table = table_for(2, 2);
        let view = table.view();
        let mut scratch = [0.0; 2];
        let mut out = [LogProb(0.0); 3];
        bundle_with(
            view.log_multinomial_coeffs(),
            &view.homozygous_alleles()[..2],
            view.genotype_allele_counts(),
            &mut scratch,
            &mut out,
        );
    }

    /// The copy-count table is the one two-dimensional buffer, and a short one would shift
    /// every genotype's row against the allele it names.
    #[test]
    #[should_panic(expected = "genotype_count × allele_count")]
    fn a_short_allele_count_table_is_refused() {
        let table = table_for(2, 2);
        let view = table.view();
        let mut scratch = [0.0; 2];
        let mut out = [LogProb(0.0); 3];
        bundle_with(
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &view.genotype_allele_counts()[..4],
            &mut scratch,
            &mut out,
        );
    }

    /// A scratch shorter than the allele count would be indexed past its end.
    #[test]
    #[should_panic(expected = "one entry per allele")]
    fn a_short_scratch_is_refused() {
        let table = table_for(2, 2);
        let view = table.view();
        let mut scratch = [0.0; 1];
        let mut out = [LogProb(0.0); 3];
        bundle_with(
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            view.genotype_allele_counts(),
            &mut scratch,
            &mut out,
        );
    }

    /// **An over-long scratch is refused too, and that is the direction that matters here.**
    /// The calling loop reuses one buffer across loci, so a scratch sized for the widest
    /// locus in a chunk and handed on untrimmed is the likeliest mis-shape at this seam. An
    /// implementation that reduced over the whole slice rather than over the first
    /// `allele_count` entries would fold another locus's values into this locus's prior and
    /// be silently wrong, which is exactly the failure class these checks exist for.
    #[test]
    #[should_panic(expected = "one entry per allele")]
    fn an_over_long_scratch_is_refused() {
        let table = table_for(2, 2);
        let view = table.view();
        let mut scratch = [0.0; 4];
        let mut out = [LogProb(0.0); 3];
        bundle_with(
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            view.genotype_allele_counts(),
            &mut scratch,
            &mut out,
        );
    }

    /// **The message names both buffers, not only the one found short.** `out` is the
    /// yardstick every genotype-count check is measured against, so a row buffer that is
    /// itself the wrong size makes a correct coefficient array look wrong; naming both is
    /// what stops an operator resizing the array that was right.
    #[test]
    fn a_length_message_names_the_yardstick_as_well_as_the_buffer() {
        let table = table_for(2, 2);
        let view = table.view();
        let panic_payload = std::panic::catch_unwind(|| {
            let mut scratch = [0.0; 2];
            let mut out = [LogProb(0.0); 6];
            bundle_with(
                view.log_multinomial_coeffs(),
                view.homozygous_alleles(),
                view.genotype_allele_counts(),
                &mut scratch,
                &mut out,
            );
        })
        .expect_err("a six-genotype row over a three-genotype table must panic");

        // `assert_eq!` with format arguments panics with a `String`; a bare `assert!` with a
        // literal panics with a `&'static str`. Read both, and say so if it is neither,
        // rather than defaulting to an empty message and blaming the assertion's wording.
        let message = panic_payload
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| panic_payload.downcast_ref::<&'static str>().copied())
            .unwrap_or("<panic payload was neither String nor &str>");
        assert!(
            message.contains("`out` is sized for 6 genotypes"),
            "the message should name the row buffer that is actually wrong, said: {message}"
        );
        assert!(
            message.contains("`log_multinomial_coeffs` holds 3"),
            "the message should also name what it was compared with, said: {message}"
        );
    }

    /// A seed cannot be built without saying where it came from, and the regime survives being
    /// copied around — it has to reach the run's output. The `assert_ne!` is what catches a
    /// hand-written equality that compares only the two numbers.
    #[test]
    fn a_spectrum_seed_carries_the_regime_that_produced_it() {
        let fitted_seed = SpectrumSeed::new(
            1.0,
            6e-4, // tomato1's fitted θ, spec §4.1
            SeedRegime::FittedSpectrum {
                regularizer_site_weight: 10.0,
                census_sites_outweigh_regularizer: true,
                spectrum_match: SpectrumMatch::Reproduced,
            },
        );
        assert_eq!(fitted_seed.alpha_ref(), 1.0);
        assert_eq!(fitted_seed.alpha_alt_total(), 6e-4);

        let neutral_seed = SpectrumSeed::new(1.0, 6e-4, SeedRegime::NeutralShape);
        let fallback_seed = SpectrumSeed::new(1.0, 1e-3, SeedRegime::FallbackDiversity);
        assert_ne!(
            neutral_seed, fitted_seed,
            "the regime is part of what a seed is"
        );
        assert_eq!(fallback_seed.regime(), SeedRegime::FallbackDiversity);
    }

    /// A fully invariant cohort has no alternative polymorphism at all, and that is a real
    /// answer rather than a bad one: the flooring that keeps `lgamma` finite happens when the
    /// total is split across a locus's alternative alleles, not here.
    #[test]
    fn a_seed_may_carry_no_alternative_concentration_at_all() {
        assert_eq!(
            SpectrumSeed::new(1.0, 0.0, SeedRegime::NeutralShape).alpha_alt_total(),
            0.0
        );
    }

    /// The reference concentration is what the whole hom-ref weight rests on, so a value that
    /// is not a concentration at all is refused rather than carried into `lgamma`.
    #[test]
    #[should_panic(expected = "reference concentration must be finite and strictly positive")]
    fn a_seed_with_no_reference_concentration_is_refused() {
        let _ = SpectrumSeed::new(0.0, 1e-3, SeedRegime::NeutralShape);
    }

    /// A negative alternative total would make `Σα` smaller than `α_ref` and every frequency
    /// the prior derives wrong in a direction nothing downstream checks.
    #[test]
    #[should_panic(expected = "alternative concentration total must be finite and non-negative")]
    fn a_seed_with_a_negative_alternative_total_is_refused() {
        let _ = SpectrumSeed::new(1.0, -1e-3, SeedRegime::NeutralShape);
    }

    /// **A first genotype whose copies sum to zero is refused in release, not only in debug.**
    ///
    /// It is the one mis-shaped count array that fails quietly: [`PriorRow::ploidy`] returns 0, the
    /// inbreeding mixture's correction becomes `lgamma(Σα) − lgamma(Σα)` — exactly zero — and the
    /// mixture silently reverts to the unscaled one that made an inbred sample's heterozygote up to
    /// 3,600 times too likely. `Ploidy` refuses zero outright and so does this.
    #[test]
    #[should_panic(expected = "cannot sum to zero")]
    fn a_first_genotype_carrying_no_copies_is_refused() {
        let table = table_for(2, 2);
        let view = table.view();
        let mut zeroed = view.genotype_allele_counts().to_vec();
        zeroed[0] = 0;
        zeroed[1] = 0;
        let mut scratch = [0.0; 2];
        let mut out = [LogProb(0.0); 3];
        bundle_with(
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &zeroed,
            &mut scratch,
            &mut out,
        );
    }

    /// **Genotypes that disagree on how many copies they carry are refused in debug.**
    ///
    /// Debug rather than release because it is O(genotypes) and because a disagreement mis-weights
    /// a prior rather than mis-shaping any output — the line this module draws for every check on
    /// the *values* in these buffers. Measured on a diploid biallelic table whose first row was
    /// edited to `[6, 0]`: the homozygous-reference log-prior moved 5.89 nats and nothing was
    /// raised before this check existed.
    #[test]
    #[should_panic(expected = "must sum to the same ploidy")]
    #[cfg(debug_assertions)]
    fn genotypes_that_disagree_on_the_copy_count_are_refused_in_debug() {
        let table = table_for(2, 2);
        let view = table.view();
        let mut disagreeing = view.genotype_allele_counts().to_vec();
        disagreeing[0] = 6;
        let mut scratch = [0.0; 2];
        let mut out = [LogProb(0.0); 3];
        bundle_with(
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &disagreeing,
            &mut scratch,
            &mut out,
        );
    }

    /// Wrap the two copy-count arrays and fill a concentration — the shape every test below uses.
    fn concentration_from(seed: &[f64], cohort: &[f64], own: &[f64]) -> Vec<f64> {
        let mut out = vec![f64::NAN; seed.len()];
        fill_sample_concentration(
            Concentration::new(seed),
            CohortAlleleCopies::new(cohort),
            SampleAlleleCopies::new(own),
            &mut out,
        );
        out
    }

    /// **At one sample the leave-one-out concentration *is* the seed, bit for bit** (spec §12
    /// test 8).
    ///
    /// With one sample the cohort's expected copies and the sample's own are the same number, so
    /// the difference is zero and the seed passes through untouched. **Bit equality rather than a
    /// tolerance, and no test of the cohort size anywhere** — the spec's rule is that one sample is
    /// reached by the arithmetic and not by a branch.
    ///
    /// **What this test does not do is pin the cohort term**, and saying so is the point: with the
    /// two copy arrays equal, an implementation that ignored them, halved them, or swapped them
    /// passes here too. Measured — all three do. The cohort term is pinned by
    /// [`raising_the_cohorts_evidence_never_lowers_an_alleles_weight`] and by
    /// [`the_leave_one_out_term_is_the_cohorts_evidence_less_the_samples_own`] below. What this one
    /// pins is that nothing scales or rounds the seed on the way through, at seed entries from 1 to
    /// 6,001 — one sample's own starting concentration up to a thousand diploid samples' worth.
    ///
    /// **Nor is the "no branch on cohort size" rule testable**, and that is a property of the rule
    /// rather than a gap: an implementation with an explicit `n == 1` branch returns bit-identical
    /// values, so no fixture can separate it. It is a shape requirement, held by review.
    #[test]
    fn at_one_sample_the_concentration_is_the_seed_bit_for_bit() {
        for allele_count in [1_usize, 2, 3, 8] {
            for reference in [1.0, 201.0, 6001.0] {
                for copies in [0.0, 1e-9, 0.5, 2.0, 2000.0] {
                    let seed: Vec<f64> = (0..allele_count)
                        .map(|a| if a == 0 { reference } else { 1e-3 * (a as f64) })
                        .collect();
                    // One sample: whatever it showed is the whole cohort's evidence.
                    let own: Vec<f64> = (0..allele_count)
                        .map(|a| copies / (1.0 + a as f64))
                        .collect();
                    let out = concentration_from(&seed, &own, &own);
                    for (allele, (got, want)) in out.iter().zip(&seed).enumerate() {
                        assert_eq!(
                            got.to_bits(),
                            want.to_bits(),
                            "{allele_count} alleles, reference {reference}, own copies {copies}, \
                             allele {allele}: got {got}, seed {want}"
                        );
                    }
                }
            }
        }
    }

    /// **The term added to the seed is the cohort's evidence less this sample's own**, entry by
    /// entry — the identity the whole function is, checked directly rather than through a
    /// monotonicity that several wrong implementations also satisfy.
    ///
    /// Measured on the flat-slice version this replaces: an implementation using the cohort's
    /// copies alone, or halving the difference, or subtracting the wrong way round, satisfies the
    /// monotonicity assertion below. None of them satisfies this one.
    #[test]
    fn the_leave_one_out_term_is_the_cohorts_evidence_less_the_samples_own() {
        let seed = [1.0, 1e-3, 1e-3, 0.5];
        let own = [1.4, 0.6, 0.0, 2.0];
        let cohort = [90.0, 0.6, 12.5, 2.0];
        let out = concentration_from(&seed, &cohort, &own);
        for (allele, ((got, seed_a), (cohort_a, own_a))) in out
            .iter()
            .zip(&seed)
            .zip(cohort.iter().zip(&own))
            .enumerate()
        {
            let want = seed_a + (cohort_a - own_a).max(0.0);
            assert!(
                (got - want).abs() < 1e-12,
                "allele {allele}: seed {seed_a} + ({cohort_a} − {own_a}) should be {want}, got \
                 {got}"
            );
        }
        // Spelled out at one allele so the identity is readable rather than only computed:
        // 90 cohort copies less this sample's 1.4 puts allele 0 at 1 + 88.6.
        assert!((out[0] - 89.6).abs() < 1e-12, "allele 0 was {}", out[0]);
    }

    /// **Raising the cohort's evidence for an allele never lowers that allele's weight for a
    /// sample that did not contribute the rise** (spec §12 test 9), and never moves any other
    /// allele's.
    #[test]
    fn raising_the_cohorts_evidence_never_lowers_an_alleles_weight() {
        let seed = [1.0, 1e-3, 1e-3];
        let own = [1.4, 0.6, 0.0];
        // Seeded below every value the loop can produce rather than at −∞, so the first of the
        // five rises carries an assertion too.
        let mut previous = vec![0.0_f64; 3];
        for extra in [0.0, 0.5, 2.0, 40.0, 4000.0] {
            let cohort = [own[0], own[1] + extra, own[2]];
            let out = concentration_from(&seed, &cohort, &own);
            assert!(
                out[1] >= previous[1],
                "raising the cohort's copies of allele 1 by {extra} lowered its weight from {} \
                 to {}",
                previous[1],
                out[1]
            );
            // The alleles whose cohort evidence did not move must not move either.
            assert_eq!(out[0].to_bits(), seed[0].to_bits());
            assert_eq!(out[2].to_bits(), seed[2].to_bits());
            previous = out;
        }
    }

    /// **The `max(0, ·)` absorbs float noise and returns the seed untouched**, rather than letting
    /// a concentration dip below the seed where [`Concentration::new`]'s floor would refuse it.
    ///
    /// A difference of `-1e-13` is the size of noise this guards; anything materially negative is
    /// a caller bug and is refused instead — the two sides are bracketed by the pair of tests
    /// below.
    #[test]
    fn float_noise_below_zero_leaves_the_seed_untouched() {
        let seed = [1.0, 1e-3];
        let own = [2.0, 0.5];
        let cohort = [2.0 - 1e-13, 0.5 - 1e-15];
        let out = concentration_from(&seed, &cohort, &own);
        assert_eq!(out[0].to_bits(), seed[0].to_bits());
        assert_eq!(out[1].to_bits(), seed[1].to_bits());
    }

    /// **The noise side of the desync threshold**: a difference of `-5e-7` is absorbed, not
    /// refused.
    ///
    /// **Written as a literal, not as a fraction of [`COUNT_PATH_DESYNC_THRESHOLD`]**, which is the
    /// whole point of the pair. A test that derives its input from the constant it means to pin
    /// moves with it and can never fail; these two brackets stay where they are, so widening the
    /// constant breaks the one below and tightening it breaks this one. Measured before they
    /// existed: the constant could be widened three orders of magnitude, to `-1e-3`, with every
    /// other test in this module still green — a drift in the direction that hides defects.
    #[test]
    fn a_difference_just_inside_the_desync_threshold_is_absorbed() {
        let seed = [1.0, 1e-3];
        let own = [2.0, 0.5];
        let cohort = [2.0 - 5e-7, 0.5];
        let out = concentration_from(&seed, &cohort, &own);
        assert_eq!(out[0].to_bits(), seed[0].to_bits());
    }

    /// **The defect side of the desync threshold**: a difference of `-2e-6` is refused. The other
    /// bracket; see the test above for why both are literals.
    #[test]
    #[should_panic(expected = "the two count paths have gone out of step")]
    #[cfg(debug_assertions)]
    fn a_difference_just_outside_the_desync_threshold_is_refused_in_debug() {
        let seed = [1.0, 1e-3];
        let own = [2.0, 0.5];
        let cohort = [2.0 - 2e-6, 0.5];
        let _ = concentration_from(&seed, &cohort, &own);
    }

    /// **A cohort total materially below the sample's own means the two count paths disagree**, and
    /// that is a caller bug rather than a small number to clamp. Debug-only, matching where
    /// production holds the same check and where this module holds every check on a *value*.
    #[test]
    #[should_panic(expected = "the two count paths have gone out of step")]
    #[cfg(debug_assertions)]
    fn a_cohort_total_below_the_samples_own_is_refused_in_debug() {
        let seed = [1.0, 1e-3];
        let own = [2.0, 0.5];
        let cohort = [1.0, 0.5];
        let _ = concentration_from(&seed, &cohort, &own);
    }

    /// **A short output buffer is refused in release.** It is the mis-shape that would otherwise
    /// leave the previous locus's entries standing in this locus's prior — the silent failure this
    /// module refuses everywhere. Neither production spelling checks this at any level; both can
    /// afford not to, one allocating its output and the other sizing its scratch once per locus.
    #[test]
    #[should_panic(expected = "one output entry per allele")]
    fn a_short_output_buffer_is_refused() {
        let seed = [1.0, 1e-3, 1e-3];
        let own = [0.0; 3];
        let cohort = [0.0; 3];
        let mut out = [f64::NAN; 2];
        fill_sample_concentration(
            Concentration::new(&seed),
            CohortAlleleCopies::new(&cohort),
            SampleAlleleCopies::new(&own),
            &mut out,
        );
    }

    /// **A short cohort count array is refused in release**, for the same reason: `zip` would stop
    /// at the shorter slice and leave the remaining alleles carrying whatever `out` held before.
    #[test]
    #[should_panic(expected = "one cohort copy count per allele")]
    fn a_short_cohort_count_array_is_refused() {
        let seed = [1.0, 1e-3, 1e-3];
        let own = [0.0; 3];
        let cohort = [0.0; 2];
        let _ = concentration_from(&seed, &cohort, &own);
    }

    /// **A short own-copies array is refused in release.** The sibling of the check above, and it
    /// was the one of the three with no test: measured, downgrading it to a `debug_assert_eq!`
    /// left this module green in both profiles while the trailing alleles silently carried the
    /// previous locus's entries.
    #[test]
    #[should_panic(expected = "one own copy count per allele")]
    fn a_short_own_count_array_is_refused() {
        let seed = [1.0, 1e-3, 1e-3];
        let own = [0.0; 2];
        let cohort = [0.0; 3];
        let _ = concentration_from(&seed, &cohort, &own);
    }

    /// **An empty seed is refused in release, and the refusal has moved into the type.** It used
    /// to be this function's own check; the seed now arrives as a [`Concentration`], which cannot
    /// be built empty, so the refusal fires one frame earlier and covers every caller rather than
    /// this one. The other three arrays keep their own wording, because a caller reusing one
    /// buffer across loci needs to know which of them is the short one.
    #[test]
    #[should_panic(expected = "a concentration needs one entry per allele")]
    fn an_empty_seed_is_refused() {
        let cohort = [0.0; 1];
        let own = [0.0; 1];
        let mut out = [f64::NAN; 1];
        fill_sample_concentration(
            Concentration::new(&[]),
            CohortAlleleCopies::new(&cohort),
            SampleAlleleCopies::new(&own),
            &mut out,
        );
    }

    /// **A `NaN` copy count is refused when the array is wrapped, in debug.** Left to the
    /// arithmetic it would pass silently: `f64::max` returns the other operand on a `NaN`, so the
    /// difference would become `0.0` and the allele would come back carrying nothing but its seed —
    /// a plausible-looking number with the cohort's evidence quietly gone.
    #[test]
    #[should_panic(expected = "must be finite and at or above zero")]
    #[cfg(debug_assertions)]
    fn a_nan_cohort_copy_count_is_refused_in_debug() {
        let _ = CohortAlleleCopies::new(&[1.0, f64::NAN]);
    }

    /// **An infinite own-copy count is refused when the array is wrapped, in debug.** It passes the
    /// desync guard — an infinite difference is not negative — and would leave an infinite
    /// concentration that [`Concentration::new`] then refuses one step later, naming the wrong
    /// buffer.
    #[test]
    #[should_panic(expected = "must be finite and at or above zero")]
    #[cfg(debug_assertions)]
    fn an_infinite_own_copy_count_is_refused_in_debug() {
        let _ = SampleAlleleCopies::new(&[1.0, f64::INFINITY]);
    }

    /// **The result is always a [`Concentration`]**, which is what lets the loop wrap it without
    /// re-checking: every entry is its seed entry plus a finite non-negative number, so a seed at
    /// or above the alternative floor cannot produce an entry below it.
    ///
    /// **The fixture gives one allele a noise-negative difference on purpose.** Without it the
    /// `max(0, ·)` is inert here — measured, deleting it left this test green — so the entry that
    /// most needs the floor was the one the test could not see.
    #[test]
    fn the_result_is_accepted_as_a_concentration() {
        let seed = [1.0, MIN_ALT_CONCENTRATION, 1e-3];
        let own = [1.0, 1e-9, 0.5];
        // Allele 1's difference is negative by float noise, which is exactly where a missing
        // `max(0, ·)` would push an entry under the floor.
        let cohort = [40.0, 1e-9 - 1e-16, 0.5];
        let out = concentration_from(&seed, &cohort, &own);
        let concentration = Concentration::new(&out);
        assert_eq!(concentration.allele_count(), 3);
        for (allele, (entry, floor)) in out.iter().zip(&seed).enumerate() {
            assert!(
                entry >= floor,
                "allele {allele}: {entry} fell below its seed {floor}"
            );
        }
        assert_eq!(out[1].to_bits(), MIN_ALT_CONCENTRATION.to_bits());
    }
}
