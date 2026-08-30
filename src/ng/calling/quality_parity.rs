//! **ng's site quality against the shipping caller's, at the same prior — and what moves when
//! the prior becomes the run's own.**
//!
//! [`score_uncorrected_site_quality`](super::quality::score_uncorrected_site_quality) is a port
//! of `src/var_calling/posterior_engine.rs`'s `compute_qual_via_exact_af` — the same collapse
//! to non-reference copies, the same linear-domain fold over the cohort's allele count, the
//! same Beta-Binomial prior on that count, the same normalisation — and what a port has to
//! prove is not that it is self-consistent but that it agrees with what it was ported from
//! (`doc/devel/ng/spec/calling_quality.md` §11 and §14's test 2).
//!
//! **It shipped without this oracle.** The module's own tests check the single-sample closed
//! form, that a locus nobody carries does not inflate as samples are added, and that the exact
//! log-domain zero term keeps a confident cohort off the ceiling. Every one of those is a
//! property of ng alone; none of them would notice the whole calculation sitting a few Phred
//! away from production's.
//!
//! # ⚠ The prior this quality reads is **not** the one the frequency loop reads
//!
//! Production carries **two** concentration pairs per record, and they are different numbers
//! from different sources:
//!
//! - `scratch.alpha`, `[1, θ̂/k, …]`, derived from the run's nucleotide diversity — the
//!   Dirichlet the **EM** uses, and the one [`loop_parity`](super::loop_parity) holds identical;
//! - `scratch.pseudocounts`, `[10, 0.01, …]`, four compiled-in constants inherited from GATK —
//!   the Beta-Binomial the **site quality** uses, and the pair spec §5.4 says ng replaces.
//!
//! `PosteriorEngineConfig::with_nucleotide_diversity` moves the first and **not** the second, so
//! a fixture that varied θ̂ and expected production's quality to follow would be testing
//! nothing. Arm one below therefore seeds ng from
//! [`DEFAULT_REF_PSEUDOCOUNT`] and [`DEFAULT_SNP_ALT_PSEUDOCOUNT`] rather than from a diversity,
//! and [`tests::both_sides_read_the_same_concentrations_and_not_two_that_happen_to_agree`]
//! moves those constants to prove that is what both sides are reading.
//!
//! # The two arms, and only one of them is parity
//!
//! **Arm one — production's prior.** Both sides are handed the same genotype log-likelihood
//! table, in the genotype order [`genotype_table_parity`](super::genotype_table_parity) already
//! pins value for value, and the same two concentrations. The numbers must agree. This is the
//! arm that says the four steps were ported correctly.
//!
//! **Arm two — ng's fitted prior.** The same table, with ng seeded from a run's own fitted
//! spectrum instead of the two GATK constants (spec §5.4). **A silent agreement here is a
//! failure rather than a pass**: it would mean the seed never reached the prior. So
//! [`tests::the_fitted_seed_moves_the_quality_and_moves_it_both_ways`] asserts a signed
//! difference in each of the two directions §5.4's table gives, and not a bound.
//!
//! Two arms that differ by a recorded decision must not be compared as though one were a bug in
//! the other — the *differential rather than parity* shape
//! [`loop_parity`](super::loop_parity) uses, and for the same reason.
//!
//! # The one place the two legitimately part
//!
//! **The ceiling is ng's alone, and it is asserted rather than avoided.** Production returns an
//! unbounded `f64` and leaves the capping to its VCF writer; ng caps inside the function, at
//! [`MAX_SITE_QUALITY`](super::quality::MAX_SITE_QUALITY), because
//! [`Phred`](crate::ng::types::Phred) refuses an infinity (spec §5.3). So the parity fixtures
//! are built below that ceiling, and
//! [`tests::a_cohort_past_the_ceiling_is_where_the_two_part_and_the_difference_is_the_cap`]
//! pins the one place they diverge.
//!
//! Nothing here changes production. Every item it names is already visible: `run_em_columnar`,
//! `EmInputs`, `MergedAllelesView`, `RecordScratch` and `RecordScratch::empty` are `pub(crate)`,
//! and `PosteriorEngineConfig`, `RecordLocus`, `AlleleSupportStats`, `MergedAllele`, the two
//! pseudocount constants and `mod backends` are `pub` — the same arrangement
//! [`loop_parity`](super::loop_parity) found.

use crate::ng::calling::genotype_prior::{SeedRegime, SpectrumSeed};
use crate::ng::calling::quality::ArtifactTestCounts;
use crate::ng::calling::quality::artifact_correction::correct_site_quality;
use crate::ng::calling::quality::score_uncorrected_site_quality;
use crate::ng::calling::{CallingScratch, CandidateAlleles, GenotypeTable};
use crate::ng::locus_generation::LocusKind;
use crate::ng::types::{AlleleId, LogProb, Phred, Ploidy};
use crate::pileup_record::AlleleSupportStats;
use crate::var_calling::per_group_merger::MergedAllele;
use crate::var_calling::posterior_engine::backends::InterpUnivariateSimdMath;
use crate::var_calling::posterior_engine::{
    DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT, EmInputs, MergedAllelesView,
    PosteriorEngineConfig, RecordLocus, RecordScratch, run_em_columnar,
};
use crate::var_calling::posterior_engine::{EmDiagnostics, PosteriorRecord};
use crate::vcf::qual_refine::refine_qual;

/// The four bases, so that each allele of a fixture is a different one byte long.
const BASES: &[u8] = b"ACGT";

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("a diploid")
}

/// **One locus, as both site qualities are handed it.** The table is `samples × genotypes`,
/// row-major — production's own layout, and ng's.
struct Locus<'a> {
    log_likelihoods: &'a [f64],
    samples: usize,
    alleles: usize,
    /// The Beta-Binomial's reference-side concentration, which production reads out of its
    /// `ref_pseudocount` and ng out of its seed.
    reference_concentration: f64,
    /// What **each** alternative allele contributes to the other side. Production sums these
    /// across the alternatives; ng's seed carries the total, so the two are tied together in
    /// [`Self::ng_at_productions_prior`] rather than typed twice.
    per_alternative_concentration: f64,
}

impl<'a> Locus<'a> {
    /// A locus at the two concentrations production ships, which is what arm one compares at.
    fn at_productions_prior(log_likelihoods: &'a [f64], samples: usize, alleles: usize) -> Self {
        Self {
            log_likelihoods,
            samples,
            alleles,
            reference_concentration: DEFAULT_REF_PSEUDOCOUNT,
            per_alternative_concentration: DEFAULT_SNP_ALT_PSEUDOCOUNT,
        }
    }

    fn genotypes(&self) -> usize {
        self.log_likelihoods.len() / self.samples
    }

    /// **The site quality production computes from this table**, in Phred.
    ///
    /// Driven through `run_em_columnar` rather than through `compute_qual_via_exact_af`, which
    /// is private — and nothing is lost by it, because the quality is a function of the
    /// *input* likelihood table and the record-static pseudocounts, neither of which the EM
    /// moves. The loop runs, and its converged frequencies do not reach this number.
    fn production_site_quality(&self) -> f64 {
        // **Every allele one base long, so nothing production does with allele *lengths* comes
        // into it** — every alternative classifies as a SNP and so draws
        // `snp_alt_pseudocount`, which is what makes the alternative side of the Beta a plain
        // multiple of one number.
        let merged: Vec<MergedAllele> = (0..self.alleles)
            .map(|allele| MergedAllele {
                seq: vec![BASES[allele % BASES.len()]],
                is_compound: false,
                constituents: Vec::new(),
            })
            .collect();
        let view = MergedAllelesView::new(&merged);
        // The per-allele read summaries production carries beside the likelihoods. The site
        // quality reads none of them — it is a fold over the table and the prior — so one
        // plausible shape repeated is enough to make a well-formed record.
        let scalars = vec![
            AlleleSupportStats {
                num_obs: 10,
                q_sum: -30.0,
                fwd: 5,
                placed_left: 5,
                placed_start: 0,
                mapq_sum: 600,
                mapq_sum_sq: 36_000,
            };
            self.samples * self.alleles
        ];
        let anchor_flags = vec![false; self.samples * self.alleles];

        let config = PosteriorEngineConfig::new()
            .with_ref_pseudocount(self.reference_concentration)
            .expect("a reference concentration in range")
            .with_snp_alt_pseudocount(self.per_alternative_concentration)
            .expect("an alternative concentration in range");
        let math = InterpUnivariateSimdMath;
        let mut scratch = RecordScratch::empty();
        run_em_columnar(
            EmInputs {
                locus: RecordLocus {
                    chrom_id: 1,
                    start: 1_000,
                    end: 1_000,
                },
                ploidy: 2,
                n_samples: self.samples,
                n_genotypes: self.genotypes(),
                alleles: &view,
                scalars: &scalars,
                log_likelihoods: self.log_likelihoods,
                chain_anchor_flags_for_validation: &anchor_flags,
            },
            &config,
            &math,
            &mut scratch,
        )
        .expect("production scores this record")
        .qual_phred
    }

    /// The candidate table the fixture's allele count spells.
    fn candidates(&self) -> CandidateAlleles {
        let mut candidates = CandidateAlleles::new(Box::from(&BASES[0..1]), LocusKind::Generic);
        for allele in 1..self.alleles {
            candidates.admit(Box::from(
                &BASES[allele % BASES.len()..allele % BASES.len() + 1],
            ));
        }
        assert_eq!(
            candidates.len(),
            self.alleles,
            "a fixture of more than four alleles would repeat a base and the table would \
             collapse it, so the two sides would be scored over different allele counts"
        );
        candidates
    }

    /// **The site quality ng computes from the same table**, under whichever prior it is given.
    fn ng_site_quality(&self, seed: SpectrumSeed) -> f64 {
        let candidates = self.candidates();
        let table = GenotypeTable::build(diploid(), self.alleles);
        let genotypes = table.view();
        assert_eq!(
            genotypes.genotype_count(),
            self.genotypes(),
            "the fixture's table is one row per genotype of this shape"
        );

        let mut scratch: CallingScratch<()> = CallingScratch::default();
        scratch.prepare_for_locus(self.samples, &candidates, &genotypes);
        for sample in 0..self.samples {
            let row =
                &self.log_likelihoods[sample * self.genotypes()..(sample + 1) * self.genotypes()];
            for (slot, &value) in scratch
                .sample_genotype_likelihoods_mut(sample)
                .iter_mut()
                .zip(row)
            {
                *slot = LogProb(value);
            }
        }
        f64::from(
            score_uncorrected_site_quality(scratch.site_quality_buffers_mut(), &genotypes, seed)
                .get(),
        )
    }

    /// ng under **production's** prior.
    ///
    /// **The same construction rather than a transcription.** Production sums its per-allele
    /// pseudocounts over the alternatives to reach the Beta's `α_alt`; this multiplies the same
    /// per-allele number by the same count. A hand-typed total would agree at two alleles and
    /// part at three, and would go on passing after a change to production's constant.
    fn ng_at_productions_prior(&self) -> f64 {
        self.ng_site_quality(SpectrumSeed::new(
            self.reference_concentration,
            self.per_alternative_concentration * (self.alleles - 1) as f64,
            SeedRegime::NeutralShape,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::quality::MAX_SITE_QUALITY;

    /// **The tolerance the two sides are compared to, and why it is not zero.**
    ///
    /// ng carries its quality as an `f32` — [`Phred`](crate::ng::types::Phred) is what reaches
    /// the `QUAL` column — and production returns an `f64`. At the 1-to-25 Phred these fixtures
    /// score, one unit in the last place of an `f32` is between `1e-7` and `2e-6`, so a
    /// tolerance below that would be a test of the narrowing conversion rather than of the
    /// arithmetic.
    ///
    /// **Measured across the six fixtures of arm one:** the largest disagreement is
    /// `1.3e-5`, at the triallelic locus; the smallest is `9.1e-8`. `1e-4` clears all six by
    /// about an order of magnitude and is still three orders below the 0.1 Phred margin arm
    /// two asserts, so it cannot swallow one of those.
    const TOLERANCE: f64 = 1e-4;

    /// One diploid sample's log-likelihoods over the three genotypes of a biallelic locus, in
    /// the table's order: homozygous reference, heterozygous, homozygous alternative.
    fn row(hom_ref: f64, het: f64, hom_alt: f64) -> [f64; 3] {
        [hom_ref, het, hom_alt]
    }

    /// One carrier of the alternative among `others` samples whose reads favour the reference —
    /// the ordinary shape of a rare variant in a cohort, and the one the prior on the count
    /// does most of its work at.
    fn one_carrier_among(others: usize) -> Vec<f64> {
        let mut table = row(-8.0, 0.0, -8.0).to_vec();
        for _ in 0..others {
            table.extend_from_slice(&row(0.0, -6.0, -12.0));
        }
        table
    }

    fn assert_agrees(locus: &Locus<'_>, what: &str) {
        let production = locus.production_site_quality();
        let ng = locus.ng_at_productions_prior();
        assert!(
            (production - ng).abs() < TOLERANCE,
            "{what}: ng's site quality is {ng} where production's is {production}, a \
             difference of {} against a tolerance of {TOLERANCE}",
            (production - ng).abs()
        );
    }

    // -----------------------------------------------------------------------------------
    // Arm one — the port, at production's prior
    // -----------------------------------------------------------------------------------

    /// **One sample, which is the case the whole calculation is checkable by hand at.** The
    /// cohort allele count runs 0, 1, 2 and the fold is three terms.
    #[test]
    fn one_sample_scores_the_same_on_both_sides() {
        let table = row(-12.0, 0.0, -12.0);
        assert_agrees(
            &Locus::at_productions_prior(&table, 1, 2),
            "one heterozygous sample",
        );
    }

    /// **A locus nobody carries.** Every sample's reads favour the reference, and both sides
    /// must come back near zero — the property §5.1's marginal exists for, checked here
    /// against production rather than against ng's own arithmetic.
    #[test]
    fn a_locus_nobody_carries_scores_the_same_on_both_sides() {
        let mut table = Vec::new();
        for _ in 0..12 {
            table.extend_from_slice(&row(0.0, -9.0, -18.0));
        }
        assert_agrees(
            &Locus::at_productions_prior(&table, 12, 2),
            "twelve reference-looking samples",
        );
    }

    /// **A cohort with one carrier**, where the fold has to spread one copy over eighty
    /// chromosomes.
    #[test]
    fn one_carrier_in_a_cohort_scores_the_same_on_both_sides() {
        let table = one_carrier_among(39);
        assert_agrees(
            &Locus::at_productions_prior(&table, 40, 2),
            "one carrier among forty",
        );
    }

    /// **A triallelic locus**, where the collapse to (reference, any-non-reference) is doing
    /// something rather than being the identity: six genotypes fold into three copy counts,
    /// and production's `α_alt` is the sum over *two* alternatives rather than one.
    #[test]
    fn a_triallelic_locus_scores_the_same_on_both_sides() {
        // The six genotypes of a diploid three-allele locus, in the table's order.
        let carrier = [-9.0, 0.0, -9.0, -3.0, -9.0, -9.0];
        let reference = [0.0, -7.0, -14.0, -7.0, -14.0, -14.0];
        let mut table = carrier.to_vec();
        for _ in 0..9 {
            table.extend_from_slice(&reference);
        }
        assert_agrees(
            &Locus::at_productions_prior(&table, 10, 3),
            "one carrier among ten at a triallelic locus",
        );
    }

    /// **A cohort whose reads say nothing at all.** Every genotype is equally likely, so what
    /// is left is the prior and the fold's own combinatorics — the one fixture where a mistake
    /// in the Beta-Binomial cannot hide behind the evidence.
    #[test]
    fn a_cohort_with_no_evidence_scores_the_prior_alone_on_both_sides() {
        let table = row(0.0, 0.0, 0.0).repeat(8);
        assert_agrees(
            &Locus::at_productions_prior(&table, 8, 2),
            "eight silent samples",
        );
    }

    /// **Both sides read the same two concentrations, and this is what says so.**
    ///
    /// Every other test in arm one runs at production's shipped constants, where a port that
    /// had hard-coded `(10, 0.01)` instead of reading the pair it is given would pass all of
    /// them. Moving both constants — the reference side down by a factor of ten, the
    /// alternative side up by a hundred — makes that port disagree, and makes a port that
    /// reads only one of the two disagree as well.
    #[test]
    fn both_sides_read_the_same_concentrations_and_not_two_that_happen_to_agree() {
        let table = one_carrier_among(19);
        let moved = Locus {
            log_likelihoods: &table,
            samples: 20,
            alleles: 2,
            reference_concentration: 1.0,
            per_alternative_concentration: 1.0,
        };
        let shipped = Locus::at_productions_prior(&table, 20, 2);
        assert!(
            (moved.production_site_quality() - shipped.production_site_quality()).abs() > 1.0,
            "the two concentration pairs have to give production materially different \
             qualities, or this fixture is not testing that either side reads them"
        );
        assert_agrees(
            &moved,
            "one carrier among twenty at a moved concentration pair",
        );
    }

    // -----------------------------------------------------------------------------------
    // Arm two — the fitted seed, where the two are meant to differ
    // -----------------------------------------------------------------------------------

    /// **The fitted seed must move the quality, and §5.4's table says which way at each end.**
    ///
    /// The prior sets how much read evidence a site must produce before its quality climbs off
    /// zero. At 63 diploid samples that toll is 16 Phred under production's constants, **23**
    /// under a neutral seed at one variant per kilobase, and **13** under one at ten per
    /// kilobase (spec §5.4). A higher toll is a lower quality for the same reads, so one
    /// carrier among forty must score *below* production's pair at human diversity and *above*
    /// it at ten times that.
    ///
    /// **Measured on this fixture:** 14.63 Phred under production's `(10, 0.01)`, **6.16**
    /// under `(1, 1e-3)` and **15.09** under `(1, 1e-2)` — so 8.5 Phred down at the neutral end
    /// and 0.5 up at the diverse one. The two are not the table's 7 and 3, and are not checked
    /// against them: §5.4 tabulates the prior toll alone, at a different cohort size, where
    /// these are qualities with a carrier's reads in them.
    ///
    /// **The margin is 0.1 Phred and it is not a claim about size.** A seed that never reached
    /// the prior would move the quality by exactly nothing; 0.1 is four orders of magnitude
    /// above the `f32` noise of [`TOLERANCE`] and far below either measured shift, so it
    /// separates *moved* from *did not move* without encoding a number this fixture cannot
    /// justify.
    #[test]
    fn the_fitted_seed_moves_the_quality_and_moves_it_both_ways() {
        let table = one_carrier_among(39);
        let locus = Locus::at_productions_prior(&table, 40, 2);
        let at_productions_pair = locus.ng_at_productions_prior();

        let human_like =
            locus.ng_site_quality(SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape));
        assert!(
            human_like < at_productions_pair - 0.1,
            "a neutral seed at one variant per kilobase asks more evidence than production's \
             constants at this cohort size, so the same reads must score lower: {human_like} \
             against {at_productions_pair}"
        );

        let diverse = locus.ng_site_quality(SpectrumSeed::new(1.0, 1e-2, SeedRegime::NeutralShape));
        assert!(
            diverse > at_productions_pair + 0.1,
            "a neutral seed at ten variants per kilobase asks less, so the same reads must \
             score higher: {diverse} against {at_productions_pair}"
        );
    }

    /// **Permuting the cohort moves the last bits and must not move the quality.**
    ///
    /// The fold walks the samples in the run's order, so reordering them genuinely changes the
    /// summation order and the answer is allowed to move in its last bits — spec §14's test 5,
    /// whose whole point is that this is the half asserted **to a tolerance**, where worker
    /// count is the half asserted bitwise. Writing the second where the first belongs is how a
    /// run-order dependency ships unnoticed; writing the first where the second belongs is a
    /// test that fails on a machine with a different sample count.
    ///
    /// **Forty samples, no two alike, and the order reversed.** A cohort of repeated rows would
    /// not test this: the fold would meet the same numbers in the same sequence whichever way
    /// the rows were listed, and an order-dependent fold would pass. Each sample here leans a
    /// different amount, so reversing them really does change what is multiplied into what.
    ///
    /// **Measured: the two orders come out identical, to the last bit of the `f32` the quality
    /// is carried in.** That is not asserted, and the difference between the two statements is
    /// the point of the test. ng's answer leaves this function as a `Phred`, whose unit in the
    /// last place at the 0.39 this fixture scores is about `3e-8`; an `f64` fold that differed
    /// by less than that rounds to the same `f32` and would pass an equality assertion having
    /// proved nothing about the fold. So the assertion stays a tolerance, which is also what
    /// spec §14's test 5 asks for — bitwise is the *worker-count* half, and that half is a
    /// run-level test this plan does not own.
    #[test]
    fn the_same_cohort_in_a_different_order_scores_the_same_within_a_tolerance() {
        let ascending: Vec<f64> = (0..40)
            .flat_map(|sample| {
                // A lean that walks from "the reads say nothing" to "firmly reference", so no
                // two rows carry the same numbers.
                let lean = f64::from(sample) * 0.37;
                row(0.0, -lean, -2.0 * lean)
            })
            .collect();
        let mut descending = Vec::with_capacity(ascending.len());
        for sample in (0..40).rev() {
            descending.extend_from_slice(&ascending[sample * 3..sample * 3 + 3]);
        }

        let up = Locus::at_productions_prior(&ascending, 40, 2).ng_at_productions_prior();
        let down = Locus::at_productions_prior(&descending, 40, 2).ng_at_productions_prior();
        assert!(
            (up - down).abs() < 1e-5,
            "the same forty samples in two orders scored {up} and {down}, a difference of {}",
            (up - down).abs()
        );
    }

    // -----------------------------------------------------------------------------------
    // Where the two part, and it is one place
    // -----------------------------------------------------------------------------------

    /// **Past ng's ceiling the two part, and the difference is the cap rather than the
    /// arithmetic.** Production returns an unbounded `f64` and leaves capping to its writer;
    /// ng caps inside the function because [`Phred`](crate::ng::types::Phred) refuses an
    /// infinity (spec §5.3). Asserted so the cap is a recorded difference rather than a
    /// surprise the day someone runs a deep cohort.
    #[test]
    fn a_cohort_past_the_ceiling_is_where_the_two_part_and_the_difference_is_the_cap() {
        let mut table = Vec::new();
        for _ in 0..60 {
            table.extend_from_slice(&row(-40.0, 0.0, -40.0));
        }
        let locus = Locus::at_productions_prior(&table, 60, 2);
        let production = locus.production_site_quality();
        let ng = locus.ng_at_productions_prior();
        assert!(
            production > f64::from(MAX_SITE_QUALITY),
            "the fixture has to drive production past ng's ceiling for this test to be about \
             the cap: production returned {production}"
        );
        assert!(
            (ng - f64::from(MAX_SITE_QUALITY)).abs() < TOLERANCE,
            "ng caps at its ceiling rather than following production upward: {ng}"
        );
    }

    // -----------------------------------------------------------------------------------
    // The artifact correction, against the code it was ported from
    // -----------------------------------------------------------------------------------

    /// **One sample's reads at a biallelic locus**, as both corrections are handed them: how
    /// many reads showed each allele, how many of those were on the forward strand, how many
    /// were placed left, and which genotype the sample was called.
    struct SampleReads {
        /// `(reads, forward, placed left)` for the reference allele.
        reference: (u32, u32, u32),
        /// The same three for the alternative.
        alternative: (u32, u32, u32),
        /// Index into [`GENOTYPES`].
        genotype: usize,
    }

    /// The three diploid genotypes of a biallelic locus, as allele indices — the table
    /// production's correction reads a sample's called genotype out of. **Supplied by the
    /// fixture rather than borrowed from either side's enumeration**, so a difference between
    /// the two corrections cannot be a difference between two genotype orders; that agreement is
    /// [`genotype_table_parity`](super::genotype_table_parity)'s to pin, and it does.
    const GENOTYPES: [[u8; 2]; 3] = [[0, 0], [0, 1], [1, 1]];

    fn stats((reads, forward, placed_left): (u32, u32, u32)) -> AlleleSupportStats {
        AlleleSupportStats {
            num_obs: reads,
            // The three moments below reach the VCF's annotations and neither correction reads
            // them; one plausible shape keeps the record well-formed.
            q_sum: -3.0 * f64::from(reads),
            fwd: forward,
            placed_left,
            placed_start: 0,
            mapq_sum: reads * 60,
            mapq_sum_sq: u64::from(reads) * 3_600,
        }
    }

    /// **The corrected quality production computes** from these reads and this baseline.
    fn production_corrected_quality(samples: &[SampleReads], baseline: f64) -> f64 {
        let alleles: Vec<MergedAllele> = b"AT"
            .iter()
            .map(|&base| MergedAllele {
                seq: vec![base],
                is_compound: false,
                constituents: Vec::new(),
            })
            .collect();
        let scalars: Vec<AlleleSupportStats> = samples
            .iter()
            .flat_map(|sample| [stats(sample.reference), stats(sample.alternative)])
            .collect();
        let record = PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 1,
                start: 1_000,
                end: 1_000,
            },
            alleles,
            ploidy: 2,
            n_samples: samples.len(),
            n_genotypes: GENOTYPES.len(),
            allele_frequencies: vec![0.5, 0.5],
            compound_frequencies: vec![None, None],
            posteriors: vec![1.0 / 3.0; samples.len() * GENOTYPES.len()],
            best_genotype: samples.iter().map(|sample| sample.genotype).collect(),
            gq_phred: vec![30.0; samples.len()],
            qual_phred: baseline,
            scalars,
            other_scalars: Vec::new(),
            chain_anchor_flags: vec![false; samples.len() * 2],
            diagnostics: EmDiagnostics {
                iterations: 3,
                final_max_delta_p: 1e-5,
                converged: true,
            },
            paralog_posterior: None,
        };
        let table: Vec<Vec<u8>> = GENOTYPES.iter().map(|genotype| genotype.to_vec()).collect();
        refine_qual(&record, &table, baseline)
    }

    /// **The nine numbers ng's correction reads**, pooled from the same samples.
    ///
    /// This repeats the definition rather than calling the worker that fills it in a run
    /// ([`summarise_condition`](super::inference::summarise_condition)), because that function
    /// needs a whole locus of evidence and a finished frequency loop to reach. What it pools is
    /// pinned by its own tests; **what this fixture has to get right is only that both sides are
    /// given the same locus**, and a disagreement caused by pooling it differently shows up here
    /// as a difference between the two corrections, which is exactly what the test is watching
    /// for.
    fn ng_summary(samples: &[SampleReads]) -> ArtifactTestCounts {
        let mut counts = ArtifactTestCounts {
            primary_alternative: AlleleId(1),
            reference_reads: 0.0,
            reference_forward_reads: 0.0,
            reference_placed_left_reads: 0.0,
            alternative_reads: 0.0,
            alternative_forward_reads: 0.0,
            alternative_placed_left_reads: 0.0,
            total_reads: 0.0,
            genotype_expected_alternative_reads: 0.0,
        };
        for sample in samples {
            let (reference_reads, reference_forward, reference_placed_left) = sample.reference;
            let (alternative_reads, alternative_forward, alternative_placed_left) =
                sample.alternative;
            counts.reference_reads += f64::from(reference_reads);
            counts.reference_forward_reads += f64::from(reference_forward);
            counts.reference_placed_left_reads += f64::from(reference_placed_left);
            counts.alternative_reads += f64::from(alternative_reads);
            counts.alternative_forward_reads += f64::from(alternative_forward);
            counts.alternative_placed_left_reads += f64::from(alternative_placed_left);
            let depth = f64::from(reference_reads) + f64::from(alternative_reads);
            counts.total_reads += depth;
            let alternative_copies = GENOTYPES[sample.genotype]
                .iter()
                .filter(|&&allele| allele == 1)
                .count();
            counts.genotype_expected_alternative_reads += (alternative_copies as f64 / 2.0) * depth;
        }
        counts
    }

    /// **The tolerance the two corrections are compared to.** ng's answer crosses into a
    /// `Phred`, whose unit in the last place at the 900 Phred these fixtures start from is about
    /// `6e-5`. **Measured: the largest disagreement across the six loci below is `2.5e-5` Phred**,
    /// at one thin heterozygote among ten at the deeper of that fixture's two scales.
    const CORRECTION_TOLERANCE: f64 = 1e-3;

    fn assert_corrections_agree(samples: &[SampleReads], baseline: f64, what: &str) {
        let production = production_corrected_quality(samples, baseline);
        let (ng, _) = correct_site_quality(
            Phred::try_new(baseline as f32).expect("a baseline quality"),
            &ng_summary(samples),
        );
        let ng = f64::from(ng.get());
        assert!(
            (production - ng).abs() < CORRECTION_TOLERANCE,
            "{what}: ng corrects to {ng} where production corrects to {production}, a \
             difference of {}",
            (production - ng).abs()
        );
    }

    /// **Production's ramp endpoints must be the compiled-in ones for any of this to mean
    /// anything.** It reads `PVC_BIAS_RAMP` once into a `OnceLock`, so a value set in the
    /// environment would silently give the two sides different ramps and every comparison below
    /// would be measuring that instead. ng has no such knob — its endpoints are typed constants,
    /// which is the point of §3.5 — so the check has only one side to make.
    #[test]
    fn productions_ramp_is_not_overridden_from_the_environment() {
        assert!(
            std::env::var("PVC_BIAS_RAMP").is_err(),
            "PVC_BIAS_RAMP is set, so production's strand ramp is not the (3, 7) ng compiles \
             in and the differential below is comparing two different tests"
        );
    }

    /// **A clean heterozygote: both sides take nothing off.** Ten samples, each half reference
    /// and half alternative reads, evenly drawn on both axes, every one called heterozygous.
    #[test]
    fn a_clean_cohort_is_corrected_the_same_by_both() {
        let samples: Vec<SampleReads> = (0..10)
            .map(|_| SampleReads {
                reference: (20, 10, 10),
                alternative: (20, 10, 10),
                genotype: 1,
            })
            .collect();
        assert_corrections_agree(&samples, 900.0, "ten clean heterozygotes");
    }

    /// **An allele-balance deficit, at two depths.** One sample in ten is called heterozygous
    /// and the alternative reads are a fifth of what that implies — the artifact shape — and the
    /// same locus again at ten times the depth, where the penalty is an order of magnitude
    /// larger and any disagreement in the tail would be too.
    #[test]
    fn an_allele_balance_deficit_is_corrected_the_same_at_two_depths() {
        for scale in [1_u32, 10] {
            let mut samples: Vec<SampleReads> = (0..9)
                .map(|_| SampleReads {
                    reference: (20 * scale, 10 * scale, 10 * scale),
                    alternative: (0, 0, 0),
                    genotype: 0,
                })
                .collect();
            samples.push(SampleReads {
                reference: (16 * scale, 8 * scale, 8 * scale),
                alternative: (4 * scale, 2 * scale, 2 * scale),
                genotype: 1,
            });
            assert_corrections_agree(&samples, 900.0, "one thin heterozygote among ten");
        }
    }

    /// **A strand-piled artifact above the ramp and below it**, which is the one place the two
    /// corrections could differ by a whole penalty rather than by a last bit. Forty alternative
    /// reads all on one strand is the full charge; three is charged nothing, and a ramp
    /// transcribed with its endpoints the other way round would give that locus the full charge
    /// on one side and nothing on the other.
    #[test]
    fn a_strand_piled_artifact_is_corrected_the_same_above_and_below_the_ramp() {
        let deep = vec![SampleReads {
            reference: (60, 30, 30),
            alternative: (40, 40, 20),
            genotype: 1,
        }];
        assert_corrections_agree(&deep, 900.0, "forty alternative reads on one strand");

        let thin = vec![SampleReads {
            reference: (60, 30, 30),
            alternative: (3, 3, 2),
            genotype: 1,
        }];
        assert_corrections_agree(&thin, 900.0, "three alternative reads on one strand");

        let midway = vec![SampleReads {
            reference: (60, 30, 30),
            alternative: (5, 5, 3),
            genotype: 1,
        }];
        assert_corrections_agree(
            &midway,
            900.0,
            "five alternative reads, halfway up the ramp",
        );
    }

    /// **A homozygous-variant cohort, which the allele-balance guard skips on both sides.** Its
    /// few reference reads are sequencing error, and a correction that had transcribed the 0.9
    /// guard as a strict inequality the other way would charge this locus where the other does
    /// not.
    #[test]
    fn a_homozygous_variant_cohort_is_skipped_by_both() {
        let samples: Vec<SampleReads> = (0..10)
            .map(|_| SampleReads {
                reference: (1, 1, 0),
                alternative: (39, 20, 19),
                genotype: 2,
            })
            .collect();
        assert_corrections_agree(&samples, 900.0, "ten homozygous-variant samples");
    }

    /// **A weak baseline that both corrections floor at zero.** The penalties exceed what there
    /// is to take them from, and neither side goes negative.
    #[test]
    fn a_baseline_smaller_than_its_penalties_floors_at_zero_on_both_sides() {
        let samples = vec![SampleReads {
            reference: (400, 200, 200),
            alternative: (100, 100, 50),
            genotype: 1,
        }];
        assert_corrections_agree(
            &samples,
            20.0,
            "a weak baseline against two large penalties",
        );
        assert_eq!(production_corrected_quality(&samples, 20.0), 0.0);
    }

    /// **A locus no read reached an alternative at keeps its baseline**, which is production's
    /// early return and ng's `None` — the two spellings of the same answer.
    ///
    /// ng's correction is never called here: the worker hands the output stage no summary at
    /// all, and that stage passes the baseline through (§3.4). So what this asserts is the half
    /// that can be asserted — that production agrees there is nothing to charge — and that ng's
    /// producer returns `None` is [`summarise_condition`](super::inference::summarise_condition)'s
    /// own test.
    #[test]
    fn a_locus_with_no_alternative_reads_keeps_its_baseline_on_both_sides() {
        let samples: Vec<SampleReads> = (0..10)
            .map(|_| SampleReads {
                reference: (20, 10, 10),
                alternative: (0, 0, 0),
                genotype: 0,
            })
            .collect();
        assert_eq!(production_corrected_quality(&samples, 743.25), 743.25);
    }
}
