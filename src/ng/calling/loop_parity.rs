//! **ng's frequency loop against the shipping caller's, genotype for genotype.**
//!
//! The loop in [`inference::summarise_condition`](super::inference::summarise_condition) is a
//! port of `src/var_calling/posterior_engine.rs`'s — the same expectation-maximization over the
//! same Dirichlet-multinomial prior with the same leave-one-out sharing — and what a port has to
//! prove is not that it is self-consistent but that it agrees with what it was ported from
//! (`doc/devel/ng/spec/calling_em_loop.md` §10, §13 test 8).
//!
//! **What is held identical.** Both sides are handed the same genotype log-likelihood table, in
//! the same genotype order — an order [`genotype_table_parity`](super::genotype_table_parity)
//! already pins value for value — the same prior concentration, the same inbreeding coefficient
//! and the same allele count. Neither side is given the other's output, and **neither side
//! computes the table**: it is the fixture's, so a difference between the two emissions cannot be
//! mistaken for a difference between the two loops.
//!
//! **Three things both loops read are not tied together here, and saying so is the point.** The
//! convergence threshold and the pass cap are two *independent* constant pairs, `1e-3` and `50` on
//! each side, equal today because two edits agree rather than because anything holds them
//! together — so [`tests::the_two_loops_run_under_the_same_threshold_and_the_same_cap`] asserts
//! the four values, and a fixture that set one side's would be hiding the other's. **The
//! transcendental backend is genuinely different and is left that way**: production interpolates
//! `ln` and `exp` from a table and ng calls the library's, which is a difference of units in the
//! last place — far below the nat-scale margins every fixture here is built on, and not something
//! a parity oracle should paper over by making ng use production's approximation.
//!
//! **The concentrations are the same construction rather than two numbers that happen to
//! agree.** Production derives `α = [ALPHA_REF, θ/k, …]` from `nucleotide_diversity`
//! ([`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs), `run_em_columnar`'s
//! prologue); ng derives it from a [`SpectrumSeed`](super::genotype_prior::SpectrumSeed)'s
//! reference concentration and alternative total, floored the same way. Seeding ng at
//! `(ALPHA_REF, θ)` therefore makes the two arrays equal entry for entry, and both read the same
//! `ALPHA_REF` out of [`crate::genetics`]. **A test that set ng's seed to a hand-typed pair would
//! be comparing two transcriptions**, and a change to production's constant would leave it
//! passing.
//!
//! # Where the two are expected to differ, and both places are pinned rather than excused
//!
//! Spec §10 asks for parity **or a difference traced to a decision one of the calling documents
//! records**. Two such differences exist, this oracle found both by failing, and each is now
//! asserted in the direction its record predicts — so losing either correction breaks a test
//! rather than quietly restoring "parity".
//!
//! **The inbreeding mixture, which is the one place the port departs on purpose.** Production
//! mixes its random-mating and identical-by-descent branches on two different scales — the first
//! is a log-prior up to a shared constant, the second a true probability — so the coefficient does
//! a fraction of the work it should; ng adds the missing term and it does all of it. **The
//! decision is recorded in a calling document, which is what §10 requires of a traced
//! difference**: `doc/devel/ng/spec/calling_priors.md` §3.2, *"ng corrects it deliberately (owner,
//! 2026-08-22), and it is the one place the port departs from what it was ported from"* — and the
//! measurement behind it is in [`MarginalizedDirichletPrior`]'s own documentation. At `F = 0` the
//! branch short-circuits away on both sides and every fixture here agrees exactly, which is why
//! production's default hides it. At `F = 0.9` the two call a thin sample differently, and
//! `the_two_loops_part_at_an_inbred_panel_and_ng_is_the_corrected_one` fixes which way round.
//!
//! **The pass count, which is a difference in what is counted and not in where either stopped.**
//! Both loops begin with one E-step on the reads alone, before any prior; production counts it as
//! its first iteration, ng does not, because its `passes` counts the passes that had a prior. So
//! production's count is ng's plus one, and every fixture here asserts that rather than ignoring
//! it: a real difference in stopping point would break it by more than one.
//!
//! # One sample and a cohort are two different production code paths, and both are compared
//!
//! **An earlier draft of this module excluded one sample, and the argument was wrong.** It said
//! the two loops have both already finished there — true of the *stopping rule*, since production
//! tests `max|p̂ − p̂'|` at one sample
//! ([`posterior_engine.rs:2724`](../../../../src/var_calling/posterior_engine.rs)) and its own
//! comment explains that a single-sample E-step is `p̂`-independent and reaches its fixed point in
//! one iteration. **It is not an argument about genotypes**, and it hid something worth knowing:
//! at one sample production dispatches to a different E-step entirely — the record-static
//! `e_step`, not the `e_step_cohort_loo` every cohort fixture here exercises. Excluding one sample
//! meant the oracle never ran that function at all, on a caller whose committed range starts
//! there (`CLAUDE.md`).
//!
//! So both are compared, and what differs between them is only *which* checks apply: the pass-count
//! relation below is a cohort fact and is asserted through [`Comparison::assert_parity`], while the
//! single-sample fixture asserts genotypes and convergence through
//! [`Comparison::assert_same_calls`] and says separately what its counts were.
//!
//! # This is ng's test; production is only the yardstick
//!
//! Nothing here changes production and nothing in production names ng, exactly as
//! [`genotype_table_parity`](super::genotype_table_parity) and
//! [`crate::ng::scanner_parity`] are arranged. **It needed no edit to the frozen tree at all.**
//! Every production item named here is already visible: `run_em_columnar`, `EmInputs`,
//! `MergedAllelesView`, `RecordScratch` and `RecordScratch::empty` are `pub(crate)`, and
//! `PosteriorEngineConfig`, `RecordLocus`, `EmDiagnostics`, `AlleleSupportStats`, `MergedAllele`
//! and `mod backends` are `pub`. The one earlier parity module cost production a `mod shape`
//! declaration; this one costs nothing.

use std::collections::BTreeMap;

use crate::genetics::ALPHA_REF;
use crate::ng::calling::genotype_prior::seed_generic::{VariantClass, fill_locus_concentration};
use crate::ng::calling::genotype_prior::{MarginalizedDirichletPrior, SeedRegime, SpectrumSeed};
use crate::ng::calling::inference::RunnableCallingLoopConfig;
use crate::ng::calling::inference::summarise_condition::{
    run_frequency_loop, summarise_final_pass,
};
use crate::ng::calling::{
    CallingScratch, CandidateAlleles, FrozenParameters, GenericLocusSample, GenericSampleEvidence,
    GenotypeIdx, GenotypeTable, LocusEvidence, LocusInference, ReadGroupCalibration,
};
use crate::ng::locus_generation::LocusKind;
use crate::ng::parameter_estimation::Provenance;
use crate::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::ng::types::{ContigId, GenomeRegion, InbreedingF, LogProb, Ploidy, Position};
use crate::pileup_record::AlleleSupportStats;
use crate::var_calling::per_group_merger::MergedAllele;
use crate::var_calling::posterior_engine::backends::InterpUnivariateSimdMath;
use crate::var_calling::posterior_engine::{
    EmDiagnostics, EmInputs, MergedAllelesView, PosteriorEngineConfig, RecordLocus, RecordScratch,
    run_em_columnar,
};

/// The four bases, so that each allele of a fixture is a different one byte long.
const BASES: &[u8] = b"ACGT";

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("a diploid")
}

fn region() -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(1),
        start: Position(1_000),
        end: Position(1_000),
    }
}

/// **One locus, as both loops are handed it.** The table is `samples × genotypes`, row-major —
/// production's own layout, and ng's.
struct Locus<'a> {
    log_likelihoods: &'a [f64],
    samples: usize,
    alleles: usize,
    /// The population's diversity, which becomes both sides' alternative concentration.
    diversity: f64,
    /// The inbreeding coefficient of each sample, in the run's sample order.
    ///
    /// **A slice rather than one value**, so a fixture can give two samples different
    /// coefficients: at one shared value a loop that read every row's coefficient off row 0 gives
    /// the same answer as one that reads each row's own, and no fixture could tell them apart.
    inbreeding: &'a [f64],
}

impl Locus<'_> {
    fn genotypes(&self) -> usize {
        self.log_likelihoods.len() / self.samples
    }

    /// **The genotypes production calls**, one index per sample into the genotype table's order,
    /// with what its loop reports about how it got there.
    fn production_calls(&self) -> (Vec<usize>, EmDiagnostics) {
        // **Every allele one base long, so nothing production does with allele *lengths* comes
        // into it.** The likelihood table is the input under test; the bytes only have to make a
        // well-formed record whose reference is entry 0.
        let merged: Vec<MergedAllele> = (0..self.alleles)
            .map(|allele| MergedAllele {
                seq: vec![BASES[allele % BASES.len()]],
                is_compound: false,
                constituents: Vec::new(),
            })
            .collect();
        let view = MergedAllelesView::new(&merged);
        // **The per-allele read summaries production carries beside the likelihoods**, which its
        // loop reads none of: they reach the reported record and its quality, not the E-step or
        // the M-step. One plausible shape repeated, rather than a second encoding of evidence
        // that would then have to agree with the table.
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
            .with_nucleotide_diversity(self.diversity)
            .expect("a diversity in range")
            .with_fixation_index_default(self.inbreeding[0])
            .expect("an inbreeding coefficient in range")
            .with_fixation_index_overrides(Some(self.inbreeding.to_vec()))
            .expect("every coefficient in range");
        let math = InterpUnivariateSimdMath;
        let mut scratch = RecordScratch::empty();
        let outputs = run_em_columnar(
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
        .expect("production calls this record");
        (outputs.best_genotype, outputs.diagnostics)
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
            "a fixture of more than four alleles would repeat a base and the table would collapse \
             it, so the two sides would be called over different allele counts"
        );
        candidates
    }

    /// **The genotypes ng calls from the same table.**
    ///
    /// The loop is driven at the two functions the shipped driver drives — the frequency loop and
    /// the final pass — rather than through `call_locus`, because `call_locus` builds the
    /// likelihood table from the evidence and this fixture's whole point is that the table is
    /// supplied to both sides.
    fn ng_calls(&self) -> LocusInference {
        let candidates = self.candidates();
        let table = GenotypeTable::build(diploid(), self.alleles);
        let genotypes = table.view();
        assert_eq!(
            genotypes.genotype_count(),
            self.genotypes(),
            "the fixture's table is one row per genotype of this shape"
        );

        // **The evidence carries no observations**, and that is what makes this a test of the
        // loop: nothing here scores a read. What it still has to say is how many samples there
        // are, that none of them was set aside, and where the locus is — all three of which the
        // final pass reads.
        let covering: Vec<GenericLocusSample<'_>> = (0..self.samples)
            .map(|_| GenericLocusSample {
                evidence: GenericSampleEvidence {
                    supported: &[],
                    unmatched_q_sum: 0.0,
                    partials: &[],
                },
                genotype_must_be_missing: false,
            })
            .collect();
        let evidence = LocusEvidence::generic(region(), &covering);

        let calibration = [ReadGroupCalibration::defaulted()];
        let inbreeding_by_sample: Vec<InbreedingF> = self
            .inbreeding
            .iter()
            .map(|&f| InbreedingF::try_new(f).expect("an inbreeding coefficient in range"))
            .collect();
        assert_eq!(
            inbreeding_by_sample.len(),
            self.samples,
            "one coefficient per sample of the run"
        );
        let strata = StratumFits::over(&[], BTreeMap::new());
        let substitution = BTreeMap::new();
        let parameters = FrozenParameters::uncontaminated(
            &calibration,
            &inbreeding_by_sample,
            // **The same construction production uses**, not a transcription of its output.
            SpectrumSeed::new(ALPHA_REF, self.diversity, SeedRegime::NeutralShape),
            &strata,
            &substitution,
            diploid(),
        );

        let mut scratch: CallingScratch<()> = CallingScratch::default();
        scratch.prepare_for_locus(self.samples, &candidates, &genotypes);
        for (sample, &coefficient) in inbreeding_by_sample.iter().enumerate() {
            scratch.claim_row_for(sample, coefficient);
        }
        let _ = fill_locus_concentration(
            parameters.prior_seed(),
            VariantClass::Substitution,
            candidates.len(),
            scratch.seed_concentration_mut(),
        );
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

        let prior = MarginalizedDirichletPrior;
        let outcome = run_frequency_loop(
            &mut scratch,
            &genotypes,
            &prior,
            &RunnableCallingLoopConfig::default(),
            None,
        );
        summarise_final_pass(
            &mut scratch,
            &genotypes,
            &evidence,
            &parameters,
            &prior,
            candidates,
            outcome,
            Provenance::FittedHere,
            None,
        )
    }

    /// One sample's called genotype as allele copies, in allele order — the shape production's
    /// genotype index resolves to through the shared table.
    fn ng_copies(&self, inference: &LocusInference, sample: usize) -> Vec<u32> {
        let genotype = inference.per_sample[sample]
            .genotype()
            .unwrap_or_else(|| panic!("sample {sample} was called missing"));
        let mut copies = vec![0_u32; self.alleles];
        for allele in genotype.alleles() {
            copies[usize::from(allele.get())] += 1;
        }
        copies
    }

    /// **Both sides' calls, as allele copies**, ready to be compared without either side's index
    /// convention standing in the way.
    fn both(&self) -> Comparison {
        let table = GenotypeTable::build(diploid(), self.alleles);
        let genotypes = table.view();
        let (indices, diagnostics) = self.production_calls();
        let production = indices
            .into_iter()
            .map(|index| {
                genotypes
                    .allele_counts_of(GenotypeIdx(
                        u32::try_from(index).expect("a genotype index fits a u32"),
                    ))
                    .expect("production's winner is an index into this shape's table")
                    .to_vec()
            })
            .collect();
        let inference = self.ng_calls();
        let ng = (0..self.samples)
            .map(|sample| self.ng_copies(&inference, sample))
            .collect();
        Comparison {
            production,
            production_iterations: diagnostics.iterations,
            production_converged: diagnostics.converged,
            ng,
            ng_passes: inference.passes,
            ng_converged: inference.converged,
        }
    }

    /// **What each sample's reads alone say**, as allele copies — the argmax of its likelihood
    /// row, with no prior at all.
    ///
    /// A fixture whose calls equal this everywhere has not tested a prior: the reads decided it,
    /// and both loops would agree however either one weighted its seed.
    fn reads_alone(&self) -> Vec<Vec<u32>> {
        let table = GenotypeTable::build(diploid(), self.alleles);
        let genotypes = table.view();
        (0..self.samples)
            .map(|sample| {
                let row = &self.log_likelihoods
                    [sample * self.genotypes()..(sample + 1) * self.genotypes()];
                let best = row
                    .iter()
                    .enumerate()
                    .max_by(|left, right| left.1.total_cmp(right.1))
                    .expect("every sample has at least one genotype")
                    .0;
                genotypes
                    .allele_counts_of(GenotypeIdx(
                        u32::try_from(best).expect("a genotype index fits a u32"),
                    ))
                    .expect("an index into this shape's table")
                    .to_vec()
            })
            .collect()
    }
}

/// **What the two loops said, side by side**, with what each reported about its own run.
struct Comparison {
    production: Vec<Vec<u32>>,
    production_iterations: u32,
    production_converged: bool,
    ng: Vec<Vec<u32>>,
    ng_passes: u32,
    ng_converged: bool,
}

impl Comparison {
    /// **The genotypes agree sample for sample, and both loops say they settled** — everything
    /// parity means that is not about counting passes.
    ///
    /// **Convergence is asserted on both sides rather than only compared**: two loops that had
    /// each run out of passes would agree on that too, and their genotypes would then be two
    /// unfinished answers rather than one right one.
    fn assert_same_calls(&self) {
        assert_eq!(
            self.ng, self.production,
            "ng called {:?} where production called {:?}",
            self.ng, self.production
        );
        assert!(
            self.ng_converged && self.production_converged,
            "one of the loops ran out of passes — ng converged: {}, production converged: {}",
            self.ng_converged,
            self.production_converged
        );
    }

    /// **The calls agree and the two pass counts differ by exactly one** — the cohort check.
    ///
    /// Both loops begin with an E-step on the reads alone before any prior; production counts it
    /// as its first iteration and ng does not, because `passes` counts the passes that had a
    /// prior. Asserted on every cohort fixture rather than on the one that takes 35 passes, so
    /// that the relation is an invariant of the port rather than an observation about a locus.
    ///
    /// **Not asserted at one sample**, where production takes a different E-step and its own
    /// stopping rule; that fixture reports its counts instead of constraining them.
    fn assert_parity(&self) {
        self.assert_same_calls();
        assert_eq!(
            self.production_iterations,
            self.ng_passes + 1,
            "ng took {} passes and production {} iterations: the two stopped in different \
             places, which is more than the one pass their two definitions differ by",
            self.ng_passes,
            self.production_iterations
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Ten outbred samples, which is at least as many as any fixture here has.
    const OUTBRED: [f64; 10] = [0.0; 10];
    /// The same ten at a strongly selfing panel's coefficient — tomato's order.
    const SELFING: [f64; 10] = [0.9; 10];

    /// A two-allele diploid table's three genotypes, in the order both sides enumerate them:
    /// `0/0`, `0/1`, `1/1`.
    fn row(hom_ref: f64, het: f64, hom_alt: f64) -> [f64; 3] {
        [hom_ref, het, hom_alt]
    }

    /// **Three samples whose reads say three different things, and both loops say the same three
    /// things back.**
    ///
    /// The baseline case, and the one that proves least: **ten to twelve nats separate each
    /// sample's top two genotypes**, so the reads decide every sample and no prior either loop
    /// could hold would move one. It is
    /// here so that a wiring mistake — a table read in the wrong order, a sample's row handed to
    /// its neighbour — fails somewhere obvious, and the fixtures below are what test the loop.
    #[test]
    fn the_two_loops_call_the_same_genotypes_at_a_biallelic_site() {
        let mut table = Vec::new();
        table.extend_from_slice(&row(-30.0, -10.0, 0.0));
        table.extend_from_slice(&row(-12.0, 0.0, -12.0));
        table.extend_from_slice(&row(0.0, -10.0, -30.0));
        let locus = Locus {
            log_likelihoods: &table,
            samples: 3,
            alleles: 2,
            diversity: 1e-3,
            inbreeding: &OUTBRED[..3],
        };
        let comparison = locus.both();
        comparison.assert_parity();
        assert_eq!(
            comparison.production,
            vec![vec![0, 2], vec![1, 1], vec![2, 0]],
            "the reads' own answer, which is what a fixture this loud must give"
        );
    }

    /// **A cohort where the prior decides, and the two loops decide it the same way.**
    ///
    /// Nine samples carry the alternative allele unambiguously and one — the last — has reads
    /// that put its heterozygote **0.7 nats** above its homozygous reference, which is inside
    /// what a prior can move. The nine make the alternative common, the leave-one-out share hands
    /// that to the tenth, and both loops call it heterozygous **where its own reads alone say
    /// heterozygous too** — so this fixture is checked the other way round, against a prior that
    /// would have overturned the reads if it were the seed's rather than the cohort's.
    ///
    /// **The check that this fixture tests anything is the assertion on the first sample of the
    /// pair below**, not this one. Here what is asserted is that ten samples, a moving frequency
    /// and more than one pass produce the same answer on both sides.
    #[test]
    fn the_two_loops_agree_where_the_cohort_carries_a_thin_sample() {
        let mut table = Vec::new();
        for _ in 0..9 {
            table.extend_from_slice(&row(-25.0, -8.0, 0.0));
        }
        table.extend_from_slice(&row(-0.7, 0.0, -6.0));
        let locus = Locus {
            log_likelihoods: &table,
            samples: 10,
            alleles: 2,
            diversity: 1e-3,
            inbreeding: &OUTBRED,
        };
        let comparison = locus.both();
        comparison.assert_parity();
        assert!(
            comparison.ng_passes > 1,
            "a locus whose frequency moves takes more than one pass; this one took {}",
            comparison.ng_passes
        );
    }

    /// **The prior overturns the reads at one sample, and both loops overturn them the same
    /// way** — which is the fixture that makes this oracle worth running.
    ///
    /// Nine samples are homozygous reference beyond doubt, so the cohort's alternative frequency
    /// goes to nearly nothing. The tenth's reads put its **heterozygote 1.2 nats above** its
    /// homozygous reference — the reads' own answer is `0/1`. Both loops call it `0/0`: the
    /// leave-one-out prior from nine reference homozygotes is worth more than 1.2 nats.
    ///
    /// **So the call is not the reads' answer**, and that is asserted rather than described. A
    /// loop that ignored its prior entirely would pass every other fixture in this file and fail
    /// this one.
    #[test]
    fn both_loops_let_the_cohort_overturn_a_thin_samples_reads() {
        let mut table = Vec::new();
        for _ in 0..9 {
            table.extend_from_slice(&row(0.0, -25.0, -50.0));
        }
        table.extend_from_slice(&row(-1.2, 0.0, -20.0));
        let locus = Locus {
            log_likelihoods: &table,
            samples: 10,
            alleles: 2,
            diversity: 1e-3,
            inbreeding: &OUTBRED,
        };
        let comparison = locus.both();
        comparison.assert_parity();

        let reads_alone = locus.reads_alone();
        assert_eq!(
            reads_alone[9],
            vec![1, 1],
            "the tenth sample's own reads say heterozygous"
        );
        assert_eq!(
            comparison.ng[9],
            vec![2, 0],
            "and both loops call it homozygous reference, on the cohort's evidence rather than \
             its own"
        );
    }

    /// **The one place the port departs from what it was ported from, pinned rather than
    /// excused** — and it is invisible until the panel is inbred.
    ///
    /// The two priors mix the same two branches, and production's mixes them **on two different
    /// scales**: its random-mating branch is a log-prior up to a shared additive constant, its
    /// identical-by-descent branch a true probability, so the first is inflated by `Σα(Σα + 1)`
    /// and the inbreeding coefficient does a fraction of the work it should. ng adds the missing
    /// `lgamma(Σα + m) − lgamma(Σα)` and the coefficient does all of it (owner, 2026-08-22;
    /// [`MarginalizedDirichletPrior`]'s own documentation carries the measurement).
    ///
    /// **Where it shows and where it cannot.** At `F = 0` the branch short-circuits away on both
    /// sides, which is why every other fixture in this file agrees exactly and why production's
    /// own default hides this. At `F = 0.9` over ten samples it shows: the thin sample's reads put
    /// its heterozygote **1.0 nat** above its homozygous reference, ng's corrected coefficient is
    /// worth more than that and calls `0/0`, and production's is not and calls `0/1`.
    ///
    /// **How much of its work production's coefficient does depends on the cohort, and this
    /// fixture's ten samples are near the favourable end.** The inflation is `Σα(Σα + 1)` and `Σα`
    /// is the leave-one-out concentration, which grows with the cohort — so a *smaller* cohort
    /// leaves production's coefficient further along, and the sibling measurement puts fifty
    /// samples at 3.6% of the way and a thousand at 0.09%. Ten samples get more than 3.6% of the
    /// way and it is still not enough to move one nat, which is what makes this the fixture rather
    /// than a panel of sixty.
    ///
    /// **So this test asserts the difference**, in both directions: agreement at `F = 0`, and one
    /// named sample differing at `F = 0.9` with ng on the side its correction predicts. A parity
    /// oracle with an escape clause would have skipped the inbred case; this one fails if the
    /// departure ever stops happening, which is what would tell somebody the correction had been
    /// lost.
    #[test]
    fn the_two_loops_part_at_an_inbred_panel_and_ng_is_the_corrected_one() {
        let mut table = Vec::new();
        for _ in 0..5 {
            table.extend_from_slice(&row(-25.0, -8.0, 0.0));
        }
        for _ in 0..4 {
            table.extend_from_slice(&row(0.0, -8.0, -25.0));
        }
        table.extend_from_slice(&row(-1.0, 0.0, -8.0));

        let outbred = Locus {
            log_likelihoods: &table,
            samples: 10,
            alleles: 2,
            diversity: 1e-3,
            inbreeding: &OUTBRED,
        };
        let selfing = Locus {
            inbreeding: &SELFING,
            ..outbred
        };

        // Outbred: the branch is not taken, and the two loops agree sample for sample.
        let without = outbred.both();
        without.assert_parity();
        assert_eq!(
            without.ng[9],
            vec![1, 1],
            "outbred, the thin sample's own reads decide it"
        );

        // Inbred: they part, and only at the sample whose reads left a nat of room.
        let with = selfing.both();
        // **Both still finished, and both still counted the same passes.** Without this the
        // fixture could be showing two unfinished answers rather than one departure — the failure
        // `assert_same_calls` exists to prevent, which the genotype assertions below cannot make
        // because they are asserting a *difference*.
        assert!(
            with.ng_converged && with.production_converged,
            "one of the loops ran out of passes at F = 0.9, so its genotypes are unfinished \
             rather than corrected"
        );
        assert_eq!(
            with.production_iterations,
            with.ng_passes + 1,
            "the inbred run stopped in different places, which is a second difference this \
             fixture does not claim"
        );
        assert_eq!(
            with.ng[..9],
            with.production[..9],
            "the nine samples whose reads decide them are called alike either way"
        );
        assert_eq!(
            with.ng[9],
            vec![2, 0],
            "ng's corrected coefficient is worth more than the reads' one nat"
        );
        assert_eq!(
            with.production[9],
            vec![1, 1],
            "production's is worth a fraction of that and leaves the reads standing"
        );
        assert_ne!(
            with.ng, without.ng,
            "and the coefficient moves ng at all — without this the fixture would show a \
             departure that was really a prior doing nothing on either side"
        );
    }

    /// **Three alleles: the concentration is split rather than handed to one allele, and the
    /// split is asserted by value because the genotypes cannot see it.**
    ///
    /// Production divides `θ` among the alternatives and so does ng, and at two alleles the
    /// division is by one — so every other fixture would pass against a loop that never divided.
    /// **So would this one's genotype comparison**, and that was a defect in an earlier draft
    /// which claimed otherwise: deleting the division from `fill_locus_concentration` changes no
    /// call and no pass count here. It cannot: at a diversity of 1 in 1,000 the difference between
    /// `θ` and `θ/2` is a ten-thousandth of a chromosome against a leave-one-out term carrying
    /// whole ones, and every sample's reads are decided by six nats or more.
    ///
    /// **What can fail is the concentration itself**, so that is what is checked: `[1, θ/2, θ/2]`,
    /// entry for entry, which is production's `α` prologue written out. Beside it the two loops
    /// are still compared at three alleles, which is the shape the M-step and the table are
    /// exercised at and which no other fixture reaches.
    ///
    /// *(The stronger check — ng's concentration against production's own array rather than
    /// against its formula — is not available: `RecordScratch::alpha` is private to
    /// `posterior_engine`, and widening the frozen tree for a test is not a trade this module
    /// makes.)*
    #[test]
    fn the_two_loops_split_the_alternative_concentration_the_same_way() {
        // Six genotypes at three alleles, in table order: 0/0, 0/1, 1/1, 0/2, 1/2, 2/2.
        let table = vec![
            // four samples carrying allele 1
            -20.0, -6.0, 0.0, -20.0, -6.0, -20.0, //
            -20.0, -6.0, 0.0, -20.0, -6.0, -20.0, //
            -20.0, -6.0, 0.0, -20.0, -6.0, -20.0, //
            -20.0, -6.0, 0.0, -20.0, -6.0, -20.0, //
            // two carrying allele 2
            -20.0, -20.0, -20.0, -6.0, -6.0, 0.0, //
            -20.0, -20.0, -20.0, -6.0, -6.0, 0.0, //
            // and one whose reads put its heterozygote a nat ahead
            -1.0, 0.0, -8.0, -1.5, -9.0, -12.0,
        ];
        let diversity = 1e-3;
        let locus = Locus {
            log_likelihoods: &table,
            samples: 7,
            alleles: 3,
            diversity,
            inbreeding: &OUTBRED[..7],
        };
        locus.both().assert_parity();

        // **The split, by value.** Production's prologue writes `α = [ALPHA_REF, θ/k, …]`; this is
        // ng's, at `k = 2`. A loop that handed each alternative the whole total would read
        // `[1, θ, θ]` and fail here — which is the check the genotype comparison above cannot make.
        let candidates = locus.candidates();
        let table_view = GenotypeTable::build(diploid(), 3);
        let genotypes = table_view.view();
        let mut scratch: CallingScratch<()> = CallingScratch::default();
        scratch.prepare_for_locus(1, &candidates, &genotypes);
        let _ = fill_locus_concentration(
            SpectrumSeed::new(ALPHA_REF, diversity, SeedRegime::NeutralShape),
            VariantClass::Substitution,
            3,
            scratch.seed_concentration_mut(),
        );
        assert_eq!(
            scratch.seed_concentration(),
            &[ALPHA_REF, diversity / 2.0, diversity / 2.0],
            "the alternative total is divided among the alternatives, not given to each of them"
        );
    }

    /// **One sample, which is a different function in production and the bottom of this
    /// caller's committed range.**
    ///
    /// At `n_samples == 1` production dispatches to its record-static `e_step` rather than the
    /// `e_step_cohort_loo` every other fixture here exercises, and it stops on `max|p̂ − p̂'|`
    /// rather than on expected copies over chromosomes. ng runs the same code at every cohort
    /// size. **So this is the fixture that runs production's other E-step**, and an earlier draft
    /// of this module left it out on an argument about the stopping rule that says nothing about
    /// genotypes.
    ///
    /// **The pass relation is reported rather than asserted here**: the `+1` holds because
    /// production counts its prior-free first iteration, but at one sample its loop is also
    /// stopping on a different quantity, so tying the two counts together would be pinning a
    /// coincidence. What is asserted is the call and that both finished.
    #[test]
    fn the_two_loops_call_one_sample_alike_on_production_s_other_e_step() {
        let table = row(-1.2, 0.0, -20.0).to_vec();
        let locus = Locus {
            log_likelihoods: &table,
            samples: 1,
            alleles: 2,
            diversity: 1e-3,
            inbreeding: &OUTBRED[..1],
        };
        let comparison = locus.both();
        comparison.assert_same_calls();
        assert_eq!(
            comparison.ng[0],
            vec![2, 0],
            "one sample, no cohort to borrow from: the seed's pull toward the reference is worth \
             more than the 1.2 nats its reads put on the heterozygote"
        );
        assert_ne!(
            comparison.ng[0],
            locus.reads_alone()[0],
            "and the call is not the reads' own answer, so the prior decided it here too"
        );
    }

    /// **Two samples with identical reads and different inbreeding coefficients: the departure
    /// appears at one of them and not the other.**
    ///
    /// Every other fixture gives every sample the same coefficient, which makes *this row's
    /// coefficient* and *row 0's coefficient* the same number — so a loop that read one for the
    /// other would pass all of them. Here samples 6 and 7 have the same likelihood row and
    /// coefficients of 0 and 0.95, and production is given the same vector through the per-sample
    /// override its own pipeline uses when the diversity estimator has fitted them.
    ///
    /// **This cannot assert parity, and the reason is the recorded departure itself**: production
    /// and ng part at any inbred sample. What it asserts instead is the *pattern* — identical at
    /// every outbred sample, different at the inbred one.
    ///
    /// **What it catches, measured rather than argued: the final pass reading sample 0's
    /// coefficient for every sample.** That is where a genotype is decided, so replacing the
    /// per-sample coefficient there makes sample 7's call the outbred one and this fixture fails.
    ///
    /// **What it does not catch is the same substitution inside the frequency loop**, and the
    /// reason is worth stating rather than leaving for somebody to rediscover: the loop's
    /// coefficient moves only the cohort's frequency trajectory, and over eight samples that
    /// movement does not reach a call. `summarise_condition`'s own
    /// `each_sample_is_scored_against_its_own_inbreeding_coefficient` is what pins that one, by
    /// reading the posterior row rather than the genotype.
    #[test]
    fn a_samples_own_inbreeding_coefficient_is_what_reaches_its_row() {
        let mut table = Vec::new();
        for _ in 0..4 {
            table.extend_from_slice(&row(-25.0, -8.0, 0.0));
        }
        for _ in 0..2 {
            table.extend_from_slice(&row(0.0, -8.0, -25.0));
        }
        // The two whose reads leave a nat of room, and whose coefficients differ.
        table.extend_from_slice(&row(-1.0, 0.0, -8.0));
        table.extend_from_slice(&row(-1.0, 0.0, -8.0));

        let coefficients = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.95];
        let locus = Locus {
            log_likelihoods: &table,
            samples: 8,
            alleles: 2,
            diversity: 1e-3,
            inbreeding: &coefficients,
        };
        let comparison = locus.both();

        assert_eq!(
            comparison.ng[..7],
            comparison.production[..7],
            "the seven outbred samples take the branch that short-circuits, so the two loops are \
             the same arithmetic there"
        );
        assert_eq!(
            comparison.ng[6],
            vec![1, 1],
            "sample 6 is outbred and its reads decide it"
        );
        assert_eq!(
            comparison.ng[7],
            vec![2, 0],
            "sample 7 has the identical reads and F = 0.95, and ng's corrected coefficient \
             overturns them"
        );
        assert_ne!(
            comparison.ng[7], comparison.production[7],
            "and this is where the two part — a loop that read row 0's coefficient for every row \
             would agree with production here and this assertion would fail"
        );
    }

    /// **The two loops run under the same threshold and the same cap, and nothing but this
    /// asserts it.**
    ///
    /// They are two independent pairs of constants — production's in `posterior_engine`, ng's in
    /// `inference` — equal today because two edits agree. Every fixture here runs both at their
    /// defaults, so a change to one side alone would turn every parity assertion in this file
    /// into a comparison of two differently-configured loops, and nothing would say so.
    #[test]
    fn the_two_loops_run_under_the_same_threshold_and_the_same_cap() {
        use crate::ng::calling::inference::{DEFAULT_CONVERGENCE_THRESHOLD, DEFAULT_MAX_PASSES};
        use crate::var_calling::posterior_engine::{
            DEFAULT_CONVERGENCE_THRESHOLD as PRODUCTION_THRESHOLD,
            DEFAULT_MAX_ITERATIONS as PRODUCTION_CAP,
        };

        assert_eq!(
            DEFAULT_CONVERGENCE_THRESHOLD, PRODUCTION_THRESHOLD,
            "the two loops stop at different movements, so every fixture in this file is \
             comparing two differently-configured loops"
        );
        assert_eq!(
            DEFAULT_MAX_PASSES, PRODUCTION_CAP,
            "the two loops give up after different numbers of passes"
        );
    }

    /// **A locus the reads leave undecided is decided identically, and the two pass counts differ
    /// by exactly one for a reason both loops write down.**
    ///
    /// Every sample's three genotypes are within a nat of each other, so the prior carries the
    /// locus outright and any difference in either loop's arithmetic has nothing to hide behind.
    /// It takes 35 passes to settle, which is what makes the counting worth comparing here and
    /// nowhere else: on a locus that settles in two, an off-by-one is most of the count.
    ///
    /// **The two loops run the same number of scoring passes and report different numbers.**
    /// Both begin with one E-step on the reads alone, before any prior — production's
    /// `EmStepPhase::FirstIteration`, ng's initialisation. Production counts it as iteration 1;
    /// ng does not count it at all, because `passes` counts the passes that had a prior
    /// (`run_frequency_loop`'s own comment). So production's count is ng's plus one, always, and
    /// that is asserted rather than tolerated: a real difference in where the two stopped would
    /// break this by more than one, and a change to either definition would break it too.
    #[test]
    fn a_locus_the_reads_leave_open_is_decided_identically_and_counts_its_passes_one_apart() {
        let mut table = Vec::new();
        for sample in 0..8 {
            let tilt = f64::from(sample) * 0.05;
            table.extend_from_slice(&row(-0.3 + tilt, 0.0, -0.6 - tilt));
        }
        let locus = Locus {
            log_likelihoods: &table,
            samples: 8,
            alleles: 2,
            diversity: 1e-3,
            inbreeding: &OUTBRED[..8],
        };
        let comparison = locus.both();
        comparison.assert_parity();
        assert!(
            comparison.ng_passes > 10,
            "a locus that settled in a couple of passes would make the count comparison below \
             almost vacuous; this one took {}",
            comparison.ng_passes
        );
        assert_eq!(
            comparison.production_iterations,
            comparison.ng_passes + 1,
            "production counts the prior-free first E-step and ng does not, so the two differ by \
             exactly that one pass — anything else means they stopped in different places"
        );
    }
}
