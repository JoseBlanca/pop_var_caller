use super::*;
use crate::calling::likelihood::TEST_TRACT_LEN;
use crate::parameter_estimation::joint::stratum_fits::AlleleLengthVariationRung;

use proptest::prelude::*;

use crate::calling::GenericObservation;
use crate::calling::genotype_prior::{MarginalizedDirichletPrior, SeedRegime, SpectrumSeed};
use crate::calling::likelihood::ssr_emission::SsrCandidate;
use crate::calling::quality::MAX_GENOTYPE_QUALITY;
use crate::calling::{
    CandidateAlleles, ExpectedAlleleCopies, FrozenParameters, GenotypeIdx, GenotypeTable,
    LocusEvidence, LocusInference, ReadGroupCalibration, SampleGenotypeCall, SsrSampleEvidence,
    UNWRITTEN_SCRATCH_VALUE,
};
use crate::cohort_observation::PartialObservation;
use crate::locus_generation::{
    LocusKind, ReadWitness, SequenceObservation, SsrDetail, WitnessedLocusPositions,
};
use crate::parameter_estimation::Provenance;
use crate::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::types::{
    AlleleId, ContigId, GenomeRegion, Genotype, Motif, Phred, Position, ReadGroupId, SummedLogError,
};
use std::num::NonZeroU32;
use std::sync::Arc;

use super::super::{CallingLoopConfig, DEFAULT_CONVERGENCE_THRESHOLD};
use crate::calling::MIN_CONTAMINANT_FREQUENCY;
use crate::calling::tests::one_batch;
use crate::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::read::input::read_groups::ReadGroups;
use crate::types::BatchId;

use crate::calling::likelihood::ssr_emission::{
    StutterSubstitutionEmission, StutterSubstitutionScratch,
};
use crate::calling::{ContaminationView, EmissionCost};
use crate::parameter_estimation::joint::contamination::ContaminationSource;

mod contamination;
mod convergence;
mod cost;
mod driver;
mod final_pass;
mod frequency_loop;
mod repeat_tract;

/// The repeat unit every repeat tract in this file is built on — **a dinucleotide**, so that
/// a whole repeat and a single base are different steps and the part-repeat half of the
/// stutter distribution is not hidden.
const TRACT_MOTIF: &[u8] = b"AT";

/// The two candidate lengths every repeat tract here is called over, in candidate-table
/// order: **6 whole repeats and 7**. Entry 0 is the reference allele's, which is what the
/// prior's allele length variation measures its offsets from.
///
/// **Different from each other**, so a lookup that read one candidate's repeat count for
/// another lands in a different stratum rather than the same one. **Neither is what a
/// reader could derive from the other's bases**, which
/// `an_interrupted_candidates_supplied_repeat_count_is_what_picks_its_stratum` is about.
const TRACT_CANDIDATE_REPEATS: [u32; 2] = [6, 7];

/// **How many libraries a tract fixture's run has, and it is not two.** The parameter table
/// a tract is scored from is `read groups × candidates`, filled read-group-major; at an
/// equal shape a table filled the other way round is the same length and the same set of
/// cells, so a transposition passes every shape check. At three libraries against two
/// candidates the two orders disagree cell by cell.
///
/// Only the first library sends any read below. The other two are in the run, which is what
/// makes them visible to a tract's warrant fold — that fold runs over the run's read-group
/// axis rather than over the reads that arrived.
const TRACT_READ_GROUPS: usize = 3;

/// The two candidates' repeat counts. **Supplied, not selected** — repeat-tract candidate
/// selection is unwritten, so a tract fixture states its candidates' repeat counts and a
/// reader must not take them for a step's output.
fn tract_repeat_counts() -> Vec<NonZeroU32> {
    TRACT_CANDIDATE_REPEATS
        .iter()
        .map(|count| NonZeroU32::new(*count).expect("a candidate always holds a repeat"))
        .collect()
}

/// The bases of a tract holding `repeats` whole copies of the motif.
fn tract_bases(repeats: u32) -> Vec<u8> {
    TRACT_MOTIF.repeat(repeats as usize)
}

fn tract_detail() -> SsrDetail {
    SsrDetail {
        motif: Motif::new(TRACT_MOTIF).expect("a dinucleotide motif"),
        left_flank: Box::from(b"CCCGGG".as_slice()),
        right_flank: Box::from(b"TTTAAA".as_slice()),
    }
}

/// The candidate table the two lengths above spell, with the shorter as the reference.
fn tract_alleles() -> CandidateAlleles {
    let mut alleles = CandidateAlleles::new(
        tract_bases(TRACT_CANDIDATE_REPEATS[0]).into_boxed_slice(),
        LocusKind::Ssr(tract_detail()),
    );
    alleles.admit(tract_bases(TRACT_CANDIDATE_REPEATS[1]).into_boxed_slice());
    alleles
}

/// `reads` reads that spanned the whole tract and showed `repeats` whole copies of it, from
/// the run's first library.
///
/// **The per-read error sum is set to a plausible value and nothing on this path reads
/// it**, which is worth saying because it is what a reader would assume decides these
/// calls. A tract's row is scored from a stutter model and a per-base substitution rate;
/// `q_sum` is the SNP/indel emission's charge and neither `likelihood::ssr` nor its
/// emission model touches it. Zeroing it here would change no number below.
fn tract_reads(repeats: u32, reads: u32) -> SequenceObservation {
    SequenceObservation {
        flank_deletion: crate::locus_generation::FlankDeletion::NONE,
        held_flank: None,
        bases: tract_bases(repeats).into_boxed_slice(),
        read_witness: ReadWitness::Complete,
        read_group: ReadGroupId(0),
        num_obs: reads,
        num_fwd: reads / 2,
        q_sum: SummedLogError::from_nats(-10.0 * f64::from(reads)),
        mapq_sum: 60 * reads,
        mapq_sum_sq: u64::from(reads) * 3_600,
        placed_left: reads / 2,
        num_low_quality: 0,
        chain_ids: Vec::new(),
    }
}

/// **`reads` reads that ran out inside the tract**, having seen `repeats` whole copies of it
/// before they did — a lower bound on the tract's length rather than a measurement of it.
///
/// **Over half the reads overlapping a 60-base tract are these on real data**, and they are
/// scored by a different term of the row from a read that spanned the whole of it. A fixture
/// of complete reads only is the one shape in which counting the observations and counting
/// the *complete* observations are the same number.
fn tract_reads_that_ran_out(repeats: u32, reads: u32) -> SequenceObservation {
    let bases = tract_bases(repeats);
    let covered = u16::try_from(bases.len()).expect("a short fixture tract");
    SequenceObservation {
        flank_deletion: crate::locus_generation::FlankDeletion::NONE,
        held_flank: None,
        read_witness: ReadWitness::from_left(
            covered,
            crate::locus_generation::LocusLen::from_positions(u64::from(covered) + 4),
        )
        .expect("a reach that covers part of the locus"),
        bases: bases.into_boxed_slice(),
        read_group: ReadGroupId(0),
        num_obs: reads,
        num_fwd: reads / 2,
        q_sum: SummedLogError::from_nats(-10.0 * f64::from(reads)),
        mapq_sum: 60 * reads,
        mapq_sum_sq: u64::from(reads) * 3_600,
        placed_left: reads / 2,
        num_low_quality: 0,
        chain_ids: Vec::new(),
    }
}

/// **A repeat fit that reached every stratum this file's tracts are called over**, each
/// with its own slippage numbers and its own allele length variation.
///
/// **Three sources of coincidence are removed on purpose**, each because a mutation
/// survived it:
///
/// - **the two strata share no number.** Level, direction split and fall-off all differ, so
///   a lookup that answered from the wrong stratum changes a score rather than only the one
///   number a reader happens to check;
/// - **the run's libraries do not share one slippage group.** Library 0 sits in group 0 and
///   libraries 1 and 2 in group 1, whose numbers differ — so a parameter table filled
///   candidate-major instead of read-group-major hands library 0's reads another group's
///   polymerase;
/// - **no allele length variation is a palindrome, no two of its classes share a weight, and no two
///   adjacent pairs of classes are in the same ratio.** Each is heavier on the contraction
///   side, as every real repeat fit is, so reversing one is a different prior; reading class
///   *i* for class *j* is a different number; and — the repair a review forced — a spectrum
///   read one repeat off centre gives the candidates a *different* pair of shares, where a
///   geometric tail would have given the same pair.
///
/// Every warrant here is the cell's own, so a tract scored under it is
/// [`Provenance::FittedHere`] — which is what makes the calibration's own `Defaulted`
/// visible if it ever entered the tract's fold.
fn tract_strata() -> StratumFits {
    tract_strata_describing(std::collections::BTreeMap::from([
        (ReadGroupId(0), 0),
        (ReadGroupId(1), 1),
        (ReadGroupId(2), 1),
    ]))
}

/// The same fit, over whichever of the run's libraries it claims to describe — what the
/// fixture about a library the fit never saw needs.
fn tract_strata_describing(
    slippage_group_of_each_library: std::collections::BTreeMap<ReadGroupId, u32>,
) -> StratumFits {
    use crate::parameter_estimation::joint::census::Stratum;
    use crate::parameter_estimation::joint::share_curve::ShareSource;
    use crate::parameter_estimation::joint::slippage_curve::LevelSource;
    use crate::parameter_estimation::joint::ssr_fit::{
        LevelProvenance, ShareProvenance, SharesProvenance, Slippage, StratumFit, StratumOutcome,
    };

    let level = LevelProvenance {
        source: LevelSource::Cell,
        curve: None,
        reach: None,
        slipped_reads: Some(400.0),
    };
    let share = ShareProvenance {
        source: ShareSource::Stratum,
        curve: None,
        reach: None,
    };
    // One entry per slippage group: group 0 is library 0's, group 1 is libraries 1 and 2's.
    let slippage_of_each_group = |level_in_group_0: f64| {
        vec![
            Some(Slippage {
                level: level_in_group_0,
                shorter_share: 0.83,
                fall_off: 0.25,
            }),
            Some(Slippage {
                level: level_in_group_0 * 1.75,
                shorter_share: 0.62,
                fall_off: 0.41,
            }),
        ]
    };
    let fitted_stratum = |repeats: u64,
                          level_in_group_0: f64,
                          allele_length_variation: Vec<f64>,
                          concentration: f64| {
        StratumOutcome::Fitted(Box::new(StratumFit {
            standard_errors: None,
            stratum: Stratum {
                period: 2,
                reference_repeats: repeats,
            },
            slippage: slippage_of_each_group(level_in_group_0),
            allele_length_variation,
            concentration: Some(concentration),
            samples_fitted_on: None,
            tracts_of_its_own: 40,
            reads_crossing: 400,
            level_provenance: vec![Some(level), Some(level)],
            shares_provenance: vec![
                Some(SharesProvenance {
                    slipped_reads: Some(400.0),
                    shorter_share: share,
                    fall_off: share,
                }),
                Some(SharesProvenance {
                    slipped_reads: Some(400.0),
                    shorter_share: share,
                    fall_off: share,
                }),
            ],
        }))
    };
    StratumFits::over(
        &[
            // The weights run over offsets −2 … +2 whole repeats from the stratum's own
            // reference tract length.
            //
            // **No two adjacent pairs of these weights are in the same ratio**, which is a
            // deliberate repair: the first version ran `0.10, 0.25, 0.45, 0.15, 0.05`,
            // whose upper tail falls by a factor of three at each step — so a spectrum read
            // one repeat off centre gave the two candidates the *same* pair of shares,
            // 0.75 and 0.25, and the re-centring this fixture exists to catch was
            // invisible. Here offsets {0, +1} give 0.8 and 0.2 where {+1, +2} give 0.6875
            // and 0.3125.
            fitted_stratum(4, 0.03, vec![0.12, 0.27, 0.42, 0.14, 0.05], 16.0),
            fitted_stratum(6, 0.04, vec![0.10, 0.30, 0.44, 0.11, 0.05], 20.0),
            fitted_stratum(7, 0.06, vec![0.09, 0.21, 0.43, 0.19, 0.08], 25.0),
        ],
        slippage_group_of_each_library,
    )
}

/// **The record a called repeat tract carries**, which every tract has and no ordinary site
/// does — so a call to this on a SNP/indel locus's result is a test asking the wrong path.
fn tract_record(inference: &LocusInference) -> RepeatTractProvenance {
    inference
        .repeat_tract
        .expect("every called repeat tract carries what its parameters rested on")
}

/// Which rung of the tract ladder that record's prior shape came from.
fn rung_of(inference: &LocusInference) -> AlleleLengthVariationRung {
    tract_record(inference).allele_length_variation_rung()
}

/// The run's libraries, none of them calibrated — a tract's row never reads a calibration,
/// so this says nothing about the tract and everything about whether the warrant fold has
/// wandered onto the wrong list.
fn tract_libraries() -> Vec<ReadGroupCalibration> {
    vec![ReadGroupCalibration::unscaled(); TRACT_READ_GROUPS]
}

/// **A fitted per-base substitution rate for every `(library, stratum)` cell a tract's
/// parameter table covers** — every library of the run against every stratum this file's
/// tracts reach.
///
/// **Every cell, not only the ones whose reads arrived**, because the table is built on the
/// run's read-group axis: a run that fitted only the library that sent reads still has
/// cells for the other two, and those cells fall to the stated constant and pull the
/// locus's warrant down with them. `a_library_that_sent_this_tract_no_read_still_weakens_it`
/// is that case; this is its opposite.
///
/// The rate differs per stratum, so a lookup keyed by the wrong repeat count is a different
/// number rather than the same one.
fn tract_substitution_rates() -> std::collections::BTreeMap<
    crate::parameter_estimation::repeat_strata::StratumKey,
    crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
> {
    substitution_rates_over(&[4, 6, 7], TRACT_READ_GROUPS)
}

/// The same, over the strata and libraries named — so a fixture can leave a cell out.
fn substitution_rates_over(
    strata: &[u32],
    read_groups: usize,
) -> std::collections::BTreeMap<
    crate::parameter_estimation::repeat_strata::StratumKey,
    crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
> {
    use crate::parameter_estimation::Estimate;
    use crate::parameter_estimation::repeat_strata::{
        RepeatCount, Stratum as SsrStratum, StratumKey,
    };
    use crate::types::{ErrorRate, SsrPeriod};

    let period = SsrPeriod::try_new(2).expect("a dinucleotide");
    let mut rates = std::collections::BTreeMap::new();
    for library in 0..read_groups {
        for repeats in strata {
            rates.insert(
                StratumKey {
                    read_group: ReadGroupId(u32::try_from(library).expect("a small library index")),
                    stratum: SsrStratum::new(period, RepeatCount(*repeats)),
                    ploidy: diploid(),
                },
                Estimate {
                    // A different rate per stratum, so a lookup keyed by the wrong repeat
                    // count is a different number.
                    value: ErrorRate::try_new(0.001 * f64::from(*repeats)).expect("a probability"),
                    provenance: Provenance::FittedHere,
                    observations: 4_000,
                    standard_error: None,
                },
            );
        }
    }
    rates
}

/// A genotype prior that writes numbers the test chose.
///
/// **A stub is what makes the E-step's own arithmetic checkable by hand.** Both shipped
/// priors go through `lgamma`, so a hand-computed expectation for one of them would be
/// a transcription of a calculator rather than a check on this function; the real ones
/// are used below on the properties they can be held to instead. The concentration this
/// sees is not thrown away — it is the scratch's own buffer, so the test reads it back
/// and checks it there.
struct FixedLogPriors(Vec<f64>);

impl GenotypePriorModel for FixedLogPriors {
    fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, _inbreeding: InbreedingF) {
        let (_workspace, out) = row.scratch_and_out();
        assert_eq!(
            out.len(),
            self.0.len(),
            "the fixture's row is the wrong width"
        );
        for (slot, &value) in out.iter_mut().zip(&self.0) {
            *slot = LogProb(value);
        }
    }
}

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("a diploid")
}

/// A run whose repeat-tract substitution rates were never fitted — the empty map, which is
/// what a fixture calling no tract needs, and what
/// `FrozenParameters::ssr_substitution_rate_at` answers `None` from.
///
/// A `static` rather than a function, so that a call site can borrow it for as long as the
/// parameters live: `BTreeMap::new` is a `const fn`, and a temporary would be freed at the
/// end of the statement that built the view.
static NO_SUBSTITUTION_RATES: std::collections::BTreeMap<
    crate::parameter_estimation::repeat_strata::StratumKey,
    crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
> = std::collections::BTreeMap::new();

/// A site quality standing in for one the worker computed, where what the test is about
/// is something else.
///
/// **Deliberately not zero.** Zero is a real answer — a locus at which nobody carries
/// anything — so a test whose fixture wrote zero would pass against an implementation
/// that lost the value on the way in.
fn a_worker_written_site_quality() -> Phred {
    Phred::try_new(37.0).expect("a legal quality, and below the site-quality ceiling")
}

fn outbred() -> InbreedingF {
    InbreedingF::try_new(0.0).expect("an outbred sample")
}

/// A SNP locus over `alternatives + 1` alleles, with its diploid genotype table.
fn generic_locus(alternatives: usize) -> (CandidateAlleles, Arc<GenotypeTable>) {
    let mut alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
    for base in b"TCG".iter().take(alternatives) {
        alleles.admit(Box::from(&[*base][..]));
    }
    let table = GenotypeTable::build(diploid(), alleles.len());
    (alleles, table)
}

/// A configuration that ships, with the pass cap moved.
fn capped_at(max_passes: u32) -> RunnableCallingLoopConfig {
    CallingLoopConfig {
        max_passes: NonZeroU32::new(max_passes).expect("a cap of at least one pass"),
        ..CallingLoopConfig::DEFAULT
    }
    .validate()
    .expect("only the cap moved, and it is not a value validation refuses")
}

/// Run the frequency loop over a cohort whose genotype likelihoods are `likelihoods`, one
/// row per sample, every sample outbred — and hand back the scratch so a test can read
/// what the loop left in it: the cohort's expected copies, the previous pass's row, and
/// the last sample scored.
///
/// For a cohort whose samples differ in inbreeding, use [`loop_over_inbred`].
fn loop_over(
    likelihoods: &[Vec<f64>],
    seed_concentration: &[f64],
    alleles: &CandidateAlleles,
    view: &GenotypeTableView<'_>,
    config: &RunnableCallingLoopConfig,
) -> (FrequencyLoopOutcome, CallingScratch<()>) {
    let outbred_cohort = vec![outbred(); likelihoods.len()];
    loop_over_inbred(
        likelihoods,
        &outbred_cohort,
        seed_concentration,
        alleles,
        view,
        config,
    )
}

/// [`loop_over`] with an inbreeding coefficient chosen per sample, in the run's sample
/// order.
fn loop_over_inbred(
    likelihoods: &[Vec<f64>],
    inbreeding_by_sample: &[InbreedingF],
    seed_concentration: &[f64],
    alleles: &CandidateAlleles,
    view: &GenotypeTableView<'_>,
    config: &RunnableCallingLoopConfig,
) -> (FrequencyLoopOutcome, CallingScratch<()>) {
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(likelihoods.len(), alleles, view);
    // Every sample of this fixture is callable, so the rows are the run's samples
    // one for one — which is what a locus with nothing set aside looks like.
    for (sample, &inbreeding) in inbreeding_by_sample.iter().enumerate() {
        scratch.claim_row_for(sample, inbreeding);
    }
    scratch
        .seed_concentration_mut()
        .copy_from_slice(seed_concentration);
    for (sample, row) in likelihoods.iter().enumerate() {
        for (slot, &value) in scratch
            .sample_genotype_likelihoods_mut(sample)
            .iter_mut()
            .zip(row)
        {
            *slot = LogProb(value);
        }
    }
    let outcome = run_frequency_loop(
        &mut scratch,
        view,
        &MarginalizedDirichletPrior,
        config,
        None,
    );
    (outcome, scratch)
}

fn locus_region() -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(3),
        start: Position(940),
        end: Position(940),
    }
}

/// One `(allele, read group)` row of a sample's evidence, with the two counts the
/// artifact summary reads.
///
/// **`q_sum` is filled with something that would be wrong if anything read it**, rather
/// than with zero: the final pass never looks at the error mass, and a fixture of zeroes
/// would not say so.
fn observation(
    allele: u16,
    read_group: u32,
    num_reads: u32,
    forward_reads: u32,
    placed_left_reads: u32,
) -> GenericObservation {
    GenericObservation {
        allele: AlleleId(allele),
        read_group: ReadGroupId(read_group),
        num_reads,
        q_sum: -3.0 * f64::from(num_reads),
        forward_reads: f64::from(forward_reads),
        placed_left_reads: f64::from(placed_left_reads),
        low_quality_reads: 0,
    }
}

/// A sample whose reads are `rows` and which the candidate step did not set aside.
fn called_sample<'a>(rows: &'a [GenericObservation]) -> GenericLocusSample<'a> {
    GenericLocusSample {
        evidence: GenericSampleEvidence::new(rows, 0.0, &[]),
        genotype_must_be_missing: false,
    }
}

/// A neutral panel's fitted spectrum at one variant per kilobase — human diversity, and
/// the seed the site quality's own tests use.
fn human_like_seed() -> SpectrumSeed {
    SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape)
}

/// **The whole of one locus over a hand-built likelihood table**: the frequency loop, then
/// the final pass over the same scratch — the order and the sharing the design requires,
/// since the pass reads the posterior row the loop leaves and the likelihood table beside
/// it. **The table is supplied rather than assembled here**, so nothing in this fixture
/// re-assembles it; the driver's own tests are where a contaminated locus's per-pass
/// assembly is exercised.
///
/// `likelihoods` holds **one row per sample of the run**, and a sample the candidate step
/// ruled uncallable simply does not get a scratch row — so its row here is never read,
/// which is what the driver does with it too.
fn call_locus(
    likelihoods: &[Vec<f64>],
    evidence: &LocusEvidence<'_>,
    seed_concentration: &[f64],
    alleles: &CandidateAlleles,
    view: &GenotypeTableView<'_>,
    model: &dyn GenotypePriorModel,
) -> LocusInference {
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = vec![outbred(); likelihoods.len()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = FrozenParameters::uncontaminated(
        &calibration,
        &inbreeding,
        human_like_seed(),
        &strata,
        &NO_SUBSTITUTION_RATES,
        view.ploidy(),
    );
    let callable: Vec<usize> = (0..evidence.sample_count())
        .filter(|run_sample| is_callable(evidence, *run_sample))
        .collect();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(callable.len(), alleles, view);
    for &run_sample in &callable {
        scratch.claim_row_for(run_sample, outbred());
    }
    scratch
        .seed_concentration_mut()
        .copy_from_slice(seed_concentration);
    for (row, &run_sample) in callable.iter().enumerate() {
        for (slot, &value) in scratch
            .sample_genotype_likelihoods_mut(row)
            .iter_mut()
            .zip(&likelihoods[run_sample])
        {
            *slot = LogProb(value);
        }
    }
    let outcome = run_frequency_loop(
        &mut scratch,
        view,
        model,
        &RunnableCallingLoopConfig::default(),
        None,
    );
    summarise_final_pass(
        &mut scratch,
        view,
        evidence,
        &parameters,
        model,
        alleles.clone(),
        outcome,
        Provenance::FittedHere,
        // **The fixture's own, because this helper drives the final pass and not the
        // loop**: nothing here gathers a tract's scoring parameters, so there is no filled
        // table to read the counts off. A tract must carry a record and an ordinary site
        // must not, so the locus's kind decides which.
        (!matches!(alleles.kind(), LocusKind::Generic))
            .then(a_tract_record_this_fixture_did_not_gather),
    )
}

/// **A stand-in for what a tract's parameters rested on**, for a fixture that never
/// gathered any — the counts are the fixture's, not a measurement, and no test of this
/// helper's may assert them.
///
/// The driver's own tests are where the counts are checked against a gather.
fn a_tract_record_this_fixture_did_not_gather() -> RepeatTractProvenance {
    RepeatTractProvenance::new(AlleleLengthVariationRung::StatedFlat, 0, 0, 0, 0, false)
}

/// The call this fixture expects, unwrapped — a `Missing` here is the test's failure and
/// not something to match on.
fn called(inference: &LocusInference, sample: usize) -> (&Genotype, Phred) {
    match &inference.per_sample[sample] {
        SampleGenotypeCall::Called {
            genotype,
            genotype_quality,
            ..
        } => (genotype, *genotype_quality),
        SampleGenotypeCall::Missing => {
            panic!("sample {sample} was called, so it is not missing")
        }
    }
}
/// The shipped arm: repeat tracts scored by the shipped emission model, genotypes under
/// the marginalized Dirichlet prior.
fn shipped_arm() -> SummariseConditionLoop<StutterSubstitutionEmission, MarginalizedDirichletPrior>
{
    SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior)
}

/// A worker's scratch, typed by the shipped emission model.
fn worker_scratch() -> CallingScratch<StutterSubstitutionScratch> {
    CallingScratch::default()
}

/// One run's frozen parameters over `samples` outbred samples and one read group, with
/// nothing contaminated.
fn uncontaminated_run<'a>(
    calibration: &'a [ReadGroupCalibration],
    inbreeding: &'a [InbreedingF],
    strata: &'a StratumFits,
    substitution: &'a std::collections::BTreeMap<
        crate::parameter_estimation::repeat_strata::StratumKey,
        crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
    >,
    ploidy: Ploidy,
) -> FrozenParameters<'a> {
    FrozenParameters::uncontaminated(
        calibration,
        inbreeding,
        human_like_seed(),
        strata,
        substitution,
        ploidy,
    )
}

/// The genotype a sample was called, as allele ids.
fn called_alleles(inference: &LocusInference, sample: usize) -> Vec<u16> {
    inference.per_sample[sample]
        .genotype()
        .unwrap_or_else(|| panic!("sample {sample} was called missing"))
        .alleles()
        .iter()
        .map(|allele| allele.get())
        .collect()
}

/// A read group whose fit found the fraction below, measured on real evidence.
fn contaminated_at(fraction: f64) -> ContaminationView {
    ContaminationView {
        fraction,
        markers_with_reads: 400,
        reads_on_markers: 1_000,
        source: ContaminationSource::ThisReadGroupsReads,
    }
}

/// The frozen parameters for a contaminated run, over the default batching — one batch
/// holding all of it, which is what a run that declared nothing gets.
fn contaminated_run<'a>(
    calibration: &'a [ReadGroupCalibration],
    contamination: &'a [ContaminationView],
    batching: &'a SequencingBatches,
    inbreeding: &'a [InbreedingF],
    strata: &'a StratumFits,
) -> FrozenParameters<'a> {
    FrozenParameters::new(
        calibration,
        contamination,
        batching,
        inbreeding,
        human_like_seed(),
        strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    )
}
