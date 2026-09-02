//! **ng calls genotypes** — per-sample reads in at one cohort locus, called genotypes out, with
//! the answers derived by hand. Both paths: an ordinary SNP or indel, and a repeat tract.
//!
//! Every earlier test of this loop hands it a likelihood table, an evidence view or a candidate
//! set built for the occasion. This one hands it a cohort locus and takes genotypes out:
//!
//! ```text
//! SNP/indel: per-sample observations → ClosedLocus → CohortObservation::over → select_generic
//!                                   → shape_generic_locus → call_locus → LocusInference
//! repeat tract: per-sample observations + supplied candidates
//!                                   → shape_ssr_locus      → call_locus → LocusInference
//! ```
//!
//! # What this fixture supplies and what it runs, said exactly
//!
//! **Supplied**, and each for a reason:
//!
//! - **the per-sample observations** — what each sample's reads showed at the locus, as the
//!   locus generator emits them. Turning aligned reads into observations is step 5's
//!   and outside this plan (`calling_loop.md`'s Scope); a fixture that ran it would be testing
//!   that step here;
//! - **the `ClosedLocus`** the merge is handed, on the SNP/indel path. The chaining walk that
//!   groups overlapping observations into loci — `LocusCloser` — is not run. What [`merge`] does
//!   reproduce is the **keep rule**, `MinAltReads::DEFAULT` asked of each sample separately, so a
//!   fixture the real walk would have discarded as too quiet cannot pass here unnoticed;
//! - **at a repeat tract, the candidates themselves and each candidate's repeat count.** The
//!   repeat-tract half of candidate selection is not written, so there is no step to run — and a
//!   later reader must not take these candidate sets for a step's output. The tract section
//!   below says it again where the fixtures are.
//! - the run's frozen parameters and the loop's configuration, which are a run's inputs rather
//!   than a locus's.
//!
//! **Run rather than supplied**: on the SNP/indel path the merge's allele unification and read
//! attribution (`CohortObservation::over`), candidate selection (`select_generic`) and the input
//! edge (`shape_generic_locus`); at a tract the input edge (`shape_ssr_locus`); and on both, the
//! loop (`call_locus`).
//!
//! # Why a test binary rather than a module inside the library
//!
//! It imports only what the crate exports, so it is also a check that both paths' seams are
//! `pub` — the 56 items it names, not every seam in the middle of them. Half of those are the
//! repeat fit's own output types, which a fixture has to build by hand because the pre-pass that
//! produces them is a different subsystem.

use pop_var_caller::ng::calling::allele_candidates::generic::select_generic;
use pop_var_caller::ng::calling::allele_candidates::{CandidateSelectionConfig, SelectionScratch};
use pop_var_caller::ng::calling::evidence_shaping::{
    GenericEvidenceScratch, shape_generic_locus, shape_ssr_locus,
};
use pop_var_caller::ng::calling::genotype_prior::{
    MarginalizedDirichletPrior, SeedRegime, SpectrumSeed,
};
use pop_var_caller::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use pop_var_caller::ng::calling::inference::{CallingLoopConfig, LocusGenotyper};
use pop_var_caller::ng::calling::likelihood::ssr_emission::{
    StutterSubstitutionEmission, StutterSubstitutionScratch,
};
use pop_var_caller::ng::calling::{
    CallingScratch, CandidateAlleles, ContaminationView, FrozenParameters, GenericLocusSample,
    LocusInference, ReadGroupCalibration, RepeatTractProvenance, SsrSampleEvidence,
};
use pop_var_caller::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation, SsrDetail,
};
use pop_var_caller::ng::parameter_estimation::Provenance;
use pop_var_caller::ng::parameter_estimation::joint::census::Stratum as FitStratum;
use pop_var_caller::ng::parameter_estimation::joint::contamination::ContaminationSource;
use pop_var_caller::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use pop_var_caller::ng::parameter_estimation::joint::share_curve::ShareSource;
use pop_var_caller::ng::parameter_estimation::joint::slippage_curve::LevelSource;
use pop_var_caller::ng::parameter_estimation::joint::ssr_fit::{
    LevelProvenance, ShareProvenance, SharesProvenance, Slippage, StratumFit, StratumOutcome,
};
use pop_var_caller::ng::parameter_estimation::joint::stratum_fits::{
    LengthSpectrumRung, StratumFits,
};
use pop_var_caller::ng::parameter_estimation::ssr::{
    RepeatCount, Stratum as SsrStratum, StratumKey,
};
use pop_var_caller::ng::run::cohort_merge::MinAltReads;
use pop_var_caller::ng::run::cohort_merge::build::CohortObservation;
use pop_var_caller::ng::run::cohort_merge::close::{ClosedLocus, SampleMembers, Verdict};
use pop_var_caller::ng::types::{
    AlleleId, ContigId, ErrorRate, GenomeRegion, InbreedingF, Motif, Ploidy, Position, ReadGroupId,
    SsrPeriod, SummedLogError,
};

use pop_var_caller::ng::parameter_estimation::Estimate;
use std::collections::BTreeMap;
use std::num::NonZeroU32;

/// The one position every locus here sits at.
fn region() -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(1),
        start: Position(1_000),
        end: Position(1_000),
    }
}

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("a diploid")
}

/// **A sample's reads at the locus, as the generator emits them**: how many showed each base.
///
/// The per-read error sum is what the SNP/indel emission charges an observation, so it is set
/// from a real per-read quality rather than left at zero — `q_sum` of `−3 · reads` is about
/// Phred 13 a read, low enough to be ordinary and high enough that twenty reads decide a
/// genotype.
fn showed(bases: &[u8], reads: u32) -> SequenceObservation {
    showed_from(bases, reads, ReadGroupId(0))
}

/// The same, from a named library — what the warrant tests need, since a locus's warrant is
/// folded over the read groups whose reads reached it.
fn showed_from(bases: &[u8], reads: u32, read_group: ReadGroupId) -> SequenceObservation {
    SequenceObservation {
        bases: bases.to_vec(),
        read_witness: ReadWitness::Complete,
        read_group,
        num_obs: reads,
        num_fwd: reads / 2,
        q_sum: SummedLogError::from_nats(-3.0 * f64::from(reads)),
        mapq_sum: 60 * reads,
        mapq_sum_sq: u64::from(reads) * 3_600,
        placed_left: reads / 2,
        chain_ids: Vec::new(),
    }
}

/// One sample's locus, over the reference base `A`.
fn sample_locus(observations: Vec<SequenceObservation>) -> SampleLocusObservations {
    SampleLocusObservations {
        region: region(),
        reference_bases: b"A".to_vec(),
        observations,
        reads_without_observation: 0,
        reads_discarded_by_cap: 0,
        kind: LocusKind::Generic,
    }
}

/// **The merge's allele unification and read attribution, run** over a closed locus built here.
///
/// **A sample with no observations gets no `SampleMembers` entry**, which is that type's own
/// contract — *"a sample with nothing here has no `SampleMembers` at all"* — and is what makes
/// the covering-sample list shorter than the run's. The first draft pushed an entry per sample
/// unconditionally, and the silent-sample test below then exercised nothing: with every sample
/// covering, the merge's index and the run's are the same number, and an input edge that
/// confused the two passed all seven tests.
///
/// **The keep rule is reproduced rather than assumed.** `LocusCloser`'s chaining walk is not run
/// here, but its verdict rests on one question asked of each sample separately — at least
/// `MinAltObs::DEFAULT` non-reference reads, and at least `MinAltReadShare::DEFAULT` of that
/// sample's own compared reads. A fixture no sample reaches it with is one the real walk would
/// discard as too quiet, so this refuses it rather than handing the loop a locus that could not
/// arise.
///
/// # Panics
///
/// If no sample covers the locus, or if no sample reaches the keep rule.
fn merge(per_sample: &[SampleLocusObservations]) -> CohortObservation {
    let members: Vec<SampleMembers<'_>> = per_sample
        .iter()
        .enumerate()
        .filter(|(_, locus)| !locus.observations.is_empty())
        .map(|(sample, locus)| SampleMembers {
            sample,
            observations: std::slice::from_ref(locus),
        })
        .collect();
    assert!(
        !members.is_empty(),
        "a locus nobody covered is not a locus the merge is handed"
    );

    let non_reference_of = |locus: &SampleLocusObservations| -> u32 {
        locus
            .observations
            .iter()
            .filter(|observation| observation.bases.as_slice() != b"A")
            .map(|observation| observation.num_obs)
            .sum()
    };
    let some_sample_reached_the_rule = per_sample.iter().any(|locus| {
        let compared: u32 = locus
            .observations
            .iter()
            .map(|observation| observation.num_obs)
            .sum();
        MinAltReads::DEFAULT.reached_by(non_reference_of(locus), compared)
    });
    assert!(
        some_sample_reached_the_rule,
        "no sample reaches the merge's keep rule at this locus, so the walk would call it too \
         quiet and never build it — a fixture asserting `Build` past that would be testing a \
         locus no run can produce"
    );

    CohortObservation::over(&ClosedLocus {
        region: region(),
        members,
        non_reference_reads: per_sample.iter().map(non_reference_of).sum(),
        verdict: Verdict::Build,
    })
}

/// A run of `samples` outbred diploids sequenced from one library, with nothing contaminated
/// and no repeat tract anywhere.
struct Run {
    calibration: Vec<ReadGroupCalibration>,
    inbreeding: Vec<InbreedingF>,
    strata: StratumFits,
    substitution: BTreeMap<StratumKey, Estimate<ErrorRate>>,
}

impl Run {
    fn of(samples: usize, calibration: &[ReadGroupCalibration]) -> Self {
        Self {
            calibration: calibration.to_vec(),
            inbreeding: vec![InbreedingF::try_new(0.0).expect("an outbred sample"); samples],
            // **A gather over no outcomes rather than nothing at all**, so a lookup answers
            // *no such stratum* — which is what a run with no repeat tracts has.
            strata: StratumFits::over(&[], BTreeMap::new()),
            substitution: BTreeMap::new(),
        }
    }

    fn parameters(&self) -> FrozenParameters<'_> {
        FrozenParameters::uncontaminated(
            &self.calibration,
            &self.inbreeding,
            // A human-like spectrum: one chromosome of reference belief and a thousandth of
            // that for the alternative, on the neutral shape.
            SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
            &self.strata,
            &self.substitution,
            diploid(),
        )
    }
}

/// Run the whole path over one cohort locus, under one library whose calibration was never
/// fitted, and hand back what came out.
fn call(per_sample: &[SampleLocusObservations]) -> LocusInference {
    call_with_calibration(per_sample, &[ReadGroupCalibration::defaulted()])
}

/// The same, over a run whose libraries carry the calibrations given.
fn call_with_calibration(
    per_sample: &[SampleLocusObservations],
    calibration: &[ReadGroupCalibration],
) -> LocusInference {
    let merged = merge(per_sample);
    let mut selection_scratch = SelectionScratch::new();
    let selection = select_generic(
        &merged,
        &CandidateSelectionConfig::DEFAULT,
        &mut selection_scratch,
    );

    let mut shaping = GenericEvidenceScratch::default();
    let mut views: Vec<GenericLocusSample<'_>> = Vec::new();
    let evidence = shape_generic_locus(
        &mut shaping,
        &merged,
        &selection,
        per_sample.len(),
        &mut views,
    );

    let run = Run::of(per_sample.len(), calibration);
    let parameters = run.parameters();
    let arm = SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior);
    let config = CallingLoopConfig::DEFAULT
        .validate()
        .expect("the shipped configuration");
    let mut scratch: CallingScratch<StutterSubstitutionScratch> = CallingScratch::default();
    arm.call_locus(
        &evidence,
        &parameters,
        selection.alleles().clone(),
        &config,
        &mut scratch,
    )
}

/// The genotype a sample was called, as allele ids — `[0, 0]` for the homozygous reference.
fn genotype_of(inference: &LocusInference, sample: usize) -> Vec<u16> {
    inference.per_sample[sample]
        .genotype()
        .unwrap_or_else(|| panic!("sample {sample} was called missing"))
        .alleles()
        .iter()
        .map(|allele| allele.get())
        .collect()
}

/// **What a repeat tract's call rested on**, which every tract carries and no SNP or indel
/// does — so calling this on an ordinary site's result is a test asking the wrong path.
fn tract_record(inference: &LocusInference) -> RepeatTractProvenance {
    inference
        .repeat_tract
        .expect("every called repeat tract carries what its parameters rested on")
}

/// Which rung of the tract ladder that record's prior shape came from.
fn rung_of(inference: &LocusInference) -> LengthSpectrumRung {
    tract_record(inference).length_spectrum_rung()
}

/// **ng calls genotypes at a SNP, and the three a reader can derive by hand come out.**
///
/// Three samples at one position whose reference base is `A`. The first showed 20 reads of `T`
/// and nothing else, the second 10 of each, the third 20 of `A`. At Phred 13 a read those are
/// not close calls, so the answers do not turn on the prior's strength: **`1/1`, `0/1`, `0/0`**.
///
/// **The candidates are selected, not supplied.** Selection asks its bar of each sample
/// separately rather than of a cohort pool, and sample 0's 20 `T` reads of 20 clear its own bar
/// — at least 2 reads and at least 10% of that sample's compared reads — which admits `T` for
/// the whole cohort. So the locus is called over two alleles and the diploid genotype table
/// holds three.
///
/// **The loop settles in one pass here**, asserted below, because at twenty reads a sample the
/// emission decides every genotype and the first pass's frequencies are already the answer. That
/// is why `at_three_reads_a_sample_the_loop_iterates_and_the_calls_hold` exists: without it
/// nothing in this file would notice a convergence rule that never ran.
#[test]
fn three_samples_at_a_snp_are_called_from_their_reads() {
    let per_sample = vec![
        sample_locus(vec![showed(b"T", 20)]),
        sample_locus(vec![showed(b"A", 10), showed(b"T", 10)]),
        sample_locus(vec![showed(b"A", 20)]),
    ];
    let inference = call(&per_sample);

    assert_eq!(inference.region, region());
    assert_eq!(
        inference.alleles().len(),
        2,
        "the reference and one alternative"
    );
    assert_eq!(inference.alleles().reference(), b"A");
    assert_eq!(inference.per_sample.len(), 3);

    assert_eq!(genotype_of(&inference, 0), vec![1, 1], "20 reads of T");
    assert_eq!(genotype_of(&inference, 1), vec![0, 1], "10 of each");
    assert_eq!(genotype_of(&inference, 2), vec![0, 0], "20 reads of A");

    assert!(inference.converged, "a locus this one-sided settles");
    assert_eq!(
        inference.passes, 1,
        "twenty reads a sample decide it on the first pass"
    );
}

/// **The cohort's expected allele copies are what the three calls carry**, to within what the
/// prior is still worth at twenty reads a sample.
///
/// Two chromosomes each from `1/1`, `0/1` and `0/0` put 3 copies of the alternative and 3 of the
/// reference among the cohort's six. Measured, the loop settles at **3.0000019 reference copies
/// against 2.9999981 alternative** — two parts in a million off an even split, which is what is
/// left of the prior's pull toward the reference once each sample has twenty reads.
#[test]
fn the_cohorts_expected_copies_are_the_three_calls_chromosomes() {
    let per_sample = vec![
        sample_locus(vec![showed(b"T", 20)]),
        sample_locus(vec![showed(b"A", 10), showed(b"T", 10)]),
        sample_locus(vec![showed(b"A", 20)]),
    ];
    let inference = call(&per_sample);

    let copies = inference.cohort_expected_copies().copies();
    assert_eq!(copies.len(), 2);
    assert!(
        (copies[0] - 3.000_001_9).abs() < 1e-6 && (copies[1] - 2.999_998_1).abs() < 1e-6,
        "three calls' six chromosomes, split as their reads split them: {copies:?}"
    );
    assert!(
        (copies[0] + copies[1] - 6.0).abs() < 1e-9,
        "every called sample contributes its ploidy: {copies:?}"
    );
}

/// **One sample is called on its own reads and nothing else** — the small end of the cohort-size
/// range this caller commits to (`CLAUDE.md`), where there is no panel to draw a frequency
/// from.
///
/// The single sample shows 20 reads of `T` and none of `A`, and comes out `1/1`. It is the
/// hardest case for the prior, not for the loop: with one sample the leave-one-out term is
/// empty, so the concentration is the seed alone.
#[test]
fn one_sample_is_called_on_its_own_reads() {
    let per_sample = vec![sample_locus(vec![showed(b"T", 20)])];
    let inference = call(&per_sample);

    assert_eq!(inference.per_sample.len(), 1);
    assert_eq!(genotype_of(&inference, 0), vec![1, 1]);
    assert!(inference.converged);
    assert_eq!(inference.passes, 1);
}

/// **A sample that covered nothing keeps its place in the run's sample order, and is called
/// heterozygous — which is not what a reader would guess, and the margin is two parts in ten
/// thousand.**
///
/// Three samples, of which the middle one has no reads at the locus at all. **So the merge's
/// covering-sample list holds two entries and the loop's holds three**, and the join between them
/// is the one no type enforces (`evidence_shaping`'s own doc): a defect there gives sample 2 the
/// call sample 1 should have had. Asserted below on both, and on the merge's list itself, because
/// the first draft of this fixture gave the silent sample a covering entry — and with every
/// sample covering, the two indices are the same number and the join cannot be tested at all.
///
/// **What the silent sample is called, measured: `0/1`.** Its reads score every genotype alike,
/// so what decides it is its own prior — the seed plus what the *other* samples showed, which the
/// leave-one-out subtraction is. Its two neighbours are called `1/1` and `0/0`, so they contribute
/// exactly `[2.0, 2.0]` copies, and against a seed of `[1.0, 0.001]` its concentration is
/// `[3.0000019, 2.000998]`. **The heterozygote wins by two parts in ten thousand** — posterior
/// `[0.39985, 0.40005, 0.20009]`, a genotype quality of about 2 Phred — so this is a knife-edge
/// rather than a confident call, and the assertion below is on the genotype the code produces
/// rather than on one anybody derived first.
///
/// **This is close to a question the owner has open and is not the same one.** C1's ⚑ is that a
/// silent sample's *flat pass* inflates the cohort estimate **other** samples see; here the
/// leave-one-out subtraction removes its own vote from its own prior, so what decides it is
/// genuinely its two neighbours. What the two share is the underlying fact — a sample with no
/// reads is not a sample with no answer — and at three reads a position roughly one sample in
/// twenty is silent at any given position.
#[test]
fn a_sample_that_covered_nothing_keeps_its_place_and_is_called_from_its_neighbours() {
    let per_sample = vec![
        sample_locus(vec![showed(b"T", 20)]),
        sample_locus(Vec::new()),
        sample_locus(vec![showed(b"A", 20)]),
    ];

    // **The join, before the loop ever runs**: two covering samples against three of the run.
    let merged = merge(&per_sample);
    assert_eq!(
        merged.per_sample.len(),
        2,
        "the merge lists only the samples that covered the locus"
    );
    assert_eq!(
        merged
            .per_sample
            .iter()
            .map(|support| support.sample)
            .collect::<Vec<_>>(),
        vec![0, 2],
        "and it names which samples of the run they are"
    );

    let inference = call(&per_sample);
    assert_eq!(inference.per_sample.len(), 3);
    assert!(
        !inference.per_sample[1].is_missing(),
        "a sample with no reads is called, not set aside — it is the candidate step that sets a \
         sample aside, and it ruled nobody out here"
    );
    assert_eq!(genotype_of(&inference, 0), vec![1, 1], "20 reads of T");
    assert_eq!(
        genotype_of(&inference, 1),
        vec![0, 1],
        "no reads, so its own prior decides — see this test's doc comment for the margin"
    );
    assert_eq!(genotype_of(&inference, 2), vec![0, 0], "20 reads of A");

    let copies = inference.cohort_expected_copies().copies();
    assert!(
        (copies[0] - 3.199_763_6).abs() < 1e-6,
        "the cohort after the silent sample's own posterior is folded in: {copies:?}"
    );
}

/// **The weakest warrant behind the locus reaches the output**, and here it is the one the run's
/// calibration carries.
///
/// The fixture's read group has a defaulted calibration — nothing was fitted for it — so the
/// locus's warrant is `Defaulted`, not `FittedHere`. That is the whole point of the field: a run
/// that trusted its instrument and one that measured it are otherwise indistinguishable in the
/// output (`doc/devel/ng/spec/read_likelihoods.md` §3.2).
#[test]
fn the_locus_carries_the_weakest_warrant_of_the_parameters_that_reached_it() {
    let per_sample = vec![
        sample_locus(vec![showed(b"T", 20)]),
        sample_locus(vec![showed(b"A", 20)]),
    ];
    let inference = call(&per_sample);
    assert_eq!(inference.weakest_provenance, Provenance::Defaulted);

    // And with the calibration fitted, the same locus claims the stronger warrant — so the
    // field is reporting what reached it rather than a constant.
    let merged = merge(&per_sample);
    let mut selection_scratch = SelectionScratch::new();
    let selection = select_generic(
        &merged,
        &CandidateSelectionConfig::DEFAULT,
        &mut selection_scratch,
    );
    let mut shaping = GenericEvidenceScratch::default();
    let mut views: Vec<GenericLocusSample<'_>> = Vec::new();
    let evidence = shape_generic_locus(
        &mut shaping,
        &merged,
        &selection,
        per_sample.len(),
        &mut views,
    );
    let calibration = vec![ReadGroupCalibration {
        scale: 1.0,
        provenance: Provenance::FittedHere,
    }];
    let inbreeding = vec![InbreedingF::try_new(0.0).expect("an outbred sample"); 2];
    let strata = StratumFits::over(&[], BTreeMap::new());
    let substitution = BTreeMap::new();
    let parameters = FrozenParameters::uncontaminated(
        &calibration,
        &inbreeding,
        SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
        &strata,
        &substitution,
        diploid(),
    );
    let arm = SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior);
    let config = CallingLoopConfig::DEFAULT
        .validate()
        .expect("the shipped configuration");
    let mut scratch: CallingScratch<StutterSubstitutionScratch> = CallingScratch::default();
    let fitted = arm.call_locus(
        &evidence,
        &parameters,
        selection.alleles().clone(),
        &config,
        &mut scratch,
    );
    assert_eq!(fitted.weakest_provenance, Provenance::FittedHere);
}

/// **A locus the merge builds and selection then calls over the reference alone.**
///
/// **The two bars are different and this is the gap between them**, which is the only reason
/// such a locus exists: the merge keeps a locus a sample showed at least 2 non-reference reads
/// and at least **2%** of its own at, and selection admits a sequence a sample showed at least 2
/// reads and at least **10%** of its own at. So one sample showing **3 `T` reads out of 100**
/// clears the merge (3 against 2) and fails selection (3 against 10). The candidate table holds
/// one allele, the genotype table one genotype, and both samples are called `0/0`.
///
/// **A locus of one allele is the narrow end of the loop's own range** — one genotype, whose
/// posterior is 1 whatever the prior — and nothing else in this file reaches it through
/// selection.
///
/// *(The first version of this test used one stray read of 21, which fails **both** bars: the
/// merge would have called that locus too quiet and never built it. `merge`'s keep-rule check is
/// what caught it.)*
#[test]
fn a_locus_the_merge_builds_can_still_be_called_over_the_reference_alone() {
    let per_sample = vec![
        sample_locus(vec![showed(b"A", 97), showed(b"T", 3)]),
        sample_locus(vec![showed(b"A", 100)]),
    ];
    let inference = call(&per_sample);

    assert_eq!(
        inference.alleles().len(),
        1,
        "3 reads in 100 clears the merge's 2% and not selection's 10%"
    );
    assert_eq!(genotype_of(&inference, 0), vec![0, 0]);
    assert_eq!(genotype_of(&inference, 1), vec![0, 0]);
    assert_eq!(inference.passes, 1, "one allele, one genotype, one pass");
}

/// **A three-allele locus is called over all three, and the candidate ids are not the merge's.**
///
/// Five samples over a reference `A`. The merge interns sequences in the order it meets them
/// walking samples in ascending order, so its table is **`[A, T, G, C]`**; `G` is one sample's 3
/// reads in 100, which clears the merge's 2% bar and fails selection's 10%, so it is dropped and
/// the candidates are **`[A, T, C]`**. **`C` is therefore merge index 3 and candidate id 2**, and
/// the assertions below are on the literal ids rather than on a lookup by bases — which is what
/// makes them able to fail.
///
/// *(The first version of this test had nothing dropped, so the remapping was the identity and
/// the ids were resolved by looking each sequence up in the output. Deleting the remapping
/// entirely left it green.)*
#[test]
fn a_locus_with_two_alternatives_is_called_over_all_three_and_the_ids_are_selections() {
    let per_sample = vec![
        sample_locus(vec![showed(b"T", 20)]),
        sample_locus(vec![showed(b"A", 97), showed(b"G", 3)]),
        sample_locus(vec![showed(b"C", 20)]),
        sample_locus(vec![showed(b"A", 10), showed(b"T", 10)]),
        sample_locus(vec![showed(b"A", 20)]),
    ];

    let merged = merge(&per_sample);
    assert_eq!(
        merged
            .alleles
            .iter()
            .map(|allele| allele.as_ref())
            .collect::<Vec<&[u8]>>(),
        vec![b"A".as_slice(), b"T", b"G", b"C"],
        "the merge interns in the order it meets each sequence"
    );

    let inference = call(&per_sample);
    assert_eq!(inference.alleles().len(), 3, "G is below selection's bar");
    let bases_of = |candidate: u16| {
        inference
            .alleles()
            .bases_of(AlleleId(candidate))
            .expect("an allele of this table")
            .to_vec()
    };
    assert_eq!(bases_of(0), b"A".to_vec(), "the reference is candidate 0");
    assert_eq!(bases_of(1), b"T".to_vec());
    assert_eq!(
        bases_of(2),
        b"C".to_vec(),
        "C is merge index 3 and candidate 2 — the remapping closed the gap G left"
    );

    assert_eq!(genotype_of(&inference, 0), vec![1, 1], "20 reads of T");
    assert_eq!(
        genotype_of(&inference, 1),
        vec![0, 0],
        "3 G in 100 is not a candidate"
    );
    assert_eq!(genotype_of(&inference, 2), vec![2, 2], "20 reads of C");
    assert_eq!(genotype_of(&inference, 3), vec![0, 1], "10 A and 10 T");
    assert_eq!(genotype_of(&inference, 4), vec![0, 0], "20 reads of A");
}

/// **The loop iterates at three reads a sample, which is where this project's own cohort sits.**
///
/// Every other fixture here is at twenty reads a sample, where the emission swamps the prior and
/// the frequency loop settles in a single pass — so nothing else in this file would notice a
/// convergence rule that never ran. **At three reads the same three genotypes come out and the
/// loop takes four passes**, which is the regime the tomato cohort is in (about three reads a
/// position, 63 accessions) and the one `CLAUDE.md` names as the hardest.
///
/// The cohort's expected copies land at **3.3044 reference against 2.6956 alternative** out of
/// six — further from an even split than the twenty-read fixture's two parts in a million,
/// because at three reads the prior is still worth something against the reads.
#[test]
fn at_three_reads_a_sample_the_loop_iterates_and_the_calls_hold() {
    let per_sample = vec![
        sample_locus(vec![showed(b"T", 3)]),
        sample_locus(vec![showed(b"A", 2), showed(b"T", 1)]),
        sample_locus(vec![showed(b"A", 3)]),
    ];
    let inference = call(&per_sample);

    assert_eq!(inference.alleles().len(), 2);
    assert_eq!(genotype_of(&inference, 0), vec![1, 1]);
    assert_eq!(genotype_of(&inference, 1), vec![0, 1]);
    assert_eq!(genotype_of(&inference, 2), vec![0, 0]);
    assert_eq!(
        inference.passes, 4,
        "at three reads the frequency loop has work to do"
    );
    assert!(inference.converged);

    let copies = inference.cohort_expected_copies().copies();
    assert!(
        (copies[0] - 3.304_434_8).abs() < 1e-6,
        "the prior still moves the cohort at three reads: {copies:?}"
    );
}

/// **The locus's warrant is the weakest of the read groups whose reads reached it**, and this is
/// the only test in the repository that combines two different ones.
///
/// One sample sequenced from two libraries — 10 reads of `T` from each — where the first
/// library's calibration was fitted and the second's was not. The locus's warrant is `Defaulted`,
/// because a call resting on one fitted parameter and one defaulted one is a defaulted call
/// (`spec/read_likelihoods.md` §4.4).
///
/// **Asserted in both orders**, because the fold is over a list: with the defaulted calibration
/// first and with it second. A fold that kept the last warrant instead of the weaker of the two
/// passes one order and fails the other — measured, and until this test existed it passed all
/// 4,815 of the library's tests as well.
#[test]
fn the_locus_takes_the_weakest_warrant_of_two_read_groups() {
    for defaulted_library in [0_u32, 1] {
        let per_sample = vec![
            sample_locus(vec![
                showed_from(b"T", 10, ReadGroupId(0)),
                showed_from(b"T", 10, ReadGroupId(1)),
            ]),
            sample_locus(vec![showed_from(b"A", 20, ReadGroupId(0))]),
        ];
        let calibration: Vec<ReadGroupCalibration> = (0..2)
            .map(|library| ReadGroupCalibration {
                scale: 1.0,
                provenance: if library == defaulted_library {
                    Provenance::Defaulted
                } else {
                    Provenance::FittedHere
                },
            })
            .collect();
        let inference = call_with_calibration(&per_sample, &calibration);
        assert_eq!(
            inference.weakest_provenance,
            Provenance::Defaulted,
            "with library {defaulted_library} defaulted, the locus is a defaulted call"
        );
        assert_eq!(genotype_of(&inference, 0), vec![1, 1]);
    }
}

/// **And with every read group fitted, the same locus claims the stronger warrant** — so the
/// field reports what reached it rather than a constant.
#[test]
fn a_locus_whose_read_groups_were_all_fitted_says_so() {
    let per_sample = vec![
        sample_locus(vec![
            showed_from(b"T", 10, ReadGroupId(0)),
            showed_from(b"T", 10, ReadGroupId(1)),
        ]),
        sample_locus(vec![showed_from(b"A", 20, ReadGroupId(0))]),
    ];
    let calibration = vec![
        ReadGroupCalibration {
            scale: 1.0,
            provenance: Provenance::FittedHere,
        };
        2
    ];
    let inference = call_with_calibration(&per_sample, &calibration);
    assert_eq!(inference.weakest_provenance, Provenance::FittedHere);
}

// ────────────────────────────────────────────────────────────────────────────────────────
// The repeat-tract path
// ────────────────────────────────────────────────────────────────────────────────────────
//
// **What is supplied here that is chosen above.** On the SNP/indel path the candidates come out
// of `select_generic`, run by this file. **At a repeat tract they are supplied by the fixture,
// and so are their repeat counts** — the repeat-tract half of candidate selection is not
// written, so there is no step to run. A reader must not take these candidate sets for a
// step's output: nothing here chose them.
//
// **Why the repeat counts travel at all.** How many whole repeats a candidate holds is not its
// byte length divided by the motif's: an interrupted tract — one whose repeat is broken by a
// substitution — holds fewer. The count is what picks the stratum whose slippage numbers score
// the candidate, so it has to be measured rather than derived. The locus generator measures it;
// until selection exists, a caller supplies it.
//
// Everything after the evidence is run: `shape_ssr_locus`, and the loop.

/// A dinucleotide `AT` tract, with the flanks the read model anchors its alignment on.
fn tract_detail() -> SsrDetail {
    SsrDetail {
        motif: Motif::new(b"AT").expect("a dinucleotide motif"),
        left_flank: Box::from(b"CCCGGG".as_slice()),
        right_flank: Box::from(b"TTTAAA".as_slice()),
    }
}

/// The two lengths every tract here is called over, in candidate-table order: **6 whole repeats,
/// which is the reference, and 7**.
const TRACT_CANDIDATE_REPEATS: [u32; 2] = [6, 7];

/// **How many libraries a tract fixture's run has, and it is deliberately not two.** The
/// parameter table a tract is scored from is `read groups × candidates`, filled
/// read-group-major; at an equal shape a table filled the other way round is the same length and
/// the same set of cells, so a transposition passes every shape check. Only the first library
/// sends a read.
const TRACT_READ_GROUPS: usize = 3;

fn tract_bases(repeats: u32) -> Vec<u8> {
    b"AT".repeat(repeats as usize)
}

/// The supplied candidate table — the reference tract and one longer allele.
fn tract_alleles() -> CandidateAlleles {
    let mut alleles = CandidateAlleles::new(
        tract_bases(TRACT_CANDIDATE_REPEATS[0]).into_boxed_slice(),
        LocusKind::Ssr(tract_detail()),
    );
    alleles.admit(tract_bases(TRACT_CANDIDATE_REPEATS[1]).into_boxed_slice());
    alleles
}

/// The supplied repeat counts, parallel to that table.
fn tract_repeat_counts() -> Vec<NonZeroU32> {
    TRACT_CANDIDATE_REPEATS
        .iter()
        .map(|count| NonZeroU32::new(*count).expect("a candidate always holds a repeat"))
        .collect()
}

/// `reads` reads that spanned the whole tract and showed `repeats` whole copies of it, from the
/// run's first library.
///
/// **The per-read error sum is filled in and nothing on this path reads it.** A tract's row is
/// scored from a stutter model and a per-base substitution rate; `q_sum` is the SNP/indel
/// emission's charge, and neither the repeat-tract row nor its emission model touches it.
fn tract_reads(repeats: u32, reads: u32) -> SequenceObservation {
    SequenceObservation {
        bases: tract_bases(repeats),
        read_witness: ReadWitness::Complete,
        read_group: ReadGroupId(0),
        num_obs: reads,
        num_fwd: reads / 2,
        q_sum: SummedLogError::from_nats(-10.0 * f64::from(reads)),
        mapq_sum: 60 * reads,
        mapq_sum_sq: u64::from(reads) * 3_600,
        placed_left: reads / 2,
        chain_ids: Vec::new(),
    }
}

/// **A repeat fit that reached both of this tract's candidate strata**, each with its own
/// slippage numbers and its own length spectrum.
///
/// The length spectrum is one share per whole-repeat offset from the reference tract length, two
/// repeats either way. **No two of its classes share a weight and neither spectrum is a
/// palindrome**, so reading one class for another, or reading a spectrum backwards, is a
/// different prior rather than the same one — and each leans toward contraction, as every real
/// repeat fit does.
///
/// **The run's three libraries do not share one slippage group.** Library 0 sits in group 0 and
/// libraries 1 and 2 in group 1, whose numbers differ, so a parameter table filled
/// candidate-major rather than read-group-major hands library 0's reads another group's
/// polymerase.
fn tract_strata() -> StratumFits {
    tract_strata_describing(BTreeMap::from([
        (ReadGroupId(0), 0),
        (ReadGroupId(1), 1),
        (ReadGroupId(2), 1),
    ]))
}

/// The same fit, over whichever of the run's libraries it claims to describe — what the fixture
/// about a library the fit never saw needs.
fn tract_strata_describing(
    slippage_group_of_each_library: BTreeMap<ReadGroupId, u32>,
) -> StratumFits {
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
    let fitted_stratum =
        |repeats: u64, level_in_group_0: f64, length_spectrum: Vec<f64>, concentration: f64| {
            StratumOutcome::Fitted(Box::new(StratumFit {
                stratum: FitStratum {
                    period: 2,
                    reference_repeats: repeats,
                },
                slippage: vec![
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
                ],
                length_spectrum,
                concentration,
                log_likelihood_a_tract: -1.5,
                tracts_fitted: 40,
                borrowed: Vec::new(),
                converged: true,
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
            // **No two adjacent pairs of these weights are in the same ratio.** The first
            // version's upper tail fell by a factor of three at each step, so a spectrum read
            // one repeat off centre gave the two candidates the same pair of shares and the
            // re-centring was invisible.
            fitted_stratum(6, 0.04, vec![0.10, 0.30, 0.44, 0.11, 0.05], 20.0),
            fitted_stratum(7, 0.06, vec![0.09, 0.21, 0.43, 0.19, 0.08], 25.0),
        ],
        slippage_group_of_each_library,
    )
}

/// **A fitted per-base substitution rate for every `(library, stratum)` cell the tract's
/// parameter table covers** — every library of the run against both candidate strata, so no cell
/// falls to the stated constant.
///
/// **Every library, not only the one whose reads arrive**, because the table is built on the
/// run's read-group axis and the locus's warrant is folded over all of it. The rate differs per
/// stratum, so a lookup keyed by the wrong repeat count is a different number.
fn tract_substitution_rates() -> BTreeMap<StratumKey, Estimate<ErrorRate>> {
    tract_substitution_rates_over(TRACT_READ_GROUPS)
}

/// The same, fitted for the first `libraries` of the run and for none of the rest.
fn tract_substitution_rates_over(libraries: usize) -> BTreeMap<StratumKey, Estimate<ErrorRate>> {
    let period = SsrPeriod::try_new(2).expect("a dinucleotide");
    let mut rates = BTreeMap::new();
    for library in 0..libraries {
        for repeats in TRACT_CANDIDATE_REPEATS {
            rates.insert(
                StratumKey {
                    read_group: ReadGroupId(u32::try_from(library).expect("a small index")),
                    stratum: SsrStratum::new(period, RepeatCount(repeats)),
                    ploidy: diploid(),
                },
                Estimate {
                    value: ErrorRate::try_new(0.001 * f64::from(repeats)).expect("a probability"),
                    provenance: Provenance::FittedHere,
                    observations: 4_000,
                },
            );
        }
    }
    rates
}

/// Run the tract path over one cohort's observations and hand back what came out.
fn call_tract(
    observations_of_each_sample: &[Vec<SequenceObservation>],
    strata: &StratumFits,
    substitution: &BTreeMap<StratumKey, Estimate<ErrorRate>>,
) -> LocusInference {
    let detail = tract_detail();
    let per_run_sample: Vec<&[SequenceObservation]> = observations_of_each_sample
        .iter()
        .map(Vec::as_slice)
        .collect();
    let repeat_counts = tract_repeat_counts();
    let mut views: Vec<SsrSampleEvidence<'_>> = Vec::new();
    let evidence = shape_ssr_locus(
        region(),
        &per_run_sample,
        &detail,
        &repeat_counts,
        &mut views,
    );

    let calibration = vec![ReadGroupCalibration::defaulted(); TRACT_READ_GROUPS];
    let inbreeding =
        vec![InbreedingF::try_new(0.0).expect("an outbred sample"); per_run_sample.len()];
    let parameters = FrozenParameters::uncontaminated(
        &calibration,
        &inbreeding,
        SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
        strata,
        substitution,
        diploid(),
    );
    let arm = SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior);
    let config = CallingLoopConfig::DEFAULT
        .validate()
        .expect("the shipped configuration");
    let mut scratch: CallingScratch<StutterSubstitutionScratch> = CallingScratch::default();
    arm.call_locus(
        &evidence,
        &parameters,
        tract_alleles(),
        &config,
        &mut scratch,
    )
}

/// The same, in a run whose parameter fit found `fraction` of each library's reads to have come
/// from another individual.
fn call_contaminated_tract(
    observations_of_each_sample: &[Vec<SequenceObservation>],
    fraction: f64,
) -> LocusInference {
    let detail = tract_detail();
    let per_run_sample: Vec<&[SequenceObservation]> = observations_of_each_sample
        .iter()
        .map(Vec::as_slice)
        .collect();
    let repeat_counts = tract_repeat_counts();
    let mut views: Vec<SsrSampleEvidence<'_>> = Vec::new();
    let evidence = shape_ssr_locus(
        region(),
        &per_run_sample,
        &detail,
        &repeat_counts,
        &mut views,
    );

    let calibration = vec![ReadGroupCalibration::defaulted(); TRACT_READ_GROUPS];
    let contamination = vec![
        ContaminationView {
            fraction,
            markers_with_reads: 400,
            reads_on_markers: 1_000,
            source: ContaminationSource::ThisReadGroupsReads,
        };
        TRACT_READ_GROUPS
    ];
    let inbreeding =
        vec![InbreedingF::try_new(0.0).expect("an outbred sample"); per_run_sample.len()];
    let batching = SequencingBatches::all_together_over(TRACT_READ_GROUPS, per_run_sample.len());
    let strata = tract_strata();
    let substitution = tract_substitution_rates();
    let parameters = FrozenParameters::new(
        &calibration,
        &contamination,
        &batching,
        &inbreeding,
        SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
        &strata,
        &substitution,
        diploid(),
    );
    let arm = SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior);
    let config = CallingLoopConfig::DEFAULT
        .validate()
        .expect("the shipped configuration");
    let mut scratch: CallingScratch<StutterSubstitutionScratch> = CallingScratch::default();
    arm.call_locus(
        &evidence,
        &parameters,
        tract_alleles(),
        &config,
        &mut scratch,
    )
}

/// **ng calls genotypes at a repeat tract** — reads in, genotypes out, over the run's own fitted
/// slippage numbers and its own fitted length spectrum.
///
/// Three samples at a dinucleotide `AT` tract called over 6 whole repeats — the reference — and
/// 7. The first showed 20 reads of the 7-repeat tract, the second 10 of each length, the third
/// 20 of the 6-repeat one. At a slippage level of 4 in 100, stutter cannot manufacture ten
/// reads at a wrong length, so the three come out **`1/1`, `0/1`, `0/0`**.
///
/// **The candidates are supplied, not selected** — see this section's own note above.
#[test]
fn three_samples_at_a_repeat_tract_are_called_from_their_reads() {
    let per_sample = vec![
        vec![tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)],
        vec![
            tract_reads(TRACT_CANDIDATE_REPEATS[0], 10),
            tract_reads(TRACT_CANDIDATE_REPEATS[1], 10),
        ],
        vec![tract_reads(TRACT_CANDIDATE_REPEATS[0], 20)],
    ];
    let inference = call_tract(&per_sample, &tract_strata(), &tract_substitution_rates());

    assert_eq!(inference.region, region());
    assert_eq!(inference.alleles().len(), 2, "the two supplied lengths");
    assert_eq!(inference.per_sample.len(), 3);

    assert_eq!(
        genotype_of(&inference, 0),
        vec![1, 1],
        "20 reads at 7 repeats"
    );
    assert_eq!(genotype_of(&inference, 1), vec![0, 1], "10 reads at each");
    assert_eq!(
        genotype_of(&inference, 2),
        vec![0, 0],
        "20 reads at 6 repeats"
    );
    assert!(inference.converged);

    assert_eq!(
        rung_of(&inference),
        LengthSpectrumRung::StratumsOwnFit,
        "the tract's prior came from its own stratum's fitted length spectrum"
    );
    assert_eq!(
        inference.weakest_provenance,
        Provenance::FittedHere,
        "every slippage number and every substitution rate this tract read was the fit's"
    );
}

/// **A run whose repeat fit reached nothing still calls its tracts**, from the bottom rung of
/// the ladder — and the record says so rather than claiming a measurement.
///
/// The same three samples and the same reads, against an empty fit. The stutter model falls to
/// the one HipSTR ships and the substitution rate to a stated constant, so the locus's warrant
/// is `Defaulted`; the prior's shape falls to a flat spectrum at a stated concentration, so the
/// rung is `StatedFlat`. **The genotypes are unchanged**, which is the point worth having: at
/// twenty reads a sample the reads decide, and what the fit buys is not this locus's calls but
/// the confidence attached to them.
///
/// **⚖ This fixture is the run the owner ruled on, 2026-08-27.**
/// `population_diversity.md` §5 used to want a tract in a run carrying no repeat-tract
/// parameters refused by name, where §4.4 wants the ladder to always answer; the two meet only
/// at a run that fitted nothing, which is this one. **It is called, and the rung says so**,
/// because refusing would turn a whole class of runs into a hard failure for a condition the
/// record already states.
#[test]
fn a_tract_in_a_run_that_fitted_nothing_is_still_called_and_says_what_it_rested_on() {
    let per_sample = vec![
        vec![tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)],
        vec![
            tract_reads(TRACT_CANDIDATE_REPEATS[0], 10),
            tract_reads(TRACT_CANDIDATE_REPEATS[1], 10),
        ],
        vec![tract_reads(TRACT_CANDIDATE_REPEATS[0], 20)],
    ];
    let inference = call_tract(
        &per_sample,
        &StratumFits::over(&[], BTreeMap::new()),
        &BTreeMap::new(),
    );

    assert_eq!(genotype_of(&inference, 0), vec![1, 1]);
    assert_eq!(genotype_of(&inference, 1), vec![0, 1]);
    assert_eq!(genotype_of(&inference, 2), vec![0, 0]);
    assert_eq!(
        rung_of(&inference),
        LengthSpectrumRung::StatedFlat,
        "the ladder always answers, and says from how far down"
    );
    assert_eq!(inference.weakest_provenance, Provenance::Defaulted);
}

/// **One sample's tract is called on its own reads** — the small end of the cohort-size range
/// this caller commits to (`CLAUDE.md`), where there is no panel to draw a length frequency
/// from.
///
/// **This is the case the prior's old construction refused outright.** A tract's prior belief
/// used to be built as a decay away from the cohort's commonest length and scaled to reproduce a
/// measured diversity, which one outbred sample cannot supply. It is now read from the fit's own
/// per-stratum length spectrum, which is fitted across *tracts* — and a single genome carries the
/// same tracts a panel does.
#[test]
fn one_sample_is_called_at_a_tract_on_its_own_reads() {
    let per_sample = vec![vec![tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)]];
    let inference = call_tract(&per_sample, &tract_strata(), &tract_substitution_rates());

    assert_eq!(inference.per_sample.len(), 1);
    assert_eq!(genotype_of(&inference, 0), vec![1, 1]);
    assert_eq!(
        rung_of(&inference),
        LengthSpectrumRung::StratumsOwnFit,
        "one sample reaches the same rung as a panel: the spectrum is fitted across tracts"
    );
}

/// **A tract carries no artifact summary**, where a SNP/indel locus carries one.
///
/// The artifact tests weigh strand and read-position imbalance between the reference and one
/// alternative allele. **At a tract what goes wrong is slippage**, which is already inside the
/// read likelihood, and a tract's own site quality is left to a document that is not written
/// (`doc/devel/ng/spec/calling_quality.md` §8). So the summary is absent rather than computed
/// over quantities that do not mean there what they mean at a SNP.
#[test]
fn a_repeat_tract_carries_no_artifact_summary_where_a_snp_carries_one() {
    let at_a_tract = call_tract(
        &[vec![tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)]],
        &tract_strata(),
        &tract_substitution_rates(),
    );
    assert!(at_a_tract.artifact_test_counts().is_none());

    let at_a_snp = call(&[sample_locus(vec![showed(b"T", 20)])]);
    assert!(at_a_snp.artifact_test_counts().is_some());
}

/// **A contaminated library's reads are not called as a second allele at a repeat tract** — the
/// failure the third term of the tract's read-likelihood mixture exists to prevent, end to end.
///
/// Three samples. The middle one carries two copies of the 6-repeat tract and shows twenty reads
/// of its own plus four at the 7-repeat length that came from another individual's DNA. **With
/// no fraction fitted it is called `0/1`**. The other two explanations cannot carry four reads:
/// slippage to exactly one repeat longer runs about **one read in two hundred** at this
/// stratum's fitted numbers, and the outlier term — reads no allele explains — is spread flat
/// over every length the tract can reach and is smaller still. **With the fraction the pre-pass
/// measured, 8 in 100, it is called `0/0`**, which is what it is.
///
/// **How this differs from the same correction at a SNP.** At an ordinary site the contaminating
/// population's frequency for the allele an observation shows is the cohort's own estimate, so it
/// moves as the loop iterates. At a tract there is no such number per length that is fixed before
/// calling — so the mixture uses the genotype prior's own belief about which lengths this tract
/// can be, which is the joint repeat fit's length spectrum for this tract's stratum. It is
/// specific to the locus, because it is indexed from this tract's own reference length, and it is
/// frozen, because the fit produced it before calling started.
#[test]
fn a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele() {
    let per_sample = vec![
        vec![tract_reads(TRACT_CANDIDATE_REPEATS[0], 20)],
        vec![
            tract_reads(TRACT_CANDIDATE_REPEATS[0], 20),
            tract_reads(TRACT_CANDIDATE_REPEATS[1], 4),
        ],
        vec![tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)],
    ];

    let clean = call_tract(&per_sample, &tract_strata(), &tract_substitution_rates());
    assert_eq!(
        genotype_of(&clean, 1),
        vec![0, 1],
        "with no fraction fitted, four reads at another length are a second allele"
    );

    let contaminated = call_contaminated_tract(&per_sample, 0.08);
    assert_eq!(
        genotype_of(&contaminated, 1),
        vec![0, 0],
        "with the fitted fraction they are somebody else's DNA"
    );

    // **The fraction's own value has to do work, not merely its existence.** At the same four
    // reads a fitted 5 in 100 is not enough mass to beat a heterozygote that must also account
    // for twenty reference reads, so it still calls `0/1`. A model that read `c` as a flag would
    // pass the two assertions above and fail this one.
    let barely = call_contaminated_tract(&per_sample, 0.05);
    assert_eq!(
        genotype_of(&barely, 1),
        vec![0, 1],
        "a smaller fitted fraction cannot explain the same four reads"
    );

    for sample in [0, 2] {
        assert_eq!(
            genotype_of(&clean, sample),
            genotype_of(&contaminated, sample),
            "sample {sample} is unambiguous at twenty reads and is called the same either way"
        );
    }
    // **The warrant does not cover the fraction, and this assertion does not claim it does.**
    // A tract's warrant is folded over the stutter and substitution parameters its row reads
    // per `(read group, candidate)`; the contamination fraction carries where it was fitted
    // from, and nothing folds that in. Reporting it is a separate job.
    assert_eq!(
        contaminated.weakest_provenance,
        Provenance::FittedHere,
        "every slippage number and every substitution rate this tract read was the fit's"
    );
}

/// **What a tract's call rested on comes out beside the genotypes, and a SNP's does not.**
///
/// This is the end-to-end half of what the run has to be able to say about itself. A genotype
/// scored under a fitted length spectrum and one scored under a stated flat shape are different
/// claims; so are one whose reads were shared out with a contaminant and one whose were not.
/// Nothing in the called genotype says which, so the record carries it.
///
/// **The numbers are worth writing out, and the middle run is arranged so that no two of them
/// are the same.** The run has three libraries and the tract two candidates, so the parameter
/// table is six cells — over every library of the run rather than the one whose reads arrived,
/// because that is the axis the read likelihood's context table is indexed on.
///
/// - **Fully fitted**: the fit names every library and every stratum these candidates reach, so
///   none of the six falls back.
/// - **Partly fitted**: the slippage fit names library 0 only and the substitution rates are
///   fitted for libraries 0 and 1, so the four numbers come out **6, 4, 4 and 2** — and the two
///   that coincide are the ones that must, since every one of library 1's and 2's cells is
///   defaulted for the same reason.
/// - **Fitted nowhere**: every cell falls back, by the reason that means the parameters and the
///   reads came from different runs.
#[test]
fn a_tract_says_what_its_parameters_rested_on_and_a_snp_says_nothing() {
    let per_sample = vec![
        vec![tract_reads(TRACT_CANDIDATE_REPEATS[0], 20)],
        vec![tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)],
    ];

    let clean = call_tract(&per_sample, &tract_strata(), &tract_substitution_rates());
    let record = tract_record(&clean);
    assert_eq!(
        record.scoring_cells(),
        TRACT_READ_GROUPS * 2,
        "three libraries over two candidates"
    );
    assert_eq!(record.cells_with_no_fitted_slippage(), 0);
    assert_eq!(record.cells_whose_read_group_the_fit_does_not_describe(), 0);
    assert_eq!(record.cells_with_no_fitted_substitution_rate(), 0);
    assert!(
        !record.contaminant_term_was_built(),
        "this run's fit found no contamination, so the mixture has two terms"
    );

    // **A fit that reaches part of the run, so that the four counts are four different-sized
    // answers rather than all-or-nothing.** Without this arm every count here is either 0 or the
    // cell total, and any two of the record's fields could be swapped unnoticed.
    let partly = call_tract(
        &per_sample,
        &tract_strata_describing(BTreeMap::from([(ReadGroupId(0), 0)])),
        &tract_substitution_rates_over(2),
    );
    let partial = tract_record(&partly);
    assert_eq!(partial.scoring_cells(), 6);
    assert_eq!(
        partial.cells_with_no_fitted_slippage(),
        4,
        "libraries 1 and 2 over both candidates"
    );
    assert_eq!(
        partial.cells_whose_read_group_the_fit_does_not_describe(),
        4,
        "all four, because both candidates' strata are fitted and only the libraries are not"
    );
    assert_eq!(
        partial.cells_with_no_fitted_substitution_rate(),
        2,
        "library 2's two cells; libraries 0 and 1 have a rate at both strata"
    );

    // **A run that fitted nothing falls back in every cell, and says so** — the same tract, the
    // same reads, over a fit that describes no library and no stratum.
    let unfitted = call_tract(
        &per_sample,
        &StratumFits::over(&[], BTreeMap::new()),
        &BTreeMap::new(),
    );
    let fell_back = tract_record(&unfitted);
    assert_eq!(fell_back.scoring_cells(), TRACT_READ_GROUPS * 2);
    assert_eq!(
        fell_back.cells_with_no_fitted_slippage(),
        TRACT_READ_GROUPS * 2
    );
    assert_eq!(
        fell_back.cells_whose_read_group_the_fit_does_not_describe(),
        TRACT_READ_GROUPS * 2,
        "a fit naming no library at all is the absence that means the parameters and the reads \
         came from different runs"
    );
    assert_eq!(
        fell_back.cells_with_no_fitted_substitution_rate(),
        TRACT_READ_GROUPS * 2
    );

    // **And the contaminant term is reported where it was built.**
    let contaminated = call_contaminated_tract(&per_sample, 0.08);
    assert!(tract_record(&contaminated).contaminant_term_was_built());

    // A SNP/indel locus carries no such record: its prior comes from the population's frequency
    // spectrum, whose ladder has different rungs, and it has no per-library stutter table at all.
    let snp = call(&[
        sample_locus(vec![showed(b"T", 20)]),
        sample_locus(vec![showed(b"A", 20)]),
    ]);
    assert_eq!(snp.repeat_tract, None);
}
