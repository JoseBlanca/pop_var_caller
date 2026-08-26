//! **ng calls genotypes** — the merge's output, through candidate selection and the input edge,
//! into the calling loop, with the genotypes derived by hand.
//!
//! Every earlier test of this loop hands it a likelihood table, an evidence view or a candidate
//! set built for the occasion. This one hands it a cohort locus and takes genotypes out, with
//! nothing between the two supplied by the test:
//!
//! ```text
//! per-sample observations → ClosedLocus → CohortObservation::over → select_generic
//!                        → shape_generic_locus → call_locus → LocusInference
//! ```
//!
//! # What this fixture supplies and what it runs, said exactly
//!
//! **Supplied**, and each for a reason:
//!
//! - **the per-sample observations** — what each sample's reads showed at the locus, as the
//!   SNP/indel locus generator emits them. Turning aligned reads into observations is step 5's
//!   and outside this plan (`calling_loop.md`'s Scope); a fixture that ran it would be testing
//!   that step here;
//! - **the `ClosedLocus`** the merge is handed. The chaining walk that groups overlapping
//!   observations into loci — `LocusCloser` — is not run. What [`merge`] does reproduce is the
//!   **keep rule**, `MinAltReads::DEFAULT` asked of each sample separately, so a fixture the
//!   real walk would have discarded as too quiet cannot pass here unnoticed;
//! - the run's frozen parameters and the loop's configuration, which are a run's inputs rather
//!   than a locus's.
//!
//! **Run rather than supplied**: the merge's allele unification and read attribution
//! (`CohortObservation::over`), candidate selection (`select_generic`), the input edge
//! (`shape_generic_locus`), and the loop (`call_locus`).
//!
//! # Why a test binary rather than a module inside the library
//!
//! It imports only what the crate exports, so it is also a check that this path's own seams are
//! `pub` — the twenty-odd items it names, not every seam in the middle of it.
//!
//! **The repeat-tract path is not here.** A tract needs its genotype prior's seed, which is
//! keyed by repeat count and takes two run-level numbers the parameter pre-pass does not emit
//! yet; that half is the plan's E3b.

use pop_var_caller::ng::calling::allele_candidates::generic::select_generic;
use pop_var_caller::ng::calling::allele_candidates::{CandidateSelectionConfig, SelectionScratch};
use pop_var_caller::ng::calling::evidence_shaping::{GenericEvidenceScratch, shape_generic_locus};
use pop_var_caller::ng::calling::genotype_prior::{
    MarginalizedDirichletPrior, SeedRegime, SpectrumSeed,
};
use pop_var_caller::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use pop_var_caller::ng::calling::inference::{CallingLoopConfig, LocusGenotyper};
use pop_var_caller::ng::calling::likelihood::ssr_emission::{
    StutterSubstitutionEmission, StutterSubstitutionScratch,
};
use pop_var_caller::ng::calling::{
    CallingScratch, FrozenParameters, GenericLocusSample, LocusInference, ReadGroupCalibration,
};
use pop_var_caller::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
};
use pop_var_caller::ng::parameter_estimation::Provenance;
use pop_var_caller::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use pop_var_caller::ng::parameter_estimation::ssr::StratumKey;
use pop_var_caller::ng::run::cohort_merge::build::CohortObservation;
use pop_var_caller::ng::run::cohort_merge::MinAltReads;
use pop_var_caller::ng::run::cohort_merge::close::{ClosedLocus, SampleMembers, Verdict};
use pop_var_caller::ng::types::{
    AlleleId, ContigId, ErrorRate, GenomeRegion, InbreedingF, Ploidy, Position, ReadGroupId,
};

use pop_var_caller::ng::parameter_estimation::Estimate;
use std::collections::BTreeMap;

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
        bases: bases.to_vec().into_boxed_slice(),
        read_witness: ReadWitness::Complete,
        read_group,
        num_obs: reads,
        num_fwd: reads / 2,
        q_sum: -3.0 * f64::from(reads),
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
        reference_bases: Box::from(b"A".as_slice()),
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
            .filter(|observation| observation.bases.as_ref() != b"A")
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

/// **One sample is called on its own reads and nothing else** — the cohort end of the range
/// this caller commits to (`CLAUDE.md`), where there is no panel to draw a frequency from.
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
    assert_eq!(genotype_of(&inference, 1), vec![0, 0], "3 G in 100 is not a candidate");
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
