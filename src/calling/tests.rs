use super::*;
use crate::calling::genotype_prior::SpectrumSeed;
use crate::locus_generation::SsrDetail;
use crate::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::types::{AlleleId, ContigId, InbreedingF, Motif, Ploidy, Position};

pub(super) fn generic_reference() -> CandidateAlleles {
    CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic)
}

/// A site quality standing in for one the worker computed, where what the test is
/// about is something else.
///
/// **Deliberately not zero.** Zero is a real answer — a locus at which nobody carries
/// anything — so a test whose fixture wrote zero would pass against an
/// implementation that lost the value on the way in.
fn a_worker_written_site_quality() -> Phred {
    Phred::try_new(37.0).expect("a legal quality, and below the site-quality ceiling")
}

/// A run whose repeat-tract substitution rates were never fitted — the empty map, which is
/// what a fixture calling no tract needs, and what
/// `FrozenParameters::ssr_substitution_rate_at` answers `None` from.
///
/// A `static` rather than a function, so that a call site can borrow it for as long as the
/// parameters live: `BTreeMap::new` is a `const fn`, and a temporary would be freed at the
/// end of the statement that built the view.
pub(super) static NO_SUBSTITUTION_RATES: std::collections::BTreeMap<
    crate::parameter_estimation::repeat_strata::StratumKey,
    crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
> = std::collections::BTreeMap::new();

/// **What a repeat tract's parameters rested on, for a fixture that is testing something
/// else** — a tract scored over six cells, one of which found no fitted slippage and none
/// of which was defaulted by a library the fit had never seen.
///
/// **The four counts are all different from one another**, so that the record is a different
/// value under any transposition — which matters not here, where `LocusInference::new` moves
/// the record whole and never reads inside it, but for anything that later does. **The
/// fixture that really discriminates the counts is the driver's**,
/// `a_tract_reports_how_many_of_its_cells_fell_back_and_why`, which reads them off a gather.
fn a_tract_record(rung: AlleleLengthVariationRung) -> RepeatTractProvenance {
    RepeatTractProvenance::new(rung, 6, 1, 0, 2, false)
}

fn diploid_call(first: u16, second: u16, quality: f32) -> SampleGenotypeCall {
    SampleGenotypeCall::Called {
        genotype: Genotype::new(vec![AlleleId(first), AlleleId(second)]),
        genotype_quality: Phred::try_new(quality).expect("a legal quality"),
        // These fixtures are about what a called locus holds, not about what the reads
        // said, so they use the case where the reads said something.
        reads_were_uninformative: false,
        genotype_rows: placeholder_genotype_rows(),
    }
}

pub(super) fn two_allele_locus() -> (CandidateAlleles, ExpectedAlleleCopies) {
    let mut alleles = generic_reference();
    alleles.admit(Box::from(b"T".as_slice()));
    let copies = ExpectedAlleleCopies::new(vec![2.6, 1.4], &alleles);
    (alleles, copies)
}

/// A repeat tract with two candidate lengths — the path the gene-diversity seed
/// marker belongs to, and the only one on which a locus may carry it.
pub(super) fn str_two_allele_locus() -> (CandidateAlleles, ExpectedAlleleCopies) {
    let mut alleles = CandidateAlleles::new(
        Box::from(b"ATAT".as_slice()),
        LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        }),
    );
    alleles.admit(Box::from(b"ATATAT".as_slice()));
    let copies = ExpectedAlleleCopies::new(vec![2.6, 1.4], &alleles);
    (alleles, copies)
}

pub(super) fn region() -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(3),
        start: Position(940),
        end: Position(940),
    }
}

/// A settled locus keeps everything a consumer needs and cannot reconstruct: the
/// alleles finally called over, one call per sample in run order, the cohort's
/// expected copies, and the four pieces of evidence for how the answer was reached.
#[test]
fn a_locus_carries_its_calls_and_the_evidence_for_how_they_were_reached() {
    let (alleles, copies) = two_allele_locus();
    let inference = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 0, 40.0), diploid_call(0, 1, 25.0)],
        copies,
        true,
        4,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );

    assert_eq!(inference.region, region());
    assert_eq!(inference.alleles().len(), 2);
    assert_eq!(inference.cohort_expected_copies().copies(), [2.6, 1.4]);
    assert_eq!(inference.passes, 4);
    assert!(inference.converged);
    assert_eq!(inference.weakest_provenance, Provenance::FittedHere);
    assert_eq!(inference.repeat_tract, None);

    // The per-sample calls are a sequence, in the run's sample order — the second
    // sample is the heterozygote, and swapping them would be a different record.
    assert_eq!(inference.per_sample.len(), 2);
    assert_eq!(
        inference.per_sample[0]
            .genotype()
            .expect("a called sample")
            .alleles(),
        [AlleleId(0), AlleleId(0)]
    );
    assert_eq!(
        inference.per_sample[1]
            .genotype()
            .expect("a called sample")
            .alleles(),
        [AlleleId(0), AlleleId(1)]
    );
    assert_eq!(
        inference.per_sample[1]
            .score_best_genotype()
            .expect("a called sample")
            .get(),
        25.0
    );
}

/// A locus that ran out of passes is emitted with the flag set, never dropped and
/// never fatal — one hard locus must not kill a cohort run. The flag is the whole
/// point: a genotype from a loop that did not settle is a weaker claim, and nothing
/// downstream can tell it from a settled one otherwise.
///
/// Built on the **repeat** path, because the tract ladder's rung this also carries
/// cannot arise on the SNP/indel one — a fixture that set it there would be pinning a
/// state its own field documents as impossible.
#[test]
fn a_locus_that_ran_out_of_passes_is_emitted_with_the_flag_set() {
    let (alleles, copies) = str_two_allele_locus();
    let capped = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 1, 3.0)],
        copies,
        false,
        50,
        Provenance::Defaulted,
        Some(a_tract_record(
            AlleleLengthVariationRung::PeriodsPooledTracts,
        )),
        a_worker_written_site_quality(),
        None,
    );

    assert!(!capped.converged, "the loop hit its cap");
    assert_eq!(capped.passes, 50);
    // The call is still there. Nothing about not converging removes it.
    assert_eq!(capped.per_sample.len(), 1);
    assert_eq!(
        capped.per_sample[0]
            .genotype()
            .expect("a called sample")
            .alleles(),
        [AlleleId(0), AlleleId(1)]
    );
    // And the two other warrants travel independently of it.
    assert_eq!(capped.weakest_provenance, Provenance::Defaulted);
    assert_eq!(
        capped.repeat_tract,
        Some(a_tract_record(
            AlleleLengthVariationRung::PeriodsPooledTracts
        ))
    );
}

/// One pass is a legitimate answer and must be constructible, because the loop's
/// first pass already has something to settle against: the reads-only estimate made
/// before iteration begins. A bound that drifted from `> 0` to `> 1` would reject
/// real single-pass loci at run time — and a cohort of one at three reads a position,
/// this caller's hardest committed case, is where a loop settles fastest.
#[test]
fn a_locus_that_settled_on_its_first_pass_is_a_locus() {
    let (alleles, copies) = two_allele_locus();
    let quick = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 0, 60.0)],
        copies,
        true,
        1,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );
    assert_eq!(quick.passes, 1);
    assert!(quick.converged);
}

/// Copies of a width the allele table does not have, in **both** directions.
///
/// Narrower is the obvious slip. **Wider is the one the pipeline produces**: the
/// final prune shrinks the allele table, so a copies vector that was not re-cut
/// alongside it is left long, and since every consumer indexes by `AlleleId` the
/// trailing entries would ride along unread — a pruned-away allele's expected copies
/// travelling in the output with nothing to notice.
#[test]
fn a_locus_cannot_carry_copies_of_a_width_its_alleles_do_not_have() {
    let build = |allele_count: usize, copy_count: usize| {
        let mut alleles = generic_reference();
        for _ in 1..allele_count {
            alleles.admit(Box::from(b"T".as_slice()));
        }
        let mut copies_alleles = generic_reference();
        for _ in 1..copy_count {
            copies_alleles.admit(Box::from(b"T".as_slice()));
        }
        let copies = ExpectedAlleleCopies::new(vec![1.0; copy_count], &copies_alleles);
        LocusInference::new(
            region(),
            alleles,
            vec![diploid_call(0, 0, 30.0)],
            copies,
            true,
            2,
            Provenance::FittedHere,
            None,
            a_worker_written_site_quality(),
            None,
        )
    };

    // Copies narrower than the alleles.
    let narrower = std::panic::catch_unwind(|| build(3, 2));
    assert!(narrower.is_err(), "3 alleles cannot carry 2 copy counts");

    // Copies wider than the alleles — what an un-re-cut vector looks like after the
    // prune, and the direction no fixture reached before.
    let wider = std::panic::catch_unwind(|| build(2, 3));
    assert!(wider.is_err(), "2 alleles cannot carry 3 copy counts");
}

/// Zero passes is not a locus that converged instantly — the loop always runs its
/// body at least once — it is a counter that was never incremented, and it would
/// distort the distribution the pass cap is to be set from.
#[test]
#[should_panic(expected = "at least one pass")]
fn a_locus_cannot_report_no_passes_at_all() {
    let (alleles, copies) = two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 0, 30.0)],
        copies,
        true,
        0,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );
}

/// A locus naming no call at all has lost them: a cohort has at least one sample and
/// every sample is called, so an empty list is a dropped result rather than a locus
/// with nothing to say.
#[test]
#[should_panic(expected = "has lost them")]
fn a_locus_cannot_name_no_sample_at_all() {
    let (alleles, copies) = two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles,
        Vec::new(),
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );
}

/// The tract ladder belongs to the repeat prior. A SNP/indel locus seeds from a
/// *frequency* spectrum, whose own ladder has different rungs, so a locus carrying a
/// tract record has had one path's ladder wired onto the other — which is exactly what an
/// implementation slip in the seed's routing would look like.
#[test]
#[should_panic(expected = "SNP/indel locus carries a repeat-tract record")]
fn a_snp_locus_cannot_carry_what_a_repeat_tract_rested_on() {
    let (alleles, copies) = two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 0, 30.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        Some(a_tract_record(
            AlleleLengthVariationRung::PeriodsPooledTracts,
        )),
        a_worker_written_site_quality(),
        None,
    );
}

/// **And the other direction, which is the newer half of the same check.** A tract whose
/// record went missing reports no rung and no defaulted-cell counts, which downstream reads
/// as *this call rested on nothing worth stating* rather than as a dropped field — so the
/// record is refused as absent at a tract exactly as it is refused as present at a SNP.
///
/// This direction is what a driver that built the record on the wrong branch would produce,
/// and the fixture is a tract in every other respect.
#[test]
#[should_panic(expected = "repeat-tract locus carries no repeat-tract record")]
fn a_repeat_tract_cannot_be_called_without_saying_what_it_rested_on() {
    let (alleles, copies) = str_two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 1, 30.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );
}

/// **A tract's record cannot claim more fallbacks than it has cells**, in either parameter.
///
/// The count and its denominator are read off one table, so a larger count is two tracts'
/// numbers in one record — and it would print as an answer rather than fail.
#[test]
#[should_panic(expected = "fall back on the slippage fit")]
fn a_tract_record_cannot_default_more_slippage_cells_than_it_has() {
    let _ =
        RepeatTractProvenance::new(AlleleLengthVariationRung::StratumsOwnFit, 6, 7, 0, 0, false);
}

/// The same on the other parameter, checked apart because a run of one candidate per library
/// cannot tell the two denominators from each other.
#[test]
#[should_panic(expected = "fall back on the substitution rate")]
fn a_tract_record_cannot_default_more_substitution_cells_than_it_has() {
    let _ =
        RepeatTractProvenance::new(AlleleLengthVariationRung::StratumsOwnFit, 6, 0, 0, 7, false);
}

/// **And the share cannot exceed what it is a share of.** The cells defaulted by a library
/// the fit never named are counted inside the arm that counts every defaulted cell, so a
/// larger one is the two counters swapped — the one mistake a reader of these numbers could
/// not detect, since both are plausible counts of the same table.
#[test]
#[should_panic(expected = "the two counters swapped")]
fn a_tract_record_cannot_blame_more_cells_on_an_unknown_library_than_it_defaulted() {
    let _ =
        RepeatTractProvenance::new(AlleleLengthVariationRung::StratumsOwnFit, 6, 2, 3, 0, false);
}

/// A region that runs backwards is not a locus. `GenomeRegion` holds no constructor
/// of its own and says a caller needing `start <= end` must say so; a called locus
/// needs it, because the region reaches the output's position column.
#[test]
#[should_panic(expected = "cannot run backwards")]
fn a_locus_cannot_cover_a_region_that_runs_backwards() {
    let (alleles, copies) = two_allele_locus();
    let _ = LocusInference::new(
        GenomeRegion {
            contig: ContigId(3),
            start: Position(1_100),
            end: Position(940),
        },
        alleles,
        vec![diploid_call(0, 0, 30.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );
}

/// A cohort of one is in scope, and is where this caller is weakest, so it must be
/// representable: one sample, one call, and the cohort's expected copies are that
/// sample's own, since it is the whole cohort.
#[test]
fn a_cohort_of_one_sample_is_a_locus_like_any_other() {
    let (alleles, copies) = two_allele_locus();
    let single = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 1, 12.0)],
        copies,
        true,
        2,
        Provenance::Borrowed,
        None,
        a_worker_written_site_quality(),
        None,
    );
    assert_eq!(single.per_sample.len(), 1);
    assert_eq!(single.cohort_expected_copies().copies(), [2.6, 1.4]);
}

/// The nine pooled counts as the final pass builds them at a biallelic locus whose
/// alternative drew 6 of the 8 reads — a summary of the right shape, for the two checks
/// this type makes on it.
fn artifact_counts_naming(primary_alternative: AlleleId) -> ArtifactTestCounts {
    ArtifactTestCounts {
        primary_alternative,
        reference_reads: 2.0,
        reference_forward_reads: 1.0,
        reference_placed_left_reads: 1.0,
        alternative_reads: 6.0,
        alternative_forward_reads: 5.0,
        alternative_placed_left_reads: 2.0,
        reference_low_quality_reads: 0.0,
        alternative_low_quality_reads: 0.0,
        total_reads: 8.0,
        genotype_expected_alternative_reads: 4.0,
    }
}

/// **Both of the quality plan's fields travel onto the record**, and the site quality is
/// read back through the one accessor there is.
///
/// `doc/devel/ng/spec/calling_quality.md` §3.5 makes that single accessor the design: the
/// worker writes the baseline, the artifact-correction stage overwrites the same field
/// with the corrected value, and **there is never a second quality field** for a gate and
/// a column to disagree about. Production shipped the other shape for sixteen days and
/// emitted 40 false positives at 30× on GIAB HG002 carrying a written `QUAL` of 0.
#[test]
fn the_site_quality_and_the_artifact_summary_travel_onto_the_record() {
    let (alleles, copies) = two_allele_locus();
    let inference = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 1, 25.0)],
        copies,
        true,
        3,
        Provenance::FittedHere,
        None,
        Phred::try_new(431.5).expect("a legal quality"),
        Some(artifact_counts_naming(AlleleId(1))),
    );

    assert_eq!(inference.uncorrected_site_quality().get(), 431.5);
    let counts = inference
        .artifact_test_counts()
        .expect("this locus carries a summary");
    assert_eq!(counts.primary_alternative, AlleleId(1));
    assert_eq!(counts.total_reads, 8.0);
}

/// A locus with no summary carries `None` rather than a summary of zeroes — the two are
/// different claims, and only one of them can be told from *the tests found no bias*.
#[test]
fn a_locus_with_nothing_for_the_artifact_tests_to_weigh_carries_no_summary() {
    let (alleles, copies) = two_allele_locus();
    let inference = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 0, 25.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );
    assert!(inference.artifact_test_counts().is_none());
}

/// A site quality above the ceiling did not come from the function that owns the ceiling.
///
/// `quality::score_uncorrected_site_quality` caps at `MAX_SITE_QUALITY`, which is also
/// its answer to a probability of exactly zero, so nothing it returns can be above it.
/// The field this is going into is the one the correction stage overwrites in place, and
/// a number from somewhere else in it is the shape §3.5 exists to prevent.
#[test]
#[should_panic(expected = "above the ceiling")]
fn a_site_quality_above_the_ceiling_is_refused() {
    let (alleles, copies) = two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 0, 25.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        Phred::try_new(quality::MAX_SITE_QUALITY + 1.0).expect("a legal Phred"),
        None,
    );
}

/// **The ceiling itself is accepted**, and the boundary is not academic: it is exactly
/// what `quality::score_uncorrected_site_quality` returns for a probability of no variant
/// of zero — a locus whose reads exclude the reference outright. A `<`/`<=` slip in the
/// check above would panic the worker on the loci the caller is most confident about.
#[test]
fn a_site_quality_at_the_ceiling_is_accepted() {
    let (alleles, copies) = two_allele_locus();
    let inference = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(1, 1, 25.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        Phred::try_new(quality::MAX_SITE_QUALITY).expect("a legal Phred"),
        None,
    );
    assert_eq!(
        inference.uncorrected_site_quality().get(),
        quality::MAX_SITE_QUALITY
    );
}

/// An artifact summary that names the **reference** as its primary alternative is
/// refused: both tests weigh one alternative against the reference, so a locus with
/// nothing to weigh carries no summary at all rather than one pointing at allele 0.
#[test]
#[should_panic(expected = "has no comparison to make")]
fn an_artifact_summary_naming_the_reference_as_its_alternative_is_refused() {
    let (alleles, copies) = two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 1, 25.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        Some(artifact_counts_naming(AlleleId::REFERENCE)),
    );
}

/// An artifact summary naming an allele this locus was not called over is refused, for
/// the reason a stale genotype is: the final prune renumbers every id above the one it
/// drops, so a summary not re-cut alongside the table points somewhere else.
#[test]
#[should_panic(expected = "pooled against a different allele table")]
fn an_artifact_summary_naming_an_allele_the_locus_lacks_is_refused() {
    let (alleles, copies) = two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles, // two alleles: ids 0 and 1
        vec![diploid_call(0, 1, 25.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        Some(artifact_counts_naming(AlleleId(4))),
    );
}

// ────────────────────────────────────────────────────────────────────────────────
// The three types the calling loop takes and gives back, and the missing genotype.
// ────────────────────────────────────────────────────────────────────────────────

use std::collections::BTreeMap;

pub(super) fn no_strata() -> StratumFits {
    StratumFits::over(&[], BTreeMap::new())
}

pub(super) fn diploid_ploidy() -> Ploidy {
    Ploidy::try_new(2).expect("a diploid")
}

pub(super) fn outbred_samples(samples: usize) -> Vec<InbreedingF> {
    vec![InbreedingF::try_new(0.0).expect("a legal coefficient"); samples]
}

pub(super) fn one_read_group() -> Vec<ReadGroupCalibration> {
    vec![ReadGroupCalibration::unscaled()]
}

pub(super) fn neutral_seed() -> SpectrumSeed {
    SpectrumSeed::new(1.0, 1e-3, genotype_prior::SeedRegime::NeutralShape)
}

/// The frozen parameters for a run of diploids, routed to whichever constructor the
/// contamination list calls for — an empty one is the run where the fit identified
/// none, which has its own named door.
///
/// **The contaminated branch takes the batching as an argument**, because it is borrowed
/// for as long as the parameters are and a helper that built one inside would hand back a
/// borrow of something already dropped. A run with no mixture has nothing to read it, so
/// the uncontaminated branch never looks at it.
pub(super) fn frozen_parameters<'a>(
    calibration: &'a [ReadGroupCalibration],
    contamination: &'a [ContaminationView],
    batching: &'a SequencingBatches,
    inbreeding: &'a [InbreedingF],
    strata: &'a StratumFits,
    substitution: &'a std::collections::BTreeMap<
        crate::parameter_estimation::repeat_strata::StratumKey,
        crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
    >,
) -> FrozenParameters<'a> {
    if contamination.is_empty() {
        FrozenParameters::uncontaminated(
            calibration,
            inbreeding,
            neutral_seed(),
            strata,
            substitution,
            diploid_ploidy(),
        )
    } else {
        FrozenParameters::new(
            calibration,
            contamination,
            batching,
            inbreeding,
            neutral_seed(),
            strata,
            substitution,
            diploid_ploidy(),
        )
    }
}

/// The default batching — one batch holding a run of `libraries` read groups over
/// `samples` samples, which is what a run that declared nothing gets.
pub(crate) fn one_batch(libraries: usize, samples: usize) -> SequencingBatches {
    let names: Vec<(String, String)> = (0..libraries)
        .map(|library| {
            (
                format!("rg{library}"),
                format!("s{}", library.min(samples.saturating_sub(1))),
            )
        })
        .collect();
    let borrowed: Vec<(&str, &str)> = names
        .iter()
        .map(|(id, sample)| (id.as_str(), sample.as_str()))
        .collect();
    let groups = crate::read::input::read_groups::ReadGroups::of_libraries(&borrowed);
    assert_eq!(
        groups.read_groups_per_sample().len(),
        samples,
        "the fixture asked for {samples} samples over {libraries} libraries"
    );
    SequencingBatches::all_together(&groups)
}

/// A missing genotype has no genotype and no quality — the two are absent together,
/// which is what an enum buys over a struct with an optional field.
///
/// **A quality beside a missing genotype would be the failure.** Emission has to write
/// this sample's `GT` as missing because no genotype was scored for it, not because a
/// scored one came out weak, and the two must not be conflated in the output.
#[test]
fn a_missing_call_carries_neither_a_genotype_nor_a_quality() {
    let missing = SampleGenotypeCall::Missing;
    assert!(missing.is_missing());
    assert_eq!(missing.genotype(), None);
    assert_eq!(missing.score_best_genotype(), None);

    let called = diploid_call(0, 1, 25.0);
    assert!(!called.is_missing());
    assert_eq!(
        called.genotype().expect("a called sample").alleles(),
        [AlleleId(0), AlleleId(1)]
    );
    assert_eq!(
        called.score_best_genotype().expect("a called sample").get(),
        25.0
    );
}

/// A locus can hand back a missing call beside called ones, in the run's sample order —
/// which is what a cohort whose cap cut one sample's earned allele produces.
#[test]
fn a_locus_carries_a_missing_call_beside_called_ones_in_run_order() {
    let (alleles, copies) = two_allele_locus();
    let inference = LocusInference::new(
        region(),
        alleles,
        vec![
            diploid_call(0, 0, 40.0),
            SampleGenotypeCall::Missing,
            diploid_call(0, 1, 25.0),
        ],
        copies,
        true,
        3,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );

    assert!(!inference.per_sample[0].is_missing());
    assert!(inference.per_sample[1].is_missing());
    assert!(!inference.per_sample[2].is_missing());
}

/// A repeat tract sets no sample aside, so a missing call at one is the SNP/indel path's
/// ruling wired onto the wrong path.
///
/// **The mirror of the tract-ladder rung's check**, which refuses the repeat path's
/// ruling at a SNP/indel locus. What this one catches is a sample silently dropped from
/// a tract's output: such a sample has no expected copies at all, so it also leaves the
/// cohort's expected-copies denominator, and the locus's allele frequencies come out
/// wrong with a well-formed record beside them.
#[test]
#[should_panic(expected = "a repeat tract sets no sample aside")]
fn a_repeat_tract_locus_cannot_carry_a_missing_call() {
    let (alleles, copies) = str_two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles,
        vec![diploid_call(0, 1, 20.0), SampleGenotypeCall::Missing],
        copies,
        true,
        2,
        Provenance::FittedHere,
        Some(a_tract_record(AlleleLengthVariationRung::StratumsOwnFit)),
        a_worker_written_site_quality(),
        None,
    );
}

/// A call naming an allele the locus was not called over is refused.
///
/// **The prune is what produces one.** Dropping an allele renumbers every id above it,
/// so a call minted before the prune and not remapped alongside it names an id past the
/// end — which reaches the `GT` column as an index past the `ALT` list, an unparseable
/// record rather than a crash.
#[test]
#[should_panic(expected = "names an allele this locus was not called over")]
fn a_call_naming_an_allele_the_locus_lost_is_refused() {
    let (alleles, copies) = two_allele_locus();
    let _ = LocusInference::new(
        region(),
        alleles, // two alleles: ids 0 and 1
        vec![diploid_call(0, 2, 20.0)],
        copies,
        true,
        2,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );
}
