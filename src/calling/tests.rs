use super::*;
use crate::locus_generation::SsrDetail;
use crate::types::{AlleleId, ContigId, Motif, Position};

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

use crate::parameter_estimation::joint::contamination::ContaminationSource;
use std::collections::BTreeMap;

pub(super) fn no_strata() -> StratumFits {
    StratumFits::over(&[], BTreeMap::new())
}

/// A run whose fit produced an allele length variation at one dinucleotide stratum and at one
/// trinucleotide stratum, with **different shapes and different concentrations** — so that a
/// lookup answering from the wrong stratum, or from the wrong period, is a different number
/// rather than the same one.
fn strata_with_allele_length_variations() -> StratumFits {
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
    let one = |period: u8, repeats: u64, allele_length_variation: Vec<f64>, concentration: f64| {
        StratumOutcome::Fitted(Box::new(StratumFit {
            standard_errors: None,
            stratum: Stratum {
                period,
                reference_repeats: repeats,
            },
            slippage: vec![Some(Slippage {
                level: 0.05,
                shorter_share: 0.83,
                fall_off: 0.25,
            })],
            allele_length_variation,
            concentration: Some(concentration),
            samples_fitted_on: None,
            tracts_of_its_own: 40,
            reads_crossing: 400,
            level_provenance: vec![Some(level)],
            shares_provenance: vec![Some(SharesProvenance {
                slipped_reads: Some(400.0),
                shorter_share: share,
                fall_off: share,
            })],
        }))
    };
    StratumFits::over(
        &[
            one(2, 10, vec![0.6, 0.3, 0.1], 4.0),
            one(3, 10, vec![0.1, 0.3, 0.6], 25.0),
        ],
        BTreeMap::from([(ReadGroupId(0), 0)]),
    )
}

/// **The run's frozen parameters answer a tract's prior shape from the tract's own stratum**,
/// and the accessor is keyed by the tract's reference repeat count where its two neighbours
/// on this type are keyed by the candidate's.
///
/// **Three fixtures' worth of coincidence is removed on purpose.** The two strata carry
/// different motif periods, different shapes and different concentrations, and the repeat
/// count asked for is not one either stratum would answer by accident — so a lookup off by
/// one repeat, or reading the wrong period, is a visibly different answer.
#[test]
fn the_run_answers_a_tracts_prior_shape_from_its_own_stratum() {
    use crate::parameter_estimation::joint::stratum_fits::AlleleLengthVariationRung;
    use crate::parameter_estimation::repeat_strata::RepeatCount;

    let strata = strata_with_allele_length_variations();
    let calibration = one_read_group();
    let inbreeding = outbred_samples(1);
    let batching = SequencingBatches::all_together_over(1, 1);
    let parameters = frozen_parameters(
        &calibration,
        &[],
        &batching,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
    );

    let dinucleotide = parameters.ssr_allele_length_variation_at(
        SsrPeriod::try_new(2).expect("a dinucleotide"),
        RepeatCount(10),
    );
    assert_eq!(
        dinucleotide.rung(),
        AlleleLengthVariationRung::StratumsOwnFit
    );
    assert_eq!(dinucleotide.fitted_weights(), Some(&[0.6, 0.3, 0.1][..]));
    assert_eq!(dinucleotide.concentration(), 4.0);

    // Same repeat count, other motif period: a different stratum and a different answer.
    let trinucleotide = parameters.ssr_allele_length_variation_at(
        SsrPeriod::try_new(3).expect("a trinucleotide"),
        RepeatCount(10),
    );
    assert_eq!(trinucleotide.fitted_weights(), Some(&[0.1, 0.3, 0.6][..]));
    assert_eq!(trinucleotide.concentration(), 25.0);

    // A repeat count neither stratum holds falls down the ladder rather than answering.
    let elsewhere = parameters.ssr_allele_length_variation_at(
        SsrPeriod::try_new(2).expect("a dinucleotide"),
        RepeatCount(11),
    );
    assert_eq!(elsewhere.rung(), AlleleLengthVariationRung::StatedFlat);
    // The median is taken on the between-locus share 1 / (1 + c): 4 and 25 have shares 1/5
    // and 1/26, whose mean is 31/260, a concentration of (1 − 31/260) / (31/260) = 229/31 —
    // not the mean of the two concentrations, 14.5.
    let expected = 229.0 / 31.0;
    assert!(
        (elsewhere.concentration() - expected).abs() < 1e-12,
        "the run fitted 4.0 and 25.0, whose median share converts back to 229/31 ≈ 7.39; \
         got {}",
        elsewhere.concentration()
    );
}

fn diploid_ploidy() -> Ploidy {
    Ploidy::try_new(2).expect("a diploid")
}

pub(super) fn outbred_samples(samples: usize) -> Vec<InbreedingF> {
    vec![InbreedingF::try_new(0.0).expect("a legal coefficient"); samples]
}

pub(super) fn one_read_group() -> Vec<ReadGroupCalibration> {
    vec![ReadGroupCalibration::unscaled()]
}

fn measured_contamination() -> ContaminationView {
    ContaminationView {
        fraction: 0.03,
        markers_with_reads: 400,
        reads_on_markers: 1_200,
        source: ContaminationSource::ThisReadGroupsReads,
    }
}

fn neutral_seed() -> SpectrumSeed {
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

/// A two-allele diploid locus: its allele table and the genotype table over it. Three
/// genotypes against two alleles, so a scratch that swapped the two counts is caught.
fn biallelic_locus() -> (CandidateAlleles, std::sync::Arc<GenotypeTable>) {
    let (alleles, _) = two_allele_locus();
    let table = GenotypeTable::build(diploid_ploidy(), alleles.len());
    (alleles, table)
}

/// Contamination is fitted per read group, so a run whose calibration covers one group
/// and whose contamination covers two has one of the two indexed by something else.
///
/// **The check has to be at construction, not at the first read.** The row looks a
/// fraction up by read group id, so a mismatched pair is found at whichever locus first
/// carries a read from the group past the end — or never, if no read ever comes from
/// one, in which case every genotype of the run is scored under somebody else's
/// contamination.
#[test]
#[should_panic(expected = "contamination is fitted per read group")]
fn frozen_parameters_refuse_a_contamination_list_of_another_read_group_count() {
    let calibration = one_read_group();
    let contamination = vec![measured_contamination(), measured_contamination()];
    let inbreeding = outbred_samples(2);
    let strata = no_strata();
    // A batching that matches neither list, so that what refuses this is the check the
    // test is about rather than a fixture that happened to disagree first.
    let _ = frozen_parameters(
        &calibration,
        &contamination,
        &one_batch(2, 2),
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
    );
}

/// No contamination anywhere has a **named** door, and it is a different claim from a
/// fitted zero: at one sample there is no panel to be surprised by, so the fraction is
/// not estimable at all rather than estimated and found clean.
#[test]
fn a_run_with_no_contamination_fitted_says_so_by_name() {
    let calibration = one_read_group();
    let inbreeding = outbred_samples(1);
    let strata = no_strata();
    let parameters = FrozenParameters::uncontaminated(
        &calibration,
        &inbreeding,
        neutral_seed(),
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid_ploidy(),
    );

    assert!(parameters.contamination_is_absent());
    assert!(parameters.contamination_by_read_group().is_empty());
    assert_eq!(parameters.read_group_count(), 1);
    assert_eq!(parameters.sample_count(), 1);
    assert_eq!(parameters.ploidy().get(), 2);
    assert_eq!(parameters.prior_seed().alpha_ref(), 1.0);
    assert_eq!(parameters.ssr_slippage_fits().strata(), 0);
}

/// The shortest thing that compiles is refused, so a caller reaches the decision.
///
/// **The two spellings are not equivalent and that is the whole point.** A consumer
/// handed an empty list writes `contamination.get(group).map(|v| v.fraction)
/// .unwrap_or(0.0)` and has silently turned *not estimable* into *estimated and found
/// clean* — which are different claims about the sample. The sibling half of the same
/// mixture, `ContaminationMixture::new`, refuses the empty spelling for this reason,
/// and this is the frozen half spelled to match.
#[test]
#[should_panic(expected = "is spelled `FrozenParameters::uncontaminated`")]
fn an_empty_contamination_list_is_refused_in_favour_of_the_named_constructor() {
    let calibration = one_read_group();
    let inbreeding = outbred_samples(1);
    let strata = no_strata();
    let no_declared_batching = one_batch(1, 1);
    let _ = FrozenParameters::new(
        &calibration,
        &[],
        &no_declared_batching,
        &inbreeding,
        neutral_seed(),
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid_ploidy(),
    );
}

/// The two axes the type calls "not interchangeable" are told apart by a fixture where
/// they differ: one library, three samples.
///
/// **Every fixture with one of each passes a swapped implementation.** `sample_count()`
/// is what the evidence check compares against, so `read_group_count()` and
/// `sample_count()` reading each other's field turns the run-order guard into a
/// read-group-count guard — invisible on a panel of 63 accessions with one library
/// each, and wrong on any run where a sample is sequenced twice.
#[test]
fn read_group_count_and_sample_count_are_two_different_axes() {
    let calibration = one_read_group();
    let inbreeding = outbred_samples(3);
    let strata = no_strata();
    let no_declared_batching = one_batch(1, 1);
    let parameters = frozen_parameters(
        &calibration,
        &[],
        &no_declared_batching,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
    );

    assert_eq!(parameters.read_group_count(), 1);
    assert_eq!(parameters.sample_count(), 3);
    assert_eq!(parameters.calibration_by_read_group().len(), 1);
    assert_eq!(parameters.inbreeding_coefficient_by_sample().len(), 3);
}

/// A run whose read-group axis went missing is refused at construction, rather than
/// producing an out-of-range index at whichever locus first carries a read.
#[test]
#[should_panic(expected = "read-group axis went missing")]
fn frozen_parameters_refuse_a_run_with_no_read_groups() {
    let inbreeding = outbred_samples(2);
    let strata = no_strata();
    let no_declared_batching = one_batch(1, 1);
    let _ = frozen_parameters(
        &[],
        &[],
        &no_declared_batching,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
    );
}

/// **A batching minted over a different set of libraries is refused where the run is
/// assembled**, not at whichever locus first carries a read from the library it misses.
///
/// The batching is declared by the user and the parameters are fitted, so the two can
/// disagree about how many libraries a run has without either being wrong on its own — and
/// the symptom without this is a mixture scoring some library against the neighbours of
/// another.
#[test]
#[should_panic(expected = "the batching covers 2 read groups")]
fn a_batching_over_another_runs_read_groups_is_refused() {
    let calibration = one_read_group();
    let contamination = vec![measured_contamination()];
    let inbreeding = outbred_samples(1);
    let strata = no_strata();
    let two_libraries = one_batch(2, 2);
    let _ = FrozenParameters::new(
        &calibration,
        &contamination,
        &two_libraries,
        &inbreeding,
        neutral_seed(),
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid_ploidy(),
    );
}

/// **And the other axis**, which is the one a run of one library per sample cannot tell
/// apart from the first: the sample-keyed batching is read by the run's own sample index, so
/// a shorter one drops the last samples out of every batch's copies.
#[test]
#[should_panic(expected = "the batching covers 1 samples")]
fn a_batching_over_another_runs_samples_is_refused() {
    let calibration = one_read_group();
    let contamination = vec![measured_contamination()];
    let inbreeding = outbred_samples(2);
    let strata = no_strata();
    let one_library_one_sample = one_batch(1, 1);
    let _ = FrozenParameters::new(
        &calibration,
        &contamination,
        &one_library_one_sample,
        &inbreeding,
        neutral_seed(),
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid_ploidy(),
    );
}

/// A run whose sample order went missing is refused, rather than producing a locus with
/// no calls in it.
#[test]
#[should_panic(expected = "an empty list is a run whose sample order went missing")]
fn frozen_parameters_refuse_a_run_with_no_samples() {
    let calibration = one_read_group();
    let strata = no_strata();
    let no_declared_batching = one_batch(1, 1);
    let _ = frozen_parameters(
        &calibration,
        &[],
        &no_declared_batching,
        &[],
        &strata,
        &NO_SUBSTITUTION_RATES,
    );
}

/// Each sample reads and writes its own row of the two flat tables, and the rows do not
/// overlap.
///
/// **The flat tables are where a slip is silent.** A caller slicing `lg_table` itself
/// with the allele count instead of the genotype count would read a window straddling
/// two samples' rows: every entry a real log-likelihood, none of them this sample's.
/// The fixture writes each sample's index into its own row so a wrong window names the
/// sample it actually came from.
#[test]
fn each_sample_reads_and_writes_its_own_row_of_the_scratch_tables() {
    let (alleles, table) = biallelic_locus();
    let view = table.view();
    assert_eq!(view.genotype_count(), 3);
    assert_eq!(view.allele_count(), 2);

    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(3, &alleles, &view);
    assert_eq!(scratch.row_count(), 3);
    assert_eq!(scratch.genotype_count(), 3);
    assert_eq!(scratch.allele_count(), 2);

    for sample in 0..3 {
        for (genotype, slot) in scratch
            .sample_genotype_likelihoods_mut(sample)
            .iter_mut()
            .enumerate()
        {
            *slot = LogProb(-(sample as f64) - 0.1 * genotype as f64);
        }
        for (allele, slot) in scratch
            .sample_expected_copies_mut(sample)
            .iter_mut()
            .enumerate()
        {
            *slot = sample as f64 + 0.5 * allele as f64;
        }
    }

    for sample in 0..3 {
        let row = scratch.sample_genotype_likelihoods(sample);
        assert_eq!(row.len(), 3);
        assert_eq!(row[0].get(), -(sample as f64));
        assert_eq!(row[2].get(), -(sample as f64) - 0.2);

        let copies = scratch.sample_expected_copies(sample);
        assert_eq!(copies.len(), 2);
        assert_eq!(copies, [sample as f64, sample as f64 + 0.5]);
    }
}

/// **[`CallingScratch::sample_scoring_buffers_mut`] hands back that sample's two rows**
/// of the flat sample-major tables — its likelihood row and its expected-copies row —
/// and not another sample's.
///
/// **The scoring path's own tests cannot see this**, which is why the check lives here.
/// Measured: with the method's two `sample` arguments replaced by `0`, so every sample
/// is handed sample 0's rows, the whole of `calling` still passes. The test that
/// looks like the one to catch it — the three-sample scoring test, whose `assert_ne!`
/// is commented *"sample 1 was scored on sample 0's likelihood row"* — passes under
/// that mutation, because it also redirects the expected-copies row, which the scorer
/// **overwrites** on every call; so each sample gets a different leave-one-out term and
/// therefore a different posterior, from the same reads.
///
/// The two rows are 3 and 2 entries wide here, so a wrong row is a legal slice of the
/// right table — the hazard [`row_range`](CallingScratch::row_range)'s own doc names.
#[test]
fn the_scoring_buffers_hand_back_that_samples_own_rows() {
    let (alleles, table) = biallelic_locus();
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(3, &alleles, &view);
    for sample in 0..3 {
        for (genotype, slot) in scratch
            .sample_genotype_likelihoods_mut(sample)
            .iter_mut()
            .enumerate()
        {
            *slot = LogProb(10.0 * sample as f64 + genotype as f64);
        }
        scratch
            .sample_expected_copies_mut(sample)
            .copy_from_slice(&[sample as f64, 2.0 - sample as f64]);
    }

    let buffers = scratch.sample_scoring_buffers_mut(1);
    assert_eq!(buffers.sample, 1);
    assert_eq!(
        buffers.genotype_likelihoods,
        &[LogProb(10.0), LogProb(11.0), LogProb(12.0)],
        "sample 1's likelihood row"
    );
    assert_eq!(
        buffers.sample_expected_copies,
        &[1.0, 1.0],
        "sample 1's own expected copies"
    );
}

/// **The buffers a repeat tract's rows are scored from refuse a scratch whose rows were
/// never claimed**, and refuse it at the door rather than inside the walk.
///
/// The walk over a tract's rows runs inside that one borrow, where no accessor of this
/// scratch can be called — so the map from a row back to the run's sample has to be checked
/// on the way in. Without the check, a scratch prepared for two rows and claimed for one
/// reads row 1 of a map one entry long, which is an index panic naming a slice where the
/// reader needs the cohort.
#[test]
#[should_panic(expected = "of them were claimed")]
fn the_tracts_locus_buffers_refuse_a_scratch_whose_rows_were_never_claimed() {
    let (alleles, table) = biallelic_locus();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &table.view());
    scratch.claim_row_for(0, InbreedingF::try_new(0.0).expect("an outbred sample"));
    let _ = scratch.tract_locus_buffers_mut();
}

/// A sample past the count the scratch was prepared for is refused by name, rather than
/// reading the next sample's row.
#[test]
#[should_panic(expected = "sample 2 is past the 2 this scratch was prepared for")]
fn a_sample_past_the_prepared_count_is_refused() {
    let (alleles, table) = biallelic_locus();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &table.view());
    let _ = scratch.sample_genotype_likelihoods(2);
}

/// **Preparing a locus overwrites every entry, even when the shape has not changed.**
///
/// This is the failure `Vec::resize` alone would leave: at two loci of the same shape it
/// keeps the leading entries, so the second locus would be scored against the first
/// one's likelihoods and priors with nothing failing. Every buffer comes back `NaN`, so
/// a value that survived is one that reaches an arithmetic check rather than a genotype.
#[test]
fn preparing_a_locus_overwrites_the_previous_locus_of_the_same_shape() {
    let (alleles, table) = biallelic_locus();
    let view = table.view();

    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &view);
    scratch.sample_genotype_likelihoods_mut(0)[0] = LogProb(-7.0);
    scratch.sample_genotype_likelihoods_mut(1)[2] = LogProb(-9.0);
    scratch.sample_expected_copies_mut(0)[1] = 1.75;
    scratch.cohort_expected_copies_mut()[0] = 3.5;
    scratch.seed_concentration_mut()[0] = 1.0;

    scratch.prepare_for_locus(2, &alleles, &view);

    assert!(
        scratch
            .sample_genotype_likelihoods(0)
            .iter()
            .all(|entry| entry.get().is_nan())
    );
    assert!(
        scratch
            .sample_genotype_likelihoods(1)
            .iter()
            .all(|entry| entry.get().is_nan())
    );
    assert!(
        scratch
            .sample_expected_copies(0)
            .iter()
            .all(|copy| copy.is_nan())
    );
    assert!(
        scratch
            .cohort_expected_copies()
            .iter()
            .all(|copy| copy.is_nan())
    );
    assert!(scratch.seed_concentration().iter().all(|a| a.is_nan()));
}

/// Advancing hands the previous pass's copies to the convergence test and gives back a
/// buffer for the next pass — a swap, so no allele-length copy is paid per pass.
#[test]
fn advancing_makes_the_current_copies_the_previous_ones() {
    let (alleles, table) = biallelic_locus();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &table.view());

    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[1.25, 0.75]);
    let next = scratch.advance_cohort_expected_copies();
    // Handed back unwritten, not holding the pass before last's real numbers: a pass
    // that skipped an allele would otherwise leave a value no pass wrote for the
    // convergence test to compare.
    assert!(next.iter().all(|copy| copy.is_nan()), "{next:?}");
    next.copy_from_slice(&[1.5, 0.5]);

    assert_eq!(scratch.previous_cohort_expected_copies(), [1.25, 0.75]);
    assert_eq!(scratch.cohort_expected_copies(), [1.5, 0.5]);
}

/// A locus prepared for no rows at all is a run whose sample order went missing. **Rows
/// rather than samples**: the scratch is sized for the samples the locus is *called on*,
/// so zero here is a locus nobody can be called at rather than a cohort of nobody.
#[test]
#[should_panic(expected = "a locus prepared for no rows")]
fn a_scratch_prepared_for_no_rows_is_refused() {
    let (alleles, table) = biallelic_locus();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(0, &alleles, &table.view());
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

/// A scratch that was never sized for a locus is refused by name, rather than handing a
/// pass an empty buffer to fold nothing into.
///
/// **The unwritten-value fill cannot cover this door**, and that is why the check
/// exists: an unprepared buffer has no slots to fill, so a fold over it runs zero
/// iterations, writes nothing, and leaves the cohort's expected copies summing to a
/// plausible `0.0` — which is exactly the value the fill exists to keep out of a
/// genotype. Step B2's oracle is a bitwise comparison of those sums, and an all-zero
/// sum passes against another all-zero sum.
#[test]
#[should_panic(expected = "has not been prepared for a locus")]
fn an_unprepared_scratch_is_refused() {
    let mut scratch = CallingScratch::<()>::default();
    for slot in scratch.cohort_expected_copies_mut() {
        *slot = 99.0;
    }
}

/// A genotype table built for a different allele set is refused where the locus's shape
/// is fixed, rather than at the far end of the locus.
///
/// **This is the first of the three caller bugs the spec names as assertions**, and a
/// discovery round admitting an allele is exactly how the two come apart. With a table
/// one allele narrow, every per-allele buffer is sized for the old count, and the first
/// thing to notice is the copies-versus-alleles check when the locus is finished —
/// thousands of arithmetic operations later, naming the output rather than the table.
#[test]
#[should_panic(expected = "the table was built for a different allele set")]
fn a_genotype_table_built_for_another_allele_set_is_refused() {
    let (mut alleles, table) = biallelic_locus();
    // A discovery round admits a third allele; the table is still the two-allele one.
    alleles.admit(Box::from(b"G".as_slice()));

    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &table.view());
}

/// The previous pass's copies are poisoned by the next locus too, not only the current
/// pass's.
///
/// **Without this, locus *n*'s first convergence comparison runs against locus *n−1*'s
/// final allele copies** — a locus declared settled on the previous locus's numbers,
/// with nothing failing. The fixture writes the previous-pass buffer between the two
/// loci, which is the one state the same-shape test does not reach.
#[test]
fn preparing_a_locus_poisons_the_previous_passs_copies_as_well() {
    let (alleles, table) = biallelic_locus();
    let view = table.view();

    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &view);
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[1.25, 0.75]);
    scratch
        .advance_cohort_expected_copies()
        .copy_from_slice(&[1.5, 0.5]);
    assert_eq!(scratch.previous_cohort_expected_copies(), [1.25, 0.75]);

    // The next locus, same shape.
    scratch.prepare_for_locus(2, &alleles, &view);
    assert!(
        scratch
            .previous_cohort_expected_copies()
            .iter()
            .all(|copy| copy.is_nan()),
        "the previous pass's copies survived into the next locus: {:?}",
        scratch.previous_cohort_expected_copies()
    );
}

/// The four buffers no accessor reached are sized on the axis they belong to.
///
/// **Three genotypes against two alleles is what makes this discriminating.** A buffer
/// sized on the wrong axis is still a legal length, so at a locus where the two counts
/// were equal every wrong sizing would pass; here the prior and posterior rows must be
/// 3 and the two per-allele buffers must be 2.
#[test]
fn the_prior_and_posterior_buffers_are_sized_on_their_own_axes() {
    let (alleles, table) = biallelic_locus();
    let view = table.view();

    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &view);

    assert_eq!(scratch.prior_row().len(), view.genotype_count());
    assert_eq!(scratch.posterior_row().len(), view.genotype_count());
    assert_eq!(scratch.sample_concentration().len(), view.allele_count());
    assert_eq!(
        scratch.prior_per_allele_workspace().len(),
        view.allele_count()
    );
    assert_ne!(view.genotype_count(), view.allele_count());

    // And they arrive unwritten, like every other buffer.
    assert!(scratch.prior_row().iter().all(|entry| entry.get().is_nan()));
    assert!(scratch.posterior_row().iter().all(|entry| entry.is_nan()));
    assert!(scratch.sample_concentration().iter().all(|a| a.is_nan()));
    assert!(
        scratch
            .prior_per_allele_workspace()
            .iter()
            .all(|a| a.is_nan())
    );

    // Written, then overwritten by the next locus of the same shape.
    scratch.prior_row_mut()[0] = LogProb(-1.5);
    scratch.posterior_row_mut()[0] = 0.25;
    scratch.sample_concentration_mut()[0] = 1.0;
    scratch.prior_per_allele_workspace_mut()[0] = 2.0;
    scratch.prepare_for_locus(2, &alleles, &view);
    assert!(scratch.prior_row()[0].get().is_nan());
    assert!(scratch.posterior_row()[0].is_nan());
    assert!(scratch.sample_concentration()[0].is_nan());
    assert!(scratch.prior_per_allele_workspace()[0].is_nan());
}

/// The configuration a run actually uses — the shipped repeat-tract emission model's
/// own scratch — is built here and nowhere else in the module's tests.
///
/// **Every other fixture uses `()` as the model scratch**, which compiles whether or
/// not the shipped model's `Scratch` still satisfies what the type needs. This also
/// gives the three sub-scratch accessors their only call site, and pins that a scratch
/// can cross a thread boundary, which "allocated once per worker" requires and nothing
/// else checked.
#[test]
fn the_shipped_emission_scratch_builds_and_can_cross_a_worker_boundary() {
    fn assert_send<T: Send>() {}
    assert_send::<CallingScratch<likelihood::ssr_emission::StutterSubstitutionScratch>>();

    let (alleles, table) = biallelic_locus();
    let mut scratch =
        CallingScratch::<likelihood::ssr_emission::StutterSubstitutionScratch>::default();
    scratch.prepare_for_locus(2, &alleles, &table.view());
    assert_eq!(scratch.row_count(), 2);

    // The three sub-scratches are reachable, and are the worker's rather than a locus's.
    let _ = scratch.candidate_selection_mut();
    let _ = scratch.generic_row_mut(0);
    let _ = scratch.ssr_row_mut();
}

/// **The contaminant tables refuse to be read before they are sized**, and that refusal is
/// what the sizing's own separateness costs.
///
/// `prepare_for_locus` un-sizes them and only a run that fitted a fraction re-sizes them, so
/// a driver that forgot the second call would otherwise sum every batch's copies into a
/// zero-length table: the fill writes nothing, returns having done nothing wrong, and the
/// mixture is built over an empty slice.
#[test]
#[should_panic(expected = "contaminant tables were not prepared")]
fn the_contaminant_tables_refuse_to_be_read_before_they_are_sized() {
    let (alleles, table) = biallelic_locus();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &table.view());
    scratch.claim_row_for(0, InbreedingF::try_new(0.0).expect("an outbred sample"));
    let _ = scratch.batch_copy_buffers_mut();
}

/// A locus prepared for no sequencing batches is a run whose batching went missing rather
/// than a run without one — the default is one batch holding all of it.
#[test]
#[should_panic(expected = "at least one sequencing batch")]
fn contaminant_tables_for_no_batches_are_refused() {
    let (alleles, table) = biallelic_locus();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &table.view());
    scratch.claim_row_for(0, InbreedingF::try_new(0.0).expect("an outbred sample"));
    scratch.prepare_contaminant_tables(0, 1);
}

/// **A row naming a sample the run does not have is not this run's locus**, and the copies
/// are scattered onto the run's axis by that name — so it would put one sample's copies
/// outside the table they are indexed in.
///
/// **The row count is not what decides it.** A one-row locus whose row names sample 5 of a
/// two-sample run has fewer rows than the run has samples and is still wrong, which is what
/// this fixture is: one row, a run of two, and the row names sample 5.
#[test]
#[should_panic(expected = "rows name sample 5 and the run has 2")]
fn a_row_naming_a_sample_the_run_does_not_have_is_refused() {
    let (alleles, table) = biallelic_locus();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &table.view());
    scratch.claim_row_for(5, InbreedingF::try_new(0.0).expect("an outbred sample"));
    scratch.prepare_contaminant_tables(1, 2);
}

/// One emission cache per scratch row, and a row past the locus's own count is refused —
/// which is what a driver walking one row too far would otherwise fill, leaving the last
/// row's cache holding another sample's observations.
#[test]
#[should_panic(expected = "row 2 is past the 2 this scratch was prepared for")]
fn an_emission_cache_past_the_prepared_rows_is_a_caller_bug() {
    let (alleles, table) = biallelic_locus();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &table.view());
    let _ = scratch.generic_row_mut(2);
}

/// **The contaminant tables are un-sized by every `prepare_for_locus`**, so a worker that
/// called a contaminated locus and then an uncontaminated one cannot read the first one's
/// frequencies at the second.
///
/// Sizing them and then preparing another locus is the whole of the fixture: `resize` alone
/// would have left a table of the right length holding the last locus's numbers, and only
/// the count coming back to zero says the tables were released rather than reused.
#[test]
fn preparing_a_locus_releases_the_last_ones_contaminant_tables() {
    let (alleles, table) = biallelic_locus();
    let outbred = InbreedingF::try_new(0.0).expect("an outbred sample");
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &table.view());
    scratch.claim_row_for(0, outbred);
    scratch.prepare_contaminant_tables(3, 4);
    assert_eq!(scratch.contaminant_batch_count(), 3);

    scratch.prepare_for_locus(1, &alleles, &table.view());
    assert_eq!(
        scratch.contaminant_batch_count(),
        0,
        "the next locus has to size them again, and an uncontaminated one never will"
    );
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
