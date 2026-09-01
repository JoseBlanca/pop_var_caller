//! What the mapper does with a called locus.
//!
//! **These test settled rules against a provisional interface.** If the stream hands over a
//! different shape than [`LocusEvidenceForOutput`], these fixtures change and their assertions
//! do not: what a no-call means, which samples reach `AN`, and where the counts land are
//! decisions with documents behind them.

use super::*;
use crate::ng::calling::{CandidateAlleles, ExpectedAlleleCopies};
use crate::ng::locus_generation::LocusKind;
use crate::ng::parameter_estimation::Provenance;
use crate::ng::types::{
    AlleleId, ContigId, GenomeRegion, Genotype, Motif, Phred, Ploidy, Position,
};
use crate::ng::vcf::{format_keys, info_column, sample_columns};

fn quality(phred: f32) -> Phred {
    Phred::try_new(phred).expect("a non-negative finite quality")
}

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("two copies is a genome")
}

fn region() -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(0),
        start: Position(1_000),
        end: Position(1_000),
    }
}

/// A two-allele SNP table.
fn alleles() -> CandidateAlleles {
    let mut table = CandidateAlleles::new(b"A".to_vec().into_boxed_slice(), LocusKind::Generic);
    table.admit(b"T".to_vec().into_boxed_slice());
    table
}

fn called(alleles: &[u16], genotype_quality: f32) -> SampleGenotypeCall {
    called_saying(alleles, genotype_quality, false)
}

/// The same, with **whether the sample's reads said anything** named — the bit that turns a
/// called sample into a `./.`, which since 2026-09-01 the calling loop mints onto the call
/// rather than the emission step being handed it separately.
fn called_saying(
    alleles: &[u16],
    genotype_quality: f32,
    reads_were_uninformative: bool,
) -> SampleGenotypeCall {
    SampleGenotypeCall::Called {
        genotype: Genotype::new(alleles.iter().copied().map(AlleleId).collect()),
        genotype_quality: quality(genotype_quality),
        reads_were_uninformative,
    }
}

/// A called locus over the two-allele table, with the given per-sample calls.
fn locus(calls: Vec<SampleGenotypeCall>, copies: Vec<f64>, converged: bool) -> LocusInference {
    let table = alleles();
    let expected = ExpectedAlleleCopies::new(copies, &table);
    LocusInference::new(
        region(),
        table,
        calls,
        expected,
        converged,
        3,
        Provenance::FittedHere,
        None,
        quality(0.0),
        None,
    )
}

fn evidence(
    samples: Vec<SampleEvidenceForOutput>,
    filter: FilterVerdict,
) -> LocusEvidenceForOutput {
    let mapq = (0..2)
        .map(|allele| {
            let reads: u64 = samples
                .iter()
                .map(|sample| u64::from(sample.allele_reads[allele]))
                .sum();
            MapqPool {
                reads,
                mapq_sum: reads * 60,
            }
        })
        .collect();
    LocusEvidenceForOutput {
        samples,
        allele_mapq: mapq,
        padding_base: None,
        corrected_site_quality: quality(120.0),
        artifact_penalties: None,
        repeat_tract: None,
        filter,
    }
}

fn sample(allele_reads: Vec<u32>, unexplained: u32) -> SampleEvidenceForOutput {
    SampleEvidenceForOutput {
        allele_reads,
        reads_no_written_allele_explains: unexplained,
    }
}

// ---------------------------------------------------------------------------
// The no-call rule — the settled part
// ---------------------------------------------------------------------------

#[test]
fn a_sample_whose_reads_said_nothing_is_written_as_a_no_call() {
    // **The rule the owner settled.** The loop calls such a sample — with a flat likelihood the
    // prior decides alone — and the file must not repeat that as a genotype.
    // **The bit is on the call now**, minted by the loop that scored the sample — which is
    // why the second sample here is built with `called_saying(..., true)` rather than by a
    // flag handed to the emission step beside it.
    let locus = locus(
        vec![called(&[0, 1], 40.0), called_saying(&[0, 0], 3.0, true)],
        vec![3.0, 1.0],
        true,
    );
    let record = assemble_record(
        &locus,
        evidence(
            vec![sample(vec![10, 9], 0), sample(vec![0, 0], 0)],
            FilterVerdict::Pass,
        ),
    );

    let columns = sample_columns(&record, diploid());
    assert_eq!(columns, "0/1:40:19:10,9\t./.:.:0:0,0");
    // The loop had given the second sample a genotype; the file does not repeat it.
    assert!(matches!(
        locus.per_sample[1],
        SampleGenotypeCall::Called { .. }
    ));
}

#[test]
fn a_sample_the_caller_declined_to_call_is_written_as_a_no_call_too() {
    // The other route, and a different fact: candidate selection cut an allele this sample's
    // reads had earned. One spelling, two reasons.
    let locus = locus(
        vec![called(&[0, 1], 40.0), SampleGenotypeCall::Missing],
        vec![3.0, 1.0],
        true,
    );
    let record = assemble_record(
        &locus,
        evidence(
            vec![sample(vec![10, 9], 0), sample(vec![2, 2], 5)],
            FilterVerdict::Pass,
        ),
    );

    // Its evidence is still written beside the missing call.
    assert_eq!(
        sample_columns(&record, diploid()),
        "0/1:40:19:10,9\t./.:.:9:2,2"
    );
}

#[test]
fn a_sample_with_reads_that_spoke_keeps_the_genotype_the_loop_gave_it() {
    let locus = locus(vec![called(&[1, 1], 55.0)], vec![0.2, 1.8], true);
    let record = assemble_record(
        &locus,
        evidence(vec![sample(vec![1, 14], 2)], FilterVerdict::Pass),
    );
    assert_eq!(sample_columns(&record, diploid()), "1/1:55:17:1,14");
}

#[test]
fn a_sample_written_as_a_no_call_leaves_an_and_ac_alone() {
    // The consequence that reaches the site's own numbers: `AN` counts written genotypes, so a
    // sample the loop scored but the file no-calls is out of it — while its copies stay in the
    // fit `AF` is taken from. The two denominators differ, and that is correct.
    let locus = locus(
        vec![called(&[0, 1], 40.0), called_saying(&[0, 0], 2.0, true)],
        vec![3.0, 1.0],
        true,
    );
    let record = assemble_record(
        &locus,
        evidence(
            vec![sample(vec![6, 6], 0), sample(vec![0, 0], 0)],
            FilterVerdict::Pass,
        ),
    );

    let info = info_column(&record);
    // AN is 2, not 4. AF is 1.0 of 4.0 copies — a quarter, not a half.
    assert!(
        info.starts_with("AF=0.250000;AC=1;AN=2;DP=12;"),
        "got {info}"
    );
}

// ---------------------------------------------------------------------------
// What the mapper reads off the real called locus
// ---------------------------------------------------------------------------

#[test]
fn the_alleles_and_the_fitted_copies_come_from_the_locus_itself() {
    let locus = locus(vec![called(&[0, 1], 40.0)], vec![1.25, 0.75], true);
    let record = assemble_record(
        &locus,
        evidence(vec![sample(vec![5, 5], 0)], FilterVerdict::Pass),
    );

    assert_eq!(record.reference(), b"A");
    assert_eq!(&*record.alternatives()[0], b"T");
    assert_eq!(record.expected_copies(), [1.25, 0.75]);
    assert_eq!(record.region(), region());
}

#[test]
fn an_unconverged_locus_reaches_the_file_on_its_filter() {
    let locus = locus(vec![called(&[0, 1], 40.0)], vec![1.0, 1.0], false);
    let record = assemble_record(
        &locus,
        evidence(vec![sample(vec![5, 5], 0)], FilterVerdict::EmDidNotConverge),
    );
    assert_eq!(record.filter(), FilterVerdict::EmDidNotConverge);
}

#[test]
fn a_tract_annotation_makes_it_a_tract_record() {
    let locus = locus(vec![called(&[0, 1], 40.0)], vec![1.0, 1.0], true);
    let mut inputs = evidence(vec![sample(vec![5, 5], 0)], FilterVerdict::Pass);
    inputs.repeat_tract = Some(TractAnnotation::new(
        Motif::new(b"AT").expect("a two-base motif"),
    ));
    let record = assemble_record(&locus, inputs);

    assert!(record.is_repeat_tract());
    assert_eq!(format_keys(&record), "GT:GQ:DP:AD:REPCN");
}

// ---------------------------------------------------------------------------
// The joins that could silently lie
// ---------------------------------------------------------------------------

#[test]
#[should_panic(expected = "gathered over different cohorts")]
fn evidence_for_a_different_number_of_samples_is_refused() {
    let locus = locus(
        vec![called(&[0, 1], 40.0), called(&[0, 0], 40.0)],
        vec![3.0, 1.0],
        true,
    );
    let _ = assemble_record(
        &locus,
        evidence(vec![sample(vec![5, 5], 0)], FilterVerdict::Pass),
    );
}

#[test]
#[should_panic(expected = "one of them was set by hand")]
fn a_filter_that_disagrees_with_the_loop_about_convergence_is_refused() {
    // The one verdict derivable from the locus itself. If the two could disagree, the file
    // would state a convergence the loop never reported.
    let locus = locus(vec![called(&[0, 1], 40.0)], vec![1.0, 1.0], false);
    let _ = assemble_record(
        &locus,
        evidence(vec![sample(vec![5, 5], 0)], FilterVerdict::Pass),
    );
}

#[test]
#[should_panic(expected = "AD is written one entry per allele")]
fn per_sample_counts_of_the_wrong_width_are_refused() {
    let locus = locus(vec![called(&[0, 1], 40.0)], vec![1.0, 1.0], true);
    let mut inputs = evidence(vec![sample(vec![5, 5], 0)], FilterVerdict::Pass);
    inputs.samples[0].allele_reads = vec![5, 5, 5];
    let _ = assemble_record(&locus, inputs);
}
