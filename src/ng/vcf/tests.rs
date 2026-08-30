//! What a record can and cannot be made to say.
//!
//! The two halves of Milestone A's checkpoint: every shape the format has to express is built
//! here (both kinds of locus, a no-call, a refused tract locus, a full-tract deletion at either
//! end of a contig), and every shape it forbids is shown to be refused at construction or
//! impossible to spell at all.

use super::*;
use crate::ng::types::{AlleleId, ContigId};

fn region(start: u64, end: u64) -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(0),
        start: Position(start),
        end: Position(end),
    }
}

fn quality(phred: f32) -> Phred {
    Phred::try_new(phred).expect("a non-negative finite quality")
}

fn called(alleles: &[u16], genotype_quality: f32) -> SampleCall {
    SampleCall::Called {
        genotype: Genotype::new(alleles.iter().copied().map(AlleleId).collect()),
        genotype_quality: quality(genotype_quality),
    }
}

fn allele(bases: &[u8]) -> Box<[u8]> {
    bases.to_vec().into_boxed_slice()
}

fn motif(bases: &[u8]) -> Motif {
    Motif::new(bases).expect("a motif of one to six bases")
}

/// Pool `mapq_sum` over `reads` reads — the identity `VcfRecord::new` enforces means the read
/// count has to match what the samples' `AD` attributes to that allele.
fn pool(reads: u64, mean_mapq: u64) -> MapqPool {
    MapqPool {
        reads,
        mapq_sum: reads * mean_mapq,
    }
}

/// A biallelic SNP: two samples, one heterozygous and one homozygous reference.
fn a_snp_record() -> VcfRecord {
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"A"), allele(b"T")],
        vec![1.4, 0.6],
        vec![
            SampleColumn {
                call: called(&[0, 1], 42.0),
                read_counts: SampleReadCounts::new(vec![10, 9], 0),
            },
            SampleColumn {
                call: called(&[0, 0], 60.0),
                read_counts: SampleReadCounts::new(vec![20, 0], 0),
            },
        ],
        vec![pool(30, 60), pool(9, 60)],
        None,
        quality(300.0),
        Some(ArtifactPenalties {
            allele_balance: Phred::ZERO,
            strand_and_read_position: quality(2.5),
        }),
        FilterVerdict::Pass,
        None,
    )
}

/// A repeat tract of `CA` units: the reference at eight copies, one alternative at six.
fn a_tract_record() -> VcfRecord {
    VcfRecord::new(
        region(2_000, 2_015),
        vec![allele(b"CACACACACACACACA"), allele(b"CACACACACACA")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: called(&[0, 1], 35.0),
            read_counts: SampleReadCounts::new(vec![12, 11], 7),
        }],
        vec![pool(12, 58), pool(11, 55)],
        None,
        quality(150.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(motif(b"CA"))),
    )
}

// ---------------------------------------------------------------------------
// The shapes the format has to express
// ---------------------------------------------------------------------------

#[test]
fn both_kinds_of_locus_are_one_record_shape() {
    let snp = a_snp_record();
    let tract = a_tract_record();

    // The only structural difference between the two is the tract annotation.
    assert!(!snp.is_repeat_tract());
    assert!(tract.is_repeat_tract());
    assert_eq!(snp.reference(), b"A");
    assert_eq!(tract.reference(), b"CACACACACACACACA");
    assert_eq!(snp.alternatives().len(), 1);
    assert_eq!(tract.alternatives().len(), 1);
}

#[test]
fn a_sample_can_have_reads_and_no_call() {
    // The distinction the format exists to keep: a sample no-called on 7 reads and one
    // no-called on none are different facts, and both are spelled `./.`.
    let record = VcfRecord::new(
        region(10, 10),
        vec![allele(b"G"), allele(b"C")],
        vec![2.0, 0.0],
        vec![
            SampleColumn {
                call: SampleCall::NoCall,
                read_counts: SampleReadCounts::new(vec![4, 3], 0),
            },
            SampleColumn {
                call: SampleCall::NoCall,
                read_counts: SampleReadCounts::new(vec![0, 0], 0),
            },
        ],
        vec![pool(4, 60), pool(3, 44)],
        None,
        quality(0.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    let with_reads = &record.sample_columns()[0];
    let without = &record.sample_columns()[1];
    assert!(with_reads.call.is_no_call() && without.call.is_no_call());
    assert_eq!(with_reads.read_counts.depth(), 7);
    assert_eq!(without.read_counts.depth(), 0);
    // No genotype means no genotype quality — absent, not zero.
    assert!(with_reads.call.genotype_quality().is_none());
}

#[test]
fn a_refused_tract_locus_is_expressible() {
    // The shape spec §8 fixes for a locus the caller looked at and could not call: the
    // reference alone, every sample no-called, quality zero, and the filter saying why.
    let record = VcfRecord::new(
        region(500, 511),
        vec![allele(b"ATATATATATAT")],
        vec![0.0],
        vec![SampleColumn {
            call: SampleCall::NoCall,
            read_counts: SampleReadCounts::new(vec![0], 2),
        }],
        vec![MapqPool::default()],
        None,
        quality(0.0),
        None,
        FilterVerdict::LowDepth,
        Some(TractAnnotation::new(motif(b"AT"))),
    );

    assert!(record.alternatives().is_empty());
    assert_eq!(record.filter().as_str(), "lowDepth");
    assert_eq!(record.site_quality.get(), 0.0);
    // Two reads were seen and no allele explains them — the counts still say so.
    assert_eq!(
        record.sample_columns()[0].read_counts.unexplained_reads(),
        2
    );
}

#[test]
fn every_filter_value_has_its_written_spelling() {
    // Five values in one namespace, spelled as the two production writers spell them.
    assert_eq!(FilterVerdict::Pass.as_str(), "PASS");
    assert_eq!(FilterVerdict::EmDidNotConverge.as_str(), "EMNoConv");
    assert_eq!(FilterVerdict::NotPeriodic.as_str(), "notPeriodic");
    assert_eq!(FilterVerdict::TooManyAlleles.as_str(), "tooManyAlleles");
    assert_eq!(FilterVerdict::LowDepth.as_str(), "lowDepth");
}

/// A full-tract deletion away from the contig's start: padded from the left.
fn a_full_tract_deletion() -> VcfRecord {
    VcfRecord::new(
        region(700, 705),
        vec![allele(b"ATATAT"), allele(b"")],
        vec![0.4, 1.6],
        vec![SampleColumn {
            call: called(&[1, 1], 30.0),
            read_counts: SampleReadCounts::new(vec![0, 14], 0),
        }],
        vec![MapqPool::default(), pool(14, 57)],
        Some(PaddingBase::Left(b'C')),
        quality(90.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(motif(b"AT"))),
    )
}

#[test]
fn an_empty_alternative_is_carried_with_the_base_it_will_be_padded_with() {
    let record = a_full_tract_deletion();
    assert_eq!(record.alternatives()[0].len(), 0);
    assert_eq!(record.padding_base(), Some(PaddingBase::Left(b'C')));
    assert_eq!(record.padding_base().expect("a padding base").base(), b'C');
}

#[test]
fn a_deletion_at_the_first_base_of_a_contig_is_padded_from_the_right() {
    // There is no base to the left of position 1, so the rule appends the base to the right
    // and POS does not move — where production's tract writer invents an `N` instead.
    let record = VcfRecord::new(
        region(1, 4),
        vec![allele(b"ATAT"), allele(b"")],
        vec![0.5, 1.5],
        vec![SampleColumn {
            call: called(&[1, 1], 25.0),
            read_counts: SampleReadCounts::new(vec![0, 8], 0),
        }],
        vec![MapqPool::default(), pool(8, 50)],
        Some(PaddingBase::Right(b'G')),
        quality(40.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(motif(b"AT"))),
    );

    assert_eq!(record.padding_base(), Some(PaddingBase::Right(b'G')));
    assert_eq!(record.region().start, Position(1));
}

#[test]
fn a_record_with_no_empty_allele_carries_no_padding_base() {
    assert_eq!(a_snp_record().padding_base(), None);
    assert_eq!(a_tract_record().padding_base(), None);
}

// ---------------------------------------------------------------------------
// Derivations that cannot drift
// ---------------------------------------------------------------------------

#[test]
fn the_str_flag_is_the_annotation_and_cannot_disagree_with_it() {
    // There is no separate flag to fall out of step with `RU` and `PERIOD`: one is the other.
    let tract = a_tract_record();
    assert_eq!(tract.is_repeat_tract(), tract.repeat_tract().is_some());

    let annotation = tract.repeat_tract().expect("a tract record");
    assert_eq!(annotation.motif(), b"CA");
    assert_eq!(annotation.period(), annotation.motif().len());
}

#[test]
fn repeat_copies_come_from_the_allele_the_record_writes() {
    let tract = a_tract_record();
    let annotation = tract.repeat_tract().expect("a tract record");

    // Sixteen bases of a two-base unit is eight copies; twelve is six.
    assert_eq!(annotation.repeat_copies_of(tract.reference()), 8);
    assert_eq!(annotation.repeat_copies_of(&tract.alternatives()[0]), 6);
}

#[test]
fn a_partial_final_unit_reports_the_whole_copies_it_has() {
    let annotation = TractAnnotation::new(motif(b"CAG"));
    // Eight bases of a three-base unit: two whole copies and two bases over.
    assert_eq!(annotation.repeat_copies_of(b"CAGCAGCA"), 2);
}

#[test]
fn a_deleted_tract_holds_no_repeat_copies() {
    let record = a_full_tract_deletion();
    let annotation = record.repeat_tract().expect("a tract record");
    assert_eq!(annotation.repeat_copies_of(&record.alternatives()[0]), 0);
}

#[test]
fn unexplained_reads_are_the_depth_the_alleles_do_not_account_for() {
    // The per-sample artifact signal: 30 reads seen, 23 attributed, 7 explained by nothing
    // this record writes.
    let tract = a_tract_record();
    let counts = &tract.sample_columns()[0].read_counts;
    assert_eq!(counts.depth(), 30);
    assert_eq!(counts.allele_reads(), [12, 11]);
    assert_eq!(counts.unexplained_reads(), 7);
}

#[test]
fn every_read_explained_leaves_nothing_unexplained() {
    let counts = SampleReadCounts::new(vec![10, 9], 0);
    assert_eq!(counts.unexplained_reads(), 0);
    assert_eq!(counts.depth(), 19);
}

#[test]
fn the_depth_is_derived_and_so_cannot_contradict_the_allele_counts() {
    // `DP` is not stored: it is `ΣAD` plus the reads no written allele explains, so the
    // failure the old shape allowed — a depth passed in below the counts it has to cover,
    // which no VCF parser would reject — cannot be spelled at all.
    let counts = SampleReadCounts::new(vec![3, 4, 5], 6);
    assert_eq!(counts.depth(), 18);
    assert_eq!(
        counts.depth(),
        counts.allele_reads().iter().sum::<u32>() + 6
    );
}

#[test]
fn a_sample_whose_reads_all_miss_the_written_alleles_still_has_depth() {
    // Every read explained by an allele candidate selection dropped: `AD` is all zeroes and
    // `DP` is not, which is the per-sample signal the difference exists to publish.
    let counts = SampleReadCounts::new(vec![0, 0], 11);
    assert_eq!(counts.depth(), 11);
    assert_eq!(counts.unexplained_reads(), 11);
}

#[test]
fn a_pooled_mapping_quality_of_no_reads_has_no_mean() {
    assert_eq!(MapqPool::default().mean(), None);
    assert_eq!(pool(4, 60).mean(), Some(60.0));
}

#[test]
fn the_expected_copies_run_parallel_to_the_alleles() {
    // What `AF` is written from — the loop's fit, one entry per allele, carried rather than
    // recomputed from the genotypes.
    let snp = a_snp_record();
    assert_eq!(snp.expected_copies().len(), snp.alleles().len());
    assert_eq!(snp.expected_copies(), [1.4, 0.6]);
}

// ---------------------------------------------------------------------------
// The shapes the format forbids
// ---------------------------------------------------------------------------

#[test]
#[should_panic(expected = "AD is written one entry per allele")]
fn a_read_count_vector_wider_than_the_allele_table_is_refused() {
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: called(&[0, 1], 20.0),
            read_counts: SampleReadCounts::new(vec![5, 5, 5], 0),
        }],
        vec![pool(5, 60), pool(5, 60)],
        None,
        quality(10.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "is a corrupt count rather than a deep locus")]
fn a_read_total_too_large_for_the_depth_column_is_refused() {
    // The one arithmetic failure the derived depth still admits: two counts that cannot both
    // be real, caught rather than wrapped into a small depth.
    SampleReadCounts::new(vec![u32::MAX, u32::MAX], 1);
}

#[test]
#[should_panic(expected = "one entry per allele, reference first")]
fn a_mapping_quality_pool_of_the_wrong_width_is_refused() {
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: called(&[0, 1], 20.0),
            read_counts: SampleReadCounts::new(vec![5, 5], 0),
        }],
        vec![pool(5, 60)],
        None,
        quality(10.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "mean two different pools")]
fn mapping_qualities_pooled_over_other_reads_than_ad_counts_are_refused() {
    // The identity production keeps by construction — its `AD` and its MQ denominator are one
    // field — asserted here, because two totals mean an `MQDIFF` over reads the `AD` beside it
    // does not describe.
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: called(&[0, 1], 20.0),
            read_counts: SampleReadCounts::new(vec![5, 5], 0),
        }],
        vec![pool(5, 60), pool(4, 60)],
        None,
        quality(10.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "so that AF names the alleles this record holds")]
fn expected_copies_of_the_wrong_width_are_refused() {
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0],
        vec![SampleColumn {
            call: called(&[0, 1], 20.0),
            read_counts: SampleReadCounts::new(vec![5, 5], 0),
        }],
        vec![pool(5, 60), pool(5, 60)],
        None,
        quality(10.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "names an allele this record does not hold")]
fn a_call_naming_an_allele_past_the_table_is_refused() {
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: called(&[1, 2], 20.0),
            read_counts: SampleReadCounts::new(vec![5, 5], 0),
        }],
        vec![pool(5, 60), pool(5, 60)],
        None,
        quality(10.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "has lost the table rather than having none")]
fn a_record_with_no_alleles_is_refused() {
    VcfRecord::new(
        region(1_000, 1_000),
        Vec::new(),
        Vec::new(),
        vec![SampleColumn {
            call: SampleCall::NoCall,
            read_counts: SampleReadCounts::new(Vec::new(), 0),
        }],
        Vec::new(),
        None,
        quality(0.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "reaches the REF column as an unparseable record")]
fn an_empty_reference_allele_is_refused() {
    // The gate `CandidateAlleles::new` keeps upstream, named after this very column — the
    // writer must not be the weaker one.
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b""), allele(b"T")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: called(&[0, 1], 20.0),
            read_counts: SampleReadCounts::new(vec![5, 5], 0),
        }],
        vec![pool(5, 60), pool(5, 60)],
        Some(PaddingBase::Left(b'C')),
        quality(10.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "describes a different stretch of reference")]
fn a_reference_allele_that_does_not_span_the_region_is_refused() {
    // Three bases of REF over a one-base span: the record's POS would claim ground the
    // sequence does not cover, and the line would still parse.
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"ATG"), allele(b"T")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: called(&[0, 1], 20.0),
            read_counts: SampleReadCounts::new(vec![5, 5], 0),
        }],
        vec![pool(5, 60), pool(5, 60)],
        None,
        quality(10.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "a record naming no sample has lost them")]
fn a_record_with_no_samples_is_refused() {
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"A")],
        vec![0.0],
        Vec::new(),
        vec![MapqPool::default()],
        None,
        quality(0.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "cannot run backwards")]
fn a_backwards_region_is_refused() {
    VcfRecord::new(
        region(100, 99),
        vec![allele(b"A")],
        vec![0.0],
        vec![SampleColumn {
            call: SampleCall::NoCall,
            read_counts: SampleReadCounts::new(vec![0], 0),
        }],
        vec![MapqPool::default()],
        None,
        quality(0.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "a padding base is carried exactly when some allele is empty")]
fn an_empty_allele_without_a_padding_base_is_refused() {
    // Without the flanking base the encoder cannot write this record at all: the base is
    // outside the span, and nothing downstream holds the reference to fetch it.
    VcfRecord::new(
        region(700, 705),
        vec![allele(b"ATATAT"), allele(b"")],
        vec![0.4, 1.6],
        vec![SampleColumn {
            call: called(&[1, 1], 30.0),
            read_counts: SampleReadCounts::new(vec![0, 14], 0),
        }],
        vec![MapqPool::default(), pool(14, 57)],
        None,
        quality(90.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(motif(b"AT"))),
    );
}

#[test]
#[should_panic(expected = "a padding base is carried exactly when some allele is empty")]
fn a_padding_base_with_no_empty_allele_is_refused() {
    VcfRecord::new(
        region(1_000, 1_000),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: called(&[0, 1], 20.0),
            read_counts: SampleReadCounts::new(vec![5, 5], 0),
        }],
        vec![pool(5, 60), pool(5, 60)],
        Some(PaddingBase::Left(b'C')),
        quality(10.0),
        None,
        FilterVerdict::Pass,
        None,
    );
}

#[test]
#[should_panic(expected = "the padding base is the one to the left of the span")]
fn a_right_hand_padding_base_away_from_the_contig_start_is_refused() {
    VcfRecord::new(
        region(700, 705),
        vec![allele(b"ATATAT"), allele(b"")],
        vec![0.4, 1.6],
        vec![SampleColumn {
            call: called(&[1, 1], 30.0),
            read_counts: SampleReadCounts::new(vec![0, 14], 0),
        }],
        vec![MapqPool::default(), pool(14, 57)],
        Some(PaddingBase::Right(b'G')),
        quality(90.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(motif(b"AT"))),
    );
}

#[test]
#[should_panic(expected = "the padding base is the one to the left of the span")]
fn a_left_hand_padding_base_at_the_contig_start_is_refused() {
    // There is no base to the left of position 1 to have taken.
    VcfRecord::new(
        region(1, 4),
        vec![allele(b"ATAT"), allele(b"")],
        vec![0.5, 1.5],
        vec![SampleColumn {
            call: called(&[1, 1], 25.0),
            read_counts: SampleReadCounts::new(vec![0, 8], 0),
        }],
        vec![MapqPool::default(), pool(8, 50)],
        Some(PaddingBase::Left(b'C')),
        quality(40.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(motif(b"AT"))),
    );
}

#[test]
fn a_motif_cannot_be_empty_so_a_tract_annotation_cannot_be() {
    // The refusal that used to live here is gone because the state is now unrepresentable:
    // `Motif` is one to six bases by construction, so no `PERIOD=0` can be written and
    // `repeat_copies_of` can never divide by zero.
    assert!(Motif::new(b"").is_err());
    assert!(Motif::new(b"ACGTACGT").is_err());
}
