//! The fixed columns, asserted byte for byte.
//!
//! **Golden strings rather than field-by-field checks**, because the failure this step exists to
//! catch is a record that writes cleanly and says something different from what it holds — a
//! shifted `POS`, a padding base on the wrong side, a quality rendered to a different precision.
//! Every one of those parses.

use super::*;
use crate::ng::calling::quality::artifact_correction::ArtifactPenalties;
use crate::ng::types::{AlleleId, ContigId, GenomeRegion, Genotype, Motif, Phred, Position};
use crate::ng::vcf::{
    FilterVerdict, MapqPool, SampleCall, SampleColumn, SampleReadCounts, TractAnnotation,
};

fn contigs() -> Vec<HeaderContig> {
    vec![
        HeaderContig {
            name: "chr1".to_string(),
            length: 248_956_422,
            md5: None,
        },
        HeaderContig {
            name: "chr2".to_string(),
            length: 242_193_529,
            md5: None,
        },
    ]
}

fn region_on(contig: u32, start: u64, end: u64) -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(contig),
        start: Position(start),
        end: Position(end),
    }
}

fn quality(phred: f32) -> Phred {
    Phred::try_new(phred).expect("a non-negative finite quality")
}

fn allele(bases: &[u8]) -> Box<[u8]> {
    bases.to_vec().into_boxed_slice()
}

fn one_sample(allele_reads: Vec<u32>, unexplained: u32, genotype: &[u16]) -> Vec<SampleColumn> {
    vec![SampleColumn {
        call: SampleCall::Called {
            genotype: Genotype::new(genotype.iter().copied().map(AlleleId).collect()),
            genotype_quality: quality(40.0),
        },
        read_counts: SampleReadCounts::new(allele_reads, unexplained),
    }]
}

/// Pools matching what the samples' `AD` attributes, which the record type requires.
fn pools(reads: &[u64]) -> Vec<MapqPool> {
    reads
        .iter()
        .map(|count| MapqPool {
            reads: *count,
            mapq_sum: count * 60,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// The ordinary shapes
// ---------------------------------------------------------------------------

#[test]
fn a_snp_writes_its_seven_columns() {
    let record = VcfRecord::new(
        region_on(0, 1_000, 1_000),
        vec![allele(b"A"), allele(b"T")],
        vec![1.4, 0.6],
        one_sample(vec![10, 9], 0, &[0, 1]),
        pools(&[10, 9]),
        None,
        quality(300.0),
        Some(ArtifactPenalties {
            allele_balance: Phred::ZERO,
            strand_and_read_position: quality(2.5),
        }),
        FilterVerdict::Pass,
        None,
    );

    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr1\t1000\t.\tA\tT\t300.0\tPASS"
    );
}

#[test]
fn the_contig_is_named_from_the_header_table_the_record_indexes() {
    let record = VcfRecord::new(
        region_on(1, 55, 55),
        vec![allele(b"G"), allele(b"C")],
        vec![1.0, 1.0],
        one_sample(vec![5, 5], 0, &[0, 1]),
        pools(&[5, 5]),
        None,
        quality(12.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr2\t55\t.\tG\tC\t12.0\tPASS"
    );
}

#[test]
fn several_alternatives_are_comma_joined_in_table_order() {
    let record = VcfRecord::new(
        region_on(0, 20, 20),
        vec![allele(b"A"), allele(b"C"), allele(b"G")],
        vec![1.0, 0.6, 0.4],
        one_sample(vec![8, 6, 4], 0, &[1, 2]),
        pools(&[8, 6, 4]),
        None,
        quality(99.5),
        None,
        FilterVerdict::Pass,
        None,
    );

    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr1\t20\t.\tA\tC,G\t99.5\tPASS"
    );
}

#[test]
fn a_record_with_no_alternative_writes_a_missing_alt() {
    // A tract the caller looked at and could not call: the reference alone, on its filter.
    let record = VcfRecord::new(
        region_on(0, 500, 511),
        vec![allele(b"ATATATATATAT")],
        vec![0.0],
        vec![SampleColumn {
            call: SampleCall::NoCall,
            read_counts: SampleReadCounts::new(vec![0], 2),
        }],
        pools(&[0]),
        None,
        quality(0.0),
        None,
        FilterVerdict::LowDepth,
        Some(TractAnnotation::new(
            Motif::new(b"AT").expect("a two-base motif"),
        )),
    );

    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr1\t500\t.\tATATATATATAT\t.\t0.0\tlowDepth"
    );
}

#[test]
fn an_unconverged_locus_writes_its_filter() {
    let record = VcfRecord::new(
        region_on(0, 7, 7),
        vec![allele(b"T"), allele(b"A")],
        vec![1.0, 1.0],
        one_sample(vec![3, 3], 0, &[0, 1]),
        pools(&[3, 3]),
        None,
        quality(31.25),
        None,
        FilterVerdict::EmDidNotConverge,
        None,
    );

    // 31.25 renders to one decimal, and the tie rounds to even: 31.2, not 31.3.
    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr1\t7\t.\tT\tA\t31.2\tEMNoConv"
    );
}

// ---------------------------------------------------------------------------
// The padding rule — this step's silent-failure site
// ---------------------------------------------------------------------------

#[test]
fn a_deletion_takes_the_base_to_its_left_and_moves_back_one() {
    // The ordinary VCF deletion: REF gains the flanking base, ALT becomes that base alone, and
    // POS moves from 700 to 699. A POS that failed to move would parse and name a variant one
    // base to the right of the real one.
    let record = VcfRecord::new(
        region_on(0, 700, 705),
        vec![allele(b"ATATAT"), allele(b"")],
        vec![0.4, 1.6],
        one_sample(vec![0, 14], 0, &[1, 1]),
        pools(&[0, 14]),
        Some(PaddingBase::Left(b'C')),
        quality(90.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(
            Motif::new(b"AT").expect("a two-base motif"),
        )),
    );

    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr1\t699\t.\tCATATAT\tC\t90.0\tPASS"
    );
}

#[test]
fn a_deletion_at_the_first_base_of_a_contig_takes_the_base_to_its_right_and_stays() {
    // There is no base to the left of position 1. VCF 4.4 appends the base to the right and
    // leaves POS alone — where production's repeat-tract writer invents the letter `N` at an
    // unshifted position, a base the reference does not contain. This is the case that rule
    // exists for, and the output is what it should have been.
    let record = VcfRecord::new(
        region_on(0, 1, 4),
        vec![allele(b"ATAT"), allele(b"")],
        vec![0.5, 1.5],
        one_sample(vec![0, 8], 0, &[1, 1]),
        pools(&[0, 8]),
        Some(PaddingBase::Right(b'G')),
        quality(40.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(
            Motif::new(b"AT").expect("a two-base motif"),
        )),
    );

    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr1\t1\t.\tATATG\tG\t40.0\tPASS"
    );
}

#[test]
fn every_allele_is_padded_and_not_only_the_empty_one() {
    // VCF states a deletion by giving all alleles a shared flanking base. Padding only the
    // empty allele would describe a different variant and would still parse.
    let record = VcfRecord::new(
        region_on(0, 300, 303),
        vec![allele(b"GGGG"), allele(b""), allele(b"GG")],
        vec![0.6, 0.9, 0.5],
        one_sample(vec![4, 5, 3], 0, &[1, 2]),
        pools(&[4, 5, 3]),
        Some(PaddingBase::Left(b'T')),
        quality(55.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr1\t299\t.\tTGGGG\tT,TGG\t55.0\tPASS"
    );
}

#[test]
fn a_record_with_no_empty_allele_is_written_unpadded_at_its_own_position() {
    // The negative of the two tests above: nothing is prefixed and POS does not move.
    let record = VcfRecord::new(
        region_on(0, 42, 43),
        vec![allele(b"AC"), allele(b"AG")],
        vec![1.0, 1.0],
        one_sample(vec![6, 6], 0, &[0, 1]),
        pools(&[6, 6]),
        None,
        quality(70.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    assert_eq!(
        fixed_columns(&record, &contigs()),
        "chr1\t42\t.\tAC\tAG\t70.0\tPASS"
    );
}

// ---------------------------------------------------------------------------
// Formatting details that are part of the format
// ---------------------------------------------------------------------------

#[test]
fn a_quality_is_written_to_one_decimal_whatever_its_size() {
    let cases = [
        (0.0_f32, "0.0"),
        (9.0, "9.0"),
        (30.0, "30.0"),
        (0.04, "0.0"),
        (9999.0, "9999.0"),
    ];
    for (phred, expected) in cases {
        let record = VcfRecord::new(
            region_on(0, 5, 5),
            vec![allele(b"A"), allele(b"T")],
            vec![1.0, 1.0],
            one_sample(vec![2, 2], 0, &[0, 1]),
            pools(&[2, 2]),
            None,
            quality(phred),
            None,
            FilterVerdict::Pass,
            None,
        );
        let columns = fixed_columns(&record, &contigs());
        let quality_column = columns.split('\t').nth(5).expect("a QUAL column");
        assert_eq!(quality_column, expected, "for a Phred of {phred}");
    }
}

#[test]
fn the_id_column_is_always_missing() {
    // Nothing populates it: neither production writer does, and no consumer has asked.
    let record = VcfRecord::new(
        region_on(0, 5, 5),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        one_sample(vec![2, 2], 0, &[0, 1]),
        pools(&[2, 2]),
        None,
        quality(1.0),
        None,
        FilterVerdict::Pass,
        None,
    );
    assert_eq!(
        fixed_columns(&record, &contigs())
            .split('\t')
            .nth(2)
            .expect("an ID column"),
        MISSING_FIELD
    );
}

#[test]
fn the_columns_are_tab_separated_with_no_trailing_separator() {
    let record = VcfRecord::new(
        region_on(0, 5, 5),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        one_sample(vec![2, 2], 0, &[0, 1]),
        pools(&[2, 2]),
        None,
        quality(1.0),
        None,
        FilterVerdict::Pass,
        None,
    );
    let columns = fixed_columns(&record, &contigs());
    assert_eq!(columns.split('\t').count(), 7);
    assert!(!columns.ends_with('\t'));
    assert!(!columns.contains('\n'));
}

#[test]
#[should_panic(expected = "built from different references")]
fn a_contig_the_header_does_not_hold_is_refused() {
    let record = VcfRecord::new(
        region_on(9, 5, 5),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        one_sample(vec![2, 2], 0, &[0, 1]),
        pools(&[2, 2]),
        None,
        quality(1.0),
        None,
        FilterVerdict::Pass,
        None,
    );
    // Bound rather than discarded: `fixed_columns` is `#[must_use]`, correctly — throwing away
    // an encoding is a defect everywhere except here, where the panic is the subject.
    let _ = fixed_columns(&record, &contigs());
}
