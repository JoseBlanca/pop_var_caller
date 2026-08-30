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

// ---------------------------------------------------------------------------
// INFO
// ---------------------------------------------------------------------------

#[test]
fn a_snp_writes_every_site_annotation_in_order() {
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

    // One sample, called 0/1: AN 2, AC 1. Copies 1.4 and 0.6 of a total 2.0, so AF 0.3.
    // Both pools mean 60, so MQDIFF is 0.
    assert_eq!(
        info_column(&record),
        "AF=0.300000;AC=1;AN=2;DP=19;ABPEN=0.0;SPPEN=2.5;MQREF=60.00;MQALT=60.00;MQDIFF=0.00"
    );
}

#[test]
fn a_tract_record_carries_the_str_flag_beside_its_motif_and_period() {
    let record = VcfRecord::new(
        region_on(0, 2_000, 2_015),
        vec![allele(b"CACACACACACACACA"), allele(b"CACACACACACA")],
        vec![1.0, 1.0],
        one_sample(vec![12, 11], 7, &[0, 1]),
        pools(&[12, 11]),
        None,
        quality(150.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(
            Motif::new(b"CA").expect("a two-base motif"),
        )),
    );

    let info = info_column(&record);
    assert!(info.ends_with(";STR;RU=CA;PERIOD=2"), "got {info}");
    // The three travel together or not at all — there is no record with one and not the others.
    assert_eq!(info.matches("STR").count(), 1);
}

#[test]
fn a_snp_record_carries_no_tract_annotation() {
    let info = info_column(&a_biallelic_snp());
    assert!(!info.contains("STR"));
    assert!(!info.contains("RU="));
    assert!(!info.contains("PERIOD="));
}

/// A plain two-allele SNP with one heterozygous sample.
fn a_biallelic_snp() -> VcfRecord {
    VcfRecord::new(
        region_on(0, 30, 30),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        one_sample(vec![5, 5], 0, &[0, 1]),
        pools(&[5, 5]),
        None,
        quality(50.0),
        None,
        FilterVerdict::Pass,
        None,
    )
}

#[test]
fn the_frequencies_are_normalised_over_the_copies_and_not_over_an() {
    // **The distinction this test exists for.** Two samples take part: one is called 0/1, the
    // other is written as a no-call because its reads said nothing (spec §7.1) — but the loop
    // scored it, so its copies are in the fit. AN counts only the called sample's two alleles,
    // while the copies total 4. Normalising over AN would give frequencies summing to 2.
    let record = VcfRecord::new(
        region_on(0, 90, 90),
        vec![allele(b"A"), allele(b"T")],
        vec![3.0, 1.0],
        vec![
            SampleColumn {
                call: SampleCall::Called {
                    genotype: Genotype::new(vec![AlleleId(0), AlleleId(1)]),
                    genotype_quality: quality(30.0),
                },
                read_counts: SampleReadCounts::new(vec![6, 6], 0),
            },
            SampleColumn {
                call: SampleCall::NoCall,
                read_counts: SampleReadCounts::new(vec![0, 0], 0),
            },
        ],
        pools(&[6, 6]),
        None,
        quality(80.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    let info = info_column(&record);
    // 1.0 of 4.0 copies is 0.25 — not 1.0/2 = 0.5, which normalising over AN would give.
    assert!(info.starts_with("AF=0.250000;AC=1;AN=2;"), "got {info}");
}

#[test]
fn a_no_called_sample_is_in_neither_ac_nor_an_but_its_reads_are_in_dp() {
    let record = VcfRecord::new(
        region_on(0, 90, 90),
        vec![allele(b"A"), allele(b"T")],
        vec![2.0, 0.0],
        vec![
            SampleColumn {
                call: SampleCall::Called {
                    genotype: Genotype::new(vec![AlleleId(0), AlleleId(0)]),
                    genotype_quality: quality(30.0),
                },
                read_counts: SampleReadCounts::new(vec![8, 0], 0),
            },
            SampleColumn {
                call: SampleCall::NoCall,
                read_counts: SampleReadCounts::new(vec![2, 1], 4),
            },
        ],
        pools(&[10, 1]),
        None,
        quality(20.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    let info = info_column(&record);
    // AN is 2, not 4: the no-called sample contributes no called alleles. Its seven reads are
    // still in DP — the evidence exists even where the call does not.
    assert!(info.contains(";AC=0;AN=2;DP=15;"), "got {info}");
}

#[test]
fn an_allele_no_read_reached_writes_a_missing_entry_rather_than_a_zero() {
    // A mean over no reads is absent, not zero — zero would claim every read mapped as badly
    // as possible. `MQALT` keeps its slot and writes `.`; `MQDIFF` follows it.
    let record = VcfRecord::new(
        region_on(0, 12, 12),
        vec![allele(b"A"), allele(b"T"), allele(b"G")],
        vec![1.5, 0.5, 0.0],
        one_sample(vec![7, 3, 0], 0, &[0, 1]),
        vec![
            MapqPool {
                reads: 7,
                mapq_sum: 420,
            },
            MapqPool {
                reads: 3,
                mapq_sum: 150,
            },
            MapqPool::default(),
        ],
        None,
        quality(45.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    let info = info_column(&record);
    // REF mean 60, first ALT mean 50 (difference -10), second ALT unreached.
    assert!(
        info.contains("MQREF=60.00;MQALT=50.00,.;MQDIFF=-10.00,."),
        "got {info}"
    );
}

#[test]
fn a_reference_no_read_reached_omits_the_key_rather_than_writing_a_missing_value() {
    // The other half of the rule: an undefined *key* disappears, an undefined *entry* stays
    // and writes `.`. A parser distinguishes the two.
    let record = VcfRecord::new(
        region_on(0, 12, 12),
        vec![allele(b"A"), allele(b"T")],
        vec![0.0, 2.0],
        one_sample(vec![0, 9], 0, &[1, 1]),
        pools(&[0, 9]),
        None,
        quality(45.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    let info = info_column(&record);
    assert!(!info.contains("MQREF"), "got {info}");
    // MQALT still has a value; MQDIFF cannot be computed without the reference, so it is `.`.
    assert!(info.contains("MQALT=60.00;MQDIFF=."), "got {info}");
}

#[test]
fn a_record_with_no_alternative_omits_every_per_alternative_key() {
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

    let info = info_column(&record);
    for key in ["AF=", "AC=", "MQALT=", "MQDIFF="] {
        assert!(!info.contains(key), "{key} should be absent, got {info}");
    }
    // What remains is the site's own counts and the tract annotation.
    assert_eq!(info, "AN=0;DP=2;STR;RU=AT;PERIOD=2");
}

#[test]
fn the_penalties_are_absent_where_the_correction_did_not_run() {
    let info = info_column(&a_biallelic_snp());
    assert!(!info.contains("ABPEN"), "got {info}");
    assert!(!info.contains("SPPEN"), "got {info}");
}

#[test]
fn the_uncorrected_quality_is_the_written_one_plus_its_two_penalties() {
    // Spec §6's reason for publishing them: the correction stays recoverable. Both are written
    // at the quality's own precision so the addition is exact in the file's own digits.
    let record = VcfRecord::new(
        region_on(0, 60, 60),
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        one_sample(vec![5, 5], 0, &[0, 1]),
        pools(&[5, 5]),
        None,
        quality(112.5),
        Some(ArtifactPenalties {
            allele_balance: quality(3.5),
            strand_and_read_position: quality(4.0),
        }),
        FilterVerdict::Pass,
        None,
    );

    let info = info_column(&record);
    assert!(info.contains("ABPEN=3.5;SPPEN=4.0"), "got {info}");
    let written = fixed_columns(&record, &contigs());
    let quality: f64 = written
        .split('\t')
        .nth(5)
        .expect("a QUAL column")
        .parse()
        .expect("a numeric quality");
    assert!((quality + 3.5 + 4.0 - 120.0).abs() < 1e-9);
}

#[test]
fn several_alternatives_give_one_entry_per_alternative_in_every_per_alternative_key() {
    let record = VcfRecord::new(
        region_on(0, 77, 77),
        vec![allele(b"A"), allele(b"C"), allele(b"G")],
        vec![2.0, 1.0, 1.0],
        vec![
            SampleColumn {
                call: SampleCall::Called {
                    genotype: Genotype::new(vec![AlleleId(0), AlleleId(1)]),
                    genotype_quality: quality(30.0),
                },
                read_counts: SampleReadCounts::new(vec![4, 4, 0], 0),
            },
            SampleColumn {
                call: SampleCall::Called {
                    genotype: Genotype::new(vec![AlleleId(0), AlleleId(2)]),
                    genotype_quality: quality(30.0),
                },
                read_counts: SampleReadCounts::new(vec![4, 0, 4], 0),
            },
        ],
        pools(&[8, 4, 4]),
        None,
        quality(60.0),
        None,
        FilterVerdict::Pass,
        None,
    );

    let info = info_column(&record);
    // Copies 1.0 and 1.0 of a total 4.0; one copy of each ALT called; AN 4.
    assert!(
        info.starts_with("AF=0.250000,0.250000;AC=1,1;AN=4;DP=16;"),
        "got {info}"
    );
    for key in ["MQALT=", "MQDIFF="] {
        let field = info
            .split(';')
            .find(|field| field.starts_with(key))
            .unwrap_or_else(|| panic!("{key} present"));
        assert_eq!(field.split(',').count(), 2, "{key} has one entry per ALT");
    }
}
