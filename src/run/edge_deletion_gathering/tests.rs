use super::*;
use crate::types::{ContigId, Motif};
use crate::vcf::{FilterVerdict, TractAnnotation};

fn region(start: u64, end: u64) -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(0),
        start: Position(start),
        end: Position(end),
    }
}

fn called(alleles: &[u16], quality: f32) -> SampleCall {
    SampleCall::Called {
        genotype: Genotype::new(alleles.iter().copied().map(AlleleId).collect()),
        genotype_quality: Phred::try_new(quality).expect("a quality"),
    }
}

/// A record from its columns, the mapping qualities pooled over their `AD` at 60 each.
fn record_of(
    at: GenomeRegion,
    alleles: &[&[u8]],
    columns: Vec<SampleColumn>,
    padding_base: Option<PaddingBase>,
    repeat_tract: Option<TractAnnotation>,
) -> VcfRecord {
    let pools = (0..alleles.len())
        .map(|allele| {
            let reads: u64 = columns
                .iter()
                .map(|column| u64::from(column.read_counts.allele_reads()[allele]))
                .sum();
            MapqPool {
                reads,
                mapq_sum: reads * 60,
            }
        })
        .collect();
    VcfRecord::new(
        at,
        alleles.iter().map(|bases| Box::from(*bases)).collect(),
        vec![0.5; alleles.len()],
        columns,
        pools,
        padding_base,
        Phred::try_new(50.0).expect("a quality"),
        None,
        FilterVerdict::Pass,
        repeat_tract,
    )
}

/// A record with each sample's call and one read per allele copy it names.
fn record(
    at: GenomeRegion,
    alleles: &[&[u8]],
    calls: &[(&[u16], f32)],
    repeat_tract: Option<TractAnnotation>,
) -> VcfRecord {
    let columns = calls
        .iter()
        .map(|(genotype, quality)| {
            let mut reads = vec![0_u32; alleles.len()];
            for &allele in *genotype {
                reads[usize::from(allele)] += 1;
            }
            SampleColumn {
                call: called(genotype, *quality),
                read_counts: SampleReadCounts::new(reads, 0),
            }
        })
        .collect();
    record_of(at, alleles, columns, None, repeat_tract)
}

fn widened(motif: &[u8], widening: TractWidening) -> Option<TractAnnotation> {
    Some(TractAnnotation::new(Motif::new(motif).expect("a motif")).widened(widening))
}

/// `GTT` ×3 at 101–109, widened one base right over the flank's `A` at 110, with a deletion
/// of `TTA` taking it (allele 1, reaching one base right) and an expansion (allele 2).
fn a_tract_widened_right(calls: &[(&[u16], f32)]) -> VcfRecord {
    record(
        region(101, 110),
        &[b"GTTGTTGTTA", b"GTTGTTG", b"GTTGTTGTTGTTA"],
        calls,
        widened(
            b"GTT",
            TractWidening {
                left: 0,
                right: 1,
                flank_bases_per_allele: vec![1, 0, 1],
                reach_per_allele: vec![(0, 0), (0, 1), (0, 0)],
                base_before: Some(b'C'),
            },
        ),
    )
}

/// A SNP at the flank's `A`, 110: `A` or `G`.
fn a_snp_at_110(calls: &[(&[u16], f32)]) -> VcfRecord {
    record(region(110, 110), &[b"A", b"G"], calls, None)
}

fn written(gathering: &mut EdgeDeletionGathering, records: Vec<VcfRecord>) -> Vec<VcfRecord> {
    let mut out = Vec::new();
    for record in records {
        out.extend(
            gathering
                .accept(record, Vec::new())
                .into_iter()
                .map(|(r, _)| r),
        );
    }
    out.extend(gathering.finish().into_iter().map(|(r, _)| r));
    out
}

fn genotypes(record: &VcfRecord) -> Vec<Option<Vec<u16>>> {
    record
        .sample_columns()
        .iter()
        .map(|column| {
            column
                .call
                .genotype()
                .map(|genotype| genotype.alleles().iter().map(|a| a.get()).collect())
        })
        .collect()
}

fn alleles(record: &VcfRecord) -> Vec<&[u8]> {
    record.alleles().iter().map(|allele| &**allele).collect()
}

fn no_record_holds_a_star(records: &[VcfRecord]) {
    for record in records {
        assert!(
            record.alleles().iter().all(|allele| &**allele != b"*"),
            "{record:?}"
        );
    }
}

/// **The flank half the mapper wrote as a substitution is the deletion's own, and the merged
/// record is the STR record** (the stuttering benchmark's `seg02581`, 2026-10-11). The sample
/// carries the deletion taking the flank's `A`, and the generic path called that `A` as a `G`
/// from the reads aligned with a mismatch: the deletion's copy takes the `G`, nothing is left
/// for the generic record to add, and one record is written.
#[test]
fn a_flank_half_written_as_a_substitution_is_absorbed_into_the_str_record() {
    let mut gathering = EdgeDeletionGathering::default();
    let tract = a_tract_widened_right(&[(&[1, 2], 30.0)]);
    let out = written(
        &mut gathering,
        vec![tract.clone(), a_snp_at_110(&[(&[0, 1], 40.0)])],
    );
    assert_eq!(out.len(), 1, "one record: {out:?}");
    assert_eq!(alleles(&out[0]), alleles(&tract));
    assert_eq!(genotypes(&out[0]), vec![Some(vec![1, 2])]);
    assert_eq!(
        out[0].sample_columns()[0].call.genotype_quality(),
        Phred::try_new(30.0).ok()
    );
    let counts = gathering.counts();
    assert_eq!(counts.samples_in_agreement, 1);
    assert_eq!(counts.generic_records_merged, 1);
}

/// **Another sample's own SNP becomes an allele of the STR record, and the STR record wins
/// where the generic call disagrees.** Sample 0 carries the deletion and its generic call is
/// `0/1`: the deletion takes the `G`, the expansion pairs with the reference. Sample 1 carries
/// the deletion but the generic path called it `0/0` with `GQ` 50: the STR record wins, and
/// `GQ` is the STR record's 12. Sample 2 has no deletion and a real `G`: its reference STR
/// allele carries the `G` on one copy, a new allele.
#[test]
fn a_real_snp_becomes_an_allele_and_a_disagreeing_call_is_overridden() {
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(
        &mut gathering,
        vec![
            a_tract_widened_right(&[(&[1, 2], 30.0), (&[0, 1], 12.0), (&[0, 0], 40.0)]),
            a_snp_at_110(&[(&[0, 1], 40.0), (&[0, 0], 50.0), (&[0, 1], 45.0)]),
        ],
    );
    no_record_holds_a_star(&out);
    assert_eq!(out.len(), 1);
    let merged = &out[0];
    assert_eq!(merged.region(), region(101, 110));
    assert_eq!(
        alleles(merged),
        vec![
            &b"GTTGTTGTTA"[..],
            b"GTTGTTG",
            b"GTTGTTGTTGTTA",
            b"GTTGTTGTTG"
        ]
    );
    assert_eq!(
        genotypes(merged),
        vec![Some(vec![1, 2]), Some(vec![0, 1]), Some(vec![0, 3])]
    );
    let qualities: Vec<Option<Phred>> = merged
        .sample_columns()
        .iter()
        .map(|column| column.call.genotype_quality())
        .collect();
    assert_eq!(
        qualities,
        vec![
            Phred::try_new(30.0).ok(),
            Phred::try_new(12.0).ok(),
            Phred::try_new(40.0).ok()
        ],
        "the lower GQ where the generic call adds an allele or was overridden"
    );
    assert_eq!(gathering.counts().samples_overridden, 1);
    assert_eq!(gathering.counts().samples_in_agreement, 1);
}

/// **A record widened to the left absorbs the generic record that runs past its start**, and
/// the anchor base that record carried, shared by every allele, is trimmed: the record starts
/// where the STR record did, after a SNP further left that it does not reach.
#[test]
fn a_left_widened_record_absorbs_the_flank_half_and_trims_its_anchor() {
    // Flank ... C(96) A(97) G(98) C(99) | tract `CA` ×3 at 100–105. The STR record is widened
    // two bases left, to 98, over a deletion of `GC` + `CA` (allele 1). The generic path wrote
    // the flank half as `AGC`→`A` at 97–99.
    let tract = record(
        region(98, 105),
        &[b"GCCACACA", b"CACA"],
        &[(&[0, 1], 30.0)],
        widened(
            b"CA",
            TractWidening {
                left: 2,
                right: 0,
                flank_bases_per_allele: vec![2, 0],
                reach_per_allele: vec![(0, 0), (2, 0)],
                base_before: Some(b'A'),
            },
        ),
    );
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(
        &mut gathering,
        vec![
            record(region(96, 96), &[b"C", b"T"], &[(&[0, 1], 40.0)], None),
            record(region(97, 99), &[b"AGC", b"A"], &[(&[0, 1], 40.0)], None),
            tract,
        ],
    );
    no_record_holds_a_star(&out);
    let starts: Vec<u64> = out.iter().map(|r| r.region().start.get()).collect();
    assert_eq!(starts, vec![96, 98], "the SNP, then the merged record");
    assert_eq!(alleles(&out[1]), vec![&b"GCCACACA"[..], b"CACA"]);
    assert_eq!(genotypes(&out[1]), vec![Some(vec![0, 1])]);
    assert_eq!(gathering.counts().generic_records_merged, 1);
}

/// **A sample heterozygous at both records takes the pairing the cohort's settled samples
/// show** (spec §6.2, step 3). Sample 0 is homozygous at the STR and carries the `G` with the
/// reference allele; sample 1 is homozygous for the expansion with the reference flank. Sample
/// 2 is `0/2` at the STR and `A/G` in the flank: the `G` goes with the reference allele, and its
/// `GQ` says how sure that is.
#[test]
fn a_sample_heterozygous_at_both_records_is_paired_from_the_cohort() {
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(
        &mut gathering,
        vec![
            a_tract_widened_right(&[(&[0, 0], 30.0), (&[2, 2], 30.0), (&[0, 2], 30.0)]),
            a_snp_at_110(&[(&[0, 1], 40.0), (&[0, 0], 40.0), (&[0, 1], 40.0)]),
        ],
    );
    assert_eq!(out.len(), 1);
    let merged = &out[0];
    assert_eq!(
        alleles(merged),
        vec![&b"GTTGTTGTTA"[..], b"GTTGTTGTTGTTA", b"GTTGTTGTTG"]
    );
    assert_eq!(genotypes(merged)[2], Some(vec![1, 2]));
    // Pairing (A with allele 0, G with allele 2) weighs 1.5 × 0.5; the other, 1.5 × 2.5.
    let p = (1.5 * 2.5) / (1.5 * 2.5 + 1.5 * 0.5);
    let expected = -10.0 * crate::float::log10(1.0 - (1.0 - 1e-3) * p);
    let gq = merged.sample_columns()[2]
        .call
        .genotype_quality()
        .expect("called")
        .get();
    assert!(
        (f64::from(gq) - expected).abs() < 1e-3,
        "{gq} against {expected}"
    );
    assert_eq!(gathering.counts().samples_paired_from_the_cohort, 1);
}

/// **At one sample the cohort has nothing to offer**: the first pairing by allele order, and a
/// `GQ` of a coin's chance of being wrong.
#[test]
fn at_one_sample_the_pairing_is_the_first_and_its_quality_is_three() {
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(
        &mut gathering,
        vec![
            a_tract_widened_right(&[(&[0, 2], 30.0)]),
            a_snp_at_110(&[(&[0, 1], 40.0)]),
        ],
    );
    let merged = &out[0];
    assert_eq!(
        alleles(merged),
        vec![&b"GTTGTTGTTA"[..], b"GTTGTTGTTGTTG"],
        "the reference allele with the A, the expansion with the G"
    );
    assert_eq!(genotypes(merged), vec![Some(vec![0, 1])]);
    let gq = merged.sample_columns()[0]
        .call
        .genotype_quality()
        .expect("called")
        .get();
    assert!(gq > 2.9 && gq <= 3.02, "{gq}");
}

/// **A generic no-call keeps the STR genotype; an STR no-call is a no-call of the merged
/// record, its reads still written, and a variant the generic record called for it is lost and
/// counted.**
#[test]
fn a_no_call_at_either_record() {
    let tract = record_of(
        region(101, 110),
        &[b"GTTGTTGTTA", b"GTTGTTG"],
        vec![
            SampleColumn {
                call: called(&[0, 1], 30.0),
                read_counts: SampleReadCounts::new(vec![4, 5], 1),
            },
            SampleColumn {
                call: SampleCall::NoCall,
                read_counts: SampleReadCounts::new(vec![3, 0], 0),
            },
        ],
        None,
        widened(
            b"GTT",
            TractWidening {
                left: 0,
                right: 1,
                flank_bases_per_allele: vec![1, 0],
                reach_per_allele: vec![(0, 0), (0, 1)],
                base_before: Some(b'C'),
            },
        ),
    );
    let snp = record_of(
        region(110, 110),
        &[b"A", b"G"],
        vec![
            SampleColumn {
                call: SampleCall::NoCall,
                read_counts: SampleReadCounts::new(vec![0, 0], 0),
            },
            SampleColumn {
                call: called(&[0, 1], 20.0),
                read_counts: SampleReadCounts::new(vec![1, 2], 0),
            },
        ],
        None,
        None,
    );
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(&mut gathering, vec![tract, snp]);
    assert_eq!(out.len(), 1);
    let merged = &out[0];
    assert_eq!(alleles(merged), vec![&b"GTTGTTGTTA"[..], b"GTTGTTG"]);
    assert_eq!(genotypes(merged), vec![Some(vec![0, 1]), None]);
    assert_eq!(
        merged.sample_columns()[0].call.genotype_quality(),
        Phred::try_new(30.0).ok()
    );
    assert_eq!(
        merged.sample_columns()[1].read_counts.allele_reads(),
        &[3, 0]
    );
    assert_eq!(gathering.counts().generic_calls_lost_to_an_str_no_call, 1);
}

/// **An STR allele's reads cannot be split between two merged alleles built on it**, so at a
/// sample naming both they go to the reads no written allele explains; and the two alleles,
/// differing only in the flank, have one `REPCN`.
#[test]
fn reads_of_an_str_allele_split_over_two_alleles_are_unexplained_and_repcn_is_shared() {
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(
        &mut gathering,
        vec![
            a_tract_widened_right(&[(&[0, 0], 30.0)]),
            a_snp_at_110(&[(&[0, 1], 40.0)]),
        ],
    );
    let merged = &out[0];
    assert_eq!(alleles(merged), vec![&b"GTTGTTGTTA"[..], b"GTTGTTGTTG"]);
    let column = &merged.sample_columns()[0];
    assert_eq!(column.read_counts.allele_reads(), &[0, 0]);
    assert_eq!(column.read_counts.depth(), 2);
    let tract = merged.repeat_tract().expect("an STR record");
    assert_eq!(tract.repeat_copies_of_allele(0, &merged.alleles()[0]), 3);
    assert_eq!(tract.repeat_copies_of_allele(1, &merged.alleles()[1]), 3);
}

/// **A haplotype that spells another STR allele's sequence takes that allele's reads and fitted
/// copies with it.** STR allele 1 changes the flank's `A` to `G` itself; a sample called 0/0 at
/// the STR record and `G/G` at the SNP spells the same sequence, so it is written 1/1 — and its
/// reads, which the STR record put on allele 0, go to allele 1, as does allele 0's `AF` share.
#[test]
fn a_haplotype_spelling_another_str_alleles_sequence_takes_its_reads_along() {
    let tract = record(
        region(101, 110),
        &[b"GTTGTTGTTA", b"GTTGTTGTTG"],
        &[(&[1, 1], 30.0), (&[0, 0], 30.0)],
        widened(
            b"GTT",
            TractWidening {
                left: 0,
                right: 1,
                flank_bases_per_allele: vec![1, 1],
                reach_per_allele: vec![(0, 0), (0, 1)],
                base_before: Some(b'C'),
            },
        ),
    );
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(
        &mut gathering,
        vec![tract, a_snp_at_110(&[(&[1, 1], 40.0), (&[1, 1], 40.0)])],
    );
    assert_eq!(out.len(), 1);
    let merged = &out[0];
    assert_eq!(alleles(merged), vec![&b"GTTGTTGTTA"[..], b"GTTGTTGTTG"]);
    assert_eq!(genotypes(merged), vec![Some(vec![1, 1]), Some(vec![1, 1])]);
    let column = &merged.sample_columns()[1];
    assert_eq!(column.read_counts.allele_reads(), &[0, 2]);
    assert_eq!(merged.expected_copies(), &[0.0, 1.0]);
    assert_eq!(merged.allele_mapq()[1].reads, 4);
}

/// **A generic record on each flank of one STR record**: both merge, one after the other, the
/// second paired against the haplotypes the first left.
#[test]
fn generic_records_on_both_flanks_merge_into_one_record() {
    // G(98) C(99) | `CA` ×3 at 100–105 | T(106). Allele 1 deletes `GC` and reaches two bases
    // left; allele 2 deletes the `T` and reaches one right.
    let tract = record(
        region(98, 106),
        &[b"GCCACACAT", b"CACACAT", b"GCCACACA"],
        &[(&[0, 0], 30.0)],
        widened(
            b"CA",
            TractWidening {
                left: 2,
                right: 1,
                flank_bases_per_allele: vec![3, 1, 2],
                reach_per_allele: vec![(0, 0), (2, 0), (0, 1)],
                base_before: Some(b'A'),
            },
        ),
    );
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(
        &mut gathering,
        vec![
            record(region(98, 98), &[b"G", b"A"], &[(&[0, 1], 40.0)], None),
            tract,
            record(region(106, 106), &[b"T", b"C"], &[(&[0, 1], 40.0)], None),
        ],
    );
    no_record_holds_a_star(&out);
    assert_eq!(out.len(), 1, "{out:?}");
    let merged = &out[0];
    assert_eq!(merged.region(), region(98, 106));
    assert_eq!(alleles(merged), vec![&b"GCCACACAT"[..], b"ACCACACAC"]);
    assert_eq!(genotypes(merged), vec![Some(vec![0, 1])]);
    assert_eq!(gathering.counts().generic_records_merged, 2);
    assert_eq!(gathering.counts().samples_paired_from_the_cohort, 1);
}

/// **An allele that comes out empty is padded with the base the widening kept** — here the STR
/// record's deletion of everything it spans, at a sample whose other copy carries a flank SNP.
#[test]
fn an_empty_allele_is_padded_with_the_base_before_the_widened_span() {
    let tract = record_of(
        region(98, 105),
        &[b"GCCACACA", b""],
        vec![SampleColumn {
            call: called(&[0, 1], 30.0),
            read_counts: SampleReadCounts::new(vec![1, 1], 0),
        }],
        Some(PaddingBase::Left(b'A')),
        widened(
            b"CA",
            TractWidening {
                left: 2,
                right: 0,
                flank_bases_per_allele: vec![2, 0],
                reach_per_allele: vec![(0, 0), (2, 0)],
                base_before: Some(b'A'),
            },
        ),
    );
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(
        &mut gathering,
        vec![
            record(region(99, 99), &[b"C", b"T"], &[(&[0, 1], 40.0)], None),
            tract,
        ],
    );
    assert_eq!(out.len(), 1, "{out:?}");
    let merged = &out[0];
    assert_eq!(alleles(merged), vec![&b"GCCACACA"[..], b"", b"GTCACACA"]);
    assert_eq!(genotypes(merged), vec![Some(vec![1, 2])]);
    assert_eq!(merged.padding_base(), Some(PaddingBase::Left(b'A')));
    assert_eq!(gathering.counts().samples_overridden, 1);
}

/// **A widened record nothing overlaps is written as it came.**
#[test]
fn a_widened_record_alone_is_unchanged() {
    let tract = a_tract_widened_right(&[(&[1, 2], 30.0)]);
    let mut gathering = EdgeDeletionGathering::default();
    let out = written(&mut gathering, vec![tract.clone()]);
    assert_eq!(out, vec![tract]);
    assert_eq!(gathering.counts(), EdgeDeletionGatheringCounts::default());
}
