use super::*;
use crate::read_phase::record_reads::RecordReadAlleles;
use crate::types::{AlleleId, GenomeRegion, Genotype, Phred, Position};
use crate::vcf::{FilterVerdict, MapqPool, SampleCall, SampleColumn, SampleReadCounts};

/// A biallelic SNP record over two samples at `position`, carrying `reads` as `(sample, read,
/// allele)`.
fn snp_at(contig: u32, position: u64, reads: &[(usize, ChainId, u32)]) -> VcfRecord {
    let column = || SampleColumn {
        genotype_rows: None,
        call: SampleCall::Called {
            genotype: Genotype::new(vec![AlleleId(0), AlleleId(1)]),
            genotype_quality: Phred::try_new(30.0).expect("a quality"),
        },
        read_counts: SampleReadCounts::new(vec![5, 5], 0),
    };
    let mut record = VcfRecord::new(
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(position),
            end: Position(position),
        },
        vec![
            b"A".to_vec().into_boxed_slice(),
            b"T".to_vec().into_boxed_slice(),
        ],
        vec![2.0, 2.0],
        vec![column(), column()],
        vec![
            MapqPool {
                reads: 10,
                mapq_sum: 600,
            };
            2
        ],
        None,
        Phred::try_new(50.0).expect("a quality"),
        None,
        FilterVerdict::Pass,
        None,
    );
    let mut by_sample: Vec<(usize, Vec<ReadAllele>)> = Vec::new();
    for &(sample, read, allele) in reads {
        if by_sample.last().is_none_or(|(last, _)| *last != sample) {
            by_sample.push((sample, Vec::new()));
        }
        by_sample
            .last_mut()
            .expect("just pushed")
            .1
            .push(ReadAllele {
                read,
                allele,
                error_phred: 25,
            });
    }
    record.replace_read_alleles(RecordReadAlleles::from_merge(by_sample, Some));
    record
}

/// One released record as `(position, [(sample, fragment, allele)])`.
type Released = (u64, Vec<(u32, u32, u16)>);

/// Every record through a keeper, as `(position, [(sample, fragment, allele)])`.
fn kept(records: Vec<VcfRecord>) -> (Vec<Released>, KeeperCounts) {
    let mut keeper = FragmentKeeper::default();
    let mut out = Vec::new();
    for record in records {
        out.extend(keeper.accept(record, ()));
    }
    out.extend(keeper.finish());
    let shown = out
        .into_iter()
        .map(|(record, (), fragments)| {
            (
                record.region().start.get(),
                fragments
                    .entries()
                    .iter()
                    .map(|entry| (entry.sample, entry.fragment, entry.allele))
                    .collect(),
            )
        })
        .collect();
    (shown, keeper.counts())
}

/// **A read seen at two records within reach is kept at both; one seen once is dropped; one more
/// than `R` past its read's first is dropped and counted** (§2.4). Read 1 is at 100 and 600, read
/// 2 at 100 only, read 3 at 100, 600 and 1,200 — the last 1,100 bases past its first.
#[test]
fn a_read_at_two_records_is_kept_one_at_one_is_not_and_one_beyond_reach_is_dropped() {
    let (shown, counts) = kept(vec![
        snp_at(0, 100, &[(0, 1, 1), (0, 2, 0), (0, 3, 1)]),
        snp_at(0, 600, &[(0, 1, 0), (0, 3, 0)]),
        snp_at(0, 1_200, &[(0, 3, 1)]),
    ]);
    assert_eq!(
        shown,
        vec![
            (100, vec![(0, 0, 1), (0, 2, 1)]),
            (600, vec![(0, 0, 0), (0, 2, 0)]),
            (1_200, vec![]),
        ]
    );
    assert_eq!(
        counts,
        KeeperCounts {
            read_alleles_in: 6,
            read_alleles_kept: 4,
            read_alleles_alone: 1,
            read_alleles_beyond_reach: 1,
        }
    );
}

/// **A record is held until a record more than `R` past it arrives, and records leave in the order
/// they came**: the record at 100 comes out when 1,101 arrives and not when 1,100 does.
#[test]
fn a_record_is_held_until_the_stream_is_more_than_the_reach_past_it() {
    let mut keeper: FragmentKeeper<()> = FragmentKeeper::default();
    assert!(keeper.accept(snp_at(0, 100, &[]), ()).is_empty());
    assert!(keeper.accept(snp_at(0, 1_100, &[]), ()).is_empty());
    let released = keeper.accept(snp_at(0, 1_101, &[]), ());
    assert_eq!(released.len(), 1);
    assert_eq!(released[0].0.region().start.get(), 100);
    let rest: Vec<u64> = keeper
        .finish()
        .into_iter()
        .map(|(record, (), _)| record.region().start.get())
        .collect();
    assert_eq!(rest, vec![1_100, 1_101]);
}

/// **A contig change releases everything and numbers fragments from 0 again**; and **a chain id is
/// per sample**, so read 7 of sample 0 and read 7 of sample 1 are two fragments, each sample
/// numbering its own.
#[test]
fn fragments_are_numbered_per_sample_and_restart_on_a_new_contig() {
    let (shown, _) = kept(vec![
        snp_at(0, 100, &[(0, 7, 1), (1, 7, 0), (1, 9, 1)]),
        snp_at(0, 150, &[(0, 7, 0), (1, 7, 1), (1, 9, 1)]),
        snp_at(1, 100, &[(1, 9, 1)]),
        snp_at(1, 120, &[(1, 9, 0)]),
    ]);
    assert_eq!(
        shown,
        vec![
            (100, vec![(0, 0, 1), (1, 0, 0), (1, 1, 1)]),
            (150, vec![(0, 0, 0), (1, 0, 1), (1, 1, 1)]),
            // Read 9 again on the next contig: fragment 0 there, the numbering having restarted.
            (100, vec![(1, 0, 1)]),
            (120, vec![(1, 0, 0)]),
        ]
    );
}

/// **A read met again after it was forgotten is a new fragment**: read 4 at 100 and 200, then —
/// its first record now more than `R` behind everything held — at 3,000 and 3,100. Both pairs are
/// kept, as two fragments; nothing joins the two halves.
#[test]
fn a_read_met_again_after_it_was_forgotten_is_a_new_fragment() {
    let (shown, counts) = kept(vec![
        snp_at(0, 100, &[(0, 4, 1)]),
        snp_at(0, 200, &[(0, 4, 1)]),
        snp_at(0, 3_000, &[(0, 4, 0)]),
        snp_at(0, 3_100, &[(0, 4, 0)]),
    ]);
    assert_eq!(
        shown,
        vec![
            (100, vec![(0, 0, 1)]),
            (200, vec![(0, 0, 1)]),
            (3_000, vec![(0, 1, 0)]),
            (3_100, vec![(0, 1, 0)]),
        ]
    );
    assert_eq!(counts.read_alleles_beyond_reach, 0);
}
