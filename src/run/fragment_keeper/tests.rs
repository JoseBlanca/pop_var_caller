use super::*;
use crate::read_phase::record_reads::{RecordReadAlleles, SampleReadAllele};
use crate::types::{AlleleId, GenomeRegion, Genotype, Phred, Position};
use crate::vcf::{
    FilterVerdict, MapqPool, PaddingBase, SampleCall, SampleColumn, SampleReadCounts,
};

/// A biallelic SNP record over two samples at `position`, carrying `reads` as `(sample, read,
/// allele)`.
fn snp_at(contig: u32, position: u64, reads: &[(usize, ChainId, u32)]) -> VcfRecord {
    record_at(contig, position, None, reads)
}

/// The same, or with a padding base a one-base deletion, which `PaddingBase::Left` writes one base
/// before `position`.
fn record_at(
    contig: u32,
    position: u64,
    padding: Option<PaddingBase>,
    reads: &[(usize, ChainId, u32)],
) -> VcfRecord {
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
            // A padded record is a deletion: one allele empty, which is what padding is for.
            if padding.is_some() {
                b"".to_vec()
            } else {
                b"T".to_vec()
            }
            .into_boxed_slice(),
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
        padding,
        Phred::try_new(50.0).expect("a quality"),
        None,
        FilterVerdict::Pass,
        None,
    );
    let mut entries: Vec<SampleReadAllele> = reads
        .iter()
        .map(|&(sample, read, allele)| SampleReadAllele {
            sample: sample as u32,
            read_allele: ReadAllele {
                read,
                allele,
                error_phred: 25,
            },
        })
        .collect();
    entries.sort_unstable();
    record.replace_read_alleles(RecordReadAlleles::from_merge(entries, Some));
    record
}

/// One released record as `(position, [(sample, fragment, allele)])`.
type Released = (u64, Vec<(u32, u32, u16)>);

/// Every record through a keeper, as `(position, [(sample, fragment, allele)])`.
fn kept(records: Vec<VcfRecord>) -> (Vec<Released>, KeeperCounts) {
    let mut keeper = FragmentKeeper::default();
    let mut out = Vec::new();
    for record in records {
        keeper.accept(record, (), &mut out);
    }
    keeper.finish(&mut out);
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
            read_alleles_over_the_cap: 0,
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
    let mut released = Vec::new();
    keeper.accept(snp_at(0, 100, &[]), (), &mut released);
    keeper.accept(snp_at(0, 1_100, &[]), (), &mut released);
    assert!(released.is_empty());
    keeper.accept(snp_at(0, 1_101, &[]), (), &mut released);
    assert_eq!(released.len(), 1);
    assert_eq!(released[0].0.region().start.get(), 100);
    let mut rest = Vec::new();
    keeper.finish(&mut rest);
    let rest: Vec<u64> = rest
        .into_iter()
        .map(|(record, (), _)| record.region().start.get())
        .collect();
    assert_eq!(rest, vec![1_100, 1_101]);
}

/// **A read allele exactly `R` past its read's first is kept**: the rule drops only one more than
/// `R` past it. Kills `>` written `>=` in the reach test.
#[test]
fn a_read_exactly_the_reach_past_its_first_is_kept() {
    let (shown, counts) = kept(vec![
        snp_at(0, 100, &[(0, 1, 1)]),
        snp_at(0, 100 + READ_REACH, &[(0, 1, 0)]),
    ]);
    assert_eq!(
        shown,
        vec![(100, vec![(0, 0, 1)]), (100 + READ_REACH, vec![(0, 0, 0)])]
    );
    assert_eq!(counts.read_alleles_beyond_reach, 0);
}

/// **A read is not forgotten while its first lies exactly `R` behind the oldest record held**:
/// read 4 at 100, an empty record at 1,100, then read 4 at 2,100 and 2,200. When 2,100 arrives the
/// record at 1,100 is the oldest held, exactly `R` past read 4's first, so read 4 is remembered
/// and its read allele at 2,100, more than `R` past its first, is dropped. (By 2,200 the oldest
/// held is 2,100 and read 4 is forgotten, so there it starts a new fragment seen once.) Forgotten
/// one record early, read 4 would be a new fragment at 2,100 and 2,200 both, and kept at both.
/// Kills `<` written `<=` in the forgetting test.
#[test]
fn a_read_whose_first_is_exactly_the_reach_behind_the_oldest_held_is_remembered() {
    let (shown, counts) = kept(vec![
        snp_at(0, 100, &[(0, 4, 1)]),
        snp_at(0, 1_100, &[]),
        snp_at(0, 2_100, &[(0, 4, 0)]),
        snp_at(0, 2_200, &[(0, 4, 0)]),
    ]);
    assert_eq!(shown[2], (2_100, vec![]));
    assert_eq!(shown[3], (2_200, vec![]));
    assert_eq!(counts.read_alleles_beyond_reach, 1);
}

/// **Records are ordered by written position, as the writer orders them** (review M1): a generic
/// record padded on the left, its span starting at 101 and written at 100, then an STR record
/// starting at 100, is the writer's one legal tie. Keyed on region starts the keeper saw 101 then
/// 100 and panicked; keyed on written positions it sees 100 twice, and read 3 links the two.
#[test]
fn a_record_after_a_left_padded_one_at_the_same_written_position_is_accepted() {
    let (shown, _) = kept(vec![
        record_at(0, 101, Some(PaddingBase::Left(b'C')), &[(0, 3, 1)]),
        snp_at(0, 100, &[(0, 3, 0)]),
    ]);
    assert_eq!(shown, vec![(101, vec![(0, 0, 1)]), (100, vec![(0, 0, 0)])]);
}

/// **A record written before the one handed in before it is a wiring defect**, and stops the run.
#[test]
#[should_panic(expected = "records out of order")]
fn a_record_out_of_order_panics() {
    let mut keeper: FragmentKeeper<()> = FragmentKeeper::default();
    let mut released = Vec::new();
    keeper.accept(snp_at(0, 200, &[]), (), &mut released);
    keeper.accept(snp_at(0, 199, &[]), (), &mut released);
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

/// **What the cap removed at capture reaches the keeper's counts**, summed over records (review
/// Mi8): the cap's work is otherwise invisible in the log at the deep end of the depth range.
#[test]
fn the_read_alleles_the_cap_removed_are_counted() {
    let mut keeper: FragmentKeeper<()> = FragmentKeeper::default();
    let mut released = Vec::new();
    for (position, over_the_cap) in [(100, 3), (200, 5)] {
        let mut record = snp_at(0, position, &[]);
        let read_alleles = record.replace_read_alleles(Default::default());
        record.replace_read_alleles(read_alleles.with_dropped_by_the_cap(over_the_cap));
        keeper.accept(record, (), &mut released);
    }
    assert_eq!(keeper.counts().read_alleles_over_the_cap, 8);
}

/// **A keeper for a run whose output reads no fragments holds nothing and keeps nothing**: each
/// record leaves the moment it arrives, in the order it came, its read alleles taken off it and no
/// fragment beside it — the same three records the keeping keeper holds and links in the first
/// test above. Its counts stay at zero, so the driver has nothing to log.
#[test]
fn a_passing_through_keeper_releases_each_record_at_once_with_no_fragments() {
    let mut keeper: FragmentKeeper<u32> = FragmentKeeper::for_a_run(false);
    assert!(!keeper.keeps_fragments());
    let mut released = Vec::new();
    for (payload, record) in [
        snp_at(0, 100, &[(0, 1, 1), (0, 2, 0), (0, 3, 1)]),
        snp_at(0, 600, &[(0, 1, 0), (0, 3, 0)]),
        snp_at(0, 1_200, &[(0, 3, 1)]),
    ]
    .into_iter()
    .enumerate()
    {
        keeper.accept(record, payload as u32, &mut released);
        assert_eq!(
            released.len(),
            payload + 1,
            "the record leaves as it arrives, with nothing held behind it"
        );
    }
    keeper.finish(&mut released);
    assert_eq!(released.len(), 3, "finishing releases nothing more");
    for (index, (record, payload, fragments)) in released.iter().enumerate() {
        assert_eq!(*payload as usize, index, "in the order they came");
        assert!(
            fragments.is_empty(),
            "no fragment at {}",
            record.region().start.get()
        );
        assert!(
            record.read_alleles().is_empty(),
            "the read alleles are taken off, as the keeping keeper takes them"
        );
    }
    assert_eq!(keeper.counts(), KeeperCounts::default());

    let keeping: FragmentKeeper<()> = FragmentKeeper::for_a_run(true);
    assert!(
        keeping.keeps_fragments(),
        "a run that reads fragments keeps them"
    );
}

/// **Sharing the samples among the threads changes nothing a sample is given.** A 300-sample
/// stream spans five groups of [`SAMPLES_A_TASK`], the last one short; each sample's fragments out
/// of it must be the ones a keeper handed that sample's reads alone gives, and the counts their
/// sums. Each sample's reads differ — which records it has reads at, how many, how far they
/// reach — so a read allele handed to the wrong sample's state, or a group boundary off by one,
/// shows as a sample whose fragments differ from its own run.
#[test]
fn each_sample_is_kept_as_it_would_be_alone_across_groups_of_samples() {
    const SAMPLES: usize = 300;
    assert!(SAMPLES > 4 * SAMPLES_A_TASK && SAMPLES % SAMPLES_A_TASK != 0);
    // A cheap, fixed scramble, so the stream is the same on every run.
    let mix = |a: u64, b: u64| {
        (a.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ b).wrapping_mul(0xBF58_476D_1CE4_E5B9) >> 33
    };
    // (contig, position) of each record: two contigs, records 37 to 400 bases apart, so reads
    // reach some neighbours and not others, and the contig change resets every sample.
    let mut places = Vec::new();
    for contig in 0..2u32 {
        let mut position = 1;
        for step in 0..120u64 {
            position += 37 + mix(u64::from(contig), step) % 364;
            places.push((contig, position));
        }
    }
    // Sample `s` at record `r`: none, or up to three reads, each a read the sample started at one
    // of its last few records, so reads recur across records within and beyond the reach.
    let reads_of = |sample: usize, record: usize| -> Vec<(ChainId, u32)> {
        let roll = mix(sample as u64, record as u64);
        if roll % 4 == 0 {
            return Vec::new();
        }
        let count = 1 + roll % 3;
        (0..count)
            .map(|k| {
                let started = record.saturating_sub((mix(roll, k) % 4) as usize);
                let read: ChainId = sample as u64 * 1_000_000 + started as u64 * 4 + k;
                (read, (mix(read, record as u64) % 2) as u32)
            })
            .collect()
    };
    let stream = |samples: &[usize]| -> Vec<VcfRecord> {
        places
            .iter()
            .enumerate()
            .map(|(record, &(contig, position))| {
                let reads: Vec<(usize, ChainId, u32)> = samples
                    .iter()
                    .enumerate()
                    .flat_map(|(slot, &sample)| {
                        reads_of(sample, record)
                            .into_iter()
                            .map(move |(read, allele)| (slot, read, allele))
                    })
                    .collect();
                snp_at(contig, position, &reads)
            })
            .collect()
    };

    let all: Vec<usize> = (0..SAMPLES).collect();
    let (together, counts) = kept(stream(&all));
    let mut summed = KeeperCounts::default();
    for sample in 0..SAMPLES {
        let (alone, alone_counts) = kept(stream(&[sample]));
        let picked: Vec<Released> = together
            .iter()
            .map(|(position, entries)| {
                (
                    *position,
                    entries
                        .iter()
                        .filter(|entry| entry.0 as usize == sample)
                        .map(|&(_, fragment, allele)| (0, fragment, allele))
                        .collect(),
                )
            })
            .collect();
        assert_eq!(picked, alone, "sample {sample}");
        summed.read_alleles_in += alone_counts.read_alleles_in;
        summed.read_alleles_kept += alone_counts.read_alleles_kept;
        summed.read_alleles_alone += alone_counts.read_alleles_alone;
        summed.read_alleles_beyond_reach += alone_counts.read_alleles_beyond_reach;
    }
    assert_eq!(counts, summed);
    // The stream exercises every outcome, or the comparison above proves less than it says.
    assert!(counts.read_alleles_kept > 0);
    assert!(counts.read_alleles_alone > 0);
    assert!(counts.read_alleles_beyond_reach > 0);
}
