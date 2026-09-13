//! What the patch must not do, mostly: move a column the verdict does not touch.
//!
//! The standing oracle for the whole filter is that a run at a target no record reaches, with
//! the two `INFO` keys stripped, is byte for byte the run with the filter off. That holds only
//! if the eight columns before `FORMAT` and everything after it come back exactly as they went
//! in, so most of what follows compares columns rather than checking a substring is present.

use super::{LinePatchError, rewrite_filter_and_info};

/// A record as ng's encoder writes it: eight fixed columns, `FORMAT`, then two samples.
fn a_record_line() -> Vec<u8> {
    b"SL4.0ch04\t1000000\t.\tA\tG\t42.5\tPASS\tAF=0.5;AC=1\tGT:GQ:AD\t0/1:30:5,5\t0/0:99:8,0"
        .to_vec()
}

/// The columns of a line, so a test can say which ones moved.
fn columns(line: &[u8]) -> Vec<String> {
    line.split(|&byte| byte == b'\t')
        .map(|column| String::from_utf8(column.to_vec()).expect("the fixtures are text"))
        .collect()
}

/// Patch, and hand back the columns.
fn patched_columns(
    line: &[u8],
    filter_to_add: Option<&str>,
    info_to_add: &[String],
) -> Vec<String> {
    let patched =
        rewrite_filter_and_info(line, filter_to_add, info_to_add).expect("a well-formed line");
    columns(&patched)
}

/// The two fields pass three appends to a scored record.
fn the_two_info_fields() -> Vec<String> {
    vec![
        "PARALOG_LR=4.2117".to_string(),
        "PARALOG_POST=0.9931".to_string(),
    ]
}

#[test]
fn a_line_with_nothing_to_add_comes_back_byte_for_byte() {
    // The one that matters most: an unscored record goes out unchanged, and it goes out
    // through the same splice as every other call rather than around it.
    let line = a_record_line();

    let patched = rewrite_filter_and_info(&line, None, &[]).expect("a well-formed line");

    assert_eq!(patched, line);
}

#[test]
fn only_the_filter_and_info_columns_ever_move() {
    let line = a_record_line();
    let before = columns(&line);

    let after = patched_columns(&line, Some("hiddenParalog"), &the_two_info_fields());

    assert_eq!(after.len(), before.len(), "a column was added or lost");
    for (index, (was, now)) in before.iter().zip(&after).enumerate() {
        if index == 6 || index == 7 {
            continue;
        }
        assert_eq!(
            was, now,
            "column {index} moved, and the verdict may only touch seven and eight"
        );
    }
}

#[test]
fn a_passing_record_carries_the_filter_in_place_of_pass() {
    let after = patched_columns(&a_record_line(), Some("hiddenParalog"), &[]);

    assert_eq!(after[6], "hiddenParalog");
}

#[test]
fn a_record_that_already_carries_a_filter_keeps_it_and_gains_the_new_one() {
    // Spec §6 trap 7: a record the calling loop marked `EMNoConv` and the filter then flags
    // carries both, `;`-separated. Replacing would erase what the loop said.
    let line = b"SL4.0ch04\t1000000\t.\tA\tG\t12.0\tEMNoConv\tAF=0.5\tGT\t0/1\t0/0".to_vec();

    let after = patched_columns(&line, Some("hiddenParalog"), &[]);

    assert_eq!(after[6], "EMNoConv;hiddenParalog");
}

#[test]
fn a_missing_filter_is_replaced_rather_than_joined_to_though_ng_never_writes_one() {
    // `.;hiddenParalog` would say "no filter was applied, and this one was".
    //
    // **Defensive.** `FilterVerdict` has five values and all are words, so this column is never
    // `.` on a line ng emits — `the_encoder_writes_neither_a_missing_nor_an_empty_filter_or_info`
    // says so against the encoder. The *reachable* half of spec §6 trap 7, a real verdict being
    // joined to, is pinned in `encoded_lines`.
    let line = b"SL4.0ch04\t1000000\t.\tA\tG\t12.0\t.\tAF=0.5\tGT\t0/1".to_vec();

    let after = patched_columns(&line, Some("hiddenParalog"), &[]);

    assert_eq!(after[6], "hiddenParalog");
}

#[test]
fn the_info_fields_are_appended_in_the_order_they_are_given() {
    let after = patched_columns(&a_record_line(), None, &the_two_info_fields());

    assert_eq!(
        after[7],
        "AF=0.5;AC=1;PARALOG_LR=4.2117;PARALOG_POST=0.9931"
    );
}

#[test]
fn an_info_that_says_nothing_is_replaced_rather_than_appended_to_though_ng_never_writes_one() {
    // Both spellings: `.`, which is VCF's missing, and an empty column.
    //
    // **Defensive.** `info_column` pushes `AN=` and `DP=` before any conditional field, so the
    // thinnest `INFO` ng writes is `AN=0;DP=0` — never `.`, never empty. Appending to either
    // would put a leading `;` on the column, so the rule is kept for a line from anywhere else.
    let with_a_dot = b"SL4.0ch04\t1000000\t.\tA\t.\t0\tPASS\t.\tGT\t./.".to_vec();
    let with_nothing = b"SL4.0ch04\t1000000\t.\tA\t.\t0\tPASS\t\tGT\t./.".to_vec();

    for line in [with_a_dot, with_nothing] {
        let after = patched_columns(&line, None, &the_two_info_fields());
        assert_eq!(after[7], "PARALOG_LR=4.2117;PARALOG_POST=0.9931");
    }
}

#[test]
fn the_filter_and_the_info_are_independent() {
    let line = a_record_line();

    let filter_only = patched_columns(&line, Some("hiddenParalog"), &[]);
    let info_only = patched_columns(&line, None, &the_two_info_fields());

    assert_eq!(
        filter_only[7], "AF=0.5;AC=1",
        "the INFO moved on a filter-only patch"
    );
    assert_eq!(
        info_only[6], "PASS",
        "the FILTER moved on an INFO-only patch"
    );
}

#[test]
fn a_record_with_one_sample_and_one_with_a_thousand_patch_the_same_two_columns() {
    // The plan's case: the split stops after the eighth tab, so the sample columns are one
    // piece however many of them there are.
    let head = "SL4.0ch04\t1000000\t.\tA\tG\t42.5\tPASS\tAF=0.5\tGT:AD";
    let one = format!("{head}\t0/1:5,5");
    let thousand = format!("{head}{}", "\t0/1:5,5".repeat(1000));

    for line in [one.as_bytes(), thousand.as_bytes()] {
        let before = columns(line);
        let after = patched_columns(line, Some("hiddenParalog"), &the_two_info_fields());

        assert_eq!(after.len(), before.len());
        assert_eq!(after[6], "hiddenParalog");
        assert_eq!(after[7], "AF=0.5;PARALOG_LR=4.2117;PARALOG_POST=0.9931");
        assert_eq!(
            after[8..],
            before[8..],
            "a sample column moved on a {}-sample record",
            before.len() - 9
        );
    }
}

#[test]
fn a_tab_run_inside_the_sample_columns_survives_untouched() {
    // The columns after the eighth tab are never separated, so whatever is in them — including
    // adjacent tabs from a sample with an empty column — comes back as it went in.
    let line = b"SL4.0ch04\t1000000\t.\tA\tG\t42.5\tPASS\tAF=0.5\tGT\t\t0/1".to_vec();

    let patched = rewrite_filter_and_info(&line, None, &[]).expect("a well-formed line");

    assert_eq!(patched, line);
}

#[test]
fn a_line_with_too_few_columns_is_refused_rather_than_patched_in_the_wrong_place() {
    // Seven columns: `FILTER` is there but `INFO` is not, so column seven is the last thing on
    // the line and patching it would append to a record's `FILTER` with no `INFO` behind it.
    let line = b"SL4.0ch04\t1000000\t.\tA\tG\t42.5\tPASS".to_vec();

    match rewrite_filter_and_info(&line, Some("hiddenParalog"), &[]) {
        Err(LinePatchError::TooFewColumns { columns }) => assert_eq!(columns, 7),
        other => panic!("expected the line to be refused, got {other:?}"),
    }
}

#[test]
fn a_line_with_exactly_the_eight_fixed_columns_is_refused() {
    // Eight columns and nothing after them is not a record ng writes — every record has a
    // `FORMAT` and at least one sample — and the split cannot tell it from a truncated line.
    let line = b"SL4.0ch04\t1000000\t.\tA\tG\t42.5\tPASS\tAF=0.5".to_vec();

    assert!(matches!(
        rewrite_filter_and_info(&line, None, &[]),
        Err(LinePatchError::TooFewColumns { columns: 8 })
    ));
}

#[test]
fn an_empty_line_is_refused() {
    assert!(matches!(
        rewrite_filter_and_info(b"", None, &[]),
        Err(LinePatchError::TooFewColumns { columns: 1 })
    ));
}

#[test]
fn an_empty_filter_is_replaced_rather_than_joined_to_though_ng_never_writes_one() {
    // Defensive, like the `.` case above: `FilterVerdict`'s five values are all words, so ng's
    // encoder never leaves this column empty. Joining to it would put a leading `;` on the
    // column.
    let line = b"SL4.0ch04\t1000000\t.\tA\tG\t12.0\t\tAF=0.5\tGT\t0/1".to_vec();

    let after = patched_columns(&line, Some("hiddenParalog"), &[]);

    assert_eq!(after[6], "hiddenParalog");
}

#[test]
fn a_missing_info_with_nothing_to_add_comes_back_byte_for_byte() {
    // The identity law at the one `INFO` spelling `a_line_with_nothing_to_add_…` does not use.
    // Appending nothing must not turn `.` into an empty column.
    let line = b"SL4.0ch04\t1000000\t.\tA\tG\t12.0\tPASS\t.\tGT\t0/1".to_vec();

    let patched = rewrite_filter_and_info(&line, None, &[]).expect("a well-formed line");

    assert_eq!(patched, line);
}

#[test]
fn patching_a_line_twice_adds_the_verdict_twice() {
    // Pinned, not fixed: one call per line is the design, and this says what a second call
    // costs so a pass-three retry does not discover it in a cohort VCF.
    let once = rewrite_filter_and_info(
        &a_record_line(),
        Some("hiddenParalog"),
        &["PARALOG_LR=4.2117".to_string()],
    )
    .expect("a well-formed line");

    let twice = patched_columns(
        &once,
        Some("hiddenParalog"),
        &["PARALOG_LR=4.2117".to_string()],
    );

    assert_eq!(twice[6], "hiddenParalog;hiddenParalog");
    assert_eq!(twice[7], "AF=0.5;AC=1;PARALOG_LR=4.2117;PARALOG_LR=4.2117");
}

#[test]
fn the_capacity_is_never_an_under_estimate() {
    // `added_byte_count` sizes the buffer and nothing links it to what the two append functions
    // write. An under-estimate costs a realloc and a copy of a line that runs to tens of
    // kilobytes at three thousand samples, and no other test can see it.
    let heads = [
        "SL4.0ch04\t1000000\t.\tA\tG\t42.5\tPASS\tAF=0.5;AC=1",
        "SL4.0ch04\t1000000\t.\tA\tG\t42.5\t.\t.",
        "SL4.0ch04\t1000000\t.\tA\tG\t42.5\tEMNoConv\t",
        "SL4.0ch04\t1000000\t.\tA\tG\t42.5\t\tAF=0.5",
    ];
    let filters = [None, Some("hiddenParalog"), Some("x")];
    let infos: [Vec<String>; 3] = [
        vec![],
        vec!["PARALOG_LR=4.2117".to_string()],
        the_two_info_fields(),
    ];

    for head in heads {
        let line = format!("{head}\tGT:AD\t0/1:5,5\t0/0:8,0").into_bytes();
        for filter in filters {
            for info in &infos {
                let capacity = line.len() + super::added_byte_count(filter, info);
                let patched =
                    rewrite_filter_and_info(&line, filter, info).expect("a well-formed line");
                assert!(
                    patched.len() <= capacity,
                    "under-estimate: {} bytes written into a capacity of {capacity}",
                    patched.len()
                );
            }
        }
    }
}

/// **Lines the encoder actually made**, rather than hand-typed ones.
///
/// Everything above patches a fixture a person wrote, and a person writing one gets it slightly
/// wrong: the head of this file spells `AF` at two decimals where
/// [`encode`](crate::vcf::encode) fixes six, and a `FORMAT` of `GT:GQ:AD` where the encoder
/// writes `GT:GQ:DP:AD`. So those tests say the splice is self-consistent and say nothing about
/// the splice against the real thing.
///
/// **This module is the oracle the step exists to protect** (spec §10): a run with the filter on
/// at a target no record reaches, with the two `INFO` keys stripped, must equal the filter-off
/// run byte for byte. Here that is checked one record at a time — build a `VcfRecord`, encode it,
/// patch it with nothing to add, and require the bytes back — over every shape ng emits.
mod encoded_lines {
    use super::{columns, the_two_info_fields};
    use crate::calling::quality::artifact_correction::ArtifactPenalties;
    use crate::run::paralog_filter::rewrite_filter_and_info;
    use crate::types::{
        AlleleId, ContigId, GenomeRegion, Genotype, Motif, Phred, Ploidy, Position,
    };
    use crate::vcf::encode::record_line;
    use crate::vcf::{
        FilterVerdict, HeaderContig, MapqPool, PaddingBase, SampleCall, SampleColumn,
        SampleReadCounts, TractAnnotation, VcfRecord,
    };

    /// Two contigs, so a shape can sit somewhere other than the first.
    fn contigs() -> Vec<HeaderContig> {
        vec![
            HeaderContig {
                name: "chr1".to_string(),
                length: 1_000_000,
                md5: None,
            },
            HeaderContig {
                name: "chr2".to_string(),
                length: 900_000,
                md5: None,
            },
        ]
    }

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("two copies")
    }

    fn quality(phred: f32) -> Phred {
        Phred::try_new(phred).expect("a quality")
    }

    fn allele(bases: &[u8]) -> Box<[u8]> {
        bases.to_vec().into_boxed_slice()
    }

    fn pools(columns: &[SampleColumn], alleles: usize) -> Vec<MapqPool> {
        (0..alleles)
            .map(|index| {
                let reads: u64 = columns
                    .iter()
                    .map(|column| u64::from(column.read_counts.allele_reads()[index]))
                    .sum();
                MapqPool {
                    reads,
                    mapq_sum: reads * 60,
                }
            })
            .collect()
    }

    fn called(alleles: &[u16], counts: Vec<u32>) -> SampleColumn {
        SampleColumn {
            call: SampleCall::Called {
                genotype: Genotype::new(alleles.iter().map(|id| AlleleId(*id)).collect()),
                genotype_quality: quality(30.0),
            },
            read_counts: SampleReadCounts::new(counts, 0),
        }
    }

    fn no_call(counts: Vec<u32>) -> SampleColumn {
        SampleColumn {
            call: SampleCall::NoCall,
            read_counts: SampleReadCounts::new(counts, 0),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn record(
        contig: u32,
        start: u64,
        end: u64,
        alleles: Vec<Box<[u8]>>,
        columns: Vec<SampleColumn>,
        padding: Option<PaddingBase>,
        filter: FilterVerdict,
        tract: Option<TractAnnotation>,
        penalties: Option<ArtifactPenalties>,
    ) -> VcfRecord {
        let mapq = pools(&columns, alleles.len());
        let copies = vec![1.0; alleles.len()];
        VcfRecord::new(
            GenomeRegion {
                contig: ContigId(contig),
                start: Position(start),
                end: Position(end),
            },
            alleles,
            copies,
            columns,
            mapq,
            padding,
            quality(50.0),
            penalties,
            filter,
            tract,
        )
    }

    /// Every shape ng emits, named — every `FilterVerdict` spelling, both `FORMAT` strings, both
    /// padding sides, `ALT .`, a no-called sample, and a cohort at the top of the declared range.
    fn every_shape() -> Vec<(&'static str, VcfRecord)> {
        let mut a_big_cohort: Vec<SampleColumn> =
            (0..2_999).map(|_| called(&[0, 1], vec![5, 5])).collect();
        a_big_cohort.push(no_call(vec![0, 0]));

        vec![
            (
                "a biallelic SNP",
                record(
                    0,
                    100,
                    100,
                    vec![allele(b"A"), allele(b"T")],
                    vec![called(&[0, 1], vec![5, 5])],
                    None,
                    FilterVerdict::Pass,
                    None,
                    None,
                ),
            ),
            (
                "a multiallelic SNP over three samples, one of them a no-call",
                record(
                    0,
                    200,
                    200,
                    vec![allele(b"A"), allele(b"T"), allele(b"G")],
                    vec![
                        called(&[0, 1], vec![5, 5, 0]),
                        called(&[1, 2], vec![0, 4, 4]),
                        no_call(vec![0, 0, 0]),
                    ],
                    None,
                    FilterVerdict::Pass,
                    None,
                    None,
                ),
            ),
            (
                "an insertion, no padding base",
                record(
                    0,
                    300,
                    300,
                    vec![allele(b"A"), allele(b"ACCG")],
                    vec![called(&[0, 1], vec![7, 3])],
                    None,
                    FilterVerdict::Pass,
                    None,
                    None,
                ),
            ),
            (
                "a deletion with a left-hand padding base",
                record(
                    0,
                    400,
                    401,
                    vec![allele(b"AT"), allele(b"")],
                    vec![called(&[1, 1], vec![0, 9])],
                    Some(PaddingBase::Left(b'C')),
                    FilterVerdict::Pass,
                    None,
                    None,
                ),
            ),
            (
                "a deletion at the contig's first base, right-hand padding",
                record(
                    0,
                    1,
                    2,
                    vec![allele(b"AT"), allele(b"")],
                    vec![called(&[1, 1], vec![0, 9])],
                    Some(PaddingBase::Right(b'G')),
                    FilterVerdict::Pass,
                    None,
                    None,
                ),
            ),
            (
                "a repeat tract",
                record(
                    0,
                    500,
                    503,
                    vec![allele(b"ATAT"), allele(b"AT")],
                    vec![called(&[1, 1], vec![3, 3])],
                    None,
                    FilterVerdict::NotPeriodic,
                    Some(TractAnnotation::new(Motif::new(b"AT").expect("a motif"))),
                    None,
                ),
            ),
            (
                "a repeat tract whose sample is a no-call, so REPCN is `.`",
                record(
                    0,
                    600,
                    603,
                    vec![allele(b"ATAT"), allele(b"AT")],
                    vec![no_call(vec![0, 0])],
                    None,
                    FilterVerdict::TooManyAlleles,
                    Some(TractAnnotation::new(Motif::new(b"AT").expect("a motif"))),
                    None,
                ),
            ),
            (
                "a record with no alternatives, every sample a no-call",
                record(
                    0,
                    700,
                    700,
                    vec![allele(b"A")],
                    vec![no_call(vec![0]), no_call(vec![0])],
                    None,
                    FilterVerdict::LowDepth,
                    None,
                    None,
                ),
            ),
            (
                "a record carrying artifact penalties and a non-PASS filter",
                record(
                    0,
                    800,
                    800,
                    vec![allele(b"A"), allele(b"T")],
                    vec![called(&[0, 1], vec![5, 5])],
                    None,
                    FilterVerdict::EmDidNotConverge,
                    None,
                    Some(ArtifactPenalties {
                        allele_balance: quality(3.5),
                        strand_and_read_position: quality(1.25),
                    }),
                ),
            ),
            (
                "a record whose reference allele nobody's reads reached",
                record(
                    0,
                    900,
                    900,
                    vec![allele(b"A"), allele(b"T")],
                    vec![called(&[1, 1], vec![0, 8])],
                    None,
                    FilterVerdict::Pass,
                    None,
                    None,
                ),
            ),
            (
                "a record on the second contig",
                record(
                    1,
                    1000,
                    1000,
                    vec![allele(b"A"), allele(b"T")],
                    vec![called(&[0, 1], vec![5, 5])],
                    None,
                    FilterVerdict::Pass,
                    None,
                    None,
                ),
            ),
            (
                "a three-thousand-sample cohort, the top of the declared range",
                record(
                    0,
                    1100,
                    1100,
                    vec![allele(b"A"), allele(b"T")],
                    a_big_cohort,
                    None,
                    FilterVerdict::Pass,
                    None,
                    None,
                ),
            ),
        ]
    }

    /// The line the encoder writes for a record.
    fn encoded(record: &VcfRecord) -> Vec<u8> {
        record_line(record, &contigs(), diploid()).into_bytes()
    }

    #[test]
    fn every_encoded_shape_with_nothing_to_add_comes_back_byte_for_byte() {
        // Spec §10's oracle, one record at a time. A splice that rebuilt any column from its
        // parts, or that mis-counted the tabs on any shape, fails here.
        for (name, record) in every_shape() {
            let line = encoded(&record);

            let patched = rewrite_filter_and_info(&line, None, &[])
                .unwrap_or_else(|error| panic!("{name} was refused: {error}"));

            assert_eq!(
                String::from_utf8_lossy(&patched),
                String::from_utf8_lossy(&line),
                "{name} did not come back byte for byte"
            );
        }
    }

    #[test]
    fn only_the_filter_and_info_columns_move_on_every_encoded_shape() {
        for (name, record) in every_shape() {
            let line = encoded(&record);
            let before = columns(&line);

            let patched =
                rewrite_filter_and_info(&line, Some("hiddenParalog"), &the_two_info_fields())
                    .unwrap_or_else(|error| panic!("{name} was refused: {error}"));
            let after = columns(&patched);

            assert_eq!(
                after.len(),
                before.len(),
                "{name}: a column was added or lost"
            );
            for (index, (was, now)) in before.iter().zip(&after).enumerate() {
                if index == 6 || index == 7 {
                    continue;
                }
                assert_eq!(was, now, "{name}: column {index} moved");
            }
        }
    }

    #[test]
    fn the_filter_and_info_columns_of_an_encoded_line_gain_exactly_the_verdict() {
        // The two tests above skip columns seven and eight, so without this one nothing says
        // what the function is *for*: an implementation that put the ratio in `FILTER` and the
        // tag in `INFO`, or that dropped the record's own annotations, passes both of them.
        //
        // This is also where spec §6 trap 7 is pinned on lines the encoder made: the four
        // non-`PASS` verdicts are joined to, not replaced.
        let added = vec!["PARALOG_LR=4.2117".to_string()];

        for (name, record) in every_shape() {
            let line = encoded(&record);
            let before = columns(&line);

            let patched = rewrite_filter_and_info(&line, Some("hiddenParalog"), &added)
                .unwrap_or_else(|error| panic!("{name} was refused: {error}"));
            let after = columns(&patched);

            let expected_filter = if before[6] == "PASS" {
                "hiddenParalog".to_string()
            } else {
                format!("{};hiddenParalog", before[6])
            };
            assert_eq!(after[6], expected_filter, "{name}: the FILTER column");
            assert_eq!(
                after[7],
                format!("{};PARALOG_LR=4.2117", before[7]),
                "{name}: the INFO column"
            );
        }
    }

    #[test]
    fn the_encoder_writes_neither_a_missing_nor_an_empty_filter_or_info() {
        // Why the `.`-and-empty branches of the patch are defensive rather than reachable, said
        // once against the encoder instead of assumed: `info_column` pushes `AN=` and `DP=`
        // before any conditional field, and every `FilterVerdict` is a word.
        for (name, record) in every_shape() {
            let columns = columns(&encoded(&record));

            assert!(!columns[6].is_empty(), "{name}: FILTER is empty");
            assert_ne!(columns[6], ".", "{name}: FILTER is `.`");
            assert!(!columns[7].is_empty(), "{name}: INFO is empty");
            assert_ne!(columns[7], ".", "{name}: INFO is `.`");
        }
    }
}

/// **The two laws, over a domain the fixtures cannot span.**
///
/// The fixtures above are about a dozen lines a person chose. These say the same two things over
/// arbitrary column contents and column counts: with nothing to add the function is the identity,
/// and with something to add it moves the seventh and eighth columns and no others. A splice that
/// depended on what a column *contains*, or that mishandled an empty or trailing piece, fails
/// here and passes every fixture.
mod laws {
    use super::{columns, the_two_info_fields};
    use crate::run::paralog_filter::rewrite_filter_and_info;

    /// A line's worth of columns: never a tab or a newline inside one, since those are what
    /// separates columns and lines, and the patch's doc says it does not check for them.
    const A_COLUMN: &str = "[^\t\n]{0,12}";

    proptest::proptest! {
        #[test]
        fn nothing_to_add_is_the_identity(
            fields in proptest::collection::vec(A_COLUMN, 10..40)
        ) {
            let line = fields.join("\t").into_bytes();

            let patched = rewrite_filter_and_info(&line, None, &[]).expect("ten columns or more");

            proptest::prop_assert_eq!(patched, line);
        }

        #[test]
        fn only_the_filter_and_info_columns_move(
            fields in proptest::collection::vec(A_COLUMN, 10..40)
        ) {
            let line = fields.join("\t").into_bytes();

            let patched =
                rewrite_filter_and_info(&line, Some("hiddenParalog"), &the_two_info_fields())
                    .expect("ten columns or more");

            let before = columns(&line);
            let after = columns(&patched);
            proptest::prop_assert_eq!(after.len(), before.len());
            for (index, (was, now)) in before.iter().zip(&after).enumerate() {
                if index == 6 || index == 7 {
                    continue;
                }
                proptest::prop_assert_eq!(was, now, "column {} moved", index);
            }
        }
    }
}
