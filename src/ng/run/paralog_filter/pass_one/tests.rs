//! What pass one has to get right about a record before anything can score it.
//!
//! Two of these matter more than the rest. **The entry's position is the *written* one**, so
//! pass three's ordering check runs against a position the file actually contains — an entry
//! carrying the span start would order a left-padded deletion against a base one further on than
//! the `POS` column says. And **the row shape follows the record's tract flag**, because that is
//! what decides whether the scorer ever sees an allele split (spec §3.2).

use super::{CalledRecordSink, entry_for};
use crate::ng::calling::quality::artifact_correction::ArtifactPenalties;
use crate::ng::run::paralog_filter::SpilledSamples;
use crate::ng::types::{
    AlleleId, ContigId, GenomeRegion, Genotype, Motif, Phred, Ploidy, Position,
};
use crate::ng::vcf::encode::record_line;
use crate::ng::vcf::{
    FilterVerdict, HeaderContig, MapqPool, PaddingBase, SampleCall, SampleColumn, SampleReadCounts,
    TractAnnotation, VcfRecord,
};
use crate::ng::window_coverage::WindowCoverage;

/// Two contigs, so a fixture can sit somewhere other than the first.
fn two_contigs() -> Vec<HeaderContig> {
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

fn contigs() -> Vec<HeaderContig> {
    vec![HeaderContig {
        name: "chr1".to_string(),
        length: 1_000_000,
        md5: None,
    }]
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

fn called(alleles: &[u16], counts: Vec<u32>) -> SampleColumn {
    called_with_unexplained(alleles, counts, 0)
}

/// A called sample that also saw reads **no written allele explains** — the `DP` column counts
/// them, the per-allele counts do not. Every other fixture here sets that to zero, which is what
/// makes `Σ AD[1..]` and `DP − AD[0]` indistinguishable; on real reads it is routinely non-zero.
fn called_with_unexplained(alleles: &[u16], counts: Vec<u32>, unexplained: u32) -> SampleColumn {
    SampleColumn {
        call: SampleCall::Called {
            genotype: Genotype::new(alleles.iter().map(|id| AlleleId(*id)).collect()),
            genotype_quality: quality(30.0),
        },
        read_counts: SampleReadCounts::new(counts, unexplained),
    }
}

fn no_call(counts: Vec<u32>) -> SampleColumn {
    SampleColumn {
        call: SampleCall::NoCall,
        read_counts: SampleReadCounts::new(counts, 0),
    }
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

#[allow(clippy::too_many_arguments)]
fn record(
    start: u64,
    end: u64,
    alleles: Vec<Box<[u8]>>,
    columns: Vec<SampleColumn>,
    padding: Option<PaddingBase>,
    tract: Option<TractAnnotation>,
) -> VcfRecord {
    let mapq = pools(&columns, alleles.len());
    let copies = vec![1.0; alleles.len()];
    VcfRecord::new(
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        },
        alleles,
        copies,
        columns,
        mapq,
        padding,
        quality(50.0),
        None::<ArtifactPenalties>,
        FilterVerdict::Pass,
        tract,
    )
}

/// A scratch directory inside the project, per `CLAUDE.md` — never the system temp.
fn scratch(name: &str) -> std::path::PathBuf {
    let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tmp")
        .join("paralog_pass_one_tests")
        .join(name);
    std::fs::create_dir_all(&directory).expect("a scratch directory");
    directory
}

fn a_window(mean_depth: f32) -> WindowCoverage {
    WindowCoverage {
        gc_fraction: 0.4,
        mean_depth,
    }
}

#[test]
fn a_left_padded_deletion_is_parked_at_the_position_it_is_written_at() {
    // **The trap this test exists for.** A deletion that pads left is written one base before its
    // span, and pass three re-runs the writer's ordering check from the entry's head fields — so
    // an entry carrying the span start would order the file against a base the file does not
    // contain, and the two would disagree by exactly one on every padded deletion.
    let deletion = record(
        400,
        401,
        vec![allele(b"AT"), allele(b"")],
        vec![called(&[1, 1], vec![0, 9])],
        Some(PaddingBase::Left(b'C')),
        None,
    );

    let entry = entry_for(&deletion, &[a_window(5.0)], &contigs(), diploid());

    assert_eq!(
        entry.position,
        Position(399),
        "the span starts at 400 and the padding moves the written position to 399"
    );
    // And the entry agrees with the line beside it, which is the thing that matters.
    let line = record_line(&deletion, &contigs(), diploid());
    let written_pos = line.split('\t').nth(1).expect("a POS column");
    assert_eq!(written_pos, "399");
    assert_eq!(entry.position.get().to_string(), written_pos);
}

#[test]
fn a_records_line_is_parked_exactly_as_the_encoder_wrote_it() {
    let snp = record(
        100,
        100,
        vec![allele(b"A"), allele(b"T")],
        vec![called(&[0, 1], vec![5, 5])],
        None,
        None,
    );

    let entry = entry_for(&snp, &[a_window(5.0)], &contigs(), diploid());

    assert_eq!(
        entry.line,
        record_line(&snp, &contigs(), diploid()).into_bytes(),
        "the parked line must be the encoder's, since pass three writes it back unchanged"
    );
}

#[test]
fn a_repeat_tract_is_parked_with_windows_and_no_read_counts() {
    // Spec §3.2's rule, at the point it is applied: the tract flag chooses the row shape, and a
    // tract's rows carry no allele split because slippage has made it unreadable.
    let tract = record(
        500,
        503,
        vec![allele(b"ATAT"), allele(b"AT")],
        vec![called(&[1, 1], vec![3, 3]), called(&[0, 0], vec![6, 0])],
        None,
        Some(TractAnnotation::new(Motif::new(b"AT").expect("a motif"))),
    );

    let entry = entry_for(
        &tract,
        &[a_window(5.0), a_window(6.0)],
        &contigs(),
        diploid(),
    );

    assert!(entry.is_repeat_tract);
    match entry.samples {
        SpilledSamples::RepeatTract(rows) => {
            assert_eq!(rows.len(), 2);
            assert_eq!(rows[0].window.mean_depth, 5.0);
            assert_eq!(rows[1].window.mean_depth, 6.0);
        }
        SpilledSamples::GenericLocus(_) => {
            panic!("a repeat tract was parked with an allele split it cannot report")
        }
    }
}

#[test]
fn a_multiallelic_site_is_parked_with_its_alternatives_summed() {
    // Pooling is what lets a multiallelic site carry the allele signal at all: the story asks
    // what share of the reads is non-reference, not which alternative each one carried.
    let multiallelic = record(
        200,
        200,
        vec![allele(b"A"), allele(b"T"), allele(b"G")],
        vec![called(&[1, 2], vec![2, 3, 4])],
        None,
        None,
    );

    let entry = entry_for(&multiallelic, &[a_window(5.0)], &contigs(), diploid());

    match entry.samples {
        SpilledSamples::GenericLocus(rows) => {
            assert_eq!(rows[0].ref_reads, 2);
            assert_eq!(
                rows[0].alt_reads, 7,
                "3 + 4, not just the first alternative"
            );
        }
        SpilledSamples::RepeatTract(_) => panic!("a SNP was parked as a tract"),
    }
}

#[test]
fn a_record_with_no_alternative_parks_a_zero_alternative_count() {
    // A refused locus written `ALT .` has one allele and no split to report. It is still parked,
    // because pass three has to write it back.
    let refused = record(
        700,
        700,
        vec![allele(b"A")],
        vec![called(&[0, 0], vec![4])],
        None,
        None,
    );

    let entry = entry_for(&refused, &[a_window(5.0)], &contigs(), diploid());

    match entry.samples {
        SpilledSamples::GenericLocus(rows) => {
            assert_eq!(rows[0].ref_reads, 4);
            assert_eq!(rows[0].alt_reads, 0);
        }
        SpilledSamples::RepeatTract(_) => panic!("a generic locus was parked as a tract"),
    }
}

#[test]
fn the_windows_reach_the_entry_in_the_runs_sample_order() {
    // The slice is dense over the run's samples and so are the record's columns; parking them
    // out of step would give every sample its neighbour's coverage, which no later step can
    // detect.
    let snp = record(
        100,
        100,
        vec![allele(b"A"), allele(b"T")],
        vec![
            called(&[0, 1], vec![5, 5]),
            called(&[0, 0], vec![9, 0]),
            called(&[1, 1], vec![0, 7]),
        ],
        None,
        None,
    );

    let entry = entry_for(
        &snp,
        &[a_window(1.0), a_window(2.0), a_window(3.0)],
        &contigs(),
        diploid(),
    );

    match entry.samples {
        SpilledSamples::GenericLocus(rows) => {
            assert_eq!(
                rows.iter().map(|r| r.window.mean_depth).collect::<Vec<_>>(),
                vec![1.0, 2.0, 3.0]
            );
            assert_eq!(
                rows.iter().map(|r| r.ref_reads).collect::<Vec<_>>(),
                vec![5, 9, 0],
                "a sample was given another's read counts"
            );
        }
        SpilledSamples::RepeatTract(_) => panic!("a SNP was parked as a tract"),
    }
}

#[test]
fn the_entry_a_record_parks_can_be_written_back_where_the_record_would_have_gone() {
    // The whole point of the head fields: pass three rebuilds the writer's ordering key from the
    // entry, without the record. This is that key, checked against the one the writer builds
    // from the record itself.
    let deletion = record(
        400,
        401,
        vec![allele(b"AT"), allele(b"")],
        vec![called(&[1, 1], vec![0, 9])],
        Some(PaddingBase::Left(b'C')),
        None,
    );
    let tract = record(
        500,
        503,
        vec![allele(b"ATAT"), allele(b"AT")],
        vec![called(&[1, 1], vec![3, 3])],
        None,
        Some(TractAnnotation::new(Motif::new(b"AT").expect("a motif"))),
    );

    for source in [&deletion, &tract] {
        let entry = entry_for(source, &[a_window(5.0)], &contigs(), diploid());
        let from_the_entry = crate::ng::vcf::RecordPlace::from(&entry);

        assert_eq!(from_the_entry.at.contig, source.region().contig);
        assert_eq!(
            from_the_entry.is_repeat_tract,
            source.is_repeat_tract(),
            "the tract flag decides the one legal tie, so it has to survive the spill"
        );
        assert_eq!(
            from_the_entry.at.position.get().to_string(),
            record_line(source, &contigs(), diploid())
                .split('\t')
                .nth(1)
                .expect("a POS column"),
            "the place pass three checks must be the position the line says"
        );
    }
}

#[test]
fn a_sink_that_parks_records_does_not_create_the_vcf() {
    // Spec §2: no record's verdict is known until every record has been scored, so a VCF written
    // during pass one would be a file someone could mistake for the finished output.
    let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tmp")
        .join("paralog_pass_one_tests")
        .join("parking_does_not_create_the_vcf");
    std::fs::create_dir_all(&directory).expect("a scratch directory");
    let output = directory.join("cohort.vcf");
    let _ = std::fs::remove_file(&output);

    let mut sink = CalledRecordSink::ParkedOnTheSpill(super::SpillingSink::beside(
        &output,
        contigs(),
        diploid(),
    ));
    let snp = record(
        100,
        100,
        vec![allele(b"A"), allele(b"T")],
        vec![called(&[0, 1], vec![5, 5])],
        None,
        None,
    );

    sink.accept(&snp, &[a_window(5.0)]).expect("the record");

    assert!(
        !output.exists(),
        "pass one created the VCF, which cannot be written until every record has been scored"
    );
    if let CalledRecordSink::ParkedOnTheSpill(sink) = &sink {
        assert!(sink.spill().path().exists(), "the spill was not created");
    }
}

#[test]
fn the_alternative_count_is_the_written_alleles_and_not_the_depth_minus_the_reference() {
    // **The two rules agree on every other fixture in this file and disagree on real reads.**
    // A sample's `DP` counts reads no written allele explains — reads that reached the locus and
    // supported nothing the record lists — and `Σ AD[1..]` does not. Deriving the alternative
    // count as `DP − AD[0]` would fold those into the non-reference share the scorer sees, at
    // every locus in every sample, with nothing to notice it.
    //
    // Here: 4 reference, 5 alternative, and 6 reads explained by neither. The answer is 5.
    let snp = record(
        100,
        100,
        vec![allele(b"A"), allele(b"T")],
        vec![called_with_unexplained(&[0, 1], vec![4, 5], 6)],
        None,
        None,
    );

    let entry = entry_for(&snp, &[a_window(5.0)], &contigs(), diploid());

    match entry.samples {
        SpilledSamples::GenericLocus(rows) => {
            assert_eq!(rows[0].ref_reads, 4);
            assert_eq!(
                rows[0].alt_reads, 5,
                "the six reads no allele explains are not evidence for the alternative"
            );
        }
        SpilledSamples::RepeatTract(_) => panic!("a SNP was parked as a tract"),
    }
}

#[test]
fn a_record_is_parked_on_the_contig_it_was_called_on() {
    // Every other fixture here sits on contig 0, so hardcoding the contig would pass all of them
    // — and a spill whose entries all claim one contig orders the whole file against the wrong
    // sequence at pass three.
    let on_the_second = VcfRecord::new(
        GenomeRegion {
            contig: ContigId(1),
            start: Position(100),
            end: Position(100),
        },
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        vec![called(&[0, 1], vec![5, 5])],
        pools(&[called(&[0, 1], vec![5, 5])], 2),
        None,
        quality(50.0),
        None::<ArtifactPenalties>,
        FilterVerdict::Pass,
        None,
    );

    let entry = entry_for(&on_the_second, &[a_window(5.0)], &two_contigs(), diploid());

    assert_eq!(entry.contig, ContigId(1));
    assert!(
        entry.line.starts_with(b"chr2\t"),
        "the parked line names a different contig than the entry does"
    );
}

#[test]
fn the_runs_ploidy_reaches_the_parked_line() {
    // The ploidy shows in how a *no-call* is spelled — `./.` at two copies, `./././.` at four —
    // so a fixture with no no-call cannot tell a forwarded ploidy from a hardcoded one, and every
    // other fixture here has none.
    let with_a_no_call = record(
        100,
        100,
        vec![allele(b"A"), allele(b"T")],
        vec![no_call(vec![0, 0])],
        None,
        None,
    );

    let diploid_line = entry_for(&with_a_no_call, &[a_window(5.0)], &contigs(), diploid()).line;
    let tetraploid = Ploidy::try_new(4).expect("four copies");
    let tetraploid_line = entry_for(&with_a_no_call, &[a_window(5.0)], &contigs(), tetraploid).line;

    assert!(String::from_utf8_lossy(&diploid_line).contains("./."));
    assert!(String::from_utf8_lossy(&tetraploid_line).contains("./././."));
    assert_ne!(
        diploid_line, tetraploid_line,
        "the ploidy the run was given did not reach the parked line"
    );
}

#[test]
#[should_panic(expected = "disagree in length")]
fn a_record_and_a_window_slice_that_disagree_are_refused_rather_than_zipped_short() {
    // The zip would take the shorter of the two and park a prefix of the cohort. It cannot
    // mis-pair — both are the same sample order — but a lost suffix is a record scored on fewer
    // samples than the run has, which pass two would report as *its* problem.
    let snp = record(
        100,
        100,
        vec![allele(b"A"), allele(b"T")],
        vec![called(&[0, 1], vec![5, 5]), called(&[0, 0], vec![9, 0])],
        None,
        None,
    );

    let _ = entry_for(&snp, &[a_window(5.0)], &contigs(), diploid());
}

#[test]
fn parked_records_are_readable_once_the_sink_is_finished() {
    // `finish_parking` is what makes the bytes durable, and nothing exercised it: a sink that
    // never flushed would leave passes two and three reading a file shorter than the run wrote,
    // which the spill's own completeness check reports as a truncation rather than as a
    // forgotten flush.
    let directory = scratch("parked_records_are_readable");
    let output = directory.join("cohort.vcf");

    let mut sink = CalledRecordSink::ParkedOnTheSpill(super::SpillingSink::beside(
        &output,
        contigs(),
        diploid(),
    ));
    let snp = record(
        100,
        100,
        vec![allele(b"A"), allele(b"T")],
        vec![called(&[0, 1], vec![5, 5])],
        None,
        None,
    );
    let tract = record(
        500,
        503,
        vec![allele(b"ATAT"), allele(b"AT")],
        vec![called(&[1, 1], vec![3, 3])],
        None,
        Some(TractAnnotation::new(Motif::new(b"AT").expect("a motif"))),
    );

    sink.accept(&snp, &[a_window(5.0)]).expect("the SNP");
    sink.accept(&tract, &[a_window(6.0)]).expect("the tract");

    let CalledRecordSink::ParkedOnTheSpill(mut parked) = sink else {
        panic!("the sink changed variant")
    };
    parked.finish_parking().expect("the parked records");

    let read: Vec<_> = parked
        .spill()
        .read()
        .expect("the parked file")
        .map(|entry| entry.expect("every entry the sink wrote reads back"))
        .collect();

    assert_eq!(read.len(), 2, "both records must survive the flush");
    assert!(!read[0].is_repeat_tract);
    assert!(read[1].is_repeat_tract);
    assert_eq!(
        read[0].line,
        record_line(&snp, &contigs(), diploid()).into_bytes()
    );
}

#[test]
fn a_sink_that_cannot_write_says_which_half_of_pass_one_failed() {
    // The two error variants are the two paths, and telling them apart is the whole reason there
    // are two: a VCF that cannot be written and a spill that cannot be appended to send an
    // operator to different places.
    let directory = scratch("a_sink_that_cannot_write");
    let output = directory.join("nowhere").join("deeper").join("cohort.vcf");

    let mut sink = CalledRecordSink::ParkedOnTheSpill(super::SpillingSink::beside(
        &output,
        contigs(),
        diploid(),
    ));
    let snp = record(
        100,
        100,
        vec![allele(b"A"), allele(b"T")],
        vec![called(&[0, 1], vec![5, 5])],
        None,
        None,
    );

    let refused = sink
        .accept(&snp, &[a_window(5.0)])
        .expect_err("there is no directory to park in");

    assert!(
        matches!(refused, super::PassOneError::Spill(_)),
        "a spill failure was reported as a VCF failure: {refused:?}"
    );
}
