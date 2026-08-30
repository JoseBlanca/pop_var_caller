//! The two things the writer owns: the order records may arrive in, and a file that appears
//! whole or not at all.

use std::fs;

use super::*;
use crate::ng::types::{AlleleId, GenomeRegion, Genotype, Motif, Phred, Position};
use crate::ng::vcf::{
    FilterVerdict, HeaderContig, MapqPool, SampleCall, SampleColumn, SampleReadCounts,
    TractAnnotation, VcfHeaderMetadata,
};

/// A scratch directory inside the project, per `CLAUDE.md` — never the system temp.
fn scratch(name: &str) -> PathBuf {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tmp")
        .join("vcf_writer_tests")
        .join(name);
    fs::create_dir_all(&directory).expect("a scratch directory");
    directory
}

fn metadata() -> VcfHeaderMetadata {
    VcfHeaderMetadata::try_new(
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
        ],
        vec!["one".to_string()],
        "ng call".to_string(),
        "/genomes/ref.fa".to_string(),
        "run.parameters.toml".to_string(),
    )
    .expect("well-formed header metadata")
}

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("two copies is a genome")
}

fn quality(phred: f32) -> Phred {
    Phred::try_new(phred).expect("a non-negative finite quality")
}

fn allele(bases: &[u8]) -> Box<[u8]> {
    bases.to_vec().into_boxed_slice()
}

/// A one-base SNP on `contig` at `position`, called heterozygous.
fn snp_at(contig: u32, position: u64) -> VcfRecord {
    VcfRecord::new(
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(position),
            end: Position(position),
        },
        vec![allele(b"A"), allele(b"T")],
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: SampleCall::Called {
                genotype: Genotype::new(vec![AlleleId(0), AlleleId(1)]),
                genotype_quality: quality(30.0),
            },
            read_counts: SampleReadCounts::new(vec![5, 5], 0),
        }],
        vec![
            MapqPool {
                reads: 5,
                mapq_sum: 300,
            },
            MapqPool {
                reads: 5,
                mapq_sum: 300,
            },
        ],
        None,
        quality(50.0),
        None,
        FilterVerdict::Pass,
        None,
    )
}

/// A repeat tract whose span starts at `span_start`; with `padded`, a full-tract deletion whose
/// written position is one base earlier.
fn tract_at(contig: u32, span_start: u64, padded: bool) -> VcfRecord {
    let alleles = if padded {
        vec![allele(b"ATAT"), allele(b"")]
    } else {
        vec![allele(b"ATAT"), allele(b"AT")]
    };
    let counts = if padded { vec![0, 6] } else { vec![3, 3] };
    VcfRecord::new(
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(span_start),
            end: Position(span_start + 3),
        },
        alleles,
        vec![1.0, 1.0],
        vec![SampleColumn {
            call: SampleCall::Called {
                genotype: Genotype::new(vec![AlleleId(1), AlleleId(1)]),
                genotype_quality: quality(30.0),
            },
            read_counts: SampleReadCounts::new(counts.clone(), 0),
        }],
        counts
            .iter()
            .map(|reads| MapqPool {
                reads: u64::from(*reads),
                mapq_sum: u64::from(*reads) * 60,
            })
            .collect(),
        padded.then_some(PaddingBase::Left(b'C')),
        quality(50.0),
        None,
        FilterVerdict::Pass,
        Some(TractAnnotation::new(
            Motif::new(b"AT").expect("a two-base motif"),
        )),
    )
}

// ---------------------------------------------------------------------------
// Order
// ---------------------------------------------------------------------------

#[test]
fn records_in_genome_order_are_accepted() {
    let path = scratch("in_order").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    for position in [100, 200, 300] {
        writer.write_record(&snp_at(0, position)).expect("in order");
    }
    // A later contig starts over at a lower position, and that is forward, not backward.
    writer.write_record(&snp_at(1, 5)).expect("a later contig");
    assert_eq!(writer.records_written(), 4);
    writer.finish().expect("the file finishes");
}

#[test]
fn a_record_that_runs_backwards_is_refused() {
    let path = scratch("backwards").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer.write_record(&snp_at(0, 200)).expect("the first");
    let refused = writer.write_record(&snp_at(0, 100));
    assert!(matches!(
        refused,
        Err(VcfWriteError::OutOfOrder {
            previous_position: 200,
            position: 100,
            ..
        })
    ));
}

#[test]
fn an_earlier_contig_after_a_later_one_is_refused() {
    let path = scratch("contig_backwards").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer.write_record(&snp_at(1, 100)).expect("the first");
    assert!(matches!(
        writer.write_record(&snp_at(0, 100)),
        Err(VcfWriteError::OutOfOrder { .. })
    ));
}

#[test]
fn a_generic_locus_and_a_tract_padded_onto_it_may_share_a_position() {
    // **The one legal tie.** The tract's span starts at 101; its deletion pads left, so it is
    // written at 100 — the position of the SNP that owns the anchor base. The two describe
    // different bases, so both belong in the file.
    let path = scratch("legal_tie").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer
        .write_record(&snp_at(0, 100))
        .expect("the generic one");
    writer
        .write_record(&tract_at(0, 101, true))
        .expect("the tract padded onto it");
    assert_eq!(writer.records_written(), 2);
    writer.finish().expect("the file finishes");
}

#[test]
fn a_tract_followed_by_a_generic_locus_at_one_position_is_refused() {
    // The order within the tie is fixed: the generic record's span genuinely starts there and
    // the tract's starts one base later.
    let path = scratch("tie_wrong_way").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer
        .write_record(&tract_at(0, 101, true))
        .expect("the tract first");
    assert!(matches!(
        writer.write_record(&snp_at(0, 100)),
        Err(VcfWriteError::IllegalTie {
            position: 100,
            previous_was_repeat_tract: true,
            is_repeat_tract: false,
            ..
        })
    ));
}

#[test]
fn two_generic_records_at_one_position_are_refused() {
    let path = scratch("two_generic").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer.write_record(&snp_at(0, 100)).expect("the first");
    assert!(matches!(
        writer.write_record(&snp_at(0, 100)),
        Err(VcfWriteError::IllegalTie { .. })
    ));
}

#[test]
fn a_third_record_at_one_position_is_refused() {
    // The tie is admitted once: after the legal pair, anything else at that position is a bug.
    let path = scratch("third_at_one_position").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer
        .write_record(&snp_at(0, 100))
        .expect("the generic one");
    writer
        .write_record(&tract_at(0, 101, true))
        .expect("the tract padded onto it");
    assert!(matches!(
        writer.write_record(&tract_at(0, 101, true)),
        Err(VcfWriteError::IllegalTie {
            previous_was_repeat_tract: true,
            ..
        })
    ));
}

#[test]
fn the_order_is_checked_on_the_written_position_not_the_span_start() {
    // A tract whose span starts at 101 is *written* at 100 when it pads left. A writer that
    // checked the span start would accept a SNP at 100 after it, producing a file whose
    // positions run backwards.
    let path = scratch("written_position").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer
        .write_record(&tract_at(0, 101, true))
        .expect("written at 100");
    // 100 would be a legal *span* successor to 101-padded-to-100 only if the check used the
    // span; against the written position it is a tie in the wrong order.
    assert!(writer.write_record(&snp_at(0, 100)).is_err());
}

// ---------------------------------------------------------------------------
// The file
// ---------------------------------------------------------------------------

#[test]
fn the_finished_file_holds_the_header_and_every_record() {
    let path = scratch("whole_file").join("out.vcf");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer.write_record(&snp_at(0, 100)).expect("a record");
    writer.write_record(&snp_at(0, 200)).expect("a record");
    writer.finish().expect("the file finishes");

    let written = fs::read_to_string(&path).expect("the finished file");
    assert!(written.starts_with("##fileformat=VCFv4.4\n"));
    let records: Vec<&str> = written
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(records.len(), 2);
    assert!(records[0].starts_with("chr1\t100\t"));
    assert!(records[1].starts_with("chr1\t200\t"));
    // Every line ends with a newline, including the last.
    assert!(written.ends_with('\n'));
}

#[test]
fn nothing_appears_at_the_output_path_until_the_writer_finishes() {
    // A crash before `finish` leaves the in-flight file and no output, so nobody mistakes a
    // half-written VCF for a finished run.
    let directory = scratch("only_on_finish");
    let path = directory.join("out.vcf");
    let _ = fs::remove_file(&path);

    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer.write_record(&snp_at(0, 100)).expect("a record");
    assert!(!path.exists(), "the output must not exist yet");
    assert!(
        directory.join("out.vcf.tmp").exists(),
        "the in-flight file should be there instead"
    );

    writer.finish().expect("the file finishes");
    assert!(path.exists(), "the output appears on finish");
    assert!(
        !directory.join("out.vcf.tmp").exists(),
        "and the in-flight file is gone"
    );
}

#[test]
fn a_gzipped_name_gets_a_bgzf_file_ending_in_the_marker_htslib_requires() {
    let path = scratch("bgzf").join("out.vcf.gz");
    let mut writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    writer.write_record(&snp_at(0, 100)).expect("a record");
    writer.finish().expect("the file finishes");

    let bytes = fs::read(&path).expect("the finished file");
    assert!(bytes.ends_with(BGZF_EOF), "bgzf files end with the marker");
    // And it is not plain text: the gzip magic is the first two bytes.
    assert_eq!(&bytes[..2], &[0x1f, 0x8b]);
}

#[test]
fn the_suffix_is_matched_whatever_its_case() {
    assert!(path_is_bgzf(Path::new("a/out.vcf.gz")));
    assert!(path_is_bgzf(Path::new("a/OUT.VCF.GZ")));
    assert!(path_is_bgzf(Path::new("a/out.vcf.bgz")));
    assert!(!path_is_bgzf(Path::new("a/out.vcf")));
    assert!(!path_is_bgzf(Path::new("a/out.txt")));
}

#[test]
fn a_run_that_called_nothing_still_writes_a_header() {
    // An empty cohort VCF is a legitimate answer — no locus passed — and a consumer should
    // meet a well-formed file saying so rather than an absent one.
    let path = scratch("no_records").join("out.vcf");
    let writer = VcfWriter::create(&path, metadata(), diploid()).expect("the writer opens");
    assert_eq!(writer.records_written(), 0);
    writer.finish().expect("the file finishes");

    let written = fs::read_to_string(&path).expect("the finished file");
    assert!(written.contains("#CHROM\t"));
    assert_eq!(
        written
            .lines()
            .filter(|line| !line.starts_with('#'))
            .count(),
        0
    );
}
