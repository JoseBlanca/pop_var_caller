//! Write the interleaved fixture VCF that an outside parser judges.
//!
//! **This is Milestone C's gate** (`doc/devel/ng/impl_plan/vcf_output.md` C3): a single file
//! holding every shape ng's format can write — SNPs and repeat tracts interleaved, a
//! multi-allelic site, both padding cases, the one legal position tie, a locus the caller
//! refused, an unconverged locus, and no-called samples beside called ones — so that `bcftools`
//! can be asked whether any of it is real VCF.
//!
//! It exists because the conventions ng inherits from production's repeat-tract writer **were
//! never pushed through an external parser**, and because the interleaving itself is new. A
//! defect found here costs a day; the same defect found after the differential runs costs the
//! differential.
//!
//! Usage: `ng_vcf_fixture <output directory>` — writes `fixture.vcf` and `fixture.vcf.gz`.

use std::path::PathBuf;

use pop_var_caller::ng::types::{
    AlleleId, ContigId, GenomeRegion, Genotype, Motif, Phred, Ploidy, Position,
};
use pop_var_caller::ng::vcf::{
    FilterVerdict, HeaderContig, MapqPool, PaddingBase, SampleCall, SampleColumn, SampleReadCounts,
    TractAnnotation, VcfHeaderMetadata, VcfRecord, VcfWriter,
};

fn main() {
    let mut arguments = std::env::args().skip(1);
    let directory = PathBuf::from(
        arguments
            .next()
            .unwrap_or_else(|| panic!("usage: ng_vcf_fixture <output directory>")),
    );
    std::fs::create_dir_all(&directory).expect("the output directory");

    for name in ["fixture.vcf", "fixture.vcf.gz"] {
        let path = directory.join(name);
        let mut writer = VcfWriter::create(&path, metadata(), diploid())
            .unwrap_or_else(|error| panic!("opening {}: {error}", path.display()));
        for record in every_shape() {
            writer
                .write_record(&record)
                .unwrap_or_else(|error| panic!("writing to {}: {error}", path.display()));
        }
        let written = writer.records_written();
        writer
            .finish()
            .unwrap_or_else(|error| panic!("finishing {}: {error}", path.display()));
        println!("{} — {written} records", path.display());
    }
}

fn metadata() -> VcfHeaderMetadata {
    VcfHeaderMetadata::try_new(
        vec![HeaderContig {
            name: "chr1".to_string(),
            length: 1_000_000,
            md5: Some([
                0x6a, 0xef, 0x89, 0x7c, 0x3d, 0x6f, 0xf0, 0xc7, 0x8a, 0xff, 0x06, 0xac, 0x18, 0x91,
                0x78, 0xdd,
            ]),
        }],
        vec!["HG002".to_string(), "HG003".to_string()],
        "ng_vcf_fixture".to_string(),
        "/genomes/ref.fa".to_string(),
        "run.parameters.toml".to_string(),
    )
    .expect("well-formed header metadata")
}

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("two copies is a genome")
}

fn phred(value: f32) -> Phred {
    Phred::try_new(value).expect("a non-negative finite quality")
}

fn allele(bases: &[u8]) -> Box<[u8]> {
    bases.to_vec().into_boxed_slice()
}

fn region(start: u64, end: u64) -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(0),
        start: Position(start),
        end: Position(end),
    }
}

fn called(alleles: &[u16], quality: f32, reads: Vec<u32>, unexplained: u32) -> SampleColumn {
    SampleColumn {
        call: SampleCall::Called {
            genotype: Genotype::new(alleles.iter().copied().map(AlleleId).collect()),
            genotype_quality: phred(quality),
        },
        read_counts: SampleReadCounts::new(reads, unexplained),
    }
}

fn no_call(reads: Vec<u32>, unexplained: u32) -> SampleColumn {
    SampleColumn {
        call: SampleCall::NoCall,
        read_counts: SampleReadCounts::new(reads, unexplained),
    }
}

/// Mapping-quality pools that match what the columns' `AD` attributes, which the record type
/// requires — every allele gets a slightly different mean so a swap would be visible.
fn pools(columns: &[SampleColumn], allele_count: usize) -> Vec<MapqPool> {
    (0..allele_count)
        .map(|allele| {
            let reads: u64 = columns
                .iter()
                .map(|column| u64::from(column.read_counts.allele_reads()[allele]))
                .sum();
            let mean = 60 - u64::try_from(allele).expect("a small allele index");
            MapqPool {
                reads,
                mapq_sum: reads * mean,
            }
        })
        .collect()
}

#[allow(
    clippy::too_many_arguments,
    reason = "it mirrors `VcfRecord::new`'s own list, minus the mapping-quality pools it \
              derives; giving the fixture a shape of its own would make the records harder to \
              read against the type they build"
)]
fn record(
    region: GenomeRegion,
    alleles: Vec<Box<[u8]>>,
    expected_copies: Vec<f64>,
    columns: Vec<SampleColumn>,
    padding: Option<PaddingBase>,
    quality: f32,
    filter: FilterVerdict,
    tract: Option<TractAnnotation>,
) -> VcfRecord {
    let mapq = pools(&columns, alleles.len());
    VcfRecord::new(
        region,
        alleles,
        expected_copies,
        columns,
        mapq,
        padding,
        phred(quality),
        None,
        filter,
        tract,
    )
}

fn at_motif(bases: &[u8]) -> TractAnnotation {
    TractAnnotation::new(Motif::new(bases).expect("a motif of one to six bases"))
}

/// Every shape the format can write, in genome order.
fn every_shape() -> Vec<VcfRecord> {
    vec![
        // A full-tract deletion at the contig's very first base: padded from the right, POS
        // unmoved. The case production's tract writer fills with an invented `N`.
        record(
            region(1, 4),
            vec![allele(b"ATAT"), allele(b"")],
            vec![1.0, 3.0],
            vec![
                called(&[1, 1], 40.0, vec![0, 9], 0),
                called(&[0, 1], 25.0, vec![4, 5], 1),
            ],
            Some(PaddingBase::Right(b'G')),
            88.0,
            FilterVerdict::Pass,
            Some(at_motif(b"AT")),
        ),
        // A biallelic SNP, with the second sample no-called on reads it did have.
        record(
            region(100, 100),
            vec![allele(b"A"), allele(b"T")],
            vec![2.6, 1.4],
            vec![
                called(&[0, 1], 45.0, vec![12, 9], 0),
                no_call(vec![1, 1], 0),
            ],
            None,
            310.0,
            FilterVerdict::Pass,
            None,
        ),
        // A repeat tract padded onto that SNP's position — the one legal tie.
        record(
            region(101, 106),
            vec![allele(b"ATATAT"), allele(b"")],
            vec![1.2, 2.8],
            vec![
                called(&[1, 1], 38.0, vec![0, 11], 2),
                called(&[0, 1], 22.0, vec![3, 4], 3),
            ],
            Some(PaddingBase::Left(b'C')),
            120.5,
            FilterVerdict::Pass,
            Some(at_motif(b"AT")),
        ),
        // A multi-allelic site: two samples carrying different alternatives.
        record(
            region(200, 200),
            vec![allele(b"A"), allele(b"C"), allele(b"G")],
            vec![2.0, 1.0, 1.0],
            vec![
                called(&[0, 1], 35.0, vec![6, 6, 0], 0),
                called(&[0, 2], 33.0, vec![7, 0, 5], 1),
            ],
            None,
            250.0,
            FilterVerdict::Pass,
            None,
        ),
        // An ordinary repeat tract: no padding, two lengths.
        record(
            region(300, 315),
            vec![allele(b"CACACACACACACACA"), allele(b"CACACACACACA")],
            vec![2.2, 1.8],
            vec![
                called(&[0, 1], 41.0, vec![10, 8], 4),
                called(&[1, 1], 29.0, vec![1, 12], 3),
            ],
            None,
            175.0,
            FilterVerdict::Pass,
            Some(at_motif(b"CA")),
        ),
        // A tract the caller looked at and could not call: reference alone, every sample a
        // no-call, quality zero, on its filter.
        record(
            region(500, 511),
            vec![allele(b"ATATATATATAT")],
            vec![0.0],
            vec![no_call(vec![0], 2), no_call(vec![0], 0)],
            None,
            0.0,
            FilterVerdict::LowDepth,
            Some(at_motif(b"AT")),
        ),
        // A deletion away from the contig start: padded from the left, written at 699.
        record(
            region(700, 705),
            vec![allele(b"GTGTGT"), allele(b"")],
            vec![0.6, 3.4],
            vec![
                called(&[1, 1], 44.0, vec![0, 15], 0),
                called(&[1, 1], 36.0, vec![0, 13], 2),
            ],
            Some(PaddingBase::Left(b'T')),
            99.9,
            FilterVerdict::Pass,
            Some(at_motif(b"GT")),
        ),
        // A locus whose calling loop ran out of passes: written, on its filter.
        record(
            region(800, 800),
            vec![allele(b"C"), allele(b"G")],
            vec![2.1, 1.9],
            vec![
                called(&[0, 1], 12.0, vec![3, 3], 0),
                called(&[0, 1], 11.0, vec![2, 2], 1),
            ],
            None,
            31.25,
            FilterVerdict::EmDidNotConverge,
            None,
        ),
    ]
}
