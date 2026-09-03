//! **Who genotypes an indel that sits on the line between a repeat tract and its
//! flank** — ownership checked over both calling paths at once, end to end.
//!
//! The partition gives a tract's own span to the repeat-tract path and its flanks to
//! the SNP/indel path (`doc/devel/ng/spec/typed_regions.md` §2.2). An indel *event*
//! can sit exactly on that line: a mapper anchors an insertion between the last flank
//! base and the tract's first base **on the flank base**, which is SNP/indel
//! territory, while the ruling is that the event belongs to the tract. Each half of
//! that ruling has a unit test of its own — `open_record.rs` pins that an insertion
//! anchored on a generic region's last base contributes no allele (the `b6309954`
//! rule), and `ssr.rs` pins that a left-junction insertion whose last base differs
//! from the flank base beside it is spelled into the tract's observation (the L4
//! junction convention). What neither can check, and these tests do, is the two
//! paths **together**: that the event reaches the output genotyped exactly once, by
//! the tract path; and that a deletion reaching across the line is genotyped twice —
//! once per path, each over exactly its own bases.
//!
//! **Everything here is deliberately self-contained** — its own reference, its own
//! ground, local copies of three tiny fixture helpers from `callers.rs`'s tests —
//! so the file compiles and runs without reaching into the test modules of files
//! under active change.
//!
//! **The mirror case is known-broken and deliberately not pinned here**: an
//! insertion anchored on the *tract's last* base whose spelling the junction
//! convention hands to the right flank is expelled by the tract path and never seen
//! by the generic walk — no path owns it
//! (`doc/devel/ng/research/tract_accuracy_program_report.md` §L4, "right→silent").
//! A test for it belongs with the fix.

use std::io::Write as _;
use std::path::PathBuf;
use std::sync::Arc;

use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
use tempfile::TempDir;

use crate::ng::calling::LocusInference;
use crate::ng::calling::allele_candidates::CandidateSelectionConfig;
use crate::ng::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use crate::ng::calling::inference::CallingLoopConfig;
use crate::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use crate::ng::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use crate::ng::calling::parameters_file::DeclaredInbreeding;
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::locus_generation::LocusKind;
use crate::ng::locus_generation::pileup::PileupGeneratorConfig;
use crate::ng::read::filtering::ReadFilterConfig;
use crate::ng::read::input::read_groups::build_read_groups;
use crate::ng::read::input::reference::OpenReference;
use crate::ng::read::input::test_fixtures::{
    header, indexed_named_bam, matching_contigs, read_named_with_length,
};
use crate::ng::reference_info::{ReferenceSource, read_reference_info};
use crate::ng::region_typing::segment_criteria::SsrSegment;
use crate::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion};
use crate::ng::repeat_catalog::{RepeatCatalogHeader, StrRepeatCriteria};
use crate::ng::run::callers::{AlignedFilesVariantCaller, AlignmentInputs, MergeParameters};
use crate::ng::run::segments::Segmentation;
use crate::ng::tandem_repeat::ScanParams;
use crate::ng::types::{AlleleId, ContigId, GenomeRegion, Motif, Ploidy, Position};
use crate::regions::ContigBounds;

// ---------------------------------------------------------------------------
// The ground: a real tract with distinct flanks.
//
// `callers.rs`'s tract fixture runs over a reference of a hundred identical `A`s,
// which is exactly wrong for a junction test: an indel there is ambiguous under
// left-alignment and the declared AT tract spells no repeat. This reference makes
// every junction event unambiguous by construction: the flank's last base (`C` at
// 40) differs from the tract's first (`A` at 41), from the tract's alphabet, and
// from every inserted base a test uses, so no representation of any fixture indel
// can slide across the boundary.
// ---------------------------------------------------------------------------

/// chr1:1–25, sequence before the flank a read starts in.
const FILLER: &[u8] = b"TGCATGCCTAGGATCCTGAACGTCG";
/// chr1:26–40 — the fifteen flank bases the tract generator fetches on the left.
/// Ends `...GC`, so position 38 is `C`, 39 `G`, 40 `C`.
const LEFT_FLANK: &[u8] = b"CTGGACGTCTACCGC";
/// chr1:41–52 — the tract itself: `AT` × 6.
const TRACT: &[u8] = b"ATATATATATAT";
/// chr1:53–100. Starts with `C` so the repeat ends where the segment says it does.
const RIGHT: &[u8] = b"CGGCTAGGCTCCAGTTCGACCTGGACGTTCAGGCTACCGGATCGGCAC";

/// Where the tract sits, 1-based inclusive — the segment the ground declares.
const TRACT_SPAN: (u64, u64) = (41, 52);

/// The whole of chr1, positions 1–100.
fn chr1_bases() -> Vec<u8> {
    let bases = [FILLER, LEFT_FLANK, TRACT, RIGHT].concat();
    assert_eq!(bases.len(), 100, "the four blocks tile chr1 exactly");
    bases
}

/// Reference bases over 1-based inclusive positions `start..=end` of chr1.
fn chr1_window(start: usize, end: usize) -> Vec<u8> {
    chr1_bases()[start - 1..end].to_vec()
}

/// The reference on disk, opened the way a real run holds it: geometry from the
/// `.fai`, the FASTA's path kept so the walk can reach the bases. The same shape as
/// `test_fixtures::fixture_reference_from_its_index`, with chr1 carrying this
/// module's sequence instead of a hundred `A`s. chr2 stays 200 `A`s so the FASTA
/// agrees with the shared header fixture's contig list.
fn junction_reference() -> (TempDir, OpenReference) {
    let dir = tempfile::tempdir().expect("a tempdir for the reference");
    let fasta = dir.path().join("ref.fa");
    let fai = dir.path().join("ref.fa.fai");
    let contigs: [(&str, Vec<u8>); 2] = [("chr1", chr1_bases()), ("chr2", vec![b'A'; 200])];

    let mut fa = std::fs::File::create(&fasta).expect("the FASTA opens for writing");
    let mut index = std::fs::File::create(&fai).expect("the index opens for writing");
    let mut offset: u64 = 0;
    for (name, bases) in &contigs {
        let heading = format!(">{name}\n");
        fa.write_all(heading.as_bytes()).expect("heading written");
        offset += heading.len() as u64;
        fa.write_all(bases).expect("bases written");
        fa.write_all(b"\n").expect("newline written");
        let line_bases = bases.len() as u64;
        writeln!(
            index,
            "{name}\t{line_bases}\t{offset}\t{line_bases}\t{}",
            line_bases + 1
        )
        .expect("index row written");
        offset += line_bases + 1;
    }

    let from_the_index =
        read_reference_info(ReferenceSource::Fai(fai)).expect("the index reads back");
    (
        dir,
        OpenReference::new(Arc::new(crate::ng::reference_info::ReferenceInfo {
            fasta_path: Some(fasta),
            md5: from_the_index.md5,
            contigs: from_the_index.contigs.clone(),
        })),
    )
}

/// The analysed ground: ordinary stretch, the AT tract, ordinary stretch — the same
/// three-segment shape as `callers.rs`'s `ground_with_a_tract_interleaved`, over
/// this module's reference.
fn ground_with_the_tract() -> Segmentation {
    let chr1 = |start: u64, end: u64| GenomeRegion {
        contig: ContigId(0),
        start: Position(start),
        end: Position(end),
    };
    let bounds = [ContigBounds {
        name: "chr1",
        length: 100,
    }];
    let (tract_start, tract_end) = TRACT_SPAN;
    let segments = vec![
        TypedRegion {
            region: chr1(1, tract_start - 1),
            kind: RegionKind::Generic,
        },
        TypedRegion {
            region: chr1(tract_start, tract_end),
            kind: RegionKind::SsrSegment(
                SsrSegment::new(
                    "chr1".into(),
                    tract_start,
                    tract_end,
                    Motif::new(b"AT").expect("AT is a motif"),
                    1.0,
                )
                .expect("a twelve-base AT tract inside the contig"),
            ),
        },
        TypedRegion {
            region: chr1(tract_end + 1, 100),
            kind: RegionKind::Generic,
        },
    ];
    Segmentation::build(
        segments.into_iter().map(Ok),
        GenomeRegions::whole_contigs(&bounds),
        RepeatCatalogHeader {
            contigs: Vec::new(),
            reference_md5: [7; 16],
            built_under: StrRepeatCriteria::default(),
            scan: ScanParams::default(),
            tool_version: "test".to_string(),
            longest_tract_bp: Vec::new(),
        },
        StrRepeatCriteria::default(),
        PathBuf::from("/genomes/test.catalog.parquet"),
    )
    .expect("a clean stream builds")
}

// ---------------------------------------------------------------------------
// Reads and the caller.
// ---------------------------------------------------------------------------

/// A read spanning chr1:25–69 on the reference — the tract with the fifteen
/// fetched flank bases either side and one to spare — with the cigar and sequence
/// a test states. Four of these clear the merge's two-carrying-reads floor twice
/// over.
fn junction_read(name: &str, cigar: &[(Kind, usize)], sequence: Vec<u8>) -> RecordBuf {
    let mut record = read_named_with_length(name, 0, 25, 45);
    *record.cigar_mut() = cigar
        .iter()
        .map(|(kind, len)| Op::new(*kind, *len))
        .collect();
    *record.quality_scores_mut() = QualityScores::from(vec![30u8; sequence.len()]);
    *record.sequence_mut() = Sequence::from(sequence);
    record
}

/// One sample, four copies of the same read shape, indexed on disk.
fn a_sample_of(read: impl Fn(String) -> RecordBuf, file_name: &str) -> (TempDir, PathBuf) {
    let records: Vec<RecordBuf> = (0..4).map(|n| read(format!("junction-r{n}"))).collect();
    indexed_named_bam(
        &header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg1", Some("junction"))],
        ),
        &records,
        file_name,
    )
}

/// Everything the run's loci came back through, with the shipped settings — the
/// same opening `callers.rs`'s tract fixture uses, over this module's ground.
fn called_over(bam: &PathBuf, reference: &OpenReference) -> Vec<LocusInference> {
    let paths = std::slice::from_ref(bam);
    let read_groups = build_read_groups(paths).expect("the fixture declares a read group");
    let parameters = RunParameters::of_defaults(
        &read_groups,
        Ploidy::try_new(2).expect("a diploid"),
        &DeclaredInbreeding::nothing_said(),
    );
    AlignedFilesVariantCaller::open(
        AlignmentInputs {
            read_groups: &read_groups,
            reference,
            read_filters: ReadFilterConfig::default(),
            build_index_if_missing: false,
            locus_generator_settings: PileupGeneratorConfig::default(),
            reference_with_checksums: reference.info(),
        },
        ground_with_the_tract(),
        parameters,
        CallingLoopConfig::DEFAULT
            .validate()
            .expect("the shipped calling-loop settings are runnable"),
        CandidateSelectionConfig::DEFAULT,
        MergeParameters::DEFAULT,
    )
    .expect("one readable sample over a readable reference opens")
    .call_cohort(&SummariseConditionLoop::new(
        StutterSubstitutionEmission,
        MarginalizedDirichletPrior,
    ))
    .expect("the fixture cohort calls")
    .called_loci
}

/// The loci that carry any allele beyond the reference — the ones a genotype
/// could reach the output through.
fn variant_loci(called: &[LocusInference]) -> Vec<&LocusInference> {
    called
        .iter()
        .filter(|locus| locus.alleles().len() > 1)
        .collect()
}

/// The one sample's called genotype at `locus`, as (both copies same allele,
/// that allele is not the reference) — the homozygous-alternative shape every
/// fixture here expects, since all four reads carry the event.
fn assert_called_homozygous_alternative(locus: &LocusInference, what: &str) {
    let genotype = locus.per_sample[0]
        .genotype()
        .unwrap_or_else(|| panic!("{what}: the sample was called, not set aside"));
    let alleles = genotype.alleles();
    assert_eq!(
        alleles.len(),
        2,
        "{what}: a diploid call carries two copies"
    );
    assert_eq!(
        alleles[0], alleles[1],
        "{what}: all four reads carry the event, so both copies are the same allele",
    );
    assert_ne!(
        alleles[0],
        AlleleId::REFERENCE,
        "{what}: the called allele is the event, not the reference",
    );
}

/// Lengths of the non-reference alleles at `locus`, in admission order.
fn alternative_lengths(locus: &LocusInference) -> Vec<usize> {
    locus.alleles().iter().skip(1).map(<[u8]>::len).collect()
}

// ---------------------------------------------------------------------------
// The tests.
// ---------------------------------------------------------------------------

/// **An insertion between the last flank base and the tract's first base, spelling
/// bases foreign to both, is genotyped by the tract path alone.**
///
/// The read's cigar anchors the two inserted `G`s on chr1:40 — the last base of
/// the generic region, SNP/indel territory. The ruling is that the event is the
/// tract's: the SNP/indel path must not open a locus for it (`open_record.rs`'s
/// `b6309954` rule), and the tract path must carry it in its own observation
/// (`ssr.rs`'s left-junction convention: the inserted run's last base `G` differs
/// from the flank base `C` beside it, so left-alignment cannot carry it into the
/// flank). End to end that means **exactly one locus with a non-reference allele:
/// the tract, two bases longer than its reference** — a run that genotyped the
/// insertion on both paths would show a second variant locus ending at 40, and a
/// run that dropped it on both would show none.
#[test]
fn a_foreign_insertion_at_the_tracts_first_base_is_genotyped_by_the_tract_path_alone() {
    let (_reference_dir, reference) = junction_reference();
    let sequence = [chr1_window(25, 40), b"GG".to_vec(), chr1_window(41, 69)].concat();
    let (_bam_dir, bam) = a_sample_of(
        |name| {
            junction_read(
                &name,
                &[(Kind::Match, 16), (Kind::Insertion, 2), (Kind::Match, 29)],
                sequence.clone(),
            )
        },
        "foreign_insertion.bam",
    );

    let called = called_over(&bam, &reference);

    let variants = variant_loci(&called);
    assert_eq!(
        variants.len(),
        1,
        "the insertion is genotyped exactly once, and the run called variants at {:?}",
        variants
            .iter()
            .map(|locus| locus.region)
            .collect::<Vec<_>>(),
    );
    let tract = variants[0];
    assert_eq!(
        (tract.region.start.get(), tract.region.end.get()),
        TRACT_SPAN,
        "the one variant locus is the tract's own span",
    );
    assert!(
        matches!(tract.alleles().kind(), LocusKind::Ssr(_)),
        "and it was called through the tract path, its candidates saying {:?}",
        tract.alleles().kind(),
    );
    assert_eq!(
        alternative_lengths(tract),
        vec![TRACT.len() + 2],
        "the tract's alternative carries the two inserted bases",
    );
    assert_called_homozygous_alternative(tract, "the tract locus");
    assert!(
        called
            .iter()
            .all(|locus| !(locus.region.start.get() <= 40 && 40 <= locus.region.end.get())),
        "no called locus claims chr1:40, the flank base the insertion is anchored on — \
         the loci are {:?}",
        called.iter().map(|locus| locus.region).collect::<Vec<_>>(),
    );
}

/// **The same adjudication for the event a tract actually produces**: an inserted
/// repeat unit (`AT`) anchored on the flank's last base — how a mapper spells a
/// read that crossed a tract one unit longer than the reference's. The inserted
/// run's last base `T` differs from the flank base `C`, so the junction convention
/// keeps it in the tract; the SNP/indel path refuses the anchor. One variant
/// locus, the tract, one unit longer.
#[test]
fn a_repeat_unit_insertion_anchored_on_the_flanks_last_base_is_the_tracts_alone() {
    let (_reference_dir, reference) = junction_reference();
    let sequence = [chr1_window(25, 40), b"AT".to_vec(), chr1_window(41, 69)].concat();
    let (_bam_dir, bam) = a_sample_of(
        |name| {
            junction_read(
                &name,
                &[(Kind::Match, 16), (Kind::Insertion, 2), (Kind::Match, 29)],
                sequence.clone(),
            )
        },
        "unit_insertion.bam",
    );

    let called = called_over(&bam, &reference);

    let variants = variant_loci(&called);
    assert_eq!(
        variants.len(),
        1,
        "the extra unit is genotyped exactly once, and the run called variants at {:?}",
        variants
            .iter()
            .map(|locus| locus.region)
            .collect::<Vec<_>>(),
    );
    let tract = variants[0];
    assert_eq!(
        (tract.region.start.get(), tract.region.end.get()),
        TRACT_SPAN,
        "the one variant locus is the tract's own span",
    );
    assert!(
        matches!(tract.alleles().kind(), LocusKind::Ssr(_)),
        "and it was called through the tract path, its candidates saying {:?}",
        tract.alleles().kind(),
    );
    assert_eq!(
        alternative_lengths(tract),
        vec![TRACT.len() + 2],
        "the tract's alternative is one repeat unit longer than its reference",
    );
    assert_called_homozygous_alternative(tract, "the tract locus");
}

/// **A deletion reaching across the flank–tract line becomes two loci, each
/// genotyping exactly its own bases.**
///
/// The reads delete chr1:39–42 — the flank's `GC` and the tract's first `AT`. The
/// ideal the design names for this shape is two variable loci, one per path, each
/// with its own indel allele: the SNP/indel path's deletion stops at the flank's
/// last base (its footprint is clamped at the region edge, so the tract's bases
/// are not its to call), and the tract path reads the same alignment as a tract
/// one unit short. Neither path sees the other's half, and together they spell
/// the whole event.
#[test]
fn a_deletion_across_flank_and_tract_is_two_loci_each_with_its_own_allele() {
    let (_reference_dir, reference) = junction_reference();
    let sequence = [chr1_window(25, 38), chr1_window(43, 69)].concat();
    let (_bam_dir, bam) = a_sample_of(
        |name| {
            junction_read(
                &name,
                &[(Kind::Match, 14), (Kind::Deletion, 4), (Kind::Match, 27)],
                sequence.clone(),
            )
        },
        "spanning_deletion.bam",
    );

    let called = called_over(&bam, &reference);

    let mut variants = variant_loci(&called);
    variants.sort_by_key(|locus| locus.region.start.get());
    assert_eq!(
        variants.len(),
        2,
        "one locus per path, and the run called variants at {:?}",
        variants
            .iter()
            .map(|locus| locus.region)
            .collect::<Vec<_>>(),
    );

    let (flank, tract) = (variants[0], variants[1]);

    assert!(
        matches!(flank.alleles().kind(), LocusKind::Generic),
        "the flank's half is called through the SNP/indel path, its candidates saying {:?}",
        flank.alleles().kind(),
    );
    assert!(
        flank.region.end.get() < TRACT_SPAN.0,
        "the SNP/indel locus stops before the tract begins — the tract's bases are not its \
         to call — and its span is {:?}",
        flank.region,
    );
    let flank_reference_len = flank.alleles().reference().len();
    assert_eq!(
        alternative_lengths(flank),
        vec![flank_reference_len - 2],
        "its one alternative deletes exactly the two flank bases, chr1:39–40",
    );
    assert_called_homozygous_alternative(flank, "the flank locus");

    assert_eq!(
        (tract.region.start.get(), tract.region.end.get()),
        TRACT_SPAN,
        "the tract's half is the tract's own span",
    );
    assert!(
        matches!(tract.alleles().kind(), LocusKind::Ssr(_)),
        "called through the tract path, its candidates saying {:?}",
        tract.alleles().kind(),
    );
    assert_eq!(
        alternative_lengths(tract),
        vec![TRACT.len() - 2],
        "its one alternative is the tract one repeat unit short — the deleted first unit",
    );
    assert_called_homozygous_alternative(tract, "the tract locus");
}
