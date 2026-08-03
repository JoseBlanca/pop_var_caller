//! Fixtures shared by this module's tests: a small multi-contig reference, a
//! header builder, and an indexed BAM on disk.
//!
//! They live in their own file because both the gate's tests and the region
//! query's need them, and a `#[cfg(test)] mod tests` block is private to its
//! own module — the alternative is two copies that drift.

use std::fs::File;
use std::num::NonZero;
use std::path::PathBuf;

use noodles_bam as bam;
use noodles_core::Position as RecordPosition;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::{Flags, MappingQuality};
use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
use sam::header::record::value::Map;
use sam::header::record::value::map::{ReadGroup, ReferenceSequence};
use tempfile::TempDir;

use crate::bam::index_preflight::preflight_alignment_indexes;
use crate::ng::read::filtering::ReadFilterCounts;
use crate::ng::read::input::read_groups::ReadGroupResolution;
use crate::ng::read::input::reference::OpenReference;
use crate::ng::reference_info::{ReferenceSource, read_reference_info};
use crate::ng::types::ReadGroupId;
use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};

/// The fixture reference: two contigs of **different lengths**, so a
/// permutation of an `@SQ` list is detectable on `name` and a re-labelling on
/// `length`.
pub(crate) const FIXTURE_CONTIGS: [(&str, usize); 2] = [("chr1", 100), ("chr2", 200)];

/// An all-`A` accessor over [`FIXTURE_CONTIGS`] — **named**, so its contig table is the one
/// the fixture files declare in their `@SQ`.
///
/// **Use this rather than `InMemoryRefSeq::from_contigs` for anything that opens a cursor.**
/// `from_contigs` names its contigs `contig0`, `contig1`, … , and every cursor fixture used it
/// until 2026-08-03 — so every one of them was handing `AlignmentFile::cursor` an accessor
/// whose table disagreed with the file's on **every name**. Nothing noticed, because the check
/// of the day fetched a zero-length window per contig and a window resolves whatever the
/// contig is called. B1 replaced that with a comparison of the two tables and 23 tests failed
/// at once; this is what they were fixed to.
pub(crate) fn fixture_reference_bases() -> crate::ng::ref_seq::InMemoryRefSeq {
    crate::ng::ref_seq::InMemoryRefSeq::from_named_contigs(
        FIXTURE_CONTIGS
            .iter()
            .map(|(name, length)| ((*name).to_string(), vec![b'A'; *length]))
            .collect(),
    )
}

/// The same, over [`BIG_FIXTURE_CONTIG`] — for the tests that need a contig larger than BAI's
/// 16 kb finest bin.
pub(crate) fn big_fixture_reference_bases() -> crate::ng::ref_seq::InMemoryRefSeq {
    crate::ng::ref_seq::InMemoryRefSeq::from_named_contigs(vec![(
        BIG_FIXTURE_CONTIG.0.to_string(),
        vec![b'A'; BIG_FIXTURE_CONTIG.1],
    )])
}

/// One `@RG` for a test header: its id, plus whichever tags the test is about.
/// A `None` tag is **omitted entirely**, which is how the "no `SM`" and "no
/// `LB`" inputs are built.
pub(crate) struct FixtureReadGroup<'a> {
    pub(crate) id: &'a str,
    pub(crate) sample: Option<&'a str>,
    pub(crate) library: Option<&'a str>,
    pub(crate) platform: Option<&'a str>,
}

impl<'a> FixtureReadGroup<'a> {
    /// A read group carrying only `SM` — what a test that is about neither the
    /// library nor the platform wants.
    pub(crate) fn new(id: &'a str, sample: Option<&'a str>) -> Self {
        Self {
            id,
            sample,
            library: None,
            platform: None,
        }
    }

    pub(crate) fn with_library(mut self, library: &'a str) -> Self {
        self.library = Some(library);
        self
    }

    pub(crate) fn with_platform(mut self, platform: &'a str) -> Self {
        self.platform = Some(platform);
        self
    }
}

/// A header builder taking `(sort order, contigs, read groups)`, so each test
/// states only the field it is about. `None` omits the tag entirely.
///
/// Read groups are `(id, SM)` pairs — the shape almost every test wants. A test
/// that is about the library or the platform builds its read groups as
/// [`FixtureReadGroup`]s and calls [`header_with_read_groups`] instead.
pub(crate) fn header(
    sort_order_value: Option<&str>,
    contigs: &[(&str, usize, Option<&str>)],
    read_groups: &[(&str, Option<&str>)],
) -> sam::Header {
    let read_groups: Vec<FixtureReadGroup<'_>> = read_groups
        .iter()
        .map(|(id, sample)| FixtureReadGroup::new(id, *sample))
        .collect();
    header_with_read_groups(sort_order_value, contigs, &read_groups)
}

/// [`header`], with read groups that can carry `LB` and `PL` as well as `SM`.
pub(crate) fn header_with_read_groups(
    sort_order_value: Option<&str>,
    contigs: &[(&str, usize, Option<&str>)],
    read_groups: &[FixtureReadGroup<'_>],
) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::{LIBRARY, PLATFORM, SAMPLE};
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let mut hd = Map::<sam::header::record::value::map::Header>::default();
    if let Some(value) = sort_order_value {
        hd.other_fields_mut()
            .insert(SORT_ORDER, value.as_bytes().into());
    }

    let mut builder = sam::Header::builder().set_header(hd);

    for (name, length, md5) in contigs {
        let mut sq = Map::<ReferenceSequence>::new(NonZero::new(*length).unwrap());
        if let Some(md5) = md5 {
            sq.other_fields_mut()
                .insert(MD5_CHECKSUM, md5.as_bytes().into());
        }
        builder = builder.add_reference_sequence(*name, sq);
    }

    for read_group in read_groups {
        let mut rg = Map::<ReadGroup>::default();
        for (tag, value) in [
            (SAMPLE, read_group.sample),
            (LIBRARY, read_group.library),
            (PLATFORM, read_group.platform),
        ] {
            if let Some(value) = value {
                rg.other_fields_mut().insert(tag, value.as_bytes().into());
            }
        }
        builder = builder.add_read_group(read_group.id, rg);
    }

    builder.build()
}

/// A header whose `@SQ` list is `contigs`, `SO:coordinate`, and one read group
/// naming `NA12878` — the shape a file has to have to get past the gate.
pub(crate) fn bam_header(contigs: &[(&str, usize, Option<&str>)]) -> sam::Header {
    header(Some("coordinate"), contigs, &[("rg1", Some("NA12878"))])
}

/// How a fixture file's records resolve: every fixture header declares exactly
/// one `@RG`, so every record is that one and no record's `RG` is read.
///
/// The identifier is `0` because a test opens one file and the run's table would
/// have minted `0` for its only read group. A test that cares which identifier
/// reaches a read builds the table with `build_read_groups` instead; one about
/// per-record resolution builds a `PerRecord` resolution itself.
pub(crate) fn fixture_read_group() -> ReadGroupResolution {
    ReadGroupResolution::Sole(ReadGroupId(0))
}

/// [`FIXTURE_CONTIGS`] in the `@SQ` shape, with no `M5` tags.
pub(crate) fn matching_contigs() -> Vec<(&'static str, usize, Option<&'static str>)> {
    FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| (*name, *length, None))
        .collect()
}

/// An [`OpenReference`] over [`FIXTURE_CONTIGS`] — the shape every open takes.
///
/// Hands back the run-scoped handle rather than the bare [`ReferenceInfo`]
/// because that is what a caller really holds: one reference, shared by every
/// file it opens. A test that wants the description alone asks it for
/// [`info()`](OpenReference::info).
///
/// `with_digests` picks the arm: the `Fasta` arm reads the genome and carries
/// real per-contig MD5s, while the `Fai` arm cannot, so its digests are `None`
/// and the MD5 half of reconciliation is a no-op. Tests that care about the
/// digest comparison need both. The `Fai` arm also carries no `fasta_path`, so
/// it has no bases either — a CRAM cannot be opened against it.
pub(crate) fn fixture_reference(with_digests: bool) -> (TempDir, OpenReference) {
    let specs: Vec<ContigSpec> = FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| ContigSpec {
            name: (*name).to_string(),
            length: *length as u64,
        })
        .collect();
    let (dir, fasta) = build_fasta(&specs).expect("build fasta");

    let source = if with_digests {
        ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        }
    } else {
        ReferenceSource::Fai(crate::ng::reference_info::sibling_fai_path(&fasta))
    };
    (
        dir,
        OpenReference::from(read_reference_info(source).expect("read reference")),
    )
}

/// **A contig long enough that the index has resolution**, for the tests that need it.
///
/// BAI's finest bins are 16 kb and a BGZF block is 64 kB, so a fixture smaller than that
/// resolves *every* region to the same single chunk — and a test on one cannot tell a reader
/// that positions from one that bounds. That is not hypothetical: the first version of the
/// cursor's differential oracle ran on a 100-base contig and passed with the bounding defect
/// in place.
pub(crate) const BIG_FIXTURE_CONTIG: (&str, usize) = ("chrBig", 200_000);

/// An [`OpenReference`] over [`BIG_FIXTURE_CONTIG`], from a `.fai`-only read.
pub(crate) fn big_fixture_reference() -> (TempDir, OpenReference) {
    let (dir, fasta) = build_fasta(&[ContigSpec {
        name: BIG_FIXTURE_CONTIG.0.to_string(),
        length: BIG_FIXTURE_CONTIG.1 as u64,
    }])
    .expect("build fasta");
    (
        dir,
        OpenReference::from(
            read_reference_info(ReferenceSource::Fai(
                crate::ng::reference_info::sibling_fai_path(&fasta),
            ))
            .expect("read reference"),
        ),
    )
}

/// [`BIG_FIXTURE_CONTIG`] in the `@SQ` shape.
pub(crate) fn big_contig_specs() -> Vec<(&'static str, usize, Option<&'static str>)> {
    vec![(BIG_FIXTURE_CONTIG.0, BIG_FIXTURE_CONTIG.1, None)]
}

/// Reads spread across [`BIG_FIXTURE_CONTIG`], 30 bases each starting every `stride`, which
/// is enough data to span several BGZF blocks and many index bins.
pub(crate) fn big_spread_of_reads(stride: usize) -> Vec<RecordBuf> {
    let mut records = Vec::new();
    let mut start = 1;
    while start + 30 < BIG_FIXTURE_CONTIG.1 {
        records.push(read_named_with_length(&format!("r{start}"), 0, start, 30));
        start += stride;
    }
    records
}

/// A 10 bp perfectly-matching read at `start` on `reference_sequence_id`.
///
/// `qname` matters for the region-query oracle, which identifies reads by name
/// to compare two independently-produced streams.
pub(crate) fn read_named(qname: &str, reference_sequence_id: usize, start: usize) -> RecordBuf {
    read_named_with_length(qname, reference_sequence_id, start, 10)
}

/// A [`read_named`] record carrying the `RG:Z` tag naming the read group it
/// belongs to.
///
/// Separate from [`read_named`] because tagging every fixture record would
/// change what the existing tests are about: a file that declares one `@RG` is
/// read without the tag being consulted at all, so an untagged record is the
/// **normal** input, not a degenerate one.
pub(crate) fn read_named_in_read_group(
    qname: &str,
    reference_sequence_id: usize,
    start: usize,
    read_group_id: &str,
) -> RecordBuf {
    read_named_with_length_in_read_group(qname, reference_sequence_id, start, 10, read_group_id)
}

/// [`read_named_in_read_group`], at a chosen length.
///
/// **A test that runs the real filter needs at least `DEFAULT_MIN_READ_LENGTH`
/// (30).** The 10 bp default is below it, so a filtered stream drops such reads
/// entirely — which looks exactly like a read-group resolution that returned
/// nothing.
pub(crate) fn read_named_with_length_in_read_group(
    qname: &str,
    reference_sequence_id: usize,
    start: usize,
    length: usize,
    read_group_id: &str,
) -> RecordBuf {
    use sam::alignment::record::data::field::Tag;
    use sam::alignment::record_buf::data::field::Value;

    let mut record = read_named_with_length(qname, reference_sequence_id, start, length);
    record.data_mut().insert(
        Tag::READ_GROUP,
        Value::String(read_group_id.as_bytes().to_vec().into()),
    );
    record
}

pub(crate) fn read_named_with_length(
    qname: &str,
    reference_sequence_id: usize,
    start: usize,
    length: usize,
) -> RecordBuf {
    RecordBuf::builder()
        .set_name(qname.as_bytes())
        .set_reference_sequence_id(reference_sequence_id)
        // Explicit, because `RecordBuf`'s default flags are `UNMAPPED` — a
        // fixture that left them alone would be silently dropped by filter #1
        // in any test that runs the real filter.
        .set_flags(Flags::empty())
        .set_mapping_quality(MappingQuality::new(60).expect("mapq in range"))
        .set_alignment_start(RecordPosition::try_from(start).unwrap())
        .set_cigar([Op::new(Kind::Match, length)].into_iter().collect())
        .set_sequence(Sequence::from(vec![b'A'; length]))
        .set_quality_scores(QualityScores::from(vec![30u8; length]))
        .build()
}

/// Write a BAM holding `records` and build its index beside it.
///
/// Returns the `TempDir` as well as the path: bind it, or the directory is
/// removed the moment the call returns and the handle points at nothing.
pub(crate) fn indexed_bam(header: &sam::Header, records: &[RecordBuf]) -> (TempDir, PathBuf) {
    indexed_named_bam(header, records, "sample.bam")
}

/// [`indexed_bam`] under a chosen file name.
///
/// **Two fixture files opened together need different names.** A fixture header
/// declares no `LB`, so a read group's library name is synthesized from its
/// sample, its `@RG ID` and its file's name — and two files agreeing on all
/// three are indistinguishable, which the read-group pre-pass rejects. Real
/// inputs opened together differ in name for the same reason.
pub(crate) fn indexed_named_bam(
    header: &sam::Header,
    records: &[RecordBuf],
    file_name: &str,
) -> (TempDir, PathBuf) {
    let (dir, path) = named_bam(header, records, file_name);
    preflight_alignment_indexes(std::slice::from_ref(&path), true).expect("build index");
    (dir, path)
}

/// The same BAM with **no** index beside it — for the gate's index check.
pub(crate) fn unindexed_bam(header: &sam::Header, records: &[RecordBuf]) -> (TempDir, PathBuf) {
    named_bam(header, records, "sample.bam")
}

/// An unindexed BAM written under a chosen file name.
///
/// The name is a test input in its own right: a synthesized library name is
/// built from the file's name, so two files that share one — the same name in
/// different directories — are what makes those names collide.
pub(crate) fn named_bam(
    header: &sam::Header,
    records: &[RecordBuf],
    file_name: &str,
) -> (TempDir, PathBuf) {
    let dir = TempDir::new().expect("tempdir");
    let path = dir.path().join(file_name);

    let mut writer = bam::io::Writer::new(File::create(&path).expect("create bam"));
    writer.write_header(header).expect("write header");
    for record in records {
        writer
            .write_alignment_record(header, record)
            .expect("write record");
    }
    writer.try_finish().expect("finish");

    (dir, path)
}

/// The same records written as a **CRAM**, with a `.crai` beside it, over the
/// fixture reference. Returns the CRAM's dir and path plus the FASTA's dir and
/// path — a CRAM cannot be decoded without the reference, so its `ReferenceInfo`
/// has to come from the `Fasta` arm.
///
/// **`records` must all be on one contig.** Not a choice: `noodles_cram::fs::index`
/// decodes *multi-reference* slices with `fasta::Repository::default()` — an
/// empty repository, marked `// TODO` in noodles 0.93 (`src/fs/index.rs:137`) —
/// so indexing a CRAM whose reads are stored as differences from the reference
/// panics with "invalid reference sequence name". Reads spanning two contigs
/// land in one slice, which is exactly that case. Porting noodles' indexer with
/// a real repository is not possible from outside the crate:
/// `ReferenceSequenceContext` and `Slice::header` are private.
///
/// Production's CRAM fixtures are single-contig too, which is why this has not
/// bitten before. The `.crai` contig walk is covered instead by hand-built
/// indexes, which need no file at all: the grouping in `open_bam`'s tests and
/// the positioning in `aligned_reads_reader::cram`'s.
pub(crate) fn indexed_cram(records: &[RecordBuf]) -> (TempDir, PathBuf, TempDir, PathBuf) {
    indexed_cram_declaring(records, &[("rg1", Some("NA12878"))])
}

/// [`indexed_cram`], with a chosen set of `@RG` records.
///
/// A CRAM declaring several read groups is what makes the CRAM sources resolve
/// each record's own `RG` rather than take the file's only one — a path the
/// single-`@RG` fixture cannot reach at all, and where a stamp taken from the
/// wrong record would be silent.
pub(crate) fn indexed_cram_declaring(
    records: &[RecordBuf],
    read_groups: &[(&str, Option<&str>)],
) -> (TempDir, PathBuf, TempDir, PathBuf) {
    debug_assert!(
        records
            .iter()
            .all(|record| record.reference_sequence_id() == Some(0)),
        "the CRAM fixture must stay on one contig — see this function's docs"
    );
    use crate::pileup::per_sample::cram_files::{HeaderOverrides, build_cram};

    let specs: Vec<ContigSpec> = FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| ContigSpec {
            name: (*name).to_string(),
            length: *length as u64,
        })
        .collect();
    let (fasta_dir, fasta) = build_fasta(&specs).expect("build fasta");
    let (cram_dir, cram_path) = build_cram(
        &fasta,
        &specs,
        &HeaderOverrides {
            read_groups: read_groups
                .iter()
                .map(|(id, sample)| ((*id).to_string(), sample.map(str::to_string)))
                .collect(),
            ..HeaderOverrides::default()
        },
        records,
    )
    .expect("build cram");

    let index = noodles_cram::fs::index(&cram_path).expect(
        "noodles can index a single-reference CRAM; see the doc above for why \
         the fixture must stay single-contig",
    );
    let crai_path = PathBuf::from(format!("{}.crai", cram_path.display()));
    noodles_cram::crai::fs::write(&crai_path, &index).expect("write crai");

    (cram_dir, cram_path, fasta_dir, fasta)
}

/// A single-contig CRAM over a **long** contig, with enough reads to fill
/// several containers.
///
/// noodles writes 10240 records per container, so a fixture that stays under
/// that produces one container and one `.crai` entry — which exercises the
/// container decode but *none* of the `.crai` walk: not the multi-entry loop,
/// not the container-level early stop, not the span skip. Returns the CRAM, its
/// FASTA, and the contig length.
pub(crate) fn multi_container_cram(
    contig_length: usize,
    read_count: usize,
) -> (TempDir, PathBuf, TempDir, PathBuf) {
    use crate::pileup::per_sample::cram_files::{HeaderOverrides, build_cram};

    let specs = vec![ContigSpec {
        name: "chr1".to_string(),
        length: contig_length as u64,
    }];
    let (fasta_dir, fasta) = build_fasta(&specs).expect("build fasta");

    // Spread the reads evenly along the contig and keep them in coordinate
    // order, so successive containers cover successive stretches — which is
    // what makes the container-level early stop observable.
    let step = (contig_length - 40) / read_count.max(1);
    let records: Vec<RecordBuf> = (0..read_count)
        .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * step.max(1), 30))
        .collect();

    let (cram_dir, cram_path) = build_cram(
        &fasta,
        &specs,
        &HeaderOverrides {
            read_groups: vec![("rg1".to_string(), Some("NA12878".to_string()))],
            ..HeaderOverrides::default()
        },
        &records,
    )
    .expect("build cram");

    let index = noodles_cram::fs::index(&cram_path).expect("index a single-reference CRAM");
    let crai_path = PathBuf::from(format!("{}.crai", cram_path.display()));
    noodles_cram::crai::fs::write(&crai_path, &index).expect("write crai");

    (cram_dir, cram_path, fasta_dir, fasta)
}

/// One 10 bp read at the start of the first contig — enough to make a file
/// non-empty for tests that never read it.
pub(crate) fn one_read() -> Vec<RecordBuf> {
    vec![read_named("read-1", 0, 1)]
}

/// The fixtures test themselves, because the read-group tests that use them
/// cannot tell "the parser ignored the tag" from "the fixture never wrote it".
/// A silently-absent `LB` would make a "missing library" test pass for the
/// wrong reason.
mod tests {
    use super::*;

    #[test]
    fn a_read_group_carries_every_tag_it_was_given() {
        use sam::header::record::value::map::read_group::tag::{LIBRARY, PLATFORM, SAMPLE};

        let header = header_with_read_groups(
            Some("coordinate"),
            &matching_contigs(),
            &[FixtureReadGroup::new("rg1", Some("NA12878"))
                .with_library("lib-A")
                .with_platform("ILLUMINA")],
        );

        let (_, rg) = header
            .read_groups()
            .first()
            .expect("the header has one read group");
        let tag = |tag| rg.other_fields().get(&tag).map(|v| v.to_vec());
        assert_eq!(tag(SAMPLE).as_deref(), Some(&b"NA12878"[..]));
        assert_eq!(tag(LIBRARY).as_deref(), Some(&b"lib-A"[..]));
        assert_eq!(tag(PLATFORM).as_deref(), Some(&b"ILLUMINA"[..]));
    }

    /// The omissions matter as much as the values: an absent tag must be absent
    /// from the header, not present and empty.
    #[test]
    fn an_untagged_read_group_carries_nothing() {
        use sam::header::record::value::map::read_group::tag::{LIBRARY, PLATFORM, SAMPLE};

        let header = header(Some("coordinate"), &matching_contigs(), &[("rg1", None)]);

        let (_, rg) = header
            .read_groups()
            .first()
            .expect("the header has one read group");
        for tag in [SAMPLE, LIBRARY, PLATFORM] {
            assert!(rg.other_fields().get(&tag).is_none(), "{tag:?} is absent");
        }
    }

    #[test]
    fn a_tagged_record_names_its_read_group() {
        use sam::alignment::record::data::field::Tag;
        use sam::alignment::record_buf::data::field::Value;

        let record = read_named_in_read_group("read-1", 0, 1, "rg2");
        assert!(matches!(
            record.data().get(&Tag::READ_GROUP),
            Some(Value::String(id)) if id == "rg2"
        ));
        assert!(
            read_named("read-1", 0, 1)
                .data()
                .get(&Tag::READ_GROUP)
                .is_none(),
            "the untagged fixture stays untagged"
        );
    }
}

/// The one read group's tally, for a fixture file that declares exactly one.
///
/// Asserting the count rather than taking the first entry is the point: a
/// fixture that quietly grew a second read group would otherwise have half its
/// drops silently ignored by whatever test called this.
pub(crate) fn only_tally(
    counts: &[(Option<crate::ng::types::ReadGroupId>, ReadFilterCounts)],
) -> ReadFilterCounts {
    assert_eq!(
        counts.len(),
        1,
        "this fixture is expected to meet exactly one read group"
    );
    counts[0].1.clone()
}
