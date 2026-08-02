//! **A reference and one sample's reads on disk, built from nothing** — the fixture
//! ng's generic-walk measurements run on, in one body instead of two.
//!
//! Not an example itself: cargo discovers `examples/*.rs` and `examples/*/main.rs`, so a
//! plain file in a subdirectory is only compiled where a target asks for it by
//! `#[path = "shared/synthetic_alignment.rs"] mod synthetic_alignment;` (from an example)
//! or `#[path = "../examples/shared/synthetic_alignment.rs"] …` (from a bench).
//!
//! # Why it is shared
//!
//! `benches/ng_generic_pileup_perf.rs` built this contig and this BAM to time the generic
//! locus generator; `examples/ng_generic_walk_probe.rs` needs the same pair to *test* the
//! walk it measures. Two copies of a seeded contig and a BAM writer would be two things
//! that can drift while both still pass — and the drift would be invisible, because the
//! bench reports a time and the probe reports a count, neither of which says which fixture
//! produced it. One body means a bench number and a probe assertion describe the same
//! reads.
//!
//! # What the fixture is, and what it deliberately is not
//!
//! A single contig of pseudo-random bases, tiled by fixed-length reads that each carry one
//! mismatch, written as a coordinate-sorted single-read-group BAM. It is **deterministic**:
//! the bases come from a fixed-seed generator and the read starts from an integer rule, so
//! two runs on one commit compare, and a run today compares with a run next month.
//!
//! It is **not** a realistic genome. Nothing here places a repeat, a variant, a duplicate or
//! an unmapped mate. What it must do is (a) vary — a homopolymer contig gives the
//! left-aligner and the walker's fold an unrepresentatively easy time — and (b) cover its
//! span exactly, so *every* reference position yields one locus and a walk that quietly
//! generated nothing is visible as a wrong count rather than as a fast run.

use std::path::{Path, PathBuf};

use tempfile::TempDir;

/// How much reference a synthetic sample spans, how long its reads are, and how deep they
/// lie on it.
///
/// A struct rather than three arguments because two of the three are base-pair counts and
/// would transpose silently.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SyntheticGeometry {
    /// Reference length, in base pairs. Also the number of loci a full walk must produce:
    /// the reads tile `[1, span]`.
    pub span: u64,
    /// Read length, in base pairs. Illumina-shaped at 150.
    pub read_len: u64,
    /// Fold coverage: how many reads lie over an average position. `span * coverage /
    /// read_len` reads are written, so at span 5,000 / read length 150 / coverage 10 that
    /// is 333 reads, not 10.
    pub coverage: u64,
}

impl SyntheticGeometry {
    /// The reads this geometry writes — the fixture's own definition of depth, and the
    /// floor an admission check holds a walk to.
    pub fn num_reads(&self) -> u64 {
        self.span * self.coverage / self.read_len
    }

    /// The widest gap between two consecutive read starts.
    ///
    /// The reads are spread evenly over the startable positions `[1, span - read_len + 1]`,
    /// so this is that range divided by the gaps between reads, rounded **up** — integer
    /// spacing puts the largest step at the ceiling. It has to stay at or below `read_len`
    /// or the reads stop tiling; [`SyntheticSample::build`] refuses a geometry where it
    /// does not.
    pub fn stride(&self) -> u64 {
        let startable = self.span - self.read_len;
        match self.num_reads() {
            0 | 1 => startable,
            reads => startable.div_ceil(reads - 1),
        }
    }
}

/// A reference FASTA and a sample BAM in a temporary directory, held together so the
/// directory outlives the paths into it.
///
/// It does not carry the geometry it was built from: every caller passes one in and still
/// holds it, and a second copy here is a second thing that can be read after the first has
/// moved on.
pub struct SyntheticSample {
    /// Deleted on drop, which is why the paths below are only valid while this is alive.
    _dir: TempDir,
    pub fasta: PathBuf,
    pub bam: PathBuf,
}

impl SyntheticSample {
    /// Write the reference and the reads for one geometry.
    ///
    /// Untimed by construction — a bench body takes this already built, and a test builds
    /// it once per case.
    pub fn build(geometry: SyntheticGeometry) -> Self {
        // The three fields used to be two file-level constants in the bench, which made a
        // read longer than its reference unrepresentable. They are public now, so the
        // shape a caller can build is wider than the shape this can serve — and without
        // this the failure is `attempt to subtract with overflow` inside a private helper
        // in debug, and a wrapped length reaching a slice in release.
        assert!(
            geometry.read_len >= 1 && geometry.read_len <= geometry.span,
            "a read must fit the reference it tiles: read_len {} against span {}",
            geometry.read_len,
            geometry.span,
        );
        assert!(
            geometry.num_reads() >= 1,
            "coverage {} over span {} at read length {} writes no reads at all",
            geometry.coverage,
            geometry.span,
            geometry.read_len,
        );
        // **The claim in this module's header, enforced instead of hoped for.** The reads
        // are spread evenly over the startable positions, so at low depth consecutive
        // starts can be further apart than a read is long and the tiling develops holes —
        // and a hole is invisible in the only place it shows up, as a *slightly* smaller
        // locus count. Measured: span 5,000, read length 150, coverage 1 gives 33 reads at
        // a stride of 152 and covers **4,950 of 5,000** positions.
        assert!(
            geometry.stride() <= geometry.read_len,
            "coverage {} over span {} writes {} reads at a stride of {}, which leaves holes \
             between reads {} bases long — the fixture's whole point is that it tiles",
            geometry.coverage,
            geometry.span,
            geometry.num_reads(),
            geometry.stride(),
            geometry.read_len,
        );
        let bases = synthetic_contig(geometry.span as usize);
        let dir = TempDir::new().expect("a temporary directory for the fixture");
        let fasta = dir.path().join("ref.fa");
        let bam = dir.path().join("sample.bam");
        write_fasta(&fasta, &bases);
        write_bam(&bam, bases.len(), &synthetic_reads(&bases, geometry));
        Self {
            _dir: dir,
            fasta,
            bam,
        }
    }
}

/// `len` deterministic bases from a fixed-seed linear congruential generator.
///
/// Not a real sequence and not meant to be. Seeded rather than `rand` so the same commit
/// gives the same contig on every run.
fn synthetic_contig(len: usize) -> Vec<u8> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            BASES[(state >> 33) as usize % BASES.len()]
        })
        .collect()
}

/// The base that is not this one — used to plant a mismatch.
fn other_base(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

/// Reads tiling `[1, span]` at approximately `coverage`×, each carrying **one** mismatch.
///
/// The mismatch is what makes the walker's fold do the work it exists for: a run of reads
/// matching the reference everywhere produces one row per locus and never exercises a
/// second bucket, a split row or a chain id. One per read, at a position that moves with
/// the read index, keeps every locus's row count small and realistic rather than turning
/// the fixture into a variant-caller stress test.
///
/// Starts ascend, which is both what a coordinate-sorted BAM needs and what the walker
/// requires of its input.
fn synthetic_reads(
    contig: &[u8],
    geometry: SyntheticGeometry,
) -> Vec<noodles_sam::alignment::RecordBuf> {
    use noodles_core::Position as CorePosition;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::{Flags, MappingQuality};
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

    // `coverage` is read just below, through `num_reads()`. Named rather than `..` so a
    // field added to the geometry — a mismatch rate, an indel rate — cannot slip past the
    // one function that lays reads down.
    let SyntheticGeometry {
        span,
        read_len,
        coverage: _,
    } = geometry;
    let last_start = span - read_len + 1;
    let num_reads = geometry.num_reads();
    let mut reads = Vec::with_capacity(num_reads as usize);
    for i in 0..num_reads {
        // Even spacing over the startable positions; several reads share a start once
        // `num_reads` exceeds them, which is what real depth looks like.
        //
        // The divisor is `num_reads - 1`, not `num_reads`, so the last read starts at
        // `last_start` and its final base is `span`. With `num_reads` the walk stopped
        // fifteen positions short of the span at 10× — invisible to a timing number, which
        // is why the checks downstream assert on the loci.
        let start = 1 + (i * (last_start - 1) / (num_reads - 1).max(1));
        let mut seq = contig[start as usize - 1..(start + read_len - 1) as usize].to_vec();
        let at = (i % read_len) as usize;
        seq[at] = other_base(seq[at]);
        reads.push(
            noodles_sam::alignment::RecordBuf::builder()
                .set_name(format!("r{i}").as_bytes())
                .set_reference_sequence_id(0)
                .set_flags(Flags::empty())
                .set_mapping_quality(MappingQuality::new(60).expect("60 is a mapping quality"))
                .set_alignment_start(
                    CorePosition::try_from(start as usize).expect("starts are one-based"),
                )
                .set_cigar(
                    [Op::new(Kind::Match, read_len as usize)]
                        .into_iter()
                        .collect(),
                )
                .set_sequence(Sequence::from(seq))
                .set_quality_scores(QualityScores::from(vec![40u8; read_len as usize]))
                .build(),
        );
    }
    reads
}

/// Write `chr1` as a one-line FASTA. The `.fai` is created by the reference read.
fn write_fasta(path: &Path, bases: &[u8]) {
    use std::io::Write as _;
    let mut file = std::fs::File::create(path).expect("the fixture's reference is writable");
    writeln!(file, ">chr1").expect("writing the FASTA header");
    file.write_all(bases).expect("writing the FASTA bases");
    writeln!(file).expect("terminating the FASTA");
}

/// Write a coordinate-sorted single-contig, single-read-group BAM.
fn write_bam(path: &Path, contig_len: usize, reads: &[noodles_sam::alignment::RecordBuf]) {
    use bstr::BString;
    use noodles_bam as bam;
    use noodles_sam as sam;
    use sam::alignment::io::Write as _;
    use sam::header::record::value::Map;
    use sam::header::record::value::map::header::Version;
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};
    use std::num::NonZero;

    let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
    hd.other_fields_mut()
        .insert(SORT_ORDER, BString::from("coordinate"));
    let sq = Map::<ReferenceSequence>::new(NonZero::new(contig_len).expect("a non-empty contig"));
    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut()
        .insert(SAMPLE, BString::from("sample0"));
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(b"chr1".to_vec(), sq)
        .add_read_group(b"rg0".to_vec(), rg)
        .build();

    let mut writer = bam::io::Writer::new(
        std::fs::File::create(path).expect("the fixture's alignment file is writable"),
    );
    writer
        .write_header(&header)
        .expect("writing the BAM header");
    for record in reads {
        writer
            .write_alignment_record(&header, record)
            .expect("writing a fixture read");
    }
    writer.try_finish().expect("finishing the BAM");
}
