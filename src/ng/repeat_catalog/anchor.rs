//! The whole catalog stack, end to end, on a real multi-contig FASTA.
//!
//! **The port anchor.** Everything else in this module is a unit test standing next to the
//! code it checks, handing a function rows it built by hand or bases it wrote to make a
//! point. This module drives what actually ships, in order: a **real multi-contig FASTA on
//! disk**, through the reference pass, through [`RepeatCatalogBuilder`] into a Parquet file,
//! back out through [`RepeatCatalog`], and into typed regions. Nothing here reaches inside.
//!
//! It exists to fail when the pieces stop fitting together, which no unit test can: the
//! FASTA reader, the observer seam the builder hangs off, the per-contig scan, the row
//! encoding, the file's own footer and page index, the reader's criteria check, the widening
//! a region read needs, and admission itself — all at once, on sequence nobody wrote to make
//! a point.
//!
//! **`#[cfg(test)]` in-crate, not `tests/`, for one reason:** the golden catalog's reader is
//! `pub(crate)` and production is frozen (`typed_regions.md`, Revision), so an out-of-crate
//! test cannot open the oracle at all. It still consumes only what a caller could.
//!
//! What it asserts:
//!
//! - **`.cat` parity** — the file reproduces the committed trf-mod-built golden catalog;
//! - the **partition invariant** — contiguous, non-overlapping, complete, maximal;
//! - **region-invariance** — asking for part of a contig chooses what you are shown, never
//!   what things are;
//! - the **edge cases** the spec lists: a tract at position 1, a repeat-free contig, and a
//!   tract at one contig's end abutting one at the next's start;
//! - the **tally** describes the regions that came out.
//!
//! # What used to be here and is not
//!
//! An earlier version of this file drove a windowed walk of the reference and asserted that
//! its window size changed nothing. **The catalog's builder scans a contig whole**
//! (`builder.rs`), so there is no window to be invariant to — a satellite of any length comes
//! out as one row, and the three things a windowed scan needs (a margin carried across each
//! chunk, a rule for which side a straddling detection belongs to, and a cap on the repeat
//! length it can promise to catch) do not exist here. The knob went with the walk.

use crate::ng::reference_info::{ReferenceInfo, ReferenceSource, read_reference_info_observing};
use crate::ng::region_typing::segment_criteria::{MinCopies, SsrSegmentCriteria};
use crate::ng::region_typing::{RegionKind, TypedRegion, TypedRegionConfig};
use crate::ng::repeat_catalog::criteria::StrRepeatCriteria;
use crate::ng::repeat_catalog::{ReadScope, RepeatCatalog, RepeatCatalogBuilder};
use crate::ng::tandem_repeat::{PeriodRange, ScanParams};
use crate::ng::types::{Bp, ContigId, GenomeRegion, Position};
use std::io::Write;
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------
// Fixtures on disk
// ---------------------------------------------------------------------

/// Write `contigs` as a FASTA and return where it went.
///
/// No `.fai`: the reference pass streams the file from the start and the catalog reader
/// never opens it, so nothing here seeks.
fn write_fasta(dir: &Path, contigs: &[(String, Vec<u8>)]) -> PathBuf {
    let path = dir.join("ref.fa");
    let mut file = std::fs::File::create(&path).expect("create the FASTA");
    for (name, bases) in contigs {
        writeln!(file, ">{name}").expect("header");
        for chunk in bases.chunks(60) {
            file.write_all(chunk).expect("bases");
            file.write_all(b"\n").expect("newline");
        }
    }
    path
}

/// The floors the anchor's catalogs are built at.
///
/// **Below the golden catalog's on every period**, which is what lets the parity test read
/// the file at the golden catalog's own settings: a reader asking for a floor beneath the
/// built one is refused, and the golden catalog's period-4 and period-5 floors are 3 where
/// the shipped catalog default is 4.
fn built_criteria() -> StrRepeatCriteria {
    StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            periods: PeriodRange::new(1, 6).expect("1..=6 is a valid period range"),
            min_copies: MinCopies::new([3, 3, 3, 3, 3, 3], 3),
            ..SsrSegmentCriteria::default()
        },
        ..StrRepeatCriteria::default()
    }
}

/// A FASTA on disk, its catalog beside it, and the reference the pass computed.
struct BuiltCatalog {
    /// Kept so the directory outlives the file paths taken from it.
    _dir: tempfile::TempDir,
    path: PathBuf,
    reference: ReferenceInfo,
}

/// Write `contigs`, run the reference pass with the builder attached, and finish the file.
fn build(contigs: &[(String, Vec<u8>)]) -> BuiltCatalog {
    std::fs::create_dir_all("tmp").expect("project-local scratch (CLAUDE.md: never /tmp)");
    let dir = tempfile::tempdir_in("tmp").expect("a scratch directory");
    let fasta = write_fasta(dir.path(), contigs);
    let path = dir.path().join("ref.fa.repeats.parquet");
    let mut builder = RepeatCatalogBuilder::create(&path, built_criteria(), ScanParams::default())
        .expect("the builder opens its output");
    let reference =
        read_reference_info_observing(ReferenceSource::Fasta { fasta, fai: None }, &mut builder)
            .expect("the reference pass runs");
    builder.finish(&reference).expect("the file is finished");
    BuiltCatalog {
        _dir: dir,
        path,
        reference,
    }
}

impl BuiltCatalog {
    /// The typed regions the file gives for `wanted`, at `criteria`.
    fn segments(&self, criteria: &StrRepeatCriteria, wanted: &[GenomeRegion]) -> Vec<TypedRegion> {
        let catalog = RepeatCatalog::open_checking_against_reference(&self.path, &self.reference)
            .expect("the catalog describes the reference it was built from");
        catalog
            .genome_segments(criteria, ReadScope::Regions(wanted))
            .expect("the policy is servable from this file")
            .map(|region| region.expect("the file reads cleanly"))
            .collect()
    }

    /// Every contig, end to end, as the region list the catalog is asked with.
    fn whole_reference(&self) -> Vec<GenomeRegion> {
        self.reference
            .contigs
            .iter()
            .enumerate()
            .map(|(index, contig)| GenomeRegion {
                contig: ContigId(index as u32),
                start: Position(1),
                end: Position(contig.length),
            })
            .collect()
    }
}

/// The committed golden reference, as `(name, bases)` — real sequence, and the same file
/// the golden `.cat` catalog was built from.
fn golden_contigs() -> Vec<(String, Vec<u8>)> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("tandem_repeat")
        .join("synthetic_ref.fa");
    let file = std::fs::File::open(path).expect("the committed golden reference");
    let mut reader = noodles_fasta::io::Reader::new(std::io::BufReader::new(file));
    reader
        .records()
        .map(|record| {
            let record = record.expect("a FASTA record");
            (
                String::from_utf8_lossy(record.name()).into_owned(),
                record.sequence().as_ref().to_vec(),
            )
        })
        .collect()
}

// ---------------------------------------------------------------------
// The invariant, stated here rather than borrowed
// ---------------------------------------------------------------------

/// **The partition invariant** (`typed_regions.md` §2.3), over one contig's slice of the
/// output: contiguous, non-overlapping, complete over `[1, contig_len]`, and maximal.
///
/// Written out again rather than shared with the unit tests, deliberately: this file is the
/// outside view, and an anchor that imports the thing it is anchoring proves less. One
/// property — *concatenating the regions reconstructs the contig, exactly.*
#[track_caller]
fn assert_partitions(regions: &[TypedRegion], contig: ContigId, contig_len: u64, case: &str) {
    assert!(
        !regions.is_empty(),
        "{case}: a non-empty contig has regions"
    );
    let mut expect = 1u64;
    let mut prev: Option<std::mem::Discriminant<RegionKind>> = None;
    for region in regions {
        assert_eq!(region.region.contig, contig, "{case}: contig");
        assert_eq!(
            region.region.start.get(),
            expect,
            "{case}: gap or overlap at {}",
            region.region.start.get()
        );
        assert!(
            region.region.end >= region.region.start,
            "{case}: empty region"
        );
        let kind = std::mem::discriminant(&region.kind);
        assert_ne!(
            Some(kind),
            prev,
            "{case}: two consecutive regions share a kind at {} — MAXIMALITY. For Generic \
             that is a correctness bug: the pileup mints loci INSIDE a generic region, so a \
             split run makes an indel across the join callable by neither half",
            region.region.start.get()
        );
        prev = Some(kind);
        expect = region.region.end.get() + 1;
    }
    assert_eq!(
        expect - 1,
        contig_len,
        "{case}: the partition must cover exactly [1, {contig_len}] — COMPLETENESS"
    );
}

fn kinds(regions: &[TypedRegion]) -> Vec<&'static str> {
    regions
        .iter()
        .map(|region| match &region.kind {
            RegionKind::SsrSegment(_) => "locus",
            RegionKind::SsrBundle { .. } => "bundle",
            RegionKind::Generic => "generic",
            RegionKind::Satellite => "satellite",
        })
        .collect()
}

fn on_contig(regions: &[TypedRegion], contig: ContigId) -> Vec<TypedRegion> {
    regions
        .iter()
        .filter(|region| region.region.contig == contig)
        .cloned()
        .collect()
}

// ---------------------------------------------------------------------
// The real reference, end to end
// ---------------------------------------------------------------------

/// **`.cat` parity through the whole stack** — the anchor's headline.
///
/// The file read at the golden catalog's settings must reproduce the golden catalog: every
/// golden locus is present, **or** absent *and* inside a satellite run, **or** absent *and*
/// within 15 bases of a contig's end, which is the one thing the file deliberately does not
/// hold (`repeat_catalog.md` §4.1). A strict subset otherwise, and that shape is earned by
/// the spec's ordering — the satellite cap applies to the *cleaned* coverage, after
/// classification, so the difference can only go one way.
///
/// The oracle is the committed **trf-mod-built** golden catalog: a different detector, a
/// different code path, nothing ng touched. Overlap matching, inherited from
/// `scanner_parity`, because the detector difference is characterised (±1–2 bp of
/// boundary/phase wobble) and is a yardstick rather than a confound.
///
/// The unit tests assert this too, from bases held in memory. What this adds is everything
/// between: the FASTA on disk, the streaming reference pass, the builder's per-contig scan,
/// the Parquet encoding, and the reader.
#[test]
fn the_catalog_reproduces_the_golden_catalog_through_the_shipping_stack() {
    use crate::ssr::catalog::io::CatalogReader;

    let fixture = |name: &str| {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("data")
            .join("tandem_repeat")
            .join(name)
    };
    let mut reader =
        CatalogReader::new(std::fs::File::open(fixture("golden.ssr_catalog.bed.gz")).unwrap())
            .unwrap();
    let cat_params = reader.header().params.clone();
    let golden = reader.read_all().unwrap();
    assert!(!golden.is_empty(), "the golden catalog must have loci");

    let contigs = golden_contigs();
    assert!(
        contigs.len() > 1,
        "the anchor walks a MULTI-contig reference"
    );
    let built = build(&contigs);

    // The golden catalog's settings, pinned explicitly rather than to whatever `Default`
    // is — or this starts failing the first time someone moves a floor and reads a *result*
    // as a bug. `min_purity` / `min_score` / `bundle_threshold` still come from the `.cat`
    // header; the period scope, the copy floors and the satellite cap are the catalog's
    // build settings (periods 2..=6, copy floors `[10,5,4,3,3,3]`, satellite cap 1 kb).
    let criteria = StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            min_purity: cat_params.min_purity,
            min_score: cat_params.min_score,
            bundle_threshold: u64::from(cat_params.bundle_threshold),
            periods: PeriodRange::new(2, 6).expect("2..=6 is a valid period range"),
            min_copies: MinCopies::new([10, 5, 4, 3, 3, 3], 3),
        },
        max_str_len_bp: Bp(1000),
        ..StrRepeatCriteria::default()
    };
    let regions = built.segments(&criteria, &built.whole_reference());

    // The partition holds on every contig of real sequence.
    for (index, (name, bases)) in contigs.iter().enumerate() {
        let contig = ContigId(index as u32);
        assert_partitions(
            &on_contig(&regions, contig),
            contig,
            bases.len() as u64,
            &format!("golden contig {name}"),
        );
    }

    // Present, or absent AND inside a satellite, or absent AND at a contig's very end.
    let name_of = |id: ContigId| contigs[id.get() as usize].0.clone();
    let ours: Vec<(String, u64, u64)> = regions
        .iter()
        .filter_map(|region| match &region.kind {
            RegionKind::SsrSegment(locus) => {
                Some((name_of(region.region.contig), locus.start(), locus.end()))
            }
            _ => None,
        })
        .collect();
    let satellites: Vec<(String, u64, u64)> = regions
        .iter()
        .filter(|region| matches!(region.kind, RegionKind::Satellite))
        .map(|region| {
            (
                name_of(region.region.contig),
                region.region.start.get(),
                region.region.end.get(),
            )
        })
        .collect();
    assert!(!ours.is_empty(), "the file must hold loci");

    let length_of = |name: &str| {
        contigs
            .iter()
            .find(|(contig, _)| contig == name)
            .map(|(_, bases)| bases.len() as u64)
            .expect("the golden catalog names contigs this reference has")
    };
    // Production's `Locus` is 0-based half-open, ng's 1-based inclusive: `[s, e)` is
    // `[s + 1, e]` (`typed_regions.md` §4).
    let overlaps =
        |a: &(String, u64, u64), b: &(String, u64, u64)| a.0 == b.0 && a.1 <= b.2 && b.1 <= a.2;
    let flank = built_criteria().min_flank_bp.get();
    let (mut missed, mut at_an_edge) = (Vec::new(), 0u64);
    for locus in &golden {
        let it = (
            locus.chrom().to_string(),
            u64::from(locus.start()) + 1,
            u64::from(locus.end()),
        );
        if ours.iter().any(|one| overlaps(&it, one)) || satellites.iter().any(|s| overlaps(&it, s))
        {
            continue;
        }
        if it.1 <= flank || it.2 + flank > length_of(&it.0) {
            at_an_edge += 1;
            continue;
        }
        missed.push(format!("{}:{}-{}", it.0, it.1, it.2));
    }
    assert!(
        missed.is_empty(),
        "every golden locus must be present, or absent AND inside a satellite run, or \
         absent AND within {flank} bases of a contig's end. At the golden catalog's \
         settings, a locus missing for any other reason is a machinery bug. Missing: \
         {missed:#?}"
    );
    // The edge exemption is stated so it stays visible: if it ever swallows most of the
    // catalog, the parity claim above has stopped meaning anything.
    assert!(
        at_an_edge * 20 < golden.len() as u64,
        "{at_an_edge} of the golden catalog's {} loci were excused for sitting within \
         {flank} bases of a contig's end — that exemption is supposed to be a handful",
        golden.len()
    );
}

/// **Region-invariance through the shipping stack** (`typed_regions.md` §2.5): asking for
/// part of a contig chooses what you are shown, never what things are.
///
/// Every base of every region a subset read returns is compared against the whole-reference
/// read of the same file. The span is chosen to cut through the middle of the sequence, not
/// to land tidily between features.
#[test]
fn a_region_subset_does_not_change_what_things_are() {
    let contigs = golden_contigs();
    let built = build(&contigs);
    let criteria = StrRepeatCriteria::from(&TypedRegionConfig::default());

    let whole = built.segments(&criteria, &built.whole_reference());

    let length = contigs[0].1.len() as u64;
    assert!(length > 800, "the fixture contig must be worth subsetting");
    let wanted = [GenomeRegion {
        contig: ContigId(0),
        start: Position(401),
        end: Position(800),
    }];
    let subset = built.segments(&criteria, &wanted);
    assert!(!subset.is_empty(), "the span must contain something");

    for region in &subset {
        for position in [region.region.start.get(), region.region.end.get()] {
            // A finding straddling the edge is emitted whole, so only ask about the bases
            // actually requested.
            if !(401..=800).contains(&position) {
                continue;
            }
            let truth = whole
                .iter()
                .find(|one| {
                    one.region.contig == region.region.contig
                        && one.region.contains(Position(position))
                })
                .expect("the whole-reference read covers every base");
            assert_eq!(
                std::mem::discriminant(&region.kind),
                std::mem::discriminant(&truth.kind),
                "base {position} is {:?} inside a region subset and {:?} without one",
                region.kind,
                truth.kind
            );
        }
    }
}

// ---------------------------------------------------------------------
// The edge cases the spec lists
// ---------------------------------------------------------------------

/// Aperiodic filler — **not** a homopolymer, which would be a period-1 tract and make these
/// fixtures pass for the wrong reason.
fn filler(n: usize) -> Vec<u8> {
    b"ACGTTGCAAGCTTGCA"
        .iter()
        .copied()
        .cycle()
        .take(n)
        .collect()
}

fn tract(copies: usize) -> Vec<u8> {
    b"AT".iter().copied().cycle().take(copies * 2).collect()
}

/// **The three edge cases the spec names, on one multi-contig reference** — and the third is
/// the one only a multi-contig reference can reach.
///
/// | contig | shape | what it pins |
/// |---|---|---|
/// | `at_start` | a tract at base 1, then filler | no room beside it at the CONTIG's start → `Generic`, and the partition still starts at 1 |
/// | `empty_of_repeats` | aperiodic filler only | **one** `Generic`, not a run of them (maximality) |
/// | `ends_with_tract` / `starts_with_tract` | a tract at one contig's end abutting one at the next's start | contigs are independent: neither tract borrows the other's flank, and nothing carries across the seam |
/// | `flanked` | **the control**: the same tract, flanks both sides | it IS a locus |
///
/// **The control is what makes the rest mean anything.** Three of these assert *"not a
/// locus"*, and a tract that was never admissible would satisfy all three for free.
/// `flanked` is the same `tract(10)` built by the same helper at the same settings: it comes
/// back a locus, so the absences above are the contig **ends** doing their job and not the
/// tract being unremarkable.
#[test]
fn the_edge_cases_hold_on_a_multi_contig_reference() {
    let mut at_start = tract(10);
    at_start.extend(filler(200));

    let repeat_free = filler(300);

    let mut ends_with = filler(200);
    ends_with.extend(tract(10));

    let mut starts_with = tract(10);
    starts_with.extend(filler(200));

    let mut flanked = filler(200);
    flanked.extend(tract(10));
    flanked.extend(filler(200));

    let contigs: Vec<(String, Vec<u8>)> = vec![
        ("at_start".to_string(), at_start.clone()),
        ("empty_of_repeats".to_string(), repeat_free.clone()),
        ("ends_with_tract".to_string(), ends_with.clone()),
        ("starts_with_tract".to_string(), starts_with.clone()),
        ("flanked".to_string(), flanked.clone()),
    ];
    let built = build(&contigs);
    let regions = built.segments(
        &StrRepeatCriteria::from(&TypedRegionConfig::default()),
        &built.whole_reference(),
    );

    for (index, (name, bases)) in contigs.iter().enumerate() {
        assert_partitions(
            &on_contig(&regions, ContigId(index as u32)),
            ContigId(index as u32),
            bases.len() as u64,
            name,
        );
    }

    // A tract at base 1 has no room beside it, so it is not a locus — and the partition
    // still tiles from base 1.
    let first = on_contig(&regions, ContigId(0));
    assert_eq!(first[0].region.start, Position(1));
    assert!(
        !kinds(&first).contains(&"locus"),
        "a tract at base 1 has nothing beside it: {:?}",
        kinds(&first)
    );

    // A repeat-free contig is exactly ONE Generic region. Not many.
    let second = on_contig(&regions, ContigId(1));
    assert_eq!(kinds(&second), vec!["generic"]);
    assert_eq!(second[0].region.end, Position(repeat_free.len() as u64));

    // **The seam.** A tract at one contig's end and one at the next's start are 20 bp apart
    // *in the file* and on different chromosomes: neither may borrow the other's flank,
    // bundle with it, or carry a coverage run across. Both are dropped for the same reason —
    // each abuts its own contig's end — and the two partitions are independent.
    let third = on_contig(&regions, ContigId(2));
    let fourth = on_contig(&regions, ContigId(3));
    assert!(
        !kinds(&third).contains(&"locus"),
        "the tract at this contig's END has nothing beyond it: {:?}",
        kinds(&third)
    );
    assert!(
        !kinds(&fourth).contains(&"locus"),
        "the tract at this contig's START has nothing before it: {:?}",
        kinds(&fourth)
    );
    assert!(
        !kinds(&third).contains(&"bundle") && !kinds(&fourth).contains(&"bundle"),
        "and they must not bundle with each other ACROSS the contig seam"
    );
    assert_eq!(
        third.last().unwrap().region.end,
        Position(ends_with.len() as u64),
        "the third contig's partition ends at its own last base"
    );
    assert_eq!(
        fourth[0].region.start,
        Position(1),
        "and the fourth's starts at its own first base"
    );

    // **The control.** The same tract, the same settings, flanks either side: a locus. So
    // every "not a locus" above is the contig's end doing its job, and not a tract that was
    // never going to be classified.
    let fifth = on_contig(&regions, ContigId(4));
    assert_eq!(
        kinds(&fifth),
        vec!["generic", "locus", "generic"],
        "the SAME tract, given flanks, is a locus"
    );
}

/// The running tally describes what came out, on the real reference — *no silent caps*,
/// checked against the regions themselves rather than against literals.
#[test]
fn the_counts_describe_the_read_of_a_real_reference() {
    let contigs = golden_contigs();
    let built = build(&contigs);
    let criteria = StrRepeatCriteria::from(&TypedRegionConfig::default());
    let wanted = built.whole_reference();

    let catalog = RepeatCatalog::open_checking_against_reference(&built.path, &built.reference)
        .expect("the catalog describes its reference");
    let mut segments = catalog
        .genome_segments(&criteria, ReadScope::Regions(&wanted))
        .expect("the policy is servable");
    let mut regions = Vec::new();
    for region in segments.by_ref() {
        regions.push(region.expect("the file reads cleanly"));
    }
    let counts = segments.counts();

    assert_eq!(counts.spans, wanted.len() as u64, "one span per contig");
    let count = |kind: &str| kinds(&regions).iter().filter(|k| **k == kind).count() as u64;
    assert_eq!(counts.ssr_loci, count("locus"));
    assert_eq!(counts.ssr_bundles, count("bundle"));
    assert_eq!(counts.generic, count("generic"));
    assert_eq!(counts.satellites, count("satellite"));
    assert!(counts.ssr_loci > 0, "real sequence, real loci");

    // Repeat coverage that yielded no locus is a *subset* of the repeat coverage, so it
    // cannot exceed the bases the read typed as repeat at all.
    let repeat_bp: u64 = regions
        .iter()
        .filter(|region| !matches!(region.kind, RegionKind::Generic))
        .map(|region| region.region.len())
        .sum();
    assert!(
        counts.repeat_bp_with_no_locus <= repeat_bp + counts.ssr_loci,
        "the no-locus gap ({}) cannot exceed the repeat coverage it is part of ({repeat_bp})",
        counts.repeat_bp_with_no_locus
    );
}
