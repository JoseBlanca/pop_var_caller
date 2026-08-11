//! The catalog's north-star test: **the segmentation derived from the file equals the one a
//! live scan produces**, region for region, at several policies.
//!
//! This is what the whole design rests on (spec §5.1, §10.1). Every other test says a part
//! behaves; this one says the parts compose into the same answer `partition_resident` gives
//! — which is the claim a consumer relies on when it stops opening the FASTA.
//!
//! **One stated exception, and it is the file's flank floor.** Step 3 drops a locus only
//! when its flank clamps to zero, so a live scan keeps tracts 1 to 14 bases from a contig
//! end that the catalog never recorded. The comparison therefore runs at a reader flank of
//! 15 bp or more, which every calling policy already satisfies, and the tests below assert
//! that the difference is confined to the contig edges rather than assuming it.

use std::io::Write;
use std::path::{Path, PathBuf};

use pop_var_caller::ng::reference_info::{
    ReferenceInfo, ReferenceSource, read_reference_info_observing,
};
use pop_var_caller::ng::region_typing::segment_criteria::{MinCopies, SsrSegmentCriteria};
use pop_var_caller::ng::region_typing::{
    RegionKind, TypedRegion, TypedRegionConfig, partition_resident,
};
use pop_var_caller::ng::repeat_catalog::criteria::{
    CATALOG_MIN_COPIES, CATALOG_MIN_COPIES_BEYOND_TABLE,
};
use pop_var_caller::ng::repeat_catalog::{
    ReadScope, RepeatCatalog, RepeatCatalogBuilder, StrRepeatCriteria,
};
use pop_var_caller::ng::tandem_repeat::{PeriodRange, ScanParams};
use pop_var_caller::ng::types::{Bp, ContigId, GenomeRegion, Position};

/// A region covering one contig whole — what "just this chromosome" looks like now that the
/// read surface speaks regions.
fn whole(contig: ContigId) -> GenomeRegion {
    GenomeRegion {
        contig,
        start: Position(1),
        end: Position(u64::MAX),
    }
}

/// Filler with no tandem structure of its own at periods 1..=6.
fn filler(len: usize) -> String {
    const CYCLE: &[u8] = b"ACGTTGCAAGGTCCAT";
    (0..len).map(|i| CYCLE[i % CYCLE.len()] as char).collect()
}

fn write_fasta(dir: &Path, contigs: &[(&str, String)]) -> PathBuf {
    let path = dir.join("ref.fa");
    let mut file = std::fs::File::create(&path).expect("create");
    for (name, seq) in contigs {
        writeln!(file, ">{name}").expect("header");
        for chunk in seq.as_bytes().chunks(60) {
            file.write_all(chunk).expect("bases");
            file.write_all(b"\n").expect("newline");
        }
    }
    path
}

/// Build a catalog at the catalog's own (permissive) settings and return where it went,
/// with the reference the pass computed.
fn build_catalog(dir: &Path, contigs: &[(&str, String)]) -> (PathBuf, ReferenceInfo, Vec<String>) {
    let fasta = write_fasta(dir, contigs);
    let catalog_path = dir.join("ref.fa.repeats.parquet");
    let mut builder = RepeatCatalogBuilder::create(
        &catalog_path,
        StrRepeatCriteria::default(),
        ScanParams::default(),
    )
    .expect("builder");
    let reference =
        read_reference_info_observing(ReferenceSource::Fasta { fasta, fai: None }, &mut builder)
            .expect("the reference pass runs");
    builder.finish(&reference).expect("finish");
    let names = reference.contigs.iter().map(|c| c.name.clone()).collect();
    (catalog_path, reference, names)
}

/// The live scan, at the same policy the catalog is read with.
fn scanned(
    seq: &str,
    chrom: &str,
    contig: ContigId,
    criteria: &StrRepeatCriteria,
) -> Vec<TypedRegion> {
    let config = TypedRegionConfig {
        criteria: criteria.classification.clone(),
        max_str_len: criteria.max_str_len_bp,
        ..TypedRegionConfig::default()
    };
    partition_resident(chrom, contig, seq.as_bytes(), &config)
}

/// The same segmentation, derived from the file.
fn derived(
    catalog_path: &Path,
    reference: &ReferenceInfo,
    contig: ContigId,
    criteria: &StrRepeatCriteria,
) -> Vec<TypedRegion> {
    let catalog =
        RepeatCatalog::open_checking_against_reference(catalog_path, reference).expect("opens");
    catalog
        .genome_segments(criteria, ReadScope::Regions(&[whole(contig)]))
        .expect("the policy is servable")
        .map(|r| r.expect("a region"))
        .collect()
}

/// A tract nearer a contig end than the catalog's 15 bp flank floor is the one place the
/// file holds less than a scan (spec §5.1). Drop those regions from both sides before
/// comparing, and report how many were dropped so a test can assert the exception is small
/// and edge-confined rather than a general disagreement.
fn away_from_contig_edges(
    regions: &[TypedRegion],
    contig_len: u64,
    flank: u64,
) -> Vec<TypedRegion> {
    regions
        .iter()
        .filter(|r| r.region.start.get() > flank && r.region.end.get() + flank <= contig_len)
        .cloned()
        .collect()
}

/// The policies the comparison runs at: the catalog's own, then ones differing from it on
/// every bounded axis (a higher copy floor, a narrower period range) and on the unbounded
/// ones (purity, satellite cap).
fn policies() -> Vec<(&'static str, StrRepeatCriteria)> {
    let base = StrRepeatCriteria::default();

    let calling_floors = StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            min_copies: MinCopies::default(),
            ..base.classification.clone()
        },
        ..base.clone()
    };

    let narrow_periods = StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            periods: PeriodRange::new(2, 4).expect("valid"),
            ..base.classification.clone()
        },
        ..base.clone()
    };

    let strict_purity = StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            min_purity: 0.95,
            ..base.classification.clone()
        },
        ..base.clone()
    };

    let small_satellite_cap = StrRepeatCriteria {
        max_str_len_bp: Bp(100),
        ..base.clone()
    };

    let higher_floor_one_period = StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            min_copies: MinCopies::new(
                {
                    let mut floors = CATALOG_MIN_COPIES;
                    floors[2] = 9; // period 3 only
                    floors
                },
                CATALOG_MIN_COPIES_BEYOND_TABLE,
            ),
            ..base.classification.clone()
        },
        ..base.clone()
    };

    vec![
        ("the catalog's own settings", base),
        ("the calling copy floors", calling_floors),
        ("a narrower period range", narrow_periods),
        ("a stricter purity floor", strict_purity),
        ("a smaller satellite cap", small_satellite_cap),
        ("a higher floor at period 3", higher_floor_one_period),
    ]
}

/// A reference carrying the structures that make classification interesting: clean tracts
/// of several periods, two tracts close enough to bundle, a long array, an impure tract,
/// and tracts at both contig ends.
fn fixture() -> Vec<(&'static str, String)> {
    let clean = format!(
        "{}{}{}{}{}{}{}",
        filler(200),
        "CAG".repeat(10),
        filler(200),
        "AT".repeat(14),
        filler(200),
        "A".repeat(12),
        filler(200),
    );

    // Two tracts 10 bp apart — inside the 15 bp bundle radius, so they cluster.
    let bundled = format!(
        "{}{}{}{}{}",
        filler(200),
        "GATA".repeat(8),
        filler(10),
        "TCTC".repeat(8),
        filler(200),
    );

    // A 1.2 kb array, over the calling satellite cap of 100 and under the catalog's 500.
    let array = format!("{}{}{}", filler(200), "AT".repeat(600), filler(200));

    // An interrupted tract whose purity lands **between** the two policies' floors: three
    // clean copies then one broken one, three times over — 33 of 36 bases match a perfect
    // `CAG` tiling, so 0.92. A single interruption is not enough; the earlier fixture scored
    // 0.96 and the strict-purity policy discriminated nothing (see
    // `the_fixture_drives_the_purity_gate_and_the_satellite_cap`).
    let impure = format!("{}{}{}", filler(200), "CAGCAGCAGCTG".repeat(3), filler(200));

    // Tracts hard against both contig ends — the stated exception's own fixture.
    let edges = format!("{}{}{}", "CAG".repeat(10), filler(400), "CAG".repeat(10));

    // A tract 5 bases from the contig's end: **a locus to a live scan and absent from the
    // file**, because the catalog stores nothing with less than 15 bases beside it. Distinct
    // from `chr_edges`, whose tracts abut the very first and last base and so are rejected by
    // a live scan too. This is the contig that makes the difference measurable instead of
    // hypothetical (`the_tally_from_the_file_matches_the_walks`).
    let near_edge = format!("{}{}{}", filler(400), "CAG".repeat(10), filler(5));

    // 180 bases: a locus under the catalog's 500 bp satellite cap and a satellite under the
    // calling one of 100, which is what makes the cap a knob this fixture can move.
    let mid_array = format!("{}{}{}", filler(200), "AT".repeat(90), filler(200));

    vec![
        ("chr_clean", clean),
        ("chr_bundled", bundled),
        ("chr_array", array),
        ("chr_mid_array", mid_array),
        ("chr_impure", impure),
        ("chr_edges", edges),
        ("chr_near_edge", near_edge),
    ]
}

/// **The differential.** Every contig, every policy: what the file says equals what a scan
/// says, away from the contig edges the file deliberately omits.
#[test]
fn the_derived_segmentation_equals_the_scanned_one_at_every_policy() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let (catalog_path, reference, names) = build_catalog(dir.path(), &contigs);
    assert_eq!(names.len(), contigs.len());

    for (index, (name, seq)) in contigs.iter().enumerate() {
        let contig = ContigId(index as u32);
        for (label, criteria) in policies() {
            let flank = criteria.min_flank_bp.get();
            let live = scanned(seq, name, contig, &criteria);
            let from_file = derived(&catalog_path, &reference, contig, &criteria);

            let live_inner = away_from_contig_edges(&live, seq.len() as u64, flank);
            let file_inner = away_from_contig_edges(&from_file, seq.len() as u64, flank);

            assert_eq!(
                file_inner, live_inner,
                "{name} under {label}: the catalog and a live scan disagree\n\
                 from the file: {from_file:#?}\nfrom a scan: {live:#?}"
            );
        }
    }
}

/// The loci are what a consumer actually asks for, so compare them by name rather than
/// letting a generic-stretch coincidence carry the test.
#[test]
fn the_derived_loci_are_the_scanned_loci() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let (catalog_path, reference, _) = build_catalog(dir.path(), &contigs);

    let mut compared = 0usize;
    for (index, (name, seq)) in contigs.iter().enumerate() {
        let contig = ContigId(index as u32);
        let criteria = StrRepeatCriteria::default();
        let flank = criteria.min_flank_bp.get();

        let loci = |regions: &[TypedRegion]| -> Vec<(u64, u64, String, f32)> {
            away_from_contig_edges(regions, seq.len() as u64, flank)
                .iter()
                .filter_map(|r| match &r.kind {
                    RegionKind::SsrSegment(s) => Some((
                        r.region.start.get(),
                        r.region.end.get(),
                        String::from_utf8_lossy(s.motif().as_bytes()).into_owned(),
                        s.purity_fraction(),
                    )),
                    _ => None,
                })
                .collect()
        };

        let live = loci(&scanned(seq, name, contig, &criteria));
        let from_file = loci(&derived(&catalog_path, &reference, contig, &criteria));
        compared += live.len();
        assert_eq!(from_file, live, "{name}: loci differ");
    }
    assert!(
        compared >= 4,
        "the fixture must actually produce loci, or the comparison proves nothing ({compared})"
    );
}

/// **The fixture must drive every gate the comparison claims to test**, or the differential
/// passes for the wrong reason. Found by mutation: removing the purity floor from the
/// derived path left every test green, because nothing in an earlier fixture produced a
/// locus whose purity sat between the two policies' floors.
#[test]
fn the_fixture_drives_the_purity_gate_and_the_satellite_cap() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let (catalog_path, reference, _) = build_catalog(dir.path(), &contigs);

    let base = StrRepeatCriteria::default();
    let strict = StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            min_purity: 0.95,
            ..base.classification.clone()
        },
        ..base.clone()
    };

    let mut purities = Vec::new();
    let mut satellites_at_small_cap = 0usize;
    let mut satellites_at_large_cap = 0usize;
    let small_cap = StrRepeatCriteria {
        max_str_len_bp: Bp(100),
        ..base.clone()
    };

    for index in 0..contigs.len() {
        let contig = ContigId(index as u32);
        for region in derived(&catalog_path, &reference, contig, &base) {
            if let RegionKind::SsrSegment(s) = &region.kind {
                purities.push(s.purity_fraction());
            }
        }
        satellites_at_large_cap += derived(&catalog_path, &reference, contig, &base)
            .iter()
            .filter(|r| matches!(r.kind, RegionKind::Satellite))
            .count();
        satellites_at_small_cap += derived(&catalog_path, &reference, contig, &small_cap)
            .iter()
            .filter(|r| matches!(r.kind, RegionKind::Satellite))
            .count();
    }

    let between = purities.iter().filter(|p| (0.8..0.95).contains(*p)).count();
    assert!(
        between >= 1,
        "no locus has a purity between the two policies' floors, so the strict-purity policy \
         discriminates nothing. Purities seen: {purities:?}"
    );

    // And the strict policy must actually drop one of them.
    let strict_loci: usize = (0..contigs.len())
        .map(|i| {
            derived(&catalog_path, &reference, ContigId(i as u32), &strict)
                .iter()
                .filter(|r| matches!(r.kind, RegionKind::SsrSegment(_)))
                .count()
        })
        .sum();
    assert!(
        strict_loci < purities.len(),
        "the strict purity floor kept every locus ({strict_loci} of {}), so it tests nothing",
        purities.len()
    );

    assert!(
        satellites_at_small_cap > satellites_at_large_cap,
        "the satellite cap must move something: {satellites_at_small_cap} at 100 bp vs \
         {satellites_at_large_cap} at 500 bp"
    );
}

/// The per-stratum tally read from the file equals one taken by scanning the reference
/// directly (spec §10.3). This is the number the pre-pass reweights by, so an error here
/// propagates into a diversity estimate rather than into a crash.
#[test]
fn the_strata_tally_matches_a_direct_count() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let (catalog_path, reference, _) = build_catalog(dir.path(), &contigs);
    let catalog =
        RepeatCatalog::open_checking_against_reference(&catalog_path, &reference).expect("opens");
    let criteria = StrRepeatCriteria::default();

    // Counted by scanning, exactly as a caller without a catalog would have to.
    let mut scanned_counts: std::collections::BTreeMap<(u8, u64), u64> = Default::default();
    for (index, (name, seq)) in contigs.iter().enumerate() {
        for region in scanned(seq, name, ContigId(index as u32), &criteria) {
            if let RegionKind::SsrSegment(locus) = &region.kind {
                let period = locus.motif().period() as u8;
                let copies = (locus.end() - locus.start() + 1) / u64::from(period);
                *scanned_counts.entry((period, copies)).or_insert(0) += 1;
            }
        }
    }

    // Counted from the file.
    let from_file = catalog
        .count_loci_per_stratum(&criteria, ReadScope::WholeReference)
        .expect("servable");

    // The one difference the file admits to: loci within the flank floor of a contig's end.
    // Count those separately and require the rest to agree exactly.
    let mut edge_loci = 0u64;
    for (index, (name, seq)) in contigs.iter().enumerate() {
        let len = seq.len() as u64;
        let flank = criteria.min_flank_bp.get();
        for region in scanned(seq, name, ContigId(index as u32), &criteria) {
            if matches!(region.kind, RegionKind::SsrSegment(_))
                && (region.region.start.get() <= flank || region.region.end.get() + flank > len)
            {
                edge_loci += 1;
            }
        }
    }

    let scanned_total: u64 = scanned_counts.values().sum();
    assert!(
        scanned_total > 0,
        "the fixture must produce loci for this comparison to mean anything"
    );
    assert_eq!(
        from_file.total() + edge_loci,
        scanned_total,
        "from the file: {:?}\nfrom a scan: {scanned_counts:?}\nat contig edges: {edge_loci}",
        from_file.iter_sorted()
    );

    // And stratum by stratum, away from the edges.
    for ((period, copies), count) in from_file.iter_sorted() {
        let scanned_here = scanned_counts.get(&(period, copies)).copied().unwrap_or(0);
        assert!(
            count <= scanned_here,
            "period {period}, {copies} copies: the file holds {count}, a scan finds {scanned_here}"
        );
    }
}

/// A sample drawn from the file is stable under its seed and never larger than the cap, and
/// the counts beside it see every locus (spec §5.3).
#[test]
fn the_sample_is_capped_stable_and_counted_in_full() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let (catalog_path, reference, _) = build_catalog(dir.path(), &contigs);
    let catalog =
        RepeatCatalog::open_checking_against_reference(&catalog_path, &reference).expect("opens");
    let criteria = StrRepeatCriteria::default();

    let (counts, sample) = catalog
        .sample_loci_per_stratum(&criteria, ReadScope::WholeReference, 1, 42)
        .expect("servable");
    let (counts_again, sample_again) = catalog
        .sample_loci_per_stratum(&criteria, ReadScope::WholeReference, 1, 42)
        .expect("servable");

    assert_eq!(counts, counts_again);
    assert_eq!(
        sample.iter_sorted(),
        sample_again.iter_sorted(),
        "the same seed keeps the identical loci"
    );
    for ((period, copies), loci) in sample.iter_sorted() {
        assert!(
            loci.len() <= 1,
            "period {period}, {copies} copies: over cap"
        );
        assert!(counts.get(period, copies) >= loci.len() as u64);
    }
    assert_eq!(
        counts.total(),
        catalog
            .str_loci(&criteria, ReadScope::WholeReference)
            .expect("servable")
            .count() as u64,
        "the counts see every locus, sampled or not"
    );
}

/// **Wiring step 3 to the catalog needs a region subset**, since every consumer of the walk
/// asks for one. What the walk emits for a set of spans and what the catalog emits for the
/// same spans must be the same regions — including the rule that a locus overlapping an edge
/// comes out whole while a generic stretch is clipped.
#[test]
fn a_region_subset_from_the_file_equals_the_walk_over_the_same_spans() {
    use pop_var_caller::ng::WindowedRefSeq;
    use pop_var_caller::ng::region_typing::{GenomeRegions, TypedRegionIterator};
    use pop_var_caller::regions::ContigBounds;

    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let fasta = write_fasta(dir.path(), &contigs);
    let catalog_path = dir.path().join("ref.fa.repeats.parquet");
    let mut builder = RepeatCatalogBuilder::create(
        &catalog_path,
        StrRepeatCriteria::default(),
        ScanParams::default(),
    )
    .expect("builder");
    let reference = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        },
        &mut builder,
    )
    .expect("the pass runs");
    builder.finish(&reference).expect("finish");

    let criteria = StrRepeatCriteria::default();
    let bounds: Vec<ContigBounds> = reference
        .contigs
        .iter()
        .map(|c| ContigBounds {
            name: &c.name,
            length: c.length as u32,
        })
        .collect();

    // Two spans on one contig, chosen to cut through the middle of the sequence rather than
    // to land tidily between features. A BED on disk because that is the only way to build a
    // region set — `GenomeRegions` has no in-memory constructor.
    let clean_len = contigs[0].1.len() as u64;
    let bed = dir.path().join("wanted.bed");
    std::fs::write(
        &bed,
        format!(
            "{name}\t149\t450\n{name}\t599\t{end}\n",
            name = contigs[0].0,
            end = clean_len - 50
        ),
    )
    .expect("write bed");
    let spans = GenomeRegions::from_bed_path(&bed, &bounds).expect("valid spans");

    let walk_config = TypedRegionConfig {
        criteria: criteria.classification.clone(),
        max_str_len: criteria.max_str_len_bp,
        ..TypedRegionConfig::default()
    };
    // The walk seeks by index, so it needs a `.fai`; the catalog path needs none.
    pop_var_caller::ng::reference_info::write_fai(
        &reference.contigs,
        &pop_var_caller::ng::reference_info::sibling_fai_path(&fasta),
    )
    .expect("write fai");
    let walked: Vec<TypedRegion> = TypedRegionIterator::over_regions(
        WindowedRefSeq::new(fasta.clone(), reference.contig_list()),
        spans.clone(),
        walk_config,
    )
    .expect("the walk starts")
    .map(|r| r.expect("a region"))
    .collect();

    let catalog =
        RepeatCatalog::open_checking_against_reference(&catalog_path, &reference).expect("opens");
    let wanted: Vec<GenomeRegion> = spans.iter().collect();
    let from_file: Vec<TypedRegion> = catalog
        .genome_segments(&criteria, ReadScope::Regions(&wanted))
        .expect("servable")
        .map(|r| r.expect("a region"))
        .collect();

    assert!(!walked.is_empty(), "the spans must contain something");
    assert_eq!(
        from_file, walked,
        "the catalog's region subset and the walk's disagree\nfile: {from_file:#?}\nwalk: {walked:#?}"
    );
}

/// **Classification is not local, and this is what proves the file's answer accounts for it.**
///
/// A read that wants a stretch must also read every row that can reach *into* that stretch: a
/// tract just outside it bundles with one inside, and a long array overlapping it swallows a
/// locus inside. The catalog widens each requested span by the contig's longest stored tract
/// plus the reader's bundle radius, which is what `longest_tract_bp` is in the header for.
///
/// Nothing held either half of that until now. Deleting the widening, or dropping the
/// longest-tract term from it, each left the whole suite green while turning a bundle into a
/// locus and a satellite into a locus — which is the difference between skipping an array and
/// calling it.
///
/// Both windows here start *after* the feature they must still see, so the widening is the
/// only thing that can bring it in.
#[test]
fn a_feature_outside_the_window_still_decides_what_is_inside_it() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let sides = both_sides(dir.path(), &contigs);
    let criteria = StrRepeatCriteria::default();
    let catalog =
        RepeatCatalog::open_checking_against_reference(&sides.catalog_path, &sides.reference)
            .expect("opens");

    let index_of = |name: &str| {
        sides
            .reference
            .contigs
            .iter()
            .position(|c| c.name == name)
            .expect("the fixture has this contig") as u32
    };

    // `chr_bundled` carries two tracts 10 bases apart. A window over the second alone must
    // still bundle them, because the first is within the bundle radius of it.
    let inside_the_second_tract = [GenomeRegion {
        contig: ContigId(index_of("chr_bundled")),
        start: Position(245),
        end: Position(260),
    }];
    let kinds: Vec<&'static str> = catalog
        .genome_segments(&criteria, ReadScope::Regions(&inside_the_second_tract))
        .expect("servable")
        .map(|r| match r.expect("a region").kind {
            RegionKind::SsrBundle { .. } => "bundle",
            RegionKind::SsrSegment(_) => "locus",
            RegionKind::Satellite => "satellite",
            RegionKind::Generic => "generic",
        })
        .collect();
    assert!(
        kinds.contains(&"bundle"),
        "a window over one member of a bundle must still see the other: {kinds:?}"
    );
    assert!(
        !kinds.contains(&"locus"),
        "the bundled tracts must not come back as clean loci: {kinds:?}"
    );

    // `chr_array` carries a 1.2 kb array running from base 201. A window over its last few
    // bases must call the whole thing a satellite and hand it back **whole**: a satellite's
    // extent is its claim, and 60 bases of "array too long to be a microsatellite"
    // contradicts the cap that produced the label.
    let inside_the_array = [GenomeRegion {
        contig: ContigId(index_of("chr_array")),
        start: Position(1_390),
        end: Position(1_450),
    }];
    let regions: Vec<TypedRegion> = catalog
        .genome_segments(&criteria, ReadScope::Regions(&inside_the_array))
        .expect("servable")
        .map(|r| r.expect("a region"))
        .collect();
    let satellite = regions
        .iter()
        .find(|r| matches!(r.kind, RegionKind::Satellite))
        .unwrap_or_else(|| {
            panic!("a window inside a 1.2 kb array must be satellite: {regions:#?}")
        });

    // **And it comes out whole, not clipped to the window.** A satellite's extent is its
    // claim — a 50-base "array too long to be a microsatellite" contradicts the cap that
    // produced the label — so only generic stretches clip at a requested edge, which is the
    // walk's own rule.
    assert!(
        satellite.region.start.get() < 1_390,
        "the satellite must be emitted whole, not cut to the window: {:?}",
        satellite.region
    );
}

/// Sequence with no structure of its own, from a fixed seed — the same bases every run.
fn pseudorandom(len: usize, seed: u64) -> String {
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[((state >> 33) % 4) as usize] as char
        })
        .collect()
}

/// [`fixture`] plus 200 kb of sequence with no designed structure.
///
/// **The hand-built contigs drive no rejection at all**, which would make the tally
/// comparison assert four zeroes against four zeroes. Two of admission's gates need a
/// haystack rather than a planted tract: over 200 kb of random bases the detector emits a few
/// dozen tracts that lose a copy to the whole-motif cut (the copy floor) and a few dozen with
/// no whole-motif boundaries to cut back to at all. The other two gates — an impure tract and
/// a compound motif — a scan does not reach even over 2 Mb, and they are pinned instead by
/// `segments.rs`'s own tests, which hand admission the rows directly.
fn tally_fixture() -> Vec<(&'static str, String)> {
    let mut contigs = fixture();
    contigs.push(("chr_random", pseudorandom(200_000, 5)));
    contigs
}

/// Everything a consumer needs to run both sides of the tally comparison over one contig
/// subset: the FASTA (which the walk reads), the catalog, and the reference the pass
/// computed.
struct BothSides {
    fasta: PathBuf,
    catalog_path: PathBuf,
    reference: ReferenceInfo,
}

fn both_sides(dir: &Path, contigs: &[(&str, String)]) -> BothSides {
    let fasta = write_fasta(dir, contigs);
    let catalog_path = dir.join("ref.fa.repeats.parquet");
    let mut builder = RepeatCatalogBuilder::create(
        &catalog_path,
        StrRepeatCriteria::default(),
        ScanParams::default(),
    )
    .expect("builder");
    let reference = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        },
        &mut builder,
    )
    .expect("the reference pass runs");
    builder.finish(&reference).expect("finish");
    // The walk seeks by index, so it needs a `.fai`; the catalog path opens no reference.
    pop_var_caller::ng::reference_info::write_fai(
        &reference.contigs,
        &pop_var_caller::ng::reference_info::sibling_fai_path(&fasta),
    )
    .expect("write fai");
    BothSides {
        fasta,
        catalog_path,
        reference,
    }
}

/// The named contigs, whole, as the region list both sides are asked with.
fn whole_contigs_named(reference: &ReferenceInfo, names: &[&str]) -> Vec<GenomeRegion> {
    reference
        .contigs
        .iter()
        .enumerate()
        .filter(|(_, info)| names.contains(&info.name.as_str()))
        .map(|(index, info)| GenomeRegion {
            contig: ContigId(index as u32),
            start: Position(1),
            end: Position(info.length),
        })
        .collect()
}

/// What the walk tallies over `spans`, and what the catalog tallies over the same spans.
fn both_tallies(
    sides: &BothSides,
    spans: &[GenomeRegion],
    criteria: &StrRepeatCriteria,
) -> (
    pop_var_caller::ng::region_typing::TypedRegionCounts,
    pop_var_caller::ng::repeat_catalog::CatalogRegionCounts,
    Vec<TypedRegion>,
) {
    use pop_var_caller::ng::WindowedRefSeq;
    use pop_var_caller::ng::region_typing::{GenomeRegions, TypedRegionIterator};
    use pop_var_caller::regions::ContigBounds;

    let bounds: Vec<ContigBounds> = sides
        .reference
        .contigs
        .iter()
        .map(|c| ContigBounds {
            name: &c.name,
            length: c.length as u32,
        })
        .collect();
    // The walk takes its spans through `GenomeRegions`, which is built from whole contigs or
    // a BED and nothing else — so the subset goes to disk as a BED.
    let bed = sides.fasta.with_extension("wanted.bed");
    let lines: String = spans
        .iter()
        .map(|s| {
            format!(
                "{}\t{}\t{}\n",
                bounds[s.contig.get() as usize].name,
                s.start.get() - 1,
                s.end.get()
            )
        })
        .collect();
    std::fs::write(&bed, lines).expect("write bed");
    let regions = GenomeRegions::from_bed_path(&bed, &bounds).expect("valid spans");

    let walk_config = TypedRegionConfig {
        criteria: criteria.classification.clone(),
        max_str_len: criteria.max_str_len_bp,
        ..TypedRegionConfig::default()
    };
    let mut walk = TypedRegionIterator::over_regions(
        WindowedRefSeq::new(sides.fasta.clone(), sides.reference.contig_list()),
        regions,
        walk_config,
    )
    .expect("the walk starts");
    let walked: Vec<TypedRegion> = walk.by_ref().map(|r| r.expect("a region")).collect();

    let catalog =
        RepeatCatalog::open_checking_against_reference(&sides.catalog_path, &sides.reference)
            .expect("opens");
    let mut from_file = catalog
        .genome_segments(criteria, ReadScope::Regions(spans))
        .expect("servable");
    let derived: Vec<TypedRegion> = from_file.by_ref().map(|r| r.expect("a region")).collect();

    assert_eq!(
        derived, walked,
        "the regions themselves must agree before their tallies mean anything"
    );
    (*walk.counts(), *from_file.counts(), walked)
}

/// **The tally a consumer keeps when it stops walking the reference.** Over contigs whose
/// repeats all sit clear of the ends, every counter the catalog has holds the walk's number.
///
/// The five contigs excluded are the two carrying structure within 15 bases of a contig end;
/// `the_tally_differs_only_at_the_contig_ends` is where those are accounted for, and it says
/// by how much.
#[test]
fn the_tally_from_the_file_matches_the_walks() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = tally_fixture();
    let sides = both_sides(dir.path(), &contigs);
    let criteria = StrRepeatCriteria::default();

    let spans = whole_contigs_named(
        &sides.reference,
        &[
            "chr_clean",
            "chr_bundled",
            "chr_array",
            "chr_mid_array",
            "chr_impure",
            "chr_random",
        ],
    );
    let (walk, file, walked) = both_tallies(&sides, &spans, &criteria);

    // The fixture has to exercise what the comparison claims to compare.
    assert!(
        walk.ssr_loci > 0 && walk.ssr_bundles > 0 && walk.satellites > 0 && walk.generic > 0,
        "these contigs must produce all four kinds of region, or the comparison proves \
         little: {walk:?}"
    );
    assert!(
        walk.rejected_by_reason.copy_floor > 0 && walk.rejected_by_reason.no_clean_trim > 0,
        "two of admission's gates must actually fire, or comparing them compares zero \
         against zero: {:?}",
        walk.rejected_by_reason
    );
    assert_eq!(
        walk.rejected_by_reason.flank_clamped, 0,
        "these contigs were chosen for having nothing at a contig end; if that stops being \
         true the exact comparison below is no longer the right test"
    );

    assert_eq!(walk.spans, file.spans);
    assert_eq!(walk.ssr_loci, file.ssr_loci);
    assert_eq!(walk.ssr_bundles, file.ssr_bundles);
    assert_eq!(walk.ssr_bundle_bp, file.ssr_bundle_bp);
    assert_eq!(walk.generic, file.generic);
    assert_eq!(walk.satellites, file.satellites);
    assert_eq!(walk.satellite_bp, file.satellite_bp);
    assert_eq!(
        walk.repeat_bp_with_no_locus, file.repeat_bp_with_no_locus,
        "repeat coverage that yielded no locus"
    );

    let reasons = walk.rejected_by_reason;
    assert_eq!(reasons.copy_floor, file.rejected_by_reason.copy_floor);
    assert_eq!(reasons.purity, file.rejected_by_reason.purity);
    assert_eq!(reasons.compound, file.rejected_by_reason.compound);
    assert_eq!(reasons.no_clean_trim, file.rejected_by_reason.no_clean_trim);

    // And the regions counted are the regions emitted, on both sides.
    assert_eq!(
        walk.generic + walk.ssr_loci + walk.ssr_bundles + walk.satellites,
        walked.len() as u64
    );
}

/// **The two differences, named and measured** — the file is not short by accident, and this
/// is where it says by how much.
///
/// `chr_edges` carries a 30-base tract hard against each contig end: a live scan sees both,
/// finds no sequence to anchor them against, and charges them to the contig-end rejection.
/// `chr_near_edge` carries one 5 bases from the end: a live scan makes it a **locus**. The
/// catalog holds neither, because the file stores nothing with less than 15 bases beside it.
#[test]
fn the_tally_differs_only_at_the_contig_ends() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let sides = both_sides(dir.path(), &contigs);

    let criteria = StrRepeatCriteria::default();

    let spans = whole_contigs_named(&sides.reference, &["chr_edges", "chr_near_edge"]);

    // The regions themselves differ here, so `both_tallies`' region check would fire — run
    // the two sides directly instead.
    use pop_var_caller::ng::WindowedRefSeq;
    use pop_var_caller::ng::region_typing::{GenomeRegions, TypedRegionIterator};
    use pop_var_caller::regions::ContigBounds;

    let bounds: Vec<ContigBounds> = sides
        .reference
        .contigs
        .iter()
        .map(|c| ContigBounds {
            name: &c.name,
            length: c.length as u32,
        })
        .collect();
    let bed = dir.path().join("edges.bed");
    let lines: String = spans
        .iter()
        .map(|s| {
            format!(
                "{}\t0\t{}\n",
                bounds[s.contig.get() as usize].name,
                s.end.get()
            )
        })
        .collect();
    std::fs::write(&bed, lines).expect("bed");
    let regions = GenomeRegions::from_bed_path(&bed, &bounds).expect("spans");
    let walk_config = TypedRegionConfig {
        criteria: criteria.classification.clone(),
        max_str_len: criteria.max_str_len_bp,
        ..TypedRegionConfig::default()
    };
    let mut walk = TypedRegionIterator::over_regions(
        WindowedRefSeq::new(sides.fasta.clone(), sides.reference.contig_list()),
        regions,
        walk_config,
    )
    .expect("walk");
    for region in walk.by_ref() {
        region.expect("a region");
    }
    let walk_counts = *walk.counts();

    let catalog =
        RepeatCatalog::open_checking_against_reference(&sides.catalog_path, &sides.reference)
            .expect("opens");
    let mut segments = catalog
        .genome_segments(&criteria, ReadScope::Regions(&spans))
        .expect("servable");
    for region in segments.by_ref() {
        region.expect("a region");
    }
    let file_counts = *segments.counts();

    // 1. Two 30-base tracts abut a contig end. The walk charges 60 bases to the contig-end
    //    rejection; the catalog has no counter for them at all, which is the point — a `0`
    //    here would say "this genome has none", and it has two.
    assert_eq!(
        walk_counts.rejected_by_reason.flank_clamped, 60,
        "chr_edges' two 30-base tracts, one at each end"
    );

    // 2. Those same 60 bases are repeat coverage that yielded no locus on the walk's side and
    //    coverage the file never saw on the catalog's, so the totals differ by exactly them.
    assert_eq!(
        walk_counts.repeat_bp_with_no_locus - file_counts.repeat_bp_with_no_locus,
        60,
        "walk {} vs file {}",
        walk_counts.repeat_bp_with_no_locus,
        file_counts.repeat_bp_with_no_locus
    );

    // 3. `chr_near_edge`'s tract is a locus to the walk and absent from the file — the one
    //    place the file holds fewer loci, and it is 5 bases from a contig end.
    assert_eq!(
        walk_counts.ssr_loci - file_counts.ssr_loci,
        1,
        "walk {} loci vs file {}",
        walk_counts.ssr_loci,
        file_counts.ssr_loci
    );
    // And the locus splits a generic stretch in two on the walk's side — generic, locus,
    // generic — where the file has one uninterrupted stretch.
    assert_eq!(walk_counts.generic - file_counts.generic, 1);

    // Everything that is not at a contig end still agrees.
    assert_eq!(walk_counts.spans, file_counts.spans);
    assert_eq!(walk_counts.ssr_bundles, file_counts.ssr_bundles);
    assert_eq!(walk_counts.satellites, file_counts.satellites);
    assert_eq!(
        walk_counts.rejected_by_reason.copy_floor,
        file_counts.rejected_by_reason.copy_floor
    );
    assert_eq!(
        walk_counts.rejected_by_reason.purity,
        file_counts.rejected_by_reason.purity
    );
    assert_eq!(
        walk_counts.rejected_by_reason.compound,
        file_counts.rejected_by_reason.compound
    );
    assert_eq!(
        walk_counts.rejected_by_reason.no_clean_trim,
        file_counts.rejected_by_reason.no_clean_trim
    );
}

/// A reader more permissive than the file is refused **before any row is read**, naming the
/// axis and both values (spec §4.3, §10.2).
#[test]
fn a_more_permissive_reader_is_refused_eagerly() {
    let dir = tempfile::tempdir().expect("tmp");
    let contigs = fixture();
    let (catalog_path, reference, _) = build_catalog(dir.path(), &contigs);
    let catalog =
        RepeatCatalog::open_checking_against_reference(&catalog_path, &reference).expect("opens");

    let base = StrRepeatCriteria::default();

    let lower_floor = StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            min_copies: MinCopies::new([2; 6], CATALOG_MIN_COPIES_BEYOND_TABLE),
            ..base.classification.clone()
        },
        ..base.clone()
    };
    let err = catalog
        .genome_segments(&lower_floor, ReadScope::WholeReference)
        .err()
        .expect("a lower copy floor is refused");
    let message = format!("{err}");
    assert!(message.contains("period 1"), "names the axis: {message}");
    assert!(message.contains('5'), "names the built floor: {message}");

    let smaller_flank = StrRepeatCriteria {
        min_flank_bp: Bp(5),
        ..base.clone()
    };
    assert!(
        catalog
            .genome_segments(&smaller_flank, ReadScope::WholeReference)
            .is_err()
    );

    // And the mirror case: the unbounded axes are served, not refused.
    let permissive_elsewhere = StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            min_purity: 0.0,
            bundle_threshold: 40,
            ..base.classification.clone()
        },
        max_str_len_bp: Bp(1_000_000),
        ..base
    };
    assert!(
        catalog
            .genome_segments(&permissive_elsewhere, ReadScope::WholeReference)
            .is_ok()
    );
}
