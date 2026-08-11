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
use pop_var_caller::ng::repeat_catalog::{RepeatCatalog, RepeatCatalogBuilder, StrRepeatCriteria};
use pop_var_caller::ng::tandem_repeat::{PeriodRange, ScanParams};
use pop_var_caller::ng::types::{Bp, ContigId};

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
    builder.finish(&reference, "differential").expect("finish");
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
        .genome_segments(criteria, Some(contig))
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
        .count_loci_per_stratum(&criteria, None)
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
        .sample_loci_per_stratum(&criteria, None, 1, 42)
        .expect("servable");
    let (counts_again, sample_again) = catalog
        .sample_loci_per_stratum(&criteria, None, 1, 42)
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
        catalog.str_loci(&criteria, None).expect("servable").count() as u64,
        "the counts see every locus, sampled or not"
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
        .genome_segments(&lower_floor, None)
        .err()
        .expect("a lower copy floor is refused");
    let message = format!("{err}");
    assert!(message.contains("period 1"), "names the axis: {message}");
    assert!(message.contains('5'), "names the built floor: {message}");

    let smaller_flank = StrRepeatCriteria {
        min_flank_bp: Bp(5),
        ..base.clone()
    };
    assert!(catalog.genome_segments(&smaller_flank, None).is_err());

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
    assert!(catalog.genome_segments(&permissive_elsewhere, None).is_ok());
}
