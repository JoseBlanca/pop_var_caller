//! Is a duplicated locus a class worth carrying, and can window coverage find one?
//!
//! `parameter_prepass_joint_fit.md` §2.2 adopts a third class of site — a locus the
//! *sample* carries more copies of than the reference does, collecting two copies' reads at
//! one position and so reading about half alternative wherever the copies differ — and
//! makes its discriminator **local relative coverage**, never the site's own depth. The
//! object that costs something is
//! [`parameter_prepass_joint_records.md`](../doc/devel/ng/spec/parameter_prepass_joint_records.md)
//! §4's per-sample coverage-by-window summary, 1.6 MB per sample on tomato. Both documents
//! gate it on one measurement, and this program is that measurement:
//!
//! > on one sample, what fraction of positions sit in windows near two copies **and** show a
//! > near-half alternative fraction.
//!
//! One walk over one alignment, no cohort and no truth set.
//!
//! # What it does
//!
//! One pass over the reference's generic regions inside the BED gives each fixed window its
//! **denominator** — how many generic positions it holds — and its **GC fraction**, over
//! exactly those positions. One walk over the alignment then gives every generic position
//! its observation depth and, at single-base loci, its alternative-read count.
//!
//! Coverage tracks GC content, so a window near an extreme of it reads high or low for a
//! reason that has nothing to do with copy number. The curve is fitted from this sample:
//! the median window depth in each GC bin, falling back to the global median where a bin is
//! thin. A window's **relative coverage** is its mean depth over the curve's value at its
//! own GC, rescaled so the median window sits at 1.0 — so 2.0 is a window carrying twice
//! the sample's normal depth.
//!
//! Both the raw and the GC-corrected classification are printed, because what the correction
//! is worth on this sample is itself unknown.
//!
//! # The three discriminators it compares
//!
//! The window summary is the expensive object: a genome-wide depth accumulator, a GC curve
//! fitted per sample, denominators taken from the reference, and a second full pass over
//! every sample's stored pileup when the estimate runs from stored data. The cheap
//! alternative is already in the records — **the position's own read count**, five bits per
//! kept position. So the probe scores the same positions three ways and prints the three
//! side by side:
//!
//! 1. **the window's relative coverage** — the sample's mean depth in the window that holds
//!    the position, over what the sample's own depth-against-GC curve predicts there;
//! 2. **the position's own relative depth** — its read count over what one copy would give
//!    at that position, with and without the same GC correction, and once more after the
//!    read count has been through the depth ladder the records actually store it in;
//! 3. **both at once** — a position only flagged when the window and the read count agree.
//!
//! A fourth print, the near-half rate in each cell of window-band × position-band, says
//! whether either discriminator adds anything to the other.
//!
//! # What it cannot say
//!
//! **There is no truth set here**, so nothing below establishes that a window near two
//! copies *is* a duplication. What it can establish is whether the population the class was
//! invented for exists at all, how big it is, and whether relative coverage separates it
//! from the rest — which is what decides whether the accumulator is worth building. Every
//! number the arms are compared on is *enrichment over independence*, a comparison between
//! discriminators and never an accuracy.
//!
//! Depth is *observation* depth, the same quantity the generic path's accumulator sees: a
//! read that covered a position but anchored no border is not counted
//! (`SampleLocusObservations::num_obs_along_locus`).
//!
//! # The per-window dump
//!
//! With a fifth argument the probe writes one row per window — its coordinates, its mean
//! depth, its GC, its relative coverage and how many of its positions read near half. That
//! is what answers a question the gate number cannot: **are the two-copy windows the same
//! windows in every sample?** If they are, the amplification is a property of the reference
//! and the aligner, a per-locus class would hold it, and
//! `parameter_prepass_joint_records.md` §4's per-sample object is not needed. If they differ
//! between samples, the grain in `parameter_prepass_joint_fit.md` §2.2 — the (locus, sample)
//! pair — is the right one.
//!
//! ```text
//! ng_joint_duplication_probe <reference.fa> <catalog.parquet> <alignment> <regions.bed> \
//!     [window-bp] [windows.tsv]
//! ```

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, LocusKind, SampleLocusObservationsIterator, UnhandledReason,
};
use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::{RefSeq, WindowedRefSeq};
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{
    GenomeRegions, RegionKind, TypedRegion, TypedRegionConfig,
};
use pop_var_caller::ng::repeat_catalog::{
    ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria,
};
use pop_var_caller::ng::types::GenomeRegion;
use pop_var_caller::regions::ContigBounds;

/// One fixed window of the reference, as this sample fills it.
#[derive(Default, Clone)]
struct Window {
    /// Which contig, and which window of it — kept so the dump can name the window.
    contig: u32,
    /// The window's first base, 1-based.
    start: u64,
    /// Generic positions the window holds — the denominator, read off the typed regions and
    /// not off the walk, so a window whose reads are missing still divides by its full size.
    generic_positions: u32,
    /// Positions the walk actually emitted a locus for.
    visited_positions: u32,
    /// Σ observation depth over the emitted positions.
    depth_sum: u64,
    /// Reference `G` or `C` over the generic positions.
    gc: u32,
    /// Reference `A` or `T` over the generic positions.
    at: u32,
}

impl Window {
    fn mean_depth(&self) -> f64 {
        if self.generic_positions == 0 {
            0.0
        } else {
            self.depth_sum as f64 / f64::from(self.generic_positions)
        }
    }

    fn gc_fraction(&self) -> Option<f64> {
        let called = self.gc + self.at;
        (called > 0).then(|| f64::from(self.gc) / f64::from(called))
    }
}

/// One generic position, as the walk saw it. Twelve bytes, because there are millions.
#[derive(Copy, Clone)]
struct PositionObservation {
    /// Which coverage window holds it — the index the window arm reads.
    window: u32,
    /// Which fixed 500 bp stretch of reference holds it. Separate from `window` because the
    /// coverage window's width is a command-line argument and the GC a position sits in is
    /// not: GC bias acts at about a fragment's length however wide the coverage window is.
    gc_window: u32,
    depth: u16,
    alt: u16,
}

/// How much of a window must be generic, and covered, before its mean depth means anything.
const MIN_GENERIC_SHARE: f64 = 0.5;

/// GC bins for the depth-against-GC curve: 2 percentage points wide.
const GC_BINS: usize = 51;

/// Below this many windows a GC bin's own median is scatter, and the global median is used.
const MIN_WINDOWS_PER_GC_BIN: usize = 50;

/// The stretch of reference a single position's GC content is read over.
///
/// **Fixed, and not the coverage window's width.** GC bias enters through the fragments that
/// were sequenced, so the GC a position's own depth should be corrected against is the GC of
/// a fragment-sized neighbourhood — the same 500 bp the summary stores at, and unchanged
/// when the coverage arm is run at 5 kb.
const GC_WINDOW_BP: u64 = 500;

/// Below this many positions a GC bin's own mean depth is scatter, and the global mean is
/// used. Positions, not windows: this curve is fitted over single positions.
const MIN_POSITIONS_PER_GC_BIN: u64 = 10_000;

/// The same floor for the sparse curve, which has a few hundred positions a bin to work
/// with rather than a few hundred thousand. Set at the point where a GC bin's mean depth is
/// worth about a percent at tomato's depth: 100 positions at three reads each has a standard
/// error of 17% of one read, which is under 6% of the 1.79-fold swing being corrected for.
const MIN_POSITIONS_PER_SPARSE_GC_BIN: u64 = 100;

/// The deepest position the depth histograms count exactly. Anything deeper is counted in
/// the overflow, which only the maximum and the above-cap share ever read.
const DEPTH_HISTOGRAM_LIMIT: usize = 4_096;

/// One walked position in this many feeds the sparse depth-against-GC curve — the order of
/// how many positions the joint records keep, so the curve it produces is the one a fit that
/// had only the records could have fitted.
const SPARSE_CURVE_STRIDE: u64 = 300;

/// Which entry of [`BANDS`] is "two copies". Named once because four reports index it.
const TWO_COPY_BAND: usize = 3;

fn main() {
    let mut args = std::env::args().skip(1);
    let usage = "usage: <reference.fa> <catalog.parquet> <alignment> <regions.bed> [window-bp]";
    let fasta = PathBuf::from(args.next().expect(usage));
    let catalog_path = PathBuf::from(args.next().expect(usage));
    let reads = PathBuf::from(args.next().expect(usage));
    let bed = PathBuf::from(args.next().expect(usage));
    let window_bp: u64 = args
        .next()
        .map_or(500, |a| a.parse().expect("a window size"));
    let dump = args.next().map(PathBuf::from);

    let started = Instant::now();

    // ---- the reference, the catalog, the regions -------------------------------------
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verification) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.clone(),
        ReferenceCheck::VerifyAgainstIndex,
    )
    .expect("the reference is readable and has (or can derive) a .fai");
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(&fasta).expect("the .fai beside the reference reads");
    let catalog = RepeatCatalog::open_checking_against_reference(&catalog_path, &info)
        .expect("the catalog is this reference's");

    let bounds: Vec<ContigBounds> = contigs
        .entries
        .iter()
        .map(|entry| ContigBounds {
            name: &entry.name,
            length: u32::try_from(entry.length).expect("a contig shorter than 4 Gb"),
        })
        .collect();
    let spans: Vec<GenomeRegion> = GenomeRegions::from_bed_path(&bed, &bounds)
        .expect("the BED resolves against this reference's contigs")
        .iter()
        .collect();
    let bed_bases: u64 = spans.iter().map(|s| s.len()).sum();

    let criteria = StrRepeatCriteria::from(&TypedRegionConfig::default());
    let typed: Vec<TypedRegion> = catalog
        .genome_segments(&criteria, ReadScope::Regions(&spans))
        .expect("the BED's spans name contigs this catalog holds")
        .map(|item| item.expect("the catalog reads through the whole of the BED"))
        .collect();
    let generic: Vec<GenomeRegion> = typed
        .iter()
        .filter(|region| region.kind == RegionKind::Generic)
        .map(|region| region.region)
        .collect();
    let generic_bases: u64 = generic.iter().map(|r| r.len()).sum();

    println!("reference        {}", fasta.display());
    println!("catalog          {}", catalog_path.display());
    println!("alignment        {}", reads.display());
    println!(
        "regions          {} — {} spans, {} bases",
        bed.display(),
        spans.len(),
        bed_bases
    );
    println!(
        "generic regions  {} spans, {} bases ({:.1}% of the BED)",
        generic.len(),
        generic_bases,
        100.0 * generic_bases as f64 / bed_bases as f64
    );
    println!("window           {window_bp} bp");

    // ---- pass 1: every window's denominator and its GC ---------------------------------
    let accessor = WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone());
    let mut window_index: HashMap<u64, u32> = HashMap::new();
    let mut windows: Vec<Window> = Vec::new();
    // The fixed 500 bp stretches a single position's GC is read over — `(G+C, A+T)` each,
    // and nothing else, because GC comes off the reference and never off a depth
    // accumulator.
    let mut gc_window_index: HashMap<u64, u32> = HashMap::new();
    let mut gc_windows: Vec<(u32, u32)> = Vec::new();
    let mut bases = Vec::new();
    for region in &generic {
        accessor
            .fetch_into(region.contig, region.start.get(), region.len(), &mut bases)
            .expect("a generic region reads from the reference");
        for (offset, base) in bases.iter().enumerate() {
            let position = region.start.get() + offset as u64;
            let ordinal = position / window_bp;
            let key = window_key(region.contig.0, ordinal);
            let slot = *window_index.entry(key).or_insert_with(|| {
                windows.push(Window {
                    contig: region.contig.0,
                    start: ordinal * window_bp,
                    ..Window::default()
                });
                (windows.len() - 1) as u32
            });
            let window = &mut windows[slot as usize];
            window.generic_positions += 1;
            let gc_key = window_key(region.contig.0, position / GC_WINDOW_BP);
            let gc_slot = *gc_window_index.entry(gc_key).or_insert_with(|| {
                gc_windows.push((0, 0));
                (gc_windows.len() - 1) as u32
            });
            match base.to_ascii_uppercase() {
                b'G' | b'C' => {
                    window.gc += 1;
                    gc_windows[gc_slot as usize].0 += 1;
                }
                b'A' | b'T' => {
                    window.at += 1;
                    gc_windows[gc_slot as usize].1 += 1;
                }
                _ => {}
            }
        }
    }
    println!(
        "windows          {} touched, {:.1} s to here",
        windows.len(),
        started.elapsed().as_secs_f64()
    );

    // ---- pass 2: the walk ---------------------------------------------------------------
    let read_groups = build_read_groups(std::slice::from_ref(&reads))
        .expect("the header declares its read groups");
    let sample = match read_groups.read_groups_per_sample() {
        [only] => only.clone(),
        other => panic!(
            "{} holds {} samples; this probe is per sample",
            reads.display(),
            other.len()
        ),
    };
    let reference = OpenReference::new(info.clone());
    let sample_reads = SampleReads::open(
        &sample,
        &read_groups,
        &reference,
        ReadFilterConfig::default(),
        true,
    )
    .expect("the alignment file opens against this reference");
    let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
        fasta.clone(),
        contigs.clone(),
        index.clone(),
    ));
    let generator = {
        let fasta = fasta.clone();
        let contigs = contigs.clone();
        let index = index.clone();
        #[allow(
            clippy::arc_with_non_send_sync,
            reason = "file-backed and single-threaded, as in real_alignments.rs"
        )]
        let reference = Arc::new(WindowedRefSeq::with_shared_index(
            fasta.clone(),
            contigs.clone(),
            index.clone(),
        ));
        PileupGenerator::new(
            reference,
            move || {
                WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone())
            },
            preparer,
            PileupGeneratorConfig::default(),
        )
        .expect("the generic generator builds against this reference")
    };
    let generators = GeneratorSet::new(
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        GeneratorSlot::Generator(Box::new(generator)),
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
    );
    let regions: Vec<Result<TypedRegion, RepeatCatalogError>> = typed.into_iter().map(Ok).collect();
    let mut stream =
        SampleLocusObservationsIterator::new(regions.into_iter(), sample_reads, generators);

    let mut observations: Vec<PositionObservation> = Vec::new();
    let mut loci = 0_u64;
    let mut wide_loci = 0_u64;
    let mut positions_outside_windows = 0_u64;
    let mut total_depth = 0_u64;
    // Every walked position's depth, once overall and once per GC bin: the scale the
    // position's own read count is read against, and the curve it is corrected by.
    let mut depth_histogram = vec![0_u64; DEPTH_HISTOGRAM_LIMIT + 1];
    let mut depth_by_gc = vec![vec![0_u32; DEPTH_HISTOGRAM_LIMIT + 1]; GC_BINS];
    // The same curve fitted from one walked position in `SPARSE_CURVE_STRIDE` — the order of
    // how many positions the records keep. If this curve gives the cheap arm what the dense
    // one does, correcting a position's depth for GC needs no genome-wide pass.
    let mut depth_by_gc_sparse = vec![vec![0_u32; DEPTH_HISTOGRAM_LIMIT + 1]; GC_BINS];
    let mut walked = 0_u64;
    let mut deepest_position = 0_u32;
    for locus in &mut stream {
        let locus = locus.expect("the walk runs to completion on a well-formed alignment");
        if locus.kind != LocusKind::Generic {
            continue;
        }
        loci += 1;
        let depths = locus.num_obs_along_locus();
        let single_base = locus.region.len() == 1;
        if !single_base {
            wide_loci += 1;
        }
        // The alternative count is only meaningful where the locus is one position: at a
        // wider one an observation's bases span the whole region and "the alternative at
        // this position" is not a question the record answers. Their depth still counts
        // towards the window, which is what the window is for.
        let alt = if single_base {
            locus
                .complete_observations()
                .filter(|obs| *obs.bases != *locus.reference_bases)
                .map(|obs| obs.num_obs)
                .sum::<u32>()
        } else {
            0
        };
        for (offset, depth) in depths.iter().enumerate() {
            let position = locus.region.start.get() + offset as u64;
            let key = window_key(locus.region.contig.0, position / window_bp);
            let Some(&slot) = window_index.get(&key) else {
                positions_outside_windows += 1;
                continue;
            };
            let gc_key = window_key(locus.region.contig.0, position / GC_WINDOW_BP);
            let gc_slot = *gc_window_index
                .get(&gc_key)
                .expect("a walked generic position sits in a GC stretch pass 1 filled");
            let window = &mut windows[slot as usize];
            window.visited_positions += 1;
            window.depth_sum += u64::from(*depth);
            total_depth += u64::from(*depth);

            deepest_position = deepest_position.max(*depth);
            let cell = (*depth as usize).min(DEPTH_HISTOGRAM_LIMIT);
            depth_histogram[cell] += 1;
            if let Some(fraction) = gc_fraction_of(gc_windows[gc_slot as usize]) {
                depth_by_gc[gc_bin(fraction)][cell] += 1;
                if walked.is_multiple_of(SPARSE_CURVE_STRIDE) {
                    depth_by_gc_sparse[gc_bin(fraction)][cell] += 1;
                }
            }
            walked += 1;

            if single_base {
                observations.push(PositionObservation {
                    window: slot,
                    gc_window: gc_slot,
                    depth: u16::try_from(*depth).unwrap_or(u16::MAX),
                    alt: u16::try_from(alt).unwrap_or(u16::MAX),
                });
            }
        }
    }
    if let Some(handle) = verification {
        handle
            .join()
            .expect("the .fai beside the reference describes it");
    }

    let visited: u64 = windows.iter().map(|w| u64::from(w.visited_positions)).sum();
    println!(
        "\nwalk             {loci} generic loci, {visited} positions, {:.1} s",
        started.elapsed().as_secs_f64()
    );
    println!(
        "  of them        {} loci wider than one base (their positions carry depth but no \
         alternative fraction)",
        wide_loci
    );
    println!(
        "  never emitted  {} of {} generic positions ({:.2}%) — no locus reached them",
        generic_bases - visited,
        generic_bases,
        100.0 * (generic_bases - visited) as f64 / generic_bases as f64
    );
    if positions_outside_windows > 0 {
        println!("  outside window {positions_outside_windows} (a bug if non-zero)");
    }
    println!(
        "  mean depth     {:.2} reads a position over the generic bases",
        total_depth as f64 / generic_bases as f64
    );
    println!(
        "  single-base    {} positions carry an alternative fraction",
        observations.len()
    );

    // ---- the windows, and what GC does to them ------------------------------------------
    let usable: Vec<u32> = (0..windows.len() as u32)
        .filter(|&i| {
            let w = &windows[i as usize];
            f64::from(w.generic_positions) >= MIN_GENERIC_SHARE * window_bp as f64
                && w.gc_fraction().is_some()
        })
        .collect();
    println!(
        "\nwindows          {} of {} are at least {:.0}% generic and callable",
        usable.len(),
        windows.len(),
        100.0 * MIN_GENERIC_SHARE
    );

    let mean_depths: Vec<f64> = usable
        .iter()
        .map(|&i| windows[i as usize].mean_depth())
        .collect();
    let global_median = median(&mut mean_depths.clone());
    println!("  window depth   {}", percentile_line(&mean_depths));
    println!("  median         {global_median:.3} reads a position");

    // The depth-against-GC curve, this sample's own.
    let mut per_bin: Vec<Vec<f64>> = vec![Vec::new(); GC_BINS];
    for &i in &usable {
        let w = &windows[i as usize];
        per_bin[gc_bin(w.gc_fraction().expect("a callable window"))].push(w.mean_depth());
    }
    let curve: Vec<f64> = per_bin
        .iter()
        .map(|depths| {
            if depths.len() >= MIN_WINDOWS_PER_GC_BIN {
                median(&mut depths.clone())
            } else {
                global_median
            }
        })
        .collect();
    println!("\n  depth against GC — the curve the correction divides by");
    println!("  GC%    windows     median depth   ratio to global");
    for (bin, depths) in per_bin.iter().enumerate() {
        if depths.len() < MIN_WINDOWS_PER_GC_BIN {
            continue;
        }
        println!(
            "  {:>3}  {:>9}  {:>13.3}  {:>15.3}",
            bin * 2,
            depths.len(),
            curve[bin],
            curve[bin] / global_median
        );
    }

    // Relative coverage, raw and corrected, each rescaled so the median window sits at 1.
    let raw: Vec<f64> = usable
        .iter()
        .map(|&i| windows[i as usize].mean_depth() / global_median)
        .collect();
    let corrected: Vec<f64> = usable
        .iter()
        .map(|&i| {
            let w = &windows[i as usize];
            let expected = curve[gc_bin(w.gc_fraction().expect("a callable window"))];
            if expected > 0.0 {
                w.mean_depth() / expected
            } else {
                0.0
            }
        })
        .collect();
    let raw = rescale_to_unit_median(raw);
    let corrected = rescale_to_unit_median(corrected);

    println!("\n  relative coverage per window, median rescaled to 1.0");
    println!("  raw            {}", percentile_line(&raw));
    println!("  GC-corrected   {}", percentile_line(&corrected));

    // Per-window relative coverage, indexed the way the observations are.
    let mut relative = vec![f64::NAN; windows.len()];
    for (slot, &i) in usable.iter().enumerate() {
        relative[i as usize] = corrected[slot];
    }
    let mut relative_raw = vec![f64::NAN; windows.len()];
    for (slot, &i) in usable.iter().enumerate() {
        relative_raw[i as usize] = raw[slot];
    }

    println!("\n  windows by band, and the positions in them");
    print_band_table(&windows, &relative, "GC-corrected");
    print_band_table(&windows, &relative_raw, "raw");

    // ---- the joint question --------------------------------------------------------------
    for &min_depth in &[2_u16, 4, 6, 8, 12, 16] {
        report_joint(&observations, &relative, min_depth, "GC-corrected");
    }
    report_joint(&observations, &relative_raw, 8, "raw");

    // ---- the cheap arm: the position's own read count -------------------------------------
    //
    // Everything below asks whether the five bits already in every kept record do what the
    // window summary does. The read count is read three ways — as the record stores it
    // (through the depth ladder, after the per-position cap), as the walk saw it, and without
    // the GC correction — so that what the encoding costs and what the correction is worth
    // are separable from what the arm is worth.
    let edges = DepthBinEdges::for_census();
    let visited_positions: u64 = depth_histogram.iter().sum();
    let mean_position_depth = total_depth as f64 / visited_positions as f64;
    let median_position_depth = median_of_histogram(&depth_histogram);

    println!("\n\n================ the position's own read count ================");
    println!(
        "\nwalked positions {visited_positions}, mean depth {mean_position_depth:.2}, \
         median {median_position_depth}, deepest {deepest_position}"
    );

    // Where the cap ends the arm. The record holds a bin, so a doubled position is only
    // visible while the bin its doubled count lands in is far enough above the bin one
    // copy's count lands in — and above the cap both are written as the top bin.
    let above_cap: u64 = depth_histogram[(edges.max_depth() as usize + 1)..]
        .iter()
        .sum();
    let (dulls_at, blind_at) = blinding_depths(&edges);
    println!(
        "  ladder cap {} reads a position; {above_cap} walked positions are over it \
         ({:.4}% of them)",
        edges.max_depth(),
        100.0 * above_cap as f64 / visited_positions as f64
    );
    println!(
        "  a doubled position's recorded count stops reaching 1.6 times an ordinary one's \
         from {dulls_at} reads a position,"
    );
    println!(
        "  and is written as the very same code from {blind_at} reads a position — so above \
         {blind_at} the arm cannot see a duplication at all."
    );
    println!(
        "  this sample sits at {median_position_depth} reads a position, which is {}",
        if median_position_depth >= u64::from(dulls_at) {
            "past the point where the cap limits the arm"
        } else {
            "below it, so the cap does not fire on this sample"
        }
    );

    // What one copy would give at this position's GC, fitted from the walk's own depths. The
    // mean and not the median, because a single position's depth is a small integer and its
    // median is a staircase — 2, 2, 3, 3 across the whole GC range at tomato's depth.
    let (exact_curve, _) = gc_curve(&depth_by_gc, MIN_POSITIONS_PER_GC_BIN, f64::from);
    let (recorded_curve, recorded_mean) =
        gc_curve(&depth_by_gc, MIN_POSITIONS_PER_GC_BIN, |depth| {
            recorded_depth(&edges, depth)
        });
    let (sparse_curve, _) = gc_curve(
        &depth_by_gc_sparse,
        MIN_POSITIONS_PER_SPARSE_GC_BIN,
        |depth| recorded_depth(&edges, depth),
    );
    let sparse_positions: u64 = depth_by_gc_sparse
        .iter()
        .flat_map(|counts| counts.iter().map(|&count| u64::from(count)))
        .sum();
    let flat_curve = vec![recorded_mean; GC_BINS];
    println!("\n  depth against GC over single positions — the cheap arm's own curve");
    println!("  GC%    positions    mean depth   ratio to global");
    for (bin, counts) in depth_by_gc.iter().enumerate() {
        let positions: u64 = counts.iter().map(|&count| u64::from(count)).sum();
        if positions < MIN_POSITIONS_PER_GC_BIN {
            continue;
        }
        println!(
            "  {:>3}  {:>10}  {:>12.3}  {:>15.3}",
            bin * 2,
            positions,
            exact_curve[bin],
            exact_curve[bin] / mean_position_depth
        );
    }

    let recorded = position_relative(&observations, &gc_windows, &recorded_curve, |depth| {
        recorded_depth(&edges, depth)
    });
    let exact = position_relative(&observations, &gc_windows, &exact_curve, f64::from);
    let sparse = position_relative(&observations, &gc_windows, &sparse_curve, |depth| {
        recorded_depth(&edges, depth)
    });
    let uncorrected = position_relative(&observations, &gc_windows, &flat_curve, |depth| {
        recorded_depth(&edges, depth)
    });
    let window_here: Vec<f32> = observations
        .iter()
        .map(|observation| relative[observation.window as usize] as f32)
        .collect();

    println!("\n  relative depth per position, median rescaled to 1.0");
    println!(
        "  recorded       {}",
        percentile_line(&sample_every(&recorded, 37))
    );
    println!(
        "  exact          {}",
        percentile_line(&sample_every(&exact, 37))
    );

    // Above the cap a position's alternative reads are thinned by the same ratio as its
    // depth, so the fraction should survive it. This counts the positions where rounding
    // moves that fraction across a band edge anyway.
    let thinning_flips = observations
        .iter()
        .filter(|observation| {
            let depth = u32::from(observation.depth);
            let alt = thin_to_cap(u32::from(observation.alt), depth, edges.max_depth());
            near_half_of(alt, depth.min(edges.max_depth())) != near_half(observation)
        })
        .count();
    println!(
        "  positions whose near-half status changes when the cap thins their reads: \
         {thinning_flips}"
    );

    // ---- the arms, side by side -----------------------------------------------------------
    for &min_depth in &[4_u16, 8] {
        report_arms(
            &observations,
            &window_here,
            min_depth,
            &[
                ("window relative coverage", &window_here),
                ("position, recorded (binned, capped)", &recorded),
                ("position, exact read count", &exact),
                ("position, recorded, no GC correction", &uncorrected),
                ("position, recorded, GC curve from 1 in 300", &sparse),
            ],
            &[
                (
                    "window and recorded position, both",
                    &window_here,
                    &recorded,
                ),
                ("window and exact position, both", &window_here, &exact),
            ],
        );
    }

    println!(
        "\n  the sparse GC curve was fitted from {sparse_positions} positions, 1 in \
         {SPARSE_CURVE_STRIDE} of the walk"
    );

    // Does either discriminator add anything to the other? The near-half rate in each cell of
    // window band × position band answers it without a threshold argument.
    println!("\nnear-half rate by window band and position band, positions at depth 4 or more");
    print_cross_table(&observations, &window_here, &recorded, 4);

    // The shape behind the gate: how the alternative fraction is distributed inside each
    // band. A duplication class shows as a bump at a half that the one-copy band does not
    // have.
    println!("\nalternative-read fraction by band, positions at depth 8 or more");
    print_alt_fraction_shape(&observations, &relative, 8);

    if let Some(dump) = dump {
        write_windows(
            &dump,
            &windows,
            &relative,
            &relative_raw,
            &observations,
            &contigs,
            8,
        );
        println!("\nper-window dump  {}", dump.display());
    }

    println!(
        "\ntotal            {:.1} s",
        started.elapsed().as_secs_f64()
    );
}

/// One row per window: what this sample made of it. **The columns two samples are compared
/// on are `relative` and `near_half`** — the rest are there so a disagreement can be traced
/// to a depth, a GC fraction or a thin denominator rather than left as a disagreement.
fn write_windows(
    path: &std::path::Path,
    windows: &[Window],
    relative: &[f64],
    relative_raw: &[f64],
    observations: &[PositionObservation],
    contigs: &pop_var_caller::fasta::ContigList,
    min_depth: u16,
) {
    let mut half_count = vec![0_u32; windows.len()];
    let mut scored = vec![0_u32; windows.len()];
    for observation in observations {
        if observation.depth < min_depth {
            continue;
        }
        scored[observation.window as usize] += 1;
        if near_half(observation) {
            half_count[observation.window as usize] += 1;
        }
    }
    let mut out = String::from(
        "contig\tstart\tgeneric_positions\tmean_depth\tgc\trelative\trelative_raw\tscored\tnear_half\n",
    );
    for (i, window) in windows.iter().enumerate() {
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\n",
            contigs.entries[window.contig as usize].name,
            window.start,
            window.generic_positions,
            window.mean_depth(),
            window.gc_fraction().unwrap_or(f64::NAN),
            relative[i],
            relative_raw[i],
            scored[i],
            half_count[i],
        ));
    }
    std::fs::write(path, out).expect("the dump path is writable");
}

fn window_key(contig: u32, window: u64) -> u64 {
    (u64::from(contig) << 40) | window
}

fn gc_bin(fraction: f64) -> usize {
    ((fraction * (GC_BINS - 1) as f64).round() as usize).min(GC_BINS - 1)
}

fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).expect("no NaN among window depths"));
    values[values.len() / 2]
}

fn rescale_to_unit_median(values: Vec<f64>) -> Vec<f64> {
    let centre = median(&mut values.clone());
    if centre <= 0.0 {
        return values;
    }
    values.into_iter().map(|v| v / centre).collect()
}

fn percentile_line(values: &[f64]) -> String {
    if values.is_empty() {
        return "no windows".to_string();
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    let at = |q: f64| sorted[((sorted.len() - 1) as f64 * q).round() as usize];
    format!(
        "p1 {:.3}  p10 {:.3}  p50 {:.3}  p90 {:.3}  p99 {:.3}  p99.9 {:.3}  max {:.3}",
        at(0.01),
        at(0.10),
        at(0.50),
        at(0.90),
        at(0.99),
        at(0.999),
        sorted[sorted.len() - 1]
    )
}

/// The bands relative coverage is read in. `[low, high)` on the corrected value.
const BANDS: [(f64, f64, &str); 6] = [
    (0.0, 0.6, "below 0.6 — half a copy or less"),
    (0.6, 1.4, "0.6 to 1.4 — one copy"),
    (1.4, 1.6, "1.4 to 1.6 — between"),
    (1.6, 2.4, "1.6 to 2.4 — two copies"),
    (2.4, 3.5, "2.4 to 3.5 — three copies"),
    (3.5, f64::INFINITY, "3.5 and above"),
];

fn band_of(value: f64) -> Option<usize> {
    if value.is_nan() {
        return None;
    }
    BANDS
        .iter()
        .position(|&(low, high, _)| value >= low && value < high)
}

fn print_band_table(windows: &[Window], relative: &[f64], label: &str) {
    let mut window_count = [0_u64; BANDS.len()];
    let mut position_count = [0_u64; BANDS.len()];
    for (i, window) in windows.iter().enumerate() {
        if let Some(band) = band_of(relative[i]) {
            window_count[band] += 1;
            position_count[band] += u64::from(window.generic_positions);
        }
    }
    let windows_total: u64 = window_count.iter().sum();
    let positions_total: u64 = position_count.iter().sum();
    println!("  {label}");
    println!("    band                             windows          %      positions          %");
    for (band, (_, _, name)) in BANDS.iter().enumerate() {
        println!(
            "    {name:<30}  {:>10}  {:>8.3}  {:>12}  {:>8.3}",
            window_count[band],
            100.0 * window_count[band] as f64 / windows_total as f64,
            position_count[band],
            100.0 * position_count[band] as f64 / positions_total as f64
        );
    }
}

/// Is this position reading about half alternative? Two reads at least, so one stray read
/// at low depth is not a half.
fn near_half(observation: &PositionObservation) -> bool {
    near_half_of(u32::from(observation.alt), u32::from(observation.depth))
}

fn near_half_of(alt: u32, depth: u32) -> bool {
    if alt < 2 || depth == 0 {
        return false;
    }
    let fraction = f64::from(alt) / f64::from(depth);
    (0.35..=0.65).contains(&fraction)
}

/// The GC fraction of a `(G+C, A+T)` count, `None` where the stretch is all unknown bases.
fn gc_fraction_of((gc, at): (u32, u32)) -> Option<f64> {
    let called = gc + at;
    (called > 0).then(|| f64::from(gc) / f64::from(called))
}

/// The depth at which half the counted positions sit at or below.
fn median_of_histogram(counts: &[u64]) -> u64 {
    let total: u64 = counts.iter().sum();
    let mut seen = 0_u64;
    for (depth, &count) in counts.iter().enumerate() {
        seen += count;
        if seen * 2 >= total {
            return depth as u64;
        }
    }
    0
}

/// **What the record says a position's depth was.** The stored value is a five-bit code
/// standing for a range of read counts, so the best a reader can do is a point inside that
/// range; the midpoint is that point. Below nine reads the range is a single depth and this
/// is exact.
fn recorded_depth(edges: &DepthBinEdges, depth: u32) -> f64 {
    let range = edges.depth_range(edges.bin_for(depth.min(edges.max_depth())));
    f64::from(*range.start() + *range.end()) / 2.0
}

/// How many of `reads` survive when a position's `depth` reads are thinned to `cap` — the
/// same proportional rule the record builder applies, copied so the probe can ask what the
/// cap does to an alternative fraction.
fn thin_to_cap(reads: u32, depth: u32, cap: u32) -> u32 {
    if depth <= cap || depth == 0 {
        return reads;
    }
    let scaled = f64::from(reads) * f64::from(cap) / f64::from(depth);
    (scaled.round() as u32).min(cap).max(u32::from(reads > 0))
}

/// The two depths at which the stored code stops separating one copy from two.
///
/// Returns the one-copy depth from which a doubled position's recorded count no longer
/// reaches 1.6 times an ordinary one's — the band edge the two-copy class is read at — and
/// the one-copy depth from which the two are written as the identical code. Both are
/// properties of the ladder and the cap alone, so they are the same number for every sample;
/// what changes between samples is whether the sample's own depth has reached them.
fn blinding_depths(edges: &DepthBinEdges) -> (u32, u32) {
    let ratio = |depth: u32| recorded_depth(edges, 2 * depth) / recorded_depth(edges, depth);
    let ceiling = 4 * edges.max_depth();
    let dulls_at = (1..=ceiling)
        .rev()
        .find(|&depth| ratio(depth) >= 1.6)
        .map_or(1, |depth| depth + 1);
    let blind_at = (1..=ceiling)
        .rev()
        .find(|&depth| ratio(depth) > 1.0)
        .map_or(1, |depth| depth + 1);
    (dulls_at, blind_at)
}

/// One copy's depth at each GC content, and over the whole sample.
///
/// `counts[bin][depth]` is how many positions of that GC sat at that depth; `read_as` says
/// how the reader of a depth sees it — exactly, or through the ladder. The mean rather than
/// the median because a single position's depth is a small integer, and a median over it is
/// a staircase that would flatten the very bias the correction exists for. A bin holding too
/// few positions falls back to the sample's own mean.
fn gc_curve<F: Fn(u32) -> f64>(
    counts: &[Vec<u32>],
    min_positions: u64,
    read_as: F,
) -> (Vec<f64>, f64) {
    let read: Vec<f64> = (0..=DEPTH_HISTOGRAM_LIMIT as u32).map(&read_as).collect();
    let mut global_sum = 0.0;
    let mut global_positions = 0_u64;
    let mut curve = vec![0.0; counts.len()];
    for (bin, depths) in counts.iter().enumerate() {
        let mut sum = 0.0;
        let mut positions = 0_u64;
        for (depth, &count) in depths.iter().enumerate() {
            sum += read[depth] * f64::from(count);
            positions += u64::from(count);
        }
        global_sum += sum;
        global_positions += positions;
        curve[bin] = if positions >= min_positions {
            sum / positions as f64
        } else {
            f64::NAN
        };
    }
    let global = if global_positions == 0 {
        1.0
    } else {
        global_sum / global_positions as f64
    };
    for value in &mut curve {
        if value.is_nan() || *value <= 0.0 {
            *value = global;
        }
    }
    (curve, global)
}

/// Each position's read count over what one copy would give at its GC, rescaled so the
/// median position sits at 1.0 — the cheap arm's answer to the window's relative coverage.
fn position_relative<F: Fn(u32) -> f64>(
    observations: &[PositionObservation],
    gc_windows: &[(u32, u32)],
    curve: &[f64],
    read_as: F,
) -> Vec<f32> {
    let mut values: Vec<f32> = observations
        .iter()
        .map(|observation| {
            let expected = gc_fraction_of(gc_windows[observation.gc_window as usize])
                .map_or(f64::NAN, |fraction| curve[gc_bin(fraction)]);
            if expected > 0.0 {
                (read_as(u32::from(observation.depth)) / expected) as f32
            } else {
                f32::NAN
            }
        })
        .collect();
    let centre = median(&mut sample_every(&values, 37));
    if centre > 0.0 {
        for value in &mut values {
            *value /= centre as f32;
        }
    }
    values
}

/// Every `stride`-th value, as `f64` — enough of a few million to find a median with.
fn sample_every(values: &[f32], stride: usize) -> Vec<f64> {
    values
        .iter()
        .step_by(stride)
        .map(|&value| f64::from(value))
        .filter(|value| !value.is_nan())
        .collect()
}

/// The arms, on one set of positions, in one table.
///
/// **Every arm scores the same positions** — those at `min_depth` or deeper whose window has
/// a relative coverage at all — so the columns are comparable. `enrichment` is the count of
/// positions an arm flags *and* that read near half, over what independence between the two
/// would give: a comparison between arms, and never a detection rate, because no truth set
/// says which positions are duplicated.
fn report_arms(
    observations: &[PositionObservation],
    scorable: &[f32],
    min_depth: u16,
    arms: &[(&str, &[f32])],
    combinations: &[(&str, &[f32], &[f32])],
) {
    let scored: Vec<usize> = (0..observations.len())
        .filter(|&i| observations[i].depth >= min_depth && !scorable[i].is_nan())
        .collect();
    let half: Vec<bool> = scored
        .iter()
        .map(|&i| near_half(&observations[i]))
        .collect();
    let half_total = half.iter().filter(|&&is_half| is_half).count();
    println!(
        "\n\nthe arms at depth {min_depth} or more — {} positions scored, {half_total} of \
         them near half",
        scored.len()
    );
    println!(
        "  arm                                          flagged        %   near half   \
         rate in flag   rate elsewhere   enrichment"
    );
    let row = |name: &str, flags: &dyn Fn(usize) -> bool| {
        let mut flagged = 0_u64;
        let mut both = 0_u64;
        for (slot, &i) in scored.iter().enumerate() {
            if flags(i) {
                flagged += 1;
                if half[slot] {
                    both += 1;
                }
            }
        }
        let elsewhere = half_total as u64 - both;
        let outside = scored.len() as u64 - flagged;
        let expected = flagged as f64 * half_total as f64 / scored.len() as f64;
        println!(
            "  {name:<42}  {flagged:>9}  {:>7.3}  {both:>9}  {:>12.4}%  {:>14.4}%  {:>10.2}×",
            100.0 * flagged as f64 / scored.len() as f64,
            100.0 * both as f64 / flagged.max(1) as f64,
            100.0 * elsewhere as f64 / outside.max(1) as f64,
            both as f64 / expected.max(1e-9),
        );
    };
    for (name, values) in arms {
        row(name, &|i| {
            band_of(f64::from(values[i])) == Some(TWO_COPY_BAND)
        });
    }
    for (name, first, second) in combinations {
        row(name, &|i| {
            band_of(f64::from(first[i])) == Some(TWO_COPY_BAND)
                && band_of(f64::from(second[i])) == Some(TWO_COPY_BAND)
        });
    }
}

/// The near-half rate in each cell of window band × position band.
///
/// **This is the question the enrichment table cannot answer**: whether a discriminator adds
/// anything to the other one. If the position's own count carried what the window carries,
/// the rate would rise along the position axis inside a fixed window band; if it carries
/// nothing, the rate is flat along that axis and set entirely by the window.
fn print_cross_table(
    observations: &[PositionObservation],
    window: &[f32],
    position: &[f32],
    min_depth: u16,
) {
    const CELLS: [(f64, f64, &str); 3] = [
        (0.6, 1.4, "one copy"),
        (1.4, 1.6, "between"),
        (1.6, 2.4, "two copies"),
    ];
    let cell_of = |value: f32| {
        CELLS
            .iter()
            .position(|&(low, high, _)| f64::from(value) >= low && f64::from(value) < high)
    };
    let mut counts = [[0_u64; CELLS.len()]; CELLS.len()];
    let mut halves = [[0_u64; CELLS.len()]; CELLS.len()];
    for (i, observation) in observations.iter().enumerate() {
        if observation.depth < min_depth {
            continue;
        }
        let (Some(row), Some(column)) = (cell_of(window[i]), cell_of(position[i])) else {
            continue;
        };
        counts[row][column] += 1;
        if near_half(observation) {
            halves[row][column] += 1;
        }
    }
    print!("  window \\ position           ");
    for (_, _, name) in CELLS {
        print!("{name:>24}");
    }
    println!();
    for (row, (_, _, window_name)) in CELLS.iter().enumerate() {
        print!("  {window_name:<28}");
        for column in 0..CELLS.len() {
            print!(
                "{:>24}",
                format!(
                    "{:.4}% of {}",
                    100.0 * halves[row][column] as f64 / counts[row][column].max(1) as f64,
                    counts[row][column]
                )
            );
        }
        println!();
    }
}

fn report_joint(
    observations: &[PositionObservation],
    relative: &[f64],
    min_depth: u16,
    label: &str,
) {
    let mut in_band = 0_u64;
    let mut half = 0_u64;
    let mut both = 0_u64;
    let mut scored = 0_u64;
    let two_copy = 3_usize; // the "1.6 to 2.4" band
    let mut half_by_band = [0_u64; BANDS.len()];
    let mut scored_by_band = [0_u64; BANDS.len()];
    for observation in observations {
        if observation.depth < min_depth {
            continue;
        }
        let Some(band) = band_of(relative[observation.window as usize]) else {
            continue;
        };
        scored += 1;
        scored_by_band[band] += 1;
        let is_half = near_half(observation);
        if is_half {
            half += 1;
            half_by_band[band] += 1;
        }
        if band == two_copy {
            in_band += 1;
            if is_half {
                both += 1;
            }
        }
    }
    if scored == 0 {
        println!("\nat depth {min_depth} or more ({label}): no position scored");
        return;
    }
    let p_band = in_band as f64 / scored as f64;
    let p_half = half as f64 / scored as f64;
    let expected = p_band * p_half * scored as f64;
    println!(
        "\nat depth {min_depth} or more ({label}) — {scored} positions scored, {:.1}% of the \
         single-base ones",
        100.0 * scored as f64 / observations.len() as f64
    );
    println!(
        "  in a two-copy window                 {in_band} ({:.4}% of scored)",
        100.0 * p_band
    );
    println!(
        "  near half alternative                {half} ({:.4}% of scored)",
        100.0 * p_half
    );
    println!(
        "  both — the gate                      {both} ({:.4}% of scored, 1 in {:.0})",
        100.0 * both as f64 / scored as f64,
        if both > 0 {
            scored as f64 / both as f64
        } else {
            f64::INFINITY
        }
    );
    println!(
        "  expected if the two were independent {expected:.1} — observed is {:.2}× that",
        both as f64 / expected.max(1e-9)
    );
    println!("  near-half rate inside each band");
    for (band, (_, _, name)) in BANDS.iter().enumerate() {
        if scored_by_band[band] == 0 {
            continue;
        }
        println!(
            "    {name:<30}  {:>10} of {:>10}  {:>8.4}%",
            half_by_band[band],
            scored_by_band[band],
            100.0 * half_by_band[band] as f64 / scored_by_band[band] as f64
        );
    }
}

fn print_alt_fraction_shape(
    observations: &[PositionObservation],
    relative: &[f64],
    min_depth: u16,
) {
    const CELLS: usize = 10;
    let mut counts = vec![[0_u64; CELLS]; BANDS.len()];
    let mut totals = [0_u64; BANDS.len()];
    for observation in observations {
        if observation.depth < min_depth {
            continue;
        }
        let Some(band) = band_of(relative[observation.window as usize]) else {
            continue;
        };
        let fraction = f64::from(observation.alt) / f64::from(observation.depth);
        let cell = ((fraction * CELLS as f64) as usize).min(CELLS - 1);
        counts[band][cell] += 1;
        totals[band] += 1;
    }
    print!("  band                            ");
    for cell in 0..CELLS {
        print!(
            "{:>9}",
            format!("{:.1}-{:.1}", cell as f64 / 10.0, (cell + 1) as f64 / 10.0)
        );
    }
    println!();
    for (band, (_, _, name)) in BANDS.iter().enumerate() {
        if totals[band] == 0 {
            continue;
        }
        print!("  {name:<30}  ");
        for cell in 0..CELLS {
            print!(
                "{:>9.4}",
                100.0 * counts[band][cell] as f64 / totals[band] as f64
            );
        }
        println!("   ({} positions)", totals[band]);
    }
    println!("  (each row is a percentage of that band's positions, summing to 100)");
}
