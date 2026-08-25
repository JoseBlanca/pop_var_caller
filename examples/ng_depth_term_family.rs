//! Is a position's read count Poisson, and does a likelihood ratio separate better than a
//! threshold did?
//!
//! `parameter_prepass_joint_fit.md` §2.2 carries a third class of site — a stretch of genome
//! the sample holds twice while the reference holds it once, so both copies' reads land on one
//! position and the position reads about half non-reference wherever the copies differ. With
//! twenty-five samples or more the cohort recognises such a site for free. Below that, and at
//! one sample, the only evidence left is that the position carries about twice the reads the
//! sample normally has there, and the plan is to write that evidence into the fit as
//!
//! ```text
//!     ln P(d | 2m) − ln P(d | m)
//! ```
//!
//! where `d` is the position's read count and `m` is the depth that sample is expected to
//! reach at that position. **Two things about that formula are unsettled and this program
//! measures both.**
//!
//! # 1. Which family is `P`?
//!
//! Writing `P` as a Poisson asserts the variance of depth equals its mean. If depth is
//! overdispersed the Poisson ratio is confidently wrong in the direction that costs: it hands
//! strong evidence of *two copies* to a position that is merely deep. The measurement needs no
//! model — inside one GC bin of one sample, take **variance ÷ mean** over single positions.
//! Near 1 and Poisson is right; at 3 it is not.
//!
//! Three things the measurement has to get right, and each is a column below:
//!
//! - **the true read count, never the record's five-bit code.** The stored ladder is exact to
//!   eight reads and geometric above, so a dispersion measured on it is a property of the
//!   ladder. This program walks the alignment and holds the exact counts.
//! - **GC content held fixed.** Per-position depth on one accession runs from 11.7 reads at
//!   18% GC to 29.0 at 34% — a factor of 2.5, larger than the doubling the term exists to
//!   detect. So the variance is measured inside a GC bin and never across them.
//! - **how many positions stand behind each ratio.** A thin bin returns noise, so every bin
//!   prints its count and its own standard error, and bins under
//!   [`MIN_POSITIONS_PER_GC_BIN`] are dropped.
//!
//! Two families fit an overdispersed count and they disagree about what happens as depth
//! grows: a **negative binomial** (variance `m + m²/r`, so variance ÷ mean rises with depth)
//! and a **quasi-Poisson** (variance `φ·m`, so variance ÷ mean is one number at every depth).
//! GC bins span a 2.5-fold range of `m` within one sample, which is enough to tell them apart,
//! and the program fits both and prints which describes the bins better.
//!
//! **Mappability is not GC.** A position in a repeat-rich neighbourhood reads deep for a
//! reason neither the GC curve nor the copy-number term models. If that arrives as a few very
//! deep positions rather than as a uniformly wider spread it changes *which family* is right,
//! not just its parameter — so each bin also prints what share of its variance comes from its
//! deepest 1% of positions, beside what a Poisson and a fitted negative binomial would put
//! there.
//!
//! # 2. Does the log-ratio separate as well as a threshold did?
//!
//! Everything measured so far scored a position by a threshold — *is this about two copies?* —
//! and reported **enrichment**: how often a flagged position reads near half, over what the
//! flagging rate and the near-half rate together predict. The log-ratio is not a threshold, so
//! it is scored here as its own arm and reported on the same two quantities: enrichment, and
//! the share of scored positions it puts above zero.
//!
//! The arms are scored on identical positions so the columns can be compared. An arm that
//! flags more positions can buy enrichment with them, so the table also holds two arms scored
//! at the flagged share the threshold arm chose — that is the comparison at equal cost.
//!
//! **There is no truth set.** Nobody has a validated list of the stretches these accessions
//! carry twice. Enrichment compares discriminators on the same positions; it is not a
//! detection rate, and a flagged position is not thereby a duplication.
//!
//! ```text
//! ng_depth_term_family <reference.fa> <catalog.parquet> <alignment> <regions.bed>
//! ```
//!
//! Lines beginning `TSV` carry the same numbers in a machine-readable form, so a run over
//! several accessions can be assembled into one table.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::genetics::lgamma;
use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, LocusKind, SampleLocusObservationsIterator, UnhandledReason,
};
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

/// How wide a GC bin is, in percentage points, unless the command line says otherwise.
///
/// A bin still holds a range of GC contents, and depth varies inside it, so some of what the
/// measurement reads as *spread at fixed GC* is really the bin's own width. Narrowing the bin
/// is the check: if variance ÷ mean falls when the bin halves, the width was manufacturing it.
const DEFAULT_GC_BIN_WIDTH: usize = 2;

/// The stretch of reference a single position's GC content is read over. GC bias enters
/// through the fragments that were sequenced, so the neighbourhood a position's depth is
/// corrected against is a fragment's length — the same 500 bp the coverage summary stores at.
const GC_WINDOW_BP: u64 = 500;

/// The deepest position counted exactly; anything deeper joins the top cell. Tomato's deepest
/// single position across eight accessions is 234, so nothing reaches this.
const DEPTH_LIMIT: usize = 4_096;

/// Below this many positions a GC bin's variance is mostly its own noise and the bin is
/// dropped. At 5,000 positions the standard error of variance ÷ mean is about 2% of it, which
/// is small against the difference between 1 and 3 the measurement exists to resolve. Every
/// surviving bin still prints its own standard error, computed from its fourth moment rather
/// than assumed.
const MIN_POSITIONS_PER_GC_BIN: u64 = 5_000;

/// Positions shallower than this are not scored by any arm — the floor the earlier threshold
/// measurement used, kept so the two are comparable.
const MIN_DEPTH_SCORED: u16 = 4;

/// The alternative-read fraction that counts as *near half*, and the least number of
/// alternative reads for the fraction to mean anything. Both as in
/// `locus_depth_vs_window_2026-08-13.md`.
const NEAR_HALF: std::ops::RangeInclusive<f64> = 0.35..=0.65;
const MIN_ALT_READS: u32 = 2;

/// The threshold arm's band on relative depth, `[low, high)`: *about two copies*.
const TWO_COPY_BAND: (f64, f64) = (1.6, 2.4);

/// The strengths of evidence the report counts positions above, as the log-ratio gives them —
/// in nats, so one factor of ten in the odds is `ln 10`. Even odds, then 10, 100 and 1,000 to 1.
const EVIDENCE_CUTS: [f64; 4] = [
    0.0,
    std::f64::consts::LN_10,
    2.0 * std::f64::consts::LN_10,
    3.0 * std::f64::consts::LN_10,
];

/// One generic position, as the walk saw it. Six bytes, because there are millions.
#[derive(Copy, Clone)]
struct PositionObservation {
    /// Which GC bin the 500 bp around it falls in.
    gc_bin: u8,
    /// The exact read count — not the record's five-bit code.
    depth: u16,
    /// Reads disagreeing with the reference.
    alt: u16,
}

impl PositionObservation {
    fn near_half(&self) -> bool {
        if u32::from(self.alt) < MIN_ALT_READS || self.depth == 0 {
            return false;
        }
        NEAR_HALF.contains(&(f64::from(self.alt) / f64::from(self.depth)))
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let usage =
        "usage: <reference.fa> <catalog.parquet> <alignment> <regions.bed> [gc-bin-width-%]";
    let fasta = PathBuf::from(args.next().expect(usage));
    let catalog_path = PathBuf::from(args.next().expect(usage));
    let reads = PathBuf::from(args.next().expect(usage));
    let bed = PathBuf::from(args.next().expect(usage));
    let gc_bin_width: usize = args.next().map_or(DEFAULT_GC_BIN_WIDTH, |a| {
        a.parse()
            .expect("a GC bin width in whole percentage points")
    });
    assert!(
        gc_bin_width > 0 && 100_usize.is_multiple_of(gc_bin_width),
        "a GC bin width has to divide 100 percentage points"
    );
    let gc_bins = 100 / gc_bin_width + 1;
    let sample_name = reads
        .file_stem()
        .map_or_else(|| "sample".to_string(), |s| s.to_string_lossy().to_string());

    let started = Instant::now();

    // ---- the reference, the catalog, the regions ---------------------------------------
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

    println!("sample           {sample_name}");
    println!("alignment        {}", reads.display());
    println!("GC bin           {gc_bin_width} percentage points wide");
    println!(
        "generic regions  {} spans, {generic_bases} bases",
        generic.len()
    );

    // ---- pass 1: each 500 bp stretch's GC content, off the reference ---------------------
    let accessor = WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone());
    let mut gc_window_index: HashMap<u64, u32> = HashMap::new();
    let mut gc_windows: Vec<(u32, u32)> = Vec::new();
    // How many generic positions each stretch holds. **The walk emits no locus at a position no
    // read reached**, so without this the depth distribution is the one with its zeros removed
    // — which is narrower than the distribution itself, and at three reads a position narrower
    // by enough to read as under-dispersion. The denominator comes off the reference, so the
    // missing positions can be put back.
    let mut generic_positions: Vec<u32> = Vec::new();
    let mut bases = Vec::new();
    for region in &generic {
        accessor
            .fetch_into(region.contig, region.start.get(), region.len(), &mut bases)
            .expect("a generic region reads from the reference");
        for (offset, base) in bases.iter().enumerate() {
            let position = region.start.get() + offset as u64;
            let key = gc_key(region.contig.0, position / GC_WINDOW_BP);
            let slot = *gc_window_index.entry(key).or_insert_with(|| {
                gc_windows.push((0, 0));
                generic_positions.push(0);
                (gc_windows.len() - 1) as u32
            });
            generic_positions[slot as usize] += 1;
            match base.to_ascii_uppercase() {
                b'G' | b'C' => gc_windows[slot as usize].0 += 1,
                b'A' | b'T' => gc_windows[slot as usize].1 += 1,
                _ => {}
            }
        }
    }

    // ---- pass 2: the walk ----------------------------------------------------------------
    let read_groups = build_read_groups(std::slice::from_ref(&reads))
        .expect("the header declares its read groups");
    let sample = match read_groups.read_groups_per_sample() {
        [only] => only.clone(),
        other => panic!(
            "{} holds {} samples; this program is per sample",
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
            reason = "file-backed and single-threaded, as in the duplication probe"
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
    let mut walked = 0_u64;
    let mut total_depth = 0_u64;
    // Every generic position the walk reached, whatever the locus's width — the number the
    // reference's own count is compared against to find the positions no read reached.
    let mut reached: Vec<u32> = vec![0; gc_windows.len()];
    for locus in &mut stream {
        let locus = locus.expect("the walk runs to completion on a well-formed alignment");
        if locus.kind != LocusKind::Generic {
            continue;
        }
        for offset in 0..locus.region.len() {
            let key = gc_key(
                locus.region.contig.0,
                (locus.region.start.get() + offset) / GC_WINDOW_BP,
            );
            if let Some(&slot) = gc_window_index.get(&key) {
                reached[slot as usize] += 1;
            }
        }
        // Only a one-position locus answers "what fraction of the reads here disagree" — at a
        // wider one an observation's bases span the whole region. Wider loci are 1 in 1,100 of
        // the walk on tomato.
        if locus.region.len() != 1 {
            continue;
        }
        let alt = locus
            .complete_observations()
            .filter(|obs| *obs.bases != *locus.reference_bases)
            .map(|obs| obs.num_obs)
            .sum::<u32>();
        let depth = locus.num_obs_along_locus()[0];
        let key = gc_key(
            locus.region.contig.0,
            locus.region.start.get() / GC_WINDOW_BP,
        );
        let slot = *gc_window_index
            .get(&key)
            .expect("a walked generic position sits in a stretch pass 1 filled");
        let Some(fraction) = gc_fraction_of(gc_windows[slot as usize]) else {
            continue;
        };
        walked += 1;
        total_depth += u64::from(depth);
        observations.push(PositionObservation {
            gc_bin: gc_bin(fraction, gc_bins) as u8,
            depth: u16::try_from(depth).unwrap_or(u16::MAX),
            alt: u16::try_from(alt).unwrap_or(u16::MAX),
        });
    }
    if let Some(handle) = verification {
        handle
            .join()
            .expect("the .fai beside the reference describes it");
    }
    let mean_depth = total_depth as f64 / walked as f64;
    println!(
        "walk             {walked} single-base generic positions, mean depth {mean_depth:.2} \
         reads a position, {:.1} s",
        started.elapsed().as_secs_f64()
    );

    // ---- measurement 1: is depth Poisson? --------------------------------------------------
    let all = histograms(&observations, gc_bins, |_| true);
    let ordinary = histograms(&observations, gc_bins, |observation| {
        !observation.near_half()
    });
    // The positions no read reached, put back at depth zero. Without them the low-coverage
    // accessions read as *under*-dispersed, which is what a Poisson looks like once its zeros
    // have been removed and not a property of the data.
    let mut with_zeros = all.clone();
    let mut unreached = 0_u64;
    for (slot, &here) in generic_positions.iter().enumerate() {
        let Some(fraction) = gc_fraction_of(gc_windows[slot]) else {
            continue;
        };
        let missing = here.saturating_sub(reached[slot]);
        unreached += u64::from(missing);
        with_zeros[gc_bin(fraction, gc_bins)][0] += missing;
    }
    println!(
        "                 {unreached} generic positions no read reached ({:.2}% of them), put \
         back at depth zero for the variance",
        100.0 * unreached as f64 / generic_bases as f64
    );
    let bins: Vec<usize> = (0..gc_bins)
        .filter(|&bin| positions_in(&all[bin]) >= MIN_POSITIONS_PER_GC_BIN)
        .collect();
    let stats: HashMap<usize, BinStatistics> = bins
        .iter()
        .map(|&bin| (bin, BinStatistics::of(&all[bin])))
        .collect();
    let stats_with_zeros: HashMap<usize, BinStatistics> = bins
        .iter()
        .map(|&bin| (bin, BinStatistics::of(&with_zeros[bin])))
        .collect();

    println!("\n================ 1. variance against mean, inside a GC bin ================\n");
    println!(
        "  Bins under {MIN_POSITIONS_PER_GC_BIN} positions are dropped. `± ` is the standard \
         error of variance ÷ mean,"
    );
    println!(
        "  computed from the bin's own fourth moment. `ordinary` drops the near-half positions \
         — the ones"
    );
    println!(
        "  a duplication would put there — so the spread is not being read off the artefact \
         itself."
    );
    println!(
        "  `with zeros` adds back the positions no read reached; `ordinary` drops the near-half \
         ones."
    );
    println!(
        "\n  GC%   positions    mean d   variance   var/mean          with zeros   ordinary  \
         implied r   deepest 1% of positions"
    );
    println!(
        "                                                                                     \
                hold this share of the variance"
    );
    println!(
        "                                                                                     \
             observed   Poisson   neg-binom"
    );
    for &bin in &bins {
        let s = &stats[&bin];
        let ordinary_ratio = BinStatistics::of(&ordinary[bin]).variance_to_mean();
        let with_zero_ratio = stats_with_zeros[&bin].variance_to_mean();
        let implied = s.implied_negative_binomial_size();
        let (poisson_tail, nb_tail) = s.variance_share_of_deepest_hundredth_under_models();
        println!(
            "  {:>3}  {:>10}  {:>8.3}  {:>9.3}  {:>6.3} ± {:>5.3}  {:>10.3}  {:>9.3}  {:>9}  \
             {:>8.1}% {:>8.1}%  {:>9.1}%",
            bin * gc_bin_width,
            s.positions,
            s.mean,
            s.variance,
            s.variance_to_mean(),
            s.variance_to_mean_standard_error(),
            with_zero_ratio,
            ordinary_ratio,
            implied.map_or("—".to_string(), |r| format!("{r:.1}")),
            100.0 * s.variance_share_of_deepest_hundredth(),
            100.0 * poisson_tail,
            100.0 * nb_tail,
        );
        println!(
            "TSV\tbin\t{sample_name}\t{mean_depth:.4}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t\
             {:.4}\t{:.4}\t{}\t{:.4}\t{:.4}\t{:.4}",
            bin * gc_bin_width,
            s.positions,
            s.mean,
            s.variance,
            s.variance_to_mean(),
            s.variance_to_mean_standard_error(),
            with_zero_ratio,
            ordinary_ratio,
            implied.unwrap_or(f64::NAN),
            s.variance_share_of_deepest_hundredth(),
            poisson_tail,
            nb_tail,
        );
    }

    // Which of the two overdispersed families describes the bins? A negative binomial says
    // variance ÷ mean rises with depth as `1 + m/r`; a quasi-Poisson says it is one number at
    // every depth. The GC bins span a 2.5-fold range of mean depth inside one sample, which is
    // the lever that tells them apart.
    let fit = FamilyFit::over(&bins, &stats);
    fit.report("positions a read reached", &sample_name, mean_depth);
    FamilyFit::over(&bins, &stats_with_zeros).report(
        "the same positions with the unreached ones back at zero",
        &sample_name,
        mean_depth,
    );

    // What each family predicts for the thing that matters: how often an ordinary position
    // reads as deep as a doubled one. This is the "confidently wrong" quantity — a Poisson that
    // under-predicts the 2× share is a Poisson that will call those positions duplicated.
    println!("\n  how often an ordinary position reads deep, observed against each family");
    println!(
        "  (share of the bin's positions at or above 1.5×, 2× and 3× its mean depth, in \
         positions per 10,000)"
    );
    println!(
        "  GC%    mean d        at 1.5× mean                  at 2× mean                    \
         at 3× mean"
    );
    println!(
        "                 observed Poisson neg-bin      observed Poisson neg-bin      \
         observed Poisson neg-bin"
    );
    for &bin in &bins {
        let s = &stats[&bin];
        let r = fit.negative_binomial_size;
        print!("  {:>3}  {:>8.2}  ", bin * gc_bin_width, s.mean);
        for multiple in [1.5, 2.0, 3.0] {
            let threshold = (s.mean * multiple).ceil() as u32;
            print!(
                "{:>9.1} {:>7.1} {:>7.1}      ",
                10_000.0 * s.share_at_or_above(threshold),
                10_000.0 * poisson_tail_share(s.mean, threshold),
                10_000.0 * negative_binomial_tail_share(s.mean, r, threshold),
            );
        }
        println!();
        println!(
            "TSV\ttail\t{sample_name}\t{}\t{:.4}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t\
             {:.6}\t{:.6}\t{:.6}",
            bin * gc_bin_width,
            s.mean,
            s.share_at_or_above((s.mean * 1.5).ceil() as u32),
            poisson_tail_share(s.mean, (s.mean * 1.5).ceil() as u32),
            negative_binomial_tail_share(s.mean, r, (s.mean * 1.5).ceil() as u32),
            s.share_at_or_above((s.mean * 2.0).ceil() as u32),
            poisson_tail_share(s.mean, (s.mean * 2.0).ceil() as u32),
            negative_binomial_tail_share(s.mean, r, (s.mean * 2.0).ceil() as u32),
            s.share_at_or_above((s.mean * 3.0).ceil() as u32),
            poisson_tail_share(s.mean, (s.mean * 3.0).ceil() as u32),
            negative_binomial_tail_share(s.mean, r, (s.mean * 3.0).ceil() as u32),
        );
    }

    // ---- measurement 2: the log-ratio against the threshold --------------------------------
    //
    // `m` is what one copy is expected to give at this position's GC content. Two readings of
    // it are carried: the bin's mean, which is what an estimator with no truth set can compute,
    // and the bin's median, which a heavy upper tail moves far less. They differ by exactly the
    // amount the tail pulls the mean, so printing both says whether the choice matters.
    let expected_from_mean: Vec<f64> = (0..gc_bins)
        .map(|bin| stats.get(&bin).map_or(f64::NAN, |s| s.mean))
        .collect();
    let expected_from_median: Vec<f64> = (0..gc_bins)
        .map(|bin| stats.get(&bin).map_or(f64::NAN, |s| s.median))
        .collect();

    println!("\n\n================ 2. the log-ratio against the threshold ================");
    scored_report(
        &observations,
        &expected_from_mean,
        fit.negative_binomial_size,
        &sample_name,
        mean_depth,
        "one copy's depth read as the GC bin's mean",
        true,
    );
    scored_report(
        &observations,
        &expected_from_median,
        fit.negative_binomial_size,
        &sample_name,
        mean_depth,
        "one copy's depth read as the GC bin's median",
        false,
    );

    println!(
        "\ntotal            {:.1} s",
        started.elapsed().as_secs_f64()
    );
}

/// Every arm, on one set of positions, with the two quantities the earlier threshold
/// measurement reported: enrichment, and the share of scored positions the arm flags.
///
/// `verbose` adds the evidence-weight table, which is a property of the family and not of the
/// `m` reading, so it is printed once.
fn scored_report(
    observations: &[PositionObservation],
    expected: &[f64],
    negative_binomial_size: f64,
    sample_name: &str,
    mean_depth: f64,
    label: &str,
    verbose: bool,
) {
    let scored: Vec<usize> = (0..observations.len())
        .filter(|&i| {
            observations[i].depth >= MIN_DEPTH_SCORED
                && expected[observations[i].gc_bin as usize].is_finite()
        })
        .collect();
    let near_half: Vec<bool> = scored
        .iter()
        .map(|&i| observations[i].near_half())
        .collect();
    let half_total = near_half.iter().filter(|&&is| is).count();

    // The three scores, one per scored position.
    let relative: Vec<f64> = scored
        .iter()
        .map(|&i| f64::from(observations[i].depth) / expected[observations[i].gc_bin as usize])
        .collect();
    // The threshold arm read relative depth after rescaling the *median position of the whole
    // walk* to 1.0 — not the median of the positions it went on to score, which at three reads
    // a position is a deep subset and would move the band by nearly a factor of two. Rescaling
    // the same way is what makes this a reproduction of
    // `locus_depth_vs_window_2026-08-13.md` rather than something that resembles it.
    let centre = median_of(
        &observations
            .iter()
            .map(|observation| f64::from(observation.depth) / expected[observation.gc_bin as usize])
            .filter(|value| value.is_finite())
            .collect::<Vec<f64>>(),
    );
    let rescaled: Vec<f64> = relative.iter().map(|value| value / centre).collect();
    println!(
        "\n  ({label}; relative depth rescaled by {centre:.4} so the median position of the \
         whole walk sits at 1.0)"
    );
    let poisson: Vec<f64> = scored
        .iter()
        .map(|&i| {
            poisson_log_ratio(
                f64::from(observations[i].depth),
                expected[observations[i].gc_bin as usize],
            )
        })
        .collect();
    let negative_binomial: Vec<f64> = scored
        .iter()
        .map(|&i| {
            negative_binomial_log_ratio(
                f64::from(observations[i].depth),
                expected[observations[i].gc_bin as usize],
                negative_binomial_size,
                negative_binomial_size,
            )
        })
        .collect();
    // The other reading of a doubled position: if the two copies' read counts are drawn
    // independently, the doubled position is a negative binomial of twice the mean *and* twice
    // the size, which is a narrower distribution than one of the same size.
    let negative_binomial_independent: Vec<f64> = scored
        .iter()
        .map(|&i| {
            negative_binomial_log_ratio(
                f64::from(observations[i].depth),
                expected[observations[i].gc_bin as usize],
                negative_binomial_size,
                2.0 * negative_binomial_size,
            )
        })
        .collect();

    println!(
        "\n--- {label} — {} positions scored, {half_total} near half",
        scored.len()
    );
    println!(
        "  arm                                                flagged        %   near half    \
         rate in flag   rate elsewhere   enrichment"
    );
    let mut band_flagged = 0_usize;
    let mut rows: Vec<(String, Vec<bool>)> = Vec::new();
    rows.push((
        format!(
            "threshold, {:.1} ≤ depth/m < {:.1}",
            TWO_COPY_BAND.0, TWO_COPY_BAND.1
        ),
        rescaled
            .iter()
            .map(|&v| v >= TWO_COPY_BAND.0 && v < TWO_COPY_BAND.1)
            .collect(),
    ));
    rows.push((
        format!("threshold, depth/m ≥ {:.1}, no upper edge", TWO_COPY_BAND.0),
        rescaled.iter().map(|&v| v >= TWO_COPY_BAND.0).collect(),
    ));
    // The same band without the rescaling, which is what the log-ratio arms work on: `m` used
    // as the expected depth it is, and nothing recentred. The gap between this row and the
    // first says how much of the threshold arm's answer came from the recentring.
    rows.push((
        format!(
            "threshold, {:.1} ≤ depth/m < {:.1}, no rescaling",
            TWO_COPY_BAND.0, TWO_COPY_BAND.1
        ),
        relative
            .iter()
            .map(|&v| v >= TWO_COPY_BAND.0 && v < TWO_COPY_BAND.1)
            .collect(),
    ));
    rows.push((
        "Poisson log-ratio above zero".to_string(),
        poisson.iter().map(|&v| v > 0.0).collect(),
    ));
    rows.push((
        "negative-binomial log-ratio above zero".to_string(),
        negative_binomial.iter().map(|&v| v > 0.0).collect(),
    ));
    rows.push((
        "negative-binomial, copies drawn independently".to_string(),
        negative_binomial_independent
            .iter()
            .map(|&v| v > 0.0)
            .collect(),
    ));
    for (index, (name, flags)) in rows.iter().enumerate() {
        let outcome = enrichment_of(flags, &near_half, half_total);
        if index == 0 {
            band_flagged = outcome.flagged;
        }
        print_arm(name, &outcome, scored.len());
        println!(
            "TSV\tarm\t{sample_name}\t{mean_depth:.4}\t{label}\t{name}\t{}\t{}\t{:.6}\t{:.6}",
            outcome.flagged,
            outcome.both,
            outcome.flagged as f64 / scored.len() as f64,
            outcome.enrichment
        );
    }

    // The comparison at equal cost: give each log-ratio exactly as many positions as the
    // threshold arm flagged, taking its highest-scoring ones. An arm that flags more positions
    // can always buy enrichment with them, and this removes that.
    for (name, score) in [
        ("Poisson log-ratio, top scores", &poisson),
        (
            "negative-binomial log-ratio, top scores",
            &negative_binomial,
        ),
    ] {
        let flags = top_k(score, band_flagged);
        let outcome = enrichment_of(&flags, &near_half, half_total);
        print_arm(
            &format!("{name} ({band_flagged} of them, as the threshold arm)"),
            &outcome,
            scored.len(),
        );
        println!(
            "TSV\tarm\t{sample_name}\t{mean_depth:.4}\t{label}\t{name}, matched share\t{}\t{}\t\
             {:.6}\t{:.6}",
            outcome.flagged,
            outcome.both,
            outcome.flagged as f64 / scored.len() as f64,
            outcome.enrichment
        );
    }

    if !verbose {
        return;
    }

    // What the term would actually add to the fit's likelihood. The flagging share says which
    // positions an arm picks; this says how loudly it argues for them, and it is where the two
    // families part company.
    println!(
        "\n  how much evidence of two copies each family hands out, in positions per 10,000 \
         scored"
    );
    println!(
        "  odds of two copies over one:            better than even   10 to 1   100 to 1   \
         1,000 to 1"
    );
    for (name, score) in [
        ("Poisson", &poisson),
        ("negative binomial, same spread", &negative_binomial),
        (
            "negative binomial, copies independent",
            &negative_binomial_independent,
        ),
    ] {
        print!("  {name:<38}");
        for cut in EVIDENCE_CUTS {
            let share = score.iter().filter(|&&value| value > cut).count() as f64
                / score.len().max(1) as f64;
            print!("{:>18.1}", 10_000.0 * share);
        }
        println!();
        println!(
            "TSV\tevidence\t{sample_name}\t{mean_depth:.4}\t{name}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            share_above(score, EVIDENCE_CUTS[0]),
            share_above(score, EVIDENCE_CUTS[1]),
            share_above(score, EVIDENCE_CUTS[2]),
            share_above(score, EVIDENCE_CUTS[3]),
        );
    }

    // The same weights, separated by whether the position reads near half. The difference is
    // what the fit sees; the level of the second column is what it costs everywhere else.
    println!(
        "\n  mean evidence, in nats, at a near-half position and at an ordinary one \
         (positive argues for two copies)"
    );
    for (name, score) in [
        ("Poisson", &poisson),
        ("negative binomial, same spread", &negative_binomial),
        (
            "negative binomial, copies independent",
            &negative_binomial_independent,
        ),
    ] {
        let (mut half_sum, mut half_count, mut rest_sum, mut rest_count) = (0.0, 0_u64, 0.0, 0_u64);
        for (slot, &value) in score.iter().enumerate() {
            if near_half[slot] {
                half_sum += value;
                half_count += 1;
            } else {
                rest_sum += value;
                rest_count += 1;
            }
        }
        println!(
            "  {name:<38}near half {:>8.2}   ordinary {:>8.2}   difference {:>8.2}",
            half_sum / half_count.max(1) as f64,
            rest_sum / rest_count.max(1) as f64,
            half_sum / half_count.max(1) as f64 - rest_sum / rest_count.max(1) as f64,
        );
    }
}

fn share_above(score: &[f64], cut: f64) -> f64 {
    score.iter().filter(|&&value| value > cut).count() as f64 / score.len().max(1) as f64
}

/// What one arm did: how many positions it flagged and how much more often those read near
/// half than independence would give.
struct ArmOutcome {
    flagged: usize,
    both: usize,
    rate_in_flag: f64,
    rate_elsewhere: f64,
    enrichment: f64,
}

fn enrichment_of(flags: &[bool], near_half: &[bool], half_total: usize) -> ArmOutcome {
    let mut flagged = 0_usize;
    let mut both = 0_usize;
    for (slot, &is_flagged) in flags.iter().enumerate() {
        if is_flagged {
            flagged += 1;
            if near_half[slot] {
                both += 1;
            }
        }
    }
    let scored = flags.len();
    let outside = scored - flagged;
    let expected = flagged as f64 * half_total as f64 / scored.max(1) as f64;
    ArmOutcome {
        flagged,
        both,
        rate_in_flag: both as f64 / flagged.max(1) as f64,
        rate_elsewhere: (half_total - both) as f64 / outside.max(1) as f64,
        enrichment: both as f64 / expected.max(1e-9),
    }
}

fn print_arm(name: &str, outcome: &ArmOutcome, scored: usize) {
    println!(
        "  {name:<48}  {:>9}  {:>7.3}  {:>9}  {:>12.4}%  {:>14.4}%  {:>10.2}×",
        outcome.flagged,
        100.0 * outcome.flagged as f64 / scored as f64,
        outcome.both,
        100.0 * outcome.rate_in_flag,
        100.0 * outcome.rate_elsewhere,
        outcome.enrichment,
    );
}

/// The `k` highest-scoring positions, as a flag per position.
fn top_k(score: &[f64], k: usize) -> Vec<bool> {
    let mut order: Vec<usize> = (0..score.len()).collect();
    order.sort_unstable_by(|&a, &b| score[b].total_cmp(&score[a]));
    let mut flags = vec![false; score.len()];
    for &index in order.iter().take(k) {
        flags[index] = true;
    }
    flags
}

fn median_of(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 1.0;
    }
    let mut sorted: Vec<f64> = values.iter().step_by(37).copied().collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    sorted[sorted.len() / 2]
}

/// `ln P(d | 2m) − ln P(d | m)` under a Poisson: the log-ratio reduces to `d ln 2 − m`, so it
/// crosses zero at `d = m / ln 2 = 1.44 m` and gains `ln 2` for every further read.
fn poisson_log_ratio(depth: f64, mean: f64) -> f64 {
    depth * std::f64::consts::LN_2 - mean
}

/// The same ratio under a negative binomial of mean `m` and size `r` against one of mean `2m`
/// and size `doubled_size`.
///
/// The size a doubled position takes is a modelling choice and both are measured here. If the
/// overdispersion is a local property the two copies share — mappability, a GC residual — the
/// doubled position keeps the same size, and its spread doubles with its mean. If the two
/// copies' read counts are drawn independently, the doubled position has twice the size, which
/// is a relatively narrower distribution.
fn negative_binomial_log_ratio(depth: f64, mean: f64, size: f64, doubled_size: f64) -> f64 {
    negative_binomial_log_pmf(depth, 2.0 * mean, doubled_size)
        - negative_binomial_log_pmf(depth, mean, size)
}

fn negative_binomial_log_pmf(depth: f64, mean: f64, size: f64) -> f64 {
    lgamma(depth + size) - lgamma(size) - lgamma(depth + 1.0)
        + size * (size / (size + mean)).ln()
        + depth * (mean / (size + mean)).ln()
}

fn poisson_log_pmf(depth: f64, mean: f64) -> f64 {
    depth * mean.ln() - mean - lgamma(depth + 1.0)
}

/// `P(d ≥ threshold)` under a Poisson of this mean.
fn poisson_tail_share(mean: f64, threshold: u32) -> f64 {
    let ceiling = ((mean + 30.0 * mean.sqrt()).ceil() as u32).max(threshold + 1);
    (threshold..=ceiling)
        .map(|depth| poisson_log_pmf(f64::from(depth), mean).exp())
        .sum()
}

/// `P(d ≥ threshold)` under a negative binomial of this mean and size.
fn negative_binomial_tail_share(mean: f64, size: f64, threshold: u32) -> f64 {
    let sd = (mean + mean * mean / size).sqrt();
    let ceiling = ((mean + 60.0 * sd).ceil() as u32).max(threshold + 1);
    (threshold..=ceiling)
        .map(|depth| negative_binomial_log_pmf(f64::from(depth), mean, size).exp())
        .sum()
}

/// One depth histogram per GC bin, over the positions the filter admits.
fn histograms<F: Fn(&PositionObservation) -> bool>(
    observations: &[PositionObservation],
    gc_bins: usize,
    admit: F,
) -> Vec<Vec<u32>> {
    let mut counts = vec![vec![0_u32; DEPTH_LIMIT + 1]; gc_bins];
    for observation in observations.iter().filter(|o| admit(o)) {
        let cell = (observation.depth as usize).min(DEPTH_LIMIT);
        counts[observation.gc_bin as usize][cell] += 1;
    }
    counts
}

fn positions_in(counts: &[u32]) -> u64 {
    counts.iter().map(|&count| u64::from(count)).sum()
}

/// What one GC bin's depths look like — computed from its histogram, so every moment is exact
/// rather than sampled.
struct BinStatistics {
    positions: u64,
    mean: f64,
    median: f64,
    variance: f64,
    /// The fourth central moment, which is what the variance's own standard error needs.
    fourth_moment: f64,
    /// `Σ(d − mean)²` over the deepest 1% of the bin's positions, over the whole of it.
    variance_share_of_deepest_hundredth: f64,
    counts: Vec<u32>,
}

impl BinStatistics {
    fn of(counts: &[u32]) -> Self {
        let positions = positions_in(counts);
        if positions == 0 {
            return Self {
                positions: 0,
                mean: f64::NAN,
                median: f64::NAN,
                variance: f64::NAN,
                fourth_moment: f64::NAN,
                variance_share_of_deepest_hundredth: f64::NAN,
                counts: counts.to_vec(),
            };
        }
        let total: f64 = counts
            .iter()
            .enumerate()
            .map(|(depth, &count)| depth as f64 * f64::from(count))
            .sum();
        let mean = total / positions as f64;
        let (mut second, mut fourth) = (0.0, 0.0);
        for (depth, &count) in counts.iter().enumerate() {
            let deviation = depth as f64 - mean;
            second += deviation * deviation * f64::from(count);
            fourth += deviation.powi(4) * f64::from(count);
        }
        let variance = second / positions as f64;
        // Where the deepest hundredth of the positions begins, and what share of the total
        // squared deviation sits at or above it.
        let cut = quantile_of(counts, 0.99);
        let tail: f64 = counts
            .iter()
            .enumerate()
            .skip(cut as usize)
            .map(|(depth, &count)| {
                let deviation = depth as f64 - mean;
                deviation * deviation * f64::from(count)
            })
            .sum();
        Self {
            positions,
            mean,
            median: f64::from(quantile_of(counts, 0.5)),
            variance,
            fourth_moment: fourth / positions as f64,
            variance_share_of_deepest_hundredth: tail / second.max(1e-9),
            counts: counts.to_vec(),
        }
    }

    fn variance_to_mean(&self) -> f64 {
        self.variance / self.mean
    }

    /// The standard error of variance ÷ mean, from the bin's own fourth moment — so a bin whose
    /// ratio is noise says so rather than being assumed adequate.
    fn variance_to_mean_standard_error(&self) -> f64 {
        let n = self.positions as f64;
        let variance_of_variance = (self.fourth_moment - self.variance * self.variance) / n;
        variance_of_variance.max(0.0).sqrt() / self.mean
    }

    /// The negative-binomial size this bin implies on its own: `m² / (variance − m)`.
    fn implied_negative_binomial_size(&self) -> Option<f64> {
        (self.variance > self.mean).then(|| self.mean * self.mean / (self.variance - self.mean))
    }

    fn variance_share_of_deepest_hundredth(&self) -> f64 {
        self.variance_share_of_deepest_hundredth
    }

    /// The same share a Poisson and a negative binomial of this bin's own mean and variance
    /// would put in their own deepest hundredth. **A trimmed variance is not comparable across
    /// families without this**: trimming lowers a Poisson's variance too, so only the
    /// difference between observed and predicted says the tail is heavier than the family
    /// allows.
    fn variance_share_of_deepest_hundredth_under_models(&self) -> (f64, f64) {
        let poisson = model_tail_share(self.mean, |depth| poisson_log_pmf(depth, self.mean));
        let size = self.implied_negative_binomial_size();
        let negative_binomial = size.map_or(f64::NAN, |size| {
            model_tail_share(self.mean, |depth| {
                negative_binomial_log_pmf(depth, self.mean, size)
            })
        });
        (poisson, negative_binomial)
    }

    fn share_at_or_above(&self, threshold: u32) -> f64 {
        let above: u64 = self
            .counts
            .iter()
            .skip(threshold as usize)
            .map(|&count| u64::from(count))
            .sum();
        above as f64 / self.positions as f64
    }
}

/// Under a distribution given by `log_pmf`, what share of its variance sits in its deepest 1%.
fn model_tail_share<F: Fn(f64) -> f64>(mean: f64, log_pmf: F) -> f64 {
    let ceiling = (mean * 40.0).ceil() as u32 + 200;
    let probabilities: Vec<f64> = (0..=ceiling)
        .map(|depth| log_pmf(f64::from(depth)).exp())
        .collect();
    let total_variance: f64 = probabilities
        .iter()
        .enumerate()
        .map(|(depth, &p)| p * (depth as f64 - mean).powi(2))
        .sum();
    // Walk down from the top until 1% of the mass is accounted for.
    let mut mass = 0.0;
    let mut tail_variance = 0.0;
    for (depth, &p) in probabilities.iter().enumerate().rev() {
        if mass >= 0.01 {
            break;
        }
        let take = p.min(0.01 - mass);
        tail_variance += take * (depth as f64 - mean).powi(2);
        mass += take;
    }
    tail_variance / total_variance.max(1e-12)
}

/// The depth at or below which this share of the bin's positions sit.
fn quantile_of(counts: &[u32], share: f64) -> u32 {
    let total = positions_in(counts) as f64;
    let mut seen = 0.0;
    for (depth, &count) in counts.iter().enumerate() {
        seen += f64::from(count);
        if seen >= share * total {
            return depth as u32;
        }
    }
    counts.len() as u32 - 1
}

/// Which of the two overdispersed families the GC bins agree with, and its one parameter.
struct FamilyFit {
    lowest_mean: f64,
    highest_mean: f64,
    lowest_ratio: f64,
    highest_ratio: f64,
    /// The one size a negative binomial would use for the whole sample: `variance = m + m²/r`.
    negative_binomial_size: f64,
    /// How far that leaves each bin's variance ÷ mean, averaged over bins and weighted by how
    /// many positions each bin holds.
    negative_binomial_residual: f64,
    /// The one multiplier a quasi-Poisson would use: `variance = φ·m` at every depth.
    quasi_poisson_multiplier: f64,
    quasi_poisson_residual: f64,
}

impl FamilyFit {
    fn report(&self, over: &str, sample_name: &str, mean_depth: f64) {
        println!("\n  over {over}:");
        println!(
            "    mean depth across the surviving bins runs {:.2} to {:.2} reads a position, a \
             factor of {:.2}",
            self.lowest_mean,
            self.highest_mean,
            self.highest_mean / self.lowest_mean.max(1e-9)
        );
        println!(
            "    variance ÷ mean over the same bins runs {:.2} to {:.2}",
            self.lowest_ratio, self.highest_ratio
        );
        println!(
            "    negative binomial, one size for the sample:  r = {:.2}  — misses each bin's \
             ratio by {:.3} on average",
            self.negative_binomial_size, self.negative_binomial_residual
        );
        println!(
            "    quasi-Poisson, one multiplier:               φ = {:.2}  — misses each bin's \
             ratio by {:.3} on average",
            self.quasi_poisson_multiplier, self.quasi_poisson_residual
        );
        println!(
            "    the family whose shape fits these bins: {}",
            if self.negative_binomial_residual < self.quasi_poisson_residual {
                "negative binomial"
            } else {
                "quasi-Poisson"
            }
        );
        println!(
            "TSV\tfamily\t{sample_name}\t{mean_depth:.4}\t{over}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t\
             {:.4}\t{:.4}\t{:.4}\t{:.4}",
            self.lowest_mean,
            self.highest_mean,
            self.lowest_ratio,
            self.highest_ratio,
            self.negative_binomial_size,
            self.negative_binomial_residual,
            self.quasi_poisson_multiplier,
            self.quasi_poisson_residual,
        );
    }

    fn over(bins: &[usize], stats: &HashMap<usize, BinStatistics>) -> Self {
        let usable: Vec<&BinStatistics> = bins.iter().map(|bin| &stats[bin]).collect();
        // Weighted least squares of (variance − mean) against mean², the negative binomial's
        // shape; and of variance against mean, the quasi-Poisson's. Weight by positions.
        let (mut numerator, mut denominator) = (0.0, 0.0);
        let (mut variance_sum, mut mean_sum) = (0.0, 0.0);
        for s in &usable {
            let weight = s.positions as f64;
            numerator += weight * s.mean.powi(4);
            denominator += weight * s.mean.powi(2) * (s.variance - s.mean);
            variance_sum += weight * s.variance;
            mean_sum += weight * s.mean;
        }
        // A sample whose bins are not overdispersed at all leaves nothing for the size to fit.
        // The negative binomial's answer there is *a Poisson*, and a size this large is one to
        // every digit that prints, so the arms below degrade to the Poisson arm rather than to
        // an infinity that makes every score `NaN`.
        const POISSON_SIZE: f64 = 1e9;
        let size = if denominator > 0.0 {
            (numerator / denominator).min(POISSON_SIZE)
        } else {
            POISSON_SIZE
        };
        let multiplier = variance_sum / mean_sum.max(1e-9);
        let residual = |predict: &dyn Fn(f64) -> f64| {
            let mut weighted = 0.0;
            let mut weight_total = 0.0;
            for s in &usable {
                let weight = s.positions as f64;
                weighted += weight * (s.variance_to_mean() - predict(s.mean)).abs();
                weight_total += weight;
            }
            weighted / weight_total.max(1e-9)
        };
        Self {
            lowest_mean: usable.iter().map(|s| s.mean).fold(f64::INFINITY, f64::min),
            highest_mean: usable
                .iter()
                .map(|s| s.mean)
                .fold(f64::NEG_INFINITY, f64::max),
            lowest_ratio: usable
                .iter()
                .map(|s| s.variance_to_mean())
                .fold(f64::INFINITY, f64::min),
            highest_ratio: usable
                .iter()
                .map(|s| s.variance_to_mean())
                .fold(f64::NEG_INFINITY, f64::max),
            negative_binomial_size: size,
            negative_binomial_residual: residual(&|mean| 1.0 + mean / size),
            quasi_poisson_multiplier: multiplier,
            quasi_poisson_residual: residual(&|_| multiplier),
        }
    }
}

fn gc_key(contig: u32, window: u64) -> u64 {
    (u64::from(contig) << 40) | window
}

fn gc_bin(fraction: f64, bins: usize) -> usize {
    ((fraction * (bins - 1) as f64).round() as usize).min(bins - 1)
}

fn gc_fraction_of((gc, at): (u32, u32)) -> Option<f64> {
    let called = gc + at;
    (called > 0).then(|| f64::from(gc) / f64::from(called))
}
