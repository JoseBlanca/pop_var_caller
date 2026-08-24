//! **How wide the allele table is, and what a bar and a cap would do to it** — the
//! measurement candidate selection (step 6) has no design without.
//!
//! The cohort merge unifies every sequence any sample showed into one allele table
//! (`spec/cohort_merge.md` §4.2) and narrows nothing. Selection has to decide what survives,
//! and two numbers that nobody has measured decide whether the decision matters at all: **how
//! many alternatives a real locus carries**, and **how many of them a support bar removes**.
//! This walks a real cohort's reads, assembles the loci the merge would assemble, and prints
//! both distributions.
//!
//! ```text
//! ng_candidate_selection_probe <reference.fa> <cram-dir> <regions.bed>
//! ```
//!
//! `NG_REAL_SAMPLES=n` walks only the first `n` alignment files of the directory in name order;
//! `NG_REAL_REGIONS=n` only the first `n` intervals of the BED. Both default to everything, and
//! both cost memory: every sample's observations over every interval are held at once, which is
//! what the merge consumes.
//!
//! **The bar asked here is the merge's own keep rule, one level down.** The merge asks each
//! sample whether its *non-reference reads* reach `max(floor, ceil(share × its compared
//! reads))` and builds the locus if any one does (`MinAltReads`,
//! `ng::run::cohort_merge`). This asks the identical question of each *alternative allele*
//! separately: an alternative survives if some single sample lent it that many reads. Nothing
//! new is introduced — the floor and the share are the merge's, so a sweep here and a sweep
//! there move the same two knobs.
//!
//! **The cap is applied after the bar and ranks by the largest within-sample share** — the
//! allele's reads in one sample divided by that sample's compared reads, maximised over
//! samples. Production ranks by the cohort's raw read total
//! (`var_calling::per_group_merger::enforce_max_alleles`), which at large cohorts truncates
//! away the private alleles first; this probe prints how often the two rankings disagree, so
//! the difference is a number rather than an argument.
//!
//! **Every allele the bar or the cap removes is counted with its summed error mass**, because
//! the SNP/indel read likelihood needs that pool and nothing upstream produces it
//! (`spec/read_likelihoods.md` §3.3's `q_sum_other`). Printing it here shows what it is made of
//! before anything consumes it.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, SampleLocusObservations, SampleLocusObservationsIterator,
    UnhandledReason,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{RegionKind, TypedRegion};
use pop_var_caller::ng::run::cohort_merge::build::CohortObservation;
use pop_var_caller::ng::run::cohort_merge::close::{LocusCloser, Verdict};
use pop_var_caller::ng::run::cohort_merge::{
    MaxCohortLocusSpan, MinAltObs, MinAltReadShare, MinAltReads,
};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

#[path = "shared/reference_check.rs"]
mod reference_check_knob;
use reference_check_knob::reference_check_from_env;

/// Allele index 0 is the reference, by the merge's own contract.
const REFERENCE_ALLELE: usize = 0;

/// The caps swept. 6 is production's `DEFAULT_MAX_ALLELES_PER_RECORD` and GATK's
/// `--max-alternate-alleles` default; the rest bracket it.
const CAPS: [usize; 4] = [3, 4, 6, 8];

/// What one setting of the bar did to one cohort's loci.
#[derive(Default)]
struct Tally {
    loci: u64,
    /// Alternatives in the table the merge handed over, and after the bar.
    alts_before: u64,
    alts_after: u64,
    /// Loci whose alternatives all failed the bar — nothing left to call.
    emptied: u64,
    /// How many loci carry each alternative count, before and after (index = count, last
    /// bucket is "that many or more").
    hist_before: [u64; 12],
    hist_after: [u64; 12],
    /// Reads and summed error mass the bar sent to the pool.
    dropped_reads: u64,
    dropped_q_sum: f64,
    /// Reads that stayed on a surviving allele — the pool's denominator.
    kept_reads: u64,
    /// Loci over each cap, after the bar, and the alternatives each cap would truncate.
    over_cap: [u64; CAPS.len()],
    truncated: [u64; CAPS.len()],
    /// Loci where ranking by largest within-sample share and by cohort read total keep
    /// different alleles at a cap of 6.
    ranking_disagrees: u64,
}

impl Tally {
    fn bucket(hist: &mut [u64; 12], count: usize) {
        hist[count.min(11)] += 1;
    }
}

/// One alternative allele's support, folded across the cohort — what a bar and a cap read.
struct AlleleSupportSummary {
    allele: usize,
    /// The largest share of one sample's compared reads this allele took.
    best_within_sample_share: f64,
    /// Reads over the whole cohort — production's ranking key.
    cohort_reads: u64,
    /// Samples whose reads reached the bar.
    samples_clearing: u32,
    /// Reads and summed error mass, over every sample, for the pool.
    reads: u64,
    q_sum: f64,
    /// Whether some single sample's reads reached the bar.
    passes: bool,
}

/// Fold one locus's table into per-allele summaries, asking the bar of each sample.
fn summarise(observation: &CohortObservation, bar: MinAltReads) -> Vec<AlleleSupportSummary> {
    summarise_over(observation, bar, usize::MAX)
}

/// The same fold, asking only the first `samples` samples of the run — how many alternatives
/// clear the bar at a fixed allele table as the cohort grows.
fn summarise_over(
    observation: &CohortObservation,
    bar: MinAltReads,
    samples: usize,
) -> Vec<AlleleSupportSummary> {
    let n_alleles = observation.alleles.len();
    let mut out: Vec<AlleleSupportSummary> = (0..n_alleles)
        .map(|allele| AlleleSupportSummary {
            allele,
            best_within_sample_share: 0.0,
            cohort_reads: 0,
            samples_clearing: 0,
            reads: 0,
            q_sum: 0.0,
            passes: false,
        })
        .collect();
    for sample in observation.per_sample.iter().filter(|s| s.sample < samples) {
        // A sample's compared reads at this locus are its reads across the whole table: the
        // merge admits only reads that spanned the locus, and every one lands on some allele.
        let compared: u32 = sample
            .supported
            .iter()
            .map(|supported| supported.support.num_reads)
            .sum();
        for supported in &sample.supported {
            let entry = &mut out[supported.allele];
            let reads = supported.support.num_reads;
            entry.cohort_reads += u64::from(reads);
            entry.reads += u64::from(reads);
            entry.q_sum += supported.support.q_sum;
            if compared > 0 {
                let share = f64::from(reads) / f64::from(compared);
                if share > entry.best_within_sample_share {
                    entry.best_within_sample_share = share;
                }
            }
            if bar.reached_by(reads, compared) {
                entry.samples_clearing += 1;
                entry.passes = true;
            }
        }
    }
    out[REFERENCE_ALLELE].passes = true; // the reference is always present
    out
}

/// The alleles a cap keeps, ranked by the largest within-sample share, then by how many
/// samples cleared the bar, then by cohort reads. The reference is never prunable.
fn keep_by_share(survivors: &[&AlleleSupportSummary], cap: usize) -> Vec<usize> {
    let mut ranked: Vec<&&AlleleSupportSummary> = survivors
        .iter()
        .filter(|s| s.allele != REFERENCE_ALLELE)
        .collect();
    ranked.sort_by(|a, b| {
        b.best_within_sample_share
            .total_cmp(&a.best_within_sample_share)
            .then(b.samples_clearing.cmp(&a.samples_clearing))
            .then(b.cohort_reads.cmp(&a.cohort_reads))
            .then(a.allele.cmp(&b.allele))
    });
    ranked.truncate(cap.saturating_sub(1));
    let mut kept: Vec<usize> = ranked.iter().map(|s| s.allele).collect();
    kept.sort_unstable();
    kept
}

/// The same cap, ranked production's way: the cohort's raw read total.
fn keep_by_cohort_reads(survivors: &[&AlleleSupportSummary], cap: usize) -> Vec<usize> {
    let mut ranked: Vec<&&AlleleSupportSummary> = survivors
        .iter()
        .filter(|s| s.allele != REFERENCE_ALLELE)
        .collect();
    ranked.sort_by(|a, b| {
        b.cohort_reads
            .cmp(&a.cohort_reads)
            .then(a.allele.cmp(&b.allele))
    });
    ranked.truncate(cap.saturating_sub(1));
    let mut kept: Vec<usize> = ranked.iter().map(|s| s.allele).collect();
    kept.sort_unstable();
    kept
}

fn fold_locus(tally: &mut Tally, observation: &CohortObservation, bar: MinAltReads) {
    let summaries = summarise(observation, bar);
    let alts_before = observation.alleles.len().saturating_sub(1);
    tally.loci += 1;
    tally.alts_before += alts_before as u64;
    Tally::bucket(&mut tally.hist_before, alts_before);

    let survivors: Vec<&AlleleSupportSummary> =
        summaries.iter().filter(|entry| entry.passes).collect();
    let alts_after = survivors.len().saturating_sub(1);
    tally.alts_after += alts_after as u64;
    Tally::bucket(&mut tally.hist_after, alts_after);
    if alts_after == 0 && alts_before > 0 {
        tally.emptied += 1;
    }

    for entry in &summaries {
        if entry.passes {
            tally.kept_reads += entry.reads;
        } else {
            tally.dropped_reads += entry.reads;
            tally.dropped_q_sum += entry.q_sum;
        }
    }

    for (at, cap) in CAPS.iter().enumerate() {
        if survivors.len() > *cap {
            tally.over_cap[at] += 1;
            tally.truncated[at] += (survivors.len() - cap) as u64;
        }
    }
    if survivors.len() > 6 && keep_by_share(&survivors, 6) != keep_by_cohort_reads(&survivors, 6) {
        tally.ranking_disagrees += 1;
    }
}

fn report(label: &str, bar: MinAltReads, tally: &Tally) {
    println!(
        "\n## bar: {} reads or {:.0}% of a sample's compared reads, whichever is more — {label}",
        bar.floor.get(),
        100.0 * bar.share.get()
    );
    println!("built loci: {}", tally.loci);
    let n = tally.loci.max(1) as f64;
    println!(
        "alternatives per locus: {:.2} from the merge, {:.2} after the bar ({:.1}% removed)",
        tally.alts_before as f64 / n,
        tally.alts_after as f64 / n,
        100.0 * (tally.alts_before.saturating_sub(tally.alts_after)) as f64
            / tally.alts_before.max(1) as f64,
    );
    println!(
        "loci left with no alternative at all: {} ({:.1}%)",
        tally.emptied,
        100.0 * tally.emptied as f64 / n
    );
    print!("alternatives per locus, before the bar:");
    for (count, loci) in tally.hist_before.iter().enumerate() {
        if *loci > 0 {
            print!(" {count}:{loci}");
        }
    }
    println!();
    print!("alternatives per locus, after the bar: ");
    for (count, loci) in tally.hist_after.iter().enumerate() {
        if *loci > 0 {
            print!(" {count}:{loci}");
        }
    }
    println!();
    println!(
        "reads the bar sent to the pool: {} of {} ({:.2}%), summed error mass {:.0} nats",
        tally.dropped_reads,
        tally.dropped_reads + tally.kept_reads,
        100.0 * tally.dropped_reads as f64 / (tally.dropped_reads + tally.kept_reads).max(1) as f64,
        tally.dropped_q_sum,
    );
    for (at, cap) in CAPS.iter().enumerate() {
        println!(
            "cap {cap}: {} loci over it ({:.3}%), {} alternatives truncated",
            tally.over_cap[at],
            100.0 * tally.over_cap[at] as f64 / n,
            tally.truncated[at],
        );
    }
    println!(
        "loci where the two rankings keep different alleles at a cap of 6: {}",
        tally.ranking_disagrees
    );
}

/// The BED's intervals as the analysed regions, in the reference's contig order.
/// One-based inclusive, where a BED is zero-based half-open.
fn analysed_regions_of(
    bed: &Path,
    contig_index: impl Fn(&str) -> Option<u32>,
    limit: Option<usize>,
) -> Result<Vec<GenomeRegion>, Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(bed)?;
    let mut regions = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let (Some(contig), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
        else {
            return Err(format!("a BED line has fewer than three fields: {line}").into());
        };
        let Some(contig) = contig_index(contig) else {
            return Err(format!("the reference has no contig named {contig}").into());
        };
        regions.push(GenomeRegion {
            contig: ContigId(contig),
            start: Position(start.parse::<u64>()? + 1),
            end: Position(end.parse::<u64>()?),
        });
    }
    regions.sort_by_key(|region| (region.contig.0, region.start.0));
    if let Some(limit) = limit {
        regions.truncate(limit);
    }
    Ok(regions)
}

/// Walk one sample's reads over `analysed` and keep every observation, in coordinate order —
/// `ng_cohort_merge_real_cost`'s pipeline, unchanged.
fn walk_one_sample(
    fasta: &Path,
    cram: &Path,
    analysed: &[GenomeRegion],
    cache: &Arc<ReferenceInfoCache>,
) -> Result<Vec<SampleLocusObservations>, Box<dyn std::error::Error>> {
    let check = reference_check_from_env()?;
    let (info, _verify) =
        read_reference_verifying_or_creating_fai(cache, fasta.to_path_buf(), check)?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(fasta)?;
    let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
        fasta.to_path_buf(),
        contigs.clone(),
        index.clone(),
    ));

    let reference = OpenReference::new(info);
    let reads = SampleReads::open_only_sample(
        &[cram.to_path_buf()],
        &reference,
        ReadFilterConfig::default(),
        true,
    )?;

    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "PileupGenerator::new takes Arc and this accessor is file-backed and single-threaded, as in ng_cohort_merge_real_cost"
    )]
    let shared = Arc::new(WindowedRefSeq::with_shared_index(
        fasta.to_path_buf(),
        contigs.clone(),
        index.clone(),
    ));
    let make_reference = {
        let fasta = fasta.to_path_buf();
        let contigs = contigs.clone();
        let index = index.clone();
        move || WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone())
    };
    let generator = PileupGenerator::new(
        shared,
        make_reference,
        preparer,
        PileupGeneratorConfig::default(),
    )?;
    let generators = GeneratorSet::new(
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        GeneratorSlot::Generator(Box::new(generator)),
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
    );

    let regions: Vec<Result<TypedRegion, pop_var_caller::ng::repeat_catalog::RepeatCatalogError>> =
        analysed
            .iter()
            .map(|region| {
                Ok(TypedRegion {
                    region: *region,
                    kind: RegionKind::Generic,
                })
            })
            .collect();

    let mut observations = Vec::new();
    let mut stream = SampleLocusObservationsIterator::new(regions.into_iter(), reads, generators);
    for locus in &mut stream {
        observations.push(locus?);
    }
    Ok(observations)
}

fn run(fasta: &Path, crams: &Path, bed: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let limit_of = |name: &str| -> Option<usize> {
        std::env::var(name)
            .ok()
            .map(|value| value.parse().expect("a count"))
    };

    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, _verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        reference_check_from_env()?,
    )?;
    let contigs = info.contig_list();
    let analysed = analysed_regions_of(
        bed,
        |name| {
            contigs
                .entries
                .iter()
                .position(|entry| entry.name == name)
                .map(|at| at as u32)
        },
        limit_of("NG_REAL_REGIONS"),
    )?;

    let mut cram_paths: Vec<PathBuf> = std::fs::read_dir(crams)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| {
            path.extension()
                .is_some_and(|kind| kind == "cram" || kind == "bam")
        })
        .collect();
    cram_paths.sort();
    if let Some(limit) = limit_of("NG_REAL_SAMPLES") {
        cram_paths.truncate(limit);
    }
    if cram_paths.is_empty() {
        return Err(format!("no .cram or .bam under {}", crams.display()).into());
    }

    println!("# samples: {}", cram_paths.len());
    println!("# analysed intervals: {}", analysed.len());
    println!(
        "# analysed bases: {}",
        analysed.iter().map(|region| region.len()).sum::<u64>()
    );

    let mut cohort: Vec<Vec<SampleLocusObservations>> = Vec::with_capacity(cram_paths.len());
    for cram in &cram_paths {
        cohort.push(walk_one_sample(fasta, cram, &analysed, &cache)?);
    }

    // The merge's keep rule stays at its default here: this probe sweeps the *allele* bar,
    // and moving both at once would leave nothing to attribute a change to.
    let keep_rule = MinAltReads::DEFAULT;
    let all: Vec<&[SampleLocusObservations]> = cohort.iter().map(Vec::as_slice).collect();
    let built: Vec<CohortObservation> =
        LocusCloser::over(&all, MaxCohortLocusSpan::DEFAULT, keep_rule)
            .filter(|locus| locus.verdict == Verdict::Build)
            .map(|locus| CohortObservation::over(&locus))
            .collect();
    println!("# loci the merge built: {}", built.len());

    let bars: Vec<(String, MinAltReads)> =
        [(2u32, 0.0f64), (2, 0.02), (2, 0.05), (2, 0.10), (3, 0.02)]
            .iter()
            .map(|(floor, share)| {
                (
                    format!("floor {floor}, share {share}"),
                    MinAltReads {
                        floor: MinAltObs(
                            std::num::NonZeroU32::new(*floor).expect("non-zero floor"),
                        ),
                        share: MinAltReadShare::new(*share).expect("a fraction of one"),
                    },
                )
            })
            .collect();

    for (label, bar) in &bars {
        let mut tally = Tally::default();
        for observation in &built {
            fold_locus(&mut tally, observation, *bar);
        }
        report(label, *bar, &tally);
    }

    // **What the bar admitted, allele by allele, for scoring against a truth set.**
    // `NG_SELECT_DUMP=<path>` writes one row per (locus, allele): the locus span, the
    // reference bases over it and the allele's own, and whether the bar passed it. Joining
    // that against a benchmark's truth VCF — projecting each truth record onto the same span,
    // which is what the merge itself does (`spec/cohort_merge.md` §4.2) — turns "the bar
    // removed 64% of the table" into "the bar lost this many true variants".
    if let Ok(path) = std::env::var("NG_SELECT_DUMP") {
        use std::io::Write;
        let mut out = std::io::BufWriter::new(std::fs::File::create(&path)?);
        writeln!(
            out,
            "contig\tstart\tend\tref\talt\tpassed\tbest_share\tcohort_reads"
        )?;
        // `NG_SELECT_DUMP_FLOOR` / `NG_SELECT_DUMP_SHARE` set the bar the dump's `passed`
        // column reports, so the recall curve is one run per point rather than one guess.
        let bar = MinAltReads {
            floor: match std::env::var("NG_SELECT_DUMP_FLOOR").ok() {
                Some(raw) => MinAltObs(
                    std::num::NonZeroU32::new(raw.parse().expect("a number")).expect("non-zero"),
                ),
                None => MinAltObs::DEFAULT,
            },
            share: match std::env::var("NG_SELECT_DUMP_SHARE").ok() {
                Some(raw) => {
                    MinAltReadShare::new(raw.parse().expect("a number")).expect("a fraction of one")
                }
                None => MinAltReadShare::DEFAULT,
            },
        };
        for observation in &built {
            let summaries = summarise(observation, bar);
            let reference =
                String::from_utf8_lossy(&observation.alleles[REFERENCE_ALLELE]).into_owned();
            for entry in summaries.iter().skip(1) {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
                    observation.region.contig.0,
                    observation.region.start.get(),
                    observation.region.end.get(),
                    reference,
                    String::from_utf8_lossy(&observation.alleles[entry.allele]),
                    u8::from(entry.passes),
                    entry.best_within_sample_share,
                    entry.cohort_reads,
                )?;
            }
        }
        out.flush()?;
        println!("# allele dump written to {path}");
    }

    // **Does the admitted allele count grow with the cohort?** The allele table is held
    // fixed — the whole run's — and only the number of samples the bar is asked of varies.
    // A per-sample bar asks the same question however many samples there are, but more
    // samples are more independent chances for one to reach it, so what is admitted can
    // still grow. This says by how much, which decides whether the cap is a safety valve
    // or a working part at large cohorts.
    let bar = MinAltReads::DEFAULT;
    println!(
        "\n## alternatives admitted as the cohort grows (allele table fixed, bar {} reads or {:.0}%)",
        bar.floor.get(),
        100.0 * bar.share.get()
    );
    println!(
        "{:>9}{:>14}{:>16}{:>16}{:>14}",
        "samples", "alts/locus", "loci with none", "loci over 6", "max alts"
    );
    let n_samples = cram_paths.len();
    let mut sizes: Vec<usize> = vec![1];
    while *sizes.last().expect("seeded") * 4 < n_samples {
        sizes.push(sizes.last().expect("seeded") * 4);
    }
    sizes.push(n_samples);
    for size in sizes {
        let (mut alts, mut none, mut over, mut max) = (0u64, 0u64, 0u64, 0usize);
        for observation in &built {
            let summaries = summarise_over(observation, bar, size);
            let passing = summaries.iter().filter(|e| e.passes).count();
            let a = passing.saturating_sub(1);
            alts += a as u64;
            if a == 0 {
                none += 1;
            }
            if passing > 6 {
                over += 1;
            }
            max = max.max(a);
        }
        let n = built.len().max(1) as f64;
        println!(
            "{:>9}{:>14.3}{:>13} ({:4.1}%){:>11} ({:5.3}%){:>14}",
            size,
            alts as f64 / n,
            none,
            100.0 * none as f64 / n,
            over,
            100.0 * over as f64 / n,
            max
        );
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: ng_candidate_selection_probe <reference.fa> <cram-dir> <regions.bed>");
        return ExitCode::from(2);
    }
    match run(
        Path::new(&args[1]),
        Path::new(&args[2]),
        Path::new(&args[3]),
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
