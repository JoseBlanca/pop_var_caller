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
//! **This runs the shipped rule rather than a copy of it.** Every allele question here is
//! answered by [`select_generic`] — which alternatives clear the bar, which the cap keeps, and
//! how many reads and how much error mass the dropped alleles carry — so a number printed here
//! is a number the caller itself would produce. **One rule stays the probe's own, deliberately:**
//! production's ranking by cohort read total, which is the rival the two-ranking comparison below
//! measures the shipped ranking against.
//!
//! **The bar is the merge's own keep rule, one level down.** The merge asks each sample whether
//! its *non-reference reads* reach `max(floor, ceil(share × its compared reads))` and builds the
//! locus if any one does (`MinAltReads`, `ng::run::cohort_merge`). Selection asks the identical
//! question of each *alternative allele* separately: an alternative survives if some single
//! sample lent it that many reads. Nothing new is introduced — the floor and the share are the
//! merge's, so a sweep here and a sweep there move the same two knobs.
//!
//! **The cap is applied after the bar and ranks by the largest within-sample share** — the
//! allele's reads in one sample divided by that sample's compared reads, maximised over the
//! samples that cleared the bar for it (`spec/candidate_alleles.md` §4.1). Production ranks by
//! the cohort's raw read total (`var_calling::per_group_merger::enforce_max_alleles`), which at
//! large cohorts truncates away the private alleles first; this probe prints how often the two
//! rankings keep different alleles, so the difference is a number rather than an argument.
//!
//! **Three rules differ between that copy and the shipped module, and all three differences are
//! deliberate.** The figures `spec/candidate_alleles.md` §3.3, §4.2 and §5 quote were taken with
//! the copy, so each of the three is a reason a number here could legitimately differ from them.
//!
//! 1. **The bar was asked of each `(allele, read group)` row separately; the module pools a
//!    sample's rows first.** The copy therefore applied a stricter rule to exactly the samples
//!    carrying more than one library — 157 of 1,707 in a surveyed tomato archive
//!    (`spec/read_groups.md` §1) — so its figures are a lower bound on what the module admits.
//!    **Measured: it moves nothing on either benchmark.** Every bar total, every cap count and
//!    the leftover are identical under the two rules on the 63-accession panel and on the trio.
//! 2. **The within-sample share was maximised over every sample; the module maximises over the
//!    samples that cleared the bar** (`spec/candidate_alleles.md` §4.1, the owner's decision of
//!    2026-08-24). **Measured: this is the one difference the tomato panel exercises**, and it
//!    is why the two cap rankings disagree at 19 loci where the copy found 17.
//! 3. **The cap's last tie-break was the merge table's index; the module's is the allele's
//!    bases** (`compare_best_first`). The merge interns alleles in first-seen order rather than
//!    byte order, so these are genuinely different orders. **Measured: it cannot decide anything
//!    on this panel.** The bases only speak when all three numeric keys tie, and at every
//!    cap-binding tomato locus, at all five bars swept, no two surviving alternatives even share
//!    a within-sample share and a cohort read total — so the tie-break is never reached and it
//!    contributes none of the 17-to-19 move.
//!
//! **Any difference beyond these three is a defect in one of the two implementations, and is to
//! be traced rather than accepted.**
//!
//! **The second difference also empties a column of the dump, which is a change in what the
//! column means rather than in what the bar did** — so the column is renamed rather than left to
//! change under its old name. `best_share_of_clearing_samples` is 0 for every allele the bar
//! rejected, two rows in every three of the tomato dump (73,649 of 115,329), and it is smaller
//! than the old `best_share` on admitted rows too, wherever the old maximum was lent by a sample
//! that never cleared the bar.
//!
//! **Every allele the bar or the cap removes is counted with its summed error mass**, because
//! the SNP/indel read likelihood needs that pool and nothing upstream produces it
//! (`spec/read_likelihoods.md` §3.3's `q_sum_other`). It is read straight off the module's
//! per-sample leftover, which is the value the calling loop will consume.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::ng::calling::allele_candidates::generic::select_generic;
use pop_var_caller::ng::calling::allele_candidates::{
    CandidateSelectionConfig, DEFAULT_MAX_CANDIDATE_ALLELES, LocusSelection, MaxCandidateAlleles,
    SelectionScratch, SelectionVerdict,
};
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

/// The caps swept, counting the reference. 6 is production's
/// `DEFAULT_MAX_ALLELES_PER_RECORD` ([`PRODUCTION_CAP`]) and the rest bracket it.
///
/// **It is not GATK's number, though production's comment says it is**: GATK's
/// `--max-alternate-alleles` defaults to 6 *alternates*, so seven alleles where these are six
/// (`DEFAULT_MAX_CANDIDATE_ALLELES` sets that out).
const CAPS: [usize; 4] = [3, 4, 6, 8];

/// **The module's shipping default, which is production's `DEFAULT_MAX_ALLELES_PER_RECORD`** —
/// the cap the two rankings are compared at, and the one the cohort-growth table reports loci
/// above, so both questions are asked where a run would be asking them.
///
/// **Read from the module rather than written down**, because a step whose whole point is to
/// stop copying the module's rules should not copy its number either: change the default and
/// this run's two comparisons follow it.
const PRODUCTION_CAP: usize = DEFAULT_MAX_CANDIDATE_ALLELES.get() as usize;

/// **A cap that cannot bind**, so that a selection run under it reports the bar alone.
///
/// Not an `Option`: `MaxCandidateAlleles` refuses anything below two precisely because a
/// cap is never absent, and the widest one the type can hold is 65,534 alternatives, where
/// the widest allele table either benchmark builds carries 102 — one tomato locus; the
/// trio reaches 19 at 300× and 4 at 30×.
const WIDEST_CAP: MaxCandidateAlleles = MaxCandidateAlleles::new_or_panic(u16::MAX);

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
    /// different alleles at [`PRODUCTION_CAP`].
    ranking_disagrees: u64,
}

impl Tally {
    fn bucket(hist: &mut [u64; 12], count: usize) {
        hist[count.min(11)] += 1;
    }
}

/// Narrow one locus with the shipped selection, under a given bar and cap.
///
/// **The scratch is handed back full**: it holds this locus's per-allele fold until the next
/// call, so [`SelectionScratch::best_within_sample_share_of`] and
/// [`SelectionScratch::cohort_reads_of`] answer for this locus until then. That is where the
/// ranking keys this probe reports come from.
fn select_under(
    observation: &CohortObservation,
    bar: MinAltReads,
    cap: MaxCandidateAlleles,
    scratch: &mut SelectionScratch,
) -> LocusSelection {
    let config = CandidateSelectionConfig {
        min_allele_support: bar,
        max_candidate_alleles: cap,
    };
    select_generic(observation, &config, scratch)
}

/// The cap as the selection config wants it. Every value in [`CAPS`] is at least two, which
/// is what `MaxCandidateAlleles` refuses to go below.
fn cap_of(alleles: usize) -> MaxCandidateAlleles {
    let alleles = u16::try_from(alleles).expect("a cap that fits a u16");
    MaxCandidateAlleles::new(alleles).expect("a cap of at least two alleles")
}

/// **The alternatives a selection kept, as indices into the merge's table**, ascending —
/// the same shape [`keep_by_cohort_reads`] returns, so the two can be compared directly.
fn kept_alternative_indices(selection: &LocusSelection, table_len: usize) -> Vec<usize> {
    (1..table_len)
        .filter(|&index| selection.remap().candidate_for(index).is_some())
        .collect()
}

/// **The cap ranked production's way: the cohort's raw read total**, over the alternatives
/// that already cleared the bar. Each entry is a merge table index and that allele's reads
/// across every covering sample.
///
/// **This is the one rule the probe still owns**, because it is production's and not this
/// module's — it is what the shipped ranking is being compared against
/// (`enforce_max_alleles`, `src/var_calling/per_group_merger.rs`, a stable sort on
/// `Reverse(cohort_count)`).
fn keep_by_cohort_reads(survivors: &[(usize, u64)], cap: usize) -> Vec<usize> {
    let mut ranked = survivors.to_vec();
    ranked.sort_by(|(left_allele, left_reads), (right_allele, right_reads)| {
        right_reads
            .cmp(left_reads)
            .then(left_allele.cmp(right_allele))
    });
    ranked.truncate(cap.saturating_sub(1));
    let mut kept: Vec<usize> = ranked.into_iter().map(|(allele, _)| allele).collect();
    kept.sort_unstable();
    kept
}

/// **Every read the merge attributed to an allele at this locus**, over every sample and
/// every read group — the denominator the pool is reported as a share of.
fn table_reads_of(observation: &CohortObservation) -> u64 {
    observation
        .per_sample
        .iter()
        .flat_map(|sample| sample.supported.iter())
        .map(|row| u64::from(row.support.num_reads))
        .sum()
}

/// **The same locus with only the first `samples` samples of the run covering it**, its
/// allele table untouched — how the bar's answer moves as a cohort grows, with the table
/// held fixed so that nothing but the number of samples asked can move it.
fn restrict_to_first_samples(observation: &CohortObservation, samples: usize) -> CohortObservation {
    CohortObservation {
        region: observation.region,
        alleles: observation.alleles.clone(),
        per_sample: observation
            .per_sample
            .iter()
            .filter(|sample| sample.sample < samples)
            .cloned()
            .collect(),
    }
}

fn fold_locus(
    tally: &mut Tally,
    observation: &CohortObservation,
    bar: MinAltReads,
    scratch: &mut SelectionScratch,
) {
    let table_len = observation.alleles.len();
    let alts_before = table_len.saturating_sub(1);
    tally.loci += 1;
    tally.alts_before += alts_before as u64;
    Tally::bucket(&mut tally.hist_before, alts_before);

    // The bar on its own, under a cap that cannot bind: what survives here is exactly what
    // some single sample's reads earned.
    let admitted = select_under(observation, bar, WIDEST_CAP, scratch);
    let alts_after = admitted.alternative_allele_count();
    tally.alts_after += alts_after as u64;
    Tally::bucket(&mut tally.hist_after, alts_after);
    if alts_after == 0 && alts_before > 0 {
        tally.emptied += 1;
    }

    // The pool, read off the per-sample leftover the module produced — the same value the
    // read likelihood will consume, rather than a second sum over the dropped alleles.
    let mut dropped_reads = 0_u64;
    for leftover in admitted.unmatched() {
        dropped_reads += u64::from(leftover.num_reads);
        tally.dropped_q_sum += leftover.q_sum;
    }
    tally.dropped_reads += dropped_reads;
    tally.kept_reads += table_reads_of(observation).saturating_sub(dropped_reads);

    // The survivors and their cohort read totals, taken from the fold before the capped runs
    // below overwrite it. Only the rival ranking needs them; the shipped one is asked by
    // running the cap itself.
    let survivors: Vec<(usize, u64)> = kept_alternative_indices(&admitted, table_len)
        .into_iter()
        .map(|allele| (allele, scratch.cohort_reads_of(allele)))
        .collect();

    for (at, cap) in CAPS.iter().enumerate() {
        let capped = select_under(observation, bar, cap_of(*cap), scratch);
        // `SelectionVerdict` is `#[non_exhaustive]` and the repeat-tract path adds a variant,
        // so this asks for the one it wants rather than matching exhaustively.
        let SelectionVerdict::Truncated { dropped } = capped.verdict() else {
            continue;
        };
        tally.over_cap[at] += 1;
        tally.truncated[at] += u64::from(dropped);
        if *cap == PRODUCTION_CAP
            && kept_alternative_indices(&capped, table_len)
                != keep_by_cohort_reads(&survivors, *cap)
        {
            tally.ranking_disagrees += 1;
        }
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
        "loci where the two rankings keep different alleles at a cap of {PRODUCTION_CAP}: {}",
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

    // One set of buffers for the whole run, as a worker would hold it.
    let mut scratch = SelectionScratch::new();

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
            fold_locus(&mut tally, observation, *bar, &mut scratch);
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
        // **The seventh column is renamed rather than left to change meaning under its old
        // name**, which is the one surface change this step makes deliberately. It used to be
        // the largest within-sample share over *every* covering sample; the shipped fold
        // raises the share only inside the bar's own branch, so it is now a maximum over the
        // samples that cleared the bar — identically 0 on every rejected allele, which was two
        // rows in three of the tomato dump. Under the old name a join against an older dump
        // would keep working and compare two different quantities.
        writeln!(
            out,
            "contig\tstart\tend\tref\talt\tpassed\tbest_share_of_clearing_samples\tcohort_reads"
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
            // The bar alone, so `passed` means the bar and not the cap; the two ranking keys
            // are the module's own, read from the fold it just ran.
            let admitted = select_under(observation, bar, WIDEST_CAP, &mut scratch);
            let reference =
                String::from_utf8_lossy(&observation.alleles[REFERENCE_ALLELE]).into_owned();
            // The fold answers for this locus's table in that table's own order, so a
            // sequence and its ranking keys are the same index — walked together rather
            // than indexed apart.
            for (allele, bases) in observation.alleles.iter().enumerate().skip(1) {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
                    observation.region.contig.0,
                    observation.region.start.get(),
                    observation.region.end.get(),
                    reference,
                    String::from_utf8_lossy(bases),
                    u8::from(admitted.remap().candidate_for(allele).is_some()),
                    scratch.best_within_sample_share_of(allele),
                    scratch.cohort_reads_of(allele),
                )?;
            }
        }
        out.flush()?;
        println!("# allele dump written to {path}");
    }

    // **The same loci with the sample axis kept, which the dump above folds away.**
    // `NG_SELECT_ROWS=<path>` writes one row per (locus, sample, allele): the sample's pooled
    // reads on that allele and its compared reads at the locus. **The dump above cannot serve
    // as a test fixture and this can**, because the admission rule is per sample — a checked-in
    // test that re-derives it needs each sample's own numerator and denominator, where the
    // dump carries a cohort total and one maximised share.
    //
    // **Read groups are pooled here, as the shipped fold pools them** (`one_run_per_allele`),
    // so a fixture built from this exercises the rule but not the pooling; that half is the
    // module's own unit tests'.
    if let Ok(path) = std::env::var("NG_SELECT_ROWS") {
        use std::io::Write;
        let mut out = std::io::BufWriter::new(std::fs::File::create(&path)?);
        writeln!(
            out,
            "contig\tstart\tend\tallele\tbases\tsample\treads\tcompared_reads"
        )?;
        for observation in &built {
            for sample in &observation.per_sample {
                let compared: u32 = sample
                    .supported
                    .iter()
                    .map(|row| row.support.num_reads)
                    .sum();
                let mut pooled: Vec<u32> = vec![0; observation.alleles.len()];
                for row in &sample.supported {
                    pooled[row.allele] += row.support.num_reads;
                }
                // **Every allele of the table, including the ones this sample showed no reads
                // for.** An earlier version skipped the zero-read rows, and at a homozygous site
                // that silently dropped the *reference* — every read carries the alternative, so
                // allele 0 has no row — which cost a fixture built from this 221 loci of 7,478 at
                // 300× and 279 of 4,177 at 30×, every one of them a clean homozygous truth
                // variant. A consumer that wants only the non-zero rows can filter; one that
                // needs the table's shape cannot invent it.
                for (allele, reads) in pooled.iter().enumerate() {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        observation.region.contig.0,
                        observation.region.start.get(),
                        observation.region.end.get(),
                        allele,
                        String::from_utf8_lossy(&observation.alleles[allele]),
                        sample.sample,
                        reads,
                        compared,
                    )?;
                }
            }
        }
        out.flush()?;
        println!("# per-sample allele rows written to {path}");
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
        "samples",
        "alts/locus",
        "loci with none",
        format!("loci over {PRODUCTION_CAP}"),
        "max alts"
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
            let asked = restrict_to_first_samples(observation, size);
            let admitted = select_under(&asked, bar, WIDEST_CAP, &mut scratch);
            let a = admitted.alternative_allele_count();
            alts += a as u64;
            if a == 0 {
                none += 1;
            }
            if a + 1 > PRODUCTION_CAP {
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
