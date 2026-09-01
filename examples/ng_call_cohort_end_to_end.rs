//! **ng calls genotypes from real alignment files, and says where the time went.**
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --features merge-timing \
//!     --example ng_call_cohort_end_to_end -- \
//!     <reference.fa> <catalog.parquet> <regions.bed> <cram-dir>
//! ```
//!
//! `NG_SAMPLES=n` calls only the first `n` alignment files of the directory, in name order;
//! `NG_REGIONS=n` only the first `n` intervals of the BED. `NG_COVER=serial|parallel` picks
//! the arm: `serial` (default) is `call_cohort`, the oracle, one thread throughout;
//! `parallel` is the run's own record path, whose cover sweeps the samples concurrently
//! (Milestone E1) — compare the two on one ground to see what the parallel cover buys.
//! **Both default to a handful**, not
//! to everything, because this is meant to run inside the development loop: at its defaults it
//! takes about **5 seconds**, of which nearly 3 are reading and checksumming the reference. The
//! 63 tomato accessions over all 80 regions is a different measurement and should be asked for
//! deliberately.
//!
//! # The two questions
//!
//! **First, that ng calls a cohort of real reads at all.** Every test behind
//! `AlignedFilesVariantCaller::call_cohort` is a fabricated BAM of three or four reads over a
//! hundred bases of a reference that is a single homopolymer. This is the first run over
//! sequenced DNA, a real repeat catalog, and ground that has repeat tracts, gaps and contig
//! changes in it.
//!
//! **Second, where a calling run spends its time**, which is the measurement Milestone E's
//! shape is decided from (`doc/devel/ng/spec/run_streaming.md` §11, question 7 — *nobody has
//! measured it*). Two arrangements can genotype several loci at once — the merge's own region
//! batching switched on, so each thread assembles and genotypes its own stretch of ground, or
//! the merge left on one thread handing each finished locus to workers that only genotype —
//! and which is worth building depends on how the time divides between assembling loci and
//! calling them. This prints that division:
//!
//! - **drawing the readers forward** — every sample's walk, which is the reads being decoded;
//! - **assembling the loci** — closing the cohort's loci over the ground and folding each
//!   covering sample's reads onto their alleles;
//! - **genotyping them** — candidate selection, evidence shaping and the calling loop, which
//!   is what happens to each locus after it is assembled.
//!
//! **The last two are the split that decides the milestone**, and they are only separable
//! because calling happens inside the builder: the run's own stopwatch cannot see inside a
//! merge that returns everything at once. `--features merge-timing` is what turns the counters
//! on; **without it every time below reads zero** and only the wall clock and the counts are
//! real.
//!
//! # What it measured, 2026-09-01
//!
//! Tomato accessions from `benchmarks/tomato1/crams/` over the first **two** intervals of
//! `benchmarks/tomato1/regions.bed` — 200 kb of SL4.0 at about three reads a position — in the
//! development container, release build, `--features merge-timing`:
//!
//! | samples | compressed MB | loci called | `call_cohort` | drawing the readers | assembling | genotyping |
//! |---|---|---|---|---|---|---|
//! | 3 | 123.5 | 3,291 | 0.61 s | 589.2 ms (97.3%) | 7.4 ms (1.2%) | 5.8 ms (1.0%) |
//! | 6 | 216.8 | 4,235 | 1.07 s | 1,040.7 ms (96.8%) | 14.9 ms (1.4%) | 12.9 ms (1.2%) |
//! | 12 | 383.1 | 5,675 | 1.98 s | 1,895.9 ms (96.0%) | 32.9 ms (1.7%) | 34.1 ms (1.7%) |
//! | 24 | 672.8 | 8,825 | 3.89 s | 3,671.6 ms (94.3%) | 89.6 ms (2.3%) | 103.1 ms (2.6%) |
//! | 63 | 1,840.7 | 23,450 | 20.55 s | 18,105.0 ms (88.1%) | 1,135.0 ms (5.5%) | 1,090.9 ms (5.3%) |
//!
//! **Decoding reads is `call_cohort`.** Assembling and genotyping together are 2.2% of it at
//! three samples, 4.9% at twenty-four and 10.8% at the whole cohort of sixty-three — and the
//! shape survives the ground: the same sixty-three accessions over all eighty intervals
//! (8 Mb, 1,069,772 loci, 723.2 s on one thread, measured 2026-09-01 at Milestone E) split
//! 87.9% drawing, 5.4% assembling, 5.9% genotyping. **That split is what shaped Milestone E**:
//! the parallelism went to the cover (`NG_COVER=parallel` above), and no pool of genotyping
//! workers was built, because at every measured size such a pool reaches at most the
//! genotyping share.
//!
//! **But `call_cohort` is not the whole wait.** Reading and checksumming the 795 MB reference,
//! opening the catalog and building the segments cost **2.73 to 2.81 seconds** across those
//! four runs — constant in the cohort and in the ground — so at this probe's defaults they are
//! more than half of the 4.8 seconds a person waits. They are a row of the output.
//!
//! **That 88–97% was one thread until Milestone E1 (2026-09-01), and this probe's two arms
//! are how the fix was measured.** `call_cohort` still draws every sample forward one after
//! another (`ObservationCache::cover`) — it is the oracle — while the run's record path sweeps
//! the samples concurrently. Three alternated pairs on a quiet machine, 63 samples over this
//! same 200 kb: calling went from **12.06–12.32 s serial to 6.61–6.68 s parallel** on 8 rayon
//! threads — **1.8×** — and on the whole 8 Mb benchmark ground from **473.6 s to 308.8 s** —
//! **1.5×** — identical loci everywhere. The parallel arm spends about half again the CPU
//! (22.2–22.8 s of user time against 14.7–15.0 on the slice; 13m37 against 7m59 on the whole
//! ground): per-cover scheduling, Jacobi re-sweeps, and records freed on a different thread
//! from the one that allocated them (Milestone G's territory). The gap to the naive 8-thread
//! ceiling is the cover's granularity — one 200-base building region per sweep, whose barrier
//! waits for the slowest sample's decode each time.
//!
//! # The two rates, and what the fifth row did to one of them
//!
//! **The per-locus calling rate is stable and is what a later run should be compared
//! against**: about 1 microsecond per locus per sample — 1.34, 1.09, 0.98, 0.91 and 0.74
//! across the five cohorts, flat and falling.
//!
//! **The reading rate is not.** 4.77, 4.80, 4.95 and 5.46 ms per compressed megabyte at 3–24
//! samples, then **9.84 at 63** — the fifth row broke the "linear to within 14%" claim this
//! header used to carry, and why the whole-cohort rate is near double is unmeasured. Compare a
//! later run against the range, 4.8–9.8 ms/MB, not against 5.
//!
//! **Neither is a share, and the share is what Milestone E turned on.** Calling's share grows
//! 2.2% → 4.9% → 10.8% across 3 → 24 → 63 samples, and **the whole of that growth is the
//! number of loci**, not the cost of one: more accessions segregate more sites, so the count
//! goes 3,291 to 23,450 while the cost of each falls. **That curve has to flatten** — 200 kb
//! of SL4.0 holds a finite number of segregating sites — and where it flattens is what decides
//! whether calling is a tenth or a third of a thousand-sample run. **This probe cannot see
//! that**, and no extrapolation from five cohorts of at most sixty-three should be trusted to:
//! the answer is to run it at a thousand samples, which nothing here prevents.
//!
//! **What this must not be read as saying is that reading scales sublinearly in the cohort.**
//! It does not. The first three accessions in name order happen to be 1.5 times the average
//! file size, so eight times the samples is only 5.45 times the bytes — which is why the table
//! carries the megabytes, and why the rate to compare is per megabyte rather than per sample.
//!
//! **The two per-building-region costs inside "assembling" are nothing here** — 0.0 to 0.3 ms
//! against 7 to 90 ms of assembling — so on this ground the assembling is per-locus work almost
//! entirely, and dividing the genome more finely would not cost the run much.
//!
//! # What it does not do
//!
//! **It writes no VCF.** The serial arm prints what the caller produced — called loci, and per
//! sample how many were called, set aside, or called as carrying an alternative. The parallel
//! arm prints record totals and no per-sample lines: its records are dropped at the sink, so
//! the disk stays out of the timing.
//!
//! **Every locus goes down the SNP/indel path.** Repeat-tract candidate selection is specified
//! and unbuilt, so both tract generator slots are refused as such: a repeat tract in the
//! analysed ground is counted as ground this caller has not built yet, and the per-sample walk
//! tallies say how much. A run over tract-rich ground is therefore **short, not wrong**, and
//! the number to read is `not built yet` below.
//!
//! # The catalog is given by path, not found beside the reference
//!
//! `RepeatCatalog::open_beside_reference` is the convention and cannot be used here: the
//! benchmark reference lives under `$HOME/genomes`, which the development container mounts
//! read-only, so no catalog can be written beside it. The catalog for that reference is inside
//! the project tree instead, and this takes its path — as
//! `examples/ng_open_cohort_descriptors.rs` does, for the same reason.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::calling::SampleGenotypeCall;
use pop_var_caller::ng::calling::allele_candidates::CandidateSelectionConfig;
use pop_var_caller::ng::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use pop_var_caller::ng::calling::inference::CallingLoopConfig;
use pop_var_caller::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use pop_var_caller::ng::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use pop_var_caller::ng::calling::parameters_file::DeclaredInbreeding;
use pop_var_caller::ng::calling::run_parameters::RunParameters;
use pop_var_caller::ng::locus_generation::pileup::PileupGeneratorConfig;
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::GenomeRegions;
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::run::cohort_merge::timing;
use pop_var_caller::ng::run::{
    AlignedFilesVariantCaller, AlignmentInputs, AssemblyCheckOutcome, CalledCohort,
    MergeParameters, Segmentation,
};
use pop_var_caller::ng::types::{Genotype, Ploidy};
use pop_var_caller::regions::ContigBounds;

/// How many alignment files a run takes when nothing says otherwise.
///
/// **A handful, because this is a development-loop probe.** Six samples over four of the
/// tomato benchmark's hundred-kilobase intervals is the shape that answers both questions in
/// minutes; the whole cohort is a different measurement and `NG_SAMPLES` is how it is asked
/// for.
const SAMPLES_BY_DEFAULT: usize = 6;

/// How many analysed intervals a run takes when nothing says otherwise.
const REGIONS_BY_DEFAULT: usize = 4;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [fasta, catalog, bed, crams] = args.as_slice() else {
        eprintln!(
            "usage: ng_call_cohort_end_to_end <reference.fa> <catalog.parquet> <regions.bed> \
             <cram-dir>\n\
             calls a cohort of alignment files end to end and reports where the time went.\n\
             NG_SAMPLES=n calls the first n files (default {SAMPLES_BY_DEFAULT}); \
             NG_REGIONS=n the first n BED intervals (default {REGIONS_BY_DEFAULT})."
        );
        return ExitCode::from(2);
    };
    match run(
        Path::new(fasta),
        Path::new(catalog),
        Path::new(bed),
        Path::new(crams),
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            let mut cause = error.source();
            while let Some(next) = cause {
                eprintln!("  caused by: {next}");
                cause = next.source();
            }
            ExitCode::FAILURE
        }
    }
}

/// How many of something this run takes: the environment's answer, or `fallback`.
fn how_many(name: &str, fallback: usize) -> usize {
    std::env::var(name)
        .ok()
        .map(|value| value.parse().expect("a count"))
        .unwrap_or(fallback)
}

fn run(
    fasta: &Path,
    catalog_path: &Path,
    bed: &Path,
    crams: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    // **Everything before the first read is decoded is timed too.** Reading and checksumming
    // an 800 MB reference and opening a 60 MB catalog cost seconds, are constant in the cohort
    // and in the ground, and are most of what a person waits for at this probe's defaults — so
    // a breakdown that began at `call_cohort` would describe a minority of the run.
    let setting_up = Instant::now();

    // **The reference with its FASTA read to the end**, because the assembly check and the
    // catalog check both compare checksums and a reference whose background verification has
    // not been joined carries none.
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let with_checksums = match verify {
        Some(handle) => handle.join()?,
        None => Arc::clone(&info),
    };
    let contigs: ContigList = info.contig_list();
    let reference = OpenReference::new(info);

    let bounds: Vec<ContigBounds<'_>> = contigs
        .entries
        .iter()
        .map(|entry| ContigBounds {
            name: &entry.name,
            length: entry.length as u32,
        })
        .collect();
    let scratch = tempfile::tempdir()?;
    let trimmed = first_regions_of(
        bed,
        how_many("NG_REGIONS", REGIONS_BY_DEFAULT),
        scratch.path(),
    )?;
    let analysed = GenomeRegions::from_bed_path(&trimmed, &bounds)?;
    let analysed_bases: u64 = analysed.iter().map(|region| region.len()).sum();

    let criteria = StrRepeatCriteria::default();
    let catalog = RepeatCatalog::open_checking_against_reference(catalog_path, &with_checksums)?;
    let spans: Vec<_> = analysed.iter().collect();
    let segments = catalog.genome_segments(&criteria, ReadScope::Regions(&spans))?;
    let segmentation = Segmentation::build(
        segments,
        analysed,
        catalog.header().clone(),
        criteria,
        catalog_path.to_path_buf(),
    )?;

    let mut paths: Vec<PathBuf> = std::fs::read_dir(crams)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| {
            path.extension()
                .is_some_and(|kind| kind == "cram" || kind == "bam")
        })
        .collect();
    paths.sort();
    paths.truncate(how_many("NG_SAMPLES", SAMPLES_BY_DEFAULT));
    if paths.is_empty() {
        return Err(format!("no .cram or .bam under {}", crams.display()).into());
    }

    let compressed_bytes: u64 = paths
        .iter()
        .filter_map(|path| std::fs::metadata(path).ok())
        .map(|file| file.len())
        .sum();

    let read_groups = build_read_groups(&paths)?;
    // **Defaulted parameters, and the run says so.** Nothing here has been fitted from this
    // cohort: what this probe measures is the calling path's cost and that it produces
    // genotypes at all, not what those genotypes are worth.
    let parameters = RunParameters::of_defaults(
        &read_groups,
        Ploidy::try_new(2).expect("a diploid"),
        &DeclaredInbreeding::nothing_said(),
    );

    println!("# reference: {}", fasta.display());
    println!("# analysed intervals: {}", spans.len());
    println!("# analysed bases: {analysed_bases}");
    println!("# alignment files: {}", paths.len());
    println!(
        "# compressed alignment bytes: {compressed_bytes} ({:.1} MB)",
        compressed_bytes as f64 / 1e6,
    );
    println!(
        "# segments: {} over {} analysed region(s)",
        segmentation.segments().len(),
        segmentation.analysed_regions().len(),
    );

    let opening = Instant::now();
    let caller = AlignedFilesVariantCaller::open(
        AlignmentInputs {
            read_groups: &read_groups,
            reference: &reference,
            read_filters: ReadFilterConfig::default(),
            build_index_if_missing: false,
            locus_generator_settings: PileupGeneratorConfig::default(),
            reference_with_checksums: &with_checksums,
        },
        segmentation,
        parameters,
        CallingLoopConfig::DEFAULT
            .validate()
            .expect("the shipped calling-loop settings are runnable"),
        CandidateSelectionConfig::DEFAULT,
        MergeParameters::DEFAULT,
    )?;
    let opening_seconds = opening.elapsed().as_secs_f64();
    let setup_seconds = setting_up.elapsed().as_secs_f64();
    println!("# samples: {}", caller.sample_count());
    println!("# opening every sample's files: {opening_seconds:.2} s");
    println!(
        "# everything before the first read is decoded: {setup_seconds:.2} s \
         (the reference, the catalog, the segments, and the opens above)"
    );

    // **The counters are global and the report is not reentrant**, so they are put back to
    // zero immediately before the run that is being measured.
    timing::reset();
    let genotyper =
        SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior);

    // **Two arms, chosen by `NG_COVER`.** `serial` (the default) is `call_cohort`, the
    // oracle, whose cover draws every sample on one thread; `parallel` is
    // `call_cohort_handing_each_record_over`, the path the command takes, whose cover sweeps
    // the samples concurrently (Milestone E1). Same walk, same merge, same calls per locus —
    // what the two arms measure against each other is the cover's schedule, which is why the
    // probe grew the switch.
    let cover = std::env::var("NG_COVER").unwrap_or_else(|_| "serial".to_string());
    match cover.as_str() {
        "serial" => {
            let calling = Instant::now();
            let called = caller.call_cohort(&genotyper)?;
            let calling_seconds = calling.elapsed().as_secs_f64();
            println!("# cover: serial (call_cohort, the oracle)");
            report_the_calls(&called);
            report_the_ground(&called, analysed_bases);
            report_where_the_time_went(calling_seconds, setup_seconds, compressed_bytes);
        }
        "parallel" => {
            let calling = Instant::now();
            // The records are dropped where they are handed over: this arm measures the
            // calling path, and writing a VCF would time the disk beside it. The count comes
            // back in the answer (`records_written`), so nothing here keeps one.
            let written = caller.call_cohort_handing_each_record_over(
                &genotyper,
                &mut |_record| -> Result<(), std::io::Error> { Ok(()) },
            )?;
            let calling_seconds = calling.elapsed().as_secs_f64();
            println!(
                "# cover: parallel (call_cohort_handing_each_record_over, the run's path) \
                 over {} rayon threads",
                rayon::current_num_threads(),
            );
            println!(
                "# loci called: {} — {} written as records, {} establishing no variant",
                written.loci_called(),
                written.records_written,
                written.loci_called_but_not_written,
            );
            println!(
                "# loci the merge declined to assemble for being too wide: {}",
                written.loci_too_wide_to_assemble.len(),
            );
            report_where_the_time_went(calling_seconds, setup_seconds, compressed_bytes);
        }
        other => {
            return Err(format!("NG_COVER={other}: it answers `serial` or `parallel`").into());
        }
    }
    Ok(())
}

/// The BED's first `how_many` intervals, as a BED of their own under `scratch`.
///
/// **Written back out rather than filtered in memory**, because `GenomeRegions` is built from a
/// BED path or from whole contigs and from nothing else — and adding a third constructor to a
/// production type for a probe's convenience is the wrong way round.
///
/// **Truncating rather than sampling** keeps a re-run comparable with the one before it, which
/// is what a probe read inside a development loop needs.
fn first_regions_of(
    bed: &Path,
    how_many: usize,
    scratch: &Path,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let whole = std::fs::read_to_string(bed)?;
    let kept: Vec<&str> = whole
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .take(how_many)
        .collect();
    if kept.is_empty() {
        return Err(format!("{} holds no intervals to call", bed.display()).into());
    }
    let trimmed = scratch.join("analysed.bed");
    std::fs::write(&trimmed, format!("{}\n", kept.join("\n")))?;
    Ok(trimmed)
}

/// What the cohort was called: how many loci, and what each sample got.
fn report_the_calls(called: &CalledCohort) {
    println!("# loci called: {}", called.called_loci.len());
    println!(
        "# loci the merge declined to assemble for being too wide: {}",
        called.loci_too_wide_to_assemble.len()
    );
    // **A third fact, and not a rename of either of the others.** Too wide is ground the merge
    // would not assemble; too quiet is ground no sample varied at and is counted nowhere; this
    // is ground that was assembled and where the allele cap then left nobody callable.
    //
    // **No flag is named here because none exists yet.** The candidate cap is
    // `CandidateSelectionConfig::max_candidate_alleles`, which this probe leaves at its shipped
    // value and the command surface does not expose (that is the plan's Milestone F); telling
    // an operator to raise a flag they cannot type would be worse than telling them nothing.
    println!(
        "# loci where the allele cap left nobody callable: {} (assembled, and no sample of the \
         run was left to call: every sample covered them and every one had earned a sequence \
         the cap cut)",
        called.loci_with_nobody_to_call.len()
    );
    println!(
        "# loci the merge assembled: {} (called plus nobody-callable — a locus nobody varied \
         at is counted nowhere, by design)",
        called.called_loci.len() + called.loci_with_nobody_to_call.len(),
    );
    if called.called_loci.is_empty() {
        println!("# (nothing was called, so there is nothing to break down per sample)");
        return;
    }

    let alleles: usize = called
        .called_loci
        .iter()
        .map(|locus| locus.alleles().len())
        .sum();
    println!(
        "# alleles a locus is called over, mean: {:.2} (the reference counts as one)",
        alleles as f64 / called.called_loci.len() as f64
    );

    // **The third number is a subset of the second, not a third part of the total.** Printed
    // as "of those" because three counts in a row read as a partition, and these are not one.
    println!("# sample, loci called, loci set aside, of those called: carrying an alternative");
    for (sample, walk) in called.walk.per_sample.iter().enumerate() {
        let mut called_here = 0_usize;
        let mut set_aside = 0_usize;
        let mut carrying = 0_usize;
        for locus in &called.called_loci {
            match &locus.per_sample[sample] {
                SampleGenotypeCall::Missing => set_aside += 1,
                SampleGenotypeCall::Called { genotype, .. } => {
                    called_here += 1;
                    if carries_an_alternative(genotype) {
                        carrying += 1;
                    }
                }
            }
        }
        println!(
            "{}, {called_here}, {set_aside}, {carrying}",
            walk.sample_name
        );
    }
}

/// Whether a genotype carries any allele but the reference.
fn carries_an_alternative(genotype: &Genotype) -> bool {
    genotype
        .alleles()
        .iter()
        .any(|allele| !allele.is_reference())
}

/// What each sample's walk covered, and what it could not.
///
/// **The two kinds of nothing are kept apart** and a reader must not add them: ground no
/// generator is built for **yet** is this caller's own gap, where a satellite is ground it will
/// **never** call. And the share that matters is in bases: typed regions differ in length by
/// orders of magnitude, so half the regions can be a twentieth of the ground.
fn report_the_ground(called: &CalledCohort, analysed_bases: u64) {
    println!("# {}", the_assembly_check(called));

    let first = called
        .walk
        .per_sample
        .first()
        .map(|walk| &walk.regions)
        .expect("a run has at least one sample");
    println!(
        "# typed regions this caller walked: {} of {}; not built yet: {} regions, {} bases \
         ({:.1}% of the analysed ground); never called (satellite): {} regions, {} bases",
        first.regions_handled,
        first.regions_in,
        first.unhandled_not_implemented,
        first.unhandled_not_implemented_bp,
        100.0 * first.unhandled_not_implemented_bp as f64 / analysed_bases.max(1) as f64,
        first.unhandled_out_of_scope,
        first.unhandled_out_of_scope_bp,
    );
    println!(
        "# (every sample walks the same ground, so the line above is the run's; the loci each \
         sample's walk emitted differ and are below)"
    );

    println!(
        "# sample, loci its walk emitted, reads admitted, positions short of the hold ceiling, \
         columns the per-position caps cut"
    );
    for walk in &called.walk.per_sample {
        let counted = walk.snp_indel;
        println!(
            "{}, {}, {}, {}, {}",
            walk.sample_name,
            walk.regions.loci_emitted,
            counted.map_or(0, |counts| counts.reads_admitted),
            counted.map_or(0, |counts| counts.positions_short_of_cap),
            counted.map_or(0, |counts| counts.column_depth_truncations),
        );
    }
}

/// **What the run could learn about the assembly its samples were aligned to**, as a sentence.
///
/// Not the enum's `Debug`: this is the one line of the output that answers *were these files
/// aligned to this reference*, and a reader of a run's report should not have to decode a Rust
/// struct to get it. The two cases say different things and only one is reassuring.
fn the_assembly_check(called: &CalledCohort) -> String {
    match called.walk.assembly_check {
        AssemblyCheckOutcome::EverySampleMatchedTheReference {
            alignment_files,
            checksums_compared,
            checksums_possible,
        } => format!(
            "assembly check: every one of {alignment_files} alignment file(s) matched the \
             reference; {checksums_compared} of {checksums_possible} contig checksums could be \
             compared and all agreed"
        ),
        AssemblyCheckOutcome::NothingCouldBeChecked { because } => format!(
            "assembly check: NOT ONE checksum could be compared ({because:?}), so nothing here \
             says these files were aligned to this reference"
        ),
    }
}

/// **Where a calling run's time went** — the measurement spec §11 question 7 asks for.
///
/// `calling_seconds` is the whole of `call_cohort`, timed from outside it, and `setup_seconds`
/// is everything before it. The rest comes from the merge's own counters and is **zero without
/// `--features merge-timing`**, which is said out loud rather than left to be inferred from a
/// table of zeros.
fn report_where_the_time_went(calling_seconds: f64, setup_seconds: f64, compressed_bytes: u64) {
    let counted = timing::report(rayon::current_num_threads());
    println!(
        "# what a person waits for: {:.2} s — {setup_seconds:.2} s before the first read is \
         decoded, {calling_seconds:.2} s calling the cohort",
        setup_seconds + calling_seconds,
    );
    if counted.merge_wall_ms == 0.0 {
        println!(
            "# (built without --features merge-timing, so the breakdown below is all zeros: \
             re-run with it to get the split)"
        );
    }

    // **Milliseconds and a share of `call_cohort`, not seconds.** The first run of this probe
    // printed two decimals of a second and read `0.00` for three of its four rows, which says
    // nothing about a split whose whole purpose is to be compared.
    let calling_ms = calling_seconds * 1e3;
    let share = |part: f64| {
        if calling_ms > 0.0 {
            100.0 * part / calling_ms
        } else {
            0.0
        }
    };
    println!("# inside those {calling_seconds:.2} s — milliseconds, and share of them:");
    for (what, part) in [
        ("drawing the readers forward", counted.cover_ms),
        ("evicting what the merge has passed", counted.evict_ms),
        ("assembling the loci", counted.assembling_loci_ms()),
        ("genotyping them", counted.after_assembly_ms),
    ] {
        println!("{what}, {part:.1}, {:.1}%", share(part));
    }
    // **Two of those milliseconds are charged per building region rather than per locus**, and
    // they sit inside "assembling the loci" above. Which of the two Milestone E's arrangements
    // is worth building turns on that distinction, so they are broken out rather than left for
    // somebody to discover in `timing.rs`.
    println!(
        "# of the assembling, per building region rather than per locus: {:.1} ms building the \
         per-sample windows, {:.1} ms setting each region's walk up",
        counted.window_ms, counted.walk_setup_ms,
    );
    println!(
        "# merge working windows: {} ({} held no locus); the merge's own wall clock: {:.1} ms",
        counted.regions, counted.regions_with_no_locus, counted.merge_wall_ms,
    );

    // The two rates this probe exists to produce. Both are stable where a share is not: a
    // share moves with the cohort because the locus count does, and neither of these does.
    if compressed_bytes > 0 && counted.cover_ms > 0.0 {
        // **Per megabyte of the files opened, not of the bytes decoded**, which nothing
        // counts: a run walks only the analysed ground, so this rate is comparable across
        // cohort sizes over the *same* intervals and not across different ones. Doubling the
        // regions doubles it while the files are unchanged.
        println!(
            "# reading: {:.2} ms per compressed megabyte of the files opened — comparable \
             across cohort sizes over the same intervals, not across different ground",
            counted.cover_ms / (compressed_bytes as f64 / 1e6),
        );
    }
}
