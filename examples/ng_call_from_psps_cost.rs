//! **ng calls a cohort from stored psp files, and says where the time went.**
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --features merge-timing \
//!     --example ng_call_from_psps_cost -- \
//!     <reference.fa> <catalog.parquet> <psp-dir>
//! ```
//!
//! `NG_SAMPLES=n` calls only the first `n` psp files of the directory, in name order.
//!
//! # The question
//!
//! **Nobody has ever timed a psp-mode run.** Direct mode's own probe
//! (`ng_call_cohort_end_to_end.rs`) measured that at 63 tomato accessions over 200 kb,
//! **88.1% of a 20.55 s call is drawing the readers forward** — which in direct mode is reads
//! being decoded out of CRAM — against 5.5% assembling the cohort's loci and 5.3% genotyping
//! them. psp mode replaces that first term with reading records back out of a file somebody
//! else already walked, and **what that does to the other two shares is the number this
//! prints**.
//!
//! It is asked because it gates a design. `doc/devel/ng/spec/cohort_merge_psp_path.md` proposes
//! that a psp-mode run stop building every stored observation and build only the roughly one
//! locus in a hundred the cohort keeps, deciding the rest from each record's head. The saving
//! that design can possibly deliver is bounded by decoding's share of a run, and **that share
//! is what is measured here, before any of it is built** (its plan's Milestone A, and its
//! checkpoint: if decoding plus merging is a small part of a run, most of the design is not
//! worth building).
//!
//! # What it prints, and what each row means
//!
//! The same four-way split direct mode's probe prints, so the two are read side by side:
//!
//! - **drawing the readers forward** — every sample's records coming out of its psp, which is
//!   block decompression plus building each record;
//! - **evicting what the merge has passed** — releasing records behind the merge's window;
//! - **assembling the loci** — closing the cohort's loci and folding each covering sample's
//!   reads onto their alleles;
//! - **genotyping them** — candidate selection, evidence shaping and the calling loop.
//!
//! **Without `--features merge-timing` every one of those reads zero** and only the wall clock
//! and the counts are real. The feature is what turns the merge's counters on.
//!
//! # The ground is the files' own, not this probe's
//!
//! **The analysed regions come from the psp cohort itself**, not from a BED given here. Every
//! psp records the ground it was walked over, a calling run refuses a cohort whose files
//! disagree about it, and the segmentation a run calls against must be built over that same
//! ground or `PspVariantCaller::open` refuses it. Taking the regions from the cohort makes that
//! agreement structural: there is no BED argument to get wrong, and a probe run cannot
//! accidentally measure a segmentation the files were not written under.
//!
//! What must still match is the **catalog** and the repeat criteria: the psps record the
//! segmentation inputs their typing used, so a catalog other than the one `generate-psps` was
//! given is refused by name.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::calling::allele_candidates::CandidateSelectionConfig;
use pop_var_caller::ng::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
use pop_var_caller::ng::calling::inference::CallingLoopConfig;
use pop_var_caller::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use pop_var_caller::ng::calling::likelihood::ssr_emission::StutterSubstitutionEmission;
use pop_var_caller::ng::calling::parameters_file::DeclaredInbreeding;
use pop_var_caller::ng::calling::run_parameters::RunParameters;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::DEFAULT_MAX_STR_LEN;
use pop_var_caller::ng::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use pop_var_caller::ng::run::cohort_merge::timing;
use pop_var_caller::ng::run::{
    MergeParameters, OpenPspCohort, PspVariantCaller, StoredCohortInputs,
};
use pop_var_caller::ng::types::Ploidy;
use pop_var_caller::pop_var_caller_exp::run_ground::{
    self, GroundRequest, RepeatRouting,
};

/// How many psps a run takes when nothing says otherwise.
///
/// **A handful, not everything**, for the reason direct mode's probe gives: this is meant to
/// run inside the development loop. The 63-accession cohort is a different measurement and is
/// asked for deliberately with `NG_SAMPLES=63`.
const SAMPLES_BY_DEFAULT: usize = 6;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [fasta, catalog, psps] = args.as_slice() else {
        eprintln!(
            "usage: ng_call_from_psps_cost <reference.fa> <catalog.parquet> <psp-dir>\n\
             calls a cohort of stored psp files and reports where the time went.\n\
             NG_SAMPLES=n calls the first n psps (default {SAMPLES_BY_DEFAULT}).\n\
             The analysed ground comes from the psps themselves, so there is no BED argument."
        );
        return ExitCode::from(2);
    };
    match run(Path::new(fasta), Path::new(catalog), Path::new(psps)) {
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
    psp_dir: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    // **Everything before the first record is decoded is timed too**, on the same grounds
    // direct mode's probe states: reading and checksumming an 800 MB reference and opening a
    // 60 MB catalog cost seconds and are most of what a person waits for at these defaults.
    let setting_up = Instant::now();

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
    let reference = OpenReference::new(info);

    let mut paths: Vec<PathBuf> = std::fs::read_dir(psp_dir)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|kind| kind == "psp"))
        .collect();
    paths.sort();
    paths.truncate(how_many("NG_SAMPLES", SAMPLES_BY_DEFAULT));
    if paths.is_empty() {
        return Err(format!("no .psp under {}", psp_dir.display()).into());
    }

    let stored_bytes: u64 = paths
        .iter()
        .filter_map(|path| std::fs::metadata(path).ok())
        .map(|file| file.len())
        .sum();

    // **The cohort is opened before the segmentation is built, because it is what says over
    // what ground to build it.** Every refusal comparing the files with each other fires here,
    // before a block is decoded.
    let opening = Instant::now();
    let cohort = OpenPspCohort::open(&paths)?;
    let opening_seconds = opening.elapsed().as_secs_f64();

    let analysed = cohort.analysed_regions().clone();
    let analysed_bases: u64 = analysed.iter().map(|region| region.len()).sum();

    // **The segmentation is built through the subcommands' own path, not assembled here.**
    // A run's repeat criteria come from five routing flags, and their defaults are *not*
    // `StrRepeatCriteria::default()` — that value is the catalog file's own storage floors,
    // which would route about seven times more reference to the repeat path than a run's
    // calling floors do (`run_ground::routing_criteria`). Building them any other way than
    // `generate-psps` did earns the refusal this probe met on its first run: *"written under a
    // different set of repeat-tract criteria from this run's"*. Sharing the function is what
    // makes the two agree by construction rather than by matching two lists of defaults.
    let ground = GroundRequest {
        reference: fasta,
        catalog: Some(catalog_path),
        regions: None,
        routing: RepeatRouting {
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        },
    };
    let spans: Vec<_> = analysed.iter().collect();
    let segmentation = run_ground::segments_over(&ground, &analysed, &with_checksums)?;

    // **Defaulted parameters, and the run says so.** Nothing here has been fitted from this
    // cohort: what this probe measures is the calling path's cost, not what its genotypes are
    // worth.
    let parameters = RunParameters::of_defaults(
        cohort.read_groups(),
        Ploidy::try_new(2).expect("a diploid"),
        &DeclaredInbreeding::nothing_said(),
    );

    println!("# reference: {}", fasta.display());
    println!("# psp files: {}", paths.len());
    println!(
        "# stored psp bytes: {stored_bytes} ({:.1} MB)",
        stored_bytes as f64 / 1e6,
    );
    println!("# analysed intervals: {}", spans.len());
    println!("# analysed bases: {analysed_bases}");
    println!("# segments: {}", segmentation.segments().len());
    println!("# opening every sample's psp: {opening_seconds:.2} s");

    let caller = PspVariantCaller::open(
        cohort,
        StoredCohortInputs {
            reference: &reference,
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
    let setup_seconds = setting_up.elapsed().as_secs_f64();
    println!("# samples: {}", caller.sample_count());
    println!(
        "# everything before the first record is decoded: {setup_seconds:.2} s \
         (the reference, the catalog, the segments, and the opens above)"
    );

    // **The counters are global and the report is not reentrant**, so they are put back to zero
    // immediately before the run being measured.
    timing::reset();
    let genotyper =
        SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior);

    let calling = Instant::now();
    // The records are dropped where they are handed over: this measures the calling path, and
    // writing a VCF would time the disk beside it. The counts come back in the answer.
    // **Two tallies since Milestone F1**: what the calling did, and what each stored sample
    // contributed. This probe reports the first; the second is what the subcommand's run
    // report states per sample, and nothing here needs it.
    let (written, _per_sample) = caller.call_cohort_handing_each_record_over(
        &genotyper,
        &mut |_record| -> Result<(), std::io::Error> { Ok(()) },
    )?;
    let calling_seconds = calling.elapsed().as_secs_f64();

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
    report_where_the_time_went(calling_seconds, setup_seconds, stored_bytes);
    Ok(())
}

/// **Where a psp-mode calling run's time went.**
///
/// `calling_seconds` is the whole call, timed from outside it, and `setup_seconds` is
/// everything before. The rest comes from the merge's own counters and is **zero without
/// `--features merge-timing`**, said out loud rather than left to be inferred from a table of
/// zeros. Deliberately the same rows, in the same order and units, as direct mode's probe, so
/// the two runs can be read against each other.
fn report_where_the_time_went(calling_seconds: f64, setup_seconds: f64, stored_bytes: u64) {
    let counted = timing::report(rayon::current_num_threads());
    println!(
        "# what a person waits for: {:.2} s — {setup_seconds:.2} s before the first record is \
         decoded, {calling_seconds:.2} s calling the cohort",
        setup_seconds + calling_seconds,
    );
    if counted.merge_wall_ms == 0.0 {
        println!(
            "# (built without --features merge-timing, so the breakdown below is all zeros: \
             re-run with it to get the split)"
        );
    }

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
    println!(
        "# of the assembling, per building region rather than per locus: {:.1} ms building the \
         per-sample windows, {:.1} ms setting each region's walk up",
        counted.window_ms, counted.walk_setup_ms,
    );
    println!(
        "# merge working windows: {} ({} held no locus); the merge's own wall clock: {:.1} ms",
        counted.regions, counted.regions_with_no_locus, counted.merge_wall_ms,
    );
    println!(
        "# of the drawing: {:.1} ms the samples' own drawing summed over threads, which spread \
         perfectly over {} threads would be {:.1} ms of wall, against {:.1} ms measured; \
         {} cover sweeps over {} working windows",
        counted.cover_busy_ms,
        counted.threads,
        counted.cover_busy_ms / counted.threads as f64,
        counted.cover_ms,
        counted.cover_sweeps,
        counted.regions,
    );
    // **Records drawn is what the deferred-build design would change**, so it is printed beside
    // the time: the design builds only the records a kept locus needs, and this is the count it
    // would be cutting into.
    println!(
        "# records drawn: {} carrying {} observations",
        timing::RECORDS_DRAWN.get(),
        timing::OBSERVATIONS_DRAWN.get(),
    );

    if stored_bytes > 0 && counted.cover_ms > 0.0 {
        println!(
            "# reading: {:.2} ms per stored megabyte — comparable across cohort sizes over the \
             same intervals, not across different ground",
            counted.cover_ms / (stored_bytes as f64 / 1e6),
        );
    }
}
