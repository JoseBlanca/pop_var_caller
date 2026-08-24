//! **The histogram route's error rate, per read group, on real reads — for as many samples as
//! you hand it.**
//!
//! ```text
//! cargo run --release --example ng_histogram_error_rates -- \
//!     <reference.fa> <catalog.parquet> <regions.bed> <sample.bam|cram> [sample ...]
//! ```
//!
//! # Why this exists
//!
//! The parameter pre-pass estimates a per-base sequencing error rate **twice**, by two routes
//! that read different objects, and which of the two is kept has never been measured
//! (`spec/parameter_prepass.md` §4.1; `spec/parameter_prepass_joint_fit.md` §8 sets out the
//! comparison and its third measurement is deliberately on real reads). The other route — the
//! census fit at the gather — already had a driver that runs a whole cohort of real alignments,
//! `examples/ng_joint_records_walk.rs`. **This route had none.** Its only real-data exercise was
//! one `#[ignore]`d unit test taking a single sample from three environment variables
//! (`parameter_estimation::generic::real_alignments`), so comparing the two over a 63-accession
//! cohort meant 63 invocations and scraping standard error.
//!
//! This is that test's walk, taking samples as arguments and printing one line per read group.
//!
//! # It is the real route, not a re-derivation of it
//!
//! Same typed regions from the reference's repeat catalog, same left-aligning read preparer,
//! same pileup generator, and the fit is
//! [`estimate_generic_parameters`](pop_var_caller::ng::parameter_estimation::generic::estimate_generic_parameters)
//! itself — the pre-pass's own public entry point, which builds the accumulator, pours every
//! locus into it and runs the coupled error-rate/genotype-frequency fit. Nothing here reimplements
//! a fit.
//!
//! # Three things about the run that a reader of the numbers has to know
//!
//! **The inbreeding coefficient is supplied at zero, not fitted.** The runs model refuses unless
//! 3,000 separate 100 kb windows each hold a site, and both benchmark BEDs are too scattered to
//! reach that — tomato's is 8.0 Mb in 80 spans. The error rate is fitted jointly with the sample's
//! genotype frequencies, and those frequencies were fitted under a stated `F`, so this is a
//! property of the run and not a detail. It is printed on every line.
//!
//! **Every site is scored at ploidy 2.** `ConstantPloidy` is what the pre-pass builds and nothing
//! here overrides it. Both benchmark BEDs are autosomal.
//!
//! **The rate printed is the *marginal* one** — the rate a read disagrees at a site drawn at
//! random, which is what a sample emits. Where the fit found a second, noisier class of site, the
//! clean and noisy rates behind it are printed too, because the other route reports that pair
//! rather than the marginal and folding one into the other is how the two are compared at all.
//!
//! # What a reader should look at first
//!
//! `provenance`. A rate that says `Borrowed` was not fitted from this library's reads at all — it
//! is the unweighted mean of the other libraries' rates, handed over because this one stood on
//! fewer than 10,000 sites — and `Defaulted` means not even that was available. Only `FittedHere`
//! is a measurement. Then `railed`, which says the fit stopped at an end of its ladder rather than
//! finding a maximum inside it: the one way this estimator returns a confident wrong number
//! instead of failing.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, SampleLocusObservationsIterator, UnhandledReason,
};
use pop_var_caller::ng::parameter_estimation::generic::accumulators::{
    ConstantPloidy, InbreedingMode,
};
use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::generic::error_rate_ladder;
use pop_var_caller::ng::parameter_estimation::generic::estimate::{
    GenericEstimationConfig, estimate_generic_parameters,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{
    GenomeRegions, RegionKind, TypedRegion, TypedRegionConfig,
};
use pop_var_caller::ng::repeat_catalog::{
    ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria,
};
use pop_var_caller::ng::types::{InbreedingF, Ploidy};
use pop_var_caller::regions::ContigBounds;

#[path = "shared/reference_check.rs"]
mod reference_check_knob;
use reference_check_knob::{reference_check_from_env, reference_check_label};

/// Fit one sample and print a line per read group.
fn fit_one_sample(
    fasta: &Path,
    catalog_path: &Path,
    bed: &Path,
    sample_path: &Path,
    cache: &Arc<ReferenceInfoCache>,
) -> Result<(), Box<dyn std::error::Error>> {
    let (info, verify) = read_reference_verifying_or_creating_fai(
        cache,
        fasta.to_path_buf(),
        reference_check_from_env()?,
    )?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(fasta)?;

    // **The catalog by explicit path, not beside the reference.** Both benchmark references live
    // under a read-only mount, so the file that would normally sit next to the FASTA cannot be
    // written there; the reader takes a path and checks the contig digests either way, so nothing
    // about the check is weakened by moving the file.
    let catalog = RepeatCatalog::open_checking_against_reference(catalog_path, &info)?;

    let bounds: Vec<ContigBounds<'_>> = contigs
        .entries
        .iter()
        .map(|entry| ContigBounds {
            name: &entry.name,
            length: entry.length as u32,
        })
        .collect();
    let spans: Vec<_> = GenomeRegions::from_bed_path(bed, &bounds)?.iter().collect();

    // The pre-pass's own region typing: the catalog says which stretches are repeat tracts and
    // which are ordinary sequence, and only the ordinary ones reach the accumulator.
    let criteria = StrRepeatCriteria::from(&TypedRegionConfig::default());
    let typed: Vec<Result<TypedRegion, RepeatCatalogError>> = catalog
        .genome_segments(&criteria, ReadScope::Regions(&spans))?
        .collect();
    let generic_regions = typed
        .iter()
        .filter(|item| item.as_ref().is_ok_and(|r| r.kind == RegionKind::Generic))
        .count();
    if generic_regions == 0 {
        return Err("the catalog typed no ordinary-sequence region in this BED".into());
    }

    let paths = [sample_path.to_path_buf()];
    let read_groups = build_read_groups(&paths)?;
    let sample_groups = match read_groups.read_groups_per_sample() {
        [only] => only.clone(),
        other => {
            return Err(format!("expected one sample per file, found {}", other.len()).into());
        }
    };

    let reference = OpenReference::new(info);
    let reads = SampleReads::open(
        &sample_groups,
        &read_groups,
        &reference,
        ReadFilterConfig::default(),
        true,
    )?;

    let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
        fasta.to_path_buf(),
        contigs.clone(),
        index.clone(),
    ));
    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "PileupGenerator::new takes an Arc'd accessor; this one is file-backed and this tool is single-threaded, as in the harness it copies"
    )]
    let accessor = Arc::new(WindowedRefSeq::with_shared_index(
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
        accessor,
        make_reference,
        preparer,
        PileupGeneratorConfig::default(),
    )?;
    let generators = GeneratorSet::new(
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        GeneratorSlot::Generator(Box::new(generator)),
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
    );

    let supplied = InbreedingF::try_new(0.0).expect("zero is a fraction");
    let config = GenericEstimationConfig {
        sample_name: sample_groups.sample.to_string(),
        read_groups: sample_groups.read_groups.clone(),
        ploidy: Arc::new(ConstantPloidy(
            Ploidy::try_new(2).expect("a positive copy number"),
        )),
        inbreeding: InbreedingMode::Supplied(supplied),
        fallback_error_rates: BTreeMap::new(),
        edges: Arc::new(DepthBinEdges::new()),
        read_admission: ReadFilterConfig::default(),
    };

    let stream = SampleLocusObservationsIterator::new(typed.into_iter(), reads, generators);
    let parameters = estimate_generic_parameters(stream, &config)?;

    // The ladder ascends in Phred and so descends in error rate: its first rung is the coarsest
    // rate and its last the finest.
    let ladder = error_rate_ladder();
    let coarsest = ladder.first().expect("the ladder has rungs").get();
    let finest = ladder.last().expect("the ladder has rungs").get();
    // **Which of the three numbers is shared differs between the two routes, and that is the
    // whole difficulty of comparing them.** Here the *clean* rate is per read group and the noisy
    // rate and the noisy share are per **sample**; the census route puts clean and noisy both per
    // read group and holds the share per **cohort**. What both can produce is the marginal, so the
    // marginal is what the comparison runs on — and the pieces are printed beside it so the fold
    // can be checked rather than trusted.
    let (noisy, noisy_share) = parameters.site_noise.map_or((f64::NAN, f64::NAN), |noise| {
        (noise.noisy_error_rate().get(), noise.noisy_fraction())
    });

    for (group, rate) in &parameters.error_rate {
        let value = rate.value.get();
        // **The clean rate recovered by inverting the fold, because nothing emits it.** The
        // sample emits the marginal — `(1 - share) * clean + share * noisy` — and the pair behind
        // it, so `clean` follows exactly. At a share of one there is no clean class to recover and
        // at no second class the marginal *is* the clean rate.
        let clean = if noisy_share.is_nan() {
            value
        } else if noisy_share >= 1.0 {
            f64::NAN
        } else {
            (value - noisy_share * noisy) / (1.0 - noisy_share)
        };
        println!(
            "sample={sample}\tread_group={group}\tmarginal_error_rate={value:.6e}\t\
             phred={phred:.2}\tprovenance={provenance:?}\treads={reads}\t\
             railed={railed}\tclean={clean:.6e}\tnoisy={noisy:.6e}\t\
             noisy_share={noisy_share:.4}\tsecond_class_refused={refused}\t\
             coupled_fit_converged={converged}\tinbreeding=supplied_0.0",
            sample = config.sample_name,
            group = group.get(),
            phred = -10.0 * value.log10(),
            provenance = rate.provenance,
            reads = rate.observations,
            // Two independent answers to "did the fit stop at an edge of its search": the fit's
            // own report, and a reconstruction from the ladder's ends. They must agree.
            railed = parameters.error_rate_on_a_ladder_end.contains(group)
                || value >= coarsest
                || value <= finest,
            refused = parameters.site_noise_off_the_ladder,
            converged = parameters.coupled_fit.converged,
        );
    }

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, fasta, catalog, bed, samples @ ..] = args.as_slice() else {
        eprintln!(
            "usage: ng_histogram_error_rates <reference.fa> <catalog.parquet> <regions.bed> \
             <sample.bam|cram> [sample ...]"
        );
        return ExitCode::from(2);
    };
    if samples.is_empty() {
        eprintln!("error: no sample given; at least one BAM or CRAM is needed");
        return ExitCode::from(2);
    }

    let check = match reference_check_from_env() {
        Ok(check) => check,
        Err(error) => {
            eprintln!("error: {error}");
            return ExitCode::from(2);
        }
    };
    println!("reference_check={}", reference_check_label(check));

    let (fasta, catalog, bed) = (
        PathBuf::from(fasta),
        PathBuf::from(catalog),
        PathBuf::from(bed),
    );
    let cache = Arc::new(ReferenceInfoCache::new());
    // **One sample at a time, and a failure names the sample rather than ending the run.** A
    // cohort of 63 is a long walk and one unreadable alignment should not cost the other 62.
    let mut failures = 0_u32;
    for sample in samples {
        let path = PathBuf::from(sample);
        if let Err(error) = fit_one_sample(&fasta, &catalog, &bed, &path, &cache) {
            eprintln!("error: {}: {error}", path.display());
            failures += 1;
        }
    }
    if failures > 0 {
        eprintln!("{failures} sample(s) failed");
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
