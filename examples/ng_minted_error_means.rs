//! **How far apart are the two averages of a read's error probability, on real reads?**
//!
//! ```text
//! PVC_MINTED_ERROR_CENSUS=1 cargo run --release --example ng_minted_error_means -- \
//!     <reference.fa> <regions.bed> <sample.bam|cram> [sample ...]
//! ```
//!
//! # The question
//!
//! The walk mints one error probability per read — the worse of what the instrument said about the
//! bases and what the aligner said about the placement — and then throws the reads away, keeping per
//! allele per library a read count and the sum of those probabilities' logarithms. **So what the
//! caller charges is one number per observation**, `exp(Σ ln ε / n)`, and a single factor per library
//! rescales it so the average comes out at the rate the parameter pre-pass measured. That average is
//! therefore the **geometric** mean
//! (`doc/devel/ng/spec/read_likelihoods.md` §3.2, corrected by the owner on 2026-08-24). The
//! specification first asked for the **arithmetic** mean, `Σ ε / n`, and nothing in the walk
//! carries it: the fold sums logarithms into an observation's `q_sum` and discards the reads.
//!
//! The geometric mean was chosen partly because it is the quantity the model charges, and partly
//! on the expectation that the two sit close together. **Nobody had measured the gap.** This tool
//! measures it, and it is the whole of what it does.
//!
//! Why it might not be small: the minted error is `max(ln ε_BQ, ln ε_MQ)` — the *worse* of the
//! two — and mapping quality spans a wider range than base quality does. An arithmetic mean
//! follows the badly-mapped tail; a geometric mean follows the bulk. So the ratio this prints is
//! a fact about how spread out one library's per-read errors are, and about nothing else: reads
//! that all carry the same error give a ratio of exactly 1.
//!
//! # Where the arithmetic sum comes from, since nothing stores it
//!
//! From [`minted_error_census`](pop_var_caller::ng::locus_generation::pileup::minted_error_census),
//! which sums both shapes at the one place a read's own error still exists as one read's — the
//! record's finalise, before `add_contribution` pools it into an observation. It is off unless
//! `PVC_MINTED_ERROR_CENSUS=1` is set, and **this tool refuses to run without it** rather than
//! printing zeros: `std::env::set_var` is `unsafe` in edition 2024 and this crate forbids
//! `unsafe`, so the variable has to come from the invocation.
//!
//! # The identity this tool also checks
//!
//! The census and the pre-pass's own accumulator
//! ([`calibration`](pop_var_caller::ng::parameter_estimation::generic::calibration)) are supposed
//! to see **the same reads**: complete witnesses at generic loci. This tool folds every locus
//! through the accumulator as well and prints both read counts and both log-error sums, so a run
//! where the two disagree says so in its own output instead of leaving the claim to an argument.
//! `reads_agree=yes` is the line to read.
//!
//! # What this tool's site set is, and how it differs from the pre-pass's
//!
//! Every BED interval is walked as a **generic** region, because the typed-region stream reads a
//! repeat catalog built beside the reference and neither benchmark reference has one (and
//! `$HOME/genomes` is mounted read-only, so this tool cannot write one). So repeat tracts inside
//! the intervals are walked through the generic generator here, where the real pre-pass routes
//! them to the repeat-tract generator and the calibration accumulator's `LocusKind::Generic` gate
//! excludes them. That widens the site set; it does not change what is being compared, because
//! both means are taken over the same reads whichever reads those are.
//!
//! # What it prints, one line per read group
//!
//! Tab-separated `name=value`, so a run can be pasted into a report and the columns picked off by
//! name rather than by position:
//!
//! - `reads_census` and `reads_accumulator` — the same count from the two paths, and
//!   `reads_agree` says whether they match. **A read here is a read at a position**, counted once
//!   for every position it is seen at, which is what the fitted error rate is a rate per.
//! - `geometric_mean`, `arithmetic_mean`, `ratio_arithmetic_over_geometric` — the answer, plus
//!   `phred_geometric` and `phred_arithmetic` for readers who think in Phred.
//! - `geometric_mean_from_accumulator` — the same mean from the pre-pass's fixed-point sum rather
//!   than this tool's `f64` one. It differs from `geometric_mean` only by the accumulator's
//!   documented rounding.
//! - `reads_charged_a_full_unit` and `their_share_of_the_arithmetic_mean` — how many reads carry
//!   an error probability of exactly one (a mate the walk silenced), and how much of `Σ ε` they
//!   are. **This is the column that says whether the arithmetic mean is about the chemistry.**
//! - `reads_under_the_histograms_cap`, `geometric_mean_under_that_cap` and
//!   `cap_moves_the_geometric_mean_by` — what the denominator would be if it thinned each position
//!   the way the error-rate histogram does. Exactly 1.0000 where no position was over the cap.
//! - `loci_never_ruled_on`, once per sample — loci the walk built and the generator neither kept
//!   nor discarded. Zero on every run so far; a non-zero says how many reads are missing from the
//!   answer rather than making it wrong.
//!
//! A final `COHORT` line pools every read group of every sample given.
//!
//! # It is a measurement tool and it ships as one
//!
//! Single-threaded, one sample at a time, and it takes a mutex per read while the census is
//! armed. Nothing in the caller reads it.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::ng::locus_generation::pileup::minted_error_census::{self, MintedErrorTotals};
use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, SampleLocusObservationsIterator, UnhandledReason,
};
use pop_var_caller::ng::parameter_estimation::generic::calibration::{
    MintedReadErrors, fold_into, minted_error_by_read_group,
};
use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion};
use pop_var_caller::ng::repeat_catalog::RepeatCatalogError;
use pop_var_caller::ng::types::ReadGroupId;
use pop_var_caller::regions::ContigBounds;

#[path = "shared/reference_check.rs"]
mod reference_check_knob;
use reference_check_knob::{reference_check_from_env, reference_check_label};

/// One read group's answer: the two means, and the evidence that both were taken over the same
/// reads.
///
/// **Both counts are printed, not just their agreement**, because the interesting failure is not
/// "they differ by one" but "one path saw a tenth of the reads the other did" — and a `yes`/`no`
/// hides the size of that.
struct GroupAnswer {
    name: String,
    /// From the pre-pass's accumulator: Σ `ln ε` in fixed point, and the read count.
    accumulator: MintedReadErrors,
    /// From the census: the same Σ `ln ε` in `f64`, plus Σ `ε`, plus its own read count.
    census: MintedErrorTotals,
    /// **The same fold, but with each site thinned to the depth the error-rate histogram bins
    /// at.** The histogram caps every position at
    /// [`DepthBinEdges::max_depth`](pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges::max_depth)
    /// before the rate is fitted; the calibration fold caps nothing. Per site that is harmless —
    /// the draw is on counts and never on a quality — but it changes **how much weight each site
    /// carries**: a 500-read position gets 500 votes in the denominator and 124 in the population
    /// the numerator was fitted from. This column is what says whether that matters, and it can
    /// only differ where some position was deeper than the cap.
    thinned_to_the_histograms_cap: MintedReadErrors,
}

impl GroupAnswer {
    /// The line one read group gets. Everything on it is measured; nothing is derived from an
    /// assumption about the other columns.
    fn render(&self) -> String {
        let reads = self.census.reads;
        let geometric = self.census.geometric_mean().map_or(f64::NAN, |mean| mean);
        let arithmetic = self.census.arithmetic_mean().map_or(f64::NAN, |mean| mean);
        // The accumulator's own geometric mean, from its fixed-point sum rather than the
        // census's `f64` one. If the two site sets agree this equals `geometric` to the
        // accumulator's documented 2^-21 on the mean log.
        let accumulator_geometric = self
            .accumulator
            .mean_error_probability()
            .map_or(f64::NAN, |mean| mean);
        format!(
            "read_group={name}\treads_census={reads}\treads_accumulator={accumulator_reads}\t\
             reads_agree={agree}\tgeometric_mean={geometric:.6e}\t\
             arithmetic_mean={arithmetic:.6e}\tratio_arithmetic_over_geometric={ratio:.4}\t\
             geometric_mean_from_accumulator={accumulator_geometric:.6e}\t\
             phred_geometric={phred_geometric:.2}\tphred_arithmetic={phred_arithmetic:.2}\t\
             reads_charged_a_full_unit={full_unit}\t\
             their_share_of_the_arithmetic_mean={full_unit_share:.4}\t\
             reads_under_the_histograms_cap={thinned_reads}\t\
             geometric_mean_under_that_cap={thinned_geometric:.6e}\t\
             cap_moves_the_geometric_mean_by={cap_ratio:.4}",
            name = self.name,
            accumulator_reads = self.accumulator.reads(),
            agree = if self.accumulator.reads() == reads {
                "yes"
            } else {
                "NO"
            },
            ratio = arithmetic / geometric,
            // The same two numbers as Phred scores, because that is the scale the reader's
            // qualities are in: 30 is one wrong base in a thousand.
            phred_geometric = -10.0 * geometric.log10(),
            phred_arithmetic = -10.0 * arithmetic.log10(),
            // **What the arithmetic mean is actually made of.** A read the mate-overlap rule
            // silenced carries ε = 1 exactly, so it contributes a whole unit to `Σ ε` and
            // nothing to `Σ ln ε`. If this share is near one, the arithmetic mean is a
            // measurement of how often mates overlap and not of the chemistry.
            full_unit = self.census.reads_charged_a_full_unit,
            full_unit_share = self
                .census
                .full_unit_share_of_arithmetic_mean()
                .unwrap_or(f64::NAN),
            thinned_reads = self.thinned_to_the_histograms_cap.reads(),
            thinned_geometric = self
                .thinned_to_the_histograms_cap
                .mean_error_probability()
                .unwrap_or(f64::NAN),
            // **Above one means the denominator is charging a worse average error than the
            // population the numerator was fitted from**, and the whole of the difference is
            // which sites got how many votes. Exactly one where no position was over the cap.
            cap_ratio = accumulator_geometric
                / self
                    .thinned_to_the_histograms_cap
                    .mean_error_probability()
                    .unwrap_or(f64::NAN),
        )
    }
}

/// Walk one sample and answer for each of its read groups, plus one pooled row.
fn walk_sample(
    fasta: &Path,
    bed: &Path,
    sample_path: &Path,
    cache: &Arc<ReferenceInfoCache>,
) -> Result<Vec<GroupAnswer>, Box<dyn std::error::Error>> {
    // **Reset before, not after.** A run that dies mid-sample then leaves the table dirty for the
    // next one, and the next sample's numbers would be a silent mixture.
    minted_error_census::reset();

    let (info, verify) = read_reference_verifying_or_creating_fai(
        cache,
        fasta.to_path_buf(),
        reference_check_from_env()?,
    )?;
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(fasta)?;

    let paths = [sample_path.to_path_buf()];
    let read_groups = build_read_groups(&paths)?;

    let reference = OpenReference::new(info);
    let sample =
        SampleReads::open_only_sample(&paths, &reference, ReadFilterConfig::default(), true)?;

    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "PileupGenerator::new takes an Arc'd accessor; this one is file-backed and this tool is single-threaded — same as ng_generic_walk_probe"
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
    // The same preparer the walk probe uses, so the reads this tool admits are the reads the
    // real generic walk admits — a different normalizer would change which indels are
    // left-aligned and so which reads witness what.
    let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::with_shared_index(
        fasta.to_path_buf(),
        contigs.clone(),
        index.clone(),
    ));
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

    let bounds: Vec<ContigBounds<'_>> = contigs
        .entries
        .iter()
        .map(|entry| ContigBounds {
            name: &entry.name,
            length: entry.length as u32,
        })
        .collect();
    let spans = GenomeRegions::from_bed_path(bed, &bounds)?;
    // Every interval as a **generic** region — see the module doc for what that costs and why
    // the catalog is not available here.
    // **The error type is the catalog's, because that is the one the locus stream converts
    // from** — this tool's regions come from a BED and every item is `Ok`, so the channel is
    // never used; naming the catalog's type is what lets the stream be built at all.
    let regions: Vec<Result<TypedRegion, RepeatCatalogError>> = spans
        .iter()
        .map(|region| {
            Ok(TypedRegion {
                region,
                kind: RegionKind::Generic,
            })
        })
        .collect();

    let mut accumulated: BTreeMap<ReadGroupId, MintedReadErrors> = BTreeMap::new();
    let mut thinned: BTreeMap<ReadGroupId, MintedReadErrors> = BTreeMap::new();
    let mut scratch: Vec<(ReadGroupId, MintedReadErrors)> = Vec::new();
    let mut thinned_scratch: Vec<(ReadGroupId, MintedReadErrors)> = Vec::new();
    // **The cap read off the ladder the pre-pass actually bins with**, not written down here: a
    // number copied into a measuring tool is a number that can drift away from the thing it
    // measures.
    let cap = u64::from(DepthBinEdges::new().max_depth());
    let mut stream = SampleLocusObservationsIterator::new(regions.into_iter(), sample, generators);
    for locus in &mut stream {
        let locus = locus?;
        // Exactly what `GenericAccumulators::add_locus` does with the locus, and nothing else —
        // so the numbers below are the accumulator's own and not a re-derivation of them.
        minted_error_by_read_group(&locus, &mut scratch);
        fold_into(&mut accumulated, &scratch);
        // And the same locus again with each read group's reads thinned to the cap, at that
        // group's own mean for this site. **The expectation of the real draw, not the draw**:
        // the histogram's thinning is hypergeometric on counts, so the reads it keeps have the
        // site's mean quality in expectation, and using that mean removes the sampling noise
        // rather than adding a second source of it.
        thinned_scratch.clear();
        thinned_scratch.extend(scratch.iter().map(|&(group, totals)| {
            let kept = totals.reads().min(cap);
            let mean = totals.mean_log_error().unwrap_or(0.0);
            (
                group,
                MintedReadErrors::of_observation(
                    mean * kept as f64,
                    u32::try_from(kept).unwrap_or(u32::MAX),
                ),
            )
        }));
        fold_into(&mut thinned, &thinned_scratch);
        drop(locus);
    }

    let census: BTreeMap<ReadGroupId, MintedErrorTotals> =
        minted_error_census::snapshot().into_iter().collect();
    // Loci the walk built and the generator never ruled on. Printed rather than asserted: a
    // non-zero does not make the answer wrong, it says how many reads are missing from it.
    println!(
        "sample={}\tloci_never_ruled_on={}",
        sample_path.display(),
        minted_error_census::loci_never_ruled_on(),
    );

    // Every read group either path saw, so a group present in one and missing from the other is
    // a row rather than a silence.
    let mut every_group: Vec<ReadGroupId> = accumulated.keys().copied().collect();
    for group in census.keys() {
        if !every_group.contains(group) {
            every_group.push(*group);
        }
    }
    every_group.sort_unstable();

    let mut answers: Vec<GroupAnswer> = every_group
        .iter()
        .map(|&group| GroupAnswer {
            name: format!(
                "{}/{}",
                read_groups.get(group).sample,
                read_groups.get(group).id
            ),
            accumulator: accumulated.get(&group).copied().unwrap_or_default(),
            census: census.get(&group).copied().unwrap_or_default(),
            thinned_to_the_histograms_cap: thinned.get(&group).copied().unwrap_or_default(),
        })
        .collect();

    // **The pooled row is the sample's own, not the cohort's.** Read groups within one library
    // preparation are the grain the scale is per; pooling them says what the sample as a whole
    // looks like, which is what makes a single-read-group sample and a many-group one comparable
    // on one line.
    if answers.len() > 1 {
        let mut pooled_accumulator = MintedReadErrors::default();
        let mut pooled_thinned = MintedReadErrors::default();
        let mut pooled_census = MintedErrorTotals::default();
        for answer in &answers {
            pooled_accumulator.add(answer.accumulator);
            pooled_thinned.add(answer.thinned_to_the_histograms_cap);
            pooled_census.log_error_sum += answer.census.log_error_sum;
            pooled_census.error_sum += answer.census.error_sum;
            pooled_census.reads += answer.census.reads;
            pooled_census.reads_charged_a_full_unit += answer.census.reads_charged_a_full_unit;
        }
        answers.push(GroupAnswer {
            name: format!("{}/ALL", read_groups.get(ReadGroupId(0)).sample),
            accumulator: pooled_accumulator,
            census: pooled_census,
            thinned_to_the_histograms_cap: pooled_thinned,
        });
    }

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(answers)
}

fn main() -> ExitCode {
    if !minted_error_census::enabled() {
        eprintln!(
            "error: the per-read census is off, so the arithmetic mean would come back as a \
             count of zero reads. Re-run with PVC_MINTED_ERROR_CENSUS=1 in the environment — \
             this tool cannot set it for you, because std::env::set_var is unsafe in edition \
             2024 and this crate forbids unsafe."
        );
        return ExitCode::from(2);
    }

    let args: Vec<String> = std::env::args().collect();
    let [_, fasta, bed, samples @ ..] = args.as_slice() else {
        eprintln!(
            "usage: ng_minted_error_means <reference.fa> <regions.bed> <sample.bam|cram> \
             [sample ...]"
        );
        return ExitCode::from(2);
    };
    if samples.is_empty() {
        eprintln!("error: no sample given; at least one BAM or CRAM is needed");
        return ExitCode::from(2);
    }

    let fasta = PathBuf::from(fasta);
    let bed = PathBuf::from(bed);
    let cache = Arc::new(ReferenceInfoCache::new());
    println!(
        "reference_check={}",
        reference_check_label(match reference_check_from_env() {
            Ok(check) => check,
            Err(error) => {
                eprintln!("error: {error}");
                return ExitCode::from(2);
            }
        })
    );

    // **Everything pooled across samples too**, because the question was asked about a cohort
    // and one accession's answer is one accession's.
    let mut cohort = MintedErrorTotals::default();
    for sample in samples {
        let path = PathBuf::from(sample);
        match walk_sample(&fasta, &bed, &path, &cache) {
            Ok(answers) => {
                for answer in &answers {
                    // The pooled-within-sample row is printed but not added twice.
                    if !answer.name.ends_with("/ALL") {
                        cohort.log_error_sum += answer.census.log_error_sum;
                        cohort.error_sum += answer.census.error_sum;
                        cohort.reads += answer.census.reads;
                        cohort.reads_charged_a_full_unit += answer.census.reads_charged_a_full_unit;
                    }
                    println!("{}", answer.render());
                }
            }
            Err(error) => {
                eprintln!("error: {}: {error}", path.display());
                return ExitCode::FAILURE;
            }
        }
    }

    let geometric = cohort.geometric_mean().unwrap_or(f64::NAN);
    let arithmetic = cohort.arithmetic_mean().unwrap_or(f64::NAN);
    println!(
        "COHORT\treads={reads}\tgeometric_mean={geometric:.6e}\tarithmetic_mean={arithmetic:.6e}\t\
         ratio_arithmetic_over_geometric={ratio:.4}\tphred_geometric={phred_geometric:.2}\t\
         phred_arithmetic={phred_arithmetic:.2}\treads_charged_a_full_unit={full_unit}\t\
         their_share_of_the_arithmetic_mean={full_unit_share:.4}",
        reads = cohort.reads,
        ratio = arithmetic / geometric,
        phred_geometric = -10.0 * geometric.log10(),
        phred_arithmetic = -10.0 * arithmetic.log10(),
        full_unit = cohort.reads_charged_a_full_unit,
        full_unit_share = cohort
            .full_unit_share_of_arithmetic_mean()
            .unwrap_or(f64::NAN),
    );
    ExitCode::SUCCESS
}
