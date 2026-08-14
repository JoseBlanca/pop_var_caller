//! What the joint route's records actually cost, and what a cohort of them says.
//!
//! Everything measured for this route so far was measured against truths the measuring
//! program made up itself. This one drives the real writer —
//! `parameter_estimation::joint::census::CensusWriter` — through real alignments, and it
//! answers two questions the specifications assert without evidence:
//!
//! 1. **What does the census weigh?** `parameter_prepass_joint_records.md` §6 prices it by
//!    arithmetic and §7.8 says to measure instead. The generic depth array, the sparse list of
//!    non-reference observations and the STR set are reported **separately**, because a single
//!    total would hide any one of them being wrong.
//! 2. **How much population structure does the cohort have?** Handed more than one alignment,
//!    the program holds every sample's records at once — which is the state the fit runs in —
//!    and asks how far apart the samples are. That decides which row of
//!    `parameter_prepass_joint_fit.md` §3.4.2 a cohort sits in, and so how much the
//!    contamination work matters.
//!
//! # The structure measurement, and its control
//!
//! Allele frequencies estimated from a panel are pulled towards the panel's average, and how
//! badly depends on how diverged the panel is. On a panel with no divergence a contamination
//! estimate comes back exact; at `F_st` 0.20 a genuinely 3%-contaminated accession comes back
//! at 0.5% and passes as clean (`joint_contamination_2026-08-12.md`). So the number wanted
//! here is the cohort's own divergence.
//!
//! **Divergence is measured against a null that runs the identical pipeline.** Each sample's
//! reads at a position are re-drawn from a binomial at the cohort's own allele frequency and
//! that sample's own depth — which destroys every trace of structure and keeps the depths,
//! the missingness and the frequencies exactly. The same principal-component decomposition
//! and the same split-and-measure then run on the re-drawn data. Whatever the null returns is
//! what read sampling alone manufactures at this depth and this sample count; the excess over
//! it is the structure.
//!
//! Without that control the measurement has the shape of the failure recorded in
//! `joint_route_research_narrative_2026-08-13.md` §7: a number that looks like an answer and
//! is an artefact of the estimator.
//!
//! ```text
//! ng_joint_records_walk <reference.fa> <catalog.parquet> <regions.bed> <generic-target> \
//!     <alignment> [alignment...]
//! ```

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use pop_var_caller::ng::locus_generation::ssr::{SsrGenerator, SsrGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, LocusKind, SampleLocusObservationsIterator, UnhandledReason,
};
use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    CensusWriter, DepthCap, DepthCode, ReadCap, SampleCensusEvidence,
};
use pop_var_caller::ng::parameter_estimation::joint::contamination::{
    ContaminationConfig, ContaminationEstimate, OwnCoordinates, fit_contamination,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{JointFit, JointFitConfig, fit_jointly};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLoci, ReferenceDigest, RegionSetDigest, SelectableRegions,
    SelectionTerms, UnambiguousRuns, select_kept_loci,
};
use pop_var_caller::ng::parameter_estimation::joint::ssr_fit;
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceInfo, ReferenceSource, read_reference_info_observing,
};
use pop_var_caller::ng::region_typing::{
    GenomeRegions, RegionKind, TypedRegion, TypedRegionConfig,
};
use pop_var_caller::ng::repeat_catalog::{
    ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria,
};
use pop_var_caller::ng::types::{Bp, ContigId};
use pop_var_caller::regions::ContigBounds;

/// The seed the locus selection is drawn with. One number, shared by every sample, and part
/// of what the fit refuses to pool across.
const SEED: u64 = 42;

/// The depth above which a position's reads stop being counted one by one.
///
/// **Held at 124 across the ladder's extension on purpose.** It used to be read off the
/// ladder's top rung, which made one number of two: the ladder now reaches about 1,500 and
/// this does not follow it, because what the cap buys is a bound on the sparse list of
/// non-reference reads and what the ladder buys is reach for the position's own depth.
const DEPTH_CAP: DepthCap = DepthCap::new(124);

fn main() {
    let usage = "usage: <reference.fa> <catalog.parquet> <regions.bed> <generic-target> \
                 <alignment> [alignment...]";
    let mut args = std::env::args().skip(1);
    let fasta = PathBuf::from(args.next().expect(usage));
    let catalog_path = PathBuf::from(args.next().expect(usage));
    let bed = PathBuf::from(args.next().expect(usage));
    let generic_target: u64 = args.next().expect(usage).parse().expect("a position count");
    let alignments: Vec<PathBuf> = args.map(PathBuf::from).collect();
    assert!(!alignments.is_empty(), "{usage}");

    let started = Instant::now();

    // ---- the reference, and where it is sequence at all --------------------------------
    let mut callable = UnambiguousRuns::default();
    let info = Arc::new(
        read_reference_info_observing(
            ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            },
            &mut callable,
        )
        .expect("the reference reads"),
    );
    let contigs = Arc::new(info.contig_list());
    let index = WindowedRefSeq::read_index(&fasta).expect("the .fai beside the reference reads");
    let unambiguous = callable
        .into_selectable()
        .expect("maximal runs of A/C/G/T are disjoint");
    println!(
        "reference        {} — {} contigs, read in {:.1} s",
        fasta.display(),
        contigs.entries.len(),
        started.elapsed().as_secs_f64()
    );

    // ---- the analysed regions, and the catalog ------------------------------------------
    let bounds: Vec<ContigBounds> = contigs
        .entries
        .iter()
        .map(|entry| ContigBounds {
            name: &entry.name,
            length: u32::try_from(entry.length).expect("a contig shorter than 4 Gb"),
        })
        .collect();
    let analysed = SelectableRegions::new(
        GenomeRegions::from_bed_path(&bed, &bounds)
            .expect("the BED resolves against this reference's contigs")
            .iter()
            .collect(),
    )
    .expect("a BED's merged spans are disjoint");
    let catalog = RepeatCatalog::open_checking_against_reference(&catalog_path, &info)
        .expect("the catalog is this reference's");
    let criteria = StrRepeatCriteria::from(&TypedRegionConfig::default());

    let terms = SelectionTerms {
        seed: SEED,
        reference: ReferenceDigest::of(&info).expect("the reference has a digest"),
        analysed_regions: RegionSetDigest::of(&analysed),
        catalog_built_under: CatalogBuildSettings::of(&catalog),
        ssr_criteria: criteria.clone(),
        generic_target,
        // No stratum is capped here: the tomato catalog's largest holds far fewer loci than
        // this, so the cap never fires and every STR locus is kept.
        ssr_cap: 1_000_000,
    };

    let kept = select_kept_loci(&terms, &catalog, &analysed, &unambiguous)
        .expect("the catalog serves these criteria");
    println!(
        "analysed regions {} — {} spans, {} bases",
        bed.display(),
        analysed.spans().len(),
        analysed.total_length().get()
    );
    println!(
        "kept loci        {} generic positions (target {generic_target}), {} STR loci in {} \
         strata, {:.1} s",
        kept.generic().len(),
        kept.ssr_stratum_counts().total(),
        kept.ssr_stratum_counts().iter_sorted().len(),
        started.elapsed().as_secs_f64()
    );

    // The walk and the selection run over the same stretch of genome: the analysed regions
    // with the ambiguous runs cut out.
    let masked = analysed.intersect(&unambiguous);
    let typed: Vec<TypedRegion> = catalog
        .genome_segments(&criteria, ReadScope::Regions(masked.spans()))
        .expect("the BED's spans name contigs this catalog holds")
        .map(|item| item.expect("the catalog reads through the whole of the BED"))
        .collect();
    let generic_domain = SelectableRegions::new(
        typed
            .iter()
            .filter(|region| region.kind == RegionKind::Generic)
            .map(|region| region.region)
            .collect(),
    )
    .expect("typed regions are disjoint");
    // ---- one walk per sample -------------------------------------------------------------
    let mut cohort: Vec<SampleCensusEvidence> = Vec::new();
    for alignment in &alignments {
        let at = Instant::now();
        let records = walk_one(
            &fasta,
            &info,
            &contigs,
            &index,
            alignment,
            &typed,
            &generic_domain,
            &kept,
            &terms,
        );
        println!(
            "\n=== {} — {:.1} s",
            records.sample,
            at.elapsed().as_secs_f64()
        );
        report_sizes(&records);
        cohort.push(records);
    }

    // ---- do they pool at all? --------------------------------------------------------------
    println!("\n--- the recording-terms check, which is what lets these be pooled ---");
    let first = &cohort[0];
    let mut refused = 0;
    for other in &cohort[1..] {
        if let Some(field) = first.terms.first_disagreement(&other.terms) {
            println!(
                "  {} disagrees with {} on {field}",
                other.sample, first.sample
            );
            refused += 1;
        }
    }
    println!(
        "  {} of {} samples agree on all twelve values",
        cohort.len() - refused,
        cohort.len()
    );

    if cohort.len() >= 8 {
        structure(&cohort, kept.generic().len());
    } else {
        println!(
            "\n--- structure: skipped, {} samples is too few to decompose ---",
            cohort.len()
        );
    }

    // The kept positions themselves, where a run asks for them. **This is what makes a
    // comparison against a benchmark VCF exact**: the fit's rates are means over these
    // positions and no others, so the truth has to be counted over the same list.
    if let Ok(path) = std::env::var("KEPT_POSITIONS_TSV") {
        let mut out = String::new();
        for position in kept.generic() {
            out.push_str(&format!(
                "{}\t{}\n",
                contigs.entries[position.contig.get() as usize].name,
                position.position.get()
            ));
        }
        std::fs::write(&path, out).expect("the kept-positions path is writable");
        println!("\nkept positions   {path}");
    }

    // **The duplicated class's depth term has no supplier here any more.** It used to come
    // from a per-sample coverage-by-window summary built during this walk; that summary is
    // gone (`parameter_prepass_joint_records.md` §4), and what replaces it is the position's
    // own depth, which the census already carries. Until the fit reads that for itself the
    // class is left with the cohort pattern alone — a position where the carriers all read
    // near a half and nobody is homozygous for the non-reference allele — which needs about
    // twenty-five samples before its absence means anything.
    let coverage_odds = Vec::new();

    depth_ladder_occupancy(&cohort);

    // **The ordinary-position half runs once and its answer is handed on.** The repeat-tract
    // half needs each sample's homozygote excess from it and gives nothing back, so re-fitting
    // to obtain that would be 883 seconds of the tomato cohort's time spent twice.
    if let Some(fit) = fit_the_cohort(&cohort, coverage_odds, kept.generic(), &contigs) {
        fit_the_tracts(&cohort, &kept, &contigs, &fit);
    }

    println!(
        "\ntotal            {:.1} s",
        started.elapsed().as_secs_f64()
    );
}

/// Where each sample's positions sit on the depth ladder.
///
/// **The ladder is exact to eight reads a position and a range above it**, and how much of a
/// run is in the exact part decides how much anything about the range can matter. A cohort at
/// three reads a position is almost entirely exact; the band above is where a stored code
/// stands for several depths at once.
///
/// **The last column is the one to watch, and it is about the cap rather than the ladder.**
/// The depth recorded is now the position's own, so the ladder no longer piles deep positions
/// onto one rung; what still collapses is the *allele counts*, which are thinned to the
/// per-position cap. A position above that cap has an exact depth and a thinned count beside
/// it, and the fraction it reports is the thinned one.
fn depth_ladder_occupancy(cohort: &[SampleCensusEvidence]) {
    let edges = DepthBinEdges::for_census();
    println!("\n--- where the positions sit on the depth ladder ---");
    println!(
        "  {:<24}{:>14}{:>18}{:>16}{:>12}",
        "sample", "exact (0–8)", "a range (9+)", "above the cap", "no reads"
    );
    for sample in cohort {
        let cap = sample.terms.depth_cap.get();
        let (mut exact, mut ranged, mut capped, mut unwalked) = (0_u64, 0_u64, 0_u64, 0_u64);
        for records in sample.generic.values() {
            for code in records.depth().iter() {
                match code {
                    DepthCode::NeverWalked => unwalked += 1,
                    DepthCode::Binned(bin) => {
                        let range = edges.depth_range(bin);
                        if *range.start() > cap {
                            capped += 1;
                        } else if range.start() == range.end() {
                            exact += 1;
                        } else {
                            ranged += 1;
                        }
                    }
                }
            }
        }
        let total = (exact + ranged + capped).max(1) as f64;
        println!(
            "  {:<24}{:>13.1}%{:>17.1}%{:>13.1}%{:>12}",
            sample.sample,
            100.0 * exact as f64 / total,
            100.0 * ranged as f64 / total,
            100.0 * capped as f64 / total,
            unwalked
        );
    }
}

/// Fit every parameter from the records just built, and print what came back.
///
/// **These are real reads, so there is no truth here** beyond what a benchmark VCF supplies
/// separately. What the run shows is the shape of the answer — whether the two classes of
/// position stay apart, where the population's frequency density lands, and what each sample's
/// heterozygosity comes out at with the count of positions that actually carried a read beside
/// it.
fn fit_the_cohort(
    cohort: &[SampleCensusEvidence],
    coverage_odds: Vec<Arc<[f32]>>,
    kept: &[pop_var_caller::ng::types::GenomePosition],
    contigs: &ContigList,
) -> Option<JointFit> {
    println!("\n--- the joint fit, ordinary positions ---");
    let at = Instant::now();
    let genotype_posteriors = std::env::var("GENOTYPE_POSTERIORS_TSV").ok();
    let config = JointFitConfig {
        quadrature_nodes: 12,
        coverage_odds,
        // The shipped default, unless the run says otherwise — a run comparing contamination
        // against an earlier one wants the rest of the fit held still.
        duplicated_positions: std::env::var("DUPLICATED_CLASS").as_deref() != Ok("0"),
        genotype_posteriors: genotype_posteriors.is_some(),
        pass_trace: true,
        // `DEPTH_AS_A_RANGE=0` restores the point-read, so the same reads can be fitted both
        edges: std::sync::Arc::new(DepthBinEdges::for_census()),
        // ways and the difference attributed.
        depth_as_a_range: std::env::var("DEPTH_AS_A_RANGE").as_deref() != Ok("0"),
        max_passes: std::env::var("MAX_PASSES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(200),
        ..JointFitConfig::default()
    };
    let fit = match fit_jointly(cohort, &config) {
        Ok(fit) => fit,
        Err(error) => {
            println!("  refused: {error}");
            return None;
        }
    };
    println!(
        "  {} passes, {}, log-likelihood {:.0}, {:.1} s",
        fit.passes,
        if fit.converged {
            "converged"
        } else {
            "RAN OUT OF PASSES"
        },
        fit.log_likelihood,
        at.elapsed().as_secs_f64()
    );
    for (group, estimate) in &fit.noise {
        println!(
            "  read group {group:?}: a read misreads at {:.5} at an ordinary position and \
             {:.4} at a mismapped one",
            estimate.value.clean, estimate.value.noisy
        );
    }
    println!(
        "  positions mismapped {:.4}; carrying only the reference {:.4}, only a non-reference \
         base {:.5}; the rest segregate with Beta({:.3}, {:.3})",
        fit.noisy_share,
        fit.density.value.p_invariant,
        fit.density.value.p_fixed_alt,
        fit.density.value.a,
        fit.density.value.b
    );
    match &fit.duplicated {
        Some(duplicated) => println!(
            "  positions a sample carries an extra copy of {:.5}; the share of the panel carrying \
             one is Beta({:.3}, {:.3})",
            duplicated.value.share, duplicated.value.carrier_a, duplicated.value.carrier_b
        ),
        None => println!("  the duplicated class was not fitted"),
    }
    println!(
        "  the population's expected heterozygosity {:.5} ({:.3} per kilobase)",
        fit.expected_heterozygosity,
        1_000.0 * fit.expected_heterozygosity
    );

    // **What the alternation was doing when it stopped.** A fit that ran out of passes and one
    // that settled report the same numbers, and only the trajectory tells them apart: a
    // heterozygosity still climbing at the last pass is not the estimator's answer, it is where
    // the estimator had got to.
    if !fit.trace.is_empty() {
        println!(
            "\n  {:<8}{:>16}{:>14}{:>22}{:>12}{:>10}",
            "pass",
            "log-likelihood",
            "largest move",
            "het/kb, first sample",
            "less het by",
            "Beta a"
        );
        let last = fit.trace.len();
        for entry in &fit.trace {
            let shown = entry.pass <= 5
                || entry.pass % 25 == 0
                || entry.pass as usize == last
                || entry.pass == 29;
            if !shown {
                continue;
            }
            println!(
                "  {:<8}{:>16.0}{:>14.6}{:>22.4}{:>12.3}{:>10.3}",
                entry.pass,
                entry.log_likelihood,
                entry.largest_move,
                1_000.0 * entry.heterozygous.first().copied().unwrap_or(0.0),
                entry.hom_excess.first().copied().unwrap_or(0.0),
                entry.density_a
            );
        }
    }

    // Every sample's posterior at every kept position, for a comparison against a benchmark
    // VCF that can ask *which* positions the two disagree at rather than only by how much.
    if let Some(path) = &genotype_posteriors {
        // The reads the fit saw, beside what it made of them: a posterior is not evidence, and
        // a position called heterozygous at a read share nothing like a half is a different
        // finding from one called heterozygous at a half.
        let edges = DepthBinEdges::for_census();
        let non_reference: Vec<Vec<u32>> = cohort
            .iter()
            .map(|sample| {
                let mut counts = vec![0_u32; kept.len()];
                for group in sample.generic.values() {
                    for observation in group.non_reference() {
                        counts[observation.index as usize] += u32::from(observation.reads);
                    }
                }
                counts
            })
            .collect();

        let mut out = String::new();
        out.push_str("contig\tposition");
        for sample in cohort {
            out.push_str(&format!(
                "\t{0}_het\t{0}_homalt\t{0}_depth\t{0}_nonref",
                sample.sample
            ));
        }
        out.push('\n');
        let width = cohort.len() * 2;
        for (index, position) in kept.iter().enumerate() {
            let row = &fit.genotype_posterior[index * width..][..width];
            out.push_str(&contigs.entries[position.contig.get() as usize].name);
            out.push('\t');
            out.push_str(&position.position.get().to_string());
            for (s, sample) in cohort.iter().enumerate() {
                let mut depth = 0_u32;
                for group in sample.generic.values() {
                    if let DepthCode::Binned(bin) = group.depth().get(index) {
                        let range = edges.depth_range(bin);
                        depth += (*range.start() + *range.end()) / 2;
                    }
                }
                out.push_str(&format!(
                    "\t{:.4}\t{:.4}\t{}\t{}",
                    row[s * 2],
                    row[s * 2 + 1],
                    depth,
                    non_reference[s][index]
                ));
            }
            out.push('\n');
        }
        std::fs::write(path, out).expect("the genotype-posterior path is writable");
        println!("\n  genotype posteriors  {path}");
    }
    println!(
        "\n  {:<26}{:>10}{:>12}{:>14}{:>14}",
        "sample", "het/kb", "hom-alt/kb", "less het by", "positions read"
    );
    let mut names: Vec<&String> = fit.rates.keys().collect();
    names.sort();
    for name in names {
        let rates = &fit.rates[name].value;
        println!(
            "  {:<26}{:>10.3}{:>12.3}{:>14.3}{:>13.1}%",
            name,
            1_000.0 * rates.heterozygous,
            1_000.0 * rates.homozygous_alt,
            fit.hom_excess[name].value.get(),
            100.0 * rates.positions_with_reads as f64
                / cohort[0]
                    .generic
                    .values()
                    .next()
                    .map_or(1.0, |g| g.depth().len() as f64)
        );
    }

    // **Judged against the cohort's own spread, not against a constant.** What a threshold has
    // to clear is the noise floor on the clean samples, and that floor is a property of this
    // panel's depth, marker count and structure rather than a number anyone can quote in
    // advance — so the panel's own distribution is printed beside each value.
    println!("\n  contamination — the share of a sample's reads from another plant");
    let mut estimated: Vec<(&String, f64, f64)> = Vec::new();
    let mut refused = 0;
    let mut markers = 0;
    for (name, estimate) in &fit.contamination {
        match estimate {
            ContaminationEstimate::Estimated {
                alpha,
                markers: count,
                leverage,
            } => {
                markers = *count;
                estimated.push((name, *alpha, *leverage));
            }
            ContaminationEstimate::NotIdentified { reason } => {
                println!("    {name:<24} not identified — {reason}");
                refused += 1;
            }
        }
    }
    estimated.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("no NaN"));
    println!(
        "    {} of {} samples estimated over {markers} varying positions, {refused} refused",
        estimated.len(),
        fit.contamination.len()
    );
    if !estimated.is_empty() {
        let values: Vec<f64> = estimated.iter().map(|e| e.1).collect();
        let median = values[values.len() / 2];
        println!("    the panel's own spread: median {median:.4}, and these are the highest");
        for (name, alpha, leverage) in estimated.iter().take(8) {
            println!("      {name:<24}{alpha:>8.4}   supplies {leverage:.3} of its own frequency");
        }
    }

    // **What excluding the mismapped positions is worth, on this cohort.** A position where two
    // stretches of genome the reference holds once both pile reads up puts a small share of
    // unexpected reads into every sample at once, which is the contamination signature exactly.
    // The fit says how likely each position is to be one; the arms below are that judgement
    // used and ignored, crossed with the three readings of where a sample stands on the panel's
    // axes while its own fraction is searched for. The drawn control that says which arm to
    // believe is `examples/ng_joint_contamination_control.rs`.
    let condemned = fit.noisy_posterior.iter().filter(|p| **p > 0.5).count();
    println!(
        "\n  {condemned} of {} positions are more likely mismapped than not",
        fit.noisy_posterior.len()
    );
    let error: Vec<f64> = cohort
        .iter()
        .map(|sample| {
            let group = *sample.generic.keys().next().expect("a read group");
            fit.noise[&group].value.clean
        })
        .collect();
    let excess: Vec<f64> = cohort
        .iter()
        .map(|sample| fit.hom_excess[&sample.sample].value.get())
        .collect();
    println!(
        "\n  {:<44}{:>10}{:>10}{:>10}{:>10}",
        "", "markers", "median", "highest", "lowest"
    );
    for keep_mismapped in [true, false] {
        for own in [
            OwnCoordinates::AsRead,
            OwnCoordinates::UndoneByAlpha,
            OwnCoordinates::MaximisedFreely,
        ] {
            let at = Instant::now();
            let settings = ContaminationConfig {
                max_noisy_posterior: if keep_mismapped { 1.0 } else { 0.5 },
                weight_by_posterior: !keep_mismapped,
                own_coordinates: own,
                ..ContaminationConfig::default()
            };
            let estimates = fit_contamination(
                cohort,
                &config.edges,
                &error,
                &excess,
                &fit.noisy_posterior,
                &settings,
            );
            let mut values: Vec<f64> = estimates
                .iter()
                .filter_map(|estimate| match estimate {
                    ContaminationEstimate::Estimated { alpha, .. } => Some(*alpha),
                    ContaminationEstimate::NotIdentified { .. } => None,
                })
                .collect();
            values.sort_by(f64::total_cmp);
            if values.is_empty() {
                continue;
            }
            let markers = estimates
                .iter()
                .find_map(|estimate| match estimate {
                    ContaminationEstimate::Estimated { markers, .. } => Some(*markers),
                    ContaminationEstimate::NotIdentified { .. } => None,
                })
                .unwrap_or(0);
            println!(
                "  {:<44}{markers:>10}{:>10.4}{:>10.4}{:>10.4}   {:.0} s",
                format!(
                    "{}, {}",
                    if keep_mismapped {
                        "every position"
                    } else {
                        "mismapped dropped"
                    },
                    match own {
                        OwnCoordinates::AsRead => "coordinates as read",
                        OwnCoordinates::UndoneByAlpha => "coordinates undone by α",
                        OwnCoordinates::MaximisedFreely => "coordinates free",
                    }
                ),
                values[values.len() / 2],
                values.last().copied().unwrap_or(0.0),
                values.first().copied().unwrap_or(0.0),
                at.elapsed().as_secs_f64(),
            );
        }
    }

    Some(fit)
}

/// Fit the repeat tracts, stratum by stratum, from the same records.
///
/// **This half runs after the ordinary-position one and takes one thing from it**: each
/// sample's homozygote excess, which weights the genotype drawn from a tract's own length
/// frequencies. It gives nothing back.
///
/// `SLIPPAGE_PER_READ_GROUP=1` fits the three slippage numbers separately for every read group,
/// which is the grain spec §4 names. The default pools every read group into one set, because
/// a cohort of sixty-three single-read-group samples would otherwise ask 189 numbers of a
/// stratum holding a few dozen tracts.
fn fit_the_tracts(
    cohort: &[SampleCensusEvidence],
    kept: &CensusLoci,
    contigs: &ContigList,
    fit: &JointFit,
) {
    println!("\n--- the joint fit, repeat tracts ---");

    let by_name: BTreeMap<&str, usize> = contigs
        .entries
        .iter()
        .enumerate()
        .map(|(index, entry)| (entry.name.as_str(), index))
        .collect();
    let contig_of = |name: &str| {
        by_name
            .get(name)
            .map(|index| ContigId(u32::try_from(*index).expect("a contig index")))
    };
    let strata = ssr_fit::strata_of_kept_loci(kept, &contig_of);

    // One slippage group per read group is the specified grain; pooling them is what a cohort
    // this thin can afford, and which was used is printed rather than assumed.
    let per_read_group = std::env::var("SLIPPAGE_PER_READ_GROUP").is_ok_and(|v| v == "1");
    let mut slippage_group_of: BTreeMap<_, u32> = BTreeMap::new();
    for sample in cohort {
        for group in sample.ssr.keys() {
            let next = if per_read_group {
                slippage_group_of.len() as u32
            } else {
                0
            };
            slippage_group_of.entry(*group).or_insert(next);
        }
    }
    println!(
        "  {} read groups in {} slippage group{}",
        slippage_group_of.len(),
        slippage_group_of
            .values()
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        if per_read_group {
            "s, one each"
        } else {
            ", pooled"
        }
    );

    let homozygote_excess: Vec<f64> = cohort
        .iter()
        .map(|sample| {
            fit.hom_excess
                .get(&sample.sample)
                .map_or(0.0, |estimate| estimate.value.get())
        })
        .collect();

    let at = Instant::now();
    let evidence = ssr_fit::gather_strata(cohort, &strata, &slippage_group_of);
    println!(
        "  {} strata over {} tracts, gathered in {:.1} s",
        evidence.len(),
        strata.len(),
        at.elapsed().as_secs_f64()
    );
    let over_guard: u64 = evidence.iter().map(|s| s.tracts_over_guard_threshold).sum();
    let not_crossed: u64 = evidence.iter().map(|s| s.reads_reaching_not_crossing).sum();
    let guard_reads: u64 = evidence.iter().map(|s| s.guard_reads).sum();
    let spanning: u64 = evidence.iter().map(|s| s.spanning_reads()).sum();
    println!(
        "  {spanning} reads crossed a tract; {not_crossed} reached one without crossing it; \
         {guard_reads} differed by a non-whole number of repeats; {over_guard} tracts were over \
         the guard's threshold and left out"
    );

    // The lengths the fit may place allele mass on. The specification's ±6 is the default; a
    // narrower span costs a third of the run's time and is what a first pass over a cohort can
    // afford, so which was used is printed rather than left to be guessed.
    let allele_span: i32 = std::env::var("SSR_ALLELE_SPAN")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(ssr_fit::ALLELE_SPAN);
    // `SSR_BORROWING_FLOOR=0` fits every stratum on its own tracts alone, which is the arm that
    // says what borrowing changed rather than only what it produced.
    let borrowing_floor: usize = std::env::var("SSR_BORROWING_FLOOR")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(ssr_fit::DEFAULT_BORROWING_FLOOR);
    let config = ssr_fit::SsrFitConfig {
        allele_span,
        borrowing_floor,
        ..ssr_fit::SsrFitConfig::default()
    };
    println!("  allele mass may sit within ±{allele_span} repeats of the reference length");
    let at = Instant::now();
    let outcomes = ssr_fit::fit_strata(&evidence, &homozygote_excess, &config);
    println!(
        "  fitted in {:.1} s, borrowing below {} tracts and refusing below {}",
        at.elapsed().as_secs_f64(),
        config.borrowing_floor,
        config.refusal_floor
    );

    let held: BTreeMap<_, &ssr_fit::StratumEvidence> =
        evidence.iter().map(|e| (e.stratum, e)).collect();
    println!(
        "\n  {:<8}{:>8}{:>10}{:>9}{:>9}{:>9}{:>9}{:>8}  borrowed from",
        "motif", "repeats", "tracts", "reads", "level", "shorter", "fall-off", "conc.",
    );
    for outcome in &outcomes {
        match outcome {
            ssr_fit::StratumOutcome::Fitted(fitted) => {
                let own = held[&fitted.stratum];
                // One set of slippage numbers a group; with the groups pooled there is one,
                // and with them apart the first fitted group stands for the row.
                let slippage = fitted
                    .slippage
                    .iter()
                    .flatten()
                    .next()
                    .expect("a fitted stratum has at least one group with reads");
                println!(
                    "  {:<8}{:>8}{:>10}{:>9}{:>9.4}{:>9.3}{:>9.3}{:>8.3}  {}",
                    fitted.stratum.period,
                    fitted.stratum.reference_repeats,
                    own.tracts_with_reads(),
                    own.spanning_reads(),
                    slippage.level,
                    slippage.shorter_share,
                    slippage.fall_off,
                    fitted.concentration,
                    if fitted.borrowed.is_empty() {
                        "on its own".to_string()
                    } else {
                        format!(
                            "{} strata, {} tracts",
                            fitted.borrowed.len(),
                            fitted.tracts_fitted
                        )
                    }
                );
            }
            ssr_fit::StratumOutcome::Refused {
                stratum,
                tracts,
                reason,
            } => {
                println!(
                    "  {:<8}{:>8}{:>10}{:>9}   refused: {reason:?}",
                    stratum.period,
                    stratum.reference_repeats,
                    tracts,
                    held[stratum].spanning_reads(),
                );
            }
        }
    }
}

/// One sample: one walk, filling that sample's census evidence.
#[allow(
    clippy::too_many_arguments,
    reason = "a walk needs all of it, and a config struct \
                                              would only move the list"
)]
fn walk_one(
    fasta: &Path,
    info: &Arc<ReferenceInfo>,
    contigs: &Arc<ContigList>,
    index: &Arc<noodles_fasta::fai::Index>,
    alignment: &Path,
    typed: &[TypedRegion],
    _generic_domain: &SelectableRegions,
    kept: &pop_var_caller::ng::parameter_estimation::joint::loci::CensusLoci,
    terms: &SelectionTerms,
) -> SampleCensusEvidence {
    let read_groups =
        build_read_groups(&[alignment.to_path_buf()]).expect("the header declares read groups");
    let sample = match read_groups.read_groups_per_sample() {
        [only] => only.clone(),
        other => panic!(
            "{} holds {} samples; this walk is per sample",
            alignment.display(),
            other.len()
        ),
    };
    let reference = OpenReference::new(Arc::clone(info));
    let sample_reads = SampleReads::open(
        &sample,
        &read_groups,
        &reference,
        ReadFilterConfig::default(),
        true,
    )
    .expect("the alignment file opens against this reference");

    // Owned captures, because the generator outlives this function's borrows and keeps the
    // maker to mint a fresh view of the reference whenever it evicts one.
    let accessor = {
        let fasta = fasta.to_path_buf();
        let contigs = Arc::clone(contigs);
        let index = Arc::clone(index);
        move || {
            WindowedRefSeq::with_shared_index(
                fasta.clone(),
                Arc::clone(&contigs),
                Arc::clone(&index),
            )
        }
    };
    #[allow(
        clippy::arc_with_non_send_sync,
        reason = "file-backed and single-threaded, as in real_alignments.rs"
    )]
    let shared = Arc::new(accessor());
    let generic_generator = PileupGenerator::new(
        Arc::clone(&shared),
        accessor.clone(),
        LeftAlignPreparer::with_default_normalizer(accessor()),
        PileupGeneratorConfig::default(),
    )
    .expect("the generic generator builds against this reference");
    let bundle = Bp(TypedRegionConfig::default().criteria.bundle_threshold);
    let ssr_generator = SsrGenerator::with_default_aligner(
        Arc::clone(&shared),
        {
            let shared = Arc::clone(&shared);
            move || Arc::clone(&shared)
        },
        SsrGeneratorConfig {
            flank_bp: bundle,
            ..SsrGeneratorConfig::default()
        },
        bundle,
    )
    .expect("flank within the bundle threshold");
    let generators = GeneratorSet::new(
        GeneratorSlot::Generator(Box::new(ssr_generator)),
        GeneratorSlot::Generator(Box::new(generic_generator)),
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
    );

    let contig_of = |name: &str| {
        contigs
            .entries
            .iter()
            .position(|entry| entry.name == name)
            .map(|i| ContigId(i as u32))
    };
    let mut writer = CensusWriter::new(
        sample.sample.to_string(),
        kept,
        sample.read_groups.clone(),
        &contig_of,
        terms.clone(),
        DepthBinEdges::for_census(),
        ReadCap(pop_var_caller::ng::locus_generation::ssr::DEFAULT_SSR_MAX_READS_PER_LOCUS),
        // **Not the ladder's top any more, and the two are separate knobs.** The depth code
        // records what the position held, all the way to the ladder's ceiling near 1,500; this
        // is where a position's *reads* stop being counted one by one, and the allele counts
        // are thinned to it proportionally so the fractions they showed survive.
        DEPTH_CAP,
    );

    // **Every generic stretch handed to the walk is walked**, whether or not a read reached
    // it. Without saying so, a position no read reached is indistinguishable from a region the
    // run never opened, because the generic generator emits no locus where there is no read.
    for region in typed
        .iter()
        .filter(|region| region.kind == RegionKind::Generic)
    {
        writer.mark_walked(region.region);
    }

    let regions: Vec<Result<TypedRegion, RepeatCatalogError>> =
        typed.iter().cloned().map(Ok).collect();
    let mut stream =
        SampleLocusObservationsIterator::new(regions.into_iter(), sample_reads, generators);
    let mut generic_loci = 0_u64;
    let mut ssr_loci = 0_u64;
    for locus in &mut stream {
        let locus = locus.expect("the walk runs to completion on a well-formed alignment");
        match locus.kind {
            LocusKind::Generic => generic_loci += 1,
            LocusKind::Ssr(_) => ssr_loci += 1,
            _ => {}
        }
        writer.add_locus(&locus);
    }
    println!("  walked         {generic_loci} generic loci, {ssr_loci} STR loci");
    writer.finish()
}

/// What one sample's records weigh, each object on its own line.
///
/// **Separately and not as a total**, because `parameter_prepass_joint_records.md` §6 claims
/// the STR set is the larger half and a single number would hide that being wrong.
fn report_sizes(records: &SampleCensusEvidence) {
    let positions = records
        .generic
        .values()
        .next()
        .map_or(0, |g| g.depth().len());
    let depth_bytes: usize = records
        .generic
        .values()
        .map(|g| g.depth().as_bytes().len())
        .sum();
    let sparse_entries: usize = records
        .generic
        .values()
        .map(|g| g.non_reference().len())
        .sum();
    let sparse_bytes = sparse_entries
        * std::mem::size_of::<
            pop_var_caller::ng::parameter_estimation::joint::census::AlleleObservation,
        >();

    let ssr_loci = records.ssr.values().next().map_or(0, |s| s.len());
    // A tract costs its offset buckets and one **bit** saying the walk reached it. The two
    // counts that used to sit beside them — the reads that reached without crossing, and the
    // base-comparison denominator — are now one pair per stratum, which is the second term.
    let ssr_dense_bytes: usize = records
        .ssr
        .values()
        .map(|s| {
            s.len()
                * std::mem::size_of::<
                    pop_var_caller::ng::parameter_estimation::joint::census::OffsetCounts,
                >()
                + s.walked_bits().as_bytes().len()
                + s.covering_not_crossing_by_stratum().count()
                    * 2
                    * (std::mem::size_of::<
                        pop_var_caller::ng::parameter_estimation::joint::census::Stratum,
                    >() + std::mem::size_of::<u64>())
        })
        .sum();
    let guard_entries: usize = records.ssr.values().map(|s| s.guard().len()).sum();
    let difference_entries: usize = records.ssr.values().map(|s| s.differences().len()).sum();
    let difference_bytes = difference_entries
        * std::mem::size_of::<
            pop_var_caller::ng::parameter_estimation::joint::census::TractDifference,
        >();

    // The three states the depth array must keep apart.
    let mut never_walked = 0_u64;
    let mut zero_depth = 0_u64;
    let mut covered = 0_u64;
    if let Some(first) = records.generic.values().next() {
        for code in first.depth().iter() {
            match code {
                DepthCode::NeverWalked => never_walked += 1,
                DepthCode::Binned(bin) if bin.0 == 0 => zero_depth += 1,
                DepthCode::Binned(_) => covered += 1,
            }
        }
    }

    let mb = |bytes: usize| bytes as f64 / 1e6;
    println!(
        "  read groups    {}, {positions} kept generic positions, {ssr_loci} kept STR loci",
        records.generic.len()
    );
    println!(
        "  generic depth  {:.3} MB ({:.2} bits a position a read group)",
        mb(depth_bytes),
        8.0 * depth_bytes as f64 / (positions.max(1) * records.generic.len()) as f64
    );
    println!(
        "  generic sparse {:.3} MB, {sparse_entries} entries ({:.1} per 1,000 positions)",
        mb(sparse_bytes),
        1_000.0 * sparse_entries as f64 / positions.max(1) as f64
    );
    println!(
        "  STR dense      {:.3} MB over {ssr_loci} loci ({:.1} bytes a locus a read group)",
        mb(ssr_dense_bytes),
        ssr_dense_bytes as f64 / (ssr_loci.max(1) * records.ssr.len().max(1)) as f64
    );
    println!(
        "  STR guard      {guard_entries} entries; difference list {:.3} MB, \
         {difference_entries} entries",
        mb(difference_bytes)
    );
    println!(
        "  states         {never_walked} never walked, {zero_depth} walked at zero depth, \
         {covered} with reads"
    );
    println!(
        "  TOTAL          {:.3} MB held for this sample",
        mb(depth_bytes + sparse_bytes + ssr_dense_bytes + difference_bytes)
    );
}

// ---------------------------------------------------------------------
// How far apart are the samples?
// ---------------------------------------------------------------------

/// One position the cohort is decomposed on: the alternative allele it segregates, and each
/// sample's reads there.
struct Marker {
    /// Alternative reads and total depth, per sample.
    alt: Vec<u16>,
    depth: Vec<u16>,
    /// The cohort's allele frequency, from the samples that have data.
    frequency: f64,
}

/// A position enters the decomposition only if this many samples put a read on it.
const MIN_SAMPLES_WITH_DATA: usize = 8;
/// …and only if the cohort's own alternative frequency is inside this band. Outside it the
/// position carries almost no information about who differs from whom, and the standardising
/// divisor blows up.
const MIN_FREQUENCY: f64 = 0.05;

fn structure(cohort: &[SampleCensusEvidence], positions: usize) {
    println!("\n--- how far apart the samples are ---");
    let started = Instant::now();
    let markers = markers(cohort, positions);
    println!(
        "  markers        {} of {positions} kept positions segregate with at least \
         {MIN_SAMPLES_WITH_DATA} samples covered and frequency in {MIN_FREQUENCY}–{:.2} ({:.1} s)",
        markers.len(),
        1.0 - MIN_FREQUENCY,
        started.elapsed().as_secs_f64()
    );
    if markers.len() < 100 {
        println!("  too few markers to decompose");
        return;
    }

    let measured = decompose(&markers, cohort.len(), None);
    // The control: the same depths, the same missingness, the same allele frequencies, and no
    // structure whatever. Anything the pipeline returns here it manufactures.
    let mut rng = 0x2545_F491_4F6C_DD1D_u64;
    let null_markers = redraw(&markers, &mut rng);
    let null = decompose(&null_markers, cohort.len(), None);

    println!("\n  {:<28}{:>14}{:>14}", "", "measured", "no structure");
    for axis in 0..6.min(measured.eigenvalues.len()) {
        println!(
            "  principal axis {:<13}{:>13.2}%{:>13.2}%",
            axis + 1,
            100.0 * measured.share[axis],
            100.0 * null.share[axis]
        );
    }
    println!(
        "  {:<28}{:>14.4}{:>14.4}",
        "F_st, split on axis 1", measured.fst, null.fst
    );
    println!(
        "  {:<28}{:>14.4}{:>14.4}",
        "F_st, random split", measured.fst_random, null.fst_random
    );
    println!(
        "\n  structure above the null: {:+.4} in F_st on the leading split",
        measured.fst - null.fst
    );

    // Which samples sit where, so a split that is one accession against the rest is visible
    // as such rather than reported as divergence.
    let mut order: Vec<usize> = (0..cohort.len()).collect();
    order.sort_by(|&a, &b| {
        measured.axis1[a]
            .partial_cmp(&measured.axis1[b])
            .expect("no NaN on the axis")
    });
    println!("\n  samples along axis 1, most negative first");
    for chunk in order.chunks(4) {
        let line: Vec<String> = chunk
            .iter()
            .map(|&s| format!("{} {:+.3}", cohort[s].sample, measured.axis1[s]))
            .collect();
        println!("    {}", line.join("   "));
    }
    println!(
        "\n  {} of {} samples on the negative side of axis 1",
        measured.axis1.iter().filter(|v| **v < 0.0).count(),
        cohort.len()
    );
}

/// The positions worth decomposing on, with each sample's reads at them.
fn markers(cohort: &[SampleCensusEvidence], positions: usize) -> Vec<Marker> {
    let samples = cohort.len();
    // Per position, each sample's four allele counts. Built one position at a time would be
    // a binary search per sample per position; instead each sample's sparse list is swept
    // once into a per-position column.
    let mut alt_counts = vec![[0_u32; 5]; positions];
    let mut per_sample_alt: Vec<Vec<u16>> = vec![vec![0; positions]; samples];
    let mut per_sample_depth: Vec<Vec<u16>> = vec![vec![0; positions]; samples];

    // First sweep: which non-reference allele the cohort carries at each position.
    for records in cohort {
        for group in records.generic.values() {
            for observation in group.non_reference() {
                alt_counts[observation.index as usize][observation.allele.code() as usize] +=
                    u32::from(observation.reads);
            }
        }
    }
    let major: Vec<u8> = alt_counts
        .iter()
        .map(|counts| {
            counts
                .iter()
                .enumerate()
                .max_by_key(|(_, n)| **n)
                .map(|(allele, _)| allele as u8)
                .unwrap_or(0)
        })
        .collect();

    // Second sweep: each sample's reads on that allele, and its depth.
    for (s, records) in cohort.iter().enumerate() {
        for group in records.generic.values() {
            for (index, code) in group.depth().iter().enumerate() {
                if let DepthCode::Binned(bin) = code {
                    // The stored code is a bin, and the fit reads it as one. A bin's lower
                    // edge is exact below depth 9, which is where a three-read cohort lives.
                    let depth = u32::from(bin.0).min(u32::from(u16::MAX));
                    per_sample_depth[s][index] =
                        per_sample_depth[s][index].saturating_add(depth as u16);
                }
            }
            for observation in group.non_reference() {
                if observation.allele.code() == major[observation.index as usize] {
                    per_sample_alt[s][observation.index as usize] = per_sample_alt[s]
                        [observation.index as usize]
                        .saturating_add(u16::from(observation.reads));
                }
            }
        }
    }

    let mut markers = Vec::new();
    for index in 0..positions {
        let covered = (0..samples)
            .filter(|&s| per_sample_depth[s][index] > 0)
            .count();
        if covered < MIN_SAMPLES_WITH_DATA {
            continue;
        }
        let mut sum = 0.0;
        for s in 0..samples {
            if per_sample_depth[s][index] > 0 {
                sum += f64::from(per_sample_alt[s][index]) / f64::from(per_sample_depth[s][index]);
            }
        }
        let frequency = sum / covered as f64;
        if !(MIN_FREQUENCY..=1.0 - MIN_FREQUENCY).contains(&frequency) {
            continue;
        }
        markers.push(Marker {
            alt: (0..samples).map(|s| per_sample_alt[s][index]).collect(),
            depth: (0..samples).map(|s| per_sample_depth[s][index]).collect(),
            frequency,
        });
    }
    markers
}

struct Decomposition {
    eigenvalues: Vec<f64>,
    share: Vec<f64>,
    axis1: Vec<f64>,
    fst: f64,
    fst_random: f64,
}

fn decompose(markers: &[Marker], samples: usize, _unused: Option<()>) -> Decomposition {
    // The covariance between two samples, over the markers both have reads at.
    //
    // **The diagonal is not computed and not used.** A sample against itself carries its own
    // read-sampling noise, which is not shared with anyone and would dominate a shallow
    // cohort's leading axis; every off-diagonal cell is free of it because two samples' reads
    // are drawn independently.
    let mut covariance = vec![0.0_f64; samples * samples];
    let mut pairs = vec![0.0_f64; samples * samples];
    for marker in markers {
        let p = marker.frequency;
        let scale = (p * (1.0 - p)).sqrt();
        if scale <= 0.0 {
            continue;
        }
        let z: Vec<Option<f64>> = (0..samples)
            .map(|s| {
                (marker.depth[s] > 0)
                    .then(|| (f64::from(marker.alt[s]) / f64::from(marker.depth[s]) - p) / scale)
            })
            .collect();
        for a in 0..samples {
            let Some(za) = z[a] else { continue };
            for b in (a + 1)..samples {
                let Some(zb) = z[b] else { continue };
                covariance[a * samples + b] += za * zb;
                pairs[a * samples + b] += 1.0;
            }
        }
    }
    for a in 0..samples {
        for b in (a + 1)..samples {
            let cell = if pairs[a * samples + b] > 0.0 {
                covariance[a * samples + b] / pairs[a * samples + b]
            } else {
                0.0
            };
            covariance[a * samples + b] = cell;
            covariance[b * samples + a] = cell;
        }
    }
    // The diagonal stands in as the mean of the row's own off-diagonals, which is the value a
    // noise-free self-comparison would take under one shared axis of variation.
    for a in 0..samples {
        let mut sum = 0.0;
        for b in 0..samples {
            if a != b {
                sum += covariance[a * samples + b];
            }
        }
        covariance[a * samples + a] = sum / (samples - 1) as f64;
    }

    let (values, vectors) = jacobi_eigen(&covariance, samples);
    let total: f64 = values.iter().map(|v| v.abs()).sum();
    let share: Vec<f64> = values.iter().map(|v| v.abs() / total.max(1e-12)).collect();
    let axis1: Vec<f64> = (0..samples).map(|s| vectors[s * samples]).collect();

    // Split on the leading axis and measure how far the two halves' allele frequencies sit
    // apart; then do the same on a split that knows nothing, which is the floor a split
    // chosen by the data has to beat.
    let leading: Vec<bool> = axis1.iter().map(|v| *v < 0.0).collect();
    let mut random = vec![false; samples];
    for (s, slot) in random.iter_mut().enumerate() {
        *slot = s % 2 == 0;
    }
    Decomposition {
        eigenvalues: values,
        share,
        fst: hudson_fst(markers, &leading),
        fst_random: hudson_fst(markers, &random),
        axis1,
    }
}

/// Hudson's `F_st` between the two groups a boolean split names, as a ratio of averages.
///
/// A group's allele frequency at a position is the mean of its members' own read fractions,
/// so a deep sample and a shallow one count the same — which is what makes the number about
/// the plants rather than about the sequencing.
fn hudson_fst(markers: &[Marker], left: &[bool]) -> f64 {
    let mut numerator = 0.0_f64;
    let mut denominator = 0.0_f64;
    for marker in markers {
        let mut sums = [0.0_f64; 2];
        let mut counts = [0.0_f64; 2];
        for (s, depth) in marker.depth.iter().enumerate() {
            if *depth == 0 {
                continue;
            }
            let side = usize::from(!left[s]);
            sums[side] += f64::from(marker.alt[s]) / f64::from(*depth);
            counts[side] += 1.0;
        }
        if counts[0] < 2.0 || counts[1] < 2.0 {
            continue;
        }
        let p1 = sums[0] / counts[0];
        let p2 = sums[1] / counts[1];
        // Sampled chromosomes, two per diploid with data.
        let n1 = 2.0 * counts[0];
        let n2 = 2.0 * counts[1];
        numerator +=
            (p1 - p2).powi(2) - p1 * (1.0 - p1) / (n1 - 1.0) - p2 * (1.0 - p2) / (n2 - 1.0);
        denominator += p1 * (1.0 - p2) + p2 * (1.0 - p1);
    }
    if denominator <= 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

/// The same cohort with every trace of structure removed: each sample's alternative reads are
/// re-drawn from a binomial at the position's own cohort frequency and that sample's own
/// depth.
fn redraw(markers: &[Marker], state: &mut u64) -> Vec<Marker> {
    markers
        .iter()
        .map(|marker| {
            let alt: Vec<u16> = marker
                .depth
                .iter()
                .map(|depth| {
                    let mut drawn = 0_u16;
                    for _ in 0..*depth {
                        if uniform(state) < marker.frequency {
                            drawn += 1;
                        }
                    }
                    drawn
                })
                .collect();
            Marker {
                alt,
                depth: marker.depth.clone(),
                frequency: marker.frequency,
            }
        })
        .collect()
}

fn uniform(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 11) as f64 / (1_u64 << 53) as f64
}

/// Eigenvalues and eigenvectors of a small symmetric matrix, by cyclic Jacobi rotations.
///
/// Returned sorted by descending eigenvalue; `vectors[s * n + k]` is sample `s`'s coordinate
/// on axis `k`.
fn jacobi_eigen(matrix: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut a = matrix.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
            .map(|(i, j)| a[i * n + j].powi(2))
            .sum();
        if off < 1e-18 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p * n + q].abs() < 1e-15 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * a[p * n + q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let vkp = v[k * n + p];
                    let vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        a[j * n + j]
            .abs()
            .partial_cmp(&a[i * n + i].abs())
            .expect("no NaN in a covariance matrix")
    });
    let values: Vec<f64> = order.iter().map(|&i| a[i * n + i]).collect();
    let mut vectors = vec![0.0; n * n];
    for (k, &i) in order.iter().enumerate() {
        for s in 0..n {
            vectors[s * n + k] = v[s * n + i];
        }
    }
    (values, vectors)
}

/// Unused today, kept because the STR half of the report will want it.
#[allow(dead_code)]
fn per_read_group(records: &SampleCensusEvidence) -> BTreeMap<String, usize> {
    records
        .generic
        .iter()
        .map(|(group, g)| (format!("{group:?}"), g.non_reference().len()))
        .collect()
}
