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
    CensusWriter, CohortCensusEvidence, DepthCap, DepthCode, ReadCap, SampleCensusEvidence,
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

/// Every census in this program is one a walk just built and is holding, so no scoped call has a
/// file to fail on.
const RESIDENT: &str = "a resident census has no file to fail on";

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
        // **5,000 tracts a stratum, which is the figure the design settled on and not a number
        // chosen here.** `str_stratum_size_sweep_2026-08-13.md` drew strata at a known truth and
        // fitted them from 50 tracts up to 20,000: at three reads a site — tomato's depth — 5,000
        // is the smallest count where all five fitted numbers land within a few percent of the
        // truth both on average and between draws, and `parameter_prepass_joint_loci.md` §6
        // question 1 records the question closed on it.
        //
        // **The cap is what bounds the run, and it was previously set where it could never
        // fire.** Cost is linear in how many tracts a fit reads; this line stood at 1,000,000,
        // above the largest stratum tomato has, so nothing was ever capped and the 63-accession
        // cohort extrapolated to days. At 5,000 tomato keeps 86,688 of its 462,701 tracts and 8
        // of its 141 strata are capped at all — the three fattest strata hold 79% of the loci,
        // and they are the ones whose parameters are already best determined.
        //
        // **What it does not touch is the thin strata**, and that is most of them: 68 of the 141
        // hold fewer than a hundred tracts each, are far under any cap, and reach a fittable size
        // only by borrowing from neighbouring repeat counts.
        ssr_cap: 5_000,
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
        let mut records = walk_one(
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
        report_sizes(&mut records);
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

    // **From here on the samples are one cohort.** Its own door makes the same check the loop
    // above reports on, so it cannot refuse what that loop said was fine.
    let mut cohort = CohortCensusEvidence::new(cohort).expect("every sample recorded one way");

    if cohort.len() >= 8 {
        structure(&mut cohort, kept.generic().len());
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

    depth_ladder_occupancy(&mut cohort);

    // **The ordinary-position half runs once and its answer is handed on.** The repeat-tract
    // half needs each sample's homozygote excess from it and gives nothing back, so re-fitting
    // to obtain that would be 883 seconds of the tomato cohort's time spent twice.
    if let Some(fit) = fit_the_cohort(&mut cohort, coverage_odds, kept.generic(), &contigs) {
        fit_the_tracts(&mut cohort, &kept, &contigs, &fit);
    }

    println!(
        "\ntotal            {:.1} s",
        started.elapsed().as_secs_f64()
    );
}

/// Where each sample's positions sit on the depth ladder.
///
/// **The ladder is exact to the cap of 124 reads a position and a range only above it**, so the
/// middle column should read zero for every sample and is kept as the assertion that it does. It
/// was where a run's positions piled up until 2026-08-16, when a code above eight reads stopped
/// standing for several depths at once
/// (`doc/devel/ng/reports/census_depth_resolution_2026-08-16.md`).
///
/// **The last column is the one to watch, and it is about the cap rather than the ladder.**
/// The depth recorded is now the position's own, so the ladder no longer piles deep positions
/// onto one rung; what still collapses is the *allele counts*, which are thinned to the
/// per-position cap. A position above that cap has an exact depth and a thinned count beside
/// it, and the fraction it reports is the thinned one.
fn depth_ladder_occupancy(cohort: &mut CohortCensusEvidence) {
    let edges = DepthBinEdges::for_census();
    println!("\n--- where the positions sit on the depth ladder ---");
    println!(
        "  {:<24}{:>14}{:>18}{:>16}{:>12}",
        "sample", "an exact depth", "a range", "above the cap", "no reads"
    );
    let names: Vec<String> = cohort.sample_names().map(str::to_string).collect();
    let cap = cohort.terms().map_or(0, |terms| terms.depth_cap.get());
    let groups = cohort.read_groups().to_vec();
    // **The counting happens inside the call**, because that is how long the sections are lent
    // for: the census hands them over, the closure reads them, and it takes them back.
    let rows = cohort
        .with_generic(&groups, |lent| {
            lent.iter()
                .map(|sections| {
                    let (mut exact, mut ranged, mut capped, mut unwalked) =
                        (0_u64, 0_u64, 0_u64, 0_u64);
                    for (_, records) in sections {
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
                    (exact, ranged, capped, unwalked)
                })
                .collect::<Vec<_>>()
        })
        .expect(RESIDENT);
    for (name, (exact, ranged, capped, unwalked)) in names.iter().zip(rows) {
        let total = (exact + ranged + capped).max(1) as f64;
        println!(
            "  {:<24}{:>13.1}%{:>17.1}%{:>13.1}%{:>12}",
            name,
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
    cohort: &mut CohortCensusEvidence,
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
    // **The milestone-B assertion, on real reads**: the same cohort fitted from memory and from
    // files must give the same parameters. `CENSUS_FILES=<dir>` writes each sample's census
    // there, opens them again as files and refits; without it nothing is written and the run is
    // what it was.
    if let Ok(dir) = std::env::var("CENSUS_FILES") {
        refit_from_files(cohort, Path::new(&dir), &config, &fit);
    }

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
        // Each sample's non-reference reads and depth at each kept position, read while the
        // sections are lent — one call, since the closure sees the whole cohort at once.
        let groups = cohort.read_groups().to_vec();
        let (non_reference, depths) = cohort
            .with_generic(&groups, |lent| {
                let mut non_reference: Vec<Vec<u32>> = Vec::with_capacity(lent.len());
                let mut depths: Vec<Vec<u32>> = Vec::with_capacity(lent.len());
                for sections in lent {
                    let mut counts = vec![0_u32; kept.len()];
                    let mut depth = vec![0_u32; kept.len()];
                    for (_, group) in sections {
                        for observation in group.non_reference() {
                            counts[observation.index as usize] += u32::from(observation.reads);
                        }
                        for (index, code) in group.depth().iter().enumerate() {
                            if let DepthCode::Binned(bin) = code {
                                let range = edges.depth_range(bin);
                                depth[index] += (*range.start() + *range.end()) / 2;
                            }
                        }
                    }
                    non_reference.push(counts);
                    depths.push(depth);
                }
                (non_reference, depths)
            })
            .expect(RESIDENT);

        let mut out = String::new();
        out.push_str("contig\tposition");
        for name in cohort.sample_names() {
            out.push_str(&format!(
                "\t{0}_het\t{0}_homalt\t{0}_depth\t{0}_nonref",
                name
            ));
        }
        out.push('\n');
        let width = cohort.len() * 2;
        for (index, position) in kept.iter().enumerate() {
            let row = &fit.genotype_posterior[index * width..][..width];
            out.push_str(&contigs.entries[position.contig.get() as usize].name);
            out.push('\t');
            out.push_str(&position.position.get().to_string());
            for s in 0..cohort.len() {
                out.push_str(&format!(
                    "\t{:.4}\t{:.4}\t{}\t{}",
                    row[s * 2],
                    row[s * 2 + 1],
                    depths[s][index],
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
            100.0 * rates.positions_with_reads as f64 / kept.len() as f64
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
    let group = *cohort.read_groups().first().expect("a read group");
    let error: Vec<f64> = (0..cohort.len())
        .map(|_| fit.noise[&group].value.clean)
        .collect();
    let excess: Vec<f64> = cohort
        .sample_names()
        .map(|name| fit.hom_excess[name].value.get())
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
            )
            .expect("a resident census has no file to fail on");
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
    cohort: &mut CohortCensusEvidence,
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
    for group in cohort.read_groups().to_vec() {
        {
            let next = if per_read_group {
                slippage_group_of.len() as u32
            } else {
                0
            };
            slippage_group_of.entry(group).or_insert(next);
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
        .sample_names()
        .map(|name| {
            fit.hom_excess
                .get(name)
                .map_or(0.0, |estimate| estimate.value.get())
        })
        .collect();

    let at = Instant::now();
    let mut evidence = ssr_fit::gather_strata(cohort, &strata, &slippage_group_of)
        .expect("a resident census has no file to fail on");
    println!(
        "  {} strata over {} tracts, gathered in {:.1} s",
        evidence.len(),
        strata.len(),
        at.elapsed().as_secs_f64()
    );

    // **`SSR_TRACT_CAP=<n>` keeps only the first `n` tracts of each stratum**, which is the
    // question the run time turns on: the fit's cost is linear in how many tracts it reads, and
    // nobody has measured how many it needs before the numbers it returns stop moving. The cap is
    // the selection's own `ssr_cap` applied after the walk rather than before it, so one walk can
    // be fitted at several caps and the answers compared. Tracts are kept in genome order, which
    // is neither sorted by depth nor by length, so the kept ones are not a favourable sample.
    if let Ok(cap) = std::env::var("SSR_TRACT_CAP") {
        let cap: usize = cap.parse().expect("a tract count");
        let before: usize = evidence.iter().map(|s| s.tracts.len()).sum();
        for stratum in &mut evidence {
            stratum.tracts.truncate(cap);
        }
        let after: usize = evidence.iter().map(|s| s.tracts.len()).sum();
        println!("  tract cap {cap}: {before} tracts kept as {after}");
    }
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

    // **`SSR_CELL_TABLE=<path>` writes one row a stratum**, with the evidence beside the fit.
    // The printed table below is for reading; this one is for a curve to be fitted against, and
    // it therefore carries the three counts that set how sharply a stratum determines each of
    // its numbers — tracts, spanning reads, and reads sitting off the reference length — which
    // the printed table has never carried.
    if let Ok(path) = std::env::var("SSR_CELL_TABLE") {
        write_cell_table(&path, &outcomes, &held);
    }

    // **`SSR_CELL_TABLE_NO_CURVE=<path>` fits the same strata again with no curve drawn at all.**
    // Nothing about how a stratum is fitted depends on the curve, so every cell's own level must
    // come back identical between the two tables; a cell that moves is a defect in the plumbing
    // rather than a consequence of the design (`str_slippage_level_curve.md` §12).
    if let Ok(path) = std::env::var("SSR_CELL_TABLE_NO_CURVE") {
        let plain = ssr_fit::SsrFitConfig {
            curve: ssr_fit::SlippageCurveConfig {
                draw_curves: false,
                ..ssr_fit::SlippageCurveConfig::default()
            },
            ..config.clone()
        };
        let at = Instant::now();
        let outcomes = ssr_fit::fit_strata(&evidence, &homozygote_excess, &plain);
        println!(
            "  refitted with no curve in {:.1} s",
            at.elapsed().as_secs_f64()
        );
        write_cell_table(&path, &outcomes, &held);
    }

    // **`SSR_CELL_TABLE_BORROWED=<path>` fits the same strata a second time with borrowing on**,
    // and writes that table beside the first. Both fits read one walk's census, so the two
    // tables differ in the borrowing rule and in nothing else — which is what makes "how far
    // borrowing moved this cell" a difference rather than a comparison across runs.
    if let Ok(path) = std::env::var("SSR_CELL_TABLE_BORROWED") {
        let borrowed_config = ssr_fit::SsrFitConfig {
            allele_span,
            borrowing_floor: ssr_fit::DEFAULT_BORROWING_FLOOR,
            ..ssr_fit::SsrFitConfig::default()
        };
        let at = Instant::now();
        let borrowed = ssr_fit::fit_strata(&evidence, &homozygote_excess, &borrowed_config);
        println!(
            "  refitted in {:.1} s, borrowing below {} tracts",
            at.elapsed().as_secs_f64(),
            borrowed_config.borrowing_floor
        );
        write_cell_table(&path, &borrowed, &held);
    }
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
            ssr_fit::StratumOutcome::Derived(derived) => {
                let own = held[&derived.stratum];
                let slippage = derived
                    .slippage
                    .iter()
                    .flatten()
                    .next()
                    .expect("a derived stratum has at least one group with numbers");
                println!(
                    "  {:<8}{:>8}{:>10}{:>9}{:>9.4}{:>9.3}{:>9.3}{:>8}  nothing fitted here",
                    derived.stratum.period,
                    derived.stratum.reference_repeats,
                    own.tracts_with_reads(),
                    own.spanning_reads(),
                    slippage.level,
                    slippage.shorter_share,
                    slippage.fall_off,
                    "-",
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

/// Write every sample's census, read them back as files, and fit the cohort again.
///
/// **What it prints is the difference, not the second answer.** The two fits read the same
/// evidence through the same calls; if they disagree, the file is not the census.
fn refit_from_files(
    cohort: &CohortCensusEvidence,
    dir: &Path,
    config: &JointFitConfig,
    from_memory: &JointFit,
) {
    use pop_var_caller::ng::parameter_estimation::joint::census_file::{
        bytes_read, open_census, reset_bytes_read, write_census,
    };

    println!("\n--- the same cohort, fitted from files ---");
    std::fs::create_dir_all(dir).expect("the census directory is writable");
    let at = Instant::now();
    let mut written = 0_u64;
    let mut backed = Vec::new();
    for sample in cohort.samples() {
        let path = dir.join(format!("{}.census", sample.sample));
        let mut file = std::fs::File::create(&path).expect("a census file is writable");
        write_census(sample, None, &mut file).expect("a census writes");
        drop(file);
        written += std::fs::metadata(&path).expect("the census exists").len();
        backed.push(open_census(&path).expect("this build's own census").0);
    }
    println!(
        "  wrote          {} files, {:.3} MB in total ({:.3} MB a sample), in {:.1} s",
        backed.len(),
        written as f64 / 1e6,
        written as f64 / 1e6 / backed.len().max(1) as f64,
        at.elapsed().as_secs_f64()
    );

    let mut from_files = CohortCensusEvidence::new(backed).expect("every census records one way");
    let at = Instant::now();
    reset_bytes_read();
    let refit = match fit_jointly(&mut from_files, config) {
        Ok(refit) => refit,
        Err(error) => {
            println!("  refused: {error}");
            return;
        }
    };
    let read = bytes_read();
    println!(
        "  read back      {:.3} MB of section in {:.1} s, {:.2} times what the files hold",
        read as f64 / 1e6,
        at.elapsed().as_secs_f64(),
        read as f64 / written.max(1) as f64
    );

    // The largest disagreement anywhere, in the numbers a reader of this program acts on.
    let mut worst: (f64, String) = (0.0, "nothing".to_string());
    let mut note = |gap: f64, what: String| {
        if gap.abs() > worst.0 {
            worst = (gap.abs(), what);
        }
    };
    note(
        refit.log_likelihood - from_memory.log_likelihood,
        "the log-likelihood".to_string(),
    );
    for (name, rates) in &from_memory.rates {
        note(
            refit.rates[name].value.heterozygous - rates.value.heterozygous,
            format!("{name}'s heterozygosity"),
        );
        note(
            refit.hom_excess[name].value.get() - from_memory.hom_excess[name].value.get(),
            format!("{name}'s homozygote excess"),
        );
    }
    for (group, noise) in &from_memory.noise {
        note(
            refit.noise[group].value.clean - noise.value.clean,
            format!("{group:?}'s error rate"),
        );
    }
    if worst.0 == 0.0 {
        println!("  agreement      every fitted number identical, to the last bit");
    } else {
        println!("  DISAGREEMENT   {} differs by {:e}", worst.1, worst.0);
    }
}

/// What one sample's records weigh, each object on its own line.
///
/// **Separately and not as a total**, because `parameter_prepass_joint_records.md` §6 claims
/// the STR set is the larger half and a single number would hide that being wrong.
fn report_sizes(records: &mut SampleCensusEvidence) {
    use pop_var_caller::ng::parameter_estimation::joint::census::{
        AlleleObservation, OffsetCounts, TractDifference,
    };

    let groups = records.read_groups();
    let strata = records.strata();

    // ---- the ordinary positions, read while their sections are lent ----------------------
    let (positions, depth_bytes, sparse_entries, never_walked, zero_depth, covered) = records
        .with_generic(&groups, |sections| {
            let positions = sections.first().map_or(0, |g| g.depth().len());
            let depth_bytes: usize = sections.iter().map(|g| g.depth().as_bytes().len()).sum();
            let sparse_entries: usize = sections.iter().map(|g| g.non_reference().len()).sum();
            // The three states the depth array must keep apart, counted on one read group.
            let (mut never_walked, mut zero_depth, mut covered) = (0_u64, 0_u64, 0_u64);
            if let Some(first) = sections.first() {
                for code in first.depth().iter() {
                    match code {
                        DepthCode::NeverWalked => never_walked += 1,
                        DepthCode::Binned(bin) if bin.0 == 0 => zero_depth += 1,
                        DepthCode::Binned(_) => covered += 1,
                    }
                }
            }
            (
                positions,
                depth_bytes,
                sparse_entries,
                never_walked,
                zero_depth,
                covered,
            )
        })
        .expect(RESIDENT);
    let sparse_bytes = sparse_entries * std::mem::size_of::<AlleleObservation>();

    // ---- the tracts, one read group's band of strata at a time ---------------------------
    // **A tract costs its offset buckets and one bit saying the walk reached it.** The two
    // counts that used to sit beside them — the reads that reached without crossing, and the
    // base-comparison denominator — are one number each *for the section*, since a section is
    // one stratum, so they no longer scale with the tracts at all.
    let (mut ssr_loci, mut ssr_dense_bytes, mut guard_entries, mut difference_entries) =
        (0_usize, 0_usize, 0_usize, 0_usize);
    let (mut highest_read, mut at_the_ceiling) = (0_u16, 0_usize);
    for (which, group) in groups.iter().enumerate() {
        let (loci, dense, guard, differences, highest, saturated) = records
            .with_strata(*group, &strata, |sections| {
                let mut loci = 0_usize;
                let mut dense = 0_usize;
                let mut guard = 0_usize;
                let mut differences = 0_usize;
                let (mut highest, mut saturated) = (0_u16, 0_usize);
                for section in sections {
                    loci += section.len();
                    dense += section.len() * std::mem::size_of::<OffsetCounts>()
                        + section.walked_bits().as_bytes().len()
                        + 2 * std::mem::size_of::<u64>();
                    guard += section.guard().len();
                    differences += section.differences().len();
                    for difference in section.differences() {
                        highest = highest.max(difference.read);
                        saturated += usize::from(difference.read == u16::MAX);
                    }
                }
                (loci, dense, guard, differences, highest, saturated)
            })
            .expect(RESIDENT);
        // The tract count is one read group's, as the position count above is.
        if which == 0 {
            ssr_loci = loci;
        }
        ssr_dense_bytes += dense;
        guard_entries += guard;
        difference_entries += differences;
        highest_read = highest_read.max(highest);
        at_the_ceiling += saturated;
    }
    let difference_bytes = difference_entries * std::mem::size_of::<TractDifference>();

    let mb = |bytes: usize| bytes as f64 / 1e6;
    println!(
        "  read groups    {}, {positions} kept generic positions, {ssr_loci} kept STR loci",
        groups.len()
    );
    println!(
        "  generic depth  {:.3} MB ({:.2} bits a position a read group)",
        mb(depth_bytes),
        8.0 * depth_bytes as f64 / (positions.max(1) * groups.len()) as f64
    );
    println!(
        "  generic sparse {:.3} MB, {sparse_entries} entries ({:.1} per 1,000 positions)",
        mb(sparse_bytes),
        1_000.0 * sparse_entries as f64 / positions.max(1) as f64
    );
    println!(
        "  STR dense      {:.3} MB over {ssr_loci} loci ({:.1} bytes a locus a read group)",
        mb(ssr_dense_bytes),
        ssr_dense_bytes as f64 / (ssr_loci.max(1) * groups.len().max(1)) as f64
    );
    // **How close the read a difference sits on comes to the field that holds it.** The read is
    // numbered within the locus and the field reaches 65,535, past the 1,000-read cap a locus
    // is entered under — so this should never saturate. It is reported because it did, at 255,
    // while the field was one byte: on the trio one mismatch in five came back at the ceiling,
    // several reads reported as one. Printed rather than assumed, in both directions.
    println!(
        "  STR guard      {guard_entries} entries; difference list {:.3} MB, \
         {difference_entries} entries; highest read {highest_read}, {at_the_ceiling} at the \
         field's ceiling",
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

fn structure(cohort: &mut CohortCensusEvidence, positions: usize) {
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
    let names: Vec<String> = cohort.sample_names().map(str::to_string).collect();
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
            .map(|&s| format!("{} {:+.3}", names[s], measured.axis1[s]))
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
fn markers(cohort: &mut CohortCensusEvidence, positions: usize) -> Vec<Marker> {
    let samples = cohort.len();
    // Per position, each sample's four allele counts. Built one position at a time would be
    // a binary search per sample per position; instead each sample's sparse list is swept
    // once into a per-position column.
    let mut alt_counts = vec![[0_u32; 5]; positions];
    let mut per_sample_alt: Vec<Vec<u16>> = vec![vec![0; positions]; samples];
    let mut per_sample_depth: Vec<Vec<u16>> = vec![vec![0; positions]; samples];

    // First sweep: which non-reference allele the cohort carries at each position.
    let groups = cohort.read_groups().to_vec();
    cohort
        .with_generic(&groups, |lent| {
            for sections in lent {
                for (_, group) in sections {
                    for observation in group.non_reference() {
                        alt_counts[observation.index as usize]
                            [observation.allele.code() as usize] += u32::from(observation.reads);
                    }
                }
            }
        })
        .expect(RESIDENT);
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
    cohort
        .with_generic(&groups, |lent| {
            for (s, sections) in lent.iter().enumerate() {
                for (_, group) in sections {
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
        })
        .expect(RESIDENT);

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
fn per_read_group(records: &mut SampleCensusEvidence) -> BTreeMap<String, usize> {
    let groups = records.read_groups();
    records
        .with_generic(&groups, |sections| {
            groups
                .iter()
                .zip(sections)
                .map(|(group, g)| (format!("{group:?}"), g.non_reference().len()))
                .collect()
        })
        .expect(RESIDENT)
}

/// One row a stratum: the evidence the walk gathered, then the fit if there was one.
///
/// **The printed table is for reading and this one is for fitting a curve against**, so it
/// carries the counts that set how sharply a stratum determines each of its numbers — tracts,
/// reads crossing, and reads sitting off the reference length — which the printed table never
/// carried.
fn write_cell_table(
    path: &str,
    outcomes: &[ssr_fit::StratumOutcome],
    held: &BTreeMap<ssr_fit::Stratum, &ssr_fit::StratumEvidence>,
) {
    let mut out = String::from(
        "period,repeats,tracts,spanning_reads,off_reference_reads,bases_compared,\
         mismatching_bases,substitution_rate,fitted,level,shorter_share,fall_off,concentration,\
         log_likelihood_a_tract,converged,tracts_fitted,borrowed_from,level_source,curve_weight,\
         curve_reach,curve_shape,curve_held_out_error,curve_cells\n",
    );
    for outcome in outcomes {
        let stratum = outcome.stratum();
        let own = held[&stratum];
        let substitution = own
            .substitution_rate()
            .map_or(String::new(), |rate| format!("{rate:.8}"));
        let fit = match outcome {
            ssr_fit::StratumOutcome::Fitted(fitted) => {
                let slippage = fitted
                    .slippage
                    .iter()
                    .flatten()
                    .next()
                    .expect("a fitted stratum has at least one group with reads");
                format!(
                    "1,{:.8},{:.6},{:.6},{:.6},{:.6},{},{},{}",
                    slippage.level,
                    slippage.shorter_share,
                    slippage.fall_off,
                    fitted.concentration,
                    fitted.log_likelihood_a_tract,
                    u8::from(fitted.converged),
                    fitted.tracts_fitted,
                    fitted
                        .borrowed
                        .iter()
                        .map(u64::to_string)
                        .collect::<Vec<_>>()
                        .join(" "),
                )
            }
            // **A derived stratum is fitted from the caller's point of view and from nothing
            // else's**: it carries the three numbers a read likelihood asks for and no spectrum,
            // concentration or log-likelihood, because nothing here was estimated.
            ssr_fit::StratumOutcome::Derived(derived) => {
                let slippage = derived
                    .slippage
                    .iter()
                    .flatten()
                    .next()
                    .expect("a derived stratum has at least one group with numbers");
                format!(
                    "derived,{:.8},{:.6},{:.6},,,,0,",
                    slippage.level, slippage.shorter_share, slippage.fall_off,
                )
            }
            ssr_fit::StratumOutcome::Refused { .. } => "0,,,,,,,,".to_string(),
        };
        // **Where the level came from, beside the level.** After the curve, a value fitted from
        // eight thousand slipped reads and one drawn across a gap are the same number, and these
        // six columns are what tells them apart (`str_slippage_level_curve.md` §8).
        let provenance = match outcome {
            ssr_fit::StratumOutcome::Refused { .. } => ",,,,,".to_string(),
            _ => outcome
                .level_provenance()
                .iter()
                .flatten()
                .next()
                .map(|provenance| {
                    let (source, weight) = match provenance.source {
                        ssr_fit::LevelSource::Cell => ("cell".to_string(), 0.0),
                        ssr_fit::LevelSource::Curve => ("curve".to_string(), 1.0),
                        ssr_fit::LevelSource::Blend { curve_weight } => {
                            ("blend".to_string(), curve_weight)
                        }
                    };
                    let reach = match provenance.reach {
                        Some(ssr_fit::CurveReach::Inside) => "inside",
                        Some(ssr_fit::CurveReach::BelowFitted) => "below",
                        Some(ssr_fit::CurveReach::AboveFitted) => "above",
                        None => "",
                    };
                    match provenance.curve {
                        Some(curve) => format!(
                            "{source},{weight:.6},{reach},{:.2},{:.6},{}",
                            curve.rise_shape.get(),
                            curve.held_out_error,
                            curve.cells
                        ),
                        None => format!("{source},{weight:.6},{reach},,,"),
                    }
                })
                .unwrap_or_else(|| ",,,,,".to_string()),
        };
        out.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{}\n",
            stratum.period,
            stratum.reference_repeats,
            own.tracts_with_reads(),
            own.spanning_reads(),
            own.reads_off_reference_length(),
            own.bases_compared,
            own.mismatching_bases,
            substitution,
            fit,
            provenance,
        ));
    }
    std::fs::write(path, out).expect("the cell table's directory exists");
    println!("  cell table written to {path}");
}
