//! What the joint route's records actually cost, and what a cohort of them says.
//!
//! Everything measured for this route so far was measured against truths the measuring
//! program made up itself. This one drives the real writer —
//! `parameter_estimation::joint::records::RecordWriter` — through real alignments, and it
//! answers three questions the specifications assert without evidence:
//!
//! 1. **What do the records weigh?** `parameter_prepass_joint_records.md` §6 prices them by
//!    arithmetic and §7.8 says to measure instead. The generic depth array, the sparse list of
//!    non-reference observations, the STR set and the coverage-by-window summary are reported
//!    **separately**, because a single total would hide any one of them being wrong.
//! 2. **Does the coverage-by-window summary come out where §4 says it does?** It is built in
//!    the same walk, by `parameter_estimation::joint::coverage`, and its size and its
//!    short-window list are printed rather than assumed.
//! 3. **How much population structure does the cohort have?** Handed more than one alignment,
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
use pop_var_caller::ng::parameter_estimation::joint::contamination::{
    ContaminationConfig, ContaminationEstimate, OwnCoordinates, fit_contamination,
};
use pop_var_caller::ng::parameter_estimation::joint::coverage::{
    CoverageAccumulator, CoverageGrid, GC_BINS,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{JointFit, JointFitConfig, fit_jointly};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, ReferenceDigest, RegionSetDigest, SelectableRegions, SelectionIdentity,
    UnambiguousRuns, select_kept_loci,
};
use pop_var_caller::ng::parameter_estimation::joint::records::{
    DepthCap, DepthCode, ReadCap, RecordWriter, SampleRecords,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::read::left_align::LeftAlignPreparer;
use pop_var_caller::ng::ref_seq::{RefSeq, WindowedRefSeq};
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

/// The stored coverage window, and the one the identity makes every sample agree on.
const WINDOW_BP: Bp = Bp(500);

/// The seed the locus selection is drawn with. One number, shared by every sample, and part
/// of what the fit refuses to pool across.
const SEED: u64 = 42;

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

    let identity = SelectionIdentity {
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

    let kept = select_kept_loci(&identity, &catalog, &analysed, &unambiguous)
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

    // The walk, the coverage denominator and the selection all run over the same stretch of
    // genome: the analysed regions with the ambiguous runs cut out.
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
    let grid = CoverageGrid::over(&generic_domain, WINDOW_BP);
    println!(
        "coverage grid    {} windows of {} bp over {} generic bases",
        grid.len(),
        WINDOW_BP.get(),
        generic_domain.total_length().get()
    );

    // ---- one walk per sample -------------------------------------------------------------
    let mut cohort: Vec<SampleRecords> = Vec::new();
    let mut window_gc: Vec<f32> = Vec::new();
    for alignment in &alignments {
        let at = Instant::now();
        let (records, gc) = walk_one(
            &fasta,
            &info,
            &contigs,
            &index,
            alignment,
            &typed,
            &generic_domain,
            &grid,
            &kept,
            &identity,
        );
        println!(
            "\n=== {} — {:.1} s",
            records.sample,
            at.elapsed().as_secs_f64()
        );
        report_sizes(&records);
        cohort.push(records);
        window_gc = gc;
    }

    // ---- do they pool at all? --------------------------------------------------------------
    println!("\n--- the identity check, which is what lets these be pooled ---");
    let first = &cohort[0];
    let mut refused = 0;
    for other in &cohort[1..] {
        if let Some(field) = first.identity.first_disagreement(&other.identity) {
            println!(
                "  {} disagrees with {} on {field}",
                other.sample, first.sample
            );
            refused += 1;
        }
    }
    println!(
        "  {} of {} samples agree on all thirteen values",
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

    // **The second discriminator for the duplicated class, and the only one that works at
    // three samples.** A plant carrying two copies of a stretch the reference holds once
    // collects twice the reads around it. The cohort pattern — a position where the carriers
    // are all at a half and nobody is homozygous for the non-reference allele — needs about
    // twenty-five samples before that absence means anything; a window's depth says the same
    // thing whoever else is in the run. This turns each sample's coverage summary into
    // `ln P(this depth | two copies) − ln P(this depth | one)` at every kept position.
    // **Off unless asked for**, because the GC correction below is missing: on tomato, where
    // coverage runs from 16.2 reads a base at 20% GC to 29.0 at 36%, an uncorrected reading is
    // mostly GC content. On the human trio's benchmark regions at 300 reads a position it is
    // the only discriminator there is, three samples being far below the twenty-five the
    // cohort pattern needs.
    let coverage_odds = if std::env::var("COVERAGE_DISCRIMINATOR").is_ok() {
        let odds = coverage_odds(&cohort, &grid, kept.generic());
        let loud = odds
            .first()
            .map_or(0, |sample| sample.iter().filter(|o| **o > 0.0).count());
        println!(
            "\ncoverage odds    built for {} samples; {loud} of {} positions read more like two \
             copies than one in the first sample",
            odds.len(),
            kept.generic().len()
        );
        odds
    } else {
        Vec::new()
    };

    depth_code_census(&cohort);

    // **What the duplicated class claims on real reads, and whether the reads agree.** The
    // class is on by default and on this panel it takes the median accession's heterozygosity
    // from 0.867 to 0.064 per kilobase; the audit fits the cohort both ways and asks the two
    // questions the fit itself cannot answer — do the positions it claims carry twice their
    // sample's normal depth, and how many samples carry the alternative allele there.
    match std::env::var("DUPLICATED_CLASS_AUDIT") {
        Ok(path) => duplicated_class_audit(
            &cohort,
            &grid,
            &window_gc,
            kept.generic(),
            &contigs,
            if path == "1" {
                None
            } else {
                Some(path.as_str())
            },
        ),
        Err(_) => fit_the_cohort(&cohort, coverage_odds, kept.generic(), &contigs),
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
/// three reads a position is almost entirely exact; one at three hundred is entirely in the
/// ladder's top bin, which is the cap and so exact again because a deeper position was thinned
/// down to it. The band in between is where a stored code stands for several depths at once.
fn depth_code_census(cohort: &[SampleRecords]) {
    let edges = DepthBinEdges::new();
    println!("\n--- where the positions sit on the depth ladder ---");
    println!(
        "  {:<24}{:>14}{:>18}{:>14}{:>12}",
        "sample", "exact (0–8)", "a range (9–97)", "the cap (124)", "no reads"
    );
    for sample in cohort {
        let (mut exact, mut ranged, mut capped, mut unwalked) = (0_u64, 0_u64, 0_u64, 0_u64);
        for records in sample.generic.values() {
            for code in records.depth().iter() {
                match code {
                    DepthCode::NeverWalked => unwalked += 1,
                    DepthCode::Binned(bin) => {
                        let range = edges.recorded_depths(bin);
                        if range.start() == range.end() {
                            if *range.start() == edges.max_depth() {
                                capped += 1;
                            } else {
                                exact += 1;
                            }
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
    cohort: &[SampleRecords],
    coverage_odds: Vec<Arc<[f32]>>,
    kept: &[pop_var_caller::ng::types::GenomePosition],
    contigs: &ContigList,
) {
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
            return;
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
        let edges = DepthBinEdges::new();
        let non_reference: Vec<Vec<u32>> = cohort
            .iter()
            .map(|sample| {
                let mut counts = vec![0_u32; kept.len()];
                for group in sample.generic.values() {
                    for observation in group.non_reference() {
                        counts[observation.index as usize] += observation.reads;
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
                        let range = edges.recorded_depths(bin);
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
}

// =====================================================================================
// What the duplicated class claims on real reads
// =====================================================================================

/// The GC-corrected relative depth of every kept position, for one sample.
///
/// **1.0 is this sample's own normal depth and 2.0 is twice it**, which is what a stretch of
/// genome the sample carries twice and the reference holds once collects. Three steps, all of
/// them the ones `duplicated_locus_probe_2026-08-12.md` §2 measured on real alignments:
///
/// 1. **Windows are summed until they hold about 12,000 aligned bases**, because below that a
///    window's mean depth is scatter — one 500 bp window at 25 reads a position, ten at 2.5.
///    Only windows that abut *in the genome* are summed: a run over a BED holds windows that
///    are neighbours in the summary and megabases apart on the chromosome.
/// 2. **The sum is divided by what this sample's coverage does at that stretch's GC content.**
///    On tomato the depth-against-GC curve spans a factor of 1.79 — 16.2 reads a position at
///    20% GC against 29.0 at 36% — which is larger than the doubling being looked for, so an
///    uncorrected reading is mostly GC content.
/// 3. **The result is rescaled so the median window reads 1.0.**
///
/// `f32::NAN` at a position no window holds, or where the reference has no called base to take
/// a GC fraction from.
fn relative_depth(
    records: &SampleRecords,
    grid: &CoverageGrid,
    window_gc: &[f32],
    kept: &[pop_var_caller::ng::types::GenomePosition],
) -> Vec<f32> {
    let coverage = records
        .coverage
        .as_ref()
        .expect("every walk here builds a coverage summary");
    let span = coverage.windows_to_sum();
    let curve = coverage.gc_curve();
    let mut raw = vec![f32::NAN; coverage.len()];
    for slot in 0..coverage.len() {
        let (mut first, mut last) = (slot, slot + 1);
        while last - first < span {
            let left = first > 0 && grid.adjacent(first - 1, first);
            let right = last < coverage.len() && grid.adjacent(last - 1, last);
            if !left && !right {
                break;
            }
            if left && (!right || slot - first <= last - slot) {
                first -= 1;
            } else {
                last += 1;
            }
        }
        let (mut gc_sum, mut counted) = (0.0_f64, 0_u32);
        for window in first..last {
            if window_gc[window].is_finite() {
                gc_sum += f64::from(window_gc[window]);
                counted += 1;
            }
        }
        if counted == 0 {
            continue;
        }
        let bin = ((gc_sum / f64::from(counted) * (GC_BINS - 1) as f64).round() as usize)
            .min(GC_BINS - 1);
        let expected = f64::from(curve[bin]).max(1e-6);
        raw[slot] = (f64::from(coverage.mean_depth_over(first..last)) / expected) as f32;
    }
    let mut finite: Vec<f32> = raw.iter().copied().filter(|v| v.is_finite()).collect();
    let scale = f64::from(median_of(&mut finite)).max(1e-6);
    kept.iter()
        .map(|position| {
            grid.slot_of(position.contig, position.position)
                .map_or(f32::NAN, |slot| (f64::from(raw[slot]) / scale) as f32)
        })
        .collect()
}

/// Sorts in place and hands back the middle value, or `NaN` where there is none.
fn median_of(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    values.sort_by(f32::total_cmp);
    values[values.len() / 2]
}

/// One accession's reads at every kept position, as the audit needs them.
struct SampleReadsAtPositions {
    /// Reads disagreeing with the reference.
    non_reference: Vec<u32>,
    /// Reads in total, taken as the middle of the range the stored code stands for.
    depth: Vec<u32>,
}

fn reads_at_positions(records: &SampleRecords, positions: usize) -> SampleReadsAtPositions {
    let edges = DepthBinEdges::new();
    let mut non_reference = vec![0_u32; positions];
    let mut depth = vec![0_u32; positions];
    for group in records.generic.values() {
        for observation in group.non_reference() {
            non_reference[observation.index as usize] += observation.reads;
        }
        for index in 0..positions {
            if let DepthCode::Binned(bin) = group.depth().get(index) {
                let range = edges.recorded_depths(bin);
                depth[index] += (*range.start() + *range.end()) / 2;
            }
        }
    }
    SampleReadsAtPositions {
        non_reference,
        depth,
    }
}

/// The fit's summary numbers, for one arm of the audit.
struct ArmSummary {
    median_heterozygosity: f64,
    median_hom_excess: f64,
    expected_heterozygosity: f64,
    density_a: f64,
    density_b: f64,
    p_invariant: f64,
    class_share: f64,
    carrier_a: f64,
    carrier_b: f64,
    passes: u32,
    converged: bool,
}

fn fit_one_arm(
    cohort: &[SampleRecords],
    name: &str,
    with_class: bool,
    coverage_odds: Vec<Arc<[f32]>>,
) -> (JointFit, ArmSummary) {
    let at = Instant::now();
    let config = JointFitConfig {
        quadrature_nodes: 12,
        duplicated_positions: with_class,
        coverage_odds,
        genotype_posteriors: true,
        max_passes: std::env::var("MAX_PASSES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(200),
        ..JointFitConfig::default()
    };
    let fit = fit_jointly(cohort, &config).expect("the cohort fits");
    let mut heterozygosity: Vec<f32> = fit
        .rates
        .values()
        .map(|estimate| estimate.value.heterozygous as f32)
        .collect();
    let mut excess: Vec<f32> = fit
        .hom_excess
        .values()
        .map(|estimate| estimate.value.get() as f32)
        .collect();
    let summary = ArmSummary {
        median_heterozygosity: f64::from(median_of(&mut heterozygosity)),
        median_hom_excess: f64::from(median_of(&mut excess)),
        expected_heterozygosity: fit.expected_heterozygosity,
        density_a: fit.density.value.a,
        density_b: fit.density.value.b,
        p_invariant: fit.density.value.p_invariant,
        class_share: fit.duplicated.as_ref().map_or(0.0, |d| d.value.share),
        carrier_a: fit.duplicated.as_ref().map_or(0.0, |d| d.value.carrier_a),
        carrier_b: fit.duplicated.as_ref().map_or(0.0, |d| d.value.carrier_b),
        passes: fit.passes,
        converged: fit.converged,
    };
    println!(
        "  {:<30}{:>8}{:>12}{:>14.3}{:>14.3}{:>14.3}{:>10.0} s",
        name,
        summary.passes,
        if summary.converged {
            "converged"
        } else {
            "ran out"
        },
        1_000.0 * summary.median_heterozygosity,
        summary.median_hom_excess,
        1_000.0 * summary.expected_heterozygosity,
        at.elapsed().as_secs_f64()
    );
    (fit, summary)
}

/// **Does the class claim duplications, or rare real variants?** The fit cannot answer that
/// about itself — inside its own model a claimed position is a duplication by construction —
/// so this asks two things the fit never looks at.
///
/// 1. **Read depth.** A stretch of genome a plant carries twice collects two copies' reads, so
///    the positions the class claims should sit at about twice that accession's normal depth
///    for their GC content. If they sit at ordinary depth they are not duplications. This is
///    independent of every piece of evidence the class itself uses.
/// 2. **How many accessions carry the alternative allele.** A duplication and a rare variant
///    both have few carriers, so a claim rate that is flat in the carrier count is what a
///    duplication population looks like, and one that climbs as the allele gets rarer is what
///    claiming rare variants looks like.
///
/// **Neither is a truth set**, and there is none for duplications in this cohort. Every number
/// below compares two discriminators with each other.
fn duplicated_class_audit(
    cohort: &[SampleRecords],
    grid: &CoverageGrid,
    window_gc: &[f32],
    kept: &[pop_var_caller::ng::types::GenomePosition],
    contigs: &ContigList,
    dump: Option<&str>,
) {
    let positions = kept.len();
    let samples = cohort.len();
    println!("\n--- the duplicated class on real reads: the same cohort fitted three ways ---");
    println!(
        "  {:<30}{:>8}{:>12}{:>14}{:>14}{:>14}{:>12}",
        "", "passes", "", "median het/kb", "less het by", "population/kb", "time"
    );

    // The class-off arm's per-position heterozygosity is all that is kept from it: which
    // positions the class-on arm takes its heterozygosity from is the question, and one number
    // a position answers it. The 1.5 GB of per-sample posteriors goes when the fit does.
    let (heterozygous_without_class, without) = {
        let (fit, summary) = fit_one_arm(cohort, "the class off", false, Vec::new());
        let per_sample: Vec<f32> = (0..positions * samples)
            .map(|slot| fit.genotype_posterior[slot * 3])
            .collect();
        (per_sample, summary)
    };
    let heterozygosity_without_class: Vec<f32> = (0..positions)
        .map(|index| {
            heterozygous_without_class[index * samples..][..samples]
                .iter()
                .sum()
        })
        .collect();
    let (fit, with) = fit_one_arm(cohort, "the class on", true, Vec::new());

    // **The third arm: the class handed the depth reading as well.** Everything the class knows
    // by default is the cohort's genotype composition; this arm adds, per accession per
    // position, how much more likely that accession's read depth around there is if it has two
    // copies rather than one — GC-corrected, which the stored summary alone cannot do.
    let with_coverage = if std::env::var("SKIP_COVERAGE_ARM").is_ok() {
        None
    } else {
        let odds: Vec<Arc<[f32]>> = cohort
            .iter()
            .map(|records| {
                let coverage = records.coverage.as_ref().expect("a coverage summary");
                let aligned = f64::from(coverage.median_depth())
                    * coverage.window_bp().get() as f64
                    * coverage.windows_to_sum() as f64;
                let sd = |copies: f64| {
                    (copies * copies * 0.194 * 0.194 + copies * 313.0 / aligned.max(1.0)).sqrt()
                };
                let relative = relative_depth(records, grid, window_gc, kept);
                let odds: Vec<f32> = relative
                    .iter()
                    .map(|value| {
                        if !value.is_finite() {
                            return 0.0;
                        }
                        let reading = f64::from(*value).clamp(0.02, 8.0);
                        let ln_normal = |copies: f64| {
                            let s = sd(copies);
                            -0.5 * ((reading - copies) / s).powi(2) - s.ln()
                        };
                        (ln_normal(2.0) - ln_normal(1.0)) as f32
                    })
                    .collect();
                Arc::from(odds)
            })
            .collect();
        let (_, summary) = fit_one_arm(cohort, "the class on, with depth", true, odds);
        Some(summary)
    };

    println!("\n  the population's allele-frequency density, which is what H3 asks about");
    println!(
        "  {:<26}{:>16}{:>21}{:>18}{:>14}",
        "", "carrying only", "the rest", "expected het/kb", "class weight"
    );
    println!(
        "  {:<26}{:>16}{:>21}{:>18}{:>14}",
        "", "the reference", "segregate with", "", ""
    );
    let mut arms: Vec<(&str, &ArmSummary)> =
        vec![("the class off", &without), ("the class on", &with)];
    if let Some(summary) = with_coverage.as_ref() {
        arms.push(("the class on, with depth", summary));
    }
    for (name, arm) in &arms {
        println!(
            "  {:<26}{:>16.4}{:>21}{:>18.3}{:>14.5}",
            name,
            arm.p_invariant,
            format!("Beta({:.3}, {:.3})", arm.density_a, arm.density_b),
            1_000.0 * arm.expected_heterozygosity,
            arm.class_share
        );
    }

    // ---- the cheapest check: how many positions does the class claim? ----------------------
    let class = &fit.duplicated_posterior;
    let expected_positions: f64 = class.iter().map(|p| f64::from(*p)).sum();
    let above_half = class.iter().filter(|p| **p > 0.5).count();
    let above_ninety = class.iter().filter(|p| **p > 0.9).count();
    println!("\n  --- how many positions the class claims ---");
    println!(
        "  fitted share of positions                 {:.5}  (Beta({:.3}, {:.3}) for the share \
         of the panel carrying one)",
        with.class_share, with.carrier_a, with.carrier_b
    );
    println!(
        "  expected positions in the class           {:>12.0} of {positions} — 1 in {:.0}",
        expected_positions,
        positions as f64 / expected_positions.max(1e-9)
    );
    println!(
        "  positions more likely in it than not      {above_half:>12} — 1 in {:.0}",
        positions as f64 / (above_half.max(1) as f64)
    );
    println!(
        "  positions at a posterior above 0.9        {above_ninety:>12} — 1 in {:.0}",
        positions as f64 / (above_ninety.max(1) as f64)
    );

    // ---- the per-accession pass ------------------------------------------------------------
    //
    // One accession at a time: its reads, its relative depth and its own carrier posteriors,
    // held together only for as long as it takes to count them.
    let mut carriers_at: Vec<u16> = vec![0; positions];
    let mut carriers_any: Vec<u16> = vec![0; positions];
    let mut class_carriers: Vec<f32> = vec![0.0; positions];
    let claimed: Vec<usize> = (0..positions).filter(|index| class[*index] > 0.5).collect();
    let mut claimed_alt_share = vec![0.0_f32; claimed.len()];
    let mut claimed_relative = vec![0.0_f32; claimed.len()];
    let mut claimed_carriers = vec![0.0_f32; claimed.len()];

    println!(
        "\n  --- H1, the depth test: is a claimed position at twice its accession's depth? ---"
    );
    println!(
        "  Relative depth is 1.0 at this accession's own normal for that stretch's GC content \
         and 2.0 at twice it."
    );
    println!(
        "  Beside each accession's claimed positions stand the ones the same fit calls \
         heterozygous\n  with the class switched off, which is where the class's claims came \
         from, and the whole\n  genome, which is what a position picked at random reads."
    );
    println!(
        "\n  {:<24}{:>8}{:>9}{:>10}{:>12}{:>12}{:>10}{:>12}{:>12}",
        "accession",
        "reads/pos",
        "read at",
        "claimed",
        "median rel.",
        "in 1.6–2.4",
        "het, off",
        "in 1.6–2.4",
        "background"
    );
    let mut per_accession: Vec<(String, f64, f64, usize, f64, f64, f64)> = Vec::new();
    let mut probe_cell: Vec<(String, f64, f64, u64, f64, u64, f64)> = Vec::new();
    for (s, records) in cohort.iter().enumerate() {
        let reads = reads_at_positions(records, positions);
        let relative = relative_depth(records, grid, window_gc, kept);
        let coverage = records.coverage.as_ref().expect("a coverage summary");
        let width_kb = coverage.windows_to_sum() as f64 * coverage.window_bp().get() as f64 / 1e3;

        let mut background_two_copy = 0_u64;
        let mut background_scored = 0_u64;
        let mut claimed_relatives: Vec<f32> = Vec::new();
        let mut claimed_two_copy = 0_u64;
        let mut expected_claimed = 0.0_f64;
        let (mut heterozygous_here, mut heterozygous_two_copy) = (0_u64, 0_u64);
        // The measurement the duplication probe made on real alignments, taken here on the
        // kept positions so the two are over the same list: a position reading between 35%
        // and 65% non-reference at eight reads or more, inside a window near two copies —
        // and the same near-half positions at ordinary coverage, which is the contrast that
        // says whether the class is following depth at all.
        let (mut near_half_two_copy, mut near_half_one_copy) = (0_u64, 0_u64);
        let (mut claim_in_two_copy, mut claim_in_one_copy) = (0.0_f64, 0.0_f64);

        for index in 0..positions {
            let carrier = fit.genotype_posterior[(index * samples + s) * 3 + 2];
            expected_claimed += f64::from(carrier);
            class_carriers[index] += carrier;
            if reads.non_reference[index] >= 1 {
                carriers_any[index] += 1;
            }
            if reads.non_reference[index] >= 2 {
                carriers_at[index] += 1;
            }
            let value = relative[index];
            if value.is_finite() {
                background_scored += 1;
                if (1.6..=2.4).contains(&value) {
                    background_two_copy += 1;
                }
            }
            if carrier > 0.5 && value.is_finite() {
                claimed_relatives.push(value);
                if (1.6..=2.4).contains(&value) {
                    claimed_two_copy += 1;
                }
            }
            if heterozygous_without_class[index * samples + s] > 0.5 && value.is_finite() {
                heterozygous_here += 1;
                if (1.6..=2.4).contains(&value) {
                    heterozygous_two_copy += 1;
                }
            }
            let depth = reads.depth[index];
            if depth >= 8 && value.is_finite() {
                let share = f64::from(reads.non_reference[index]) / f64::from(depth);
                if (0.35..=0.65).contains(&share) {
                    if (1.6..=2.4).contains(&value) {
                        near_half_two_copy += 1;
                        claim_in_two_copy += f64::from(carrier);
                    } else if (0.6..=1.4).contains(&value) {
                        near_half_one_copy += 1;
                        claim_in_one_copy += f64::from(carrier);
                    }
                }
            }
        }
        for (slot, index) in claimed.iter().enumerate() {
            let carrier = fit.genotype_posterior[(index * samples + s) * 3 + 2];
            if carrier > 0.5 {
                claimed_carriers[slot] += 1.0;
                claimed_relative[slot] += relative[*index];
                let depth = reads.depth[*index];
                if depth > 0 {
                    claimed_alt_share[slot] += reads.non_reference[*index] as f32 / depth as f32;
                }
            }
        }

        let claimed_here = claimed_relatives.len();
        let median_claimed = f64::from(median_of(&mut claimed_relatives));
        let share_claimed = claimed_two_copy as f64 / claimed_here.max(1) as f64;
        let share_background = background_two_copy as f64 / background_scored.max(1) as f64;
        println!(
            "  {:<24}{:>8.1}{:>7.1} kb{:>10}{:>12.2}{:>11.1}%{:>10}{:>11.1}%{:>11.1}%",
            records.sample,
            coverage.median_depth(),
            width_kb,
            claimed_here,
            median_claimed,
            100.0 * share_claimed,
            heterozygous_here,
            100.0 * heterozygous_two_copy as f64 / heterozygous_here.max(1) as f64,
            100.0 * share_background,
        );
        per_accession.push((
            records.sample.clone(),
            f64::from(coverage.median_depth()),
            expected_claimed,
            claimed_here,
            median_claimed,
            share_claimed,
            share_background,
        ));
        probe_cell.push((
            records.sample.clone(),
            f64::from(coverage.median_depth()),
            expected_claimed,
            near_half_two_copy,
            claim_in_two_copy / near_half_two_copy.max(1) as f64,
            near_half_one_copy,
            claim_in_one_copy / near_half_one_copy.max(1) as f64,
        ));
    }

    // ---- the probe's own cell, on the kept positions ---------------------------------------
    println!(
        "\n  --- the same positions the duplication probe counts, and what the class does with \
         them ---"
    );
    println!(
        "  A near-half position reads 35–65% non-reference at eight reads or more. The probe \
         counted\n  those inside windows near two copies and found 150 to 590 in two million \
         (1 in 4,000)."
    );
    println!(
        "  Beside them stands the count the class books for that accession — the sum of its own \
         posteriors\n  that the accession carries an extra copy — with both scaled to two \
         million positions."
    );
    println!(
        "\n  {:<24}{:>10}{:>12}{:>10}{:>12}{:>12}{:>12}{:>12}",
        "accession",
        "reads/pos",
        "the class",
        "per 2 M",
        "near half",
        "per 2 M",
        "taken",
        "near half"
    );
    println!(
        "  {:<24}{:>10}{:>12}{:>10}{:>12}{:>12}{:>12}{:>12}{:>12}",
        "", "", "books", "", "at 2×", "", "", "at 1×", "taken"
    );
    let scale = 2e6 / positions as f64;
    for (name, depth, books, two, claim_two, one, claim_one) in &probe_cell {
        println!(
            "  {name:<24}{depth:>10.1}{books:>12.0}{:>10.0}{two:>12}{:>12.0}{:>11.1}%{one:>12}\
             {:>11.1}%",
            books * scale,
            *two as f64 * scale,
            100.0 * claim_two,
            100.0 * claim_one,
        );
    }

    // ---- H2: the posterior against how many accessions carry the allele ---------------------
    println!("\n  --- H2, the carrier-count test ---");
    println!(
        "  An accession counts as carrying the alternative allele where two or more of its \
         reads show it."
    );
    let silent = claimed
        .iter()
        .filter(|index| carriers_any[**index] == 0)
        .count();
    println!(
        "  {silent} of the {} positions the class claims outright have no accession showing \
         even one\n  read that disagrees with the reference.",
        claimed.len()
    );
    println!(
        "\n  {:>18}{:>14}{:>18}{:>20}{:>18}{:>20}",
        "accessions carrying",
        "positions",
        "the class claims",
        "share of them claimed",
        "carriers it sees",
        "heterozygosity there"
    );
    let bands: [(u16, u16); 10] = [
        (0, 0),
        (1, 1),
        (2, 2),
        (3, 3),
        (4, 5),
        (6, 10),
        (11, 20),
        (21, 40),
        (41, 62),
        (63, u16::MAX),
    ];
    for (low, high) in bands {
        let (mut count, mut mass, mut het, mut seen) = (0_u64, 0.0_f64, 0.0_f64, 0.0_f64);
        for index in 0..positions {
            if (low..=high).contains(&carriers_at[index]) {
                count += 1;
                mass += f64::from(class[index]);
                het += f64::from(heterozygosity_without_class[index]);
                seen += f64::from(class_carriers[index]);
            }
        }
        if count == 0 {
            continue;
        }
        println!(
            "  {:>18}{count:>14}{mass:>18.0}{:>19.2}%{:>18.2}{:>20.0}",
            if low == high {
                format!("{low}")
            } else if high == u16::MAX {
                format!("{low} or more")
            } else {
                format!("{low}–{high}")
            },
            100.0 * mass / count as f64,
            seen / mass.max(1e-9),
            het
        );
    }

    // ---- where the heterozygosity went -------------------------------------------------------
    let total_het: f64 = heterozygosity_without_class
        .iter()
        .map(|v| f64::from(*v))
        .sum();
    let at_claimed: f64 = claimed
        .iter()
        .map(|index| f64::from(heterozygosity_without_class[*index]))
        .sum();
    println!("\n  --- where the class took the heterozygosity from ---");
    println!(
        "  with the class off the fit books {total_het:.0} heterozygous sample-positions; \
         {at_claimed:.0} of them\n  ({:.1}%) sit at the {} positions the class claims outright.",
        100.0 * at_claimed / total_het.max(1e-9),
        claimed.len()
    );

    // ---- the claimed positions themselves ----------------------------------------------------
    if let Some(path) = dump {
        let mut out = String::from(
            "contig\tposition\tclass_posterior\tmismapped_posterior\tcarriers\t\
             carriers_claimed\tmean_alt_share\tmean_relative_depth\theterozygosity_class_off\n",
        );
        for (slot, index) in claimed.iter().enumerate() {
            let claimed_here = claimed_carriers[slot].max(1.0);
            out.push_str(&format!(
                "{}\t{}\t{:.4}\t{:.4}\t{}\t{:.0}\t{:.3}\t{:.3}\t{:.3}\n",
                contigs.entries[kept[*index].contig.get() as usize].name,
                kept[*index].position.get(),
                class[*index],
                fit.noisy_posterior[*index],
                carriers_at[*index],
                claimed_carriers[slot],
                claimed_alt_share[slot] / claimed_here,
                claimed_relative[slot] / claimed_here,
                heterozygosity_without_class[*index],
            ));
        }
        std::fs::write(path, out).expect("the audit dump path is writable");
        println!("\n  the claimed positions   {path}");
    }

    let _ = per_accession;
}

/// `ln P(this sample's depth around here | it has two copies) − ln P(… | one copy)`, per sample
/// per kept position.
///
/// **A window's mean depth is read against this sample's own median**, never an absolute depth,
/// and windows are summed until they hold the 12,000 aligned bases the classification needs —
/// one window at 25× and ten at 2.5×.
///
/// ⚠ **The GC correction is missing, and it is not a small omission.** Coverage varies with GC
/// content by a factor of 1.79 across tomato's range, which is larger than the doubling being
/// looked for, so the reading below is that variation plus the copy number rather than the copy
/// number. `CoverageByWindow` keeps the sample's depth-against-GC curve but not each window's own
/// GC content, so the correction cannot be applied from the stored summary; the report proposes
/// storing it.
///
/// The scatter a relative coverage reading has at `copies` copies is
/// `√(copies² × 0.194² + copies × 313 / aligned bases)`, fitted to two points measured across
/// eight tomato accessions (`duplicated_locus_probe_2026-08-12.md` §4).
fn coverage_odds(
    cohort: &[SampleRecords],
    grid: &CoverageGrid,
    kept: &[pop_var_caller::ng::types::GenomePosition],
) -> Vec<Arc<[f32]>> {
    let mut out = Vec::with_capacity(cohort.len());
    for sample in cohort {
        let Some(coverage) = sample.coverage.as_ref() else {
            return Vec::new();
        };
        let span = coverage.windows_to_sum();
        let aligned =
            f64::from(coverage.median_depth()) * coverage.window_bp().get() as f64 * span as f64;
        let sd = |copies: f64| {
            (copies * copies * 0.194 * 0.194 + copies * 313.0 / aligned.max(1.0)).sqrt()
        };
        // One relative reading per window, so a position only costs a lookup.
        let mut per_window = vec![0.0_f32; coverage.len()];
        for slot in 0..coverage.len() {
            let first = slot.saturating_sub(span / 2);
            let last = (first + span).min(coverage.len());
            let mean = f64::from(coverage.mean_depth_over(first..last));
            let expected = f64::from(coverage.median_depth()).max(1e-6);
            let reading = (mean / expected).clamp(0.02, 8.0);
            let ln_normal = |copies: f64| {
                let s = sd(copies);
                -0.5 * ((reading - copies) / s).powi(2) - s.ln()
            };
            per_window[slot] = (ln_normal(2.0) - ln_normal(1.0)) as f32;
        }
        let odds: Vec<f32> = kept
            .iter()
            .map(|position| {
                grid.slot_of(position.contig, position.position)
                    .map_or(0.0, |slot| per_window[slot])
            })
            .collect();
        out.push(Arc::from(odds));
    }
    out
}

/// One sample: one walk, filling the records and the coverage summary together.
///
/// Returns the records and, beside them, **each coverage window's reference GC fraction** —
/// `f32::NAN` where the window holds no called base. That is a property of the reference and
/// the grid rather than of the sample, so every walk returns the same list; the stored summary
/// keeps the sample's depth-against-GC curve but not the window's own GC, so without this the
/// curve cannot be looked up and a relative depth reading is copy number times GC content
/// (`duplicated_locus_probe_2026-08-12.md` §2).
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
    generic_domain: &SelectableRegions,
    grid: &CoverageGrid,
    kept: &pop_var_caller::ng::parameter_estimation::joint::loci::KeptLoci,
    identity: &SelectionIdentity,
) -> (SampleRecords, Vec<f32>) {
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

    // The coverage denominator and the GC content come from the reference over the same
    // generic stretch the walk covers, not from what the walk emitted: a window whose reads
    // are missing has a low mean depth, and that is the truth about it.
    let mut coverage = CoverageAccumulator::new(grid.clone());
    let reader = accessor();
    let mut bases = Vec::new();
    for span in generic_domain.spans() {
        reader
            .fetch_into(span.contig, span.start.get(), span.len(), &mut bases)
            .expect("a generic span reads from the reference");
        coverage.observe_reference(*span, &bases);
    }

    let contig_of = |name: &str| {
        contigs
            .entries
            .iter()
            .position(|entry| entry.name == name)
            .map(|i| ContigId(i as u32))
    };
    let mut writer = RecordWriter::new(
        sample.sample.to_string(),
        kept,
        sample.read_groups.clone(),
        &contig_of,
        identity.clone(),
        DepthBinEdges::new(),
        ReadCap(pop_var_caller::ng::locus_generation::ssr::DEFAULT_SSR_MAX_READS_PER_LOCUS),
        // **The ladder's own top.** The stored code is a bin and the ladder's last bin is the
        // last there is, so a position deeper than that is thinned down to it before anything
        // is recorded — otherwise a 300× sample records a depth of about 111 beside an
        // undiminished count of alternative reads.
        DepthCap(DepthBinEdges::new().max_depth()),
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
            LocusKind::Generic => {
                generic_loci += 1;
                coverage.add_depths(
                    locus.region.contig,
                    locus.region.start,
                    &locus.num_obs_along_locus(),
                );
            }
            LocusKind::Ssr(_) => ssr_loci += 1,
            _ => {}
        }
        writer.add_locus(&locus);
    }
    println!(
        "  walked         {generic_loci} generic loci, {ssr_loci} STR loci; median window depth \
         {:.2} reads a position, {} windows summed to reach 12,000 aligned bases",
        coverage.median_depth(),
        (12_000.0 / (f64::from(coverage.median_depth()) * WINDOW_BP.get() as f64))
            .ceil()
            .max(1.0)
    );
    let window_gc: Vec<f32> = (0..coverage.grid().len())
        .map(|slot| coverage.gc_fraction(slot).unwrap_or(f32::NAN))
        .collect();
    (writer.finish(Some(coverage.finish())), window_gc)
}

/// What one sample's records weigh, each object on its own line.
///
/// **Separately and not as a total**, because `parameter_prepass_joint_records.md` §6 claims
/// the STR set is the larger half and a single number would hide that being wrong.
fn report_sizes(records: &SampleRecords) {
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
            pop_var_caller::ng::parameter_estimation::joint::records::AlleleObservation,
        >();

    let ssr_loci = records.ssr.values().next().map_or(0, |s| s.len());
    let ssr_dense_bytes: usize = records
        .ssr
        .values()
        .map(|s| {
            s.len()
                * (std::mem::size_of::<
                    pop_var_caller::ng::parameter_estimation::joint::records::OffsetCounts,
                >() + std::mem::size_of::<u16>()
                    + std::mem::size_of::<u32>()
                    + std::mem::size_of::<bool>())
        })
        .sum();
    let guard_entries: usize = records.ssr.values().map(|s| s.guard().len()).sum();
    let difference_entries: usize = records.ssr.values().map(|s| s.differences().len()).sum();
    let difference_bytes = difference_entries
        * std::mem::size_of::<
            pop_var_caller::ng::parameter_estimation::joint::records::TractDifference,
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
    if let Some(coverage) = &records.coverage {
        println!(
            "  coverage       {:.3} MB over {} windows of {} bp",
            coverage.len() as f64 / 1e6,
            coverage.len(),
            coverage.window_bp().get()
        );
    }
    println!(
        "  states         {never_walked} never walked, {zero_depth} walked at zero depth, \
         {covered} with reads"
    );
    println!(
        "  TOTAL          {:.3} MB held for this sample",
        mb(depth_bytes + sparse_bytes + ssr_dense_bytes + difference_bytes)
            + records
                .coverage
                .as_ref()
                .map_or(0.0, |c| c.len() as f64 / 1e6)
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

fn structure(cohort: &[SampleRecords], positions: usize) {
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
fn markers(cohort: &[SampleRecords], positions: usize) -> Vec<Marker> {
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
                    observation.reads;
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
                        .saturating_add(observation.reads.min(u32::from(u16::MAX)) as u16);
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
fn per_read_group(records: &SampleRecords) -> BTreeMap<String, usize> {
    records
        .generic
        .iter()
        .map(|(group, g)| (format!("{group:?}"), g.non_reference().len()))
        .collect()
}
