//! **Does excluding the mismapped positions leave a real 3% contamination still findable?**
//!
//! On the 63 tomato accessions the contamination estimate came back with a median of 6.5% and a
//! highest of 12.5%, which cannot be true of sixty-three plants from a public archive: it is a
//! floor, not a measurement. The named cause is that **a position where two stretches of genome
//! the reference holds once both pile reads up produces the contamination signature in every
//! sample at once** — a small share of reads carrying an allele the sample should not have —
//! and nothing excluded such positions from the ones contamination was measured over.
//!
//! The fit already computes, for every position, the probability that it is one of those. This
//! program is the control that says whether using it works, and it needs both halves:
//!
//! - **the floor** — a drawn panel *with* mismapped positions planted and nobody contaminated,
//!   which must come back near zero once they are excluded and does not today;
//! - **the signal** — the same panel with one sample genuinely contaminated at 3%, which must
//!   still be found. An estimator that has been made blind returns zero on both, and that is
//!   the trap this project has already fallen into once.
//!
//! Everything runs through the library: the cohort is drawn as the records a walk would have
//! written, [`fit_jointly`] fits it and produces the per-position probabilities, and
//! [`fit_contamination`] is then called once per arm on the same fit. So what is compared is
//! the shipped code under different settings, not a re-implementation of it.
//!
//! # The arms
//!
//! Two settings crossed. **Which positions are markers**: all of them, as before, against
//! dropping those the fit says are more likely mismapped than not. And **where the sample is
//! taken to stand on the panel's axes of variation** while its fraction is searched for: as the
//! decomposition read it off the sample's own reads, or divided by `1 − α` — which undoes
//! exactly the drag a fraction `α` of stray reads causes — or each axis searched freely beside
//! `α`.
//!
//! ```text
//! ng_joint_contamination_control [positions] [samples] [depth] [groups] [fst] [noisy-share] [noisy-rate] [alpha]
//! ```
//!
//! # Several libraries from one plant — added 2026-08-20
//!
//! **A plant's DNA is prepared into libraries, and a library is what a second plant's DNA gets
//! into.** Two libraries made from one plant can therefore carry different stray-read fractions —
//! one contaminated and the other clean — which one number for the whole plant cannot say. The
//! estimator produces one number per plant today, and
//! `doc/devel/ng/impl_plan/contamination_read_group_grain.md` is moving it to the library. No
//! cohort we hold can show the difference, because every sample of every benchmark here was
//! sequenced from one library, so the case has to be drawn.
//!
//! Set `LIBRARIES` to divide each plant's reads into that many libraries — the plant's total
//! depth is unchanged, so what varies is how its reads are divided, not how many it has. By
//! default the whole stray-read fraction goes into the first library and the others are clean,
//! which is the case the two grains disagree about most; `LIBRARY_ALPHAS` and `LIBRARY_DEPTHS`
//! set the fractions and the read shares by hand.
//!
//! **`LIBRARIES` unset draws exactly what it always drew**, generator call for generator call,
//! so the numbers already published from this program can be reproduced from the same command
//! line — checked, and the whole table matched.

use std::collections::BTreeMap;
use std::time::Instant;

use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, CohortCensusEvidence, DepthCap, DepthCode, DepthLadderDigest,
    GenericEvidence, ObservedAllele, PackedDepthCodes, ReadCap, RecordingTerms,
    SampleCensusEvidence, Section, SectionKey, SelectionTermsDigest,
};
use pop_var_caller::ng::parameter_estimation::joint::contamination::{
    ContaminationConfig, ContaminationEstimate, ContaminationGrain, OwnCoordinates,
    fit_contamination,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{
    FrequencyDensity, JointFit, JointFitConfig, fit_jointly,
};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::ReadGroupId;

/// How often a read misreads a base at an ordinary position.
const CLEAN: f64 = 0.003;

/// How much less heterozygous than random mating each drawn plant is — a selfing crop.
const HOM_EXCESS: f64 = 0.6;

/// One library of every drawn plant: what share of that plant's reads were prepared in it, and
/// what share of *those* reads came from another plant.
///
/// **A plant's DNA is prepared into libraries and a library is what gets contaminated**, so the
/// stray-read fraction belongs here and not on the plant. One library holding all the reads and
/// carrying the whole fraction is the shape every measurement before 2026-08-20 was made at, and
/// it is what `LIBRARIES` unset still draws.
#[derive(Copy, Clone, Debug)]
struct Library {
    /// This library's share of the plant's reads. The shares are used as given and are expected
    /// to sum to one; anything else changes the plant's depth, which is a different experiment.
    depth_share: f64,
    /// The fraction of this library's reads that came from another plant — applied to sample 0,
    /// which is the one the panel is asked to find.
    alpha: f64,
}

/// How the drawn plants are divided into libraries, from the environment.
///
/// Three knobs, all optional, because the positional arguments are this program's published
/// interface and adding to them would change what an old command line means:
///
/// - `LIBRARIES=2` — how many libraries each plant's reads are prepared into.
/// - `LIBRARY_ALPHAS=0.06,0` — what fraction of each library's reads came from another plant,
///   for sample 0. **Default: the whole fraction in the first library and none in the rest**,
///   which is the case the two grains disagree about most.
/// - `LIBRARY_DEPTHS=0.5,0.5` — each library's share of the plant's reads. Default: equal.
///
/// **The plant's total depth is the same however many libraries it has**, so what a sweep varies
/// is how one plant's reads are divided rather than how many it has.
fn libraries_from_environment(alpha: f64) -> Vec<Library> {
    let count: usize = std::env::var("LIBRARIES")
        .ok()
        .map_or(1, |v| v.parse().expect("LIBRARIES is a count"));
    assert!(count >= 1, "a plant is sequenced from at least one library");
    let list = |name: &str| -> Option<Vec<f64>> {
        std::env::var(name).ok().map(|v| {
            let parsed: Vec<f64> = v
                .split(',')
                .map(|field| field.trim().parse().expect("a number"))
                .collect();
            assert_eq!(
                parsed.len(),
                count,
                "{name} gives {} values for {count} libraries",
                parsed.len()
            );
            parsed
        })
    };
    let alphas = list("LIBRARY_ALPHAS").unwrap_or_else(|| {
        let mut all = vec![0.0; count];
        all[0] = alpha;
        all
    });
    let depths = list("LIBRARY_DEPTHS").unwrap_or_else(|| vec![1.0 / count as f64; count]);
    (0..count)
        .map(|k| Library {
            depth_share: depths[k],
            alpha: alphas[k],
        })
        .collect()
}

fn main() {
    let mut args = std::env::args().skip(1);
    let positions: usize = args.next().map_or(400_000, |a| a.parse().expect("a count"));
    let samples: usize = args.next().map_or(40, |a| a.parse().expect("a count"));
    let depth: f64 = args.next().map_or(3.0, |a| a.parse().expect("a depth"));
    let groups: usize = args.next().map_or(4, |a| a.parse().expect("a count"));
    let fst: f64 = args.next().map_or(0.20, |a| a.parse().expect("an F_st"));
    let noisy_share: f64 = args.next().map_or(0.033, |a| a.parse().expect("a share"));
    let noisy: f64 = args.next().map_or(0.024, |a| a.parse().expect("a rate"));
    let alpha: f64 = args.next().map_or(0.03, |a| a.parse().expect("a fraction"));
    let components: usize = args.next().map_or(4, |a| a.parse().expect("a count"));

    println!("positions        {positions}");
    println!("samples          {samples} in {groups} subpopulations at F_st {fst}");
    println!("depth            {depth} reads a position");
    println!(
        "mismapped        {noisy_share} of positions, where a read disagrees at {noisy} instead \
         of {CLEAN}"
    );
    println!("contaminated     sample 0, at {alpha} of its reads, by a plant drawn from the panel");

    // **One panel is one draw, and at a few thousand markers a draw is not a small thing.** The
    // fraction returned for the same planted 3% moved between 0.0121 and 0.0278 across three
    // seeds at 40,000 positions, so a comparison between two settings has to be made either at
    // the same seed or over several. `SEED` is what lets a caller do the second.
    let seed: u64 = std::env::var("SEED").ok().map_or(
        0x9E37_79B9_7F4A_7C15,
        // Hexadecimal as well as decimal, because the default is written that way and a caller
        // varying it will copy that form.
        |v| match v.strip_prefix("0x").or_else(|| v.strip_prefix("0X")) {
            Some(digits) => u64::from_str_radix(digits, 16).expect("a hexadecimal seed"),
            None => v.parse().expect("a decimal seed"),
        },
    );
    println!("seed             {seed:#x}");

    let libraries = libraries_from_environment(alpha);
    if libraries.len() > 1 {
        let weighted: f64 = libraries.iter().map(|l| l.depth_share * l.alpha).sum();
        println!(
            "libraries        {} per plant, sharing its reads {}",
            libraries.len(),
            libraries
                .iter()
                .map(|l| format!("{:.3}", l.depth_share))
                .collect::<Vec<_>>()
                .join("/")
        );
        println!(
            "                 sample 0's stray-read fraction per library {} — a per-sample \
             estimate can at best return the read-weighted mean, {weighted:.4}",
            libraries
                .iter()
                .map(|l| format!("{:.3}", l.alpha))
                .collect::<Vec<_>>()
                .join("/")
        );
    }
    println!(
        "\nBoth halves are needed. A panel with nobody contaminated says whether the floor is\n\
         gone; the same panel with one sample at {alpha} says whether the estimator can still\n\
         see anything at all."
    );

    // **`SWEEP=1` runs one scenario and one arm instead of sixteen**, and it is the one that
    // ships: mismapped positions planted and dropped, the plant's coordinates as the
    // decomposition read them, the depth summed over its range. A sweep over depths and library
    // counts wants many *panels* at the settings that ship, where the four-scenario crossing
    // above wants one panel at many settings, and the fit is the expensive half of both.
    let sweeping = std::env::var("SWEEP").as_deref() == Ok("1");
    let mut scenarios = if sweeping {
        vec![(noisy_share, alpha)]
    } else {
        vec![
            (noisy_share, alpha),
            (noisy_share, 0.0),
            (0.0, alpha),
            (0.0, 0.0),
        ]
    };
    scenarios.dedup();
    for (planted, spike) in scenarios {
        let label = format!(
            "{}, {}",
            if planted > 0.0 {
                format!("{planted} of positions mismapped")
            } else {
                "no mismapped positions at all".to_string()
            },
            if spike > 0.0 {
                format!("sample 0 contaminated at {spike}")
            } else {
                "nobody contaminated".to_string()
            }
        );
        println!("\n=== {label} ===");
        let at = Instant::now();
        // A scenario with nobody contaminated zeroes every library's fraction; one with a spike
        // uses the fractions as configured. The library *shares of the reads* never change, so
        // the two scenarios differ in the one thing they are meant to.
        let scenario_libraries: Vec<Library> = libraries
            .iter()
            .map(|library| Library {
                alpha: if spike > 0.0 { library.alpha } else { 0.0 },
                ..*library
            })
            .collect();
        let drawn = draw(
            samples,
            positions,
            depth,
            planted,
            noisy,
            groups,
            fst,
            &scenario_libraries,
            seed,
        );
        let config = JointFitConfig {
            quadrature_nodes: 12,
            // **The third class is off here**, and this harness is not the place to judge it:
            // it is measured in `ng_joint_duplicated_in_fit.rs`. Left on, every fit below pays
            // for a second frequency integral and the numbers move for a reason that has
            // nothing to do with contamination.
            duplicated_positions: false,
            ..JointFitConfig::default()
        };
        let mut cohort = CohortCensusEvidence::new(drawn.samples.clone())
            .expect("a drawn cohort records one way");
        let fit = fit_jointly(&mut cohort, &config).expect("a drawn cohort pools");
        println!(
            "  the fit: {} of positions booked mismapped against {planted} planted, error rate \
             {:.5} against {CLEAN}, {} passes, {:.0} s",
            format_args!("{:.4}", fit.noisy_share),
            fit.noise[&ReadGroupId(0)].value.clean,
            fit.passes,
            at.elapsed().as_secs_f64()
        );
        let condemned = fit.noisy_posterior.iter().filter(|p| **p > 0.5).count();
        println!(
            "  {condemned} of {} positions are more likely mismapped than not",
            fit.noisy_posterior.len()
        );

        println!(
            "\n  {:<58}{:>9}{:>9}{:>36}{:>9}{:>9}{:>9}",
            "", "grain", "markers", "sample 0's libraries", "median", "worst", "refused"
        );
        println!(
            "  {:<58}{:>9}{:>9}{:>36}{:>9}{:>9}{:>9}",
            "", "", "", "", "clean", "clean", ""
        );
        let mismapped_arms: &[bool] = if sweeping { &[false] } else { &[true, false] };
        let coordinate_arms: &[(OwnCoordinates, bool)] = if sweeping {
            &[(OwnCoordinates::AsRead, true)]
        } else {
            &[
                (OwnCoordinates::AsRead, false),
                (OwnCoordinates::AsRead, true),
                (OwnCoordinates::UndoneByAlpha, true),
                (OwnCoordinates::MaximisedFreely, true),
            ]
        };
        for &keep_mismapped in mismapped_arms {
            for &(own, integrate) in coordinate_arms {
                let label = format!(
                    "{}, {}{}",
                    if keep_mismapped {
                        "every position"
                    } else {
                        "mismapped dropped"
                    },
                    match own {
                        OwnCoordinates::AsRead => "coordinates as read",
                        OwnCoordinates::UndoneByAlpha => "coordinates undone by α",
                        OwnCoordinates::MaximisedFreely => "coordinates free",
                    },
                    if integrate { ", depth summed over" } else { "" }
                );
                // **Both grains on the same drawn panel**, which is the comparison this
                // program exists for since 2026-08-20: the reads are identical and the only
                // difference is whether one fraction is fitted from each library's reads or
                // one from all of a plant's.
                for grain in [ContaminationGrain::ReadGroup, ContaminationGrain::Sample] {
                    let settings = ContaminationConfig {
                        components,
                        max_noisy_posterior: if keep_mismapped { 1.0 } else { 0.5 },
                        weight_by_posterior: !keep_mismapped,
                        own_coordinates: own,
                        integrate_over_depth_bin: integrate,
                        // `LEAVE_SELF_OUT=1` runs every arm with each sample taken out of the
                        // frequency it is judged against. Off by default, as it ships, so the
                        // measurement that put it there can be repeated from either side.
                        leave_self_out: std::env::var("LEAVE_SELF_OUT").as_deref() == Ok("1"),
                        grain,
                    };
                    let arm = run(&drawn.samples, &fit, &settings);
                    // **A refused library is not a clean one, and it is not a zero either.**
                    // It comes back `NaN`; sorting with `total_cmp` puts those last, so a
                    // single refusal used to be printed as the worst clean library's fraction.
                    // They are counted instead, because how many a panel is too small to answer
                    // for is itself the finding.
                    let mut clean: Vec<f64> =
                        arm.clean.iter().copied().filter(|a| !a.is_nan()).collect();
                    let refused = arm.clean.len() - clean.len()
                        + arm.spiked.iter().filter(|a| a.is_nan()).count();
                    clean.sort_by(f64::total_cmp);
                    println!(
                        "  {:<58}{:>9}{:>9}{:>36}{:>9.4}{:>9.4}{refused:>9}",
                        if grain == ContaminationGrain::ReadGroup {
                            label.as_str()
                        } else {
                            ""
                        },
                        if grain == ContaminationGrain::ReadGroup {
                            "library"
                        } else {
                            "plant"
                        },
                        arm.markers,
                        arm.spiked
                            .iter()
                            .map(|a| format!("{a:.4}"))
                            .collect::<Vec<_>>()
                            .join("/"),
                        clean[clean.len() / 2],
                        clean.last().copied().unwrap_or(0.0),
                    );
                }
            }
        }
    }
}

/// One arm's answers, split into the plant the panel is asked to find and all the others.
struct Arm {
    markers: u64,
    /// Sample 0's libraries in order — **the row the two grains disagree on**, because at the
    /// plant grain they are all one number and at the library grain they need not be.
    spiked: Vec<f64>,
    /// Every library of every other plant, flattened.
    clean: Vec<f64>,
}

/// One arm: how many markers survived, and every library's fraction.
fn run(samples: &[SampleCensusEvidence], fit: &JointFit, settings: &ContaminationConfig) -> Arm {
    // **One error rate per read group**, which is the grain the fit produced them at. Library
    // slot `k` is read group `k` for every plant, so this is one rate a slot.
    let error: BTreeMap<ReadGroupId, f64> = fit
        .noise
        .iter()
        .map(|(group, estimate)| (*group, estimate.value.clean))
        .collect();
    let excess: Vec<f64> = samples
        .iter()
        .map(|sample| fit.hom_excess[&sample.sample].value.get())
        .collect();
    let mut cohort =
        CohortCensusEvidence::new(samples.to_vec()).expect("a drawn cohort records one way");
    let estimates = fit_contamination(
        &mut cohort,
        &DepthBinEdges::for_census(),
        &error,
        &excess,
        &fit.noisy_posterior,
        settings,
    )
    .expect("a resident census has no file to fail on");
    let alpha_of = |estimate: &ContaminationEstimate| match estimate {
        ContaminationEstimate::Estimated { alpha, .. } => *alpha,
        ContaminationEstimate::NotIdentified { .. } => f64::NAN,
    };
    Arm {
        markers: estimates
            .iter()
            .flatten()
            .find_map(|(_, estimate)| match estimate {
                ContaminationEstimate::Estimated { panel_markers, .. } => Some(*panel_markers),
                ContaminationEstimate::NotIdentified { .. } => None,
            })
            .unwrap_or(0),
        spiked: estimates[0]
            .iter()
            .map(|(_, estimate)| alpha_of(estimate))
            .collect(),
        clean: estimates[1..]
            .iter()
            .flatten()
            .map(|(_, estimate)| alpha_of(estimate))
            .collect(),
    }
}

struct Drawn {
    samples: Vec<SampleCensusEvidence>,
}

/// Draw a structured cohort and write it into the records the fit reads.
///
/// **The contaminant is drawn from the whole panel**, not from the contaminated sample's own
/// subpopulation, because a neighbouring library on a plate is whoever else was on the plate.
#[allow(
    clippy::too_many_arguments,
    reason = "the drawn cohort's own parameters"
)]
fn draw(
    samples: usize,
    positions: usize,
    depth: f64,
    noisy_share: f64,
    noisy: f64,
    groups: usize,
    fst: f64,
    libraries: &[Library],
    seed: u64,
) -> Drawn {
    // A Beta(0.7, 2.5) population, and the same three branches the fit's own density has: most
    // positions carry one allele, a few are fixed on a non-reference one, the rest segregate.
    let (a, b) = (0.7_f64, 2.5_f64);
    let density = FrequencyDensity {
        p_invariant: 0.97,
        p_fixed_alt: 0.0005,
        a,
        b,
    };
    let mut rng = Rng(seed);
    let edges = DepthBinEdges::for_census();
    // **One depth code list and one observation list per (plant, library)**, because that is
    // what the census stores: a section is one read group's evidence, and a plant's libraries
    // are separate sections of its record.
    let section_count = samples * libraries.len();
    let mut codes: Vec<PackedDepthCodes> = (0..section_count)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); section_count];
    let slot = |sample: usize, library: usize| sample * libraries.len() + library;
    let group_of = |s: usize| s * groups / samples;

    for index in 0..positions {
        let rate = if rng.uniform() < noisy_share {
            noisy
        } else {
            CLEAN
        };
        let branch = rng.pick(&[
            density.p_invariant,
            density.p_fixed_alt,
            1.0 - density.p_invariant - density.p_fixed_alt,
        ]);
        let allele = 1 + ((rng.uniform() * 3.0) as usize).min(2);
        // Balding–Nichols: one ancestral frequency, then a frequency per subpopulation drawn
        // around it with a spread set by F_st.
        let ancestral = match branch {
            0 => 0.0,
            1 => 1.0,
            _ => rng.beta(a, b),
        };
        let per_group: Vec<f64> = (0..groups)
            .map(|_| {
                if branch != 2 || fst <= 0.0 {
                    ancestral
                } else {
                    let scale = (1.0 - fst) / fst;
                    rng.beta(
                        (ancestral * scale).max(1e-3),
                        ((1.0 - ancestral) * scale).max(1e-3),
                    )
                }
            })
            .collect();

        let genotype = |rng: &mut Rng, f: f64| -> usize {
            let het = 2.0 * f * (1.0 - f) * (1.0 - HOM_EXCESS);
            let shift = HOM_EXCESS * f * (1.0 - f);
            rng.pick(&[(1.0 - f) * (1.0 - f) + shift, het, f * f + shift])
        };
        let genotypes: Vec<usize> = (0..samples)
            .map(|s| match branch {
                0 => 0,
                1 => 2,
                _ => genotype(&mut rng, per_group[group_of(s)]),
            })
            .collect();

        for sample in 0..samples {
            // **Each library draws its own reads at its own share of the plant's depth, and
            // stray reads enter one library at a time.** With one library holding all the reads
            // this is the draw it always was, right down to the order the generator is consumed
            // in, which is what lets a one-library run reproduce the published numbers.
            for (library, settings) in libraries.iter().enumerate() {
                let reads = rng.poisson(depth * settings.depth_share);
                let mut counts = [0_u32; 5];
                for _ in 0..reads {
                    let from =
                        if sample == 0 && settings.alpha > 0.0 && rng.uniform() < settings.alpha {
                            genotypes[(rng.uniform() * samples as f64) as usize % samples]
                        } else {
                            genotypes[sample]
                        };
                    let carried = from as f64 / 2.0;
                    let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
                    let on_reference = (1.0 - carried) * (1.0 - rate) + carried * rate / 3.0;
                    let u = rng.uniform();
                    let code = if u < on_candidate {
                        allele
                    } else if u < on_candidate + on_reference {
                        0
                    } else if u < on_candidate + on_reference + rate / 3.0 {
                        (allele % 3) + 1
                    } else {
                        ((allele + 1) % 3) + 1
                    };
                    counts[code] += 1;
                }
                let here = slot(sample, library);
                codes[here].set(index, DepthCode::Binned(edges.bin_for(reads)));
                for (code, count) in counts.iter().enumerate() {
                    if code == 0 || *count == 0 {
                        continue;
                    }
                    sparse[here].push(AlleleObservation {
                        index: index as u32,
                        allele: match code {
                            1 => ObservedAllele::C,
                            2 => ObservedAllele::G,
                            _ => ObservedAllele::T,
                        },
                        reads: u8::try_from(*count)
                            .expect("a drawn count fits the census's one-byte field"),
                    });
                }
            }
        }
    }

    let terms = RecordingTerms {
        selection: SelectionTermsDigest::of(&SelectionTerms {
            seed,
            reference: ReferenceDigest([7; 16]),
            analysed_regions: RegionSetDigest([9; 16]),
            catalog_built_under: CatalogBuildSettings {
                criteria: StrRepeatCriteria::default(),
                scan: ScanParams::default(),
                tool_version: "0.1.0".to_string(),
            },
            ssr_criteria: StrRepeatCriteria::default(),
            generic_target: positions as u64,
            ssr_cap: 1_000,
        }),
        kept_loci: CensusLociDigester::new().finish(),
        ssr_stratum_counts: Default::default(),
        read_cap: ReadCap(1_000),
        depth_ladder: DepthLadderDigest::of(&DepthBinEdges::for_census()),
        depth_cap: DepthCap::new(124),
    };
    Drawn {
        samples: (0..samples)
            .map(|s| {
                // **Library `k` of every plant is read group `k`**, shared across the panel
                // rather than one identifier per plant per library. Real read groups are unique
                // to a plant, but sharing them here keeps the error rate fitted from the whole
                // panel exactly as it is today, so a difference in the fraction is the split of
                // the fraction and not a change in how well the error rate is determined.
                let sections = (0..libraries.len())
                    .map(|k| {
                        let here = slot(s, k);
                        (
                            SectionKey::Generic(ReadGroupId(k as u32)),
                            Section::Generic(GenericEvidence::from_parts(
                                std::mem::replace(
                                    &mut codes[here],
                                    PackedDepthCodes::never_walked(0),
                                ),
                                std::mem::take(&mut sparse[here]),
                            )),
                        )
                    })
                    .collect();
                SampleCensusEvidence::resident(format!("s{s:02}"), terms.clone(), sections)
            })
            .collect(),
    }
}

/// The stream every drawn number comes from. Deterministic, so a run is reproducible.
struct Rng(u64);

impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn pick(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut cut = self.uniform() * total;
        for (index, weight) in weights.iter().enumerate() {
            cut -= weight;
            if cut <= 0.0 {
                return index;
            }
        }
        weights.len() - 1
    }

    fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a);
        let y = self.gamma(b);
        (x / (x + y)).clamp(1e-9, 1.0 - 1e-9)
    }

    fn gamma(&mut self, shape: f64) -> f64 {
        if shape < 1.0 {
            let u = self.uniform().max(1e-12);
            return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let z = self.normal();
            let v = (1.0 + c * z).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.uniform().max(1e-12);
            if u.ln() < 0.5 * z * z + d - d * v + d * v.ln() {
                return d * v;
            }
        }
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-12);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    fn poisson(&mut self, mean: f64) -> u32 {
        let limit = (-mean).exp();
        let mut product = self.uniform();
        let mut count = 0;
        while product > limit && count < 400 {
            count += 1;
            product *= self.uniform();
        }
        count
    }
}
