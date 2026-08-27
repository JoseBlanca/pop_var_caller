//! **What holding a whole cohort's pre-pass results at once costs, in bytes** — measured, at
//! cohort sizes from one sample to a thousand.
//!
//! ```text
//! ./scripts/dev.sh cargo run --profile profiling --example ng_prepass_handover_footprint \
//!     --features dhat-heap
//! ```
//!
//! # The question, and why it is asked here
//!
//! The function that turns the parameter pre-pass's outputs into what calling reads
//! (`RunParameters::from_prepass`) is handed **every sample's results at the same moment**,
//! because the object it builds covers the whole run: one calibration per library, one inbreeding
//! coefficient per sample, one substitution rate per (library, tract shape). The project commits
//! the caller to cohorts from one sample to several thousand, so what that costs at a thousand is
//! a number somebody has to have before a run driver is written. **If it is large, that is a
//! finding for the driver's plan** — where a sample's results are produced and could be released
//! — and not something the seam can fix.
//!
//! # What is measured
//!
//! Live bytes on the heap, from `dhat`'s allocator, in three parts that add up to the peak:
//!
//! 1. **the inputs** — every sample's SNP/indel parameters and repeat-tract parameters, held at
//!    once, which is what a caller of the seam must have in hand before it can call;
//! 2. **the union** — the run-wide maps the seam builds out of them: one error rate and one
//!    minted-error total per library, one substitution rate per (library, tract shape), one
//!    coefficient per sample;
//! 3. **what assembling adds** — the dense per-library vectors calling indexes. It is small
//!    because the maps above are *moved* into the result rather than copied.
//!
//! The three are summed as **the peak**, which is what the process holds at the moment the seam
//! returns: the seam does not consume its inputs, so both are live at once. Each part is the
//! difference between two readings of `curr_bytes`, so every figure is what is *live* rather than
//! what was allocated and freed.
//!
//! # The shape of a sample, and where the numbers come from
//!
//! A footprint is meaningless without saying what a sample holds, and two of the three counts are
//! measured elsewhere in this project rather than chosen here:
//!
//! - **one library a sample.** Every sample of both benchmark cohorts here — the 63-accession
//!   tomato panel and the GIAB human trio — was sequenced once.
//! - **338 repeat-tract strata a library**, from the repeat-tract fit's own report
//!   (`doc/devel/reports/implementations/ng_parameter_prepass_ssr_e5_2026-08-13.md`), which is the
//!   larger of the two counts this project has measured; the tomato SL4.00 catalogue holds 141
//!   (`census_tract_grain_b4_2026-08-14.md`).
//! - **78 allele-length genotypes a stratum**, which is what a five-copy dinucleotide tract has at
//!   two genome copies: twelve allele lengths make 78 unordered pairs. Longer tracts have more, so
//!   this is a floor rather than a typical value.
//!
//! All three are arguments, so a reader who disagrees can re-run rather than rescale by hand:
//!
//! ```text
//! cargo run … --example ng_prepass_handover_footprint --features dhat-heap -- 338 1 78
//! ```
//!
//! # What it does not measure
//!
//! Nothing here runs a pre-pass or reads a file. The values are synthesised at the shape above,
//! so this is a **footprint** measurement and not an accuracy one: no fitted number in it means
//! anything.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::collections::{BTreeMap, BTreeSet};

use smallvec::SmallVec;

use pop_var_caller::ng::calling::genotype_prior::{SeedRegime, SpectrumSeed};
use pop_var_caller::ng::calling::run_parameters::RunParameters;
use pop_var_caller::ng::parameter_estimation::generic::GenericSampleParameters;
use pop_var_caller::ng::parameter_estimation::generic::calibration::MintedReadErrors;
use pop_var_caller::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use pop_var_caller::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use pop_var_caller::ng::parameter_estimation::ssr::slippage::SlippageModel;
use pop_var_caller::ng::parameter_estimation::ssr::{
    GenotypeFrequency, RepeatCount, SsrSampleParameters, Stratum, StratumFit, StratumFitSummary,
    StratumKey, WholeRepeatOffset,
};
use pop_var_caller::ng::parameter_estimation::{Estimate, Provenance};
use pop_var_caller::ng::types::{ErrorRate, InbreedingF, Ploidy, ReadGroupId, SsrPeriod};

/// The cohort sizes reported, from the smallest run the caller must serve to the largest the
/// project commits to.
const COHORT_SIZES: [usize; 4] = [1, 10, 100, 1_000];

/// Strata a library holds, unless the command line says otherwise — see this file's header.
const STRATA_A_LIBRARY: usize = 338;

/// Libraries a sample holds, unless the command line says otherwise.
const LIBRARIES_A_SAMPLE: usize = 1;

/// Allele-length genotypes a stratum's fit weighed, unless the command line says otherwise.
const GENOTYPES_A_STRATUM: usize = 78;

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("two genome copies")
}

/// One library's SNP/indel parameters: a fitted rate and a minted-error total for each of its
/// read groups.
fn generic_parameters_of(read_groups: &[ReadGroupId]) -> GenericSampleParameters {
    GenericSampleParameters {
        error_rate: read_groups
            .iter()
            .map(|&group| {
                (
                    group,
                    Estimate {
                        value: ErrorRate::try_new(0.002).expect("a legal rate"),
                        provenance: Provenance::FittedHere,
                        observations: 1_000_000,
                    },
                )
            })
            .collect(),
        minted_errors: read_groups
            .iter()
            .map(|&group| {
                (
                    group,
                    MintedReadErrors::of_observation(-7.0 * 1_000_000.0, 1_000_000),
                )
            })
            .collect(),
        rates: BTreeMap::new(),
        inbreeding: Some(Estimate {
            value: InbreedingF::try_new(0.1).expect("a coefficient"),
            provenance: Provenance::FittedHere,
            observations: 1_000_000,
        }),
        runs_model: None,
        site_noise: None,
        site_noise_off_the_ladder: false,
        error_rate_on_a_ladder_end: BTreeSet::new(),
        coupled_fit: pop_var_caller::ng::parameter_estimation::fitting::FitTermination {
            iterations: 3,
            converged: true,
        },
    }
}

/// One sample's repeat-tract parameters at the stated shape: `strata` records for each of its
/// read groups, each weighing `genotypes` allele-length pairs.
fn repeat_tract_parameters_of(
    read_groups: &[ReadGroupId],
    strata: usize,
    genotypes: usize,
) -> SsrSampleParameters {
    let a_model = SlippageModel::try_new(0.01, 0.20, 0.065).expect("a slippage model");
    let a_start = pop_var_caller::ng::parameter_estimation::ssr::SlippageStart {
        from: a_model,
        reached: a_model,
        log_likelihood: -1.0,
    };
    let mut by_stratum = BTreeMap::new();
    for &group in read_groups {
        for index in 0..strata {
            // Repeat counts from four upward, at the six motif periods a catalogue holds, so the
            // keys are distinct and spread across periods the way a real run's are.
            let period = SsrPeriod::try_new(index % 6 + 1).expect("a legal period");
            let repeats = RepeatCount((index / 6 + 4) as u32);
            let stratum = Stratum::new(period, repeats);
            by_stratum.insert(
                StratumKey {
                    read_group: group,
                    stratum,
                    ploidy: diploid(),
                },
                StratumFit {
                    stratum,
                    slippage: Estimate {
                        value: a_model,
                        provenance: Provenance::FittedHere,
                        observations: 100_000,
                    },
                    substitution: Estimate {
                        value: ErrorRate::try_new(0.003).expect("a legal rate"),
                        provenance: Provenance::FittedHere,
                        observations: 1_000_000,
                    },
                    // Each genotype is a diploid pair of allele lengths, which is what the
                    // fit emits and what carries the `SmallVec` inside every record.
                    genotypes: (0..genotypes)
                        .map(|at| {
                            GenotypeFrequency::new(
                                [
                                    WholeRepeatOffset((at % 12) as i8 - 6),
                                    WholeRepeatOffset((at % 7) as i8 - 3),
                                ],
                                1.0 / genotypes as f64,
                            )
                        })
                        .collect(),
                    not_whole_repeat_share: 0.01,
                    unexplained_locus_share: 0.001,
                    starts_tried: SmallVec::from_slice(&[a_start]),
                    fitted_over: SmallVec::from_slice(&[stratum]),
                    shares_fitted_over: SmallVec::from_slice(&[stratum]),
                    slipped_reads: 9_000,
                },
            );
        }
    }
    SsrSampleParameters {
        by_stratum,
        summary: read_groups
            .iter()
            .map(|&group| (group, StratumFitSummary::default()))
            .collect(),
    }
}

/// Live bytes on the heap right now.
#[cfg(feature = "dhat-heap")]
fn live_bytes() -> u64 {
    dhat::HeapStats::get().curr_bytes as u64
}

#[cfg(not(feature = "dhat-heap"))]
fn live_bytes() -> u64 {
    0
}

fn main() {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    let number = |at: usize, fallback: usize| -> usize {
        arguments
            .get(at)
            .map_or(fallback, |value| value.parse().expect("a whole number"))
    };
    let strata = number(0, STRATA_A_LIBRARY);
    let libraries = number(1, LIBRARIES_A_SAMPLE);
    let genotypes = number(2, GENOTYPES_A_STRATUM);

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::builder().testing().build();
    #[cfg(not(feature = "dhat-heap"))]
    eprintln!("warning: built without --features dhat-heap, so every byte count below reads zero");

    println!(
        "shape: {libraries} librar{} a sample, {strata} strata a library, {genotypes} genotypes \
         a stratum",
        if libraries == 1 { "y" } else { "ies" }
    );
    println!(
        "{:>8}  {:>14}  {:>13}  {:>13}  {:>14}  {:>13}",
        "samples", "inputs (B)", "union (B)", "assembly (B)", "peak (B)", "peak/sample"
    );

    for cohort in COHORT_SIZES {
        let read_group_count = cohort * libraries;
        let read_groups: Vec<Vec<ReadGroupId>> = (0..cohort)
            .map(|sample| {
                (0..libraries)
                    .map(|library| ReadGroupId((sample * libraries + library) as u32))
                    .collect()
            })
            .collect();

        let before_inputs = live_bytes();
        let generic: Vec<GenericSampleParameters> = read_groups
            .iter()
            .map(|mine| generic_parameters_of(mine))
            .collect();
        let repeat_tract: Vec<SsrSampleParameters> = read_groups
            .iter()
            .map(|mine| repeat_tract_parameters_of(mine, strata, genotypes))
            .collect();
        let after_inputs = live_bytes();

        // The seam's own union, assembled the way `from_prepass` assembles it. `assemble` is
        // called rather than `from_prepass` because the latter needs the run's read-group table,
        // which is read from the alignment files' headers and cannot be built in memory — and
        // what is being measured is what is *held*, which is the same either way.
        let before_union = live_bytes();
        let mut error_rate = BTreeMap::new();
        let mut minted = BTreeMap::new();
        let mut substitution = BTreeMap::new();
        let mut inbreeding = Vec::with_capacity(cohort);
        for (sample, of_sample) in generic.iter().enumerate() {
            error_rate.extend(
                of_sample
                    .error_rate
                    .iter()
                    .map(|(&group, rate)| (group, rate.clone())),
            );
            minted.extend(of_sample.minted_errors.iter().map(|(&g, &t)| (g, t)));
            inbreeding.push(
                of_sample
                    .inbreeding
                    .as_ref()
                    .expect("every sample here has one")
                    .value,
            );
            substitution.extend(repeat_tract[sample].substitution_rate_by_stratum());
        }

        let before_assembled = live_bytes();
        let assembled = RunParameters::assemble(
            &error_rate,
            &minted,
            &BTreeMap::new(),
            SequencingBatches::all_together_over(read_group_count, cohort),
            inbreeding,
            SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
            StratumFits::over(&[], BTreeMap::new()),
            substitution,
            diploid(),
        );
        let after_assembled = live_bytes();

        let inputs = after_inputs - before_inputs;
        let union = before_assembled - before_union;
        let assembly = after_assembled - before_assembled;
        let peak = after_assembled - before_inputs;
        println!(
            "{cohort:>8}  {inputs:>14}  {union:>13}  {assembly:>13}  {peak:>14}  {:>13}",
            peak / cohort as u64
        );

        drop(assembled);
        drop(repeat_tract);
        drop(generic);
    }
}
