//! **What the repeat-tract guard leaves out of a cohort's evidence, under the test it used to make
//! and under the one the spec states.**
//!
//! The guard leaves a repeat tract out of the slippage fit when, of its reads that differ from the
//! reference tract length, more than one in ten differ by a number of bases that is not a whole
//! number of repeat units (`parameter_prepass_joint_records.md` §3.3). The spec states that as a
//! test of the tract. Until 2026-10-14 the fit made it once for each sample's each read group and
//! left the tract out if any one failed; at about three reads a sample, one stray read in one
//! sample was enough. This tool reads a cohort's tract evidence once and counts both ways:
//!
//! - **tracts left out** and **tracts kept**;
//! - **repeat-length classes with enough tracts to fit** — tracts with reads at least the climb's
//!   refusal floor (8, kept as the bench's `COUNTING_FIT_FLOOR` so these counts compare with those
//!   printed before the climb was deleted; the penalised fit fits every class with reads);
//! - **(tract, sample) pairs with reads** among the kept tracts;
//! - for the tracts the spec's test brings back, how their share of non-whole reads is spread,
//!   against the same share at the tracts both tests keep.
//!
//! It then gathers the evidence as `estimate-parameters` does, which prints the run's own
//! "repeat-tract evidence: arranged … left out by the guard" line: that number must equal this
//! tool's count under the spec's test. **On a large cohort that step holds the arranged evidence
//! — 19.4 GiB on 2,169 tomato samples — so `SKIP_GATHER=1` stops before it**; the counts above
//! are complete without it.
//!
//! ```text
//! cargo run --release --example ng_str_guard_census -- \
//!     --reference ref.fa --catalog catalog.parquet --psp psps/ --output unused.toml
//! ```
//!
//! **The flags are `estimate-parameters`'s**, and the cohort is opened and checked as that
//! command opens it; `--output` is required by them and never written. Nothing is fitted: the
//! `FIT_CLASSES` option, which fitted the classes the pooled test brings back the most tracts to
//! with the climb, went with the climb on 2026-10-16 (`doc/devel/ng/arch/str_penalised_fit.md`
//! §4.1).

use clap::Parser;

use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::census::{GuardTally, SsrLocusState, Stratum};
use pop_var_caller::parameter_estimation::joint::ssr_fit::expected_counts::COUNTING_FIT_FLOOR;
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    StratumEvidence, gather_strata, strata_of_kept_loci,
};
use pop_var_caller::run::every_read_group_pooled;

#[derive(Parser)]
#[command(about = "Count the tracts the repeat-tract guard leaves out, per section and pooled")]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
}

/// What the guard sees at one tract.
#[derive(Default, Clone, Copy)]
struct TractGuard {
    /// The guard's two counts over every sample and read group.
    pooled: GuardTally,
    /// Whether any one sample's read group failed the threshold on its own — the old test.
    any_section_over: bool,
    /// Samples with a read crossing the whole tract.
    samples_with_reads: u32,
}

impl TractGuard {
    fn kept_before(&self) -> bool {
        !self.any_section_over
    }

    fn kept_now(&self) -> bool {
        !self.pooled.is_over_threshold()
    }
}

/// One side's totals: the old test or the spec's.
#[derive(Default)]
struct Totals {
    left_out: u64,
    kept: u64,
    pairs_with_reads: u64,
    classes_enough: u64,
}

fn totals(classes: &[(Stratum, Vec<TractGuard>)], keep: impl Fn(&TractGuard) -> bool) -> Totals {
    let mut out = Totals::default();
    for (_, tracts) in classes {
        let mut with_reads: usize = 0;
        for tract in tracts {
            if keep(tract) {
                out.kept += 1;
                out.pairs_with_reads += u64::from(tract.samples_with_reads);
                with_reads += usize::from(tract.samples_with_reads > 0);
            } else {
                out.left_out += 1;
            }
        }
        out.classes_enough += u64::from(with_reads >= COUNTING_FIT_FLOOR);
    }
    out
}

/// How a set of tracts' non-whole shares are spread, as counts in fixed bands, plus the pooled
/// share over all their reads.
fn spread(tracts: &[&TractGuard]) -> String {
    let edges = [0.0, 0.02, 0.05, 0.08, 0.10, 0.20, 0.50, 1.0];
    let mut counts = vec![0_u64; edges.len()];
    let mut no_differing = 0_u64;
    let mut all = GuardTally::default();
    for tract in tracts {
        all.add(tract.pooled);
        match tract.pooled.non_whole_share() {
            None => no_differing += 1,
            Some(0.0) => counts[0] += 1,
            Some(share) => {
                let band = edges
                    .windows(2)
                    .position(|pair| share > pair[0] && share <= pair[1])
                    .expect("a share is in 0 to 1");
                counts[band + 1] += 1;
            }
        }
    }
    let mut text = format!(
        "{} tracts; {} with no differing read; share 0: {}",
        tracts.len(),
        no_differing,
        counts[0]
    );
    for (band, pair) in edges.windows(2).enumerate() {
        text.push_str(&format!(
            "; ({:.2}, {:.2}]: {}",
            pair[0],
            pair[1],
            counts[band + 1]
        ));
    }
    text.push_str(&format!(
        "; over all their reads {} non-whole of {} differing ({:.4})",
        all.non_whole_reads,
        all.non_whole_reads + all.whole_repeat_differing_reads,
        all.non_whole_share().unwrap_or(0.0)
    ));
    text
}

fn main() {
    let Cli { args } = Cli::parse();
    let OpenedCensuses {
        mut evidence, plan, ..
    } = open_the_censuses(&args).unwrap_or_else(|error| panic!("{}", with_its_causes(&error)));
    let samples = evidence.sample_names().count();
    let contig_of = contig_id_of(&plan);
    let kept = strata_of_kept_loci(&plan.loci, &contig_of);
    let mut tracts_in: std::collections::BTreeMap<Stratum, usize> = Default::default();
    for stratum in &kept {
        *tracts_in.entry(*stratum).or_insert(0) += 1;
    }
    let band: Vec<Stratum> = tracts_in.keys().copied().collect();

    let classes: Vec<(Stratum, Vec<TractGuard>)> = evidence
        .with_strata(&band, |lent| {
            tracts_in
                .iter()
                .map(|(stratum, tracts)| {
                    let mut guards = vec![TractGuard::default(); *tracts];
                    for sections in lent {
                        let mut crossed = vec![false; *tracts];
                        for (_, section_stratum, records) in sections {
                            if section_stratum != stratum {
                                continue;
                            }
                            for (locus, guard) in guards.iter_mut().enumerate() {
                                let tally = records.section_guard_tally(locus);
                                guard.pooled.add(tally);
                                guard.any_section_over |= tally.is_over_threshold();
                                crossed[locus] |= records.state(locus) == SsrLocusState::Crossed;
                            }
                        }
                        for (guard, crossed) in guards.iter_mut().zip(crossed) {
                            guard.samples_with_reads += u32::from(crossed);
                        }
                    }
                    (*stratum, guards)
                })
                .collect()
        })
        .unwrap_or_else(|error| panic!("{error}"));

    let before = totals(&classes, TractGuard::kept_before);
    let now = totals(&classes, TractGuard::kept_now);
    println!(
        "{samples} samples, {} tracts in {} repeat-length classes",
        kept.len(),
        classes.len()
    );
    for (name, side) in [
        ("each sample's read group on its own (before)", &before),
        ("pooled over the tract (spec, now)", &now),
    ] {
        println!(
            "{name}: {} tracts left out by the guard, {} kept, {} classes with at least {} \
             tracts with reads, {} (tract, sample) pairs with reads",
            side.left_out,
            side.kept,
            side.classes_enough,
            COUNTING_FIT_FLOOR,
            side.pairs_with_reads
        );
    }

    let all: Vec<&TractGuard> = classes.iter().flat_map(|(_, t)| t.iter()).collect();
    let kept_both: Vec<&TractGuard> = all
        .iter()
        .copied()
        .filter(|t| t.kept_before() && t.kept_now())
        .collect();
    let returning: Vec<&TractGuard> = all
        .iter()
        .copied()
        .filter(|t| !t.kept_before() && t.kept_now())
        .collect();
    let dropped_now: Vec<&TractGuard> = all.iter().copied().filter(|t| !t.kept_now()).collect();
    let newly_dropped = all
        .iter()
        .filter(|t| t.kept_before() && !t.kept_now())
        .count();
    println!("kept by both tests: {}", spread(&kept_both));
    println!("brought back by the pooled test: {}", spread(&returning));
    println!("still left out: {}", spread(&dropped_now));
    println!("kept before and left out now: {newly_dropped} (expected 0)");

    // The run's own gathering, which prints the evidence line this tool's counts must match.
    if std::env::var("SKIP_GATHER").is_ok_and(|value| value == "1") {
        return;
    }
    let pooled = every_read_group_pooled(&evidence);
    let gathered: Vec<StratumEvidence> =
        gather_strata(&mut evidence, &kept, &pooled).unwrap_or_else(|error| panic!("{error}"));
    let gathered_left_out: u64 = gathered
        .iter()
        .map(|stratum| stratum.tracts_over_guard_threshold)
        .sum();
    assert_eq!(
        gathered_left_out, now.left_out,
        "the gathering and this tool disagree on the tracts the pooled test leaves out"
    );
}

fn with_its_causes(error: &dyn std::error::Error) -> String {
    let mut text = error.to_string();
    let mut source = error.source();
    while let Some(cause) = source {
        text.push_str(&format!("\n  because: {cause}"));
        source = cause.source();
    }
    text
}
