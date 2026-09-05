//! **The joint parameters fit's heap profile — both halves, no CRAM, no reference.**
//!
//! Build and run on the host (the review's target machine):
//!
//! ```text
//! cargo run --profile profiling --example dhat_ng_joint_fit --features dhat-heap -- all
//! cargo run --profile profiling --example dhat_ng_joint_fit --features dhat-heap -- generic
//! cargo run --profile profiling --example dhat_ng_joint_fit --features dhat-heap -- gather
//! cargo run --profile profiling --example dhat_ng_joint_fit --features dhat-heap -- tracts
//! ```
//!
//! **Four runs and not one, because dhat reports a single `t-gmax`** — the one
//! instant the process held the most — and only the sites live at that instant
//! get exact byte attribution. The three phases have three separate peak
//! questions, so in one process the largest would hide the other two. `all` still
//! answers the *churn* question (total blocks and bytes a site), which is a sum
//! over the run and does not care which phase peaked.
//!
//! **`--profile profiling` and not `--release`**, for the reason `Cargo.toml`
//! gives at that profile: release is `lto = "fat"` with `debug =
//! "line-tables-only"`, and dhat's backtraces then collapse into frames that name
//! no source line. The `profiling` profile is release speed with `lto = false`,
//! `codegen-units = 16` and full debug info, which is what makes a site
//! attributable to `src/ng/parameter_estimation/…rs:line`. The counts themselves
//! are the same either way.
//!
//! Each run writes `dhat-heap-<mode>.json` in the working directory. Open one at
//! <https://nnethercote.github.io/dh_view/dh_view.html>, or parse it offline —
//! the stacks are deep, so attribute a site to the first
//! `src/ng/parameter_estimation/…rs:line` frame, skipping this file's own
//! allocator hook and the alloc/core/BTreeMap internals.
//!
//! **One caveat when parsing offline, and it bites.** dhat splits one source
//! line into many program points when the stacks above it differ, which under
//! rayon they always do — 48 program points for one `Vec::with_capacity` in
//! `TractLikelihoods::of`. So a site's `tb`/`tbk` (total bytes and blocks) may be
//! summed across its program points, because totals add; its `mb` (that program
//! point's own maximum) may **not**, because two maxima need not happen at the
//! same instant. For live footprint use `gb` — bytes at `t-gmax` — and arrange
//! for the phase you care about to be the one that peaks, which is what the modes
//! are for.
//!
//! ## Why this example exists
//!
//! The performance review of 2026-08-15
//! (`doc/devel/reports/reviews/perf_ng-census-joint-fit_2026-08-15.md`, §3 item 5)
//! records that **no allocation profile of this module has ever been taken**, and
//! that the footprint arithmetic in its findings L2, L3 and L9 is computed from
//! `size_of` rather than measured. The same review closed the allocator question
//! on the *time* axis — every `libsystem_malloc` symbol in its profile summed to
//! 516 of 245,018 samples, 0.21% — so **nothing here is a wall-clock claim**. What
//! it measures is footprint: how many blocks and bytes each phase takes, and how
//! that grows with the two axes `CLAUDE.md` §0 commits the caller to, cohort size
//! and evidence size.
//!
//! ## What it reaches, and what it does not
//!
//! Both halves of the fit, from data built in memory:
//!
//! | mode | phase | what it runs | size here |
//! |---|---|---|---|
//! | `generic` | 1 | `fit_jointly` — the **ordinary-position half**, contamination included | 8 and 24 samples × 30,000 kept positions |
//! | `gather` | 2 | `gather_strata` — the census→fit re-encoding of the **repeat-tract** evidence | 8 samples × 6 strata × 2,000 and 4,000 tracts |
//! | `tracts` | 3 | `fit_strata` — the **repeat-tract half**, borrowing and dedup included | 8 samples × 4 strata × 25 and 50 tracts |
//! | `tracts` | 4 | `fit_stratum` with the climb switched off — what one fat stratum's read likelihoods weigh | 8 samples × 3,000 tracts |
//!
//! Every phase runs at two sizes, because a single number cannot separate a fixed
//! cost from a per-sample or per-tract one. Phase 3 is small on purpose: the
//! review measures the repeat-tract fit at about 0.045 s a tract at 8 samples on
//! four threads, so a stratum wide enough to be interesting takes longer than a
//! profiling run should — which is why phase 4 exists to measure the wide case's
//! *footprint* with `max_rounds = 0` and one starting point, buying the width by
//! giving up the climb.
//!
//! **Not reached:** the census *walk* (it needs reads), the census *file*
//! (`census_file.rs` — a file-backed run needs a written census, which needs the
//! walk), and `loci.rs`'s selection (it needs a reference and a repeat catalog).
//! Those three are the CRAM-gated part of the module and are left to a harness
//! that has one.
//!
//! ## Where the fixtures come from
//!
//! The drawn cohort of phase 1 and the drawn stratum of phase 3 are **duplicated
//! from the two positive controls inside the module's own `#[cfg(test)]` blocks**
//! — `draw_cohort_with_duplications` (`fit.rs`) and `draw_stratum`
//! (`ssr_fit.rs`) — rather than moved out of them, so this example changes no
//! line of `src/`. They are simplified: the draw only has to produce evidence of
//! the right *shape and size* for an allocation profile, not evidence the fit
//! recovers known parameters from. Nothing here asserts a fitted value, and
//! nothing here should be read as an accuracy check.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::collections::BTreeMap;
use std::sync::Arc;

use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, CohortCensusEvidence, DepthCap, DepthCode, DepthLadderDigest,
    GenericEvidence, NamedReadGroup, ObservedAllele, OffsetCounts, PackedDepthCodes, ReadCap,
    RecordingTerms, SampleCensusEvidence, Section, SectionKey, SelectionTermsDigest, SsrEvidence,
    Stratum, WalkedBits,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{JointFitConfig, fit_jointly};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::parameter_estimation::joint::ssr_fit::{
    self, SampleTractReads, Slippage, SsrFitConfig, StratumEvidence, TractReads,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::{Ploidy, ReadGroupId};

fn main() {
    // Which phases to run. **The mode exists for one reason: dhat reports a
    // single `t-gmax`**, the instant the whole process held the most, and only
    // the sites live *then* get exact byte attribution. Three phases of this
    // module have separate peak questions, so running them in one process would
    // let the largest hide the other two. One mode a question, one file each.
    let mode = std::env::args().nth(1).unwrap_or_else(|| "all".to_string());
    // `prepared` is `tracts` without its two churn phases: phase 4 alone, which is
    // the cheap one, so a change to the quadrature or the read likelihoods can be
    // checked for **bit-identity of the fitted numbers** in seconds rather than
    // after the millions of objective evaluations the climb makes.
    let (generic_half, gather, tract_churn, prepared) = match mode.as_str() {
        "generic" => (true, false, false, false),
        "gather" => (false, true, false, false),
        "tracts" => (false, false, true, true),
        "prepared" => (false, false, false, true),
        "all" => (true, true, true, true),
        other => {
            eprintln!("unknown mode {other:?}: expected all, generic, gather, tracts or prepared");
            std::process::exit(2);
        }
    };

    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::builder()
        .file_name(format!("dhat-heap-{mode}.json"))
        .build();

    #[cfg(not(feature = "dhat-heap"))]
    {
        let _ = (generic_half, gather, tract_churn, prepared);
        eprintln!(
            "This example measures allocations and needs the dhat allocator.\n\
             Re-run with: cargo run --profile profiling --example dhat_ng_joint_fit \
             --features dhat-heap -- {mode}"
        );
    }

    #[cfg(feature = "dhat-heap")]
    {
        println!(
            "{:<46} {:>9} {:>13} {:>12} {:>12}",
            "phase", "blocks", "bytes", "live", "peak"
        );
        println!("{}", "-".repeat(96));

        // ---- phase 1: the ordinary-position half, at two cohort sizes --------
        //
        // Two sizes and not one, because the four dense `samples × positions`
        // arrays the review's L9 names are the only thing here that should grow
        // with both — everything else the pass allocates is sized once.
        for samples in [8_usize, 24]
            .into_iter()
            .take(if generic_half { 2 } else { 0 })
        {
            let positions = 30_000;
            let (cohort, _) =
                measured(&format!("build resident cohort, {samples} samples"), || {
                    drawn_generic_cohort(samples, positions, 12.0, 0x9E37_79B9_7F4A_7C15)
                });
            let mut cohort = cohort;
            let config = JointFitConfig {
                ploidy: Ploidy::try_new(2).expect("two is a ploidy"),
                // Six rather than the default 200: a pass allocates the same
                // things as the one before it, so passes past the first say
                // nothing new about the shape and only cost time.
                max_passes: 6,
                edges: Arc::new(DepthBinEdges::for_census()),
                ..JointFitConfig::default()
            };
            let (fit, _) = measured(
                &format!("fit_jointly, {samples} samples x {positions} positions"),
                || fit_jointly(&mut cohort, &config).expect("a drawn cohort fits"),
            );
            // Kept so the fit is not optimised away, and printed so a reader can
            // see the run was a real one.
            println!(
                "    (fitted {} read groups, {} samples' contamination)",
                fit.noise.len(),
                fit.contamination.len()
            );
        }

        // ---- phase 2: the census → fit re-encoding of the tract evidence -----
        //
        // Two evidence sizes, so the per-(tract, sample) slope the review's L2
        // computes from `size_of` can be measured as a difference instead — and
        // then a third run in which **every sample has a read at every tract**,
        // because the structure is linear in (tract, sample-with-reads) pairs and
        // that is the memory bound's worst case rather than its typical one. It
        // goes last so it is the run's peak, which is the only instant dhat
        // attributes live bytes at.
        for (tracts_a_stratum, one_sample_in) in [(2_000_usize, 3_u32), (4_000, 3), (4_000, 1)]
            .into_iter()
            .take(if gather { 3 } else { 0 })
        {
            let samples = 8;
            let strata_count = 6;
            let (built, _) = measured(
                &format!(
                    "build resident SSR cohort, {tracts_a_stratum} tracts, 1 sample in \
                     {one_sample_in} reads"
                ),
                || {
                    drawn_ssr_cohort(
                        samples,
                        strata_count,
                        tracts_a_stratum,
                        one_sample_in,
                        0xD1B5_4A32_D192_ED03,
                    )
                },
            );
            let (mut cohort, strata) = built;
            let groups: BTreeMap<ReadGroupId, u32> = (0..samples)
                .map(|s| (ReadGroupId(s as u32), 0_u32))
                .collect();
            let (gathered, _) = measured(
                &format!(
                    "gather_strata, {samples} samples x {} tracts",
                    strata_count as usize * tracts_a_stratum
                ),
                || ssr_fit::gather_strata(&mut cohort, &strata, &groups).expect("resident"),
            );
            // **Rows, not tracts, is the unit the footprint is linear in**: one
            // row is one (tract, sample-with-reads), and a sample with no read
            // there contributes none.
            let rows: usize = gathered
                .iter()
                .flat_map(|stratum| stratum.tracts.iter())
                .map(|tract| tract.samples.len())
                .sum();
            println!(
                "    ({} strata, {} tracts with reads, {rows} rows, {} spanning reads)",
                gathered.len(),
                gathered
                    .iter()
                    .map(StratumEvidence::tracts_with_reads)
                    .sum::<usize>(),
                gathered
                    .iter()
                    .map(StratumEvidence::spanning_reads)
                    .sum::<u64>()
            );
            drop(gathered);
        }

        // ---- phase 3: the repeat-tract half ---------------------------------
        //
        // Four thin strata of one motif length. The borrowing floor is left at
        // its default, so every stratum borrows the other three, the pooled sets
        // are identical, and `fit_strata`'s dedup collapses four fits to one —
        // which is exactly the path the review's L4 and L5 are about.
        //
        // Two tract counts, so what the fit allocates can be told from what it
        // allocates *per tract* — the review's L3 is an arithmetic claim about
        // exactly that slope.
        let samples = 8;
        for tracts in [25_usize, 50]
            .into_iter()
            .take(if tract_churn { 2 } else { 0 })
        {
            let (strata, _) = measured(&format!("draw 4 strata x {tracts} tracts"), || {
                (0..4)
                    .map(|ring| {
                        draw_stratum(
                            Slippage {
                                level: 0.08,
                                shorter_share: 0.8,
                                fall_off: 0.25,
                            },
                            &spectrum_of(13),
                            3.0,
                            0.2,
                            tracts,
                            samples,
                            6,
                            4,
                            Stratum {
                                period: 2,
                                reference_repeats: 10 + ring,
                            },
                            0xA076_1D64_78BD_642F ^ ring,
                        )
                    })
                    .collect::<Vec<_>>()
            });
            let excess = vec![0.2_f64; samples];
            let config = SsrFitConfig::default();
            let (outcomes, _) = measured(
                &format!("fit_strata, 4 strata x {tracts} tracts x {samples} samples"),
                || ssr_fit::fit_strata(&strata, &excess, &config),
            );
            println!("    ({} outcomes)", outcomes.len());
        }

        // ---- phase 4: what one fat stratum's `Prepared` actually weighs ------
        //
        // **The climb is switched off here and that is the point.** `max_rounds =
        // 0` and one starting point leave `fit_stratum` building the `Scorer`,
        // refreshing it once and scoring once — which is every allocation the
        // review's L3 is an arithmetic claim about, and none of the millions of
        // objective evaluations that would otherwise make a stratum this wide take
        // minutes. What it costs is that this phase says nothing about churn; the
        // two phases above say that.
        if prepared {
            let tracts = 3_000;
            let (evidence, _) = measured(&format!("draw one stratum x {tracts} tracts"), || {
                draw_stratum(
                    Slippage {
                        level: 0.08,
                        shorter_share: 0.8,
                        fall_off: 0.25,
                    },
                    &spectrum_of(13),
                    3.0,
                    0.2,
                    tracts,
                    samples,
                    6,
                    4,
                    Stratum {
                        period: 2,
                        reference_repeats: 10,
                    },
                    0x2545_F491_4F6C_DD1D,
                )
            });
            let excess = vec![0.2_f64; samples];
            // **One start, then the three the estimator ships with, and nothing
            // else different.** `fit_pooled` constructs its `Scorer` *inside* the
            // starting-point loop, so if the review's L3 is right that the read
            // likelihoods are built once per start, the second run's blocks at
            // `TractLikelihoods::of` must be exactly three times the first's. That
            // is a countable prediction and this is what counts it.
            for starts in [1_usize, 3] {
                let config = SsrFitConfig {
                    max_rounds: 0,
                    starting_points: ssr_fit::StartingPoint::spanning_the_monomorphic_range()
                        .into_iter()
                        .take(starts)
                        .collect(),
                    ..SsrFitConfig::default()
                };
                let (fit, _) = measured(
                    &format!(
                        "fit_stratum, {tracts} tracts x {samples} samples, no climb, \
                         {starts} start(s)"
                    ),
                    || ssr_fit::fit_stratum(&evidence, &excess, &config),
                );
                // **Printed at full precision on purpose.** A change to the
                // quadrature or the read likelihoods that claims to be
                // value-preserving has to show these bytes unchanged, and `{:?}`
                // on an `f64` is the shortest string that round-trips — so a diff
                // of two runs is a bit-identity check.
                match fit {
                    Some(fit) => println!(
                        "    (ln L a tract {:?}, concentration {:?}, slippage {:?}, \
                         spectrum[6] {:?})",
                        fit.log_likelihood_a_tract,
                        fit.concentration,
                        fit.slippage[0],
                        fit.length_spectrum[6],
                    ),
                    None => println!("    (no fit)"),
                }
            }
        }

        let stats = dhat::HeapStats::get();
        println!();
        println!("whole run:");
        println!("  total blocks {}", stats.total_blocks);
        println!("  total bytes  {}", stats.total_bytes);
        println!("  peak blocks  {}", stats.max_blocks);
        println!("  peak bytes   {}", stats.max_bytes);
        println!(
            "  live at end  {} blocks, {} bytes",
            stats.curr_blocks, stats.curr_bytes
        );
    }
}

/// Run `body`, print what it allocated, and hand back its value.
///
/// **Blocks and bytes, not seconds.** Allocation counts are identical run to run
/// whatever else the machine is doing, which is the whole reason this vehicle is
/// a heap profile rather than a timer.
/// Four numbers a phase, and they answer different questions. **Blocks** and
/// **bytes** are what the phase asked the allocator for in total, churn included;
/// **live** is what it still held when it returned, which for a phase that builds
/// a structure is that structure's own footprint; **peak** is the whole run's
/// high-water mark so far, so a phase that does not raise it is a phase whose own
/// peak sat under an earlier one's.
#[cfg(feature = "dhat-heap")]
fn measured<T>(what: &str, body: impl FnOnce() -> T) -> (T, u64) {
    let before = dhat::HeapStats::get();
    let value = body();
    let after = dhat::HeapStats::get();
    let blocks = after.total_blocks - before.total_blocks;
    let bytes = after.total_bytes - before.total_bytes;
    let live = after.curr_bytes as i64 - before.curr_bytes as i64;
    println!(
        "{what:<46} {blocks:>9} {bytes:>13} {live:>12} {:>12}",
        after.max_bytes
    );
    (value, blocks)
}

// ---------------------------------------------------------------------
// The fixtures — duplicated from the module's own positive controls
// ---------------------------------------------------------------------

/// The deterministic stream every drawn number comes from.
///
/// The same xorshift the two `#[cfg(test)]` generators use, copied so this file
/// depends on nothing behind a test gate.
struct Draw(u64);

impl Draw {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    /// A Poisson count, by inversion. Small means, so the loop is short.
    fn poisson(&mut self, mean: f64) -> u32 {
        let limit = (-mean).exp();
        let mut product = self.uniform();
        let mut count = 0;
        while product > limit && count < 1_000 {
            product *= self.uniform();
            count += 1;
        }
        count
    }

    fn pick(&mut self, weights: &[f64]) -> usize {
        let mut u = self.uniform() * weights.iter().sum::<f64>();
        for (index, weight) in weights.iter().enumerate() {
            u -= weight;
            if u <= 0.0 {
                return index;
            }
        }
        weights.len() - 1
    }
}

/// The twelve recording values every sample in one cohort has to agree on.
fn recording_terms() -> RecordingTerms {
    RecordingTerms {
        selection: SelectionTermsDigest::of(&SelectionTerms {
            seed: 42,
            reference: ReferenceDigest([7; 16]),
            analysed_regions: RegionSetDigest([9; 16]),
            catalog_built_under: CatalogBuildSettings {
                criteria: StrRepeatCriteria::default(),
                scan: ScanParams::default(),
                tool_version: "0.1.0".to_string(),
            },
            ssr_criteria: StrRepeatCriteria::default(),
            generic_target: 2_000_000,
            ssr_cap: 1_000,
        }),
        kept_loci: CensusLociDigester::new().finish(),
        ssr_stratum_counts: Default::default(),
        read_cap: ReadCap(100),
        depth_ladder: DepthLadderDigest::of(&DepthBinEdges::for_census()),
        depth_cap: DepthCap::new(124),
    }
}

/// A resident cohort of ordinary-position evidence — one read group a sample.
///
/// Simplified from `fit.rs`'s `draw_cohort_with_duplications`: nine positions in
/// ten carry no variation, the rest segregate at a frequency drawn flat, and
/// reads disagree with the reference at a fixed rate. What matters for a heap
/// profile is that the sparse non-reference list has a realistic *length* — it
/// is the one part of a generic section that is not a fixed size.
fn drawn_generic_cohort(
    samples: usize,
    positions: usize,
    mean_depth: f64,
    seed: u64,
) -> CohortCensusEvidence {
    let mut draw = Draw(seed);
    let edges = DepthBinEdges::for_census();
    let mut depth: Vec<PackedDepthCodes> = (0..samples)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
    let rate = 0.01;

    for index in 0..positions {
        let segregating = draw.uniform() < 0.10;
        let frequency = if segregating { draw.uniform() } else { 0.0 };
        let allele = 1 + (draw.uniform() * 3.0) as usize % 3;
        for sample in 0..samples {
            let reads = draw.poisson(mean_depth);
            let copies = if segregating {
                draw.pick(&[
                    (1.0 - frequency) * (1.0 - frequency),
                    2.0 * frequency * (1.0 - frequency),
                    frequency * frequency,
                ])
            } else {
                0
            };
            let carried = copies as f64 / 2.0;
            let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
            let mut counts = [0_u32; 5];
            for _ in 0..reads {
                let u = draw.uniform();
                let code = if u < on_candidate {
                    allele
                } else if u < on_candidate + rate / 3.0 {
                    (allele % 3) + 1
                } else if u < on_candidate + 2.0 * rate / 3.0 {
                    ((allele + 1) % 3) + 1
                } else {
                    0
                };
                counts[code] += 1;
            }
            depth[sample].set(index, DepthCode::Binned(edges.bin_for(reads)));
            for (code, count) in counts.iter().enumerate() {
                if code == 0 || *count == 0 {
                    continue;
                }
                sparse[sample].push(AlleleObservation {
                    index: index as u32,
                    allele: match code {
                        1 => ObservedAllele::C,
                        2 => ObservedAllele::G,
                        _ => ObservedAllele::T,
                    },
                    reads: u8::try_from(*count).unwrap_or(u8::MAX),
                });
            }
        }
    }

    let terms = recording_terms();
    let records: Vec<SampleCensusEvidence> = (0..samples)
        .map(|s| {
            SampleCensusEvidence::resident(
                format!("s{s}"),
                terms.clone(),
                NamedReadGroup::drawn_for(&format!("s{s}"), [ReadGroupId(s as u32)]),
                BTreeMap::new(),
                BTreeMap::from([(
                    SectionKey::Generic(ReadGroupId(s as u32)),
                    Section::Generic(GenericEvidence::from_parts(
                        std::mem::replace(&mut depth[s], PackedDepthCodes::never_walked(0)),
                        std::mem::take(&mut sparse[s]),
                    )),
                )]),
            )
        })
        .collect();
    CohortCensusEvidence::new(records).expect("every sample recorded the same thing")
}

/// A resident cohort of repeat-tract evidence, and the flat stratum-per-tract
/// list `gather_strata` reads the section lengths against.
///
/// **`share_with_reads` is the knob that matters here**: a sample with no read at
/// a tract contributes no row to the gathered structure, so the fraction of
/// (tract, sample) pairs that carry a read is what the footprint is linear in.
/// One in `share_with_reads` carries one, so 3 means about a third — the review
/// quotes census.rs:694 saying close to half the reads reaching a tomato tract
/// never cross it.
fn drawn_ssr_cohort(
    samples: usize,
    strata_count: u64,
    tracts_a_stratum: usize,
    share_with_reads: u32,
    seed: u64,
) -> (CohortCensusEvidence, Vec<Stratum>) {
    let mut draw = Draw(seed);
    let strata: Vec<Stratum> = (0..strata_count)
        .map(|ring| Stratum {
            period: 2,
            reference_repeats: 10 + ring,
        })
        .collect();
    let terms = recording_terms();

    let records: Vec<SampleCensusEvidence> = (0..samples)
        .map(|s| {
            let sections = strata
                .iter()
                .map(|stratum| {
                    let mut offsets = vec![OffsetCounts::default(); tracts_a_stratum];
                    let mut walked = WalkedBits::none_of(tracts_a_stratum);
                    for (tract, counts) in offsets.iter_mut().enumerate() {
                        walked.set(tract);
                        if (draw.uniform() * f64::from(share_with_reads)) as u32 != 0 {
                            continue;
                        }
                        for _ in 0..3 {
                            counts.add((draw.uniform() * 3.0) as i32 - 1, 1);
                        }
                    }
                    (
                        // One read group a sample: a library belongs to one plant.
                        SectionKey::Ssr(ReadGroupId(s as u32), *stratum),
                        Section::Ssr(SsrEvidence::from_parts(
                            offsets,
                            0,
                            walked,
                            0,
                            Vec::new(),
                            Vec::new(),
                        )),
                    )
                })
                .collect::<BTreeMap<_, _>>();
            SampleCensusEvidence::resident(
                format!("s{s}"),
                terms.clone(),
                NamedReadGroup::drawn_for(&format!("s{s}"), [ReadGroupId(s as u32)]),
                BTreeMap::new(),
                sections,
            )
        })
        .collect();

    // One entry a kept tract, in genome order — which for a drawn selection is
    // stratum by stratum.
    let per_tract: Vec<Stratum> = strata
        .iter()
        .flat_map(|stratum| std::iter::repeat_n(*stratum, tracts_a_stratum))
        .collect();
    (
        CohortCensusEvidence::new(records).expect("every sample recorded the same thing"),
        per_tract,
    )
}

/// One stratum drawn at known slippage — duplicated from `ssr_fit.rs`'s
/// `draw_stratum`, with the stratum key made an argument so several rings of one
/// motif length can be drawn and `fit_strata`'s borrowing has neighbours to take.
#[allow(
    clippy::too_many_arguments,
    reason = "the drawn stratum's own parameters"
)]
fn draw_stratum(
    slippage: Slippage,
    spectrum: &[f64],
    concentration: f64,
    homozygote_excess: f64,
    tracts: usize,
    samples: usize,
    depth: u32,
    span: i32,
    stratum: Stratum,
    seed: u64,
) -> StratumEvidence {
    let classes = spectrum.len();
    let buckets = (2 * span + 1) as usize;
    let per_allele: Vec<Vec<f64>> = (0..classes)
        .map(|class| slippage.read_probabilities(class as i32 - span, span))
        .collect();
    let mut draw = Draw(seed);
    let mut drawn = Vec::with_capacity(tracts);
    for _ in 0..tracts {
        // The test generator draws the tract's allele frequencies from a
        // Dirichlet; a flat draw over the same classes gives evidence of the
        // same shape and size, which is all a heap profile reads.
        let frequencies: Vec<f64> = spectrum
            .iter()
            .map(|weight| weight * (0.5 + draw.uniform()) * concentration.max(1e-3))
            .collect();
        let mut reads = TractReads::default();
        for sample in 0..samples {
            let first = draw.pick(&frequencies);
            let second = if draw.uniform() < homozygote_excess {
                first
            } else {
                draw.pick(&frequencies)
            };
            let mut counts = vec![0_u32; buckets];
            for _ in 0..depth {
                let allele = if draw.uniform() < 0.5 { first } else { second };
                counts[draw.pick(&per_allele[allele])] += 1;
            }
            reads.samples.push(SampleTractReads {
                sample: sample as u32,
                by_group: vec![(0, counts)],
            });
        }
        drawn.push(reads);
    }
    StratumEvidence {
        stratum,
        tracts: drawn,
        read_span: span,
        groups: 1,
        tracts_over_guard_threshold: 0,
        reads_reaching_not_crossing: 0,
        guard_reads: 0,
        // A drawn stratum has no sequence behind it, so it has no substitution rate. Zero bases
        // compared is what `substitution_rate()` returns `None` for, which is the honest answer
        // here rather than a fitted zero.
        bases_compared: 0,
        mismatching_bases: 0,
    }
}

/// Most chromosomes at the reference length, one repeat either side carrying most
/// of the rest — `ssr_fit.rs`'s `spectrum_of`.
fn spectrum_of(classes: usize) -> Vec<f64> {
    let middle = classes / 2;
    let mut spectrum: Vec<f64> = (0..classes)
        .map(|class| 0.55_f64.powi((class as i32 - middle as i32).abs()))
        .collect();
    let total: f64 = spectrum.iter().sum();
    for weight in &mut spectrum {
        *weight /= total;
    }
    spectrum
}
