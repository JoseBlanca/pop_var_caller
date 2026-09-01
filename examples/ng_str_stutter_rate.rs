//! **How often does a read at an STR tract slip, per library and per repeat count?** — the
//! estimate the copy-floor question actually needs.
//!
//! `doc/devel/ng/spec/parameter_prepass_ssr.md` §4 fits four numbers per (read group × stratum), a
//! stratum being one (motif period, reference repeat count). The first of them — **how often a read
//! shows a length other than its allele's** — is the quantity behind "at how many repeats does
//! stutter start to matter", and it is per library because stutter is a property of the chemistry.
//!
//! ## What this measures that the library survey could not
//!
//! [`ng_str_stutter_by_library.rs`](ng_str_stutter_by_library.rs) counts reads that differ from the
//! **reference** tract length. That number is slippage *plus* reads that correctly measured an
//! allele the sample genuinely carries, and it cannot separate them: it pools loci, and the
//! separation lives inside a locus. **A real allele shifts all of a locus's reads; slippage shifts
//! some of them.**
//!
//! So this program keeps each locus's reads together, and sums over the genotype rather than
//! choosing one — scoring every possible pair of allele lengths and weighting each by how common
//! that pair turns out to be, with those frequencies fitted alongside. Pool the loci first and the
//! fitted rate moves **333-fold depending only on where the search starts** (spec §4.1).
//!
//! ## How it is checked, which is the part that matters
//!
//! A survey ran for two days on this archive and could not answer the question it was built for.
//! The instrument was never asked to reproduce an answer somebody already knew. So this program
//! **runs its checks before it will report anything**, through the same component builder, kernel
//! and optimiser it uses on real reads:
//!
//! 1. **Three algebraic gates, before any fit.** The bucket probabilities sum to one for every
//!    genotype; none of them is outside `[0, 1]`; and with the slippage level at zero every read
//!    lands on its locus's own alleles. Each rejects a broken scoring rule in one line and none
//!    needs data. The gates are also run against a rule known to be *wrong* — scoring a saturating
//!    end bucket at its edge — so they are proven to bite rather than merely to pass.
//! 2. **Exact recovery of a known truth.** Every possible locus shape, weighted by its exact
//!    probability under a chosen truth, then fitted. There is no sampling noise in this, so the
//!    answer is decided rather than estimated and the bias must be **zero**.
//! 3. **The answer is converged, not merely stopped.** Refit with ten times the climb passes; if
//!    the answer moves, what the previous check reported was this program's convergence rather
//!    than a property of the estimator. It caught exactly that on its first run.
//! 4. **The score at the truth.** A correctly specified model cannot be beaten at the truth, so a
//!    fitted point scoring *above* it is a defect here rather than a finding. This is the check
//!    that a spread across starting points cannot make: a deterministic search returns the same
//!    point from every start wherever the objective is flat, which is how a fall-off 0.02 out once
//!    passed four agreeing starts. It caught a second defect — a warm-started climb that made the
//!    objective depend on the order candidates were tried in.
//! 5. **A silent world.** Loci generated with no slippage at all must fit to no slippage.
//!
//! **All of that takes about ten seconds**, so every run does it and refuses to walk on a failure.
//! A number from an unchecked instrument is what this exercise exists to stop producing.
//!
//! **`--drawn` adds a sixth**, and it is opt-in because it costs minutes rather than seconds:
//! recovery from tens of thousands of *drawn* loci at the depth this archive has, which carries
//! sampling error and so answers a different question — not "is the estimator biased" but "is
//! there enough data here to see". Worth running when the model or the fit changes; not worth
//! paying before every walk. An earlier version bundled it into `--self-check` and left the reader
//! in front of a silent terminal for twenty minutes.
//!
//! ```text
//! ng_str_stutter_rate --self-check [--drawn]
//! ng_str_stutter_rate [--contigs a,b] [--regions r.bed] [--min-loci 100]
//!     [--max-repeats 40] <reference.fa> <sample.cram> [more...]
//! ```
//!
//! Output: one row per (library, period, reference repeat count) with the fitted slippage rate, the
//! direction split, the fall-off, the spread across starting points, and the evidence behind it.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::locus_generation::ssr::{SsrGenerator, SsrGeneratorConfig};
use pop_var_caller::ng::locus_generation::{
    LocusGenerator, LocusKind, ReadWitness, SampleLocusObservations,
};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::read::input::read_groups::build_read_groups;
use pop_var_caller::ng::read::input::reference::OpenReference;
use pop_var_caller::ng::ref_seq::WindowedRefSeq;
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::region_typing::segment_criteria::SsrSegment;
use pop_var_caller::ng::region_typing::{GenomeRegions, RegionKind, TypedRegionConfig};
use pop_var_caller::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use pop_var_caller::ng::types::GenomeRegion;
use pop_var_caller::ng::types::{Bp, ContigId, ReadGroupId};
use pop_var_caller::regions::ContigBounds;

#[path = "shared/stutter_table.rs"]
mod stutter_table;

use stutter_table::*;

#[path = "shared/catalog_regions.rs"]
mod catalog_regions;

/// Offsets are recorded over `±this`, the end buckets absorbing everything beyond. Measured to
/// matter far less than it looks: with the ends scored by their marginal, a range of ±1 against
/// alleles reaching ±3 still returns the rate to within 0.05% (spec §4.1).
const OFFSET_HALF_RANGE: i32 = 4;

/// How far from the reference length the fit may place an allele — **the width that is
/// load-bearing**, and a threshold to clear rather than a number to tune (spec §8.1).
const ALLELE_OFFSET_LIMIT: i32 = 6;

/// Reads entered from one locus. A deeper locus is entered from a subsample down to this, seeded
/// from the locus's position so a region-sharded walk and a single-threaded one keep the same reads.
const MAX_LOCUS_READS: u32 = 12;

/// Above this share of the reads that differ from the reference, this noise model does not describe
/// the stratum and a rate fitted from it is mostly mis-modelled ordinary indel (spec §5).
const GUARD_SHARE_LIMIT: f64 = 0.10;

/// Expectation-maximization passes the inner climb over the genotype frequencies is allowed.
///
/// **A knob because a climb that stops short tilts the outer search, and that is indistinguishable
/// from bias from the outside.** The self-check runs at exactly this value and then again at ten
/// times it; if the answer moves, this is too low and the check says so.
///
/// **200 was too low and the check caught it, but the cap was the symptom rather than the cause.**
/// The exact control recovered a 2.0000% truth as 1.9817% at 200 passes and 1.9977% at 2,000, so
/// most of the apparent bias was this number rather than the estimator — the failure the research
/// note records for the mode-centred arm (§6.3.1), arriving from a different direction.
///
/// What made the climb slow was the control's own key. Expectation-maximization converges at a rate
/// set by how much the genotypes overlap, and the three-bucket key chosen to keep the enumeration
/// cheap cannot tell fifteen genotypes apart. Widening it to five buckets converges in far fewer
/// passes and recovers the truth exactly, so the cap now has margin rather than being the limit.
///
/// **It costs little where it is not needed**, because the climb stops as soon as the frequencies
/// settle; a high cap binds only where convergence is genuinely slow.
const CLIMB_ITERATIONS: usize = 2_000;

fn keying() -> Keying {
    Keying {
        half_range: OFFSET_HALF_RANGE,
        edges: EdgeScoring::Marginal,
    }
}

// ---------------------------------------------------------------------------
// Turning one locus into one entry
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Stratum {
    period: u8,
    repeats: u32,
}

/// The stratum a **typed segment** belongs to, read off the reference tract before any read is
/// touched.
///
/// This is what lets a full stratum's segments be skipped without paying for them. The expensive
/// half of the walk is fetching a locus's reads and aligning them, and that happens inside
/// `next_locus`; the stratum is a property of the reference, so the decision can be taken before
/// the generator is asked for anything.
fn stratum_of_segment(segment: &SsrSegment) -> Option<Stratum> {
    let period = segment.period();
    let length = segment.tract_len();
    if period == 0 || !length.is_multiple_of(period as u64) {
        return None;
    }
    Some(Stratum {
        period: period as u8,
        repeats: (length / period as u64) as u32,
    })
}

/// The stratum a locus belongs to, from the **reference** tract alone — a pure function of the
/// reference, so every sample stratifies identically and a cohort can be compared.
fn stratum_of(locus: &SampleLocusObservations) -> Option<Stratum> {
    let LocusKind::Ssr(detail) = &locus.kind else {
        return None;
    };
    let period = detail.motif.period();
    if period == 0 || !locus.reference_bases.len().is_multiple_of(period) {
        return None;
    }
    Some(Stratum {
        period: period as u8,
        repeats: (locus.reference_bases.len() / period) as u32,
    })
}

/// A small deterministic hash, so a locus's subsample is a function of where it is and nothing
/// else — the property that makes a sharded walk and a single-threaded one agree.
fn seed_at(contig: u32, start: u64) -> u64 {
    let mut x = (u64::from(contig) << 40) ^ start.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    x ^= x >> 33;
    x = x.wrapping_mul(0xFF51_AFD7_ED55_8CCD);
    x ^= x >> 33;
    x
}

/// Reduce one locus to one shape per read group.
///
/// **Complete witnesses only.** A partial witness saw part of the tract, so its length is a lower
/// bound; scoring it as a length would read as a read that lost repeats, which is a direct bias in
/// the direction split — the parameter §3 exists to protect.
fn shapes_of(locus: &SampleLocusObservations) -> BTreeMap<ReadGroupId, LocusShape> {
    let buckets = keying().buckets();
    let middle = buckets / 2;
    let reference_len = locus.reference_bases.len() as i64;
    let Some(stratum) = stratum_of(locus) else {
        return BTreeMap::new();
    };
    let period = i64::from(stratum.period);

    let mut by_group: BTreeMap<ReadGroupId, LocusShape> = BTreeMap::new();
    for observation in locus.observations.iter() {
        if observation.read_witness != ReadWitness::Complete {
            continue;
        }
        let entry = by_group
            .entry(observation.read_group)
            .or_insert_with(|| LocusShape {
                counts: vec![0; buckets],
                not_whole_repeat: 0,
            });
        let difference = observation.bases.len() as i64 - reference_len;
        let reads = observation.num_obs;
        if difference != 0 && !difference.rem_euclid(period).eq(&0) {
            entry.not_whole_repeat = entry.not_whole_repeat.saturating_add(reads as u16);
            continue;
        }
        let offset = (difference / period) as i32;
        let bucket = offset.clamp(-OFFSET_HALF_RANGE, OFFSET_HALF_RANGE) + OFFSET_HALF_RANGE;
        let _ = middle;
        entry.counts[bucket as usize] = entry.counts[bucket as usize].saturating_add(reads as u16);
    }

    // Cap each group's depth by thinning uniformly, seeded from the locus position. A uniform
    // subsample is exact rather than approximate: it leaves the bucket counts distributed exactly
    // as they would be at the lower depth. What it costs is precision.
    let mut seed = seed_at(locus.region.contig.0, locus.region.start.0);
    for shape in by_group.values_mut() {
        let depth = shape.scored_depth() + u32::from(shape.not_whole_repeat);
        if depth <= MAX_LOCUS_READS {
            continue;
        }
        let mut kept = vec![0u16; shape.counts.len()];
        let mut kept_guard = 0u16;
        let mut pool: Vec<usize> = Vec::with_capacity(depth as usize);
        for (bucket, &count) in shape.counts.iter().enumerate() {
            for _ in 0..count {
                pool.push(bucket);
            }
        }
        for _ in 0..shape.not_whole_repeat {
            pool.push(usize::MAX);
        }
        // Partial Fisher-Yates against the seeded stream.
        for index in 0..MAX_LOCUS_READS as usize {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let pick = index + (seed >> 33) as usize % (pool.len() - index);
            pool.swap(index, pick);
            if pool[index] == usize::MAX {
                kept_guard += 1;
            } else {
                kept[pool[index]] += 1;
            }
        }
        shape.counts = kept;
        shape.not_whole_repeat = kept_guard;
    }
    by_group.retain(|_, shape| shape.scored_depth() > 0);
    by_group
}

// ---------------------------------------------------------------------------
// The checks
// ---------------------------------------------------------------------------

/// Every way `total` reads split across `buckets` buckets — the enumeration the exact control needs
/// and the real path deliberately avoids.
fn for_each_split(total: u32, buckets: usize, visit: &mut impl FnMut(&[u16])) {
    let mut counts = vec![0u16; buckets];
    fn recurse(counts: &mut Vec<u16>, slot: usize, left: u32, visit: &mut impl FnMut(&[u16])) {
        if slot + 1 == counts.len() {
            counts[slot] = left as u16;
            visit(counts);
            counts[slot] = 0;
            return;
        }
        for take in 0..=left {
            counts[slot] = take as u16;
            recurse(counts, slot + 1, left - take, visit);
        }
        counts[slot] = 0;
    }
    recurse(&mut counts, 0, total, visit);
}

/// A truth to recover: the slippage, the allele spectrum the loci are drawn from, and the depth.
struct World {
    slip: Slip,
    alleles: Vec<i32>,
    allele_probs: Vec<f64>,
    mean_depth: f64,
}

impl World {
    fn genotype_probs(&self) -> Vec<f64> {
        genotype_pairs(self.alleles.len())
            .into_iter()
            .map(|(i, j)| {
                let (pi, pj) = (self.allele_probs[i], self.allele_probs[j]);
                if i == j { pi * pi } else { 2.0 * pi * pj }
            })
            .collect()
    }

    /// Depth per locus: Poisson, conditioned on at least one read and truncated.
    fn depth_probs(&self, max_depth: u32) -> Vec<f64> {
        let mut probs = vec![0.0; (max_depth + 1) as usize];
        for (n, slot) in probs.iter_mut().enumerate().skip(1) {
            let n = n as u32;
            *slot =
                (-self.mean_depth + f64::from(n) * self.mean_depth.ln() - ln_factorial(n)).exp();
        }
        let total: f64 = probs.iter().sum();
        for p in &mut probs {
            *p /= total;
        }
        probs
    }
}

/// A table whose entries carry **exact probabilities** rather than integer locus counts.
///
/// The same shapes and the same component builder the real path uses; only the weights differ, and
/// that difference is what removes the sampling noise so bias is decided rather than estimated.
fn exact_table(world: &World, max_depth: u32, key: &Keying) -> (StratumTable, Vec<f64>) {
    let genotype_probs = world.genotype_probs();
    let genotypes = genotype_pairs(world.alleles.len());
    let depth_probs = world.depth_probs(max_depth);
    let mut weight_by_shape: BTreeMap<LocusShape, f64> = BTreeMap::new();

    for (genotype, &g_prob) in genotypes.iter().zip(&genotype_probs) {
        if g_prob == 0.0 {
            continue;
        }
        let bucket_probs = genotype_bucket_probs(&world.slip, &world.alleles, *genotype, key);
        for (depth, &d_prob) in depth_probs.iter().enumerate() {
            if d_prob == 0.0 {
                continue;
            }
            let depth = depth as u32;
            for_each_split(depth, key.buckets(), &mut |counts| {
                let mut ln = ln_factorial(depth);
                for (&count, &p) in counts.iter().zip(&bucket_probs) {
                    ln -= ln_factorial(u32::from(count));
                    if count > 0 {
                        ln += f64::from(count) * p.max(1e-300).ln();
                    }
                }
                let shape = LocusShape {
                    counts: counts.to_vec(),
                    not_whole_repeat: 0,
                };
                *weight_by_shape.entry(shape).or_insert(0.0) += g_prob * d_prob * ln.exp();
            });
        }
    }

    // Rebuild as a StratumTable so the real path's own accessors are exercised, then hand the exact
    // weights alongside — the integer counter cannot hold a probability.
    let mut table = StratumTable::default();
    let mut weights = Vec::with_capacity(weight_by_shape.len());
    for (shape, weight) in &weight_by_shape {
        table.add_locus(shape.clone());
        weights.push(*weight);
    }
    (table, weights)
}

/// Draw `loci` loci from the world and tally their shapes — the realistic case, with sampling noise.
fn drawn_table(world: &World, loci: u32, max_depth: u32, key: &Keying, seed: u64) -> StratumTable {
    let genotypes = genotype_pairs(world.alleles.len());
    let genotype_probs = world.genotype_probs();
    let depth_probs = world.depth_probs(max_depth);
    let mut state = seed | 1;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / ((1u64 << 53) as f64)
    };
    let pick = |u: f64, probs: &[f64]| {
        let mut acc = 0.0;
        for (i, &p) in probs.iter().enumerate() {
            acc += p;
            if u <= acc {
                return i;
            }
        }
        probs.len() - 1
    };

    let mut table = StratumTable::default();
    for _ in 0..loci {
        let genotype = genotypes[pick(next(), &genotype_probs)];
        let depth = pick(next(), &depth_probs) as u32;
        if depth == 0 {
            continue;
        }
        let bucket_probs = genotype_bucket_probs(&world.slip, &world.alleles, genotype, key);
        let mut counts = vec![0u16; key.buckets()];
        for _ in 0..depth {
            counts[pick(next(), &bucket_probs)] += 1;
        }
        table.add_locus(LocusShape {
            counts,
            not_whole_repeat: 0,
        });
    }
    table
}

/// The truth every check is run against: tomato's dinucleotides at six or more repeats — 2 reads in
/// 100 slip, 17 in 100 of those gain rather than lose, 9 in 100 take a second step (spec §3, §5).
fn measured_slip() -> Slip {
    Slip {
        level: 0.02,
        up_share: 0.17,
        falloff: 0.09,
    }
}

struct CheckOutcome {
    name: String,
    detail: String,
    passed: bool,
}

/// The checks, printed **as each one finishes** rather than collected and printed at the end.
///
/// A buffered report is fine for a check that takes a second and useless for one that takes
/// minutes: the reader sits in front of a silent terminal with no way to tell a slow check from a
/// hung one, and a run they interrupt shows nothing at all — including the results that had already
/// passed. That happened, so the collector prints on push.
#[derive(Default)]
struct Checks(Vec<CheckOutcome>);

impl Checks {
    fn push(&mut self, outcome: CheckOutcome) {
        eprintln!(
            "  [{}] {}\n        {}",
            if outcome.passed { "ok" } else { "FAILED" },
            outcome.name,
            outcome.detail
        );
        self.0.push(outcome);
    }
    fn iter(&self) -> std::slice::Iter<'_, CheckOutcome> {
        self.0.iter()
    }
    fn len(&self) -> usize {
        self.0.len()
    }
}

/// Run the checks. `include_drawn` adds the two that draw tens of thousands of loci.
///
/// **The split is about cost, not about importance.** Everything below the drawn pair runs on a
/// couple of hundred enumerated entries and takes under a second, so a real run pays for it every
/// time. The drawn checks build tables the size of a real stratum and take minutes — worth running
/// whenever the model or the fit changes, and not worth paying before every walk. `--self-check`
/// runs everything.
fn run_self_check(report: &mut String, include_drawn: bool) -> bool {
    let mut outcomes = Checks::default();
    eprintln!("# self-check\n");
    let alleles = vec![-2, -1, 0, 1, 2];
    let truth = measured_slip();

    // --- 1. The gates, which need no data and no fit ---------------------------------------
    let good = Keying {
        half_range: 2,
        edges: EdgeScoring::Marginal,
    };
    let sums = gate_sums_to_one(&truth, &alleles, &good);
    outcomes.push(CheckOutcome {
        name: "gate: bucket probabilities sum to one".into(),
        detail: format!("worst deviation {sums:.2e}"),
        passed: sums < 1e-12,
    });
    let range = gate_probabilities_in_range(&truth, &alleles, &good);
    outcomes.push(CheckOutcome {
        name: "gate: every bucket probability is in [0,1]".into(),
        detail: format!("worst violation {range:.2e}"),
        passed: range == 0.0,
    });
    let silent = gate_silent_kernel(&alleles, &good);
    outcomes.push(CheckOutcome {
        name: "gate: at a zero level every read sits on its own alleles".into(),
        detail: format!("worst stray mass {silent:.2e}"),
        passed: silent < 1e-12,
    });
    // The gate must reject a rule known to be wrong, or it is not a gate.
    let bad = Keying {
        half_range: 2,
        edges: EdgeScoring::PlugAtEdge,
    };
    let bad_sum = gate_sums_to_one(&truth, &alleles, &bad);
    outcomes.push(CheckOutcome {
        name: "gate bites: the edge plug-in is rejected by it".into(),
        detail: format!("plug-in deviates {bad_sum:.4} from one"),
        passed: bad_sum > 1e-6,
    });

    // --- 2. Exact recovery: no sampling noise, so the bias is decided -----------------------
    // A narrow recorded range keeps the enumerated shape space affordable, which §6.4 measured
    // costs nothing: the marginal rule is exactly unbiased at every range tried.
    // **Five buckets rather than three, and the width is chosen against the climb rather than
    // against the enumeration.** The shape space grows steeply with the bucket count, so the
    // cheapest control is the narrowest — but expectation-maximization converges at a rate set by
    // how much the genotypes overlap, and three buckets cannot tell fifteen genotypes apart, so the
    // climb crawls and the control ends up measuring its own convergence instead of the estimator.
    // Five buckets at a lower depth keeps the enumeration small (461 shapes) and the genotypes
    // separable.
    let narrow = Keying {
        half_range: 2,
        edges: EdgeScoring::Marginal,
    };
    let world = World {
        slip: truth,
        alleles: alleles.clone(),
        allele_probs: vec![0.05, 0.15, 0.60, 0.15, 0.05],
        mean_depth: 4.0,
    };
    let (exact, weights) = exact_table(&world, 6, &narrow);
    let exact_fit = fit_exact(&exact, &weights, &alleles, &narrow, CLIMB_ITERATIONS);
    let level_bias = 100.0 * (exact_fit.slip.level - truth.level) / truth.level;
    outcomes.push(CheckOutcome {
        name: "exact: the slippage rate is recovered without bias".into(),
        detail: format!(
            "{:.4}% against a truth of {:.4}% — bias {level_bias:+.3}%, {} entries, starts agree to {:.3}x",
            100.0 * exact_fit.slip.level,
            100.0 * truth.level,
            exact.entries(),
            exact_fit.level_spread
        ),
        passed: level_bias.abs() < 1.0 && exact_fit.level_spread < 1.06,
    });

    // **Is the residual the estimator's, or the inner climb stopping short?** The two are
    // indistinguishable from the outside, and telling them apart is the check that caught a
    // retracted finding: a fall-off 0.02 above its truth, from four starting points that agreed to
    // four decimal places, turned out to be the climb rather than the accumulator (spec §4.2). So
    // refit with ten times the passes and see whether the answer moves. If it does, the number
    // above is this program's convergence and not a property of the estimator, and
    // `CLIMB_ITERATIONS` is too low.
    let patient = fit_exact(&exact, &weights, &alleles, &narrow, CLIMB_ITERATIONS * 10);
    let moved = 100.0 * (patient.slip.level - exact_fit.slip.level).abs() / exact_fit.slip.level;
    outcomes.push(CheckOutcome {
        name: "exact: the answer is converged, not merely stopped".into(),
        detail: format!(
            "{} passes gives {:.4}%, {} gives {:.4}% — moves {moved:.3}%",
            CLIMB_ITERATIONS,
            100.0 * exact_fit.slip.level,
            CLIMB_ITERATIONS * 10,
            100.0 * patient.slip.level,
        ),
        passed: moved < 0.5,
    });
    let share_bias = (exact_fit.slip.up_share - truth.up_share).abs();
    outcomes.push(CheckOutcome {
        name: "exact: the direction split is recovered without bias".into(),
        detail: format!(
            "{:.4} against a truth of {:.4} — out by {share_bias:.4}",
            exact_fit.slip.up_share, truth.up_share
        ),
        passed: share_bias < 0.01,
    });

    // --- 3. The score at the truth, which a correct model cannot be beaten at ---------------
    let at_truth = score_exact(&exact, &weights, &truth, &alleles, &narrow);
    let at_fit = score_exact(&exact, &weights, &exact_fit.slip, &alleles, &narrow);
    outcomes.push(CheckOutcome {
        name: "exact: the fit does not beat the truth".into(),
        detail: format!("fitted − truth = {:+.3e} nats", at_fit - at_truth),
        passed: at_fit - at_truth < 1e-6,
    });

    // --- 4. A silent world: no slippage in, no slippage out ---------------------------------
    let quiet = World {
        slip: Slip {
            level: 0.0,
            up_share: 0.5,
            falloff: 0.1,
        },
        alleles: alleles.clone(),
        allele_probs: vec![0.05, 0.15, 0.60, 0.15, 0.05],
        mean_depth: 6.0,
    };
    let (quiet_table, quiet_weights) = exact_table(&quiet, 6, &narrow);
    let quiet_fit = fit_exact(
        &quiet_table,
        &quiet_weights,
        &alleles,
        &narrow,
        CLIMB_ITERATIONS,
    );
    outcomes.push(CheckOutcome {
        name: "control: a world with no slippage fits to none".into(),
        detail: format!("{:.5}% recovered", 100.0 * quiet_fit.slip.level),
        passed: quiet_fit.slip.level < 1e-3,
    });

    // --- 5. Recovery from drawn loci, at the depth this archive has -------------------------
    if include_drawn {
        let wide = keying();
        let realistic = World {
            slip: truth,
            alleles: (-4..=4).collect(),
            allele_probs: {
                let mut p: Vec<f64> = (-4i32..=4).map(|d| 0.6f64.powi(d.abs())).collect();
                let total: f64 = p.iter().sum();
                for v in &mut p {
                    *v /= total;
                }
                p
            },
            mean_depth: 8.0,
        };
        for &loci in &[20_000u32, 60_000] {
            let drawn = drawn_table(&realistic, loci, MAX_LOCUS_READS, &wide, 0x5EED_1234);
            let fit = fit_stratum(
                &drawn,
                12,
                ALLELE_OFFSET_LIMIT,
                &wide,
                CLIMB_ITERATIONS,
                SearchPrecision::fast(),
            );
            let bias = 100.0 * (fit.slip.level - truth.level) / truth.level;
            outcomes.push(CheckOutcome {
                name: format!("drawn: {loci} loci at 8 reads each recover the rate"),
                detail: format!(
                    "{:.4}% against {:.4}% — {bias:+.1}%, {} entries, starts agree to {:.2}x",
                    100.0 * fit.slip.level,
                    100.0 * truth.level,
                    drawn.entries(),
                    fit.level_spread
                ),
                // Sampling noise is real here, so this is a tolerance and not a zero. At 20,000 loci
                // and 8 reads the count of slipped reads is a few thousand, so a few percent is the
                // resolution; 25% is loose enough not to fire on noise and tight enough that the
                // failures this program exists to prevent — a factor of two or a sign flip — cannot pass.
                passed: bias.abs() < 25.0,
            });
        }
    }

    let passed = outcomes.iter().all(|o| o.passed);
    let _ = writeln!(
        report,
        "  {} of {} checks passed",
        outcomes.iter().filter(|o| o.passed).count(),
        outcomes.len()
    );
    passed
}

/// Fit a table whose weights are exact probabilities rather than locus counts.
fn fit_exact(
    table: &StratumTable,
    weights: &[f64],
    alleles: &[i32],
    key: &Keying,
    climb_iterations: usize,
) -> MultistartFit {
    let shapes = table.shapes();
    // Cold every time — see `fit_stratum`: a warm start makes the objective depend on the order
    // candidates were tried in, and the outer search cannot compare two such numbers.
    let objective = |slip: &Slip| {
        let component = component_matrix(&shapes, slip, alleles, key);
        let freqs = climb_genotype_frequencies(&component, weights, climb_iterations, None);
        let score = score_at(&component, weights, &freqs);
        (score, freqs)
    };
    fit_from_starts(objective, &starting_points(0.02), SearchPrecision::fine())
}

fn score_exact(
    table: &StratumTable,
    weights: &[f64],
    slip: &Slip,
    alleles: &[i32],
    key: &Keying,
) -> f64 {
    let shapes = table.shapes();
    let component = component_matrix(&shapes, slip, alleles, key);
    let freqs = climb_genotype_frequencies(&component, weights, CLIMB_ITERATIONS, None);
    score_at(&component, weights, &freqs)
}

// ---------------------------------------------------------------------------
// The walk
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn run(
    fasta: &Path,
    alignments: &[PathBuf],
    contig_filter: &[String],
    regions_bed: Option<&Path>,
    min_loci: u64,
    max_loci_per_stratum: u64,
    max_repeats: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        fasta.to_path_buf(),
        ReferenceCheck::VerifyAgainstIndex,
    )?;
    let contigs: ContigList = info.contig_list();
    // **The typed regions come from the catalog beside the reference**, checked against what
    // the pass just reported. No catalog, no run: the error names the command that writes one.
    let catalog = RepeatCatalog::open_beside_reference(fasta, &info)?;
    let reference = OpenReference::new(info);
    let read_groups = build_read_groups(alignments)?;
    let samples: Vec<SampleReads> = read_groups
        .read_groups_per_sample()
        .iter()
        .map(|entry| {
            SampleReads::open(
                entry,
                &read_groups,
                &reference,
                ReadFilterConfig::default(),
                true,
            )
        })
        .collect::<Result<_, _>>()?;
    eprintln!(
        "  {} sample(s), {} read group(s)",
        samples.len(),
        read_groups.iter().count()
    );

    let walk_config = TypedRegionConfig::default();
    let criteria = StrRepeatCriteria::from(&walk_config);
    let bundle_threshold = Bp(walk_config.criteria.bundle_threshold);
    let shared_reference = Arc::new(WindowedRefSeq::new(fasta.to_path_buf(), contigs.clone()));
    let mut generators = samples
        .iter()
        .map(|_| {
            SsrGenerator::with_default_aligner(
                Arc::clone(&shared_reference),
                {
                    let reference = Arc::clone(&shared_reference);
                    move || Arc::clone(&reference)
                },
                SsrGeneratorConfig::default(),
                bundle_threshold,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;

    let bed_spans = regions_bed
        .map(|bed| {
            let bounds: Vec<ContigBounds<'_>> = contigs
                .entries
                .iter()
                .map(|e| ContigBounds {
                    name: &e.name,
                    length: e.length as u32,
                })
                .collect();
            GenomeRegions::from_bed_path(bed, &bounds)
        })
        .transpose()?;

    let wanted_contig = |contig: ContigId| {
        contig_filter.is_empty()
            || contigs
                .entries
                .get(contig.0 as usize)
                .is_some_and(|e| contig_filter.iter().any(|n| n == &e.name))
    };

    // **Built on demand rather than once, because the region walk is made twice**: a
    // reference-only pass to count each stratum's segments, then the real pass that fetches reads.
    // A `TypedRegionIterator` is consumed by walking it, so the second pass needs its own.
    // The stretches to ask the catalog for, labelled. Built once and reused by both passes:
    // a reader borrows the catalog, so what is collected is the spans, not the readers.
    let batches: Vec<(String, Vec<GenomeRegion>)> = {
        let mut walks = Vec::new();
        match bed_spans.clone() {
            Some(spans) => {
                walks.push((
                    "BED".to_string(),
                    spans.iter().collect::<Vec<GenomeRegion>>(),
                ));
            }
            None => {
                for (index, entry) in contigs.entries.iter().enumerate() {
                    if !wanted_contig(ContigId(index as u32)) {
                        continue;
                    }
                    walks.push((
                        entry.name.clone(),
                        vec![catalog_regions::whole_contig(
                            ContigId(index as u32),
                            entry.length,
                        )],
                    ));
                }
            }
        }
        walks
    };

    // **A reference-only pre-pass, so the sample can be spread instead of truncated.**
    //
    // Region typing reads the reference and nothing else, so counting how many segments each
    // stratum holds costs no read fetching and no alignment — the two things that make the real
    // walk expensive. Knowing the totals up front is what lets the cap keep an even spread across
    // the whole region rather than the first N it meets, which measurement showed is biased by a
    // quarter (see the skip below).
    let mut segments_by_stratum: BTreeMap<Stratum, u64> = BTreeMap::new();
    if max_loci_per_stratum < u64::MAX {
        for (label, spans) in &batches {
            eprintln!("  counting strata over {label} (reference only)");
            let mut walk = catalog.genome_segments(&criteria, ReadScope::Regions(spans))?;
            for region in walk.by_ref() {
                let region = region?;
                let RegionKind::SsrSegment(segment) = &region.kind else {
                    continue;
                };
                if !wanted_contig(region.region.contig) {
                    continue;
                }
                if let Some(stratum) = stratum_of_segment(segment)
                    && stratum.repeats <= max_repeats
                {
                    *segments_by_stratum.entry(stratum).or_insert(0) += 1;
                }
            }
        }
    }

    let mut tables: BTreeMap<(ReadGroupId, Stratum), StratumTable> = BTreeMap::new();
    // How many segments of each stratum have been seen, and how many were skipped by the cap.
    // Counted per **segment** and not per read group, because the skip covers every sample at
    // once — they all walk the same reference.
    let mut walked_by_stratum: BTreeMap<Stratum, u64> = BTreeMap::new();
    let mut skipped_by_stratum: BTreeMap<Stratum, u64> = BTreeMap::new();
    for (label, spans) in &batches {
        eprintln!("  walking {label}");
        let mut walk = catalog.genome_segments(&criteria, ReadScope::Regions(spans))?;
        for region in walk.by_ref() {
            let region = region?;
            let RegionKind::SsrSegment(segment) = &region.kind else {
                continue;
            };
            if !wanted_contig(region.region.contig) {
                continue;
            }

            // **Decide before fetching anything, and spread the sample rather than truncating it.**
            //
            // A stratum holding enough loci learns nothing from another one, and the reads it would
            // cost are the walk's whole expense — fetching and aligning happen inside `next_locus`
            // below, while the stratum is a property of the reference tract and needs no reads to
            // read off. So a cap can be applied here, for free.
            //
            // **But "the first N" is not "N of them", and the difference is a quarter of the
            // answer.** Measured over 20 Mb of tomato chromosome 1: capping mononucleotides at 6
            // repeats to the first 2,953 of 28,125 returns 0.0541% against the full 0.0803% —
            // **32.6% low** — with every uncapped stratum in the same run agreeing to the digit.
            // Tracts near the start of a chromosome stutter measurably less, so a prefix is a
            // biased sample of one and the bias is one-directional.
            //
            // The counts come from a reference-only pre-pass, so each stratum takes an **even
            // spread across the whole region**: segment `i` of `n` is kept when it falls on the
            // next of `cap` evenly spaced slots. Deterministic, exactly `cap` kept, and no reliance
            // on where a stratum's segments happen to sit.
            if let Some(stratum) = stratum_of_segment(segment) {
                if stratum.repeats > max_repeats {
                    continue;
                }
                let total = *segments_by_stratum.get(&stratum).unwrap_or(&0);
                let seen = walked_by_stratum.entry(stratum).or_insert(0);
                let index = *seen;
                *seen += 1;
                if total > max_loci_per_stratum {
                    let before = index * max_loci_per_stratum / total;
                    let after = (index + 1) * max_loci_per_stratum / total;
                    if before == after {
                        *skipped_by_stratum.entry(stratum).or_insert(0) += 1;
                        continue;
                    }
                }
            }

            for (sample, generator) in samples.iter().zip(generators.iter_mut()) {
                generator.begin_segment(region.region);
                while let Some(locus) = generator.next_locus(segment, sample)? {
                    let Some(stratum) = stratum_of(&locus) else {
                        continue;
                    };
                    if stratum.repeats > max_repeats {
                        continue;
                    }
                    for (group, shape) in shapes_of(&locus) {
                        tables.entry((group, stratum)).or_default().add_locus(shape);
                    }
                }
            }
        }
    }

    let stdout = std::io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());
    for (id, group) in read_groups.iter() {
        writeln!(
            out,
            "#rg\t{}\t{}\t{}\t{}\t{}",
            id.get(),
            group.id,
            group.sample,
            group.library.value,
            group.file.display(),
        )?;
    }
    writeln!(
        out,
        "read_group\tperiod\trepeats\tloci\treads\tentries\tslip_rate\tgain_share\t\
         step_decay\tstart_spread\tguard_share\theterozygosity\tidentified"
    )?;

    let key = keying();
    let total = tables.len();
    for (index, ((group, stratum), table)) in tables.iter().enumerate() {
        if table.loci() < min_loci {
            continue;
        }
        eprintln!(
            "  fitting {}/{}: read group {}, period {}, {} repeats — {} loci, {} entries",
            index + 1,
            total,
            group.get(),
            stratum.period,
            stratum.repeats,
            table.loci(),
            table.entries()
        );
        let fit = fit_stratum(
            table,
            stratum.repeats,
            ALLELE_OFFSET_LIMIT,
            &key,
            CLIMB_ITERATIONS,
            SearchPrecision::fast(),
        );
        // **Two ways a row is not an answer, and both are reported rather than hidden.** A search
        // whose starting points disagreed found where it stopped, not what the data says; and a
        // stratum whose differing reads are mostly not whole-repeat changes is one this noise model
        // does not describe, so its rate is mis-modelled indel however many loci stood behind it.
        let identified = if fit.level_spread > 1.06 {
            "starts-disagree"
        } else if fit.not_whole_repeat_share > GUARD_SHARE_LIMIT {
            "model-does-not-fit"
        } else {
            "yes"
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.4}\t{:.4}\t{:.3}\t{:.4}\t{:.4}\t{}",
            group.get(),
            stratum.period,
            stratum.repeats,
            fit.loci,
            fit.reads,
            table.entries(),
            fit.slip.level,
            fit.slip.up_share,
            fit.slip.falloff,
            fit.level_spread,
            fit.not_whole_repeat_share,
            fit.heterozygosity,
            identified,
        )?;
        out.flush()?;
    }

    if let Some(handle) = verify {
        handle.join()?;
    }
    Ok(())
}

fn main() -> ExitCode {
    let mut positional: Vec<String> = Vec::new();
    let mut contig_filter: Vec<String> = Vec::new();
    let mut regions_bed: Option<PathBuf> = None;
    let mut min_loci: u64 = 2_000;
    let mut max_repeats: u32 = 30;
    let mut max_loci_per_stratum: u64 = 20_000;
    let mut only_check = false;
    let mut include_drawn = false;
    let mut rest = std::env::args().skip(1);
    while let Some(arg) = rest.next() {
        match arg.as_str() {
            "--self-check" => only_check = true,
            "--drawn" => include_drawn = true,
            "--regions" => regions_bed = rest.next().map(PathBuf::from),
            "--contigs" => {
                contig_filter = rest
                    .next()
                    .map(|l| l.split(',').map(|s| s.trim().to_string()).collect())
                    .unwrap_or_default()
            }
            "--min-loci" => min_loci = rest.next().and_then(|v| v.parse().ok()).unwrap_or(min_loci),
            "--max-loci-per-stratum" => {
                max_loci_per_stratum = rest
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(max_loci_per_stratum)
            }
            "--max-repeats" => {
                max_repeats = rest
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(max_repeats)
            }
            _ => positional.push(arg),
        }
    }

    // **The checks run before anything else, always.** A number from an unchecked instrument is
    // what this exercise exists to stop producing, so a failure here stops the run rather than
    // printing a warning above thousands of rows nobody will scroll back through.
    let mut report = String::new();
    let checks_passed = run_self_check(&mut report, include_drawn);
    eprint!("{report}");
    if !checks_passed {
        eprintln!("\nerror: the self-check failed — refusing to measure anything");
        return ExitCode::FAILURE;
    }
    if only_check {
        return ExitCode::SUCCESS;
    }

    if positional.len() < 2 {
        eprintln!(
            "usage: ng_str_stutter_rate [--contigs a,b] [--regions r.bed] [--min-loci N] \
             [--max-repeats N] <reference.fa> <sample.bam|cram> [sample ...]\n\
             fits the per-read slippage rate for each (library, motif period, reference repeat \
             count), summing over the genotype so a real allele is not charged to stutter."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&positional[0]);
    let alignments: Vec<PathBuf> = positional[1..].iter().map(PathBuf::from).collect();
    match run(
        &fasta,
        &alignments,
        &contig_filter,
        regions_bed.as_deref(),
        min_loci,
        max_loci_per_stratum,
        max_repeats,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
