//! **What the multi-library cell key costs** — measures the proposed key against exact
//! per-library scoring, over the range the caller will actually meet.
//!
//! The question. When two libraries with different error rates cover one site, the site's
//! likelihood is one sum over genotypes with *both* libraries' reads inside it, each
//! weighed against its own rate. A cell key of total depth and total alternative count has
//! thrown that away. The proposal is to key on **which libraries the alternative reads came
//! from** (for `k <= K`) and to drop the per-library **depth** split, approximating each
//! library's depth at a site by the sample's average share.
//!
//! What this measures: the error that approximation puts into the fitted parameters, in the
//! units the design already uses to decide such things — **rungs of the error-rate ladder**,
//! a quarter-Phred each, about 6% in probability.
//!
//! Why generating from the model is legitimate here, when the reviews warn against it: we
//! are not asking whether the model is right. We are asking what one *approximation* to it
//! costs against computing the same model exactly. Both sides use the same model by
//! construction, so a known truth is what makes the comparison possible at all.
//!
//! Only the depth split is approximated. Depth binning is a separate approximation and is
//! deliberately not mixed in — total depth is exact on both sides.
//!
//! ```text
//! cargo run --release --example ng_multilib_key_sweep
//! ```

use std::collections::BTreeMap;

/// SplitMix64 — dependency-free and reproducible, as `src/ssr/cohort/sim.rs` uses.
struct SplitMix64(u64);

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn unit(&mut self) -> f64 {
        // 53 bits of mantissa, in [0, 1).
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn binomial(&mut self, n: u32, p: f64) -> u32 {
        // Direct simulation: n is small here (depths up to a few hundred) and this keeps
        // the draw exact rather than approximate.
        (0..n).filter(|_| self.unit() < p).count() as u32
    }
    fn poisson(&mut self, mean: f64) -> u32 {
        // Knuth. Means here are <= 300, so the product form is fine.
        let limit = (-mean).exp();
        let mut k = 0u32;
        let mut prod = self.unit();
        while prod > limit {
            k += 1;
            prod *= self.unit();
            if k > 10_000 {
                break;
            }
        }
        k
    }
}

/// One genotype's chance that a single read shows a non-reference base, at error rate
/// `eps` and `j` alternative copies out of `ploidy` (spec `parameter_prepass.md` §3).
fn p_alt(j: u32, ploidy: u32, eps: f64) -> f64 {
    let frac = f64::from(j) / f64::from(ploidy);
    frac * (1.0 - eps / 3.0) + (1.0 - frac) * eps
}

#[derive(Clone)]
struct Scenario {
    name: &'static str,
    /// One error rate per library.
    eps: Vec<f64>,
    /// Each library's share of the reads. Sums to 1.
    weights: Vec<f64>,
    mean_depth: f64,
    pi_het: f64,
    pi_hom_alt: f64,
}

impl Scenario {
    fn genotype_freqs(&self) -> [f64; 3] {
        [
            1.0 - self.pi_het - self.pi_hom_alt,
            self.pi_het,
            self.pi_hom_alt,
        ]
    }
}

/// One site, as the exact key would hold it: per-library depth and alternative count.
type ExactCell = Vec<(u32, u32)>;

/// Draw `sites` sites and tally them by their exact per-library breakdown.
fn simulate(scn: &Scenario, sites: u64, seed: u64) -> BTreeMap<ExactCell, u64> {
    let mut rng = SplitMix64::new(seed);
    let freqs = scn.genotype_freqs();
    let libs = scn.eps.len();
    let mut table: BTreeMap<ExactCell, u64> = BTreeMap::new();

    for _ in 0..sites {
        // Genotype.
        let u = rng.unit();
        let j = if u < freqs[0] {
            0
        } else if u < freqs[0] + freqs[1] {
            1
        } else {
            2
        };

        let depth = rng.poisson(scn.mean_depth);
        if depth == 0 {
            continue;
        }

        // Split the depth between libraries: a chain of binomials is an exact multinomial.
        let mut remaining = depth;
        let mut remaining_weight = 1.0;
        let mut cell: ExactCell = Vec::with_capacity(libs);
        for g in 0..libs {
            let n_g = if g + 1 == libs {
                remaining
            } else {
                let p = (scn.weights[g] / remaining_weight).clamp(0.0, 1.0);
                let drawn = rng.binomial(remaining, p);
                remaining -= drawn;
                remaining_weight -= scn.weights[g];
                drawn
            };
            let k_g = rng.binomial(n_g, p_alt(j, 2, scn.eps[g]));
            cell.push((n_g, k_g));
        }
        *table.entry(cell).or_insert(0) += 1;
    }
    table
}

/// Log-likelihood of one cell under **exact** per-library scoring.
fn ln_exact(cell: &ExactCell, eps: &[f64], freqs: &[f64; 3]) -> f64 {
    let mut total = 0.0;
    for (j, &pi) in freqs.iter().enumerate() {
        if pi <= 0.0 {
            continue;
        }
        let mut ln_term = pi.ln();
        for (g, &(n_g, k_g)) in cell.iter().enumerate() {
            let p = p_alt(j as u32, 2, eps[g]);
            ln_term += f64::from(k_g) * p.ln() + f64::from(n_g - k_g) * (1.0 - p).ln();
        }
        total = ln_sum_exp(total, ln_term, j == 0);
    }
    total
}

/// Log-likelihood of the same cell under the **proposed** key: the alternative reads keep
/// their library attribution (when `k <= max_attributed`), the depths do not — each
/// library's depth is taken as its average share of the total.
fn ln_proposed(
    cell: &ExactCell,
    eps: &[f64],
    weights: &[f64],
    freqs: &[f64; 3],
    max_attributed: u32,
    exact_below_depth: u32,
) -> f64 {
    let depth: u32 = cell.iter().map(|&(n, _)| n).sum();
    let alt: u32 = cell.iter().map(|&(_, k)| k).sum();

    // Shallow sites keep their whole breakdown — depths as well as attribution — because
    // the average-share approximation is worst exactly where depth is small: crediting two
    // libraries with 1.5 reads each when the real split is 3/0 or 2/1.
    if depth <= exact_below_depth {
        return ln_exact(cell, eps, freqs);
    }

    // What the key retains about the alternative reads.
    let attributed: Vec<f64> = if alt <= max_attributed {
        cell.iter().map(|&(_, k)| f64::from(k)).collect()
    } else {
        weights.iter().map(|w| w * f64::from(alt)).collect()
    };

    let mut total = 0.0;
    for (j, &pi) in freqs.iter().enumerate() {
        if pi <= 0.0 {
            continue;
        }
        let mut ln_term = pi.ln();
        for g in 0..eps.len() {
            let p = p_alt(j as u32, 2, eps[g]);
            let n_g = weights[g] * f64::from(depth);
            let k_g = attributed[g];
            ln_term += k_g * p.ln() + (n_g - k_g).max(0.0) * (1.0 - p).ln();
        }
        total = ln_sum_exp(total, ln_term, j == 0);
    }
    total
}

fn ln_sum_exp(acc: f64, term: f64, first: bool) -> f64 {
    if first {
        return term;
    }
    let (hi, lo) = if acc > term { (acc, term) } else { (term, acc) };
    hi + (lo - hi).exp().ln_1p()
}

/// The ladder the design scans: Phred 10 to 50 in quarter-Phred steps.
fn ladder() -> Vec<f64> {
    (0..=160)
        .map(|i| 10f64.powf(-(10.0 + 0.25 * f64::from(i)) / 10.0))
        .collect()
}

/// Profile one library's error rate: hold everything else at truth, scan that rate, and
/// return the index of the rung each scoring puts its maximum on.
fn profile_eps(
    table: &BTreeMap<ExactCell, u64>,
    scn: &Scenario,
    target: usize,
    max_attributed: u32,
    exact_below: u32,
) -> (usize, usize) {
    let rungs = ladder();
    let freqs = scn.genotype_freqs();
    let mut best_exact = (f64::NEG_INFINITY, 0usize);
    let mut best_prop = (f64::NEG_INFINITY, 0usize);

    for (r, &candidate) in rungs.iter().enumerate() {
        let mut eps = scn.eps.clone();
        eps[target] = candidate;
        let mut score_exact = 0.0;
        let mut score_prop = 0.0;
        for (cell, &count) in table {
            let c = count as f64;
            score_exact += c * ln_exact(cell, &eps, &freqs);
            score_prop += c * ln_proposed(
                cell,
                &eps,
                &scn.weights,
                &freqs,
                max_attributed,
                exact_below,
            );
        }
        if score_exact > best_exact.0 {
            best_exact = (score_exact, r);
        }
        if score_prop > best_prop.0 {
            best_prop = (score_prop, r);
        }
    }
    (best_exact.1, best_prop.1)
}

/// Profile the heterozygosity the same way, on a log grid, returning the relative gap
/// between where the two scorings put their maxima.
fn profile_het(
    table: &BTreeMap<ExactCell, u64>,
    scn: &Scenario,
    max_attributed: u32,
    exact_below: u32,
) -> f64 {
    // Centred on truth and fine enough that a real shift is not confused with a grid step:
    // 121 points spanning 0.2x to 3x the true value.
    let grid: Vec<f64> = (0..=120)
        .map(|i| scn.pi_het * (0.2 + 2.8 * f64::from(i) / 120.0))
        .filter(|&h| h + scn.pi_hom_alt < 1.0)
        .collect();
    let mut best_exact = (f64::NEG_INFINITY, scn.pi_het);
    let mut best_prop = (f64::NEG_INFINITY, scn.pi_het);

    for &h in &grid {
        let freqs = [1.0 - h - scn.pi_hom_alt, h, scn.pi_hom_alt];
        let mut score_exact = 0.0;
        let mut score_prop = 0.0;
        for (cell, &count) in table {
            let c = count as f64;
            score_exact += c * ln_exact(cell, &scn.eps, &freqs);
            score_prop += c * ln_proposed(
                cell,
                &scn.eps,
                &scn.weights,
                &freqs,
                max_attributed,
                exact_below,
            );
        }
        if score_exact > best_exact.0 {
            best_exact = (score_exact, h);
        }
        if score_prop > best_prop.0 {
            best_prop = (score_prop, h);
        }
    }
    (best_prop.1 - best_exact.1) / best_exact.1
}

fn scenarios() -> Vec<Scenario> {
    let mut out = Vec::new();
    // The axes that matter: error rate across the ladder, the ratio between libraries,
    // depth, heterozygosity, split skew, and library count.
    for &(name, eps_lo, ratio) in &[
        ("eps=1e-3 ratio=1", 1e-3, 1.0),
        ("eps=1e-3 ratio=4", 1e-3, 4.0),
        ("eps=1e-3 ratio=10", 1e-3, 10.0),
        ("eps=1e-2 ratio=4", 1e-2, 4.0),
        ("eps=1e-4 ratio=4", 1e-4, 4.0),
    ] {
        for &depth in &[3.0f64, 20.0, 60.0, 300.0] {
            for &(split_name, w0) in &[("even", 0.5), ("skew90", 0.9)] {
                out.push(Scenario {
                    name: Box::leak(
                        format!("{name} depth={depth:.0} split={split_name}").into_boxed_str(),
                    ),
                    eps: vec![eps_lo, eps_lo * ratio],
                    weights: vec![w0, 1.0 - w0],
                    mean_depth: depth,
                    pi_het: 1e-3,
                    pi_hom_alt: 6e-4,
                });
            }
        }
    }
    // Low heterozygosity, and four libraries, at the depth where it should matter most.
    // Heterozygosity needs its own scenarios: at one per kilobase, two million sites hold
    // only ~2,000 heterozygotes and the profile is Monte Carlo noise. Ten per kilobase
    // makes the column resolvable without changing what is being tested.
    for &depth in &[3.0f64, 20.0, 300.0] {
        out.push(Scenario {
            name: Box::leak(format!("het-resolvable ratio=4 depth={depth:.0}").into_boxed_str()),
            eps: vec![1e-3, 4e-3],
            weights: vec![0.5, 0.5],
            mean_depth: depth,
            pi_het: 1e-2,
            pi_hom_alt: 6e-3,
        });
    }
    out.push(Scenario {
        name: "four libraries depth=20",
        eps: vec![5e-4, 1e-3, 2e-3, 4e-3],
        weights: vec![0.25, 0.25, 0.25, 0.25],
        mean_depth: 20.0,
        pi_het: 1e-3,
        pi_hom_alt: 6e-4,
    });
    out
}

/// How many distinct cells the proposed key would create, which is what it costs. Counted
/// from the data rather than from the combinatorics, so it is the *observed* key count and
/// not the possible one. Depth is left unbinned here, so this is an upper bound — the real
/// design bins it, which folds the pooled arm down further.
fn distinct_keys(
    table: &BTreeMap<ExactCell, u64>,
    max_attributed: u32,
    exact_below_depth: u32,
) -> usize {
    let mut keys: std::collections::BTreeSet<(u32, u32, Vec<u32>, bool)> =
        std::collections::BTreeSet::new();
    for cell in table.keys() {
        let depth: u32 = cell.iter().map(|&(n, _)| n).sum();
        let alt: u32 = cell.iter().map(|&(_, k)| k).sum();
        if depth <= exact_below_depth {
            // The whole breakdown is the key.
            let flat: Vec<u32> = cell.iter().flat_map(|&(n, k)| [n, k]).collect();
            keys.insert((depth, alt, flat, true));
        } else {
            let attribution: Vec<u32> = if alt <= max_attributed {
                cell.iter().map(|&(_, k)| k).collect()
            } else {
                Vec::new()
            };
            keys.insert((depth, alt, attribution, false));
        }
    }
    keys.len()
}

/// How much evidence a scenario actually carries, so a reader can tell a real shift from
/// Monte Carlo noise. The heterozygosity profile needs heterozygous sites; the error-rate
/// profile needs alternative reads at homozygous-reference sites.
fn evidence(table: &BTreeMap<ExactCell, u64>) -> (u64, u64) {
    let mut sites = 0;
    let mut alt_reads = 0;
    for (cell, &count) in table {
        sites += count;
        let k: u32 = cell.iter().map(|&(_, k)| k).sum();
        alt_reads += count * u64::from(k);
    }
    (sites, alt_reads)
}

fn main() {
    const SITES: u64 = 2_000_000;
    // (K = alternative reads that keep their attribution, D = depth below which the whole
    // breakdown is kept). K=4,D=0 is the previous proposal; D>0 is what this run measures.
    let variants: [(u32, u32); 6] = [(0, 0), (4, 0), (4, 4), (4, 8), (4, 12), (4, 20)];

    println!("# Multi-library cell key: what the approximation costs");
    println!("#");
    println!("# Sites per scenario: {SITES}. Depth exact on both sides; only the per-library");
    println!("# depth split is approximated. Shift is in ladder rungs (quarter-Phred, ~6% each);");
    println!(
        "# a shift of 0 means the approximation put the fit on the same rung as exact scoring."
    );
    println!("#");
    println!("# K = how many alternative reads keep their library attribution. K=0 is the");
    println!("# pooled key we have today — the baseline the proposal has to beat.");
    println!();
    println!(
        "{:<40} {:>2} {:>3} {:>9} {:>9} {:>11} {:>7}",
        "scenario", "K", "D", "eps0", "eps1", "het rel.err", "keys"
    );

    for scn in scenarios() {
        let table = simulate(&scn, SITES, 0xC0FFEE_u64 ^ (scn.name.len() as u64));
        let (sites, alt_reads) = evidence(&table);
        let expected_hets = (sites as f64 * scn.pi_het) as u64;
        println!(
            "# {} — {} sites, {} alternative reads, ~{} heterozygous sites",
            scn.name, sites, alt_reads, expected_hets
        );
        for (k, d) in variants {
            let (e0_exact, e0_prop) = profile_eps(&table, &scn, 0, k, d);
            let (e1_exact, e1_prop) = profile_eps(&table, &scn, 1, k, d);
            let het = profile_het(&table, &scn, k, d);
            let keys = distinct_keys(&table, k, d);
            println!(
                "{:<40} {:>2} {:>3} {:>9} {:>9} {:>10.2}% {:>7}",
                scn.name,
                k,
                d,
                e0_prop as i64 - e0_exact as i64,
                e1_prop as i64 - e1_exact as i64,
                het * 100.0,
                keys
            );
        }
        println!();
    }
}
