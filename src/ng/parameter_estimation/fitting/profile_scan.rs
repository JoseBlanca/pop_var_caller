//! The scan: step through a ladder of noise parameters, climb to the best genotype
//! frequencies at each rung, and keep the rung that scores highest.
//!
//! In the standard vocabulary this is a **profile likelihood** over the noise
//! parameters — the frequencies are maximised out at every value of the parameter being
//! scanned, leaving a curve in that parameter alone. Splitting the search this way puts
//! the effort where the difficulty is: the frequencies are provably concave and never
//! needed searching, while about the noise parameters there is no proof either way
//! (`spec/parameter_prepass.md` §3.1).
//!
//! **Every rung is scored and no early exit is taken**, because nobody has shown the
//! curve has a single hump. That is the whole reason for a scan rather than a
//! one-dimensional optimiser, and it is why the cost is stated rather than avoided: the
//! generic ladder is 161 rungs, once per read group.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §4.2.

use std::collections::BTreeMap;

use smallvec::SmallVec;

use crate::ng::parameter_estimation::fitting::mixture_weights::{
    GenotypeLikelihoodTable, climb_mixture_weights,
};
use crate::ng::parameter_estimation::fitting::{NoiseModel, WeightedCell};
use crate::ng::types::{LogProb, Ploidy};

/// What one scan over a ladder of noise parameters returned.
///
/// Generic over the noise parameters, because the two paths scan different things: one
/// error rate here, three stutter parameters on the STR path.
#[derive(Clone, PartialEq, Debug)]
pub struct ScanResult<P> {
    /// The winning rung.
    pub noise: P,
    /// The genotype frequencies climbed to at that rung, **one set per ploidy the cells
    /// covered**. On the error-rate scan these are a means rather than an output — the
    /// scan is run for `noise` and they are discarded — while the sample's own rates
    /// come from a scan run for these.
    ///
    /// **Keyed by ploidy rather than the single vector the architecture sketches**
    /// (§5.2). A haploid region has two genotype classes and a diploid three,
    /// so they cannot share a weight vector and the scan climbs once per ploidy (§4.2);
    /// a single vector would mean picking one of them to report and dropping the rest,
    /// silently, in the module whose whole difficulty is that its wrong numbers have no
    /// symptom.
    ///
    /// The same quantity `MixtureWeightsFit::genotype_frequencies` carries and under the
    /// same name, because it is the same numbers: the climb returns them and the scan
    /// files them by ploidy.
    pub genotype_frequencies: BTreeMap<Ploidy, SmallVec<[f64; 3]>>,
    /// What makes "the best-scoring iterate" a defined comparison in an alternating
    /// fit. A [`LogProb`] rather than a bare `f64` because comparing is the only thing
    /// it is for, and `LogProb` carries `ln(0)` as `-∞` — the score of a rung nothing
    /// could have produced — where a linear probability would reach `0` and be
    /// indistinguishable from a value that merely got too small to represent.
    pub log_likelihood: LogProb,
    /// **Whether the answer sat on the ladder's edge, and it is not decoration.** A
    /// read group whose true rate lies outside the ladder — a bad run, heavy
    /// contamination, or any of the arithmetic failures a scan can suffer — has its
    /// answer silently clamped to an endpoint and emitted as though it were fitted,
    /// with a large observation count behind it. This one bit is what stands between a
    /// railed fit and a plausible-looking number.
    pub argmax_at_ladder_end: bool,
}

/// Step through `ladder`; at each rung climb to the genotype frequencies that best
/// explain `cells`, and return the best-scoring rung with the frequencies found there.
///
/// **Takes a slice of cells, not a histogram and not an iterator.** Not a histogram,
/// because the generic path's concrete table is not a shape the STR accumulator has —
/// that would make the "shared with the STR path" claim false. Not an iterator, because
/// the scan re-walks the cells once per rung, and a slice is re-walkable without asking
/// the caller for `Clone` on a map traversal or re-deriving the keys 161 times. The
/// caller materialises once and lends.
///
/// **Ploidy travels with each cell rather than being one argument**, because one noise
/// parameter is fitted across every ploidy its reads covered — a haploid sex chromosome
/// and the diploid autosomes were prepared by the same chemistry. Each cell is scored
/// against its own genotype set, and the frequencies are climbed **once per ploidy** on
/// that ploidy's cells, because a haploid cell has two genotypes to mix and a diploid
/// three. A rung's score is the sum over ploidies.
///
/// **Ties resolve to the last of the tied rungs**, which is the only rule a scan can
/// state without knowing which way its ladder runs. The generic ladder ascends in Phred
/// and so **descends** in error rate (`generic::error_rate_ladder`), so the last of a
/// tied run is the lowest error rate — which is what
/// `impl_plan/parameter_prepass_generic.md` D3 asks for, and
/// `a_tie_resolves_to_the_lower_error_rate_on_the_generic_ladder` is what binds the two
/// together, so that reversing the ladder cannot quietly reverse the rule.
///
/// **It reports whether its answer sat on the ladder's edge, and that is not
/// decoration.** A read group whose true rate lies outside the ladder — a bad run, heavy
/// contamination, or any of the arithmetic failures a scan can suffer — has its answer
/// silently clamped to an endpoint and emitted as though it were fitted, with a large
/// observation count behind it. That one bit is the only thing standing between a railed
/// fit and a plausible-looking number.
///
/// # Panics
///
/// If `cells` or `ladder` is empty, if any ploidy present carries no sites at all, if the
/// model declares no genotypes at a ploidy present, if it appends a number of entries
/// other than the [`NoiseModel::genotypes`] it declares, or if what it wrote is not a
/// log-likelihood — a `NaN`, a `+∞`, or a whole row of `−∞` on a cell that carries sites.
/// **The last three name the rung**, which is the part no later frame can recover: by the
/// time the climb refuses the table, the noise parameters are gone.
#[must_use]
pub fn fit_by_profile_scan<M>(
    model: &M,
    cells: &[M::Cell],
    ladder: &[M::NoiseParams],
) -> ScanResult<M::NoiseParams>
where
    M: NoiseModel,
    M::NoiseParams: Clone,
{
    assert!(!cells.is_empty(), "a scan needs at least one cell to score");
    assert!(!ladder.is_empty(), "a scan needs at least one rung to try");

    // Everything that does not move along the ladder is built once, before it: which
    // cells belong to each ploidy, what each of them weighs, and the interior start the
    // climb takes. `BTreeMap` rather than `HashMap` because the per-ploidy scores are
    // added into one total and floating-point addition is not associative, so the order
    // has to be the same on every run.
    let mut cells_of_ploidy: BTreeMap<Ploidy, Vec<usize>> = BTreeMap::new();
    for (position, cell) in cells.iter().enumerate() {
        cells_of_ploidy
            .entry(cell.ploidy())
            .or_default()
            .push(position);
    }

    let mut plans: Vec<PloidyPlan> = Vec::with_capacity(cells_of_ploidy.len());
    for (ploidy, positions) in cells_of_ploidy {
        let genotypes = model.genotypes(ploidy);
        assert!(
            genotypes > 0,
            "the model scores a cell of ploidy {ploidy} against no genotypes at all"
        );
        let cell_weights: Vec<f64> = positions
            .iter()
            .map(|&position| cells[position].sites() as f64)
            .collect();
        let total_sites: f64 = cell_weights.iter().sum();
        assert!(
            total_sites > 0.0,
            "every cell of ploidy {ploidy} holds zero sites, so there is nothing to fit"
        );
        plans.push(PloidyPlan {
            ploidy,
            genotypes,
            positions,
            cell_weights,
            uniform_start: vec![1.0 / genotypes as f64; genotypes],
        });
    }

    // Scratch, cleared and refilled once per (rung, ploidy) rather than allocated there:
    // the generic path's cell table has room for 583 cells and walks whatever it holds at
    // each of 161 rungs — the ladder's capacity is the bound, not the count, since only
    // non-empty cells are materialised and the attributed arm adds one entry per listing.
    // The model appends
    // straight into this, so what comes out is the row-major table the climb borrows —
    // no per-cell row and no copy.
    let mut ln_likelihood_row_major: Vec<f64> = Vec::new();

    let mut best: Option<ScanResult<M::NoiseParams>> = None;

    for (rung, noise) in ladder.iter().enumerate() {
        let mut rung_log_likelihood = 0.0;
        let mut genotype_frequencies: BTreeMap<Ploidy, SmallVec<[f64; 3]>> = BTreeMap::new();

        for plan in &plans {
            ln_likelihood_row_major.clear();
            for &position in &plan.positions {
                let before = ln_likelihood_row_major.len();
                model.append_genotype_likelihoods(
                    &cells[position],
                    noise,
                    plan.ploidy,
                    &mut ln_likelihood_row_major,
                );
                assert_eq!(
                    ln_likelihood_row_major.len() - before,
                    plan.genotypes,
                    "the model declares {} genotypes at ploidy {} and appended {} for \
                     cell {position}",
                    plan.genotypes,
                    plan.ploidy,
                    ln_likelihood_row_major.len() - before
                );
            }

            // **What the model wrote is checked here, where the rung is still in
            // scope.** The same faults are refused two frames down, by
            // `GenotypeLikelihoodTable::from_natural_logs` and by the climb — but
            // those messages name a row of a buffer and nothing else, and the rung is
            // the one part no later frame can recover, because the noise parameters
            // are gone by then. A model that goes wrong at one rung of 161 is exactly
            // the case worth localising.
            for ((row, &position), &cell_weight) in ln_likelihood_row_major
                .chunks_exact(plan.genotypes)
                .zip(&plan.positions)
                .zip(&plan.cell_weights)
            {
                for (genotype, &entry) in row.iter().enumerate() {
                    assert!(
                        entry.is_finite() || entry == f64::NEG_INFINITY,
                        "rung {rung}, ploidy {}: cell {position} scored {entry} under \
                         genotype {genotype}, which is not a log-likelihood",
                        plan.ploidy
                    );
                }
                // A cell that carries no sites may legally say no genotype produced
                // it — the climb never looks at it — so the weight is what makes this
                // a fault rather than a shape.
                assert!(
                    cell_weight == 0.0 || row.iter().any(|&entry| entry > f64::NEG_INFINITY),
                    "rung {rung}, ploidy {}: no genotype can have produced cell \
                     {position}, which holds {cell_weight} sites",
                    plan.ploidy
                );
            }

            let table = GenotypeLikelihoodTable::from_natural_logs(
                &ln_likelihood_row_major,
                plan.genotypes,
            );
            let climbed = climb_mixture_weights(table, &plan.cell_weights, &plan.uniform_start);
            // **A rung whose climb ran out of passes is scored below its own summit**,
            // so the rung beside it can win on that alone and the argmax stops being
            // the argmax. There is no channel to report it — `FitTermination` covers
            // the outer alternation of Milestone E and not this — and `MAX_CLIMB_PASSES`
            // is measured on four-cell fixtures where the slowest truth took 1,234
            // passes, against a real table of 583 cells climbed 161 times. Debug only,
            // because a slow climb is a data condition rather than a bug
            // (`spec/parameter_prepass.md` §3.1) and the release path must not abort on
            // one.
            debug_assert!(
                climbed.converged,
                "rung {rung}, ploidy {}: the climb used all {} passes, so this rung is \
                 scored short of its summit and the rung beside it may win on that alone",
                plan.ploidy, climbed.passes
            );
            rung_log_likelihood += climbed.log_likelihood.get();
            genotype_frequencies.insert(plan.ploidy, climbed.genotype_frequencies);
        }

        // The rung's answer is built here rather than at the end from a stored index,
        // so that the rail flag is decided where `rung` is in scope and there is no
        // second, private copy of `ScanResult` to keep in step with this one.
        let scanned = ScanResult {
            noise: noise.clone(),
            genotype_frequencies,
            log_likelihood: LogProb(rung_log_likelihood),
            argmax_at_ladder_end: rung == 0 || rung == ladder.len() - 1,
        };

        // `>=` and not `>`: a tie keeps the **later** rung. See the tie rule above — it
        // is what makes the generic ladder, which descends in error rate, resolve a tie
        // to the lower rate.
        let scores_at_least_as_well = best
            .as_ref()
            .is_none_or(|winner| scanned.log_likelihood >= winner.log_likelihood);
        if scores_at_least_as_well {
            best = Some(scanned);
        }
    }

    best.expect("the ladder is not empty, so some rung won")
}

/// Everything about one ploidy's cells that does not change along the ladder.
///
/// Built once before the rungs, because the cells, their weights and the width of their
/// genotype set are the same at every rung — only the noise parameters move.
struct PloidyPlan {
    ploidy: Ploidy,
    /// How many genotypes this ploidy's cells are scored against, as the model declares
    /// it. Not `ploidy + 1`: that is the SNP/indel path's count and not every path's.
    genotypes: usize,
    /// Where this ploidy's cells sit in the caller's slice.
    positions: Vec<usize>,
    /// One site count per entry of `positions`, in the same order.
    ///
    /// **Parallel to `positions` rather than zipped into it**, because
    /// [`climb_mixture_weights`] takes the weights as one contiguous `&[f64]`; a
    /// `Vec<(usize, f64)>` would have to be unzipped into a fresh buffer at every one of
    /// the ladder's 161 rungs. Nothing can be transposed between the two — one is a
    /// `usize` index and the other an `f64` weight — so what has to hold is only that
    /// they stay the same length and the same order, and both are built in one pass
    /// below.
    cell_weights: Vec<f64>,
    /// The interior point the climb starts from. Where it starts does not matter — the
    /// surface is concave — so it is built once rather than at every rung.
    uniform_start: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;

    use super::*;

    /// A cell of a toy table: which genotype set it belongs to, how many sites landed in
    /// it, and — for the fixture models below — which column of a hand-written matrix
    /// says how likely each genotype makes it.
    #[derive(Clone, Debug)]
    struct ToyCell {
        ploidy: Ploidy,
        sites: u64,
        row: usize,
    }

    impl WeightedCell for ToyCell {
        fn ploidy(&self) -> Ploidy {
            self.ploidy
        }
        fn sites(&self) -> u64 {
            self.sites
        }
    }

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("a positive copy number")
    }

    /// The rail flag is the field with teeth, so it is stated rather than left to a
    /// reader to notice: a scan that railed reports the same shape as one that did not.
    #[test]
    fn a_scan_result_reports_whether_its_answer_sat_on_the_ladders_edge() {
        let diploid = ploidy(2);
        let railed = ScanResult {
            noise: 0.1_f64,
            genotype_frequencies: BTreeMap::from([(
                diploid,
                SmallVec::from_slice(&[0.98, 0.015, 0.005]),
            )]),
            log_likelihood: LogProb(-1.2e9),
            argmax_at_ladder_end: true,
        };

        assert!(railed.argmax_at_ladder_end);
        assert_eq!(railed.genotype_frequencies[&diploid].len(), 3);
    }

    /// A model whose likelihoods come from a table indexed by rung and cell, so a test
    /// can state the whole likelihood surface and know what the answer must be.
    ///
    /// It also **records which rungs it was asked about**, which is how
    /// `every_rung_is_scored_and_none_is_skipped` checks the contract that no early exit
    /// is taken.
    struct TableModel {
        /// `[rung][row][genotype]` — the linear probability, not its log.
        likelihood: Vec<Vec<Vec<f64>>>,
        genotypes: usize,
        asked: RefCell<Vec<usize>>,
    }

    impl NoiseModel for TableModel {
        type Cell = ToyCell;
        type NoiseParams = usize;

        fn genotypes(&self, _ploidy: Ploidy) -> usize {
            self.genotypes
        }

        fn append_genotype_likelihoods(
            &self,
            cell: &ToyCell,
            noise: &usize,
            _ploidy: Ploidy,
            out: &mut Vec<f64>,
        ) {
            self.asked.borrow_mut().push(*noise);
            for likelihood in &self.likelihood[*noise][cell.row] {
                out.push(likelihood.ln());
            }
        }
    }

    /// Three rungs over four cells and three genotypes. Rung 1 is the truth: its columns
    /// are the ones the cell weights below are generated from.
    fn three_rung_table() -> Vec<Vec<Vec<f64>>> {
        vec![
            // rung 0 — a poor explanation of the table
            vec![
                vec![0.40, 0.30, 0.30],
                vec![0.30, 0.30, 0.20],
                vec![0.20, 0.20, 0.30],
                vec![0.10, 0.20, 0.20],
            ],
            // rung 1 — the truth
            vec![
                vec![0.70, 0.10, 0.02],
                vec![0.20, 0.40, 0.08],
                vec![0.07, 0.35, 0.30],
                vec![0.03, 0.15, 0.60],
            ],
            // rung 2 — poor in the other direction
            vec![
                vec![0.10, 0.20, 0.60],
                vec![0.20, 0.30, 0.30],
                vec![0.30, 0.30, 0.08],
                vec![0.40, 0.20, 0.02],
            ],
        ]
    }

    /// The exact cell weights an infinite genome at `truth` would produce under a set of
    /// `columns` — the same device the research harnesses use, so the answer is what the
    /// estimator converges to with no sampling noise in it.
    ///
    /// Takes the columns rather than a table and a rung index, so that the rail test can
    /// ask it for the weights under a truth that is on **no** rung of its ladder.
    fn weights_under(columns: &[Vec<f64>], truth: &[f64], sites: f64) -> Vec<u64> {
        columns
            .iter()
            .map(|row| {
                let probability: f64 = row
                    .iter()
                    .zip(truth)
                    .map(|(&likelihood, &frequency)| frequency * likelihood)
                    .sum();
                (sites * probability).round() as u64
            })
            .collect()
    }

    fn diploid_cells(sites: &[u64]) -> Vec<ToyCell> {
        sites
            .iter()
            .enumerate()
            .map(|(row, &sites)| ToyCell {
                ploidy: ploidy(2),
                sites,
                row,
            })
            .collect()
    }

    /// **The claim the step exists for: a table generated at a known rung recovers that
    /// rung.** The cell weights are the infinite-genome ones under rung 1, so rung 1 is
    /// the truth exactly, and the two rungs either side of it are wrong in opposite
    /// directions — a scan that stopped at the first improvement would take rung 0.
    #[test]
    fn a_table_generated_at_a_known_rung_recovers_it() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000_000.0));

        let scan = fit_by_profile_scan(&model, &cells, &[0, 1, 2]);

        assert_eq!(scan.noise, 1, "the middle rung is the truth");
        assert!(
            !scan.argmax_at_ladder_end,
            "the winner is interior, so nothing railed"
        );
        let frequencies = &scan.genotype_frequencies[&ploidy(2)];
        assert_eq!(frequencies.len(), 3);
        for (genotype, (&got, &want)) in frequencies.iter().zip(&[0.90, 0.07, 0.03]).enumerate() {
            assert!(
                (got - want).abs() < 1e-3,
                "genotype {genotype}: fitted {got}, truth {want}"
            );
        }
    }

    /// **Every rung is scored — no early exit.** The model records which rungs it was
    /// asked about, so this is checked directly rather than inferred from the answer.
    /// Nobody has shown the profile curve has a single hump, and a scan that stopped
    /// climbing would be sound only if it did.
    #[test]
    fn every_rung_is_scored_and_none_is_skipped() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));

        let _ = fit_by_profile_scan(&model, &cells, &[0, 1, 2]);

        let mut asked: Vec<usize> = model.asked.borrow().clone();
        asked.sort_unstable();
        asked.dedup();
        assert_eq!(asked, vec![0, 1, 2], "a rung went unscored");
    }

    /// A curve with two humps, which is the case the scan exists for and the case an
    /// optimiser that climbed would get wrong. Rung 0 is a local best — better than rung
    /// 1 beside it — and rung 3 is the global one.
    #[test]
    fn a_curve_with_two_humps_returns_the_higher_one() {
        // Four rungs. The truth is rung 3; rung 0 is made a local summit by putting a
        // dip at rung 1 between them.
        let truth = vec![
            vec![0.70, 0.10, 0.02],
            vec![0.20, 0.40, 0.08],
            vec![0.07, 0.35, 0.30],
            vec![0.03, 0.15, 0.60],
        ];
        let near = vec![
            vec![0.55, 0.15, 0.05],
            vec![0.25, 0.35, 0.10],
            vec![0.12, 0.30, 0.28],
            vec![0.08, 0.20, 0.57],
        ];
        let dip = vec![
            vec![0.25, 0.25, 0.25],
            vec![0.25, 0.25, 0.25],
            vec![0.25, 0.25, 0.25],
            vec![0.25, 0.25, 0.25],
        ];
        let table = vec![near.clone(), dip.clone(), dip, truth.clone()];
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[3], &[0.90, 0.07, 0.03], 1_000_000.0));

        // First prove the local summit is really there: over the first two rungs alone,
        // rung 0 wins. Without this the test below could pass on a table with one hump.
        let sub = fit_by_profile_scan(&model, &cells, &[0, 1]);
        assert_eq!(
            sub.noise, 0,
            "rung 0 is not a local summit, so there is one hump"
        );

        let scan = fit_by_profile_scan(&model, &cells, &[0, 1, 2, 3]);

        assert_eq!(scan.noise, 3, "the scan settled for the near rung");
        assert!(scan.argmax_at_ladder_end, "rung 3 is the ladder's last");
    }

    /// A ladder of five rungs whose columns march towards `truth` in equal steps, so the
    /// last rung is the closest the ladder gets and the truth itself is past its end.
    /// `reversed` marches the other way, so the closest rung is the first.
    fn ladder_marching_towards(truth: &[Vec<f64>], reversed: bool) -> Vec<Vec<Vec<f64>>> {
        const RUNGS: usize = 5;
        // **One over the number of *cells*, not the number of genotypes.** A column is one
        // genotype's distribution over the cells, so the flat column that sums to one is
        // `1/cells`. At `1/genotypes` the blend below stops preserving the sum and the
        // renormalisation is doing work it should not have to.
        let flat = 1.0 / truth.len() as f64;
        (0..RUNGS)
            .map(|rung| {
                let step = if reversed { RUNGS - 1 - rung } else { rung };
                // The last rung sits four fifths of the way from flat to the truth, so no
                // rung *is* the truth and the closest one is at the end — which is what
                // makes this a railed fit rather than a recovered one.
                let towards = 0.8 * step as f64 / (RUNGS - 1) as f64;
                let mut blended: Vec<Vec<f64>> = truth
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|&likelihood| flat + towards * (likelihood - flat))
                            .collect()
                    })
                    .collect();
                // Both `flat` and every column of `truth` sum to one, so a convex
                // combination of them does too and this loop is a no-op — kept because it
                // is what makes that a checked property of the fixture rather than an
                // assumption about `flat`.
                //
                // **The first version of this fixture used `1/genotypes` and skipped
                // it, and what went wrong is worth recording, because it is not what a
                // reader would guess.** Unnormalised columns did *not* hand the win to the
                // flattest rung: measured, that rung scored **worst of the five**, about
                // 160,000 nats below the winner. What the extra mass did was move the
                // argmax **one rung inward**, from 4 to 3 — a symptom far easier to read
                // as rounding than as a broken fixture.
                let genotypes = blended[0].len();
                for genotype in 0..genotypes {
                    let column: f64 = blended.iter().map(|row| row[genotype]).sum();
                    for row in &mut blended {
                        row[genotype] /= column;
                    }
                }
                blended
            })
            .collect()
    }

    /// **The rail flag, which is the only thing between a railed fit and a
    /// plausible-looking number.** A table whose truth lies past the ladder's end has
    /// its answer clamped to that end, and the answer looks exactly like a fitted one —
    /// same shape, same score, a large observation count behind it.
    ///
    /// Both ends are checked, because a flag that only watched one of them would pass a
    /// test written at the other. The contrast that gives this test teeth is
    /// `a_table_generated_at_a_known_rung_recovers_it`, where the winner is interior and
    /// the flag is **false** — without that, a flag hard-coded to `true` would pass here.
    #[test]
    fn a_table_whose_truth_lies_past_the_ladder_sets_the_rail_flag() {
        let truth_columns = vec![
            vec![0.70, 0.10, 0.02],
            vec![0.20, 0.40, 0.08],
            vec![0.07, 0.35, 0.30],
            vec![0.03, 0.15, 0.60],
        ];
        let truth = [0.90, 0.07, 0.03];

        for reversed in [false, true] {
            let table = ladder_marching_towards(&truth_columns, reversed);
            let model = TableModel {
                likelihood: table.clone(),
                genotypes: 3,
                asked: RefCell::new(Vec::new()),
            };
            // The weights come from the truth's own columns, which are on no rung.
            let sites = weights_under(&truth_columns, &truth, 1_000_000.0);

            let scan = fit_by_profile_scan(&model, &diploid_cells(&sites), &[0, 1, 2, 3, 4]);

            let expected_end = if reversed { 0 } else { 4 };
            assert_eq!(
                scan.noise, expected_end,
                "reversed = {reversed}: the ladder's closest rung should win"
            );
            assert!(
                scan.argmax_at_ladder_end,
                "reversed = {reversed}: the answer sat on the ladder's edge and said \
                 nothing"
            );
        }

        // And the degenerate ladder is all edge, which is why a one-rung ladder would
        // flag every read group — the failure `ERROR_RATE_LADDER_RUNGS` guards against.
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));
        assert!(fit_by_profile_scan(&model, &cells, &[0]).argmax_at_ladder_end);
    }

    /// **One noise parameter across every ploidy its reads covered, and one set of
    /// frequencies per ploidy.** A haploid cell has two genotypes to mix and a diploid
    /// three, so they cannot share a weight vector — the scan climbs once per ploidy and
    /// adds the two scores.
    #[test]
    fn cells_of_two_ploidies_are_scored_against_their_own_genotype_sets() {
        /// A model whose genotype count really does depend on the ploidy, so a scan that
        /// used one width for both would build a table of the wrong shape.
        struct DosageModel;
        impl NoiseModel for DosageModel {
            type Cell = ToyCell;
            type NoiseParams = f64;

            fn genotypes(&self, ploidy: Ploidy) -> usize {
                usize::from(ploidy.get()) + 1
            }

            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                noise: &f64,
                ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                // A crude but honest per-genotype likelihood: `row` stands for how many
                // reads showed the alternative allele out of four, and `noise` for the
                // per-read error rate.
                for alt_copies in 0..=ploidy.get() {
                    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
                    let p = carried * (1.0 - noise / 3.0) + (1.0 - carried) * noise;
                    let alt = cell.row as u32;
                    out.push(f64::from(alt).mul_add(p.ln(), f64::from(4 - alt) * (1.0 - p).ln()));
                }
            }
        }

        let cells: Vec<ToyCell> = (0..=4)
            .map(|row| ToyCell {
                ploidy: ploidy(2),
                sites: if row == 0 { 10_000 } else { 20 },
                row,
            })
            .chain((0..=4).map(|row| ToyCell {
                ploidy: ploidy(1),
                sites: if row == 0 { 5_000 } else { 10 },
                row,
            }))
            .collect();

        let scan = fit_by_profile_scan(&DosageModel, &cells, &[0.01, 0.001, 0.0001]);

        assert_eq!(
            scan.genotype_frequencies.len(),
            2,
            "one set of frequencies per ploidy"
        );
        assert_eq!(
            scan.genotype_frequencies[&ploidy(1)].len(),
            2,
            "a haploid has two genotype classes"
        );
        assert_eq!(
            scan.genotype_frequencies[&ploidy(2)].len(),
            3,
            "a diploid has three"
        );
        for (ploidy, weights) in &scan.genotype_frequencies {
            let total: f64 = weights.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-9,
                "ploidy {ploidy}: the frequencies sum to {total}"
            );
        }
    }

    /// **A tie resolves to the later rung, and on the generic ladder that is the lower
    /// error rate.** The two halves are asserted together on purpose: the scan cannot
    /// see which way its ladder runs, so the positional rule and the ladder's direction
    /// are one fact split across two files, and reversing either without the other would
    /// quietly reverse the answer.
    #[test]
    fn a_tie_resolves_to_the_lower_error_rate_on_the_generic_ladder() {
        use crate::ng::parameter_estimation::generic::error_rate_ladder;

        let ladder = error_rate_ladder();
        assert!(
            ladder[0].get() > ladder[ladder.len() - 1].get(),
            "the generic ladder descends in error rate, which is what makes the \
             positional tie rule the lower-rate one"
        );

        // Three rungs whose columns are identical, so every rung scores the same and the
        // whole ladder is one tie.
        let flat = vec![
            vec![0.70, 0.10, 0.02],
            vec![0.20, 0.40, 0.08],
            vec![0.07, 0.35, 0.30],
            vec![0.03, 0.15, 0.60],
        ];
        let table = vec![flat.clone(), flat.clone(), flat];
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[0], &[0.90, 0.07, 0.03], 1_000.0));

        let scan = fit_by_profile_scan(&model, &cells, &[0, 1, 2]);
        assert_eq!(scan.noise, 2, "a tie keeps the last of the tied rungs");
    }

    /// The score reported is the sum over ploidies at the winning rung, which is what
    /// makes "best-scoring iterate" a defined comparison in the coupled fit of E2.
    #[test]
    fn the_reported_score_is_the_sum_over_the_ploidies_at_the_winning_rung() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let sites = weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0);

        let diploid_only = fit_by_profile_scan(&model, &diploid_cells(&sites), &[0, 1, 2]);

        // The same cells again under a second ploidy label doubles the score, because
        // the two groups are scored independently and added.
        let mut doubled = diploid_cells(&sites);
        doubled.extend(sites.iter().enumerate().map(|(row, &sites)| ToyCell {
            ploidy: ploidy(4),
            sites,
            row,
        }));
        let both = fit_by_profile_scan(&model, &doubled, &[0, 1, 2]);

        assert_eq!(both.noise, diploid_only.noise);
        assert!(
            (both.log_likelihood.get() - 2.0 * diploid_only.log_likelihood.get()).abs() < 1e-6,
            "two identical groups scored {} against one group's {}",
            both.log_likelihood.get(),
            diploid_only.log_likelihood.get()
        );
    }

    /// A model that appends a different number of entries from the count it declares
    /// would reshape the table without changing its length, and the climb would run on
    /// transposed rows. Named rather than absorbed.
    #[test]
    #[should_panic(expected = "declares 3 genotypes at ploidy 2 and appended 2")]
    fn a_model_that_appends_a_width_it_did_not_declare_is_refused() {
        struct Liar;
        impl NoiseModel for Liar {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                3
            }

            fn append_genotype_likelihoods(
                &self,
                _cell: &ToyCell,
                _noise: &usize,
                _ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                out.extend([0.5_f64.ln(), 0.5_f64.ln()]);
            }
        }

        let _ = fit_by_profile_scan(&Liar, &diploid_cells(&[10, 20]), &[0]);
    }

    /// **The zero-sites guard is per ploidy, and this is what makes its scope
    /// testable.** With one ploidy in the scan, "this ploidy's cells" and "every cell in
    /// the slice" are the same set, so a guard summing the whole slice would pass every
    /// single-ploidy test. A read group covering a haploid contig that contributed no
    /// sites is the ordinary case — one error rate is fitted across every ploidy the
    /// group covered — and the guard is the only thing that names *which* ploidy was
    /// empty. Two frames down the message is "every cell carries zero weight", which
    /// names neither.
    #[test]
    #[should_panic(expected = "every cell of ploidy 1 holds zero sites")]
    fn a_ploidy_holding_no_sites_beside_one_that_does_names_that_ploidy() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let mut cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));
        cells.extend((0..4).map(|row| ToyCell {
            ploidy: ploidy(1),
            sites: 0,
            row,
        }));
        let _ = fit_by_profile_scan(&model, &cells, &[0, 1, 2]);
    }

    /// A model that scores a cell as `NaN` at one rung of three names **that rung**, not
    /// a row of a buffer nobody outside the loop can see. The rung is what no later frame
    /// can recover: by the time the climb refuses the table, the noise parameters are
    /// gone.
    #[test]
    #[should_panic(
        expected = "rung 2, ploidy 2: cell 3 scored NaN under genotype 0, which is not a"
    )]
    fn a_model_that_scores_a_cell_as_nan_names_the_rung_and_the_cell() {
        struct NanAtTheLastRung;
        impl NoiseModel for NanAtTheLastRung {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                3
            }

            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                noise: &usize,
                _ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                if *noise == 2 && cell.row == 3 {
                    out.extend([f64::NAN, 0.5_f64.ln(), 0.5_f64.ln()]);
                } else {
                    out.extend([0.5_f64.ln(), 0.3_f64.ln(), 0.2_f64.ln()]);
                }
            }
        }

        let _ = fit_by_profile_scan(
            &NanAtTheLastRung,
            &diploid_cells(&[10, 20, 30, 40]),
            &[0, 1, 2],
        );
    }

    /// A rung at which some weighted cell becomes impossible under every genotype is a
    /// fault in the model, and it is the rung that has to be named — the same cell is
    /// perfectly scorable at the other 160.
    #[test]
    #[should_panic(
        expected = "rung 2, ploidy 2: no genotype can have produced cell 3, which holds 40 sites"
    )]
    fn a_rung_at_which_a_weighted_cell_is_impossible_names_that_rung() {
        struct ImpossibleAtTheLastRung;
        impl NoiseModel for ImpossibleAtTheLastRung {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                3
            }

            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                noise: &usize,
                _ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                if *noise == 2 && cell.row == 3 {
                    out.extend([f64::NEG_INFINITY; 3]);
                } else {
                    out.extend([0.5_f64.ln(), 0.3_f64.ln(), 0.2_f64.ln()]);
                }
            }
        }

        let _ = fit_by_profile_scan(
            &ImpossibleAtTheLastRung,
            &diploid_cells(&[10, 20, 30, 40]),
            &[0, 1, 2],
        );
    }

    /// The companion to the two above, and what keeps the impossible-cell check from
    /// being tightened into one that refuses a legal table: a cell carrying **no** sites
    /// may say no genotype produced it, because the climb never looks at it.
    #[test]
    fn a_cell_holding_no_sites_may_say_no_genotype_produced_it() {
        struct ImpossibleWhereEmpty;
        impl NoiseModel for ImpossibleWhereEmpty {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                3
            }

            fn append_genotype_likelihoods(
                &self,
                cell: &ToyCell,
                _noise: &usize,
                _ploidy: Ploidy,
                out: &mut Vec<f64>,
            ) {
                if cell.sites == 0 {
                    out.extend([f64::NEG_INFINITY; 3]);
                } else {
                    out.extend([0.5_f64.ln(), 0.3_f64.ln(), 0.2_f64.ln()]);
                }
            }
        }

        let cells = diploid_cells(&[10, 0, 30, 40]);
        let scan = fit_by_profile_scan(&ImpossibleWhereEmpty, &cells, &[0, 1, 2]);
        assert!(scan.log_likelihood.get().is_finite());
    }

    /// A cell holding no sites is legal on its own — only the ploidy's total is
    /// checked — and it must not move the answer, because a cell nothing landed in is no
    /// evidence.
    #[test]
    fn a_cell_holding_no_sites_reaches_the_climb_without_moving_the_fit() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let sites = weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000_000.0);

        let without = fit_by_profile_scan(&model, &diploid_cells(&sites), &[0, 1, 2]);

        let mut with_empty = diploid_cells(&sites);
        with_empty.push(ToyCell {
            ploidy: ploidy(2),
            sites: 0,
            row: 0,
        });
        let with = fit_by_profile_scan(&model, &with_empty, &[0, 1, 2]);

        assert_eq!(with.noise, without.noise);
        assert_eq!(
            with.log_likelihood.get().to_bits(),
            without.log_likelihood.get().to_bits(),
            "an empty cell moved the score"
        );
    }

    /// A model that scores some ploidy against no genotypes at all is named where the
    /// ploidy is still in scope. Two frames down the message knows only that some table
    /// was empty.
    #[test]
    #[should_panic(expected = "scores a cell of ploidy 2 against no genotypes at all")]
    fn a_model_that_scores_a_ploidy_against_no_genotypes_at_all_is_refused() {
        struct NoGenotypes;
        impl NoiseModel for NoGenotypes {
            type Cell = ToyCell;
            type NoiseParams = usize;

            fn genotypes(&self, _ploidy: Ploidy) -> usize {
                0
            }

            fn append_genotype_likelihoods(
                &self,
                _cell: &ToyCell,
                _noise: &usize,
                _ploidy: Ploidy,
                _out: &mut Vec<f64>,
            ) {
            }
        }

        let _ = fit_by_profile_scan(&NoGenotypes, &diploid_cells(&[10, 20]), &[0]);
    }

    /// **The per-ploidy scores are added in a fixed order whatever order the cells
    /// arrive in**, which is what the `BTreeMap` is for: floating-point addition is not
    /// associative, so two runs that summed the ploidies in different orders would return
    /// different bits. Asserted on the bits and not against a tolerance, because a
    /// tolerance is exactly what this cannot afford — the scan compares rungs on this
    /// number and adjacent rungs can be close.
    ///
    /// **Six ploidies, and their site counts spread over four orders of magnitude.** Both
    /// are what make the test able to fail. At two groups a `HashMap` happens to iterate
    /// them in the same order often enough to pass; at six the chance is one in 720. And
    /// with equal group scores the sum is order-independent whatever the container does,
    /// so the counts are scaled ×1 to ×10⁵ to force the partial sums apart in the last
    /// bits.
    #[test]
    fn the_per_ploidy_scores_are_summed_in_a_fixed_order_whatever_the_slice_order() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let base = weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0);

        let groups: Vec<Vec<ToyCell>> = (1..=6_u8)
            .map(|copies| {
                let scale = 10_u64.pow(u32::from(copies) - 1);
                base.iter()
                    .enumerate()
                    .map(|(row, &sites)| ToyCell {
                        ploidy: ploidy(copies),
                        sites: sites * scale + 1,
                        row,
                    })
                    .collect()
            })
            .collect();

        let grouped: Vec<ToyCell> = groups.iter().flatten().cloned().collect();
        // The same cells, one from each ploidy in turn, so a container that iterated in
        // insertion order would add the six partial sums the other way round.
        let mut interleaved: Vec<ToyCell> = Vec::with_capacity(grouped.len());
        for row in 0..base.len() {
            for group in groups.iter().rev() {
                interleaved.push(group[row].clone());
            }
        }

        let first = fit_by_profile_scan(&model, &grouped, &[0, 1, 2]);
        let second = fit_by_profile_scan(&model, &interleaved, &[0, 1, 2]);

        assert_eq!(first.genotype_frequencies.len(), 6);
        assert_eq!(
            first.log_likelihood.get().to_bits(),
            second.log_likelihood.get().to_bits(),
            "the ploidies were summed in a different order: {} against {}",
            first.log_likelihood.get(),
            second.log_likelihood.get()
        );
    }

    /// A ladder carrying the same rung twice is the other way a tie arises, and the one
    /// an alternating fit re-scanning a narrowed ladder would produce. The rule is the
    /// same: the later of the two wins.
    #[test]
    fn a_ladder_with_a_repeated_rung_returns_the_later_of_the_two() {
        let table = three_rung_table();
        let model = TableModel {
            likelihood: table.clone(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let cells = diploid_cells(&weights_under(&table[1], &[0.90, 0.07, 0.03], 1_000.0));

        // Rung 1 is the truth and appears at positions 1 and 2 of the ladder.
        let scan = fit_by_profile_scan(&model, &cells, &[0, 1, 1]);

        assert_eq!(scan.noise, 1);
        assert!(
            scan.argmax_at_ladder_end,
            "the later of the two tied rungs is the ladder's last"
        );
    }

    #[test]
    #[should_panic(expected = "at least one cell")]
    fn a_scan_over_no_cells_is_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_profile_scan(&model, &[], &[0]);
    }

    #[test]
    #[should_panic(expected = "at least one rung")]
    fn a_scan_over_no_rungs_is_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_profile_scan(&model, &diploid_cells(&[10, 20]), &[]);
    }

    #[test]
    #[should_panic(expected = "every cell of ploidy 2 holds zero sites")]
    fn a_ploidy_whose_cells_hold_no_sites_is_refused() {
        let model = TableModel {
            likelihood: three_rung_table(),
            genotypes: 3,
            asked: RefCell::new(Vec::new()),
        };
        let _ = fit_by_profile_scan(&model, &diploid_cells(&[0, 0]), &[0, 1, 2]);
    }
}
