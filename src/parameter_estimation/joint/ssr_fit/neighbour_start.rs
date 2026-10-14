//! **One walk from a neighbouring class's answer**, for measuring whether it ends where the best of
//! the usual starting points ends, in fewer rounds — the owner's first idea for a faster fit
//! (2026-10-10). A measurement helper; nothing in the fit calls it.

use super::*;

/// How one walk from a neighbour's answer went.
#[derive(Debug, Clone, PartialEq)]
pub struct NeighbourWalk {
    /// How it ended and how many rounds it took.
    pub record: WalkRecord,
    /// The mean log-likelihood a tract where it stopped.
    pub log_likelihood_a_tract: f64,
    /// Its slippage numbers where it stopped, per slippage group.
    pub slippage: Vec<Slippage>,
    /// Its allele-length frequencies and their spread where it stopped.
    pub allele_length_variation: Vec<f64>,
    pub concentration: f64,
}

/// **Climbs `evidence` once, from `neighbour`'s answer**: its slippage numbers, its allele-length
/// frequencies and their spread. A group the neighbour has no numbers for starts at the first
/// usual starting point's level; if the neighbour has no concentration, the walk starts at that
/// point's concentration. `None` when no group put a read in `evidence`.
#[must_use]
pub fn walk_from_a_neighbour(
    evidence: &StratumEvidence,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    neighbour: &StratumFit,
) -> Option<NeighbourWalk> {
    let genotypes = genotype_pairs(allele_classes_of(config));
    let live_groups = evidence.groups_with_reads();
    if !live_groups.iter().any(|live| *live) {
        return None;
    }
    let fallback = Parameters::start(
        config.starting_points[0],
        evidence.groups,
        allele_classes_of(config),
    );
    let parameters = Parameters {
        slippage: (0..evidence.groups)
            .map(|group| {
                neighbour
                    .slippage
                    .get(group)
                    .copied()
                    .flatten()
                    .unwrap_or(fallback.slippage[group])
            })
            .collect(),
        allele_length_variation: neighbour.allele_length_variation.clone(),
        concentration: neighbour.concentration.unwrap_or(fallback.concentration),
    };
    let climb = climb_from_parameters(
        evidence,
        parameters,
        homozygote_excess,
        &genotypes,
        &live_groups,
        config,
        WhereTheThreadsGo::AcrossStrata,
    );
    Some(NeighbourWalk {
        record: climb.record(),
        log_likelihood_a_tract: climb.score,
        slippage: climb.parameters.slippage,
        allele_length_variation: climb.parameters.allele_length_variation,
        concentration: climb.parameters.concentration,
    })
}

/// **The mean log-likelihood a tract of `evidence` at these numbers**, its average over each tract's
/// allele frequencies taken over `points` points rather than the fit's 256 — for scoring two answers
/// with an average precise enough to tell them apart (`fit_precision_quadrature_average_2026-10-01.md`).
#[must_use]
pub fn mean_log_likelihood_at(
    evidence: &StratumEvidence,
    slippage: &[Slippage],
    allele_length_variation: &[f64],
    concentration: f64,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    points: usize,
) -> f64 {
    let genotypes = genotype_pairs(allele_classes_of(config));
    let precise = SsrFitConfig {
        quadrature_points: points,
        ..config.clone()
    };
    let mut scorer = Scorer::new(
        evidence,
        homozygote_excess,
        &genotypes,
        &precise,
        WhereTheThreadsGo::AcrossTheTractsOfOneStratum,
    );
    scorer.score(&Parameters {
        slippage: slippage.to_vec(),
        allele_length_variation: allele_length_variation.to_vec(),
        concentration,
    })
}
