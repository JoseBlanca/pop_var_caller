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
}

/// **Climbs `evidence` once, from `neighbour`'s answer**: its slippage numbers, its allele-length
/// frequencies and their spread. A group the neighbour has no numbers for starts at the first
/// usual starting point's level. `None` when no group put a read in `evidence`.
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
        length_spectrum: neighbour.length_spectrum.clone(),
        concentration: neighbour.concentration,
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
    })
}
