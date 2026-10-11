//! **Each sample's genotype likelihoods at each tract of a class, with the stutter held** — the
//! input to estimating a class's allele-length frequencies without the average over tract
//! frequencies (the owner's "route 2" research, 2026-10-11). A measurement helper; nothing in the
//! fit calls it.

use super::*;

/// One class's genotype likelihoods: a row a (tract, sample with reads), over the genotype pairs of
/// [`genotype_pairs`]`(classes)`, each row divided by its largest value.
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypeRows {
    /// How many allele classes, `2 × allele_span + 1`.
    pub classes: usize,
    /// The genotype pairs, in the order a row's values are.
    pub genotypes: Vec<(usize, usize)>,
    /// Per row, which tract of the class it is at, counting only tracts with reads.
    pub tract: Vec<u32>,
    /// Per row, the sample's homozygote excess.
    pub homozygote_excess: Vec<f64>,
    /// The rows, laid end to end: row `r` is `values[r * genotypes.len()..][..genotypes.len()]`.
    pub values: Vec<f64>,
}

/// **Every sample's genotype likelihoods at every tract of `evidence`**, its reads scored under
/// `slippage` — the stutter held where a fit left it.
#[must_use]
pub fn genotype_rows(
    evidence: &StratumEvidence,
    slippage: &[Slippage],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> GenotypeRows {
    let classes = allele_classes_of(config);
    let genotypes = genotype_pairs(classes);
    let distributions = ReadDistributionsByLevel::of(
        slippage,
        &evidence.libraries,
        config.allele_span,
        evidence.read_span,
    );
    let mut rows = GenotypeRows {
        classes,
        genotypes: genotypes.clone(),
        tract: Vec::new(),
        homozygote_excess: Vec::new(),
        values: Vec::new(),
    };
    for (index, tract) in evidence
        .tracts
        .iter()
        .filter(|tract| !tract.samples.is_empty())
        .enumerate()
    {
        let likelihoods =
            TractLikelihoods::of(tract, &distributions, &genotypes, homozygote_excess);
        rows.values.extend_from_slice(&likelihoods.scaled);
        rows.homozygote_excess
            .extend_from_slice(&likelihoods.homozygote_excess);
        rows.tract
            .extend(std::iter::repeat_n(index as u32, tract.samples.len()));
    }
    rows
}
