//! The two ways into step 4's generic path, and everything the fits need that is not in
//! the loci themselves.
//!
//! **Two entry points because the step has two callers.** A run that walks the loci for
//! step 4 alone wants one function; a run that folds step 4 into an existing consumer loop
//! drove the accumulator itself and wants only the reduction. The first is the second over
//! an accumulator fed by the stream, so the two cannot answer differently — which is the
//! property, not the convenience.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §1.1.

use std::collections::BTreeMap;
use std::sync::Arc;

use crate::ng::locus_generation::{LocusGenerationError, SampleLocusObservations};
use crate::ng::parameter_estimation::generic::accumulators::{
    GenericAccumulators, InbreedingMode, PloidyMap,
};
use crate::ng::parameter_estimation::generic::coupled_fit::{fit_coupled, library_shares};
use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use crate::ng::parameter_estimation::generic::fallback::take_supplied_inbreeding;
use crate::ng::parameter_estimation::generic::noise_model::{LibraryNoise, SampleLibraryNoise};
use crate::ng::parameter_estimation::generic::runs::{RunsModelStarts, fit_inbreeding};
use crate::ng::parameter_estimation::generic::{
    GenericSampleParameters, SampleRates, error_rate_ladder,
};
use crate::ng::parameter_estimation::{Estimate, ParameterEstimationError};
use crate::ng::read::filtering::ReadFilterConfig;
use crate::ng::types::{ErrorRate, Ploidy, ReadGroupId};

/// Everything the fits need that is not in the loci themselves.
pub struct GenericEstimationConfig {
    /// The sample's name as the alignment files declared it — only ever used to name the
    /// sample in an error message and in the emitted summary.
    ///
    /// **A `String`, not a newtype and not `SampleIdentity`**: ng has no `SampleId`, and the
    /// one identity type it does have (`read/input/mod.rs`) compares files by pointer, which
    /// answers *are these the same open sample?* rather than *what is this sample called*.
    pub sample_name: String,
    /// The sample's own read groups, which decide whether a site's alternative reads keep
    /// the library they came from. At one group no site ever takes the attributed arm.
    pub read_groups: Vec<ReadGroupId>,
    pub ploidy: Arc<dyn PloidyMap>,
    pub inbreeding: InbreedingMode,
    /// Error rates the run states rather than fits — the **third** rung of the fallback
    /// ladder, below *fitted here* and below *borrowed*.
    ///
    /// **The name is the decision** (owner, 2026-08-09). `fallback_error_rates` says these
    /// apply only where the sample's own data could not answer, which is the order
    /// `resolve_error_rates` implements. A field named `error_rate_overrides` would be a
    /// different feature and would have to sit *above* *fitted here*, not below *borrowed*
    /// — for an operator who supplies a rate precisely because this library is chemically
    /// unlike its siblings.
    pub fallback_error_rates: BTreeMap<ReadGroupId, ErrorRate>,
    /// The binning rule, shared by every accumulator of every shard so that their cells
    /// mean the same thing and `merge` can prove it by pointer identity.
    pub edges: Arc<DepthBinEdges>,
    /// The read-admission policy the loci were produced under, recorded so that every
    /// emitted `ε` says what population of reads it describes
    /// (`spec/parameter_prepass.md` §2).
    ///
    /// **Recorded and not yet read.** An error rate is a property of the reads that were
    /// admitted, so two runs at different `min_mapq` produce rates that are not comparable;
    /// what carries this into the emitted summary is the `SampleSummary` assembly, which
    /// belongs to the cohort gather and is out of this plan's scope (spec §7). It is on the
    /// config rather than waiting for that plan because the config is where a caller states
    /// it, and a field added later would be one every existing caller had to be revisited to
    /// fill.
    pub read_admission: ReadFilterConfig,
}

impl GenericEstimationConfig {
    /// A fresh accumulator for one region shard of this sample.
    ///
    /// **Every shard of a sample must come from here**, because `merge` requires the same
    /// `edges` object and proves it by `Arc::ptr_eq` rather than by comparing lengths.
    #[must_use]
    pub fn accumulators(&self) -> GenericAccumulators {
        GenericAccumulators::new(
            Arc::clone(&self.edges),
            &self.read_groups,
            Arc::clone(&self.ploidy),
            self.inbreeding,
        )
    }
}

/// Walk a sample's loci and return its generic parameters — the whole step in one call,
/// for a caller with nothing else to do with the stream.
///
/// **This is [`GenericAccumulators::estimate`] over an accumulator fed by the stream**, and
/// it is written that way rather than duplicating the reduction so that the two entry points
/// cannot diverge.
///
/// # Errors
///
/// A [`LocusGenerationError`] in the stream is **fatal and propagates**, as
/// [`ParameterEstimationError::LocusGeneration`]. The loci a walk failed to produce are
/// missing evidence, not zero evidence, and a rate fitted over a truncated genome is wrong
/// in a way nothing downstream would announce.
///
/// Otherwise as [`GenericAccumulators::estimate`].
pub fn estimate_generic_parameters(
    loci: impl Iterator<Item = Result<SampleLocusObservations, LocusGenerationError>>,
    config: &GenericEstimationConfig,
) -> Result<GenericSampleParameters, ParameterEstimationError> {
    let mut accumulators = config.accumulators();
    for locus in loci {
        let locus = locus.map_err(|failure| ParameterEstimationError::LocusGeneration {
            sample: config.sample_name.clone(),
            cause: failure.to_string(),
        })?;
        accumulators.add_locus(&locus);
    }
    accumulators.estimate(config)
}

impl GenericAccumulators {
    /// The reduction, for a caller that drove the accumulator itself — one per region shard,
    /// merged. **The half that does no I/O.**
    ///
    /// The order is the design's dependency order and not a preference: the coupled fit
    /// settles the error rates and the genotype frequencies together, and the runs model
    /// then takes **both** — each library's own rate rather than one pooled one, and the
    /// sample's frequencies as the outside state's starting point.
    ///
    /// # Errors
    ///
    /// [`ParameterEstimationError::GenotypeFrequenciesNotFittable`] when some ploidy holds
    /// too few sites. [`ParameterEstimationError::InbreedingNotFittable`],
    /// [`InbreedingStatesNotSeparated`](ParameterEstimationError::InbreedingStatesNotSeparated)
    /// or [`InbreedingStartsDisagree`](ParameterEstimationError::InbreedingStartsDisagree)
    /// when `F` was to be fitted and could not be. **None of those has a default**: `F` is
    /// the parameter that differs most between an outcrosser and a selfing landrace, and a
    /// cohort's diversity divides by `1 − F`, so a wrong constant would be amplified rather
    /// than absorbed. The answer to all three is to supply one.
    ///
    /// # Panics
    ///
    /// If the sample has no read group with reads.
    pub fn estimate(
        &self,
        config: &GenericEstimationConfig,
    ) -> Result<GenericSampleParameters, ParameterEstimationError> {
        let ladder = error_rate_ladder();
        let coupled = fit_coupled(
            &config.sample_name,
            self,
            &ladder,
            &config.fallback_error_rates,
        )?;

        let (inbreeding, runs_model) = match config.inbreeding {
            // **A supplied `F` is taken without a fit, and that is not an optimisation.**
            // `InbreedingMode::Supplied` dropped the window key at accumulation, so there is
            // no chain to walk here even if a caller wanted one.
            InbreedingMode::Supplied(_) => (
                take_supplied_inbreeding(config.inbreeding, self.covered_positions()),
                None,
            ),
            InbreedingMode::Fitted => {
                self.fit_inbreeding_if_diploid(config, &coupled.rates, &coupled.error_rate)?
            }
        };

        Ok(GenericSampleParameters {
            error_rate: coupled.error_rate,
            rates: coupled.rates,
            inbreeding,
            runs_model,
            coupled_fit: coupled.termination,
        })
    }

    /// The runs model, run only where there is a diploid region to run it on.
    ///
    /// **`None` rather than an error above and below two copies**, which is the one place
    /// this step declines to answer without that being a failure: above two, `F` needs
    /// several identity-by-descent coefficients and is deferred (spec §7); below two there
    /// are no heterozygotes to be short of. A sample whose whole genome is haploid therefore
    /// estimates its other three parameters and reports no `F`.
    fn fit_inbreeding_if_diploid(
        &self,
        config: &GenericEstimationConfig,
        rates: &BTreeMap<Ploidy, Estimate<SampleRates>>,
        settled_rates: &BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
    ) -> Result<
        (
            Option<Estimate<crate::ng::types::InbreedingF>>,
            Option<crate::ng::parameter_estimation::generic::runs::RunsModelFit>,
        ),
        ParameterEstimationError,
    > {
        let diploid = Ploidy::try_new(2).expect("two is a positive copy number");
        let Some(outside_rates) = rates.get(&diploid) else {
            return Ok((None, None));
        };

        let noise = self.library_noise(settled_rates);
        let (estimate, fit) = fit_inbreeding(
            &config.sample_name,
            self.windowed_histograms(),
            &noise,
            diploid,
            &outside_rates.value,
            &RunsModelStarts::default(),
        )?;
        Ok((Some(estimate), Some(fit)))
    }

    /// Every library's share of the sample's reads paired with the rate that library was
    /// settled at.
    ///
    /// **Paired by read group and never by position**, which is the same discipline the
    /// coupled fit's own `noise_from` keeps and for the same reason: a rule with two
    /// libraries' rates swapped between their groups is still a probability over the cell
    /// space, so none of the scoring rule's identities can see it. The only thing that
    /// prevents it is that a share and a rate reach [`LibraryNoise`] under one key.
    ///
    /// **The rates are the *settled* ones, so a borrowed or supplied rate reaches the runs
    /// model** — the alternative is a library scored against a rate nobody chose.
    ///
    /// # Panics
    ///
    /// If a library has a share of the reads but no settled rate. `resolve_error_rates`
    /// answers for every group that contributed a read, and a share exists only for such a
    /// group, so the two key sets are the same set by construction — asserted rather than
    /// assumed, because a library silently dropped here would be scored against the
    /// share-weighted mean of the others.
    fn library_noise(
        &self,
        settled_rates: &BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
    ) -> SampleLibraryNoise {
        SampleLibraryNoise::new(
            library_shares(self.read_group_histograms())
                .into_iter()
                .map(|(read_group, share_of_reads)| LibraryNoise {
                    read_group,
                    share_of_reads,
                    error_rate: settled_rates
                        .get(&read_group)
                        .unwrap_or_else(|| {
                            panic!(
                                "read group {} produced {share_of_reads} of the sample's reads \
                             and was given no error rate",
                                read_group.get()
                            )
                        })
                        .value,
                }),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::parameter_estimation::generic::MIN_SITES_TO_FIT;
    use crate::ng::parameter_estimation::generic::accumulators::ConstantPloidy;
    use crate::ng::types::{ContigId, GenomeRegion, InbreedingF, Position};

    /// One diploid site at ten reads, `alt_reads` of them showing an alternative base.
    ///
    /// A one-base locus, so `region.len()` is one and the covered-position count is the
    /// site count — which is what lets the supplied-`F` warrant below be checked against a
    /// number a reader can compute.
    fn site(start: u64, alt_reads: u32, group: ReadGroupId) -> SampleLocusObservations {
        const DEPTH: u32 = 10;
        let mut observations = Vec::new();
        if alt_reads > 0 {
            observations.push(observation(b"C", group, alt_reads));
        }
        if DEPTH > alt_reads {
            observations.push(observation(b"A", group, DEPTH - alt_reads));
        }
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(start),
                end: Position(start),
            },
            reference_bases: b"A".as_slice().into(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    fn observation(bases: &[u8], read_group: ReadGroupId, num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.into(),
            read_witness: ReadWitness::Complete,
            read_group,
            num_obs,
            num_fwd: 0,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    /// A sample of `MIN_SITES_TO_FIT` sites, one in twenty of them heterozygous.
    ///
    /// **Above the floor by one site and not by a comfortable margin**, so that a fit
    /// rejected for too little data is a visible failure here rather than a silent one two
    /// milestones later.
    fn a_samples_loci(group: ReadGroupId) -> Vec<SampleLocusObservations> {
        (0..=MIN_SITES_TO_FIT)
            .map(|index| {
                let alt_reads = if index.is_multiple_of(20) { 5 } else { 0 };
                site(index + 1, alt_reads, group)
            })
            .collect()
    }

    fn config(inbreeding: InbreedingMode) -> GenericEstimationConfig {
        GenericEstimationConfig {
            sample_name: "SL_landrace_07".to_string(),
            read_groups: vec![ReadGroupId(1)],
            ploidy: Arc::new(ConstantPloidy(
                Ploidy::try_new(2).expect("a positive copy number"),
            )),
            inbreeding,
            fallback_error_rates: BTreeMap::new(),
            edges: Arc::new(DepthBinEdges::new()),
            read_admission: ReadFilterConfig::default(),
        }
    }

    fn supplied(value: f64) -> InbreedingMode {
        InbreedingMode::Supplied(InbreedingF::try_new(value).expect("a fraction"))
    }

    /// **The property F1 exists for: the two entry points cannot answer differently.**
    ///
    /// One walks the stream itself; the other is handed an accumulator the test drove. They
    /// are asserted equal as whole `GenericSampleParameters`, not field by field, so a
    /// field added later is covered without this test being revisited.
    ///
    /// **What this cannot say**: the two share `estimate`, so it proves the stream-walking
    /// half feeds the accumulator the same loci — not that the reduction is right. What the
    /// reduction computes is F2's question.
    #[test]
    fn the_two_entry_points_answer_identically() {
        let group = ReadGroupId(1);
        let config = config(supplied(0.2));

        let walked =
            estimate_generic_parameters(a_samples_loci(group).into_iter().map(Ok), &config)
                .expect("a sample above the site floor is fittable");

        let mut driven = config.accumulators();
        for locus in a_samples_loci(group) {
            driven.add_locus(&locus);
        }
        let reduced = driven
            .estimate(&config)
            .expect("the same sample, the same fit");

        assert_eq!(walked, reduced);
    }

    /// **A failed walk stops the estimate rather than shortening it.** The loci a walk did
    /// not produce are missing evidence, not zero evidence: a rate fitted over the prefix
    /// that did arrive is a plausible number describing a genome nobody chose.
    ///
    /// The stream below yields the whole sample and *then* fails, so the prefix is above
    /// every floor and would have fitted — which is what makes this a test of the propagation
    /// rather than of the floors.
    #[test]
    fn a_failed_walk_propagates_rather_than_fitting_the_prefix() {
        let group = ReadGroupId(1);
        let config = config(supplied(0.2));
        let stream = a_samples_loci(group)
            .into_iter()
            .map(Ok)
            .chain(std::iter::once(Err(LocusGenerationError::ForeignSample {
                region: GenomeRegion {
                    contig: ContigId(0),
                    start: Position(1),
                    end: Position(1),
                },
            })));

        let error = estimate_generic_parameters(stream, &config)
            .expect_err("a walk that failed must not return a fit over what arrived first");

        assert!(
            matches!(error, ParameterEstimationError::LocusGeneration { .. }),
            "the wrong failure: {error}"
        );
        assert!(
            error.to_string().contains("SL_landrace_07"),
            "a cohort run of hundreds needs the sample named: {error}"
        );
    }

    /// **A supplied `F` is reported as supplied and no chain is walked** — the first
    /// production caller `take_supplied_inbreeding` has ever had.
    ///
    /// Its warrant is the reference positions the sample covered, which for these one-base
    /// loci is the site count. Asserted against that number rather than against itself, so a
    /// warrant taken from the wrong table would show.
    #[test]
    fn a_supplied_inbreeding_coefficient_is_taken_without_a_fit() {
        let group = ReadGroupId(1);
        let config = config(supplied(0.4));

        let parameters =
            estimate_generic_parameters(a_samples_loci(group).into_iter().map(Ok), &config)
                .expect("fittable");

        let inbreeding = parameters
            .inbreeding
            .expect("a supplied coefficient is reported");
        assert_eq!(inbreeding.value.get(), 0.4);
        assert_eq!(inbreeding.provenance, Provenance::Supplied);
        assert_eq!(inbreeding.observations, MIN_SITES_TO_FIT + 1);
        assert!(
            parameters.runs_model.is_none(),
            "no chain was walked, so there is nothing to report about one"
        );
    }

    /// **The `Supplied` rung's first production source** — `fallback_error_rates` had no
    /// caller before F1, so nothing had ever exercised it non-empty.
    ///
    /// **Reaching that rung takes a specific shape and it is worth saying why.** The ladder
    /// is *fitted here → borrowed → supplied → defaulted*, so a supplied rate is only
    /// consulted when the group's own sites are too few **and** no sibling group qualifies to
    /// lend. Three read groups each covering a disjoint third of the sites, each one site
    /// below the floor, is that shape: every group is too thin to fit and none can lend,
    /// while the whole-sample table — which is the union — is comfortably above the floor the
    /// genotype frequencies are checked against. Two groups are supplied a rate and the third
    /// is not, so the last two rungs are both reached in one sample.
    #[test]
    fn a_supplied_rate_is_used_where_no_group_can_fit_or_lend() {
        let thin = MIN_SITES_TO_FIT - 1;
        let groups = [ReadGroupId(1), ReadGroupId(2), ReadGroupId(3)];

        let mut loci = Vec::new();
        for (index, &group) in groups.iter().enumerate() {
            let first = index as u64 * thin + 1;
            for offset in 0..thin {
                let position = first + offset;
                let alt_reads = if position.is_multiple_of(20) { 5 } else { 0 };
                loci.push(site(position, alt_reads, group));
            }
        }
        loci.sort_by_key(|locus| locus.region.start.get());

        let mut config = config(supplied(0.2));
        config.read_groups = groups.to_vec();
        config.fallback_error_rates = BTreeMap::from([
            (groups[0], ErrorRate::try_new(0.002).expect("a probability")),
            (groups[1], ErrorRate::try_new(0.004).expect("a probability")),
        ]);

        let parameters = estimate_generic_parameters(loci.into_iter().map(Ok), &config)
            .expect("the union of three thin groups is above the site floor");

        assert_eq!(
            parameters.error_rate[&groups[0]].provenance,
            Provenance::Supplied
        );
        assert_eq!(parameters.error_rate[&groups[0]].value.get(), 0.002);
        assert_eq!(parameters.error_rate[&groups[1]].value.get(), 0.004);
        assert_eq!(
            parameters.error_rate[&groups[2]].provenance,
            Provenance::Defaulted,
            "the group nothing was supplied for falls through to the last rung"
        );
    }
}
