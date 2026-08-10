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
use crate::ng::parameter_estimation::generic::coupled_fit::{fit_coupled, library_shares_over};
use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use crate::ng::parameter_estimation::generic::fallback::{
    as_marginal_rates, take_supplied_inbreeding,
};
use crate::ng::parameter_estimation::generic::noise_model::{LibraryNoise, SampleLibraryNoise};
use crate::ng::parameter_estimation::generic::runs::{RunsModelStarts, fit_inbreeding};
use crate::ng::parameter_estimation::generic::{
    GenericSampleParameters, SampleRates, SiteNoise, error_rate_ladder,
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
/// [`ParameterEstimationError::LocusGeneration`], carrying the walk's own error so that a
/// caller can see *which* stage broke and over which region. The loci a walk failed to
/// produce are missing evidence, not zero evidence, and a rate fitted over a truncated
/// genome is wrong in a way nothing downstream would announce.
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
            source: failure,
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
    /// **Any of those three discards the three parameters that *were* fitted, and that is
    /// deliberate** (owner, 2026-08-09). The error rates and the genotype frequencies are
    /// settled before the runs model runs, and they are dropped with them — on the streaming
    /// entry point, after a whole-genome walk. The reason is downstream: **`F` is a prior the
    /// calling step needs**, so a sample with the other three and no `F` cannot be called
    /// anyway, and returning it would only move the failure to a place with less context. It
    /// is worth knowing what that costs before running a cohort: on a genome with no runs the
    /// refusal rate is about **nine in twenty-three** (research note §6.3), and the two floors
    /// are far apart — [`MIN_SITES_TO_FIT`](super::MIN_SITES_TO_FIT) is 10,000 sites where
    /// [`MIN_WINDOWS_TO_FIT_INBREEDING`](super::runs::MIN_WINDOWS_TO_FIT_INBREEDING) is 3,000
    /// windows, which is 300 Mb that must hold sites. A region-restricted run, and an
    /// outcrossing cohort run without a supplied `F`, will meet this. **The answer in both
    /// cases is `InbreedingMode::Supplied`**, which is per-sample.
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
            InbreedingMode::Fitted => self.fit_inbreeding_if_diploid(
                config,
                &coupled.rates,
                &coupled.error_rate,
                coupled.site_noise,
            )?,
        };

        Ok(GenericSampleParameters {
            // **The marginal is applied here and nowhere earlier.** Everything inside the
            // fit — the coupled alternation, and the runs model below — scores with the
            // *clean* rates and the second class beside them, which is the pair the
            // scoring rule takes. Folding them into one number is what a *consumer* needs,
            // and doing it sooner would hand the runs model a single rate whose site
            // classes had been averaged away, putting back inside it exactly the bias the
            // second class exists to remove.
            error_rate: as_marginal_rates(coupled.error_rate, coupled.site_noise),
            rates: coupled.rates,
            inbreeding,
            runs_model,
            site_noise: coupled.site_noise,
            site_noise_off_the_ladder: coupled.site_noise_off_the_ladder,
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
        site_noise: Option<SiteNoise>,
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

        let noise = self.library_noise(settled_rates, diploid, site_noise);
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
    /// **And the shares are restricted to the ploidy being fitted**, which pooling them
    /// across the sample is not. The runs model walks the diploid windows and nothing else,
    /// so a library that contributed only to a haploid arm must not be given a share of the
    /// reads it is scoring: measured on a fixture whose two arms come from different
    /// libraries, pooling reports 0.5 and 0.5 where the diploid arm's reads are entirely one
    /// of them, and at Phred 20 against Phred 30 that puts the share-weighted rate 5.5 times
    /// above the truth.
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
        at: Ploidy,
        site_noise: Option<SiteNoise>,
    ) -> SampleLibraryNoise {
        // **The runs model scores with the same pair the coupled fit did** — each library's
        // clean rate and the sample's second class of site — and not with the marginal the
        // emitted summary carries. A single averaged rate would put back inside the runs
        // model the tail misspecification the second class exists to remove, and `F` is
        // read off a contrast between windows that the tail moves.
        let libraries = library_shares_over(self.read_group_histograms(), |ploidy| ploidy == at)
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
            });
        match site_noise {
            None => SampleLibraryNoise::new(libraries),
            Some(site_noise) => SampleLibraryNoise::with_site_noise(libraries, site_noise),
        }
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

    /// One locus that **both** libraries covered, five reads each.
    ///
    /// The shape every other fixture here lacks: elsewhere each site belongs to one read
    /// group, so the read-group table and the windowed one hold the same number of entries
    /// and a count taken from the wrong one is invisible.
    fn shared_site(
        start: u64,
        heterozygous: bool,
        groups: [ReadGroupId; 2],
    ) -> SampleLocusObservations {
        let mut observations = Vec::new();
        for group in groups {
            if heterozygous {
                observations.push(observation(b"C", group, 2));
                observations.push(observation(b"A", group, 3));
            } else {
                observations.push(observation(b"A", group, 5));
            }
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

    /// **A site two libraries covered is one reference position, not two.**
    ///
    /// The read-group table enters that site **twice** — once per group, at each group's own
    /// depth — and the windowed table once, at their combined depth. So a warrant summed
    /// over the wrong table doubles on exactly this shape and on no other fixture in this
    /// file: everywhere else a site belongs to one group, and the two tables agree.
    /// Measured: 20,002 against the 10,001 asserted below.
    #[test]
    fn a_site_two_libraries_covered_counts_one_position_and_not_two() {
        let groups = [ReadGroupId(1), ReadGroupId(2)];
        let mut config = config(supplied(0.3));
        config.read_groups = groups.to_vec();
        let loci: Vec<SampleLocusObservations> = (0..=MIN_SITES_TO_FIT)
            .map(|index| shared_site(index + 1, index.is_multiple_of(20), groups))
            .collect();

        let parameters = estimate_generic_parameters(loci.into_iter().map(Ok), &config)
            .expect("a two-library sample above the site floor is fittable");

        assert_eq!(
            parameters.inbreeding.expect("supplied").observations,
            MIN_SITES_TO_FIT + 1,
            "each site covers one reference position, however many libraries covered it"
        );
    }

    /// **The `Fitted` arm is entered at all**, which no other test here does — every other
    /// one supplies `F` to keep the runs model out of the way, and that left the whole arm
    /// dead: replacing it with `unreachable!()` left the suite green, as did suppressing the
    /// runs model with a bare `Ok((None, None))`.
    ///
    /// The fixture's 10,001 one-base sites all fall in one 100 kb window, against a floor of
    /// 3,000, so the runs model refuses — which is the point. Reaching that refusal proves
    /// the arm runs, the ploidy-2 lookup finds its entry and `library_noise` was built. **It
    /// says nothing about whether the runs model gets the right numbers**; that needs 300 Mb
    /// of windows and is F2's question.
    #[test]
    fn a_fitted_coefficient_is_attempted_rather_than_reported_absent() {
        let config = config(InbreedingMode::Fitted);

        let error = estimate_generic_parameters(
            a_samples_loci(ReadGroupId(1)).into_iter().map(Ok),
            &config,
        )
        .expect_err("one window of evidence is far below the runs model's floor");

        assert!(
            matches!(
                error,
                ParameterEstimationError::InbreedingNotFittable { windows: 1, .. }
            ),
            "the wrong failure, so the Fitted arm may not have run: {error}"
        );
    }

    /// **How the coupled fit terminated is reported rather than invented.** It is one of five
    /// fields on the emitted parameters and the equality test above cannot see it: both sides
    /// of that comparison would carry the same constant. A fit that ran out of iterations
    /// reported as converged is the failure this closes.
    #[test]
    fn the_coupled_fits_termination_is_the_one_it_reached() {
        let config = config(supplied(0.2));

        let parameters = estimate_generic_parameters(
            a_samples_loci(ReadGroupId(1)).into_iter().map(Ok),
            &config,
        )
        .expect("fittable");

        assert!(parameters.coupled_fit.converged);
        assert_eq!(
            parameters.coupled_fit.iterations, 2,
            "a single-library sample reaches its fixed point in two rounds — one to settle \
             and one to see nothing move"
        );
    }

    /// **Two shards merged answer as one walk did** — the sharded path is the whole reason
    /// the second entry point exists and nothing had exercised it.
    ///
    /// It also puts `accumulators()`'s stated contract under test: every shard must come from
    /// the config so that all of them share one `edges` object, which `merge` proves by
    /// pointer identity. Handing each call a fresh `Arc` instead leaves every other test in
    /// this file green and panics here.
    #[test]
    fn two_shards_merged_answer_as_one_walk() {
        let group = ReadGroupId(1);
        let config = config(supplied(0.2));
        let loci = a_samples_loci(group);
        let seam = loci.len() / 2;

        let mut first = config.accumulators();
        for locus in &loci[..seam] {
            first.add_locus(locus);
        }
        let mut second = config.accumulators();
        for locus in &loci[seam..] {
            second.add_locus(locus);
        }
        first.merge(second);

        let sharded = first.estimate(&config).expect("fittable");
        let single =
            estimate_generic_parameters(loci.into_iter().map(Ok), &config).expect("fittable");

        assert_eq!(sharded, single);
    }

    /// **The runs model is handed the clean rates *and* the second class, not the
    /// marginal** — and nothing else in the suite could see it if it were not.
    ///
    /// `F` is read off a contrast between windows, and the tail of the error distribution
    /// is what moves that contrast; a single averaged rate would put back inside the runs
    /// model the misspecification the second class exists to remove. The direct assertion
    /// is here because the alternative is a fixture of three thousand windows carrying a
    /// noisy site population, and this states the property the fixture would be built to
    /// state. Measured: replacing the pair with a plain rate leaves all 329 tests green.
    #[test]
    fn the_runs_model_is_given_the_second_class_of_site_and_not_the_marginal() {
        let group = ReadGroupId(1);
        let config = config(supplied(0.0));
        let mut accumulators = config.accumulators();
        for locus in a_samples_loci(group) {
            accumulators.add_locus(&locus);
        }

        let clean = ErrorRate::try_new(1.895e-3).expect("a probability");
        let mut settled = BTreeMap::new();
        settled.insert(
            group,
            Estimate {
                value: clean,
                provenance: Provenance::FittedHere,
                observations: 1_000,
            },
        );
        let site_noise =
            SiteNoise::try_new(0.0088, ErrorRate::try_new(5.29e-2).expect("a probability"))
                .expect("a share and a rate");
        let diploid = Ploidy::try_new(2).expect("a positive copy number");

        let with = accumulators.library_noise(&settled, diploid, Some(site_noise));
        assert_eq!(
            with.site_noise(),
            Some(site_noise),
            "the runs model was not given the sample's second class of site"
        );
        assert!(
            (with.libraries()[0].error_rate.get() - clean.get()).abs() < 1e-15,
            "the runs model was given {:e} rather than the clean rate {:e} — the marginal \
             has been folded in a step too early",
            with.libraries()[0].error_rate.get(),
            clean.get()
        );

        let without = accumulators.library_noise(&settled, diploid, None);
        assert_eq!(without.site_noise(), None);
    }

    /// One diploid site at `depth`, `alt_reads` of them showing an alternative base.
    fn site_at_depth(
        start: u64,
        depth: u32,
        alt_reads: u32,
        group: ReadGroupId,
    ) -> SampleLocusObservations {
        let mut observations = Vec::new();
        if alt_reads > 0 {
            observations.push(observation(b"C", group, alt_reads));
        }
        if depth > alt_reads {
            observations.push(observation(b"A", group, depth - alt_reads));
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

    /// A sample whose sites really do come in two populations: a clean body, a heterozygous
    /// class, and **a tail of sites disagreeing with the reference far more than any single
    /// per-base rate explains** — the shape measured on HG002 and the whole reason the
    /// second class exists.
    ///
    /// **Two earlier versions of this fixture bought the second class nothing, and why is
    /// worth keeping.** The first gave every noisy site exactly two alternative reads of
    /// ten, leaving three occupied cells — fitted exactly by three genotype frequencies and
    /// a free rate. The second spread them but left a *gap* at two alternative reads, and a
    /// gap is **under**-dispersion: a mixture of binomials can only ever be more spread out
    /// than a binomial, never less, so no second class could help and the fit correctly
    /// declined one, gaining 0.044 nats against a floor of 3.
    ///
    /// Thirty reads a site, so a heterozygote sits near fifteen and the tail at three to six
    /// cannot be mistaken for one. A rate explaining the 900 sites with one alternative read
    /// predicts about two sites with three, against the 30 here.
    fn loci_with_a_noisy_population(group: ReadGroupId) -> Vec<SampleLocusObservations> {
        const DEPTH: u32 = 30;
        let mut loci = Vec::new();
        let mut at = 1u64;
        for (count, alt_reads) in [
            (9_000, 0),
            (900, 1),
            (60, 2),
            (30, 3),
            (15, 4),
            (8, 5),
            (100, 15),
        ] {
            for _ in 0..count {
                loci.push(site_at_depth(at, DEPTH, alt_reads, group));
                at += 1;
            }
        }
        assert!(loci.len() as u64 > MIN_SITES_TO_FIT, "above the site floor");
        loci
    }

    /// **End to end: a sample with two populations of site is fitted with two, and the rate
    /// it emits is the marginal of them.**
    ///
    /// The discriminator is that **a marginal is not a rung**. Every rate the ladder can
    /// fit is one of its 161 rungs; a share-weighted average of two rungs is not, except by
    /// coincidence. So an emitted rate that lands exactly on a rung means the marginal was
    /// never applied — or was applied somewhere else and the clean rate emitted in its
    /// place, which is the swap this catches and which every other test in the module
    /// leaves alive.
    #[test]
    fn a_sample_with_two_populations_of_site_emits_the_marginal_rate() {
        let group = ReadGroupId(1);
        let config = config(supplied(0.0));
        let parameters = estimate_generic_parameters(
            loci_with_a_noisy_population(group).into_iter().map(Ok),
            &config,
        )
        .expect("a sample above the site floor is fittable");

        let site_noise = parameters
            .site_noise
            .expect("a tail no single per-base rate can reach is a second class");
        assert!(
            site_noise.noisy_error_rate().get()
                > site_noise
                    .marginal_error_rate(ErrorRate::try_new(1e-4).expect("a probability"))
                    .get(),
            "the noisy class must be the worse-behaved of the two"
        );

        let emitted = parameters.error_rate[&group].value.get();
        let ladder = error_rate_ladder();
        let nearest = ladder
            .iter()
            .map(|rung| (rung.get() - emitted).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(
            nearest > 1e-12,
            "the emitted rate {emitted:e} sits exactly on a ladder rung, so it is a fitted \
             clean rate rather than the marginal of the sample's two classes of site"
        );
    }
}
