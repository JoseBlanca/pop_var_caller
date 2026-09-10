//! **The censuses of a cohort of psps, read out of the files that carry them**, and the refusal
//! for a set that cannot be fitted as one cohort.
//!
//! A census holds what one sample showed at the loci a run chose to keep, and it lives in its
//! psp's trailer (`psp_census_pair.md` §3). A fit needs every sample's at once, and needs them to
//! be comparable: recorded under the same selection, in the same units. This module reads them
//! and hands back the object the estimator reads ([`CohortCensusEvidence`]).
//!
//! # What is checked here, and what was checked before it
//!
//! **Three ways a set of samples is not one cohort**, and all three are
//! [`CohortCensusEvidence::new`]'s:
//!
//! - **The samples recorded different things.** Twelve terms say which loci were asked for, which
//!   came back, and in what units — a selection made under another seed, a different depth
//!   ladder, a different cap. Two samples that disagree on any of them hold rows that mean
//!   different things and every one of them fails silently, so the cohort is refused at the door
//!   rather than fitted.
//! - **Two samples declare one `@RG ID`.** The identifier is unique across a whole run, so a pair
//!   that shares one was not produced by one run, and pooling their libraries into a single
//!   sequencing-error rate is damage nothing downstream would report.
//! - **A sample holds evidence filed under a read group it does not declare.** Renumbering onto
//!   run-wide identifiers has nowhere to put that section, and the alternative — dropping it —
//!   has no symptom.
//!
//! **A whole class of refusal went when the census moved into the psp.** A census file could be
//! copied away from its psp, or beside another sample's, so it carried the digest of the psp it
//! was built from and the opener compared the two — *this census was built from another psp*,
//! *this census names a psp and there is none beside it*. A trailer cannot be separated from the
//! file it is the tail of, so the identity is structural and those refusals have no subject
//! (spec §3). What replaces them is the freshness judgement:
//! [`census_freshness`](super::census_freshness) says whether the trailer holds a census this
//! build reads at all.
//!
//! **The cohort's own agreement — its ground, its catalog, its criteria — is checked before this
//! runs**, when the psps are opened ([`OpenPspCohort::open`]).
//!
//! **The psps' records are never read.** What is taken from each file is the head of its trailer —
//! at most a megabyte, and less for a census shorter than that — out of which the census's header
//! and its directory are decoded, a few hundred bytes of them; every section after those is read
//! only when something asks for it (spec §5).

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::ng::parameter_estimation::joint::census::{
    ByteExtent, CensusError, CohortCensusEvidence, CohortRefusal, SampleCensusEvidence,
};
use crate::ng::parameter_estimation::joint::census_file::open_census_within;
use crate::ng::psp::PspReader;
use crate::ng::read::input::read_groups::{
    NameOrigin, NameWithOrigin, ReadGroup, ReadGroups, SampleReadGroups,
};
use crate::ng::run::psp_caller::OpenPspCohort;

/// Why a cohort's censuses could not be read as one.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CensusCohortError {
    /// A psp's trailer could not be read as a census.
    ///
    /// **A psp with no census in it reaches this too**, as a malformed one, and that is a message
    /// about bytes where the user needs a message about regenerating. Plan step C3 takes the
    /// verdicts of [`census_freshness`](super::census_freshness) first, so a stale cohort is
    /// named as stale before anything is decoded.
    #[error("the census in {} could not be read", path.display())]
    CensusNotRead {
        /// The psp carrying it.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: Box<CensusError>,
    },

    /// The samples cannot be fitted as one cohort.
    #[error("these censuses are not one cohort")]
    NotOneCohort {
        /// Which of the three disagreements it is.
        #[source]
        source: Box<CohortRefusal>,
    },
}

/// **Every census of an opened cohort, read out of the psps that carry them.**
///
/// One sample a psp, in the cohort's own order, which is the order the run's samples are in.
///
/// **Lazily**, which is the whole reason the fit reads a census by path and offset rather than
/// taking the psp's trailer whole: what is decoded here is each census's header and its directory,
/// and a section only when something asks for it (spec §5). What that costs off the disk is the
/// head read the census reader does — at most a megabyte a sample, transient, one buffer at a time
/// — against a thousand samples' trailers taken whole, which would be tens of gigabytes held.
///
/// **The two halves composed, and what production calls is the halves.**
/// `estimate-parameters` reads the censuses, judges each against its own run's settings, and
/// assembles afterwards ([`each_census_in_the_cohorts_psps`], [`the_censuses_as_one_cohort`]);
/// this is the convenience a test or a caller with nothing to judge wants. **So it has no caller
/// in a shipped command, and that is on purpose rather than an oversight** — it is the whole-cohort
/// read the tests of that assembly are written against, and deleting it would leave each of them
/// composing the two halves by hand.
///
/// # Errors
///
/// [`CensusCohortError::CensusNotRead`] naming the psp whose trailer will not decode, and
/// [`CensusCohortError::NotOneCohort`] for samples that cannot be fitted together.
pub fn every_census_in_the_cohorts_psps(
    cohort: &OpenPspCohort,
) -> Result<CohortCensusEvidence, CensusCohortError> {
    the_censuses_as_one_cohort(each_census_in_the_cohorts_psps(cohort)?)
}

/// **Each census of an opened cohort, read out of its psp and not yet assembled into a cohort** —
/// one sample a psp, in the cohort's own order, read exactly as
/// [`every_census_in_the_cohorts_psps`] reads them.
///
/// **Apart from the assembly because a command may have to judge each census before refusing the
/// cohort.** Assembling refuses two samples that recorded different settings, naming the pair
/// and nothing to do; a command that compares each one with the settings its own run records
/// under can instead name every stale sample and the fix (plan step C5), and assembles
/// afterwards with [`the_censuses_as_one_cohort`].
///
/// # Errors
///
/// [`CensusCohortError::CensusNotRead`] naming the psp whose trailer will not decode.
pub fn each_census_in_the_cohorts_psps(
    cohort: &OpenPspCohort,
) -> Result<Vec<SampleCensusEvidence>, CensusCohortError> {
    let mut opened = Vec::with_capacity(cohort.sample_count());
    for (path, psp) in cohort.each_psp_with_its_path_read_only() {
        opened.push(the_census_in_a_psp(path, psp)?);
    }
    Ok(opened)
}

/// **One psp's census, out of its own trailer** — its header and its directory decoded, and not
/// one section touched.
///
/// **What a command reads when it may not need every sample's census.** `regenerate-census` asks
/// this of the psps whose heads look fresh, to compare what each census recorded with what its own
/// run records under (spec §8's third cause). What it costs is the census reader's head read — at
/// most a mebibyte, and less for a census shorter than that — out of which a few hundred bytes are
/// decoded; **what it does not cost is the psp's records**, which is the pass that command exists
/// to avoid (plan step D2).
///
/// # Errors
///
/// [`CensusCohortError::CensusNotRead`] naming the psp whose trailer will not decode — which
/// includes a psp carrying no census at all, so a caller that means to tell those apart judges the
/// head first ([`census_freshness`](super::census_freshness)).
pub fn the_census_in_a_psp(
    path: &Path,
    psp: &PspReader,
) -> Result<SampleCensusEvidence, CensusCohortError> {
    // **The footer says where the census is**, and it was read when the psp was opened, so
    // finding it costs nothing.
    let footer = psp.footer();
    let evidence = open_census_within(
        path,
        ByteExtent::new(footer.trailer_offset, footer.trailer_bytes),
    )
    .map_err(|source| CensusCohortError::CensusNotRead {
        path: path.to_path_buf(),
        source: Box::new(source),
    })?;
    Ok(evidence)
}

/// **A cohort's censuses assembled into one**, refusing samples that cannot be fitted together.
///
/// # Errors
///
/// [`CensusCohortError::NotOneCohort`] for samples that recorded different settings or that claim
/// one read group.
pub fn the_censuses_as_one_cohort(
    censuses: Vec<SampleCensusEvidence>,
) -> Result<CohortCensusEvidence, CensusCohortError> {
    CohortCensusEvidence::new(censuses).map_err(|source| CensusCohortError::NotOneCohort {
        source: Box::new(source),
    })
}

/// **The run's read-group table, built from what the censuses declare.**
///
/// A parameters file names a read group by its `@RG ID`, its library and its sample, and the
/// table is where those three meet. The censuses carry all three, so this reads no file at all —
/// which is the whole reason the read-group names went into the census
/// (`parameter_prepass_runs.md` step C1).
///
/// **The identifiers are the cohort's, not each census's own.** Assembling the cohort renumbered
/// every sample onto run-wide identifiers; this walks the samples in that order and mints the
/// same numbering, so the table and the evidence agree by construction.
///
/// **The library's origin is recorded as synthesized**, which is the weaker of the two claims and
/// the one that cannot be false: a census records the library the walk *resolved* — `@RG LB`, or
/// the name the walk invented where the file declared none — and not which of the two it was.
/// The calling stage's own table over stored files says the same thing for the same reason.
#[must_use]
pub fn read_groups_of(cohort: &CohortCensusEvidence, psps: &[PathBuf]) -> ReadGroups {
    let mut groups = Vec::new();
    let mut per_sample = Vec::with_capacity(cohort.len());
    for (index, sample) in cohort.samples().iter().enumerate() {
        // Only for a message to be able to name a file; nothing keys on it.
        let file: Arc<Path> = psps
            .get(index)
            .map_or_else(|| Arc::from(Path::new("")), |it| Arc::from(it.as_path()));
        let mut mine = Vec::new();
        for (id, named) in sample.declared_read_groups() {
            groups.push(ReadGroup {
                file: Arc::clone(&file),
                id: named.declared_id.clone().into_boxed_str(),
                sample: sample.sample.clone().into_boxed_str(),
                library: NameWithOrigin {
                    value: named.library.clone().into_boxed_str(),
                    origin: NameOrigin::Synthesized,
                },
                // The experiment is the library copied, because nothing reads an experiment tag
                // yet — direct mode's own rule.
                experiment: NameWithOrigin {
                    value: named.library.clone().into_boxed_str(),
                    origin: NameOrigin::Synthesized,
                },
                platform: None,
            });
            mine.push(*id);
        }
        per_sample.push(SampleReadGroups {
            sample: sample.sample.clone().into_boxed_str(),
            read_groups: mine,
        });
    }
    ReadGroups::of_merged_tables(groups, per_sample)
}

#[cfg(test)]
mod tests {
    //! A cohort's censuses are read out of its psps, and a set that is not one cohort is refused.
    //!
    //! **Five tests went with the census file** (plan step C2). Three provoked ways a census file
    //! and a psp could be put wrong beside each other — a census with no psp beside it, a census
    //! built from another sample's psp, and the stem rule that paired the two — and a trailer can
    //! be none of those. The other two changed subject rather than disappearing: a cohort of no
    //! files is now the psp opener's refusal, and *a file that is not a census* is now *a psp
    //! whose trailer is not one*, below. What a psp carrying no census is **told** is plan step
    //! C3's report, over [`census_freshness`](super::super::census_freshness)'s verdicts.

    use super::*;
    use crate::pop_var_caller_exp::generate_psps::{
        GeneratePspsArgs, psp_path_for, run_generate_psps,
    };
    use crate::pop_var_caller_exp::test_fixtures::a_cohort_on_disk;

    /// The fixture cohort walked into psps, opened as the run opens it.
    fn a_walked_cohort() -> (
        crate::pop_var_caller_exp::test_fixtures::ACohortOnDisk,
        Vec<PathBuf>,
        OpenPspCohort,
    ) {
        use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
        use crate::ng::region_typing::segment_criteria::{
            DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
        };

        let cohort = a_cohort_on_disk();
        let psps = cohort.directory.path().join("psps");
        let walk = GeneratePspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            alignments: cohort.alignments.clone(),
            output_dir: psps.clone(),
            regions: None,
            force: false,
            build_index_if_missing: false,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        };
        run_generate_psps(&walk).expect("the cohort walks into psps");
        let paths = vec![psp_path_for(&psps, "alpha"), psp_path_for(&psps, "zeta")];
        let open = OpenPspCohort::open(&paths).expect("the walk wrote one cohort");
        (cohort, paths, open)
    }

    /// **A walked cohort's censuses read out of its psps**, and the cohort's read groups are as
    /// many as the samples' libraries rather than as many as one sample's.
    #[test]
    fn the_censuses_of_a_walked_cohort_are_read_from_its_psps() {
        let (_cohort, _paths, open) = a_walked_cohort();

        let evidence = every_census_in_the_cohorts_psps(&open).expect("the walk wrote them");

        assert_eq!(evidence.len(), 2);
        assert_eq!(evidence.samples()[0].sample, "alpha");
        assert_eq!(evidence.samples()[1].sample, "zeta");
        assert_eq!(
            evidence.read_groups().len(),
            2,
            "two samples' libraries, under identifiers of their own",
        );
    }

    /// **Reading a cohort's censuses takes no psp's trailer whole, decodes no section, and
    /// leaves every section where it is** — the property that lets a cohort too large to hold be
    /// opened at all (spec §5).
    ///
    /// **Three assertions, because the two counters see different halves and neither sees the
    /// third.** `trailer_bytes_read` counts what comes out of a psp through its own reader, so a
    /// version of this that took each trailer whole and decoded it resident moves it off zero.
    /// The census's `bytes_read` counts *section* reads only — not the head read that opening
    /// does — so at this point it must still be zero. And neither of them can tell a lazy census
    /// from one already held in memory: what does is asking for a section afterwards and finding
    /// that it reached the disk.
    #[test]
    fn reading_a_cohorts_censuses_takes_no_trailer_whole_and_leaves_every_section_on_disk() {
        use crate::ng::parameter_estimation::joint::census_file::{bytes_read, reset_bytes_read};

        let (_cohort, paths, open) = a_walked_cohort();
        let carried: u64 = paths
            .iter()
            .map(|path| {
                let psp = crate::ng::psp::PspReader::open(path).expect("a walked psp opens");
                psp.footer().trailer_bytes
            })
            .sum();
        assert!(carried > 0, "the walk sealed a census into each psp");
        crate::ng::psp::reset_trailer_bytes_read();
        reset_bytes_read();

        let mut evidence = every_census_in_the_cohorts_psps(&open).expect("the walk wrote them");

        assert_eq!(evidence.len(), 2);
        assert_eq!(
            crate::ng::psp::trailer_bytes_read(),
            0,
            "no trailer was taken whole, and the two carry {carried} bytes of census between them",
        );
        assert_eq!(
            bytes_read(),
            0,
            "and no section of either census was decoded"
        );

        let groups = evidence.read_groups().to_vec();
        evidence
            .with_generic(&groups, |samples| samples.len())
            .expect("the generic sections are there to be read");

        assert!(
            bytes_read() > 0,
            "asking for a section reached the file, so the census was never held whole",
        );
    }

    /// **A psp whose trailer is not a census is refused, naming the psp.**
    ///
    /// The trailer is replaced with bytes that are not a census — which is what a psp written
    /// before the census moved into it looks like, except that its trailer is empty rather than
    /// prose, and both arrive here the same way. **This is the message plan step C3 puts a
    /// regeneration report in front of**; until then it is what a stale cohort says.
    #[test]
    fn a_psp_whose_trailer_is_not_a_census_is_refused_naming_the_psp() {
        let (_cohort, paths, _open) = a_walked_cohort();
        crate::ng::psp::replace_trailer(&paths[1], b"not a census").expect("the tail rewrites");
        let open = OpenPspCohort::open(&paths).expect("the psps still open");

        let error = every_census_in_the_cohorts_psps(&open).expect_err("that is not a census");

        assert!(
            matches!(&error, CensusCohortError::CensusNotRead { path, .. } if path == &paths[1]),
            "{error:?}",
        );
    }
}
