//! **Opening a cohort of census files**, and refusing one that cannot be fitted as a cohort.
//!
//! A census holds what one sample showed at the loci a run chose to keep. A fit needs every
//! sample's at once, and needs them to be comparable: recorded under the same selection, in the
//! same units, and built from the psps they claim to be built from. This module opens the files,
//! makes those checks, and hands back the object the estimator reads
//! ([`CohortCensusEvidence`]).
//!
//! # Three ways a cohort is not one
//!
//! - **The samples recorded different things.** Twelve terms say which loci were asked for, which
//!   came back, and in what units — a selection made under another seed, a different depth
//!   ladder, a different cap. Two samples that disagree on any of them hold rows that mean
//!   different things and every one of them fails silently, so the cohort is refused at the door
//!   rather than fitted.
//! - **Two samples declare one `@RG ID`.** The identifier is unique across a whole run, so a pair
//!   that shares one was not produced by one run, and pooling their libraries into a single
//!   sequencing-error rate is damage nothing downstream would report.
//! - **A census is not the psp beside it any more.** A census names the psp it was built from;
//!   evidence from other reads is otherwise indistinguishable from this run's.
//!
//! # What the psp is needed for, and what it is not
//!
//! **The psps are not read.** What is taken from each is the digest of its header — one short
//! read — and compared against the digest the census carries. A census whose psp is not beside it
//! is refused rather than trusted, because nothing can then check it and nothing can rebuild it
//! ([`freshness_by_header`]). That function's own documentation says precisely what a
//! header-only check leaves out.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::ng::parameter_estimation::joint::census::{
    CensusError, CohortCensusEvidence, CohortRefusal,
};
use crate::ng::parameter_estimation::joint::census_file::{
    Freshness, freshness_by_header, open_census,
};
use crate::ng::psp::{self, PspReadError};
use crate::ng::read::input::read_groups::{
    NameOrigin, NameWithOrigin, ReadGroup, ReadGroups, SampleReadGroups,
};
use crate::pop_var_caller_exp::generate_census::CENSUS_FILE_EXTENSION;
use crate::pop_var_caller_exp::generate_psps::PSP_FILE_EXTENSION;

/// A cohort of censuses, opened and agreed with.
#[derive(Debug)]
#[must_use]
pub struct OpenCensusCohort {
    /// What the fit reads — every sample's census, on run-wide read-group identifiers.
    pub evidence: CohortCensusEvidence,
    /// One entry a sample, in the order the censuses were named.
    pub samples: Vec<CensusInCohort>,
}

/// One census of a cohort, and where it came from.
#[derive(Debug, Clone)]
pub struct CensusInCohort {
    /// The individual, as its census names it.
    pub sample: String,
    /// The census file.
    pub census: PathBuf,
    /// The psp it was checked against.
    pub psp: PathBuf,
}

/// Why a cohort of censuses could not be opened.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CensusCohortError {
    /// No census was named.
    #[error("a fit needs at least one census")]
    NoCensuses,

    /// A census file could not be opened or is not one.
    #[error("the census {} could not be read", path.display())]
    CensusNotRead {
        /// The file.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: Box<CensusError>,
    },

    /// The psp a census names is not beside it.
    #[error(
        "{}'s census names a psp and there is none at {}; a census cannot be checked without \
         the psp it was built from, and cannot be rebuilt without it either",
        sample,
        psp.display()
    )]
    PspNotBesideIt {
        /// The individual.
        sample: String,
        /// Where its psp was looked for.
        psp: PathBuf,
    },

    /// The psp is there and could not be read far enough to be identified.
    #[error("the psp {} could not be identified", path.display())]
    PspNotIdentified {
        /// The psp.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: Box<PspReadError>,
    },

    /// The census was built from a different psp than the one beside it.
    #[error(
        "{sample}'s census was built from another psp — {field} differs from the one at {}. \
         Rebuild it with generate-census, or put the psp it was built from beside it",
        psp.display()
    )]
    BuiltFromAnotherPsp {
        /// The individual.
        sample: String,
        /// The psp beside it now.
        psp: PathBuf,
        /// What differs, in the words the census file's own check uses.
        field: &'static str,
    },

    /// The census names no psp at all, so nothing can check it.
    #[error(
        "{sample}'s census names no psp, so nothing can say whether it was built from these \
         reads. Rebuild it with generate-census"
    )]
    NamesNoPsp {
        /// The individual.
        sample: String,
    },

    /// The samples cannot be fitted as one cohort.
    #[error("these censuses are not one cohort")]
    NotOneCohort {
        /// Which of the three disagreements it is.
        #[source]
        source: Box<CohortRefusal>,
    },
}

/// **Open every census named, check each against its psp, and assemble the cohort.**
///
/// The psp is looked for beside the census, under the same stem: `<sample>.census` is checked
/// against `<sample>.psp` in the same directory. That is where `generate-psps` leaves the pair
/// and where `generate-census --output-dir` puts them when it is pointed at the psps' own
/// directory.
///
/// # Errors
///
/// See [`CensusCohortError`]. **Every file is opened and checked before the cohort is
/// assembled**, so a run over sixty samples reports the one that is stale rather than the first
/// disagreement the merge happens to reach.
pub fn open_census_cohort(paths: &[PathBuf]) -> Result<OpenCensusCohort, CensusCohortError> {
    if paths.is_empty() {
        return Err(CensusCohortError::NoCensuses);
    }

    let mut opened = Vec::with_capacity(paths.len());
    let mut samples = Vec::with_capacity(paths.len());
    for path in paths {
        let (evidence, named) =
            open_census(path).map_err(|source| CensusCohortError::CensusNotRead {
                path: path.clone(),
                source: Box::new(source),
            })?;
        let sample = evidence.sample.clone();
        let psp = psp_beside(path);

        // **The digest, not the psp.** One short read of the file's own header bytes, which is
        // what a census names its psp by.
        let header_in_hand = match psp.exists() {
            true => Some(psp::header_digest(&psp).map_err(|source| {
                CensusCohortError::PspNotIdentified {
                    path: psp.clone(),
                    source: Box::new(source),
                }
            })?),
            false => None,
        };

        match freshness_by_header(named, header_in_hand) {
            Freshness::Fresh => {}
            Freshness::Rebuild(field) if header_in_hand.is_some() => {
                return Err(CensusCohortError::BuiltFromAnotherPsp { sample, psp, field });
            }
            Freshness::Rebuild(_) => return Err(CensusCohortError::NamesNoPsp { sample }),
            Freshness::Refused(_) if psp.exists() => {
                return Err(CensusCohortError::NamesNoPsp { sample });
            }
            Freshness::Refused(_) => {
                return Err(CensusCohortError::PspNotBesideIt { sample, psp });
            }
        }

        samples.push(CensusInCohort {
            sample,
            census: path.clone(),
            psp,
        });
        opened.push(evidence);
    }

    let evidence =
        CohortCensusEvidence::new(opened).map_err(|source| CensusCohortError::NotOneCohort {
            source: Box::new(source),
        })?;
    Ok(OpenCensusCohort { evidence, samples })
}

/// The psp a census is checked against: the same stem, in the same directory.
fn psp_beside(census: &Path) -> PathBuf {
    let mut psp = census.to_path_buf();
    // `set_extension` replaces the last one, which is what turns `zeta.census` into `zeta.psp`
    // and leaves a stem holding dots alone.
    if psp
        .extension()
        .is_some_and(|it| it == CENSUS_FILE_EXTENSION)
    {
        psp.set_extension(PSP_FILE_EXTENSION);
    } else {
        psp.as_mut_os_string()
            .push(format!(".{PSP_FILE_EXTENSION}"));
    }
    psp
}

/// **The run's read-group table, built from what the censuses declare.**
///
/// A parameters file names a read group by its `@RG ID`, its library and its sample, and the
/// table is where those three meet. The censuses carry all three, so no alignment file and no psp
/// is opened to build it — which is the whole reason the names went into the census (step C1).
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
pub fn read_groups_of(cohort: &CohortCensusEvidence, censuses: &[CensusInCohort]) -> ReadGroups {
    let mut groups = Vec::new();
    let mut per_sample = Vec::with_capacity(cohort.len());
    for (index, sample) in cohort.samples().iter().enumerate() {
        // Only for a message to be able to name a file; nothing keys on it.
        let file: Arc<Path> = censuses.get(index).map_or_else(
            || Arc::from(Path::new("")),
            |it| Arc::from(it.census.as_path()),
        );
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
    //! **Plan step C3**: a cohort of censuses opens, and every way it can fail to be one is
    //! provoked and named.

    use super::*;
    use crate::pop_var_caller_exp::generate_psps::{
        GeneratePspsArgs, census_path_for, psp_path_for, run_generate_psps,
    };
    use crate::pop_var_caller_exp::test_fixtures::a_cohort_on_disk;

    /// The fixture cohort walked into psps, with its censuses beside them — the pair a fit is
    /// meant to be handed.
    fn a_cohort_with_its_censuses() -> (
        crate::pop_var_caller_exp::test_fixtures::ACohortOnDisk,
        PathBuf,
        Vec<PathBuf>,
    ) {
        use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
        use crate::ng::region_typing::segment_criteria::{
            DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
        };

        let cohort = a_cohort_on_disk();
        let psps = cohort.directory.path().join("psps");
        run_generate_psps(&GeneratePspsArgs {
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
        })
        .expect("the cohort walks into psps");
        let censuses = vec![
            census_path_for(&psps, "alpha"),
            census_path_for(&psps, "zeta"),
        ];
        (cohort, psps, censuses)
    }

    /// **The pair a run is meant to be handed opens**, and the cohort's read groups are as many
    /// as the samples' libraries rather than as many as one sample's.
    #[test]
    fn a_cohort_of_censuses_beside_their_psps_opens() {
        let (_cohort, psps, censuses) = a_cohort_with_its_censuses();

        let open = open_census_cohort(&censuses).expect("the censuses are this cohort's");

        assert_eq!(open.samples.len(), 2);
        assert_eq!(open.samples[0].sample, "alpha");
        assert_eq!(open.samples[1].sample, "zeta");
        assert_eq!(open.samples[1].psp, psp_path_for(&psps, "zeta"));
        assert_eq!(
            open.evidence.read_groups().len(),
            2,
            "two samples' libraries, under identifiers of their own",
        );
    }

    /// **A census with no psp beside it is refused**, because nothing can then check it and
    /// nothing can rebuild it.
    #[test]
    fn a_census_whose_psp_is_not_beside_it_is_refused() {
        let (cohort, psps, _censuses) = a_cohort_with_its_censuses();
        // The census is moved somewhere with no psp, which is what naming an --output-dir of its
        // own leaves a person with.
        let alone = cohort.directory.path().join("alone");
        std::fs::create_dir_all(&alone).expect("the scratch dir is ours");
        let moved = alone.join("zeta.census");
        std::fs::copy(census_path_for(&psps, "zeta"), &moved).expect("a copy");

        let error = open_census_cohort(&[moved]).expect_err("nothing can check it");

        assert!(
            matches!(&error, CensusCohortError::PspNotBesideIt { sample, .. } if sample == "zeta"),
            "{error:?}",
        );
    }

    /// **A census built from another psp is refused, and the refusal says to rebuild it.**
    ///
    /// The pair is put wrong the way a person puts it wrong: one sample's census is copied in
    /// beside another sample's psp. **Deliberately not "walk the sample twice"** — two walks of
    /// one sample a second apart write the same header, because the only field they legitimately
    /// differ in is a timestamp at one-second resolution, so that test would pass or fail on how
    /// fast the machine is.
    #[test]
    fn a_census_built_from_another_psp_is_refused() {
        let (cohort, psps, _censuses) = a_cohort_with_its_censuses();
        let wrong = cohort.directory.path().join("wrong");
        std::fs::create_dir_all(&wrong).expect("the scratch dir is ours");
        // zeta's census, put where alpha's psp is — so the census names a psp that is not the
        // one beside it.
        std::fs::copy(census_path_for(&psps, "zeta"), wrong.join("alpha.census")).expect("a copy");
        std::fs::copy(psp_path_for(&psps, "alpha"), wrong.join("alpha.psp")).expect("a copy");

        let error = open_census_cohort(&[wrong.join("alpha.census")])
            .expect_err("this census was built from zeta's psp");

        let CensusCohortError::BuiltFromAnotherPsp { sample, field, .. } = &error else {
            panic!("a census naming another psp is its own refusal: {error:?}");
        };
        assert_eq!(sample, "zeta", "the census says which sample it holds");
        assert_eq!(*field, "the pileup's header");
        assert!(
            error.to_string().contains("generate-census"),
            "the refusal says what to do about it: {error}",
        );
    }

    /// **A run with no census named is refused rather than fitted over nothing.**
    #[test]
    fn a_cohort_of_no_censuses_is_refused() {
        let error = open_census_cohort(&[]).expect_err("a fit needs evidence");
        assert!(matches!(error, CensusCohortError::NoCensuses), "{error:?}");
    }

    /// **A file that is not a census is refused by name.**
    #[test]
    fn a_file_that_is_not_a_census_is_refused() {
        let (cohort, _psps, _censuses) = a_cohort_with_its_censuses();
        let not_a_census = cohort.directory.path().join("prose.census");
        std::fs::write(&not_a_census, b"this is not a census").expect("the scratch dir is ours");

        let error = open_census_cohort(std::slice::from_ref(&not_a_census))
            .expect_err("prose is not a census");

        assert!(
            matches!(&error, CensusCohortError::CensusNotRead { path, .. } if path == &not_a_census),
            "{error:?}",
        );
    }

    /// The psp a census is checked against is its own stem with the psp extension.
    #[test]
    fn the_psp_a_census_is_checked_against_sits_beside_it() {
        assert_eq!(
            psp_beside(Path::new("/where/zeta.census")),
            PathBuf::from("/where/zeta.psp"),
        );
        assert_eq!(
            psp_beside(Path::new("/where/SRR7279481.p1.census")),
            PathBuf::from("/where/SRR7279481.p1.psp"),
            "a sample name holding a dot keeps it",
        );
        assert_eq!(
            psp_beside(Path::new("/where/zeta")),
            PathBuf::from("/where/zeta.psp"),
            "a census named without the extension still names its psp",
        );
    }
}
