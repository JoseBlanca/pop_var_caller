//! **The second producer: a census built from a stored psp rather than from alignment files.**
//!
//! A census is a small set of loci at which every sample keeps its raw evidence instead of
//! folding it into a histogram, so that the fit can ask the same question of everybody
//! ([`parameter_prepass_census_sites.md`](../../../doc/devel/ng/spec/parameter_prepass_census_sites.md)).
//! One is written while the reads are walked ([`gatherer`](super::gatherer)); this module builds
//! the same object afterwards, from the psp that walk produced.
//!
//! # Why a second producer exists
//!
//! Three situations need one, and none of them can re-walk the reads
//! ([`parameter_prepass_joint_records.md`](../../../doc/devel/ng/spec/parameter_prepass_joint_records.md)
//! §6.1): psps written before censuses existed; a census lost, or built under settings since
//! changed; and a census wanted larger than the one on disk.
//!
//! # What makes the two agree
//!
//! **Both build their writer through [`CensusPlan::writer_for`](super::CensusPlan::writer_for)
//! and feed it the same type** — [`SampleLocusObservations`], which is what the walk yields and
//! what a psp decodes back to. So the two producers differ in where the loci come from and in
//! nothing else, which is what makes comparing their two files byte for byte a statement about
//! the psp: whether it carries everything a census needs (§7.12). A field the psp drops shows up
//! there, and a repeat tract's per-read length is the one most likely to.

use std::path::{Path, PathBuf};

use crate::ng::parameter_estimation::joint::census::{
    CensusError, DepthCode, SampleCensusEvidence,
};
use crate::ng::parameter_estimation::joint::census_file::PileupIdentity;
use crate::ng::psp::{self, PspReadError, PspReader};
use crate::ng::run::{CensusPlan, Segmentation};
use crate::ng::types::ReadGroupId;

/// One sample's census, built from its stored psp, and which psp it came from.
#[derive(Debug)]
#[must_use]
pub struct CensusOfStoredPileup {
    /// The evidence, ready to be written to a census file or handed to a fit.
    pub evidence: SampleCensusEvidence,
    /// **Which psp built it** — the digest of that file's header and how many records it holds.
    /// A fit compares this against the psp beside it and refuses a pair that has come apart.
    pub identity: PileupIdentity,
    /// The individual, as the psp's own header names it.
    pub sample: String,
    /// Every read group the psp declares, numbered as that psp numbers them.
    pub read_groups: Vec<ReadGroupId>,
}

/// **What one sample actually put into its census** — the counts a run reports about it.
///
/// **Two numbers a piece, because one alone cannot be read.** A sample with reads at 900 kept
/// positions has done well or badly depending on whether the selection kept a thousand of them
/// or two million, and a reader who is told only the first cannot tell which.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[must_use]
pub struct CensusTally {
    /// Kept ordinary positions this sample has at least one read at.
    pub positions_with_reads: u64,
    /// Kept ordinary positions the selection holds at all — **the same for every sample of the
    /// run**, since the selection is the run's and not the sample's.
    pub positions_kept: u64,
    /// Kept repeat tracts this sample has at least one read at.
    pub tracts_with_reads: u64,
    /// Kept repeat tracts the selection holds at all.
    pub tracts_kept: u64,
}

impl CensusTally {
    /// **Whether this sample's census would tell a fit anything at all.**
    ///
    /// A psp whose walk covered ground the selection kept nothing in, or whose reads reached
    /// none of what it did keep, produces a census that is entirely denominator. That is a
    /// legitimate outcome and not an error — but a run that omitted such a sample from its
    /// report would leave somebody hunting for a file that was written exactly as asked.
    #[must_use]
    pub fn contributes_nothing(self) -> bool {
        self.positions_with_reads == 0 && self.tracts_with_reads == 0
    }
}

/// Why a census could not be built from a stored psp.
///
/// **Both variants are about the file rather than about ng.** A psp that cannot be opened or
/// whose records stop part-way is data a run was handed, so neither is a panic and neither may
/// reach a caller as a half-built census.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CensusFromPspError {
    /// The file could not be opened, or its header could not be read.
    #[error("the psp {} could not be opened", path.display())]
    NotOpened {
        /// The psp.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: Box<PspReadError>,
    },
    /// The file opened and then failed part-way through its records, naming how far it got.
    #[error(
        "the psp {} failed after {records} records; a census cannot be built from part of a file",
        path.display()
    )]
    RecordsStopped {
        /// The psp.
        path: PathBuf,
        /// How many records were read before it failed.
        records: u64,
        /// What the reader said.
        #[source]
        source: Box<PspReadError>,
    },
}

impl CensusOfStoredPileup {
    /// Count what this sample put into its census.
    ///
    /// **Read off the census itself rather than counted while it was built**, so the numbers a
    /// run reports and the file it wrote cannot come to disagree. A position counts when any of
    /// the sample's read groups saw a read there; the depth ladder's bin 0 is zero reads, which
    /// is a walked position with nothing at it and not the same as a position never walked.
    ///
    /// # Errors
    ///
    /// [`CensusError`] when a section cannot be read, which for a census still resident in
    /// memory cannot happen.
    pub fn tally(&mut self) -> Result<CensusTally, CensusError> {
        let groups = self.evidence.read_groups();
        let strata = self.evidence.strata();

        let positions = self.evidence.with_generic(&groups, |lent| {
            let kept = lent.first().map_or(0, |evidence| evidence.depth().len());
            let mut with_reads = 0_u64;
            for index in 0..kept {
                let any = lent.iter().any(|evidence| {
                    matches!(evidence.depth().get(index), DepthCode::Binned(bin) if bin.get() > 0)
                });
                if any {
                    with_reads += 1;
                }
            }
            (with_reads, kept as u64)
        })?;

        let mut tracts_with_reads = 0_u64;
        let mut tracts_kept = 0_u64;
        for stratum in &strata {
            let counted =
                self.evidence
                    .with_strata(groups[0], std::slice::from_ref(stratum), |lent| {
                        lent.first().map_or(0, |section| section.len())
                    })?;
            tracts_kept += counted as u64;
            for locus in 0..counted {
                let mut reached = false;
                for group in &groups {
                    reached |= self.evidence.with_strata(
                        *group,
                        std::slice::from_ref(stratum),
                        |lent| {
                            lent.first()
                                .is_some_and(|section| section.offsets(locus).total() > 0)
                        },
                    )?;
                }
                if reached {
                    tracts_with_reads += 1;
                }
            }
        }

        Ok(CensusTally {
            positions_with_reads: positions.0,
            positions_kept: positions.1,
            tracts_with_reads,
            tracts_kept,
        })
    }
}

/// **Build one sample's census from its stored psp.**
///
/// `plan` is the run's selection — which positions and which tracts every sample of the cohort
/// keeps — and `segmentation` is the ground the psp's walk covered. Both must be the ones the
/// psp was written under: a selection made under other terms produces a census that cannot be
/// pooled with the cohort's others, and the recording terms travelling in the file are what say
/// so at the fit.
///
/// # What it reads, and what it does not
///
/// It reads every record of the file once, in order, and holds one at a time. **The record
/// count in the identity is counted here rather than taken from the file's own index**, so a
/// file whose index disagrees with its blocks cannot produce a census claiming a count the
/// records do not support.
///
/// # Errors
///
/// [`CensusFromPspError::NotOpened`] when the file is not a readable psp, and
/// [`CensusFromPspError::RecordsStopped`] when it fails part-way through, naming how many
/// records had been read.
pub fn census_from_psp(
    path: &Path,
    plan: &CensusPlan,
    segmentation: &Segmentation,
) -> Result<CensusOfStoredPileup, CensusFromPspError> {
    let not_opened = |source: PspReadError| CensusFromPspError::NotOpened {
        path: path.to_path_buf(),
        source: Box::new(source),
    };
    // **The digest of the header as it stands in the file, not of the one this process decodes.**
    // The psp writer amends the header before encoding it, so the two are not the same bytes,
    // and a census naming a header no psp carries would make every freshness check answer
    // *rebuild* for ever.
    let header = psp::header_digest(path).map_err(not_opened)?;
    let mut reader = PspReader::open(path).map_err(not_opened)?;

    let sample = reader.header().sample.clone();
    // **Walk-local numbers, which is the grain a census is keyed at.** A census belongs to one
    // sample, so the numbers its own psp gave its read groups are the right ones; renumbering
    // across a cohort is the calling stage's business and happens later, elsewhere.
    let read_groups: Vec<ReadGroupId> = reader
        .header()
        .read_groups
        .iter()
        .map(|group| group.walk_local_id)
        .collect();

    let mut writer = plan.writer_for(sample.clone(), read_groups.clone(), segmentation);
    let mut records = 0u64;
    for streamed in reader.records().map_err(not_opened)? {
        let streamed = streamed.map_err(|source| CensusFromPspError::RecordsStopped {
            path: path.to_path_buf(),
            records,
            source: Box::new(source),
        })?;
        records += 1;
        // **Every record's body is decoded here**, because a census is fed from bodies. The
        // `None` arm belongs to a selective walk that skipped a body on its head alone, which
        // this walk never asks for.
        if let Some(record) = streamed.record.as_ref() {
            writer.add_locus(record);
        }
    }

    Ok(CensusOfStoredPileup {
        evidence: writer.finish(),
        identity: PileupIdentity { header, records },
        sample,
        read_groups,
    })
}

#[cfg(test)]
mod tests {
    //! **Plan step A1**: the producer reads a stored psp and builds a census from it.
    //!
    //! What the two producers *agree* on is step A2's question and is tested there; these
    //! establish that this one runs at all, names the file it read, and refuses what it cannot
    //! read.

    use super::*;
    use crate::ng::psp::{self, PspReader};
    use crate::ng::run::test_fixtures::{a_census_plan_over, gatherer_over};
    use crate::pop_var_caller_exp::test_fixtures::a_cohort_on_disk;

    /// **The census names the psp it was read from, by that file's own bytes.**
    ///
    /// The identity is the digest of the header as it stands on disk and the number of records
    /// the file holds, and both are compared against what the walk that wrote the psp reported.
    /// **This is the check that would fail if the digest were taken from a decoded `Header`
    /// rather than from the file**: the writer records the compression level into the header
    /// before encoding it, so a value rebuilt in memory is not the value in the file, and a
    /// census naming it would send every freshness check to *rebuild* for ever.
    #[test]
    fn the_census_it_builds_names_the_psp_it_read() {
        let cohort = a_cohort_on_disk();
        let (segmentation, plan) = a_census_plan_over(&cohort.reference, &cohort.catalog);
        let psp = cohort.directory.path().join("zeta.psp");

        let (stats, _) = gatherer_over(
            &cohort.alignments[0],
            &cohort.reference,
            &segmentation,
            Some(&plan),
        )
        .write_psp(&psp, None)
        .expect("the walk writes its psp");

        let built = census_from_psp(&psp, &plan, &segmentation).expect("the psp is readable");

        assert_eq!(
            built.identity.header, stats.header_digest,
            "the digest is of the header the writer actually wrote",
        );
        assert_eq!(
            built.identity.records, stats.records,
            "the record count is the one the walk stored",
        );
    }

    /// **The digest taken from the file's bytes equals the one taken from its decoded header.**
    ///
    /// These are two routes to one number and they must not come apart: the walk-time producer
    /// takes it from the writer, this producer takes it from the file, and a fit compares the
    /// two. A format change that made a header re-encode differently from how it was written
    /// would break the pairing silently, and this is what says so.
    #[test]
    fn the_digest_off_the_file_is_the_digest_of_its_own_header() {
        let cohort = a_cohort_on_disk();
        let (segmentation, plan) = a_census_plan_over(&cohort.reference, &cohort.catalog);
        let psp = cohort.directory.path().join("zeta.psp");
        let _ = gatherer_over(
            &cohort.alignments[0],
            &cohort.reference,
            &segmentation,
            Some(&plan),
        )
        .write_psp(&psp, None)
        .expect("the walk writes its psp");

        let from_the_file = psp::header_digest(&psp).expect("the header reads");
        let re_encoded = PileupIdentity::of_header(
            &PspReader::open(&psp)
                .expect("the psp opens")
                .header()
                .encode()
                .expect("the header it was written with re-encodes"),
            0,
        );

        assert_eq!(from_the_file, re_encoded.header);
    }

    /// **The sample and its read groups come from the psp, not from anything handed in.**
    ///
    /// A census belongs to one individual and is keyed by read group, so a producer that took
    /// either from its caller could build a census under a name the evidence does not belong to.
    #[test]
    fn the_sample_and_its_read_groups_come_from_the_psp() {
        let cohort = a_cohort_on_disk();
        let (segmentation, plan) = a_census_plan_over(&cohort.reference, &cohort.catalog);
        let psp = cohort.directory.path().join("zeta.psp");
        let _ = gatherer_over(
            &cohort.alignments[0],
            &cohort.reference,
            &segmentation,
            Some(&plan),
        )
        .write_psp(&psp, None)
        .expect("the walk writes its psp");

        let stored = PspReader::open(&psp).expect("the psp opens");
        let expected_sample = stored.header().sample.clone();
        let expected_groups: Vec<ReadGroupId> = stored
            .header()
            .read_groups
            .iter()
            .map(|group| group.walk_local_id)
            .collect();
        drop(stored);

        let built = census_from_psp(&psp, &plan, &segmentation).expect("the psp is readable");

        assert_eq!(built.sample, expected_sample);
        assert_eq!(built.read_groups, expected_groups);
        assert!(
            !built.read_groups.is_empty(),
            "the fixture's sample declares at least one read group, so an empty list here \
             would mean the header's table was not read at all",
        );
    }

    /// **A file that is not a psp is refused, and the refusal names the file.**
    #[test]
    fn a_file_that_is_not_a_psp_is_refused() {
        let cohort = a_cohort_on_disk();
        let (segmentation, plan) = a_census_plan_over(&cohort.reference, &cohort.catalog);
        let not_a_psp = cohort.directory.path().join("not-a-psp.psp");
        std::fs::write(&not_a_psp, b"this is not a psp at all").expect("the scratch dir is ours");

        let error = census_from_psp(&not_a_psp, &plan, &segmentation)
            .expect_err("a file of prose is not a psp");

        assert!(
            matches!(&error, CensusFromPspError::NotOpened { path, .. } if path == &not_a_psp),
            "{error:?}",
        );
    }
}

#[cfg(test)]
mod the_two_producers_agree {
    //! **Plan step A2, and the question it answers is about the psp** (spec §7.12).
    //!
    //! One sample's census built while its reads are walked, and built again from the psp that
    //! walk wrote, must be the same file byte for byte. **That is what says the psp holds
    //! everything a census needs**: a field the record format drops survives every other test —
    //! the psp still reads back, the caller still calls, the census still writes — and shows up
    //! only here, as two files that differ.
    //!
    //! **The fixture is the varying cohort and not the plain one**, because the field most
    //! likely not to survive the round trip is a read's length at a repeat tract, and only that
    //! fixture has a tract for a read to have a length at.

    use super::*;
    use crate::ng::parameter_estimation::joint::census_file::write_census;
    use crate::ng::run::test_fixtures::{a_census_plan_over, gatherer_over};
    use crate::pop_var_caller_exp::test_fixtures::a_varying_cohort_on_disk;

    /// Build both censuses for one sample and return the two files' bytes.
    ///
    /// The walk writes its own; this then re-reads the psp it left and writes the second
    /// through the same `write_census`, so the comparison is of two census *files* rather than
    /// of two in-memory values that a writing difference could still separate.
    fn both_censuses_for(which: usize) -> (Vec<u8>, Vec<u8>) {
        let cohort = a_varying_cohort_on_disk();
        let (segmentation, plan) = a_census_plan_over(&cohort.reference, &cohort.catalog);
        let psp = cohort.directory.path().join(format!("sample{which}.psp"));
        let walked = cohort
            .directory
            .path()
            .join(format!("sample{which}.census"));

        let _ = gatherer_over(
            &cohort.alignments[which],
            &cohort.reference,
            &segmentation,
            Some(&plan),
        )
        .write_psp(&psp, Some(&walked))
        .expect("the walk writes both files");

        let rebuilt = census_from_psp(&psp, &plan, &segmentation).expect("the psp is readable");
        let mut from_the_psp = Vec::new();
        write_census(&rebuilt.evidence, Some(rebuilt.identity), &mut from_the_psp)
            .expect("a vector accepts every write");

        let from_the_walk = std::fs::read(&walked).expect("the walk's census is on disk");
        (from_the_walk, from_the_psp)
    }

    /// **The sample that carries reads: the two censuses are one file.**
    #[test]
    fn a_sample_with_reads_gets_the_same_census_either_way() {
        let (from_the_walk, from_the_psp) = both_censuses_for(0);

        assert_eq!(
            from_the_walk.len(),
            from_the_psp.len(),
            "the two censuses are different lengths, so the psp is not round-tripping some \
             field the census records",
        );
        assert!(
            from_the_walk == from_the_psp,
            "the two censuses are the same length and differ at byte {}",
            from_the_walk
                .iter()
                .zip(&from_the_psp)
                .position(|(walked, rebuilt)| walked != rebuilt)
                .expect("the lengths match and the vectors differ, so a first difference exists"),
        );
    }

    /// **The second sample too**, because the two carry their variation in different places and
    /// a producer that read only the first record of a block would pass on one of them alone.
    #[test]
    fn the_other_sample_gets_the_same_census_either_way() {
        let (from_the_walk, from_the_psp) = both_censuses_for(1);
        assert_eq!(from_the_walk, from_the_psp);
    }

    /// **The census the psp produces is not empty**, which is what stops the test above passing
    /// on two censuses that both recorded nothing.
    ///
    /// A census whose every position reads *never walked* would compare equal to another of the
    /// same shape while proving nothing about the psp at all. This asserts the fixture's walk
    /// actually reached the kept loci — the depth codes are not all the never-walked one.
    #[test]
    fn the_census_being_compared_holds_something() {
        let cohort = a_varying_cohort_on_disk();
        let (segmentation, plan) = a_census_plan_over(&cohort.reference, &cohort.catalog);
        let psp = cohort.directory.path().join("sample0.psp");
        let _ = gatherer_over(
            &cohort.alignments[0],
            &cohort.reference,
            &segmentation,
            Some(&plan),
        )
        .write_psp(&psp, None)
        .expect("the walk writes its psp");

        let built = census_from_psp(&psp, &plan, &segmentation).expect("the psp is readable");

        assert!(
            built.identity.records > 0,
            "the fixture's walk stored no records at all, so nothing was fed to either producer",
        );
        assert!(
            !built.evidence.strata().is_empty(),
            "the census holds no repeat-tract stratum, so the comparison above says nothing \
             about the field it was built to catch — a read's length at a tract",
        );
        let strata = built.evidence.strata();
        let groups = built.evidence.read_groups();
        let mut evidence = built.evidence;
        let mut reads_at_tracts: u32 = 0;
        for group in groups {
            reads_at_tracts += evidence
                .with_strata(group, &strata, |lent| {
                    lent.iter()
                        .flat_map(|section| {
                            (0..section.len()).map(|locus| section.offsets(locus).total())
                        })
                        .sum::<u32>()
                })
                .expect("the sections are resident");
        }
        assert!(
            reads_at_tracts > 0,
            "no read reached a kept tract, so a producer that dropped every tract read would \
             pass the byte-for-byte comparison",
        );
    }
}
