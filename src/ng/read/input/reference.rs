//! [`OpenReference`] — the reference every file in a run is opened against,
//! and the one question that is asked of it at open.
//!
//! The read-input layer needs two different things from a reference. It needs
//! the *description* — the contig table each file's `@SQ` is validated against,
//! which is [`ReferenceInfo`], a plain data record. And, for CRAM only, it
//! needs to know that the *bases* can be had at all, because a CRAM stores its
//! reads as differences from the reference and cannot be decoded without them.
//!
//! **"Open" is the difference between the two, and it is the same "open" as
//! [`AlignmentFile::open`](super::open_bam::AlignmentFile::open)**: a
//! `ReferenceInfo` merely *describes* a reference — names, lengths, digests, a
//! path — while an `OpenReference` has been proved readable. Holding the two
//! together is not tidiness: they must be the same reference, and separate
//! arguments could be mismatched, which would decode every read against the
//! wrong bases with nothing to catch it.
//!
//! ## What this type is not, since it used to be
//!
//! Until 2026-09-07 it also **held the bases**: one shared `fasta::Repository`
//! per run, bounded to the contig in hand, and every CRAM cursor took a handle
//! onto it. That was the answer to a real fault — a repository memoises whole
//! contigs and never evicts, so one per file costs `files × genome`, which for
//! a 51-sample tomato cohort was ~38 GiB and an OOM kill — but it was an
//! expensive answer: noodles holds a contig at about two bytes a base, so one
//! chromosome resident is 198 MB for tomato and 493 MB for human.
//!
//! The decode never needed a chromosome. It needs the span of the slice it is
//! decoding, about 9 kb on a real coordinate-sorted CRAM, and it now reads that
//! span from **a reference reader of its own**, minted by the same factory
//! every other consumer of bases in this project is served by
//! ([`alignment_cursor.md`](../../../../../doc/devel/ng/spec/alignment_cursor.md)
//! §10). So the repository, the one-contig bound and the escape hatch from it
//! are all gone, and what is left here is a description plus a proof.
//!
//! ## The proof, and what it does not cover
//!
//! [`check_bases_can_be_read`](OpenReference::check_bases_can_be_read) asks two
//! things and reads no sequence: does this reference name a FASTA at all, and
//! does that FASTA's `.fai` parse. A reference described by its `.fai` alone
//! fails the first — there is no sequence anywhere to be had — and a FASTA
//! whose sibling index is missing or malformed fails the second. Both are
//! faults at the moment a CRAM is opened rather than mysteries at the first
//! query, which is the whole point of asking here.
//!
//! **Opening the FASTA itself is proved one level later**, by the zero-length
//! fetch [`AlignmentFile::cursor`](super::open_bam::AlignmentFile::cursor)
//! makes through every reference reader it mints. That is a stronger check than
//! this one — it proves the bases of the cursor's own contig can be served, not
//! merely that the geometry parses — and it is where a reader that cannot open
//! its FASTA is caught.
//!
//! Cheap to clone: the description is behind an `Arc`, so passing one by value
//! is a pointer bump.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::ng::ref_seq::WindowedRefSeq;
use crate::ng::reference_info::ReferenceInfo;

/// Why a run's reference could not supply the bases a CRAM needs to decode.
///
/// Two distinct faults, kept apart because they call for different fixes: the
/// caller described the reference by its `.fai` alone and there are no bases to
/// be had, or there are bases but the index that says where they lie could not
/// be read.
#[derive(Debug, thiserror::Error)]
pub enum ReferenceBasesError {
    /// The reference was read from a `.fai` alone. A `.fai` describes a
    /// genome's geometry but holds no sequence, so nothing here can serve a
    /// CRAM decode.
    #[error("the reference was described by a `.fai` alone, which holds no bases")]
    NoFasta,
    /// The FASTA is named but its `.fai` could not be read — usually because
    /// the sibling index is missing.
    #[error("reference FASTA '{fasta}' cannot be read: its `.fai` did not load: {source}")]
    Build {
        /// The FASTA whose index could not be read.
        fasta: PathBuf,
        /// What went wrong reading it.
        #[source]
        source: std::io::Error,
    },
}

/// The reference a run reads against: its description, and the proof — taken
/// once per CRAM opened — that bases for it can be reached.
///
/// Built once per run and handed to every file's open (module docs). Cloning
/// shares the description rather than copying it, so passing it by value is as
/// cheap as passing it by reference.
#[derive(Clone)]
pub struct OpenReference {
    info: Arc<ReferenceInfo>,
}

impl OpenReference {
    /// Wrap a run's reference description. Nothing is read here; the check
    /// below is made when a CRAM is opened against it.
    pub fn new(info: Arc<ReferenceInfo>) -> Self {
        Self { info }
    }

    /// The reference description — the contig table, the digests, the FASTA
    /// path.
    pub fn info(&self) -> &ReferenceInfo {
        &self.info
    }

    /// The FASTA this reference was read from, when it was read from one.
    /// Private: a caller that wants it asks [`info()`](Self::info), which
    /// carries it; this exists so the check below reads the field in one place.
    fn fasta_path(&self) -> Option<&Path> {
        self.info.fasta_path.as_deref()
    }

    /// Prove that a CRAM opened against this reference has bases to decode
    /// with: the reference names a FASTA, and that FASTA's `.fai` reads.
    ///
    /// Called once per CRAM open. No sequence is read — the `.fai` of a
    /// GRCh38-shaped reference is 2,580 lines and parses in about 137 µs, which
    /// is why this is not memoised: a thousand-sample cohort pays it once per
    /// file at startup and never again.
    ///
    /// # Errors
    ///
    /// [`ReferenceBasesError::NoFasta`] for a `.fai`-only reference;
    /// [`ReferenceBasesError::Build`] when the index does not read.
    pub(crate) fn check_bases_can_be_read(&self) -> Result<(), ReferenceBasesError> {
        let fasta = self
            .fasta_path()
            .ok_or(ReferenceBasesError::NoFasta)?
            .to_path_buf();
        WindowedRefSeq::read_index(&fasta)
            .map(|_| ())
            .map_err(|source| ReferenceBasesError::Build { fasta, source })
    }
}

impl std::fmt::Debug for OpenReference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OpenReference")
            .field("contigs", &self.info.contigs.len())
            .field("fasta_path", &self.info.fasta_path)
            .finish()
    }
}

impl From<ReferenceInfo> for OpenReference {
    fn from(info: ReferenceInfo) -> Self {
        Self::new(Arc::new(info))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::input::test_fixtures::fixture_reference;

    /// A `.fai`-only reference names the fault for what it is rather than
    /// failing as a missing file.
    #[test]
    fn a_fai_only_reference_has_no_bases() {
        let info = ReferenceInfo {
            md5: None,
            contigs: Vec::new(),
            fasta_path: None,
        };
        let reference = OpenReference::from(info);
        assert!(matches!(
            reference.check_bases_can_be_read(),
            Err(ReferenceBasesError::NoFasta)
        ));
    }

    /// A reference that names a FASTA with a readable index passes.
    #[test]
    fn a_reference_with_a_fasta_and_its_index_can_be_read() {
        let (_dir, reference) = fixture_reference(true);
        reference
            .check_bases_can_be_read()
            .expect("the fixture has a FASTA and a .fai beside it");
    }

    /// **The half the `.fai`-only test cannot reach.** A FASTA is named and
    /// present, and its index is not — which is the ordinary way this fails on
    /// real input, and it must be a fault at open with the FASTA named rather
    /// than a missing-file error from somewhere inside a decode.
    #[test]
    fn a_fasta_whose_index_is_missing_is_a_fault_naming_the_fasta() {
        let (_dir, reference) = fixture_reference(true);
        let fasta = reference
            .info()
            .fasta_path
            .clone()
            .expect("the fixture names its FASTA");
        std::fs::remove_file(crate::ng::reference_info::sibling_fai_path(&fasta))
            .expect("the fixture wrote a .fai to remove");

        let error = reference
            .check_bases_can_be_read()
            .expect_err("the index is gone");
        match error {
            ReferenceBasesError::Build { fasta: named, .. } => assert_eq!(named, fasta),
            other => panic!("expected the index fault, got {other:?}"),
        }
    }
}
