//! [`RunReference`] — the reference every file in a run is opened against,
//! **including the one shared copy of its bases**.
//!
//! The read-input layer needs two different things from a reference. It needs
//! the *description* — the contig table each file's `@SQ` is validated against,
//! which is [`ReferenceInfo`], a plain data record. And, for CRAM only, it
//! needs the *bases*, because a CRAM stores its reads as differences from the
//! reference and cannot be decoded without them.
//!
//! The bases are the expensive half, and they are expensive in a way that is
//! easy to miss: noodles' [`fasta::Repository`] is a whole-contig memoising
//! cache with **no eviction**, so every contig a decode touches stays resident
//! for the repository's life. One repository per file therefore costs
//! `files × genome`, which for a 51-sample tomato cohort is ~38 GiB and an
//! OOM kill. One repository per *run* costs `genome`, once — the reference is
//! the same for every file by construction, since the `@SQ` gate has just
//! proved it.
//!
//! That is what this type exists to make true by construction: a caller holds
//! one `RunReference` for the whole run and hands it to every
//! [`AlignmentFile::open`](super::open_bam::AlignmentFile::open), so there is
//! no per-file arm in which a second repository could be built. The
//! repository is built **lazily**, on the first CRAM open, so a BAM-only run
//! against a reference with no `.fai` still works — it never asks for bases.
//!
//! Cheap to clone: the description is an `Arc` and the bases are shared behind
//! one, so a clone is a pointer bump onto *the same* cache rather than a second
//! one.
//!
//! ## What this deliberately does not do yet
//!
//! One repository still holds every contig it has touched — 746 MiB for tomato,
//! ~3 GiB for a human reference, more for a polyploid crop. `Repository::clear`
//! exists, and production's own pileup loop already calls it on contig
//! transition to keep exactly one contig resident
//! (`bam::alignment_input::build_fasta_repository`); a cohort walk is
//! region-outer and sample-inner, so every file is on the same contig at the
//! same time and that policy would cap the resident bases at the largest single
//! contig (~87 MiB for tomato) instead of the whole genome. That is a second,
//! independent factor of ~8.6 on this reference — and it is what
//! [`WindowedRefSeq`](crate::ng::ref_seq::WindowedRefSeq)'s `evict_before`
//! already does for the *walk's* view of the bases.
//!
//! It is not done here because eviction needs a signal this layer does not
//! have: only the caller knows when the walk has left a contig for good, and
//! clearing early just re-reads. Adding it means giving `RunReference` an
//! evict-on-transition entry point and a caller that calls it — a change to the
//! walk's contract, not to this type alone. Worth doing when a reference is
//! big enough for `genome` to be the problem; for tomato, `files × genome` was.

use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

use noodles_fasta as fasta;

use crate::bam::alignment_input::build_fasta_repository;
use crate::bam::errors::AlignmentInputError;
use crate::ng::reference_info::ReferenceInfo;

/// Why a run's reference could not supply the bases a CRAM needs to decode.
///
/// Two distinct faults, kept apart because they call for different fixes: the
/// caller described the reference by its `.fai` alone and there are no bases to
/// be had, or there are bases but the indexed FASTA could not be opened.
#[derive(Debug, thiserror::Error)]
pub enum ReferenceBasesError {
    /// The reference was read from a `.fai` alone. A `.fai` describes a
    /// genome's geometry but holds no sequence, so nothing here can serve a
    /// CRAM decode.
    #[error("the reference was described by a `.fai` alone, which holds no bases")]
    NoFasta,
    /// The FASTA is named but could not be opened as an indexed reader —
    /// usually a missing sibling `.fai`.
    #[error("reference FASTA '{fasta}' cannot be read: {source}")]
    Build {
        /// The FASTA that could not be opened.
        fasta: PathBuf,
        /// What went wrong opening it.
        #[source]
        source: AlignmentInputError,
    },
}

/// The reference a run reads against: its description, plus the single shared
/// cache of its bases that every CRAM in the run decodes against.
///
/// Built once per run and handed to every file's open (module docs). Cloning
/// shares — it does not copy the cache — so passing it by value is as cheap as
/// passing it by reference.
#[derive(Clone)]
pub struct RunReference {
    inner: Arc<Inner>,
}

struct Inner {
    info: Arc<ReferenceInfo>,
    /// Built on the first CRAM open and reused by every later one.
    ///
    /// A `OnceLock` rather than an eager field so a BAM-only run never touches
    /// the FASTA, and **only successes are stored**, so a transient failure to
    /// open the FASTA can be retried rather than being remembered as the
    /// answer — the same policy
    /// [`ReferenceInfoCache`](crate::ng::reference_info::ReferenceInfoCache)
    /// states for its own reads.
    bases: OnceLock<fasta::Repository>,
}

impl RunReference {
    /// Wrap a run's reference description. Nothing is read here; the bases are
    /// opened on first use.
    pub fn new(info: Arc<ReferenceInfo>) -> Self {
        Self {
            inner: Arc::new(Inner {
                info,
                bases: OnceLock::new(),
            }),
        }
    }

    /// The reference description — the contig table, the digests, the FASTA
    /// path.
    pub fn info(&self) -> &ReferenceInfo {
        &self.inner.info
    }

    /// The FASTA this reference was read from, when it was read from one.
    /// Private: a caller that wants it asks [`info()`](Self::info), which
    /// carries it; this exists so [`bases()`](Self::bases) reads the field in
    /// one place.
    fn fasta_path(&self) -> Option<&Path> {
        self.inner.info.fasta_path.as_deref()
    }

    /// The run's **one** reference-bases repository, built on first call.
    ///
    /// The returned `Repository` is a handle onto the shared cache (it is an
    /// `Arc` inside), so every caller gets the same memoised contigs. A losing
    /// racer's freshly built repository is dropped before it has cached
    /// anything, so the race costs an indexed-reader open and nothing else.
    ///
    /// **Sharing the cache shares its lock.** noodles takes the repository's
    /// *write* lock across the adapter read, so the first fetch of a contig
    /// blocks every other reader of that repository until the contig is in.
    /// That is a convoy of one whole-contig read per contig per run — where
    /// per-file repositories had each file pay that read in full, so the total
    /// work drops even as the blocking appears. Every fetch after the first
    /// takes the read lock and is concurrent. Today's callers are
    /// single-threaded, so it costs nothing yet; it is worth re-timing when one
    /// is not.
    ///
    /// # Errors
    ///
    /// [`ReferenceBasesError::NoFasta`] for a `.fai`-only reference;
    /// [`ReferenceBasesError::Build`] if the FASTA cannot be opened.
    pub fn bases(&self) -> Result<fasta::Repository, ReferenceBasesError> {
        if let Some(repository) = self.inner.bases.get() {
            return Ok(repository.clone());
        }
        let fasta = self
            .fasta_path()
            .ok_or(ReferenceBasesError::NoFasta)?
            .to_path_buf();
        let built = build_fasta_repository(&fasta)
            .map_err(|source| ReferenceBasesError::Build { fasta, source })?;
        Ok(self.inner.bases.get_or_init(|| built).clone())
    }

    /// Whether the bases have been opened yet — the observable form of "one
    /// repository per run", which `every_caller_gets_one_repository` asserts
    /// directly.
    #[cfg(test)]
    pub(crate) fn bases_opened(&self) -> bool {
        self.inner.bases.get().is_some()
    }
}

impl std::fmt::Debug for RunReference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RunReference")
            .field("contigs", &self.inner.info.contigs.len())
            .field("fasta_path", &self.inner.info.fasta_path)
            .field("bases_opened", &self.inner.bases.get().is_some())
            .finish()
    }
}

impl From<ReferenceInfo> for RunReference {
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
        let reference = RunReference::from(info);
        assert!(matches!(
            reference.bases(),
            Err(ReferenceBasesError::NoFasta)
        ));
        assert!(!reference.bases_opened());
    }

    /// **The property the type exists for.** Two `bases()` calls — which is
    /// what two CRAM opens do — yield handles onto *one* cache, so a contig
    /// fetched through the first is already resident in the second.
    #[test]
    fn every_caller_gets_one_repository() {
        let (_dir, reference) = fixture_reference(true);
        assert!(!reference.bases_opened(), "nothing is read until asked");

        let first = reference.bases().expect("the fixture has a FASTA");
        assert!(reference.bases_opened());
        let name = b"chr1";
        first.get(name).expect("chr1 is in the fixture").unwrap();
        assert_eq!(first.len(), 1);

        let second = reference.bases().expect("the fixture has a FASTA");
        assert_eq!(
            second.len(),
            1,
            "a second caller sees the first caller's cached contig, \
             which is only true if they share one repository"
        );

        // And a *clone* of the handle shares too — the trap a plain `Clone`
        // derive over a bare `OnceLock` would have set.
        let cloned = reference.clone();
        assert!(cloned.bases_opened());
        assert_eq!(cloned.bases().expect("shared").len(), 1);
    }
}
