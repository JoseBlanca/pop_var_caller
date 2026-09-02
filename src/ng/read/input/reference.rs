//! [`OpenReference`] — the reference every file in a run is opened against,
//! **including the one shared copy of its bases**.
//!
//! The read-input layer needs two different things from a reference. It needs
//! the *description* — the contig table each file's `@SQ` is validated against,
//! which is [`ReferenceInfo`], a plain data record. And, for CRAM only, it
//! needs the *bases*, because a CRAM stores its reads as differences from the
//! reference and cannot be decoded without them.
//!
//! **"Open" is the difference between the two, and it is the same "open" as
//! [`AlignmentFile::open`](super::open_bam::AlignmentFile::open)**: a
//! `ReferenceInfo` merely *describes* a reference — names, lengths, digests, a
//! path — while an `OpenReference` can be *read from*. It holds the FASTA open
//! and hands out sequence. Holding the two together is not tidiness: they must
//! be the same reference, and separate arguments could be mismatched, which
//! would decode every read against the wrong bases with nothing to catch it.
//!
//! The bases are the expensive half, and they are expensive in a way that is
//! easy to miss: noodles' [`fasta::Repository`] is a whole-contig memoising
//! cache with **no eviction**, so every contig a decode touches stays resident
//! for the repository's life. One repository per file therefore costs
//! `files × genome`, which for a 51-sample tomato cohort is ~38 GiB and an
//! OOM kill. One repository per *run* costs `genome`, once — the reference is
//! the same for every file by construction, since the `@SQ` gate has just
//! proved it — and bounding *that* to the contig in hand is the second half,
//! below.
//!
//! That is what this type exists to make true by construction: a caller holds
//! one `OpenReference` for the whole run and hands it to every
//! [`AlignmentFile::open`](super::open_bam::AlignmentFile::open), so there is
//! no per-file arm in which a second repository could be built. The
//! repository is built **lazily**, on the first CRAM open, so a BAM-only run
//! against a reference with no `.fai` still works — it never asks for bases.
//!
//! Cheap to clone: the description is an `Arc` and the bases are shared behind
//! one, so a clone is a pointer bump onto *the same* cache rather than a second
//! one.
//!
//! ## One contig resident, not the whole genome
//!
//! Holding the *genome* was never required either — that too was ours, not
//! noodles'. What noodles needs at any moment is the contig of the slice it is
//! decoding: `get_slice_reference_sequence` asks the repository for one contig
//! by name and holds it for that slice. Everything beyond that is cache policy,
//! and the policy of an unevicted `Repository` is "keep every contig ever
//! touched", which over a genome-wide walk is the whole genome.
//!
//! So [`bases_for_contig`](OpenReference::bases_for_contig) clears the cache
//! when the run moves to a new contig, keeping one chromosome resident instead
//! of all of them — 86.6 MiB rather than 746 MiB for tomato, and ~238 MiB
//! rather than ~3 GiB for a human reference. Production's pileup loop already
//! does exactly this (`bam::alignment_input::build_fasta_repository`), and it is
//! what [`WindowedRefSeq`](crate::ng::ref_seq::WindowedRefSeq)'s `evict_before`
//! does for the walk's own view of the bases.
//!
//! **It costs no extra reading in an ordered walk.** A contig is cleared only
//! when a query for a *different* contig arrives, so a walk that goes chr1 →
//! chr2 → … reads each contig exactly once, as it did with no eviction at all.
//! A cohort walk is region-outer and sample-inner, so all k files are on the
//! same contig at the same time and the transition happens once per contig for
//! the whole cohort, not once per file.
//!
//! **The assumption that buys it is contig-ordered access**, which every caller
//! in this tree makes. A caller that alternates between contigs — a
//! *contig-parallel* walk is the realistic one — would clear a chromosome its
//! neighbour is still working through and re-read it, over and over. Such a
//! caller must build its reference with
//! [`unbounded`](OpenReference::unbounded), and will want a k-slot policy (one
//! repository per worker contig) rather than this one; the shape is a map of
//! contig → `Repository` in place of the single `OnceLock` here. In-flight
//! decodes are never *broken* by a clear, only slowed: noodles holds each
//! contig by `Arc` for as long as it is decoding against it.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Arc, OnceLock};

use noodles_fasta as fasta;

use crate::bam::alignment_input::build_fasta_repository;
use crate::bam::errors::AlignmentInputError;
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::types::ContigId;

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
pub struct OpenReference {
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
    /// Which contig the cached bases are for, or [`NO_CONTIG`] before the first
    /// query — the state behind the one-contig bound (module docs).
    ///
    /// An atomic rather than a lock because it is read on every region query
    /// and the update is a single compare: two threads that both observe a
    /// transition both clear, which is idempotent and needs no mutual
    /// exclusion.
    resident_contig: AtomicU32,
    /// Whether to clear on contig transition at all. `false` restores
    /// keep-everything for a caller whose access is not contig-ordered
    /// ([`unbounded`](OpenReference::unbounded)).
    bound_to_one_contig: bool,
}

/// `resident_contig` before any query has narrowed the cache. Not a legal
/// [`ContigId`] in any reference this project can read — a contig table that
/// large would need 4 billion entries.
const NO_CONTIG: u32 = u32::MAX;

impl OpenReference {
    /// Wrap a run's reference description. Nothing is read here; the bases are
    /// opened on first use, and only one contig of them is kept resident
    /// (module docs).
    pub fn new(info: Arc<ReferenceInfo>) -> Self {
        Self::with_bound(info, true)
    }

    /// [`new`](Self::new), but keeping **every** contig the run touches.
    ///
    /// For a caller whose queries are not contig-ordered — a contig-parallel
    /// walk — where the one-contig bound would clear a chromosome another
    /// worker is still reading and re-read it repeatedly. Costs the whole
    /// genome resident, which is the price of that access pattern until a
    /// k-slot policy exists (module docs).
    pub fn unbounded(info: Arc<ReferenceInfo>) -> Self {
        Self::with_bound(info, false)
    }

    fn with_bound(info: Arc<ReferenceInfo>, bound_to_one_contig: bool) -> Self {
        Self {
            inner: Arc::new(Inner {
                info,
                bases: OnceLock::new(),
                resident_contig: AtomicU32::new(NO_CONTIG),
                bound_to_one_contig,
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

    /// The run's repository, holding **this contig's** bases and no other's.
    ///
    /// Called per region query. When the contig differs from the last query's,
    /// the cache is cleared first, so the resident bases are one chromosome
    /// rather than every chromosome the walk has passed (module docs). A query
    /// for the contig already resident does nothing but hand back the handle.
    ///
    /// `None` if the bases were never opened — a BAM-only file, which never
    /// asks. Infallible otherwise: `open` already proved they open, and the
    /// result is memoised.
    pub(crate) fn bases_for_contig(&self, contig: ContigId) -> Option<fasta::Repository> {
        let repository = self.inner.bases.get()?;
        // **The load is not an optimisation of the swap; it is the same answer.**
        // Where the stored contig already is this one, the swap writes back what it
        // read and the comparison below is false — so reading first can only skip a
        // no-op. What it buys is that the common case takes the cache line shared
        // instead of exclusive: this is called per region query by every sample's
        // draw, and a cohort walk keeps every sample on one contig at a time, so all
        // but the first query of each contig is that case.
        if self.inner.bound_to_one_contig
            && self.inner.resident_contig.load(Ordering::Relaxed) != contig.0
        {
            // `swap`, not load-then-store: two threads observing the same
            // transition both clear, which is harmless, but neither may miss it.
            let previous = self.inner.resident_contig.swap(contig.0, Ordering::Relaxed);
            if previous != contig.0 {
                // Drops this repository's *cache entries* only. A decode still
                // running against the outgoing contig holds it by `Arc` and is
                // unaffected — it simply keeps the bases alive a little longer.
                repository.clear();
            }
        }
        Some(repository.clone())
    }

    /// Whether the bases have been opened yet — the observable form of "one
    /// repository per run", which `every_caller_gets_one_repository` asserts
    /// directly.
    #[cfg(test)]
    pub(crate) fn bases_opened(&self) -> bool {
        self.inner.bases.get().is_some()
    }
}

impl std::fmt::Debug for OpenReference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OpenReference")
            .field("contigs", &self.inner.info.contigs.len())
            .field("fasta_path", &self.inner.info.fasta_path)
            .field("bases_opened", &self.inner.bases.get().is_some())
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

    /// **The one-contig bound.** Moving to a new contig drops the previous
    /// one's bases; staying on a contig keeps them. Both halves matter: the
    /// first is the memory saving, the second is what makes it free — without
    /// it every query would re-read the chromosome.
    #[test]
    fn moving_to_a_new_contig_drops_the_previous_contig_s_bases() {
        let (_dir, reference) = fixture_reference(true);
        reference.bases().expect("the fixture has a FASTA");

        let first = reference
            .bases_for_contig(ContigId(0))
            .expect("the bases are open");
        first.get(b"chr1").expect("chr1 is in the fixture").unwrap();
        first.get(b"chr2").expect("chr2 is in the fixture").unwrap();
        assert_eq!(first.len(), 2, "both fetched, both cached");

        // Still chr1: nothing is dropped, so a second query costs no re-read.
        let same = reference
            .bases_for_contig(ContigId(0))
            .expect("the bases are open");
        assert_eq!(same.len(), 2);

        // On to chr2 — and chr1's 90 Mb equivalent goes with it.
        let next = reference
            .bases_for_contig(ContigId(1))
            .expect("the bases are open");
        assert_eq!(next.len(), 0, "the cache was cleared at the transition");
    }

    /// `unbounded` is the escape hatch for a caller whose queries are not
    /// contig-ordered: it never clears, whatever contig it is asked for.
    #[test]
    fn an_unbounded_reference_keeps_every_contig() {
        let (_dir, bounded) = fixture_reference(true);
        let reference = OpenReference::unbounded(Arc::new(bounded.info().clone()));

        let bases = reference.bases().expect("the fixture has a FASTA");
        bases.get(b"chr1").expect("chr1 is in the fixture").unwrap();

        let after_transition = reference
            .bases_for_contig(ContigId(1))
            .expect("the bases are open");
        assert_eq!(after_transition.len(), 1, "nothing is ever dropped");
    }
}
