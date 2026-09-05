//! One sample's stored observations behind the merge's source interface.
//!
//! **A source answers one question: what did this sample see next?** (`doc/devel/ng/arch/run_streaming.md`
//! §2). Direct mode's answer is [`AlignmentFilesWalker`](super::walker::AlignmentFilesWalker),
//! which mints the observations from alignment files; psp mode's is this file, which decodes
//! them from a psp somebody else already walked. **Nothing above the trait can tell the two
//! apart**, and that is the whole of spec §3.1's "the two callers differ only in what a source
//! is".
//!
//! **It is an adapter and not a reader.** [`PspReader`] already streams a psp's records in
//! coordinate order and already holds nothing that grows with the file
//! ([`psp::walk`](crate::ng::psp::walk)); what this adds is the three things the merge's trait
//! asks for that the store's own walk cannot give.
//!
//! - **A failure that names the sample and how far it had got.** The merge adds nothing to a
//!   source's error and passes it through, so in a run over a thousand samples an error saying
//!   only *block 41 would not inflate* names neither the individual to look at nor the ground
//!   already called ([`RunError::SourceFailed`], spec §9). On the path a run takes
//!   ([`PspObservationSource::over`]) the name is read from the file's own header, so the name
//!   in the error and the file it came from cannot come apart.
//! - **Refusals where the merge would otherwise assert or go quiet** — see [`PspSourceError`],
//!   which says of each whether the mistake is the file's or this crate's.
//! - **The read groups renumbered from the file's terms into the run's.** Every psp numbers its
//!   own read groups from zero, so a cohort's files collide by construction; the calling stage
//!   merges their tables into one numbering and hands each source its own map (spec §6.2). This
//!   is where a record stops being one file's and becomes the run's.
//! - **A parameter it takes and drops: the offer of a spent record back for reuse.**
//!   [`ObservationSource::next_observation`] hands a source a record the merge will not read
//!   again, and this one releases it. A decoder is the reuse hook's best customer — it fills
//!   buffers per record where a walker mints whole new ones — but the reader's own interface
//!   hands records out rather than filling ones it is given, so taking the offer is a change
//!   to [`PspReader`]'s walk and not to this adapter. The psp-mode plan defers it with the
//!   rest of psp-mode performance.
//!
//! **Generic over the walk, and that is what makes the head-only refusal reachable.** A psp's
//! records can be walked two ways: [`RecordIter`], which builds every record, and
//! [`SelectiveRecordIter`](crate::ng::psp::SelectiveRecordIter), which builds only the bodies a
//! predicate asks for and hands the rest back as a head alone. **A merge fed the second would
//! be handed a cohort with observations missing and no error anywhere** — wrong genotypes, not
//! a failure. So this type takes either and refuses the moment a body is missing, rather than
//! taking only the first and leaving the trap for whoever wires the cheap-numbers pass spec §10
//! defers.

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::psp::{PspReadError, PspReader, RecordIter, StreamedRecord};
use crate::ng::types::{GenomePosition, GenomeRegion, ReadGroupId};

use super::cohort_merge::observation_cache::ObservationSource;
use super::cohort_merge::observation_cache::{Drawn, LocusSummary};
use crate::ng::psp::RecordHead;
use crate::ng::psp::record::{LocatedRecord, RecordLayout, decode_the_body_of};
use crate::ng::psp::chain_ids::LiveSet;

/// **A stored record's head is a summary** — the claim the whole deferred-build design rests
/// on, written down as a conversion so it can be tested rather than asserted.
///
/// The merge decides everything it decides before assembling a locus from three facts about an
/// observation: the reference it covers, and the keep rule's numerator and denominator
/// ([`LocusSummary`]). A psp record carries all three in the fixed head in front of its body,
/// so a reader can produce this without decompressing the evidence it describes — which is what
/// lets a run decline to build the roughly ninety-nine positions in a hundred no sample varied
/// at (`spec/cohort_merge_psp_path.md` §2).
///
/// **The head's two counts are the same numbers, not an approximation of them**: the writer
/// derives them from the record it is about to write, by the same call the merge would make,
/// and the reader checks them against the rebuilt body whenever it does build one
/// (`psp/record.rs`, `decode_the_body_of`). `a_head_summarises_a_record_exactly` pins the
/// mapping from the other end.
impl From<&RecordHead> for LocusSummary {
    fn from(head: &RecordHead) -> Self {
        Self {
            region: head.region,
            non_reference_reads: head.non_reference_reads,
            reads_compared_with_reference: head.reads_compared_with_reference,
        }
    }
}
use super::{RunError, WalkProgress};

/// What is wrong with a psp's observations themselves, once the file has decoded.
///
/// **What these share is not whose mistake they are — it is that the merge would otherwise
/// meet them as an assertion or as silence.** Each variant's own doc says where the mistake
/// is, and they differ: an out-of-order file is damaged input, a head-only record is this
/// crate's own wiring.
///
/// The out-of-order one is why the type exists. The merge's cache refuses a backwards source
/// with a release assertion, on the ground that a source going backwards is a bug in the
/// generators that mint them ([`observation_cache`](super::cohort_merge::observation_cache),
/// `draw_next`). That ground disappears the moment observations arrive from a file written by
/// another process, possibly by another build: damaged input is refused with a message rather
/// than by aborting the process. Arch §8 records this as owed the day the psp path reaches
/// that check, and this is where it lands.
///
/// **None of these names the sample**, because [`RunError::SourceFailed`] does — this travels
/// as its cause, and a sample named twice in one rendered chain reads as two samples.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PspSourceError {
    /// A record starts before the one in front of it.
    ///
    /// **What it costs if it is not caught**: the merge draws each sample forward only until
    /// its window is covered, so an observation behind the window's left edge ends a cover
    /// early and hands a builder a locus cut short — a wrong genotype rather than a failure.
    /// The store refuses such a record at the point of writing
    /// ([`PspWriter::push`](crate::ng::psp::PspWriter::push)), so a file reaching this was
    /// damaged after it was written or written by something that is not this store.
    ///
    /// **Two start positions rather than two spans, and rendered the same way**: the order
    /// check compares where records *start*, so these are the two values it compared. A span
    /// on one side and a position on the other would leave a reader working out that
    /// `contig 0:101-101` and `contig 0 position 101` are the same fact.
    #[error(
        "its stored observations are not in coordinate order: a record starting at contig {} \
         position {} follows one starting at contig {} position {}",
        offered.contig.get(),
        offered.position.get(),
        previous.contig.get(),
        previous.position.get(),
    )]
    ObservationsOutOfOrder {
        /// Where the record in front of it starts.
        previous: GenomePosition,
        /// Where the record that goes backwards starts.
        offered: GenomePosition,
    },

    /// A record arrived with its head alone, its body never built.
    ///
    /// **A wiring mistake rather than damage**, and the only one of the three: it means the
    /// walk underneath was a selective one, whose predicate declined this record's body. A
    /// merge cannot be handed a locus with no observations in it, so this refuses rather than
    /// skipping — a skip would be an observation silently missing from one sample of the
    /// cohort, which is the class of failure `run_streaming.md` §3.3 exists to keep visible.
    ///
    /// **The message names the walk and not the file**, because the file is sound: sending
    /// somebody to rebuild a psp that is not broken is the wrong instruction.
    #[error(
        "the stored observation at {at} arrived without its body, so it holds no evidence: \
         this source was given a walk that builds only some record bodies"
    )]
    ObservationBodyNotBuilt {
        /// The record whose body was declined, as its head describes it.
        at: GenomeRegion,
    },

    /// A record names a read group its own file's table does not hold.
    ///
    /// **The number in a record is walk-local**, and the run renumbers it through the table
    /// its header carries (spec §6.2) — so a number past the end of that table is a file
    /// disagreeing with itself, and there is nothing to renumber it into. Caught here rather
    /// than by the indexing, because a panic in the middle of a cohort says nothing about
    /// which file to look at.
    #[error(
        "the stored observation at {at} names read group {names}, and this sample's psp \
         declares {in_the_table}"
    )]
    ReadGroupNotInThisFilesTable {
        /// The record whose observation names it.
        at: GenomeRegion,
        /// The walk-local number the record carries.
        names: u32,
        /// How many read groups the file's own header declares.
        in_the_table: usize,
    },

    /// The source refused a record earlier in this file, so it will not go on.
    ///
    /// **A refusal has to end the source, and this is what says so on every later draw.**
    /// Without it, a consumer that swallowed one of the two refusals above and asked again
    /// would be handed the record *after* the refused one, and the refused observation would
    /// be gone from the cohort with nothing left to say so — one sample short at one locus,
    /// a wrong genotype rather than a failure. Answering `None` instead would be the same
    /// silence in a different shape: the merge reads `None` as a sample that ran out.
    ///
    /// **Nothing reaches this today**: the merge abandons its cache at the first failure and
    /// never draws again. It exists so that the day something retries, it is told.
    #[error(
        "it refused a record earlier in this file, so what follows that record is not this \
         sample's whole evidence"
    )]
    AlreadyRefused,
}

/// One sample's observations, decoded from its psp (arch §2, §5).
///
/// Constructed once per sample and advanced by the merge alone — **one source per sample for
/// the whole run**, not one per worker and not one per building region (spec §3.4). The merge
/// is its only consumer and only moves forward, so each block of the file is decoded once and
/// the backward jump the reader is capable of never happens.
///
/// **Not `Clone`**, for the walker's reason: two sources over one sample would decode the same
/// ground twice while each told the merge a different story about how far it had got.
///
/// **Deliberately not an [`Iterator`]**, also for the walker's reason: every iterator of one
/// sample's observations is already a source through the blanket implementation
/// ([`observation_cache`](super::cohort_merge::observation_cache)), so a type that was both
/// would implement the trait twice and Rust refuses the overlap. Taking the blanket
/// implementation instead was the alternative, and it cannot work here: it drops the spare and
/// passes the error through untouched, where this has to name the sample.
///
/// **It hands over every record the file holds, and knows nothing about the run's analysed
/// regions or its reach ceiling.** Both are the calling stage's: the analysed regions are what
/// `PspVariantCaller::open` compares across the cohort, and the ceiling is read from the header
/// where the merge needs it (plan steps E1 and E4). A source is asked one question and this is
/// the whole of its answer.
/// One sample's stored records read **as summaries, with their evidence kept but not built**.
///
/// This is the reading half of the deferred build (`spec/cohort_merge_psp_path.md` §3.1). It
/// walks the file's record heads, hands back the [`LocusSummary`] each one carries, and appends
/// the record's still-compressed-out body bytes to an arena of its own — so the cohort can
/// decide which loci are worth calling before anything is decoded, and the roughly ninety-nine
/// positions in a hundred no sample varied at are never built at all.
///
/// **The arena is one buffer per sample, not a box per record.** That is the measurement's
/// doing rather than tidiness: a skipping walk's speed is the speed of a walk that allocates
/// nothing, and a box per record would put an allocation back on every record in the window to
/// save one on the eighth of them that gets built (§3.4 there).
///
/// **What it does not do is decide.** A per-sample predicate cannot: a sample showing no
/// non-reference read at a position can never be the one that admits the locus, but its
/// evidence is still needed at every locus another sample admits, because thirty reference
/// reads and no reads at all are different genotypes. So this keeps every body it passes and
/// lets the cohort choose.
pub struct PspSummarySource<'a> {
    walk: RecordIter<'a>,
    kept: Vec<u8>,
    /// The individual this file holds, so a failure names them.
    sample: String,
    /// How far the walk has got — the other half of locating a failure.
    reached: WalkProgress,
    /// This file's own record layout, kept so a body can be built later without the walk.
    layout: RecordLayout,
    /// Walk-local read-group identifier `i` is this run's `read_groups[i]` (spec §6.2).
    read_groups: Vec<ReadGroupId>,
    /// The heads handed over so far, so a body can be located and described when it is built.
    heads: Vec<KeptRecord>,
    /// What this sample contributed, for the run report.
    read: StoredSampleTallies,
}

/// One stored record, summarised, with its body waiting in the source's arena.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KeptRecord {
    /// What the head said — everything the merge decides on before it assembles a locus.
    pub summary: LocusSummary,
    /// Where the body's bytes sit in the source's arena.
    pub body: core::ops::Range<usize>,
    /// **The reads live at this record, which its body cannot be built without.**
    ///
    /// One observation in a record has its read list *derived* rather than stored — the
    /// residual is this set minus every other observation's list, and the body declares what
    /// that should come to. So a body is not self-contained after all: built against the wrong
    /// set it is refused, which is how this field came to exist.
    ///
    /// **Kept per record, which is the simple answer and not yet the measured one.** The set is
    /// the reads open at one position, so it grows with depth: at three reads a position it is
    /// a handful of identifiers and at three hundred it is three hundred, over a window of a
    /// few thousand records a sample. The alternatives both trade memory for work — a snapshot
    /// every so many records with the changes replayed from it, or one per building region —
    /// and choosing between them wants the measurement that has not been taken
    /// (`spec/cohort_merge_psp_path.md` §3.2).
    pub live: LiveSet,
}

impl<'a> PspSummarySource<'a> {
    /// Walk `psp` from its first record, keeping every body.
    ///
    /// # Errors
    ///
    /// Whatever opening the walk refuses — a psp whose first block will not read.
    pub fn over(psp: &'a mut PspReader, read_groups: &[ReadGroupId]) -> Result<Self, RunError> {
        let sample = psp.header().sample.clone();
        let layout = RecordLayout::from_manifest(&psp.header().manifest)
            .map_err(|source| source_failed(&sample, WalkProgress::NothingYet, source))?;
        let read_groups = read_groups.to_vec();
        let walk = psp
            .records()
            .map_err(|source| source_failed(&sample, WalkProgress::NothingYet, source))?;
        Ok(Self {
            walk,
            kept: Vec::new(),
            sample,
            reached: WalkProgress::NothingYet,
            layout,
            read_groups,
            heads: Vec::new(),
            read: StoredSampleTallies::default(),
        })
    }

    /// A failure of this sample's walk, named and located.
    fn refuse(&self, cause: impl std::error::Error + Send + Sync + 'static) -> RunError {
        source_failed(&self.sample, self.reached, cause)
    }

    /// The individual this file holds.
    #[must_use]
    pub fn sample_name(&self) -> &str {
        &self.sample
    }

    /// What this sample contributed to the run — **counted from the heads**, which is where
    /// both numbers live, so a run that never builds a body still reports them.
    #[must_use]
    pub fn read(&self) -> StoredSampleTallies {
        self.read
    }

    /// The bytes kept so far. A body's range addresses this.
    #[must_use]
    pub fn kept(&self) -> &[u8] {
        &self.kept
    }

    /// The next record's summary, its body appended to the arena.
    ///
    /// # Errors
    ///
    /// Whatever the walk refuses — a damaged block, or a record that will not parse.
    pub fn next_summary(&mut self) -> Option<Result<KeptRecord, PspReadError>> {
        match self.walk.next_keeping(&mut self.kept) {
            Some(Ok(streamed)) => {
                let body = streamed
                    .body
                    .expect("a walk that builds nothing keeps every body");
                Some(Ok(KeptRecord {
                    summary: LocusSummary::from(&streamed.head),
                    body,
                    // **After the step, so the changes this record's head carried are in.**
                    // The walk parses a head's live-set changes and applies them before
                    // handing the record over, so this is the set as of this record and not
                    // as of the one before it.
                    live: self.walk.live_reads().clone(),
                }))
            }
            Some(Err(failed)) => Some(Err(failed)),
            None => None,
        }
    }
}

pub struct PspObservationSource<W> {
    /// The individual this file holds. **[`over`](Self::over) takes it from the psp's own
    /// header**, so on the path a run uses, the name a failure carries and the file it came
    /// from cannot come apart. [`new`](Self::new) takes whatever it is given, which is why it
    /// is not part of this module's public surface.
    sample: String,
    /// How far the source has got — the second half of locating a failure (spec §9).
    reached: WalkProgress,
    /// Where the last record handed over started, which is what the order check compares
    /// against. **Not [`reached`](Self::reached)**, which is where the last record *ended*:
    /// records may overlap — a deletion covers the bases a later record starts on — so
    /// comparing starts against ends would refuse sound files.
    last_start: Option<GenomePosition>,
    /// Set once this source has refused a record of its own — see
    /// [`PspSourceError::AlreadyRefused`] for what it stops.
    refused: bool,
    /// What this sample's walk-local read groups are called in this run: entry `i` is the
    /// run-wide identifier of the group its own walk numbered `i` (spec §6.2).
    ///
    /// **Every sample numbers its read groups from zero**, so without this every sample's
    /// first group would reach the merge as identifier 0 and the cohort would score them all
    /// against one calibration. A single-sample run's map is the identity and the work is
    /// wasted there; it is applied unconditionally anyway, because a source that sometimes
    /// renumbers is a source somebody forgets to hand a map to.
    read_groups: Vec<ReadGroupId>,
    /// **What this source drew, kept past the draw** — the only per-sample facts psp mode's
    /// run report has, because a stored file records no walk.
    read: StoredSampleTallies,
    walk: W,
}

/// **What one sample's stored file gave this run**, counted as it was decoded.
///
/// **Measured by the run that read it, not copied from the file.** A psp carries no count of
/// what its walk kept or dropped — those tallies belong to the read cursor of a walk that
/// happened in another process — so a run over stored files can only say what it drew, and
/// these two are that. They are the psp-mode counterpart of direct mode's
/// [`SampleWalkTallies`](super::callers::SampleWalkTallies), and they are deliberately not the
/// same fields: a number that cannot be measured here is absent rather than zero.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct StoredSampleTallies {
    /// How many stored loci this run drew out of the file.
    ///
    /// **Loci and not observations.** One psp record is one locus with every read that spoke
    /// at it inside, so this is several times smaller than the observations it carries — the
    /// same distinction `generate-psps` draws when it reports what it stored.
    pub loci_read: u64,
    /// **Summed over those loci, how many reads each compared with the reference** — the
    /// record head's `reads-compared-with-reference`, the keep rule's own denominator.
    ///
    /// **Not the sample's depth, and divided by [`loci_read`](Self::loci_read) not its mean
    /// depth either.** The head's own field says what it leaves out
    /// ([`RecordHead::reads_compared_with_reference`](crate::ng::psp::RecordHead::reads_compared_with_reference)):
    /// reads a filter turned away, reads the per-position cap discarded, reads that covered the
    /// locus and produced no observation, and reads whose witness stopped inside it. At a repeat
    /// tract 40 reads cover but only 22 anchor both borders of, this counts 22.
    ///
    /// It is worth carrying because it is the number the keep rule was applied with, and it is
    /// the closest a stored file comes to direct mode's *what this sample's reads did* — which
    /// no psp records at all.
    pub reads_compared_with_reference: u64,
}

impl StoredSampleTallies {
    /// **How many reads went into the comparison at one locus, on average** — or `None` where
    /// this source read none. Not depth: see
    /// [`reads_compared_with_reference`](Self::reads_compared_with_reference).
    ///
    /// `None` rather than zero, because the two are different facts about a sample: a file
    /// with no loci over this ground contributed nothing, where a mean of zero would say every
    /// locus it did hold was compared against no reads at all.
    #[must_use]
    pub fn mean_reads_a_locus(&self) -> Option<f64> {
        (self.loci_read > 0)
            .then(|| self.reads_compared_with_reference as f64 / self.loci_read as f64)
    }
}

/// **The sample's name and how far it has got, not the walk.** The sibling walker hand-writes
/// its own for the same reason ([`AlignmentFilesWalker`](super::walker::AlignmentFilesWalker)):
/// a derived one would print the decoder's buffers. It also drops the `W: Debug` bound a derive
/// would add, which a source over a mapped or filtered walk would otherwise fail.
impl<W> std::fmt::Debug for PspObservationSource<W> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("PspObservationSource")
            .field("sample", &self.sample)
            .field("reached", &self.reached)
            .field("refused", &self.refused)
            .finish_non_exhaustive()
    }
}

impl<W> PspObservationSource<W> {
    /// A source over `walk`, reporting failures under `sample`.
    ///
    /// **The general constructor, and the one a fixture uses**: `walk` is any iterator of the
    /// store's streamed records, so a test can hand over a `Vec`'s iterator where a run hands
    /// over [`PspReader::records`]. Named for its sibling: `AlignmentFilesWalker::new` is
    /// direct mode's general constructor and `::over` is what a run builds, and these two are
    /// the same pair.
    ///
    /// **Crate-private, and that is what makes the sample name trustworthy.** The name is
    /// whatever the caller passes, so a public one would let a failure name an individual the
    /// file does not hold — and at a cohort of thousands the whole value of
    /// [`RunError::SourceFailed`] is that its name identifies the file.
    pub(crate) fn new(sample: String, walk: W, read_groups: &[ReadGroupId]) -> Self {
        Self {
            sample,
            reached: WalkProgress::NothingYet,
            last_start: None,
            refused: false,
            read_groups: read_groups.to_vec(),
            read: StoredSampleTallies::default(),
            walk,
        }
    }

    /// The individual this source reads.
    #[must_use]
    pub fn sample_name(&self) -> &str {
        &self.sample
    }

    /// **What this source drew out of its file**, which is what a run over stored files has to
    /// say about the sample (plan step F1).
    #[must_use]
    pub fn read(&self) -> StoredSampleTallies {
        self.read
    }

    /// How far this source has got. `NothingYet` until the first observation is handed over,
    /// and the last base of the last one after that — **not** where the reader has decoded to,
    /// which is ahead of it and is the block stream's business.
    #[must_use]
    pub fn reached(&self) -> WalkProgress {
        self.reached
    }

    /// Fail this source, naming the sample and how far it had got.
    ///
    /// **It does not latch on its own.** The walk's own failures fuse the walk and need no
    /// flag; the two this type mints itself do, and they set it at the call site so that the
    /// two kinds of failure are visibly different in the code as well as in the doc.
    fn refuse(&self, cause: impl std::error::Error + Send + Sync + 'static) -> RunError {
        source_failed(&self.sample, self.reached, cause)
    }
}

/// **The one place a psp source's failure is shaped**, so that a failure before the source
/// exists — [`PspObservationSource::over`], where the walk itself would not start — and one
/// from a draw cannot render differently.
fn source_failed(
    sample: &str,
    reached: WalkProgress,
    cause: impl std::error::Error + Send + Sync + 'static,
) -> RunError {
    RunError::SourceFailed {
        sample: sample.to_string(),
        reached,
        source: Box::new(cause),
    }
}

impl<'a> PspObservationSource<RecordIter<'a>> {
    /// Every record of an open psp, as a source named by the file's own header.
    ///
    /// **This is what a calling run builds**, one per open psp, and it is where the file's
    /// sample name becomes the name every failure of this source carries. Its sibling is
    /// [`AlignmentFilesWalker::over`](super::walker::AlignmentFilesWalker::over), direct
    /// mode's run-shaped constructor.
    ///
    /// # Errors
    ///
    /// [`RunError::SourceFailed`] if the walk cannot be started — the file's manifest declares
    /// an encoding this build cannot read, or the first seek fails. It reports
    /// [`WalkProgress::NothingYet`], which is exact: nothing has been decoded.
    pub fn over(psp: &'a mut PspReader, read_groups: &[ReadGroupId]) -> Result<Self, RunError> {
        // **Cloned before the walk is asked for**, and not merely to keep the borrow checker
        // happy: the walk borrows the reader for as long as this source lives, so the header is
        // out of reach afterwards.
        let sample = psp.header().sample.clone();
        let walk = psp
            .records()
            // **`NothingYet` is exact rather than a placeholder**: the source does not exist
            // yet, so nothing has been decoded.
            //
            // ⚠ **This arm has no test, and no fixture can reach it through a `PspReader`.**
            // `records` fails on a seek this file's offsets have already been bounded inside
            // (`psp::walk`, whose own note records the same gap), on a manifest `open` has
            // already parsed, or on a buffer ceiling `with_a_record_buffer_ceiling` refuses
            // where it is set. Measured: replacing `&sample` here with a literal leaves all
            // fifteen tests of this module green. It is written out rather than left to `?` so
            // that the sample is named the day something does reach it — but do not read it as
            // covered.
            .map_err(|source| source_failed(&sample, WalkProgress::NothingYet, source))?;
        Ok(Self::new(sample, walk, read_groups))
    }
}


/// **Reading a stored sample the way a cohort run wants it read.**
///
/// Every draw is a [`Drawn::Kept`]: the head's summary, and the body left in this source's
/// arena. The cohort decides from the summaries and asks for evidence only where a locus
/// survived, which at about one position in a hundred is what makes the other ninety-nine free
/// (`spec/cohort_merge_psp_path.md` §3.1).
impl ObservationSource for PspSummarySource<'_> {
    type Error = RunError;

    /// **Kept and then built at once** — the shape a caller wanting every record takes.
    ///
    /// It exists because the trait requires it and because a source must still be usable by
    /// something that wants everything; the cache calls [`next_drawn`](Self::next_drawn)
    /// instead, which is the whole point of this type.
    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, RunError>> {
        match self.next_drawn(spare)? {
            Ok(Drawn::Kept { body, .. }) => Some(self.build(body)),
            Ok(Drawn::Built(record)) => Some(Ok(record)),
            Err(failed) => Some(Err(failed)),
        }
    }

    fn next_drawn(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<Drawn, RunError>> {
        drop(spare);
        let kept = match self.next_summary()? {
            Ok(kept) => kept,
            Err(failed) => return Some(Err(self.refuse(failed))),
        };
        self.reached = WalkProgress::After(kept.summary.reach_position());
        self.read.loci_read += 1;
        self.read.reads_compared_with_reference +=
            u64::from(kept.summary.reads_compared_with_reference);
        self.heads.push(kept.clone());
        Some(Ok(Drawn::Kept {
            summary: kept.summary,
            body: kept.body,
        }))
    }

    /// Build the body kept at `body`.
    ///
    /// **⚠ Decoded against an empty live set, which is right today and will not always be.**
    /// A body needs the set of reads live at its own record only to derive the one observation
    /// whose read list is stored as a residual; the encoder does not yet write chain ids at all
    /// (`psp/record.rs`, `encode_record_body`), so every file this caller produces has an empty
    /// set everywhere and the residual is trivial. When the encoding's Milestone E writes them,
    /// the set as of *this* record has to reach here — and since the head walk is what carries
    /// it, that means keeping it beside the bytes or replaying the heads to reach it. Building
    /// against the wrong set is the failure that does not announce itself: a body decoded
    /// against a plausible-but-wrong set is a plausible body.
    fn build(&self, body: core::ops::Range<usize>) -> Result<SampleLocusObservations, RunError> {
        let kept = self
            .heads
            .binary_search_by_key(&body.start, |kept| kept.body.start)
            .map(|at| &self.heads[at])
            .expect("a body range this source handed out");
        let head = RecordHead {
            region: kept.summary.region,
            non_reference_reads: kept.summary.non_reference_reads,
            reads_compared_with_reference: kept.summary.reads_compared_with_reference,
            body_bytes: u32::try_from(body.len()).expect("a body this source wrote down"),
        };
        let found = LocatedRecord {
            head,
            body: &self.kept[body.clone()],
            record_bytes: body.len(),
        };
        let mut record = decode_the_body_of(&found, &kept.live, &self.layout)
            .map_err(|source| self.refuse(source))?
            .record;
        // **The same renumbering the building source makes, for the same reason**: every
        // sample numbers its read groups from zero, so without it every sample's first group
        // would reach the merge as identifier 0 and the cohort would score them as one lane.
        for observation in &mut record.observations {
            let walk_local = observation.read_group.get() as usize;
            let Some(run_wide) = self.read_groups.get(walk_local).copied() else {
                return Err(self.refuse(PspSourceError::ReadGroupNotInThisFilesTable {
                    at: record.region,
                    names: observation.read_group.get(),
                    in_the_table: self.read_groups.len(),
                }));
            };
            observation.read_group = run_wide;
        }
        Ok(record)
    }
}

impl<W> ObservationSource for PspObservationSource<W>
where
    W: Iterator<Item = Result<StreamedRecord, PspReadError>>,
{
    type Error = RunError;

    /// The next observation stored for this sample, or `None` once the file is spent.
    ///
    /// **The spare is dropped**, which the trait permits — see the module documentation for
    /// why taking it is a change to the store's walk rather than to this type. It is released
    /// at the top rather than left to fall out of scope, so its buffers are back with the
    /// allocator before the next record asks for any.
    ///
    /// **Exhaustion is final**, which is the walk's guarantee rather than a fresh one: a
    /// [`RecordIter`] that fails yields its `Err` once and then `None`, and one that reaches
    /// the end of the blocks answers `None` for ever.
    ///
    /// **⚑ A failed source is not left live, and the trait's contract says a failure should
    /// leave it so.** This deviates on all three of its failure paths, in two different ways,
    /// and neither is reachable today — the merge abandons its cache at the first failure
    /// (`ObservationCache::draw_next` propagates without marking the source spent, and both
    /// drivers drop the cache), so nothing ever draws again. Written down because what a retry
    /// would meet is not obvious from the code:
    ///
    /// - **a failure inside the walk fuses the walk**, so the next draw answers `None` — the
    ///   walker's deviation exactly, and the silent one: a consumer that swallowed the error
    ///   would read `None` as exhaustion and build cohort loci without this sample;
    /// - **the two refusals this type mints latch**, and every later draw returns
    ///   [`PspSourceError::AlreadyRefused`] rather than `None` or the next record. Latching is
    ///   what stops the refused observation being dropped from the cohort in silence.
    ///
    /// Anything that adds a retry has to settle both.
    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, RunError>> {
        drop(spare);
        if self.refused {
            return Some(Err(self.refuse(PspSourceError::AlreadyRefused)));
        }
        // **Destructured with no `..`**, so a field added to `StreamedRecord` has to be
        // considered here rather than dropped silently — the same rule `psp::reader` states
        // where it reads the walk's own shape. Step E2's read-group remap lands in this
        // function, and it is a second reason to care what a streamed record carries.
        let StreamedRecord {
            block: _,
            // This adapter builds every body, so nothing is ever kept for later.
            body: _,
            head,
            record,
        } = match self.walk.next()? {
            Ok(streamed) => streamed,
            Err(source) => return Some(Err(self.refuse(source))),
        };
        let Some(mut record) = record else {
            self.refused = true;
            return Some(Err(self.refuse(PspSourceError::ObservationBodyNotBuilt {
                at: head.region,
            })));
        };
        // **`start_position` and not the head's region**, so that this compares the value the
        // merge's own cache orders on rather than a second reading of the same fact.
        let start = record.start_position();
        if let Some(previous) = self.last_start
            && start < previous
        {
            self.refused = true;
            return Some(Err(self.refuse(PspSourceError::ObservationsOutOfOrder {
                previous,
                offered: start,
            })));
        }
        // **Renumbered here, at the one place a record crosses from the file's terms into the
        // run's.** Doing it further up — in the merge, or at the call — would mean every
        // consumer knowing which sample a record came from in order to read its read group,
        // which is exactly what the run-wide numbering exists to remove.
        for observation in &mut record.observations {
            let walk_local = observation.read_group.get() as usize;
            let Some(run_wide) = self.read_groups.get(walk_local).copied() else {
                self.refused = true;
                return Some(Err(self.refuse(
                    PspSourceError::ReadGroupNotInThisFilesTable {
                        at: record.region,
                        names: observation.read_group.get(),
                        in_the_table: self.read_groups.len(),
                    },
                )));
            };
            observation.read_group = run_wide;
        }
        self.last_start = Some(start);
        // **Counted where the record is handed over, not where it is decoded.** A record the
        // source refuses is not a record this run read, and the two differ by exactly the
        // refused one — which is the number the report would otherwise over-state.
        self.read.loci_read += 1;
        self.read.reads_compared_with_reference += u64::from(head.reads_compared_with_reference);
        // **`reach_position` rather than `region.end`**: `reach` is `end.max(start)`, and
        // `GenomeRegion` has public fields and no constructor, so an inverted region read
        // straight off `end` would put an observation's reach before its own first base. The
        // walker and the merge's cache both key on this same call.
        self.reached = WalkProgress::After(record.reach_position());
        Some(Ok(record))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{
        LocusKind, ReadWitness, SsrDetail, WitnessedLocusPositions,
    };
    use crate::ng::psp::PspWriter;
    use crate::ng::psp::writer::tests_support::{a_header, a_record, a_sample};
    use crate::ng::types::{ContigId, Motif, Position};
    use std::path::{Path, PathBuf};

    /// **The map a single-sample run has: the identity.** The fixture header declares two read
    /// groups, so a source over one of these file renumbers 0 to 0 and 1 to 1 — which is what
    /// leaves every other test in this module about what it was about before the run-wide
    /// numbering existed.
    fn as_walked() -> Vec<ReadGroupId> {
        vec![ReadGroupId(0), ReadGroupId(1)]
    }

    /// A psp holding `records`, and the directory it lives in — dropped by the caller, which
    /// is what deletes the file.
    ///
    /// The 1 kb grid is `a_sample`'s: it cuts the twenty records of each contig into four
    /// blocks, so a walk over one of these crosses block boundaries rather than reading a
    /// single block and stopping. Take the grid as a parameter when a test needs a second one.
    fn a_psp_of(records: &[SampleLocusObservations]) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let path = dir.path().join("sample.psp");
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("the header writes");
        for record in records {
            writer.push(record).expect("the fixtures are in order");
        }
        let _ = writer.finish(b"").expect("the file seals");
        (dir, path)
    }

    /// Everything a source hands over, and the failure that ended it if one did.
    fn drain(
        source: &mut impl ObservationSource<Error = RunError>,
    ) -> (Vec<SampleLocusObservations>, Option<RunError>) {
        let mut observations = Vec::new();
        while let Some(next) = source.next_observation(None) {
            match next {
                Ok(observation) => observations.push(observation),
                Err(failed) => return (observations, Some(failed)),
            }
        }
        (observations, None)
    }

    /// The whole point of the adapter: what the walk stage stored is what the merge is handed.
    ///
    /// **Field for field and in order**, over a file whose records cross eight blocks and two
    /// contigs — a source that read one block and stopped, or that reordered on a contig
    /// change, fails here rather than at a cohort.
    #[test]
    fn the_records_a_psp_holds_come_back_as_the_observations_that_were_stored() {
        let written = a_sample();
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert!(
            psp.block_index().len() > 1,
            "the fixture is meant to cross blocks; it holds {}",
            psp.block_index().len()
        );
        let mut source =
            PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        let (read, failed) = drain(&mut source);
        assert!(failed.is_none(), "the file is sound: {failed:?}");
        assert_eq!(read, written);
    }

    /// **A repeat tract carries a motif and two flanks the generic path has no field for**, and
    /// the merge needs them: the tract model is what reads them. So the kind travels through
    /// the source whole rather than as *some locus was here*.
    #[test]
    fn a_repeat_tract_reaches_the_merge_with_its_motif_and_both_flanks() {
        let mut tract = a_record(0, 41, 6);
        tract.kind = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide is a motif"),
            // Different lengths, because a left flank is clamped at a contig's start: a source
            // that swapped them would pass if they matched.
            left_flank: b"GG".to_vec().into_boxed_slice(),
            right_flank: b"CCCCCCCCCC".to_vec().into_boxed_slice(),
        });
        let written = vec![a_record(0, 1, 1), tract, a_record(0, 101, 1)];
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let mut source =
            PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        let (read, failed) = drain(&mut source);
        assert!(failed.is_none(), "the file is sound: {failed:?}");
        assert_eq!(read, written);
    }

    /// **The sample a failure names is the file's own**, not something a caller passed beside
    /// it — which is what stops a cohort's error naming the wrong individual after a misplaced
    /// argument.
    #[test]
    fn the_sample_a_source_reports_under_is_the_one_its_header_names() {
        let (_dir, path) = a_psp_of(&a_sample());
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let in_the_header = psp.header().sample.clone();
        let source = PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        assert_eq!(source.sample_name(), in_the_header);
    }

    /// **Analysed and empty is a file, not an absence** (spec §8): a sample whose ground held
    /// no observations writes a psp with no blocks, and the merge must read it as a sample that
    /// saw nothing rather than refusing it.
    #[test]
    fn a_psp_holding_no_records_is_a_source_that_is_spent_from_the_start() {
        let (_dir, path) = a_psp_of(&[]);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert!(psp.block_index().is_empty(), "no records, so no blocks");
        let mut source =
            PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        assert!(source.next_observation(None).is_none());
        assert_eq!(source.reached(), WalkProgress::NothingYet);
    }

    /// **A walk that declined a body is refused, not handed over short.**
    ///
    /// The predicate here keeps only records with a non-reference read, and the fixture's
    /// observations all match the reference — so every record arrives as a head alone. Fed to
    /// the merge, that would be a cohort locus with one sample's evidence silently missing;
    /// what it must be instead is an error naming the sample and the locus.
    #[test]
    fn a_walk_that_skipped_a_body_is_refused_rather_than_passed_on_without_its_evidence() {
        let written = a_sample();
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let sample = psp.header().sample.clone();
        let walk = psp
            .records_where(|head| head.non_reference_reads > 0)
            .expect("the walk starts");
        let mut source = PspObservationSource::new(sample.clone(), walk, &as_walked());
        let (read, failed) = drain(&mut source);
        assert!(read.is_empty(), "no record had a body to hand over");
        let Some(RunError::SourceFailed {
            sample: named,
            reached,
            source: cause,
        }) = failed
        else {
            panic!("a head-only record must fail the source: {failed:?}");
        };
        assert_eq!(named, sample);
        assert_eq!(reached, WalkProgress::NothingYet, "nothing was handed over");
        let cause = cause
            .downcast::<PspSourceError>()
            .expect("the cause says what was wrong with the file");
        let PspSourceError::ObservationBodyNotBuilt { at } = *cause else {
            panic!("the body is what was missing: {cause:?}");
        };
        assert_eq!(
            at, written[0].region,
            "the first record is where it stopped"
        );
        // **A record whose body was never built is not a record this run read.** The counters
        // sit below this refusal, so nothing is counted here; moving them above it would put
        // loci in a run report that no genotype could ever be built from.
        assert_eq!(
            source.read(),
            StoredSampleTallies::default(),
            "nothing was handed over, so nothing is counted: {:?}",
            source.read(),
        );
    }

    /// **Backwards observations are an error, not an assertion** (arch §8).
    ///
    /// The merge's own cache aborts the process on a source that goes backwards, on the ground
    /// that its generators cannot mint one — a ground that does not hold for a file written by
    /// another process. This is the fixture that proves the psp path refuses first: the records
    /// are the ones a sound psp holds, handed over in the wrong order.
    #[test]
    fn stored_observations_that_go_backwards_are_refused_naming_the_sample() {
        let written = a_sample();
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        // Read the file honestly first, then hand the records back in the wrong order: the
        // store refuses an out-of-order `push`, so a psp on disk cannot hold one.
        let sound: Vec<_> = psp
            .records()
            .expect("the walk starts")
            .collect::<Result<Vec<_>, _>>()
            .expect("the file is sound");
        let mut backwards = vec![Ok(sound[1].clone()), Ok(sound[0].clone())];
        backwards.extend(sound.into_iter().skip(2).map(Ok));
        let mut source = PspObservationSource::new(
            "SRR7279481".to_string(),
            backwards.into_iter(),
            &as_walked(),
        );
        let (read, failed) = drain(&mut source);
        assert_eq!(
            read,
            vec![written[1].clone()],
            "the record that went backwards is refused, the one in front of it is not"
        );
        let Some(RunError::SourceFailed {
            sample,
            reached,
            source: cause,
        }) = failed
        else {
            panic!("going backwards must fail the source: {failed:?}");
        };
        assert_eq!(sample, "SRR7279481");
        assert_eq!(
            reached,
            WalkProgress::After(written[1].reach_position()),
            "how far it got is the record before the refusal, not the refused one"
        );
        let cause = cause
            .downcast::<PspSourceError>()
            .expect("the cause says what was wrong with the file");
        let PspSourceError::ObservationsOutOfOrder { previous, offered } = *cause else {
            panic!("the order is what was wrong: {cause:?}");
        };
        assert_eq!(previous, written[1].start_position());
        assert_eq!(offered, written[0].start_position());
        // **A refused record is not a record this run read**, and the counters are incremented
        // where the record is handed over rather than where it is decoded, so exactly one is
        // counted here. Counting at the decode instead would say two, and a run report would
        // then credit this sample with a locus no genotype was ever built from.
        assert_eq!(
            source.read().loci_read,
            1,
            "one record was handed over and one was refused: {:?}",
            source.read(),
        );
    }

    /// **A source that refused a record is finished, and says so on every later draw.**
    ///
    /// Without the latch the next draw hands back the record *after* the refused one and goes
    /// on succeeding, so a consumer that swallowed the error would be handed a stream that
    /// looks sound and is one observation short — one sample missing at one locus of the
    /// cohort, with nothing left to say so.
    #[test]
    fn a_source_that_refused_a_record_refuses_every_later_draw() {
        let written = a_sample();
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let sound: Vec<_> = psp
            .records()
            .expect("the walk starts")
            .collect::<Result<Vec<_>, _>>()
            .expect("the file is sound");
        let mut backwards = vec![Ok(sound[1].clone()), Ok(sound[0].clone())];
        backwards.extend(sound.into_iter().skip(2).map(Ok));
        let mut source = PspObservationSource::new(
            "SRR7279481".to_string(),
            backwards.into_iter(),
            &as_walked(),
        );
        let _ = source
            .next_observation(None)
            .expect("a record")
            .expect("the first is in order");
        let _ = source
            .next_observation(None)
            .expect("a draw")
            .expect_err("the second goes backwards");

        for draw in 0..2 {
            let failed = source
                .next_observation(None)
                .unwrap_or_else(|| panic!("draw {draw} must not be silence"))
                .expect_err("a refused source hands over nothing more");
            let RunError::SourceFailed { source: cause, .. } = failed else {
                panic!("still a source failure: {failed:?}");
            };
            let cause = cause
                .downcast::<PspSourceError>()
                .expect("the cause is this source's own");
            assert!(
                matches!(*cause, PspSourceError::AlreadyRefused),
                "draw {draw}: {cause:?}"
            );
        }
    }

    /// **`reached` is the last base the last observation covered, not the base it began on.**
    ///
    /// Every other fixture here is built from one-base records, where the two are the same
    /// number — so this is the only test that can tell them apart, and the difference is the
    /// width of a deletion or a repeat tract: exactly the loci an operator reading a failure
    /// would go and look at.
    #[test]
    fn reached_is_the_last_base_the_last_observation_covered() {
        // Six bases wide: it starts at 11 and reaches 16.
        let written = vec![a_record(0, 11, 6)];
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let mut source =
            PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        let _ = source
            .next_observation(None)
            .expect("the one record")
            .expect("the file is sound");
        assert_eq!(
            source.reached(),
            WalkProgress::After(GenomePosition {
                contig: ContigId(0),
                position: Position(16),
            })
        );
    }

    /// **A spent source answers `None` for ever** — the trait's "exhaustion is final", which
    /// the merge's cache relies on and which this type inherits from the walk rather than
    /// enforcing itself. A source that yielded `Some` after a `None` would be drawn in behind
    /// the cache's own window and so silently out of coordinate order.
    #[test]
    fn a_spent_source_answers_none_for_ever() {
        let written = a_sample();
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let mut source =
            PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        let (read, failed) = drain(&mut source);
        assert!(failed.is_none(), "the file is sound: {failed:?}");
        assert_eq!(read.len(), written.len());
        assert!(source.next_observation(None).is_none());
        assert!(source.next_observation(None).is_none());
    }

    /// **A body declined part-way through reports the ground already handed over.** The
    /// head-only test above refuses on the very first record, where `reached` reads
    /// `NothingYet` whether it is right or has never been set.
    #[test]
    fn a_body_declined_partway_through_reports_the_records_already_handed_over() {
        let written = a_sample();
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let sample = psp.header().sample.clone();
        // `a_sample`'s first four records start at 1, 101, 201 and 301, so this builds three
        // bodies and declines the fourth.
        let walk = psp
            .records_where(|head| head.region.start.get() < 300)
            .expect("the walk starts");
        let mut source = PspObservationSource::new(sample, walk, &as_walked());
        let (read, failed) = drain(&mut source);
        assert_eq!(read, written[..3]);
        let Some(RunError::SourceFailed {
            reached,
            source: cause,
            ..
        }) = failed
        else {
            panic!("a declined body must fail the source: {failed:?}");
        };
        assert_eq!(reached, WalkProgress::After(written[2].reach_position()));
        let cause = cause
            .downcast::<PspSourceError>()
            .expect("the cause is this source's own");
        let PspSourceError::ObservationBodyNotBuilt { at } = *cause else {
            panic!("the body is what was missing: {cause:?}");
        };
        assert_eq!(at, written[3].region);
    }

    /// **Two records starting on the same base are in order**, and the check must not refuse
    /// them: the writer admits them (a block's byte ceiling can close after a record whose
    /// successor starts on the same base), and a deletion makes records overlap by
    /// construction.
    #[test]
    fn records_that_start_on_the_same_base_are_not_backwards() {
        // A deletion-shaped first record covering six bases, then two records starting inside
        // it — the overlap a comparison against the previous record's *end* would refuse.
        let written = vec![a_record(0, 11, 6), a_record(0, 13, 1), a_record(0, 13, 1)];
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let mut source =
            PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        let (read, failed) = drain(&mut source);
        assert!(
            failed.is_none(),
            "these are in coordinate order: {failed:?}"
        );
        assert_eq!(read, written);
    }

    /// **A read failure names the sample and how far it had got**, which is the whole of what
    /// [`RunError::SourceFailed`] exists for: at a cohort of thousands, the file that failed
    /// and the ground already called are two different questions and an operator has both.
    #[test]
    fn a_failing_walk_reports_the_sample_and_the_last_observation_it_handed_over() {
        let written = a_sample();
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let sound: Vec<_> = psp
            .records()
            .expect("the walk starts")
            .collect::<Result<Vec<_>, _>>()
            .expect("the file is sound");
        let mut walk: Vec<Result<StreamedRecord, PspReadError>> =
            sound.into_iter().take(3).map(Ok).collect();
        walk.push(Err(PspReadError::Incomplete { path: path.clone() }));
        let mut source =
            PspObservationSource::new("SRR7279481".to_string(), walk.into_iter(), &as_walked());
        let (read, failed) = drain(&mut source);
        assert_eq!(read, written[..3]);
        let Some(RunError::SourceFailed {
            sample, reached, ..
        }) = failed
        else {
            panic!("a walk failure must reach the run as a source failure: {failed:?}");
        };
        assert_eq!(sample, "SRR7279481");
        assert_eq!(reached, WalkProgress::After(written[2].reach_position()));
    }

    /// The same, through a **damaged file** rather than a hand-made failure: the reader's own
    /// refusal is what a run meets, and this is what proves it arrives as a run error rather
    /// than as a panic or a short read.
    #[test]
    fn a_damaged_block_ends_the_source_with_the_ground_it_had_already_covered() {
        let written = a_sample();
        let (dir, path) = a_psp_of(&written);
        let block_at = {
            let psp = PspReader::open(&path).expect("a finished psp opens");
            // The third block, so records from the first two are handed over before the walk
            // meets the damage — which is what makes `reached` say something.
            let Some(third) = psp.block_index().get(2) else {
                panic!(
                    "the fixture is meant to hold at least three blocks; it holds {}",
                    psp.block_index().len()
                );
            };
            third.block_offset as usize
        };
        let mut bytes = std::fs::read(&path).expect("the file reads");
        // **Inside the compressed frame, past the four-byte length and past zstd's own frame
        // header**, so the block's declared length still points where it did and the damage is
        // met by the decompressor rather than by the framing.
        let declared = u32::from_le_bytes(
            bytes[block_at..block_at + 4]
                .try_into()
                .expect("four bytes of length"),
        ) as usize;
        let frame_at = block_at + 4;
        let damage_from = frame_at + declared / 2;
        for byte in &mut bytes[damage_from..frame_at + declared] {
            *byte ^= 0xff;
        }
        let damaged = dir.path().join("damaged.psp");
        std::fs::write(&damaged, &bytes).expect("the damaged copy writes");

        let mut psp = PspReader::open(&damaged).expect("the footer and index are untouched");
        let mut source =
            PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        let (read, failed) = drain(&mut source);
        assert!(
            !read.is_empty() && read.len() < written.len(),
            "the walk must get some way in and then stop, not all of it and not none: it \
             handed over {} of {}",
            read.len(),
            written.len()
        );
        assert_eq!(
            read,
            written[..read.len()],
            "what it did hand over is what was stored"
        );
        let Some(RunError::SourceFailed {
            sample, reached, ..
        }) = failed
        else {
            panic!("a damaged block must reach the run as a source failure: {failed:?}");
        };
        assert_eq!(sample, "SRR7279481");
        assert_eq!(
            reached,
            WalkProgress::After(written[read.len() - 1].reach_position()),
            "how far it got is the last observation it handed over"
        );
    }

    /// **A source that fails on its very first draw reports that nothing was reached**, rather
    /// than inventing a position: a run that said "failed at contig 0 position 1" when nothing
    /// had been decoded would send an operator to an innocent locus.
    #[test]
    fn a_source_that_fails_on_its_very_first_draw_reports_nothing_reached() {
        let failure = PspReadError::damaged(Path::new("sample.psp"), "a manifest".to_string());
        let mut source = PspObservationSource::new(
            "SRR7279481".to_string(),
            std::iter::once(Err::<StreamedRecord, _>(failure)),
            &as_walked(),
        );
        let Some(Err(RunError::SourceFailed {
            sample, reached, ..
        })) = source.next_observation(None)
        else {
            panic!("the first draw fails");
        };
        assert_eq!(sample, "SRR7279481");
        assert_eq!(reached, WalkProgress::NothingYet);
    }

    /// **Two samples whose read groups collide, renumbered apart.**
    ///
    /// This is the case the run-wide numbering exists for and the one a psp cannot avoid: every
    /// walk sees one sample and numbers its groups from zero, so both files here call their two
    /// groups 0 and 1. Handed to the merge unrenumbered, every sample's first group would be
    /// identifier 0 and the cohort would score four libraries against two calibrations —
    /// silently, because the numbers are all in range.
    #[test]
    fn two_samples_whose_read_groups_collide_are_renumbered_apart() {
        let written: Vec<SampleLocusObservations> = [0u32, 1]
            .iter()
            .enumerate()
            .map(|(at, group)| {
                let mut record = a_record(0, 1 + at as u64 * 100, 1);
                record.observations[0].read_group = ReadGroupId(*group);
                record
            })
            .collect();
        let (_first_dir, first) = a_psp_of(&written);
        let (_second_dir, second) = a_psp_of(&written);

        // What the calling stage's merged table hands each source: the first file's groups keep
        // their numbers, the second file's start where the first's ended.
        for (path, remap) in [
            (&first, vec![ReadGroupId(0), ReadGroupId(1)]),
            (&second, vec![ReadGroupId(2), ReadGroupId(3)]),
        ] {
            let mut psp = PspReader::open(path).expect("a finished psp opens");
            let mut source = PspObservationSource::over(&mut psp, &remap).expect("the walk starts");
            let (read, failed) = drain(&mut source);
            assert!(failed.is_none(), "the file is sound: {failed:?}");
            let groups: Vec<ReadGroupId> = read
                .iter()
                .map(|record| record.observations[0].read_group)
                .collect();
            assert_eq!(
                groups, remap,
                "each stored group reaches the merge under the number this run gave it",
            );
        }
    }

    /// **A number past the end of the file's own table is the file disagreeing with itself**,
    /// and there is nothing to renumber it into — so it is refused rather than left to panic in
    /// the middle of a cohort, where nothing would say which file to look at.
    #[test]
    fn an_observation_naming_a_read_group_the_file_does_not_declare_is_refused() {
        let mut record = a_record(0, 1, 1);
        record.observations[0].read_group = ReadGroupId(1);
        let written = vec![record];
        let (_dir, path) = a_psp_of(&written);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        // A table of one, where the record names group 1.
        let mut source =
            PspObservationSource::over(&mut psp, &[ReadGroupId(7)]).expect("the walk starts");

        let (read, failed) = drain(&mut source);
        assert!(read.is_empty());
        let Some(RunError::SourceFailed { source: cause, .. }) = failed else {
            panic!("an unknown read group must fail the source: {failed:?}");
        };
        let cause = cause
            .downcast::<PspSourceError>()
            .expect("the cause is this source's own");
        let PspSourceError::ReadGroupNotInThisFilesTable {
            at,
            names,
            in_the_table,
        } = *cause
        else {
            panic!("the read group is what was wrong: {cause:?}");
        };
        assert_eq!((at, names, in_the_table), (written[0].region, 1, 1));
    }

    /// **A source over a psp can cross a thread**, which is what the merge's parallel cover
    /// needs of it: `merge_cohort_handing_each_locus_over_covering_samples_in_parallel` sweeps
    /// the cohort's sources from a worker pool, one thread at a time each, and refuses a source
    /// that is not `Send`. The walker earned the same proof at direct mode's Milestone E
    /// (`a_run_walker_can_cross_a_thread`); this is psp mode's half, and it is a compile-time
    /// check rather than a runtime one — a `!Send` source would not build this test.
    #[test]
    fn a_psp_source_can_cross_a_thread() {
        fn only_takes_send<T: Send>(_: &T) {}
        let (_dir, path) = a_psp_of(&a_sample());
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let source = PspObservationSource::over(&mut psp, &as_walked()).expect("the walk starts");
        only_takes_send(&source);
    }

    /// **A head summarises its record exactly** — every field, on every record of a psp,
    /// including the ones where the two counts differ.
    ///
    /// This is the property the deferred build is bought with: if a head's summary can differ
    /// from the summary of the record behind it, then deciding from heads decides on something
    /// other than the evidence, and the cohort keeps or drops the wrong loci — silently, since
    /// nothing downstream ever sees the discarded records. The fixture deliberately includes a
    /// record whose reads are not all complete, because that is where numerator and denominator
    /// come apart: a read whose bases stop inside the locus is counted by neither, so a summary
    /// that took read depth for its denominator would pass this test only on the easy records.
    #[test]
    fn a_head_summarises_a_record_exactly() {
        let mut records = a_sample();
        // One record with a partial read beside its complete ones, and one whose only
        // non-reference evidence is partial — so `reads_compared_with_reference` is below the
        // reads that covered the ground, and `non_reference_reads` ignores what it cannot
        // compare.
        let mut mixed = a_record(0, 5_000, 4);
        let mut partial = mixed.observations[0].clone();
        partial.read_witness = ReadWitness::Partial {
            positions: WitnessedLocusPositions::one_run_from_offset_and_length(0, 2)
                .expect("a two-base witness inside a four-base locus"),
        };
        partial.bases = Box::from(&b"AC"[..]);
        partial.num_obs = 7;
        mixed.observations.push(partial);
        let mut varying = a_record(0, 6_000, 3);
        varying.observations[0].bases = Box::from(&b"TTT"[..]);
        records.push(mixed);
        records.push(varying);
        records.sort_by_key(|record| (record.region.contig.get(), record.region.start.get()));

        let (_dir, path) = a_psp_of(&records);
        let mut psp = PspReader::open(&path).expect("the file opens");

        let mut checked = 0usize;
        let mut differing_counts = 0usize;
        for found in psp.records().expect("the walk starts") {
            let found = found.expect("the fixture reads back");
            let record = found.record.as_ref().expect("a full walk builds every body");
            let from_the_head = LocusSummary::from(&found.head);
            assert_eq!(
                from_the_head,
                LocusSummary::of(record),
                "the head and the record disagree at {}",
                found.head.region,
            );
            if from_the_head.non_reference_reads != from_the_head.reads_compared_with_reference {
                differing_counts += 1;
            }
            checked += 1;
        }

        assert_eq!(checked, records.len(), "every record was compared");
        assert!(differing_counts > 0, "the fixture separates the two counts");

        // **The sharp case, asserted on the record built for it.** At 5,000 the fixture has
        // three complete reads and seven partial ones. A denominator taken from read depth
        // would say ten; the rule's denominator is the reads whose whole sequence over the
        // locus could be compared, which is three — and the head has to say three, or a run
        // deciding from heads applies a bar ten-thirds too high at exactly the sites where
        // reads run out mid-locus.
        let mut psp = PspReader::open(&path).expect("the file opens again");
        let mixed_head = psp
            .records()
            .expect("the walk starts")
            .map(|found| found.expect("the fixture reads back").head)
            .find(|head| head.region.start.get() == 5_000)
            .expect("the mixed record is in the file");
        assert_eq!(
            mixed_head.reads_compared_with_reference, 3,
            "the seven partial reads are in neither count"
        );
        assert_eq!(
            mixed_head.non_reference_reads, 0,
            "and the complete reads all matched the reference"
        );
    }


    /// **Reading a sample as summaries sees exactly what reading it as records sees**, and the
    /// kept bytes build the records that the building walk built.
    ///
    /// The two halves are the two claims the deferred build stands on. The first is that a run
    /// deciding from heads decides on the same numbers a run holding records would: same
    /// coordinates, same keep-rule numerator and denominator, same order, same count. The
    /// second is that nothing is lost by not building — the evidence is still there, in the
    /// arena, and comes back identical.
    ///
    /// Built here **in reverse**, because that is the order a cohort will not use and therefore
    /// the one that would hide a dependence between records: a body decoded against the wrong
    /// state is a plausible body, not a refusal.
    #[test]
    fn reading_as_summaries_sees_what_reading_as_records_sees() {
        let mut records = a_sample();
        let mut mixed = a_record(0, 5_000, 4);
        let mut partial = mixed.observations[0].clone();
        partial.read_witness = ReadWitness::Partial {
            positions: WitnessedLocusPositions::one_run_from_offset_and_length(0, 2)
                .expect("a two-base witness inside a four-base locus"),
        };
        partial.bases = Box::from(&b"AC"[..]);
        partial.num_obs = 7;
        mixed.observations.push(partial);
        records.push(mixed);
        records.sort_by_key(|record| (record.region.contig.get(), record.region.start.get()));

        let (_dir, path) = a_psp_of(&records);

        let built: Vec<_> = {
            let mut psp = PspReader::open(&path).expect("the file opens");
            psp.records()
                .expect("the walk starts")
                .map(|found| found.expect("the fixture reads back"))
                .map(|found| found.record.expect("a full walk builds every body"))
                .collect()
        };

        let mut psp = PspReader::open(&path).expect("the file opens again");
        let mut source = PspSummarySource::over(&mut psp, &as_walked()).expect("the summary walk starts");
        let mut kept = Vec::new();
        while let Some(next) = source.next_summary() {
            kept.push(next.expect("the fixture reads back"));
        }

        assert_eq!(kept.len(), built.len(), "the two walks met the same records");
        for (at, (summarised, record)) in kept.iter().zip(&built).enumerate() {
            assert_eq!(
                summarised.summary,
                LocusSummary::of(record),
                "record {at}: the summary read differs from the record's own"
            );
        }

        // And the evidence is still there. Reverse, for the reason the doc gives.
        let layout = {
            let reference = PspReader::open(&path).expect("the file opens once more");
            crate::ng::psp::record::RecordLayout::from_manifest(&reference.header().manifest)
                .expect("the file's own manifest")
        };
        let live = crate::ng::psp::chain_ids::LiveSet::default();
        for (at, kept_record) in kept.iter().enumerate().rev() {
            let head = crate::ng::psp::RecordHead {
                region: kept_record.summary.region,
                non_reference_reads: kept_record.summary.non_reference_reads,
                reads_compared_with_reference: kept_record.summary.reads_compared_with_reference,
                body_bytes: u32::try_from(kept_record.body.len()).expect("a small body"),
            };
            let found = crate::ng::psp::record::LocatedRecord {
                head,
                body: &source.kept()[kept_record.body.clone()],
                record_bytes: kept_record.body.len(),
            };
            let rebuilt = crate::ng::psp::record::decode_the_body_of(&found, &live, &layout)
                .expect("a kept body builds on its own")
                .record;
            assert_eq!(
                rebuilt, built[at],
                "record {at} built out of order differs from the one the building walk built"
            );
        }
    }


    /// **The summary source, driven as the cache drives it, yields the records the building
    /// source yields** — every draw kept, every body built afterwards and out of order.
    ///
    /// This is the deferred build end to end at the source's own level: the cohort's decisions
    /// would be made on the summaries this returns, and its evidence on the records these
    /// builds produce. If the two sources can disagree, psp mode calls something direct mode
    /// does not, and nothing downstream would say so.
    #[test]
    fn keeping_and_building_later_gives_what_building_at_once_gives() {
        let records = a_sample();
        let (_dir, path) = a_psp_of(&records);
        let groups = as_walked();

        let built: Vec<_> = {
            let mut psp = PspReader::open(&path).expect("the file opens");
            let mut source =
                PspObservationSource::over(&mut psp, &groups).expect("the walk starts");
            std::iter::from_fn(|| source.next_observation(None))
                .map(|next| next.expect("the fixture reads back"))
                .collect()
        };

        let mut psp = PspReader::open(&path).expect("the file opens again");
        let mut source = PspSummarySource::over(&mut psp, &groups).expect("the walk starts");
        let mut drawn = Vec::new();
        while let Some(next) = source.next_drawn(None) {
            match next.expect("the fixture reads back") {
                Drawn::Kept { summary, body } => drawn.push((summary, body)),
                Drawn::Built(_) => panic!("this source keeps every body"),
            }
        }

        assert_eq!(drawn.len(), built.len(), "the two sources met the same records");
        for (at, ((summary, _), record)) in drawn.iter().zip(&built).enumerate() {
            assert_eq!(
                *summary,
                LocusSummary::of(record),
                "record {at}: the kept summary differs from the built record's own"
            );
        }

        // Out of order, because that is the order a cohort will not use.
        for (at, (_, body)) in drawn.iter().enumerate().rev() {
            let rebuilt = source.build(body.clone()).expect("a kept body builds");
            assert_eq!(
                rebuilt, built[at],
                "record {at} built out of order differs from the one built at once"
            );
        }
    }

}
