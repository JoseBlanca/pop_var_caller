//! **Encoding a sample's psp records on a thread of its own, while the walk carries on.**
//!
//! # What is being moved, and why nothing downstream waits for it
//!
//! The walk yields one [`SampleLocusObservations`] per stored locus, hands it to the census,
//! and then hands it to the psp writer, which encodes it into the open block's payload. The
//! census reads the locus first and keeps its own numbers, so **the writer is the last thing
//! that wants a locus** and it can take ownership of it.
//!
//! Nothing the walk does next depends on what the writer produced: the next locus is generated
//! from reads and reference, not from anything the file holds. The only thing tying the writer
//! to the walk is the **order** records reach it in, and one writing thread with a FIFO queue
//! keeps that by construction. The bytes are the same bytes — [`PspWriter`] is the same object
//! doing the same work in the same order, so which thread ran it cannot be read out of the file.
//!
//! **Measured, on 10 Mb of tomato chromosome 1 at 103.5×**: the record encoder was 14.8% of the
//! walking thread — `write_a_record` 10.3% and `encode_record_body_reusing` 4.5% — about 4.2
//! seconds of a 28.4-second run, and the walking thread was busy 99.8% of the time, so it was
//! the run's own critical path.
//!
//! # Failures arrive a few records late, and one of them is not this module's to change
//!
//! `push` is no longer called on the walking thread, so a record the store refuses raises its
//! error at the walk's *next* handover or at [`finish`](PspWriterLine::finish) rather than at
//! the record that caused it. Nothing a reader can see changes — the file stops in the same
//! place and [`PspWriter::finish`] still refuses to seal it — and the error still names the
//! locus, because the thread holds the locus it failed on. This is the same shift the psp's
//! block compression already made when it moved to a thread.
//!
//! **After a failure the thread keeps draining and stops writing.** A thread that returned
//! early would drop its receiver, and the walk would then fail on its next `send` with nothing
//! but "the channel closed" to say — losing the real error. So the first failure is kept, every
//! later locus is read and dropped, and the kept failure is what [`finish`](PspWriterLine::finish)
//! returns.
//!
//! # ⚠ Closing the queue is not the same as asking for the file
//!
//! **A walk that fails drops this line without calling [`finish`](PspWriterLine::finish)**, and
//! the file it leaves behind must not read back as whole — `write_psp_propagates_a_walk_failure`
//! is what says so, and it is the test this module broke on its first draft. A thread that
//! sealed whenever its queue closed would seal exactly then, and an interrupted walk would
//! produce a psp with a footer over a chromosome missing its tail.
//!
//! So the queue carries [`ToWriter`], and the seal is a **message** rather than the absence of
//! one: `finish` sends [`ToWriter::Seal`] and then closes, while `drop` only closes. A thread
//! that reaches the end of its queue without having been asked returns the writer unsealed and
//! drops it, which is what the walk did itself before the writer moved off its thread.
//!
//! # The seal carries the trailer, because the trailer is known only when the walk is over
//!
//! [`PspWriter::finish`] takes the closing payload the file will carry, and for a sample walk
//! that payload is its census — accumulated from every locus as it goes past and encoded once
//! the last one has (`psp_census_pair.md` §3.1). So the bytes have to reach the writing thread,
//! and the seal is the one message that can carry them: they cross the seam **once**, moved
//! rather than copied, at the end of a walk that has already handed over millions of loci.
//! (It is encoded once: plan step A2 removed the file that used to sit beside the psp, and the
//! trailer is the only copy.)
//!
//! A walk that builds no census closes with [`PspTrailer::Nothing`], which is what every psp in
//! this tree carried before. **That is a named choice rather than an empty `Vec`** — see the
//! type — because a psp sealed with no census reads back as a whole file, and the parameters fit
//! is where the mistake would be found, a stage and a re-walk later.
//!
//! # Loci travel in batches, and that is not a tuning knob
//!
//! **Handing one locus at a time cost more than it saved.** A walk emits about 350,000 loci a
//! second at this fixture's depth, and a bounded channel operation is an atomic pair plus, when
//! the queue runs empty or full, a park and an unpark. Measured over 10 Mb: one locus per
//! message won **0.65 s of wall** and cost **13 s of user time** — 39.1 s to 51.9 s, with cycles
//! up 35.6% — which on the saturated machine a cohort actually runs on is a loss rather than a
//! win. Batching them at [`LOCI_PER_BATCH`] divides that traffic by the batch size.
//!
//! # Why the queue is bounded, and at what
//!
//! Unbounded, a walk that outruns its writer would hold every locus it has ever emitted — 9.7
//! million of them over 10 Mb, about a kilobyte each. Bounded, the walk waits, which is what it
//! did before this existed.
//!
//! [`BATCHES_AWAITING_ENCODING`] batches of [`LOCI_PER_BATCH`] is about 4 MB at this fixture's
//! depth. The writer is not the slow side — it is 4.2 s of work against the walk's remaining 24
//! — so what the queue is sized for is the **burst**: when a block closes, the writer hands it
//! to the compressing thread, and if that thread already holds its two blocks the writer waits
//! for one to come back. A block is 100 kb of genome by default, so those bursts are 100,000
//! loci apart, and the depth only has to cover one of them.

use std::path::PathBuf;

use crate::locus_generation::SampleLocusObservations;
use crate::psp::{PspWriter, WriteStats};
use crate::run::RunError;

/// **How many loci travel in one message.** See the module note: one at a time was measured and
/// cost 13 seconds of user time for 0.65 of wall.
const LOCI_PER_BATCH: usize = 256;

/// **How many batches may be waiting before the walk has to wait.** See the module note.
const BATCHES_AWAITING_ENCODING: usize = 16;

/// What the walk hands the writing thread: loci to write, or the request to seal the file.
///
/// **The seal is a message and not the end of the queue** — see the module note's ⚠.
enum ToWriter {
    Records(Vec<SampleLocusObservations>),
    /// Seal the file, closing it with this payload.
    Seal(PspTrailer),
}

/// **What a psp closes with**: the sample's census, or a deliberate nothing.
///
/// **A named type rather than a bare `Vec<u8>`, and the reason is what the mistake costs.** An
/// empty payload is a legal trailer and a psp that carries one is a whole psp — it opens, it
/// indexes, it reads back every record and reports the right count. So a caller that *meant* to
/// pass a census and passed `Vec::new()` produces a file nothing downstream refuses, and the
/// omission surfaces at the parameters fit, which is a stage and a whole re-walk away from the
/// line that caused it (`psp_census_pair.md` §3.3 names that state as the one the design exists
/// to remove). With this type, closing with nothing is a sentence the author had to write.
#[derive(Debug)]
pub enum PspTrailer {
    /// The sample's census, encoded (`psp_census_pair.md` §3.1).
    Census(Vec<u8>),
    /// No census was built, so the file closes with nothing.
    Nothing,
}

impl PspTrailer {
    /// What [`PspWriter::finish`] is handed.
    fn bytes(&self) -> &[u8] {
        match self {
            Self::Census(bytes) => bytes,
            Self::Nothing => &[],
        }
    }
}

/// The psp writer, running on a thread of its own.
///
/// Built by [`start`](Self::start), fed by [`push`](Self::push), and ended by
/// [`finish`](Self::finish) — which is the only way to get the file's [`WriteStats`] back, and
/// so the only way to learn that the file was written at all.
#[derive(Debug)]
pub struct PspWriterLine {
    /// `None` once [`finish`](Self::finish) has taken it, which is what closes the channel and
    /// lets the thread seal the file.
    to_write: Option<crossbeam_channel::Sender<ToWriter>>,
    worker: Option<std::thread::JoinHandle<Result<Option<WriteStats>, RunError>>>,
    /// The loci gathered since the last handover — see the module note on batching.
    batch: Vec<SampleLocusObservations>,
    /// Batch buffers the writing thread has emptied and given back, ready to be filled again,
    /// so a run of ten million loci allocates a few tens of these rather than one per batch.
    spare: crossbeam_channel::Receiver<Vec<SampleLocusObservations>>,
}

impl PspWriterLine {
    /// Start the writing thread, which owns `writer` from here on.
    ///
    /// The writer is **created by the caller** rather than here, so that a file that cannot be
    /// created fails where the path is, before any locus exists to be lost.
    pub fn start(writer: PspWriter, path: PathBuf) -> Self {
        let (to_write, queue) = crossbeam_channel::bounded(BATCHES_AWAITING_ENCODING);
        // **Unbounded, for the reason the block compressor's return channel is**: the walk may
        // block handing a batch over, and the writer must never block handing an empty one back,
        // or the two would wait on each other. It is bounded in fact by the sending side, which
        // never has more than `BATCHES_AWAITING_ENCODING` batches in flight.
        let (emptied, spare) = crossbeam_channel::unbounded();
        let worker = std::thread::Builder::new()
            .name("psp-record-encoding".to_string())
            .spawn(move || write_loci(writer, path, &queue, &emptied))
            .expect("a thread for encoding psp records");
        Self {
            to_write: Some(to_write),
            worker: Some(worker),
            batch: Vec::with_capacity(LOCI_PER_BATCH),
            spare,
        }
    }

    /// Hand one locus over.
    ///
    /// **A send that fails means the thread has gone**, which happens only if it panicked —
    /// a failure to write is kept and drained rather than returned early. The locus is dropped
    /// and the error comes out of [`finish`](Self::finish), which joins the thread and so
    /// re-raises the panic where the caller can see it.
    pub fn push(&mut self, locus: SampleLocusObservations) {
        self.batch.push(locus);
        if self.batch.len() >= LOCI_PER_BATCH {
            self.hand_the_batch_over();
        }
    }

    /// Send whatever is gathered, and take an emptied buffer to gather the next lot into.
    fn hand_the_batch_over(&mut self) {
        if self.batch.is_empty() {
            return;
        }
        let mut next = self.spare.try_recv().unwrap_or_default();
        next.clear();
        next.reserve(LOCI_PER_BATCH);
        let batch = std::mem::replace(&mut self.batch, next);
        if let Some(to_write) = self.to_write.as_ref() {
            let _ = to_write.send(ToWriter::Records(batch));
        }
    }

    /// Ask for the file, wait for the last record to be written, and seal it with `trailer`.
    ///
    /// **The seal is asked for and then the queue is closed**, in that order and both before the
    /// join: a join before the close would wait for a thread that is waiting for a locus, and a
    /// close before the request would be indistinguishable from the walk giving up (the module
    /// note's ⚠).
    ///
    /// **`trailer` is the file's closing payload**, handed straight to [`PspWriter::finish`] —
    /// the sample's census for a walk that built one, and [`PspTrailer::Nothing`] for one that
    /// did not. It is taken by value because it crosses the thread seam (the module note).
    pub fn finish(mut self, trailer: PspTrailer) -> Result<WriteStats, RunError> {
        self.hand_the_batch_over();
        if let Some(to_write) = self.to_write.as_ref() {
            // **A send that fails means the thread panicked** — nothing else drops the receiver
            // while this sender is alive — and the `join` below re-raises that panic. So the
            // census travelling in the message is not a loss that can pass for a success.
            let _ = to_write.send(ToWriter::Seal(trailer));
        }
        self.to_write = None;
        let worker = self
            .worker
            .take()
            .expect("the handle is taken only here, and `finish` consumes the line");
        match worker.join() {
            Ok(Ok(Some(stats))) => Ok(stats),
            Ok(Ok(None)) => unreachable!("`finish` asks for the seal before it closes the queue"),
            Ok(Err(failure)) => Err(failure),
            Err(panic) => std::panic::resume_unwind(panic),
        }
    }
}

impl Drop for PspWriterLine {
    /// **Joined on drop, so no writing thread outlives its psp.**
    ///
    /// Reached when the walk fails and `finish` is never called. Closing the channel first is
    /// what stops the join waiting forever; the thread's own result is discarded, because the
    /// failure the caller is already carrying is the one that matters. **No seal is asked for**,
    /// so the file is left without a footer and a reader refuses it as interrupted — see the
    /// module note's ⚠.
    fn drop(&mut self) {
        self.to_write = None;
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

/// The writing thread's whole body: write every locus in the order it arrives, and seal the
/// file if it is asked to.
///
/// `Ok(None)` is the queue ending without a [`ToWriter::Seal`] — the walk gave up, and the
/// writer is dropped unsealed. See the module note for that and for why a failure keeps
/// draining instead of returning.
fn write_loci(
    mut writer: PspWriter,
    path: PathBuf,
    queue: &crossbeam_channel::Receiver<ToWriter>,
    emptied: &crossbeam_channel::Sender<Vec<SampleLocusObservations>>,
) -> Result<Option<WriteStats>, RunError> {
    let mut failure: Option<RunError> = None;
    let mut sealing: Option<PspTrailer> = None;
    for message in queue {
        match message {
            ToWriter::Seal(trailer) => {
                // **One sender asks once.** `finish` is the only place a seal is sent, it
                // consumes the line, and the sender is never cloned — so a second seal would
                // mean a change elsewhere, and it would silently replace this payload with
                // another rather than fail.
                debug_assert!(
                    sealing.is_none(),
                    "the file was asked to seal twice, and the second census would win",
                );
                sealing = Some(trailer);
            }
            ToWriter::Records(mut batch) => {
                // **And nothing arrives after it**: `finish` hands the last batch over before
                // it asks for the seal. A locus arriving later would be written into the block
                // stream and sealed over, giving a whole-looking psp holding a record its own
                // census never saw.
                debug_assert!(
                    sealing.is_none(),
                    "a locus arrived after the seal was asked for",
                );
                for locus in batch.drain(..) {
                    if failure.is_some() {
                        continue;
                    }
                    if let Err(source) = writer.push(&locus) {
                        failure = Some(RunError::RecordNotWritten {
                            locus: locus.region,
                            source: Box::new(source),
                        });
                    }
                }
                // The buffer goes back to be filled again; a send that fails means the walk has
                // gone, and then nothing wants it.
                let _ = emptied.send(batch);
            }
        }
    }
    if let Some(failure) = failure {
        return Err(failure);
    }
    let Some(trailer) = sealing else {
        return Ok(None);
    };
    writer
        .finish(trailer.bytes())
        .map(Some)
        .map_err(|source| RunError::PspNotWritten {
            path,
            source: Box::new(source),
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::psp::PspWriteError;
    use crate::psp::writer::tests_support::{a_file, a_header, a_record};

    /// **A line dropped without being asked to seal leaves a file every reader refuses** — the
    /// invariant this module broke on its first draft, where the thread sealed whenever its
    /// queue closed and an interrupted walk produced a psp with a footer over a chromosome
    /// missing its tail.
    ///
    /// `write_psp_propagates_a_walk_failure` in `gatherer.rs` is the end-to-end statement of the
    /// same thing; this one names the mechanism, so a change here fails beside the change.
    #[test]
    fn a_line_dropped_without_being_asked_to_seal_leaves_no_readable_psp() {
        let (_dir, path) = a_file();
        {
            let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
            let mut line = PspWriterLine::start(writer, path.clone());
            for step in 0..10u64 {
                line.push(a_record(0, 1 + step * 1_000, 1));
            }
            // Dropped here, deliberately, without `finish` — which is what a walk that failed
            // part-way leaves behind.
        }
        let refused = PspWriter::append(&path).expect_err("nothing sealed it");
        assert!(
            matches!(refused, PspWriteError::Reopen { .. }),
            "an unsealed file is refused for append, not opened: {refused:?}",
        );
    }

    /// **Twenty lines abandoned in turn, each with loci in flight, neither hang nor pile up.**
    ///
    /// The path nothing else reaches: every other test either finishes its line or abandons one
    /// with nothing sent, so no batch is ever in the queue when the sender goes. Mirrors
    /// `writers_abandoned_with_blocks_in_flight_neither_hang_nor_pile_up`, which covers the same
    /// hazard one stage further down.
    ///
    /// **Mutation-verified twice**: a `drop` that joins the thread without first closing the
    /// queue hangs this, and ninety seconds was not enough for it to finish; a thread that seals
    /// whenever its queue closes fails it, along with
    /// `a_line_dropped_without_being_asked_to_seal_leaves_no_readable_psp` and, end to end,
    /// `write_psp_propagates_a_walk_failure`.
    #[test]
    fn many_abandoned_lines_with_loci_in_flight_neither_hang_nor_pile_up() {
        let (_dir, path) = a_file();
        for round in 0..20u64 {
            let writer = PspWriter::create(&path, a_header(1_000)).expect("a header each round");
            let mut line = PspWriterLine::start(writer, path.clone());
            // More than one batch, so the queue holds one while another is being filled.
            for step in 0..(LOCI_PER_BATCH as u64 * 3) {
                line.push(a_record(0, 1 + step + round, 1));
            }
            drop(line);
        }
        let refused = PspWriter::append(&path).expect_err("nothing sealed the last one");
        assert!(
            matches!(refused, PspWriteError::Reopen { .. }),
            "an unsealed file is refused for append, not opened: {refused:?}",
        );
    }

    /// **A record the store refuses still names its locus, though it is raised late.**
    ///
    /// The writer is no longer called on the walking thread, so the failure surfaces at
    /// [`PspWriterLine::finish`] rather than at the `push` that caused it. What must not change
    /// is *which* locus the message names, and that survives because the thread holds the locus
    /// it failed on.
    #[test]
    fn a_refused_record_is_reported_by_finish_and_names_its_locus() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let mut line = PspWriterLine::start(writer, path.clone());
        line.push(a_record(0, 1_000, 1));
        // Backwards, which the block builder refuses — and the walk is told nothing here.
        line.push(a_record(0, 900, 1));
        line.push(a_record(0, 1_100, 1));
        let failure = line
            .finish(PspTrailer::Nothing)
            .expect_err("the second record goes backwards");
        match failure {
            RunError::RecordNotWritten { locus, .. } => {
                assert_eq!(
                    locus.start.get(),
                    900,
                    "the failure names the record it was"
                );
            }
            other => panic!("expected the refused record to be named: {other:?}"),
        }
    }

    /// **A line asked to seal writes every locus it was given and reports the file.**
    #[test]
    fn a_finished_line_seals_a_readable_psp_holding_every_locus() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let mut line = PspWriterLine::start(writer, path.clone());
        // Two and a bit batches, so the last one is handed over by `finish` rather than by a
        // full batch — the path a walk that ends mid-batch takes, which is every walk.
        let loci = LOCI_PER_BATCH as u64 * 2 + 7;
        for step in 0..loci {
            line.push(a_record(0, 1 + step * 10, 1));
        }
        let stats = line.finish(PspTrailer::Nothing).expect("it seals");
        assert_eq!(
            stats.records, loci,
            "every locus handed over is in the file"
        );
        // The claim is that the call succeeds — a sealed file can be reopened for appending —
        // so the writer itself is not wanted. Bound and dropped rather than discarded with
        // `let _ =`, which would also swallow a writer this test had meant to use.
        let reopened = PspWriter::append(&path).expect("a sealed file reopens");
        drop(reopened);
    }

    /// **The bytes handed to `finish` are the bytes the sealed file carries.**
    ///
    /// This is the seam the census crosses (`psp_census_pair.md` §3.1): the payload is built on
    /// the walking thread and written by the writing one, and nothing else in the pipeline would
    /// notice if it arrived empty — a psp with an empty trailer is a well-formed psp, and every
    /// one this tree wrote before carried exactly that. So the round trip is asserted here,
    /// where the handover is, rather than only end to end in the gatherer.
    #[test]
    fn the_file_carries_the_trailer_the_line_was_asked_to_seal_with() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let mut line = PspWriterLine::start(writer, path.clone());
        line.push(a_record(0, 1, 1));
        let payload = b"a sample's census would be here".to_vec();
        let stats = line
            .finish(PspTrailer::Census(payload.clone()))
            .expect("it seals");
        assert_eq!(stats.records, 1, "the record is still in the file");

        let mut reader = crate::psp::PspReader::open(&path).expect("the psp opens");
        assert_eq!(
            reader.trailer().expect("the trailer reads"),
            payload,
            "the payload the line was given is the one the file holds",
        );
    }

    /// **A trailer the size a real census is.** The round trip above is 31 bytes, and until this
    /// change 31 bytes was the largest trailer anything in this tree had ever written — every
    /// other `finish` under test passes an empty slice or a short literal. A whole genome's
    /// census is a few megabytes of positions plus tens of megabytes of tracts
    /// (`psp_census_pair.md` §3.1), so **the size class this change moves the trailer to is the
    /// only one production will ever see, and it was the one class with no test**.
    ///
    /// What it would catch is silent: a truncation or an offset that only shows once the
    /// payload outgrows the writer's buffer leaves a psp that opens, indexes and reads back
    /// every record, and hands the fit a census that stops early — reported as a malformed
    /// census a whole pipeline stage from its cause. The bytes vary rather than repeat, because
    /// a constant payload would survive a defect that wrote one page twice.
    #[test]
    fn a_trailer_of_megabytes_round_trips_through_the_line() {
        let (_dir, path) = a_file();
        let writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        let mut line = PspWriterLine::start(writer, path.clone());
        line.push(a_record(0, 1, 1));
        let payload: Vec<u8> = (0..4_000_000u32).map(|at| (at % 251) as u8).collect();
        let stats = line
            .finish(PspTrailer::Census(payload.clone()))
            .expect("it seals");
        assert_eq!(stats.records, 1, "the record is still in the file");

        let mut reader = crate::psp::PspReader::open(&path).expect("the psp opens");
        let read_back = reader.trailer().expect("the trailer reads");
        assert_eq!(
            read_back.len(),
            payload.len(),
            "the trailer came back a different length than it went in",
        );
        assert!(
            read_back == payload,
            "the trailer first differs at byte {}",
            read_back
                .iter()
                .zip(&payload)
                .position(|(back, sent)| back != sent)
                .expect("the lengths match, so a difference has a position"),
        );
    }
}
