//! **Pass three: the verdict on each parked record, and the file that comes out.**
//!
//! Pass one parked every record the run would have written; pass two scored them all and
//! resolved the operator's target false-discovery rate to a cut. This reads the spill a second
//! time, in the same order, and writes the VCF — dropping the records the cut removes, or, when
//! the operator asked to see them, writing them on the `hiddenDup` filter instead.
//!
//! **The *i*th record read here is the *i*th ratio pass two kept**, and that pairing is the whole
//! of how a verdict reaches a record: the ratios carry no contig and no position. So the one
//! index a ratio is looked up at is the count of entries the reading thread copied before the
//! record — never a number derived from the record — and the batches come back in the order they
//! were read. A read that filtered, chunked or sharded would give every record its neighbour's
//! verdict, with nothing about the file looking wrong.
//!
//! **A record that is written comes out byte for byte as pass one encoded it, except for the two
//! columns the verdict touches.** That is what makes spec §10's oracle checkable rather than
//! hoped for: a run at a target no record reaches, with the two `INFO` keys stripped, is the
//! filter-off run exactly. The splice is [`rewrite_filter_and_info`]'s, and it is the same call
//! whether anything is added or not.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.5.

use crate::vcf::{
    HIDDEN_PARALOG_FILTER_ID, PARALOG_POSTERIOR_DECIMALS, PARALOG_RATIO_DECIMALS, RecordPlace,
    VcfWriteError, VcfWriter,
};

use std::io::BufRead;
use std::sync::mpsc::{self, Receiver, Sender, SyncSender};

use rayon::prelude::*;

use super::{
    CompressedSpillEntry, LinePatchError, ParalogVerdicts, SpillEntryDecoder, SpillError,
    SpillFile, SpillFileError, SpillReaderWithFragments, rewrite_filter_and_info,
};
use crate::read_phase::fragments::RecordFragments;
use crate::run::explain::{ParalogExplanations, ParalogVerdictAt};
use crate::types::ContigId;

/// **What the filter did to the run's records**, for the run report and spec §3.5's lines.
///
/// The four counts partition the records pass one parked: every parked record is written or
/// dropped, and every written one is tagged or untouched.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct WhatTheFilterDid {
    /// Records written to the VCF — the survivors, plus the tagged ones.
    pub written: u64,
    /// Records the cut removed and the file does not contain. Zero when the operator asked for
    /// tagging.
    pub dropped: u64,
    /// Records the cut removed but the operator asked to keep, written on the filter's id. Zero
    /// when the operator did not.
    pub tagged: u64,
    /// Records no sample could speak for. They are written untouched and carry neither `INFO`
    /// field — **and they are counted**, because a run where this is every record produced a
    /// file the filter never looked at, and nothing else about the output would say so.
    pub unscored: u64,
}

/// What can go wrong writing the calls back.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PassThreeError {
    /// The parked records could not be read back.
    #[error("the records parked for the hidden-duplication filter could not be read back")]
    Spill(#[source] SpillFileError),
    /// A parked line could not be given its verdict. **Names the record**, because the spill
    /// holds millions and a failure that named only the reason would not say which.
    #[error(
        "the verdict could not be put on the record at contig {contig} position {position}, \
         which is record {ordinal} of the spill"
    )]
    Patch {
        /// The record's contig.
        contig: u32,
        /// The record's written position.
        position: u64,
        /// Which record of the spill it is, counting from one.
        ordinal: u64,
        /// What was wrong with the line.
        #[source]
        source: LinePatchError,
    },
    /// A record could not be written.
    #[error("the record at contig {contig} position {position} could not be written")]
    Write {
        /// The record's contig.
        contig: u32,
        /// The record's written position.
        position: u64,
        /// What the writer said.
        #[source]
        source: VcfWriteError,
    },
    /// The spill holds a different number of records than pass two scored.
    ///
    /// **Not reachable from a single run**, since both passes read one file whose record count
    /// the reader checks — but the two are paired by position and nothing in the types says
    /// they came from the same run, so the pairing is checked rather than assumed.
    #[error(
        "the spill holds {records} record(s) and the scoring pass produced {ratios} ratio(s); \
         pairing them by position would give records their neighbours' verdicts"
    )]
    RatiosDoNotMatchTheSpill {
        /// How many records the spill holds.
        records: u64,
        /// How many ratios pass two produced.
        ratios: usize,
    },
}

#[cfg(test)]
/// **Write the VCF, applying each record's verdict.**
///
/// `tag_instead_of_dropping` is the operator's `--hidden-dup-tag`: with it, a record the cut
/// removes is written on the `hiddenDup` filter rather than left out, so a run that drops can
/// be audited against one that does not.
///
/// The writer is not finished here — the caller owns that, because a run has more to do after the
/// last record.
///
/// # Errors
///
/// If the spill cannot be read, if a line cannot be patched, if a record cannot be written, or if
/// the spill and the ratios do not describe the same run.
pub fn write_the_records_the_filter_kept(
    spill: &SpillFile,
    verdicts: &ParalogVerdicts,
    tag_instead_of_dropping: bool,
    writer: &mut VcfWriter,
) -> Result<WhatTheFilterDid, PassThreeError> {
    write_the_records_the_filter_kept_explaining(
        spill,
        verdicts,
        tag_instead_of_dropping,
        writer,
        None,
    )
}

/// `write_the_records_the_filter_kept`, also keeping the verdict of every record `explained`
/// asks about (`--explain-loci`, `doc/devel/ng/spec/explain_loci.md`). **What is written does
/// not depend on it**: the verdict is read where it is decided and nothing else is touched.
///
/// # Errors
///
/// As `write_the_records_the_filter_kept`.
pub fn write_the_records_the_filter_kept_explaining(
    spill: &SpillFile,
    verdicts: &ParalogVerdicts,
    tag_instead_of_dropping: bool,
    writer: &mut VcfWriter,
    mut explained: Option<&mut ParalogExplanations>,
) -> Result<WhatTheFilterDid, PassThreeError> {
    if spill.entries_written() != verdicts.ratios.len() as u64 {
        return Err(PassThreeError::RatiosDoNotMatchTheSpill {
            records: spill.entries_written(),
            ratios: verdicts.ratios.len(),
        });
    }

    if let Some(windowed) = &verdicts.windowed
        && windowed.posteriors.len() != verdicts.ratios.len()
    {
        return Err(PassThreeError::RatiosDoNotMatchTheSpill {
            records: spill.entries_written(),
            ratios: windowed.posteriors.len(),
        });
    }

    // The one reader of the fragments: they go to the writer beside each line.
    let entries = spill.read_with_fragments().map_err(PassThreeError::Spill)?;
    let mut did = WhatTheFilterDid::default();

    // **Three threads, each doing what only it can.** The spill has no index, so one thread walks
    // it — and only copies each entry's two compressed frames off the stream. A second hands each
    // batch to the pool, which unpacks the record and its fragments, reaches its verdict and puts
    // it on the line. This thread keeps what depends on the order: the counts, the explained loci
    // and the writer, whose own stage runs on the same pool. Before, this thread did all of it,
    // and the pool waited on the unpacking while the reading waited on the writer.
    //
    // **The *i*th record still meets the *i*th ratio.** The reading thread numbers every entry as
    // it copies it, and that number — the count of entries read before it — is the one index
    // the verdict is looked up at. Nothing filters, reorders or skips between the read and the
    // lookup, and the batches come back here in the order they were read, so the pairing is the
    // serial walk's.
    std::thread::scope(|scope| -> Result<(), PassThreeError> {
        let (read_batches, batches_to_prepare) = mpsc::sync_channel::<ReadBatch>(1);
        let (emptied_batches, batches_to_refill) = mpsc::channel::<Vec<CompressedSpillEntry>>();
        let (prepared_batches, batches_to_write) = mpsc::sync_channel::<PreparedBatch>(1);
        scope.spawn(move || read_ahead(entries, read_batches, batches_to_refill));
        scope.spawn(move || {
            prepare_ahead(
                spill,
                verdicts,
                tag_instead_of_dropping,
                batches_to_prepare,
                emptied_batches,
                prepared_batches,
            );
        });

        for batch in batches_to_write {
            for record in batch.records {
                let record = record?;
                write_one(
                    record,
                    tag_instead_of_dropping,
                    writer,
                    explained.as_deref_mut(),
                    &mut did,
                )?;
            }
            if let Some(source) = batch.ended_with {
                return Err(PassThreeError::Spill(spill.naming_this_file(source)));
            }
        }
        Ok(())
    })?;

    Ok(did)
}

/// How many entries the pool prepares at once. **Smaller than pass two's batch** because a prepared
/// record holds its whole line and its fragments — about 450 kB at 4,352 samples — and up to three
/// prepared batches are in flight, so about 350 MB there.
const PREPARING_BATCH: usize = 256;

/// **A batch of compressed entries, as the reading thread hands it over.** The first `filled` are
/// this batch's, the first of them the spill's entry number `first`, counting from 0; the rest are
/// buffers kept from an earlier, fuller batch. `ended_with` is the reader's failure, when the batch
/// stopped short because of one.
struct ReadBatch {
    entries: Vec<CompressedSpillEntry>,
    filled: usize,
    first: usize,
    ended_with: Option<SpillError>,
}

/// **A batch the pool has prepared**, one item a record in spill order, then the reader's
/// failure, if the batch stopped short because of one.
struct PreparedBatch {
    records: Vec<Result<PreparedRecord, PassThreeError>>,
    ended_with: Option<SpillError>,
}

/// **One record with its verdict reached and, unless it is dropped, its line patched** — all a
/// pure function of the record and its ratio, so it is done on the pool.
struct PreparedRecord {
    place: RecordPlace,
    contig: ContigId,
    written_at: u64,
    ratio: f64,
    posterior: Option<f64>,
    flagged: bool,
    /// The patched line; `None` for a record the cut drops, which is never written.
    line: Option<Vec<u8>>,
    fragments: RecordFragments,
}

/// **The reading thread**: fill batches of compressed entries from `reader`, number them, and hand
/// each on, until the spill ends, the reader fails, or nothing is listening.
///
/// It refills the batches handed back rather than allocating new ones, and never waits for one: a
/// batch not yet returned is replaced by a fresh one.
fn read_ahead<R: BufRead>(
    mut reader: SpillReaderWithFragments<R>,
    read_batches: SyncSender<ReadBatch>,
    batches_to_refill: Receiver<Vec<CompressedSpillEntry>>,
) {
    let mut first = 0;
    loop {
        let mut entries = batches_to_refill.try_recv().unwrap_or_default();
        if entries.len() < PREPARING_BATCH {
            entries.resize_with(PREPARING_BATCH, CompressedSpillEntry::default);
        }
        let mut filled = 0;
        let mut ended_with = None;
        let mut at_the_end = false;
        while filled < PREPARING_BATCH {
            match reader.next_compressed_entry(&mut entries[filled]) {
                Some(Ok(())) => filled += 1,
                Some(Err(error)) => {
                    ended_with = Some(error);
                    at_the_end = true;
                    break;
                }
                None => {
                    at_the_end = true;
                    break;
                }
            }
        }
        if filled == 0 && ended_with.is_none() {
            return;
        }
        let handed_over = read_batches.send(ReadBatch {
            entries,
            filled,
            first,
            ended_with,
        });
        first += filled;
        if handed_over.is_err() || at_the_end {
            return;
        }
    }
}

/// **The preparing thread**: unpack each read batch on the pool, reach every record's verdict and
/// patch its line, hand the prepared batch on in the order it came, and return the emptied buffers
/// to the reading thread.
fn prepare_ahead(
    spill: &SpillFile,
    verdicts: &ParalogVerdicts,
    tag_instead_of_dropping: bool,
    batches_to_prepare: Receiver<ReadBatch>,
    emptied_batches: Sender<Vec<CompressedSpillEntry>>,
    prepared_batches: SyncSender<PreparedBatch>,
) {
    for batch in batches_to_prepare {
        let ReadBatch {
            entries,
            filled,
            first,
            ended_with,
        } = batch;
        let mut records = Vec::with_capacity(filled);
        entries[..filled]
            .par_iter()
            .enumerate()
            .map_init(SpillEntryDecoder::new, |decoder, (offset, entry)| {
                prepare_one(
                    spill,
                    decoder,
                    entry,
                    first + offset,
                    verdicts,
                    tag_instead_of_dropping,
                )
            })
            .collect_into_vec(&mut records);
        // The reading thread may already have finished; a batch it no longer wants is dropped.
        let _ = emptied_batches.send(entries);
        if prepared_batches
            .send(PreparedBatch {
                records,
                ended_with,
            })
            .is_err()
        {
            return;
        }
    }
}

/// **The verdict on the spill's entry number `index`, and its line patched.**
///
/// `index` counts from 0 and is the position of the record's ratio in `verdicts`: the count of
/// entries read before it, which is the whole of the pairing (see the module's notes).
fn prepare_one(
    spill: &SpillFile,
    decoder: &mut SpillEntryDecoder,
    entry: &CompressedSpillEntry,
    index: usize,
    verdicts: &ParalogVerdicts,
    tag_instead_of_dropping: bool,
) -> Result<PreparedRecord, PassThreeError> {
    let (entry, fragments) = decoder
        .decode_with_fragments(entry)
        .map_err(|source| PassThreeError::Spill(spill.naming_this_file(source)))?;
    // PANIC-FREE: the lengths were compared with the spill's count before the walk, and the
    // reader refuses a spill holding any other number of entries, so `index` is below both.
    let ratio = verdicts.ratios[index];

    // **The verdict and the probability**: the window's where the run used one, or the run-wide
    // rate's at a window of 0. Either way a record with no finite ratio has no probability and is
    // not flagged.
    let (flagged, posterior) = match &verdicts.windowed {
        Some(windowed) => {
            let posterior = windowed.posteriors[index];
            (
                ratio.is_finite() && windowed.flags(posterior),
                (ratio.is_finite() && posterior.is_finite()).then_some(posterior),
            )
        }
        None => (
            verdicts.calibration.flags(ratio),
            if ratio.is_finite() {
                verdicts.calibration.posterior(ratio)
            } else {
                None
            },
        ),
    };

    let line = if flagged && !tag_instead_of_dropping {
        None
    } else {
        // **What decides the fields is whether the record was scored — a finite ratio** (spec
        // §3.5). An unscored record gets neither, because it has nothing to report.
        let mut info_to_add: Vec<String> = Vec::new();
        if ratio.is_finite() {
            info_to_add.push(format!("PARALOG_LR={:.*}", PARALOG_RATIO_DECIMALS, ratio));
            // **The probability can be missing where the ratio is not.** It needs the run's
            // fitted duplication rate to be a rate — strictly between none and all — and a run
            // that fitted 0 or 1 has no log-odds to add. The copied calibration answers `None`
            // there rather than a saturated 0 or 1, and the field is then simply absent.
            if let Some(posterior) = posterior {
                info_to_add.push(format!(
                    "PARALOG_POST={:.*}",
                    PARALOG_POSTERIOR_DECIMALS, posterior
                ));
            }
        }
        let filter_to_add = flagged.then_some(HIDDEN_PARALOG_FILTER_ID);
        Some(
            rewrite_filter_and_info(&entry.line, filter_to_add, &info_to_add).map_err(
                |source| PassThreeError::Patch {
                    contig: entry.contig.get(),
                    position: entry.position.get(),
                    // Counting from one, because it is a place in a file a person will look for.
                    ordinal: index as u64 + 1,
                    source,
                },
            )?,
        )
    };

    Ok(PreparedRecord {
        place: RecordPlace::from(&entry),
        contig: entry.contig,
        written_at: entry.position.get(),
        ratio,
        posterior,
        flagged,
        line,
        fragments,
    })
}

/// **What only the writing thread can do with a prepared record**, in spill order: keep its
/// verdict where `--explain-loci` asked, count it, and write it unless the cut dropped it.
fn write_one(
    record: PreparedRecord,
    tag_instead_of_dropping: bool,
    writer: &mut VcfWriter,
    explained: Option<&mut ParalogExplanations>,
    did: &mut WhatTheFilterDid,
) -> Result<(), PassThreeError> {
    let PreparedRecord {
        place,
        contig,
        written_at,
        ratio,
        posterior,
        flagged,
        line,
        fragments,
    } = record;
    if let Some(explained) = explained
        && explained.wants(contig, written_at)
    {
        explained.record(ParalogVerdictAt {
            contig,
            written_at,
            ratio,
            posterior,
            flagged,
            dropped: flagged && !tag_instead_of_dropping,
        });
    }
    let Some(line) = line else {
        did.dropped += 1;
        return Ok(());
    };
    if !ratio.is_finite() {
        did.unscored += 1;
    }
    if flagged {
        did.tagged += 1;
    }
    writer
        .write_line(place, line, fragments)
        .map_err(|source| PassThreeError::Write {
            contig: contig.get(),
            position: written_at,
            source,
        })?;
    did.written += 1;
    Ok(())
}

#[cfg(test)]
mod tests;
