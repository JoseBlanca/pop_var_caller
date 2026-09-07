//! **Pass three: the verdict on each parked record, and the file that comes out.**
//!
//! Pass one parked every record the run would have written; pass two scored them all and
//! resolved the operator's target false-discovery rate to a cut. This reads the spill a second
//! time, in the same order, and writes the VCF — dropping the records the cut removes, or, when
//! the operator asked to see them, writing them on the `hiddenParalog` filter instead.
//!
//! **The *i*th record read here is the *i*th ratio pass two kept**, and that pairing is the whole
//! of how a verdict reaches a record: the ratios carry no contig and no position. So this walks
//! the spill front to back beside the ratio vector and never indexes into it — a read that
//! filtered, chunked or sharded would give every record its neighbour's verdict, with nothing
//! about the file looking wrong.
//!
//! **A record that is written comes out byte for byte as pass one encoded it, except for the two
//! columns the verdict touches.** That is what makes spec §10's oracle checkable rather than
//! hoped for: a run at a target no record reaches, with the two `INFO` keys stripped, is the
//! filter-off run exactly. The splice is [`rewrite_filter_and_info`]'s, and it is the same call
//! whether anything is added or not.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.5.

use crate::ng::vcf::{RecordPlace, VcfWriteError, VcfWriter};

use super::{LinePatchError, ParalogVerdicts, SpillFile, SpillFileError, rewrite_filter_and_info};

/// **The filter's id, as it appears in a written record's `FILTER` column.**
///
/// Declared in the header by [`crate::ng::vcf::header`] whenever the filter ran.
pub const HIDDEN_PARALOG_FILTER_ID: &str = "hiddenParalog";

/// **How a record's likelihood ratio is written**, and it is the same number of decimals the
/// header's `lr_cut` carries — so an operator auditing why a record went can compare the two
/// as written, without wondering whether a difference is real or a rounding.
const RATIO_DECIMALS: usize = 4;

/// How a record's posterior probability is written — the precision the header gives `pi`, which
/// is the other number it is read against.
const POSTERIOR_DECIMALS: usize = 6;

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

/// **Write the VCF, applying each record's verdict.**
///
/// `tag_instead_of_dropping` is the operator's `--paralog-filter-tag`: with it, a record the cut
/// removes is written on the `hiddenParalog` filter rather than left out, so a run that drops can
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
    if spill.entries_written() != verdicts.ratios.len() as u64 {
        return Err(PassThreeError::RatiosDoNotMatchTheSpill {
            records: spill.entries_written(),
            ratios: verdicts.ratios.len(),
        });
    }

    let entries = spill.read().map_err(PassThreeError::Spill)?;
    let mut did = WhatTheFilterDid::default();
    // **Taken front to back, never indexed.** The pairing is positional and this is what makes
    // it hold: the iterator advances once per record, in step with the read.
    let mut ratios = verdicts.ratios.iter();

    for (read_so_far, entry) in entries.enumerate() {
        let entry =
            entry.map_err(|source| PassThreeError::Spill(spill.naming_this_file(source)))?;
        // Counting from one, because it is a place in a file a person will look for.
        let ordinal = read_so_far as u64 + 1;
        // PANIC-FREE: the counts were compared above and both are walked once per record.
        let ratio = *ratios
            .next()
            .expect("the spill and the ratios were counted equal before the walk");

        let flagged = verdicts.calibration.flags(ratio);
        if flagged && !tag_instead_of_dropping {
            did.dropped += 1;
            continue;
        }

        // **What decides the fields is whether the record was scored — a finite ratio** (spec
        // §3.5). An unscored record gets neither, because it has nothing to report.
        let mut info_to_add: Vec<String> = Vec::new();
        if ratio.is_finite() {
            info_to_add.push(format!("PARALOG_LR={:.*}", RATIO_DECIMALS, ratio));
            // **The probability can be missing where the ratio is not.** It needs the run's
            // fitted duplication rate to be a rate — strictly between none and all — and a run
            // that fitted 0 or 1 has no log-odds to add. The copied calibration answers `None`
            // there rather than a saturated 0 or 1, and the field is then simply absent.
            if let Some(posterior) = verdicts.calibration.posterior(ratio) {
                info_to_add.push(format!("PARALOG_POST={:.*}", POSTERIOR_DECIMALS, posterior));
            }
        } else {
            did.unscored += 1;
        }

        let filter_to_add = if flagged {
            did.tagged += 1;
            Some(HIDDEN_PARALOG_FILTER_ID)
        } else {
            None
        };

        let line = rewrite_filter_and_info(&entry.line, filter_to_add, &info_to_add).map_err(
            |source| PassThreeError::Patch {
                contig: entry.contig.get(),
                position: entry.position.get(),
                ordinal,
                source,
            },
        )?;

        writer
            .write_line(RecordPlace::from(&entry), &line)
            .map_err(|source| PassThreeError::Write {
                contig: entry.contig.get(),
                position: entry.position.get(),
                source,
            })?;
        did.written += 1;
    }

    Ok(did)
}

#[cfg(test)]
mod tests;
