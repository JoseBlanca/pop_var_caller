//! **Pass one: where a finished record goes, and the flag that decides.**
//!
//! Every record the run would write reaches one sink. With the filter **off** that sink is the
//! VCF writer, exactly as it was before this work existed — and *exactly* is the point, because
//! the standing oracle for the whole filter is that `--paralog-fdr 0` reproduces the previous
//! run byte for byte (spec §1.1 goal 7, §10). Nothing is inserted in that path: no encode and
//! no copy — only the match that chooses the sink, taken once per record.
//!
//! With the filter **on**, the record's line is parked on the spill together with what the
//! scorer will need per sample, and **the VCF is not opened at all** — it cannot be, since no
//! record's verdict is known until every record has been scored, and a half-written VCF beside a
//! finished run is a file someone can mistake for the output (spec §2, §3.4).
//!
//! **The position the entry carries is the *written* one.** A left-padded deletion is written one
//! base before its span, and pass three re-runs the writer's ordering check from the entry's head
//! fields — so an entry carrying the span start would order the file against a position the file
//! does not contain. It comes from the same function the `POS` column does
//! ([`written_position`](crate::ng::vcf::encode)), which is the only way the two cannot drift.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.4, §3.6.

use crate::ng::types::Ploidy;
use crate::ng::vcf::encode::{record_line, written_position};
use crate::ng::vcf::{HeaderContig, VcfRecord, VcfWriteError, VcfWriter};
use crate::ng::window_coverage::WindowCoverage;

use super::{
    GenericLocusSample, RepeatTractSample, SpillEntry, SpillFile, SpillFileError, SpilledSamples,
};

/// **The spill entry for one finished record**, ready to append.
///
/// `windows` is dense over the run's samples, in the run's sample order, and so are the record's
/// sample columns.
///
/// # Panics
///
/// If the two disagree in length. They come from one place — `assemble_record` builds both from
/// the same evidence under its own assertion — so a disagreement is a wiring error rather than an
/// input, and the alternative is an entry that silently carries a prefix of the cohort.
///
/// **Which shape the rows take is the record's tract flag and nothing else** (spec §3.2). A
/// repeat tract's rows carry the window alone, because slippage makes its allele split unreadable
/// and the filter declines to use it; every other record's rows carry the reference count and
/// every alternative's reads summed.
#[must_use]
pub fn entry_for(
    record: &VcfRecord,
    windows: &[WindowCoverage],
    contigs: &[HeaderContig],
    ploidy: Ploidy,
) -> SpillEntry {
    // **Both slices are the run's samples, and a disagreement is a wiring error.** The zip below
    // would take the shorter of the two — which never *mis-pairs*, since both are prefixes of the
    // same sample order, but silently loses a suffix. Pass two does catch it, for either row
    // shape: its check is on the entry's row count and answers for a tract as well as a generic
    // locus. Asserting here names the sink instead of pointing at pass two after the whole
    // calling pass has finished.
    assert_eq!(
        windows.len(),
        record.sample_columns().len(),
        "the run's window coverage and the record's sample columns disagree in length, so the \
         entry would silently carry a prefix of the cohort"
    );

    // Read once and used twice, so the flag the entry carries and the rows it carries cannot
    // come from two different answers to the same question.
    let is_repeat_tract = record.is_repeat_tract();

    let samples = if is_repeat_tract {
        SpilledSamples::RepeatTract(
            windows
                .iter()
                .map(|window| RepeatTractSample { window: *window })
                .collect(),
        )
    } else {
        SpilledSamples::GenericLocus(
            windows
                .iter()
                .zip(record.sample_columns())
                .map(|(window, column)| {
                    let reads = column.read_counts.allele_reads();
                    GenericLocusSample {
                        window: *window,
                        ref_reads: reads.first().copied().unwrap_or(0),
                        // **Summed, not `reads[1]`.** A multiallelic site takes this path too,
                        // and the collapsed-duplication story asks what share of the reads is
                        // non-reference — not which alternative each one carried (spec §3.2).
                        alt_reads: reads.iter().skip(1).copied().sum(),
                    }
                })
                .collect(),
        )
    };

    SpillEntry {
        contig: record.region().contig,
        position: crate::ng::types::Position(written_position(record)),
        is_repeat_tract,
        line: record_line(record, contigs, ploidy).into_bytes(),
        samples,
    }
}

/// **Where pass one puts each finished record**, chosen once when the run starts.
///
/// The two variants are the two runs: the one that exists today, and the one that parks its
/// records for a verdict it cannot reach yet.
pub enum CalledRecordSink {
    /// **The filter is off.** Records go straight to the VCF, through the same call the run made
    /// before this work — which is what makes the byte-identity oracle a property of the code
    /// rather than a hope.
    StraightToTheVcf(VcfWriter),
    /// **The filter is on.** Records are parked; the VCF is not opened until pass three.
    ParkedOnTheSpill(SpillingSink),
}

/// The state pass one needs to park a record: the file, and what turns a record into a line.
pub struct SpillingSink {
    spill: SpillFile,
    contigs: Vec<HeaderContig>,
    ploidy: Ploidy,
}

impl SpillingSink {
    /// Park records beside `output`, encoding their lines against `contigs` at `ploidy`.
    #[must_use]
    pub fn beside(output: &std::path::Path, contigs: Vec<HeaderContig>, ploidy: Ploidy) -> Self {
        Self {
            spill: SpillFile::beside(output),
            contigs,
            ploidy,
        }
    }

    /// The file the records are parked in, for passes two and three.
    #[must_use]
    pub fn spill(&self) -> &SpillFile {
        &self.spill
    }

    /// Stop appending and make the bytes durable, so passes two and three can read them.
    ///
    /// # Errors
    ///
    /// If the file cannot be flushed.
    pub fn finish_parking(&mut self) -> Result<(), SpillFileError> {
        self.spill.finish_writing()
    }
}

impl CalledRecordSink {
    /// **Close whatever this sink opened.**
    ///
    /// Consuming, so a forgotten close is a missing output rather than a truncated one — and it
    /// is what lets the two call sites hand the sink over without recovering the writer from
    /// inside it. An earlier version pattern-matched the off variant back out and declared the
    /// other `unreachable!`, which was sound only because the caller had refused a non-zero
    /// target ninety lines earlier: a panic waiting for the step that moves that guard.
    ///
    /// # Errors
    ///
    /// If the VCF cannot be finished, or the parked records cannot be made durable.
    pub fn finish(self) -> Result<(), PassOneError> {
        match self {
            Self::StraightToTheVcf(writer) => writer.finish().map_err(PassOneError::Vcf),
            Self::ParkedOnTheSpill(mut sink) => sink.finish_parking().map_err(PassOneError::Spill),
        }
    }

    /// Take one finished record.
    ///
    /// # Errors
    ///
    /// If the VCF write fails or the record runs backwards (off), or the spill cannot be
    /// appended to (on).
    pub fn accept(
        &mut self,
        record: &VcfRecord,
        windows: &[WindowCoverage],
    ) -> Result<(), PassOneError> {
        match self {
            // **Nothing between the record and the writer.** Not a helper, not a re-encode —
            // this is the line the run ran before the filter existed.
            Self::StraightToTheVcf(writer) => {
                writer.write_record(record).map_err(PassOneError::Vcf)
            }
            Self::ParkedOnTheSpill(sink) => {
                let entry = entry_for(record, windows, &sink.contigs, sink.ploidy);
                sink.spill.append(&entry).map_err(PassOneError::Spill)
            }
        }
    }
}

/// What can go wrong taking a record in pass one.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PassOneError {
    /// The record could not be written to the VCF — the filter-off path.
    #[error("the record could not be written")]
    Vcf(#[source] VcfWriteError),
    /// The record could not be parked on the spill — the filter-on path.
    #[error("the record could not be parked for the hidden-duplication filter")]
    Spill(#[source] SpillFileError),
}

#[cfg(test)]
mod tests;
