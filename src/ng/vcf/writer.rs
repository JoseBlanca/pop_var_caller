//! **Writing the file** — the header, then records in genome order, to a file that appears
//! whole or not at all.
//!
//! Two jobs, and they are separate on purpose. **Ordering** is a property of the record stream
//! and is checked here, because a VCF whose records run backwards is not a VCF and no consumer
//! would notice until it indexed one. **Durability** is a property of the file and is the
//! sink's: bytes go to `<output>.tmp` and are renamed into place only once they are on disk, so
//! a crash leaves no half-written VCF for anyone to mistake for a finished run.
//!
//! Encoding, by contrast, happens in [`encode`](super::encode) and cannot fail — everything
//! that could make a record unwritable was refused when the record was built. So this module is
//! where the `Result`s live, and all of them are about the file rather than about the data.

use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use thiserror::Error;

use super::encode::record_line;
use super::header::{VcfHeaderMetadata, header_text};
use super::{PaddingBase, VcfRecord};
use crate::ng::types::{ContigId, Ploidy};

/// The 28-byte empty-block marker every well-formed bgzf file ends with.
///
/// `noodles_bgzf`'s writer emits it on finish; the constant exists so a test can assert the
/// bytes are there without reaching into the library, since it is what makes `tabix` and
/// `bcftools` accept the file at all.
#[cfg(test)]
pub(crate) const BGZF_EOF: &[u8; 28] = &[
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Buffer size for the plain-text sink, so per-record writes coalesce before reaching the
/// kernel. The bgzf sink does its own block-level buffering.
const WRITE_BUFFER_BYTES: usize = 64 * 1024;

/// Where one record sits in the file's order: its contig, the position it is *written* at, and
/// whether it is a repeat tract.
///
/// **The written position, not the record's own span start.** A left-padded deletion is written
/// one base before its span (spec §5), and it is the written position a consumer sorts on — so
/// it is the one the order has to be checked against.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct RecordPlace {
    contig: ContigId,
    position: u64,
    is_repeat_tract: bool,
}

/// **The file's writer: a header, then records in genome order.**
///
/// Records must arrive non-decreasing in (contig, written position). Two records may share a
/// position, and exactly one shape of tie is legal — see [`Self::write_record`].
pub struct VcfWriter {
    sink: Sink,
    final_path: PathBuf,
    ploidy: Ploidy,
    metadata: VcfHeaderMetadata,
    last: Option<RecordPlace>,
    records_written: u64,
}

impl VcfWriter {
    /// Open the output and write the header.
    ///
    /// The path's suffix chooses the encoding: `.vcf.gz` or `.vcf.bgz` (either case) is bgzf,
    /// anything else is plain text. Bytes go to `<path>.tmp` until [`Self::finish`].
    ///
    /// # Errors
    ///
    /// If the temporary file cannot be created or the header cannot be written.
    pub fn create(
        final_path: &Path,
        metadata: VcfHeaderMetadata,
        ploidy: Ploidy,
    ) -> Result<Self, VcfWriteError> {
        let mut sink = Sink::open_tmp(final_path)?;
        let header = header_text(&metadata);
        sink.write_all(header.as_bytes())
            .map_err(|source| VcfWriteError::Write {
                tmp_path: tmp_path_for(final_path),
                source,
            })?;
        Ok(Self {
            sink,
            final_path: final_path.to_path_buf(),
            ploidy,
            metadata,
            last: None,
            records_written: 0,
        })
    }

    /// Write one record, checking it does not run backwards.
    ///
    /// **The order is non-decreasing, not strictly increasing, and the difference is exactly one
    /// case.** Production's generic writer demands strictly increasing positions and can, because
    /// it never shares a file with repeat-tract records. ng interleaves the two, and the padding
    /// rule creates one legal collision: a tract whose record moved one base left can land on the
    /// position of the generic locus that owns the anchor base. The two describe different
    /// bases — the reference partition guarantees it — so both belong in the file.
    ///
    /// So a tie is admitted **once**, and only generic-then-tract: the generic record's span
    /// genuinely starts there and the tract's starts one base later. A second record at the same
    /// position, a tract followed by a generic one, or two tracts, are refused.
    ///
    /// # Errors
    ///
    /// If the record runs backwards, if it forms a tie that is not the legal one, or if the
    /// write fails.
    pub fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfWriteError> {
        let place = place_of(record);
        if let Some(last) = self.last {
            self.check_order(last, place)?;
        }

        let line = record_line(record, self.metadata.contigs(), self.ploidy);
        self.sink
            .write_all(line.as_bytes())
            .and_then(|()| self.sink.write_all(b"\n"))
            .map_err(|source| VcfWriteError::Write {
                tmp_path: tmp_path_for(&self.final_path),
                source,
            })?;

        self.last = Some(place);
        self.records_written += 1;
        Ok(())
    }

    /// The ordering rule, in one place so the refusals cannot disagree with each other.
    fn check_order(&self, last: RecordPlace, next: RecordPlace) -> Result<(), VcfWriteError> {
        let going_backwards = (next.contig.0, next.position) < (last.contig.0, last.position);
        if going_backwards {
            return Err(VcfWriteError::OutOfOrder {
                previous_contig: last.contig.0,
                previous_position: last.position,
                contig: next.contig.0,
                position: next.position,
            });
        }
        let tied = next.contig == last.contig && next.position == last.position;
        if tied && !(!last.is_repeat_tract && next.is_repeat_tract) {
            return Err(VcfWriteError::IllegalTie {
                contig: next.contig.0,
                position: next.position,
                previous_was_repeat_tract: last.is_repeat_tract,
                is_repeat_tract: next.is_repeat_tract,
            });
        }
        Ok(())
    }

    /// **Write a whole stream of records, in the order it yields them.**
    ///
    /// This is the shape the run is built around: variants come off the caller, pass through
    /// filters and through the mappers that attach genotypes and annotations
    /// ([`assemble_record`](super::assemble::assemble_record) being the last of them), and end
    /// here. Ordering is still checked per record, so a filter or mapper that reorders the
    /// stream is refused rather than silently written.
    ///
    /// **Takes the records by value and one at a time**, so the stream can be lazy: nothing
    /// requires the whole cohort's records to exist at once, which is what keeps the writer off
    /// the memory budget at three thousand samples.
    ///
    /// # Errors
    ///
    /// The first record that runs backwards, forms an illegal tie, or cannot be written stops
    /// the stream and returns.
    pub fn write_stream(
        &mut self,
        records: impl IntoIterator<Item = VcfRecord>,
    ) -> Result<(), VcfWriteError> {
        for record in records {
            self.write_record(&record)?;
        }
        Ok(())
    }

    /// How many records have been written.
    #[must_use]
    pub fn records_written(&self) -> u64 {
        self.records_written
    }

    /// Flush everything, make it durable, and move the finished file into place.
    ///
    /// **Consumes the writer**, so a forgotten `finish` shows up as a missing output rather than
    /// as a silently truncated file.
    ///
    /// # Errors
    ///
    /// If any of the flush, the bgzf terminator, the `fsync`s or the rename fails.
    pub fn finish(self) -> Result<(), VcfWriteError> {
        self.sink.finish(&self.final_path)
    }
}

/// Where a record is written, after the padding rule has moved it.
fn place_of(record: &VcfRecord) -> RecordPlace {
    let start = record.region().start.get();
    let position = match record.padding_base() {
        // The record type refuses a left-hand padding base at position 1, so this cannot
        // underflow.
        Some(PaddingBase::Left(_)) => start - 1,
        Some(PaddingBase::Right(_)) | None => start,
    };
    RecordPlace {
        contig: record.region().contig,
        position,
        is_repeat_tract: record.is_repeat_tract(),
    }
}

/// Where the in-flight bytes live before the file is finished.
fn tmp_path_for(final_path: &Path) -> PathBuf {
    let mut name = final_path.as_os_str().to_os_string();
    name.push(".tmp");
    PathBuf::from(name)
}

/// Whether the path names a bgzf-compressed VCF, matched case-insensitively.
fn path_is_bgzf(path: &Path) -> bool {
    let name = path.to_string_lossy().to_lowercase();
    name.ends_with(".vcf.gz") || name.ends_with(".vcf.bgz")
}

/// The bytes' destination: plain text or bgzf, fixed at creation by the path's suffix.
enum Sink {
    Plain(BufWriter<File>),
    Bgzf(Box<noodles_bgzf::io::Writer<File>>),
}

impl Sink {
    fn open_tmp(final_path: &Path) -> Result<Self, VcfWriteError> {
        let tmp_path = tmp_path_for(final_path);
        let file = File::create(&tmp_path).map_err(|source| VcfWriteError::CreateTmp {
            tmp_path: tmp_path.clone(),
            source,
        })?;
        Ok(if path_is_bgzf(final_path) {
            Self::Bgzf(Box::new(noodles_bgzf::io::Writer::new(file)))
        } else {
            Self::Plain(BufWriter::with_capacity(WRITE_BUFFER_BYTES, file))
        })
    }

    /// Flush, terminate, `fsync`, rename, then `fsync` the parent directory.
    ///
    /// **The parent-directory sync is the part that is easy to leave out and hard to notice
    /// missing.** Without it a crash between the rename returning and the directory's journal
    /// reaching disk leaves the file's contents durable and the name pointing at nothing.
    fn finish(self, final_path: &Path) -> Result<(), VcfWriteError> {
        let tmp_path = tmp_path_for(final_path);

        let file = match self {
            Self::Plain(buffered) => {
                buffered
                    .into_inner()
                    .map_err(|error| VcfWriteError::Write {
                        tmp_path: tmp_path.clone(),
                        source: error.into_error(),
                    })?
            }
            Self::Bgzf(bgzf) => bgzf.finish().map_err(|source| VcfWriteError::FinishBgzf {
                tmp_path: tmp_path.clone(),
                source,
            })?,
        };

        file.sync_all().map_err(|source| VcfWriteError::Fsync {
            path: tmp_path.clone(),
            source,
        })?;
        drop(file);

        fs::rename(&tmp_path, final_path).map_err(|source| VcfWriteError::Rename {
            tmp_path: tmp_path.clone(),
            final_path: final_path.to_path_buf(),
            source,
        })?;

        sync_parent_directory(final_path)
    }
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(sink) => sink.write(buf),
            Self::Bgzf(sink) => sink.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(sink) => sink.flush(),
            Self::Bgzf(sink) => sink.flush(),
        }
    }
}

/// `fsync` the directory the output now lives in, so the rename survives a crash.
fn sync_parent_directory(final_path: &Path) -> Result<(), VcfWriteError> {
    let parent = final_path.parent().unwrap_or_else(|| Path::new("."));
    let directory = if parent.as_os_str().is_empty() {
        Path::new(".")
    } else {
        parent
    };
    File::open(directory)
        .and_then(|handle| handle.sync_all())
        .map_err(|source| VcfWriteError::Fsync {
            path: directory.to_path_buf(),
            source,
        })
}

/// What can go wrong writing the file. **All of it is about the file**, not about the records:
/// a record that could not be written was refused when it was built.
#[derive(Debug, Error)]
pub enum VcfWriteError {
    /// The temporary file could not be created.
    #[error("could not create the in-flight output `{tmp_path}`")]
    CreateTmp {
        /// The path that could not be created.
        tmp_path: PathBuf,
        /// What the filesystem said.
        source: io::Error,
    },

    /// A write to the temporary file failed.
    #[error("could not write to the in-flight output `{tmp_path}`")]
    Write {
        /// The path being written.
        tmp_path: PathBuf,
        /// What the filesystem said.
        source: io::Error,
    },

    /// The bgzf terminator could not be written.
    #[error("could not finish the bgzf stream in `{tmp_path}`")]
    FinishBgzf {
        /// The path being finished.
        tmp_path: PathBuf,
        /// What the filesystem said.
        source: io::Error,
    },

    /// A file or directory could not be made durable.
    #[error("could not flush `{path}` to disk")]
    Fsync {
        /// What was being synced.
        path: PathBuf,
        /// What the filesystem said.
        source: io::Error,
    },

    /// The finished file could not be moved into place.
    #[error("could not move `{tmp_path}` into place as `{final_path}`")]
    Rename {
        /// The in-flight path.
        tmp_path: PathBuf,
        /// The intended output path.
        final_path: PathBuf,
        /// What the filesystem said.
        source: io::Error,
    },

    /// **A record ran backwards.** A VCF's records are in genome order, and a consumer that
    /// indexes one would find the index disagreeing with the file rather than fail here.
    #[error(
        "records run in genome order and this one goes backwards: contig {contig} position \
         {position} after contig {previous_contig} position {previous_position}"
    )]
    OutOfOrder {
        /// The previous record's contig.
        previous_contig: u32,
        /// The previous record's written position.
        previous_position: u64,
        /// This record's contig.
        contig: u32,
        /// This record's written position.
        position: u64,
    },

    /// **Two records shared a position in a way the format does not allow.** The one legal tie
    /// is a generic locus followed by a repeat tract whose padding moved it onto that position.
    #[error(
        "two records share contig {contig} position {position}, and the only tie the format \
         allows is a SNP or indel followed by a repeat tract padded onto it — this is a {} \
         after a {}",
        if *is_repeat_tract { "repeat tract" } else { "SNP or indel" },
        if *previous_was_repeat_tract { "repeat tract" } else { "SNP or indel" }
    )]
    IllegalTie {
        /// The shared contig.
        contig: u32,
        /// The shared written position.
        position: u64,
        /// Whether the record already at this position was a repeat tract.
        previous_was_repeat_tract: bool,
        /// Whether the arriving record is a repeat tract.
        is_repeat_tract: bool,
    },
}

#[cfg(test)]
mod tests;
