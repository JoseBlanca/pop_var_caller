//! **The spill as a file on disk: where it lives, when it appears, and what makes it go away.**
//!
//! [`spill`](super::spill) says what an entry is and how it turns into bytes. This says where
//! those bytes go: a file beside the output, written once in pass one and read twice
//! afterwards, and **gone when the run ends whatever way it ends** — success, an error on its
//! way up, or a panic unwinding through the caller.
//!
//! **The removal is a `Drop`, and that is the whole point.** A guard called on the success path
//! is a guard that is not called on the two paths where a leftover file is most likely, and
//! this file is the size of the VCF. Nothing has to remember to clean up: the value going out
//! of scope is the cleanup, so an early `?` and an unwinding panic remove the file for the same
//! reason a normal return does.
//!
//! **A run killed from outside still leaves the file**, exactly as it leaves `<output>.tmp`
//! ([`vcf::writer`](crate::vcf::writer)); nothing inside a process runs when the process is
//! killed. That is the spec's stated position (§3.4), not an oversight.
//!
//! **The file appears on the first record and not before**, so a run whose filter is off names
//! nothing and leaves nothing. A run *with* the filter on that called no records still ends
//! pass one with an empty file, because pass one ran: that is a different thing from a spill
//! that was never opened, and it is what lets passes two and three read an empty stream rather
//! than meet a missing file.
//!
//! **The path is fixed by spec §3.4, so it can collide.** Anything already there — a leftover
//! from a killed run, or another run writing the same output — is named rather than truncated:
//! one of the two is evidence and the other is live.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.4, §5.

use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::time::Instant;

use thiserror::Error;

use super::spill::{SpillEntry, SpillError, SpillReader, SpillWriter};

/// How many of a record's bytes sit in memory before they reach the kernel, on each side.
///
/// **The two sides need it for different reasons.** The writer issues one `write_all` per
/// record, so unbuffered it would cost a write syscall per called record. The reader takes the
/// variable-length integers a byte at a time, so unbuffered it would cost a read syscall per
/// *byte*. 64 KiB is the VCF writer's buffer
/// ([`vcf::writer`](crate::vcf::writer)), and one entry fits inside it up to about six
/// thousand samples — past that a record spans two flushes, which costs nothing but a second
/// write.
const BUFFER_BYTES: usize = 64 * 1024;

/// **The spill file: `<output>.paralog-spill.tmp`, beside the output on the same filesystem.**
///
/// One of these exists for one run. Pass one appends to it; [`Self::finish_writing`] flushes it;
/// passes two and three each open their own cursor over it with [`Self::read`]. Dropping it
/// unlinks the file.
///
/// The name follows the convention the VCF writer already uses for `<output>.tmp`: beside the
/// output rather than in a scratch directory anywhere else, so that the spill lands on the
/// filesystem the operator chose for the run's output and a leftover is found next to the run
/// it belongs to. **It is not the same size as the output**: the spill holds each record's line
/// uncompressed plus about ten bytes a sample, so against a `.vcf.gz` output it is several
/// times larger.
pub struct SpillFile {
    path: PathBuf,
    /// When this spill was named, which is just before the calling loop starts — **not** when its
    /// first record arrived, since naming one creates nothing.
    ///
    /// **Plan step D2 asks for the wall each pass spent**, and the calling pass is the only one
    /// that runs outside `fit_score_and_write_the_calls`; taking the clock here rather than in
    /// each subcommand keeps one clock instead of two copies of it in the two copies of the
    /// wiring that already differ only by an error type.
    named_at: Instant,
    /// Bytes on disk when pass one finished, or `None` before that. **What it is for**: the spill
    /// holds every record's line uncompressed plus about ten bytes a sample, so it is several
    /// times a compressed output and nothing said by how much (B2 review M6).
    bytes_on_disk: Option<u64>,
    /// The open writer, from the first appended entry until [`Self::finish_writing`].
    writer: Option<SpillWriter<BufWriter<File>>>,
    /// How far through its life the file is, which is what says whether there is anything to
    /// unlink and whether another record may still be appended.
    stage: Stage,
    entries_written: u64,
}

/// Where a spill is in its one life: named, being written, or finished.
///
/// **The three are told apart by this and not by asking the filesystem.** A file at the path
/// may be a leftover from a run that was killed, which must be named rather than silently
/// truncated, and `append` runs once per called record — a `stat` per record would be a syscall
/// per record for an answer this type already knows.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Stage {
    /// Named, nothing on disk. A run with the filter off never leaves this state.
    NotYetCreated,
    /// Pass one is running: the file exists and the writer is open.
    BeingWritten,
    /// Pass one has ended: the file exists, is flushed, and takes no more records.
    Finished,
}

impl SpillFile {
    /// Name the spill that would sit beside `output`. **Creates nothing** — the file appears
    /// when the first entry is appended.
    #[must_use]
    pub fn beside(output: &Path) -> Self {
        Self {
            path: spill_path_beside(output),
            named_at: Instant::now(),
            bytes_on_disk: None,
            writer: None,
            stage: Stage::NotYetCreated,
            entries_written: 0,
        }
    }

    /// Where the file is, or would be.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// How many entries have been appended.
    #[must_use]
    pub fn entries_written(&self) -> u64 {
        self.entries_written
    }

    /// Append one entry, creating the file if this is the first.
    ///
    /// # Errors
    ///
    /// If pass one has already ended, if something is already at the path, if the file cannot
    /// be created, or if the entry cannot be written. All name the path.
    pub fn append(&mut self, entry: &SpillEntry) -> Result<(), SpillFileError> {
        let writer = match self.stage {
            Stage::NotYetCreated => {
                let opened = self.create()?;
                self.stage = Stage::BeingWritten;
                self.writer.insert(opened)
            }
            // PANIC-FREE: `BeingWritten` is set only where the writer is stored just above, and
            // only `finish_writing` takes it back out — and that moves the stage to `Finished`.
            Stage::BeingWritten => self
                .writer
                .as_mut()
                .expect("a spill being written holds its writer"),
            Stage::Finished => {
                return Err(SpillFileError::AppendedAfterPassOneEnded {
                    path: self.path.clone(),
                });
            }
        };
        writer
            .append(entry)
            .map_err(|source| SpillFileError::Write {
                path: self.path.clone(),
                source,
            })?;
        self.entries_written += 1;
        Ok(())
    }

    /// End pass one: flush what has been appended, so that [`Self::read`] sees all of it.
    ///
    /// **Idempotent, and required before reading.**
    ///
    /// **On a run that appended nothing this creates the file, empty.** Pass one ran, so there
    /// is a spill and it holds no records — which is a different thing from a spill that was
    /// never opened, and it lets [`Self::read`] hand back an empty stream rather than a missing
    /// file. A run with the filter off never reaches this and leaves nothing on disk.
    ///
    /// # Errors
    ///
    /// If the file cannot be created, or the flush fails. Both name the path.
    pub fn finish_writing(&mut self) -> Result<(), SpillFileError> {
        let writer = match self.stage {
            Stage::Finished => return Ok(()),
            Stage::NotYetCreated => self.create()?,
            // PANIC-FREE: as in `append` — the writer is present for exactly this stage.
            Stage::BeingWritten => self
                .writer
                .take()
                .expect("a spill being written holds its writer"),
        };
        // The stage moves before the flush can fail, because the file exists either way and a
        // failed flush must not leave the spill looking as though pass one never opened it.
        self.stage = Stage::Finished;
        let buffered = writer.finish().map_err(|source| SpillFileError::Flush {
            path: self.path.clone(),
            source,
        })?;
        // Dropping the `BufWriter` closes the handle. It has just been flushed, so its own
        // `Drop` has nothing left to swallow.
        drop(buffered);
        // **After the flush, so it is the whole file and not what had reached the disk.** A
        // failure here is not the run's problem — the spill is written and readable, and the size
        // is a line in a report — so it leaves `None` rather than ending the run.
        self.bytes_on_disk = fs::metadata(&self.path).ok().map(|it| it.len());
        Ok(())
    }

    /// When this spill was named — the start of the calling pass, and before its first record.
    #[must_use]
    pub fn named_at(&self) -> Instant {
        self.named_at
    }

    /// How large the spill grew, once pass one has finished; `None` before that, or where the
    /// size could not be read.
    #[must_use]
    pub fn bytes_on_disk(&self) -> Option<u64> {
        self.bytes_on_disk
    }

    /// Open a cursor over the file, from the beginning. **May be called more than once** — pass
    /// two scores the entries and pass three writes them, and each gets its own cursor.
    ///
    /// The reader knows how many entries were written, so a file that ends early — a lost tail,
    /// a partial final write — is refused instead of read as a complete, shorter spill.
    ///
    /// # Errors
    ///
    /// If [`Self::finish_writing`] has not been called — reading before the flush would see
    /// whatever had reached the disk — or if the file cannot be opened. Both name the path.
    pub fn read(&self) -> Result<SpillReader<BufReader<File>>, SpillFileError> {
        if self.stage != Stage::Finished {
            return Err(SpillFileError::PassOneHasNotEnded {
                path: self.path.clone(),
            });
        }
        let file = File::open(&self.path).map_err(|source| SpillFileError::Open {
            path: self.path.clone(),
            source,
        })?;
        Ok(SpillReader::new(
            BufReader::with_capacity(BUFFER_BYTES, file),
            self.entries_written,
        ))
    }

    /// Give a failure that came off one of this file's readers the path it happened on.
    ///
    /// The reader is generic over its stream and carries no path, by design; spec §5 asks that
    /// a spill read failure name the file, and this is where the name is.
    #[must_use]
    pub fn naming_this_file(&self, source: SpillError) -> SpillFileError {
        SpillFileError::Read {
            path: self.path.clone(),
            source,
        }
    }

    /// Create the file, refusing to overwrite anything already at the path.
    ///
    /// **`create_new` rather than `truncate`**, because everything that could be there is
    /// something this run must not destroy: a spill left by a run that was killed, or one
    /// another `SpillFile` for the same output is still writing. Truncating either loses a
    /// finished run's evidence or half of a live run's records, and returns `Ok` while doing
    /// it.
    fn create(&self) -> Result<SpillWriter<BufWriter<File>>, SpillFileError> {
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&self.path)
            .map_err(|source| {
                if source.kind() == std::io::ErrorKind::AlreadyExists {
                    SpillFileError::SomethingIsAlreadyThere {
                        path: self.path.clone(),
                    }
                } else {
                    SpillFileError::Create {
                        path: self.path.clone(),
                        source,
                    }
                }
            })?;
        Ok(SpillWriter::new(BufWriter::with_capacity(
            BUFFER_BYTES,
            file,
        )))
    }
}

impl Drop for SpillFile {
    /// Close the handle, then unlink the file.
    ///
    /// **Nothing is reported if the removal fails.** A `Drop` that runs while a panic is
    /// unwinding cannot return an error and must not panic itself, and the failure it would be
    /// reporting — a temporary file left on disk — is the same one a killed process leaves. The
    /// path is `<output>.paralog-spill.tmp`, so a leftover names the run it came from.
    fn drop(&mut self) {
        let Self {
            path,
            writer,
            stage,
            named_at: _,
            bytes_on_disk: _,
            entries_written: _,
        } = self;
        // Closing the handle before unlinking is not needed on Unix and is on Windows.
        *writer = None;
        if *stage == Stage::NotYetCreated {
            return;
        }
        if let Err(source) = fs::remove_file(&*path) {
            // Silently leaking a file the size of the VCF is worse than a line on stderr, and
            // the crate already answers this the same way in `reference_info`'s
            // `VerificationHandle`. Stay quiet while a panic is unwinding: a second message
            // there buries the first, and a `Drop` must not panic.
            if !std::thread::panicking() {
                eprintln!(
                    "warning: the paralog spill `{}` could not be removed and is still on \
                     disk: {source}",
                    path.display()
                );
            }
        }
    }
}

/// What can go wrong with the file itself, as opposed to with a record in it.
///
/// **Every variant names the path**, which is what spec §5 asks of a spill failure: a run over
/// thousands of samples that says only "the spill failed" leaves nobody anywhere to look. The
/// causes are [`SpillError`]s and `io::Error`s, and both are carried as sources so the whole
/// chain renders.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SpillFileError {
    /// The file could not be created on the first record.
    #[error("the paralog spill `{}` could not be created", path.display())]
    Create {
        /// Where it was to be created.
        path: PathBuf,
        /// What the filesystem said.
        source: std::io::Error,
    },

    /// Something is already at the path — a spill left by a run that was killed, or one
    /// another run for the same output is writing now.
    ///
    /// **Refused rather than overwritten.** The name is fixed by spec §3.4, so two runs writing
    /// one output collide on it; truncating would destroy a live run's records, and a leftover
    /// is evidence about a run that failed. Remove the file named here, or give the run its own
    /// output path.
    #[error(
        "something is already at the paralog spill's path `{}`: remove it, or give this run \
         its own output path",
        path.display()
    )]
    SomethingIsAlreadyThere {
        /// What is in the way.
        path: PathBuf,
    },

    /// A record was appended after pass one had ended.
    ///
    /// **A caller's ordering mistake**: the file is flushed and its record count is settled, so
    /// a further record would either be lost or, if the file were reopened, replace what is
    /// there.
    #[error(
        "a record was appended to the paralog spill `{}` after pass one had ended",
        path.display()
    )]
    AppendedAfterPassOneEnded {
        /// The file that is no longer taking records.
        path: PathBuf,
    },

    /// A record could not be written.
    #[error("a record could not be written to the paralog spill `{}`", path.display())]
    Write {
        /// The file being written.
        path: PathBuf,
        /// What went wrong.
        #[source]
        source: SpillError,
    },

    /// The file could not be flushed, so what it holds is not what was appended.
    #[error("the paralog spill `{}` could not be flushed to disk", path.display())]
    Flush {
        /// The file being flushed.
        path: PathBuf,
        /// What went wrong.
        #[source]
        source: SpillError,
    },

    /// The file could not be opened for reading.
    #[error("the paralog spill `{}` could not be opened for reading", path.display())]
    Open {
        /// The file being opened.
        path: PathBuf,
        /// What the filesystem said.
        source: std::io::Error,
    },

    /// A record could not be read back.
    ///
    /// The reader is generic over its stream and carries no path; this is
    /// [`SpillFile::naming_this_file`]'s doing, so that spec §5's "name the file" holds on the
    /// read side as well as the write side.
    #[error("a record could not be read from the paralog spill `{}`", path.display())]
    Read {
        /// The file being read.
        path: PathBuf,
        /// What went wrong.
        #[source]
        source: SpillError,
    },

    /// A read was asked for before pass one had ended.
    ///
    /// **This is a caller's ordering mistake, not a filesystem failure.** Reading before the
    /// flush would see whatever had reached the disk — a spill short by its last records, and
    /// short in a way nothing downstream could notice — and reading a spill that was never
    /// opened would report a missing file where the answer is "pass one has not finished".
    #[error(
        "the paralog spill `{}` has not been finished; call finish_writing before reading it",
        path.display()
    )]
    PassOneHasNotEnded {
        /// The file that is not ready to be read.
        path: PathBuf,
    },
}

/// `<output>.paralog-spill.tmp`, the convention [`vcf::writer`](crate::vcf::writer) already
/// uses for `<output>.tmp`.
///
/// Appended to the whole path rather than substituted into it, so an output with no extension,
/// several extensions, or a name that looks like a directory all give one predictable answer.
fn spill_path_beside(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_os_string();
    name.push(".paralog-spill.tmp");
    PathBuf::from(name)
}

#[cfg(test)]
mod tests;
