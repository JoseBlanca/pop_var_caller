//! Replacing a finished psp's trailer, touching neither its blocks nor its index.
//!
//! **This is why the index sits before the trailer** (spec §3, §6.5): the trailer's offset is
//! where the rewrite starts, so everything before it stays exactly as it was and only the
//! trailer and the fixed tail are written. It is the cheap operation, and it exists because the
//! trailer is where things computed *after* the records land — the per-sample summary today,
//! whatever the statistical work adds later (spec §3.4).
//!
//! **The whole file is not rewritten and the writer is not involved.** A run that has already
//! spent an hour writing blocks does not spend another hour to change a histogram.

use std::fs::{File, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
use std::path::Path;

use super::footer::{Footer, decode_footer, encode_footer};
use super::reader::read_and_check_the_footer;
use super::{PspReadError, PspWriteError, read_header_from};

/// What a failed [`replace_trailer`] left on disk.
///
/// **The two need opposite actions from the caller**, which is the whole reason this exists: one
/// says *try again*, the other says *this sample has to be written from its reads again*. Before
/// this type they differed only by a phrase inside an [`PspWriteError::Io`] variant, which a
/// caller cannot branch on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum FileAfterAFailedReplacement {
    /// **Nothing was written.** The file is byte for byte what it was — the failure happened
    /// while reading it, or while checking what would be written back — so calling again is a
    /// retry and not a repair.
    #[error("the file is exactly as it was, and the call can be made again")]
    Unchanged,
    /// **The file was cut back to its trailer and the replacement did not land.** It has no
    /// footer, so every reader refuses it, and the trailer it used to hold is gone: there is
    /// nothing to retry from. The sample has to be written again.
    ///
    /// ⚠ **A failed *truncation* is reported here too, and that is the conservative side of a
    /// judgement.** `ftruncate` does not cut a file half-way, so in that one case the file is
    /// very likely untouched — but the failure that actually reaches it is a disk error, and a
    /// caller told *retry* about a file that a failing disk has cut goes on to trust it. The
    /// mistake this way costs one sample rewritten; the other way costs a sample silently short.
    #[error("the file is cut short and has no footer; the sample has to be written again")]
    Torn,
}

/// A [`replace_trailer`] that failed, and what it left behind.
///
/// **The state comes first because it is what the caller acts on**; the cause underneath says
/// what went wrong, and the two print as separate sentences so neither repeats the other.
#[derive(Debug, thiserror::Error)]
#[error("the trailer was not replaced: {file}")]
pub struct TrailerReplacementFailure {
    /// What the file is now — *the* thing to branch on.
    pub file: FileAfterAFailedReplacement,
    /// Why it failed.
    #[source]
    pub source: PspWriteError,
}

/// Replace a finished psp's trailer with `trailer`, leaving its blocks and its index untouched.
///
/// **The file must be finished.** A psp with no valid footer is one a run was killed part-way
/// through writing, and there is nothing to say where its trailer would go — so it is refused as
/// [`PspReadError::Incomplete`], wrapped in [`PspWriteError::Reopen`] (spec §6.7).
///
/// **⚠ The old trailer is gone the moment this starts writing.** The file is truncated at the
/// trailer's offset before a byte is written, so from that instant until the new footer lands
/// there is no footer and every reader refuses the file — which is right, and is what makes an
/// interruption safe. But it means this is **not** an in-place edit of a file whose *existing*
/// trailer matters: spec §6.4 gives append the same warning and the same answer, **write to a
/// new path and rename if that matters to the caller**.
///
/// ⚠ **Truncating first is the whole of that guarantee, and the first version did not.** It
/// overwrote in place and trimmed afterwards, so until a write passed the old trailer's end the
/// old footer was still there and still consistent — and the file opened, handing back a trailer
/// that was neither the old one nor the new one. Measured on the twenty-byte fixture: replacing
/// it with `short` and stopping left a file that opened with the trailer
/// `"short-sample summary"`. **Twenty of the twenty-one torn states were accepted**; with the
/// truncation first, one is — the complete write.
///
/// **What it checks.** The fixed tail, the two rules that need the file's length (that the
/// sections end exactly where the footer begins, and that the index does not start inside the
/// header's magic) — both shared with [`super::PspReader::open`] rather than copied — and the
/// header, for the one bound only the header can give: **that the trailer does not begin inside
/// the header or the blocks.**
///
/// ⚠ **That last check is the one the G3 review found missing, and its absence destroyed
/// files.** `decode_footer` proves the index ends where the trailer begins, and nothing bounded
/// either below except a four-byte magic — so a footer claiming `trailer_offset = 4` seeked to
/// byte 4 and wrote there, and a 3,742-byte psp that `PspReader::open` **already refused**
/// became 56 bytes and returned `Ok(())`.
///
/// **It does not read the block index**, and that is still right: reading it would cost a decode
/// per call, and every field written back except the trailer's length is the file's own.
///
/// **Which failures left the file alone, and which left it torn** — and the caller is told which
/// rather than having to infer it. Everything up to and including the footer read-back leaves the
/// file byte-identical, so the call can simply be made again; from the truncation onwards the old
/// trailer is gone and the file has no footer until the last write lands, so the sample has to be
/// written again. **The two halves are two functions and each labels its own failures**, so the
/// verdict is decided where the file's state is known rather than at the join, and a line moving
/// from one half to the other takes its verdict with it. See [`FileAfterAFailedReplacement`].
pub fn replace_trailer(path: &Path, trailer: &[u8]) -> Result<(), TrailerReplacementFailure> {
    let (mut file, replaced, trailer_offset) = read_what_the_replacement_needs(path, trailer)?;
    write_the_replacement(path, &mut file, &replaced, trailer_offset, trailer)
}

/// Open the file, check it, and work out the footer that will replace its own.
///
/// **Nothing here writes**, which is the property the caller's `Unchanged` verdict rests on: the
/// handle is opened without `truncate`, and every other step reads or computes. A failure leaves
/// the file byte for byte what it was.
fn read_what_the_replacement_needs(
    path: &Path,
    trailer: &[u8],
) -> Result<(File, Footer, u64), TrailerReplacementFailure> {
    // **Each half puts its own verdict on its own failures**, rather than the caller deciding
    // from which of the two returned. A verdict chosen at the join is a judgement made where
    // nothing knows the file's state; here it is made where everything does, and there is one
    // place per half to get it wrong.
    let unchanged = |source: PspWriteError| TrailerReplacementFailure {
        file: FileAfterAFailedReplacement::Unchanged,
        source,
    };
    let reopen = move |source: PspReadError| {
        unchanged(PspWriteError::Reopen {
            path: path.to_path_buf(),
            source,
        })
    };
    let io = move |while_doing: &'static str| {
        move |source: std::io::Error| {
            unchanged(PspWriteError::Io {
                path: path.to_path_buf(),
                while_doing,
                source,
            })
        }
    };
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .open(path)
        .map_err(io("reopening the file to replace its trailer"))?;
    let file_bytes = file.metadata().map_err(io("measuring the file"))?.len();

    let footer = read_and_check_the_footer(path, &mut file, file_bytes).map_err(reopen)?;
    let trailer_offset =
        check_the_trailer_is_past_the_header(path, &file, &footer).map_err(reopen)?;

    let replaced = Footer {
        trailer_bytes: trailer.len() as u64,
        ..footer
    };
    // **Read back before believing it**, the obligation `finish` carries for the same reason: a
    // footer that reaches disk and is then rejected has cost a file that was readable before.
    // Nothing here can make it fail — the fields come from a footer that already decoded and the
    // only one that moves is a length — and it is kept because that is an argument about today's
    // fields rather than about the encoding.
    //
    // **And it belongs on this side of the split**: a footer that would not decode is caught
    // before the file is cut, so it costs a refusal rather than a file.
    decode_footer(&encode_footer(&replaced)).map_err(|source| {
        unchanged(PspWriteError::WouldNotBeReadable {
            path: path.to_path_buf(),
            reason: "the footer it would write does not decode".to_string(),
            source: Some(source.into()),
        })
    })?;
    Ok((file, replaced, trailer_offset))
}

/// Cut the file back to its trailer and write the new trailer and footer.
///
/// **Every failure in here leaves a file no reader accepts**, which is what makes this the whole
/// of the `Torn` half. The first statement is the truncation, deliberately: from that instant
/// until the footer lands the file has none, so an interruption leaves something every reader
/// refuses instead of something that opens and lies.
fn write_the_replacement(
    path: &Path,
    file: &mut File,
    replaced: &Footer,
    trailer_offset: u64,
    trailer: &[u8],
) -> Result<(), TrailerReplacementFailure> {
    // **Every failure below is `Torn`, including the truncation's**, and there is one place that
    // says so — see [`FileAfterAFailedReplacement::Torn`] for why the truncation counts.
    let io = |while_doing: &'static str| {
        move |source: std::io::Error| TrailerReplacementFailure {
            file: FileAfterAFailedReplacement::Torn,
            source: PspWriteError::Io {
                path: path.to_path_buf(),
                while_doing,
                source,
            },
        }
    };
    file.set_len(trailer_offset)
        .map_err(io("cutting the file back to its trailer"))?;
    file.seek(SeekFrom::Start(trailer_offset))
        .map_err(io("seeking to the trailer"))?;
    file.write_all(trailer).map_err(io("writing the trailer"))?;
    file.write_all(&encode_footer(replaced))
        .map_err(io("writing the footer"))?;

    // **Durable, for the reason `finish` is** (spec §6.3): a caller that is told the trailer was
    // replaced and then loses power must not find the old one.
    file.sync_all().map_err(io("syncing the file"))?;
    Ok(())
}

/// The one bound only the header can give: **the trailer does not begin inside the header or the
/// blocks**.
///
/// The footer proves the index ends exactly where the trailer begins, and `open`'s own rule
/// bounds the index below by the four-byte magic — neither is enough. A file whose footer says
/// the trailer starts at byte 4 passes both, and rewriting it would put a trailer over the
/// header. So the header is read for its length alone, which is what the blocks begin after.
///
/// **This is the one thing this operation reads that spec §6.7's table does not account for**,
/// and it earns a refusal class that table does not list for it: a file written by a newer format
/// comes back as `UnsupportedVersion` rather than being rewritten. That is the safe answer and
/// the spec's table should gain the row.
fn check_the_trailer_is_past_the_header(
    path: &Path,
    file: &File,
    footer: &Footer,
) -> Result<u64, PspReadError> {
    // **From the handle already open**, the way `open` reads its own header.
    let for_the_header = file.try_clone().map_err(|source| PspReadError::Io {
        path: path.to_path_buf(),
        while_doing: "reading the header",
        source,
    })?;
    let (_, header_bytes) = read_header_from(for_the_header, path)?;
    if footer.trailer_offset < header_bytes as u64 {
        return Err(PspReadError::damaged(
            path,
            format!(
                "the footer puts the trailer at byte {}, which is inside the {header_bytes}-byte \
                 header",
                footer.trailer_offset
            ),
        ));
    }
    Ok(footer.trailer_offset)
}

#[cfg(test)]
mod tests {
    use super::super::footer::FOOTER_BYTES;
    use super::*;
    use crate::ng::psp::PspReader;
    use crate::ng::psp::writer::PspWriter;
    use crate::ng::psp::writer::tests_support::{
        a_file, a_finished_psp, a_header, a_sample, bytes_of, footer_of, rewrite,
    };

    /// **A longer trailer, a shorter one and an empty one all read back**, and each leaves a file
    /// whose sections still account for every byte.
    ///
    /// The shorter case is the one that needs the file trimmed: without it the tail of the old
    /// trailer sits past the new footer, and a reader taking the last forty-eight bytes finds
    /// something that is not a footer.
    #[test]
    fn a_trailer_of_any_length_replaces_the_one_before_it() {
        let (_dir, path) = a_finished_psp();
        let first = {
            let mut psp = PspReader::open(&path).expect("a finished psp opens");
            psp.trailer().expect("the trailer reads")
        };
        assert_eq!(first, b"a per-sample summary");

        for payload in [
            b"a much longer per-sample summary than the one before it".as_slice(),
            b"short".as_slice(),
            b"".as_slice(),
            b"and one more".as_slice(),
        ] {
            replace_trailer(&path, payload).expect("the trailer is replaced");
            let mut psp = PspReader::open(&path).expect("the file is still a finished psp");
            assert_eq!(psp.trailer().expect("the trailer reads"), payload);
            let whole = bytes_of(&path);
            let footer = footer_of(&whole);
            assert_eq!(
                footer.trailer_offset + footer.trailer_bytes,
                (whole.len() - FOOTER_BYTES) as u64,
                "the sections end exactly where the footer begins"
            );
        }
    }

    /// **The blocks and the index are untouched, byte for byte**, which is the whole reason the
    /// index sits before the trailer (spec §3, §6.5).
    ///
    /// Compared as bytes rather than as a decode: an index that happened to decode the same way
    /// after being rewritten would pass a weaker check.
    #[test]
    fn replacing_a_trailer_leaves_the_blocks_and_the_index_byte_identical() {
        let (_dir, path) = a_finished_psp();
        let before = bytes_of(&path);
        let footer = footer_of(&before);
        let up_to_the_trailer = before[..footer.trailer_offset as usize].to_vec();

        replace_trailer(&path, b"something else entirely, and longer").expect("it is replaced");

        let after = bytes_of(&path);
        assert_eq!(
            after[..footer.trailer_offset as usize],
            up_to_the_trailer[..],
            "the header, every block and the index are the same bytes"
        );
        // And the records still come back.
        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        let read = psp
            .records()
            .expect("the walk starts")
            .map(|found| found.expect("a finished psp walks").record.expect("a body"))
            .collect::<Vec<_>>();
        assert_eq!(read, a_sample());
    }

    /// **Only the trailer's length changes in the footer.** Everything else it carries describes
    /// the blocks and the index, which this operation does not touch — so writing a different
    /// value there would be inventing one.
    #[test]
    fn only_the_trailers_length_changes_in_the_footer() {
        let (_dir, path) = a_finished_psp();
        let before = footer_of(&bytes_of(&path));
        replace_trailer(&path, b"a summary of a different length").expect("it is replaced");
        let after = footer_of(&bytes_of(&path));

        assert_eq!(after.trailer_bytes, 31);
        assert_ne!(after.trailer_bytes, before.trailer_bytes);
        assert_eq!(
            Footer {
                trailer_bytes: before.trailer_bytes,
                ..after
            },
            before,
            "every other field is the file's own"
        );
    }

    /// **A failure once the file has been cut is reported as torn, and the two halves are told
    /// apart by which function failed rather than by how far a single one got.**
    ///
    /// The two verdicts ask the caller for opposite things — call again, or write the sample
    /// again — so the one thing worth testing is that a failure on the writing side is never
    /// reported as *unchanged*. Every other test here covers the reading side, and asserts
    /// `Unchanged` beside the file's own bytes.
    ///
    /// **The write is made to fail with a handle opened read-only**, which is the same trick the
    /// walk uses in reverse to make a `read(2)` fail on a sound file: no `unsafe`, a real refusal
    /// from the kernel, and the psp underneath is a good one. **What it cannot do is reach
    /// `replace_trailer` itself**, which opens the file for writing — so it drives the writing
    /// half directly. That is enough because each half now puts the verdict on its own failures:
    /// the value asserted below is the one a caller would receive, not one built here.
    #[test]
    fn a_failure_after_the_file_is_cut_is_reported_as_torn() {
        let (_dir, path) = a_finished_psp();
        let before = bytes_of(&path);

        // A sound psp, read for everything the replacement needs — this half must succeed, or
        // the test would be proving something about the reading side.
        let (_writable, replaced, trailer_offset) =
            read_what_the_replacement_needs(&path, b"a different summary")
                .expect("the fixture is a psp this operation accepts");

        let mut read_only = File::open(&path).expect("the file opens for reading");
        let refused = write_the_replacement(
            &path,
            &mut read_only,
            &replaced,
            trailer_offset,
            b"a different summary",
        )
        .expect_err("a handle with no write permission cannot cut or write the file");

        assert_eq!(
            refused.file,
            FileAfterAFailedReplacement::Torn,
            "a failure on the writing side is never *unchanged*, whatever the caller does next"
        );
        assert!(
            matches!(refused.source, PspWriteError::Io { .. }),
            "a refused syscall, not a decode: {}",
            refused.source
        );
        assert!(
            refused.to_string().contains("written again"),
            "the sentence has to tell the caller what to do; got {refused}"
        );
        assert!(
            !refused.to_string().contains("call can be made again"),
            "and must not tell it the opposite; got {refused}"
        );

        assert_eq!(
            bytes_of(&path),
            before,
            "nothing wrote here — the file is the fixture's, and the verdict is about what a \
             failure on this side means, not about this particular failure"
        );
    }

    /// **A file with no footer has no trailer to replace**, and is refused as incomplete — the
    /// class spec §6.7 names for this operation. It is what a killed run leaves.
    #[test]
    fn a_file_with_no_footer_is_refused_as_incomplete() {
        let (_dir, path) = a_file();
        let mut writer = PspWriter::create(&path, a_header(1_000)).expect("a header");
        for record in a_sample() {
            writer.push(&record).expect("in order");
        }
        drop(writer); // killed before `finish`: header and blocks, no footer

        let refused = replace_trailer(&path, b"anything").expect_err("there is no trailer");
        assert_eq!(
            refused.file,
            FileAfterAFailedReplacement::Unchanged,
            "nothing was written, so the caller may simply call again"
        );
        match refused.source {
            PspWriteError::Reopen {
                source: PspReadError::Incomplete { .. },
                ..
            } => {}
            other => panic!("got {other}"),
        }
    }

    /// **Sections that stop short of the footer and sections that run past it are both
    /// refused**, and the file is left alone either way.
    ///
    /// ⚠ **Nudged one way only, this could not fail on the half it did not name.** Weakening the
    /// rule from an equality to a *greater-than* — which accepts every file whose sections stop
    /// short — left all 370 tests green. That is the F4 Blocker exactly, in the rule this
    /// operation shares with `open`, one commit after it was fixed there.
    #[test]
    fn sections_that_stop_short_or_run_past_the_footer_are_both_refused() {
        for (nudge, what) in [(1i64, "run past"), (-1, "stop short")] {
            let (_dir, path) = a_finished_psp();
            let mut whole = bytes_of(&path);
            let at = whole.len() - FOOTER_BYTES;
            let mut footer = footer_of(&whole);
            footer.trailer_bytes = footer.trailer_bytes.wrapping_add(nudge as u64);
            whole[at..].copy_from_slice(&encode_footer(&footer));
            rewrite(&path, &whole);
            let before = bytes_of(&path);

            let refused = replace_trailer(&path, b"anything")
                .expect_err("the sections must account for the file");
            assert_eq!(
                refused.file,
                FileAfterAFailedReplacement::Unchanged,
                "nothing was written, so the caller may simply call again"
            );
            match refused.source {
                PspWriteError::Reopen {
                    source: PspReadError::Damaged { reason, .. },
                    ..
                } => assert!(
                    reason.contains("sections end at byte"),
                    "sections that {what} the footer must be refused: {reason}"
                ),
                other => panic!("sections that {what} the footer: got {other}"),
            }
            assert_eq!(bytes_of(&path), before, "and nothing was written to it");
        }
    }

    /// **A file shorter than a footer is refused as incomplete**, which is the class spec §6.7
    /// gives this operation, and the file is left alone.
    #[test]
    fn a_file_shorter_than_a_footer_is_refused_as_incomplete() {
        let (_dir, path) = a_file();
        let too_short = vec![0x5a; FOOTER_BYTES - 1];
        rewrite(&path, &too_short);

        let refused = replace_trailer(&path, b"anything").expect_err("there is no footer");
        assert_eq!(
            refused.file,
            FileAfterAFailedReplacement::Unchanged,
            "nothing was written, so the caller may simply call again"
        );
        match refused.source {
            PspWriteError::Reopen {
                source: PspReadError::Incomplete { .. },
                ..
            } => {}
            other => panic!("got {other}"),
        }
        assert_eq!(bytes_of(&path), too_short, "the file is exactly as it was");
    }

    /// **A file that was never a psp is told apart from a killed run**, and nothing is written
    /// to it.
    ///
    /// ⚠ **The first version of this test was 47 bytes long**, one short of a footer, so it
    /// returned before the magic was ever looked at: it named the foreign-file path and
    /// exercised the short-file one.
    ///
    /// ⚠ **And it asserted the wrong class.** This operation used to map every unreadable tail
    /// to *incomplete*, reasoning that a killed run and a foreign file both mean there is no
    /// trailer to replace. Sharing `open`'s footer read brought `open`'s distinction with it,
    /// and the distinction is worth having: *you handed me the wrong file* and *rebuild this
    /// one* are different instructions, and only one of them is answered by re-running a
    /// pileup.
    #[test]
    fn a_foreign_file_longer_than_a_footer_is_told_apart_from_a_killed_run() {
        let (_dir, path) = a_file();
        let not_a_psp = b"this is a text file, it is not a psp at all, and it is comfortably \
                          longer than the forty-eight bytes a footer occupies."
            .to_vec();
        assert!(
            not_a_psp.len() > FOOTER_BYTES,
            "or it takes the short-file path"
        );
        rewrite(&path, &not_a_psp);

        let refused = replace_trailer(&path, b"anything").expect_err("that is not a psp");
        assert_eq!(
            refused.file,
            FileAfterAFailedReplacement::Unchanged,
            "nothing was written, so the caller may simply call again"
        );
        match refused.source {
            PspWriteError::Reopen {
                source: PspReadError::NotAnNgPsp { found, .. },
                ..
            } => assert_eq!(&found, b"this", "the file's own first four bytes"),
            other => panic!("got {other}"),
        }
        assert_eq!(bytes_of(&path), not_a_psp, "the file is exactly as it was");
    }

    /// **A file that is not there names itself and what was being done to it**, rather than
    /// arriving as a bare `No such file or directory`.
    #[test]
    fn a_missing_file_names_itself() {
        let (_dir, path) = a_file();
        let refused = replace_trailer(&path, b"anything").expect_err("there is no file");
        assert_eq!(
            refused.file,
            FileAfterAFailedReplacement::Unchanged,
            "nothing was written, so the caller may simply call again"
        );
        match refused.source {
            PspWriteError::Io {
                while_doing,
                path: named,
                ..
            } => {
                assert_eq!(named, path);
                assert_eq!(while_doing, "reopening the file to replace its trailer");
            }
            other => panic!("got {other}"),
        }
    }

    /// **Replacing a trailer twice leaves the second one**, and the file does not grow each
    /// time — the length is set from the trailer's offset, not from what was there before.
    #[test]
    fn replacing_a_trailer_twice_does_not_grow_the_file() {
        let (_dir, path) = a_finished_psp();
        replace_trailer(&path, b"the first replacement").expect("it is replaced");
        let once = bytes_of(&path).len();
        replace_trailer(&path, b"the first replacement").expect("it is replaced again");
        let twice = bytes_of(&path).len();
        assert_eq!(once, twice, "the same payload gives the same file length");

        let mut psp = PspReader::open(&path).expect("a finished psp opens");
        assert_eq!(
            psp.trailer().expect("the trailer reads"),
            b"the first replacement"
        );
    }

    /// **A footer that puts the block index inside the header is refused.** The footer's own
    /// decode cannot see it — 48 bytes do not know where the header ends — and this operation
    /// does not read the header, so this rule is the only thing standing between a file whose
    /// offsets are nonsense and a fresh footer blessing them.
    #[test]
    fn a_footer_that_puts_the_index_inside_the_header_is_refused() {
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let footer_at = whole.len() - FOOTER_BYTES;
        let mut footer = footer_of(&whole);
        // The index moved to byte 0 and stretched to reach the trailer, so the footer's own
        // agreements still hold and only the length-aware rule can refuse it.
        footer.index_bytes += footer.index_offset;
        footer.index_offset = 0;
        whole[footer_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &whole);

        let refused = replace_trailer(&path, b"anything").expect_err("byte 0 is the head magic");
        assert_eq!(
            refused.file,
            FileAfterAFailedReplacement::Unchanged,
            "nothing was written, so the caller may simply call again"
        );
        match refused.source {
            PspWriteError::Reopen {
                source: PspReadError::Damaged { reason, .. },
                ..
            } => assert!(reason.contains("inside the header"), "got {reason}"),
            other => panic!("got {other}"),
        }
        assert_eq!(
            bytes_of(&path).len(),
            whole.len(),
            "and nothing was written to it"
        );
    }

    /// **A footer that puts the trailer inside the header is refused, and the file survives.**
    ///
    /// ⚠ **This is the G3 review's Blocker, and it destroyed files.** `decode_footer` proves the
    /// index ends exactly where the trailer begins, and the only lower bound on either was a
    /// four-byte magic — so a footer claiming the index and the trailer both start at byte 4
    /// passes every check the footer can make about itself, and the rewrite seeked to byte 4 and
    /// wrote there. **A 3,742-byte psp that `PspReader::open` already refuses became 56 bytes,
    /// and the call returned `Ok(())`.**
    #[test]
    fn a_footer_that_puts_the_trailer_inside_the_header_is_refused() {
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let footer_at = whole.len() - FOOTER_BYTES;
        let mut footer = footer_of(&whole);
        footer.index_offset = 4;
        footer.index_bytes = 0;
        footer.trailer_offset = 4;
        footer.trailer_bytes = footer_at as u64 - 4;
        whole[footer_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &whole);
        let before = bytes_of(&path);
        assert!(
            PspReader::open(&path).is_err(),
            "the fixture is a file no reader accepts, which is the point"
        );

        let refused = replace_trailer(&path, b"tiny")
            .expect_err("byte 4 is inside the header, not a trailer");
        assert_eq!(
            refused.file,
            FileAfterAFailedReplacement::Unchanged,
            "nothing was written, so the caller may simply call again"
        );
        match refused.source {
            PspWriteError::Reopen {
                source: PspReadError::Damaged { reason, .. },
                ..
            } => assert!(reason.contains("inside the"), "got {reason}"),
            other => panic!("got {other}"),
        }
        assert_eq!(
            bytes_of(&path),
            before,
            "and every byte of it is still there"
        );
    }

    /// **An interruption between the truncation and the footer leaves a file every reader
    /// refuses**, whatever the new trailer's length.
    ///
    /// ⚠ **The first version of this operation overwrote in place and trimmed afterwards, and
    /// this test is what that could not survive.** Until a write passed the old trailer's end
    /// the old footer was still there and still consistent, so the file opened — handing back a
    /// trailer that was neither the old one nor the new one. Replacing `a per-sample summary`
    /// with `short` and stopping gave `short-sample summary`. Every stopping point is walked
    /// here, at every payload length that matters.
    #[test]
    fn a_rewrite_stopped_part_way_leaves_a_file_no_reader_accepts() {
        for payload in [
            b"".as_slice(),
            b"short".as_slice(),
            b"a per-sample summar!".as_slice(), // the same length as the one it replaces
            b"a much longer per-sample summary than before".as_slice(),
        ] {
            let footer = {
                let (_dir, path) = a_finished_psp();
                footer_of(&bytes_of(&path))
            };
            let whole_len = footer.trailer_offset as usize + payload.len() + FOOTER_BYTES;
            for stopped_after in 0..=whole_len - footer.trailer_offset as usize {
                let (_dir, path) = a_finished_psp();
                // The file as it stands the instant the rewrite is interrupted: cut back to the
                // trailer's offset, then as many of the new bytes as had landed.
                let mut torn = bytes_of(&path)[..footer.trailer_offset as usize].to_vec();
                let replaced = Footer {
                    trailer_bytes: payload.len() as u64,
                    ..footer
                };
                let rest: Vec<u8> = payload
                    .iter()
                    .copied()
                    .chain(encode_footer(&replaced))
                    .collect();
                torn.extend_from_slice(&rest[..stopped_after.min(rest.len())]);
                rewrite(&path, &torn);

                let whole = stopped_after >= rest.len();
                match PspReader::open(&path) {
                    Ok(mut psp) => {
                        assert!(
                            whole,
                            "a rewrite stopped {stopped_after} bytes in opened, with the trailer \
                             {:?}",
                            String::from_utf8_lossy(&psp.trailer().expect("it reads"))
                        );
                        assert_eq!(psp.trailer().expect("it reads"), payload);
                    }
                    Err(_) => assert!(!whole, "the finished rewrite must open"),
                }
            }
        }
    }
}
