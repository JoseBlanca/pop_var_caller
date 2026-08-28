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

use std::fs::OpenOptions;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;

use super::footer::{FOOTER_BYTES, Footer, decode_footer, encode_footer};
use super::header::HEAD_MAGIC;
use super::{PspReadError, PspWriteError};

/// Replace a finished psp's trailer with `trailer`, leaving its blocks and its index untouched.
///
/// **The file must be finished.** A psp with no valid footer is one a run was killed part-way
/// through writing, and there is nothing to say where its trailer would go — so it is refused as
/// [`PspReadError::Incomplete`], wrapped in [`PspWriteError::Reopen`] (spec §6.7).
///
/// **⚠ The old trailer is gone the moment this starts writing, and a crash in the middle leaves
/// a file every reader refuses.** That is the right outcome — a half-written footer must not
/// look finished — but it means this is not a safe in-place edit of a file whose *existing*
/// trailer matters. Spec §6.4 gives append the same warning and the same answer: **write to a new
/// path and rename if that matters to the caller.**
///
/// **What it checks, and what it deliberately does not.** It reads the fixed tail — which checks
/// its own internal agreements, the index ending exactly where the trailer begins among them —
/// and adds the two rules that need the file's length: that the sections end exactly where the
/// footer begins, and that the index does not start inside the header. **It does not read the
/// header or the block index**, and spec §6.7's table is why: the only refusals it lists for this
/// operation are a missing footer and the file's own bytes. Reading them would make the cheap
/// operation cost an index decode per call, and it could not make the file safer than it already
/// is — **every field this writes back except the trailer's length is the file's own**.
pub fn replace_trailer(path: &Path, trailer: &[u8]) -> Result<(), PspWriteError> {
    let reopen = |source: PspReadError| PspWriteError::Reopen {
        path: path.to_path_buf(),
        source,
    };
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .open(path)
        .map_err(|source| PspWriteError::Io {
            path: path.to_path_buf(),
            while_doing: "reopening the file to replace its trailer",
            source,
        })?;
    let file_bytes = file
        .metadata()
        .map_err(|source| PspWriteError::Io {
            path: path.to_path_buf(),
            while_doing: "measuring the file",
            source,
        })?
        .len();

    let footer = read_the_footer(path, &mut file, file_bytes).map_err(reopen)?;

    // **From the trailer's offset forward, and nothing before it.** The blocks and the index end
    // where the trailer begins, so this seek is the whole of what makes the operation cheap.
    file.seek(SeekFrom::Start(footer.trailer_offset))
        .map_err(|source| PspWriteError::Io {
            path: path.to_path_buf(),
            while_doing: "seeking to the trailer",
            source,
        })?;

    let replaced = Footer {
        trailer_bytes: trailer.len() as u64,
        ..footer
    };
    // **Read back before believing it**, the obligation `finish` carries for the same reason: a
    // footer that reaches disk and is then rejected has cost a file that was readable before.
    decode_footer(&encode_footer(&replaced)).map_err(|source| {
        PspWriteError::WouldNotBeReadable {
            path: path.to_path_buf(),
            reason: "the footer it would write does not decode".to_string(),
            source: Some(source.into()),
        }
    })?;

    file.write_all(trailer)
        .map_err(|source| PspWriteError::Io {
            path: path.to_path_buf(),
            while_doing: "writing the trailer",
            source,
        })?;
    file.write_all(&encode_footer(&replaced))
        .map_err(|source| PspWriteError::Io {
            path: path.to_path_buf(),
            while_doing: "writing the footer",
            source,
        })?;

    // **The file has to end here, and a shorter trailer is why.** Without this, replacing a
    // twenty-byte trailer with a five-byte one leaves fifteen bytes of the old one past the new
    // footer — and a reader takes the *last* forty-eight bytes for the footer, so it would read
    // the tail of a trailer as one and refuse the file.
    let ends_at = footer.trailer_offset + trailer.len() as u64 + FOOTER_BYTES as u64;
    file.set_len(ends_at).map_err(|source| PspWriteError::Io {
        path: path.to_path_buf(),
        while_doing: "trimming the file to its new length",
        source,
    })?;

    // **Durable, for the reason `finish` is** (spec §6.3): a caller that is told the trailer was
    // replaced and then loses power must not find the old one.
    file.sync_all().map_err(|source| PspWriteError::Io {
        path: path.to_path_buf(),
        while_doing: "syncing the file",
        source,
    })?;
    Ok(())
}

/// Read the fixed tail and check the two rules that need the file's length.
///
/// **The same rules `PspReader::open` applies**, and for the same reason: 48 bytes cannot check
/// their own offsets against a file whose length they do not know. They are repeated rather than
/// shared because `open` reads the index and the header in the same breath, and this operation
/// deliberately reads neither.
fn read_the_footer(
    path: &Path,
    file: &mut std::fs::File,
    file_bytes: u64,
) -> Result<Footer, PspReadError> {
    let Some(footer_at) = file_bytes.checked_sub(FOOTER_BYTES as u64) else {
        return Err(PspReadError::Incomplete {
            path: path.to_path_buf(),
        });
    };
    file.seek(SeekFrom::Start(footer_at))
        .map_err(|source| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing: "seeking to the footer",
            source,
        })?;
    let mut tail = [0u8; FOOTER_BYTES];
    file.read_exact(&mut tail)
        .map_err(|source| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing: "reading the footer",
            source,
        })?;

    let footer = decode_footer(&tail).map_err(|refused| match refused {
        // **No tail magic means the writer never finished**, which is the everyday case and the
        // one spec §6.7 names for this operation. Telling a foreign file apart from a killed one
        // is `open`'s job; here both answers are the same — there is no trailer to replace.
        super::footer::FooterDecodeError::NotAFooter { .. } => PspReadError::Incomplete {
            path: path.to_path_buf(),
        },
        damaged => PspReadError::damaged_by(path, "the footer does not decode", damaged.into()),
    })?;

    let sections_end = footer
        .trailer_offset
        .checked_add(footer.trailer_bytes)
        .ok_or_else(|| {
            PspReadError::damaged(path, "the trailer's end is past any address".to_string())
        })?;
    if sections_end != footer_at {
        return Err(PspReadError::damaged(
            path,
            format!(
                "the footer says the file's sections end at byte {sections_end}, but the footer \
                 itself begins at byte {footer_at}"
            ),
        ));
    }
    if footer.index_offset < HEAD_MAGIC.len() as u64 {
        return Err(PspReadError::damaged(
            path,
            format!(
                "the footer puts the block index at byte {}, which is inside the header",
                footer.index_offset
            ),
        ));
    }
    Ok(footer)
}

#[cfg(test)]
mod tests {
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
        match refused {
            PspWriteError::Reopen {
                source: PspReadError::Incomplete { .. },
                ..
            } => {}
            other => panic!("got {other}"),
        }
    }

    /// **A footer that does not describe the file it sits in is refused**, rather than being
    /// rewritten with the same wrong offsets and a fresh trailer.
    #[test]
    fn a_footer_that_does_not_describe_the_file_is_refused() {
        let (_dir, path) = a_finished_psp();
        let mut whole = bytes_of(&path);
        let footer_at = whole.len() - FOOTER_BYTES;
        let mut footer = footer_of(&whole);
        footer.trailer_bytes += 1;
        whole[footer_at..].copy_from_slice(&encode_footer(&footer));
        rewrite(&path, &whole);

        let refused = replace_trailer(&path, b"anything").expect_err("the sections do not add up");
        match refused {
            PspWriteError::Reopen {
                source: PspReadError::Damaged { .. },
                ..
            } => {}
            other => panic!("got {other}"),
        }
    }

    /// **A file that was never a psp is refused, and nothing is written to it.** The bytes are
    /// compared afterwards: a refusal that had already written the trailer would leave a file
    /// that is neither what it was nor a psp.
    #[test]
    fn a_file_that_is_not_a_psp_is_refused_and_left_alone() {
        let (_dir, path) = a_file();
        let not_a_psp = b"this is a text file, and it is not a psp at all".to_vec();
        rewrite(&path, &not_a_psp);

        let refused = replace_trailer(&path, b"anything").expect_err("that is not a psp");
        assert!(
            matches!(refused, PspWriteError::Reopen { .. }),
            "got {refused}"
        );
        assert_eq!(bytes_of(&path), not_a_psp, "the file is exactly as it was");
    }

    /// **A file that is not there names itself and what was being done to it**, rather than
    /// arriving as a bare `No such file or directory`.
    #[test]
    fn a_missing_file_names_itself() {
        let (_dir, path) = a_file();
        let refused = replace_trailer(&path, b"anything").expect_err("there is no file");
        match refused {
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
        match refused {
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
}
