//! **Rewriting a finished record's `FILTER` and `INFO`, and touching nothing else.**
//!
//! Pass three writes back the lines pass one parked, and the verdict changes at most two of
//! their columns: `FILTER` gains the filter's id where a record is tagged, and `INFO` gains the
//! ratio and the posterior where a record was scored. Every other column — the sample columns
//! above all, which are all but nine of a cohort's line — goes out as it came in.
//!
//! **That is what makes byte-identity checkable rather than hoped for.** The standing oracle is
//! that a run with the filter on at a target no record reaches, with the two `INFO` keys
//! stripped, equals the run with the filter off exactly (spec §10). A patch that rebuilt the
//! line from its parts would have to reproduce every column's spelling; this one keeps the
//! bytes and splices two of them.
//!
//! **The two rules that are easy to get backwards** (spec §6 traps 5 and 7):
//!
//! - **A `FILTER` that already says something is joined, not replaced.** A record the calling
//!   loop marked `EMNoConv` and the filter then flags carries `EMNoConv;hiddenParalog`, which
//!   the VCF grammar allows and which the header will declare once step C4 adds the
//!   `hiddenParalog` id to [`FILTER_DECLARATIONS`](crate::ng::vcf::header). Only `PASS` — and
//!   `.`, which means no filter was applied — are replaced, because joining to either would
//!   state a contradiction.
//! - **An `INFO` that says nothing is replaced, not appended to.** `.` is VCF's *missing*, and
//!   `.;PARALOG_LR=1.5` is not a field list with a missing value in it, it is a parse error
//!   waiting to happen.
//!
//! **Both of those spellings are defensive, not reachable.** ng's encoder writes neither in
//! either column: [`FilterVerdict`](crate::ng::vcf::FilterVerdict) has five values and all are
//! words, and `info_column` pushes `AN=` and `DP=` before any conditional field, so the thinnest
//! `INFO` any record gets is `AN=0;DP=0`. The rules are kept because a patch that got them
//! backwards would be silently wrong on a line from anywhere else, and because they cost two
//! comparisons a record.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.5.

use thiserror::Error;

/// How many pieces the line is split into: the six columns before `FILTER`, then `FILTER` and
/// `INFO` — the two the verdict may touch — then everything from `FORMAT` onwards in one piece.
///
/// **`FORMAT` and the sample columns are never separated**, so a cohort of three thousand costs
/// one split of nine rather than three thousand of one.
const PIECES: usize = 9;

/// VCF's *missing* value, which both columns may hold and neither may be joined to.
const MISSING: &[u8] = crate::ng::vcf::MISSING_FIELD.as_bytes();

/// The `FILTER` value that means every filter passed, which is replaced rather than joined.
const PASS: &[u8] = crate::ng::vcf::FilterVerdict::Pass.as_str().as_bytes();

/// Rewrite a record's line, adding a filter id, some `INFO` fields, or neither.
///
/// With nothing to add the line comes back byte for byte — and it comes back through the same
/// splice as every other call, so that identity is a property of the code rather than of a
/// short circuit around it.
///
/// **Call once per line.** A second call adds the id and the fields again, giving
/// `hiddenParalog;hiddenParalog` and a repeated `PARALOG_LR`, which strict readers reject.
///
/// **What this does not refuse**, because a line reaching pass three came from the encoder by way
/// of the spill and neither adds these: a `\n` anywhere in `line`, which makes the writer put two
/// record lines in the file while counting one; and a `\t` in `filter_to_add` or in an
/// `info_to_add` field, which adds a column and shifts every sample column right. The column
/// count below is the only shape this checks.
///
/// # Errors
///
/// If the line does not split into at least nine tab-separated pieces — the eight fixed columns,
/// then `FORMAT` and the sample columns together.
pub fn rewrite_filter_and_info(
    line: &[u8],
    filter_to_add: Option<&str>,
    info_to_add: &[String],
) -> Result<Vec<u8>, LinePatchError> {
    // A fixed-size array rather than a slice pattern, so the destructure's arity is checked
    // against `PIECES` at compile time: changing `PIECES` without changing the columns named
    // below is `error[E0527]`, not a run-time refusal of every real line.
    let pieces: Vec<&[u8]> = line.splitn(PIECES, |&byte| byte == b'\t').collect();
    let pieces: [&[u8]; PIECES] = match pieces.try_into() {
        Ok(pieces) => pieces,
        Err(short) => {
            return Err(LinePatchError::TooFewColumns {
                columns: Vec::<&[u8]>::len(&short),
            });
        }
    };
    let [
        chrom,
        position,
        id,
        reference,
        alternatives,
        quality,
        filter,
        info,
        format_and_samples,
    ] = pieces;

    let mut patched = Vec::with_capacity(line.len() + added_byte_count(filter_to_add, info_to_add));
    for column in [chrom, position, id, reference, alternatives, quality] {
        patched.extend_from_slice(column);
        patched.push(b'\t');
    }
    append_filter_column(&mut patched, filter, filter_to_add);
    patched.push(b'\t');
    append_info_column(&mut patched, info, info_to_add);
    patched.push(b'\t');
    patched.extend_from_slice(format_and_samples);
    Ok(patched)
}

/// The `FILTER` column: the id joined on, unless what is there is `PASS`, missing, or empty.
fn append_filter_column(out: &mut Vec<u8>, filter: &[u8], filter_to_add: Option<&str>) {
    let Some(added) = filter_to_add else {
        out.extend_from_slice(filter);
        return;
    };
    if filter != PASS && filter != MISSING && !filter.is_empty() {
        out.extend_from_slice(filter);
        out.push(b';');
    }
    out.extend_from_slice(added.as_bytes());
}

/// The `INFO` column: the fields appended, unless what is there says nothing.
fn append_info_column(out: &mut Vec<u8>, info: &[u8], info_to_add: &[String]) {
    if info_to_add.is_empty() {
        out.extend_from_slice(info);
        return;
    }
    let says_something = info != MISSING && !info.is_empty();
    if says_something {
        out.extend_from_slice(info);
    }
    let mut needs_separator = says_something;
    for field in info_to_add {
        if needs_separator {
            out.push(b';');
        }
        out.extend_from_slice(field.as_bytes());
        needs_separator = true;
    }
}

/// How much longer the patched line will be, so the buffer is sized once.
///
/// **This must be an upper bound on what [`append_filter_column`] and [`append_info_column`]
/// write**, or the `Vec` reallocates and copies a line that runs to tens of kilobytes at three
/// thousand samples. Nothing but this sentence links the three, so
/// `the_capacity_is_never_an_under_estimate` in the tests checks it over every combination of
/// `FILTER` and `INFO` spelling.
fn added_byte_count(filter_to_add: Option<&str>, info_to_add: &[String]) -> usize {
    let filter = filter_to_add.map_or(0, |added| added.len() + 1);
    let info: usize = info_to_add.iter().map(|field| field.len() + 1).sum();
    filter + info
}

/// What can go wrong patching a line.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum LinePatchError {
    /// The line does not split into the nine pieces the patch needs.
    ///
    /// **The patch needs nine**: the eight fixed columns, then `FORMAT` and the sample columns
    /// together. A line that splits into fewer did not come from the encoder, so rewriting its
    /// seventh piece would be rewriting something other than `FILTER`.
    ///
    /// **The tenth piece is not checked, and a record's line always has one.**
    /// [`VcfRecord::new`](crate::ng::vcf::VcfRecord) refuses a record with no sample columns, so
    /// the encoder cannot produce a line with `FORMAT` and nothing after it. This guard is about
    /// the pieces the patch must index, not about whether the line is a valid cohort record.
    ///
    /// **The error names no record.** This function is handed bytes and nothing else. A caller
    /// walking a spill — step C4 — should wrap it with the entry's contig and position, or a
    /// failure on a file of millions of records is untraceable.
    #[error(
        "a record's line splits into {columns} tab-separated piece(s), and rewriting FILTER and \
         INFO needs the eight fixed columns followed by FORMAT"
    )]
    #[non_exhaustive]
    TooFewColumns {
        /// How many pieces the line split into — at most eight, since nine is what is needed.
        columns: usize,
    },
}

#[cfg(test)]
mod tests;
