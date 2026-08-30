//! **A record as bytes** — the fixed columns, and the padding rule that makes an empty allele
//! writable (`doc/devel/ng/spec/vcf_output.md` §5).
//!
//! **Encoding is total: it returns text, not a `Result`.** Everything that could make a record
//! unwritable was refused when the record was built — the reference spells its span, the
//! per-allele vectors match the table, a padding base is present exactly when one is needed and
//! on the side the span allows. So there is nothing left for this stage to reject, and the
//! signature says so. What *is* fallible is the I/O that follows, and keeping the two apart is
//! what lets spec §11's rule — "formatting must be a pure function of the record" — be read off
//! the type rather than promised in a comment.

use std::fmt::Write as _;

use super::{PaddingBase, VcfRecord};
use crate::ng::vcf::header::HeaderContig;

/// What a VCF writes where it has no value: the `ID` column of every record ng emits, and the
/// `ALT` column of a record that established no alternative.
pub const MISSING_FIELD: &str = ".";

/// How many decimal places a quality is written to.
///
/// **One, matching both production writers.** The precision is part of the format rather than a
/// display choice: spec §11 requires the same run to produce the same bytes, so a quality
/// rendered two ways is two different files.
pub const QUALITY_DECIMALS: usize = 1;

/// **The seven fixed columns**, `CHROM` through `FILTER`, tab-separated and with no trailing
/// separator — `INFO`, `FORMAT` and the sample columns follow in later steps.
///
/// `contigs` is the header's contig list, which a record's contig id indexes.
///
/// # Panics
///
/// On a contig id the header does not hold, and on allele bytes that are not UTF-8. Both are
/// wiring defects rather than inputs: the record and the header are built from one run's
/// reference, and allele bases are `A`/`C`/`G`/`T`/`N`-validated at the merge boundary. A
/// writer that emitted `?` for an unknown contig would produce a file that parses and names the
/// wrong chromosome, which is the failure this refuses to have.
#[must_use]
pub fn fixed_columns(record: &VcfRecord, contigs: &[HeaderContig]) -> String {
    let contig_index = record.region().contig.0 as usize;
    let contig = contigs.get(contig_index).unwrap_or_else(|| {
        panic!(
            "this record is on contig {contig_index} and the header holds {}: the record and \
             the header were built from different references, and a CHROM written from the \
             wrong table names the wrong chromosome rather than failing",
            contigs.len()
        )
    });

    let mut out = String::new();
    let _ = write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{:.*}\t{}",
        contig.name,
        written_position(record),
        MISSING_FIELD,
        padded(record.reference(), record.padding_base()),
        alternatives_column(record),
        QUALITY_DECIMALS,
        record.site_quality.get(),
        record.filter().as_str(),
    );
    out
}

/// **The `POS` column**, after the padding rule has been applied.
///
/// A left-hand padding base is the reference base *before* the span, so it becomes the record's
/// first base and the position moves back one. A right-hand one is appended instead and the
/// position does not move — which is only ever the contig's first base, where there is nothing
/// to the left to move onto.
#[must_use]
fn written_position(record: &VcfRecord) -> u64 {
    let start = record.region().start.get();
    match record.padding_base() {
        // The record type refuses a left-hand base at position 1, so this cannot underflow.
        Some(PaddingBase::Left(_)) => start - 1,
        Some(PaddingBase::Right(_)) | None => start,
    }
}

/// The `ALT` column: every alternative, padded, comma-joined — or `.` where the locus
/// established none.
fn alternatives_column(record: &VcfRecord) -> String {
    if record.alternatives().is_empty() {
        return MISSING_FIELD.to_string();
    }
    let padding = record.padding_base();
    let mut column = String::new();
    for (index, allele) in record.alternatives().iter().enumerate() {
        if index > 0 {
            column.push(',');
        }
        column.push_str(&padded(allele, padding));
    }
    column
}

/// One allele's bases with the padding base applied, or unchanged where the record carries
/// none.
///
/// **This is what makes an empty allele writable**, and it is applied to *every* allele of the
/// record rather than to the empty one alone: VCF states a deletion by giving all alleles a
/// shared flanking base, so padding one and not the others would describe a different variant.
fn padded(allele: &[u8], padding: Option<PaddingBase>) -> String {
    let bases = std::str::from_utf8(allele)
        .expect("allele bases are A/C/G/T/N-validated at the merge boundary");
    match padding {
        None => bases.to_string(),
        Some(PaddingBase::Left(base)) => format!("{}{}", base as char, bases),
        Some(PaddingBase::Right(base)) => format!("{}{}", bases, base as char),
    }
}

#[cfg(test)]
mod tests;
