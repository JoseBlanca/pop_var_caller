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

use super::{PaddingBase, SampleCall, SampleColumn, VcfRecord};
use crate::ng::calling::quality::MAX_GENOTYPE_QUALITY;
use crate::ng::types::{Phred, Ploidy};
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

/// How many decimal places an allele frequency is written to.
///
/// **Six, because the smallest frequency the caller can state is smaller than four would
/// show.** One copy in a cohort of 3,000 diploid samples is 1 in 6,000, or 0.000167; at four
/// decimals two different singleton frequencies render identically, and at three they render as
/// zero.
pub const FREQUENCY_DECIMALS: usize = 6;

/// How many decimal places a pooled mapping quality is written to.
///
/// Two. Mapping qualities are integers 0..=60 and these are means over reads, so the fraction
/// carries the pooling and nothing finer is meaningful.
pub const MAPPING_QUALITY_DECIMALS: usize = 2;

/// How many decimal places an artifact penalty is written to.
///
/// One, the same as [`QUALITY_DECIMALS`] — the penalties are Phreds subtracted from the site
/// quality, so writing them at a different precision would stop `QUAL + ABPEN + SPPEN` from
/// recovering the uncorrected quality the way spec §6 says it does.
pub const PENALTY_DECIMALS: usize = 1;

/// **A site or genotype quality, as the file writes it** — [`QUALITY_DECIMALS`] places.
///
/// **Negative zero cannot reach here**, and that is upstream's doing rather than this
/// function's: [`Phred::try_new`] normalises the sign at its only door, because `from_log_prob`
/// produces exactly `-0.0` at a log-probability of zero and a `QUAL` column cannot carry a minus
/// sign on a certainty. Production's repeat-tract writer adds `0.0` at write time for the same
/// reason; ng does not need to, and the test table pins that it does not need to.
#[must_use]
pub fn quality_text(quality: Phred) -> String {
    format!("{:.*}", QUALITY_DECIMALS, quality.get())
}

/// An allele frequency, as the file writes it — [`FREQUENCY_DECIMALS`] places.
///
/// **A frequency below the precision renders as zero**, which is a real limit rather than a
/// defect: at six places the smallest distinguishable value is one in a million, and a cohort
/// would need 500,000 diploid samples for a singleton to fall under it.
#[must_use]
pub fn frequency_text(frequency: f64) -> String {
    format!("{frequency:.*}", FREQUENCY_DECIMALS)
}

/// A pooled mapping quality or a difference of two, as the file writes it —
/// [`MAPPING_QUALITY_DECIMALS`] places.
///
/// Signed: `MQDIFF` is negative exactly when an alternative's reads map worse than the
/// reference's, which is the multi-mapper signal the field exists to publish.
#[must_use]
pub fn mapping_quality_text(mapping_quality: f64) -> String {
    format!("{mapping_quality:.*}", MAPPING_QUALITY_DECIMALS)
}

/// An artifact penalty, as the file writes it — [`PENALTY_DECIMALS`] places, the same as a
/// quality's, so that `QUAL` plus the two penalties recovers the uncorrected quality in the
/// file's own digits.
#[must_use]
pub fn penalty_text(penalty: Phred) -> String {
    format!("{:.*}", PENALTY_DECIMALS, penalty.get())
}

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
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        contig.name,
        written_position(record),
        MISSING_FIELD,
        padded(record.reference(), record.padding_base()),
        alternatives_column(record),
        quality_text(record.site_quality),
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

/// **The `INFO` column**: the site's own annotations, `;`-separated, in the order spec §6
/// declares them.
///
/// **A key whose value is undefined for this record is omitted; a defined key whose *entry* is
/// undefined writes `.` in that slot.** The two mean different things to a parser and production
/// keeps them apart, so this does too: `MQREF` disappears when no read reached the reference,
/// while `MQALT` stays and writes `.` for the one alternative nobody's reads reached.
///
/// The per-alternative keys — `AF`, `AC`, `MQALT`, `MQDIFF`, all `Number=A` — are omitted
/// entirely at a record with no alternative, which is what a locus the caller refused comes to.
#[must_use]
pub fn info_column(record: &VcfRecord) -> String {
    let mut fields: Vec<String> = Vec::new();
    let counts = called_allele_counts(record);

    if !record.alternatives().is_empty() {
        fields.push(format!("AF={}", allele_frequencies(record)));
        fields.push(format!(
            "AC={}",
            join_with_commas(counts.per_alternative.iter().map(u64::to_string))
        ));
    }
    fields.push(format!("AN={}", counts.total));
    fields.push(format!("DP={}", total_depth(record)));

    if let Some(penalties) = record.artifact_penalties {
        fields.push(format!("ABPEN={}", penalty_text(penalties.allele_balance)));
        fields.push(format!(
            "SPPEN={}",
            penalty_text(penalties.strand_and_read_position)
        ));
    }

    let reference_mapq = record.allele_mapq()[0].mean();
    if let Some(mean) = reference_mapq {
        fields.push(format!("MQREF={}", mapping_quality_text(mean)));
    }
    if !record.alternatives().is_empty() {
        let alternative_means: Vec<Option<f64>> = record.allele_mapq()[1..]
            .iter()
            .map(|pool| pool.mean())
            .collect();
        fields.push(format!(
            "MQALT={}",
            join_with_commas(alternative_means.iter().map(optional_mapping_quality))
        ));
        fields.push(format!(
            "MQDIFF={}",
            join_with_commas(alternative_means.iter().map(|alternative| {
                optional_mapping_quality(
                    &alternative
                        .zip(reference_mapq)
                        .map(|(alternative, reference)| alternative - reference),
                )
            }))
        ));
    }

    if let Some(tract) = record.repeat_tract() {
        let motif = std::str::from_utf8(tract.motif())
            .expect("motif bases are A/C/G/T/N-validated at the catalog boundary");
        fields.push("STR".to_string());
        fields.push(format!("RU={motif}"));
        fields.push(format!("PERIOD={}", tract.period()));
    }

    fields.join(";")
}

/// The `AF` values: each alternative's share of the cohort's expected allele copies.
///
/// **Normalised over the copies' own total, not over `AN`, and the difference is real.** The
/// copies are the calling loop's converged fit and they sum to `ploidy ×` the samples the *loop*
/// scored; `AN` counts the samples the *file* writes a genotype for, which is fewer whenever a
/// sample the loop scored is written `./.` — the ordinary case for a sample whose reads said
/// nothing (spec §7.1). Dividing by `AN` would therefore make the frequencies sum to more than
/// one. `AF` is an estimate and `AC`/`AN` are counts of called genotypes; they are different
/// quantities and are allowed different denominators.
fn allele_frequencies(record: &VcfRecord) -> String {
    let total: f64 = record.expected_copies().iter().sum();
    join_with_commas(record.expected_copies()[1..].iter().map(|copies| {
        // A locus every sample was set aside at cannot reach here — the record type refuses an
        // empty cohort — but a cohort whose copies are all zero is arithmetic, not data.
        let frequency = if total > 0.0 { copies / total } else { 0.0 };
        frequency_text(frequency)
    }))
}

/// `AC` per alternative and `AN`, counted from the genotypes the file actually writes.
struct CalledAlleleCounts {
    per_alternative: Vec<u64>,
    total: u64,
}

/// Count the called genotypes' alleles.
///
/// **A no-called sample is in neither number.** `AN` is "total called allele copies", so a
/// sample the file writes as `./.` contributes nothing to it — which is what makes `AN` differ
/// from `ploidy × samples` and is the whole reason it is written at all.
fn called_allele_counts(record: &VcfRecord) -> CalledAlleleCounts {
    let mut per_alternative = vec![0u64; record.alternatives().len()];
    let mut total = 0u64;
    for column in record.sample_columns() {
        let Some(genotype) = column.call.genotype() else {
            continue;
        };
        for allele in genotype.alleles() {
            total += 1;
            let index = usize::from(allele.get());
            if index > 0 {
                // PANIC-FREE: `VcfRecord::new` refused a genotype naming an allele past the
                // table, so `index - 1` is in range for the alternatives.
                per_alternative[index - 1] += 1;
            }
        }
    }
    CalledAlleleCounts {
        per_alternative,
        total,
    }
}

/// The record's `INFO/DP`: the sum of the samples' own depths.
fn total_depth(record: &VcfRecord) -> u64 {
    record
        .sample_columns()
        .iter()
        .map(|column| u64::from(column.read_counts.depth()))
        .sum()
}

/// One `Number=A` entry that may be undefined — a mean over no reads, or a difference against
/// one.
fn optional_mapping_quality(value: &Option<f64>) -> String {
    value.map_or_else(|| MISSING_FIELD.to_string(), mapping_quality_text)
}

fn join_with_commas(values: impl Iterator<Item = String>) -> String {
    values.collect::<Vec<_>>().join(",")
}

/// **The `FORMAT` column**: which per-sample fields this record carries, in the order the
/// sample columns write them.
///
/// A repeat-tract record carries one more than a SNP or indel — `REPCN`, each called allele's
/// length in whole repeat units — and that is the only difference between the two.
#[must_use]
pub fn format_keys(record: &VcfRecord) -> &'static str {
    if record.is_repeat_tract() {
        "GT:GQ:DP:AD:REPCN"
    } else {
        "GT:GQ:DP:AD"
    }
}

/// **The sample columns**, tab-separated, one per sample of the run in the run's sample order.
///
/// `ploidy` is the run's, and it is a parameter rather than a field on the record because a
/// no-call has to be spelled `./.` at a ploidy the sample itself no longer states: a called
/// sample's ploidy is the length of its genotype, but a **refused locus writes every sample as
/// a no-call** (spec §8), so there is no sibling to read it from. One ploidy per run is what
/// the run's frozen parameters hold.
#[must_use]
pub fn sample_columns(record: &VcfRecord, ploidy: Ploidy) -> String {
    let mut out = String::new();
    for (index, column) in record.sample_columns().iter().enumerate() {
        if index > 0 {
            out.push('\t');
        }
        out.push_str(&one_sample_column(record, column, ploidy));
    }
    out
}

/// One sample's fields, `:`-separated, matching [`format_keys`].
fn one_sample_column(record: &VcfRecord, column: &SampleColumn, ploidy: Ploidy) -> String {
    let genotype = genotype_field(&column.call, ploidy);
    let quality = column.call.genotype_quality().map_or_else(
        || MISSING_FIELD.to_string(),
        |phred| written_genotype_quality(phred).to_string(),
    );
    let depth = column.read_counts.depth().to_string();
    let allele_depths =
        join_with_commas(column.read_counts.allele_reads().iter().map(u32::to_string));

    let mut fields = format!("{genotype}:{quality}:{depth}:{allele_depths}");
    if let Some(tract) = record.repeat_tract() {
        let repeat_counts = column.call.genotype().map_or_else(
            || MISSING_FIELD.to_string(),
            |genotype| {
                join_with_commas(genotype.alleles().iter().map(|allele| {
                    // PANIC-FREE: `VcfRecord::new` refused a genotype naming an allele past the
                    // table, so every id indexes the allele table.
                    let bases = &record.alleles()[usize::from(allele.get())];
                    tract.repeat_copies_of(bases).to_string()
                }))
            },
        );
        fields.push(':');
        fields.push_str(&repeat_counts);
    }
    fields
}

/// The `GT` field: the called alleles, `/`-joined and never phased — or the no-call, which is
/// one `.` per copy of the genome.
///
/// **The alleles are already ascending**, because [`Genotype`](crate::ng::types::Genotype)
/// sorts at construction and has no other way in. That is what makes `REPCN` line up with `GT`
/// for free: both walk the same slice. Production's repeat-tract writer sorts its `GT` while
/// building `REPCN` from the unsorted candidate order, so the two fields' entries need not
/// correspond — a mismatch that is invisible unless the two alleles have different repeat
/// counts.
///
/// **Unphased always.** ng computes no phasing, so a `|` would claim something no step
/// established.
fn genotype_field(call: &SampleCall, ploidy: Ploidy) -> String {
    match call.genotype() {
        Some(genotype) => join_with_separator(
            genotype
                .alleles()
                .iter()
                .map(|allele| allele.get().to_string()),
            '/',
        ),
        None => join_with_separator(
            std::iter::repeat_n(MISSING_FIELD.to_string(), usize::from(ploidy.get())),
            '/',
        ),
    }
}

/// The `GQ` field's integer value.
///
/// Rounded to the nearest whole Phred and held to [`MAX_GENOTYPE_QUALITY`], the cap GATK and
/// bcftools use and the one this project's quality module already applies — restated here
/// because the file's own field is what downstream tools read, and production clamps at write
/// time for the same reason.
fn written_genotype_quality(quality: Phred) -> u32 {
    let capped = quality.get().clamp(0.0, MAX_GENOTYPE_QUALITY);
    capped.round() as u32
}

fn join_with_separator(values: impl Iterator<Item = String>, separator: char) -> String {
    values.collect::<Vec<_>>().join(&separator.to_string())
}

#[cfg(test)]
mod tests;
