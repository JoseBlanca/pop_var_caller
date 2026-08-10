//! Deriving the genome's segmentation from the catalog, with no FASTA open.
//!
//! **Every gate below is step 3's own**, reached through the file instead of through a
//! scan: the same pre-screen, the same bundling, the same admission, the same feature
//! resolution. That is not tidiness — it is the only way the differential of
//! `catalog_differential.rs` can come out *identical* rather than merely close (spec §5.1).
//!
//! What the file replaces is exactly the three values that needed the bases — the
//! whole-motif cut, the motif and the purity — which the builder computed once and every
//! reader now inherits.

use crate::ng::region_typing::segment_criteria::{
    self, SsrSegment, SsrSegmentCriteria, is_compound, split_bundles,
};
use crate::ng::region_typing::{
    TypedRegion, TypedRegionConfig, coverage_runs, fill_generic_gaps, resolve_features,
};
use crate::ng::repeat_catalog::criteria::StrRepeatCriteria;
use crate::ng::repeat_catalog::{FoundRepeat, RepeatCatalogError, TractSpan};
use crate::ng::tandem_repeat::RepeatInterval;
use crate::ng::types::{Bp, ContigId, Motif};

/// One contig's rows, turned into the typed regions a live scan would have produced.
///
/// `rows` are the catalog's rows for `contig`, in file order; `contig_len` is that contig's
/// own length, which the flank test asks about and no row carries.
///
/// The output covers `1..=contig_len` with no gap and no overlap, exactly as
/// [`crate::ng::region_typing::partition_resident`] does — the generic stretches are the
/// spaces between the repeat features (spec §5.1).
pub fn segments_of_contig(
    chrom: &str,
    contig: ContigId,
    contig_len: Bp,
    rows: &[FoundRepeat],
    criteria: &StrRepeatCriteria,
) -> Vec<TypedRegion> {
    if contig_len.get() == 0 {
        return Vec::new();
    }
    let class = &criteria.classification;

    // 1. The detected spans, back in the detector's own coordinates, because that is what
    //    the pre-screen and the bundling read (spec §3.1).
    let detected: Vec<RepeatInterval> = rows.iter().map(detected_interval).collect();

    // 2. The pre-screen, unchanged: step 3's own `prefilter` over those intervals.
    let cleaned = segment_criteria::prefilter(&detected, class);

    // 3. Admission, gate for gate as `classify` runs it — but reading the stored cut,
    //    motif and purity instead of the bases (§5.1).
    let admitted = admit(chrom, rows, &cleaned, contig_len, class);

    // 4-5. The satellite cap over the cleaned coverage, the features, then the gaps.
    let runs = coverage_runs(&cleaned);
    let config = TypedRegionConfig {
        criteria: class.clone(),
        max_str_len: criteria.max_str_len_bp,
        ..TypedRegionConfig::default()
    };
    let features = resolve_features(&runs, admitted.loci, &admitted.bundled, contig, &config);
    fill_generic_gaps(features, contig, contig_len.get())
}

/// Every segment of every contig the catalog holds, in reference order.
///
/// `contigs` are the header's, so their order is the file's and their lengths are the
/// reference's.
pub fn segments_of_rows<I>(
    rows_by_contig: I,
    criteria: &StrRepeatCriteria,
) -> Result<Vec<TypedRegion>, RepeatCatalogError>
where
    I: IntoIterator<Item = (String, ContigId, Bp, Vec<FoundRepeat>)>,
{
    let mut out = Vec::new();
    for (chrom, contig, contig_len, rows) in rows_by_contig {
        out.extend(segments_of_contig(
            &chrom, contig, contig_len, &rows, criteria,
        ));
    }
    Ok(out)
}

/// What admission produced: the loci, and the bundle members it set aside.
struct Admitted {
    loci: Vec<SsrSegment>,
    bundled: Vec<RepeatInterval>,
}

/// `classify`'s gates, in `classify`'s order, over stored fields.
///
/// The one difference is where three values come from: the whole-motif cut, the motif and
/// the purity are read off the row rather than recomputed from bases. Everything the
/// criteria decide is decided here, now, so a reader may move any of them.
fn admit(
    chrom: &str,
    rows: &[FoundRepeat],
    cleaned: &[RepeatInterval],
    contig_len: Bp,
    class: &SsrSegmentCriteria,
) -> Admitted {
    // The pre-screen returns intervals, not rows; match them back by span and period so the
    // stored cut travels with each survivor.
    let mut survivors: Vec<&FoundRepeat> = Vec::with_capacity(cleaned.len());
    for interval in cleaned {
        if let Some(row) = rows.iter().find(|row| &detected_interval(row) == interval) {
            survivors.push(row);
        }
    }

    // Scope, score and the compound-motif gate — `classify`'s step 2, over stored fields.
    let mut kept: Vec<&FoundRepeat> = Vec::with_capacity(survivors.len());
    for row in survivors {
        if row.period < class.periods.min()
            || row.period > class.periods.max()
            || row.score < class.min_score
            || is_compound(row.motif.as_bytes())
        {
            continue;
        }
        kept.push(row);
    }

    // Bundling, on the raw pre-trim coordinates and start-sorted, as `classify` does.
    kept.sort_by_key(|row| (row.detected.start, row.detected.end));
    let intervals: Vec<RepeatInterval> = kept.iter().map(|row| detected_interval(row)).collect();
    let (kept_intervals, bundled) = split_bundles(intervals, class.bundle_threshold);

    let mut loci = Vec::with_capacity(kept_intervals.len());
    for interval in &kept_intervals {
        let Some(row) = kept.iter().find(|row| detected_interval(row) == *interval) else {
            continue;
        };
        if let Some(locus) = finish_from_row(chrom, row, contig_len, class) {
            loci.push(locus);
        }
    }

    Admitted { loci, bundled }
}

/// `finish_locus` with the trim, the motif and the purity already known.
///
/// The gates and their order are that function's: the cut must exist, the **trimmed** tract
/// must clear the copy floor, the purity must clear its floor, and the flank must not clamp
/// to nothing at either contig end.
fn finish_from_row(
    chrom: &str,
    row: &FoundRepeat,
    contig_len: Bp,
    class: &SsrSegmentCriteria,
) -> Option<SsrSegment> {
    let trimmed = row.trimmed?;
    if trimmed.repeat_count(row.period) < u64::from(class.min_copies.for_period(row.period)) {
        return None;
    }
    let purity = row.purity?;
    if purity < class.min_purity {
        return None;
    }

    // The flank test, and it is the contig's ends that answer it — never a window's.
    //
    // `.max(1)` mirrors `finish_locus` rather than reasoning about when it can fire: for a
    // catalog row it cannot, since the file's own 15 bp flank floor already excludes every
    // tract nearer the start than the default bundle radius. Mutating it away leaves the
    // differential green, which is the honest reason it is here — parity with the function
    // this one stands in for, not a case of its own.
    let tract_start = trimmed.start.get();
    let tract_end = trimmed.end.get();
    let ref_start = tract_start.saturating_sub(class.bundle_threshold).max(1);
    let ref_end = (tract_end + class.bundle_threshold).min(contig_len.get());
    if ref_start == tract_start || ref_end == tract_end {
        return None;
    }

    SsrSegment::new(
        chrom.to_string().into_boxed_str(),
        tract_start,
        tract_end,
        Motif::new(row.motif.as_bytes()).ok()?,
        purity,
    )
    .ok()
}

/// A row's detected span, back in the detector's 0-based half-open coordinates.
///
/// The inverse of the conversion in [`crate::ng::repeat_catalog::row`], and the reason the
/// round trip is lossless: the file stores what the detector said, not a reshaping of it.
fn detected_interval(row: &FoundRepeat) -> RepeatInterval {
    RepeatInterval {
        start: row.detected.start.get() - 1,
        end: row.detected.end.get(),
        period: row.period,
        score: row.score,
    }
}

/// The trimmed span of a row, for callers that want the locus rather than the detection.
#[inline]
pub fn locus_span(row: &FoundRepeat) -> Option<TractSpan> {
    row.trimmed
}
