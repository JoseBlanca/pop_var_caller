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
    self, RejectionReason, SsrSegment, SsrSegmentCriteria, is_compound, split_bundles,
};
use crate::ng::region_typing::{
    RegionKind, TypedRegion, TypedRegionConfig, coverage_runs, fill_generic_gaps, resolve_features,
};
use crate::ng::repeat_catalog::FoundRepeat;
use crate::ng::repeat_catalog::criteria::StrRepeatCriteria;
use crate::ng::repeat_catalog::tally::{CatalogRejectionCounts, ContigTally};
use crate::ng::tandem_repeat::{RepeatInterval, ScanParams};
use crate::ng::types::{Bp, ContigId, GenomeRegion, Motif, Position};

/// One contig's rows, turned into the typed regions a live scan would have produced.
///
/// `rows` are the catalog's rows for `contig`, in file order; `contig_len` is that contig's
/// own length, which the flank test asks about and no row carries.
///
/// The output covers `1..=contig_len` with no gap and no overlap, exactly as
/// [`crate::ng::region_typing::partition_resident`] does — the generic stretches are the
/// spaces between the repeat features (spec §5.1).
///
/// The [`ContigTally`] beside the regions is what cannot be read off them: this contig's
/// repeat coverage, and what admission turned down here.
pub fn segments_of_contig(
    chrom: &str,
    contig: ContigId,
    contig_len: Bp,
    rows: &[FoundRepeat],
    criteria: &StrRepeatCriteria,
) -> (Vec<TypedRegion>, ContigTally) {
    if contig_len.get() == 0 {
        return (Vec::new(), ContigTally::default());
    }
    let (features, tally) = repeat_features_of_contig(chrom, contig, contig_len, rows, criteria);
    (fill_generic_gaps(features, contig, contig_len.get()), tally)
}

/// The segments of one contig that overlap `wanted`, clipped to it.
///
/// **The whole contig is typed first, and only then clipped** — which is what step 3's walk
/// does (`region_typing/mod.rs`, `SpanWalk::requested`), and it has to be: a tract outside
/// the requested span still bundles with one inside it, and a satellite outside still
/// swallows a locus inside. Classifying only the requested stretch would answer a different
/// question and quietly answer it differently.
///
/// **Only a generic stretch is clipped.** A locus, a bundle and a satellite each come out
/// whole, because each is a claim about its own extent and half of it is a different claim:
/// half a locus is not a locus, a clipped bundle's members describe bases outside their own
/// region, and a satellite clipped to 100 bases contradicts the very cap that made it one.
/// A generic stretch is the only kind that is not a finding — *nothing more specific can be
/// said here* stays true of any part of it. That is the walk's rule, stated in
/// `region_typing/mod.rs`'s `clips_at_a_bed_edge`, and this is where the two must agree.
///
/// **The tally is the whole contig's, not the requested stretches'**, because the whole
/// contig is what was typed — the walk's own tally works the same way, since its scan span
/// is always the whole contig whatever was asked for (`region_typing/mod.rs`, `scan_set`).
pub fn segments_of_contig_in(
    chrom: &str,
    contig: ContigId,
    contig_len: Bp,
    rows: &[FoundRepeat],
    criteria: &StrRepeatCriteria,
    wanted: &[GenomeRegion],
) -> (Vec<TypedRegion>, ContigTally) {
    let (all, tally) = segments_of_contig(chrom, contig, contig_len, rows, criteria);
    if wanted.is_empty() {
        return (Vec::new(), tally);
    }

    let mut out = Vec::new();
    for region in all {
        for span in wanted.iter().filter(|s| s.contig == contig) {
            let start = region.region.start.get().max(span.start.get());
            let end = region.region.end.get().min(span.end.get());
            if start > end {
                continue;
            }
            let whole = matches!(
                region.kind,
                RegionKind::SsrSegment(_) | RegionKind::SsrBundle { .. } | RegionKind::Satellite
            );
            out.push(TypedRegion {
                region: if whole {
                    region.region
                } else {
                    GenomeRegion {
                        contig,
                        start: Position(start),
                        end: Position(end),
                    }
                },
                kind: region.kind.clone(),
            });
            if whole {
                // Emitted whole once, however many requested spans it touches.
                break;
            }
        }
    }
    (out, tally)
}

/// The **STR loci** of one contig, without the generic stretches between them.
///
/// The same admission [`segments_of_contig`] runs — a locus here is a locus there,
/// satellite absorption included — but the generic spans a caller would immediately discard
/// are never built (spec §5.3).
pub fn loci_of_contig(
    chrom: &str,
    contig: ContigId,
    contig_len: Bp,
    rows: &[FoundRepeat],
    criteria: &StrRepeatCriteria,
) -> Vec<SsrSegment> {
    if contig_len.get() == 0 {
        return Vec::new();
    }
    repeat_features_of_contig(chrom, contig, contig_len, rows, criteria)
        .0
        .into_iter()
        .filter_map(|region| match region.kind {
            RegionKind::SsrSegment(locus) => Some(locus),
            _ => None,
        })
        .collect()
}

/// The STR loci of one contig that overlap `wanted`.
///
/// A locus is inside or it is not — there is no clipping, since half a locus is not one — so
/// this is a filter, where [`segments_of_contig_in`] is a filter plus a cut for the stretches
/// between loci.
pub fn loci_of_contig_in(
    chrom: &str,
    contig: ContigId,
    contig_len: Bp,
    rows: &[FoundRepeat],
    criteria: &StrRepeatCriteria,
    wanted: &[GenomeRegion],
) -> Vec<SsrSegment> {
    loci_of_contig(chrom, contig, contig_len, rows, criteria)
        .into_iter()
        .filter(|locus| {
            wanted.iter().any(|s| {
                s.contig == contig && locus.start() <= s.end.get() && locus.end() >= s.start.get()
            })
        })
        .collect()
}

/// Everything up to but not including the generic fill: the repeat features of one contig,
/// coordinate-ordered.
///
/// **Absorption is why the loci cannot be taken before this point**: a satellite swallows a
/// locus too close to it, so a locus is only a locus once the features are resolved.
fn repeat_features_of_contig(
    chrom: &str,
    contig: ContigId,
    contig_len: Bp,
    rows: &[FoundRepeat],
    criteria: &StrRepeatCriteria,
) -> (Vec<TypedRegion>, ContigTally) {
    let class = &criteria.classification;

    // 1. The detected spans, back in the detector's own coordinates, because that is what
    //    the pre-screen and the bundling read (spec §3.1).
    let detected: Vec<RepeatInterval> = rows.iter().map(detected_interval).collect();

    // 2. The pre-screen, unchanged: step 3's own `prefilter` over those intervals.
    let cleaned = segment_criteria::prefilter(&detected, class);

    // 3. Admission, gate for gate as `classify` runs it — but reading the stored cut,
    //    motif and purity instead of the bases (§5.1).
    let admitted = admit(chrom, rows, &cleaned, contig_len, class);

    // 4-5. The satellite cap over the cleaned coverage, then the features.
    let runs = coverage_runs(&cleaned);
    // The coverage the walk charges to `repeat_bp_with_no_locus`: the merged cleaned
    // intervals, so overlapping repeats count their bases once. The walk sums the same runs
    // window by window, clipped to each core; cores tile a contig, so the totals are the
    // same number reached two ways.
    let repeat_bp = runs.iter().map(|run| run.len()).sum();
    // **Every field named, no `..default()` tail.** The walk's own construction carries the
    // same rule for the same reason: a field added to `TypedRegionConfig` later must break
    // this line rather than silently take the calling walk's value, and the differential
    // cannot catch that because it would default on both sides at once. `scan` and
    // `window_bp` are the two the catalog genuinely has no use for — there is no scan here
    // and there are no windows — and each says so.
    let config = TypedRegionConfig {
        criteria: class.clone(),
        max_str_len: criteria.max_str_len_bp,
        // Nothing is scanned: the tracts came from the file. `resolve_features` reads the
        // weights not at all, so this value cannot reach an answer.
        scan: ScanParams::default(),
        // A whole contig is resolved at once, so there is no window to size.
        window_bp: Bp(crate::ng::region_typing::DEFAULT_WINDOW_BP),
    };
    (
        resolve_features(&runs, admitted.loci, &admitted.bundled, contig, &config),
        ContigTally {
            repeat_bp,
            rejected: admitted.rejected,
        },
    )
}

/// What admission produced: the loci, the bundle members it set aside, and the repeats it
/// turned down.
struct Admitted {
    loci: Vec<SsrSegment>,
    bundled: Vec<RepeatInterval>,
    rejected: CatalogRejectionCounts,
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
    // The pre-screen and the bundler both hand back **intervals**, not rows, and what a
    // locus needs — the whole-motif cut, the motif, the purity — lives on the row. So each
    // survivor has to be matched back to the row it came from.
    //
    // **Both sequences are sorted the same way and the survivors are a subsequence**, since
    // `prefilter` only drops, so one cursor walks them together. A lookup per survivor
    // would scan a chromosome's rows once per survivor.
    let survivors: Vec<&FoundRepeat> = pair_with_rows(rows, cleaned);

    // Scope, score and the compound-motif gate — `classify`'s step 2, over stored fields.
    //
    // **Only the compound gate is charged**, exactly as `classify` charges it: a repeat
    // outside the period range or under the score floor is out of the question being asked,
    // not turned down by it, so neither is a rejection in either implementation.
    let mut rejected = CatalogRejectionCounts::default();
    let mut kept: Vec<&FoundRepeat> = Vec::with_capacity(survivors.len());
    for row in survivors {
        if row.period < class.periods.min()
            || row.period > class.periods.max()
            || row.score < class.min_score
        {
            continue;
        }
        if is_compound(row.motif.as_bytes()) {
            rejected.add(RejectionReason::Compound, row.detected.len_bp());
            continue;
        }
        kept.push(row);
    }

    // Bundling, on the raw pre-trim coordinates and start-sorted, as `classify` does.
    kept.sort_by_key(|row| (row.detected.start, row.detected.end));
    let intervals: Vec<RepeatInterval> = kept.iter().map(|row| detected_interval(row)).collect();
    let (kept_intervals, bundled) = split_bundles(intervals, class.bundle_threshold);

    // The isolated intervals are a subsequence of what went in, and `kept` is in that same
    // order, so the same cursor walk pairs them — the rows this time being `kept`, not the
    // contig's, since the bundling was fed `kept`'s order.
    let isolated_rows = pair_with_rows_by_ref(&kept, &kept_intervals);

    let mut loci = Vec::with_capacity(isolated_rows.len());
    for row in isolated_rows {
        match finish_from_row(chrom, row, contig_len, class) {
            Ok(locus) => loci.push(locus),
            // The bases charged are the **detected** length, before any trim — what the walk
            // charges (`RejectionCounts`), and what makes the two tallies comparable.
            Err(reason) => rejected.add(reason, row.detected.len_bp()),
        }
    }

    Admitted {
        loci,
        bundled,
        rejected,
    }
}

/// Pair each interval of `subsequence` with the row it came from, walking both once.
///
/// `subsequence` must be what it says: the same intervals `rows` describes, in the same
/// order, with some dropped. An interval that does not appear is skipped rather than
/// searched for — it cannot have come from these rows.
fn pair_with_rows<'a>(
    rows: &'a [FoundRepeat],
    subsequence: &[RepeatInterval],
) -> Vec<&'a FoundRepeat> {
    let mut out = Vec::with_capacity(subsequence.len());
    let mut cursor = rows.iter();
    for interval in subsequence {
        for row in cursor.by_ref() {
            if detected_interval(row) == *interval {
                out.push(row);
                break;
            }
        }
    }
    out
}

/// [`pair_with_rows`] over rows already borrowed — the second walk, whose left-hand side is
/// the survivors of the first rather than the contig's rows.
fn pair_with_rows_by_ref<'a>(
    rows: &[&'a FoundRepeat],
    subsequence: &[RepeatInterval],
) -> Vec<&'a FoundRepeat> {
    let mut out = Vec::with_capacity(subsequence.len());
    let mut cursor = rows.iter();
    for interval in subsequence {
        for row in cursor.by_ref() {
            if detected_interval(row) == *interval {
                out.push(*row);
                break;
            }
        }
    }
    out
}

/// `finish_locus` with the trim, the motif and the purity already known.
///
/// The gates and their order are that function's: the cut must exist, the **trimmed** tract
/// must clear the copy floor, the purity must clear its floor, and the flank must not clamp
/// to nothing at either contig end.
///
/// **The refusal carries its reason**, and it is the walk's own
/// [`RejectionReason`] rather than a second vocabulary — the two tallies are compared
/// against each other, so one name per gate is what makes that comparison mean anything.
fn finish_from_row(
    chrom: &str,
    row: &FoundRepeat,
    contig_len: Bp,
    class: &SsrSegmentCriteria,
) -> Result<SsrSegment, RejectionReason> {
    let trimmed = row.trimmed.ok_or(RejectionReason::NoCleanTrim)?;
    if trimmed.repeat_count(row.period) < u64::from(class.min_copies.for_period(row.period)) {
        return Err(RejectionReason::CopyFloor);
    }
    // `purity` is `Some` exactly when `trimmed` is, so this cannot be the arm that fires
    // once the cut is known to exist; it is here because the type says it can be, and
    // `NoCleanTrim` is what a missing purity would mean.
    let purity = row.purity.ok_or(RejectionReason::NoCleanTrim)?;
    if purity < class.min_purity {
        return Err(RejectionReason::Purity);
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
        return Err(RejectionReason::FlankClamped);
    }

    // Unreachable, both of them: the motif was checked when the row was written, and every
    // gate above guarantees `SsrSegment`'s invariants. They are charged to `NoCleanTrim` for
    // the reason `finish_locus` gives — of the reasons available it is the one that means
    // "this tract did not turn out to be a well-formed repeat" — and the `debug_assert` is
    // the part that matters, since a count here would be a count of arithmetic bugs.
    let motif = Motif::new(row.motif.as_bytes()).map_err(|e| {
        debug_assert!(false, "a stored motif is not a valid motif: {e}");
        RejectionReason::NoCleanTrim
    })?;
    SsrSegment::new(
        chrom.to_string().into_boxed_str(),
        tract_start,
        tract_end,
        motif,
        purity,
    )
    .map_err(|e| {
        debug_assert!(
            false,
            "the catalog built an invalid locus, so the arithmetic is wrong: {e} \
             (tract [{tract_start}, {tract_end}], contig is {} long)",
            contig_len.get()
        );
        RejectionReason::NoCleanTrim
    })
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::repeat_catalog::TractSpan;
    use crate::ng::repeat_catalog::tally::CatalogRejectionCounts;
    use crate::ng::types::Position;

    const CONTIG_LEN: Bp = Bp(10_000);

    fn span(start: u64, end: u64) -> TractSpan {
        TractSpan {
            start: Position(start),
            end: Position(end),
        }
    }

    /// A row that clears the pre-screen: ten whole copies, far from either contig end and
    /// far enough from its neighbours not to bundle.
    fn row(start: u64, period: u8, motif: &str) -> FoundRepeat {
        let end = start + u64::from(period) * 10 - 1;
        FoundRepeat {
            contig: ContigId(0),
            detected: span(start, end),
            trimmed: Some(span(start, end)),
            period,
            score: 100,
            motif: Motif::new(motif.as_bytes()).expect("a valid motif"),
            purity: Some(1.0),
        }
    }

    fn tally_of(rows: &[FoundRepeat], criteria: &StrRepeatCriteria) -> (usize, ContigTally) {
        let (features, tally) =
            repeat_features_of_contig("chr1", ContigId(0), CONTIG_LEN, rows, criteria);
        let loci = features
            .iter()
            .filter(|r| matches!(r.kind, RegionKind::SsrSegment(_)))
            .count();
        (loci, tally)
    }

    /// **Each gate charges its own counter, and charges the detected length.** Two of these
    /// four a live scan never reaches — the detector emits only primitive motifs, and it
    /// emits the purest sub-segment of an impure tract — so a fixture reference cannot drive
    /// them and the rows are handed to admission directly, the same way
    /// `segment_criteria`'s own compound-motif test does.
    #[test]
    fn every_gate_charges_its_own_counter() {
        let criteria = StrRepeatCriteria::default();

        // No whole-motif boundaries to cut back to.
        let mut no_trim = row(1_001, 3, "CAG");
        no_trim.trimmed = None;
        no_trim.purity = None;

        // Ten detected copies, three surviving the cut — under the period-3 floor of 4.
        let mut short_after_trim = row(2_001, 3, "CAG");
        short_after_trim.trimmed = Some(span(2_001, 2_009));

        // Half the bases match a perfect tiling, against a floor of 0.8.
        let mut impure = row(3_001, 3, "CAG");
        impure.purity = Some(0.5);

        // `ATAT` is `AT` twice, so the period is a lie.
        let compound = row(4_001, 4, "ATAT");

        let rows = vec![no_trim, short_after_trim, impure, compound];
        let (loci, tally) = tally_of(&rows, &criteria);

        assert_eq!(loci, 0, "every row here is turned down");
        assert_eq!(
            tally.rejected,
            CatalogRejectionCounts {
                no_clean_trim: 30,
                copy_floor: 30,
                purity: 30,
                compound: 40,
            },
            "the bases charged are each row's detected length, before any cut"
        );
    }

    /// The mirror case: a row that clears every gate is a locus and charges nothing, so the
    /// counters above are not simply charging everything they see.
    #[test]
    fn a_row_that_clears_every_gate_charges_nothing() {
        let criteria = StrRepeatCriteria::default();
        let (loci, tally) = tally_of(&[row(1_001, 3, "CAG")], &criteria);

        assert_eq!(loci, 1);
        assert_eq!(tally.rejected, CatalogRejectionCounts::default());
        assert_eq!(
            tally.repeat_bp, 30,
            "the contig's repeat coverage is the one tract"
        );
    }

    /// Repeat coverage is bases, not tracts: two tracts sharing a base cover it once.
    #[test]
    fn overlapping_rows_charge_a_shared_base_once() {
        let criteria = StrRepeatCriteria::default();
        let mut second = row(1_030, 3, "CAG");
        second.trimmed = Some(span(1_030, 1_059));

        // 1001..=1030 and 1030..=1059: sixty bases of tract over fifty-nine of contig.
        let (_, tally) = tally_of(&[row(1_001, 3, "CAG"), second], &criteria);
        assert_eq!(tally.repeat_bp, 59);
    }

    /// A repeat outside the reader's period range or under its score floor is **not** a
    /// rejection — it is out of the question being asked, exactly as `classify` treats it.
    #[test]
    fn out_of_scope_is_not_a_rejection() {
        let base = StrRepeatCriteria::default();
        let narrow = StrRepeatCriteria {
            classification: SsrSegmentCriteria {
                periods: crate::ng::tandem_repeat::PeriodRange::new(4, 6).expect("valid"),
                ..base.classification.clone()
            },
            ..base.clone()
        };
        let (loci, tally) = tally_of(&[row(1_001, 3, "CAG")], &narrow);

        assert_eq!(loci, 0, "period 3 is outside the range asked for");
        assert_eq!(tally.rejected, CatalogRejectionCounts::default());
    }
}
