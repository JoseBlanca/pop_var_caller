//! **Production's catalog post-filter, copied** — `src/ssr/catalog/postprocess.rs` at commit
//! `d9e7b076`. Promotion step C12 moved it into `scanner_parity.rs` when production's deletion would
//! otherwise have taken the golden path's second half with it, and step C24 moved it here so that
//! `region_typing::segment_criteria`'s port-fidelity differential can compare against it too.
//!
//! The code below is production's, line for line, with three changes and no others: it takes the
//! golden catalog's settings ([`GoldenCatalogSettings`]) where production took its `CatalogParams`
//! (the same four fields); production's `Motif::new` length check is written inline; and it returns
//! a [`PostFilteredLocus`] — where the tract is, what repeats in it and its recomputed purity —
//! where production built a `Locus`, whose own checks cannot fire on what this function hands it.
//! **It is frozen**: `scanner_parity` measures ng's scanner against trf-mod's path, and
//! `segment_criteria` measures ng's own port of the same policy
//! ([`classify`](crate::ng::region_typing::segment_criteria::classify)) against it, so it must not
//! follow either.

use crate::ng::golden_catalog::GoldenCatalogSettings;

/// Production's `TrfRecord`, as `TrfRecord::for_test` built it: `frac_match` 1 and an empty
/// pattern, neither of which the post-filter reads.
#[derive(Debug, Clone)]
pub(crate) struct TrfRecord {
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) period: u16,
    pub(crate) score: i32,
}

/// What the post-filter keeps of a tract: where it is, what repeats in it and how purely.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PostFilteredLocus {
    /// The contig's name.
    pub(crate) chrom: String,
    /// The trimmed tract's first base, 0-based.
    pub(crate) start: u32,
    /// One past the trimmed tract's last base.
    pub(crate) end: u32,
    /// The repeat unit, upper-cased, in the tract's own phase.
    pub(crate) motif: Vec<u8>,
    /// The fraction of the trimmed tract matching a perfect tiling of the motif.
    pub(crate) purity_fraction: f32,
}

/// Production's `ssr::types::MAX_MOTIF_LEN`.
const MAX_MOTIF_LEN: usize = 6;

/// Per-period minimum copy number a tract must reach to survive (GangSTR
/// `minimal_trim.py` `thresholds = {1:10, 2:5, 3:4, 4:3, 5:3, 6:3}`; default 3
/// for any other period). Period 1 is filtered out upstream ([`MIN_PERIOD`]),
/// so its floor is retained only for parity with the source table.
const fn copy_number_floor(period: usize) -> u32 {
    match period {
        1 => 10,
        2 => 5,
        3 => 4,
        4 => 3,
        5 => 3,
        6 => 3,
        _ => 3,
    }
}

/// The narrowest SSR period the catalog keeps. Period-1 **homopolymers are
/// excluded** (architecture §4): the standard GangSTR/HipSTR drop — error-prone
/// for STR genotyping and not the di/tri/tetra-nucleotide target — and dropping
/// them *before* bundling stops a long poly-A/T run from bundle-dropping an
/// adjacent real SSR.
const MIN_PERIOD: u16 = 2;

/// The widest SSR period the catalog keeps (architecture §4 step 1).
const MAX_PERIOD: u16 = 6;

/// Post-process one contig's TRF-mod records into start-sorted catalog loci.
/// `chrom` is the contig name; `contig_seq` is the full contig (any case — the
/// tract, motif, and `ref_seq` are upper-cased here for case-stable identity).
pub(crate) fn build_loci(
    recs: Vec<TrfRecord>,
    chrom: &str,
    contig_seq: &[u8],
    p: &GoldenCatalogSettings,
) -> Vec<PostFilteredLocus> {
    // 1. scope + score gate, then 2. compound-motif drop. Both need the tract,
    //    so we slice the (upper-cased) prefix motif here and filter on it.
    let mut kept: Vec<TrfRecord> = recs
        .into_iter()
        .filter(|r| {
            r.period >= MIN_PERIOD
                && r.period <= MAX_PERIOD
                && r.score >= p.min_score
                && r.end > r.start
                && (r.end as usize) <= contig_seq.len()
        })
        .filter(|r| {
            let period = r.period as usize;
            let tract = &contig_seq[r.start as usize..r.end as usize];
            // Need at least one full motif to form the prefix.
            if tract.len() < period {
                return false;
            }
            let motif = upper(&tract[..period]);
            !is_compound(&motif)
        })
        .collect();

    // 3. drop bundles (on the raw, pre-trim coordinates). Records must be
    //    start-sorted for the streaming clustering.
    kept.sort_by_key(|r| (r.start, r.end));
    let kept = drop_bundles(kept, p.bundle_threshold);

    // 4-5. per record: end-trim + copy floor, recompute purity + floor, embed.
    let mut out = Vec::with_capacity(kept.len());
    for r in &kept {
        if let Some(locus) = finish_locus(r, chrom, contig_seq, p) {
            out.push(locus);
        }
    }
    out
}

/// Steps 4-5 for one record: end-trim, copy-number floor, motif, purity floor,
/// `ref_seq` embed. `None` if the record fails any gate.
fn finish_locus(
    r: &TrfRecord,
    chrom: &str,
    contig_seq: &[u8],
    p: &GoldenCatalogSettings,
) -> Option<PostFilteredLocus> {
    let period = r.period as usize;
    let raw_tract = upper(&contig_seq[r.start as usize..r.end as usize]);
    let motif_bytes = raw_tract.get(..period)?.to_vec();

    // End-trim to clean whole-motif boundaries (GangSTR minimal_trim).
    let (st, en) = minimal_trim(&raw_tract, &motif_bytes)?;
    let new_start = r.start + st as u32;
    let new_end = r.start + en as u32;
    let trimmed = &raw_tract[st..en];

    // Copy-number floor — GangSTR computes copies from the ORIGINAL TRF span
    // (integer division), as an accept-gate after trimming.
    let ref_copy = (r.end - r.start) / r.period as u32;
    if ref_copy < copy_number_floor(period) {
        return None;
    }

    // After minimal_trim the tract starts on a motif boundary, so
    // `trimmed[..period] == motif_bytes`; use it as the phase-faithful motif.
    // Production's `Motif::new`, which refused an empty unit or one longer than six bases.
    if motif_bytes.is_empty() || motif_bytes.len() > MAX_MOTIF_LEN {
        return None;
    }

    // Recompute purity from the trimmed tract vs a perfect motif tiling.
    let purity = recompute_purity(trimmed, &motif_bytes);
    if purity < p.min_purity {
        return None;
    }

    // Embed ref_seq: trimmed tract + flank each side, clamped at contig ends.
    let ref_start = new_start.saturating_sub(p.flank_bp);
    let ref_end = (new_end + p.flank_bp).min(contig_seq.len() as u32);

    // Drop a locus whose flank clamped to nothing on either side — a tract
    // abutting position 0 of the contig (empty left flank) or ending on the
    // contig's last base (empty right flank). The empirical-candidate delimiter
    // (`ssr::pileup::alignment`) anchors the repeat region on *both* flank
    // junctions; a zero-length flank leaves nothing to anchor against, so the
    // tract is not genotypeable and must not reach Stage 1.
    if ref_start == new_start || ref_end == new_end {
        return None;
    }

    // Production's `Locus::new` also embedded the flanked reference bytes and refused
    // inconsistent coordinates or a purity outside [0, 1]; by construction here neither refusal
    // can fire, and the tests read only where the tract is, what repeats in it and its purity.
    Some(PostFilteredLocus {
        chrom: chrom.to_string(),
        start: new_start,
        end: new_end,
        motif: motif_bytes,
        purity_fraction: purity,
    })
}

/// ASCII upper-case copy of `bytes`.
fn upper(bytes: &[u8]) -> Vec<u8> {
    bytes.iter().map(|b| b.to_ascii_uppercase()).collect()
}

/// Count non-overlapping greedy occurrences of `motif` in `repeat` (GangSTR
/// `count_motif`): scan left-to-right, advancing by `motif.len()` on a match
/// and by 1 otherwise.
fn count_motif(repeat: &[u8], motif: &[u8]) -> usize {
    let m = motif.len();
    let mut c = 0;
    let mut s = 0;
    while s < repeat.len() {
        if s + m <= repeat.len() && &repeat[s..s + m] == motif {
            c += 1;
            s += m;
        } else {
            s += 1;
        }
    }
    c
}

/// `true` if `motif` is itself internally periodic — a non-fundamental period
/// such as `ATAT = (AT)²` (GangSTR `is_compound`, threshold 0.8). A motif whose
/// shorter prefix `sub` tiles more than 80% of it is compound.
fn is_compound(motif: &[u8]) -> bool {
    let l = motif.len();
    const THRESHOLD: f64 = 0.8;
    for i in 1..=(l / 2) {
        let sub = &motif[..i];
        if (count_motif(motif, sub) * i) as f64 > l as f64 * THRESHOLD {
            return true;
        }
    }
    false
}

/// End-trim a tract to clean whole-motif boundaries (GangSTR `minimal_trim`):
/// find the smallest `start_offset` where two motif copies (`motif*2`) appear,
/// and the largest `end_offset` (exclusive) where they appear before it.
/// Returns `(start_offset, end_offset)` into `rep`, or `None` if the tract has
/// no clean motif boundary within `min(3·|motif|, |rep|/2)` of each end.
fn minimal_trim(rep: &[u8], motif: &[u8]) -> Option<(usize, usize)> {
    let m = motif.len();
    let ll = m * 2;
    let max_trim_len = (m * 3).min(rep.len() / 2);

    let two_motifs = |slice: &[u8]| -> bool {
        slice.len() == ll && &slice[..m] == motif && &slice[m..] == motif
    };

    // Smallest start offset whose next 2·m bytes are motif*2. GangSTR checks
    // the bound first and bails (not continues) — replicate.
    let mut start = None;
    for so in 0..=max_trim_len {
        if so + ll >= rep.len() {
            return None;
        }
        if two_motifs(&rep[so..so + ll]) {
            start = Some(so);
            break;
        }
    }
    let st = start?;

    // Largest end offset (exclusive) whose preceding 2·m bytes are motif*2.
    let mut end = None;
    for eo in (rep.len().saturating_sub(max_trim_len)..=rep.len()).rev() {
        if eo < ll {
            return None;
        }
        if two_motifs(&rep[eo - ll..eo]) {
            end = Some(eo);
            break;
        }
    }
    let en = end?;
    if st >= en {
        return None;
    }
    Some((st, en))
}

/// Fraction of `tract` matching a perfect tiling of `motif` from phase 0 (spec
/// §3.2). With `motif == tract[..period]`, the first repeat always matches, so
/// a perfect tract scores 1.0 and interruptions lower it proportionally.
fn recompute_purity(tract: &[u8], motif: &[u8]) -> f32 {
    if tract.is_empty() || motif.is_empty() {
        return 0.0;
    }
    let m = motif.len();
    let matches = tract
        .iter()
        .enumerate()
        .filter(|&(i, &b)| b == motif[i % m])
        .count();
    matches as f32 / tract.len() as f32
}

/// Two records are "close" (bundle candidates) if any of their start/end
/// coordinates are within `thresh` bp (GangSTR `is_close`, `check_motif=False`;
/// chrom equality is implicit — these are one contig's records).
fn is_close(a: &TrfRecord, b: &TrfRecord, thresh: u32) -> bool {
    a.start.abs_diff(b.start) < thresh
        || a.start.abs_diff(b.end) < thresh
        || b.start.abs_diff(a.end) < thresh
        || b.end.abs_diff(a.end) < thresh
}

/// Drop bundles: discard every maximal run of consecutive close records in
/// their entirety (GangSTR `remove_bundles.py`). `recs` must be start-sorted.
fn drop_bundles(recs: Vec<TrfRecord>, thresh: u32) -> Vec<TrfRecord> {
    let mut out = Vec::new();
    let n = recs.len();
    let mut i = 0;
    while i < n {
        if i + 1 < n && is_close(&recs[i], &recs[i + 1], thresh) {
            // Advance through the whole close cluster; drop all of it.
            let mut j = i;
            while j + 1 < n && is_close(&recs[j], &recs[j + 1], thresh) {
                j += 1;
            }
            i = j + 1;
        } else {
            out.push(recs[i].clone());
            i += 1;
        }
    }
    out
}
