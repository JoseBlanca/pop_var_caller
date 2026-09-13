//! **Golden-catalog parity for ng's tandem-repeat scanner.**
//!
//! [`find_tandem_repeats`](crate::ng::tandem_repeat::find_tandem_repeats) must
//! reproduce the **final catalog** the external `trf-mod` produces, run through
//! production's post-filter (`ssr::catalog::postprocess`) — since promotion step C12 a frozen
//! copy of it, [`production_post_filter`], and the golden catalog read by ng's own
//! [`golden_catalog`](crate::ng::golden_catalog).
//!
//! The oracle is a golden catalog snapshot committed under
//! `tests/data/tandem_repeat/` — built once by the current `trf-mod` →
//! `postprocess` path (`ssr-catalog`) on a synthetic STR-diversity reference. This
//! test bypasses `trf-mod`: it runs the scanner on the same reference, bridges the
//! intervals into the post-filter's input record, and asserts the resulting locus set
//! reproduces the golden one. This only *validates* ng's scanner.
//!
//! ## Why this lives in `src/ng/`, having been written in `src/ssr/catalog/`
//!
//! It is **ng's test**: its subject is ng's scanner, and production is only the
//! yardstick. It began in `src/ssr/catalog/` because that is where the golden path
//! was, and that made `src/ssr/` — which is **frozen** (owner, 2026-07-16: *"leave
//! production as is"*) — the **only** part of the tree that `use`d `crate::ng`.
//!
//! That dependency pointed the wrong way, and it bit: typed-regions Milestone B2
//! widens ng's coordinates to `u64`, and `RepeatInterval`'s `u32`-ness was baked
//! into this file (four sites, most bindingly `TrfRecord::for_test(start: u32,
//! end: u32, …)`) — so an **ng-internal** change could not compile without editing
//! frozen production. A freeze that ng can break is not a freeze; the whole point
//! is that an experiment cannot destabilise the yardstick it is scored against.
//!
//! Moving it (owner-approved) severs that: **production now depends on nothing in
//! ng**, permanently, and this test follows ng's types as they move. The move is
//! test-only — no shipped behaviour changed, and the parity bar is the same one.
//!
//! ## The catalog pre-filter (a validated integration finding)
//!
//! `trf-mod` hands `postprocess` a **clean** candidate set — statistically
//! significant repeats, redundancy already eliminated. The raw scanner is
//! deliberately permissive (`min_copies = 2`), so in aperiodic sequence it emits
//! many low-copy noise repeats and every period-multiple of a real tract. Fed
//! straight to `build_loci`, that noise triggers `drop_bundles` (which runs
//! *before* the copy-number floor) and cascades away the real loci. So a consumer
//! must apply, **before** `build_loci`, the same two cleanups `trf-mod` bakes in —
//! the per-period copy floor and period-multiple redundancy elimination (its
//! `IsRedundant`).
//!
//! [`catalog_prefilter`] is that policy, kept here **verbatim as trf-mod's
//! shape**. Note ng has since ported it properly, with the floors as a real knob,
//! to [`crate::ng::region_typing::segment_criteria::prefilter`] — this copy stays frozen
//! against the golden path rather than tracking ng's, so that the two oracles stay
//! independent: if ng's prefilter drifts, `classification`'s own differential says so,
//! and this one keeps measuring the scanner against trf-mod.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use crate::ng::golden_catalog::{
    GoldenCatalogSettings, GoldenLocus, golden_catalog, tandem_repeat_fixture as fixture,
};
use crate::ng::tandem_repeat::{PeriodRange, RepeatInterval, ScanParams, find_tandem_repeats};
use production_post_filter::{PostFilteredLocus, TrfRecord, build_loci};

/// **Production's catalog post-filter, copied** — `src/ssr/catalog/postprocess.rs` at commit
/// `d9e7b076`, moved here at promotion step C12 when production's deletion would otherwise have
/// taken the golden path's second half with it.
///
/// The code below is production's, line for line, with three changes and no others: it takes the
/// golden catalog's settings ([`GoldenCatalogSettings`]) where production took its `CatalogParams`
/// (the same four fields); production's `Motif::new` length check is written inline; and it returns
/// a [`PostFilteredLocus`] where production built a `Locus`, whose own checks cannot fire on what
/// this function hands it. **It is frozen**: this file measures ng's scanner against trf-mod's
/// path, and it stays independent of ng's own port of the same policy
/// ([`classify`](crate::ng::region_typing::segment_criteria::classify)) for the reason the module header gives.
mod production_post_filter {
    use crate::ng::golden_catalog::GoldenCatalogSettings;

    /// Production's `TrfRecord`, as `TrfRecord::for_test` built it: `frac_match` 1 and an empty
    /// pattern, neither of which the post-filter reads.
    #[derive(Debug, Clone)]
    pub(super) struct TrfRecord {
        pub(super) start: u32,
        pub(super) end: u32,
        pub(super) period: u16,
        pub(super) score: i32,
    }

    /// What the post-filter keeps of a tract: where it is and what repeats in it.
    #[derive(Debug, Clone, PartialEq)]
    pub(super) struct PostFilteredLocus {
        pub(super) chrom: String,
        pub(super) start: u32,
        pub(super) end: u32,
        pub(super) motif: Vec<u8>,
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
    pub(super) fn build_loci(
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
        // can fire, and the test reads only where the tract is and what repeats in it.
        Some(PostFilteredLocus {
            chrom: chrom.to_string(),
            start: new_start,
            end: new_end,
            motif: motif_bytes,
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
}

/// The per-period minimum copy count `postprocess` applies (as its `copy_number_floor`),
/// used here as a **pre-`build_loci`** filter — the catalog's job, since `drop_bundles` runs
/// before the post-filter's own floor.
fn copy_floor(period: u8) -> u32 {
    match period {
        2 => 5,
        3 => 4,
        _ => 3,
    }
}

/// The catalog's pre-`build_loci` cleanup of raw scanner intervals (see the module docs):
/// keep period-2..6 intervals clearing the per-period copy floor, then drop period-multiple
/// redundancy (a higher-period interval overlapping a divisor-period one is the same tract).
fn catalog_prefilter(intervals: &[RepeatInterval]) -> Vec<RepeatInterval> {
    let mut floored: Vec<RepeatInterval> = intervals
        .iter()
        .copied()
        .filter(|iv| {
            iv.period >= 2
                && (iv.end - iv.start) / u64::from(iv.period) >= u64::from(copy_floor(iv.period))
        })
        .collect();
    // Process low periods first so a fundamental tract is kept and its multiples dropped.
    floored.sort_by_key(|iv| (iv.period, iv.start));
    let mut kept: Vec<RepeatInterval> = Vec::new();
    for iv in floored {
        let redundant = kept.iter().any(|a| {
            a.period < iv.period
                && iv.period % a.period == 0
                && iv.start < a.end
                && a.start < iv.end
        });
        if !redundant {
            kept.push(iv);
        }
    }
    kept
}

/// Run the scanner path (`find_tandem_repeats` → catalog pre-filter → `build_loci`) over the
/// reference the golden was built from, and return its loci.
fn scanner_loci(reference: &Path, params: &GoldenCatalogSettings) -> Vec<PostFilteredLocus> {
    let file = File::open(reference).unwrap();
    let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));
    let mut out = Vec::new();
    for result in reader.records() {
        let rec = result.unwrap();
        let name = String::from_utf8_lossy(rec.name()).into_owned();
        let seq = rec.sequence().as_ref();
        // The catalog scans period 1..=6; the pre-filter and post-filter drop period 1.
        let intervals =
            find_tandem_repeats(seq, PeriodRange::new(1, 6).unwrap(), &ScanParams::default());
        // **The narrowing that made this file block B2, now on ng's side of the
        // fence.** `TrfRecord` is trf-mod's parse shape and is `u32`; ng's
        // `RepeatInterval` is `u64` (spec §4). While this test lived in
        // `src/ssr/catalog/`, widening ng meant editing frozen production — which is
        // why it moved here. Now the conversion is ng's own, at ng's own boundary,
        // and production never notices.
        //
        // PANIC-FREE: the synthetic parity fixture is a few kb; a `u32` holds 4.29 Gb.
        // `expect` over `as` because a truncated coordinate would silently compare the
        // *wrong tract* against the golden catalog — a green test asserting nothing.
        let recs: Vec<TrfRecord> = catalog_prefilter(&intervals)
            .iter()
            .map(|iv| TrfRecord {
                start: u32::try_from(iv.start).expect("parity fixture coordinates fit u32"),
                end: u32::try_from(iv.end).expect("parity fixture coordinates fit u32"),
                period: u16::from(iv.period),
                score: iv.score,
            })
            .collect();
        out.extend(build_loci(recs, &name, seq, params));
    }
    out
}

/// Two loci are the *same tract* if they are on the same contig and their `[start, end)`
/// spans overlap — the parity-relevant match, tolerant of the small boundary/phase wobble
/// two different detectors place on the same repeat (spec §6.1). `(chrom, start, end)` for a
/// readable diff.
fn overlaps(a: &GoldenLocus, b: &PostFilteredLocus) -> bool {
    a.chrom == b.chrom && a.start < b.end && b.start < a.end
}

fn describe(chrom: &str, start: u32, end: u32, motif: &[u8]) -> String {
    format!("{chrom}:{start}-{end} {}", String::from_utf8_lossy(motif))
}

#[test]
fn scanner_reproduces_the_trf_mod_golden_catalog() {
    // The golden catalog (trf-mod → postprocess) and its build parameters.
    let catalog = golden_catalog();
    let params = catalog.settings;
    let golden = catalog.loci;

    // The scanner path over the same reference, with the same post-filter parameters.
    let scanner = scanner_loci(&fixture("synthetic_ref.fa"), &params);

    // Recall: every golden tract must be reproduced (overlap match). Track exact-coordinate
    // hits vs boundary/phase wobble, and scanner-only loci, so the diff stays reviewed.
    let mut exact = 0usize;
    let mut wobble = Vec::new();
    let mut missed = Vec::new();
    for g in &golden {
        match scanner.iter().find(|s| overlaps(g, s)) {
            Some(s) if (s.start, s.end, &s.motif) == (g.start, g.end, &g.motif) => {
                exact += 1;
            }
            Some(s) => wobble.push(format!(
                "{}  ~  {}",
                describe(&g.chrom, g.start, g.end, &g.motif),
                describe(&s.chrom, s.start, s.end, &s.motif)
            )),
            None => missed.push(describe(&g.chrom, g.start, g.end, &g.motif)),
        }
    }
    let extra: Vec<String> = scanner
        .iter()
        .filter(|s| !golden.iter().any(|g| overlaps(g, s)))
        .map(|s| describe(&s.chrom, s.start, s.end, &s.motif))
        .collect();

    let recovered = golden.len() - missed.len();
    let recall = recovered as f64 / golden.len() as f64;

    assert!(
        recall >= 0.99,
        "scanner recall {recall:.4} ({recovered}/{}) below the 0.99 parity bar (spec §9.2).\n\
         missed golden loci: {missed:#?}",
        golden.len(),
    );

    // Reviewed diff (spec §6.1): a handful of boundary/phase wobbles and scanner-only loci are
    // expected between two different detectors. Printed (visible with `--nocapture`), and the
    // scanner-only loci must be genuine STRs — pinned by count so the diff can't silently grow.
    eprintln!(
        "parity: {}/{} recall — {exact} exact, {} boundary/phase wobble, {} scanner-only.\n\
         wobble: {wobble:#?}\n  scanner-only (genuine STRs trf-mod's significance model rejected): {extra:#?}",
        recovered,
        golden.len(),
        wobble.len(),
        extra.len(),
    );
    assert!(
        extra.len() <= 1,
        "more scanner-only loci ({}) than the reviewed baseline of 1 (the ctg2 CCG tract): {extra:#?}",
        extra.len()
    );
}
