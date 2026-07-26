//! The STR tract delimiter's benchmark — `SsrUnitRobustAligner` (algorithm 4u), the default
//! delimiter [`SsrGenerator::next_locus`] drives once per kept read.
//!
//! **Why this exists.** Two sampling profiles taken on 2026-07-26
//! (`doc/devel/reports/reviews/perf_ng-ssr-locus-observations_2026-07-26.md`) put
//! `classify::delimit` at the top of the self-time ranking in both regimes the caller serves:
//! **36.0%** of a 51-sample shallow tomato cohort (~7 reads per covered locus) and **67.0%** of a
//! deep single sample (HG002 300×, ~412 reads per covered locus). The DP findings against it are
//! gated at ~10% *on the delimiter*, which Amdahl translates to 3.6–6.7% end to end — inside the
//! run-to-run spread that has already produced one documented thermal-drift false positive on this
//! path. So a delimiter-level instrument is the prerequisite for touching the DP at all.
//!
//! Two groups, because the findings split on which one they trade against:
//!
//! 1. `ng_ssr_delimiter/frame` — one `align` per iteration over a named (period, tract) frame.
//!    This is the group that resolves per-cell work (hoisting a per-column constant, the emission
//!    table lookup, the backpointer store width).
//! 2. `ng_ssr_delimiter/depth` — `N` reads against **one** frame with **one** shared scratch per
//!    iteration, `N ∈ {1, 7, 30, 100, 412}` spanning the two measured regimes. This is the group
//!    that resolves per-*locus* work (anything hoisted out of the per-read loop amortizes here and
//!    nowhere else), and it is reported per read via `Throughput::Elements`.
//!
//! Three things in here are load-bearing rather than decoration, each pinned to a real trap:
//!
//! - **`black_box` on both sides.** Under `lto = "fat"` a `const` frame lets the compiler hoist the
//!   geometry and flank resolution out of the timed loop — which is exactly the work some of the
//!   candidate fixes hoist, so an un-`black_box`ed bench would report those wins as already
//!   present.
//! - **A verification assertion in the timed body.** `delimit` has a silent zero-work path: it
//!   returns `None` the moment `reference_len == 0`, and `align` maps that to
//!   `RepeatSpan::Unanchored`. A fixture built with an empty or mis-sliced frame would produce a
//!   plausible-looking, meaninglessly fast number. Asserting the *measured tract length* is what
//!   makes "this measured the DP" survive the next refactor.
//! - **More than one frame per group.** `UnitRobustScratch` is grow-and-keep and deliberately not
//!   re-zeroed between reads, so a single-input microbench measures a permanently-warm L1 that real
//!   loci never see. The `frame` group therefore sweeps distinct (period, tract) shapes.
//!
//! The fixtures are built here rather than shared with `examples/ng_ssr_synthetic_bakeoff.rs`
//! (which has a richer scenario grid behind an `[[example]]` a bench cannot import): a bench needs
//! representative shapes, not the bake-off's accuracy axes, and duplicating four constructors is
//! cheaper than making test support public. If the two ever need to *provably* agree, lift the
//! bake-off's `Scenario` into a shared module and delete these.

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

use pop_var_caller::ng::alignment::ssr_unit_robust::{SsrUnitRobustAligner, UnitRobustScratch};
use pop_var_caller::ng::alignment::{
    BestPathAligner, PerQualityEmission, ReadBases, RepeatContext, RepeatGeometry, RepeatSpan,
    StutterModel,
};
use pop_var_caller::ng::types::{Bp, Motif};

/// The flank the generator fetches either side of the tract, per
/// `SsrGeneratorConfig::default()` (= `DEFAULT_BUNDLE_THRESHOLD`). The DP's reference frame is
/// `tract + 2 * FLANK`, so this constant, not the tract, sets the matrix width — at the measured
/// median tract (8 bp tomato / 14 bp HG002) about 88% of the columns are flank.
const FLANK: usize = 30;

/// Phred quality every fixture base carries. A single value keeps the emission lookup's *input*
/// constant so the bench measures the DP rather than the table; the DP's branches are exercised by
/// the base mismatches instead.
const QUAL: u8 = 40;

/// One delimiter fixture: a reference frame, a read over it, and the tract length the DP must
/// measure. `expected` is asserted inside the timed body — see the module docs.
struct Frame {
    label: String,
    reference: Vec<u8>,
    read: Vec<u8>,
    quality: Vec<u8>,
    geometry: RepeatGeometry,
    expected: Option<u64>,
}

/// Deterministic flank bases — not a repeat of the motif, so the tract boundaries are unambiguous
/// and the anchor test has real evidence to count. Two different sequences either side, so a
/// left/right transposition in any candidate fix shows up as a wrong measurement.
fn left_flank() -> Vec<u8> {
    b"GTCAGTCAGTCAGTCAGTCAGTCAGTCAGT"[..FLANK].to_vec()
}

fn right_flank() -> Vec<u8> {
    b"TTGCATTGCATTGCATTGCATTGCATTGCA"[..FLANK].to_vec()
}

/// `units` copies of `motif`, truncated to `len` bases so a tract length that is not a whole
/// multiple of the period is reachable (the common case in real catalogs).
fn tract(motif: &[u8], len: usize) -> Vec<u8> {
    motif.iter().copied().cycle().take(len).collect()
}

/// A read spanning the whole frame with `mismatches` evenly-spaced substituted bases — a complete
/// observation of the reference allele, with enough noise that the mismatch/match select in the
/// cell is not a constant branch.
fn spanning(reference: &[u8], mismatches: usize) -> Vec<u8> {
    let mut read = reference.to_vec();
    if mismatches > 0 {
        let step = read.len() / (mismatches + 1);
        for i in 0..mismatches {
            let at = step * (i + 1);
            read[at] = match read[at] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
    }
    read
}

/// The frame set the `frame` group sweeps: the two measured medians (period 1 tract 8, period 2
/// tract 14), a long tract where the per-tract-column findings should pay most, and a hexamer
/// where the runtime `%` in the slip emission is widest.
fn frames() -> Vec<Frame> {
    let mut out = Vec::new();
    for (motif_bytes, tract_len, mismatches) in [
        (&b"A"[..], 8usize, 1usize), // tomato median: mononucleotide, short
        (&b"AC"[..], 14, 2),         // HG002 median: dinucleotide
        (&b"AGC"[..], 30, 2),        // trinucleotide, long tract
        (&b"AGGCTC"[..], 42, 3),     // hexamer, longest tract: widest slip emission
    ] {
        let motif = Motif::new(motif_bytes).expect("valid period");
        let mut reference = left_flank();
        reference.extend_from_slice(&tract(motif_bytes, tract_len));
        reference.extend_from_slice(&right_flank());

        // A read that measures the reference allele exactly.
        let read = spanning(&reference, mismatches);
        let quality = vec![QUAL; read.len()];
        out.push(Frame {
            label: format!("p{}_t{}", motif_bytes.len(), tract_len),
            reference: reference.clone(),
            read,
            quality,
            geometry: RepeatGeometry {
                left_flank_len: Bp(FLANK as u64),
                right_flank_len: Bp(FLANK as u64),
                motif,
            },
            expected: Some(tract_len as u64),
        });

        // A read carrying one extra unit — the whole-unit slip route, which is the one the
        // per-tract-column emission work sits on and which a reference-only fixture never enters.
        let expanded_len = tract_len + motif_bytes.len();
        let mut expanded = left_flank();
        expanded.extend_from_slice(&tract(motif_bytes, expanded_len));
        expanded.extend_from_slice(&right_flank());
        let quality = vec![QUAL; expanded.len()];
        out.push(Frame {
            label: format!("p{}_t{}_expanded", motif_bytes.len(), tract_len),
            reference,
            read: spanning(&expanded, mismatches),
            quality,
            geometry: RepeatGeometry {
                left_flank_len: Bp(FLANK as u64),
                right_flank_len: Bp(FLANK as u64),
                motif,
            },
            expected: Some(expanded_len as u64),
        });
    }
    out
}

/// Align one read and assert the DP measured what the fixture says it should. Returns the span so
/// the caller can `black_box` it.
#[inline]
fn align_checked(
    aligner: &SsrUnitRobustAligner<PerQualityEmission>,
    frame: &Frame,
    stutter: &StutterModel,
    scratch: &mut UnitRobustScratch,
) -> RepeatSpan {
    let bases = ReadBases::try_new(black_box(&frame.read), black_box(&frame.quality))
        .expect("fixture read and quality are the same length");
    let context = RepeatContext {
        geometry: &frame.geometry,
        stutter,
    };
    let span = aligner.align(bases, black_box(&frame.reference), context, scratch);
    assert_eq!(
        span.measured_length(),
        frame.expected,
        "fixture {} did not measure its own tract — the bench is not measuring the DP",
        frame.label
    );
    span
}

/// Per-read cost by frame shape. The group every per-cell finding is judged on.
fn bench_frame(c: &mut Criterion) {
    let aligner = SsrUnitRobustAligner::new(PerQualityEmission::new());
    let stutter = StutterModel::hipstr_shipped();
    let frames = frames();

    let mut group = c.benchmark_group("ng_ssr_delimiter/frame");
    for frame in &frames {
        // One scratch per frame, built outside the timed loop: grow-and-keep is the production
        // shape (the generator holds one for the whole run).
        let mut scratch = UnitRobustScratch::new();
        group.bench_function(BenchmarkId::from_parameter(&frame.label), |b| {
            b.iter(|| black_box(align_checked(&aligner, frame, &stutter, &mut scratch)));
        });
    }
    group.finish();
}

/// Per-read cost at a locus of `depth` reads, sharing one scratch — the group every per-*locus*
/// finding is judged on, and the only one where the two measured regimes (~7 and ~412 reads per
/// covered locus) are distinguishable.
fn bench_depth(c: &mut Criterion) {
    let aligner = SsrUnitRobustAligner::new(PerQualityEmission::new());
    let stutter = StutterModel::hipstr_shipped();
    // The HG002 median shape: dinucleotide, 14 bp tract. Two reads alternate (reference allele and
    // one-unit expansion) so the locus is not a single input repeated.
    let frames = frames();
    let locus: Vec<&Frame> = frames
        .iter()
        .filter(|f| f.label.starts_with("p2_t14"))
        .collect();
    assert_eq!(locus.len(), 2, "expected the reference and expanded reads");

    let mut group = c.benchmark_group("ng_ssr_delimiter/depth");
    for depth in [1usize, 7, 30, 100, 412] {
        group.throughput(Throughput::Elements(depth as u64));
        let mut scratch = UnitRobustScratch::new();
        group.bench_function(BenchmarkId::from_parameter(depth), |b| {
            b.iter(|| {
                for index in 0..depth {
                    let frame = locus[index % locus.len()];
                    black_box(align_checked(&aligner, frame, &stutter, &mut scratch));
                }
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_frame, bench_depth);
criterion_main!(benches);
