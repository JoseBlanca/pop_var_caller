//! Does the transcription still compute production's window? A differential against the
//! accumulator it was copied from.
//!
//! **Production's side is frozen since promotion step C13.** Until then each of the 200 generated
//! streams was fed to production's accumulator beside ng's and every window compared. Production is
//! being deleted, and its 447,581 windows are too many to commit one by one, so what it emitted on
//! each stream was recorded at commit `d9e7b076` as a count and a digest
//! (`testdata/production_windows_by_seed.tsv`), and ng's stream is reduced the same way and compared.
//! **What that gives up:** a divergence still fails, at the seed it happens in, but no longer names
//! the window or the field; bisecting it needs ng's windows printed and read by hand. The prose below
//! describes the comparison as it was built.
//!
//! [`accumulator`](super::accumulator) claims to be `SlidingWindowCoverageAccumulator`
//! (`src/sample_summary/coverage.rs`) with its vocabulary changed and its arithmetic untouched.
//! The transcribed unit tests check the arithmetic against hand-computed means; this checks it
//! against the original, on streams neither test was written for — jumping positions, repeated
//! positions, `N` bases in both cases, contig changes, and window widths from one base to
//! wider than the 500 the run configures, over contigs long enough that a 500-base window
//! genuinely slides rather than swallowing the contig whole.
//!
//! **This has narrowed to the windows, and will not narrow further.** ng's accumulator has
//! both of its additions now: the floor (spec §3.3), which is switched off here so the two
//! sides stay comparable, and a depth bin width fitted to each sample (spec §3.4), which
//! production does not have and which ends the histogram half of the comparison — two
//! histograms cut on different axes are not comparable, whatever the windows behind them.
//! **What the windows themselves report is still production's, and is what this holds.** Four
//! behaviours the transcribed unit tests used to leave entirely to this file — the overflow
//! column's boundary, the GC clamp, lowercase bases, and the closing frontier — were given unit
//! tests of their own before the narrowing, which is why the narrowing costs nothing.

use super::{WindowCoverageAccumulator, WindowCoverageConfig};
use crate::types::{ContigId, Position};

/// **Production's windows on one stream, as recorded**: the width that stream drew, how many
/// windows production emitted, and the digest of all of them — see the fixture's header.
struct ProductionStream {
    window_bp: u32,
    windows: usize,
    digest: u64,
}

/// Production's recorded stream for `seed`.
fn production_stream(seed: u64) -> ProductionStream {
    let line = include_str!("testdata/production_windows_by_seed.tsv")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .find(|line| line.split('\t').next() == Some(seed.to_string().as_str()))
        .unwrap_or_else(|| panic!("no recorded production stream for seed {seed}"));
    let fields: Vec<&str> = line.split('\t').collect();
    let [_, window_bp, windows, digest] = fields[..] else {
        panic!("a recorded stream has four columns: {line}");
    };
    ProductionStream {
        window_bp: window_bp.parse().expect("a window width"),
        windows: windows.parse().expect("a window count"),
        digest: u64::from_str_radix(digest, 16).expect("a digest"),
    }
}

/// Knuth's multiplier for a 64-bit linear congruential generator, with the increment from the
/// same table (`MMIX`).
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;
const LCG_INCREMENT: u64 = 1_442_695_040_888_963_407;
/// An LCG's low bits cycle short, so values are taken from the high half.
const LCG_HIGH_BITS_SHIFT: u32 = 33;
/// Knuth's multiplicative-hash constant, used only to spread consecutive seeds apart.
const SEED_SPREAD: u64 = 2_654_435_761;

/// A linear congruential generator, so ng is fed exactly the stream production's windows were recorded on, and the
/// stream is the same on every machine and every run. Numbers, not randomness, is the point —
/// a failure has to be reproducible from its seed alone.
struct DeterministicStream(u64);

impl DeterministicStream {
    fn from_seed(seed: u64) -> Self {
        Self(seed.wrapping_mul(SEED_SPREAD).wrapping_add(12_345))
    }

    fn next_value(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(LCG_MULTIPLIER)
            .wrapping_add(LCG_INCREMENT);
        self.0 >> LCG_HIGH_BITS_SHIFT
    }

    fn next_below(&mut self, bound: u64) -> u64 {
        self.next_value() % bound
    }
}

/// How many windows the 200 streams compare. Asserted rather than described, so that a change
/// which quietly stops the streams producing windows cannot leave this test passing on nothing.
/// It is the count of non-`N` covered positions the generator draws, and so it does not move
/// when the window widths change — widening those is free coverage, and this number does not
/// notice it.
const WINDOWS_COMPARED: u64 = 447_581;

/// Every window ng emits equals production's, bit for bit — through a digest of the whole stream
/// since promotion step C13.
///
/// **Bit equality, not a tolerance.** The two implementations do the same arithmetic in the
/// same order, so any difference at all is a transcription slip rather than rounding, and a
/// tolerance would hide exactly the class of error this test exists to catch.
///
/// The mutation that proves it discriminates: pulling the window's left edge **in** by one base
/// (`centre - half_window_bp` to `centre - (half_window_bp - 1)`, which narrows the window
/// rather than widening it) fails the digest at the first seed, and fails three of the
/// eleven transcribed tests. At `window_bp = 1` that literal mutation underflows instead of
/// producing a wrong window, because the half-window is zero there.
#[test]
fn the_transcription_matches_production_on_streams_neither_test_was_written_for() {
    let mut windows_compared = 0u64;

    for seed in 0..200u64 {
        let mut stream = DeterministicStream::from_seed(seed);
        // 1 to 600 bases wide: 1 is the degenerate window that holds only its own centre, and
        // the range spans the 500 the run configures (`super::WINDOW_BP`), which the
        // transcribed unit tests never reach.
        let window_bp = 1 + stream.next_below(600) as u32;
        let config = WindowCoverageConfig {
            window_bp,
            gc_bins: 5,
            depth_bins: 17,
            // The floor at 1 is the floor switched off — every window holds at least its own
            // centre — which is what keeps this comparable with production, whose window has
            // no floor at all. The floor is ng's own and is unit-tested.
            min_window_positions: 1,
            depth_scale_windows: 10_000,
            depth_range_in_medians: 10.0,
        };
        let mut transcription = WindowCoverageAccumulator::new(config);
        let mut transcribed_windows = Vec::new();

        for contig in 0..3u32 {
            let mut position = 1u64;
            // Up to 2,000 covered positions at gaps of 0 to 8 bases, so a contig spans some
            // thousands of bases and a 500-base window slides along it instead of covering it
            // whole.
            let covered_positions = 1 + stream.next_below(2_000);
            for _ in 0..covered_positions {
                // A gap of 0 repeats the previous position — the walk emits at most one
                // record a position, but the accumulator's order rule is non-decreasing, and
                // a repeat is the boundary of it.
                position += stream.next_below(9);
                let alphabet = *b"ACGTNgcn";
                let reference_base = alphabet[stream.next_below(alphabet.len() as u64) as usize];
                let depth = stream.next_below(40) as u32;
                transcription.observe(ContigId(contig), Position(position), reference_base, depth);
                while let Some(window) = transcription.pop_ready() {
                    transcribed_windows.push(window);
                }
            }
        }

        // Production's histogram cells are not comparable — it cut its depth axis at a fixed width
        // and ng cuts at one fitted to the sample — but two things about ng's still are, and
        // neither has anything to do with the depth axis: that a stream of real windows produces a
        // histogram at all, and that every window emitted was folded into it.
        let (transcribed_tail, transcribed_histogram) = transcription.finish();
        let transcribed_histogram = transcribed_histogram.fitted().expect(
            "every window here clears a floor of 1 and the depths are positive, so a width is \
             fitted and a histogram comes back",
        );
        transcribed_windows.extend(transcribed_tail);

        let production = production_stream(seed);
        assert_eq!(
            production.window_bp, window_bp,
            "seed {seed}: the stream drew a different window width from the one recorded, so the \
             generator has changed and the recording no longer describes it",
        );
        assert_eq!(
            transcribed_windows.len(),
            production.windows,
            "seed {seed}, window_bp {window_bp}: different number of windows emitted",
        );
        // Reduced exactly as production's were: contig and centre as `u32`, then the mean depth
        // and GC fraction as the bit patterns of their `f32`s, all little-endian, in emission
        // order. Converting the centre rather than casting it, so that a wider stream fails here
        // instead of silently digesting a wrapped one.
        let mut bytes = Vec::with_capacity(transcribed_windows.len() * 16);
        for (at, window) in &transcribed_windows {
            bytes.extend(at.contig.get().to_le_bytes());
            bytes.extend(
                u32::try_from(at.position.get())
                    .expect("fixture positions stay inside u32")
                    .to_le_bytes(),
            );
            bytes.extend(window.mean_depth.to_bits().to_le_bytes());
            bytes.extend(window.gc_fraction.to_bits().to_le_bytes());
        }
        let digest = bytes
            .iter()
            .fold(0xcbf2_9ce4_8422_2325_u64, |digest, &byte| {
                (digest ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
            });
        assert_eq!(
            digest, production.digest,
            "seed {seed}, window_bp {window_bp}: some window's contig, centre, mean depth or GC \
             fraction differs from production's — print ng's windows for this seed to find which",
        );
        windows_compared += transcribed_windows.len() as u64;
        assert!(
            transcribed_windows
                .iter()
                .all(|(_, window)| !window.is_absent()),
            "seed {seed}: the floor is set to 1 here so that this stays a comparison with \
             production, whose window has no floor — an absent window means it bound",
        );
        assert_eq!(
            transcribed_histogram.windows_folded as usize,
            transcribed_windows.len(),
            "seed {seed}: every window emitted must reach the histogram, whatever axis it is \
             cut on",
        );
    }

    assert_eq!(
        windows_compared, WINDOWS_COMPARED,
        "the 200 streams compare a fixed number of windows; a change here means the fixture \
         moved, not that the accumulator is wrong",
    );
}
