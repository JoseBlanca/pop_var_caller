//! **The pass's random numbers** (`doc/devel/ng/spec/neighbour_evidence_phasing.md` §5).
//!
//! Phasing draws one state path a sample a round, and hides a few records of each sample before the
//! rounds start. Both read uniforms from a SplitMix64 stream whose seed is mixed from a constant and
//! numbers that depend on the input only — the run-wide number of the window's first core site, the
//! sample, the round — so the output is the same at any thread count and in any solving order.
//! `doc/devel/ng/research/neighbour_evidence/phaser_reference.py` implements the same stream, so the
//! Rust pass and the Python prototype can be compared exactly.

/// The constant every seed starts from. `phaser_reference.py`'s `SEED` must equal it, and its
/// purpose numbers must equal [`StreamPurpose`]'s, or the reference test fails.
const SEED: u64 = 0x6E65_6967_6862_6F75;

/// **What a stream is for**, mixed into its seed so two purposes never share one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StreamPurpose {
    /// The records of a sample hidden from the rounds for the fill check.
    HiddenRecords,
    /// The state path a sample draws in a round.
    Path,
    /// The extra paths a sample draws in a round for the phase's link test
    /// (`neighbour_evidence_phase_output.md` §2.3).
    ExtraPaths,
    /// The orientations of a sample's read-linked records along a round's paths, the round's own
    /// path first and then its extra paths (`neighbour_evidence_reads_in_draw.md` §4).
    Orientation,
}

/// **A SplitMix64 stream** of uniforms in `[0, 1)`.
#[derive(Clone, Debug)]
pub struct Uniforms {
    state: u64,
}

impl Uniforms {
    /// The stream for `purpose`, keyed by `window` (the window's key,
    /// `refine::window_key`), `sample` and `round` (0 for [`StreamPurpose::HiddenRecords`]).
    #[must_use]
    pub fn new(purpose: StreamPurpose, window: u64, sample: usize, round: usize) -> Self {
        let purpose = match purpose {
            StreamPurpose::HiddenRecords => 1,
            StreamPurpose::Path => 2,
            StreamPurpose::ExtraPaths => 3,
            StreamPurpose::Orientation => 4,
        };
        let mut state = SEED;
        for value in [purpose, window, sample as u64, round as u64] {
            state = finalise(state ^ value);
        }
        Self { state }
    }

    /// The next uniform: the top 53 bits of the next SplitMix64 output, over 2⁵³.
    pub fn next_uniform(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        (finalise(self.state) >> 11) as f64 / (1_u64 << 53) as f64
    }
}

/// SplitMix64's output function.
fn finalise(value: u64) -> u64 {
    let mut z = value;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// **The index a uniform `u` picks from `weights`**, in their order: the first whose cumulative
/// weight exceeds `u` times the total. Weights are non-negative and at least one is positive.
///
/// # Panics
///
/// If every weight is zero.
#[must_use]
pub fn pick(weights: &[f64], u: f64) -> usize {
    let total: f64 = weights.iter().sum();
    assert!(total > 0.0, "a draw needs a positive weight");
    let target = u * total;
    let mut cumulative = 0.0;
    let mut last_positive = 0;
    for (index, &weight) in weights.iter().enumerate() {
        if weight > 0.0 {
            last_positive = index;
        }
        cumulative += weight;
        if cumulative > target {
            return index;
        }
    }
    last_positive
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The stream is a function of its key**, and different keys give different streams.
    #[test]
    fn the_stream_depends_on_its_key_only() {
        let mut one = Uniforms::new(StreamPurpose::Path, 7, 3, 2);
        let mut again = Uniforms::new(StreamPurpose::Path, 7, 3, 2);
        let mut other = Uniforms::new(StreamPurpose::Path, 7, 3, 3);
        let first: Vec<f64> = (0..5).map(|_| one.next_uniform()).collect();
        let second: Vec<f64> = (0..5).map(|_| again.next_uniform()).collect();
        let third: Vec<f64> = (0..5).map(|_| other.next_uniform()).collect();
        assert_eq!(first, second);
        assert_ne!(first, third);
        assert!(first.iter().all(|&u| (0.0..1.0).contains(&u)));
    }

    /// **A pick follows the cumulative weights**, and a zero weight is never picked.
    #[test]
    fn a_pick_follows_the_cumulative_weights() {
        let weights = [0.0, 1.0, 0.0, 3.0];
        assert_eq!(pick(&weights, 0.0), 1);
        assert_eq!(pick(&weights, 0.24), 1);
        assert_eq!(pick(&weights, 0.26), 3);
        assert_eq!(pick(&weights, 0.999_999), 3);
    }
}
