//! **Research switches for the GIAB alpha study (branch `research/giab-alpha`), not for main.**
//!
//! Each switch is an environment variable read once per process. Unset, every switch leaves the
//! caller exactly as main builds it, so a run without them is the baseline.

use std::sync::OnceLock;

fn read_f64(name: &str) -> Option<f64> {
    std::env::var(name).ok().map(|v| {
        v.parse::<f64>()
            .unwrap_or_else(|_| panic!("{name}={v} is not a number"))
    })
}

/// `PVC_EXP_OFF_GRID_SHARE`: the repeat-tract periodicity gate's share (main: 0.10).
pub fn off_grid_share() -> Option<f64> {
    static V: OnceLock<Option<f64>> = OnceLock::new();
    *V.get_or_init(|| read_f64("PVC_EXP_OFF_GRID_SHARE"))
}
