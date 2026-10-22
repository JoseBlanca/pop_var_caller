//! **Research switches for the GIAB alpha study (branch `research/giab-alpha`), not for main.**
//!
//! Each switch is an environment variable read once per process. Unset, every switch leaves the
//! caller exactly as main builds it, so a run without them is the baseline.
//!
//! The first, `PVC_EXP_OFF_GRID_SHARE`, set the repeat-tract periodicity check's share; it went
//! with the check itself on 2026-10-22 (`doc/devel/ng/spec/candidate_alleles_ssr.md` §7).
//!
//! **The windowed duplication rate** ([`crate::paralog::windowed_prior`]) is switched on by the
//! command line's `--hidden-dup-window`, not by an environment variable. Two of its knobs stay
//! here:
//!
//! - `PVC_EXP_PARALOG_SHRINK` — how hard each window's rate is pulled towards the run-wide one,
//!   in records. Default `10`.
//! - `PVC_EXP_PARALOG_CUT` — `fdr` (default) or `posterior`; see
//!   [`WindowedCut`](crate::paralog::windowed_prior::WindowedCut).
//!
//! Neither is read on a run whose window is `0`. A value that does not parse stops the run rather
//! than being read as unset, so a typo cannot quietly produce a different run.

use std::sync::OnceLock;

use crate::paralog::windowed_prior::{WindowedCut, WindowedPriorConfig};

fn read_var(name: &str) -> Option<String> {
    std::env::var(name).ok().filter(|value| !value.is_empty())
}

/// The windowed duplication rate's knobs for a window of `window_bases`, or `None` where the
/// window is `0` — the one run-wide rate.
///
/// # Panics
///
/// If a switch is set to something it cannot be: a pull that is not a finite number at or above
/// zero, or a cut other than `fdr` or `posterior`.
pub fn windowed_paralog_prior(window_bases: u64) -> Option<WindowedPriorConfig> {
    static KNOBS: OnceLock<(f64, WindowedCut)> = OnceLock::new();
    if window_bases == 0 {
        return None;
    }
    let (shrink_records, cut) = *KNOBS.get_or_init(|| {
        let shrink_records = read_var("PVC_EXP_PARALOG_SHRINK").map_or(10.0, |v| {
            v.parse::<f64>()
                .ok()
                .filter(|k| k.is_finite() && *k >= 0.0)
                .unwrap_or_else(|| {
                    panic!("PVC_EXP_PARALOG_SHRINK={v} is not a finite number of records >= 0")
                })
        });
        let cut = match read_var("PVC_EXP_PARALOG_CUT").as_deref() {
            None | Some("fdr") => WindowedCut::TailFdr,
            Some("posterior") => WindowedCut::Posterior,
            Some(other) => panic!("PVC_EXP_PARALOG_CUT={other} is neither fdr nor posterior"),
        };
        (shrink_records, cut)
    });
    Some(WindowedPriorConfig {
        window_bases,
        shrink_records,
        cut,
    })
}
