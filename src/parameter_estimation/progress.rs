//! **Where parameter estimation has got to**, printed to stderr.
//!
//! Estimating a whole-genome cohort's parameters takes hours and writes nothing until the
//! parameters file at the end, so a person watching it cannot tell a long fit from a hung one.
//! Each stage says when it starts and when it ends, and the stages that loop — the SNP/indel
//! fit's passes, the repeat-tract fit's walks — say where they are **at most once a minute**:
//!
//! ```text
//! estimating: SNP/indel fit, start 1 of 3, pass 17 of at most 200: log-likelihood moved 3.1e3 (stops below 2.5e2), largest parameter move 2.1e-3 (stops below 1e-4); 34s a pass, 9m38s into this stage
//! ```
//!
//! **A minute rather than a count of items**, because what an item costs spans orders of
//! magnitude across the inputs this caller takes: a pass over one sample's small panel takes
//! milliseconds and one over a thousand whole genomes takes minutes. A count that suits one end
//! floods or starves the other; a clock suits both, and a small run prints only its stage lines.
//!
//! Nothing here touches what the estimation computes.

use std::sync::Mutex;
use std::time::{Duration, Instant};

/// The shortest time between two "where it is" lines of one stage.
const BETWEEN_LINES: Duration = Duration::from_secs(60);

/// One stage's clock: when it began and when it last said where it was.
///
/// **`Sync`**, so the walks of a parallel stage can share one and still print at most once a
/// minute between them.
pub(crate) struct StageProgress {
    started: Instant,
    last_line: Mutex<Instant>,
}

impl StageProgress {
    /// Start a stage, announcing it with `line`.
    pub(crate) fn begin(line: impl std::fmt::Display) -> Self {
        eprintln!("estimating: {line}");
        let now = Instant::now();
        Self {
            started: now,
            last_line: Mutex::new(now),
        }
    }

    /// Say where the stage is, if a minute has passed since its last line. `line` is only
    /// formatted when it is printed, and is given the time since the stage began.
    pub(crate) fn now_and_then(&self, line: impl FnOnce(&str) -> String) {
        let mut last = self
            .last_line
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if last.elapsed() < BETWEEN_LINES {
            return;
        }
        *last = Instant::now();
        drop(last);
        eprintln!("estimating: {}", line(&self.time_into_the_stage()));
    }

    /// Say something unconditionally — a milestone inside the stage, such as a start that ended.
    pub(crate) fn always(&self, line: impl FnOnce(&str) -> String) {
        eprintln!("estimating: {}", line(&self.time_into_the_stage()));
    }

    /// The time since the stage began, as `9m38s into this stage`.
    fn time_into_the_stage(&self) -> String {
        format!("{} into this stage", duration(self.started.elapsed()))
    }
}

/// A duration as `1h12m`, `12m05s`, `42s` or `350ms`.
pub(crate) fn duration(took: Duration) -> String {
    let seconds = took.as_secs();
    let (hours, minutes, seconds) = (seconds / 3600, seconds / 60 % 60, seconds % 60);
    if hours > 0 {
        format!("{hours}h{minutes:02}m")
    } else if minutes > 0 {
        format!("{minutes}m{seconds:02}s")
    } else if seconds > 0 {
        format!("{seconds}s")
    } else {
        format!("{}ms", took.as_millis())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn durations_are_read_at_the_largest_unit() {
        assert_eq!(duration(Duration::from_millis(350)), "350ms");
        assert_eq!(duration(Duration::from_secs(42)), "42s");
        assert_eq!(duration(Duration::from_secs(12 * 60 + 5)), "12m05s");
        assert_eq!(duration(Duration::from_secs(3600 + 12 * 60 + 59)), "1h12m");
    }

    /// A stage that has only just begun prints no "where it is" line: the minute runs from its
    /// start, so a fast stage says only that it began and ended.
    #[test]
    fn a_line_waits_a_minute_from_the_stage_start() {
        let stage = StageProgress::begin("a test stage");
        let mut formatted = false;
        stage.now_and_then(|_| {
            formatted = true;
            String::new()
        });
        assert!(!formatted);
    }
}
