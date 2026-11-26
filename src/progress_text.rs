//! **How a run's messages write numbers, sizes and times**, and the target that marks a progress
//! line.
//!
//! Every stage that reports where it has got to — calling, parameter estimation, the command
//! that wrote the run — formats with these, so a count, a size or a time reads the same way in
//! every message. They sit at the crate root, below the library and the command line both, so a
//! library module can state a progress event without importing [`crate::cli`].
//!
//! **[`PROGRESS`] is a contract with the code that draws the line.** A library module emits
//! `info!(target: PROGRESS, …)`; the terminal layer in [`crate::cli::run_output`] recognises the
//! target and draws that event in place over the last one. The two sides live apart, and each
//! names the other.
//!
//! **Two duration styles, on purpose.** [`elapsed_text`] stops at whole seconds, which is what a
//! run's closing lines print; [`elapsed_text_to_the_millisecond`] goes below a second, which is
//! what parameter estimation's stage lines print for a stage that took `350ms`. Merging them would
//! change log text, which is the owner's decision.

use std::time::Duration;

/// The target a progress event carries, so the terminal can draw it in place:
/// `info!(target: PROGRESS, "calling: …")`. [`crate::cli::run_output`] is what reads it.
pub const PROGRESS: &str = "pvc::progress";

/// `1234567` as `1,234,567`.
#[must_use]
pub fn with_thousands(number: u64) -> String {
    let digits = number.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(digit);
    }
    grouped
}

/// A duration as a person reads it: `1h12m`, `12m05s` or `42s` — to the second, since a
/// progress line comes every few minutes at the least. Below a second it reads `0s`.
#[must_use]
pub fn elapsed_text(took: Duration) -> String {
    above_a_second(took).unwrap_or_else(|| format!("{}s", took.as_secs()))
}

/// A duration as `1h12m`, `12m05s`, `42s` or `350ms`: [`elapsed_text`], except that below a
/// second it reads in milliseconds.
#[must_use]
pub fn elapsed_text_to_the_millisecond(took: Duration) -> String {
    above_a_second(took).unwrap_or_else(|| format!("{}ms", took.as_millis()))
}

/// `took` at its largest unit, or `None` below one second.
fn above_a_second(took: Duration) -> Option<String> {
    let seconds = took.as_secs();
    let (hours, minutes, seconds) = (seconds / 3600, seconds / 60 % 60, seconds % 60);
    if hours > 0 {
        Some(format!("{hours}h{minutes:02}m"))
    } else if minutes > 0 {
        Some(format!("{minutes}m{seconds:02}s"))
    } else if seconds > 0 {
        Some(format!("{seconds}s"))
    } else {
        None
    }
}

/// A size in bytes as a person reads it: whole mebibytes below one gibibyte (`850 MiB`), where one
/// decimal of a gibibyte would read `0.0` for a cohort of a few hundred samples, and gibibytes to
/// one decimal above (`3.2 GiB`).
#[must_use]
pub fn size_in_bytes(bytes: u64) -> String {
    // Precision loss is irrelevant: a size to read, to one decimal.
    #[allow(clippy::cast_precision_loss)]
    let as_float = bytes as f64;
    if bytes >= 1 << 30 {
        format!("{:.1} GiB", as_float / f64::from(1u32 << 30))
    } else {
        format!("{:.0} MiB", as_float / f64::from(1u32 << 20))
    }
}

/// One memory field of the text of `/proc/self/status`, in kibibytes: `field` is the line's
/// label with its colon, such as `"VmHWM:"`. `None` where the line is missing or unreadable.
///
/// **Linux only has that file**, so on any other system the caller has no text to hand here and
/// prints its message without memory.
#[must_use]
pub fn kibibytes_in_status(status: &str, field: &str) -> Option<u64> {
    status
        .lines()
        .find_map(|line| line.strip_prefix(field))?
        .trim()
        .trim_end_matches("kB")
        .trim()
        .parse()
        .ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thousands_are_marked_with_commas() {
        assert_eq!(with_thousands(0), "0");
        assert_eq!(with_thousands(999), "999");
        assert_eq!(with_thousands(1_000), "1,000");
        assert_eq!(with_thousands(1_234_567), "1,234,567");
    }

    /// **The two duration styles differ below a second only.**
    #[test]
    fn durations_are_read_at_the_largest_unit() {
        for style in [elapsed_text, elapsed_text_to_the_millisecond] {
            assert_eq!(style(Duration::from_secs(42)), "42s");
            assert_eq!(style(Duration::from_secs(12 * 60 + 5)), "12m05s");
            assert_eq!(style(Duration::from_secs(3600 + 12 * 60 + 59)), "1h12m");
        }
        assert_eq!(elapsed_text(Duration::from_millis(350)), "0s");
        assert_eq!(
            elapsed_text_to_the_millisecond(Duration::from_millis(350)),
            "350ms"
        );
    }

    /// Whole mebibytes below a gibibyte, gibibytes to one decimal above.
    #[test]
    fn sizes_are_read_in_mib_and_gib() {
        assert_eq!(size_in_bytes(850 << 20), "850 MiB");
        assert_eq!(size_in_bytes(3 << 30), "3.0 GiB");
        assert_eq!(size_in_bytes(364_000_000), "347 MiB");
        assert_eq!(size_in_bytes(7_895_160_000), "7.4 GiB");
    }

    #[test]
    fn a_memory_field_is_read_from_the_status_text() {
        let status = "Name:\tpvc\nVmHWM:\t 9542041 kB\nVmRSS:\t 7864320 kB\n";
        assert_eq!(kibibytes_in_status(status, "VmHWM:"), Some(9_542_041));
        assert_eq!(kibibytes_in_status(status, "VmRSS:"), Some(7_864_320));
        assert_eq!(kibibytes_in_status(status, "VmPeak:"), None);
        assert_eq!(kibibytes_in_status("VmHWM:\t lots kB\n", "VmHWM:"), None);
    }
}
