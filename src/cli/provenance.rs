//! **What a run records about itself**: the command line that invoked it and the time it ran, for
//! the `##commandline` VCF header line and a psp's provenance block.
//!
//! Copied from production's `pop_var_caller::common` at promotion step C15, unchanged: ng's
//! commands had been calling production's copies, and production is being deleted.

use std::ffi::OsString;
use std::time::{SystemTime, UNIX_EPOCH};

/// Reconstruct the invoking process's command line as a single
/// space-joined string suitable for `WriterProvenance.command_line`
/// / VCF `##commandline=`. Lossy-UTF-8 for non-Unicode argv entries.
pub(crate) fn current_command_line() -> String {
    std::env::args_os()
        .map(|a: OsString| a.to_string_lossy().into_owned())
        .collect::<Vec<_>>()
        .join(" ")
}

/// Format the current UTC time as a TOML-compatible RFC3339 string
/// (`YYYY-MM-DDTHH:MM:SSZ`). Done by hand so we don't pull in a date
/// crate just for this one call; uses Howard Hinnant's civil-date
/// algorithm ([`civil_from_days`]) to turn days-since-epoch into
/// (year, month, day).
pub(crate) fn rfc3339_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let sod = secs % 86_400;
    let (y, m, d) = civil_from_days(days);
    let h = sod / 3600;
    let min = (sod % 3600) / 60;
    let s = sod % 60;
    format!("{y:04}-{m:02}-{d:02}T{h:02}:{min:02}:{s:02}Z")
}

/// Howard Hinnant's `civil_from_days` — translate days-since-epoch
/// (1970-01-01 = 0) into a proleptic Gregorian (year, month, day).
/// Verified by external reference; see
/// <http://howardhinnant.github.io/date_algorithms.html>.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365; // [0, 399]
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32; // [1, 12]
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn civil_from_days_matches_known_dates() {
        // Unix epoch
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        // 1970-12-31 (day 364, 1970 is not a leap year)
        assert_eq!(civil_from_days(364), (1970, 12, 31));
        // 2000-01-01 (day 10957 since epoch, common Y2K test)
        assert_eq!(civil_from_days(10957), (2000, 1, 1));
        // 2024-02-29 — leap day
        assert_eq!(civil_from_days(19782), (2024, 2, 29));
    }

    #[test]
    fn rfc3339_now_parses_as_toml_datetime() {
        let s = rfc3339_now();
        let _: toml::value::Datetime = s.parse().expect("must parse as toml Datetime");
    }
}
