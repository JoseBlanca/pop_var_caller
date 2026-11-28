//! The read filters a walk applied, written into its psp header as provenance parameters, and
//! the key prefix a calling run uses to pick them back out (spec `run_streaming.md` §6.1).
//!
//! **Here, in `run`, because both ends of these keys are here**: the walk stage
//! ([`gatherer`](super::gatherer)) records them and the psp calling stage
//! ([`psp_caller`](super::psp_caller)) reads them back by prefix. Beside the filter's own config
//! in `read`, they made `read` import `psp` for this one function; in `psp`, the file-format
//! module would have had to know the read filter's settings. `run` already depends on both.

use crate::psp::ParameterValue;
use crate::psp::header::wire_float_of;
use crate::read::ReadFilterConfig;

#[cfg(test)]
/// Every provenance key [`read_filter_parameters`] can write, so a
/// caller replacing an earlier recording can clear the whole family first — the key set
/// shrinks when the mismatch filter turns off, and a stale floor beside an `"off"` is
/// exactly the misreading the conditional key exists to prevent.
pub const READ_FILTER_PROVENANCE_KEYS: [&str; 7] = [
    "read-filter-min-mapq",
    "read-filter-min-read-length-bp",
    "read-filter-drop-qc-fail",
    "read-filter-drop-duplicates",
    "read-filter-max-read-mismatch-fraction",
    "read-filter-mismatch-bq-floor",
    "read-filter-max-deletion-len-bp",
];

/// **The prefix every one of those keys begins with**, so a reader can pick the whole family
/// out of a header it did not write.
///
/// Both spellings are needed and they answer different questions.
/// `READ_FILTER_PROVENANCE_KEYS` is *what this build records* — a writer clearing its own
/// keys, or a test holding the set exhaustive, wants exactly those six. This is *what a read
/// filter's key looks like*, which is what a reader comparing two psps written by different
/// builds needs: a key a later ng added is a read-filter setting, and dropping it because this
/// build's list does not name it would make two differently-filtered walks compare equal.
pub const READ_FILTER_PROVENANCE_PREFIX: &str = "read-filter-";

/// How an off filter is spelled in the header: a value, not an absent key — an absent
/// key would read as "unrecorded", which is a different fact.
const OFF: &str = "off";

/// A read-filter policy as psp-header provenance parameters — recorded, never compared
/// (spec `run_streaming.md` §6.1). **Recorded so that a run can say when a cohort's
/// files were walked under different filters**, which is the one thing done with
/// them: `OpenPspCohort::open` does not compare them and does not refuse over them,
/// and the calling run's report prints a line naming the ones that disagree
/// ([`StoredSample::read_filters_the_walk_applied`](crate::run::StoredSample)).
/// A header that omitted them would let a cohort assembled from files walked under
/// different filters be called without a word about it.
///
/// One key per **configurable** filter, each readable by eye; an off filter is the
/// string `"off"`. The one absent-when-meaningless key is the mismatch base-quality
/// floor, which exists only while the mismatch-fraction filter is on. The
/// unconditional drops — secondary, supplementary and unmapped records, a read
/// whose CIGAR will not decode, and a read with no base qualities — are not settings,
/// so no key records them.
///
/// **The longest deletion a kept read may carry is recorded although no option sets
/// it** (`read-filter-max-deletion-len-bp`, added 2026-11-28). It follows the walk's
/// record cap, and a build that changes the cap changes which reads a walk kept. A psp
/// written before the key existed has no such key, and a cohort mixing the two gets a
/// report line naming the files that do not record it.
pub fn read_filter_parameters(filters: &ReadFilterConfig) -> Vec<(String, ParameterValue)> {
    // Exhaustive, so a filter added to the policy fails to compile here rather than going
    // unrecorded — which would make two differently-filtered walks digest identically in the
    // census, silently. It still holds away from the struct because every field is `pub`.
    let ReadFilterConfig {
        min_mapq,
        min_read_length,
        drop_qc_fail,
        drop_duplicate,
        max_read_mismatch_fraction,
        mismatch_bq_floor,
        max_deletion_len,
    } = *filters;

    let mut entries = vec![
        (
            "read-filter-min-mapq".to_string(),
            match min_mapq {
                Some(quality) => ParameterValue::Integer(i64::from(quality.get())),
                None => ParameterValue::String(OFF.to_string()),
            },
        ),
        (
            "read-filter-min-read-length-bp".to_string(),
            match min_read_length {
                // A length is written as a number while one fits TOML's signed
                // integer, and as its own digits in a string past that —
                // unreachable from any real configuration, but a writer must never
                // build a header its own reader refuses.
                Some(length) => match i64::try_from(length.get()) {
                    Ok(length) => ParameterValue::Integer(length),
                    Err(_) => ParameterValue::String(length.get().to_string()),
                },
                None => ParameterValue::String(OFF.to_string()),
            },
        ),
        (
            "read-filter-drop-qc-fail".to_string(),
            ParameterValue::Boolean(drop_qc_fail),
        ),
        (
            "read-filter-drop-duplicates".to_string(),
            ParameterValue::Boolean(drop_duplicate),
        ),
    ];
    match max_read_mismatch_fraction {
        Some(fraction) => {
            entries.push((
                "read-filter-max-read-mismatch-fraction".to_string(),
                ParameterValue::Float(wire_float_of(fraction.get())),
            ));
            entries.push((
                "read-filter-mismatch-bq-floor".to_string(),
                ParameterValue::Integer(i64::from(mismatch_bq_floor.get())),
            ));
        }
        None => entries.push((
            "read-filter-max-read-mismatch-fraction".to_string(),
            ParameterValue::String(OFF.to_string()),
        )),
    }
    entries.push((
        "read-filter-max-deletion-len-bp".to_string(),
        // The same fallback as the read length's: a string of digits past TOML's integer.
        match i64::try_from(max_deletion_len.get()) {
            Ok(length) => ParameterValue::Integer(length),
            Err(_) => ParameterValue::String(max_deletion_len.get().to_string()),
        },
    ));
    entries
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{BaseQual, Bp, MapQual, MismatchFraction};

    /// Every configurable filter lands under its own key with its **exact** value, from
    /// a policy whose two booleans differ — so transposing the qc-fail and duplicate
    /// sources, or nudging any number, fails here rather than passing everywhere (a
    /// mutation pass found exactly that: four value mutations survived a test that
    /// pinned only two of the six).
    #[test]
    fn read_filter_parameters_pin_every_filter_value() {
        let policy = ReadFilterConfig {
            min_mapq: Some(MapQual(17)),
            min_read_length: Some(Bp(31)),
            drop_qc_fail: true,
            drop_duplicate: false,
            max_read_mismatch_fraction: Some(
                MismatchFraction::try_new(0.25).expect("a quarter is in [0, 1]"),
            ),
            mismatch_bq_floor: BaseQual(13),
            max_deletion_len: Bp(4_321),
        };
        let entries: std::collections::BTreeMap<String, ParameterValue> =
            read_filter_parameters(&policy).into_iter().collect();

        assert_eq!(
            entries.get("read-filter-min-mapq"),
            Some(&ParameterValue::Integer(17))
        );
        assert_eq!(
            entries.get("read-filter-min-read-length-bp"),
            Some(&ParameterValue::Integer(31))
        );
        assert_eq!(
            entries.get("read-filter-drop-qc-fail"),
            Some(&ParameterValue::Boolean(true))
        );
        assert_eq!(
            entries.get("read-filter-drop-duplicates"),
            Some(&ParameterValue::Boolean(false))
        );
        assert_eq!(
            entries.get("read-filter-max-read-mismatch-fraction"),
            Some(&ParameterValue::Float(0.25))
        );
        assert_eq!(
            entries.get("read-filter-mismatch-bq-floor"),
            Some(&ParameterValue::Integer(13))
        );
        assert_eq!(
            entries.get("read-filter-max-deletion-len-bp"),
            Some(&ParameterValue::Integer(4_321))
        );
        assert_eq!(entries.len(), 7, "seven keys, nothing extra");
        for key in entries.keys() {
            assert!(
                READ_FILTER_PROVENANCE_KEYS.contains(&key.as_str()),
                "{key} is missing from the published key family"
            );
        }
    }

    /// **Every published key begins with the published prefix**, which is what lets a reader
    /// pick the whole family out of a header written by another build of ng.
    ///
    /// The two constants say different things — the list is *what this build records*, the
    /// prefix is *what a read filter's key looks like* — and a reader comparing two psps takes
    /// the prefix, so that a key a later ng added is still recognised as a read-filter setting
    /// rather than dropped. If they drift apart, that reader silently compares fewer settings
    /// than it reports.
    #[test]
    fn every_published_provenance_key_carries_the_published_prefix() {
        for key in READ_FILTER_PROVENANCE_KEYS {
            assert!(
                key.starts_with(READ_FILTER_PROVENANCE_PREFIX),
                "{key} does not begin with {READ_FILTER_PROVENANCE_PREFIX:?}, so a reader \
                 matching on the prefix would not compare it",
            );
        }
    }

    /// An off filter is spelled `"off"`, never omitted — an absent key would read as
    /// "unrecorded", a different fact — and the one deliberate omission is the
    /// base-quality floor beside an off mismatch filter, which would invite misreading.
    #[test]
    fn read_filter_parameters_spell_an_off_filter_rather_than_omitting_it() {
        let all_off = ReadFilterConfig {
            min_mapq: None,
            min_read_length: None,
            max_read_mismatch_fraction: None,
            ..ReadFilterConfig::default()
        };
        let entries: std::collections::BTreeMap<String, ParameterValue> =
            read_filter_parameters(&all_off).into_iter().collect();
        for key in [
            "read-filter-min-mapq",
            "read-filter-min-read-length-bp",
            "read-filter-max-read-mismatch-fraction",
        ] {
            assert_eq!(
                entries.get(key),
                Some(&ParameterValue::String("off".to_string())),
                "{key} must say off, not vanish"
            );
        }
        assert!(
            !entries.contains_key("read-filter-mismatch-bq-floor"),
            "a floor beside an off mismatch filter would invite misreading"
        );
        assert_eq!(entries.len(), 6, "six keys while the floor is meaningless");
    }

    /// A length past TOML's signed integer is recorded as its own digits in a string —
    /// unreachable from any real policy, but a wrapping-cast regression here would
    /// record a negative length silently.
    #[test]
    fn a_read_length_past_toml_integers_is_recorded_as_its_own_digits() {
        let policy = ReadFilterConfig {
            min_read_length: Some(Bp(u64::MAX)),
            ..ReadFilterConfig::default()
        };
        let entries: std::collections::BTreeMap<String, ParameterValue> =
            read_filter_parameters(&policy).into_iter().collect();
        assert_eq!(
            entries.get("read-filter-min-read-length-bp"),
            Some(&ParameterValue::String("18446744073709551615".to_string()))
        );
    }
}
