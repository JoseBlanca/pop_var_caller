//! **The fixed order a repeat-tract stratum takes its samples in** (`fit_precision.md` §4.4).
//!
//! A large cohort's repeat-tract strata are fitted on a subset of samples — the first ones in this
//! order — grown until the stratum's slippage level is measured precisely enough. The order is a
//! ranking of the samples by a hash of each sample's **name** and a fixed seed, so:
//!
//! - it is the same on every machine and in every run;
//! - it does not depend on the order the psps were listed in, since a sample's place depends on its
//!   name alone;
//! - it is chosen without looking at depth, genotypes or anything else in the data.
//!
//! One order serves every stratum, so a sample is in every stratum's first subset or in none
//! (spec §4.4, decided).
//!
//! The hash is the census's, xxh3 with a seed ([`hash_position`](super::loci::hash_position)),
//! taken over the name's UTF-8 bytes and nothing else. The census feeds its names through Rust's
//! `Hash` for `str`, which also appends a byte whose choice the standard library has not settled;
//! here the bytes go to xxh3 directly, so the order cannot move with a toolchain upgrade.

use xxhash_rust::xxh3::xxh3_64_with_seed;

/// **The seed the sample order is drawn with.** Any fixed value would do; this one is fixed so that
/// the same cohort takes the same subset in every run. Changing it changes which samples every
/// large cohort's repeat-tract strata are fitted on.
pub const SAMPLE_ORDER_SEED: u64 = 0x5EED_0F5A_4D50_1E5D;

/// The value a sample is ranked by: a hash of its name and `seed`.
fn hash_sample(name: &str, seed: u64) -> u64 {
    xxh3_64_with_seed(name.as_bytes(), seed)
}

/// **The samples in their fixed order**: the indices of `names`, ranked by the hash of each name
/// and `seed`, lowest first.
///
/// Two different names whose hashes are equal are ranked by name, so the order of the names it
/// returns is the same whatever order `names` arrives in. A run never holds two samples of one name
/// (`RunError::SampleAppearsTwice` refuses the cohort).
pub fn sample_order<S: AsRef<str>>(names: &[S], seed: u64) -> Vec<usize> {
    let mut keyed: Vec<(u64, &str, usize)> = names
        .iter()
        .enumerate()
        .map(|(index, name)| (hash_sample(name.as_ref(), seed), name.as_ref(), index))
        .collect();
    keyed.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
    keyed.into_iter().map(|(_, _, index)| index).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The names `sample_order` puts first to last.
    fn ordered_names<'a>(names: &[&'a str], seed: u64) -> Vec<&'a str> {
        sample_order(names, seed)
            .into_iter()
            .map(|index| names[index])
            .collect()
    }

    /// **The order of the names is the same whatever order they arrive in**: every rotation and the
    /// reversal of thirty names give the same sequence of names.
    #[test]
    fn the_order_does_not_depend_on_the_order_the_names_arrive_in() {
        let owned: Vec<String> = (0..30).map(|i| format!("SRR50799{i:02}")).collect();
        let names: Vec<&str> = owned.iter().map(String::as_str).collect();
        let expected = ordered_names(&names, SAMPLE_ORDER_SEED);
        for shift in 1..names.len() {
            let mut rotated = names.clone();
            rotated.rotate_left(shift);
            assert_eq!(
                ordered_names(&rotated, SAMPLE_ORDER_SEED),
                expected,
                "shift {shift}"
            );
        }
        let mut reversed = names.clone();
        reversed.reverse();
        assert_eq!(ordered_names(&reversed, SAMPLE_ORDER_SEED), expected);
    }

    /// **The order is a permutation of the samples**, and not the order they arrived in: thirty
    /// names come back once each, shuffled.
    #[test]
    fn every_sample_is_ranked_once() {
        let owned: Vec<String> = (0..30).map(|i| format!("sample_{i}")).collect();
        let order = sample_order(&owned, SAMPLE_ORDER_SEED);
        let mut seen = order.clone();
        seen.sort_unstable();
        assert_eq!(seen, (0..30).collect::<Vec<_>>());
        assert_ne!(order, (0..30).collect::<Vec<_>>());
        assert!(sample_order::<&str>(&[], SAMPLE_ORDER_SEED).is_empty());
    }

    /// **A different seed gives a different order** of the same names.
    #[test]
    fn the_seed_sets_the_order() {
        let owned: Vec<String> = (0..30).map(|i| format!("sample_{i}")).collect();
        assert_ne!(
            sample_order(&owned, SAMPLE_ORDER_SEED),
            sample_order(&owned, SAMPLE_ORDER_SEED + 1)
        );
    }

    /// **The order is pinned**: four tomato accessions and four more names, as the hash ranks them.
    /// The order must be the same on every platform; a change of hash or of how a name reaches it —
    /// another crate version, another toolchain — would move every large cohort's subsets, and fails
    /// here first.
    #[test]
    fn the_order_of_known_names_is_pinned() {
        let names = [
            "SRR5079906",
            "SRR5079864",
            "SRR5079878",
            "SRR5079876",
            "HG002",
            "HG003",
            "HG004",
            "TS-1",
        ];
        assert_eq!(
            ordered_names(&names, SAMPLE_ORDER_SEED),
            [
                "SRR5079876",
                "SRR5079878",
                "HG004",
                "SRR5079906",
                "HG002",
                "HG003",
                "TS-1",
                "SRR5079864",
            ]
        );
    }
}
