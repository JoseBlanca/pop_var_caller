//! **ng's, not a copy** — the one number that decides which reads survive a cap.
//!
//! Two places in the walk have to drop reads and cannot drop all of them: the
//! per-position depth cap ([`WalkerConfig::max_snp_column_depth`] /
//! `max_indel_column_depth`), and the ceiling on how many reads the walk holds open
//! at once (`PileupGeneratorConfig::max_active_reads`). Both used to answer the
//! question the same way — *keep whichever came first* — and that answer has two
//! faults the owner ruled out:
//!
//! - **It depends on something that is not the data.** The per-position cap kept a
//!   prefix of the active set's storage order, which `swap_remove` scrambles; a change
//!   to how the set stores reads moved 88,351 of 341,094 emitted rows on a 300× walk
//!   without any read changing. The evidence at a capped position is supposed to be a
//!   fact about the reads, not about the container.
//! - **It is biased towards early alignment starts.** Reads arrive sorted by start, so
//!   "first come" and "starts leftmost" are the same rule wherever arrival order is
//!   preserved. A subsample chosen that way over-represents reads whose left end is
//!   near the position and under-represents reads that reach it from the left, which
//!   is a systematic tilt in `placed_left` and in witness extent.
//!
//! # The rule
//!
//! Every read gets a 64-bit **sampling key** derived from its identity alone —
//! its query name and which mate of the pair it is — and a cap keeps the reads with
//! the **smallest** keys. Nothing else enters: not the alignment start, not the order
//! of arrival, not the read's position in any container.
//!
//! Three properties follow, and the third only approximately:
//!
//! 1. **Deterministic across runs and across processes.** The hash below is written
//!    out here rather than taken from `ahash` or `DefaultHasher` precisely because
//!    those are free to seed themselves per process or to change between compiler
//!    versions. `parity::ng_emits_the_same_bytes_in_a_second_process` is the test that
//!    would catch a regression, and it should never have to.
//! 2. **Unbiased with respect to alignment start** — exactly, not approximately. The
//!    key is a function of the query name, which carries no positional information, so
//!    the kept subsample is a uniform random subset of the covering reads with respect
//!    to every alignment property.
//! 3. **Stable between adjacent positions — approximately.** Because the kept set at a
//!    position is *the k smallest keys among the reads covering it*, two neighbouring
//!    positions covered by the same reads keep the same reads. The set changes only
//!    when a read enters or leaves the covering set, and then by one read at a time:
//!    a read that ends drops out and the next-smallest key takes its slot. This is the
//!    bottom-k sketch's stability property. What it does **not** give is stability
//!    across a change in the cap itself, or across the ceiling evicting a read that a
//!    later position would have wanted — those are the approximations, and the second
//!    only bites in a region deep enough to fill the walk.
//!
//! # Why the query name and not the chain id
//!
//! The chain id is already on every contributor and would need no lookup, and it is
//! deterministic. It is still the wrong input: it is minted in **admission order**, so
//! it is a fact about the walk's history — how many reads preceded this one, which
//! regions were walked, whether an earlier read was shed — rather than about the read.
//! Hashing it would hide the ordering dependence rather than remove it. The query name
//! is what a BAM says the read *is*.
//!
//! [`WalkerConfig::max_snp_column_depth`]: crate::pileup::walker::WalkerConfig::max_snp_column_depth

use super::{MateRole, PreparedRead};

/// FNV-1a's 64-bit offset basis and prime. Spelled out because this hash is a
/// **format**: it decides which reads reach a `.psp`, so it may not change silently
/// with a dependency bump.
const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

/// The read's sampling key: **smaller survives a cap**.
///
/// FNV-1a over the query name, one byte of mate role folded in so the two mates of a
/// pair are two draws rather than one, then a splitmix64 finaliser. FNV-1a alone
/// avalanches poorly in its high bits, and "smallest wins" reads the high bits first —
/// so without the finaliser the rule would degenerate towards "whichever name starts
/// with the earliest byte", which is a property of the sequencer's naming scheme.
pub(super) fn sampling_key(read: &PreparedRead) -> u64 {
    let mut hash = FNV_OFFSET_BASIS;
    for &byte in read.qname.as_bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash ^= match read.mate_role {
        MateRole::Solo => 0,
        MateRole::FirstOfPair => 1,
        MateRole::SecondOfPair => 2,
    };
    hash = hash.wrapping_mul(FNV_PRIME);
    splitmix64(hash)
}

/// The finaliser from splitmix64 — a bijection, so it cannot create collisions FNV-1a
/// did not already have, and it spreads every input bit across all 64 output bits.
fn splitmix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::ng::types::ReadGroupId;
    use crate::pileup::walker::CigarOp;

    fn read(qname: &str, alignment_start: u32, mate_role: MateRole) -> PreparedRead {
        PreparedRead {
            chrom_id: 0,
            alignment_start,
            alignment_end: alignment_start + 9,
            cigar: vec![CigarOp::Match(10)],
            seq: vec![b'A'; 10],
            bq_baq: vec![30; 10],
            mq_log_err: -3.0,
            mapq: 60,
            is_reverse_strand: false,
            qname: Arc::from(qname),
            mate_role,
            adaptor_boundary: None,
            read_group: ReadGroupId(0),
        }
    }

    /// **The key is the read's identity and nothing else.** Two reads with the same
    /// name and mate role placed a megabase apart draw the same number, so no cap can
    /// prefer one for where it landed.
    #[test]
    fn the_key_ignores_where_the_read_aligned() {
        assert_eq!(
            sampling_key(&read("q", 1, MateRole::Solo)),
            sampling_key(&read("q", 1_000_001, MateRole::Solo)),
        );
    }

    /// The two mates of one fragment are two draws. Sharing a key would make a cap keep
    /// or drop a pair together, which is a correlation the sampling does not want and
    /// nothing asks for.
    #[test]
    fn the_two_mates_of_a_pair_draw_separately() {
        let first = sampling_key(&read("q", 1, MateRole::FirstOfPair));
        let second = sampling_key(&read("q", 1, MateRole::SecondOfPair));
        let solo = sampling_key(&read("q", 1, MateRole::Solo));
        assert_ne!(first, second);
        assert_ne!(first, solo);
        assert_ne!(second, solo);
    }

    /// **The literal values, pinned.** This hash decides which reads reach a `.psp`,
    /// so a change to it is a change to the output on every capped position — it has to
    /// arrive as a failing test, not as a quiet re-baseline. Recomputing these three
    /// numbers is the price of changing the rule deliberately.
    #[test]
    fn the_hash_is_pinned_to_these_values() {
        assert_eq!(
            sampling_key(&read("", 1, MateRole::Solo)),
            2_737_183_428_366_584_608
        );
        assert_eq!(
            sampling_key(&read("read1", 1, MateRole::FirstOfPair)),
            15_641_639_245_897_738_950
        );
        assert_eq!(
            sampling_key(&read("read1", 1, MateRole::SecondOfPair)),
            2_625_164_446_083_436_970
        );
    }

    /// The rule is only unbiased if the keys spread. Ten thousand names in the shape a
    /// sequencer emits, split at the median of the whole draw: a hash that leaked the
    /// counter into the high bits would put the low names on one side.
    #[test]
    fn keys_spread_across_the_range_for_sequencer_style_names() {
        let mut keys: Vec<u64> = (0..10_000)
            .map(|index| {
                sampling_key(&read(
                    &format!("A00123:45:HXXXX:1:1101:{index}:1000"),
                    1,
                    MateRole::FirstOfPair,
                ))
            })
            .collect();
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(
            keys.len(),
            10_000,
            "the keys must not collide at this scale"
        );
        // The smallest tenth of the key space should hold about a tenth of the draws.
        let tenth = u64::MAX / 10;
        let below = keys.iter().filter(|&&key| key < tenth).count();
        assert!(
            (800..=1200).contains(&below),
            "{below} of 10,000 keys fell in the smallest tenth of the range; a uniform \
             hash gives about 1,000"
        );
    }
}
