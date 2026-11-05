//! **The splitmix64 finaliser**, shared by the two places that must thin reads by the same key: the
//! psp read cap (`run::psp_source::read_cap`) and the read-phase cap (`read_phase::capture`).

/// The splitmix64 finaliser: a bijection that spreads consecutive identifiers across the whole
/// range, so "lowest hash" is not "lowest identifier" — which would keep the reads that arrived
/// first, a sample tilted towards early alignment starts.
#[must_use]
pub fn splitmix64(value: u64) -> u64 {
    let mut z = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}
