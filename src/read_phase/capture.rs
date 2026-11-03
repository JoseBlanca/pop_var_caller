//! **One sample's read alleles at one cohort locus** — which allele each read showed, captured in
//! the merge where the reads' chain ids and the composed alleles meet, and nowhere later
//! (read phase spec §2.2).

use crate::splitmix::splitmix64;
use crate::types::ChainId;

/// **The most read alleles one sample keeps at one record** (§2.4), a starting value adopted by the
/// owner on 2026-11-03. It bounds what a deep sample costs from the merge to the phase writer
/// whatever its depth; WhatsHap caps at 15 fragments a site for the same reason.
pub const READ_ALLELES_KEPT: usize = 32;

/// The largest error Phred a read allele stores. A geometric mean of read errors can be far below
/// any one read's (§2.2), and the link test floors every error at 0.01 anyway (§2.5), so nothing
/// above this changes a weight.
pub const LARGEST_ERROR_PHRED: u8 = 60;

/// **One read's allele at one locus**: the read, the allele it showed and how likely that is a
/// sequencing or mapping error.
///
/// Ordered by read first, which is the order a sample's list is kept in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReadAllele {
    /// The read's chain id. **Meaningful only within its own sample**: two samples' ids are two
    /// numberings (§6, first trap).
    pub read: ChainId,
    /// The allele, by its index in whichever table the list is currently numbered against: the
    /// merge's allele table at capture, the record's alleles once calling has written one.
    pub allele: u32,
    /// The read's error as a Phred integer, `0..=LARGEST_ERROR_PHRED`: the mean natural-log error
    /// of the reads behind its observation, converted and rounded once, here, so every later step
    /// works on the same integer (§2.2).
    pub error_phred: u8,
}

/// **A mean natural-log read error as the stored Phred byte**: `−10·log₁₀ e`, rounded, at least 0
/// and at most [`LARGEST_ERROR_PHRED`].
///
/// The mean is of logarithms, so it is a *geometric* mean of the errors, which sits well below the
/// arithmetic one (§2.2); the link test's error floor is what answers for that. A mean that is not
/// a number reads as Phred 0, an error of one, which the link test weighs at exactly nothing.
#[must_use]
pub fn error_phred_of_mean_ln_error(mean_ln_error: f64) -> u8 {
    let phred = -mean_ln_error * (10.0 / std::f64::consts::LN_10);
    if phred.is_nan() || phred <= 0.0 {
        return 0;
    }
    phred.round().min(f64::from(LARGEST_ERROR_PHRED)) as u8
}

/// **What decides which reads a deep sample keeps**: the read's rank in its walked region, hashed
/// (§2.4). The rank, `read − chain_id_base`, is the same whether the region was walked alone or
/// after other ground, where the raw id is not; `splitmix64` is the hash the psp read cap already
/// thins by, and it is a bijection, so two reads never share a key.
#[must_use]
pub fn cap_key(read: ChainId, chain_id_base: u64) -> u64 {
    splitmix64(read.wrapping_sub(chain_id_base))
}

/// **Settle one sample's captured list**, in place: sorted by read, one entry a read, at most
/// [`READ_ALLELES_KEPT`] of them (§2.2).
///
/// - **A read listed twice with one allele keeps one entry**, the one with the lower Phred, so one
///   read never counts twice.
/// - **A read listed with two alleles is removed.** That is a pair whose mates overlap and
///   disagree, which the merge emits as two sequences of one minted record under one id; like
///   HapCUT2, a molecule that showed two things at a site says nothing there.
/// - **Above the cap, the reads with the smallest [`cap_key`] are kept**, and the list is put back
///   in read order.
///
/// `chain_id_base` is the walked region's first id, shared by every minted record of one sample
/// inside one cohort locus: a locus never crosses a segment boundary.
pub fn settle(read_alleles: &mut Vec<ReadAllele>, chain_id_base: u64) {
    read_alleles.sort_unstable();
    read_alleles.dedup_by(|later, kept| later.read == kept.read && later.allele == kept.allele);
    // After that, a read still listed twice showed two alleles.
    let mut kept = 0;
    let mut at = 0;
    while at < read_alleles.len() {
        let read = read_alleles[at].read;
        let mut end = at + 1;
        while end < read_alleles.len() && read_alleles[end].read == read {
            end += 1;
        }
        if end - at == 1 {
            read_alleles[kept] = read_alleles[at];
            kept += 1;
        }
        at = end;
    }
    read_alleles.truncate(kept);
    if read_alleles.len() > READ_ALLELES_KEPT {
        read_alleles.select_nth_unstable_by_key(READ_ALLELES_KEPT - 1, |read_allele| {
            cap_key(read_allele.read, chain_id_base)
        });
        read_alleles.truncate(READ_ALLELES_KEPT);
        read_alleles.sort_unstable();
    }
}

#[cfg(test)]
mod tests;
