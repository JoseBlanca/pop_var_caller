//! **Tract lengths the model is explaining as slippage, which some sample showed too often for
//! that to be true** — the discovery round of `doc/devel/ng/spec/calling_em_loop.md` §4.1.
//!
//! # What the round is for
//!
//! A true allele can hide *under* stutter. Every read carrying it is booked as a slip product of
//! a called length, so its repeat count never surfaces as a candidate and no sample can be
//! genotyped for it. Selection cannot find it by construction: **a sample nominates at most
//! `ploidy` lengths** — its best-supported ones — so a third length in that sample is not put
//! forward however many reads it has
//! ([`promote_rungs_for_sample`](crate::ng::calling::allele_candidates::ssr::promote_rungs_for_sample)
//! truncates to `ploidy`). A heterozygote carrying a long allele and a short one is exactly the
//! case: the long allele's own contraction slips can outnumber the short allele's reads, so the
//! sample's two peaks are the long allele and its slip, and the short allele is third.
//!
//! # ng's retrace is the evidence it already holds
//!
//! HipSTR retraces each read's maximum-likelihood alignment and, where the trace says the read
//! slipped, counts the tract sequence the trace implies (spec §4.1). **ng has that answer
//! already**: the tract locus generator realigns every read against the tract before the caller
//! sees it, so a [`SequenceObservation`](crate::ng::locus_generation::SequenceObservation)'s
//! `bases` *are* the implied tract sequence. So the retrace is a walk over observations rather
//! than a second alignment, and what is left to decide is which of them the candidate set does
//! not already hold.
//!
//! **An observation whose bases are no candidate's is, by construction, what the model is
//! explaining as slippage or as junk** — those are the only two terms in the tract row that can
//! account for a read matching no candidate (`doc/devel/ng/spec/read_likelihoods.md` §4.5). That
//! is what makes the eligible set exactly "the sequences some sample showed that are not
//! candidates", with no posterior to consult.
//!
//! # The bar, and which reads it is a share of
//!
//! HipSTR's, inherited and soft: a sequence needs **at least
//! [`DEFAULT_DISCOVERY_MIN_READS`](super::DEFAULT_DISCOVERY_MIN_READS) reads and at least
//! [`DEFAULT_DISCOVERY_MIN_SPANNING_READ_SHARE`](super::DEFAULT_DISCOVERY_MIN_SPANNING_READ_SHARE)
//! of one sample's tract-spanning reads**, and one sample clearing it admits the sequence for the
//! cohort. **Spanning reads, not all reads**: a read that ran out inside the tract says the tract
//! is *at least* this long and cannot say which length it is, so counting it in the denominator
//! would make the share easier to clear exactly where the evidence is weakest.
//!
//! **The two halves bind at opposite ends of the depth range and that is why they are both here**
//! (spec §4.1): below about 13 reads a position, 2 reads already clears 15%, so the count is the
//! only constraint — and two reads is what a single stutter product looks like. Above it the
//! share takes over. A run at three reads a position is therefore admitting on "2 of this
//! sample's 3 reads", which is the corner this mechanism is most dangerous in, and it is why the
//! numbers are soft and swept rather than fixed.
//!
//! # Why a round can only fire once on ng's evidence
//!
//! **This module's own consequence, and it is worth stating because the spec's cost argument
//! assumes otherwise.** The eligible set is a function of the observations and the candidate
//! table alone — no posterior enters it. So a second round over the same evidence, with the
//! first round's admissions now in the table, finds the same sequences and admits none of them.
//! Rounds therefore stop at two: one that admits and one that establishes there is nothing left.
//! The round cap ([`DEFAULT_DISCOVERY_MAX_ROUNDS`](super::DEFAULT_DISCOVERY_MAX_ROUNDS)) is a
//! guard that cannot bind, and the "locus that keeps finding one more allele" §4.1 worries about
//! cannot arise here. That is a property of ng's retrace being structural, not a claim about
//! HipSTR's.
//!
//! **What that does *not* make free is the round's cost**, which is the second convergence and
//! is paid at every tract whether or not anything is found.

use std::num::NonZeroU32;

use crate::ng::calling::CandidateAlleles;
use crate::ng::calling::allele_candidates::ssr::repeat_count_of_bases;
use crate::ng::calling::likelihood::SsrSampleEvidence;
use crate::ng::run::cohort_merge::MinAltReads;
use crate::ng::types::Motif;

/// **One tract sequence a discovery round would add, and what earned it.**
///
/// The bases rather than the repeat count, because ng's candidate table is a table of sequences:
/// two spellings of one length are two alleles and the tract prior lays its mass over the rung
/// they share (`doc/devel/ng/spec/calling_priors.md` §5.1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscoveredAllele {
    /// The tract bases, exactly as some sample's reads showed them.
    pub bases: Box<[u8]>,
    /// Whole motif copies the sequence carries — **non-zero**, because the stutter ladder is
    /// written in whole repeats and has no rung below its first. A sequence shorter than one
    /// copy of the unit is not offered at all rather than offered and refused downstream.
    pub repeats: NonZeroU32,
    /// The most reads any one sample showed it with — **the count that cleared the bar**, not
    /// the cohort's total. Two samples showing it twice each is not four reads for this purpose:
    /// the bar is per sample and one sample has to clear it alone.
    pub best_sample_reads: u32,
    /// Which sample that was, in the run's sample order — so a report can say who found it.
    pub best_sample: usize,
}

/// Reusable buffers for a discovery round, so a round allocates nothing per locus.
///
/// **Two tallies and not one.** The per-sample one is rebuilt for each sample and holds that
/// sample's distinct sequences; the cohort one accumulates across samples and is what the round
/// returns. Merging them would lose the per-sample denominator the bar is a share of.
#[derive(Debug, Default)]
pub struct DiscoveryScratch {
    /// Per sample: for each distinct sequence among its complete observations, the index of the
    /// first observation showing it and how many reads showed it in total.
    ///
    /// **Indices rather than bases**, so a sample's tally borrows nothing and copies nothing —
    /// the bases are read back through the evidence when a sequence clears the bar.
    per_sample: Vec<(usize, u32)>,
    /// The sequences admitted so far, across every sample walked.
    admitted: Vec<DiscoveredAllele>,
}

impl DiscoveryScratch {
    /// A round's buffers, empty.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

/// **What a discovery round would admit at one repeat tract**, best-supported first.
///
/// Walks each sample's complete observations, tallies the sequences the candidate table does not
/// hold, and keeps those a sample showed with enough reads to clear `bar`. The result is capped
/// at `room`, which is how many more alleles the locus's table may carry.
///
/// **The order is the order they are admitted in and it is deliberate**: best-supported first,
/// so that a cap cutting the list keeps the sequences one sample showed most. Ties break on the
/// bases, so two sequences with equal support cannot swap between runs — a tract's candidate
/// table has to be a function of its evidence and nothing else (spec §8's determinism).
///
/// # What it does not do
///
/// **It does not touch the candidate table**, and it does not read a posterior. Growing the
/// table, rebuilding what a wider table changes and re-running the loop are the caller's, in
/// [`summarise_condition`](super::summarise_condition); this function is the decision alone so
/// that its own tests can reach it without a loop.
///
/// # Panics
///
/// On a candidate table that is not a repeat tract's — the caller dispatches on the locus kind,
/// so reaching here with a SNP or indel locus is a routing bug and admitting a "tract sequence"
/// at one would put a length change into a table that has no motif to measure it against.
pub fn discover_tract_alleles<'a>(
    per_sample: &[SsrSampleEvidence<'a>],
    candidates: &CandidateAlleles,
    motif: &Motif,
    bar: MinAltReads,
    room: usize,
    scratch: &mut DiscoveryScratch,
) -> &'a [DiscoveredAllele] {
    let _ = per_sample;
    let _ = candidates;
    let _ = motif;
    let _ = bar;
    let _ = room;
    let _ = scratch;
    unimplemented!("filled in below")
}
