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
//! **And it is why the wiring did not land here.** A decision that reads no posterior has
//! nothing to wait for, so discovery runs as a **pre-pass inside candidate selection** —
//! [`select_ssr`](crate::ng::calling::allele_candidates::ssr::select_ssr)'s
//! [`nominate_discovered_sequences`](crate::ng::calling::allele_candidates::ssr::nominate_discovered_sequences),
//! over the merge's rows, whose spellings are these same realigned observations — and the
//! second convergence the spec's cost argument priced in is never paid: the loop runs once,
//! over a table already widened (`doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md`
//! §6.5). This module remains the decision half, and its tests are the decision's oracle.

use std::num::NonZeroU32;

use crate::ng::calling::CandidateAlleles;
use crate::ng::calling::allele_candidates::ssr::repeat_count_of_bases;
use crate::ng::calling::likelihood::SsrSampleEvidence;
use crate::ng::locus_generation::LocusKind;
use crate::ng::run::cohort_merge::MinAltReads;

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
/// **It does not touch the candidate table**, and it does not read a posterior. This function
/// is the decision alone, so that its own tests can reach it without a loop; the shipped
/// wiring applies the same decision to the merge's rows inside candidate selection
/// ([`nominate_discovered_sequences`](crate::ng::calling::allele_candidates::ssr::nominate_discovered_sequences)),
/// where admission, the cap and the per-sample leftover already live.
///
/// # Panics
///
/// On a candidate table that is not a repeat tract's — the caller dispatches on the locus kind,
/// so reaching here with a SNP or indel locus is a routing bug and admitting a "tract sequence"
/// at one would put a length change into a table that has no motif to measure it against.
pub fn discover_tract_alleles<'s>(
    per_sample: &[SsrSampleEvidence<'_>],
    candidates: &CandidateAlleles,
    bar: MinAltReads,
    room: usize,
    scratch: &'s mut DiscoveryScratch,
) -> &'s [DiscoveredAllele] {
    let LocusKind::Ssr(detail) = candidates.kind() else {
        panic!(
            "a discovery round was handed a {:?} candidate table: the round is defined on \
             stutter attribution and has no motif to measure a length change against, so \
             reaching here off the repeat path is a routing bug",
            candidates.kind()
        );
    };
    let motif = &detail.motif;
    scratch.admitted.clear();

    for (sample, evidence) in per_sample.iter().enumerate() {
        // **The denominator is the sample's spanning reads alone.** A read that ran out
        // inside the tract says the tract is *at least* this long and cannot say which
        // length it is, so counting it here would make the share easier to clear exactly
        // where the evidence is weakest.
        let spanning: u32 = evidence
            .complete_observations()
            .map(|(_, observation)| observation.num_obs)
            .fold(0, u32::saturating_add);

        // One entry per distinct sequence this sample showed, holding where to read its
        // bases back from and how many reads carried it. Linear, because a sample's distinct
        // sequences at one tract are a handful and a map would allocate per locus.
        scratch.per_sample.clear();
        for (position, observation) in evidence.complete_observations() {
            let bases = &*observation.bases;
            match scratch
                .per_sample
                .iter_mut()
                .find(|(seen, _)| &*evidence.observations[*seen].bases == bases)
            {
                Some((_, reads)) => *reads = reads.saturating_add(observation.num_obs),
                None => scratch.per_sample.push((position, observation.num_obs)),
            }
        }

        for &(position, reads) in &scratch.per_sample {
            let bases = &*evidence.observations[position].bases;
            if candidates.iter().any(|candidate| candidate == bases) {
                continue;
            }
            if !bar.reached_by(reads, spanning) {
                continue;
            }
            // **A sequence below one whole copy of the unit has no rung and is not offered.**
            // The tract model's ladder is written in whole repeats, so admitting it would
            // reach the evidence's `NonZeroU32` count and stop the locus one step later —
            // refusing it here keeps the refusal where the reason is.
            let Some(repeats) = NonZeroU32::new(repeat_count_of_bases(bases, motif)) else {
                continue;
            };
            match scratch
                .admitted
                .iter_mut()
                .find(|held| &*held.bases == bases)
            {
                // **The best sample's count, not the cohort's sum.** The bar is per sample
                // and one sample has to clear it alone, so two samples showing a sequence
                // twice each is not four reads for this purpose.
                Some(held) => {
                    if reads > held.best_sample_reads {
                        held.best_sample_reads = reads;
                        held.best_sample = sample;
                    }
                }
                None => scratch.admitted.push(DiscoveredAllele {
                    bases: Box::from(bases),
                    repeats,
                    best_sample_reads: reads,
                    best_sample: sample,
                }),
            }
        }
    }

    // Best-supported first so a cap keeps the sequences a sample showed most, and ties broken
    // on the bases so the table is a function of the evidence rather than of the walk order.
    scratch.admitted.sort_unstable_by(|left, right| {
        right
            .best_sample_reads
            .cmp(&left.best_sample_reads)
            .then_with(|| left.bases.cmp(&right.bases))
    });
    scratch.admitted.truncate(room);
    &scratch.admitted
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::inference::DEFAULT_DISCOVERY_BAR;
    use crate::ng::locus_generation::{
        ReadWitness, SequenceObservation, SsrDetail, WitnessedLocusPositions,
    };
    use crate::ng::run::cohort_merge::{MinAltObs, MinAltReadShare};
    use crate::ng::types::{Motif, ReadGroupId, SummedLogError};

    /// **A dinucleotide**, so that one whole repeat and one base are different steps — a
    /// homopolymer fixture cannot tell a sequence one repeat long from one base long, and every
    /// off-by-a-period mistake in this module would pass on one.
    const MOTIF: &[u8] = b"AT";

    fn tract(repeats: u32) -> Vec<u8> {
        MOTIF.repeat(repeats as usize)
    }

    fn detail() -> SsrDetail {
        SsrDetail {
            motif: Motif::new(MOTIF).expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        }
    }

    /// A candidate table over the tract lengths given, the first as the reference.
    fn table(repeats: &[u32]) -> CandidateAlleles {
        let mut alleles = CandidateAlleles::new(
            tract(repeats[0]).into_boxed_slice(),
            LocusKind::Ssr(detail()),
        );
        for &length in &repeats[1..] {
            alleles.admit(tract(length).into_boxed_slice());
        }
        alleles
    }

    fn spanning(bases: &[u8], num_obs: u32) -> SequenceObservation {
        observation(bases, ReadWitness::Complete, num_obs)
    }

    /// A read that ran out inside the tract — it says the tract is *at least* this long.
    fn ran_out(bases: &[u8], num_obs: u32) -> SequenceObservation {
        let witness = ReadWitness::Partial {
            positions: WitnessedLocusPositions::one_run_from_offset_and_length(0, 3)
                .expect("a run from offset zero over three positions is a witness"),
        };
        observation(bases, witness, num_obs)
    }

    fn observation(bases: &[u8], witness: ReadWitness, num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.into(),
            read_witness: witness,
            read_group: ReadGroupId(0),
            num_obs,
            num_fwd: num_obs / 2,
            q_sum: SummedLogError::from_nats(-13.5),
            mapq_sum: 60 * u32::from(num_obs as u16),
            mapq_sum_sq: 3_600 * u64::from(num_obs),
            placed_left: 1,
            chain_ids: Vec::new(),
        }
    }

    /// The shipped bar: 2 reads **and** 15 in 100 of the sample's spanning reads.
    fn shipped_bar() -> MinAltReads {
        DEFAULT_DISCOVERY_BAR
    }

    /// A bar spelled out, so a test can put the two halves where it wants them.
    fn bar_of(floor: u32, share: f64) -> MinAltReads {
        MinAltReads {
            floor: MinAltObs(NonZeroU32::new(floor).expect("a floor of at least one read")),
            share: MinAltReadShare::new_or_panic(share),
        }
    }

    fn found<'a>(
        observations: &[Vec<SequenceObservation>],
        candidates: &CandidateAlleles,
        bar: MinAltReads,
        room: usize,
        detail: &'a SsrDetail,
        scratch: &'a mut DiscoveryScratch,
    ) -> Vec<(Vec<u8>, u32, usize)> {
        let evidence: Vec<SsrSampleEvidence<'_>> = observations
            .iter()
            .map(|rows| SsrSampleEvidence::new(rows, detail))
            .collect();
        discover_tract_alleles(&evidence, candidates, bar, room, scratch)
            .iter()
            .map(|one| (one.bases.to_vec(), one.best_sample_reads, one.best_sample))
            .collect()
    }

    /// **The case the whole mechanism exists for**, and it is the one selection cannot reach: a
    /// diploid sample nominates at most two lengths, so a third it carries is never put forward
    /// however many reads it has.
    ///
    /// Here the sample shows 12 reads at 8 repeats, 10 at 7 — the reference and its slip — and 8
    /// at 5. Selection's two peaks are 8 and 7; the 5 is invisible to it. The round finds it.
    #[test]
    fn a_third_length_one_sample_carries_is_found_where_selection_took_its_two_peaks() {
        let candidates = table(&[8, 7]);
        let rows = vec![vec![
            spanning(&tract(8), 12),
            spanning(&tract(7), 10),
            spanning(&tract(5), 8),
        ]];
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        assert_eq!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch),
            vec![(tract(5), 8, 0)],
            "the length neither peak covered, with the reads that earned it"
        );
    }

    /// A sequence the table already holds is not offered again, whatever its support.
    #[test]
    fn a_length_already_on_the_table_is_not_offered_again() {
        let candidates = table(&[8, 5]);
        let rows = vec![vec![spanning(&tract(8), 12), spanning(&tract(5), 20)]];
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        assert!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch).is_empty(),
            "both sequences are candidates already"
        );
    }

    /// **The read floor binds where the share cannot**: at 6 spanning reads, 15 in 100 is one
    /// read, so the floor of two is the only thing standing between a single stray read and a
    /// minted allele. This is the low-depth corner spec §4.1 names as the dangerous one.
    #[test]
    fn one_read_never_mints_an_allele_however_few_the_sample_has() {
        let candidates = table(&[8]);
        let rows = vec![vec![spanning(&tract(8), 5), spanning(&tract(5), 1)]];
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        assert!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch).is_empty(),
            "one read is below the floor of two, and the share of six reads asks for one"
        );
        assert_eq!(
            found(
                &[vec![spanning(&tract(8), 5), spanning(&tract(5), 2)]],
                &candidates,
                shipped_bar(),
                4,
                &detail,
                &mut scratch
            ),
            vec![(tract(5), 2, 0)],
            "and two reads clears both halves at this depth"
        );
    }

    /// **The share binds where the floor cannot**: at 40 spanning reads, 15 in 100 asks for six,
    /// so a sequence with four reads is refused although it clears the floor of two. This is the
    /// high-depth end of the same bar.
    #[test]
    fn the_share_refuses_a_sequence_the_read_floor_would_admit_at_depth() {
        let candidates = table(&[8]);
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        let rows = vec![vec![spanning(&tract(8), 36), spanning(&tract(5), 4)]];
        assert!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch).is_empty(),
            "four reads of forty is 10 in 100, below the bar's 15"
        );
        let rows = vec![vec![spanning(&tract(8), 34), spanning(&tract(5), 6)]];
        assert_eq!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch),
            vec![(tract(5), 6, 0)],
            "six of forty is 15 in 100 and clears it"
        );
    }

    /// **A read that ran out inside the tract counts in neither half of the bar.** It says the
    /// tract is *at least* this long and cannot say which length it is, so counting it in the
    /// numerator would mint an allele from evidence that names none, and counting it in the
    /// denominator would make the share easier to clear exactly where the reads are weakest.
    #[test]
    fn a_read_that_ran_out_inside_the_tract_counts_in_neither_half_of_the_bar() {
        let candidates = table(&[8]);
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();

        // Numerator: five partial reads at a length nothing spanned admit nothing.
        let rows = vec![vec![spanning(&tract(8), 10), ran_out(&tract(5), 5)]];
        assert!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch).is_empty(),
            "a partial names no length, so it cannot earn one"
        );

        // Denominator: two spanning reads of ten spanning is 20 in 100 and clears the bar —
        // and stays cleared however many partials sit beside them. Counted in, thirty
        // partials would put the share at 2 in 32 and refuse it.
        let rows = vec![vec![
            spanning(&tract(8), 8),
            spanning(&tract(5), 2),
            ran_out(&tract(9), 30),
        ]];
        assert_eq!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch),
            vec![(tract(5), 2, 0)],
            "the share is of the spanning reads alone"
        );
    }

    /// **One sample has to clear the bar alone.** Three samples showing a length once each is
    /// three reads across the cohort and one read anywhere, and one read is a stutter product.
    #[test]
    fn a_cohort_sum_does_not_clear_a_bar_no_single_sample_clears() {
        let candidates = table(&[8]);
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        let rows = vec![
            vec![spanning(&tract(8), 5), spanning(&tract(5), 1)],
            vec![spanning(&tract(8), 5), spanning(&tract(5), 1)],
            vec![spanning(&tract(8), 5), spanning(&tract(5), 1)],
        ];
        assert!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch).is_empty(),
            "three reads in three samples is one read in each"
        );
    }

    /// The count reported is the **best sample's**, and it names that sample.
    #[test]
    fn the_reported_count_is_the_best_samples_and_it_names_that_sample() {
        let candidates = table(&[8]);
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        let rows = vec![
            vec![spanning(&tract(8), 6), spanning(&tract(5), 2)],
            vec![spanning(&tract(8), 3), spanning(&tract(5), 7)],
        ];
        assert_eq!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch),
            vec![(tract(5), 7, 1)],
            "sample 1 showed it seven times; sample 0's two are not added to them"
        );
    }

    /// **The cap keeps the best-supported**, because a cap that kept the first-found would make
    /// the table depend on the order the samples happen to sit in.
    #[test]
    fn the_cap_keeps_the_best_supported_and_cuts_the_rest() {
        let candidates = table(&[8]);
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        let rows = vec![vec![
            spanning(&tract(8), 4),
            spanning(&tract(3), 5),
            spanning(&tract(5), 9),
            spanning(&tract(7), 7),
        ]];
        assert_eq!(
            found(
                &rows,
                &candidates,
                bar_of(2, 0.05),
                2,
                &detail,
                &mut scratch
            ),
            vec![(tract(5), 9, 0), (tract(7), 7, 0)],
            "nine and seven survive a cap of two; five does not"
        );
    }

    /// A sequence below one whole copy of the unit has no rung on the stutter ladder, so it is
    /// not offered — the refusal sits here, where the reason is, rather than at the conversion
    /// that would fail one step later.
    #[test]
    fn a_sequence_shorter_than_one_whole_repeat_is_not_offered() {
        let candidates = table(&[8]);
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        let rows = vec![vec![spanning(&tract(8), 10), spanning(b"A", 6)]];
        assert!(
            found(
                &rows,
                &candidates,
                bar_of(2, 0.05),
                4,
                &detail,
                &mut scratch
            )
            .is_empty(),
            "one base of a two-base unit floors to zero repeats and has no rung"
        );
    }

    /// **Ties break on the bases, so the table is a function of the evidence.** Two sequences
    /// with equal support must come back in one order however the samples were walked.
    #[test]
    fn two_sequences_with_equal_support_come_back_in_one_order() {
        let candidates = table(&[8]);
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        let forwards = vec![vec![
            spanning(&tract(8), 10),
            spanning(&tract(5), 4),
            spanning(&tract(6), 4),
        ]];
        let backwards = vec![vec![
            spanning(&tract(8), 10),
            spanning(&tract(6), 4),
            spanning(&tract(5), 4),
        ]];
        let one = found(
            &forwards,
            &candidates,
            bar_of(2, 0.05),
            4,
            &detail,
            &mut scratch,
        );
        let other = found(
            &backwards,
            &candidates,
            bar_of(2, 0.05),
            4,
            &detail,
            &mut scratch,
        );
        assert_eq!(one, other, "the walk order must not reach the answer");
        assert_eq!(
            one.first().map(|(bases, _, _)| bases.clone()),
            Some(tract(5)),
            "the shorter spelling sorts first"
        );
    }

    /// **A second round over the same evidence admits nothing**, which is this module's own
    /// consequence and the reason the round cap cannot bind: the eligible set is a function of
    /// the observations and the table alone, so once the first round's finds are on the table
    /// there is nothing left to find.
    #[test]
    fn a_second_round_over_the_same_evidence_finds_nothing() {
        let mut candidates = table(&[8, 7]);
        let rows = vec![vec![
            spanning(&tract(8), 12),
            spanning(&tract(7), 10),
            spanning(&tract(5), 8),
        ]];
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        let first = found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch);
        assert_eq!(first.len(), 1, "the first round finds the hidden length");
        for (bases, _, _) in &first {
            candidates.admit(bases.clone().into_boxed_slice());
        }
        assert!(
            found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch).is_empty(),
            "and the second finds nothing, so the rounds stop at two"
        );
    }

    /// A discovery round at a SNP or indel locus is a routing bug, and a crash is the right
    /// answer: the alternative is measuring a length change against a motif that does not exist.
    #[test]
    #[should_panic(expected = "a discovery round was handed a")]
    fn a_generic_candidate_table_is_refused() {
        let candidates = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
        let detail = detail();
        let mut scratch = DiscoveryScratch::new();
        let rows = vec![vec![spanning(b"AT", 10)]];
        let _ = found(&rows, &candidates, shipped_bar(), 4, &detail, &mut scratch);
    }
}
