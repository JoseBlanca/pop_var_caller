//! One STR locus, reduced to the entry a stratum's table is keyed on.
//!
//! This is the only place that decides what a read's repeat count *is*, which the locus
//! type cannot answer: a locus carries its reference bases, its motif and the sequence each
//! read showed, and turning those into "this read sits two copies short of the reference"
//! is a division that has to be done somewhere and done once.
//!
//! Its own file because the shaping of data and the mathematics on it never live together
//! (`arch/parameter_prepass_ssr.md` §Module home) — [`super::slippage`] is the mathematics.
//!
//! The reductions land one at a time through Milestone C; this file's first is the one that
//! decides which stratum a locus belongs to at all.

use crate::ng::locus_generation::{LocusKind, SampleLocusObservations};
use crate::ng::types::{Bp, SsrPeriod};

use super::{RepeatCount, Stratum};

/// Which stratum a locus belongs to, or why none does.
///
/// **Three answers and not two**, where the architecture sketches an `Option`: the two ways a
/// locus can fail to have a stratum mean opposite things to the accumulator that asks. A locus
/// that is not one repeat tract is the SNP/indel path's business, so it is passed over in
/// silence; a repeat tract whose reference length is *not a whole number of motif copies* is a
/// delimiting fault in whatever admitted it, so it is counted and reported. The counter it
/// feeds — `SsrAccumulationCounts::loci_without_whole_repeat_reference`, which
/// `arch/parameter_prepass_ssr.md` §4 says should be near zero — is named after this variant
/// for that reason. An `Option` would collapse the two, and the counter would then have to
/// re-derive from the locus what this function had already decided.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum LocusStratum {
    /// The stratum this locus's **reference** tract puts it in.
    Stratified(Stratum),
    /// Not one delimited repeat tract, so this path does not fit it: a SNP/indel candidate, or
    /// a repeat cluster with no clean flanks. Nothing to record — the STR accumulator passes it
    /// over and the generic one has its own.
    NotOneRepeatTract,
    /// A repeat tract whose reference length is not a whole number of motif copies, so no
    /// stratum holds it: 13 bases of a 3-base motif is four copies and a base, and rounding
    /// either way would file it with tracts it is not comparable to.
    ///
    /// **Counted and skipped, never rounded**, and a run where this is not near zero is a bug
    /// report against the classification that admitted the locus rather than something this
    /// unit absorbs. Named to match the counter it feeds.
    WithoutWholeRepeatReference {
        /// How long the reference tract is.
        tract_len: Bp,
        /// The motif period it failed to divide by.
        period: SsrPeriod,
    },
}

/// Which stratum a locus belongs to — its reference tract's motif period and how many whole
/// copies of that motif the tract holds.
///
/// **A pure function of the reference**, which is the property the whole design rests on: every
/// sample files a locus under the same stratum whatever its own reads showed, so one sample's
/// stutter can be compared with another's and a cohort can pool them
/// (`spec/parameter_prepass_ssr.md` §4.1). A sample whose alleles differ from the reference does
/// not move between strata; its reads land at an offset instead.
///
/// Both halves come from the locus itself — `reference_bases` and the motif the STR generator
/// put in the locus's own detail — so nothing here re-reads the reference or the catalog.
///
/// **A tract of no bases is filed at zero copies rather than reported.** Zero divides by every
/// period, and zero copies is a stratum such tracts could be compared inside, where a
/// fractional count is comparable with nothing. Region typing emits none — the catalog's
/// loosest floor is three copies — so this is an answer to arithmetic rather than to data.
///
/// # Panics
///
/// If a tract holds more than `u32::MAX` whole copies. Unreachable twice over: the reference
/// bases are held in memory, so that would be a 4.29-billion-base allocation, and the catalog
/// refuses to serve any tract past `CATALOG_MAX_STR_LEN_BP` — 500 bases, against step 3's
/// calling default of 100 — so a tract that arrives here holds at most a few hundred copies.
/// It panics rather than saturating because saturating would file the locus under a stratum it
/// does not belong to and say nothing, which is the one thing this function is written not to
/// do.
#[must_use]
pub fn stratum_of(locus: &SampleLocusObservations) -> LocusStratum {
    // Matched kind by kind rather than through a wildcard, so that a locus kind added later
    // is a compile error here and its routing is a decision rather than a default: the one
    // outcome this function leaves no record of is the locus it passes over.
    let detail = match &locus.kind {
        LocusKind::Ssr(detail) => detail,
        LocusKind::Generic | LocusKind::SsrBundle => return LocusStratum::NotOneRepeatTract,
    };

    let period = detail.motif.ssr_period();
    let tract_bases = locus.reference_bases.len();
    let motif_bases = usize::from(period.get());

    if !tract_bases.is_multiple_of(motif_bases) {
        return LocusStratum::WithoutWholeRepeatReference {
            tract_len: Bp(tract_bases as u64),
            period,
        };
    }

    let repeats = u32::try_from(tract_bases / motif_bases)
        .expect("a tract holds at most a few hundred copies; see this function's panic note");
    LocusStratum::Stratified(Stratum::new(period, RepeatCount(repeats)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{ReadWitness, SequenceObservation, SsrDetail};
    use crate::ng::types::{ContigId, GenomeRegion, MAX_MOTIF_LEN, Motif, Position, ReadGroupId};

    /// A locus over `reference_bases`, tiled by `motif`, carrying whatever reads a test wants.
    fn ssr_locus(
        reference_bases: &[u8],
        motif: &[u8],
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1_000),
                end: Position(1_000 + reference_bases.len().max(1) as u64 - 1),
            },
            reference_bases: Box::from(reference_bases),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Ssr(SsrDetail {
                motif: Motif::new(motif).expect("a motif inside the STR period range"),
                left_flank: Box::from(&b"CCCGGG"[..]),
                right_flank: Box::from(&b"TTTAAA"[..]),
            }),
        }
    }

    /// One read's worth of observation, at whatever length it showed.
    fn observation(bases: &[u8], num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: Box::from(bases),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs,
            num_fwd: num_obs / 2,
            q_sum: -10.0,
            mapq_sum: 60 * num_obs,
            mapq_sum_sq: 3_600 * u64::from(num_obs),
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn stratum(period: u8, repeats: u32) -> Stratum {
        Stratum::new(
            SsrPeriod::try_new(usize::from(period)).expect("a period inside the STR scope"),
            RepeatCount(repeats),
        )
    }

    /// The tract's length over its motif's period, which is the whole of the arithmetic: twenty
    /// bases of a dinucleotide is a stratum of ten repeats.
    #[test]
    fn a_tract_is_stratified_by_its_reference_length_over_its_motif_period() {
        let locus = ssr_locus(b"ATATATATATATATATATAT", b"AT", Vec::new());

        assert_eq!(stratum_of(&locus), LocusStratum::Stratified(stratum(2, 10)));
    }

    /// **The stratum is a property of the reference alone**, which is what lets a cohort compare
    /// one sample's stutter with another's. Four samples that saw completely different things at
    /// the same tract — every read at the reference length, every read two copies short, every
    /// read two copies long, and nothing at all — file it under the same stratum.
    ///
    /// **The expansion is the sample a hand-written fixture forgets**, and the one this path's
    /// own measurements make rarest: a read is about five times as likely to have lost a copy as
    /// gained one, so a length rule reading the *longest* thing it can see would pass every
    /// other fixture here. A sample carrying an expansion that migrated to another stratum would
    /// be compared against nothing.
    #[test]
    fn every_sample_files_the_same_tract_under_the_same_stratum_whatever_its_reads_showed() {
        let reference = b"ATATATATATATATATATAT";
        let at_reference = ssr_locus(reference, b"AT", vec![observation(reference, 9)]);
        let two_copies_short =
            ssr_locus(reference, b"AT", vec![observation(b"ATATATATATATATAT", 9)]);
        let two_copies_long = ssr_locus(
            reference,
            b"AT",
            vec![observation(b"ATATATATATATATATATATATAT", 9)],
        );
        let saw_nothing = ssr_locus(reference, b"AT", Vec::new());

        let expected = LocusStratum::Stratified(stratum(2, 10));
        assert_eq!(stratum_of(&at_reference), expected);
        assert_eq!(stratum_of(&two_copies_short), expected);
        assert_eq!(stratum_of(&two_copies_long), expected);
        assert_eq!(stratum_of(&saw_nothing), expected);
    }

    /// **The divisor is the motif's own period at every period the STR scope admits**, not only
    /// at the mononucleotides and dinucleotides the other fixtures use. Periods 1 to 6 are all
    /// ordinary loci, and a rule that clamped the divisor — at 2, say — would file a
    /// sixteen-base tetranucleotide tract at eight copies instead of four, pooling it with a
    /// stratum whose reads are measured from a different origin, with nothing to say so.
    #[test]
    fn a_tract_is_stratified_by_its_own_period_at_every_period_the_str_scope_admits() {
        for (bases, motif, period, repeats) in [
            (&b"AAAAAAAA"[..], &b"A"[..], 1u8, 8u32),
            (&b"ATATATATATATATATATAT"[..], &b"AT"[..], 2, 10),
            (&b"CAGCAGCAGCAGCAG"[..], &b"CAG"[..], 3, 5),
            (&b"AAAGAAAGAAAGAAAG"[..], &b"AAAG"[..], 4, 4),
            (&b"AAAAGAAAAGAAAAG"[..], &b"AAAAG"[..], 5, 3),
            (&b"AAAAAGAAAAAGAAAAAG"[..], &b"AAAAAG"[..], 6, 3),
        ] {
            assert_eq!(
                stratum_of(&ssr_locus(bases, motif, Vec::new())),
                LocusStratum::Stratified(stratum(period, repeats)),
                "{} bases of a {period}-base motif",
                bases.len()
            );
        }
    }

    /// A tract that is exactly its own motif is one copy — not zero, and not "not a repeat".
    /// The catalog's floors emit none this short, but the arithmetic answers for them, and the
    /// mononucleotide loop below reaches one copy only at period 1, where a tract's length and
    /// its repeat count are the same number and so cannot separate an off-by-one from a
    /// division.
    #[test]
    fn a_tract_of_exactly_one_motif_copy_is_a_stratum_of_one_repeat() {
        for (motif, period) in [
            (&b"AT"[..], 2u8),
            (&b"CAG"[..], 3),
            (&b"AAAG"[..], 4),
            (&b"AAAAG"[..], 5),
            (&b"AAAAAG"[..], 6),
        ] {
            assert_eq!(
                stratum_of(&ssr_locus(motif, motif, Vec::new())),
                LocusStratum::Stratified(stratum(period, 1)),
                "one copy of a {period}-base motif"
            );
        }
    }

    /// **A tract that is not a whole number of copies is counted and skipped, not rounded.**
    /// Thirteen bases of a three-base motif is four copies and a base; filing it at four would
    /// pool it with clean four-copy tracts, and at five with clean five-copy ones, and its
    /// reads' offsets would be measured from a length no allele has.
    #[test]
    fn a_tract_whose_length_is_not_a_whole_number_of_copies_is_reported_rather_than_rounded() {
        let locus = ssr_locus(b"CAGCAGCAGCAGC", b"CAG", Vec::new());

        assert_eq!(
            stratum_of(&locus),
            LocusStratum::WithoutWholeRepeatReference {
                tract_len: Bp(13),
                period: SsrPeriod::try_new(3).expect("a trinucleotide"),
            },
            "13 bases of a 3-base motif is four copies and a base"
        );
    }

    /// The report names the two numbers a reader needs to find the fault upstream: how long the
    /// tract was, and what it failed to divide by. A count alone would say a locus was skipped
    /// and nothing about which.
    ///
    /// **Every remainder, not only one base over.** A rule that reported a tract only when it
    /// sat one base past a whole number of copies would floor all the others silently — and the
    /// counter that is supposed to catch exactly that would read zero while it happened.
    #[test]
    fn the_report_names_the_tract_length_and_the_period_it_did_not_divide_by() {
        for (bases, motif, tract_len, period) in [
            (&b"CAGCAGCAGCAGC"[..], &b"CAG"[..], 13, 3),
            (&b"ATATATATA"[..], &b"AT"[..], 9, 2),
            (&b"AAAGAAAGAAAGA"[..], &b"AAAG"[..], 13, 4),
            (&b"CAGCAGCAGCAGCA"[..], &b"CAG"[..], 14, 3),
            (&b"AAAGAAAGAAAGAA"[..], &b"AAAG"[..], 14, 4),
            (&b"AAAAAGAAAAAGAAAA"[..], &b"AAAAAG"[..], 16, 6),
            (&b"AAAAAGAAAAAGAAAAA"[..], &b"AAAAAG"[..], 17, 6),
        ] {
            assert_eq!(
                stratum_of(&ssr_locus(bases, motif, Vec::new())),
                LocusStratum::WithoutWholeRepeatReference {
                    tract_len: Bp(tract_len),
                    period: SsrPeriod::try_new(period).expect("a period inside the STR scope"),
                }
            );
        }
    }

    /// **At period 1 every length divides**, so a mononucleotide tract can never be reported as
    /// a non-whole number of copies — which matters because mononucleotides are the bulk of what
    /// this path fits: the copy floors admit them at 8 copies where a dinucleotide needs 6.
    ///
    /// What this cannot check, and what the every-period test above exists for: at period 1 a
    /// tract's length and its repeat count are the same number, so nothing here separates
    /// dividing by the period from not dividing at all.
    #[test]
    fn a_mononucleotide_tract_is_always_a_whole_number_of_copies() {
        for length in 1..40usize {
            let locus = ssr_locus(&b"A".repeat(length), b"A", Vec::new());

            assert_eq!(
                stratum_of(&locus),
                LocusStratum::Stratified(stratum(1, u32::try_from(length).expect("small"))),
                "a run of {length} A's"
            );
        }
    }

    /// A locus of any other kind is passed over in silence rather than counted: the SNP/indel
    /// path has its own accumulator, and a repeat *cluster* is not one delimited tract, so
    /// neither is a fault to report.
    #[test]
    fn a_locus_that_is_not_one_repeat_tract_has_no_stratum_and_is_not_a_fault() {
        let mut generic = ssr_locus(b"ATATATATAT", b"AT", Vec::new());
        generic.kind = LocusKind::Generic;
        assert_eq!(stratum_of(&generic), LocusStratum::NotOneRepeatTract);

        let mut bundle = ssr_locus(b"ATATATATAT", b"AT", Vec::new());
        bundle.kind = LocusKind::SsrBundle;
        assert_eq!(stratum_of(&bundle), LocusStratum::NotOneRepeatTract);
    }

    /// A tract of no bases at all is not something the catalog emits — its floors admit nothing
    /// below three copies — but the arithmetic has to answer for it, and zero copies is what
    /// zero bases is. Pinned so that the answer is a decision rather than an accident.
    ///
    /// It also pins **where the length comes from**: zero is the only length at which this
    /// fixture's region and its reference bases disagree, the helper giving every locus a region
    /// at least one base long, so a rule reading `region.len()` instead of the tract answers one
    /// copy here and agrees everywhere else.
    #[test]
    fn a_reference_tract_of_no_bases_is_a_stratum_of_no_copies() {
        let locus = ssr_locus(b"", b"AT", Vec::new());

        assert_eq!(stratum_of(&locus), LocusStratum::Stratified(stratum(2, 0)));
    }

    proptest::proptest! {
        /// **A tract is stratified exactly when its period divides its length, and its count is
        /// then the quotient** — stated as a law over every period the scope admits and every
        /// remainder each of them can leave, so no period and no remainder goes untried. The
        /// tables above are the readable pins; this is the one that cannot be satisfied by a
        /// rule that happens to fit them.
        #[test]
        fn a_tract_is_stratified_exactly_when_its_period_divides_its_length(
            period in 1usize..=MAX_MOTIF_LEN,
            tract_bases in 0usize..200,
        ) {
            let motif = &b"ACGTAC"[..period];
            let bases: Vec<u8> = motif.iter().copied().cycle().take(tract_bases).collect();

            let expected = if tract_bases.is_multiple_of(period) {
                LocusStratum::Stratified(stratum(
                    u8::try_from(period).expect("a period inside the STR scope"),
                    u32::try_from(tract_bases / period).expect("under 200 bases"),
                ))
            } else {
                LocusStratum::WithoutWholeRepeatReference {
                    tract_len: Bp(tract_bases as u64),
                    period: SsrPeriod::try_new(period).expect("a period inside the STR scope"),
                }
            };

            proptest::prop_assert_eq!(stratum_of(&ssr_locus(&bases, motif, Vec::new())), expected);
        }
    }
}
