//! Step 6 — **choosing the short list of sequences a locus is called over.**
//!
//! The cohort merge collects every sequence any sample's reads showed over a stretch of
//! genome and unifies them into one table; it narrows nothing
//! (`doc/devel/ng/spec/cohort_merge.md` §4.2). Everything downstream is defined over the
//! narrowed list — the read likelihood scores each observation against each candidate,
//! the genotype prior lays its mass over the genotypes those candidates make, and the
//! VCF's `ALT` column is what survived. This module does the narrowing
//! (`doc/devel/ng/spec/candidate_alleles.md`, and `doc/devel/ng/arch/candidate_alleles.md`
//! for the shapes).
//!
//! **A folder rather than a file, and no trait.** Two paths — the ordinary SNP/indel one
//! in `generic.rs` and the repeat tract's — take different evidence and return different
//! extras, and which runs is decided by the locus's kind rather than by a swappable
//! recipe. Two functions, not two impls of one seam. What both share lives here: the
//! config, the verdict, the leftover, the remapping and the ranking.
//!
//! **The rule, in one line.** An alternative survives if *some single sample's* reads
//! lent it at least `max(2 reads, 5 in 100 of that sample's reads at the locus)`, and a
//! locus is called over at most six alleles counting the reference. **No term of the bar
//! reads the cohort** — one sample reaching it admits the sequence for everyone —
//! because otherwise a sample's candidate list would depend on who else is in the run
//! (spec §3.2).

use crate::ng::run::cohort_merge::{MinAltObs, MinAltReadShare, MinAltReads};

/// **The support one sample must lend one sequence for it to be called over, and the cap
/// on how many sequences a locus is called over at all** — the two halves of the
/// narrowing, and the run's only knobs on it (spec §3, §4).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct CandidateSelectionConfig {
    /// `max(floor, ceil(share × that sample's reads at the locus))` non-reference reads,
    /// asked of **each sample separately against that sample's own reads**; one sample
    /// reaching it admits the sequence for the whole cohort (spec §3.2).
    ///
    /// **The merge's own type, reused rather than copied**, so the rule that decides
    /// whether a locus is built and the rule that decides which of its alleles are
    /// called over cannot drift apart — a sweep of one is a sweep of both. Only the
    /// number differs; see [`DEFAULT_ALLELE_SUPPORT`].
    pub support: MinAltReads,
    /// How many alleles a locus may be called over, **counting the reference**. Above it
    /// the list is cut to the best-ranked and the locus is still called; it is never
    /// refused (spec §4.1).
    ///
    /// **A bare `u16`, so 0 and 1 are representable and neither is a legal cap** — at
    /// either value the reference is the only survivor and every alternative becomes a
    /// truncation, which is refusal under another name and is what spec §4.1 rules out.
    /// Nothing validates it here, and that is an obligation this field hands to the fold:
    /// **`select_generic` asserts a cap of at least 2** when it lands (plan step C2). The
    /// alternative — a newtype refusing anything below 2, the shape
    /// [`MinAltObs`](crate::ng::run::cohort_merge::MinAltObs) and its neighbours already
    /// have — changes what arch §2.1 declares and is raised at Checkpoint A rather than
    /// taken here.
    pub max_candidate_alleles: u16,
}

impl CandidateSelectionConfig {
    /// [`DEFAULT_ALLELE_SUPPORT`] and [`DEFAULT_MAX_CANDIDATE_ALLELES`] — both soft, and
    /// both carrying their source in their own documentation.
    pub const DEFAULT: Self = Self {
        support: DEFAULT_ALLELE_SUPPORT,
        max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES,
    };
}

impl Default for CandidateSelectionConfig {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// **Two reads, or 5 in 100 of that sample's reads at the locus, whichever is more.**
///
/// **The floor is the merge's own** ([`MinAltObs::DEFAULT`], production's number) **and
/// should stay at 2**: measured against the GIAB trio's v4.2.1 truth set over 572 kb on
/// 2026-08-24, **at 30×** raising it from 2 to 3 loses five true alternative alleles,
/// where raising the share to 10 in 100 loses two for the same reduction in table size —
/// 1,539 alternatives kept against 1,601. The floor is the expensive knob (spec §3.3).
///
/// **The share is 5 in 100 where the merge's keep rule uses 2**, because the
/// allele-level question tolerates a stricter share at depth. On the same trio **at
/// 300×** it cuts the merge's 15,474 alternatives to 2,308, where a bar of 2 reads alone
/// keeps 10,793 — and it loses the same two true alleles the 2-in-100 share loses.
///
/// **It is inert below 41 compared reads a sample**, which is the arithmetic rather than
/// a measurement: `ceil(0.05 × 40) = 2` is the floor, and 41 is the first count at which
/// the share asks for more. So a tomato-depth run — about 11 compared reads a sample at a
/// locus — sees the identical rule it would have seen at 2 in 100. What was measured on
/// that panel is the neighbouring comparison: turning the share off entirely against the
/// merge's 2 in 100 moves 4 loci in 53,935 (spec §3.3).
///
/// **Soft.** Measured on one human trio over 572 kb (spec §11, Q3); what would move it is
/// the same scoring on a second high-depth cohort.
pub const DEFAULT_ALLELE_SUPPORT: MinAltReads = MinAltReads {
    floor: MinAltObs::DEFAULT,
    share: MinAltReadShare::new_const(0.05),
};

/// **Six alleles including the reference** — production's `DEFAULT_MAX_ALLELES_PER_RECORD`
/// (`src/var_calling/per_group_merger.rs`) and GATK's `--max-alternate-alleles` default,
/// inherited and declared inherited.
///
/// **Measured to be a safety valve rather than a working part at the cohort sizes we
/// have** (spec §4.2, 2026-08-24): it binds at 23 of 53,935 tomato loci — one in 2,300 —
/// and at none of the GIAB trio's 4,177 loci at 30× or 7,478 at 300×. What it guards
/// against grows with the cohort, which is why it is here: holding the tomato allele
/// table fixed and asking the bar of 1, 4, 16 and 63 samples gives 0, 0, 3 and 23 loci
/// above six alleles.
///
/// **Those counts were taken with the merge's 2-in-100 share, not the 5 in 100 this
/// module ships**, which spec §4.2 states in its own header and a reader of this constant
/// would otherwise not know. The direction is safe: at tomato depth the two shares are
/// provably the same rule (see [`DEFAULT_ALLELE_SUPPORT`]), and everywhere else 5 in 100
/// admits fewer alternatives, so the cap binds no more often than these numbers say.
///
/// **Soft, and never measured at its own value.** Whether it becomes load-bearing past a
/// few hundred samples is an extrapolation from that table, not a measurement
/// (spec §11, Q2).
pub const DEFAULT_MAX_CANDIDATE_ALLELES: u16 = 6;

#[cfg(test)]
mod tests {
    use super::*;

    /// The two numbers themselves, and the floor's **coupling** to the merge's constant
    /// rather than to the digit 2 — the doc comment says the floor *is* the merge's own,
    /// and only the second assertion holds that.
    #[test]
    fn the_default_bar_is_two_reads_or_five_in_a_hundred() {
        assert_eq!(DEFAULT_ALLELE_SUPPORT.floor.get(), 2);
        assert_eq!(DEFAULT_ALLELE_SUPPORT.floor, MinAltObs::DEFAULT);
        assert_eq!(DEFAULT_ALLELE_SUPPORT.share.get(), 0.05);
    }

    /// The two ends of the committed depth range, as spec §3 states them: at 3 compared
    /// reads the rule asks 2, and at 300 it asks 15.
    ///
    /// The third count is there because the first two cannot see the rounding: `0.05 ×
    /// 300` is exactly 15, so rounding the share *down* would answer 15 as well. At 301
    /// the share is 15.05, and up and down are 16 and 15.
    #[test]
    fn the_floor_decides_at_three_reads_and_the_share_at_three_hundred() {
        assert_eq!(DEFAULT_ALLELE_SUPPORT.required_of(3), 2);
        assert_eq!(DEFAULT_ALLELE_SUPPORT.required_of(300), 15);
        assert_eq!(DEFAULT_ALLELE_SUPPORT.required_of(301), 16);
    }

    /// The share is stricter than the merge's own at depth and **indistinguishable from
    /// it below 41 compared reads** — the claim [`DEFAULT_ALLELE_SUPPORT`]'s
    /// documentation makes, held against the merge's constant rather than against a
    /// number retyped here.
    ///
    /// **40 and 41 are the fixture, and stopping short of them is what makes this test
    /// vacuous.** With the equality arm ending at 20 compared reads, a share of 10 in 100
    /// also passes — `ceil(0.10 × 20) = 2` is still the floor — so the test would admit a
    /// bar twice as strict as the one shipped, which spec §3.3 measures losing two more
    /// true alleles at 300×. Carrying the pair to 40 and 41 pins the share to **more than
    /// 2/41 — about 0.0488 — and no more than 0.05**, which is the narrowest window the
    /// rule's own integer arithmetic can express.
    #[test]
    fn the_allele_share_binds_only_above_forty_compared_reads() {
        for compared_reads in [1_u32, 3, 11, 20, 40] {
            assert_eq!(
                DEFAULT_ALLELE_SUPPORT.required_of(compared_reads),
                MinAltReads::DEFAULT.required_of(compared_reads),
                "at {compared_reads} compared reads the floor decides for both rules"
            );
        }
        assert!(
            DEFAULT_ALLELE_SUPPORT.required_of(41) > MinAltReads::DEFAULT.required_of(41),
            "41 compared reads is where the allele rule first asks for more than the merge's"
        );
        assert!(
            DEFAULT_ALLELE_SUPPORT.required_of(300) > MinAltReads::DEFAULT.required_of(300),
            "at 300 compared reads the allele rule must ask for more than the merge's"
        );
    }

    #[test]
    fn the_cap_default_is_six_and_the_config_carries_it() {
        assert_eq!(DEFAULT_MAX_CANDIDATE_ALLELES, 6);
        assert_eq!(
            CandidateSelectionConfig::default().max_candidate_alleles,
            DEFAULT_MAX_CANDIDATE_ALLELES
        );
    }

    /// **The config a run gets by default is the two announced constants, and its support
    /// rule is the allele one rather than the merge's.**
    ///
    /// Both names are in scope in this file and both are a `MinAltReads`, so writing
    /// `MinAltReads::DEFAULT` here is a one-token slip that nothing else would catch: the
    /// type is right, the floor is right, and only the share moves — from 5 in 100 to 2 in
    /// 100. It is invisible at tomato depth, where the two are the same rule, and on the
    /// GIAB trio at 300× it is the difference between keeping 2,308 alternatives and
    /// keeping 5,596 (spec §3.3), each of which the genotype prior divides its
    /// concentration by. The last assertion is what makes the test about the *number* and
    /// not merely about the type.
    #[test]
    fn the_default_config_is_the_two_announced_constants_and_not_the_merges_rule() {
        let config = CandidateSelectionConfig::default();
        assert_eq!(config, CandidateSelectionConfig::DEFAULT);
        assert_eq!(config.support, DEFAULT_ALLELE_SUPPORT);
        assert_eq!(config.max_candidate_alleles, DEFAULT_MAX_CANDIDATE_ALLELES);
        assert_ne!(
            config.support,
            MinAltReads::DEFAULT,
            "the allele rule and the merge's keep rule share a type, not a share"
        );
    }
}
