//! What a segmentation derived from the catalog contained, and what admission turned down
//! getting there.
//!
//! **This is the catalog's answer to a scan's own running tally**
//! ([`crate::ng::region_typing::TypedRegionCounts`]): a consumer that stops scanning the
//! reference and reads the file instead must not lose its report. Every counter there has a
//! counterpart here and holds the same number, with **one deliberate omission**.
//!
//! ## The one question the file cannot answer
//!
//! A scan counts repeats it turns down for **touching a contig's very first or last
//! base** — a tract with no sequence on one side has nothing for a read to anchor against,
//! so it is not an STR (`RejectionReason::FlankClamped`).
//!
//! **The catalog cannot count those, and there is no field here for them.** The file only
//! holds a repeat with at least 15 bases of sequence beside it on each side (spec §4.1), so
//! a repeat sitting on base 1 was never written down and no reader can see it. A counter
//! here would therefore report `0` for every genome — zero by construction, not zero in the
//! genome — and nothing downstream would know the difference. A field that cannot exist
//! cannot be misread; a field holding `0` will be.
//!
//! Two consequences a reader of these numbers has to know, and neither is large:
//!
//! - a repeat 1 to 14 bases from a contig's end is a **locus** to a scan and is absent
//!   from the file altogether, so [`ssr_loci`](CatalogRegionCounts::ssr_loci) and
//!   [`generic`](CatalogRegionCounts::generic) can differ there too. Measured over whole
//!   genomes: 7 such tracts in tomato SL4.00 and 857 in GRCh38;
//! - everything else agrees exactly, because both sides run the same pre-screen and the same
//!   admission gates — the file merely supplies the whole-motif cut, the motif and the purity
//!   instead of the bases.

use crate::ng::region_typing::segment_criteria::RejectionReason;

/// The running tally of a segmentation derived from the catalog — readable mid-read,
/// complete once the iterator is exhausted
/// ([`GenomeSegments::counts`](crate::ng::repeat_catalog::GenomeSegments::counts)).
///
/// Field for field a scan's [`TypedRegionCounts`](crate::ng::region_typing::TypedRegionCounts),
/// minus the contig-end rejection the module doc explains.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CatalogRegionCounts {
    /// Regions the reader asked for. One per contig when the whole reference was asked for,
    /// a zero-length contig contributing none — which is what
    /// [`GenomeRegions::whole_contigs`](crate::ng::region_typing::GenomeRegions::whole_contigs)
    /// produces.
    pub spans: u64,
    /// STR loci emitted.
    pub ssr_loci: u64,
    /// Clusters of repeats emitted, none of whose members has clean flanks.
    pub ssr_bundles: u64,
    /// Bases those clusters cover.
    pub ssr_bundle_bp: u64,
    /// Generic stretches emitted — the spaces between the repeat features.
    pub generic: u64,
    /// Tandem arrays too long to be microsatellites.
    pub satellites: u64,
    /// Bases those arrays cover.
    pub satellite_bp: u64,
    /// Repeat coverage that yielded **no locus**, in bases: every base of cleaned repeat
    /// coverage that did not come out as an STR locus — because it was bundled, capped as a
    /// satellite, or turned down by one of admission's gates.
    ///
    /// In bases rather than per repeat, for the reason a scan gives: admission trims every
    /// survivor, so a repeat that yields one locus and sheds 200 bases would contribute
    /// nothing to a per-repeat counter.
    ///
    /// **Mid-read this leads.** A contig's whole repeat coverage is charged when that contig
    /// is read, and the loci that cancel it are subtracted as they are handed out. It is
    /// exact once the iterator is exhausted, which is what this type promises.
    pub repeat_bp_with_no_locus: u64,
    /// [`Self::repeat_bp_with_no_locus`] broken out by the gate that turned the repeat down.
    ///
    /// **It does not partition the total**, for the same two reasons a scan's does not:
    /// overlapping rejected repeats are both charged, and bases with no locus for a reason
    /// that is not a rejection (bundled, capped as a satellite, out of the period range) are
    /// in the total and under no reason here. A diagnosis of admission's gates, not an
    /// account of the genome.
    pub rejected_by_reason: CatalogRejectionCounts,
}

impl CatalogRegionCounts {
    /// Charge one contig's repeat coverage and its rejections.
    pub(crate) fn add_contig(&mut self, tally: &ContigTally) {
        self.repeat_bp_with_no_locus += tally.repeat_bp;
        let r = &mut self.rejected_by_reason;
        r.copy_floor += tally.rejected.copy_floor;
        r.purity += tally.rejected.purity;
        r.compound += tally.rejected.compound;
        r.no_clean_trim += tally.rejected.no_clean_trim;
    }
}

/// Repeat bases admission turned down, by reason — **four of a scan's five**, and the
/// module doc says which is missing and why.
///
/// The bases charged are the repeat's **detected** length, before any trim, exactly as a
/// scan charges them.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CatalogRejectionCounts {
    /// Fewer whole motif copies in the cut tract than the reader's floor for that period.
    pub copy_floor: u64,
    /// Below the reader's purity floor once measured against a perfect motif tiling.
    pub purity: u64,
    /// The motif is itself a repeat (`ATAT` is `AT` twice), so the period is a lie.
    pub compound: u64,
    /// No whole-motif boundaries to cut back to — nothing in the tract tiles.
    pub no_clean_trim: u64,
}

impl CatalogRejectionCounts {
    /// Charge `bp` to `reason`.
    ///
    /// [`RejectionReason::FlankClamped`] is **not chargeable and cannot arrive.** The clamp
    /// fires only for a tract that abuts a contig's very first or last base, and the builder
    /// stores no tract with less than its flank floor beside it — which the command refuses
    /// to set below 1 for exactly this reason (`repeat-catalog`'s `NoFlankFloor`). So one
    /// base of flank is what the property needs; the default of 15 is far above it. A row
    /// that did clamp would be a builder bug rather than a policy rejection, which is why it
    /// trips a `debug_assert` instead of quietly landing in a counter that would then mean
    /// two different things.
    pub(crate) fn add(&mut self, reason: RejectionReason, bp: u64) {
        let slot = match reason {
            RejectionReason::CopyFloor => &mut self.copy_floor,
            RejectionReason::Purity => &mut self.purity,
            RejectionReason::Compound => &mut self.compound,
            RejectionReason::NoCleanTrim => &mut self.no_clean_trim,
            RejectionReason::FlankClamped => {
                debug_assert!(
                    false,
                    "a catalog row cannot clamp its flank — the file's 15 bp flank floor \
                     excludes every tract that could ({bp} bp charged)"
                );
                return;
            }
        };
        *slot += bp;
    }

    /// Bases turned down, over every reason this type holds.
    pub fn total(&self) -> u64 {
        self.copy_floor + self.purity + self.compound + self.no_clean_trim
    }
}

/// One contig's contribution to [`CatalogRegionCounts`] — the parts that cannot be read off
/// the regions themselves, because they are about repeats that never became one.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ContigTally {
    /// Cleaned repeat coverage on this contig, in bases. The whole contig's, not the
    /// requested stretches' — a region read still types the contig whole, exactly as the
    /// walk does, because a repeat outside the request still bundles with one inside it.
    pub repeat_bp: u64,
    /// What admission turned down on this contig.
    pub rejected: CatalogRejectionCounts,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_reason_lands_in_its_own_counter() {
        let mut counts = CatalogRejectionCounts::default();
        counts.add(RejectionReason::CopyFloor, 10);
        counts.add(RejectionReason::Purity, 20);
        counts.add(RejectionReason::Compound, 30);
        counts.add(RejectionReason::NoCleanTrim, 40);

        assert_eq!(
            counts,
            CatalogRejectionCounts {
                copy_floor: 10,
                purity: 20,
                compound: 30,
                no_clean_trim: 40,
            }
        );
        assert_eq!(counts.total(), 100);
    }

    /// A contig's coverage and its rejections both roll up; the contig-end reason has no
    /// slot to roll up into.
    #[test]
    fn a_contigs_tally_rolls_up_into_the_run_total() {
        let mut counts = CatalogRegionCounts::default();
        counts.add_contig(&ContigTally {
            repeat_bp: 500,
            rejected: CatalogRejectionCounts {
                copy_floor: 7,
                ..Default::default()
            },
        });
        counts.add_contig(&ContigTally {
            repeat_bp: 300,
            rejected: CatalogRejectionCounts {
                purity: 9,
                ..Default::default()
            },
        });

        assert_eq!(counts.repeat_bp_with_no_locus, 800);
        assert_eq!(counts.rejected_by_reason.copy_floor, 7);
        assert_eq!(counts.rejected_by_reason.purity, 9);
        assert_eq!(counts.rejected_by_reason.total(), 16);
    }
}
