//! One record's head: the fixed fields at the front of every record that let a reader
//! decide whether it wants the record without building it.
//!
//! ```text
//! record = position_offset | reference_span | non_reference_reads | body_bytes | body
//!          └─────────────────────── the head ──────────────────────────┘   └─ skip ─┘
//! ```
//!
//! A reader takes the head, decides, and either builds the body or advances `body_bytes`
//! past it; nothing else in the block has to be touched to make that decision. **Measured
//! on a tomato accession at three reads a position, 7.69 M records: a walk keeping one
//! record in a hundred takes 0.141 s against 0.29 s for one that builds every record —
//! 2.06× faster** (spec §4.3).
//!
//! **The head is not free, and most of its price is not the length field.** It costs 9.2 %
//! of the file at three reads a position and 5.8 % at 279, of which the length field alone
//! is 1.4 % and 3.3 %. The rest is what skippability forces on the body: a record's
//! coverage and its chain ids would otherwise be coded as differences from the previous
//! record, and a reader that skips a body never sees those differences — so both restart at
//! every record instead (spec §4.3).

use crate::ng::types::GenomeRegion;

/// What a reader learns about a record before deciding to build it.
///
/// Every field is read from the record's head; none requires touching the body.
/// **`body_bytes` is what makes skipping possible at all** — the encoded bytes carry no
/// separators, so without it a reader that wants to reach the next record must decode
/// every variable-length integer in this one to find where it ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordHead {
    /// Where the record sits. **Absolute, and rebuilt by the reader** from the block's
    /// first position and the difference the head encodes — the difference restarts at
    /// every block boundary (spec §3.2).
    ///
    /// It is a region and not a position because a record widened by a deletion covers
    /// more than one base, so a reader indexed by position cannot work out what a record
    /// reaches from its start alone. The cohort merge names this as one of two things it
    /// asks of the format.
    pub region: GenomeRegion,
    /// Reads at this locus that supported something other than the reference, summed over
    /// the observations. **Zero and "nothing varies here" are the same condition**, which
    /// is what the cohort's first pass filters on.
    ///
    /// A count of reads rather than of alternative alleles: it answers *does anything vary
    /// here* just as well, and it also lets a reader apply a threshold.
    pub non_reference_reads: u32,
    /// Length of the body that follows, in bytes.
    pub body_bytes: u32,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Position};

    /// A head is read once per record on a path that runs at about twenty million records
    /// a second, and a reader that had to clone one — or that paid a pointer chase to
    /// reach it — would be paying that on every record it skips. `Copy` is the property
    /// that says it does not.
    #[test]
    fn a_head_is_copied_not_cloned_and_stays_small() {
        let head = RecordHead {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1_000),
                end: Position(1_002),
            },
            non_reference_reads: 3,
            body_bytes: 47,
        };
        let copied = head;
        assert_eq!(copied, head, "a head is Copy, so this is not a move");
        assert!(
            std::mem::size_of::<RecordHead>() <= 32,
            "the head is {} bytes; it is read once per record and holds no allocation",
            std::mem::size_of::<RecordHead>()
        );
    }

    /// The span is a field of the region rather than something a reader derives from the
    /// next record's start: a record widened by a deletion reaches past the record that
    /// follows it, so the distance between two starts is not a span.
    #[test]
    fn a_head_carries_a_span_wider_than_one_base() {
        let deletion = RecordHead {
            region: GenomeRegion {
                contig: ContigId(1),
                start: Position(90_667_287),
                end: Position(90_667_293),
            },
            non_reference_reads: 162,
            body_bytes: 231,
        };
        assert_eq!(deletion.region.len(), 7);
    }
}
