//! **Reading a decoded record's fields without a trait object per field.**
//!
//! `Record` reaches its caller through `sam::alignment::Record`, whose accessors each return a
//! `Box<dyn …>`: eight heap allocations for every record, before a byte is copied. A caller
//! that already knows it is holding a CRAM record does not need the dynamic dispatch, and a
//! caller filling its own flat buffers does not need the intermediate `RecordBuf` either.
//!
//! Every method here writes into a buffer the caller owns and clears it first, so a caller can
//! keep one set of buffers for a whole container.

use noodles_sam::alignment::record::cigar::Op;

use super::{Record, Sequence};
use crate::io::reader::container::slice::ReferenceSequence;

impl Record<'_> {
    /// Which `@RG` in the header this record declares, as the index the CRAM stores.
    pub fn read_group_index(&self) -> Option<usize> {
        self.read_group_id
    }

    /// The read's name, as stored — `None` where the file did not preserve names.
    pub fn name_bytes(&self) -> Option<&[u8]> {
        self.name.as_deref()
    }

    /// The read's bases, reconstructed against the reference where the CRAM stores only the
    /// differences.
    pub fn write_bases_into(&self, dst: &mut Vec<u8>) {
        dst.clear();

        if self.bam_flags.is_unmapped() || self.cram_flags.sequence_is_missing() {
            dst.extend_from_slice(&self.sequence);
            return;
        }

        let (reference_sequence, alignment_start) = match self.reference_sequence.as_ref() {
            // The caller's window and the file's embedded reference are indexed the same way:
            // both are bases starting at a position that is not 1.
            Some(
                ReferenceSequence::Embedded {
                    reference_start,
                    sequence,
                }
                | ReferenceSequence::Window {
                    reference_start,
                    sequence,
                },
            ) => {
                let alignment_start = usize::from(self.alignment_start.unwrap());
                let offset = usize::from(*reference_start);
                let offset_alignment_start =
                    noodles_core::Position::new(alignment_start - offset + 1).unwrap();
                (Some(*sequence), offset_alignment_start)
            }
            Some(ReferenceSequence::External { sequence, .. }) => {
                (Some((**sequence).as_ref()), self.alignment_start.unwrap())
            }
            None => (None, noodles_core::Position::MIN),
        };

        dst.reserve(self.read_length);
        let sequence = Sequence::new(
            reference_sequence,
            self.substitution_matrix.clone(),
            &self.features,
            alignment_start,
            self.read_length,
        );
        dst.extend(sequence.iter_concrete());
    }

    /// The read's quality scores, as the raw 0–93 values — not offset by 33.
    pub fn write_quality_scores_into(&self, dst: &mut Vec<u8>) {
        use noodles_sam::alignment::record::QualityScores as _;

        dst.clear();
        dst.reserve(self.read_length);

        if self.bam_flags.is_unmapped() || self.cram_flags.quality_scores_are_stored_as_array() {
            dst.extend_from_slice(&self.quality_scores);
        } else {
            let scores = super::QualityScores::new(&self.features, self.read_length);
            for score in scores.iter() {
                dst.push(score.unwrap_or(0xff));
            }
        }
    }

    /// The read's alignment, as CIGAR operations.
    pub fn write_cigar_into(&self, dst: &mut Vec<Op>) {
        use noodles_sam::alignment::record::Cigar as _;

        dst.clear();

        let cigar = super::Cigar::new(
            &self.features,
            self.bam_flags.is_unmapped(),
            self.read_length,
        );
        for op in cigar.iter() {
            match op {
                Ok(op) => dst.push(op),
                Err(_) => return,
            }
        }
    }
}
