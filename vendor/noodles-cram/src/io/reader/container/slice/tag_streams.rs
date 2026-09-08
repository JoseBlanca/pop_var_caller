//! **Can a slice's auxiliary tags be left unread altogether?**
//!
//! A caller that reads no tags would like not to decode them. Whether it may is a property of
//! the file, not a choice: CRAM's data series and tags are decoded out of the same set of
//! streams, and a stream is a cursor. Skipping a read that some other read depends on having
//! happened leaves every value after it wrong, silently.
//!
//! So this decides the question from the compression header, before any record is read, and
//! answers *no* whenever it cannot see that skipping is safe. Two conditions, and both must
//! hold for **every** tag the container declares:
//!
//! 1. **The tag's encoding must read nothing from the core data stream.** That stream is one
//!    bit reader shared by every data series, so a tag that consumes bits from it is part of
//!    the same sequence as the fields around it and cannot be stepped over.
//! 2. **The external blocks the tag reads must be read by no data series.** An external block
//!    is its own cursor, so leaving one untouched costs nothing — unless something else reads
//!    the same block, in which case not advancing it desynchronises that reader instead.
//!
//! In practice both hold for CRAMs written by samtools, which gives every tag its own external
//! block. They do not hold for a writer that Huffman-codes a tag into the core stream, and for
//! such a file the answer here is *no* and the tags are read as they always were.

use std::collections::HashSet;

use crate::container::{
    block,
    compression_header::{
        CompressionHeader, Encoding,
        encoding::codec::{Byte, ByteArray, Integer},
    },
};

/// Whether every tag in this container can be left undecoded.
///
/// `false` is always a safe answer and is the one returned whenever anything is unrecognised.
pub(super) fn tags_can_be_left_unread(compression_header: &CompressionHeader) -> bool {
    tag_only_blocks(compression_header).is_some()
}

/// The external blocks **only the tags read**, or `None` when the tags may not be left unread.
///
/// `Some(set)` is the same answer [`tags_can_be_left_unread`] gives as `true`, carrying the
/// blocks that answer makes dead: a caller reading no tags will never look inside any of them,
/// so their bytes need not be decompressed at all. The set is disjoint from every external
/// block a data series reads — that disjointness *is* the condition — so dropping them cannot
/// starve a reader that is still consulted.
pub(super) fn tag_only_blocks(
    compression_header: &CompressionHeader,
) -> Option<HashSet<block::ContentId>> {
    let mut data_series_blocks = HashSet::new();
    if !collect_data_series_blocks(compression_header, &mut data_series_blocks) {
        // A data series reads the core stream, which is expected and fine; what is collected
        // here is only the set of external blocks. `collect_data_series_blocks` returns false
        // when it meets an encoding it does not recognise, and then nothing is skipped.
        return None;
    }

    let mut tag_blocks = HashSet::new();
    for encoding in compression_header.tag_encodings().values() {
        if !collect_byte_array_blocks(encoding, &mut tag_blocks) {
            return None;
        }
    }

    tag_blocks
        .is_disjoint(&data_series_blocks)
        .then_some(tag_blocks)
}

/// The external blocks every data series reads. `false` if an encoding is not one this
/// function knows how to walk, in which case its blocks are unknown and nothing may be skipped.
fn collect_data_series_blocks(
    compression_header: &CompressionHeader,
    blocks: &mut HashSet<block::ContentId>,
) -> bool {
    let encodings = compression_header.data_series_encodings();

    let integers = [
        encodings.bam_flags(),
        encodings.cram_flags(),
        encodings.reference_sequence_ids(),
        encodings.read_lengths(),
        encodings.alignment_starts(),
        encodings.read_group_ids(),
        encodings.mate_flags(),
        encodings.mate_reference_sequence_ids(),
        encodings.mate_alignment_starts(),
        encodings.template_lengths(),
        encodings.mate_distances(),
        encodings.tag_set_ids(),
        encodings.feature_counts(),
        encodings.feature_position_deltas(),
        encodings.deletion_lengths(),
        encodings.reference_skip_lengths(),
        encodings.padding_lengths(),
        encodings.hard_clip_lengths(),
        encodings.mapping_qualities(),
    ];
    for encoding in integers.into_iter().flatten() {
        collect_integer_blocks(encoding, blocks);
    }

    let bytes = [
        encodings.feature_codes(),
        encodings.base_substitution_codes(),
        encodings.bases(),
        encodings.quality_scores(),
    ];
    for encoding in bytes.into_iter().flatten() {
        collect_byte_blocks(encoding, blocks);
    }

    let byte_arrays = [
        encodings.names(),
        encodings.stretches_of_bases(),
        encodings.stretches_of_quality_scores(),
        encodings.insertion_bases(),
        encodings.soft_clip_bases(),
    ];
    for encoding in byte_arrays.into_iter().flatten() {
        // A data series may read the core stream; that is not a reason to refuse, because a
        // data series is read either way. Only its external blocks are wanted here, so the
        // return value is deliberately ignored.
        collect_byte_array_blocks(encoding, blocks);
    }

    true
}

/// `false` when the encoding reads the core data stream, which makes it unskippable.
fn collect_byte_array_blocks(
    encoding: &Encoding<ByteArray>,
    blocks: &mut HashSet<block::ContentId>,
) -> bool {
    match encoding.get() {
        ByteArray::ByteArrayStop {
            block_content_id, ..
        } => {
            blocks.insert(*block_content_id);
            true
        }
        ByteArray::ByteArrayLength {
            len_encoding,
            value_encoding,
        } => {
            // **Both sides are collected before they are combined, and `&&` would not do
            // that.** It short-circuits, so a length encoding that reads the core stream —
            // which returns `false` and is perfectly ordinary — would leave the *value*
            // encoding's external block uninserted. On the tag side that is harmless: the
            // `false` propagates and nothing is skipped. On the data-series side it is not,
            // because `collect_data_series_blocks` deliberately ignores this return value for
            // its `byte_arrays` encodings — so such a series' external block would be missing
            // from the set a tag block is tested for disjointness against, and a tag sharing
            // it would look safe to skip. Skipping it desynchronises the reader that is still
            // consulted, which is wrong values rather than an error. It does not fire on the
            // files here (samtools writes read names as `ByteArrayStop`, the other arm).
            let length_recognised = collect_integer_blocks(len_encoding, blocks);
            let value_recognised = collect_byte_blocks(value_encoding, blocks);
            length_recognised && value_recognised
        }
    }
}

/// `false` when the encoding reads the core data stream.
fn collect_integer_blocks(
    encoding: &Encoding<Integer>,
    blocks: &mut HashSet<block::ContentId>,
) -> bool {
    match encoding.get() {
        Integer::External { block_content_id } => {
            blocks.insert(*block_content_id);
            true
        }
        // **A one-symbol Huffman alphabet reads no bits**, and it is not a curiosity: it is
        // how a writer stores a value that is the same in every record, which is what the
        // length of a fixed-width tag is. Refusing it would refuse most real files.
        Integer::Huffman { alphabet, .. } if alphabet.len() == 1 => true,
        // Golomb, Huffman over more than one symbol, Beta, Subexp, GolombRice and Gamma all
        // read bits from the core stream. Listed rather than matched with a wildcard so that a
        // codec added later is a compile error here instead of a silent "skippable".
        Integer::Golomb { .. }
        | Integer::Huffman { .. }
        | Integer::Beta { .. }
        | Integer::Subexp { .. }
        | Integer::GolombRice { .. }
        | Integer::Gamma { .. } => false,
    }
}

/// `false` when the encoding reads the core data stream.
fn collect_byte_blocks(encoding: &Encoding<Byte>, blocks: &mut HashSet<block::ContentId>) -> bool {
    match encoding.get() {
        Byte::External { block_content_id } => {
            blocks.insert(*block_content_id);
            true
        }
        // As above: one symbol is a constant, and reads nothing.
        Byte::Huffman { alphabet, .. } if alphabet.len() == 1 => true,
        Byte::Huffman { .. } => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::container::compression_header::{
        DataSeriesEncodings, PreservationMap, TagEncodings, data_series_encodings::DataSeries,
    };

    /// A container with the standard data-series encodings and one tag.
    fn header_with(
        tag_block: block::ContentId,
        tag_encoding: Encoding<ByteArray>,
    ) -> CompressionHeader {
        let mut tag_encodings = TagEncodings::default();
        tag_encodings.insert(tag_block, tag_encoding);

        CompressionHeader::new(
            PreservationMap::default(),
            DataSeriesEncodings::init(),
            tag_encodings,
        )
    }

    /// The content id a real tag gets: `NM` of type `c`, packed as the specification says.
    const NM_C: block::ContentId = ((b'N' as i32) << 16) | ((b'M' as i32) << 8) | (b'c' as i32);

    #[test]
    fn a_tag_in_its_own_external_block_can_be_left_unread() {
        let header = header_with(
            NM_C,
            Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0,
                block_content_id: NM_C,
            }),
        );
        assert!(tags_can_be_left_unread(&header));
    }

    #[test]
    fn a_tag_that_reads_the_core_stream_may_not_be_left_unread() {
        // Two symbols, so the Huffman decoder consumes bits from the shared core reader.
        let header = header_with(
            NM_C,
            Encoding::new(ByteArray::ByteArrayLength {
                len_encoding: Encoding::new(Integer::Huffman {
                    alphabet: vec![1, 2],
                    bit_lens: vec![1, 1],
                }),
                value_encoding: Encoding::new(Byte::External {
                    block_content_id: NM_C,
                }),
            }),
        );
        assert!(!tags_can_be_left_unread(&header));
    }

    #[test]
    fn a_one_symbol_huffman_length_reads_nothing_and_is_allowed() {
        let header = header_with(
            NM_C,
            Encoding::new(ByteArray::ByteArrayLength {
                len_encoding: Encoding::new(Integer::Huffman {
                    alphabet: vec![4],
                    bit_lens: vec![0],
                }),
                value_encoding: Encoding::new(Byte::External {
                    block_content_id: NM_C,
                }),
            }),
        );
        assert!(tags_can_be_left_unread(&header));
    }

    #[test]
    fn a_tag_sharing_a_block_with_a_data_series_may_not_be_left_unread() {
        let header = header_with(
            NM_C,
            Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0,
                block_content_id: block::ContentId::from(DataSeries::BamFlags),
            }),
        );
        assert!(!tags_can_be_left_unread(&header));
    }
}
