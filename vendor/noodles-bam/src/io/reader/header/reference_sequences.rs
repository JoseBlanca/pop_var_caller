mod reference_sequence;

use std::io::{self, Read};

use noodles_sam::header::ReferenceSequences;

use self::reference_sequence::read_reference_sequence;
use crate::io::reader::num::read_u32_le;

pub(super) fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
where
    R: Read,
{
    let n_ref = read_u32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    // pvc fork (FORK.md, change 2): the declared count no longer sizes the map — 2^31 of them
    // asked for 38 GB before the first entry was read. A real header's few thousand are reserved
    // as they come.
    let mut reference_sequences = ReferenceSequences::with_capacity(n_ref.min(1 << 12));

    for _ in 0..n_ref {
        let (name, reference_sequence) = read_reference_sequence(reader)?;
        reference_sequences.insert(name, reference_sequence);
    }

    Ok(reference_sequences)
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use bstr::BString;
    use noodles_sam::header::record::value::{Map, map::ReferenceSequence};

    use super::*;

    /// pvc fork, change 2: a reference-sequence count the stream cannot fill is an end-of-file
    /// error, reached without reserving room for every declared entry.
    #[test]
    fn test_read_reference_sequences_with_a_count_past_the_end() {
        let data = [
            0xff, 0xff, 0xff, 0x7f, // n_ref = 2^31 - 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];
        let mut reader = &data[..];
        let error = read_reference_sequences(&mut reader).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_read_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let data = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequences(&mut reader)?;

        let expected: ReferenceSequences = [(
            BString::from("sq0"),
            Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
        )]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);

        Ok(())
    }
}
