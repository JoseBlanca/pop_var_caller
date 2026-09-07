use std::io;

use super::{read_states, state_renormalize, table};
use crate::{
    codecs::rans_4x8::ALPHABET_SIZE,
    io::reader::num::{read_itf8_as, read_u8},
};

type Frequencies = [u16; ALPHABET_SIZE]; // F

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    #[cfg(feature = "perf-counters")]
    let table_started = std::time::Instant::now();

    let frequencies = read_frequencies(src)?;

    // A fixed-size table, so the inner loop's index — a state masked to twelve bits — needs no
    // bounds check.
    let mut slots = Box::new([0; table::TOTAL_FREQUENCY]);
    table::fill_context(&frequencies, slots.as_mut_slice());

    #[cfg(feature = "perf-counters")]
    let table_nanos = table_started.elapsed().as_nanos() as u64;
    #[cfg(feature = "perf-counters")]
    let decode_started = std::time::Instant::now();

    let mut states = read_states(src)?;

    for chunk in dst.chunks_mut(states.len()) {
        for (d, state) in chunk.iter_mut().zip(states.iter_mut()) {
            let slot = slots[(*state & 0x0fff) as usize];

            *d = table::symbol(slot);

            *state = table::frequency(slot) * (*state >> 12) + (*state & 0x0fff)
                - table::range_start(slot);
            *state = state_renormalize(*state, src)?;
        }
    }

    #[cfg(feature = "perf-counters")]
    crate::perf::record_rans(
        false,
        dst.len(),
        table_nanos,
        decode_started.elapsed().as_nanos() as u64,
    );

    Ok(())
}

pub(super) fn read_frequencies(src: &mut &[u8]) -> io::Result<Frequencies> {
    const NUL: u8 = 0x00;

    let mut frequencies = [0; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        let f = read_itf8_as(src)?;
        frequencies[usize::from(sym)] = f;

        sym = read_u8(src)?;

        if sym == NUL {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                let f = read_itf8_as(src)?;
                frequencies[usize::from(sym)] = f;
                sym += 1;
            }
        }

        prev_sym = sym;
    }

    Ok(frequencies)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_frequencies() -> io::Result<()> {
        let src = [
            b'a', // symbol = 'a'
            0x05, // frequencies['a'] = 5
            b'b', // symbol = 'b'
            0x02, // run length = 2
            0x02, // frequencies['b'] = 2
            0x01, // frequencies['c'] = 1
            0x01, // frequencies['d'] = 1
            b'r', // symbol = 'r'
            0x02, // frequencies['r'] = 2
            0x00, // EOF
        ];

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(read_frequencies(&mut &src[..])?, expected);

        Ok(())
    }

    /// The same frequencies the test above reads, laid out as slots: every slot of a
    /// symbol's range answers with that symbol, its frequency, and where its range starts.
    #[test]
    fn a_slot_table_hands_each_slot_to_the_symbol_whose_range_covers_it() {
        let mut frequencies = [0; ALPHABET_SIZE];
        frequencies[usize::from(b'a')] = 5;
        frequencies[usize::from(b'b')] = 2;
        frequencies[usize::from(b'c')] = 1;
        frequencies[usize::from(b'd')] = 1;
        frequencies[usize::from(b'r')] = 2;

        let mut slots = vec![0; table::TOTAL_FREQUENCY];
        table::fill_context(&frequencies, &mut slots);

        // 'a' owns 0..5, 'b' 5..7, 'c' 7..8, 'd' 8..9, 'r' 9..11.
        for (slot, (symbol, frequency, range_start)) in [
            (0, (b'a', 5, 0)),
            (4, (b'a', 5, 0)),
            (5, (b'b', 2, 5)),
            (6, (b'b', 2, 5)),
            (7, (b'c', 1, 7)),
            (8, (b'd', 1, 8)),
            (9, (b'r', 2, 9)),
            (10, (b'r', 2, 9)),
        ] {
            assert_eq!(table::symbol(slots[slot]), symbol, "slot {slot}");
            assert_eq!(table::frequency(slots[slot]), frequency, "slot {slot}");
            assert_eq!(table::range_start(slots[slot]), range_start, "slot {slot}");
        }
    }

    /// **A context with one symbol gives it every slot**, so its frequency is 4,096 — one more
    /// than twelve bits hold, which is why the packing stores the frequency less one.
    #[test]
    fn a_single_symbol_context_packs_its_full_frequency() {
        let mut frequencies = [0; ALPHABET_SIZE];
        frequencies[usize::from(b'N')] = table::TOTAL_FREQUENCY as u16;

        let mut slots = vec![0; table::TOTAL_FREQUENCY];
        table::fill_context(&frequencies, &mut slots);

        for slot in [0, 1, table::TOTAL_FREQUENCY - 1] {
            assert_eq!(table::symbol(slots[slot]), b'N');
            assert_eq!(table::frequency(slots[slot]), table::TOTAL_FREQUENCY as u32);
            assert_eq!(table::range_start(slots[slot]), 0);
        }
    }
}
