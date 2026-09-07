use std::io;

use super::{order_0, read_states, state_renormalize, table};
use crate::{
    codecs::rans_4x8::{ALPHABET_SIZE, STATE_COUNT},
    io::reader::num::read_u8,
};

type Frequencies = [[u16; ALPHABET_SIZE]; ALPHABET_SIZE]; // F

/// The slot tables of the contexts this block actually uses, back to back.
///
/// The order-1 model keys on the previous symbol, so it declares 256 contexts. A block of DNA
/// bases uses about five of them and a block of quality scores a few dozen; the rest have an
/// all-zero frequency row that no state can ever reach, because the decoder indexes context
/// *i* only after emitting symbol *i*. Holding only the used ones is why a block allocates
/// tens of kilobytes here rather than the four megabytes 256 packed contexts would take.
struct DecodeTables {
    /// One fixed-size table per used context. Fixed-size rather than a flat `Vec` and an
    /// offset **so that the inner loop's slot index needs no bounds check**: the index is a
    /// state masked to twelve bits, which the compiler can see is below the row's length.
    rows: Vec<[table::Slot; table::TOTAL_FREQUENCY]>,
    /// Which row each context uses. An unused context is never indexed, and points at row zero
    /// so that a corrupt stream reads a wrong byte rather than panicking in the inner loop.
    row_of_context: [u16; ALPHABET_SIZE],
}

impl DecodeTables {
    fn build(frequencies: &Frequencies) -> Self {
        let mut rows = Vec::new();
        let mut row_of_context = [0u16; ALPHABET_SIZE];

        for (context, frequencies) in frequencies.iter().enumerate() {
            if frequencies.iter().all(|&frequency| frequency == 0) {
                continue;
            }

            row_of_context[context] = rows.len() as u16;
            rows.push([0; table::TOTAL_FREQUENCY]);
            let row = rows.last_mut().expect("a row was just pushed");
            table::fill_context(frequencies, row);
        }

        if rows.is_empty() {
            rows.push([0; table::TOTAL_FREQUENCY]);
        }

        Self {
            rows,
            row_of_context,
        }
    }
}

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    #[cfg(feature = "perf-counters")]
    let table_started = std::time::Instant::now();

    let frequencies = read_frequencies(src)?;
    let tables = DecodeTables::build(&frequencies);

    #[cfg(feature = "perf-counters")]
    let table_nanos = table_started.elapsed().as_nanos() as u64;
    #[cfg(feature = "perf-counters")]
    let decode_started = std::time::Instant::now();

    let mut states = read_states(src)?;
    let mut prev_syms = [0; STATE_COUNT];

    let [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4] = split_chunks(dst);
    let chunks = chunk_0.iter_mut().zip(chunk_1).zip(chunk_2).zip(chunk_3);

    for (((d0, d1), d2), d3) in chunks {
        let dsts = [d0, d1, d2, d3];

        for (state, (prev_sym, d)) in states.iter_mut().zip(prev_syms.iter_mut().zip(dsts)) {
            let row = &tables.rows[usize::from(tables.row_of_context[usize::from(*prev_sym)])];
            let slot = row[(*state & 0x0fff) as usize];
            let sym = table::symbol(slot);

            *d = sym;

            *state = table::frequency(slot) * (*state >> 12) + (*state & 0x0fff)
                - table::range_start(slot);
            *state = state_renormalize(*state, src)?;

            *prev_sym = sym;
        }
    }

    let mut state = states[3];
    let mut prev_sym = prev_syms[3];

    for d in chunk_4 {
        let row = &tables.rows[usize::from(tables.row_of_context[usize::from(prev_sym)])];
        let slot = row[(state & 0x0fff) as usize];
        let sym = table::symbol(slot);

        *d = sym;

        state = table::frequency(slot) * (state >> 12) + (state & 0x0fff)
            - table::range_start(slot);
        state = state_renormalize(state, src)?;

        prev_sym = sym;
    }

    #[cfg(feature = "perf-counters")]
    crate::perf::record_rans(
        true,
        dst.len(),
        table_nanos,
        decode_started.elapsed().as_nanos() as u64,
    );

    Ok(())
}

fn read_frequencies(src: &mut &[u8]) -> io::Result<Frequencies> {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        let f = order_0::read_frequencies(src)?;
        frequencies[usize::from(sym)] = f;

        sym = read_u8(src)?;

        if sym == 0 {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                let f = order_0::read_frequencies(src)?;
                frequencies[usize::from(sym)] = f;
                sym += 1;
            }
        }

        prev_sym = sym;
    }

    Ok(frequencies)
}

fn split_chunks(dst: &mut [u8]) -> [&mut [u8]; 5] {
    let chunk_size = dst.len() / STATE_COUNT;

    let (left_chunk, right_chunk) = dst.split_at_mut(2 * chunk_size);
    let (chunk_0, chunk_1) = left_chunk.split_at_mut(chunk_size);
    let (chunk_2, chunk_3_4) = right_chunk.split_at_mut(chunk_size);
    let (chunk_3, chunk_4) = chunk_3_4.split_at_mut(chunk_size);

    [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_frequencies() -> io::Result<()> {
        const NUL: u8 = 0x00;

        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07):
        // "abracadabraabracadabraabracadabraabracadabrad".
        let src = [
            0x00, // symbols[0] = '\0' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['\0']['a'] = 4
            0x00, // }
            b'a', // symbols[0] = 'a' {
            b'a', //   symbols[1] = 'a'
            0x03, //     frequencies['a']['a'] = 3
            b'b', //   symbols[1] = 'b'
            0x02, //     run length = 2
            0x08, //     frequencies['a']['b'] = 8
            0x04, //     frequencies['a']['c'] = 4
            0x05, //     frequencies['a']['d'] = 5
            0x00, // }
            b'b', // symbols[0] = 'b' {
            0x02, //   run length = 2
            b'r', //   symbols[1] = 'r'
            0x08, //     frequencies['b']['r'] = 8
            0x00, // }
            //    // symbols[0] = 'c' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['c']['a'] = 4
            0x00, // }
            //    // symbols[0] = 'd' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['d']['a'] = 4
            0x00, // }
            b'r', // symbols[0] = 'r' {
            b'a', //   symbols[1] = 'a'
            0x08, //     frequencies['r']['a'] = 8
            0x00, // }
            0x00, // EOF
        ];

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 4;
        expected[usize::from(b'a')][usize::from(b'a')] = 3;
        expected[usize::from(b'a')][usize::from(b'b')] = 8;
        expected[usize::from(b'a')][usize::from(b'c')] = 4;
        expected[usize::from(b'a')][usize::from(b'd')] = 5;
        expected[usize::from(b'b')][usize::from(b'r')] = 8;
        expected[usize::from(b'c')][usize::from(b'a')] = 4;
        expected[usize::from(b'd')][usize::from(b'a')] = 4;
        expected[usize::from(b'r')][usize::from(b'a')] = 8;

        assert_eq!(read_frequencies(&mut &src[..])?, expected);

        Ok(())
    }
}
