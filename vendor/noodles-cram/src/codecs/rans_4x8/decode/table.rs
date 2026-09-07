//! **One lookup a byte instead of three.**
//!
//! Decoding a byte needs three numbers keyed on the same slot of the same context: which symbol
//! the slot belongs to, that symbol's frequency, and where its range starts. They were held in
//! three separate arrays, so the inner loop did three loads and two of them depended on the
//! result of the first. Here they are one 32-bit word, so it does one.
//!
//! ```text
//!   bits 31..20   where the symbol's range starts   (0..4095)
//!   bits 19..8    the symbol's frequency, less one  (1..4096 stored as 0..4095)
//!   bits  7..0    the symbol
//! ```
//!
//! **The frequency is stored less one because it does not fit otherwise.** Frequencies are
//! normalised to sum to 4,096, so a context with one symbol gives that symbol a frequency of
//! 4,096 — thirteen bits. Its range then starts at 0, and every other frequency is at most
//! 4,095, so subtracting one makes all three fields fit in exactly 32 bits with nothing spare.

use super::super::ALPHABET_SIZE;

/// How many slots a context's range is divided into — the sum every frequency table is
/// normalised to.
pub(super) const TOTAL_FREQUENCY: usize = 1 << 12;

/// One slot's answer: the symbol, its frequency less one, and where its range starts.
pub(super) type Slot = u32;

pub(super) fn pack(symbol: u8, frequency: u16, range_start: u16) -> Slot {
    (u32::from(range_start) << 20) | (u32::from(frequency - 1) << 8) | u32::from(symbol)
}

pub(super) fn symbol(slot: Slot) -> u8 {
    slot as u8
}

pub(super) fn frequency(slot: Slot) -> u32 {
    ((slot >> 8) & 0xfff) + 1
}

pub(super) fn range_start(slot: Slot) -> u32 {
    slot >> 20
}

/// Fill one context's 4,096 slots from its frequencies.
///
/// A symbol with a frequency of zero occupies no slot, which is what makes a slot's symbol
/// well defined: the slots are handed out in symbol order, `frequency` of them each.
pub(super) fn fill_context(frequencies: &[u16; ALPHABET_SIZE], slots: &mut [Slot]) {
    debug_assert_eq!(slots.len(), TOTAL_FREQUENCY);

    let mut range_start = 0u16;

    for (symbol, &frequency) in frequencies.iter().enumerate() {
        if frequency == 0 {
            continue;
        }

        let end = usize::from(range_start) + usize::from(frequency);
        let end = end.min(TOTAL_FREQUENCY);
        let packed = pack(symbol as u8, frequency, range_start);

        for slot in &mut slots[usize::from(range_start)..end] {
            *slot = packed;
        }

        range_start = range_start.saturating_add(frequency);
    }
}
