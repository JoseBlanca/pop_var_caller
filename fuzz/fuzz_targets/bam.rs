//! **A BAM.** The input is the decompressed BAM stream — magic, header, records — and the target
//! compresses it into BGZF blocks itself, so the fuzzer's mutations land on header text, record
//! fields and CIGARs rather than on deflate bytes that would only fail their checksum. The file is
//! called with the fixture's second sample, its index built from it by `--build-index-if-missing`.
//!
//! **Inputs declaring a size no real file has are skipped**, and that is a known gap, not
//! coverage. noodles-bam allocates what a BAM's reference-sequence count, a reference name's
//! length, or a record's size declares before reading the bytes (0.89.0, and 0.98.0 still does
//! for the last two), so such an input makes the process ask for gigabytes. Without a memory
//! cap that ends in a clean error; under one it aborts (fuzz/README.md). Every such input would
//! re-find that one behaviour and stop the campaign, so the filter lets the fuzzer reach the
//! code past it.
#![no_main]

use std::io::Write;

use libfuzzer_sys::fuzz_target;
use pop_var_caller_fuzz::{CallFromAlignments, remove_scratch, scratch_file};

/// More reference sequences than any assembly has: GRCh38 with decoys and alts is about 3,400.
const MOST_REFERENCE_SEQUENCES: u32 = 1 << 20;
/// A reference name or a record longer than this is not from a real file.
const LONGEST_FIELD: u32 = 1 << 26;

fn u32_at(data: &[u8], at: usize) -> Option<u32> {
    Some(u32::from_le_bytes(data.get(at..at + 4)?.try_into().ok()?))
}

/// Whether `data` declares a size past the bounds above, anywhere a reader would believe it.
fn declares_an_impossible_size(data: &[u8]) -> bool {
    let Some(text) = u32_at(data, 4) else {
        return false;
    };
    if text > LONGEST_FIELD {
        return true;
    }
    let mut at = 8 + text as usize;
    let Some(references) = u32_at(data, at) else {
        return false;
    };
    if references > MOST_REFERENCE_SEQUENCES {
        return true;
    }
    at += 4;
    for _ in 0..references {
        let Some(name) = u32_at(data, at) else {
            return false;
        };
        if name > LONGEST_FIELD {
            return true;
        }
        at += 4 + name as usize + 4;
    }
    while let Some(record) = u32_at(data, at) {
        if record > LONGEST_FIELD {
            return true;
        }
        at += 4 + record as usize;
    }
    false
}

fuzz_target!(|data: &[u8]| {
    if declares_an_impossible_size(data) {
        return;
    }
    let path = scratch_file("sample.bam");
    let file = std::fs::File::create(&path).expect("the scratch BAM can be created");
    let mut writer = noodles_bgzf::io::Writer::new(file);
    writer
        .write_all(data)
        .expect("the scratch BAM can be written");
    writer.finish().expect("the scratch BAM can be closed");
    remove_scratch("sample.bam.bai");
    remove_scratch("sample.bam.csi");
    CallFromAlignments {
        alignments: vec![path, pop_var_caller_fuzz::fixture("HG003.bam")],
        extra: &["--build-index-if-missing"],
        ..Default::default()
    }
    .run();
});
