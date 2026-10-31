//! **A BAM.** The input is the decompressed BAM stream — magic, header, records — and the target
//! compresses it into BGZF blocks itself, so the fuzzer's mutations land on header text, record
//! fields and CIGARs rather than on deflate bytes that would only fail their checksum. The file is
//! called with the fixture's second sample, its index built from it by `--build-index-if-missing`.
#![no_main]

use std::io::Write;

use libfuzzer_sys::fuzz_target;
use pop_var_caller_fuzz::{CallFromAlignments, remove_scratch, scratch_file};

fuzz_target!(|data: &[u8]| {
    let path = scratch_file("sample.bam");
    let file = std::fs::File::create(&path).expect("the scratch BAM can be created");
    let mut writer = noodles_bgzf::io::Writer::new(file);
    writer.write_all(data).expect("the scratch BAM can be written");
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
