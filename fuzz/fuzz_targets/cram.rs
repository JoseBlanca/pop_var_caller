//! **A CRAM**, as raw bytes: CRAM's own containers and blocks are what the reader parses, and only
//! some of them are compressed. Called against the fixture's reference with the second sample,
//! its index built by `--build-index-if-missing`.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pop_var_caller_fuzz::{CallFromAlignments, remove_scratch, write_scratch};

fuzz_target!(|data: &[u8]| {
    let path = write_scratch("sample.cram", data);
    remove_scratch("sample.cram.crai");
    CallFromAlignments {
        alignments: vec![path, pop_var_caller_fuzz::fixture("HG003.bam")],
        extra: &["--build-index-if-missing"],
        ..Default::default()
    }
    .run();
});
