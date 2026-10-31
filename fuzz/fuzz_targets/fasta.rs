//! **A reference FASTA and its `.fai`.** The input is the two files joined by a line holding only
//! `@@FAI`; without that line there is no index, and pvc writes one. The reference is read twice,
//! the two ways a user's reference is: by `repeat-catalog`, which scans every base, and by
//! `call-from-alignments`, which checks it against the alignments' headers and then reads it a
//! window at a time, here with the catalog just built from it.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pop_var_caller_fuzz::{CallFromAlignments, pvc, remove_scratch, scratch_file, write_scratch};

const SEPARATOR: &[u8] = b"\n@@FAI\n";

fuzz_target!(|data: &[u8]| {
    let split = data.windows(SEPARATOR.len()).position(|window| window == SEPARATOR);
    let (fasta, fai) = match split {
        Some(at) => (&data[..=at], Some(&data[at + SEPARATOR.len()..])),
        None => (data, None),
    };
    let reference = write_scratch("reference.fa", fasta);
    match fai {
        Some(fai) => {
            write_scratch("reference.fa.fai", fai);
        }
        None => remove_scratch("reference.fa.fai"),
    }
    let catalog = scratch_file("catalog.parquet");
    let log = scratch_file("catalog.log");
    remove_scratch("catalog.parquet");
    pvc(&[
        "repeat-catalog", "--reference", &reference, "--output", &catalog, "--log", &log,
        "--force", "--quiet",
    ]);
    if std::path::Path::new(&catalog).exists() {
        CallFromAlignments { reference, catalog, ..Default::default() }.run();
    }
});
