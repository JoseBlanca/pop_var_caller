//! **A psp read by `estimate-parameters`**, beside the fixture's second sample. Slow — the fit
//! takes seconds a run — so run it with a long `-timeout`; its own paths are the census sections
//! of the psp's trailer, which only the fit reads.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pop_var_caller_fuzz::{fixture, pvc, scratch_file, write_scratch};

fuzz_target!(|data: &[u8]| {
    let psp = write_scratch("sample.psp", data);
    let (reference, catalog) = (fixture("reference.fa"), fixture("catalog.parquet"));
    let other = fixture("psps/HG003.psp");
    let (toml, log) = (scratch_file("fit.toml"), scratch_file("fit.log"));
    pvc(&[
        "estimate-parameters", "--reference", &reference, "--catalog", &catalog, "--psp", &psp,
        "--psp", &other, "--threads", "1", "--quiet", "--force", "--output", &toml, "--log",
        &log,
    ]);
});
