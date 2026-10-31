//! **A psp**, read by `call-from-psps` beside the fixture's second sample, and by
//! `pvc dev inspect-psp`. The third reader, `estimate-parameters`, is `psp_fit`: its fit takes
//! seconds a run, which would slow these two a hundredfold.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pop_var_caller_fuzz::{fixture, pvc, scratch_file, write_scratch};

fuzz_target!(|data: &[u8]| {
    let psp = write_scratch("sample.psp", data);
    let (reference, catalog) = (fixture("reference.fa"), fixture("catalog.parquet"));
    let (other, parameters) = (fixture("psps/HG003.psp"), fixture("parameters.toml"));
    let (vcf, log, table) =
        (scratch_file("psp.vcf"), scratch_file("psp.log"), scratch_file("psp.tsv"));
    pvc(&[
        "call-from-psps", "--reference", &reference, "--catalog", &catalog, "--psp", &psp,
        "--psp", &other, "--parameters", &parameters, "--threads", "1", "--quiet", "--force",
        "--output", &vcf, "--log", &log,
    ]);
    pvc(&[
        "dev", "inspect-psp", "--psp", &psp, "--regions", &fixture("regions.bed"), "--output",
        &table,
    ]);
});
