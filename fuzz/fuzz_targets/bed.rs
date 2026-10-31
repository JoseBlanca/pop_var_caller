//! **A BED of regions** (`--regions`), called over the fixture.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pop_var_caller_fuzz::{CallFromAlignments, write_scratch};

fuzz_target!(|data: &[u8]| {
    let regions = write_scratch("regions.bed", data);
    CallFromAlignments { regions, ..Default::default() }.run();
});
