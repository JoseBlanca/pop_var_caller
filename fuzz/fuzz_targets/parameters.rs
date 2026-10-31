//! **A parameters file** (`--parameters`), called over the fixture. The seed is the file
//! `estimate-parameters` fitted from the fixture's own psps, so its reference digest and samples
//! match and a mutated copy is read all the way to the calling.
#![no_main]

use libfuzzer_sys::fuzz_target;
use pop_var_caller_fuzz::{CallFromAlignments, write_scratch};

fuzz_target!(|data: &[u8]| {
    let parameters = write_scratch("parameters.toml", data);
    CallFromAlignments { parameters: Some(parameters), ..Default::default() }.run();
});
