//! `pop_var_caller` binary entry point. Parses the top-level command line and dispatches to the
//! subcommand's driver; the logic lives in `src/cli/`, so this file is intentionally thin. Errors
//! are rendered through the library's `format_error_chain` (spec T7a). This was the
//! `pop_var_caller_exp` binary, `src/main_exp.rs`, until promotion step E2 gave it the name the
//! deleted production binary had.

use std::process;

// The `mimalloc` global allocator. A `#[global_allocator]` is per *binary*, not per crate, so
// the `alloc-mimalloc` default feature does nothing for a binary that does not
// declare one: without this line every `call-from-alignments` run — and every
// number measured from one — used the system allocator while every probe in
// `examples/` used mimalloc.
//
// It is worth its line here: a calling run frees far more blocks than it
// allocates on the merge thread, because the observations it walks were
// allocated by the sample sweeps and released as it passes them, and a
// system allocator takes a lock per cross-thread free.
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

use clap::Parser;
use pop_var_caller::cli::{
    Cli, PopVarCallerCommand, run_call_from_alignments, run_call_from_psps,
    run_estimate_contamination, run_estimate_parameters, run_generate_psps, run_regenerate_census,
    run_repeat_catalog, run_typed_regions,
};
use pop_var_caller::error_render::format_error_chain;

fn main() {
    let cli = Cli::parse();
    let result = match cli.cmd {
        PopVarCallerCommand::TypeRegions(args) => {
            run_typed_regions(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::RepeatCatalog(args) => {
            run_repeat_catalog(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::CallFromAlignments(args) => {
            run_call_from_alignments(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::CallFromPsps(args) => {
            run_call_from_psps(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::GeneratePsps(args) => {
            run_generate_psps(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::RegenerateCensus(args) => {
            run_regenerate_census(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::EstimateParameters(args) => {
            run_estimate_parameters(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::EstimateContamination(args) => {
            run_estimate_contamination(&args).map_err(|e| format_error_chain(&e))
        }
    };
    if let Err(msg) = result {
        eprintln!("error: {msg}");
        process::exit(1);
    }
}
