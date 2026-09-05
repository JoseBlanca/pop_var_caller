//! `pop_var_caller_exp` binary entry point — ng's experiment command
//! surface. Parses the top-level CLI and dispatches to the subcommand
//! driver; the logic lives in `pop_var_caller_exp/`, so this file is
//! intentionally thin, the way `src/main.rs` is. Both binaries render
//! errors through the shared `format_error_chain` (spec T7a).

use std::process;

// The `mimalloc` global allocator, the same one `src/main.rs` installs and for
// the same reason. A `#[global_allocator]` is per *binary*, not per crate, so
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
use pop_var_caller::error_render::format_error_chain;
use pop_var_caller::pop_var_caller_exp::{
    Cli, PopVarCallerExpCommand, run_call_from_alignments, run_call_from_psps,
    run_estimate_contamination, run_estimate_parameters, run_generate_census, run_generate_psps,
    run_repeat_catalog, run_typed_regions,
};

fn main() {
    let cli = Cli::parse();
    let result = match cli.cmd {
        PopVarCallerExpCommand::TypeRegions(args) => {
            run_typed_regions(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerExpCommand::RepeatCatalog(args) => {
            run_repeat_catalog(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerExpCommand::CallFromAlignments(args) => {
            run_call_from_alignments(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerExpCommand::CallFromPsps(args) => {
            run_call_from_psps(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerExpCommand::GeneratePsps(args) => {
            run_generate_psps(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerExpCommand::GenerateCensus(args) => {
            run_generate_census(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerExpCommand::EstimateParameters(args) => {
            run_estimate_parameters(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerExpCommand::EstimateContamination(args) => {
            run_estimate_contamination(&args).map_err(|e| format_error_chain(&e))
        }
    };
    if let Err(msg) = result {
        eprintln!("error: {msg}");
        process::exit(1);
    }
}
