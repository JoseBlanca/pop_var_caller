//! `pvc` binary entry point. Parses the top-level command line and dispatches to the
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

use clap::{CommandFactory, FromArgMatches};
use pop_var_caller::cli::command_line::settings_of;
use pop_var_caller::cli::run_output;
use pop_var_caller::cli::{
    Cli, DevCommand, PopVarCallerCommand, run_call_from_alignments, run_call_from_psps,
    run_estimate_contamination, run_estimate_parameters, run_generate_psps, run_inspect_psp,
    run_repeat_catalog, run_typed_regions,
};
use pop_var_caller::error_render::format_error_chain;

fn main() {
    let matches = Cli::command().get_matches();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|error| error.exit());
    run_output::install(cli.quiet, settings_of(&matches));
    // **A panic is recorded in the log before the process stops**: a release build aborts on a
    // panic, so nothing after it runs. Rust's own message still goes to the terminal.
    let print_the_panic = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        print_the_panic(info);
        run_output::close_log_after_a_panic(&format!("pvc stopped on an internal error: {info}"));
    }));
    let result = match cli.cmd {
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
        PopVarCallerCommand::EstimateParameters(args) => {
            run_estimate_parameters(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::Dev(DevCommand::InspectPsp(args)) => {
            run_inspect_psp(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::Dev(DevCommand::TypeRegions(args)) => {
            run_typed_regions(&args).map_err(|e| format_error_chain(&e))
        }
        PopVarCallerCommand::Dev(DevCommand::EstimateContamination(args)) => {
            run_estimate_contamination(&args).map_err(|e| format_error_chain(&e))
        }
    };
    // **The error is the run's last word, so it goes to the log too**, which is why the log is
    // closed here and not when the subcommand returns.
    if let Err(msg) = result {
        tracing::error!("{msg}");
        run_output::close_log(1);
        process::exit(1);
    }
    run_output::close_log(0);
}
