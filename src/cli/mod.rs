//! The `pop_var_caller` binary's command line: [`command_line`] owns the top-level `Parser` and the
//! subcommand enum, [`parsers`] the value parsers several subcommands share, and one module per
//! subcommand owns its `Args`, its `run_*` and its `#[non_exhaustive]` error enum.
//!
//! This was `pop_var_caller_exp`, ng's experiment binary, kept apart from the production
//! `pop_var_caller` binary so that ng's experiment knobs never grew it
//! (`doc/devel/ng/spec/typed_regions_cli.md` §2). Promotion step E2 renamed the binary and moved
//! this module from `src/pop_var_caller_exp/`; the subcommands kept their names.

pub mod call_from_alignments;
pub mod call_from_psps;
pub mod calling_run;
pub mod command_line;
mod cross_platform_digests;
pub mod estimate_contamination;
pub mod estimate_parameters;
pub mod generate_psps;
pub mod mode_equivalence;
pub mod parsers;
pub(crate) mod provenance;
pub mod psp_inputs;
pub mod regenerate_census;
pub mod repeat_catalog;
pub mod run_ground;
#[cfg(test)]
pub(crate) mod test_fixtures;
pub mod typed_regions;

pub use call_from_alignments::{
    CallFromAlignmentsArgs, CallFromAlignmentsCliError, run_call_from_alignments,
};
pub use call_from_psps::{CallFromPspsArgs, CallFromPspsCliError, run_call_from_psps};
pub use command_line::{Cli, PopVarCallerCommand};
pub use estimate_contamination::{
    EstimateContaminationArgs, EstimateContaminationCliError, run_estimate_contamination,
};
pub use estimate_parameters::{
    EstimateParametersArgs, EstimateParametersCliError, run_estimate_parameters,
};
pub use generate_psps::{GeneratePspsArgs, GeneratePspsCliError, run_generate_psps};
pub use regenerate_census::{
    CensusReport, RegenerateCensusArgs, RegenerateCensusCliError, SampleCensusOutcome, SkippedPsp,
    run_regenerate_census,
};
pub use repeat_catalog::{RepeatCatalogArgs, RepeatCatalogCliError, run_repeat_catalog};
pub use typed_regions::{TypedRegionsArgs, TypedRegionsCliError, run_typed_regions};
