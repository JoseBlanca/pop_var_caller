//! `pop_var_caller_exp` binary namespace — ng's experiment command
//! surface. A second binary in the same crate, kept apart from the
//! production `pop_var_caller` CLI so ng's experiment knobs never grow it
//! (see `doc/devel/ng/spec/typed_regions_cli.md` §2). The library both
//! binaries link is the same one; only what a user can *invoke* is split.
//!
//! Layout mirrors [`crate::pop_var_caller`]: [`cli`] owns the top-level
//! `Parser` plus the subcommand enum, and one module per subcommand owns
//! its `Args`, its `run_*`, and its `#[non_exhaustive]` error enum.

pub mod call_from_alignments;
pub mod call_from_psps;
pub mod calling_run;
pub mod cli;
pub mod estimate_contamination;
pub mod estimate_parameters;
pub mod generate_census;
pub mod generate_psps;
pub mod mode_equivalence;
pub mod repeat_catalog;
pub mod run_ground;
#[cfg(test)]
pub(crate) mod test_fixtures;
pub mod typed_regions;

pub use call_from_alignments::{
    CallFromAlignmentsArgs, CallFromAlignmentsCliError, run_call_from_alignments,
};
pub use call_from_psps::{CallFromPspsArgs, CallFromPspsCliError, run_call_from_psps};
pub use cli::{Cli, PopVarCallerExpCommand};
pub use estimate_contamination::{
    EstimateContaminationArgs, EstimateContaminationCliError, run_estimate_contamination,
};
pub use estimate_parameters::{
    EstimateParametersArgs, EstimateParametersCliError, run_estimate_parameters,
};
pub use generate_census::{
    CensusReport, GenerateCensusArgs, GenerateCensusCliError, SampleCensusOutcome,
    run_generate_census,
};
pub use generate_psps::{GeneratePspsArgs, GeneratePspsCliError, run_generate_psps};
pub use repeat_catalog::{RepeatCatalogArgs, RepeatCatalogCliError, run_repeat_catalog};
pub use typed_regions::{TypedRegionsArgs, TypedRegionsCliError, run_typed_regions};
