//! Top-level CLI for the `pop_var_caller_exp` binary: the `Parser` and the
//! subcommand enum. Shape copied from [`crate::pop_var_caller::cli`].

use clap::{Parser, Subcommand};

use super::estimate_contamination::EstimateContaminationArgs;
use super::repeat_catalog::RepeatCatalogArgs;
use super::typed_regions::TypedRegionsArgs;

pub mod parsers;

/// Top-level CLI for the `pop_var_caller_exp` binary.
#[derive(Debug, Parser)]
#[command(name = "pop_var_caller_exp", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub cmd: PopVarCallerExpCommand,
}

/// The exp binary's subcommands. Each kebab-cases to its command name, as
/// `SsrCatalog` → `ssr-catalog`.
#[derive(Debug, Subcommand)]
pub enum PopVarCallerExpCommand {
    /// Run step 3's walk over a reference and write the typed-region
    /// partition to a file (contig, span, kind, and STR detail per region).
    TypeRegions(TypedRegionsArgs),

    /// Scan a reference for tandem repeats once and write the catalog beside
    /// it, so that every later run reads the file instead of re-scanning the
    /// genome (doc/devel/ng/spec/repeat_catalog.md).
    RepeatCatalog(RepeatCatalogArgs),

    /// Estimate, for each sample in a panel of alignments, what share of its
    /// reads came from another individual, and write the answers as JSON.
    ///
    /// **A side tool, not a step of the caller.** ng fits contamination
    /// internally before calling; this exposes the same estimator to somebody
    /// comparing methods. It needs a panel — about a dozen samples — because
    /// the allele frequencies it judges each sample against are fitted from
    /// the run itself and no outside panel.
    EstimateContamination(EstimateContaminationArgs),
}
