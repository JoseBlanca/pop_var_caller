//! Top-level CLI for the `pop_var_caller_exp` binary: the `Parser` and the
//! subcommand enum. Shape copied from [`crate::pop_var_caller::cli`].

use clap::{Parser, Subcommand};

use super::call_from_alignments::CallFromAlignmentsArgs;
use super::call_from_psps::CallFromPspsArgs;
use super::estimate_contamination::EstimateContaminationArgs;
use super::generate_psps::GeneratePspsArgs;
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

    /// Call a cohort of alignment files and write a VCF, in one process.
    ///
    /// Every sample's file is held open for the whole run, every sample's reads
    /// are walked at one shared frontier, and each locus is genotyped where it
    /// is built. Nothing is written between the alignments and the VCF.
    ///
    /// One limit to know before you read the output: the run decodes every
    /// sample's reads across the machine's cores, but assembles and genotypes
    /// on one thread, and says nothing while it works.
    ///
    /// Repeat tracts ARE called, through their own stutter model, and their
    /// records carry STR, RU, PERIOD and each called allele's REPCN. What is
    /// still set aside is a repeat cluster too tangled to have clean flanks;
    /// its ground is charged to `not built yet` and reported at the end.
    CallFromAlignments(CallFromAlignmentsArgs),

    /// Walk each sample's alignment files once and store what they showed as a psp,
    /// one file per sample, calling nothing.
    ///
    /// A psp holds one sample's evidence: what its reads showed at every position of the
    /// ground you asked for — the alleles, their support, which reads carried them — in
    /// the form the caller reads. This is psp mode's first stage, and a sample walked here
    /// can join any later cohort without being read again: adding a sample re-walks that
    /// sample only, and a failed sample is one sample to re-run.
    ///
    /// Samples are walked one at a time, in the order given. There is no thread knob:
    /// each sample's walk is independent, so a cohort is spread by running this command
    /// once per sample rather than by threading one invocation.
    ///
    /// Each psp is named for the sample its reads declare, inside --output-dir. Two things
    /// this stage does not do yet: it writes no census beside the psp (the file a
    /// parameters fit reads), and it does not refuse to overwrite a psp already there.
    GeneratePsps(GeneratePspsArgs),

    /// Call a cohort of stored psps and write a VCF — psp mode's second stage.
    ///
    /// Same cohort, same parameters, same VCF as `call-from-alignments`: what differs is
    /// that the evidence is read from files a `generate-psps` run already wrote, so no
    /// alignment file is opened and no read is decoded again.
    ///
    /// There is no --regions here. A psp records the ground its walk covered, the cohort is
    /// refused unless every file agrees about it, and that agreed ground is what this run
    /// calls over. To call over less, walk less.
    ///
    /// What it can say about each sample is narrower than direct mode's, and for a reason no
    /// flag can fix: a psp carries no count of the reads its walk kept or dropped. So the
    /// report says how many stored loci this run read out of each file and how deep they
    /// were — and it names any file whose walk applied different read filters from the rest,
    /// which nothing else in the pipeline checks.
    CallFromPsps(CallFromPspsArgs),

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
