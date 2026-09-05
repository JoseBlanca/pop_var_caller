//! Top-level CLI for the `pop_var_caller_exp` binary: the `Parser` and the
//! subcommand enum. Shape copied from [`crate::pop_var_caller::cli`].

use clap::{Parser, Subcommand};

use super::call_from_alignments::CallFromAlignmentsArgs;
use super::call_from_psps::CallFromPspsArgs;
use super::estimate_contamination::EstimateContaminationArgs;
use super::estimate_parameters::EstimateParametersArgs;
use super::generate_census::GenerateCensusArgs;
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
    /// Each psp is named for the sample its reads declare, inside --output-dir, and beside
    /// it goes that sample's census — the smaller file a parameters fit reads — from the
    /// same single pass over the reads. A psp already in --output-dir is refused before
    /// anything is walked; --force replaces it.
    GeneratePsps(GeneratePspsArgs),

    /// Build each stored psp's census, without re-reading a single alignment file.
    ///
    /// A census is the small file a parameters fit reads: what one sample showed at a
    /// fixed set of positions and repeat tracts chosen for the whole run, so the fit can
    /// ask the same question of every sample and compare their answers.
    ///
    /// `generate-psps` already writes one beside each psp. This is for the cases that
    /// cannot re-walk the reads: psps written before censuses existed, a census lost or
    /// built under settings since changed, and a census wanted larger than the one on
    /// disk.
    ///
    /// There is no --regions, for the reason call-from-psps has none. The psps record the
    /// ground they were walked over and the cohort is refused unless they agree about it;
    /// that agreed ground is what the positions are chosen from. Choosing them over other
    /// ground would produce censuses the cohort cannot be fitted from.
    GenerateCensus(GenerateCensusArgs),

    /// Fit a cohort's parameters from its censuses and write them as a parameters file.
    ///
    /// This is the file a calling run scores with. Without one a run has two choices and
    /// neither is a fit: the constants compiled into the binary, or a file somebody hands
    /// it.
    ///
    /// It reads the censuses. It does not read the psps — but each census's psp must be
    /// beside it, because a census names the psp it was built from and evidence from other
    /// reads is otherwise indistinguishable from this run's. What is taken from each psp is
    /// its header: one short read.
    ///
    /// The reference and the catalog are needed because a census stores a repeat tract by
    /// its index within its stratum and nothing else, so the selection has to be rebuilt —
    /// and the rebuild is checked against a digest every census carries.
    ///
    /// The inbreeding coefficient is declared rather than fitted here; --inbreeding states
    /// one for the cohort and the file records it as supplied.
    EstimateParameters(EstimateParametersArgs),

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
