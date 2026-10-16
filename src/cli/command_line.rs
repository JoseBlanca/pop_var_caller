//! Top-level CLI for the `pvc` binary: the `Parser` and the
//! subcommand enum. Shape first copied from the deleted production binary's `cli` module.

use clap::parser::ValueSource;
use clap::{ArgMatches, CommandFactory, Parser, Subcommand};

use super::call_from_alignments::CallFromAlignmentsArgs;
use super::call_from_psps::CallFromPspsArgs;
use super::estimate_contamination::EstimateContaminationArgs;
use super::estimate_parameters::EstimateParametersArgs;
use super::generate_psps::GeneratePspsArgs;
use super::inspect_psp::InspectPspArgs;
use super::regenerate_census::RegenerateCensusArgs;
use super::repeat_catalog::RepeatCatalogArgs;
use super::typed_regions::TypedRegionsArgs;

/// Top-level CLI for the `pvc` binary.
#[derive(Debug, Parser)]
#[command(name = "pvc", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub cmd: PopVarCallerCommand,

    /// Print only warnings and errors. The run's log, where a subcommand keeps one, is complete
    /// either way.
    #[arg(long, global = true)]
    pub quiet: bool,
}

/// **Every option of the subcommand `matches` ran, with the value the run used** — typed, taken
/// from its default, or neither — for the run's log
/// (`doc/devel/ng/spec/command_line.md` §7.3).
///
/// Read from the parsed command line rather than from each subcommand's `Args`, so a new option
/// is in the log without anyone remembering to add it.
#[must_use]
pub fn settings_of(matches: &ArgMatches) -> Vec<super::run_output::Setting> {
    let Some((name, ran)) = matches.subcommand() else {
        return Vec::new();
    };
    let command = Cli::command();
    let Some(subcommand) = command.find_subcommand(name) else {
        return Vec::new();
    };
    subcommand
        .get_arguments()
        .filter(|argument| !matches!(argument.get_id().as_str(), "help" | "version"))
        .map(|argument| {
            let id = argument.get_id().as_str();
            let source = ran.value_source(id);
            let value = source.and_then(|_| {
                ran.get_raw(id).map(|values| {
                    values
                        .map(|value| value.to_string_lossy().into_owned())
                        .collect::<Vec<_>>()
                        .join(" ")
                })
            });
            super::run_output::Setting {
                name: argument.get_long().unwrap_or(id).to_owned(),
                value,
                typed: source == Some(ValueSource::CommandLine),
            }
        })
        .collect()
}

/// The binary's subcommands. Each kebab-cases to its command name, as
/// `SsrCatalog` → `ssr-catalog`.
#[derive(Debug, Subcommand)]
pub enum PopVarCallerCommand {
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
    /// records carry STR, RU, PERIOD and each called allele's REPCN. An STR
    /// bundle (repeats too close together for any to have clean flanks) is
    /// called as one locus over its span, with stutter in each of its runs:
    /// its record carries one RU and PERIOD value per run and, per sample,
    /// each called allele's run counts in REPCNB (none of the four when the
    /// bundle holds no run).
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
    /// Each psp is named for the sample its reads declare, inside --output-dir. A sample is
    /// one file: the census a parameters fit reads is sealed into that psp's tail, from the
    /// same single pass over the reads, so there is nothing to keep beside it and nothing to
    /// lose. A psp already in --output-dir is refused before anything is walked; --force
    /// replaces it.
    GeneratePsps(GeneratePspsArgs),

    /// Rebuild each psp's census from the records it already holds, without re-reading a
    /// single alignment file.
    ///
    /// A census is what a parameters fit reads: what one sample showed at a fixed set of
    /// positions and repeat tracts chosen for the whole run, so the fit can ask the same
    /// question of every sample and compare their answers. It lives in the psp's trailer,
    /// written there by the walk.
    ///
    /// This is the repair for the psps a fit refuses — one written before the census moved
    /// into the psp, one whose tail is not a census at all, one written by another build, and
    /// one recorded under settings the run does not use.
    /// estimate-parameters names those samples and tells you to run this. Each named psp's
    /// trailer is replaced; its header, blocks and index are the bytes they were.
    ///
    /// There is no flag here for what counts as a repeat, and no --regions: the ground the
    /// walk covered and the criteria it cut that ground with are in every psp's header, and
    /// the cohort is refused unless the files agree about them. The reference and the
    /// catalog are needed to choose the positions again, and both are checked against the
    /// psps before anything is written — a census rebuilt against another reference is one
    /// every fit refuses.
    RegenerateCensus(RegenerateCensusArgs),

    /// Fit a cohort's parameters from its psps and write them as a parameters file.
    ///
    /// This is the file a calling run scores with. Without one a run has two choices and
    /// neither is a fit: the constants compiled into the binary, or a file somebody hands
    /// it.
    ///
    /// It reads each psp's census — the small object in the file's tail — and none of its
    /// records. From each header it takes the ground the walk covered and the repeat
    /// criteria it cut that ground with, so there is no flag here that says what a repeat
    /// is: the psps do. A cohort whose files disagree about either is refused.
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

    /// Dump what one psp file stored over the regions of a BED, as a TSV with one row per
    /// locus and observation.
    ///
    /// For tracing a missing variant: whether its reads reached the psp, with how many reads,
    /// from which read group, and whether they covered the whole locus. Everything a calling
    /// run decides afterwards is not here. The column names are stable, and the header names
    /// the psp format version and the column layout version.
    InspectPsp(InspectPspArgs),
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The log tells a typed value from a default, and records an option given neither** — the
    /// three cases a person reading it needs apart to reproduce the run.
    #[test]
    fn every_option_is_recorded_typed_defaulted_or_not_given() {
        let matches = Cli::command()
            .try_get_matches_from([
                "pvc",
                "call-from-alignments",
                "--reference",
                "ref.fa",
                "--alignment",
                "one.cram",
                "--output",
                "calls.vcf",
                "--defaults",
                "--max-candidate-alleles",
                "4",
            ])
            .expect("the command line parses");
        let settings = settings_of(&matches);
        let named = |name: &str| {
            settings
                .iter()
                .find(|setting| setting.name == name)
                .unwrap_or_else(|| panic!("{name} is recorded"))
                .clone()
        };

        let typed = named("max-candidate-alleles");
        assert_eq!((typed.value.as_deref(), typed.typed), (Some("4"), true));
        let defaulted = named("max-reads-per-position");
        assert_eq!(
            (defaulted.value.as_deref(), defaulted.typed),
            (Some("1000"), false)
        );
        let not_given = named("regions");
        assert_eq!((not_given.value, not_given.typed), (None, false));
        assert!(settings.iter().all(|setting| setting.name != "help"));
    }
}
