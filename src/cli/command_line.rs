//! Top-level CLI for the `pvc` binary: the `Parser`, the subcommand enum, and `pvc dev`'s group.
//! Shape first copied from the deleted production binary's `cli` module.

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

/// The `pvc` command line: the top-level parser, whose help is the first thing a user reads.
#[derive(Debug, Parser)]
#[command(
    name = "pvc",
    version,
    about = "Call SNPs, indels and short tandem repeats (STRs) from aligned reads, for one sample \
             or a cohort, and write one VCF.",
    long_about = "Call SNPs, indels and short tandem repeats (STRs) from aligned reads, for one \
                  sample or a cohort, and write one VCF.

Two ways to run it:

  Direct mode — from alignments straight to a VCF, with the model's built-in numbers:
    pvc repeat-catalog --reference genome.fa
    pvc call-from-alignments --reference genome.fa --alignment s1.bam --alignment s2.bam \\
        --output calls.vcf.gz --defaults

  psp mode — store each sample's evidence once, fit the cohort's parameters, then call:
    pvc repeat-catalog --reference genome.fa
    pvc generate-psps --reference genome.fa --alignment s1.bam --output-dir psps/
    pvc estimate-parameters --reference genome.fa --psp psps/ --output cohort.parameters.toml
    pvc call-from-psps --reference genome.fa --psp psps/ --parameters cohort.parameters.toml \\
        --output calls.vcf.gz

A psp file holds what one sample's reads showed at every position walked, in the form the caller \
reads, so a cohort can be fitted and called — and a sample added later — without reading the \
alignments again."
)]
pub struct Cli {
    #[command(subcommand)]
    pub cmd: PopVarCallerCommand,

    /// Print only warnings and errors. The run's log, where a subcommand keeps one, is complete
    /// either way.
    #[arg(long, global = true)]
    pub quiet: bool,
}

/// The binary's subcommands, **in the order a user runs them**. Each kebab-cases to its command
/// name, as `RepeatCatalog` → `repeat-catalog`.
#[derive(Debug, Subcommand)]
pub enum PopVarCallerCommand {
    /// Build the reference's repeat catalogue, once per reference.
    ///
    /// Scans the reference for tandem repeats and writes the catalogue beside it, as
    /// `<reference>.repeats.parquet`. Every other command reads it from there.
    RepeatCatalog(RepeatCatalogArgs),

    /// Call samples from their BAM/CRAM files and write a VCF.
    ///
    /// Direct mode: every sample's alignments are read once and called together, with nothing
    /// stored in between. Without a fitted parameters file it calls with the numbers built into
    /// pvc (`--defaults`); fitted numbers come from psp mode, through `estimate-parameters`.
    CallFromAlignments(CallFromAlignmentsArgs),

    /// Store each sample's evidence as a psp file.
    ///
    /// Reads each sample's alignments once and writes one `<sample>.psp`, calling nothing.
    /// Samples given together are written one after another; to spread a cohort over cores or
    /// machines, run one command per sample, all with the same `--regions`.
    GeneratePsps(GeneratePspsArgs),

    /// Fit the model's parameters from psp files.
    ///
    /// Writes the parameters file `call-from-psps` scores with: base-quality calibration,
    /// contamination, inbreeding, population diversity and repeat stutter, fitted from the
    /// cohort itself.
    EstimateParameters(EstimateParametersArgs),

    /// Call samples from their psp files and write a VCF.
    ///
    /// psp mode's calling: the same calls as `call-from-alignments` on the same cohort and
    /// parameters, read from the psps instead of the alignments. The ground called and what
    /// counts as a repeat are the ones the psps were written with.
    CallFromPsps(CallFromPspsArgs),

    /// Rebuild the parameter-fitting data in psps written by older builds.
    RegenerateCensus(RegenerateCensusArgs),

    /// Tools for developing and debugging the caller.
    #[command(subcommand)]
    Dev(DevCommand),
}

/// `pvc dev`: the tools a developer of the caller uses, kept out of a user's way
/// (`doc/devel/ng/spec/command_line.md` §5.3).
#[derive(Debug, Subcommand)]
pub enum DevCommand {
    /// What one psp stored over a BED, as a table.
    ///
    /// For tracing whether a missing variant's reads reached the psp, with how many reads and
    /// from which read group.
    InspectPsp(InspectPspArgs),

    /// How the reference is divided into repeat and non-repeat ground, as a table.
    TypeRegions(TypedRegionsArgs),

    /// Contamination estimated on its own, for comparing methods.
    ///
    /// Each sample's share of reads from another individual, from a panel of alignments, written
    /// as JSON. The caller fits contamination itself in `estimate-parameters`; this exposes the
    /// same estimator alone.
    EstimateContamination(EstimateContaminationArgs),
}

/// **Every option of the subcommand `matches` ran, with the value the run used** — typed, taken
/// from its default, or neither — for the run's log
/// (`doc/devel/ng/spec/command_line.md` §7.3).
///
/// Read from the parsed command line rather than from each subcommand's `Args`, so a new option
/// is in the log without anyone remembering to add it.
#[must_use]
pub fn settings_of(matches: &ArgMatches) -> Vec<super::run_output::Setting> {
    // **Down to the subcommand that ran**, through `pvc dev`'s group.
    let mut command = Cli::command();
    let Some((mut name, mut ran)) = matches.subcommand() else {
        return Vec::new();
    };
    while let Some((inner, inner_matches)) = ran.subcommand() {
        let Some(group) = command.find_subcommand(name) else {
            return Vec::new();
        };
        command = group.clone();
        (name, ran) = (inner, inner_matches);
    }
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
                "--max-generic-alleles",
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

        let typed = named("max-generic-alleles");
        assert_eq!((typed.value.as_deref(), typed.typed), (Some("4"), true));
        let defaulted = named("max-depth");
        assert_eq!(
            (defaulted.value.as_deref(), defaulted.typed),
            (Some("1000"), false)
        );
        let not_given = named("max-reads-per-position");
        assert_eq!((not_given.value, not_given.typed), (None, false));
        assert!(settings.iter().all(|setting| setting.name != "help"));
    }

    /// **A `pvc dev` tool's options are found through the group**, not read as the group's own,
    /// which has none.
    #[test]
    fn a_dev_tools_options_are_recorded() {
        let matches = Cli::command()
            .try_get_matches_from([
                "pvc",
                "dev",
                "inspect-psp",
                "--psp",
                "a.psp",
                "--regions",
                "r.bed",
            ])
            .expect("the command line parses");
        let settings = settings_of(&matches);
        assert!(
            settings
                .iter()
                .any(|setting| setting.name == "psp" && setting.value.as_deref() == Some("a.psp")),
            "the tool's own option is recorded, and got: {settings:?}",
        );
    }
}
