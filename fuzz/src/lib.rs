//! **What every fuzz target shares**: the fixture, a scratch directory, and a way to run one of
//! `pvc`'s subcommands in-process.
//!
//! Each target writes the fuzzer's bytes over one of the fixture's files and runs the subcommands
//! that read that file, exactly as `pvc` would, with the rest of the fixture as the other inputs.
//! **An error is the right answer to malformed input; a panic or an abort is the finding.** So
//! the harness throws every `Err` away and lets a panic reach libFuzzer, which records the input.

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use clap::Parser;
use pop_var_caller::cli::{
    Cli, DevCommand, PopVarCallerCommand, run_call_from_alignments, run_call_from_psps,
    run_estimate_contamination, run_estimate_parameters, run_generate_psps, run_inspect_psp,
    run_repeat_catalog, run_typed_regions,
};

/// A file of the fixture `make_fixture.sh` built.
pub fn fixture(name: &str) -> String {
    format!("{}/fixture/{name}", env!("CARGO_MANIFEST_DIR"))
}

/// This process's scratch directory, under the repository's `tmp/` (never the system `/tmp`), so
/// that fuzzers running side by side do not overwrite each other's inputs.
pub fn scratch() -> &'static Path {
    static SCRATCH: OnceLock<PathBuf> = OnceLock::new();
    SCRATCH.get_or_init(|| {
        let directory = PathBuf::from(format!(
            "{}/../tmp/fuzz-scratch/{}",
            env!("CARGO_MANIFEST_DIR"),
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).expect("the scratch directory can be made");
        directory
    })
}

/// A path in this process's scratch directory.
pub fn scratch_file(name: &str) -> String {
    scratch().join(name).display().to_string()
}

/// Write `bytes` to the scratch file `name` and return its path.
pub fn write_scratch(name: &str, bytes: &[u8]) -> String {
    let path = scratch_file(name);
    std::fs::write(&path, bytes).expect("the scratch file can be written");
    path
}

/// Remove a scratch file if it is there: an index left by the last input would otherwise be read
/// against this one.
pub fn remove_scratch(name: &str) {
    let _ = std::fs::remove_file(scratch().join(name));
}

/// Run `pvc <arguments>` in-process and discard its result. The arguments are the command line
/// without the program name; one that clap refuses is a harness mistake, so it panics.
pub fn pvc(arguments: &[&str]) {
    let cli = Cli::try_parse_from(std::iter::once("pvc").chain(arguments.iter().copied()))
        .unwrap_or_else(|error| panic!("the harness built a command pvc refuses: {error}"));
    let _ = match cli.cmd {
        PopVarCallerCommand::RepeatCatalog(args) => run_repeat_catalog(&args).map_err(drop),
        PopVarCallerCommand::CallFromAlignments(args) => {
            run_call_from_alignments(&args).map_err(drop)
        }
        PopVarCallerCommand::CallFromPsps(args) => run_call_from_psps(&args).map_err(drop),
        PopVarCallerCommand::GeneratePsps(args) => run_generate_psps(&args).map_err(drop),
        PopVarCallerCommand::EstimateParameters(args) => {
            run_estimate_parameters(&args).map_err(drop)
        }
        PopVarCallerCommand::Dev(DevCommand::InspectPsp(args)) => {
            run_inspect_psp(&args).map_err(drop)
        }
        PopVarCallerCommand::Dev(DevCommand::TypeRegions(args)) => {
            run_typed_regions(&args).map_err(drop)
        }
        PopVarCallerCommand::Dev(DevCommand::EstimateContamination(args)) => {
            run_estimate_contamination(&args).map_err(drop)
        }
    };
}

/// `pvc call-from-alignments` over the fixture, with any of its inputs replaced. One thread and
/// `--quiet`, so a run is as cheap as it can be.
pub struct CallFromAlignments<'a> {
    pub reference: String,
    pub catalog: String,
    pub alignments: Vec<String>,
    pub regions: String,
    /// `None` runs with `--defaults`.
    pub parameters: Option<String>,
    pub extra: &'a [&'a str],
}

impl Default for CallFromAlignments<'_> {
    fn default() -> Self {
        Self {
            reference: fixture("reference.fa"),
            catalog: fixture("catalog.parquet"),
            alignments: vec![fixture("HG002.bam"), fixture("HG003.bam")],
            regions: fixture("regions.bed"),
            parameters: None,
            extra: &[],
        }
    }
}

impl CallFromAlignments<'_> {
    pub fn run(&self) {
        let output = scratch_file("calls.vcf");
        let log = scratch_file("calls.log");
        let mut arguments = vec![
            "call-from-alignments",
            "--reference",
            &self.reference,
            "--catalog",
            &self.catalog,
            "--regions",
            &self.regions,
            "--threads",
            "1",
            "--quiet",
            "--force",
            "--output",
            &output,
            "--log",
            &log,
        ];
        for alignment in &self.alignments {
            arguments.extend(["--alignment", alignment.as_str()]);
        }
        match &self.parameters {
            Some(parameters) => arguments.extend(["--parameters", parameters.as_str()]),
            None => arguments.push("--defaults"),
        }
        arguments.extend(self.extra);
        pvc(&arguments);
    }
}
