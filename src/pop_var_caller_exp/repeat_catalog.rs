//! The `repeat-catalog` subcommand: scan a reference for tandem repeats once and write
//! what was found beside it.
//!
//! **This is the only place the scan happens** (spec §2.5): every other run reads the file
//! this writes. The scan rides the reference pass that already computes the digests, so a
//! build costs one read of the FASTA and not two.

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::ng::reference_info::{
    ReferenceInfoError, ReferenceSource, read_reference_info_observing,
};
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_BUNDLE_THRESHOLD, DEFAULT_MIN_PURITY, DEFAULT_MIN_SCORE, MAX_MOTIF_LEN, MinCopies,
    SsrSegmentCriteria,
};
use crate::ng::repeat_catalog::criteria::{
    CATALOG_MAX_PERIOD, CATALOG_MAX_STR_LEN_BP, CATALOG_MIN_FLANK_BP, CATALOG_MIN_PERIOD,
};
use crate::ng::repeat_catalog::{
    BuildTally, RepeatCatalogBuilder, RepeatCatalogError, StrRepeatCriteria,
};
use crate::ng::tandem_repeat::{
    DEFAULT_MATCH_REWARD, DEFAULT_MISMATCH_PENALTY, PeriodRange, PeriodRangeError, ScanParams,
};
use crate::ng::types::Bp;
use crate::pop_var_caller_exp::cli::parsers::parse_min_copies;

/// What the file is called when `--output` is not given: a sibling of the reference, the
/// way a `.fai` is, so a later run finds it without being told where it is.
pub const CATALOG_SUFFIX: &str = ".repeats.parquet";

/// The permissive floor handed to the scanner itself.
///
/// The per-period floors are applied when a row is written, on the **detected** span; the
/// scanner is given the table's minimum so nothing a reader could ask for is dropped before
/// it is seen (spec §4.1).
const SCANNER_MIN_COPIES: u32 = 2;

/// `repeat-catalog` arguments.
///
/// The defaults are the **catalog's**, not step 3's calling defaults: this file is built to
/// outlive any one calling policy, so its floors sit below every one of them (spec §4.1).
#[derive(Debug, Args, Clone)]
pub struct RepeatCatalogArgs {
    /// Reference FASTA. Read from byte zero — the digests and the scan come from the same
    /// pass.
    #[arg(long)]
    pub reference: PathBuf,

    /// Where to write the catalog. Defaults to `<reference>.repeats.parquet`.
    #[arg(long)]
    pub output: Option<PathBuf>,

    /// How many contigs to scan at once. A speed knob only: the file is byte-identical at
    /// every value.
    #[arg(long, default_value_t = 1)]
    pub threads: usize,

    /// Overwrite an existing catalog.
    ///
    /// Without it an existing file is left alone: a run reading a catalog while another
    /// rebuilds it in place fails with no error message, so replacing one is said out loud.
    #[arg(long)]
    pub force: bool,

    /// Smallest motif length scanned, in bases.
    #[arg(long, default_value_t = CATALOG_MIN_PERIOD)]
    pub min_period: u8,

    /// Largest motif length scanned, in bases.
    #[arg(long, default_value_t = CATALOG_MAX_PERIOD)]
    pub max_period: u8,

    /// The fewest motif copies a tract needs to be written down, one per period 1..=6.
    ///
    /// Measured on the span the scanner reports, as classification's own pre-screen
    /// measures it. Below every floor a caller routes on, so a caller can move its floor by
    /// filtering this file rather than re-scanning the genome.
    #[arg(long, value_parser = parse_min_copies, default_value = "5,5,4,4,4,3")]
    pub min_copies: MinCopies,

    /// Sequence required on each side of a tract, to the contig's own end, in bases. A
    /// tract with less than this beside it cannot be anchored by a read.
    #[arg(long, default_value_t = CATALOG_MIN_FLANK_BP)]
    pub min_flank_bp: u64,

    /// Score added when two bases one period apart agree.
    #[arg(long, default_value_t = DEFAULT_MATCH_REWARD, help_heading = "Scoring")]
    pub match_reward: i32,

    /// Score subtracted when they disagree.
    #[arg(long, default_value_t = DEFAULT_MISMATCH_PENALTY, help_heading = "Scoring")]
    pub mismatch_penalty: i32,
}

/// Everything that can stop a build, rendered for a person at a terminal.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum RepeatCatalogCliError {
    /// `--min-period` / `--max-period` do not describe a range.
    #[error("--min-period {min} / --max-period {max}")]
    PeriodRange {
        /// What was asked for.
        min: u8,
        /// What was asked for.
        max: u8,
        /// Why it is not a range.
        #[source]
        source: PeriodRangeError,
    },

    /// A catalog is already there and `--force` was not given.
    #[error(
        "a catalog already exists at {path}; pass --force to replace it \
         (a run reading it while this one rewrites it would fail silently)"
    )]
    Exists {
        /// The file that was left alone.
        path: PathBuf,
    },

    /// The reference could not be read.
    #[error("reading the reference {path}")]
    Reference {
        /// The FASTA.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: ReferenceInfoError,
    },

    /// The catalog could not be written.
    #[error("building the catalog")]
    Build {
        /// What the builder said.
        #[source]
        source: RepeatCatalogError,
    },
}

/// Scan `--reference` and write its catalog. Prints the tally when it finishes.
pub fn run_repeat_catalog(args: &RepeatCatalogArgs) -> Result<(), RepeatCatalogCliError> {
    // The criteria depend on nothing but the args, so a bad flag pair costs no genome read.
    let criteria = catalog_criteria(args)?;
    let scan = ScanParams {
        match_reward: args.match_reward,
        mismatch_penalty: args.mismatch_penalty,
        min_copies: SCANNER_MIN_COPIES,
    };

    let output = args.output.clone().unwrap_or_else(|| {
        let mut path = args.reference.clone().into_os_string();
        path.push(CATALOG_SUFFIX);
        PathBuf::from(path)
    });
    if output.exists() && !args.force {
        return Err(RepeatCatalogCliError::Exists { path: output });
    }

    eprintln!(
        "repeat-catalog: reference={} output={} threads={}",
        args.reference.display(),
        output.display(),
        args.threads
    );

    let mut builder = RepeatCatalogBuilder::create(&output, criteria, scan)
        .map_err(|source| RepeatCatalogCliError::Build { source })?
        .with_threads(args.threads);

    let reference = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: args.reference.clone(),
            fai: None,
        },
        &mut builder,
    )
    .map_err(|source| RepeatCatalogCliError::Reference {
        path: args.reference.clone(),
        source,
    })?;

    let tally = builder
        .finish(&reference, env!("CARGO_PKG_VERSION"))
        .map_err(|source| RepeatCatalogCliError::Build { source })?;

    report(&tally, reference.contigs.len());
    Ok(())
}

/// The criteria a build applies, from the args.
fn catalog_criteria(args: &RepeatCatalogArgs) -> Result<StrRepeatCriteria, RepeatCatalogCliError> {
    let periods = PeriodRange::new(args.min_period, args.max_period).map_err(|source| {
        RepeatCatalogCliError::PeriodRange {
            min: args.min_period,
            max: args.max_period,
            source,
        }
    })?;

    Ok(StrRepeatCriteria {
        classification: SsrSegmentCriteria {
            periods,
            min_copies: args.min_copies,
            // Not applied by the builder — a reader chooses them freely — but recorded in
            // the header as what this run was configured with (spec §4.2).
            min_purity: DEFAULT_MIN_PURITY,
            min_score: DEFAULT_MIN_SCORE,
            bundle_threshold: DEFAULT_BUNDLE_THRESHOLD,
        },
        min_flank_bp: Bp(args.min_flank_bp),
        max_str_len_bp: Bp(CATALOG_MAX_STR_LEN_BP),
    })
}

/// What the run wrote, and what it dropped.
///
/// The rejections are printed beside the rows because a count of what was kept cannot say
/// whether the floors did what they were meant to (spec §4.4) — and the per-period tally is
/// the measurement open question 1 asks for, so the command produces it rather than a
/// separate harness.
fn report(tally: &BuildTally, contigs: usize) {
    println!("contigs\t{contigs}");
    println!("repeats\t{}", tally.rows.total());
    for period in 1..=MAX_MOTIF_LEN as u8 {
        println!("period_{period}\t{}", tally.rows.for_period(period));
    }
    println!("below_copy_floor\t{}", tally.below_copy_floor);
    println!("too_close_to_contig_end\t{}", tally.too_close_to_contig_end);
    if tally.malformed > 0 {
        println!("malformed\t{}", tally.malformed);
    }
}

/// The default catalog floors, as the flag's default string spells them.
///
/// A test rather than a constant expression because clap needs a literal; this is what
/// keeps the literal and [`CATALOG_MIN_COPIES`] from drifting apart.
#[cfg(test)]
fn default_min_copies_flag() -> String {
    crate::ng::repeat_catalog::criteria::CATALOG_MIN_COPIES
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::repeat_catalog::criteria::CATALOG_MIN_COPIES;
    use clap::Parser;

    use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};

    fn args_of(argv: &[&str]) -> RepeatCatalogArgs {
        match Cli::parse_from(argv).cmd {
            PopVarCallerExpCommand::RepeatCatalog(args) => args,
            other => panic!("expected repeat-catalog, got {other:?}"),
        }
    }

    #[test]
    fn the_flag_default_is_the_catalogs_own_floor_table() {
        let args = args_of(&[
            "pop_var_caller_exp",
            "repeat-catalog",
            "--reference",
            "r.fa",
        ]);
        let criteria = catalog_criteria(&args).expect("valid");
        for (index, expected) in CATALOG_MIN_COPIES.iter().enumerate() {
            let period = index as u8 + 1;
            assert_eq!(
                criteria.classification.min_copies.for_period(period),
                *expected,
                "period {period}: the --min-copies default string and CATALOG_MIN_COPIES \
                 have drifted apart (default string is {})",
                default_min_copies_flag()
            );
        }
        assert_eq!(criteria.min_flank_bp, Bp(CATALOG_MIN_FLANK_BP));
        assert_eq!(criteria.classification.periods.min(), CATALOG_MIN_PERIOD);
        assert_eq!(criteria.classification.periods.max(), CATALOG_MAX_PERIOD);
    }

    #[test]
    fn the_output_defaults_to_a_sibling_of_the_reference() {
        let args = args_of(&[
            "pop_var_caller_exp",
            "repeat-catalog",
            "--reference",
            "/data/ref.fa",
        ]);
        assert!(args.output.is_none());
        // The path the driver builds — kept in step with `run_repeat_catalog` by being the
        // same expression.
        let mut expected = args.reference.clone().into_os_string();
        expected.push(CATALOG_SUFFIX);
        assert_eq!(
            PathBuf::from(expected),
            PathBuf::from("/data/ref.fa.repeats.parquet")
        );
    }

    #[test]
    fn an_impossible_period_range_is_an_error_not_a_scan() {
        let args = args_of(&[
            "pop_var_caller_exp",
            "repeat-catalog",
            "--reference",
            "r.fa",
            "--min-period",
            "6",
            "--max-period",
            "2",
        ]);
        assert!(matches!(
            catalog_criteria(&args),
            Err(RepeatCatalogCliError::PeriodRange { min: 6, max: 2, .. })
        ));
    }

    /// End to end: the command scans a reference, writes a catalog beside it, and that
    /// catalog opens against the very reference it was built from.
    #[test]
    fn the_command_writes_a_catalog_that_opens_against_its_reference() {
        use crate::ng::reference_info::{ReferenceSource, read_reference_info};
        use crate::ng::repeat_catalog::{ReadScope, RepeatCatalog};
        use std::io::Write;

        let dir = tempfile::tempdir().expect("tmp");
        let fasta = dir.path().join("ref.fa");
        let mut file = std::fs::File::create(&fasta).expect("create");
        // 200 bases of filler, a clean (CAG)10, 200 more.
        let cycle = b"ACGTTGCAAGGTCCAT";
        let filler: String = (0..200).map(|i| cycle[i % cycle.len()] as char).collect();
        let seq = format!("{filler}{}{filler}", "CAG".repeat(10));
        writeln!(file, ">chr1").expect("header");
        for chunk in seq.as_bytes().chunks(60) {
            file.write_all(chunk).expect("bases");
            file.write_all(b"\n").expect("newline");
        }
        drop(file);

        let args = args_of(&[
            "pop_var_caller_exp",
            "repeat-catalog",
            "--reference",
            fasta.to_str().expect("utf8"),
        ]);
        run_repeat_catalog(&args).expect("the build succeeds");

        let mut catalog_path = fasta.clone().into_os_string();
        catalog_path.push(CATALOG_SUFFIX);
        let catalog_path = PathBuf::from(catalog_path);
        assert!(catalog_path.exists(), "the catalog lands beside the FASTA");

        let reference = read_reference_info(ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        })
        .expect("reads");
        let catalog = RepeatCatalog::open_checking_against_reference(&catalog_path, &reference)
            .expect("the catalog describes the reference it was built from");
        assert_eq!(catalog.contigs().len(), 1);
        assert_eq!(
            catalog.header().tool_version,
            env!("CARGO_PKG_VERSION"),
            "the header records what wrote it"
        );
        assert!(
            catalog
                .repeats(ReadScope::WholeReference)
                .expect("rows")
                .iter()
                .any(|row| row.motif.as_bytes() == b"CAG"),
            "the planted tract is in the file"
        );

        // A second run leaves it alone unless told otherwise — and says so.
        let err = run_repeat_catalog(&args).expect_err("an existing catalog is not replaced");
        assert!(matches!(err, RepeatCatalogCliError::Exists { .. }));
        assert!(format!("{err}").contains("--force"));

        let forced = args_of(&[
            "pop_var_caller_exp",
            "repeat-catalog",
            "--reference",
            fasta.to_str().expect("utf8"),
            "--force",
        ]);
        run_repeat_catalog(&forced).expect("--force replaces it");
    }

    #[test]
    fn a_misspelled_flag_does_not_parse() {
        assert!(
            Cli::try_parse_from([
                "pop_var_caller_exp",
                "repeat-catalog",
                "--reference",
                "r.fa",
                "--min-flank",
                "15",
            ])
            .is_err(),
            "a misspelled knob must fail rather than fall back to a default"
        );
    }
}
