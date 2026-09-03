//! `generate-psps` — psp mode's walk stage at the command line.
//!
//! **A psp is one sample's evidence, stored.** It holds what that sample's reads showed at
//! every position of the analysed ground — the alleles, their support, which reads carried
//! them — in the form the caller consumes, so that calling can happen later without touching
//! the alignment files again.
//!
//! **One sample, one walk, one file.** Each sample's alignment files are read once and its
//! observations stored as a psp; nothing is called and no cohort is assembled
//! (`doc/devel/ng/spec/run_streaming.md` §2, §5.2). What that buys is the reason the psp
//! exists: a sample can be added to a cohort later without re-walking the others, a failed
//! sample is one sample to re-run, and the alignments are decoded exactly once instead of
//! once per cohort they take part in.
//!
//! **The samples are walked one at a time, in the order given** (owner's ruling,
//! 2026-09-03). There is no concurrency knob here and that is deliberate: each sample's
//! generation is independent of every other's, so a cohort is parallelised by running
//! invocations — typically one sample each, on as many machines or cores as there are. That
//! independence is the difference from direct mode, which must hold every sample open at one
//! shared frontier and therefore has to parallelise inside one process.
//!
//! **What this command does *not* write is the census** — the second file spec §2 gives the
//! walk stage, which the parameters fit reads. It lands with Milestone G of
//! `run_driver_psp_mode.md`; until then a psp written here is a complete record of the
//! sample's observations and the fit still has to be fed from elsewhere.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use clap::Args;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::locus_generation::LocusCounts;
use crate::ng::locus_generation::pileup::PileupGeneratorConfig;
use crate::ng::psp::{WriteStats, WriterProvenance};
use crate::ng::read::ReadFilterConfig;
use crate::ng::read::input::read_groups::{
    ReadGroupError, ReadGroups, SampleReadGroups, build_read_groups,
};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, ReferenceInfoError,
    read_reference_verifying_or_creating_fai,
};
use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use crate::ng::run::report::{describe, plural, share_of};
use crate::ng::run::{RunError, SampleObservationGatherer, SampleWalkInputs};
use crate::ng::types::MAX_MOTIF_LEN;
use crate::pop_var_caller::common::{current_command_line, rfc3339_now};
use crate::pop_var_caller_exp::run_ground::{self, GroundError};

#[cfg(test)]
mod tests;

/// What this subcommand is called on the command line, and what every psp it writes records
/// as the subcommand that produced it.
///
/// **One constant for both**, because they must agree: a psp naming a subcommand nobody can
/// type is a dead end for whoever finds the file. `the_recorded_subcommand_is_the_one_a_person_types`
/// ties this to what clap derives from the enum variant.
pub const SUBCOMMAND: &str = "generate-psps";

/// Walk each sample's alignment files once and store its observations as a psp.
#[derive(Debug, Args)]
pub struct GeneratePspsArgs {
    /// Reference FASTA — the one every alignment file was made against. A `.fai` is built
    /// beside it if there is none.
    #[arg(long)]
    pub reference: PathBuf,

    /// The tandem-repeat catalog, which says where the repeat tracts are. Build it first with
    /// `pop_var_caller_exp repeat-catalog --reference <reference>`; it is not optional.
    ///
    /// Defaults to `<reference>.repeats.parquet`, which is where that command writes it. A
    /// catalog built on another reference is refused: its coordinates would put every tract in
    /// the wrong place, genome-wide, with nothing to notice.
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// One alignment file per sample (BAM or CRAM, indexed). Repeat the flag.
    ///
    /// A sample is named by its files' `@RG SM` tag, and **files sharing one `SM` are one
    /// sample and are walked together into one psp**. Naming several samples here is allowed
    /// and writes one psp each, in the order the samples were first seen — but each sample's
    /// walk is independent, so running this command once per sample is the way to spread a
    /// cohort over cores or machines.
    #[arg(long = "alignment", required = true, num_args = 1..)]
    pub alignments: Vec<PathBuf>,

    /// The directory the psp files are written into, one `<sample>.psp` per sample. Created
    /// if it does not exist.
    #[arg(long)]
    pub output_dir: PathBuf,

    /// BED of the stretch of genome to walk. Without it, every base of every contig.
    ///
    /// **This is recorded in every psp written**, and a later calling run compares it across
    /// the cohort: samples analysed over different ground are not comparable, so a cohort
    /// whose files disagree about it is refused rather than called over the ground they
    /// happen to share (spec §6.2). Walking a cohort in several invocations therefore means
    /// giving every one of them the same `--regions`.
    #[arg(long)]
    pub regions: Option<PathBuf>,

    /// Overwrite psps that are already in `--output-dir`.
    ///
    /// Without it a run refuses as soon as it finds one, before walking anything: a psp is
    /// hours of decoding, and silently replacing one because a command was re-typed is not a
    /// thing a person can undo. **The refusal comes before the first sample is walked**, so a
    /// cohort is never left half-replaced.
    #[arg(long)]
    pub force: bool,

    /// Build a `.bai`/`.crai` beside any alignment file that has none.
    #[arg(long, help_heading = "Advanced")]
    pub build_index_if_missing: bool,

    /// The fewest motif copies a tract needs before this run treats it as a repeat: six
    /// comma-separated numbers, one per period 1 to 6. Any other count is refused.
    ///
    /// The default is the copy count at which each period starts to stutter, measured over a
    /// tomato archive on 2026-08-10. Below its floor, a tract is ordinary sequence and the
    /// SNP/indel path handles it.
    #[arg(
        long,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_copies,
        default_value = "8,6,6,6,5,4",
        help_heading = "What counts as a repeat"
    )]
    pub min_copies: MinCopies,

    /// The shortest repeat unit this run treats as a repeat. 1 puts homopolymers on the
    /// repeat path.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub min_period: u8,

    /// The longest repeat unit this run treats as a repeat. Six is the longest the catalog
    /// holds, and the longest a motif can be.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "What counts as a repeat"
    )]
    pub max_period: u8,

    /// A tract longer than this many bases is a satellite: no generator speaks for it, so no
    /// observation of it reaches the psp.
    ///
    /// A round number at the read-length limit rather than a measured one — with 150 bp reads
    /// a read spans a tract plus an anchor each side only up to about 90 bp, so past 100 the
    /// repeat path has nothing to offer.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_STR_LEN,
        help_heading = "What counts as a repeat"
    )]
    pub max_str_len: u64,

    /// How much of a tract must match a perfect tiling of its motif, from 0 to 1. Below this
    /// the tract is ordinary sequence.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PURITY,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_purity,
        help_heading = "What counts as a repeat"
    )]
    pub min_purity: f32,
}

/// Everything that can stop a walk, rendered for a person at a terminal.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum GeneratePspsCliError {
    /// The reference could not be read.
    #[error("reading the reference {}", path.display())]
    Reference {
        /// The FASTA.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: ReferenceInfoError,
    },

    /// The reference's FASTA could not be verified against its index.
    ///
    /// **Named apart from the read itself** because the read succeeds and the verification
    /// runs on a second thread: a run reaches this after its contig table is already in hand.
    #[error("verifying the reference {} against its index", path.display())]
    ReferenceVerification {
        /// The FASTA.
        path: PathBuf,
        /// What the verification said.
        #[source]
        source: ReferenceInfoError,
    },

    /// The ground this walk would cover could not be worked out.
    ///
    /// **Transparent**, so the sentence is the shared one every mode renders: these refusals
    /// are `call-from-alignments`' too, and the same mistake must not read differently
    /// depending on which command was typed.
    #[error(transparent)]
    Ground(#[from] GroundError),

    /// The alignment files' read groups could not be read.
    #[error("the read groups of this walk's alignment files could not be read")]
    ReadGroups {
        /// What the read said.
        #[source]
        source: ReadGroupError,
    },

    /// The directory the psp files would go in could not be made ready.
    #[error("the output directory {} could not be prepared", path.display())]
    OutputDir {
        /// The directory.
        path: PathBuf,
        /// What the filesystem said.
        #[source]
        source: std::io::Error,
    },

    /// A sample is named something that cannot be a file name of its own.
    ///
    /// **`@RG SM` is free header text**, and this command derives each psp's name from it, so
    /// a sample called `../elsewhere` would write outside `--output-dir` and one called
    /// `lane/1` would fail at the write after its whole walk had been decoded. Checked for
    /// every sample before the reference is read, because it costs nothing and the
    /// alternative is finding out at the end.
    #[error(
        "the sample name {sample:?}, read from an @RG SM tag, cannot be a file name; \
         a psp is written as <sample>.psp inside --output-dir"
    )]
    SampleNameNotAFileName {
        /// The name the alignment file declared.
        sample: String,
    },

    /// The moment this run started could not be written down.
    ///
    /// **Not silently substituted**, because the timestamp is the one header field two
    /// otherwise identical runs legitimately differ in — the field spec §12.1's byte-identity
    /// oracle is written to exempt. A constant standing in for it would make two runs
    /// identical for the wrong reason. The two other places in this tree that build the same
    /// provenance refuse here too.
    #[error("this run's start time {stamp:?} is not a timestamp a psp header can carry")]
    Timestamp {
        /// What the clock formatted to.
        stamp: String,
    },

    /// A psp for this sample is already in the output directory.
    ///
    /// **Refused before anything is walked**, and never overwritten by accident: the file is
    /// a walk that already happened, and re-typing a command is not a reason to spend it.
    /// `--force` is how a person says they mean it.
    #[error(
        "{sample} already has a psp at {}; pass --force to walk it again and replace it",
        path.display()
    )]
    PspAlreadyThere {
        /// The individual whose psp is already there.
        sample: String,
        /// The file that would be replaced.
        path: PathBuf,
    },

    /// One sample's walk stopped, and that sample has no psp.
    ///
    /// **Names the sample and the file it would have been**, because the loop walks samples
    /// one at a time and the ones before this one finished: what a person needs to know is
    /// which sample to re-run. **The named path holds no partial file** — a walk writes
    /// beside it and renames only once the psp is whole — so a stopped re-walk leaves any
    /// earlier psp at that path untouched.
    #[error("walking {sample} into {} stopped", path.display())]
    Walk {
        /// The individual whose walk stopped.
        sample: String,
        /// The psp it would have been.
        path: PathBuf,
        /// What stopped it.
        #[source]
        source: RunError,
    },
}

/// What one sample's walk produced.
pub struct SampleWalkOutcome {
    /// The individual walked.
    pub sample: String,
    /// The psp now on disk.
    pub psp: PathBuf,
    /// What the store wrote: **loci**, blocks, bytes. `WriteStats::records` counts psp
    /// records, and one record is a locus with its observations inside it — so this is not a
    /// count of observations, which is several times larger.
    pub stats: WriteStats,
    /// What the walk met: segments dispatched, handled, and the two kinds of refused.
    pub counts: LocusCounts,
}

/// **What the walk stage has to say about itself when it finishes.**
///
/// Held as a value rather than printed where it is computed, so that what a run says is
/// something a test can hold — `call-from-alignments` splits its own report the same way,
/// having found that the printing was the one part of it a mutation could change with the
/// suite still green.
pub struct WalkReport {
    /// The ground every sample was walked over, named as a person would name it —
    /// `chr1:1-100`, or how many intervals and the first and last of them.
    pub ground: String,
    /// How many bases of it were asked for.
    pub analysed_bases: u64,
    /// One entry per sample, in the order they were walked.
    pub samples: Vec<SampleWalkOutcome>,
}

impl SampleWalkOutcome {
    /// **The bases this walk was actually handed**, which is not always the bases asked for.
    ///
    /// A repeat tract is typed and walked whole even where a BED cuts one (spec §4.2), so the
    /// three parts can sum past the ask — measured in the sibling report's own review, a BED
    /// of 120 bases inside two tracts charged 240 and dividing by the 120 printed *200.0%*.
    /// Every share is of this sum, so the parts add to a hundred by construction.
    #[must_use]
    pub fn bases_walked(&self) -> u64 {
        self.counts.regions_handled_bp
            + self.counts.unhandled_not_implemented_bp
            + self.counts.unhandled_out_of_scope_bp
    }

    /// The one line that says what this sample's walk produced — **shared by the progress
    /// note printed as the sample finishes and by the report at the end**, so the two cannot
    /// come to say different things about one walk.
    #[must_use]
    pub fn line(&self) -> String {
        let counts = &self.counts;
        let walked = self.bases_walked();
        let mut line = String::new();
        let _ = write!(
            line,
            "{}: {} loci stored, {} bytes at {}",
            self.sample,
            self.stats.records,
            self.stats.bytes,
            self.psp.display(),
        );
        let _ = write!(
            line,
            "; spoke for {} of {} typed region{} ({} of {} bases walked, {})",
            counts.regions_handled,
            counts.regions_in,
            plural(counts.regions_in),
            counts.regions_handled_bp,
            walked,
            share_of(counts.regions_handled_bp, walked),
        );
        if counts.unhandled_not_implemented_bp > 0 {
            let _ = write!(
                line,
                "; not stored — clusters of repeats too close together to have clean flanks: {} bases ({})",
                counts.unhandled_not_implemented_bp,
                share_of(counts.unhandled_not_implemented_bp, walked),
            );
        }
        if counts.unhandled_out_of_scope_bp > 0 {
            let _ = write!(
                line,
                "; not stored — tandem arrays longer than this run types as callable: {} bases ({})",
                counts.unhandled_out_of_scope_bp,
                share_of(counts.unhandled_out_of_scope_bp, walked),
            );
        }
        line
    }
}

impl WalkReport {
    /// The lines a person reads at the end of a run.
    ///
    /// **Two things every line is about**: how much of the ground this walk could speak for,
    /// and how much it stored. A segment the walk met but has no generator for is *not*
    /// silence — it is analysed ground this build cannot yet describe, and a psp that omits
    /// it looks exactly like one where nothing was there. So what a sample could not cover is
    /// named by its two kinds rather than left to a subtraction, and every share is of the
    /// ground the walk was handed rather than of the ground asked for
    /// ([`SampleWalkOutcome::bases_walked`]).
    #[must_use]
    pub fn lines(&self) -> Vec<String> {
        let mut lines = Vec::with_capacity(self.samples.len() + 3);
        lines.push(format!(
            "walked {} sample{} over {} — {} bases asked for",
            self.samples.len(),
            plural(self.samples.len() as u64),
            self.ground,
            self.analysed_bases,
        ));
        for outcome in &self.samples {
            let walked = outcome.bases_walked();
            if walked != self.analysed_bases {
                lines.push(format!(
                    "  (a typed region is walked whole, so {} spoke for {walked} bases; its shares are of that)",
                    outcome.sample,
                ));
            }
            lines.push(format!("  {}", outcome.line()));
        }
        lines.push(format!(
            "{} psp{}, {} bytes in total",
            self.samples.len(),
            plural(self.samples.len() as u64),
            self.samples
                .iter()
                .map(|outcome| outcome.stats.bytes)
                .sum::<u64>(),
        ));
        lines
    }
}

/// This command's flags, as the shared ground assembly asks for them.
fn ground_request(args: &GeneratePspsArgs) -> run_ground::GroundRequest<'_> {
    run_ground::GroundRequest {
        reference: &args.reference,
        catalog: args.catalog.as_deref(),
        regions: args.regions.as_deref(),
        routing: run_ground::RepeatRouting {
            min_copies: args.min_copies,
            min_period: args.min_period,
            max_period: args.max_period,
            max_str_len: args.max_str_len,
            min_purity: args.min_purity,
        },
    }
}

/// Walk every sample this command was given and write one psp each.
///
/// # Errors
///
/// In order, and all of it before a read is decoded: the output directory must be usable,
/// the run's own start time must be a timestamp a header can carry, the alignment files'
/// read groups must be readable and every sample they name must be a usable file name, the
/// reference must read and verify, and the ground and the catalog must answer. After that,
/// the first sample whose walk stops ends the run naming that sample — with every earlier
/// sample's psp finished on disk, and nothing left at the stopped sample's own path.
pub fn run_generate_psps(args: &GeneratePspsArgs) -> Result<(), GeneratePspsCliError> {
    let report = walk_every_sample(args)?;
    for line in report.lines() {
        println!("{line}");
    }
    Ok(())
}

/// The walk itself: every sample of this run, one at a time, and what each produced.
///
/// **Split from the printing** so what a run says about itself is a value a test can hold,
/// rather than something only a terminal ever sees — the same split
/// `call-from-alignments` makes, and for the reason its own note records: the report was the
/// one part of that command a mutation could change with the whole suite still green.
///
/// # Errors
///
/// As [`run_generate_psps`].
fn walk_every_sample(args: &GeneratePspsArgs) -> Result<WalkReport, GeneratePspsCliError> {
    // **Everything a person typed is judged before a byte is read**, the same order direct
    // mode uses: reading a reference and opening a cohort of CRAMs is minutes, and a path or
    // a number that was never going to work should cost none of them. So the output
    // directory and the run's own timestamp are settled here, at the door — direct mode
    // refuses its output before its reference read for exactly this reason.
    std::fs::create_dir_all(&args.output_dir).map_err(|source| {
        GeneratePspsCliError::OutputDir {
            path: args.output_dir.clone(),
            source,
        }
    })?;
    let provenance = provenance()?;

    // The read-group table is what says how many samples there are and what they are called,
    // and a name this command cannot turn into a file is refused here rather than at the
    // write, after that sample's reads have been decoded.
    let read_groups = build_read_groups(&args.alignments)
        .map_err(|source| GeneratePspsCliError::ReadGroups { source })?;
    for sample in read_groups.read_groups_per_sample() {
        let name = sample.sample.as_ref();
        refuse_a_sample_name_that_is_not_a_file_name(name)?;
        // **Every sample's psp is checked before any sample is walked**, so a cohort whose
        // second psp is already there does not lose its first sample's hours before saying so.
        let psp = psp_path_for(&args.output_dir, name);
        // `try_exists` rather than `exists`: the second answers "no" both when there is no
        // file and when this process cannot tell, and replacing a psp because its directory
        // was briefly unreadable is the failure this check exists to prevent.
        let already_there = psp
            .try_exists()
            .map_err(|source| GeneratePspsCliError::OutputDir {
                path: psp.clone(),
                source,
            })?;
        if !args.force && already_there {
            return Err(GeneratePspsCliError::PspAlreadyThere {
                sample: name.to_string(),
                path: psp,
            });
        }
    }

    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, verify) = read_reference_verifying_or_creating_fai(
        &cache,
        args.reference.clone(),
        ReferenceCheck::VerifyAgainstIndex,
    )
    .map_err(|source| GeneratePspsCliError::Reference {
        path: args.reference.clone(),
        source,
    })?;
    let with_checksums = match verify {
        Some(handle) => {
            handle
                .join()
                .map_err(|source| GeneratePspsCliError::ReferenceVerification {
                    path: args.reference.clone(),
                    source,
                })?
        }
        None => Arc::clone(&info),
    };
    let contigs: ContigList = info.contig_list();
    let reference = OpenReference::new(info);

    let ground = ground_request(args);
    let analysed = run_ground::analysed_regions(&ground, &contigs)?;
    // **One segmentation, shared by every sample's walk** — the same object each gatherer
    // holds a handle on, which is what makes "these psps were written over one ground" true
    // by construction rather than by comparison (spec §4.2).
    let segmentation = Arc::new(run_ground::segments_over(
        &ground,
        &analysed,
        &with_checksums,
    )?);

    // **The samples are walked one at a time, in the order given** (spec §5.2). Nothing here
    // holds two samples' files open at once, which is the property psp mode exists for.
    let mut walked: Vec<SampleWalkOutcome> =
        Vec::with_capacity(read_groups.read_groups_per_sample().len());
    for sample in read_groups.read_groups_per_sample() {
        let files = alignment_files_of(sample, &read_groups);
        let psp_path = psp_path_for(&args.output_dir, sample.sample.as_ref());
        // **Written beside the psp and renamed once it is whole.** `PspWriter::create`
        // truncates what it finds, and this command's advertised repair is *re-run the one
        // sample that failed* — so writing straight to the final path would destroy a good
        // psp at the first byte and leave a refused stump if the re-walk stopped too. The
        // writer's own doc hands that choice to the caller ("write to a new path and rename
        // if that matters"); it matters here.
        // **Unique to this process**, because this command's own advice is to walk a cohort
        // by running one invocation per sample — and two invocations that named the same
        // sample would otherwise interleave into one scratch file and both rename it into
        // place. With a name of its own each writes its own file, and the rename is atomic,
        // so the loser replaces the winner's psp with an equally whole one instead of a
        // shredded one. It is not a lock: the existence check below is advisory, and two
        // invocations racing on one sample is still a thing not to do.
        let while_writing = psp_path.with_extension(format!("psp.{}.partial", std::process::id()));
        let walk = || -> Result<(WriteStats, LocusCounts), RunError> {
            let gatherer = SampleObservationGatherer::open(
                SampleWalkInputs {
                    alignments: &files,
                    reference: &reference,
                    read_filters: ReadFilterConfig::default(),
                    locus_generator_settings: PileupGeneratorConfig::default(),
                    build_index_if_missing: args.build_index_if_missing,
                },
                Arc::clone(&segmentation),
                provenance.clone(),
            )?;
            gatherer.write_psp(&while_writing)
        };
        let produced = match walk() {
            Ok(produced) => produced,
            Err(source) => {
                // The stump says nothing a reader wants and its name is not a psp's; leaving
                // it would put a file in the output directory no later run should ever open.
                let _ = std::fs::remove_file(&while_writing);
                return Err(GeneratePspsCliError::Walk {
                    sample: sample.sample.to_string(),
                    path: psp_path,
                    source,
                });
            }
        };
        std::fs::rename(&while_writing, &psp_path).map_err(|source| {
            GeneratePspsCliError::OutputDir {
                path: psp_path.clone(),
                source,
            }
        })?;
        let (stats, counts) = produced;
        let outcome = SampleWalkOutcome {
            sample: sample.sample.to_string(),
            psp: psp_path,
            stats,
            counts,
        };
        // **Said as each sample finishes, not only at the end**: a cohort of sixty is an hour
        // or more, and a command that says nothing until it is done cannot be told from a
        // hung one. **To stderr**, so a shell capturing the report gets the report and a
        // person watching gets the progress — and in the same words, because both come from
        // `SampleWalkOutcome::line`.
        eprintln!("{}", outcome.line());
        walked.push(outcome);
    }
    Ok(WalkReport {
        ground: describe(&analysed, &contigs),
        analysed_bases: analysed.iter().map(|region| region.len()).sum(),
        samples: walked,
    })
}

/// Refuse a sample whose name cannot be the file name of its own psp.
///
/// **One normal path component and nothing else**: no separator, no `.` or `..`, not empty,
/// not a root or a prefix. `@RG SM` is free header text, so without this a sample could name
/// a path outside `--output-dir` or one that cannot be created at all.
fn refuse_a_sample_name_that_is_not_a_file_name(sample: &str) -> Result<(), GeneratePspsCliError> {
    let refused = || GeneratePspsCliError::SampleNameNotAFileName {
        sample: sample.to_string(),
    };
    let mut components = Path::new(sample).components();
    match (components.next(), components.next()) {
        (Some(std::path::Component::Normal(only)), None) if only == sample => Ok(()),
        _ => Err(refused()),
    }
}

/// The distinct alignment files this sample's read groups live in, in the order the table
/// holds them.
///
/// **Read from the sample's own read-group list**, which the table already grouped, rather
/// than re-scanning every read group and comparing sample names: the by-sample view exists
/// precisely so a caller does not have to do that.
///
/// **Distinct, because a file declaring several of one sample's read groups appears once per
/// group** and handing the same path twice to one walk would open it twice. The order is the
/// table's, which is the order the files were named on the command line — so the psp's own
/// read-group numbering follows what a person typed.
fn alignment_files_of(sample: &SampleReadGroups, read_groups: &ReadGroups) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = Vec::new();
    for id in &sample.read_groups {
        let file = read_groups.get(*id).file.as_ref();
        if !files.iter().any(|seen| seen == file) {
            files.push(file.to_path_buf());
        }
    }
    files
}

/// What only this command can know about the files it writes.
///
/// The gatherer fills in the rest — the input basenames from the files it actually opens, and
/// the read filters it applied — so what is stated here is the program, the subcommand, the
/// command line and the moment (spec §6.1). **The timestamp is the one field two otherwise
/// identical runs legitimately differ in**, which is why §12.1's byte-identity oracle is
/// written to exempt it.
fn provenance() -> Result<WriterProvenance, GeneratePspsCliError> {
    let stamp = rfc3339_now();
    let created = stamp.parse().map_err(|_| GeneratePspsCliError::Timestamp {
        stamp: stamp.clone(),
    })?;
    Ok(WriterProvenance {
        tool: "ng".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        subcommand: SUBCOMMAND.to_string(),
        // Both overwritten by the gatherer from the files it opens.
        input_alignments: Vec::new(),
        input_reference: String::new(),
        command_line: current_command_line(),
        parameters: std::collections::BTreeMap::new(),
        created,
    })
}

/// Where this command writes `sample`'s psp, given the directory it was told to write into.
///
/// **One name, derived, not a flag**: a cohort of psps is opened by naming files, and a
/// per-sample output flag would let two samples be written to one path with nothing noticing.
#[must_use]
pub fn psp_path_for(output_dir: &Path, sample: &str) -> PathBuf {
    output_dir.join(format!("{sample}.psp"))
}
