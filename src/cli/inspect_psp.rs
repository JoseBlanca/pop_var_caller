//! The `inspect-psp` subcommand: what one psp file stored over the regions of a BED, as a TSV
//! with one row per locus and observation.
//!
//! **What it is for.** A psp holds one sample's evidence as the caller reads it. When a variant
//! is missing from a VCF, the first question is whether its reads reached the psp at all, and
//! with how many reads; this answers that question without anybody having to parse the binary
//! format, which changes. Everything after the psp — which candidates were kept, how the
//! genotypes were scored — is the calling run's and is not here.
//!
//! **The output is a contract.** Column names and their order are fixed, and a change to either
//! changes [`COLUMNS_VERSION`]. The header records the psp's own format version, so a script
//! can tell which file and which layout it is reading.
//!
//! # What a row is
//!
//! One row per observation: one distinct sequence a set of reads showed at the locus, from one
//! read group, with whether those reads covered the whole locus. A locus no read showed a
//! sequence at still gets one row, with the observation columns empty, so that every stored
//! locus in the regions appears. A record too large for the reader to hold gets one row from
//! its head alone.
//!
//! # What the counts are and are not
//!
//! **They are the counts the walk stored**, after its own per-position read cap
//! (`reads_discarded_by_cap` says how many that cap dropped at the locus). `call-from-psps`
//! thins each record again to `--max-reads-per-position` as it reads it, so a calling run can
//! have used fewer reads than these.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use crate::psp::PspReadError;
use crate::psp::block::StreamedRecord;
use crate::psp::header::Header;
use crate::psp::reader::PspReader;
use crate::psp::record::RecordHead;
use crate::regions::{BedError, ContigBounds, RegionSet};
use crate::types::{ContigId, GenomePosition, GenomeRegion, Position};

/// The layout of the TSV this writes. **Changes whenever a column is added, removed, renamed or
/// moved**, so a script reading it can refuse a layout it does not know.
pub const COLUMNS_VERSION: u32 = 1;

/// The columns, in order. The names are the contract; see the module note.
pub const COLUMNS: [&str; 22] = [
    "contig",
    "start",
    "end",
    "locus_kind",
    "motif",
    "left_flank",
    "right_flank",
    "reference_bases",
    "non_reference_reads",
    "reads_compared_with_reference",
    "reads_without_observation",
    "reads_discarded_by_cap",
    "record_status",
    "observed_bases",
    "witness",
    "witnessed_positions",
    "read_group",
    "num_obs",
    "num_fwd",
    "q_sum",
    "mapq_mean",
    "mapq_sd",
];

/// How far before a region's first base the walk starts looking for records that reach into
/// it. A psp's index says where each block starts and nothing about how far its records reach,
/// so a record starting before the region and spanning into it — a deletion, a repeat tract —
/// is found only by starting earlier. **Ten thousand bases is far wider than any record the
/// walk writes**: an ordinary locus is bounded by the header's observation reach ceiling, a
/// tract by the repeat-tract length cap, and both are hundreds of bases at most.
const LOOK_BEHIND_BP: u64 = 10_000;

/// Dump what one psp stored over the regions of a BED, one row per locus and observation.
#[derive(Debug, Args)]
pub struct InspectPspArgs {
    /// The psp file to read.
    #[arg(long)]
    pub psp: PathBuf,

    /// BED of the regions to dump. Every stored locus overlapping one is written, including a
    /// locus that starts before the region and reaches into it.
    #[arg(long)]
    pub regions: PathBuf,

    /// Where to write the TSV. Standard output when not given.
    #[arg(long)]
    pub output: Option<PathBuf>,
}

/// Why `inspect-psp` could not finish.
#[derive(Debug, Error)]
pub enum InspectPspCliError {
    /// The psp could not be opened or read.
    #[error("reading the psp {path}")]
    Psp {
        /// The file.
        path: PathBuf,
        /// What the reader hit.
        #[source]
        source: PspReadError,
    },
    /// The BED could not be read against the psp's contigs.
    #[error("reading the regions in {path}")]
    Regions {
        /// The BED.
        path: PathBuf,
        /// What the BED reader hit.
        #[source]
        source: BedError,
    },
    /// A contig in the psp's header is longer than a BED coordinate can say.
    #[error("contig {name} is {length} bases, longer than a BED coordinate can address")]
    ContigTooLong {
        /// The contig.
        name: String,
        /// Its length.
        length: u64,
    },
    /// The output could not be written.
    #[error("writing the output")]
    Write(#[source] io::Error),
}

/// Run `inspect-psp`.
///
/// # Errors
///
/// The psp or the BED cannot be read, or the output cannot be written.
pub fn run_inspect_psp(args: &InspectPspArgs) -> Result<(), InspectPspCliError> {
    let psp_error = |source| InspectPspCliError::Psp {
        path: args.psp.clone(),
        source,
    };
    let mut reader = PspReader::open(&args.psp).map_err(psp_error)?;
    let header = reader.header().clone();
    let bounds = contig_bounds(&header)?;
    let regions = RegionSet::from_bed_path(&args.regions, &bounds).map_err(|source| {
        InspectPspCliError::Regions {
            path: args.regions.clone(),
            source,
        }
    })?;

    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(
            File::create(path).map_err(InspectPspCliError::Write)?,
        )),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };
    write_header(&mut out, &header).map_err(InspectPspCliError::Write)?;

    // **A record is written once even where two regions both reach it**: the BED's regions
    // are sorted and merged, so a record overlapping the previous region was written by that
    // region's walk.
    let mut previous: Option<GenomeRegion> = None;
    for region in regions.regions() {
        let wanted = GenomeRegion {
            contig: ContigId(region.chrom_id),
            start: Position(u64::from(region.start)),
            end: Position(u64::from(region.end)),
        };
        let look_from = GenomePosition {
            contig: wanted.contig,
            position: Position(wanted.start.get().saturating_sub(LOOK_BEHIND_BP).max(1)),
        };
        let mut records = reader.records_from(look_from).map_err(psp_error)?;
        records.skipping_records_too_large_to_hold(true);
        for streamed in records {
            let streamed = streamed.map_err(psp_error)?;
            let at = streamed.head.region;
            if at.contig > wanted.contig || (at.contig == wanted.contig && at.start > wanted.end) {
                break;
            }
            if at.contig < wanted.contig || at.end < wanted.start {
                continue;
            }
            if previous.is_some_and(|before| overlaps(at, before)) {
                continue;
            }
            write_record(&mut out, &header, &streamed).map_err(InspectPspCliError::Write)?;
        }
        previous = Some(wanted);
    }
    out.flush().map_err(InspectPspCliError::Write)
}

/// Whether two stretches of one genome share a base.
fn overlaps(a: GenomeRegion, b: GenomeRegion) -> bool {
    a.contig == b.contig && a.start <= b.end && b.start <= a.end
}

/// The psp's contigs as the BED reader wants them.
fn contig_bounds(header: &Header) -> Result<Vec<ContigBounds<'_>>, InspectPspCliError> {
    header
        .contigs
        .iter()
        .map(|contig| {
            Ok(ContigBounds {
                name: &contig.name,
                length: u32::try_from(contig.length).map_err(|_| {
                    InspectPspCliError::ContigTooLong {
                        name: contig.name.clone(),
                        length: contig.length,
                    }
                })?,
            })
        })
        .collect()
}

/// The `#` lines that say what the file is, then the column names.
fn write_header(out: &mut impl Write, header: &Header) -> io::Result<()> {
    let (major, minor) = header.format_version;
    writeln!(out, "#inspect_psp_columns_version={COLUMNS_VERSION}")?;
    writeln!(out, "#psp_format_version={major}.{minor}")?;
    writeln!(out, "#sample={}", header.sample)?;
    writeln!(
        out,
        "#coordinates=1-based, both ends included, as in a VCF; a locus is written when it \
         overlaps a region of the BED"
    )?;
    writeln!(
        out,
        "#counts=as the walk stored them, after its own read cap; call-from-psps thins each \
         record again to --max-reads-per-position"
    )?;
    writeln!(
        out,
        "#witness=complete: the reads covered every position of the locus; partial: they ran \
         out inside it, and witnessed_positions lists the locus positions they covered, \
         0-based, as start-end runs with the end excluded"
    )?;
    writeln!(
        out,
        "#q_sum=sum over the observation's reads of the natural log of each read's base-error \
         probability; mapq_mean and mapq_sd are over the same reads"
    )?;
    writeln!(out, "{}", COLUMNS.join("\t"))
}

/// One record's rows: one per observation, or one with the observation columns empty.
fn write_record(
    out: &mut impl Write,
    header: &Header,
    streamed: &StreamedRecord,
) -> io::Result<()> {
    let head = &streamed.head;
    let Some(record) = &streamed.record else {
        // Too large to hold: the head is all there is.
        write_locus_columns(out, header, head, None)?;
        return writeln!(out, "\ttoo_large_to_hold\t\t\t\t\t\t\t\t\t");
    };
    if record.observations.is_empty() {
        write_locus_columns(out, header, head, Some(record))?;
        return writeln!(out, "\tdecoded\t\t\t\t\t\t\t\t\t");
    }
    for observation in &record.observations {
        write_locus_columns(out, header, head, Some(record))?;
        let (witness, positions) = match &observation.read_witness {
            ReadWitness::Complete => ("complete", String::new()),
            ReadWitness::Partial { positions } => (
                "partial",
                positions
                    .runs()
                    .map(|(start, end)| format!("{start}-{end}"))
                    .collect::<Vec<_>>()
                    .join(","),
            ),
        };
        let read_group = header
            .read_groups
            .get(observation.read_group.0 as usize)
            .map_or_else(
                || observation.read_group.0.to_string(),
                |group| group.id.clone(),
            );
        let (mapq_mean, mapq_sd) = mapq_mean_and_sd(
            observation.num_obs,
            observation.mapq_sum,
            observation.mapq_sum_sq,
        );
        writeln!(
            out,
            "\tdecoded\t{}\t{witness}\t{positions}\t{read_group}\t{}\t{}\t{:.4}\t{}\t{}",
            String::from_utf8_lossy(&observation.bases),
            observation.num_obs,
            observation.num_fwd,
            observation.q_sum.nats(),
            format_optional(mapq_mean),
            format_optional(mapq_sd),
        )?;
    }
    Ok(())
}

/// The twelve columns every row of a locus shares, without a line end.
fn write_locus_columns(
    out: &mut impl Write,
    header: &Header,
    head: &RecordHead,
    record: Option<&SampleLocusObservations>,
) -> io::Result<()> {
    let region = head.region;
    let contig = header
        .contigs
        .get(region.contig.0 as usize)
        .map_or("?", |contig| contig.name.as_str());
    let (kind, motif, left_flank, right_flank) = match record.map(|record| &record.kind) {
        None => ("", String::new(), String::new(), String::new()),
        Some(LocusKind::Generic) => ("generic", String::new(), String::new(), String::new()),
        Some(LocusKind::Ssr(detail)) => (
            "repeat_tract",
            String::from_utf8_lossy(detail.motif.as_bytes()).into_owned(),
            String::from_utf8_lossy(&detail.left_flank).into_owned(),
            String::from_utf8_lossy(&detail.right_flank).into_owned(),
        ),
        Some(LocusKind::SsrBundle) => {
            ("repeat_bundle", String::new(), String::new(), String::new())
        }
    };
    let reference_bases = record.map_or_else(String::new, |record| {
        String::from_utf8_lossy(&record.reference_bases).into_owned()
    });
    let reads_without_observation = record.map_or_else(String::new, |record| {
        record.reads_without_observation.to_string()
    });
    write!(
        out,
        "{contig}\t{}\t{}\t{kind}\t{motif}\t{left_flank}\t{right_flank}\t{reference_bases}\t{}\t{}\t{reads_without_observation}\t{}",
        region.start.get(),
        region.end.get(),
        head.non_reference_reads,
        head.reads_compared_with_reference,
        head.reads_discarded_by_cap,
    )
}

/// The mean and the sample standard deviation of an observation's MAPQs, from the sums the
/// psp stores. `None` where they are undefined: no reads for the mean, fewer than two for the
/// deviation.
fn mapq_mean_and_sd(reads: u32, sum: u32, sum_of_squares: u64) -> (Option<f64>, Option<f64>) {
    if reads == 0 {
        return (None, None);
    }
    let n = f64::from(reads);
    let mean = f64::from(sum) / n;
    if reads < 2 {
        return (Some(mean), None);
    }
    // Σx² − n·mean², divided by n − 1; floored at zero against rounding.
    let variance = ((sum_of_squares as f64) - n * mean * mean).max(0.0) / (n - 1.0);
    (Some(mean), Some(variance.sqrt()))
}

fn format_optional(value: Option<f64>) -> String {
    value.map_or_else(String::new, |value| format!("{value:.2}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn mapq_moments_come_back_from_the_sums() {
        // MAPQs 60, 60, 30: mean 50, sample variance ((100 + 100 + 400) / 2) = 300.
        let (mean, sd) = mapq_mean_and_sd(3, 150, 60 * 60 * 2 + 30 * 30);
        assert_eq!(mean, Some(50.0));
        assert!((sd.unwrap() - 300.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn mapq_moments_are_absent_where_undefined() {
        assert_eq!(mapq_mean_and_sd(0, 0, 0), (None, None));
        assert_eq!(mapq_mean_and_sd(1, 60, 3600), (Some(60.0), None));
    }

    /// The fixture cohort walked into psps; returns the cohort (which keeps the files alive) and
    /// the path of `zeta`'s psp, the sample carrying reads.
    fn zetas_psp() -> (crate::cli::test_fixtures::ACohortOnDisk, PathBuf) {
        use crate::cli::generate_psps::{GeneratePspsArgs, psp_path_for, run_generate_psps};
        use crate::region_typing::DEFAULT_MAX_STR_LEN;
        use crate::region_typing::segment_criteria::{
            DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
        };
        let cohort = crate::cli::test_fixtures::a_cohort_on_disk();
        let psps = cohort.directory.path().join("psps");
        run_generate_psps(&GeneratePspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            alignments: cohort.alignments.clone(),
            output_dir: psps.clone(),
            regions: None,
            force: false,
            build_index_if_missing: false,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
            max_reads_per_position: crate::locus_generation::pileup::DEFAULT_MAX_SNP_COLUMN_DEPTH,
        })
        .expect("the cohort walks into psps");
        let zeta = psp_path_for(&psps, "zeta");
        (cohort, zeta)
    }

    fn dump(psp: &Path, bed: &str, directory: &Path) -> String {
        let regions = directory.join("regions.bed");
        let output = directory.join("dump.tsv");
        std::fs::write(&regions, bed).expect("the BED is written");
        run_inspect_psp(&InspectPspArgs {
            psp: psp.to_path_buf(),
            regions,
            output: Some(output.clone()),
        })
        .expect("the dump runs");
        std::fs::read_to_string(output).expect("the dump is read back")
    }

    /// **The command end to end**, over a psp another command wrote: the header names the
    /// versions, the column line is the contract, and every row is as wide as it.
    #[test]
    fn a_dump_names_its_versions_and_every_row_fills_every_column() {
        let (cohort, zeta) = zetas_psp();
        let text = dump(&zeta, "chr1\t0\t100\n", cohort.directory.path());

        assert!(text.contains(&format!("#inspect_psp_columns_version={COLUMNS_VERSION}\n")));
        assert!(text.contains("#psp_format_version=1.1\n"));
        assert!(text.contains("#sample=zeta\n"));
        let mut lines = text.lines().filter(|line| !line.starts_with('#'));
        assert_eq!(lines.next(), Some(COLUMNS.join("\t").as_str()));
        let rows: Vec<&str> = lines.collect();
        assert!(!rows.is_empty(), "zeta's reads put loci on chr1");
        for row in &rows {
            assert_eq!(row.split('\t').count(), COLUMNS.len(), "{row}");
            assert!(row.starts_with("chr1\t"), "{row}");
        }
    }

    /// **A locus is written once however the BED cuts around it**, and a region no locus
    /// overlaps writes nothing.
    #[test]
    fn two_regions_write_each_locus_once_and_an_empty_region_writes_nothing() {
        let (cohort, zeta) = zetas_psp();
        let directory = cohort.directory.path();
        let whole = dump(&zeta, "chr1\t0\t100\n", directory);
        let split = dump(&zeta, "chr1\t0\t30\nchr1\t30\t100\n", directory);
        let rows = |text: &str| text.lines().filter(|line| !line.starts_with('#')).count();
        assert_eq!(rows(&split), rows(&whole));

        // chr2's first ten bases: zeta's chr2 read starts at 40.
        let empty = dump(&zeta, "chr2\t0\t10\n", directory);
        assert_eq!(rows(&empty), 1, "the column line alone");
    }

    #[test]
    fn every_row_has_as_many_columns_as_the_header_names() {
        // The two head-only rows: 12 locus columns, then the rest padded.
        let locus_columns = 12;
        let padded = "\ttoo_large_to_hold\t\t\t\t\t\t\t\t\t";
        assert_eq!(locus_columns + padded.matches('\t').count(), COLUMNS.len());
    }
}
