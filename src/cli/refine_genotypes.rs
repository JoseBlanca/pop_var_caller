//! **`pvc refine-genotypes`**: the neighbour-evidence pass over a finished VCF
//! (`doc/devel/ng/spec/neighbour_evidence_pass.md` §9 Q4).
//!
//! The same pass a calling run puts in front of its writer, run on a VCF pvc wrote earlier: each
//! sample's genotypes re-called from its own `PL` and the samples it shares stretches of chromosome
//! with. Each sample's starting inbreeding coefficient comes from the run's parameters file. The
//! header is the input's, with the `##neighbourEvidence=` line and `GP`'s declaration saying the
//! pass ran.
//!
//! **A VCF the pass has already refined is refused.** Its `GP` would be the pass's own, and using
//! them as templates would count the cohort's evidence a second time (§9 Q4, §10).
//!
//! **So is a VCF that is already phased** (`neighbour_evidence_phase_output.md` §3) — one whose
//! header declares `PS`, or with a `|` in a `GT`: the pass writes its own `PS`, and the columns it
//! leaves unchanged would keep the input's blocks beside the pass's under the same convention.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;
use std::time::Instant;

use clap::Args;
use md5::{Digest, Md5};
use thiserror::Error;

use crate::cli::calling_run;
use crate::cli::provenance::current_command_line;
use crate::cli::run_output::{self, LogFileError};
use crate::neighbour_evidence::stage::{NeighbourEvidenceSwitch, decide};
use crate::types::{ContigId, GenomePosition, Ploidy, Position};
use crate::vcf::header::{
    HeaderContig, HeaderMetadataError, VcfHeaderMetadata, phased_genotype_declarations,
};
use crate::vcf::{RecordPlace, VcfWriteError, VcfWriter};

/// `GP`'s declaration as calling writes it, and as the pass rewrites it (§4.1).
const CALLING_GP: &str = r#"##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posterior probabilities, VCF genotype order">"#;
const REFINED_GP: &str = r#"##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posterior probabilities after the neighbour-evidence pass, VCF genotype order">"#;

/// `pvc refine-genotypes`'s options.
#[derive(Debug, Args, Clone)]
pub struct RefineGenotypesArgs {
    /// The VCF to refine: one pvc wrote, `.vcf` or `.vcf.gz`.
    #[arg(long)]
    pub vcf: PathBuf,

    /// The parameters file the VCF was called with, for each sample's inbreeding coefficient.
    #[arg(long)]
    pub parameters: PathBuf,

    /// Where to write the refined VCF; `.vcf.gz` or `.vcf.bgz` writes bgzf.
    #[arg(long)]
    pub output: PathBuf,

    /// Threads to share each window's samples among. 0 uses every core.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Replace the output and the log if they are already there.
    #[arg(long)]
    pub force: bool,

    /// Where to write the run's log. Defaults to `pvc_refine_genotypes.log` beside the output.
    #[arg(long)]
    pub log: Option<PathBuf>,

    /// The share of the pass's drawn paths that must agree on how two neighbouring heterozygous
    /// sites are phased for the VCF to phase them together. 1, all of them, is the strictest.
    #[arg(
        long,
        default_value = "1",
        value_parser = calling_run::parse_path_agreement
    )]
    pub phase_min_path_agreement: crate::neighbour_evidence::phase::PathAgreement,
}

/// What stops `pvc refine-genotypes`.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum RefineGenotypesCliError {
    /// The output is there and `--force` was not given.
    #[error("`{path}` is already there; pass --force to replace it")]
    OutputAlreadyThere { path: PathBuf },
    /// The log could not be opened.
    #[error(transparent)]
    Log(#[from] LogFileError),
    /// The input could not be read.
    #[error("could not read `{path}`")]
    Read {
        path: PathBuf,
        source: std::io::Error,
    },
    /// The input is not a VCF this command can refine.
    #[error("`{path}` is not a VCF pvc can refine: {reason}")]
    NotRefinable { path: PathBuf, reason: String },
    /// The input says the pass already ran.
    #[error(
        "`{path}` was already refined by the neighbour-evidence pass (its header says `ran=yes`); \
         refining it again would count the cohort's evidence twice — refine the VCF calling wrote"
    )]
    AlreadyRefined { path: PathBuf },
    /// The input already carries a phase.
    #[error(
        "`{path}` is already phased ({what}); the pass writes its own phase sets — strip the \
         input's first, e.g. `bcftools annotate -x FORMAT/PS` and `bcftools +setGT -- -t a -n u`"
    )]
    AlreadyPhased { path: PathBuf, what: String },
    /// The parameters file could not be read or lacks what the pass needs.
    #[error("the parameters file `{path}`: {reason}")]
    Parameters { path: PathBuf, reason: String },
    /// The cohort is one the pass cannot run on.
    #[error("the neighbour-evidence pass cannot run on `{path}`: {why}")]
    CannotRun { path: PathBuf, why: String },
    /// The output's header metadata was refused.
    #[error(transparent)]
    Header(#[from] HeaderMetadataError),
    /// The output could not be written.
    #[error("could not write `{path}`")]
    Write {
        path: PathBuf,
        source: VcfWriteError,
    },
}

/// The input VCF's header, as read.
struct InputHeader {
    /// Every `##` line, in order.
    meta_lines: Vec<String>,
    /// The `#CHROM` line.
    column_line: String,
    contigs: Vec<HeaderContig>,
    samples: Vec<String>,
}

/// **Runs `pvc refine-genotypes`.**
///
/// # Errors
///
/// See [`RefineGenotypesCliError`].
pub fn run_refine_genotypes(args: &RefineGenotypesArgs) -> Result<(), RefineGenotypesCliError> {
    if args.threads > 0 {
        // A pool already built means a second call in one process; the result is unaffected.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }
    if !args.force && args.output.exists() {
        return Err(RefineGenotypesCliError::OutputAlreadyThere {
            path: args.output.clone(),
        });
    }
    let log = args
        .log
        .clone()
        .unwrap_or_else(|| run_output::beside(&args.output, "pvc_refine_genotypes.log"));
    run_output::open_log(&log, args.force)?;
    let started = Instant::now();

    let read_error = |source| RefineGenotypesCliError::Read {
        path: args.vcf.clone(),
        source,
    };
    let mut lines = open_text(&args.vcf).map_err(read_error)?.lines();
    let header = read_header(&mut lines, args)?;
    if header
        .meta_lines
        .iter()
        .any(|line| line.starts_with("##neighbourEvidence=<ran=yes"))
    {
        return Err(RefineGenotypesCliError::AlreadyRefined {
            path: args.vcf.clone(),
        });
    }
    if header
        .meta_lines
        .iter()
        .any(|line| line.starts_with("##FORMAT=<ID=PS,"))
    {
        return Err(RefineGenotypesCliError::AlreadyPhased {
            path: args.vcf.clone(),
            what: "its header declares PS".to_string(),
        });
    }
    let (inbreeding, ploidy) = read_parameters(args, &header)?;
    let decision = decide(NeighbourEvidenceSwitch::On, &inbreeding, ploidy)
        .with_path_agreement(args.phase_min_path_agreement);
    if !decision.ran {
        return Err(RefineGenotypesCliError::CannotRun {
            path: args.vcf.clone(),
            why: decision.why,
        });
    }
    tracing::info!("neighbour-evidence pass: runs ({})", decision.why);

    let mut out_header: Vec<String> = Vec::with_capacity(header.meta_lines.len() + 2);
    let mut placed = false;
    let mut declared_phase = false;
    for line in &header.meta_lines {
        if line.starts_with("##neighbourEvidence=") {
            continue;
        }
        if !placed
            && (line.starts_with("##INFO=")
                || line.starts_with("##FILTER=")
                || line.starts_with("##FORMAT="))
        {
            out_header.push(format!(
                "##neighbourEvidence={}",
                decision.as_header_value()
            ));
            out_header.push(format!(
                "##refineGenotypesCommandline={}",
                current_command_line()
            ));
            placed = true;
        }
        // `GT` now says it is phased where `PS` is set, and `PS` is declared beside it.
        if line.starts_with("##FORMAT=<ID=GT,") {
            out_header.extend(phased_genotype_declarations());
            declared_phase = true;
            continue;
        }
        out_header.push(if line == CALLING_GP {
            REFINED_GP.to_string()
        } else {
            line.clone()
        });
    }
    if !placed {
        out_header.push(format!(
            "##neighbourEvidence={}",
            decision.as_header_value()
        ));
        out_header.push(format!(
            "##refineGenotypesCommandline={}",
            current_command_line()
        ));
    }
    if !declared_phase {
        out_header.extend(phased_genotype_declarations());
    }
    out_header.push(header.column_line.clone());
    let mut header_text = out_header.join("\n");
    header_text.push('\n');

    let contig_ids: HashMap<String, ContigId> = header
        .contigs
        .iter()
        .enumerate()
        .map(|(index, contig)| {
            (
                contig.name.clone(),
                ContigId(u32::try_from(index).expect("fewer than 2^32 contigs")),
            )
        })
        .collect();
    let metadata = VcfHeaderMetadata::try_new(
        header.contigs,
        header.samples,
        current_command_line(),
        String::new(),
        String::new(),
        String::new(),
    )?
    .the_neighbour_evidence_decision(decision);
    let write_error = |source| RefineGenotypesCliError::Write {
        path: args.output.clone(),
        source,
    };
    let mut writer = VcfWriter::create_with_header(&args.output, metadata, ploidy, &header_text)
        .map_err(write_error)?;
    for line in lines {
        let line = line.map_err(read_error)?;
        if line.is_empty() {
            continue;
        }
        if has_a_phased_genotype(&line) {
            let mut fields = line.split('\t');
            let (contig, position) = (fields.next().unwrap_or("?"), fields.next().unwrap_or("?"));
            return Err(RefineGenotypesCliError::AlreadyPhased {
                path: args.vcf.clone(),
                what: format!("the record at {contig}:{position} has a GT with `|`"),
            });
        }
        let place = place_of_line(&line, &contig_ids).ok_or_else(|| {
            RefineGenotypesCliError::NotRefinable {
                path: args.vcf.clone(),
                reason: format!(
                    "a record names a contig the header does not, or has no position: {}",
                    line.chars().take(80).collect::<String>()
                ),
            }
        })?;
        writer
            .write_line(place, line.as_bytes())
            .map_err(write_error)?;
    }
    let finished = writer.finish().map_err(write_error)?;
    calling_run::report_neighbour_evidence(&finished, &args.output);
    tracing::info!(
        "done: {} records written to {} in {}",
        run_output::with_thousands(finished.records_written),
        calling_run::where_the_calls_went(&args.output),
        run_output::elapsed_text(started.elapsed())
    );
    Ok(())
}

/// Whether a record has a `GT` written with `|`.
fn has_a_phased_genotype(line: &str) -> bool {
    let mut fields = line.split('\t');
    let Some(format) = fields.nth(8) else {
        return false;
    };
    let Some(genotype_at) = format.split(':').position(|key| key == "GT") else {
        return false;
    };
    fields.any(|column| {
        column
            .split(':')
            .nth(genotype_at)
            .is_some_and(|genotype| genotype.contains('|'))
    })
}

/// A text reader over a plain or bgzf-compressed file.
fn open_text(path: &std::path::Path) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    let mut probe = [0_u8; 2];
    let mut peek = File::open(path)?;
    let compressed = peek.read(&mut probe)? == 2 && probe == [0x1f, 0x8b];
    Ok(if compressed {
        Box::new(BufReader::new(noodles_bgzf::io::Reader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

fn read_header(
    lines: &mut std::io::Lines<Box<dyn BufRead>>,
    args: &RefineGenotypesArgs,
) -> Result<InputHeader, RefineGenotypesCliError> {
    let not_refinable = |reason: &str| RefineGenotypesCliError::NotRefinable {
        path: args.vcf.clone(),
        reason: reason.to_string(),
    };
    let mut meta_lines = Vec::new();
    let mut contigs = Vec::new();
    for line in lines.by_ref() {
        let line = line.map_err(|source| RefineGenotypesCliError::Read {
            path: args.vcf.clone(),
            source,
        })?;
        if let Some(column_line) = line.strip_prefix("#CHROM") {
            let samples: Vec<String> = column_line
                .split('\t')
                .skip(9)
                .map(str::to_string)
                .collect();
            if samples.is_empty() {
                return Err(not_refinable("it has no samples"));
            }
            return Ok(InputHeader {
                meta_lines,
                column_line: line,
                contigs,
                samples,
            });
        }
        if !line.starts_with("##") {
            return Err(not_refinable("a record comes before the #CHROM line"));
        }
        if let Some(contig) = line.strip_prefix("##contig=<") {
            let field = |key: &str| {
                contig
                    .trim_end_matches('>')
                    .split(',')
                    .find_map(|entry| entry.strip_prefix(key)?.strip_prefix('='))
            };
            let name = field("ID").ok_or_else(|| not_refinable("a ##contig line has no ID"))?;
            let length = field("length")
                .and_then(|length| length.parse().ok())
                .ok_or_else(|| not_refinable("a ##contig line has no length"))?;
            contigs.push(HeaderContig {
                name: name.to_string(),
                length,
                md5: None,
            });
        }
        meta_lines.push(line);
    }
    Err(not_refinable("it has no #CHROM line"))
}

/// Each sample's inbreeding coefficient, in the VCF's sample order, and the run's ploidy, from
/// the parameters file. A file whose checksum differs from the one the VCF's header records is
/// warned about, not refused: the coefficients only start each sample's fit.
fn read_parameters(
    args: &RefineGenotypesArgs,
    header: &InputHeader,
) -> Result<(Vec<f64>, Ploidy), RefineGenotypesCliError> {
    let refused = |reason: String| RefineGenotypesCliError::Parameters {
        path: args.parameters.clone(),
        reason,
    };
    let text =
        std::fs::read_to_string(&args.parameters).map_err(|error| refused(error.to_string()))?;
    let recorded = header.meta_lines.iter().find_map(|line| {
        let value = line.strip_prefix("##pvcParameters=")?;
        let digest = value.split("(md5 ").nth(1)?;
        Some(digest.trim_end_matches(')').to_string())
    });
    let digest: String = Md5::digest(text.as_bytes())
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    if let Some(recorded) = recorded
        && recorded != digest
    {
        tracing::warn!(
            "the parameters file's checksum is {digest}; the VCF was called with one whose \
             checksum is {recorded}"
        );
    }
    let table: toml::Table = toml::from_str(&text).map_err(|error| refused(error.to_string()))?;
    let ploidy = table
        .get("ploidy")
        .and_then(toml::Value::as_integer)
        .and_then(|ploidy| u8::try_from(ploidy).ok())
        .and_then(|ploidy| Ploidy::try_new(ploidy).ok())
        .ok_or_else(|| refused("it states no ploidy".to_string()))?;
    let by_sample: HashMap<&str, f64> = table
        .get("inbreeding")
        .and_then(|inbreeding| inbreeding.get("by_sample"))
        .and_then(toml::Value::as_array)
        .ok_or_else(|| refused("it has no [inbreeding] by_sample".to_string()))?
        .iter()
        .filter_map(|row| {
            Some((
                row.get("sample")?.as_str()?,
                row.get("inbreeding_coefficient")?
                    .get("value")?
                    .as_float()?,
            ))
        })
        .collect();
    let inbreeding = header
        .samples
        .iter()
        .map(|sample| {
            by_sample.get(sample.as_str()).copied().ok_or_else(|| {
                refused(format!(
                    "it has no inbreeding coefficient for sample `{sample}`"
                ))
            })
        })
        .collect::<Result<Vec<f64>, _>>()?;
    Ok((inbreeding, ploidy))
}

/// Where a record line sits: its contig, its `POS`, and whether it is a repeat tract (`STR` in
/// `INFO`), which is what the writer's order check admits one tie for.
fn place_of_line(line: &str, contig_ids: &HashMap<String, ContigId>) -> Option<RecordPlace> {
    let mut fields = line.split('\t');
    let contig = *contig_ids.get(fields.next()?)?;
    let position: u64 = fields.next()?.parse().ok()?;
    let info = fields.nth(5)?;
    Some(RecordPlace {
        at: GenomePosition {
            contig,
            position: Position(position),
        },
        is_repeat_tract: info.split(';').any(|entry| entry == "STR"),
    })
}

#[cfg(test)]
mod tests;
