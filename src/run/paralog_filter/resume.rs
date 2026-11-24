//! **What a run with the hidden-duplication filter keeps beside its spill, so that a later run
//! can score and write from the spill without calling again** (`pvc resume-calling`).
//!
//! The filter's passes two and three read the spill, and besides it only what calling handed
//! them: each sample's coverage histogram, each sample's inbreeding coefficient, what the operator
//! asked of the filter, the VCF header's metadata and the ploidy
//! ([`fit_score_and_write_the_calls_explaining`](super::fit_score_and_write_the_calls_explaining)).
//! Everything but the spill lives only in memory, and the histograms cannot be rebuilt without
//! reading every sample's reads again — so a run that stopped in pass two or three, after days of
//! calling, could only be started over. **This keeps them in the run's scratch directory as soon
//! as calling has ended**, where the spill is, and reads them back.
//!
//! **Two files.** The histograms go to [`HISTOGRAMS_FILE_NAME`], one line a sample in the form
//! [`write_the_histogram`] writes, every cell and the bin width by its bit pattern; the rest goes
//! to [`RESUME_FILE_NAME`], a JSON manifest, with every float as its bit pattern. **The manifest is
//! written last**, each file under a temporary name and renamed into place, so a directory holding
//! the manifest holds a complete resume point.
//!
//! **What it is checked against.** The spill has no version of its own — it is written and read
//! inside one run ([`super::spill`]) — so the manifest records [`SPILL_FORMAT`] and a build whose
//! spill differs refuses to read it rather than misread it. The entry count is checked by the
//! readers as they go, as in a run that never stopped.
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use thiserror::Error;

use super::finish::WhatTheOperatorAskedFor;
use super::pass_two::TargetFdr;
use super::spill::SPILL_FORMAT;
use super::spill_file::{SpillFile, SpillFileError};
use crate::neighbour_evidence::phase::PathAgreement;
use crate::neighbour_evidence::stage::NeighbourEvidenceDecision;
use crate::neighbour_evidence::windows::{MarginThinning, WindowSizes};
use crate::run::cohort_merge::recorded_windows::{read_the_histogram, write_the_histogram};
use crate::types::{InbreedingF, Ploidy};
use crate::vcf::header::{HeaderContig, VcfHeaderMetadata};
use crate::window_coverage::SampleHistogram;

/// The manifest's name in the scratch directory.
pub const RESUME_FILE_NAME: &str = "resume.json";

/// The histograms' name in the scratch directory.
pub const HISTOGRAMS_FILE_NAME: &str = "coverage-histograms.tsv";

/// **The manifest's own layout**, raised whenever a field is added, removed or read differently.
const RESUME_FORMAT: u32 = 1;

/// **Everything a resumed run needs besides the spill**, read back from a scratch directory.
pub struct ResumePoint {
    /// The finished spill, with the count of entries calling appended.
    pub spill: SpillFile,
    /// Each sample's coverage histogram, in the run's sample order.
    pub histograms: Vec<SampleHistogram>,
    /// Each sample's inbreeding coefficient, in the run's sample order.
    pub inbreeding: Vec<InbreedingF>,
    /// What the operator asked of the filter.
    pub asked_for: WhatTheOperatorAskedFor,
    /// The VCF header's metadata, with the neighbour-evidence decision the run made.
    pub metadata: VcfHeaderMetadata,
    /// The run's ploidy.
    pub ploidy: Ploidy,
}

/// **Whether `directory` holds a complete resume point** — its manifest, written last.
#[must_use]
pub fn is_resumable(directory: &Path) -> bool {
    directory.join(RESUME_FILE_NAME).is_file()
}

/// **Writes the resume point beside `spill`**, in the run's scratch directory: the histograms, then
/// the manifest. A spill with no directory of its own — one named beside an output, in tests —
/// keeps none.
///
/// # Errors
///
/// If either file cannot be written.
pub fn save(
    spill: &SpillFile,
    histograms: &[SampleHistogram],
    inbreeding: &[InbreedingF],
    asked_for: WhatTheOperatorAskedFor,
    metadata: &VcfHeaderMetadata,
    ploidy: Ploidy,
) -> Result<(), ResumeError> {
    let Some(directory) = spill.directory() else {
        return Ok(());
    };
    write_then_rename(&directory.join(HISTOGRAMS_FILE_NAME), |out| {
        let mut line = String::new();
        for (sample, histogram) in histograms.iter().enumerate() {
            line.clear();
            write_the_histogram(sample, histogram, &mut line);
            out.write_all(line.as_bytes())?;
        }
        Ok(())
    })?;
    let manifest = Manifest {
        resume_format: RESUME_FORMAT,
        spill_format: SPILL_FORMAT,
        spill_entries: spill.entries_written(),
        ploidy: ploidy.get(),
        inbreeding_bits: inbreeding.iter().map(|f| f.get().to_bits()).collect(),
        hidden_dup: HiddenDupAsked {
            target_fdr_bits: asked_for.target_fdr.get().to_bits(),
            tag_instead_of_dropping: asked_for.tag_instead_of_dropping,
            window_bases: asked_for.window_bases,
        },
        header: HeaderParts::of(metadata),
        refinement: metadata.neighbour_evidence().map(RefinementParts::of),
    };
    write_then_rename(&directory.join(RESUME_FILE_NAME), |out| {
        serde_json::to_writer_pretty(&mut *out, &manifest).map_err(std::io::Error::other)?;
        out.write_all(b"\n")
    })
}

/// **Reads the resume point in `directory`** back into what passes two and three take.
///
/// # Errors
///
/// If the manifest or the histograms are missing, unreadable or malformed, if either was written
/// by a build with another layout of the manifest or of the spill, if a value fails its type's
/// check, or if there is no spill.
pub fn load(directory: &Path) -> Result<ResumePoint, ResumeError> {
    let manifest_path = directory.join(RESUME_FILE_NAME);
    let manifest: Manifest = {
        let file = File::open(&manifest_path).map_err(|source| ResumeError::Read {
            path: manifest_path.clone(),
            source,
        })?;
        serde_json::from_reader(BufReader::new(file)).map_err(|source| ResumeError::Manifest {
            path: manifest_path.clone(),
            source,
        })?
    };
    if manifest.resume_format != RESUME_FORMAT {
        return Err(ResumeError::OtherFormat {
            what: "resume manifest",
            found: manifest.resume_format,
            expected: RESUME_FORMAT,
        });
    }
    if manifest.spill_format != SPILL_FORMAT {
        return Err(ResumeError::OtherFormat {
            what: "spill",
            found: manifest.spill_format,
            expected: SPILL_FORMAT,
        });
    }
    let invalid = |what: &str, why: String| ResumeError::Invalid {
        what: what.to_string(),
        why,
    };
    let ploidy = Ploidy::try_new(manifest.ploidy)
        .map_err(|why| invalid("the ploidy", format!("{why:?}")))?;
    let inbreeding = manifest
        .inbreeding_bits
        .iter()
        .map(|&bits| {
            InbreedingF::try_new(f64::from_bits(bits))
                .map_err(|why| invalid("an inbreeding coefficient", format!("{why:?}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let asked_for = WhatTheOperatorAskedFor {
        target_fdr: TargetFdr::try_new(f64::from_bits(manifest.hidden_dup.target_fdr_bits))
            .map_err(|why| invalid("the target FDR", why.to_string()))?,
        tag_instead_of_dropping: manifest.hidden_dup.tag_instead_of_dropping,
        window_bases: manifest.hidden_dup.window_bases,
    };
    let samples = manifest.header.sample_names.len();
    if inbreeding.len() != samples {
        return Err(invalid(
            "the inbreeding coefficients",
            format!("{} for {samples} samples", inbreeding.len()),
        ));
    }
    let metadata = manifest
        .header
        .into_metadata(manifest.refinement, &inbreeding)?;
    let histograms = read_histograms(&directory.join(HISTOGRAMS_FILE_NAME), samples)?;
    let spill = SpillFile::finished_in(directory, manifest.spill_entries)
        .map_err(|source| ResumeError::Spill { source })?;
    Ok(ResumePoint {
        spill,
        histograms,
        inbreeding,
        asked_for,
        metadata,
        ploidy,
    })
}

/// The histograms' file, one line a sample, every sample once and in order.
fn read_histograms(path: &Path, samples: usize) -> Result<Vec<SampleHistogram>, ResumeError> {
    let file = File::open(path).map_err(|source| ResumeError::Read {
        path: path.to_path_buf(),
        source,
    })?;
    let mut histograms = Vec::with_capacity(samples);
    for (index, line) in BufReader::with_capacity(1 << 20, file).lines().enumerate() {
        let line = line.map_err(|source| ResumeError::Read {
            path: path.to_path_buf(),
            source,
        })?;
        let malformed = |why: String| ResumeError::Histograms {
            path: path.to_path_buf(),
            line: index + 1,
            why,
        };
        let (sample, histogram) = read_the_histogram(&line).map_err(malformed)?;
        if sample != index {
            return Err(malformed(format!(
                "sample {sample} where sample {index} was due"
            )));
        }
        histograms.push(histogram);
    }
    if histograms.len() != samples {
        return Err(ResumeError::Histograms {
            path: path.to_path_buf(),
            line: histograms.len(),
            why: format!("{} histograms for {samples} samples", histograms.len()),
        });
    }
    Ok(histograms)
}

/// Writes `path` under a temporary name beside it, then renames it into place.
fn write_then_rename(
    path: &Path,
    write: impl FnOnce(&mut BufWriter<File>) -> std::io::Result<()>,
) -> Result<(), ResumeError> {
    let mut temporary = path.as_os_str().to_owned();
    temporary.push(".partial");
    let temporary = PathBuf::from(temporary);
    let written = (|| {
        let mut out = BufWriter::new(File::create(&temporary)?);
        write(&mut out)?;
        out.into_inner()
            .map_err(|error| error.into_error())?
            .sync_all()?;
        std::fs::rename(&temporary, path)
    })();
    written.map_err(|source| ResumeError::Write {
        path: path.to_path_buf(),
        source,
    })
}

/// **The manifest, as written.** Plain data, every float as its bit pattern; read back through
/// the types' own checks.
#[derive(Serialize, Deserialize)]
struct Manifest {
    resume_format: u32,
    spill_format: u32,
    spill_entries: u64,
    ploidy: u8,
    inbreeding_bits: Vec<u64>,
    hidden_dup: HiddenDupAsked,
    header: HeaderParts,
    refinement: Option<RefinementParts>,
}

#[derive(Serialize, Deserialize)]
struct HiddenDupAsked {
    target_fdr_bits: u64,
    tag_instead_of_dropping: bool,
    window_bases: u64,
}

#[derive(Serialize, Deserialize)]
struct ContigParts {
    name: String,
    length: u64,
    md5: Option<[u8; 16]>,
}

#[derive(Serialize, Deserialize)]
struct HeaderParts {
    contigs: Vec<ContigParts>,
    sample_names: Vec<String>,
    command_line: String,
    reference_path: String,
    log_file_name: String,
    parameters_source: String,
    low_base_quality_test: bool,
}

impl HeaderParts {
    fn of(metadata: &VcfHeaderMetadata) -> Self {
        Self {
            contigs: metadata
                .contigs()
                .iter()
                .map(|contig| ContigParts {
                    name: contig.name.clone(),
                    length: contig.length,
                    md5: contig.md5,
                })
                .collect(),
            sample_names: metadata.sample_names().to_vec(),
            command_line: metadata.command_line().to_string(),
            reference_path: metadata.reference_path().to_string(),
            log_file_name: metadata.log_file_name().to_string(),
            parameters_source: metadata.parameters_source().to_string(),
            low_base_quality_test: metadata.low_base_quality_test(),
        }
    }

    fn into_metadata(
        self,
        refinement: Option<RefinementParts>,
        inbreeding: &[InbreedingF],
    ) -> Result<VcfHeaderMetadata, ResumeError> {
        let metadata = VcfHeaderMetadata::try_new(
            self.contigs
                .into_iter()
                .map(|contig| HeaderContig {
                    name: contig.name,
                    length: contig.length,
                    md5: contig.md5,
                })
                .collect(),
            self.sample_names,
            self.command_line,
            self.reference_path,
            self.log_file_name,
            self.parameters_source,
        )
        .map_err(|why| ResumeError::Invalid {
            what: "the header".to_string(),
            why: why.to_string(),
        })?;
        let metadata = if self.low_base_quality_test {
            metadata.the_low_base_quality_test_ran()
        } else {
            metadata
        };
        Ok(match refinement {
            Some(parts) => {
                metadata.the_neighbour_evidence_decision(parts.into_decision(inbreeding)?)
            }
            None => metadata,
        })
    }
}

#[derive(Serialize, Deserialize)]
struct RefinementParts {
    ran: bool,
    why: String,
    core_bases: u64,
    margin_bases: u64,
    margin_floor_bases: u64,
    /// The shortest core (`--refine-core-floor`); absent in a manifest written before cores could
    /// be shorter than the margin floor, which is then the core floor too.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    core_floor_bases: Option<u64>,
    budget_bytes: u64,
    /// The far margins' walked records a kb, in thousandths (`--refine-margin-density`); absent
    /// where every record is walked, as in a manifest written before margins could be thinned.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    margin_density_thousandths: Option<u32>,
    path_agreement_bits: u64,
    reads_in_draw: bool,
}

impl RefinementParts {
    fn of(decision: &NeighbourEvidenceDecision) -> Self {
        Self {
            ran: decision.ran,
            why: decision.why.clone(),
            core_bases: decision.sizes.core_bases,
            margin_bases: decision.sizes.margin_bases,
            margin_floor_bases: decision.sizes.margin_floor_bases,
            core_floor_bases: Some(decision.sizes.core_floor_bases),
            budget_bytes: decision.sizes.budget_bytes,
            margin_density_thousandths: decision
                .sizes
                .thinning
                .map(|thinning| thinning.thousandths_per_kb.get()),
            path_agreement_bits: decision.path_agreement.get().to_bits(),
            reads_in_draw: decision.reads_in_draw,
        }
    }

    /// The decision back; its inbreeding is the run's, which the manifest keeps once.
    fn into_decision(
        self,
        inbreeding: &[InbreedingF],
    ) -> Result<NeighbourEvidenceDecision, ResumeError> {
        let path_agreement = PathAgreement::try_new(f64::from_bits(self.path_agreement_bits))
            .ok_or_else(|| ResumeError::Invalid {
                what: "the phase path agreement".to_string(),
                why: "not a share in (0, 1]".to_string(),
            })?;
        Ok(NeighbourEvidenceDecision {
            ran: self.ran,
            why: self.why,
            inbreeding_by_sample: inbreeding.iter().map(|f| f.get()).collect(),
            sizes: WindowSizes {
                core_bases: self.core_bases,
                margin_bases: self.margin_bases,
                margin_floor_bases: self.margin_floor_bases,
                core_floor_bases: self.core_floor_bases.unwrap_or(self.margin_floor_bases),
                budget_bytes: self.budget_bytes,
                thinning: self
                    .margin_density_thousandths
                    .and_then(std::num::NonZeroU32::new)
                    .map(|thousandths_per_kb| MarginThinning { thousandths_per_kb }),
            },
            path_agreement,
            reads_in_draw: self.reads_in_draw,
        })
    }
}

/// **Why a resume point could not be written or read.**
#[derive(Debug, Error)]
pub enum ResumeError {
    /// A file of the resume point could not be written.
    #[error("the resume point's `{}` could not be written", path.display())]
    Write {
        /// The file.
        path: PathBuf,
        /// What the file system said.
        #[source]
        source: std::io::Error,
    },
    /// A file of the resume point could not be read.
    #[error(
        "the resume point's `{}` could not be read: was the run stopped before calling ended?",
        path.display()
    )]
    Read {
        /// The file.
        path: PathBuf,
        /// What the file system said.
        #[source]
        source: std::io::Error,
    },
    /// The manifest is not one this build wrote.
    #[error("the resume manifest `{}` could not be read", path.display())]
    Manifest {
        /// The manifest.
        path: PathBuf,
        /// What the JSON reader said.
        #[source]
        source: serde_json::Error,
    },
    /// The manifest or the spill was written by a build with another layout.
    #[error(
        "the {what} was written in layout {found} and this pvc reads layout {expected}: resume \
         with the pvc that made the run"
    )]
    OtherFormat {
        /// The manifest or the spill.
        what: &'static str,
        /// The layout found.
        found: u32,
        /// The layout this build reads.
        expected: u32,
    },
    /// A histogram line is malformed.
    #[error("line {line} of `{}`: {why}", path.display())]
    Histograms {
        /// The histograms' file.
        path: PathBuf,
        /// The line, from 1.
        line: usize,
        /// What was wrong.
        why: String,
    },
    /// A value failed its type's check.
    #[error("{what} in the resume manifest is not valid: {why}")]
    Invalid {
        /// Which value.
        what: String,
        /// Why.
        why: String,
    },
    /// The spill is not in the directory.
    #[error("the run's spill is not in the scratch directory")]
    Spill {
        /// What opening it said.
        #[source]
        source: SpillFileError,
    },
}

#[cfg(test)]
mod tests;
