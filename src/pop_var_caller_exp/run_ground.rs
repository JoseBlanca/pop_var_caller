//! **What ground a run walks, and how it is cut into segments** — the assembly every mode
//! does before a read is decoded.
//!
//! A run needs two things from the reference and the catalog before it can start: the
//! *analysed regions*, which is the stretch of genome it undertakes to speak for, and the
//! *segmentation*, which is that stretch cut into the pieces each locus generator owns
//! (`doc/devel/ng/spec/run_streaming.md` §4). Neither depends on what the run then does with
//! the observations, so both belong to every mode rather than to one.
//!
//! **Lifted out of `call-from-alignments` when psp mode's walk needed it**
//! (`run_driver_psp_mode.md` step C1, which says to reuse these). The alternative — a second
//! copy in `generate-psps` — would let the two modes drift on what counts as a repeat, and
//! the whole psp-mode design rests on a psp being walkable ground that a calling run can
//! check its own segmentation against (§6.2). One copy is what makes that check meaningful.
//!
//! Every refusal here renders exactly as it did inside the direct-mode command: the messages
//! moved verbatim, and both subcommands' error types carry [`GroundError`] transparently, so
//! a person sees the same sentence whichever command they typed.

use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::region_typing::segment_criteria::{MinCopies, SsrSegmentCriteria};
use crate::ng::region_typing::{GenomeRegions, TypedRegionConfig};
use crate::ng::repeat_catalog::{
    CriteriaRefusal, ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria,
    sibling_catalog_path,
};
use crate::ng::run::{RunError, Segmentation};
use crate::ng::tandem_repeat::{PeriodRange, PeriodRangeError};
use crate::ng::types::Bp;
use crate::regions::{BedError, ContigBounds};

/// **What this run counts as a repeat** — the five flags that say so, as both subcommands
/// spell them.
///
/// Grouped because they travel together into one question: which of the catalog's tracts
/// this run treats as repeat ground and which it leaves to the SNP/indel path.
#[derive(Debug, Clone)]
pub struct RepeatRouting {
    /// The fewest motif copies a tract needs, per period 1 to 6.
    pub min_copies: MinCopies,
    /// The shortest repeat unit treated as a repeat.
    pub min_period: u8,
    /// The longest repeat unit treated as a repeat.
    pub max_period: u8,
    /// A tract longer than this many bases is a satellite.
    pub max_str_len: u64,
    /// How much of a tract must match a perfect tiling of its motif, 0 to 1.
    pub min_purity: f32,
}

/// The files and settings the ground is computed from.
///
/// **Borrowed rather than owned**, because a subcommand's `Args` already holds every one of
/// these and this is read once at startup.
pub struct GroundRequest<'a> {
    /// The reference FASTA, named so a refusal can say which assembly it is about.
    pub reference: &'a Path,
    /// The repeat catalog, or `None` for the file `repeat-catalog` writes beside the
    /// reference.
    pub catalog: Option<&'a Path>,
    /// The BED naming the ground, or `None` for every base of every contig.
    pub regions: Option<&'a Path>,
    /// What this run counts as a repeat. **Owned**, because it is five small values a
    /// subcommand builds from its flags rather than something it already holds.
    pub routing: RepeatRouting,
}

impl GroundRequest<'_> {
    /// The catalog this request names — the flag's path, or the one beside the reference.
    #[must_use]
    pub fn catalog_path(&self) -> PathBuf {
        self.catalog
            .map(Path::to_path_buf)
            .unwrap_or_else(|| sibling_catalog_path(self.reference))
    }
}

/// Everything that can stop a run before it has any ground to walk.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum GroundError {
    /// A contig is longer than the run's coordinates can address.
    ///
    /// **The analysed ground is resolved against contig lengths a `u32` carries**
    /// ([`ContigBounds`]), so a reference whose contig is longer than that cannot be resolved
    /// against — and narrowing it silently would compute the run's ground against a wrong
    /// length. No assembly in use has one; the refusal is what makes that a fact rather than
    /// an assumption.
    #[error(
        "contig {name} of {} is {length} bases, longer than the {limit} a run can address",
        reference.display()
    )]
    ContigTooLong {
        /// The reference.
        reference: PathBuf,
        /// The contig.
        name: String,
        /// Its length.
        length: u64,
        /// The most a run can address.
        limit: u64,
    },

    /// The BED naming the ground could not be read.
    #[error("the regions file {} could not be read", path.display())]
    Bed {
        /// The BED.
        path: PathBuf,
        /// What the reader said.
        #[source]
        source: BedError,
    },

    /// There is no repeat catalog to route this run's ground with.
    #[error(
        "no repeat catalog at {}; build one with `pop_var_caller_exp repeat-catalog --reference {}`",
        path.display(),
        reference.display()
    )]
    MissingCatalog {
        /// Where one was looked for.
        path: PathBuf,
        /// The reference it would be built from.
        reference: PathBuf,
    },

    /// The run's period range runs the wrong way round.
    ///
    /// Both ends are bounded to 1..=6 by clap, so this is reachable only by asking for a
    /// narrowest period wider than the widest.
    #[error("--min-period and --max-period do not make a range")]
    PeriodRange {
        /// Which way the range is wrong.
        #[source]
        source: PeriodRangeError,
    },

    /// The catalog could not be read, or was built on another reference.
    #[error("the repeat catalog {} could not be read", path.display())]
    Catalog {
        /// The catalog.
        path: PathBuf,
        /// What the catalog said.
        #[source]
        source: RepeatCatalogError,
    },

    /// This run asked for repeats the catalog was never built to hold.
    ///
    /// **Not a policy refusal**: the rows below the file's own floors were never written, so
    /// the request cannot be served at all. Either move the named flag back up, or build a
    /// catalog at floors low enough to answer it.
    #[error(
        "{flag} asks for repeats the catalog {} does not hold; raise it, or rebuild the catalog \
         at lower floors",
        path.display()
    )]
    RoutingBelowCatalog {
        /// The flag to move.
        flag: &'static str,
        /// The catalog.
        path: PathBuf,
        /// What the catalog said.
        #[source]
        source: RepeatCatalogError,
    },

    /// The catalog's segment stream failed part-way through building the segmentation.
    #[error("the run's segments could not be built")]
    Segmentation {
        /// What the stream said.
        #[source]
        source: RunError,
    },
}

/// The ground this run walks: the BED it was given, or every base of every contig.
///
/// # Errors
///
/// [`GroundError::ContigTooLong`] for a reference a run cannot address, and
/// [`GroundError::Bed`] for a regions file that will not read.
pub fn analysed_regions(
    request: &GroundRequest<'_>,
    contigs: &ContigList,
) -> Result<GenomeRegions, GroundError> {
    // **The narrowing is refused, not taken.** `ContigBounds` carries a `u32`, and a contig
    // past that would have its ground resolved against a wrong length with nothing to notice —
    // the rule `typed_regions.rs`'s `ContigTooLong` records, and the reason the one precedent
    // in this tree that casts is a test.
    let mut bounds: Vec<ContigBounds<'_>> = Vec::with_capacity(contigs.entries.len());
    for entry in &contigs.entries {
        let length = u32::try_from(entry.length).map_err(|_| GroundError::ContigTooLong {
            reference: request.reference.to_path_buf(),
            name: entry.name.clone(),
            length: entry.length,
            limit: u64::from(u32::MAX),
        })?;
        bounds.push(ContigBounds {
            name: &entry.name,
            length,
        });
    }
    match request.regions {
        Some(bed) => {
            GenomeRegions::from_bed_path(bed, &bounds).map_err(|source| GroundError::Bed {
                path: bed.to_path_buf(),
                source,
            })
        }
        None => Ok(GenomeRegions::whole_contigs(&bounds)),
    }
}

/// **What this run counts as a repeat**, built from the five flags that say so.
///
/// The catalog is a source of *candidates*: it is deliberately built below every calling
/// floor, so that a caller can put its own line anywhere inside that gap by filtering the
/// file rather than re-scanning the genome (`repeat_catalog.md` §4.1). Asking the file with
/// [`StrRepeatCriteria::default()`] instead — which *is* the file's own storage floors —
/// makes everything the file holds an STR locus of the run, and on the human benchmark that
/// routed about seven times more reference to the repeat path than ng's calling floors would
/// (`run_ssr_observations.md` §2).
///
/// **The flank floor is not a flag**, because it is the one axis a reader cannot move
/// downwards: the rows below the file's 15 bp were never written, so the request could not be
/// served. It comes from the conversion, which is where that reasoning lives
/// ([`StrRepeatCriteria::from`]).
///
/// The scan half of the [`TypedRegionConfig`] built here is unread: the conversion takes only
/// the classification rules and the satellite cap, and a run that reads a catalog detects
/// nothing itself.
///
/// # Errors
///
/// [`GroundError::PeriodRange`] when the two period flags do not make a range.
pub fn routing_criteria(routing: &RepeatRouting) -> Result<StrRepeatCriteria, GroundError> {
    // Both ends are already bounded to 1..=6 by clap, so the only way left to fail is a
    // range the wrong way round.
    let periods = PeriodRange::new(routing.min_period, routing.max_period)
        .map_err(|source| GroundError::PeriodRange { source })?;
    Ok(StrRepeatCriteria::from(&TypedRegionConfig {
        max_str_len: Bp(routing.max_str_len),
        criteria: SsrSegmentCriteria {
            periods,
            min_copies: routing.min_copies,
            min_purity: routing.min_purity,
            // Not a flag: the score floor gates the *scanner*'s output, and a catalog reader
            // has no scanner. `SsrSegmentCriteria::default()`'s 0 rejects nothing.
            ..SsrSegmentCriteria::default()
        },
        ..TypedRegionConfig::default()
    }))
}

/// The run's segments: the analysed ground cut into the stretches each generator owns, drawn
/// from the catalog.
///
/// # Errors
///
/// [`GroundError::MissingCatalog`] when there is no catalog to read,
/// [`GroundError::Catalog`] when it will not read or was built on another reference,
/// [`GroundError::RoutingBelowCatalog`] when this run asked for repeats it does not hold, and
/// [`GroundError::Segmentation`] when its segment stream fails part-way.
pub fn segments_over(
    request: &GroundRequest<'_>,
    analysed: &GenomeRegions,
    with_checksums: &ReferenceInfo,
) -> Result<Segmentation, GroundError> {
    let path = request.catalog_path();
    if !path.exists() {
        return Err(GroundError::MissingCatalog {
            path,
            reference: request.reference.to_path_buf(),
        });
    }
    let criteria = routing_criteria(&request.routing)?;
    let catalog = RepeatCatalog::open_checking_against_reference(&path, with_checksums).map_err(
        |source| GroundError::Catalog {
            path: path.clone(),
            source,
        },
    )?;
    let spans: Vec<_> = analysed.iter().collect();
    let segments = catalog
        .genome_segments(&criteria, ReadScope::Regions(&spans))
        .map_err(|source| catalog_error_naming_the_flag(source, &path))?;
    Segmentation::build(
        segments,
        analysed.clone(),
        catalog.header().clone(),
        criteria,
        path,
    )
    .map_err(|source| GroundError::Segmentation { source })
}

/// Render a catalog failure, naming the flag to move when the failure is that this run asked
/// for more than the file holds.
///
/// **The refusal is real and it is not policy** — the rows below the file's floors were never
/// written, so the request cannot be served (`run_ssr_observations.md` §2.3) — but on its own
/// it names two numbers and no way to change either. A person who typed `--min-copies 3,3,3,3,3,3`
/// should be told that flag, not left to infer which of five knobs produced *"period 1: catalog
/// holds tracts of 5 copies and up, reader asked for 3"*.
///
/// **The flank floor has no flag and so has no arm here**: a run pins it at the catalog's
/// own, so a catalog built at a wider flank than 15 bp is a file that has to be rebuilt, which
/// is what the general catalog error already says.
fn catalog_error_naming_the_flag(source: RepeatCatalogError, path: &Path) -> GroundError {
    // Exhaustive on the refusal on purpose: a new bounded axis must not silently inherit the
    // no-flag answer and leave a person hunting five knobs for the one they moved.
    let flag = match &source {
        RepeatCatalogError::CriteriaTooPermissive(refusal) => match refusal {
            CriteriaRefusal::CopyFloor { .. } => Some("--min-copies"),
            // Whichever end reaches outside what was built; `serves` checks the low end
            // first, so a range outside at both ends names `--min-period`.
            CriteriaRefusal::PeriodRange {
                built_min,
                wanted_min,
                ..
            } if wanted_min < built_min => Some("--min-period"),
            CriteriaRefusal::PeriodRange { .. } => Some("--max-period"),
            CriteriaRefusal::MinFlank { .. } => None,
        },
        _ => None,
    };
    match flag {
        Some(flag) => GroundError::RoutingBelowCatalog {
            flag,
            path: path.to_path_buf(),
            source,
        },
        None => GroundError::Catalog {
            path: path.to_path_buf(),
            source,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **A run asking for repeats the file does not hold is told which flag to move.**
    ///
    /// The catalog's own refusal names two numbers and no knob: *"period 1: catalog holds
    /// tracts of 5 copies and up, reader asked for 3"* leaves a person to work out which of
    /// five flags produced it. Each bounded axis maps to the flag that moves it.
    ///
    /// Moved here with the function it tests when the ground assembly was lifted out of
    /// `call-from-alignments` for `generate-psps` to share.
    #[test]
    fn a_request_the_catalog_cannot_serve_names_the_flag_that_made_it() {
        let path = Path::new("ref.fa.repeats.parquet");
        let named = |refusal: CriteriaRefusal| match catalog_error_naming_the_flag(
            RepeatCatalogError::CriteriaTooPermissive(refusal),
            path,
        ) {
            GroundError::RoutingBelowCatalog { flag, .. } => Some(flag),
            GroundError::Catalog { .. } => None,
            other => panic!("expected a catalog refusal, got {other:?}"),
        };

        assert_eq!(
            named(CriteriaRefusal::CopyFloor {
                period: 1,
                built: 5,
                wanted: 3
            }),
            Some("--min-copies"),
        );
        assert_eq!(
            named(CriteriaRefusal::PeriodRange {
                built_min: 2,
                built_max: 6,
                wanted_min: 1,
                wanted_max: 6,
            }),
            Some("--min-period"),
            "the low end is the one outside what was built",
        );
        assert_eq!(
            named(CriteriaRefusal::PeriodRange {
                built_min: 1,
                built_max: 4,
                wanted_min: 1,
                wanted_max: 6,
            }),
            Some("--max-period"),
        );
        assert_eq!(
            named(CriteriaRefusal::MinFlank {
                built: 30,
                wanted: 15
            }),
            None,
            "no flag moves the flank floor, so this is a catalog to rebuild and says so",
        );
    }
}
